!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_svnorms --- Implements Singular Vector Norms
!
! !INTERFACE:
!
module m_svnorms
!
! !USES:
#if defined( SPMD )
      use mod_comm,  only : gid
#define CPP_PRT_PREFIX  if(gid == 0)
#else
#define CPP_PRT_PREFIX
#endif

  use precision

  use m_const, only : c_p      => cpm
  use m_const, only : r_dryair => rgas
  use m_const, only : pstd
  use m_const, only : tref     => tstd
  use m_const, only : kappa
  use m_const, only : zvir
  use m_const, only : alhl       ! latent heat condensation

  use m_sVecDef, only : eps_eer  ! eps from Ehrendorfer, Errico and Raeder (1999)

  use m_parDOT, only : parDOT

  use prognostics, only : dyn_prog
  use prognostics, only : jfirst
  use prognostics, only : jlast
  use prognostics, only : imr
  use prognostics, only : jnp
  use prognostics, only : nl

  use m_mpif90,only : MP_comm_rank
  use mod_comm, only : mp_barrier, mp_reduce_sum, mp_recv_s, mp_send_n

  use m_stdio,   only: stdout
  use m_stdio,   only: stderr
  use m_chars,   only: lowercase

  use m_mpout,only : mpout_log
  use m_die,     only: MP_die

  implicit none

#include "lapack.H"

!
! !PUBLIC MEMBER FUNCTIONS:
!
  PRIVATE

  PUBLIC  norm_zvecsize  
  PUBLIC  norm_zz
  PUBLIC  norm_svecA
  PUBLIC  norm_svecBC
  PUBLIC  norm_svecD

  interface norm_zvecsize ; module procedure zvecsize_    ; end interface 
  interface norm_svecA    ; module procedure norm_svecA_  ; end interface
  interface norm_svecBC   ; module procedure norm_svecBC_ ; end interface
  interface norm_svecD    ; module procedure norm_svecD_  ; end interface
  interface norm_zz       ; module procedure norm_zz_     ; end interface

  interface p_for_norms   ; module procedure p_for_norms_ ; end interface
!
! !DESCRIPTION: This module implements a variety of possible norms
!               to apply to the singular-vectors calculation.
!
! !REVISION HISTORY:
!
!  24Oct2002  Todling   Initial design/interfaces.
!  24Jan2002  Errico    Development of correct TE and KE norms
!  11Apr2005  Elena N.  Changed Hard-coded constants to point to m_const
!  05Jun2007  KIM       MPI Parallelization: norm routines
!
!EOP
!-------------------------------------------------------------------------
  character(len=*), parameter :: myname = 'm_SVnorms'
  real(r8), parameter :: pref = 100.d0*pstd

  integer, dimension(5), parameter :: rc = (/99, &  ! Invalid option
                                             98, &  ! Already initialized
                                             97, &  ! Not yet implemented
                                             96, &  ! Undefined normalization
                                             95 /)  ! Packg. not initialized
  character(len=32), dimension(5), parameter :: rcmsg = (/ &
                     'Invalid option                  ',   &
                     'Already initialized             ',   &
                     'Not yet implemented             ',   &
                     'Undefined normalization         ',   &
                     'Package not initialized         ' /)

CONTAINS

! 
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: zvecsize_ --- determine size of z vector
!
! !INTERFACE:

subroutine zvecsize_ ( comm, ROOT, mynorm, nzvecsize, stat )

! !USES:
! determine size of z vector for Lanczos routine based on norm type

  implicit none

! !INPUT PARAMETERS:

  character(len=*), intent(in) :: mynorm

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:

  integer,  intent(in)  ::  comm
  integer,  intent(in)  ::  ROOT

  integer,  intent(out) ::  nzvecsize
  integer,  intent(out), optional ::  stat

! !DESCRIPTION: Initialized weights and whatever else for norm utilization.
!
! !REVISION HISTORY:
!
!  21Jan2003 Errico    Initial code.
!  17Mar2003 Todling   Generalized character usage.
!  13Dec2007 Todling   Augmented to handle q-component in total energy
!
!EOP
!-------------------------------------------------------------------------
  character(len=len(mynorm)) :: usrnorm
  integer  jt, jj, js, je
  integer  myID, ierr

  call MP_comm_rank(comm,myID,ierr)
      if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

  usrnorm = lowercase ( mynorm )
  if (usrnorm=='te' .or. usrnorm=='we') then
      nzvecsize = 0
          ! start counting for U
      js = jfirst; if ( jfirst==1 ) js = 2
      je = jlast

      jt=je-js+1

      nzvecsize = nzvecsize + imr * jt * nl  ! u
           ! start counting for V
      js = jfirst; if ( jfirst==1  ) js = 2
      je = jlast ; if ( jlast ==jnp) je = jnp-1
      jt=je-js+1

      nzvecsize = nzvecsize + imr * jt * nl  ! v
           ! start counting for T
      js = jfirst; if ( jfirst==1  ) js = 2
      je = jlast ; if ( jlast ==jnp) je = jnp-1
      jt=je-js+1

      nzvecsize = nzvecsize + imr * jt * nl  ! t
      if (jfirst == 1)   nzvecsize = nzvecsize + nl 
      if (jlast  == jnp) nzvecsize = nzvecsize + nl 

           ! start counting for PS
      js = jfirst; if ( jfirst==1  ) js = 2
      je = jlast ; if ( jlast ==jnp) je = jnp-1
      jt=je-js+1

      nzvecsize = nzvecsize + imr * jt       ! ps
      if (jfirst == 1)   nzvecsize = nzvecsize + 1 
      if (jlast  == jnp) nzvecsize = nzvecsize + 1 

      if ( usrnorm == 'we' ) then
           nzvecsize = nzvecsize + imr * jt * nl  ! q
           if (jfirst == 1)   nzvecsize = nzvecsize + nl 
           if (jlast  == jnp) nzvecsize = nzvecsize + nl 
           write(stdout,'(2(a,i8))') ' Control vector size for WE-norm: ', nzvecsize, ' on PE: ', myID
      else
           write(stdout,'(2(a,i8))') ' Control vector size for TE-norm: ', nzvecsize, ' on PE: ', myID
      endif

  else if (usrnorm=='ke') then
          ! start counting for U
      nzvecsize = 0
      js = jfirst; if ( jfirst==1 ) js = 2
      je = jlast
      jt=je-js+1

      nzvecsize = nzvecsize + imr * jt * nl  ! u
           ! start counting for V
      js = jfirst; if ( jfirst==1  ) js = 2
      je = jlast ; if ( jlast ==jnp) je = jnp-1
      jt=je-js+1

      nzvecsize = nzvecsize + imr * jt * nl  ! v
      write(stdout,'(2(a,i8))') ' Control vector size for KE-norm: ', nzvecsize, ' on PE: ', myID
  else
      if(present(stat)) then
        write(stdout,'(a)') ' Unknown norm ... '
        stat = 1
      endif
  endif

end subroutine zvecsize_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Norm_sVecA_ --- Select norm for step A of norm application
!
! !INTERFACE:

subroutine norm_svecA_ ( comm, ROOT, mynorm, x, xpert, nzvecsize, zvector, stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:

  integer,  intent(in) :: comm
  integer,  intent(in) :: ROOT
  character(len=*), intent(in) :: mynorm      ! name of norm
  type(dyn_prog), intent(in) :: x             ! reference state stored as vector
  integer,  intent(in) :: nzvecsize           ! size of zvector
  real(r8), intent(in) :: zvector(nzvecsize)  ! z Lanczos vector

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:

  type(dyn_prog)                  :: xpert    ! perturbation corresponding to x
  integer,  optional, intent(out) :: stat     ! status indicator 

! !DESCRIPTION: This routine calls a specific norm routine given the choice
!               defined in the rc file, for use as norm application step A.
!
!               The norm is applied as possibly 5 distinct operators:
!
!      step A: This is the operator $S_1^{1/2}$ acting on the Lanczos vector z to
!              create the perturbation input xpert for the TLM.  The linear operator
!              $S_1$ defines the norm for the initial constraint, such that 
!              $z^T~z$ corresponds to $xpert^T~S_1~xpert$.
!
!      step BC: This is the operator $S_2^{T~1/2}~S_2^{1/2}$ acting on the 
!               output $ypert$ of the TLM (after any local projection operator $L_p$ is
!               applied, to produce an input to the adm (or the adjoint of the $L_p$
!               operator if any is used). The linear operator $S_2$ defines 
!               the end time norm.  This routine also determines the value of the 
!               end time norm $s= ypert^T~S_2~ypert$.
!               Note that if $S_2$ and $L_p$ are both really diagonal, then $L_p$ need 
!               only be called before $S_2$ is applied.  The nomenclature `BC` 
!               signifies the combine operators  $S_2^{T~1/2}$ and $~S_2^{1/2}$. 
!
!      step D: This is the operator $S_1^{T~1/2}$ acting on the adjoint (dual) 
!              vector corresponding to the initial reference state x, called pert
!              here, to produce the vector $z^T$ which is the output of the 
!              Olsedec operator. xpert is the output from the adm. This step 
!              is the adjoint of step A.  
!
! !REMARKS: see also m\_sVecDef
!
! !REVISION HISTORY:
!
!  24Jan2003  Errico    Initial code 
!  17Mar2003  Todling   Generalized character usage.
!  16May2007  Todling   Introduced dyn_prog.
!  13Dec2007  Todling   Augmented to handle q-component in total energy
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::norm_svecA_'
  character(len=len(mynorm)) :: usrnorm

  usrnorm = lowercase ( mynorm )
  if(present(stat)) stat=0

  if (usrnorm=='te' .or. usrnorm=='ke' .or. usrnorm=='we') then
      call E_normA_ ( comm, ROOT, mynorm, x, xpert, nzvecsize, zvector, stat )
  else
      if(present(stat)) then
        stat = rc(1)
        return
      else
        call MP_die(myname_,trim(rcmsg(1)),rc(1)) 
      endif
  endif

end subroutine norm_svecA_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Norm_sVecBC_ --- Select norm for step C*B of norm application
!
! !INTERFACE:

subroutine norm_svecBC_ ( comm, ROOT, mynorm, x, xpert, ynorm, flag, stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:

  integer, intent(in) :: comm
  integer, intent(in) :: ROOT
  character(len=*), intent(in) :: mynorm   ! name of norm
  character(len=*), intent(in) :: flag     ! flag for usage:'BC' means rescale pert, 
                                           ! otherwise only compute ynorm (pert unchanged)  
  type(dyn_prog), intent(in)   :: x        ! reference state stored as vector

! !INPUT/OUTPUT PARAMETERS:

  type(dyn_prog)               :: xpert    ! perturbation corresponding to y
                                           ! On input, is result of tlm after local 
                                           ! projection operator is applied. On output
                                           ! is result to initialize adm before the
                                           ! adjoint of the local projection operator
                                           ! is applied.                                                

! !OUTPUT PARAMETERS:

  real(r8), intent(out)           ::  ynorm ! value of norm applied to y
  integer,  intent(out), optional ::  stat  ! status indicator 

! !DESCRIPTION: This routine calls a specific norm routine given the choice
!               defined in the rc file, for use as norm application step BC.
!
!      step BC: This is the operator $S_2^{T~1/2}~S_2^{1/2}$ acting on the 
!               output $ypert$ of the TLM (after any local projection operator $L_p$ is
!               applied, to produce an input to the adm (or the adjoint of the $L_p$
!               operator if any is used). The linear operator $S_2$ defines 
!               the end time norm.  This routine also determines the value of the 
!               end time norm $ynorm = ypert^T~S_2~ypert$.
!               Note that if $S_2$ and $L_p$ are both really diagonal, then $L_p$ need 
!               only be called before $S_2$ is applied.  The nomenclature `BC` 
!               signifies the combine operators  $C=S_2^{T~1/2}$ and $B=S_2^{1/2}$. 
!
! !REMARKS: see routine norm\_svecA\_ for further description
!
! !REVISION HISTORY:
!
!  24Jan2003  Errico    Initial code 
!  17Mar2003  Todling   Generalized character usage.
!  16May2007  Todling   Introduced dyn_prog.
!  13Dec2007  Todling   Augmented to handle q-component in total energy
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::norm_svecBC_'
  character(len=len(mynorm)) :: usrnorm

  usrnorm = lowercase ( mynorm )
  if(present(stat)) stat=0

  if (usrnorm=='te' .or. usrnorm=='ke' .or. usrnorm=='we') then
     call E_normBC_ ( comm, ROOT, mynorm, x, xpert, ynorm, flag, stat )
  else
      if(present(stat)) then
        stat = rc(1)
        return
      else
        call MP_die(myname_,trim(rcmsg(1)),rc(1)) 
      endif
  endif

end subroutine norm_svecBC_



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Norm_sVecD_ --- Select norm for step D of norm application
!
! !INTERFACE:

subroutine norm_svecD_ ( comm, ROOT, mynorm, x, xpert, nzvecsize, zvector, stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:

  integer, intent(in) :: comm
  integer, intent(in) :: ROOT
  character(len=*), intent(in) :: mynorm    ! name of norm
  type(dyn_prog),   intent(in) :: x         ! reference state stored as vector
  integer,          intent(in) :: nzvecsize ! size of zvector
  type(dyn_prog)               :: xpert     ! result of adm corresponding to x

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:

  real(r8), intent(out) :: zvector(nzvecsize)  ! z Lanczos vector
  integer,  intent(out), optional ::  stat     ! status indicator 

! !DESCRIPTION: This routine calls a specific norm routine given the choice
!               defined in the rc file, for use as norm application step D.
!
!      step D: This is the operator $S_1^{T~1/2}$ acting on the adjoint (dual) 
!              vector corresponding to the initial reference state x, called pert
!              here, to produce the vector $z^T$ which is the output of the 
!              Olsedec operator. xpert is the output from the adm. This step 
!              is the adjoint of step A.  
!
! !REMARKS: see routine norm\_svecA\_ for further description
!
! !REVISION HISTORY:
!
!  24Jan2003  Errico    Initial code 
!  16May2007  Todling   Introduced dyn_prog.
!  13Dec2007  Todling   Augmented to handle q-component in total energy
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::norm_svecD_'
  character(len=len(mynorm)) :: usrnorm

  usrnorm = lowercase ( mynorm )
  if(present(stat)) stat=0

  if (usrnorm=='te' .or. mynorm=='ke' .or. usrnorm=='we') then
      call e_normD_ ( comm, ROOT, mynorm, x, xpert, nzvecsize, zvector, stat )
  else
      if(present(stat)) then
        stat = rc(1)
        return
      else
        call MP_die(myname_,trim(rcmsg(1)),rc(1)) 
      endif
  endif

end subroutine norm_svecD_



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Norm_zz_ --- Compute inner product of two z vectors
!
! !INTERFACE:

subroutine norm_zz_ ( comm, ROOT, nzvecsize, z1, z2, product, stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:

  integer, intent (in) :: ROOT
  integer, intent (in) :: comm
  integer,  intent(in) :: nzvecsize       ! size of zvector
  real(r8), intent(in) :: z1(nzvecsize)   ! z1
  real(r8), intent(in) :: z2(nzvecsize)   ! z2

! !INPUT/OUTPUT PARAMETERS:


! !OUTPUT PARAMETERS:

  real(r8), intent(out)           :: product  ! inner product
  integer,  intent(out), optional ::  stat    ! status indicator 

! !DESCRIPTION: This routine computes the inner product of two vectors
!
!
! !REVISION HISTORY:
!
!  22Jan2003  Errico     Initial routine
!  25Aug2006  Todling    Replaced by BLAS call
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::norm_zz_'

  real(r8) _SDOT
  integer :: n   

  if(present(stat)) stat=0

#ifdef SPMD

  product = parDOT ( z1, z2, comm )

#else

  product = _SDOT ( nzvecsize, z1, 1, z2, 1 )

#endif

end subroutine norm_zz_



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: E_NormA_ --- Energy Norm: z to xpert transform
!
! !INTERFACE:

subroutine E_normA_ ( comm, ROOT, mynorm, x, xpert, nzvecsize, zvector, stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:

  integer,          intent(in) :: ROOT
  integer,          intent(in) :: comm
  character(len=*), intent(in) :: mynorm             ! name of norm
  type(dyn_prog),   intent(in) :: x                  ! reference state
  integer,          intent(in) :: nzvecsize          ! size of zvector
  real(r8),         intent(in) :: zvector(nzvecsize) ! z Lanczos vector

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:

  type(dyn_prog)                   :: xpert    ! perturbation corresponding to x
  integer, optional, intent(out)   :: stat     ! status indicator 

! !DESCRIPTION: This routine implements step A of either the 
!               total energy norm or the KE norm.
!               It copys z to xpert and scales xpert by the the inverse of 
!               the sqrt of the 
!               energy coefficients.  The energy is defined following the
!               routine benergy except for a corrected grid dimension 
!               normalization constant and replacing PE by perturbation APE
!               For scalar fields that are defined at the poles, the input 
!               zvector contains only one value (for i=1), since the other 
!               i=2,...,imr values are constrained to be identical to it.
!
!               See norm\_svecA\_ for descripion of the nomenclature A,BC,D
!
! !REVISION HISTORY:
!
!  24Nov2002  Todling    Setting routine's template.
!  22Jan2003  Errico     Designing working routine
!  17Mar2003  Todling    Generalized character usage.
!  15Apr2005  Elena N.   Changed te-norm calculation to count for p_ref
!                        Changed te-norm calculation not to count for q_ref 
!                           (MOIST NORM vs DRY NORM difference) 
!  16May2007  Todling    Introduced dyn_prog.
!  21May2007  Todling    Fix aux_q(norm=ke) memory leak.
!  05Jun2007  Kim        MPI Parallelization
!  13Dec2007  Todling    Add q-component to total energy
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::E_normA_'


  real(r8)   :: scalej                     ! not used
  real(r8)   :: scalei                     ! not used  
  real(r8)   :: pi                         ! pi
  real(r8)   :: dblejnp                    ! r8 version of jnp
  real(r8)   :: dbleimr                    ! r8 version of imr
  real(r8)   :: ptop                       ! specified pressure at model top
  real(r8)   :: pint                       ! not used here  
  real(r8)   :: polefac                    ! gridfac*polefac=area of polar cap
  real(r8)   :: gridfac                    ! area factor for all points
  real(r8)   :: scalefac                   ! scaling for z to x
  real(r8)   :: tfactor                    ! c_p/tref
  real(r8)   :: qfactor                    ! one
  real(r8)   :: pfactor                    ! r_dryair*tref
  real(r8)   :: bksum                      ! sum of bk vertical corrdinates
  real(r8), allocatable ::  sine(:)
  real(r8), allocatable ::  cosp(:)
  real(r8), allocatable ::  sinp(:)
  real(r8), allocatable ::  cose(:) 
  real(r8), allocatable ::  ak(:)
  real(r8), allocatable ::  bk(:)
  real(r8), allocatable ::  dsigma(:,:,:)  ! reference state (delta p)/ps
  real(r8), allocatable ::  pkz(:,:,:)     ! reference state (p/po)**kappa
  real(r8), allocatable ::  ps(:,:)        ! reference state ps (surface p)
  real(r8), allocatable ::  pspert(:,:)    ! perturbed ps
  real(r8), allocatable ::  aux_q(:,:,:)   ! work array for q_ref to convert TV to T


  integer  i        ! longitudinal index
  integer  j        ! latitudinal index
  integer  k        ! vertical index
  integer  nlp1     ! number of levels plus 1
  integer  ks       ! not used
  integer  ierr     ! error indicator
  integer  izcount  ! counter for elements of the z vector
  integer  myID     ! processor ID
  integer  js, je

  character(len=len(mynorm)) :: usrnorm

  call MP_comm_rank(comm,myID,ierr)
      if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

  usrnorm = lowercase ( mynorm )

  if(present(stat)) stat=0

! specify required horizontal grid area-weighting factors

  allocate( sine(jnp), cosp(jnp), sinp(jnp), cose(jnp), stat=ierr )
	if(ierr/=0) call MP_die(myname_,'Alloc(sine,...)',ierr)
  call setrig(imr, jnp, scalej, scalei, cosp, cose, sinp, sine)
  pi = 4.d0 * datan(1.d0)
  dblejnp = dble(jnp)
  dbleimr = dble(imr)
  polefac = dbleimr*(dblejnp-1)*(1.d0+sine(2))/pi
  gridfac = pi/(4.d0*dbleimr*(dblejnp-1.d0))
  deallocate( sine, cose, sinp, stat=ierr )
  if(ierr/=0) call MP_die(myname_,'Dealloc(sine, etc)',ierr)

! specify required vertical grid weighting constant factors

  nlp1=nl+1
  allocate( ak(nlp1), bk(nlp1), stat=ierr )
	if(ierr/=0) call MP_die(myname_,'Alloc(sine,...)',ierr)
  call set_eta(nl, ks, ptop, pint, ak, bk)

! specify grid weights that depend on reference pressure 

  if ( usrnorm == 'ke' ) then  
    allocate( dsigma(imr,jfirst-1:jlast,nl), stat=ierr )
  	if(ierr/=0) call MP_die(myname_,'Alloc(dsigma)',ierr)
    call p_for_norms ( x%delp, r_dryair, c_p, ptop, dsigma, ierr)  
	if(ierr/=0) call MP_die(myname_,'p_for_norms 1',ierr)
  else ! TE and WE norm also need pkz and ps
    allocate( dsigma(imr,jfirst-1:jlast,nl), pkz(imr,jfirst:jlast,nl), ps(imr,jfirst:jlast), stat=ierr )
  	if(ierr/=0) call MP_die(myname_,'Alloc(sine, etc.)',ierr)
    call p_for_norms ( x%delp, r_dryair, c_p, ptop, dsigma, ierr, pkz=pkz, ps=ps )
	if(ierr/=0) call MP_die(myname_,'p_for_norms 2',ierr)
  endif 

! initialize the zvector counter and allocate auxillary array
 
  izcount=0

  allocate( aux_q(imr,jfirst:jlast,nl), stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Alloc(aux_q',ierr)   
!
! zero moisture perturbations since this is a dry norm
    xpert%q(:,jfirst:jlast,:,1) = 0.d0

!!
! map the u-component of z to the model grid, scale it, then map to xpert
!!  
!   zero extraneous grid values
    if (jfirst == 1) xpert%u(:,1,:) = 0.d0 

!   u points next to South Pole first
  if (jfirst == 1) then
      do k=1,nl     
        do i=1,imr 
           izcount=izcount+1
           scalefac=gridfac * ( 0.5d0 * dsigma(i,2,k) * cosp(2) + &
                                     dsigma(1,1,k) * polefac / dbleimr )
           xpert%u(i,2,k)=zvector(izcount) / sqrt(scalefac)
        enddo
      enddo
  endif

!   u points not adjacent to either pole
  js = jfirst
  if (jfirst == 1) js = 3

  je = jlast
  if (jlast == jnp) je = jnp-1

  do j=js,je
    do k=1,nl
      do i=1,imr
        izcount=izcount+1
        scalefac=0.5d0 * gridfac * ( dsigma(i,  j,k) * cosp(  j) +  &
                                     dsigma(i,j-1,k) * cosp(j-1)   )
        xpert%u(i,j,k)=zvector(izcount) / sqrt(scalefac)
      enddo
    enddo
  enddo

!   u points next to North Pole last
  if (jlast == jnp) then
      do k=1,nl     
        do i=1,imr 
           izcount=izcount+1
           scalefac=gridfac * ( 0.5d0 * dsigma(i,jnp-1,k) * cosp(jnp-1) + &
                                        dsigma(1,jnp,k) * polefac / dbleimr )
           xpert%u(i,jnp,k)=zvector(izcount) / sqrt(scalefac)
        enddo
      enddo
  endif



! end of mappings and calculations for u

!!
! map the v-component of z to the model grid, scale it, then map to xpert
!!  
!   zero extraneous points

    if (jfirst == 1)  xpert%v(:,  1,:) = 0.d0 
    if (jlast == jnp) xpert%v(:,jnp,:) = 0.d0 

!   v points next to South Pole first
  if (jfirst == 1) then
      do k=1,nl     
         scalefac=gridfac * ( dsigma(1,1,k) * polefac / dbleimr +  &
                  0.5d0 * cosp(2) * ( dsigma(1,2,k) + dsigma(imr,2,k) ) )
         xpert%v(1,2,k)=zvector(izcount+1) / sqrt(scalefac)
         do i=2,imr 
            scalefac=gridfac * ( dsigma(1,1,k) * polefac / dbleimr +  &
                     0.5d0 * cosp(2) * ( dsigma(i,2,k) + dsigma(i-1,2,k) ) )

            xpert%v(i,2,k)=zvector(izcount+i) / sqrt(scalefac)
         enddo
         izcount=izcount+imr
      enddo
  endif

!   v points not adjacent to either pole
  js = jfirst
  if (jfirst == 1) js = 3

  je = jlast
  if (jlast == jnp) je = jnp-2

  do j=js,je
    do k=1,nl
      scalefac=gridfac *                                             &
              0.5d0 * cosp(j) * ( dsigma(1,j,k) + dsigma(imr,j,k) ) 
      xpert%v(1,j,k)=zvector(izcount+1) / sqrt(scalefac)
      do i=2,imr
        scalefac=gridfac *                                           & 
              0.5d0 * cosp(j) * ( dsigma(i,j,k) + dsigma(i-1,j,k) ) 
        xpert%v(i,j,k)=zvector(izcount+i) / sqrt(scalefac)
      enddo
      izcount=izcount+imr
    enddo
  enddo

!   v points next to North Pole last
  if (jlast == jnp) then  
      do k=1,nl
         scalefac=gridfac * ( dsigma(1,jnp,k) * polefac / dbleimr +  &
                  0.5d0 * cosp(jnp-1) * ( dsigma(1,jnp-1,k) + dsigma(imr,jnp-1,k) ) )
         xpert%v(1,jnp-1,k)=zvector(izcount+1) / sqrt(scalefac)
         do i=2,imr 
            scalefac=gridfac * ( dsigma(1,jnp,k) * polefac / dbleimr +  &
                  0.5d0 * cosp(jnp-1) * ( dsigma(i,jnp-1,k) + dsigma(i-1,jnp-1,k) ) )
            xpert%v(i,jnp-1,k)=zvector(izcount+i) / sqrt(scalefac)
         enddo
         izcount=izcount+imr
      enddo
  endif


! end of mappings and calculations for v

!!
! If KE norm, then finish up here
!!
  if ( usrnorm == 'ke' ) then  
       if (izcount /= nzvecsize) then 
         call MP_die(myname_,'izcount /= nzvecsize',ierr)
       endif
       xpert%pt   = 0.d0
       xpert%delp = 0.d0
       deallocate( dsigma, ak, bk, cosp, aux_q, stat=ierr )
         if(ierr/=0) call MP_die(myname_,'Dealloc(dsigma, etc.)',ierr)
       return
  endif

!!
! map the pt-component of z to the model grid, scale it, then map to xpert
!!  
    aux_q(:,jfirst:jlast,:) = x%q(:,jfirst:jlast,:,1)

    tfactor = c_p / tref

!   theta points at the South Pole first
  if (jfirst == 1) then    
      do k=1,nl 
         izcount=izcount+1
         scalefac=gridfac * tfactor * dsigma(1,1,k) * polefac &
                 * (pkz(1,1,k)/(1.+zvir*aux_q(1,1,k)))**2   
         xpert%pt(1,1,k)=zvector(izcount) / sqrt(scalefac)
         xpert%pt(2:imr,1,k)=xpert%pt(1,1,k)  ! copy extraneous pole elements
      enddo
  endif

!   theta points not at either pole
  js = jfirst
  if (jfirst == 1) js = 2

  je = jlast
  if (jlast == jnp) je = jnp-1

  do j=js,je
    do k=1,nl
      do i=1,imr
        izcount=izcount+1
        scalefac=gridfac * tfactor * dsigma(i,j,k) * cosp(j) * (pkz(i,j,k)/(1.+zvir*aux_q(i,j,k)))**2 
        xpert%pt(i,j,k)=zvector(izcount) / sqrt(scalefac)
      enddo
    enddo
  enddo

!   theta points at the North Pole last
  if (jlast == jnp) then    
      do k=1,nl 
         izcount=izcount+1
         scalefac=gridfac * tfactor * dsigma(1,jnp,k) * polefac &
                 * (pkz(1,jnp,k)/(1.+zvir*aux_q(1,jnp,k)))**2 
         xpert%pt(1,jnp,k)=zvector(izcount) / sqrt(scalefac)
         xpert%pt(2:imr,jnp,k)=xpert%pt(1,jnp,k)  ! copy extraneous pole elements
      enddo
  endif

! end of mappings and calculations for pt

! If so, account for q-component in total energy norm

  if ( usrnorm == 'we' ) then ! in case of doing total wet-energy
!!
! map the q-component of z to the model grid, scale it, then map to xpert
!!  
    qfactor = eps_eer * alhl * alhl / (c_p * tref)

!   theta points at the South Pole first
  if (jfirst == 1) then    
      do k=1,nl 
         izcount=izcount+1
         scalefac=gridfac * qfactor * dsigma(1,1,k) * polefac
         xpert%q(1,1,k,1)=zvector(izcount) / sqrt(scalefac)
         xpert%q(2:imr,1,k,1)=xpert%q(1,1,k,1)  ! copy extraneous pole elements
      enddo
  endif

!   theta points not at either pole
  js = jfirst
  if (jfirst == 1) js = 2

  je = jlast
  if (jlast == jnp) je = jnp-1

  do j=js,je
    do k=1,nl
      do i=1,imr
        izcount=izcount+1
        scalefac=gridfac * qfactor * dsigma(i,j,k) * cosp(j)
        xpert%q(i,j,k,1)=zvector(izcount) / sqrt(scalefac)
      enddo
    enddo
  enddo

!   theta points at the North Pole last
  if (jlast == jnp) then    
      do k=1,nl 
         izcount=izcount+1
         scalefac=gridfac * qfactor * dsigma(1,jnp,k) * polefac
         xpert%q(1,jnp,k,1)=zvector(izcount) / sqrt(scalefac)
         xpert%q(2:imr,jnp,k,1)=xpert%q(1,jnp,k,1)  ! copy extraneous pole elements
      enddo
  endif

! end of mappings and calculations for q
!!
  endif ! < wet energy >

! deallocate some reference-state arrays no longer needed
!!
  deallocate( pkz, dsigma, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(pkz,sigma)',ierr)
  deallocate( aux_q, stat=ierr ) 
    if(ierr/=0) call MP_die(myname_,'Dealloc(aux_q)',ierr)
 
!!
! map the ps-component of z to the model grid, scale it, 
! then compute corresponding perturbations of delp and finally map to xpert
!!  

  pfactor = r_dryair * tref 

  allocate( pspert(imr,jfirst:jlast), stat=ierr )
	if(ierr/=0) call MP_die(myname_,'Alloc(pspert)',ierr)

!   ps points at the South Pole first
  if ( jfirst == 1 ) then
       izcount=izcount+1
       scalefac=gridfac * pfactor * polefac /pref**2
       pspert(1,1)=zvector(izcount) / sqrt(scalefac)
       pspert(2:imr,1)=pspert(1,1)  ! copy pole values to extraneous array elements
  endif


!   ps points not at either pole
  js = jfirst
  if (jfirst == 1) js = 2

  je = jlast
  if (jlast == jnp) je = jnp-1

  do j=js,je
     scalefac=gridfac * pfactor * cosp(j) / pref**2
    do i=1,imr
      izcount=izcount+1
      pspert(i,j)=zvector(izcount) / sqrt(scalefac)
    enddo
  enddo

!   ps points at the North Pole last
  if (jlast == jnp) then
      izcount=izcount+1
      scalefac=gridfac * pfactor * polefac /pref**2
      pspert(1,jnp)=zvector(izcount) / sqrt(scalefac)
      pspert(2:imr,jnp)=pspert(1,jnp)  ! copy pole values to extra array elements 
  endif


!  distribute ps values in the vertical

  do j=jfirst,jlast
    do k=1,nl
      do i=1,imr
        xpert%delp(i,j,k) = pspert(i,j) * ( bk(k+1) - bk(k) ) / ( bk(nl+1) - bk(1) )
      enddo
    enddo
  enddo


! end of mappings and calculations for ps

!!
!  deallocate remaining local arrays
!!

  deallocate( pspert, ps, ak, bk, cosp, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(pspert, etc.)',ierr)
!
!! check that entire z vector is accounted for

  if (izcount /= nzvecsize) then 
       write(stderr,'(2(a,i8))') 'iz = ', izcount, 'nz = ', nzvecsize
       call MP_die(myname_,'izcount /= nzvecsize',ierr)
  endif

end subroutine e_normA_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: E_NormBC_ --- Total Energy Norm operator applied to pert
!
! !INTERFACE:

subroutine e_normBC_ ( comm, ROOT, mynorm, x, xpert, E, opt, stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:
! it is assumed here that the output y vector from the tlm has the 
! same structure as its input x. x and xpert refer to y and ypert here

  integer, intent (in) :: ROOT
  integer, intent (in) :: comm
  character(len=*), intent(in) :: mynorm   ! name of norm
  character(len=*), intent(in) :: opt      ! if = 'BC' xpert -> E     * xpert
                                           ! if = 'SQ' xpert -> E^1/2 * xpert
                                           ! otherwise xpert is left unchanged

  type(dyn_prog),   intent(in) :: x        ! reference state stored as vector

! !INPUT/OUTPUT PARAMETERS:

  type(dyn_prog)               :: xpert    ! perturbation variable corresponding to x

! !OUTPUT PARAMETERS:

  real(r8), intent(out)           ::  E       ! total energy summed
  integer,  intent(out), optional ::  stat    ! status indicator 

! !DESCRIPTION: This routine implements step BC of the either the total energy 
!  or kinetic energy norms, or is used to compute TE or KE from xpert.
!
!               It is the operator $S^{T~1/2}~S^{1/2}$ applied to the output 
!               from the TLM (after any local projection operator is applied).
!               For this u and v, this operator is equivalent to operating with 
!               $S$, but for T and delp involves some other operations accounting
!               for peculiar storage of pole values and that the delp contribution 
!               to TE is determined by the vertical sums of delp alone.
!               If opt = 'BC' then its is assumed that this routine is being used 
!               to calculate SVs, in which case the input xpert is modified before 
!               output.  Otherwise, xpert is not modified and it is assumed that the 
!               purpose here is to only compute TE.  Also used for KE norm.
!
!
! !REVISION HISTORY:
!
!  24Nov2002  Todling    Setting routine's template.
!  22Jan2003  Errico     Designing working routine
!  17Mar2003  Todling    Generalized character usage.
!  15Apr2005  Elena N.   Changed te-norm calculation to count for p_ref
!                        Changed te-norm calculation not to count for q_ref (MOIST NORM vs DRY NORM difference) 
!  28Aug2006  Todling    Added option SQ to return E^1/2 * xpert  
!  16May2007  Todling    Introduced dyn_prog.
!  05Jun2007  Kim        MPI Parallelization
!  13Dec2007  Todling    Add q-component to total energy
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::e_normBC_'

  real(r8)   :: scalej                     ! not used
  real(r8)   :: scalei                     ! not used  
  real(r8)   :: pi                         ! pi
  real(r8)   :: dblejnp                    ! r8 version of jnp
  real(r8)   :: dbleimr                    ! r8 version of imr
  real(r8)   :: dblenl                     ! r8 version of nl
  real(r8)   :: ptop                       ! specified pressure at model top
  real(r8)   :: pint                       ! not used here  
  real(r8)   :: polefac                    ! gridfac*polefac=area of polar cap
  real(r8)   :: gridfac                    ! area factor for all points
  real(r8)   :: scalefac                   ! scaling for z to x
  real(r8)   :: tfactor                    ! c_p/tref
  real(r8)   :: qfactor                    ! one
  real(r8)   :: pfactor                    ! r_dryair*tref
  real(r8)   :: KEU                        ! kinetic energy (u)
  real(r8)   :: KEV                        ! kinetic energy (v)
  real(r8)   :: APET                       ! available potential energy (T)
  real(r8)   :: AMEQ                       ! available moist energy (Q)
  real(r8)   :: APEP                       ! available potential energy (ps)
  real(r8), allocatable ::  sine(:)
  real(r8), allocatable ::  cosp(:)
  real(r8), allocatable ::  sinp(:)
  real(r8), allocatable ::  cose(:) 
  real(r8), allocatable ::  ak(:)
  real(r8), allocatable ::  bk(:)
  real(r8), allocatable ::  dsigma(:,:,:)  ! reference state (delta p)/ps
  real(r8), allocatable ::  pkz(:,:,:)     ! reference state (p/po)**kappa
  real(r8), allocatable ::  ps(:,:)        ! reference state ps (surface p)
  real(r8), allocatable ::  pspert(:,:)    ! perturbed ps
  real(r8), allocatable :: aux  (:,:,:)
  real(r8), allocatable :: aux_q(:,:,:)    ! work array for q_ref to convert TV to T

  integer  i        ! longitudinal index
  integer  j        ! latitudinal index
  integer  k        ! vertical index
  integer  nlp1     ! number of levels plus 1
  integer  ks       ! defines levels which are pure p coordinate ???
  integer  ierr     ! error indicator
  integer  myID     ! processor ID
  integer  js,je

  character(len=len(mynorm)) :: usrnorm

  call MP_comm_rank(comm,myID,ierr)
      if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

  usrnorm = lowercase ( mynorm )

  if(present(stat)) stat=0


! specify required horizontal grid area-weighting factors

  allocate( sine(jnp), cosp(jnp), sinp(jnp), cose(jnp), stat=ierr )
	if(ierr/=0) call MP_die(myname_,'Alloc(sine,...)',ierr)
  call setrig(imr, jnp, scalej, scalei, cosp, cose, sinp, sine)
  pi = 4.d0 * datan(1.d0)
  dblejnp = dble(jnp)
  dbleimr = dble(imr)
  dblenl  = dble(nl)
  polefac = dbleimr*(dblejnp-1)*(1.d0+sine(2))/pi
  gridfac = pi/(4.d0*dbleimr*(dblejnp-1.d0))
  deallocate( sine, cose, sinp, stat=ierr )
  if(ierr/=0) call MP_die(myname_,'Dealloc(sine, etc)',ierr)

! specify required vertical grid weighting constant factors

  nlp1=nl+1
  allocate( ak(nlp1), bk(nlp1), stat=ierr )
	if(ierr/=0) call MP_die(myname_,'Alloc(sine,...)',ierr)
  call set_eta(nl, ks, ptop, pint, ak, bk)

! specify grid weights that depend on reference pressure 

   if ( usrnorm == 'ke' ) then  
    allocate( dsigma(imr,jfirst-1:jlast,nl), stat=ierr )
  	if(ierr/=0) call MP_die(myname_,'Alloc(dsigma)',ierr)
    call p_for_norms ( x%delp, r_dryair, c_p, ptop, dsigma, ierr )  
	if(ierr/=0) call MP_die(myname_,'p_for_norms 1',ierr)
  else ! TE and WE norm also need pkz and ps
    allocate( dsigma(imr,jfirst-1:jlast,nl), pkz(imr,jfirst:jlast,nl), ps(imr,jfirst:jlast), stat=ierr )
  	if(ierr/=0) call MP_die(myname_,'Alloc(sine, etc.)',ierr)
    call p_for_norms ( x%delp, r_dryair, c_p, ptop, dsigma, ierr, pkz=pkz, ps=ps )
	if(ierr/=0) call MP_die(myname_,'p_for_norms 2',ierr)
  endif 

! allocate auxillary array

 
  allocate( aux(imr,jfirst:jlast,nl), stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Alloc(aux)',ierr)
  allocate( aux_q(imr,jfirst:jlast,nl), stat=ierr )		   
    if(ierr/=0) call MP_die(myname_,'Alloc(aux_q)',ierr)  

!
! zero moisture perturbations since this is a dry norm
    if (opt == 'BC' .or. opt == 'SQ') then
      xpert%q = 0.d0
    endif

  aux_q(:,jfirst:jlast,:) = x%q(:,jfirst:jlast,:,1)
!!
! map the u-component of pert to the model grid, scale it, add to E
!!  

  KEU = 0.d0

!   map xpert to aux  

  aux(1:imr,jfirst:jlast,1:nl) = xpert%u(1:imr,jfirst:jlast,1:nl)

!   u points next to South Pole first
  if (jfirst == 1) then
      j=2
      do k=1,nl     
         do i=1,imr 
            scalefac=gridfac * ( 0.5d0 * dsigma(i,2,k) * cosp(2) + &
                                         dsigma(1,1,k) * polefac / dbleimr )
            KEU = KEU + aux(i,j,k) * aux(i,j,k) * scalefac
            aux(i,j,k) = aux(i,j,k) * mysqrt(scalefac)
         enddo
      enddo
  endif

!   u points not adjacent to either pole
  js=jfirst
  if (jfirst == 1) js=3

  je=jlast
  if (jlast == jnp) je=jnp-1

  do j=js,je
    do k=1,nl
      do i=1,imr
        scalefac=0.5d0 * gridfac * ( dsigma(i,  j,k) * cosp(  j) +  &
                                     dsigma(i,j-1,k) * cosp(j-1)   )
        KEU = KEU + aux(i,j,k) * aux(i,j,k) * scalefac
        aux(i,j,k) = aux(i,j,k) * mysqrt(scalefac)
      enddo
    enddo
  enddo

!   u points next to North Pole last
  if (jlast == jnp) then
      j=jnp
      do k=1,nl     
         do i=1,imr 
            scalefac=gridfac * ( 0.5d0 * dsigma(i,jnp-1,k) * cosp(jnp-1) + &
                                         dsigma(1,jnp,k) * polefac / dbleimr )
            KEU = KEU + aux(i,j,k) * aux(i,j,k) * scalefac
            aux(i,j,k) = aux(i,j,k) * mysqrt(scalefac)
         enddo
      enddo
  endif
 
  call mp_reduce_sum(KEU)
!   map aux to xpert  
  if (opt=='BC' .or. opt=='SQ') then
!   zero extraneous points
    if (jfirst == 1) aux(:,  1,:) = 0.d0 
    xpert%u(1:imr,jfirst:jlast,1:nl) = aux(1:imr,jfirst:jlast,1:nl)
  endif

! end of mappings and calculations for u


!!
! map the v-component of pert to the model grid, scale it, then map to z
!!  

  KEV = 0.d0

!   map xpert to aux  

  aux(1:imr,jfirst:jlast,1:nl) = xpert%v(1:imr,jfirst:jlast,1:nl)

!   v points next to South Pole first
  if (jfirst == 1) then
      j=2
      do k=1,nl     
         scalefac=gridfac * ( dsigma(1,1,k) * polefac / dbleimr +  &
                  0.5d0 * cosp(2) * ( dsigma(1,2,k) + dsigma(imr,2,k) ) )
         KEV = KEV + aux(1,j,k) * aux(1,j,k) * scalefac
         aux(1,j,k) = aux(1,j,k) * mysqrt(scalefac)
         do i=2,imr 
            scalefac=gridfac * ( dsigma(1,1,k) * polefac / dbleimr +  &
                     0.5d0 * cosp(2) * ( dsigma(i,2,k) + dsigma(i-1,2,k) ) )
            KEV = KEV + aux(i,j,k) * aux(i,j,k) *scalefac
            aux(i,j,k) = aux(i,j,k) * mysqrt(scalefac)
         enddo
      enddo
  endif

!   v points not adjacent to either pole
  js=jfirst
  if (jfirst == 1) js=3

  je=jlast
  if (jlast == jnp) je=jnp-2

  do j=js,je
    do k=1,nl
      scalefac=gridfac *                                             &
              0.5d0 * cosp(j) * ( dsigma(1,j,k) + dsigma(imr,j,k) ) 
      KEV = KEV + aux(1,j,k) * aux(1,j,k) * scalefac
      aux(1,j,k) = aux(1,j,k) * mysqrt(scalefac)
      do i=2,imr
        scalefac=gridfac *                                           & 
              0.5d0 * cosp(j) * ( dsigma(i,j,k) + dsigma(i-1,j,k) ) 
        KEV = KEV + aux(i,j,k) * aux(i,j,k) * scalefac
        aux(i,j,k) = aux(i,j,k) * mysqrt(scalefac)
      enddo
    enddo
  enddo

!   v points next to North Pole last
  if (jlast == jnp) then
      j=jnp-1
      do k=1,nl  
         scalefac=gridfac * ( dsigma(1,jnp,k) * polefac / dbleimr +  &
                  0.5d0 * cosp(jnp-1) * ( dsigma(1,jnp-1,k) + dsigma(imr,jnp-1,k) ) )
         KEV = KEV + aux(1,j,k) * aux(1,j,k) * scalefac
         aux(1,j,k) = aux(1,j,k) * mysqrt(scalefac)
         do i=2,imr 
            scalefac=gridfac * ( dsigma(1,jnp,k) * polefac / dbleimr +  0.5d0 &
                    * cosp(jnp-1) * ( dsigma(i,jnp-1,k) + dsigma(i-1,jnp-1,k) ) )
            KEV = KEV + aux(i,j,k) * aux(i,j,k) * scalefac
            aux(i,j,k) = aux(i,j,k) * mysqrt(scalefac)
         enddo
      enddo
  endif


  call mp_reduce_sum(KEV)
!   map aux to xpert  

  if (opt=='BC' .or. opt=='SQ') then
!   zero extraneous points
    if (jfirst == 1) aux(:,  1,:) = 0.d0 
    if (jlast == jnp) aux(:,jnp,:) = 0.d0 
    xpert%v(1:imr,jfirst:jlast,1:nl) = aux(1:imr,jfirst:jlast,1:nl)
  endif

! end of mappings and calculations for v

!!
! If KE norm, then finish up here
!!

  if ( usrnorm == 'ke' ) then  
       E = KEU + KEV
       if (opt=='BC' .or. opt=='SQ') then
         xpert%pt  (1:imr,jfirst:jlast,1:nl) = 0.d0
         xpert%delp(1:imr,jfirst:jlast,1:nl) = 0.d0
       endif
       deallocate( dsigma, ak, bk, aux, cosp, stat=ierr )
         if(ierr/=0) call MP_die(myname_,'Dealloc(dsigma, etc.)',ierr)
       return
  endif

!!
! map the pt-component of pert to the model grid, scale it, compute the pt 
! contribution to the enorm, and create the ypert input to the adjoint model.
!
! Average all the i=1,...,imr the pole values of the TLM output, and then set 
! all pole values equal to that average. This should only have the effect of
! removing roundoff differences, since all TLM pole values should already be 
! identical when input to this routine. The energy at the pole is then computed
! as though only a single pole point is defined, with area that of the polar cap.
! The adjoint then distributes the scaled value at the single pole point equally
! to all i=1,...,imr pole points, but with each 1/imr the original value, as 
! though this value is distributed to imr equal sectors at the pole. This is 
! equivalent to 1) setting each TLM point equal to the average but defined on 
! imr sectors; 2) summing the energy contributed by imr sectors, each having 
! 1/imr of the area of the polar cap, so that the sum is identical to considering
! one pole value on a single polar cap grid area, and (3) setting the input pole
! values for the adjoint equal to the average of the sector values, which is 
! the adjoint of what is done for step 1 (setting all points to their the average 
! is a self-adjoint calculation). 
!!  
 
  APET = 0.d0

!   map xpert to aux  

  aux(1:imr,jfirst:jlast,1:nl) = xpert%pt(1:imr,jfirst:jlast,1:nl)
 
    tfactor = c_p / tref

!   theta points at the South Pole first
  if (jfirst == 1) then
      j=1
      do k=1,nl    
         do i=2,imr
            aux(1,j,k) = aux(1,j,k) + aux(i,j,k) 
         enddo
         aux(1,j,k) = aux(1,j,k) / dble(imr)
         scalefac=gridfac * tfactor * dsigma(1,j,k) * polefac &
                 * (pkz(1,j,k)/(1.+zvir*aux_q(1,j,k)))**2     
         APET = APET + aux(1,j,k) * aux(1,j,k) * scalefac 
         aux(1,j,k) = aux(1,j,k) * mysqrt(scalefac,dbleimr) / dbleimr
         aux(2:imr,j,k)=aux(1,j,k)
      enddo
  endif

!   theta points not at either pole
  js=jfirst
  if (jfirst == 1) js=2
  je=jlast
  if (jlast == jnp) je=jnp-1

  do j=js,je
    do k=1,nl
      do i=1,imr
       scalefac=gridfac * tfactor * dsigma(i,j,k) * cosp(j) * (pkz(i,j,k)/(1.+zvir*aux_q(i,j,k)))**2
        APET = APET + aux(i,j,k) * aux(i,j,k) * scalefac
        aux(i,j,k) = aux(i,j,k) * mysqrt(scalefac)
      enddo
    enddo
  enddo

!   theta points at the North Pole last
  if (jlast == jnp) then
      j=jnp
      do k=1,nl     
         do i=2,imr
            aux(1,j,k) = aux(1,j,k) + aux(i,j,k) 
         enddo
         aux(1,j,k) = aux(1,j,k) / dble(imr)
         scalefac=gridfac * tfactor * dsigma(1,j,k) * polefac &
                 * (pkz(1,j,k)/(1.+zvir*aux_q(1,j,k)))**2
         APET = APET + aux(1,j,k) * aux(1,j,k) * scalefac
         aux(1,j,k) = aux(1,j,k) * mysqrt(scalefac,dbleimr) / dbleimr
         aux(2:imr,j,k)=aux(1,j,k)
      enddo
  endif

  call mp_reduce_sum(APET)
!   map aux to xpert  

  if (opt=='BC' .or. opt=='SQ') then
    xpert%pt(1:imr,jfirst:jlast,1:nl) = aux(1:imr,jfirst:jlast,1:nl)
  endif 

! end of mappings and calculations for pt

! If so, account for q-component in total energy norm

  AMEQ = 0.d0

  if ( usrnorm == 'we' ) then

!   map xpert to aux  

  aux(1:imr,jfirst:jlast,1:nl) = xpert%q(1:imr,jfirst:jlast,1:nl,1)
 
    qfactor = eps_eer * alhl * alhl / (c_p * tref)

!   theta points at the South Pole first
  if (jfirst == 1) then
      j=1
      do k=1,nl    
         do i=2,imr
            aux(1,j,k) = aux(1,j,k) + aux(i,j,k) 
         enddo
         aux(1,j,k) = aux(1,j,k) / dble(imr)
         scalefac=gridfac * tfactor * dsigma(1,j,k) * polefac
         AMEQ = AMEQ + aux(1,j,k) * aux(1,j,k) * scalefac 
         aux(1,j,k) = aux(1,j,k) * mysqrt(scalefac,dbleimr) / dbleimr
         aux(2:imr,j,k)=aux(1,j,k)
      enddo
  endif

!   theta points not at either pole
  js=jfirst
  if (jfirst == 1) js=2
  je=jlast
  if (jlast == jnp) je=jnp-1

  do j=js,je
    do k=1,nl
      do i=1,imr
       scalefac=gridfac * qfactor * dsigma(i,j,k) * cosp(j)
        AMEQ = AMEQ + aux(i,j,k) * aux(i,j,k) * scalefac
        aux(i,j,k) = aux(i,j,k) * mysqrt(scalefac)
      enddo
    enddo
  enddo

!   theta points at the North Pole last
  if (jlast == jnp) then
      j=jnp
      do k=1,nl     
         do i=2,imr
            aux(1,j,k) = aux(1,j,k) + aux(i,j,k) 
         enddo
         aux(1,j,k) = aux(1,j,k) / dble(imr)
         scalefac=gridfac * tfactor * dsigma(1,j,k) * polefac
         AMEQ = AMEQ + aux(1,j,k) * aux(1,j,k) * scalefac
         aux(1,j,k) = aux(1,j,k) * mysqrt(scalefac,dbleimr) / dbleimr
         aux(2:imr,j,k)=aux(1,j,k)
      enddo
  endif

  call mp_reduce_sum(AMEQ)

!   map aux to xpert  

  if (opt=='BC' .or. opt=='SQ') then
    xpert%q(1:imr,jfirst:jlast,1:nl,1) = aux(1:imr,jfirst:jlast,1:nl)
  endif 

! end of mappings and calculations for q

  endif ! < wet energy >

!!
! deallocate some reference-state arrays no longer needed
!!
  deallocate( pkz, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(pkz)',ierr)
  deallocate( dsigma, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(dsigma)',ierr)
  deallocate( aux_q, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(dsigma)',ierr)
!!
! map the delp-component of xpert to the model grid, compute ps, scale it, 
! compute APEP
!!  

  APEP =0.d0

!   map xpert to aux  

  aux(1:imr,jfirst:jlast,1:nl) = xpert%delp(1:imr,jfirst:jlast,1:nl)

  pfactor = r_dryair * tref 

  allocate( pspert(imr,jfirst:jlast), stat=ierr )
	if(ierr/=0) call MP_die(myname_,'Alloc(pspert)',ierr)

!  compute pspert  (pspert is simply the sum of the delp) 
!  see the pt calculations above for explanation of treatment of the poles.

  pspert(:,jfirst:jlast) = 0.d0

  do j=jfirst,jlast
    do k=1,nl
      do i=1,imr
        pspert(i,j) = pspert(i,j) + aux(i,j,k) 
      enddo
    enddo
  enddo

!   ps points at the South Pole first
  if (jfirst == 1) then
      j=1
      do i=2,imr
         pspert(1,j) = pspert(1,j) + pspert(i,j) 
      enddo
      pspert(1,j) = pspert(1,j) / dble(imr)
      scalefac=gridfac * pfactor * polefac / pref**2
      APEP = APEP + pspert(1,j) * pspert(1,j) * scalefac
      pspert(1,j) = pspert(1,j) * mysqrt(scalefac,dbleimr/dblenl) / dbleimr
      pspert(2:imr,j)=pspert(1,j)  ! copy pole values to extra array elements 
  endif

!   ps points not at either pole
  js=jfirst
  if (jfirst == 1) js=2
  je=jlast
  if (jlast == jnp) je=jnp-1

  do j=js,je
      scalefac=gridfac * pfactor * cosp(j) / pref**2
    do i=1,imr
      APEP = APEP + pspert(i,j) * pspert(i,j) * scalefac
      pspert(i,j) = pspert(i,j) * mysqrt(scalefac,1.0/dblenl) 
    enddo
  enddo

!   ps points at the North Pole last
  if (jlast == jnp) then
      j=jnp
      do i=2,imr
         pspert(1,j) = pspert(1,j) + pspert(i,j) 
      enddo
      pspert(1,j) = pspert(1,j) / dble(imr)
      scalefac=gridfac * pfactor * polefac / pref**2
      APEP = APEP + pspert(1,j) * pspert(1,j) * scalefac
      pspert(1,j) = pspert(1,j) * mysqrt(scalefac,dbleimr/dblenl) / dbleimr
      pspert(2:imr,j)=pspert(1,j)  ! copy pole values to extra array elements 
  endif

  call mp_reduce_sum(APEP)
!  distribute ps values in the vertical 
  do k=1,nl
    aux(:,jfirst:jlast,k) = pspert(:,jfirst:jlast) 
  enddo

!   map aux to xpert  

  if (opt=='BC' .or. opt=='SQ') then
    xpert%delp(1:imr,jfirst:jlast,1:nl) = aux(1:imr,jfirst:jlast,1:nl)
  endif

! end of mappings and calculations for ps

!!
!  deallocate remaining local arrays
!!

  deallocate( pspert, ps, ak, bk, aux, cosp, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(pspert, etc.)',ierr)

!!
!    Compute total energy
!!
  
  E = KEU + KEV + APET + APEP + AMEQ

! should print out KEU, KEV, APET, APEP, AMEQ if opt /= 'BC' somehow  
  if ( myID==ROOT ) then
       write(stdout,'(1x,a)') 'e_normBC:          KEU     |        KEV       |        APET      |         APEP        |         AMEQ        '
       write(stdout,'(1x,a)') 'e_normBC: ---------------------------------------------------------------------------------------------------'
       write(stdout,'(1x,a,5(1p,3x,e15.9),/)') 'e_normBC: ',  KEU, KEV,               APET,               APEP,                 AMEQ
  endif

contains

  real function mysqrt (a,myim)
     real a
     real, optional :: myim
     mysqrt = a
     if (opt=='SQ') then
         if(present(myim)) then
            mysqrt = sqrt(a*myim)
         else
            mysqrt = sqrt(a) 
         endif
     endif
  end function mysqrt


end subroutine e_normBC_




!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: E_NormD_ --- Total Energy or Kinetic Energy Norm: xpert to z transform
!
! !INTERFACE:

subroutine e_normD_ ( comm, ROOT, mynorm, x, xpert, nzvecsize, zvector, stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:
  
  integer, intent (in) :: ROOT
  integer, intent (in) :: comm
  character(len=*), intent(in) :: mynorm     ! name of norm
  type(dyn_prog), intent(in)  :: x           ! reference state stored as vector
  integer,  intent(in)  :: nzvecsize         ! size of zvector
  type(dyn_prog), intent(in)  :: xpert       ! adjoint variable corresponding to x

! !INPUT/OUTPUT PARAMETERS:


! !OUTPUT PARAMETERS:

  real(r8), intent(out) :: zvector(nzvecsize) ! z Lanczos vector
  integer,  intent(out), optional ::  stat    ! status indicator 

! !DESCRIPTION: This routine implements step D of the total energy norm.
!               It is the adjoint of the routine e\_normA\_
!               It copys xpert to z after scaling xpert by the inverse of the 
!               sqrt of the
!               energy coefficients.  The energy is defined following the
!               routine benergy except for a corrected grid dimension 
!               normalization constant and replacing PE by perturbation APE
!               This adjoint routine need not operate backwards compared to
!               step A since the operator is diagonal, except for relation
!               between delp and ps.  Also doubles as KE norm operator. 
!
!
! !REVISION HISTORY:
!
!  24Nov2002  Todling    Setting routine's template.
!  22Jan2003  Errico     Designing working routine
!  17Mar2003  Todling    Generalized character usage.
!  15Apr2005  Elena N.   Changed te-norm calculation to count for p_ref
!                        Changed te-norm calculation not to count for q_ref (MOIST NORM vs DRY NORM difference)
!  16May2007  Todling    Introduced dyn_prog.
!  21May2007  Todling    Fix aux_q(norm=ke) memory leak.
!  05Jun2007  Kim        MPI Parallelization
!  13Dec2007  Todling    Add q-component to total energy
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::e_normD_'

  real(r8)   :: scalej                     ! not used
  real(r8)   :: scalei                     ! not used  
  real(r8)   :: pi                         ! pi
  real(r8)   :: dblejnp                    ! r8 version of jnp
  real(r8)   :: dbleimr                    ! r8 version of imr
  real(r8)   :: ptop                       ! specified pressure at model top
  real(r8)   :: pint                       ! not used here  
  real(r8)   :: polefac                    ! gridfac*polefac=area of polar cap
  real(r8)   :: gridfac                    ! area factor for all points
  real(r8)   :: scalefac                   ! scaling for z to x
  real(r8)   :: tfactor                    ! c_p/tref
  real(r8)   :: qfactor                    ! one
  real(r8)   :: pfactor                    ! r_dryair*tref
  real(r8), allocatable ::  sine(:)
  real(r8), allocatable ::  cosp(:)
  real(r8), allocatable ::  sinp(:)
  real(r8), allocatable ::  cose(:) 
  real(r8), allocatable ::  ak(:)
  real(r8), allocatable ::  bk(:)
  real(r8), allocatable ::  dsigma(:,:,:)  ! reference state (delta p)/ps
  real(r8), allocatable ::  pkz(:,:,:)     ! reference state (p/po)**kappa
  real(r8), allocatable ::  ps(:,:)        ! reference state ps (surface p)
  real(r8), allocatable ::  pspert(:,:)    ! perturbed ps
  real(r8), allocatable ::  aux_q(:,:,:)   ! work array q_ref to convert TV to T


  integer  i        ! longitudinal index
  integer  j        ! latitudinal index
  integer  k        ! vertical index
  integer  nlp1     ! number of vertical levels plus 1
  integer  ks       ! not used
  integer  ierr     ! error indicator
  integer  izcount  ! counter for elements of the z vector
  integer  myID     ! processor ID
  integer  js,je

  character(len=len(mynorm)) :: usrnorm

  call MP_comm_rank(comm,myID,ierr)
      if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

  usrnorm = lowercase ( mynorm )

  if(present(stat)) stat=0

! specify required horizontal grid area-weighting factors

  allocate( sine(jnp), cosp(jnp), sinp(jnp), cose(jnp), stat=ierr )
	if(ierr/=0) call MP_die(myname_,'Alloc(sine,...)',ierr)
  call setrig(imr, jnp, scalej, scalei, cosp, cose, sinp, sine)
  pi = 4.d0 * datan(1.d0)
  dblejnp = dble(jnp)
  dbleimr = dble(imr)
  polefac = dbleimr*(dblejnp-1)*(1.d0+sine(2))/pi
  gridfac = pi/(4.d0*dbleimr*(dblejnp-1.d0))
  deallocate( sine, cose, sinp, stat=ierr )
  if(ierr/=0) call MP_die(myname_,'Dealloc(sine, etc)',ierr)

! specify required vertical grid weighting constant factors

  nlp1=nl+1
  allocate( ak(nlp1), bk(nlp1), stat=ierr )
	if(ierr/=0) call MP_die(myname_,'Alloc(sine,...)',ierr)
  call set_eta(nl, ks, ptop, pint, ak, bk)

! specify grid weights that depend on reference pressure 

   if ( usrnorm == 'ke' ) then  
    allocate( dsigma(imr,jfirst-1:jlast,nl), stat=ierr )
  	if(ierr/=0) call MP_die(myname_,'Alloc(dsigma)',ierr)
    call p_for_norms_ ( x%delp, r_dryair, c_p, ptop, dsigma, ierr )  
	if(ierr/=0) call MP_die(myname_,'p_for_norms 1',ierr)
  else ! TE norm also need pkz and ps
    allocate( dsigma(imr,jfirst-1:jlast,nl), pkz(imr,jfirst:jlast,nl), ps(imr,jfirst:jlast), stat=ierr )
  	if(ierr/=0) call MP_die(myname_,'Alloc(sine, etc.)',ierr)
    call p_for_norms_ ( x%delp, r_dryair, c_p, ptop, dsigma, ierr, pkz=pkz, ps=ps )
	if(ierr/=0) call MP_die(myname_,'p_for_norms 2',ierr)
  endif 

! initialize the zvector counter and allocate auxillary array
 
  izcount=0
  allocate( aux_q(imr,jfirst:jlast,nl), stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Alloc(aux_q)',ierr)
!!
! map the u-component of pert to the model grid, scale it, then map to z
!!  


  aux_q(:,jfirst:jlast,:) = x%q(:,jfirst:jlast,:,1)

!   u points next to South Pole first
  if (jfirst == 1) then
  do k=1,nl     
    do i=1,imr 
        izcount=izcount+1
        scalefac=gridfac * ( 0.5d0 * dsigma(i,2,k) * cosp(2) + &
                                     dsigma(1,1,k) * polefac / dbleimr )
        zvector(izcount) = xpert%u(i,2,k) / sqrt(scalefac)
    enddo
  enddo
  endif

!   u points not adjacent to either pole
  js=jfirst
  if (jfirst == 1) js=3

  je=jlast
  if (jlast == jnp) je=jnp-1

  do j=js,je
    do k=1,nl
      do i=1,imr
        izcount=izcount+1
        scalefac=0.5d0 * gridfac * ( dsigma(i,  j,k) * cosp(  j) +  &
                                     dsigma(i,j-1,k) * cosp(j-1)   )
        zvector(izcount) = xpert%u(i,j,k) / sqrt(scalefac)
      enddo
    enddo
  enddo

!   u points next to North Pole last
if (jlast == jnp) then
  do k=1,nl     
    do i=1,imr 
        izcount=izcount+1
        scalefac=gridfac * ( 0.5d0 * dsigma(i,jnp-1,k) * cosp(jnp-1) + &
                                     dsigma(1,jnp,k) * polefac / dbleimr )
        zvector(izcount) = xpert%u(i,jnp,k) / sqrt(scalefac)
    enddo
  enddo
  endif

! end of mappings and calculations for u

!!
! map the v-component of xpert to the model grid, scale it, then map to z
!!  


!   v points next to South Pole first
  if (jfirst == 1) then
      do k=1,nl     
         scalefac=gridfac * ( dsigma(1,1,k) * polefac / dbleimr +  &
                  0.5d0 * cosp(2) * ( dsigma(1,2,k) + dsigma(imr,2,k) ) )
         zvector(izcount+1) = xpert%v(1,2,k) / sqrt(scalefac)
         do i=2,imr 
            scalefac=gridfac * ( dsigma(1,1,k) * polefac / dbleimr +  &
                    0.5d0 * cosp(2) * ( dsigma(i,2,k) + dsigma(i-1,2,k) ) )
            zvector(izcount+i) = xpert%v(i,2,k) / sqrt(scalefac)
         enddo
         izcount=izcount+imr
      enddo
  endif
!   v points not adjacent to either pole
  js=jfirst
  if (jfirst == 1) js=3

  je=jlast
  if (jlast == jnp) je=jnp-2

  do j=js,je
    do k=1,nl
      scalefac=gridfac *                                             &
              0.5d0 * cosp(j) * ( dsigma(1,j,k) + dsigma(imr,j,k) ) 
        zvector(izcount+1) = xpert%v(1,j,k)/ sqrt(scalefac)
      do i=2,imr
        scalefac=gridfac *                                           & 
              0.5d0 * cosp(j) * ( dsigma(i,j,k) + dsigma(i-1,j,k) ) 
        zvector(izcount+i) = xpert%v(i,j,k)/ sqrt(scalefac)
      enddo
      izcount=izcount+imr
    enddo
  enddo

!   v points next to North Pole last
  if (jlast == jnp) then
      do k=1,nl  
         scalefac=gridfac * ( dsigma(1,jnp,k) * polefac / dbleimr +  &
                  0.5d0 * cosp(jnp-1) * ( dsigma(1,jnp-1,k) + dsigma(imr,jnp-1,k) ) )
         zvector(izcount+1) = xpert%v(1,jnp-1,k) / sqrt(scalefac)
         do i=2,imr 
            scalefac=gridfac * ( dsigma(1,jnp,k) * polefac / dbleimr + 0.5d0 &
                    * cosp(jnp-1) * ( dsigma(i,jnp-1,k) + dsigma(i-1,jnp-1,k) ) )
            zvector(izcount+i) = xpert%v(i,jnp-1,k) / sqrt(scalefac)
         enddo
         izcount=izcount+imr
      enddo
  endif

! end of mappings and calculations for v

!!
! If KE norm, then finish up here
!!
   if ( usrnorm == 'ke' ) then
       if (izcount /= nzvecsize) then 
         call MP_die(myname_,'izcount /= izvecsize',ierr)
       endif
       deallocate( dsigma, ak, bk, cosp, aux_q, stat=ierr )
         if(ierr/=0) call MP_die(myname_,'Dealloc(dsigma, etc.)',ierr)
       return
  endif

!!
! map the pt-component of xpert (=adjoint output) to the model grid, scale it, 
! then map to the z vector.
! Assumes that the adjoint variables at the pole have not been accumulated
! (summed over i) into a single value that represents the actual single pole
! point.  The summing is done here.  The z vector assumes that the pole is 
! defined as a single point.
!!  


    tfactor = c_p / tref

!   theta points at the South Pole first
  if (jfirst == 1) then
      do k=1,nl     
         izcount=izcount+1
         scalefac=gridfac * tfactor * dsigma(1,1,k) * polefac &
                 * (pkz(1,1,k)/(1.+zvir*aux_q(1,1,k)))**2  
         zvector(izcount)=0.d0
         do i=1,imr
            zvector(izcount) = zvector(izcount) + xpert%pt(i,1,k) / sqrt(scalefac)
         enddo
      enddo
  endif

!   theta points not at either pole
  js=jfirst
  if (jfirst == 1 ) js=2

  je=jlast
  if (jlast == jnp ) je=jnp-1

  do j=js,je
    do k=1,nl
      do i=1,imr
        izcount=izcount+1
        scalefac=gridfac * tfactor * dsigma(i,j,k) * cosp(j) &
                * (pkz(i,j,k)/(1.+zvir*aux_q(i,j,k)))**2 
        zvector(izcount) = xpert%pt(i,j,k) / sqrt(scalefac)
      enddo
    enddo
  enddo

!   theta points at the North Pole last
  if (jlast == jnp) then
      do k=1,nl     
         izcount=izcount+1
         scalefac=gridfac * tfactor * dsigma(1,jnp,k) * polefac &
                 * (pkz(1,jnp,k)/(1.+zvir*aux_q(1,jnp,k)))**2 
         zvector(izcount)=0.d0
         do i=1,imr
            zvector(izcount) = zvector(izcount) + xpert%pt(i,jnp,k) / sqrt(scalefac)
         enddo
      enddo
  endif

! end of mappings and calculations for pt

! If so, acccount for q-component in total energy

  if ( usrnorm == 'we' ) then

    qfactor = eps_eer * alhl * alhl / (c_p * tref)

!   theta points at the South Pole first
  if (jfirst == 1) then
      do k=1,nl     
         izcount=izcount+1
         scalefac=gridfac * qfactor * dsigma(1,1,k) * polefac
         zvector(izcount)=0.d0
         do i=1,imr
            zvector(izcount) = zvector(izcount) + xpert%q(i,1,k,1) / sqrt(scalefac)
         enddo
      enddo
  endif

!   theta points not at either pole
  js=jfirst
  if (jfirst == 1 ) js=2

  je=jlast
  if (jlast == jnp ) je=jnp-1

  do j=js,je
    do k=1,nl
      do i=1,imr
        izcount=izcount+1
        scalefac=gridfac * qfactor * dsigma(i,j,k) * cosp(j)
        zvector(izcount) = xpert%q(i,j,k,1) / sqrt(scalefac)
      enddo
    enddo
  enddo

!   theta points at the North Pole last
  if (jlast == jnp) then
      do k=1,nl     
         izcount=izcount+1
         scalefac=gridfac * qfactor * dsigma(1,jnp,k) * polefac
         zvector(izcount)=0.d0
         do i=1,imr
            zvector(izcount) = zvector(izcount) + xpert%q(i,jnp,k,1) / sqrt(scalefac)
         enddo
      enddo
  endif

! end of mappings and calculations for q

  endif  ! < wet energy >

!!
! deallocate some reference-state arrays no longer needed
!!
  deallocate( aux_q, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(aux_q)',ierr)
  deallocate( pkz, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(pkz)',ierr)
  deallocate( dsigma, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(dsigma)',ierr)

!!
! map the delp-component of xpert to the model grid, compute ps, scale it, 
! then map to the ps component of z vector.
! Assumes that the adjoint variables at the pole have not been accumulated
! (summed over i) into a single value that represents the actual single pole
! point.  The summing is done here.  The z vector assumes that the pole is 
! defined as a single point.
!!  


  pfactor = r_dryair * tref 

  allocate( pspert(imr,jfirst:jlast), stat=ierr )
	if(ierr/=0) call MP_die(myname_,'Alloc(pspert)',ierr)

!  adjoint of the functional relationship between delp and ps

  pspert(:,jfirst:jlast) = 0.d0

  do j=jfirst,jlast
    do k=1,nl
      do i=1,imr
        pspert(i,j) = pspert(i,j) + xpert%delp(i,j,k) * ( bk(k+1) - bk(k) ) / &
                                                        ( bk(nl+1) - bk(1) )
      enddo
    enddo
  enddo

!   ps points at the South Pole first
  if (jfirst == 1) then
      izcount=izcount+1
      scalefac=gridfac * pfactor * polefac / pref**2
      zvector(izcount)=0.d0
      do i=1,imr
         zvector(izcount) = zvector(izcount) + pspert(i,1) / sqrt(scalefac)
      enddo
  endif

!   ps points not at either pole
  js=jfirst
  if (jfirst == 1) js=2

  je=jlast
  if (jlast == jnp) je=jnp-1

  do j=js,je
    scalefac=gridfac * pfactor * cosp(j) / pref**2
    do i=1,imr
      izcount=izcount+1
      zvector(izcount) = pspert(i,j) / sqrt(scalefac)
    enddo
  enddo

!   ps points at the North Pole last
  if (jlast == jnp) then
      izcount=izcount+1
      scalefac=gridfac * pfactor * polefac / pref**2
      zvector(izcount)=0.d0
      do i=1,imr
         zvector(izcount) = zvector(izcount) + pspert(i,jnp) / sqrt(scalefac)
      enddo
  endif

! end of mappings and calculations for ps

!!
!  deallocate remaining local arrays
!!

  deallocate( pspert, ps, ak, bk, cosp, stat=ierr )
    if(ierr/=0) call MP_die(myname_,'Dealloc(pspert, etc.)',ierr)

!
!! check that entire z vector is accounted for

  if (izcount /= nzvecsize) then 
       call MP_die(myname_,'izcount /= nzvecsize',ierr)
  endif

end subroutine e_normD_






!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: p_for_norms_ --- Compute delta sigma and p/po **kappa for norms
!
! !INTERFACE:

subroutine p_for_norms_ ( xdelp, r_dryair, c_p, ptop, dsigma, stat, &
                          pkz, ps )   ! optionals

! !USES:

  implicit none


! !INPUT PARAMETERS:

  real(r8), intent(in) :: xdelp(imr,jfirst:jlast,nl)  ! complete refererence state vector
  real(r8), intent(in) :: r_dryair      ! gas constant for dry air
  real(r8), intent(in) :: c_p           ! specific heat of dry air
  real(r8), intent(in) :: ptop          ! pressure at model top

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:

  real,    intent(out)             ::  dsigma(imr,jfirst-1:jlast,nl)
  integer, intent(out)             ::  stat
  real,    intent(out), optional   ::  pkz(imr,jfirst:jlast,nl)
  real,    intent(out), optional   ::  ps(imr,jfirst:jlast)

! !DESCRIPTION: This routine determines delta sigma and/or (p/po)**kappa
!               and/or p\_surface from the reference state vectors x or y 
!               for use in construction of coefficients used to define 
!               norms in the SV software. It is intended to be called by
!               a norm scaling or computing  routine such as te\_norm\_  
!               It essentially duplicates algorithms found in various 
!               places in the FVGCM.  
!
! To do:
!  \begin{itemize}
!      \item MPI part
!      \item double check this is really correct
!  \end{itemize}
!
! !REVISION HISTORY:
!
!  22Jan2003 Errico   initial construction
!  16May2007 Todling  Introduced dyn_prog.
!  22May2007 Todling  This prog was calling pk what really is pkz
!  05Jun2007 KIM      MPI Parallelization
!
!EOP
!-------------------------------------------------------------------------

  real, parameter :: p00=1.  ! reference p where theta=T (Note = 1 Pa !!!)
  character(len=*), parameter :: myname_ = myname//'::p_for_norms_'

  real(r8), allocatable ::  pe(:,:,:)   ! p at interfaces
  real(r8), allocatable ::  pke(:,:)    ! (p/po)**kappa at interfaces
  real(r8), allocatable ::  logpe(:,:)  ! log(p) at interfaces
  real(r8)              ::  cappa       ! R/C_p
  real(r8)              ::  ptopk       ! ptop^cappa
 
  logical  failed
  integer  nlp1
  integer  k, i, j
  integer  ierr

  stat=0

  cappa=kappa

  nlp1=nl+1

! Consistency checks
! ------------------
  failed = .false.
  if(size(dsigma,1)/=imr           ) failed = .true.
  if(size(dsigma,2)/=jlast-jfirst+2) failed = .true.
  if(size(dsigma,3)/=nl            ) failed = .true.
    if(failed)then
       stat = 91
       call mpout_log(myname_,' dimensions(dsigma)')
       return
    endif

  if(size(xdelp,1)/=imr           ) failed = .true.
  if(size(xdelp,2)/=jlast-jfirst+1) failed = .true.
  if(size(xdelp,3)/=nl            ) failed = .true.
    if(failed)then
       stat = 92
       call mpout_log(myname_,' dimensions(xdelp)')
       return
    endif
  
! compute p at interfaces
! assume here that array is complete (no special pole treatment)
! NOTE: pe is defined here transpose to what the model defines it

  allocate( pe(imr,jfirst:jlast,nlp1), stat=ierr )
    if(ierr/=0)then
       stat=93
       call mpout_log(myname_,'Alloc(pe)')
       return
    endif
 
   do j=jfirst,jlast
      do i=1,imr
         pe(i,j,1) = ptop
      enddo
      do k=2,nl+1
         do i=1,imr
            pe(i,j,k) = pe(i,j,k-1) + xdelp(i,j,k-1)
         enddo
      enddo
   enddo

! compute dsigma if requested where pe(:,:,nlp1)=ps

  do k=1,nl
     dsigma(:,jfirst:jlast,k) = xdelp(:,jfirst:jlast,k) &
                              / (pe(:,jfirst:jlast,nlp1) - ptop)
  enddo

#if defined( SPMD )
!-------------------------------------
! Send/recv dsigma from south to north 
!-------------------------------------
      call mp_send_n(imr, jnp, jfirst, jlast, 1, nl, 1, 0, dsigma)
      call mp_barrier
      call mp_recv_s(imr, jnp, jfirst, jlast, 1, nl, 1, 0, dsigma)
#endif

! copy ps for output if requested 

  if(present(ps)) then  
     if(size(ps,1)/=imr           ) failed = .true.
     if(size(ps,2)/=jlast-jfirst+1) failed = .true.
       if(failed)then
          stat = 94
          call mpout_log(myname_,' dimensions(ps)')
          return
       endif
     ps(:,jfirst:jlast)= pe(:,jfirst:jlast,nlp1)
  endif

! compute (p/po)**kappa if requested

  if(present(pkz)) then  

    if(size(pkz,1)/=imr           ) failed = .true.
    if(size(pkz,2)/=jlast-jfirst+1) failed = .true.
    if(size(pkz,3)/=nl            ) failed = .true.
      if(failed)then
         stat = 95
         call mpout_log(myname_,' dimensions(pkz)')
         return
      endif
  
    allocate( pke(imr,nlp1), logpe(imr,nlp1), stat=ierr )
      if(ierr/=0)then
         stat=96
         call mpout_log(myname_,'alloc(pke,logpe)')
         return
      endif
    
    do j=jfirst,jlast
      pke(:,:)=(pe(:,j,:)/p00)**cappa
      logpe(:,:)=log(pe(:,j,:))
      do k=1,nl
        pkz(:,j,k)=(pke(:,k+1)-pke(:,k)) &
                  / ( cappa * (logpe(:,k+1)-logpe(:,k)) )  
      enddo
    enddo

    deallocate( pke, logpe, stat=ierr )
      if(ierr/=0)then
         stat=97
         call mpout_log(myname_,'dealloc(pke,logpe)')
         return
      endif

  endif   ! test on present(pkz)

  deallocate( pe, stat=ierr )
      if(ierr/=0)then
         stat=98
         call mpout_log(myname_,'dealloc(pe)')
         return
      endif

end subroutine p_for_norms_

end module m_svnorms

