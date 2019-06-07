!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_admtlm --- Apply Oseledec's operator
!
! !INTERFACE:

  module m_admtlm

! !USES:

  use precision

  use prognostics, only : prognostics_initial
  use prognostics, only : prognostics_final
  use prognostics, only : dyn_prog

  use stepon,  only : stepon_set
  use stepon,  only : nymd
  use stepon,  only : nhms

  use m_model_ad, only : amodel_ad
  use m_model_tl, only : amodel_tl

  use m_sVecDef, only : test_norm  ! Determines test levs for norms
  use m_sVecDef, only : svecnormI  ! Initial state norm
  use m_sVecDef, only : svecnormF  ! Final state norm

  use m_svnorms, only : norm_svecA
  use m_svnorms, only : norm_svecBC
  use m_svnorms, only : norm_svecD
  use m_svnorms, only : norm_zz
  use m_svprojs, only : proj_svec

  use m_mpif90,only : MP_comm_rank

  use m_die, only : die
  use m_die, only: MP_die

  use timingModule
  use m_stdio, only : stdout


  implicit NONE
 
! !PUBLIC MEMBER FUNCTIONS:
 
  PRIVATE

  PUBLIC admtlm_init
  PUBLIC admtlm_set
  PUBLIC admtlm
  PUBLIC admtlm_clean

  interface admtlm_init ; module procedure init_  ; end interface
  interface admtlm_set  ; module procedure set_   ; end interface
  interface admtlm      ; module procedure &
                          admtlm0_, &
                          admtlm1_
  end interface
  interface admtlm_clean; module procedure clean_ ; end interface

  integer,  save                            :: nn  ! dim of input  state vector
  integer,  save                            :: mm  ! dim of output state vector
  integer,  save                            :: nnz ! dim of pert  vector
  real(r8), save, allocatable, dimension(:) :: zzi ! input  perturbation array
  real(r8), save, allocatable, dimension(:) :: zzo ! output perturbation array

  type(dyn_prog), save :: xx                       ! input  state vector


!
! !DESCRIPTION: Apply Oseledec's operator
!
! !REVISION HISTORY:
!
!  26Apr2004  Todling   Modularized.
!  14May2007  Todling   Introduced dyn_prog; global change.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'm_admtlm'
  logical,          parameter :: debug  = .false.

  CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: admtlm0_ --- Apply Oseledec's operator to perturbation vector
!
! !INTERFACE:

subroutine admtlm0_( comm, ROOT, &
                     nzvecsize, zvector_in, zvector_out,  &
                     nymd_b, nhms_b, nymd_e, nhms_e )  !optionals

  implicit none

! !INPUT PARAMETERS:

  integer, intent(in) :: ROOT
  integer, intent(in) :: comm
  integer, intent(in) :: nzvecsize             ! size of z vector
  real(r8),intent(in) :: zvector_in(nzvecsize) ! z used to create initial pert


! !OUTPUT PARAMETERS:

  integer,  intent(out), optional  :: nymd_b
  integer,  intent(out), optional  :: nhms_b
  integer,  intent(out), optional  :: nymd_e
  integer,  intent(out), optional  :: nhms_e
  real(r8), intent(out)            :: zvector_out(nzvecsize) ! result of this operator*z

! !DESCRIPTION:  This routine calculates the result of applying the operator
!  \begin{equation*}
!    O=S_1^{T~-1/2} L_M^T P^T S_2^{T~1/2} S_2^{1/2} P L_M S_1^{-1/2}
!  \end{equation*}
! to the input vector $z=zvector_in$ from the vector-space considered by the 
! Lanczos algorithm, where
! \begin{enumerate}
!  \item $S_1^{-1/2}$ defines the transform from $z$ to pert $x$.  The
!   constraining applied to the initial time is $x^T~S_1~x=z^T~z=1$
!  \item $L_M$ is the tangent linear model with the trajectory initialized
!    by x and perturbation input pert
!  \item $P$ is the local projection operator that specifies a geographical
!    region over which to calculate the end-time norm
!  \item $S_2$ defines the end-time norm to be applied to the locally 
!    projected TLM output (called $g_y$ here).  It also creates the vector 
!    $g_y = S_2 ypert$ for input into the adjoint routines. 
!  \item superscript $T$ denotes the usual adjoint, defined with respect 
!    to the identity norm in the space of the $x$ or $y$ vectors.
! \end{enumerate}
!
!  Tests:
!     test\_norm=0  : perform no tests of norm
!     test\_norm=1  : compare ynorm $=y^T~S_2~y$ with znorm $=z^T~O~z$
!                    These should be identical except for roundoff. Values of
!                    diff\_norm=ynorm-znorm that are 10*(-10) smaller than
!                    ynorm are acceptable
!     test\_norm=2  : like test\_norm=1, except ADM and TLM replaced by identity
!
! !REVISION HISTORY:
!
!  15Oct2002 FastOpt - Initial code.
!  25Oct2002 Todling - Added prologue; added hook for norms
!  20Nov2000 Todling - Added call to local projection operator
!  24Jan2003 Errico  - Adding corrrect TE norm and checking of operators
!  11Mar2003 Todling - Edited prologue and comments
!  16May2007 Todling - Introduced dyn_prog
!
!EOP
!-------------------------------------------------------------------------
!BOC

  character(len=*), parameter :: myname_ = myname//'::admtlm0_'

  type(dyn_prog) :: xpert    ! local perturbations

  real(r8)   ynorm           ! value of norm $S_2$
  real(r8)   znorm           ! value of norm $z^T~Oz$ where $O$ is entire operator 
  real(r8)   diff_norm       ! difference in value of norms z/ynorm (should=0)
  integer    myid,ierr

  call MP_comm_rank(comm,myID,ierr)
      if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)

! For debugging purposes, make matrix equals identity
! ---------------------------------------------------
  if ( debug ) then
       zvector_out = zvector_in
       return
  endif

  call prognostics_initial (  xpert )

! Apply the operator S^{-1/2} to zvector_in
! -----------------------------------------
	call timing_on ( 'normA' )

  call norm_svecA  ( comm, ROOT, svecnormI, xx, xpert, nzvecsize, zvector_in, ierr )
     if (ierr/=0) call die(myname_,'Error from norm_svecA',ierr)
	call timing_off ( 'normA' )

! Apply the TLM operator L_M to xpert to create ypert
! ---------------------------------------------------
  if (test_norm < 2) then
    	call timing_on ( 'func_TL' )
        call amodel_tl ( xpert )
    	call timing_off ( 'func_TL' )
    if(present(nymd_e)) nymd_e = nymd	! Tap on final date/time
    if(present(nhms_e)) nhms_e = nhms
  endif
 
! Apply local projection operator P to ypert
! ------------------------------------------
 	call timing_on ( 'ProjSvec' )
  call proj_svec ( comm, ROOT, xpert, ierr )   ! Apply local projection operator
        if(ierr/=0) call die(myname_,'Error from proj_svec',ierr)
 	call timing_off ( 'ProjSvec' )
 
! Apply the operator S^2=S^(T/2)S^{1/2)
! -------------------------------------
	call timing_on ( 'normBC' )
  call norm_svecBC ( comm, ROOT, svecnormF, xx, xpert, ynorm, 'BC', ierr )
	if(ierr/=0) call die(myname_,'Error from norm_svecBC',ierr)
	call timing_off ( 'normBC' )
 
! Apply the operator P^T=P (assumed such that P^T=P here)
! -------------------------------------------------------
  	call timing_on ( 'ProjSvec' )
  call proj_svec ( comm, ROOT, xpert, ierr )   ! Apply local projection operator
	if(ierr/=0) call die(myname_,'Error from proj_svec',ierr)
  	call timing_off ( 'ProjSvec' )
 
! Apply the ADM operator L_M^T to ypert adjoint to create xpert adjoint
! ---------------------------------------------------------------------
  if (test_norm < 2) then
     call timing_on ( 'func_AD' )
         call amodel_ad( xpert )
     call timing_off ( 'func_AD' )
     if(present(nymd_b)) nymd_b = nymd	! Tap on initial date/time
     if(present(nhms_b)) nhms_b = nhms
  endif
 
! Apply the operator S^{-T/2} to create zvector_out
! -------------------------------------------------
	call timing_on ( 'normD' )
  call norm_svecD ( comm, ROOT, svecnormI, xx, xpert, nzvecsize, zvector_out, ierr)
	if(ierr/=0) call die(myname_,'Error from norm_svecD',ierr)
	call timing_off ( 'normD' )

! If testing, we are done
! -----------------------
  if (test_norm == 0) return  

! Apply norm
! ----------
  call norm_zz (comm, ROOT, nzvecsize, zvector_in, zvector_out, znorm, ierr)


  diff_norm = ynorm - znorm
  if(myid==root) write(stdout,*) myname_, ': y,z,diff norms =',ynorm, znorm, diff_norm

  call prognostics_final (  xpert )

end subroutine admtlm0_
!EOC

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_ --- Initialize internal state vector
!
! !INTERFACE:

subroutine init_ ( n, m, nz )
  implicit none

! !INPUT PARAMETERS:

  integer, intent(in)  :: n   ! dim of initial state
  integer, intent(in)  :: m   ! dim of final state
  integer, intent(in)  :: nz  ! dim of perturbation

! !DESCRIPTION:  Initilization of internal state vector.
!
! !REVISION HISTORY: 
!
!   26Mar2004  Todling    Initial code.
!   16May2007  Todling    Updated to handle dyn_prog-like vectors 
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::init_'

  integer  i, ier

  if(m/=n) call die (myname_,'dim error')
  nn  = n
  mm  = m
  nnz = nz
  allocate ( zzi(nz), zzo(nz), stat=ier )
    if(ier/=0) call die (myname_,'Alloc(zz)')

  call prognostics_initial ( xx )

end subroutine init_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_ --- Set input arrays in state vector
!
! !INTERFACE:

subroutine set_ ( x )

  implicit none

! !INPUT PARAMETERS:

  type(dyn_prog), intent(in) :: x

! !DESCRIPTION:  Sets up input state vector.
!
! !REVISION HISTORY:
!
!   26Mar2004  Todling    Initial code.
!   16May2007  Todling    Updated to handle dyn_prog-like vectors 
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::set_'

  integer  ier

  call stepon_set ( xx )

end subroutine set_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: admtlm1_ --- Alternative interface to Oseledec's operator
!
! !INTERFACE:

subroutine admtlm1_ ( z, nz )

! !USES:
                                                                                        
  use m_sVecDef, only : ncalls        ! counter for lanczos iterations
  use m_sVecDef, only : lanunit       ! unit for lanczos restart data

  use m_mpif90, only : comm => MP_comm_world

  use m_ioutil, only : luavail

  use timingModule

  implicit none

! !INPUT PARAMETERS: 
 
  integer, intent(in)    :: nz      ! dimension of z

! !INPUT/OUTPUT PARAMETERS:
 
  real(r8), intent(out)  :: z(nz)   ! result of this operator*z

! !OUTPUT PARAMETERS:

! !DESCRIPTION:  Interface to admtlm for both ARPACK and NAG-lib solver packages. 
!                Also controls restart capability for Lanczos iterations for 
!                either package. To use this interface for admtlm, the following 
!                should be the sequence of calls. See admtlm\_ for details.
!
!       call admtlm\_init ( n, m, nz )
!       call admtlm\_set  ( x, z )
!       call admtlm       ( z, nz [,m, y] )
!       call admtlm\_clean( )
!                
!
! !REVISION HISTORY: 
!
!   26Mar2004  Todling    Initial code.
!   21Dec2005  Gelaro     Modified for Lanczos restart capability
!   16May2007  Todling    Updated to handle dyn_prog-like vectors
!   09Jun2007  Kim        Added local Lanczos restart flag: iflag_restart
!   23Aug2007  Todling    Renamed routine (overload admtlm)
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::admtlm1_'

  integer  :: ROOT = 0 ! this is bad, but I cannot pass ROOT as argument
  integer  :: myid,ierr
  real(r4) :: z4(nz)
  integer  :: iflag_restart
  character(len=*), parameter :: fname = 'lancdat.bin'

  if(nz/=nnz) call die (myname_,'dim error(nz)')

  call MP_comm_rank(comm,myID,ierr)
      if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

! Open/create file with Lanczos restart data
! ------------------------------------------
  if ( ncalls == 0 ) then
    lanunit = luavail()
    open ( unit=lanunit, file=trim(fname), form='unformatted', access='direct', &
        recl=4*nz,status='old',iostat=ierr )
    if (ierr == 0) iflag_restart=1
    if (ierr /= 0) then
        lanunit = luavail()
        open ( unit=lanunit, file=fname, form='unformatted', access='direct', &
               recl=4*nz,status='new',iostat=ierr )
        iflag_restart=0
    endif
  end if


! If restart data exist, read and return, else call adm_tlm
! ---------------------------------------------------------
  ncalls = ncalls + 1
  if(myid==root) print *, 'ncalls= ', ncalls

  if (iflag_restart == 1) then
      read ( lanunit, rec=ncalls, iostat=ierr ) z4   ! data saved at r4
      if ( ierr == 0 ) then
           print *, 'read zvector for ncalls= ', ncalls
           z = z4
           return
      endif
  end if

  zzi = z

! At this point x and zi must have been set
! -----------------------------------------
	call timing_on ( 'AdmTlm' )
  call admtlm0_( comm, ROOT, nz, zzi, zzo )
	call timing_off ( 'AdmTlm' )

  z = zzo

! Write vector to restart data file for this iteration
! ----------------------------------------------------
  z4 = z
  write ( lanunit, rec=ncalls, iostat=ierr ) z4
  if(myid==root .and. ierr==0) print *, 'wrote zvector for ncalls= ', ncalls

end subroutine admtlm1_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: clean_ --- Deallocate arrays memory
!
! !INTERFACE:

subroutine clean_ ( )
  implicit none

! !INPUT PARAMETERS:

! !DESCRIPTION:  Deallocate state vector
!
! !REVISION HISTORY:
!
!   26Mar2004  Todling    Initial code.
!   16May2007  Todling    Updated to handle dyn_prog-like vectors 
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::clean_'

  integer  i, ier

  call prognostics_final ( xx )

  deallocate ( zzi, zzo, stat=ier )
    if(ier/=0) call die (myname_,'Dealloc(xx)')

end subroutine clean_

end module m_admtlm
