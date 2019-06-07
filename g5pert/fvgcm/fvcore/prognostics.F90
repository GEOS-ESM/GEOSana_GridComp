
module prognostics
!BOP
!
! !MODULE: prognostics --- dynamics-physics coupling module
!
   use precision
   use m_die, only : die
   use mod_comm, only : gid,mp_add1d
   use mod_comm, only : imr, jnp, nl, nc
   use m_random, only :  zufalli, normalen

   implicit none

! !PUBLIC MEMBER FUNCTIONS:

!  PRIVATE
   SAVE

   public prognostics_allocate
   public prognostics_deallocate
   public prognostics_initial
   public prognostics_zero
   public prognostics_dotp
   public prognostics_dup
   public prognostics_scal
   public prognostics_axpy
   public prognostics_cnst
   public prognostics_rand

   public dyn_prog

   type dyn_prog
!
! Ultimately the fields will be allocatable arrays, e.g., 
!
      real(r8), pointer ::    u(:,:,:)   ! zonal wind on D-grid
      real(r8), pointer ::    v(:,:,:)   ! meridional wind
      real(r8), pointer ::   pt(:,:,:)   ! virtual potential temperature
      real(r8), pointer :: delp(:,:,:)   ! pressure thickness (pascal)
      real(r8), pointer ::    q(:,:,:,:) ! specific humidity & tracer mixing ratios

   end type dyn_prog

! Dependent variables (everything that can be derived from prog. variables)

      real(r8), save, allocatable ::   pe(:,:,:)   ! pressure at layer edges
      real(r8), save, allocatable ::   pk(:,:,:)   ! pe**kappa
      real(r8), save, allocatable ::  pkz(:,:,:)   ! layer-mean value of pk
      real(r8), save, allocatable ::   ps(:,:)     ! surface pressure (pascal)
      real(r8), save, allocatable :: peln(:,:,:)   ! log pressure (pe) at layer edges

      real(r8), save, allocatable :: phis(:,:)     ! surface geopotential
      real(r8), save, allocatable ::  sgh(:,:)     ! std dev. of topg.

      integer  jfirst,jlast
      integer  kfirst,klast
      integer  klastp
      integer  iord
      integer  jord
      integer  kord
      integer  nq                            ! Total number of tracers

      integer :: ng_d_save
      integer :: ng_s_save

! !DESCRIPTION:
!
!   {\bf Purpose:} Prognostic variables held in-core for convenient 
!   access. q3 is specific humidity (water vapor) and other 
!   constituents. pcnst is advected constituents, pnats is non-advected.
! 
! !REVISION HISTORY:
!
!   00.08.15     Lin        Modifications
!   00.12.14     Sawyer     SPMD bug: do j=1,plat => beglat,endlat
!   01.03.26     Sawyer     Added ProTeX documentation
!   12May2007    Todling    Introduced dyn_prog data type
!   21Nov2007    Todling    Dynamic memory code
!
!EOP
!-----------------------------------------------------------------------

interface prognostics_initial; module procedure &
          prognostics_initial_, &
          prognostics_initial1_ 
end interface

interface prognostics_rand; module procedure &
          prognostics_rand_
end interface

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_allocate --- allocate prognostic variables
!
! !INTERFACE:
subroutine prognostics_allocate ( prog )

implicit none

  type(dyn_prog) :: prog

! !DESCRIPTION:
!
!   Initialize the prognostic variables
!
! !REVISION HISTORY:
!   00.10.15  Lin      modified to declare only one-time-level; 
!                      may want to declare them as 3D arrays eventually
!   01.03.26  Sawyer   Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

   character(len=*), parameter :: myname_ = 'prognostics_allocate_'
   integer :: i, j, k, m, ierr
   integer :: klastp

      klastp = klast
      if (klast .eq. nl) klastp = klast + 1


!
! Allocate variables
!
      call prognostics_initial_ ( prog )

      allocate(  pkz(imr,jfirst:jlast,kfirst:klast), stat=ierr )          ! layer-mean value of pk
         if(ierr/=0) call die(myname_,'Alloc(pkz)',ierr)
      allocate(  pk(imr,jfirst:jlast,kfirst:klastp), stat=ierr )         ! pe**kappa
         if(ierr/=0) call die(myname_,'Alloc(pk)',ierr)
      allocate(  pe(imr,kfirst:klastp,jfirst:jlast), stat=ierr )         ! pressure at layer edges
         if(ierr/=0) call die(myname_,'Alloc(pe)',ierr)
      allocate(  peln(imr,kfirst:klastp,jfirst:jlast) )
         if(ierr/=0) call die(myname_,'Alloc(peln)',ierr)
      allocate(  ps(imr,jfirst:jlast), stat=ierr )
         if(ierr/=0) call die(myname_,'Alloc(ps)',ierr)

      allocate( phis(imr,jfirst:jlast) )
      allocate(  sgh(imr,jfirst:jlast) )

#if defined (ALT_PBL)
      allocate(  tke(imr,jfirst:jlast,nl) )                    ! turbulent kinetic energy
#endif

!EOC
 end subroutine prognostics_allocate
!-----------------------------------------------------------------------

subroutine prognostics_initial_ ( prog )
implicit none
type(dyn_prog) :: prog

   character(len=*), parameter :: myname_ = 'prognostics_initial_'
   integer :: ng_s, ng_d, ierr

#if defined( SPMD )
        ng_d = max(2,min(abs(jord),3)) 
        ng_s = 3
#else
        ng_d = 0
        ng_s = 0
#endif
      ng_d_save = ng_d
      ng_s_save = ng_s
      klastp = klast
      nq = nc
      if (klast .eq. nl) klastp = klast + 1
!
! Allocate variables
!
      allocate(  prog%u(imr,jfirst-ng_d:jlast+ng_s,kfirst:klast), stat=ierr )
         if(ierr/=0) call die(myname_,'Alloc(u)',ierr)
      allocate(  prog%v(imr,jfirst-ng_s:jlast+ng_d,kfirst:klast), stat=ierr )
         if(ierr/=0) call die(myname_,'Alloc(v)',ierr)
      allocate( prog%pt(imr,jfirst-ng_d:jlast+ng_d,kfirst:klast), stat=ierr )  ! virtual potential temperature
         if(ierr/=0) call die(myname_,'Alloc(pt)',ierr)
      allocate( prog%delp(imr,jfirst:jlast,kfirst:klast), stat=ierr )          ! pressure thickness (pascal)
         if(ierr/=0) call die(myname_,'Alloc(delp)',ierr)
      allocate( prog%q(imr,jfirst-ng_d:jlast+ng_d,kfirst:klast,nc), stat=ierr )! specific humidity & tracer
         if(ierr/=0) call die(myname_,'Alloc(q)',ierr)

      call prognostics_zero ( prog )


end subroutine prognostics_initial_

subroutine prognostics_initial1_ ( prog, imr_, jnp_, lm_, nc_ )
implicit none
! 16May 2007 Todling This interface added in support on single-pe mpi applications
type(dyn_prog) :: prog
                                                                                                                       
   character(len=*), parameter :: myname_ = 'prognostics_initial1_'
   integer, intent(in) :: imr_, jnp_, lm_, nc_

   integer :: ierr
                                                                                                                       
!
! Allocate variables
!
      allocate(  prog%u(imr_,jnp_,lm_), stat=ierr )
         if(ierr/=0) call die(myname_,'Alloc(u)',ierr)
      allocate(  prog%v(imr_,jnp_,lm_), stat=ierr )
         if(ierr/=0) call die(myname_,'Alloc(v)',ierr)
      allocate( prog%pt(imr_,jnp_,lm_), stat=ierr )            ! virtual potential temperature
         if(ierr/=0) call die(myname_,'Alloc(pt)',ierr)
      allocate( prog%delp(imr_,jnp_,lm_), stat=ierr )          ! pressure thickness (pascal)
         if(ierr/=0) call die(myname_,'Alloc(delp)',ierr)
      allocate( prog%q(imr_,jnp_,lm_,nc_), stat=ierr )         ! specific humidity & tracer
         if(ierr/=0) call die(myname_,'Alloc(q)',ierr)

      call prognostics_zero ( prog )
                                                                                                                       
                                                                                                                       
end subroutine prognostics_initial1_


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_deallocate --- deallocate all prognostic variables
!
! !INTERFACE:
subroutine prognostics_deallocate (prog)      ! FastOpt
  implicit none
  type(dyn_prog) :: prog

! !DESCRIPTION: deallocate all prognostic variables
!
! !REVISION HISTORY:
!   01.08.21  Ralf Giering, FastOpt : initial version
!
!EOP
!-----------------------------------------------------------------------
!BOC
   character(len=*), parameter :: myname_ = 'prognostics_deallocate'
   integer  ierr

      call prognostics_final ( prog )

      deallocate( pkz  )
      deallocate( pk   )
      deallocate( pe   )
      deallocate( peln )
      deallocate( ps   )

      deallocate( phis, stat=ierr )
         if(ierr/=0) call die(myname_,'Dealloc(phis)',ierr)
      deallocate(  sgh, stat=ierr )
         if(ierr/=0) call die(myname_,'Dealloc(sgh)',ierr)


#if defined (ALT_PBL)
      allocate(  tke, stat=ierr )     ! turbulent kinetic energy
         if(ierr/=0) call die(myname_,'Dealloc(tke)',ierr)
#endif

!EOC
 end subroutine prognostics_deallocate       ! FastOpt
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_final --- deallocate prognostic variables
!
! !INTERFACE:
subroutine prognostics_final ( prog )
implicit none
! !INPUT PARAMETERS:

type(dyn_prog) :: prog

! !DESCRIPTION: deallocate prognostic variables
!
! !REVISION HISTORY:
!   12May2007 Todling Initial code
!
!EOP
!-----------------------------------------------------------------------
!BOC
   character(len=*), parameter :: myname_ = 'prognostics_final'
   integer  ierr

   if ( associated(prog%u)     ) deallocate(prog%u)
   if ( associated(prog%v)     ) deallocate(prog%v)
   if ( associated(prog%pt)    ) deallocate(prog%pt)
   if ( associated(prog%delp)  ) deallocate(prog%delp)
   if ( associated(prog%q)     ) deallocate(prog%q)

      nullify( prog%u    )
      nullify( prog%v    )
      nullify( prog%pt   )
      nullify( prog%delp )
      nullify( prog%q    )


end subroutine prognostics_final

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_zero --- zero out prognostic variables
!
! !INTERFACE:
subroutine prognostics_zero ( prog )
      implicit none
! !INPUT/OUTPUT PARAMETERS:
      type(dyn_prog) :: prog
! !DESCRIPTION: zero out prognostic variables
!
! !REVISION HISTORY:
!   12May2007 Todling Initial code
!
!EOP
!-----------------------------------------------------------------------
!BOC

      prog%u    = 0.0d0
      prog%v    = 0.0d0
      prog%pt   = 0.0d0
      prog%delp = 0.0d0
      prog%q    = 0.0d0

end subroutine prognostics_zero 

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_scal --- scale dyn_prog by a constant
!
! !INTERFACE:
subroutine prognostics_scal ( a, prog )
      implicit none
! !INPUT PARAMETERS:
      real(r8), intent(in) :: a
! !INPUT/OUTPUT PARAMETERS:
      type(dyn_prog) :: prog
! !DESCRIPTION: zero out prognostic variables
!
! !REVISION HISTORY:
!   12May2007 Todling Initial code
!
!EOP
!-----------------------------------------------------------------------
!BOC

      prog%u    = a * prog%u
      prog%v    = a * prog%v
      prog%pt   = a * prog%pt
      prog%delp = a * prog%delp
      prog%q    = a * prog%q

end subroutine prognostics_scal

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_cnst --- sets vector to a constant value
!
! !INTERFACE:
subroutine prognostics_cnst ( a, prog )
      implicit none
! !INPUT PARAMETERS:
      real(r8), intent(in) :: a
! !INPUT/OUTPUT PARAMETERS:
      type(dyn_prog) :: prog
! !DESCRIPTION: zero out prognostic variables
!
! !REVISION HISTORY:
!   12May2007 Todling Initial code
!
!EOP
!-----------------------------------------------------------------------
!BOC

      prog%u    = a
      prog%v    = a
      prog%pt   = a
      prog%delp = a
      prog%q    = a

end subroutine prognostics_cnst

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_rand_ --- random dyn_prog
!
! !INTERFACE:
subroutine prognostics_rand_ ( prog, seed )
      implicit none
! !INPUT PARAMETERS:
      integer, intent(in),optional :: seed
! !INPUT/OUTPUT PARAMETERS:
      type(dyn_prog) :: prog
! !DESCRIPTION: generate random dyn_prog vector
!
! !REVISION HISTORY:
!   14Aug2007 Todling Initial code.
!   31Jul2008 Oloso   Fix for proper distribution of random numbers.
!   08Dec2009 Todling Fix Oloso's fix!
!   20Jan2010 Todling Reset default: using m_random
!
!EOP
!-----------------------------------------------------------------------
!BOC

      character(len=*), parameter :: myname_ = 'prognostics_rand'
      integer seed_,n,isize,ierr,nn,j
      integer ng_d,ng_s
      real(r8),allocatable::aux(:)
      integer,allocatable:: iseed(:)
      logical  inter

      inter=.false.
      ng_d = ng_d_save
      ng_s = ng_s_save
      seed_ = 0
      if (present(seed)) then
          seed_ = seed
          if(seed_<0) seed_=0
          inter=.false.
      endif

      n = imr*(jlast-jfirst+1)*nl

      if (inter) then
          call random_seed(size=isize)
          allocate(iseed(isize))
          iseed(1)=seed_+gid
          call random_seed(put=iseed)
          deallocate(iseed)
      else
          call zufalli ( seed_+gid )
      endif

      j = (jlast+ng_s)-(jfirst-ng_d)+1
      n = imr*j*nl
      allocate(aux(n),stat=ierr)
         if(ierr/=0) call die(myname_,'Alloc(u-aux)',ierr)
      if (inter) then
          call random_number(aux)
      else
          call normalen(n,aux)
      endif
      prog%u(:,jfirst-ng_d:jlast+ng_s,:)=reshape(aux,(/imr,j,nl/))
      deallocate(aux,stat=ierr)
         if(ierr/=0) call die(myname_,'Dealloc(u-aux)',ierr)

      j = (jlast+ng_d)-(jfirst-ng_s)+1
      n = imr*j*nl
      allocate(aux(n),stat=ierr)
         if(ierr/=0) call die(myname_,'Alloc(v-aux)',ierr)
      if (inter) then
          call random_number(aux)
      else
          call normalen(n,aux)
      endif
      prog%v(:,jfirst-ng_s:jlast+ng_d,:)=reshape(aux,(/imr,j,nl/))
      deallocate(aux,stat=ierr)
         if(ierr/=0) call die(myname_,'Dealloc(v-aux)',ierr)

      j = jlast-jfirst+1
      n = imr*j*nl
      allocate(aux(n),stat=ierr)
         if(ierr/=0) call die(myname_,'Alloc(delp-aux)',ierr)
      if (inter) then
          call random_number(aux)
      else
          call normalen(n,aux)
      endif
      prog%delp(:,jfirst:jlast,:)=reshape(aux,(/imr,j,nl/))
      deallocate(aux,stat=ierr)
         if(ierr/=0) call die(myname_,'Dealloc(delp-aux)',ierr)

      j = (jlast+ng_d)-(jfirst-ng_d)+1
      n = imr*j*nl
      allocate(aux(n),stat=ierr)
         if(ierr/=0) call die(myname_,'Alloc(pt-aux)',ierr)
      if (inter) then
          call random_number(aux)
      else
          call normalen(n,aux)
      endif
      prog%pt(:,jfirst-ng_d:jlast+ng_d,:)=reshape(aux,(/imr,j,nl/))
      deallocate(aux,stat=ierr)
         if(ierr/=0) call die(myname_,'Dealloc(pt-aux)',ierr)

      n = imr*j*nl*nc
      allocate(aux(n),stat=ierr)
         if(ierr/=0) call die(myname_,'Alloc(q-aux)',ierr)
      if (inter) then
          call random_number(aux)
      else
          call normalen(n,aux)
      endif
      prog%q(:,jfirst-ng_d:jlast+ng_d,:,:)=reshape(aux,(/imr,j,nl,nc/))
      deallocate(aux,stat=ierr)
         if(ierr/=0) call die(myname_,'Dealloc(q-aux)',ierr)

end subroutine prognostics_rand_

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_dup --- duplicate state
!
! !INTERFACE:
subroutine prognostics_dup ( prog_in, prog_out )
      implicit none
! !INPUT PARAMETERS:
      type(dyn_prog), intent(in)  :: prog_in
! !INPUT/OUTPUT PARAMETERS:
      type(dyn_prog), intent(out) :: prog_out
! !DESCRIPTION: duplicate state vector (could be acheived w/ _scal,
!               but this saves the price of multiplications.
!
! !REVISION HISTORY:
!   12May2007 Todling Initial code
!
!EOP
!-----------------------------------------------------------------------

      prog_out%u    = prog_in%u
      prog_out%v    = prog_in%v
      prog_out%pt   = prog_in%pt
      prog_out%delp = prog_in%delp
      prog_out%q    = prog_in%q

end subroutine prognostics_dup

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_axpy --- mimic blas routine to mult vect and add
!
! !INTERFACE:
subroutine prognostics_axpy ( a, x, y ) 
! !USES:
implicit none
! !INPUT PARAMETERS:
real(r8), intent(in) :: a
type(dyn_prog) :: x
! !INPUT/OUTPUT PARAMETERS:
type(dyn_prog) :: y
! !DESCRIPTION: mimics BLAS call to multiply a vector by a constant
!               and add to another vector
!
! !REVISION HISTORY:
!   12May2007 Todling Initial code
!
!EOP
!-----------------------------------------------------------------------

          y%u    =  a * x%u    + y%u
          y%v    =  a * x%v    + y%v
          y%pt   =  a * x%pt   + y%pt
          y%delp =  a * x%delp + y%delp
          y%q    =  a * x%q    + y%q

end subroutine prognostics_axpy 

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prognostics_dopt --- dot product for prog-vector
!
! !INTERFACE:
real(r8) function prognostics_dotp ( x, y )
! !USES:
implicit none
! !INPUT PARAMETERS:
type(dyn_prog) :: x
type(dyn_prog) :: y
! !OUTPUT PARAMETERS:
! !DESCRIPTION: implements a dot product for the prog-vector
!
! !REVISION HISTORY:
!   30May2007 Todling Initial code
!
!EOP
!-----------------------------------------------------------------------
integer, parameter :: ng = 1
real(r8),allocatable :: sum(:)
real(r8) sgp(ng)
integer i,j,k,n,ns
integer ntot,npptot

ns = nc+4
allocate (sum(ns))

ntot= imr*jnp*klast*(nc+4)

sum = 0._r8
sgp = 0._r8
do k = kfirst, klast
   do j = jfirst,jlast
      do i = 1, imr
         sum(1) = sum(1) + y%delp (i,j,k)   * x%delp(i,j,k)
         sum(2) = sum(2) + y%pt   (i,j,k)   * x%pt  (i,j,k)
         sum(3) = sum(3) + y%u    (i,j,k)   * x%u   (i,j,k)
         sum(4) = sum(4) + y%v    (i,j,k)   * x%v   (i,j,k)
         sgp(1) = sgp(1) + 4
      enddo
   enddo
enddo
do n = 1, nc
   do k = kfirst, klast
      do j = jfirst,jlast
         do i = 1, imr
            sum(n+4) = sum(n+4) + y%q(i,j,k,n) * x%q(i,j,k,n)
            sgp(1) = sgp(1) + 1
         enddo
      enddo
   enddo
enddo

  call mp_add1d ( ns, sum )
  call mp_add1d ( ng, sgp )
 
  prognostics_dotp=0._r8
  do i=1,ns
    prognostics_dotp = prognostics_dotp + sum(i)
  enddo
  npptot = nint(sgp(1))
  if(ntot/=npptot)then
     print *, '            total # grid points*(nc+4): ', ntot
     print *, 'counting pp total # grid points*(nc+4): ', sgp(1)
     call die('prognostics_dopt',' inconsistent values ... aborting!')
  endif

  deallocate(sum)

end function prognostics_dotp

end module prognostics
