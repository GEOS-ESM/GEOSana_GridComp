!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  geos_gsi_util --- utilities for GEOS-GSI coupler
!
! !INTERFACE:
!
module geos_gsi_util
use mpeu_util, only: perr,die

character(len=*), parameter :: myname='geos_cplr_util'
public :: geos_gsi_divvor
public :: geos_gsi_getqs
  interface geos_gsi_divvor
    module procedure divvor8_
    module procedure divvor4_
  end interface
  interface geos_gsi_getqs
    module procedure getqs_
  end interface

! !DESCRIPTION: general utilites needed by GEOS-GSI coupler routines.
!
! !REVISION HISTORY:
!
!  13Nov2011  Todling   Initial code.
!
!EOP
!-----------------------------------------------------------------------
contains
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: divvor8_:  Calculate voricity/divergence using TRANSF
!
! !INTERFACE:

   subroutine divvor8_ (grid,rlats,u,v,div,vor)

! !USES:
   use kinds, only: i_kind,r_double
   use mpimod, only: mype,mpi_comm_world
   use constants, only: zero

   use general_sub2grid_mod, only: sub2grid_info

   use m_ggGradient,only : ggGradient
   use m_ggGradient,only : ggGradient_init,clean
   use m_ggGradient,only : ggDivo

   implicit none

! !INPUT PARAMETERS:

   type(sub2grid_info),intent(in) :: grid
   real(r_double),intent(in)   ,dimension(:)    ::rlats
   real(r_double),intent(in)   ,dimension(:,:,:)::u,v

! !OUTPUT PARAMETERS:

   real(r_double),intent(inout),dimension(:,:,:)::div,vor

! !DESCRIPTION: double precision interface to calc of vor/div using
!               Jing's TRANSF routines.
!
! !REVISION HISTORY:
!
!  13Nov2011  Todling   Initial code.
!
!EOP
!-----------------------------------------------------------------------
   character(len=*), parameter :: myname_=myname//'*divvor8_'
   type(ggGradient) :: gr
   integer :: iGdim,iGloc,iGlen
   integer :: jGdim,jGloc,jGlen
   real(r_double),dimension(:,:,:), allocatable :: up8
   real(r_double),dimension(:,:,:), allocatable :: vp8

   vor=zero
   div=zero
   iGloc=grid%jstart(mype+1)
   iGlen=grid%jlon1 (mype+1)
   jGloc=grid%istart(mype+1)
   jGlen=grid%ilat1 (mype+1)

   call ggGradient_init(gr,grid%nlon,rlats,   &
        iGloc,iGlen,jGloc,jGlen,grid%nsig, mpi_comm_world)

   call ggDivo(gr,  &
        u,          &       ! in: u
        v,          &       ! in: v
        div,        &       ! div(u,v)
        vor,        &       ! vor(u,v)
        mpi_comm_world)     ! communicator

   call clean(gr)
   end subroutine divvor8_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: divvor4_:  Calculate voricity/divergence using TRANSF
!
! !INTERFACE:

   subroutine divvor4_ (grid,rlats,u,v,div,vor)

! !USES:
   use kinds, only: i_kind,r_single,r_double
   use mpimod, only: mype,mpi_comm_world
   use constants, only: zero

   use general_sub2grid_mod, only: sub2grid_info

   use m_ggGradient,only : ggGradient
   use m_ggGradient,only : ggGradient_init,clean
   use m_ggGradient,only : ggDivo

   implicit none

! !INPUT PARAMETERS:

   type(sub2grid_info),intent(in) :: grid
   real(r_double),intent(in)   ,dimension(:)    ::rlats
   real(r_single),intent(in),   dimension(:,:,:)::u,v

! !OUTPUT PARAMETERS:

   real(r_single),intent(inout),dimension(:,:,:)::div,vor

! !DESCRIPTION: double precision interface to calc of vor/div using
!               Jing's TRANSF routines.
!
! !REVISION HISTORY:
!
!  13Nov2011  Todling   Initial code.
!
!EOP
!-----------------------------------------------------------------------
   character(len=*), parameter :: myname_=myname//'*divvor4_'
   type(ggGradient) :: gr
   integer    ier
   integer :: iGdim,iGloc,iGlen
   integer :: jGdim,jGloc,jGlen
   real(r_double),dimension(:,:,:), allocatable :: up8
   real(r_double),dimension(:,:,:), allocatable :: vp8
   real(r_double),dimension(:,:,:), allocatable :: dv8
   real(r_double),dimension(:,:,:), allocatable :: vr8

   iGloc=grid%jstart(mype+1)
   iGlen=grid%jlon1 (mype+1)
   jGloc=grid%istart(mype+1)
   jGlen=grid%ilat1 (mype+1)

   allocate(up8 (grid%lat2,grid%lon2,grid%nsig), &
            vp8 (grid%lat2,grid%lon2,grid%nsig), &
            dv8 (grid%lat2,grid%lon2,grid%nsig), &
            vr8 (grid%lat2,grid%lon2,grid%nsig), &
            stat=ier)
        if(ier/=0) then
          call perr(myname_,'alloc()')
          call die(myname_)
        endif

   up8=u
   vp8=v
   dv8=zero
   vr8=zero

   call ggGradient_init(gr,grid%nlon,rlats,   &
        iGloc,iGlen,jGloc,jGlen,grid%nsig, mpi_comm_world)

   call ggDivo(gr,  &
        up8,        &       ! in: u
        vp8,        &       ! in: v
        dv8,        &       ! div(u,v)
        vr8,        &       ! vor(u,v)
        mpi_comm_world)     ! communicator

   div=dv8
   vor=vr8

   call clean(gr)
   deallocate(up8, vp8, dv8, vr8, stat=ier)
        if(ier/=0) then
          call perr(myname_,'dealloc()')
          call die(myname_)
        endif

   end subroutine divvor4_

! the following is sort of redundant w/ get_gefs_ensperts_dualres.f90
! TBD: need to reconcile
   subroutine getqs_ (ps,q,tv,im,jm,km,qs)
 
   use kinds, only: i_kind,r_kind
   use gridmod, only: idsl5
   use constants, only: zero,one,half,fv,rd_over_cp
   implicit none

   integer(i_kind),intent(in)  :: im,jm,km
   real(r_kind),   intent(in)  :: ps(im,jm)
   real(r_kind),   intent(in)  :: tv(im,jm,km)
   real(r_kind),   intent(in)  :: q(im,jm,km)
   real(r_kind),   intent(out) :: qs(im,jm,km)

! locals
   integer(i_kind) i,j,k,iderivative
   logical ice
   real(r_kind) bar_norm,sig_norm,kapr,kap1,rh
   real(r_kind),allocatable,dimension(:,:,:) :: tsen,prsl,pri

  kap1=rd_over_cp+one
  kapr=one/rd_over_cp

! Compute RH
! Get 3d pressure field now on interfaces
    allocate(pri(im,jm,km+1))
    call general_getprs_glb(ps,tv,pri)
    allocate(prsl(im,jm,km),tsen(im,jm,km))
! Get sensible temperature and 3d layer pressure
    if (idsl5 /= 2) then
!$omp parallel do schedule(dynamic,1) private(k,j,i)
      do k=1,km
        do j=1,jm
          do i=1,im
            prsl(i,j,k)=((pri(i,j,k)**kap1-pri(i,j,k+1)**kap1)/&
                           (kap1*(pri(i,j,k)-pri(i,j,k+1))))**kapr
            tsen(i,j,k)= tv(i,j,k)/(one+fv*max(zero,q(i,j,k)))
          end do
        end do
      end do
    else
!$omp parallel do schedule(dynamic,1) private(k,j,i)
      do k=1,km
        do j=1,jm
          do i=1,im
            prsl(i,j,k)=(pri(i,j,k)+pri(i,j,k+1))*half
            tsen(i,j,k)= tv(i,j,k)/(one+fv*max(zero,q(i,j,k)))
          end do
        end do
      end do
    end if
    deallocate(pri)

    ice=.true.
    iderivative=0
    call genqsat(qs,tsen,prsl,im,jm,km,ice,iderivative)
    deallocate(tsen,prsl)

   end subroutine getqs_

end module geos_gsi_util
