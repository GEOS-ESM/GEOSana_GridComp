module stptgasmod

!$$$ module documentation block
!           .        .    .                                       .
! module:   stptgasmod    module for stptgas and its tangent linear stptgas_tl
!  prgmmr:
!
! abstract: module for stptgas and its tangent linear stptgas_tl
!
! program history log:
!   2014-08-22  weir    - initial code, based on stpozmod
!   2021-05-10  guo     - replaced ob_type by polymorphic obsNode with type casting
!
! subroutines included:
!   sub stptgas
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

private
public :: stptgas

contains

subroutine stptgas(tgashead, rval, sval, out, sges, nstep)
!$$$  subprogram documentation block
!                .     .    .                                       .
! subprogram:    stptgas    compute contribution to penalty and stepsize for
!                           trace gases using nonlinear qc
!   prgmmr: weir            org: gmao                 date: 2014-08-22
!
! abstract: The routine computes the contribution to the penalty from trace
!           gas observations.  The routine also computes the contribution of
!           trace gas observations to the step size.  This version includes
!           nonlinear qc.
!
!  warning: I have no idea what this subroutine is supposed to compute, it
!           was merely translated from other stp routines and tested a
!           couple of times. It could contain serious errors (weir).
!
!
! program history log:
!   2014-08-22  weir    - initial code, based on stpozlay_
!
!   input argument list:
!     tgashead - trace gas obs type pointer to obs structure
!     rval     - search direction for trace gas
!     sval     - input trace gas correction field
!     sges     - step size estimates (nstep)
!     nstep    - number of stepsize estimates (== 0 means use outer iteration
!                value)
!
!   output argument list:
!     out(1:nstep) - contribution of trace gas data to penalty sges(1:nstep)
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

  use m_obsNode,  only: obsNode
  use m_tgasNode, only: tgasNode
  use m_tgasNode, only: tgasNode_typecast
  use m_tgasNode, only: tgasNode_nextcast

  use kinds,      only: r_kind, i_kind, r_quad
  use constants,  only: zero, zero_quad, one, r3600,                           &
                        varlen => max_varname_length
  use gridmod,    only: lat2, lon2, nsig
! use jfunc,      only: l_foto, xhat_dt, dhat_dt
  use tgasinfo,   only: ntgas, tgnames

  use gsi_bundlemod,     only: gsi_bundle, gsi_bundlegetpointer
  use gsi_chemguess_mod, only: gsi_chemguess_get

  implicit none

! Declare passed variables
  class(obsNode),     pointer, intent(in   ) :: tgashead
  integer(i_kind),             intent(in   ) :: nstep
  real(r_quad),                intent(inout) :: out(max(1,nstep))
  type(gsi_bundle),            intent(in   ) :: rval, sval
  real(r_kind),                intent(in   ) :: sges(max(1,nstep))

! Declare local variables
  integer(i_kind) :: i, j, k, l, n, ij, is, ier, istatus
  integer(i_kind) :: ij1, ij2, ij3, ij4
! integer(i_kind) :: ijx, ij1x, ij2x, ij3x, ij4x
  integer(i_kind) :: nchanl, navg

  real(r_kind) :: w1, w2, w3, w4, time, tgas
  real(r_quad) :: val, val1, val_lay, val1_lay, valr, vals

  character(len=256) :: obstype

  real(r_kind), dimension(max(1,nstep)) :: pen

  real(r_kind), allocatable, dimension(:,:,:) :: stg,    rtg
  real(r_kind), pointer,     dimension(:,:,:) :: stgptr, rtgptr
! real(r_kind), allocatable, dimension(:,:)   :: xtg,    dtg
! real(r_kind), pointer,     dimension(:)     :: xtgptr, dtgptr

  type(tgasNode), pointer :: tgasptr

! Return if no trace gas observations
  if (.not. associated(tgashead)) return

! Can't do rank-3 pointer into rank-3, therefore, allocate work space
  allocate(stg(lat2*lon2,nsig,ntgas))
  allocate(rtg(lat2*lon2,nsig,ntgas))
! allocate(xtg(lat2*lon2*nsig,ntgas))
! allocate(dtg(lat2*lon2*nsig,ntgas))

! Initialize trace gas arrays
  do n = 1,ntgas
!    Retrieve pointers and return if any pointer not found
     istatus = 0
     call gsi_bundlegetpointer(sval, trim(tgnames(n)), stgptr, ier)
     istatus = istatus + ier
     call gsi_bundlegetpointer(rval, trim(tgnames(n)), rtgptr, ier)
     istatus = istatus + ier

!    Probably could be a little more informative (fixme)
     if (istatus /= 0) return

     do k = 1,nsig
        ij = 0
        do j = 1,lon2
           do i = 1,lat2
              ij = ij+1
              stg(ij,k,n) = stgptr(i,j,k)
              rtg(ij,k,n) = rtgptr(i,j,k)
           end do
        end do
     end do

!     if (l_foto) then
!        call gsi_bundlegetpointer(xhat_dt, trim(tgnames(n)), xtgptr, ier)
!        istatus = istatus + ier
!        call gsi_bundlegetpointer(dhat_dt, trim(tgnames(n)), dtgptr, ier)
!        istatus = istatus + ier
!
!!       Probably could be a little more informative (fixme)
!        if (istatus /= 0) return
!
!        ijx = 0
!        do k = 1,nsig
!           do j = 1,lon2
!              do i = 1,lat2
!                 ijx = ijx+1
!                 xtg(ijx,n) = xtgptr(ijx)
!                 dtg(ijx,n) = dtgptr(ijx)
!              end do
!           end do
!        end do
!     end if
  end do

! Loop over trace gas observations
  tgasptr => tgasNode_typecast(tgashead)
  do while (associated(tgasptr))
     if (tgasptr%luse) then
!       Get level and constituent info
        obstype = tgasptr%obstype
        nchanl  = tgasptr%nchanl
        navg    = tgasptr%navg

        if (0 < nstep) then
!          Get location
           ij1 = tgasptr%ij(1)
           ij2 = tgasptr%ij(2)
           ij3 = tgasptr%ij(3)
           ij4 = tgasptr%ij(4)

           w1 = tgasptr%wij(1)
           w2 = tgasptr%wij(2)
           w3 = tgasptr%wij(3)
           w4 = tgasptr%wij(4)

!          if (l_foto) time = tgasptr%time*r3600
        end if

        do k = 1,nchanl
           n = tgasptr%itgas(k)   ! index of the trace gas observed

           if (0 < nstep) then
              val  =  zero_quad
              val1 = -tgasptr%res(k)

              do j = 1,navg
                 val_lay  = zero_quad
                 val1_lay = zero_quad

                 do l = 1,nsig
                    valr = w1*rtg(ij1,l,n) + w2*rtg(ij2,l,n) + &
                           w3*rtg(ij3,l,n) + w4*rtg(ij4,l,n)
                    vals = w1*stg(ij1,l,n) + w2*stg(ij2,l,n) + &
                           w3*stg(ij3,l,n) + w4*stg(ij4,l,n)

                    val_lay  = val_lay  + tgasptr%avgwgt(j,l)*valr
                    val1_lay = val1_lay + tgasptr%avgwgt(j,l)*vals

!                    if (l_foto) then
!                       ij1x = ij1 + (l-1)*lat2*lon2
!                       ij2x = ij2 + (l-1)*lat2*lon2
!                       ij3x = ij3 + (l-1)*lat2*lon2
!                       ij4x = ij4 + (l-1)*lat2*lon2
!
!                       valr = w1*dtg(ij1x,n) + w2*dtg(ij2x,n) + &
!                              w3*dtg(ij3x,n) + w4*dtg(ij4x,n)
!                       vals = w1*xtg(ij1x,n) + w2*xtg(ij2x,n) + &
!                              w3*xtg(ij3x,n) + w4*xtg(ij4x,n)
!
!                       val_lay  = val_lay  + time*tgasptr%avgwgt(j,l)*valr
!                       val1_lay = val1_lay + time*tgasptr%avgwgt(j,l)*vals
!                    end if
                 end do

                 val  = val  + tgasptr%avgker(k,j)*val_lay
                 val1 = val1 + tgasptr%avgker(k,j)*val1_lay
              end do

              do is = 1,nstep
                 tgas = val1 + sges(is)*val
                 pen(is) = tgasptr%err2(k)*tgas*tgas
              end do
           else
              pen(1) = tgasptr%res(k)*tgasptr%res(k)*tgasptr%err2(k)
           end if

           out(1) = out(1) + pen(1)*tgasptr%raterr2(k)
           do is = 2,nstep
              out(is) = out(is) + (pen(is) - pen(1))*tgasptr%raterr2(k)
           end do
        end do
     end if ! luse

     tgasptr => tgasNode_nextcast(tgasptr)
  end do ! loop over observations

! Clean 
  deallocate(stg, rtg)
! deallocate(xtg, dtg)

  return
end subroutine stptgas

end module stptgasmod
