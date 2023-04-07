module inttgasmod
!$$$ module documentation block
!           .      .    .                                       .
!   module: inttgasmod  module for inttgas and its tangent linear inttgas_tl
!   prgmmr: bweir       org: gmao                date: 2014-04-19
!
! abstract: module for inttgas and its tangent linear inttgas_tl
!
! program history log:
!   2014-04-19  bweir    - initial code, based on carbon monoxide
!
! subroutines included:
!   sub inttgas
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use m_obsNode, only: obsNode
use m_tgasNode, only: tgasNode
use m_tgasNode, only: tgasNode_typecast
use m_tgasNode, only: tgasNode_nextcast
use m_obsdiagNode, only: obsdiagNode_set
use m_obsdiagNode, only: obsdiagNode_get

implicit none

private
public :: inttgas

! Interface for other languages?
interface inttgas
   module procedure inttgas_
end interface

contains

subroutine inttgas_(tgashead, rval, sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    inttgas     apply nonlin qc obs operator for trace gases
!     prgmmr:    bweir       org: gmao                date: 2014-04-19
!
! abstract:      This routine applies the observation operator (forward
!                model) and adjoint of this operator for trace gas
!                observations with the addition of nonlinear qc.
!
! program history log:
!   2014-04-19  weir     - initial code, based on carbon monoxide
!   2014-11-10  weir     - simplified valx calculation
!
!   input argument list:
!     tgashead - trace gas obs type pointer to obs structure
!     sval     - increment in grid space
!
!   output argument list:
!     rval     - results from observation operator 
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use kinds,      only: r_kind, i_kind
  use obsmod,     only: lsaveobsens, l_do_adjoint
  use gridmod,    only: lat2, lon2, nsig
  use jfunc,      only: jiter
  use constants,  only: zero, varlen => max_varname_length
  use gsi_4dvar,  only: ladtest_obs
  use tgasinfo,   only: ntgas, tgnames

  use gsi_bundlemod,     only: gsi_bundle, gsi_bundlegetpointer
  use gsi_chemguess_mod, only: gsi_chemguess_get

  implicit none

! Declare passed variables
  class(obsNode),     pointer, intent(in   ) :: tgashead
  type(gsi_bundle),            intent(in   ) :: sval
  type(gsi_bundle),            intent(inout) :: rval

! Declare local variables
  real(r_kind)    :: val1, valx, w1, w2, w3, w4
  integer(i_kind) :: i, j, k, l, n, ij, ij1, ij2, ij3, ij4, ier, istatus
  integer(i_kind) :: nchanl, navg

  character(len=256) :: obstype

  real(r_kind), allocatable, dimension(:,:)   :: vals, vali, val_ret, tgasak
  real(r_kind), allocatable, dimension(:,:,:) :: stg,    rtg
  real(r_kind), pointer,     dimension(:,:,:) :: stgptr, rtgptr

  type(tgasNode), pointer :: tgasptr

! Return if no trace gas observations
  if (.not. associated(tgashead)) return

! Can't do rank-3 pointer into rank-3, therefore, allocate work space
  allocate(stg(lat2*lon2,nsig,ntgas), rtg(lat2*lon2,nsig,ntgas))

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
  end do

! Loop over trace gas observations
  tgasptr => tgasNode_typecast(tgashead)
  do while (associated(tgasptr))
!    Get level and constituent info
     obstype = tgasptr%obstype
     nchanl  = tgasptr%nchanl
     navg    = tgasptr%navg

!    Set location
     ij1 = tgasptr%ij(1)
     ij2 = tgasptr%ij(2)
     ij3 = tgasptr%ij(3)
     ij4 = tgasptr%ij(4)

     w1 = tgasptr%wij(1)
     w2 = tgasptr%wij(2)
     w3 = tgasptr%wij(3)
     w4 = tgasptr%wij(4)

     allocate(vals(nsig,ntgas), vali(navg,ntgas), val_ret(nchanl,ntgas),       &
              tgasak(navg,ntgas))

!    Accumulate contribution (for all tracers) from model levels
     do n = 1,ntgas
        do l = 1,nsig
           vals(l,n) = w1*stg(ij1,l,n) + w2*stg(ij2,l,n) +                     &
                       w3*stg(ij3,l,n) + w4*stg(ij4,l,n)
        end do

        do j = 1,navg
           vali(j,n) = zero
           do l = 1,nsig
              vali(j,n) = vali(j,n) + tgasptr%avgwgt(j,l)*vals(l,n)
           end do
        end do
     end do

!    Apply the averaging kernel
     do k = 1,nchanl
        n = tgasptr%itgas(k)   ! index of the trace gas observed

        val1 = zero
        do j = 1,navg
           val1 = val1 + tgasptr%avgker(k,j)*vali(j,n)
        end do

        if (lsaveobsens) then
           tgasptr%diags(k)%ptr%obssen(jiter) = val1 * tgasptr%err2(k)         &
                                                     * tgasptr%raterr2(k)
        else
           if (tgasptr%luse) tgasptr%diags(k)%ptr%tldepart(jiter) = val1
        end if

        if (l_do_adjoint) then
           if (lsaveobsens) then
              valx = tgasptr%diags(k)%ptr%obssen(jiter)
           else
              if (ladtest_obs) then
                 valx = val1
              else
                 valx = (val1 - tgasptr%res(k)) * tgasptr%err2(k)              &
                                                * tgasptr%raterr2(k)
              end if
           end if
           val_ret(k,n) = valx
        end if
     end do

!    Spread values to averaging kernel contribution levels
     if (l_do_adjoint) then
        tgasak = zero
        do k = 1,navg
           do j = 1,nchanl
              n = tgasptr%itgas(j)   ! index of the trace gas observed
!             Contribution to kth obs level from jth retrieval channel
              tgasak(k,n) = tgasak(k,n) + tgasptr%avgker(j,k)*val_ret(j,n)
           end do
        end do

!       Adjoint of interpolation
!       Spread each obs level to interpolant gridpoints
        do n = 1,ntgas
           do l = nsig,1,-1
              do j = 1,navg
                 rtg(ij1,l,n) = rtg(ij1,l,n) + w1*tgasptr%avgwgt(j,l)*tgasak(j,n)
                 rtg(ij2,l,n) = rtg(ij2,l,n) + w2*tgasptr%avgwgt(j,l)*tgasak(j,n)
                 rtg(ij3,l,n) = rtg(ij3,l,n) + w3*tgasptr%avgwgt(j,l)*tgasak(j,n)
                 rtg(ij4,l,n) = rtg(ij4,l,n) + w4*tgasptr%avgwgt(j,l)*tgasak(j,n)
              end do
           end do
        end do
     end if ! l_do_adjoint

     deallocate(tgasak, vals, vali, val_ret)

     tgasptr => tgasNode_nextcast(tgasptr)
  end do ! loop over observations

! Copy output and clean up 
  do n = 1,ntgas
!    Retrieve pointers and return if any pointer not found
     call gsi_bundlegetpointer(rval, trim(tgnames(n)), rtgptr, ier)

     if (ier /= 0) return

     do k = 1,nsig
        ij = 0
        do j = 1,lon2
           do i = 1,lat2
              ij = ij+1
              rtgptr(i,j,k) = rtg(ij,k,n)
           end do
        end do
     end do
  end do

  deallocate(stg, rtg)

  return
end subroutine inttgas_

end module inttgasmod
