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

!                This routine splits up the co and co2 arrays to allow for
!                assimilating either or both. The solution is a little hacky,
!                but it works for now. It is important to leave open the
!                possibility of observation operators that depend upon multiple
!                gases or cross-correlations between gases.
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

  use kinds,         only: r_kind, i_kind
  use obsmod,        only: lsaveobsens, l_do_adjoint, luse_obsdiag
  use gridmod,       only: lat2, lon2, nsig
  use jfunc,         only: jiter
  use constants,     only: zero
  use gsi_bundlemod, only: gsi_bundle, gsi_bundlegetpointer
  use gsi_4dvar,     only: ladtest_obs
  use tgasinfo,      only: ihave_co, ihave_co2

  implicit none

! Declare passed variables
  class(obsNode)  , pointer, intent(in   ) :: tgashead
  type(gsi_bundle),            intent(in   ) :: sval
  type(gsi_bundle),            intent(inout) :: rval

! Declare local variables
  real(r_kind)    :: val1, valx, w1, w2, w3, w4
  integer(i_kind) :: i, j, k, l, ij, ij1, ij2, ij3, ij4, ier, istatus
  integer(i_kind) :: nchanl, npro

  character(len=256) :: obstype

  real(r_kind), allocatable, dimension(:)     :: tgasak, vals, vali, val_ret
  real(r_kind), allocatable, dimension(:,:)   :: sco, rco, sco2, rco2
  real(r_kind), pointer,     dimension(:,:,:) :: scoptr,  rcoptr
  real(r_kind), pointer,     dimension(:,:,:) :: sco2ptr, rco2ptr

  type(tgasNode), pointer :: tgasptr

! Return if no trace gas observations
  if (.not. associated(tgashead)) return

! Initialize carbon monoxide arrays
  if (ihave_co) then
!    Retrieve pointers and return if any pointer not found
     ier = 0

     call gsi_bundlegetpointer(sval, 'co', scoptr, istatus)
     ier = istatus + ier
     call gsi_bundlegetpointer(rval, 'co', rcoptr, istatus)
     ier = istatus + ier

     if (ier /= 0) return

!    Can't do rank-2 pointer into rank-2, therefore, allocate work space
     allocate(sco(lat2*lon2,nsig), rco(lat2*lon2,nsig))

     do k = 1,nsig
        ij = 0
        do j = 1,lon2
           do i = 1,lat2
              ij = ij+1
              sco(ij,k) = scoptr(i,j,k)
              rco(ij,k) = rcoptr(i,j,k)
           end do
        end do
     end do
  end if

! Initialize carbon dioxide arrays
  if (ihave_co2) then
!    Retrieve pointers and return if any pointer not found
     ier = 0

     call gsi_bundlegetpointer(sval, 'co2', sco2ptr, istatus)
     ier = istatus + ier
     call gsi_bundlegetpointer(rval, 'co2', rco2ptr, istatus)
     ier = istatus + ier

     if (ier /= 0) return

!    Can't do rank-2 pointer into rank-2, therefore, allocate work space
     allocate(sco2(lat2*lon2,nsig), rco2(lat2*lon2,nsig))

     do k = 1,nsig
        ij = 0
        do j = 1,lon2
           do i = 1,lat2
              ij = ij+1
              sco2(ij,k) = sco2ptr(i,j,k)
              rco2(ij,k) = rco2ptr(i,j,k)
           end do
        end do
     end do
  end if

! Loop over trace gas observations
  !tgasptr => tgashead
  tgasptr => tgasNode_typecast(tgashead)
  do while (associated(tgasptr))
!    Get level and constituent info
     obstype = tgasptr%obstype
     nchanl  = tgasptr%nchanl
     npro    = tgasptr%npro

!    Set location
     ij1 = tgasptr%ij(1)
     ij2 = tgasptr%ij(2)
     ij3 = tgasptr%ij(3)
     ij4 = tgasptr%ij(4)

     w1 = tgasptr%wij(1)
     w2 = tgasptr%wij(2)
     w3 = tgasptr%wij(3)
     w4 = tgasptr%wij(4)

     allocate(tgasak(npro), vals(nsig), vali(npro), val_ret(nchanl))

!    Accumulate contribution from model levels
     if (trim(obstype) == 'mopitt' .and. ihave_co) then
        do l = 1,nsig
           vals(l) = w1*sco(ij1,l) + w2*sco(ij2,l) +   &
                     w3*sco(ij3,l) + w4*sco(ij4,l)
        end do

        do j = 1,npro
           vali(j) = zero
           do l = 1,nsig
              vali(j) = vali(j) + tgasptr%avgwgt(j,l)*vals(l)
           end do
        end do

     else if (trim(obstype) == 'acos' .and. ihave_co2) then
        do l = 1,nsig
           vals(l) = w1*sco2(ij1,l) + w2*sco2(ij2,l) + &
                     w3*sco2(ij3,l) + w4*sco2(ij4,l)
        end do

        do j = 1,npro
           vali(j) = zero
           do l = 1,nsig
              vali(j) = vali(j) + tgasptr%avgwgt(j,l)*vals(l)
           end do
        end do

     else if (trim(obstype) == 'flask') then
         vali = zero
         if (ihave_co) then
            vali(1) = w1*sco(ij1,1)  + w2*sco(ij2,1) +  &
                      w3*sco(ij3,1)  + w4*sco(ij4,1)
         end if
         if (ihave_co2) then
            vali(2) = w1*sco2(ij1,1) + w2*sco2(ij2,1) + &
                      w3*sco2(ij3,1) + w4*sco2(ij4,1)
         end if
     end if

!    Apply the averaging kernel
     do k = 1,nchanl
        val1 = zero
        do j = 1,npro
           val1 = val1 + tgasptr%avgker(k,j)*vali(j)
        end do

       if(luse_obsdiag) then
         if (lsaveobsens) then
           valx = val1 * tgasptr%err2(k) * tgasptr%raterr2(k)
           call obsdiagNode_set(tgasptr%diags(k)%ptr,jiter=jiter, &
             obssen = valx )
         else
           if (tgasptr%luse) call obsdiagNode_set(tgasptr%diags(k)%ptr,jiter=jiter, &
             tldepart = val1 )
         end if
       endif

        if (l_do_adjoint) then
           if (lsaveobsens) then
              ! valx = tgasptr%diags(k)%ptr%obssen(jiter)
              call obsdiagNode_get(tgasptr%diags(k)%ptr,jiter=jiter,obssen=valx)
           else
              if (ladtest_obs) then
                 valx = val1
              else
                 valx = (val1 - tgasptr%res(k)) * tgasptr%err2(k)              &
                                                * tgasptr%raterr2(k)
              end if
           end if
           val_ret(k) = valx
        end if
     end do

!    Spread values to averaging kernel contribution levels
     if (l_do_adjoint) then
        do k = 1,npro
           tgasak(k) = zero
           do j = 1,nchanl
!             Contribution to kth obs level from jth retrieval channel
              tgasak(k) = tgasak(k) + tgasptr%avgker(j,k)*val_ret(j)
           end do
        end do

!       Adjoint of interpolation
!       Spread each obs level to interpolant gridpoints
!       Maybe could set this up so it just checks ihave_co and ihave_co2?
        if (trim(obstype) == 'mopitt' .and. ihave_co) then
           do l = nsig,1,-1
              do j = 1,npro
                 rco(ij1,l) = rco(ij1,l) + w1*tgasptr%avgwgt(j,l)*tgasak(j)
                 rco(ij2,l) = rco(ij2,l) + w2*tgasptr%avgwgt(j,l)*tgasak(j)
                 rco(ij3,l) = rco(ij3,l) + w3*tgasptr%avgwgt(j,l)*tgasak(j)
                 rco(ij4,l) = rco(ij4,l) + w4*tgasptr%avgwgt(j,l)*tgasak(j)
              end do
           end do

        else if (trim(obstype) == 'acos' .and. ihave_co2) then
           do l = nsig,1,-1
              do j = 1,npro
                 rco2(ij1,l) = rco2(ij1,l) + w1*tgasptr%avgwgt(j,l)*tgasak(j)
                 rco2(ij2,l) = rco2(ij2,l) + w2*tgasptr%avgwgt(j,l)*tgasak(j)
                 rco2(ij3,l) = rco2(ij3,l) + w3*tgasptr%avgwgt(j,l)*tgasak(j)
                 rco2(ij4,l) = rco2(ij4,l) + w4*tgasptr%avgwgt(j,l)*tgasak(j)
              end do
           end do

        else if (trim(obstype) == 'flask') then
           if (ihave_co) then
              rco(ij1,1) = rco(ij1,1) + w1*tgasak(1)
              rco(ij2,1) = rco(ij2,1) + w2*tgasak(1)
              rco(ij3,1) = rco(ij3,1) + w3*tgasak(1)
              rco(ij4,1) = rco(ij4,1) + w4*tgasak(1)
           end if

           if (ihave_co2) then
              rco2(ij1,1) = rco2(ij1,1) + w1*tgasak(2)
              rco2(ij2,1) = rco2(ij2,1) + w2*tgasak(2)
              rco2(ij3,1) = rco2(ij3,1) + w3*tgasak(2)
              rco2(ij4,1) = rco2(ij4,1) + w4*tgasak(2)
           end if
        end if
     end if ! l_do_adjoint

     deallocate(tgasak, vals, vali, val_ret)

     !tgasptr => tgasptr%llpoint
     tgasptr => tgasNode_nextcast(tgasptr)

  end do ! loop over observations

! Copy output and clean up 
  do k = 1,nsig
     ij = 0
     do j = 1,lon2
        do i = 1,lat2
           ij = ij+1
           if (ihave_co)  rcoptr(i,j,k)  = rco(ij,k)
           if (ihave_co2) rco2ptr(i,j,k) = rco2(ij,k)
        end do
     end do
  end do

  if (ihave_co)  deallocate(sco,  rco)
  if (ihave_co2) deallocate(sco2, rco2)

  return
end subroutine inttgas_

end module inttgasmod
