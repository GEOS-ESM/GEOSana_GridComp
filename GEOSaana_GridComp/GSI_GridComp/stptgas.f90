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
  use kinds,         only: r_kind, i_kind, r_quad
  use obsmod,        only: tgas_ob_type
  use constants,     only: zero, zero_quad, one, r3600
  use gridmod,       only: lat2, lon2, nsig
  use jfunc,         only: l_foto, xhat_dt, dhat_dt
  use gsi_bundlemod, only: gsi_bundle, gsi_bundlegetpointer
  use tgasinfo,      only: ihave_co, ihave_co2

  implicit none

! Declare passed variables
  type(tgas_ob_type), pointer, intent(in   ) :: tgashead
  integer(i_kind),             intent(in   ) :: nstep
  real(r_quad),                intent(inout) :: out(max(1,nstep))
  type(gsi_bundle),            intent(in   ) :: rval, sval
  real(r_kind),                intent(in   ) :: sges(max(1,nstep))

! Declare local variables
  integer(i_kind) :: i, j, k, l, ij, is, ier, istatus
  integer(i_kind) :: ij1, ij2, ij3, ij4, ij1x, ij2x, ij3x, ij4x
  integer(i_kind) :: nchanl, npro

  real(r_kind) :: w1, w2, w3, w4, time, tgas
  real(r_quad) :: val, val1, val_lay, val1_lay, valr, vals

  character(len=256) :: obstype

  real(r_kind), dimension(max(1,nstep)) :: pen

  real(r_kind), pointer,     dimension(:)     :: xhat_dt_co,  dhat_dt_co
  real(r_kind), pointer,     dimension(:)     :: xhat_dt_co2, dhat_dt_co2
  real(r_kind), allocatable, dimension(:,:)   :: sco, rco, sco2, rco2
  real(r_kind), pointer,     dimension(:,:,:) :: scoptr,  rcoptr
  real(r_kind), pointer,     dimension(:,:,:) :: sco2ptr, rco2ptr

  type(tgas_ob_type), pointer :: tgasptr

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

     if (l_foto) then
        call gsi_bundlegetpointer(xhat_dt, 'co', xhat_dt_co, istatus)
        ier = istatus + ier
        call gsi_bundlegetpointer(dhat_dt, 'co', dhat_dt_co, istatus)
        ier = istatus + ier
     end if

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

     if (l_foto) then
        call gsi_bundlegetpointer(xhat_dt, 'co2', xhat_dt_co2, istatus)
        ier = istatus + ier
        call gsi_bundlegetpointer(dhat_dt, 'co2', dhat_dt_co2, istatus)
        ier = istatus + ier
     end if

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
  tgasptr => tgashead
  do while (associated(tgasptr))
     if (tgasptr%luse) then
!       Get level and constituent info
        obstype = tgasptr%obstype
        nchanl  = tgasptr%nchanl
        npro    = tgasptr%npro

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

           if (l_foto) time = tgasptr%time*r3600
        end if

        do k = 1,nchanl
           if (0 < nstep) then
              val  =  zero_quad
              val1 = -tgasptr%res(k)

              if (trim(obstype) == 'mopitt' .and. ihave_co) then
                 do j = 1,npro
                    val_lay  = zero_quad
                    val1_lay = zero_quad

                    do l = 1,nsig
                       valr = w1*rco(ij1,l) + w2*rco(ij2,l) + &
                              w3*rco(ij3,l) + w4*rco(ij4,l)
                       vals = w1*sco(ij1,l) + w2*sco(ij2,l) + &
                              w3*sco(ij3,l) + w4*sco(ij4,l)

                       val_lay  = val_lay  + tgasptr%avgwgt(j,l)*valr
                       val1_lay = val1_lay + tgasptr%avgwgt(j,l)*vals

                       if (l_foto) then
                          ij1x = ij1 + (l-1)*lat2*lon2
                          ij2x = ij2 + (l-1)*lat2*lon2
                          ij3x = ij3 + (l-1)*lat2*lon2
                          ij4x = ij4 + (l-1)*lat2*lon2
  
                          valr = w1*dhat_dt_co(ij1x) + w2*dhat_dt_co(ij2x) + &
                                 w3*dhat_dt_co(ij3x) + w4*dhat_dt_co(ij4x)
                          vals = w1*xhat_dt_co(ij1x) + w2*xhat_dt_co(ij2x) + &
                                 w3*xhat_dt_co(ij3x) + w4*xhat_dt_co(ij4x)
 
                          val_lay  = val_lay  + time*tgasptr%avgwgt(j,l)*valr
                          val1_lay = val1_lay + time*tgasptr%avgwgt(j,l)*vals
                       end if
                    end do

                    val  = val  + tgasptr%avgker(k,j)*val_lay
                    val1 = val1 + tgasptr%avgker(k,j)*val1_lay
                 end do

              else if (trim(obstype) == 'acos' .and. ihave_co2) then
                 do j = 1,npro
                    val_lay  = zero_quad
                    val1_lay = zero_quad

                    do l = 1,nsig
                       valr = w1*rco2(ij1,l) + w2*rco2(ij2,l) + &
                              w3*rco2(ij3,l) + w4*rco2(ij4,l)
                       vals = w1*sco2(ij1,l) + w2*sco2(ij2,l) + &
                              w3*sco2(ij3,l) + w4*sco2(ij4,l)

                       val_lay  = val_lay  + tgasptr%avgwgt(j,l)*valr
                       val1_lay = val1_lay + tgasptr%avgwgt(j,l)*vals

                       if (l_foto) then
                          ij1x = ij1 + (l-1)*lat2*lon2
                          ij2x = ij2 + (l-1)*lat2*lon2
                          ij3x = ij3 + (l-1)*lat2*lon2
                          ij4x = ij4 + (l-1)*lat2*lon2
  
                          valr = w1*dhat_dt_co2(ij1x) + w2*dhat_dt_co2(ij2x) + &
                                 w3*dhat_dt_co2(ij3x) + w4*dhat_dt_co2(ij4x)
                          vals = w1*xhat_dt_co2(ij1x) + w2*xhat_dt_co2(ij2x) + &
                                 w3*xhat_dt_co2(ij3x) + w4*xhat_dt_co2(ij4x)
 
                          val_lay  = val_lay  + time*tgasptr%avgwgt(j,l)*valr
                          val1_lay = val1_lay + time*tgasptr%avgwgt(j,l)*vals
                       end if
                    end do

                    val  = val  + tgasptr%avgker(k,j)*val_lay
                    val1 = val1 + tgasptr%avgker(k,j)*val1_lay
                 end do

              else if (trim(obstype) == 'flask') then
                 if (ihave_co) then
                    val_lay  = w1*rco(ij1,1) + w2*rco(ij2,1) &
                             + w3*rco(ij3,1) + w4*rco(ij4,1)
                    val1_lay = w1*sco(ij1,1) + w2*sco(ij2,1) &
                             + w3*sco(ij3,1) + w4*sco(ij4,1)
   
                    if (l_foto) then
                       val_lay  = val_lay  + time*(  w1*dhat_dt_co(ij1) &
                                                   + w2*dhat_dt_co(ij2) &
                                                   + w3*dhat_dt_co(ij3) &
                                                   + w4*dhat_dt_co(ij4) )
                       val1_lay = val1_lay + time*(  w1*xhat_dt_co(ij1) &
                                                   + w2*xhat_dt_co(ij2) &
                                                   + w3*xhat_dt_co(ij3) &
                                                   + w4*xhat_dt_co(ij4) )
                    end if

                    val  = val  + tgasptr%avgker(k,1)*val_lay
                    val1 = val1 + tgasptr%avgker(k,1)*val1_lay
                 end if

                 if (ihave_co2) then
                    val_lay  = w1*rco2(ij1,1) + w2*rco2(ij2,1) &
                             + w3*rco2(ij3,1) + w4*rco2(ij4,1)
                    val1_lay = w1*sco2(ij1,1) + w2*sco2(ij2,1) &
                             + w3*sco2(ij3,1) + w4*sco2(ij4,1)
                    if (l_foto) then
                       val_lay  = val_lay  + time*(  w1*dhat_dt_co2(ij1) &
                                                   + w2*dhat_dt_co2(ij2) &
                                                   + w3*dhat_dt_co2(ij3) &
                                                   + w4*dhat_dt_co2(ij4) )
                       val1_lay = val1_lay + time*(  w1*xhat_dt_co2(ij1) &
                                                   + w2*xhat_dt_co2(ij2) &
                                                   + w3*xhat_dt_co2(ij3) &
                                                   + w4*xhat_dt_co2(ij4) )
                    end if

                    val  = val  + tgasptr%avgker(k,2)*val_lay
                    val1 = val1 + tgasptr%avgker(k,2)*val1_lay
                 end if

              end if

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
     end if

     tgasptr => tgasptr%llpoint
  end do ! loop over observations

! Clean 
  if (ihave_co)  deallocate(sco,  rco)
  if (ihave_co2) deallocate(sco2, rco2)

  return

end subroutine stptgas

end module stptgasmod
