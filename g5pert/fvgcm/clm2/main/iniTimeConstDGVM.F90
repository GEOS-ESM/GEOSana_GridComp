#include <misc.h>
#include <preproc.h>

subroutine iniTimeConstDGVM

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize time invariant clm variables
! 
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use infnan
  use clm_varder
  use clm_varpar, only : nlevsoi
  use clm_varmap, only : begpatch, endpatch
  use clm_varcon, only : spval
  use pft_varcon, only : roota_par, rootb_par, z0mr, displar, &
                         dleaf, rhol, rhos, taul, taus, xl,   &
                         qe25, vcmx25, mp, c3psn
  implicit none

! ------------------------ arguments ---------------------------------
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  integer  :: i,j,k,l,m,ib                      !indices
  integer  :: ivt                               !vegetation type index
  real(r8) :: bd              !bulk density of dry soil material [kg/m^3]
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! Initialize root fraction (computing from surface, d is depth in meter):
! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with beta & d_obs
! given in Zeng et al. (1998).
! --------------------------------------------------------------------

  do k = begpatch, endpatch
     if (.not. clm(k)%lakpoi) then
        ivt = clm(k)%itypveg
        do j = 1, nlevsoi-1
           clm(k)%rootfr(j) = .5*( exp(-roota_par(ivt)*clm(k)%zi(j-1))  &
                                 + exp(-rootb_par(ivt)*clm(k)%zi(j-1))  &
                                 - exp(-roota_par(ivt)*clm(k)%zi(j  ))  &
                                 - exp(-rootb_par(ivt)*clm(k)%zi(j  )) )
        end do
        clm(k)%rootfr(nlevsoi) = .5*( exp(-roota_par(ivt)*clm(k)%zi(nlevsoi-1))  &
                                    + exp(-rootb_par(ivt)*clm(k)%zi(nlevsoi-1)) )
     else
        clm(k)%rootfr(1:nlevsoi) = spval
     end if
  end do

! --------------------------------------------------------------------
! Initialize clm derived type components from pft_varcon to avoid
! indirect addressing and be compatible with offline CLM code
! --------------------------------------------------------------------

  do k = begpatch, endpatch
     ivt = clm(k)%itypveg
     clm(k)%z0mr    = z0mr(ivt)
     clm(k)%displar = displar(ivt)
     clm(k)%dleaf  = dleaf(ivt)
     clm(k)%xl     = xl(ivt)
     do ib = 1,numrad 
        clm(k)%rhol(ib) = rhol(ivt,ib)
        clm(k)%rhos(ib) = rhos(ivt,ib)
        clm(k)%taul(ib) = taul(ivt,ib)
        clm(k)%taus(ib) = taus(ivt,ib)
     end do
     clm(k)%qe25   = qe25(ivt)      ! quantum efficiency at 25c (umol co2 / umol photon)
     clm(k)%vcmx25 = vcmx25(ivt)    ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
     clm(k)%mp     = mp(ivt)        ! slope for conductance-to-photosynthesis relationship
     clm(k)%c3psn  = c3psn(ivt)     ! photosynthetic pathway: 0. = c4, 1. = c3
  end do

  return
end subroutine iniTimeConstDGVM
