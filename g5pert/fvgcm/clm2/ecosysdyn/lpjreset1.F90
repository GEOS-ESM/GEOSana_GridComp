#include <misc.h>
#include <preproc.h>

subroutine lpjreset1 (caldayp1, eccen, obliqr, lambm0, mvelpp)
 
  use precision
  use clmtype
  use clm_varder
  use clm_varmap, only : begpatch,endpatch
  implicit none

! ----------------------------------------------------------------------
! purpose           : to reset variables related to lpj
! date first created: November 2000 - lsm version 2 
! by whom           : Sam Levis
! date last revised : 
! by whom           : 
! ----------------------------------------------------------------------
! $Id$
! ----------------------------------------------------------------------

! -------------------------- arguments ---------------------------------
  real(r8), intent(in) :: caldayp1 !calendar day at Greenwich (1.00, ..., 365.99) for nstep+1
  real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox long. of perihelion + pi (radians)
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
  integer  :: k                    !indices 
! ----------------------------------------------------------------------

! reset a few variables here at the very end of the year

  do k = begpatch,endpatch
     clm(k)%annpsn     = 0.
     clm(k)%annpsnpot  = 0.
     clm(k)%bm_inc     = 0.
     clm(k)%afmicr     = 0.
     clm(k)%firelength = 0.
     clm(k)%agddtw     = 0.
     clm(k)%agdd       = 0.
     clm(k)%t10min     = 1.0e+36
     clm(k)%t_mo_min   = 1.0e+36
  end do

! call EcosystemDyn because need info for first timestep of next year

  do k = begpatch,endpatch

     call EcosystemDyn (clm(k), .false., .true.)

     call SurfaceAlbedo (clm(k), caldayp1, eccen, obliqr, lambm0, mvelpp)

  end do

  call iniTimeConstDGVM ()

  return
end subroutine lpjreset1
