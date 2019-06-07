#include <misc.h>
#include <preproc.h>

subroutine FireSeason (clm)

!----------------------------------------------------------------------- 
! 
! Purpose: Calculate length of fire season in a year
! 
! Method: Orig. code was called once per day.
!         slevis adapted to call every tstep.
!         Orig. code operated on a grid cell basis.
!         slevis adapted to operate on a patch basis.
! 
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine fire)
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use clmtype
  use clm_varcon, only : tfrz
  use pft_varcon, only : pftpar
  use shr_const_mod, only : SHR_CONST_PI,SHR_CONST_CDAY
  implicit none

! ------------------------ arguments ------------------------------
  type (clm1d), intent(inout) :: clm   !CLM 1-D Module
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  real(r8) :: flam
  real(r8) :: fire_prob
  real(r8) :: pi
! -----------------------------------------------------------------

  pi = SHR_CONST_PI

! Calculate the length of the fire season (in days)

! Calculate today's fire probability, fire_prob
! Assume fire is only possible when temperature is above tfrz
! slevis: *wf is top 0.5 m soil water as a fraction of the whc
!         *divide fire_prob (days) by tsteps/day to get fire_prob (tsteps)
!         *else need daily avg t_ref2m and wf to calc. fire_prob

  flam = pftpar(clm%itypveg,6)

  if (clm%t_ref2m > tfrz .and. clm%litterag > 0.0) then
    fire_prob=EXP((-pi/4.0) * (max(0.0,clm%wf)/flam)**2) * clm%dtime / SHR_CONST_CDAY
  else
    fire_prob=0.0
  endif

  clm%firelength = clm%firelength + fire_prob !reset 1/yr in subroutine lpj

  return
end subroutine FireSeason
