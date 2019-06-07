#include <params.h>
      subroutine ccmslp(prsfac, phis  ,ps   ,t1   ,t   ,psl )
C-----------------------------------------------------------------------
C
C Hybrid coord version:  Compute sea level pressure for a latitude line
C
C Algorithm: using the two lowest levels of free atmosphere temperature,
C   define a temperature profile that is linear in the log of hyam +
C   hybm.  Warning: this profile definition may not be reasonable if the
C   hybrid coordinate base pressure is very different than the surface
C   pressure.  Using this profile, extrapolate linearly to the surface
C   and then assume a moist adiabatic lapse rate between surface and sea
C   level to determine a 'sea-level temperature' for evaluating the
C   integral of the hydrostatic equation between surface and sea level.
C
C-------------------------Code History----------------------------------
C
C Original version:  Stolen from the Processor by John Truesdale 
C Standardized:      J. Truesdale, Mar 1994
C                    L. Buja,      Feb 1996
C 
C Modified by SJL Feb 8, 1998
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
      implicit none
#include <pmgrid.h>
#include <comcon.h>
 
C Input
 
      real phis(plond)      ! Surface geopotential height
      real ps(plond)        ! Surface pressure
      real t1(plond)        ! Temperature at plev-1
      real t (plond)        ! Temperature at plev
      real prsfac
C
C Output
C
      real psl(plond)           ! Sea level pressure
C
C---------------------------Local variables-----------------------------
C
      real econst               ! Constant of integration for hydostatic eq.
      real gmoist               ! Moist lapse rate
      real pconst               ! Lapse rate/grav
      integer i                 ! Longitude index
C
C     rair = 287.04, grav = 9.80616
C
      gmoist = 6.5e-3
      econst = gravit / (rair * gmoist)
      pconst = gmoist * rga
 
C Calculate sea level pressure for a latitude
 
      do i=1,plon
        psl(i) = ps(i) * (1. + pconst * phis(i) /
     &          (t(i) + prsfac * (t1(i) - t(i))) ) ** econst
      end do
 
      return
      end
