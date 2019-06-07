#include <params.h>
      subroutine zenith(calday  ,dodiavg ,clat    ,coszrs  )
C-----------------------------------------------------------------------
C
C Modified for generalized orbit
C 19 November 1996    Bruce P. Briegleb
C
C Compute cosine of solar zenith angle for albedo and radiation 
C   computations.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski, May 1994
C Standardized:      L. Buja,     Feb 1996
C Reviewed:          J. Kiehl,    Apr 1996    
C 
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <crdcon.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real calday              ! Calendar day, including fraction
      logical dodiavg          ! true => do diurnal averaging
      real clat                ! Current latitude (radians)
C
C Output arguments
C
      real coszrs(plond)       ! Cosine solar zenith angle
C
C---------------------------Local variables-----------------------------
C
      integer i     ! Longitude loop index
      real phi      ! Greenwich calendar day + local time + long offset
      real delta    ! Solar declination angle  in radians
      real eccf     ! Earth orbit eccentricity factor
      real sinc     ! Sine   of latitude
      real cosc     ! Cosine of latitude
      real sind     ! Sine   of declination
      real cosd     ! Cosine of declination
      real frac     ! Daylight fraction
      real arg      ! Cosine of the half day length (P&P eq. 3.21)
      real tsun     ! Temporary term in diurnal averaging
      real coszrsu  ! Uniform cosine zenith solar angle 
C
C-----------------------------------------------------------------------
C
      call orbit(calday,delta,eccf)
C
C Compute local cosine solar zenith angle,
C
      sinc = sin(clat)
      sind = sin(delta)
      cosc = cos(clat)
C SJL
      if(cosc .lt. 1.E-6) cosc = 0.
C SJL

      cosd = cos(delta)
C
C If using diurnal averaging, then compute the average local cosine solar 
C zenith angle using formulas from paltridge and platt 1976  p. 57, p. 62,63.
C
      if (dodiavg) then
        arg = -(sinc/cosc)*(sind/cosd)
        if (arg .lt. -1.) then
          frac = 1.0
        else if (arg .gt. 1.) then
          frac = 0.0
        else
          frac = (1./pie)*acos(arg)
        endif
        tsun = pie*frac
        if (tsun .gt. 0.) then
          coszrsu =  sinc*sind + (cosc*cosd*sin(tsun))/tsun
        else
          coszrsu = 0.0
        endif
        do i=1,plon
         coszrs(i) = coszrsu
        end do
      else                       ! No diurnal averaging
C
C Calday is the calender day for Greenwich, including fraction
C of day; the fraction of the day represents a local time at
C Greenwich; to adjust this to produce a true instantaneous time
C For other longitudes, we must correct for the local time change:
C local time based on the longitude and day of year
C then compute the local cosine solar zenith angle
C
! SJL: calday is negative if leap year
! Apr 11, 2000
        do i=1,plon
! SJL: note that this code segment will not work for non-uniform dx
!         phi       = calday + (real(i-1)/real(plon))
          phi       = abs(calday) + (real(i-1)/real(plon))
          coszrs(i) = sinc*sind - cosc*cosd*cos(2.*pie*phi)
        end do
      end if
C
      return
      end
