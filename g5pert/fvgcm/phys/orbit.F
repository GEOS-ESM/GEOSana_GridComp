      subroutine orbit(calday,delta,eccf)
C
C Compute earth/orbit parameters using formula from ccm2
C
C---------------------------Code history--------------------------------
C
C Original version:     B. Briegleb, December 1996
C Modified:             B. Briegleb, January 1997
C (to correct 'calday' to 'calday-1.')
C
C-----------------------------------------------------------------------
C
      implicit none
C
C------------------------------Arguments--------------------------------
C Input arguments
C
C SJL: calday negative if leap year [-1, -367]; [1,366] otherwise
      real calday   ! Calendar day, including fraction
C
C Output arguments
C
      real delta,   ! Solar declination angle in radians
     $     eccf     ! Earth-sun distance factor ( i.e. (1/r)**2 )
C
C---------------------------Local variables-----------------------------
C
      real dayspy     ! Number of days per 1 year
      real pie        ! 3.14.....
      real theta      ! earth orbit mean angle
C
C-----------------------------------------------------------------------
C
C Compute earth orbit angle
C
      dayspy = 365.
      pie    = 4.*atan(1.)

      if( calday .lt. 0.) then
          theta  = 2.*pie*(abs(calday)-1.) / 366.
      else
          theta  = 2.*pie*(calday-1.) / dayspy
      endif
C
C Solar declination in radians:
C
      delta = .006918 - .399912*cos(theta) + .070257*sin(theta) -
     $        .006758*cos(2.*theta) + .000907*sin(2.*theta) -
     $        .002697*cos(3.*theta) + .001480*sin(3.*theta)
C
C Compute eccentricity factor (sun-earth distance factor)
C
      eccf  = 1.000110 + .034221*cos(theta) + .001280*sin(theta) +
     $         .000719*cos(2.*theta) + .000077*sin(2.*theta)
C
      return
C
      end
