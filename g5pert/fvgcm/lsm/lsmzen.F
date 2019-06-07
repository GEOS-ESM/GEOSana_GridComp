#include <preproc.h>
      subroutine lsmzen (npt, calday, lon, lat, cosz)
C
C Modified for generalized orbit
C 19 November 1996    Bruce P. Briegleb
C
 
#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       lsmzen.F
* purpose:           cosine of solar zenith angle
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* cosine solar zenith angle from:
*    o day (1.x to 365.x), where x=0 (e.g. 213.0) denotes 00:00 at greenwich
*    o latitude,  where SH = - and NH = + 
*    o longitude, where WH = - and EH = +

* the solar declination must match that used in the atmospheric model.
* for ccm2, this code matches the ccm2 cosz to within +/- 0.0001.

* this discrepancy between lsm cosz and atm cosz causes a problem.
* lsm cosz may be <= zero (sun below horizon), in which case albedos
* equal zero, but atm cosz may be > zero (sun above horizon), in which
* case atm model needs albedos. no problem if atm model has sun below 
* horizon but lsm has sun above horizon so long as atm solar fluxes 
* equal zero. a possible solution then is to reset points with sun 
* slightly below horizon to slightly above horizon. 

* in practice this error is not very large. e.g., if albedo error is 
* 0.0001 (atm cosz = 0.0001, lsm cosz = 0) absorbed solar radiation 
* error is incident flux * 0.0001. since incident flux is itself 
* multiplied by atm cosz, incident flux is small. hence, error is small.
* in fact the error is smaller than the difference between atm net solar 
* radiation at surface and lsm net solar radiation at surface, which arises
* due to the different surface radiation parameterizations.

* however, reset points as discussed above just in case atm model may
* blow up if albedos equal zero when atm cosz > 0.
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer npt              !number of points to process
      real calday              !calendar day + fraction (1.xx -> 365.xx)
      real lat(npt)            !latitude  (radians): + = NH
      real lon(npt)            !longitude (radians): + = EH

* output
      real cosz(npt)           !cosine zenith angle 
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      real pi                  !pi
      real delta               !solar declination angle  in radians
      real eccf                !earth orbit eccentricity factor
      real sind                !sine   of declination
      real cosd                !cosine of declination
      real phi                 !greenwich calendar day + longitude offset
      real loctim              !local time (hour)
      real hrang               !solar hour angle, 24 hour periodicity (radians)

      integer k                !index
      integer mcsec            !current seconds in day (0, ..., 86400)
! SJL
      real caldayloc
* -----------------------------------------------------------------

      pi = 4.*atan(1.)

* solar declination: 

      call orbit(calday,delta,eccf)

      sind = sin(delta)
      cosd = cos(delta)

      do k = 1, npt

* local time

!
! SJL: modifications for leap years
!
!        mcsec = nint((calday - int(calday)) * 86400.)
!        phi = calday + lon(k)/(2.*pi)
!        loctim = (mcsec + (phi-calday)*86400.) / 3600.

         caldayloc  = abs(calday)
         mcsec = nint((caldayloc - int(caldayloc)) * 86400.)
         loctim = (mcsec + lon(k)/(2.*pi)*86400.) / 3600.
! SJL
 
* hour angle

         hrang = 360./24. * (loctim-12.) * pi/180.

* cosine solar zenith angle. reset points with sun slightly below horizon 
* to slightly above horizon, as discussed in notes.

         cosz(k) = sin(lat(k))*sind + cos(lat(k))*cosd*cos(hrang)
         if (cosz(k) .ge. -0.001 .and. cosz(k).le. 0.) cosz(k) = 0.001

      end do

      return
      end
 
