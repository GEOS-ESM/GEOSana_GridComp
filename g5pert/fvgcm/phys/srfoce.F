#include <params.h>
      subroutine srfoce(oro     ,ubot    ,vbot    ,tbot    ,qbot    ,
     $                  thbot   ,zbot    ,pmidm1  ,tssub   ,qflx    ,
     $                  taux    ,tauy    ,ts      ,shflx   ,lhflx   ,
     $                  lwup    )
C-----------------------------------------------------------------------
C
C Compute ocean to atmosphere sensible/latent heat fluxes and 
C stress components
C 
C---------------------------Code history--------------------------------
C
C Original version:      CCM1
C Standardized:          L. Buja, Jun 1992, Feb 1996
C Split by surface type: L. Bath, Dec 1994
C Reviewed:              B. Briegleb, J. Hack, Aug 1992, March 1996
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
#include <compbl.h>
C-----------------------------------------------------------------------
#include <comtsc.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real oro(plond)           ! Land/ocean/seaice flag
      real ubot(plond)          ! Bottom level u wind
      real vbot(plond)          ! Bottom level v wind
      real tbot(plond)          ! Bottom level temperature
      real qbot(plond)          ! Bottom level specific humidity
      real thbot(plond)         ! Bottom level potential temperature
      real zbot(plond)          ! Bottom level height above surface
      real pmidm1(plond)        ! Bottom level pressure
C
C Input/output arguments
C
      real tssub(plond,plevmx)  ! Surface/sub-surface temperatures
C
C Output arguments
C
      real qflx(plond,pcnst)    ! Constituent flux (kg/m2/s)
      real taux(plond)          ! X surface stress (N/m2)
      real tauy(plond)          ! Y surface stress (N/m2)
      real ts(plond)            ! surface temperature (K)
      real shflx(plond)         ! Surface sensible heat flux (J/m2/s)
      real lhflx(plond)         ! Surface latent   heat flux (J/m2/s)
      real lwup(plond)          ! surface longwave up flux (W/m2)
C
C---------------------------Local variables-----------------------------
C
      integer i,ii              ! Longitude indices
      integer m                 ! constituent index
      integer indx(plond)       ! Longitude index array (land)
      integer npts              ! Number of land points

      real ltheat(plond)        ! Latent heat for given sfc conditions
C
C--------------------------Statement Functions--------------------------
C
      logical ocean
      ocean(i) = nint(oro(i)).eq.0
C
C-----------------------------------------------------------------------
C
C Set up index array of ocean surfaces
C
      npts = 0
      do i=1,plon
        if (ocean(i)) then
          npts = npts + 1
          indx(npts) = i
        end if
      end do
      if (npts.eq.0) return
C
C Determine latent heat
C
      do ii=1,npts
         i = indx(ii)
         ltheat(i) = latvap
      end do
C
C Compute surface fluxes, derviatives, and exchange coefficiants
C
      call flxoce(indx    ,npts    ,pmidm1  ,ubot    ,vbot      ,
     $            tbot    ,qbot    ,thbot   ,zbot    ,tssub(1,1),
     $            ltheat  ,shflx   ,lhflx   ,taux    ,tauy      ,
     $            lwup    )
C
C Update surface temperature
C
      do ii=1,npts
         i = indx(ii)
         ts(i) = tssub(i,1)
      end do
C
C Evaluate contituent fluxes
C
      do ii=1,npts
         i = indx(ii)
         qflx(i,1) = lhflx(i)/ltheat(i)
      end do
c
c Set non-water constituent fluxes to zero
c
      do m=2,pcnst
        do ii=1,npts
          i = indx(ii)
          qflx(i,m) = 0.
        end do
      end do
C
      return
      end
 
