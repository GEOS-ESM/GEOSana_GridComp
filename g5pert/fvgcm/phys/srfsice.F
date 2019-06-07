#include <params.h>
      subroutine srfsice(oro     ,snowh   ,ubot    ,vbot    ,tbot    ,
     $                   qbot    ,thbot   ,zbot    ,pmidm1  ,srfrad  ,
     $                   tssub   ,qflx    ,taux    ,tauy    ,ts      ,
     $                   shflx   ,lhflx   ,lwup    )
C-----------------------------------------------------------------------
C Compute sea ice to atmosphere surface fluxes; then compute
C sea ice temperature change.
C 
C Temperatures over sea-ice surfaces are specified in 'plevmx' layers of
C fixed thickness and thermal properties.  The forecast temperatures are
C determined from a backward/implicit diffusion calculation using
C linearized sensible/latent heat fluxes. The bottom ocean temperature
C is fixed at -2C, allowing heat flux exchange with underlying ocean.
C 
C Temperature over sea ice is not allowed to exceed melting temperature.
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
      real snowh(plond)         ! Snow depth (liquid water equivalent)
      real ubot(plond)          ! Bottom level u wind
      real vbot(plond)          ! Bottom level v wind
      real tbot(plond)          ! Bottom level temperature
      real qbot(plond)          ! Bottom level specific humidity
      real thbot(plond)         ! Bottom level potential temperature
      real zbot(plond)          ! Bottom level height above surface
      real pmidm1(plond)        ! Bottom level pressure
      real srfrad(plond)        ! Srf solar abs flux plus down longwave
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
      integer indx(plond)       ! Longitude index array (land)
      integer m                 ! constituent index
      integer npts              ! Number of land points
      integer i,ii              ! Longitude indices
      integer k                 ! Sub-surface level index
      integer isrftyp(plond)    ! Integer type for surface
      real tsbsf(plond,plevmx)  ! Non-adjusted srfc/sub-srfc temperatures
      real fnt(plond)           ! Net surface flux for input conditions
      real dfntdt(plond)        ! ts partial derivative of net srf flx
      real ltheat(plond)        ! Latent heat for given sfc conditions
C
C--------------------------Statement Functions--------------------------
C
      logical seaice
      seaice(i) = nint(oro(i)).eq.2
C
C-----------------------------------------------------------------------
C
C Find sea ice points
C
      npts = 0
      do i=1,plon
        if (seaice(i)) then
          npts = npts + 1
          indx(npts) = i
        end if
      end do
      if (npts.eq.0) return
C
C Set latent heat for evaporation sea ice surface
C
      do ii=1,npts
         i = indx(ii)
         isrftyp(i) = 2
         ltheat(i) = latvap + latice
      end do
C
C Compute surface fluxes, derviatives, and exchange coefficiants
C
      call flxsice(indx    ,npts    ,pmidm1  ,ubot    ,vbot    ,
     $             tbot    ,qbot    ,thbot   ,zbot    ,srfrad  ,
     $             tssub(1,1),ltheat,fnt     ,dfntdt  ,shflx   ,
     $             lhflx   ,taux    ,tauy    )
C
C Initialize surface/subsurface temperatures for srftsb
C
      do k=1,plevmx
         do ii=1,npts
            i = indx(ii)
            tsbsf(i,k) = tssub(i,k)
         end do
      end do
C
C Diffusion calculation for temperature
C
      call srftsb(isrftyp ,indx    ,npts    ,fnt     ,dfntdt  ,
     $            snowh   ,tsbsf   )
C 
C Modification to updated surface temperatures
C Reset temperature to melting point.
C
      do ii=1,npts
         i = indx(ii)
         do k=1,plevmx
            tsbsf(i,k) = min(tsbsf(i,k),tmelt)
         end do
      end do
C
C Update surface and sub-surface temperatures
C
      do k=1,plevmx
         do ii=1,npts
            i = indx(ii)
            tssub(i,k) = tsbsf(i,k)
         end do
      end do
      do ii=1,npts
         i = indx(ii)
         ts(i) = tssub(i,1)
         lwup(i) = stebol * ts(i)**4
      end do
C
C Evaluate constituent fluxes
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
 
