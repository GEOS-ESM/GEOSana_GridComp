#include <params.h>
      subroutine flxsice(indx    ,npts    ,pmidm1  ,ubot    ,vbot    ,
     $                   tbot    ,qbot    ,thbot   ,zbot    ,srfrad  ,
     $                   ts      ,ltheat  ,fnt     ,dfntdt  ,shf     ,
     $                   lhf     ,taux    ,tauy    )
C-----------------------------------------------------------------------
C
C Compute sea ice to atmosphere surface fluxes of sensible, latent heat
C and stress components
C
C Follows the same basic parameterizations as for ocean surfaces
C
C--------------------------Code History---------------------------------
C
C Original:      Bill Large/M.Vertenstein, Sep. 1995
C Standardized:  L. Buja,     Feb 1996
C Reviewed:      B. Briegleb, March 1996
C
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
#include <parpbl.h>
C------------------------------Commons----------------------------------
#include <comcon.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer indx(plond) ! Longitude index array (land)
      integer npts        ! Number of land points
C                         
      real pmidm1(plond)  ! Bottom level pressure
      real ubot(plond)    ! Bottom level u wind
      real vbot(plond)    ! Bottom level v wind
      real tbot(plond)    ! Bottom level temperature
      real qbot(plond)    ! Bottom level specific humidity
      real thbot(plond)   ! Bottom level potential temperature
      real zbot(plond)    ! Bottom level height above surface
      real srfrad(plond)  ! Solar absorbed plus down longwave flux
      real ts(plond)      ! Surface temperature
      real ltheat(plond)  ! Latent heat for given srf conditions
C                         
C Output arguments        
C                         
      real fnt(plond)     ! Net surface flux for input conditions (W/m2)
      real dfntdt(plond)  ! Net surface flux ts partial derivative (W/m2)
      real shf(plond)     ! Initial sensible heat flux (W/m2)
      real lhf(plond)     ! Initial latent heat flux (W/m2)
      real taux(plond)    ! X surface stress (N/m2)
      real tauy(plond)    ! Y surface stress (N/m2)
C
C---------------------------Local variables-----------------------------
C
      integer i,ii        ! Longitude indices
      real ssq            ! Surface saturation specific humidity
      real ustar          ! ustar
      real tstar          ! tstar
      real qstar          ! qstar
      real vmag           ! Surface wind magnitude
      real thvbot         ! Bottom level virtual potential temperature
      real delt           ! potential T difference (k)
      real delq           ! specific humidity difference (kg/kg)
      real rdn            ! sqrt of neutral exchange coeff (momentum)
      real rhn            ! sqrt of neutral exchg coeff (heat & tracers)
      real ren            ! sqrt of neutral exchange coeff (water)
      real rd             ! sqrt of exchange coefficient (momentum)
      real rh             ! sqrt of exchange coefficient(heat & tracers)
      real re             ! sqrt of exchange coefficient (water) 
      real hol            ! Ref height (10m) / monin-obukhov length
      real xsq            ! Temporary variable
      real xqq            ! Temporary variable
      real alz            ! ln(zbot/10)
      real cp             ! Specific heat of moist air
      real tau            ! Reference height stress
      real psimh          ! Stability function at ref lev (momentum)
      real psixh          ! Stability function at ref lev (heat & tracers) 
      real stable         ! Stability factor
      real rbot(plond)    ! Density at bottom model level
      real dssqdt         ! Derivate of qs wrt surface temperature
      real dshf(plond)    ! Ts partial derivative for sensible heat flux
      real dlhf(plond)    ! Ts partial derivative for latent heat flux
      real lwup(plond)    ! Long wave up surface flux
C
C------------------------------Functions--------------------------------
C
      real psimhu         ! Unstable part of psimh
      real psixhu         ! Unstable part of psixh
      real qsat           ! Saturation humidty of air
      real dqsatdt        ! Derivative of qsat wrt surface temperature
      real xd             ! Dummy argument
      real Tk             ! Temperature (K)
C
C      qsat(Tk)    = 640380. / exp(5107.4/Tk)
C      dqsatdt(Tk) = (5107.4 / Tk**2) * 640380. / exp(5107.4/Tk)
      qsat(Tk)    = 11637800. / exp(5897.8/Tk)
      dqsatdt(Tk) = (5897.8 / Tk**2) * 11637800. / exp(5897.8/Tk)
      psimhu(xd)  = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.)
     $               - 2.*atan(xd) + 1.571
      psixhu(xd)  = 2. * log((1. + xd*xd)/2.)
C
C-----------------------------------------------------------------------
C
      do ii=1,npts
         i = indx(ii)
C
C-----------------------------------------------------------------------
C Determine some necessary variables
C-----------------------------------------------------------------------
C
         rbot(i)= pmidm1(i)/(rair*tbot(i))
         vmag   = max(umin, sqrt(ubot(i)**2+vbot(i)**2))
         thvbot = thbot(i) * (1.0 + zvir*qbot(i))
C         ssq    = 0.98 * qsat(ts(i)) / rbot(i)
C         dssqdt = 0.98 * dqsatdt(ts(i)) / rbot(i) 
         ssq    =  qsat(ts(i)) / rbot(i)
         dssqdt =  dqsatdt(ts(i)) / rbot(i) 
         delt   = thbot(i) - ts(i)
         delq   = qbot(i) - ssq 
         alz    = log(zbot(i)/zref) 
         cp     = cpair*(1. + cpvir*ssq) 
C
C---------------------------------------------------------------
C First iteration to converge on Z/L and hence the fluxes
C---------------------------------------------------------------
C
C Determine roots of neutral exchange coefficients
C
         rdn = xkar/log(zref/zzsice)
         rhn = rdn
         ren = rdn
C
C Determine initial guess of ustar,tstar and qstar
C
         ustar = rdn*vmag
         tstar = rhn*delt
         qstar = ren*delq
C
C Compute stability and evaluate all stability functions
C Stable if (thbot > ts or hol > 0 )
C
         hol = xkar * gravit * zbot(i)
     $         * (tstar/thvbot + qstar/(1./zvir+qbot(i)))
     $         / ustar**2
         hol = sign( min(abs(hol),10.), hol )
         stable = 0.5 + sign(0.5 , hol)
         xsq   = max(sqrt(abs(1. - 16.*hol)) , 1.)
         xqq   = sqrt(xsq)
         psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
         psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
C
C Shift all coeffs to measurement height and stability
C
         rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
         rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
         re = ren / (1.+ren/xkar*(alz-psixh))
C
C Update ustar, tstar, qstar using updated, shifted coeffs 
C
         ustar = rd * vmag 
         tstar = rh * delt 
         qstar = re * delq 
C
C---------------------------------------------------------------
C Second iteration to converge on Z/L and hence the fluxes
C---------------------------------------------------------------
C
C Recompute stability & evaluate all stability functions  
C Stable if (thbot > ts or hol > 0 )
C 
         hol = xkar * gravit * zbot(i)
     $        * (tstar/thvbot + qstar/(1./zvir+qbot(i)))
     $        / ustar**2
         hol = sign( min(abs(hol),10.), hol )
         stable = 0.5 + sign(0.5 , hol)
         xsq   = max(sqrt(abs(1. - 16.*hol)) , 1.)
         xqq   = sqrt(xsq)
         psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
         psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
C
C Shift all coeffs to measurement height and stability
C
         rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
         rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
         re = ren / (1.+ren/xkar*(alz-psixh)) 
C
C Update ustar, tstar, qstar using updated, shifted coeffs 
C
         ustar = rd * vmag 
         tstar = rh * delt 
         qstar = re * delq 
C
C Compute surface stress components
C
         tau     =  rbot(i) * ustar * ustar 
         taux(i) = -tau * ubot(i) / vmag 
         tauy(i) = -tau * vbot(i) / vmag 
C
C Compute heat flux components at current surface temperature
C (Define positive latent and sensible heat as upwards into the atm)
C
         shf(i) = -cp * tau * tstar / ustar 
         lhf(i) = -ltheat(i) * tau * qstar / ustar
         lwup(i) = stebol * ts(i)**4 
C
C Compute net surface flux surface temperature derivative at the current
C surface temperature (ignore the variation of the exchange coefficients
C with temperature).
C
         dshf(i) = cp * rbot(i) * rd*rh * vmag
         dlhf(i) = ltheat(i) * rbot(i) * rd*re * vmag * dssqdt
C
C Compute net surface flux at current surface temperature
C (Define positive net flux as downwards into surface)
C
         fnt(i) = srfrad(i) - lwup(i) - shf(i) - lhf(i)
C
C Compute derivate of net surface flux (ignore changes due to radiation)
C
         dfntdt(i) = -(dshf(i) + dlhf(i)) - stebol * 4.*ts(i)**3 
C
      end do
C
      return
      end
 
