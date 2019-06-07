#include <params.h>
      subroutine trefoce(indx    ,npts    ,pmidm1  ,ubot    ,vbot    ,
     $                   tbot    ,qbot    ,thbot   ,zbot    ,ts      ,
     $                   tref    ) 

C-----------------------------------------------------------------------
C
C Computes 10m reference height temperature over ocean
C
C ASSUME:
c   1) Neutral 10m drag coeff: 
c         cdn = .0027/U10N + .000142 + .0000764 U10N
c   2) Neutral 10m stanton number: 
c         ctn = .0327 sqrt(cdn), unstable
c         ctn = .0180 sqrt(cdn), stable
c   3) Neutral 10m dalton number:  
c         cen = .0346 sqrt(cdn)
c   4) The saturation humidity of air at T(K): 
c         qsat(T)  (kg/m^3)
C
c NOTE:
c   1) here, tstar = <WT>/U*, and qstar = <WQ>/U*.
c   2) wind speeds should all be above a minimum speed, say 1.0 m/s
c
C---------------------------Code history--------------------------------
C
C From Bil. Large code:  M.Vertenstein, Sep. 1995
C Standardized:          T. Acker, March 1996
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
C-----------------------------------------------------------------------
#include <parpbl.h>
C------------------------------Commons----------------------------------
#include <comcon.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer indx(plond)      ! Longitude index array (land)
      integer npts             ! Number of land points
C
      real pmidm1(plond)       ! Bottom level pressure
      real ubot(plond)         ! Bottom level u wind
      real vbot(plond)         ! Bottom level v wind
      real tbot(plond)         ! Bottom level temperature
      real qbot(plond)         ! Bottom level specific humidity
      real thbot(plond)        ! Bottom level potential temperature
      real zbot(plond)         ! Bottom level height above surface
      real ts(plond)           ! Surface temperature
C
C Output arguments
C
      real tref(plond)         ! 10m reference temperature over ocean
C
C---------------------------Local variables-----------------------------
C
C
      integer i,ii       ! Longitude indices
      real ssq           ! Surface saturation specific humidity
      real ustar         ! ustar
      real tstar         ! tstar
      real qstar         ! qstar
      real u10n          ! neutral 10 m wind speed over ocean
      real vmag          ! Surface wind magnitude
      real thvbot        ! Bottom level virtual potential temperature
      real delt          ! potential T difference (k)
      real delq          ! specific humidity difference (kg/kg)
      real rdn           ! sqrt of neutral exchange coeff (momentum)
      real rhn           ! sqrt of neutral exchange coefficient (heat)
      real ren           ! sqrt of neutral exchange coeff (tracers)
      real rd            ! sqrt of exchange coefficient (momentum)
      real rh            ! sqrt of exchange coefficient (heat)
      real re            ! sqrt of exchange coefficient (tracers)
      real hol           ! reference height (10m) / monin-obukhov length
      real xsq           ! temporary variable
      real xqq           ! temporary variable
      real alz           ! ln(zbot/z10)
      real psimh         ! Stability function at ref lev (momentum)
      real psixh         ! Stab function at ref lev (heat & tracers)
      real stable        ! Stability factor
      real rbot(plond)   ! Density at bottom model level
      real bn            ! exchange coef funct for interpolation
      real bh            ! exchange coef funct for interpolation
      real fac           ! interpolation factor
      real ln0           ! log factor for interpolation
      real ln3           ! log factor for interpolation
C
C--------------------------Statement Functions--------------------------
C
      real psimhu          ! unstable part of psimh
      real psixhu          ! unstable part of psixh
      real qsat            ! saturation humidty of air
      real cdn             ! function: neutral drag coeff at btm level
      real xd              ! dummy argument
      real Tk              ! temperature (K)
      real Umps            ! wind velocity (m/sec)

C
      qsat(Tk) = 640380. / exp(5107.4/Tk)
C
      cdn(Umps) = 0.0027 / Umps + .000142 + .0000764 * Umps
C
      psimhu(xd) =  log((1.+xd*(2.+xd))*(1.+xd*xd)/8.)
     $             - 2.*atan(xd) + 1.571
C
      psixhu(xd) =  2. * log((1. + xd*xd)/2.)
C
C-----------------------------------------------------------------------
C
      do ii=1,npts
         i = indx(ii)
C
C---------------------------------------------------------------
C Set up necessary variables
C---------------------------------------------------------------
C
         rbot(i)= pmidm1(i)/(rair*tbot(i))
         vmag   = max(umin, sqrt(ubot(i)**2+vbot(i)**2))
         thvbot = thbot(i) * (1.0 + zvir*qbot(i))
         ssq    = 0.98 * qsat(ts(i)) / rbot(i)
         delt   = thbot(i) - ts(i)
         delq   = qbot(i) - ssq  
         alz    = log(zbot(i)/zref) 
C
C---------------------------------------------------------------
C First iteration to converge on Z/L and hence the fluxes
C---------------------------------------------------------------
C
C Initial guess for roots of neutral exchange coefficients, 
C assume z/L=0. and u10n is approximated by vmag.
C Stable if (thbot > ts ).
C
         stable = 0.5 + sign(0.5 , delt)
         rdn  = sqrt(cdn(vmag))
         rhn  = (1.-stable) * 0.0327 + stable * 0.018 
         ren  = 0.0346 
C
C Initial guess of ustar,tstar and qstar
C
         ustar = rdn*vmag
         tstar = rhn*delt
         qstar = ren*delq
C
C Compute stability and evaluate all stability functions
C Stable if (thbot > ts or hol > 0 )
C
         hol = xkar *gravit *zbot(i)
     $         * (tstar/thvbot + qstar/(1./zvir+qbot(i)))
     $         / ustar**2
         hol = sign( min(abs(hol),10.), hol )
         stable = 0.5 + sign(0.5 , hol)
         xsq = max(sqrt(abs(1. - 16.*hol)) , 1.)
         xqq = sqrt(xsq)
         psimh = -5.*hol*stable + (1.-stable)*psimhu(xqq)
         psixh = -5.*hol*stable + (1.-stable)*psixhu(xqq)
C
C Shift 10m neutral wind speed using old rdn coefficient
C
         rd = rdn / (1.+rdn/xkar*(alz-psimh))
         u10n = vmag * rd / rdn
C
C Update the neutral transfer coefficients at 10m and neutral stability
C
         rdn = sqrt(cdn(u10n))
         ren = 0.0346
         rhn = (1.-stable) * 0.0327 + stable * 0.018 
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
         xsq = max(sqrt(abs(1. - 16.*hol)) , 1.)
         xqq = sqrt(xsq)
         psimh = -5.*hol*stable + (1.-stable)*psimhu(xqq)
         psixh = -5.*hol*stable + (1.-stable)*psixhu(xqq)
C
C Shift 10m neutral wind speed using old rdn coefficient
C
         rd = rdn / (1.+rdn/xkar*(alz-psimh))
         u10n = vmag * rd / rdn
C
C Update the neutral transfer coefficients at 10m and neutral stability
C
        rdn = sqrt(cdn(u10n))
        ren = 0.0346
        rhn = (1.-stable) * 0.0327 + stable * 0.018 
C
C Shift all coeffs to measurement height and stability
C
        rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
        rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
        re = ren / (1.+ren/xkar*(alz-psixh))
C
C---------------------------------------------------------------
C Following Geleyn(1988), interpolate ts to a fixed height of 
C 10m in the surface layer. 
C---------------------------------------------------------------
C
C Compute function of exchange coefficients. Assume that 
C cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
C 1/sqrt(cn(i))=1/rdn and sqrt(cm(i))/ch(i)=1/rh 
C
         bn = xkar/rdn
         bh = xkar/rh
C
C Interpolation factor for stable and unstable cases
C
         ln0 = log(1.0 + (10.0/zbot(i))*(exp(bn) - 1.0))
         ln3 = log(1.0 + (10.0/zbot(i))*(exp(bn - bh) - 1.0))
C
C Interpolation factor for stable and unstable cases
C
         fac = (ln0 - 10.0/zbot(i)*(bn - bh))/bh * stable
     $       + (ln0 - ln3)/bh * (1.-stable)
         fac = min(max(fac,0.),1.)
C
C Actual interpolation
C
         tref(i) = ts(i) + (tbot(i) - ts(i))*fac
C
      end do
C
      return
C
      end

 
