#include <params.h>
      subroutine radctl(hbuf    ,clat    ,coslat  ,lat     ,lwup    ,
     $                  pmid    ,pint    ,pmln    ,piln    ,t       ,
     $                  h2ommr  ,cld     ,effcld  ,clwp    ,coszrs  ,
     $                  asdir   ,asdif   ,aldir   ,aldif   ,fsns    ,
     $                  qrs     ,qrl     ,flwds   ,rel     ,rei     ,
     $                  fice    ,sols    ,soll    ,solsd   ,solld   ,
     $                  absnxt  ,abstot  ,emstot  , hattr)
C-----------------------------------------------------------------------
C
C Driver for radiation computation.
C
C Radiation uses cgs units, so conversions must be done from
C model fields to radiation fields.
C
C Calling sequence:
C
C     radinp      Converts units of model fields and computes ozone
C                 mixing ratio for solar scheme
C
C     radcsw      Performs solar computation
C       radalb    Computes surface albedos
C       radded    Computes delta-Eddington solution
C       radclr    Computes diagnostic clear sky fluxes
C
C     radclw      Performs longwave computation
C
C       radtpl    Computes path quantities
C       radems    Computes emissivity
C       radabs    Computes absorptivity
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C Modified:          B. Briegleb, March 1995 to add aerosol
C                    to shortwave code
C
C Reviewed:          J. Kiehl, April 1996
C Reviewed:          B. Briegleb, May 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
#include <ptrrgrid.h>
#include <pagrid.h>
C------------------------------Commons----------------------------------
C SJL     #include <comhst.h>
C-----------------------------------------------------------------------
#include <comtim.h>
C------------------------------Diagnostics------------------------------
#include <diag.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
C SJL
#if (defined R4BYTE)
C 4-byte
      real*4 absnxt(plond,plev,4)      ! Nearest layer absorptivities
      real*4 abstot(plond,plevp,plevp) ! Non-adjacent layer absorptivites
      real*4 emstot(plond,plevp)       ! Total emissivity

      real*4 n2o(plond,plev)         ! nitrous oxide mass mixing ratio
      real*4 ch4(plond,plev)         ! methane mass mixing ratio
      real*4 cfc11(plond,plev)       ! cfc11 mass mixing ratio
      real*4 cfc12(plond,plev)       ! cfc12 mass mixing ratio
#else
      real   absnxt(plond,plev,4)      ! Nearest layer absorptivities
      real   abstot(plond,plevp,plevp) ! Non-adjacent layer absorptivites
      real   emstot(plond,plevp)       ! Total emissivity

      real n2o(plond,plev)         ! nitrous oxide mass mixing ratio
      real ch4(plond,plev)         ! methane mass mixing ratio
      real cfc11(plond,plev)       ! cfc11 mass mixing ratio
      real cfc12(plond,plev)       ! cfc12 mass mixing ratio
#endif

C SJL
      integer hattr(5,pdiag)    ! History field attributes
      real*4  hbuf(*)           ! History buffer
      integer lat               ! Latitude row index
      real lwup(plond)          ! Longwave up flux at surface
      real pmid(plond,plev)     ! Model level pressures
      real pint(plond,plevp)    ! Model interface pressures
      real pmln(plond,plev)     ! Natural log of pmid
      real rel(plond,plev)      ! liquid effective drop size (microns)
      real rei(plond,plev)      ! ice effective drop size (microns)
      real fice(plond,plev)     ! fractional ice content within cloud
      real piln(plond,plevp)    ! Natural log of pint
      real t(plond,plev)        ! Model level temperatures
      real h2ommr(plond,plev)   ! Model level specific humidity
      real cld(plond,plevp)     ! Fractional cloud cover
      real effcld(plond,plevp)  ! Effective fractional cloud cover
      real clwp(plond,plev)     ! Cloud liquid water path
      real coszrs(plond)        ! Cosine solar zenith angle
      real asdir(plond)         ! albedo shortwave direct
      real asdif(plond)         ! albedo shortwave diffuse
      real aldir(plond)         ! albedo longwave direct
      real aldif(plond)         ! albedo longwave diffuse
      real clat                 ! current latitude(radians)
      real coslat               ! cosine latitude
C
C Output solar arguments
C
      real fsns(plond)          ! Surface absorbed solar flux
      real sols(plond)          ! Downward solar rad onto surface (sw direct)
      real soll(plond)          ! Downward solar rad onto surface (lw direct)
      real solsd(plond)         ! Downward solar rad onto surface (sw diffuse)
      real solld(plond)         ! Downward solar rad onto surface (lw diffuse)
      real qrs(plond,plev)      ! Solar heating rate
C
C Output longwave arguments
C
      real qrl(plond,plev)      ! Longwave cooling rate
      real flwds(plond)         ! Surface down longwave flux
C
C---------------------------Local variables-----------------------------
C
      integer i                 ! index

      real solin(plond)         ! Solar incident flux
      real fsds(plond)          ! Flux Shortwave Downwelling Surface
      real fsnt(plond)          ! Net column abs solar flux at model top
      real fsntc(plond)         ! Clear sky total column abs solar flux
      real fsnsc(plond)         ! Clear sky surface abs solar flux
      real flnt(plond)          ! Net outgoing lw flux at model top
      real flns(plond)          ! Srf longwave cooling (up-down) flux
      real flntc(plond)         ! Clear sky lw flux at model top
      real flnsc(plond)         ! Clear sky lw flux at srf (up-down)
      real pbr(plond,plevr)     ! Model mid-level pressures (dynes/cm2)
      real pnm(plond,plevrp)    ! Model interface pressures (dynes/cm2)
      real o3vmr(plond,plevr)   ! Ozone volume mixing ratio
      real o3mmr(plond,plevr)   ! Ozone mass mixing ratio
      real plco2(plond,plevrp)  ! Prs weighted CO2 path
      real plh2o(plond,plevrp)  ! Prs weighted H2O path
      real tclrsf(plond,plevrp) ! Total clear sky fraction  level to space
      real aermmr(plond,plevr)  ! level aerosol mass mixing ratio
      real rh(plond,plevr)      ! level relative humidity (fraction)
      real lwupcgs(plond)       ! Upward longwave flux in cgs units
C JDC ADD
      real pardif(plond)        ! diffuse photosythetically active radiation (0.35-0.70 um)
      real pardir(plond)        ! direct photosythetically active radiation (0.35-0.70 um)

      real eccf                 ! Earth/sun distance factor
C
C Declare local arrays to which model input arrays are interpolated here.
C Current default is none since radiation grid = model grid.
C
C--------------------------------------------------------------------------
C
C Interpolate model input arrays to radiation vertical grid.  Currently this 
C is a do-nothing routine because radiation grid = model grid.
C
C SJL
c     call torgrid(pmid    ,pint    ,pmln    ,piln    ,t       ,
c    $             h2ommr  ,cld     ,effcld  ,clwp    ,
c    $             pmid    ,pint    ,pmln    ,piln    ,t       ,
c    $             h2ommr  ,cld     ,effcld  ,clwp    )
C SJL

C Interpolate ozone volume mixing ratio to model levels

      call radozn(lat     ,pmid    ,o3vmr   )

      call outfld(hattr(1,iO3VMR),o3vmr ,plond, lat, hbuf)
C
C Set latitude dependent radiation input
C
      call radinp(pmid    ,pint    ,h2ommr  ,cld     ,o3vmr   ,
     $            pbr     ,pnm     ,plco2   ,plh2o   ,tclrsf  ,
     $            eccf    ,o3mmr   )
C
C Solar radiation computation
C
      if (dosw) then
C
C Specify aerosol mass mixing ratio
C
        call aermix(pnm     ,aermmr  ,rh      )
        call t_startf('radcsw')
        call radcsw(pnm     ,h2ommr  ,o3mmr   ,aermmr  ,rh      ,
     $              cld     ,clwp    ,rel     ,rei     ,fice    ,
     $              eccf    ,coszrs  ,asdir   ,asdif   ,aldir   ,
     $              aldif   ,solin   ,qrs     ,fsns    ,fsnt    ,
     $              fsds    ,fsnsc   ,fsntc   ,sols    ,soll    ,
     $              solsd   ,solld   ,pardif  ,pardir  )
        call t_stopf('radcsw')
C
C Convert units of shortwave fields needed by rest of model from CGS to MKS
C
        do i=1,plon
          solin(i) = solin(i)*1.e-3
          fsds(i)  = fsds(i)*1.e-3
          fsnt(i)  = fsnt(i) *1.e-3
          fsns(i)  = fsns(i) *1.e-3
          fsntc(i) = fsntc(i)*1.e-3
          fsnsc(i) = fsnsc(i)*1.e-3
        end do
C
C Dump shortwave radiation information to history tape buffer (diagnostics)
C
        call outfld(hattr(1,iSOLIN), solin ,plond,lat,hbuf)
        call outfld(hattr(1,iFSDS), fsds  ,plond,lat,hbuf)
        call outfld(hattr(1,iFSNT), fsnt  ,plond,lat,hbuf)
        call outfld(hattr(1,iFSNS), fsns  ,plond,lat,hbuf)
        call outfld(hattr(1,iFSNTC),fsntc ,plond,lat,hbuf)
        call outfld(hattr(1,iFSNSC),fsnsc ,plond,lat,hbuf)
C JDC ADD
        call outfld(hattr(1,iPARDIF), pardif ,plond,lat,hbuf)
        call outfld(hattr(1,iPARDIR), pardir ,plond,lat,hbuf)
      end if
C
C Longwave radiation computation
C
      if (dolw) then
c
c Convert upward longwave flux units to CGS
c
        do i=1,plon
          lwupcgs(i) = lwup(i)*1000.
        end do
c
c Specify trace gas mixing ratios    
c
        call trcmix(pmid, clat, coslat, n2o, ch4, cfc11, cfc12)
c
        call radclw(lat     ,lwupcgs ,t       ,h2ommr  ,o3vmr   ,
     $              pbr     ,pnm     ,pmln    ,piln    ,plco2   ,
     $              plh2o   ,n2o     ,ch4     ,cfc11   ,cfc12   ,
     $              effcld  ,tclrsf  ,qrl     ,flns    ,flnt    ,
     $              flnsc   ,flntc   ,flwds   ,
     $              absnxt  ,abstot  ,emstot    )

C
C Convert units of longwave fields needed by rest of model from CGS to MKS
C
        do i=1,plon
          flnt(i)  = flnt(i)*1.e-3
          flns(i)  = flns(i)*1.e-3
          flntc(i) = flntc(i)*1.e-3
          flnsc(i) = flnsc(i)*1.e-3
          flwds(i) = flwds(i)*1.e-3
        end do
C
C Dump longwave radiation information to history tape buffer (diagnostics)
C
        call outfld(hattr(1,iFLNT), flnt  ,plond,lat,hbuf)
        call outfld(hattr(1,iFLNTC),flntc ,plond,lat,hbuf)
        call outfld(hattr(1,iFLNS), flns  ,plond,lat,hbuf)
        call outfld(hattr(1,iFLNSC),flnsc ,plond,lat,hbuf)
      end if
C
C Interpolate radiation output arrays to model vertical grid.  Currently this 
C is a do-nothing routine because radiation grid = model grid.

C SJL
c     call fmrgrid(qrs     ,qrl     ,
c    $             qrs     ,qrl     )

      return
      end
 
