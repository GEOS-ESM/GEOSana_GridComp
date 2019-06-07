! 1st: things that are needed either on restart file and/or initial
! dataset: make global even in SPMD case
 
      common/comsrf/asdir   ,asdif   ,aldir   ,aldif   ,lwup    
      common/comsrf/oro     ,ts      ,tssub   ,sicthk  ,snowh   
      common/comsrf/flwds   ,sols    ,soll    ,solsd   ,solld
 
      real asdir(plond,plat)       ! albedo: shortwave, direct
      real asdif(plond,plat)       ! albedo: shortwave, diffuse
      real aldir(plond,plat)       ! albedo: longwave, direct
      real aldif(plond,plat)       ! albedo: longwave, diffuse
      real lwup(plond,plat)        ! longwave up radiative flux
      real oro(plond,plat)         ! land/ocean/sea ice flag
      real ts(plond,plat)          ! sfc temp (merged w/ocean if coupled)
      real tssub(plond,plevmx,plat)! ccm surface/subsurface temperatures 
      real sicthk(plond,plat)      ! ccm sea-ice thickness (m)
      real snowh(plond,plat)       ! ccm snow depth (liquid water)
      real flwds(plond,plat)       ! downward longwave radiation at surface
      real sols(plond,plat)        ! direct beam solar radiation onto srf
      real soll(plond,plat)        ! direct beam solar radiation onto srf
      real solsd(plond,plat)       ! diffuse solar radiation onto srf (sw)
      real solld(plond,plat)       ! diffuse solar radiation onto srf (lw)

      common/comsrf/lhf     ,shf     ,cflx    ,wsx     ,wsy     
      common/comsrf/tref    ,tbot    ,zbot    ,ubot    ,vbot    
      common/comsrf/qbot    ,pbot    ,precc   ,precl   ,thbot   
      common/comsrf/srfrad  ,prcsnw  ,z0m     ,z0h     ,zpd
 
      real lhf(plond,beglat:endlat)    ! latent heat flux
      real shf(plond,beglat:endlat)    ! sensible heat flux
      real cflx(plond,pcnst,beglat:endlat) ! constituent flux (evap)
      real wsx(plond,beglat:endlat)    ! surface u-stress (N)
      real wsy(plond,beglat:endlat)    ! surface v-stress (N)
      real tref(plond,beglat:endlat)   ! ref height surface air temp
! JDC ADDED
      real z0m(plond,beglat:endlat)    !roughness length, momentum (m)
      real z0h(plond,beglat:endlat)    !roughness length, sensible heat (m)
      real zpd(plond,beglat:endlat)    !displacement height (m)
 
! Atmosphere quantities: needed *from* atmosphere
 
      real tbot(plond,beglat:endlat)   ! bottom level temperature 
      real zbot(plond,beglat:endlat)   ! bottom level height above surface
      real ubot(plond,beglat:endlat)   ! bottom level u wind
      real vbot(plond,beglat:endlat)   ! bottom level v wind
      real qbot(plond,beglat:endlat)   ! bottom level specific humidity
      real pbot(plond,beglat:endlat)   ! bottom level pressure
      real precc(plond,beglat:endlat)  ! convective precipitation rate
      real precl(plond,beglat:endlat)  ! large-scale precipitation rate
      real thbot(plond,beglat:endlat)  ! bottom level potential temperature
      real srfrad(plond,beglat:endlat) ! surface net radiative flux
      real prcsnw(plond,beglat:endlat) ! ccm total snow precip
