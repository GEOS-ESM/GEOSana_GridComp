      module m_comsrf
      use mod_comm, only: plat,pcnst
      use m_pmgrid, only: plond,plevmx,beglat,endlat

      PUBLIC

      real,allocatable:: asdir(:,:)       ! albedo: shortwave, direct
      real,allocatable:: asdif(:,:)       ! albedo: shortwave, diffuse
      real,allocatable:: aldir(:,:)       ! albedo: longwave, direct
      real,allocatable:: aldif(:,:)       ! albedo: longwave, diffuse
      real,allocatable:: lwup(:,:)        ! longwave up radiative flux
      real,allocatable:: oro(:,:)         ! land/ocean/sea ice flag
      real,allocatable:: ts(:,:)          ! sfc temp (merged w/ocean if coupled)
      real,allocatable:: tssub(:,:,:)     ! ccm surface/subsurface temperatures
      real,allocatable:: sicthk(:,:)      ! ccm sea-ice thickness (m)
      real,allocatable:: snowh(:,:)       ! ccm snow depth (liquid water)
      real,allocatable:: flwds(:,:)       ! downward longwave radiation at surface
      real,allocatable:: sols(:,:)        ! direct beam solar radiation onto srf
      real,allocatable:: soll(:,:)        ! direct beam solar radiation onto srf
      real,allocatable:: solsd(:,:)       ! diffuse solar radiation onto srf (sw)
      real,allocatable:: solld(:,:)       ! diffuse solar radiation onto srf (lw)
                                                                                                                                            
      real,allocatable:: lhf(:,:)    ! latent heat flux
      real,allocatable:: shf(:,:)    ! sensible heat flux
      real,allocatable:: cflx(:,:,:) ! constituent flux (evap)
      real,allocatable:: wsx(:,:)    ! surface u-stress (N)
      real,allocatable:: wsy(:,:)    ! surface v-stress (N)
      real,allocatable:: tref(:,:)   ! ref height surface air temp
! JDC ADDED
      real,allocatable:: z0m(:,:)    !roughness length, momentum (m)
      real,allocatable:: z0h(:,:)    !roughness length, sensible heat (m)
      real,allocatable:: zpd(:,:)    !displacement height (m)
                                                                                                                                            
! Atmosphere quantities: needed *from* atmosphere
                                                                                                                                            
      real,allocatable:: tbot(:,:)   ! bottom level temperature
      real,allocatable:: zbot(:,:)   ! bottom level height above surface
      real,allocatable:: ubot(:,:)   ! bottom level u wind
      real,allocatable:: vbot(:,:)   ! bottom level v wind
      real,allocatable:: qbot(:,:)   ! bottom level specific humidity
      real,allocatable:: pbot(:,:)   ! bottom level pressure
      real,allocatable:: precc(:,:)  ! convective precipitation rate
      real,allocatable:: precl(:,:)  ! large-scale precipitation rate
      real,allocatable:: thbot(:,:)  ! bottom level potential temperature
      real,allocatable:: srfrad(:,:) ! surface net radiative flux
      real,allocatable:: prcsnw(:,:) ! ccm total snow precip

      contains
 
      subroutine comsrf_init 

! 1st: things that are needed either on restart file and/or initial
! dataset: make global even in SPMD case
 
!     common/comsrf/asdir   ,asdif   ,aldir   ,aldif   ,lwup    
!     common/comsrf/oro     ,ts      ,tssub   ,sicthk  ,snowh   
!     common/comsrf/flwds   ,sols    ,soll    ,solsd   ,solld
 
      allocate( asdir(plond,plat) )      ! albedo: shortwave, direct
      allocate( asdif(plond,plat) )      ! albedo: shortwave, diffuse
      allocate( aldir(plond,plat) )      ! albedo: longwave, direct
      allocate( aldif(plond,plat) )      ! albedo: longwave, diffuse
      allocate( lwup(plond,plat) )       ! longwave up radiative flux
      allocate( oro(plond,plat) )        ! land/ocean/sea ice flag
      allocate( ts(plond,plat) )         ! sfc temp (merged w/ocean if coupled)
      allocate( tssub(plond,plevmx,plat))! ccm surface/subsurface temperatures 
      allocate( sicthk(plond,plat))      ! ccm sea-ice thickness (m)
      allocate( snowh(plond,plat) )      ! ccm snow depth (liquid water)
      allocate( flwds(plond,plat) )      ! downward longwave radiation at surface
      allocate( sols(plond,plat)  )      ! direct beam solar radiation onto srf
      allocate( soll(plond,plat)  )      ! direct beam solar radiation onto srf
      allocate( solsd(plond,plat) )      ! diffuse solar radiation onto srf (sw)
      allocate( solld(plond,plat) )      ! diffuse solar radiation onto srf (lw)

!     common/comsrf/lhf     ,shf     ,cflx    ,wsx     ,wsy     
!     common/comsrf/tref    ,tbot    ,zbot    ,ubot    ,vbot    
!     common/comsrf/qbot    ,pbot    ,precc   ,precl   ,thbot   
!     common/comsrf/srfrad  ,prcsnw  ,z0m     ,z0h     ,zpd
 
      allocate( lhf(plond,beglat:endlat) )   ! latent heat flux
      allocate( shf(plond,beglat:endlat) )   ! sensible heat flux
      allocate( cflx(plond,pcnst,beglat:endlat)) ! constituent flux (evap)
      allocate( wsx(plond,beglat:endlat) )   ! surface u-stress (N)
      allocate( wsy(plond,beglat:endlat) )   ! surface v-stress (N)
      allocate( tref(plond,beglat:endlat))   ! ref height surface air temp
! JDC ADDED
      allocate( z0m(plond,beglat:endlat) )   !roughness length, momentum (m)
      allocate( z0h(plond,beglat:endlat) )   !roughness length, sensible heat (m)
      allocate( zpd(plond,beglat:endlat) )   !displacement height (m)
 
! Atmosphere quantities: needed *from* atmosphere
 
      allocate( tbot(plond,beglat:endlat))   ! bottom level temperature 
      allocate( zbot(plond,beglat:endlat))   ! bottom level height above surface
      allocate( ubot(plond,beglat:endlat))   ! bottom level u wind
      allocate( vbot(plond,beglat:endlat))   ! bottom level v wind
      allocate( qbot(plond,beglat:endlat))   ! bottom level specific humidity
      allocate( pbot(plond,beglat:endlat))   ! bottom level pressure
      allocate( precc(plond,beglat:endlat))  ! convective precipitation rate
      allocate( precl(plond,beglat:endlat))  ! large-scale precipitation rate
      allocate( thbot(plond,beglat:endlat))  ! bottom level potential temperature
      allocate( srfrad(plond,beglat:endlat)) ! surface net radiative flux
      allocate( prcsnw(plond,beglat:endlat)) ! ccm total snow precip

      end subroutine comsrf_init 
      end module m_comsrf
