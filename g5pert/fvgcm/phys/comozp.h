#include <fvgcm.h>
      integer pozlon, pozlev, pozlat, poztim
c     parameter (pozlon=1, pozlev=23, pozlat=FVGCM_LAT, poztim=12)
C SJL 03/25/99
      parameter (pozlon=1, pozlev=55, pozlat=FVGCM_LAT, poztim=12) 

      common/ozone/oznbdy(pozlon,pozlev,pozlat,2), ozlat(pozlat),
     &             ozmixm(pozlon,pozlev,plat,2),
     &             ozmix(pozlev,plat), pin(pozlev),
     &             cdayozm, cdayozp, cplos, cplol
      common/ozone/nm, np, oznid, date_oz(poztim), sec_oz(poztim), 
     &             lonsiz, levsiz, latsiz, timesiz, np1

      real oznbdy     ! O3 values on boundary dataset
      real ozlat      ! Latitude array for bdy dataset values
      real ozmixm     ! O3 mixing ratios interp. in latitude
      real ozmix      ! O3 mixing ratios interp. in time
      real pin        ! O3 pressure values (pascals)
      real cdayozm    ! Calendar day for prv. month O3 values read in
      real cdayozp    ! Calendar day for nxt. month O3 values read in
      real cplos      ! Const for ozone path length calculation
      real cplol      ! Const for pressure-weighted o3 path length calc.
      
      integer nm,np   ! Array indices for prv., nxt month ozone data
      integer oznid   ! netcdf id for ozone variable
      integer date_oz ! Date on ozone dataset (YYYYMMDD)
      integer sec_oz  ! seconds of date on ozone dataset (0-86399)
      integer lonsiz  ! size of longitude dimension on ozone dataset
      integer levsiz  ! size of level dimension on ozone dataset
      integer latsiz  ! size of latitude dimension on ozone dataset
      integer timesiz ! size of time dimension on ozone dataset
      integer np1     ! current forward time index of ozone dataset
 
