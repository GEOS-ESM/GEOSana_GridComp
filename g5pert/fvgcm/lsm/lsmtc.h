* ------------------------ code history ---------------------------
* source file:       lsmtc.h
* purpose:           lsmtc common block for time constant variables 
*                    without dynamic memory allocation 
*                    initialized in lsmini.F 
* date last revised: August 1996 (M. Vertenstein)
* author:            M. Vertenstein
* -----------------------------------------------------------------

* surface characteristics - integer values


      integer surf2d(lsmlon,lsmlat) !surface type for lsmlon x lsmlat grid: 0=ocean, >0 = land
      integer numlon(lsmlat)        !number of longitude points for each latitude strip
      integer soic2d(lsmlon,lsmlat) !soil color

      common /lsmtcnd_i/ surf2d, numlon, soic2d

* surface characteristics - real values

      real latixy(lsmlon,lsmlat) !latitude  on lsmlon x lsmlat grid
      real longxy(lsmlon,lsmlat) !longitude on lsmlon x lsmlat grid
      real fland (lsmlon,lsmlat) !fractional land on lsmlon x lsmlat grid
      real sand2d(lsmlon,lsmlat) !percent sand
      real clay2d(lsmlon,lsmlat) !percent clay
      real silt2d(lsmlon,lsmlat) !percent silt
      real pctlak(lsmlon,lsmlat) !percent lake
      real pctwet(lsmlon,lsmlat) !percent wetland

      common /lsmtcnd_r/ latixy, longxy, fland  , sand2d, clay2d, 
     &                   silt2d, pctlak, pctwet

* ------------------------ end lsmtc.h ----------------------------

 
