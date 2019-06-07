#include <preproc.h>
      subroutine landtyp (fin      , igeo_i   , nlon_i, nlat_i, igeo_o, 
     $                    numlon_o , mxovr_i2o, fland , ndiag , lsm_o )

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       landtyp.F
* purpose:           generate land types for lsm grid
* date last revised: April 1996
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ parameters -----------------------------
* define grid for output data

      integer nlon_o    !maximum number of output grid longitude points
      integer nlat_o    !number of output grid latitude  points
      parameter( nlon_o = lsmlon, nlat_o = lsmlat)

* -----------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <lsmctl.h>
* ----------------------------------------------------------------------

* ------------------------ input variables ------------------------

* igeo_i and igeo_o, in "mn" format, define the grid types using:
*   o  m = 1: regular grid
*      m = 2: gaussian grid
*   o  n = 1: grid starts at dateline.  western edge ON dateline
*   o  n = 2: grid starts at greenwich. western edge ON greenwich
*   o  n = 3: grid starts at greenwich. is centered  ON greenwich

      character*(*) fin            !input file name
      
      integer igeo_i               !input grid type in "mn" format
      integer nlon_i               !input grid number of longitude points
      integer nlat_i               !input grid number of latitude  points
      integer igeo_o               !output grid type in "mn" format (see above)
      integer numlon_o(nlat_o)     !output grid number of longitude points at each lat
      integer mxovr_i2o            !max number of input cells that overlap output cell
      real fland (nlon_o,nlat_o)   !output grid fractional land 
      integer ndiag                !unit number for diagnostic output  
* -----------------------------------------------------------------

* ------------------------ output variables -----------------------
      integer lsm_o(nlon_o,nlat_o) !lsm surface type on output grid
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer ii                   !longitude index for olson grid
      integer ji                   !latitude  index for olson grid
      integer ji_t                 !transformed latitude index for olson grid
      integer k                    !temporary olson or lsm surface type
      integer miss                 !missing data indicator
      data miss /99999/            
      integer ols2lsm(100)           !lsm surface type for each olson type
      integer olson_o(nlon_i,nlat_i) !olson surface type: original
      integer olson_m(nlon_i,nlat_i) !olson surface type: modified
      integer lsm_i(nlon_i,nlat_i)   !lsm surface type on olson grid

      integer iovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lon index of overlap input cell
      integer jovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lat index of overlap input cell
      real    wovr_i2o(nlon_o,nlat_o,mxovr_i2o) !weight    of overlap input cell

      integer io                  !output grid longitude index 
      integer jo                  !output grid latitude  index 
      integer n                   !loop index

      integer nlsm                !number of lsm surface types
      parameter (nlsm = 28)
      real wst(0:nlsm)            !overlap weights, by surface type
      integer wst_1               !index to largest value in wst vector
      integer wst_2               !index to 2nd largest value in wst vector
      real gast_i(0:nlsm)         !input grid : global area, by surface type
      real gast_o(0:nlsm)         !output grid: global area, by surface type

      integer numlon_i(nlat_i)    !input grid number of longitude points at each lat
      real lon_i(nlon_i+1,nlat_i) !input grid cell longitude, western edge (degrees)
      real lat_i(nlat_i+1)        !input grid cell latitude, south edge    (degrees)
      real area_i(nlon_i,nlat_i)  !input grid cell area
      real fld_i (nlon_i,nlat_i)  !input grid dummy field

      real lon_o(nlon_o+1,nlat_o) !output grid cell longitude, western edge  (degrees)
      real lat_o(nlat_o+1)        !output grid cell latitude, south edge (degree)
      real area_o(nlon_o,nlat_o)  !output grid cell area
      real fld_o (nlon_o,nlat_o)  !output grid dummy field

      integer ii_w  !grid cell to west
      integer ii_e  !grid cell to east
      integer ji_n  !grid cell to north
      integer ji_s  !grid cell to south
      integer nbad
      integer iisv
      integer jisv

      integer getavu             !next available unit number
      integer ndata              !unit number of input data

      character*35 veg(0:nlsm)    !name of each surface type
      save veg

      data veg( 0)/'ocean                              '/
      data veg( 1)/'land ice                           '/
      data veg( 2)/'desert                             '/
      data veg( 3)/'cool needleleaf evergreen forest   '/
      data veg( 4)/'cool needleleaf deciduous forest   '/
      data veg( 5)/'cool broadleaf deciduous forest    '/
      data veg( 6)/'cool mixed ne+bd forest            '/
      data veg( 7)/'warm needleleaf evergreen forest   '/
      data veg( 8)/'warm broadleaf deciduous forest    '/
      data veg( 9)/'warm mixed ne+bd forest            '/
      data veg(10)/'tropical broadleaf evergreen forest'/
      data veg(11)/'tropical seasonal deciduous forest '/
      data veg(12)/'savanna                            '/
      data veg(13)/'evergreen forest tundra            '/
      data veg(14)/'deciduous forest tundra            '/
      data veg(15)/'cool forest crop                   '/
      data veg(16)/'warm forest crop                   '/
      data veg(17)/'cool grassland                     '/
      data veg(18)/'warm grassland                     '/
      data veg(19)/'tundra                             '/
      data veg(20)/'evergreen shrubland                '/
      data veg(21)/'deciduous shrubland                '/
      data veg(22)/'semi-desert                        '/
      data veg(23)/'cool irrigated crop                '/
      data veg(24)/'cool crop                          '/
      data veg(25)/'warm irrigated crop                '/
      data veg(26)/'warm crop                          '/
      data veg(27)/'forest wetland                     '/
      data veg(28)/'non-forest wetland                 '/
      integer lencs              !length of character string
      character*256 fulpath      !full pathname of input dataset 

* -----------------------------------------------------------------

      write (6,*)
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*) 'Making Surface Types'
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*)

* -----------------------------------------------------------------
* define the lsm surface types: 0 to nlsm
* -----------------------------------------------------------------


* -----------------------------------------------------------------
* read in Olson's input data
* -----------------------------------------------------------------

* olson data are 1/2 x 1/2 degree, stored in latitude bands, 
* from north to south. in a given latitude band, data begin 
* at the date line (180w) and proceed eastward. so first data  
* point (olson(1,1)) is a box centered at 89.75n, 179.75w. 
*
*   90.0N  ---------------------  
*          |         |         |
*          |    x    |    x    |
*          |  (1,1)  |  (2,1)  |
*          |         |         |
*   89.5N  --------------------- 
*        180.0W    179.5W    179.0W

      
      fulpath = srfpath(1:lencs(srfpath)) //'/'// fin(1:lencs(fin))
      call getfil (fulpath, fin)
      ndata = getavu()
      call opnfil (fin, ndata, 'f')
      read (ndata,*) olson_o
      call relavu (ndata) 

* modify data (similar to BATS data set code) and
* transform lat index FROM np=1 and sp=nlat_i TO sp=1 and np=nlat_i
 
      do ji = 1, nlat_i
         ji_t = nlat_i - ji + 1  

         do ii = 1, nlon_i

            k = olson_o(ii,ji)

* there are several values olson=2, olson=6, olson=8 in the data set that
* are not defined. use neighboring cells in this order: west, east, north, south

            ii_w = max(     1,ii-1)
            ii_e = min(nlon_i,ii+1)
            ji_n = min(     1,ji-1)
            ji_s = max(nlat_i,ji+1)

            if (k.eq.2 .or. k.eq.6 .or. k.eq.8) k = olson_o(ii_w,ji  )
            if (k.eq.2 .or. k.eq.6 .or. k.eq.8) k = olson_o(ii_e,ji  )
            if (k.eq.2 .or. k.eq.6 .or. k.eq.8) k = olson_o(ii  ,ji_n)
            if (k.eq.2 .or. k.eq.6 .or. k.eq.8) k = olson_o(ii  ,ji_s)

* split antarctica (17) into polar desert (69) and ice (70)
 
            if (k .eq. 17) then
               if (ji .le. 313) then
                  k = 69
               else
                  k = 70
               end if
            end if

* olson=61 (eastern south taiga) will be classified as needleleaf
* deciduous tree. change olson=61 to olson=20 (main taiga =
* needleleaf evergreen tree) based on longitude

            if (k.eq.61 .and. ii.le.576) k = 20 

* olson=61 (eastern south taiga) will be classified needleleaf deciduous tree. 
* create additional needleleaf deciduous tree from olson=21 (main taiga) and
* olson=60 (southern taiga) based on longitude

            if (k.eq.21 .and. ii.ge.555) k = 61  
            if (k.eq.60 .and. ii.ge.582) k = 61   

* change olson=26 (warm mixed) to broad-leaved humid forest based on latitude

            if (k.eq.26 .and. ji.ge.113) k = 29

* split forest tundra (62, 63) into needleleaf evergreen forest tundra (62) 
* and needleleaf deciduous forest tundra (63) based on longitude

            if (k.eq.63) k = 62
            if (k.eq.62 .and. ii.ge.490) k = 63   

* assign modified olson type to transformed grid

            olson_m(ii,ji_t) = k

         end do
      end do

* error checks

      nbad = 0
      do ji = 1, nlat_i
         do ii = 1, nlon_i
            if (olson_m(ii,ji).gt.100 .or. olson_m(ii,ji).lt.0) then
               nbad = nbad + 1
               iisv = ii
               jisv = ji
            end if
         end do
      end do
      if (nbad.gt.0) then
         write (6,*) nbad,' landtyp errors. Last one is olson surface',
     &' type = ',olson_m(iisv,jisv),' is undef for lon,lat = ',iisv,jisv
         call stoprun
      end if

* -----------------------------------------------------------------
* convert olson surface types to lsm surface types
* -----------------------------------------------------------------

* initialize all cells on input grid to missing value 

      do ji = 1, nlat_i                                                  
         do ii = 1, nlon_i                                                  
            lsm_i(ii,ji) = miss 
         end do
      end do

* assign each of the olson surface types to an lsm surface type.
* this mapping from olson to lsm is based on the BATS data set code.
* note: ols2lsm(i) = olson type i

      do ii=1,19
         ols2lsm(ii) = miss
      end do
      ols2lsm(20) = 3                                                     
      ols2lsm(21) = 3                                                     
      ols2lsm(22) = 3                                                     
      ols2lsm(23) = 6                                                    
      ols2lsm(24) = 8                                                    
      ols2lsm(25) = 9                                                     
      ols2lsm(26) = 9                                                     
      ols2lsm(27) = 7             
      ols2lsm(28) = 10                                                     
      ols2lsm(29) = 10                                                     
      ols2lsm(30) = 24                                                     
      ols2lsm(31) = 26                                                     
      ols2lsm(32) = 12                                                     
      ols2lsm(33) = 10                                                     
      ols2lsm(34) = miss                                                 
      ols2lsm(35) = miss                                                 
      ols2lsm(36) = 28                                                    
      ols2lsm(37) = 25                                                    
      ols2lsm(38) = 23                                                    
      ols2lsm(39) = 23                                                    
      ols2lsm(40) = 17                                                     
      ols2lsm(41) = 18                                                     
      ols2lsm(42) = 17                                                     
      ols2lsm(43) = 12                                                     
      ols2lsm(44) = 28                                                    
      ols2lsm(45) = 28                                                    
      ols2lsm(46) = 20                                                    
      ols2lsm(47) = 20                                                    
      ols2lsm(48) = 20         
      ols2lsm(49) = 22                                                    
      ols2lsm(50) = 2                                                     
      ols2lsm(51) = 22                                                    
      ols2lsm(52) = 22         
      ols2lsm(53) = 19                                                     
      ols2lsm(54) = 19                                                     
      ols2lsm(55) = 15          
      ols2lsm(56) = 16         
      ols2lsm(57) = 15         
      ols2lsm(58) = 16          
      ols2lsm(59) = 21                                                    
      ols2lsm(60) = 6                                                    
      ols2lsm(61) = 4                                                     
      ols2lsm(62) = 13          
      ols2lsm(63) = 14                                                    
      ols2lsm(64) = 20                                                    
      ols2lsm(65) = 0         
      ols2lsm(66) = 0         
      ols2lsm(67) = 0         
      ols2lsm(68) = 0         
      ols2lsm(69) = 2                                                     
      ols2lsm(70) = 1                                                   
      ols2lsm(71) = 22                                                    
      ols2lsm(72) = 27                                                    
      ols2lsm(73) = 0         
      do ii=74,100
         ols2lsm(ii) = miss
      end do

* transform olson surface types to lsm surface types

      nbad = 0
      do ji = 1 , nlat_i                                                 
         do ii = 1, nlon_i                                                 

            if (olson_m(ii,ji) .eq. 0) then
               lsm_i(ii,ji) = 0
            else
               lsm_i(ii,ji) = ols2lsm(olson_m(ii,ji))
            end if

            if (lsm_i(ii,ji).gt.nlsm .or. lsm_i(ii,ji).lt.0) then
               nbad = nbad + 1
               iisv = ii
               jisv = ji
            end if

         end do
      end do
      if (nbad.gt.0) then
         write (6,*) nbad,' landtyp errors: Last one is lsm surface',
     &' type = ',lsm_i(iisv,jisv),' is undef for lon,lat = ',iisv,jisv
         call stoprun
      end if

* -----------------------------------------------------------------
* map input grid to output grid
* -----------------------------------------------------------------

      do ji=1,nlat_i
         numlon_i(ji) = nlon_i
      end do

      call map_i (igeo_i  , nlon_i  , nlat_i  , numlon_i ,
     &            igeo_o  , nlon_o  , nlat_o  , numlon_o ,
     &            iovr_i2o, jovr_i2o, wovr_i2o, mxovr_i2o,
     &            lon_i   , lat_i   , area_i  , lon_o    ,
     &            lat_o   , area_o  , fld_o   , fld_i    )

* -----------------------------------------------------------------
* process each cell on output grid
* -----------------------------------------------------------------

      do jo = 1, nlat_o
         do io = 1, numlon_o(jo)

* sum overlap weights by surface type

            do k = 0, nlsm
               wst(k) = 0.
            end do

            do n = 1, mxovr_i2o        !overlap cell index
               ii = iovr_i2o(io,jo,n)  !lon index (input grid) of overlap cell
               ji = jovr_i2o(io,jo,n)  !lat index (input grid) of overlap cell
               k = lsm_i(ii,ji)        !surface type (input grid)
               wst(k) = wst(k) + wovr_i2o(io,jo,n)
            end do

* rank non-zero weights by surface type. wst_1 is the most extensive
* surface type. wst_2 is the second most extensive surface type

            call rank (nlsm, wst, miss, wst_1, wst_2)

* set surface type as:

* o if fractional land = 0: cell = ocean

* o if fractional land > 0: cell = land
*   a. use most frequent surface type based on area of overlap unless
*      this is ocean. in this case, input grid says ocean but output
*      grid wants land
*   b. so use next most extensive surface type so long as is not ocean
*   c. if this is ocean or if there is none, use swamp/marsh surface type

            if (fland(io,jo) .eq. 0.) then        !ocean
               lsm_o(io,jo) = 0
            else                                  !land
               if (wst_1 .ne. 0) then
                  lsm_o(io,jo) = wst_1 
               else
                  if (wst_2 .eq. 0 .or. wst_2 .eq. miss) then
                     lsm_o(io,jo) = 28
                  else
                     lsm_o(io,jo) = wst_2
                  end if
               end if
            end if

* error checks

            if (lsm_o(io,jo).lt.0 .or. lsm_o(io,jo).gt.nlsm) then
               write (6,*) 'landtyp error: lsm surface type = ',
     &         lsm_o(io,jo),' is not valid for lon,lat = ',io,jo
               stop
            end if

            if (fland(io,jo).gt.0. .and. lsm_o(io,jo).eq.0) then
               write (6,*) 'landtyp error: fland = land. lsm = ',
     &         lsm_o(io,jo),' for lon,lat = ',io,jo
               stop
            end if

            if (fland(io,jo).eq.0. .and. lsm_o(io,jo).ne.0) then
               write (6,*) 'landtyp error: fland = ocean. lsm = ',
     &         lsm_o(io,jo),' for lon,lat = ',io,jo
               stop
            end if

         end do
      end do

* -----------------------------------------------------------------
* compare global area of each lsm surface type on input and output grids
* -----------------------------------------------------------------

* input grid

      do k = 0, nlsm
         gast_i(k) = 0.
      end do

      do ji = 1, nlat_i
         do ii = 1, nlon_i
            k = lsm_i(ii,ji)
            gast_i(k) = gast_i(k) + area_i(ii,ji)
         end do
      end do

* output grid

      do k = 0, nlsm
         gast_o(k) = 0.
      end do

      do jo = 1, nlat_o
         do io = 1, numlon_o(jo)
            k = lsm_o(io,jo)
            gast_o(k) = gast_o(k) + area_o(io,jo)
         end do
      end do

* area comparison

      write (ndiag,*)
      write (ndiag,'(1x,70a1)') ('=',k=1,70)
      write (ndiag,*) 'Surface Type Output'
      write (ndiag,'(1x,70a1)') ('=',k=1,70)

      write (ndiag,*)
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,1001)
1001  format (1x,'vegetation type',20x,' input grid area',
     &         ' output grid area',/
     &         1x,33x,'     10**6 km**2','      10**6 km**2')
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,*)

      do k = 0, nlsm
         write (ndiag,1002) veg(k),gast_i(k)*1.e-6,gast_o(k)*1.e-6
1002     format (1x,a35,f16.3,f17.3)
      end do

      write (ndiag,*)
      write (ndiag,*) 'For reference the area on the output grid of a',
     & ' cell near the equator is: '
      write (ndiag,'(f10.3,a14)')area_o(1,nlat_o/2)*1.e-06,
     &     ' x 10**6 km**2'
      write (ndiag,*)

      write (6,*) 'successfully made surface types'
      write (6,*)

      return
      end
 
