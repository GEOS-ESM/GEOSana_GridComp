#include <preproc.h>
      subroutine soilcol (fin     , igeo_i   , nlon_i, nlat_i, igeo_o,
     $                    numlon_o, mxovr_i2o, lon_c , lat_c , veg_o ,
     $                    ndiag   , lsm_o    )

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       soilcol.F
* purpose:           generate soil colors for lsm grid
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

      character*(*) fin             !input file name
      integer igeo_i                !input grid type in "mn" format
      integer nlon_i                !input grid number of longitude points
      integer nlat_i                !input grid number of latitude  points
      integer igeo_o                !output grid type in "mn" format (see above)
      integer numlon_o(nlat_o)      !output grid number of longitude points at each lat
      integer mxovr_i2o             !max number of input cells that overlap output cell
      real lon_c(nlon_o,nlat_o)     !output grid cell center longitudes (degrees)
      real lat_c(nlat_o)            !output grid cell center latitudes  (degrees)
      integer veg_o(nlon_o,nlat_o)  !output grid land type 
      integer ndiag                 !unit number for diagnostic output  
* -----------------------------------------------------------------

* ------------------------ output variables -----------------------
      integer lsm_o(nlon_o,nlat_o)  !lsm color class on output grid
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer ii                  !input grid longitude index 
      integer ji                  !input grid latitude  index
      integer k                   !temporary bats or lsm color type
      integer bats(nlon_i,nlat_i) !input grid bats soil color

      integer iovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lon indx of overlap input cell
      integer jovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lat indx of overlap input cell
      real    wovr_i2o(nlon_o,nlat_o,mxovr_i2o) !weight   of overlap input cell

      integer io                  !longitude index for lsm grid
      integer jo                  !latitude  index for lsm grid
      integer n                   !loop index

      integer nlsm                !number of lsm color classes
      parameter (nlsm = 9)
      character*35 col(0:nlsm)    !name of each color
      real wst(0:nlsm)            !overlap weights, by color type
      integer wst_1               !index to largest value in wst vector
      integer wst_2               !index to 2nd largest value in wst vector
      real gast_i(0:nlsm)         !input grid: global area, by color type
      real gast_o(0:nlsm)         !output grid: global area, by color type

      integer numlon_i(nlat_i)    !input grid number of lon points at each lat
      real lon_i(nlon_i+1,nlat_i) !input grid cell longitude, w. edge (degrees)
      real lat_i(nlat_i+1)        !input grid cell latitude, s. edge (degrees)
      real area_i(nlon_i,nlat_i)  !input grid cell area
      real fld_i (nlon_i,nlat_i)  !input grid dummy field

      real lon_o(nlon_o+1,nlat_o) !output grid cell longitude, w. edge (deg)
      real lat_o(nlat_o+1)        !output grid cell latitude, s. edge (deg)
      real area_o(nlon_o,nlat_o)  !output grid cell area
      real fld_o (nlon_o,nlat_o)  !output grid dummy field

      integer miss                !missing data indicator
      data miss /99999/

      integer getavu              !next available unit number
      integer ndata               !unit number of input data
      integer lencs               !length of character string
      character*256 fulpath       !full pathname of input dataset 

* -----------------------------------------------------------------

      write (6,*)
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*) 'Making Soil Colors'
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*)

* -----------------------------------------------------------------
* define the lsm color classes: 0 to nlsm
* -----------------------------------------------------------------

      col(0) = 'no soil                            '
      col(1) = 'class 1: light                     '
      col(2) = 'class 2:                           '
      col(3) = 'class 3:                           '
      col(4) = 'class 4:                           '
      col(5) = 'class 5:                           '
      col(6) = 'class 6:                           '
      col(7) = 'class 7:                           '
      col(8) = 'class 8: dark                      '
      col(9) = 'class 9: very light North Africa   '

* -----------------------------------------------------------------
* read in BATS input data
* -----------------------------------------------------------------

* BATS data are on t42 gaussian grid, approximately 2.8 x 2.8 degrees,
* stored in latitude bands, from south to north. in a given latitude band, 
* data begin at greenwich, centered on greenwich, and proceed eastward. 

      fulpath = srfpath(1:lencs(srfpath)) //'/'// fin(1:lencs(fin))
      call getfil (fulpath, fin)
      ndata = getavu()
      call opnfil (fin, ndata, 'f')
      do ji = 1, nlat_i
         do ii = 1, nlon_i
            read (ndata,*) bats(ii,ji)
            if (bats(ii,ji).lt.0 .or. bats(ii,ji).gt.8) then
               write (6,*) 'soilcol error: bats soil color = ',
     &         bats(ii,ji),' is not valid for lon,lat = ',ii,ji
               stop
            end if
         end do
      end do
      call relavu (ndata) 

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

* sum overlap weights by color type

            do k = 0, nlsm
               wst(k) = 0.
            end do

            do n = 1, mxovr_i2o         !overlap cell index
               ii = iovr_i2o(io,jo,n)  !lon index (input grid) of overlap cell
               ji = jovr_i2o(io,jo,n)  !lat index (input grid) of overlap cell
               k = bats(ii,ji)         !color type (input grid)
               wst(k) = wst(k) + wovr_i2o(io,jo,n)
            end do

* rank non-zero weights by color type. wst_1 is the most extensive
* color type. wst_2 is the second most extensive color type

            call rank (nlsm, wst, miss, wst_1, wst_2)
            lsm_o(io,jo) = wst_1

         end do
      end do

* -----------------------------------------------------------------
* corrections and error checks
* -----------------------------------------------------------------

      do jo = 1, nlat_o
         do io = 1, numlon_o(jo)

* if land but no color, set color to 4

            if (veg_o(io,jo).gt.0 .and. lsm_o(io,jo).eq.0)
     &      lsm_o(io,jo) = 4

* set ocean colors to zero

            if (veg_o(io,jo) .eq. 0) lsm_o(io,jo) = 0

* set glacier colors to zero

            if (veg_o(io,jo) .eq. 1) lsm_o(io,jo) = 0

* increase desert and semi-desert albedos for north africa and arabian 
* peninsula. this region is defined by the lat and lon at the center
* of the grid cell. 

            if (lat_c(jo).ge.15. .and. lat_c(jo).le.35.) then
               if (lon_c(io,jo).ge.-17. .and. lon_c(io,jo).le.57.) then
                  if (veg_o(io,jo).eq. 2) lsm_o(io,jo) = 9
                  if (veg_o(io,jo).eq.22) lsm_o(io,jo) = 9
               end if
            end if

* error checks

            if (lsm_o(io,jo).lt.0 .or. lsm_o(io,jo).gt.nlsm) then
               write (6,*) 'soilcol error: lsm soil color = ',
     &         lsm_o(io,jo),' is not valid for lon,lat = ',io,jo
               stop
            end if

         end do
      end do

* -----------------------------------------------------------------
* compare global area of each soil color on input and output grids
* -----------------------------------------------------------------

* input grid

      do k = 0, nlsm
         gast_i(k) = 0.
      end do

      do ji = 1, nlat_i
         do ii = 1, nlon_i
            k = bats(ii,ji)
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
      write (ndiag,*) 'Soil Color Output'
      write (ndiag,'(1x,70a1)') ('=',k=1,70)

      write (ndiag,*)
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,1001)
1001  format (1x,'soil color type',20x,' input grid area',
     &         ' output grid area',/
     &         1x,33x,'     10**6 km**2','      10**6 km**2')
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,*)

      do k = 0, nlsm
         write (ndiag,1002) col(k),gast_i(k)*1.e-6,gast_o(k)*1.e-6
1002     format (1x,a35,f16.3,f17.3)
      end do

      write (ndiag,*)
      write (ndiag,*) 'For reference the area on the output grid of a',
     & ' cell near the equator is: '
      write (ndiag,'(f10.3,a14)')area_o(1,nlat_o/2)*1.e-06,
     &     ' x 10**6 km**2'
      write (ndiag,*)

      write (6,*) 'successfully made soil color types'
      write (6,*)

      return
      end
 
