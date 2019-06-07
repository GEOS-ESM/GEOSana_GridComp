#include <preproc.h>
      subroutine basinmat (fin       , igeo_i   , nlon_i    , nlat_i ,
     $                     igeo_o    , numlon_o , mxovr_i2o , ndiag  ,
     $                     drnbasn_o , drnfrac_o, drnarea_o )

#include <vartyp.h>
#include <lsmpar.h>
#include <basinpar.h>

* ------------------------ code history ---------------------------
* source file:       basinmat.F
* purpose:           generate basin drainage matrix 
* date last revised: March 1997
* author:            Mariana Vertenstein
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ parameters -----------------------------
* define grid for output data
      integer nlon_o           !maximum number of output grid longitude points
      integer nlat_o           !number of output grid latitude  points
      parameter( nlon_o = lsmlon, nlat_o = lsmlat)
* -----------------------------------------------------------------

* ------------------------ common block variables --------------------
#include <lsmctl.h>
#include <basin.h>
#include <netcdf.inc>
* -----------------------------------------------------------------

* ------------------------ input variables ------------------------
* igeo_i and igeo_o, in "mn" format, define the grid types using:
*   o  m = 1: regular grid
*      m = 2: gaussian grid
*   o  n = 1: grid starts at dateline.  western edge ON dateline
*   o  n = 2: grid starts at greenwich. western edge ON greenwich
*   o  n = 3: grid starts at greenwich. is centered  ON greenwich

      character*(*) fin          !input file name
      integer igeo_i             !input grid type in "mn" format
      integer nlon_i             !input grid number of longitude points
      integer nlat_i             !input grid number of latitude  points
      integer igeo_o             !output grid type in "mn" format (see above)
      integer mxovr_i2o          !max number of input cells that overlap output cell
      integer numlon_o(nlat_o)   !output grid number of longitude points at each lat
      integer ndiag              !unit number for diagnostic output  
* -----------------------------------------------------------------

* ------------------------ output variables -----------------------
      integer drnbasn_o(ndrnmax,nlon_o,nlat_o) !output drainage basin matrix
      real drnfrac_o   (ndrnmax,nlon_o,nlat_o) !output fractional drainage matrix
      real drnarea_o           (nlon_o,nlat_o) !output grid cell area
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
    
      integer ii               !input grid longitude index 
      integer ji               !input grid latitude  index
      integer i,k,nb           !loop indices

      integer iovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lon indx of overlap input cell
      integer jovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lat indx of overlap input cell
      real    wovr_i2o(nlon_o,nlat_o,mxovr_i2o) !weight   of overlap input cell

      integer io               !longitude index for lsm grid
      integer jo               !latitude  index for lsm grid
      integer n                !loop index

      real wgt(ndrn)           !overlap weights, by drainage basin type
      real gast_i(ndrn)        !input  grid: global area, by drainage basin type
      real gast_o(ndrn)        !output grid: global area, by drainage basin type
      real wt                  !cumulative weight on output area
      real norm                !normalization factor

      integer drnbasn_i(nlon_i,nlat_i) !input grid drainage basin data

      integer numlon_i(nlat_i)    !input grid number of lon points at each lat
      real area_i(nlon_i,nlat_i)  !input grid cell area
      real lon_i(nlon_i+1,nlat_i) !input grid cell longitude, w. edge (degrees)
      real lat_i(nlat_i+1)        !input grid cell latitude , s. edge (degrees)
      real fld_i (nlon_i,nlat_i)  !input grid dummy field

      real lon_o(nlon_o+1,nlat_o) !output grid cell longitude, w. edge (deg)
      real lat_o(nlat_o+1)        !output grid cell latitude , s. edge (deg)
      real fld_o(nlon_o,nlat_o)   !output grid dummy field

      integer miss                !missing data indicator
      data miss /99999/

      real eps                    !permissable error for normalization
      data eps  /1.e-6/

      integer ncid                !netcdf id
      integer loni_id             !netcdf input longitude dimension id
      integer lati_id             !netcdf input latitude dimension id
      integer matr_id             !netcdf input basin matrix id
      integer nlon_ncd            !netcdf number of longitudes 
      integer nlat_ncd            !netcdf number of latitudes 
      integer rcode               !netcdf error code
      integer start(2)            !netcdf hyperslab corner
      integer count(2)            !netcdf hyperslab edge lengths (from hyperslab corner)
      character*80 lon_nam        !netcdf longitude coordinate variable name
      character*80 lat_nam        !netcdf latitude coordinate variable name
      integer lencs               !length of character string
      character*256 fulpath       !full pathname of input dataset
      
* -----------------------------------------------------------------

      write (6,*)
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*) 'Making Drainage Basin Matrix'
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*)

* -----------------------------------------------------------------
* read in drainage basin input data
* -----------------------------------------------------------------

* obtain and open netcdf file

      fulpath = srfpath(1:lencs(srfpath)) //'/'// fin(1:lencs(fin))
      call getfil (fulpath, fin)
      ncid = ncopn(fin, NCNOWRIT, rcode)

* get longitude, latitude id's. Check that file has same dimensions
* as those contained in parameter declaration.

      loni_id = ncdid(ncid, 'longitude', rcode) 
      lati_id = ncdid(ncid, 'latitude' , rcode)                 

      call ncdinq (ncid, loni_id, lon_nam, nlon_ncd, rcode) 
      call ncdinq (ncid, lati_id, lat_nam, nlat_ncd, rcode) 

      if (nlon_ncd .ne. nlon_i) then
         write(6,*)'(DRAINBAS): nlon_i and nlon_ncd do not match'
         write(6,*)' nlon_i = ',nlon_i,' nlon_ncd = ',nlon_ncd
         call stoprun
      endif

      if (nlat_ncd .ne. nlat_i) then
         write(6,*)'(DRAINBAS): nlat_i and nlat do not match'
         write(6,*)' nlat_i = ',nlat_i,' nlat_ncd = ',nlat_ncd
         call stoprun
      endif

* obtain drainage basin matrix on high resolution grid

      start(1) = 1
      count(1) = nlon_i
      start(2) = 1
      count(2) = nlat_i
      matr_id  = ncvid (ncid, 'DRNMAT' , rcode)
      call ncvgt(ncid, matr_id , start, count, drnbasn_i , rcode)

* close netcf file

      call ncclos(ncid, rcode)

* -----------------------------------------------------------------
* map input grid to output grid
* -----------------------------------------------------------------

      do ji=1,nlat_i
         numlon_i(ji) = nlon_i
      end do

      call map_i (igeo_i  , nlon_i   , nlat_i  , numlon_i ,
     &            igeo_o  , nlon_o   , nlat_o  , numlon_o ,
     &            iovr_i2o, jovr_i2o , wovr_i2o, mxovr_i2o,
     &            lon_i   , lat_i    , area_i  , lon_o    ,
     &            lat_o   , drnarea_o, fld_o   , fld_i    )

* -----------------------------------------------------------------
* process each cell on output grid
* -----------------------------------------------------------------

* loop over output grid points

      do jo = 1, nlat_o
         do io = 1, numlon_o(jo)

* for each output cell, sum input overlap weights by drainage type
            
            do nb = 1, ndrn
               wgt(nb) = 0.
            end do
            
            do n = 1, mxovr_i2o                       !overlap input cell index
               ii = iovr_i2o(io,jo,n)                 !lon index  (input grid) of overlap cell
               ji = jovr_i2o(io,jo,n)                 !lat index  (input grid) of overlap cell
               nb = drnbasn_i(ii,ji)                  !basin type (input grid)
               wgt(nb) = wgt(nb) + wovr_i2o(io,jo,n)  !overlap area at io,jo for basin type nb
            end do
            
* sort overalp areas in descending order - keep only 
* top "ndrnmax" elements. drnbasn_o(1,io,jo) is the drainage 
* basin type with largest overlap area, drnbasn_o(2,io,jo) is
* the next largest overlap area, and so on.

            call rankn (ndrn, wgt, ndrnmax, drnbasn_o(1,io,jo))

* normalize these top basin drainage fractions to 1.
* Note that wgt goes from 1 to ndrn, whereas drnfrac_o goes 
* from 1 to ndrnmax for a given [io,jo].

            wt = 0              ! determine normalization weight
            do nb = 1, ndrnmax
               i = drnbasn_o(nb,io,jo)
               wt = wt + wgt(i)
            end do
            if (wt .eq. 0) then
               write(6,*)'(DRAINBAS): total weight of zero encountered'
               write(6,*)'(DRAINBAS): io,jo = ',io,jo
               stop
            endif

            do nb = 1, ndrnmax  ! normalize basin fractions
               i = drnbasn_o(nb,io,jo)
               wgt(i) = wgt(i)/wt 
               drnfrac_o(nb,io,jo) = wgt(i)
            end do

            wt = 0              ! check that total weights sum to 0
            do nb = 1, ndrnmax
               i = drnbasn_o(nb,io,jo)
               wt = wt + drnfrac_o(nb,io,jo)
            end do
            if (wt .lt. 1-eps .or. wt .gt. 1+eps) then
               write(6,*)'(DRAINBAS): normalization incorrect'
               write(6,*)'(DRAINBAS): io,jo = ',io,jo
               stop
            endif

         end do                 ! end of jo loop
      end do                    ! end of io loop

* -----------------------------------------------------------------
* compare global area of each basin type on input and output grids
* -----------------------------------------------------------------

* input grid

      do nb = 1, ndrn
         gast_i(nb) = 0.
      end do

      do ji = 1, nlat_i
         do ii = 1, nlon_i
            nb = drnbasn_i(ii,ji)
            gast_i(nb) = gast_i(nb) + area_i(ii,ji)
         end do
      end do

* output grid

      do nb = 1, ndrn
         gast_o(nb) = 0.
      end do

      do nb = 1, ndrn
         do jo = 1, nlat_o
            do io = 1, numlon_o(jo)
               gast_o(nb) = gast_o(nb) + wgt(nb)*drnarea_o(io,jo)
            end do
         end do
      end do

* area comparison

      write (ndiag,*)
      write (ndiag,'(1x,70a1)') ('=',k=1,70)
      write (ndiag,*) 'Drainage Basin Output'
      write (ndiag,'(1x,70a1)') ('=',k=1,70)

      write (ndiag,*)
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,1001)
1001  format (1x,'drainage basin type',20x,' input grid area',
     &     ' output grid area',/
     &     1x,33x,'     10**6 km**2','      10**6 km**2')
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,*)

      do k = 1, ndrn
         write (ndiag,1002) k,gast_i(k)*1.e-6,gast_o(k)*1.e-6
1002     format (1x,i35,f16.3,f17.3)
      end do

      write (ndiag,*)
      write (ndiag,*) 'For reference the area on the output grid of a',
     & ' cell near the equator is: '
      write (ndiag,'(f10.3,a14)')drnarea_o(1,nlat_o/2)*1.e-06,
     &     ' x 10**6 km**2'
      write (ndiag,*)

      write (6,*) 'successfully made drainage basin types'
      write (6,*)

      return
      end
 
