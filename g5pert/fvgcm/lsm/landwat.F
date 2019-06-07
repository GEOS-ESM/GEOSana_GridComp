#include <preproc.h>
      subroutine landwat (fin1  , fin2    , igeo_i   , nlon_i, nlat_i,
     $                    igeo_o, numlon_o, mxovr_i2o, surtyp, ndiag ,
     $                    lake_o, swmp_o)

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       landwat.F
* purpose:           generate %lake and %wetland for lsm grid
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
      character*(*) fin1            !input file name
      character*(*) fin2            !input file name
      integer igeo_i                !input grid type in "mn" format
      integer nlon_i                !input grid number of longitude points
      integer nlat_i                !input grid number of latitude  points
      integer igeo_o                !output grid type in "mn" format (see above)
      integer numlon_o(nlat_o)      !output grid number of longitude points for each latitude 
      integer mxovr_i2o             !max number of input cells that overlap output cell
      integer surtyp(nlon_o,nlat_o) !output grid surface type 
      integer ndiag                 !unit number for diagnostic output  
* -----------------------------------------------------------------

* ------------------------ output variables -----------------------
      real lake_o(nlon_o,nlat_o) !percent lake on output grid
      real swmp_o(nlon_o,nlat_o) !percent wetland on output grid
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer ii                    !longitude index for input grid
      integer ji                    !latitude  index for input grid
      integer k                     !loop index 

      integer lake_i(nlon_i,nlat_i) !input grid: percent lake
      integer swmp_i(nlon_i,nlat_i) !input grid: percent wetland

      integer iovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lon index of overlap input cell
      integer jovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lat index of overlap input cell
      real    wovr_i2o(nlon_o,nlat_o,mxovr_i2o) !weight    of overlap input cell

      integer io                   !longitude index for lsm grid
      integer jo                   !latitude  index for lsm grid
      integer n

      real wt                      !overlap weight
      real w_sum                   !sum of %lake and %wetland

      real glake_o                 !output grid: global area lakes
      real gswmp_o                 !output grid: global area wetlands
      real garea_o                 !output grid: global area
      real glake_i                 !input grid: global area lakes
      real gswmp_i                 !input grid: global area wetlands
      real garea_i                 !input grid: global area

      integer numlon_i(nlat_i)     !input grid number of longitude points at each lat
      real lon_i(nlon_i+1,nlat_i)  !input grid cell longitude, western edge (degrees)
      real lat_i(nlat_i+1)         !input grid cell latitude, south edge (degrees)
      real area_i(nlon_i,nlat_i)   !input grid cell area
      real fld_i (nlon_i,nlat_i)   !input grid dummy field

      real lon_o(nlon_o+1,nlat_o)  !output grid cell longitude, western edge  (degrees)
      real lat_o(nlat_o+1)         !output grid cell latitude, south edge (degree)
      real area_o(nlon_o,nlat_o)   !output grid cell area
      real fld_o (nlon_o,nlat_o)   !output grid dummy field

      integer getavu               !next available unit number
      integer ndata                !unit number of input data
      integer lencs                !length of character string
      character*256 fulpath        !full pathname of input dataset 
* -----------------------------------------------------------------
      write (6,*)
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*) 'Making Inland Water (%lake, %wetland)'
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*)

* -----------------------------------------------------------------
* read in Cogley's input data
* -----------------------------------------------------------------

* data are 1 x 1 degree, stored in latitude bands,
* from south to north. in a given latitude band, data begin
* at greenwich and proceed eastward. so first data
* point (x(1,1)) is a box centered at 89.5s, 0.5e
*
*   89.0S  ---------------------
*          |         |         |
*          |    x    |    x    |
*          |  (1,1)  |  (2,1)  |
*          |         |         |
*   90.0S  ---------------------
*          0.       1.0E      2.0E
 
* percent lake

      fulpath = srfpath(1:lencs(srfpath)) //'/'// fin1(1:lencs(fin1))
      call getfil (fulpath, fin1)
      ndata = getavu()
      call opnfil (fin1, ndata, 'f')
      read (ndata,*) lake_i
      call relavu (ndata) 

* percent wetland
 
      fulpath = srfpath(1:lencs(srfpath)) //'/'// fin2(1:lencs(fin2))
      call getfil (fulpath, fin2)
      ndata = getavu()
      call opnfil (fin2, ndata, 'f')
      read (ndata,*) swmp_i
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

* area average

            lake_o(io,jo) = 0.
            swmp_o(io,jo) = 0.

            do n = 1, mxovr_i2o         !overlap cell index
               ii = iovr_i2o(io,jo,n)  !lon index (input grid) of overlap cell
               ji = jovr_i2o(io,jo,n)  !lat index (input grid) of overlap cell
               wt = wovr_i2o(io,jo,n)  !overlap weight
               lake_o(io,jo) = lake_o(io,jo) + lake_i(ii,ji) * wt
               swmp_o(io,jo) = swmp_o(io,jo) + swmp_i(ii,ji) * wt
            end do

* corrections: set oceans to zero. exclude areas less than 5% of cell

            if (surtyp(io,jo) .eq. 0) then
               lake_o(io,jo) = 0.
               swmp_o(io,jo) = 0.
            else
               if (lake_o(io,jo) .lt. 5.) lake_o(io,jo) = 0.
               if (swmp_o(io,jo) .lt. 5.) swmp_o(io,jo) = 0.
            end if

* error checks

            w_sum = lake_o(io,jo) + swmp_o(io,jo) 
            if (w_sum .gt. 100.) then
               write (6,*) 'landwat error: lake = ',lake_o(io,jo),
     &         ' and wetland = ',swmp_o(io,jo),
     &         ' are greater than 100 for lon,lat = ',io,jo
               stop
            end if

         end do
      end do

* -----------------------------------------------------------------
* compare global areas on input and output grids
* -----------------------------------------------------------------

* input grid

      glake_i = 0.
      gswmp_i = 0.
      garea_i = 0.

      do ji = 1, nlat_i
         do ii = 1, nlon_i
            garea_i = garea_i + area_i(ii,ji)
            glake_i = glake_i + lake_i(ii,ji)*area_i(ii,ji)/100.
            gswmp_i = gswmp_i + swmp_i(ii,ji)*area_i(ii,ji)/100.
         end do
      end do

* output grid

      glake_o = 0.
      gswmp_o = 0.
      garea_o = 0.

      do jo = 1, nlat_o
         do io = 1, numlon_o(jo)
            garea_o = garea_o + area_o(io,jo)
            glake_o = glake_o + lake_o(io,jo)*area_o(io,jo)/100.
            gswmp_o = gswmp_o + swmp_o(io,jo)*area_o(io,jo)/100.
         end do
      end do

* comparison

      write (ndiag,*)
      write (ndiag,'(1x,70a1)') ('=',k=1,70)
      write (ndiag,*) 'Inland Water Output'
      write (ndiag,'(1x,70a1)') ('=',k=1,70)

      write (ndiag,*)
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,2001)
2001  format (1x,'surface type   input grid area  output grid area'/
     &        1x,'                 10**6 km**2      10**6 km**2   ')
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,*)
      write (ndiag,2002) glake_i*1.e-06,glake_o*1.e-06
      write (ndiag,2003) gswmp_i*1.e-06,gswmp_o*1.e-06
      write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002  format (1x,'lakes       ',f14.3,f17.3)
2003  format (1x,'wetlands    ',f14.3,f17.3)
2004  format (1x,'all surface ',f14.3,f17.3)

      write (ndiag,*)
      write (ndiag,*) 'For reference the area on the output grid of a',
     & ' cell near the equator is: '
      write (ndiag,'(f10.3,a14)')area_o(1,nlat_o/2)*1.e-06,
     &     ' x 10**6 km**2'
      write (ndiag,*)

      write (6,*) 'successfully made inland water'
      write (6,*)

      return
      end
 
