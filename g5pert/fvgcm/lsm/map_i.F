#include <preproc.h>
      subroutine map_i (igeo_i  , nlon_i  , nlat_i  , numlon_i ,
     &                  igeo_o  , nlon_o  , nlat_o  , numlon_o ,
     &                  iovr_i2o, jovr_i2o, wovr_i2o, mxovr_i2o,
     $                  lon_i   , lat_i   , area_i  , lon_o    ,
     $                  lat_o   , area_o  , fld_o   , fld_i    )

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       map_i.F
* purpose:           driver for area averaging initialization
* date last revised: November 1996
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* o get indices and weights for area-averaging:

*   from input surface grid to output lsm grid

* o input surface and output lsm grids can be any resolution BUT:

*   both grids must be oriented south to north, i.e., cell(lat+1)
*   must be north of cell(lat). the southern edge of the first row
*   must be -90 (south pole) and the northern edge of the last row
*   must be +90 (north pole)

*   both grids must be oriented eastwards, i.e., cell(lon+1) must be
*   east of cell(lon). but the two grids do not have to start at the
*   same longitude, i.e., one grid can start at dateline and go east;
*   the other grid can start at greenwich and go east. longitudes for
*   the western edge of the cells must increase continuously and span
*   360 degrees. examples
*   dateline            :  -180 to 180         (- longitudes west of greenwich)
*   greenwich           :     0 to 360
*   greenwich (centered): -dx/2 to -dx/2 + 360 (- longitudes west of greenwich)

* o igeo_o and igeo_i, in "mn" format, define the lsm and srf grids:

*   m = 1: regular grid
*   m = 2: gaussian grid

*   n = 1: grid starts at dateline.  western edge ON dateline
*   n = 2: grid starts at greenwich. western edge ON greenwich
*   n = 3: grid starts at greenwich. is centered  ON greenwich

* o field values fld_i on an  input grid with dimensions nlon_i and nlat_i =>
*   field values fld_o on an output grid with dimensions nlon_o and nlat_o as

*   fld_o(io,jo) =
*   fld_i(i_ovr(io,jo,     1 ),j_ovr(io,jo,     1 )) * w_ovr(io,jo,     1 ) +
*                             ... + ... +
*   fld_i(i_ovr(io,jo,mxovr_i),j_ovr(io,jo,mxovr_i)) * w_ovr(io,jo,mxovr_i)

* o error checks:

*   overlap weights of input cells sum to 1 for each output cell

*   global sums of dummy fields are conserved for srf => lsm area-averaging
* -----------------------------------------------------------------

* ------------------- input variables -----------------------------
      integer igeo_o      !output lsm grid type in "mn" format
      integer igeo_i      !input srf grid type in "mn" format

      integer nlon_i      !input srf grid : max number of longitude points
      integer nlat_i      !input srf grid : number of latitude  points
      integer nlon_o      !output lsm grid: max number of longitude points
      integer nlat_o      !output lsm grid: number of latitude  points

      integer numlon_i(nlat_i)  !input grid number of longitude points at each lat
      integer numlon_o(nlat_o)  !input grid number of longitude points at each lat

      integer mxovr_i2o   !maximum number of srf cells that overlap lsm cell
* -----------------------------------------------------------------

* ------------------- output variables ----------------------------
      integer iovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lon index of overlap srf cell
      integer jovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lat index of overlap srf cell
      real    wovr_i2o(nlon_o,nlat_o,mxovr_i2o) !weight    of overlap srf cell
* -----------------------------------------------------------------

* ------------------- local variables -----------------------------
      real lon_o(nlon_o+1,nlat_o) !lsm grid cell longitude, west edge  (degrees)
      real lat_o(nlat_o+1)        !lsm grid cell latitude, south edge (degrees)
      real area_o(nlon_o,nlat_o)  !lsm grid cell area
      real  fld_o(nlon_o,nlat_o)  !dummy lsm grid field

      real lon_i(nlon_i+1,nlat_i) !srf grid cell longitude, west edge (degrees)
      real lat_i(nlat_i+1)        !srf grid cell latitude, south edge (degrees)
      real area_i(nlon_i,nlat_i)  !srf grid cell area
      real  fld_i(nlon_i,nlat_i)  !dummy surface grid field

      real re                     !radius of earth
      real sum_fldo               !global sum of dummy lsm field
      real sum_fldi               !global sum of dummy srf field
      real relerr                 !relative error for error checks

      integer io,ii               !lsm and surface longitude loop indices
      integer jo,ji               !lsm and surface latitude  loop indices
* -----------------------------------------------------------------

      relerr = 0.000001

* -----------------------------------------------------------------
* set up input surface grid southern and western edges. get cell areas
* -----------------------------------------------------------------

      call cell_edge (igeo_i , nlat_i, nlon_i, numlon_i, lon_i, lat_i) 

      call cell_area (nlat_i, nlon_i, numlon_i, lon_i, lat_i, re,
     $                area_i)

* -----------------------------------------------------------------
* set up output lsm grid southern and western edges. get cell areas
* -----------------------------------------------------------------

      call cell_edge (igeo_o, nlat_o, nlon_o, numlon_o, lon_o, lat_o) 

      call cell_area (nlat_o, nlon_o, numlon_o, lon_o, lat_o, re,
     $                area_o)

* -----------------------------------------------------------------
* get indices and weights for mapping from surface grid to lsm grid
* -----------------------------------------------------------------

      call ao_i (nlon_i   , nlat_i   , numlon_i, lon_i    , lat_i ,
     &           nlon_o   , nlat_o   , numlon_o, lon_o    , lat_o ,
     &           mxovr_i2o, iovr_i2o , jovr_i2o, wovr_i2o , re    ,
     &           area_o   , relerr   )

* -----------------------------------------------------------------
* error check: global sum fld_o = global sum fld_i
* -----------------------------------------------------------------

* make dummy surface field and sum globally

      sum_fldi = 0.
      do ji = 1, nlat_i      
         do ii = 1, numlon_i(ji)
            fld_i(ii,ji) = (ji-1)*nlon_i + ii
            sum_fldi = sum_fldi + area_i(ii,ji)*fld_i(ii,ji)
         end do
      end do

* area-average lsm field from surface field

      call area_ave (nlat_i   , nlon_i   , numlon_i ,fld_i    ,
     &               nlat_o   , nlon_o   , numlon_o ,fld_o    ,
     &               iovr_i2o , jovr_i2o , wovr_i2o , mxovr_i2o)

* global sum of lsm field

      sum_fldo = 0.
      do jo = 1, nlat_o
         do io = 1, numlon_o(jo)
            sum_fldo = sum_fldo + area_o(io,jo)*fld_o(io,jo)
         end do
      end do

* check for conservation

      if ( abs(sum_fldo/sum_fldi-1.) .gt. relerr ) then
         write (6,*) 'map_i error srf => lsm: srf field not conserved'
         write (6,'(a23,e20.10)') 'global sum lsm field = ',sum_fldo
         write (6,'(a23,e20.10)') 'global sum srf field = ',sum_fldi
         call stoprun
      end if

      return
      end
 
