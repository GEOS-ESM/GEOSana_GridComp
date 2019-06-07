#include <preproc.h>
      subroutine cell_edge (igeo, nlat, nlon, numlon, lon_w, lat_s)

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       cell_edge.F
* purpose:           southern and western grid cell edges 
*                    of grid cells
* date last revised: November 1996
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------- input variables -----------------------------
      integer igeo         !grid type in "mn" format
      integer nlat         !number of latitude points
      integer nlon         !maximum number of longitude points
      integer numlon(nlat) !number of longitude points for each latitude 
* -----------------------------------------------------------------

* ------------------- output variables ----------------------------
      real lon_w(nlon+1,nlat) !grid cell longitude, western edge  (degrees)
      real lat_s(nlat+1)      !grid cell latitude , southern edge (degrees)
* -----------------------------------------------------------------

* ------------------- local variables -----------------------------
      integer i            !longitude index
      integer j            !latitude index
      integer m            !"m" part of igeo = mn
      integer n            !"n" part of igeo = mn
      real deg2rad         !pi/180
      real dx              !cell width
      real dy              !cell length
      real gauwt(nlat)     !grid cell gaussian weights
      real lat_c(nlat)     !grid cell latitude center (degrees)

* -----------------------------------------------------------------
* get grid type from igeo = mn

*   m = 1: regular grid 
*   m = 2: gaussian grid
*   m = 2: regular grid with pole caps

*   n = 1: grid starts at dateline.  western edge ON dateline
*   n = 2: grid starts at greenwich. western edge ON greenwich
*   n = 3: grid starts at greenwich. is centered  ON greenwich
* -----------------------------------------------------------------

      m = igeo/10 
      if (m.lt.1 .or. m.gt.3) then
         write (6,*)'CELL_EDGE error: lat grid type ',m,' not supported'
         call stoprun
      end if

      n = mod(igeo,10)
      if (n.lt.1 .or. n.gt.3 ) then
         write (6,*)'CELL_EDGE error: lon grid type ',n,' not supported'
         call stoprun
      end if

* -----------------------------------------------------------------
* latitudes: grid starts at south pole (southern edge = -90).
* proceedes northwards to north pole (northern edge = 90)
* -----------------------------------------------------------------

      if (m .eq. 1) then

* regular grid in increments of dy

         dy = 180./nlat
         do j = 1, nlat+1
            lat_s(j) = -90.0 + (j-1)*dy
         end do

      else if (m .eq. 2) then

* gaussian grid with non-uniform latitudes

         deg2rad = (4.*atan(1.)) / 180.               !degrees => radians
         call gaulat(lat_c, gauwt, nlat)              !latitudes at cell center
         do j = 1, nlat                               !degrees, j=1 south pole
            lat_c(j) = -asin(lat_c(j))/deg2rad
         end do
         do j = 1, nlat+1                             !southern edges
            if (j .eq. 1) then                        !south pole
               lat_s(j) = -90.
            else if (j .eq. nlat+1) then              !north pole
               lat_s(j) = 90.
            else                                      !edge = average latitude
               lat_s(j) = (lat_c(j-1) + lat_c(j)) / 2.
            end if
         end do
      else if (m .eq. 3) then
* regular grid in increments of dy

         dy = 180./(nlat-1)
         lat_s(1) = -90.
         lat_s(2) = -90.+dy/2.
         lat_s(nlat) = 90. -dy/2.
         lat_s(nlat+1) = 90.
         do j = 3, nlat-1
            lat_s(j) = lat_s(j-1)+dy
         end do
      end if

* -----------------------------------------------------------------
* longitudes. longitudes for the western edge of the cells must 
* increase continuously and span 360 degrees. three types of grids:

* dateline            :  -180 to 180         (- longitudes west of greenwich)
* greenwich           :     0 to 360
* greenwich (centered): -dx/2 to -dx/2 + 360 (- longitudes west of greenwich)
* -----------------------------------------------------------------

      do j = 1,nlat
         dx = 360./(numlon(j))
         do i = 1, nlon+1
            if (n .eq. 1) lon_w(i,j) = -180.  + (i-1)*dx
            if (n .eq. 2) lon_w(i,j) =    0.  + (i-1)*dx
            if (n .eq. 3) lon_w(i,j) = -dx/2. + (i-1)*dx
         end do
      end do

      return
      end
 
