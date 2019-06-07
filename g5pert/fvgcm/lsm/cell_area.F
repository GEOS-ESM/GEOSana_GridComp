#include <preproc.h>
      subroutine cell_area (nlat, nlon, numlon, lon_w, lat_s, re, area)

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       cell_area.F
* purpose:           area of grid cells
* date last revised: March 1996
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------- input variables -----------------------------
      integer nlat            !number of latitude points
      integer nlon            !maximum number of longitude points
      integer numlon(nlat)    !number of longitude points for each latitude 
      real lon_w(nlon+1,nlat) !grid cell longitude, western edge (degrees)
      real lat_s(nlat+1)      !grid cell latitude, southern edge (degrees)
* -----------------------------------------------------------------

* ------------------- output variables ----------------------------
      real re              !radius of earth (km)
      real area(nlon,nlat) !cell area (km**2)
* -----------------------------------------------------------------

* ------------------- local variables -----------------------------
      integer i            !longitude index
      integer j            !latitude index
      real dx              !cell width
      real dy              !cell length
      real deg2rad         !pi/180
* -----------------------------------------------------------------

      deg2rad = (4.*atan(1.)) / 180.
      re = 6371.227709

      do j = 1, nlat
         do i = 1, numlon(j)
            dx = (lon_w(i+1,j)-lon_w(i,j)) * deg2rad
            dy = sin(lat_s(j+1)*deg2rad) - sin(lat_s(j)*deg2rad)
            area(i,j) = dx*dy*re*re
         end do
      end do

      return
      end
 
