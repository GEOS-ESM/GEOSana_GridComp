#include <preproc.h>
      subroutine max_ovr (igeo_i  , nlon_i  , nlat_i , nlon_o ,
     $                    nlat_o  , numlon_o, lon_o  , lat_o  ,
     $                    novr_max)

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       ao.F
* purpose:           determine maximum number of overlapping cells 
*                    input and output grids
* date last revised: March 1996 
* author:            Gordon Bonan
* standardized: 
* reviewed:    
* -----------------------------------------------------------------

* ------------------- input variables -----------------------------
      integer igeo_i              !grid type in "mn" format
      integer nlon_i              !number of input  longitude points
      integer nlat_i              !number of input  latitude points
      integer nlon_o              !number of output longitude points
      integer nlat_o              !number of output latitude points
      integer numlon_o(nlat_o)    !number of longitude points for each latitude 
      real lon_o(nlon_o+1,nlat_o) !output grid cell longitude, western edge  (degrees)
      real lat_o(nlat_o+1)        !output grid cell latitude , southern edge (degrees)
* -----------------------------------------------------------------

* ------------------- output variables -----------------------------
      integer novr_max            !maximum number of overlapping input cells
* -----------------------------------------------------------------

* ------------------- local variables -----------------------------
      integer novr                !number of overlapping input cells

      integer io,ii               !output and input grids longitude loop index
      integer jo,ji               !output and input grids latitude  loop index

      integer numlon_i(nlat_i)    !input grid cell number of longitude points for each lat
      real lon_i(nlon_i+1,nlat_i) !input grid cell longitude, western edge (degrees)
      real lat_i(nlat_i+1)        !input grid cell latitude, southern edge (degrees)

* -----------------------------------------------------------------

* -----------------------------------------------------------------
* for each output grid cell: find overlapping input grid cell and area of
* input grid cell that overlaps with output grid cell. cells overlap if:
*
* southern edge of input grid < northern edge of output grid AND
* northern edge of input grid > southern edge of output grid
*
* western edge of input grid < eastern edge of output grid AND
* eastern edge of input grid > western edge of output grid
*
*           lon_o(io,jo)       lon_o(io+1,jo)
*
*              |                   |
*              --------------------- lat_o(jo+1)
*              |                   |
*              |                   |
*    xxxxxxxxxxxxxxx lat_i(ji+1)   |
*    x         |   x               |
*    x  input  |   x   output      |
*    x  cell   |   x    cell       |
*    x  ii,ji  |   x   io,jo       |
*    x         |   x               |
*    x         ----x---------------- lat_o(jo  )
*    x             x
*    xxxxxxxxxxxxxxx lat_i(ji  )
*    x             x
* lon_i(ii,ji) lon_i(ii+1,ji)
* -----------------------------------------------------------------

* determine surface cell edges for input grid
        
      do ji=1,nlat_i
         numlon_i(ji) = nlon_i
      end do

      call cell_edge (igeo_i, nlat_i, nlon_i, numlon_i, lon_i, lat_i) 

* determine maximum number of overlapping cells
* loop through all input grid cells to find overlap with output grid.
* code does not vectorize but is only called during initialization.

      novr_max = 0
      do jo = 1, nlat_o
         do io = 1, numlon_o(jo)
            novr = 0
            do ji = 1, nlat_i                            
               if ( lat_i(ji  ).lt.lat_o(jo+1) .and. 
     &              lat_i(ji+1).gt.lat_o(jo  ) ) then                !lat ok
                  do ii = 1, numlon_i(ji)
                     if ( lon_i(ii  ,ji).lt.lon_o(io+1,jo) .and. 
     &                    lon_i(ii+1,ji).gt.lon_o(io  ,jo) ) then    !lon okay
                        novr = novr + 1    ! increment number of ovrlap cells for io,jo
                     end if
                  end do
               end if
            end do
            if (novr .gt. novr_max) then
               novr_max = novr
            endif
         end do
      end do
      
      return
      end
 
