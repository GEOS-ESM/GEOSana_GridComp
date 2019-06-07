#include <preproc.h>
      subroutine soiltex (fin     , igeo_i   , nlon_i, nlat_i, igeo_o,
     $                    numlon_o, mxovr_i2o, surtyp, ndiag , sand_o,
     $                    silt_o  , clay_o   )

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       soiltex.F
* purpose:           generate %sand, %silt, %clay for lsm grid
* date last revised: September 1996
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
      integer surtyp(nlon_o,nlat_o) !output grid surface type 
      integer ndiag                 !unit number for diagnostic output  
* -----------------------------------------------------------------

* ------------------------ output variables -----------------------
      real sand_o(nlon_o,nlat_o)    !percent sand on output grid
      real silt_o(nlon_o,nlat_o)    !percent silt on output grid
      real clay_o(nlon_o,nlat_o)    !percent clay on output grid
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer ii                   !longitude index for input grid
      integer ji                   !latitude  index for input grid
      integer k                    !loop index

      real sand_i(nlon_i,nlat_i)   !input grid: percent sand
      real silt_i(nlon_i,nlat_i)   !input grid: percent silt
      real clay_i(nlon_i,nlat_i)   !input grid: percent clay
      integer land_i(nlon_i,nlat_i)!input grid: land mask

      integer iovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lon index of overlap input cell
      integer jovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lat index of overlap input cell
      real    wovr_i2o(nlon_o,nlat_o,mxovr_i2o) !weight    of overlap input cell

      integer io                   !longitude index for lsm grid
      integer jo                   !latitude  index for lsm grid
      integer n                    !loop index

      real w_max  

      integer nlsm                 !number of soil textures
      parameter (nlsm = 4)
      character*35 soil(0:nlsm)    !name of each soil texture
      character*35 typ             !soil texture based on %sand, silt, clay
      real gast_i(0:nlsm)          !input grid : global area, by texture type
      real gast_o(0:nlsm)          !output grid: global area, by texture type

      integer numlon_i(nlat_i)     !input grid number of longitude points at each lat
      real lon_i(nlon_i+1,nlat_i)  !input grid cell longitude, western edge (degrees)
      real lat_i(nlat_i+1)         !input grid cell latitude, south edge (degrees)
      real area_i(nlon_i,nlat_i)   !input grid cell area
      real fld_i (nlon_i,nlat_i)   !input grid dummy field

      real lon_o(nlon_o+1,nlat_o)  !output grid cell longitude, western edge  (degrees)
      real lat_o(nlat_o+1)         !output grid cell latitude, south edge (degree)
      real area_o(nlon_o,nlat_o)   !output grid cell area
      real fld_o(nlon_o,nlat_o)    !output grid dummy field

      integer getavu               !next available unit number
      integer ndata                !unit number of input data
      integer lencs                !length of character string
      character*256 fulpath        !full pathname of input dataset 
      real del_max,del_tst         !delta for maximum test

* -----------------------------------------------------------------

      write (6,*)
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*) 'Making Soil Textures (%sand, %silt, %clay)'
      write (6,'(1x,70a1)') ('=',k=1,70)
      write (6,*)

* -----------------------------------------------------------------
* define the lsm surface types: 0 to nlsm
* -----------------------------------------------------------------

      soil( 0) = 'no soil: ocean, glacier            '
      soil( 1) = 'clays                              '
      soil( 2) = 'sands                              '
      soil( 3) = 'loams                              '
      soil( 4) = 'silts                              '

* -----------------------------------------------------------------
* read in Webb's input data
* -----------------------------------------------------------------

* data are 1 x 1 degree, stored in latitude bands,
* from south to north. in a given latitude band, data begin
* at the date line (180w) and proceed eastward. so first data
* point (x(1,1)) is a box centered at 89.5s, 179.5w.
*
*   89.0S  ---------------------
*          |         |         |
*          |    x    |    x    |
*          |  (1,1)  |  (2,1)  |
*          |         |         |
*   90.0S  ---------------------
*        180.0W    179.0W    178.0W
 
      fulpath = srfpath(1:lencs(srfpath)) //'/'// fin(1:lencs(fin))
      call getfil (fulpath, fin)
      ndata = getavu()
      call opnfil (fin, ndata, 'f')
      do ji = 1, nlat_i
         do ii = 1, nlon_i
            read(ndata,'(3f10.4)')sand_i(ii,ji),silt_i(ii,ji),
     &           clay_i(ii,ji)
            if (sand_i(ii,ji).eq.0. .and. silt_i(ii,ji).eq.0.
     &      .and. clay_i(ii,ji).eq.0.) then
               land_i(ii,ji) = 0
            else
               land_i(ii,ji) = 1
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
* process each cell on output grid: assign output cell the %sand,
* %silt, and %clay of the input cell with the largest overlap.
* can not simply area average %sand, %silt, %clay because this
* will create new soil textures: e.g.

*  .... input cells .............
*  weight  sand silt clay
*   0.50    92   5    3    = sand
*   0.50    22  20   58    = clay

*  .... output cell .............
*  weight  sand silt clay
*   1.00    57  12   31    = sand clay loam

* so whereas 50% of the output cell is in reality sand
* and 50% is clay, the area averaged output cell is a new
* texture, sandy clay loam
* -----------------------------------------------------------------

      do jo = 1, nlat_o
         do io = 1, numlon_o(jo)

* loop through overlapping input cells to find one with greatest 
* area of overlap. note that the following test includes del_max 
* so only 6 significant figures are retained in the comparison.

            w_max = -99999.
            del_max = w_max/1.e6
            k = 0
            do n = 1, mxovr_i2o   
               if ((wovr_i2o(io,jo,n) - del_max) .gt. w_max) then
                  w_max = wovr_i2o(io,jo,n)
                  del_max = w_max/1.e6
                  k = n
               else
                  del_tst = w_max/1.e7
                  if ((wovr_i2o(io,jo,n) - del_tst) .gt. w_max) then
                     write(6,*)'soiltex: WARNING maximum test failed',
     $                    ' at 10^6 and passed at 10^7'
                  endif
               end if
            end do
            if (k .eq. 0) then
               write (6,*) 'soiltex error: did not find input cell ',
     &         'with maximum overlap'
               stop
            end if


* assign cell on output grid the soil of most extensive land input cell

            sand_o(io,jo) = 0.
            silt_o(io,jo) = 0.
            clay_o(io,jo) = 0.

            n = k                      !overlap cell index
            ii = iovr_i2o(io,jo,n)     !lon index (input grid) of overlap cell
            ji = jovr_i2o(io,jo,n)     !lat index (input grid) of overlap cell

            if (land_i(ii,ji) .eq. 1) then
               sand_o(io,jo) = sand_i(ii,ji)
               silt_o(io,jo) = silt_i(ii,ji) 
               clay_o(io,jo) = clay_i(ii,ji) 
            end if

* corrections: 
* o ocean and glaciers have no soil
* o use loam if lsm = land but no soil from input data

            if (surtyp(io,jo).eq.0 .or. surtyp(io,jo).eq.1) then

               sand_o(io,jo) = 0.
               silt_o(io,jo) = 0.
               clay_o(io,jo) = 0. 

            else

               if (sand_o(io,jo).eq.0. .and. silt_o(io,jo).eq.0.
     &         .and. clay_o(io,jo).eq.0.) then
                  sand_o(io,jo) = 43.
                  silt_o(io,jo) = 39.
                  clay_o(io,jo) = 18.
               end if

            end if

* error checks

            if (surtyp(io,jo) .gt. 1) then
               w_max = sand_o(io,jo) + silt_o(io,jo) + clay_o(io,jo)
               if (abs(w_max-100.) .gt. 0.001) then
                  write (6,*) 'soiltex error: sand = ',sand_o(io,jo),
     &            ' silt = ',silt_o(io,jo),' and clay = ',clay_o(io,jo),
     &            ' do not sum to 100 for lon,lat = ',io,jo
                  stop
               end if
            end if

         end do
      end do

* -----------------------------------------------------------------
* compare global area of each soil type on input and output grids
* -----------------------------------------------------------------

* input grid: global areas by texture class

      do k = 0, nlsm
         gast_i(k) = 0.
      end do

      do ji = 1, nlat_i
         do ii = 1, nlon_i
            call soiltyp (sand_i(ii,ji),silt_i(ii,ji),clay_i(ii,ji),typ)
            do k = 0, nlsm
               if (typ .eq. soil(k)) go to 101
            end do
            write (6,*) 'soiltex error: sand = ',sand_i(ii,ji),
     &      ' silt = ',silt_i(ii,ji),' and clay = ',clay_i(ii,ji),
     &      ' not assigned to soil type for input grid lon,lat = ',ii,ji
            stop
101         continue
            gast_i(k) = gast_i(k) + area_i(ii,ji)
         end do
      end do

* output grid: global areas by texture class

      do k = 0, nlsm
         gast_o(k) = 0.
      end do

      do jo = 1, nlat_o
         do io = 1, numlon_o(jo)
            call soiltyp (sand_o(io,jo),silt_o(io,jo),clay_o(io,jo),typ)
            do k = 0, nlsm
               if (typ .eq. soil(k)) go to 102
            end do
            write (6,*) 'soiltex error: sand = ',sand_o(io,jo),
     &     ' silt = ',silt_o(io,jo),' and clay = ',clay_o(io,jo),
     &     ' not assigned to soil type for output grid lon,lat = ',io,jo
            stop
102         continue
            gast_o(k) = gast_o(k) + area_o(io,jo)
         end do
      end do

* comparisons

      write (ndiag,*)
      write (ndiag,'(1x,70a1)') ('=',k=1,70)
      write (ndiag,*) 'Soil Texture Output'
      write (ndiag,'(1x,70a1)') ('=',k=1,70)
      write (ndiag,*)

      write (ndiag,*) 'The following table of soil texture classes is ',
     & 'for comparison only.'
      write (ndiag,*) 'The actual data is continuous %sand, %silt ',
     & 'and %clay not textural classes'
      write (ndiag,*)

      write (ndiag,*)
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,1001)
1001  format (1x,'soil texture class',17x,' input grid area',
     &         ' output grid area',/
     &         1x,33x,'     10**6 km**2','      10**6 km**2')
      write (ndiag,'(1x,70a1)') ('.',k=1,70)
      write (ndiag,*)

      do k = 0, nlsm
         write (ndiag,1002) soil(k),gast_i(k)*1.e-6,gast_o(k)*1.e-6
1002     format (1x,a35,f16.3,f17.3)
      end do

      write (ndiag,*)
      write (ndiag,*) 'For reference the area on the output grid of a',
     & ' cell near the equator is: '
      write (ndiag,'(f10.3,a14)')area_o(1,nlat_o/2)*1.e-06,
     &     ' x 10**6 km**2'
      write (ndiag,*)

      write (6,*) 'successfully made soil textures'
      write (6,*)

      return
      end

      subroutine soiltyp (sand, silt, clay, typ)
      implicit none
      real sand, silt, clay     !percent sand, silt, clay
      character*(*) typ         !soil type

      if (clay.eq.0. .and. sand.eq.0. .and. silt.eq.0.) then
         typ = 'no soil: ocean, glacier'
         return
      end if

      if (clay .ge. 40.) then
         typ = 'clays'
      else if (sand .ge. 50.) then
         typ = 'sands'
      else if (silt .ge. 50.) then
         typ = 'silts'
      else
         typ = 'loams'
      end if
 
      return
      end
 
