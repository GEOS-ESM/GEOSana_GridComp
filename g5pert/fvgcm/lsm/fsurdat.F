#include <preproc.h>
      subroutine fsurdat(igeo_o  , numlon_o, fland   , surtyp, soicol,
     $                   sand    , silt    , clay    , pctlak, pctwet,
     $                   drngbasn, drngarea, drngfrac) 
#include <vartyp.h>
#include <lsmpar.h>
#include <basinpar.h>

* ------------------------ code history ---------------------------
* source file:       fsurdat.F
* purpose:           determines lsm surface dataset
* date last revised: November 1996
* authors:           Gordon Bonan, Mariana Vertenstein
* standardized:
* reviewed:
* -----------------------------------------------------------------
* ------------------------parameters------------------------------------
      integer nlon_o                !maximum number of output grid longitude points
      integer nlat_o                !number of output grid latitude  points
      parameter( nlon_o = lsmlon, nlat_o = lsmlat)
* ----------------------------------------------------------------------

* ------------------------input/output variables------------------------
* input variables
      integer igeo_o                !output grid type in "mn" format
      integer numlon_o(nlat_o)      !output grid number of longitude points for each lat
      real    fland (nlon_o,nlat_o) !output grid fractional land 

* output variables
      integer surtyp(nlon_o,nlat_o) !lsm surface type on output grid
      integer soicol(nlon_o,nlat_o) !lsm soil color on output grid
      real      sand(nlon_o,nlat_o) !percent sand on output grid
      real      silt(nlon_o,nlat_o) !percent silt on output grid
      real      clay(nlon_o,nlat_o) !percent clay on output grid
      real    pctlak(nlon_o,nlat_o) !percent lake on output grid
      real    pctwet(nlon_o,nlat_o) !percent wetland on output grid

      integer drngbasn(ndrnmax,nlon_o,nlat_o)  !output drainage basin matrix
      real    drngfrac(ndrnmax,nlon_o,nlat_o)  !output fracional drainage matrix
      real    drngarea        (nlon_o,nlat_o)  !output grid cell area
 
* ----------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <lsmhis.h>
#include <lsmctl.h>
* ----------------------------------------------------------------------

* ------------------------ local miscellaneous variables ---------------
      integer i,j                 !indices

      real lat_o(nlat_o+1)        !output grid cell latitude , s edge (degrees)
      real lon_o(nlon_o+1,nlat_o) !output grid cell longitude, w edge (degrees)
      real lat_c(nlat_o)          !output grid cell latitude , center (degrees)
      real lon_c(nlon_o  ,nlat_o) !output grid cell longitude, center (degrees)

      character*80 fin_vegt       !file name: vegetation type
      character*80 fin_soic       !file name: soil color
      character*80 fin_soit       !file name: soil texture
      character*80 fin_lake       !file name: lakes
      character*80 fin_swmp       !file name: swamps
      character*80 fin_basin      !file name: basin drainage matrix 
      character*80 msscom         !mass store comment

      integer getavu              !next available unit number
      integer ndiag               !unit number for diagnostic summary output
      integer nsurf               !unit number for diagnostic surface output

      integer mxovr_i2o           !max # input cells that overlap output cell
      integer igeo_i              !input grid type in "mn" format
      integer nlon_i              !input grid max number of longitude points
      integer nlat_i              !input grid number of latitude  points
* ----------------------------------------------------------------------

* -----------------------------------------------------------------
* open diagnostic output log file
* -----------------------------------------------------------------

      ndiag = getavu()
      call opnfil ('lsm_diagnos.logfile', ndiag, 'f')

* -----------------------------------------------------------------
* adjust fractional land values if necessary 
* -----------------------------------------------------------------

#ifdef COUP_CCM      
* values less than cutoff value are set to zero

      do j = 1, nlat_o
         do i = 1, numlon_o(j)
            if (fland(i,j) .lt. cutoff) fland(i,j) = 0.
         end do
      end do

      write (ndiag,*)
      write (ndiag,*) 'using fractional land data '
      write (ndiag,'(1x,a14,f7.4,a29)') 'with cutoff = ',cutoff,
     &     ' to create the output dataset'
#else
* adjustment done in lsmini via land mask sent by flux coupler

      write (ndiag,*)
      write (ndiag,*) 'using fractional land data and land mask ',
     &     'from flux coupler to create output dataset'
#endif

* -----------------------------------------------------------------
* determine output grid southern and western cell edges and 
* latitudes and longitudes at cell centers 
* -----------------------------------------------------------------

      call cell_box (igeo_o, nlat_o, nlon_o, numlon_o, lon_o, lat_o,
     &               lon_c , lat_c)

* -----------------------------------------------------------------
* make land types
* -----------------------------------------------------------------

* determine maximum number of overlapping cells

      nlon_i = 720
      nlat_i = 360
      igeo_i = 11

      call max_ovr (igeo_i   , nlon_i   , nlat_i  , nlon_o ,
     &              nlat_o   , numlon_o , lon_o   , lat_o  , 
     &              mxovr_i2o)

* determine land type

      fin_vegt = 'olson.data'
      call landtyp (fin_vegt, igeo_i   , nlon_i  , nlat_i, igeo_o,
     &              numlon_o, mxovr_i2o, fland   , ndiag , surtyp)

* -----------------------------------------------------------------
* make soil color types
* -----------------------------------------------------------------

* determine maximum number of overlapping cells

      nlon_i = 128
      nlat_i = 64
      igeo_i = 23

      call max_ovr (igeo_i   , nlon_i   , nlat_i  , nlon_o ,
     &              nlat_o   , numlon_o , lon_o   , lat_o  , 
     &              mxovr_i2o)

* determine soil color types

      fin_soic = 'bats.data'
      call soilcol (fin_soic, igeo_i   , nlon_i, nlat_i, igeo_o,
     &              numlon_o, mxovr_i2o, lon_c , lat_c , surtyp,
     &              ndiag   , soicol  )

* -----------------------------------------------------------------
* make soil textures: %sand, %silt, %clay
* -----------------------------------------------------------------

* determine maximum number of overlapping cells

      nlon_i = 360
      nlat_i = 180
      igeo_i = 11

      call max_ovr (igeo_i   , nlon_i   , nlat_i  , nlon_o ,
     &              nlat_o   , numlon_o , lon_o   , lat_o  , 
     &              mxovr_i2o)

* determine soil textures

      fin_soit = 'webb.1x1'
      call soiltex (fin_soit, igeo_i   , nlon_i, nlat_i, igeo_o,
     &              numlon_o, mxovr_i2o, surtyp, ndiag , sand  ,
     &              silt    , clay     )

* -----------------------------------------------------------------
* make inland water: %lake, %wetland
* -----------------------------------------------------------------

* determine maximum number of overlapping cells

      nlon_i = 360
      nlat_i = 180
      igeo_i = 12

      call max_ovr (igeo_i   , nlon_i   , nlat_i  , nlon_o ,
     &              nlat_o   , numlon_o , lon_o   , lat_o  , 
     &              mxovr_i2o)

* determine inland water

      fin_lake = 'flak.1x1'
      fin_swmp = 'swmp.1x1'
      call landwat (fin_lake, fin_swmp, igeo_i   , nlon_i, nlat_i,
     &              igeo_o  , numlon_o, mxovr_i2o, surtyp, ndiag ,
     &              pctlak  , pctwet  ) 

* -----------------------------------------------------------------
* Determine basin drainage matrix
* -----------------------------------------------------------------

      nlon_i = 720
      nlat_i = 360
      igeo_i = 13

* determine maximum number of overlapping cells

      call max_ovr (igeo_i   , nlon_i   , nlat_i  , nlon_o ,
     &              nlat_o   , numlon_o , lon_o   , lat_o  , 
     &              mxovr_i2o)

      
* determine drainage basin matrix on output grid

      fin_basin = 'drainage.p5xp5.nc'
      call basinmat (fin_basin, igeo_i  , nlon_i   , nlat_i,
     $               igeo_o   , numlon_o, mxovr_i2o, ndiag ,
     $               drngbasn , drngfrac, drngarea )

#ifdef COUP_CCM      
* -----------------------------------------------------------------
* extend Antarctica: make Ross ice shelf glacier (***SEE JEFF KEIHL***).
* set all cells south of 70S to glacier. this has the effect of
* making the Ross ice shelf a glacier
* -----------------------------------------------------------------

      if (antartica) then
         write(6,*)'(fsurdat) extending Antartica: making Ross ice',
     &        ' shelf glacier'
         do j = 1, nlat_o
            if (lat_c(j) .lt. -79.) then
               do i = 1, numlon_o(j)
                  fland(i,j)  = 1.
                  surtyp(i,j) = 1
                  soicol(i,j) = 0
                  sand(i,j)   = 0.
                  silt(i,j)   = 0.
                  clay(i,j)   = 0.
                  pctlak(i,j) = 0.
                  pctwet(i,j) = 0.
               end do
            end if
         end do
      endif
#endif

* -----------------------------------------------------------------
* truncate all percentage fields on output grid
* -----------------------------------------------------------------

         do j = 1, nlat_o
            do i = 1, numlon_o(j)
               sand(i,j)   = float(nint(sand(i,j)  ))
               silt(i,j)   = float(nint(silt(i,j)  ))
               clay(i,j)   = float(nint(clay(i,j)  ))
               pctlak(i,j) = float(nint(pctlak(i,j)))
               pctwet(i,j) = float(nint(pctwet(i,j)))
            end do
         end do

* -----------------------------------------------------------------
* write out diagnostic surface data
* -----------------------------------------------------------------

      nsurf = getavu()
      call opnfil ('lsm_diagnos.fsurdat', nsurf, 'f')

      write (nsurf,*)'---------------------------------------'
      write (nsurf,*)' The surface dataset contains the following:'
      write (nsurf,*)' 1st col: latitude  at center of grid cell ',
     &     '(degrees)'
      write (nsurf,*)' 2nd col: longitude at center of grid cell ',
     &     '(degrees)'
      write (nsurf,*)' 3rd col: fractional land at center of grid cell'
      write (nsurf,*)' 4th col: surface type (0 to 28).',
     &     'see vegtypi.F for definitions'
      write (nsurf,*)' 5th col: soil color (1 to 9) ',
     &     'for use with soil albedos'
      write (nsurf,*)' 6th col: soil texture, %sand, ',
     &     'for thermal and hydraulic properties '
      write (nsurf,*)' 7th col: soil texture, %silt,',
     &     'for thermal and hydraulic properties '
      write (nsurf,*)' 8th col: soil texture, %clay, ',
     &     'for thermal and hydraulic properties'
      write (nsurf,*)' 9th col: % of cell covered by lake. ',
     &     'for use as multiple subgrid point'
      write (nsurf,*)'10th col: % of cell covered by wetland. ',
     &     'for use as multiple subgrid point'
      write (nsurf,*)'---------------------------------------'

      do j = 1, nlat_o
         do i = 1, numlon_o(j)
            write (nsurf,10)lat_c(j),lon_c(i,j),fland(i,j),surtyp(i,j),
     &                   soicol(i,j),sand(i,j) ,silt(i,j) ,clay(i,j)  ,
     &                   pctlak(i,j),pctwet(i,j)
   10       format (3f10.4,i3,i3,3f5.0,2f5.0)
         end do
      end do

      call relavu (nsurf) 
      call relavu (ndiag) 

* -----------------------------------------------------------------
* terminate program
* -----------------------------------------------------------------

      write (6,*) '---------------------------------------'

C SJL Changed f5.1 to f5.2
C The following is wrong !!!!
C Meridional resolution should be 180 / (nlat_o - 1)

c     write (6,'(a46,f5.2,a4,f5.2,a5)') 
c    &     'LSM surface data set successfully created for ',
c    &     360./nlon_o,' by ',180./nlat_o,' grid'

      write (6,'(a46,f5.2,a4,f5.2,a5)') 
     &     'LSM surface data set successfully created for ',
     &     360./nlon_o,' by ',180./(nlat_o-1),' grid'
C SJL
      write (6,*)

      msscom = 'diagnostic output: surface data file'
      call putfil ('lsm_diagnos.fsurdat', msspnh, irt, nswrps, msscom)
      msscom = 'diagnostic output: summary of surface data'
      call putfil ('lsm_diagnos.logfile', msspnh, irt, nswrps, msscom)
 
      write (6,*) 'Diagnostic output file = lsm_diagnos.fsurdat. '
      write (6,*) '   This file contains the surface data'
      write (6,*) 'Diagnostic log file = lsm_diagnos.logfile. '
      write (6,*) '   See this file for a summary of the dataset'

      write (6,*) '---------------------------------------'
      write (6,*)

      return
      end
 
