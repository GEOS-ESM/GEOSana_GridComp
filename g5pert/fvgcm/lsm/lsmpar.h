* ------------------------ code history ---------------------------
* source file:       lsmpar.h
* purpose:           land surface model array dimensions
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* -----------------------------------------------------------------

* define land surface 2-d grid 

      integer lsmlon  !maximum number of longitude points on lsm grid
      integer lsmlat  !number of latitude points on lsm grid
      integer msub    !maximum number of subgrid points
 
      parameter (lsmlon=LSMLON, lsmlat=LSMLAT, msub=5)

* define miscellaneous lsm parameters

      integer mband   !number of solar radiation bands: vis, nir
      integer msl     !number of soil layers
      integer mst     !number of "soil" types (soil, ice, 2 lakes, wetland)
      integer mvt     !number of plant types
      integer msc     !number of soil color types

      parameter (mband=2, msl=6, mst=5, mvt=14, msc=9)

* define number of single-level equivalent variables for restart

      integer mlsf    !number of land surface fields to save for each point
      parameter (mlsf = 21 + 2*msl + 9*mband)

* define maximum number of history fields

      integer mslflds   !maximum number of active single-level fields
      integer mmlflds   !maximum number of active multi-level fields
      integer malflds   !maximum number of active fields (all levels)

      parameter (mslflds=60, mmlflds=3, malflds=mslflds+mmlflds)

* ------------------------ end lsmpar.h ---------------------------
 
