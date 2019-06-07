#include <preproc.h>
      block data vegtypi

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       vegtypi.F
* purpose:           subgrid plant type and fractional area for surface types
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

#include <vegtyp.h>

      integer i          ! loop index

* there are 29 land surface types: 0 = ocean, 1 to 28 = land. each
* land point has up to three vegetation types, ranging in value from
* 1 to 14. [plant] contains the vegetation type of the 3 subgrid points 
* for each surface type. [cover] contains the fractional area of the 3 
* subgrid points for each surface type.
 
      data (plant(i,1),i=1,29) /   0,
     &                            14,  14,   1,   2,   4,   1  , 1, 
     &                             4,   1,   3,   5,  13,   1,   2,
     &                            11,  11,   6,  13,   9,   7,   8,
     &                             8,  12,  11,  12,  11,   3,  14/
      data (cover(i,1),i=1,29) /0.00,
     &                          1.00,1.00,0.75,0.50,0.75,0.37,0.75,
     &                          0.75,0.37,0.95,0.75,0.70,0.25,0.25,
     &                          0.40,0.40,0.60,0.60,0.30,0.80,0.80,
     &                          0.10,0.85,0.85,0.85,0.85,0.80,1.00/
      data (plant(i,2),i=1,29) /   0,
     &                            14,  14,  14,  14,  14,   4  ,14,
     &                            14,   4,  14,  14,   5,  10,  10,
     &                             4,   4,  13,   6,  10,  14,  14,
     &                            14,  14,  14,  14,  14,  14,  14/
      data (cover(i,2),i=1,29) /0.00,
     &                          0.00,0.00,0.25,0.50,0.25,0.37,0.25,
     &                          0.25,0.37,0.05,0.25,0.30,0.25,0.25,
     &                          0.30,0.30,0.20,0.20,0.30,0.20,0.20,
     &                          0.90,0.15,0.15,0.15,0.15,0.20,0.00/
      data (plant(i,3),i=1,29) /   0,
     &                            14,  14,  14,  14,  14,  14,  14,
     &                            14,  14,  14,  14,  14,  14,  14,
     &                             1,   1,  14,  14,  14,  14,  14,
     &                            14,  14,  14,  14,  14,  14,  14/
      data (cover(i,3),i=1,29) /0.00,
     &                          0.00,0.00,0.00,0.00,0.00,0.26,0.00,
     &                          0.00,0.26,0.00,0.00,0.00,0.50,0.50,
     &                          0.30,0.30,0.20,0.20,0.40,0.00,0.00,
     &                          0.00,0.00,0.00,0.00,0.00,0.00,0.00/

* ------------------------------------------------------------------
* description of the 29 surface types 
* ------------------------------------------------------------------

* no vegetation
* -------------
*  0 ocean                                 
*  1 land ice (glacier)                             
*  2 desert                                

* forest vegetation
* -----------------
*  3 cool needleleaf evergreen tree           
*  4 cool needleleaf deciduous tree           
*  5 cool broadleaf  deciduous tree           
*  6 cool mixed needleleaf evergreen and broadleaf deciduous tree    
*  7 warm needleleaf evergreen tree           
*  8 warm broadleaf  deciduous tree            
*  9 warm mixed needleleaf evergreen and broadleaf deciduous tree    
* 10 tropical broadleaf evergreen tree  
* 11 tropical seasonal deciduous tree         

* interrupted woods
* ----------------
* 12 savanna                               
* 13 evergreen forest tundra               
* 14 deciduous forest tundra               
* 15 cool forest crop                           
* 16 warm forest crop                           

* non-woods
* ---------
* 17 cool grassland                             
* 18 warm grassland                            
* 19 tundra                              
* 20 evergreen shrub                   
* 21 deciduous shrub                 
* 22 semi-desert                     
* 23 cool irrigated crop                
* 24 cool non-irrigated crop                
* 25 warm irrigated crop               
* 26 warm non-irrigated crop               

* wetlands
* --------
* 27 forest (mangrove)                    
* 28 non-forest                          
* ------------------------------------------------------------------

* ------------------------------------------------------------------
* description of the 14 plant types. see vegconi.F for parameters
* that depend on vegetation type
* ------------------------------------------------------------------

*  1 = needleleaf evergreen tree
*  2 = needleleaf deciduous tree
*  3 = broadleaf evergreen tree
*  4 = broadleaf deciduous tree
*  5 = tropical seasonal tree
*  6 = cool grass (c3)
*  7 = evergreen shrub
*  8 = deciduous shrub
*  9 = arctic deciduous shrub
* 10 = arctic grass
* 11 = crop
* 12 = irrigated crop
* 13 = warm grass (c4)
* 14 = not vegetated

      end
 
