* ------------------------ code history ---------------------------
* source file:       lsmtcdyn.h
* purpose:           lsmtcdyn common block for time constant 
*                    variables with dynamic memory allocation 
*                    in lsmini.F 
* date last revised: August 1996 (M. Vertenstein)
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* -----------------------------------------------------------------

* lpt and kpt will be determined at run time
* set up pointers so that can dynamically allocate memory for
* arrays dependent on these lengths

      pointer (pixy   ,ixy   )           
      pointer (pjxy   ,jxy   )          
      pointer (pkvec  ,kvec  )
      pointer (pklnd  ,klnd  )          
      pointer (pwsg2g ,wsg2g )
      pointer (pivt   ,ivt   ) 
      pointer (pist   ,ist   ) 
      pointer (pisc   ,isc   ) 
      pointer (pwatsat,watsat) 
      pointer (phksat ,hksat )  
      pointer (psmpsat,smpsat) 
      pointer (pbch   ,bch   )    
      pointer (pwatdry,watdry) 
      pointer (pwatopt,watopt) 
      pointer (pcsol  ,csol  )   
      pointer (ptksol ,tksol )  
      pointer (ptkdry ,tkdry )  
      pointer (pdzsoi ,dzsoi )
      pointer (pzsoi  ,zsoi  )
      pointer (proot  ,root  )
      pointer (psand  ,sand  )   
      pointer (pclay  ,clay  )   
      pointer (plati  ,lati  )   
      pointer (plong  ,long  )
      pointer (pbegkpt,begkpt)
      pointer (pnumkpt,numkpt)

* the land surface model works by gathering all the land points on a
* [lsmlon] x [lsmlat] grid into a vector of [lpt] land points. this
* is then expanded into a vector of [kpt] subgrid points, allowing
* for up to [msub] subgrid points per land point. [ixy], [jxy], [kvec],
* and [klnd] are indices for the mapping: [lsmlon] x [lsmlat] grid <->
* [lpt] vector of land points <-> [kpt] vector of subgrid points. [wsg2g]
* are the weights to obtain the land average from the subgrid points.
* [numlon] <= [lsmlon] allows for for variable longitudinal resolution
* for each latitude strip

      integer ixy(lpt)        !longitude index of lpt land points: 1 to lsmlon
      integer jxy(lpt)        !latitude index of lpt land points: 1 to lsmlat
      integer kvec(lpt,msub)  !msub subgrid vector indices for lpt land points: 1 to kpt
      real   wsg2g(lpt,msub)  !msub subgrid weights for lpt land points
      integer klnd(kpt)       !land index for kpt subgrid points: 1 to lpt 

      common /lsmtc_p/ pixy, pjxy, pkvec, pklnd, pwsg2g

* the "big" vectors of [kpt] subgrid points are processed as [numlv]-1
* "little" vectors of [lvec] points and one "little" vector of [mpt]
* points. for each "little" vector: [begkpt] (<= [kpt]) is the first
* point in the "big" [kpt] vector to process. [numkpt] is the number of
* points in the "big" [kpt] vector to process

      integer begkpt(numlv)  !first point in "big" vector to process
      integer numkpt(numlv)  !number of points in "big" kpt vector to process

      common /lsmtc_p/ pbegkpt, pnumkpt

* time-invariant boundary data for the [kpt] subgrid points

      integer ivt(kpt)     !vegetation type: see vegconi.F for definitions 
      integer ist(kpt)     !"soil" type: see soiconi.F for definitions 
      integer isc(kpt)     !color classes for soil albedos
      real watsat(kpt)     !volumetric soil water content at saturation (porosity)
      real hksat(kpt)      !hydraulic conductivity at saturation (mm h2o /s) 
      real smpsat(kpt)     !soil matrix potential at saturation (mm) 
      real bch(kpt)        !clapp and hornberger "b"
      real watdry(kpt)     !water content when evapotranspiration stops
      real watopt(kpt)     !optimal water content for evapotranspiration 
      real csol(kpt)       !specific heat capacity, soil solids (j/m**3/kelvin)
      real tksol(kpt)      !thermal conductivity, soil solids (w/m/kelvin)
      real tkdry(kpt)      !thermal conductivity, dry soil (w/m/kelvin)
      real dzsoi(msl,kpt)  !soil node thickness (m)
      real zsoi(msl,kpt)   !soil node depth (m)
      real root(msl,kpt)   !root fraction
      real sand(kpt)       !percent sand
      real clay(kpt)       !percent clay

      common /lsmtc_p/ pivt    ,pist   ,pisc    ,pwatsat ,phksat , 
     &                 psmpsat ,pbch   ,pwatdry ,pwatopt ,pcsol  , 
     &                 ptksol  ,ptkdry ,pdzsoi  ,pzsoi   ,proot  ,
     &                 psand   ,pclay

* latitudes and longitudes

      real lati(kpt)             !latitude  for subgrid vector
      real long(kpt)             !longitude for subgrid vector

      common /lsmtc_p/ plati, plong

* ------------------------ end lsmtcdyn.h -------------------------


 
