* basin drainage - real values

      real    drngfrac(ndrnmax,lsmlon,lsmlat) !fractional drainage basin matrix
      real    drngarea        (lsmlon,lsmlat) !grid cell area

      common /lsmbas_r/ drngfrac, drngarea

* basin drainage - integer values

      integer drngbasn(ndrnmax,lsmlon,lsmlat) !basin type drainage matrix

      common /lsmtbas_i/ drngbasn

* drainage vector 

      real drnvec(nbasmax)         ! total runoff(kg/sec) to each of ndrn basins

      common / lsmbas_r/ drnvec

* drainage vector send to history file

      logical ncbasin              ! true => put runoff on lsm history file

      common / lsmbas_l/ ncbasin

      integer dim_bas_id           ! id for basin vector dimension
      integer var_bas_id           ! id for basin vector variable
      integer bascnt               ! accumulator index

      common / lsmbas_i/ bascnt, dim_bas_id, var_bas_id

      real bashist(ndrn)           ! runoff to each of ndrn basins to history file

      common / lsmbas_r/ bashist

 
