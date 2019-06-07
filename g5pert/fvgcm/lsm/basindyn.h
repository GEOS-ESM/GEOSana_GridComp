* drainage basin dynamic memory allocation
      
      pointer (pdrnbasn,drnbasn)
      pointer (pdrnfrac,drnfrac)
      pointer (pdrnarea,drnarea)
      pointer (prunoff ,runoff )

      integer drnbasn(ndrnmax,lpt) !basin drainage index array
      real    drnfrac(ndrnmax,lpt) !fractional drainage matrix
      real    drnarea(lpt)         !drainage basin area of grid cell 
      real    runoff (lpt)         !total runoff (kg/sec) at each land point

      common /lsmbas_p/ pdrnbasn, pdrnfrac, pdrnarea, prunoff
 
