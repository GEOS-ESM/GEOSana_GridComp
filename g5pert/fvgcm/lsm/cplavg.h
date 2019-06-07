* ------------------------ code history ---------------------------
* source file:       lsmavg.h
* purpose:           maintain time averages for flux coupled run
* date last revised: Nov 1996 - lsm version 1
* author:            Mariana Vertenstein
* standardized:
* reviewed:
* -----------------------------------------------------------------

      logical dosend       !true if data is to be sent on current time step
      logical dorecv       !true if data is to be recd on current time step
      common /cplavg_l/ dosend, dorecv

      integer ncnt         !number of steps over which to average output fluxes 
      integer icnt         !step counter for flux averager 

      common /cplavg_i/ icnt, ncnt


      real tauxa (lsmlon,lsmlat)   ! accumulator for land flux
      real tauya (lsmlon,lsmlat)   ! accumulator for land flux
      real lhflxa(lsmlon,lsmlat)   ! accumulator for land flux
      real shflxa(lsmlon,lsmlat)   ! accumulator for land flux
      real lwupa (lsmlon,lsmlat)   ! accumulator for land flux
      real qflxa (lsmlon,lsmlat)   ! accumulator for land flux
      real drnveca(nbasmax)        ! total runoff(kg/sec) to each of ndrn basins

      common /cplavg_r/ tauxa, tauya, lhflxa, shflxa, lwupa, qflxa,
     $                  drnveca


 
