* ------------------------ code history ---------------------------
* source file:       lsmctl.h
* purpose:           run control variables 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

      common /lsmctl_i/ nsrest, nestep, nelapse, nbstep,
     &                  mdbase, msbase, mbdate , mbsec ,
     &                  irad  , lsmgeo, ncpday

      integer nsrest  !0: initial run. 1: restart: 3: branch
      integer nestep  !end of run. overrides nelapse (iteration)
      integer nelapse !elapsed time of run (iteration)
      integer nbstep  !start of run (iteration)
      integer mdbase  !base day of run (e.g., 0)
      integer msbase  !base seconds of base day (e.g., 0)
      integer mbdate  !base date of run (yyyymmdd format) (e.g., 00000901)
      integer mbsec   !base seconds of base date (e.g., 0)
      integer irad    !solar radiation frequency (iterations)
      integer lsmgeo  !output grid type in "mn" format where
                      !m = 1: regular grid
                      !m = 2: gaussian grid
                      !n = 1: grid starts at dateline.  western edge ON dateline
                      !n = 2: grid starts at greenwich. western edge ON greenwich
                      !n = 3: grid starts at greenwich. is centered  ON greenwich
      integer ncpday  !number of send/recv calls to flux coupler

      common /lsmctl_c/ flondat, flnddat, finidat, srfpath, nrevsn  

      character*80 finidat  !initial conditions file name
      character*80 flondat  !number of longitudes per latitude band filename 
      character*80 flnddat  !fractional land data file name
      character*80 srfpath  !local or mss path for initial surface datasets
      character*80 nrevsn   !fulll restart pathname for branch run
 
      common /lsmctl_r/ dtlsm, dtsoi, cutoff

      real dtlsm      !main lsm time step = atmospheric model time step (s)
      real dtsoi      !soil hydrology time step (s), <= lsm time step
      real cutoff     !minimum fractional land for land

      common /lsmctl_l/ hydro, pergro, conchk, antartica, flxave

      logical hydro      !true if using prognostic hydrology
      logical pergro     !true if random perturbation growth test
      logical conchk     !true if want error energy and water conservation checks
      logical antartica  !if true, extend Antartica for Ross Ice Shelf->glacier
      logical flxave     !if true, only communicate with the flux coupler on
                         !time steps where albedos are calculated

* ------------------------ end lsmctl.h ----------------------------


 
