* ------------------------ code history ---------------------------
* source file:       lsmhis.h
* purpose:           land surface model history and restart file common block
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

      common /lsmhis_c/ locpnh, locpnr, msspnh, msspnr, 
     &                  locfnh, locfnr, ctitle, timcom, rpnt_dir, 
     &                  rpnt_fil,rest_pfile,home_dir,case_nam,
     &                  caseid,logid

      character*80 locpnh          !directory for local history files
      character*80 locpnr          !directory for local restart files
      character*80 msspnh          !mass store path name for history file 
      character*80 msspnr          !mass store path name for restart file
      character*80 locfnh          !local history file name
      character*80 locfnr          !local restart file name
      character*80 ctitle          !case title
      character*80 timcom          !comment: start and end of history interval
      character*80 rpnt_dir        !directory name for local restart pointer file
      character*80 rpnt_fil        !file name for local restart pointer file
      character*80 rest_pfile !full UNIX pathname of restart pointer file
      character*80 home_dir   !full UNIX filepath name of home directory
      character*80 case_nam   !case name
      character*256 caseid         !case title
      character*256 logid          !logname

      common /lsmhis_c/ slfdes, mlfdes 

      character*40 slfdes(mslflds) !description of single-level fields
      character*40 mlfdes(mmlflds) !description of multi-level fields

      common /lsmhis_c/ slfnam, slfuni, slftyp, 
     &                  mlfnam, mlfuni, mlftyp, 
     &                  naver , nmaxi , nmini , 
     &                  ninst , excl  , nswrps, chntyp

      character*8 slfnam(mslflds)  !single-level field: name
      character*8 slfuni(mslflds)  !single-level field: units
      character*8 slftyp(mslflds)  !sin-lev field type: inst, maxi, mini, aver
      character*8 mlfnam(mmlflds)  !multi-level field : name
      character*8 mlfuni(mmlflds)  !multi-level field : units
      character*8 mlftyp(mmlflds)  !mul-lev field type: inst, maxi, mini, aver
      character*8 naver            !average field over history file interval
      character*8 nmaxi            !max field value over history file interval
      character*8 nmini            !min field value over history file interval
      character*8 ninst            !instantaneous field value
      character*8 excl(malflds)    !names of fields to exclude
      character*8 nswrps           !mass store write password for output files
      character*8 chntyp(2,malflds)!paired field names and type for namelist

      common /lsmhis_c/ ninavg

      character*1 ninavg  !equals "q" or "Q" if using monthly averager

      common /lsmhis_c/ sl1dfdes, ml1dfdes 

      character*40 sl1dfdes(mslflds) !description of single-level 1d fields
      character*40 ml1dfdes(mmlflds) !description of multi-level 1d fields

      common /lsmhis_c/ flds1d  ,
     &                  sl1dfnam, sl1dfuni, sl1dftyp,
     &                  ml1dfnam, ml1dfuni, ml1dftyp
 
      character*8  flds1d(malflds)    !names of fields to produce xy output 
      character*10 sl1dfnam(mslflds)  !single-level 1d field: name
      character*8  sl1dfuni(mslflds)  !single-level 1d field: units
      character*8  sl1dftyp(mslflds)  !sin-lev field 1d type: inst, maxi, mini, aver
      character*10 ml1dfnam(mmlflds)  !multi-level 1d field : name
      character*8  ml1dfuni(mmlflds)  !multi-level 1d field : units
      character*8  ml1dftyp(mmlflds)  !mul-lev 1d field type: inst, maxi, mini, aver

      common /lsmhis_i/ nslflds , nmlflds, irt     , nlocfnh,
     &                  nhtfrq  , mfilt  , ntim    , nfil   ,
     &                  mdcur_f , mscur_f, mcdate_f, mcsec_f,
     $                  mdcur   , mscur  , mcdate  , mcsec  ,
     $                  mdcur_t , mscur_t, mcdate_t, mcsec_t,
     $                  nmon    , nyr    , slfcnt  , mlfcnt

      integer nslflds   !number of active single-level fields
      integer nmlflds   !number of active multi-level fields
      integer irt       !mass store retention period
      integer nlocfnh   !unit number for opened local history file
      integer nhtfrq    !history interval (iterations)
      integer mfilt     !max number of time samples per history file
      integer ntim      !number of time samples written to history file
      integer nfil      !number of history files
      integer mdcur_f   !current day when 1st time sample written to hist file
      integer mscur_f   !current seconds of current day for ...
      integer mcdate_f  !current date (yyyymmdd format) for ...
      integer mcsec_f   !current secs or current date for ...
      integer mdcur     !current day for current nstep
      integer mscur     !current seconds of current day for ...
      integer mcdate    !current date (yyyymmdd format) for ...
      integer mcsec     !current seconds of current date for ...
      integer mdcur_t   !current day at start of history interval
      integer mscur_t   !current seconds of current day for ...
      integer mcdate_t  !current date (yyyymmdd format) for ...
      integer mcsec_t   !current secs or current date for ...
      integer nmon      !month (1 to 12)
      integer nyr       !year (00, ...)
      integer slfcnt(mslflds) !number accumulations, sing-lev history field
      integer mlfcnt(mmlflds) !number accumulations, mult-lev history field

      common /lsmhis_i/ nsl1dflds , nml1dflds

      integer nsl1dflds !number of active single-level fields
      integer nml1dflds !number of active multi-level fields

      common /lsmhis_l/  ehi, nlend, ncgetvid, ncopnfil

      logical ehi      !true => current time step is end of history interval
      logical nlend    !true => end of run
      logical ncgetvid !true => need to get netcdf variables id's
      logical ncopnfil !true => netcdf file is open

* ------------------------ end lsmhis.h ---------------------------
 
