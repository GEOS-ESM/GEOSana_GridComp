#include <preproc.h>
      subroutine lsmctli (xmdbase, xmsbase, xmbdate, xmbsec,
     &                    xnsrest, xnstep , xnestep, dtatm ,
     &                    xnhtfrq, xmfilt , xninavg, xirt  ,
     &                    xirad  )

#if defined ( SPMD ) 
      use mod_comm, only: gid
#endif

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>

* ------------------------ code history --------------------------------
* source file:       lsmctli.F
* purpose:           initialize run control variables in lsmctl and 
*                    lsmhis common blocks
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* ----------------------------------------------------------------------

* ------------------------ notes ---------------------------------------
* when coupled to ccm: 
* base calendar info, nstep, nestep, nsrest, time step, mfilt, nhtfrq, 
* ninavg are input to lsm from ccm. the values in the lsmexp 
* namelist are not used. the minimum lsmexp namelist parameters are:
*    o finidat
*
* when uncoupled to ccm:
* the fractional land mask is provided by the flux coupler so flnddat 
* is not used. note that flondat is not currently used.  
* it will be used when the lsm reduced grid is implemented. 
* the minimum lsmexp namelist parameters are:
*    o flondat
*    o finidat
*    o nsrest
*    o nestep or nelapse (set to large negative integer like -9999) 
*    o nnbdat, nnbsec, nndbas, nnsbas
*    o dtime
*
* note: when coupled to the ccm, the user may specify that the fractional 
* land data may be read in from a specified dataset, rather than be 
* obtained as input to lsmini from the atm. If this is done, the user must 
* be sure that the fractional land data is consistent with other atm 
* surface fields (such as the geopotential height, etc.)
* ----------------------------------------------------------------------

* ------------------------ input variables -----------------------------
      integer xmdbase      !base day of run (e.g., 0)
      integer xmsbase      !base seconds of base day (e.g., 0)
      integer xmbdate      !base date of run (yymmdd format) (e.g., 000901)
      integer xmbsec       !base seconds of base date (e.g., 0)
      integer xnstep       !current time index
      integer xnestep      !ending time index
      integer xnsrest      !0 = initial run, > 0 = continuation run
      real dtatm           !time step (s)
      integer xnhtfrq      !history interval(iterations)
      integer xmfilt       !max number of time samples per history file
      integer xirt         !mss retention time
      integer xirad        !solar radiation frequency (iterations)
      character*1 xninavg  !equals "q" or "Q" if using monthly averager
* ----------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <dynmem.h>
#include <lsmhis.h>
#include <lsmctl.h>
* ----------------------------------------------------------------------

* ------------------------ local variables -----------------------------
      character*7 runtyp(4)   !run type
      integer i               !loop index
      integer iundef          !integer undefined 
      real rundef             !real undefined
      integer lencs           !integer function to get length of string
      integer nndbas          !base day of run
      integer nnsbas          !base seconds of base day
      integer nnbdat          !base date of run (yymmdd format)
      integer nnbsec          !base seconds of base date
      real    dtime           !time step
      character*256 text      !text for shell commands
      logical mkdir           !if true, make directory
      integer len             !index for length of string
      integer ierr            !error value

* ----------------------------------------------------------------------

* ------------------------ namelist variables --------------------------
      namelist /lsmexp/ ctitle    , nsrest   , nestep   , nelapse  ,
     &                  flnddat   , flondat  , finidat  , nrevsn   ,
     $                  nswrps    , irt      , nhtfrq   , ninavg   ,
     $                  mfilt     , flds1d   , excl     , lsmgeo   ,
     $                  locpnh    , locpnr   , srfpath  , 
     &                  nnbdat    , nnbsec   , nndbas   , nnsbas   , 
     &                  dtime     , hydro    , pergro   , conchk   , 
     &                  irad      , chntyp   , cutoff   , antartica,
     &                  flxave    , rest_pfile

* === define run ===
*    o ctitle  = 80 character case title
*    o nsrest  = integer flag. 0: initial run. 1: restart: 3: branch
*    o nestep  = integer end of run. supercedes nelapse (+ = iteration. - = day)
*    o nelapse = integer elapsed time of run (+ = iteration. - = day)

* === generation of surface data === 
*    o cutoff    = real minimum fractional land for land
*    o flondat   = 80 character input file name: number of lon points per latitude
*                  (only to be used for flux coupled runs and not currently used)
*    o flnddat   = 80 character input file name: fractional land 
*    o antartica = if true, extend Antartica for Ross Ice Shelf->glacier 
*    o lsmgeo    = integer, output grid type in "mn" format where
*                  m = 1: regular grid
*                  m = 2: gaussian grid
*                  n = 1: grid starts at dateline.  western edge ON dateline
*                  n = 2: grid starts at greenwich. western edge ON greenwich
*                  n = 3: grid starts at greenwich. is centered  ON greenwich
*    o srfpath   = 80 character local or mss path where input surface datasets reside 

* === history and restart files ===
*    o finidat    = 80 character initial conditions file name for use with initial run
*    o nhtfrq     = integer history interval (+ = iterations. - = hours)
*    o ninavg     = 'q' or 'Q' if using monthly average history option
*    o mfilt      = integer number of time samples per history file
*    o flds1d     = 8 character names of fields to output as 1d (up to malflds names)
*    o excl       = 8 character names of fields to exclude (up to malflds names)
*    o chntyp     = paired 8 character field name and field type. overrides
*                   default settings in subroutine histlst. e.g., 'TV','maximum'
*    o locpnh     = 80 character directory for local history files
*    o locpnr     = 80 character directory for local restart files
*    o rest_pfile = 80 character full UNIX pathname of the local restart pointer file.
*                   This file must exist when the model is restarted. This file is 
*                   overwritten and updated every time new restart data files are output. 
*    o nrevsn     = 80 character restart file name for use with branch run
*    o nswrps     = 8 character mass store write password for output data sets
*    o irt        = integer mass store retention period (days)

* === model time ===
*    o nnbdat  = integer base date of run (yymmdd format) (e.g., 000901)
*    o nnbsec  = integer base seconds of base date (e.g., 0)
*    o nndbas  = integer base day of run (e.g., 0)
*    o nnsbas  = integer base seconds of base day (e.g., 0)
*    o dtime   = real model time step (s)

* === model physics ===
*    o hydro   = true if using prognostic hydrology
*    o pergro  = true if random perturbation growth test
*    o conchk  = true if want error energy and water conservation checks
*    o irad    = integer solar radiation frequency (+ = iterations. - = hours)

* === flux averaging ===
*    o flxave  = true if flux averaging is to be performed, only used when
*                uncoupled to ccm 
* ----------------------------------------------------------------------

#if ( !defined SPMD ) 
      integer gid
      gid = 0
#endif

      runtyp(0 + 1) = 'initial'
      runtyp(1 + 1) = 'restart'
      runtyp(3 + 1) = 'branch '
      iundef = -9999999
      rundef = -9999999.

* ----------------------------------------------------------------------
* default values
* ----------------------------------------------------------------------

* default values for namelist variables

      ctitle  = ' '
      nsrest  = iundef
      nestep  = iundef
      nelapse = iundef

      nrevsn  = ' '

      lsmgeo    =  33
#ifdef COUP_CCM
      antartica = .true.
#else
      antartica = .false.
#endif
      cutoff    = .50

      finidat   = 'arbitrary initialization'
#ifdef COUP_CCM
      flnddat   = 'fractional land data obtained from ccm model'
#else
      flnddat   = 'fractional land data obtained from flux coupler'
      flondat   = ' '
#endif

      locpnh = '.'
      locpnr = '.'
C SJL
c     srfpath = 'not_specified'
C Modified by SJL to use the current working directory.
      srfpath = '.'
C SJL

      nswrps = '        '
      irt = 0

      nhtfrq = -24
      ninavg = ' '
      mfilt = 1
      do i = 1, malflds
         flds1d(i)   = ' '
         excl(i)     = ' '
         chntyp(1,i) = ' '
         chntyp(2,i) = ' '
      end do

      do i = 1,80
         home_dir(i:i)   = ' '
         case_nam(i:i)   = ' '
         rest_pfile(i:i) = ' '
      end do
      call getenv ('HOME',home_dir)
      call getenv ('CASE',case_nam) 
      nnbdat  = iundef
      nnbsec  = iundef
      nndbas  = iundef
      nnsbas  = iundef
      dtime   = rundef

      hydro   = .true.
      pergro  = .false.
      conchk  = .false.
      irad    = -1

      flxave  = .true.

* default values for other variables

      nbstep  = 0

* ----------------------------------------------------------------------
* read namelist from standard input. override if coupled to ccm
* ----------------------------------------------------------------------

      if (masterproc) then
         read (12,lsmexp)
      endif	

C SJL --   from namelist
c     rest_pfile = home_dir(1:lencs(home_dir))
c    &     //'/lnd.'// case_nam(1:lencs(case_nam)) //'.rpointer'

#if ( defined SPMD )
! SJL note: lsmdistnl broadcasts nemalist vars to all processes
!     call lsmdistnl(nnbdat  , nnbsec   ,nndbas  , nnsbas ,dtime )
#endif

      nsrest  = xnsrest
      nestep  = xnestep
      nnbdat  = xmbdate    
      nnbsec  = xmbsec      
      nndbas  = xmdbase  
      nnsbas  = xmsbase   
      dtime   = dtatm         
      nbstep  = xnstep
      nhtfrq  = xnhtfrq
      mfilt   = xmfilt
      ninavg  = xninavg
      irt     = xirt
      irad    = xirad

* ----------------------------------------------------------------------
* main model time step (dtlsm) and time step for soil hydrology (dtsoi)
* if time step <= 10 min, do soil hydrology with     dtlsm time step
* if time step  = 20 min, do soil hydrology with 1/2 dtlsm time step (10 min)
* if time step  = 30 min, do soil hydrology with 1/3 dtlsm time step (10 min)
* ----------------------------------------------------------------------

      if (dtime  .eq. rundef) then
         if (masterproc) write (6,*) 'lsmctli error: must set dtime'
         call stoprun
      end if

      dtlsm = dtime

      if (dtlsm <= 600.) then
         dtsoi = dtlsm
      else if (dtlsm >  600. .and. dtlsm <= 1200.) then
         dtsoi = dtlsm/2.
      else if (dtlsm > 1200. .and. dtlsm <= 1800.) then
         dtsoi = dtlsm/3.
      else if (dtlsm > 1800. .and. dtlsm <= 2400.) then
         dtsoi = dtlsm/4.
      else
         dtsoi = dtlsm/6.
      end if

      if (masterproc) then
         write (6,*) 'model time steps:'
         write (6,*) '   main time step (s)           = ',dtlsm
         write (6,*) '   soil hydrology time step (s) = ',dtsoi
      end if	

* model calendar

      if (nnbdat .eq. iundef) then
         if (masterproc) write (6,*) 'lsmctli error: must set nnbdat'
         call stoprun
      end if
      if (nnbsec .eq. iundef) then
         if (masterproc) write (6,*) 'lsmctli error: must set nnbsec'
         call stoprun
      end if
      if (nndbas .eq. iundef) then
         if (masterproc) write (6,*) 'lsmctli error: must set nndbas'
         call stoprun
      end if
      if (nnsbas .eq. iundef) then
         if (masterproc) write (6,*) 'lsmctli error: must set nnsbas'
         call stoprun
      end if

      mdbase = nndbas           !base day of run
      msbase = nnsbas           !base seconds of base day
      mbdate = nnbdat           !base date of run (yymmdd format)
      mbsec  = nnbsec           !base seconds of base date

      if (masterproc) then
      write (6,*) 'model calendar:'
      write (6,*) '   base date of run (yymmdd format) = ',mbdate
      write (6,*) '   base seconds of base date        = ',mbsec
      write (6,*) '   base day of run                  = ',mdbase
      write (6,*) '   base seconds of base day         = ',msbase

* define run 

      write (6,*) 'define run:'
      write (6,*) '   case title            = ',ctitle(1:lencs(ctitle))
      end if

      if (nsrest .eq. iundef) then
         if (masterproc) write (6,*) 'lsmctli error: must set nsrest'
         call stoprun
      end if
      if (masterproc) 
     $     write (6,*) '   run type              = ',runtyp(nsrest+1)

      if (nestep .ne. iundef) then
         nelapse = 0
         if (nestep  .lt. 0) nestep  = nint(-nestep *86400./dtime)
         if (masterproc) 
     $        write (6,*) '   last time step (iter) = ',nestep
      else if (nelapse .ne. iundef) then
         nestep = 0
         if (nelapse .lt. 0) nelapse = nint(-nelapse*86400./dtime)
         if (masterproc) 
     $        write (6,*) '   elapsed time (iter)   = ',nelapse
      else
         if (masterproc) 
     $        write (6,*) 'lsmctli error: must set nestep or nelapse'
         call stoprun
      end if
 
* for initial run and for restart, override nrevsn, 
* for print out purposes only, so that the print out is
* consistent with what the model does. 

      if (nsrest .eq. 0) nrevsn = ' '
      if (nsrest .eq. 1) nrevsn = 'set by restart pointer file'
      if (nsrest .eq. 3 .and. nrevsn .eq. ' ') then
         if (masterproc) 
     $        write (6,*) 'lsmctli error: need to set restart',
     $        ' data file name' 
         call stoprun
      end if

      if (masterproc) then 
         write (6,*)'input data files:'
         write (6,*)'  initial  data         = ',
     $        finidat(1:lencs(finidat))
         if (nsrest .ne. 0) then
            write (6,*)'  restart  data         = ',
     $           nrevsn(1:lencs(nrevsn))
         endif
      end if

* initial surface data

      if (nsrest.eq.0 .and. finidat.eq.'arbitrary initialization') then
         if (masterproc) then 
            write (6,*) '  fractional land data  = ',
     &           flnddat(1:lencs(flnddat))
            if (srfpath .eq. 'not_specified') then
               write (6,*) 'lsmctli error: need to specify srfpath'
               call stoprun
            else
               write (6,*) '  local or mss path for input surface ',
     $              'data = ',srfpath(1:lencs(srfpath))
            endif
         end if
      endif

* history and restart files 

      if (nhtfrq .lt. 0) nhtfrq = nint(-nhtfrq*3600./dtime)
      irt = min0(irt,1825)
         if (masterproc)
     $     write (6,*) 'history and restart parameters:'

      if (ninavg.eq.'q' .or. ninavg.eq.'Q') then
         nhtfrq = 0
         mfilt = 1
         if (masterproc) 
     $      write(6,*)'   history time interval              = monthly'
      else
         if (masterproc)
     $      write(6,*)'   history time interval (iterations) = ',nhtfrq
      end if
      if (masterproc) then 
         write(6,*)'   time samples per history file      = ',mfilt
         write (6,*) '   directory for local history files  = ',
     &        locpnh(1:lencs(locpnh))
         write (6,*) '   directory for local restart files  = ',
     &        locpnr(1:lencs(locpnr))
         write (6,*) '   mass store retention (days) = ',irt
         write (6,*) '   mass store write password = ',nswrps
      end if
* restart pointer file
*  o split the full pathname of the restart pointer file into a 
*    directory name and a file name
*  o check if the directory exists and if not, make it 

      do i=1,80
         rpnt_dir(i:i) = ' '
         rpnt_fil(i:i) = ' '
      end do

      do i=1,256
         text(i:i)     = ' '
      end do

      do len = lencs(rest_pfile),0,-1     
         if (rest_pfile(len:len).eq. '/') then
            rpnt_dir = rest_pfile(1:len-1)
            rpnt_fil = rest_pfile(len+1:lencs(rest_pfile))
            go to 100
         endif
      enddo
      rpnt_dir = "."                        ! no "/" found, set path = "."
      rpnt_fil = rest_pfile                 ! no "/" found, use whole input string.
  100 continue

      if ( masterproc ) then
         mkdir = .false.
         inquire (file=rpnt_dir, exist=mkdir)
         if (.not. mkdir) then  ! if directory does not exist
            text = "mkdir -p " // rpnt_dir(1:lencs(rpnt_dir))
            call   system(text(1:lencs(text))) ! Non-cray shell call
         endif

         if(masterproc) write (6,*) '   directory for restart pointer files  = ',
     &        rpnt_dir(1:lencs(rpnt_dir))
         if(masterproc) write (6,*) '   file name for restart pointer files  = ',
     &        rpnt_fil(1:lencs(rpnt_fil))
      end if
* model physics

      if (irad .lt. 0) irad = nint(-irad*3600./dtime)
      if (masterproc) then 
        write(6,*)'model physics parameters:'
        write(6,*)'   prognostic hydrology                   = ',hydro
        write(6,*)'   flag for random perturbation test      = ',pergro
        write(6,*)'   energy and water conservation checks   = ',conchk
        write(6,*)'   solar radiation frequency (iterations) = ',irad
        write(6,*)
      end if

* convert nelapse to nestep: this is done again in subroutine
* lsmres for a continuation run, in which case nbstep > 1

      if (nelapse .gt. 0) nestep = nbstep + nelapse - 1

      return
      end
 
