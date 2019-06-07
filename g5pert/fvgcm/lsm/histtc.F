#include <preproc.h>
      subroutine histtc(beglatspmd)

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>

* ------------------------ code history ---------------------------
* source file:       histtc.F
* purpose:           write time-const data [lsmtc common block] to history file
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* the history file variables (lsmhis common block) set in this
* subroutine are overwritten by a subsequent call to histini.
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmtc.h>
#include <lsmtcdyn.h>
#include <lsmhis.h>
#include <ncd.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      real, allocatable ::  rivt(:)        !real ivt      
      real, allocatable ::  rist(:)        !real ist      
      real, allocatable ::  risc(:)        !real isc      
      integer lencs         !integer function: length of character string
      character*80 msscom   !mass store comment
      integer i,k           !loop index
      integer npt           !lenght of "little" vector
      integer kpti          !first point in "big" kpt vec for loopi "little" vec
      integer nstep         !dummy time index
      integer beglatspmd    !beg lat on proc, degenerates to 1 for nonspmd
* -----------------------------------------------------------------

!     allocate(rivt(begkptspmd(beglatspmd):endkptspmd(endlatspmd))) 
!     allocate(rist(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
!     allocate(risc(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))

      allocate(rivt(begkptspmd(beglatspmd):begkptspmd(beglatspmd)+numks-1)) 
      allocate(rist(begkptspmd(beglatspmd):begkptspmd(beglatspmd)+numks-1)) 
      allocate(risc(begkptspmd(beglatspmd):begkptspmd(beglatspmd)+numks-1)) 

      if (masterproc) then
         write (6,*) '---------------------------------------'
         write (6,*) 'creating lsm time-constant history file'
      endif	

* dummy time info

      nstep = 0
      mcdate = 0
      mcsec = 0
      mdcur = 0 
      mscur = 0
      timcom = 'time constant data in lsmtc common block'
      msscom = timcom
      ntim = 1      

* convert integer data to real

!     do k = begkptspmd(beglatspmd), endkptspmd(endlatspmd)
      do k = begkptspmd(beglatspmd), begkptspmd(beglatspmd)+numks-1
         rivt(k) = float(ivt(k))
         rist(k) = float(ist(k))
         risc(k) = float(isc(k))
      end do

* set up single-level history fields

      ninst = 'instant ' 
      nslflds = 0
      nslflds = nslflds + 1
      slfnam(nslflds)  =  'IVT     '
      slfuni(nslflds)  =  '        '
      slftyp(nslflds)  =  ninst
      slfdes(nslflds)  =  'plant type'
 
      nslflds = nslflds + 1
      slfnam(nslflds)  =  'IST     '
      slfuni(nslflds)  =  '        '
      slftyp(nslflds)  =  ninst
      slfdes(nslflds)  =  'soil type'

      nslflds = nslflds + 1
      slfnam(nslflds)  =  'ISC     '
      slfuni(nslflds)  =  '        '
      slftyp(nslflds)  =  ninst
      slfdes(nslflds)  =  'soil color class'

      nslflds = nslflds + 1
      slfnam(nslflds)  =  'WATSAT  '
      slfuni(nslflds)  =  'MM3/MM3 '
      slftyp(nslflds)  =  ninst
      slfdes(nslflds)  =  'saturated soil water content (porosity)'

      nslflds = nslflds + 1
      slfnam(nslflds)  =  'HKSAT   '
      slfuni(nslflds)  =  'MM/S    '
      slftyp(nslflds)  =  ninst
      slfdes(nslflds)  =  'saturated hydraulic conductivity'

      nslflds = nslflds + 1
      slfnam(nslflds)  =  'SMPSAT  '
      slfuni(nslflds)  =  'MM      '
      slftyp(nslflds)  =  ninst
      slfdes(nslflds)  =  'saturated soil matrix potential'

      nslflds = nslflds + 1
      slfnam(nslflds)  =  'BCH     '
      slfuni(nslflds)  =  '        '
      slftyp(nslflds)  =  ninst
      slfdes(nslflds)  =  'clapp and hornberger "b"'

      nslflds = nslflds + 1
      slfnam(nslflds)  =  'WATDRY  '
      slfuni(nslflds)  =  'MM3/MM3 '
      slftyp(nslflds)  =  ninst
      slfdes(nslflds)  =  'soil water when et stops'

      nslflds = nslflds + 1
      slfnam(nslflds)  =  'WATOPT  '
      slfuni(nslflds)  =  'MM3/MM3 '
      slftyp(nslflds)  =  ninst
      slfdes(nslflds)  =  'optimal soil water content for et'

      nslflds = nslflds + 1
      slfnam(nslflds) =  'CSOL    '
      slfuni(nslflds) =  'J/M3K   '
      slftyp(nslflds) =  ninst
      slfdes(nslflds) =  'heat capacity, soil solids'

      nslflds = nslflds + 1
      slfnam(nslflds) =  'TKSOL   '
      slfuni(nslflds) =  'W/MK    '
      slftyp(nslflds) =  ninst
      slfdes(nslflds) =  'thermal conductivity, soil solids'

      nslflds = nslflds + 1
      slfnam(nslflds) =  'TKDRY   '
      slfuni(nslflds) =  'W/MK    '
      slftyp(nslflds) =  ninst
      slfdes(nslflds) =  'thermal conductivity, dry soil'

      nslflds = nslflds + 1
      slfnam(nslflds) =  'SAND    '
      slfuni(nslflds) =  'percent '
      slftyp(nslflds) =  ninst
      slfdes(nslflds) =  'percent sand'

      nslflds = nslflds + 1
      slfnam(nslflds) =  'CLAY    '
      slfuni(nslflds) =  'percent '
      slftyp(nslflds) =  ninst
      slfdes(nslflds) =  'percent clay'

      if (nslflds .gt. mslflds) then
         if (masterproc) then
            write (6,*) 'histtc error: number single-level fields',
     &           ' > maximum allowed'
	 end if
         call stoprun
      end if

* set up multi-level history fields

      nmlflds = 0
      nmlflds = nmlflds + 1
      mlfnam(nmlflds) =  'DZSOI   '
      mlfuni(nmlflds) =  'M       '
      mlftyp(nmlflds) =  ninst
      mlfdes(nmlflds) =  'soil layer thickness'

      nmlflds = nmlflds + 1
      mlfnam(nmlflds) =  'ZSOI    '
      mlfuni(nmlflds) =  'M       '
      mlftyp(nmlflds) =  ninst
      mlfdes(nmlflds) =  'soil layer depth'

      nmlflds = nmlflds + 1
      mlfnam(nmlflds) =  'ROOT    '
      mlfuni(nmlflds) =  '        '
      mlftyp(nmlflds) =  ninst
      mlfdes(nmlflds) =  'relative root abundance'

      if (nmlflds .gt. mmlflds) then
         if (masterproc) then
            write (6,*) 'histtc error: number multi-level fields',
     &           ' > maximum allowed'
         endif
         call stoprun
      end if

* echo active fields

      if (masterproc) then
         write (6,*) 'active fields: single-level'
         write (6,1002)
         do i = 1, nslflds
            write (6,1003) i,slfnam(i),slfuni(i),slftyp(i),slfdes(i)
         end do
         write (6,*)
         
         write (6,*) 'active fields: multi-level'
         write (6,1002)
         do i = 1, nmlflds
            write (6,1003) i,mlfnam(i),mlfuni(i),mlftyp(i),mlfdes(i)
         end do
         write (6,*)
      end if
1002  format(' no',' name    ',' units   ',' type    ',' description')
1003  format((1x,i2),(1x,a8),(1x,a8),(1x,a8),(1x,a40))

* set control variables for calls to histslf and histmlf

      ehi = .true.
      do i = 1, nslflds
         slfcnt(i) = 0
      end do
      do i = 1, nmlflds
         mlfcnt(i) = 0
      end do

* accumulate field values over history time interval
!  Which is better small vectors or big vectors
!     do i = 1, 1
!        kpti = begkptspmd(beglatspmd)
!        npt = numks
      do i = 1, numlv
         kpti = begkpt(i)
         npt = numkpt(i)
         call histslf ('IVT     ',npt, rivt(kpti)    , kpti, kpti+npt-1)
         call histslf ('IST     ',npt, rist(kpti)    , kpti, kpti+npt-1)
         call histslf ('ISC     ',npt, risc(kpti)    , kpti, kpti+npt-1)
         call histslf ('WATSAT  ',npt, watsat(kpti)  , kpti, kpti+npt-1)
         call histslf ('HKSAT   ',npt, hksat(kpti)   , kpti, kpti+npt-1)
         call histslf ('SMPSAT  ',npt, smpsat(kpti)  , kpti, kpti+npt-1)
         call histslf ('BCH     ',npt, bch(kpti)     , kpti, kpti+npt-1)
         call histslf ('WATDRY  ',npt, watdry(kpti)  , kpti, kpti+npt-1)
         call histslf ('WATOPT  ',npt, watopt(kpti)  , kpti, kpti+npt-1)
         call histslf ('CSOL    ',npt, csol(kpti)    , kpti, kpti+npt-1)
         call histslf ('TKSOL   ',npt, tksol(kpti)   , kpti, kpti+npt-1)
         call histslf ('TKDRY   ',npt, tkdry(kpti)   , kpti, kpti+npt-1)
         call histslf ('SAND    ',npt, sand(kpti)    , kpti, kpti+npt-1)
         call histslf ('CLAY    ',npt, clay(kpti)    , kpti, kpti+npt-1)
         call histmlf ('DZSOI   ',npt, dzsoi(1,kpti) , kpti, kpti+npt-1)
         call histmlf ('ZSOI    ',npt, zsoi(1,kpti)  , kpti, kpti+npt-1)
         call histmlf ('ROOT    ',npt, root(1,kpti)  , kpti, kpti+npt-1)
      end do

* set 1d fields to zero
      
      nsl1dflds = 0
      nml1dflds = 0

* create time-invariant netcdf history file

      locfnh = locpnh(1:lencs(locpnh))//'/lsmh_timcon.nc'
      if (masterproc) call histcrt

* write history file

      call histwrt (nstep,beglatspmd)

* close history file and dispose to mass store

      if (masterproc) then
         write (6,*) 'closing local netcdf history file ',
     $        locfnh(1:lencs(locfnh))
      endif 
      if (masterproc) call ncclos(ncid, rcode)
      if (masterproc ) then
         call putfil (locfnh, msspnh, irt, nswrps, msscom) 
         write (6,*) msscom(1:lencs(msscom))
         write (6,*) '---------------------------------------'
         write (6,*)
      end if                    !masterproc

! Added by SJL

      deallocate(rivt)
      deallocate(rist)
      deallocate(risc)

      return 
      end
 
