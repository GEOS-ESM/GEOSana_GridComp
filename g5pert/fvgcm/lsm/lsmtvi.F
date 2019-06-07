#include <preproc.h>
      subroutine lsmtvi (finidat, nio, calday, nstep, kmo, kda)

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>
#include <basinpar.h>

* ------------------------ code history ------------------------------
* source file:       lsmtvi.F
* purpose:           initialize time varying variables for lsmtv common block
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      
* reviewed:        
* --------------------------------------------------------------------

* ------------------------ input variables ---------------------------
      character*80 finidat     !initial conditions file name
      integer nio              !unit number (not used if rdlsf is false)
      real calday              !current gmt calendar day (1.00 -> 365.00)
      integer nstep            !time step index
      integer kmo              !current month (1 -> 12)
      integer kda              !current day (1 -> 31)
* --------------------------------------------------------------------

* ------------------------ common block variables --------------------
#include <dynmem.h>
#include <lsmtc.h>
#include <lsmtcdyn.h>
#include <lsmtvdyn.h>
#include <soicon.h>
#include <snocon.h>
#include <lsmhis.h>
#include <basin.h>
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      real latdeg              !latitude,  + = NH (degree)
      real londeg              !longitude, + = EH (degree)
      integer i,k,m            !indices
      logical doalb            !true if albedo time step
      character*80 msscom      !mass store comment for arbitrary init
      real pi                  !pi
      integer nfini            !unit number for init conditions output
      integer getavu           !next available unit number
* --------------------------------------------------------------------
* ----------------------------------------------------------------------
* initialize water and temperature based on:
*
*    o finidat  = 'arbitrary initialization' : arbitrary initialization
*    o finidat != 'arbitrary initialization' : read initial data set 
*
* if an initial data set is read, the file is opened in lsmini.
* the time invariant data part of the file is read in lsmini and the 
* remainder of the dataset is read here.
* ----------------------------------------------------------------------

      if (finidat .ne. 'arbitrary initialization') then

         if (masterproc) then
* read data
         read (nio) h2osno
         read (nio) h2ocan
         read (nio) h2osoi
         read (nio) tv
         read (nio) tg
         read (nio) tsoi

* close and release unit

         call relavu (nio) 
         endif
      else

* set water and temperatures to constant values

         call setarr (h2osno, kpt    ,   0.)
         call setarr (h2ocan, kpt    ,   0.)
         call setarr (h2osoi, msl*kpt,  0.3)
         call setarr (tv    , kpt    , 283.)
         call setarr (tg    , kpt    , 283.)
         call setarr (tsoi  , msl*kpt, 283.)

* adjust h2osoi and tsoi for land ice, lakes, and wetlands.
* make sure h2osoi <= watsat

         do m = 1, msl
            do k = 1, kpt
               if (ist(k) .eq. istice) then
                  h2osoi(m,k) = 1.0
                  tsoi(m,k) = 273.0
               else if (ist(k) .eq. istdlak) then
                  h2osoi(m,k) = 1.0
                  tsoi(m,k) = 277.0
               else if (ist(k) .eq. istslak) then
                  h2osoi(m,k) = 1.0
                  tsoi(m,k) = 277.0
               else if (ist(k) .eq. istwet) then
                  h2osoi(m,k) = 1.0
                  tsoi(m,k) = 277.0
               end if
               h2osoi(m,k) = min(h2osoi(m,k),watsat(k))
            end do
         end do

         if (masterproc) then
* write initial data (including surface data) out to file

         nfini = getavu()
         call opnfil ('lsmi_arbini', nfini, 'u')
         write (nfini) numlon
         write (nfini) latixy
         write (nfini) longxy
         write (nfini) surf2d
         write (nfini) soic2d
         write (nfini) sand2d
         write (nfini) silt2d
         write (nfini) clay2d
         write (nfini) pctlak
         write (nfini) pctwet
         write (nfini) drngarea
         write (nfini) drngfrac
         write (nfini) drngbasn
         write (nfini) h2osno
         write (nfini) h2ocan
         write (nfini) h2osoi
         write (nfini) tv
         write (nfini) tg
         write (nfini) tsoi
         call relavu (nfini)
         write (6,*) 
         write (6,*) 'successfully wrote arbitrary initialization file',
     $        ' lsmi_arbini'
         msscom = 'initial data file: arbitrary initialization'
         call putfil ('lsmi_arbini', msspnr, irt, nswrps, msscom)

         end if
      end if

* ----------------------------------------------------------------------
* the remaining variables are not part of the initial data set because:

*    o soot  : time-invariant
*    o moz   : need only for first iteration in [surtem]. set to neutral
*    o eah   : need for first call to [stomata]. set to arbitray high 
*              value because [stomata] limits eah to <= saturation 
*              vapor pressure. result is that vapor pressure does not 
*              affect stomatal resistance on first time step
*    o hsno  : from h2osno
*    o fsno  : from h2osno
*    o fwet  : from h2ocan
* ----------------------------------------------------------------------

* moz, eah

      call setarr (moz , kpt,      0.)
      call setarr (eah , kpt,  10000.)

* soot. antarctica: soot=0. siberia: soot=0.2e-06. rest: soot=5.0e-06

      call setarr (soot, kpt, 5.0e-06)
      pi = 4.*atan(1.) 
      do k = 1, kpt
         latdeg = lati(k) * 180./pi
         londeg = long(k) * 180./pi
         if (latdeg.lt.-60. .and. ist(k).eq.istice) soot(k) = 0. !antarctica
         if (latdeg.gt.40. .and. (londeg.gt.45. .and.            !siberia
     &         londeg.lt.165.)) soot(k) = 0.2e-06
      end do

* hsno, fsno: see subroutine snoh2o

      do k = 1, kpt
         hsno(k) = h2osno(k)/bdsno            
         fsno(k) = min( hsno(k)/hsnoc, 1. )
      end do

* ----------------------------------------------------------------------
* remaining variables are initialized by calls to ecosystem dynamics and
* albedo subroutines. mimic the call to subroutine lsm in subroutine
* lsmdrv: process the "big" vectors of [kpt] points as [numlv] "little"
* vectors of [numkpt] points. [begkpt] is the starting location of the
* [numkpt] points for the "little" vector in the "big" [kpt] vector.
* ----------------------------------------------------------------------

      doalb = .true.
      do i = 1, numlv
         k = begkpt(i)

         call ecodyn (numkpt(i),kmo     ,kda      ,ivt(k)  ,hsno(k) ,
     &                lati(k)  ,tv(k)   ,htop(k)  ,elai(k) ,esai(k) ,
     &                tlai(k)  ,tsai(k) ,igs(k)   ,stemb(k),rootb(k),
     &                foln(k)  ,soilc(k),h2ocan(k),fwet(k) ,doalb   ,
     &                nstep    ,k       ,i        )

* note: need tam for first call to suralb. use tv instead.
* does not affect restart so long as suralb not called on restart

         call suralb (numkpt(i)  ,ivt(k)     ,ist(k)     ,isc(k)    ,
     &                lati(k)    ,long(k)    ,calday     ,soot(k)   ,
     &                elai(k)    ,esai(k)    ,tv(k)      ,tg(k)     ,
     &                tv(k)      ,h2osno(k)  ,h2osoi(1,k),fwet(k)   ,
     &                albgrd(1,k),albgri(1,k),albd(1,k)  ,albi(1,k) ,
     &                fabd(1,k)  ,fabi(1,k)  ,ftdd(1,k)  ,ftid(1,k) ,
     &                ftii(1,k)  ,fsun(k)    ,fsno(k)    ,k         ,
     &                i          ,nstep      )
      end do

      return
      end

 
