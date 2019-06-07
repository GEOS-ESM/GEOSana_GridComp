#include <preproc.h>
      subroutine lsmres(nio)

#if defined (SPMD)
      use mod_comm, only: mp_bcst_int, mp_scat3d, mp_scat3d_int
#endif

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>
#include <basinpar.h>

* ------------------------ code history ---------------------------
* source file:       lsmres.F
* purpose:           read a lsm restart file. position history file if needed
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* this code reads the lsmr_xxxx restart file, set by either the
* lsmr_rstrt file (if restart) or nrevsn (if branch run). this file
* was opened in lsmini and the surface data at the beginning of 
* the file was read. the remaining restart variables are read here 
* and the file is then closed. if a restart and the current history 
* file is not full, this code also opens the current history file 
* and positions to the end of the file so that subsequent time 
* samples are added until the file is full.
* [note that when using the mass store, this history file, if needed, 
* is copied to the current working directory not to the locpnh directory].
* a new history file is used if it is a branch run. 
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmtvdyn.h>
#include <lsmctl.h>
#include <lsmhis.h>
#include <lsmhisdyn.h>
#include <cplavg.h>
#include <basin.h>
* -----------------------------------------------------------------

* ------------------------ input variables ------------------------
      integer nio              !unit number
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i,j
      integer nstep            !time index
#ifndef COUP_CCM
      logical flxave_res       !flux averaging flag read from restart file
#endif
* -----------------------------------------------------------------

* -----------------------------------------------------------------
* read restart file
* -----------------------------------------------------------------

      if (masterproc) then
         write (6,*) 'attempting to read restart data'

* read land surface variables needed for restart and that are part of 
* dataset that can be used as initial data for new run

         read (nio) h2osno
         read (nio) h2ocan
         read (nio) h2osoi
         read (nio) tv    
         read (nio) tg    
         read (nio) tsoi  

* read other main land surface variables needed for restart

         read (nio) moz, eah, soot, hsno, fsno, fwet

* read vegetation for next time step

         read (nio) htop, tlai  ,tsai,  elai , esai,  
     $              foln, stemb ,rootb, soilc, igs

* read albedo for next time step

         read (nio ) albd ,albi ,albgrd ,albgri ,fabd ,
     &               fabi ,ftdd ,ftid   ,ftii   ,fsun    

* read history file related info

         read (nio) nstep,nfil,ntim
         read (nio) nmon,nyr
         read (nio) mcdate_f,mcsec_f,mdcur_f,mscur_f
         read (nio) mcdate_t,mcsec_t,mdcur_t,mscur_t  
         read (nio) slfval,mlfval
         read (nio) slfcnt,mlfcnt
         read (nio) bashist
         read (nio) bascnt
         read (nio) flondat,flnddat,finidat

#ifndef COUP_CCM
* read in flux averaging flag from restart dataset - this will
* only be used for error checking  in restart runs

         read (nio) flxave_res

* restart run: check to see that flux averaging is specified the
* same way in both namelist and restart data. if not, then print 
* warning and chose restart flag and read in remaining restart data.

         if (nsrest.eq.1) then
            if ((      flxave_res .and. .not. flxave) .or.
     &           (.not. flxave_res .and.       flxave)) then
               flxave = flxave_res
               write(6,*)'(lsmres): restart run WARNING '
               write(6,*)'(lsmres): flxave from namelist input  = ',
     &              flxave
               write(6,*)'(lsmres): flxave from restart dataset = ',
     &              flxave_res
               write(6,*)'(lsmres): flxave from restart dataset ',
     &              'will be used.'
            endif
            if (flxave) then
               read (nio) dosend
               read (nio) icnt
               read (nio) tauxa, tauya, lhflxa, shflxa, lwupa, qflxa,
     &                    drnveca
            endif
         endif

* branch run: if flux averaged run is branched from a non flux-averaged 
* run, check that the time step obtained from the restart dataset
* is an albedo calculation time step and then set the send flag to true.

         if (nsrest.eq.3) then
            if (flxave) then
               read (nio,end=900) dosend
               read (nio) icnt
               read (nio) tauxa, tauya, lhflxa, shflxa, lwupa, qflxa
               go to 910
 900           write(6,*)'(lsmres): EOF reached in restart dataset: ',
     &              'branching a flux averaged run from a non',
     &              'flux-averaged run'
               if (mod(nstep,irad) .ne. 0) then
               write(6,*)'(lsmres): ERROR, branching a flux averaged ',
     &            'run from a non-flux averaged run must be done on ',
     &            'an albedo calculation time step'
                  call stoprun
               endif
               dosend = .true.
 910           continue
            endif
         endif
#endif

         call relavu (nio)
      end if                    !if masterproc

* -----------------------------------------------------------------
* set beginning and ending times
* -----------------------------------------------------------------

      nbstep = nstep+1
      if (nelapse .gt. 0) nestep = nbstep + nelapse - 1

#if defined(SPMD)
* ----------------------------------------------------------------------
* end lsm initialization
* ----------------------------------------------------------------------

!
! Some vars are too large for comm buffer and are looped to send by length=kpt
!

      call mp_scat3d(h2osno, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)
      call mp_scat3d(h2ocan, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)
      call mp_scat3d(h2osoi, msl, kpt, 1, 1,msl, 1,kpt, 1,1, 0)
      call mp_scat3d(tsoi, msl, kpt, 1, 1,msl, 1,kpt, 1,1, 0)
      call mp_scat3d(tv, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)
      call mp_scat3d(tg, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)
      call mp_scat3d(moz, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(eah, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(soot, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(hsno, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(fsno, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(fwet, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(htop, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(tlai, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(tsai, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(elai, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(esai, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(foln, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(stemb, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(rootb, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(soilc, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(igs, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)      
      call mp_scat3d(albd, mband, kpt, 1, 1,mband, 1,kpt, 1,1, 0)     
      call mp_scat3d(albi, mband, kpt, 1, 1,mband, 1,kpt, 1,1, 0)
      call mp_scat3d(albgrd, mband, kpt, 1, 1,mband, 1,kpt, 1,1, 0)
      call mp_scat3d(albgri, mband, kpt, 1, 1,mband, 1,kpt, 1,1, 0)
      call mp_scat3d(fabd, mband, kpt, 1, 1,mband, 1,kpt, 1,1, 0)
      call mp_scat3d(fabi, mband, kpt, 1, 1,mband, 1,kpt, 1,1, 0)
      call mp_scat3d(ftdd, mband, kpt, 1, 1,mband, 1,kpt, 1,1, 0)
      call mp_scat3d(ftid, mband, kpt, 1, 1,mband, 1,kpt, 1,1, 0)
      call mp_scat3d(ftii, mband, kpt, 1, 1,mband, 1,kpt, 1,1, 0)
      call mp_scat3d(fsun, 1, kpt, 1, 1,1, 1,kpt, 1,1, 0)     
      call mp_bcst_int(nstep)
      call mp_bcst_int(nfil)
      call mp_bcst_int(ntim)
      call mp_bcst_int(nmon)
      call mp_bcst_int(nyr)
      call mp_bcst_int(mcdate_f)
      call mp_bcst_int(mcsec_f)
      call mp_bcst_int(mdcur_f)
      call mp_bcst_int(mscur_f)
      call mp_bcst_int(mcdate_t)
      call mp_bcst_int(mcsec_t)
      call mp_bcst_int(mdcur_t)
      call mp_bcst_int(mscur_t)
      call mp_bcst_int(bascnt)
      call mp_scat3d(bashist, ndrn, 1, 1, 1,ndrn, 1,1, 1,1, 0)          
      call mp_scat3d(slfval, kpt, mslflds, 1, 1,kpt, 1,mslflds, 1,1, 0)           
      call mp_scat3d(mlfval, kpt, msl, mmlflds, 1,kpt, 1,msl, 1,mmlflds, 0)        
      call mp_scat3d_int(slfcnt, mslflds, 1, 1, 1,mslflds, 1,1, 1,1, 0)           
      call mp_scat3d_int(mlfcnt, mmlflds, 1, 1, 1,mmlflds, 1,1, 1,1, 0)        
#endif

      if (masterproc) then
         write (6,*) 'successfully read restart data '
         write (6,*)
      end if
      return
      end
 
