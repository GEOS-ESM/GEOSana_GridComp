#include <preproc.h>
      subroutine rst_lsm (fileout, nio, nstep)

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>
#include <basinpar.h>

* ------------------------ input variables ------------------------
      character*(*)	::	fileout		! output filename
      integer 		:: 	nio        	! logical unit number for file
      integer 		:: 	nstep          	! current time step

* ------------------------ common block variables -----------------
c
cjr Moved dynmem.h before lsmtcdyn.h because some compilers complain about
cjr arrays dimensioned some size which is defined later in the code
c
#include <dynmem.h>
#include <lsmtcdyn.h>
#include <lsmhis.h>
#include <lsmhisdyn.h>
#include <lsmctl.h>
#include <lsmtvdyn.h>
#include <lsmtc.h>
#include <cplavg.h>
#include <basin.h>
* ------------------------ local variables ------------------------
      integer 		:: 	ierr        	! error flag from open
* -----------------------------------------------------------------
* open main restart file. write data. close. 

      open (nio, file=fileout,access='sequential',
     .              form='unformatted',iostat=ierr)

      if (ierr .ne. 0) then
        write(6,*) ' lsmrst subroutine error opening ',fileout
        write(6,*) ' Will continue with run without output.'

      else

* write land surface variables needed for restart and that can be used 
* as initial data for new run

         write (nio) numlon
         write (nio) latixy
         write (nio) longxy
         write (nio) surf2d
         write (nio) soic2d
         write (nio) sand2d
         write (nio) silt2d
         write (nio) clay2d
         write (nio) pctlak
         write (nio) pctwet
         write (nio) drngarea
         write (nio) drngfrac
         write (nio) drngbasn
         write (nio) h2osno
         write (nio) h2ocan
         write (nio) h2osoi
         write (nio) tv    
         write (nio) tg    
         write (nio) tsoi  

* write other main land surface variables needed for restart

         write (nio) moz, eah, soot, hsno, fsno, fwet

* write vegetation for next time step

         write (nio) htop, tlai  ,tsai,  elai , esai,  
     $               foln, stemb ,rootb, soilc, igs

* write albedo for next time step

         write (nio ) albd ,albi ,albgrd ,albgri ,fabd ,
     $                fabi ,ftdd ,ftid   ,ftii   ,fsun    

* write history file related info

         write (nio) nstep,nfil,ntim
         write (nio) nmon,nyr
         write (nio) mcdate_f,mcsec_f,mdcur_f,mscur_f
         write (nio) mcdate_t,mcsec_t,mdcur_t,mscur_t  
         write (nio) slfval,mlfval
         write (nio) slfcnt,mlfcnt
         write (nio) bashist
         write (nio) bascnt
         write (nio) flondat,flnddat, finidat

#ifndef COUP_CCM
* write flux average info
         write (nio) flxave
         if (flxave) then
            write (nio) dosend
            write (nio) icnt
            write (nio) tauxa, tauya, lhflxa, shflxa, lwupa, qflxa,
     $                  drnveca
         endif
#endif

         close(nio)
  
      endif

      return
      end


 
