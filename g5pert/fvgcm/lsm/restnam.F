#include <preproc.h>
      subroutine restnam 

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       restnam.F
* purpose:           generate lsm restart file names 
* date last revised: April 1997
* author:            Mariana Vertenstein
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* generate restart file names  of the form lsmr_yyyymmdd_sssss, where
* yyyy is the current year, mm is the current month, dd is the 
* current day and sssss is the seconds into the current day.
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <lsmhis.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      character*19 filext !unique file extension: lsmr_yyyymmdd_sssss
      integer lencs       !integer function to get length of character string
* -----------------------------------------------------------------

      filext = 'lsmr_yyyymmdd_sssss'
      write(filext( 6:13),'(i8.8)') mcdate
      write(filext(15:19),'(i5.5)') mcsec

* add extension to local directory to get full file name

      locfnr = locpnr(1:lencs(locpnr)) // '/'
     $     // filext(1:lencs(filext)) 

      return
      end
 
