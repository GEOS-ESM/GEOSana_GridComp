#include <preproc.h>
      subroutine getfil (fulpath, locfn)

#include <vartyp.h>

* ------------------------ code history ------------------------------
* source file:       getfil.F
* purpose:           get local copy [locfn] using full pathname [fulpath]. 
*                    use local file if it exists in working directory.
*                    otherwise check for permanent disk pathname and then
*                    mass store
* date last revised: April 1997
* author:            Gordon Bonan
* revision           Mariana Vertenstein 
* standardized:
* reviewed:
* --------------------------------------------------------------------

* ------------------------ input/output variables --------------------
      character*(*) fulpath   !mass store or permanent disk full pathname
      character*(*) locfn     !output local file name
* --------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <lsmctl.h>
#include <lsmpar.h>
#include <lsmspmd.h>
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      integer i               !loop index
      integer lencs           !length of character string
      integer klen            !length of fulpath character string
      integer ioe             !error return from msread
      logical lexist          !true if local file exists
      character*256 text      !mswrite command
#if (defined( CRAY ) || defined( CRAY_T3E ))
      integer ishell
      external ishell       ! System routine, execute shell command
#endif
* --------------------------------------------------------------------
* get local file name from full name: start at end. look for first "/"

      klen = lencs(fulpath)
      do i = klen, 1, -1
         if (fulpath(i:i).eq.'/') go to 100
      end do
      i = 0
  100 locfn = fulpath(i+1:klen)
      if( masterproc ) write(6,*)'(GETFIL): attempting to find local file ',
     $     locfn(1:lencs(locfn))

* first check if file is in current working directory.

      inquire (file=locfn,exist=lexist)
      if (lexist) then
         if(masterproc) write (6,*) '(GETFIL): using ',locfn(1:lencs(locfn)),
     &        ' in current working directory'
         return
      endif

* second check on permanent disk 

      inquire(file=fulpath,exist=lexist)
      if (lexist) then
         locfn = fulpath(1:lencs(fulpath))
         if(masterproc) write(6,*)'(GETFIL): using ',fulpath(1:lencs(fulpath))
         return
      endif

* finally check mass store
      
      text='msread '//locfn(1:lencs(locfn))//' '
     $     //fulpath(1:lencs(fulpath))
      if(masterproc) write(6,*)'(GETFIL): Calling msread as follows:  ',
     $     text(1:lencs(text))
      call system(text(1:lencs(text)))
      ioe = 0

      return
      end


 
