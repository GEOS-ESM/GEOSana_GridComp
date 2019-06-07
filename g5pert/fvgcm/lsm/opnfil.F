#include <preproc.h>
      subroutine opnfil (locfn, iun, form)

#include <vartyp.h>

* ------------------------ code history ------------------------------
* source file:       opnfil.F
* purpose:           open file locfn to unit iun 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* --------------------------------------------------------------------

* ------------------------ input variables ---------------------------
      character*(*) locfn     !file name
      integer iun             !fortran unit number
      character*1 form        !file format: u = unformatted. f = formatted
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      integer ioe             !error return from fortran open
      integer lencs           !integer function: length of character string
      integer len             !length of character string
      character*11 ft         !format type: formatted. unformatted
* --------------------------------------------------------------------

      len = lencs(locfn)
      if (form.eq.'u' .or. form.eq.'U') then
         ft = 'unformatted'
      else
         ft = 'formatted  '
      end if
      open (unit=iun,file=locfn,status='unknown',form=ft,iostat=ioe)
      if (ioe .ne. 0) then
         write (6,*) 'opnfil: failed to open file ',locfn(1:len),
     &   ' to unit ',iun,' ierr=',ioe
         call stoprun
      end if
 
      return
      end
 
