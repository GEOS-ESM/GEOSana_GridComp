#include <preproc.h>
      subroutine putfil (locfn, msspn, irt, nswrps, icomnt)

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       putfil.F
* purpose:           copy local file locfn to mass store file mssfn
*                    If on NCAR Cray, Fortran mass store interface will be 
*                    used if cpp token SHELL_MSS is not defined.  For other
*                    platforms, this token needs to be defined.
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ input variables ------------------------
      character*(*) locfn  !local file name
      character*(*) msspn  !mass store path
      integer irt          !mass store retention period
      character*(*) nswrps !mass store write password
      character*(*) icomnt !mass store comment 
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      character*256 mssfn  !mass store file name   
      character*80 rwpass  !read and write passwords
      character*256 opts   !options for mswrite
      character*80 fname   !local file name without directory name  
      integer lencs        !integer function: length of character string
      integer kloc         !length of character string
      integer i            !loop index
      integer ier          !mswrite error
#if ( defined SHELL_MSS )
      character*256 cmd     !
      character*256 cmdtem  !
      character crt*4       ! Retention time as characters
#if (defined( CRAY ) || defined( CRAY_T3E ))
      integer ishell
      external ishell       ! System routine, execute shell command
#endif
#endif
* -----------------------------------------------------------------

      if (irt .eq. 0) return

      rwpass = ','//nswrps(1:lencs(nswrps))
* set mass store file name. first remove directory name from local file 
* name by starting at end and looking for first "/". then add this
* name to the mass store path name

      kloc = lencs(locfn)
      do i = kloc, 1, -1
         if (locfn(i:i).eq.'/') go to 100
      end do
      i = 0
  100 fname = locfn(i+1:kloc)
      mssfn = msspn(1:lencs(msspn))//'/'//fname(1:lencs(fname))

#if ( defined SHELL_MSS )
      write (crt,'(i4)') irt
      if (lencs(icomnt).gt.0) then
        write(cmd,'(100a)') 'mswrite -c "',icomnt(1:lencs(icomnt)),
     $     '" -t ',crt
      else
        write(cmd,'(100a)') 'mswrite ',' -t ',crt
      end if
      if (lencs(rwpass).gt.1) then
         cmd=cmd(1:lencs(cmd))//' -w '//rwpass(1:lencs(rwpass))
     $      //' '//locfn(1:lencs(locfn))//' '//mssfn(1:lencs(mssfn))
      else
         cmd=cmd(1:lencs(cmd))//' '//locfn(1:lencs(locfn))
     $      //' '//mssfn(1:lencs(mssfn))
      end if

* Put mswrite cmd in background for asynch. behavior.

      cmdtem = '('//cmd(1:lencs(cmd))//'; /bin/rm '
     $         //locfn(1:lencs(locfn))//' )&'
      write(6,*)'(PUTFIL): Issuing shell cmd:',cmdtem(1:lencs(cmdtem))
#if (defined( CRAY ) || defined( CRAY_T3E ))
      ier = ishell(cmdtem(1:lencs(cmdtem)))
      if (ier.ne.0) then
        write(6,*) '(PUTFIL): ishell(mswrite(...)) returns ', ier
        call endrun
      end if
#else
      call system(cmdtem(1:lencs(cmdtem)))
      ier = 0
#endif

#else

* set options for mswrite

      if (icomnt(1:1) .ne. ' ') then
         opts = 'CMNT='//icomnt(1:lencs(icomnt))
     &        //',FORMAT=TR,MVN=CTPUBLIC,NOWAIT,NOMAIL'
      else
         opts = 'FORMAT=TR,MVN=CTPUBLIC,NOWAIT,NOMAIL'
      end if

* set read and write passwords

      write (6,*) '(PUTFIL): disposing local lsm file ',
     &     locfn(1:lencs(locfn))
      write (6,*) '        mass store file = ',mssfn(1:lencs(mssfn))
      write (6,*) '        read,write passwords = ',
     &     rwpass(1:lencs(rwpass))
      write (6,*) '        retention time (days) = ',irt

      call mswrite (ier, locfn, mssfn, rwpass, irt, opts)

      if (ier .ne. 0) then
         write (6,*) '        unsuccessful mass store write'
         call stoprun
      end if

* remove disk copy

      call unlink(locfn(1:lencs(locfn)))
#endif

      return
      end
 
