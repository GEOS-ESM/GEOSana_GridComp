#include <params.h>
       subroutine iostop(iostat  ,nunit   ,nrec    ,clabel)
C-----------------------------------------------------------------------
C
C Explain the CRAY FORTRAN I/O error, then call endrun
C  
C---------------------------Code history--------------------------------
C
C Original version:  L. Buja
C Standardized:      L. Buja, Jun 1992, Feb 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer iostat         ! I/O error number from err=iostat option
      integer nrec           ! Number of current record (ignored if <=0)
      integer nunit          ! I/O Unit number
      character*(*) clabel   ! Users written diagnostic
C
C---------------------------Local variables-----------------------------
C
      integer i              ! Do loop counter
#if ( defined CRAY )
      integer iret           ! Return code for the ishell call
      character*16 iolabel   ! String to store error explanation
C
C------------------------------Externals--------------------------------
C
      integer ishell
      external ishell,endrun
#endif
C
C-----------------------------------------------------------------------
C
      if (iostat.ne.0) then
        write (6,*) 'IOSTOP:',('*',i=1,30),'  I/O ERROR  ',('*',i=1,29)
        write (6,*) '       ',clabel 
        if (nrec.ge.1) then
          write (6,*) 'I/O Unit = ',nunit,'   Record number = ',nrec
     $               ,'    Error number = ',iostat
        else
          write (6,*) 'I/O Unit = ',nunit,'   Error number = ',iostat
        end if
        if (iostat.gt.0) then
#if ( defined CRAY )
          write(iolabel(1:16),'(a12,i4)') 'explain lib ',iostat
          write (6,*) iolabel(1:16)
          iret = ishell(iolabel(1:16))
#endif  
        else
          write (6,*) 'End Of File (EOF) was encountered.'
        end if
        call endrun
      end if
C
      return
      end
 
