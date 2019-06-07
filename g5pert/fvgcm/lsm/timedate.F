#include <preproc.h>
      subroutine timedate(cdate,ctime)

C-----------------------------------------------------------------------
C
C  Generic Date and Time routine
C
C-----------------------------------------------------------------------
      implicit none
C-----------------------------Arguments---------------------------------
C
C Output
C
      character*8 cdate       ! Output:  Current date
      character*8 ctime       ! Output:  Current time
C
C--------------------------Local Variables------------------------------

      integer i
      integer mon,day,yr,hr,min,sec

      call get_date(mon,day,yr,hr,min,sec)

      write(cdate,'(i2,a1,i2,a1,i2)')mon,'/',day,'/',yr
      write(ctime,'(i2,a1,i2,a1,i2)')hr,':',min,':',sec
      do i=1,8
         if (cdate(i:i).eq.' ') cdate(i:i) = '0'
         if (ctime(i:i).eq.' ') ctime(i:i) = '0'
      end do

      return
      end
 
