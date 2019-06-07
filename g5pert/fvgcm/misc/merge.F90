      subroutine merge (pd2d, pd3d, qname, fldname, inx)
!
! ... Merge 2 sorted tables into a single sorted table
!
!     Input: 
!            1st sorted table:  qname[    1:p2d]
!            2nd sorted table:  qname[p2d+1:p2d+p3d]
!     Output:    
!            Sorted table:      fldname[1:p2d+p3d]
!            Indices wrt qname:     inx[1:p2d+p3d]
!
      implicit none

      integer        pd2d, pd3d
      integer        i2d, i3d, i, n
      integer        inx(pd2d+pd3d)

      character*8    qname(pd2d+pd3d)
      character*8    fldname(pd2d+pd3d)

      i2d = 1
      i3d = 1
      n = 1
      do while (i2d .le. pd2d .and. i3d .le. pd3d)
        if (qname(i2d) .lt. qname(i3d+pd2d)) then
          fldname(n) = qname(i2d)
          inx(n) = i2d
          i2d = i2d + 1
        else
          fldname(n) = qname(i3d+pd2d)
          inx(n) = i3d + pd2d
          i3d = i3d + 1
        endif
        n = n + 1
      enddo

      if (i2d .gt. pd2d) then
        do i = pd2d + i3d, pd2d + pd3d
          fldname(n) = qname(i)
          inx(n) = i
          n = n + 1
        enddo
      else if (i3d .gt. pd3d) then
        do i = i2d, pd2d
          fldname(n) = qname(i)
          inx(n) = i
          n = n + 1
        enddo
      endif

      return
      end
