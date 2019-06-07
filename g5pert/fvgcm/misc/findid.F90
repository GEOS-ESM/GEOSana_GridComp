      integer function findid (n, fldtbl, fld, head, tail)

      implicit none

      integer     n, i
      integer     head, tail, mid
      character*8 fldtbl(n), fld
      logical     found

      findid = 0
      found = .false.
      
      do while (.not. found .and. tail .ge. head)
        mid = (head + tail) / 2
        if (fld .gt. fldtbl(mid)) then
          head = mid + 1
        else if (fld .lt. fldtbl(mid)) then
          tail = mid - 1
        else
          found  = .true.
          findid = mid
        endif
      enddo

      return
      end
