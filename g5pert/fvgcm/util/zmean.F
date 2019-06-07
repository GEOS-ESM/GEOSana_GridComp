      program          zmean
c
c Compiling -
c
c SGI: f77 -O -n32 -mips4 -bytereclen -o zmean zmean.F
c PGI: pgf90 -O -byteswapio -tp p6 -o zmean zmean.F
c
      implicit         none

      integer          imr, jnp, nl, nc
      parameter        (nc  = 1)

      integer          maxfld
      parameter        (maxfld = 100)

      integer          iargc
      integer          max
      integer          len_trim

      integer          inx(maxfld+1)
      integer          ndim(maxfld)
      integer          ictl, nt, nvar, frame
      integer          yy, mo, dd, hh, mn, tinc
      integer          ir, iw, m, n

      real             undef

      character*80     ctlflnm
      character*80     iflnm, oflnm
      character*8      varname(maxfld)
      character*8      start_rec
      character*2      tunit

      n = iargc()
      if (n .ne. 1) then
        print *
        print *, 'Usage: zmean ctl_flnm'
        print *
        stop
      endif
      call getarg (1, ctlflnm)

      ictl = 35
      ir   = 36
      iw   = 37

      call readctl (ictl, ctlflnm, iflnm, imr, jnp, nl, nt, yy, mo,
     .              dd, hh, mn, tinc, tunit, nvar, inx, varname,
     .              undef)

      if (nvar .gt. maxfld) then
        print *, 'Number of fields exceeds maxfld, stop'
        stop
      endif

      n = len_trim(iflnm)+1
      m = index(iflnm, '.bin')
      if (m .ne. 0) n = m
      m = index(iflnm, '.grd')
      if (m .ne. 0) n = m
     
      oflnm = iflnm(:n-1) // '_zm' // iflnm(n:)

      do n = 1, nvar
        ndim(n) = inx(n+1) - inx(n)
      enddo

c     print *, 'yy,mo,dd,hh,mn,tinc,tunit=',
c    .          yy,mo,dd,hh,mn,tinc,tunit

      frame = inx(nvar+1) - 1
c     print *, 'frame=', frame
c     print *, 'iflnm=(', iflnm(:len_trim(iflnm)), ')'
c     print *, 'oflnm=(', oflnm(:len_trim(oflnm)), ')'
c     print *, 'nt=', nt
c     print *, 'nvar=', nvar
c     print *, 'undef=', undef

c     do n = 1, nvar
c       print '(a,a8,a,3x,i6)', 'varname=(', varname(n), ')', ndim(n)
c     enddo

      open (ir, file=iflnm, form='unformatted', status='unknown',
     .      access='direct', recl=imr*jnp*frame*4)

      open (iw, file=oflnm, form='unformatted', status='unknown',
     .      access='direct', recl=jnp*frame*4)

      call mean (ir, iw, imr, jnp, frame, nt, undef)

      close (ir)
      close (iw)

      stop
      end


      subroutine mean (ir, iw, imr, jnp, frame, nt, undef)

      implicit         none

      integer          malloc
      integer          ir, iw, imr, jnp, frame
      integer          nt
      integer          n, i, j, k
      integer          count

      real             undef
      real             buf(imr,jnp,frame)
      real             sum(jnp,frame)
      pointer          (xBUF, buf)
      pointer          (xSUM, sum)

      xBUF = malloc(imr*jnp*frame*4)
      xSUM = malloc(jnp*frame*4)

      do k = 1, frame
        do j = 1, jnp
          do i = 1, imr
            buf(i,j,k) = 0.
          enddo
        enddo
      enddo

      do n = 1, nt
        read (ir, rec=n) buf
        do k = 1, frame
          do j = 1, jnp
            sum(j,k) = 0.
            count = 0
            do i = 1, imr
              if (abs(buf(i,j,k)) .lt. 1.e-25) buf(i,j,k) = 0.
              if (buf(i,j,k) .ne. undef) then
                sum(j,k) = sum(j,k) + buf(i,j,k)
                count = count + 1
              endif
            enddo
            if (count .ne. 0) then
              sum(j,k) = sum(j,k) / float(count)
            else
              sum(j,k) = undef
            endif
          enddo
        enddo
        write (iw, rec=n) sum
      enddo

      call free (xBUF)
      call free (xSUM)

      return
      end


      subroutine readctl (ictl, ctlflnm, iflnm, m, n, l, nt, yy, mo,
     .                    dd, hh, mn, tinc, tunit, nvar, inx, varname,
     .                    undef)
c
c **********************************************************************
c *                                                                    *
c *   Purpose:                                                         *
c *                Parse GrADS control file (.ctl) to find out the     *
c *                information for the binary data structure.          *
c *                                                                    *
c *   On entry:                                                        *
c *                                                                    *
c *        ictl    logic unit number of the control file               *
c *     ctlflnm    filename of the control file                        *
c *                                                                    *
c *   On  exit:                                                        *
c *                                                                    *
c *           m    dimension in x direction                            *
c *           n    dimension in y direction                            *
c *           l    dimension in z direction                            *
c *          nt    number of times                                     *
c *          yy    4-digit year                                        *
c *          mo    2-digit month                                       *
c *          dd    2-digit date                                        *
c *          hh    2-digit hour                                        *
c *          mn    2-digit minute                                      *
c *        nvar    number of variables                                 *
c *         inx    array of size nvar+1 holding variable indices       *
c *     varname    array of size nvar holding variables names          *
c *       undef    value for undefined data                            *
c *                                                                    *
c *   Calling:                                                         *
c *                parse                                               *
c *                len_trim                                            *
c *   Called by:                                                       *
c *                track                                               *
c *                                                                    *
c *   Last Modified:  Thu Mar 13 16:01:58 EST 1997                     *
c *                                                                    *
c **********************************************************************

      implicit         none

      integer          ns, ntoken 
      parameter        (ns = 2, ntoken = 20)

      integer          index
      integer          len_trim
      integer          inx(*)
      integer          ictl, ntk, nt, nvar, ios
      integer          m, n, l, nlev, i, j
      integer          yy, mo, dd, mn, hh, tinc
      integer          zloc, monloc, colloc
      integer          slen

      real             undef

      logical          done, found

      character*(*)    ctlflnm, iflnm
      character*128    line
      character*80     fmt
      character*80     token(ntoken), keyword, buf
      character*8      varname(*)
      character*3      month(12)
      character*2      tunit
      character        separator(ns)
      character        blank, tab
c
      data             month /'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
     .                        'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/
c
      blank     = ' '
      tab       = char(9)
      separator(1) = blank
      separator(2) = tab
      ntk = ntoken
c
c --- Open GrADS control file
c
      open (ictl, file=ctlflnm, form='formatted', status='unknown')
c
c --- Parse the control file
c
      done = .false.
      do while (.not. done)
        read (ictl, '(a128)', iostat=ios) line
        if (ios .eq. 0) then
          call parse (line, ns, separator, ntk, token)
          keyword = token(1)
          call upper (keyword)
          if (keyword(1:4) .eq. 'DSET')then
c           print *, 'DEST=', token(2)
            iflnm = token(2)
          else if (keyword(1:4) .eq. 'XDEF') then
            read (token(2), '(i3)') m
          else if (keyword(1:4) .eq. 'YDEF') then
            read (token(2), '(i3)') n
          else if (keyword(1:4) .eq. 'ZDEF') then
            read (token(2), '(i3)') l
          else if (keyword(1:4) .eq. 'TDEF') then
            read (token(2), '(i3)') nt
            buf = token(4)
            call upper (buf)
            dd = 1
            mn = 0
            hh = 0
            zloc = index(buf, 'Z')
            colloc = index(buf, ':')
            i = 1
            mo = 0
            found = .false.
            do while (.not. found .and. i .le. 12)
              monloc = index(buf, month(i))
              if (monloc .ne. 0) then 
                mo = i
                found = .true. 
              endif
              i = i + 1
            end do
            read (buf(monloc+3:), '(i4)') yy
c           if (yy .lt. 100) then
c             if (yy .ge. 50) then
c               yy = yy + 1900
c             else
c               yy = yy + 2000
c             endif
c           endif
            if (zloc .ne. 0) then
              read (buf(zloc+1:monloc-1), '(i2)') dd
              if (colloc .ne. 0) then
                read (buf(colloc+1:zloc-1), '(i2)') mn
                read (buf(1:colloc-1), '(i2)') hh
              endif
            endif
            buf = token(5)
            call upper (buf)
            slen = len_trim(buf)
            tunit = buf(slen-1:slen)
            read (buf(1:slen-2), '(i2)')  tinc
c           print *, 'yy,mo,dd,mn,hh,tinc,tunit=',
c    .                yy,mo,dd,mn,hh,tinc,tunit
          else if (keyword(1:5) .eq. 'UNDEF') then
            buf = token(2)
            call upper (buf)
            i = index(buf, '.')
            j = index(buf, 'E')
            if (j .eq. 0) j = len_trim(buf) + 1
            write (fmt, '(a2,i2.2,a1,i2.2,a1)')
     .             '(g', len_trim(buf), '.', j-i-1, ')'
            read (token(2), fmt) undef
          else if (keyword(1:4) .eq. 'VARS') then
            read (token(2), '(i3)') nvar
c
c --- If VARS is found, start reading fields
c
            inx(1) = 1
            do i = 1, nvar
              read (ictl, '(a128)') line
              call parse (line, ns, separator, ntk, token)
              keyword = token(1)
              call upper (keyword)
              varname(i) = keyword(:len_trim(keyword))
              read (token(2), '(i3)') nlev
              if (nlev .eq. 0) nlev = 1
              inx(i+1) = inx(i) + nlev
            end do
          endif
        else if (ios .lt. 0) then
          done = .true.
c         print '(a16,a)', 'End of file --> ', ctlflnm
        else
          done = .true.
          print '(a21,a)', 'Error in reading --> ', ctlflnm
          print '(a14,i5)', 'Return code = ', ios
        endif
      end do

c     print *, 'm,n,l,nt,nvar=', m,n,l,nt,nvar
c     do i = 1, nvar
c       print *, 'varname(',i,')=', varname(i), inx(i)
c     end do
c
c --- Close the control file
c
      close (ictl)

      return
      end


      subroutine parse (string, ns, separator, nt, token)
c
c **********************************************************************
c *                                                                    *
c *   Purpose:                                                         *
c *                Parse a line of text to a number of tokens          *
c *                separated by the separator(s).                      *
c *                                                                    *
c *   On entry:                                                        *
c *                                                                    *
c *      string    input string                                        *
c *   separator    token separator(s)                                  *
c *          ns    number of token separator(s)                        *
c *                                                                    *
c *   On  exit:                                                        *
c *                                                                    *
c *          nt    number of tokens                                    *
c *       token    array of tokens                                     *
c *                                                                    *
c *   Calling:                                                         *
c *                len_trim                                            *
c *                locate                                              *
c *   Called by:                                                       *
c *                main                                                *
c *                                                                    *
c *   Last Modified:  Thu Mar 13 17:47:34 EST 1997                     *
c *                                                                    *
c **********************************************************************

      implicit         none

      integer          len_trim, locate
      integer          ns, nt
      character*(*)    string, token(*)
      character*80     fmt
      character        separator(*)
      character        comment, blank, tab
      character        delimiter
      logical          done, found, skip
      integer          k, n, ptr, head, tail

      comment   = '!'
      blank     = ' '
      tab       = char(9)
c
c --- Locate the first non-blank character
c
      head = 0
      found = .false.
      do while (.not. found )
        head = head  + 1
        if (string(head:head) .ne. blank .and.
     .      string(head:head) .ne. tab) found = .true.
      end do 
c
c --- Trim the trailing blank character and comment
c
      tail = locate(string, 1, comment)
      if (tail .eq. 0) then
        tail = len_trim(string)
      else
        found = .false.
        tail = tail + 1
        do while (.not. found )
          tail = tail - 1
          if (string(tail:tail) .ne. blank .and.
     .        string(tail:tail) .ne. tab) found = .true.
        end do 
      endif
c
c --- Search for the delimiter (first available separator)
c
      found = .false.
      ptr = head - 1
      do while (ptr .lt. tail .and. .not. found)
        ptr = ptr + 1
        n = 1
        do while (n .le. ns .and. .not. found)
          if (string(ptr:ptr) .eq. separator(n)) then
            if (string(ptr:ptr) .ne. blank .and.
     .          string(ptr:ptr) .ne. tab) then
              found = .true.
              delimiter = separator(n)
            endif
          endif
          n = n + 1
        end do
      end do
      if (.not. found) delimiter = blank
c     print *, 'found=', found
c     print *, 'delimiter=(', delimiter, ')'
c     print *, 'string(:tail)=', string(:tail)
c
c --- Trim the extra trailing non-blank separator
c
      if (string(tail:tail) .eq. delimiter) tail = tail - 1

c     print *, 'string=(', string(head:tail), ')'

      k = 0
      done = .false.
c
c --- Loop for tokens
c
      do while (.not. done)

        ptr = head
        found = .false.

        do while (.not. found .and. ptr .lt. tail)
          ptr = ptr + 1
c
c ...     look for the delimiter
c
          n = 1
          do while (n .le. ns .and. .not. found)
            if (string(ptr:ptr) .eq. separator(n)) then
              if (separator(n) .eq. delimiter) then
                found = .true.
c               print *, 'found, ptr=', ptr
              else
                if (string(ptr:ptr) .ne. blank .and.
     .              string(ptr:ptr) .ne. tab) then
                  print *, 'Syntax Error : delimiter "', 
     .                      delimiter, '" expected.'
                  print '(/a)', string(:tail)
                  write (fmt, '(a1, i2.2, a6)')
     .                   '(', ptr - 1, 'x,a1/)'
                  print fmt, '^'
                  stop
                endif
              endif
            endif
            n = n + 1
          end do
c
c ...     found a token, save it
c
          if (found) then
            k = k + 1
            token(k) = string(head:ptr-1)
c           token(k) = token(k)(:len_trim(token(k)))
c           print *, 'token(',k,')=',token(k)
c
c ...       skip to next non-separator character
c
            ptr = ptr + 1
            skip = .true.
            do while (ptr .le. tail .and. skip)
              skip = .false.
              if (string(ptr:ptr) .eq. blank .or. 
     .            string(ptr:ptr) .eq. tab) then 
                skip = .true.
c               print *, 'ptr=', ptr
                ptr = ptr + 1
              endif
            end do
            head = ptr
          endif

c         print *, 'k, head, tail, ptr, found=', 
c    .              k, head, tail, ptr, found
        end do

        if (ptr .eq. tail) then
          k = k + 1
          token(k) = string(head:ptr)
c         print *, 'Token(',k,')=',token(k)
          done = .true.
        endif

c       print *, 'K, HEAD, TAIL, PTR, FOUND=', 
c    .            k, head, tail, ptr, found

      end do

      nt = k

c     print *, 'NT =', nt
c     do n = 1, nt
c       print *, 'token(',n,')=[',token(n)(:len_trim(token(n))),']'
c     end do

      return
      end



      integer function locate (string, n, chars)
c
c **********************************************************************
c *                                                                    *
c *   Purpose:                                                         *
c *                Locate certain character(s) in the input character  *
c *                string and return the position of the character as  *
c *                the function value.                                 *
c *                                                                    *
c *   On entry:                                                        *
c *                                                                    *
c *      string    input character string                              *
c *           n    dimension of char                                   *
c *       chars    a set of character to be located in 'string'        *
c *                                                                    *
c *   On  exit:                                                        *
c *                                                                    *
c *      locate    position of first available char found in 'string'  *
c *                                                                    *
c *   Calling subroutines:                                             *
c *                                                                    *
c *                None                                                *
c *   Called by:                                                       *
c *                main                                                *
c *                                                                    *
c *   Last modified: Tue Jun 11 02:03:03 EDT 1996                      *
c *                                                                    *
c **********************************************************************

      implicit         none

      integer          len
      integer          n, i, k
      character*(*)    string
      character        chars(n)
      logical          found

      found = .false. 
      locate = 0

      i = 0
      do while (i .lt. len(string) .and. .not. found)
        i = i + 1
        k = 1
        do while (k .le. n .and. .not. found)
          if (string(i:i) .eq. chars(k)) found = .true.
          k = k + 1
        end do
      end do

      if (found) locate = i

      return
      end


      integer function len_trim (string)

************************************************************************
*                                                                      * 
*     len_trim.f - return the length of string without the trailing    * 
*                  blanks and tabs                                     * 
*                                                                      * 
*     Last Modified:  Thu Mar 13 17:47:34 EST 1997                     * 
*                                                                      * 
************************************************************************

      implicit         none

      integer          tab, blank
      parameter        (tab   =  9)
      parameter        (blank = 32)

      character*(*)    string
      integer          n, code

      len_trim = 0
      n = len(string)
      if (n .eq. 0) return

      do while (n .gt. 0)
        code = ichar(string(n:n))
        if (code .eq. blank  .or. code .eq. tab) then
          len_trim = n - 1
	  n = n - 1
        else
	  len_trim = n
          n = 0
        end if
      end do

      return 
      end


      subroutine upper(string)

************************************************************************
*                                                                      * 
*     upper.f - change lower case letter to upper case letter          *
*                                                                      *
*     George Lai Tue Jun 28 16:37:00 1994                              *
*                                                                      * 
************************************************************************

      implicit         none

      character*(*) string
      integer i, n
      integer a, z, dist

      a = ichar('a')
      z = ichar('z')
      n = len(string)
      dist = ichar('A') - a

      do i = 1,n
        if (ichar(string(i:i)) .ge. a .and.
     .      ichar(string(i:i)) .le. z) then
          string(i:i) = char(ichar(string(i:i))+dist)
        endif
      end do  

      return 
      end

