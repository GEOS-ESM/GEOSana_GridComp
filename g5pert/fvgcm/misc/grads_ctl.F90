      subroutine grads_ctl (nctl, ctlflnm, oflnm, title, ix, jx, kx,  &
                            nt, undef, date, time, idt,               &
                            pstart, pend, maxdiag,                    & 
                            fldname, pick, vcoord, vdim, unit, label) 
!
! **********************************************************************
! *                                                                    *
! *   Purpose:                                                         *
! *                Create a GrADS control/description file.            *
! *                                                                    *
! *   On entry:                                                        *
! *                                                                    *
! *        nctl    logic unit number for control file                  *
! *     ctlflnm    control/description filename                        *
! *       oflnm    output data filename                                *
! *       title    description of data set                             *
! *          ix    dimension in x direction                            *
! *          jx    dimension in y direction                            *
! *          kx    dimension in z direction                            *
! *          nt    dimension in time axis (# of time steps)            *
! *       ndiag    number of selected diagnostic variables             *
! *     maxdiag    maximum number of diagnostic variables              *
! *       undef    undefined number in the data set (missing data)     *
! *        slat    latitude of the south boundary (degree)             *
! *        dlat    latitude increment (degree)                         *
! *        wlon    longitude of the west boundary (degree)             *
! *        dlon    longitude increment (degree)                        *
! *        date    initial date (YYYYMMDD)                             *
! *        time    initial time (HHMMSS)                               *
! *         idt    time increment (in seconds)                         *
! *     fldname    variable names abbreviation                         *
! *        pick    flag for selected fields                            *
! *        vdim    vertical dimension for fields                       *
! *        unit    unit of the variables                               *
! *       label    description of the variables                        *
! *                                                                    *
! *   On  exit:                                                        *
! *                none                                                *
! *                                                                    *
! *   Calling subroutines:                                             *
! *                                                                    *
! *                len_trim                                            *
! *                trimleft                                            *
! *   Called by:                                                       *
! *                fvgcm                                               *
! *                                                                    *
! *   Last Modified: Tue Dec  5 18:18:26 EST 1995                      *
! *                  WS 00.07.08  Use precision module                 *
! **********************************************************************

      use precision
      implicit          none

      character*(*)     ctlflnm, oflnm, title
      character*(*)     fldname(*)
      character*(*)     label(*)
      character*(*)     unit(*)
      character*256     line, zline
      character*80      fmt
      character*80      desc
      character*8       buf 
      character*3       month(12)

      external          len_trim

      integer(i4)       len_trim
      integer(i4)       nctl, ix, jx, kx, nt
      integer(i4)       pstart, pend
      integer(i4)       maxdiag, ndiag
      integer(i4)       pick(*)
      integer(i4)       vdim(*), verdim
      integer(i4)       date, time, idt
      integer(i4)       yy, mm, dd, hh, mn, ss
      integer(i4)       len1, len2
      integer(i4)       k, n, inx, nc

      real(r8)          vcoord(*)
      real(r8)          undef
      real(r8)          wlon, dlon, slat, dlat

      data              month /'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',  &
                               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/

      open (nctl, file=ctlflnm, form='formatted', status='unknown')

      nc = len_trim(oflnm)
      write (fmt, '(a5,i3.3,a1)') '(a7,a', nc, ')'
      write (line, fmt)  'DSET   ', oflnm(:nc)
      write (nctl, '(a)') line(:len_trim(line))

      nc = len_trim(title)
      write (fmt, '(a5,i3.3,a1)') '(a7,a', nc, ')'
      write (line, fmt)  'TITLE  ', title(:nc)
      write (nctl, '(a)') line(:len_trim(line))

      write (line, '(a8,a10)') 'OPTIONS ', 'big_endian'
      write (nctl, '(a)') line(:len_trim(line))

      write (line, '(a7,e12.6)') 'UNDEF  ', undef
      write (nctl, '(a)') line(:len_trim(line))

      wlon = 0.
      dlon = 360. / real(ix)
      write (line, '(a7,i4,a8,2(1x,f12.6))')             &
             'XDEF   ', ix, ' LINEAR ', wlon, dlon
      write (nctl, '(a)') line(:len_trim(line))

      slat = -90.
      dlat = 180. / real(jx - 1)
      write (line, '(a7,i4,a8,2(1x,f12.6))')             &
             'YDEF   ', jx, ' LINEAR ', slat, dlat
      write (nctl, '(a)') line(:len_trim(line))

!     write (line, '(a7,i4,a8,2(1x,i2))')                &
!           'ZDEF   ', kx, ' LINEAR ', 1, 1
! SJL     write (line, '(a7,i3,a8))')
      write (line, '(a7,i4,a8)')                         &
            'ZDEF   ', kx, ' LEVELS '
      write (zline, '(a18,6(1x,f9.4))') line, (vcoord(k), k=1,min(6,kx))
      write (nctl, '(a)') zline(:len_trim(zline))
      if (kx .gt. 6) then
        do n = 7, kx, 6
          write (zline, '(18x,6(1x,f9.4))') (vcoord(k), k=n,min(n+5,kx))
          write (nctl, '(a)') zline(:len_trim(zline))
        end do
      endif

      write (line, '(i8.8)')  date
      read  (line, '(i4.4,2i2.2)') yy, mm, dd
      write (line, '(i6.6)')  time
      read  (line, '(3i2.2)') hh, mn, ss

      if (idt .ne. 0) then
        if (mod(idt,86400) .eq. 0) then
          write (line,                                   &
                 '(a7,i4,a8,i2.2,a1,i2.2,a1,i2.2,a3,i4.4,3x,i4,a2)')  &
                 'TDEF   ', nt, ' LINEAR ', hh, ':', mn, &
                 'Z', dd, month(mm), yy, idt/86400, 'dy'
        else if (mod(idt,3600) .eq. 0) then
          write (line,                                   &
                 '(a7,i4,a8,i2.2,a1,i2.2,a1,i2.2,a3,i4.4,3x,i2,a2)') &
                 'TDEF   ', nt, ' LINEAR ', hh, ':', mn, &
                 'Z', dd, month(mm), yy, idt/3600, 'hr'                         
        else if (mod(idt,60) .eq. 0) then
          write (line,                                   &
                 '(a7,i4,a8,i2.2,a1,i2.2,a1,i2.2,a3,i4.4,3x,i2,a2)') &
                 'TDEF   ', nt, ' LINEAR ', hh, ':', mn, &
                 'Z', dd, month(mm), yy, idt/60, 'mn'
        else

! ... GrADS doesn't take zero time increment

          write (line,                                   &
                 '(a7,i4,a8,i2.2,a1,i2.2,a1,i2.2,a3,i4.4,3x,i1,a2)') &
                 'TDEF   ', nt, ' LINEAR ', hh, ':', mn, &
                 'Z', dd, month(mm), yy, 1, 'mn'
        endif
      else
        write (line,                                     &
               '(a7,i4,a8,i2.2,a1,i2.2,a1,i2.2,a3,i4.4,3x,i1,a2)') &
               'TDEF   ', nt, ' LINEAR ', hh, ':', mn,   &
               'Z', dd, month(mm), yy, 1, 'mo'
      endif
      write (nctl, '(a)') line(:len_trim(line))

      ndiag = 0
      do n = pstart, pend
        if (pick(n) .eq. 1) ndiag = ndiag + 1
      end do
      write (line, '(a10,i3)') 'VARS      ', ndiag
      write (nctl, '(a)') line(:len_trim(line))

      do n = pstart, pend
        if (pick(n) .eq. 1) then
          len1 = len_trim(label(n))
          len2 = max(1,len_trim(unit(n)))
          write (fmt, '(a2,i2.2,a5,i2.2,a4)')           &
                '(a', len1, ',a2,a', min(len2,77-len1), ',a1)'
          write (desc, fmt) label(n), ' [', unit(n), ']'
          if (vdim(n) .eq. 1) then
            verdim = 0
          else
            verdim = vdim(n)
          endif
          write (line, '(a8,2x,i3,a4,a)')               &
                 fldname(n), verdim, ' 99 ', desc(:len_trim(desc))
          write (nctl, '(a)') line(:len_trim(line))
        endif
      end do

      write (nctl, '(a7)') 'ENDVARS'

      close (nctl)

      return
      end 
