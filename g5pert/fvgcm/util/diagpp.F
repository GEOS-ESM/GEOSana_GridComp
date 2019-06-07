      program diagpp


c
c ... Postprocess diagnostics fields from fvgcm
c
c     1. vertical interpolateion to NCEP standard pressure levels
c     2. time average
c
      implicit         none

#_RT     include <fvgcm.h>
#     include <diag.h>

      integer          imr, jnp, nlayr

      integer          plev

#if (FVGCM_LEV >= 55) 
      parameter        (plev = 26)        ! NCEP data levels
      real             pm(plev)
      data pm /0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 3., 5.,
     .          10.   , 20.   , 30.  , 50.  , 70.   ,100.   ,150.   ,
     .         200.   ,250.   ,300.  ,400.  ,500.   ,600.   ,700.   ,
     .         850.   ,925.   ,1000.  /
#else
      parameter        (plev = 17)        ! ERA-15 data levels
      real             pm(plev)
      data pm / 10.   , 30.  , 50.  , 70.   ,100.   ,150.   ,
     .         200.   ,250.   ,300.  ,400.  ,500.   ,600.   ,700.   ,
     .         775.   ,850.   ,925.   ,1000.  /
! NCEP
!     data pm / 10.   , 20.   , 30.  , 50.  , 70.   ,100.   ,150.   ,
!    .         200.   ,250.   ,300.  ,400.  ,500.   ,600.   ,700.   ,
!    .         850.   ,925.   ,1000.  /

#endif
      real             vcoord(plev)

      real             undef
      parameter        (undef = 1.e25)

      integer          iargc
      integer          len_trim

      integer          inx(pdiag+1)
      integer          ndim(pdiag)
      integer          pdim(pdiag)
      integer          nl
      integer          ictl, nt, nvar
      integer          frame, pframe
      integer          yy, mo, dd, hh, mn, tinc
      integer          ymd(2), hms(2)
      integer          ir, iw1, iw2, n, k
      integer          nstr

      character*80     ictlflnm
      character*80     octlflnm1
      character*80     octlflnm2
      character*80     iflnm
      character*80     oflnm1
      character*80     oflnm2
      character*40     case
      character*40     type
      character*8      flag
      character*8      varname(pdiag)
      character*2      tunit
      
      logical          saveprs
      logical          exist


      n = iargc()
      if (n .lt. 1 .or. n .gt. 2) then
        print *
        print *, 'Usage: diagpp [-saveprs] filename'
        print *
        stop
      endif
c
c     iflnm: Input filename in the form of
c
c       CASE_TYPE_YYYYMMDD.hhmmss-YYYYMMDD.hhmmss
c
c      CASE: case name
c      TYPE: rout or diag
c
      saveprs = .false.

      if (n .eq. 2) then
        call Getarg (1, flag)
        if (index(flag, '-saveprs') .ne. 0) then
          saveprs = .true.
        endif
      endif

      ictl  = 35
      ir    = 36
      iw1   = 37
      iw2   = 38

      call Getarg (n, iflnm)
      inquire (file=iflnm, exist=exist)

      if (.not. exist) then
        print *, 'Input file does not exist -->  ',
     .            iflnm(:len_trim(iflnm))
        stop
      endif

      call name_decode (iflnm, case, type, ymd, hms)

      nstr = len_trim(iflnm)
      write (ictlflnm, '(a,a4)') iflnm(:nstr), '.ctl'

      call read_ctl (ictl, ictlflnm, iflnm, imr, jnp, nlayr, nt, yy, mo,
     .               dd, hh, mn, tinc, tunit, nvar, inx, varname)

      do n = 1, nvar
        ndim(n) = inx(n+1) - inx(n)
        if (ndim(n) .eq. 1) then
          pdim(n) = 1
        else
          pdim(n) = plev
        endif
      end do

      frame  = 0
      pframe = 0
      nl = ndim(1)
      do n = 1, nvar
        nl = max(nl,ndim(n))
        frame = frame + ndim(n)
        pframe = pframe + pdim(n)
      end do

      open (ir, file=iflnm, form='unformatted', status='unknown',
     .      access='direct', recl=imr*jnp*4)
c
c ... oflnm1: time-averaged  fields on pressure coordinates
c ... oflnm2: time-dependent fields on pressure coordinates (optional)
c
      write (oflnm1, '(a,a1,a,a8,a)')
     .       case(:len_trim(case)), '_', type(:len_trim(type)),
     .       '_prs_tm_', iflnm(nstr-30:nstr)
      open (iw1, file=oflnm1, form='unformatted', status='unknown',
     .      access='direct', recl=imr*jnp*pframe*4)
      write (octlflnm1, '(a,a4)') oflnm1(:len_trim(oflnm1)), '.ctl'
      do k = 1, plev
        vcoord(k) = pm(plev-k+1)
      end do
      call write_ctl (ictlflnm, octlflnm1, oflnm1, undef,
     .                plev, vcoord, 1)

      if (saveprs) then
        write (oflnm2, '(a,a1,a,a5,a)')
     .         case(:len_trim(case)), '_', type(:len_trim(type)),
     .         '_prs_', iflnm(nstr-30:nstr)
        open (iw2, file=oflnm2, form='unformatted', status='unknown',
     .        access='direct', recl=imr*jnp*4)
        write (octlflnm2, '(a,a4)') oflnm2(:len_trim(oflnm2)), '.ctl'
        call write_ctl (ictlflnm, octlflnm2, oflnm2, undef,
     .                plev, vcoord, nt)
      endif

      call postproc (ir, iw1, iw2, imr, jnp, nl, nvar, nt, ndim,
     .               frame, varname, plev, pdim, pframe, pm, undef,
     .               saveprs)

      close (ir)
      close (iw1)
      if (saveprs) then
        close (iw2)
      endif

      stop
      end

      subroutine name_decode (flnm, case, type, ymd, hms)

      implicit none

      character*(*)     flnm
      character*(*)     case
      character*(*)     type

      integer           index
      integer           len_trim
      integer           ymd(2)
      integer           hms(2)
      integer           nc, n
c
      nc = len_trim(flnm)
c
c     flnm: Filename in the form of
c
c           CASE_TYPE_YYYYMMDD.hhmmss-YYYYMMDD.hhmmss
c
      read(flnm(nc-5:nc),     '(i6.6)') hms(2)
      read(flnm(nc-14:nc-7),  '(i8.8)') ymd(2)
      read(flnm(nc-21:nc-16), '(i6.6)') hms(1)
      read(flnm(nc-30:nc-23), '(i8.8)') ymd(1)

      n = index(flnm(:nc-32), '_')

      case = flnm(:n-1)
      type = flnm(n+1:nc-32)
      
      return
      end

      subroutine postproc (ir, iw1, iw2, imr, jnp, nl, nvar, nt, ndim,
     .            frame, varname, plev, pdim, pframe, p_out, undef,
     .            saveprs)

      use precision

      implicit         none

      integer          memget
      integer          nvar, nt
      integer          imr, jnp, nl
      integer          frame, pframe
      integer          ndim(nvar)
      integer          pdim(nvar)
      integer          plev
      integer          ir, iw1, iw2
      integer          ks
      integer          n, i, j, k, l
      integer          offset
      integer          recno
      integer          iPS
      integer          iPTR
      integer          idx, idx2
      integer          wgt(imr,jnp,plev)
      integer          wgt2d(imr,jnp,frame)

      character*(*)    varname(nvar)

      real             undef

C Grid info
      real             pe(imr,nl+1,jnp)            ! could be p, p**kappa, or log(p)
      real             ak(nl+1)
      real             bk(nl+1)
      real(r8)	       akr8(nl+1)
      real(r8)	       bkr8(nl+1)
      real(r8)	       pintr8
      real(r8)	       ptopr8

      pointer          (xPE,   pe)
      pointer          (xAK,   ak)
      pointer          (xBK,   bk)

      real             sum(imr,jnp,pframe)
      real             ps(imr,jnp)
      real             f2d(imr,jnp)                   ! 2D data (horizontal)
      real             f3d(imr,nl,jnp)                ! 3D data on eta coordinates
      real             f2dp(imr,plev)                 ! 2D data (vertical)
      real             f3dp(imr,jnp,plev)             ! 3D data on pressure coordinates

      pointer          (xSUM,  sum)
      pointer          (xWGT,  wgt)
      pointer          (xWGT2d,  wgt2d)
      pointer          (xPS,   ps)
      pointer          (xF2D,  f2d)
      pointer          (xF3D,  f3d)
      pointer          (xF2DP, f2dp)
      pointer          (xF3DP, f3dp)

C Local work arrays
      real             dp1(imr,nl)
      real             dq1(imr,nl)
      real             dc(imr,nl)
      real             qr(imr,nl)
      real             ql(imr,nl)
      real             q6(imr,nl)

      pointer          (xDP1,  dp1)
      pointer          (xDQ1,  dq1)
      pointer          (xDC,   dc)
      pointer          (xQR,   qr)
      pointer          (xQL,   ql)
      pointer          (xQ6,   q6)

      real             p_out(plev)
      real             ptop, pint
      real             fmax, fmin

      logical          found
      logical          saveprs
      logical          savewgt

      xPE  = memget(imr*jnp*(nl+1)*4)
      xAK  = memget((nl+1)*4)
      xBK  = memget((nl+1)*4)

      xSUM = memget(imr*jnp*pframe*4)
      xWGT = memget(imr*jnp*plev*4)
      xWGT2d = memget(imr*jnp*frame*4)
      xF2D = memget(imr*jnp*4)
      xF3D = memget(imr*jnp*nl*4)
      xF2DP= memget(imr*plev*4)
      xF3DP= memget(imr*jnp*plev*4)

      xDP1 = memget(imr*jnp*4)
      xDQ1 = memget(imr*jnp*4)
      xDC  = memget(imr*jnp*4)
      xQR  = memget(imr*jnp*4)
      xQL  = memget(imr*jnp*4)
      xQ6  = memget(imr*jnp*4)

      xPS = xF2D

      do l = 1, pframe
        do j = 1, jnp
          do i = 1, imr
            sum(i,j,l) = 0.
          enddo
        enddo
      enddo
      do k = 1, plev
        do j = 1, jnp
          do i = 1, imr
            wgt(i,j,k) = 0
          enddo
        enddo
      enddo
      do k = 1, pframe
        do j = 1, jnp
          do i = 1, imr
            wgt2d(i,j,k) = 0
          enddo
        enddo
      enddo

      do k = 1, plev 
        p_out(k) = log(p_out(k) * 1.e2)
      end do

c ... Set up the Eta-coordinate

c  Quick fix, output from set_eta kind (r8)
      call set_eta (nl, ks, ptopr8, pintr8, akr8, bkr8)
      ptop = ptopr8
      pint = pintr8
      do k = 1, nl+1
        ak(k) = akr8(k)
        bk(k) = bkr8(k)
      enddo

      found = .false.
      iPS = 0
      l = 1
      do while (.not. found .and. l .le. nvar)
        if (varname(l) .eq. 'SURFP') then
          found = .true.
          iPS = l
        endif
        l = l + 1
      end do
      if (iPS .eq. 0) then
        print *, 'Cannot find surface pressure in the dataset, stop.'
        stop
      endif

      idx2 = 1

      do n = 1, nt

c ... To perform time average on data which may contain
c     values, a "time weight" should be accumulated for each
c     grid point in the 3-D space. Since the undefined data 
c     are only created in the vertical interpolation for grid
c     points under terrain, these undefined points are the 
c     same for all 3-D fields. Therefore, saving the "time 
c     weight" from the first 3-D field is sufficient.

        savewgt = .true.
c
c ... Find surface pressure in the eta dataset first to be used for
c     vertical interpolation to pressure coordinates.

        offset = (n - 1) * frame
        read (ir, rec=offset+iPS) ps            ! SURFP

c ... Compute pe

        do j = 1, jnp
          do k = 1, nl + 1
            do i = 1, imr
              pe(i,k,j) = ak(k) + bk(k) * ps(i,j)
c
c ... Use log-P for vertical interpolation
c
              pe(i,k,j) = log(pe(i,k,j))
            enddo
          enddo
        enddo


        iPTR = offset + 1 
        idx = 1

        do l = 1, nvar                    ! loop over variables
          if (ndim(l) .gt. 1) then
            do k = 1, nl
              read (ir, rec=iPTR) f2d
              do j = 1, jnp
                do i = 1, imr
c
c ... Invert vertical eta coordinates from bottom-up to top-down
c
                  f3d(i,nl-k+1,j) = f2d(i,j)
                end do
              end do
              iPTR = iPTR + 1
            end do

            do j = 1, jnp

c ... Compute coordinate increments : thickness

              do k = 1, nl
                do i = 1, imr
                  dp1(i,k) =  pe(i,k+1,j) -  pe(i,k,j)
                  dq1(i,k) = f3d(i,k+1,j) - f3d(i,k,j)
                end do
              end do
c
c ... Vertical interpolation from eta to pressure levels
c     (both eta and pressure levels are top-down)
c
              call int_ppm (imr, nl, plev, 1, imr, pe(1,1,j),
     .                      p_out, f3d(1,1,j), f2dp, 3, dp1,
     .                      dq1, dc, ql, qr, q6, undef)

              do k = 1, plev
                do i = 1, imr
                  f3dp(i,j,k) = f2dp(i,plev-k+1)     ! bottom-up
                  if (savewgt .and. f3dp(i,j,k) .ne. undef) then
                    wgt(i,j,k) = wgt(i,j,k) + 1
                  endif
                end do
              end do
            end do
            savewgt = .false.
          else
            read (ir, rec=iPTR) f2d
            do j = 1, jnp
              do i = 1, imr
                f3dp(i,j,1) = f2d(i,j)
              end do
            end do
            iPTR = iPTR + 1
          endif
c
c ... Back to bottom-up in the output pressure coordinates
c
          do k = 1, pdim(l)
            do j = 1, jnp
              do i = 1, imr
                if (f3dp(i,j,k) .ne. undef) then
                  sum(i,j,idx) = sum(i,j,idx) + f3dp(i,j,k)
                  wgt2d(i,j,idx) = wgt2d(i,j,idx) + 1
                endif
              end do
            end do
            if (saveprs) then
              write (iw2, rec=idx2) ((f3dp(i,j,k), i=1,imr), j=1,jnp)
              idx2 = idx2 + 1
            endif
            idx  = idx  + 1
          end do

        end do            ! End variable loop

      end do              ! End time loop
      
      idx = 1
      do l = 1, nvar
        if (pdim(l) .gt. 1) then
          do k = 1, pdim(l)
            do j = 1, jnp
              do i = 1, imr
                if (wgt(i,j,k) .ne. 0.) then
                  sum(i,j,idx) = sum(i,j,idx) / real(wgt(i,j,k))
                  if (abs(sum(i,j,idx)) .lt. 1.e-25) sum(i,j,idx) = 0.
                  if (abs(sum(i,j,idx)) .gt. 1.e+15) sum(i,j,idx) = undef
                else
                  sum(i,j,idx) = undef
                endif
              enddo
            enddo
            idx  = idx  + 1
          enddo
        else
          do j = 1, jnp
            do i = 1, imr
              if (wgt2d(i,j,idx) .ne. 0.) then
                sum(i,j,idx) = sum(i,j,idx) / real(wgt2d(i,j,idx))
                if (abs(sum(i,j,idx)) .gt. 1.e+15) sum(i,j,idx) = undef
              else
                sum(i,j,idx) = undef
              endif
            enddo
          enddo
          idx  = idx  + 1
        endif
      enddo

      write (iw1, rec=1) sum

      call memfree (xPE)
      call memfree (xAK)
      call memfree (xBK)

      call memfree (xSUM)
      call memfree (xWGT)
      call memfree (xWGT2d)
      call memfree (xF2D)
      call memfree (xF3D)
      call memfree (xF2DP)
      call memfree (xF3DP)

      call memfree (xDP1)
      call memfree (xDQ1)
      call memfree (xDC)
      call memfree (xQR)
      call memfree (xQL)
      call memfree (xQ6)

      return
      end


      subroutine read_ctl (ictl, ctlflnm, iflnm, m, n, l, nt, yy, mo,
     .                     dd, hh, mn, tinc, tunit, nvar, inx, varname)
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
c *                                                                    *
c *   Calling:                                                         *
c *                parse                                               *
c *                len_trim                                            *
c *   Called by:                                                       *
c *                diagpp                                              *
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
      integer          m, n, l, nlev, i
      integer          yy, mo, dd, mn, hh, tinc
      integer          zloc, monloc, colloc
      integer          slen

      logical          done, found

      character*(*)    ctlflnm, iflnm
      character*128    line
      character*128    buf
      character*80     token(ntoken), keyword
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
          else if (keyword(1:4) .eq. 'VARS') then
            read (token(2), '(i3)') nvar
c
c --- If VARS is found, start reading fields
c
            inx(1) = 1
            do i = 1, nvar
c  ... Read the field lines below 'VARS', parse it to tokens
              read (ictl, '(a128)') line
              call parse (line, ns, separator, ntk, token)
c  ... The 1st token is field name
              keyword = token(1)
              call upper (keyword)
              varname(i) = keyword(:len_trim(keyword))
c  ... The 2nd token is number of levels
              read (token(2), '(i3)') nlev
              if (nlev .eq. 0) nlev = 1
              inx(i+1) = inx(i) + nlev
              read (token(2), '(i3)') nlev
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


      subroutine write_ctl (ictlflnm, octlflnm, oflnm, undef,
     .                      nlev, vcoord, nt)
c
c **********************************************************************
c *                                                                    *
c *   Purpose:                                                         *
c *                Create GrADS control file for diagnostics output    *
c *                at pressure levels (with or without time averaging) *
c *                by modifying the original control file.             *
c *                                                                    *
c *   On entry:                                                        *
c *                                                                    *
c *    ictlflnm    original control filename                           *
c *    octlflnm    new control filename                                *
c *       oflnm    filename for the new dataset                        *
c *       undef    undefined value                                     *
c *        nlev    number of levels in the new vertical coordinates    *
c *      vcoord    vertical coordinates of the new dataset             *
c *          nt    number of times in the new dataset                  *
c *                                                                    *
c *   On  exit:                                                        *
c *                none                                                *
c *                                                                    *
c *   Calling:                                                         *
c *                parse                                               *
c *                len_trim                                            *
c *   Called by:                                                       *
c *                diagpp                                              *
c *                                                                    *
c *   Last Modified:  Thu Mar 13 16:01:58 EST 1997                     *
c *                                                                    *
c **********************************************************************

      implicit         none

      integer          ns, ntoken 
      parameter        (ns = 2, ntoken = 20)

      integer          index
      integer          len_trim
      integer          ictl, octl
      integer          sloc, slen
      integer          ntk, nt, nlev, nvar
      integer          ios, i, n, k

      real             undef
      real             vcoord(*)

      logical          done, found

      character*(*)    ictlflnm, octlflnm, oflnm
      character*128    line
      character*80     token(ntoken), keyword
      character        separator(ns)
      character        blank, tab
c
      blank     = ' '
      tab       = char(9)
      separator(1) = blank
      separator(2) = tab
      ntk = ntoken
c
c --- Open GrADS control files
c
      ictl = 51
      octl = 52
      open (ictl, file=ictlflnm, form='formatted', status='unknown')
      open (octl, file=octlflnm, form='formatted', status='unknown')
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
          if (keyword(1:4) .eq. 'DSET') then
            write (octl, '(a7,a)') 'DSET   ', oflnm(:len_trim(oflnm))
          else if (keyword(1:5) .eq. 'UNDEF') then
            write (octl, '(a7,e12.6)') 'UNDEF  ', undef
          else if (keyword(1:4) .eq. 'ZDEF') then
            read (token(2), '(i3)') n
c  ... skip all the lines containing vertical coordinates in the input
            n = n - ntk + 3
            do while (n .gt. 0)
              read (ictl, '(a128)') line
              call parse (line, ns, separator, ntk, token)
              n = n - ntk
            end do
c  ... write the new vertical coordinates
            write (octl, '(a7,i3,a8,6(1x,f9.4))')
     .             'ZDEF   ', nlev, ' LEVELS ',
     .             (vcoord(k), k=1,min(6,nlev))
            if (nlev .gt. 6) then
              do n = 7, nlev, 6
                write (octl, '(18x,6(1x,f9.4))')
     .                (vcoord(k), k=n,min(n+5,nlev))
              end do
            endif
          else if (keyword(1:4) .eq. 'TDEF') then
            slen = len_trim(token(2))
            sloc = index(line, token(2)(:slen))
            write (octl, '(a7,i3,1x,a)') 'TDEF   ', nt,
     .             line(sloc+slen+1:len_trim(line))
          else if (keyword(1:4) .eq. 'VARS') then
            write (octl, '(a)') line(:len_trim(line))
            read (token(2), '(i3)') nvar
c
c --- VARS is found, start reading fields
c
            do i = 1, nvar
c  ... Read the field lines below 'VARS', parse it to tokens
              read (ictl, '(a128)') line
              call parse (line, ns, separator, ntk, token)
              if (token(2) .eq. '0') then
                write (octl, '(a)') line(:len_trim(line))
              else
                slen = len_trim(token(2))
                sloc = index(line, token(2)(:slen))
                write (octl, '(a8,2x,i3,1x,a)') token(1), nlev,
     .                 line(sloc+slen+1:len_trim(line))
              endif
            end do
          else
c  ... Unchanged lines
            write (octl, '(a)') line(:len_trim(line))
          endif
        else if (ios .lt. 0) then
          done = .true.
c         print '(a16,a)', 'End of file --> ', ictlflnm
        else
          done = .true.
          print '(a21,a)', 'Error in reading --> ', ictlflnm
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
      close (octl)

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


      subroutine int_ppm(im,km,kn,i1,i2,pe1,pp,q1,q2,
     &                   lmt,dp1,dq1,dc,ql,qr,q6,undef)

C This routine applies a high order monotonic interpolation algorithm
C for "inerpolating" cell-mena values to prescribed coordinate
C Input:
      real pe1(im,km+1)      !coordinate value at original grid edges
      real q1(im,km)         !original field
      real pp(kn)            !output coordinate
      real undef
C Input work arrays
      real dp1(im,km), qr(im,km), ql(im,km), q6(im,km)
      real dq1(im,km), dc(im,km)
C
C Output:
      real q2(im,kn)            ! output fields
C Local variabls:
      real qmax, qmin, qtmp
      real s

c     do k=1,km
c       do i=i1,i2
c         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
c         dq1(i,k) =  q1(i,k+1) -  q1(i,k)
c       end do
c     end do

C Construct the subgrid distribution of the original data
C accordingto the chosen lmt value
      call ppm2d (q1,ql,qr,q6,dp1,dq1,dc,im,km,i1,i2,lmt)

      do 555 i=i1,i2
         k0 = 1
      do 500 k=1,kn

C Start searching
      if(pp(k) .lt. pe1(i,1) .or. pp(k) .gt. pe1(i,km+1) ) then
         q2(i,k) = undef
      else
        do 45 l=k0,km
        if(pp(k) .le. pe1(i,L+1) .and. pp(k) .ge. pe1(i,l)) then

         if(L .eq. 1) then
            qmax = max(q1(i,L), q1(i,L+1))
            qmin = min(q1(i,L), q1(i,L+1))
         elseif(L .eq. km) then
            qmax = max(q1(i,L), q1(i,L-1))
            qmin = min(q1(i,L), q1(i,L-1))
         else
            qmax = max(q1(i,L-1),q1(i,L), q1(i,L+1))
            qmin = min(q1(i,L-1),q1(i,L), q1(i,L+1))
         endif

           s  = (pp(k)-pe1(i,L)) / dp1(i,L)
         qtmp = ql(i,l) + s*(qr(i,l) - ql(i,l) +
     &                       q6(i,l)*(1.-s)     )

         qtmp    = max(qtmp, qmin)
         q2(i,k) = min(qtmp, qmax)
         k0 = l
         goto 500
        endif
45      continue
      endif
500   continue
555   continue

      return
      end

      subroutine fmap_ppm(im,km,kn,i1,i2,pe1,prk,q1,q2,
     &          dp1, dq1, dc, ql, qr, q6, fac, undef)

C This routine is currently not used; it is for the interface between
C the analysis and the model
C Partial mapping using the piecewise parabolic finite-volume method

      parameter (R3 = 1./3., R23 = 2./3.)
C Input
      real pe1(im,km+1)      ! Cell edges of the roiginal coordinate
      real q1(im,km)         ! original cell-mena values 
      real prk(kn)           ! cell center coordinate
      real undef
      real fac               ! output cell thickness ( 0 < fac < 0.5)

C Output
      real q2(im,kn)         ! output partially mapped data
      real DPK

C Local work arrays
      real dp1(im,km),qr(im,km), ql(im,km), q6(im,km)
      real dq1(im,km),dc(im,km)
C
      do k=1,kn
      do i=i1,i2
         q2(i,k) = undef
      enddo
      enddo

C fac should be < 0.5

c     do k=1,km
c       do i=i1,i2
c         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
c         dq1(i,k) =  q1(i,k+1) -  q1(i,k)
c       end do
c     end do

C Construct the monotonic subgrid distribution of the original data
      call ppm2d (q1,ql,qr,q6,dp1,dq1,dc,im,km,i1,i2,0)

      do 555 i=i1,i2
         k0 = 1
      do 555 k=1,kn

      if(prk(k) .gt. pe1(i,km+1)) go to 555  

         if(k .eq. 1) then
            dpk =  fac*min( prk(1)-pe1(i,1), prk(2)-prk(1) )
         elseif(k .eq. kn) then
            dpk = fac*min(prk(kn)-prk(kn-1), pe1(i,km) - prk(kn))
         else
            dpk = fac*min( prk(k)-prk(k-1), prk(k+1)-prk(k) )
         endif

         p_hi = prk(k) - DPK
         p_lo = prk(k) + DPK

      if(p_lo .gt. pe1(i,km+1)) then
         p_lo = pe1(i,km+1)
	   p_hi = prk(k) - max(p_lo - prk(k), 0.1)
      endif
         dp2 = p_lo - p_hi

         do 45 L=k0,km
         if( p_hi.ge.pe1(i,L) .and. p_hi.le. pe1(i,L+1)) then
             k0 = L
             PL = (p_hi-pe1(i,L)) / dp1(i,L)

              
             if(p_lo .le. pe1(i,L+1)) then
C within one cell
               PR = (p_lo-pe1(i,L)) / dp1(i,L)
               TT = R3*(PR*(PR+PL)+PL**2)
               q2(i,k) = QL(i,L) + 0.5*(Q6(i,L)+QR(i,L)
     &                   - QL(i,L))*(PR+PL) - Q6(i,L)*TT
              GOTO 555
             else
C
C Fractional area...
C
              delp = pe1(i,L+1) - p_hi
              TT   = R3*(1.+PL*(1.+PL))
              qsum = delp*(QL(i,L)+0.5*(Q6(i,L)+
     &               QR(i,L)-QL(i,L))*(1.+PL)-Q6(i,L)*TT)
              k1 = L + 1
             GOTO 111
             endif
         endif
45       continue

111      continue
         do 55 L=k1,km
         if( p_lo .gt. pe1(i,L+1) ) then
C
C Whole layer..
C
            qsum = qsum + dp1(i,L)*q1(i,L)
         else
           delp = p_lo - pe1(i,L)
           esl  = delp / dp1(i,L)
           qsum = qsum + delp * ( QL(i,L)+0.5*esl*
     &            (QR(i,L)-QL(i,L)+Q6(i,L)*(1.-R23*esl)) )
           k0 = L
           GOTO 123
         endif
55       continue
123   q2(i,k) = qsum / dp2
555   continue

      return
      end

C****6***0*********0*********0*********0*********0*********0**********72
      subroutine ppm2d(p,al,ar,a6,delp,delq,dc,im,km,i1,i2,lmt)
C****6***0*********0*********0*********0*********0*********0**********72
      real p(im,km),al(im,km),ar(im,km),a6(im,km),
     &     delp(im,km)
      real dc(im,km), delq(im,km)
 
      km1 = km - 1
 
      do 500 k=2,km
      do 500 i=i1,i2
500   A6(i,k) = delp(i,k-1) + delp(i,k)

c     do 1000 k=1,km1
c     do 1000 i=i1,i2
c     delq(i,k) = P(i,k+1) - P(i,k)
c1000  continue
 
      DO 1220 k=2,km1
      DO 1220 i=i1,i2
      c1 = (delp(i,k-1)+0.5*delp(i,k))/A6(i,k+1)
      c2 = (delp(i,k+1)+0.5*delp(i,k))/A6(i,k)
      tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /
     &                              (A6(i,k)+delp(i,k+1))
      Qmax = max(P(i,k-1),P(i,k),P(i,k+1)) - P(i,k)
      Qmin = P(i,k) - min(P(i,k-1),P(i,k),P(i,k+1))
      dc(i,k) = sign(min(abs(tmp),Qmax,Qmin), tmp)
1220  continue
 
C****6***0*********0*********0*********0*********0*********0**********72
C 4th order interpolation of the provisional cell edge value
C****6***0*********0*********0*********0*********0*********0**********72
 
      DO 12 k=3,km1
      DO 12 i=i1,i2
      c1 = delq(i,k-1)*delp(i,k-1) / A6(i,k)
      A1 = A6(i,k-1) / (A6(i,k) + delp(i,k-1))
      A2 = A6(i,k+1) / (A6(i,k) + delp(i,k))
      AL(i,k) = P(i,k-1) + c1 + 2./(A6(i,k-1)+A6(i,k+1)) *
     &          ( delp(i,k  )*(c1*(A1 - A2)+A2*DC(i,k-1)) -
     &                          delp(i,k-1)*A1*DC(i,k  ) )
12    CONTINUE
 
 
C three-cell parabolic subgrid distribution at model top
 
      DO 10 i=i1,i2
C three-cell PP-distribution
C Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
C a3 = a / 3
C b2 = b / 2
      S1 = delp(i,1)
      S2 = delp(i,2) + S1
C
      S3 = delp(i,2) + delp(i,3)
      S4 = S3 + delp(i,4)
      SS3 =  S3 + S1
      S32 = S3*S3
      S42 = S4*S4
      S34 = S3*S4
C model top
      a3 = (delq(i,2) - delq(i,1)*S3/S2) / (S3*SS3)
C
      if(abs(a3) .gt. 1.E-14) then
         b2 =  delq(i,1)/S2 - a3*(S1+S2)
         SC = -b2/(3.*a3)
         if(SC .lt. 0. .or. SC .gt. S1) then
             AL(i,1) = P(i,1) - S1*(a3*S1 + b2)
         else
             AL(i,1) = P(i,1) - delq(i,1)*S1/S2
         endif
      else
C Linear
         AL(i,1) = P(i,1) - delq(i,1)*S1/S2
      endif
      DC(i,1) = P(i,1) - AL(i,1)
C compute coef. for the off-centered area preserving cubic poly.
      DM = delp(i,1) / (S34*SS3*(delp(i,2)+S3)*(S4+delp(i,1)))
      F1 = delp(i,2)*S34 / ( S2*SS3*(S4+delp(i,1)) )
      F2 = (delp(i,2)+S3) * (SS3*(delp(i,2)*S3+S34+delp(i,2)*S4)
     &      + S42*(delp(i,2)+S3+S32/S2))
      F3 = -delp(i,2)*( SS3*(S32*(S3+S4)/(S4-delp(i,2))
     &      + (delp(i,2)*S3+S34+delp(i,2)*S4))
     &      + S42*(delp(i,2)+S3) )
      F4 = SS3*delp(i,2)*S32*(delp(i,2)+S3) / (S4-delp(i,2))
      AL(i,2) = F1*P(i,1)+(F2*P(i,2)+F3*P(i,3)+F4*P(i,4))*DM
C****6***0*********0*********0*********0*********0*********0**********72
C No over- and undershoot condition
      Cmax = max(P(i,1), P(i,2))
      Cmin = min(P(i,1), P(i,2))
      AL(i,2) = max(Cmin,AL(i,2))
      AL(i,2) = min(Cmax,AL(i,2))
10    continue

C****6***0*********0*********0*********0*********0*********0**********72

C Bottom
C Area preserving cubic with 2nd deriv. = 0 at the surface
      DO 15 i=i1,i2
      d1 = delp(i,km )
      d2 = delp(i,km1)
      qm = (d2*P(i,km)+d1*P(i,km1)) / (d1+d2)
      dq = 2.*(P(i,km1)-P(i,km)) / (d1+d2)
      c1 = (AL(i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
      AL(i,km) = qm - c1*d1*d2*(d2+3.*d1)
      AR(i,km) = d1*(8.*c1*d1**2-c3) + AL(i,km)
      DC(i,km) = AR(i,km) -  P(i,km)
C****6***0*********0*********0*********0*********0*********0**********72
C No over- and undershoot condition
      Cmax = max(P(i,km), P(i,km1))
      Cmin = min(P(i,km), P(i,km1))
      AL(i,km) = max(Cmin,AL(i,km))
      AL(i,km) = min(Cmax,AL(i,km))
C****6***0*********0*********0*********0*********0*********0**********72
15    continue

      do 20 k=1,km1
      do 20 i=i1,i2
      AR(i,k) = AL(i,k+1)
20    continue
C
C f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
C
      do 30 k=1,km
      do 30 i=i1,i2
      A6(i,k) = 3.*(P(i,k)+P(i,k) - (AL(i,k)+AR(i,k)))
30    continue
 
      if(LMT.LE.2) then
	len = i2-i1+1
         do k=1,km
         call lmppm(DC(i1,k),A6(i1,k),AR(i1,k),AL(i1,k),
     &               P(i1,k),len,LMT)
         enddo
      endif
      return
      end

C****6***0*********0*********0*********0*********0*********0**********72
      subroutine lmppm(DM,A6,AR,AL,P,IM,LMT)
C****6***0*********0*********0*********0*********0*********0**********72
      parameter ( R12 = 1./12. )
      REAL A6(IM),AR(IM),AL(IM),P(IM),DM(IM)
 
C LMT = 0: full monotonicity
C LMT = 1: semi-monotonic constraint (no undershoot)
C LMT = 2: positive-definite constraint
 
      if(LMT.eq.0) then
C Full constraint
      do 100 i=1,IM
      if(DM(i).eq.0.) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      else
      da1  = AR(i) - AL(i)
      da2  = da1**2
      A6DA = A6(i)*da1
      if(A6DA .lt. -da2) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      elseif(A6DA .gt. da2) then
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
      endif
100   continue
      elseif(LMT.eq.1) then
C Semi-monotonic constraint
      do 150 i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) go to 150
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
150   continue
      elseif(LMT.eq.2) then
C Positive definite constraint
      do 250 i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) go to 250
      fmin = p(i) + 0.25*(AR(i)-AL(i))**2/A6(i) + A6(i)*R12
      if(fmin.ge.0.) go to 250
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
250   continue
      endif
      return
      end


#if defined( CRAY_T3E )

      subroutine Getarg (argnum, arg)

      implicit none

      integer       :: argnum
      character*128 :: arg

      integer :: ierror
      integer :: ilen

      call Pxfgetarg (argnum, arg, ilen, ierror)

      if (ierror /= 0) then
        Write (6,*) 'Problem in Getarg/CRAY_T3E.'
        Stop
      end if

      return

      end

#endif

