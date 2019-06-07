      program  interpter
!
! Compiling -
!
!=======
! Linux:
!=======
! To compile using Lahey/Fujitsu f95:
!    lf95 -fw  -I/usr/local/netcdf/include -o surf surf.F -L/usr/local/netcdf/lib -lnetcdf
! To run: surf -Wl,-T
!

! Programmer: S.-J. Lin

      implicit none

      integer im, jm

      character*80     oflnm
      character*80     hgtflnm
 
! Get dimension

      write(*,*) 'IM (E-W dimension) =?', im
      read(*,*) im
      write(*,*) 'JM (N-S dimension) =?', jm
      read(*,*) jm

      if(im .eq. 1000) then
         oflnm    = 'surf_r.data_1000x721'
         hgtflnm  = 'topo_r.1000x721.bin'                 ! Grads readable data
      elseif(im .eq. 576) then
         oflnm    = 'surf_r.data_576x361'
         hgtflnm  = 'topo_r.576x361.bin'                 ! Grads readable data
      elseif(im .eq. 288) then
         oflnm    = 'surf_r.data_288x181'
         hgtflnm  = 'topo_r.288x181.bin'
      elseif(im .eq. 144) then
         oflnm    = 'surf_r.data_144x91'
         hgtflnm  = 'topo_r.144x91.bin'
      elseif(im .eq. 120) then
         oflnm    = 'surf_r.data_120x73'
         hgtflnm  = 'topo_r.120x73.bin'
      elseif(im .eq. 100) then
         oflnm    = 'surf_r.data_100x61'
         hgtflnm  = 'topo_r.100x61.bin'
      elseif(im .eq.  72) then
         oflnm    = 'surf_r.data_72x46'
         hgtflnm  = 'topo_r.72x46.bin'
      elseif(im .eq.  24) then
         oflnm    = 'surf_r.data_24x19'
         hgtflnm  = 'topo_r.24x19.bin'
      else
         write(6,*) 'Please specify file names for the
     &               chosen output resolution'
         stop
      endif


      call surfdrv(im, jm, oflnm, hgtflnm)
      end



      subroutine surfdrv(mlon, mlat, oflnm, hgtflnm)


      implicit         none
#include <netcdf.inc>
      integer          nlon
      integer          nlat
      parameter       (nlon = 2160)
      parameter       (nlat = 1080)

      integer          mlon
      integer          mlat
      integer  mphis, msgh
      real fmin, fmax, vmax

      double precision undef
      parameter       (undef = 1.d30)

      integer          ncid
      integer          lonid
      integer          latid
      integer          ftopoid
      integer          htopoid
      integer          londim
      integer          latdim
      integer          status
      integer          n
      integer          i
      integer          j

      real             ftopo(nlon,nlat)
      real             htopo(nlon,nlat)
      real             xn(nlon)
      real             yn(nlat)
      real             buf(mlon,mlat)

      double precision sgh(mlon,mlat)
      double precision phis(mlon,mlat)
      double precision oro(mlon,mlat)
      double precision xm(mlon)
      double precision ym(mlat)

      character*80     iflnm
      character*80     oflnm
      character*80     topoflnm
      character*80     hgtflnm

      logical          check_orig
      logical          check

      real grav

! Output the original 10 min NAVY data in grads readable format
      data             check_orig /.false./
      data             check      /.true./

! **** tunable parameters *****   
      mphis = min(12, mlon/12)
      msgh  = min(12, mlon/12)
! ******************************

      iflnm    = 'topo.nc'
      topoflnm = 'topo.bin'

! ... Open input netCDF file

      status = nf_open (iflnm, NF_NOWRITE, ncid)
      if (status .ne. NF_NOERR) call handle_err(status)

      status = nf_inq_dimid (ncid, 'lon', lonid)
      if (status .ne. NF_NOERR) call handle_err(status)
      status = nf_inq_dimlen (ncid, lonid, londim)
      if (status .ne. NF_NOERR) call handle_err(status)

      status = nf_inq_dimid (ncid, 'lat', latid)
      if (status .ne. NF_NOERR) call handle_err(status)
      status = nf_inq_dimlen (ncid, latid, latdim)
      if (status .ne. NF_NOERR) call handle_err(status)

      status = nf_inq_varid (ncid, 'ftopo', ftopoid)
      if (status .ne. NF_NOERR) call handle_err(status)
      status = nf_get_var_real (ncid, ftopoid, ftopo)
      if (status .ne. NF_NOERR) call handle_err(status)

      status = nf_inq_varid (ncid, 'htopo', htopoid)
      if (status .ne. NF_NOERR) call handle_err(status)
      status = nf_get_var_real (ncid, htopoid, htopo)
      if (status .ne. NF_NOERR) call handle_err(status)

      status = nf_close (ncid)
      if (status .ne. NF_NOERR) call handle_err(status)
c
c ... Set check_orig=.true. to output original 10-minute
c         real*4 data (GrADS format)
c
      if (check_orig) then
        open (31, file=topoflnm, form='unformatted',
     .      status='unknown', access='direct', recl=nlon*nlat*4)
        write (31, rec=1) ftopo
        write (31, rec=2) htopo
        close (31)
      endif

        write(6,*) 'check original 10 min data ..'
        fmax =  vmax(htopo,fmin,nlon*nlat)
	write(6,*) 'hmax=', fmax
	write(6,*) 'hmin=', fmin

      call get_phis (nlon, nlat, xn, yn, ftopo, htopo, mlon, mlat, 
     .               xm, ym, undef, phis, oro )

      if (check) then
        open (32, file=hgtflnm, form='unformatted', status='unknown',
     .        access='direct', recl=mlon*mlat*4)

C Output un-smoothed data
            grav =  9.80616d0

        do j = 1, mlat
          do i = 1, mlon
            buf(i,j) = phis(i,j) * grav
          end do
        end do
        write (32, rec=1) buf
       endif
C SJL
c
c --> Smooth phis 

      write (*,*) 'Computing sgh'

      call get_sgh (nlon, nlat, xn, yn, ftopo, htopo, mlon, mlat, 
     .              xm, ym, undef, phis, sgh )

! *****************************   
      call sm2(mlon, mlat, phis, mphis, 0.25D0)
! *****************************   


      do j = 1, mlat      
        do i = 1, mlon      
          phis(i,j) = phis(i,j) * grav
          if (oro(i,j) .ge. .5d0) then
            oro(i,j) = 1.d0
          else
            oro(i,j) = 0.d0
          endif
C SJL
          if (sgh(i,j) .gt. .0d0) then  
            sgh(i,j) = sqrt(sgh(i,j))
          else
            sgh(i,j) = 0.d0
          endif

        end do
      end do

c ... Set check=.true. to output real*4 data (GrADS format)
c
      if (check) then
        do j = 1, mlat
          do i = 1, mlon
            buf(i,j) = sgh(i,j)
          end do
        end do
        write (32, rec=2) buf
        do j = 1, mlat
          do i = 1, mlon
            buf(i,j) = oro(i,j)
          end do
        end do
        write (32, rec=3) buf
      endif
c
! *****************************   
      call sm2(mlon, mlat, sgh,  msgh, 0.25D0)
! *****************************   
c
      if (check) then
        do j = 1, mlat
          do i = 1, mlon
            buf(i,j) = phis(i,j)
          end do
        end do
        write (32, rec=4) buf

        fmax =  vmax(buf,fmin,mlon*mlat)
	write(6,*) 'hmax=', fmax / grav
	write(6,*) 'hmin=', fmin / grav

        do j = 1, mlat
          do i = 1, mlon
            buf(i,j) = sgh(i,j)
          end do
        end do
        write (32, rec=5) buf
        close (32)
        fmax =  vmax(buf,fmin,mlon*mlat)
	write(6,*) 'STDmax=', fmax
	write(6,*) 'STDmin=', fmin
      endif

! Safety check
      do j = 1, mlat      
        do i = 1, mlon      
            sgh(i,j) = max(0., sgh(i,j))
        end do
      end do
c ... Open output file (unformatter sequential format)
c
      open (33, file=oflnm, form='unformatted', status='unknown')
      write (33) phis
      write (33) sgh
c     write (33) oro
      close (33)

      open (34, file='oro.init', form='unformatted', status='unknown')
      write (34) oro                          ! sjl
      close (34)

      return
      end

      subroutine get_phis(nlon, nlat, xn, yn, ftopo, htopo, mlon, mlat,
     .                    xm, ym, undef, phis, oro)

      implicit         none

      integer          nlon
      integer          nlat
      integer          mlon
      integer          mlat
      integer          i
      integer          j
      integer          i1
      integer          i2
      integer          icn
      integer          icp
      integer          j1
      integer          j2
      integer          ii
      integer          jj
      integer          count

      real             ftopo(nlon,nlat)
      real             htopo(nlon,nlat)
      real             xn(nlon)
      real             yn(nlat)

      double precision undef
      double precision phis(mlon,mlat)
      double precision oro(mlon,mlat)
      double precision xm(mlon)
      double precision ym(mlat)
      double precision dxm
      double precision dym
      double precision dxn
      double precision dyn
      double precision x1
      double precision x2
      double precision y1
      double precision y2
      double precision fsum
      double precision hsum
      double precision ssum


      dxn = 1.d0 / 6.d0
      dyn = 1.d0 / 6.d0
      dxm = 360.d0 / dble(mlon)
      dym = 180.d0 / dble(mlat - 1)

      do i = 1, nlon
        xn(i) = real(i - 1) * dxn
      end do
      do j = 1, nlat
        yn(j) = -90.d0 + .5d0 * dyn +
     .          real(j - 1) * dyn
      end do

      do i = 1, mlon
        xm(i) = real(i - 1) * dxm
      end do
      do j = 1, mlat
        ym(j) = -90.d0 + real(j - 1) * dym
      end do

      do j = 2, mlat - 1
        y1 = ym(j) - .5d0 * dym
        y2 = y1 + dym
        j1 = int((y1 + 90.d0) / dyn) + 1
        j2 = int((y2 + 90.d0) / dyn)
c       print '(6(1x,f12.4))', yn(j1)-dyn,y1,yn(j1),
c    .                         yn(j2), y2, yn(j2)+dyn
        do i = 1, mlon
          x1 = xm(i) - .5d0 * dxm
          x2 = x1 + dxm
          icn = mod(int((xm(i) + 360.d0) / dxn) - 1, nlon) + 1
          icp = mod(icn + 1, nlon)
          i1 = mod(int((x1 + 360.d0) / dxn), nlon) + 2
          i2 = int(x2 / dxn) + 1
c         print '(8(1x,f12.4))', xn(i1)-dxn,x1,xn(i1),xn(icn),
c    .                          xn(icp),xn(i2),x2,xn(i2)+dxn
c         if (xn(i1) .le. x1 .or. xn(i2) .ge. x2) then
c           print *, ' === x1,xn(i1),xn(i2),x2=', x1,xn(i1),xn(i2),x2
c         endif

          fsum = 0.
          hsum = 0.
          count = 0
          do jj = j1, j2
            do ii = i1, icn
              count = count + 1
              hsum = hsum + htopo(ii,jj)
              fsum = fsum + ftopo(ii,jj)
            end do
            do ii = icp, i2
              count = count + 1
              hsum = hsum + htopo(ii,jj)
              fsum = fsum + ftopo(ii,jj)
            end do
          end do

          if (count .ne. 0) then
            phis(i,j) = hsum / dble(count)
            oro(i,j)  = fsum / dble(count)
          else
            phis(i,j) = undef
            oro(i,j)  = undef
          endif
        end do
      end do
c
c ... South pole
c
      j = 1
      j1 = 1
      y2 = ym(1) + .5d0 * dym
      j2 = int((y2 + 90.d0) / dyn)
      do i = 1, mlon
        x1 = xm(i) - .5d0 * dxm
        x2 = x1 + dxm
        icn = mod(int((xm(i) + 360.d0) / dxn) - 1, nlon) + 1
        icp = mod(icn + 1, nlon)
        i1 = mod(int((x1 + 360.d0) / dxn), nlon) + 2
        i2 = int(x2 / dxn) + 1

        fsum = 0.
        hsum = 0.
        count = 0
        do jj = j1, j2
          do ii = i1, icn
            count = count + 1
            hsum = hsum + htopo(ii,jj)
            fsum = fsum + ftopo(ii,jj)
          end do
          do ii = icp, i2
            count = count + 1
            hsum = hsum + htopo(ii,jj)
            fsum = fsum + ftopo(ii,jj)
          end do
        end do

        if (count .ne. 0) then
          phis(i,j) = hsum / dble(count)
          oro(i,j)  = fsum / dble(count)
        else
          phis(i,j) = undef
          oro(i,j)  = undef
        endif
      end do
c
c ... North pole
c
      j = mlat
      y1 = ym(mlat) - .5d0 * dym
      j1 = int((y1 + 90.d0) / dyn)
      j2 = nlat
      do i = 1, mlon
        x1 = xm(i) - .5d0 * dxm
        x2 = x1 + dxm
        icn = mod(int((xm(i) + 360.d0) / dxn) - 1, nlon) + 1
        icp = mod(icn + 1, nlon)
        i1 = mod(int((x1 + 360.d0) / dxn), nlon) + 2
        i2 = int(x2 / dxn) + 1

        fsum = 0.
        hsum = 0.
        count = 0
        do jj = j1, j2
          do ii = i1, icn
            count = count + 1
            hsum = hsum + htopo(ii,jj)
            fsum = fsum + ftopo(ii,jj)
          end do
          do ii = icp, i2
            count = count + 1
            hsum = hsum + htopo(ii,jj)
            fsum = fsum + ftopo(ii,jj)
          end do
        end do

        if (count .ne. 0) then
          phis(i,j) = hsum / dble(count)
          oro(i,j)  = fsum / dble(count)
        else
          phis(i,j) = undef
          oro(i,j)  = undef
        endif
      end do

      return
      end 

      subroutine get_sgh (nlon, nlat, xn, yn, ftopo, htopo, mlon, mlat,
     .                   xm, ym, undef, phis, sgh)

      implicit         none

      integer          nlon
      integer          nlat
      integer          mlon
      integer          mlat
      integer          i
      integer          j
      integer          i1
      integer          i2
      integer          icn
      integer          icp
      integer          j1
      integer          j2
      integer          ii
      integer          jj
      integer          count

      real             ftopo(nlon,nlat)
      real             htopo(nlon,nlat)
      real             xn(nlon)
      real             yn(nlat)

      double precision undef
C SJL
C Phis must be input !!!
      double precision phis(mlon,mlat)

      double precision sgh(mlon,mlat)
      double precision xm(mlon)
      double precision ym(mlat)
      double precision dxm
      double precision dym
      double precision dxn
      double precision dyn
      double precision x1
      double precision x2
      double precision y1
      double precision y2
      double precision fsum
      double precision hsum
      double precision ssum


      dxn = 1.d0 / 6.d0
      dyn = 1.d0 / 6.d0
      dxm = 360.d0 / dble(mlon)
      dym = 180.d0 / dble(mlat - 1)

      do i = 1, nlon
        xn(i) = real(i - 1) * dxn
      end do
      do j = 1, nlat
        yn(j) = -90.d0 + .5d0 * dyn +
     .          real(j - 1) * dyn
      end do

      do i = 1, mlon
        xm(i) = real(i - 1) * dxm
      end do
      do j = 1, mlat
        ym(j) = -90.d0 + real(j - 1) * dym
      end do

      do j = 2, mlat - 1
        y1 = ym(j) - .5d0 * dym
        y2 = y1 + dym
        j1 = int((y1 + 90.d0) / dyn) + 1
        j2 = int((y2 + 90.d0) / dyn)
c       print '(6(1x,f12.4))', yn(j1)-dyn,y1,yn(j1),
c    .                         yn(j2), y2, yn(j2)+dyn
        do i = 1, mlon
          x1 = xm(i) - .5d0 * dxm
          x2 = x1 + dxm
          icn = mod(int((xm(i) + 360.d0) / dxn) - 1, nlon) + 1
          icp = mod(icn + 1, nlon)
          i1 = mod(int((x1 + 360.d0) / dxn), nlon) + 2
          i2 = int(x2 / dxn) + 1
c         print '(8(1x,f12.4))', xn(i1)-dxn,x1,xn(i1),xn(icn),
c    .                          xn(icp),xn(i2),x2,xn(i2)+dxn
c         if (xn(i1) .le. x1 .or. xn(i2) .ge. x2) then
c           print *, ' === x1,xn(i1),xn(i2),x2=', x1,xn(i1),xn(i2),x2
c         endif

          fsum = 0.
          hsum = 0.
          count = 0
          do jj = j1, j2
            do ii = i1, icn
              count = count + 1
              hsum = hsum + htopo(ii,jj)
              fsum = fsum + ftopo(ii,jj)
            end do
            do ii = icp, i2
              count = count + 1
              hsum = hsum + htopo(ii,jj)
              fsum = fsum + ftopo(ii,jj)
            end do
          end do

          if (count .ne. 0) then
c           phis(i,j) = hsum / dble(count)

            ssum = 0.d0
            do jj = j1, j2
              do ii = i1, icn
                ssum = ssum +
     .                 (htopo(ii,jj) - phis(i,j)) ** 2
              end do
              do ii = icp, i2
                ssum = ssum +
     .                 (htopo(ii,jj) - phis(i,j)) ** 2
              end do
            end do
            sgh(i,j) = ssum / dble(count)
          else
c           phis(i,j) = undef
            sgh(i,j)  = undef
          endif
        end do
      end do
c
c ... South pole
c
      j = 1
      j1 = 1
      y2 = ym(1) + .5d0 * dym
      j2 = int((y2 + 90.d0) / dyn)
      do i = 1, mlon
        x1 = xm(i) - .5d0 * dxm
        x2 = x1 + dxm
        icn = mod(int((xm(i) + 360.d0) / dxn) - 1, nlon) + 1
        icp = mod(icn + 1, nlon)
        i1 = mod(int((x1 + 360.d0) / dxn), nlon) + 2
        i2 = int(x2 / dxn) + 1

        fsum = 0.
        hsum = 0.
        count = 0
        do jj = j1, j2
          do ii = i1, icn
            count = count + 1
            hsum = hsum + htopo(ii,jj)
            fsum = fsum + ftopo(ii,jj)
          end do
          do ii = icp, i2
            count = count + 1
            hsum = hsum + htopo(ii,jj)
            fsum = fsum + ftopo(ii,jj)
          end do
        end do

        if (count .ne. 0) then
c         phis(i,j) = hsum / dble(count)
          ssum = 0.d0
          do jj = j1, j2
            do ii = i1, icn
              ssum = ssum +
     .               (htopo(ii,jj) - phis(i,j)) ** 2
            end do
            do ii = icp, i2
              ssum = ssum +
     .               (htopo(ii,jj) - phis(i,j)) ** 2
            end do
          end do
          sgh(i,j) = ssum / dble(count)
        else
c         phis(i,j) = undef
          sgh(i,j)  = undef
        endif
      end do
c
c ... North pole
c
      j = mlat
      y1 = ym(mlat) - .5d0 * dym
      j1 = int((y1 + 90.d0) / dyn)
      j2 = nlat
      do i = 1, mlon
        x1 = xm(i) - .5d0 * dxm
        x2 = x1 + dxm
        icn = mod(int((xm(i) + 360.d0) / dxn) - 1, nlon) + 1
        icp = mod(icn + 1, nlon)
        i1 = mod(int((x1 + 360.d0) / dxn), nlon) + 2
        i2 = int(x2 / dxn) + 1

        fsum = 0.
        hsum = 0.
        count = 0
        do jj = j1, j2
          do ii = i1, icn
            count = count + 1
            hsum = hsum + htopo(ii,jj)
            fsum = fsum + ftopo(ii,jj)
          end do
          do ii = icp, i2
            count = count + 1
            hsum = hsum + htopo(ii,jj)
            fsum = fsum + ftopo(ii,jj)
          end do
        end do

        if (count .ne. 0) then
c         phis(i,j) = hsum / dble(count)
          ssum = 0.d0
          do jj = j1, j2
            do ii = i1, icn
              ssum = ssum +
     .               (htopo(ii,jj) - phis(i,j)) ** 2
            end do
            do ii = icp, i2
              ssum = ssum +
     .               (htopo(ii,jj) - phis(i,j)) ** 2
            end do
          end do
          sgh(i,j) = ssum / dble(count)
        else
c         phis(i,j) = undef
          sgh(i,j)  = undef
        endif
      end do

      return
      end 


      subroutine handle_err(status)

      implicit         none

#     include          <netcdf.inc>

      integer          status

      if (status .ne. nf_noerr) then
        print *, nf_strerror(status)
        stop 'Stopped'
      endif

      return
      end

      real function vmax(a,pmin,n)
      implicit none
      integer n, i
      real pmin, pmax
      real a(*)
      pmax =a(n)
      pmin =a(n)
      do 10 i=1,n-1
      pmax = max(pmax,a(i))
      pmin = min(pmin,a(i))
10    continue
      vmax = pmax
      return
      end

! **********
! Module sm2
! **********

      subroutine sm2(im, jm, ht, itmax, c)
!
! Del-2 diffusion on the sphere
!
      implicit none

! Input:
      integer im               ! e-w dimension (eg, 144 for 2.5 deg resolution)
      integer jm               ! n-s doemsnion (eg, 91 for 2 deg resolution)
      integer itmax            ! iteration count
      real*8  c                ! filter coefficient

! Input/Output
      real*8 ht(im,jm)         ! array to be filtered

! Local
      real*8 dg(im,jm)   ! del2 of h
      real*8 cose(jm), cosp(jm), sinp(jm), sine(jm)
      real*8 dl
      real*8 dp
      real*8 fmin, fmax
      integer jm1
      integer mnk, mxk
      integer ndeg
      integer it, i, j
      real*8 s1, s2

      jm1 = jm-1

      call setrig(im, jm, dp, DL, cosp, cose, sinp, sine)

      call pmnx(ht, im, jm, fmin, fmax, mnk, mxk)
      write(6,*) 'hmax=', fmax,' at j= ',mxk
      write(6,*) 'hmin=', fmin,' at j= ',mnk

      ndeg = 60                    ! starting latitude for the monotonicity
                                   ! preserving polar filter

      call pmnx(ht,im,jm,fmin,fmax,mnk,mxk)
      write(6,*) 'hmax=', fmax,' at j= ',mxk
      write(6,*) 'hmin=', fmin,' at j= ',mnk

! Apply Monotonicity preserving polar filter
      call plft2d(im, jm, ht, 2, jm1, ndeg)
      call avgp2(ht, sine, im, jm)

      do it=1,itmax
	call del2(ht, im, jm, dg, cosp, cose, sine, DL, dp, ndeg)
        call plft2d(im, jm, dg, 2, jm1, ndeg)

	do j=1,jm
	   do i=1,im
             ht(i,j) = ht(i,j) + c*dg(i,j)
           enddo
        enddo
      enddo

! Final polar filter
      call plft2d(im, jm, ht, 2, jm1, ndeg)

      return
      end

      subroutine del2(h, im, jm, dg, cosp, cose, sine, dL, dp, ndeg)
      implicit none

! AE = 1  (unit radius)
! Input:
        integer im
        integer jm
        integer ndeg
! Input-output

	real*8 h(im,jm)
	real*8 dg(im,jm)              ! del2 of h
	real*8 cose(jm),cosp(jm)
	real*8 sine(jm)
	real*8 PI, ycrit, coszc, CD
	real*8 DL, dp

! Local
	real*8 fx(im,jm)   ! e-w fluxes
	real*8 fy(im,jm)   ! n-s fluxes
        integer i, j
 
	call grad(h, im, jm, fx, fy, cosp, dl, dp)

        PI = 4. * ATAN(1.)
        ycrit = float(ndeg)*PI/180.
        coszc = cos(ycrit)

 	CD = 0.25*DL*DP*coszc**2
!	CD = 0.25*DL*DP*cosp(2)**2

        do j=2,jm-1
           do i=1,im
              fx(i,j) = fx(i,j) * CD
           enddo
        enddo

        do j=2,jm
	   do i=1,im
	      fy(i,j) = fy(i,j) * CD
           enddo
        enddo

	call divg(im,jm,fx,fy,DG,cosp,cose,sine, dl, dp)

	return
	end

      subroutine divg(im, jm, fx, fy, dg, cosp, cose, sine, dl, dp)
      implicit none

      integer im
      integer jm
      real*8 fx(im,jm)   ! e-w fluxes
      real*8 fy(im,jm)   ! n-s fluxes
      real*8 DG(im,jm)   ! del2 of h
      real*8 wk(im,jm)
      real*8 cosp(jm),  cose(jm), sine(jm)
      real*8 rdx
      real*8 dl, dp, CDP, sum1, sum2
      integer i,j

	do j=2,jm-1

	rdx = 1./ (cosp(j)*DL)

           do i=1,im-1
              DG(i,j) = (fx(i+1,j) - fx(i,j)) * rdx
           enddo
              DG(im,j) = (fx(1,j) - fx(im,j)) * rdx
        enddo

	do j=2,jm
           do i=1,im
              wk(i,j) = fy(i,j) * cose(j)
           enddo
        enddo

	do j=2,jm-1
           CDP = 1./ (DP*cosp(j))
           do i=1,im
              DG(i,j) = DG(i,j) + (wk(i,j+1) - wk(i,j)) * CDP
           enddo
	enddo

C Poles;

	sum1 = wk(im, 2)
	sum2 = wk(im,jm)

	do i=1,im-1
           sum1 = sum1 + wk(i, 2)
           sum2 = sum2 + wk(i,jm)
	enddo

	sum1 =  sum1 / ( float(im)*(1.+sine(2)) )
	sum2 = -sum2 / ( float(im)*(1.+sine(2)) )

	do i=1,im
           DG(i, 1) =  sum1
           DG(i,jm) =  sum2
	enddo

	return
	end

      subroutine grad(h, im, jm, fx, fy, cosp, DL, DP)
      implicit none
      integer im
      integer jm
      real*8 h(im,jm)
      real*8 fx(im,jm)   ! e-w fluxes
      real*8 fy(im,jm)   ! n-s fluxes
      real*8 cosp(jm)
      real*8 RDP, DL, DP, rdx
      integer i, j

	RDP = 1./ DP

      do j=2,jm
         do i=1,im
            fy(i,j) = (h(i,j) - h(i,j-1)) * RDP
         enddo
      enddo

      do j=2,jm-1

	rdx = 1./ (cosp(j)*DL)
        fx(1,j) = (h(1,j) - h(im,j)) * rdx
        do i=2,im
           fx(i,j) = (h(i,j) - h(i-1,j)) * rdx
        enddo
      enddo

      return
      end

      subroutine avgp2(p, sine, im, jm)
      implicit none
      integer im, jm
      real*8 p(im,jm)
      real*8 sine(jm)
      real*8 sum1, sum2
      real*8 sum3, sum4
      real*8 rim
      integer i
      integer j
      integer jm1

      jm1 = jm-1
      rim = 1./ float(im)

      call sump2(p(1,1),p(1,jm),IM,sum1,sum2)
      sum1 = sum1*(1.+sine(2))
      sum2 = sum2*(1.+sine(2))

      call sump2(p(1,2),p(1,jm1),IM,sum3,sum4)
      sum1 = rim * ( sum1 + sum3*(sine(3)-sine(2)) ) / (1.+sine(3))
      sum2 = rim * ( sum2 + sum4*(sine(3)-sine(2)) ) / (1.+sine(3))

      do i=1,im
      P(i,  1) = sum1
      P(i,  2) = sum1
      P(i,jm1) = sum2
      P(i, jm) = sum2
      enddo
      return
      end

      subroutine sump2(p1,p2,im,s1,s2)
      implicit none
      integer im,i
      real*8 s1,s2
      real*8 p1(*),p2(*)
 
         s1 =  p1(im)
         s2 =  p2(im)
 
      do i=1,im-1
         s1 =  s1 + p1(i)
         s2 =  s2 + p2(i)
      enddo
      return
      end

      subroutine plft2d(im, jm, p, JS, JN, ndeg)
!
! This is a weak LOCAL polar filter.
! Developer: Shian-Jiann Lin

      implicit none

      integer im
      integer jm
      integer js, jn, ndeg
      real*8 p(im,jm)

      integer i, j, n, ideg, jj, jc
      real*8 cosp(jm),cose(jm)
      real*8 a(0:im/2+1)

      real*8 sine(jm),sinp(jm)
      real*8, allocatable, save :: se(:), sc(:)

      real*8 pi, dp, dl, e0, ycrit, coszc, smax, rn, rn2, esl, tmp

      data IDEG /0/

      if(IDEG .ne. ndeg) then
	IDEG = ndeg
C  (e0 = 2.6)
      e0 = 0.5 * sqrt(27.)
      PI = 4. * ATAN(1.)

      allocate( sc(jm), se(jm))

      call setrig(im, jm, DP, DL, cosp, cose, sinp, sine)

      ycrit = IDEG*PI/180.
      coszc = cos(ycrit)

      smax = (jm-1)/2
      write(6,*) 'Critical latitude in local pft = ',ndeg

      a(0) = 1.
      do n=1,im/2+1
      rn = n
      rn2 = 2*n
      a(n) = sqrt(rn2+1.) * ((rn2+1.)/rn2)**rn
      enddo

      do j=2,jm-1
      sc(j) = coszc / cosp(j)

      IF(sc(j) .gt.1. .and. sc(j) .le. 1.5 ) THEN
         esl = 1./ sc(j)
         sc(j) =  1. +  (1.-esl) / (1.+esl)
      ELSEIF(sc(j) .gt.1.5 .and. sc(j) .le. e0 ) THEN
         esl = 1./ sc(j)
         sc(j) =  1. + 2./ (27.*esl**2 - 2.)
      ELSEIF(sc(j).gt. e0) THEN
C Search
      do jj=1,im/2
      if(sc(j).le. a(jj)) then
      jc = jj
      write(*,*) 'jc=', jc
      goto 111
      endif
      enddo

      jc = im/2 + 1
111   continue
      tmp = ((sc(j) - a(jc-1))/(a(jc) - a(jc-1)))**0.25
      sc(j) =  jc + min(1., tmp)
c     sc(j) =  min(smax,sc(j))
      ENDIF
      enddo
C ====================================================
      do j=2,jm
      se(j) = coszc / cose(j)
      IF(se(j) .gt.1. .and. se(j) .le. 1.5 ) THEN
         esl = 1./ se(j)
         se(j) =  1. + (1.-esl) / (1.+esl)
      ELSEIF(se(j) .gt.1.5 .and. se(j) .le. e0 ) THEN
         esl = 1./ se(j)
         se(j) =  1. + 2./ (27.*esl**2 - 2.)
      ELSEIF(se(j).gt. e0) THEN
C Search
      do jj=1,im/2
      if(se(j) .le. a(jj)) then
      jc = jj
      goto 222
      endif
      enddo

      jc = im/2 + 1
222   continue
      tmp = ((se(j) - a(jc-1))/(a(jc) - a(jc-1)))**0.25
      se(j) =  jc + min(1., tmp)
c     se(j) =  min(smax,se(j))
      ENDIF
      enddo

      do i=1,im
        se( 2) = sc(2)
        se(jm) = sc(jm-1)
      enddo

      do j=2,jm-1
         write(*,*) j,sc(j)
      enddo
      ENDIF

      if(JN .eq. (jm-1) ) then
C Cell-centered variables
         call lpft(im, jm, p, 2, jm-1, Sc)
      else
C Cell-edge variables
         call lpft(im, jm, p, 2, jm, Se)
      endif
      return
      end


      subroutine lpft(im, jm, p, j1, j2, s)
      implicit none

      integer im, jm, j1, j2
      real*8  p(im,jm)
      real*8  s(jm)

! Local
      integer i, j, n, nt

      real*8 ptmp(0:im+1)
      real*8 q(0:im+1)
      real*8 frac, rsc, bt

      do 2500 j=j1,j2
      if(s(j) .gt. 1.02) then

        NT  = INT(S(j))
        frac = S(j) - NT
        NT = NT-1

        rsc = 1. / (1.+frac)
        bt = 0.5 * frac

        do i=1,im
           ptmp(i) = p(i,j)
        enddo

        ptmp(0)    = p(im,j)
        ptmp(im+1) = p(1 ,j)

        if( NT .lt. 1 ) then
          do i=1,im
             p(i,j) = rsc * (ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)))
          enddo
        else
          do i=1,im
             q(i) = rsc * (ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)))
          enddo

           do 500 N=1,NT
              q(0) = q(im)
              do i=1,im
                 ptmp(i) = q(i) + q(i-1)
              enddo
                 ptmp(im+1) = ptmp(1)

             if ( n .eq. nt ) then
               do i=1,im
                  p(i,j) = 0.25*(ptmp(i) + ptmp(i+1))
               enddo
             else
               do i=1,im
                  q(i) = 0.25*(ptmp(i) + ptmp(i+1))
               enddo
             endif
500        continue
        endif
      endif
2500  continue

      return
      end

      subroutine pmnx(a,nx,ny,fmin,fmax,mnk,mxk)
      implicit none
      integer nx
      integer ny
      integer mnk
      integer mxk
      real*8 a(nx,*)
      real*8 fmax, fmin, temp
      integer i,j

      fmax = a(1,1)
      fmin = a(1,1)
      mnk = 1
      mxk = 1

      do j=1,ny
        do i=1,nx
           temp = a(i,j)
           if(temp.gt.fmax) then
              fmax = temp
              mxk = j
           elseif(temp .lt. fmin) then
              fmin = temp
              mnk = j
           endif
        enddo
      enddo

      return
      end

      subroutine setrig(im, jm, dp, dl, cosp, cose, sinp, sine)

      implicit none

      integer im, jm
      integer j, jm1
      real*8 sine(jm),cosp(jm),sinp(jm),cose(jm)
      real*8 dp, dl
      double precision pi, ph5

      jm1 = jm - 1
      pi  = 4.d0 * datan(1.d0)
      dl  = (pi+pi)/dble(im)
      dp  = pi/dble(jm1)

      do 10 j=2,jm
         ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(jm1))
10    sine(j) = dsin(ph5)

      cosp( 1) =  0.
      cosp(jm) =  0.

      do 80 j=2,jm1
80    cosp(j) = (sine(j+1)-sine(j)) / dp

c Define cosine at edges..

      do 90 j=2,jm
90    cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      cose(1) = cose(2)

      sinp( 1) = -1.
      sinp(jm) =  1.

      do 100 j=2,jm1
100   sinp(j) = 0.5 * (sine(j) + sine(j+1))

      return
      end
