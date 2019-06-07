      program cmp_drst
      implicit none
      
      integer    imr, jmr, nlayr, nc
      parameter  (imr   = 144)
      parameter  (jmr   =  90)
      parameter  (nlayr =  55)
      parameter  (nc    =   1)

      integer  jnp, nl
      parameter(jnp= jmr+1, nl=nlayr)

      real ps1(imr,jnp)
      real delp1(imr,jnp,nl)
      real u1(imr,jnp,nl)
      real v1(imr,jnp,nl)
      real pt1(imr,jnp,nl)
      real q1(imr,jnp,nl,nc)

      real ps2(imr,jnp)
      real delp2(imr,jnp,nl)
      real u2(imr,jnp,nl)
      real v2(imr,jnp,nl)
      real pt2(imr,jnp,nl)
      real q2(imr,jnp,nl,nc)
      real qsum

      real dif, dmax, dmin

      integer nstep, nymd, nhms
      integer ir1, ir2
      integer i, j, k
      integer imax, jmax, kmax
      integer count

      ir1 = 31
      open (ir1, file='d_rst1', form='unformatted', status='unknown')
      read (ir1) nstep, nymd, nhms
      read (ir1) ps1, delp1, u1, v1, pt1
      if(nc .ne. 0) read(ir1) q1

      ir2 = 32
      open (ir2, file='d_rst2', form='unformatted', status='unknown')
      read (ir2) nstep, nymd, nhms
      read (ir2) ps2, delp2, u2, v2, pt2
      if(nc .ne. 0) read(ir2) q2

      write(*,*) ' '
      write(*,*) '<<< Max and min >>>'
      call minmax ('ps1', imr, jnp, 1, ps1)
      call minmax ('ps2', imr, jnp, 1, ps2)
      call minmax ('delp1', imr, jnp, nl, delp1)
      call minmax ('delp2', imr, jnp, nl, delp2)
      call minmax ('u1', imr, jnp, nl, u1)
      call minmax ('u2', imr, jnp, nl, u2)
      call minmax ('v1', imr, jnp, nl, v1)
      call minmax ('v2', imr, jnp, nl, v2)
      call minmax ('pt1', imr, jnp, nl, pt1)
      call minmax ('pt2', imr, jnp, nl, pt2)
      call minmax ('q1', imr, jnp, nl, q1(1,1,1,1))
      call minmax ('q2', imr, jnp, nl, q2(1,1,1,1))


      count = 0
      dmax = 1.E-25
              imax = 1
              jmax = 1
              kmax = 1

      do j = 1, jnp
        do i = 1, imr
          dif = ps2(i,j) - ps1(i,j)
          if (dif .ne. 0.) then
            count = count + 1
            if (abs(dif/ps1(i,j)) .gt. dmax) then
              dmax = abs(dif/ps1(i,j))
              imax = i
              jmax = j
            endif
          endif
        end do
      end do

      dmax = dmax*100.

      write(*,*) ' '
      write(*,*) '<<< Relative Diff (%) >>> '
      print *, '  ps:', count, dmax
      print *, 'psmax:', imax, jmax,
     .                  ps1(imax,jmax), ps2(imax,jmax)

      count = 0
      dmax = 1.E-25
              imax = 1
              jmax = 1

      do k = 1, nl
        do j = 1, jnp
          do i = 1, imr
            dif = delp2(i,j,k) - delp1(i,j,k)
            if (dif .ne. 0.) then
              count = count + 1
              if (abs(dif/delp1(i,j,k)) .gt. dmax) then
                dmax = abs(dif/delp1(i,j,k))
                imax = i
                jmax = j
                kmax = k
              endif
            endif
          end do
        end do
      end do
      write(*,*) ' '
      print *, 'delp:', count, dmax
      print *, 'delpmax:', imax, jmax, kmax,
     .                  delp1(imax,jmax,kmax), delp2(imax,jmax,kmax)

      count = 0
      imax = 0
      jmax = 0
      kmax = 0
      dmax = 1.e-25
      do k = 1, nl
        do j = 1, jnp
          do i = 1, imr
            dif = u2(i,j,k) - u1(i,j,k)
            if (dif .ne. 0.) then
              count = count + 1
              if (abs(dif/u1(i,j,k)) .gt. dmax) then
                dmax = abs(dif/u1(i,j,k))
                imax = i
                jmax = j
                kmax = k
              endif
            endif
          end do
        end do
      end do
      write(*,*) ' '
      print *, '   u:', count, dmax
      print *, 'umax:', imax, jmax, kmax,
     .                  u1(imax,jmax,kmax), u2(imax,jmax,kmax)

      count = 0
      dmax = 1.e-25
      imax = 0
      jmax = 0
      kmax = 0
      do k = 1, nl
        do j = 1, jnp
          do i = 1, imr
            dif = v2(i,j,k) - v1(i,j,k)
            if (dif .ne. 0.) then
              count = count + 1
              if (abs(dif/v1(i,j,k)) .gt. dmax) then
                dmax = abs(dif/v1(i,j,k))
                imax = i
                jmax = j
                kmax = k
              endif
            endif
          end do
        end do
      end do
      write(*,*) ' '
      print *, '   v:', count, dmax
      print *, 'vmax:', imax, jmax, kmax,
     .                  v1(imax,jmax,kmax), v2(imax,jmax,kmax)

      count = 0
      dmax = 1.e-25
      do k = 1, nl
        do j = 1, jnp
          do i = 1, imr
            dif = pt2(i,j,k) - pt1(i,j,k)
            if (dif .ne. 0.) then
              count = count + 1
              if (abs(dif/pt1(i,j,k)) .gt. dmax) then
                dmax = abs(dif/pt1(i,j,k))
                imax = i
                jmax = j
                kmax = k
              endif
            endif
          end do
        end do
      end do
      write(*,*) ' '
      print *, '  pt:', count, dmax
      print *, 'ptmax:', imax, jmax, kmax,
     .                  pt1(imax,jmax,kmax), pt2(imax,jmax,kmax)

      count = 0
      dmax = 0
      imax = 1
      jmax = 1
      kmax = 1
      do k = 1, nl
        do j = 1, jnp
          do i = 1, imr
            dif = q2(i,j,k,1) - q1(i,j,k,1)
            if (dif .ne. 0.) then
              count = count + 1
              if (abs(dif) .gt. dmax) then
                dmax = abs(dif)
                imax = i
                jmax = j
                kmax = k
              endif
            endif
          end do
        end do
      end do
      write(*,*) ' '
      print *, '   q:', count, dmax
      print *, 'qmax:', imax, jmax, kmax,
     .                  q1(imax,jmax,kmax,1), q2(imax,jmax,kmax,1)

      close (ir1)
      close (ir2)

      stop
      end

      subroutine minmax (name, m, n, l, f)
      
      implicit         none

      character*(*)      name

      integer          m, n, l
      integer          i, j, k

      real             f(m,n,l)
      real             fmax
      real             fmin

      fmax = f(1,1,1)
      fmin = f(1,1,1)

      do k = 1, l
        do j = 1, n
          do i = 1, m
            fmax = max(fmax,f(i,j,k))
            fmin = min(fmin,f(i,j,k))
          end do
        end do
      end do

      write(6,*) name, ' max = ', fmax, ' min = ', fmin
c     print '(a8, 2(1x,a7,e16.8))',
c    .      name, ' max = ', fmax, ' min = ', fmin

      return
      end
