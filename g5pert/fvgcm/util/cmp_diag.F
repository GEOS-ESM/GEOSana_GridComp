      program cmp_diag
c
c ... To compile: f77 -bytereclen -mips4 -n32 -O -o cmp_diag cmp_diag.F
c                 f90 -cpp -mips4 -n32 -O -o cmp_diag cmp_diag.F
c
      implicit none

      integer imr, jnp, nl
      parameter  (imr = 144)
      parameter  (jnp =  91)
      parameter  (nl  =  55)

      integer n2d, n3d
      parameter  (n2d = 36)
      parameter  (n3d = 24)
      integer nslice
      parameter  (nslice = n2d + n3d * nl )

      integer inx1(n2d+n3d)
      integer inx2(n2d+n3d)
      integer count(n2d+n3d)
      integer nt                    ! time samples
      parameter (nt = 4)
      integer n1, n2
      integer ir1, ir2
      integer i, j, k, n, iv

      logical pick1(n2d+n3d)
      logical pick2(n2d+n3d)

      real buf1(imr,jnp)
      real buf2(imr,jnp)
      
      ir1 =11
      ir2 =12

      do n = 1, n2d + n3d
        pick1(n) = .true.
        pick2(n) = .true.
        inx1(n) = 0
        inx2(n) = 0
        count(n) = 0
      end do

c     pick1(27) = .false.   ! SURFP
c     pick2(27) = .false.   ! SURFP

      n1 = 0
      n2 = 0

      do n = 1, n2d
        if (pick1(n)) then
         inx1(n) = n1 + 1
         n1 = n1 + 1
        endif
        if (pick2(n)) then
         inx2(n) = n2 + 1
         n2 = n2 + 1
        endif
      end do
      do n = n2d+1, n2d+n3d
        if (pick1(n)) then
          inx1(n) = n1 + 1
          n1 = n1 + nl
        endif
        if (pick2(n)) then
          inx2(n) = n2 + 1
          n2 = n2 + nl
        endif
      end do

      do n = 1, n2d + n3d
        print *, n, inx1(n), inx2(n)
      end do
      print *,  '    ', n1, n2

      open (ir1, file='diag1', status='unknown', form='unformatted',
     .      access='direct', recl=imr*jnp*4)
      open (ir2, file='diag2', status='unknown', form='unformatted',
     .      access='direct', recl=imr*jnp*4)

      do n = 1, nt
        print *, '******** Time = ', n, ' *******'
        do iv = 1, n2d
          if (pick1(iv) .and. pick2(iv)) then
            read (ir1, rec=inx1(iv)) buf1
            read (ir2, rec=inx2(iv)) buf2
            do j = 1, jnp
              do i = 1, imr
                if (buf1(i,j) .ne. buf2(i,j)) then
                  count(iv) = count(iv) + 1
                endif
              end do
            end do
          endif
          print *, iv, count(iv)
        end do

        do iv = n2d+1, n2d+n3d
          if (pick1(iv) .and. pick2(iv)) then
            read (ir1, rec=inx1(iv)) buf1
            read (ir2, rec=inx2(iv)) buf2
            do k = 1, nl
              do j = 1, jnp
                do i = 1, imr
                  if (buf1(i,j) .ne. buf2(i,j)) then
                    count(iv) = count(iv) + 1
                  endif
                end do
              end do
            end do
          endif
          print *, iv, count(iv)
        end do
      end do

      close (ir1)
      close (ir2)

      stop
      end
