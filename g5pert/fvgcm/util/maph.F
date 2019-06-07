      program mapxy
 
! Program to map fvgcm dynamical core restart file: d_rst from
! one horizontal to another
! S.-J. Lin: April 1, 2000

      implicit none

      integer im, jm                    ! original dimension
      integer in, jn                    ! target dimension
      integer km
      integer nc
      integer hold
      integer hnew

      integer iuic, iout

      write(6,*) '**************************************************'
      write(6,*) 'This program is for the horizontal mapping of fvgcm'
      write(6,*) 'dynamical restart file: d_rst'
      write(6,*) '**************************************************'

      write(*,*) 'Original data resolution?'
      write(*,*) '0: 4x5; 1: 2x2.5; 2: 1x1.25; 3: 0.5x0.625'
      read(*,*) hold

      write(*,*) 'Target data resolution?'
      write(*,*) '0: 4x5; 1: 2x2.5; 2: 1x1.25; 3: 0.5x0.625'
      read(*,*) hnew

! Original resolution
      if( hold .eq. 0 ) then
          im = 72
          jm = 46
      elseif( hold .eq. 1) then
          im = 144
          jm =  91
      elseif( hold .eq. 2) then
          im = 288
          jm = 181
      elseif( hold .eq. 3) then
          im = 576
          jm = 361
      else
          write(*,*) 'No suitable resolution chosen'
          stop
      endif

! New resolution
      if( hnew .eq. 0 ) then
          in = 72
          jn = 46
      elseif( hnew .eq. 1) then
          in = 144
          jn =  91
      elseif( hnew .eq. 2) then
          in = 288
          jn = 181
      elseif( hnew .eq. 3) then
          in = 576
          jn = 361
      else
          write(*,*) 'No suitable resolution chosen'
          stop
      endif

      write(6,*) 'Vertical dimension km=?'
      read(*,*) km

      write(6,*) 'Total number of tracers including water vapor=?'
      read(*,*) nc

      write(6,*) ' '
      write(6,*) 'Input file name is assumed to be d_rst'
      write(6,*) 'Output file will be d_rst_new'

      iuic = 71
      iout = 81

      open (unit=iuic,file='d_rst',form='unformatted',status='old')
      open (unit=iout,file='d_rst_new',form='unformatted',
     &      status='unknown')

      call h_map_drst(iuic, iout, im, jm, km, nc, in, jn)

      end


      subroutine h_map_drst(iuic, iout, im, jm, km, nc, in, jn)
!
! Purpose: perform horizontal mapping of d_rst
!
      implicit none

      integer im, jm, km, nc
      integer in, jn
      integer iuic, iout

! Original
      real    psm(im, jm)     
      real     um(im,jm,km)        !zonal wind on D-grid
      real     vm(im,jm,km)        !meridional wind
      real    ptm(im,jm,km)        !scaled potential temperature
      real     qm(im,jm,km,nc)     !specific humidity & mixing ratios
      real  delpm(im,jm,km)        ! pressure thickness

! New data 
      real    psn(in,jn)     
      real     un(in,jn,km)        !zonal wind on D-grid
      real     vn(in,jn,km)        !meridional wind
      real    ptn(in,jn,km)        !scaled potential temperature
      real     qn(in,jn,km,nc)     !specific humidity & mixing ratios
      real  delpn(in,jn,km)        ! pressure thickness

      integer nymd, nhms, nstep
      real pmax, pmin, mean
      real  gmean

      real dx1, dx2
      real dy1, dy2

      real lon1(im+1)
      real lon2(in+1)

! Pole to pole
      real sin1(jm+1)
      real sin2(jn+1)

      integer i, j
      integer ks
      real pint, ptop
      real ak(km+1), bk(km+1)
      integer ic
      real dl, dp
      real cosp(jm), cose(jm), sine(jm), sinp(jm)
      real sinlon(im), coslon(im), cosl5(im), sinl5(im)
      integer imh, k
      double precision pi, zamda, zam5

      imh = im/2

      pi = 4.d0 * datan(1.d0)
      call setrig(im,jm,dp,dl,cosp,cose,sinp,sine)

        do i=1,imh
           zam5          = (dble(i-1)-0.5d0) * dl
           cosl5(i)      = -dcos(zam5)
           cosl5(i+imh)  = -cosl5(i)
           sinl5(i)      = -dsin(zam5)
           sinl5(i+imh)  = -sinl5(i)
           zamda         = dble(i-1)*dl
           coslon(i)     = -dcos(zamda)
           coslon(i+imh) = -coslon(i)
           sinlon(i)     = -dsin(zamda)
           sinlon(i+imh) = -sinlon(i)
        enddo

      dx1 = 360./float(im)
      dx2 = 360./float(in)

      dy1 = pi/float(jm-1)
      dy2 = pi/float(jn-1)

! get original data with km layers

      call rst_dyn(im, jm, km, 1 , iuic , delpm, um, vm, ptm, nc,
     &             qm   , psm ,  nymd,  nhms, nstep )
 
      call pmaxmin('PS_old', psm, pmin, pmax, im, jm, 0.01 )
      mean = gmean(im,jm,1,jm,psm)
      write(*,*) 'Mean ps_old=',mean

!**************
! Map surface pressure
!**************

      do i=1,im+1
         lon1(i) = dx1 * (-0.5 + (i-1) )
      enddo

      do i=1,in+1
         lon2(i) = dx2 * (-0.5 + (i-1) )
      enddo

         sin1(   1) = -1.
         sin1(jm+1) =  1.

         sin2(   1) = -1.
         sin2(jn+1) =  1.

           do j=2,jm
              sin1(j) = sin( -0.5*pi + dy1*(-0.5+float(j-1)) )
           enddo

           do j=2,jn
              sin2(j) = sin( -0.5*pi + dy2*(-0.5+float(j-1)) )
           enddo

      call map_a2a(1 , im, jm, lon1, sin1, psm,
     &                 in, jn, lon2, sin2, psn, 0, 0)

      call pmaxmin('PS_new', psn, pmin, pmax, in, jn, 0.01 )
      mean = gmean(in,jn,1,jn,psn)
      write(*,*) 'Mean ps_new=',mean

!**************
! Map delp
!**************
      call map_a2a(km , im, jm, lon1, sin1, delpm,
     &                  in, jn, lon2, sin2, delpn, 0, 0)

!**************
! Map pt
!**************
      call map_a2a(km , im, jm, lon1, sin1, ptm,
     &                  in, jn, lon2, sin2, ptn, 0, 0)
      call pmaxmin('PT_new',ptn, pmin, pmax,in*jn, km, 1. )

!**************
! Constituents:
!**************
      do ic=1,nc
          call map_a2a(km, im, jm, lon1, sin1, qm(1,1,1,ic),
     &                     in, jn, lon2, sin2, qn(1,1,1,ic), 0, 0)
          call pmaxmin('Q_new', qn(1,1,1,ic), pmin, pmax,in*jn, km, 1. )
      enddo


! get D-grid v-wind at poles:
      do k=1,km
         call vpol5(um(1,1,k) ,vm(1,1,k), im, jm, coslon, sinlon,
     &              cosl5, sinl5, 1, jm)
      enddo

!**************
! Map u-wind
!**************
! ig = 1                   ! data defined "half dy" from pole
         sin1( 1) = -1.
         sin1(jm) =  1.

         sin2( 1) = -1.
         sin2(jn) =  1.

           do j=2,jm-1
              sin1(j) = sin( -0.5*pi + dy1*float(j-1) )
           enddo

           do j=2,jn-1
              sin2(j) = sin( -0.5*pi + dy2*float(j-1) )
           enddo

! Note: u(i,j=1,k) are NOT defined

      call map_a2a(km, im, jm, lon1, sin1, um,
     &                 in, jn, lon2, sin2, un, 1, 1)
      call pmaxmin('U_new', un, pmin, pmax, in*jn, km, 1. )

!**************
! Map v-wind
!**************
! v loated half-dx west cell center

      do i=1,im+1
         lon1(i) = dx1 * float(i-2)
      enddo

      do i=1,in+1
         lon2(i) = dx2 * float(i-2)
      enddo

         sin1(   1) = -1.
         sin1(jm+1) =  1.

         sin2(   1) = -1.
         sin2(jn+1) =  1.

      do j=2,jm
         sin1(j) = sin( -0.5*pi + dy1*(-0.5+float(j-1)) )
      enddo

      do j=2,jn
         sin2(j) = sin( -0.5*pi + dy2*(-0.5+float(j-1)) )
      enddo

      call map_a2a(km , im, jm, lon1, sin1, vm,
     &                  in, jn, lon2, sin2, vn, 0, 1)
      call pmaxmin('V_new', vn, pmin, pmax, in*jn, km, 1. )

      write(6,*) ' '
      write(6,*) 'Done horizontal mapping'

! The user may wish to modify nymd, nhms, and nstep here for "branch run"

      write(6,*) ' '
      write(6,*) 'NYMD=',nymd,' NHMS=',nhms
      write(6,*) 'nstep=', nstep

C Fix Dry air mass
      if ( km .eq. 18 ) then
           call set_18(km, ks, ptop, pint, ak, bk)
      elseif( km .eq. 32 ) then
           call set_32(km, ks, ptop, pint, ak, bk)
      elseif( km .eq. 55 ) then
           call set_55(km, ks, ptop, pint, ak, bk)
      else
         write(6,*) 'No suitable set_eta specified'
         stop
      endif

!**************
! Adjust dry mass
!**************
!     call drymadj( in, jn, km, 1, jn,
!    &              .true., ak, bk, psn, delpn, nc, qn, 1 )

!**************
! write new dyn rest file.
!**************
      call rst_dyn(in, jn, km, -1 , iout , delpn, un, vn, ptn, nc,
     &             qn, psn,  nymd,  nhms, nstep )

      return
      end

      subroutine map_a2a(km, im, jm, lon1, sin1, q1,
     &                   in, jn, lon2, sin2, q2, ig, iv)

! Horizontal arbitrary grid to arbitrary grid conservative high-order mapping
! S.-J. Lin

      implicit none
      integer km
      integer im, jm
      integer in, jn
      integer ig                ! ig=0: pole to pole; ig=1 j=1 is half-dy
                                ! north of south pole
      integer iv                ! iv=0: scalar; iv=1: vector
      real q1(im,jm,km)

      real lon1(im+1)
      real lon2(in+1)

      real sin1(jm+1)
      real sin2(jn+1)

! Output
      real q2(in,jn,km)
! local
      integer i,j,k
      real qtmp(in,jm)

      do k=1, km
        if( im .eq. in ) then
            do j=1,jm-ig
               do i=1,im
                  qtmp(i,j+ig) = q1(i,j+ig,k)
               enddo
            enddo
        else
! Mapping in the E-W direction
         call xmap(im, jm-ig, lon1, q1(1,1+ig,k),
     &             in, lon2, qtmp(1,1+ig) )
        endif

        if( jm .eq. jn ) then
            do j=1,jm-ig
               do i=1,in
                  q2(i,j+ig,k) = qtmp(i,j+ig)
               enddo
            enddo
        else
! Mapping in the N-S direction
         call ymap(in, jm, sin1, qtmp(1,1+ig),
     &                 jn, sin2, q2(1,1+ig,k), ig, iv)
        endif
      enddo

      return
      end

      subroutine ymap(im, jm, sin1, q1, jn, sin2, q2, ig, iv)

! Routine to perform area preserving mapping in N-S from an arbitrary
! resolution to another.
!
! sin1 (1) = -1 must be south pole; sin1(jm+1)=1 must be N pole.
!
! sin1(1) < sin1(2) < sin1(3) < ... < sin1(jm) < sin1(jm+1)
! sin2(1) < sin2(2) < sin2(3) < ... < sin2(jn) < sin2(jn+1)
!
! Developer: S.-J. Lin
! First version: piece-wise constant mapping
! Apr 1, 2000
! Last modified:

      implicit none

! Input
      integer im              ! original E-W dimension
      integer jm              ! original N-S dimension
      integer jn              ! Target N-S dimension
      integer ig              ! ig=0: scalars from S.P. to N. P.
                              ! D-grid v-wind is also ig 0
                              ! ig=1: D-grid u-wind
      integer iv              ! iv=0 scalar; iv=1: vector
      real    sin1(jm+1-ig)   ! original southern edge of the cell
                              ! sin(lat1)
!     real    q1(im,jm-ig)    ! original data at center of the cell
      real    q1(im,jm)       ! original data at center of the cell
      real    sin2(jn+1-ig)   ! Target cell's southern edge
                              ! sin(lat2)

! Output
!     real    q2(im,jn-ig)    ! Mapped data at the target resolution
      real    q2(im,jn)       ! Mapped data at the target resolution

! Local
      integer i, j0, m, mm
      integer j

! PPM related arrays
      real   al(im,jm)
      real   ar(im,jm)
      real   a6(im,jm)
      real  dy1(jm)

      real  r3, r23
      parameter ( r3 = 1./3., r23 = 2./3. )
      real pl, pr, qsum, esl
      real dy, sum

      do j=1,jm-ig
         dy1(j) = sin1(j+1) - sin1(j)
      enddo

! ***********************
! Area preserving mapping
! ***********************

! Construct subgrid PP distribution
      call ppm_lat(im, jm, ig, q1, al, ar, a6, 3, iv)

      do 1000 i=1,im
         j0 = 1
      do 555 j=1,jn-ig
      do 100 m=j0,jm-ig
!
! locate the southern edge: sin2(i)
!
      if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then
         pl = (sin2(j)-sin1(m)) / dy1(m)
         if(sin2(j+1) .le. sin1(m+1)) then
! entire new cell is within the original cell
            pr = (sin2(j+1)-sin1(m)) / dy1(m)
            q2(i,j) = al(i,m) + 0.5*(a6(i,m)+ar(i,m)-al(i,m))
     .                    *(pr+pl)-a6(i,m)*r3*(pr*(pr+pl)+pl**2)
               j0 = m
               goto 555
          else
! South most fractional area
            qsum = (sin1(m+1)-sin2(j))*(al(i,m)+0.5*(a6(i,m)+
     .              ar(i,m)-al(i,m))*(1.+pl)-a6(i,m)*
     .               (r3*(1.+pl*(1.+pl))))
              do mm=m+1,jm-ig
! locate the eastern edge: sin2(j+1)
                 if(sin2(j+1) .gt. sin1(mm+1) ) then
! Whole layer
                     qsum = qsum + dy1(mm)*q1(i,mm)
                 else
! North most fractional area
                     dy = sin2(j+1)-sin1(mm)
                    esl = dy / dy1(mm)
                   qsum = qsum + dy*(al(i,mm)+0.5*esl*
     .                   (ar(i,mm)-al(i,mm)+a6(i,mm)*(1.-r23*esl)))
                     j0 = mm
                     goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
555   continue
1000  continue

! Final processing for poles

      if ( ig .eq. 0 .and. iv .eq. 0 ) then

! South pole
           sum = 0.
         do i=1,im
           sum = sum + q2(i,1)
         enddo

           sum = sum / float(im)
         do i=1,im
           q2(i,1) = sum
         enddo

! North pole:
           sum = 0.
         do i=1,im
           sum = sum + q2(i,jn)
         enddo

           sum = sum / float(im)
         do i=1,im
           q2(i,jn) = sum
         enddo

      endif

      return
      end

      subroutine ppm_lat(im, jm, ig, q, al, ar, a6, jord, iv)
      implicit none

!INPUT
! ig=0: scalar pole to pole
! ig=1: D-grid u-wind; not defined at poles because of staggering

      integer im, jm                      !  Dimensions
      integer ig
      real  q(im,jm-ig)
      real al(im,jm-ig)
      real ar(im,jm-ig)
      real a6(im,jm-ig)
      integer jord
      integer iv                             ! iv=0 scalar
                                             ! iv=1 vector
! Local
      real dm(im,jm-ig)
      real    r3
      parameter ( r3 = 1./3. )
      integer i, j, im2, iop, jm1
      real tmp, qmax, qmin
      real qop

! Compute dm: linear slope

      do j=2,jm-1-ig
         do i=1,im
            dm(i,j) = 0.25*(q(i,j+1) - q(i,j-1))
            qmax = max(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
            qmin = q(i,j) - min(q(i,j-1),q(i,j),q(i,j+1))
            dm(i,j) = sign(min(abs(dm(i,j)),qmin,qmax),dm(i,j))
         enddo
      enddo

      im2 = im/2
      jm1 = jm - 1

!Poles:

      if (iv .eq. 1 ) then
!
!*********
! u-wind (ig=1)
! v-wind (ig=0)
!*********
!
! SP
          do i=1,im
              if( i .le. im2) then
                  qop = -q(i+im2,2-ig)
              else
                  qop = -q(i-im2,2-ig)
              endif
              tmp = 0.25*(q(i,2) - qop)
              qmax = max(q(i,2),q(i,1), qop) - q(i,1)
              qmin = q(i,1) - min(q(i,2),q(i,1), qop)
              dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
           enddo
! NP
           do i=1,im
              if( i .le. im2) then
                  qop = -q(i+im2,jm1)
              else
                  qop = -q(i-im2,jm1)
              endif
              tmp = 0.25*(qop - q(i,jm1-ig))
              qmax = max(qop,q(i,jm-ig), q(i,jm1-ig)) - q(i,jm-ig)
              qmin = q(i,jm-ig) - min(qop,q(i,jm-ig), q(i,jm1-ig))
              dm(i,jm-ig) = sign(min(abs(tmp),qmax,qmin),tmp)
           enddo
      else
!
!*********
! Scalar:
!*********
! This code segment currently works only if ig=0
! SP
          do i=1,im2
            tmp = 0.25*(q(i,2)-q(i+im2,2))
            qmax = max(q(i,2),q(i,1), q(i+im2,2)) - q(i,1)
            qmin = q(i,1) - min(q(i,2),q(i,1), q(i+im2,2))
            dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i, 1) =  - dm(i-im2, 1)
          enddo
! NP
          do i=1,im2
            tmp = 0.25*(q(i+im2,jm1)-q(i,jm1))
            qmax = max(q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
            qmin = q(i,jm) - min(q(i+im2,jm1),q(i,jm), q(i,jm1))
            dm(i,jm) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i,jm) =  - dm(i-im2,jm)
          enddo
      endif

      do j=2,jm-ig
        do i=1,im
          al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
        enddo
      enddo

      do j=1,jm-1-ig
        do i=1,im
          ar(i,j) = al(i,j+1)
        enddo
      enddo

      if ( iv .eq. 1 ) then
! Vector:
! ig=0
        if ( ig .eq. 0 ) then
          do i=1,im2
            al(i,    1) = -al(i+im2,2)
            al(i+im2,1) = -al(i,    2)
          enddo

          do i=1,im2
            ar(i,    jm) = -ar(i+im2,jm1)
            ar(i+im2,jm) = -ar(i,    jm1)
          enddo
        else
! ig=1
! SP
          do i=1,im
             if( i .le. im2) then
                 iop = i+im2
             else
                 iop = i-im2
             endif
             al(i,1) = 0.5*(q(i,1)-q(iop,1)) - r3*(dm(iop,1) + dm(i,1))
          enddo
! NP
          do i=1,im
             if( i .le. im2) then
                 iop = i+im2
             else
                 iop = i-im2
             endif
             ar(i,jm1) = 0.5*(q(i,jm1)-q(iop,jm1)) -
     &                   r3*(dm(iop,jm1) + dm(i,jm1))
          enddo
        endif
      else
! Scalar (works for ig=0 only):
          do i=1,im2
            al(i,    1) = al(i+im2,2)
            al(i+im2,1) = al(i,    2)
          enddo

          do i=1,im2
            ar(i,    jm) = ar(i+im2,jm1)
            ar(i+im2,jm) = ar(i,    jm1)
          enddo
      endif

      do j=1,jm-ig
        do i=1,im
          a6(i,j) = 3.*(q(i,j)+q(i,j) - (al(i,j)+ar(i,j)))
        enddo
        call lmppm(dm(1,j), a6(1,j), ar(1,j), 
     &             al(1,j),  q(1,j), im, jord-3)
      enddo

      return
      end

      subroutine xmap(im, jm, lon1, q1, in, lon2, q2)

! Routine to perform area preserving mapping in E-W from an arbitrary
! resolution to another.
! Periodic domain will be assumed, i.e., the eastern wall bounding cell
! im is lon1(im+1) = lon1(1); Note the equal sign is true geographysically.
!
! lon1(1) < lon1(2) < lon1(3) < ... < lon1(im) < lon1(im+1)
! lon2(1) < lon2(2) < lon2(3) < ... < lon2(in) < lon2(in+1)
!
! Developer: S.-J. Lin
! First version: piece-wise constant mapping
! Apr 1, 2000
! Last modified:

      implicit none

! Input
      integer im              ! original E-W dimension
      integer in              ! Target E-W dimension
      integer jm              ! original N-S dimension
      real    lon1(im+1)      ! original western edge of the cell
      real    q1(im,jm)       ! original data at center of the cell
      real    lon2(in+1)      ! Target cell's western edge

! Output
      real    q2(in,jm)       ! Mapped data at the target resolution

! Local
      integer i1, i2
      integer i, i0, m, mm
      integer j

! PPM related arrays
      real qtmp(-im:im+im)
      real   al(-im:im+im)
      real   ar(-im:im+im)
      real   a6(-im:im+im)
      real   x1(-im:im+im+1)
      real  dx1(-im:im+im)
      real  r3, r23
      parameter ( r3 = 1./3., r23 = 2./3. )
      real pl, pr, qsum, esl
      real dx
      integer iord
      data iord /3/
      logical found

      do i=1,im+1
         x1(i) = lon1(i)
      enddo

      do i=1,im
         dx1(i) = x1(i+1) - x1(i)
      enddo

! check to see if ghosting is necessary

!**************
! Western edge:
!**************
          found = .false.
          i1 = 1
      do while ( .not. found )
         if( lon2(1) .ge. x1(i1) ) then
             found = .true.
         else
                  i1 = i1 - 1
             if (i1 .lt. -im) then
                 write(6,*) 'failed in xmap'
                 stop
             else
                 x1(i1) = x1(i1+1) - dx1(im+i1)
                dx1(i1) = dx1(im+i1)
             endif
         endif
      enddo

!**************
! Eastern edge:
!**************
          found = .false.
          i2 = im+1
      do while ( .not. found )
         if( lon2(in+1) .le. x1(i2) ) then
             found = .true.
         else
                  i2 = i2 + 1
             if (i2 .gt. 2*im) then
                 write(6,*) 'failed in xmap'
                 stop
             else
                dx1(i2-1) = dx1(i2-1-im)
                 x1(i2) = x1(i2-1) + dx1(i2-1)
             endif
         endif
      enddo

!     write(6,*) 'i1,i2=',i1,i2

      do 1000 j=1,jm

! ***********************
! Area preserving mapping
! ***********************

! Construct subgrid PP distribution
      call ppm_cycle(im, q1(1,j), al(1), ar(1), a6(1), qtmp(0), iord)

! check to see if ghosting is necessary

! Western edge
          if ( i1 .le. 0 ) then
               do i=i1,0
                  qtmp(i) = qtmp(im+i)
                    al(i) = al(im+i)
                    ar(i) = ar(im+i)
                    a6(i) = a6(im+i)
               enddo
          endif

! Eastern edge:
          if ( i2 .gt. im+1 ) then
             do i=im+1,i2-1
                qtmp(i) = qtmp(i-im)
                  al(i) =   al(i-im)
                  ar(i) =   ar(i-im)
                  a6(i) =   a6(i-im)
             enddo
          endif

         i0 = i1

      do 555 i=1,in
      do 100 m=i0,i2-1
!
! locate the western edge: lon2(i)
!
      if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then
         pl = (lon2(i)-x1(m)) / dx1(m)
         if(lon2(i+1) .le. x1(m+1)) then
! entire new grid is within the original grid
            pr = (lon2(i+1)-x1(m)) / dx1(m)
            q2(i,j) = al(m) + 0.5*(a6(m)+ar(m)-al(m))
     .                    *(pr+pl)-a6(m)*r3*(pr*(pr+pl)+pl**2)
               i0 = m
               goto 555
          else
! Left most fractional area
            qsum = (x1(m+1)-lon2(i))*(al(m)+0.5*(a6(m)+
     .              ar(m)-al(m))*(1.+pl)-a6(m)*
     .               (r3*(1.+pl*(1.+pl))))
              do mm=m+1,i2-1
! locate the eastern edge: lon2(i+1)
                 if(lon2(i+1) .gt. x1(mm+1) ) then
! Whole layer
                     qsum = qsum + dx1(mm)*qtmp(mm)
                 else
! Right most fractional area
                     dx = lon2(i+1)-x1(mm)
                    esl = dx / dx1(mm)
                   qsum = qsum + dx*(al(mm)+0.5*esl*
     .                   (ar(mm)-al(mm)+a6(mm)*(1.-r23*esl)))
                     i0 = mm
                     goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555   continue
1000  continue

      return
      end

      subroutine ppm_cycle(im, q, al, ar, a6, p, iord)
      implicit none

      real r3
      parameter ( r3 = 1./3. )

! Input 
      integer im, iord
      real  q(1)
! Output
      real al(1)
      real ar(1)
      real a6(1)
      real  p(0:im+1)

! local
      real  dm(0:im)
      integer i, lmt
      real tmp, qmax, qmin

         p(0) = q(im)
      do i=1,im
         p(i) = q(i)
      enddo
         p(im+1) = q(1)

! 2nd order slope
      do i=1,im
         tmp = 0.25*(p(i+1) - p(i-1))
         qmax = max(p(i-1), p(i), p(i+1)) - p(i)
         qmin = p(i) - min(p(i-1), p(i), p(i+1))
         dm(i) = sign(min(abs(tmp),qmax,qmin), tmp)
      enddo
         dm(0) = dm(im)

      do i=1,im
         al(i) = 0.5*(p(i-1)+p(i)) + (dm(i-1) - dm(i))*r3
      enddo

      do i=1,im-1
         ar(i) = al(i+1)
      enddo
         ar(im) = al(1)

      if(iord .le. 6) then
         do i=1,im
            a6(i) = 3.*(p(i)+p(i)  - (al(i)+ar(i)))
         enddo
         lmt = iord - 3
         if(lmt.le.2) call lmppm(dm(1),a6(1),ar(1),al(1),p(1),im,lmt)
      else
         call huynh(im, ar(1), al(1), p(1), a6(1), dm(1))
      endif

      return
      end

      subroutine lmppm(dm, a6, ar, al, p, im, lmt)
      implicit none
      real r12
      parameter ( r12 = 1./12. )

      integer im, lmt
      integer i
      real a6(im),ar(im),al(im),p(im),dm(im)
      real da1, da2, fmin, a6da

c LMT = 0: full monotonicity
c LMT = 1: semi-monotonic constraint (no undershoot)
c LMT = 2: positive-definite constraint

      if(lmt.eq.0) then

c Full constraint
      do 100 i=1,im
      if(dm(i) .eq. 0.) then
         ar(i) = p(i)
         al(i) = p(i)
         a6(i) = 0.
      else
         da1  = ar(i) - al(i)
         da2  = da1**2
         a6da = a6(i)*da1
         if(a6da .lt. -da2) then
            a6(i) = 3.*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
         elseif(a6da .gt. da2) then
            a6(i) = 3.*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
         endif
      endif
100   continue

      elseif(lmt.eq.1) then
c Semi-monotonic constraint
      do 150 i=1,im
      if(abs(ar(i)-al(i)) .ge. -a6(i)) go to 150
      if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
            ar(i) = p(i)
            al(i) = p(i)
            a6(i) = 0.
      elseif(ar(i) .gt. al(i)) then
            a6(i) = 3.*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
      else
            a6(i) = 3.*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
      endif
150   continue
      elseif(lmt.eq.2) then
c Positive definite constraint
      do 250 i=1,im
      if(abs(ar(i)-al(i)) .ge. -a6(i)) go to 250
      fmin = p(i) + 0.25*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
      if(fmin.ge.0.) go to 250
      if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
            ar(i) = p(i)
            al(i) = p(i)
            a6(i) = 0.
      elseif(ar(i) .gt. al(i)) then
            a6(i) = 3.*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
      else
            a6(i) = 3.*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
      endif
250   continue
      endif
      return
      end

      subroutine huynh(im, ar, al, p, d2, d1)

! Enforce Huynh's 2nd constraint in 1D periodic domain

      implicit none
      integer im, i
      real ar(im)
      real al(im)
      real  p(im)
      real d2(im)
      real d1(im)

! Local scalars:
      real pmp
      real lac
      real pmin
      real pmax

! Compute d1 and d2
         d1(1) = p(1) - p(im)
      do i=2,im
         d1(i) = p(i) - p(i-1)
      enddo

      do i=1,im-1
         d2(i) = d1(i+1) - d1(i)
      enddo
         d2(im) = d1(1) - d1(im)

! Constraint for AR
!            i = 1
         pmp   = p(1) + 2.0 * d1(1)
         lac   = p(1) + 0.5 * (d1(1)+d2(im)) + d2(im) 
         pmin  = min(p(1), pmp, lac)
         pmax  = max(p(1), pmp, lac)
         ar(1) = min(pmax, max(ar(1), pmin))

      do i=2, im
         pmp   = p(i) + 2.0*d1(i)
         lac   = p(i) + 0.5*(d1(i)+d2(i-1)) + d2(i-1)
         pmin  = min(p(i), pmp, lac)
         pmax  = max(p(i), pmp, lac)
         ar(i) = min(pmax, max(ar(i), pmin))
      enddo

! Constraint for AL
      do i=1, im-1
         pmp   = p(i) - 2.0*d1(i+1)
         lac   = p(i) + 0.5*(d2(i+1)-d1(i+1)) + d2(i+1)
         pmin  = min(p(i), pmp, lac)
         pmax  = max(p(i), pmp, lac)
         al(i) = min(pmax, max(al(i), pmin))
      enddo

! i=im
         i = im
         pmp    = p(im) - 2.0*d1(1)
         lac    = p(im) + 0.5*(d2(1)-d1(1)) + d2(1)
         pmin   = min(p(im), pmp, lac)
         pmax   = max(p(im), pmp, lac)
         al(im) = min(pmax, max(al(im), pmin))

! compute A6 (d2)
      do i=1, im
         d2(i) = 3.*(p(i)+p(i)  - (al(i)+ar(i)))
      enddo
      return
      end

! Parallelized utility routine for computing/printing
! max/min of an input array
!
      subroutine pmaxmin( qname, a, pmin, pmax, im, jm, fac )

#define CPP_PRT_PREFIX

      implicit none

      character*(*)  qname
      integer im, jm
      integer i, j
      real a(im,jm)

      real qmin(jm), qmax(jm)
      real pmax, pmin
      real fac                     ! multiplication factor

#if ( defined OpenMP )
!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j, pmax, pmin)
#endif

#if ( defined IRIX64 )
c$doacross   local(i,j),
c$&          share(im,jm, pmax, pmin)
#endif

      do j=1,jm
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,im
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
            pmax = qmax(1)
            pmin = qmin(1)
         do j=2,jm
            pmax = max(pmax, qmax(j))
            pmin = min(pmin, qmin(j))
         enddo

      CPP_PRT_PREFIX write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

      return
      end

      subroutine set_18(km,ks,ptop,pint,ak,bk)
      implicit none
      integer nl, ks, k, km
      parameter (nl = 18) 
      real ak(km+1),bk(km+1)
      real a(nl+1),b(nl+1)
      real ptop, pint

      data a /291.70,  792.92,  2155.39,  4918.34,  8314.25,
     &       7993.08, 7577.38,  7057.52,  6429.63,  5698.38,
     &       4879.13, 3998.95,  3096.31,  2219.02,  1420.39,
     &        754.13,  268.38,   0.0000,   0.0000 /

      data b /0.0000,  0.0000,   0.0000,   0.0000,   0.0000,
     &        0.0380541, 0.0873088, 0.1489307, 0.2232996,
     &        0.3099406, 0.4070096, 0.5112977, 0.6182465,
     &        0.7221927, 0.8168173, 0.8957590, 0.9533137,
     &        0.9851122, 1.0  /

      if( km .ne. nl) then
          write(6,*) 'Stop in set_18()'
          stop
      endif
 
      ptop =  a(1)
      ks = 4
      pint = a(ks+1) 

      do k=1,nl+1
      ak(k) = a(k)
      bk(k) = b(k)
      enddo

      return
      end

      subroutine set_32(km, ks, ptop, pint, ak, bk)

! New 32-layer setup with top at 0.4 mb for high horizontal resolution runs. 
! pint = 176.93 mb

      implicit none
      integer nl, ks, k, km
      parameter (nl = 32) 
      real ak(km+1)
      real bk(km+1)
      real a(nl+1),b(nl+1)
      real ptop, pint

!GrITAS*Vertical_Levels: 1000 850 700 500 400 300 250 200 150 100 70 50 30 10 5 2 1 0.4

      data a/40.00000,     100.00000,     200.00000,
     &      370.00000,     630.00000,    1000.00000,
     &     1510.00000,    2160.00000,    2900.00000,
     &     3680.00000,    4535.00000,    5505.00000,
     &     6607.26750,    7851.22980,    9236.56610,
     &    10866.34270,   12783.70000,   15039.30000,
     &    17693.00000,   20119.20876,   21686.49129,
     &    22436.28749,   22388.46844,   21541.75227,
     &    19873.78342,   17340.31831,   13874.44006,
     &    10167.16551,    6609.84274,    3546.59643,
     &     1270.49390,       0.00000,       0.00000   /

      data b/0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00696,       0.02801,
     &       0.06372,       0.11503,       0.18330,
     &       0.27033,       0.37844,       0.51046,
     &       0.64271,       0.76492,       0.86783,
     &       0.94329,       0.98511,       1.00000   /

      if( km .ne. nl) then
          write(6,*) 'Stop in set_eta(km=32)'
          stop
      endif
 
      ks = 18
      ptop =  a(1)
      pint = a(ks+1)

      do k=1,nl+1
         ak(k) = a(k)
         bk(k) = b(k)
      enddo

      return
      end

      subroutine set_55(km, ks, ptop, pint, ak, bk)

! Revised 55-level eta with pint at 176.93 mb
! SJL: 2000-03-20

      implicit none
      integer nl, ks, k, km
      parameter (nl = 55)
      real ak(km+1),bk(km+1)
      real a(nl+1)
      real b(nl+1)
      real ptop, pint

      data a/ 1.00000,       2.00000,       3.27000,
     &        4.75850,       6.60000,       8.93450,
     &       11.97030,      15.94950,      21.13490,
     &       27.85260,      36.50410,      47.58060,
     &       61.67790,      79.51340,     101.94420,
     &      130.05080,     165.07920,     208.49720,
     &      262.02120,     327.64330,     407.65670,
     &      504.68050,     621.68000,     761.98390,
     &      929.29430,    1127.68880,    1364.33920,
     &     1645.70720,    1979.15540,    2373.03610,
     &     2836.78160,    3380.99550,    4017.54170,
     &     4764.39320,    5638.79380,    6660.33770,
     &     7851.22980,    9236.56610,   10866.34270,
     &    12783.70000,   15039.30000,   17693.00000,
     &    20119.20876,   21686.49129,   22436.28749,
     &    22388.46844,   21541.75227,   19873.78342,
     &    17340.31831,   13874.44006,   10167.16551,
     &     6609.84274,    3546.59643,    1270.49390,
     &        0.00000,       0.00000   /

      data b / 0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00696,       0.02801,       0.06372,
     &         0.11503,       0.18330,       0.27033,
     &         0.37844,       0.51046,       0.64271,
     &         0.76492,       0.86783,       0.94329,
     &         0.98511,       1.00000  /


      if( km .ne. nl) then
          write(6,*) 'Stop in set_eta(km=55)'
          stop
      endif

         ks = 41
         ptop =  a(1)
         pint = a(ks+1)

      do k=1,NL+1
         ak(k) = a(k)
         bk(k) = b(k)
      enddo

      return
      end

      subroutine drymadj( im, jm, km, jfirst, jlast,
     &                    moun, ak, bk, ps, delp, nq, q, id )

! !USES:

      implicit   none

! !INPUT PARAMETERS:
      integer im, jm, km     ! Dimensions
      integer jfirst, jlast  ! Latitude strip
      integer       id       ! 0:  checking total dry mass only
                             ! 1:  checking total dry mass and adjust
      integer nq             ! Number of tracers         
      logical moun
      real    ak(km+1)
      real    bk(km+1)

! !INPUT/OUTPUT PARAMETERS:
      real delp(im,jfirst:jlast,km)     ! pressure thickness
      real    q(im,jfirst:jlast,km,nq) 
      real   ps(im,jfirst:jlast)        ! surface pressure

! !DESCRIPTION:
!  Perform adjustment of the total dry-air-mass while preserving total
!  tracer mass
!  Developer: S.-J. Lin
!
! !REVISION HISTORY:
!   WS  99.10.26:  Revision; documentation; removed fvcore.h
!   WS  99.11.02:  Limited to jfirst:jlast
!   SJL 00.03.20:
!
!EOP
!---------------------------------------------------------------------
!BOC

      real psd(im,jfirst:jlast)     ! surface pressure  due to dry air mass

      real drym                     ! dry air mass in pascals
      parameter ( drym = 98222. )

      integer   i, j, k
      real  psmo
      real  psdry
      real  gmean
      real  dpd
      integer ic

! Check global maximum/minimum
      PSMO = gmean(im,jm,jfirst,jlast,ps)
      write(6,*) 
     &        'Total (moist) surface pressure before adjustment = ',
     &         PSMO

#if ( defined OpenMP )
!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j,k)
#endif

#if ( defined CRAY )
cmic$*       private(i,j,k)
#endif

#if ( defined IRIX64 ) || ( defined PGI )
c$doacross local(i,j,k)
#endif

      do 1000 j=jfirst,jlast
        do i=1,im
          psd(i,j) = ak(1)
        enddo

        if(nq .ne. 0) then
          do k=1, km
            do i=1,im
              psd(i,j) = psd(i,j) +  delp(i,j,k) * ( 1. - q(i,j,k,1) )
            enddo
          enddo
        else
          do k=1, km
            do i=1,im
              psd(i,j) = psd(i,j) +  delp(i,j,k)
            enddo
          enddo
        endif
1000  continue

      psdry = gmean(im,jm,jfirst,jlast,psd(1,jfirst))
      write(6,*) 'mean (dry) surface pressure = ', psdry
      if( id .eq. 0) return

      if(moun) then
        dpd = drym - psdry
      else
        dpd = 1000.*100. - psdry
      endif

      write(6,*) 'dry mass to be added (pascals) =', dpd

#if ( defined OpenMP )
!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j)
#endif

#if defined( CRAY ) && !defined( CRAY_T3E )
cmic$ do all private(i,j,ic)
#endif

#if ( defined IRIX64 )
c$doacross local(i,j, ic)
#endif

      do 2000 j=jfirst,jlast
        do i=1,im
           ps(i,j) = ps(i,j) + dpd
        enddo

        do ic=1,nq
           do i=1,im
              q(i,j,km,ic) =  q(i,j,km,ic)*delp(i,j,km) /
     &                      ( delp(i,j,km) + dpd )
           enddo
        enddo

        do i=1,im
           delp(i,j,km) = delp(i,j,km) + dpd
        enddo
2000  continue

      PSMO = gmean(im,jm,jfirst,jlast,ps)
      write(6,*) 
     &            'Total (moist) surface pressure after adjustment = ',
     &            PSMO

      return
!EOC
      end

c****6***0*********0*********0*********0*********0*********0**********72
      subroutine rst_dyn(im, jm, km, itp , iuic , delp, u, v, pt, nq,
     &                   q   , ps ,  nymd,  nhms, nstep )
c****6***0*********0*********0*********0*********0*********0**********72

      implicit none
      integer  im, jm
      integer  km

      integer nymd, nhms, nstep
      integer nq, iuic
      integer itp
      integer i, j, k

      real    u(im,jm,km),    v(im,jm,km),  pt(im,jm,km),
     &        delp(im,jm,km),
     &        q(im,jm,km,nq), ps(im,jm)

      real pmax, pmin

      if(itp.ge.0) then

      read(iuic) nstep, nymd, nhms
      read(iuic) ps, delp, u, v, pt

      if(itp .eq. 0) then
         nstep = 0
         write(6,*) 'nstep reset to zero in rst_dyn()'
      endif

      call pmaxmin('U', u, pmin, pmax, im*jm, km, 1. )
      call pmaxmin('V', v, pmin, pmax, im*jm, km, 1. )
      call pmaxmin('PT', pt, pmin, pmax, im*jm, km, 1. )

      if(nq .ne. 0) then
	 read(iuic) q
         call pmaxmin('Q', q, pmin, pmax, im*jm, km, 1. )
      endif

      write(6,*) 'done reading ic for the dycore',nymd, nhms, nstep

	if(nstep .eq. 0) then
	do k=1,km

 	do j=1,jm
 	do i=1,im
 	   q(i,j,k,1) = max(1.2e-12, q(i,j,k,1))
 	enddo
 	enddo

	enddo
	endif

      else

c Write

      rewind iuic
      write(iuic) nstep, nymd, nhms
      write(iuic) ps, delp,u,v,pt
      if(nq .ne. 0) write(iuic) q

      endif

      return
      end

!****6***0*********0*********0*********0*********0*********0**********72
      real function gmean(im,jm,jfirst,jlast,q)
!****6***0*********0*********0*********0*********0*********0**********72
      implicit none

! !INPUT PARAMETERS:
      integer im, jm                       ! Horizontal dimensions
      integer jfirst, jlast                ! Latitude strip
      real q(im,jfirst:jlast)              ! 2D field 

      integer i, j, j2, jmm1
      real dp, acap, area, gsum, xsum
      real sine(jm),cosp(jm),sinp(jm),cose(jm)
      real dl
 
      j2   = max( 2, jfirst )
      jmm1 = min( jm-1, jlast )
      call setrig(im,jm,dp,dl,cosp,cose,sinp,sine)

c scaled polar cap area.
        acap = im*(1.+sine(2)) / dp

        area = 2.*acap
        do j=2,jm-1
          area = area + cosp(j)*im
        enddo
 
      gsum = 0.0
      if ( jfirst .eq. 1  ) gsum = gsum + q(1,1)*acap
      if ( jlast  .eq. jm ) gsum = gsum + q(1,jm)*acap

      do j=j2,jmm1
        xsum = 0.
        do i=1,im
          xsum = xsum + q(i,j)
        enddo
        gsum = gsum + xsum*cosp(j)
      enddo
      gmean = gsum/area

      return
      end


c****6***0*********0*********0*********0*********0*********0**********72
      subroutine setrig(im,jm,dp,dl,cosp,cose,sinp,sine)
c****6***0*********0*********0*********0*********0*********0**********72
      implicit none

      integer im, jm
      integer j, jm1
      real sine(jm),cosp(jm),sinp(jm),cose(jm)
      real dp, dl
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

      write(6,*) 'SETRIG called. ',im,jm
      return
      end

c****6***0*********0*********0*********0*********0*********0**********72
      subroutine vpol5(u,v,im,jm,coslon,sinlon,cosl5,sinl5,jfirst,jlast)
c****6***0*********0*********0*********0*********0*********0**********72
      implicit none

! WS 99.05.05 :  Restructured to split up SP and NP

      integer im, jm, jfirst, jlast
      real u(im,jfirst:jlast),v(im,jfirst:jlast),coslon(im),
     &     sinlon(im),cosl5(im),sinl5(im)

! Local:

      integer i, imh
      real uanp(im), uasp(im), vanp(im), vasp(im)
      real un, vn, us, vs, r2im

! WS 99.05.25 :  Replaced conversions of IMR with IM
      r2im = 0.5d0/dble(im)
      imh  = im / 2

! WS 990726 :  Added condition to decide if poles are on this processor

      IF ( jfirst .EQ. 1 ) THEN

!
! Treat SP
!
      do i=1,im
      uasp(i) = u(i,  2) + u(i,  3)
      enddo

      do i=1,im-1
      vasp(i) = v(i,  2) + v(i+1,  2)
      enddo
      vasp(im) = v(im,  2) + v(1,  2)

c Projection at SP

      us = 0.
      vs = 0.

      do i=1,imh
      us = us + (uasp(i+imh)-uasp(i))*sinlon(i)
     &        + (vasp(i)-vasp(i+imh))*coslon(i)
      vs = vs + (uasp(i+imh)-uasp(i))*coslon(i)
     &        + (vasp(i+imh)-vasp(i))*sinlon(i)
      enddo
      us = us*r2im
      vs = vs*r2im

c get V-wind at SP

      do i=1,imh
      v(i,  1) =  us*cosl5(i) - vs*sinl5(i)
      v(i+imh,  1) = -v(i,  1)
      enddo

      ENDIF

      IF ( jlast .EQ. jm ) THEN

!
! Treat NP
!
      do i=1,im
      uanp(i) = u(i,jm-1) + u(i,jm)
      enddo

      do i=1,im-1
      vanp(i) = v(i,jm-1) + v(i+1,jm-1)
      enddo
      vanp(im) = v(im,jm-1) + v(1,jm-1)

c Projection at NP

      un = 0.
      vn = 0.
      do i=1,imh
      un = un + (uanp(i+imh)-uanp(i))*sinlon(i)
     &        + (vanp(i+imh)-vanp(i))*coslon(i)
      vn = vn + (uanp(i)-uanp(i+imh))*coslon(i)
     &        + (vanp(i+imh)-vanp(i))*sinlon(i)
      enddo
      un = un*r2im
      vn = vn*r2im

c get V-wind at NP

      do i=1,imh
      v(i,jm) = -un*cosl5(i) - vn*sinl5(i)
      v(i+imh,jm) = -v(i,jm)
      enddo

      ENDIF

      return
      end
