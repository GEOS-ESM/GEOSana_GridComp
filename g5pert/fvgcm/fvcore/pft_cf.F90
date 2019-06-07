      subroutine pft_cf(im, jm, js2g0, jn2g0, jn1g1, sc, se, dc, de,   &
                        cosp, cose, ycrit)
      use precision
      implicit none

! Compute coefficients for the 3-point algebraic and the FFT
! polar filters.

! Original version: S.-J. Lin
! Modified for MPI by W. Sawyer and S.J. Lin
 
! Input:
      integer im
      integer jm
      integer js2g0, jn2g0, jn1g1
      real(r8)   cosp(jm), cose(jm)
      real(r8)   ycrit

! Algebric filter
      real(r8)   sc(js2g0:jn2g0)
      real(r8)   se(js2g0:jn1g1)

! FFT filter
      real(r8)   dc(im,js2g0:jn2g0)
      real(r8)   de(im,js2g0:jn1g1)

! Local:
      integer i, j
      real(r8)   dl, coszc, cutoff, phi, damp
      real(r8)  pi
!     real(r8) deg

      pi = 4.d0 * datan(1.d0)

      coszc = cos(ycrit*pi/180.)

! INIT fft polar coefficients:
      dl = pi/dble(im)
      cutoff = 1.e-20

      do j=js2g0,jn2g0
         do i=1,im
            dc(i,j) = 1.
         enddo
      enddo

      do j=js2g0,jn1g1
         do i=1,im
            de(i,j) = 1.
         enddo
      enddo

!     write(6,*) '3-point polar filter coefficients:'

!************
! Cell center
!************
!     do 100 j=2,jm-1
      do 100 j=js2g0,jn2g0
            sc(j) = (coszc/cosp(j))**2

#if defined ( ALT_PFT )
         if( sc(j) .gt. 1. ) then
#else
         if(sc(j) .gt.1. .and. sc(j) .le. 2.0) then
            sc(j) =  1. +  (sc(j)-1.)/(sc(j)+1.)
         elseif(sc(j) .gt. 2. .and. sc(j) .le. 4.) then
            sc(j) =  1. +  sc(j)/(8.-sc(j))
         elseif(sc(j) .gt. 4. ) then
#endif

! FFT filter
            do i=1,im/2
               phi = dl * i
               damp = min((cosp(j)/coszc)/sin(phi), 1._r8)**2
               if(damp .lt. cutoff) damp = 0.
               dc(2*i-1,j) = damp
               dc(2*i  ,j) = damp
            enddo

         endif
!     deg = -90. + (j-1)*180./(jm-1)
!     write(6,*) j, deg, sc(j)
100   continue

!************
! Cell edges
!************
!     do 200 j=2,jm
      do 200 j=js2g0,jn1g1
            se(j) = (coszc/cose(j))**2

#if defined ( ALT_PFT )
         if( se(j) .gt. 1. ) then
#else
         if(se(j) .gt.1. .and. se(j) .le. 2.0 ) then
            se(j) =  1. +  (se(j)-1.)/(se(j)+1.)
         elseif(se(j) .gt. 2. .and. se(j) .le. 4.) then
            se(j) =  1. +  se(j)/(8.-se(j))
         elseif(se(j) .gt. 4. ) then
#endif
! FFT
            do i=1,im/2
               phi = dl * i
               damp = min((cose(j)/coszc)/sin(phi), 1._r8)**2
               if(damp .lt. cutoff) damp = 0.
               de(2*i-1,j) = damp
               de(2*i  ,j) = damp
            enddo
         endif
200   continue

      return
      end
