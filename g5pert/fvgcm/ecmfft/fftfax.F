      subroutine fftfax (n, ifax, trigs)

      use precision

#if defined( SGI_FFT )
!
! WS 99.11.24 : Added SGI wrappers
!
      implicit none
      real(r8)    trigs(1)
      integer ifax(*)
      integer n
! local
      integer*4 nn

      nn=n
      call dzfftm1dui (nn,trigs)
#else
c      dimension ifax(13),trigs(3*n/2+1)
       integer ifax(13)
       real(r8) trigs(3*n/2+1)
c
c MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
c TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
c DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
c WAS WRITTEN.
c
      data mode /3/
      call fax (ifax, n, mode)
      i = ifax(1)
      if (ifax(i+1) .gt. 5 .or. n .le. 4) ifax(1) = -99
      if (ifax(1) .le. 0 ) then
        write(6,*) ' set99 -- invalid n'
        stop'set99'
      endif
      call fftrig (trigs, n, mode)
#endif
      return
      end
