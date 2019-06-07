!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  gmean4 --- Calculate the mean of a 2D field
!
! !INTERFACE:
      real*4 function gmean4(im, jm, jfirst, jlast, q)
! !USES:
      use precision
      implicit none

! !INPUT PARAMETERS:
      integer  im, jm                        ! Horizontal dimensions
      integer  jfirst, jlast                 ! Latitude strip
      real(r4) q(im,jfirst:jlast)            ! 2D field 

! !DESCRIPTION:
!     Calculate the mean of a 2D field
!
!EOP
!---------------------------------------------------------------------
!BOC

      integer i, j
      real(r8) sine(jm),cosp(jm),sinp(jm),cose(jm)
      real(r8), allocatable, save :: gw(:)
      real(r8) dl, dp, xsum(jfirst:jlast)
      real(r8) gmean
      logical first
      data first /.true./

      if(first) then
         call setrig(im,jm,dp,dl,cosp,cose,sinp,sine)

         allocate( gw(jm) )
            gw(1) = 1 + sine(2)
         do j=2,jm-1
            gw(j) = sine(j+1) - sine(j)
         enddo
            gw(jm) = 1 - sine(jm)
 
        first = .false.
      endif

      do j=jfirst,jlast
         xsum(j)= 0.
        do i=1,im
          xsum(j) = xsum(j) + q(i,j)
        enddo
          xsum(j) = xsum(j)*gw(j)
      enddo

      call par_vecsum( jm, jfirst, jlast, xsum, gmean)
      gmean4 = gmean / (2*im)

      return
!EOC
      end
!---------------------------------------------------------------------
