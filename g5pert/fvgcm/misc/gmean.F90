module gmean

      use precision
      implicit none

      real(r8), allocatable :: gw(:)

  public    :: gmean_initialize, gmean_finalize, gmean_of

contains

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  gmean_initialize --- intialize gmean
!
! !INTERFACE:
   subroutine gmean_initialize(im, jm)
! !USES:


! !INPUT PARAMETERS:
  integer  im, jm                        ! Horizontal dimensions

! !DESCRIPTION:
!
!EOP
!---------------------------------------------------------------------
!BOC
  real(r8) sine(jm), cosp(jm), sinp(jm), cose(jm)
  real(r8) dl, dp
  integer j

  call setrig(im,jm,dp,dl,cosp,cose,sinp,sine)

  allocate( gw(jm) )
  gw(1) = 1 + sine(2)
  do j=2,jm-1
     gw(j) = sine(j+1) - sine(j)
  enddo
  gw(jm) = 1 - sine(jm)
  return

!EOC
  end subroutine gmean_initialize
!---------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  gmean_finalize --- finalize gmean
!
! !INTERFACE:
   subroutine gmean_finalize
! !USES:

!EOC

  deallocate( gw )

  return
!EOC

  end subroutine gmean_finalize
!---------------------------------------------------------------------






!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  gmean_of --- Calculate the mean of a 2D field
!
! !INTERFACE:
      function gmean_of(im, jm, jfirst, jlast, q)
! !USES:

! FastOpt fixed bug gmean is used as real(kind=r8)
      real(r8) gmean_of

! !INPUT PARAMETERS:
      integer  im, jm                        ! Horizontal dimensions
      integer  jfirst, jlast                 ! Latitude strip
      real(r8) q(im,jfirst:jlast)            ! 2D field 

! !DESCRIPTION:
!     Calculate the mean of a 2D field
!
!EOP
!---------------------------------------------------------------------
!BOC

      integer i, j
      real(r8) xsum(jfirst:jlast)

      xsum = 0.
      do j=jfirst,jlast
        do i=1,im
          xsum(j) = xsum(j) + q(i,j)
        enddo
          xsum(j) = xsum(j)*gw(j)
      enddo

      call par_vecsum( jm, jfirst, jlast, xsum, gmean_of)
      gmean_of = gmean_of / (2*im)

      return
!EOC
      end function gmean_of
!---------------------------------------------------------------------


end module gmean
