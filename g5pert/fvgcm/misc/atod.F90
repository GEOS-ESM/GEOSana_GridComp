!-----------------------------------------------------------------------
!BOP
! !IROUTINE: atod --- Convert wind field from A to D grid
!
! !INTERFACE:
      subroutine atod(ua,va,u,v,im,jm,coslon,sinlon)

! !USES:
      use precision
      implicit none

! !INPUT PARAMETERS:
      integer(i4), intent(in) :: im          ! X dimension
      integer(i4), intent(in) :: jm          ! Y dimension
      real(r8), intent(in)    :: coslon(im)  ! cosine of longitude
      real(r8), intent(in)    :: sinlon(im)  ! sine of longitude

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: ua(im,jm)
      real(r8), intent(inout) :: va(im,jm)

! !OUTPUT PARAMETERS:
      real(r8), intent(out)   :: u(im,jm)
      real(r8), intent(out)   :: v(im,jm)

! !DESCRIPTION:
!     Convert winds from A grid to D grid
!
!
! !REVISION HISTORY:
!   99.09.30   Lin     Creation
!   00.07.07   Sawyer  Removed fvcore.h; use precision module; cleaning
!
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      real(r8) un, vn, us, vs
      integer(i4)   imh, i, j
!
! **********************************************
! Get consistent A-grid winds at poles.
!
      imh = im / 2
! NP
      un = 0.
      vn = 0.
! SP
      us = 0.
      vs = 0.
      j=jm-1
      do i=1,imh
        un = un + (ua(i+imh,j)-ua(i,j))*sinlon(i)  &
                + (va(i+imh,j)-va(i,j))*coslon(i)  
        vn = vn + (ua(i,j)-ua(i+imh,j))*coslon(i)  &
                + (va(i+imh,j)-va(i,j))*sinlon(i)
        us = us + (ua(i+imh,2)-ua(i,2))*sinlon(i)  &
                + (va(i,2)-va(i+imh,2))*coslon(i)
        vs = vs + (ua(i+imh,2)-ua(i,2))*coslon(i)  &
                + (va(i+imh,2)-va(i,2))*sinlon(i)
      enddo
!
      un = un/im
      vn = vn/im
      us = us/im
      vs = vs/im
!
      do i=1,imh
! SP
        ua(i,1)   = -us*sinlon(i) - vs*coslon(i)
        va(i,1)   =  us*coslon(i) - vs*sinlon(i)
        ua(i+imh,1)   = -ua(i,1)
        va(i+imh,1)   = -va(i,1)
! NP
        ua(i,jm) = -un*sinlon(i) + vn*coslon(i)
        va(i,jm) = -un*coslon(i) - vn*sinlon(i)
        ua(i+imh,jm) = -ua(i,jm)
        va(i+imh,jm) = -va(i,jm)
      enddo
! **********************************************
!
      do j=2,jm
        do i=1,im
          u(i,j) = 0.5*(ua(i,j) + ua(i,j-1))
        enddo
      enddo
!
      do j=2,jm-1
        do i=2,im
          v(i,j) = 0.5*(va(i,j) + va(i-1,j))
        enddo
      enddo
!
      do j=2,jm-1
        v(1,j) = 0.5*(va(im,j) + va(1,j))
      enddo
      return
!EOC
      end subroutine atod
!-----------------------------------------------------------------------

