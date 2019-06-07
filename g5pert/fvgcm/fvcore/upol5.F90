!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  upol5 --- Treat U winds at poles
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine upol5(u,v,im,jm,coslon,sinlon,cosl5,sinl5,jfirst,jlast)
!****6***0*********0*********0*********0*********0*********0**********72

! !USES:
      use precision
      implicit none

! !INPUT PARAMETERS:
      integer im                       ! Total longitudes
      integer jm                       ! Total latitudes
      integer jfirst                   ! First PE latitude (no ghosting)
      integer jlast                    ! Last  PE latitude (no ghosting)
      real(r8) coslon(im), sinlon(im)
      real(r8) cosl5(im),sinl5(im)

      real(r8) v(im,jfirst:jlast)      ! Winds in Y (C-grid)

! !INPUT/OUTPUT PARAMETERS:
      real(r8) u(im,jfirst:jlast)      ! Winds in X (C-grid)

!
! !DESCRIPTION:
!
!   Treat the U winds at the poles.  This requires an average 
!   of the U- and V-winds, weighted by their angles of incidence
!   at the pole points.     
!
! !REVISION HISTORY:
!   98.06.01 Lin       Creation
!   99.05.05 Sawyer    Restructured to split up SP and NP
!   99.05.25 Sawyer    Replaced conversions of IMR, JMR ; removed fvcore.h
!   99.07.26 Sawyer    Added condition to decide if poles are on this PE
!   01.03.26 Sawyer    Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      integer i, imh
      real(r8) uanp(im), uasp(im), vanp(im), vasp(im)
      real(r8) un, vn, us, vs, r2im

      imh = im/2
      r2im = 0.5d0/dble(im)

      if ( jfirst .EQ. 1 ) then
!
! Treat SP
!
      do i=1,im-1
      uasp(i) = u(i,  2) + u(i+1,2)
      enddo
      uasp(im) = u(im,  2) + u(1,2)

      do i=1,im
      vasp(i) = v(i,  2) + v(i,  3)
      enddo

! Projection at SP

      us = 0.
      vs = 0.
      do i=1,imh
      us = us + (uasp(i+imh)-uasp(i))*sinlon(i)     &
              + (vasp(i)-vasp(i+imh))*coslon(i)
      vs = vs + (uasp(i+imh)-uasp(i))*coslon(i)     &
              + (vasp(i+imh)-vasp(i))*sinlon(i)
      enddo
      us = us*r2im
      vs = vs*r2im

! get U-wind at SP

      do i=1,imh
      u(i,  1) = -us*sinl5(i) - vs*cosl5(i)
      u(i+imh,  1) = -u(i,  1)
      enddo

      endif

      if ( jlast .EQ. jm ) then
!
! Treat NP
!
      do i=1,im-1
      uanp(i) = u(i,jm-1) + u(i+1,jm-1)
      enddo
      uanp(im) = u(im,jm-1) + u(1,jm-1)

      do i=1,im
      vanp(i) = v(i,jm-1) + v(i,jm)
      enddo

! Projection at NP

      un = 0.
      vn = 0.
      do i=1,imh
      un = un + (uanp(i+imh)-uanp(i))*sinlon(i)  &
              + (vanp(i+imh)-vanp(i))*coslon(i)
      vn = vn + (uanp(i)-uanp(i+imh))*coslon(i)  &
              + (vanp(i+imh)-vanp(i))*sinlon(i)
      enddo
      un = un*r2im
      vn = vn*r2im

! get U-wind at NP

      do i=1,imh
      u(i,jm) = -un*sinl5(i) + vn*cosl5(i)
      u(i+imh,jm) = -u(i,jm)
      enddo

      endif

      return
!EOC
      end subroutine upol5
!----------------------------------------------------------------------- 
