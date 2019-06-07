!-----------------------------------------------------------------------
!BOP
! !ROUTINE: get_srf --- Read in surface data
!
! !INTERFACE:
      subroutine get_srf(im, jm, jfirst, jlast, iuhs, srfdata,
     &                   moun, phis, sgh, grav)

      implicit none

! !INPUT PARAMETERS:
      integer im, jm              ! Dimensions
      integer jfirst, jlast
      character*80 srfdata        ! Fsle name
      integer  iuhs               ! Unit number
      real grav                   ! Gravitation constant
      logical moun                ! Read real data?
                                  ! If false, set output to 0

! !OUTPUT PARAMETERS:
      real  phis(im,jfirst:jlast)     ! surface geopotential
      real   sgh(im,jfirst:jlast)     ! std dev. of topg. (m) C

! !DESCRIPTION:
!     Read in the surface data (replicated on all PEs in MP mode)
!
! !REVISION HISTORY:
!   SJL 99.01.01:  Delivery
!   WS  99.09.07:  Parallel interface, now calls readsrf
!   WS  99.11.19:  Documentation; removed imr, jnp
!   SJL 99.12.31:  Use pmaxmin
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local variables

      integer  i,  j, jt
      real pmax, pmin
      real gz(im,jfirst:jlast)
      logical exist

      jt = jlast - jfirst + 1

      if(moun) then

! Pass in global (jfirst=1, jlast=jm) for now
         call readsrf(iuhs, srfdata, im, jm, jfirst, jlast, phis, gz)

! Standard deviation is multiplied by 2 according to John T. @ NCAR
! No need to do this in ccm3.6
        do j=jfirst,jlast
          do i=1,im
             sgh(i,j) = 2.0*gz(i,j)
          enddo
        enddo

      else
        do j=jfirst,jlast
          do i=1,im
            phis(i,j) = 0.
             sgh(i,j) = 0.
          enddo
        enddo
      endif

! Check initial topography
      call pmaxmin('Hs  (m)', phis, pmin, pmax, im, jt, 1./grav)
      call pmaxmin('SGH (m)',  sgh, pmin, pmax, im, jt, 1.)

      return
      end
