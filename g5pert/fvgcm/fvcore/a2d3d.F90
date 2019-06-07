!-----------------------------------------------------------------------
!BOP
! !ROUTINE: a2d3d -- Second order A-to-D grid transformation (3D)
!
! !INTERFACE:
      subroutine a2d3d(ua, va, u, v, im, jm, km,            &
                       jfirst, jlast, ng_d, ng_s,           &
                       coslon, sinlon)

      use precision
#if defined( SPMD )
      use mod_comm, only : mp_send_s, mp_recv_n, mp_barrier, mp_send2_n, mp_recv2_s
#endif
      implicit none
! !INPUT PARAMETERS:
      integer im, jm, km                  ! Dimensions
      integer jfirst, jlast               ! Latitude strip
      integer ng_d, ng_s                  ! Ghosting parameters
      real(r8) ua(im,jfirst-1:jlast,km)   ! U-Wind on A Grid
      real(r8) va(im,jfirst:jlast,km)     ! V-Wind on A Grid
      real(r8) coslon(im),sinlon(im)      ! Sine and cosine in longitude

! !OUTPUT PARAMETERS:
      real(r8) u(im,jfirst-ng_d:jlast+ng_s,km)
      real(r8) v(im,jfirst-ng_s:jlast+ng_d,km) 

! !REVISION HISTORY:
!
!    13Jun2007  Todling  Based on Errico/Oloso's code.
!
!EOP
!-----------------------------------------------------------------------

      integer i, j, j1, k

!     Compute u at jfirst later; need ua(j=jfirst-1)
!     ----------------------------------------------
#if defined( SPMD )
      call mp_barrier
      call mp_send2_n(im,jm,jfirst,jlast,1,km,1,0,ua,ua)
      call mp_barrier
      call mp_recv2_s(im,jm,jfirst,jlast,1,km,1,0,ua,ua)
#endif

      j1=max(2,jfirst)
      do j=j1,jlast
          u(1:im,j,1:km) = 0.5d0 * (ua(1:im,j,1:km)+ua(1:im,j-1,1:km))
      enddo
                                                                                                                    
!     Add change to v by remapping changes from D grid
!     ------------------------------------------------
      do j=max(2,jfirst),min(jm-1,jlast)
         v(1,j,1:km) = 0.5d0 * (va(1,j,1:km)+va(im,j,1:km))
         do i=2,im
            v(i,j,1:km) = 0.5d0 * (va(i,j,1:km)+va(i-1,j,1:km))
         enddo
      enddo

      return
      end subroutine a2d3d
