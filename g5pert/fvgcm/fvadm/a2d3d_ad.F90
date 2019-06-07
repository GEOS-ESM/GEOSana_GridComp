!-----------------------------------------------------------------------
!BOP
! !ROUTINE: a2d3d_ad -- Adjoint of second order A-to-D grid transformation (3D)
!
! !INTERFACE:

      subroutine a2d3d_ad( ua_ad, va_ad, u_ad, v_ad, im, jm, km, &
                           jfirst, jlast, ng_d, ng_s )

      use precision
      use m_mpif, only:mpi_double_precision,mpi_comm_world
      use mod_comm, only : gid, mp_gather4d

      implicit none

! !INPUT PARAMETERS:

      integer im, jm, km                ! Dimensions
      integer jfirst, jlast             ! Latitude strip
      integer ng_d, ng_s                ! as defined in fvgcm.F
      real(kind=r8), intent(inout) :: u_ad(im,jfirst-ng_d:jlast+ng_s,km)
      real(kind=r8), intent(inout) :: v_ad(im,jfirst-ng_s:jlast+ng_d,km)

! !OUTPUT PARAMETERS:

      real(r8) ua_ad(im,jfirst-ng_d:jlast,km)
      real(r8) va_ad(im,jfirst:jlast,km) 

! !REVISION HISTORY:
!
!    13Jun2007  Todling  Based on Errico/Oloso's code
!
!EOP
!-----------------------------------------------------------------------

    integer, parameter :: ntopfl = 1  ! will pass this as argument
    integer i,j,ierr
    real(r8) :: utotal(im,jm,km)

#ifdef SPMD
    call mp_gather4d(u_ad,utotal,im,jm,km,1,jfirst,jlast, 1,km,ng_d,ng_d,0)
    call mpi_bcast(utotal,im*jm*km,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
#else
    utotal = u_ad
#endif
      do j=max(2,jfirst),min(jm-1,jlast)
        ua_ad(1:im,j,ntopfl:km)=0.5d0*(utotal(1:im,j,ntopfl:km)+utotal(1:im,j+1,ntopfl:km))
      enddo
      if (jfirst == 1) then
        ua_ad(1:im,1,ntopfl:km)=0.5d0*u_ad(1:im,2,ntopfl:km)
      endif
      if (jlast == jm) then
        ua_ad(1:im,jm,ntopfl:km)=0.5d0*u_ad(1:im,jm,ntopfl:km)
      endif
                                                                                                                            
! Adjoint of remapping of v changes back to D grid
! ------------------------------------------------
      va_ad(:,jfirst,:) = 0.0d0
      va_ad(:,jlast ,:) = 0.0d0
      do j=max(2,jfirst),min(jm-1,jlast)
        va_ad(im,j,ntopfl:km)=0.5d0*(v_ad(im,j,ntopfl:km)+v_ad(1,j,ntopfl:km))
        do i=1,im-1
          va_ad(i,j,ntopfl:km)=0.5d0*(v_ad(i,j,ntopfl:km)+v_ad(i+1,j,ntopfl:km))
        enddo
      enddo


      end subroutine a2d3d_ad
