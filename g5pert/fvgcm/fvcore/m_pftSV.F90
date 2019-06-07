

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_pftSV --- Module for polar filter of initial-time SVs
!
! !INTERFACE:
!
module m_pftSV

!
! !PUBLIC MEMBER FUNCTIONS:
!

  PRIVATE
  PUBLIC  pftSV

!
! !DESCRIPTION: Polar filter initial-time SVs by applying a damping 
!               factor that depends on latitude and zonal wavenumber.
!               The filtering is performed in zonal spectral space. 
!               See the routine pftSV_damp for the damping factors.
!
!   This routine shopuld be applied both at the beginning of the TLM
!   and end of the ADM, when (iouter.eq.1 .and. iinner.eq.1) 
!
! !REVISION HISTORY:
!
!     Ronald Errico  April 15, 2004: Initial code
!     Nathan Winslow May 1,    2004: Insertion into FVGCM
!     Ronald Errico  August 20, 2004: Entirely new algoritm and code

!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
  subroutine pftSV (imr,jnp,nl,jfirst,jlast,cosp,cose,pt,delp,u,v)

!
!  Sequence routines to filter all fields
!  Set lfilter=.false. if no filtering is desired
!  Note the filter is global: filtering can occur at all latitudes, 
!    depending on how the array damp is set.
!
      use precision
      implicit none

! Input/Output
      integer, intent(in) :: imr,jnp,nl    
      integer, intent(in) :: jfirst,jlast
      real(kind=r8), intent(in)    :: cosp(jfirst:jlast)
      real(kind=r8), intent(in)    :: cose(jfirst:jlast)
      real(kind=r8), intent(inout) :: pt(imr,jfirst:jlast,nl)
      real(kind=r8), intent(inout) :: delp(imr,jfirst:jlast,nl)
      real(kind=r8), intent(inout) :: u(imr,jfirst:jlast,nl)
      real(kind=r8), intent(inout) :: v(imr,jfirst:jlast,nl)
!
! Local variables

      integer       :: k         ! vertical index
      integer       :: j1        ! first lat to be processed (never=1)
      integer       :: jp        ! last lat on p grid processed (never =jnp)
      integer       :: je        ! last lat on edge grid processed 
      integer       :: jnump     ! number of lats on p grid to process
      integer       :: jnume     ! number of lats on edge grid to process
      integer       :: ifax(13)  ! prime factorization of imr for FFT
      real(kind=r8) :: trigs(3*imr/2+1) ! trig factors used by FFT
      real(kind=r8), allocatable  :: dampp(:,:) ! damping coefs for p grid
      real(kind=r8), allocatable  :: dampe(:,:) ! damping coefs for edge grid

      logical, parameter :: lfilter=.true.  ! true if filter is on 
                                            ! otherwise routine does nothing
!
!
      if (lfilter) then  ! perform filter
!
! Determine which latitudes to act upon (not pole points)
        j1=max(2,jfirst)
        jp=min(jnp-1,jlast)
        je=jlast
        jnump=jp-j1+1
        jnume=je-j1+1
!
! Determine FFT transformation factors and damping factors
        call fftfax (imr,ifax,trigs)
        allocate ( dampp(imr,j1:jp) )
        allocate ( dampe(imr,j1:je) )
        call pftSV_damp (imr,j1,jp,cosp(j1:jp),dampp(:,j1:jp))
        call pftSV_damp (imr,j1,je,cose(j1:je),dampe(:,j1:je))
!
! Filter each field on each level separately
        do k=1,nl
          call pft2dSV (delp(:,j1:jp,k),dampp(:,j1:jp), &
                        imr,jnump,ifax,trigs)
          call pft2dSV (  pt(:,j1:jp,k),dampp(:,j1:jp), &
                        imr,jnump,ifax,trigs)
          call pft2dSV (   u(:,j1:je,k),dampe(:,j1:je), &
                        imr,jnume,ifax,trigs)
          call pft2dSV (   v(:,j1:jp,k),dampp(:,j1:jp), &
                        imr,jnump,ifax,trigs)
        enddo  ! loop over levels
!
        deallocate ( dampp )
        deallocate ( dampe )
!
    endif     ! test on whether to perform filter

  end subroutine pftSV

!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
      subroutine pftSV_damp (imr,j1,j2,coslat,damp)
      use precision
!
! Compute coefficients for Fourier filter of initial fields in the
! tangent linear model or final fields in the adjoint model.
!
! The filtering is done based on a function of the zonal wavenumber k,
! accounting for the cosine factor in the length of the latitude circle.
! The function here is:
!
!      damp = 1.                           for k <= kmin
!           = (kmin/k)**3                  for k > kmin
!
! where kmin is the k corresponding to a 2*delta x wave at the equator
! and k is a wavenumber at the latitude being filtered.
!
! The indexes for the array damp(i,j) are j for latitude and i for
! dimensionless zonal wavenumber m=int((i+1)/2). Odd and even i correspond
! respectively to real and imaginary components of the Fourier
! coefficients (those damping coefs should be identical). m is the 
! wavenumber measured in multiples of the fundamental wavenumber (that 
! for the largest zonal wavelength at each particular latitude).
!    
! damp varies between 0 and 1. Smaller values indicate greater filtering.
! At the precise poles, there is no filtering, since conditions there are 
! constrained by using other conditions.
!   
! This routine should never be called for the poles, since coslat=0 then.
!
!
! Input:
      integer,       intent(in)  :: imr
      integer,       intent(in)  :: j1,j2
      real(kind=r8), intent(in)  :: coslat(j1:j2)    ! cosine of latitiude

! Output:
      real(kind=r8), intent(out) :: damp(imr,j1:j2)  ! filter factor 

! Local:
      real,    parameter :: power=2.     ! tuning parameter 
      integer     :: i     ! wavenumber or lomgitude index
      integer     :: j     ! index for latitude
      real(r8)    :: kmin  ! wavenumber corr. to 2 delta x at equator
      real(r8)    :: kfund ! wavenumber 1 at current latitude
      real(r8)    :: k
      real(r8)    :: fk
!
! k here is the zonal wavenumber, made dimensionless by multiplication
! by the Earth's radius a.
!
!  k(dimensional)   = imr/(i*a*cos(j)), where i is number of grid points
!                     corresponding to wavelength
!
      kmin=real(imr)/2.d0
!
      do j=j1,j2
        kfund=1.d0/coslat(j)  ! corresp. to wavenumber 1 at latitude j
        do i=2,imr,2
          k=kfund*i/2 ! divide by 2 here because i/2 is wave number
          if (k > kmin) then
            fk=(kmin/k)**power
            damp(i,j)=fk
          else      
            damp(i,j)=1.d0
          endif
          damp(i-1,j)=damp(i,j)
        enddo  ! loop over wavenumber index
      enddo    ! loop over latitudes 
!
      end subroutine pftSV_damp

!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
 subroutine pft2dSV (p,damp,im,jp,ifax,trigs)

! Apply particular polar filter for calculation of singular vectors to
! a single 2-d (lat-lon) field.
! All filtering performed by applying zonal FFT to fields, multiplying 
! zonal spectral coefficients by a supplied factor 0 <= damp <= 1, and 
! projecting back onto the grid, replacing the original field.
!

 use precision
 implicit none

! Input
     integer,  intent(in) :: im  ! number of points in zonal dirction
     integer,  intent(in) :: jp  ! number of latitudes to process
     integer,  intent(in) :: ifax(13)        ! factorization of im 
     real(r8), intent(in) :: trigs(3*im/2+1) ! trig factors for FFT
     real(r8), intent(in) :: damp(im,jp)     ! spectral damping coefs

! Input/Output
     real(r8), intent(inout) :: p(im,jp)   ! Field to be polar filtered

! Local
     integer   :: i, j
     real(r8)  :: q1(im+2,jp) ! array for copying fields
     real(r8)  :: q2(im+2,jp) ! work array for FFT

! Copy field into an array with two extra elements for each j
! as required for the FFT routine

      do j=1,jp
        do i=1,im
          q1(i,j) = p(i,j)
        enddo
        q1(im+1,j) = 0.
        q1(im+2,j) = 0.
      enddo
!
! Replace field values in q1 by Fourier coefficients
! These coefs are ordered real, imaginary, real, imaginary, ...
! for dimensionless wavenumber 0,0, 1,1, 2,2, ..., im/2, im/2)
! 
      call rfftmlt(q1, q2, trigs, ifax, 1, im+2, im, jp, -1)
!
! Damp spectral coefficients.
! Note that the array damp is presumed to be ordered the same as 
! the spectral coefficients except that the values start corresponding 
! to wavenumber 1, rather than at wavenumber 0 as for q1.

      do j=1,jp
        do i=3,im+2
          q1(i,j) = q1(i,j) * damp(i-2,j)
        enddo
      enddo
!
! Compute filtered fields from filtered coefficients
      call rfftmlt(q1, q2, trigs, ifax, 1, im+2, im, jp, 1)
!
! Copy fields from array q1 back to p
      do j=1,jp
        do i=1,im
          p(i,j) = q1(i,j)
        enddo
      enddo

 end subroutine pft2dSV

 end module m_pftSV







