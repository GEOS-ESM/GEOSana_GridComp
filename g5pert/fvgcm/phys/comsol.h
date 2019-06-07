!
!	Solar radiation
!
!	$Id$
!
! Visible optical depth
!
      real(r8) tauvis     ! Visible optical depth

      common /comvis/ tauvis
!
! Solar constant
!
      real(r8) scon       ! Solar constant

      common /comsol/ scon
!
! Earth orbital characteristics
!	
      real(r8) eccen       ! eccentricity factor (unitless) (typically 0 to 0.1)
      real(r8) obliq       ! obliquity angle (degrees) (-90 to +90) (typically 22-26)
      real(r8) mvelp       ! moving vernal equinox at perhelion (degrees) (0 to 360.0)
      integer iyear_AD     ! Year (AD) to simulate above orbital parameters for
!
! Orbital information after processed by orbit_params
!
      real(r8) obliqr      ! obliquity in radians
      real(r8) lambm0      ! Mean longitude of perihelion at the 
!                          ! vernal equinox (radians)
      real(r8) mvelpp      ! moving vernal equinox longitude
!                          ! of perihelion plus pi (radians)
!
      common /comorb/ eccen   , obliq   , mvelp   , obliqr  
      common /comorb/ lambm0  , mvelpp  , iyear_AD

