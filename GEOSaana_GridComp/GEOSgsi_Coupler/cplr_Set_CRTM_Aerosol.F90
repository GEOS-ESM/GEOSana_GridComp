!-----------------------------------------------------------------------------
!BOP
 
! !IROUTINE: set_CRTM_Aerosol --- Set CRTM Aerosols (FORTRAN-77 Interface)
!
! !INTERFACE:
!
subroutine Set_CRTM_Aerosol_ ( km, na, na_dummy, aero_name, aero_conc, rh, aerosol)

! USES:

  use kinds, only: i_kind,r_kind
  use CRTM_Aerosol_Define, only: CRTM_Aerosol_type
  use crtm_aerosol, only: SetAerosol

implicit none

! !ARGUMENTS:

  integer(i_kind) , intent(in)    :: km                ! number of levels
  integer(i_kind) , intent(in)    :: na_dummy          ! to cope w/ NCAR''s change: REVISIT
  integer(i_kind) , intent(in)    :: na                ! number of aerosols
  character(len=*), intent(in)    :: aero_name(na)     ! [na]    GOCART aerosol names: du0001, etc.
  real(r_kind),     intent(in)    :: aero_conc(km,na)  ! [km,na] aerosol concentration (Kg/m2)
  real(r_kind),     intent(in)    :: rh(km)            ! [km]    relative humdity [0,1]

  type(CRTM_Aerosol_type), intent(inout) :: aerosol(na)! [na]   CRTM Aerosol object

! !DESCRIPTION: Set the CRTM Aerosol object given GOCART aerosol properties.
!               This version based on the GEOS-5 implementation of GOCART.
!
! !REVISION HISTORY:
!
! 23feb2011  da Silva  Initial version, FORTRAN-77 interface for GSI.
! 07apr2012  todling   temporarily add na_dummy to cope w/ NCAR''s odd change.
!
!EOP
!-----------------------------------------------------------------------------
 
  call setAerosol (aero_name, aero_conc, rh, aerosol)

end subroutine Set_CRTM_Aerosol_
