!----------------------------------------------------------------------------
!BOP
!  
! !MODULE:  GSI_CRTMcoupler_AeroMod ---
!
! !INTERFACE:
!
! !DESCRIPTION: This module provides general interface for
!               Aerosol fields seen by CRTM.
!
! !REVISION HISTORY:
!
!  2026-02-04 Todling - Revive this code (unsure how this got lost over time)
!
!EOP
!-------------------------------------------------------------------------

module GSI_CRTMCoupler_AeroMod

implicit none
private

!
! !PUBLIC MEMBER FUNCTIONS:
!
public Set_CRTM_Aerosol

interface Set_CRTM_Aerosol
  subroutine Set_CRTM_Aerosol_ ( km, na, na_crtm, aero_name, aero_conc, rh, aerosol)
    use kinds,         only: r_kind, i_kind
    use CRTM_Aerosol_Define, only: CRTM_Aerosol_type
    implicit none

    integer(i_kind) , intent(in)    :: km                ! number of levels
    integer(i_kind) , intent(in)    :: na                ! number of aerosols
    integer(i_kind) , intent(in)    :: na_crtm           ! number of aerosols seen by CRTM
    character(len=*), intent(in)    :: aero_name(na)     ! [na]    GOCART aerosol names
    real(r_kind),     intent(inout) :: aero_conc(km,na)  ! [km,na] aerosol concentration (Kg/m2)
    real(r_kind),     intent(in)    :: rh(km)            ! [km]    relative humdity [0,1]

    type(CRTM_Aerosol_type), intent(inout) :: aerosol(na_crtm)! [na]   CRTM Aerosol object
     
  end subroutine Set_CRTM_Aerosol_
end interface

end module GSI_CRTMCoupler_AeroMod

