!-----------------------------------------------------------------------------
!BOP

! !MODULE: crtm_aerosol --- Implements the GOCART-CRTM Aerosol Interface
!
! !INTERFACE:
!

module crtm_aerosol

! !USES:

  use kinds, only: i_kind,r_kind
  use CRTM_Aerosol_Define, only: CRTM_Aerosol_type
  use CRTM_Parameters,     only: DUST_AEROSOL, SEASALT_SSAM_AEROSOL, &
                                 SEASALT_SSCM1_AEROSOL, SEASALT_SSCM2_AEROSOL, &
                                 SEASALT_SSCM3_AEROSOL, SEASALT_SSCM3_AEROSOL, &
                                 BLACK_CARBON_AEROSOL, ORGANIC_CARBON_AEROSOL, &
                                 SULFATE_AEROSOL


  use Chem_RegistryMod, only: Chem_Registry, Chem_RegistryCreate, Chem_RegistryDestroy
  use Chem_MieMod,      only: Chem_Mie, Chem_MieCreate, Chem_MieDestroy, &
                              Chem_MieQueryIdx, Chem_MieQuery
  use m_chars,          only: lowercase
  use m_die,            only: die

  implicit none

! !PUBLIC METHODS:

  public setAerosol

! !PUBLIC DATA MEMBERS:
! --------------------
  type(Chem_Registry), pointer :: aerReg => null()  ! Aerosol Registry
  type(Chem_Mie),      pointer :: Mie => null()     ! Mie tables

! !REVISION HISTORY:
!
! 23feb2011  da Silva  Initial version.
!
!EOP
!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------
!BOP

! !IROUTINE: setAerosol --- Set CRTM Aerosol Object
!
! !INTERFACE:
!

  subroutine setAerosol (aero_name, aero_conc, rh, aerosol)

! !ARGUMENTS:

  character(len=*), intent(in)    :: aero_name(:)      ! [na]    GOCART aerosol names: du0001, etc.
  real(r_kind),     intent(in)    :: aero_conc(:,:)    ! [km,na] aerosol concentration (Kg/m2)
  real(r_kind),     intent(in)    :: rh(:)             ! [km]    relative humdity [0,1]

  type(CRTM_Aerosol_type), intent(inout) :: aerosol(:) ! [na]   CRTM Aerosol object

! !DESCRIPTION: Set the CRTM Aerosol object given GOCART aerosol properties.
!               This version based on the GEOS-5 implementation of GOCART.
!
! !REVISION HISTORY:
!
! 23feb2011  da Silva  Initial version.
!
!EOP
!-----------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'setAerosol'
  integer :: na, km, n, k, rc, iTable
  real :: rh_, rEff_

  km = size(aero_conc,1)
  na = size(aero_conc,2)

! Load Mie Tables the first time around
! -------------------------------------
  if ( .not. associated(aerReg) ) then
     allocate(aerReg,Mie,stat=rc)
     if ( rc /= 0 ) &
          call die(myname,'cannot allocate Mie, Registry')
     aerReg = Chem_RegistryCreate(rc, 'Chem_AerRegistry.rc')
     if ( rc /= 0 ) &
          call die(myname,'cannot create aerosol registry; make sure <Chem_AerRegistry.rc> is available.')
     Mie = Chem_MieCreate('Chem_Mie-550nm.rc', rc, aerReg )
     if ( rc /= 0 ) &
        call die(myname,"cannot load Mie Tables; make sure <Chem_Mie-550nm.rc> is available and that ExtData symlink is in place.")
  end if

! Loop over species...
! --------------------
  do n = 1, na

!    Map GOCART names into CRTM Aerosol indices
!    ------------------------------------------
     Aerosol(n)%Type = AeroType_(aero_name(n))

!    Concentration
!    -------------
     Aerosol(n)%Concentration(:) = aero_conc(:,n)

!    Effective radius
!    ----------------
     iTable = Chem_MieQueryIdx(Mie, aero_name(n), rc)
     if(iTable == -1 .OR. rc/=0 ) &
        call die(myname,"cannot get Mie index for "//aero_name(n) )
     do k = 1, km
        rh_ = rh(k)
        call Chem_MieQuery(Mie, iTable, 1.0, 0.0, rh_, &
                           rEff=rEff_ )
        Aerosol(n)%Effective_radius(k) = rEff_
     end do
  end do

CONTAINS

  function AeroType_(name) Result(atype)
    character(len=*) :: name  ! GOCART name
    integer(i_kind)  :: atype ! CRTM aerosol type

    if ( lowercase(name(1:2)) == 'du' ) then
       atype = DUST_AEROSOL

    else if ( trim(lowercase(name)) == 'ss001' ) then
       atype = SEASALT_SSAM_AEROSOL
    else if ( trim(lowercase(name)) == 'ss002' ) then
       atype = SEASALT_SSCM1_AEROSOL
    else if ( trim(lowercase(name)) == 'ss003' ) then
       atype = SEASALT_SSCM2_AEROSOL
    else if ( trim(lowercase(name)) == 'ss004' ) then
       atype = SEASALT_SSCM3_AEROSOL
    else if ( trim(lowercase(name)) == 'ss005' ) then
       atype = SEASALT_SSCM3_AEROSOL

    else if ( lowercase(name(1:2))  ==    'bc' ) then
       atype = BLACK_CARBON_AEROSOL

    else if ( lowercase(name(1:2))  ==    'oc' ) then
       atype = ORGANIC_CARBON_AEROSOL

    else if ( trim(lowercase(name)) ==   'so4' ) then
       atype = SULFATE_AEROSOL

    else
       call die(myname,"cannot recognize aerosol name <"//trim(name)//">")
    end if

  end function AeroType_

end subroutine setAerosol

end module crtm_aerosol
