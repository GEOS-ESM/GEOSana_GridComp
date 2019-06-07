! $Id$


#include "MAPL_Generic.h"

module AeroOptPropTableMod
  use ESMF
  use MAPL_Mod
  use Chem_MieMod

  implicit none
  private

  type(Chem_Mie),  save :: MieTable
  logical       ,  save :: MieTableCreated=.false.

  public Get_AeroOptProp
  public Get_AeroOptPropAllAero
  public Get_AeroIndex
  public Get_AeroIndexTableProps
  public Get_AeroTauList
  public Create_AeroOptProp

  interface Get_AeroOptPropAllAero
     module procedure Get_AeroOptPropAllAero3D
     module procedure Get_AeroOptPropAllAero4D
  end interface Get_AeroOptPropAllAero

contains

!=================================================================

  subroutine Get_AeroOptProp(idx,ib,qaero,rh,tau,ssa,asy,rc)

    integer,           intent(IN ) :: idx
    integer,           intent(IN ) :: ib
    real,              intent(IN ) :: qaero
    real,              intent(IN ) :: rh
    real,              intent(OUT) :: tau
    real,              intent(OUT) :: ssa
    real,              intent(OUT) :: asy
    integer, optional, intent(OUT) :: rc

    integer :: status
    character(len=ESMF_MAXSTR) :: Iam='Get_AeroOptProp' 

    call Chem_MieQuery(MieTable,idx,real(ib),qaero,rh, &
         tau,ssa,asy,rc=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Get_AeroOptProp

!=================================================================

  subroutine Get_AeroOptPropAllAero3D(na,aerosol,nb,offset,qaero,rh,tau,ssa,asy,rc)

    integer, intent(in)            :: na         ! number of aerosols
    character(len=*),  intent(IN ) :: aerosol(:) ! list of aerosols
    integer,           intent(IN ) :: nb
    integer,           intent(IN ) :: offset
    real,              intent(IN ) :: qaero(:,:,:)
    real,              intent(IN ) :: rh(:,:)
    real,              intent(OUT) :: tau(:,:,:)
    real,              intent(OUT) :: ssa(:,:,:)
    real,              intent(OUT) :: asy(:,:,:)
    integer, optional, intent(OUT) :: rc

    integer :: status
    character(len=ESMF_MAXSTR) :: Iam='Get_AeroOptPropAllAero3D' 

    integer :: l, idx

    real*8 :: tau8(size(tau,1),size(tau,2),size(tau,3))
    real*8 :: ssa8(size(tau,1),size(tau,2),size(tau,3))
    real*8 :: asy8(size(tau,1),size(tau,2),size(tau,3))

    tau8 = 0.0
    ssa8 = 0.0
    asy8 = 0.0

    do l = 1, na
       idx = Chem_MieQueryIdx( MieTable, aerosol(l), STATUS)
       VERIFY_(STATUS)

       call Chem_MieQueryAllBand3D(MieTable,idx,nb,offset,qaero(:,:,l),rh, &
            tau,ssa,asy,rc=STATUS)
       VERIFY_(STATUS)

       tau8 = tau8 +              tau
       ssa8 = ssa8 +       (ssa * tau)
       asy8 = asy8 + asy * (ssa * tau)
    end do

    tau = tau8
    ssa = ssa8
    asy = asy8

    RETURN_(ESMF_SUCCESS)
  end subroutine Get_AeroOptPropAllAero3D

!=================================================================

  subroutine Get_AeroOptPropAllAero4D(na,aerosol,nb,offset,qaero,rh,tau,ssa,asy,rc)

    integer,           intent(IN ) :: na         ! number of aerosols
    character(len=*),  intent(IN ) :: aerosol(:) ! list of aerosols
    integer,           intent(IN ) :: nb
    integer,           intent(IN ) :: offset
    real,              intent(IN ) :: qaero(:,:,:,:)
    real,              intent(IN ) :: rh(:,:,:)
    real,              intent(OUT) :: tau(:,:,:,:)
    real,              intent(OUT) :: ssa(:,:,:,:)
    real,              intent(OUT) :: asy(:,:,:,:)
    integer, optional, intent(OUT) :: rc

    integer :: status
    character(len=ESMF_MAXSTR) :: Iam='Get_AeroOptPropAllAero4D' 

    integer :: l, idx

    real*8 :: tau8(size(tau,1),size(tau,2),size(tau,3),size(tau,4))
    real*8 :: ssa8(size(tau,1),size(tau,2),size(tau,3),size(tau,4))
    real*8 :: asy8(size(tau,1),size(tau,2),size(tau,3),size(tau,4))

    tau8 = 0.0
    ssa8 = 0.0
    asy8 = 0.0

    do l = 1, na
       idx = Chem_MieQueryIdx( MieTable, aerosol(l), STATUS)
       VERIFY_(STATUS)

       call Chem_MieQueryAllBand4D(MieTable,idx,nb,offset,qaero(:,:,:,l),rh, &
            tau,ssa,asy,rc=STATUS)
       VERIFY_(STATUS)

       tau8 = tau8 +              tau
       ssa8 = ssa8 +       (ssa * tau)
       asy8 = asy8 + asy * (ssa * tau)
    end do

    tau = tau8
    ssa = ssa8
    asy = asy8

    RETURN_(ESMF_SUCCESS)
  end subroutine Get_AeroOptPropAllAero4D

!=================================================================

  subroutine Get_AeroTauList(idx,ib,nn,qaero,rh,tau,rc)

    integer,           intent(IN ) :: idx
    integer,           intent(IN ) :: ib,nn
    real,              intent(IN ) :: qaero(nn)
    real,              intent(IN ) :: rh(nn)
    real,              intent(OUT) :: tau(nn)
    integer, optional, intent(OUT) :: rc

    integer :: status
    character(len=ESMF_MAXSTR) :: Iam='Get_AeroTauList' 


    call Chem_MieQueryTauList(MieTable,idx,real(ib),qaero,rh, &
         tau,rc=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Get_AeroTauList

!=================================================================

  integer function Get_AeroIndex(aerosol, rc ) result (idx)
    character (len=*), intent(IN ) :: aerosol
    integer, optional, intent(OUT) :: rc

    integer :: status
    character(len=ESMF_MAXSTR) :: Iam='Get_AeroIndex' 


    idx = Chem_MieQueryIdx ( MieTable, AEROSOL, STATUS ) 
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end function Get_AeroIndex

  subroutine Get_AeroIndexTableProps(na,aerosol,rh,routine,&
                                     idx,nrh,bext,bsca,gasy,arh,irh,rc)
    
    integer, intent(in) :: na                  ! number of aerosols
    character(len=*), intent(in) :: aerosol(:) ! list of aerosols
    real, intent(in) :: rh(:,:)                ! relative humidity
    character(len=5), intent(in) :: routine    ! name of routine

    integer, intent(out) :: idx(:)     ! subindex of aerosol table
    integer, intent(out) :: nrh(:)     ! number of RH levels in table
    real, intent(out) :: bext(:,:,:,:) ! mass extinction efficiency
    real, intent(out) :: bsca(:,:,:,:) ! mass scattering efficiency
    real, intent(out) :: gasy(:,:,:,:) ! asymmetry parameter
    real, intent(out) :: arh(:,:,:)    ! mapped slope from hi-res RH map
    integer, intent(out) :: irh(:,:,:) ! mapped pointer from hi-res RH map
    integer, optional, intent(out) :: rc

    integer :: status, i,j,k, ii, jj, isnap
    integer :: num_radii, num_rh
    character(len=ESMF_MAXSTR) :: Iam='Get_AeroIndexTableProps' 

    ii = size(rh,1)
    jj = size(rh,2)

    ! Initialize outgoing matrices
    idx = 0
    nrh = 0
    bext = 0.0
    bsca = 0.0
    gasy = 0.0
    arh = 0.0
    irh = 0

    do k = 1, na

      idx(k) = Chem_MieQueryIdx( MieTable, aerosol(k), STATUS)
      VERIFY_(STATUS)

      nrh(k) = MieTable%vtableUse%nrh

      num_rh = size(MieTable%vtableUse%bext,2)
      num_radii = size(MieTable%vtableUse%bext,3)

      select case(routine)
      case('SOLAR') ! Use only the first 8 tables for SOLAR
         bext(1:8,1:num_rh,1:num_radii,k) = MieTable%vtableUse%bext(1:8,:,:)
         bsca(1:8,1:num_rh,1:num_radii,k) = MieTable%vtableUse%bsca(1:8,:,:)
         gasy(1:8,1:num_rh,1:num_radii,k) = MieTable%vtableUse%g(   1:8,:,:)
      case('IRRAD') ! Use only the last 10 tables for IRRAD
         bext(1:10,1:num_rh,1:num_radii,k) = MieTable%vtableUse%bext(9:18,:,:)
         bsca(1:10,1:num_rh,1:num_radii,k) = MieTable%vtableUse%bsca(9:18,:,:)
         gasy(1:10,1:num_rh,1:num_radii,k) = MieTable%vtableUse%g(   9:18,:,:)
      case default
         ASSERT_(.FALSE.) ! We MUST have a case
      end select

      do j = 1, jj
         do i = 1, ii
            isnap = int(((min(max(rh(i,j),0.0),0.99))+0.001)*1000)

            arh(i,j,k) = MieTable%vtableUse%rha(isnap)
            irh(i,j,k) = MieTable%vtableUse%rhi(isnap)
         end do
      end do

    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine Get_AeroIndexTableProps

  subroutine Create_AeroOptProp(cf,rc)

    type(ESMF_Config), intent(IN ) :: cf
    integer, optional, intent(OUT) :: rc

    integer :: status
    character(len=ESMF_MAXSTR) :: Iam='Create_AeroOptProp' 

    MieTable = Chem_MieCreate(cf,status)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Create_AeroOptProp


end module AeroOptPropTableMod
