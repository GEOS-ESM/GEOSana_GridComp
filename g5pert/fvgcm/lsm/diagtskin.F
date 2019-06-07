!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
! !IROUTINE:  diagtskin --- Diagnostics of surface skin temperature
!
! !INTERFACE:
!

      subroutine diagtskin (npt     ,veg    ,emv  ,emg  ,tv  ,tg  ,
     &                      tskin4  ,trad4  )
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
!
      integer, intent(in) :: npt            !number of points
      logical, intent(in) :: veg(npt)       !true if vegetated surface
      real,    intent(in) :: emv(npt)       !vegetation emissivity
      real,    intent(in) :: emg(npt)       !ground emissivity
      real,    intent(in) :: tv(npt)        !vegetation temperature (kelvin)
      real,    intent(in) :: tg(npt)        !ground temperature (kelvin)
!
! !OUTPUT PARAMETERS:
!
      real,    intent(out)  :: tskin4(npt)  ! 4th power of skin temperature
      real,    intent(out)  :: trad4(npt)   ! 4th power of radiative temperature
!
! !DESCRIPTION: compute skin and radiative temperature for fvDAS
!
! !REVISION HISTORY:
!
! 28Jan2002  Chern     Initial code
!
!EOP
!-------------------------------------------------------------------------
!  local variables 
      integer k
      real    emveff          ! effective vegetation emissivity
      real    emgeff          ! effective ground emissivity
      real    tmp1, tmp2      ! temporary working space
      real    irv             ! 4th power of radiative temperature from vegetation
      real    irg             ! 4th power of radiative temperature from ground
      real    fv              ! fraction skin temperature contribution from vegetation
      
      do k = 1, npt
        if( veg(k) ) then
          emveff = emv(k) * ( 1.0 + (1.0 - emv(k)) * (1.0 - emg(k)) )
          emgeff = emg(k) * ( 1.0 -emv(k) )
          tmp1   = tv(k) ** 4
          tmp2   = tg(k) ** 4
          irv    = emveff * tmp1
          irg    = emgeff * tmp2
        else
          emveff = 0.0
          emgeff = emg(k)
          tmp1   = 0.0
          tmp2   = tg(k) ** 4
          irv    = 0.0
          irg    = emgeff * tmp2
        end if
        trad4(k) = irv + irg
        fv       = irv / (irv + irg)
        tskin4(k) = fv * tmp1 + (1.0 - fv) * tmp2
      end do

      return
      end subroutine diagtskin
