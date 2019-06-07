#include <params.h>
      subroutine diagbma(     
     &  lchnk,  pcols,                               
     &  thbot,  qbot,   ps,   ts,   oro,   shf,   lhf,   wet1,  bma  )

!----------------------------------------------------------------------- 
! 
! Purpose:  diagnostics of bulk moisture availability 
!
! Method:  Retrieve bulk moisture availability fromm known 
!          fluxes and model lowest level fields based on similarity theory
! 
! Author: Jiundar Chern, Jan. 2001
! 
!-----------------------------------------------------------------------
      use precision

      implicit none

#include <pmgrid.h>

#include <comcon.h>

!
!-----------------------------------------------------------------------
! note: Constants cpvir,cpair,latvap, and latice are 
!       defined in comcon.h
!  cpvir    = cpwv/cpair -1
!  cpwv     = Specific heat of water vapor
!  cpair    = Specific heat of dry air
!  latvap   = Latent heat of vaporization
!  latice   = Latent heat of fusion
!
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer,  intent(in) :: lchnk            ! chunk indentifer
      integer,  intent(in) :: pcols            ! longitude dimension

      real(r8), intent(in) :: thbot(pcols)     ! Bottom level potential temperature
      real(r8), intent(in) :: qbot(pcols)      ! Bottom level specific humidity
      real(r8), intent(in) :: ps(pcols)        ! surface pressure (Pa)
      real(r8), intent(in) :: ts(pcols)        ! surface temperature (K)
      real(r8), intent(in) :: oro(pcols)       ! land-sea mask 
      real(r8), intent(in) :: shf(pcols)       ! sensible heat flux (W/m2)
      real(r8), intent(in) :: lhf(pcols)       ! latent heat flux (W/m2)
      real(r8), intent(in) :: wet1(pcols)      ! first soil layer wetness (fraction)
!
! Output arguments
!
      real(r8), intent(out) :: bma(pcols)      ! bulk moisture availability (fraction)
!
!---------------------------Local variable-----------------------------
!
      integer i                                ! Longitude indices
      real(r8) qsat(pcols)                     ! saturation specific humidity at surface
      real(r8) est(pcols)                      ! saturation vapor pressure at surface
      real(r8) cp                              ! Specific heat of moist air
      real(r8) ltheat                          ! Latent heat for given srf conditions
      real(r8) tmp1,tmp2                       ! working space
!
!---------------------------------------------------------------
! Compute saturation specific humidity at surface
!---------------------------------------------------------------
!
         call aqsat(ts, ps, est, qsat, pcols, pcols, 1, 1, 1 )
!
! Loop over Longitude
!
      do i = 1, pcols
!
!---------------------------------------------------------------
! Set up necessary constants
!---------------------------------------------------------------
!
         cp      = cpair * ( 1. + cpvir * qbot(i) ) 
         if( ts(i) > 273.15 ) then
           ltheat = latvap
         else
           ltheat = latvap + latice
         end if
!
!---------------------------------------------------------------
! Compute bulk moisture avaliability (bma)
!---------------------------------------------------------------
! set bma to ocean value
!
         bma(i) = 0.99
         if( oro(i) == 1 ) then
           bma(i) = wet1(i)
           tmp1   = shf(i) * ( thbot(i) - ts(i) )
           tmp2   = lhf(i) * ( qbot(i) - qsat(i) )
           if( tmp1 < 0.0 .and. tmp2 < 0.0 ) then
             bma(i) =  cp * ( thbot(i) - ts(i) ) * lhf(i)
     &              / ( ltheat * ( qbot(i) - qsat(i) ) * shf(i) )
             if( bma(i) > 0.9 .or. bma(i) < 0.1 ) bma(i) = wet1(i)
           end if
         else
           if( oro(i) == 2 ) bma(i) = 1.0
         end if
      end do

      return
      end subroutine diagbma
