#include <params.h>
      subroutine ocez0(     
     &  lchnk  ,pcols  ,                               
     &  oro    ,tbot   ,pbot  ,taux   ,tauy   ,
     &  z0m    ,z0h    ,zpd   ) 

C----------------------------------------------------------------------- 
C 
C Purpose:  diagnostics of z0m, z0h, and zpd over ocean and seaice
C
C Method:  diagnostic z0 from Charnock's relation (1953) over ocean
C          fluxes and model lowest level fields based on similarity 
C          theory
C Assume:
C   1) Assume z0m = z0h
C   2) Assume zpd = 0.0
C 
C Author: Jiundar Chern, Jan. 2001
C 
C-----------------------------------------------------------------------
      use precision
      implicit none
#include <pmgrid.h>
#include <comcon.h>
!
!-----------------------------------------------------------------------
! note: Constants gravit and rair are defined in comcon.h
!  gravit   = Gravitational acceleration
!  rair     = Gas constant for dry air
!
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer,  intent(in) :: lchnk            ! chunk indentifer
      integer,  intent(in) :: pcols            ! longitude dimension

      real(r8), intent(in) :: oro(pcols)       ! land-sea mask 
      real(r8), intent(in) :: tbot(pcols)      ! Bottom level temperature
      real(r8), intent(in) :: pbot(pcols)      ! Bottom level pressure

      real(r8), intent(in) :: taux(pcols)      ! x surface stress (zonal)
      real(r8), intent(in) :: tauy(pcols)      ! y surface stress (meridional)
!
! Output arguments
!
      real(r8), intent(out) :: z0m(pcols)      ! roughness length, momentum (m)
      real(r8), intent(out) :: z0h(pcols)      ! roughness length, sensible heat (m)
      real(r8), intent(out) :: zpd(pcols)      ! displacement height (m)
!
!---------------------------Local variable-----------------------------
!
      integer i                                ! Longitude indices
      real(r8) rbot                            ! Density at bottom model level
      real(r8) tau                             ! Surface stress
      real(r8) ustar                           ! u*
!
! Loop over Longitude
!
      do i = 1, pcols
        if( oro (i) .eq. 0 ) then
          rbot    = pbot(i) / (rair*tbot(i))
          tau     = sqrt( taux(i) * taux(i) + tauy(i) * tauy(i) )
          ustar   = sqrt( tau / rbot )
          z0m(i)  = 0.016 * ustar * ustar / gravit
          z0h(i)  = z0m(i)
          zpd(i)  = 0.0
        else
! Use ccm value over sea ice
          if( oro(i) .eq. 2 ) then
            z0m(i) = 0.04
            z0h(i) = 0.04
            zpd(i) = 0.0
          end if
        end if
      end do

      return
      end subroutine ocez0
