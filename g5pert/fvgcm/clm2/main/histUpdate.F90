#include <misc.h>
#include <preproc.h>

subroutine histUpdate()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! update history file fields
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use clmtype
  use clm_varder
  use clm_varmap  , only : begpatch, endpatch
  use clm_varcon  , only : tfrz, istdlak, istslak
#if (defined DGVM)
  use shr_const_mod, only: SHR_CONST_CDAY
  use accumulMod  , only : accslf, accext
  use pft_varcon  , only : pftpar
  use time_manager, only : get_step_size, get_nstep
#endif
  use histFileMod , only : histslf, histmlf, spval
  implicit none

! ------------------- local variables -----------------------------
  integer j,k   !loop indices
  real(r8) tmpmlev(begpatch:endpatch,1:nlevsoi)
  real(r8) tmpslev(begpatch:endpatch)
#if (defined DGVM)
  integer :: dtime               !timestep size [seconds]
  integer :: nstep               !timestep number
#endif
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! accumulate field values over history time interval
! -----------------------------------------------------------------

! atmospheric input fields

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%forc_rain
  end do
  call histslf ('RAIN    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%forc_snow
  end do
  call histslf ('SNOW    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%forc_t
  end do
  call histslf ('TBOT    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%forc_th
  end do
  call histslf ('THBOT   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%forc_q
  end do
  call histslf ('QBOT    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%forc_hgt
  end do
  call histslf ('ZBOT    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%forc_lwrad
  end do
  call histslf ('FLDS    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) = sqrt(clm(k)%forc_u**2+clm(k)%forc_v**2)
  end do
  call histslf ('WIND    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) = clm(k)%forc_solad(1) + clm(k)%forc_solai(1) &
                 +clm(k)%forc_solad(2) + clm(k)%forc_solai(2)
  end do
  call histslf ('FSDS    ', tmpslev)

! snow properties

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%snowdp
  end do
  call histslf ('SNOWDP  ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%snowage
  end do
  call histslf ('SNOWAGE ', tmpslev)

! soil properties

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevsoi
        if (clm(k)%itypwat==istdlak .or. clm(k)%itypwat==istslak) then
           tmpmlev(k,j) = spval
        else
           tmpmlev(k,j) = clm(k)%z(j)
        endif
     end do
  end do
  call histmlf ('ZSOI   ', tmpmlev, nlevsoi)

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevsoi
        if (clm(k)%itypwat==istdlak .or. clm(k)%itypwat==istslak) then
           tmpmlev(k,j) = spval
        else
           tmpmlev(k,j) = clm(k)%dz(j)
        endif
     end do
  end do
  call histmlf ('DZSOI  ', tmpmlev, nlevsoi)

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevsoi
        tmpmlev(k,j) = clm(k)%watsat(j)
     end do
  end do
  call histmlf ('WATSAT ', tmpmlev, nlevsoi)

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevsoi
        tmpmlev(k,j) = clm(k)%sucsat(j)
     end do
  end do
  call histmlf ('SUCSAT ', tmpmlev, nlevsoi)  

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevsoi
        tmpmlev(k,j) = clm(k)%bsw(j)
     end do
  end do
  call histmlf ('BSW    ', tmpmlev, nlevsoi)

! water content

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%h2osno 
  end do
  call histslf ('H2OSNO  ', tmpslev)  

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%h2ocan
  end do
  call histslf ('H2OCAN  ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevsoi
        tmpmlev(k,j) = clm(k)%h2osoi_vol(j)
     end do
  end do
  call histmlf ('H2OSOI  ', tmpmlev, nlevsoi)

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevsoi
        tmpmlev(k,j) = clm(k)%h2osoi_liq(j)
     end do
  end do
  call histmlf ('SOILLIQ ', tmpmlev, nlevsoi)

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevsoi
        tmpmlev(k,j) = clm(k)%h2osoi_ice(j)
     end do
  end do
  call histmlf ('SOILICE ', tmpmlev, nlevsoi)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%snowliq
  end do
  call histslf ('SNOWLIQ ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%snowice
  end do
  call histslf ('SNOWICE ', tmpslev)

! temperatures

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%t_veg
  end do
  call histslf ('TV      ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%t_grnd
  end do
  call histslf ('TG      ', tmpslev)  

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevsoi
        tmpmlev(k,j) = clm(k)%t_soisno(j)
     end do
  end do
  call histmlf ('TSOI    ', tmpmlev, nlevsoi)

!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch,endpatch
     do j = 1,nlevlak
        tmpmlev(k,j) = clm(k)%t_lake(j)
     end do
  end do
  call histmlf ('TLAKE   ', tmpmlev, nlevlak)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%t_ref2m
  end do
  call histslf ('TSA     ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%t_snow
  end do
  call histslf ('TSNOW   ', tmpslev)

! canopy physiology 

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%rssun
  end do
  call histslf ('RSSUN   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%rssha
  end do
  call histslf ('RSSHA   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%btran
  end do
  call histslf ('BTRAN   ', tmpslev)

! vegetation phenology 

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%elai
  end do
  call histslf ('ELAI    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%esai
  end do
  call histslf ('ESAI    ', tmpslev)

! surface solar radation

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%fsa
  end do
  call histslf ('FSA     ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%fsr
  end do
  call histslf ('FSR     ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%ndvi
  end do
  call histslf ('NDVI    ', tmpslev)

! energy fluxes

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%eflx_lh_vegt
  end do
  call histslf ('FCTR    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%eflx_lh_vege
  end do
  call histslf ('FCEV    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%eflx_lh_grnd
  end do
  call histslf ('FGEV    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%eflx_sh_tot
  end do
  call histslf ('FSH     ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%eflx_snomelt
  end do
  call histslf ('FSM     ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%eflx_soil_grnd
  end do
  call histslf ('FGR     ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%eflx_lwrad_net
  end do
  call histslf ('FIRA    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%eflx_lwrad_out
  end do
  call histslf ('FIRE    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%taux
  end do
  call histslf ('TAUX    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%tauy
  end do
  call histslf ('TAUY    ', tmpslev)

! water fluxes

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_infl
  end do
  call histslf ('QINFL   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_surf
  end do
  call histslf ('QOVER   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_qrgwl
  end do
  call histslf ('QRGWL   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_drain
  end do
  call histslf ('QDRAI   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_evap_soi
  end do
  call histslf ('QSOIL   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_evap_veg - clm(k)%qflx_tran_veg
  end do
  call histslf ('QVEGE   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_tran_veg
  end do
  call histslf ('QVEGT   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_prec_grnd
  end do
  call histslf ('QDRIP   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_prec_intr
  end do
  call histslf ('QINTR   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%qflx_snomelt
  end do
  call histslf ('QMELT   ', tmpslev)

! conservation checks

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%errsoi
  end do
  call histslf ('ERRSOI  ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%errseb
  end do
  call histslf ('ERRSEB  ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  abs(clm(k)%errseb)
  end do
  call histslf ('ERRSEBMX', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%errsol
  end do
  call histslf ('ERRSOL  ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%errh2o
  end do
  call histslf ('ERRH2O  ', tmpslev)

#if (defined BGC)
! Biogeochemistry

!  call histslf ('FMICR   ', clm(begpatch:endpatch)%fmicr)
!  call histslf ('FPSN    ', clm(begpatch:endpatch)%fpsn )
!  call histslf ('FRMF    ', clm(begpatch:endpatch)%frmf )
!  call histslf ('FRMS    ', clm(begpatch:endpatch)%frms )
!  call histslf ('FRMR    ', clm(begpatch:endpatch)%frmr )
!  call histslf ('FRG     ', clm(begpatch:endpatch)%frg  )
!  call histslf ('FCO2    ', clm(begpatch:endpatch)%fco2 )
!  call histslf ('DMI     ', clm(begpatch:endpatch)%dmi  )

!$OMP PARALLEL DO PRIVATE (K)
   do k = begpatch,endpatch
      tmpslev(k) = 0._r8
      do j = 1, ndst
         tmpslev(k) = tmpslev(k) + clm(k)%flx_mss_vrt_dst(j) 
      end do
   end do
   call histslf ('DSTFLXT ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
   do k = begpatch,endpatch
      tmpslev(k) = 0._r8
      do j = 1, nvoc
         tmpslev(k) = tmpslev(k) + clm(k)%vocflx(j)
      end do
   end do
   call histslf ('VOCFLXT ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
   do k = begpatch,endpatch
      tmpslev(k) = clm(k)%vocflx(1)
   end do
   call histslf ('ISOPRENE', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
   do k = begpatch,endpatch
      tmpslev(k) = clm(k)%vocflx(2)
   end do
   call histslf ('MONOTERP', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
   do k = begpatch,endpatch
      tmpslev(k) = clm(k)%vocflx(3)
   end do
   call histslf ('OVOC    ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
   do k = begpatch,endpatch
      tmpslev(k) = clm(k)%vocflx(4)
   end do
   call histslf ('ORVOC   ', tmpslev)

!$OMP PARALLEL DO PRIVATE (K)
   do k = begpatch,endpatch
      tmpslev(k) = clm(k)%vocflx(5)
   end do
   call histslf ('BIOGENCO', tmpslev)
#endif

#if (defined DGVM)
! -----------------------------------------------------------------
! internally accumulate fields values over history time interval
! -----------------------------------------------------------------

! determine time step size and index

  dtime = get_step_size()
  nstep = get_nstep()

! Examples of various field accumulations for use within model
! TDA = accumulates TBOT as 30-day average (orig. was daily avg of TSA)
! T10 = accumulates TSA as 10-day running mean
! FNPSN10 = accumulates fpsn-frmf as 10-day running mean
! PREC365 = accumulates total precipitation as 365-day running mean

  call accslf ('TDA     ', clm(begpatch:endpatch)%forc_t , nstep)
  call accslf ('T10     ', clm(begpatch:endpatch)%t_ref2m, nstep)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%fpsn - clm(k)%frmf
  end do
  call accslf ('FNPSN10 ', tmpslev, nstep)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  clm(k)%forc_rain + clm(k)%forc_snow
  end do
  call accslf ('PREC365 ', tmpslev, nstep)

! Examples of how to accumulate growing degree days based on 
! 10-day running mean temperature. First, extract t10 from
! accumulated variables. Then accumulate GDD above 0C and -5C.
! The trigger to reset the accumulated values to zero is -99999.

  call accext ('T10     ', clm(begpatch:endpatch)%t10, nstep)

! slevis: in orig. tmpslev = (clm%t_da - ...) but t_da = 1.e36 between days!

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  (clm(k)%t10 - tfrz) * dtime / SHR_CONST_CDAY
     if (tmpslev(k) < 0.) tmpslev(k) = -99999.
  end do
  call accslf ('AGDD0   ', tmpslev, nstep)

! slevis: in orig. tmpslev = (clm%t_da - ...) but t_da = 1.e36 between days!

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) =  (clm(k)%t10 - (tfrz - 5.0)) * dtime / SHR_CONST_CDAY
     if (tmpslev(k) < 0.) tmpslev(k) = -99999.
  end do
  call accslf ('AGDD5   ', tmpslev, nstep)

! slevis: in orig. tmpslev = (clm%t_da - ...) but t_da = 1.e36 between days!
!         also agddtw was accumulated only when warmest month temperature
!         exceeded twmax = pftpar(clm%itypveg,31)
!         agddtw is currently reset at the end of each year in subr. lpj

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) = max(0.0, (clm(k)%t10 - (tfrz + pftpar(clm(k)%itypveg,31))) * dtime / SHR_CONST_CDAY)
  end do
  call accslf ('AGDDTW  ', tmpslev, nstep)

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     tmpslev(k) = max(0.0, (clm(k)%t_ref2m - (tfrz + 5.0)) * dtime / SHR_CONST_CDAY)
  end do
  call accslf ('AGDD    ', tmpslev, nstep)

! extract accumulated variables to send to history file

  call accext  ('TDA     ', clm(begpatch:endpatch)%t_mo   , nstep)
  call histslf ('TDA     ', clm(begpatch:endpatch)%t_mo   )

  call accext  ('AGDD0   ', clm(begpatch:endpatch)%agdd0  , nstep)
  call histslf ('AGDD0   ', clm(begpatch:endpatch)%agdd0  )

  call accext  ('AGDD5   ', clm(begpatch:endpatch)%agdd5  , nstep)
  call histslf ('AGDD5   ', clm(begpatch:endpatch)%agdd5  )

  call accext  ('FNPSN10 ', clm(begpatch:endpatch)%fnpsn10, nstep)
  call accext  ('PREC365 ', clm(begpatch:endpatch)%prec365, nstep)
  call accext  ('AGDDTW  ', clm(begpatch:endpatch)%agddtw , nstep)
  call accext  ('AGDD    ', clm(begpatch:endpatch)%agdd   , nstep)

  call histslf ('T10     ', clm(begpatch:endpatch)%t10    )

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch,endpatch
     clm(k)%t_mo_min = min(clm(k)%t_mo_min, clm(k)%t_mo)
  end do

#endif

  return
end subroutine histUpdate

!=======================================================================
