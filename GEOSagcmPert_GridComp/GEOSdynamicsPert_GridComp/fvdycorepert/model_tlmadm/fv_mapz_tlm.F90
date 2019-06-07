!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
! SJL: Apr 12, 2012
! This revision may actually produce rounding level differences due to the elimination of KS to compute
! pressure level for remapping.
module fv_mapz_tlm_mod

  use constants_mod,     only: radius, pi=>pi_8, rvgas, rdgas, grav, hlv, hlf, cp_air, cp_vapor
  use tracer_manager_mod,only: get_tracer_index
  use field_manager_mod, only: MODEL_ATMOS
  use fv_grid_utils_mod, only: ptop_min
  use fv_grid_utils_tlm_mod, only: g_sum
  use fv_grid_utils_tlm_mod, only: g_sum_tlm
  use fv_fill_mod,       only: fillz
  use mpp_domains_mod,   only: mpp_update_domains, domain2d
  use mpp_mod,           only: FATAL, mpp_error, get_unit, mpp_root_pe, mpp_pe
  use fv_arrays_mod,     only: fv_grid_type, fv_flags_type
  use fv_timing_mod,     only: timing_on, timing_off
  use fv_mp_mod,         only: is_master
  use fv_cmp_mod,        only: qs_init, fv_sat_adj
  use fv_diagnostics_mod, only: prt_mxm
  use fv_arrays_nlm_mod, only: fpp

  implicit none
  real, parameter:: consv_min= 0.001   ! below which no correction applies
  real, parameter:: t_min= 184.   ! below which applies stricter constraint
  real, parameter:: r2=1./2., r0=0.0
  real, parameter:: r3 = 1./3., r23 = 2./3., r12 = 1./12.
  real, parameter:: cv_vap = 3.*rvgas  ! 1384.5
  real, parameter:: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
! real, parameter:: c_ice = 2106.           ! heat capacity of ice at 0.C
  real, parameter:: c_ice = 1972.           ! heat capacity of ice at -15.C
  real, parameter:: c_liq = 4.1855e+3    ! GFS: heat capacity of water at 0C
! real, parameter:: c_liq = 4218.        ! ECMWF-IFS
  real, parameter:: cp_vap = cp_vapor   ! 1846.
  real, parameter:: tice = 273.16

  real :: E_Flux = 0.
  private

  public compute_total_energy, Lagrangian_to_Eulerian, moist_cv, moist_cp,   &
         rst_remap, mappm, E_Flux, map1_q2
  public compute_total_energy_tlm, Lagrangian_to_Eulerian_tlm,  &
         map1_q2_tlm

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of lagrangian_to_eulerian in forward (tangent) mode:
!   variations   of useful results: peln q u v w delp ua delz omga
!                te0_2d pkz pe pk ps pt te
!   with respect to varying inputs: ws peln q u v w delp ua delz
!                omga te0_2d pkz pe pk ps pt te
  SUBROUTINE LAGRANGIAN_TO_EULERIAN_TLM(last_step, consv, ps, ps_tl, pe&
&   , pe_tl, delp, delp_tl, pkz, pkz_tl, pk, pk_tl, mdt, pdt, km, is, ie&
&   , js, je, isd, ied, jsd, jed, nq, nwat, sphum, q_con, u, u_tl, v, &
&   v_tl, w, w_tl, delz, delz_tl, pt, pt_tl, q, q_tl, hs, r_vir, cp, &
&   akap, cappa, kord_mt, kord_wz, kord_tr, kord_tm, peln, peln_tl, &
&   te0_2d, te0_2d_tl, ng, ua, ua_tl, va, omga, omga_tl, te, te_tl, ws, &
&   ws_tl, fill, reproduce_sum, out_dt, dtdt, ptop, ak, bk, pfull, &
&   flagstruct, gridstruct, domain, do_sat_adj, hydrostatic, hybrid_z, &
&   do_omega, adiabatic, do_adiabatic_init, mfx, mfy, remap_option, &
&   kord_mt_pert, kord_wz_pert, kord_tr_pert, kord_tm_pert)
    IMPLICIT NONE
!$OMP end parallel
    LOGICAL, INTENT(IN) :: last_step
! remap time step
    REAL, INTENT(IN) :: mdt
! phys time step
    REAL, INTENT(IN) :: pdt
    INTEGER, INTENT(IN) :: km
! number of tracers (including h2o)
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: nwat
! index for water vapor (specific humidity)
    INTEGER, INTENT(IN) :: sphum
    INTEGER, INTENT(IN) :: ng
! starting & ending X-Dir index
    INTEGER, INTENT(IN) :: is, ie, isd, ied
! starting & ending Y-Dir index
    INTEGER, INTENT(IN) :: js, je, jsd, jed
! Mapping order for the vector winds
    INTEGER, INTENT(IN) :: kord_mt
! Mapping order/option for w
    INTEGER, INTENT(IN) :: kord_wz
! Mapping order for tracers
    INTEGER, INTENT(IN) :: kord_tr(nq)
! Mapping order for thermodynamics
    INTEGER, INTENT(IN) :: kord_tm
! Mapping order for the vector winds
    INTEGER, INTENT(IN) :: kord_mt_pert
! Mapping order/option for w
    INTEGER, INTENT(IN) :: kord_wz_pert
! Mapping order for tracers
    INTEGER, INTENT(IN) :: kord_tr_pert(nq)
! Mapping order for thermodynamics
    INTEGER, INTENT(IN) :: kord_tm_pert
! factor for TE conservation
    REAL, INTENT(IN) :: consv
    REAL, INTENT(IN) :: r_vir
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: akap
! surface geopotential
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: te0_2d(is:ie, js:je)
    REAL, INTENT(INOUT) :: te0_2d_tl(is:ie, js:je)
    REAL, INTENT(IN) :: ws(is:ie, js:je)
    REAL, INTENT(IN) :: ws_tl(is:ie, js:je)
    LOGICAL, INTENT(IN) :: do_sat_adj
! fill negative tracers
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
    LOGICAL, INTENT(IN) :: do_omega, adiabatic, do_adiabatic_init
    REAL, INTENT(IN) :: ptop
    REAL, INTENT(IN) :: ak(km+1)
    REAL, INTENT(IN) :: bk(km+1)
    REAL, INTENT(IN) :: pfull(km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(INOUT) :: flagstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! !INPUT/OUTPUT
! pe to the kappa
    REAL, INTENT(INOUT) :: pk(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: pk_tl(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed, km, nq)
! pressure thickness
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: delp_tl(isd:ied, jsd:jed, km)
! pressure at layer edges
    REAL, INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(INOUT) :: pe_tl(is-1:ie+1, km+1, js-1:je+1)
! surface pressure
    REAL, INTENT(INOUT) :: ps(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: ps_tl(isd:ied, jsd:jed)
! u-wind will be ghosted one latitude to the north upon exit
! u-wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: u_tl(isd:ied, jsd:jed+1, km)
! v-wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_tl(isd:ied+1, jsd:jed, km)
! vertical velocity (m/s)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: w_tl(isd:ied, jsd:jed, km)
! cp*virtual potential temperature
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: pt_tl(isd:ied, jsd:jed, km)
! as input; output: temperature
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: delz
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: delz_tl
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: q_con, cappa
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: hybrid_z
    LOGICAL, INTENT(IN) :: out_dt
! u-wind (m/s) on physics grid
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: ua_tl(isd:ied, jsd:jed, km)
! v-wind (m/s) on physics grid
    REAL, INTENT(INOUT) :: va(isd:ied, jsd:jed, km)
! vertical press. velocity (pascal/sec)
    REAL, INTENT(INOUT) :: omga(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: omga_tl(isd:ied, jsd:jed, km)
! log(pe)
    REAL, INTENT(INOUT) :: peln(is:ie, km+1, js:je)
    REAL, INTENT(INOUT) :: peln_tl(is:ie, km+1, js:je)
    REAL, INTENT(INOUT) :: dtdt(is:ie, js:je, km)
! layer-mean pk for converting t to pt
    REAL, INTENT(OUT) :: pkz(is:ie, js:je, km)
    REAL, INTENT(OUT) :: pkz_tl(is:ie, js:je, km)
    REAL, INTENT(OUT) :: te(isd:ied, jsd:jed, km)
    REAL, INTENT(OUT) :: te_tl(isd:ied, jsd:jed, km)
! Mass fluxes
! X-dir Mass Flux
    REAL, OPTIONAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, km)
! Y-dir Mass Flux
    REAL, OPTIONAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, km)
! 0: remap  T in logP
    INTEGER, INTENT(IN) :: remap_option
! 1: remap PT in P
! 3: remap TE in logP with GMAO cubic
! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
    REAL, DIMENSION(is:ie, js:je) :: te_2d, zsum0, zsum1, dpln
    REAL, DIMENSION(is:ie, js:je) :: te_2d_tl, zsum0_tl, zsum1_tl
    REAL, DIMENSION(is:ie, km) :: q2, dp2
    REAL, DIMENSION(is:ie, km) :: q2_tl, dp2_tl
    REAL, DIMENSION(is:ie, km+1) :: pe1, pe2, pk1, pk2, pn2, phis
    REAL, DIMENSION(is:ie, km+1) :: pe1_tl, pe2_tl, pk1_tl, pk2_tl, &
&   pn2_tl, phis_tl
    REAL, DIMENSION(is:ie+1, km+1) :: pe0, pe3
    REAL, DIMENSION(is:ie+1, km+1) :: pe0_tl, pe3_tl
    REAL, DIMENSION(is:ie) :: gz, cvm, qv
    REAL, DIMENSION(is:ie) :: gz_tl
    REAL :: rcp, rg, tmp, tpe, rrg, bkh, dtmp, k1k, dlnp
    REAL :: tmp_tl, tpe_tl, dtmp_tl, dlnp_tl
    LOGICAL :: fast_mp_consv
    INTEGER :: i, j, k
    INTEGER :: nt, liq_wat, ice_wat, rainwat, snowwat, cld_amt, graupel&
&   , iq, n, kmp, kp, k_next
    LOGICAL :: remap_t, remap_pt, remap_te
    INTEGER :: abs_kord_tm, abs_kord_tm_pert
    INTEGER :: iep1, jep1, iedp1, jedp1
    INTRINSIC ABS
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC PRESENT
    REAL :: abs0
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: arg2
    REAL :: arg2_tl
    REAL :: result1
    REAL :: result1_tl

    REAL :: pt_tj(isd:ied, jsd:jed, km)
    REAL :: q_tj(isd:ied, jsd:jed, km, nq)
    REAL :: q2_tj(is:ie, km)
    REAL :: delz_tj(isd:ied, jsd:jed, km)
    REAL :: u_tj(isd:ied, jsd:jed+1, km)
    REAL :: v_tj(isd:ied+1, jsd:jed, km)
    REAL :: w_tj(isd:ied, jsd:jed, km)

    IF (kord_tm .GE. 0.) THEN
      abs_kord_tm = kord_tm
    ELSE
      abs_kord_tm = -kord_tm
    END IF
    IF (kord_tm_pert .GE. 0.) THEN
      abs_kord_tm_pert = kord_tm_pert
    ELSE
      abs_kord_tm_pert = -kord_tm_pert
    END IF
    iep1 = ie + 1
    jep1 = je + 1
    iedp1 = ied + 1
    jedp1 = jed + 1
    remap_t = .false.
    remap_pt = .false.
    remap_te = .false.
    SELECT CASE  (remap_option) 
    CASE (0) 
      remap_t = .true.
    CASE (1) 
      remap_pt = .true.
    CASE (2) 
      remap_te = .true.
    CASE DEFAULT
      PRINT*, ' INVALID REMAPPING OPTION '
      STOP
    END SELECT
    IF (IS_MASTER() .AND. flagstruct%fv_debug) THEN
      PRINT*, ''
      SELECT CASE  (remap_option) 
      CASE (0) 
        PRINT*, ' REMAPPING  T in logP '
      CASE (1) 
        PRINT*, ' REMAPPING PT in P'
      CASE (2) 
        PRINT*, ' REMAPPING TE in logP with GMAO cubic'
      END SELECT
      PRINT*, ' REMAPPING CONSV:     ', consv
      PRINT*, ' REMAPPING CONSV_MIN: ', consv_min
      PRINT*, ''
    END IF
    IF (flagstruct%fv_debug) CALL PRT_MXM('remap-0  PT', pt, is, ie, js&
&                                   , je, ng, km, 1., gridstruct%area_64&
&                                   , domain)
! akap / (1.-akap) = rg/Cv=0.4
    k1k = rdgas/cv_air
    rg = rdgas
    rcp = 1./cp
    rrg = -(rdgas/grav)
    IF (fpp%fpp_mapl_mode) THEN
      liq_wat = 2
      ice_wat = 3
      rainwat = -1
      snowwat = -1
      graupel = -1
      cld_amt = -1
    ELSE
      liq_wat = GET_TRACER_INDEX(model_atmos, 'liq_wat')
      ice_wat = GET_TRACER_INDEX(model_atmos, 'ice_wat')
      rainwat = GET_TRACER_INDEX(model_atmos, 'rainwat')
      snowwat = GET_TRACER_INDEX(model_atmos, 'snowwat')
      graupel = GET_TRACER_INDEX(model_atmos, 'graupel')
      cld_amt = GET_TRACER_INDEX(model_atmos, 'cld_amt')
    END IF
    IF (do_sat_adj) THEN
      fast_mp_consv = .NOT.do_adiabatic_init .AND. consv .GT. consv_min
      DO k=1,km
        kmp = k
        IF (pfull(k) .GT. 10.e2) EXIT
      END DO
      CALL QS_INIT(kmp)
      phis_tl = 0.0
      pe0_tl = 0.0
      pe1_tl = 0.0
      pe2_tl = 0.0
      pe3_tl = 0.0
      dp2_tl = 0.0
      q2_tl = 0.0
      pn2_tl = 0.0
      pk1_tl = 0.0
      pk2_tl = 0.0
    ELSE
      phis_tl = 0.0
      pe0_tl = 0.0
      pe1_tl = 0.0
      pe2_tl = 0.0
      pe3_tl = 0.0
      dp2_tl = 0.0
      q2_tl = 0.0
      pn2_tl = 0.0
      pk1_tl = 0.0
      pk2_tl = 0.0
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,km,pe,ptop,kord_tm,hydrostatic, &
!$OMP                                  pt,pk,rg,peln,q,nwat,liq_wat,rainwat,ice_wat,snowwat,    &
!$OMP                                  graupel,q_con,sphum,cappa,r_vir,rcp,k1k,delp, &
!$OMP                                  delz,akap,pkz,te,u,v,ps, gridstruct, last_step, &
!$OMP                                  ak,bk,nq,isd,ied,jsd,jed,kord_tr,fill, adiabatic, &
!$OMP                                  hs,w,ws,kord_wz,do_omega,omga,rrg,kord_mt,ua)    &
!$OMP                          private(qv,gz,cvm,kp,k_next,bkh,dp2,   &
!$OMP                                  pe0,pe1,pe2,pe3,pk1,pk2,pn2,phis,q2)
    DO j=js,je+1
      DO k=1,km+1
        DO i=is,ie
          pe1_tl(i, k) = pe_tl(i, k, j)
          pe1(i, k) = pe(i, k, j)
        END DO
      END DO
      DO i=is,ie
        pe2_tl(i, 1) = 0.0
        pe2(i, 1) = ptop
        pe2_tl(i, km+1) = pe_tl(i, km+1, j)
        pe2(i, km+1) = pe(i, km+1, j)
      END DO
!(j < je+1)
      IF (j .NE. je + 1) THEN
        IF (remap_t) THEN
! hydro test
! Remap T in logP
! Note: pt at this stage is Theta_v
          IF (hydrostatic) THEN
! Transform virtual pt to virtual Temp
            DO k=1,km
              DO i=is,ie
                pt_tl(i, j, k) = ((pt_tl(i, j, k)*(pk(i, j, k+1)-pk(i, j&
&                 , k))+pt(i, j, k)*(pk_tl(i, j, k+1)-pk_tl(i, j, k)))*&
&                 akap*(peln(i, k+1, j)-peln(i, k, j))-pt(i, j, k)*(pk(i&
&                 , j, k+1)-pk(i, j, k))*akap*(peln_tl(i, k+1, j)-&
&                 peln_tl(i, k, j)))/(akap*(peln(i, k+1, j)-peln(i, k, j&
&                 )))**2
                pt(i, j, k) = pt(i, j, k)*(pk(i, j, k+1)-pk(i, j, k))/(&
&                 akap*(peln(i, k+1, j)-peln(i, k, j)))
              END DO
            END DO
          ELSE
! Transform "density pt" to "density temp"
            DO k=1,km
              DO i=is,ie
                arg1_tl = (rrg*delp_tl(i, j, k)*delz(i, j, k)-rrg*delp(i&
&                 , j, k)*delz_tl(i, j, k))*pt(i, j, k)/delz(i, j, k)**2&
&                 + rrg*delp(i, j, k)*pt_tl(i, j, k)/delz(i, j, k)
                arg1 = rrg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
                arg2_tl = k1k*arg1_tl/arg1
                arg2 = k1k*LOG(arg1)
                pt_tl(i, j, k) = pt_tl(i, j, k)*EXP(arg2) + pt(i, j, k)*&
&                 arg2_tl*EXP(arg2)
                pt(i, j, k) = pt(i, j, k)*EXP(arg2)
              END DO
            END DO
          END IF
        ELSE IF (.NOT.remap_pt) THEN
! Using dry pressure for the definition of the virtual potential temperature
!                    pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*    &
!                                              pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
! Remap PT in P
! pt is already virtual PT
          IF (remap_te) THEN
! Remap TE in logP
! Transform virtual pt to total energy
            CALL PKEZ_TLM(km, is, ie, js, je, j, pe, pk, pk_tl, akap, &
&                   peln, peln_tl, pkz, pkz_tl, ptop)
! Compute cp*T + KE
            DO k=1,km
              DO i=is,ie
                te_tl(i, j, k) = 0.25*gridstruct%rsin2(i, j)*(2*u(i, j, &
&                 k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i&
&                 , j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl(i+1, j, k)-&
&                 gridstruct%cosa_s(i, j)*((u_tl(i, j, k)+u_tl(i, j+1, k&
&                 ))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i, j+1, k))&
&                 *(v_tl(i, j, k)+v_tl(i+1, j, k)))) + cp_air*(pt_tl(i, &
&                 j, k)*pkz(i, j, k)+pt(i, j, k)*pkz_tl(i, j, k))
                te(i, j, k) = 0.25*gridstruct%rsin2(i, j)*(u(i, j, k)**2&
&                 +u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j&
&                 , k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&                 gridstruct%cosa_s(i, j)) + cp_air*pt(i, j, k)*pkz(i, j&
&                 , k)
              END DO
            END DO
          END IF
        END IF
        IF (.NOT.hydrostatic) THEN
          DO k=1,km
            DO i=is,ie
! ="specific volume"/grav
              delz_tl(i, j, k) = -((delz_tl(i, j, k)*delp(i, j, k)-delz(&
&               i, j, k)*delp_tl(i, j, k))/delp(i, j, k)**2)
              delz(i, j, k) = -(delz(i, j, k)/delp(i, j, k))
            END DO
          END DO
        END IF
! update ps
        DO i=is,ie
          ps_tl(i, j) = pe1_tl(i, km+1)
          ps(i, j) = pe1(i, km+1)
        END DO
!
! Hybrid sigma-P coordinate:
!
        DO k=2,km
          DO i=is,ie
            pe2_tl(i, k) = bk(k)*pe_tl(i, km+1, j)
            pe2(i, k) = ak(k) + bk(k)*pe(i, km+1, j)
          END DO
        END DO
        DO k=1,km
          DO i=is,ie
            dp2_tl(i, k) = pe2_tl(i, k+1) - pe2_tl(i, k)
            dp2(i, k) = pe2(i, k+1) - pe2(i, k)
          END DO
        END DO
!------------
! update delp
!------------
        DO k=1,km
          DO i=is,ie
            delp_tl(i, j, k) = dp2_tl(i, k)
            delp(i, j, k) = dp2(i, k)
          END DO
        END DO
!------------------
! Compute p**Kappa
!------------------
        DO k=1,km+1
          DO i=is,ie
            pk1_tl(i, k) = pk_tl(i, j, k)
            pk1(i, k) = pk(i, j, k)
          END DO
        END DO
        DO i=is,ie
          pn2_tl(i, 1) = peln_tl(i, 1, j)
          pn2(i, 1) = peln(i, 1, j)
          pn2_tl(i, km+1) = peln_tl(i, km+1, j)
          pn2(i, km+1) = peln(i, km+1, j)
          pk2_tl(i, 1) = pk1_tl(i, 1)
          pk2(i, 1) = pk1(i, 1)
          pk2_tl(i, km+1) = pk1_tl(i, km+1)
          pk2(i, km+1) = pk1(i, km+1)
        END DO
        DO k=2,km
          DO i=is,ie
            pn2_tl(i, k) = pe2_tl(i, k)/pe2(i, k)
            pn2(i, k) = LOG(pe2(i, k))
            pk2_tl(i, k) = akap*pn2_tl(i, k)*EXP(akap*pn2(i, k))
            pk2(i, k) = EXP(akap*pn2(i, k))
          END DO
        END DO
        IF (remap_t) THEN
!----------------------------------
! Map t using logp
!----------------------------------
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            CALL MAP_SCALAR_TLM(km, peln(is:ie, 1:km+1, j), peln_tl(&
&                            is:ie, 1:km+1, j), gz, km, pn2, pn2_tl, pt&
&                            , pt_tl, is, ie, j, isd, ied, jsd, jed, 1, &
&                            abs_kord_tm, t_min)
          ELSE
pt_tj = pt
            CALL MAP_SCALAR_TLM(km, peln(is:ie, 1:km+1, j), peln_tl(is:&
&                         ie, 1:km+1, j), gz, km, pn2, pn2_tl, pt_tj, pt_tl&
&                         , is, ie, j, isd, ied, jsd, jed, 1, &
&                         abs_kord_tm_pert, t_min)
      call map_scalar(km,  peln(is:ie,1:km+1,j),      gz,   &
                      km,  pn2,           pt,              &
                      is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm, t_min)
          END IF
        ELSE IF (remap_pt) THEN
!----------------------------------
! Map pt using pe
!----------------------------------
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            gz_tl = 0.0
            CALL MAP1_PPM_TLM(km, pe1, pe1_tl, gz, gz_tl, km, pe2, &
&                          pe2_tl, pt, pt_tl, is, ie, j, isd, ied, jsd, &
&                          jed, 1, abs_kord_tm)
          ELSE
            gz_tl = 0.0
pt_tj = pt
            CALL MAP1_PPM_TLM(km, pe1, pe1_tl, gz, gz_tl, km, pe2, &
&                       pe2_tl, pt_tj, pt_tl, is, ie, j, isd, ied, jsd, jed&
&                       , 1, abs_kord_tm_pert)
      call map1_ppm (km,  pe1,       gz,       &
                     km,  pe2,  pt,                  &
                     is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm)
          END IF
        ELSE IF (remap_te) THEN
!----------------------------------
! map Total Energy using GMAO cubic
!----------------------------------
          DO i=is,ie
            phis_tl(i, km+1) = 0.0
            phis(i, km+1) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              phis_tl(i, k) = phis_tl(i, k+1) + cp_air*(pt_tl(i, j, k)*(&
&               pk1(i, k+1)-pk1(i, k))+pt(i, j, k)*(pk1_tl(i, k+1)-&
&               pk1_tl(i, k)))
              phis(i, k) = phis(i, k+1) + cp_air*pt(i, j, k)*(pk1(i, k+1&
&               )-pk1(i, k))
            END DO
          END DO
          DO k=1,km+1
            DO i=is,ie
              phis_tl(i, k) = phis_tl(i, k)*pe1(i, k) + phis(i, k)*&
&               pe1_tl(i, k)
              phis(i, k) = phis(i, k)*pe1(i, k)
            END DO
          END DO
          DO k=1,km
            DO i=is,ie
              te_tl(i, j, k) = te_tl(i, j, k) + ((phis_tl(i, k+1)-&
&               phis_tl(i, k))*(pe1(i, k+1)-pe1(i, k))-(phis(i, k+1)-&
&               phis(i, k))*(pe1_tl(i, k+1)-pe1_tl(i, k)))/(pe1(i, k+1)-&
&               pe1(i, k))**2
              te(i, j, k) = te(i, j, k) + (phis(i, k+1)-phis(i, k))/(pe1&
&               (i, k+1)-pe1(i, k))
            END DO
          END DO
! Map te using log P in GMAO cubic
          CALL MAP1_CUBIC_TLM(km, pe1, pe1_tl, km, pe2, pe2_tl, te, &
&                       te_tl, is, ie, j, isd, ied, jsd, jed, akap, &
&                       t_var=1, conserv=.true.)
        END IF
!----------------
! Map constituents
!----------------
        IF (nq .GT. 5) THEN
          IF (kord_tr(1) .EQ. kord_tr_pert(1)) THEN
            CALL MAPN_TRACER_TLM(nq, km, pe1, pe1_tl, pe2, pe2_tl, q&
&                             , q_tl, dp2, dp2_tl, kord_tr, j, is, ie, &
&                             isd, ied, jsd, jed, 0., fill)
          ELSE
q_tj = q
            CALL MAPN_TRACER_TLM(nq, km, pe1, pe1_tl, pe2, pe2_tl, q_tj, &
&                          q_tl, dp2, dp2_tl, kord_tr_pert, j, is, ie, &
&                          isd, ied, jsd, jed, 0., fill)
      call mapn_tracer(nq, km, pe1, pe2, q, dp2, kord_tr, j,     &
                       is, ie, isd, ied, jsd, jed, 0., fill)
          END IF
        ELSE IF (nq .GT. 0) THEN
! Remap one tracer at a time
          DO iq=1,nq
            IF (kord_tr(iq) .EQ. kord_tr_pert(iq)) THEN
              CALL MAP1_Q2_TLM(km, pe1, pe1_tl, q(isd:ied, jsd:jed, 1&
&                           :km, iq), q_tl(isd:ied, jsd:jed, 1:km, iq), &
&                           km, pe2, pe2_tl, q2, q2_tl, dp2, dp2_tl, is&
&                           , ie, 0, kord_tr(iq), j, isd, ied, jsd, jed&
&                           , 0.)
            ELSE
q2_tj = q2
              CALL MAP1_Q2_TLM(km, pe1, pe1_tl, q(isd:ied, jsd:jed, 1:km&
&                        , iq), q_tl(isd:ied, jsd:jed, 1:km, iq), km, &
&                        pe2, pe2_tl, q2_tj, q2_tl, dp2, dp2_tl, is, ie, 0&
&                        , kord_tr_pert(iq), j, isd, ied, jsd, jed, 0.)
      call map1_q2(km, pe1, q(isd:ied,jsd:jed,1:km,iq),     &
                   km, pe2, q2, dp2,             &
                   is, ie, 0, kord_tr(iq), j, isd, ied, jsd, jed, 0.)
            END IF
            IF (fill) CALL FILLZ(ie - is + 1, km, 1, q2, dp2)
            DO k=1,km
              DO i=is,ie
                q_tl(i, j, k, iq) = q2_tl(i, k)
                q(i, j, k, iq) = q2(i, k)
              END DO
            END DO
          END DO
        END IF
        IF (.NOT.hydrostatic) THEN
! Remap vertical wind:
          IF (kord_wz .EQ. kord_wz_pert) THEN
            CALL MAP1_PPM_TLM(km, pe1, pe1_tl, ws(is:ie, j), ws_tl(is&
&                          :ie, j), km, pe2, pe2_tl, w, w_tl, is, ie, j&
&                          , isd, ied, jsd, jed, -2, kord_wz)
          ELSE
w_tj = w
            CALL MAP1_PPM_TLM(km, pe1, pe1_tl, ws(is:ie, j), ws_tl(is:ie&
&                       , j), km, pe2, pe2_tl, w_tj, w_tl, is, ie, j, isd, &
&                       ied, jsd, jed, -2, kord_wz_pert)
      call map1_ppm (km,   pe1,      ws(is:ie,j),   &
                     km,   pe2,  w,              &
                     is, ie, j, isd, ied, jsd, jed, -2, kord_wz)
          END IF
! Remap delz for hybrid sigma-p coordinate
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            gz_tl = 0.0
            CALL MAP1_PPM_TLM(km, pe1, pe1_tl, gz, gz_tl, km, pe2, &
&                          pe2_tl, delz, delz_tl, is, ie, j, isd, ied, &
&                          jsd, jed, 1, abs_kord_tm)
          ELSE
            gz_tl = 0.0
delz_tj = delz
            CALL MAP1_PPM_TLM(km, pe1, pe1_tl, gz, gz_tl, km, pe2, &
&                       pe2_tl, delz_tj, delz_tl, is, ie, j, isd, ied, jsd&
&                       , jed, 1, abs_kord_tm_pert)
      call map1_ppm (km,   pe1,        gz,   &
                     km,   pe2, delz,              &
                     is, ie, j, isd,  ied,  jsd,  jed,  1, abs_kord_tm)
          END IF
          DO k=1,km
            DO i=is,ie
              delz_tl(i, j, k) = -(delz_tl(i, j, k)*dp2(i, k)+delz(i, j&
&               , k)*dp2_tl(i, k))
              delz(i, j, k) = -(delz(i, j, k)*dp2(i, k))
            END DO
          END DO
        END IF
!----------
! Update pk
!----------
        DO k=1,km+1
          DO i=is,ie
            pk_tl(i, j, k) = pk2_tl(i, k)
            pk(i, j, k) = pk2(i, k)
          END DO
        END DO
!----------------
        IF (do_omega) THEN
! Start do_omega
! Copy omega field to pe3
          DO i=is,ie
            pe3_tl(i, 1) = 0.0
            pe3(i, 1) = 0.
          END DO
          DO k=2,km+1
            DO i=is,ie
              pe3_tl(i, k) = omga_tl(i, j, k-1)
              pe3(i, k) = omga(i, j, k-1)
            END DO
          END DO
        END IF
        DO k=1,km+1
          DO i=is,ie
            pe0_tl(i, k) = peln_tl(i, k, j)
            pe0(i, k) = peln(i, k, j)
            peln_tl(i, k, j) = pn2_tl(i, k)
            peln(i, k, j) = pn2(i, k)
          END DO
        END DO
!------------
! Compute pkz
!------------
        IF (hydrostatic) THEN
          DO k=1,km
            DO i=is,ie
              pkz_tl(i, j, k) = ((pk2_tl(i, k+1)-pk2_tl(i, k))*akap*(&
&               peln(i, k+1, j)-peln(i, k, j))-(pk2(i, k+1)-pk2(i, k))*&
&               akap*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(akap*(peln(&
&               i, k+1, j)-peln(i, k, j)))**2
              pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1&
&               , j)-peln(i, k, j)))
            END DO
          END DO
        ELSE IF (remap_te) THEN
! WMP: note that this is where TE remapping non-hydrostatic is invalid and cannot be run
          GOTO 120
        ELSE IF (remap_t) THEN
! Note: pt at this stage is T_v or T_m
          DO k=1,km
            DO i=is,ie
              arg1_tl = (rrg*delp_tl(i, j, k)*delz(i, j, k)-rrg*delp(i, &
&               j, k)*delz_tl(i, j, k))*pt(i, j, k)/delz(i, j, k)**2 + &
&               rrg*delp(i, j, k)*pt_tl(i, j, k)/delz(i, j, k)
              arg1 = rrg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
              arg2_tl = akap*arg1_tl/arg1
              arg2 = akap*LOG(arg1)
              pkz_tl(i, j, k) = arg2_tl*EXP(arg2)
              pkz(i, j, k) = EXP(arg2)
            END DO
          END DO
        ELSE
! Using dry pressure for the definition of the virtual potential temperature
!           pkz(i,j,k) = exp(akap*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
! Note: pt at this stage is Theta_v
          DO k=1,km
            DO i=is,ie
              arg1_tl = (rrg*delp_tl(i, j, k)*delz(i, j, k)-rrg*delp(i, &
&               j, k)*delz_tl(i, j, k))*pt(i, j, k)/delz(i, j, k)**2 + &
&               rrg*delp(i, j, k)*pt_tl(i, j, k)/delz(i, j, k)
              arg1 = rrg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
              arg2_tl = k1k*arg1_tl/arg1
              arg2 = k1k*LOG(arg1)
              pkz_tl(i, j, k) = arg2_tl*EXP(arg2)
              pkz(i, j, k) = EXP(arg2)
            END DO
          END DO
        END IF
! end do_omega
! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
        IF (do_omega) THEN
          DO k=1,km
            DO i=is,ie
              dp2_tl(i, k) = 0.5*(peln_tl(i, k, j)+peln_tl(i, k+1, j))
              dp2(i, k) = 0.5*(peln(i, k, j)+peln(i, k+1, j))
            END DO
          END DO
          DO i=is,ie
            k_next = 1
            DO 110 n=1,km
              kp = k_next
              DO k=kp,km
                IF (dp2(i, n) .LE. pe0(i, k+1) .AND. dp2(i, n) .GE. pe0(&
&                   i, k)) GOTO 100
              END DO
              GOTO 110
 100          omga_tl(i, j, n) = pe3_tl(i, k) + (((pe3_tl(i, k+1)-pe3_tl&
&               (i, k))*(dp2(i, n)-pe0(i, k))+(pe3(i, k+1)-pe3(i, k))*(&
&               dp2_tl(i, n)-pe0_tl(i, k)))*(pe0(i, k+1)-pe0(i, k))-(pe3&
&               (i, k+1)-pe3(i, k))*(dp2(i, n)-pe0(i, k))*(pe0_tl(i, k+1&
&               )-pe0_tl(i, k)))/(pe0(i, k+1)-pe0(i, k))**2
              omga(i, j, n) = pe3(i, k) + (pe3(i, k+1)-pe3(i, k))*(dp2(i&
&               , n)-pe0(i, k))/(pe0(i, k+1)-pe0(i, k))
              k_next = k
 110        CONTINUE
          END DO
        END IF
      END IF
      DO i=is,ie+1
        pe0_tl(i, 1) = pe_tl(i, 1, j)
        pe0(i, 1) = pe(i, 1, j)
      END DO
!------
! map u
!------
      DO k=2,km+1
        DO i=is,ie
          pe0_tl(i, k) = 0.5*(pe_tl(i, k, j-1)+pe1_tl(i, k))
          pe0(i, k) = 0.5*(pe(i, k, j-1)+pe1(i, k))
        END DO
      END DO
      DO k=1,km+1
        bkh = 0.5*bk(k)
        DO i=is,ie
          pe3_tl(i, k) = bkh*(pe_tl(i, km+1, j-1)+pe1_tl(i, km+1))
          pe3(i, k) = ak(k) + bkh*(pe(i, km+1, j-1)+pe1(i, km+1))
        END DO
      END DO
      IF (kord_mt .EQ. kord_mt_pert) THEN
        gz_tl = 0.0
        CALL MAP1_PPM_TLM(km, pe0(is:ie, :), pe0_tl(is:ie, :), gz, &
&                      gz_tl, km, pe3(is:ie, :), pe3_tl(is:ie, :), u, &
&                      u_tl, is, ie, j, isd, ied, jsd, jedp1, -1, &
&                      kord_mt)
      ELSE
        gz_tl = 0.0
u_tj = u
        CALL MAP1_PPM_TLM(km, pe0(is:ie, :), pe0_tl(is:ie, :), gz, gz_tl&
&                   , km, pe3(is:ie, :), pe3_tl(is:ie, :), u_tj, u_tl, is, &
&                   ie, j, isd, ied, jsd, jedp1, -1, kord_mt_pert)
      call map1_ppm( km, pe0(is:ie,:),        gz,   &
                     km, pe3(is:ie,:),   u,               &
                     is, ie, j, isd, ied, jsd, jedp1, -1, kord_mt)
      END IF
      IF (PRESENT(mfy)) CALL MAP1_PPM(km, pe0(is:ie, :), gz, km, pe3(is:&
&                               ie, :), mfy, is, ie, j, is, ie, js, jep1&
&                               , -1, kord_mt)
! (j < je+1)
      IF (j .LT. je + 1) THEN
!------
! map v
!------
        DO i=is,ie+1
          pe3_tl(i, 1) = 0.0
          pe3(i, 1) = ak(1)
        END DO
        DO k=2,km+1
          bkh = 0.5*bk(k)
          DO i=is,ie+1
            pe0_tl(i, k) = 0.5*(pe_tl(i-1, k, j)+pe_tl(i, k, j))
            pe0(i, k) = 0.5*(pe(i-1, k, j)+pe(i, k, j))
            pe3_tl(i, k) = bkh*(pe_tl(i-1, km+1, j)+pe_tl(i, km+1, j))
            pe3(i, k) = ak(k) + bkh*(pe(i-1, km+1, j)+pe(i, km+1, j))
          END DO
        END DO
        IF (kord_mt .EQ. kord_mt_pert) THEN
          gz_tl = 0.0
          CALL MAP1_PPM_TLM(km, pe0, pe0_tl, gz, gz_tl, km, pe3, &
&                        pe3_tl, v, v_tl, is, iep1, j, isd, iedp1, jsd, &
&                        jed, -1, kord_mt)
        ELSE
          gz_tl = 0.0
v_tj = v
          CALL MAP1_PPM_TLM(km, pe0, pe0_tl, gz, gz_tl, km, pe3, pe3_tl&
&                     , v_tj, v_tl, is, iep1, j, isd, iedp1, jsd, jed, -1, &
&                     kord_mt_pert)
      call map1_ppm (km, pe0,     gz,    &
                     km, pe3,  v, is, ie+1,    &
                     j, isd, iedp1, jsd, jed, -1, kord_mt)
        END IF
        IF (PRESENT(mfx)) CALL MAP1_PPM(km, pe0, gz, km, pe3, mfx, is, &
&                                 iep1, j, is, iep1, js, je, -1, kord_mt&
&                                )
      END IF
      DO k=1,km
        DO i=is,ie
          ua_tl(i, j, k) = pe2_tl(i, k+1)
          ua(i, j, k) = pe2(i, k+1)
        END DO
      END DO
    END DO
!$OMP parallel default(none) shared(is,ie,js,je,km,kmp,ptop,u,v,pe,ua,isd,ied,jsd,jed,kord_mt, &
!$OMP                               te_2d,te,delp,hydrostatic,hs,rg,pt,peln, adiabatic, &
!$OMP                               cp,delz,nwat,rainwat,liq_wat,ice_wat,snowwat,       &
!$OMP                               graupel,q_con,r_vir,sphum,w,pk,pkz,last_step,consv, &
!$OMP                               do_adiabatic_init,zsum1,zsum0,te0_2d,domain,        &
!$OMP                               ng,gridstruct,E_Flux,pdt,dtmp,reproduce_sum,q,      &
!$OMP                               mdt,cld_amt,cappa,dtdt,out_dt,rrg,akap,do_sat_adj,  &
!$OMP                               fast_mp_consv,kord_tm) &
!$OMP                       private(pe0,pe1,pe2,pe3,qv,cvm,gz,phis,tpe,tmp, dpln)
!$OMP do
    DO k=2,km
      DO j=js,je
        DO i=is,ie
          pe_tl(i, k, j) = ua_tl(i, j, k-1)
          pe(i, k, j) = ua(i, j, k-1)
        END DO
      END DO
    END DO
    IF (flagstruct%fv_debug) THEN
      IF (kord_tm .LT. 0) THEN
        CALL PRT_MXM('remap-1  TV', pt, is, ie, js, je, ng, km, 1., &
&              gridstruct%area_64, domain)
      ELSE
        CALL PRT_MXM('remap-1  PT', pt, is, ie, js, je, ng, km, 1., &
&              gridstruct%area_64, domain)
      END IF
    END IF
    dtmp = 0.
! end last_step check
    IF (last_step .AND. (.NOT.do_adiabatic_init)) THEN
! end consv check
      IF (consv .GT. consv_min) THEN
        gz_tl = 0.0
        te_2d_tl = 0.0
        zsum0_tl = 0.0
        zsum1_tl = 0.0
!$OMP do
        DO j=js,je
          IF (remap_t) THEN
! end non-hydro
            IF (hydrostatic) THEN
              DO i=is,ie
                gz_tl(i) = 0.0
                gz(i) = hs(i, j)
                DO k=1,km
                  gz_tl(i) = gz_tl(i) + rg*(pt_tl(i, j, k)*(peln(i, k+1&
&                   , j)-peln(i, k, j))+pt(i, j, k)*(peln_tl(i, k+1, j)-&
&                   peln_tl(i, k, j)))
                  gz(i) = gz(i) + rg*pt(i, j, k)*(peln(i, k+1, j)-peln(i&
&                   , k, j))
                END DO
              END DO
              DO i=is,ie
                te_2d_tl(i, j) = hs(i, j)*pe_tl(i, km+1, j) - pe_tl(i, 1&
&                 , j)*gz(i) - pe(i, 1, j)*gz_tl(i)
                te_2d(i, j) = pe(i, km+1, j)*hs(i, j) - pe(i, 1, j)*gz(i&
&                 )
              END DO
              DO k=1,km
                DO i=is,ie
                  te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(cp&
&                   *pt(i, j, k)+0.25*gridstruct%rsin2(i, j)*(u(i, j, k)&
&                   **2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u&
&                   (i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&                   gridstruct%cosa_s(i, j))) + delp(i, j, k)*(cp*pt_tl(&
&                   i, j, k)+0.25*gridstruct%rsin2(i, j)*(2*u(i, j, k)*&
&                   u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i, &
&                   j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl(i+1, j, k)-&
&                   gridstruct%cosa_s(i, j)*((u_tl(i, j, k)+u_tl(i, j+1&
&                   , k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i, j+1&
&                   , k))*(v_tl(i, j, k)+v_tl(i+1, j, k)))))
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cp*pt(i, j&
&                   , k)+0.25*gridstruct%rsin2(i, j)*(u(i, j, k)**2+u(i&
&                   , j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, &
&                   k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&                   gridstruct%cosa_s(i, j)))
                END DO
              END DO
            ELSE
              DO i=is,ie
                te_2d_tl(i, j) = 0.0
                te_2d(i, j) = 0.
                phis_tl(i, km+1) = 0.0
                phis(i, km+1) = hs(i, j)
              END DO
              DO k=km,1,-1
                DO i=is,ie
                  phis_tl(i, k) = phis_tl(i, k+1) - grav*delz_tl(i, j, k&
&                   )
                  phis(i, k) = phis(i, k+1) - grav*delz(i, j, k)
                END DO
              END DO
              DO k=1,km
                DO i=is,ie
                  te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(&
&                   cv_air*pt(i, j, k)/(1.+r_vir*q(i, j, k, sphum))+0.5*&
&                   (phis(i, k)+phis(i, k+1)+w(i, j, k)**2+0.5*&
&                   gridstruct%rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**&
&                   2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1&
&                   , k))*(v(i, j, k)+v(i+1, j, k))*gridstruct%cosa_s(i&
&                   , j)))) + delp(i, j, k)*((cv_air*pt_tl(i, j, k)*(1.+&
&                   r_vir*q(i, j, k, sphum))-cv_air*pt(i, j, k)*r_vir*&
&                   q_tl(i, j, k, sphum))/(1.+r_vir*q(i, j, k, sphum))**&
&                   2+0.5*(phis_tl(i, k)+phis_tl(i, k+1)+2*w(i, j, k)*&
&                   w_tl(i, j, k)+0.5*gridstruct%rsin2(i, j)*(2*u(i, j, &
&                   k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(&
&                   i, j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl(i+1, j, k&
&                   )-gridstruct%cosa_s(i, j)*((u_tl(i, j, k)+u_tl(i, j+&
&                   1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i, j+&
&                   1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k))))))
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv_air*pt(i&
&                   , j, k)/(1.+r_vir*q(i, j, k, sphum))+0.5*(phis(i, k)&
&                   +phis(i, k+1)+w(i, j, k)**2+0.5*gridstruct%rsin2(i, &
&                   j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+&
&                   1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(&
&                   i+1, j, k))*gridstruct%cosa_s(i, j))))
                END DO
              END DO
            END IF
          ELSE IF (remap_pt) THEN
! k-loop
            IF (hydrostatic) THEN
              DO i=is,ie
                gz_tl(i) = 0.0
                gz(i) = hs(i, j)
                DO k=1,km
                  gz_tl(i) = gz_tl(i) + cp_air*(pt_tl(i, j, k)*(pk(i, j&
&                   , k+1)-pk(i, j, k))+pt(i, j, k)*(pk_tl(i, j, k+1)-&
&                   pk_tl(i, j, k)))
                  gz(i) = gz(i) + cp_air*pt(i, j, k)*(pk(i, j, k+1)-pk(i&
&                   , j, k))
                END DO
              END DO
              DO i=is,ie
                te_2d_tl(i, j) = hs(i, j)*pe_tl(i, km+1, j) - pe_tl(i, 1&
&                 , j)*gz(i) - pe(i, 1, j)*gz_tl(i)
                te_2d(i, j) = pe(i, km+1, j)*hs(i, j) - pe(i, 1, j)*gz(i&
&                 )
              END DO
              DO k=1,km
                DO i=is,ie
                  te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(&
&                   cp_air*pt(i, j, k)*pkz(i, j, k)+0.25*gridstruct%&
&                   rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k&
&                   )**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i&
&                   , j, k)+v(i+1, j, k))*gridstruct%cosa_s(i, j))) + &
&                   delp(i, j, k)*(cp_air*(pt_tl(i, j, k)*pkz(i, j, k)+&
&                   pt(i, j, k)*pkz_tl(i, j, k))+0.25*gridstruct%rsin2(i&
&                   , j)*(2*u(i, j, k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl&
&                   (i, j+1, k)+2*v(i, j, k)*v_tl(i, j, k)+2*v(i+1, j, k&
&                   )*v_tl(i+1, j, k)-gridstruct%cosa_s(i, j)*((u_tl(i, &
&                   j, k)+u_tl(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(&
&                   i, j, k)+u(i, j+1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k&
&                   )))))
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cp_air*pt(i&
&                   , j, k)*pkz(i, j, k)+0.25*gridstruct%rsin2(i, j)*(u(&
&                   i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, &
&                   k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j&
&                   , k))*gridstruct%cosa_s(i, j)))
                END DO
              END DO
            ELSE
!-----------------
! Non-hydrostatic:
!-----------------
              DO i=is,ie
                phis_tl(i, km+1) = 0.0
                phis(i, km+1) = hs(i, j)
                DO k=km,1,-1
                  phis_tl(i, k) = phis_tl(i, k+1) - grav*delz_tl(i, j, k&
&                   )
                  phis(i, k) = phis(i, k+1) - grav*delz(i, j, k)
                END DO
              END DO
              DO i=is,ie
                te_2d_tl(i, j) = 0.0
                te_2d(i, j) = 0.
              END DO
              DO k=1,km
                DO i=is,ie
                  te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(&
&                   cv_air*pt(i, j, k)/(1.+r_vir*q(i, j, k, sphum))+0.5*&
&                   (phis(i, k)+phis(i, k+1)+w(i, j, k)**2+0.5*&
&                   gridstruct%rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**&
&                   2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1&
&                   , k))*(v(i, j, k)+v(i+1, j, k))*gridstruct%cosa_s(i&
&                   , j)))) + delp(i, j, k)*((cv_air*pt_tl(i, j, k)*(1.+&
&                   r_vir*q(i, j, k, sphum))-cv_air*pt(i, j, k)*r_vir*&
&                   q_tl(i, j, k, sphum))/(1.+r_vir*q(i, j, k, sphum))**&
&                   2+0.5*(phis_tl(i, k)+phis_tl(i, k+1)+2*w(i, j, k)*&
&                   w_tl(i, j, k)+0.5*gridstruct%rsin2(i, j)*(2*u(i, j, &
&                   k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(&
&                   i, j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl(i+1, j, k&
&                   )-gridstruct%cosa_s(i, j)*((u_tl(i, j, k)+u_tl(i, j+&
&                   1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i, j+&
&                   1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k))))))
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv_air*pt(i&
&                   , j, k)/(1.+r_vir*q(i, j, k, sphum))+0.5*(phis(i, k)&
&                   +phis(i, k+1)+w(i, j, k)**2+0.5*gridstruct%rsin2(i, &
&                   j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+&
&                   1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(&
&                   i+1, j, k))*gridstruct%cosa_s(i, j))))
                END DO
              END DO
            END IF
          ELSE IF (remap_te) THEN
            DO i=is,ie
              te_2d_tl(i, j) = te_tl(i, j, 1)*delp(i, j, 1) + te(i, j, 1&
&               )*delp_tl(i, j, 1)
              te_2d(i, j) = te(i, j, 1)*delp(i, j, 1)
            END DO
            DO k=2,km
              DO i=is,ie
                te_2d_tl(i, j) = te_2d_tl(i, j) + te_tl(i, j, k)*delp(i&
&                 , j, k) + te(i, j, k)*delp_tl(i, j, k)
                te_2d(i, j) = te_2d(i, j) + te(i, j, k)*delp(i, j, k)
              END DO
            END DO
          END IF
          DO i=is,ie
            te_2d_tl(i, j) = te0_2d_tl(i, j) - te_2d_tl(i, j)
            te_2d(i, j) = te0_2d(i, j) - te_2d(i, j)
            zsum1_tl(i, j) = pkz_tl(i, j, 1)*delp(i, j, 1) + pkz(i, j, 1&
&             )*delp_tl(i, j, 1)
            zsum1(i, j) = pkz(i, j, 1)*delp(i, j, 1)
          END DO
          DO k=2,km
            DO i=is,ie
              zsum1_tl(i, j) = zsum1_tl(i, j) + pkz_tl(i, j, k)*delp(i, &
&               j, k) + pkz(i, j, k)*delp_tl(i, j, k)
              zsum1(i, j) = zsum1(i, j) + pkz(i, j, k)*delp(i, j, k)
            END DO
          END DO
          IF (hydrostatic) THEN
            DO i=is,ie
              zsum0_tl(i, j) = ptop*(pk_tl(i, j, 1)-pk_tl(i, j, km+1)) +&
&               zsum1_tl(i, j)
              zsum0(i, j) = ptop*(pk(i, j, 1)-pk(i, j, km+1)) + zsum1(i&
&               , j)
            END DO
          END IF
        END DO
! j-loop
!$OMP single
        result1_tl = G_SUM_TLM(domain, te_2d, te_2d_tl, is, ie, js, je, &
&         ng, gridstruct%area_64, 0, reproduce=.true., g_sum=result1)
        tpe_tl = consv*result1_tl
        tpe = consv*result1
! unit: W/m**2
        e_flux = tpe/(grav*pdt*4.*pi*radius**2)
! Note pdt is "phys" time step
        IF (hydrostatic) THEN
          result1_tl = G_SUM_TLM(domain, zsum0, zsum0_tl, is, ie, js, je&
&           , ng, gridstruct%area_64, 0, reproduce=.true., g_sum=result1&
&           )
          dtmp_tl = (tpe_tl*cp*result1-tpe*cp*result1_tl)/(cp*result1)**&
&           2
          dtmp = tpe/(cp*result1)
        ELSE
          result1_tl = G_SUM_TLM(domain, zsum1, zsum1_tl, is, ie, js, je&
&           , ng, gridstruct%area_64, 0, reproduce=.true., g_sum=result1&
&           )
          dtmp_tl = (tpe_tl*cv_air*result1-tpe*cv_air*result1_tl)/(&
&           cv_air*result1)**2
          dtmp = tpe/(cv_air*result1)
        END IF
      ELSE IF (consv .LT. -consv_min) THEN
!$OMP end single
        zsum0_tl = 0.0
        zsum1_tl = 0.0
!$OMP do
        DO j=js,je
          DO i=is,ie
            zsum1_tl(i, j) = pkz_tl(i, j, 1)*delp(i, j, 1) + pkz(i, j, 1&
&             )*delp_tl(i, j, 1)
            zsum1(i, j) = pkz(i, j, 1)*delp(i, j, 1)
          END DO
          DO k=2,km
            DO i=is,ie
              zsum1_tl(i, j) = zsum1_tl(i, j) + pkz_tl(i, j, k)*delp(i, &
&               j, k) + pkz(i, j, k)*delp_tl(i, j, k)
              zsum1(i, j) = zsum1(i, j) + pkz(i, j, k)*delp(i, j, k)
            END DO
          END DO
          IF (hydrostatic) THEN
            DO i=is,ie
              zsum0_tl(i, j) = ptop*(pk_tl(i, j, 1)-pk_tl(i, j, km+1)) +&
&               zsum1_tl(i, j)
              zsum0(i, j) = ptop*(pk(i, j, 1)-pk(i, j, km+1)) + zsum1(i&
&               , j)
            END DO
          END IF
        END DO
        e_flux = consv
!$OMP single
        IF (hydrostatic) THEN
          result1_tl = G_SUM_TLM(domain, zsum0, zsum0_tl, is, ie, js, je&
&           , ng, gridstruct%area_64, 0, reproduce=.true., g_sum=result1&
&           )
          dtmp_tl = -(e_flux*grav*pdt*4.*pi*radius**2*cp*result1_tl/(cp*&
&           result1)**2)
          dtmp = e_flux*(grav*pdt*4.*pi*radius**2)/(cp*result1)
          gz_tl = 0.0
        ELSE
          result1_tl = G_SUM_TLM(domain, zsum1, zsum1_tl, is, ie, js, je&
&           , ng, gridstruct%area_64, 0, reproduce=.true., g_sum=result1&
&           )
          dtmp_tl = -(e_flux*grav*pdt*4.*pi*radius**2*cv_air*result1_tl/&
&           (cv_air*result1)**2)
          dtmp = e_flux*(grav*pdt*4.*pi*radius**2)/(cv_air*result1)
          gz_tl = 0.0
        END IF
      ELSE
        gz_tl = 0.0
        dtmp_tl = 0.0
      END IF
    ELSE
      gz_tl = 0.0
      dtmp_tl = 0.0
    END IF
!$OMP end single
! do_sat_adj
! Note: pt at this stage is T_v
    IF (remap_t .AND. (.NOT.do_adiabatic_init) .AND. do_sat_adj) THEN
! if ( do_sat_adj ) then
      CALL TIMING_ON('sat_adj2')
!$OMP do
      DO k=kmp,km
        DO j=js,je
          DO i=is,ie
            dpln(i, j) = peln(i, k+1, j) - peln(i, k, j)
          END DO
        END DO
        IF (mdt .GE. 0.) THEN
          abs0 = mdt
        ELSE
          abs0 = -mdt
        END IF
        CALL FV_SAT_ADJ(abs0, r_vir, is, ie, js, je, ng, hydrostatic, &
&                 fast_mp_consv, te(isd:ied, jsd:jed, k), q(isd:ied, jsd&
&                 :jed, k, sphum), q(isd:ied, jsd:jed, k, liq_wat), q(&
&                 isd:ied, jsd:jed, k, ice_wat), q(isd:ied, jsd:jed, k, &
&                 rainwat), q(isd:ied, jsd:jed, k, snowwat), q(isd:ied, &
&                 jsd:jed, k, graupel), dpln, delz(isd:ied, jsd:jed, k)&
&                 , pt(isd:ied, jsd:jed, k), delp(isd:ied, jsd:jed, k), &
&                 q_con(isd:ied, jsd:jed, k), cappa(isd:ied, jsd:jed, k)&
&                 , gridstruct%area_64, dtdt(is:ie, js:je, k), out_dt, &
&                 last_step, cld_amt .GT. 0, q(isd:ied, jsd:jed, k, &
&                 cld_amt))
        IF (.NOT.hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              arg1_tl = (rrg*delp_tl(i, j, k)*delz(i, j, k)-rrg*delp(i, &
&               j, k)*delz_tl(i, j, k))*pt(i, j, k)/delz(i, j, k)**2 + &
&               rrg*delp(i, j, k)*pt_tl(i, j, k)/delz(i, j, k)
              arg1 = rrg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
              arg2_tl = akap*arg1_tl/arg1
              arg2 = akap*LOG(arg1)
              pkz_tl(i, j, k) = arg2_tl*EXP(arg2)
              pkz(i, j, k) = EXP(arg2)
            END DO
          END DO
        END IF
      END DO
! OpenMP k-loop
      IF (fast_mp_consv) THEN
!$OMP do
        DO j=js,je
          DO i=is,ie
            DO k=kmp,km
              te0_2d_tl(i, j) = te0_2d_tl(i, j) + te_tl(i, j, k)
              te0_2d(i, j) = te0_2d(i, j) + te(i, j, k)
            END DO
          END DO
        END DO
      END IF
      CALL TIMING_OFF('sat_adj2')
    END IF
! last_step
    IF (last_step) THEN
! Output temperature if last_step
      IF (remap_t) THEN
!$OMP do
        DO k=1,km
          DO j=js,je
            IF (.NOT.adiabatic) THEN
              DO i=is,ie
                pt_tl(i, j, k) = ((pt_tl(i, j, k)+dtmp_tl*pkz(i, j, k)+&
&                 dtmp*pkz_tl(i, j, k))*(1.+r_vir*q(i, j, k, sphum))-(pt&
&                 (i, j, k)+dtmp*pkz(i, j, k))*r_vir*q_tl(i, j, k, sphum&
&                 ))/(1.+r_vir*q(i, j, k, sphum))**2
                pt(i, j, k) = (pt(i, j, k)+dtmp*pkz(i, j, k))/(1.+r_vir*&
&                 q(i, j, k, sphum))
              END DO
            END IF
          END DO
        END DO
      ELSE IF (remap_pt) THEN
! j-loop
! k-loop
!$OMP do
        DO k=1,km
          DO j=js,je
            DO i=is,ie
              pt_tl(i, j, k) = (((pt_tl(i, j, k)+dtmp_tl)*pkz(i, j, k)+(&
&               pt(i, j, k)+dtmp)*pkz_tl(i, j, k))*(1.+r_vir*q(i, j, k, &
&               sphum))-(pt(i, j, k)+dtmp)*pkz(i, j, k)*r_vir*q_tl(i, j&
&               , k, sphum))/(1.+r_vir*q(i, j, k, sphum))**2
              pt(i, j, k) = (pt(i, j, k)+dtmp)*pkz(i, j, k)/(1.+r_vir*q(&
&               i, j, k, sphum))
            END DO
          END DO
        END DO
      ELSE IF (remap_te) THEN
!$OMP do
        DO j=js,je
          DO i=is,ie
            gz_tl(i) = 0.0
            gz(i) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              tpe_tl = te_tl(i, j, k) - gz_tl(i) - 0.25*gridstruct%rsin2&
&               (i, j)*(2*u(i, j, k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i&
&               , j+1, k)+2*v(i, j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl&
&               (i+1, j, k)-gridstruct%cosa_s(i, j)*((u_tl(i, j, k)+u_tl&
&               (i, j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i, &
&               j+1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k))))
              tpe = te(i, j, k) - gz(i) - 0.25*gridstruct%rsin2(i, j)*(u&
&               (i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)&
&               **2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&               gridstruct%cosa_s(i, j))
              dlnp_tl = rg*(peln_tl(i, k+1, j)-peln_tl(i, k, j))
              dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
              tmp_tl = (tpe_tl*(cp-pe(i, k, j)*dlnp/delp(i, j, k))*(1.+&
&               r_vir*q(i, j, k, sphum))-tpe*((cp-pe(i, k, j)*dlnp/delp(&
&               i, j, k))*r_vir*q_tl(i, j, k, sphum)-((pe_tl(i, k, j)*&
&               dlnp+pe(i, k, j)*dlnp_tl)*delp(i, j, k)-pe(i, k, j)*dlnp&
&               *delp_tl(i, j, k))*(1.+r_vir*q(i, j, k, sphum))/delp(i, &
&               j, k)**2))/((cp-pe(i, k, j)*dlnp/delp(i, j, k))*(1.+&
&               r_vir*q(i, j, k, sphum)))**2
              tmp = tpe/((cp-pe(i, k, j)*dlnp/delp(i, j, k))*(1.+r_vir*q&
&               (i, j, k, sphum)))
              pt_tl(i, j, k) = tmp_tl + ((dtmp_tl*pkz(i, j, k)+dtmp*&
&               pkz_tl(i, j, k))*(1.+r_vir*q(i, j, k, sphum))-dtmp*pkz(i&
&               , j, k)*r_vir*q_tl(i, j, k, sphum))/(1.+r_vir*q(i, j, k&
&               , sphum))**2
              pt(i, j, k) = tmp + dtmp*pkz(i, j, k)/(1.+r_vir*q(i, j, k&
&               , sphum))
              gz_tl(i) = gz_tl(i) + (dlnp_tl*tmp+dlnp*tmp_tl)*(1.+r_vir*&
&               q(i, j, k, sphum)) + dlnp*tmp*r_vir*q_tl(i, j, k, sphum)
              gz(i) = gz(i) + dlnp*tmp*(1.+r_vir*q(i, j, k, sphum))
            END DO
          END DO
        END DO
      END IF
! end k-loop
      IF (flagstruct%fv_debug) CALL PRT_MXM('remap-3  TA', pt, is, ie, &
&                                     js, je, ng, km, 1., gridstruct%&
&                                     area_64, domain)
    ELSE
! not last_step
      IF (remap_t) THEN
!$OMP do
        DO k=1,km
          DO j=js,je
            DO i=is,ie
              pt_tl(i, j, k) = (pt_tl(i, j, k)*pkz(i, j, k)-pt(i, j, k)*&
&               pkz_tl(i, j, k))/pkz(i, j, k)**2
              pt(i, j, k) = pt(i, j, k)/pkz(i, j, k)
            END DO
          END DO
        END DO
      ELSE IF (remap_te) THEN
!$OMP do
        DO j=js,je
          DO i=is,ie
            gz_tl(i) = 0.0
            gz(i) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              tpe_tl = te_tl(i, j, k) - gz_tl(i) - 0.25*gridstruct%rsin2&
&               (i, j)*(2*u(i, j, k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i&
&               , j+1, k)+2*v(i, j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl&
&               (i+1, j, k)-gridstruct%cosa_s(i, j)*((u_tl(i, j, k)+u_tl&
&               (i, j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i, &
&               j+1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k))))
              tpe = te(i, j, k) - gz(i) - 0.25*gridstruct%rsin2(i, j)*(u&
&               (i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)&
&               **2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&               gridstruct%cosa_s(i, j))
              dlnp_tl = rg*(peln_tl(i, k+1, j)-peln_tl(i, k, j))
              dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
              tmp_tl = (tpe_tl*(cp-pe(i, k, j)*dlnp/delp(i, j, k))+tpe*(&
&               (pe_tl(i, k, j)*dlnp+pe(i, k, j)*dlnp_tl)*delp(i, j, k)-&
&               pe(i, k, j)*dlnp*delp_tl(i, j, k))/delp(i, j, k)**2)/(cp&
&               -pe(i, k, j)*dlnp/delp(i, j, k))**2
              tmp = tpe/(cp-pe(i, k, j)*dlnp/delp(i, j, k))
              pt_tl(i, j, k) = (tmp_tl*pkz(i, j, k)-tmp*pkz_tl(i, j, k))&
&               /pkz(i, j, k)**2 + dtmp_tl
              pt(i, j, k) = tmp/pkz(i, j, k) + dtmp
              gz_tl(i) = gz_tl(i) + dlnp_tl*tmp + dlnp*tmp_tl
              gz(i) = gz(i) + dlnp*tmp
            END DO
          END DO
        END DO
      END IF
! end k-loop
      IF (flagstruct%fv_debug) CALL PRT_MXM('remap-3  PT', pt, is, ie, &
&                                     js, je, ng, km, 1., gridstruct%&
&                                     area_64, domain)
    END IF
    GOTO 130
 120 PRINT*, 'TE remapping non-hydrostatic is invalid and cannot be run'
    STOP
 130 CONTINUE
  END SUBROUTINE LAGRANGIAN_TO_EULERIAN_TLM
  SUBROUTINE LAGRANGIAN_TO_EULERIAN(last_step, consv, ps, pe, delp, pkz&
&   , pk, mdt, pdt, km, is, ie, js, je, isd, ied, jsd, jed, nq, nwat, &
&   sphum, q_con, u, v, w, delz, pt, q, hs, r_vir, cp, akap, cappa, &
&   kord_mt, kord_wz, kord_tr, kord_tm, peln, te0_2d, ng, ua, va, omga, &
&   te, ws, fill, reproduce_sum, out_dt, dtdt, ptop, ak, bk, pfull, &
&   flagstruct, gridstruct, domain, do_sat_adj, hydrostatic, hybrid_z, &
&   do_omega, adiabatic, do_adiabatic_init, mfx, mfy, remap_option, &
&   kord_mt_pert, kord_wz_pert, kord_tr_pert, kord_tm_pert)
    IMPLICIT NONE
!$OMP end parallel
    LOGICAL, INTENT(IN) :: last_step
! remap time step
    REAL, INTENT(IN) :: mdt
! phys time step
    REAL, INTENT(IN) :: pdt
    INTEGER, INTENT(IN) :: km
! number of tracers (including h2o)
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: nwat
! index for water vapor (specific humidity)
    INTEGER, INTENT(IN) :: sphum
    INTEGER, INTENT(IN) :: ng
! starting & ending X-Dir index
    INTEGER, INTENT(IN) :: is, ie, isd, ied
! starting & ending Y-Dir index
    INTEGER, INTENT(IN) :: js, je, jsd, jed
! Mapping order for the vector winds
    INTEGER, INTENT(IN) :: kord_mt
! Mapping order/option for w
    INTEGER, INTENT(IN) :: kord_wz
! Mapping order for tracers
    INTEGER, INTENT(IN) :: kord_tr(nq)
! Mapping order for thermodynamics
    INTEGER, INTENT(IN) :: kord_tm
! Mapping order for the vector winds
    INTEGER, INTENT(IN) :: kord_mt_pert
! Mapping order/option for w
    INTEGER, INTENT(IN) :: kord_wz_pert
! Mapping order for tracers
    INTEGER, INTENT(IN) :: kord_tr_pert(nq)
! Mapping order for thermodynamics
    INTEGER, INTENT(IN) :: kord_tm_pert
! factor for TE conservation
    REAL, INTENT(IN) :: consv
    REAL, INTENT(IN) :: r_vir
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: akap
! surface geopotential
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: te0_2d(is:ie, js:je)
    REAL, INTENT(IN) :: ws(is:ie, js:je)
    LOGICAL, INTENT(IN) :: do_sat_adj
! fill negative tracers
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
    LOGICAL, INTENT(IN) :: do_omega, adiabatic, do_adiabatic_init
    REAL, INTENT(IN) :: ptop
    REAL, INTENT(IN) :: ak(km+1)
    REAL, INTENT(IN) :: bk(km+1)
    REAL, INTENT(IN) :: pfull(km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(INOUT) :: flagstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! !INPUT/OUTPUT
! pe to the kappa
    REAL, INTENT(INOUT) :: pk(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, nq)
! pressure thickness
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, km)
! pressure at layer edges
    REAL, INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
! surface pressure
    REAL, INTENT(INOUT) :: ps(isd:ied, jsd:jed)
! u-wind will be ghosted one latitude to the north upon exit
! u-wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
! v-wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
! vertical velocity (m/s)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, km)
! cp*virtual potential temperature
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
! as input; output: temperature
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: delz
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: q_con, cappa
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: hybrid_z
    LOGICAL, INTENT(IN) :: out_dt
! u-wind (m/s) on physics grid
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
! v-wind (m/s) on physics grid
    REAL, INTENT(INOUT) :: va(isd:ied, jsd:jed, km)
! vertical press. velocity (pascal/sec)
    REAL, INTENT(INOUT) :: omga(isd:ied, jsd:jed, km)
! log(pe)
    REAL, INTENT(INOUT) :: peln(is:ie, km+1, js:je)
    REAL, INTENT(INOUT) :: dtdt(is:ie, js:je, km)
! layer-mean pk for converting t to pt
    REAL, INTENT(OUT) :: pkz(is:ie, js:je, km)
    REAL, INTENT(OUT) :: te(isd:ied, jsd:jed, km)
! Mass fluxes
! X-dir Mass Flux
    REAL, OPTIONAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, km)
! Y-dir Mass Flux
    REAL, OPTIONAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, km)
! 0: remap  T in logP
    INTEGER, INTENT(IN) :: remap_option
! 1: remap PT in P
! 3: remap TE in logP with GMAO cubic
! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
    REAL, DIMENSION(is:ie, js:je) :: te_2d, zsum0, zsum1, dpln
    REAL, DIMENSION(is:ie, km) :: q2, dp2
    REAL, DIMENSION(is:ie, km+1) :: pe1, pe2, pk1, pk2, pn2, phis
    REAL, DIMENSION(is:ie+1, km+1) :: pe0, pe3
    REAL, DIMENSION(is:ie) :: gz, cvm, qv
    REAL :: rcp, rg, tmp, tpe, rrg, bkh, dtmp, k1k, dlnp
    LOGICAL :: fast_mp_consv
    INTEGER :: i, j, k
    INTEGER :: nt, liq_wat, ice_wat, rainwat, snowwat, cld_amt, graupel&
&   , iq, n, kmp, kp, k_next
    LOGICAL :: remap_t, remap_pt, remap_te
    INTEGER :: abs_kord_tm, abs_kord_tm_pert
    INTEGER :: iep1, jep1, iedp1, jedp1
    INTRINSIC ABS
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC PRESENT
    REAL :: abs0
    REAL :: arg1
    REAL :: arg2
    REAL :: result1
    IF (kord_tm .GE. 0.) THEN
      abs_kord_tm = kord_tm
    ELSE
      abs_kord_tm = -kord_tm
    END IF
    IF (kord_tm_pert .GE. 0.) THEN
      abs_kord_tm_pert = kord_tm_pert
    ELSE
      abs_kord_tm_pert = -kord_tm_pert
    END IF
    iep1 = ie + 1
    jep1 = je + 1
    iedp1 = ied + 1
    jedp1 = jed + 1
    remap_t = .false.
    remap_pt = .false.
    remap_te = .false.
    SELECT CASE  (remap_option) 
    CASE (0) 
      remap_t = .true.
    CASE (1) 
      remap_pt = .true.
    CASE (2) 
      remap_te = .true.
    CASE DEFAULT
      PRINT*, ' INVALID REMAPPING OPTION '
      STOP
    END SELECT
    IF (IS_MASTER() .AND. flagstruct%fv_debug) THEN
      PRINT*, ''
      SELECT CASE  (remap_option) 
      CASE (0) 
        PRINT*, ' REMAPPING  T in logP '
      CASE (1) 
        PRINT*, ' REMAPPING PT in P'
      CASE (2) 
        PRINT*, ' REMAPPING TE in logP with GMAO cubic'
      END SELECT
      PRINT*, ' REMAPPING CONSV:     ', consv
      PRINT*, ' REMAPPING CONSV_MIN: ', consv_min
      PRINT*, ''
    END IF
    IF (flagstruct%fv_debug) CALL PRT_MXM('remap-0  PT', pt, is, ie, js&
&                                   , je, ng, km, 1., gridstruct%area_64&
&                                   , domain)
! akap / (1.-akap) = rg/Cv=0.4
    k1k = rdgas/cv_air
    rg = rdgas
    rcp = 1./cp
    rrg = -(rdgas/grav)
    IF (fpp%fpp_mapl_mode) THEN
      liq_wat = 2
      ice_wat = 3
      rainwat = -1
      snowwat = -1
      graupel = -1
      cld_amt = -1
    ELSE
      liq_wat = GET_TRACER_INDEX(model_atmos, 'liq_wat')
      ice_wat = GET_TRACER_INDEX(model_atmos, 'ice_wat')
      rainwat = GET_TRACER_INDEX(model_atmos, 'rainwat')
      snowwat = GET_TRACER_INDEX(model_atmos, 'snowwat')
      graupel = GET_TRACER_INDEX(model_atmos, 'graupel')
      cld_amt = GET_TRACER_INDEX(model_atmos, 'cld_amt')
    END IF
    IF (do_sat_adj) THEN
      fast_mp_consv = .NOT.do_adiabatic_init .AND. consv .GT. consv_min
      DO k=1,km
        kmp = k
        IF (pfull(k) .GT. 10.e2) GOTO 100
      END DO
 100  CALL QS_INIT(kmp)
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,km,pe,ptop,kord_tm,hydrostatic, &
!$OMP                                  pt,pk,rg,peln,q,nwat,liq_wat,rainwat,ice_wat,snowwat,    &
!$OMP                                  graupel,q_con,sphum,cappa,r_vir,rcp,k1k,delp, &
!$OMP                                  delz,akap,pkz,te,u,v,ps, gridstruct, last_step, &
!$OMP                                  ak,bk,nq,isd,ied,jsd,jed,kord_tr,fill, adiabatic, &
!$OMP                                  hs,w,ws,kord_wz,do_omega,omga,rrg,kord_mt,ua)    &
!$OMP                          private(qv,gz,cvm,kp,k_next,bkh,dp2,   &
!$OMP                                  pe0,pe1,pe2,pe3,pk1,pk2,pn2,phis,q2)
    DO j=js,je+1
      DO k=1,km+1
        DO i=is,ie
          pe1(i, k) = pe(i, k, j)
        END DO
      END DO
      DO i=is,ie
        pe2(i, 1) = ptop
        pe2(i, km+1) = pe(i, km+1, j)
      END DO
!(j < je+1)
      IF (j .NE. je + 1) THEN
        IF (remap_t) THEN
! hydro test
! Remap T in logP
! Note: pt at this stage is Theta_v
          IF (hydrostatic) THEN
! Transform virtual pt to virtual Temp
            DO k=1,km
              DO i=is,ie
                pt(i, j, k) = pt(i, j, k)*(pk(i, j, k+1)-pk(i, j, k))/(&
&                 akap*(peln(i, k+1, j)-peln(i, k, j)))
              END DO
            END DO
          ELSE
! Transform "density pt" to "density temp"
            DO k=1,km
              DO i=is,ie
                arg1 = rrg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
                arg2 = k1k*LOG(arg1)
                pt(i, j, k) = pt(i, j, k)*EXP(arg2)
              END DO
            END DO
          END IF
        ELSE IF (.NOT.remap_pt) THEN
! Using dry pressure for the definition of the virtual potential temperature
!                    pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*    &
!                                              pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
! Remap PT in P
! pt is already virtual PT
          IF (remap_te) THEN
! Remap TE in logP
! Transform virtual pt to total energy
            CALL PKEZ(km, is, ie, js, je, j, pe, pk, akap, peln, pkz, &
&               ptop)
! Compute cp*T + KE
            DO k=1,km
              DO i=is,ie
                te(i, j, k) = 0.25*gridstruct%rsin2(i, j)*(u(i, j, k)**2&
&                 +u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j&
&                 , k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&                 gridstruct%cosa_s(i, j)) + cp_air*pt(i, j, k)*pkz(i, j&
&                 , k)
              END DO
            END DO
          END IF
        END IF
        IF (.NOT.hydrostatic) THEN
          DO k=1,km
            DO i=is,ie
! ="specific volume"/grav
              delz(i, j, k) = -(delz(i, j, k)/delp(i, j, k))
            END DO
          END DO
        END IF
! update ps
        DO i=is,ie
          ps(i, j) = pe1(i, km+1)
        END DO
!
! Hybrid sigma-P coordinate:
!
        DO k=2,km
          DO i=is,ie
            pe2(i, k) = ak(k) + bk(k)*pe(i, km+1, j)
          END DO
        END DO
        DO k=1,km
          DO i=is,ie
            dp2(i, k) = pe2(i, k+1) - pe2(i, k)
          END DO
        END DO
!------------
! update delp
!------------
        DO k=1,km
          DO i=is,ie
            delp(i, j, k) = dp2(i, k)
          END DO
        END DO
!------------------
! Compute p**Kappa
!------------------
        DO k=1,km+1
          DO i=is,ie
            pk1(i, k) = pk(i, j, k)
          END DO
        END DO
        DO i=is,ie
          pn2(i, 1) = peln(i, 1, j)
          pn2(i, km+1) = peln(i, km+1, j)
          pk2(i, 1) = pk1(i, 1)
          pk2(i, km+1) = pk1(i, km+1)
        END DO
        DO k=2,km
          DO i=is,ie
            pn2(i, k) = LOG(pe2(i, k))
            pk2(i, k) = EXP(akap*pn2(i, k))
          END DO
        END DO
        IF (remap_t) THEN
!----------------------------------
! Map t using logp
!----------------------------------
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            CALL MAP_SCALAR(km, peln(is:ie, 1:km+1, j), gz, km, pn2, &
&                        pt, is, ie, j, isd, ied, jsd, jed, 1, &
&                        abs_kord_tm, t_min)
          ELSE
      call map_scalar(km,  peln(is:ie,1:km+1,j),      gz,   &
                      km,  pn2,           pt,              &
                      is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm, t_min)
          END IF
        ELSE IF (remap_pt) THEN
!----------------------------------
! Map pt using pe
!----------------------------------
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            CALL MAP1_PPM(km, pe1, gz, km, pe2, pt, is, ie, j, isd, &
&                      ied, jsd, jed, 1, abs_kord_tm)
          ELSE
      call map1_ppm (km,  pe1,       gz,       &
                     km,  pe2,  pt,                  &
                     is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm)
          END IF
        ELSE IF (remap_te) THEN
!----------------------------------
! map Total Energy using GMAO cubic
!----------------------------------
          DO i=is,ie
            phis(i, km+1) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              phis(i, k) = phis(i, k+1) + cp_air*pt(i, j, k)*(pk1(i, k+1&
&               )-pk1(i, k))
            END DO
          END DO
          DO k=1,km+1
            DO i=is,ie
              phis(i, k) = phis(i, k)*pe1(i, k)
            END DO
          END DO
          DO k=1,km
            DO i=is,ie
              te(i, j, k) = te(i, j, k) + (phis(i, k+1)-phis(i, k))/(pe1&
&               (i, k+1)-pe1(i, k))
            END DO
          END DO
! Map te using log P in GMAO cubic
          CALL MAP1_CUBIC(km, pe1, km, pe2, te, is, ie, j, isd, ied, jsd&
&                   , jed, akap, 1, .true.)
        END IF
!----------------
! Map constituents
!----------------
        IF (nq .GT. 5) THEN
          IF (kord_tr(1) .EQ. kord_tr_pert(1)) THEN
            CALL MAPN_TRACER(nq, km, pe1, pe2, q, dp2, kord_tr, j, is&
&                         , ie, isd, ied, jsd, jed, 0., fill)
          ELSE
      call mapn_tracer(nq, km, pe1, pe2, q, dp2, kord_tr, j,     &
                       is, ie, isd, ied, jsd, jed, 0., fill)
          END IF
        ELSE IF (nq .GT. 0) THEN
! Remap one tracer at a time
          DO iq=1,nq
            IF (kord_tr(iq) .EQ. kord_tr_pert(iq)) THEN
              CALL MAP1_Q2(km, pe1, q(isd:ied, jsd:jed, 1:km, iq), km&
&                       , pe2, q2, dp2, is, ie, 0, kord_tr(iq), j, isd, &
&                       ied, jsd, jed, 0.)
            ELSE
      call map1_q2(km, pe1, q(isd:ied,jsd:jed,1:km,iq),     &
                   km, pe2, q2, dp2,             &
                   is, ie, 0, kord_tr(iq), j, isd, ied, jsd, jed, 0.)
            END IF
            IF (fill) CALL FILLZ(ie - is + 1, km, 1, q2, dp2)
            DO k=1,km
              DO i=is,ie
                q(i, j, k, iq) = q2(i, k)
              END DO
            END DO
          END DO
        END IF
        IF (.NOT.hydrostatic) THEN
! Remap vertical wind:
          IF (kord_wz .EQ. kord_wz_pert) THEN
            CALL MAP1_PPM(km, pe1, ws(is:ie, j), km, pe2, w, is, ie, &
&                      j, isd, ied, jsd, jed, -2, kord_wz)
          ELSE
      call map1_ppm (km,   pe1,      ws(is:ie,j),   &
                     km,   pe2,  w,              &
                     is, ie, j, isd, ied, jsd, jed, -2, kord_wz)
          END IF
! Remap delz for hybrid sigma-p coordinate
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            CALL MAP1_PPM(km, pe1, gz, km, pe2, delz, is, ie, j, isd&
&                      , ied, jsd, jed, 1, abs_kord_tm)
          ELSE
      call map1_ppm (km,   pe1,        gz,   &
                     km,   pe2, delz,              &
                     is, ie, j, isd,  ied,  jsd,  jed,  1, abs_kord_tm)
          END IF
          DO k=1,km
            DO i=is,ie
              delz(i, j, k) = -(delz(i, j, k)*dp2(i, k))
            END DO
          END DO
        END IF
!----------
! Update pk
!----------
        DO k=1,km+1
          DO i=is,ie
            pk(i, j, k) = pk2(i, k)
          END DO
        END DO
!----------------
        IF (do_omega) THEN
! Start do_omega
! Copy omega field to pe3
          DO i=is,ie
            pe3(i, 1) = 0.
          END DO
          DO k=2,km+1
            DO i=is,ie
              pe3(i, k) = omga(i, j, k-1)
            END DO
          END DO
        END IF
        DO k=1,km+1
          DO i=is,ie
            pe0(i, k) = peln(i, k, j)
            peln(i, k, j) = pn2(i, k)
          END DO
        END DO
!------------
! Compute pkz
!------------
        IF (hydrostatic) THEN
          DO k=1,km
            DO i=is,ie
              pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1&
&               , j)-peln(i, k, j)))
            END DO
          END DO
        ELSE IF (remap_te) THEN
! WMP: note that this is where TE remapping non-hydrostatic is invalid and cannot be run
          PRINT*, &
&         'TE remapping non-hydrostatic is invalid and cannot be run'
          STOP
        ELSE IF (remap_t) THEN
! Note: pt at this stage is T_v or T_m
          DO k=1,km
            DO i=is,ie
              arg1 = rrg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
              arg2 = akap*LOG(arg1)
              pkz(i, j, k) = EXP(arg2)
            END DO
          END DO
        ELSE
! Using dry pressure for the definition of the virtual potential temperature
!           pkz(i,j,k) = exp(akap*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
! Note: pt at this stage is Theta_v
          DO k=1,km
            DO i=is,ie
              arg1 = rrg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
              arg2 = k1k*LOG(arg1)
              pkz(i, j, k) = EXP(arg2)
            END DO
          END DO
        END IF
! end do_omega
! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
        IF (do_omega) THEN
          DO k=1,km
            DO i=is,ie
              dp2(i, k) = 0.5*(peln(i, k, j)+peln(i, k+1, j))
            END DO
          END DO
          DO i=is,ie
            k_next = 1
            DO 110 n=1,km
              kp = k_next
              DO k=kp,km
                IF (dp2(i, n) .LE. pe0(i, k+1) .AND. dp2(i, n) .GE. pe0(&
&                   i, k)) THEN
                  omga(i, j, n) = pe3(i, k) + (pe3(i, k+1)-pe3(i, k))*(&
&                   dp2(i, n)-pe0(i, k))/(pe0(i, k+1)-pe0(i, k))
                  k_next = k
                  GOTO 110
                END IF
              END DO
 110        CONTINUE
          END DO
        END IF
      END IF
      DO i=is,ie+1
        pe0(i, 1) = pe(i, 1, j)
      END DO
!------
! map u
!------
      DO k=2,km+1
        DO i=is,ie
          pe0(i, k) = 0.5*(pe(i, k, j-1)+pe1(i, k))
        END DO
      END DO
      DO k=1,km+1
        bkh = 0.5*bk(k)
        DO i=is,ie
          pe3(i, k) = ak(k) + bkh*(pe(i, km+1, j-1)+pe1(i, km+1))
        END DO
      END DO
      IF (kord_mt .EQ. kord_mt_pert) THEN
        CALL MAP1_PPM(km, pe0(is:ie, :), gz, km, pe3(is:ie, :), u, is&
&                  , ie, j, isd, ied, jsd, jedp1, -1, kord_mt)
      ELSE
      call map1_ppm( km, pe0(is:ie,:),        gz,   &
                     km, pe3(is:ie,:),   u,               &
                     is, ie, j, isd, ied, jsd, jedp1, -1, kord_mt)
      END IF
      IF (PRESENT(mfy)) CALL MAP1_PPM(km, pe0(is:ie, :), gz, km, pe3(is:&
&                               ie, :), mfy, is, ie, j, is, ie, js, jep1&
&                               , -1, kord_mt)
! (j < je+1)
      IF (j .LT. je + 1) THEN
!------
! map v
!------
        DO i=is,ie+1
          pe3(i, 1) = ak(1)
        END DO
        DO k=2,km+1
          bkh = 0.5*bk(k)
          DO i=is,ie+1
            pe0(i, k) = 0.5*(pe(i-1, k, j)+pe(i, k, j))
            pe3(i, k) = ak(k) + bkh*(pe(i-1, km+1, j)+pe(i, km+1, j))
          END DO
        END DO
        IF (kord_mt .EQ. kord_mt_pert) THEN
          CALL MAP1_PPM(km, pe0, gz, km, pe3, v, is, iep1, j, isd, &
&                    iedp1, jsd, jed, -1, kord_mt)
        ELSE
      call map1_ppm (km, pe0,     gz,    &
                     km, pe3,  v, is, ie+1,    &
                     j, isd, iedp1, jsd, jed, -1, kord_mt)
        END IF
        IF (PRESENT(mfx)) CALL MAP1_PPM(km, pe0, gz, km, pe3, mfx, is, &
&                                 iep1, j, is, iep1, js, je, -1, kord_mt&
&                                )
      END IF
      DO k=1,km
        DO i=is,ie
          ua(i, j, k) = pe2(i, k+1)
        END DO
      END DO
    END DO
!$OMP parallel default(none) shared(is,ie,js,je,km,kmp,ptop,u,v,pe,ua,isd,ied,jsd,jed,kord_mt, &
!$OMP                               te_2d,te,delp,hydrostatic,hs,rg,pt,peln, adiabatic, &
!$OMP                               cp,delz,nwat,rainwat,liq_wat,ice_wat,snowwat,       &
!$OMP                               graupel,q_con,r_vir,sphum,w,pk,pkz,last_step,consv, &
!$OMP                               do_adiabatic_init,zsum1,zsum0,te0_2d,domain,        &
!$OMP                               ng,gridstruct,E_Flux,pdt,dtmp,reproduce_sum,q,      &
!$OMP                               mdt,cld_amt,cappa,dtdt,out_dt,rrg,akap,do_sat_adj,  &
!$OMP                               fast_mp_consv,kord_tm) &
!$OMP                       private(pe0,pe1,pe2,pe3,qv,cvm,gz,phis,tpe,tmp, dpln)
!$OMP do
    DO k=2,km
      DO j=js,je
        DO i=is,ie
          pe(i, k, j) = ua(i, j, k-1)
        END DO
      END DO
    END DO
    IF (flagstruct%fv_debug) THEN
      IF (kord_tm .LT. 0) THEN
        CALL PRT_MXM('remap-1  TV', pt, is, ie, js, je, ng, km, 1., &
&              gridstruct%area_64, domain)
      ELSE
        CALL PRT_MXM('remap-1  PT', pt, is, ie, js, je, ng, km, 1., &
&              gridstruct%area_64, domain)
      END IF
    END IF
    dtmp = 0.
! end last_step check
    IF (last_step .AND. (.NOT.do_adiabatic_init)) THEN
! end consv check
      IF (consv .GT. consv_min) THEN
!$OMP do
        DO j=js,je
          IF (remap_t) THEN
! end non-hydro
            IF (hydrostatic) THEN
              DO i=is,ie
                gz(i) = hs(i, j)
                DO k=1,km
                  gz(i) = gz(i) + rg*pt(i, j, k)*(peln(i, k+1, j)-peln(i&
&                   , k, j))
                END DO
              END DO
              DO i=is,ie
                te_2d(i, j) = pe(i, km+1, j)*hs(i, j) - pe(i, 1, j)*gz(i&
&                 )
              END DO
              DO k=1,km
                DO i=is,ie
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cp*pt(i, j&
&                   , k)+0.25*gridstruct%rsin2(i, j)*(u(i, j, k)**2+u(i&
&                   , j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, &
&                   k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&                   gridstruct%cosa_s(i, j)))
                END DO
              END DO
            ELSE
              DO i=is,ie
                te_2d(i, j) = 0.
                phis(i, km+1) = hs(i, j)
              END DO
              DO k=km,1,-1
                DO i=is,ie
                  phis(i, k) = phis(i, k+1) - grav*delz(i, j, k)
                END DO
              END DO
              DO k=1,km
                DO i=is,ie
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv_air*pt(i&
&                   , j, k)/(1.+r_vir*q(i, j, k, sphum))+0.5*(phis(i, k)&
&                   +phis(i, k+1)+w(i, j, k)**2+0.5*gridstruct%rsin2(i, &
&                   j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+&
&                   1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(&
&                   i+1, j, k))*gridstruct%cosa_s(i, j))))
                END DO
              END DO
            END IF
          ELSE IF (remap_pt) THEN
! k-loop
            IF (hydrostatic) THEN
              DO i=is,ie
                gz(i) = hs(i, j)
                DO k=1,km
                  gz(i) = gz(i) + cp_air*pt(i, j, k)*(pk(i, j, k+1)-pk(i&
&                   , j, k))
                END DO
              END DO
              DO i=is,ie
                te_2d(i, j) = pe(i, km+1, j)*hs(i, j) - pe(i, 1, j)*gz(i&
&                 )
              END DO
              DO k=1,km
                DO i=is,ie
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cp_air*pt(i&
&                   , j, k)*pkz(i, j, k)+0.25*gridstruct%rsin2(i, j)*(u(&
&                   i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, &
&                   k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j&
&                   , k))*gridstruct%cosa_s(i, j)))
                END DO
              END DO
            ELSE
!-----------------
! Non-hydrostatic:
!-----------------
              DO i=is,ie
                phis(i, km+1) = hs(i, j)
                DO k=km,1,-1
                  phis(i, k) = phis(i, k+1) - grav*delz(i, j, k)
                END DO
              END DO
              DO i=is,ie
                te_2d(i, j) = 0.
              END DO
              DO k=1,km
                DO i=is,ie
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv_air*pt(i&
&                   , j, k)/(1.+r_vir*q(i, j, k, sphum))+0.5*(phis(i, k)&
&                   +phis(i, k+1)+w(i, j, k)**2+0.5*gridstruct%rsin2(i, &
&                   j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+&
&                   1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(&
&                   i+1, j, k))*gridstruct%cosa_s(i, j))))
                END DO
              END DO
            END IF
          ELSE IF (remap_te) THEN
            DO i=is,ie
              te_2d(i, j) = te(i, j, 1)*delp(i, j, 1)
            END DO
            DO k=2,km
              DO i=is,ie
                te_2d(i, j) = te_2d(i, j) + te(i, j, k)*delp(i, j, k)
              END DO
            END DO
          END IF
          DO i=is,ie
            te_2d(i, j) = te0_2d(i, j) - te_2d(i, j)
            zsum1(i, j) = pkz(i, j, 1)*delp(i, j, 1)
          END DO
          DO k=2,km
            DO i=is,ie
              zsum1(i, j) = zsum1(i, j) + pkz(i, j, k)*delp(i, j, k)
            END DO
          END DO
          IF (hydrostatic) THEN
            DO i=is,ie
              zsum0(i, j) = ptop*(pk(i, j, 1)-pk(i, j, km+1)) + zsum1(i&
&               , j)
            END DO
          END IF
        END DO
! j-loop
!$OMP single
        result1 = G_SUM(domain, te_2d, is, ie, js, je, ng, gridstruct%&
&         area_64, 0, .true.)
        tpe = consv*result1
! unit: W/m**2
        e_flux = tpe/(grav*pdt*4.*pi*radius**2)
! Note pdt is "phys" time step
        IF (hydrostatic) THEN
          result1 = G_SUM(domain, zsum0, is, ie, js, je, ng, gridstruct%&
&           area_64, 0, .true.)
          dtmp = tpe/(cp*result1)
        ELSE
          result1 = G_SUM(domain, zsum1, is, ie, js, je, ng, gridstruct%&
&           area_64, 0, .true.)
          dtmp = tpe/(cv_air*result1)
        END IF
      ELSE IF (consv .LT. -consv_min) THEN
!$OMP end single
!$OMP do
        DO j=js,je
          DO i=is,ie
            zsum1(i, j) = pkz(i, j, 1)*delp(i, j, 1)
          END DO
          DO k=2,km
            DO i=is,ie
              zsum1(i, j) = zsum1(i, j) + pkz(i, j, k)*delp(i, j, k)
            END DO
          END DO
          IF (hydrostatic) THEN
            DO i=is,ie
              zsum0(i, j) = ptop*(pk(i, j, 1)-pk(i, j, km+1)) + zsum1(i&
&               , j)
            END DO
          END IF
        END DO
        e_flux = consv
!$OMP single
        IF (hydrostatic) THEN
          result1 = G_SUM(domain, zsum0, is, ie, js, je, ng, gridstruct%&
&           area_64, 0, .true.)
          dtmp = e_flux*(grav*pdt*4.*pi*radius**2)/(cp*result1)
        ELSE
          result1 = G_SUM(domain, zsum1, is, ie, js, je, ng, gridstruct%&
&           area_64, 0, .true.)
          dtmp = e_flux*(grav*pdt*4.*pi*radius**2)/(cv_air*result1)
        END IF
      END IF
    END IF
!$OMP end single
! do_sat_adj
! Note: pt at this stage is T_v
    IF (remap_t .AND. (.NOT.do_adiabatic_init) .AND. do_sat_adj) THEN
! if ( do_sat_adj ) then
      CALL TIMING_ON('sat_adj2')
!$OMP do
      DO k=kmp,km
        DO j=js,je
          DO i=is,ie
            dpln(i, j) = peln(i, k+1, j) - peln(i, k, j)
          END DO
        END DO
        IF (mdt .GE. 0.) THEN
          abs0 = mdt
        ELSE
          abs0 = -mdt
        END IF
        CALL FV_SAT_ADJ(abs0, r_vir, is, ie, js, je, ng, hydrostatic, &
&                 fast_mp_consv, te(isd:ied, jsd:jed, k), q(isd:ied, jsd&
&                 :jed, k, sphum), q(isd:ied, jsd:jed, k, liq_wat), q(&
&                 isd:ied, jsd:jed, k, ice_wat), q(isd:ied, jsd:jed, k, &
&                 rainwat), q(isd:ied, jsd:jed, k, snowwat), q(isd:ied, &
&                 jsd:jed, k, graupel), dpln, delz(isd:ied, jsd:jed, k)&
&                 , pt(isd:ied, jsd:jed, k), delp(isd:ied, jsd:jed, k), &
&                 q_con(isd:ied, jsd:jed, k), cappa(isd:ied, jsd:jed, k)&
&                 , gridstruct%area_64, dtdt(is:ie, js:je, k), out_dt, &
&                 last_step, cld_amt .GT. 0, q(isd:ied, jsd:jed, k, &
&                 cld_amt))
        IF (.NOT.hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              arg1 = rrg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
              arg2 = akap*LOG(arg1)
              pkz(i, j, k) = EXP(arg2)
            END DO
          END DO
        END IF
      END DO
! OpenMP k-loop
      IF (fast_mp_consv) THEN
!$OMP do
        DO j=js,je
          DO i=is,ie
            DO k=kmp,km
              te0_2d(i, j) = te0_2d(i, j) + te(i, j, k)
            END DO
          END DO
        END DO
      END IF
      CALL TIMING_OFF('sat_adj2')
    END IF
! last_step
    IF (last_step) THEN
! Output temperature if last_step
      IF (remap_t) THEN
!$OMP do
        DO k=1,km
          DO j=js,je
            IF (.NOT.adiabatic) THEN
              DO i=is,ie
                pt(i, j, k) = (pt(i, j, k)+dtmp*pkz(i, j, k))/(1.+r_vir*&
&                 q(i, j, k, sphum))
              END DO
            END IF
          END DO
        END DO
      ELSE IF (remap_pt) THEN
! j-loop
! k-loop
!$OMP do
        DO k=1,km
          DO j=js,je
            DO i=is,ie
              pt(i, j, k) = (pt(i, j, k)+dtmp)*pkz(i, j, k)/(1.+r_vir*q(&
&               i, j, k, sphum))
            END DO
          END DO
        END DO
      ELSE IF (remap_te) THEN
!$OMP do
        DO j=js,je
          DO i=is,ie
            gz(i) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              tpe = te(i, j, k) - gz(i) - 0.25*gridstruct%rsin2(i, j)*(u&
&               (i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)&
&               **2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&               gridstruct%cosa_s(i, j))
              dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
              tmp = tpe/((cp-pe(i, k, j)*dlnp/delp(i, j, k))*(1.+r_vir*q&
&               (i, j, k, sphum)))
              pt(i, j, k) = tmp + dtmp*pkz(i, j, k)/(1.+r_vir*q(i, j, k&
&               , sphum))
              gz(i) = gz(i) + dlnp*tmp*(1.+r_vir*q(i, j, k, sphum))
            END DO
          END DO
        END DO
      END IF
! end k-loop
      IF (flagstruct%fv_debug) CALL PRT_MXM('remap-3  TA', pt, is, ie, &
&                                     js, je, ng, km, 1., gridstruct%&
&                                     area_64, domain)
    ELSE
! not last_step
      IF (remap_t) THEN
!$OMP do
        DO k=1,km
          DO j=js,je
            DO i=is,ie
              pt(i, j, k) = pt(i, j, k)/pkz(i, j, k)
            END DO
          END DO
        END DO
      ELSE IF (remap_te) THEN
!$OMP do
        DO j=js,je
          DO i=is,ie
            gz(i) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              tpe = te(i, j, k) - gz(i) - 0.25*gridstruct%rsin2(i, j)*(u&
&               (i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)&
&               **2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&               gridstruct%cosa_s(i, j))
              dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
              tmp = tpe/(cp-pe(i, k, j)*dlnp/delp(i, j, k))
              pt(i, j, k) = tmp/pkz(i, j, k) + dtmp
              gz(i) = gz(i) + dlnp*tmp
            END DO
          END DO
        END DO
      END IF
! end k-loop
      IF (flagstruct%fv_debug) CALL PRT_MXM('remap-3  PT', pt, is, ie, &
&                                     js, je, ng, km, 1., gridstruct%&
&                                     area_64, domain)
    END IF
  END SUBROUTINE LAGRANGIAN_TO_EULERIAN
!  Differentiation of compute_total_energy in forward (tangent) mode:
!   variations   of useful results: teq te_2d
!   with respect to varying inputs: qc peln q u v w delp delz pe
!                pt
  SUBROUTINE COMPUTE_TOTAL_ENERGY_TLM(is, ie, js, je, isd, ied, jsd, jed&
&   , km, u, u_tl, v, v_tl, w, w_tl, delz, delz_tl, pt, pt_tl, delp, &
&   delp_tl, q, q_tl, qc, qc_tl, pe, pe_tl, peln, peln_tl, hs, rsin2_l, &
&   cosa_s_l, r_vir, cp, rg, hlv, te_2d, te_2d_tl, ua, va, teq, teq_tl, &
&   moist_phys, nwat, sphum, liq_wat, rainwat, ice_wat, snowwat, graupel&
&   , hydrostatic, id_te)
    IMPLICIT NONE
!------------------------------------------------------
! Compute vertically integrated total energy per column
!------------------------------------------------------
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: km, is, ie, js, je, isd, ied, jsd, jed, id_te
    INTEGER, INTENT(IN) :: sphum, liq_wat, ice_wat, rainwat, snowwat, &
&   graupel, nwat
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt, delp
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt_tl, delp_tl
    REAL, DIMENSION(isd:ied, jsd:jed, km, *), INTENT(IN) :: q
    REAL, DIMENSION(isd:ied, jsd:jed, km, *), INTENT(IN) :: q_tl
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: qc
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: qc_tl
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: u_tl(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_tl(isd:ied+1, jsd:jed, km)
! vertical velocity (m/s)
    REAL, INTENT(IN) :: w(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: w_tl(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delz(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delz_tl(isd:ied, jsd:jed, km)
! surface geopotential
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
! pressure at layer edges
    REAL, INTENT(IN) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(IN) :: pe_tl(is-1:ie+1, km+1, js-1:je+1)
! log(pe)
    REAL, INTENT(IN) :: peln(is:ie, km+1, js:je)
    REAL, INTENT(IN) :: peln_tl(is:ie, km+1, js:je)
    REAL, INTENT(IN) :: cp, rg, r_vir, hlv
    REAL, INTENT(IN) :: rsin2_l(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: cosa_s_l(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: moist_phys, hydrostatic
! Output:
! vertically integrated TE
    REAL, INTENT(OUT) :: te_2d(is:ie, js:je)
    REAL, INTENT(OUT) :: te_2d_tl(is:ie, js:je)
! Moist TE
    REAL, INTENT(OUT) :: teq(is:ie, js:je)
    REAL, INTENT(OUT) :: teq_tl(is:ie, js:je)
! Local
    REAL, DIMENSION(is:ie, km) :: tv
    REAL, DIMENSION(is:ie, km) :: tv_tl
    REAL :: phiz(is:ie, km+1)
    REAL :: phiz_tl(is:ie, km+1)
    REAL :: cvm(is:ie), qd(is:ie)
    INTEGER :: i, j, k
    te_2d_tl = 0.0
    phiz_tl = 0.0
    tv_tl = 0.0
!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, flagstruct%c2l_ord)
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,km,hydrostatic,hs,pt,qc,rg,peln,te_2d, &
!$OMP                                  pe,delp,cp,rsin2_l,u,v,cosa_s_l,delz,moist_phys,w, &
!$OMP                                  q,nwat,liq_wat,rainwat,ice_wat,snowwat,graupel,sphum)   &
!$OMP                          private(phiz, tv, cvm, qd)
    DO j=js,je
      IF (hydrostatic) THEN
        DO i=is,ie
          phiz_tl(i, km+1) = 0.0
          phiz(i, km+1) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            tv_tl(i, k) = pt_tl(i, j, k)*(1.+qc(i, j, k)) + pt(i, j, k)*&
&             qc_tl(i, j, k)
            tv(i, k) = pt(i, j, k)*(1.+qc(i, j, k))
            phiz_tl(i, k) = phiz_tl(i, k+1) + rg*(tv_tl(i, k)*(peln(i, k&
&             +1, j)-peln(i, k, j))+tv(i, k)*(peln_tl(i, k+1, j)-peln_tl&
&             (i, k, j)))
            phiz(i, k) = phiz(i, k+1) + rg*tv(i, k)*(peln(i, k+1, j)-&
&             peln(i, k, j))
          END DO
        END DO
        DO i=is,ie
          te_2d_tl(i, j) = pe_tl(i, km+1, j)*phiz(i, km+1) + pe(i, km+1&
&           , j)*phiz_tl(i, km+1) - pe_tl(i, 1, j)*phiz(i, 1) - pe(i, 1&
&           , j)*phiz_tl(i, 1)
          te_2d(i, j) = pe(i, km+1, j)*phiz(i, km+1) - pe(i, 1, j)*phiz(&
&           i, 1)
        END DO
        DO k=1,km
          DO i=is,ie
            te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(cp*tv(i&
&             , k)+0.25*rsin2_l(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i&
&             , j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i&
&             , j, k)+v(i+1, j, k))*cosa_s_l(i, j))) + delp(i, j, k)*(cp&
&             *tv_tl(i, k)+0.25*rsin2_l(i, j)*(2*u(i, j, k)*u_tl(i, j, k&
&             )+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i, j, k)*v_tl(i, j, k&
&             )+2*v(i+1, j, k)*v_tl(i+1, j, k)-cosa_s_l(i, j)*((u_tl(i, &
&             j, k)+u_tl(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, &
&             k)+u(i, j+1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k)))))
            te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cp*tv(i, k)+0.25*&
&             rsin2_l(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2&
&             +v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i&
&             +1, j, k))*cosa_s_l(i, j)))
          END DO
        END DO
      ELSE
!-----------------
! Non-hydrostatic:
!-----------------
        DO i=is,ie
          phiz_tl(i, km+1) = 0.0
          phiz(i, km+1) = hs(i, j)
          DO k=km,1,-1
            phiz_tl(i, k) = phiz_tl(i, k+1) - grav*delz_tl(i, j, k)
            phiz(i, k) = phiz(i, k+1) - grav*delz(i, j, k)
          END DO
        END DO
        DO i=is,ie
          te_2d_tl(i, j) = 0.0
          te_2d(i, j) = 0.
        END DO
        IF (moist_phys) THEN
          DO k=1,km
            DO i=is,ie
              te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(cv_air&
&               *pt(i, j, k)+0.5*(phiz(i, k)+phiz(i, k+1)+w(i, j, k)**2+&
&               0.5*rsin2_l(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j&
&               , k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, &
&               j, k)+v(i+1, j, k))*cosa_s_l(i, j)))) + delp(i, j, k)*(&
&               cv_air*pt_tl(i, j, k)+0.5*(phiz_tl(i, k)+phiz_tl(i, k+1)&
&               +2*w(i, j, k)*w_tl(i, j, k)+0.5*rsin2_l(i, j)*(2*u(i, j&
&               , k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i&
&               , j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl(i+1, j, k)-&
&               cosa_s_l(i, j)*((u_tl(i, j, k)+u_tl(i, j+1, k))*(v(i, j&
&               , k)+v(i+1, j, k))+(u(i, j, k)+u(i, j+1, k))*(v_tl(i, j&
&               , k)+v_tl(i+1, j, k))))))
              te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv_air*pt(i, j&
&               , k)+0.5*(phiz(i, k)+phiz(i, k+1)+w(i, j, k)**2+0.5*&
&               rsin2_l(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)&
&               **2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k&
&               )+v(i+1, j, k))*cosa_s_l(i, j))))
            END DO
          END DO
        ELSE
          DO k=1,km
            DO i=is,ie
              te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(cv_air&
&               *pt(i, j, k)+0.5*(phiz(i, k)+phiz(i, k+1)+w(i, j, k)**2+&
&               0.5*rsin2_l(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j&
&               , k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, &
&               j, k)+v(i+1, j, k))*cosa_s_l(i, j)))) + delp(i, j, k)*(&
&               cv_air*pt_tl(i, j, k)+0.5*(phiz_tl(i, k)+phiz_tl(i, k+1)&
&               +2*w(i, j, k)*w_tl(i, j, k)+0.5*rsin2_l(i, j)*(2*u(i, j&
&               , k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i&
&               , j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl(i+1, j, k)-&
&               cosa_s_l(i, j)*((u_tl(i, j, k)+u_tl(i, j+1, k))*(v(i, j&
&               , k)+v(i+1, j, k))+(u(i, j, k)+u(i, j+1, k))*(v_tl(i, j&
&               , k)+v_tl(i+1, j, k))))))
              te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv_air*pt(i, j&
&               , k)+0.5*(phiz(i, k)+phiz(i, k+1)+w(i, j, k)**2+0.5*&
&               rsin2_l(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)&
&               **2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k&
&               )+v(i+1, j, k))*cosa_s_l(i, j))))
            END DO
          END DO
        END IF
      END IF
    END DO
!-------------------------------------
! Diganostics computation for moist TE
!-------------------------------------
    IF (id_te .GT. 0) THEN
      teq_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,teq,te_2d,moist_phys,km,hlv,sphum,q,delp)
      DO j=js,je
        DO i=is,ie
          teq_tl(i, j) = te_2d_tl(i, j)
          teq(i, j) = te_2d(i, j)
        END DO
        IF (moist_phys) THEN
          DO k=1,km
            DO i=is,ie
              teq_tl(i, j) = teq_tl(i, j) + hlv*(q_tl(i, j, k, sphum)*&
&               delp(i, j, k)+q(i, j, k, sphum)*delp_tl(i, j, k))
              teq(i, j) = teq(i, j) + hlv*q(i, j, k, sphum)*delp(i, j, k&
&               )
            END DO
          END DO
        END IF
      END DO
    ELSE
      teq_tl = 0.0
    END IF
  END SUBROUTINE COMPUTE_TOTAL_ENERGY_TLM
  SUBROUTINE COMPUTE_TOTAL_ENERGY(is, ie, js, je, isd, ied, jsd, jed, km&
&   , u, v, w, delz, pt, delp, q, qc, pe, peln, hs, rsin2_l, cosa_s_l, &
&   r_vir, cp, rg, hlv, te_2d, ua, va, teq, moist_phys, nwat, sphum, &
&   liq_wat, rainwat, ice_wat, snowwat, graupel, hydrostatic, id_te)
    IMPLICIT NONE
!------------------------------------------------------
! Compute vertically integrated total energy per column
!------------------------------------------------------
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: km, is, ie, js, je, isd, ied, jsd, jed, id_te
    INTEGER, INTENT(IN) :: sphum, liq_wat, ice_wat, rainwat, snowwat, &
&   graupel, nwat
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt, delp
    REAL, DIMENSION(isd:ied, jsd:jed, km, *), INTENT(IN) :: q
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: qc
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
! vertical velocity (m/s)
    REAL, INTENT(IN) :: w(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delz(isd:ied, jsd:jed, km)
! surface geopotential
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
! pressure at layer edges
    REAL, INTENT(IN) :: pe(is-1:ie+1, km+1, js-1:je+1)
! log(pe)
    REAL, INTENT(IN) :: peln(is:ie, km+1, js:je)
    REAL, INTENT(IN) :: cp, rg, r_vir, hlv
    REAL, INTENT(IN) :: rsin2_l(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: cosa_s_l(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: moist_phys, hydrostatic
! Output:
! vertically integrated TE
    REAL, INTENT(OUT) :: te_2d(is:ie, js:je)
! Moist TE
    REAL, INTENT(OUT) :: teq(is:ie, js:je)
! Local
    REAL, DIMENSION(is:ie, km) :: tv
    REAL :: phiz(is:ie, km+1)
    REAL :: cvm(is:ie), qd(is:ie)
    INTEGER :: i, j, k
!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, flagstruct%c2l_ord)
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,km,hydrostatic,hs,pt,qc,rg,peln,te_2d, &
!$OMP                                  pe,delp,cp,rsin2_l,u,v,cosa_s_l,delz,moist_phys,w, &
!$OMP                                  q,nwat,liq_wat,rainwat,ice_wat,snowwat,graupel,sphum)   &
!$OMP                          private(phiz, tv, cvm, qd)
    DO j=js,je
      IF (hydrostatic) THEN
        DO i=is,ie
          phiz(i, km+1) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            tv(i, k) = pt(i, j, k)*(1.+qc(i, j, k))
            phiz(i, k) = phiz(i, k+1) + rg*tv(i, k)*(peln(i, k+1, j)-&
&             peln(i, k, j))
          END DO
        END DO
        DO i=is,ie
          te_2d(i, j) = pe(i, km+1, j)*phiz(i, km+1) - pe(i, 1, j)*phiz(&
&           i, 1)
        END DO
        DO k=1,km
          DO i=is,ie
            te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cp*tv(i, k)+0.25*&
&             rsin2_l(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2&
&             +v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i&
&             +1, j, k))*cosa_s_l(i, j)))
          END DO
        END DO
      ELSE
!-----------------
! Non-hydrostatic:
!-----------------
        DO i=is,ie
          phiz(i, km+1) = hs(i, j)
          DO k=km,1,-1
            phiz(i, k) = phiz(i, k+1) - grav*delz(i, j, k)
          END DO
        END DO
        DO i=is,ie
          te_2d(i, j) = 0.
        END DO
        IF (moist_phys) THEN
          DO k=1,km
            DO i=is,ie
              te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv_air*pt(i, j&
&               , k)+0.5*(phiz(i, k)+phiz(i, k+1)+w(i, j, k)**2+0.5*&
&               rsin2_l(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)&
&               **2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k&
&               )+v(i+1, j, k))*cosa_s_l(i, j))))
            END DO
          END DO
        ELSE
          DO k=1,km
            DO i=is,ie
              te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv_air*pt(i, j&
&               , k)+0.5*(phiz(i, k)+phiz(i, k+1)+w(i, j, k)**2+0.5*&
&               rsin2_l(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)&
&               **2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k&
&               )+v(i+1, j, k))*cosa_s_l(i, j))))
            END DO
          END DO
        END IF
      END IF
    END DO
!-------------------------------------
! Diganostics computation for moist TE
!-------------------------------------
    IF (id_te .GT. 0) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,teq,te_2d,moist_phys,km,hlv,sphum,q,delp)
      DO j=js,je
        DO i=is,ie
          teq(i, j) = te_2d(i, j)
        END DO
        IF (moist_phys) THEN
          DO k=1,km
            DO i=is,ie
              teq(i, j) = teq(i, j) + hlv*q(i, j, k, sphum)*delp(i, j, k&
&               )
            END DO
          END DO
        END IF
      END DO
    END IF
  END SUBROUTINE COMPUTE_TOTAL_ENERGY
!  Differentiation of pkez in forward (tangent) mode:
!   variations   of useful results: peln pkz
!   with respect to varying inputs: peln pkz pk
  SUBROUTINE PKEZ_TLM(km, ifirst, ilast, jfirst, jlast, j, pe, pk, pk_tl&
&   , akap, peln, peln_tl, pkz, pkz_tl, ptop)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: km, j
! Latitude strip
    INTEGER, INTENT(IN) :: ifirst, ilast
! Latitude strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    REAL, INTENT(IN) :: akap
    REAL, INTENT(IN) :: pe(ifirst-1:ilast+1, km+1, jfirst-1:jlast+1)
    REAL, INTENT(IN) :: pk(ifirst:ilast, jfirst:jlast, km+1)
    REAL, INTENT(IN) :: pk_tl(ifirst:ilast, jfirst:jlast, km+1)
    REAL, INTENT(IN) :: ptop
! !OUTPUT
    REAL, INTENT(OUT) :: pkz(ifirst:ilast, jfirst:jlast, km)
    REAL, INTENT(OUT) :: pkz_tl(ifirst:ilast, jfirst:jlast, km)
! log (pe)
    REAL, INTENT(INOUT) :: peln(ifirst:ilast, km+1, jfirst:jlast)
    REAL, INTENT(INOUT) :: peln_tl(ifirst:ilast, km+1, jfirst:jlast)
! Local
    REAL :: pk2(ifirst:ilast, km+1)
    REAL :: pk2_tl(ifirst:ilast, km+1)
    REAL :: pek
    REAL :: pek_tl
    REAL :: lnp
    REAL :: ak1
    INTEGER :: i, k
    INTRINSIC LOG
    ak1 = (akap+1.)/akap
    pek_tl = pk_tl(ifirst, j, 1)
    pek = pk(ifirst, j, 1)
    pk2_tl = 0.0
    DO i=ifirst,ilast
      pk2_tl(i, 1) = pek_tl
      pk2(i, 1) = pek
    END DO
    DO k=2,km+1
      DO i=ifirst,ilast
!             peln(i,k,j) =  log(pe(i,k,j))
        pk2_tl(i, k) = pk_tl(i, j, k)
        pk2(i, k) = pk(i, j, k)
      END DO
    END DO
!---- GFDL modification
    IF (ptop .LT. ptop_min) THEN
      DO i=ifirst,ilast
        peln_tl(i, 1, j) = peln_tl(i, 2, j)
        peln(i, 1, j) = peln(i, 2, j) - ak1
      END DO
    ELSE
      lnp = LOG(ptop)
      DO i=ifirst,ilast
        peln_tl(i, 1, j) = 0.0
        peln(i, 1, j) = lnp
      END DO
    END IF
!---- GFDL modification
    DO k=1,km
      DO i=ifirst,ilast
        pkz_tl(i, j, k) = ((pk2_tl(i, k+1)-pk2_tl(i, k))*akap*(peln(i, k&
&         +1, j)-peln(i, k, j))-(pk2(i, k+1)-pk2(i, k))*akap*(peln_tl(i&
&         , k+1, j)-peln_tl(i, k, j)))/(akap*(peln(i, k+1, j)-peln(i, k&
&         , j)))**2
        pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1, j)-&
&         peln(i, k, j)))
      END DO
    END DO
  END SUBROUTINE PKEZ_TLM
  SUBROUTINE PKEZ(km, ifirst, ilast, jfirst, jlast, j, pe, pk, akap, &
&   peln, pkz, ptop)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: km, j
! Latitude strip
    INTEGER, INTENT(IN) :: ifirst, ilast
! Latitude strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    REAL, INTENT(IN) :: akap
    REAL, INTENT(IN) :: pe(ifirst-1:ilast+1, km+1, jfirst-1:jlast+1)
    REAL, INTENT(IN) :: pk(ifirst:ilast, jfirst:jlast, km+1)
    REAL, INTENT(IN) :: ptop
! !OUTPUT
    REAL, INTENT(OUT) :: pkz(ifirst:ilast, jfirst:jlast, km)
! log (pe)
    REAL, INTENT(INOUT) :: peln(ifirst:ilast, km+1, jfirst:jlast)
! Local
    REAL :: pk2(ifirst:ilast, km+1)
    REAL :: pek
    REAL :: lnp
    REAL :: ak1
    INTEGER :: i, k
    INTRINSIC LOG
    ak1 = (akap+1.)/akap
    pek = pk(ifirst, j, 1)
    DO i=ifirst,ilast
      pk2(i, 1) = pek
    END DO
    DO k=2,km+1
      DO i=ifirst,ilast
!             peln(i,k,j) =  log(pe(i,k,j))
        pk2(i, k) = pk(i, j, k)
      END DO
    END DO
!---- GFDL modification
    IF (ptop .LT. ptop_min) THEN
      DO i=ifirst,ilast
        peln(i, 1, j) = peln(i, 2, j) - ak1
      END DO
    ELSE
      lnp = LOG(ptop)
      DO i=ifirst,ilast
        peln(i, 1, j) = lnp
      END DO
    END IF
!---- GFDL modification
    DO k=1,km
      DO i=ifirst,ilast
        pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1, j)-&
&         peln(i, k, j)))
      END DO
    END DO
  END SUBROUTINE PKEZ
  SUBROUTINE REMAP_Z(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Method order
    INTEGER, INTENT(IN) :: kord
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
    INTEGER, INTENT(IN) :: iv
! height at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
! (from model top to bottom surface)
! hieght at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
! (from model top to bottom surface)
! Field input
    REAL, INTENT(IN) :: q1(i1:i2, km)
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(i1:i2, kn)
! !LOCAL VARIABLES:
    REAL :: qs(i1:i2)
    REAL :: dp1(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: pl, pr, qsum, delp, esl
    INTEGER :: i, k, l, m, k0
    DO k=1,km
      DO i=i1,i2
! negative
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q1(i, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE(qs, q4, dp1, km, i1, i2, iv, kord)
    ELSE
      CALL PPM_PROFILE(q4, dp1, km, i1, i2, iv, kord)
    END IF
! Mapping
    DO i=i1,i2
      k0 = 1
      DO k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .LE. pe1(i, l) .AND. pe2(i, k) .GE. pe1(i, l+1)&
&         ) THEN
            pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
            IF (pe2(i, k+1) .GE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
              pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
              q2(i, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2&
&               , i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
              k0 = l
              GOTO 555
            ELSE
! Fractional area...
              qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, &
&               l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+&
&               pl*(1.+pl))))
              DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                IF (pe2(i, k+1) .LT. pe1(i, m+1)) THEN
! Whole layer..
                  qsum = qsum + dp1(i, m)*q4(1, i, m)
                ELSE
                  delp = pe2(i, k+1) - pe1(i, m)
                  esl = delp/dp1(i, m)
                  qsum = qsum + delp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-&
&                   q4(2, i, m)+q4(4, i, m)*(1.-r23*esl)))
                  k0 = m
                  GOTO 123
                END IF
              END DO
              GOTO 123
            END IF
          END IF
        END DO
 123    q2(i, k) = qsum/(pe2(i, k+1)-pe2(i, k))
 555    CONTINUE
      END DO
    END DO
  END SUBROUTINE REMAP_Z
  SUBROUTINE MAP_SCALAR(km, pe1, qs, kn, pe2, q2, i1, i2, j, ibeg, iend&
&   , jbeg, jend, iv, kord, q_min)
    IMPLICIT NONE
! iv=1
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Mode: 0 == constituents  1 == temp
    INTEGER, INTENT(IN) :: iv
!       2 == remap temp with cs scheme
! Method order
    INTEGER, INTENT(IN) :: kord
! Current latitude
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! bottom BC
    REAL, INTENT(IN) :: qs(i1:i2)
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(IN) :: q_min
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL :: dp1(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    INTEGER :: i, k, l, m, k0
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL SCALAR_PROFILE(qs, q4, dp1, km, i1, i2, iv, kord, q_min)
    ELSE
      CALL PPM_PROFILE(q4, dp1, km, i1, i2, iv, kord)
    END IF
    DO i=i1,i2
      k0 = 1
      DO k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
            IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
              pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
              q2(i, j, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-&
&               q4(2, i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
              k0 = l
              GOTO 555
            ELSE
! Fractional area...
              qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, &
&               l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+&
&               pl*(1.+pl))))
              DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
                  qsum = qsum + dp1(i, m)*q4(1, i, m)
                ELSE
                  dp = pe2(i, k+1) - pe1(i, m)
                  esl = dp/dp1(i, m)
                  qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(&
&                   2, i, m)+q4(4, i, m)*(1.-r23*esl)))
                  k0 = m
                  GOTO 123
                END IF
              END DO
              GOTO 123
            END IF
          END IF
        END DO
 123    q2(i, j, k) = qsum/(pe2(i, k+1)-pe2(i, k))
 555    CONTINUE
      END DO
    END DO
  END SUBROUTINE MAP_SCALAR
  SUBROUTINE MAP1_PPM(km, pe1, qs, kn, pe2, q2, i1, i2, j, ibeg, iend, &
&   jbeg, jend, iv, kord)
    IMPLICIT NONE
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Mode: 0 == constituents  1 == ???
    INTEGER, INTENT(IN) :: iv
!       2 == remap temp with cs scheme
! Method order
    INTEGER, INTENT(IN) :: kord
! Current latitude
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! bottom BC
    REAL, INTENT(IN) :: qs(i1:i2)
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL :: dp1(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    INTEGER :: i, k, l, m, k0
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE(qs, q4, dp1, km, i1, i2, iv, kord)
    ELSE
      CALL PPM_PROFILE(q4, dp1, km, i1, i2, iv, kord)
    END IF
    DO i=i1,i2
      k0 = 1
      DO k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
            IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
              pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
              q2(i, j, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-&
&               q4(2, i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
              k0 = l
              GOTO 555
            ELSE
! Fractional area...
              qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, &
&               l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+&
&               pl*(1.+pl))))
              DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
                  qsum = qsum + dp1(i, m)*q4(1, i, m)
                ELSE
                  dp = pe2(i, k+1) - pe1(i, m)
                  esl = dp/dp1(i, m)
                  qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(&
&                   2, i, m)+q4(4, i, m)*(1.-r23*esl)))
                  k0 = m
                  GOTO 123
                END IF
              END DO
              GOTO 123
            END IF
          END IF
        END DO
 123    q2(i, j, k) = qsum/(pe2(i, k+1)-pe2(i, k))
 555    CONTINUE
      END DO
    END DO
  END SUBROUTINE MAP1_PPM
  SUBROUTINE MAPN_TRACER(nq, km, pe1, pe2, q1, dp2, kord, j, i1, i2, isd&
&   , ied, jsd, jed, q_min, fill)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! vertical dimension
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: j, nq, i1, i2
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: kord(nq)
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
    REAL, INTENT(IN) :: dp2(i1:i2, km)
    REAL, INTENT(IN) :: q_min
    LOGICAL, INTENT(IN) :: fill
! Field input
    REAL, INTENT(INOUT) :: q1(isd:ied, jsd:jed, km, nq)
! !LOCAL VARIABLES:
    REAL :: q4(4, i1:i2, km, nq)
! Field output
    REAL :: q2(i1:i2, km, nq)
    REAL :: qsum(nq)
    REAL :: dp1(i1:i2, km)
    REAL :: qs(i1:i2)
    REAL :: pl, pr, dp, esl, fac1, fac2
    INTEGER :: i, k, l, m, k0, iq
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
      END DO
    END DO
    DO iq=1,nq
      DO k=1,km
        DO i=i1,i2
          q4(1, i, k, iq) = q1(i, j, k, iq)
        END DO
      END DO
      CALL SCALAR_PROFILE(qs, q4(1:4, i1:i2, 1:km, iq), dp1, km, i1, i2&
&                   , 0, kord(iq), q_min)
    END DO
! Mapping
    DO i=i1,i2
      k0 = 1
      DO k=1,km
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
            IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
              pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
              fac1 = pr + pl
              fac2 = r3*(pr*fac1+pl*pl)
              fac1 = 0.5*fac1
              DO iq=1,nq
                q2(i, k, iq) = q4(2, i, l, iq) + (q4(4, i, l, iq)+q4(3, &
&                 i, l, iq)-q4(2, i, l, iq))*fac1 - q4(4, i, l, iq)*fac2
              END DO
              k0 = l
              GOTO 555
            ELSE
! Fractional area...
              dp = pe1(i, l+1) - pe2(i, k)
              fac1 = 1. + pl
              fac2 = r3*(1.+pl*fac1)
              fac1 = 0.5*fac1
              DO iq=1,nq
                qsum(iq) = dp*(q4(2, i, l, iq)+(q4(4, i, l, iq)+q4(3, i&
&                 , l, iq)-q4(2, i, l, iq))*fac1-q4(4, i, l, iq)*fac2)
              END DO
              DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
                  DO iq=1,nq
                    qsum(iq) = qsum(iq) + dp1(i, m)*q4(1, i, m, iq)
                  END DO
                ELSE
                  dp = pe2(i, k+1) - pe1(i, m)
                  esl = dp/dp1(i, m)
                  fac1 = 0.5*esl
                  fac2 = 1. - r23*esl
                  DO iq=1,nq
                    qsum(iq) = qsum(iq) + dp*(q4(2, i, m, iq)+fac1*(q4(3&
&                     , i, m, iq)-q4(2, i, m, iq)+q4(4, i, m, iq)*fac2))
                  END DO
                  k0 = m
                  GOTO 123
                END IF
              END DO
              GOTO 123
            END IF
          END IF
        END DO
 123    CONTINUE
        DO iq=1,nq
          q2(i, k, iq) = qsum(iq)/dp2(i, k)
        END DO
 555    CONTINUE
      END DO
    END DO
    IF (fill) CALL FILLZ(i2 - i1 + 1, km, nq, q2, dp2)
    DO iq=1,nq
!    if (fill) call fillz(i2-i1+1, km, 1, q2(i1,1,iq), dp2)
      DO k=1,km
        DO i=i1,i2
          q1(i, j, k, iq) = q2(i, k, iq)
        END DO
      END DO
    END DO
  END SUBROUTINE MAPN_TRACER
  SUBROUTINE MAP1_Q2(km, pe1, q1, kn, pe2, q2, dp2, i1, i2, iv, kord, j&
&   , ibeg, iend, jbeg, jend, q_min)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: i1, i2
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Mode: 0 ==  constituents 1 == ???
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! Field input
    REAL, INTENT(IN) :: q1(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(IN) :: dp2(i1:i2, kn)
    REAL, INTENT(IN) :: q_min
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(i1:i2, kn)
! !LOCAL VARIABLES:
    REAL :: qs(i1:i2)
    REAL :: dp1(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    INTEGER :: i, k, l, m, k0
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q1(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL SCALAR_PROFILE(qs, q4, dp1, km, i1, i2, iv, kord, q_min)
    ELSE
      CALL PPM_PROFILE(q4, dp1, km, i1, i2, iv, kord)
    END IF
! Mapping
    DO i=i1,i2
      k0 = 1
      DO k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
            IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
              pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
              q2(i, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2&
&               , i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
              k0 = l
              GOTO 555
            ELSE
! Fractional area...
              qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, &
&               l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+&
&               pl*(1.+pl))))
              DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
                  qsum = qsum + dp1(i, m)*q4(1, i, m)
                ELSE
                  dp = pe2(i, k+1) - pe1(i, m)
                  esl = dp/dp1(i, m)
                  qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(&
&                   2, i, m)+q4(4, i, m)*(1.-r23*esl)))
                  k0 = m
                  GOTO 123
                END IF
              END DO
              GOTO 123
            END IF
          END IF
        END DO
 123    q2(i, k) = qsum/dp2(i, k)
 555    CONTINUE
      END DO
    END DO
  END SUBROUTINE MAP1_Q2
  SUBROUTINE SCALAR_PROFILE(qs, a4, delp, km, i1, i2, iv, kord, qmin)
    IMPLICIT NONE
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
    INTEGER, INTENT(IN) :: i1, i2
! vertical dimension
    INTEGER, INTENT(IN) :: km
! iv =-1: winds
    INTEGER, INTENT(IN) :: iv
! iv = 0: positive definite scalars
! iv = 1: others
    INTEGER, INTENT(IN) :: kord
    REAL, INTENT(IN) :: qs(i1:i2)
! layer pressure thickness
    REAL, INTENT(IN) :: delp(i1:i2, km)
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(IN) :: qmin
!-----------------------------------------------------------------------
    LOGICAL, DIMENSION(i1:i2, km) :: extm, ext6
    REAL :: gam(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    INTEGER :: i, k, im
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: abs0
    INTEGER :: abs1
    REAL :: abs2
    INTEGER :: abs3
    INTEGER :: abs4
    REAL :: abs5
    INTEGER :: abs6
    REAL :: abs7
    INTEGER :: abs8
    REAL :: abs9
    INTEGER :: abs10
    INTEGER :: abs11
    INTEGER :: abs12
    REAL :: abs13
    REAL :: abs14
    REAL :: abs15
    REAL :: abs16
    REAL :: x12
    REAL :: x11
    REAL :: y29
    REAL :: x10
    REAL :: y28
    REAL :: y27
    REAL :: y26
    REAL :: y25
    REAL :: y24
    REAL :: y23
    REAL :: y22
    REAL :: y21
    REAL :: y20
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: x1
    REAL :: y19
    REAL :: y18
    REAL :: y17
    REAL :: y16
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
    REAL :: y32
    REAL :: y31
    REAL :: y30
    REAL :: y9
    REAL :: y8
    REAL :: y7
    REAL :: y6
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: y2
    REAL :: y1
    IF (iv .EQ. -2) THEN
      DO i=i1,i2
        gam(i, 2) = 0.5
        q(i, 1) = 1.5*a4(1, i, 1)
      END DO
      DO k=2,km-1
        DO i=i1,i2
          grat = delp(i, k-1)/delp(i, k)
          bet = 2. + grat + grat - gam(i, k)
          q(i, k) = (3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))/bet
          gam(i, k+1) = grat/bet
        END DO
      END DO
      DO i=i1,i2
        grat = delp(i, km-1)/delp(i, km)
        q(i, km) = (3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-&
&         1))/(2.+grat+grat-gam(i, km))
        q(i, km+1) = qs(i)
      END DO
      DO k=km-1,1,-1
        DO i=i1,i2
          q(i, k) = q(i, k) - gam(i, k+1)*q(i, k+1)
        END DO
      END DO
    ELSE
      DO i=i1,i2
! grid ratio
        grat = delp(i, 2)/delp(i, 1)
        bet = grat*(grat+0.5)
        q(i, 1) = ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))/bet
        gam(i, 1) = (1.+grat*(grat+1.5))/bet
      END DO
      DO k=2,km
        DO i=i1,i2
          d4(i) = delp(i, k-1)/delp(i, k)
          bet = 2. + d4(i) + d4(i) - gam(i, k-1)
          q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
          gam(i, k) = d4(i)/bet
        END DO
      END DO
      DO i=i1,i2
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
        END DO
      END DO
    END IF
    IF (kord .GE. 0.) THEN
      abs0 = kord
    ELSE
      abs0 = -kord
    END IF
!----- Perfectly linear scheme --------------------------------
    IF (abs0 .GT. 16) THEN
      DO k=1,km
        DO i=i1,i2
          a4(2, i, k) = q(i, k)
          a4(3, i, k) = q(i, k+1)
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
      END DO
      RETURN
    ELSE
!----- Perfectly linear scheme --------------------------------
!------------------
! Apply constraints
!------------------
      im = i2 - i1 + 1
! Apply *large-scale* constraints
      DO i=i1,i2
        IF (a4(1, i, 1) .LT. a4(1, i, 2)) THEN
          y1 = a4(1, i, 2)
        ELSE
          y1 = a4(1, i, 1)
        END IF
        IF (q(i, 2) .GT. y1) THEN
          q(i, 2) = y1
        ELSE
          q(i, 2) = q(i, 2)
        END IF
        IF (a4(1, i, 1) .GT. a4(1, i, 2)) THEN
          y2 = a4(1, i, 2)
        ELSE
          y2 = a4(1, i, 1)
        END IF
        IF (q(i, 2) .LT. y2) THEN
          q(i, 2) = y2
        ELSE
          q(i, 2) = q(i, 2)
        END IF
      END DO
      DO k=2,km
        DO i=i1,i2
          gam(i, k) = a4(1, i, k) - a4(1, i, k-1)
        END DO
      END DO
! Interior:
      DO k=3,km-1
        DO i=i1,i2
          IF (gam(i, k-1)*gam(i, k+1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y3 = a4(1, i, k)
            ELSE
              y3 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .GT. y3) THEN
              q(i, k) = y3
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y4 = a4(1, i, k)
            ELSE
              y4 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .LT. y4) THEN
              q(i, k) = y4
            ELSE
              q(i, k) = q(i, k)
            END IF
          ELSE IF (gam(i, k-1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y5 = a4(1, i, k)
            ELSE
              y5 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .LT. y5) THEN
              q(i, k) = y5
            ELSE
              q(i, k) = q(i, k)
            END IF
          ELSE
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y6 = a4(1, i, k)
            ELSE
              y6 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .GT. y6) THEN
              q(i, k) = y6
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (iv .EQ. 0) THEN
              IF (0. .LT. q(i, k)) THEN
                q(i, k) = q(i, k)
              ELSE
                q(i, k) = 0.
              END IF
            END IF
          END IF
        END DO
      END DO
! Bottom:
      DO i=i1,i2
        IF (a4(1, i, km-1) .LT. a4(1, i, km)) THEN
          y7 = a4(1, i, km)
        ELSE
          y7 = a4(1, i, km-1)
        END IF
        IF (q(i, km) .GT. y7) THEN
          q(i, km) = y7
        ELSE
          q(i, km) = q(i, km)
        END IF
        IF (a4(1, i, km-1) .GT. a4(1, i, km)) THEN
          y8 = a4(1, i, km)
        ELSE
          y8 = a4(1, i, km-1)
        END IF
        IF (q(i, km) .LT. y8) THEN
          q(i, km) = y8
        ELSE
          q(i, km) = q(i, km)
        END IF
      END DO
      DO k=1,km
        DO i=i1,i2
          a4(2, i, k) = q(i, k)
          a4(3, i, k) = q(i, k+1)
        END DO
      END DO
      DO k=1,km
        IF (k .EQ. 1 .OR. k .EQ. km) THEN
          DO i=i1,i2
            extm(i, k) = (a4(2, i, k)-a4(1, i, k))*(a4(3, i, k)-a4(1, i&
&             , k)) .GT. 0.
          END DO
        ELSE
          DO i=i1,i2
            extm(i, k) = gam(i, k)*gam(i, k+1) .LT. 0.
          END DO
        END IF
        IF (kord .GE. 0.) THEN
          abs1 = kord
        ELSE
          abs1 = -kord
        END IF
        IF (abs1 .EQ. 16) THEN
          DO i=i1,i2
            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
            IF (a4(4, i, k) .GE. 0.) THEN
              abs2 = a4(4, i, k)
            ELSE
              abs2 = -a4(4, i, k)
            END IF
            IF (a4(2, i, k) - a4(3, i, k) .GE. 0.) THEN
              abs13 = a4(2, i, k) - a4(3, i, k)
            ELSE
              abs13 = -(a4(2, i, k)-a4(3, i, k))
            END IF
            ext6(i, k) = abs2 .GT. abs13
          END DO
        END IF
      END DO
!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping
      IF (iv .EQ. 0) THEN
        DO i=i1,i2
          IF (0. .LT. a4(2, i, 1)) THEN
            a4(2, i, 1) = a4(2, i, 1)
          ELSE
            a4(2, i, 1) = 0.
          END IF
        END DO
      ELSE IF (iv .EQ. -1) THEN
        DO i=i1,i2
          IF (a4(2, i, 1)*a4(1, i, 1) .LE. 0.) a4(2, i, 1) = 0.
        END DO
      ELSE IF (iv .EQ. 2) THEN
        DO i=i1,i2
          a4(2, i, 1) = a4(1, i, 1)
          a4(3, i, 1) = a4(1, i, 1)
          a4(4, i, 1) = 0.
        END DO
      END IF
      IF (iv .NE. 2) THEN
        DO i=i1,i2
          a4(4, i, 1) = 3.*(2.*a4(1, i, 1)-(a4(2, i, 1)+a4(3, i, 1)))
        END DO
        CALL CS_LIMITERS(im, extm(i1, 1), a4(1, i1, 1), 1)
      END IF
! k=2
      DO i=i1,i2
        a4(4, i, 2) = 3.*(2.*a4(1, i, 2)-(a4(2, i, 2)+a4(3, i, 2)))
      END DO
      CALL CS_LIMITERS(im, extm(i1, 2), a4(1, i1, 2), 2)
!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
      DO k=3,km-2
        IF (kord .GE. 0.) THEN
          abs3 = kord
        ELSE
          abs3 = -kord
        END IF
        IF (abs3 .LT. 9) THEN
          DO i=i1,i2
! Left  edges
            pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
            lac_1 = pmp_1 + 1.5*gam(i, k+2)
            IF (a4(1, i, k) .GT. pmp_1) THEN
              IF (pmp_1 .GT. lac_1) THEN
                y21 = lac_1
              ELSE
                y21 = pmp_1
              END IF
            ELSE IF (a4(1, i, k) .GT. lac_1) THEN
              y21 = lac_1
            ELSE
              y21 = a4(1, i, k)
            END IF
            IF (a4(2, i, k) .LT. y21) THEN
              x1 = y21
            ELSE
              x1 = a4(2, i, k)
            END IF
            IF (a4(1, i, k) .LT. pmp_1) THEN
              IF (pmp_1 .LT. lac_1) THEN
                y9 = lac_1
              ELSE
                y9 = pmp_1
              END IF
            ELSE IF (a4(1, i, k) .LT. lac_1) THEN
              y9 = lac_1
            ELSE
              y9 = a4(1, i, k)
            END IF
            IF (x1 .GT. y9) THEN
              a4(2, i, k) = y9
            ELSE
              a4(2, i, k) = x1
            END IF
! Right edges
            pmp_2 = a4(1, i, k) + 2.*gam(i, k)
            lac_2 = pmp_2 - 1.5*gam(i, k-1)
            IF (a4(1, i, k) .GT. pmp_2) THEN
              IF (pmp_2 .GT. lac_2) THEN
                y22 = lac_2
              ELSE
                y22 = pmp_2
              END IF
            ELSE IF (a4(1, i, k) .GT. lac_2) THEN
              y22 = lac_2
            ELSE
              y22 = a4(1, i, k)
            END IF
            IF (a4(3, i, k) .LT. y22) THEN
              x2 = y22
            ELSE
              x2 = a4(3, i, k)
            END IF
            IF (a4(1, i, k) .LT. pmp_2) THEN
              IF (pmp_2 .LT. lac_2) THEN
                y10 = lac_2
              ELSE
                y10 = pmp_2
              END IF
            ELSE IF (a4(1, i, k) .LT. lac_2) THEN
              y10 = lac_2
            ELSE
              y10 = a4(1, i, k)
            END IF
            IF (x2 .GT. y10) THEN
              a4(3, i, k) = y10
            ELSE
              a4(3, i, k) = x2
            END IF
            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
          END DO
        ELSE
          IF (kord .GE. 0.) THEN
            abs4 = kord
          ELSE
            abs4 = -kord
          END IF
          IF (abs4 .EQ. 9) THEN
            DO i=i1,i2
              IF (extm(i, k) .AND. extm(i, k-1)) THEN
! grid-scale 2-delta-z wave detected
                a4(2, i, k) = a4(1, i, k)
                a4(3, i, k) = a4(1, i, k)
                a4(4, i, k) = 0.
              ELSE IF (extm(i, k) .AND. extm(i, k+1)) THEN
! grid-scale 2-delta-z wave detected
                a4(2, i, k) = a4(1, i, k)
                a4(3, i, k) = a4(1, i, k)
                a4(4, i, k) = 0.
              ELSE IF (extm(i, k) .AND. a4(1, i, k) .LT. qmin) THEN
! grid-scale 2-delta-z wave detected
                a4(2, i, k) = a4(1, i, k)
                a4(3, i, k) = a4(1, i, k)
                a4(4, i, k) = 0.
              ELSE
                a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k&
&                 )))
                IF (a4(4, i, k) .GE. 0.) THEN
                  abs5 = a4(4, i, k)
                ELSE
                  abs5 = -a4(4, i, k)
                END IF
                IF (a4(2, i, k) - a4(3, i, k) .GE. 0.) THEN
                  abs14 = a4(2, i, k) - a4(3, i, k)
                ELSE
                  abs14 = -(a4(2, i, k)-a4(3, i, k))
                END IF
! Check within the smooth region if subgrid profile is non-monotonic
                IF (abs5 .GT. abs14) THEN
                  pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                  lac_1 = pmp_1 + 1.5*gam(i, k+2)
                  IF (a4(1, i, k) .GT. pmp_1) THEN
                    IF (pmp_1 .GT. lac_1) THEN
                      y23 = lac_1
                    ELSE
                      y23 = pmp_1
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                    y23 = lac_1
                  ELSE
                    y23 = a4(1, i, k)
                  END IF
                  IF (a4(2, i, k) .LT. y23) THEN
                    x3 = y23
                  ELSE
                    x3 = a4(2, i, k)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_1) THEN
                    IF (pmp_1 .LT. lac_1) THEN
                      y11 = lac_1
                    ELSE
                      y11 = pmp_1
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                    y11 = lac_1
                  ELSE
                    y11 = a4(1, i, k)
                  END IF
                  IF (x3 .GT. y11) THEN
                    a4(2, i, k) = y11
                  ELSE
                    a4(2, i, k) = x3
                  END IF
                  pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                  lac_2 = pmp_2 - 1.5*gam(i, k-1)
                  IF (a4(1, i, k) .GT. pmp_2) THEN
                    IF (pmp_2 .GT. lac_2) THEN
                      y24 = lac_2
                    ELSE
                      y24 = pmp_2
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                    y24 = lac_2
                  ELSE
                    y24 = a4(1, i, k)
                  END IF
                  IF (a4(3, i, k) .LT. y24) THEN
                    x4 = y24
                  ELSE
                    x4 = a4(3, i, k)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_2) THEN
                    IF (pmp_2 .LT. lac_2) THEN
                      y12 = lac_2
                    ELSE
                      y12 = pmp_2
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                    y12 = lac_2
                  ELSE
                    y12 = a4(1, i, k)
                  END IF
                  IF (x4 .GT. y12) THEN
                    a4(3, i, k) = y12
                  ELSE
                    a4(3, i, k) = x4
                  END IF
                  a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i&
&                   , k)))
                END IF
              END IF
            END DO
          ELSE
            IF (kord .GE. 0.) THEN
              abs6 = kord
            ELSE
              abs6 = -kord
            END IF
            IF (abs6 .EQ. 10) THEN
              DO i=i1,i2
                IF (extm(i, k)) THEN
                  IF ((a4(1, i, k) .LT. qmin .OR. extm(i, k-1)) .OR. &
&                     extm(i, k+1)) THEN
! grid-scale 2-delta-z wave detected; or q is too small -> ehance vertical mixing
                    a4(2, i, k) = a4(1, i, k)
                    a4(3, i, k) = a4(1, i, k)
                    a4(4, i, k) = 0.
                  ELSE
! True local extremum
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                  END IF
                ELSE
! not a local extremum
                  a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i&
&                   , k))
                  IF (a4(4, i, k) .GE. 0.) THEN
                    abs7 = a4(4, i, k)
                  ELSE
                    abs7 = -a4(4, i, k)
                  END IF
                  IF (a4(2, i, k) - a4(3, i, k) .GE. 0.) THEN
                    abs15 = a4(2, i, k) - a4(3, i, k)
                  ELSE
                    abs15 = -(a4(2, i, k)-a4(3, i, k))
                  END IF
! Check within the smooth region if subgrid profile is non-monotonic
                  IF (abs7 .GT. abs15) THEN
                    pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                    lac_1 = pmp_1 + 1.5*gam(i, k+2)
                    IF (a4(1, i, k) .GT. pmp_1) THEN
                      IF (pmp_1 .GT. lac_1) THEN
                        y25 = lac_1
                      ELSE
                        y25 = pmp_1
                      END IF
                    ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                      y25 = lac_1
                    ELSE
                      y25 = a4(1, i, k)
                    END IF
                    IF (a4(2, i, k) .LT. y25) THEN
                      x5 = y25
                    ELSE
                      x5 = a4(2, i, k)
                    END IF
                    IF (a4(1, i, k) .LT. pmp_1) THEN
                      IF (pmp_1 .LT. lac_1) THEN
                        y13 = lac_1
                      ELSE
                        y13 = pmp_1
                      END IF
                    ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                      y13 = lac_1
                    ELSE
                      y13 = a4(1, i, k)
                    END IF
                    IF (x5 .GT. y13) THEN
                      a4(2, i, k) = y13
                    ELSE
                      a4(2, i, k) = x5
                    END IF
                    pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                    lac_2 = pmp_2 - 1.5*gam(i, k-1)
                    IF (a4(1, i, k) .GT. pmp_2) THEN
                      IF (pmp_2 .GT. lac_2) THEN
                        y26 = lac_2
                      ELSE
                        y26 = pmp_2
                      END IF
                    ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                      y26 = lac_2
                    ELSE
                      y26 = a4(1, i, k)
                    END IF
                    IF (a4(3, i, k) .LT. y26) THEN
                      x6 = y26
                    ELSE
                      x6 = a4(3, i, k)
                    END IF
                    IF (a4(1, i, k) .LT. pmp_2) THEN
                      IF (pmp_2 .LT. lac_2) THEN
                        y14 = lac_2
                      ELSE
                        y14 = pmp_2
                      END IF
                    ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                      y14 = lac_2
                    ELSE
                      y14 = a4(1, i, k)
                    END IF
                    IF (x6 .GT. y14) THEN
                      a4(3, i, k) = y14
                    ELSE
                      a4(3, i, k) = x6
                    END IF
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                  END IF
                END IF
              END DO
            ELSE
              IF (kord .GE. 0.) THEN
                abs8 = kord
              ELSE
                abs8 = -kord
              END IF
              IF (abs8 .EQ. 12) THEN
                DO i=i1,i2
                  IF (extm(i, k)) THEN
                    a4(2, i, k) = a4(1, i, k)
                    a4(3, i, k) = a4(1, i, k)
                    a4(4, i, k) = 0.
                  ELSE
! not a local extremum
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                    IF (a4(4, i, k) .GE. 0.) THEN
                      abs9 = a4(4, i, k)
                    ELSE
                      abs9 = -a4(4, i, k)
                    END IF
                    IF (a4(2, i, k) - a4(3, i, k) .GE. 0.) THEN
                      abs16 = a4(2, i, k) - a4(3, i, k)
                    ELSE
                      abs16 = -(a4(2, i, k)-a4(3, i, k))
                    END IF
! Check within the smooth region if subgrid profile is non-monotonic
                    IF (abs9 .GT. abs16) THEN
                      pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                      lac_1 = pmp_1 + 1.5*gam(i, k+2)
                      IF (a4(1, i, k) .GT. pmp_1) THEN
                        IF (pmp_1 .GT. lac_1) THEN
                          y27 = lac_1
                        ELSE
                          y27 = pmp_1
                        END IF
                      ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                        y27 = lac_1
                      ELSE
                        y27 = a4(1, i, k)
                      END IF
                      IF (a4(2, i, k) .LT. y27) THEN
                        x7 = y27
                      ELSE
                        x7 = a4(2, i, k)
                      END IF
                      IF (a4(1, i, k) .LT. pmp_1) THEN
                        IF (pmp_1 .LT. lac_1) THEN
                          y15 = lac_1
                        ELSE
                          y15 = pmp_1
                        END IF
                      ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                        y15 = lac_1
                      ELSE
                        y15 = a4(1, i, k)
                      END IF
                      IF (x7 .GT. y15) THEN
                        a4(2, i, k) = y15
                      ELSE
                        a4(2, i, k) = x7
                      END IF
                      pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                      lac_2 = pmp_2 - 1.5*gam(i, k-1)
                      IF (a4(1, i, k) .GT. pmp_2) THEN
                        IF (pmp_2 .GT. lac_2) THEN
                          y28 = lac_2
                        ELSE
                          y28 = pmp_2
                        END IF
                      ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                        y28 = lac_2
                      ELSE
                        y28 = a4(1, i, k)
                      END IF
                      IF (a4(3, i, k) .LT. y28) THEN
                        x8 = y28
                      ELSE
                        x8 = a4(3, i, k)
                      END IF
                      IF (a4(1, i, k) .LT. pmp_2) THEN
                        IF (pmp_2 .LT. lac_2) THEN
                          y16 = lac_2
                        ELSE
                          y16 = pmp_2
                        END IF
                      ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                        y16 = lac_2
                      ELSE
                        y16 = a4(1, i, k)
                      END IF
                      IF (x8 .GT. y16) THEN
                        a4(3, i, k) = y16
                      ELSE
                        a4(3, i, k) = x8
                      END IF
                      a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(&
&                       3, i, k))
                    END IF
                  END IF
                END DO
              ELSE
                IF (kord .GE. 0.) THEN
                  abs10 = kord
                ELSE
                  abs10 = -kord
                END IF
                IF (abs10 .EQ. 13) THEN
                  DO i=i1,i2
                    IF (extm(i, k)) THEN
                      IF (extm(i, k-1) .AND. extm(i, k+1)) THEN
! grid-scale 2-delta-z wave detected
                        a4(2, i, k) = a4(1, i, k)
                        a4(3, i, k) = a4(1, i, k)
                        a4(4, i, k) = 0.
                      ELSE
! Left  edges
                        pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                        lac_1 = pmp_1 + 1.5*gam(i, k+2)
                        IF (a4(1, i, k) .GT. pmp_1) THEN
                          IF (pmp_1 .GT. lac_1) THEN
                            y29 = lac_1
                          ELSE
                            y29 = pmp_1
                          END IF
                        ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                          y29 = lac_1
                        ELSE
                          y29 = a4(1, i, k)
                        END IF
                        IF (a4(2, i, k) .LT. y29) THEN
                          x9 = y29
                        ELSE
                          x9 = a4(2, i, k)
                        END IF
                        IF (a4(1, i, k) .LT. pmp_1) THEN
                          IF (pmp_1 .LT. lac_1) THEN
                            y17 = lac_1
                          ELSE
                            y17 = pmp_1
                          END IF
                        ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                          y17 = lac_1
                        ELSE
                          y17 = a4(1, i, k)
                        END IF
                        IF (x9 .GT. y17) THEN
                          a4(2, i, k) = y17
                        ELSE
                          a4(2, i, k) = x9
                        END IF
! Right edges
                        pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                        lac_2 = pmp_2 - 1.5*gam(i, k-1)
                        IF (a4(1, i, k) .GT. pmp_2) THEN
                          IF (pmp_2 .GT. lac_2) THEN
                            y30 = lac_2
                          ELSE
                            y30 = pmp_2
                          END IF
                        ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                          y30 = lac_2
                        ELSE
                          y30 = a4(1, i, k)
                        END IF
                        IF (a4(3, i, k) .LT. y30) THEN
                          x10 = y30
                        ELSE
                          x10 = a4(3, i, k)
                        END IF
                        IF (a4(1, i, k) .LT. pmp_2) THEN
                          IF (pmp_2 .LT. lac_2) THEN
                            y18 = lac_2
                          ELSE
                            y18 = pmp_2
                          END IF
                        ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                          y18 = lac_2
                        ELSE
                          y18 = a4(1, i, k)
                        END IF
                        IF (x10 .GT. y18) THEN
                          a4(3, i, k) = y18
                        ELSE
                          a4(3, i, k) = x10
                        END IF
                        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4&
&                         (3, i, k)))
                      END IF
                    ELSE
                      a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3&
&                       , i, k)))
                    END IF
                  END DO
                ELSE
                  IF (kord .GE. 0.) THEN
                    abs11 = kord
                  ELSE
                    abs11 = -kord
                  END IF
                  IF (abs11 .EQ. 14) THEN
                    DO i=i1,i2
                      a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3&
&                       , i, k)))
                    END DO
                  ELSE
                    IF (kord .GE. 0.) THEN
                      abs12 = kord
                    ELSE
                      abs12 = -kord
                    END IF
                    IF (abs12 .EQ. 16) THEN
                      DO i=i1,i2
                        IF (ext6(i, k)) THEN
                          IF (extm(i, k-1) .OR. extm(i, k+1)) THEN
! Left  edges
                            pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                            lac_1 = pmp_1 + 1.5*gam(i, k+2)
                            IF (a4(1, i, k) .GT. pmp_1) THEN
                              IF (pmp_1 .GT. lac_1) THEN
                                y31 = lac_1
                              ELSE
                                y31 = pmp_1
                              END IF
                            ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                              y31 = lac_1
                            ELSE
                              y31 = a4(1, i, k)
                            END IF
                            IF (a4(2, i, k) .LT. y31) THEN
                              x11 = y31
                            ELSE
                              x11 = a4(2, i, k)
                            END IF
                            IF (a4(1, i, k) .LT. pmp_1) THEN
                              IF (pmp_1 .LT. lac_1) THEN
                                y19 = lac_1
                              ELSE
                                y19 = pmp_1
                              END IF
                            ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                              y19 = lac_1
                            ELSE
                              y19 = a4(1, i, k)
                            END IF
                            IF (x11 .GT. y19) THEN
                              a4(2, i, k) = y19
                            ELSE
                              a4(2, i, k) = x11
                            END IF
! Right edges
                            pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                            lac_2 = pmp_2 - 1.5*gam(i, k-1)
                            IF (a4(1, i, k) .GT. pmp_2) THEN
                              IF (pmp_2 .GT. lac_2) THEN
                                y32 = lac_2
                              ELSE
                                y32 = pmp_2
                              END IF
                            ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                              y32 = lac_2
                            ELSE
                              y32 = a4(1, i, k)
                            END IF
                            IF (a4(3, i, k) .LT. y32) THEN
                              x12 = y32
                            ELSE
                              x12 = a4(3, i, k)
                            END IF
                            IF (a4(1, i, k) .LT. pmp_2) THEN
                              IF (pmp_2 .LT. lac_2) THEN
                                y20 = lac_2
                              ELSE
                                y20 = pmp_2
                              END IF
                            ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                              y20 = lac_2
                            ELSE
                              y20 = a4(1, i, k)
                            END IF
                            IF (x12 .GT. y20) THEN
                              a4(3, i, k) = y20
                            ELSE
                              a4(3, i, k) = x12
                            END IF
                            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k&
&                             )+a4(3, i, k)))
                          END IF
                        END IF
                      END DO
                    ELSE
! kord = 11, 13
                      DO i=i1,i2
                        IF (extm(i, k) .AND. ((extm(i, k-1) .OR. extm(i&
&                           , k+1)) .OR. a4(1, i, k) .LT. qmin)) THEN
! Noisy region:
                          a4(2, i, k) = a4(1, i, k)
                          a4(3, i, k) = a4(1, i, k)
                          a4(4, i, k) = 0.
                        ELSE
                          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+&
&                           a4(3, i, k)))
                        END IF
                      END DO
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
! Additional constraint to ensure positivity
        IF (iv .EQ. 0) CALL CS_LIMITERS(im, extm(i1, k), a4(1, i1, k), 0&
&                                )
      END DO
! k-loop
!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
      IF (iv .EQ. 0) THEN
        DO i=i1,i2
          IF (0. .LT. a4(3, i, km)) THEN
            a4(3, i, km) = a4(3, i, km)
          ELSE
            a4(3, i, km) = 0.
          END IF
        END DO
      ELSE IF (iv .EQ. -1) THEN
        DO i=i1,i2
          IF (a4(3, i, km)*a4(1, i, km) .LE. 0.) a4(3, i, km) = 0.
        END DO
      END IF
      DO k=km-1,km
        DO i=i1,i2
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
        IF (k .EQ. km - 1) CALL CS_LIMITERS(im, extm(i1, k), a4(1, i1, k&
&                                     ), 2)
        IF (k .EQ. km) CALL CS_LIMITERS(im, extm(i1, k), a4(1, i1, k), 1&
&                                )
      END DO
    END IF
  END SUBROUTINE SCALAR_PROFILE
!  Differentiation of cs_limiters in forward (tangent) mode:
!   variations   of useful results: a4
!   with respect to varying inputs: a4
  SUBROUTINE CS_LIMITERS_TLM(im, extm, a4, a4_tl, iv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    INTEGER, INTENT(IN) :: iv
    LOGICAL, INTENT(IN) :: extm(im)
! PPM array
    REAL, INTENT(INOUT) :: a4(4, im)
    REAL, INTENT(INOUT) :: a4_tl(4, im)
! !LOCAL VARIABLES:
    REAL :: da1, da2, a6da
    INTEGER :: i
    INTRINSIC ABS
    REAL :: abs0
    IF (iv .EQ. 0) THEN
! Positive definite constraint
      DO i=1,im
        IF (a4(1, i) .LE. 0.) THEN
          a4_tl(2, i) = a4_tl(1, i)
          a4(2, i) = a4(1, i)
          a4_tl(3, i) = a4_tl(1, i)
          a4(3, i) = a4(1, i)
          a4_tl(4, i) = 0.0
          a4(4, i) = 0.
        ELSE
          IF (a4(3, i) - a4(2, i) .GE. 0.) THEN
            abs0 = a4(3, i) - a4(2, i)
          ELSE
            abs0 = -(a4(3, i)-a4(2, i))
          END IF
          IF (abs0 .LT. -a4(4, i)) THEN
            IF (a4(1, i) + 0.25*(a4(3, i)-a4(2, i))**2/a4(4, i) + a4(4, &
&               i)*r12 .LT. 0.) THEN
! local minimum is negative
              IF (a4(1, i) .LT. a4(3, i) .AND. a4(1, i) .LT. a4(2, i)) &
&             THEN
                a4_tl(3, i) = a4_tl(1, i)
                a4(3, i) = a4(1, i)
                a4_tl(2, i) = a4_tl(1, i)
                a4(2, i) = a4(1, i)
                a4_tl(4, i) = 0.0
                a4(4, i) = 0.
              ELSE IF (a4(3, i) .GT. a4(2, i)) THEN
                a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
                a4(4, i) = 3.*(a4(2, i)-a4(1, i))
                a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
                a4(3, i) = a4(2, i) - a4(4, i)
              ELSE
                a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
                a4(4, i) = 3.*(a4(3, i)-a4(1, i))
                a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
                a4(2, i) = a4(3, i) - a4(4, i)
              END IF
            END IF
          END IF
        END IF
      END DO
    ELSE IF (iv .EQ. 1) THEN
      DO i=1,im
        IF ((a4(1, i)-a4(2, i))*(a4(1, i)-a4(3, i)) .GE. 0.) THEN
          a4_tl(2, i) = a4_tl(1, i)
          a4(2, i) = a4(1, i)
          a4_tl(3, i) = a4_tl(1, i)
          a4(3, i) = a4(1, i)
          a4_tl(4, i) = 0.0
          a4(4, i) = 0.
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
            a4(3, i) = a4(2, i) - a4(4, i)
          ELSE IF (a6da .GT. da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
            a4(2, i) = a4(3, i) - a4(4, i)
          END IF
        END IF
      END DO
    ELSE
! Standard PPM constraint
      DO i=1,im
        IF (extm(i)) THEN
          a4_tl(2, i) = a4_tl(1, i)
          a4(2, i) = a4(1, i)
          a4_tl(3, i) = a4_tl(1, i)
          a4(3, i) = a4(1, i)
          a4_tl(4, i) = 0.0
          a4(4, i) = 0.
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
            a4(3, i) = a4(2, i) - a4(4, i)
          ELSE IF (a6da .GT. da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
            a4(2, i) = a4(3, i) - a4(4, i)
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE CS_LIMITERS_TLM
!  Differentiation of ppm_profile in forward (tangent) mode:
!   variations   of useful results: a4
!   with respect to varying inputs: delp a4
  SUBROUTINE PPM_PROFILE_TLM(a4, a4_tl, delp, delp_tl, km, i1, i2, iv, &
&   kord)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! iv =-1: winds
    INTEGER, INTENT(IN) :: iv
! iv = 0: positive definite scalars
! iv = 1: others
! iv = 2: w (iv=-2)
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! vertical dimension
    INTEGER, INTENT(IN) :: km
! Order (or more accurately method no.):
    INTEGER, INTENT(IN) :: kord
!
! layer pressure thickness
    REAL, INTENT(IN) :: delp(i1:i2, km)
    REAL, INTENT(IN) :: delp_tl(i1:i2, km)
! !INPUT/OUTPUT PARAMETERS:
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_tl(4, i1:i2, km)
! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
!
! !REVISION HISTORY:
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
    REAL :: dc(i1:i2, km)
    REAL :: dc_tl(i1:i2, km)
    REAL :: h2(i1:i2, km)
    REAL :: h2_tl(i1:i2, km)
    REAL :: delq(i1:i2, km)
    REAL :: delq_tl(i1:i2, km)
    REAL :: df2(i1:i2, km)
    REAL :: df2_tl(i1:i2, km)
    REAL :: d4(i1:i2, km)
    REAL :: d4_tl(i1:i2, km)
! local scalars:
    INTEGER :: i, k, km1, lmt, it
    REAL :: fac
    REAL :: a1, a2, c1, c2, c3, d1, d2
    REAL :: a1_tl, a2_tl, c1_tl, c2_tl, c3_tl, d1_tl, d2_tl
    REAL :: qm, dq, lac, qmp, pmp
    REAL :: qm_tl, dq_tl, lac_tl, qmp_tl, pmp_tl
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min1_tl
    INTEGER :: abs0
    REAL :: max1
    REAL :: max1_tl
    REAL :: min2
    REAL :: min2_tl
    REAL :: x1_tl
    REAL :: y1_tl
    REAL :: z1_tl
    REAL :: y2_tl
    REAL :: y3_tl
    REAL :: y4_tl
    REAL :: y5_tl
    REAL :: y8_tl
    REAL :: x2_tl
    REAL :: y6_tl
    REAL :: y9_tl
    REAL :: x3_tl
    REAL :: y7_tl
    REAL :: x3
    REAL :: x2
    REAL :: x1
    REAL :: z1
    REAL :: y9
    REAL :: y8
    REAL :: y7
    REAL :: y6
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: y2
    REAL :: y1
    km1 = km - 1
    it = i2 - i1 + 1
    delq_tl = 0.0
    d4_tl = 0.0
    DO k=2,km
      DO i=i1,i2
        delq_tl(i, k-1) = a4_tl(1, i, k) - a4_tl(1, i, k-1)
        delq(i, k-1) = a4(1, i, k) - a4(1, i, k-1)
        d4_tl(i, k) = delp_tl(i, k-1) + delp_tl(i, k)
        d4(i, k) = delp(i, k-1) + delp(i, k)
      END DO
    END DO
    df2_tl = 0.0
    dc_tl = 0.0
    DO k=2,km1
      DO i=i1,i2
        c1_tl = ((delp_tl(i, k-1)+0.5*delp_tl(i, k))*d4(i, k+1)-(delp(i&
&         , k-1)+0.5*delp(i, k))*d4_tl(i, k+1))/d4(i, k+1)**2
        c1 = (delp(i, k-1)+0.5*delp(i, k))/d4(i, k+1)
        c2_tl = ((delp_tl(i, k+1)+0.5*delp_tl(i, k))*d4(i, k)-(delp(i, k&
&         +1)+0.5*delp(i, k))*d4_tl(i, k))/d4(i, k)**2
        c2 = (delp(i, k+1)+0.5*delp(i, k))/d4(i, k)
        df2_tl(i, k) = ((delp_tl(i, k)*(c1*delq(i, k)+c2*delq(i, k-1))+&
&         delp(i, k)*(c1_tl*delq(i, k)+c1*delq_tl(i, k)+c2_tl*delq(i, k-&
&         1)+c2*delq_tl(i, k-1)))*(d4(i, k)+delp(i, k+1))-delp(i, k)*(c1&
&         *delq(i, k)+c2*delq(i, k-1))*(d4_tl(i, k)+delp_tl(i, k+1)))/(&
&         d4(i, k)+delp(i, k+1))**2
        df2(i, k) = delp(i, k)*(c1*delq(i, k)+c2*delq(i, k-1))/(d4(i, k)&
&         +delp(i, k+1))
        IF (df2(i, k) .GE. 0.) THEN
          x1_tl = df2_tl(i, k)
          x1 = df2(i, k)
        ELSE
          x1_tl = -df2_tl(i, k)
          x1 = -df2(i, k)
        END IF
        IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
          IF (a4(1, i, k) .LT. a4(1, i, k+1)) THEN
            max1_tl = a4_tl(1, i, k+1)
            max1 = a4(1, i, k+1)
          ELSE
            max1_tl = a4_tl(1, i, k)
            max1 = a4(1, i, k)
          END IF
        ELSE IF (a4(1, i, k-1) .LT. a4(1, i, k+1)) THEN
          max1_tl = a4_tl(1, i, k+1)
          max1 = a4(1, i, k+1)
        ELSE
          max1_tl = a4_tl(1, i, k-1)
          max1 = a4(1, i, k-1)
        END IF
        y1_tl = max1_tl - a4_tl(1, i, k)
        y1 = max1 - a4(1, i, k)
        IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
          IF (a4(1, i, k) .GT. a4(1, i, k+1)) THEN
            min2_tl = a4_tl(1, i, k+1)
            min2 = a4(1, i, k+1)
          ELSE
            min2_tl = a4_tl(1, i, k)
            min2 = a4(1, i, k)
          END IF
        ELSE IF (a4(1, i, k-1) .GT. a4(1, i, k+1)) THEN
          min2_tl = a4_tl(1, i, k+1)
          min2 = a4(1, i, k+1)
        ELSE
          min2_tl = a4_tl(1, i, k-1)
          min2 = a4(1, i, k-1)
        END IF
        z1_tl = a4_tl(1, i, k) - min2_tl
        z1 = a4(1, i, k) - min2
        IF (x1 .GT. y1) THEN
          IF (y1 .GT. z1) THEN
            min1_tl = z1_tl
            min1 = z1
          ELSE
            min1_tl = y1_tl
            min1 = y1
          END IF
        ELSE IF (x1 .GT. z1) THEN
          min1_tl = z1_tl
          min1 = z1
        ELSE
          min1_tl = x1_tl
          min1 = x1
        END IF
        dc_tl(i, k) = min1_tl*SIGN(1.d0, min1*df2(i, k))
        dc(i, k) = SIGN(min1, df2(i, k))
      END DO
    END DO
!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------
    DO k=3,km1
      DO i=i1,i2
        c1_tl = ((delq_tl(i, k-1)*delp(i, k-1)+delq(i, k-1)*delp_tl(i, k&
&         -1))*d4(i, k)-delq(i, k-1)*delp(i, k-1)*d4_tl(i, k))/d4(i, k)&
&         **2
        c1 = delq(i, k-1)*delp(i, k-1)/d4(i, k)
        a1_tl = (d4_tl(i, k-1)*(d4(i, k)+delp(i, k-1))-d4(i, k-1)*(d4_tl&
&         (i, k)+delp_tl(i, k-1)))/(d4(i, k)+delp(i, k-1))**2
        a1 = d4(i, k-1)/(d4(i, k)+delp(i, k-1))
        a2_tl = (d4_tl(i, k+1)*(d4(i, k)+delp(i, k))-d4(i, k+1)*(d4_tl(i&
&         , k)+delp_tl(i, k)))/(d4(i, k)+delp(i, k))**2
        a2 = d4(i, k+1)/(d4(i, k)+delp(i, k))
        a4_tl(2, i, k) = a4_tl(1, i, k-1) + c1_tl + 2.*(delp_tl(i, k)*(&
&         c1*(a1-a2)+a2*dc(i, k-1))+delp(i, k)*(c1_tl*(a1-a2)+c1*(a1_tl-&
&         a2_tl)+a2_tl*dc(i, k-1)+a2*dc_tl(i, k-1))-delp_tl(i, k-1)*a1*&
&         dc(i, k)-delp(i, k-1)*(a1_tl*dc(i, k)+a1*dc_tl(i, k)))/(d4(i, &
&         k-1)+d4(i, k+1)) - 2.*(d4_tl(i, k-1)+d4_tl(i, k+1))*(delp(i, k&
&         )*(c1*(a1-a2)+a2*dc(i, k-1))-delp(i, k-1)*a1*dc(i, k))/(d4(i, &
&         k-1)+d4(i, k+1))**2
        a4(2, i, k) = a4(1, i, k-1) + c1 + 2./(d4(i, k-1)+d4(i, k+1))*(&
&         delp(i, k)*(c1*(a1-a2)+a2*dc(i, k-1))-delp(i, k-1)*a1*dc(i, k)&
&         )
      END DO
    END DO
!     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)
! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
    DO i=i1,i2
      d1_tl = delp_tl(i, 1)
      d1 = delp(i, 1)
      d2_tl = delp_tl(i, 2)
      d2 = delp(i, 2)
      qm_tl = ((d2_tl*a4(1, i, 1)+d2*a4_tl(1, i, 1)+d1_tl*a4(1, i, 2)+d1&
&       *a4_tl(1, i, 2))*(d1+d2)-(d2*a4(1, i, 1)+d1*a4(1, i, 2))*(d1_tl+&
&       d2_tl))/(d1+d2)**2
      qm = (d2*a4(1, i, 1)+d1*a4(1, i, 2))/(d1+d2)
      dq_tl = (2.*(a4_tl(1, i, 2)-a4_tl(1, i, 1))*(d1+d2)-2.*(a4(1, i, 2&
&       )-a4(1, i, 1))*(d1_tl+d2_tl))/(d1+d2)**2
      dq = 2.*(a4(1, i, 2)-a4(1, i, 1))/(d1+d2)
      c1_tl = (4.*(a4_tl(2, i, 3)-qm_tl-d2_tl*dq-d2*dq_tl)*d2*(2.*d2*d2+&
&       d1*(d2+3.*d1))-4.*(a4(2, i, 3)-qm-d2*dq)*(d2_tl*(2.*d2*d2+d1*(d2&
&       +3.*d1))+d2*(2.*(d2_tl*d2+d2*d2_tl)+d1_tl*(d2+3.*d1)+d1*(d2_tl+&
&       3.*d1_tl))))/(d2*(2.*d2*d2+d1*(d2+3.*d1)))**2
      c1 = 4.*(a4(2, i, 3)-qm-d2*dq)/(d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3_tl = dq_tl - 0.5*(c1_tl*(d2*(5.*d1+d2)-3.*d1*d1)+c1*(d2_tl*(5.*&
&       d1+d2)+d2*(5.*d1_tl+d2_tl)-3.*(d1_tl*d1+d1*d1_tl)))
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      a4_tl(2, i, 2) = qm_tl - 0.25*(((c1_tl*d1+c1*d1_tl)*d2+c1*d1*d2_tl&
&       )*(d2+3.*d1)+c1*d1*d2*(d2_tl+3.*d1_tl))
      a4(2, i, 2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
      a4_tl(2, i, 1) = d1_tl*(2.*c1*d1**2-c3) + d1*(2.*(c1_tl*d1**2+c1*2&
&       *d1*d1_tl)-c3_tl) + a4_tl(2, i, 2)
      a4(2, i, 1) = d1*(2.*c1*d1**2-c3) + a4(2, i, 2)
      IF (a4(1, i, 1) .GT. a4(1, i, 2)) THEN
        y2_tl = a4_tl(1, i, 2)
        y2 = a4(1, i, 2)
      ELSE
        y2_tl = a4_tl(1, i, 1)
        y2 = a4(1, i, 1)
      END IF
      IF (a4(2, i, 2) .LT. y2) THEN
        a4_tl(2, i, 2) = y2_tl
        a4(2, i, 2) = y2
      ELSE
        a4(2, i, 2) = a4(2, i, 2)
      END IF
      IF (a4(1, i, 1) .LT. a4(1, i, 2)) THEN
        y3_tl = a4_tl(1, i, 2)
        y3 = a4(1, i, 2)
      ELSE
        y3_tl = a4_tl(1, i, 1)
        y3 = a4(1, i, 1)
      END IF
      IF (a4(2, i, 2) .GT. y3) THEN
        a4_tl(2, i, 2) = y3_tl
        a4(2, i, 2) = y3
      ELSE
        a4(2, i, 2) = a4(2, i, 2)
      END IF
      dc_tl(i, 1) = 0.5*(a4_tl(2, i, 2)-a4_tl(1, i, 1))
      dc(i, 1) = 0.5*(a4(2, i, 2)-a4(1, i, 1))
    END DO
! Enforce monotonicity  within the top layer
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(2, i, 1)) THEN
          a4(2, i, 1) = a4(2, i, 1)
        ELSE
          a4_tl(2, i, 1) = 0.0
          a4(2, i, 1) = 0.
        END IF
        IF (0. .LT. a4(2, i, 2)) THEN
          a4(2, i, 2) = a4(2, i, 2)
        ELSE
          a4_tl(2, i, 2) = 0.0
          a4(2, i, 2) = 0.
        END IF
      END DO
    ELSE IF (iv .EQ. -1) THEN
      DO i=i1,i2
        IF (a4(2, i, 1)*a4(1, i, 1) .LE. 0.) THEN
          a4_tl(2, i, 1) = 0.0
          a4(2, i, 1) = 0.
        END IF
      END DO
    ELSE
      IF (iv .GE. 0.) THEN
        abs0 = iv
      ELSE
        abs0 = -iv
      END IF
      IF (abs0 .EQ. 2) THEN
        DO i=i1,i2
          a4_tl(2, i, 1) = a4_tl(1, i, 1)
          a4(2, i, 1) = a4(1, i, 1)
          a4_tl(3, i, 1) = a4_tl(1, i, 1)
          a4(3, i, 1) = a4(1, i, 1)
        END DO
      END IF
    END IF
! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
    DO i=i1,i2
      d1_tl = delp_tl(i, km)
      d1 = delp(i, km)
      d2_tl = delp_tl(i, km1)
      d2 = delp(i, km1)
      qm_tl = ((d2_tl*a4(1, i, km)+d2*a4_tl(1, i, km)+d1_tl*a4(1, i, km1&
&       )+d1*a4_tl(1, i, km1))*(d1+d2)-(d2*a4(1, i, km)+d1*a4(1, i, km1)&
&       )*(d1_tl+d2_tl))/(d1+d2)**2
      qm = (d2*a4(1, i, km)+d1*a4(1, i, km1))/(d1+d2)
      dq_tl = (2.*(a4_tl(1, i, km1)-a4_tl(1, i, km))*(d1+d2)-2.*(a4(1, i&
&       , km1)-a4(1, i, km))*(d1_tl+d2_tl))/(d1+d2)**2
      dq = 2.*(a4(1, i, km1)-a4(1, i, km))/(d1+d2)
      c1_tl = ((a4_tl(2, i, km1)-qm_tl-d2_tl*dq-d2*dq_tl)*d2*(2.*d2*d2+&
&       d1*(d2+3.*d1))-(a4(2, i, km1)-qm-d2*dq)*(d2_tl*(2.*d2*d2+d1*(d2+&
&       3.*d1))+d2*(2.*(d2_tl*d2+d2*d2_tl)+d1_tl*(d2+3.*d1)+d1*(d2_tl+3.&
&       *d1_tl))))/(d2*(2.*d2*d2+d1*(d2+3.*d1)))**2
      c1 = (a4(2, i, km1)-qm-d2*dq)/(d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3_tl = dq_tl - 2.0*(c1_tl*(d2*(5.*d1+d2)-3.*d1*d1)+c1*(d2_tl*(5.*&
&       d1+d2)+d2*(5.*d1_tl+d2_tl)-3.*(d1_tl*d1+d1*d1_tl)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      a4_tl(2, i, km) = qm_tl - ((c1_tl*d1+c1*d1_tl)*d2+c1*d1*d2_tl)*(d2&
&       +3.*d1) - c1*d1*d2*(d2_tl+3.*d1_tl)
      a4(2, i, km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
      a4_tl(3, i, km) = d1_tl*(8.*c1*d1**2-c3) + d1*(8.*(c1_tl*d1**2+c1*&
&       2*d1*d1_tl)-c3_tl) + a4_tl(2, i, km)
      a4(3, i, km) = d1*(8.*c1*d1**2-c3) + a4(2, i, km)
      IF (a4(1, i, km) .GT. a4(1, i, km1)) THEN
        y4_tl = a4_tl(1, i, km1)
        y4 = a4(1, i, km1)
      ELSE
        y4_tl = a4_tl(1, i, km)
        y4 = a4(1, i, km)
      END IF
      IF (a4(2, i, km) .LT. y4) THEN
        a4_tl(2, i, km) = y4_tl
        a4(2, i, km) = y4
      ELSE
        a4(2, i, km) = a4(2, i, km)
      END IF
      IF (a4(1, i, km) .LT. a4(1, i, km1)) THEN
        y5_tl = a4_tl(1, i, km1)
        y5 = a4(1, i, km1)
      ELSE
        y5_tl = a4_tl(1, i, km)
        y5 = a4(1, i, km)
      END IF
      IF (a4(2, i, km) .GT. y5) THEN
        a4_tl(2, i, km) = y5_tl
        a4(2, i, km) = y5
      ELSE
        a4(2, i, km) = a4(2, i, km)
      END IF
      dc_tl(i, km) = 0.5*(a4_tl(1, i, km)-a4_tl(2, i, km))
      dc(i, km) = 0.5*(a4(1, i, km)-a4(2, i, km))
    END DO
! Enforce constraint on the "slope" at the surface
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(2, i, km)) THEN
          a4(2, i, km) = a4(2, i, km)
        ELSE
          a4_tl(2, i, km) = 0.0
          a4(2, i, km) = 0.
        END IF
        IF (0. .LT. a4(3, i, km)) THEN
          a4(3, i, km) = a4(3, i, km)
        ELSE
          a4_tl(3, i, km) = 0.0
          a4(3, i, km) = 0.
        END IF
      END DO
    ELSE IF (iv .LT. 0) THEN
      DO i=i1,i2
        IF (a4(1, i, km)*a4(3, i, km) .LE. 0.) THEN
          a4_tl(3, i, km) = 0.0
          a4(3, i, km) = 0.
        END IF
      END DO
    END IF
    DO k=1,km1
      DO i=i1,i2
        a4_tl(3, i, k) = a4_tl(2, i, k+1)
        a4(3, i, k) = a4(2, i, k+1)
      END DO
    END DO
!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
    DO k=1,2
      DO i=i1,i2
        a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3, i&
&         , k))
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
      CALL PPM_LIMITERS_TLM(dc(i1, k), dc_tl(i1, k), a4(1, i1, k), a4_tl&
&                     (1, i1, k), it, 0)
    END DO
    IF (kord .GE. 7) THEN
      h2_tl = 0.0
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      DO k=2,km1
        DO i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
          h2_tl(i, k) = (2.*((dc_tl(i, k+1)*delp(i, k+1)-dc(i, k+1)*&
&           delp_tl(i, k+1))/delp(i, k+1)**2-(dc_tl(i, k-1)*delp(i, k-1)&
&           -dc(i, k-1)*delp_tl(i, k-1))/delp(i, k-1)**2)*(delp(i, k)+&
&           0.5*(delp(i, k-1)+delp(i, k+1)))-2.*(dc(i, k+1)/delp(i, k+1)&
&           -dc(i, k-1)/delp(i, k-1))*(delp_tl(i, k)+0.5*(delp_tl(i, k-1&
&           )+delp_tl(i, k+1))))*delp(i, k)**2/(delp(i, k)+0.5*(delp(i, &
&           k-1)+delp(i, k+1)))**2 + 2.*(dc(i, k+1)/delp(i, k+1)-dc(i, k&
&           -1)/delp(i, k-1))*2*delp(i, k)*delp_tl(i, k)/(delp(i, k)+0.5&
&           *(delp(i, k-1)+delp(i, k+1)))
          h2(i, k) = 2.*(dc(i, k+1)/delp(i, k+1)-dc(i, k-1)/delp(i, k-1)&
&           )/(delp(i, k)+0.5*(delp(i, k-1)+delp(i, k+1)))*delp(i, k)**2
        END DO
      END DO
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
! original quasi-monotone
      fac = 1.5
      DO k=3,km-2
        DO i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
          pmp_tl = 2.*dc_tl(i, k)
          pmp = 2.*dc(i, k)
          qmp_tl = a4_tl(1, i, k) + pmp_tl
          qmp = a4(1, i, k) + pmp
          lac_tl = a4_tl(1, i, k) + fac*h2_tl(i, k-1) + dc_tl(i, k)
          lac = a4(1, i, k) + fac*h2(i, k-1) + dc(i, k)
          IF (a4(1, i, k) .GT. qmp) THEN
            IF (qmp .GT. lac) THEN
              y8_tl = lac_tl
              y8 = lac
            ELSE
              y8_tl = qmp_tl
              y8 = qmp
            END IF
          ELSE IF (a4(1, i, k) .GT. lac) THEN
            y8_tl = lac_tl
            y8 = lac
          ELSE
            y8_tl = a4_tl(1, i, k)
            y8 = a4(1, i, k)
          END IF
          IF (a4(3, i, k) .LT. y8) THEN
            x2_tl = y8_tl
            x2 = y8
          ELSE
            x2_tl = a4_tl(3, i, k)
            x2 = a4(3, i, k)
          END IF
          IF (a4(1, i, k) .LT. qmp) THEN
            IF (qmp .LT. lac) THEN
              y6_tl = lac_tl
              y6 = lac
            ELSE
              y6_tl = qmp_tl
              y6 = qmp
            END IF
          ELSE IF (a4(1, i, k) .LT. lac) THEN
            y6_tl = lac_tl
            y6 = lac
          ELSE
            y6_tl = a4_tl(1, i, k)
            y6 = a4(1, i, k)
          END IF
          IF (x2 .GT. y6) THEN
            a4_tl(3, i, k) = y6_tl
            a4(3, i, k) = y6
          ELSE
            a4_tl(3, i, k) = x2_tl
            a4(3, i, k) = x2
          END IF
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
          qmp_tl = a4_tl(1, i, k) - pmp_tl
          qmp = a4(1, i, k) - pmp
          lac_tl = a4_tl(1, i, k) + fac*h2_tl(i, k+1) - dc_tl(i, k)
          lac = a4(1, i, k) + fac*h2(i, k+1) - dc(i, k)
          IF (a4(1, i, k) .GT. qmp) THEN
            IF (qmp .GT. lac) THEN
              y9_tl = lac_tl
              y9 = lac
            ELSE
              y9_tl = qmp_tl
              y9 = qmp
            END IF
          ELSE IF (a4(1, i, k) .GT. lac) THEN
            y9_tl = lac_tl
            y9 = lac
          ELSE
            y9_tl = a4_tl(1, i, k)
            y9 = a4(1, i, k)
          END IF
          IF (a4(2, i, k) .LT. y9) THEN
            x3_tl = y9_tl
            x3 = y9
          ELSE
            x3_tl = a4_tl(2, i, k)
            x3 = a4(2, i, k)
          END IF
          IF (a4(1, i, k) .LT. qmp) THEN
            IF (qmp .LT. lac) THEN
              y7_tl = lac_tl
              y7 = lac
            ELSE
              y7_tl = qmp_tl
              y7 = qmp
            END IF
          ELSE IF (a4(1, i, k) .LT. lac) THEN
            y7_tl = lac_tl
            y7 = lac
          ELSE
            y7_tl = a4_tl(1, i, k)
            y7 = a4(1, i, k)
          END IF
          IF (x3 .GT. y7) THEN
            a4_tl(2, i, k) = y7_tl
            a4(2, i, k) = y7
          ELSE
            a4_tl(2, i, k) = x3_tl
            a4(2, i, k) = x3
          END IF
!-------------
! Recompute A6
!-------------
          a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3&
&           , i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
! Additional constraint to ensure positivity when kord=7
        IF (iv .EQ. 0 .AND. kord .GE. 6) CALL PPM_LIMITERS_TLM(dc(i1, k)&
&                                                        , dc_tl(i1, k)&
&                                                        , a4(1, i1, k)&
&                                                        , a4_tl(1, i1, &
&                                                        k), it, 2)
      END DO
    ELSE
      lmt = kord - 3
      IF (0 .LT. lmt) THEN
        lmt = lmt
      ELSE
        lmt = 0
      END IF
      IF (iv .EQ. 0) THEN
        IF (2 .GT. lmt) THEN
          lmt = lmt
        ELSE
          lmt = 2
        END IF
      END IF
      DO k=3,km-2
        IF (kord .NE. 4) THEN
          DO i=i1,i2
            a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(&
&             3, i, k))
            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
          END DO
        END IF
        IF (kord .NE. 6) CALL PPM_LIMITERS_TLM(dc(i1, k), dc_tl(i1, k), &
&                                        a4(1, i1, k), a4_tl(1, i1, k), &
&                                        it, lmt)
      END DO
    END IF
    DO k=km1,km
      DO i=i1,i2
        a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3, i&
&         , k))
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
      CALL PPM_LIMITERS_TLM(dc(i1, k), dc_tl(i1, k), a4(1, i1, k), a4_tl&
&                     (1, i1, k), it, 0)
    END DO
  END SUBROUTINE PPM_PROFILE_TLM
!  Differentiation of ppm_limiters in forward (tangent) mode:
!   variations   of useful results: a4
!   with respect to varying inputs: dm a4
  SUBROUTINE PPM_LIMITERS_TLM(dm, dm_tl, a4, a4_tl, itot, lmt)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! the linear slope
    REAL, INTENT(IN) :: dm(*)
    REAL, INTENT(IN) :: dm_tl(*)
! Total Longitudes
    INTEGER, INTENT(IN) :: itot
! 0: Standard PPM constraint
    INTEGER, INTENT(IN) :: lmt
! 1: Improved full monotonicity constraint (Lin)
! 2: Positive definite constraint
! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
! PPM array
    REAL, INTENT(INOUT) :: a4(4, *)
    REAL, INTENT(INOUT) :: a4_tl(4, *)
! AA <-- a4(1,i)
! AL <-- a4(2,i)
! AR <-- a4(3,i)
! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
    REAL :: qmp
    REAL :: qmp_tl
    REAL :: da1, da2, a6da
    REAL :: fmin
    INTEGER :: i
    INTRINSIC ABS
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min1_tl
    REAL :: min2
    REAL :: min2_tl
    REAL :: abs0
    REAL :: x1_tl
    REAL :: y1_tl
    REAL :: x2_tl
    REAL :: y2_tl
    REAL :: x2
    REAL :: x1
    REAL :: y2
    REAL :: y1
! Developer: S.-J. Lin
    IF (lmt .EQ. 3) THEN
      RETURN
    ELSE IF (lmt .EQ. 0) THEN
! Standard PPM constraint
      DO i=1,itot
        IF (dm(i) .EQ. 0.) THEN
          a4_tl(2, i) = a4_tl(1, i)
          a4(2, i) = a4(1, i)
          a4_tl(3, i) = a4_tl(1, i)
          a4(3, i) = a4(1, i)
          a4_tl(4, i) = 0.0
          a4(4, i) = 0.
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
            a4(3, i) = a4(2, i) - a4(4, i)
          ELSE IF (a6da .GT. da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
            a4(2, i) = a4(3, i) - a4(4, i)
          END IF
        END IF
      END DO
    ELSE IF (lmt .EQ. 1) THEN
! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      DO i=1,itot
        qmp_tl = 2.*dm_tl(i)
        qmp = 2.*dm(i)
        IF (qmp .GE. 0.) THEN
          x1_tl = qmp_tl
          x1 = qmp
        ELSE
          x1_tl = -qmp_tl
          x1 = -qmp
        END IF
        IF (a4(2, i) - a4(1, i) .GE. 0.) THEN
          y1_tl = a4_tl(2, i) - a4_tl(1, i)
          y1 = a4(2, i) - a4(1, i)
        ELSE
          y1_tl = -(a4_tl(2, i)-a4_tl(1, i))
          y1 = -(a4(2, i)-a4(1, i))
        END IF
        IF (x1 .GT. y1) THEN
          min1_tl = y1_tl
          min1 = y1
        ELSE
          min1_tl = x1_tl
          min1 = x1
        END IF
        a4_tl(2, i) = a4_tl(1, i) - min1_tl*SIGN(1.d0, min1*qmp)
        a4(2, i) = a4(1, i) - SIGN(min1, qmp)
        IF (qmp .GE. 0.) THEN
          x2_tl = qmp_tl
          x2 = qmp
        ELSE
          x2_tl = -qmp_tl
          x2 = -qmp
        END IF
        IF (a4(3, i) - a4(1, i) .GE. 0.) THEN
          y2_tl = a4_tl(3, i) - a4_tl(1, i)
          y2 = a4(3, i) - a4(1, i)
        ELSE
          y2_tl = -(a4_tl(3, i)-a4_tl(1, i))
          y2 = -(a4(3, i)-a4(1, i))
        END IF
        IF (x2 .GT. y2) THEN
          min2_tl = y2_tl
          min2 = y2
        ELSE
          min2_tl = x2_tl
          min2 = x2
        END IF
        a4_tl(3, i) = a4_tl(1, i) + min2_tl*SIGN(1.d0, min2*qmp)
        a4(3, i) = a4(1, i) + SIGN(min2, qmp)
        a4_tl(4, i) = 3.*(2.*a4_tl(1, i)-a4_tl(2, i)-a4_tl(3, i))
        a4(4, i) = 3.*(2.*a4(1, i)-(a4(2, i)+a4(3, i)))
      END DO
    ELSE IF (lmt .EQ. 2) THEN
! Positive definite constraint
      DO i=1,itot
        IF (a4(3, i) - a4(2, i) .GE. 0.) THEN
          abs0 = a4(3, i) - a4(2, i)
        ELSE
          abs0 = -(a4(3, i)-a4(2, i))
        END IF
        IF (abs0 .LT. -a4(4, i)) THEN
          fmin = a4(1, i) + 0.25*(a4(3, i)-a4(2, i))**2/a4(4, i) + a4(4&
&           , i)*r12
          IF (fmin .LT. 0.) THEN
            IF (a4(1, i) .LT. a4(3, i) .AND. a4(1, i) .LT. a4(2, i)) &
&           THEN
              a4_tl(3, i) = a4_tl(1, i)
              a4(3, i) = a4(1, i)
              a4_tl(2, i) = a4_tl(1, i)
              a4(2, i) = a4(1, i)
              a4_tl(4, i) = 0.0
              a4(4, i) = 0.
            ELSE IF (a4(3, i) .GT. a4(2, i)) THEN
              a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
              a4(4, i) = 3.*(a4(2, i)-a4(1, i))
              a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
              a4(3, i) = a4(2, i) - a4(4, i)
            ELSE
              a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
              a4(4, i) = 3.*(a4(3, i)-a4(1, i))
              a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
              a4(2, i) = a4(3, i) - a4(4, i)
            END IF
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE PPM_LIMITERS_TLM
  SUBROUTINE STEEPZ(i1, i2, km, a4, df2, dm, dq, dp, d4)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km, i1, i2
! grid size
    REAL, INTENT(IN) :: dp(i1:i2, km)
! backward diff of q
    REAL, INTENT(IN) :: dq(i1:i2, km)
! backward sum:  dp(k)+ dp(k-1)
    REAL, INTENT(IN) :: d4(i1:i2, km)
! first guess mismatch
    REAL, INTENT(IN) :: df2(i1:i2, km)
! monotonic mismatch
    REAL, INTENT(IN) :: dm(i1:i2, km)
! !INPUT/OUTPUT PARAMETERS:
! first guess/steepened
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
! !LOCAL VARIABLES:
    INTEGER :: i, k
    REAL :: alfa(i1:i2, km)
    REAL :: f(i1:i2, km)
    REAL :: rat(i1:i2, km)
    REAL :: dg2
    INTRINSIC MIN
    INTRINSIC MAX
    REAL :: y1
! Compute ratio of dq/dp
    DO k=2,km
      DO i=i1,i2
        rat(i, k) = dq(i, k-1)/d4(i, k)
      END DO
    END DO
! Compute F
    DO k=2,km-1
      DO i=i1,i2
        f(i, k) = (rat(i, k+1)-rat(i, k))/(dp(i, k-1)+dp(i, k)+dp(i, k+1&
&         ))
      END DO
    END DO
    DO k=3,km-2
      DO i=i1,i2
        IF (f(i, k+1)*f(i, k-1) .LT. 0. .AND. df2(i, k) .NE. 0.) THEN
          dg2 = (f(i, k+1)-f(i, k-1))*((dp(i, k+1)-dp(i, k-1))**2+d4(i, &
&           k)*d4(i, k+1))
          IF (0.5 .GT. -(0.1875*dg2/df2(i, k))) THEN
            y1 = -(0.1875*dg2/df2(i, k))
          ELSE
            y1 = 0.5
          END IF
          IF (0. .LT. y1) THEN
            alfa(i, k) = y1
          ELSE
            alfa(i, k) = 0.
          END IF
        ELSE
          alfa(i, k) = 0.
        END IF
      END DO
    END DO
    DO k=4,km-2
      DO i=i1,i2
        a4(2, i, k) = (1.-alfa(i, k-1)-alfa(i, k))*a4(2, i, k) + alfa(i&
&         , k-1)*(a4(1, i, k)-dm(i, k)) + alfa(i, k)*(a4(1, i, k-1)+dm(i&
&         , k-1))
      END DO
    END DO
  END SUBROUTINE STEEPZ
  SUBROUTINE RST_REMAP(km, kn, is, ie, js, je, isd, ied, jsd, jed, nq, &
&   ntp, delp_r, u_r, v_r, w_r, delz_r, pt_r, q_r, qdiag_r, delp, u, v, &
&   w, delz, pt, q, qdiag, ak_r, bk_r, ptop, ak, bk, hydrostatic, &
&   make_nh, domain, square_domain)
    IMPLICIT NONE
!------------------------------------
! Assuming hybrid sigma-P coordinate:
!------------------------------------
! !INPUT PARAMETERS:
! Restart z-dimension
    INTEGER, INTENT(IN) :: km
! Run time dimension
    INTEGER, INTENT(IN) :: kn
! number of tracers (including h2o)
    INTEGER, INTENT(IN) :: nq, ntp
! starting & ending X-Dir index
    INTEGER, INTENT(IN) :: is, ie, isd, ied
! starting & ending Y-Dir index
    INTEGER, INTENT(IN) :: js, je, jsd, jed
    LOGICAL, INTENT(IN) :: hydrostatic, make_nh, square_domain
    REAL, INTENT(IN) :: ptop
    REAL, INTENT(IN) :: ak_r(km+1)
    REAL, INTENT(IN) :: bk_r(km+1)
    REAL, INTENT(IN) :: ak(kn+1)
    REAL, INTENT(IN) :: bk(kn+1)
! pressure thickness
    REAL, INTENT(IN) :: delp_r(is:ie, js:je, km)
! u-wind (m/s)
    REAL, INTENT(IN) :: u_r(is:ie, js:je+1, km)
! v-wind (m/s)
    REAL, INTENT(IN) :: v_r(is:ie+1, js:je, km)
    REAL, INTENT(INOUT) :: pt_r(is:ie, js:je, km)
    REAL, INTENT(IN) :: w_r(is:ie, js:je, km)
    REAL, INTENT(IN) :: q_r(is:ie, js:je, km, ntp)
    REAL, INTENT(IN) :: qdiag_r(is:ie, js:je, km, ntp+1:nq)
    REAL, INTENT(INOUT) :: delz_r(is:ie, js:je, km)
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Output:
! pressure thickness
    REAL, INTENT(OUT) :: delp(isd:ied, jsd:jed, kn)
! u-wind (m/s)
    REAL, INTENT(OUT) :: u(isd:ied, jsd:jed+1, kn)
! v-wind (m/s)
    REAL, INTENT(OUT) :: v(isd:ied+1, jsd:jed, kn)
! vertical velocity (m/s)
    REAL, INTENT(OUT) :: w(isd:, jsd:, :)
! temperature
    REAL, INTENT(OUT) :: pt(isd:ied, jsd:jed, kn)
    REAL, INTENT(OUT) :: q(isd:ied, jsd:jed, kn, ntp)
    REAL, INTENT(OUT) :: qdiag(isd:ied, jsd:jed, kn, ntp+1:nq)
! delta-height (m)
    REAL, INTENT(OUT) :: delz(isd:, jsd:, :)
!-----------------------------------------------------------------------
    REAL :: r_vir, rgrav
! surface pressure
    REAL :: ps(isd:ied, jsd:jed)
    REAL :: pe1(is:ie, km+1)
    REAL :: pe2(is:ie, kn+1)
    REAL :: pv1(is:ie+1, km+1)
    REAL :: pv2(is:ie+1, kn+1)
    INTEGER :: i, j, k, iq
    INTEGER, PARAMETER :: kord=4
    INTRINSIC LOG
    r_vir = rvgas/rdgas - 1.
    rgrav = 1./grav
!$OMP parallel do default(none) shared(is,ie,js,je,ps,ak_r)
    DO j=js,je
      DO i=is,ie
        ps(i, j) = ak_r(1)
      END DO
    END DO
! this OpenMP do-loop setup cannot work in it's current form....
!$OMP parallel do default(none) shared(is,ie,js,je,km,ps,delp_r)
    DO j=js,je
      DO k=1,km
        DO i=is,ie
          ps(i, j) = ps(i, j) + delp_r(i, j, k)
        END DO
      END DO
    END DO
! only one cell is needed
    IF (square_domain) THEN
      CALL MPP_UPDATE_DOMAINS(ps, domain, complete=.true., whalo=1, &
&                       ehalo=1, shalo=1, nhalo=1)
    ELSE
      CALL MPP_UPDATE_DOMAINS(ps, domain, complete=.true.)
    END IF
! Compute virtual Temp
!$OMP parallel do default(none) shared(is,ie,js,je,km,pt_r,r_vir,q_r)
    DO k=1,km
      DO j=js,je
        DO i=is,ie
          pt_r(i, j, k) = pt_r(i, j, k)*(1.+r_vir*q_r(i, j, k, 1))
        END DO
      END DO
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,km,ak_r,bk_r,ps,kn,ak,bk,u_r,u,delp, &
!$OMP                                  ntp,nq,hydrostatic,make_nh,w_r,w,delz_r,delp_r,delz, &
!$OMP                                  pt_r,pt,v_r,v,q,q_r,qdiag,qdiag_r) &
!$OMP                          private(pe1,  pe2, pv1, pv2)
    DO j=js,je+1
!------
! map u
!------
      DO k=1,km+1
        DO i=is,ie
          pe1(i, k) = ak_r(k) + 0.5*bk_r(k)*(ps(i, j-1)+ps(i, j))
        END DO
      END DO
      DO k=1,kn+1
        DO i=is,ie
          pe2(i, k) = ak(k) + 0.5*bk(k)*(ps(i, j-1)+ps(i, j))
        END DO
      END DO
      CALL REMAP_2D(km, pe1, u_r(is:ie, j:j, 1:km), kn, pe2, u(is:ie, j:&
&             j, 1:kn), is, ie, -1, kord)
!(j < je+1)
      IF (j .NE. je + 1) THEN
!---------------
! Hybrid sigma-p
!---------------
        DO k=1,km+1
          DO i=is,ie
            pe1(i, k) = ak_r(k) + bk_r(k)*ps(i, j)
          END DO
        END DO
        DO k=1,kn+1
          DO i=is,ie
            pe2(i, k) = ak(k) + bk(k)*ps(i, j)
          END DO
        END DO
!-------------
! Compute delp
!-------------
        DO k=1,kn
          DO i=is,ie
            delp(i, j, k) = pe2(i, k+1) - pe2(i, k)
          END DO
        END DO
!----------------
! Map constituents
!----------------
        IF (nq .NE. 0) THEN
          DO iq=1,ntp
            CALL REMAP_2D(km, pe1, q_r(is:ie, j:j, 1:km, iq:iq), kn, pe2&
&                   , q(is:ie, j:j, 1:kn, iq:iq), is, ie, 0, kord)
          END DO
          DO iq=ntp+1,nq
            CALL REMAP_2D(km, pe1, qdiag_r(is:ie, j:j, 1:km, iq:iq), kn&
&                   , pe2, qdiag(is:ie, j:j, 1:kn, iq:iq), is, ie, 0, &
&                   kord)
          END DO
        END IF
        IF (.NOT.hydrostatic .AND. (.NOT.make_nh)) THEN
! Remap vertical wind:
          CALL REMAP_2D(km, pe1, w_r(is:ie, j:j, 1:km), kn, pe2, w(is:ie&
&                 , j:j, 1:kn), is, ie, -1, kord)
! Remap delz for hybrid sigma-p coordinate
          DO k=1,km
            DO i=is,ie
! ="specific volume"/grav
              delz_r(i, j, k) = -(delz_r(i, j, k)/delp_r(i, j, k))
            END DO
          END DO
          CALL REMAP_2D(km, pe1, delz_r(is:ie, j:j, 1:km), kn, pe2, delz&
&                 (is:ie, j:j, 1:kn), is, ie, 1, kord)
          DO k=1,kn
            DO i=is,ie
              delz(i, j, k) = -(delz(i, j, k)*delp(i, j, k))
            END DO
          END DO
        END IF
! Geopotential conserving remap of virtual temperature:
        DO k=1,km+1
          DO i=is,ie
            pe1(i, k) = LOG(pe1(i, k))
          END DO
        END DO
        DO k=1,kn+1
          DO i=is,ie
            pe2(i, k) = LOG(pe2(i, k))
          END DO
        END DO
        CALL REMAP_2D(km, pe1, pt_r(is:ie, j:j, 1:km), kn, pe2, pt(is:ie&
&               , j:j, 1:kn), is, ie, 1, kord)
!------
! map v
!------
        DO k=1,km+1
          DO i=is,ie+1
            pv1(i, k) = ak_r(k) + 0.5*bk_r(k)*(ps(i-1, j)+ps(i, j))
          END DO
        END DO
        DO k=1,kn+1
          DO i=is,ie+1
            pv2(i, k) = ak(k) + 0.5*bk(k)*(ps(i-1, j)+ps(i, j))
          END DO
        END DO
        CALL REMAP_2D(km, pv1, v_r(is:ie+1, j:j, 1:km), kn, pv2, v(is:ie&
&               +1, j:j, 1:kn), is, ie + 1, -1, kord)
      END IF
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,kn,pt,r_vir,q)
    DO k=1,kn
      DO j=js,je
        DO i=is,ie
          pt(i, j, k) = pt(i, j, k)/(1.+r_vir*q(i, j, k, 1))
        END DO
      END DO
    END DO
  END SUBROUTINE RST_REMAP
  SUBROUTINE REMAP_2D(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1, i2
! Mode: 0 ==  constituents 1 ==others
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! Field input
    REAL, INTENT(IN) :: q1(i1:i2, km)
! Field output
    REAL, INTENT(OUT) :: q2(i1:i2, kn)
! !LOCAL VARIABLES:
    REAL :: qs(i1:i2)
    REAL :: dp1(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    INTEGER :: i, k, l, m, k0
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q1(i, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE(qs, q4, dp1, km, i1, i2, iv, kord)
    ELSE
      CALL PPM_PROFILE(q4, dp1, km, i1, i2, iv, kord)
    END IF
    DO i=i1,i2
      k0 = 1
      DO k=1,kn
        IF (pe2(i, k) .LE. pe1(i, 1)) THEN
! above old ptop:
          q2(i, k) = q1(i, 1)
        ELSE
          DO l=k0,km
! locate the top edge: pe2(i,k)
            IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1&
&               )) THEN
              pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
              IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
                pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
                q2(i, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4&
&                 (2, i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
                k0 = l
                GOTO 555
              ELSE
! Fractional area...
                qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i&
&                 , l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*&
&                 (1.+pl*(1.+pl))))
                DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                  IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
                    qsum = qsum + dp1(i, m)*q4(1, i, m)
                  ELSE
                    dp = pe2(i, k+1) - pe1(i, m)
                    esl = dp/dp1(i, m)
                    qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-&
&                     q4(2, i, m)+q4(4, i, m)*(1.-r23*esl)))
                    k0 = m
                    GOTO 123
                  END IF
                END DO
                GOTO 123
              END IF
            END IF
          END DO
 123      q2(i, k) = qsum/(pe2(i, k+1)-pe2(i, k))
        END IF
 555    CONTINUE
      END DO
    END DO
  END SUBROUTINE REMAP_2D
  SUBROUTINE MAPPM(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord, ptop)
    IMPLICIT NONE
! IV = 0: constituents
! IV = 1: potential temp
! IV =-1: winds
! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
    INTEGER, INTENT(IN) :: i1, i2, km, kn, kord, iv
    REAL, INTENT(IN) :: pe1(i1:i2, km+1), pe2(i1:i2, kn+1)
    REAL, INTENT(IN) :: q1(i1:i2, km)
    REAL, INTENT(OUT) :: q2(i1:i2, kn)
    REAL, INTENT(IN) :: ptop
! local
    REAL :: qs(i1:i2)
    REAL :: dp1(i1:i2, km)
    REAL :: a4(4, i1:i2, km)
    INTEGER :: i, k, l
    INTEGER :: k0, k1
    REAL :: pl, pr, tt, delp, qsum, dpsum, esl
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        a4(1, i, k) = q1(i, k)
      END DO
    END DO
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE(qs, a4, dp1, km, i1, i2, iv, kord)
    ELSE
      CALL PPM_PROFILE(a4, dp1, km, i1, i2, iv, kord)
    END IF
!------------------------------------
! Lowest layer: constant distribution
!------------------------------------
    DO i=i1,i2
      k0 = 1
      DO k=1,kn
        IF (pe2(i, k) .LE. pe1(i, 1)) THEN
! above old ptop
          q2(i, k) = q1(i, 1)
        ELSE IF (pe2(i, k) .GE. pe1(i, km+1)) THEN
! Entire grid below old ps
          q2(i, k) = q1(i, km)
        ELSE
          DO l=k0,km
! locate the top edge at pe2(i,k)
            IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1&
&               )) THEN
              k0 = l
              pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
              IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
                pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
                tt = r3*(pr*(pr+pl)+pl**2)
                q2(i, k) = a4(2, i, l) + 0.5*(a4(4, i, l)+a4(3, i, l)-a4&
&                 (2, i, l))*(pr+pl) - a4(4, i, l)*tt
                GOTO 555
              ELSE
! Fractional area...
                delp = pe1(i, l+1) - pe2(i, k)
                tt = r3*(1.+pl*(1.+pl))
                qsum = delp*(a4(2, i, l)+0.5*(a4(4, i, l)+a4(3, i, l)-a4&
&                 (2, i, l))*(1.+pl)-a4(4, i, l)*tt)
                dpsum = delp
                k1 = l + 1
                GOTO 111
              END IF
            END IF
          END DO
 111      CONTINUE
          DO l=k1,km
            IF (pe2(i, k+1) .GT. pe1(i, l+1)) THEN
! Whole layer..
              qsum = qsum + dp1(i, l)*q1(i, l)
              dpsum = dpsum + dp1(i, l)
            ELSE
              delp = pe2(i, k+1) - pe1(i, l)
              esl = delp/dp1(i, l)
              qsum = qsum + delp*(a4(2, i, l)+0.5*esl*(a4(3, i, l)-a4(2&
&               , i, l)+a4(4, i, l)*(1.-r23*esl)))
              dpsum = dpsum + delp
              k0 = l
              GOTO 123
            END IF
          END DO
          delp = pe2(i, k+1) - pe1(i, km+1)
          IF (delp .GT. 0.) THEN
! Extended below old ps
            qsum = qsum + delp*q1(i, km)
            dpsum = dpsum + delp
          END IF
 123      q2(i, k) = qsum/dpsum
        END IF
 555    CONTINUE
      END DO
    END DO
  END SUBROUTINE MAPPM
  SUBROUTINE CS_PROFILE(qs, a4, delp, km, i1, i2, iv, kord)
    IMPLICIT NONE
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
    INTEGER, INTENT(IN) :: i1, i2
! vertical dimension
    INTEGER, INTENT(IN) :: km
! iv =-1: winds
    INTEGER, INTENT(IN) :: iv
! iv = 0: positive definite scalars
! iv = 1: others
    INTEGER, INTENT(IN) :: kord
    REAL, INTENT(IN) :: qs(i1:i2)
! layer pressure thickness
    REAL, INTENT(IN) :: delp(i1:i2, km)
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
!-----------------------------------------------------------------------
    LOGICAL :: extm(i1:i2, km)
    REAL :: gam(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    INTEGER :: i, k, im
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: abs0
    INTEGER :: abs1
    INTEGER :: abs2
    REAL :: abs3
    INTEGER :: abs4
    REAL :: abs5
    INTEGER :: abs6
    REAL :: abs7
    INTEGER :: abs8
    INTEGER :: abs9
    REAL :: abs10
    REAL :: abs11
    REAL :: abs12
    REAL :: x10
    REAL :: y28
    REAL :: y27
    REAL :: y26
    REAL :: y25
    REAL :: y24
    REAL :: y23
    REAL :: y22
    REAL :: y21
    REAL :: y20
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: x1
    REAL :: y19
    REAL :: y18
    REAL :: y17
    REAL :: y16
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
    REAL :: y9
    REAL :: y8
    REAL :: y7
    REAL :: y6
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: y2
    REAL :: y1
    IF (iv .EQ. -2) THEN
      DO i=i1,i2
        gam(i, 2) = 0.5
        q(i, 1) = 1.5*a4(1, i, 1)
      END DO
      DO k=2,km-1
        DO i=i1,i2
          grat = delp(i, k-1)/delp(i, k)
          bet = 2. + grat + grat - gam(i, k)
          q(i, k) = (3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))/bet
          gam(i, k+1) = grat/bet
        END DO
      END DO
      DO i=i1,i2
        grat = delp(i, km-1)/delp(i, km)
        q(i, km) = (3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-&
&         1))/(2.+grat+grat-gam(i, km))
        q(i, km+1) = qs(i)
      END DO
      DO k=km-1,1,-1
        DO i=i1,i2
          q(i, k) = q(i, k) - gam(i, k+1)*q(i, k+1)
        END DO
      END DO
    ELSE
      DO i=i1,i2
! grid ratio
        grat = delp(i, 2)/delp(i, 1)
        bet = grat*(grat+0.5)
        q(i, 1) = ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))/bet
        gam(i, 1) = (1.+grat*(grat+1.5))/bet
      END DO
      DO k=2,km
        DO i=i1,i2
          d4(i) = delp(i, k-1)/delp(i, k)
          bet = 2. + d4(i) + d4(i) - gam(i, k-1)
          q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
          gam(i, k) = d4(i)/bet
        END DO
      END DO
      DO i=i1,i2
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
        END DO
      END DO
    END IF
    IF (kord .GE. 0.) THEN
      abs0 = kord
    ELSE
      abs0 = -kord
    END IF
!----- Perfectly linear scheme --------------------------------
    IF (abs0 .GT. 16) THEN
      DO k=1,km
        DO i=i1,i2
          a4(2, i, k) = q(i, k)
          a4(3, i, k) = q(i, k+1)
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
      END DO
      RETURN
    ELSE
!----- Perfectly linear scheme --------------------------------
!------------------
! Apply constraints
!------------------
      im = i2 - i1 + 1
! Apply *large-scale* constraints
      DO i=i1,i2
        IF (a4(1, i, 1) .LT. a4(1, i, 2)) THEN
          y1 = a4(1, i, 2)
        ELSE
          y1 = a4(1, i, 1)
        END IF
        IF (q(i, 2) .GT. y1) THEN
          q(i, 2) = y1
        ELSE
          q(i, 2) = q(i, 2)
        END IF
        IF (a4(1, i, 1) .GT. a4(1, i, 2)) THEN
          y2 = a4(1, i, 2)
        ELSE
          y2 = a4(1, i, 1)
        END IF
        IF (q(i, 2) .LT. y2) THEN
          q(i, 2) = y2
        ELSE
          q(i, 2) = q(i, 2)
        END IF
      END DO
      DO k=2,km
        DO i=i1,i2
          gam(i, k) = a4(1, i, k) - a4(1, i, k-1)
        END DO
      END DO
! Interior:
      DO k=3,km-1
        DO i=i1,i2
          IF (gam(i, k-1)*gam(i, k+1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y3 = a4(1, i, k)
            ELSE
              y3 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .GT. y3) THEN
              q(i, k) = y3
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y4 = a4(1, i, k)
            ELSE
              y4 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .LT. y4) THEN
              q(i, k) = y4
            ELSE
              q(i, k) = q(i, k)
            END IF
          ELSE IF (gam(i, k-1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y5 = a4(1, i, k)
            ELSE
              y5 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .LT. y5) THEN
              q(i, k) = y5
            ELSE
              q(i, k) = q(i, k)
            END IF
          ELSE
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y6 = a4(1, i, k)
            ELSE
              y6 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .GT. y6) THEN
              q(i, k) = y6
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (iv .EQ. 0) THEN
              IF (0. .LT. q(i, k)) THEN
                q(i, k) = q(i, k)
              ELSE
                q(i, k) = 0.
              END IF
            END IF
          END IF
        END DO
      END DO
! Bottom:
      DO i=i1,i2
        IF (a4(1, i, km-1) .LT. a4(1, i, km)) THEN
          y7 = a4(1, i, km)
        ELSE
          y7 = a4(1, i, km-1)
        END IF
        IF (q(i, km) .GT. y7) THEN
          q(i, km) = y7
        ELSE
          q(i, km) = q(i, km)
        END IF
        IF (a4(1, i, km-1) .GT. a4(1, i, km)) THEN
          y8 = a4(1, i, km)
        ELSE
          y8 = a4(1, i, km-1)
        END IF
        IF (q(i, km) .LT. y8) THEN
          q(i, km) = y8
        ELSE
          q(i, km) = q(i, km)
        END IF
      END DO
      DO k=1,km
        DO i=i1,i2
          a4(2, i, k) = q(i, k)
          a4(3, i, k) = q(i, k+1)
        END DO
      END DO
      DO k=1,km
        IF (k .EQ. 1 .OR. k .EQ. km) THEN
          DO i=i1,i2
            extm(i, k) = (a4(2, i, k)-a4(1, i, k))*(a4(3, i, k)-a4(1, i&
&             , k)) .GT. 0.
          END DO
        ELSE
          DO i=i1,i2
            extm(i, k) = gam(i, k)*gam(i, k+1) .LT. 0.
          END DO
        END IF
      END DO
!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping
      IF (iv .EQ. 0) THEN
        DO i=i1,i2
          IF (0. .LT. a4(2, i, 1)) THEN
            a4(2, i, 1) = a4(2, i, 1)
          ELSE
            a4(2, i, 1) = 0.
          END IF
        END DO
      ELSE IF (iv .EQ. -1) THEN
        DO i=i1,i2
          IF (a4(2, i, 1)*a4(1, i, 1) .LE. 0.) a4(2, i, 1) = 0.
        END DO
      ELSE IF (iv .EQ. 2) THEN
        DO i=i1,i2
          a4(2, i, 1) = a4(1, i, 1)
          a4(3, i, 1) = a4(1, i, 1)
          a4(4, i, 1) = 0.
        END DO
      END IF
      IF (iv .NE. 2) THEN
        DO i=i1,i2
          a4(4, i, 1) = 3.*(2.*a4(1, i, 1)-(a4(2, i, 1)+a4(3, i, 1)))
        END DO
        CALL CS_LIMITERS(im, extm(i1, 1), a4(1, i1, 1), 1)
      END IF
! k=2
      DO i=i1,i2
        a4(4, i, 2) = 3.*(2.*a4(1, i, 2)-(a4(2, i, 2)+a4(3, i, 2)))
      END DO
      CALL CS_LIMITERS(im, extm(i1, 2), a4(1, i1, 2), 2)
!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
      DO k=3,km-2
        IF (kord .GE. 0.) THEN
          abs1 = kord
        ELSE
          abs1 = -kord
        END IF
        IF (abs1 .LT. 9) THEN
          DO i=i1,i2
! Left  edges
            pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
            lac_1 = pmp_1 + 1.5*gam(i, k+2)
            IF (a4(1, i, k) .GT. pmp_1) THEN
              IF (pmp_1 .GT. lac_1) THEN
                y19 = lac_1
              ELSE
                y19 = pmp_1
              END IF
            ELSE IF (a4(1, i, k) .GT. lac_1) THEN
              y19 = lac_1
            ELSE
              y19 = a4(1, i, k)
            END IF
            IF (a4(2, i, k) .LT. y19) THEN
              x1 = y19
            ELSE
              x1 = a4(2, i, k)
            END IF
            IF (a4(1, i, k) .LT. pmp_1) THEN
              IF (pmp_1 .LT. lac_1) THEN
                y9 = lac_1
              ELSE
                y9 = pmp_1
              END IF
            ELSE IF (a4(1, i, k) .LT. lac_1) THEN
              y9 = lac_1
            ELSE
              y9 = a4(1, i, k)
            END IF
            IF (x1 .GT. y9) THEN
              a4(2, i, k) = y9
            ELSE
              a4(2, i, k) = x1
            END IF
! Right edges
            pmp_2 = a4(1, i, k) + 2.*gam(i, k)
            lac_2 = pmp_2 - 1.5*gam(i, k-1)
            IF (a4(1, i, k) .GT. pmp_2) THEN
              IF (pmp_2 .GT. lac_2) THEN
                y20 = lac_2
              ELSE
                y20 = pmp_2
              END IF
            ELSE IF (a4(1, i, k) .GT. lac_2) THEN
              y20 = lac_2
            ELSE
              y20 = a4(1, i, k)
            END IF
            IF (a4(3, i, k) .LT. y20) THEN
              x2 = y20
            ELSE
              x2 = a4(3, i, k)
            END IF
            IF (a4(1, i, k) .LT. pmp_2) THEN
              IF (pmp_2 .LT. lac_2) THEN
                y10 = lac_2
              ELSE
                y10 = pmp_2
              END IF
            ELSE IF (a4(1, i, k) .LT. lac_2) THEN
              y10 = lac_2
            ELSE
              y10 = a4(1, i, k)
            END IF
            IF (x2 .GT. y10) THEN
              a4(3, i, k) = y10
            ELSE
              a4(3, i, k) = x2
            END IF
            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
          END DO
        ELSE
          IF (kord .GE. 0.) THEN
            abs2 = kord
          ELSE
            abs2 = -kord
          END IF
          IF (abs2 .EQ. 9) THEN
            DO i=i1,i2
              IF (extm(i, k) .AND. extm(i, k-1)) THEN
! c90_mp122
! grid-scale 2-delta-z wave detected
                a4(2, i, k) = a4(1, i, k)
                a4(3, i, k) = a4(1, i, k)
                a4(4, i, k) = 0.
              ELSE IF (extm(i, k) .AND. extm(i, k+1)) THEN
! c90_mp122
! grid-scale 2-delta-z wave detected
                a4(2, i, k) = a4(1, i, k)
                a4(3, i, k) = a4(1, i, k)
                a4(4, i, k) = 0.
              ELSE
                a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i, &
&                 k))
                IF (a4(4, i, k) .GE. 0.) THEN
                  abs3 = a4(4, i, k)
                ELSE
                  abs3 = -a4(4, i, k)
                END IF
                IF (a4(2, i, k) - a4(3, i, k) .GE. 0.) THEN
                  abs10 = a4(2, i, k) - a4(3, i, k)
                ELSE
                  abs10 = -(a4(2, i, k)-a4(3, i, k))
                END IF
! Check within the smooth region if subgrid profile is non-monotonic
                IF (abs3 .GT. abs10) THEN
                  pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                  lac_1 = pmp_1 + 1.5*gam(i, k+2)
                  IF (a4(1, i, k) .GT. pmp_1) THEN
                    IF (pmp_1 .GT. lac_1) THEN
                      y21 = lac_1
                    ELSE
                      y21 = pmp_1
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                    y21 = lac_1
                  ELSE
                    y21 = a4(1, i, k)
                  END IF
                  IF (a4(2, i, k) .LT. y21) THEN
                    x3 = y21
                  ELSE
                    x3 = a4(2, i, k)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_1) THEN
                    IF (pmp_1 .LT. lac_1) THEN
                      y11 = lac_1
                    ELSE
                      y11 = pmp_1
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                    y11 = lac_1
                  ELSE
                    y11 = a4(1, i, k)
                  END IF
                  IF (x3 .GT. y11) THEN
                    a4(2, i, k) = y11
                  ELSE
                    a4(2, i, k) = x3
                  END IF
                  pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                  lac_2 = pmp_2 - 1.5*gam(i, k-1)
                  IF (a4(1, i, k) .GT. pmp_2) THEN
                    IF (pmp_2 .GT. lac_2) THEN
                      y22 = lac_2
                    ELSE
                      y22 = pmp_2
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                    y22 = lac_2
                  ELSE
                    y22 = a4(1, i, k)
                  END IF
                  IF (a4(3, i, k) .LT. y22) THEN
                    x4 = y22
                  ELSE
                    x4 = a4(3, i, k)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_2) THEN
                    IF (pmp_2 .LT. lac_2) THEN
                      y12 = lac_2
                    ELSE
                      y12 = pmp_2
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                    y12 = lac_2
                  ELSE
                    y12 = a4(1, i, k)
                  END IF
                  IF (x4 .GT. y12) THEN
                    a4(3, i, k) = y12
                  ELSE
                    a4(3, i, k) = x4
                  END IF
                  a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i&
&                   , k))
                END IF
              END IF
            END DO
          ELSE
            IF (kord .GE. 0.) THEN
              abs4 = kord
            ELSE
              abs4 = -kord
            END IF
            IF (abs4 .EQ. 10) THEN
              DO i=i1,i2
                IF (extm(i, k)) THEN
                  IF (extm(i, k-1) .OR. extm(i, k+1)) THEN
! grid-scale 2-delta-z wave detected
                    a4(2, i, k) = a4(1, i, k)
                    a4(3, i, k) = a4(1, i, k)
                    a4(4, i, k) = 0.
                  ELSE
! True local extremum
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                  END IF
                ELSE
! not a local extremum
                  a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i&
&                   , k))
                  IF (a4(4, i, k) .GE. 0.) THEN
                    abs5 = a4(4, i, k)
                  ELSE
                    abs5 = -a4(4, i, k)
                  END IF
                  IF (a4(2, i, k) - a4(3, i, k) .GE. 0.) THEN
                    abs11 = a4(2, i, k) - a4(3, i, k)
                  ELSE
                    abs11 = -(a4(2, i, k)-a4(3, i, k))
                  END IF
! Check within the smooth region if subgrid profile is non-monotonic
                  IF (abs5 .GT. abs11) THEN
                    pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                    lac_1 = pmp_1 + 1.5*gam(i, k+2)
                    IF (a4(1, i, k) .GT. pmp_1) THEN
                      IF (pmp_1 .GT. lac_1) THEN
                        y23 = lac_1
                      ELSE
                        y23 = pmp_1
                      END IF
                    ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                      y23 = lac_1
                    ELSE
                      y23 = a4(1, i, k)
                    END IF
                    IF (a4(2, i, k) .LT. y23) THEN
                      x5 = y23
                    ELSE
                      x5 = a4(2, i, k)
                    END IF
                    IF (a4(1, i, k) .LT. pmp_1) THEN
                      IF (pmp_1 .LT. lac_1) THEN
                        y13 = lac_1
                      ELSE
                        y13 = pmp_1
                      END IF
                    ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                      y13 = lac_1
                    ELSE
                      y13 = a4(1, i, k)
                    END IF
                    IF (x5 .GT. y13) THEN
                      a4(2, i, k) = y13
                    ELSE
                      a4(2, i, k) = x5
                    END IF
                    pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                    lac_2 = pmp_2 - 1.5*gam(i, k-1)
                    IF (a4(1, i, k) .GT. pmp_2) THEN
                      IF (pmp_2 .GT. lac_2) THEN
                        y24 = lac_2
                      ELSE
                        y24 = pmp_2
                      END IF
                    ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                      y24 = lac_2
                    ELSE
                      y24 = a4(1, i, k)
                    END IF
                    IF (a4(3, i, k) .LT. y24) THEN
                      x6 = y24
                    ELSE
                      x6 = a4(3, i, k)
                    END IF
                    IF (a4(1, i, k) .LT. pmp_2) THEN
                      IF (pmp_2 .LT. lac_2) THEN
                        y14 = lac_2
                      ELSE
                        y14 = pmp_2
                      END IF
                    ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                      y14 = lac_2
                    ELSE
                      y14 = a4(1, i, k)
                    END IF
                    IF (x6 .GT. y14) THEN
                      a4(3, i, k) = y14
                    ELSE
                      a4(3, i, k) = x6
                    END IF
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                  END IF
                END IF
              END DO
            ELSE
              IF (kord .GE. 0.) THEN
                abs6 = kord
              ELSE
                abs6 = -kord
              END IF
              IF (abs6 .EQ. 12) THEN
                DO i=i1,i2
                  IF (extm(i, k)) THEN
! grid-scale 2-delta-z wave detected
                    a4(2, i, k) = a4(1, i, k)
                    a4(3, i, k) = a4(1, i, k)
                    a4(4, i, k) = 0.
                  ELSE
! not a local extremum
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                    IF (a4(4, i, k) .GE. 0.) THEN
                      abs7 = a4(4, i, k)
                    ELSE
                      abs7 = -a4(4, i, k)
                    END IF
                    IF (a4(2, i, k) - a4(3, i, k) .GE. 0.) THEN
                      abs12 = a4(2, i, k) - a4(3, i, k)
                    ELSE
                      abs12 = -(a4(2, i, k)-a4(3, i, k))
                    END IF
! Check within the smooth region if subgrid profile is non-monotonic
                    IF (abs7 .GT. abs12) THEN
                      pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                      lac_1 = pmp_1 + 1.5*gam(i, k+2)
                      IF (a4(1, i, k) .GT. pmp_1) THEN
                        IF (pmp_1 .GT. lac_1) THEN
                          y25 = lac_1
                        ELSE
                          y25 = pmp_1
                        END IF
                      ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                        y25 = lac_1
                      ELSE
                        y25 = a4(1, i, k)
                      END IF
                      IF (a4(2, i, k) .LT. y25) THEN
                        x7 = y25
                      ELSE
                        x7 = a4(2, i, k)
                      END IF
                      IF (a4(1, i, k) .LT. pmp_1) THEN
                        IF (pmp_1 .LT. lac_1) THEN
                          y15 = lac_1
                        ELSE
                          y15 = pmp_1
                        END IF
                      ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                        y15 = lac_1
                      ELSE
                        y15 = a4(1, i, k)
                      END IF
                      IF (x7 .GT. y15) THEN
                        a4(2, i, k) = y15
                      ELSE
                        a4(2, i, k) = x7
                      END IF
                      pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                      lac_2 = pmp_2 - 1.5*gam(i, k-1)
                      IF (a4(1, i, k) .GT. pmp_2) THEN
                        IF (pmp_2 .GT. lac_2) THEN
                          y26 = lac_2
                        ELSE
                          y26 = pmp_2
                        END IF
                      ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                        y26 = lac_2
                      ELSE
                        y26 = a4(1, i, k)
                      END IF
                      IF (a4(3, i, k) .LT. y26) THEN
                        x8 = y26
                      ELSE
                        x8 = a4(3, i, k)
                      END IF
                      IF (a4(1, i, k) .LT. pmp_2) THEN
                        IF (pmp_2 .LT. lac_2) THEN
                          y16 = lac_2
                        ELSE
                          y16 = pmp_2
                        END IF
                      ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                        y16 = lac_2
                      ELSE
                        y16 = a4(1, i, k)
                      END IF
                      IF (x8 .GT. y16) THEN
                        a4(3, i, k) = y16
                      ELSE
                        a4(3, i, k) = x8
                      END IF
                      a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(&
&                       3, i, k))
                    END IF
                  END IF
                END DO
              ELSE
                IF (kord .GE. 0.) THEN
                  abs8 = kord
                ELSE
                  abs8 = -kord
                END IF
                IF (abs8 .EQ. 13) THEN
                  DO i=i1,i2
                    IF (extm(i, k)) THEN
                      IF (extm(i, k-1) .AND. extm(i, k+1)) THEN
! grid-scale 2-delta-z wave detected
                        a4(2, i, k) = a4(1, i, k)
                        a4(3, i, k) = a4(1, i, k)
                        a4(4, i, k) = 0.
                      ELSE
! Left  edges
                        pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                        lac_1 = pmp_1 + 1.5*gam(i, k+2)
                        IF (a4(1, i, k) .GT. pmp_1) THEN
                          IF (pmp_1 .GT. lac_1) THEN
                            y27 = lac_1
                          ELSE
                            y27 = pmp_1
                          END IF
                        ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                          y27 = lac_1
                        ELSE
                          y27 = a4(1, i, k)
                        END IF
                        IF (a4(2, i, k) .LT. y27) THEN
                          x9 = y27
                        ELSE
                          x9 = a4(2, i, k)
                        END IF
                        IF (a4(1, i, k) .LT. pmp_1) THEN
                          IF (pmp_1 .LT. lac_1) THEN
                            y17 = lac_1
                          ELSE
                            y17 = pmp_1
                          END IF
                        ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                          y17 = lac_1
                        ELSE
                          y17 = a4(1, i, k)
                        END IF
                        IF (x9 .GT. y17) THEN
                          a4(2, i, k) = y17
                        ELSE
                          a4(2, i, k) = x9
                        END IF
! Right edges
                        pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                        lac_2 = pmp_2 - 1.5*gam(i, k-1)
                        IF (a4(1, i, k) .GT. pmp_2) THEN
                          IF (pmp_2 .GT. lac_2) THEN
                            y28 = lac_2
                          ELSE
                            y28 = pmp_2
                          END IF
                        ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                          y28 = lac_2
                        ELSE
                          y28 = a4(1, i, k)
                        END IF
                        IF (a4(3, i, k) .LT. y28) THEN
                          x10 = y28
                        ELSE
                          x10 = a4(3, i, k)
                        END IF
                        IF (a4(1, i, k) .LT. pmp_2) THEN
                          IF (pmp_2 .LT. lac_2) THEN
                            y18 = lac_2
                          ELSE
                            y18 = pmp_2
                          END IF
                        ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                          y18 = lac_2
                        ELSE
                          y18 = a4(1, i, k)
                        END IF
                        IF (x10 .GT. y18) THEN
                          a4(3, i, k) = y18
                        ELSE
                          a4(3, i, k) = x10
                        END IF
                        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4&
&                         (3, i, k)))
                      END IF
                    ELSE
                      a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3&
&                       , i, k)))
                    END IF
                  END DO
                ELSE
                  IF (kord .GE. 0.) THEN
                    abs9 = kord
                  ELSE
                    abs9 = -kord
                  END IF
                  IF (abs9 .EQ. 14) THEN
                    DO i=i1,i2
                      a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3&
&                       , i, k)))
                    END DO
                  ELSE
! kord = 11
                    DO i=i1,i2
                      IF (extm(i, k) .AND. (extm(i, k-1) .OR. extm(i, k+&
&                         1))) THEN
! Noisy region:
                        a4(2, i, k) = a4(1, i, k)
                        a4(3, i, k) = a4(1, i, k)
                        a4(4, i, k) = 0.
                      ELSE
                        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4&
&                         (3, i, k)))
                      END IF
                    END DO
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
! Additional constraint to ensure positivity
        IF (iv .EQ. 0) CALL CS_LIMITERS(im, extm(i1, k), a4(1, i1, k), 0&
&                                )
      END DO
! k-loop
!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
      IF (iv .EQ. 0) THEN
        DO i=i1,i2
          IF (0. .LT. a4(3, i, km)) THEN
            a4(3, i, km) = a4(3, i, km)
          ELSE
            a4(3, i, km) = 0.
          END IF
        END DO
      ELSE IF (iv .EQ. -1) THEN
        DO i=i1,i2
          IF (a4(3, i, km)*a4(1, i, km) .LE. 0.) a4(3, i, km) = 0.
        END DO
      END IF
      DO k=km-1,km
        DO i=i1,i2
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
        IF (k .EQ. km - 1) CALL CS_LIMITERS(im, extm(i1, k), a4(1, i1, k&
&                                     ), 2)
        IF (k .EQ. km) CALL CS_LIMITERS(im, extm(i1, k), a4(1, i1, k), 1&
&                                )
      END DO
    END IF
  END SUBROUTINE CS_PROFILE
  SUBROUTINE CS_LIMITERS(im, extm, a4, iv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    INTEGER, INTENT(IN) :: iv
    LOGICAL, INTENT(IN) :: extm(im)
! PPM array
    REAL, INTENT(INOUT) :: a4(4, im)
! !LOCAL VARIABLES:
    REAL :: da1, da2, a6da
    INTEGER :: i
    INTRINSIC ABS
    REAL :: abs0
    IF (iv .EQ. 0) THEN
! Positive definite constraint
      DO i=1,im
        IF (a4(1, i) .LE. 0.) THEN
          a4(2, i) = a4(1, i)
          a4(3, i) = a4(1, i)
          a4(4, i) = 0.
        ELSE
          IF (a4(3, i) - a4(2, i) .GE. 0.) THEN
            abs0 = a4(3, i) - a4(2, i)
          ELSE
            abs0 = -(a4(3, i)-a4(2, i))
          END IF
          IF (abs0 .LT. -a4(4, i)) THEN
            IF (a4(1, i) + 0.25*(a4(3, i)-a4(2, i))**2/a4(4, i) + a4(4, &
&               i)*r12 .LT. 0.) THEN
! local minimum is negative
              IF (a4(1, i) .LT. a4(3, i) .AND. a4(1, i) .LT. a4(2, i)) &
&             THEN
                a4(3, i) = a4(1, i)
                a4(2, i) = a4(1, i)
                a4(4, i) = 0.
              ELSE IF (a4(3, i) .GT. a4(2, i)) THEN
                a4(4, i) = 3.*(a4(2, i)-a4(1, i))
                a4(3, i) = a4(2, i) - a4(4, i)
              ELSE
                a4(4, i) = 3.*(a4(3, i)-a4(1, i))
                a4(2, i) = a4(3, i) - a4(4, i)
              END IF
            END IF
          END IF
        END IF
      END DO
    ELSE IF (iv .EQ. 1) THEN
      DO i=1,im
        IF ((a4(1, i)-a4(2, i))*(a4(1, i)-a4(3, i)) .GE. 0.) THEN
          a4(2, i) = a4(1, i)
          a4(3, i) = a4(1, i)
          a4(4, i) = 0.
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            a4(3, i) = a4(2, i) - a4(4, i)
          ELSE IF (a6da .GT. da2) THEN
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            a4(2, i) = a4(3, i) - a4(4, i)
          END IF
        END IF
      END DO
    ELSE
! Standard PPM constraint
      DO i=1,im
        IF (extm(i)) THEN
          a4(2, i) = a4(1, i)
          a4(3, i) = a4(1, i)
          a4(4, i) = 0.
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            a4(3, i) = a4(2, i) - a4(4, i)
          ELSE IF (a6da .GT. da2) THEN
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            a4(2, i) = a4(3, i) - a4(4, i)
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE CS_LIMITERS
  SUBROUTINE PPM_PROFILE(a4, delp, km, i1, i2, iv, kord)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! iv =-1: winds
    INTEGER, INTENT(IN) :: iv
! iv = 0: positive definite scalars
! iv = 1: others
! iv = 2: w (iv=-2)
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! vertical dimension
    INTEGER, INTENT(IN) :: km
! Order (or more accurately method no.):
    INTEGER, INTENT(IN) :: kord
!
! layer pressure thickness
    REAL, INTENT(IN) :: delp(i1:i2, km)
! !INPUT/OUTPUT PARAMETERS:
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
!
! !REVISION HISTORY:
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
    REAL :: dc(i1:i2, km)
    REAL :: h2(i1:i2, km)
    REAL :: delq(i1:i2, km)
    REAL :: df2(i1:i2, km)
    REAL :: d4(i1:i2, km)
! local scalars:
    INTEGER :: i, k, km1, lmt, it
    REAL :: fac
    REAL :: a1, a2, c1, c2, c3, d1, d2
    REAL :: qm, dq, lac, qmp, pmp
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL :: min1
    INTEGER :: abs0
    REAL :: max1
    REAL :: min2
    REAL :: x3
    REAL :: x2
    REAL :: x1
    REAL :: z1
    REAL :: y9
    REAL :: y8
    REAL :: y7
    REAL :: y6
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: y2
    REAL :: y1
    km1 = km - 1
    it = i2 - i1 + 1
    DO k=2,km
      DO i=i1,i2
        delq(i, k-1) = a4(1, i, k) - a4(1, i, k-1)
        d4(i, k) = delp(i, k-1) + delp(i, k)
      END DO
    END DO
    DO k=2,km1
      DO i=i1,i2
        c1 = (delp(i, k-1)+0.5*delp(i, k))/d4(i, k+1)
        c2 = (delp(i, k+1)+0.5*delp(i, k))/d4(i, k)
        df2(i, k) = delp(i, k)*(c1*delq(i, k)+c2*delq(i, k-1))/(d4(i, k)&
&         +delp(i, k+1))
        IF (df2(i, k) .GE. 0.) THEN
          x1 = df2(i, k)
        ELSE
          x1 = -df2(i, k)
        END IF
        IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
          IF (a4(1, i, k) .LT. a4(1, i, k+1)) THEN
            max1 = a4(1, i, k+1)
          ELSE
            max1 = a4(1, i, k)
          END IF
        ELSE IF (a4(1, i, k-1) .LT. a4(1, i, k+1)) THEN
          max1 = a4(1, i, k+1)
        ELSE
          max1 = a4(1, i, k-1)
        END IF
        y1 = max1 - a4(1, i, k)
        IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
          IF (a4(1, i, k) .GT. a4(1, i, k+1)) THEN
            min2 = a4(1, i, k+1)
          ELSE
            min2 = a4(1, i, k)
          END IF
        ELSE IF (a4(1, i, k-1) .GT. a4(1, i, k+1)) THEN
          min2 = a4(1, i, k+1)
        ELSE
          min2 = a4(1, i, k-1)
        END IF
        z1 = a4(1, i, k) - min2
        IF (x1 .GT. y1) THEN
          IF (y1 .GT. z1) THEN
            min1 = z1
          ELSE
            min1 = y1
          END IF
        ELSE IF (x1 .GT. z1) THEN
          min1 = z1
        ELSE
          min1 = x1
        END IF
        dc(i, k) = SIGN(min1, df2(i, k))
      END DO
    END DO
!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------
    DO k=3,km1
      DO i=i1,i2
        c1 = delq(i, k-1)*delp(i, k-1)/d4(i, k)
        a1 = d4(i, k-1)/(d4(i, k)+delp(i, k-1))
        a2 = d4(i, k+1)/(d4(i, k)+delp(i, k))
        a4(2, i, k) = a4(1, i, k-1) + c1 + 2./(d4(i, k-1)+d4(i, k+1))*(&
&         delp(i, k)*(c1*(a1-a2)+a2*dc(i, k-1))-delp(i, k-1)*a1*dc(i, k)&
&         )
      END DO
    END DO
!     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)
! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
    DO i=i1,i2
      d1 = delp(i, 1)
      d2 = delp(i, 2)
      qm = (d2*a4(1, i, 1)+d1*a4(1, i, 2))/(d1+d2)
      dq = 2.*(a4(1, i, 2)-a4(1, i, 1))/(d1+d2)
      c1 = 4.*(a4(2, i, 3)-qm-d2*dq)/(d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      a4(2, i, 2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
      a4(2, i, 1) = d1*(2.*c1*d1**2-c3) + a4(2, i, 2)
      IF (a4(1, i, 1) .GT. a4(1, i, 2)) THEN
        y2 = a4(1, i, 2)
      ELSE
        y2 = a4(1, i, 1)
      END IF
      IF (a4(2, i, 2) .LT. y2) THEN
        a4(2, i, 2) = y2
      ELSE
        a4(2, i, 2) = a4(2, i, 2)
      END IF
      IF (a4(1, i, 1) .LT. a4(1, i, 2)) THEN
        y3 = a4(1, i, 2)
      ELSE
        y3 = a4(1, i, 1)
      END IF
      IF (a4(2, i, 2) .GT. y3) THEN
        a4(2, i, 2) = y3
      ELSE
        a4(2, i, 2) = a4(2, i, 2)
      END IF
      dc(i, 1) = 0.5*(a4(2, i, 2)-a4(1, i, 1))
    END DO
! Enforce monotonicity  within the top layer
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(2, i, 1)) THEN
          a4(2, i, 1) = a4(2, i, 1)
        ELSE
          a4(2, i, 1) = 0.
        END IF
        IF (0. .LT. a4(2, i, 2)) THEN
          a4(2, i, 2) = a4(2, i, 2)
        ELSE
          a4(2, i, 2) = 0.
        END IF
      END DO
    ELSE IF (iv .EQ. -1) THEN
      DO i=i1,i2
        IF (a4(2, i, 1)*a4(1, i, 1) .LE. 0.) a4(2, i, 1) = 0.
      END DO
    ELSE
      IF (iv .GE. 0.) THEN
        abs0 = iv
      ELSE
        abs0 = -iv
      END IF
      IF (abs0 .EQ. 2) THEN
        DO i=i1,i2
          a4(2, i, 1) = a4(1, i, 1)
          a4(3, i, 1) = a4(1, i, 1)
        END DO
      END IF
    END IF
! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
    DO i=i1,i2
      d1 = delp(i, km)
      d2 = delp(i, km1)
      qm = (d2*a4(1, i, km)+d1*a4(1, i, km1))/(d1+d2)
      dq = 2.*(a4(1, i, km1)-a4(1, i, km))/(d1+d2)
      c1 = (a4(2, i, km1)-qm-d2*dq)/(d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      a4(2, i, km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
      a4(3, i, km) = d1*(8.*c1*d1**2-c3) + a4(2, i, km)
      IF (a4(1, i, km) .GT. a4(1, i, km1)) THEN
        y4 = a4(1, i, km1)
      ELSE
        y4 = a4(1, i, km)
      END IF
      IF (a4(2, i, km) .LT. y4) THEN
        a4(2, i, km) = y4
      ELSE
        a4(2, i, km) = a4(2, i, km)
      END IF
      IF (a4(1, i, km) .LT. a4(1, i, km1)) THEN
        y5 = a4(1, i, km1)
      ELSE
        y5 = a4(1, i, km)
      END IF
      IF (a4(2, i, km) .GT. y5) THEN
        a4(2, i, km) = y5
      ELSE
        a4(2, i, km) = a4(2, i, km)
      END IF
      dc(i, km) = 0.5*(a4(1, i, km)-a4(2, i, km))
    END DO
! Enforce constraint on the "slope" at the surface
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(2, i, km)) THEN
          a4(2, i, km) = a4(2, i, km)
        ELSE
          a4(2, i, km) = 0.
        END IF
        IF (0. .LT. a4(3, i, km)) THEN
          a4(3, i, km) = a4(3, i, km)
        ELSE
          a4(3, i, km) = 0.
        END IF
      END DO
    ELSE IF (iv .LT. 0) THEN
      DO i=i1,i2
        IF (a4(1, i, km)*a4(3, i, km) .LE. 0.) a4(3, i, km) = 0.
      END DO
    END IF
    DO k=1,km1
      DO i=i1,i2
        a4(3, i, k) = a4(2, i, k+1)
      END DO
    END DO
!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
    DO k=1,2
      DO i=i1,i2
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
      CALL PPM_LIMITERS(dc(i1, k), a4(1, i1, k), it, 0)
    END DO
    IF (kord .GE. 7) THEN
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      DO k=2,km1
        DO i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
          h2(i, k) = 2.*(dc(i, k+1)/delp(i, k+1)-dc(i, k-1)/delp(i, k-1)&
&           )/(delp(i, k)+0.5*(delp(i, k-1)+delp(i, k+1)))*delp(i, k)**2
        END DO
      END DO
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
! original quasi-monotone
      fac = 1.5
      DO k=3,km-2
        DO i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
          pmp = 2.*dc(i, k)
          qmp = a4(1, i, k) + pmp
          lac = a4(1, i, k) + fac*h2(i, k-1) + dc(i, k)
          IF (a4(1, i, k) .GT. qmp) THEN
            IF (qmp .GT. lac) THEN
              y8 = lac
            ELSE
              y8 = qmp
            END IF
          ELSE IF (a4(1, i, k) .GT. lac) THEN
            y8 = lac
          ELSE
            y8 = a4(1, i, k)
          END IF
          IF (a4(3, i, k) .LT. y8) THEN
            x2 = y8
          ELSE
            x2 = a4(3, i, k)
          END IF
          IF (a4(1, i, k) .LT. qmp) THEN
            IF (qmp .LT. lac) THEN
              y6 = lac
            ELSE
              y6 = qmp
            END IF
          ELSE IF (a4(1, i, k) .LT. lac) THEN
            y6 = lac
          ELSE
            y6 = a4(1, i, k)
          END IF
          IF (x2 .GT. y6) THEN
            a4(3, i, k) = y6
          ELSE
            a4(3, i, k) = x2
          END IF
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
          qmp = a4(1, i, k) - pmp
          lac = a4(1, i, k) + fac*h2(i, k+1) - dc(i, k)
          IF (a4(1, i, k) .GT. qmp) THEN
            IF (qmp .GT. lac) THEN
              y9 = lac
            ELSE
              y9 = qmp
            END IF
          ELSE IF (a4(1, i, k) .GT. lac) THEN
            y9 = lac
          ELSE
            y9 = a4(1, i, k)
          END IF
          IF (a4(2, i, k) .LT. y9) THEN
            x3 = y9
          ELSE
            x3 = a4(2, i, k)
          END IF
          IF (a4(1, i, k) .LT. qmp) THEN
            IF (qmp .LT. lac) THEN
              y7 = lac
            ELSE
              y7 = qmp
            END IF
          ELSE IF (a4(1, i, k) .LT. lac) THEN
            y7 = lac
          ELSE
            y7 = a4(1, i, k)
          END IF
          IF (x3 .GT. y7) THEN
            a4(2, i, k) = y7
          ELSE
            a4(2, i, k) = x3
          END IF
!-------------
! Recompute A6
!-------------
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
! Additional constraint to ensure positivity when kord=7
        IF (iv .EQ. 0 .AND. kord .GE. 6) CALL PPM_LIMITERS(dc(i1, k), a4&
&                                                    (1, i1, k), it, 2)
      END DO
    ELSE
      lmt = kord - 3
      IF (0 .LT. lmt) THEN
        lmt = lmt
      ELSE
        lmt = 0
      END IF
      IF (iv .EQ. 0) THEN
        IF (2 .GT. lmt) THEN
          lmt = lmt
        ELSE
          lmt = 2
        END IF
      END IF
      DO k=3,km-2
        IF (kord .NE. 4) THEN
          DO i=i1,i2
            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
          END DO
        END IF
        IF (kord .NE. 6) CALL PPM_LIMITERS(dc(i1, k), a4(1, i1, k), it, &
&                                    lmt)
      END DO
    END IF
    DO k=km1,km
      DO i=i1,i2
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
      CALL PPM_LIMITERS(dc(i1, k), a4(1, i1, k), it, 0)
    END DO
  END SUBROUTINE PPM_PROFILE
  SUBROUTINE PPM_LIMITERS(dm, a4, itot, lmt)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! the linear slope
    REAL, INTENT(IN) :: dm(*)
! Total Longitudes
    INTEGER, INTENT(IN) :: itot
! 0: Standard PPM constraint
    INTEGER, INTENT(IN) :: lmt
! 1: Improved full monotonicity constraint (Lin)
! 2: Positive definite constraint
! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
! PPM array
    REAL, INTENT(INOUT) :: a4(4, *)
! AA <-- a4(1,i)
! AL <-- a4(2,i)
! AR <-- a4(3,i)
! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
    REAL :: qmp
    REAL :: da1, da2, a6da
    REAL :: fmin
    INTEGER :: i
    INTRINSIC ABS
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min2
    REAL :: abs0
    REAL :: x2
    REAL :: x1
    REAL :: y2
    REAL :: y1
! Developer: S.-J. Lin
    IF (lmt .EQ. 3) THEN
      RETURN
    ELSE IF (lmt .EQ. 0) THEN
! Standard PPM constraint
      DO i=1,itot
        IF (dm(i) .EQ. 0.) THEN
          a4(2, i) = a4(1, i)
          a4(3, i) = a4(1, i)
          a4(4, i) = 0.
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            a4(3, i) = a4(2, i) - a4(4, i)
          ELSE IF (a6da .GT. da2) THEN
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            a4(2, i) = a4(3, i) - a4(4, i)
          END IF
        END IF
      END DO
    ELSE IF (lmt .EQ. 1) THEN
! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      DO i=1,itot
        qmp = 2.*dm(i)
        IF (qmp .GE. 0.) THEN
          x1 = qmp
        ELSE
          x1 = -qmp
        END IF
        IF (a4(2, i) - a4(1, i) .GE. 0.) THEN
          y1 = a4(2, i) - a4(1, i)
        ELSE
          y1 = -(a4(2, i)-a4(1, i))
        END IF
        IF (x1 .GT. y1) THEN
          min1 = y1
        ELSE
          min1 = x1
        END IF
        a4(2, i) = a4(1, i) - SIGN(min1, qmp)
        IF (qmp .GE. 0.) THEN
          x2 = qmp
        ELSE
          x2 = -qmp
        END IF
        IF (a4(3, i) - a4(1, i) .GE. 0.) THEN
          y2 = a4(3, i) - a4(1, i)
        ELSE
          y2 = -(a4(3, i)-a4(1, i))
        END IF
        IF (x2 .GT. y2) THEN
          min2 = y2
        ELSE
          min2 = x2
        END IF
        a4(3, i) = a4(1, i) + SIGN(min2, qmp)
        a4(4, i) = 3.*(2.*a4(1, i)-(a4(2, i)+a4(3, i)))
      END DO
    ELSE IF (lmt .EQ. 2) THEN
! Positive definite constraint
      DO i=1,itot
        IF (a4(3, i) - a4(2, i) .GE. 0.) THEN
          abs0 = a4(3, i) - a4(2, i)
        ELSE
          abs0 = -(a4(3, i)-a4(2, i))
        END IF
        IF (abs0 .LT. -a4(4, i)) THEN
          fmin = a4(1, i) + 0.25*(a4(3, i)-a4(2, i))**2/a4(4, i) + a4(4&
&           , i)*r12
          IF (fmin .LT. 0.) THEN
            IF (a4(1, i) .LT. a4(3, i) .AND. a4(1, i) .LT. a4(2, i)) &
&           THEN
              a4(3, i) = a4(1, i)
              a4(2, i) = a4(1, i)
              a4(4, i) = 0.
            ELSE IF (a4(3, i) .GT. a4(2, i)) THEN
              a4(4, i) = 3.*(a4(2, i)-a4(1, i))
              a4(3, i) = a4(2, i) - a4(4, i)
            ELSE
              a4(4, i) = 3.*(a4(3, i)-a4(1, i))
              a4(2, i) = a4(3, i) - a4(4, i)
            END IF
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE PPM_LIMITERS
  SUBROUTINE MOIST_CV(is, ie, isd, ied, jsd, jed, km, j, k, nwat, sphum&
&   , liq_wat, rainwat, ice_wat, snowwat, graupel, q, qd, cvm, t1)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, isd, ied, jsd, jed, km, nwat, j, k
    INTEGER, INTENT(IN) :: sphum, liq_wat, rainwat, ice_wat, snowwat, &
&   graupel
    REAL, DIMENSION(isd:ied, jsd:jed, km, nwat), INTENT(IN) :: q
    REAL, DIMENSION(is:ie), INTENT(OUT) :: cvm, qd
    REAL, INTENT(IN), OPTIONAL :: t1(is:ie)
!
    REAL, PARAMETER :: t_i0=15.
    REAL, DIMENSION(is:ie) :: qv, ql, qs
    INTEGER :: i
    INTRINSIC PRESENT
    INTRINSIC MAX
    SELECT CASE  (nwat) 
    CASE (2) 
      IF (PRESENT(t1)) THEN
! Special case for GFS physics
        DO i=is,ie
          IF (0. .LT. q(i, j, k, liq_wat)) THEN
            qd(i) = q(i, j, k, liq_wat)
          ELSE
            qd(i) = 0.
          END IF
          IF (t1(i) .GT. tice) THEN
            qs(i) = 0.
          ELSE IF (t1(i) .LT. tice - t_i0) THEN
            qs(i) = qd(i)
          ELSE
            qs(i) = qd(i)*(tice-t1(i))/t_i0
          END IF
          ql(i) = qd(i) - qs(i)
          IF (0. .LT. q(i, j, k, sphum)) THEN
            qv(i) = q(i, j, k, sphum)
          ELSE
            qv(i) = 0.
          END IF
          cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*&
&           c_liq + qs(i)*c_ice
        END DO
      ELSE
        DO i=is,ie
          IF (0. .LT. q(i, j, k, sphum)) THEN
            qv(i) = q(i, j, k, sphum)
          ELSE
            qv(i) = 0.
          END IF
          IF (0. .LT. q(i, j, k, liq_wat)) THEN
            qs(i) = q(i, j, k, liq_wat)
          ELSE
            qs(i) = 0.
          END IF
          qd(i) = qs(i)
          cvm(i) = (1.-qv(i))*cv_air + qv(i)*cv_vap
        END DO
      END IF
    CASE (3) 
      DO i=is,ie
        qv(i) = q(i, j, k, sphum)
        ql(i) = q(i, j, k, liq_wat)
        qs(i) = q(i, j, k, ice_wat)
        qd(i) = ql(i) + qs(i)
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq &
&         + qs(i)*c_ice
      END DO
    CASE (4) 
! K_warm_rain with fake ice
      DO i=is,ie
        qv(i) = q(i, j, k, sphum)
        qd(i) = q(i, j, k, liq_wat) + q(i, j, k, rainwat)
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + qd(i)*c_liq
      END DO
    CASE (6) 
      DO i=is,ie
        qv(i) = q(i, j, k, sphum)
        ql(i) = q(i, j, k, liq_wat) + q(i, j, k, rainwat)
        qs(i) = q(i, j, k, ice_wat) + q(i, j, k, snowwat) + q(i, j, k, &
&         graupel)
        qd(i) = ql(i) + qs(i)
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq &
&         + qs(i)*c_ice
      END DO
    CASE DEFAULT
      DO i=is,ie
        qd(i) = 0.
        cvm(i) = cv_air
      END DO
    END SELECT
  END SUBROUTINE MOIST_CV
  SUBROUTINE MOIST_CP(is, ie, isd, ied, jsd, jed, km, j, k, nwat, sphum&
&   , liq_wat, rainwat, ice_wat, snowwat, graupel, q, qd, cpm, t1)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, isd, ied, jsd, jed, km, nwat, j, k
    INTEGER, INTENT(IN) :: sphum, liq_wat, rainwat, ice_wat, snowwat, &
&   graupel
    REAL, DIMENSION(isd:ied, jsd:jed, km, nwat), INTENT(IN) :: q
    REAL, DIMENSION(is:ie), INTENT(OUT) :: cpm, qd
    REAL, INTENT(IN), OPTIONAL :: t1(is:ie)
!
    REAL, PARAMETER :: t_i0=15.
    REAL, DIMENSION(is:ie) :: qv, ql, qs
    INTEGER :: i
    INTRINSIC PRESENT
    INTRINSIC MAX
    SELECT CASE  (nwat) 
    CASE (2) 
      IF (PRESENT(t1)) THEN
! Special case for GFS physics
        DO i=is,ie
          IF (0. .LT. q(i, j, k, liq_wat)) THEN
            qd(i) = q(i, j, k, liq_wat)
          ELSE
            qd(i) = 0.
          END IF
          IF (t1(i) .GT. tice) THEN
            qs(i) = 0.
          ELSE IF (t1(i) .LT. tice - t_i0) THEN
            qs(i) = qd(i)
          ELSE
            qs(i) = qd(i)*(tice-t1(i))/t_i0
          END IF
          ql(i) = qd(i) - qs(i)
          IF (0. .LT. q(i, j, k, sphum)) THEN
            qv(i) = q(i, j, k, sphum)
          ELSE
            qv(i) = 0.
          END IF
          cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*&
&           c_liq + qs(i)*c_ice
        END DO
      ELSE
        DO i=is,ie
          IF (0. .LT. q(i, j, k, sphum)) THEN
            qv(i) = q(i, j, k, sphum)
          ELSE
            qv(i) = 0.
          END IF
          IF (0. .LT. q(i, j, k, liq_wat)) THEN
            qs(i) = q(i, j, k, liq_wat)
          ELSE
            qs(i) = 0.
          END IF
          qd(i) = qs(i)
          cpm(i) = (1.-qv(i))*cp_air + qv(i)*cp_vapor
        END DO
      END IF
    CASE (3) 
      DO i=is,ie
        qv(i) = q(i, j, k, sphum)
        ql(i) = q(i, j, k, liq_wat)
        qs(i) = q(i, j, k, ice_wat)
        qd(i) = ql(i) + qs(i)
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*&
&         c_liq + qs(i)*c_ice
      END DO
    CASE (4) 
! K_warm_rain scheme with fake ice
      DO i=is,ie
        qv(i) = q(i, j, k, sphum)
        qd(i) = q(i, j, k, liq_wat) + q(i, j, k, rainwat)
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + qd(i)*&
&         c_liq
      END DO
    CASE (6) 
      DO i=is,ie
        qv(i) = q(i, j, k, sphum)
        ql(i) = q(i, j, k, liq_wat) + q(i, j, k, rainwat)
        qs(i) = q(i, j, k, ice_wat) + q(i, j, k, snowwat) + q(i, j, k, &
&         graupel)
        qd(i) = ql(i) + qs(i)
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*&
&         c_liq + qs(i)*c_ice
      END DO
    CASE DEFAULT
      DO i=is,ie
        qd(i) = 0.
        cpm(i) = cp_air
      END DO
    END SELECT
  END SUBROUTINE MOIST_CP
!  Differentiation of map1_cubic in forward (tangent) mode:
!   variations   of useful results: q2
!   with respect to varying inputs: pe1 pe2 q2
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  map1_cubic --- Cubic Interpolation for vertical re-mapping
!
! !INTERFACE:
  SUBROUTINE MAP1_CUBIC_TLM(km, pe1, pe1_tl, kn, pe2, pe2_tl, q2, q2_tl&
&   , i1, i2, j, ibeg, iend, jbeg, jend, akap, t_var, conserv)
    IMPLICIT NONE
!EOC
! !INPUT PARAMETERS:
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
    REAL, INTENT(IN) :: akap
! Thermodynamic variable to remap
    INTEGER, INTENT(IN) :: t_var
!     1:TE  2:T  3:PT
    LOGICAL, INTENT(IN) :: conserv
! Current latitude
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL, INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL, INTENT(IN) :: pe2_tl(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
!      real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(INOUT) :: q2_tl(ibeg:iend, jbeg:jend, kn)
! !DESCRIPTION:
!
!     Perform Cubic Interpolation a given latitude
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY:
!    2005.11.14   Takacs    Initial Code
!    2016.07.20   Putman    Modified to make genaric for any thermodynamic variable
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    REAL :: qx(i1:i2, km)
    REAL :: qx_tl(i1:i2, km)
    REAL :: logpl1(i1:i2, km)
    REAL :: logpl1_tl(i1:i2, km)
    REAL :: logpl2(i1:i2, kn)
    REAL :: logpl2_tl(i1:i2, kn)
    REAL :: dlogp1(i1:i2, km)
    REAL :: dlogp1_tl(i1:i2, km)
    REAL :: vsum1(i1:i2)
    REAL :: vsum1_tl(i1:i2)
    REAL :: vsum2(i1:i2)
    REAL :: vsum2_tl(i1:i2)
    REAL :: am2, am1, ap0, ap1, p, plp1, plp0, plm1, plm2, dlp0, dlm1, &
&   dlm2
    REAL :: am2_tl, am1_tl, ap0_tl, ap1_tl, p_tl, plp1_tl, plp0_tl, &
&   plm1_tl, plm2_tl, dlp0_tl, dlm1_tl, dlm2_tl
    INTEGER :: i, k, lm2, lm1, lp0, lp1
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    INTRINSIC MIN
    REAL, DIMENSION(i2-i1+1) :: arg1
    REAL, DIMENSION(i2-i1+1) :: arg1_tl
    REAL, DIMENSION(i1:i2) :: arg2
    REAL, DIMENSION(i1:i2) :: arg2_tl
! Initialization
! --------------
    SELECT CASE  (t_var) 
    CASE (1) 
      qx_tl = 0.0
      logpl1_tl = 0.0
! Total Energy Remapping in Log(P)
      DO k=1,km
        qx_tl(:, k) = q2_tl(i1:i2, j, k)
        qx(:, k) = q2(i1:i2, j, k)
        logpl1_tl(:, k) = (pe1_tl(:, k)+pe1_tl(:, k+1))/(pe1(:, k)+pe1(:&
&         , k+1))
        logpl1(:, k) = LOG(r2*(pe1(:, k)+pe1(:, k+1)))
      END DO
      logpl2_tl = 0.0
      DO k=1,kn
        logpl2_tl(:, k) = (pe2_tl(:, k)+pe2_tl(:, k+1))/(pe2(:, k)+pe2(:&
&         , k+1))
        logpl2(:, k) = LOG(r2*(pe2(:, k)+pe2(:, k+1)))
      END DO
      dlogp1_tl = 0.0
      DO k=1,km-1
        dlogp1_tl(:, k) = logpl1_tl(:, k+1) - logpl1_tl(:, k)
        dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
      END DO
    CASE (2) 
      qx_tl = 0.0
      logpl1_tl = 0.0
! Temperature Remapping in Log(P)
      DO k=1,km
        qx_tl(:, k) = q2_tl(i1:i2, j, k)
        qx(:, k) = q2(i1:i2, j, k)
        logpl1_tl(:, k) = (pe1_tl(:, k)+pe1_tl(:, k+1))/(pe1(:, k)+pe1(:&
&         , k+1))
        logpl1(:, k) = LOG(r2*(pe1(:, k)+pe1(:, k+1)))
      END DO
      logpl2_tl = 0.0
      DO k=1,kn
        logpl2_tl(:, k) = (pe2_tl(:, k)+pe2_tl(:, k+1))/(pe2(:, k)+pe2(:&
&         , k+1))
        logpl2(:, k) = LOG(r2*(pe2(:, k)+pe2(:, k+1)))
      END DO
      dlogp1_tl = 0.0
      DO k=1,km-1
        dlogp1_tl(:, k) = logpl1_tl(:, k+1) - logpl1_tl(:, k)
        dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
      END DO
    CASE (3) 
      qx_tl = 0.0
      logpl1_tl = 0.0
! Potential Temperature Remapping in P^KAPPA
      DO k=1,km
        qx_tl(:, k) = q2_tl(i1:i2, j, k)
        qx(:, k) = q2(i1:i2, j, k)
        arg1_tl(:) = r2*(pe1_tl(:, k)+pe1_tl(:, k+1))
        arg1(:) = r2*(pe1(:, k)+pe1(:, k+1))
        arg2_tl(:) = akap*arg1_tl(:)/arg1(:)
        arg2(:) = akap*LOG(arg1(:))
        logpl1_tl(:, k) = arg2_tl(:)*EXP(arg2(:))
        logpl1(:, k) = EXP(arg2(:))
      END DO
      logpl2_tl = 0.0
      DO k=1,kn
        arg1_tl(:) = r2*(pe2_tl(:, k)+pe2_tl(:, k+1))
        arg1(:) = r2*(pe2(:, k)+pe2(:, k+1))
        arg2_tl(:) = akap*arg1_tl(:)/arg1(:)
        arg2(:) = akap*LOG(arg1(:))
        logpl2_tl(:, k) = arg2_tl(:)*EXP(arg2(:))
        logpl2(:, k) = EXP(arg2(:))
      END DO
      dlogp1_tl = 0.0
      DO k=1,km-1
        dlogp1_tl(:, k) = logpl1_tl(:, k+1) - logpl1_tl(:, k)
        dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
      END DO
    CASE DEFAULT
      qx_tl = 0.0
      logpl1_tl = 0.0
      logpl2_tl = 0.0
      dlogp1_tl = 0.0
    END SELECT
    IF (conserv) THEN
! Compute vertical integral of Input TE
! -------------------------------------
      vsum1(:) = r0
      vsum1_tl = 0.0
      DO i=i1,i2
        DO k=1,km
          vsum1_tl(i) = vsum1_tl(i) + qx_tl(i, k)*(pe1(i, k+1)-pe1(i, k)&
&           ) + qx(i, k)*(pe1_tl(i, k+1)-pe1_tl(i, k))
          vsum1(i) = vsum1(i) + qx(i, k)*(pe1(i, k+1)-pe1(i, k))
        END DO
        vsum1_tl(i) = (vsum1_tl(i)*(pe1(i, km+1)-pe1(i, 1))-vsum1(i)*(&
&         pe1_tl(i, km+1)-pe1_tl(i, 1)))/(pe1(i, km+1)-pe1(i, 1))**2
        vsum1(i) = vsum1(i)/(pe1(i, km+1)-pe1(i, 1))
      END DO
    ELSE
      vsum1_tl = 0.0
    END IF
! Interpolate TE onto target Pressures
! ------------------------------------
    DO i=i1,i2
      DO k=1,kn
        lm1 = 1
        lp0 = 1
        DO WHILE (lp0 .LE. km)
          IF (logpl1(i, lp0) .LT. logpl2(i, k)) THEN
            lp0 = lp0 + 1
          ELSE
            EXIT
          END IF
        END DO
        IF (lp0 - 1 .LT. 1) THEN
          lm1 = 1
        ELSE
          lm1 = lp0 - 1
        END IF
        IF (lp0 .GT. km) THEN
          lp0 = km
        ELSE
          lp0 = lp0
        END IF
! Extrapolate Linearly in LogP above first model level
! ----------------------------------------------------
        IF (lm1 .EQ. 1 .AND. lp0 .EQ. 1) THEN
          q2_tl(i, j, k) = qx_tl(i, 1) + (((qx_tl(i, 2)-qx_tl(i, 1))*(&
&           logpl2(i, k)-logpl1(i, 1))+(qx(i, 2)-qx(i, 1))*(logpl2_tl(i&
&           , k)-logpl1_tl(i, 1)))*(logpl1(i, 2)-logpl1(i, 1))-(qx(i, 2)&
&           -qx(i, 1))*(logpl2(i, k)-logpl1(i, 1))*(logpl1_tl(i, 2)-&
&           logpl1_tl(i, 1)))/(logpl1(i, 2)-logpl1(i, 1))**2
          q2(i, j, k) = qx(i, 1) + (qx(i, 2)-qx(i, 1))*(logpl2(i, k)-&
&           logpl1(i, 1))/(logpl1(i, 2)-logpl1(i, 1))
! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
        ELSE IF (lm1 .EQ. km .AND. lp0 .EQ. km) THEN
          q2_tl(i, j, k) = qx_tl(i, km) + (((qx_tl(i, km)-qx_tl(i, km-1)&
&           )*(logpl2(i, k)-logpl1(i, km))+(qx(i, km)-qx(i, km-1))*(&
&           logpl2_tl(i, k)-logpl1_tl(i, km)))*(logpl1(i, km)-logpl1(i, &
&           km-1))-(qx(i, km)-qx(i, km-1))*(logpl2(i, k)-logpl1(i, km))*&
&           (logpl1_tl(i, km)-logpl1_tl(i, km-1)))/(logpl1(i, km)-logpl1&
&           (i, km-1))**2
          q2(i, j, k) = qx(i, km) + (qx(i, km)-qx(i, km-1))*(logpl2(i, k&
&           )-logpl1(i, km))/(logpl1(i, km)-logpl1(i, km-1))
! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
        ELSE IF (lm1 .EQ. 1 .OR. lp0 .EQ. km) THEN
          q2_tl(i, j, k) = qx_tl(i, lp0) + (((qx_tl(i, lm1)-qx_tl(i, lp0&
&           ))*(logpl2(i, k)-logpl1(i, lp0))+(qx(i, lm1)-qx(i, lp0))*(&
&           logpl2_tl(i, k)-logpl1_tl(i, lp0)))*(logpl1(i, lm1)-logpl1(i&
&           , lp0))-(qx(i, lm1)-qx(i, lp0))*(logpl2(i, k)-logpl1(i, lp0)&
&           )*(logpl1_tl(i, lm1)-logpl1_tl(i, lp0)))/(logpl1(i, lm1)-&
&           logpl1(i, lp0))**2
          q2(i, j, k) = qx(i, lp0) + (qx(i, lm1)-qx(i, lp0))*(logpl2(i, &
&           k)-logpl1(i, lp0))/(logpl1(i, lm1)-logpl1(i, lp0))
! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
        ELSE
          lp1 = lp0 + 1
          lm2 = lm1 - 1
          p_tl = logpl2_tl(i, k)
          p = logpl2(i, k)
          plp1_tl = logpl1_tl(i, lp1)
          plp1 = logpl1(i, lp1)
          plp0_tl = logpl1_tl(i, lp0)
          plp0 = logpl1(i, lp0)
          plm1_tl = logpl1_tl(i, lm1)
          plm1 = logpl1(i, lm1)
          plm2_tl = logpl1_tl(i, lm2)
          plm2 = logpl1(i, lm2)
          dlp0_tl = dlogp1_tl(i, lp0)
          dlp0 = dlogp1(i, lp0)
          dlm1_tl = dlogp1_tl(i, lm1)
          dlm1 = dlogp1(i, lm1)
          dlm2_tl = dlogp1_tl(i, lm2)
          dlm2 = dlogp1(i, lm2)
          ap1_tl = ((((p_tl-plp0_tl)*(p-plm1)+(p-plp0)*(p_tl-plm1_tl))*(&
&           p-plm2)+(p-plp0)*(p-plm1)*(p_tl-plm2_tl))*dlp0*(dlp0+dlm1)*(&
&           dlp0+dlm1+dlm2)-(p-plp0)*(p-plm1)*(p-plm2)*((dlp0_tl*(dlp0+&
&           dlm1)+dlp0*(dlp0_tl+dlm1_tl))*(dlp0+dlm1+dlm2)+dlp0*(dlp0+&
&           dlm1)*(dlp0_tl+dlm1_tl+dlm2_tl)))/(dlp0*(dlp0+dlm1)*(dlp0+&
&           dlm1+dlm2))**2
          ap1 = (p-plp0)*(p-plm1)*(p-plm2)/(dlp0*(dlp0+dlm1)*(dlp0+dlm1+&
&           dlm2))
          ap0_tl = ((((plp1_tl-p_tl)*(p-plm1)+(plp1-p)*(p_tl-plm1_tl))*(&
&           p-plm2)+(plp1-p)*(p-plm1)*(p_tl-plm2_tl))*dlp0*dlm1*(dlm1+&
&           dlm2)-(plp1-p)*(p-plm1)*(p-plm2)*((dlp0_tl*dlm1+dlp0*dlm1_tl&
&           )*(dlm1+dlm2)+dlp0*dlm1*(dlm1_tl+dlm2_tl)))/(dlp0*dlm1*(dlm1&
&           +dlm2))**2
          ap0 = (plp1-p)*(p-plm1)*(p-plm2)/(dlp0*dlm1*(dlm1+dlm2))
          am1_tl = ((((plp1_tl-p_tl)*(plp0-p)+(plp1-p)*(plp0_tl-p_tl))*(&
&           p-plm2)+(plp1-p)*(plp0-p)*(p_tl-plm2_tl))*dlm1*dlm2*(dlp0+&
&           dlm1)-(plp1-p)*(plp0-p)*(p-plm2)*((dlm1_tl*dlm2+dlm1*dlm2_tl&
&           )*(dlp0+dlm1)+dlm1*dlm2*(dlp0_tl+dlm1_tl)))/(dlm1*dlm2*(dlp0&
&           +dlm1))**2
          am1 = (plp1-p)*(plp0-p)*(p-plm2)/(dlm1*dlm2*(dlp0+dlm1))
          am2_tl = ((((plp1_tl-p_tl)*(plp0-p)+(plp1-p)*(plp0_tl-p_tl))*(&
&           plm1-p)+(plp1-p)*(plp0-p)*(plm1_tl-p_tl))*dlm2*(dlm1+dlm2)*(&
&           dlp0+dlm1+dlm2)-(plp1-p)*(plp0-p)*(plm1-p)*((dlm2_tl*(dlm1+&
&           dlm2)+dlm2*(dlm1_tl+dlm2_tl))*(dlp0+dlm1+dlm2)+dlm2*(dlm1+&
&           dlm2)*(dlp0_tl+dlm1_tl+dlm2_tl)))/(dlm2*(dlm1+dlm2)*(dlp0+&
&           dlm1+dlm2))**2
          am2 = (plp1-p)*(plp0-p)*(plm1-p)/(dlm2*(dlm1+dlm2)*(dlp0+dlm1+&
&           dlm2))
          q2_tl(i, j, k) = ap1_tl*qx(i, lp1) + ap1*qx_tl(i, lp1) + &
&           ap0_tl*qx(i, lp0) + ap0*qx_tl(i, lp0) + am1_tl*qx(i, lm1) + &
&           am1*qx_tl(i, lm1) + am2_tl*qx(i, lm2) + am2*qx_tl(i, lm2)
          q2(i, j, k) = ap1*qx(i, lp1) + ap0*qx(i, lp0) + am1*qx(i, lm1)&
&           + am2*qx(i, lm2)
        END IF
      END DO
    END DO
    IF (conserv) THEN
! Compute vertical integral of Output TE
! --------------------------------------
      vsum2(:) = r0
      vsum2_tl = 0.0
      DO i=i1,i2
        DO k=1,kn
          vsum2_tl(i) = vsum2_tl(i) + q2_tl(i, j, k)*(pe2(i, k+1)-pe2(i&
&           , k)) + q2(i, j, k)*(pe2_tl(i, k+1)-pe2_tl(i, k))
          vsum2(i) = vsum2(i) + q2(i, j, k)*(pe2(i, k+1)-pe2(i, k))
        END DO
        vsum2_tl(i) = (vsum2_tl(i)*(pe2(i, kn+1)-pe2(i, 1))-vsum2(i)*(&
&         pe2_tl(i, kn+1)-pe2_tl(i, 1)))/(pe2(i, kn+1)-pe2(i, 1))**2
        vsum2(i) = vsum2(i)/(pe2(i, kn+1)-pe2(i, 1))
      END DO
! Adjust Final TE to conserve
! ---------------------------
      DO i=i1,i2
        DO k=1,kn
          q2_tl(i, j, k) = q2_tl(i, j, k) + vsum1_tl(i) - vsum2_tl(i)
          q2(i, j, k) = q2(i, j, k) + vsum1(i) - vsum2(i)
        END DO
      END DO
    END IF
!          q2(i,j,k) = q2(i,j,k) * vsum1(i)/vsum2(i)
    RETURN
  END SUBROUTINE MAP1_CUBIC_TLM
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  map1_cubic --- Cubic Interpolation for vertical re-mapping
!
! !INTERFACE:
  SUBROUTINE MAP1_CUBIC(km, pe1, kn, pe2, q2, i1, i2, j, ibeg, iend, &
&   jbeg, jend, akap, t_var, conserv)
    IMPLICIT NONE
!EOC
! !INPUT PARAMETERS:
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
    REAL, INTENT(IN) :: akap
! Thermodynamic variable to remap
    INTEGER, INTENT(IN) :: t_var
!     1:TE  2:T  3:PT
    LOGICAL, INTENT(IN) :: conserv
! Current latitude
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
!      real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
! !DESCRIPTION:
!
!     Perform Cubic Interpolation a given latitude
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY:
!    2005.11.14   Takacs    Initial Code
!    2016.07.20   Putman    Modified to make genaric for any thermodynamic variable
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    REAL :: qx(i1:i2, km)
    REAL :: logpl1(i1:i2, km)
    REAL :: logpl2(i1:i2, kn)
    REAL :: dlogp1(i1:i2, km)
    REAL :: vsum1(i1:i2)
    REAL :: vsum2(i1:i2)
    REAL :: am2, am1, ap0, ap1, p, plp1, plp0, plm1, plm2, dlp0, dlm1, &
&   dlm2
    INTEGER :: i, k, lm2, lm1, lp0, lp1
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    INTRINSIC MIN
    REAL, DIMENSION(i2-i1+1) :: arg1
    REAL, DIMENSION(i1:i2) :: arg2
! Initialization
! --------------
    SELECT CASE  (t_var) 
    CASE (1) 
! Total Energy Remapping in Log(P)
      DO k=1,km
        qx(:, k) = q2(i1:i2, j, k)
        logpl1(:, k) = LOG(r2*(pe1(:, k)+pe1(:, k+1)))
      END DO
      DO k=1,kn
        logpl2(:, k) = LOG(r2*(pe2(:, k)+pe2(:, k+1)))
      END DO
      DO k=1,km-1
        dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
      END DO
    CASE (2) 
! Temperature Remapping in Log(P)
      DO k=1,km
        qx(:, k) = q2(i1:i2, j, k)
        logpl1(:, k) = LOG(r2*(pe1(:, k)+pe1(:, k+1)))
      END DO
      DO k=1,kn
        logpl2(:, k) = LOG(r2*(pe2(:, k)+pe2(:, k+1)))
      END DO
      DO k=1,km-1
        dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
      END DO
    CASE (3) 
! Potential Temperature Remapping in P^KAPPA
      DO k=1,km
        qx(:, k) = q2(i1:i2, j, k)
        arg1(:) = r2*(pe1(:, k)+pe1(:, k+1))
        arg2(:) = akap*LOG(arg1(:))
        logpl1(:, k) = EXP(arg2(:))
      END DO
      DO k=1,kn
        arg1(:) = r2*(pe2(:, k)+pe2(:, k+1))
        arg2(:) = akap*LOG(arg1(:))
        logpl2(:, k) = EXP(arg2(:))
      END DO
      DO k=1,km-1
        dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
      END DO
    END SELECT
    IF (conserv) THEN
! Compute vertical integral of Input TE
! -------------------------------------
      vsum1(:) = r0
      DO i=i1,i2
        DO k=1,km
          vsum1(i) = vsum1(i) + qx(i, k)*(pe1(i, k+1)-pe1(i, k))
        END DO
        vsum1(i) = vsum1(i)/(pe1(i, km+1)-pe1(i, 1))
      END DO
    END IF
! Interpolate TE onto target Pressures
! ------------------------------------
    DO i=i1,i2
      DO k=1,kn
        lm1 = 1
        lp0 = 1
        DO WHILE (lp0 .LE. km)
          IF (logpl1(i, lp0) .LT. logpl2(i, k)) THEN
            lp0 = lp0 + 1
          ELSE
            GOTO 100
          END IF
        END DO
 100    IF (lp0 - 1 .LT. 1) THEN
          lm1 = 1
        ELSE
          lm1 = lp0 - 1
        END IF
        IF (lp0 .GT. km) THEN
          lp0 = km
        ELSE
          lp0 = lp0
        END IF
! Extrapolate Linearly in LogP above first model level
! ----------------------------------------------------
        IF (lm1 .EQ. 1 .AND. lp0 .EQ. 1) THEN
          q2(i, j, k) = qx(i, 1) + (qx(i, 2)-qx(i, 1))*(logpl2(i, k)-&
&           logpl1(i, 1))/(logpl1(i, 2)-logpl1(i, 1))
! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
        ELSE IF (lm1 .EQ. km .AND. lp0 .EQ. km) THEN
          q2(i, j, k) = qx(i, km) + (qx(i, km)-qx(i, km-1))*(logpl2(i, k&
&           )-logpl1(i, km))/(logpl1(i, km)-logpl1(i, km-1))
! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
        ELSE IF (lm1 .EQ. 1 .OR. lp0 .EQ. km) THEN
          q2(i, j, k) = qx(i, lp0) + (qx(i, lm1)-qx(i, lp0))*(logpl2(i, &
&           k)-logpl1(i, lp0))/(logpl1(i, lm1)-logpl1(i, lp0))
! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
        ELSE
          lp1 = lp0 + 1
          lm2 = lm1 - 1
          p = logpl2(i, k)
          plp1 = logpl1(i, lp1)
          plp0 = logpl1(i, lp0)
          plm1 = logpl1(i, lm1)
          plm2 = logpl1(i, lm2)
          dlp0 = dlogp1(i, lp0)
          dlm1 = dlogp1(i, lm1)
          dlm2 = dlogp1(i, lm2)
          ap1 = (p-plp0)*(p-plm1)*(p-plm2)/(dlp0*(dlp0+dlm1)*(dlp0+dlm1+&
&           dlm2))
          ap0 = (plp1-p)*(p-plm1)*(p-plm2)/(dlp0*dlm1*(dlm1+dlm2))
          am1 = (plp1-p)*(plp0-p)*(p-plm2)/(dlm1*dlm2*(dlp0+dlm1))
          am2 = (plp1-p)*(plp0-p)*(plm1-p)/(dlm2*(dlm1+dlm2)*(dlp0+dlm1+&
&           dlm2))
          q2(i, j, k) = ap1*qx(i, lp1) + ap0*qx(i, lp0) + am1*qx(i, lm1)&
&           + am2*qx(i, lm2)
        END IF
      END DO
    END DO
    IF (conserv) THEN
! Compute vertical integral of Output TE
! --------------------------------------
      vsum2(:) = r0
      DO i=i1,i2
        DO k=1,kn
          vsum2(i) = vsum2(i) + q2(i, j, k)*(pe2(i, k+1)-pe2(i, k))
        END DO
        vsum2(i) = vsum2(i)/(pe2(i, kn+1)-pe2(i, 1))
      END DO
! Adjust Final TE to conserve
! ---------------------------
      DO i=i1,i2
        DO k=1,kn
          q2(i, j, k) = q2(i, j, k) + vsum1(i) - vsum2(i)
        END DO
      END DO
    END IF
!          q2(i,j,k) = q2(i,j,k) * vsum1(i)/vsum2(i)
    RETURN
  END SUBROUTINE MAP1_CUBIC
!  Differentiation of map_scalar in forward (tangent) mode:
!   variations   of useful results: q2
!   with respect to varying inputs: pe1 pe2 q2
!-----------------------------------------------------------------------
  SUBROUTINE MAP_SCALAR_TLM(km, pe1, pe1_tl, qs, kn, pe2, pe2_tl, q2&
&   , q2_tl, i1, i2, j, ibeg, iend, jbeg, jend, iv, kord, q_min)
    IMPLICIT NONE
! iv=1
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Mode: 0 == constituents  1 == temp
    INTEGER, INTENT(IN) :: iv
!       2 == remap temp with cs scheme
! Method order
    INTEGER, INTENT(IN) :: kord
! Current latitude
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! bottom BC
    REAL, INTENT(IN) :: qs(i1:i2)
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL, INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL, INTENT(IN) :: pe2_tl(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(INOUT) :: q2_tl(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(IN) :: q_min
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_tl(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: q4_tl(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    REAL :: pl_tl, pr_tl, qsum_tl, dp_tl, esl_tl
    INTEGER :: i, k, l, m, k0
    dp1_tl = 0.0
    q4_tl = 0.0
    DO k=1,km
      DO i=i1,i2
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4_tl(1, i, k) = q2_tl(i, j, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL SCALAR_PROFILE_TLM(qs, q4, q4_tl, dp1, dp1_tl, km, i1, i2&
&                          , iv, kord, q_min)
!else
!     call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
      qsum_tl = 0.0
    ELSE
      qsum_tl = 0.0
    END IF
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) GOTO 100
        END DO
        GOTO 123
 100    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2_tl(i, j, k) = q4_tl(2, i, l) + 0.5*((q4_tl(4, i, l)+q4_tl(3&
&           , i, l)-q4_tl(2, i, l))*(pr+pl)+(q4(4, i, l)+q4(3, i, l)-q4(&
&           2, i, l))*(pr_tl+pl_tl)) - r3*(q4_tl(4, i, l)*(pr*(pr+pl)+pl&
&           **2)+q4(4, i, l)*(pr_tl*(pr+pl)+pr*(pr_tl+pl_tl)+2*pl*pl_tl)&
&           )
          q2(i, j, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2&
&           , i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          qsum_tl = (pe1_tl(i, l+1)-pe2_tl(i, k))*(q4(2, i, l)+0.5*(q4(4&
&           , i, l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.&
&           +pl*(1.+pl)))) + (pe1(i, l+1)-pe2(i, k))*(q4_tl(2, i, l)+0.5&
&           *((q4_tl(4, i, l)+q4_tl(3, i, l)-q4_tl(2, i, l))*(1.+pl)+(q4&
&           (4, i, l)+q4(3, i, l)-q4(2, i, l))*pl_tl)-r3*(q4_tl(4, i, l)&
&           *(1.+pl*(1.+pl))+q4(4, i, l)*(pl_tl*(1.+pl)+pl*pl_tl)))
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
              qsum_tl = qsum_tl + dp1_tl(i, m)*q4(1, i, m) + dp1(i, m)*&
&               q4_tl(1, i, m)
              qsum = qsum + dp1(i, m)*q4(1, i, m)
            ELSE
              GOTO 110
            END IF
          END DO
          GOTO 123
 110      dp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          dp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (dp_tl*dp1(i, m)-dp*dp1_tl(i, m))/dp1(i, m)**2
          esl = dp/dp1(i, m)
          qsum_tl = qsum_tl + dp_tl*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4&
&           (2, i, m)+q4(4, i, m)*(1.-r23*esl))) + dp*(q4_tl(2, i, m)+&
&           0.5*(esl_tl*(q4(3, i, m)-q4(2, i, m)+q4(4, i, m)*(1.-r23*esl&
&           ))+esl*(q4_tl(3, i, m)-q4_tl(2, i, m)+q4_tl(4, i, m)*(1.-r23&
&           *esl)-q4(4, i, m)*r23*esl_tl)))
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
        END IF
 123    q2_tl(i, j, k) = (qsum_tl*(pe2(i, k+1)-pe2(i, k))-qsum*(pe2_tl(i&
&         , k+1)-pe2_tl(i, k)))/(pe2(i, k+1)-pe2(i, k))**2
        q2(i, j, k) = qsum/(pe2(i, k+1)-pe2(i, k))
 555  CONTINUE
    END DO
  END SUBROUTINE MAP_SCALAR_TLM
!-----------------------------------------------------------------------
!  Differentiation of map1_ppm in forward (tangent) mode:
!   variations   of useful results: q2
!   with respect to varying inputs: pe1 pe2 qs q2
  SUBROUTINE MAP1_PPM_TLM(km, pe1, pe1_tl, qs, qs_tl, kn, pe2, pe2_tl&
&   , q2, q2_tl, i1, i2, j, ibeg, iend, jbeg, jend, iv, kord)
    IMPLICIT NONE
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Mode: 0 == constituents  1 == ???
    INTEGER, INTENT(IN) :: iv
!       2 == remap temp with cs scheme
! Method order
    INTEGER, INTENT(IN) :: kord
! Current latitude
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! bottom BC
    REAL, INTENT(IN) :: qs(i1:i2)
    REAL, INTENT(IN) :: qs_tl(i1:i2)
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL, INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL, INTENT(IN) :: pe2_tl(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(INOUT) :: q2_tl(ibeg:iend, jbeg:jend, kn)
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_tl(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: q4_tl(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    REAL :: pl_tl, pr_tl, qsum_tl, dp_tl, esl_tl
    INTEGER :: i, k, l, m, k0
    dp1_tl = 0.0
    q4_tl = 0.0
    DO k=1,km
      DO i=i1,i2
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4_tl(1, i, k) = q2_tl(i, j, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE_TLM(qs, qs_tl, q4, q4_tl, dp1, dp1_tl, km, i1, &
&                      i2, iv, kord)
!else
!     call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
      qsum_tl = 0.0
    ELSE
      qsum_tl = 0.0
    END IF
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) GOTO 100
        END DO
        GOTO 123
 100    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2_tl(i, j, k) = q4_tl(2, i, l) + 0.5*((q4_tl(4, i, l)+q4_tl(3&
&           , i, l)-q4_tl(2, i, l))*(pr+pl)+(q4(4, i, l)+q4(3, i, l)-q4(&
&           2, i, l))*(pr_tl+pl_tl)) - r3*(q4_tl(4, i, l)*(pr*(pr+pl)+pl&
&           **2)+q4(4, i, l)*(pr_tl*(pr+pl)+pr*(pr_tl+pl_tl)+2*pl*pl_tl)&
&           )
          q2(i, j, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2&
&           , i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          qsum_tl = (pe1_tl(i, l+1)-pe2_tl(i, k))*(q4(2, i, l)+0.5*(q4(4&
&           , i, l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.&
&           +pl*(1.+pl)))) + (pe1(i, l+1)-pe2(i, k))*(q4_tl(2, i, l)+0.5&
&           *((q4_tl(4, i, l)+q4_tl(3, i, l)-q4_tl(2, i, l))*(1.+pl)+(q4&
&           (4, i, l)+q4(3, i, l)-q4(2, i, l))*pl_tl)-r3*(q4_tl(4, i, l)&
&           *(1.+pl*(1.+pl))+q4(4, i, l)*(pl_tl*(1.+pl)+pl*pl_tl)))
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
              qsum_tl = qsum_tl + dp1_tl(i, m)*q4(1, i, m) + dp1(i, m)*&
&               q4_tl(1, i, m)
              qsum = qsum + dp1(i, m)*q4(1, i, m)
            ELSE
              GOTO 110
            END IF
          END DO
          GOTO 123
 110      dp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          dp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (dp_tl*dp1(i, m)-dp*dp1_tl(i, m))/dp1(i, m)**2
          esl = dp/dp1(i, m)
          qsum_tl = qsum_tl + dp_tl*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4&
&           (2, i, m)+q4(4, i, m)*(1.-r23*esl))) + dp*(q4_tl(2, i, m)+&
&           0.5*(esl_tl*(q4(3, i, m)-q4(2, i, m)+q4(4, i, m)*(1.-r23*esl&
&           ))+esl*(q4_tl(3, i, m)-q4_tl(2, i, m)+q4_tl(4, i, m)*(1.-r23&
&           *esl)-q4(4, i, m)*r23*esl_tl)))
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
        END IF
 123    q2_tl(i, j, k) = (qsum_tl*(pe2(i, k+1)-pe2(i, k))-qsum*(pe2_tl(i&
&         , k+1)-pe2_tl(i, k)))/(pe2(i, k+1)-pe2(i, k))**2
        q2(i, j, k) = qsum/(pe2(i, k+1)-pe2(i, k))
 555  CONTINUE
    END DO
  END SUBROUTINE MAP1_PPM_TLM
!  Differentiation of mapn_tracer in forward (tangent) mode:
!   variations   of useful results: q1
!   with respect to varying inputs: pe1 pe2 dp2 q1
  SUBROUTINE MAPN_TRACER_TLM(nq, km, pe1, pe1_tl, pe2, pe2_tl, q1, &
&   q1_tl, dp2, dp2_tl, kord, j, i1, i2, isd, ied, jsd, jed, q_min, fill&
& )
    IMPLICIT NONE
! !INPUT PARAMETERS:
! vertical dimension
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: j, nq, i1, i2
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: kord(nq)
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL, INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, km+1)
    REAL, INTENT(IN) :: pe2_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
    REAL, INTENT(IN) :: dp2(i1:i2, km)
    REAL, INTENT(IN) :: dp2_tl(i1:i2, km)
    REAL, INTENT(IN) :: q_min
    LOGICAL, INTENT(IN) :: fill
! Field input
    REAL, INTENT(INOUT) :: q1(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: q1_tl(isd:ied, jsd:jed, km, nq)
! !LOCAL VARIABLES:
    REAL :: q4(4, i1:i2, km, nq)
    REAL :: q4_tl(4, i1:i2, km, nq)
! Field output
    REAL :: q2(i1:i2, km, nq)
    REAL :: q2_tl(i1:i2, km, nq)
    REAL :: qsum(nq)
    REAL :: qsum_tl(nq)
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_tl(i1:i2, km)
    REAL :: qs(i1:i2)
    REAL :: pl, pr, dp, esl, fac1, fac2
    REAL :: pl_tl, pr_tl, dp_tl, esl_tl, fac1_tl, fac2_tl
    INTEGER :: i, k, l, m, k0, iq
    dp1_tl = 0.0
    DO k=1,km
      DO i=i1,i2
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
      END DO
    END DO
    q4_tl = 0.0
    DO iq=1,nq
      DO k=1,km
        DO i=i1,i2
          q4_tl(1, i, k, iq) = q1_tl(i, j, k, iq)
          q4(1, i, k, iq) = q1(i, j, k, iq)
        END DO
      END DO
      CALL SCALAR_PROFILE_TLM(qs, q4(1:4, i1:i2, 1:km, iq), q4_tl(1:4&
&                          , i1:i2, 1:km, iq), dp1, dp1_tl, km, i1, i2, &
&                          0, kord(iq), q_min)
    END DO
    qsum_tl = 0.0
    q2_tl = 0.0
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,km
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) GOTO 110
        END DO
        GOTO 123
 110    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          fac1_tl = pr_tl + pl_tl
          fac1 = pr + pl
          fac2_tl = r3*(pr_tl*fac1+pr*fac1_tl+pl_tl*pl+pl*pl_tl)
          fac2 = r3*(pr*fac1+pl*pl)
          fac1_tl = 0.5*fac1_tl
          fac1 = 0.5*fac1
          DO iq=1,nq
            q2_tl(i, k, iq) = q4_tl(2, i, l, iq) + (q4_tl(4, i, l, iq)+&
&             q4_tl(3, i, l, iq)-q4_tl(2, i, l, iq))*fac1 + (q4(4, i, l&
&             , iq)+q4(3, i, l, iq)-q4(2, i, l, iq))*fac1_tl - q4_tl(4, &
&             i, l, iq)*fac2 - q4(4, i, l, iq)*fac2_tl
            q2(i, k, iq) = q4(2, i, l, iq) + (q4(4, i, l, iq)+q4(3, i, l&
&             , iq)-q4(2, i, l, iq))*fac1 - q4(4, i, l, iq)*fac2
          END DO
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          dp_tl = pe1_tl(i, l+1) - pe2_tl(i, k)
          dp = pe1(i, l+1) - pe2(i, k)
          fac1_tl = pl_tl
          fac1 = 1. + pl
          fac2_tl = r3*(pl_tl*fac1+pl*fac1_tl)
          fac2 = r3*(1.+pl*fac1)
          fac1_tl = 0.5*fac1_tl
          fac1 = 0.5*fac1
          DO iq=1,nq
            qsum_tl(iq) = dp_tl*(q4(2, i, l, iq)+(q4(4, i, l, iq)+q4(3, &
&             i, l, iq)-q4(2, i, l, iq))*fac1-q4(4, i, l, iq)*fac2) + dp&
&             *(q4_tl(2, i, l, iq)+(q4_tl(4, i, l, iq)+q4_tl(3, i, l, iq&
&             )-q4_tl(2, i, l, iq))*fac1+(q4(4, i, l, iq)+q4(3, i, l, iq&
&             )-q4(2, i, l, iq))*fac1_tl-q4_tl(4, i, l, iq)*fac2-q4(4, i&
&             , l, iq)*fac2_tl)
            qsum(iq) = dp*(q4(2, i, l, iq)+(q4(4, i, l, iq)+q4(3, i, l, &
&             iq)-q4(2, i, l, iq))*fac1-q4(4, i, l, iq)*fac2)
          END DO
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
              DO iq=1,nq
                qsum_tl(iq) = qsum_tl(iq) + dp1_tl(i, m)*q4(1, i, m, iq)&
&                 + dp1(i, m)*q4_tl(1, i, m, iq)
                qsum(iq) = qsum(iq) + dp1(i, m)*q4(1, i, m, iq)
              END DO
            ELSE
              GOTO 120
            END IF
          END DO
          GOTO 123
 120      dp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          dp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (dp_tl*dp1(i, m)-dp*dp1_tl(i, m))/dp1(i, m)**2
          esl = dp/dp1(i, m)
          fac1_tl = 0.5*esl_tl
          fac1 = 0.5*esl
          fac2_tl = -(r23*esl_tl)
          fac2 = 1. - r23*esl
          DO iq=1,nq
            qsum_tl(iq) = qsum_tl(iq) + dp_tl*(q4(2, i, m, iq)+fac1*(q4(&
&             3, i, m, iq)-q4(2, i, m, iq)+q4(4, i, m, iq)*fac2)) + dp*(&
&             q4_tl(2, i, m, iq)+fac1_tl*(q4(3, i, m, iq)-q4(2, i, m, iq&
&             )+q4(4, i, m, iq)*fac2)+fac1*(q4_tl(3, i, m, iq)-q4_tl(2, &
&             i, m, iq)+q4_tl(4, i, m, iq)*fac2+q4(4, i, m, iq)*fac2_tl)&
&             )
            qsum(iq) = qsum(iq) + dp*(q4(2, i, m, iq)+fac1*(q4(3, i, m, &
&             iq)-q4(2, i, m, iq)+q4(4, i, m, iq)*fac2))
          END DO
          k0 = m
        END IF
 123    DO iq=1,nq
          q2_tl(i, k, iq) = (qsum_tl(iq)*dp2(i, k)-qsum(iq)*dp2_tl(i, k)&
&           )/dp2(i, k)**2
          q2(i, k, iq) = qsum(iq)/dp2(i, k)
        END DO
 555  CONTINUE
    END DO
    IF (fill) CALL FILLZ(i2 - i1 + 1, km, nq, q2, dp2)
    DO iq=1,nq
!    if (fill) call fillz(i2-i1+1, km, 1, q2(i1,1,iq), dp2)
      DO k=1,km
        DO i=i1,i2
          q1_tl(i, j, k, iq) = q2_tl(i, k, iq)
          q1(i, j, k, iq) = q2(i, k, iq)
        END DO
      END DO
    END DO
  END SUBROUTINE MAPN_TRACER_TLM
!  Differentiation of map1_q2 in forward (tangent) mode:
!   variations   of useful results: q2
!   with respect to varying inputs: pe1 pe2 dp2 q1 q2
  SUBROUTINE MAP1_Q2_TLM(km, pe1, pe1_tl, q1, q1_tl, kn, pe2, pe2_tl&
&   , q2, q2_tl, dp2, dp2_tl, i1, i2, iv, kord, j, ibeg, iend, jbeg, &
&   jend, q_min)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: i1, i2
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Mode: 0 ==  constituents 1 == ???
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL, INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL, INTENT(IN) :: pe2_tl(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! Field input
    REAL, INTENT(IN) :: q1(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(IN) :: q1_tl(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(IN) :: dp2(i1:i2, kn)
    REAL, INTENT(IN) :: dp2_tl(i1:i2, kn)
    REAL, INTENT(IN) :: q_min
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(i1:i2, kn)
    REAL, INTENT(INOUT) :: q2_tl(i1:i2, kn)
! !LOCAL VARIABLES:
    REAL :: qs(i1:i2)
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_tl(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: q4_tl(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    REAL :: pl_tl, pr_tl, qsum_tl, dp_tl, esl_tl
    INTEGER :: i, k, l, m, k0
    dp1_tl = 0.0
    q4_tl = 0.0
    DO k=1,km
      DO i=i1,i2
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4_tl(1, i, k) = q1_tl(i, j, k)
        q4(1, i, k) = q1(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL SCALAR_PROFILE_TLM(qs, q4, q4_tl, dp1, dp1_tl, km, i1, i2&
&                          , iv, kord, q_min)
!else
!call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
      qsum_tl = 0.0
    ELSE
      qsum_tl = 0.0
    END IF
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) GOTO 110
        END DO
        GOTO 123
 110    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2_tl(i, k) = q4_tl(2, i, l) + 0.5*((q4_tl(4, i, l)+q4_tl(3, i&
&           , l)-q4_tl(2, i, l))*(pr+pl)+(q4(4, i, l)+q4(3, i, l)-q4(2, &
&           i, l))*(pr_tl+pl_tl)) - r3*(q4_tl(4, i, l)*(pr*(pr+pl)+pl**2&
&           )+q4(4, i, l)*(pr_tl*(pr+pl)+pr*(pr_tl+pl_tl)+2*pl*pl_tl))
          q2(i, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i&
&           , l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          qsum_tl = (pe1_tl(i, l+1)-pe2_tl(i, k))*(q4(2, i, l)+0.5*(q4(4&
&           , i, l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.&
&           +pl*(1.+pl)))) + (pe1(i, l+1)-pe2(i, k))*(q4_tl(2, i, l)+0.5&
&           *((q4_tl(4, i, l)+q4_tl(3, i, l)-q4_tl(2, i, l))*(1.+pl)+(q4&
&           (4, i, l)+q4(3, i, l)-q4(2, i, l))*pl_tl)-r3*(q4_tl(4, i, l)&
&           *(1.+pl*(1.+pl))+q4(4, i, l)*(pl_tl*(1.+pl)+pl*pl_tl)))
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
              qsum_tl = qsum_tl + dp1_tl(i, m)*q4(1, i, m) + dp1(i, m)*&
&               q4_tl(1, i, m)
              qsum = qsum + dp1(i, m)*q4(1, i, m)
            ELSE
              GOTO 120
            END IF
          END DO
          GOTO 123
 120      dp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          dp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (dp_tl*dp1(i, m)-dp*dp1_tl(i, m))/dp1(i, m)**2
          esl = dp/dp1(i, m)
          qsum_tl = qsum_tl + dp_tl*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4&
&           (2, i, m)+q4(4, i, m)*(1.-r23*esl))) + dp*(q4_tl(2, i, m)+&
&           0.5*(esl_tl*(q4(3, i, m)-q4(2, i, m)+q4(4, i, m)*(1.-r23*esl&
&           ))+esl*(q4_tl(3, i, m)-q4_tl(2, i, m)+q4_tl(4, i, m)*(1.-r23&
&           *esl)-q4(4, i, m)*r23*esl_tl)))
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
        END IF
 123    q2_tl(i, k) = (qsum_tl*dp2(i, k)-qsum*dp2_tl(i, k))/dp2(i, k)**2
        q2(i, k) = qsum/dp2(i, k)
 555  CONTINUE
    END DO
  END SUBROUTINE MAP1_Q2_TLM
!  Differentiation of scalar_profile in forward (tangent) mode:
!   variations   of useful results: a4
!   with respect to varying inputs: delp a4
  SUBROUTINE SCALAR_PROFILE_TLM(qs, a4, a4_tl, delp, delp_tl, km, i1&
&   , i2, iv, kord, qmin)
    IMPLICIT NONE
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
    INTEGER, INTENT(IN) :: i1, i2
! vertical dimension
    INTEGER, INTENT(IN) :: km
! iv =-1: winds
    INTEGER, INTENT(IN) :: iv
! iv = 0: positive definite scalars
! iv = 1: others
    INTEGER, INTENT(IN) :: kord
    REAL, INTENT(IN) :: qs(i1:i2)
! layer pressure thickness
    REAL, INTENT(IN) :: delp(i1:i2, km)
    REAL, INTENT(IN) :: delp_tl(i1:i2, km)
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_tl(4, i1:i2, km)
    REAL, INTENT(IN) :: qmin
!-----------------------------------------------------------------------
    LOGICAL, DIMENSION(i1:i2, km) :: extm, ext6
    REAL :: gam(i1:i2, km)
    REAL :: gam_tl(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: q_tl(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: d4_tl(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: bet_tl, a_bot_tl, grat_tl
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    INTEGER :: i, k, im
    INTRINSIC ABS
    INTEGER :: abs0
    IF (iv .EQ. -2) THEN
      q_tl = 0.0
      DO i=i1,i2
        gam(i, 2) = 0.5
        q_tl(i, 1) = 1.5*a4_tl(1, i, 1)
        q(i, 1) = 1.5*a4(1, i, 1)
      END DO
      gam_tl = 0.0
      DO k=2,km-1
        DO i=i1,i2
          grat_tl = (delp_tl(i, k-1)*delp(i, k)-delp(i, k-1)*delp_tl(i, &
&           k))/delp(i, k)**2
          grat = delp(i, k-1)/delp(i, k)
          bet_tl = 2*grat_tl - gam_tl(i, k)
          bet = 2. + grat + grat - gam(i, k)
          q_tl(i, k) = ((3.*(a4_tl(1, i, k-1)+a4_tl(1, i, k))-q_tl(i, k-&
&           1))*bet-(3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))*bet_tl)/&
&           bet**2
          q(i, k) = (3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))/bet
          gam_tl(i, k+1) = (grat_tl*bet-grat*bet_tl)/bet**2
          gam(i, k+1) = grat/bet
        END DO
      END DO
      DO i=i1,i2
        grat_tl = (delp_tl(i, km-1)*delp(i, km)-delp(i, km-1)*delp_tl(i&
&         , km))/delp(i, km)**2
        grat = delp(i, km-1)/delp(i, km)
        q_tl(i, km) = ((3.*(a4_tl(1, i, km-1)+a4_tl(1, i, km))-qs(i)*&
&         grat_tl-q_tl(i, km-1))*(2.+grat+grat-gam(i, km))-(3.*(a4(1, i&
&         , km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-1))*(2*grat_tl-gam_tl&
&         (i, km)))/(2.+grat+grat-gam(i, km))**2
        q(i, km) = (3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-&
&         1))/(2.+grat+grat-gam(i, km))
        q_tl(i, km+1) = 0.0
        q(i, km+1) = qs(i)
      END DO
      DO k=km-1,1,-1
        DO i=i1,i2
          q_tl(i, k) = q_tl(i, k) - gam_tl(i, k+1)*q(i, k+1) - gam(i, k+&
&           1)*q_tl(i, k+1)
          q(i, k) = q(i, k) - gam(i, k+1)*q(i, k+1)
        END DO
      END DO
    ELSE
      q_tl = 0.0
      gam_tl = 0.0
      DO i=i1,i2
! grid ratio
        grat_tl = (delp_tl(i, 2)*delp(i, 1)-delp(i, 2)*delp_tl(i, 1))/&
&         delp(i, 1)**2
        grat = delp(i, 2)/delp(i, 1)
        bet_tl = grat_tl*(grat+0.5) + grat*grat_tl
        bet = grat*(grat+0.5)
        q_tl(i, 1) = (((2*grat_tl*(grat+1.)+(grat+grat)*grat_tl)*a4(1, i&
&         , 1)+(grat+grat)*(grat+1.)*a4_tl(1, i, 1)+a4_tl(1, i, 2))*bet-&
&         ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))*bet_tl)/bet**2
        q(i, 1) = ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))/bet
        gam_tl(i, 1) = ((grat_tl*(grat+1.5)+grat*grat_tl)*bet-(1.+grat*(&
&         grat+1.5))*bet_tl)/bet**2
        gam(i, 1) = (1.+grat*(grat+1.5))/bet
      END DO
      d4_tl = 0.0
      DO k=2,km
        DO i=i1,i2
          d4_tl(i) = (delp_tl(i, k-1)*delp(i, k)-delp(i, k-1)*delp_tl(i&
&           , k))/delp(i, k)**2
          d4(i) = delp(i, k-1)/delp(i, k)
          bet_tl = 2*d4_tl(i) - gam_tl(i, k-1)
          bet = 2. + d4(i) + d4(i) - gam(i, k-1)
          q_tl(i, k) = ((3.*(a4_tl(1, i, k-1)+d4_tl(i)*a4(1, i, k)+d4(i)&
&           *a4_tl(1, i, k))-q_tl(i, k-1))*bet-(3.*(a4(1, i, k-1)+d4(i)*&
&           a4(1, i, k))-q(i, k-1))*bet_tl)/bet**2
          q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
          gam_tl(i, k) = (d4_tl(i)*bet-d4(i)*bet_tl)/bet**2
          gam(i, k) = d4(i)/bet
        END DO
      END DO
      DO i=i1,i2
        a_bot_tl = d4_tl(i)*(d4(i)+1.5) + d4(i)*d4_tl(i)
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        q_tl(i, km+1) = ((2.*((d4_tl(i)*(d4(i)+1.)+d4(i)*d4_tl(i))*a4(1&
&         , i, km)+d4(i)*(d4(i)+1.)*a4_tl(1, i, km))+a4_tl(1, i, km-1)-&
&         a_bot_tl*q(i, km)-a_bot*q_tl(i, km))*(d4(i)*(d4(i)+0.5)-a_bot*&
&         gam(i, km))-(2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))*(d4_tl(i)*(d4(i)+0.5)+d4(i)*d4_tl(i)-a_bot_tl*&
&         gam(i, km)-a_bot*gam_tl(i, km)))/(d4(i)*(d4(i)+0.5)-a_bot*gam(&
&         i, km))**2
        q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          q_tl(i, k) = q_tl(i, k) - gam_tl(i, k)*q(i, k+1) - gam(i, k)*&
&           q_tl(i, k+1)
          q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
        END DO
      END DO
    END IF
    IF (kord .GE. 0.) THEN
      abs0 = kord
    ELSE
      abs0 = -kord
    END IF
!----- Perfectly linear scheme --------------------------------
    IF (abs0 .GT. 16) THEN
      DO k=1,km
        DO i=i1,i2
          a4_tl(2, i, k) = q_tl(i, k)
          a4(2, i, k) = q(i, k)
          a4_tl(3, i, k) = q_tl(i, k+1)
          a4(3, i, k) = q(i, k+1)
          a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3&
&           , i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
      END DO
      RETURN
    END IF
  END SUBROUTINE SCALAR_PROFILE_TLM
!  Differentiation of cs_profile in forward (tangent) mode:
!   variations   of useful results: a4
!   with respect to varying inputs: qs delp a4
  SUBROUTINE CS_PROFILE_TLM(qs, qs_tl, a4, a4_tl, delp, delp_tl, km, &
&   i1, i2, iv, kord)
    IMPLICIT NONE
!----- Perfectly linear scheme --------------------------------
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
    INTEGER, INTENT(IN) :: i1, i2
! vertical dimension
    INTEGER, INTENT(IN) :: km
! iv =-1: winds
    INTEGER, INTENT(IN) :: iv
! iv = 0: positive definite scalars
! iv = 1: others
    INTEGER, INTENT(IN) :: kord
    REAL, INTENT(IN) :: qs(i1:i2)
    REAL, INTENT(IN) :: qs_tl(i1:i2)
! layer pressure thickness
    REAL, INTENT(IN) :: delp(i1:i2, km)
    REAL, INTENT(IN) :: delp_tl(i1:i2, km)
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_tl(4, i1:i2, km)
!-----------------------------------------------------------------------
    LOGICAL :: extm(i1:i2, km)
    REAL :: gam(i1:i2, km)
    REAL :: gam_tl(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: q_tl(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: d4_tl(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: bet_tl, a_bot_tl, grat_tl
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    INTEGER :: i, k, im
    INTRINSIC ABS
    INTEGER :: abs0
    IF (iv .EQ. -2) THEN
      q_tl = 0.0
      DO i=i1,i2
        gam(i, 2) = 0.5
        q_tl(i, 1) = 1.5*a4_tl(1, i, 1)
        q(i, 1) = 1.5*a4(1, i, 1)
      END DO
      gam_tl = 0.0
      DO k=2,km-1
        DO i=i1,i2
          grat_tl = (delp_tl(i, k-1)*delp(i, k)-delp(i, k-1)*delp_tl(i, &
&           k))/delp(i, k)**2
          grat = delp(i, k-1)/delp(i, k)
          bet_tl = 2*grat_tl - gam_tl(i, k)
          bet = 2. + grat + grat - gam(i, k)
          q_tl(i, k) = ((3.*(a4_tl(1, i, k-1)+a4_tl(1, i, k))-q_tl(i, k-&
&           1))*bet-(3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))*bet_tl)/&
&           bet**2
          q(i, k) = (3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))/bet
          gam_tl(i, k+1) = (grat_tl*bet-grat*bet_tl)/bet**2
          gam(i, k+1) = grat/bet
        END DO
      END DO
      DO i=i1,i2
        grat_tl = (delp_tl(i, km-1)*delp(i, km)-delp(i, km-1)*delp_tl(i&
&         , km))/delp(i, km)**2
        grat = delp(i, km-1)/delp(i, km)
        q_tl(i, km) = ((3.*(a4_tl(1, i, km-1)+a4_tl(1, i, km))-grat_tl*&
&         qs(i)-grat*qs_tl(i)-q_tl(i, km-1))*(2.+grat+grat-gam(i, km))-(&
&         3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-1))*(2*&
&         grat_tl-gam_tl(i, km)))/(2.+grat+grat-gam(i, km))**2
        q(i, km) = (3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-&
&         1))/(2.+grat+grat-gam(i, km))
        q_tl(i, km+1) = qs_tl(i)
        q(i, km+1) = qs(i)
      END DO
      DO k=km-1,1,-1
        DO i=i1,i2
          q_tl(i, k) = q_tl(i, k) - gam_tl(i, k+1)*q(i, k+1) - gam(i, k+&
&           1)*q_tl(i, k+1)
          q(i, k) = q(i, k) - gam(i, k+1)*q(i, k+1)
        END DO
      END DO
    ELSE
      q_tl = 0.0
      gam_tl = 0.0
      DO i=i1,i2
! grid ratio
        grat_tl = (delp_tl(i, 2)*delp(i, 1)-delp(i, 2)*delp_tl(i, 1))/&
&         delp(i, 1)**2
        grat = delp(i, 2)/delp(i, 1)
        bet_tl = grat_tl*(grat+0.5) + grat*grat_tl
        bet = grat*(grat+0.5)
        q_tl(i, 1) = (((2*grat_tl*(grat+1.)+(grat+grat)*grat_tl)*a4(1, i&
&         , 1)+(grat+grat)*(grat+1.)*a4_tl(1, i, 1)+a4_tl(1, i, 2))*bet-&
&         ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))*bet_tl)/bet**2
        q(i, 1) = ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))/bet
        gam_tl(i, 1) = ((grat_tl*(grat+1.5)+grat*grat_tl)*bet-(1.+grat*(&
&         grat+1.5))*bet_tl)/bet**2
        gam(i, 1) = (1.+grat*(grat+1.5))/bet
      END DO
      d4_tl = 0.0
      DO k=2,km
        DO i=i1,i2
          d4_tl(i) = (delp_tl(i, k-1)*delp(i, k)-delp(i, k-1)*delp_tl(i&
&           , k))/delp(i, k)**2
          d4(i) = delp(i, k-1)/delp(i, k)
          bet_tl = 2*d4_tl(i) - gam_tl(i, k-1)
          bet = 2. + d4(i) + d4(i) - gam(i, k-1)
          q_tl(i, k) = ((3.*(a4_tl(1, i, k-1)+d4_tl(i)*a4(1, i, k)+d4(i)&
&           *a4_tl(1, i, k))-q_tl(i, k-1))*bet-(3.*(a4(1, i, k-1)+d4(i)*&
&           a4(1, i, k))-q(i, k-1))*bet_tl)/bet**2
          q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
          gam_tl(i, k) = (d4_tl(i)*bet-d4(i)*bet_tl)/bet**2
          gam(i, k) = d4(i)/bet
        END DO
      END DO
      DO i=i1,i2
        a_bot_tl = d4_tl(i)*(d4(i)+1.5) + d4(i)*d4_tl(i)
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        q_tl(i, km+1) = ((2.*((d4_tl(i)*(d4(i)+1.)+d4(i)*d4_tl(i))*a4(1&
&         , i, km)+d4(i)*(d4(i)+1.)*a4_tl(1, i, km))+a4_tl(1, i, km-1)-&
&         a_bot_tl*q(i, km)-a_bot*q_tl(i, km))*(d4(i)*(d4(i)+0.5)-a_bot*&
&         gam(i, km))-(2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))*(d4_tl(i)*(d4(i)+0.5)+d4(i)*d4_tl(i)-a_bot_tl*&
&         gam(i, km)-a_bot*gam_tl(i, km)))/(d4(i)*(d4(i)+0.5)-a_bot*gam(&
&         i, km))**2
        q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          q_tl(i, k) = q_tl(i, k) - gam_tl(i, k)*q(i, k+1) - gam(i, k)*&
&           q_tl(i, k+1)
          q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
        END DO
      END DO
    END IF
    IF (kord .GE. 0.) THEN
      abs0 = kord
    ELSE
      abs0 = -kord
    END IF
!----- Perfectly linear scheme --------------------------------
    IF (abs0 .GT. 16) THEN
      DO k=1,km
        DO i=i1,i2
          a4_tl(2, i, k) = q_tl(i, k)
          a4(2, i, k) = q(i, k)
          a4_tl(3, i, k) = q_tl(i, k+1)
          a4(3, i, k) = q(i, k+1)
          a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3&
&           , i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
      END DO
      RETURN
    END IF
  END SUBROUTINE CS_PROFILE_TLM
end module fv_mapz_tlm_mod
