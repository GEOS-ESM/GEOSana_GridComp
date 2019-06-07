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
module fv_mapz_adm_mod

  use constants_mod,     only: radius, pi=>pi_8, rvgas, rdgas, grav, hlv, hlf, cp_air, cp_vapor
  use tracer_manager_mod,only: get_tracer_index
  use field_manager_mod, only: MODEL_ATMOS
  use fv_grid_utils_mod, only: ptop_min
  use fv_grid_utils_adm_mod, only: g_sum
  use fv_grid_utils_adm_mod, only: g_sum_adm
  use fv_fill_mod,       only: fillz
  use mpp_domains_mod,   only: mpp_update_domains, domain2d
  use mpp_mod,           only: FATAL, mpp_error, get_unit, mpp_root_pe, mpp_pe
  use fv_arrays_mod,     only: fv_grid_type, fv_flags_type
  use fv_timing_mod,     only: timing_on, timing_off
  use fv_mp_mod,         only: is_master
  use fv_cmp_mod,        only: qs_init, fv_sat_adj
  use fv_diagnostics_mod, only: prt_mxm

  use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                           pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm
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
  public compute_total_energy_fwd, Lagrangian_to_Eulerian_fwd,  &
         map1_q2_fwd
  public compute_total_energy_bwd, Lagrangian_to_Eulerian_bwd,  &
         map1_q2_bwd

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of lagrangian_to_eulerian in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4_f
!b a2b_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_
!pe dyn_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dy
!n_core_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_
!mod.Rayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils
!_mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv
!_mapz_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mo
!d.remap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.p
!pm_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map
!1_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_m
!od.fv_subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz
!_d nh_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver
! nh_utils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile 
!nh_utils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw
!_core_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_cor
!e_mod.ytp_v sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d
!_fb tp_core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner f
!v_grid_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: ws peln q u v w delp ua delz
!                omga te0_2d pkz pe pk ps pt te
!   with respect to varying inputs: ws peln q u v w delp ua delz
!                omga te0_2d pkz pe pk ps pt te
  SUBROUTINE LAGRANGIAN_TO_EULERIAN_FWD(last_step, consv, ps, pe, delp, &
&   pkz, pk, mdt, pdt, km, is, ie, js, je, isd, ied, jsd, jed, nq, nwat&
&   , sphum, q_con, u, v, w, delz, pt, q, hs, r_vir, cp, akap, cappa, &
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
    REAL :: pkz(is:ie, js:je, km)
    REAL :: te(isd:ied, jsd:jed, km)
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
    INTEGER :: arg1
    REAL :: result1
    LOGICAL :: arg10
    LOGICAL :: res
    INTEGER :: ad_count
    INTEGER :: ad_from
    INTEGER :: ad_from0
    INTEGER :: ad_from1
    INTEGER :: ad_from2
    INTEGER :: ad_from3
    INTEGER :: ad_from4
    INTEGER :: ad_from5
    INTEGER :: ad_from6
    INTEGER :: ad_from7
    INTEGER :: ad_from8
    INTEGER :: ad_from9
    INTEGER :: ad_from10
    INTEGER :: ad_from11
    INTEGER :: ad_from12
    INTEGER :: ad_from13
    INTEGER :: ad_from14
    INTEGER :: ad_from15
    INTEGER :: ad_from16
    INTEGER :: ad_from17
    INTEGER :: ad_from18
    INTEGER :: ad_from19
    INTEGER :: ad_from20
    INTEGER :: ad_from21
    INTEGER :: ad_from22
    INTEGER :: ad_from23
    INTEGER :: ad_from24
    INTEGER :: ad_from25
    INTEGER :: ad_from26
    INTEGER :: ad_count0
    INTEGER :: ad_from27
    INTEGER :: ad_from28
    INTEGER :: ad_from29
    INTEGER :: ad_from30
    INTEGER :: ad_from31
    INTEGER :: ad_from32
    INTEGER :: ad_from33
    INTEGER :: ad_count1

    te_2d = 0.0
    zsum0 = 0.0
    zsum1 = 0.0
    dpln = 0.0
    q2 = 0.0
    dp2 = 0.0
    pe1 = 0.0
    pe2 = 0.0
    pk1 = 0.0
    pk2 = 0.0
    pn2 = 0.0
    phis = 0.0
    pe0 = 0.0
    pe3 = 0.0
    gz = 0.0
    cvm = 0.0
    qv = 0.0
    rcp = 0.0
    rg = 0.0
    tmp = 0.0
    tpe = 0.0
    rrg = 0.0
    bkh = 0.0
    dtmp = 0.0
    k1k = 0.0
    dlnp = 0.0
    result1 = 0.0
    abs0 = 0.0
    nt = 0
    liq_wat = 0
    ice_wat = 0
    rainwat = 0
    snowwat = 0
    cld_amt = 0
    graupel = 0
    iq = 0
    n = 0
    kmp = 0
    kp = 0
    k_next = 0
    abs_kord_tm = 0
    abs_kord_tm_pert = 0
    iep1 = 0
    jep1 = 0
    iedp1 = 0
    jedp1 = 0
    arg1 = 0
    ad_count = 0
    ad_from = 0
    ad_from0 = 0
    ad_from1 = 0
    ad_from2 = 0
    ad_from3 = 0
    ad_from4 = 0
    ad_from5 = 0
    ad_from6 = 0
    ad_from7 = 0
    ad_from8 = 0
    ad_from9 = 0
    ad_from10 = 0
    ad_from11 = 0
    ad_from12 = 0
    ad_from13 = 0
    ad_from14 = 0
    ad_from15 = 0
    ad_from16 = 0
    ad_from17 = 0
    ad_from18 = 0
    ad_from19 = 0
    ad_from20 = 0
    ad_from21 = 0
    ad_from22 = 0
    ad_from23 = 0
    ad_from24 = 0
    ad_from25 = 0
    ad_from26 = 0
    ad_count0 = 0
    ad_from27 = 0
    ad_from28 = 0
    ad_from29 = 0
    ad_from30 = 0
    ad_from31 = 0
    ad_from32 = 0
    ad_from33 = 0
    ad_count1 = 0

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
      STOP
    END SELECT
    res = IS_MASTER()
    IF (res .AND. flagstruct%fv_debug) THEN
      PRINT*, ''
      SELECT CASE  (remap_option) 
      CASE (0) 
        CALL PUSHCONTROL(1,0)
        PRINT*, ' REMAPPING  T in logP '
      CASE (1) 
        CALL PUSHCONTROL(1,0)
        PRINT*, ' REMAPPING PT in P'
      CASE (2) 
        CALL PUSHCONTROL(1,0)
        PRINT*, ' REMAPPING TE in logP with GMAO cubic'
      CASE DEFAULT
        CALL PUSHCONTROL(1,0)
      END SELECT
      PRINT*, ' REMAPPING CONSV:     ', consv
      PRINT*, ' REMAPPING CONSV_MIN: ', consv_min
      PRINT*, ''
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
! akap / (1.-akap) = rg/Cv=0.4
    k1k = rdgas/cv_air
    rg = rdgas
    rrg = -(rdgas/grav)
    IF (do_sat_adj) THEN
      fast_mp_consv = .NOT.do_adiabatic_init .AND. consv .GT. consv_min
      ad_count = 1
      DO k=1,km
        kmp = k
        IF (pfull(k) .GT. 10.e2) THEN
          GOTO 100
        ELSE
          ad_count = ad_count + 1
        END IF
      END DO
      CALL PUSHCONTROL(1,0)
      CALL PUSHINTEGER(ad_count)
      CALL PUSHCONTROL(1,1)
      GOTO 110
 100  CALL PUSHCONTROL(1,1)
      CALL PUSHINTEGER(ad_count)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
 110 CALL PUSHINTEGER(j)
    ad_count1 = 1
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
        ad_from = is
        DO i=ad_from,ie
          CALL PUSHREALARRAY(pe1(i, k))
          pe1(i, k) = pe(i, k, j)
        END DO
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from)
      END DO
      CALL PUSHINTEGER(k - 1)
      ad_from0 = is
      DO i=ad_from0,ie
        CALL PUSHREALARRAY(pe2(i, 1))
        pe2(i, 1) = ptop
        CALL PUSHREALARRAY(pe2(i, km+1))
        pe2(i, km+1) = pe(i, km+1, j)
      END DO
      CALL PUSHINTEGER(i - 1)
      CALL PUSHINTEGER(ad_from0)
!(j < je+1)
      IF (j .NE. je + 1) THEN
        IF (remap_t) THEN
! hydro test
! Remap T in logP
! Note: pt at this stage is Theta_v
          IF (hydrostatic) THEN
! Transform virtual pt to virtual Temp
            DO k=1,km
              ad_from1 = is
              DO i=ad_from1,ie
                CALL PUSHREALARRAY(pt(i, j, k))
                pt(i, j, k) = pt(i, j, k)*(pk(i, j, k+1)-pk(i, j, k))/(&
&                 akap*(peln(i, k+1, j)-peln(i, k, j)))
              END DO
              CALL PUSHINTEGER(i - 1)
              CALL PUSHINTEGER(ad_from1)
            END DO
            CALL PUSHINTEGER(k - 1)
            CALL PUSHCONTROL(3,0)
          ELSE
! Transform "density pt" to "density temp"
            DO k=1,km
              ad_from2 = is
              DO i=ad_from2,ie
                CALL PUSHREALARRAY(pt(i, j, k))
                pt(i, j, k) = pt(i, j, k)*EXP(k1k*LOG(rrg*delp(i, j, k)/&
&                 delz(i, j, k)*pt(i, j, k)))
              END DO
              CALL PUSHINTEGER(i - 1)
              CALL PUSHINTEGER(ad_from2)
            END DO
            CALL PUSHINTEGER(k - 1)
            CALL PUSHCONTROL(3,1)
          END IF
        ELSE IF (remap_pt) THEN
! Using dry pressure for the definition of the virtual potential temperature
!                    pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*    &
!                                              pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
          CALL PUSHCONTROL(3,2)
        ELSE IF (remap_te) THEN
! Remap PT in P
! pt is already virtual PT
! Remap TE in logP
! Transform virtual pt to total energy
          CALL PKEZ_FWD(km, is, ie, js, je, j, pe, pk, akap, peln, pkz, &
&                 ptop)
! Compute cp*T + KE
          DO k=1,km
            ad_from3 = is
            DO i=ad_from3,ie
              CALL PUSHREALARRAY(te(i, j, k))
              te(i, j, k) = 0.25*gridstruct%rsin2(i, j)*(u(i, j, k)**2+u&
&               (i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)&
&               +u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*gridstruct%&
&               cosa_s(i, j)) + cp_air*pt(i, j, k)*pkz(i, j, k)
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from3)
          END DO
          CALL PUSHINTEGER(k - 1)
          CALL PUSHCONTROL(3,3)
        ELSE
          CALL PUSHCONTROL(3,4)
        END IF
        IF (.NOT.hydrostatic) THEN
          DO k=1,km
            ad_from4 = is
            DO i=ad_from4,ie
! ="specific volume"/grav
              CALL PUSHREALARRAY(delz(i, j, k))
              delz(i, j, k) = -(delz(i, j, k)/delp(i, j, k))
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from4)
          END DO
          CALL PUSHINTEGER(k - 1)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        ad_from5 = is
! update ps
        DO i=ad_from5,ie
          ps(i, j) = pe1(i, km+1)
        END DO
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from5)
!
! Hybrid sigma-P coordinate:
!
        DO k=2,km
          ad_from6 = is
          DO i=ad_from6,ie
            CALL PUSHREALARRAY(pe2(i, k))
            pe2(i, k) = ak(k) + bk(k)*pe(i, km+1, j)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from6)
        END DO
        CALL PUSHINTEGER(k - 1)
        DO k=1,km
          ad_from7 = is
          DO i=ad_from7,ie
            CALL PUSHREALARRAY(dp2(i, k))
            dp2(i, k) = pe2(i, k+1) - pe2(i, k)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from7)
        END DO
        CALL PUSHINTEGER(k - 1)
!------------
! update delp
!------------
        DO k=1,km
          ad_from8 = is
          DO i=ad_from8,ie
            CALL PUSHREALARRAY(delp(i, j, k))
            delp(i, j, k) = dp2(i, k)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from8)
        END DO
        CALL PUSHINTEGER(k - 1)
!------------------
! Compute p**Kappa
!------------------
        DO k=1,km+1
          ad_from9 = is
          DO i=ad_from9,ie
            CALL PUSHREALARRAY(pk1(i, k))
            pk1(i, k) = pk(i, j, k)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from9)
        END DO
        CALL PUSHINTEGER(k - 1)
        ad_from10 = is
        DO i=ad_from10,ie
          CALL PUSHREALARRAY(pn2(i, 1))
          pn2(i, 1) = peln(i, 1, j)
          CALL PUSHREALARRAY(pn2(i, km+1))
          pn2(i, km+1) = peln(i, km+1, j)
          CALL PUSHREALARRAY(pk2(i, 1))
          pk2(i, 1) = pk1(i, 1)
          CALL PUSHREALARRAY(pk2(i, km+1))
          pk2(i, km+1) = pk1(i, km+1)
        END DO
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from10)
        DO k=2,km
          ad_from11 = is
          DO i=ad_from11,ie
            CALL PUSHREALARRAY(pn2(i, k))
            pn2(i, k) = LOG(pe2(i, k))
            CALL PUSHREALARRAY(pk2(i, k))
            pk2(i, k) = EXP(akap*pn2(i, k))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from11)
        END DO
        CALL PUSHINTEGER(k - 1)
        IF (remap_t) THEN
!----------------------------------
! Map t using logp
!----------------------------------
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            CALL MAP_SCALAR_FWD(km, peln(is:ie, 1:km+1, j), gz, km, &
&                            pn2, pt, is, ie, j, isd, ied, jsd, jed, 1, &
&                            abs_kord_tm, t_min)
            CALL PUSHCONTROL(3,0)
          ELSE
            CALL PUSHREALARRAY(pt, (ied-isd+1)*(jed-jsd+1)*km)
            CALL MAP_SCALAR(km, peln(is:ie, 1:km+1, j), gz, km, pn2, pt&
&                     , is, ie, j, isd, ied, jsd, jed, 1, &
&                     abs_kord_tm, t_min)
            CALL PUSHCONTROL(3,1)
          END IF
        ELSE IF (remap_pt) THEN
!----------------------------------
! Map pt using pe
!----------------------------------
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            CALL MAP1_PPM_FWD(km, pe1, gz, km, pe2, pt, is, ie, j, &
&                          isd, ied, jsd, jed, 1, abs_kord_tm)
            CALL PUSHCONTROL(3,2)
          ELSE
            CALL PUSHREALARRAY(pt, (ied-isd+1)*(jed-jsd+1)*km)
            CALL MAP1_PPM(km, pe1, gz, km, pe2, pt, is, ie, j, isd, ied&
&                   , jsd, jed, 1, abs_kord_tm)
            CALL PUSHCONTROL(3,3)
          END IF
        ELSE IF (remap_te) THEN
          ad_from12 = is
!----------------------------------
! map Total Energy using GMAO cubic
!----------------------------------
          DO i=ad_from12,ie
            CALL PUSHREALARRAY(phis(i, km+1))
            phis(i, km+1) = hs(i, j)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from12)
          ad_from14 = km
          DO k=ad_from14,1,-1
            ad_from13 = is
            DO i=ad_from13,ie
              CALL PUSHREALARRAY(phis(i, k))
              phis(i, k) = phis(i, k+1) + cp_air*pt(i, j, k)*(pk1(i, k+1&
&               )-pk1(i, k))
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from13)
          END DO
          CALL PUSHINTEGER(ad_from14)
          DO k=1,km+1
            ad_from15 = is
            DO i=ad_from15,ie
              CALL PUSHREALARRAY(phis(i, k))
              phis(i, k) = phis(i, k)*pe1(i, k)
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from15)
          END DO
          CALL PUSHINTEGER(k - 1)
          DO k=1,km
            ad_from16 = is
            DO i=ad_from16,ie
              CALL PUSHREALARRAY(te(i, j, k))
              te(i, j, k) = te(i, j, k) + (phis(i, k+1)-phis(i, k))/(pe1&
&               (i, k+1)-pe1(i, k))
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from16)
          END DO
          CALL PUSHINTEGER(k - 1)
! Map te using log P in GMAO cubic
          CALL MAP1_CUBIC_FWD(km, pe1, km, pe2, te, is, ie, j, isd, ied&
&                       , jsd, jed, akap, t_var=1, conserv=.true.)
          CALL PUSHCONTROL(3,4)
        ELSE
          CALL PUSHCONTROL(3,5)
        END IF
!----------------
! Map constituents
!----------------
        IF (nq .GT. 5) THEN
          IF (kord_tr(1) .EQ. kord_tr_pert(1)) THEN
            CALL MAPN_TRACER_FWD(nq, km, pe1, pe2, q, dp2, kord_tr, j&
&                             , is, ie, isd, ied, jsd, jed, 0., fill)
            CALL PUSHCONTROL(2,0)
          ELSE
            CALL PUSHREALARRAY(q, (ied-isd+1)*(jed-jsd+1)*km*nq)
            CALL MAPN_TRACER(nq, km, pe1, pe2, q, dp2, kord_tr, j, &
&                      is, ie, isd, ied, jsd, jed, 0., fill)
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE IF (nq .GT. 0) THEN
! Remap one tracer at a time
          DO iq=1,nq
            IF (kord_tr(iq) .EQ. kord_tr_pert(iq)) THEN
              CALL MAP1_Q2_FWD(km, pe1, q(isd:ied, jsd:jed, 1:km, iq)&
&                           , km, pe2, q2, dp2, is, ie, 0, kord_tr(iq), &
&                           j, isd, ied, jsd, jed, 0.)
              CALL PUSHCONTROL(1,0)
            ELSE
              CALL MAP1_Q2(km, pe1, q(isd:ied, jsd:jed, 1:km, iq), km, &
&                    pe2, q2, dp2, is, ie, 0, kord_tr(iq), j, isd, &
&                    ied, jsd, jed, 0.)
              CALL PUSHCONTROL(1,1)
            END IF
            DO k=1,km
              ad_from17 = is
              DO i=ad_from17,ie
                CALL PUSHREALARRAY(q(i, j, k, iq))
                q(i, j, k, iq) = q2(i, k)
              END DO
              CALL PUSHINTEGER(i - 1)
              CALL PUSHINTEGER(ad_from17)
            END DO
            CALL PUSHINTEGER(k - 1)
          END DO
          CALL PUSHINTEGER(iq - 1)
          CALL PUSHCONTROL(2,2)
        ELSE
          CALL PUSHCONTROL(2,3)
        END IF
        IF (.NOT.hydrostatic) THEN
! Remap vertical wind:
          IF (kord_wz .EQ. kord_wz_pert) THEN
            CALL MAP1_PPM_FWD(km, pe1, ws(is:ie, j), km, pe2, w, is, &
&                          ie, j, isd, ied, jsd, jed, -2, kord_wz)
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHREALARRAY(w, (ied-isd+1)*(jed-jsd+1)*km)
            CALL MAP1_PPM(km, pe1, ws(is:ie, j), km, pe2, w, is, ie, j, &
&                   isd, ied, jsd, jed, -2, kord_wz)
            CALL PUSHCONTROL(1,1)
          END IF
! Remap delz for hybrid sigma-p coordinate
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            CALL MAP1_PPM_FWD(km, pe1, gz, km, pe2, delz, is, ie, j, &
&                          isd, ied, jsd, jed, 1, abs_kord_tm)
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHREALARRAY(delz, (ied-isd+1)*(jed-jsd+1)*km)
            CALL MAP1_PPM(km, pe1, gz, km, pe2, delz, is, ie, j, isd, &
&                   ied, jsd, jed, 1, abs_kord_tm)
            CALL PUSHCONTROL(1,0)
          END IF
          DO k=1,km
            ad_from18 = is
            DO i=ad_from18,ie
              CALL PUSHREALARRAY(delz(i, j, k))
              delz(i, j, k) = -(delz(i, j, k)*dp2(i, k))
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from18)
          END DO
          CALL PUSHINTEGER(k - 1)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
!----------
! Update pk
!----------
        DO k=1,km+1
          ad_from19 = is
          DO i=ad_from19,ie
            CALL PUSHREALARRAY(pk(i, j, k))
            pk(i, j, k) = pk2(i, k)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from19)
        END DO
        CALL PUSHINTEGER(k - 1)
!----------------
        IF (do_omega) THEN
          ad_from20 = is
! Start do_omega
! Copy omega field to pe3
          DO i=ad_from20,ie
            CALL PUSHREALARRAY(pe3(i, 1))
            pe3(i, 1) = 0.
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from20)
          DO k=2,km+1
            ad_from21 = is
            DO i=ad_from21,ie
              CALL PUSHREALARRAY(pe3(i, k))
              pe3(i, k) = omga(i, j, k-1)
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from21)
          END DO
          CALL PUSHINTEGER(k - 1)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        DO k=1,km+1
          ad_from22 = is
          DO i=ad_from22,ie
            CALL PUSHREALARRAY(pe0(i, k))
            pe0(i, k) = peln(i, k, j)
            CALL PUSHREALARRAY(peln(i, k, j))
            peln(i, k, j) = pn2(i, k)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from22)
        END DO
        CALL PUSHINTEGER(k - 1)
!------------
! Compute pkz
!------------
        IF (hydrostatic) THEN
          DO k=1,km
            ad_from23 = is
            DO i=ad_from23,ie
              CALL PUSHREALARRAY(pkz(i, j, k))
              pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1&
&               , j)-peln(i, k, j)))
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from23)
          END DO
          CALL PUSHINTEGER(k - 1)
          CALL PUSHCONTROL(2,0)
        ELSE IF (remap_te) THEN
! WMP: note that this is where TE remapping non-hydrostatic is invalid and cannot be run
          GOTO 140
        ELSE IF (remap_t) THEN
! Note: pt at this stage is T_v or T_m
          DO k=1,km
            ad_from24 = is
            DO i=ad_from24,ie
              CALL PUSHREALARRAY(pkz(i, j, k))
              pkz(i, j, k) = EXP(akap*LOG(rrg*delp(i, j, k)/delz(i, j, k&
&               )*pt(i, j, k)))
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from24)
          END DO
          CALL PUSHINTEGER(k - 1)
          CALL PUSHCONTROL(2,1)
        ELSE
! Using dry pressure for the definition of the virtual potential temperature
!           pkz(i,j,k) = exp(akap*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
! Note: pt at this stage is Theta_v
          DO k=1,km
            ad_from25 = is
            DO i=ad_from25,ie
              CALL PUSHREALARRAY(pkz(i, j, k))
              pkz(i, j, k) = EXP(k1k*LOG(rrg*delp(i, j, k)/delz(i, j, k)&
&               *pt(i, j, k)))
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from25)
          END DO
          CALL PUSHINTEGER(k - 1)
          CALL PUSHCONTROL(2,2)
        END IF
! end do_omega
! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
        IF (do_omega) THEN
          DO k=1,km
            ad_from26 = is
            DO i=ad_from26,ie
              CALL PUSHREALARRAY(dp2(i, k))
              dp2(i, k) = 0.5*(peln(i, k, j)+peln(i, k+1, j))
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from26)
          END DO
          CALL PUSHINTEGER(k - 1)
          ad_from27 = is
          DO i=ad_from27,ie
            k_next = 1
            DO 130 n=1,km
              kp = k_next
              CALL PUSHINTEGER(k)
              ad_count0 = 1
              DO k=kp,km
                IF (dp2(i, n) .LE. pe0(i, k+1) .AND. dp2(i, n) .GE. pe0(&
&                   i, k)) THEN
                  GOTO 120
                ELSE
                  CALL PUSHINTEGER(k)
                  ad_count0 = ad_count0 + 1
                END IF
              END DO
              CALL PUSHCONTROL(1,0)
              CALL PUSHINTEGER(ad_count0)
              CALL PUSHCONTROL(1,1)
              GOTO 130
 120          CALL PUSHCONTROL(1,1)
              CALL PUSHINTEGER(ad_count0)
              CALL PUSHREALARRAY(omga(i, j, n))
              omga(i, j, n) = pe3(i, k) + (pe3(i, k+1)-pe3(i, k))*(dp2(i&
&               , n)-pe0(i, k))/(pe0(i, k+1)-pe0(i, k))
              k_next = k
              CALL PUSHCONTROL(1,0)
 130        CONTINUE
            CALL PUSHINTEGER(n - 1)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from27)
          CALL PUSHCONTROL(2,2)
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,0)
      END IF
      ad_from28 = is
      DO i=ad_from28,ie+1
        CALL PUSHREALARRAY(pe0(i, 1))
        pe0(i, 1) = pe(i, 1, j)
      END DO
      CALL PUSHINTEGER(i - 1)
      CALL PUSHINTEGER(ad_from28)
      CALL PUSHINTEGER(k)
!------
! map u
!------
      DO k=2,km+1
        ad_from29 = is
        DO i=ad_from29,ie
          CALL PUSHREALARRAY(pe0(i, k))
          pe0(i, k) = 0.5*(pe(i, k, j-1)+pe1(i, k))
        END DO
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from29)
      END DO
      CALL PUSHINTEGER(k - 1)
      DO k=1,km+1
        CALL PUSHREALARRAY(bkh)
        bkh = 0.5*bk(k)
        ad_from30 = is
        DO i=ad_from30,ie
          CALL PUSHREALARRAY(pe3(i, k))
          pe3(i, k) = ak(k) + bkh*(pe(i, km+1, j-1)+pe1(i, km+1))
        END DO
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from30)
      END DO
      CALL PUSHINTEGER(k - 1)
      IF (kord_mt .EQ. kord_mt_pert) THEN
        CALL MAP1_PPM_FWD(km, pe0(is:ie, :), gz, km, pe3(is:ie, :), u&
&                      , is, ie, j, isd, ied, jsd, jedp1, -1, kord_mt)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHREALARRAY(u, (ied-isd+1)*(jed-jsd+2)*km)
        CALL MAP1_PPM(km, pe0(is:ie, :), gz, km, pe3(is:ie, :), u, is, &
&               ie, j, isd, ied, jsd, jedp1, -1, kord_mt)
        CALL PUSHCONTROL(1,1)
      END IF
! (j < je+1)
      IF (j .LT. je + 1) THEN
        ad_from31 = is
!------
! map v
!------
        DO i=ad_from31,ie+1
          CALL PUSHREALARRAY(pe3(i, 1))
          pe3(i, 1) = ak(1)
        END DO
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from31)
        DO k=2,km+1
          CALL PUSHREALARRAY(bkh)
          bkh = 0.5*bk(k)
          ad_from32 = is
          DO i=ad_from32,ie+1
            CALL PUSHREALARRAY(pe0(i, k))
            pe0(i, k) = 0.5*(pe(i-1, k, j)+pe(i, k, j))
            CALL PUSHREALARRAY(pe3(i, k))
            pe3(i, k) = ak(k) + bkh*(pe(i-1, km+1, j)+pe(i, km+1, j))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from32)
        END DO
        CALL PUSHINTEGER(k - 1)
        IF (kord_mt .EQ. kord_mt_pert) THEN
          CALL MAP1_PPM_FWD(km, pe0, gz, km, pe3, v, is, iep1, j, isd&
&                        , iedp1, jsd, jed, -1, kord_mt)
          CALL PUSHCONTROL(2,1)
        ELSE
          CALL PUSHREALARRAY(v, (ied-isd+2)*(jed-jsd+1)*km)
          CALL MAP1_PPM(km, pe0, gz, km, pe3, v, is, iep1, j, isd, iedp1&
&                 , jsd, jed, -1, kord_mt)
          CALL PUSHCONTROL(2,2)
        END IF
      ELSE
        CALL PUSHCONTROL(2,0)
      END IF
      DO k=1,km
        ad_from33 = is
        DO i=ad_from33,ie
          CALL PUSHREALARRAY(ua(i, j, k))
          ua(i, j, k) = pe2(i, k+1)
        END DO
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from33)
      END DO
      CALL PUSHINTEGER(k - 1)
      CALL PUSHINTEGER(j)
      ad_count1 = ad_count1 + 1
    END DO
    CALL PUSHCONTROL(1,0)
    CALL PUSHINTEGER(ad_count1)
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
          CALL PUSHREALARRAY(pe(i, k, j))
          pe(i, k, j) = ua(i, j, k-1)
        END DO
      END DO
    END DO
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
                CALL PUSHREALARRAY(gz(i))
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
              CALL PUSHCONTROL(3,5)
            ELSE
              DO i=is,ie
                te_2d(i, j) = 0.
                CALL PUSHREALARRAY(phis(i, km+1))
                phis(i, km+1) = hs(i, j)
              END DO
              DO k=km,1,-1
                DO i=is,ie
                  CALL PUSHREALARRAY(phis(i, k))
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
              CALL PUSHCONTROL(3,4)
            END IF
          ELSE IF (remap_pt) THEN
! k-loop
            IF (hydrostatic) THEN
              DO i=is,ie
                CALL PUSHREALARRAY(gz(i))
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
              CALL PUSHCONTROL(3,3)
            ELSE
!-----------------
! Non-hydrostatic:
!-----------------
              DO i=is,ie
                CALL PUSHREALARRAY(phis(i, km+1))
                phis(i, km+1) = hs(i, j)
                DO k=km,1,-1
                  CALL PUSHREALARRAY(phis(i, k))
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
              CALL PUSHCONTROL(3,2)
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
            CALL PUSHCONTROL(3,1)
          ELSE
            CALL PUSHCONTROL(3,0)
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
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
! j-loop
!$OMP single
        result1 = G_SUM(domain, te_2d, is, ie, js, je, ng, gridstruct%&
&         area_64, 0, reproduce=.true.)
        tpe = consv*result1
! unit: W/m**2
! Note pdt is "phys" time step
        IF (hydrostatic) THEN
          result1 = G_SUM(domain, zsum0, is, ie, js, je, ng, gridstruct%&
&           area_64, 0, reproduce=.true.)
          dtmp = tpe/(cp*result1)
          CALL PUSHCONTROL(3,0)
        ELSE
          result1 = G_SUM(domain, zsum1, is, ie, js, je, ng, gridstruct%&
&           area_64, 0, reproduce=.true.)
          dtmp = tpe/(cv_air*result1)
          CALL PUSHCONTROL(3,1)
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
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        CALL PUSHREALARRAY(e_flux)
        e_flux = consv
!$OMP single
        IF (hydrostatic) THEN
          result1 = G_SUM(domain, zsum0, is, ie, js, je, ng, gridstruct%&
&           area_64, 0, reproduce=.true.)
          dtmp = e_flux*(grav*pdt*4.*pi*radius**2)/(cp*result1)
          CALL PUSHCONTROL(3,2)
        ELSE
          result1 = G_SUM(domain, zsum1, is, ie, js, je, ng, gridstruct%&
&           area_64, 0, reproduce=.true.)
          dtmp = e_flux*(grav*pdt*4.*pi*radius**2)/(cv_air*result1)
          CALL PUSHCONTROL(3,3)
        END IF
      ELSE
        CALL PUSHCONTROL(3,4)
      END IF
    ELSE
      CALL PUSHCONTROL(3,5)
    END IF
!$OMP end single
! do_sat_adj
! Note: pt at this stage is T_v
    IF (remap_t .AND. (.NOT.do_adiabatic_init) .AND. do_sat_adj) THEN
!$OMP do
      DO k=kmp,km
        IF (.NOT.hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(pkz(i, j, k))
              pkz(i, j, k) = EXP(akap*LOG(rrg*delp(i, j, k)/delz(i, j, k&
&               )*pt(i, j, k)))
            END DO
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
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
        CALL PUSHCONTROL(2,0)
      ELSE
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHCONTROL(2,2)
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
                CALL PUSHREALARRAY(pt(i, j, k))
                pt(i, j, k) = (pt(i, j, k)+dtmp*pkz(i, j, k))/(1.+r_vir*&
&                 q(i, j, k, sphum))
              END DO
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
        END DO
        CALL PUSHREALARRAY(pk2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(bkh)
        CALL PUSHREALARRAY(pk1, (ie-is+1)*(km+1))
        CALL PUSHINTEGER(liq_wat)
        CALL PUSHREALARRAY(zsum1, (ie-is+1)*(je-js+1))
        CALL PUSHREALARRAY(pn2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(zsum0, (ie-is+1)*(je-js+1))
        CALL PUSHREALARRAY(rrg)
        CALL PUSHINTEGER(ice_wat)
        CALL PUSHREALARRAY(result1)
        CALL PUSHINTEGER(rainwat)
        CALL PUSHREALARRAY(te_2d, (ie-is+1)*(je-js+1))
        CALL PUSHINTEGER(abs_kord_tm_pert)
        CALL PUSHINTEGER(cld_amt)
        CALL PUSHREALARRAY(dtmp)
        CALL PUSHINTEGER(snowwat)
        CALL PUSHREALARRAY(k1k)
        CALL PUSHREALARRAY(tpe)
        CALL PUSHREALARRAY(rg)
        CALL PUSHINTEGER(iep1)
        CALL PUSHREALARRAY(dp2, (ie-is+1)*km)
        CALL PUSHINTEGER(kmp)
        CALL PUSHREALARRAY(gz, ie - is + 1)
        CALL PUSHREALARRAY(pe3, (ie-is+2)*(km+1))
        CALL PUSHREALARRAY(pe2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(pe1, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(pe0, (ie-is+2)*(km+1))
        CALL PUSHREALARRAY(phis, (ie-is+1)*(km+1))
        CALL PUSHINTEGER(graupel)
        CALL PUSHCONTROL(3,0)
      ELSE IF (remap_pt) THEN
! j-loop
! k-loop
!$OMP do
        DO k=1,km
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(pt(i, j, k))
              pt(i, j, k) = (pt(i, j, k)+dtmp)*pkz(i, j, k)/(1.+r_vir*q(&
&               i, j, k, sphum))
            END DO
          END DO
        END DO
        CALL PUSHREALARRAY(pk2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(bkh)
        CALL PUSHREALARRAY(pk1, (ie-is+1)*(km+1))
        CALL PUSHINTEGER(liq_wat)
        CALL PUSHREALARRAY(zsum1, (ie-is+1)*(je-js+1))
        CALL PUSHREALARRAY(pn2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(zsum0, (ie-is+1)*(je-js+1))
        CALL PUSHREALARRAY(rrg)
        CALL PUSHINTEGER(ice_wat)
        CALL PUSHREALARRAY(result1)
        CALL PUSHINTEGER(rainwat)
        CALL PUSHREALARRAY(te_2d, (ie-is+1)*(je-js+1))
        CALL PUSHINTEGER(abs_kord_tm_pert)
        CALL PUSHINTEGER(cld_amt)
        CALL PUSHREALARRAY(dtmp)
        CALL PUSHINTEGER(snowwat)
        CALL PUSHREALARRAY(k1k)
        CALL PUSHREALARRAY(tpe)
        CALL PUSHREALARRAY(rg)
        CALL PUSHINTEGER(iep1)
        CALL PUSHREALARRAY(dp2, (ie-is+1)*km)
        CALL PUSHINTEGER(kmp)
        CALL PUSHREALARRAY(gz, ie - is + 1)
        CALL PUSHREALARRAY(pe3, (ie-is+2)*(km+1))
        CALL PUSHREALARRAY(pe2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(pe1, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(pe0, (ie-is+2)*(km+1))
        CALL PUSHREALARRAY(phis, (ie-is+1)*(km+1))
        CALL PUSHINTEGER(graupel)
        CALL PUSHCONTROL(3,1)
      ELSE IF (remap_te) THEN
!$OMP do
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(gz(i))
            gz(i) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              CALL PUSHREALARRAY(tpe)
              tpe = te(i, j, k) - gz(i) - 0.25*gridstruct%rsin2(i, j)*(u&
&               (i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)&
&               **2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&               gridstruct%cosa_s(i, j))
              dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
              CALL PUSHREALARRAY(tmp)
              tmp = tpe/((cp-pe(i, k, j)*dlnp/delp(i, j, k))*(1.+r_vir*q&
&               (i, j, k, sphum)))
              CALL PUSHREALARRAY(pt(i, j, k))
              pt(i, j, k) = tmp + dtmp*pkz(i, j, k)/(1.+r_vir*q(i, j, k&
&               , sphum))
              CALL PUSHREALARRAY(gz(i))
              gz(i) = gz(i) + dlnp*tmp*(1.+r_vir*q(i, j, k, sphum))
            END DO
          END DO
        END DO
        CALL PUSHREALARRAY(pk2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(bkh)
        CALL PUSHREALARRAY(pk1, (ie-is+1)*(km+1))
        CALL PUSHINTEGER(liq_wat)
        CALL PUSHREALARRAY(zsum1, (ie-is+1)*(je-js+1))
        CALL PUSHREALARRAY(pn2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(zsum0, (ie-is+1)*(je-js+1))
        CALL PUSHREALARRAY(rrg)
        CALL PUSHINTEGER(ice_wat)
        CALL PUSHREALARRAY(result1)
        CALL PUSHINTEGER(rainwat)
        CALL PUSHREALARRAY(te_2d, (ie-is+1)*(je-js+1))
        CALL PUSHINTEGER(abs_kord_tm_pert)
        CALL PUSHINTEGER(cld_amt)
        CALL PUSHREALARRAY(tmp)
        CALL PUSHREALARRAY(dtmp)
        CALL PUSHINTEGER(snowwat)
        CALL PUSHREALARRAY(k1k)
        CALL PUSHREALARRAY(tpe)
        CALL PUSHREALARRAY(rg)
        CALL PUSHINTEGER(iep1)
        CALL PUSHREALARRAY(dp2, (ie-is+1)*km)
        CALL PUSHINTEGER(kmp)
        CALL PUSHREALARRAY(gz, ie - is + 1)
        CALL PUSHREALARRAY(pe3, (ie-is+2)*(km+1))
        CALL PUSHREALARRAY(pe2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(pe1, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(pe0, (ie-is+2)*(km+1))
        CALL PUSHREALARRAY(phis, (ie-is+1)*(km+1))
        CALL PUSHINTEGER(graupel)
        CALL PUSHCONTROL(3,2)
      ELSE
        CALL PUSHREALARRAY(pk2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(bkh)
        CALL PUSHREALARRAY(pk1, (ie-is+1)*(km+1))
        CALL PUSHINTEGER(liq_wat)
        CALL PUSHREALARRAY(zsum1, (ie-is+1)*(je-js+1))
        CALL PUSHREALARRAY(pn2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(zsum0, (ie-is+1)*(je-js+1))
        CALL PUSHREALARRAY(rrg)
        CALL PUSHINTEGER(ice_wat)
        CALL PUSHREALARRAY(result1)
        CALL PUSHINTEGER(rainwat)
        CALL PUSHREALARRAY(te_2d, (ie-is+1)*(je-js+1))
        CALL PUSHINTEGER(abs_kord_tm_pert)
        CALL PUSHINTEGER(cld_amt)
        CALL PUSHINTEGER(snowwat)
        CALL PUSHREALARRAY(k1k)
        CALL PUSHREALARRAY(tpe)
        CALL PUSHREALARRAY(rg)
        CALL PUSHINTEGER(iep1)
        CALL PUSHREALARRAY(dp2, (ie-is+1)*km)
        CALL PUSHINTEGER(kmp)
        CALL PUSHREALARRAY(gz, ie - is + 1)
        CALL PUSHREALARRAY(pe3, (ie-is+2)*(km+1))
        CALL PUSHREALARRAY(pe2, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(pe1, (ie-is+1)*(km+1))
        CALL PUSHREALARRAY(pe0, (ie-is+2)*(km+1))
        CALL PUSHREALARRAY(phis, (ie-is+1)*(km+1))
        CALL PUSHINTEGER(graupel)
        CALL PUSHCONTROL(3,3)
      END IF
    ELSE IF (remap_t) THEN
! not last_step
!$OMP do
      DO k=1,km
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(pt(i, j, k))
            pt(i, j, k) = pt(i, j, k)/pkz(i, j, k)
          END DO
        END DO
      END DO
      CALL PUSHREALARRAY(pk2, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(bkh)
      CALL PUSHREALARRAY(pk1, (ie-is+1)*(km+1))
      CALL PUSHINTEGER(liq_wat)
      CALL PUSHREALARRAY(zsum1, (ie-is+1)*(je-js+1))
      CALL PUSHREALARRAY(pn2, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(zsum0, (ie-is+1)*(je-js+1))
      CALL PUSHREALARRAY(rrg)
      CALL PUSHINTEGER(ice_wat)
      CALL PUSHREALARRAY(result1)
      CALL PUSHINTEGER(rainwat)
      CALL PUSHREALARRAY(te_2d, (ie-is+1)*(je-js+1))
      CALL PUSHINTEGER(abs_kord_tm_pert)
      CALL PUSHINTEGER(cld_amt)
      CALL PUSHINTEGER(snowwat)
      CALL PUSHREALARRAY(k1k)
      CALL PUSHREALARRAY(tpe)
      CALL PUSHREALARRAY(rg)
      CALL PUSHINTEGER(iep1)
      CALL PUSHREALARRAY(dp2, (ie-is+1)*km)
      CALL PUSHINTEGER(kmp)
      CALL PUSHREALARRAY(gz, ie - is + 1)
      CALL PUSHREALARRAY(pe3, (ie-is+2)*(km+1))
      CALL PUSHREALARRAY(pe2, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(pe1, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(pe0, (ie-is+2)*(km+1))
      CALL PUSHREALARRAY(phis, (ie-is+1)*(km+1))
      CALL PUSHINTEGER(graupel)
      CALL PUSHCONTROL(3,4)
    ELSE IF (remap_te) THEN
!$OMP do
      DO j=js,je
        DO i=is,ie
          CALL PUSHREALARRAY(gz(i))
          gz(i) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            CALL PUSHREALARRAY(tpe)
            tpe = te(i, j, k) - gz(i) - 0.25*gridstruct%rsin2(i, j)*(u(i&
&             , j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(&
&             u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*&
&             gridstruct%cosa_s(i, j))
            dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
            tmp = tpe/(cp-pe(i, k, j)*dlnp/delp(i, j, k))
            CALL PUSHREALARRAY(pt(i, j, k))
            pt(i, j, k) = tmp/pkz(i, j, k) + dtmp
            CALL PUSHREALARRAY(gz(i))
            gz(i) = gz(i) + dlnp*tmp
          END DO
        END DO
      END DO
      CALL PUSHREALARRAY(pk2, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(bkh)
      CALL PUSHREALARRAY(pk1, (ie-is+1)*(km+1))
      CALL PUSHINTEGER(liq_wat)
      CALL PUSHREALARRAY(zsum1, (ie-is+1)*(je-js+1))
      CALL PUSHREALARRAY(pn2, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(zsum0, (ie-is+1)*(je-js+1))
      CALL PUSHREALARRAY(rrg)
      CALL PUSHINTEGER(ice_wat)
      CALL PUSHREALARRAY(result1)
      CALL PUSHINTEGER(rainwat)
      CALL PUSHREALARRAY(te_2d, (ie-is+1)*(je-js+1))
      CALL PUSHINTEGER(abs_kord_tm_pert)
      CALL PUSHINTEGER(cld_amt)
      CALL PUSHINTEGER(snowwat)
      CALL PUSHREALARRAY(k1k)
      CALL PUSHREALARRAY(tpe)
      CALL PUSHREALARRAY(rg)
      CALL PUSHINTEGER(iep1)
      CALL PUSHREALARRAY(dp2, (ie-is+1)*km)
      CALL PUSHINTEGER(kmp)
      CALL PUSHREALARRAY(gz, ie - is + 1)
      CALL PUSHREALARRAY(pe3, (ie-is+2)*(km+1))
      CALL PUSHREALARRAY(pe2, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(pe1, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(pe0, (ie-is+2)*(km+1))
      CALL PUSHREALARRAY(phis, (ie-is+1)*(km+1))
      CALL PUSHINTEGER(graupel)
      CALL PUSHCONTROL(3,5)
    ELSE
      CALL PUSHREALARRAY(pk2, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(bkh)
      CALL PUSHREALARRAY(pk1, (ie-is+1)*(km+1))
      CALL PUSHINTEGER(liq_wat)
      CALL PUSHREALARRAY(zsum1, (ie-is+1)*(je-js+1))
      CALL PUSHREALARRAY(pn2, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(zsum0, (ie-is+1)*(je-js+1))
      CALL PUSHREALARRAY(rrg)
      CALL PUSHINTEGER(ice_wat)
      CALL PUSHREALARRAY(result1)
      CALL PUSHINTEGER(rainwat)
      CALL PUSHREALARRAY(te_2d, (ie-is+1)*(je-js+1))
      CALL PUSHINTEGER(abs_kord_tm_pert)
      CALL PUSHINTEGER(cld_amt)
      CALL PUSHINTEGER(snowwat)
      CALL PUSHREALARRAY(k1k)
      CALL PUSHREALARRAY(tpe)
      CALL PUSHREALARRAY(rg)
      CALL PUSHINTEGER(iep1)
      CALL PUSHREALARRAY(dp2, (ie-is+1)*km)
      CALL PUSHINTEGER(kmp)
      CALL PUSHREALARRAY(gz, ie - is + 1)
      CALL PUSHREALARRAY(pe3, (ie-is+2)*(km+1))
      CALL PUSHREALARRAY(pe2, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(pe1, (ie-is+1)*(km+1))
      CALL PUSHREALARRAY(pe0, (ie-is+2)*(km+1))
      CALL PUSHREALARRAY(phis, (ie-is+1)*(km+1))
      CALL PUSHINTEGER(graupel)
      CALL PUSHCONTROL(3,6)
    END IF
    GOTO 150
 140 CALL PUSHCONTROL(1,1)
    CALL PUSHINTEGER(ad_count1)
    STOP
 150 CONTINUE
  END SUBROUTINE LAGRANGIAN_TO_EULERIAN_FWD
!  Differentiation of lagrangian_to_eulerian in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4_
!fb a2b_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv
!_pe dyn_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update d
!yn_core_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics
!_mod.Rayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_util
!s_mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez f
!v_mapz_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_m
!od.remap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.
!ppm_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.ma
!p1_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_
!mod.fv_subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_d
!z_d nh_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solve
!r nh_utils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile
! nh_utils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest s
!w_core_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_co
!re_mod.ytp_v sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2
!d_fb tp_core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner 
!fv_grid_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: ws peln q u v w delp ua delz
!                omga te0_2d pkz pe pk ps pt te
!   with respect to varying inputs: ws peln q u v w delp ua delz
!                omga te0_2d pkz pe pk ps pt te
  SUBROUTINE LAGRANGIAN_TO_EULERIAN_BWD(last_step, consv, ps, ps_ad, pe&
&   , pe_ad, delp, delp_ad, pkz, pkz_ad, pk, pk_ad, mdt, pdt, km, is, ie&
&   , js, je, isd, ied, jsd, jed, nq, nwat, sphum, q_con, u, u_ad, v, &
&   v_ad, w, w_ad, delz, delz_ad, pt, pt_ad, q, q_ad, hs, r_vir, cp, &
&   akap, cappa, kord_mt, kord_wz, kord_tr, kord_tm, peln, peln_ad, &
&   te0_2d, te0_2d_ad, ng, ua, ua_ad, va, omga, omga_ad, te, te_ad, ws, &
&   ws_ad, fill, reproduce_sum, out_dt, dtdt, ptop, ak, bk, pfull, &
&   flagstruct, gridstruct, domain, do_sat_adj, hydrostatic, hybrid_z, &
&   do_omega, adiabatic, do_adiabatic_init, mfx, mfy, remap_option, &
&   kord_mt_pert, kord_wz_pert, kord_tr_pert, kord_tm_pert)
    IMPLICIT NONE
!$OMP end parallel
    LOGICAL, INTENT(IN) :: last_step
    REAL, INTENT(IN) :: mdt
    REAL, INTENT(IN) :: pdt
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: nwat
    INTEGER, INTENT(IN) :: sphum
    INTEGER, INTENT(IN) :: ng
    INTEGER, INTENT(IN) :: is, ie, isd, ied
    INTEGER, INTENT(IN) :: js, je, jsd, jed
    INTEGER, INTENT(IN) :: kord_mt
    INTEGER, INTENT(IN) :: kord_wz
    INTEGER, INTENT(IN) :: kord_tr(nq)
    INTEGER, INTENT(IN) :: kord_tm
    INTEGER, INTENT(IN) :: kord_mt_pert
    INTEGER, INTENT(IN) :: kord_wz_pert
    INTEGER, INTENT(IN) :: kord_tr_pert(nq)
    INTEGER, INTENT(IN) :: kord_tm_pert
    REAL, INTENT(IN) :: consv
    REAL, INTENT(IN) :: r_vir
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: akap
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: te0_2d(is:ie, js:je)
    REAL, INTENT(INOUT) :: te0_2d_ad(is:ie, js:je)
    REAL, INTENT(IN) :: ws(is:ie, js:je)
    REAL :: ws_ad(is:ie, js:je)
    LOGICAL, INTENT(IN) :: do_sat_adj
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
    REAL, INTENT(INOUT) :: pk(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: pk_ad(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: q_ad(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: delp_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(INOUT) :: pe_ad(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(INOUT) :: ps(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: ps_ad(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: u_ad(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_ad(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: w_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: pt_ad(isd:ied, jsd:jed, km)
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: delz
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: delz_ad
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: q_con, cappa
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: hybrid_z
    LOGICAL, INTENT(IN) :: out_dt
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: ua_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: va(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: omga(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: omga_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: peln(is:ie, km+1, js:je)
    REAL, INTENT(INOUT) :: peln_ad(is:ie, km+1, js:je)
    REAL, INTENT(INOUT) :: dtdt(is:ie, js:je, km)
    REAL :: pkz(is:ie, js:je, km)
    REAL :: pkz_ad(is:ie, js:je, km)
    REAL :: te(isd:ied, jsd:jed, km)
    REAL :: te_ad(isd:ied, jsd:jed, km)
    REAL, OPTIONAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, km)
    REAL, OPTIONAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, km)
    INTEGER, INTENT(IN) :: remap_option
    REAL, DIMENSION(is:ie, js:je) :: te_2d, zsum0, zsum1, dpln
    REAL, DIMENSION(is:ie, js:je) :: te_2d_ad, zsum0_ad, zsum1_ad
    REAL, DIMENSION(is:ie, km) :: q2, dp2
    REAL, DIMENSION(is:ie, km) :: q2_ad, dp2_ad
    REAL, DIMENSION(is:ie, km+1) :: pe1, pe2, pk1, pk2, pn2, phis
    REAL, DIMENSION(is:ie, km+1) :: pe1_ad, pe2_ad, pk1_ad, pk2_ad, &
&   pn2_ad, phis_ad
    REAL, DIMENSION(is:ie+1, km+1) :: pe0, pe3
    REAL, DIMENSION(is:ie+1, km+1) :: pe0_ad, pe3_ad
    REAL, DIMENSION(is:ie) :: gz, cvm, qv
    REAL, DIMENSION(is:ie) :: gz_ad
    REAL :: rcp, rg, tmp, tpe, rrg, bkh, dtmp, k1k, dlnp
    REAL :: tmp_ad, tpe_ad, dtmp_ad, dlnp_ad
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
    INTEGER :: arg1
    REAL :: result1
    REAL :: result1_ad
    LOGICAL :: arg10
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp2
    REAL :: temp3
    REAL :: temp4
    REAL :: temp5
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp6
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp7
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp8
    REAL :: temp9
    REAL :: temp10
    REAL :: temp_ad12
    REAL :: temp11
    REAL :: temp12
    REAL :: temp13
    REAL :: temp_ad13
    REAL :: temp14
    REAL :: temp15
    REAL :: temp16
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: temp_ad16
    REAL :: temp17
    REAL :: temp18
    REAL :: temp19
    REAL :: temp_ad17
    REAL :: temp_ad18
    REAL :: temp_ad19
    REAL :: temp20
    REAL :: temp21
    REAL :: temp22
    REAL :: temp23
    REAL :: temp24
    REAL :: temp_ad20
    REAL :: temp_ad21
    REAL :: temp_ad22
    REAL :: temp_ad23
    REAL :: temp_ad24
    REAL :: temp_ad25
    REAL :: temp25
    REAL :: temp26
    REAL :: temp27
    REAL :: temp_ad26
    REAL :: temp_ad27
    REAL :: temp_ad28
    REAL :: temp28
    REAL :: temp29
    REAL :: temp30
    REAL :: temp31
    REAL :: temp32
    REAL :: temp_ad29
    REAL :: temp_ad30
    REAL :: temp_ad31
    REAL :: temp_ad32
    REAL :: temp_ad33
    REAL :: temp_ad34
    REAL :: temp_ad35
    REAL :: temp33
    REAL :: temp34
    REAL :: temp35
    REAL :: temp_ad36
    REAL :: temp36
    REAL :: temp_ad37
    REAL :: temp37
    REAL :: temp38
    REAL :: temp39
    REAL :: temp_ad38
    REAL :: temp40
    REAL :: temp41
    REAL :: temp42
    REAL :: temp43
    REAL :: temp44
    REAL :: temp45
    REAL :: temp_ad39
    REAL :: temp_ad40
    REAL :: temp_ad41
    REAL :: temp_ad42
    REAL :: temp_ad43
    REAL :: temp_ad44
    REAL :: temp_ad45
    REAL :: temp_ad46
    REAL :: temp_ad47
    REAL :: temp46
    REAL :: temp47
    REAL :: temp48
    REAL :: temp_ad48
    REAL :: temp_ad49
    REAL :: temp_ad50
    REAL :: temp_ad51
    REAL :: temp_ad52
    REAL :: temp_ad53
    REAL :: temp_ad54
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_to0
    INTEGER :: ad_from0
    INTEGER :: ad_to1
    INTEGER :: ad_from1
    INTEGER :: ad_to2
    INTEGER :: ad_to3
    INTEGER :: ad_from2
    INTEGER :: ad_to4
    INTEGER :: ad_to5
    INTEGER :: ad_from3
    INTEGER :: ad_to6
    INTEGER :: ad_to7
    INTEGER :: ad_from4
    INTEGER :: ad_to8
    INTEGER :: ad_to9
    INTEGER :: ad_from5
    INTEGER :: ad_to10
    INTEGER :: ad_from6
    INTEGER :: ad_to11
    INTEGER :: ad_to12
    INTEGER :: ad_from7
    INTEGER :: ad_to13
    INTEGER :: ad_to14
    INTEGER :: ad_from8
    INTEGER :: ad_to15
    INTEGER :: ad_to16
    INTEGER :: ad_from9
    INTEGER :: ad_to17
    INTEGER :: ad_to18
    INTEGER :: ad_from10
    INTEGER :: ad_to19
    INTEGER :: ad_from11
    INTEGER :: ad_to20
    INTEGER :: ad_to21
    INTEGER :: ad_from12
    INTEGER :: ad_to22
    INTEGER :: ad_from13
    INTEGER :: ad_to23
    INTEGER :: ad_from14
    INTEGER :: ad_from15
    INTEGER :: ad_to24
    INTEGER :: ad_to25
    INTEGER :: ad_from16
    INTEGER :: ad_to26
    INTEGER :: ad_to27
    INTEGER :: ad_from17
    INTEGER :: ad_to28
    INTEGER :: ad_to29
    INTEGER :: ad_to30
    INTEGER :: ad_from18
    INTEGER :: ad_to31
    INTEGER :: ad_to32
    INTEGER :: ad_from19
    INTEGER :: ad_to33
    INTEGER :: ad_to34
    INTEGER :: ad_from20
    INTEGER :: ad_to35
    INTEGER :: ad_from21
    INTEGER :: ad_to36
    INTEGER :: ad_to37
    INTEGER :: ad_from22
    INTEGER :: ad_to38
    INTEGER :: ad_to39
    INTEGER :: ad_from23
    INTEGER :: ad_to40
    INTEGER :: ad_to41
    INTEGER :: ad_from24
    INTEGER :: ad_to42
    INTEGER :: ad_to43
    INTEGER :: ad_from25
    INTEGER :: ad_to44
    INTEGER :: ad_to45
    INTEGER :: ad_from26
    INTEGER :: ad_to46
    INTEGER :: ad_to47
    INTEGER :: ad_count0
    INTEGER :: i1
    INTEGER :: ad_to48
    INTEGER :: ad_from27
    INTEGER :: ad_to49
    INTEGER :: ad_from28
    INTEGER :: ad_to50
    INTEGER :: ad_from29
    INTEGER :: ad_to51
    INTEGER :: ad_to52
    INTEGER :: ad_from30
    INTEGER :: ad_to53
    INTEGER :: ad_to54
    INTEGER :: ad_from31
    INTEGER :: ad_to55
    INTEGER :: ad_from32
    INTEGER :: ad_to56
    INTEGER :: ad_to57
    INTEGER :: ad_from33
    INTEGER :: ad_to58
    INTEGER :: ad_to59
    INTEGER :: ad_count1
    INTEGER :: i2

    te_2d = 0.0
    zsum0 = 0.0
    zsum1 = 0.0
    dpln = 0.0
    q2 = 0.0
    dp2 = 0.0
    pe1 = 0.0
    pe2 = 0.0
    pk1 = 0.0
    pk2 = 0.0
    pn2 = 0.0
    phis = 0.0
    pe0 = 0.0
    pe3 = 0.0
    gz = 0.0
    cvm = 0.0
    qv = 0.0
    rcp = 0.0
    rg = 0.0
    tmp = 0.0
    tpe = 0.0
    rrg = 0.0
    bkh = 0.0
    dtmp = 0.0
    k1k = 0.0
    dlnp = 0.0
    result1 = 0.0
    abs0 = 0.0
    nt = 0
    liq_wat = 0
    ice_wat = 0
    rainwat = 0
    snowwat = 0
    cld_amt = 0
    graupel = 0
    iq = 0
    n = 0
    kmp = 0
    kp = 0
    k_next = 0
    abs_kord_tm = 0
    abs_kord_tm_pert = 0
    iep1 = 0
    jep1 = 0
    iedp1 = 0
    jedp1 = 0
    arg1 = 0
    ad_count = 0
    ad_to = 0
    ad_to0 = 0
    ad_to1 = 0
    ad_to2 = 0
    ad_to3 = 0
    ad_to4 = 0
    ad_to5 = 0
    ad_to6 = 0
    ad_to7 = 0
    ad_to8 = 0
    ad_to9 = 0
    ad_to10 = 0
    ad_to11 = 0
    ad_to12 = 0
    ad_to13 = 0
    ad_to14 = 0
    ad_to15 = 0
    ad_to16 = 0
    ad_to17 = 0
    ad_to18 = 0
    ad_to19 = 0
    ad_to20 = 0
    ad_to21 = 0
    ad_to22 = 0
    ad_to23 = 0
    ad_to24 = 0
    ad_to25 = 0
    ad_to26 = 0
    ad_to27 = 0
    ad_to28 = 0
    ad_to29 = 0
    ad_to30 = 0
    ad_to31 = 0
    ad_to32 = 0
    ad_to33 = 0
    ad_from = 0
    ad_from0 = 0
    ad_from1 = 0
    ad_from2 = 0
    ad_from3 = 0
    ad_from4 = 0
    ad_from5 = 0
    ad_from6 = 0
    ad_from7 = 0
    ad_from8 = 0
    ad_from9 = 0
    ad_from10 = 0
    ad_from11 = 0
    ad_from12 = 0
    ad_from13 = 0
    ad_from14 = 0
    ad_from15 = 0
    ad_from16 = 0
    ad_from17 = 0
    ad_from18 = 0
    ad_from19 = 0
    ad_from20 = 0
    ad_from21 = 0
    ad_from22 = 0
    ad_from23 = 0
    ad_from24 = 0
    ad_from25 = 0
    ad_from26 = 0
    ad_count0 = 0
    ad_from27 = 0
    ad_from28 = 0
    ad_from29 = 0
    ad_from30 = 0
    ad_from31 = 0
    ad_from32 = 0
    ad_from33 = 0
    ad_count1 = 0
    branch = 0

    iep1 = ie + 1
    iedp1 = ied + 1
    jedp1 = jed + 1
    CALL POPCONTROL(3,branch)
    IF (branch .LT. 3) THEN
      IF (branch .EQ. 0) THEN
        CALL POPINTEGER(graupel)
        CALL POPREALARRAY(phis, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe0, (ie-is+2)*(km+1))
        CALL POPREALARRAY(pe1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe3, (ie-is+2)*(km+1))
        CALL POPREALARRAY(gz, ie - is + 1)
        CALL POPINTEGER(kmp)
        CALL POPREALARRAY(dp2, (ie-is+1)*km)
        CALL POPINTEGER(iep1)
        CALL POPREALARRAY(rg)
        CALL POPREALARRAY(tpe)
        CALL POPREALARRAY(k1k)
        CALL POPINTEGER(snowwat)
        CALL POPREALARRAY(dtmp)
        CALL POPINTEGER(cld_amt)
        CALL POPINTEGER(abs_kord_tm_pert)
        CALL POPREALARRAY(te_2d, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(rainwat)
        CALL POPREALARRAY(result1)
        CALL POPINTEGER(ice_wat)
        CALL POPREALARRAY(rrg)
        CALL POPREALARRAY(zsum0, (ie-is+1)*(je-js+1))
        CALL POPREALARRAY(pn2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(zsum1, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(liq_wat)
        CALL POPREALARRAY(pk1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(bkh)
        CALL POPREALARRAY(pk2, (ie-is+1)*(km+1))
        dtmp_ad = 0.0
        DO k=km,1,-1
          DO j=je,js,-1
            CALL POPCONTROL(1,branch)
            IF (branch .NE. 0) THEN
              DO i=ie,is,-1
                CALL POPREALARRAY(pt(i, j, k))
                temp36 = r_vir*q(i, j, k, sphum) + 1.
                temp_ad37 = pt_ad(i, j, k)/temp36
                dtmp_ad = dtmp_ad + pkz(i, j, k)*temp_ad37
                pkz_ad(i, j, k) = pkz_ad(i, j, k) + dtmp*temp_ad37
                q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) - (pt(i, j, &
&                 k)+dtmp*pkz(i, j, k))*r_vir*temp_ad37/temp36
                pt_ad(i, j, k) = temp_ad37
              END DO
            END IF
          END DO
        END DO
        gz_ad = 0.0
      ELSE IF (branch .EQ. 1) THEN
        CALL POPINTEGER(graupel)
        CALL POPREALARRAY(phis, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe0, (ie-is+2)*(km+1))
        CALL POPREALARRAY(pe1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe3, (ie-is+2)*(km+1))
        CALL POPREALARRAY(gz, ie - is + 1)
        CALL POPINTEGER(kmp)
        CALL POPREALARRAY(dp2, (ie-is+1)*km)
        CALL POPINTEGER(iep1)
        CALL POPREALARRAY(rg)
        CALL POPREALARRAY(tpe)
        CALL POPREALARRAY(k1k)
        CALL POPINTEGER(snowwat)
        CALL POPREALARRAY(dtmp)
        CALL POPINTEGER(cld_amt)
        CALL POPINTEGER(abs_kord_tm_pert)
        CALL POPREALARRAY(te_2d, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(rainwat)
        CALL POPREALARRAY(result1)
        CALL POPINTEGER(ice_wat)
        CALL POPREALARRAY(rrg)
        CALL POPREALARRAY(zsum0, (ie-is+1)*(je-js+1))
        CALL POPREALARRAY(pn2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(zsum1, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(liq_wat)
        CALL POPREALARRAY(pk1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(bkh)
        CALL POPREALARRAY(pk2, (ie-is+1)*(km+1))
        dtmp_ad = 0.0
        DO k=km,1,-1
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(pt(i, j, k))
              temp39 = r_vir*q(i, j, k, sphum) + 1.
              temp_ad38 = pt_ad(i, j, k)/temp39
              temp38 = pkz(i, j, k)
              temp37 = pt(i, j, k) + dtmp
              dtmp_ad = dtmp_ad + temp38*temp_ad38
              pkz_ad(i, j, k) = pkz_ad(i, j, k) + temp37*temp_ad38
              q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) - temp37*&
&               temp38*r_vir*temp_ad38/temp39
              pt_ad(i, j, k) = temp38*temp_ad38
            END DO
          END DO
        END DO
        gz_ad = 0.0
      ELSE
        CALL POPINTEGER(graupel)
        CALL POPREALARRAY(phis, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe0, (ie-is+2)*(km+1))
        CALL POPREALARRAY(pe1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe3, (ie-is+2)*(km+1))
        CALL POPREALARRAY(gz, ie - is + 1)
        CALL POPINTEGER(kmp)
        CALL POPREALARRAY(dp2, (ie-is+1)*km)
        CALL POPINTEGER(iep1)
        CALL POPREALARRAY(rg)
        CALL POPREALARRAY(tpe)
        CALL POPREALARRAY(k1k)
        CALL POPINTEGER(snowwat)
        CALL POPREALARRAY(dtmp)
        CALL POPREALARRAY(tmp)
        CALL POPINTEGER(cld_amt)
        CALL POPINTEGER(abs_kord_tm_pert)
        CALL POPREALARRAY(te_2d, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(rainwat)
        CALL POPREALARRAY(result1)
        CALL POPINTEGER(ice_wat)
        CALL POPREALARRAY(rrg)
        CALL POPREALARRAY(zsum0, (ie-is+1)*(je-js+1))
        CALL POPREALARRAY(pn2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(zsum1, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(liq_wat)
        CALL POPREALARRAY(pk1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(bkh)
        CALL POPREALARRAY(pk2, (ie-is+1)*(km+1))
        gz_ad = 0.0
        dtmp_ad = 0.0
        DO j=je,js,-1
          DO k=1,km,1
            DO i=ie,is,-1
              dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
              CALL POPREALARRAY(gz(i))
              temp_ad39 = (r_vir*q(i, j, k, sphum)+1.)*gz_ad(i)
              tmp_ad = pt_ad(i, j, k) + dlnp*temp_ad39
              CALL POPREALARRAY(pt(i, j, k))
              temp45 = r_vir*q(i, j, k, sphum) + 1.
              temp_ad41 = pt_ad(i, j, k)/temp45
              dtmp_ad = dtmp_ad + pkz(i, j, k)*temp_ad41
              pkz_ad(i, j, k) = pkz_ad(i, j, k) + dtmp*temp_ad41
              pt_ad(i, j, k) = 0.0
              temp44 = r_vir*q(i, j, k, sphum) + 1.
              temp43 = delp(i, j, k)
              temp42 = pe(i, k, j)
              temp41 = temp42*dlnp/temp43
              temp40 = (cp-temp41)*temp44
              temp_ad42 = -(tpe*tmp_ad/temp40**2)
              q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) + (cp-temp41)*&
&               r_vir*temp_ad42 - dtmp*pkz(i, j, k)*r_vir*temp_ad41/&
&               temp45 + dlnp*tmp*r_vir*gz_ad(i)
              temp_ad40 = -(temp44*temp_ad42/temp43)
              dlnp_ad = temp42*temp_ad40 + tmp*temp_ad39
              CALL POPREALARRAY(tmp)
              tpe_ad = tmp_ad/temp40
              pe_ad(i, k, j) = pe_ad(i, k, j) + dlnp*temp_ad40
              delp_ad(i, j, k) = delp_ad(i, j, k) - temp41*temp_ad40
              peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + rg*dlnp_ad
              peln_ad(i, k, j) = peln_ad(i, k, j) - rg*dlnp_ad
              CALL POPREALARRAY(tpe)
              temp_ad43 = -(gridstruct%rsin2(i, j)*0.25*tpe_ad)
              temp_ad44 = -(gridstruct%cosa_s(i, j)*temp_ad43)
              temp_ad45 = (v(i, j, k)+v(i+1, j, k))*temp_ad44
              temp_ad46 = (u(i, j, k)+u(i, j+1, k))*temp_ad44
              te_ad(i, j, k) = te_ad(i, j, k) + tpe_ad
              gz_ad(i) = gz_ad(i) - tpe_ad
              u_ad(i, j, k) = u_ad(i, j, k) + temp_ad45 + 2*u(i, j, k)*&
&               temp_ad43
              u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp_ad45 + 2*u(i, j+1&
&               , k)*temp_ad43
              v_ad(i, j, k) = v_ad(i, j, k) + temp_ad46 + 2*v(i, j, k)*&
&               temp_ad43
              v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp_ad46 + 2*v(i+1, j&
&               , k)*temp_ad43
            END DO
          END DO
          DO i=ie,is,-1
            CALL POPREALARRAY(gz(i))
            gz_ad(i) = 0.0
          END DO
        END DO
      END IF
    ELSE IF (branch .LT. 5) THEN
      IF (branch .EQ. 3) THEN
        CALL POPINTEGER(graupel)
        CALL POPREALARRAY(phis, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe0, (ie-is+2)*(km+1))
        CALL POPREALARRAY(pe1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe3, (ie-is+2)*(km+1))
        CALL POPREALARRAY(gz, ie - is + 1)
        CALL POPINTEGER(kmp)
        CALL POPREALARRAY(dp2, (ie-is+1)*km)
        CALL POPINTEGER(iep1)
        CALL POPREALARRAY(rg)
        CALL POPREALARRAY(tpe)
        CALL POPREALARRAY(k1k)
        CALL POPINTEGER(snowwat)
        CALL POPINTEGER(cld_amt)
        CALL POPINTEGER(abs_kord_tm_pert)
        CALL POPREALARRAY(te_2d, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(rainwat)
        CALL POPREALARRAY(result1)
        CALL POPINTEGER(ice_wat)
        CALL POPREALARRAY(rrg)
        CALL POPREALARRAY(zsum0, (ie-is+1)*(je-js+1))
        CALL POPREALARRAY(pn2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(zsum1, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(liq_wat)
        CALL POPREALARRAY(pk1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(bkh)
        CALL POPREALARRAY(pk2, (ie-is+1)*(km+1))
        gz_ad = 0.0
        dtmp_ad = 0.0
      ELSE
        CALL POPINTEGER(graupel)
        CALL POPREALARRAY(phis, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe0, (ie-is+2)*(km+1))
        CALL POPREALARRAY(pe1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(pe3, (ie-is+2)*(km+1))
        CALL POPREALARRAY(gz, ie - is + 1)
        CALL POPINTEGER(kmp)
        CALL POPREALARRAY(dp2, (ie-is+1)*km)
        CALL POPINTEGER(iep1)
        CALL POPREALARRAY(rg)
        CALL POPREALARRAY(tpe)
        CALL POPREALARRAY(k1k)
        CALL POPINTEGER(snowwat)
        CALL POPINTEGER(cld_amt)
        CALL POPINTEGER(abs_kord_tm_pert)
        CALL POPREALARRAY(te_2d, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(rainwat)
        CALL POPREALARRAY(result1)
        CALL POPINTEGER(ice_wat)
        CALL POPREALARRAY(rrg)
        CALL POPREALARRAY(zsum0, (ie-is+1)*(je-js+1))
        CALL POPREALARRAY(pn2, (ie-is+1)*(km+1))
        CALL POPREALARRAY(zsum1, (ie-is+1)*(je-js+1))
        CALL POPINTEGER(liq_wat)
        CALL POPREALARRAY(pk1, (ie-is+1)*(km+1))
        CALL POPREALARRAY(bkh)
        CALL POPREALARRAY(pk2, (ie-is+1)*(km+1))
        DO k=km,1,-1
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(pt(i, j, k))
              temp_ad47 = pt_ad(i, j, k)/pkz(i, j, k)
              pkz_ad(i, j, k) = pkz_ad(i, j, k) - pt(i, j, k)*temp_ad47/&
&               pkz(i, j, k)
              pt_ad(i, j, k) = temp_ad47
            END DO
          END DO
        END DO
        gz_ad = 0.0
        dtmp_ad = 0.0
      END IF
    ELSE IF (branch .EQ. 5) THEN
      CALL POPINTEGER(graupel)
      CALL POPREALARRAY(phis, (ie-is+1)*(km+1))
      CALL POPREALARRAY(pe0, (ie-is+2)*(km+1))
      CALL POPREALARRAY(pe1, (ie-is+1)*(km+1))
      CALL POPREALARRAY(pe2, (ie-is+1)*(km+1))
      CALL POPREALARRAY(pe3, (ie-is+2)*(km+1))
      CALL POPREALARRAY(gz, ie - is + 1)
      CALL POPINTEGER(kmp)
      CALL POPREALARRAY(dp2, (ie-is+1)*km)
      CALL POPINTEGER(iep1)
      CALL POPREALARRAY(rg)
      CALL POPREALARRAY(tpe)
      CALL POPREALARRAY(k1k)
      CALL POPINTEGER(snowwat)
      CALL POPINTEGER(cld_amt)
      CALL POPINTEGER(abs_kord_tm_pert)
      CALL POPREALARRAY(te_2d, (ie-is+1)*(je-js+1))
      CALL POPINTEGER(rainwat)
      CALL POPREALARRAY(result1)
      CALL POPINTEGER(ice_wat)
      CALL POPREALARRAY(rrg)
      CALL POPREALARRAY(zsum0, (ie-is+1)*(je-js+1))
      CALL POPREALARRAY(pn2, (ie-is+1)*(km+1))
      CALL POPREALARRAY(zsum1, (ie-is+1)*(je-js+1))
      CALL POPINTEGER(liq_wat)
      CALL POPREALARRAY(pk1, (ie-is+1)*(km+1))
      CALL POPREALARRAY(bkh)
      CALL POPREALARRAY(pk2, (ie-is+1)*(km+1))
      gz_ad = 0.0
      dtmp_ad = 0.0
      DO j=je,js,-1
        DO k=1,km,1
          DO i=ie,is,-1
            temp_ad49 = pt_ad(i, j, k)/pkz(i, j, k)
            dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
            tmp = tpe/(cp-pe(i, k, j)*dlnp/delp(i, j, k))
            CALL POPREALARRAY(gz(i))
            tmp_ad = temp_ad49 + dlnp*gz_ad(i)
            CALL POPREALARRAY(pt(i, j, k))
            pkz_ad(i, j, k) = pkz_ad(i, j, k) - tmp*temp_ad49/pkz(i, j, &
&             k)
            dtmp_ad = dtmp_ad + pt_ad(i, j, k)
            pt_ad(i, j, k) = 0.0
            temp48 = delp(i, j, k)
            temp47 = pe(i, k, j)
            temp46 = temp47*dlnp/temp48
            temp_ad50 = tmp_ad/(cp-temp46)
            temp_ad48 = tpe*temp_ad50/((cp-temp46)*temp48)
            dlnp_ad = temp47*temp_ad48 + tmp*gz_ad(i)
            tpe_ad = temp_ad50
            pe_ad(i, k, j) = pe_ad(i, k, j) + dlnp*temp_ad48
            delp_ad(i, j, k) = delp_ad(i, j, k) - temp46*temp_ad48
            peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + rg*dlnp_ad
            peln_ad(i, k, j) = peln_ad(i, k, j) - rg*dlnp_ad
            CALL POPREALARRAY(tpe)
            temp_ad51 = -(gridstruct%rsin2(i, j)*0.25*tpe_ad)
            temp_ad52 = -(gridstruct%cosa_s(i, j)*temp_ad51)
            temp_ad53 = (v(i, j, k)+v(i+1, j, k))*temp_ad52
            temp_ad54 = (u(i, j, k)+u(i, j+1, k))*temp_ad52
            te_ad(i, j, k) = te_ad(i, j, k) + tpe_ad
            gz_ad(i) = gz_ad(i) - tpe_ad
            u_ad(i, j, k) = u_ad(i, j, k) + temp_ad53 + 2*u(i, j, k)*&
&             temp_ad51
            u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp_ad53 + 2*u(i, j+1, &
&             k)*temp_ad51
            v_ad(i, j, k) = v_ad(i, j, k) + temp_ad54 + 2*v(i, j, k)*&
&             temp_ad51
            v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp_ad54 + 2*v(i+1, j, &
&             k)*temp_ad51
          END DO
        END DO
        DO i=ie,is,-1
          CALL POPREALARRAY(gz(i))
          gz_ad(i) = 0.0
        END DO
      END DO
    ELSE
      CALL POPINTEGER(graupel)
      CALL POPREALARRAY(phis, (ie-is+1)*(km+1))
      CALL POPREALARRAY(pe0, (ie-is+2)*(km+1))
      CALL POPREALARRAY(pe1, (ie-is+1)*(km+1))
      CALL POPREALARRAY(pe2, (ie-is+1)*(km+1))
      CALL POPREALARRAY(pe3, (ie-is+2)*(km+1))
      CALL POPREALARRAY(gz, ie - is + 1)
      CALL POPINTEGER(kmp)
      CALL POPREALARRAY(dp2, (ie-is+1)*km)
      CALL POPINTEGER(iep1)
      CALL POPREALARRAY(rg)
      CALL POPREALARRAY(tpe)
      CALL POPREALARRAY(k1k)
      CALL POPINTEGER(snowwat)
      CALL POPINTEGER(cld_amt)
      CALL POPINTEGER(abs_kord_tm_pert)
      CALL POPREALARRAY(te_2d, (ie-is+1)*(je-js+1))
      CALL POPINTEGER(rainwat)
      CALL POPREALARRAY(result1)
      CALL POPINTEGER(ice_wat)
      CALL POPREALARRAY(rrg)
      CALL POPREALARRAY(zsum0, (ie-is+1)*(je-js+1))
      CALL POPREALARRAY(pn2, (ie-is+1)*(km+1))
      CALL POPREALARRAY(zsum1, (ie-is+1)*(je-js+1))
      CALL POPINTEGER(liq_wat)
      CALL POPREALARRAY(pk1, (ie-is+1)*(km+1))
      CALL POPREALARRAY(bkh)
      CALL POPREALARRAY(pk2, (ie-is+1)*(km+1))
      gz_ad = 0.0
      dtmp_ad = 0.0
    END IF
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        DO i=ie,is,-1
          DO k=km,kmp,-1
            te_ad(i, j, k) = te_ad(i, j, k) + te0_2d_ad(i, j)
          END DO
        END DO
      END DO
    ELSE IF (branch .NE. 1) THEN
      GOTO 100
    END IF
    rrg = -(rdgas/grav)
    DO k=km,kmp,-1
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(pkz(i, j, k))
            temp35 = delz(i, j, k)
            temp34 = delp(i, j, k)*pt(i, j, k)
            temp33 = temp34/temp35
            temp_ad36 = akap*EXP(akap*LOG(rrg*temp33))*pkz_ad(i, j, k)/(&
&             temp33*temp35)
            delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*temp_ad36
            pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*temp_ad36
            delz_ad(i, j, k) = delz_ad(i, j, k) - temp33*temp_ad36
            pkz_ad(i, j, k) = 0.0
          END DO
        END DO
      END IF
    END DO
 100 CALL POPCONTROL(3,branch)
    IF (branch .LT. 3) THEN
      IF (branch .EQ. 0) THEN
        temp_ad34 = dtmp_ad/(cp*result1)
        tpe_ad = temp_ad34
        result1_ad = -(tpe*temp_ad34/result1)
        CALL G_SUM_ADM(domain, zsum0, zsum0_ad, is, ie, js, je, ng, &
&                gridstruct%area_64, 0, reproduce=.true., g_sum_ad=&
&                result1_ad)
        zsum1_ad = 0.0
      ELSE IF (branch .EQ. 1) THEN
        temp_ad35 = dtmp_ad/(cv_air*result1)
        tpe_ad = temp_ad35
        result1_ad = -(tpe*temp_ad35/result1)
        CALL G_SUM_ADM(domain, zsum1, zsum1_ad, is, ie, js, je, ng, &
&                gridstruct%area_64, 0, reproduce=.true., g_sum_ad=&
&                result1_ad)
        zsum0_ad = 0.0
      ELSE
        result1_ad = -(e_flux*4.*pi*grav*pdt*radius**2*dtmp_ad/(cp*&
&         result1**2))
        CALL G_SUM_ADM(domain, zsum0, zsum0_ad, is, ie, js, je, ng, &
&                gridstruct%area_64, 0, reproduce=.true., g_sum_ad=&
&                result1_ad)
        zsum1_ad = 0.0
        GOTO 110
      END IF
      result1_ad = consv*tpe_ad
      CALL G_SUM_ADM(domain, te_2d, te_2d_ad, is, ie, js, je, ng, &
&              gridstruct%area_64, 0, reproduce=.true., g_sum_ad=&
&              result1_ad)
      phis_ad = 0.0
      DO j=je,js,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          DO i=ie,is,-1
            pk_ad(i, j, 1) = pk_ad(i, j, 1) + ptop*zsum0_ad(i, j)
            pk_ad(i, j, km+1) = pk_ad(i, j, km+1) - ptop*zsum0_ad(i, j)
            zsum1_ad(i, j) = zsum1_ad(i, j) + zsum0_ad(i, j)
            zsum0_ad(i, j) = 0.0
          END DO
        END IF
        DO k=km,2,-1
          DO i=ie,is,-1
            pkz_ad(i, j, k) = pkz_ad(i, j, k) + delp(i, j, k)*zsum1_ad(i&
&             , j)
            delp_ad(i, j, k) = delp_ad(i, j, k) + pkz(i, j, k)*zsum1_ad(&
&             i, j)
          END DO
        END DO
        DO i=ie,is,-1
          pkz_ad(i, j, 1) = pkz_ad(i, j, 1) + delp(i, j, 1)*zsum1_ad(i, &
&           j)
          delp_ad(i, j, 1) = delp_ad(i, j, 1) + pkz(i, j, 1)*zsum1_ad(i&
&           , j)
          zsum1_ad(i, j) = 0.0
          te0_2d_ad(i, j) = te0_2d_ad(i, j) + te_2d_ad(i, j)
          te_2d_ad(i, j) = -te_2d_ad(i, j)
        END DO
        CALL POPCONTROL(3,branch)
        IF (branch .LT. 3) THEN
          IF (branch .NE. 0) THEN
            IF (branch .EQ. 1) THEN
              DO k=km,2,-1
                DO i=ie,is,-1
                  te_ad(i, j, k) = te_ad(i, j, k) + delp(i, j, k)*&
&                   te_2d_ad(i, j)
                  delp_ad(i, j, k) = delp_ad(i, j, k) + te(i, j, k)*&
&                   te_2d_ad(i, j)
                END DO
              END DO
              DO i=ie,is,-1
                te_ad(i, j, 1) = te_ad(i, j, 1) + delp(i, j, 1)*te_2d_ad&
&                 (i, j)
                delp_ad(i, j, 1) = delp_ad(i, j, 1) + te(i, j, 1)*&
&                 te_2d_ad(i, j)
                te_2d_ad(i, j) = 0.0
              END DO
            ELSE
              DO k=km,1,-1
                DO i=ie,is,-1
                  temp32 = v(i, j, k) + v(i+1, j, k)
                  temp31 = u(i, j, k) + u(i, j+1, k)
                  temp30 = 0.5*gridstruct%rsin2(i, j)
                  temp29 = r_vir*q(i, j, k, sphum) + 1.
                  temp28 = pt(i, j, k)/temp29
                  temp_ad29 = delp(i, j, k)*te_2d_ad(i, j)
                  temp_ad30 = cv_air*temp_ad29/temp29
                  temp_ad31 = 0.5*temp_ad29
                  temp_ad32 = temp30*temp_ad31
                  temp_ad33 = -(gridstruct%cosa_s(i, j)*temp_ad32)
                  delp_ad(i, j, k) = delp_ad(i, j, k) + (cv_air*temp28+&
&                   0.5*(phis(i, k)+phis(i, k+1)+w(i, j, k)**2+temp30*(u&
&                   (i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j&
&                   , k)**2-gridstruct%cosa_s(i, j)*(temp31*temp32))))*&
&                   te_2d_ad(i, j)
                  pt_ad(i, j, k) = pt_ad(i, j, k) + temp_ad30
                  q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) - temp28*&
&                   r_vir*temp_ad30
                  phis_ad(i, k) = phis_ad(i, k) + temp_ad31
                  phis_ad(i, k+1) = phis_ad(i, k+1) + temp_ad31
                  w_ad(i, j, k) = w_ad(i, j, k) + 2*w(i, j, k)*temp_ad31
                  u_ad(i, j, k) = u_ad(i, j, k) + temp32*temp_ad33 + 2*u&
&                   (i, j, k)*temp_ad32
                  u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp32*temp_ad33 +&
&                   2*u(i, j+1, k)*temp_ad32
                  v_ad(i, j, k) = v_ad(i, j, k) + temp31*temp_ad33 + 2*v&
&                   (i, j, k)*temp_ad32
                  v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp31*temp_ad33 +&
&                   2*v(i+1, j, k)*temp_ad32
                END DO
              END DO
              DO i=ie,is,-1
                te_2d_ad(i, j) = 0.0
              END DO
              DO i=ie,is,-1
                DO k=1,km,1
                  CALL POPREALARRAY(phis(i, k))
                  phis_ad(i, k+1) = phis_ad(i, k+1) + phis_ad(i, k)
                  delz_ad(i, j, k) = delz_ad(i, j, k) - grav*phis_ad(i, &
&                   k)
                  phis_ad(i, k) = 0.0
                END DO
                CALL POPREALARRAY(phis(i, km+1))
                phis_ad(i, km+1) = 0.0
              END DO
            END IF
          END IF
        ELSE IF (branch .EQ. 3) THEN
          DO k=km,1,-1
            DO i=ie,is,-1
              temp27 = v(i, j, k) + v(i+1, j, k)
              temp26 = u(i, j, k) + u(i, j+1, k)
              temp25 = 0.25*gridstruct%rsin2(i, j)
              temp_ad26 = delp(i, j, k)*te_2d_ad(i, j)
              temp_ad27 = temp25*temp_ad26
              temp_ad28 = -(gridstruct%cosa_s(i, j)*temp_ad27)
              delp_ad(i, j, k) = delp_ad(i, j, k) + (cp_air*(pt(i, j, k)&
&               *pkz(i, j, k))+temp25*(u(i, j, k)**2+u(i, j+1, k)**2+v(i&
&               , j, k)**2+v(i+1, j, k)**2-gridstruct%cosa_s(i, j)*(&
&               temp26*temp27)))*te_2d_ad(i, j)
              pt_ad(i, j, k) = pt_ad(i, j, k) + cp_air*pkz(i, j, k)*&
&               temp_ad26
              pkz_ad(i, j, k) = pkz_ad(i, j, k) + cp_air*pt(i, j, k)*&
&               temp_ad26
              u_ad(i, j, k) = u_ad(i, j, k) + temp27*temp_ad28 + 2*u(i, &
&               j, k)*temp_ad27
              u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp27*temp_ad28 + 2*u&
&               (i, j+1, k)*temp_ad27
              v_ad(i, j, k) = v_ad(i, j, k) + temp26*temp_ad28 + 2*v(i, &
&               j, k)*temp_ad27
              v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp26*temp_ad28 + 2*v&
&               (i+1, j, k)*temp_ad27
            END DO
          END DO
          DO i=ie,is,-1
            pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + hs(i, j)*te_2d_ad(i&
&             , j)
            pe_ad(i, 1, j) = pe_ad(i, 1, j) - gz(i)*te_2d_ad(i, j)
            gz_ad(i) = gz_ad(i) - pe(i, 1, j)*te_2d_ad(i, j)
            te_2d_ad(i, j) = 0.0
          END DO
          DO i=ie,is,-1
            DO k=km,1,-1
              temp_ad25 = cp_air*pt(i, j, k)*gz_ad(i)
              pt_ad(i, j, k) = pt_ad(i, j, k) + cp_air*(pk(i, j, k+1)-pk&
&               (i, j, k))*gz_ad(i)
              pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp_ad25
              pk_ad(i, j, k) = pk_ad(i, j, k) - temp_ad25
            END DO
            CALL POPREALARRAY(gz(i))
            gz_ad(i) = 0.0
          END DO
        ELSE IF (branch .EQ. 4) THEN
          DO k=km,1,-1
            DO i=ie,is,-1
              temp24 = v(i, j, k) + v(i+1, j, k)
              temp23 = u(i, j, k) + u(i, j+1, k)
              temp22 = 0.5*gridstruct%rsin2(i, j)
              temp21 = r_vir*q(i, j, k, sphum) + 1.
              temp20 = pt(i, j, k)/temp21
              temp_ad20 = delp(i, j, k)*te_2d_ad(i, j)
              temp_ad21 = cv_air*temp_ad20/temp21
              temp_ad22 = 0.5*temp_ad20
              temp_ad23 = temp22*temp_ad22
              temp_ad24 = -(gridstruct%cosa_s(i, j)*temp_ad23)
              delp_ad(i, j, k) = delp_ad(i, j, k) + (cv_air*temp20+0.5*(&
&               phis(i, k)+phis(i, k+1)+w(i, j, k)**2+temp22*(u(i, j, k)&
&               **2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-&
&               gridstruct%cosa_s(i, j)*(temp23*temp24))))*te_2d_ad(i, j&
&               )
              pt_ad(i, j, k) = pt_ad(i, j, k) + temp_ad21
              q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) - temp20*r_vir&
&               *temp_ad21
              phis_ad(i, k) = phis_ad(i, k) + temp_ad22
              phis_ad(i, k+1) = phis_ad(i, k+1) + temp_ad22
              w_ad(i, j, k) = w_ad(i, j, k) + 2*w(i, j, k)*temp_ad22
              u_ad(i, j, k) = u_ad(i, j, k) + temp24*temp_ad24 + 2*u(i, &
&               j, k)*temp_ad23
              u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp24*temp_ad24 + 2*u&
&               (i, j+1, k)*temp_ad23
              v_ad(i, j, k) = v_ad(i, j, k) + temp23*temp_ad24 + 2*v(i, &
&               j, k)*temp_ad23
              v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp23*temp_ad24 + 2*v&
&               (i+1, j, k)*temp_ad23
            END DO
          END DO
          DO k=1,km,1
            DO i=ie,is,-1
              CALL POPREALARRAY(phis(i, k))
              phis_ad(i, k+1) = phis_ad(i, k+1) + phis_ad(i, k)
              delz_ad(i, j, k) = delz_ad(i, j, k) - grav*phis_ad(i, k)
              phis_ad(i, k) = 0.0
            END DO
          END DO
          DO i=ie,is,-1
            CALL POPREALARRAY(phis(i, km+1))
            phis_ad(i, km+1) = 0.0
            te_2d_ad(i, j) = 0.0
          END DO
        ELSE
          DO k=km,1,-1
            DO i=ie,is,-1
              temp19 = v(i, j, k) + v(i+1, j, k)
              temp18 = u(i, j, k) + u(i, j+1, k)
              temp17 = 0.25*gridstruct%rsin2(i, j)
              temp_ad17 = delp(i, j, k)*te_2d_ad(i, j)
              temp_ad18 = temp17*temp_ad17
              temp_ad19 = -(gridstruct%cosa_s(i, j)*temp_ad18)
              delp_ad(i, j, k) = delp_ad(i, j, k) + (cp*pt(i, j, k)+&
&               temp17*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+&
&               1, j, k)**2-gridstruct%cosa_s(i, j)*(temp18*temp19)))*&
&               te_2d_ad(i, j)
              pt_ad(i, j, k) = pt_ad(i, j, k) + cp*temp_ad17
              u_ad(i, j, k) = u_ad(i, j, k) + temp19*temp_ad19 + 2*u(i, &
&               j, k)*temp_ad18
              u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp19*temp_ad19 + 2*u&
&               (i, j+1, k)*temp_ad18
              v_ad(i, j, k) = v_ad(i, j, k) + temp18*temp_ad19 + 2*v(i, &
&               j, k)*temp_ad18
              v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp18*temp_ad19 + 2*v&
&               (i+1, j, k)*temp_ad18
            END DO
          END DO
          DO i=ie,is,-1
            pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + hs(i, j)*te_2d_ad(i&
&             , j)
            pe_ad(i, 1, j) = pe_ad(i, 1, j) - gz(i)*te_2d_ad(i, j)
            gz_ad(i) = gz_ad(i) - pe(i, 1, j)*te_2d_ad(i, j)
            te_2d_ad(i, j) = 0.0
          END DO
          DO i=ie,is,-1
            DO k=km,1,-1
              temp_ad16 = rg*pt(i, j, k)*gz_ad(i)
              pt_ad(i, j, k) = pt_ad(i, j, k) + rg*(peln(i, k+1, j)-peln&
&               (i, k, j))*gz_ad(i)
              peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad16
              peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad16
            END DO
            CALL POPREALARRAY(gz(i))
            gz_ad(i) = 0.0
          END DO
        END IF
      END DO
      GOTO 130
    ELSE IF (branch .EQ. 3) THEN
      result1_ad = -(e_flux*4.*pi*grav*pdt*radius**2*dtmp_ad/(cv_air*&
&       result1**2))
      CALL G_SUM_ADM(domain, zsum1, zsum1_ad, is, ie, js, je, ng, &
&              gridstruct%area_64, 0, reproduce=.true., g_sum_ad=&
&              result1_ad)
      zsum0_ad = 0.0
    ELSE IF (branch .EQ. 4) THEN
      GOTO 120
    ELSE
      phis_ad = 0.0
      GOTO 130
    END IF
 110 CALL POPREALARRAY(e_flux)
    DO j=je,js,-1
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        DO i=ie,is,-1
          pk_ad(i, j, 1) = pk_ad(i, j, 1) + ptop*zsum0_ad(i, j)
          pk_ad(i, j, km+1) = pk_ad(i, j, km+1) - ptop*zsum0_ad(i, j)
          zsum1_ad(i, j) = zsum1_ad(i, j) + zsum0_ad(i, j)
          zsum0_ad(i, j) = 0.0
        END DO
      END IF
      DO k=km,2,-1
        DO i=ie,is,-1
          pkz_ad(i, j, k) = pkz_ad(i, j, k) + delp(i, j, k)*zsum1_ad(i, &
&           j)
          delp_ad(i, j, k) = delp_ad(i, j, k) + pkz(i, j, k)*zsum1_ad(i&
&           , j)
        END DO
      END DO
      DO i=ie,is,-1
        pkz_ad(i, j, 1) = pkz_ad(i, j, 1) + delp(i, j, 1)*zsum1_ad(i, j)
        delp_ad(i, j, 1) = delp_ad(i, j, 1) + pkz(i, j, 1)*zsum1_ad(i, j&
&         )
        zsum1_ad(i, j) = 0.0
      END DO
    END DO
 120 phis_ad = 0.0
 130 DO k=km,2,-1
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(pe(i, k, j))
          ua_ad(i, j, k-1) = ua_ad(i, j, k-1) + pe_ad(i, k, j)
          pe_ad(i, k, j) = 0.0
        END DO
      END DO
    END DO
    CALL POPINTEGER(ad_count1)
    DO i2=1,ad_count1
      IF (i2 .EQ. 1) THEN
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          pe0_ad = 0.0
          pe1_ad = 0.0
          pe2_ad = 0.0
          pe3_ad = 0.0
          dp2_ad = 0.0
          q2_ad = 0.0
          pn2_ad = 0.0
          pk1_ad = 0.0
          pk2_ad = 0.0
          GOTO 150
        ELSE
          ws_ad = 0.0
          peln_ad = 0.0
          q_ad = 0.0
          u_ad = 0.0
          v_ad = 0.0
          w_ad = 0.0
          delp_ad = 0.0
          ua_ad = 0.0
          delz_ad = 0.0
          omga_ad = 0.0
          te0_2d_ad = 0.0
          pkz_ad = 0.0
          pe_ad = 0.0
          pk_ad = 0.0
          ps_ad = 0.0
          pt_ad = 0.0
          te_ad = 0.0
          phis_ad = 0.0
          pe0_ad = 0.0
          pe1_ad = 0.0
          pe2_ad = 0.0
          pe3_ad = 0.0
          dp2_ad = 0.0
          q2_ad = 0.0
          pn2_ad = 0.0
          pk1_ad = 0.0
          pk2_ad = 0.0
        END IF
      ELSE
        CALL POPINTEGER(ad_to59)
        DO k=ad_to59,1,-1
          CALL POPINTEGER(ad_from33)
          CALL POPINTEGER(ad_to58)
          DO i=ad_to58,ad_from33,-1
            CALL POPREALARRAY(ua(i, j, k))
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + ua_ad(i, j, k)
            ua_ad(i, j, k) = 0.0
          END DO
        END DO
        CALL POPCONTROL(2,branch)
        IF (branch .NE. 0) THEN
          IF (branch .EQ. 1) THEN
            gz_ad = 0.0
            CALL MAP1_PPM_BWD(km, pe0, pe0_ad, gz, gz_ad, km, pe3, &
&                          pe3_ad, v, v_ad, is, iep1, j, isd, iedp1, jsd&
&                          , jed, -1, kord_mt)
          ELSE
            CALL POPREALARRAY(v, (ied-isd+2)*(jed-jsd+1)*km)
            gz_ad = 0.0
            CALL MAP1_PPM_ADM(km, pe0, pe0_ad, gz, gz_ad, km, pe3, &
&                       pe3_ad, v, v_ad, is, iep1, j, isd, iedp1, jsd, &
&                       jed, -1, kord_mt_pert)
          END IF
          CALL POPINTEGER(ad_to57)
          DO k=ad_to57,2,-1
            CALL POPINTEGER(ad_from32)
            CALL POPINTEGER(ad_to56)
            DO i=ad_to56,ad_from32,-1
              CALL POPREALARRAY(pe3(i, k))
              pe_ad(i-1, km+1, j) = pe_ad(i-1, km+1, j) + bkh*pe3_ad(i, &
&               k)
              pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + bkh*pe3_ad(i, k)
              pe3_ad(i, k) = 0.0
              CALL POPREALARRAY(pe0(i, k))
              pe_ad(i-1, k, j) = pe_ad(i-1, k, j) + 0.5*pe0_ad(i, k)
              pe_ad(i, k, j) = pe_ad(i, k, j) + 0.5*pe0_ad(i, k)
              pe0_ad(i, k) = 0.0
            END DO
            CALL POPREALARRAY(bkh)
          END DO
          CALL POPINTEGER(ad_from31)
          CALL POPINTEGER(ad_to55)
          DO i=ad_to55,ad_from31,-1
            CALL POPREALARRAY(pe3(i, 1))
            pe3_ad(i, 1) = 0.0
          END DO
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          gz_ad = 0.0
          CALL MAP1_PPM_BWD(km, pe0(is:ie, :), pe0_ad(is:ie, :), gz, &
&                        gz_ad, km, pe3(is:ie, :), pe3_ad(is:ie, :), u, &
&                        u_ad, is, ie, j, isd, ied, jsd, jedp1, -1, &
&                        kord_mt)
        ELSE
          CALL POPREALARRAY(u, (ied-isd+1)*(jed-jsd+2)*km)
          gz_ad = 0.0
          CALL MAP1_PPM_ADM(km, pe0(is:ie, :), pe0_ad(is:ie, :), gz, &
&                     gz_ad, km, pe3(is:ie, :), pe3_ad(is:ie, :), u, &
&                     u_ad, is, ie, j, isd, ied, jsd, jedp1, -1, &
&                     kord_mt_pert)
        END IF
        CALL POPINTEGER(ad_to54)
        DO k=ad_to54,1,-1
          CALL POPINTEGER(ad_from30)
          CALL POPINTEGER(ad_to53)
          DO i=ad_to53,ad_from30,-1
            CALL POPREALARRAY(pe3(i, k))
            pe_ad(i, km+1, j-1) = pe_ad(i, km+1, j-1) + bkh*pe3_ad(i, k)
            pe1_ad(i, km+1) = pe1_ad(i, km+1) + bkh*pe3_ad(i, k)
            pe3_ad(i, k) = 0.0
          END DO
          CALL POPREALARRAY(bkh)
        END DO
        CALL POPINTEGER(ad_to52)
        DO k=ad_to52,2,-1
          CALL POPINTEGER(ad_from29)
          CALL POPINTEGER(ad_to51)
          DO i=ad_to51,ad_from29,-1
            CALL POPREALARRAY(pe0(i, k))
            pe_ad(i, k, j-1) = pe_ad(i, k, j-1) + 0.5*pe0_ad(i, k)
            pe1_ad(i, k) = pe1_ad(i, k) + 0.5*pe0_ad(i, k)
            pe0_ad(i, k) = 0.0
          END DO
        END DO
        CALL POPINTEGER(k)
        CALL POPINTEGER(ad_from28)
        CALL POPINTEGER(ad_to50)
        DO i=ad_to50,ad_from28,-1
          CALL POPREALARRAY(pe0(i, 1))
          pe_ad(i, 1, j) = pe_ad(i, 1, j) + pe0_ad(i, 1)
          pe0_ad(i, 1) = 0.0
        END DO
        CALL POPCONTROL(2,branch)
        IF (branch .EQ. 0) THEN
          GOTO 140
        ELSE
          IF (branch .NE. 1) THEN
            CALL POPINTEGER(ad_from27)
            CALL POPINTEGER(ad_to49)
            DO i=ad_to49,ad_from27,-1
              CALL POPINTEGER(ad_to48)
              DO n=ad_to48,1,-1
                CALL POPCONTROL(1,branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY(omga(i, j, n))
                  temp14 = pe0(i, k+1) - pe0(i, k)
                  temp_ad14 = omga_ad(i, j, n)/temp14
                  temp16 = dp2(i, n) - pe0(i, k)
                  temp15 = pe3(i, k+1) - pe3(i, k)
                  temp_ad15 = -(temp15*temp16*temp_ad14/temp14)
                  pe3_ad(i, k) = pe3_ad(i, k) + omga_ad(i, j, n) - &
&                   temp16*temp_ad14
                  pe3_ad(i, k+1) = pe3_ad(i, k+1) + temp16*temp_ad14
                  dp2_ad(i, n) = dp2_ad(i, n) + temp15*temp_ad14
                  pe0_ad(i, k) = pe0_ad(i, k) - temp_ad15 - temp15*&
&                   temp_ad14
                  pe0_ad(i, k+1) = pe0_ad(i, k+1) + temp_ad15
                  omga_ad(i, j, n) = 0.0
                END IF
                CALL POPINTEGER(ad_count0)
                DO i1=1,ad_count0
                  IF (i1 .EQ. 1) CALL POPCONTROL(1,branch)
                  CALL POPINTEGER(k)
                END DO
              END DO
            END DO
            CALL POPINTEGER(ad_to47)
            DO k=ad_to47,1,-1
              CALL POPINTEGER(ad_from26)
              CALL POPINTEGER(ad_to46)
              DO i=ad_to46,ad_from26,-1
                CALL POPREALARRAY(dp2(i, k))
                peln_ad(i, k, j) = peln_ad(i, k, j) + 0.5*dp2_ad(i, k)
                peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + 0.5*dp2_ad(i, &
&                 k)
                dp2_ad(i, k) = 0.0
              END DO
            END DO
          END IF
          CALL POPCONTROL(2,branch)
          IF (branch .EQ. 0) THEN
            CALL POPINTEGER(ad_to41)
            DO k=ad_to41,1,-1
              CALL POPINTEGER(ad_from23)
              CALL POPINTEGER(ad_to40)
              DO i=ad_to40,ad_from23,-1
                CALL POPREALARRAY(pkz(i, j, k))
                temp7 = akap*(peln(i, k+1, j)-peln(i, k, j))
                temp_ad10 = pkz_ad(i, j, k)/temp7
                temp_ad11 = -((pk2(i, k+1)-pk2(i, k))*akap*temp_ad10/&
&                 temp7)
                pk2_ad(i, k+1) = pk2_ad(i, k+1) + temp_ad10
                pk2_ad(i, k) = pk2_ad(i, k) - temp_ad10
                peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad11
                peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad11
                pkz_ad(i, j, k) = 0.0
              END DO
            END DO
          ELSE IF (branch .EQ. 1) THEN
            CALL POPINTEGER(ad_to43)
            DO k=ad_to43,1,-1
              CALL POPINTEGER(ad_from24)
              CALL POPINTEGER(ad_to42)
              DO i=ad_to42,ad_from24,-1
                CALL POPREALARRAY(pkz(i, j, k))
                temp10 = delz(i, j, k)
                temp9 = delp(i, j, k)*pt(i, j, k)
                temp8 = temp9/temp10
                temp_ad12 = akap*EXP(akap*LOG(rrg*temp8))*pkz_ad(i, j, k&
&                 )/(temp8*temp10)
                delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*&
&                 temp_ad12
                pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*&
&                 temp_ad12
                delz_ad(i, j, k) = delz_ad(i, j, k) - temp8*temp_ad12
                pkz_ad(i, j, k) = 0.0
              END DO
            END DO
          ELSE
            CALL POPINTEGER(ad_to45)
            DO k=ad_to45,1,-1
              CALL POPINTEGER(ad_from25)
              CALL POPINTEGER(ad_to44)
              DO i=ad_to44,ad_from25,-1
                CALL POPREALARRAY(pkz(i, j, k))
                temp13 = delz(i, j, k)
                temp12 = delp(i, j, k)*pt(i, j, k)
                temp11 = temp12/temp13
                temp_ad13 = k1k*EXP(k1k*LOG(rrg*temp11))*pkz_ad(i, j, k)&
&                 /(temp11*temp13)
                delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*&
&                 temp_ad13
                pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*&
&                 temp_ad13
                delz_ad(i, j, k) = delz_ad(i, j, k) - temp11*temp_ad13
                pkz_ad(i, j, k) = 0.0
              END DO
            END DO
          END IF
        END IF
      END IF
      CALL POPINTEGER(ad_to39)
      DO k=ad_to39,1,-1
        CALL POPINTEGER(ad_from22)
        CALL POPINTEGER(ad_to38)
        DO i=ad_to38,ad_from22,-1
          CALL POPREALARRAY(peln(i, k, j))
          pn2_ad(i, k) = pn2_ad(i, k) + peln_ad(i, k, j)
          peln_ad(i, k, j) = pe0_ad(i, k)
          CALL POPREALARRAY(pe0(i, k))
          pe0_ad(i, k) = 0.0
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        CALL POPINTEGER(ad_to37)
        DO k=ad_to37,2,-1
          CALL POPINTEGER(ad_from21)
          CALL POPINTEGER(ad_to36)
          DO i=ad_to36,ad_from21,-1
            CALL POPREALARRAY(pe3(i, k))
            omga_ad(i, j, k-1) = omga_ad(i, j, k-1) + pe3_ad(i, k)
            pe3_ad(i, k) = 0.0
          END DO
        END DO
        CALL POPINTEGER(ad_from20)
        CALL POPINTEGER(ad_to35)
        DO i=ad_to35,ad_from20,-1
          CALL POPREALARRAY(pe3(i, 1))
          pe3_ad(i, 1) = 0.0
        END DO
      END IF
      CALL POPINTEGER(ad_to34)
      DO k=ad_to34,1,-1
        CALL POPINTEGER(ad_from19)
        CALL POPINTEGER(ad_to33)
        DO i=ad_to33,ad_from19,-1
          CALL POPREALARRAY(pk(i, j, k))
          pk2_ad(i, k) = pk2_ad(i, k) + pk_ad(i, j, k)
          pk_ad(i, j, k) = 0.0
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        CALL POPINTEGER(ad_to32)
        DO k=ad_to32,1,-1
          CALL POPINTEGER(ad_from18)
          CALL POPINTEGER(ad_to31)
          DO i=ad_to31,ad_from18,-1
            CALL POPREALARRAY(delz(i, j, k))
            dp2_ad(i, k) = dp2_ad(i, k) - delz(i, j, k)*delz_ad(i, j, k)
            delz_ad(i, j, k) = -(dp2(i, k)*delz_ad(i, j, k))
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(delz, (ied-isd+1)*(jed-jsd+1)*km)
          gz_ad = 0.0
          CALL MAP1_PPM_ADM(km, pe1, pe1_ad, gz, gz_ad, km, pe2, pe2_ad&
&                     , delz, delz_ad, is, ie, j, isd, ied, jsd, jed, 1&
&                     , abs_kord_tm_pert)
        ELSE
          gz_ad = 0.0
          CALL MAP1_PPM_BWD(km, pe1, pe1_ad, gz, gz_ad, km, pe2, &
&                        pe2_ad, delz, delz_ad, is, ie, j, isd, ied, jsd&
&                        , jed, 1, abs_kord_tm)
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL MAP1_PPM_BWD(km, pe1, pe1_ad, ws(is:ie, j), ws_ad(is:&
&                        ie, j), km, pe2, pe2_ad, w, w_ad, is, ie, j, &
&                        isd, ied, jsd, jed, -2, kord_wz)
        ELSE
          CALL POPREALARRAY(w, (ied-isd+1)*(jed-jsd+1)*km)
          CALL MAP1_PPM_ADM(km, pe1, pe1_ad, ws(is:ie, j), ws_ad(is:ie, &
&                     j), km, pe2, pe2_ad, w, w_ad, is, ie, j, isd, ied&
&                     , jsd, jed, -2, kord_wz_pert)
        END IF
      END IF
      CALL POPCONTROL(2,branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          CALL MAPN_TRACER_BWD(nq, km, pe1, pe1_ad, pe2, pe2_ad, q, &
&                           q_ad, dp2, dp2_ad, kord_tr, j, is, ie, isd, &
&                           ied, jsd, jed, 0., fill)
        ELSE
          CALL POPREALARRAY(q, (ied-isd+1)*(jed-jsd+1)*km*nq)
          CALL MAPN_TRACER_ADM(nq, km, pe1, pe1_ad, pe2, pe2_ad, q, q_ad&
&                        , dp2, dp2_ad, kord_tr_pert, j, is, ie, isd, &
&                        ied, jsd, jed, 0., fill)
        END IF
      ELSE IF (branch .EQ. 2) THEN
        CALL POPINTEGER(ad_to30)
        DO iq=ad_to30,1,-1
          CALL POPINTEGER(ad_to29)
          DO k=ad_to29,1,-1
            CALL POPINTEGER(ad_from17)
            CALL POPINTEGER(ad_to28)
            DO i=ad_to28,ad_from17,-1
              CALL POPREALARRAY(q(i, j, k, iq))
              q2_ad(i, k) = q2_ad(i, k) + q_ad(i, j, k, iq)
              q_ad(i, j, k, iq) = 0.0
            END DO
          END DO
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL MAP1_Q2_BWD(km, pe1, pe1_ad, q(isd:ied, jsd:jed, 1:&
&                         km, iq), q_ad(isd:ied, jsd:jed, 1:km, iq), km&
&                         , pe2, pe2_ad, q2, q2_ad, dp2, dp2_ad, is, ie&
&                         , 0, kord_tr(iq), j, isd, ied, jsd, jed, 0.)
          ELSE
            CALL MAP1_Q2_ADM(km, pe1, pe1_ad, q(isd:ied, jsd:jed, 1:km, &
&                      iq), q_ad(isd:ied, jsd:jed, 1:km, iq), km, pe2, &
&                      pe2_ad, q2, q2_ad, dp2, dp2_ad, is, ie, 0, &
&                      kord_tr_pert(iq), j, isd, ied, jsd, jed, 0.)
          END IF
        END DO
      END IF
      CALL POPCONTROL(3,branch)
      IF (branch .LT. 3) THEN
        IF (branch .EQ. 0) THEN
          CALL MAP_SCALAR_BWD(km, peln(is:ie, 1:km+1, j), peln_ad(is:&
&                          ie, 1:km+1, j), gz, km, pn2, pn2_ad, pt, &
&                          pt_ad, is, ie, j, isd, ied, jsd, jed, 1, &
&                          abs_kord_tm, t_min)
        ELSE IF (branch .EQ. 1) THEN
          CALL POPREALARRAY(pt, (ied-isd+1)*(jed-jsd+1)*km)
          CALL MAP_SCALAR_ADM(km, peln(is:ie, 1:km+1, j), peln_ad(is:ie&
&                       , 1:km+1, j), gz, km, pn2, pn2_ad, pt, pt_ad, is&
&                       , ie, j, isd, ied, jsd, jed, 1, abs_kord_tm_pert&
&                       , t_min)
        ELSE
          gz_ad = 0.0
          CALL MAP1_PPM_BWD(km, pe1, pe1_ad, gz, gz_ad, km, pe2, &
&                        pe2_ad, pt, pt_ad, is, ie, j, isd, ied, jsd, &
&                        jed, 1, abs_kord_tm)
        END IF
      ELSE IF (branch .EQ. 3) THEN
        CALL POPREALARRAY(pt, (ied-isd+1)*(jed-jsd+1)*km)
        gz_ad = 0.0
        CALL MAP1_PPM_ADM(km, pe1, pe1_ad, gz, gz_ad, km, pe2, pe2_ad, &
&                   pt, pt_ad, is, ie, j, isd, ied, jsd, jed, 1, &
&                   abs_kord_tm_pert)
      ELSE IF (branch .EQ. 4) THEN
        CALL MAP1_CUBIC_BWD(km, pe1, pe1_ad, km, pe2, pe2_ad, te, te_ad&
&                     , is, ie, j, isd, ied, jsd, jed, akap, t_var=1, &
&                     conserv=.true.)
        CALL POPINTEGER(ad_to27)
        DO k=ad_to27,1,-1
          CALL POPINTEGER(ad_from16)
          CALL POPINTEGER(ad_to26)
          DO i=ad_to26,ad_from16,-1
            CALL POPREALARRAY(te(i, j, k))
            temp6 = pe1(i, k+1) - pe1(i, k)
            temp_ad8 = te_ad(i, j, k)/temp6
            temp_ad9 = -((phis(i, k+1)-phis(i, k))*temp_ad8/temp6)
            phis_ad(i, k+1) = phis_ad(i, k+1) + temp_ad8
            phis_ad(i, k) = phis_ad(i, k) - temp_ad8
            pe1_ad(i, k+1) = pe1_ad(i, k+1) + temp_ad9
            pe1_ad(i, k) = pe1_ad(i, k) - temp_ad9
          END DO
        END DO
        CALL POPINTEGER(ad_to25)
        DO k=ad_to25,1,-1
          CALL POPINTEGER(ad_from15)
          CALL POPINTEGER(ad_to24)
          DO i=ad_to24,ad_from15,-1
            CALL POPREALARRAY(phis(i, k))
            pe1_ad(i, k) = pe1_ad(i, k) + phis(i, k)*phis_ad(i, k)
            phis_ad(i, k) = pe1(i, k)*phis_ad(i, k)
          END DO
        END DO
        CALL POPINTEGER(ad_from14)
        DO k=1,ad_from14,1
          CALL POPINTEGER(ad_from13)
          CALL POPINTEGER(ad_to23)
          DO i=ad_to23,ad_from13,-1
            CALL POPREALARRAY(phis(i, k))
            temp_ad7 = cp_air*pt(i, j, k)*phis_ad(i, k)
            phis_ad(i, k+1) = phis_ad(i, k+1) + phis_ad(i, k)
            pt_ad(i, j, k) = pt_ad(i, j, k) + cp_air*(pk1(i, k+1)-pk1(i&
&             , k))*phis_ad(i, k)
            pk1_ad(i, k+1) = pk1_ad(i, k+1) + temp_ad7
            pk1_ad(i, k) = pk1_ad(i, k) - temp_ad7
            phis_ad(i, k) = 0.0
          END DO
        END DO
        CALL POPINTEGER(ad_from12)
        CALL POPINTEGER(ad_to22)
        DO i=ad_to22,ad_from12,-1
          CALL POPREALARRAY(phis(i, km+1))
          phis_ad(i, km+1) = 0.0
        END DO
      END IF
      CALL POPINTEGER(ad_to21)
      DO k=ad_to21,2,-1
        CALL POPINTEGER(ad_from11)
        CALL POPINTEGER(ad_to20)
        DO i=ad_to20,ad_from11,-1
          CALL POPREALARRAY(pk2(i, k))
          pn2_ad(i, k) = pn2_ad(i, k) + EXP(akap*pn2(i, k))*akap*pk2_ad(&
&           i, k)
          pk2_ad(i, k) = 0.0
          CALL POPREALARRAY(pn2(i, k))
          pe2_ad(i, k) = pe2_ad(i, k) + pn2_ad(i, k)/pe2(i, k)
          pn2_ad(i, k) = 0.0
        END DO
      END DO
      CALL POPINTEGER(ad_from10)
      CALL POPINTEGER(ad_to19)
      DO i=ad_to19,ad_from10,-1
        CALL POPREALARRAY(pk2(i, km+1))
        pk1_ad(i, km+1) = pk1_ad(i, km+1) + pk2_ad(i, km+1)
        pk2_ad(i, km+1) = 0.0
        CALL POPREALARRAY(pk2(i, 1))
        pk1_ad(i, 1) = pk1_ad(i, 1) + pk2_ad(i, 1)
        pk2_ad(i, 1) = 0.0
        CALL POPREALARRAY(pn2(i, km+1))
        peln_ad(i, km+1, j) = peln_ad(i, km+1, j) + pn2_ad(i, km+1)
        pn2_ad(i, km+1) = 0.0
        CALL POPREALARRAY(pn2(i, 1))
        peln_ad(i, 1, j) = peln_ad(i, 1, j) + pn2_ad(i, 1)
        pn2_ad(i, 1) = 0.0
      END DO
      CALL POPINTEGER(ad_to18)
      DO k=ad_to18,1,-1
        CALL POPINTEGER(ad_from9)
        CALL POPINTEGER(ad_to17)
        DO i=ad_to17,ad_from9,-1
          CALL POPREALARRAY(pk1(i, k))
          pk_ad(i, j, k) = pk_ad(i, j, k) + pk1_ad(i, k)
          pk1_ad(i, k) = 0.0
        END DO
      END DO
      CALL POPINTEGER(ad_to16)
      DO k=ad_to16,1,-1
        CALL POPINTEGER(ad_from8)
        CALL POPINTEGER(ad_to15)
        DO i=ad_to15,ad_from8,-1
          CALL POPREALARRAY(delp(i, j, k))
          dp2_ad(i, k) = dp2_ad(i, k) + delp_ad(i, j, k)
          delp_ad(i, j, k) = 0.0
        END DO
      END DO
      CALL POPINTEGER(ad_to14)
      DO k=ad_to14,1,-1
        CALL POPINTEGER(ad_from7)
        CALL POPINTEGER(ad_to13)
        DO i=ad_to13,ad_from7,-1
          CALL POPREALARRAY(dp2(i, k))
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp2_ad(i, k)
          pe2_ad(i, k) = pe2_ad(i, k) - dp2_ad(i, k)
          dp2_ad(i, k) = 0.0
        END DO
      END DO
      CALL POPINTEGER(ad_to12)
      DO k=ad_to12,2,-1
        CALL POPINTEGER(ad_from6)
        CALL POPINTEGER(ad_to11)
        DO i=ad_to11,ad_from6,-1
          CALL POPREALARRAY(pe2(i, k))
          pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + bk(k)*pe2_ad(i, k)
          pe2_ad(i, k) = 0.0
        END DO
      END DO
      CALL POPINTEGER(ad_from5)
      CALL POPINTEGER(ad_to10)
      DO i=ad_to10,ad_from5,-1
        pe1_ad(i, km+1) = pe1_ad(i, km+1) + ps_ad(i, j)
        ps_ad(i, j) = 0.0
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        CALL POPINTEGER(ad_to9)
        DO k=ad_to9,1,-1
          CALL POPINTEGER(ad_from4)
          CALL POPINTEGER(ad_to8)
          DO i=ad_to8,ad_from4,-1
            CALL POPREALARRAY(delz(i, j, k))
            temp_ad6 = -(delz_ad(i, j, k)/delp(i, j, k))
            delp_ad(i, j, k) = delp_ad(i, j, k) - delz(i, j, k)*temp_ad6&
&             /delp(i, j, k)
            delz_ad(i, j, k) = temp_ad6
          END DO
        END DO
      END IF
      CALL POPCONTROL(3,branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          CALL POPINTEGER(ad_to3)
          DO k=ad_to3,1,-1
            CALL POPINTEGER(ad_from1)
            CALL POPINTEGER(ad_to2)
            DO i=ad_to2,ad_from1,-1
              CALL POPREALARRAY(pt(i, j, k))
              temp1 = akap*(peln(i, k+1, j)-peln(i, k, j))
              temp_ad = pt_ad(i, j, k)/temp1
              temp0 = pk(i, j, k+1) - pk(i, j, k)
              temp = pt(i, j, k)
              temp_ad0 = -(temp*temp0*akap*temp_ad/temp1)
              pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp*temp_ad
              pk_ad(i, j, k) = pk_ad(i, j, k) - temp*temp_ad
              peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad0
              peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad0
              pt_ad(i, j, k) = temp0*temp_ad
            END DO
          END DO
        ELSE
          CALL POPINTEGER(ad_to5)
          DO k=ad_to5,1,-1
            CALL POPINTEGER(ad_from2)
            CALL POPINTEGER(ad_to4)
            DO i=ad_to4,ad_from2,-1
              CALL POPREALARRAY(pt(i, j, k))
              temp5 = delz(i, j, k)
              temp4 = delp(i, j, k)*pt(i, j, k)
              temp2 = temp4/temp5
              temp3 = k1k*LOG(rrg*temp2)
              temp_ad1 = k1k*EXP(temp3)*pt(i, j, k)*pt_ad(i, j, k)/(&
&               temp2*temp5)
              delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*temp_ad1
              delz_ad(i, j, k) = delz_ad(i, j, k) - temp2*temp_ad1
              pt_ad(i, j, k) = delp(i, j, k)*temp_ad1 + EXP(temp3)*pt_ad&
&               (i, j, k)
            END DO
          END DO
        END IF
      ELSE IF (branch .NE. 2) THEN
        IF (branch .EQ. 3) THEN
          CALL POPINTEGER(ad_to7)
          DO k=ad_to7,1,-1
            CALL POPINTEGER(ad_from3)
            CALL POPINTEGER(ad_to6)
            DO i=ad_to6,ad_from3,-1
              CALL POPREALARRAY(te(i, j, k))
              temp_ad2 = gridstruct%rsin2(i, j)*0.25*te_ad(i, j, k)
              temp_ad3 = -(gridstruct%cosa_s(i, j)*temp_ad2)
              temp_ad4 = (v(i, j, k)+v(i+1, j, k))*temp_ad3
              temp_ad5 = (u(i, j, k)+u(i, j+1, k))*temp_ad3
              u_ad(i, j, k) = u_ad(i, j, k) + temp_ad4 + 2*u(i, j, k)*&
&               temp_ad2
              u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp_ad4 + 2*u(i, j+1&
&               , k)*temp_ad2
              v_ad(i, j, k) = v_ad(i, j, k) + temp_ad5 + 2*v(i, j, k)*&
&               temp_ad2
              v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp_ad5 + 2*v(i+1, j&
&               , k)*temp_ad2
              pt_ad(i, j, k) = pt_ad(i, j, k) + cp_air*pkz(i, j, k)*&
&               te_ad(i, j, k)
              pkz_ad(i, j, k) = pkz_ad(i, j, k) + cp_air*pt(i, j, k)*&
&               te_ad(i, j, k)
              te_ad(i, j, k) = 0.0
            END DO
          END DO
          CALL PKEZ_BWD(km, is, ie, js, je, j, pe, pk, pk_ad, akap, peln&
&                 , peln_ad, pkz, pkz_ad, ptop)
        END IF
      END IF
 140  CALL POPINTEGER(ad_from0)
      CALL POPINTEGER(ad_to1)
      DO i=ad_to1,ad_from0,-1
        CALL POPREALARRAY(pe2(i, km+1))
        pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + pe2_ad(i, km+1)
        pe2_ad(i, km+1) = 0.0
        CALL POPREALARRAY(pe2(i, 1))
        pe2_ad(i, 1) = 0.0
      END DO
      CALL POPINTEGER(ad_to0)
      DO k=ad_to0,1,-1
        CALL POPINTEGER(ad_from)
        CALL POPINTEGER(ad_to)
        DO i=ad_to,ad_from,-1
          CALL POPREALARRAY(pe1(i, k))
          pe_ad(i, k, j) = pe_ad(i, k, j) + pe1_ad(i, k)
          pe1_ad(i, k) = 0.0
        END DO
      END DO
 150  CALL POPINTEGER(j)
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      CALL POPINTEGER(ad_count)
      DO i0=1,ad_count
        IF (i0 .EQ. 1) CALL POPCONTROL(1,branch)
      END DO
    END IF
    CALL POPCONTROL(1,branch)
  END SUBROUTINE LAGRANGIAN_TO_EULERIAN_BWD
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
    INTEGER :: arg1
    REAL :: result1
    LOGICAL :: arg10
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
                pt(i, j, k) = pt(i, j, k)*EXP(k1k*LOG(rrg*delp(i, j, k)/&
&                 delz(i, j, k)*pt(i, j, k)))
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
            CALL MAP_SCALAR(km, peln(is:ie, 1:km+1, j), gz, km, pn2, pt&
&                     , is, ie, j, isd, ied, jsd, jed, 1, &
&                     abs_kord_tm_pert, t_min)
          END IF
        ELSE IF (remap_pt) THEN
!----------------------------------
! Map pt using pe
!----------------------------------
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            CALL MAP1_PPM(km, pe1, gz, km, pe2, pt, is, ie, j, isd, &
&                      ied, jsd, jed, 1, abs_kord_tm)
          ELSE
            CALL MAP1_PPM(km, pe1, gz, km, pe2, pt, is, ie, j, isd, ied&
&                   , jsd, jed, 1, abs_kord_tm_pert)
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
            CALL MAPN_TRACER(nq, km, pe1, pe2, q, dp2, kord_tr_pert, j, &
&                      is, ie, isd, ied, jsd, jed, 0., fill)
          END IF
        ELSE IF (nq .GT. 0) THEN
! Remap one tracer at a time
          DO iq=1,nq
            IF (kord_tr(iq) .EQ. kord_tr_pert(iq)) THEN
              CALL MAP1_Q2(km, pe1, q(isd:ied, jsd:jed, 1:km, iq), km&
&                       , pe2, q2, dp2, is, ie, 0, kord_tr(iq), j, isd, &
&                       ied, jsd, jed, 0.)
            ELSE
              CALL MAP1_Q2(km, pe1, q(isd:ied, jsd:jed, 1:km, iq), km, &
&                    pe2, q2, dp2, is, ie, 0, kord_tr_pert(iq), j, isd, &
&                    ied, jsd, jed, 0.)
            END IF
            IF (fill) THEN
              arg1 = ie - is + 1
              CALL FILLZ(arg1, km, 1, q2, dp2)
            END IF
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
            CALL MAP1_PPM(km, pe1, ws(is:ie, j), km, pe2, w, is, ie, j, &
&                   isd, ied, jsd, jed, -2, kord_wz_pert)
          END IF
! Remap delz for hybrid sigma-p coordinate
          IF (abs_kord_tm .EQ. abs_kord_tm_pert) THEN
            CALL MAP1_PPM(km, pe1, gz, km, pe2, delz, is, ie, j, isd&
&                      , ied, jsd, jed, 1, abs_kord_tm)
          ELSE
            CALL MAP1_PPM(km, pe1, gz, km, pe2, delz, is, ie, j, isd, &
&                   ied, jsd, jed, 1, abs_kord_tm_pert)
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
              pkz(i, j, k) = EXP(akap*LOG(rrg*delp(i, j, k)/delz(i, j, k&
&               )*pt(i, j, k)))
            END DO
          END DO
        ELSE
! Using dry pressure for the definition of the virtual potential temperature
!           pkz(i,j,k) = exp(akap*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
! Note: pt at this stage is Theta_v
          DO k=1,km
            DO i=is,ie
              pkz(i, j, k) = EXP(k1k*LOG(rrg*delp(i, j, k)/delz(i, j, k)&
&               *pt(i, j, k)))
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
        CALL MAP1_PPM(km, pe0(is:ie, :), gz, km, pe3(is:ie, :), u, is, &
&               ie, j, isd, ied, jsd, jedp1, -1, kord_mt_pert)
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
          CALL MAP1_PPM(km, pe0, gz, km, pe3, v, is, iep1, j, isd, iedp1&
&                 , jsd, jed, -1, kord_mt_pert)
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
        arg10 = cld_amt .GT. 0
        CALL FV_SAT_ADJ(abs0, r_vir, is, ie, js, je, ng, hydrostatic, &
&                 fast_mp_consv, te(isd:ied, jsd:jed, k), q(isd:ied, jsd&
&                 :jed, k, sphum), q(isd:ied, jsd:jed, k, liq_wat), q(&
&                 isd:ied, jsd:jed, k, ice_wat), q(isd:ied, jsd:jed, k, &
&                 rainwat), q(isd:ied, jsd:jed, k, snowwat), q(isd:ied, &
&                 jsd:jed, k, graupel), dpln, delz(isd:ied, jsd:jed, k)&
&                 , pt(isd:ied, jsd:jed, k), delp(isd:ied, jsd:jed, k), &
&                 q_con(isd:ied, jsd:jed, k), cappa(isd:ied, jsd:jed, k)&
&                 , gridstruct%area_64, dtdt(is:ie, js:je, k), out_dt, &
&                 last_step, arg10, q(isd:ied, jsd:jed, k, cld_amt))
        IF (.NOT.hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              pkz(i, j, k) = EXP(akap*LOG(rrg*delp(i, j, k)/delz(i, j, k&
&               )*pt(i, j, k)))
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
!  Differentiation of compute_total_energy in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 
!a2b_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe
! dyn_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_
!core_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mo
!d.Rayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_m
!od.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_m
!apz_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.
!remap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm
!_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_
!cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod
!.fv_subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d
! nh_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver n
!h_utils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh
!_utils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_c
!ore_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_
!mod.ytp_v sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d_f
!b tp_core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_
!grid_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: qc peln q u v w teq delp delz
!                te_2d pe pt
!   with respect to varying inputs: qc peln q u v w delp delz pe
!                pt
  SUBROUTINE COMPUTE_TOTAL_ENERGY_FWD(is, ie, js, je, isd, ied, jsd, jed&
&   , km, u, v, w, delz, pt, delp, q, qc, pe, peln, hs, rsin2_l, &
&   cosa_s_l, r_vir, cp, rg, hlv, te_2d, ua, va, teq, moist_phys, nwat, &
&   sphum, liq_wat, rainwat, ice_wat, snowwat, graupel, hydrostatic, &
&   id_te)
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
    REAL :: te_2d(is:ie, js:je)
! Moist TE
    REAL :: teq(is:ie, js:je)
! Local
    REAL, DIMENSION(is:ie, km) :: tv
    REAL :: phiz(is:ie, km+1)
    REAL :: cvm(is:ie), qd(is:ie)
    INTEGER :: i, j, k

    tv = 0.0
    phiz = 0.0
    cvm = 0.0
    qd = 0.0

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
          CALL PUSHREALARRAY(phiz(i, km+1))
          phiz(i, km+1) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            CALL PUSHREALARRAY(tv(i, k))
            tv(i, k) = pt(i, j, k)*(1.+qc(i, j, k))
            CALL PUSHREALARRAY(phiz(i, k))
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
        CALL PUSHCONTROL(2,2)
      ELSE
!-----------------
! Non-hydrostatic:
!-----------------
        DO i=is,ie
          CALL PUSHREALARRAY(phiz(i, km+1))
          phiz(i, km+1) = hs(i, j)
          DO k=km,1,-1
            CALL PUSHREALARRAY(phiz(i, k))
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
          CALL PUSHCONTROL(2,1)
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
          CALL PUSHCONTROL(2,0)
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
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHREALARRAY(tv, (ie-is+1)*km)
      CALL PUSHREALARRAY(phiz, (ie-is+1)*(km+1))
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHREALARRAY(tv, (ie-is+1)*km)
      CALL PUSHREALARRAY(phiz, (ie-is+1)*(km+1))
      CALL PUSHCONTROL(1,0)
    END IF
  END SUBROUTINE COMPUTE_TOTAL_ENERGY_FWD
!  Differentiation of compute_total_energy in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4
! a2b_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_p
!e dyn_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn
!_core_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_m
!od.Rayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_
!mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_
!mapz_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod
!.remap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.pp
!m_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1
!_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mo
!d.fv_subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_
!d nh_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver 
!nh_utils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile n
!h_utils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_
!core_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core
!_mod.ytp_v sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d_
!fb tp_core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv
!_grid_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: qc peln q u v w teq delp delz
!                te_2d pe pt
!   with respect to varying inputs: qc peln q u v w delp delz pe
!                pt
  SUBROUTINE COMPUTE_TOTAL_ENERGY_BWD(is, ie, js, je, isd, ied, jsd, jed&
&   , km, u, u_ad, v, v_ad, w, w_ad, delz, delz_ad, pt, pt_ad, delp, &
&   delp_ad, q, q_ad, qc, qc_ad, pe, pe_ad, peln, peln_ad, hs, rsin2_l, &
&   cosa_s_l, r_vir, cp, rg, hlv, te_2d, te_2d_ad, ua, va, teq, teq_ad, &
&   moist_phys, nwat, sphum, liq_wat, rainwat, ice_wat, snowwat, graupel&
&   , hydrostatic, id_te)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km, is, ie, js, je, isd, ied, jsd, jed, id_te
    INTEGER, INTENT(IN) :: sphum, liq_wat, ice_wat, rainwat, snowwat, &
&   graupel, nwat
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt, delp
    REAL, DIMENSION(isd:ied, jsd:jed, km) :: pt_ad, delp_ad
    REAL, DIMENSION(isd:ied, jsd:jed, km, *), INTENT(IN) :: q
    REAL, DIMENSION(isd:ied, jsd:jed, km, *) :: q_ad
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: qc
    REAL, DIMENSION(isd:ied, jsd:jed, km) :: qc_ad
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: u_ad(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_ad(isd:ied+1, jsd:jed, km)
    REAL, INTENT(IN) :: w(isd:ied, jsd:jed, km)
    REAL :: w_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delz(isd:ied, jsd:jed, km)
    REAL :: delz_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL :: pe_ad(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(IN) :: peln(is:ie, km+1, js:je)
    REAL :: peln_ad(is:ie, km+1, js:je)
    REAL, INTENT(IN) :: cp, rg, r_vir, hlv
    REAL, INTENT(IN) :: rsin2_l(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: cosa_s_l(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: moist_phys, hydrostatic
    REAL :: te_2d(is:ie, js:je)
    REAL :: te_2d_ad(is:ie, js:je)
    REAL :: teq(is:ie, js:je)
    REAL :: teq_ad(is:ie, js:je)
    REAL, DIMENSION(is:ie, km) :: tv
    REAL, DIMENSION(is:ie, km) :: tv_ad
    REAL :: phiz(is:ie, km+1)
    REAL :: phiz_ad(is:ie, km+1)
    REAL :: cvm(is:ie), qd(is:ie)
    INTEGER :: i, j, k
    REAL :: temp_ad
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp2
    REAL :: temp3
    REAL :: temp4
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    INTEGER :: branch

    tv = 0.0
    phiz = 0.0
    cvm = 0.0
    qd = 0.0
    branch = 0

    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(phiz, (ie-is+1)*(km+1))
      CALL POPREALARRAY(tv, (ie-is+1)*km)
    ELSE
      CALL POPREALARRAY(phiz, (ie-is+1)*(km+1))
      CALL POPREALARRAY(tv, (ie-is+1)*km)
      DO j=je,js,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          DO k=km,1,-1
            DO i=ie,is,-1
              q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) + hlv*delp(i, &
&               j, k)*teq_ad(i, j)
              delp_ad(i, j, k) = delp_ad(i, j, k) + hlv*q(i, j, k, sphum&
&               )*teq_ad(i, j)
            END DO
          END DO
        END IF
        DO i=ie,is,-1
          te_2d_ad(i, j) = te_2d_ad(i, j) + teq_ad(i, j)
          teq_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    phiz_ad = 0.0
    tv_ad = 0.0
    DO 100 j=je,js,-1
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        DO k=km,1,-1
          DO i=ie,is,-1
            temp7 = v(i, j, k) + v(i+1, j, k)
            temp6 = u(i, j, k) + u(i, j+1, k)
            temp5 = 0.5*rsin2_l(i, j)
            temp_ad7 = delp(i, j, k)*te_2d_ad(i, j)
            temp_ad8 = 0.5*temp_ad7
            temp_ad9 = temp5*temp_ad8
            temp_ad10 = -(cosa_s_l(i, j)*temp_ad9)
            delp_ad(i, j, k) = delp_ad(i, j, k) + (cv_air*pt(i, j, k)+&
&             0.5*(phiz(i, k)+phiz(i, k+1)+w(i, j, k)**2+temp5*(u(i, j, &
&             k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-&
&             cosa_s_l(i, j)*(temp6*temp7))))*te_2d_ad(i, j)
            pt_ad(i, j, k) = pt_ad(i, j, k) + cv_air*temp_ad7
            phiz_ad(i, k) = phiz_ad(i, k) + temp_ad8
            phiz_ad(i, k+1) = phiz_ad(i, k+1) + temp_ad8
            w_ad(i, j, k) = w_ad(i, j, k) + 2*w(i, j, k)*temp_ad8
            u_ad(i, j, k) = u_ad(i, j, k) + temp7*temp_ad10 + 2*u(i, j, &
&             k)*temp_ad9
            u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp7*temp_ad10 + 2*u(i&
&             , j+1, k)*temp_ad9
            v_ad(i, j, k) = v_ad(i, j, k) + temp6*temp_ad10 + 2*v(i, j, &
&             k)*temp_ad9
            v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp6*temp_ad10 + 2*v(i+&
&             1, j, k)*temp_ad9
          END DO
        END DO
      ELSE IF (branch .EQ. 1) THEN
        DO k=km,1,-1
          DO i=ie,is,-1
            temp4 = v(i, j, k) + v(i+1, j, k)
            temp3 = u(i, j, k) + u(i, j+1, k)
            temp2 = 0.5*rsin2_l(i, j)
            temp_ad3 = delp(i, j, k)*te_2d_ad(i, j)
            temp_ad4 = 0.5*temp_ad3
            temp_ad5 = temp2*temp_ad4
            temp_ad6 = -(cosa_s_l(i, j)*temp_ad5)
            delp_ad(i, j, k) = delp_ad(i, j, k) + (cv_air*pt(i, j, k)+&
&             0.5*(phiz(i, k)+phiz(i, k+1)+w(i, j, k)**2+temp2*(u(i, j, &
&             k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-&
&             cosa_s_l(i, j)*(temp3*temp4))))*te_2d_ad(i, j)
            pt_ad(i, j, k) = pt_ad(i, j, k) + cv_air*temp_ad3
            phiz_ad(i, k) = phiz_ad(i, k) + temp_ad4
            phiz_ad(i, k+1) = phiz_ad(i, k+1) + temp_ad4
            w_ad(i, j, k) = w_ad(i, j, k) + 2*w(i, j, k)*temp_ad4
            u_ad(i, j, k) = u_ad(i, j, k) + temp4*temp_ad6 + 2*u(i, j, k&
&             )*temp_ad5
            u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp4*temp_ad6 + 2*u(i, &
&             j+1, k)*temp_ad5
            v_ad(i, j, k) = v_ad(i, j, k) + temp3*temp_ad6 + 2*v(i, j, k&
&             )*temp_ad5
            v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp3*temp_ad6 + 2*v(i+1&
&             , j, k)*temp_ad5
          END DO
        END DO
      ELSE
        DO k=km,1,-1
          DO i=ie,is,-1
            temp1 = v(i, j, k) + v(i+1, j, k)
            temp0 = u(i, j, k) + u(i, j+1, k)
            temp = 0.25*rsin2_l(i, j)
            temp_ad0 = delp(i, j, k)*te_2d_ad(i, j)
            temp_ad1 = temp*temp_ad0
            temp_ad2 = -(cosa_s_l(i, j)*temp_ad1)
            delp_ad(i, j, k) = delp_ad(i, j, k) + (cp*tv(i, k)+temp*(u(i&
&             , j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-&
&             cosa_s_l(i, j)*(temp0*temp1)))*te_2d_ad(i, j)
            tv_ad(i, k) = tv_ad(i, k) + cp*temp_ad0
            u_ad(i, j, k) = u_ad(i, j, k) + temp1*temp_ad2 + 2*u(i, j, k&
&             )*temp_ad1
            u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp1*temp_ad2 + 2*u(i, &
&             j+1, k)*temp_ad1
            v_ad(i, j, k) = v_ad(i, j, k) + temp0*temp_ad2 + 2*v(i, j, k&
&             )*temp_ad1
            v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp0*temp_ad2 + 2*v(i+1&
&             , j, k)*temp_ad1
          END DO
        END DO
        DO i=ie,is,-1
          pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + phiz(i, km+1)*te_2d_ad&
&           (i, j)
          phiz_ad(i, km+1) = phiz_ad(i, km+1) + pe(i, km+1, j)*te_2d_ad(&
&           i, j)
          pe_ad(i, 1, j) = pe_ad(i, 1, j) - phiz(i, 1)*te_2d_ad(i, j)
          phiz_ad(i, 1) = phiz_ad(i, 1) - pe(i, 1, j)*te_2d_ad(i, j)
          te_2d_ad(i, j) = 0.0
        END DO
        DO k=1,km,1
          DO i=ie,is,-1
            CALL POPREALARRAY(phiz(i, k))
            temp_ad = rg*tv(i, k)*phiz_ad(i, k)
            phiz_ad(i, k+1) = phiz_ad(i, k+1) + phiz_ad(i, k)
            tv_ad(i, k) = tv_ad(i, k) + rg*(peln(i, k+1, j)-peln(i, k, j&
&             ))*phiz_ad(i, k)
            peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad
            peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad
            phiz_ad(i, k) = 0.0
            CALL POPREALARRAY(tv(i, k))
            pt_ad(i, j, k) = pt_ad(i, j, k) + (qc(i, j, k)+1.)*tv_ad(i, &
&             k)
            qc_ad(i, j, k) = qc_ad(i, j, k) + pt(i, j, k)*tv_ad(i, k)
            tv_ad(i, k) = 0.0
          END DO
        END DO
        DO i=ie,is,-1
          CALL POPREALARRAY(phiz(i, km+1))
          phiz_ad(i, km+1) = 0.0
        END DO
        GOTO 100
      END IF
      DO i=ie,is,-1
        te_2d_ad(i, j) = 0.0
      END DO
      DO i=ie,is,-1
        DO k=1,km,1
          CALL POPREALARRAY(phiz(i, k))
          phiz_ad(i, k+1) = phiz_ad(i, k+1) + phiz_ad(i, k)
          delz_ad(i, j, k) = delz_ad(i, j, k) - grav*phiz_ad(i, k)
          phiz_ad(i, k) = 0.0
        END DO
        CALL POPREALARRAY(phiz(i, km+1))
        phiz_ad(i, km+1) = 0.0
      END DO
 100 CONTINUE
  END SUBROUTINE COMPUTE_TOTAL_ENERGY_BWD
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
!  Differentiation of pkez in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b
!_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_
!grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp 
!dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super
! fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_g
!rid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z 
!fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz
!_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_map
!z_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart
!_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z ma
!in_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Ri
!em_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3
!p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_
!halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_ve
!ct sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_
!core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.co
!py_corners_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.g
!reat_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: peln pkz pk
!   with respect to varying inputs: peln pkz pk
  SUBROUTINE PKEZ_FWD(km, ifirst, ilast, jfirst, jlast, j, pe, pk, akap&
&   , peln, pkz, ptop)
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
    REAL :: pkz(ifirst:ilast, jfirst:jlast, km)
! log (pe)
    REAL, INTENT(INOUT) :: peln(ifirst:ilast, km+1, jfirst:jlast)
! Local
    REAL :: pk2(ifirst:ilast, km+1)
    REAL :: pek
    REAL :: lnp
    REAL :: ak1
    INTEGER :: i, k
    INTRINSIC LOG

    pk2 = 0.0
    pek = 0.0
    lnp = 0.0
    ak1 = 0.0

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
        CALL PUSHREALARRAY(peln(i, 1, j))
        peln(i, 1, j) = peln(i, 2, j) - ak1
      END DO
      CALL PUSHCONTROL(1,1)
    ELSE
      lnp = LOG(ptop)
      DO i=ifirst,ilast
        CALL PUSHREALARRAY(peln(i, 1, j))
        peln(i, 1, j) = lnp
      END DO
      CALL PUSHCONTROL(1,0)
    END IF
!---- GFDL modification
    DO k=1,km
      DO i=ifirst,ilast
        CALL PUSHREALARRAY(pkz(i, j, k))
        pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1, j)-&
&         peln(i, k, j)))
      END DO
    END DO
    CALL PUSHREALARRAY(pk2, (ilast-ifirst+1)*(km+1))
  END SUBROUTINE PKEZ_FWD
!  Differentiation of pkez in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2
!b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p
!_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp
! dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Supe
!r fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_
!grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z
! fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_map
!z_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_ma
!pz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restar
!t_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z m
!ain_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.R
!iem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM
!3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest
!_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_v
!ect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw
!_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.c
!opy_corners_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.
!great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: peln pkz pk
!   with respect to varying inputs: peln pkz pk
  SUBROUTINE PKEZ_BWD(km, ifirst, ilast, jfirst, jlast, j, pe, pk, pk_ad&
&   , akap, peln, peln_ad, pkz, pkz_ad, ptop)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km, j
    INTEGER, INTENT(IN) :: ifirst, ilast
    INTEGER, INTENT(IN) :: jfirst, jlast
    REAL, INTENT(IN) :: akap
    REAL, INTENT(IN) :: pe(ifirst-1:ilast+1, km+1, jfirst-1:jlast+1)
    REAL, INTENT(IN) :: pk(ifirst:ilast, jfirst:jlast, km+1)
    REAL :: pk_ad(ifirst:ilast, jfirst:jlast, km+1)
    REAL, INTENT(IN) :: ptop
    REAL :: pkz(ifirst:ilast, jfirst:jlast, km)
    REAL :: pkz_ad(ifirst:ilast, jfirst:jlast, km)
    REAL, INTENT(INOUT) :: peln(ifirst:ilast, km+1, jfirst:jlast)
    REAL, INTENT(INOUT) :: peln_ad(ifirst:ilast, km+1, jfirst:jlast)
    REAL :: pk2(ifirst:ilast, km+1)
    REAL :: pk2_ad(ifirst:ilast, km+1)
    REAL :: pek
    REAL :: pek_ad
    REAL :: lnp
    REAL :: ak1
    INTEGER :: i, k
    INTRINSIC LOG
    REAL :: temp
    REAL :: temp_ad
    REAL :: temp_ad0
    INTEGER :: branch

    pk2 = 0.0
    pek = 0.0
    lnp = 0.0
    ak1 = 0.0
    branch = 0

    CALL POPREALARRAY(pk2, (ilast-ifirst+1)*(km+1))
    pk2_ad = 0.0
    DO k=km,1,-1
      DO i=ilast,ifirst,-1
        CALL POPREALARRAY(pkz(i, j, k))
        temp = akap*(peln(i, k+1, j)-peln(i, k, j))
        temp_ad = pkz_ad(i, j, k)/temp
        temp_ad0 = -((pk2(i, k+1)-pk2(i, k))*akap*temp_ad/temp)
        pk2_ad(i, k+1) = pk2_ad(i, k+1) + temp_ad
        pk2_ad(i, k) = pk2_ad(i, k) - temp_ad
        peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad0
        peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad0
        pkz_ad(i, j, k) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO i=ilast,ifirst,-1
        CALL POPREALARRAY(peln(i, 1, j))
        peln_ad(i, 1, j) = 0.0
      END DO
    ELSE
      DO i=ilast,ifirst,-1
        CALL POPREALARRAY(peln(i, 1, j))
        peln_ad(i, 2, j) = peln_ad(i, 2, j) + peln_ad(i, 1, j)
        peln_ad(i, 1, j) = 0.0
      END DO
    END IF
    DO k=km+1,2,-1
      DO i=ilast,ifirst,-1
        pk_ad(i, j, k) = pk_ad(i, j, k) + pk2_ad(i, k)
        pk2_ad(i, k) = 0.0
      END DO
    END DO
    pek_ad = 0.0
    DO i=ilast,ifirst,-1
      pek_ad = pek_ad + pk2_ad(i, 1)
      pk2_ad(i, 1) = 0.0
    END DO
    pk_ad(ifirst, j, 1) = pk_ad(ifirst, j, 1) + pek_ad
  END SUBROUTINE PKEZ_BWD
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
!  Differentiation of map_scalar in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn
!_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dy
!n_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_
!mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynam
!ics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils
!_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_m
!od.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scal
!ar_profile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.ste
!epz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_
!setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.co
!mpute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver
!_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver
! nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh s
!w_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_cor
!e_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.
!compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corner
!s_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circ
!le_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 q2
!   with respect to varying inputs: pe1 pe2 q2
  SUBROUTINE MAP_SCALAR_ADM(km, pe1, pe1_ad, qs, kn, pe2, pe2_ad, q2, &
&   q2_ad, i1, i2, j, ibeg, iend, jbeg, jend, iv, kord, q_min)
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
    REAL :: pe1_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL :: pe2_ad(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(INOUT) :: q2_ad(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(IN) :: q_min
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_ad(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: q4_ad(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    REAL :: pl_ad, pr_ad, qsum_ad, dp_ad, esl_ad
    INTEGER :: i, k, l, m, k0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp0
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp1
    REAL :: temp_ad11
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL PUSHREALARRAY_ADM(q4, 4*(i2-i1+1)*km)
      CALL SCALAR_PROFILE(qs, q4, dp1, km, i1, i2, iv, kord, q_min)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PPM_PROFILE_FWD(q4, dp1, km, i1, i2, iv, kord)
      CALL PUSHCONTROL1B(0)
    END IF
    DO i=i1,i2
      k0 = 1
      DO 120 k=1,kn
        CALL PUSHINTEGER4(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER4(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
        CALL PUSHINTEGER4(ad_count)
        CALL PUSHCONTROL2B(2)
        GOTO 123
 100    CALL PUSHCONTROL1B(1)
        CALL PUSHINTEGER4(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          k0 = l
          CALL PUSHCONTROL1B(0)
          GOTO 120
        ELSE
! Fractional area...
          CALL PUSHREALARRAY_ADM(qsum)
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          CALL PUSHINTEGER4(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
              qsum = qsum + dp1(i, m)*q4(1, i, m)
              CALL PUSHINTEGER4(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL1B(0)
          CALL PUSHINTEGER4(ad_count0)
          CALL PUSHCONTROL2B(1)
          GOTO 123
 110      CALL PUSHCONTROL1B(1)
          CALL PUSHINTEGER4(ad_count0)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
          CALL PUSHCONTROL2B(0)
        END IF
 123    CALL PUSHCONTROL1B(1)
 120  CONTINUE
    END DO
    dp1_ad = 0.0
    qsum_ad = 0.0
    q4_ad = 0.0
    DO i=i2,i1,-1
      DO k=kn,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          temp_ad0 = 0.5*(pr+pl)*q2_ad(i, j, k)
          temp_ad1 = 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i, l))*q2_ad(i, &
&           j, k)
          temp_ad2 = -(r3*q4(4, i, l)*q2_ad(i, j, k))
          q4_ad(2, i, l) = q4_ad(2, i, l) + q2_ad(i, j, k) - temp_ad0
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad0 - r3*(pr*(pr+pl)+pl&
&           **2)*q2_ad(i, j, k)
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad0
          pr_ad = (2*pr+pl)*temp_ad2 + temp_ad1
          pl_ad = (2*pl+pr)*temp_ad2 + temp_ad1
          q2_ad(i, j, k) = 0.0
          temp_ad3 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad3
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad3
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad3&
&           /dp1(i, l)
        ELSE
          temp1 = pe2(i, k+1) - pe2(i, k)
          temp_ad11 = -(qsum*q2_ad(i, j, k)/temp1**2)
          qsum_ad = qsum_ad + q2_ad(i, j, k)/temp1
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad11
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad11
          q2_ad(i, j, k) = 0.0
          CALL POPCONTROL2B(branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            temp0 = q4(3, i, m) - q4(2, i, m) + q4(4, i, m)*(-(r23*esl)+&
&             1.)
            temp_ad8 = dp*qsum_ad
            temp_ad9 = 0.5*esl*temp_ad8
            q4_ad(2, i, m) = q4_ad(2, i, m) + temp_ad8 - temp_ad9
            esl_ad = 0.5*temp0*temp_ad8 - q4(4, i, m)*r23*temp_ad9
            q4_ad(3, i, m) = q4_ad(3, i, m) + temp_ad9
            q4_ad(4, i, m) = q4_ad(4, i, m) + (1.-r23*esl)*temp_ad9
            temp_ad10 = esl_ad/dp1(i, m)
            dp_ad = temp_ad10 + (q4(2, i, m)+0.5*(esl*temp0))*qsum_ad
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad10/dp1(i, m)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 130
          END IF
          CALL POPINTEGER4(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL1B(branch)
            ELSE
              dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m)*qsum_ad
              q4_ad(1, i, m) = q4_ad(1, i, m) + dp1(i, m)*qsum_ad
            END IF
            CALL POPINTEGER4(m)
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY_ADM(qsum)
          temp = q4(4, i, l) + q4(3, i, l) - q4(2, i, l)
          temp_ad4 = (q4(2, i, l)+0.5*(temp*(pl+1.))-r3*(q4(4, i, l)*(pl&
&           *(pl+1.)+1.)))*qsum_ad
          temp_ad5 = (pe1(i, l+1)-pe2(i, k))*qsum_ad
          temp_ad6 = 0.5*(pl+1.)*temp_ad5
          temp_ad7 = -(r3*q4(4, i, l)*temp_ad5)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + temp_ad4
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad4
          q4_ad(2, i, l) = q4_ad(2, i, l) + temp_ad5 - temp_ad6
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad6 - r3*(pl*(pl+1.)+1.&
&           )*temp_ad5
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad6
          pl_ad = (2*pl+1.)*temp_ad7 + 0.5*temp*temp_ad5
          qsum_ad = 0.0
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 130    CALL POPINTEGER4(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL1B(branch)
          CALL POPINTEGER4(l)
        END DO
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL PPM_PROFILE_BWD(q4, q4_ad, dp1, dp1_ad, km, i1, i2, iv, kord)
    ELSE
      CALL POPREALARRAY_ADM(q4, 4*(i2-i1+1)*km)
      CALL SCALAR_PROFILE_ADM(qs, q4, q4_ad, dp1, dp1_ad, km, i1, i2, iv&
&                       , kord, q_min)
    END IF
    DO k=km,1,-1
      DO i=i2,i1,-1
        q2_ad(i, j, k) = q2_ad(i, j, k) + q4_ad(1, i, k)
        q4_ad(1, i, k) = 0.0
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE MAP_SCALAR_ADM
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
!  Differentiation of map1_ppm in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_c
!ore_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn_
!core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_mo
!d.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynamic
!s_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_m
!od.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod
!.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scalar
!_profile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steep
!z fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_se
!tup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.comp
!ute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_c
! nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver n
!h_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw_
!core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core_
!mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.co
!mpute_divergence_damping_fb sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corners_
!fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circle
!_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 qs q2
!   with respect to varying inputs: pe1 pe2 qs q2
  SUBROUTINE MAP1_PPM_ADM(km, pe1, pe1_ad, qs, qs_ad, kn, pe2, pe2_ad, &
&   q2, q2_ad, i1, i2, j, ibeg, iend, jbeg, jend, iv, kord)
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
    REAL :: qs_ad(i1:i2)
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL :: pe1_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL :: pe2_ad(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(INOUT) :: q2_ad(ibeg:iend, jbeg:jend, kn)
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_ad(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: q4_ad(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    REAL :: pl_ad, pr_ad, qsum_ad, dp_ad, esl_ad
    INTEGER :: i, k, l, m, k0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp0
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp1
    REAL :: temp_ad11
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL PUSHREALARRAY_ADM(q4, 4*(i2-i1+1)*km)
      CALL CS_PROFILE(qs, q4, dp1, km, i1, i2, iv, kord)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PPM_PROFILE_FWD(q4, dp1, km, i1, i2, iv, kord)
      CALL PUSHCONTROL1B(0)
    END IF
    DO i=i1,i2
      k0 = 1
      DO 120 k=1,kn
        CALL PUSHINTEGER4(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER4(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
        CALL PUSHINTEGER4(ad_count)
        CALL PUSHCONTROL2B(2)
        GOTO 123
 100    CALL PUSHCONTROL1B(1)
        CALL PUSHINTEGER4(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          k0 = l
          CALL PUSHCONTROL1B(0)
          GOTO 120
        ELSE
! Fractional area...
          CALL PUSHREALARRAY_ADM(qsum)
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          CALL PUSHINTEGER4(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
              qsum = qsum + dp1(i, m)*q4(1, i, m)
              CALL PUSHINTEGER4(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL1B(0)
          CALL PUSHINTEGER4(ad_count0)
          CALL PUSHCONTROL2B(1)
          GOTO 123
 110      CALL PUSHCONTROL1B(1)
          CALL PUSHINTEGER4(ad_count0)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
          CALL PUSHCONTROL2B(0)
        END IF
 123    CALL PUSHCONTROL1B(1)
 120  CONTINUE
    END DO
    dp1_ad = 0.0
    qsum_ad = 0.0
    q4_ad = 0.0
    DO i=i2,i1,-1
      DO k=kn,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          temp_ad0 = 0.5*(pr+pl)*q2_ad(i, j, k)
          temp_ad1 = 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i, l))*q2_ad(i, &
&           j, k)
          temp_ad2 = -(r3*q4(4, i, l)*q2_ad(i, j, k))
          q4_ad(2, i, l) = q4_ad(2, i, l) + q2_ad(i, j, k) - temp_ad0
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad0 - r3*(pr*(pr+pl)+pl&
&           **2)*q2_ad(i, j, k)
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad0
          pr_ad = (2*pr+pl)*temp_ad2 + temp_ad1
          pl_ad = (2*pl+pr)*temp_ad2 + temp_ad1
          q2_ad(i, j, k) = 0.0
          temp_ad3 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad3
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad3
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad3&
&           /dp1(i, l)
        ELSE
          temp1 = pe2(i, k+1) - pe2(i, k)
          temp_ad11 = -(qsum*q2_ad(i, j, k)/temp1**2)
          qsum_ad = qsum_ad + q2_ad(i, j, k)/temp1
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad11
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad11
          q2_ad(i, j, k) = 0.0
          CALL POPCONTROL2B(branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            temp0 = q4(3, i, m) - q4(2, i, m) + q4(4, i, m)*(-(r23*esl)+&
&             1.)
            temp_ad8 = dp*qsum_ad
            temp_ad9 = 0.5*esl*temp_ad8
            q4_ad(2, i, m) = q4_ad(2, i, m) + temp_ad8 - temp_ad9
            esl_ad = 0.5*temp0*temp_ad8 - q4(4, i, m)*r23*temp_ad9
            q4_ad(3, i, m) = q4_ad(3, i, m) + temp_ad9
            q4_ad(4, i, m) = q4_ad(4, i, m) + (1.-r23*esl)*temp_ad9
            temp_ad10 = esl_ad/dp1(i, m)
            dp_ad = temp_ad10 + (q4(2, i, m)+0.5*(esl*temp0))*qsum_ad
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad10/dp1(i, m)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 130
          END IF
          CALL POPINTEGER4(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL1B(branch)
            ELSE
              dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m)*qsum_ad
              q4_ad(1, i, m) = q4_ad(1, i, m) + dp1(i, m)*qsum_ad
            END IF
            CALL POPINTEGER4(m)
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY_ADM(qsum)
          temp = q4(4, i, l) + q4(3, i, l) - q4(2, i, l)
          temp_ad4 = (q4(2, i, l)+0.5*(temp*(pl+1.))-r3*(q4(4, i, l)*(pl&
&           *(pl+1.)+1.)))*qsum_ad
          temp_ad5 = (pe1(i, l+1)-pe2(i, k))*qsum_ad
          temp_ad6 = 0.5*(pl+1.)*temp_ad5
          temp_ad7 = -(r3*q4(4, i, l)*temp_ad5)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + temp_ad4
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad4
          q4_ad(2, i, l) = q4_ad(2, i, l) + temp_ad5 - temp_ad6
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad6 - r3*(pl*(pl+1.)+1.&
&           )*temp_ad5
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad6
          pl_ad = (2*pl+1.)*temp_ad7 + 0.5*temp*temp_ad5
          qsum_ad = 0.0
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 130    CALL POPINTEGER4(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL1B(branch)
          CALL POPINTEGER4(l)
        END DO
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL PPM_PROFILE_BWD(q4, q4_ad, dp1, dp1_ad, km, i1, i2, iv, kord)
    ELSE
      CALL POPREALARRAY_ADM(q4, 4*(i2-i1+1)*km)
      CALL CS_PROFILE_ADM(qs, qs_ad, q4, q4_ad, dp1, dp1_ad, km, i1, i2&
&                   , iv, kord)
    END IF
    DO k=km,1,-1
      DO i=i2,i1,-1
        q2_ad(i, j, k) = q2_ad(i, j, k) + q4_ad(1, i, k)
        q4_ad(1, i, k) = 0.0
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE MAP1_PPM_ADM
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
!  Differentiation of mapn_tracer in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dy
!n_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c d
!yn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core
!_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dyna
!mics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_util
!s_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_
!mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.sca
!lar_profile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.st
!eepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c
!_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.c
!ompute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solve
!r_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solve
!r nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh 
!sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_co
!re_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod
!.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corne
!rs_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_cir
!cle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 dp2 q1
!   with respect to varying inputs: pe1 pe2 dp2 q1
  SUBROUTINE MAPN_TRACER_ADM(nq, km, pe1, pe1_ad, pe2, pe2_ad, q1, q1_ad&
&   , dp2, dp2_ad, kord, j, i1, i2, isd, ied, jsd, jed, q_min, fill)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! vertical dimension
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: j, nq, i1, i2
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: kord(nq)
! pressure at layer edges
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL :: pe1_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, km+1)
    REAL :: pe2_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
    REAL, INTENT(IN) :: dp2(i1:i2, km)
    REAL :: dp2_ad(i1:i2, km)
    REAL, INTENT(IN) :: q_min
    LOGICAL, INTENT(IN) :: fill
! Field input
    REAL, INTENT(INOUT) :: q1(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: q1_ad(isd:ied, jsd:jed, km, nq)
! !LOCAL VARIABLES:
    REAL :: q4(4, i1:i2, km, nq)
    REAL :: q4_ad(4, i1:i2, km, nq)
! Field output
    REAL :: q2(i1:i2, km, nq)
    REAL :: q2_ad(i1:i2, km, nq)
    REAL :: qsum(nq)
    REAL :: qsum_ad(nq)
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_ad(i1:i2, km)
    REAL :: qs(i1:i2)
    REAL :: pl, pr, dp, esl, fac1, fac2
    REAL :: pl_ad, pr_ad, dp_ad, esl_ad, fac1_ad, fac2_ad
    INTEGER :: i, k, l, m, k0, iq
    INTEGER :: arg1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp0
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3
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
      CALL PUSHREALARRAY_ADM(q4(:, :, :, iq), 4*(i2-i1+1)*km)
      CALL SCALAR_PROFILE(qs, q4(1:4, i1:i2, 1:km, iq), dp1, km, i1, i2&
&                   , 0, kord(iq), q_min)
    END DO
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 130 k=1,km
        CALL PUSHINTEGER4(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER4(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
        CALL PUSHINTEGER4(ad_count)
        CALL PUSHCONTROL2B(2)
        GOTO 120
 100    CALL PUSHCONTROL1B(1)
        CALL PUSHINTEGER4(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          CALL PUSHREALARRAY_ADM(fac1)
          fac1 = pr + pl
          CALL PUSHREALARRAY_ADM(fac2)
          fac2 = r3*(pr*fac1+pl*pl)
          CALL PUSHREALARRAY_ADM(fac1)
          fac1 = 0.5*fac1
          k0 = l
          CALL PUSHCONTROL1B(0)
          GOTO 130
        ELSE
! Fractional area...
          CALL PUSHREALARRAY_ADM(dp)
          dp = pe1(i, l+1) - pe2(i, k)
          CALL PUSHREALARRAY_ADM(fac1)
          fac1 = 1. + pl
          CALL PUSHREALARRAY_ADM(fac2)
          fac2 = r3*(1.+pl*fac1)
          CALL PUSHREALARRAY_ADM(fac1)
          fac1 = 0.5*fac1
          DO iq=1,nq
            CALL PUSHREALARRAY_ADM(qsum(iq))
            qsum(iq) = dp*(q4(2, i, l, iq)+(q4(4, i, l, iq)+q4(3, i, l, &
&             iq)-q4(2, i, l, iq))*fac1-q4(4, i, l, iq)*fac2)
          END DO
          CALL PUSHINTEGER4(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
              DO iq=1,nq
                CALL PUSHREALARRAY_ADM(qsum(iq))
                qsum(iq) = qsum(iq) + dp1(i, m)*q4(1, i, m, iq)
              END DO
              CALL PUSHINTEGER4(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL1B(0)
          CALL PUSHINTEGER4(ad_count0)
          CALL PUSHCONTROL2B(1)
          GOTO 120
 110      CALL PUSHCONTROL1B(1)
          CALL PUSHINTEGER4(ad_count0)
          CALL PUSHREALARRAY_ADM(dp)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          CALL PUSHREALARRAY_ADM(fac1)
          fac1 = 0.5*esl
          CALL PUSHREALARRAY_ADM(fac2)
          fac2 = 1. - r23*esl
          DO iq=1,nq
            CALL PUSHREALARRAY_ADM(qsum(iq))
            qsum(iq) = qsum(iq) + dp*(q4(2, i, m, iq)+fac1*(q4(3, i, m, &
&             iq)-q4(2, i, m, iq)+q4(4, i, m, iq)*fac2))
          END DO
          k0 = m
          CALL PUSHCONTROL2B(0)
        END IF
 120    CALL PUSHCONTROL1B(1)
 130  CONTINUE
    END DO
    q2_ad = 0.0
    DO iq=nq,1,-1
      DO k=km,1,-1
        DO i=i2,i1,-1
          q2_ad(i, k, iq) = q2_ad(i, k, iq) + q1_ad(i, j, k, iq)
          q1_ad(i, j, k, iq) = 0.0
        END DO
      END DO
    END DO
    dp1_ad = 0.0
    qsum_ad = 0.0
    q4_ad = 0.0
    DO i=i2,i1,-1
      DO k=km,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          fac1_ad = 0.0
          fac2_ad = 0.0
          DO iq=nq,1,-1
            temp_ad2 = fac1*q2_ad(i, k, iq)
            q4_ad(2, i, l, iq) = q4_ad(2, i, l, iq) + q2_ad(i, k, iq) - &
&             temp_ad2
            q4_ad(4, i, l, iq) = q4_ad(4, i, l, iq) + temp_ad2 - fac2*&
&             q2_ad(i, k, iq)
            q4_ad(3, i, l, iq) = q4_ad(3, i, l, iq) + temp_ad2
            fac1_ad = fac1_ad + (q4(4, i, l, iq)+q4(3, i, l, iq)-q4(2, i&
&             , l, iq))*q2_ad(i, k, iq)
            fac2_ad = fac2_ad - q4(4, i, l, iq)*q2_ad(i, k, iq)
            q2_ad(i, k, iq) = 0.0
          END DO
          temp_ad0 = r3*fac2_ad
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY_ADM(fac1)
          fac1_ad = pr*temp_ad0 + 0.5*fac1_ad
          CALL POPREALARRAY_ADM(fac2)
          pr_ad = fac1_ad + fac1*temp_ad0
          pl_ad = fac1_ad + 2*pl*temp_ad0
          CALL POPREALARRAY_ADM(fac1)
          temp_ad1 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad1
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad1
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad1&
&           /dp1(i, l)
        ELSE
          DO iq=nq,1,-1
            temp_ad8 = q2_ad(i, k, iq)/dp2(i, k)
            qsum_ad(iq) = qsum_ad(iq) + temp_ad8
            dp2_ad(i, k) = dp2_ad(i, k) - qsum(iq)*temp_ad8/dp2(i, k)
            q2_ad(i, k, iq) = 0.0
          END DO
          CALL POPCONTROL2B(branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            fac1 = 0.5*esl
            fac2 = 1. - r23*esl
            dp_ad = 0.0
            fac1_ad = 0.0
            fac2_ad = 0.0
            DO iq=nq,1,-1
              CALL POPREALARRAY_ADM(qsum(iq))
              temp0 = q4(3, i, m, iq) - q4(2, i, m, iq) + q4(4, i, m, iq&
&               )*fac2
              temp_ad6 = dp*qsum_ad(iq)
              temp_ad7 = fac1*temp_ad6
              dp_ad = dp_ad + (q4(2, i, m, iq)+fac1*temp0)*qsum_ad(iq)
              q4_ad(2, i, m, iq) = q4_ad(2, i, m, iq) + temp_ad6 - &
&               temp_ad7
              fac1_ad = fac1_ad + temp0*temp_ad6
              q4_ad(3, i, m, iq) = q4_ad(3, i, m, iq) + temp_ad7
              q4_ad(4, i, m, iq) = q4_ad(4, i, m, iq) + fac2*temp_ad7
              fac2_ad = fac2_ad + q4(4, i, m, iq)*temp_ad7
            END DO
            CALL POPREALARRAY_ADM(fac2)
            esl_ad = 0.5*fac1_ad - r23*fac2_ad
            CALL POPREALARRAY_ADM(fac1)
            temp_ad5 = esl_ad/dp1(i, m)
            dp_ad = dp_ad + temp_ad5
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad5/dp1(i, m)
            CALL POPREALARRAY_ADM(dp)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 140
          END IF
          CALL POPINTEGER4(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL1B(branch)
            ELSE
              DO iq=nq,1,-1
                CALL POPREALARRAY_ADM(qsum(iq))
                dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m, iq)*qsum_ad(iq&
&                 )
                q4_ad(1, i, m, iq) = q4_ad(1, i, m, iq) + dp1(i, m)*&
&                 qsum_ad(iq)
              END DO
            END IF
            CALL POPINTEGER4(m)
          END DO
          dp_ad = 0.0
          fac1_ad = 0.0
          fac2_ad = 0.0
          DO iq=nq,1,-1
            CALL POPREALARRAY_ADM(qsum(iq))
            temp = q4(4, i, l, iq) + q4(3, i, l, iq) - q4(2, i, l, iq)
            temp_ad3 = dp*qsum_ad(iq)
            temp_ad4 = fac1*temp_ad3
            dp_ad = dp_ad + (q4(2, i, l, iq)+temp*fac1-q4(4, i, l, iq)*&
&             fac2)*qsum_ad(iq)
            q4_ad(2, i, l, iq) = q4_ad(2, i, l, iq) + temp_ad3 - &
&             temp_ad4
            q4_ad(4, i, l, iq) = q4_ad(4, i, l, iq) + temp_ad4 - fac2*&
&             temp_ad3
            q4_ad(3, i, l, iq) = q4_ad(3, i, l, iq) + temp_ad4
            fac1_ad = fac1_ad + temp*temp_ad3
            fac2_ad = fac2_ad - q4(4, i, l, iq)*temp_ad3
            qsum_ad(iq) = 0.0
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY_ADM(fac1)
          fac1_ad = r3*pl*fac2_ad + 0.5*fac1_ad
          CALL POPREALARRAY_ADM(fac2)
          pl_ad = fac1_ad + r3*fac1*fac2_ad
          CALL POPREALARRAY_ADM(fac1)
          CALL POPREALARRAY_ADM(dp)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + dp_ad
          pe2_ad(i, k) = pe2_ad(i, k) - dp_ad
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 140    CALL POPINTEGER4(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL1B(branch)
          CALL POPINTEGER4(l)
        END DO
      END DO
    END DO
    DO iq=nq,1,-1
      CALL POPREALARRAY_ADM(q4(:, :, :, iq), 4*(i2-i1+1)*km)
      CALL SCALAR_PROFILE_ADM(qs, q4(1:4, i1:i2, 1:km, iq), q4_ad(1:4, &
&                       i1:i2, 1:km, iq), dp1, dp1_ad, km, i1, i2, 0, &
&                       kord(iq), q_min)
      DO k=km,1,-1
        DO i=i2,i1,-1
          q1_ad(i, j, k, iq) = q1_ad(i, j, k, iq) + q4_ad(1, i, k, iq)
          q4_ad(1, i, k, iq) = 0.0
        END DO
      END DO
    END DO
    DO k=km,1,-1
      DO i=i2,i1,-1
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE MAPN_TRACER_ADM
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
    INTEGER :: arg1
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
    IF (fill) THEN
      arg1 = i2 - i1 + 1
      CALL FILLZ(arg1, km, nq, q2, dp2)
    END IF
    DO iq=1,nq
!    if (fill) call fillz(i2-i1+1, km, 1, q2(i1,1,iq), dp2)
      DO k=1,km
        DO i=i1,i2
          q1(i, j, k, iq) = q2(i, k, iq)
        END DO
      END DO
    END DO
  END SUBROUTINE MAPN_TRACER
!  Differentiation of map1_q2 in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_co
!re_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn_c
!ore_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_mod
!.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynamics
!_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mo
!d.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.
!map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scalar_
!profile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz
! fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_set
!up fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.compu
!te_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_c 
!nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver nh
!_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw_c
!ore_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core_m
!od.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.com
!pute_divergence_damping_fb sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corners_f
!b tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circle_
!dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 dp2 q1 q2
!   with respect to varying inputs: pe1 pe2 dp2 q1 q2
  SUBROUTINE MAP1_Q2_ADM(km, pe1, pe1_ad, q1, q1_ad, kn, pe2, pe2_ad, q2&
&   , q2_ad, dp2, dp2_ad, i1, i2, iv, kord, j, ibeg, iend, jbeg, jend, &
&   q_min)
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
    REAL :: pe1_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL :: pe2_ad(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! Field input
    REAL, INTENT(IN) :: q1(ibeg:iend, jbeg:jend, km)
    REAL :: q1_ad(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(IN) :: dp2(i1:i2, kn)
    REAL :: dp2_ad(i1:i2, kn)
    REAL, INTENT(IN) :: q_min
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(i1:i2, kn)
    REAL, INTENT(INOUT) :: q2_ad(i1:i2, kn)
! !LOCAL VARIABLES:
    REAL :: qs(i1:i2)
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_ad(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: q4_ad(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    REAL :: pl_ad, pr_ad, qsum_ad, dp_ad, esl_ad
    INTEGER :: i, k, l, m, k0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp0
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q1(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL PUSHREALARRAY_ADM(q4, 4*(i2-i1+1)*km)
      CALL SCALAR_PROFILE(qs, q4, dp1, km, i1, i2, iv, kord, q_min)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PPM_PROFILE_FWD(q4, dp1, km, i1, i2, iv, kord)
      CALL PUSHCONTROL1B(0)
    END IF
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 120 k=1,kn
        CALL PUSHINTEGER4(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER4(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
        CALL PUSHINTEGER4(ad_count)
        CALL PUSHCONTROL2B(2)
        GOTO 123
 100    CALL PUSHCONTROL1B(1)
        CALL PUSHINTEGER4(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          k0 = l
          CALL PUSHCONTROL1B(0)
          GOTO 120
        ELSE
! Fractional area...
          CALL PUSHREALARRAY_ADM(qsum)
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          CALL PUSHINTEGER4(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
              qsum = qsum + dp1(i, m)*q4(1, i, m)
              CALL PUSHINTEGER4(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL1B(0)
          CALL PUSHINTEGER4(ad_count0)
          CALL PUSHCONTROL2B(1)
          GOTO 123
 110      CALL PUSHCONTROL1B(1)
          CALL PUSHINTEGER4(ad_count0)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
          CALL PUSHCONTROL2B(0)
        END IF
 123    CALL PUSHCONTROL1B(1)
 120  CONTINUE
    END DO
    dp1_ad = 0.0
    qsum_ad = 0.0
    q4_ad = 0.0
    DO i=i2,i1,-1
      DO k=kn,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          temp_ad0 = 0.5*(pr+pl)*q2_ad(i, k)
          temp_ad1 = 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i, l))*q2_ad(i, &
&           k)
          temp_ad2 = -(r3*q4(4, i, l)*q2_ad(i, k))
          q4_ad(2, i, l) = q4_ad(2, i, l) + q2_ad(i, k) - temp_ad0
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad0 - r3*(pr*(pr+pl)+pl&
&           **2)*q2_ad(i, k)
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad0
          pr_ad = (2*pr+pl)*temp_ad2 + temp_ad1
          pl_ad = (2*pl+pr)*temp_ad2 + temp_ad1
          q2_ad(i, k) = 0.0
          temp_ad3 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad3
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad3
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad3&
&           /dp1(i, l)
        ELSE
          temp_ad11 = q2_ad(i, k)/dp2(i, k)
          qsum_ad = qsum_ad + temp_ad11
          dp2_ad(i, k) = dp2_ad(i, k) - qsum*temp_ad11/dp2(i, k)
          q2_ad(i, k) = 0.0
          CALL POPCONTROL2B(branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            temp0 = q4(3, i, m) - q4(2, i, m) + q4(4, i, m)*(-(r23*esl)+&
&             1.)
            temp_ad8 = dp*qsum_ad
            temp_ad9 = 0.5*esl*temp_ad8
            q4_ad(2, i, m) = q4_ad(2, i, m) + temp_ad8 - temp_ad9
            esl_ad = 0.5*temp0*temp_ad8 - q4(4, i, m)*r23*temp_ad9
            q4_ad(3, i, m) = q4_ad(3, i, m) + temp_ad9
            q4_ad(4, i, m) = q4_ad(4, i, m) + (1.-r23*esl)*temp_ad9
            temp_ad10 = esl_ad/dp1(i, m)
            dp_ad = temp_ad10 + (q4(2, i, m)+0.5*(esl*temp0))*qsum_ad
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad10/dp1(i, m)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 130
          END IF
          CALL POPINTEGER4(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL1B(branch)
            ELSE
              dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m)*qsum_ad
              q4_ad(1, i, m) = q4_ad(1, i, m) + dp1(i, m)*qsum_ad
            END IF
            CALL POPINTEGER4(m)
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY_ADM(qsum)
          temp = q4(4, i, l) + q4(3, i, l) - q4(2, i, l)
          temp_ad4 = (q4(2, i, l)+0.5*(temp*(pl+1.))-r3*(q4(4, i, l)*(pl&
&           *(pl+1.)+1.)))*qsum_ad
          temp_ad5 = (pe1(i, l+1)-pe2(i, k))*qsum_ad
          temp_ad6 = 0.5*(pl+1.)*temp_ad5
          temp_ad7 = -(r3*q4(4, i, l)*temp_ad5)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + temp_ad4
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad4
          q4_ad(2, i, l) = q4_ad(2, i, l) + temp_ad5 - temp_ad6
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad6 - r3*(pl*(pl+1.)+1.&
&           )*temp_ad5
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad6
          pl_ad = (2*pl+1.)*temp_ad7 + 0.5*temp*temp_ad5
          qsum_ad = 0.0
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 130    CALL POPINTEGER4(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL1B(branch)
          CALL POPINTEGER4(l)
        END DO
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL PPM_PROFILE_BWD(q4, q4_ad, dp1, dp1_ad, km, i1, i2, iv, kord)
    ELSE
      CALL POPREALARRAY_ADM(q4, 4*(i2-i1+1)*km)
      CALL SCALAR_PROFILE_ADM(qs, q4, q4_ad, dp1, dp1_ad, km, i1, i2, iv&
&                       , kord, q_min)
    END IF
    DO k=km,1,-1
      DO i=i2,i1,-1
        q1_ad(i, j, k) = q1_ad(i, j, k) + q4_ad(1, i, k)
        q4_ad(1, i, k) = 0.0
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE MAP1_Q2_ADM
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
!  Differentiation of scalar_profile in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2
! dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_
!c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_c
!ore_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_d
!ynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_u
!tils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_ma
!pz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.
!scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod
!.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.
!d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mo
!d.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_So
!lver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_so
!lver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_
!nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw
!_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_
!mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_co
!rners_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_
!circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: delp a4
!   with respect to varying inputs: delp a4
  SUBROUTINE SCALAR_PROFILE_ADM(qs, a4, a4_ad, delp, delp_ad, km, i1, i2&
&   , iv, kord, qmin)
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
    REAL :: delp_ad(i1:i2, km)
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_ad(4, i1:i2, km)
    REAL, INTENT(IN) :: qmin
!-----------------------------------------------------------------------
    LOGICAL, DIMENSION(i1:i2, km) :: extm, ext6
    REAL :: gam(i1:i2, km)
    REAL :: gam_ad(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: q_ad(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: d4_ad(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: bet_ad, a_bot_ad, grat_ad
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: pmp_1_ad, lac_1_ad, pmp_2_ad, lac_2_ad
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
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp0
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp1
    REAL :: temp2
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: y1_ad
    REAL :: y2_ad
    REAL :: y3_ad
    REAL :: y4_ad
    REAL :: y5_ad
    REAL :: y6_ad
    REAL :: y7_ad
    REAL :: y8_ad
    REAL :: temp_ad15
    REAL :: temp_ad16
    REAL :: temp_ad17
    REAL :: y21_ad
    REAL :: x1_ad
    REAL :: y9_ad
    REAL :: y22_ad
    REAL :: x2_ad
    REAL :: y10_ad
    REAL :: temp_ad18
    REAL :: temp_ad19
    REAL :: y23_ad
    REAL :: x3_ad
    REAL :: y11_ad
    REAL :: y24_ad
    REAL :: x4_ad
    REAL :: y12_ad
    REAL :: temp_ad20
    REAL :: y25_ad
    REAL :: x5_ad
    REAL :: y13_ad
    REAL :: y26_ad
    REAL :: x6_ad
    REAL :: y14_ad
    REAL :: y27_ad
    REAL :: x7_ad
    REAL :: y15_ad
    REAL :: y28_ad
    REAL :: x8_ad
    REAL :: y16_ad
    REAL :: y29_ad
    REAL :: x9_ad
    REAL :: y17_ad
    REAL :: y30_ad
    REAL :: x10_ad
    REAL :: y18_ad
    REAL :: temp_ad21
    REAL :: temp_ad22
    REAL :: temp_ad23
    REAL :: y31_ad
    REAL :: x11_ad
    REAL :: y19_ad
    REAL :: y32_ad
    REAL :: x12_ad
    REAL :: y20_ad
    REAL :: temp_ad24
    REAL :: temp_ad25
    REAL :: temp_ad26
    INTEGER :: branch
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
          CALL PUSHREALARRAY_ADM(bet)
          bet = 2. + grat + grat - gam(i, k)
          CALL PUSHREALARRAY_ADM(q(i, k))
          q(i, k) = (3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))/bet
          gam(i, k+1) = grat/bet
        END DO
      END DO
      DO i=i1,i2
        grat = delp(i, km-1)/delp(i, km)
        CALL PUSHREALARRAY_ADM(q(i, km))
        q(i, km) = (3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-&
&         1))/(2.+grat+grat-gam(i, km))
        CALL PUSHREALARRAY_ADM(q(i, km+1))
        q(i, km+1) = qs(i)
      END DO
      DO k=km-1,1,-1
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(q(i, k))
          q(i, k) = q(i, k) - gam(i, k+1)*q(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL1B(1)
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
          CALL PUSHREALARRAY_ADM(d4(i))
          d4(i) = delp(i, k-1)/delp(i, k)
          CALL PUSHREALARRAY_ADM(bet)
          bet = 2. + d4(i) + d4(i) - gam(i, k-1)
          CALL PUSHREALARRAY_ADM(q(i, k))
          q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
          gam(i, k) = d4(i)/bet
        END DO
      END DO
      DO i=i1,i2
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        CALL PUSHREALARRAY_ADM(q(i, km+1))
        q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(q(i, k))
          q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL1B(0)
    END IF
    IF (kord .GE. 0.) THEN
      abs0 = kord
    ELSE
      abs0 = -kord
    END IF
!----- Perfectly linear scheme --------------------------------
    IF (abs0 .GT. 16) THEN
      q_ad = 0.0
      DO k=km,1,-1
        DO i=i2,i1,-1
          temp_ad14 = 3.*a4_ad(4, i, k)
          a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad14
          a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad14
          a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad14
          a4_ad(4, i, k) = 0.0
          q_ad(i, k+1) = q_ad(i, k+1) + a4_ad(3, i, k)
          a4_ad(3, i, k) = 0.0
          q_ad(i, k) = q_ad(i, k) + a4_ad(2, i, k)
          a4_ad(2, i, k) = 0.0
        END DO
      END DO
      gam_ad = 0.0
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
          CALL PUSHCONTROL1B(0)
        ELSE
          y1 = a4(1, i, 1)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (q(i, 2) .GT. y1) THEN
          CALL PUSHREALARRAY_ADM(q(i, 2))
          q(i, 2) = y1
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREALARRAY_ADM(q(i, 2))
          q(i, 2) = q(i, 2)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (a4(1, i, 1) .GT. a4(1, i, 2)) THEN
          y2 = a4(1, i, 2)
          CALL PUSHCONTROL1B(0)
        ELSE
          y2 = a4(1, i, 1)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (q(i, 2) .LT. y2) THEN
          q(i, 2) = y2
          CALL PUSHCONTROL1B(0)
        ELSE
          q(i, 2) = q(i, 2)
          CALL PUSHCONTROL1B(1)
        END IF
      END DO
      DO k=2,km
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(gam(i, k))
          gam(i, k) = a4(1, i, k) - a4(1, i, k-1)
        END DO
      END DO
! Interior:
      DO k=3,km-1
        DO i=i1,i2
          IF (gam(i, k-1)*gam(i, k+1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y3 = a4(1, i, k)
              CALL PUSHCONTROL1B(0)
            ELSE
              y3 = a4(1, i, k-1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, k) .GT. y3) THEN
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = y3
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = q(i, k)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y4 = a4(1, i, k)
              CALL PUSHCONTROL1B(0)
            ELSE
              y4 = a4(1, i, k-1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, k) .LT. y4) THEN
              q(i, k) = y4
              CALL PUSHCONTROL3B(5)
            ELSE
              q(i, k) = q(i, k)
              CALL PUSHCONTROL3B(6)
            END IF
          ELSE IF (gam(i, k-1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y5 = a4(1, i, k)
              CALL PUSHCONTROL1B(0)
            ELSE
              y5 = a4(1, i, k-1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, k) .LT. y5) THEN
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = y5
              CALL PUSHCONTROL3B(3)
            ELSE
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = q(i, k)
              CALL PUSHCONTROL3B(4)
            END IF
          ELSE
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y6 = a4(1, i, k)
              CALL PUSHCONTROL1B(0)
            ELSE
              y6 = a4(1, i, k-1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, k) .GT. y6) THEN
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = y6
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = q(i, k)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (iv .EQ. 0) THEN
              IF (0. .LT. q(i, k)) THEN
                CALL PUSHCONTROL3B(0)
                q(i, k) = q(i, k)
              ELSE
                q(i, k) = 0.
                CALL PUSHCONTROL3B(2)
              END IF
            ELSE
              CALL PUSHCONTROL3B(1)
            END IF
          END IF
        END DO
      END DO
! Bottom:
      DO i=i1,i2
        IF (a4(1, i, km-1) .LT. a4(1, i, km)) THEN
          y7 = a4(1, i, km)
          CALL PUSHCONTROL1B(0)
        ELSE
          y7 = a4(1, i, km-1)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (q(i, km) .GT. y7) THEN
          CALL PUSHREALARRAY_ADM(q(i, km))
          q(i, km) = y7
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREALARRAY_ADM(q(i, km))
          q(i, km) = q(i, km)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (a4(1, i, km-1) .GT. a4(1, i, km)) THEN
          y8 = a4(1, i, km)
          CALL PUSHCONTROL1B(0)
        ELSE
          y8 = a4(1, i, km-1)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (q(i, km) .LT. y8) THEN
          q(i, km) = y8
          CALL PUSHCONTROL1B(0)
        ELSE
          q(i, km) = q(i, km)
          CALL PUSHCONTROL1B(1)
        END IF
      END DO
      DO k=1,km
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(a4(2, i, k))
          a4(2, i, k) = q(i, k)
          CALL PUSHREALARRAY_ADM(a4(3, i, k))
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
            CALL PUSHREALARRAY_ADM(a4(4, i, k))
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
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
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
            CALL PUSHREALARRAY_ADM(a4(2, i, 1))
            a4(2, i, 1) = a4(2, i, 1)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREALARRAY_ADM(a4(2, i, 1))
            a4(2, i, 1) = 0.
            CALL PUSHCONTROL1B(1)
          END IF
        END DO
        CALL PUSHCONTROL2B(0)
      ELSE IF (iv .EQ. -1) THEN
        DO i=i1,i2
          IF (a4(2, i, 1)*a4(1, i, 1) .LE. 0.) THEN
            CALL PUSHREALARRAY_ADM(a4(2, i, 1))
            a4(2, i, 1) = 0.
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL2B(1)
      ELSE IF (iv .EQ. 2) THEN
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(a4(2, i, 1))
          a4(2, i, 1) = a4(1, i, 1)
          CALL PUSHREALARRAY_ADM(a4(3, i, 1))
          a4(3, i, 1) = a4(1, i, 1)
          CALL PUSHREALARRAY_ADM(a4(4, i, 1))
          a4(4, i, 1) = 0.
        END DO
        CALL PUSHCONTROL2B(2)
      ELSE
        CALL PUSHCONTROL2B(3)
      END IF
      IF (iv .NE. 2) THEN
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(a4(4, i, 1))
          a4(4, i, 1) = 3.*(2.*a4(1, i, 1)-(a4(2, i, 1)+a4(3, i, 1)))
        END DO
        CALL CS_LIMITERS_FWD(im, extm(i1, 1), a4(1, i1, 1), 1)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
! k=2
      DO i=i1,i2
        CALL PUSHREALARRAY_ADM(a4(4, i, 2))
        a4(4, i, 2) = 3.*(2.*a4(1, i, 2)-(a4(2, i, 2)+a4(3, i, 2)))
      END DO
      CALL CS_LIMITERS_FWD(im, extm(i1, 2), a4(1, i1, 2), 2)
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
                CALL PUSHCONTROL2B(0)
              ELSE
                y21 = pmp_1
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (a4(1, i, k) .GT. lac_1) THEN
              y21 = lac_1
              CALL PUSHCONTROL2B(2)
            ELSE
              y21 = a4(1, i, k)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (a4(2, i, k) .LT. y21) THEN
              x1 = y21
              CALL PUSHCONTROL1B(0)
            ELSE
              x1 = a4(2, i, k)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (a4(1, i, k) .LT. pmp_1) THEN
              IF (pmp_1 .LT. lac_1) THEN
                y9 = lac_1
                CALL PUSHCONTROL2B(0)
              ELSE
                y9 = pmp_1
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (a4(1, i, k) .LT. lac_1) THEN
              y9 = lac_1
              CALL PUSHCONTROL2B(2)
            ELSE
              y9 = a4(1, i, k)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (x1 .GT. y9) THEN
              CALL PUSHREALARRAY_ADM(a4(2, i, k))
              a4(2, i, k) = y9
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(a4(2, i, k))
              a4(2, i, k) = x1
              CALL PUSHCONTROL1B(1)
            END IF
! Right edges
            pmp_2 = a4(1, i, k) + 2.*gam(i, k)
            lac_2 = pmp_2 - 1.5*gam(i, k-1)
            IF (a4(1, i, k) .GT. pmp_2) THEN
              IF (pmp_2 .GT. lac_2) THEN
                y22 = lac_2
                CALL PUSHCONTROL2B(0)
              ELSE
                y22 = pmp_2
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (a4(1, i, k) .GT. lac_2) THEN
              y22 = lac_2
              CALL PUSHCONTROL2B(2)
            ELSE
              y22 = a4(1, i, k)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (a4(3, i, k) .LT. y22) THEN
              x2 = y22
              CALL PUSHCONTROL1B(0)
            ELSE
              x2 = a4(3, i, k)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (a4(1, i, k) .LT. pmp_2) THEN
              IF (pmp_2 .LT. lac_2) THEN
                y10 = lac_2
                CALL PUSHCONTROL2B(0)
              ELSE
                y10 = pmp_2
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (a4(1, i, k) .LT. lac_2) THEN
              y10 = lac_2
              CALL PUSHCONTROL2B(2)
            ELSE
              y10 = a4(1, i, k)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (x2 .GT. y10) THEN
              CALL PUSHREALARRAY_ADM(a4(3, i, k))
              a4(3, i, k) = y10
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(a4(3, i, k))
              a4(3, i, k) = x2
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREALARRAY_ADM(a4(4, i, k))
            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
          END DO
          CALL PUSHCONTROL3B(0)
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
                CALL PUSHREALARRAY_ADM(a4(2, i, k))
                a4(2, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(3, i, k))
                a4(3, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(4, i, k))
                a4(4, i, k) = 0.
                CALL PUSHCONTROL3B(4)
              ELSE IF (extm(i, k) .AND. extm(i, k+1)) THEN
! grid-scale 2-delta-z wave detected
                CALL PUSHREALARRAY_ADM(a4(2, i, k))
                a4(2, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(3, i, k))
                a4(3, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(4, i, k))
                a4(4, i, k) = 0.
                CALL PUSHCONTROL3B(3)
              ELSE IF (extm(i, k) .AND. a4(1, i, k) .LT. qmin) THEN
! grid-scale 2-delta-z wave detected
                CALL PUSHREALARRAY_ADM(a4(2, i, k))
                a4(2, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(3, i, k))
                a4(3, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(4, i, k))
                a4(4, i, k) = 0.
                CALL PUSHCONTROL3B(2)
              ELSE
                CALL PUSHREALARRAY_ADM(a4(4, i, k))
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
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y23 = pmp_1
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                    y23 = lac_1
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y23 = a4(1, i, k)
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (a4(2, i, k) .LT. y23) THEN
                    x3 = y23
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    x3 = a4(2, i, k)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_1) THEN
                    IF (pmp_1 .LT. lac_1) THEN
                      y11 = lac_1
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y11 = pmp_1
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                    y11 = lac_1
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y11 = a4(1, i, k)
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (x3 .GT. y11) THEN
                    CALL PUSHREALARRAY_ADM(a4(2, i, k))
                    a4(2, i, k) = y11
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREALARRAY_ADM(a4(2, i, k))
                    a4(2, i, k) = x3
                    CALL PUSHCONTROL1B(1)
                  END IF
                  pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                  lac_2 = pmp_2 - 1.5*gam(i, k-1)
                  IF (a4(1, i, k) .GT. pmp_2) THEN
                    IF (pmp_2 .GT. lac_2) THEN
                      y24 = lac_2
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y24 = pmp_2
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                    y24 = lac_2
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y24 = a4(1, i, k)
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (a4(3, i, k) .LT. y24) THEN
                    x4 = y24
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    x4 = a4(3, i, k)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_2) THEN
                    IF (pmp_2 .LT. lac_2) THEN
                      y12 = lac_2
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y12 = pmp_2
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                    y12 = lac_2
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y12 = a4(1, i, k)
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (x4 .GT. y12) THEN
                    CALL PUSHREALARRAY_ADM(a4(3, i, k))
                    a4(3, i, k) = y12
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREALARRAY_ADM(a4(3, i, k))
                    a4(3, i, k) = x4
                    CALL PUSHCONTROL1B(1)
                  END IF
                  CALL PUSHREALARRAY_ADM(a4(4, i, k))
                  a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i&
&                   , k)))
                  CALL PUSHCONTROL3B(1)
                ELSE
                  CALL PUSHCONTROL3B(0)
                END IF
              END IF
            END DO
            CALL PUSHCONTROL3B(1)
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
                    CALL PUSHREALARRAY_ADM(a4(2, i, k))
                    a4(2, i, k) = a4(1, i, k)
                    CALL PUSHREALARRAY_ADM(a4(3, i, k))
                    a4(3, i, k) = a4(1, i, k)
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
                    a4(4, i, k) = 0.
                    CALL PUSHCONTROL2B(3)
                  ELSE
! True local extremum
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                    CALL PUSHCONTROL2B(2)
                  END IF
                ELSE
! not a local extremum
                  CALL PUSHREALARRAY_ADM(a4(4, i, k))
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
                        CALL PUSHCONTROL2B(0)
                      ELSE
                        y25 = pmp_1
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                      y25 = lac_1
                      CALL PUSHCONTROL2B(2)
                    ELSE
                      y25 = a4(1, i, k)
                      CALL PUSHCONTROL2B(3)
                    END IF
                    IF (a4(2, i, k) .LT. y25) THEN
                      x5 = y25
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      x5 = a4(2, i, k)
                      CALL PUSHCONTROL1B(1)
                    END IF
                    IF (a4(1, i, k) .LT. pmp_1) THEN
                      IF (pmp_1 .LT. lac_1) THEN
                        y13 = lac_1
                        CALL PUSHCONTROL2B(0)
                      ELSE
                        y13 = pmp_1
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                      y13 = lac_1
                      CALL PUSHCONTROL2B(2)
                    ELSE
                      y13 = a4(1, i, k)
                      CALL PUSHCONTROL2B(3)
                    END IF
                    IF (x5 .GT. y13) THEN
                      CALL PUSHREALARRAY_ADM(a4(2, i, k))
                      a4(2, i, k) = y13
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      CALL PUSHREALARRAY_ADM(a4(2, i, k))
                      a4(2, i, k) = x5
                      CALL PUSHCONTROL1B(1)
                    END IF
                    pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                    lac_2 = pmp_2 - 1.5*gam(i, k-1)
                    IF (a4(1, i, k) .GT. pmp_2) THEN
                      IF (pmp_2 .GT. lac_2) THEN
                        y26 = lac_2
                        CALL PUSHCONTROL2B(0)
                      ELSE
                        y26 = pmp_2
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                      y26 = lac_2
                      CALL PUSHCONTROL2B(2)
                    ELSE
                      y26 = a4(1, i, k)
                      CALL PUSHCONTROL2B(3)
                    END IF
                    IF (a4(3, i, k) .LT. y26) THEN
                      x6 = y26
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      x6 = a4(3, i, k)
                      CALL PUSHCONTROL1B(1)
                    END IF
                    IF (a4(1, i, k) .LT. pmp_2) THEN
                      IF (pmp_2 .LT. lac_2) THEN
                        y14 = lac_2
                        CALL PUSHCONTROL2B(0)
                      ELSE
                        y14 = pmp_2
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                      y14 = lac_2
                      CALL PUSHCONTROL2B(2)
                    ELSE
                      y14 = a4(1, i, k)
                      CALL PUSHCONTROL2B(3)
                    END IF
                    IF (x6 .GT. y14) THEN
                      CALL PUSHREALARRAY_ADM(a4(3, i, k))
                      a4(3, i, k) = y14
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      CALL PUSHREALARRAY_ADM(a4(3, i, k))
                      a4(3, i, k) = x6
                      CALL PUSHCONTROL1B(1)
                    END IF
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                    CALL PUSHCONTROL2B(1)
                  ELSE
                    CALL PUSHCONTROL2B(0)
                  END IF
                END IF
              END DO
              CALL PUSHCONTROL3B(2)
            ELSE
              IF (kord .GE. 0.) THEN
                abs8 = kord
              ELSE
                abs8 = -kord
              END IF
              IF (abs8 .EQ. 12) THEN
                DO i=i1,i2
                  IF (extm(i, k)) THEN
                    CALL PUSHREALARRAY_ADM(a4(2, i, k))
                    a4(2, i, k) = a4(1, i, k)
                    CALL PUSHREALARRAY_ADM(a4(3, i, k))
                    a4(3, i, k) = a4(1, i, k)
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
                    a4(4, i, k) = 0.
                    CALL PUSHCONTROL2B(2)
                  ELSE
! not a local extremum
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
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
                          CALL PUSHCONTROL2B(0)
                        ELSE
                          y27 = pmp_1
                          CALL PUSHCONTROL2B(1)
                        END IF
                      ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                        y27 = lac_1
                        CALL PUSHCONTROL2B(2)
                      ELSE
                        y27 = a4(1, i, k)
                        CALL PUSHCONTROL2B(3)
                      END IF
                      IF (a4(2, i, k) .LT. y27) THEN
                        x7 = y27
                        CALL PUSHCONTROL1B(0)
                      ELSE
                        x7 = a4(2, i, k)
                        CALL PUSHCONTROL1B(1)
                      END IF
                      IF (a4(1, i, k) .LT. pmp_1) THEN
                        IF (pmp_1 .LT. lac_1) THEN
                          y15 = lac_1
                          CALL PUSHCONTROL2B(0)
                        ELSE
                          y15 = pmp_1
                          CALL PUSHCONTROL2B(1)
                        END IF
                      ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                        y15 = lac_1
                        CALL PUSHCONTROL2B(2)
                      ELSE
                        y15 = a4(1, i, k)
                        CALL PUSHCONTROL2B(3)
                      END IF
                      IF (x7 .GT. y15) THEN
                        CALL PUSHREALARRAY_ADM(a4(2, i, k))
                        a4(2, i, k) = y15
                        CALL PUSHCONTROL1B(0)
                      ELSE
                        CALL PUSHREALARRAY_ADM(a4(2, i, k))
                        a4(2, i, k) = x7
                        CALL PUSHCONTROL1B(1)
                      END IF
                      pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                      lac_2 = pmp_2 - 1.5*gam(i, k-1)
                      IF (a4(1, i, k) .GT. pmp_2) THEN
                        IF (pmp_2 .GT. lac_2) THEN
                          y28 = lac_2
                          CALL PUSHCONTROL2B(0)
                        ELSE
                          y28 = pmp_2
                          CALL PUSHCONTROL2B(1)
                        END IF
                      ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                        y28 = lac_2
                        CALL PUSHCONTROL2B(2)
                      ELSE
                        y28 = a4(1, i, k)
                        CALL PUSHCONTROL2B(3)
                      END IF
                      IF (a4(3, i, k) .LT. y28) THEN
                        x8 = y28
                        CALL PUSHCONTROL1B(0)
                      ELSE
                        x8 = a4(3, i, k)
                        CALL PUSHCONTROL1B(1)
                      END IF
                      IF (a4(1, i, k) .LT. pmp_2) THEN
                        IF (pmp_2 .LT. lac_2) THEN
                          y16 = lac_2
                          CALL PUSHCONTROL2B(0)
                        ELSE
                          y16 = pmp_2
                          CALL PUSHCONTROL2B(1)
                        END IF
                      ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                        y16 = lac_2
                        CALL PUSHCONTROL2B(2)
                      ELSE
                        y16 = a4(1, i, k)
                        CALL PUSHCONTROL2B(3)
                      END IF
                      IF (x8 .GT. y16) THEN
                        CALL PUSHREALARRAY_ADM(a4(3, i, k))
                        a4(3, i, k) = y16
                        CALL PUSHCONTROL1B(0)
                      ELSE
                        CALL PUSHREALARRAY_ADM(a4(3, i, k))
                        a4(3, i, k) = x8
                        CALL PUSHCONTROL1B(1)
                      END IF
                      CALL PUSHREALARRAY_ADM(a4(4, i, k))
                      a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(&
&                       3, i, k))
                      CALL PUSHCONTROL2B(1)
                    ELSE
                      CALL PUSHCONTROL2B(0)
                    END IF
                  END IF
                END DO
                CALL PUSHCONTROL3B(3)
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
                        CALL PUSHREALARRAY_ADM(a4(2, i, k))
                        a4(2, i, k) = a4(1, i, k)
                        CALL PUSHREALARRAY_ADM(a4(3, i, k))
                        a4(3, i, k) = a4(1, i, k)
                        CALL PUSHREALARRAY_ADM(a4(4, i, k))
                        a4(4, i, k) = 0.
                        CALL PUSHCONTROL2B(2)
                      ELSE
! Left  edges
                        pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                        lac_1 = pmp_1 + 1.5*gam(i, k+2)
                        IF (a4(1, i, k) .GT. pmp_1) THEN
                          IF (pmp_1 .GT. lac_1) THEN
                            y29 = lac_1
                            CALL PUSHCONTROL2B(0)
                          ELSE
                            y29 = pmp_1
                            CALL PUSHCONTROL2B(1)
                          END IF
                        ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                          y29 = lac_1
                          CALL PUSHCONTROL2B(2)
                        ELSE
                          y29 = a4(1, i, k)
                          CALL PUSHCONTROL2B(3)
                        END IF
                        IF (a4(2, i, k) .LT. y29) THEN
                          x9 = y29
                          CALL PUSHCONTROL1B(0)
                        ELSE
                          x9 = a4(2, i, k)
                          CALL PUSHCONTROL1B(1)
                        END IF
                        IF (a4(1, i, k) .LT. pmp_1) THEN
                          IF (pmp_1 .LT. lac_1) THEN
                            y17 = lac_1
                            CALL PUSHCONTROL2B(0)
                          ELSE
                            y17 = pmp_1
                            CALL PUSHCONTROL2B(1)
                          END IF
                        ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                          y17 = lac_1
                          CALL PUSHCONTROL2B(2)
                        ELSE
                          y17 = a4(1, i, k)
                          CALL PUSHCONTROL2B(3)
                        END IF
                        IF (x9 .GT. y17) THEN
                          CALL PUSHREALARRAY_ADM(a4(2, i, k))
                          a4(2, i, k) = y17
                          CALL PUSHCONTROL1B(0)
                        ELSE
                          CALL PUSHREALARRAY_ADM(a4(2, i, k))
                          a4(2, i, k) = x9
                          CALL PUSHCONTROL1B(1)
                        END IF
! Right edges
                        pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                        lac_2 = pmp_2 - 1.5*gam(i, k-1)
                        IF (a4(1, i, k) .GT. pmp_2) THEN
                          IF (pmp_2 .GT. lac_2) THEN
                            y30 = lac_2
                            CALL PUSHCONTROL2B(0)
                          ELSE
                            y30 = pmp_2
                            CALL PUSHCONTROL2B(1)
                          END IF
                        ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                          y30 = lac_2
                          CALL PUSHCONTROL2B(2)
                        ELSE
                          y30 = a4(1, i, k)
                          CALL PUSHCONTROL2B(3)
                        END IF
                        IF (a4(3, i, k) .LT. y30) THEN
                          x10 = y30
                          CALL PUSHCONTROL1B(0)
                        ELSE
                          x10 = a4(3, i, k)
                          CALL PUSHCONTROL1B(1)
                        END IF
                        IF (a4(1, i, k) .LT. pmp_2) THEN
                          IF (pmp_2 .LT. lac_2) THEN
                            y18 = lac_2
                            CALL PUSHCONTROL2B(0)
                          ELSE
                            y18 = pmp_2
                            CALL PUSHCONTROL2B(1)
                          END IF
                        ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                          y18 = lac_2
                          CALL PUSHCONTROL2B(2)
                        ELSE
                          y18 = a4(1, i, k)
                          CALL PUSHCONTROL2B(3)
                        END IF
                        IF (x10 .GT. y18) THEN
                          CALL PUSHREALARRAY_ADM(a4(3, i, k))
                          a4(3, i, k) = y18
                          CALL PUSHCONTROL1B(0)
                        ELSE
                          CALL PUSHREALARRAY_ADM(a4(3, i, k))
                          a4(3, i, k) = x10
                          CALL PUSHCONTROL1B(1)
                        END IF
                        CALL PUSHREALARRAY_ADM(a4(4, i, k))
                        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4&
&                         (3, i, k)))
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE
                      CALL PUSHREALARRAY_ADM(a4(4, i, k))
                      a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3&
&                       , i, k)))
                      CALL PUSHCONTROL2B(0)
                    END IF
                  END DO
                  CALL PUSHCONTROL3B(4)
                ELSE
                  IF (kord .GE. 0.) THEN
                    abs11 = kord
                  ELSE
                    abs11 = -kord
                  END IF
                  IF (abs11 .EQ. 14) THEN
                    DO i=i1,i2
                      CALL PUSHREALARRAY_ADM(a4(4, i, k))
                      a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3&
&                       , i, k)))
                    END DO
                    CALL PUSHCONTROL3B(5)
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
                                CALL PUSHCONTROL2B(0)
                              ELSE
                                y31 = pmp_1
                                CALL PUSHCONTROL2B(1)
                              END IF
                            ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                              y31 = lac_1
                              CALL PUSHCONTROL2B(2)
                            ELSE
                              y31 = a4(1, i, k)
                              CALL PUSHCONTROL2B(3)
                            END IF
                            IF (a4(2, i, k) .LT. y31) THEN
                              x11 = y31
                              CALL PUSHCONTROL1B(0)
                            ELSE
                              x11 = a4(2, i, k)
                              CALL PUSHCONTROL1B(1)
                            END IF
                            IF (a4(1, i, k) .LT. pmp_1) THEN
                              IF (pmp_1 .LT. lac_1) THEN
                                y19 = lac_1
                                CALL PUSHCONTROL2B(0)
                              ELSE
                                y19 = pmp_1
                                CALL PUSHCONTROL2B(1)
                              END IF
                            ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                              y19 = lac_1
                              CALL PUSHCONTROL2B(2)
                            ELSE
                              y19 = a4(1, i, k)
                              CALL PUSHCONTROL2B(3)
                            END IF
                            IF (x11 .GT. y19) THEN
                              CALL PUSHREALARRAY_ADM(a4(2, i, k))
                              a4(2, i, k) = y19
                              CALL PUSHCONTROL1B(0)
                            ELSE
                              CALL PUSHREALARRAY_ADM(a4(2, i, k))
                              a4(2, i, k) = x11
                              CALL PUSHCONTROL1B(1)
                            END IF
! Right edges
                            pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                            lac_2 = pmp_2 - 1.5*gam(i, k-1)
                            IF (a4(1, i, k) .GT. pmp_2) THEN
                              IF (pmp_2 .GT. lac_2) THEN
                                y32 = lac_2
                                CALL PUSHCONTROL2B(0)
                              ELSE
                                y32 = pmp_2
                                CALL PUSHCONTROL2B(1)
                              END IF
                            ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                              y32 = lac_2
                              CALL PUSHCONTROL2B(2)
                            ELSE
                              y32 = a4(1, i, k)
                              CALL PUSHCONTROL2B(3)
                            END IF
                            IF (a4(3, i, k) .LT. y32) THEN
                              x12 = y32
                              CALL PUSHCONTROL1B(0)
                            ELSE
                              x12 = a4(3, i, k)
                              CALL PUSHCONTROL1B(1)
                            END IF
                            IF (a4(1, i, k) .LT. pmp_2) THEN
                              IF (pmp_2 .LT. lac_2) THEN
                                y20 = lac_2
                                CALL PUSHCONTROL2B(0)
                              ELSE
                                y20 = pmp_2
                                CALL PUSHCONTROL2B(1)
                              END IF
                            ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                              y20 = lac_2
                              CALL PUSHCONTROL2B(2)
                            ELSE
                              y20 = a4(1, i, k)
                              CALL PUSHCONTROL2B(3)
                            END IF
                            IF (x12 .GT. y20) THEN
                              CALL PUSHREALARRAY_ADM(a4(3, i, k))
                              a4(3, i, k) = y20
                              CALL PUSHCONTROL1B(0)
                            ELSE
                              CALL PUSHREALARRAY_ADM(a4(3, i, k))
                              a4(3, i, k) = x12
                              CALL PUSHCONTROL1B(1)
                            END IF
                            CALL PUSHREALARRAY_ADM(a4(4, i, k))
                            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k&
&                             )+a4(3, i, k)))
                            CALL PUSHCONTROL2B(2)
                          ELSE
                            CALL PUSHCONTROL2B(1)
                          END IF
                        ELSE
                          CALL PUSHCONTROL2B(0)
                        END IF
                      END DO
                      CALL PUSHCONTROL3B(6)
                    ELSE
! kord = 11, 13
                      DO i=i1,i2
                        IF (extm(i, k) .AND. ((extm(i, k-1) .OR. extm(i&
&                           , k+1)) .OR. a4(1, i, k) .LT. qmin)) THEN
! Noisy region:
                          CALL PUSHREALARRAY_ADM(a4(2, i, k))
                          a4(2, i, k) = a4(1, i, k)
                          CALL PUSHREALARRAY_ADM(a4(3, i, k))
                          a4(3, i, k) = a4(1, i, k)
                          CALL PUSHREALARRAY_ADM(a4(4, i, k))
                          a4(4, i, k) = 0.
                          CALL PUSHCONTROL1B(1)
                        ELSE
                          CALL PUSHREALARRAY_ADM(a4(4, i, k))
                          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+&
&                           a4(3, i, k)))
                          CALL PUSHCONTROL1B(0)
                        END IF
                      END DO
                      CALL PUSHCONTROL3B(7)
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
! Additional constraint to ensure positivity
        IF (iv .EQ. 0) THEN
          CALL CS_LIMITERS_FWD(im, extm(i1, k), a4(1, i1, k), 0)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
! k-loop
!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
      IF (iv .EQ. 0) THEN
        DO i=i1,i2
          IF (0. .LT. a4(3, i, km)) THEN
            CALL PUSHREALARRAY_ADM(a4(3, i, km))
            a4(3, i, km) = a4(3, i, km)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREALARRAY_ADM(a4(3, i, km))
            a4(3, i, km) = 0.
            CALL PUSHCONTROL1B(1)
          END IF
        END DO
        CALL PUSHCONTROL2B(2)
      ELSE IF (iv .EQ. -1) THEN
        DO i=i1,i2
          IF (a4(3, i, km)*a4(1, i, km) .LE. 0.) THEN
            CALL PUSHREALARRAY_ADM(a4(3, i, km))
            a4(3, i, km) = 0.
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
      DO k=km-1,km
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(a4(4, i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
        IF (k .EQ. km - 1) THEN
          CALL CS_LIMITERS_FWD(im, extm(i1, k), a4(1, i1, k), 2)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (k .EQ. km) THEN
          CALL CS_LIMITERS_FWD(im, extm(i1, k), a4(1, i1, k), 1)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
      DO k=km,km-1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) CALL CS_LIMITERS_BWD(im, extm(i1, k), a4(1, &
&                                         i1, k), a4_ad(1, i1, k), 1)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) CALL CS_LIMITERS_BWD(im, extm(i1, k), a4(1, &
&                                         i1, k), a4_ad(1, i1, k), 2)
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(a4(4, i, k))
          temp_ad26 = 3.*a4_ad(4, i, k)
          a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad26
          a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad26
          a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad26
          a4_ad(4, i, k) = 0.0
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .NE. 0) THEN
        IF (branch .EQ. 1) THEN
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) THEN
              CALL POPREALARRAY_ADM(a4(3, i, km))
              a4_ad(3, i, km) = 0.0
            END IF
          END DO
        ELSE
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(a4(3, i, km))
            ELSE
              CALL POPREALARRAY_ADM(a4(3, i, km))
              a4_ad(3, i, km) = 0.0
            END IF
          END DO
        END IF
      END IF
      gam_ad = 0.0
      DO k=km-2,3,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) CALL CS_LIMITERS_BWD(im, extm(i1, k), a4(1, &
&                                         i1, k), a4_ad(1, i1, k), 0)
        CALL POPCONTROL3B(branch)
        IF (branch .LT. 4) THEN
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              DO i=i2,i1,-1
                CALL POPREALARRAY_ADM(a4(4, i, k))
                temp_ad18 = 3.*a4_ad(4, i, k)
                a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad18
                a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad18
                a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad18
                a4_ad(4, i, k) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  y10_ad = a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  x2_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  x2_ad = a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  y10_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = y10_ad
                    pmp_2_ad = 0.0
                  ELSE
                    pmp_2_ad = y10_ad
                    lac_2_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_2_ad = y10_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y10_ad
                    lac_2_ad = 0.0
                  END IF
                  pmp_2_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y22_ad = x2_ad
                ELSE
                  a4_ad(3, i, k) = a4_ad(3, i, k) + x2_ad
                  y22_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = lac_2_ad + y22_ad
                  ELSE
                    pmp_2_ad = pmp_2_ad + y22_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_2_ad = lac_2_ad + y22_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y22_ad
                END IF
                pmp_2_ad = pmp_2_ad + lac_2_ad
                gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  y9_ad = a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  x1_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  x1_ad = a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  y9_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = y9_ad
                    pmp_1_ad = 0.0
                  ELSE
                    pmp_1_ad = y9_ad
                    lac_1_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_1_ad = y9_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y9_ad
                    lac_1_ad = 0.0
                  END IF
                  pmp_1_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y21_ad = x1_ad
                ELSE
                  a4_ad(2, i, k) = a4_ad(2, i, k) + x1_ad
                  y21_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = lac_1_ad + y21_ad
                  ELSE
                    pmp_1_ad = pmp_1_ad + y21_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_1_ad = lac_1_ad + y21_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y21_ad
                END IF
                pmp_1_ad = pmp_1_ad + lac_1_ad
                gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
              END DO
            ELSE
              DO i=i2,i1,-1
                CALL POPCONTROL3B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .NE. 0) THEN
                    CALL POPREALARRAY_ADM(a4(4, i, k))
                    temp_ad20 = 3.*a4_ad(4, i, k)
                    a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad20
                    a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad20
                    a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad20
                    a4_ad(4, i, k) = 0.0
                    CALL POPCONTROL1B(branch)
                    IF (branch .EQ. 0) THEN
                      CALL POPREALARRAY_ADM(a4(3, i, k))
                      y12_ad = a4_ad(3, i, k)
                      a4_ad(3, i, k) = 0.0
                      x4_ad = 0.0
                    ELSE
                      CALL POPREALARRAY_ADM(a4(3, i, k))
                      x4_ad = a4_ad(3, i, k)
                      a4_ad(3, i, k) = 0.0
                      y12_ad = 0.0
                    END IF
                    CALL POPCONTROL2B(branch)
                    IF (branch .LT. 2) THEN
                      IF (branch .EQ. 0) THEN
                        lac_2_ad = y12_ad
                        pmp_2_ad = 0.0
                      ELSE
                        pmp_2_ad = y12_ad
                        lac_2_ad = 0.0
                      END IF
                    ELSE
                      IF (branch .EQ. 2) THEN
                        lac_2_ad = y12_ad
                      ELSE
                        a4_ad(1, i, k) = a4_ad(1, i, k) + y12_ad
                        lac_2_ad = 0.0
                      END IF
                      pmp_2_ad = 0.0
                    END IF
                    CALL POPCONTROL1B(branch)
                    IF (branch .EQ. 0) THEN
                      y24_ad = x4_ad
                    ELSE
                      a4_ad(3, i, k) = a4_ad(3, i, k) + x4_ad
                      y24_ad = 0.0
                    END IF
                    CALL POPCONTROL2B(branch)
                    IF (branch .LT. 2) THEN
                      IF (branch .EQ. 0) THEN
                        lac_2_ad = lac_2_ad + y24_ad
                      ELSE
                        pmp_2_ad = pmp_2_ad + y24_ad
                      END IF
                    ELSE IF (branch .EQ. 2) THEN
                      lac_2_ad = lac_2_ad + y24_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y24_ad
                    END IF
                    pmp_2_ad = pmp_2_ad + lac_2_ad
                    gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                    a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                    gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                    CALL POPCONTROL1B(branch)
                    IF (branch .EQ. 0) THEN
                      CALL POPREALARRAY_ADM(a4(2, i, k))
                      y11_ad = a4_ad(2, i, k)
                      a4_ad(2, i, k) = 0.0
                      x3_ad = 0.0
                    ELSE
                      CALL POPREALARRAY_ADM(a4(2, i, k))
                      x3_ad = a4_ad(2, i, k)
                      a4_ad(2, i, k) = 0.0
                      y11_ad = 0.0
                    END IF
                    CALL POPCONTROL2B(branch)
                    IF (branch .LT. 2) THEN
                      IF (branch .EQ. 0) THEN
                        lac_1_ad = y11_ad
                        pmp_1_ad = 0.0
                      ELSE
                        pmp_1_ad = y11_ad
                        lac_1_ad = 0.0
                      END IF
                    ELSE
                      IF (branch .EQ. 2) THEN
                        lac_1_ad = y11_ad
                      ELSE
                        a4_ad(1, i, k) = a4_ad(1, i, k) + y11_ad
                        lac_1_ad = 0.0
                      END IF
                      pmp_1_ad = 0.0
                    END IF
                    CALL POPCONTROL1B(branch)
                    IF (branch .EQ. 0) THEN
                      y23_ad = x3_ad
                    ELSE
                      a4_ad(2, i, k) = a4_ad(2, i, k) + x3_ad
                      y23_ad = 0.0
                    END IF
                    CALL POPCONTROL2B(branch)
                    IF (branch .LT. 2) THEN
                      IF (branch .EQ. 0) THEN
                        lac_1_ad = lac_1_ad + y23_ad
                      ELSE
                        pmp_1_ad = pmp_1_ad + y23_ad
                      END IF
                    ELSE IF (branch .EQ. 2) THEN
                      lac_1_ad = lac_1_ad + y23_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y23_ad
                    END IF
                    pmp_1_ad = pmp_1_ad + lac_1_ad
                    gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                    a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                    gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
                  END IF
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  temp_ad19 = 3.*a4_ad(4, i, k)
                  a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad19
                  a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad19
                  a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad19
                  a4_ad(4, i, k) = 0.0
                ELSE IF (branch .EQ. 2) THEN
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(4, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                ELSE IF (branch .EQ. 3) THEN
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(4, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(4, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                END IF
              END DO
            END IF
          ELSE IF (branch .EQ. 2) THEN
            DO i=i2,i1,-1
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .NE. 0) THEN
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                  a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(4, i, k) = 0.0
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    y14_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    x6_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    x6_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    y14_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = y14_ad
                      pmp_2_ad = 0.0
                    ELSE
                      pmp_2_ad = y14_ad
                      lac_2_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_2_ad = y14_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y14_ad
                      lac_2_ad = 0.0
                    END IF
                    pmp_2_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y26_ad = x6_ad
                  ELSE
                    a4_ad(3, i, k) = a4_ad(3, i, k) + x6_ad
                    y26_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = lac_2_ad + y26_ad
                    ELSE
                      pmp_2_ad = pmp_2_ad + y26_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_2_ad = lac_2_ad + y26_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y26_ad
                  END IF
                  pmp_2_ad = pmp_2_ad + lac_2_ad
                  gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                  gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    y13_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    x5_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    x5_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    y13_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = y13_ad
                      pmp_1_ad = 0.0
                    ELSE
                      pmp_1_ad = y13_ad
                      lac_1_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_1_ad = y13_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y13_ad
                      lac_1_ad = 0.0
                    END IF
                    pmp_1_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y25_ad = x5_ad
                  ELSE
                    a4_ad(2, i, k) = a4_ad(2, i, k) + x5_ad
                    y25_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = lac_1_ad + y25_ad
                    ELSE
                      pmp_1_ad = pmp_1_ad + y25_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_1_ad = lac_1_ad + y25_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y25_ad
                  END IF
                  pmp_1_ad = pmp_1_ad + lac_1_ad
                  gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                  gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
                END IF
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(4, i, k) = 0.0
              ELSE IF (branch .EQ. 2) THEN
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(4, i, k) = 0.0
              ELSE
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(4, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(3, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                a4_ad(3, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(2, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                a4_ad(2, i, k) = 0.0
              END IF
            END DO
          ELSE
            DO 100 i=i2,i1,-1
              CALL POPCONTROL2B(branch)
              IF (branch .NE. 0) THEN
                IF (branch .EQ. 1) THEN
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                  a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(4, i, k) = 0.0
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    y16_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    x8_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    x8_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    y16_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = y16_ad
                      pmp_2_ad = 0.0
                    ELSE
                      pmp_2_ad = y16_ad
                      lac_2_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_2_ad = y16_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y16_ad
                      lac_2_ad = 0.0
                    END IF
                    pmp_2_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y28_ad = x8_ad
                  ELSE
                    a4_ad(3, i, k) = a4_ad(3, i, k) + x8_ad
                    y28_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = lac_2_ad + y28_ad
                    ELSE
                      pmp_2_ad = pmp_2_ad + y28_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_2_ad = lac_2_ad + y28_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y28_ad
                  END IF
                  pmp_2_ad = pmp_2_ad + lac_2_ad
                  gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                  gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    y15_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    x7_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    x7_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    y15_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = y15_ad
                      pmp_1_ad = 0.0
                    ELSE
                      pmp_1_ad = y15_ad
                      lac_1_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_1_ad = y15_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y15_ad
                      lac_1_ad = 0.0
                    END IF
                    pmp_1_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y27_ad = x7_ad
                  ELSE
                    a4_ad(2, i, k) = a4_ad(2, i, k) + x7_ad
                    y27_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = lac_1_ad + y27_ad
                    ELSE
                      pmp_1_ad = pmp_1_ad + y27_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_1_ad = lac_1_ad + y27_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y27_ad
                  END IF
                  pmp_1_ad = pmp_1_ad + lac_1_ad
                  gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                  gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
                ELSE
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(4, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  GOTO 100
                END IF
              END IF
              CALL POPREALARRAY_ADM(a4(4, i, k))
              a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
              a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
              a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
              a4_ad(4, i, k) = 0.0
 100        CONTINUE
          END IF
        ELSE IF (branch .LT. 6) THEN
          IF (branch .EQ. 4) THEN
            DO i=i2,i1,-1
              CALL POPCONTROL2B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(a4(4, i, k))
                temp_ad22 = 3.*a4_ad(4, i, k)
                a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad22
                a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad22
                a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad22
                a4_ad(4, i, k) = 0.0
              ELSE IF (branch .EQ. 1) THEN
                CALL POPREALARRAY_ADM(a4(4, i, k))
                temp_ad21 = 3.*a4_ad(4, i, k)
                a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad21
                a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad21
                a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad21
                a4_ad(4, i, k) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  y18_ad = a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  x10_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  x10_ad = a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  y18_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = y18_ad
                    pmp_2_ad = 0.0
                  ELSE
                    pmp_2_ad = y18_ad
                    lac_2_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_2_ad = y18_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y18_ad
                    lac_2_ad = 0.0
                  END IF
                  pmp_2_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y30_ad = x10_ad
                ELSE
                  a4_ad(3, i, k) = a4_ad(3, i, k) + x10_ad
                  y30_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = lac_2_ad + y30_ad
                  ELSE
                    pmp_2_ad = pmp_2_ad + y30_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_2_ad = lac_2_ad + y30_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y30_ad
                END IF
                pmp_2_ad = pmp_2_ad + lac_2_ad
                gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  y17_ad = a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  x9_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  x9_ad = a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  y17_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = y17_ad
                    pmp_1_ad = 0.0
                  ELSE
                    pmp_1_ad = y17_ad
                    lac_1_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_1_ad = y17_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y17_ad
                    lac_1_ad = 0.0
                  END IF
                  pmp_1_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y29_ad = x9_ad
                ELSE
                  a4_ad(2, i, k) = a4_ad(2, i, k) + x9_ad
                  y29_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = lac_1_ad + y29_ad
                  ELSE
                    pmp_1_ad = pmp_1_ad + y29_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_1_ad = lac_1_ad + y29_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y29_ad
                END IF
                pmp_1_ad = pmp_1_ad + lac_1_ad
                gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
              ELSE
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(4, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(3, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                a4_ad(3, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(2, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                a4_ad(2, i, k) = 0.0
              END IF
            END DO
          ELSE
            DO i=i2,i1,-1
              CALL POPREALARRAY_ADM(a4(4, i, k))
              temp_ad23 = 3.*a4_ad(4, i, k)
              a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad23
              a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad23
              a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad23
              a4_ad(4, i, k) = 0.0
            END DO
          END IF
        ELSE IF (branch .EQ. 6) THEN
          DO i=i2,i1,-1
            CALL POPCONTROL2B(branch)
            IF (branch .NE. 0) THEN
              IF (branch .NE. 1) THEN
                CALL POPREALARRAY_ADM(a4(4, i, k))
                temp_ad24 = 3.*a4_ad(4, i, k)
                a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad24
                a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad24
                a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad24
                a4_ad(4, i, k) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  y20_ad = a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  x12_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  x12_ad = a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  y20_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = y20_ad
                    pmp_2_ad = 0.0
                  ELSE
                    pmp_2_ad = y20_ad
                    lac_2_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_2_ad = y20_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y20_ad
                    lac_2_ad = 0.0
                  END IF
                  pmp_2_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y32_ad = x12_ad
                ELSE
                  a4_ad(3, i, k) = a4_ad(3, i, k) + x12_ad
                  y32_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = lac_2_ad + y32_ad
                  ELSE
                    pmp_2_ad = pmp_2_ad + y32_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_2_ad = lac_2_ad + y32_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y32_ad
                END IF
                pmp_2_ad = pmp_2_ad + lac_2_ad
                gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  y19_ad = a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  x11_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  x11_ad = a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  y19_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = y19_ad
                    pmp_1_ad = 0.0
                  ELSE
                    pmp_1_ad = y19_ad
                    lac_1_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_1_ad = y19_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y19_ad
                    lac_1_ad = 0.0
                  END IF
                  pmp_1_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y31_ad = x11_ad
                ELSE
                  a4_ad(2, i, k) = a4_ad(2, i, k) + x11_ad
                  y31_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = lac_1_ad + y31_ad
                  ELSE
                    pmp_1_ad = pmp_1_ad + y31_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_1_ad = lac_1_ad + y31_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y31_ad
                END IF
                pmp_1_ad = pmp_1_ad + lac_1_ad
                gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
              END IF
            END IF
          END DO
        ELSE
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(a4(4, i, k))
              temp_ad25 = 3.*a4_ad(4, i, k)
              a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad25
              a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad25
              a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad25
              a4_ad(4, i, k) = 0.0
            ELSE
              CALL POPREALARRAY_ADM(a4(4, i, k))
              a4_ad(4, i, k) = 0.0
              CALL POPREALARRAY_ADM(a4(3, i, k))
              a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
              a4_ad(3, i, k) = 0.0
              CALL POPREALARRAY_ADM(a4(2, i, k))
              a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
              a4_ad(2, i, k) = 0.0
            END IF
          END DO
        END IF
      END DO
      CALL CS_LIMITERS_BWD(im, extm(i1, 2), a4(1, i1, 2), a4_ad(1, i1, 2&
&                    ), 2)
      DO i=i2,i1,-1
        CALL POPREALARRAY_ADM(a4(4, i, 2))
        temp_ad17 = 3.*a4_ad(4, i, 2)
        a4_ad(1, i, 2) = a4_ad(1, i, 2) + 2.*temp_ad17
        a4_ad(2, i, 2) = a4_ad(2, i, 2) - temp_ad17
        a4_ad(3, i, 2) = a4_ad(3, i, 2) - temp_ad17
        a4_ad(4, i, 2) = 0.0
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        CALL CS_LIMITERS_BWD(im, extm(i1, 1), a4(1, i1, 1), a4_ad(1, i1&
&                      , 1), 1)
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(a4(4, i, 1))
          temp_ad16 = 3.*a4_ad(4, i, 1)
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + 2.*temp_ad16
          a4_ad(2, i, 1) = a4_ad(2, i, 1) - temp_ad16
          a4_ad(3, i, 1) = a4_ad(3, i, 1) - temp_ad16
          a4_ad(4, i, 1) = 0.0
        END DO
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(a4(2, i, 1))
            ELSE
              CALL POPREALARRAY_ADM(a4(2, i, 1))
              a4_ad(2, i, 1) = 0.0
            END IF
          END DO
        ELSE
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) THEN
              CALL POPREALARRAY_ADM(a4(2, i, 1))
              a4_ad(2, i, 1) = 0.0
            END IF
          END DO
        END IF
      ELSE IF (branch .EQ. 2) THEN
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(a4(4, i, 1))
          a4_ad(4, i, 1) = 0.0
          CALL POPREALARRAY_ADM(a4(3, i, 1))
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + a4_ad(3, i, 1)
          a4_ad(3, i, 1) = 0.0
          CALL POPREALARRAY_ADM(a4(2, i, 1))
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + a4_ad(2, i, 1)
          a4_ad(2, i, 1) = 0.0
        END DO
      END IF
      DO k=km,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          DO i=i2,i1,-1
            CALL POPREALARRAY_ADM(a4(4, i, k))
            temp_ad15 = 3.*a4_ad(4, i, k)
            a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad15
            a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad15
            a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad15
            a4_ad(4, i, k) = 0.0
          END DO
        END IF
      END DO
      q_ad = 0.0
      DO k=km,1,-1
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(a4(3, i, k))
          q_ad(i, k+1) = q_ad(i, k+1) + a4_ad(3, i, k)
          a4_ad(3, i, k) = 0.0
          CALL POPREALARRAY_ADM(a4(2, i, k))
          q_ad(i, k) = q_ad(i, k) + a4_ad(2, i, k)
          a4_ad(2, i, k) = 0.0
        END DO
      END DO
      DO i=i2,i1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y8_ad = q_ad(i, km)
          q_ad(i, km) = 0.0
        ELSE
          y8_ad = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          a4_ad(1, i, km) = a4_ad(1, i, km) + y8_ad
        ELSE
          a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + y8_ad
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY_ADM(q(i, km))
          y7_ad = q_ad(i, km)
          q_ad(i, km) = 0.0
        ELSE
          CALL POPREALARRAY_ADM(q(i, km))
          y7_ad = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          a4_ad(1, i, km) = a4_ad(1, i, km) + y7_ad
        ELSE
          a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + y7_ad
        END IF
      END DO
      DO k=km-1,3,-1
        DO 120 i=i2,i1,-1
          CALL POPCONTROL3B(branch)
          IF (branch .NE. 0) THEN
            IF (branch .LT. 4) THEN
              IF (branch .EQ. 1) THEN
                GOTO 110
              ELSE IF (branch .EQ. 2) THEN
                q_ad(i, k) = 0.0
                GOTO 110
              ELSE
                CALL POPREALARRAY_ADM(q(i, k))
                y5_ad = q_ad(i, k)
                q_ad(i, k) = 0.0
              END IF
            ELSE IF (branch .EQ. 4) THEN
              CALL POPREALARRAY_ADM(q(i, k))
              y5_ad = 0.0
            ELSE
              IF (branch .EQ. 5) THEN
                y4_ad = q_ad(i, k)
                q_ad(i, k) = 0.0
              ELSE
                y4_ad = 0.0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                a4_ad(1, i, k) = a4_ad(1, i, k) + y4_ad
              ELSE
                a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + y4_ad
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(q(i, k))
                y3_ad = q_ad(i, k)
                q_ad(i, k) = 0.0
              ELSE
                CALL POPREALARRAY_ADM(q(i, k))
                y3_ad = 0.0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                a4_ad(1, i, k) = a4_ad(1, i, k) + y3_ad
              ELSE
                a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + y3_ad
              END IF
              GOTO 120
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              a4_ad(1, i, k) = a4_ad(1, i, k) + y5_ad
            ELSE
              a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + y5_ad
            END IF
            GOTO 120
          END IF
 110      CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY_ADM(q(i, k))
            y6_ad = q_ad(i, k)
            q_ad(i, k) = 0.0
          ELSE
            CALL POPREALARRAY_ADM(q(i, k))
            y6_ad = 0.0
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            a4_ad(1, i, k) = a4_ad(1, i, k) + y6_ad
          ELSE
            a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + y6_ad
          END IF
 120    CONTINUE
      END DO
      DO k=km,2,-1
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(gam(i, k))
          a4_ad(1, i, k) = a4_ad(1, i, k) + gam_ad(i, k)
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) - gam_ad(i, k)
          gam_ad(i, k) = 0.0
        END DO
      END DO
      DO i=i2,i1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y2_ad = q_ad(i, 2)
          q_ad(i, 2) = 0.0
        ELSE
          y2_ad = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          a4_ad(1, i, 2) = a4_ad(1, i, 2) + y2_ad
        ELSE
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + y2_ad
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY_ADM(q(i, 2))
          y1_ad = q_ad(i, 2)
          q_ad(i, 2) = 0.0
        ELSE
          CALL POPREALARRAY_ADM(q(i, 2))
          y1_ad = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          a4_ad(1, i, 2) = a4_ad(1, i, 2) + y1_ad
        ELSE
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + y1_ad
        END IF
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO k=1,km,1
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(q(i, k))
          gam_ad(i, k) = gam_ad(i, k) - q(i, k+1)*q_ad(i, k)
          q_ad(i, k+1) = q_ad(i, k+1) - gam(i, k)*q_ad(i, k)
        END DO
      END DO
      d4_ad = 0.0
      DO i=i2,i1,-1
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        CALL POPREALARRAY_ADM(q(i, km+1))
        temp2 = d4(i)*(d4(i)+0.5) - a_bot*gam(i, km)
        temp_ad11 = q_ad(i, km+1)/temp2
        temp1 = d4(i)*(d4(i)+1.)
        temp_ad12 = 2.*a4(1, i, km)*temp_ad11
        temp_ad13 = -((2.*(temp1*a4(1, i, km))+a4(1, i, km-1)-a_bot*q(i&
&         , km))*temp_ad11/temp2)
        a4_ad(1, i, km) = a4_ad(1, i, km) + 2.*temp1*temp_ad11
        a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + temp_ad11
        a_bot_ad = -(gam(i, km)*temp_ad13) - q(i, km)*temp_ad11
        d4_ad(i) = d4_ad(i) + (2*d4(i)+1.5)*a_bot_ad + (2*d4(i)+0.5)*&
&         temp_ad13 + (2*d4(i)+1.)*temp_ad12
        q_ad(i, km) = q_ad(i, km) - a_bot*temp_ad11
        gam_ad(i, km) = gam_ad(i, km) - a_bot*temp_ad13
        q_ad(i, km+1) = 0.0
      END DO
      DO k=km,2,-1
        DO i=i2,i1,-1
          temp_ad9 = q_ad(i, k)/bet
          temp_ad8 = 3.*temp_ad9
          CALL POPREALARRAY_ADM(q(i, k))
          bet_ad = -((3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))*&
&           temp_ad9/bet) - d4(i)*gam_ad(i, k)/bet**2
          d4_ad(i) = d4_ad(i) + a4(1, i, k)*temp_ad8 + 2*bet_ad + gam_ad&
&           (i, k)/bet
          gam_ad(i, k) = 0.0
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + temp_ad8
          a4_ad(1, i, k) = a4_ad(1, i, k) + d4(i)*temp_ad8
          q_ad(i, k-1) = q_ad(i, k-1) - temp_ad9
          q_ad(i, k) = 0.0
          CALL POPREALARRAY_ADM(bet)
          gam_ad(i, k-1) = gam_ad(i, k-1) - bet_ad
          CALL POPREALARRAY_ADM(d4(i))
          temp_ad10 = d4_ad(i)/delp(i, k)
          delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad10
          delp_ad(i, k) = delp_ad(i, k) - delp(i, k-1)*temp_ad10/delp(i&
&           , k)
          d4_ad(i) = 0.0
        END DO
      END DO
      DO i=i2,i1,-1
        grat = delp(i, 2)/delp(i, 1)
        bet = grat*(grat+0.5)
        temp_ad4 = gam_ad(i, 1)/bet
        gam_ad(i, 1) = 0.0
        temp_ad6 = q_ad(i, 1)/bet
        temp_ad5 = a4(1, i, 1)*temp_ad6
        temp0 = 2*grat*(grat+1.)
        bet_ad = -((temp0*a4(1, i, 1)+a4(1, i, 2))*temp_ad6/bet) - (grat&
&         *(grat+1.5)+1.)*temp_ad4/bet
        grat_ad = (4*grat+2*1.)*temp_ad5 + (2*grat+0.5)*bet_ad + (2*grat&
&         +1.5)*temp_ad4
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + temp0*temp_ad6
        a4_ad(1, i, 2) = a4_ad(1, i, 2) + temp_ad6
        q_ad(i, 1) = 0.0
        temp_ad7 = grat_ad/delp(i, 1)
        delp_ad(i, 2) = delp_ad(i, 2) + temp_ad7
        delp_ad(i, 1) = delp_ad(i, 1) - delp(i, 2)*temp_ad7/delp(i, 1)
      END DO
    ELSE
      DO k=1,km-1,1
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(q(i, k))
          gam_ad(i, k+1) = gam_ad(i, k+1) - q(i, k+1)*q_ad(i, k)
          q_ad(i, k+1) = q_ad(i, k+1) - gam(i, k+1)*q_ad(i, k)
        END DO
      END DO
      DO i=i2,i1,-1
        CALL POPREALARRAY_ADM(q(i, km+1))
        q_ad(i, km+1) = 0.0
        grat = delp(i, km-1)/delp(i, km)
        CALL POPREALARRAY_ADM(q(i, km))
        temp = 2*grat - gam(i, km) + 2.
        temp_ad1 = q_ad(i, km)/temp
        temp_ad2 = -((3.*(a4(1, i, km-1)+a4(1, i, km))-qs(i)*grat-q(i, &
&         km-1))*temp_ad1/temp)
        a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + 3.*temp_ad1
        a4_ad(1, i, km) = a4_ad(1, i, km) + 3.*temp_ad1
        grat_ad = 2*temp_ad2 - qs(i)*temp_ad1
        q_ad(i, km-1) = q_ad(i, km-1) - temp_ad1
        gam_ad(i, km) = gam_ad(i, km) - temp_ad2
        q_ad(i, km) = 0.0
        temp_ad3 = grat_ad/delp(i, km)
        delp_ad(i, km-1) = delp_ad(i, km-1) + temp_ad3
        delp_ad(i, km) = delp_ad(i, km) - delp(i, km-1)*temp_ad3/delp(i&
&         , km)
      END DO
      DO k=km-1,2,-1
        DO i=i2,i1,-1
          temp_ad = q_ad(i, k)/bet
          CALL POPREALARRAY_ADM(q(i, k))
          grat = delp(i, k-1)/delp(i, k)
          bet_ad = -((3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))*temp_ad/&
&           bet) - grat*gam_ad(i, k+1)/bet**2
          grat_ad = 2*bet_ad + gam_ad(i, k+1)/bet
          gam_ad(i, k+1) = 0.0
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + 3.*temp_ad
          a4_ad(1, i, k) = a4_ad(1, i, k) + 3.*temp_ad
          q_ad(i, k-1) = q_ad(i, k-1) - temp_ad
          q_ad(i, k) = 0.0
          CALL POPREALARRAY_ADM(bet)
          gam_ad(i, k) = gam_ad(i, k) - bet_ad
          temp_ad0 = grat_ad/delp(i, k)
          delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad0
          delp_ad(i, k) = delp_ad(i, k) - delp(i, k-1)*temp_ad0/delp(i, &
&           k)
        END DO
      END DO
      DO i=i2,i1,-1
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + 1.5*q_ad(i, 1)
        q_ad(i, 1) = 0.0
      END DO
    END IF
  END SUBROUTINE SCALAR_PROFILE_ADM
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
!  Differentiation of cs_profile in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn
!_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dy
!n_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_
!mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynam
!ics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils
!_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_m
!od.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scal
!ar_profile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.ste
!epz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_
!setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.co
!mpute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver
!_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver
! nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh s
!w_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_cor
!e_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.
!compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corner
!s_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circ
!le_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: qs delp a4
!   with respect to varying inputs: qs delp a4
  SUBROUTINE CS_PROFILE_ADM(qs, qs_ad, a4, a4_ad, delp, delp_ad, km, i1&
&   , i2, iv, kord)
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
    REAL :: qs_ad(i1:i2)
! layer pressure thickness
    REAL, INTENT(IN) :: delp(i1:i2, km)
    REAL :: delp_ad(i1:i2, km)
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_ad(4, i1:i2, km)
!-----------------------------------------------------------------------
    LOGICAL :: extm(i1:i2, km)
    REAL :: gam(i1:i2, km)
    REAL :: gam_ad(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: q_ad(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: d4_ad(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: bet_ad, a_bot_ad, grat_ad
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: pmp_1_ad, lac_1_ad, pmp_2_ad, lac_2_ad
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
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp0
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp1
    REAL :: temp2
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: y1_ad
    REAL :: y2_ad
    REAL :: y3_ad
    REAL :: y4_ad
    REAL :: y5_ad
    REAL :: y6_ad
    REAL :: y7_ad
    REAL :: y8_ad
    REAL :: temp_ad15
    REAL :: temp_ad16
    REAL :: y19_ad
    REAL :: x1_ad
    REAL :: y9_ad
    REAL :: y20_ad
    REAL :: x2_ad
    REAL :: y10_ad
    REAL :: temp_ad17
    REAL :: y21_ad
    REAL :: x3_ad
    REAL :: y11_ad
    REAL :: y22_ad
    REAL :: x4_ad
    REAL :: y12_ad
    REAL :: y23_ad
    REAL :: x5_ad
    REAL :: y13_ad
    REAL :: y24_ad
    REAL :: x6_ad
    REAL :: y14_ad
    REAL :: y25_ad
    REAL :: x7_ad
    REAL :: y15_ad
    REAL :: y26_ad
    REAL :: x8_ad
    REAL :: y16_ad
    REAL :: y27_ad
    REAL :: x9_ad
    REAL :: y17_ad
    REAL :: y28_ad
    REAL :: x10_ad
    REAL :: y18_ad
    REAL :: temp_ad18
    REAL :: temp_ad19
    REAL :: temp_ad20
    REAL :: temp_ad21
    REAL :: temp_ad22
    INTEGER :: branch
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
          CALL PUSHREALARRAY_ADM(bet)
          bet = 2. + grat + grat - gam(i, k)
          CALL PUSHREALARRAY_ADM(q(i, k))
          q(i, k) = (3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))/bet
          gam(i, k+1) = grat/bet
        END DO
      END DO
      DO i=i1,i2
        grat = delp(i, km-1)/delp(i, km)
        CALL PUSHREALARRAY_ADM(q(i, km))
        q(i, km) = (3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-&
&         1))/(2.+grat+grat-gam(i, km))
        CALL PUSHREALARRAY_ADM(q(i, km+1))
        q(i, km+1) = qs(i)
      END DO
      DO k=km-1,1,-1
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(q(i, k))
          q(i, k) = q(i, k) - gam(i, k+1)*q(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL1B(1)
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
          CALL PUSHREALARRAY_ADM(d4(i))
          d4(i) = delp(i, k-1)/delp(i, k)
          CALL PUSHREALARRAY_ADM(bet)
          bet = 2. + d4(i) + d4(i) - gam(i, k-1)
          CALL PUSHREALARRAY_ADM(q(i, k))
          q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
          gam(i, k) = d4(i)/bet
        END DO
      END DO
      DO i=i1,i2
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        CALL PUSHREALARRAY_ADM(q(i, km+1))
        q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(q(i, k))
          q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL1B(0)
    END IF
    IF (kord .GE. 0.) THEN
      abs0 = kord
    ELSE
      abs0 = -kord
    END IF
!----- Perfectly linear scheme --------------------------------
    IF (abs0 .GT. 16) THEN
      q_ad = 0.0
      DO k=km,1,-1
        DO i=i2,i1,-1
          temp_ad14 = 3.*a4_ad(4, i, k)
          a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad14
          a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad14
          a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad14
          a4_ad(4, i, k) = 0.0
          q_ad(i, k+1) = q_ad(i, k+1) + a4_ad(3, i, k)
          a4_ad(3, i, k) = 0.0
          q_ad(i, k) = q_ad(i, k) + a4_ad(2, i, k)
          a4_ad(2, i, k) = 0.0
        END DO
      END DO
      gam_ad = 0.0
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
          CALL PUSHCONTROL1B(0)
        ELSE
          y1 = a4(1, i, 1)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (q(i, 2) .GT. y1) THEN
          CALL PUSHREALARRAY_ADM(q(i, 2))
          q(i, 2) = y1
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREALARRAY_ADM(q(i, 2))
          q(i, 2) = q(i, 2)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (a4(1, i, 1) .GT. a4(1, i, 2)) THEN
          y2 = a4(1, i, 2)
          CALL PUSHCONTROL1B(0)
        ELSE
          y2 = a4(1, i, 1)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (q(i, 2) .LT. y2) THEN
          q(i, 2) = y2
          CALL PUSHCONTROL1B(0)
        ELSE
          q(i, 2) = q(i, 2)
          CALL PUSHCONTROL1B(1)
        END IF
      END DO
      DO k=2,km
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(gam(i, k))
          gam(i, k) = a4(1, i, k) - a4(1, i, k-1)
        END DO
      END DO
! Interior:
      DO k=3,km-1
        DO i=i1,i2
          IF (gam(i, k-1)*gam(i, k+1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y3 = a4(1, i, k)
              CALL PUSHCONTROL1B(0)
            ELSE
              y3 = a4(1, i, k-1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, k) .GT. y3) THEN
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = y3
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = q(i, k)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y4 = a4(1, i, k)
              CALL PUSHCONTROL1B(0)
            ELSE
              y4 = a4(1, i, k-1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, k) .LT. y4) THEN
              q(i, k) = y4
              CALL PUSHCONTROL3B(5)
            ELSE
              q(i, k) = q(i, k)
              CALL PUSHCONTROL3B(6)
            END IF
          ELSE IF (gam(i, k-1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y5 = a4(1, i, k)
              CALL PUSHCONTROL1B(0)
            ELSE
              y5 = a4(1, i, k-1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, k) .LT. y5) THEN
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = y5
              CALL PUSHCONTROL3B(3)
            ELSE
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = q(i, k)
              CALL PUSHCONTROL3B(4)
            END IF
          ELSE
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y6 = a4(1, i, k)
              CALL PUSHCONTROL1B(0)
            ELSE
              y6 = a4(1, i, k-1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, k) .GT. y6) THEN
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = y6
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(q(i, k))
              q(i, k) = q(i, k)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (iv .EQ. 0) THEN
              IF (0. .LT. q(i, k)) THEN
                CALL PUSHCONTROL3B(0)
                q(i, k) = q(i, k)
              ELSE
                q(i, k) = 0.
                CALL PUSHCONTROL3B(2)
              END IF
            ELSE
              CALL PUSHCONTROL3B(1)
            END IF
          END IF
        END DO
      END DO
! Bottom:
      DO i=i1,i2
        IF (a4(1, i, km-1) .LT. a4(1, i, km)) THEN
          y7 = a4(1, i, km)
          CALL PUSHCONTROL1B(0)
        ELSE
          y7 = a4(1, i, km-1)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (q(i, km) .GT. y7) THEN
          CALL PUSHREALARRAY_ADM(q(i, km))
          q(i, km) = y7
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREALARRAY_ADM(q(i, km))
          q(i, km) = q(i, km)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (a4(1, i, km-1) .GT. a4(1, i, km)) THEN
          y8 = a4(1, i, km)
          CALL PUSHCONTROL1B(0)
        ELSE
          y8 = a4(1, i, km-1)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (q(i, km) .LT. y8) THEN
          q(i, km) = y8
          CALL PUSHCONTROL1B(0)
        ELSE
          q(i, km) = q(i, km)
          CALL PUSHCONTROL1B(1)
        END IF
      END DO
      DO k=1,km
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(a4(2, i, k))
          a4(2, i, k) = q(i, k)
          CALL PUSHREALARRAY_ADM(a4(3, i, k))
          a4(3, i, k) = q(i, k+1)
        END DO
      END DO
      DO k=1,km
        IF (k .EQ. 1 .OR. k .EQ. km) THEN
          CALL PUSHCONTROL1B(1)
          DO i=i1,i2
            extm(i, k) = (a4(2, i, k)-a4(1, i, k))*(a4(3, i, k)-a4(1, i&
&             , k)) .GT. 0.
          END DO
        ELSE
          CALL PUSHCONTROL1B(0)
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
            CALL PUSHREALARRAY_ADM(a4(2, i, 1))
            a4(2, i, 1) = a4(2, i, 1)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREALARRAY_ADM(a4(2, i, 1))
            a4(2, i, 1) = 0.
            CALL PUSHCONTROL1B(1)
          END IF
        END DO
        CALL PUSHCONTROL2B(0)
      ELSE IF (iv .EQ. -1) THEN
        DO i=i1,i2
          IF (a4(2, i, 1)*a4(1, i, 1) .LE. 0.) THEN
            CALL PUSHREALARRAY_ADM(a4(2, i, 1))
            a4(2, i, 1) = 0.
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL2B(1)
      ELSE IF (iv .EQ. 2) THEN
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(a4(2, i, 1))
          a4(2, i, 1) = a4(1, i, 1)
          CALL PUSHREALARRAY_ADM(a4(3, i, 1))
          a4(3, i, 1) = a4(1, i, 1)
          CALL PUSHREALARRAY_ADM(a4(4, i, 1))
          a4(4, i, 1) = 0.
        END DO
        CALL PUSHCONTROL2B(2)
      ELSE
        CALL PUSHCONTROL2B(3)
      END IF
      IF (iv .NE. 2) THEN
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(a4(4, i, 1))
          a4(4, i, 1) = 3.*(2.*a4(1, i, 1)-(a4(2, i, 1)+a4(3, i, 1)))
        END DO
        CALL CS_LIMITERS_FWD(im, extm(i1, 1), a4(1, i1, 1), 1)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
! k=2
      DO i=i1,i2
        CALL PUSHREALARRAY_ADM(a4(4, i, 2))
        a4(4, i, 2) = 3.*(2.*a4(1, i, 2)-(a4(2, i, 2)+a4(3, i, 2)))
      END DO
      CALL CS_LIMITERS_FWD(im, extm(i1, 2), a4(1, i1, 2), 2)
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
                CALL PUSHCONTROL2B(0)
              ELSE
                y19 = pmp_1
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (a4(1, i, k) .GT. lac_1) THEN
              y19 = lac_1
              CALL PUSHCONTROL2B(2)
            ELSE
              y19 = a4(1, i, k)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (a4(2, i, k) .LT. y19) THEN
              x1 = y19
              CALL PUSHCONTROL1B(0)
            ELSE
              x1 = a4(2, i, k)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (a4(1, i, k) .LT. pmp_1) THEN
              IF (pmp_1 .LT. lac_1) THEN
                y9 = lac_1
                CALL PUSHCONTROL2B(0)
              ELSE
                y9 = pmp_1
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (a4(1, i, k) .LT. lac_1) THEN
              y9 = lac_1
              CALL PUSHCONTROL2B(2)
            ELSE
              y9 = a4(1, i, k)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (x1 .GT. y9) THEN
              CALL PUSHREALARRAY_ADM(a4(2, i, k))
              a4(2, i, k) = y9
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(a4(2, i, k))
              a4(2, i, k) = x1
              CALL PUSHCONTROL1B(1)
            END IF
! Right edges
            pmp_2 = a4(1, i, k) + 2.*gam(i, k)
            lac_2 = pmp_2 - 1.5*gam(i, k-1)
            IF (a4(1, i, k) .GT. pmp_2) THEN
              IF (pmp_2 .GT. lac_2) THEN
                y20 = lac_2
                CALL PUSHCONTROL2B(0)
              ELSE
                y20 = pmp_2
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (a4(1, i, k) .GT. lac_2) THEN
              y20 = lac_2
              CALL PUSHCONTROL2B(2)
            ELSE
              y20 = a4(1, i, k)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (a4(3, i, k) .LT. y20) THEN
              x2 = y20
              CALL PUSHCONTROL1B(0)
            ELSE
              x2 = a4(3, i, k)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (a4(1, i, k) .LT. pmp_2) THEN
              IF (pmp_2 .LT. lac_2) THEN
                y10 = lac_2
                CALL PUSHCONTROL2B(0)
              ELSE
                y10 = pmp_2
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (a4(1, i, k) .LT. lac_2) THEN
              y10 = lac_2
              CALL PUSHCONTROL2B(2)
            ELSE
              y10 = a4(1, i, k)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (x2 .GT. y10) THEN
              CALL PUSHREALARRAY_ADM(a4(3, i, k))
              a4(3, i, k) = y10
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(a4(3, i, k))
              a4(3, i, k) = x2
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREALARRAY_ADM(a4(4, i, k))
            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
          END DO
          CALL PUSHCONTROL3B(0)
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
                CALL PUSHREALARRAY_ADM(a4(2, i, k))
                a4(2, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(3, i, k))
                a4(3, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(4, i, k))
                a4(4, i, k) = 0.
                CALL PUSHCONTROL2B(3)
              ELSE IF (extm(i, k) .AND. extm(i, k+1)) THEN
! c90_mp122
! grid-scale 2-delta-z wave detected
                CALL PUSHREALARRAY_ADM(a4(2, i, k))
                a4(2, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(3, i, k))
                a4(3, i, k) = a4(1, i, k)
                CALL PUSHREALARRAY_ADM(a4(4, i, k))
                a4(4, i, k) = 0.
                CALL PUSHCONTROL2B(2)
              ELSE
                CALL PUSHREALARRAY_ADM(a4(4, i, k))
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
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y21 = pmp_1
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                    y21 = lac_1
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y21 = a4(1, i, k)
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (a4(2, i, k) .LT. y21) THEN
                    x3 = y21
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    x3 = a4(2, i, k)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_1) THEN
                    IF (pmp_1 .LT. lac_1) THEN
                      y11 = lac_1
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y11 = pmp_1
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                    y11 = lac_1
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y11 = a4(1, i, k)
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (x3 .GT. y11) THEN
                    CALL PUSHREALARRAY_ADM(a4(2, i, k))
                    a4(2, i, k) = y11
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREALARRAY_ADM(a4(2, i, k))
                    a4(2, i, k) = x3
                    CALL PUSHCONTROL1B(1)
                  END IF
                  pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                  lac_2 = pmp_2 - 1.5*gam(i, k-1)
                  IF (a4(1, i, k) .GT. pmp_2) THEN
                    IF (pmp_2 .GT. lac_2) THEN
                      y22 = lac_2
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y22 = pmp_2
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                    y22 = lac_2
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y22 = a4(1, i, k)
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (a4(3, i, k) .LT. y22) THEN
                    x4 = y22
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    x4 = a4(3, i, k)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_2) THEN
                    IF (pmp_2 .LT. lac_2) THEN
                      y12 = lac_2
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y12 = pmp_2
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                    y12 = lac_2
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y12 = a4(1, i, k)
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (x4 .GT. y12) THEN
                    CALL PUSHREALARRAY_ADM(a4(3, i, k))
                    a4(3, i, k) = y12
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREALARRAY_ADM(a4(3, i, k))
                    a4(3, i, k) = x4
                    CALL PUSHCONTROL1B(1)
                  END IF
                  CALL PUSHREALARRAY_ADM(a4(4, i, k))
                  a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i&
&                   , k))
                  CALL PUSHCONTROL2B(1)
                ELSE
                  CALL PUSHCONTROL2B(0)
                END IF
              END IF
            END DO
            CALL PUSHCONTROL3B(1)
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
                    CALL PUSHREALARRAY_ADM(a4(2, i, k))
                    a4(2, i, k) = a4(1, i, k)
                    CALL PUSHREALARRAY_ADM(a4(3, i, k))
                    a4(3, i, k) = a4(1, i, k)
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
                    a4(4, i, k) = 0.
                    CALL PUSHCONTROL2B(3)
                  ELSE
! True local extremum
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                    CALL PUSHCONTROL2B(2)
                  END IF
                ELSE
! not a local extremum
                  CALL PUSHREALARRAY_ADM(a4(4, i, k))
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
                        CALL PUSHCONTROL2B(0)
                      ELSE
                        y23 = pmp_1
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                      y23 = lac_1
                      CALL PUSHCONTROL2B(2)
                    ELSE
                      y23 = a4(1, i, k)
                      CALL PUSHCONTROL2B(3)
                    END IF
                    IF (a4(2, i, k) .LT. y23) THEN
                      x5 = y23
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      x5 = a4(2, i, k)
                      CALL PUSHCONTROL1B(1)
                    END IF
                    IF (a4(1, i, k) .LT. pmp_1) THEN
                      IF (pmp_1 .LT. lac_1) THEN
                        y13 = lac_1
                        CALL PUSHCONTROL2B(0)
                      ELSE
                        y13 = pmp_1
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                      y13 = lac_1
                      CALL PUSHCONTROL2B(2)
                    ELSE
                      y13 = a4(1, i, k)
                      CALL PUSHCONTROL2B(3)
                    END IF
                    IF (x5 .GT. y13) THEN
                      CALL PUSHREALARRAY_ADM(a4(2, i, k))
                      a4(2, i, k) = y13
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      CALL PUSHREALARRAY_ADM(a4(2, i, k))
                      a4(2, i, k) = x5
                      CALL PUSHCONTROL1B(1)
                    END IF
                    pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                    lac_2 = pmp_2 - 1.5*gam(i, k-1)
                    IF (a4(1, i, k) .GT. pmp_2) THEN
                      IF (pmp_2 .GT. lac_2) THEN
                        y24 = lac_2
                        CALL PUSHCONTROL2B(0)
                      ELSE
                        y24 = pmp_2
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                      y24 = lac_2
                      CALL PUSHCONTROL2B(2)
                    ELSE
                      y24 = a4(1, i, k)
                      CALL PUSHCONTROL2B(3)
                    END IF
                    IF (a4(3, i, k) .LT. y24) THEN
                      x6 = y24
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      x6 = a4(3, i, k)
                      CALL PUSHCONTROL1B(1)
                    END IF
                    IF (a4(1, i, k) .LT. pmp_2) THEN
                      IF (pmp_2 .LT. lac_2) THEN
                        y14 = lac_2
                        CALL PUSHCONTROL2B(0)
                      ELSE
                        y14 = pmp_2
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                      y14 = lac_2
                      CALL PUSHCONTROL2B(2)
                    ELSE
                      y14 = a4(1, i, k)
                      CALL PUSHCONTROL2B(3)
                    END IF
                    IF (x6 .GT. y14) THEN
                      CALL PUSHREALARRAY_ADM(a4(3, i, k))
                      a4(3, i, k) = y14
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      CALL PUSHREALARRAY_ADM(a4(3, i, k))
                      a4(3, i, k) = x6
                      CALL PUSHCONTROL1B(1)
                    END IF
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
                    a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3&
&                     , i, k))
                    CALL PUSHCONTROL2B(1)
                  ELSE
                    CALL PUSHCONTROL2B(0)
                  END IF
                END IF
              END DO
              CALL PUSHCONTROL3B(2)
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
                    CALL PUSHREALARRAY_ADM(a4(2, i, k))
                    a4(2, i, k) = a4(1, i, k)
                    CALL PUSHREALARRAY_ADM(a4(3, i, k))
                    a4(3, i, k) = a4(1, i, k)
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
                    a4(4, i, k) = 0.
                    CALL PUSHCONTROL2B(2)
                  ELSE
! not a local extremum
                    CALL PUSHREALARRAY_ADM(a4(4, i, k))
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
                          CALL PUSHCONTROL2B(0)
                        ELSE
                          y25 = pmp_1
                          CALL PUSHCONTROL2B(1)
                        END IF
                      ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                        y25 = lac_1
                        CALL PUSHCONTROL2B(2)
                      ELSE
                        y25 = a4(1, i, k)
                        CALL PUSHCONTROL2B(3)
                      END IF
                      IF (a4(2, i, k) .LT. y25) THEN
                        x7 = y25
                        CALL PUSHCONTROL1B(0)
                      ELSE
                        x7 = a4(2, i, k)
                        CALL PUSHCONTROL1B(1)
                      END IF
                      IF (a4(1, i, k) .LT. pmp_1) THEN
                        IF (pmp_1 .LT. lac_1) THEN
                          y15 = lac_1
                          CALL PUSHCONTROL2B(0)
                        ELSE
                          y15 = pmp_1
                          CALL PUSHCONTROL2B(1)
                        END IF
                      ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                        y15 = lac_1
                        CALL PUSHCONTROL2B(2)
                      ELSE
                        y15 = a4(1, i, k)
                        CALL PUSHCONTROL2B(3)
                      END IF
                      IF (x7 .GT. y15) THEN
                        CALL PUSHREALARRAY_ADM(a4(2, i, k))
                        a4(2, i, k) = y15
                        CALL PUSHCONTROL1B(0)
                      ELSE
                        CALL PUSHREALARRAY_ADM(a4(2, i, k))
                        a4(2, i, k) = x7
                        CALL PUSHCONTROL1B(1)
                      END IF
                      pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                      lac_2 = pmp_2 - 1.5*gam(i, k-1)
                      IF (a4(1, i, k) .GT. pmp_2) THEN
                        IF (pmp_2 .GT. lac_2) THEN
                          y26 = lac_2
                          CALL PUSHCONTROL2B(0)
                        ELSE
                          y26 = pmp_2
                          CALL PUSHCONTROL2B(1)
                        END IF
                      ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                        y26 = lac_2
                        CALL PUSHCONTROL2B(2)
                      ELSE
                        y26 = a4(1, i, k)
                        CALL PUSHCONTROL2B(3)
                      END IF
                      IF (a4(3, i, k) .LT. y26) THEN
                        x8 = y26
                        CALL PUSHCONTROL1B(0)
                      ELSE
                        x8 = a4(3, i, k)
                        CALL PUSHCONTROL1B(1)
                      END IF
                      IF (a4(1, i, k) .LT. pmp_2) THEN
                        IF (pmp_2 .LT. lac_2) THEN
                          y16 = lac_2
                          CALL PUSHCONTROL2B(0)
                        ELSE
                          y16 = pmp_2
                          CALL PUSHCONTROL2B(1)
                        END IF
                      ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                        y16 = lac_2
                        CALL PUSHCONTROL2B(2)
                      ELSE
                        y16 = a4(1, i, k)
                        CALL PUSHCONTROL2B(3)
                      END IF
                      IF (x8 .GT. y16) THEN
                        CALL PUSHREALARRAY_ADM(a4(3, i, k))
                        a4(3, i, k) = y16
                        CALL PUSHCONTROL1B(0)
                      ELSE
                        CALL PUSHREALARRAY_ADM(a4(3, i, k))
                        a4(3, i, k) = x8
                        CALL PUSHCONTROL1B(1)
                      END IF
                      CALL PUSHREALARRAY_ADM(a4(4, i, k))
                      a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(&
&                       3, i, k))
                      CALL PUSHCONTROL2B(1)
                    ELSE
                      CALL PUSHCONTROL2B(0)
                    END IF
                  END IF
                END DO
                CALL PUSHCONTROL3B(3)
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
                        CALL PUSHREALARRAY_ADM(a4(2, i, k))
                        a4(2, i, k) = a4(1, i, k)
                        CALL PUSHREALARRAY_ADM(a4(3, i, k))
                        a4(3, i, k) = a4(1, i, k)
                        CALL PUSHREALARRAY_ADM(a4(4, i, k))
                        a4(4, i, k) = 0.
                        CALL PUSHCONTROL2B(2)
                      ELSE
! Left  edges
                        pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                        lac_1 = pmp_1 + 1.5*gam(i, k+2)
                        IF (a4(1, i, k) .GT. pmp_1) THEN
                          IF (pmp_1 .GT. lac_1) THEN
                            y27 = lac_1
                            CALL PUSHCONTROL2B(0)
                          ELSE
                            y27 = pmp_1
                            CALL PUSHCONTROL2B(1)
                          END IF
                        ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                          y27 = lac_1
                          CALL PUSHCONTROL2B(2)
                        ELSE
                          y27 = a4(1, i, k)
                          CALL PUSHCONTROL2B(3)
                        END IF
                        IF (a4(2, i, k) .LT. y27) THEN
                          x9 = y27
                          CALL PUSHCONTROL1B(0)
                        ELSE
                          x9 = a4(2, i, k)
                          CALL PUSHCONTROL1B(1)
                        END IF
                        IF (a4(1, i, k) .LT. pmp_1) THEN
                          IF (pmp_1 .LT. lac_1) THEN
                            y17 = lac_1
                            CALL PUSHCONTROL2B(0)
                          ELSE
                            y17 = pmp_1
                            CALL PUSHCONTROL2B(1)
                          END IF
                        ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                          y17 = lac_1
                          CALL PUSHCONTROL2B(2)
                        ELSE
                          y17 = a4(1, i, k)
                          CALL PUSHCONTROL2B(3)
                        END IF
                        IF (x9 .GT. y17) THEN
                          CALL PUSHREALARRAY_ADM(a4(2, i, k))
                          a4(2, i, k) = y17
                          CALL PUSHCONTROL1B(0)
                        ELSE
                          CALL PUSHREALARRAY_ADM(a4(2, i, k))
                          a4(2, i, k) = x9
                          CALL PUSHCONTROL1B(1)
                        END IF
! Right edges
                        pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                        lac_2 = pmp_2 - 1.5*gam(i, k-1)
                        IF (a4(1, i, k) .GT. pmp_2) THEN
                          IF (pmp_2 .GT. lac_2) THEN
                            y28 = lac_2
                            CALL PUSHCONTROL2B(0)
                          ELSE
                            y28 = pmp_2
                            CALL PUSHCONTROL2B(1)
                          END IF
                        ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                          y28 = lac_2
                          CALL PUSHCONTROL2B(2)
                        ELSE
                          y28 = a4(1, i, k)
                          CALL PUSHCONTROL2B(3)
                        END IF
                        IF (a4(3, i, k) .LT. y28) THEN
                          x10 = y28
                          CALL PUSHCONTROL1B(0)
                        ELSE
                          x10 = a4(3, i, k)
                          CALL PUSHCONTROL1B(1)
                        END IF
                        IF (a4(1, i, k) .LT. pmp_2) THEN
                          IF (pmp_2 .LT. lac_2) THEN
                            y18 = lac_2
                            CALL PUSHCONTROL2B(0)
                          ELSE
                            y18 = pmp_2
                            CALL PUSHCONTROL2B(1)
                          END IF
                        ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                          y18 = lac_2
                          CALL PUSHCONTROL2B(2)
                        ELSE
                          y18 = a4(1, i, k)
                          CALL PUSHCONTROL2B(3)
                        END IF
                        IF (x10 .GT. y18) THEN
                          CALL PUSHREALARRAY_ADM(a4(3, i, k))
                          a4(3, i, k) = y18
                          CALL PUSHCONTROL1B(0)
                        ELSE
                          CALL PUSHREALARRAY_ADM(a4(3, i, k))
                          a4(3, i, k) = x10
                          CALL PUSHCONTROL1B(1)
                        END IF
                        CALL PUSHREALARRAY_ADM(a4(4, i, k))
                        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4&
&                         (3, i, k)))
                        CALL PUSHCONTROL2B(1)
                      END IF
                    ELSE
                      CALL PUSHREALARRAY_ADM(a4(4, i, k))
                      a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3&
&                       , i, k)))
                      CALL PUSHCONTROL2B(0)
                    END IF
                  END DO
                  CALL PUSHCONTROL3B(4)
                ELSE
                  IF (kord .GE. 0.) THEN
                    abs9 = kord
                  ELSE
                    abs9 = -kord
                  END IF
                  IF (abs9 .EQ. 14) THEN
                    DO i=i1,i2
                      CALL PUSHREALARRAY_ADM(a4(4, i, k))
                      a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3&
&                       , i, k)))
                    END DO
                    CALL PUSHCONTROL3B(5)
                  ELSE
! kord = 11
                    DO i=i1,i2
                      IF (extm(i, k) .AND. (extm(i, k-1) .OR. extm(i, k+&
&                         1))) THEN
! Noisy region:
                        CALL PUSHREALARRAY_ADM(a4(2, i, k))
                        a4(2, i, k) = a4(1, i, k)
                        CALL PUSHREALARRAY_ADM(a4(3, i, k))
                        a4(3, i, k) = a4(1, i, k)
                        CALL PUSHREALARRAY_ADM(a4(4, i, k))
                        a4(4, i, k) = 0.
                        CALL PUSHCONTROL1B(1)
                      ELSE
                        CALL PUSHREALARRAY_ADM(a4(4, i, k))
                        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4&
&                         (3, i, k)))
                        CALL PUSHCONTROL1B(0)
                      END IF
                    END DO
                    CALL PUSHCONTROL3B(6)
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
! Additional constraint to ensure positivity
        IF (iv .EQ. 0) THEN
          CALL CS_LIMITERS_FWD(im, extm(i1, k), a4(1, i1, k), 0)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
! k-loop
!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
      IF (iv .EQ. 0) THEN
        DO i=i1,i2
          IF (0. .LT. a4(3, i, km)) THEN
            CALL PUSHREALARRAY_ADM(a4(3, i, km))
            a4(3, i, km) = a4(3, i, km)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREALARRAY_ADM(a4(3, i, km))
            a4(3, i, km) = 0.
            CALL PUSHCONTROL1B(1)
          END IF
        END DO
        CALL PUSHCONTROL2B(2)
      ELSE IF (iv .EQ. -1) THEN
        DO i=i1,i2
          IF (a4(3, i, km)*a4(1, i, km) .LE. 0.) THEN
            CALL PUSHREALARRAY_ADM(a4(3, i, km))
            a4(3, i, km) = 0.
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
      DO k=km-1,km
        DO i=i1,i2
          CALL PUSHREALARRAY_ADM(a4(4, i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
        IF (k .EQ. km - 1) THEN
          CALL CS_LIMITERS_FWD(im, extm(i1, k), a4(1, i1, k), 2)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (k .EQ. km) THEN
          CALL CS_LIMITERS_FWD(im, extm(i1, k), a4(1, i1, k), 1)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
      DO k=km,km-1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) CALL CS_LIMITERS_BWD(im, extm(i1, k), a4(1, &
&                                         i1, k), a4_ad(1, i1, k), 1)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) CALL CS_LIMITERS_BWD(im, extm(i1, k), a4(1, &
&                                         i1, k), a4_ad(1, i1, k), 2)
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(a4(4, i, k))
          temp_ad22 = 3.*a4_ad(4, i, k)
          a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad22
          a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad22
          a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad22
          a4_ad(4, i, k) = 0.0
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .NE. 0) THEN
        IF (branch .EQ. 1) THEN
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) THEN
              CALL POPREALARRAY_ADM(a4(3, i, km))
              a4_ad(3, i, km) = 0.0
            END IF
          END DO
        ELSE
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(a4(3, i, km))
            ELSE
              CALL POPREALARRAY_ADM(a4(3, i, km))
              a4_ad(3, i, km) = 0.0
            END IF
          END DO
        END IF
      END IF
      gam_ad = 0.0
      DO k=km-2,3,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) CALL CS_LIMITERS_BWD(im, extm(i1, k), a4(1, &
&                                         i1, k), a4_ad(1, i1, k), 0)
        CALL POPCONTROL3B(branch)
        IF (branch .LT. 3) THEN
          IF (branch .EQ. 0) THEN
            DO i=i2,i1,-1
              CALL POPREALARRAY_ADM(a4(4, i, k))
              temp_ad17 = 3.*a4_ad(4, i, k)
              a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad17
              a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad17
              a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad17
              a4_ad(4, i, k) = 0.0
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(a4(3, i, k))
                y10_ad = a4_ad(3, i, k)
                a4_ad(3, i, k) = 0.0
                x2_ad = 0.0
              ELSE
                CALL POPREALARRAY_ADM(a4(3, i, k))
                x2_ad = a4_ad(3, i, k)
                a4_ad(3, i, k) = 0.0
                y10_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_2_ad = y10_ad
                  pmp_2_ad = 0.0
                ELSE
                  pmp_2_ad = y10_ad
                  lac_2_ad = 0.0
                END IF
              ELSE
                IF (branch .EQ. 2) THEN
                  lac_2_ad = y10_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y10_ad
                  lac_2_ad = 0.0
                END IF
                pmp_2_ad = 0.0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                y20_ad = x2_ad
              ELSE
                a4_ad(3, i, k) = a4_ad(3, i, k) + x2_ad
                y20_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_2_ad = lac_2_ad + y20_ad
                ELSE
                  pmp_2_ad = pmp_2_ad + y20_ad
                END IF
              ELSE IF (branch .EQ. 2) THEN
                lac_2_ad = lac_2_ad + y20_ad
              ELSE
                a4_ad(1, i, k) = a4_ad(1, i, k) + y20_ad
              END IF
              pmp_2_ad = pmp_2_ad + lac_2_ad
              gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
              a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
              gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(a4(2, i, k))
                y9_ad = a4_ad(2, i, k)
                a4_ad(2, i, k) = 0.0
                x1_ad = 0.0
              ELSE
                CALL POPREALARRAY_ADM(a4(2, i, k))
                x1_ad = a4_ad(2, i, k)
                a4_ad(2, i, k) = 0.0
                y9_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_1_ad = y9_ad
                  pmp_1_ad = 0.0
                ELSE
                  pmp_1_ad = y9_ad
                  lac_1_ad = 0.0
                END IF
              ELSE
                IF (branch .EQ. 2) THEN
                  lac_1_ad = y9_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y9_ad
                  lac_1_ad = 0.0
                END IF
                pmp_1_ad = 0.0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                y19_ad = x1_ad
              ELSE
                a4_ad(2, i, k) = a4_ad(2, i, k) + x1_ad
                y19_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_1_ad = lac_1_ad + y19_ad
                ELSE
                  pmp_1_ad = pmp_1_ad + y19_ad
                END IF
              ELSE IF (branch .EQ. 2) THEN
                lac_1_ad = lac_1_ad + y19_ad
              ELSE
                a4_ad(1, i, k) = a4_ad(1, i, k) + y19_ad
              END IF
              pmp_1_ad = pmp_1_ad + lac_1_ad
              gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
              a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
              gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
            END DO
          ELSE IF (branch .EQ. 1) THEN
            DO i=i2,i1,-1
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .NE. 0) THEN
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                  a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(4, i, k) = 0.0
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    y12_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    x4_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    x4_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    y12_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = y12_ad
                      pmp_2_ad = 0.0
                    ELSE
                      pmp_2_ad = y12_ad
                      lac_2_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_2_ad = y12_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y12_ad
                      lac_2_ad = 0.0
                    END IF
                    pmp_2_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y22_ad = x4_ad
                  ELSE
                    a4_ad(3, i, k) = a4_ad(3, i, k) + x4_ad
                    y22_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = lac_2_ad + y22_ad
                    ELSE
                      pmp_2_ad = pmp_2_ad + y22_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_2_ad = lac_2_ad + y22_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y22_ad
                  END IF
                  pmp_2_ad = pmp_2_ad + lac_2_ad
                  gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                  gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    y11_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    x3_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    x3_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    y11_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = y11_ad
                      pmp_1_ad = 0.0
                    ELSE
                      pmp_1_ad = y11_ad
                      lac_1_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_1_ad = y11_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y11_ad
                      lac_1_ad = 0.0
                    END IF
                    pmp_1_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y21_ad = x3_ad
                  ELSE
                    a4_ad(2, i, k) = a4_ad(2, i, k) + x3_ad
                    y21_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = lac_1_ad + y21_ad
                    ELSE
                      pmp_1_ad = pmp_1_ad + y21_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_1_ad = lac_1_ad + y21_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y21_ad
                  END IF
                  pmp_1_ad = pmp_1_ad + lac_1_ad
                  gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                  gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
                END IF
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(4, i, k) = 0.0
              ELSE IF (branch .EQ. 2) THEN
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(4, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(3, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                a4_ad(3, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(2, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                a4_ad(2, i, k) = 0.0
              ELSE
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(4, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(3, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                a4_ad(3, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(2, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                a4_ad(2, i, k) = 0.0
              END IF
            END DO
          ELSE
            DO i=i2,i1,-1
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .NE. 0) THEN
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                  a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(4, i, k) = 0.0
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    y14_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    x6_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    x6_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    y14_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = y14_ad
                      pmp_2_ad = 0.0
                    ELSE
                      pmp_2_ad = y14_ad
                      lac_2_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_2_ad = y14_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y14_ad
                      lac_2_ad = 0.0
                    END IF
                    pmp_2_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y24_ad = x6_ad
                  ELSE
                    a4_ad(3, i, k) = a4_ad(3, i, k) + x6_ad
                    y24_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = lac_2_ad + y24_ad
                    ELSE
                      pmp_2_ad = pmp_2_ad + y24_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_2_ad = lac_2_ad + y24_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y24_ad
                  END IF
                  pmp_2_ad = pmp_2_ad + lac_2_ad
                  gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                  gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    y13_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    x5_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    x5_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    y13_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = y13_ad
                      pmp_1_ad = 0.0
                    ELSE
                      pmp_1_ad = y13_ad
                      lac_1_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_1_ad = y13_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y13_ad
                      lac_1_ad = 0.0
                    END IF
                    pmp_1_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y23_ad = x5_ad
                  ELSE
                    a4_ad(2, i, k) = a4_ad(2, i, k) + x5_ad
                    y23_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = lac_1_ad + y23_ad
                    ELSE
                      pmp_1_ad = pmp_1_ad + y23_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_1_ad = lac_1_ad + y23_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y23_ad
                  END IF
                  pmp_1_ad = pmp_1_ad + lac_1_ad
                  gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                  gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
                END IF
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(4, i, k) = 0.0
              ELSE IF (branch .EQ. 2) THEN
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                a4_ad(4, i, k) = 0.0
              ELSE
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(4, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(3, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                a4_ad(3, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(2, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                a4_ad(2, i, k) = 0.0
              END IF
            END DO
          END IF
        ELSE IF (branch .LT. 5) THEN
          IF (branch .EQ. 3) THEN
            DO 100 i=i2,i1,-1
              CALL POPCONTROL2B(branch)
              IF (branch .NE. 0) THEN
                IF (branch .EQ. 1) THEN
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
                  a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
                  a4_ad(4, i, k) = 0.0
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    y16_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    x8_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(3, i, k))
                    x8_ad = a4_ad(3, i, k)
                    a4_ad(3, i, k) = 0.0
                    y16_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = y16_ad
                      pmp_2_ad = 0.0
                    ELSE
                      pmp_2_ad = y16_ad
                      lac_2_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_2_ad = y16_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y16_ad
                      lac_2_ad = 0.0
                    END IF
                    pmp_2_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y26_ad = x8_ad
                  ELSE
                    a4_ad(3, i, k) = a4_ad(3, i, k) + x8_ad
                    y26_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_2_ad = lac_2_ad + y26_ad
                    ELSE
                      pmp_2_ad = pmp_2_ad + y26_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_2_ad = lac_2_ad + y26_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y26_ad
                  END IF
                  pmp_2_ad = pmp_2_ad + lac_2_ad
                  gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                  gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    y15_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    x7_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(a4(2, i, k))
                    x7_ad = a4_ad(2, i, k)
                    a4_ad(2, i, k) = 0.0
                    y15_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = y15_ad
                      pmp_1_ad = 0.0
                    ELSE
                      pmp_1_ad = y15_ad
                      lac_1_ad = 0.0
                    END IF
                  ELSE
                    IF (branch .EQ. 2) THEN
                      lac_1_ad = y15_ad
                    ELSE
                      a4_ad(1, i, k) = a4_ad(1, i, k) + y15_ad
                      lac_1_ad = 0.0
                    END IF
                    pmp_1_ad = 0.0
                  END IF
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    y25_ad = x7_ad
                  ELSE
                    a4_ad(2, i, k) = a4_ad(2, i, k) + x7_ad
                    y25_ad = 0.0
                  END IF
                  CALL POPCONTROL2B(branch)
                  IF (branch .LT. 2) THEN
                    IF (branch .EQ. 0) THEN
                      lac_1_ad = lac_1_ad + y25_ad
                    ELSE
                      pmp_1_ad = pmp_1_ad + y25_ad
                    END IF
                  ELSE IF (branch .EQ. 2) THEN
                    lac_1_ad = lac_1_ad + y25_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y25_ad
                  END IF
                  pmp_1_ad = pmp_1_ad + lac_1_ad
                  gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                  a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                  gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
                ELSE
                  CALL POPREALARRAY_ADM(a4(4, i, k))
                  a4_ad(4, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  GOTO 100
                END IF
              END IF
              CALL POPREALARRAY_ADM(a4(4, i, k))
              a4_ad(1, i, k) = a4_ad(1, i, k) + 6.*a4_ad(4, i, k)
              a4_ad(2, i, k) = a4_ad(2, i, k) - 3.*a4_ad(4, i, k)
              a4_ad(3, i, k) = a4_ad(3, i, k) - 3.*a4_ad(4, i, k)
              a4_ad(4, i, k) = 0.0
 100        CONTINUE
          ELSE
            DO i=i2,i1,-1
              CALL POPCONTROL2B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(a4(4, i, k))
                temp_ad19 = 3.*a4_ad(4, i, k)
                a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad19
                a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad19
                a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad19
                a4_ad(4, i, k) = 0.0
              ELSE IF (branch .EQ. 1) THEN
                CALL POPREALARRAY_ADM(a4(4, i, k))
                temp_ad18 = 3.*a4_ad(4, i, k)
                a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad18
                a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad18
                a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad18
                a4_ad(4, i, k) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  y18_ad = a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  x10_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(a4(3, i, k))
                  x10_ad = a4_ad(3, i, k)
                  a4_ad(3, i, k) = 0.0
                  y18_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = y18_ad
                    pmp_2_ad = 0.0
                  ELSE
                    pmp_2_ad = y18_ad
                    lac_2_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_2_ad = y18_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y18_ad
                    lac_2_ad = 0.0
                  END IF
                  pmp_2_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y28_ad = x10_ad
                ELSE
                  a4_ad(3, i, k) = a4_ad(3, i, k) + x10_ad
                  y28_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = lac_2_ad + y28_ad
                  ELSE
                    pmp_2_ad = pmp_2_ad + y28_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_2_ad = lac_2_ad + y28_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y28_ad
                END IF
                pmp_2_ad = pmp_2_ad + lac_2_ad
                gam_ad(i, k-1) = gam_ad(i, k-1) - 1.5*lac_2_ad
                a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_2_ad
                gam_ad(i, k) = gam_ad(i, k) + 2.*pmp_2_ad
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  y17_ad = a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  x9_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(a4(2, i, k))
                  x9_ad = a4_ad(2, i, k)
                  a4_ad(2, i, k) = 0.0
                  y17_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = y17_ad
                    pmp_1_ad = 0.0
                  ELSE
                    pmp_1_ad = y17_ad
                    lac_1_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_1_ad = y17_ad
                  ELSE
                    a4_ad(1, i, k) = a4_ad(1, i, k) + y17_ad
                    lac_1_ad = 0.0
                  END IF
                  pmp_1_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y27_ad = x9_ad
                ELSE
                  a4_ad(2, i, k) = a4_ad(2, i, k) + x9_ad
                  y27_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = lac_1_ad + y27_ad
                  ELSE
                    pmp_1_ad = pmp_1_ad + y27_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_1_ad = lac_1_ad + y27_ad
                ELSE
                  a4_ad(1, i, k) = a4_ad(1, i, k) + y27_ad
                END IF
                pmp_1_ad = pmp_1_ad + lac_1_ad
                gam_ad(i, k+2) = gam_ad(i, k+2) + 1.5*lac_1_ad
                a4_ad(1, i, k) = a4_ad(1, i, k) + pmp_1_ad
                gam_ad(i, k+1) = gam_ad(i, k+1) - 2.*pmp_1_ad
              ELSE
                CALL POPREALARRAY_ADM(a4(4, i, k))
                a4_ad(4, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(3, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
                a4_ad(3, i, k) = 0.0
                CALL POPREALARRAY_ADM(a4(2, i, k))
                a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
                a4_ad(2, i, k) = 0.0
              END IF
            END DO
          END IF
        ELSE IF (branch .EQ. 5) THEN
          DO i=i2,i1,-1
            CALL POPREALARRAY_ADM(a4(4, i, k))
            temp_ad20 = 3.*a4_ad(4, i, k)
            a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad20
            a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad20
            a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad20
            a4_ad(4, i, k) = 0.0
          END DO
        ELSE
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(a4(4, i, k))
              temp_ad21 = 3.*a4_ad(4, i, k)
              a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad21
              a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad21
              a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad21
              a4_ad(4, i, k) = 0.0
            ELSE
              CALL POPREALARRAY_ADM(a4(4, i, k))
              a4_ad(4, i, k) = 0.0
              CALL POPREALARRAY_ADM(a4(3, i, k))
              a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(3, i, k)
              a4_ad(3, i, k) = 0.0
              CALL POPREALARRAY_ADM(a4(2, i, k))
              a4_ad(1, i, k) = a4_ad(1, i, k) + a4_ad(2, i, k)
              a4_ad(2, i, k) = 0.0
            END IF
          END DO
        END IF
      END DO
      CALL CS_LIMITERS_BWD(im, extm(i1, 2), a4(1, i1, 2), a4_ad(1, i1, 2&
&                    ), 2)
      DO i=i2,i1,-1
        CALL POPREALARRAY_ADM(a4(4, i, 2))
        temp_ad16 = 3.*a4_ad(4, i, 2)
        a4_ad(1, i, 2) = a4_ad(1, i, 2) + 2.*temp_ad16
        a4_ad(2, i, 2) = a4_ad(2, i, 2) - temp_ad16
        a4_ad(3, i, 2) = a4_ad(3, i, 2) - temp_ad16
        a4_ad(4, i, 2) = 0.0
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        CALL CS_LIMITERS_BWD(im, extm(i1, 1), a4(1, i1, 1), a4_ad(1, i1&
&                      , 1), 1)
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(a4(4, i, 1))
          temp_ad15 = 3.*a4_ad(4, i, 1)
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + 2.*temp_ad15
          a4_ad(2, i, 1) = a4_ad(2, i, 1) - temp_ad15
          a4_ad(3, i, 1) = a4_ad(3, i, 1) - temp_ad15
          a4_ad(4, i, 1) = 0.0
        END DO
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(a4(2, i, 1))
            ELSE
              CALL POPREALARRAY_ADM(a4(2, i, 1))
              a4_ad(2, i, 1) = 0.0
            END IF
          END DO
        ELSE
          DO i=i2,i1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) THEN
              CALL POPREALARRAY_ADM(a4(2, i, 1))
              a4_ad(2, i, 1) = 0.0
            END IF
          END DO
        END IF
      ELSE IF (branch .EQ. 2) THEN
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(a4(4, i, 1))
          a4_ad(4, i, 1) = 0.0
          CALL POPREALARRAY_ADM(a4(3, i, 1))
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + a4_ad(3, i, 1)
          a4_ad(3, i, 1) = 0.0
          CALL POPREALARRAY_ADM(a4(2, i, 1))
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + a4_ad(2, i, 1)
          a4_ad(2, i, 1) = 0.0
        END DO
      END IF
      DO k=km,1,-1
        CALL POPCONTROL1B(branch)
      END DO
      q_ad = 0.0
      DO k=km,1,-1
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(a4(3, i, k))
          q_ad(i, k+1) = q_ad(i, k+1) + a4_ad(3, i, k)
          a4_ad(3, i, k) = 0.0
          CALL POPREALARRAY_ADM(a4(2, i, k))
          q_ad(i, k) = q_ad(i, k) + a4_ad(2, i, k)
          a4_ad(2, i, k) = 0.0
        END DO
      END DO
      DO i=i2,i1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y8_ad = q_ad(i, km)
          q_ad(i, km) = 0.0
        ELSE
          y8_ad = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          a4_ad(1, i, km) = a4_ad(1, i, km) + y8_ad
        ELSE
          a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + y8_ad
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY_ADM(q(i, km))
          y7_ad = q_ad(i, km)
          q_ad(i, km) = 0.0
        ELSE
          CALL POPREALARRAY_ADM(q(i, km))
          y7_ad = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          a4_ad(1, i, km) = a4_ad(1, i, km) + y7_ad
        ELSE
          a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + y7_ad
        END IF
      END DO
      DO k=km-1,3,-1
        DO 120 i=i2,i1,-1
          CALL POPCONTROL3B(branch)
          IF (branch .NE. 0) THEN
            IF (branch .LT. 4) THEN
              IF (branch .EQ. 1) THEN
                GOTO 110
              ELSE IF (branch .EQ. 2) THEN
                q_ad(i, k) = 0.0
                GOTO 110
              ELSE
                CALL POPREALARRAY_ADM(q(i, k))
                y5_ad = q_ad(i, k)
                q_ad(i, k) = 0.0
              END IF
            ELSE IF (branch .EQ. 4) THEN
              CALL POPREALARRAY_ADM(q(i, k))
              y5_ad = 0.0
            ELSE
              IF (branch .EQ. 5) THEN
                y4_ad = q_ad(i, k)
                q_ad(i, k) = 0.0
              ELSE
                y4_ad = 0.0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                a4_ad(1, i, k) = a4_ad(1, i, k) + y4_ad
              ELSE
                a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + y4_ad
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(q(i, k))
                y3_ad = q_ad(i, k)
                q_ad(i, k) = 0.0
              ELSE
                CALL POPREALARRAY_ADM(q(i, k))
                y3_ad = 0.0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                a4_ad(1, i, k) = a4_ad(1, i, k) + y3_ad
              ELSE
                a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + y3_ad
              END IF
              GOTO 120
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              a4_ad(1, i, k) = a4_ad(1, i, k) + y5_ad
            ELSE
              a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + y5_ad
            END IF
            GOTO 120
          END IF
 110      CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY_ADM(q(i, k))
            y6_ad = q_ad(i, k)
            q_ad(i, k) = 0.0
          ELSE
            CALL POPREALARRAY_ADM(q(i, k))
            y6_ad = 0.0
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            a4_ad(1, i, k) = a4_ad(1, i, k) + y6_ad
          ELSE
            a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + y6_ad
          END IF
 120    CONTINUE
      END DO
      DO k=km,2,-1
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(gam(i, k))
          a4_ad(1, i, k) = a4_ad(1, i, k) + gam_ad(i, k)
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) - gam_ad(i, k)
          gam_ad(i, k) = 0.0
        END DO
      END DO
      DO i=i2,i1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y2_ad = q_ad(i, 2)
          q_ad(i, 2) = 0.0
        ELSE
          y2_ad = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          a4_ad(1, i, 2) = a4_ad(1, i, 2) + y2_ad
        ELSE
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + y2_ad
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY_ADM(q(i, 2))
          y1_ad = q_ad(i, 2)
          q_ad(i, 2) = 0.0
        ELSE
          CALL POPREALARRAY_ADM(q(i, 2))
          y1_ad = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          a4_ad(1, i, 2) = a4_ad(1, i, 2) + y1_ad
        ELSE
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + y1_ad
        END IF
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO k=1,km,1
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(q(i, k))
          gam_ad(i, k) = gam_ad(i, k) - q(i, k+1)*q_ad(i, k)
          q_ad(i, k+1) = q_ad(i, k+1) - gam(i, k)*q_ad(i, k)
        END DO
      END DO
      d4_ad = 0.0
      DO i=i2,i1,-1
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        CALL POPREALARRAY_ADM(q(i, km+1))
        temp2 = d4(i)*(d4(i)+0.5) - a_bot*gam(i, km)
        temp_ad11 = q_ad(i, km+1)/temp2
        temp1 = d4(i)*(d4(i)+1.)
        temp_ad12 = 2.*a4(1, i, km)*temp_ad11
        temp_ad13 = -((2.*(temp1*a4(1, i, km))+a4(1, i, km-1)-a_bot*q(i&
&         , km))*temp_ad11/temp2)
        a4_ad(1, i, km) = a4_ad(1, i, km) + 2.*temp1*temp_ad11
        a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + temp_ad11
        a_bot_ad = -(gam(i, km)*temp_ad13) - q(i, km)*temp_ad11
        d4_ad(i) = d4_ad(i) + (2*d4(i)+1.5)*a_bot_ad + (2*d4(i)+0.5)*&
&         temp_ad13 + (2*d4(i)+1.)*temp_ad12
        q_ad(i, km) = q_ad(i, km) - a_bot*temp_ad11
        gam_ad(i, km) = gam_ad(i, km) - a_bot*temp_ad13
        q_ad(i, km+1) = 0.0
      END DO
      DO k=km,2,-1
        DO i=i2,i1,-1
          temp_ad9 = q_ad(i, k)/bet
          temp_ad8 = 3.*temp_ad9
          CALL POPREALARRAY_ADM(q(i, k))
          bet_ad = -((3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))*&
&           temp_ad9/bet) - d4(i)*gam_ad(i, k)/bet**2
          d4_ad(i) = d4_ad(i) + a4(1, i, k)*temp_ad8 + 2*bet_ad + gam_ad&
&           (i, k)/bet
          gam_ad(i, k) = 0.0
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + temp_ad8
          a4_ad(1, i, k) = a4_ad(1, i, k) + d4(i)*temp_ad8
          q_ad(i, k-1) = q_ad(i, k-1) - temp_ad9
          q_ad(i, k) = 0.0
          CALL POPREALARRAY_ADM(bet)
          gam_ad(i, k-1) = gam_ad(i, k-1) - bet_ad
          CALL POPREALARRAY_ADM(d4(i))
          temp_ad10 = d4_ad(i)/delp(i, k)
          delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad10
          delp_ad(i, k) = delp_ad(i, k) - delp(i, k-1)*temp_ad10/delp(i&
&           , k)
          d4_ad(i) = 0.0
        END DO
      END DO
      DO i=i2,i1,-1
        grat = delp(i, 2)/delp(i, 1)
        bet = grat*(grat+0.5)
        temp_ad4 = gam_ad(i, 1)/bet
        gam_ad(i, 1) = 0.0
        temp_ad6 = q_ad(i, 1)/bet
        temp_ad5 = a4(1, i, 1)*temp_ad6
        temp0 = 2*grat*(grat+1.)
        bet_ad = -((temp0*a4(1, i, 1)+a4(1, i, 2))*temp_ad6/bet) - (grat&
&         *(grat+1.5)+1.)*temp_ad4/bet
        grat_ad = (4*grat+2*1.)*temp_ad5 + (2*grat+0.5)*bet_ad + (2*grat&
&         +1.5)*temp_ad4
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + temp0*temp_ad6
        a4_ad(1, i, 2) = a4_ad(1, i, 2) + temp_ad6
        q_ad(i, 1) = 0.0
        temp_ad7 = grat_ad/delp(i, 1)
        delp_ad(i, 2) = delp_ad(i, 2) + temp_ad7
        delp_ad(i, 1) = delp_ad(i, 1) - delp(i, 2)*temp_ad7/delp(i, 1)
      END DO
    ELSE
      DO k=1,km-1,1
        DO i=i2,i1,-1
          CALL POPREALARRAY_ADM(q(i, k))
          gam_ad(i, k+1) = gam_ad(i, k+1) - q(i, k+1)*q_ad(i, k)
          q_ad(i, k+1) = q_ad(i, k+1) - gam(i, k+1)*q_ad(i, k)
        END DO
      END DO
      DO i=i2,i1,-1
        CALL POPREALARRAY_ADM(q(i, km+1))
        qs_ad(i) = qs_ad(i) + q_ad(i, km+1)
        q_ad(i, km+1) = 0.0
        grat = delp(i, km-1)/delp(i, km)
        CALL POPREALARRAY_ADM(q(i, km))
        temp = 2*grat - gam(i, km) + 2.
        temp_ad1 = q_ad(i, km)/temp
        temp_ad2 = -((3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, &
&         km-1))*temp_ad1/temp)
        a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + 3.*temp_ad1
        a4_ad(1, i, km) = a4_ad(1, i, km) + 3.*temp_ad1
        grat_ad = 2*temp_ad2 - qs(i)*temp_ad1
        qs_ad(i) = qs_ad(i) - grat*temp_ad1
        q_ad(i, km-1) = q_ad(i, km-1) - temp_ad1
        gam_ad(i, km) = gam_ad(i, km) - temp_ad2
        q_ad(i, km) = 0.0
        temp_ad3 = grat_ad/delp(i, km)
        delp_ad(i, km-1) = delp_ad(i, km-1) + temp_ad3
        delp_ad(i, km) = delp_ad(i, km) - delp(i, km-1)*temp_ad3/delp(i&
&         , km)
      END DO
      DO k=km-1,2,-1
        DO i=i2,i1,-1
          temp_ad = q_ad(i, k)/bet
          CALL POPREALARRAY_ADM(q(i, k))
          grat = delp(i, k-1)/delp(i, k)
          bet_ad = -((3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))*temp_ad/&
&           bet) - grat*gam_ad(i, k+1)/bet**2
          grat_ad = 2*bet_ad + gam_ad(i, k+1)/bet
          gam_ad(i, k+1) = 0.0
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + 3.*temp_ad
          a4_ad(1, i, k) = a4_ad(1, i, k) + 3.*temp_ad
          q_ad(i, k-1) = q_ad(i, k-1) - temp_ad
          q_ad(i, k) = 0.0
          CALL POPREALARRAY_ADM(bet)
          gam_ad(i, k) = gam_ad(i, k) - bet_ad
          temp_ad0 = grat_ad/delp(i, k)
          delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad0
          delp_ad(i, k) = delp_ad(i, k) - delp(i, k-1)*temp_ad0/delp(i, &
&           k)
        END DO
      END DO
      DO i=i2,i1,-1
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + 1.5*q_ad(i, 1)
        q_ad(i, 1) = 0.0
      END DO
    END IF
  END SUBROUTINE CS_PROFILE_ADM
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
!  Differentiation of cs_limiters in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
!mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core
!_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.
!mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleig
!h_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_or
!d4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.r
!emap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d 
!fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters
! fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_
!restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgr
!id_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils
!_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_m
!od.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mo
!d.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d
!2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
!_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core
!_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_util
!s_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: a4
!   with respect to varying inputs: a4
  SUBROUTINE CS_LIMITERS_FWD(im, extm, a4, iv)
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
          CALL PUSHREALARRAY(a4(2, i))
          a4(2, i) = a4(1, i)
          CALL PUSHREALARRAY(a4(3, i))
          a4(3, i) = a4(1, i)
          CALL PUSHREALARRAY(a4(4, i))
          a4(4, i) = 0.
          CALL PUSHCONTROL(3,5)
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
                CALL PUSHREALARRAY(a4(3, i))
                a4(3, i) = a4(1, i)
                CALL PUSHREALARRAY(a4(2, i))
                a4(2, i) = a4(1, i)
                CALL PUSHREALARRAY(a4(4, i))
                a4(4, i) = 0.
                CALL PUSHCONTROL(3,4)
              ELSE IF (a4(3, i) .GT. a4(2, i)) THEN
                CALL PUSHREALARRAY(a4(4, i))
                a4(4, i) = 3.*(a4(2, i)-a4(1, i))
                CALL PUSHREALARRAY(a4(3, i))
                a4(3, i) = a4(2, i) - a4(4, i)
                CALL PUSHCONTROL(3,3)
              ELSE
                CALL PUSHREALARRAY(a4(4, i))
                a4(4, i) = 3.*(a4(3, i)-a4(1, i))
                CALL PUSHREALARRAY(a4(2, i))
                a4(2, i) = a4(3, i) - a4(4, i)
                CALL PUSHCONTROL(3,2)
              END IF
            ELSE
              CALL PUSHCONTROL(3,1)
            END IF
          ELSE
            CALL PUSHCONTROL(3,0)
          END IF
        END IF
      END DO
      CALL PUSHCONTROL(2,0)
    ELSE IF (iv .EQ. 1) THEN
      DO i=1,im
        IF ((a4(1, i)-a4(2, i))*(a4(1, i)-a4(3, i)) .GE. 0.) THEN
          CALL PUSHREALARRAY(a4(2, i))
          a4(2, i) = a4(1, i)
          CALL PUSHREALARRAY(a4(3, i))
          a4(3, i) = a4(1, i)
          CALL PUSHREALARRAY(a4(4, i))
          a4(4, i) = 0.
          CALL PUSHCONTROL(2,3)
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            CALL PUSHREALARRAY(a4(4, i))
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            CALL PUSHREALARRAY(a4(3, i))
            a4(3, i) = a4(2, i) - a4(4, i)
            CALL PUSHCONTROL(2,2)
          ELSE IF (a6da .GT. da2) THEN
            CALL PUSHREALARRAY(a4(4, i))
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            CALL PUSHREALARRAY(a4(2, i))
            a4(2, i) = a4(3, i) - a4(4, i)
            CALL PUSHCONTROL(2,1)
          ELSE
            CALL PUSHCONTROL(2,0)
          END IF
        END IF
      END DO
      CALL PUSHCONTROL(2,1)
    ELSE
! Standard PPM constraint
      DO i=1,im
        IF (extm(i)) THEN
          CALL PUSHREALARRAY(a4(2, i))
          a4(2, i) = a4(1, i)
          CALL PUSHREALARRAY(a4(3, i))
          a4(3, i) = a4(1, i)
          CALL PUSHREALARRAY(a4(4, i))
          a4(4, i) = 0.
          CALL PUSHCONTROL(2,3)
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            CALL PUSHREALARRAY(a4(4, i))
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            CALL PUSHREALARRAY(a4(3, i))
            a4(3, i) = a4(2, i) - a4(4, i)
            CALL PUSHCONTROL(2,2)
          ELSE IF (a6da .GT. da2) THEN
            CALL PUSHREALARRAY(a4(4, i))
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            CALL PUSHREALARRAY(a4(2, i))
            a4(2, i) = a4(3, i) - a4(4, i)
            CALL PUSHCONTROL(2,1)
          ELSE
            CALL PUSHCONTROL(2,0)
          END IF
        END IF
      END DO
      CALL PUSHCONTROL(2,2)
    END IF
  END SUBROUTINE CS_LIMITERS_FWD
!  Differentiation of cs_limiters in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
!_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_cor
!e_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod
!.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Raylei
!gh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_o
!rd4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.
!remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d
! fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiter
!s fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv
!_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subg
!rid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_util
!s_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_
!mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_m
!od.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.
!d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_
!v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_cor
!e_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_uti
!ls_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: a4
!   with respect to varying inputs: a4
  SUBROUTINE CS_LIMITERS_BWD(im, extm, a4, a4_ad, iv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    INTEGER, INTENT(IN) :: iv
    LOGICAL, INTENT(IN) :: extm(im)
    REAL, INTENT(INOUT) :: a4(4, im)
    REAL, INTENT(INOUT) :: a4_ad(4, im)
    REAL :: da1, da2, a6da
    INTEGER :: i
    INTRINSIC ABS
    REAL :: abs0
    INTEGER :: branch
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      DO i=im,1,-1
        CALL POPCONTROL(3,branch)
        IF (branch .LT. 3) THEN
          IF (branch .NE. 0) THEN
            IF (branch .NE. 1) THEN
              CALL POPREALARRAY(a4(2, i))
              a4_ad(3, i) = a4_ad(3, i) + a4_ad(2, i)
              a4_ad(4, i) = a4_ad(4, i) - a4_ad(2, i)
              a4_ad(2, i) = 0.0
              CALL POPREALARRAY(a4(4, i))
              a4_ad(3, i) = a4_ad(3, i) + 3.*a4_ad(4, i)
              a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
              a4_ad(4, i) = 0.0
            END IF
          END IF
        ELSE IF (branch .EQ. 3) THEN
          CALL POPREALARRAY(a4(3, i))
          a4_ad(2, i) = a4_ad(2, i) + a4_ad(3, i)
          a4_ad(4, i) = a4_ad(4, i) - a4_ad(3, i)
          a4_ad(3, i) = 0.0
          CALL POPREALARRAY(a4(4, i))
          a4_ad(2, i) = a4_ad(2, i) + 3.*a4_ad(4, i)
          a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
          a4_ad(4, i) = 0.0
        ELSE IF (branch .EQ. 4) THEN
          CALL POPREALARRAY(a4(4, i))
          a4_ad(4, i) = 0.0
          CALL POPREALARRAY(a4(2, i))
          a4_ad(1, i) = a4_ad(1, i) + a4_ad(2, i)
          a4_ad(2, i) = 0.0
          CALL POPREALARRAY(a4(3, i))
          a4_ad(1, i) = a4_ad(1, i) + a4_ad(3, i)
          a4_ad(3, i) = 0.0
        ELSE
          CALL POPREALARRAY(a4(4, i))
          a4_ad(4, i) = 0.0
          CALL POPREALARRAY(a4(3, i))
          a4_ad(1, i) = a4_ad(1, i) + a4_ad(3, i)
          a4_ad(3, i) = 0.0
          CALL POPREALARRAY(a4(2, i))
          a4_ad(1, i) = a4_ad(1, i) + a4_ad(2, i)
          a4_ad(2, i) = 0.0
        END IF
      END DO
    ELSE IF (branch .EQ. 1) THEN
      DO i=im,1,-1
        CALL POPCONTROL(2,branch)
        IF (branch .LT. 2) THEN
          IF (branch .NE. 0) THEN
            CALL POPREALARRAY(a4(2, i))
            a4_ad(3, i) = a4_ad(3, i) + a4_ad(2, i)
            a4_ad(4, i) = a4_ad(4, i) - a4_ad(2, i)
            a4_ad(2, i) = 0.0
            CALL POPREALARRAY(a4(4, i))
            a4_ad(3, i) = a4_ad(3, i) + 3.*a4_ad(4, i)
            a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
            a4_ad(4, i) = 0.0
          END IF
        ELSE IF (branch .EQ. 2) THEN
          CALL POPREALARRAY(a4(3, i))
          a4_ad(2, i) = a4_ad(2, i) + a4_ad(3, i)
          a4_ad(4, i) = a4_ad(4, i) - a4_ad(3, i)
          a4_ad(3, i) = 0.0
          CALL POPREALARRAY(a4(4, i))
          a4_ad(2, i) = a4_ad(2, i) + 3.*a4_ad(4, i)
          a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
          a4_ad(4, i) = 0.0
        ELSE
          CALL POPREALARRAY(a4(4, i))
          a4_ad(4, i) = 0.0
          CALL POPREALARRAY(a4(3, i))
          a4_ad(1, i) = a4_ad(1, i) + a4_ad(3, i)
          a4_ad(3, i) = 0.0
          CALL POPREALARRAY(a4(2, i))
          a4_ad(1, i) = a4_ad(1, i) + a4_ad(2, i)
          a4_ad(2, i) = 0.0
        END IF
      END DO
    ELSE
      DO i=im,1,-1
        CALL POPCONTROL(2,branch)
        IF (branch .LT. 2) THEN
          IF (branch .NE. 0) THEN
            CALL POPREALARRAY(a4(2, i))
            a4_ad(3, i) = a4_ad(3, i) + a4_ad(2, i)
            a4_ad(4, i) = a4_ad(4, i) - a4_ad(2, i)
            a4_ad(2, i) = 0.0
            CALL POPREALARRAY(a4(4, i))
            a4_ad(3, i) = a4_ad(3, i) + 3.*a4_ad(4, i)
            a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
            a4_ad(4, i) = 0.0
          END IF
        ELSE IF (branch .EQ. 2) THEN
          CALL POPREALARRAY(a4(3, i))
          a4_ad(2, i) = a4_ad(2, i) + a4_ad(3, i)
          a4_ad(4, i) = a4_ad(4, i) - a4_ad(3, i)
          a4_ad(3, i) = 0.0
          CALL POPREALARRAY(a4(4, i))
          a4_ad(2, i) = a4_ad(2, i) + 3.*a4_ad(4, i)
          a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
          a4_ad(4, i) = 0.0
        ELSE
          CALL POPREALARRAY(a4(4, i))
          a4_ad(4, i) = 0.0
          CALL POPREALARRAY(a4(3, i))
          a4_ad(1, i) = a4_ad(1, i) + a4_ad(3, i)
          a4_ad(3, i) = 0.0
          CALL POPREALARRAY(a4(2, i))
          a4_ad(1, i) = a4_ad(1, i) + a4_ad(2, i)
          a4_ad(2, i) = 0.0
        END IF
      END DO
    END IF
  END SUBROUTINE CS_LIMITERS_BWD
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
!  Differentiation of ppm_profile in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
!mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core
!_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.
!mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleig
!h_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_or
!d4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.r
!emap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d 
!fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters
! fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_
!restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgr
!id_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils
!_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_m
!od.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mo
!d.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d
!2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
!_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core
!_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_util
!s_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: delp a4
!   with respect to varying inputs: delp a4
  SUBROUTINE PPM_PROFILE_FWD(a4, delp, km, i1, i2, iv, kord)
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
          CALL PUSHCONTROL(1,0)
        ELSE
          x1 = -df2(i, k)
          CALL PUSHCONTROL(1,1)
        END IF
        IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
          IF (a4(1, i, k) .LT. a4(1, i, k+1)) THEN
            max1 = a4(1, i, k+1)
            CALL PUSHCONTROL(2,0)
          ELSE
            max1 = a4(1, i, k)
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE IF (a4(1, i, k-1) .LT. a4(1, i, k+1)) THEN
          max1 = a4(1, i, k+1)
          CALL PUSHCONTROL(2,2)
        ELSE
          max1 = a4(1, i, k-1)
          CALL PUSHCONTROL(2,3)
        END IF
        y1 = max1 - a4(1, i, k)
        IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
          IF (a4(1, i, k) .GT. a4(1, i, k+1)) THEN
            min2 = a4(1, i, k+1)
            CALL PUSHCONTROL(2,0)
          ELSE
            min2 = a4(1, i, k)
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE IF (a4(1, i, k-1) .GT. a4(1, i, k+1)) THEN
          min2 = a4(1, i, k+1)
          CALL PUSHCONTROL(2,2)
        ELSE
          min2 = a4(1, i, k-1)
          CALL PUSHCONTROL(2,3)
        END IF
        z1 = a4(1, i, k) - min2
        IF (x1 .GT. y1) THEN
          IF (y1 .GT. z1) THEN
            CALL PUSHREALARRAY(min1)
            min1 = z1
            CALL PUSHCONTROL(2,0)
          ELSE
            CALL PUSHREALARRAY(min1)
            min1 = y1
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE IF (x1 .GT. z1) THEN
          CALL PUSHREALARRAY(min1)
          min1 = z1
          CALL PUSHCONTROL(2,2)
        ELSE
          CALL PUSHREALARRAY(min1)
          min1 = x1
          CALL PUSHCONTROL(2,3)
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
      CALL PUSHREALARRAY(qm)
      qm = (d2*a4(1, i, 1)+d1*a4(1, i, 2))/(d1+d2)
      CALL PUSHREALARRAY(dq)
      dq = 2.*(a4(1, i, 2)-a4(1, i, 1))/(d1+d2)
      CALL PUSHREALARRAY(c1)
      c1 = 4.*(a4(2, i, 3)-qm-d2*dq)/(d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      CALL PUSHREALARRAY(a4(2, i, 2))
      a4(2, i, 2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
      CALL PUSHREALARRAY(a4(2, i, 1))
      a4(2, i, 1) = d1*(2.*c1*d1**2-c3) + a4(2, i, 2)
      IF (a4(1, i, 1) .GT. a4(1, i, 2)) THEN
        y2 = a4(1, i, 2)
        CALL PUSHCONTROL(1,0)
      ELSE
        y2 = a4(1, i, 1)
        CALL PUSHCONTROL(1,1)
      END IF
      IF (a4(2, i, 2) .LT. y2) THEN
        CALL PUSHREALARRAY(a4(2, i, 2))
        a4(2, i, 2) = y2
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHREALARRAY(a4(2, i, 2))
        a4(2, i, 2) = a4(2, i, 2)
        CALL PUSHCONTROL(1,1)
      END IF
      IF (a4(1, i, 1) .LT. a4(1, i, 2)) THEN
        y3 = a4(1, i, 2)
        CALL PUSHCONTROL(1,0)
      ELSE
        y3 = a4(1, i, 1)
        CALL PUSHCONTROL(1,1)
      END IF
      IF (a4(2, i, 2) .GT. y3) THEN
        CALL PUSHREALARRAY(a4(2, i, 2))
        a4(2, i, 2) = y3
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHREALARRAY(a4(2, i, 2))
        a4(2, i, 2) = a4(2, i, 2)
        CALL PUSHCONTROL(1,1)
      END IF
      CALL PUSHREALARRAY(dc(i, 1))
      dc(i, 1) = 0.5*(a4(2, i, 2)-a4(1, i, 1))
    END DO
! Enforce monotonicity  within the top layer
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(2, i, 1)) THEN
          CALL PUSHREALARRAY(a4(2, i, 1))
          a4(2, i, 1) = a4(2, i, 1)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(a4(2, i, 1))
          a4(2, i, 1) = 0.
          CALL PUSHCONTROL(1,1)
        END IF
        IF (0. .LT. a4(2, i, 2)) THEN
          CALL PUSHREALARRAY(a4(2, i, 2))
          a4(2, i, 2) = a4(2, i, 2)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(a4(2, i, 2))
          a4(2, i, 2) = 0.
          CALL PUSHCONTROL(1,1)
        END IF
      END DO
      CALL PUSHCONTROL(2,3)
    ELSE IF (iv .EQ. -1) THEN
      DO i=i1,i2
        IF (a4(2, i, 1)*a4(1, i, 1) .LE. 0.) THEN
          CALL PUSHREALARRAY(a4(2, i, 1))
          a4(2, i, 1) = 0.
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHCONTROL(2,2)
    ELSE
      IF (iv .GE. 0.) THEN
        abs0 = iv
      ELSE
        abs0 = -iv
      END IF
      IF (abs0 .EQ. 2) THEN
        DO i=i1,i2
          CALL PUSHREALARRAY(a4(2, i, 1))
          a4(2, i, 1) = a4(1, i, 1)
          CALL PUSHREALARRAY(a4(3, i, 1))
          a4(3, i, 1) = a4(1, i, 1)
        END DO
        CALL PUSHCONTROL(2,1)
      ELSE
        CALL PUSHCONTROL(2,0)
      END IF
    END IF
! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
    DO i=i1,i2
      d1 = delp(i, km)
      d2 = delp(i, km1)
      CALL PUSHREALARRAY(qm)
      qm = (d2*a4(1, i, km)+d1*a4(1, i, km1))/(d1+d2)
      CALL PUSHREALARRAY(dq)
      dq = 2.*(a4(1, i, km1)-a4(1, i, km))/(d1+d2)
      CALL PUSHREALARRAY(c1)
      c1 = (a4(2, i, km1)-qm-d2*dq)/(d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      CALL PUSHREALARRAY(a4(2, i, km))
      a4(2, i, km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
      CALL PUSHREALARRAY(a4(3, i, km))
      a4(3, i, km) = d1*(8.*c1*d1**2-c3) + a4(2, i, km)
      IF (a4(1, i, km) .GT. a4(1, i, km1)) THEN
        y4 = a4(1, i, km1)
        CALL PUSHCONTROL(1,0)
      ELSE
        y4 = a4(1, i, km)
        CALL PUSHCONTROL(1,1)
      END IF
      IF (a4(2, i, km) .LT. y4) THEN
        CALL PUSHREALARRAY(a4(2, i, km))
        a4(2, i, km) = y4
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHREALARRAY(a4(2, i, km))
        a4(2, i, km) = a4(2, i, km)
        CALL PUSHCONTROL(1,1)
      END IF
      IF (a4(1, i, km) .LT. a4(1, i, km1)) THEN
        y5 = a4(1, i, km1)
        CALL PUSHCONTROL(1,0)
      ELSE
        y5 = a4(1, i, km)
        CALL PUSHCONTROL(1,1)
      END IF
      IF (a4(2, i, km) .GT. y5) THEN
        CALL PUSHREALARRAY(a4(2, i, km))
        a4(2, i, km) = y5
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHREALARRAY(a4(2, i, km))
        a4(2, i, km) = a4(2, i, km)
        CALL PUSHCONTROL(1,1)
      END IF
      CALL PUSHREALARRAY(dc(i, km))
      dc(i, km) = 0.5*(a4(1, i, km)-a4(2, i, km))
    END DO
! Enforce constraint on the "slope" at the surface
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(2, i, km)) THEN
          CALL PUSHREALARRAY(a4(2, i, km))
          a4(2, i, km) = a4(2, i, km)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(a4(2, i, km))
          a4(2, i, km) = 0.
          CALL PUSHCONTROL(1,1)
        END IF
        IF (0. .LT. a4(3, i, km)) THEN
          CALL PUSHREALARRAY(a4(3, i, km))
          a4(3, i, km) = a4(3, i, km)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(a4(3, i, km))
          a4(3, i, km) = 0.
          CALL PUSHCONTROL(1,1)
        END IF
      END DO
      CALL PUSHCONTROL(2,2)
    ELSE IF (iv .LT. 0) THEN
      DO i=i1,i2
        IF (a4(1, i, km)*a4(3, i, km) .LE. 0.) THEN
          CALL PUSHREALARRAY(a4(3, i, km))
          a4(3, i, km) = 0.
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHCONTROL(2,1)
    ELSE
      CALL PUSHCONTROL(2,0)
    END IF
    DO k=1,km1
      DO i=i1,i2
        CALL PUSHREALARRAY(a4(3, i, k))
        a4(3, i, k) = a4(2, i, k+1)
      END DO
    END DO
!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
    DO k=1,2
      DO i=i1,i2
        CALL PUSHREALARRAY(a4(4, i, k))
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
      CALL PPM_LIMITERS_FWD(dc(i1, k), a4(1, i1, k), it, 0)
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
              CALL PUSHCONTROL(2,0)
            ELSE
              y8 = qmp
              CALL PUSHCONTROL(2,1)
            END IF
          ELSE IF (a4(1, i, k) .GT. lac) THEN
            y8 = lac
            CALL PUSHCONTROL(2,2)
          ELSE
            y8 = a4(1, i, k)
            CALL PUSHCONTROL(2,3)
          END IF
          IF (a4(3, i, k) .LT. y8) THEN
            x2 = y8
            CALL PUSHCONTROL(1,0)
          ELSE
            x2 = a4(3, i, k)
            CALL PUSHCONTROL(1,1)
          END IF
          IF (a4(1, i, k) .LT. qmp) THEN
            IF (qmp .LT. lac) THEN
              y6 = lac
              CALL PUSHCONTROL(2,0)
            ELSE
              y6 = qmp
              CALL PUSHCONTROL(2,1)
            END IF
          ELSE IF (a4(1, i, k) .LT. lac) THEN
            y6 = lac
            CALL PUSHCONTROL(2,2)
          ELSE
            y6 = a4(1, i, k)
            CALL PUSHCONTROL(2,3)
          END IF
          IF (x2 .GT. y6) THEN
            CALL PUSHREALARRAY(a4(3, i, k))
            a4(3, i, k) = y6
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHREALARRAY(a4(3, i, k))
            a4(3, i, k) = x2
            CALL PUSHCONTROL(1,1)
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
              CALL PUSHCONTROL(2,0)
            ELSE
              y9 = qmp
              CALL PUSHCONTROL(2,1)
            END IF
          ELSE IF (a4(1, i, k) .GT. lac) THEN
            y9 = lac
            CALL PUSHCONTROL(2,2)
          ELSE
            y9 = a4(1, i, k)
            CALL PUSHCONTROL(2,3)
          END IF
          IF (a4(2, i, k) .LT. y9) THEN
            x3 = y9
            CALL PUSHCONTROL(1,0)
          ELSE
            x3 = a4(2, i, k)
            CALL PUSHCONTROL(1,1)
          END IF
          IF (a4(1, i, k) .LT. qmp) THEN
            IF (qmp .LT. lac) THEN
              y7 = lac
              CALL PUSHCONTROL(2,0)
            ELSE
              y7 = qmp
              CALL PUSHCONTROL(2,1)
            END IF
          ELSE IF (a4(1, i, k) .LT. lac) THEN
            y7 = lac
            CALL PUSHCONTROL(2,2)
          ELSE
            y7 = a4(1, i, k)
            CALL PUSHCONTROL(2,3)
          END IF
          IF (x3 .GT. y7) THEN
            CALL PUSHREALARRAY(a4(2, i, k))
            a4(2, i, k) = y7
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHREALARRAY(a4(2, i, k))
            a4(2, i, k) = x3
            CALL PUSHCONTROL(1,1)
          END IF
!-------------
! Recompute A6
!-------------
          CALL PUSHREALARRAY(a4(4, i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
! Additional constraint to ensure positivity when kord=7
        IF (iv .EQ. 0 .AND. kord .GE. 6) THEN
          CALL PPM_LIMITERS_FWD(dc(i1, k), a4(1, i1, k), it, 2)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHCONTROL(1,1)
    ELSE
      lmt = kord - 3
      IF (0 .LT. lmt) THEN
        lmt = lmt
      ELSE
        lmt = 0
      END IF
      IF (iv .EQ. 0) THEN
        IF (2 .GT. lmt) THEN
          CALL PUSHCONTROL(1,1)
          lmt = lmt
        ELSE
          CALL PUSHCONTROL(1,1)
          lmt = 2
        END IF
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
      DO k=3,km-2
        IF (kord .NE. 4) THEN
          DO i=i1,i2
            CALL PUSHREALARRAY(a4(4, i, k))
            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (kord .NE. 6) THEN
          CALL PPM_LIMITERS_FWD(dc(i1, k), a4(1, i1, k), it, lmt)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHCONTROL(1,0)
    END IF
    DO k=km1,km
      DO i=i1,i2
        CALL PUSHREALARRAY(a4(4, i, k))
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
      CALL PPM_LIMITERS_FWD(dc(i1, k), a4(1, i1, k), it, 0)
    END DO
    CALL PUSHREALARRAY(dc, (i2-i1+1)*km)
    CALL PUSHREALARRAY(min1)
    CALL PUSHINTEGER(it)
    CALL PUSHINTEGER(km1)
    CALL PUSHREALARRAY(d4, (i2-i1+1)*km)
    CALL PUSHREALARRAY(fac)
    CALL PUSHREALARRAY(delq, (i2-i1+1)*km)
    CALL PUSHREALARRAY(c1)
    CALL PUSHREALARRAY(dq)
    CALL PUSHREALARRAY(df2, (i2-i1+1)*km)
    CALL PUSHREALARRAY(qm)
  END SUBROUTINE PPM_PROFILE_FWD
!  Differentiation of ppm_profile in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
!_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_cor
!e_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod
!.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Raylei
!gh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_o
!rd4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.
!remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d
! fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiter
!s fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv
!_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subg
!rid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_util
!s_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_
!mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_m
!od.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.
!d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_
!v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_cor
!e_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_uti
!ls_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: delp a4
!   with respect to varying inputs: delp a4
  SUBROUTINE PPM_PROFILE_BWD(a4, a4_ad, delp, delp_ad, km, i1, i2, iv, &
&   kord)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: i1
    INTEGER, INTENT(IN) :: i2
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: kord
    REAL, INTENT(IN) :: delp(i1:i2, km)
    REAL :: delp_ad(i1:i2, km)
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_ad(4, i1:i2, km)
    REAL :: dc(i1:i2, km)
    REAL :: dc_ad(i1:i2, km)
    REAL :: h2(i1:i2, km)
    REAL :: h2_ad(i1:i2, km)
    REAL :: delq(i1:i2, km)
    REAL :: delq_ad(i1:i2, km)
    REAL :: df2(i1:i2, km)
    REAL :: df2_ad(i1:i2, km)
    REAL :: d4(i1:i2, km)
    REAL :: d4_ad(i1:i2, km)
    INTEGER :: i, k, km1, lmt, it
    REAL :: fac
    REAL :: a1, a2, c1, c2, c3, d1, d2
    REAL :: a1_ad, a2_ad, c1_ad, c2_ad, c3_ad, d1_ad, d2_ad
    REAL :: qm, dq, lac, qmp, pmp
    REAL :: qm_ad, dq_ad, lac_ad, qmp_ad, pmp_ad
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min1_ad
    INTEGER :: abs0
    REAL :: max1
    REAL :: max1_ad
    REAL :: min2
    REAL :: min2_ad
    REAL :: temp
    REAL :: temp0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: x1_ad
    REAL :: y1_ad
    REAL :: z1_ad
    REAL :: temp1
    REAL :: temp2
    REAL :: temp3
    REAL :: temp4
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp8
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: temp_ad16
    REAL :: temp_ad17
    REAL :: temp_ad18
    REAL :: temp_ad19
    REAL :: temp_ad20
    REAL :: temp_ad21
    REAL :: y2_ad
    REAL :: y3_ad
    REAL :: temp9
    REAL :: temp_ad22
    REAL :: temp_ad23
    REAL :: temp_ad24
    REAL :: temp_ad25
    REAL :: temp_ad26
    REAL :: temp_ad27
    REAL :: temp_ad28
    REAL :: temp_ad29
    REAL :: temp_ad30
    REAL :: temp_ad31
    REAL :: temp_ad32
    REAL :: y4_ad
    REAL :: y5_ad
    REAL :: temp_ad33
    REAL :: temp10
    REAL :: temp11
    REAL :: temp12
    REAL :: temp13
    REAL :: temp14
    REAL :: temp15
    REAL :: temp16
    REAL :: temp_ad34
    REAL :: temp_ad35
    REAL :: temp_ad36
    REAL :: y8_ad
    REAL :: x2_ad
    REAL :: y6_ad
    REAL :: y9_ad
    REAL :: x3_ad
    REAL :: y7_ad
    REAL :: temp_ad37
    REAL :: temp_ad38
    REAL :: temp_ad39
    INTEGER :: branch
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
    CALL POPREALARRAY(qm)
    CALL POPREALARRAY(df2, (i2-i1+1)*km)
    CALL POPREALARRAY(dq)
    CALL POPREALARRAY(c1)
    CALL POPREALARRAY(delq, (i2-i1+1)*km)
    CALL POPREALARRAY(fac)
    CALL POPREALARRAY(d4, (i2-i1+1)*km)
    CALL POPINTEGER(km1)
    CALL POPINTEGER(it)
    CALL POPREALARRAY(min1)
    CALL POPREALARRAY(dc, (i2-i1+1)*km)
    dc_ad = 0.0
    DO k=km,km1,-1
      CALL PPM_LIMITERS_BWD(dc(i1, k), dc_ad(i1, k), a4(1, i1, k), a4_ad&
&                     (1, i1, k), it, 0)
      DO i=i2,i1,-1
        CALL POPREALARRAY(a4(4, i, k))
        temp_ad39 = 3.*a4_ad(4, i, k)
        a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad39
        a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad39
        a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad39
        a4_ad(4, i, k) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO k=km-2,3,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) CALL PPM_LIMITERS_BWD(dc(i1, k), dc_ad(i1, k)&
&                                          , a4(1, i1, k), a4_ad(1, i1, &
&                                          k), it, lmt)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO i=i2,i1,-1
            CALL POPREALARRAY(a4(4, i, k))
            temp_ad38 = 3.*a4_ad(4, i, k)
            a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad38
            a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad38
            a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad38
            a4_ad(4, i, k) = 0.0
          END DO
        END IF
      END DO
      CALL POPCONTROL(1,branch)
    ELSE
      h2_ad = 0.0
      DO k=km-2,3,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) CALL PPM_LIMITERS_BWD(dc(i1, k), dc_ad(i1, k)&
&                                          , a4(1, i1, k), a4_ad(1, i1, &
&                                          k), it, 2)
        DO i=i2,i1,-1
          CALL POPREALARRAY(a4(4, i, k))
          temp_ad37 = 3.*a4_ad(4, i, k)
          a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad37
          a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad37
          a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad37
          a4_ad(4, i, k) = 0.0
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(a4(2, i, k))
            y7_ad = a4_ad(2, i, k)
            a4_ad(2, i, k) = 0.0
            x3_ad = 0.0
          ELSE
            CALL POPREALARRAY(a4(2, i, k))
            x3_ad = a4_ad(2, i, k)
            a4_ad(2, i, k) = 0.0
            y7_ad = 0.0
          END IF
          CALL POPCONTROL(2,branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              lac_ad = y7_ad
              qmp_ad = 0.0
            ELSE
              qmp_ad = y7_ad
              lac_ad = 0.0
            END IF
          ELSE
            IF (branch .EQ. 2) THEN
              lac_ad = y7_ad
            ELSE
              a4_ad(1, i, k) = a4_ad(1, i, k) + y7_ad
              lac_ad = 0.0
            END IF
            qmp_ad = 0.0
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            y9_ad = x3_ad
          ELSE
            a4_ad(2, i, k) = a4_ad(2, i, k) + x3_ad
            y9_ad = 0.0
          END IF
          CALL POPCONTROL(2,branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              lac_ad = lac_ad + y9_ad
            ELSE
              qmp_ad = qmp_ad + y9_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            lac_ad = lac_ad + y9_ad
          ELSE
            a4_ad(1, i, k) = a4_ad(1, i, k) + y9_ad
          END IF
          a4_ad(1, i, k) = a4_ad(1, i, k) + qmp_ad + lac_ad
          h2_ad(i, k+1) = h2_ad(i, k+1) + fac*lac_ad
          dc_ad(i, k) = dc_ad(i, k) - lac_ad
          pmp_ad = -qmp_ad
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(a4(3, i, k))
            y6_ad = a4_ad(3, i, k)
            a4_ad(3, i, k) = 0.0
            x2_ad = 0.0
          ELSE
            CALL POPREALARRAY(a4(3, i, k))
            x2_ad = a4_ad(3, i, k)
            a4_ad(3, i, k) = 0.0
            y6_ad = 0.0
          END IF
          CALL POPCONTROL(2,branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              lac_ad = y6_ad
              qmp_ad = 0.0
            ELSE
              qmp_ad = y6_ad
              lac_ad = 0.0
            END IF
          ELSE
            IF (branch .EQ. 2) THEN
              lac_ad = y6_ad
            ELSE
              a4_ad(1, i, k) = a4_ad(1, i, k) + y6_ad
              lac_ad = 0.0
            END IF
            qmp_ad = 0.0
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            y8_ad = x2_ad
          ELSE
            a4_ad(3, i, k) = a4_ad(3, i, k) + x2_ad
            y8_ad = 0.0
          END IF
          CALL POPCONTROL(2,branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              lac_ad = lac_ad + y8_ad
            ELSE
              qmp_ad = qmp_ad + y8_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            lac_ad = lac_ad + y8_ad
          ELSE
            a4_ad(1, i, k) = a4_ad(1, i, k) + y8_ad
          END IF
          pmp_ad = pmp_ad + qmp_ad
          a4_ad(1, i, k) = a4_ad(1, i, k) + qmp_ad + lac_ad
          h2_ad(i, k-1) = h2_ad(i, k-1) + fac*lac_ad
          dc_ad(i, k) = dc_ad(i, k) + 2.*pmp_ad + lac_ad
        END DO
      END DO
      DO k=km1,2,-1
        DO i=i2,i1,-1
          temp11 = delp(i, k) + 0.5*(delp(i, k-1)+delp(i, k+1))
          temp16 = delp(i, k)**2
          temp10 = temp16/temp11
          temp15 = delp(i, k-1)
          temp14 = dc(i, k-1)/temp15
          temp13 = delp(i, k+1)
          temp12 = dc(i, k+1)/temp13
          temp_ad34 = 2.*temp10*h2_ad(i, k)
          temp_ad35 = (temp12-temp14)*2.*h2_ad(i, k)/temp11
          temp_ad36 = -(temp10*temp_ad35)
          dc_ad(i, k+1) = dc_ad(i, k+1) + temp_ad34/temp13
          delp_ad(i, k+1) = delp_ad(i, k+1) + 0.5*temp_ad36 - temp12*&
&           temp_ad34/temp13
          dc_ad(i, k-1) = dc_ad(i, k-1) - temp_ad34/temp15
          delp_ad(i, k-1) = delp_ad(i, k-1) + 0.5*temp_ad36 + temp14*&
&           temp_ad34/temp15
          delp_ad(i, k) = delp_ad(i, k) + temp_ad36 + 2*delp(i, k)*&
&           temp_ad35
          h2_ad(i, k) = 0.0
        END DO
      END DO
    END IF
    DO k=2,1,-1
      CALL PPM_LIMITERS_BWD(dc(i1, k), dc_ad(i1, k), a4(1, i1, k), a4_ad&
&                     (1, i1, k), it, 0)
      DO i=i2,i1,-1
        CALL POPREALARRAY(a4(4, i, k))
        temp_ad33 = 3.*a4_ad(4, i, k)
        a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad33
        a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad33
        a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad33
        a4_ad(4, i, k) = 0.0
      END DO
    END DO
    DO k=km1,1,-1
      DO i=i2,i1,-1
        CALL POPREALARRAY(a4(3, i, k))
        a4_ad(2, i, k+1) = a4_ad(2, i, k+1) + a4_ad(3, i, k)
        a4_ad(3, i, k) = 0.0
      END DO
    END DO
    CALL POPCONTROL(2,branch)
    IF (branch .NE. 0) THEN
      IF (branch .EQ. 1) THEN
        DO i=i2,i1,-1
          CALL POPCONTROL(1,branch)
          IF (branch .NE. 0) THEN
            CALL POPREALARRAY(a4(3, i, km))
            a4_ad(3, i, km) = 0.0
          END IF
        END DO
      ELSE
        DO i=i2,i1,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(a4(3, i, km))
          ELSE
            CALL POPREALARRAY(a4(3, i, km))
            a4_ad(3, i, km) = 0.0
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(a4(2, i, km))
          ELSE
            CALL POPREALARRAY(a4(2, i, km))
            a4_ad(2, i, km) = 0.0
          END IF
        END DO
      END IF
    END IF
    DO i=i2,i1,-1
      CALL POPREALARRAY(dc(i, km))
      a4_ad(1, i, km) = a4_ad(1, i, km) + 0.5*dc_ad(i, km)
      a4_ad(2, i, km) = a4_ad(2, i, km) - 0.5*dc_ad(i, km)
      dc_ad(i, km) = 0.0
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(a4(2, i, km))
        y5_ad = a4_ad(2, i, km)
        a4_ad(2, i, km) = 0.0
      ELSE
        CALL POPREALARRAY(a4(2, i, km))
        y5_ad = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        a4_ad(1, i, km1) = a4_ad(1, i, km1) + y5_ad
      ELSE
        a4_ad(1, i, km) = a4_ad(1, i, km) + y5_ad
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(a4(2, i, km))
        y4_ad = a4_ad(2, i, km)
        a4_ad(2, i, km) = 0.0
      ELSE
        CALL POPREALARRAY(a4(2, i, km))
        y4_ad = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        a4_ad(1, i, km1) = a4_ad(1, i, km1) + y4_ad
      ELSE
        a4_ad(1, i, km) = a4_ad(1, i, km) + y4_ad
      END IF
      d1 = delp(i, km)
      d2 = delp(i, km1)
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      CALL POPREALARRAY(a4(3, i, km))
      temp_ad22 = d1*a4_ad(3, i, km)
      d1_ad = c1*8.*2*d1*temp_ad22 + (8.*(c1*d1**2)-c3)*a4_ad(3, i, km)
      c3_ad = -temp_ad22
      a4_ad(2, i, km) = a4_ad(2, i, km) + a4_ad(3, i, km)
      a4_ad(3, i, km) = 0.0
      CALL POPREALARRAY(a4(2, i, km))
      temp_ad23 = -((d2+3.*d1)*a4_ad(2, i, km))
      c1_ad = d2*d1*temp_ad23 - 2.0*(d2*(5.*d1+d2)-3.*d1**2)*c3_ad + 8.*&
&       d1**2*temp_ad22
      temp_ad24 = -(c1*d1*d2*a4_ad(2, i, km))
      temp_ad26 = -(2.0*c1*c3_ad)
      temp9 = 2.*d2**2 + d1*(d2+3.*d1)
      temp_ad25 = c1_ad/(d2*temp9)
      qm_ad = a4_ad(2, i, km) - temp_ad25
      a4_ad(2, i, km) = 0.0
      dq_ad = c3_ad - d2*temp_ad25
      temp_ad31 = -((a4(2, i, km1)-qm-d2*dq)*temp_ad25/(d2*temp9))
      temp_ad30 = d2*temp_ad31
      a4_ad(2, i, km1) = a4_ad(2, i, km1) + temp_ad25
      temp_ad32 = 2.*dq_ad/(d1+d2)
      temp_ad27 = -((a4(1, i, km1)-a4(1, i, km))*temp_ad32/(d1+d2))
      a4_ad(1, i, km1) = a4_ad(1, i, km1) + temp_ad32
      temp_ad29 = qm_ad/(d1+d2)
      a4_ad(1, i, km) = a4_ad(1, i, km) + d2*temp_ad29 - temp_ad32
      temp_ad28 = -((d2*a4(1, i, km)+d1*a4(1, i, km1))*temp_ad29/(d1+d2)&
&       )
      d1_ad = d1_ad + (d2*5.-3.*2*d1)*temp_ad26 + temp_ad27 + temp_ad28 &
&       + a4(1, i, km1)*temp_ad29 + (d1*3.+d2+3.*d1)*temp_ad30 + 3.*&
&       temp_ad24 + d2*c1*temp_ad23
      d2_ad = (2*d2+5.*d1)*temp_ad26 + temp_ad27 + temp_ad28 + a4(1, i, &
&       km)*temp_ad29 + (d1+2.*2*d2)*temp_ad30 + temp9*temp_ad31 - dq*&
&       temp_ad25 + temp_ad24 + c1*d1*temp_ad23
      CALL POPREALARRAY(c1)
      CALL POPREALARRAY(dq)
      CALL POPREALARRAY(qm)
      a4_ad(1, i, km1) = a4_ad(1, i, km1) + d1*temp_ad29
      delp_ad(i, km1) = delp_ad(i, km1) + d2_ad
      delp_ad(i, km) = delp_ad(i, km) + d1_ad
    END DO
    CALL POPCONTROL(2,branch)
    IF (branch .LT. 2) THEN
      IF (branch .NE. 0) THEN
        DO i=i2,i1,-1
          CALL POPREALARRAY(a4(3, i, 1))
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + a4_ad(3, i, 1)
          a4_ad(3, i, 1) = 0.0
          CALL POPREALARRAY(a4(2, i, 1))
          a4_ad(1, i, 1) = a4_ad(1, i, 1) + a4_ad(2, i, 1)
          a4_ad(2, i, 1) = 0.0
        END DO
      END IF
    ELSE IF (branch .EQ. 2) THEN
      DO i=i2,i1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          CALL POPREALARRAY(a4(2, i, 1))
          a4_ad(2, i, 1) = 0.0
        END IF
      END DO
    ELSE
      DO i=i2,i1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(a4(2, i, 2))
        ELSE
          CALL POPREALARRAY(a4(2, i, 2))
          a4_ad(2, i, 2) = 0.0
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(a4(2, i, 1))
        ELSE
          CALL POPREALARRAY(a4(2, i, 1))
          a4_ad(2, i, 1) = 0.0
        END IF
      END DO
    END IF
    DO i=i2,i1,-1
      CALL POPREALARRAY(dc(i, 1))
      a4_ad(2, i, 2) = a4_ad(2, i, 2) + 0.5*dc_ad(i, 1)
      a4_ad(1, i, 1) = a4_ad(1, i, 1) - 0.5*dc_ad(i, 1)
      dc_ad(i, 1) = 0.0
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(a4(2, i, 2))
        y3_ad = a4_ad(2, i, 2)
        a4_ad(2, i, 2) = 0.0
      ELSE
        CALL POPREALARRAY(a4(2, i, 2))
        y3_ad = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        a4_ad(1, i, 2) = a4_ad(1, i, 2) + y3_ad
      ELSE
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + y3_ad
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(a4(2, i, 2))
        y2_ad = a4_ad(2, i, 2)
        a4_ad(2, i, 2) = 0.0
      ELSE
        CALL POPREALARRAY(a4(2, i, 2))
        y2_ad = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        a4_ad(1, i, 2) = a4_ad(1, i, 2) + y2_ad
      ELSE
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + y2_ad
      END IF
      d1 = delp(i, 1)
      d2 = delp(i, 2)
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      CALL POPREALARRAY(a4(2, i, 1))
      temp_ad11 = d1*a4_ad(2, i, 1)
      d1_ad = c1*2.*2*d1*temp_ad11 + (2.*(c1*d1**2)-c3)*a4_ad(2, i, 1)
      c3_ad = -temp_ad11
      a4_ad(2, i, 2) = a4_ad(2, i, 2) + a4_ad(2, i, 1)
      a4_ad(2, i, 1) = 0.0
      CALL POPREALARRAY(a4(2, i, 2))
      temp_ad12 = -(0.25*(d2+3.*d1)*a4_ad(2, i, 2))
      c1_ad = d2*d1*temp_ad12 - 0.5*(d2*(5.*d1+d2)-3.*d1**2)*c3_ad + 2.*&
&       d1**2*temp_ad11
      temp_ad13 = -(0.25*c1*d1*d2*a4_ad(2, i, 2))
      temp_ad15 = -(0.5*c1*c3_ad)
      temp8 = 2.*d2**2 + d1*(d2+3.*d1)
      temp_ad14 = 4.*c1_ad/(d2*temp8)
      qm_ad = a4_ad(2, i, 2) - temp_ad14
      a4_ad(2, i, 2) = 0.0
      dq_ad = c3_ad - d2*temp_ad14
      temp_ad20 = -((a4(2, i, 3)-qm-d2*dq)*temp_ad14/(d2*temp8))
      temp_ad19 = d2*temp_ad20
      a4_ad(2, i, 3) = a4_ad(2, i, 3) + temp_ad14
      temp_ad21 = 2.*dq_ad/(d1+d2)
      temp_ad16 = -((a4(1, i, 2)-a4(1, i, 1))*temp_ad21/(d1+d2))
      a4_ad(1, i, 2) = a4_ad(1, i, 2) + temp_ad21
      temp_ad18 = qm_ad/(d1+d2)
      a4_ad(1, i, 1) = a4_ad(1, i, 1) + d2*temp_ad18 - temp_ad21
      temp_ad17 = -((d2*a4(1, i, 1)+d1*a4(1, i, 2))*temp_ad18/(d1+d2))
      d1_ad = d1_ad + (d2*5.-3.*2*d1)*temp_ad15 + temp_ad16 + temp_ad17 &
&       + a4(1, i, 2)*temp_ad18 + (d1*3.+d2+3.*d1)*temp_ad19 + 3.*&
&       temp_ad13 + d2*c1*temp_ad12
      d2_ad = (2*d2+5.*d1)*temp_ad15 + temp_ad16 + temp_ad17 + a4(1, i, &
&       1)*temp_ad18 + (d1+2.*2*d2)*temp_ad19 + temp8*temp_ad20 - dq*&
&       temp_ad14 + temp_ad13 + c1*d1*temp_ad12
      CALL POPREALARRAY(c1)
      CALL POPREALARRAY(dq)
      CALL POPREALARRAY(qm)
      a4_ad(1, i, 2) = a4_ad(1, i, 2) + d1*temp_ad18
      delp_ad(i, 2) = delp_ad(i, 2) + d2_ad
      delp_ad(i, 1) = delp_ad(i, 1) + d1_ad
    END DO
    delq_ad = 0.0
    d4_ad = 0.0
    DO k=km1,3,-1
      DO i=i2,i1,-1
        temp2 = d4(i, k)
        temp1 = delq(i, k-1)/temp2
        temp4 = d4(i, k) + delp(i, k)
        c1 = delq(i, k-1)*delp(i, k-1)/d4(i, k)
        a1 = d4(i, k-1)/(d4(i, k)+delp(i, k-1))
        a2 = d4(i, k+1)/(d4(i, k)+delp(i, k))
        temp7 = d4(i, k-1) + d4(i, k+1)
        temp6 = a1*dc(i, k)
        temp5 = c1*(a1-a2) + a2*dc(i, k-1)
        temp_ad4 = 2.*a4_ad(2, i, k)/temp7
        temp_ad5 = delp(i, k)*temp_ad4
        temp_ad6 = -(delp(i, k-1)*temp_ad4)
        temp_ad7 = -((delp(i, k)*temp5-delp(i, k-1)*temp6)*temp_ad4/&
&         temp7)
        a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + a4_ad(2, i, k)
        c1_ad = (a1-a2)*temp_ad5 + a4_ad(2, i, k)
        delp_ad(i, k) = delp_ad(i, k) + temp5*temp_ad4
        a1_ad = dc(i, k)*temp_ad6 + c1*temp_ad5
        a2_ad = (dc(i, k-1)-c1)*temp_ad5
        dc_ad(i, k-1) = dc_ad(i, k-1) + a2*temp_ad5
        delp_ad(i, k-1) = delp_ad(i, k-1) - temp6*temp_ad4
        dc_ad(i, k) = dc_ad(i, k) + a1*temp_ad6
        d4_ad(i, k-1) = d4_ad(i, k-1) + temp_ad7
        d4_ad(i, k+1) = d4_ad(i, k+1) + a2_ad/temp4 + temp_ad7
        a4_ad(2, i, k) = 0.0
        temp_ad8 = -(d4(i, k+1)*a2_ad/temp4**2)
        d4_ad(i, k) = d4_ad(i, k) + temp_ad8
        delp_ad(i, k) = delp_ad(i, k) + temp_ad8
        temp3 = d4(i, k) + delp(i, k-1)
        temp_ad9 = -(d4(i, k-1)*a1_ad/temp3**2)
        d4_ad(i, k-1) = d4_ad(i, k-1) + a1_ad/temp3
        delp_ad(i, k-1) = delp_ad(i, k-1) + temp1*c1_ad + temp_ad9
        temp_ad10 = delp(i, k-1)*c1_ad/temp2
        d4_ad(i, k) = d4_ad(i, k) + temp_ad9 - temp1*temp_ad10
        delq_ad(i, k-1) = delq_ad(i, k-1) + temp_ad10
      END DO
    END DO
    df2_ad = 0.0
    DO k=km1,2,-1
      DO i=i2,i1,-1
        min1_ad = SIGN(1.d0, min1*df2(i, k))*dc_ad(i, k)
        dc_ad(i, k) = 0.0
        CALL POPCONTROL(2,branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(min1)
            z1_ad = min1_ad
            y1_ad = 0.0
          ELSE
            CALL POPREALARRAY(min1)
            y1_ad = min1_ad
            z1_ad = 0.0
          END IF
          x1_ad = 0.0
        ELSE
          IF (branch .EQ. 2) THEN
            CALL POPREALARRAY(min1)
            z1_ad = min1_ad
            x1_ad = 0.0
          ELSE
            CALL POPREALARRAY(min1)
            x1_ad = min1_ad
            z1_ad = 0.0
          END IF
          y1_ad = 0.0
        END IF
        a4_ad(1, i, k) = a4_ad(1, i, k) + z1_ad
        min2_ad = -z1_ad
        CALL POPCONTROL(2,branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            a4_ad(1, i, k+1) = a4_ad(1, i, k+1) + min2_ad
          ELSE
            a4_ad(1, i, k) = a4_ad(1, i, k) + min2_ad
          END IF
        ELSE IF (branch .EQ. 2) THEN
          a4_ad(1, i, k+1) = a4_ad(1, i, k+1) + min2_ad
        ELSE
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + min2_ad
        END IF
        max1_ad = y1_ad
        a4_ad(1, i, k) = a4_ad(1, i, k) - y1_ad
        CALL POPCONTROL(2,branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            a4_ad(1, i, k+1) = a4_ad(1, i, k+1) + max1_ad
          ELSE
            a4_ad(1, i, k) = a4_ad(1, i, k) + max1_ad
          END IF
        ELSE IF (branch .EQ. 2) THEN
          a4_ad(1, i, k+1) = a4_ad(1, i, k+1) + max1_ad
        ELSE
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + max1_ad
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          df2_ad(i, k) = df2_ad(i, k) + x1_ad
        ELSE
          df2_ad(i, k) = df2_ad(i, k) - x1_ad
        END IF
        c1 = (delp(i, k-1)+0.5*delp(i, k))/d4(i, k+1)
        c2 = (delp(i, k+1)+0.5*delp(i, k))/d4(i, k)
        temp0 = d4(i, k) + delp(i, k+1)
        temp = delp(i, k)/temp0
        temp_ad = temp*df2_ad(i, k)
        temp_ad0 = (c1*delq(i, k)+c2*delq(i, k-1))*df2_ad(i, k)/temp0
        temp_ad1 = -(temp*temp_ad0)
        c1_ad = delq(i, k)*temp_ad
        delq_ad(i, k) = delq_ad(i, k) + c1*temp_ad
        c2_ad = delq(i, k-1)*temp_ad
        delq_ad(i, k-1) = delq_ad(i, k-1) + c2*temp_ad
        delp_ad(i, k) = delp_ad(i, k) + temp_ad0
        df2_ad(i, k) = 0.0
        temp_ad2 = c2_ad/d4(i, k)
        d4_ad(i, k) = d4_ad(i, k) + temp_ad1 - (delp(i, k+1)+0.5*delp(i&
&         , k))*temp_ad2/d4(i, k)
        delp_ad(i, k+1) = delp_ad(i, k+1) + temp_ad2 + temp_ad1
        delp_ad(i, k) = delp_ad(i, k) + 0.5*temp_ad2
        temp_ad3 = c1_ad/d4(i, k+1)
        delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad3
        delp_ad(i, k) = delp_ad(i, k) + 0.5*temp_ad3
        d4_ad(i, k+1) = d4_ad(i, k+1) - (delp(i, k-1)+0.5*delp(i, k))*&
&         temp_ad3/d4(i, k+1)
      END DO
    END DO
    DO k=km,2,-1
      DO i=i2,i1,-1
        delp_ad(i, k-1) = delp_ad(i, k-1) + d4_ad(i, k)
        delp_ad(i, k) = delp_ad(i, k) + d4_ad(i, k)
        d4_ad(i, k) = 0.0
        a4_ad(1, i, k) = a4_ad(1, i, k) + delq_ad(i, k-1)
        a4_ad(1, i, k-1) = a4_ad(1, i, k-1) - delq_ad(i, k-1)
        delq_ad(i, k-1) = 0.0
      END DO
    END DO
  END SUBROUTINE PPM_PROFILE_BWD
!  Differentiation of ppm_limiters in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
!_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_cor
!e_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod
!.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Raylei
!gh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_o
!rd4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.
!remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d
! fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiter
!s fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv
!_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subg
!rid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_util
!s_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_
!mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_m
!od.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.
!d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_
!v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_cor
!e_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_uti
!ls_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: dm a4
!   with respect to varying inputs: dm a4
  SUBROUTINE PPM_LIMITERS_FWD(dm, a4, itot, lmt)
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
      CALL PUSHCONTROL(3,0)
    ELSE IF (lmt .EQ. 0) THEN
! Standard PPM constraint
      DO i=1,itot
        IF (dm(i) .EQ. 0.) THEN
          CALL PUSHREALARRAY(a4(2, i))
          a4(2, i) = a4(1, i)
          CALL PUSHREALARRAY(a4(3, i))
          a4(3, i) = a4(1, i)
          CALL PUSHREALARRAY(a4(4, i))
          a4(4, i) = 0.
          CALL PUSHCONTROL(2,3)
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            CALL PUSHREALARRAY(a4(4, i))
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            CALL PUSHREALARRAY(a4(3, i))
            a4(3, i) = a4(2, i) - a4(4, i)
            CALL PUSHCONTROL(2,2)
          ELSE IF (a6da .GT. da2) THEN
            CALL PUSHREALARRAY(a4(4, i))
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            CALL PUSHREALARRAY(a4(2, i))
            a4(2, i) = a4(3, i) - a4(4, i)
            CALL PUSHCONTROL(2,1)
          ELSE
            CALL PUSHCONTROL(2,0)
          END IF
        END IF
      END DO
      CALL PUSHCONTROL(3,1)
    ELSE IF (lmt .EQ. 1) THEN
! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      DO i=1,itot
        qmp = 2.*dm(i)
        IF (qmp .GE. 0.) THEN
          x1 = qmp
          CALL PUSHCONTROL(1,0)
        ELSE
          x1 = -qmp
          CALL PUSHCONTROL(1,1)
        END IF
        IF (a4(2, i) - a4(1, i) .GE. 0.) THEN
          y1 = a4(2, i) - a4(1, i)
          CALL PUSHCONTROL(1,0)
        ELSE
          y1 = -(a4(2, i)-a4(1, i))
          CALL PUSHCONTROL(1,1)
        END IF
        IF (x1 .GT. y1) THEN
          CALL PUSHREALARRAY(min1)
          min1 = y1
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(min1)
          min1 = x1
          CALL PUSHCONTROL(1,1)
        END IF
        CALL PUSHREALARRAY(a4(2, i))
        a4(2, i) = a4(1, i) - SIGN(min1, qmp)
        IF (qmp .GE. 0.) THEN
          x2 = qmp
          CALL PUSHCONTROL(1,0)
        ELSE
          x2 = -qmp
          CALL PUSHCONTROL(1,1)
        END IF
        IF (a4(3, i) - a4(1, i) .GE. 0.) THEN
          y2 = a4(3, i) - a4(1, i)
          CALL PUSHCONTROL(1,0)
        ELSE
          y2 = -(a4(3, i)-a4(1, i))
          CALL PUSHCONTROL(1,1)
        END IF
        IF (x2 .GT. y2) THEN
          CALL PUSHREALARRAY(min2)
          min2 = y2
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(min2)
          min2 = x2
          CALL PUSHCONTROL(1,1)
        END IF
        CALL PUSHREALARRAY(a4(3, i))
        a4(3, i) = a4(1, i) + SIGN(min2, qmp)
        CALL PUSHREALARRAY(a4(4, i))
        a4(4, i) = 3.*(2.*a4(1, i)-(a4(2, i)+a4(3, i)))
      END DO
      CALL PUSHREALARRAY(min2)
      CALL PUSHREALARRAY(min1)
      CALL PUSHCONTROL(3,2)
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
              CALL PUSHREALARRAY(a4(3, i))
              a4(3, i) = a4(1, i)
              CALL PUSHREALARRAY(a4(2, i))
              a4(2, i) = a4(1, i)
              CALL PUSHREALARRAY(a4(4, i))
              a4(4, i) = 0.
              CALL PUSHCONTROL(3,4)
            ELSE IF (a4(3, i) .GT. a4(2, i)) THEN
              CALL PUSHREALARRAY(a4(4, i))
              a4(4, i) = 3.*(a4(2, i)-a4(1, i))
              CALL PUSHREALARRAY(a4(3, i))
              a4(3, i) = a4(2, i) - a4(4, i)
              CALL PUSHCONTROL(3,3)
            ELSE
              CALL PUSHREALARRAY(a4(4, i))
              a4(4, i) = 3.*(a4(3, i)-a4(1, i))
              CALL PUSHREALARRAY(a4(2, i))
              a4(2, i) = a4(3, i) - a4(4, i)
              CALL PUSHCONTROL(3,2)
            END IF
          ELSE
            CALL PUSHCONTROL(3,1)
          END IF
        ELSE
          CALL PUSHCONTROL(3,0)
        END IF
      END DO
      CALL PUSHCONTROL(3,4)
    ELSE
      CALL PUSHCONTROL(3,3)
    END IF
  END SUBROUTINE PPM_LIMITERS_FWD
!  Differentiation of ppm_limiters in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
!e_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_co
!re_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mo
!d.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayle
!igh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_
!ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod
!.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2
!d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limite
!rs fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic f
!v_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_sub
!grid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_uti
!ls_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils
!_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_
!mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod
!.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp
!_v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_co
!re_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_ut
!ils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: dm a4
!   with respect to varying inputs: dm a4
  SUBROUTINE PPM_LIMITERS_BWD(dm, dm_ad, a4, a4_ad, itot, lmt)
    IMPLICIT NONE
    REAL, INTENT(IN) :: dm(*)
    REAL :: dm_ad(*)
    INTEGER, INTENT(IN) :: itot
    INTEGER, INTENT(IN) :: lmt
    REAL, INTENT(INOUT) :: a4(4, *)
    REAL, INTENT(INOUT) :: a4_ad(4, *)
    REAL :: qmp
    REAL :: qmp_ad
    REAL :: da1, da2, a6da
    REAL :: fmin
    INTEGER :: i
    INTRINSIC ABS
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min1_ad
    REAL :: min2
    REAL :: min2_ad
    REAL :: abs0
    REAL :: x1_ad
    REAL :: y1_ad
    REAL :: x2_ad
    REAL :: y2_ad
    REAL :: temp_ad
    INTEGER :: branch
    REAL :: x2
    REAL :: x1
    REAL :: y2
    REAL :: y1
    CALL POPCONTROL(3,branch)
    IF (branch .LT. 2) THEN
      IF (branch .NE. 0) THEN
        DO i=itot,1,-1
          CALL POPCONTROL(2,branch)
          IF (branch .LT. 2) THEN
            IF (branch .NE. 0) THEN
              CALL POPREALARRAY(a4(2, i))
              a4_ad(3, i) = a4_ad(3, i) + a4_ad(2, i)
              a4_ad(4, i) = a4_ad(4, i) - a4_ad(2, i)
              a4_ad(2, i) = 0.0
              CALL POPREALARRAY(a4(4, i))
              a4_ad(3, i) = a4_ad(3, i) + 3.*a4_ad(4, i)
              a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
              a4_ad(4, i) = 0.0
            END IF
          ELSE IF (branch .EQ. 2) THEN
            CALL POPREALARRAY(a4(3, i))
            a4_ad(2, i) = a4_ad(2, i) + a4_ad(3, i)
            a4_ad(4, i) = a4_ad(4, i) - a4_ad(3, i)
            a4_ad(3, i) = 0.0
            CALL POPREALARRAY(a4(4, i))
            a4_ad(2, i) = a4_ad(2, i) + 3.*a4_ad(4, i)
            a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
            a4_ad(4, i) = 0.0
          ELSE
            CALL POPREALARRAY(a4(4, i))
            a4_ad(4, i) = 0.0
            CALL POPREALARRAY(a4(3, i))
            a4_ad(1, i) = a4_ad(1, i) + a4_ad(3, i)
            a4_ad(3, i) = 0.0
            CALL POPREALARRAY(a4(2, i))
            a4_ad(1, i) = a4_ad(1, i) + a4_ad(2, i)
            a4_ad(2, i) = 0.0
          END IF
        END DO
      END IF
    ELSE IF (branch .EQ. 2) THEN
      CALL POPREALARRAY(min1)
      CALL POPREALARRAY(min2)
      DO i=itot,1,-1
        CALL POPREALARRAY(a4(4, i))
        temp_ad = 3.*a4_ad(4, i)
        a4_ad(1, i) = a4_ad(1, i) + 2.*temp_ad
        a4_ad(2, i) = a4_ad(2, i) - temp_ad
        a4_ad(3, i) = a4_ad(3, i) - temp_ad
        a4_ad(4, i) = 0.0
        qmp = 2.*dm(i)
        CALL POPREALARRAY(a4(3, i))
        a4_ad(1, i) = a4_ad(1, i) + a4_ad(3, i)
        min2_ad = SIGN(1.d0, min2*qmp)*a4_ad(3, i)
        a4_ad(3, i) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(min2)
          y2_ad = min2_ad
          x2_ad = 0.0
        ELSE
          CALL POPREALARRAY(min2)
          x2_ad = min2_ad
          y2_ad = 0.0
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          a4_ad(3, i) = a4_ad(3, i) + y2_ad
          a4_ad(1, i) = a4_ad(1, i) - y2_ad
        ELSE
          a4_ad(1, i) = a4_ad(1, i) + y2_ad
          a4_ad(3, i) = a4_ad(3, i) - y2_ad
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          qmp_ad = x2_ad
        ELSE
          qmp_ad = -x2_ad
        END IF
        CALL POPREALARRAY(a4(2, i))
        a4_ad(1, i) = a4_ad(1, i) + a4_ad(2, i)
        min1_ad = -(SIGN(1.d0, min1*qmp)*a4_ad(2, i))
        a4_ad(2, i) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(min1)
          y1_ad = min1_ad
          x1_ad = 0.0
        ELSE
          CALL POPREALARRAY(min1)
          x1_ad = min1_ad
          y1_ad = 0.0
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          a4_ad(2, i) = a4_ad(2, i) + y1_ad
          a4_ad(1, i) = a4_ad(1, i) - y1_ad
        ELSE
          a4_ad(1, i) = a4_ad(1, i) + y1_ad
          a4_ad(2, i) = a4_ad(2, i) - y1_ad
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          qmp_ad = qmp_ad + x1_ad
        ELSE
          qmp_ad = qmp_ad - x1_ad
        END IF
        dm_ad(i) = dm_ad(i) + 2.*qmp_ad
      END DO
    ELSE IF (branch .NE. 3) THEN
      DO i=itot,1,-1
        CALL POPCONTROL(3,branch)
        IF (branch .GE. 2) THEN
          IF (branch .EQ. 2) THEN
            CALL POPREALARRAY(a4(2, i))
            a4_ad(3, i) = a4_ad(3, i) + a4_ad(2, i)
            a4_ad(4, i) = a4_ad(4, i) - a4_ad(2, i)
            a4_ad(2, i) = 0.0
            CALL POPREALARRAY(a4(4, i))
            a4_ad(3, i) = a4_ad(3, i) + 3.*a4_ad(4, i)
            a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
            a4_ad(4, i) = 0.0
          ELSE IF (branch .EQ. 3) THEN
            CALL POPREALARRAY(a4(3, i))
            a4_ad(2, i) = a4_ad(2, i) + a4_ad(3, i)
            a4_ad(4, i) = a4_ad(4, i) - a4_ad(3, i)
            a4_ad(3, i) = 0.0
            CALL POPREALARRAY(a4(4, i))
            a4_ad(2, i) = a4_ad(2, i) + 3.*a4_ad(4, i)
            a4_ad(1, i) = a4_ad(1, i) - 3.*a4_ad(4, i)
            a4_ad(4, i) = 0.0
          ELSE
            CALL POPREALARRAY(a4(4, i))
            a4_ad(4, i) = 0.0
            CALL POPREALARRAY(a4(2, i))
            a4_ad(1, i) = a4_ad(1, i) + a4_ad(2, i)
            a4_ad(2, i) = 0.0
            CALL POPREALARRAY(a4(3, i))
            a4_ad(1, i) = a4_ad(1, i) + a4_ad(3, i)
            a4_ad(3, i) = 0.0
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE PPM_LIMITERS_BWD
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
    INTEGER :: arg1
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
        arg1 = ie + 1
        CALL REMAP_2D(km, pv1, v_r(is:ie+1, j:j, 1:km), kn, pv2, v(is:ie&
&               +1, j:j, 1:kn), is, arg1, -1, kord)
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
!  Differentiation of map1_cubic in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_m
!od.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_
!mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.m
!ix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh
!_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord
!4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.re
!map_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d f
!v_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters 
!fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_r
!estart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgri
!d_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_
!mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mo
!d.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod
!.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2
!a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v_
!fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_
!mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils
!_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 q2
!   with respect to varying inputs: pe1 pe2 q2
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  map1_cubic --- Cubic Interpolation for vertical re-mapping
!
! !INTERFACE:
  SUBROUTINE MAP1_CUBIC_FWD(km, pe1, kn, pe2, q2, i1, i2, j, ibeg, iend&
&   , jbeg, jend, akap, t_var, conserv)
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
    INTEGER :: ad_count

    qx = 0.0
    logpl1 = 0.0
    logpl2 = 0.0
    dlogp1 = 0.0
    vsum1 = 0.0
    vsum2 = 0.0
    am2 = 0.0
    am1 = 0.0
    ap0 = 0.0
    ap1 = 0.0
    p = 0.0
    plp1 = 0.0
    plp0 = 0.0
    plm1 = 0.0
    plm2 = 0.0
    dlp0 = 0.0
    dlm1 = 0.0
    dlm2 = 0.0
    i = 0
    k = 0
    lm2 = 0
    lm1 = 0
    lp0 = 0
    lp1 = 0
    ad_count = 0

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
      CALL PUSHCONTROL(2,1)
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
      CALL PUSHCONTROL(2,2)
    CASE (3) 
! Potential Temperature Remapping in P^KAPPA
      DO k=1,km
        qx(:, k) = q2(i1:i2, j, k)
        logpl1(:, k) = EXP(akap*LOG(r2*(pe1(:, k)+pe1(:, k+1))))
      END DO
      DO k=1,kn
        logpl2(:, k) = EXP(akap*LOG(r2*(pe2(:, k)+pe2(:, k+1))))
      END DO
      DO k=1,km-1
        dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
      END DO
      CALL PUSHCONTROL(2,3)
    CASE DEFAULT
      CALL PUSHCONTROL(2,0)
    END SELECT
    IF (conserv) THEN
! Compute vertical integral of Input TE
! -------------------------------------
      vsum1(:) = r0
      DO i=i1,i2
        DO k=1,km
          vsum1(i) = vsum1(i) + qx(i, k)*(pe1(i, k+1)-pe1(i, k))
        END DO
        CALL PUSHREALARRAY(vsum1(i))
        vsum1(i) = vsum1(i)/(pe1(i, km+1)-pe1(i, 1))
      END DO
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
! Interpolate TE onto target Pressures
! ------------------------------------
    DO i=i1,i2
      DO k=1,kn
        CALL PUSHINTEGER(lp0)
        lp0 = 1
        ad_count = 1
        DO WHILE (lp0 .LE. km)
          IF (logpl1(i, lp0) .LT. logpl2(i, k)) THEN
            lp0 = lp0 + 1
            ad_count = ad_count + 1
          ELSE
            GOTO 100
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
        CALL PUSHINTEGER(ad_count)
        GOTO 110
 100    CALL PUSHCONTROL(1,1)
        CALL PUSHINTEGER(ad_count)
 110    IF (lp0 - 1 .LT. 1) THEN
          CALL PUSHINTEGER(lm1)
          lm1 = 1
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHINTEGER(lm1)
          lm1 = lp0 - 1
          CALL PUSHCONTROL(1,1)
        END IF
        IF (lp0 .GT. km) THEN
          lp0 = km
        ELSE
          lp0 = lp0
        END IF
! Extrapolate Linearly in LogP above first model level
! ----------------------------------------------------
        IF (lm1 .EQ. 1 .AND. lp0 .EQ. 1) THEN
          CALL PUSHREALARRAY(q2(i, j, k))
          q2(i, j, k) = qx(i, 1) + (qx(i, 2)-qx(i, 1))*(logpl2(i, k)-&
&           logpl1(i, 1))/(logpl1(i, 2)-logpl1(i, 1))
! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
          CALL PUSHCONTROL(2,3)
        ELSE IF (lm1 .EQ. km .AND. lp0 .EQ. km) THEN
          CALL PUSHREALARRAY(q2(i, j, k))
          q2(i, j, k) = qx(i, km) + (qx(i, km)-qx(i, km-1))*(logpl2(i, k&
&           )-logpl1(i, km))/(logpl1(i, km)-logpl1(i, km-1))
! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
          CALL PUSHCONTROL(2,2)
        ELSE IF (lm1 .EQ. 1 .OR. lp0 .EQ. km) THEN
          CALL PUSHREALARRAY(q2(i, j, k))
          q2(i, j, k) = qx(i, lp0) + (qx(i, lm1)-qx(i, lp0))*(logpl2(i, &
&           k)-logpl1(i, lp0))/(logpl1(i, lm1)-logpl1(i, lp0))
! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
          CALL PUSHCONTROL(2,1)
        ELSE
          CALL PUSHINTEGER(lp1)
          lp1 = lp0 + 1
          CALL PUSHINTEGER(lm2)
          lm2 = lm1 - 1
          CALL PUSHREALARRAY(p)
          p = logpl2(i, k)
          plp1 = logpl1(i, lp1)
          plp0 = logpl1(i, lp0)
          plm1 = logpl1(i, lm1)
          plm2 = logpl1(i, lm2)
          CALL PUSHREALARRAY(dlp0)
          dlp0 = dlogp1(i, lp0)
          CALL PUSHREALARRAY(dlm1)
          dlm1 = dlogp1(i, lm1)
          CALL PUSHREALARRAY(dlm2)
          dlm2 = dlogp1(i, lm2)
          CALL PUSHREALARRAY(ap1)
          ap1 = (p-plp0)*(p-plm1)*(p-plm2)/(dlp0*(dlp0+dlm1)*(dlp0+dlm1+&
&           dlm2))
          CALL PUSHREALARRAY(ap0)
          ap0 = (plp1-p)*(p-plm1)*(p-plm2)/(dlp0*dlm1*(dlm1+dlm2))
          CALL PUSHREALARRAY(am1)
          am1 = (plp1-p)*(plp0-p)*(p-plm2)/(dlm1*dlm2*(dlp0+dlm1))
          CALL PUSHREALARRAY(am2)
          am2 = (plp1-p)*(plp0-p)*(plm1-p)/(dlm2*(dlm1+dlm2)*(dlp0+dlm1+&
&           dlm2))
          CALL PUSHREALARRAY(q2(i, j, k))
          q2(i, j, k) = ap1*qx(i, lp1) + ap0*qx(i, lp0) + am1*qx(i, lm1)&
&           + am2*qx(i, lm2)
          CALL PUSHCONTROL(2,0)
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
        CALL PUSHREALARRAY(vsum2(i))
        vsum2(i) = vsum2(i)/(pe2(i, kn+1)-pe2(i, 1))
      END DO
! Adjust Final TE to conserve
! ---------------------------
      DO i=i1,i2
        DO k=1,kn
          CALL PUSHREALARRAY(q2(i, j, k))
          q2(i, j, k) = q2(i, j, k) + vsum1(i) - vsum2(i)
        END DO
      END DO
      CALL PUSHREALARRAY(logpl2, (i2-i1+1)*kn)
      CALL PUSHREALARRAY(am2)
      CALL PUSHREALARRAY(logpl1, (i2-i1+1)*km)
      CALL PUSHREALARRAY(am1)
      CALL PUSHREALARRAY(ap1)
      CALL PUSHREALARRAY(ap0)
      CALL PUSHREALARRAY(qx, (i2-i1+1)*km)
      CALL PUSHREALARRAY(dlm2)
      CALL PUSHINTEGER(lm2)
      CALL PUSHREALARRAY(dlm1)
      CALL PUSHINTEGER(lm1)
      CALL PUSHREALARRAY(p)
      CALL PUSHINTEGER(lp1)
      CALL PUSHREALARRAY(dlp0)
      CALL PUSHINTEGER(lp0)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHREALARRAY(logpl2, (i2-i1+1)*kn)
      CALL PUSHREALARRAY(am2)
      CALL PUSHREALARRAY(logpl1, (i2-i1+1)*km)
      CALL PUSHREALARRAY(am1)
      CALL PUSHREALARRAY(ap1)
      CALL PUSHREALARRAY(ap0)
      CALL PUSHREALARRAY(qx, (i2-i1+1)*km)
      CALL PUSHREALARRAY(dlm2)
      CALL PUSHINTEGER(lm2)
      CALL PUSHREALARRAY(dlm1)
      CALL PUSHINTEGER(lm1)
      CALL PUSHREALARRAY(p)
      CALL PUSHINTEGER(lp1)
      CALL PUSHREALARRAY(dlp0)
      CALL PUSHINTEGER(lp0)
      CALL PUSHCONTROL(1,1)
    END IF
  END SUBROUTINE MAP1_CUBIC_FWD
!  Differentiation of map1_cubic in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
!mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core
!_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.
!mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleig
!h_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_or
!d4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.r
!emap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d 
!fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters
! fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_
!restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgr
!id_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils
!_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_m
!od.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mo
!d.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d
!2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
!_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core
!_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_util
!s_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 q2
!   with respect to varying inputs: pe1 pe2 q2
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  map1_cubic --- Cubic Interpolation for vertical re-mapping
!
! !INTERFACE:
  SUBROUTINE MAP1_CUBIC_BWD(km, pe1, pe1_ad, kn, pe2, pe2_ad, q2, q2_ad&
&   , i1, i2, j, ibeg, iend, jbeg, jend, akap, t_var, conserv)
    IMPLICIT NONE
!EOC
    INTEGER, INTENT(IN) :: i1
    INTEGER, INTENT(IN) :: i2
    REAL, INTENT(IN) :: akap
    INTEGER, INTENT(IN) :: t_var
    LOGICAL, INTENT(IN) :: conserv
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: kn
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL :: pe1_ad(i1:i2, km+1)
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL :: pe2_ad(i1:i2, kn+1)
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(INOUT) :: q2_ad(ibeg:iend, jbeg:jend, kn)
    REAL :: qx(i1:i2, km)
    REAL :: qx_ad(i1:i2, km)
    REAL :: logpl1(i1:i2, km)
    REAL :: logpl1_ad(i1:i2, km)
    REAL :: logpl2(i1:i2, kn)
    REAL :: logpl2_ad(i1:i2, kn)
    REAL :: dlogp1(i1:i2, km)
    REAL :: dlogp1_ad(i1:i2, km)
    REAL :: vsum1(i1:i2)
    REAL :: vsum1_ad(i1:i2)
    REAL :: vsum2(i1:i2)
    REAL :: vsum2_ad(i1:i2)
    REAL :: am2, am1, ap0, ap1, p, plp1, plp0, plm1, plm2, dlp0, dlm1, &
&   dlm2
    REAL :: am2_ad, am1_ad, ap0_ad, ap1_ad, p_ad, plp1_ad, plp0_ad, &
&   plm1_ad, plm2_ad, dlp0_ad, dlm1_ad, dlm2_ad
    INTEGER :: i, k, lm2, lm1, lp0, lp1
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    INTRINSIC MIN
    REAL, DIMENSION(i2-i1+1) :: temp_ad
    REAL, DIMENSION(i2-i1+1) :: temp_ad0
    REAL, DIMENSION(i2-i1+1) :: temp_ad1
    REAL, DIMENSION(i2-i1+1) :: temp_ad2
    REAL, DIMENSION(i2-i1+1) :: temp
    REAL, DIMENSION(i2-i1+1) :: temp_ad3
    REAL, DIMENSION(i2-i1+1) :: temp0
    REAL, DIMENSION(i2-i1+1) :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp1
    REAL :: temp_ad6
    REAL :: temp2
    REAL :: temp3
    REAL :: temp4
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp8
    REAL :: temp9
    REAL :: temp10
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp11
    REAL :: temp12
    REAL :: temp13
    REAL :: temp14
    REAL :: temp15
    REAL :: temp16
    REAL :: temp17
    REAL :: temp18
    REAL :: temp19
    REAL :: temp20
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: temp_ad16
    REAL :: temp_ad17
    REAL :: temp_ad18
    REAL :: temp_ad19
    REAL :: temp_ad20
    REAL :: temp_ad21
    REAL :: temp_ad22
    REAL :: temp_ad23
    REAL :: temp_ad24
    REAL :: temp_ad25
    REAL :: temp_ad26
    REAL :: temp_ad27
    REAL :: temp_ad28
    REAL :: temp_ad29
    REAL :: temp21
    REAL :: temp_ad30
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch

    qx = 0.0
    logpl1 = 0.0
    logpl2 = 0.0
    dlogp1 = 0.0
    vsum1 = 0.0
    vsum2 = 0.0
    am2 = 0.0
    am1 = 0.0
    ap0 = 0.0
    ap1 = 0.0
    p = 0.0
    plp1 = 0.0
    plp0 = 0.0
    plm1 = 0.0
    plm2 = 0.0
    dlp0 = 0.0
    dlm1 = 0.0
    dlm2 = 0.0
    i = 0
    k = 0
    lm2 = 0
    lm1 = 0
    lp0 = 0
    lp1 = 0
    ad_count = 0
    branch = 0

    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(lp0)
      CALL POPREALARRAY(dlp0)
      CALL POPINTEGER(lp1)
      CALL POPREALARRAY(p)
      CALL POPINTEGER(lm1)
      CALL POPREALARRAY(dlm1)
      CALL POPINTEGER(lm2)
      CALL POPREALARRAY(dlm2)
      CALL POPREALARRAY(qx, (i2-i1+1)*km)
      CALL POPREALARRAY(ap0)
      CALL POPREALARRAY(ap1)
      CALL POPREALARRAY(am1)
      CALL POPREALARRAY(logpl1, (i2-i1+1)*km)
      CALL POPREALARRAY(am2)
      CALL POPREALARRAY(logpl2, (i2-i1+1)*kn)
      vsum1_ad = 0.0
      vsum2_ad = 0.0
      DO i=i2,i1,-1
        DO k=kn,1,-1
          CALL POPREALARRAY(q2(i, j, k))
          vsum1_ad(i) = vsum1_ad(i) + q2_ad(i, j, k)
          vsum2_ad(i) = vsum2_ad(i) - q2_ad(i, j, k)
        END DO
      END DO
      DO i=i2,i1,-1
        CALL POPREALARRAY(vsum2(i))
        temp21 = pe2(i, kn+1) - pe2(i, 1)
        temp_ad30 = -(vsum2(i)*vsum2_ad(i)/temp21**2)
        pe2_ad(i, kn+1) = pe2_ad(i, kn+1) + temp_ad30
        pe2_ad(i, 1) = pe2_ad(i, 1) - temp_ad30
        vsum2_ad(i) = vsum2_ad(i)/temp21
        DO k=kn,1,-1
          temp_ad29 = q2(i, j, k)*vsum2_ad(i)
          q2_ad(i, j, k) = q2_ad(i, j, k) + (pe2(i, k+1)-pe2(i, k))*&
&           vsum2_ad(i)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad29
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad29
        END DO
      END DO
    ELSE
      CALL POPINTEGER(lp0)
      CALL POPREALARRAY(dlp0)
      CALL POPINTEGER(lp1)
      CALL POPREALARRAY(p)
      CALL POPINTEGER(lm1)
      CALL POPREALARRAY(dlm1)
      CALL POPINTEGER(lm2)
      CALL POPREALARRAY(dlm2)
      CALL POPREALARRAY(qx, (i2-i1+1)*km)
      CALL POPREALARRAY(ap0)
      CALL POPREALARRAY(ap1)
      CALL POPREALARRAY(am1)
      CALL POPREALARRAY(logpl1, (i2-i1+1)*km)
      CALL POPREALARRAY(am2)
      CALL POPREALARRAY(logpl2, (i2-i1+1)*kn)
      vsum1_ad = 0.0
    END IF
    qx_ad = 0.0
    logpl1_ad = 0.0
    logpl2_ad = 0.0
    dlogp1_ad = 0.0
    DO i=i2,i1,-1
      DO k=kn,1,-1
        CALL POPCONTROL(2,branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(q2(i, j, k))
            ap1_ad = qx(i, lp1)*q2_ad(i, j, k)
            qx_ad(i, lp1) = qx_ad(i, lp1) + ap1*q2_ad(i, j, k)
            ap0_ad = qx(i, lp0)*q2_ad(i, j, k)
            qx_ad(i, lp0) = qx_ad(i, lp0) + ap0*q2_ad(i, j, k)
            am1_ad = qx(i, lm1)*q2_ad(i, j, k)
            qx_ad(i, lm1) = qx_ad(i, lm1) + am1*q2_ad(i, j, k)
            am2_ad = qx(i, lm2)*q2_ad(i, j, k)
            qx_ad(i, lm2) = qx_ad(i, lm2) + am2*q2_ad(i, j, k)
            q2_ad(i, j, k) = 0.0
            plp0 = logpl1(i, lp0)
            plp1 = logpl1(i, lp1)
            plm1 = logpl1(i, lm1)
            CALL POPREALARRAY(am2)
            temp20 = dlm2*(dlm1+dlm2)
            temp19 = temp20*(dlp0+dlm1+dlm2)
            temp_ad13 = am2_ad/temp19
            temp_ad14 = (plm1-p)*temp_ad13
            temp18 = (plp1-p)*(plp0-p)
            temp_ad15 = -(temp18*(plm1-p)*temp_ad13/temp19)
            temp_ad16 = (dlp0+dlm1+dlm2)*temp_ad15
            temp_ad17 = temp20*temp_ad15
            plm2 = logpl1(i, lm2)
            CALL POPREALARRAY(am1)
            temp17 = dlm1*dlm2
            temp_ad20 = am1_ad/(temp17*(dlp0+dlm1))
            temp_ad18 = (p-plm2)*temp_ad20
            temp16 = (plp1-p)*(plp0-p)
            temp_ad24 = -(temp16*(p-plm2)*temp_ad20/(temp17*(dlp0+dlm1))&
&             )
            CALL POPREALARRAY(ap0)
            temp15 = dlp0*dlm1
            temp_ad23 = ap0_ad/(temp15*(dlm1+dlm2))
            temp_ad19 = (p-plm2)*temp_ad23
            plp1_ad = (plp0-p)*temp_ad18 + (p-plm1)*temp_ad19 + (plp0-p)&
&             *temp_ad14
            temp14 = (plp1-p)*(p-plm1)
            temp_ad26 = -(temp14*(p-plm2)*temp_ad23/(temp15*(dlm1+dlm2))&
&             )
            CALL POPREALARRAY(ap1)
            temp13 = dlp0*(dlp0+dlm1)
            temp12 = temp13*(dlp0+dlm1+dlm2)
            temp_ad22 = ap1_ad/temp12
            temp_ad21 = (p-plm2)*temp_ad22
            plp0_ad = (plp1-p)*temp_ad18 - (p-plm1)*temp_ad21 + (plp1-p)&
&             *temp_ad14
            plm1_ad = temp18*temp_ad13 - (p-plp0)*temp_ad21 - (plp1-p)*&
&             temp_ad19
            temp11 = (p-plp0)*(p-plm1)
            p_ad = (2*p-plp1-plp0)*temp_ad18 + temp16*temp_ad20 + (2*p-&
&             plp0-plm1)*temp_ad21 + temp11*temp_ad22 + temp14*temp_ad23&
&             + (plp1-2*p+plm1)*temp_ad19 - temp18*temp_ad13 + (2*p-plp1&
&             -plp0)*temp_ad14
            plm2_ad = -(temp14*temp_ad23) - temp11*temp_ad22 - temp16*&
&             temp_ad20
            temp_ad28 = -(temp11*(p-plm2)*temp_ad22/temp12)
            temp_ad27 = (dlp0+dlm1+dlm2)*temp_ad28
            temp_ad25 = temp13*temp_ad28
            dlm2_ad = (dlp0+dlm1)*dlm1*temp_ad24 + temp_ad25 + temp15*&
&             temp_ad26 + temp_ad17 + (2*dlm2+dlm1)*temp_ad16
            dlm1_ad = (temp17+(dlp0+dlm1)*dlm2)*temp_ad24 + dlp0*&
&             temp_ad27 + temp_ad25 + (temp15+(dlm1+dlm2)*dlp0)*&
&             temp_ad26 + temp_ad17 + dlm2*temp_ad16
            dlp0_ad = temp17*temp_ad24 + (2*dlp0+dlm1)*temp_ad27 + &
&             temp_ad25 + (dlm1+dlm2)*dlm1*temp_ad26 + temp_ad17
            CALL POPREALARRAY(dlm2)
            dlogp1_ad(i, lm2) = dlogp1_ad(i, lm2) + dlm2_ad
            CALL POPREALARRAY(dlm1)
            dlogp1_ad(i, lm1) = dlogp1_ad(i, lm1) + dlm1_ad
            CALL POPREALARRAY(dlp0)
            dlogp1_ad(i, lp0) = dlogp1_ad(i, lp0) + dlp0_ad
            logpl1_ad(i, lm2) = logpl1_ad(i, lm2) + plm2_ad
            logpl1_ad(i, lm1) = logpl1_ad(i, lm1) + plm1_ad
            logpl1_ad(i, lp0) = logpl1_ad(i, lp0) + plp0_ad
            logpl1_ad(i, lp1) = logpl1_ad(i, lp1) + plp1_ad
            CALL POPREALARRAY(p)
            logpl2_ad(i, k) = logpl2_ad(i, k) + p_ad
            CALL POPINTEGER(lm2)
            CALL POPINTEGER(lp1)
          ELSE
            CALL POPREALARRAY(q2(i, j, k))
            temp8 = logpl1(i, lm1) - logpl1(i, lp0)
            temp_ad11 = q2_ad(i, j, k)/temp8
            temp10 = logpl2(i, k) - logpl1(i, lp0)
            temp9 = qx(i, lm1) - qx(i, lp0)
            temp_ad12 = -(temp9*temp10*temp_ad11/temp8)
            qx_ad(i, lp0) = qx_ad(i, lp0) + q2_ad(i, j, k) - temp10*&
&             temp_ad11
            qx_ad(i, lm1) = qx_ad(i, lm1) + temp10*temp_ad11
            logpl2_ad(i, k) = logpl2_ad(i, k) + temp9*temp_ad11
            logpl1_ad(i, lp0) = logpl1_ad(i, lp0) - temp_ad12 - temp9*&
&             temp_ad11
            logpl1_ad(i, lm1) = logpl1_ad(i, lm1) + temp_ad12
            q2_ad(i, j, k) = 0.0
          END IF
        ELSE IF (branch .EQ. 2) THEN
          CALL POPREALARRAY(q2(i, j, k))
          temp5 = logpl1(i, km) - logpl1(i, km-1)
          temp_ad9 = q2_ad(i, j, k)/temp5
          temp7 = logpl2(i, k) - logpl1(i, km)
          temp6 = qx(i, km) - qx(i, km-1)
          temp_ad10 = -(temp6*temp7*temp_ad9/temp5)
          qx_ad(i, km) = qx_ad(i, km) + temp7*temp_ad9 + q2_ad(i, j, k)
          qx_ad(i, km-1) = qx_ad(i, km-1) - temp7*temp_ad9
          logpl2_ad(i, k) = logpl2_ad(i, k) + temp6*temp_ad9
          logpl1_ad(i, km) = logpl1_ad(i, km) + temp_ad10 - temp6*&
&           temp_ad9
          logpl1_ad(i, km-1) = logpl1_ad(i, km-1) - temp_ad10
          q2_ad(i, j, k) = 0.0
        ELSE
          CALL POPREALARRAY(q2(i, j, k))
          temp3 = logpl1(i, 2) - logpl1(i, 1)
          temp4 = qx(i, 2) - qx(i, 1)
          temp2 = temp4/temp3
          temp_ad7 = (logpl2(i, k)-logpl1(i, 1))*q2_ad(i, j, k)/temp3
          temp_ad8 = -(temp2*temp_ad7)
          qx_ad(i, 1) = qx_ad(i, 1) + q2_ad(i, j, k) - temp_ad7
          logpl2_ad(i, k) = logpl2_ad(i, k) + temp2*q2_ad(i, j, k)
          logpl1_ad(i, 1) = logpl1_ad(i, 1) - temp_ad8 - temp2*q2_ad(i, &
&           j, k)
          qx_ad(i, 2) = qx_ad(i, 2) + temp_ad7
          logpl1_ad(i, 2) = logpl1_ad(i, 2) + temp_ad8
          q2_ad(i, j, k) = 0.0
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPINTEGER(lm1)
        ELSE
          CALL POPINTEGER(lm1)
        END IF
        CALL POPINTEGER(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL(1,branch)
        END DO
        CALL POPINTEGER(lp0)
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      DO i=i2,i1,-1
        CALL POPREALARRAY(vsum1(i))
        temp1 = pe1(i, km+1) - pe1(i, 1)
        temp_ad6 = -(vsum1(i)*vsum1_ad(i)/temp1**2)
        pe1_ad(i, km+1) = pe1_ad(i, km+1) + temp_ad6
        pe1_ad(i, 1) = pe1_ad(i, 1) - temp_ad6
        vsum1_ad(i) = vsum1_ad(i)/temp1
        DO k=km,1,-1
          temp_ad5 = qx(i, k)*vsum1_ad(i)
          qx_ad(i, k) = qx_ad(i, k) + (pe1(i, k+1)-pe1(i, k))*vsum1_ad(i&
&           )
          pe1_ad(i, k+1) = pe1_ad(i, k+1) + temp_ad5
          pe1_ad(i, k) = pe1_ad(i, k) - temp_ad5
        END DO
      END DO
    END IF
    CALL POPCONTROL(2,branch)
    IF (branch .LT. 2) THEN
      IF (branch .NE. 0) THEN
        DO k=km-1,1,-1
          logpl1_ad(:, k+1) = logpl1_ad(:, k+1) + dlogp1_ad(:, k)
          logpl1_ad(:, k) = logpl1_ad(:, k) - dlogp1_ad(:, k)
          dlogp1_ad(:, k) = 0.0
        END DO
        DO k=kn,1,-1
          temp_ad0 = logpl2_ad(:, k)/(pe2(:, k)+pe2(:, k+1))
          pe2_ad(:, k) = pe2_ad(:, k) + temp_ad0
          pe2_ad(:, k+1) = pe2_ad(:, k+1) + temp_ad0
          logpl2_ad(:, k) = 0.0
        END DO
        DO k=km,1,-1
          temp_ad = logpl1_ad(:, k)/(pe1(:, k)+pe1(:, k+1))
          pe1_ad(:, k) = pe1_ad(:, k) + temp_ad
          pe1_ad(:, k+1) = pe1_ad(:, k+1) + temp_ad
          logpl1_ad(:, k) = 0.0
          q2_ad(i1:i2, j, k) = q2_ad(i1:i2, j, k) + qx_ad(:, k)
          qx_ad(:, k) = 0.0
        END DO
      END IF
    ELSE IF (branch .EQ. 2) THEN
      DO k=km-1,1,-1
        logpl1_ad(:, k+1) = logpl1_ad(:, k+1) + dlogp1_ad(:, k)
        logpl1_ad(:, k) = logpl1_ad(:, k) - dlogp1_ad(:, k)
        dlogp1_ad(:, k) = 0.0
      END DO
      DO k=kn,1,-1
        temp_ad2 = logpl2_ad(:, k)/(pe2(:, k)+pe2(:, k+1))
        pe2_ad(:, k) = pe2_ad(:, k) + temp_ad2
        pe2_ad(:, k+1) = pe2_ad(:, k+1) + temp_ad2
        logpl2_ad(:, k) = 0.0
      END DO
      DO k=km,1,-1
        temp_ad1 = logpl1_ad(:, k)/(pe1(:, k)+pe1(:, k+1))
        pe1_ad(:, k) = pe1_ad(:, k) + temp_ad1
        pe1_ad(:, k+1) = pe1_ad(:, k+1) + temp_ad1
        logpl1_ad(:, k) = 0.0
        q2_ad(i1:i2, j, k) = q2_ad(i1:i2, j, k) + qx_ad(:, k)
        qx_ad(:, k) = 0.0
      END DO
    ELSE
      DO k=km-1,1,-1
        logpl1_ad(:, k+1) = logpl1_ad(:, k+1) + dlogp1_ad(:, k)
        logpl1_ad(:, k) = logpl1_ad(:, k) - dlogp1_ad(:, k)
        dlogp1_ad(:, k) = 0.0
      END DO
      DO k=kn,1,-1
        temp0 = r2*(pe2(:, k)+pe2(:, k+1))
        temp_ad4 = akap*EXP(akap*LOG(temp0))*r2*logpl2_ad(:, k)/temp0
        pe2_ad(:, k) = pe2_ad(:, k) + temp_ad4
        pe2_ad(:, k+1) = pe2_ad(:, k+1) + temp_ad4
        logpl2_ad(:, k) = 0.0
      END DO
      DO k=km,1,-1
        temp = r2*(pe1(:, k)+pe1(:, k+1))
        temp_ad3 = akap*EXP(akap*LOG(temp))*r2*logpl1_ad(:, k)/temp
        pe1_ad(:, k) = pe1_ad(:, k) + temp_ad3
        pe1_ad(:, k+1) = pe1_ad(:, k+1) + temp_ad3
        logpl1_ad(:, k) = 0.0
        q2_ad(i1:i2, j, k) = q2_ad(i1:i2, j, k) + qx_ad(:, k)
        qx_ad(:, k) = 0.0
      END DO
    END IF
  END SUBROUTINE MAP1_CUBIC_BWD
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
        logpl1(:, k) = EXP(akap*LOG(r2*(pe1(:, k)+pe1(:, k+1))))
      END DO
      DO k=1,kn
        logpl2(:, k) = EXP(akap*LOG(r2*(pe2(:, k)+pe2(:, k+1))))
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
!  Differentiation of map_scalar in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
!e_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_co
!re_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mo
!d.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayle
!igh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_
!ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod
!.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2
!d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limite
!rs fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic f
!v_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_sub
!grid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_uti
!ls_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils
!_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_
!mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod
!.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp
!_v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_co
!re_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_ut
!ils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 q2
!   with respect to varying inputs: pe1 pe2 q2
!-----------------------------------------------------------------------
  SUBROUTINE MAP_SCALAR_FWD(km, pe1, qs, kn, pe2, q2, i1, i2, j, ibeg&
&   , iend, jbeg, jend, iv, kord, q_min)
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
    INTEGER :: ad_count
    INTEGER :: ad_count0

    dp1 = 0.0
    q4 = 0.0
    pl = 0.0
    pr = 0.0
    qsum = 0.0
    dp = 0.0
    esl = 0.0
    i = 0
    k = 0
    l = 0
    m = 0
    k0 = 0
    ad_count = 0
    ad_count0 = 0

    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL SCALAR_PROFILE_FWD(qs, q4, dp1, km, i1, i2, iv, kord, &
&                          q_min)
!else
!     call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
    DO i=i1,i2
      k0 = 1
      DO 120 k=1,kn
        CALL PUSHINTEGER(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
        CALL PUSHINTEGER(ad_count)
        CALL PUSHCONTROL(2,2)
        GOTO 123
 100    CALL PUSHCONTROL(1,1)
        CALL PUSHINTEGER(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          CALL PUSHREALARRAY(q2(i, j, k))
          q2(i, j, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2&
&           , i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
          k0 = l
          CALL PUSHCONTROL(1,0)
          GOTO 120
        ELSE
! Fractional area...
          CALL PUSHREALARRAY(qsum)
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          CALL PUSHINTEGER(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
              qsum = qsum + dp1(i, m)*q4(1, i, m)
              CALL PUSHINTEGER(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL(1,0)
          CALL PUSHINTEGER(ad_count0)
          CALL PUSHCONTROL(2,1)
          GOTO 123
 110      CALL PUSHCONTROL(1,1)
          CALL PUSHINTEGER(ad_count0)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
          CALL PUSHCONTROL(2,0)
        END IF
 123    CALL PUSHREALARRAY(q2(i, j, k))
        q2(i, j, k) = qsum/(pe2(i, k+1)-pe2(i, k))
        CALL PUSHCONTROL(1,1)
 120  CONTINUE
    END DO
    CALL PUSHREALARRAY(q4, 4*(i2-i1+1)*km)
    CALL PUSHREALARRAY(qsum)
    CALL PUSHREALARRAY(dp1, (i2-i1+1)*km)
    CALL PUSHINTEGER(m)
    CALL PUSHINTEGER(l)
  END SUBROUTINE MAP_SCALAR_FWD
!  Differentiation of map_scalar in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
!ge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_c
!ore_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_m
!od.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayl
!eigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l
!_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mo
!d.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_
!2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limit
!ers fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic 
!fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_su
!bgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_ut
!ils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_util
!s_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils
!_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mo
!d.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.yt
!p_v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_c
!ore_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_u
!tils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 q2
!   with respect to varying inputs: pe1 pe2 q2
!-----------------------------------------------------------------------
  SUBROUTINE MAP_SCALAR_BWD(km, pe1, pe1_ad, qs, kn, pe2, pe2_ad, q2&
&   , q2_ad, i1, i2, j, ibeg, iend, jbeg, jend, iv, kord, q_min)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1
    INTEGER, INTENT(IN) :: i2
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: kn
    REAL, INTENT(IN) :: qs(i1:i2)
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL :: pe1_ad(i1:i2, km+1)
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL :: pe2_ad(i1:i2, kn+1)
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(INOUT) :: q2_ad(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(IN) :: q_min
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_ad(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: q4_ad(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    REAL :: pl_ad, pr_ad, qsum_ad, dp_ad, esl_ad
    INTEGER :: i, k, l, m, k0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp0
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp1
    REAL :: temp_ad11
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3

    dp1 = 0.0
    q4 = 0.0
    pl = 0.0
    pr = 0.0
    qsum = 0.0
    dp = 0.0
    esl = 0.0
    i = 0
    k = 0
    l = 0
    m = 0
    k0 = 0
    ad_count = 0
    ad_count0 = 0
    branch = 0

    CALL POPINTEGER(l)
    CALL POPINTEGER(m)
    CALL POPREALARRAY(dp1, (i2-i1+1)*km)
    CALL POPREALARRAY(qsum)
    CALL POPREALARRAY(q4, 4*(i2-i1+1)*km)
    dp1_ad = 0.0
    qsum_ad = 0.0
    q4_ad = 0.0
    DO i=i2,i1,-1
      DO k=kn,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY(q2(i, j, k))
          temp_ad0 = 0.5*(pr+pl)*q2_ad(i, j, k)
          temp_ad1 = 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i, l))*q2_ad(i, &
&           j, k)
          temp_ad2 = -(r3*q4(4, i, l)*q2_ad(i, j, k))
          q4_ad(2, i, l) = q4_ad(2, i, l) + q2_ad(i, j, k) - temp_ad0
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad0 - r3*(pr*(pr+pl)+pl&
&           **2)*q2_ad(i, j, k)
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad0
          pr_ad = (2*pr+pl)*temp_ad2 + temp_ad1
          pl_ad = (2*pl+pr)*temp_ad2 + temp_ad1
          q2_ad(i, j, k) = 0.0
          temp_ad3 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad3
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad3
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad3&
&           /dp1(i, l)
        ELSE
          CALL POPREALARRAY(q2(i, j, k))
          temp1 = pe2(i, k+1) - pe2(i, k)
          temp_ad11 = -(qsum*q2_ad(i, j, k)/temp1**2)
          qsum_ad = qsum_ad + q2_ad(i, j, k)/temp1
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad11
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad11
          q2_ad(i, j, k) = 0.0
          CALL POPCONTROL(2,branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            temp0 = q4(3, i, m) - q4(2, i, m) + q4(4, i, m)*(-(r23*esl)+&
&             1.)
            temp_ad8 = dp*qsum_ad
            temp_ad9 = 0.5*esl*temp_ad8
            q4_ad(2, i, m) = q4_ad(2, i, m) + temp_ad8 - temp_ad9
            esl_ad = 0.5*temp0*temp_ad8 - q4(4, i, m)*r23*temp_ad9
            q4_ad(3, i, m) = q4_ad(3, i, m) + temp_ad9
            q4_ad(4, i, m) = q4_ad(4, i, m) + (1.-r23*esl)*temp_ad9
            temp_ad10 = esl_ad/dp1(i, m)
            dp_ad = temp_ad10 + (q4(2, i, m)+0.5*(esl*temp0))*qsum_ad
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad10/dp1(i, m)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 100
          END IF
          CALL POPINTEGER(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL(1,branch)
            ELSE
              dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m)*qsum_ad
              q4_ad(1, i, m) = q4_ad(1, i, m) + dp1(i, m)*qsum_ad
            END IF
            CALL POPINTEGER(m)
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY(qsum)
          temp = q4(4, i, l) + q4(3, i, l) - q4(2, i, l)
          temp_ad4 = (q4(2, i, l)+0.5*(temp*(pl+1.))-r3*(q4(4, i, l)*(pl&
&           *(pl+1.)+1.)))*qsum_ad
          temp_ad5 = (pe1(i, l+1)-pe2(i, k))*qsum_ad
          temp_ad6 = 0.5*(pl+1.)*temp_ad5
          temp_ad7 = -(r3*q4(4, i, l)*temp_ad5)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + temp_ad4
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad4
          q4_ad(2, i, l) = q4_ad(2, i, l) + temp_ad5 - temp_ad6
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad6 - r3*(pl*(pl+1.)+1.&
&           )*temp_ad5
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad6
          pl_ad = (2*pl+1.)*temp_ad7 + 0.5*temp*temp_ad5
          qsum_ad = 0.0
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 100    CALL POPINTEGER(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL(1,branch)
          CALL POPINTEGER(l)
        END DO
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) CALL SCALAR_PROFILE_BWD(qs, q4, q4_ad, dp1, &
&                                           dp1_ad, km, i1, i2, iv, kord&
&                                           , q_min)
    DO k=km,1,-1
      DO i=i2,i1,-1
        q2_ad(i, j, k) = q2_ad(i, j, k) + q4_ad(1, i, k)
        q4_ad(1, i, k) = 0.0
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE MAP_SCALAR_BWD
!-----------------------------------------------------------------------
!  Differentiation of map1_ppm in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
!mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core
!_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.
!mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleig
!h_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_or
!d4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.r
!emap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d 
!fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters
! fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_
!restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgr
!id_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils
!_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_m
!od.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mo
!d.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d
!2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
!_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core
!_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_util
!s_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 qs q2
!   with respect to varying inputs: pe1 pe2 qs q2
  SUBROUTINE MAP1_PPM_FWD(km, pe1, qs, kn, pe2, q2, i1, i2, j, ibeg, &
&   iend, jbeg, jend, iv, kord)
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
    INTEGER :: ad_count
    INTEGER :: ad_count0

    dp1 = 0.0
    q4 = 0.0
    pl = 0.0
    pr = 0.0
    qsum = 0.0
    dp = 0.0
    esl = 0.0
    i = 0
    k = 0
    l = 0
    m = 0
    k0 = 0
    ad_count = 0
    ad_count0 = 0

    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE_FWD(qs, q4, dp1, km, i1, i2, iv, kord)
!else
!     call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
    DO i=i1,i2
      k0 = 1
      DO 120 k=1,kn
        CALL PUSHINTEGER(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
        CALL PUSHINTEGER(ad_count)
        CALL PUSHCONTROL(2,2)
        GOTO 123
 100    CALL PUSHCONTROL(1,1)
        CALL PUSHINTEGER(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          CALL PUSHREALARRAY(q2(i, j, k))
          q2(i, j, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2&
&           , i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
          k0 = l
          CALL PUSHCONTROL(1,0)
          GOTO 120
        ELSE
! Fractional area...
          CALL PUSHREALARRAY(qsum)
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          CALL PUSHINTEGER(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
              qsum = qsum + dp1(i, m)*q4(1, i, m)
              CALL PUSHINTEGER(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL(1,0)
          CALL PUSHINTEGER(ad_count0)
          CALL PUSHCONTROL(2,1)
          GOTO 123
 110      CALL PUSHCONTROL(1,1)
          CALL PUSHINTEGER(ad_count0)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
          CALL PUSHCONTROL(2,0)
        END IF
 123    CALL PUSHREALARRAY(q2(i, j, k))
        q2(i, j, k) = qsum/(pe2(i, k+1)-pe2(i, k))
        CALL PUSHCONTROL(1,1)
 120  CONTINUE
    END DO
    CALL PUSHREALARRAY(q4, 4*(i2-i1+1)*km)
    CALL PUSHREALARRAY(qsum)
    CALL PUSHREALARRAY(dp1, (i2-i1+1)*km)
    CALL PUSHINTEGER(m)
    CALL PUSHINTEGER(l)
  END SUBROUTINE MAP1_PPM_FWD
!  Differentiation of map1_ppm in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
!_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_cor
!e_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod
!.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Raylei
!gh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_o
!rd4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.
!remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d
! fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiter
!s fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv
!_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subg
!rid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_util
!s_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_
!mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_m
!od.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.
!d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_
!v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_cor
!e_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_uti
!ls_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 qs q2
!   with respect to varying inputs: pe1 pe2 qs q2
  SUBROUTINE MAP1_PPM_BWD(km, pe1, pe1_ad, qs, qs_ad, kn, pe2, pe2_ad&
&   , q2, q2_ad, i1, i2, j, ibeg, iend, jbeg, jend, iv, kord)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1
    INTEGER, INTENT(IN) :: i2
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: kn
    REAL, INTENT(IN) :: qs(i1:i2)
    REAL :: qs_ad(i1:i2)
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL :: pe1_ad(i1:i2, km+1)
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL :: pe2_ad(i1:i2, kn+1)
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, kn)
    REAL, INTENT(INOUT) :: q2_ad(ibeg:iend, jbeg:jend, kn)
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_ad(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: q4_ad(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    REAL :: pl_ad, pr_ad, qsum_ad, dp_ad, esl_ad
    INTEGER :: i, k, l, m, k0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp0
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp1
    REAL :: temp_ad11
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3

    dp1 = 0.0
    q4 = 0.0
    pl = 0.0
    pr = 0.0
    qsum = 0.0
    dp = 0.0
    esl = 0.0
    i = 0
    k = 0
    l = 0
    m = 0
    k0 = 0
    ad_count = 0
    ad_count0 = 0
    branch = 0

    CALL POPINTEGER(l)
    CALL POPINTEGER(m)
    CALL POPREALARRAY(dp1, (i2-i1+1)*km)
    CALL POPREALARRAY(qsum)
    CALL POPREALARRAY(q4, 4*(i2-i1+1)*km)
    dp1_ad = 0.0
    qsum_ad = 0.0
    q4_ad = 0.0
    DO i=i2,i1,-1
      DO k=kn,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY(q2(i, j, k))
          temp_ad0 = 0.5*(pr+pl)*q2_ad(i, j, k)
          temp_ad1 = 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i, l))*q2_ad(i, &
&           j, k)
          temp_ad2 = -(r3*q4(4, i, l)*q2_ad(i, j, k))
          q4_ad(2, i, l) = q4_ad(2, i, l) + q2_ad(i, j, k) - temp_ad0
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad0 - r3*(pr*(pr+pl)+pl&
&           **2)*q2_ad(i, j, k)
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad0
          pr_ad = (2*pr+pl)*temp_ad2 + temp_ad1
          pl_ad = (2*pl+pr)*temp_ad2 + temp_ad1
          q2_ad(i, j, k) = 0.0
          temp_ad3 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad3
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad3
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad3&
&           /dp1(i, l)
        ELSE
          CALL POPREALARRAY(q2(i, j, k))
          temp1 = pe2(i, k+1) - pe2(i, k)
          temp_ad11 = -(qsum*q2_ad(i, j, k)/temp1**2)
          qsum_ad = qsum_ad + q2_ad(i, j, k)/temp1
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad11
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad11
          q2_ad(i, j, k) = 0.0
          CALL POPCONTROL(2,branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            temp0 = q4(3, i, m) - q4(2, i, m) + q4(4, i, m)*(-(r23*esl)+&
&             1.)
            temp_ad8 = dp*qsum_ad
            temp_ad9 = 0.5*esl*temp_ad8
            q4_ad(2, i, m) = q4_ad(2, i, m) + temp_ad8 - temp_ad9
            esl_ad = 0.5*temp0*temp_ad8 - q4(4, i, m)*r23*temp_ad9
            q4_ad(3, i, m) = q4_ad(3, i, m) + temp_ad9
            q4_ad(4, i, m) = q4_ad(4, i, m) + (1.-r23*esl)*temp_ad9
            temp_ad10 = esl_ad/dp1(i, m)
            dp_ad = temp_ad10 + (q4(2, i, m)+0.5*(esl*temp0))*qsum_ad
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad10/dp1(i, m)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 100
          END IF
          CALL POPINTEGER(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL(1,branch)
            ELSE
              dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m)*qsum_ad
              q4_ad(1, i, m) = q4_ad(1, i, m) + dp1(i, m)*qsum_ad
            END IF
            CALL POPINTEGER(m)
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY(qsum)
          temp = q4(4, i, l) + q4(3, i, l) - q4(2, i, l)
          temp_ad4 = (q4(2, i, l)+0.5*(temp*(pl+1.))-r3*(q4(4, i, l)*(pl&
&           *(pl+1.)+1.)))*qsum_ad
          temp_ad5 = (pe1(i, l+1)-pe2(i, k))*qsum_ad
          temp_ad6 = 0.5*(pl+1.)*temp_ad5
          temp_ad7 = -(r3*q4(4, i, l)*temp_ad5)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + temp_ad4
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad4
          q4_ad(2, i, l) = q4_ad(2, i, l) + temp_ad5 - temp_ad6
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad6 - r3*(pl*(pl+1.)+1.&
&           )*temp_ad5
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad6
          pl_ad = (2*pl+1.)*temp_ad7 + 0.5*temp*temp_ad5
          qsum_ad = 0.0
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 100    CALL POPINTEGER(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL(1,branch)
          CALL POPINTEGER(l)
        END DO
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) CALL CS_PROFILE_BWD(qs, qs_ad, q4, q4_ad, dp1&
&                                       , dp1_ad, km, i1, i2, iv, kord)
    DO k=km,1,-1
      DO i=i2,i1,-1
        q2_ad(i, j, k) = q2_ad(i, j, k) + q4_ad(1, i, k)
        q4_ad(1, i, k) = 0.0
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE MAP1_PPM_BWD
!  Differentiation of mapn_tracer in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
!ge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_c
!ore_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_m
!od.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayl
!eigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l
!_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mo
!d.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_
!2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limit
!ers fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic 
!fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_su
!bgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_ut
!ils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_util
!s_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils
!_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mo
!d.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.yt
!p_v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_c
!ore_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_u
!tils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 dp2 q1
!   with respect to varying inputs: pe1 pe2 dp2 q1
  SUBROUTINE MAPN_TRACER_FWD(nq, km, pe1, pe2, q1, dp2, kord, j, i1, &
&   i2, isd, ied, jsd, jed, q_min, fill)
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
    INTEGER :: arg1
    INTEGER :: ad_count
    INTEGER :: ad_count0

    q4 = 0.0
    q2 = 0.0
    qsum = 0.0
    dp1 = 0.0
    qs = 0.0
    pl = 0.0
    pr = 0.0
    dp = 0.0
    esl = 0.0
    fac1 = 0.0
    fac2 = 0.0
    i = 0
    k = 0
    l = 0
    m = 0
    k0 = 0
    iq = 0
    arg1 = 0
    ad_count = 0
    ad_count0 = 0

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
      CALL SCALAR_PROFILE_FWD(qs, q4(1:4, i1:i2, 1:km, iq), dp1, km, &
&                          i1, i2, 0, kord(iq), q_min)
    END DO
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 130 k=1,km
        CALL PUSHINTEGER(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
        CALL PUSHINTEGER(ad_count)
        CALL PUSHCONTROL(2,2)
        GOTO 120
 100    CALL PUSHCONTROL(1,1)
        CALL PUSHINTEGER(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          CALL PUSHREALARRAY(fac1)
          fac1 = pr + pl
          CALL PUSHREALARRAY(fac2)
          fac2 = r3*(pr*fac1+pl*pl)
          CALL PUSHREALARRAY(fac1)
          fac1 = 0.5*fac1
          DO iq=1,nq
            q2(i, k, iq) = q4(2, i, l, iq) + (q4(4, i, l, iq)+q4(3, i, l&
&             , iq)-q4(2, i, l, iq))*fac1 - q4(4, i, l, iq)*fac2
          END DO
          k0 = l
          CALL PUSHCONTROL(1,0)
          GOTO 130
        ELSE
! Fractional area...
          CALL PUSHREALARRAY(dp)
          dp = pe1(i, l+1) - pe2(i, k)
          CALL PUSHREALARRAY(fac1)
          fac1 = 1. + pl
          CALL PUSHREALARRAY(fac2)
          fac2 = r3*(1.+pl*fac1)
          CALL PUSHREALARRAY(fac1)
          fac1 = 0.5*fac1
          DO iq=1,nq
            CALL PUSHREALARRAY(qsum(iq))
            qsum(iq) = dp*(q4(2, i, l, iq)+(q4(4, i, l, iq)+q4(3, i, l, &
&             iq)-q4(2, i, l, iq))*fac1-q4(4, i, l, iq)*fac2)
          END DO
          CALL PUSHINTEGER(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
              DO iq=1,nq
                CALL PUSHREALARRAY(qsum(iq))
                qsum(iq) = qsum(iq) + dp1(i, m)*q4(1, i, m, iq)
              END DO
              CALL PUSHINTEGER(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL(1,0)
          CALL PUSHINTEGER(ad_count0)
          CALL PUSHCONTROL(2,1)
          GOTO 120
 110      CALL PUSHCONTROL(1,1)
          CALL PUSHINTEGER(ad_count0)
          CALL PUSHREALARRAY(dp)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          CALL PUSHREALARRAY(fac1)
          fac1 = 0.5*esl
          CALL PUSHREALARRAY(fac2)
          fac2 = 1. - r23*esl
          DO iq=1,nq
            CALL PUSHREALARRAY(qsum(iq))
            qsum(iq) = qsum(iq) + dp*(q4(2, i, m, iq)+fac1*(q4(3, i, m, &
&             iq)-q4(2, i, m, iq)+q4(4, i, m, iq)*fac2))
          END DO
          k0 = m
          CALL PUSHCONTROL(2,0)
        END IF
 120    DO iq=1,nq
          q2(i, k, iq) = qsum(iq)/dp2(i, k)
        END DO
        CALL PUSHCONTROL(1,1)
 130  CONTINUE
    END DO
    DO iq=1,nq
!    if (fill) call fillz(i2-i1+1, km, 1, q2(i1,1,iq), dp2)
      DO k=1,km
        DO i=i1,i2
          CALL PUSHREALARRAY(q1(i, j, k, iq))
          q1(i, j, k, iq) = q2(i, k, iq)
        END DO
      END DO
    END DO
    CALL PUSHREALARRAY(q4, 4*(i2-i1+1)*km*nq)
    CALL PUSHREALARRAY(qsum, nq)
    CALL PUSHREALARRAY(dp1, (i2-i1+1)*km)
    CALL PUSHREALARRAY(qs, i2 - i1 + 1)
    CALL PUSHREALARRAY(fac2)
    CALL PUSHREALARRAY(fac1)
    CALL PUSHREALARRAY(dp)
    CALL PUSHINTEGER(m)
    CALL PUSHINTEGER(l)
  END SUBROUTINE MAPN_TRACER_FWD
!  Differentiation of mapn_tracer in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_e
!dge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_
!core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_
!mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Ray
!leigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2
!l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_m
!od.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap
!_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limi
!ters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic
! fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_s
!ubgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_u
!tils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_uti
!ls_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_util
!s_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_m
!od.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.y
!tp_v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_
!core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_
!utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 dp2 q1
!   with respect to varying inputs: pe1 pe2 dp2 q1
  SUBROUTINE MAPN_TRACER_BWD(nq, km, pe1, pe1_ad, pe2, pe2_ad, q1, &
&   q1_ad, dp2, dp2_ad, kord, j, i1, i2, isd, ied, jsd, jed, q_min, fill&
& )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: j, nq, i1, i2
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: kord(nq)
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL :: pe1_ad(i1:i2, km+1)
    REAL, INTENT(IN) :: pe2(i1:i2, km+1)
    REAL :: pe2_ad(i1:i2, km+1)
    REAL, INTENT(IN) :: dp2(i1:i2, km)
    REAL :: dp2_ad(i1:i2, km)
    REAL, INTENT(IN) :: q_min
    LOGICAL, INTENT(IN) :: fill
    REAL, INTENT(INOUT) :: q1(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: q1_ad(isd:ied, jsd:jed, km, nq)
    REAL :: q4(4, i1:i2, km, nq)
    REAL :: q4_ad(4, i1:i2, km, nq)
    REAL :: q2(i1:i2, km, nq)
    REAL :: q2_ad(i1:i2, km, nq)
    REAL :: qsum(nq)
    REAL :: qsum_ad(nq)
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_ad(i1:i2, km)
    REAL :: qs(i1:i2)
    REAL :: pl, pr, dp, esl, fac1, fac2
    REAL :: pl_ad, pr_ad, dp_ad, esl_ad, fac1_ad, fac2_ad
    INTEGER :: i, k, l, m, k0, iq
    INTEGER :: arg1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp0
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3

    q4 = 0.0
    q2 = 0.0
    qsum = 0.0
    dp1 = 0.0
    qs = 0.0
    pl = 0.0
    pr = 0.0
    dp = 0.0
    esl = 0.0
    fac1 = 0.0
    fac2 = 0.0
    i = 0
    k = 0
    l = 0
    m = 0
    k0 = 0
    iq = 0
    arg1 = 0
    ad_count = 0
    ad_count0 = 0
    branch = 0

    CALL POPINTEGER(l)
    CALL POPINTEGER(m)
    CALL POPREALARRAY(dp)
    CALL POPREALARRAY(fac1)
    CALL POPREALARRAY(fac2)
    CALL POPREALARRAY(qs, i2 - i1 + 1)
    CALL POPREALARRAY(dp1, (i2-i1+1)*km)
    CALL POPREALARRAY(qsum, nq)
    CALL POPREALARRAY(q4, 4*(i2-i1+1)*km*nq)
    q2_ad = 0.0
    DO iq=nq,1,-1
      DO k=km,1,-1
        DO i=i2,i1,-1
          CALL POPREALARRAY(q1(i, j, k, iq))
          q2_ad(i, k, iq) = q2_ad(i, k, iq) + q1_ad(i, j, k, iq)
          q1_ad(i, j, k, iq) = 0.0
        END DO
      END DO
    END DO
    dp1_ad = 0.0
    qsum_ad = 0.0
    q4_ad = 0.0
    DO i=i2,i1,-1
      DO k=km,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          fac1_ad = 0.0
          fac2_ad = 0.0
          DO iq=nq,1,-1
            temp_ad2 = fac1*q2_ad(i, k, iq)
            q4_ad(2, i, l, iq) = q4_ad(2, i, l, iq) + q2_ad(i, k, iq) - &
&             temp_ad2
            q4_ad(4, i, l, iq) = q4_ad(4, i, l, iq) + temp_ad2 - fac2*&
&             q2_ad(i, k, iq)
            q4_ad(3, i, l, iq) = q4_ad(3, i, l, iq) + temp_ad2
            fac1_ad = fac1_ad + (q4(4, i, l, iq)+q4(3, i, l, iq)-q4(2, i&
&             , l, iq))*q2_ad(i, k, iq)
            fac2_ad = fac2_ad - q4(4, i, l, iq)*q2_ad(i, k, iq)
            q2_ad(i, k, iq) = 0.0
          END DO
          temp_ad0 = r3*fac2_ad
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY(fac1)
          fac1_ad = pr*temp_ad0 + 0.5*fac1_ad
          CALL POPREALARRAY(fac2)
          pr_ad = fac1_ad + fac1*temp_ad0
          pl_ad = fac1_ad + 2*pl*temp_ad0
          CALL POPREALARRAY(fac1)
          temp_ad1 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad1
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad1
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad1&
&           /dp1(i, l)
        ELSE
          DO iq=nq,1,-1
            temp_ad8 = q2_ad(i, k, iq)/dp2(i, k)
            qsum_ad(iq) = qsum_ad(iq) + temp_ad8
            dp2_ad(i, k) = dp2_ad(i, k) - qsum(iq)*temp_ad8/dp2(i, k)
            q2_ad(i, k, iq) = 0.0
          END DO
          CALL POPCONTROL(2,branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            fac1 = 0.5*esl
            fac2 = 1. - r23*esl
            dp_ad = 0.0
            fac1_ad = 0.0
            fac2_ad = 0.0
            DO iq=nq,1,-1
              CALL POPREALARRAY(qsum(iq))
              temp0 = q4(3, i, m, iq) - q4(2, i, m, iq) + q4(4, i, m, iq&
&               )*fac2
              temp_ad6 = dp*qsum_ad(iq)
              temp_ad7 = fac1*temp_ad6
              dp_ad = dp_ad + (q4(2, i, m, iq)+fac1*temp0)*qsum_ad(iq)
              q4_ad(2, i, m, iq) = q4_ad(2, i, m, iq) + temp_ad6 - &
&               temp_ad7
              fac1_ad = fac1_ad + temp0*temp_ad6
              q4_ad(3, i, m, iq) = q4_ad(3, i, m, iq) + temp_ad7
              q4_ad(4, i, m, iq) = q4_ad(4, i, m, iq) + fac2*temp_ad7
              fac2_ad = fac2_ad + q4(4, i, m, iq)*temp_ad7
            END DO
            CALL POPREALARRAY(fac2)
            esl_ad = 0.5*fac1_ad - r23*fac2_ad
            CALL POPREALARRAY(fac1)
            temp_ad5 = esl_ad/dp1(i, m)
            dp_ad = dp_ad + temp_ad5
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad5/dp1(i, m)
            CALL POPREALARRAY(dp)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 100
          END IF
          CALL POPINTEGER(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL(1,branch)
            ELSE
              DO iq=nq,1,-1
                CALL POPREALARRAY(qsum(iq))
                dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m, iq)*qsum_ad(iq&
&                 )
                q4_ad(1, i, m, iq) = q4_ad(1, i, m, iq) + dp1(i, m)*&
&                 qsum_ad(iq)
              END DO
            END IF
            CALL POPINTEGER(m)
          END DO
          dp_ad = 0.0
          fac1_ad = 0.0
          fac2_ad = 0.0
          DO iq=nq,1,-1
            CALL POPREALARRAY(qsum(iq))
            temp = q4(4, i, l, iq) + q4(3, i, l, iq) - q4(2, i, l, iq)
            temp_ad3 = dp*qsum_ad(iq)
            temp_ad4 = fac1*temp_ad3
            dp_ad = dp_ad + (q4(2, i, l, iq)+temp*fac1-q4(4, i, l, iq)*&
&             fac2)*qsum_ad(iq)
            q4_ad(2, i, l, iq) = q4_ad(2, i, l, iq) + temp_ad3 - &
&             temp_ad4
            q4_ad(4, i, l, iq) = q4_ad(4, i, l, iq) + temp_ad4 - fac2*&
&             temp_ad3
            q4_ad(3, i, l, iq) = q4_ad(3, i, l, iq) + temp_ad4
            fac1_ad = fac1_ad + temp*temp_ad3
            fac2_ad = fac2_ad - q4(4, i, l, iq)*temp_ad3
            qsum_ad(iq) = 0.0
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY(fac1)
          fac1_ad = r3*pl*fac2_ad + 0.5*fac1_ad
          CALL POPREALARRAY(fac2)
          pl_ad = fac1_ad + r3*fac1*fac2_ad
          CALL POPREALARRAY(fac1)
          CALL POPREALARRAY(dp)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + dp_ad
          pe2_ad(i, k) = pe2_ad(i, k) - dp_ad
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 100    CALL POPINTEGER(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL(1,branch)
          CALL POPINTEGER(l)
        END DO
      END DO
    END DO
    DO iq=nq,1,-1
      CALL SCALAR_PROFILE_BWD(qs, q4(1:4, i1:i2, 1:km, iq), q4_ad(1:4&
&                          , i1:i2, 1:km, iq), dp1, dp1_ad, km, i1, i2, &
&                          0, kord(iq), q_min)
      DO k=km,1,-1
        DO i=i2,i1,-1
          q1_ad(i, j, k, iq) = q1_ad(i, j, k, iq) + q4_ad(1, i, k, iq)
          q4_ad(1, i, k, iq) = 0.0
        END DO
      END DO
    END DO
    DO k=km,1,-1
      DO i=i2,i1,-1
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE MAPN_TRACER_BWD
!  Differentiation of map1_q2 in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_m
!od.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_
!mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.m
!ix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh
!_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord
!4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.re
!map_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d f
!v_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters 
!fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_r
!estart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgri
!d_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_
!mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mo
!d.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod
!.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2
!a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v_
!fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_
!mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils
!_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 dp2 q1 q2
!   with respect to varying inputs: pe1 pe2 dp2 q1 q2
  SUBROUTINE MAP1_Q2_FWD(km, pe1, q1, kn, pe2, q2, dp2, i1, i2, iv, &
&   kord, j, ibeg, iend, jbeg, jend, q_min)
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
    INTEGER :: ad_count
    INTEGER :: ad_count0

    qs = 0.0
    dp1 = 0.0
    q4 = 0.0
    pl = 0.0
    pr = 0.0
    qsum = 0.0
    dp = 0.0
    esl = 0.0
    i = 0
    k = 0
    l = 0
    m = 0
    k0 = 0
    ad_count = 0
    ad_count0 = 0

    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q1(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL SCALAR_PROFILE_FWD(qs, q4, dp1, km, i1, i2, iv, kord, &
&                          q_min)
!else
!call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 120 k=1,kn
        CALL PUSHINTEGER(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
        CALL PUSHINTEGER(ad_count)
        CALL PUSHCONTROL(2,2)
        GOTO 123
 100    CALL PUSHCONTROL(1,1)
        CALL PUSHINTEGER(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2(i, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i&
&           , l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
          k0 = l
          CALL PUSHCONTROL(1,0)
          GOTO 120
        ELSE
! Fractional area...
          CALL PUSHREALARRAY(qsum)
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          CALL PUSHINTEGER(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
              qsum = qsum + dp1(i, m)*q4(1, i, m)
              CALL PUSHINTEGER(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL(1,0)
          CALL PUSHINTEGER(ad_count0)
          CALL PUSHCONTROL(2,1)
          GOTO 123
 110      CALL PUSHCONTROL(1,1)
          CALL PUSHINTEGER(ad_count0)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
          CALL PUSHCONTROL(2,0)
        END IF
 123    q2(i, k) = qsum/dp2(i, k)
        CALL PUSHCONTROL(1,1)
 120  CONTINUE
    END DO
    CALL PUSHREALARRAY(q4, 4*(i2-i1+1)*km)
    CALL PUSHREALARRAY(qsum)
    CALL PUSHREALARRAY(dp1, (i2-i1+1)*km)
    CALL PUSHREALARRAY(qs, i2 - i1 + 1)
    CALL PUSHINTEGER(m)
    CALL PUSHINTEGER(l)
  END SUBROUTINE MAP1_Q2_FWD
!  Differentiation of map1_q2 in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
!mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core
!_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.
!mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleig
!h_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_or
!d4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.r
!emap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d 
!fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters
! fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_
!restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgr
!id_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils
!_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_m
!od.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mo
!d.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d
!2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
!_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core
!_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_util
!s_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pe1 pe2 dp2 q1 q2
!   with respect to varying inputs: pe1 pe2 dp2 q1 q2
  SUBROUTINE MAP1_Q2_BWD(km, pe1, pe1_ad, q1, q1_ad, kn, pe2, pe2_ad&
&   , q2, q2_ad, dp2, dp2_ad, i1, i2, iv, kord, j, ibeg, iend, jbeg, &
&   jend, q_min)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: i1, i2
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: kn
    REAL, INTENT(IN) :: pe1(i1:i2, km+1)
    REAL :: pe1_ad(i1:i2, km+1)
    REAL, INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL :: pe2_ad(i1:i2, kn+1)
    REAL, INTENT(IN) :: q1(ibeg:iend, jbeg:jend, km)
    REAL :: q1_ad(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(IN) :: dp2(i1:i2, kn)
    REAL :: dp2_ad(i1:i2, kn)
    REAL, INTENT(IN) :: q_min
    REAL, INTENT(INOUT) :: q2(i1:i2, kn)
    REAL, INTENT(INOUT) :: q2_ad(i1:i2, kn)
    REAL :: qs(i1:i2)
    REAL :: dp1(i1:i2, km)
    REAL :: dp1_ad(i1:i2, km)
    REAL :: q4(4, i1:i2, km)
    REAL :: q4_ad(4, i1:i2, km)
    REAL :: pl, pr, qsum, dp, esl
    REAL :: pl_ad, pr_ad, qsum_ad, dp_ad, esl_ad
    INTEGER :: i, k, l, m, k0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp0
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3

    qs = 0.0
    dp1 = 0.0
    q4 = 0.0
    pl = 0.0
    pr = 0.0
    qsum = 0.0
    dp = 0.0
    esl = 0.0
    i = 0
    k = 0
    l = 0
    m = 0
    k0 = 0
    ad_count = 0
    ad_count0 = 0
    branch = 0

    CALL POPINTEGER(l)
    CALL POPINTEGER(m)
    CALL POPREALARRAY(qs, i2 - i1 + 1)
    CALL POPREALARRAY(dp1, (i2-i1+1)*km)
    CALL POPREALARRAY(qsum)
    CALL POPREALARRAY(q4, 4*(i2-i1+1)*km)
    dp1_ad = 0.0
    qsum_ad = 0.0
    q4_ad = 0.0
    DO i=i2,i1,-1
      DO k=kn,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          temp_ad0 = 0.5*(pr+pl)*q2_ad(i, k)
          temp_ad1 = 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i, l))*q2_ad(i, &
&           k)
          temp_ad2 = -(r3*q4(4, i, l)*q2_ad(i, k))
          q4_ad(2, i, l) = q4_ad(2, i, l) + q2_ad(i, k) - temp_ad0
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad0 - r3*(pr*(pr+pl)+pl&
&           **2)*q2_ad(i, k)
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad0
          pr_ad = (2*pr+pl)*temp_ad2 + temp_ad1
          pl_ad = (2*pl+pr)*temp_ad2 + temp_ad1
          q2_ad(i, k) = 0.0
          temp_ad3 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad3
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad3
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad3&
&           /dp1(i, l)
        ELSE
          temp_ad11 = q2_ad(i, k)/dp2(i, k)
          qsum_ad = qsum_ad + temp_ad11
          dp2_ad(i, k) = dp2_ad(i, k) - qsum*temp_ad11/dp2(i, k)
          q2_ad(i, k) = 0.0
          CALL POPCONTROL(2,branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            temp0 = q4(3, i, m) - q4(2, i, m) + q4(4, i, m)*(-(r23*esl)+&
&             1.)
            temp_ad8 = dp*qsum_ad
            temp_ad9 = 0.5*esl*temp_ad8
            q4_ad(2, i, m) = q4_ad(2, i, m) + temp_ad8 - temp_ad9
            esl_ad = 0.5*temp0*temp_ad8 - q4(4, i, m)*r23*temp_ad9
            q4_ad(3, i, m) = q4_ad(3, i, m) + temp_ad9
            q4_ad(4, i, m) = q4_ad(4, i, m) + (1.-r23*esl)*temp_ad9
            temp_ad10 = esl_ad/dp1(i, m)
            dp_ad = temp_ad10 + (q4(2, i, m)+0.5*(esl*temp0))*qsum_ad
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad10/dp1(i, m)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 100
          END IF
          CALL POPINTEGER(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL(1,branch)
            ELSE
              dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m)*qsum_ad
              q4_ad(1, i, m) = q4_ad(1, i, m) + dp1(i, m)*qsum_ad
            END IF
            CALL POPINTEGER(m)
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREALARRAY(qsum)
          temp = q4(4, i, l) + q4(3, i, l) - q4(2, i, l)
          temp_ad4 = (q4(2, i, l)+0.5*(temp*(pl+1.))-r3*(q4(4, i, l)*(pl&
&           *(pl+1.)+1.)))*qsum_ad
          temp_ad5 = (pe1(i, l+1)-pe2(i, k))*qsum_ad
          temp_ad6 = 0.5*(pl+1.)*temp_ad5
          temp_ad7 = -(r3*q4(4, i, l)*temp_ad5)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + temp_ad4
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad4
          q4_ad(2, i, l) = q4_ad(2, i, l) + temp_ad5 - temp_ad6
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad6 - r3*(pl*(pl+1.)+1.&
&           )*temp_ad5
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad6
          pl_ad = (2*pl+1.)*temp_ad7 + 0.5*temp*temp_ad5
          qsum_ad = 0.0
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 100    CALL POPINTEGER(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL(1,branch)
          CALL POPINTEGER(l)
        END DO
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) CALL SCALAR_PROFILE_BWD(qs, q4, q4_ad, dp1, &
&                                           dp1_ad, km, i1, i2, iv, kord&
&                                           , q_min)
    DO k=km,1,-1
      DO i=i2,i1,-1
        q1_ad(i, j, k) = q1_ad(i, j, k) + q4_ad(1, i, k)
        q4_ad(1, i, k) = 0.0
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE MAP1_Q2_BWD
!  Differentiation of scalar_profile in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b
!_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dy
!n_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_cor
!e_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.R
!ayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.
!c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz
!_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.rem
!ap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_li
!miters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cub
!ic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv
!_subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh
!_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_u
!tils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_ut
!ils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core
!_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod
!.ytp_v sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d t
!p_core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_gri
!d_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: delp a4
!   with respect to varying inputs: delp a4
  SUBROUTINE SCALAR_PROFILE_FWD(qs, a4, delp, km, i1, i2, iv, kord, &
&   qmin)
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
    INTEGER :: abs0

    gam = 0.0
    q = 0.0
    d4 = 0.0
    bet = 0.0
    a_bot = 0.0
    grat = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    i = 0
    k = 0
    im = 0
    abs0 = 0

    IF (iv .EQ. -2) THEN
      DO i=i1,i2
        gam(i, 2) = 0.5
        q(i, 1) = 1.5*a4(1, i, 1)
      END DO
      DO k=2,km-1
        DO i=i1,i2
          grat = delp(i, k-1)/delp(i, k)
          CALL PUSHREALARRAY(bet)
          bet = 2. + grat + grat - gam(i, k)
          CALL PUSHREALARRAY(q(i, k))
          q(i, k) = (3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))/bet
          gam(i, k+1) = grat/bet
        END DO
      END DO
      DO i=i1,i2
        grat = delp(i, km-1)/delp(i, km)
        CALL PUSHREALARRAY(q(i, km))
        q(i, km) = (3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-&
&         1))/(2.+grat+grat-gam(i, km))
        CALL PUSHREALARRAY(q(i, km+1))
        q(i, km+1) = qs(i)
      END DO
      DO k=km-1,1,-1
        DO i=i1,i2
          CALL PUSHREALARRAY(q(i, k))
          q(i, k) = q(i, k) - gam(i, k+1)*q(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
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
          CALL PUSHREALARRAY(d4(i))
          d4(i) = delp(i, k-1)/delp(i, k)
          CALL PUSHREALARRAY(bet)
          bet = 2. + d4(i) + d4(i) - gam(i, k-1)
          CALL PUSHREALARRAY(q(i, k))
          q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
          gam(i, k) = d4(i)/bet
        END DO
      END DO
      DO i=i1,i2
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        CALL PUSHREALARRAY(q(i, km+1))
        q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          CALL PUSHREALARRAY(q(i, k))
          q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
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
          CALL PUSHREALARRAY(a4(2, i, k))
          a4(2, i, k) = q(i, k)
          CALL PUSHREALARRAY(a4(3, i, k))
          a4(3, i, k) = q(i, k+1)
          CALL PUSHREALARRAY(a4(4, i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
      END DO
      CALL PUSHREALARRAY(gam, (i2-i1+1)*km)
      CALL PUSHREALARRAY(d4, i2 - i1 + 1)
      CALL PUSHREALARRAY(bet)
      CALL PUSHREALARRAY(q, (i2-i1+1)*(km+1))
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHREALARRAY(gam, (i2-i1+1)*km)
      CALL PUSHREALARRAY(d4, i2 - i1 + 1)
      CALL PUSHREALARRAY(bet)
      CALL PUSHREALARRAY(q, (i2-i1+1)*(km+1))
      CALL PUSHCONTROL(1,0)
    END IF
  END SUBROUTINE SCALAR_PROFILE_FWD
!  Differentiation of scalar_profile in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2
!b_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe d
!yn_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_co
!re_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.
!Rayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod
!.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_map
!z_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.re
!map_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_l
!imiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cu
!bic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.f
!v_subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d n
!h_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_
!utils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_u
!tils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_cor
!e_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mo
!d.ytp_v sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d 
!tp_core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_gr
!id_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: delp a4
!   with respect to varying inputs: delp a4
  SUBROUTINE SCALAR_PROFILE_BWD(qs, a4, a4_ad, delp, delp_ad, km, i1&
&   , i2, iv, kord, qmin)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1, i2
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
    REAL, INTENT(IN) :: qs(i1:i2)
    REAL, INTENT(IN) :: delp(i1:i2, km)
    REAL :: delp_ad(i1:i2, km)
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_ad(4, i1:i2, km)
    REAL, INTENT(IN) :: qmin
    LOGICAL, DIMENSION(i1:i2, km) :: extm, ext6
    REAL :: gam(i1:i2, km)
    REAL :: gam_ad(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: q_ad(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: d4_ad(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: bet_ad, a_bot_ad, grat_ad
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    INTEGER :: i, k, im
    INTRINSIC ABS
    INTEGER :: abs0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp0
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp1
    REAL :: temp2
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp_ad14
    INTEGER :: branch

    gam = 0.0
    q = 0.0
    d4 = 0.0
    bet = 0.0
    a_bot = 0.0
    grat = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    i = 0
    k = 0
    im = 0
    abs0 = 0
    branch = 0

    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(q, (i2-i1+1)*(km+1))
      CALL POPREALARRAY(bet)
      CALL POPREALARRAY(d4, i2 - i1 + 1)
      CALL POPREALARRAY(gam, (i2-i1+1)*km)
      q_ad = 0.0
    ELSE
      CALL POPREALARRAY(q, (i2-i1+1)*(km+1))
      CALL POPREALARRAY(bet)
      CALL POPREALARRAY(d4, i2 - i1 + 1)
      CALL POPREALARRAY(gam, (i2-i1+1)*km)
      q_ad = 0.0
      DO k=km,1,-1
        DO i=i2,i1,-1
          CALL POPREALARRAY(a4(4, i, k))
          temp_ad14 = 3.*a4_ad(4, i, k)
          a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad14
          a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad14
          a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad14
          a4_ad(4, i, k) = 0.0
          CALL POPREALARRAY(a4(3, i, k))
          q_ad(i, k+1) = q_ad(i, k+1) + a4_ad(3, i, k)
          a4_ad(3, i, k) = 0.0
          CALL POPREALARRAY(a4(2, i, k))
          q_ad(i, k) = q_ad(i, k) + a4_ad(2, i, k)
          a4_ad(2, i, k) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      gam_ad = 0.0
      DO k=1,km,1
        DO i=i2,i1,-1
          CALL POPREALARRAY(q(i, k))
          gam_ad(i, k) = gam_ad(i, k) - q(i, k+1)*q_ad(i, k)
          q_ad(i, k+1) = q_ad(i, k+1) - gam(i, k)*q_ad(i, k)
        END DO
      END DO
      d4_ad = 0.0
      DO i=i2,i1,-1
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        CALL POPREALARRAY(q(i, km+1))
        temp2 = d4(i)*(d4(i)+0.5) - a_bot*gam(i, km)
        temp_ad11 = q_ad(i, km+1)/temp2
        temp1 = d4(i)*(d4(i)+1.)
        temp_ad12 = 2.*a4(1, i, km)*temp_ad11
        temp_ad13 = -((2.*(temp1*a4(1, i, km))+a4(1, i, km-1)-a_bot*q(i&
&         , km))*temp_ad11/temp2)
        a4_ad(1, i, km) = a4_ad(1, i, km) + 2.*temp1*temp_ad11
        a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + temp_ad11
        a_bot_ad = -(gam(i, km)*temp_ad13) - q(i, km)*temp_ad11
        d4_ad(i) = d4_ad(i) + (2*d4(i)+1.5)*a_bot_ad + (2*d4(i)+0.5)*&
&         temp_ad13 + (2*d4(i)+1.)*temp_ad12
        q_ad(i, km) = q_ad(i, km) - a_bot*temp_ad11
        gam_ad(i, km) = gam_ad(i, km) - a_bot*temp_ad13
        q_ad(i, km+1) = 0.0
      END DO
      DO k=km,2,-1
        DO i=i2,i1,-1
          temp_ad9 = q_ad(i, k)/bet
          temp_ad8 = 3.*temp_ad9
          CALL POPREALARRAY(q(i, k))
          bet_ad = -((3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))*&
&           temp_ad9/bet) - d4(i)*gam_ad(i, k)/bet**2
          d4_ad(i) = d4_ad(i) + a4(1, i, k)*temp_ad8 + 2*bet_ad + gam_ad&
&           (i, k)/bet
          gam_ad(i, k) = 0.0
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + temp_ad8
          a4_ad(1, i, k) = a4_ad(1, i, k) + d4(i)*temp_ad8
          q_ad(i, k-1) = q_ad(i, k-1) - temp_ad9
          q_ad(i, k) = 0.0
          CALL POPREALARRAY(bet)
          gam_ad(i, k-1) = gam_ad(i, k-1) - bet_ad
          CALL POPREALARRAY(d4(i))
          temp_ad10 = d4_ad(i)/delp(i, k)
          delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad10
          delp_ad(i, k) = delp_ad(i, k) - delp(i, k-1)*temp_ad10/delp(i&
&           , k)
          d4_ad(i) = 0.0
        END DO
      END DO
      DO i=i2,i1,-1
        grat = delp(i, 2)/delp(i, 1)
        bet = grat*(grat+0.5)
        temp_ad4 = gam_ad(i, 1)/bet
        gam_ad(i, 1) = 0.0
        temp_ad6 = q_ad(i, 1)/bet
        temp_ad5 = a4(1, i, 1)*temp_ad6
        temp0 = 2*grat*(grat+1.)
        bet_ad = -((temp0*a4(1, i, 1)+a4(1, i, 2))*temp_ad6/bet) - (grat&
&         *(grat+1.5)+1.)*temp_ad4/bet
        grat_ad = (4*grat+2*1.)*temp_ad5 + (2*grat+0.5)*bet_ad + (2*grat&
&         +1.5)*temp_ad4
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + temp0*temp_ad6
        a4_ad(1, i, 2) = a4_ad(1, i, 2) + temp_ad6
        q_ad(i, 1) = 0.0
        temp_ad7 = grat_ad/delp(i, 1)
        delp_ad(i, 2) = delp_ad(i, 2) + temp_ad7
        delp_ad(i, 1) = delp_ad(i, 1) - delp(i, 2)*temp_ad7/delp(i, 1)
      END DO
    ELSE
      gam_ad = 0.0
      DO k=1,km-1,1
        DO i=i2,i1,-1
          CALL POPREALARRAY(q(i, k))
          gam_ad(i, k+1) = gam_ad(i, k+1) - q(i, k+1)*q_ad(i, k)
          q_ad(i, k+1) = q_ad(i, k+1) - gam(i, k+1)*q_ad(i, k)
        END DO
      END DO
      DO i=i2,i1,-1
        CALL POPREALARRAY(q(i, km+1))
        q_ad(i, km+1) = 0.0
        grat = delp(i, km-1)/delp(i, km)
        CALL POPREALARRAY(q(i, km))
        temp = 2*grat - gam(i, km) + 2.
        temp_ad1 = q_ad(i, km)/temp
        temp_ad2 = -((3.*(a4(1, i, km-1)+a4(1, i, km))-qs(i)*grat-q(i, &
&         km-1))*temp_ad1/temp)
        a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + 3.*temp_ad1
        a4_ad(1, i, km) = a4_ad(1, i, km) + 3.*temp_ad1
        grat_ad = 2*temp_ad2 - qs(i)*temp_ad1
        q_ad(i, km-1) = q_ad(i, km-1) - temp_ad1
        gam_ad(i, km) = gam_ad(i, km) - temp_ad2
        q_ad(i, km) = 0.0
        temp_ad3 = grat_ad/delp(i, km)
        delp_ad(i, km-1) = delp_ad(i, km-1) + temp_ad3
        delp_ad(i, km) = delp_ad(i, km) - delp(i, km-1)*temp_ad3/delp(i&
&         , km)
      END DO
      DO k=km-1,2,-1
        DO i=i2,i1,-1
          temp_ad = q_ad(i, k)/bet
          CALL POPREALARRAY(q(i, k))
          grat = delp(i, k-1)/delp(i, k)
          bet_ad = -((3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))*temp_ad/&
&           bet) - grat*gam_ad(i, k+1)/bet**2
          grat_ad = 2*bet_ad + gam_ad(i, k+1)/bet
          gam_ad(i, k+1) = 0.0
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + 3.*temp_ad
          a4_ad(1, i, k) = a4_ad(1, i, k) + 3.*temp_ad
          q_ad(i, k-1) = q_ad(i, k-1) - temp_ad
          q_ad(i, k) = 0.0
          CALL POPREALARRAY(bet)
          gam_ad(i, k) = gam_ad(i, k) - bet_ad
          temp_ad0 = grat_ad/delp(i, k)
          delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad0
          delp_ad(i, k) = delp_ad(i, k) - delp(i, k-1)*temp_ad0/delp(i, &
&           k)
        END DO
      END DO
      DO i=i2,i1,-1
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + 1.5*q_ad(i, 1)
        q_ad(i, 1) = 0.0
      END DO
    END IF
  END SUBROUTINE SCALAR_PROFILE_BWD
!  Differentiation of cs_profile in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
!e_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_co
!re_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mo
!d.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayle
!igh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_
!ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod
!.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2
!d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limite
!rs fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic f
!v_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_sub
!grid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_uti
!ls_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils
!_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_
!mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod
!.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp
!_v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_co
!re_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_ut
!ils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: qs delp a4
!   with respect to varying inputs: qs delp a4
  SUBROUTINE CS_PROFILE_FWD(qs, a4, delp, km, i1, i2, iv, kord)
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
    INTEGER :: abs0

    gam = 0.0
    q = 0.0
    d4 = 0.0
    bet = 0.0
    a_bot = 0.0
    grat = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    i = 0
    k = 0
    im = 0
    abs0 = 0

    IF (iv .EQ. -2) THEN
      DO i=i1,i2
        gam(i, 2) = 0.5
        q(i, 1) = 1.5*a4(1, i, 1)
      END DO
      DO k=2,km-1
        DO i=i1,i2
          grat = delp(i, k-1)/delp(i, k)
          CALL PUSHREALARRAY(bet)
          bet = 2. + grat + grat - gam(i, k)
          CALL PUSHREALARRAY(q(i, k))
          q(i, k) = (3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))/bet
          gam(i, k+1) = grat/bet
        END DO
      END DO
      DO i=i1,i2
        grat = delp(i, km-1)/delp(i, km)
        CALL PUSHREALARRAY(q(i, km))
        q(i, km) = (3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, km-&
&         1))/(2.+grat+grat-gam(i, km))
        CALL PUSHREALARRAY(q(i, km+1))
        q(i, km+1) = qs(i)
      END DO
      DO k=km-1,1,-1
        DO i=i1,i2
          CALL PUSHREALARRAY(q(i, k))
          q(i, k) = q(i, k) - gam(i, k+1)*q(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
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
          CALL PUSHREALARRAY(d4(i))
          d4(i) = delp(i, k-1)/delp(i, k)
          CALL PUSHREALARRAY(bet)
          bet = 2. + d4(i) + d4(i) - gam(i, k-1)
          CALL PUSHREALARRAY(q(i, k))
          q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
          gam(i, k) = d4(i)/bet
        END DO
      END DO
      DO i=i1,i2
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        CALL PUSHREALARRAY(q(i, km+1))
        q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&         a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          CALL PUSHREALARRAY(q(i, k))
          q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
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
          CALL PUSHREALARRAY(a4(2, i, k))
          a4(2, i, k) = q(i, k)
          CALL PUSHREALARRAY(a4(3, i, k))
          a4(3, i, k) = q(i, k+1)
          CALL PUSHREALARRAY(a4(4, i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
      END DO
      CALL PUSHREALARRAY(gam, (i2-i1+1)*km)
      CALL PUSHREALARRAY(d4, i2 - i1 + 1)
      CALL PUSHREALARRAY(bet)
      CALL PUSHREALARRAY(q, (i2-i1+1)*(km+1))
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHREALARRAY(gam, (i2-i1+1)*km)
      CALL PUSHREALARRAY(d4, i2 - i1 + 1)
      CALL PUSHREALARRAY(bet)
      CALL PUSHREALARRAY(q, (i2-i1+1)*(km+1))
      CALL PUSHCONTROL(1,0)
    END IF
  END SUBROUTINE CS_PROFILE_FWD
!  Differentiation of cs_profile in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
!ge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_c
!ore_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_m
!od.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayl
!eigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l
!_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mo
!d.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_
!2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limit
!ers fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic 
!fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_su
!bgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_ut
!ils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_util
!s_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils
!_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mo
!d.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.yt
!p_v_fb sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_c
!ore_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_u
!tils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: qs delp a4
!   with respect to varying inputs: qs delp a4
  SUBROUTINE CS_PROFILE_BWD(qs, qs_ad, a4, a4_ad, delp, delp_ad, km, &
&   i1, i2, iv, kord)
    IMPLICIT NONE
!----- Perfectly linear scheme --------------------------------
    INTEGER, INTENT(IN) :: i1, i2
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
    REAL, INTENT(IN) :: qs(i1:i2)
    REAL :: qs_ad(i1:i2)
    REAL, INTENT(IN) :: delp(i1:i2, km)
    REAL :: delp_ad(i1:i2, km)
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_ad(4, i1:i2, km)
    LOGICAL :: extm(i1:i2, km)
    REAL :: gam(i1:i2, km)
    REAL :: gam_ad(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: q_ad(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: d4_ad(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: bet_ad, a_bot_ad, grat_ad
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    INTEGER :: i, k, im
    INTRINSIC ABS
    INTEGER :: abs0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp0
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp1
    REAL :: temp2
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp_ad14
    INTEGER :: branch

    gam = 0.0
    q = 0.0
    d4 = 0.0
    bet = 0.0
    a_bot = 0.0
    grat = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    i = 0
    k = 0
    im = 0
    abs0 = 0
    branch = 0

    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(q, (i2-i1+1)*(km+1))
      CALL POPREALARRAY(bet)
      CALL POPREALARRAY(d4, i2 - i1 + 1)
      CALL POPREALARRAY(gam, (i2-i1+1)*km)
      q_ad = 0.0
    ELSE
      CALL POPREALARRAY(q, (i2-i1+1)*(km+1))
      CALL POPREALARRAY(bet)
      CALL POPREALARRAY(d4, i2 - i1 + 1)
      CALL POPREALARRAY(gam, (i2-i1+1)*km)
      q_ad = 0.0
      DO k=km,1,-1
        DO i=i2,i1,-1
          CALL POPREALARRAY(a4(4, i, k))
          temp_ad14 = 3.*a4_ad(4, i, k)
          a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad14
          a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad14
          a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad14
          a4_ad(4, i, k) = 0.0
          CALL POPREALARRAY(a4(3, i, k))
          q_ad(i, k+1) = q_ad(i, k+1) + a4_ad(3, i, k)
          a4_ad(3, i, k) = 0.0
          CALL POPREALARRAY(a4(2, i, k))
          q_ad(i, k) = q_ad(i, k) + a4_ad(2, i, k)
          a4_ad(2, i, k) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      gam_ad = 0.0
      DO k=1,km,1
        DO i=i2,i1,-1
          CALL POPREALARRAY(q(i, k))
          gam_ad(i, k) = gam_ad(i, k) - q(i, k+1)*q_ad(i, k)
          q_ad(i, k+1) = q_ad(i, k+1) - gam(i, k)*q_ad(i, k)
        END DO
      END DO
      d4_ad = 0.0
      DO i=i2,i1,-1
        a_bot = 1. + d4(i)*(d4(i)+1.5)
        CALL POPREALARRAY(q(i, km+1))
        temp2 = d4(i)*(d4(i)+0.5) - a_bot*gam(i, km)
        temp_ad11 = q_ad(i, km+1)/temp2
        temp1 = d4(i)*(d4(i)+1.)
        temp_ad12 = 2.*a4(1, i, km)*temp_ad11
        temp_ad13 = -((2.*(temp1*a4(1, i, km))+a4(1, i, km-1)-a_bot*q(i&
&         , km))*temp_ad11/temp2)
        a4_ad(1, i, km) = a4_ad(1, i, km) + 2.*temp1*temp_ad11
        a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + temp_ad11
        a_bot_ad = -(gam(i, km)*temp_ad13) - q(i, km)*temp_ad11
        d4_ad(i) = d4_ad(i) + (2*d4(i)+1.5)*a_bot_ad + (2*d4(i)+0.5)*&
&         temp_ad13 + (2*d4(i)+1.)*temp_ad12
        q_ad(i, km) = q_ad(i, km) - a_bot*temp_ad11
        gam_ad(i, km) = gam_ad(i, km) - a_bot*temp_ad13
        q_ad(i, km+1) = 0.0
      END DO
      DO k=km,2,-1
        DO i=i2,i1,-1
          temp_ad9 = q_ad(i, k)/bet
          temp_ad8 = 3.*temp_ad9
          CALL POPREALARRAY(q(i, k))
          bet_ad = -((3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))*&
&           temp_ad9/bet) - d4(i)*gam_ad(i, k)/bet**2
          d4_ad(i) = d4_ad(i) + a4(1, i, k)*temp_ad8 + 2*bet_ad + gam_ad&
&           (i, k)/bet
          gam_ad(i, k) = 0.0
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + temp_ad8
          a4_ad(1, i, k) = a4_ad(1, i, k) + d4(i)*temp_ad8
          q_ad(i, k-1) = q_ad(i, k-1) - temp_ad9
          q_ad(i, k) = 0.0
          CALL POPREALARRAY(bet)
          gam_ad(i, k-1) = gam_ad(i, k-1) - bet_ad
          CALL POPREALARRAY(d4(i))
          temp_ad10 = d4_ad(i)/delp(i, k)
          delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad10
          delp_ad(i, k) = delp_ad(i, k) - delp(i, k-1)*temp_ad10/delp(i&
&           , k)
          d4_ad(i) = 0.0
        END DO
      END DO
      DO i=i2,i1,-1
        grat = delp(i, 2)/delp(i, 1)
        bet = grat*(grat+0.5)
        temp_ad4 = gam_ad(i, 1)/bet
        gam_ad(i, 1) = 0.0
        temp_ad6 = q_ad(i, 1)/bet
        temp_ad5 = a4(1, i, 1)*temp_ad6
        temp0 = 2*grat*(grat+1.)
        bet_ad = -((temp0*a4(1, i, 1)+a4(1, i, 2))*temp_ad6/bet) - (grat&
&         *(grat+1.5)+1.)*temp_ad4/bet
        grat_ad = (4*grat+2*1.)*temp_ad5 + (2*grat+0.5)*bet_ad + (2*grat&
&         +1.5)*temp_ad4
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + temp0*temp_ad6
        a4_ad(1, i, 2) = a4_ad(1, i, 2) + temp_ad6
        q_ad(i, 1) = 0.0
        temp_ad7 = grat_ad/delp(i, 1)
        delp_ad(i, 2) = delp_ad(i, 2) + temp_ad7
        delp_ad(i, 1) = delp_ad(i, 1) - delp(i, 2)*temp_ad7/delp(i, 1)
      END DO
    ELSE
      gam_ad = 0.0
      DO k=1,km-1,1
        DO i=i2,i1,-1
          CALL POPREALARRAY(q(i, k))
          gam_ad(i, k+1) = gam_ad(i, k+1) - q(i, k+1)*q_ad(i, k)
          q_ad(i, k+1) = q_ad(i, k+1) - gam(i, k+1)*q_ad(i, k)
        END DO
      END DO
      DO i=i2,i1,-1
        CALL POPREALARRAY(q(i, km+1))
        qs_ad(i) = qs_ad(i) + q_ad(i, km+1)
        q_ad(i, km+1) = 0.0
        grat = delp(i, km-1)/delp(i, km)
        CALL POPREALARRAY(q(i, km))
        temp = 2*grat - gam(i, km) + 2.
        temp_ad1 = q_ad(i, km)/temp
        temp_ad2 = -((3.*(a4(1, i, km-1)+a4(1, i, km))-grat*qs(i)-q(i, &
&         km-1))*temp_ad1/temp)
        a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + 3.*temp_ad1
        a4_ad(1, i, km) = a4_ad(1, i, km) + 3.*temp_ad1
        grat_ad = 2*temp_ad2 - qs(i)*temp_ad1
        qs_ad(i) = qs_ad(i) - grat*temp_ad1
        q_ad(i, km-1) = q_ad(i, km-1) - temp_ad1
        gam_ad(i, km) = gam_ad(i, km) - temp_ad2
        q_ad(i, km) = 0.0
        temp_ad3 = grat_ad/delp(i, km)
        delp_ad(i, km-1) = delp_ad(i, km-1) + temp_ad3
        delp_ad(i, km) = delp_ad(i, km) - delp(i, km-1)*temp_ad3/delp(i&
&         , km)
      END DO
      DO k=km-1,2,-1
        DO i=i2,i1,-1
          temp_ad = q_ad(i, k)/bet
          CALL POPREALARRAY(q(i, k))
          grat = delp(i, k-1)/delp(i, k)
          bet_ad = -((3.*(a4(1, i, k-1)+a4(1, i, k))-q(i, k-1))*temp_ad/&
&           bet) - grat*gam_ad(i, k+1)/bet**2
          grat_ad = 2*bet_ad + gam_ad(i, k+1)/bet
          gam_ad(i, k+1) = 0.0
          a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + 3.*temp_ad
          a4_ad(1, i, k) = a4_ad(1, i, k) + 3.*temp_ad
          q_ad(i, k-1) = q_ad(i, k-1) - temp_ad
          q_ad(i, k) = 0.0
          CALL POPREALARRAY(bet)
          gam_ad(i, k) = gam_ad(i, k) - bet_ad
          temp_ad0 = grat_ad/delp(i, k)
          delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad0
          delp_ad(i, k) = delp_ad(i, k) - delp(i, k-1)*temp_ad0/delp(i, &
&           k)
        END DO
      END DO
      DO i=i2,i1,-1
        a4_ad(1, i, 1) = a4_ad(1, i, 1) + 1.5*q_ad(i, 1)
        q_ad(i, 1) = 0.0
      END DO
    END IF
  END SUBROUTINE CS_PROFILE_BWD
end module fv_mapz_adm_mod
