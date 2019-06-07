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
module fv_dynamics_adm_mod
   use constants_mod,           only: grav, pi=>pi_8, radius, hlv, rdgas, omega, rvgas, cp_vapor
   use dyn_core_adm_mod,        only: dyn_core, del2_cubed, init_ijk_mem
   use dyn_core_adm_mod,        only: dyn_core_fwd, dyn_core_bwd, del2_cubed_fwd, del2_cubed_bwd
   use fv_mapz_adm_mod,         only: compute_total_energy, Lagrangian_to_Eulerian, moist_cv, moist_cp
   use fv_mapz_adm_mod,         only: compute_total_energy_fwd, Lagrangian_to_Eulerian_fwd
   use fv_mapz_adm_mod,         only: compute_total_energy_bwd, Lagrangian_to_Eulerian_bwd
   use fv_tracer2d_adm_mod,     only: tracer_2d, tracer_2d_1L, tracer_2d_nested
   use fv_tracer2d_adm_mod,     only: tracer_2d_fwd, tracer_2d_1L_fwd, tracer_2d_nested_fwd
   use fv_tracer2d_adm_mod,     only: tracer_2d_bwd, tracer_2d_1L_bwd, tracer_2d_nested_bwd
   use fv_grid_utils_adm_mod,   only: cubed_to_latlon, c2l_ord2, g_sum
   use fv_grid_utils_adm_mod,   only: c2l_ord2_fwd
   use fv_grid_utils_adm_mod,   only: c2l_ord2_bwd, g_sum_adm
   use fv_fill_mod,             only: fill2D
   use fv_mp_mod,               only: is_master
   use fv_mp_mod,               only: group_halo_update_type
   use fv_mp_adm_mod,           only: start_group_halo_update, complete_group_halo_update
   use fv_mp_adm_mod,           only: start_group_halo_update_adm
   use fv_timing_mod,           only: timing_on, timing_off
   use diag_manager_mod,        only: send_data
   use fv_diagnostics_mod,      only: fv_time, prt_mxm, range_check, prt_minmax
   use mpp_domains_mod,         only: mpp_update_domains, DGRID_NE, CGRID_NE, domain2D
   use fv_mp_adm_mod,           only: mpp_update_domains_adm
   use mpp_mod,                 only: mpp_pe
   use field_manager_mod,       only: MODEL_ATMOS
   use tracer_manager_mod,      only: get_tracer_index
   use fv_sg_mod,               only: neg_adj3
   use fv_nesting_adm_mod,      only: setup_nested_grid_BCs
   use fv_nesting_adm_mod,      only: setup_nested_grid_BCs_adm
   use boundary_adm_mod,        only: nested_grid_BC_apply_intT
   use boundary_adm_mod,        only: nested_grid_BC_apply_intT_adm
   use fv_arrays_mod,           only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_grid_bounds_type
   use fv_arrays_mod,           only: R_GRID
   use fv_nwp_nudge_mod,        only: do_adiabatic_init
#ifdef MAPL_MODE
   use fv_control_mod,          only: dyn_timer, comm_timer
#endif
  use fv_arrays_mod,            only: fvprc

  use tapenade_iter,            only: pushcontrol, popcontrol, pushinteger, popinteger, &
                                      pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

  use fv_arrays_nlm_mod,        only: fv_flags_pert_type, fpp

implicit none

#ifdef MAPL_MODE
  ! Include the MPI library definitons:
  include 'mpif.h'
#endif

   logical :: RF_initialized = .false.
   logical :: pt_initialized = .false.
   logical :: bad_range = .false.
   real, allocatable ::  rf(:)
   integer :: kmax=1
   real :: agrav
   logical, public, save :: IdealTest=.false. 
#ifdef HIWPP
   real, allocatable:: u00(:,:,:), v00(:,:,:)
#endif
private
public :: fv_dynamics, fv_dynamics_fwd, fv_dynamics_bwd

!---- version number -----
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of fv_dynamics in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: peln q u v w delp ua delz va
!                pkz pe pt
!   with respect to varying inputs: peln q u v w delp delz pkz
!                pe pk pt
!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
  SUBROUTINE FV_DYNAMICS_FWD(npx, npy, npz, nq_tot, ng, bdt, consv_te, &
&   fill, reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split, &
&   q_split, u, v, w, delz, hydrostatic, pt, delp, q, ps, pe, pk, peln, &
&   pkz, phis, q_con, omga, ua, va, uc, vc, ak, bk, mfx, mfy, cx, cy, &
&   ze0, hybrid_z, gridstruct, flagstruct, flagstructp, neststruct, &
&   idiag, bd, parent_grid, domain, time_total)
    IMPLICIT NONE
! n_map loop
! Large time-step
    REAL, INTENT(IN) :: bdt
    REAL, INTENT(IN) :: consv_te
    REAL, INTENT(IN) :: kappa, cp_air
    REAL, INTENT(IN) :: zvir, ptop
    REAL, INTENT(IN), OPTIONAL :: time_total
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! transported tracers
    INTEGER, INTENT(IN) :: nq_tot
    INTEGER, INTENT(IN) :: ng
    INTEGER, INTENT(IN) :: ks
    INTEGER, INTENT(IN) :: ncnst
! small-step horizontal dynamics
    INTEGER, INTENT(IN) :: n_split
! tracer
    INTEGER, INTENT(IN) :: q_split
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
    LOGICAL, INTENT(IN) :: hydrostatic
! Using hybrid_z for remapping
    LOGICAL, INTENT(IN) :: hybrid_z
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
!  W (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! specific humidity and constituents
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! height at edges (m); non-hydrostatic
    REAL, INTENT(INOUT) :: ze0(bd%is:bd%is, bd%js:bd%js, 1)
! ze0 no longer used
!-----------------------------------------------------------------------
! Auxilliary pressure arrays:
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
! Surface pressure (pascal)
    REAL, INTENT(INOUT) :: ps(bd%isd:bd%ied, bd%jsd:bd%jed)
! edge pressure (pascal)
    REAL, INTENT(INOUT) :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
! pe**kappa
    REAL, INTENT(INOUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
! ln(pe)
    REAL, INTENT(INOUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
! finite-volume mean pk
    REAL, INTENT(INOUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: q_con(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! (uc,vc) mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua, va
    REAL, DIMENSION(npz+1), INTENT(IN) :: ak, bk
! Accumulated Mass flux arrays: the "Flux Capacitor"
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(INOUT) :: flagstruct
    TYPE(FV_FLAGS_PERT_TYPE), INTENT(INOUT) :: flagstructp
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(FV_DIAG_TYPE), INTENT(IN) :: idiag
! Local Arrays
    REAL :: ws(bd%is:bd%ie, bd%js:bd%je)
    REAL :: te_2d(bd%is:bd%ie, bd%js:bd%je)
    REAL :: teq(bd%is:bd%ie, bd%js:bd%je)
    REAL :: ps2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: m_fac(bd%is:bd%ie, bd%js:bd%je)
    REAL :: pfull(npz)
    REAL, DIMENSION(bd%is:bd%ie) :: cvm
    REAL :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz), dtdt_m(bd%is:bd%ie, &
&   bd%js:bd%je, npz), cappa(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
#ifdef OVERLOAD_R4
    REAL(kind=4) :: psx(bd%isd:bd%ied, bd%jsd:bd%jed)
#else
    REAL(kind=8) :: psx(bd%isd:bd%ied, bd%jsd:bd%jed)
#endif
    REAL(kind=8) :: dpx(bd%is:bd%ie, bd%js:bd%je)
    REAL :: akap, rdg, ph1, ph2, mdt, gam, amdt, u0
    INTEGER :: kord_tracer(ncnst), kord_mt, kord_wz, kord_tm
    INTEGER :: kord_tracer_pert(ncnst), kord_mt_pert, kord_wz_pert, &
&   kord_tm_pert
    INTEGER :: i, j, k, n, iq, n_map, nq, nwat, k_split
! GFDL physics
    INTEGER :: sphum
    INTEGER, SAVE :: liq_wat=-999
    INTEGER, SAVE :: ice_wat=-999
    INTEGER, SAVE :: rainwat=-999
    INTEGER, SAVE :: snowwat=-999
    INTEGER, SAVE :: graupel=-999
    INTEGER, SAVE :: cld_amt=-999
    INTEGER, SAVE :: theta_d=-999
    LOGICAL :: used, last_step, do_omega
    INTEGER, PARAMETER :: max_packs=12
    TYPE(GROUP_HALO_UPDATE_TYPE), SAVE :: i_pack(max_packs)
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL :: dt2
    REAL(kind=8) :: t1, t2
    INTEGER :: status
    REAL :: rf(npz)
    REAL :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pkc(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: ptc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: cry(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: divgd(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL :: delpc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: ut(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: zh(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pk3(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    INTRINSIC ANY
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC ABS
    INTRINSIC REAL
    INTRINSIC COS
    REAL :: abs0
    REAL :: abs1
    INTEGER :: arg1
    LOGICAL :: arg10
    REAL*8 :: arg11
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: arg12
    REAL :: result1
    LOGICAL :: res

    ws = 0.0
    te_2d = 0.0
    teq = 0.0
    ps2 = 0.0
    m_fac = 0.0
    pfull = 0.0
    cvm = 0.0
    dp1 = 0.0
    dtdt_m = 0.0
    cappa = 0.0
    akap = 0.0
    rdg = 0.0
    ph1 = 0.0
    ph2 = 0.0
    mdt = 0.0
    gam = 0.0
    amdt = 0.0
    u0 = 0.0
    kord_tracer = 0
    kord_mt = 0
    kord_wz = 0
    kord_tm = 0
    kord_tracer_pert = 0
    kord_mt_pert = 0
    kord_wz_pert = 0
    kord_tm_pert = 0
    iq = 0
    n_map = 0
    nq = 0
    nwat = 0
    k_split = 0
    sphum = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    dt2 = 0.0
    t1 = 0.0_8
    t2 = 0.0_8
    rf = 0.0
    gz = 0.0
    pkc = 0.0
    ptc = 0.0
    crx = 0.0
    xfx = 0.0
    cry = 0.0
    yfx = 0.0
    divgd = 0.0
    delpc = 0.0
    ut = 0.0
    vt = 0.0
    zh = 0.0
    pk3 = 0.0
    du = 0.0
    dv = 0.0
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
!     cv_air =  cp_air - rdgas
    agrav = 1./grav
    dt2 = 0.5*bdt
    k_split = flagstruct%k_split
    nwat = flagstruct%nwat
    nq = nq_tot - flagstruct%dnats
    rdg = -(rdgas*agrav)
!allocate ( dp1(isd:ied, jsd:jed, 1:npz) )
! Begin Dynamics timer for GEOS history processing
!t1 = MPI_Wtime(status)
!allocate ( cappa(isd:isd,jsd:jsd,1) )
!We call this BEFORE converting pt to virtual potential temperature,
!since we interpolate on (regular) temperature rather than theta.
    IF (gridstruct%nested .OR. ANY(neststruct%child_grids)) THEN
      CALL SETUP_NESTED_GRID_BCS(npx, npy, npz, zvir, ncnst, u, v, w, pt&
&                          , delp, delz, q, uc, vc, pkz, neststruct%&
&                          nested, flagstruct%inline_q, flagstruct%&
&                          make_nh, ng, gridstruct, flagstruct, &
&                          neststruct, neststruct%nest_timestep, &
&                          neststruct%tracer_nest_timestep, domain, bd, &
&                          nwat)
      IF (gridstruct%nested) THEN
!Correct halo values have now been set up for BCs; we can go ahead and apply them too...
        CALL NESTED_GRID_BC_APPLY_INTT(pt, 0, 0, npx, npy, npz, bd, 1., &
&                                1., neststruct%pt_bc, bctype=neststruct&
&                                %nestbctype)
        CALL PUSHCONTROL(2,0)
      ELSE
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHCONTROL(2,2)
    END IF
    IF (flagstruct%no_dycore) THEN
      IF (nwat .EQ. 2 .AND. (.NOT.hydrostatic)) THEN
        CALL PUSHCONTROL(1,0)
        sphum = GET_TRACER_INDEX(model_atmos, 'sphum')
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
!goto 911
    IF (fpp%fpp_mapl_mode) THEN
      SELECT CASE  (nwat) 
      CASE (0) 
        CALL PUSHCONTROL(1,0)
        sphum = 1
! to cause trouble if (mis)used
        cld_amt = -1
      CASE (1) 
        CALL PUSHCONTROL(1,0)
        sphum = 1
! to cause trouble if (mis)used
! to cause trouble if (mis)used
! to cause trouble if (mis)used
! to cause trouble if (mis)used
! to cause trouble if (mis)used
! to cause trouble if (mis)used
        cld_amt = -1
! to cause trouble if (mis)used
        theta_d = -1
      CASE (3) 
        CALL PUSHCONTROL(1,0)
        sphum = 1
! to cause trouble if (mis)used
! to cause trouble if (mis)used
! to cause trouble if (mis)used
! to cause trouble if (mis)used
        cld_amt = -1
! to cause trouble if (mis)used
        theta_d = -1
      CASE DEFAULT
        CALL PUSHCONTROL(1,0)
      END SELECT
    ELSE
      IF (nwat .EQ. 0) THEN
        CALL PUSHCONTROL(1,1)
        sphum = 1
! to cause trouble if (mis)used
        cld_amt = -1
      ELSE
        CALL PUSHCONTROL(1,1)
        sphum = GET_TRACER_INDEX(model_atmos, 'sphum')
        cld_amt = GET_TRACER_INDEX(model_atmos, 'cld_amt')
      END IF
      theta_d = GET_TRACER_INDEX(model_atmos, 'theta_d')
    END IF
    akap = kappa
!$OMP parallel do default(none) shared(npz,ak,bk,flagstruct,pfull) &
!$OMP                          private(ph1, ph2)
    DO k=1,npz
      ph1 = ak(k) + bk(k)*flagstruct%p_ref
      ph2 = ak(k+1) + bk(k+1)*flagstruct%p_ref
      pfull(k) = (ph2-ph1)/LOG(ph2/ph1)
    END DO
    IF (hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,zvir,nwat,q,q_con,sphum,liq_wat, &
!$OMP      rainwat,ice_wat,snowwat,graupel) private(cvm)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dp1(i, j, k) = zvir*q(i, j, k, sphum)
          END DO
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,zvir,q,q_con,sphum,liq_wat, &
!$OMP                                  rainwat,ice_wat,snowwat,graupel,pkz,flagstruct, &
!$OMP                                  cappa,kappa,rdg,delp,pt,delz,nwat)              &
!$OMP                          private(cvm)
      DO k=1,npz
        IF (flagstruct%moist_phys) THEN
          DO j=js,je
            DO i=is,ie
              dp1(i, j, k) = zvir*q(i, j, k, sphum)
              CALL PUSHREALARRAY(pkz(i, j, k))
              pkz(i, j, k) = EXP(kappa*LOG(rdg*delp(i, j, k)*pt(i, j, k)&
&               *(1.+dp1(i, j, k))/delz(i, j, k)))
            END DO
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
! Using dry pressure for the definition of the virtual potential temperature
!              pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!                                      (1.-q(i,j,k,sphum))/delz(i,j,k)) )
          DO j=js,je
            DO i=is,ie
              dp1(i, j, k) = 0.
              CALL PUSHREALARRAY(pkz(i, j, k))
              pkz(i, j, k) = EXP(kappa*LOG(rdg*delp(i, j, k)*pt(i, j, k)&
&               /delz(i, j, k)))
            END DO
          END DO
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHCONTROL(1,1)
    END IF
    IF (flagstruct%fv_debug) THEN
      IF (.NOT.hydrostatic) THEN
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
!---------------------
! Compute Total Energy
!---------------------
    IF (consv_te .GT. 0. .AND. (.NOT.do_adiabatic_init)) THEN
      CALL COMPUTE_TOTAL_ENERGY_FWD(is, ie, js, je, isd, ied, jsd, jed, &
&                             npz, u, v, w, delz, pt, delp, q, dp1, pe, &
&                             peln, phis, gridstruct%rsin2, gridstruct%&
&                             cosa_s, zvir, cp_air, rdgas, hlv, te_2d, &
&                             ua, va, teq, flagstruct%moist_phys, nwat, &
&                             sphum, liq_wat, rainwat, ice_wat, snowwat&
&                             , graupel, hydrostatic, idiag%id_te)
      IF (idiag%id_te .GT. 0) THEN
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF ((flagstruct%consv_am .OR. idiag%id_amdt .GT. 0) .AND. (.NOT.&
&       do_adiabatic_init)) THEN
      CALL COMPUTE_AAM_FWD(npz, is, ie, js, je, isd, ied, jsd, jed, &
&                    gridstruct, bd, ptop, ua, va, u, v, delp, teq, ps2&
&                    , m_fac)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (flagstruct%tau .GT. 0.) THEN
      IF (gridstruct%grid_type .LT. 4) THEN
        IF (bdt .GE. 0.) THEN
          abs0 = bdt
        ELSE
          abs0 = -bdt
        END IF
        arg10 = .NOT.neststruct%nested
        CALL RAYLEIGH_SUPER_FWD(abs0, npx, npy, npz, ks, pfull, phis, &
&                         flagstruct%tau, u, v, w, pt, ua, va, delz, &
&                         gridstruct%agrid, cp_air, rdgas, ptop, &
&                         hydrostatic, arg10, flagstruct%rf_cutoff, rf, &
&                         gridstruct, domain, bd)
        CALL PUSHCONTROL(2,0)
      ELSE
        IF (bdt .GE. 0.) THEN
          abs1 = bdt
        ELSE
          abs1 = -bdt
        END IF
        CALL RAYLEIGH_FRICTION_FWD(abs1, npx, npy, npz, ks, pfull, &
&                            flagstruct%tau, u, v, w, pt, ua, va, delz, &
&                            cp_air, rdgas, ptop, hydrostatic, .true., &
&                            flagstruct%rf_cutoff, rf, gridstruct, &
&                            domain, bd)
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHCONTROL(2,2)
    END IF
! Convert pt to virtual potential temperature on the first timestep
    IF (flagstruct%adiabatic) THEN
!$OMP parallel do default(none) shared(theta_d,is,ie,js,je,npz,pt,pkz,q)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(pt(i, j, k))
            pt(i, j, k) = pt(i, j, k)/pkz(i, j, k)
          END DO
        END DO
        IF (theta_d .GT. 0) THEN
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(q(i, j, k, theta_d))
              q(i, j, k, theta_d) = pt(i, j, k)
            END DO
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,dp1,pkz,q_con)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(pt(i, j, k))
            pt(i, j, k) = pt(i, j, k)*(1.+dp1(i, j, k))/pkz(i, j, k)
          END DO
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    END IF
    last_step = .false.
    mdt = bdt/REAL(k_split)
!DryMassRoundoffControl
!allocate(psx(isd:ied,jsd:jed),dpx(is:ie,js:je))
    IF (fpp%fpp_overload_r4) THEN
      DO j=js,je
        DO i=is,ie
          psx(i, j) = pe(i, npz+1, j)
          dpx(i, j) = 0.0
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
! first level of time-split
    DO n_map=1,k_split
      CALL PUSHREALARRAY(delp, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE(i_pack(1), delp, domain, complete=&
&                            .true.)
      CALL PUSHREALARRAY(pt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE(i_pack(1), pt, domain, complete=&
&                            .true.)
      CALL PUSHREALARRAY(v, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1)*npz)
      CALL PUSHREALARRAY(u, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2)*npz)
      CALL START_GROUP_HALO_UPDATE(i_pack(8), u, v, domain, gridtype=&
&                            dgrid_ne)
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,dp1,delp)
      DO k=1,npz
        DO j=jsd,jed
          DO i=isd,ied
            CALL PUSHREALARRAY(dp1(i, j, k))
            dp1(i, j, k) = delp(i, j, k)
          END DO
        END DO
      END DO
      IF (n_map .EQ. k_split) last_step = .true.
      CALL DYN_CORE_FWD(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir&
&                 , cp_air, akap, cappa, grav, hydrostatic, u, v, w, &
&                 delz, pt, q, delp, pe, pk, phis, ws, omga, ptop, pfull&
&                 , ua, va, uc, vc, mfx, mfy, cx, cy, pkz, peln, q_con, &
&                 ak, bk, dpx, ks, gridstruct, flagstruct, flagstructp, &
&                 neststruct, idiag, bd, domain, arg10, i_pack, &
&                 last_step, gz, pkc, ptc, crx, xfx, cry, yfx, divgd, &
&                 delpc, ut, vt, zh, pk3, du, dv, time_total)
!DryMassRoundoffControl
      IF (last_step) THEN
        IF (fpp%fpp_overload_r4) THEN
          DO j=js,je
            DO i=is,ie
              psx(i, j) = psx(i, j) + dpx(i, j)
            END DO
          END DO
          CALL MPP_UPDATE_DOMAINS(psx, domain)
          DO j=js-1,je+1
            DO i=is-1,ie+1
              CALL PUSHREALARRAY(pe(i, npz+1, j))
              pe(i, npz+1, j) = psx(i, j)
            END DO
          END DO
          CALL PUSHCONTROL(2,0)
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,2)
      END IF
!deallocate(psx,dpx)
      IF (.NOT.flagstruct%inline_q .AND. nq .NE. 0) THEN
!--------------------------------------------------------
! Perform large-time-step scalar transport using the accumulated CFL and
! mass fluxes
!!! CLEANUP: merge these two calls?
        IF (gridstruct%nested) THEN
          CALL TRACER_2D_NESTED_FWD(q, dp1, mfx, mfy, cx, cy, gridstruct&
&                             , bd, domain, npx, npy, npz, nq, &
&                             flagstruct%hord_tr, q_split, mdt, idiag%&
&                             id_divg, i_pack(10), flagstruct%nord_tr, &
&                             flagstruct%trdm2, k_split, neststruct, &
&                             parent_grid, flagstructp%hord_tr_pert, &
&                             flagstructp%nord_tr_pert, flagstructp%&
&                             trdm2_pert, flagstructp%split_damp_tr)
          CALL PUSHCONTROL(2,0)
        ELSE IF (flagstruct%z_tracer) THEN
          CALL TRACER_2D_1L_FWD(q, dp1, mfx, mfy, cx, cy, gridstruct, bd&
&                         , domain, npx, npy, npz, nq, flagstruct%&
&                         hord_tr, q_split, mdt, idiag%id_divg, i_pack(&
&                         10), flagstruct%nord_tr, flagstruct%trdm2, &
&                         flagstructp%hord_tr_pert, flagstructp%&
&                         nord_tr_pert, flagstructp%trdm2_pert, &
&                         flagstructp%split_damp_tr)
          CALL PUSHCONTROL(2,1)
        ELSE
          CALL TRACER_2D_FWD(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, &
&                      domain, npx, npy, npz, nq, flagstruct%hord_tr, &
&                      q_split, mdt, idiag%id_divg, i_pack(10), &
&                      flagstruct%nord_tr, flagstruct%trdm2, flagstructp&
&                      %hord_tr_pert, flagstructp%nord_tr_pert, &
&                      flagstructp%trdm2_pert, flagstructp%split_damp_tr&
&                     )
          CALL PUSHCONTROL(2,2)
        END IF
      ELSE
        CALL PUSHCONTROL(2,3)
      END IF
      IF (npz .GT. 4) THEN
!------------------------------------------------------------------------
! Perform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.
!------------------------------------------------------------------------
        DO iq=1,nq
          kord_tracer(iq) = flagstruct%kord_tr
! monotonic
          IF (iq .EQ. cld_amt) kord_tracer(iq) = 9
          CALL PUSHINTEGER(kord_tracer_pert(iq))
          kord_tracer_pert(iq) = flagstructp%kord_tr_pert
! linear
          IF (iq .EQ. cld_amt) THEN
            CALL PUSHCONTROL(1,1)
            kord_tracer_pert(iq) = 17
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        do_omega = hydrostatic .AND. last_step
        kord_mt = flagstruct%kord_mt
        kord_wz = flagstruct%kord_wz
        kord_tm = flagstruct%kord_tm
        kord_mt_pert = flagstructp%kord_mt_pert
        kord_wz_pert = flagstructp%kord_wz_pert
        kord_tm_pert = flagstructp%kord_tm_pert
        IF (n_map .EQ. k_split) THEN
          kord_mt = kord_mt_pert
          kord_wz = kord_wz_pert
          kord_tm = kord_tm_pert
          kord_tracer = kord_tracer_pert
        END IF
        CALL LAGRANGIAN_TO_EULERIAN_FWD(last_step, consv_te, ps, pe, &
&                                 delp, pkz, pk, mdt, bdt, npz, is, ie, &
&                                 js, je, isd, ied, jsd, jed, nq, nwat, &
&                                 sphum, q_con, u, v, w, delz, pt, q, &
&                                 phis, zvir, cp_air, akap, cappa, &
&                                 kord_mt, kord_wz, kord_tracer, kord_tm&
&                                 , peln, te_2d, ng, ua, va, omga, dp1, &
&                                 ws, fill, reproduce_sum, arg10, dtdt_m&
&                                 , ptop, ak, bk, pfull, flagstruct, &
&                                 gridstruct, domain, flagstruct%&
&                                 do_sat_adj, hydrostatic, hybrid_z, &
&                                 do_omega, flagstruct%adiabatic, &
&                                 do_adiabatic_init, mfx, mfy, &
&                                 flagstruct%remap_option, kord_mt_pert&
&                                 , kord_wz_pert, kord_tracer_pert, &
&                                 kord_tm_pert)
        IF (last_step) THEN
          IF (.NOT.hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga,delp,delz,w)
            DO k=1,npz
              DO j=js,je
                DO i=is,ie
                  CALL PUSHREALARRAY(omga(i, j, k))
                  omga(i, j, k) = delp(i, j, k)/delz(i, j, k)*w(i, j, k)
                END DO
              END DO
            END DO
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
!--------------------------
! Filter omega for physics:
!--------------------------
          IF (flagstruct%nf_omega .GT. 0) THEN
            arg11 = 0.18*gridstruct%da_min
            CALL DEL2_CUBED_FWD(omga, arg11, gridstruct, domain, npx, &
&                         npy, npz, flagstruct%nf_omega, bd)
            CALL PUSHCONTROL(2,3)
          ELSE
            CALL PUSHCONTROL(2,2)
          END IF
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,0)
      END IF
    END DO
    IF (nwat .EQ. 6) THEN
      IF (flagstruct%fv_debug) THEN
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (((flagstruct%consv_am .OR. idiag%id_amdt .GT. 0) .OR. idiag%&
&       id_aam .GT. 0) .AND. (.NOT.do_adiabatic_init)) THEN
      CALL COMPUTE_AAM_FWD(npz, is, ie, js, je, isd, ied, jsd, jed, &
&                    gridstruct, bd, ptop, ua, va, u, v, delp, te_2d, ps&
&                    , m_fac)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF ((flagstruct%consv_am .OR. idiag%id_amdt .GT. 0) .AND. (.NOT.&
&       do_adiabatic_init)) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,te_2d,teq,dt2,ps2,ps,idiag)
      DO j=js,je
        DO i=is,ie
! Note: the mountain torque computation contains also numerical error
! The numerical error is mostly from the zonal gradient of the terrain (zxg)
          te_2d(i, j) = te_2d(i, j) - teq(i, j) + dt2*(ps2(i, j)+ps(i, j&
&           ))*idiag%zxg(i, j)
        END DO
      END DO
      IF (flagstruct%consv_am .OR. prt_minmax) THEN
        amdt = G_SUM(domain, te_2d, is, ie, js, je, ng, gridstruct%&
&         area_64, 0, reproduce=.true.)
        result1 = G_SUM(domain, m_fac, is, ie, js, je, ng, gridstruct%&
&         area_64, 0, reproduce=.true.)
        u0 = -(radius*amdt/result1)
        res = IS_MASTER()
        IF (res .AND. prt_minmax) THEN
          CALL PUSHCONTROL(1,0)
          WRITE(6, *) 'Dynamic AM tendency (Hadleys)=', amdt/(bdt*1.e18)&
&         , 'del-u (per day)=', u0*86400./bdt
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
!  consv_am
      IF (flagstruct%consv_am) THEN
        CALL PUSHCONTROL(2,0)
      ELSE
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHCONTROL(2,2)
    END IF
    CALL PUSHINTEGER(jed)
    CALL PUSHREALARRAY(cry, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    CALL PUSHREALARRAY(crx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
    CALL PUSHREALARRAY(amdt)
    CALL PUSHINTEGER(sphum)
    CALL PUSHREALARRAY(result1)
    CALL PUSHREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL PUSHREALARRAY(te_2d, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL PUSHINTEGER(isd)
    CALL PUSHREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    CALL PUSHREALARRAY(m_fac, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL PUSHINTEGER(kord_tracer_pert, ncnst)
    CALL PUSHREALARRAY(delpc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL PUSHREALARRAY(pkc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(npz+1)&
&                )
    CALL PUSHREALARRAY(ut, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL PUSHINTEGER(ied)
    CALL PUSHREALARRAY(rf, npz)
    CALL PUSHINTEGER(jsd)
    CALL PUSHREALARRAY(pfull, npz)
    CALL PUSHREALARRAY(ptc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL PUSHREALARRAY(dp1, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL PUSHREALARRAY(gz, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(npz+1))
    CALL PUSHREALARRAY(ws, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL PUSHREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
    CALL PUSHREALARRAY(pk3, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(npz+1)&
&                )
  END SUBROUTINE FV_DYNAMICS_FWD
!  Differentiation of fv_dynamics in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: peln q u v w delp ua delz va
!                pkz pe pt
!   with respect to varying inputs: peln q u v w delp delz pkz
!                pe pk pt
!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
  SUBROUTINE FV_DYNAMICS_BWD(npx, npy, npz, nq_tot, ng, bdt, consv_te, &
&   fill, reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split, &
&   q_split, u, u_ad, v, v_ad, w, w_ad, delz, delz_ad, hydrostatic, pt, &
&   pt_ad, delp, delp_ad, q, q_ad, ps, ps_ad, pe, pe_ad, pk, pk_ad, peln&
&   , peln_ad, pkz, pkz_ad, phis, q_con, omga, omga_ad, ua, ua_ad, va, &
&   va_ad, uc, uc_ad, vc, vc_ad, ak, bk, mfx, mfx_ad, mfy, mfy_ad, cx, &
&   cx_ad, cy, cy_ad, ze0, hybrid_z, gridstruct, flagstruct, flagstructp&
&   , neststruct, idiag, bd, parent_grid, domain, time_total)
    IMPLICIT NONE
    REAL, INTENT(IN) :: bdt
    REAL, INTENT(IN) :: consv_te
    REAL, INTENT(IN) :: kappa, cp_air
    REAL, INTENT(IN) :: zvir, ptop
    REAL, INTENT(IN), OPTIONAL :: time_total
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: nq_tot
    INTEGER, INTENT(IN) :: ng
    INTEGER, INTENT(IN) :: ks
    INTEGER, INTENT(IN) :: ncnst
    INTEGER, INTENT(IN) :: n_split
    INTEGER, INTENT(IN) :: q_split
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: hybrid_z
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v_ad
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: w_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst&
&   )
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delz_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ze0(bd%is:bd%is, bd%js:bd%js, 1)
    REAL, INTENT(INOUT) :: ps(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: ps_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(INOUT) :: pe_ad(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1&
&   )
    REAL, INTENT(INOUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL, INTENT(INOUT) :: pk_ad(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL, INTENT(INOUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: peln_ad(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: pkz_ad(bd%is:bd%ie, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: q_con(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: omga(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: omga_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: uc_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: vc_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua_ad, va_ad
    REAL, DIMENSION(npz+1), INTENT(IN) :: ak, bk
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_ad(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_ad(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(INOUT) :: flagstruct
    TYPE(FV_FLAGS_PERT_TYPE), INTENT(INOUT) :: flagstructp
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(FV_DIAG_TYPE), INTENT(IN) :: idiag
    REAL :: ws(bd%is:bd%ie, bd%js:bd%je)
    REAL :: ws_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL :: te_2d(bd%is:bd%ie, bd%js:bd%je)
    REAL :: te_2d_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL :: teq(bd%is:bd%ie, bd%js:bd%je)
    REAL :: teq_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL :: ps2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: ps2_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: m_fac(bd%is:bd%ie, bd%js:bd%je)
    REAL :: m_fac_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL :: pfull(npz)
    REAL, DIMENSION(bd%is:bd%ie) :: cvm
    REAL :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz), dtdt_m(bd%is:bd%ie, &
&   bd%js:bd%je, npz), cappa(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: dp1_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
#ifdef OVERLOAD_R4
    REAL(kind=4) :: psx(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL(kind=4) :: psx_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
#else
    REAL(kind=8) :: psx(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL(kind=8) :: psx_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
#endif
    REAL(kind=8) :: dpx(bd%is:bd%ie, bd%js:bd%je)
    REAL(kind=8) :: dpx_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL :: akap, rdg, ph1, ph2, mdt, gam, amdt, u0
    REAL :: amdt_ad, u0_ad
    INTEGER :: kord_tracer(ncnst), kord_mt, kord_wz, kord_tm
    INTEGER :: kord_tracer_pert(ncnst), kord_mt_pert, kord_wz_pert, &
&   kord_tm_pert
    INTEGER :: i, j, k, n, iq, n_map, nq, nwat, k_split
    INTEGER :: sphum
    INTEGER, SAVE :: liq_wat=-999
    INTEGER, SAVE :: ice_wat=-999
    INTEGER, SAVE :: rainwat=-999
    INTEGER, SAVE :: snowwat=-999
    INTEGER, SAVE :: graupel=-999
    INTEGER, SAVE :: cld_amt=-999
    INTEGER, SAVE :: theta_d=-999
    LOGICAL :: used, last_step, do_omega
    INTEGER, PARAMETER :: max_packs=12
    TYPE(GROUP_HALO_UPDATE_TYPE), SAVE :: i_pack(max_packs)
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL :: dt2
    REAL(kind=8) :: t1, t2
    INTEGER :: status
    REAL :: rf(npz)
    REAL :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: gz_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pkc(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pkc_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: ptc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: ptc_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: crx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: cry(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cry_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: divgd(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL :: divgd_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL :: delpc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: delpc_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: ut(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: ut_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: vt_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: zh(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: zh_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pk3(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pk3_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL :: du_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL :: dv_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    INTRINSIC ANY
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC ABS
    INTRINSIC REAL
    INTRINSIC COS
    REAL :: abs0
    REAL :: abs1
    INTEGER :: arg1
    LOGICAL :: arg10
    REAL*8 :: arg11
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: arg12
    REAL :: result1
    REAL :: result1_ad
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp2
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp3
    REAL :: temp4
    REAL :: temp5
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp6
    REAL :: temp7
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    INTEGER :: branch

    ws = 0.0
    te_2d = 0.0
    teq = 0.0
    ps2 = 0.0
    m_fac = 0.0
    pfull = 0.0
    cvm = 0.0
    dp1 = 0.0
    dtdt_m = 0.0
    cappa = 0.0
    akap = 0.0
    rdg = 0.0
    ph1 = 0.0
    ph2 = 0.0
    mdt = 0.0
    gam = 0.0
    amdt = 0.0
    u0 = 0.0
    kord_tracer = 0
    kord_mt = 0
    kord_wz = 0
    kord_tm = 0
    kord_tracer_pert = 0
    kord_mt_pert = 0
    kord_wz_pert = 0
    kord_tm_pert = 0
    iq = 0
    n_map = 0
    nq = 0
    nwat = 0
    k_split = 0
    sphum = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    dt2 = 0.0
    t1 = 0.0_8
    t2 = 0.0_8
    rf = 0.0
    gz = 0.0
    pkc = 0.0
    ptc = 0.0
    crx = 0.0
    xfx = 0.0
    cry = 0.0
    yfx = 0.0
    divgd = 0.0
    delpc = 0.0
    ut = 0.0
    vt = 0.0
    zh = 0.0
    pk3 = 0.0
    du = 0.0
    dv = 0.0
    abs0 = 0.0
    abs1 = 0.0
    arg1 = 0
    arg11 = 0.0_r_grid
    arg12 = 0.0
    result1 = 0.0
    branch = 0

    CALL POPREALARRAY(pk3, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(npz+1))
    CALL POPREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
    CALL POPREALARRAY(ws, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL POPREALARRAY(gz, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(npz+1))
    CALL POPREALARRAY(dp1, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL POPREALARRAY(ptc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL POPREALARRAY(pfull, npz)
    CALL POPINTEGER(jsd)
    CALL POPREALARRAY(rf, npz)
    CALL POPINTEGER(ied)
    CALL POPREALARRAY(ut, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL POPREALARRAY(pkc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(npz+1))
    CALL POPREALARRAY(delpc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL POPINTEGER(kord_tracer_pert, ncnst)
    CALL POPREALARRAY(m_fac, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL POPREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    CALL POPINTEGER(isd)
    CALL POPREALARRAY(te_2d, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL POPREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL POPREALARRAY(result1)
    CALL POPINTEGER(sphum)
    CALL POPREALARRAY(amdt)
    CALL POPREALARRAY(crx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
    CALL POPREALARRAY(cry, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    CALL POPINTEGER(jed)
    js = bd%js
    ie = bd%ie
    is = bd%is
    je = bd%je
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      u0_ad = 0.0
      DO k=npz,1,-1
        DO j=je,js,-1
          DO i=ie+1,is,-1
            u0_ad = u0_ad + gridstruct%l2c_v(i, j)*v_ad(i, j, k)
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ie,is,-1
            u0_ad = u0_ad + gridstruct%l2c_u(i, j)*u_ad(i, j, k)
          END DO
        END DO
      END DO
    ELSE IF (branch .EQ. 1) THEN
      u0_ad = 0.0
    ELSE
      ps_ad = 0.0
      teq_ad = 0.0
      ps2_ad = 0.0
      m_fac_ad = 0.0
      te_2d_ad = 0.0
      GOTO 100
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      temp_ad6 = -(radius*u0_ad/result1)
      amdt_ad = temp_ad6
      result1_ad = -(amdt*temp_ad6/result1)
      CALL G_SUM_ADM(domain, m_fac, m_fac_ad, is, ie, js, je, ng, &
&              gridstruct%area_64, 0, reproduce=.true., g_sum_ad=&
&              result1_ad)
      CALL G_SUM_ADM(domain, te_2d, te_2d_ad, is, ie, js, je, ng, &
&              gridstruct%area_64, 0, reproduce=.true., g_sum_ad=amdt_ad&
&             )
    ELSE
      m_fac_ad = 0.0
      te_2d_ad = 0.0
    END IF
    dt2 = 0.5*bdt
    ps_ad = 0.0
    teq_ad = 0.0
    ps2_ad = 0.0
    DO j=je,js,-1
      DO i=ie,is,-1
        temp_ad5 = dt2*idiag%zxg(i, j)*te_2d_ad(i, j)
        teq_ad(i, j) = teq_ad(i, j) - te_2d_ad(i, j)
        ps2_ad(i, j) = ps2_ad(i, j) + temp_ad5
        ps_ad(i, j) = ps_ad(i, j) + temp_ad5
      END DO
    END DO
 100 CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      jsd = bd%jsd
      ied = bd%ied
      isd = bd%isd
      jed = bd%jed
      CALL COMPUTE_AAM_BWD(npz, is, ie, js, je, isd, ied, jsd, jed, &
&                    gridstruct, bd, ptop, ua, ua_ad, va, va_ad, u, u_ad&
&                    , v, v_ad, delp, delp_ad, te_2d, te_2d_ad, ps, &
&                    ps_ad, m_fac, m_fac_ad)
    END IF
    CALL POPCONTROL(1,branch)
    nq = nq_tot - flagstruct%dnats
    k_split = flagstruct%k_split
    akap = kappa
    uc_ad = 0.0
    omga_ad = 0.0
    vc_ad = 0.0
    pk_ad = 0.0
    pk3_ad = 0.0
    xfx_ad = 0.0
    ws_ad = 0.0
    gz_ad = 0.0
    du_ad = 0.0
    psx_ad = 0.0_8
    dv_ad = 0.0
    dp1_ad = 0.0
    ptc_ad = 0.0
    ut_ad = 0.0
    divgd_ad = 0.0
    pkc_ad = 0.0
    delpc_ad = 0.0
    yfx_ad = 0.0
    vt_ad = 0.0
    zh_ad = 0.0
    dpx_ad = 0.0_8
    crx_ad = 0.0
    cry_ad = 0.0
    DO n_map=k_split,1,-1
      CALL POPCONTROL(2,branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) GOTO 110
      ELSE
        IF (branch .NE. 2) THEN
          arg11 = 0.18*gridstruct%da_min
          CALL DEL2_CUBED_BWD(omga, omga_ad, arg11, gridstruct, domain, &
&                       npx, npy, npz, flagstruct%nf_omega, bd)
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO k=npz,1,-1
            DO j=je,js,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(omga(i, j, k))
                temp_ad4 = omga_ad(i, j, k)/delz(i, j, k)
                delp_ad(i, j, k) = delp_ad(i, j, k) + w(i, j, k)*&
&                 temp_ad4
                w_ad(i, j, k) = w_ad(i, j, k) + delp(i, j, k)*temp_ad4
                delz_ad(i, j, k) = delz_ad(i, j, k) - delp(i, j, k)*w(i&
&                 , j, k)*temp_ad4/delz(i, j, k)
                omga_ad(i, j, k) = 0.0
              END DO
            END DO
          END DO
        END IF
      END IF
      kord_mt_pert = flagstructp%kord_mt_pert
      kord_wz_pert = flagstructp%kord_wz_pert
      CALL LAGRANGIAN_TO_EULERIAN_BWD(last_step, consv_te, ps, ps_ad, pe&
&                               , pe_ad, delp, delp_ad, pkz, pkz_ad, pk&
&                               , pk_ad, mdt, bdt, npz, is, ie, js, je, &
&                               isd, ied, jsd, jed, nq, nwat, sphum, &
&                               q_con, u, u_ad, v, v_ad, w, w_ad, delz, &
&                               delz_ad, pt, pt_ad, q, q_ad, phis, zvir&
&                               , cp_air, akap, cappa, kord_mt, kord_wz&
&                               , kord_tracer, kord_tm, peln, peln_ad, &
&                               te_2d, te_2d_ad, ng, ua, ua_ad, va, omga&
&                               , omga_ad, dp1, dp1_ad, ws, ws_ad, fill&
&                               , reproduce_sum, arg10, dtdt_m, ptop, ak&
&                               , bk, pfull, flagstruct, gridstruct, &
&                               domain, flagstruct%do_sat_adj, &
&                               hydrostatic, hybrid_z, do_omega, &
&                               flagstruct%adiabatic, do_adiabatic_init&
&                               , mfx, mfy, flagstruct%remap_option, &
&                               kord_mt_pert, kord_wz_pert, &
&                               kord_tracer_pert, kord_tm_pert)
      DO iq=nq,1,-1
        CALL POPCONTROL(1,branch)
        CALL POPINTEGER(kord_tracer_pert(iq))
      END DO
 110  CALL POPCONTROL(2,branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          CALL TRACER_2D_NESTED_BWD(q, q_ad, dp1, dp1_ad, mfx, mfx_ad, &
&                             mfy, mfy_ad, cx, cx_ad, cy, cy_ad, &
&                             gridstruct, bd, domain, npx, npy, npz, nq&
&                             , flagstruct%hord_tr, q_split, mdt, idiag%&
&                             id_divg, i_pack(10), flagstruct%nord_tr, &
&                             flagstruct%trdm2, k_split, neststruct, &
&                             parent_grid, flagstructp%hord_tr_pert, &
&                             flagstructp%nord_tr_pert, flagstructp%&
&                             trdm2_pert, flagstructp%split_damp_tr)
        ELSE
          CALL TRACER_2D_1L_BWD(q, q_ad, dp1, dp1_ad, mfx, mfx_ad, mfy, &
&                         mfy_ad, cx, cx_ad, cy, cy_ad, gridstruct, bd, &
&                         domain, npx, npy, npz, nq, flagstruct%hord_tr&
&                         , q_split, mdt, idiag%id_divg, i_pack(10), &
&                         flagstruct%nord_tr, flagstruct%trdm2, &
&                         flagstructp%hord_tr_pert, flagstructp%&
&                         nord_tr_pert, flagstructp%trdm2_pert, &
&                         flagstructp%split_damp_tr)
        END IF
      ELSE IF (branch .EQ. 2) THEN
        CALL TRACER_2D_BWD(q, q_ad, dp1, dp1_ad, mfx, mfx_ad, mfy, &
&                    mfy_ad, cx, cx_ad, cy, cy_ad, gridstruct, bd, &
&                    domain, npx, npy, npz, nq, flagstruct%hord_tr, &
&                    q_split, mdt, idiag%id_divg, i_pack(10), flagstruct&
&                    %nord_tr, flagstruct%trdm2, flagstructp%&
&                    hord_tr_pert, flagstructp%nord_tr_pert, flagstructp&
&                    %trdm2_pert, flagstructp%split_damp_tr)
      ELSE
        mfx_ad = 0.0
        mfy_ad = 0.0
        cx_ad = 0.0
        cy_ad = 0.0
      END IF
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        DO j=je+1,js-1,-1
          DO i=ie+1,is-1,-1
            CALL POPREALARRAY(pe(i, npz+1, j))
            psx_ad(i, j) = psx_ad(i, j) + pe_ad(i, npz+1, j)
            pe_ad(i, npz+1, j) = 0.0
          END DO
        END DO
        CALL MPP_UPDATE_DOMAINS_ADM(psx, psx_ad, domain)
        DO j=je,js,-1
          DO i=ie,is,-1
            dpx_ad(i, j) = dpx_ad(i, j) + psx_ad(i, j)
          END DO
        END DO
      END IF
      CALL DYN_CORE_BWD(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir&
&                 , cp_air, akap, cappa, grav, hydrostatic, u, u_ad, v, &
&                 v_ad, w, w_ad, delz, delz_ad, pt, pt_ad, q, q_ad, delp&
&                 , delp_ad, pe, pe_ad, pk, pk_ad, phis, ws, ws_ad, omga&
&                 , omga_ad, ptop, pfull, ua, ua_ad, va, va_ad, uc, &
&                 uc_ad, vc, vc_ad, mfx, mfx_ad, mfy, mfy_ad, cx, cx_ad&
&                 , cy, cy_ad, pkz, pkz_ad, peln, peln_ad, q_con, ak, bk&
&                 , dpx, dpx_ad, ks, gridstruct, flagstruct, flagstructp&
&                 , neststruct, idiag, bd, domain, arg10, i_pack, &
&                 last_step, gz, gz_ad, pkc, pkc_ad, ptc, ptc_ad, crx, &
&                 crx_ad, xfx, xfx_ad, cry, cry_ad, yfx, yfx_ad, divgd, &
&                 divgd_ad, delpc, delpc_ad, ut, ut_ad, vt, vt_ad, zh, &
&                 zh_ad, pk3, pk3_ad, du, du_ad, dv, dv_ad, time_total)
      DO k=npz,1,-1
        DO j=jed,jsd,-1
          DO i=ied,isd,-1
            CALL POPREALARRAY(dp1(i, j, k))
            delp_ad(i, j, k) = delp_ad(i, j, k) + dp1_ad(i, j, k)
            dp1_ad(i, j, k) = 0.0
          END DO
        END DO
      END DO
      CALL POPREALARRAY(u, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2)*npz)
      CALL POPREALARRAY(v, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE_ADM(i_pack(8), u, u_ad, v, v_ad, &
&                                domain, gridtype=dgrid_ne)
      CALL POPREALARRAY(pt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE_ADM(i_pack(1), pt, pt_ad, domain, &
&                                complete=.true.)
      CALL POPREALARRAY(delp, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE_ADM(i_pack(1), delp, delp_ad, domain&
&                                , complete=.true.)
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        DO i=ie,is,-1
          pe_ad(i, npz+1, j) = pe_ad(i, npz+1, j) + psx_ad(i, j)
          psx_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO k=npz,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(q(i, j, k, theta_d))
              pt_ad(i, j, k) = pt_ad(i, j, k) + q_ad(i, j, k, theta_d)
              q_ad(i, j, k, theta_d) = 0.0
            END DO
          END DO
        END IF
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(pt(i, j, k))
            temp_ad2 = pt_ad(i, j, k)/pkz(i, j, k)
            pkz_ad(i, j, k) = pkz_ad(i, j, k) - pt(i, j, k)*temp_ad2/pkz&
&             (i, j, k)
            pt_ad(i, j, k) = temp_ad2
          END DO
        END DO
      END DO
    ELSE
      DO k=npz,1,-1
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(pt(i, j, k))
            temp7 = pkz(i, j, k)
            temp6 = pt(i, j, k)/temp7
            temp_ad3 = (dp1(i, j, k)+1.)*pt_ad(i, j, k)/temp7
            dp1_ad(i, j, k) = dp1_ad(i, j, k) + temp6*pt_ad(i, j, k)
            pkz_ad(i, j, k) = pkz_ad(i, j, k) - temp6*temp_ad3
            pt_ad(i, j, k) = temp_ad3
          END DO
        END DO
      END DO
    END IF
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      CALL RAYLEIGH_SUPER_BWD(abs0, npx, npy, npz, ks, pfull, phis, &
&                       flagstruct%tau, u, u_ad, v, v_ad, w, w_ad, pt, &
&                       pt_ad, ua, ua_ad, va, va_ad, delz, gridstruct%&
&                       agrid, cp_air, rdgas, ptop, hydrostatic, arg10, &
&                       flagstruct%rf_cutoff, rf, gridstruct, domain, bd&
&                      )
    ELSE IF (branch .EQ. 1) THEN
      CALL RAYLEIGH_FRICTION_BWD(abs1, npx, npy, npz, ks, pfull, &
&                          flagstruct%tau, u, u_ad, v, v_ad, w, w_ad, pt&
&                          , pt_ad, ua, ua_ad, va, va_ad, delz, delz_ad&
&                          , cp_air, rdgas, ptop, hydrostatic, .true., &
&                          flagstruct%rf_cutoff, rf, gridstruct, domain&
&                          , bd)
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) CALL COMPUTE_AAM_BWD(npz, is, ie, js, je, isd, &
&                                     ied, jsd, jed, gridstruct, bd, &
&                                     ptop, ua, ua_ad, va, va_ad, u, &
&                                     u_ad, v, v_ad, delp, delp_ad, teq&
&                                     , teq_ad, ps2, ps2_ad, m_fac, &
&                                     m_fac_ad)
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) CALL COMPUTE_TOTAL_ENERGY_BWD(is, ie, js, je, isd&
&                                              , ied, jsd, jed, npz, u, &
&                                              u_ad, v, v_ad, w, w_ad, &
&                                              delz, delz_ad, pt, pt_ad&
&                                              , delp, delp_ad, q, q_ad&
&                                              , dp1, dp1_ad, pe, pe_ad&
&                                              , peln, peln_ad, phis, &
&                                              gridstruct%rsin2, &
&                                              gridstruct%cosa_s, zvir, &
&                                              cp_air, rdgas, hlv, te_2d&
&                                              , te_2d_ad, ua, va, teq, &
&                                              teq_ad, flagstruct%&
&                                              moist_phys, nwat, sphum, &
&                                              liq_wat, rainwat, ice_wat&
&                                              , snowwat, graupel, &
&                                              hydrostatic, idiag%id_te)
    CALL POPCONTROL(1,branch)
    rdg = -(rdgas*agrav)
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO k=npz,1,-1
        DO j=je,js,-1
          DO i=ie,is,-1
            q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) + zvir*dp1_ad(i&
&             , j, k)
            dp1_ad(i, j, k) = 0.0
          END DO
        END DO
      END DO
    ELSE
      DO k=npz,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(pkz(i, j, k))
              temp5 = delz(i, j, k)
              temp4 = delp(i, j, k)*pt(i, j, k)
              temp3 = temp4/temp5
              temp_ad1 = kappa*EXP(kappa*LOG(rdg*temp3))*pkz_ad(i, j, k)&
&               /(temp3*temp5)
              delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*temp_ad1
              pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*temp_ad1
              delz_ad(i, j, k) = delz_ad(i, j, k) - temp3*temp_ad1
              pkz_ad(i, j, k) = 0.0
              dp1_ad(i, j, k) = 0.0
            END DO
          END DO
        ELSE
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(pkz(i, j, k))
              temp2 = delz(i, j, k)
              temp = (dp1(i, j, k)+1.)/temp2
              temp1 = delp(i, j, k)*pt(i, j, k)
              temp0 = rdg*temp1*temp
              temp_ad = kappa*EXP(kappa*LOG(temp0))*rdg*pkz_ad(i, j, k)/&
&               temp0
              temp_ad0 = temp1*temp_ad/temp2
              delp_ad(i, j, k) = delp_ad(i, j, k) + temp*pt(i, j, k)*&
&               temp_ad
              pt_ad(i, j, k) = pt_ad(i, j, k) + temp*delp(i, j, k)*&
&               temp_ad
              dp1_ad(i, j, k) = dp1_ad(i, j, k) + temp_ad0
              delz_ad(i, j, k) = delz_ad(i, j, k) - temp*temp_ad0
              pkz_ad(i, j, k) = 0.0
              q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) + zvir*dp1_ad(&
&               i, j, k)
              dp1_ad(i, j, k) = 0.0
            END DO
          END DO
        END IF
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    CALL POPCONTROL(1,branch)
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      CALL NESTED_GRID_BC_APPLY_INTT_ADM(pt, pt_ad, 0, 0, npx, npy, npz&
&                                  , bd, 1., 1., neststruct%pt_bc, &
&                                  neststruct%nestbctype)
    ELSE IF (branch .NE. 1) THEN
      GOTO 120
    END IF
    CALL SETUP_NESTED_GRID_BCS_ADM(npx, npy, npz, zvir, ncnst, u, u_ad, &
&                            v, v_ad, w, pt, delp, delz, q, uc, uc_ad, &
&                            vc, vc_ad, pkz, neststruct%nested, &
&                            flagstruct%inline_q, flagstruct%make_nh, ng&
&                            , gridstruct, flagstruct, neststruct, &
&                            neststruct%nest_timestep, neststruct%&
&                            tracer_nest_timestep, domain, bd, nwat)
 120 CONTINUE
  END SUBROUTINE FV_DYNAMICS_BWD
!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
  SUBROUTINE FV_DYNAMICS(npx, npy, npz, nq_tot, ng, bdt, consv_te, fill&
&   , reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split, &
&   q_split, u, v, w, delz, hydrostatic, pt, delp, q, ps, pe, pk, peln, &
&   pkz, phis, q_con, omga, ua, va, uc, vc, ak, bk, mfx, mfy, cx, cy, &
&   ze0, hybrid_z, gridstruct, flagstruct, flagstructp, neststruct, &
&   idiag, bd, parent_grid, domain, time_total)
    IMPLICIT NONE
! Large time-step
    REAL, INTENT(IN) :: bdt
    REAL, INTENT(IN) :: consv_te
    REAL, INTENT(IN) :: kappa, cp_air
    REAL, INTENT(IN) :: zvir, ptop
    REAL, INTENT(IN), OPTIONAL :: time_total
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! transported tracers
    INTEGER, INTENT(IN) :: nq_tot
    INTEGER, INTENT(IN) :: ng
    INTEGER, INTENT(IN) :: ks
    INTEGER, INTENT(IN) :: ncnst
! small-step horizontal dynamics
    INTEGER, INTENT(IN) :: n_split
! tracer
    INTEGER, INTENT(IN) :: q_split
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
    LOGICAL, INTENT(IN) :: hydrostatic
! Using hybrid_z for remapping
    LOGICAL, INTENT(IN) :: hybrid_z
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
!  W (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! specific humidity and constituents
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! height at edges (m); non-hydrostatic
    REAL, INTENT(INOUT) :: ze0(bd%is:bd%is, bd%js:bd%js, 1)
! ze0 no longer used
!-----------------------------------------------------------------------
! Auxilliary pressure arrays:
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
! Surface pressure (pascal)
    REAL, INTENT(INOUT) :: ps(bd%isd:bd%ied, bd%jsd:bd%jed)
! edge pressure (pascal)
    REAL, INTENT(INOUT) :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
! pe**kappa
    REAL, INTENT(INOUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
! ln(pe)
    REAL, INTENT(INOUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
! finite-volume mean pk
    REAL, INTENT(INOUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: q_con(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! (uc,vc) mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua, va
    REAL, DIMENSION(npz+1), INTENT(IN) :: ak, bk
! Accumulated Mass flux arrays: the "Flux Capacitor"
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(INOUT) :: flagstruct
    TYPE(FV_FLAGS_PERT_TYPE), INTENT(INOUT) :: flagstructp
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(FV_DIAG_TYPE), INTENT(IN) :: idiag
! Local Arrays
    REAL :: ws(bd%is:bd%ie, bd%js:bd%je)
    REAL :: te_2d(bd%is:bd%ie, bd%js:bd%je)
    REAL :: teq(bd%is:bd%ie, bd%js:bd%je)
    REAL :: ps2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: m_fac(bd%is:bd%ie, bd%js:bd%je)
    REAL :: pfull(npz)
    REAL, DIMENSION(bd%is:bd%ie) :: cvm
    REAL :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz), dtdt_m(bd%is:bd%ie, &
&   bd%js:bd%je, npz), cappa(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL(kind=8) :: psx(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL(kind=8) :: dpx(bd%is:bd%ie, bd%js:bd%je)
    REAL :: akap, rdg, ph1, ph2, mdt, gam, amdt, u0
    INTEGER :: kord_tracer(ncnst), kord_mt, kord_wz, kord_tm
    INTEGER :: kord_tracer_pert(ncnst), kord_mt_pert, kord_wz_pert, &
&   kord_tm_pert
    INTEGER :: i, j, k, n, iq, n_map, nq, nwat, k_split
! GFDL physics
    INTEGER :: sphum
    INTEGER, SAVE :: liq_wat=-999
    INTEGER, SAVE :: ice_wat=-999
    INTEGER, SAVE :: rainwat=-999
    INTEGER, SAVE :: snowwat=-999
    INTEGER, SAVE :: graupel=-999
    INTEGER, SAVE :: cld_amt=-999
    INTEGER, SAVE :: theta_d=-999
    LOGICAL :: used, last_step, do_omega
    INTEGER, PARAMETER :: max_packs=12
    TYPE(GROUP_HALO_UPDATE_TYPE), SAVE :: i_pack(max_packs)
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL :: dt2
    REAL(kind=8) :: t1, t2
    INTEGER :: status
    REAL :: rf(npz)
    REAL :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pkc(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: ptc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: cry(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: divgd(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL :: delpc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: ut(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: zh(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pk3(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    INTRINSIC ANY
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC ABS
    INTRINSIC REAL
    INTRINSIC COS
    REAL :: abs0
    REAL :: abs1
    INTEGER :: arg1
    LOGICAL :: arg10
    REAL*8 :: arg11
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: arg12
    REAL :: result1
    gz = 0.0
    pkc = 0.0
    ptc = 0.0
    crx = 0.0
    xfx = 0.0
    cry = 0.0
    yfx = 0.0
    divgd = 0.0
    delpc = 0.0
    ut = 0.0
    vt = 0.0
    zh = 0.0
    pk3 = 0.0
    du = 0.0
    dv = 0.0
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    dyn_timer = 0
    comm_timer = 0
!     cv_air =  cp_air - rdgas
    agrav = 1./grav
    dt2 = 0.5*bdt
    k_split = flagstruct%k_split
    nwat = flagstruct%nwat
    nq = nq_tot - flagstruct%dnats
    rdg = -(rdgas*agrav)
!allocate ( dp1(isd:ied, jsd:jed, 1:npz) )
! Begin Dynamics timer for GEOS history processing
!t1 = MPI_Wtime(status)
    t1 = 0.0
    t2 = 0.0
!allocate ( cappa(isd:isd,jsd:jsd,1) )
    cappa = 0.
!We call this BEFORE converting pt to virtual potential temperature,
!since we interpolate on (regular) temperature rather than theta.
    IF (gridstruct%nested .OR. ANY(neststruct%child_grids)) THEN
      CALL TIMING_ON('NEST_BCs')
      CALL SETUP_NESTED_GRID_BCS(npx, npy, npz, zvir, ncnst, u, v, w, pt&
&                          , delp, delz, q, uc, vc, pkz, neststruct%&
&                          nested, flagstruct%inline_q, flagstruct%&
&                          make_nh, ng, gridstruct, flagstruct, &
&                          neststruct, neststruct%nest_timestep, &
&                          neststruct%tracer_nest_timestep, domain, bd, &
&                          nwat)
      IF (gridstruct%nested) CALL NESTED_GRID_BC_APPLY_INTT(pt, 0, 0, &
&                                                     npx, npy, npz, bd&
&                                                     , 1., 1., &
&                                                     neststruct%pt_bc, &
&                                                     neststruct%&
&                                                     nestbctype)
!Correct halo values have now been set up for BCs; we can go ahead and apply them too...
      CALL TIMING_OFF('NEST_BCs')
    END IF
    IF (flagstruct%no_dycore) THEN
      IF (nwat .EQ. 2 .AND. (.NOT.hydrostatic)) sphum = GET_TRACER_INDEX&
&         (model_atmos, 'sphum')
    END IF
!goto 911
    IF (fpp%fpp_mapl_mode) THEN
      SELECT CASE  (nwat) 
      CASE (0) 
        sphum = 1
! to cause trouble if (mis)used
        cld_amt = -1
      CASE (1) 
        sphum = 1
! to cause trouble if (mis)used
        liq_wat = -1
! to cause trouble if (mis)used
        ice_wat = -1
! to cause trouble if (mis)used
        rainwat = -1
! to cause trouble if (mis)used
        snowwat = -1
! to cause trouble if (mis)used
        graupel = -1
! to cause trouble if (mis)used
        cld_amt = -1
! to cause trouble if (mis)used
        theta_d = -1
      CASE (3) 
        sphum = 1
        liq_wat = 2
        ice_wat = 3
! to cause trouble if (mis)used
        rainwat = -1
! to cause trouble if (mis)used
        snowwat = -1
! to cause trouble if (mis)used
        graupel = -1
! to cause trouble if (mis)used
        cld_amt = -1
! to cause trouble if (mis)used
        theta_d = -1
      END SELECT
    ELSE
      IF (nwat .EQ. 0) THEN
        sphum = 1
! to cause trouble if (mis)used
        cld_amt = -1
      ELSE
        sphum = GET_TRACER_INDEX(model_atmos, 'sphum')
        liq_wat = GET_TRACER_INDEX(model_atmos, 'liq_wat')
        ice_wat = GET_TRACER_INDEX(model_atmos, 'ice_wat')
        rainwat = GET_TRACER_INDEX(model_atmos, 'rainwat')
        snowwat = GET_TRACER_INDEX(model_atmos, 'snowwat')
        graupel = GET_TRACER_INDEX(model_atmos, 'graupel')
        cld_amt = GET_TRACER_INDEX(model_atmos, 'cld_amt')
      END IF
      theta_d = GET_TRACER_INDEX(model_atmos, 'theta_d')
    END IF
    akap = kappa
!$OMP parallel do default(none) shared(npz,ak,bk,flagstruct,pfull) &
!$OMP                          private(ph1, ph2)
    DO k=1,npz
      ph1 = ak(k) + bk(k)*flagstruct%p_ref
      ph2 = ak(k+1) + bk(k+1)*flagstruct%p_ref
      pfull(k) = (ph2-ph1)/LOG(ph2/ph1)
    END DO
    IF (hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,zvir,nwat,q,q_con,sphum,liq_wat, &
!$OMP      rainwat,ice_wat,snowwat,graupel) private(cvm)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dp1(i, j, k) = zvir*q(i, j, k, sphum)
          END DO
        END DO
      END DO
    ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,zvir,q,q_con,sphum,liq_wat, &
!$OMP                                  rainwat,ice_wat,snowwat,graupel,pkz,flagstruct, &
!$OMP                                  cappa,kappa,rdg,delp,pt,delz,nwat)              &
!$OMP                          private(cvm)
      DO k=1,npz
        IF (flagstruct%moist_phys) THEN
          DO j=js,je
            DO i=is,ie
              dp1(i, j, k) = zvir*q(i, j, k, sphum)
              pkz(i, j, k) = EXP(kappa*LOG(rdg*delp(i, j, k)*pt(i, j, k)&
&               *(1.+dp1(i, j, k))/delz(i, j, k)))
            END DO
          END DO
        ELSE
! Using dry pressure for the definition of the virtual potential temperature
!              pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!                                      (1.-q(i,j,k,sphum))/delz(i,j,k)) )
          DO j=js,je
            DO i=is,ie
              dp1(i, j, k) = 0.
              pkz(i, j, k) = EXP(kappa*LOG(rdg*delp(i, j, k)*pt(i, j, k)&
&               /delz(i, j, k)))
            END DO
          END DO
        END IF
      END DO
    END IF
    IF (flagstruct%fv_debug) THEN
      CALL PRT_MXM('PS', ps, is, ie, js, je, ng, 1, 0.01, gridstruct%&
&            area_64, domain)
      CALL PRT_MXM('T_dyn_b', pt, is, ie, js, je, ng, npz, 1., &
&            gridstruct%area_64, domain)
      IF (.NOT.hydrostatic) CALL PRT_MXM('delz', delz, is, ie, js, je, &
&                                  ng, npz, 1., gridstruct%area_64, &
&                                  domain)
      CALL PRT_MXM('delp_b ', delp, is, ie, js, je, ng, npz, 0.01, &
&            gridstruct%area_64, domain)
      arg1 = npz + 1
      CALL PRT_MXM('pk_b', pk, is, ie, js, je, 0, arg1, 1., gridstruct%&
&            area_64, domain)
      CALL PRT_MXM('pkz_b', pkz, is, ie, js, je, 0, npz, 1., gridstruct%&
&            area_64, domain)
    END IF
!---------------------
! Compute Total Energy
!---------------------
    IF (consv_te .GT. 0. .AND. (.NOT.do_adiabatic_init)) THEN
      CALL COMPUTE_TOTAL_ENERGY(is, ie, js, je, isd, ied, jsd, jed, npz&
&                         , u, v, w, delz, pt, delp, q, dp1, pe, peln, &
&                         phis, gridstruct%rsin2, gridstruct%cosa_s, &
&                         zvir, cp_air, rdgas, hlv, te_2d, ua, va, teq, &
&                         flagstruct%moist_phys, nwat, sphum, liq_wat, &
&                         rainwat, ice_wat, snowwat, graupel, &
&                         hydrostatic, idiag%id_te)
      IF (idiag%id_te .GT. 0) used = SEND_DATA(idiag%id_te, teq, fv_time&
&         )
!              te_den=1.E-9*g_sum(teq, is, ie, js, je, ng, area, 0)/(grav*4.*pi*radius**2)
!              if(is_master())  write(*,*) 'Total Energy Density (Giga J/m**2)=',te_den
    END IF
    IF ((flagstruct%consv_am .OR. idiag%id_amdt .GT. 0) .AND. (.NOT.&
&       do_adiabatic_init)) CALL COMPUTE_AAM(npz, is, ie, js, je, isd, &
&                                      ied, jsd, jed, gridstruct, bd, &
&                                      ptop, ua, va, u, v, delp, teq, &
&                                      ps2, m_fac)
    IF (flagstruct%tau .GT. 0.) THEN
      IF (gridstruct%grid_type .LT. 4) THEN
        IF (bdt .GE. 0.) THEN
          abs0 = bdt
        ELSE
          abs0 = -bdt
        END IF
        arg10 = .NOT.neststruct%nested
        CALL RAYLEIGH_SUPER(abs0, npx, npy, npz, ks, pfull, phis, &
&                     flagstruct%tau, u, v, w, pt, ua, va, delz, &
&                     gridstruct%agrid, cp_air, rdgas, ptop, hydrostatic&
&                     , arg10, flagstruct%rf_cutoff, rf, gridstruct, &
&                     domain, bd)
      ELSE
        IF (bdt .GE. 0.) THEN
          abs1 = bdt
        ELSE
          abs1 = -bdt
        END IF
        CALL RAYLEIGH_FRICTION(abs1, npx, npy, npz, ks, pfull, &
&                        flagstruct%tau, u, v, w, pt, ua, va, delz, &
&                        cp_air, rdgas, ptop, hydrostatic, .true., &
&                        flagstruct%rf_cutoff, rf, gridstruct, domain, &
&                        bd)
      END IF
    END IF
! Convert pt to virtual potential temperature on the first timestep
    IF (flagstruct%adiabatic) THEN
!$OMP parallel do default(none) shared(theta_d,is,ie,js,je,npz,pt,pkz,q)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            pt(i, j, k) = pt(i, j, k)/pkz(i, j, k)
          END DO
        END DO
        IF (theta_d .GT. 0) THEN
          DO j=js,je
            DO i=is,ie
              q(i, j, k, theta_d) = pt(i, j, k)
            END DO
          END DO
        END IF
      END DO
    ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,dp1,pkz,q_con)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            pt(i, j, k) = pt(i, j, k)*(1.+dp1(i, j, k))/pkz(i, j, k)
          END DO
        END DO
      END DO
    END IF
    last_step = .false.
    mdt = bdt/REAL(k_split)
    IF (idiag%id_mdt .GT. 0 .AND. (.NOT.do_adiabatic_init)) THEN
!allocate ( dtdt_m(is:ie,js:je,npz) )
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dtdt_m)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dtdt_m(i, j, k) = 0.
          END DO
        END DO
      END DO
    END IF
!DryMassRoundoffControl
!allocate(psx(isd:ied,jsd:jed),dpx(is:ie,js:je))
    IF (fpp%fpp_overload_r4) THEN
      DO j=js,je
        DO i=is,ie
          psx(i, j) = pe(i, npz+1, j)
          dpx(i, j) = 0.0
        END DO
      END DO
    END IF
    CALL TIMING_ON('FV_DYN_LOOP')
! first level of time-split
    DO n_map=1,k_split
      CALL TIMING_ON('COMM_TOTAL')
      CALL START_GROUP_HALO_UPDATE(i_pack(1), delp, domain, complete=&
&                            .true.)
      CALL START_GROUP_HALO_UPDATE(i_pack(1), pt, domain, complete=&
&                            .true.)
      CALL START_GROUP_HALO_UPDATE(i_pack(8), u, v, domain, gridtype=&
&                            dgrid_ne)
      CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,dp1,delp)
      DO k=1,npz
        DO j=jsd,jed
          DO i=isd,ied
            dp1(i, j, k) = delp(i, j, k)
          END DO
        END DO
      END DO
      IF (n_map .EQ. k_split) last_step = .true.
      CALL TIMING_ON('DYN_CORE')
      arg10 = n_map .EQ. 1
      CALL DYN_CORE(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir, &
&             cp_air, akap, cappa, grav, hydrostatic, u, v, w, delz, pt&
&             , q, delp, pe, pk, phis, ws, omga, ptop, pfull, ua, va, uc&
&             , vc, mfx, mfy, cx, cy, pkz, peln, q_con, ak, bk, dpx, ks&
&             , gridstruct, flagstruct, flagstructp, neststruct, idiag, &
&             bd, domain, arg10, i_pack, last_step, gz, pkc, ptc, crx, &
&             xfx, cry, yfx, divgd, delpc, ut, vt, zh, pk3, du, dv, &
&             time_total)
      CALL TIMING_OFF('DYN_CORE')
!DryMassRoundoffControl
      IF (last_step) THEN
        IF (fpp%fpp_overload_r4) THEN
          DO j=js,je
            DO i=is,ie
              psx(i, j) = psx(i, j) + dpx(i, j)
            END DO
          END DO
          CALL TIMING_ON('COMM_TOTAL')
          CALL MPP_UPDATE_DOMAINS(psx, domain)
          CALL TIMING_OFF('COMM_TOTAL')
          DO j=js-1,je+1
            DO i=is-1,ie+1
              pe(i, npz+1, j) = psx(i, j)
            END DO
          END DO
        END IF
      END IF
!deallocate(psx,dpx)
      IF (.NOT.flagstruct%inline_q .AND. nq .NE. 0) THEN
!--------------------------------------------------------
! Perform large-time-step scalar transport using the accumulated CFL and
! mass fluxes
        CALL TIMING_ON('tracer_2d')
!!! CLEANUP: merge these two calls?
        IF (gridstruct%nested) THEN
          CALL TRACER_2D_NESTED(q, dp1, mfx, mfy, cx, cy, gridstruct, bd&
&                         , domain, npx, npy, npz, nq, flagstruct%&
&                         hord_tr, q_split, mdt, idiag%id_divg, i_pack(&
&                         10), flagstruct%nord_tr, flagstruct%trdm2, &
&                         k_split, neststruct, parent_grid, flagstructp%&
&                         hord_tr_pert, flagstructp%nord_tr_pert, &
&                         flagstructp%trdm2_pert, flagstructp%&
&                         split_damp_tr)
        ELSE IF (flagstruct%z_tracer) THEN
          CALL TRACER_2D_1L(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, &
&                     domain, npx, npy, npz, nq, flagstruct%hord_tr, &
&                     q_split, mdt, idiag%id_divg, i_pack(10), &
&                     flagstruct%nord_tr, flagstruct%trdm2, flagstructp%&
&                     hord_tr_pert, flagstructp%nord_tr_pert, &
&                     flagstructp%trdm2_pert, flagstructp%split_damp_tr)
        ELSE
          CALL TRACER_2D(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, &
&                  domain, npx, npy, npz, nq, flagstruct%hord_tr, &
&                  q_split, mdt, idiag%id_divg, i_pack(10), flagstruct%&
&                  nord_tr, flagstruct%trdm2, flagstructp%hord_tr_pert, &
&                  flagstructp%nord_tr_pert, flagstructp%trdm2_pert, &
&                  flagstructp%split_damp_tr)
        END IF
        CALL TIMING_OFF('tracer_2d')
        IF (flagstruct%hord_tr .LT. 8 .AND. flagstruct%moist_phys) THEN
          CALL TIMING_ON('Fill2D')
          IF (liq_wat .GT. 0) CALL FILL2D(is, ie, js, je, ng, npz, q(isd&
&                                   :ied, jsd:jed, 1, liq_wat), delp, &
&                                   gridstruct%area, domain, neststruct%&
&                                   nested, npx, npy)
          IF (rainwat .GT. 0) CALL FILL2D(is, ie, js, je, ng, npz, q(isd&
&                                   :ied, jsd:jed, 1, rainwat), delp, &
&                                   gridstruct%area, domain, neststruct%&
&                                   nested, npx, npy)
          IF (ice_wat .GT. 0) CALL FILL2D(is, ie, js, je, ng, npz, q(isd&
&                                   :ied, jsd:jed, 1, ice_wat), delp, &
&                                   gridstruct%area, domain, neststruct%&
&                                   nested, npx, npy)
          IF (snowwat .GT. 0) CALL FILL2D(is, ie, js, je, ng, npz, q(isd&
&                                   :ied, jsd:jed, 1, snowwat), delp, &
&                                   gridstruct%area, domain, neststruct%&
&                                   nested, npx, npy)
          IF (graupel .GT. 0) CALL FILL2D(is, ie, js, je, ng, npz, q(isd&
&                                   :ied, jsd:jed, 1, graupel), delp, &
&                                   gridstruct%area, domain, neststruct%&
&                                   nested, npx, npy)
          CALL TIMING_OFF('Fill2D')
        END IF
        IF (last_step .AND. idiag%id_divg .GT. 0) THEN
          used = SEND_DATA(idiag%id_divg, dp1, fv_time)
          IF (flagstruct%fv_debug) CALL PRT_MXM('divg', dp1, is, ie, js&
&                                         , je, 0, npz, 1., gridstruct%&
&                                         area_64, domain)
        END IF
      END IF
      IF (npz .GT. 4) THEN
!------------------------------------------------------------------------
! Perform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.
!------------------------------------------------------------------------
        DO iq=1,nq
          kord_tracer(iq) = flagstruct%kord_tr
! monotonic
          IF (iq .EQ. cld_amt) kord_tracer(iq) = 9
          kord_tracer_pert(iq) = flagstructp%kord_tr_pert
! linear
          IF (iq .EQ. cld_amt) kord_tracer_pert(iq) = 17
        END DO
        do_omega = hydrostatic .AND. last_step
        CALL TIMING_ON('Remapping')
        kord_mt = flagstruct%kord_mt
        kord_wz = flagstruct%kord_wz
        kord_tm = flagstruct%kord_tm
        kord_mt_pert = flagstructp%kord_mt_pert
        kord_wz_pert = flagstructp%kord_wz_pert
        kord_tm_pert = flagstructp%kord_tm_pert
        IF (n_map .EQ. k_split) THEN
          kord_mt = kord_mt_pert
          kord_wz = kord_wz_pert
          kord_tm = kord_tm_pert
          kord_tracer = kord_tracer_pert
        END IF
        arg10 = idiag%id_mdt .GT. 0
        CALL LAGRANGIAN_TO_EULERIAN(last_step, consv_te, ps, pe, delp, &
&                             pkz, pk, mdt, bdt, npz, is, ie, js, je, &
&                             isd, ied, jsd, jed, nq, nwat, sphum, q_con&
&                             , u, v, w, delz, pt, q, phis, zvir, cp_air&
&                             , akap, cappa, kord_mt, kord_wz, &
&                             kord_tracer, kord_tm, peln, te_2d, ng, ua&
&                             , va, omga, dp1, ws, fill, reproduce_sum, &
&                             arg10, dtdt_m, ptop, ak, bk, pfull, &
&                             flagstruct, gridstruct, domain, flagstruct&
&                             %do_sat_adj, hydrostatic, hybrid_z, &
&                             do_omega, flagstruct%adiabatic, &
&                             do_adiabatic_init, mfx, mfy, flagstruct%&
&                             remap_option, kord_mt_pert, kord_wz_pert, &
&                             kord_tracer_pert, kord_tm_pert)
        CALL TIMING_OFF('Remapping')
        IF (last_step) THEN
          IF (.NOT.hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga,delp,delz,w)
            DO k=1,npz
              DO j=js,je
                DO i=is,ie
                  omga(i, j, k) = delp(i, j, k)/delz(i, j, k)*w(i, j, k)
                END DO
              END DO
            END DO
          END IF
!--------------------------
! Filter omega for physics:
!--------------------------
          IF (flagstruct%nf_omega .GT. 0) THEN
            arg11 = 0.18*gridstruct%da_min
            CALL DEL2_CUBED(omga, arg11, gridstruct, domain, npx, npy, &
&                     npz, flagstruct%nf_omega, bd)
          END IF
        END IF
      END IF
    END DO
! n_map loop
    CALL TIMING_OFF('FV_DYN_LOOP')
    IF (idiag%id_mdt .GT. 0 .AND. (.NOT.do_adiabatic_init)) THEN
! Output temperature tendency due to inline moist physics:
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dtdt_m,bdt)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dtdt_m(i, j, k) = dtdt_m(i, j, k)/bdt*86400.
          END DO
        END DO
      END DO
!      call prt_mxm('Fast DTDT (deg/Day)', dtdt_m, is, ie, js, je, 0, npz, 1., gridstruct%area_64, domain)
      used = SEND_DATA(idiag%id_mdt, dtdt_m, fv_time)
!deallocate ( dtdt_m )
    END IF
    IF (nwat .EQ. 6) THEN
      IF (cld_amt .GT. 0) THEN
        CALL NEG_ADJ3(is, ie, js, je, ng, npz, flagstruct%hydrostatic, &
&               peln, delz, pt, delp, q(isd:ied, jsd:jed, 1, sphum), q(&
&               isd:ied, jsd:jed, 1, liq_wat), q(isd:ied, jsd:jed, 1, &
&               rainwat), q(isd:ied, jsd:jed, 1, ice_wat), q(isd:ied, &
&               jsd:jed, 1, snowwat), q(isd:ied, jsd:jed, 1, graupel), q&
&               (isd:ied, jsd:jed, 1, cld_amt), flagstruct%&
&               check_negative)
      ELSE
        CALL NEG_ADJ3(is, ie, js, je, ng, npz, flagstruct%hydrostatic, &
&               peln, delz, pt, delp, q(isd:ied, jsd:jed, 1, sphum), q(&
&               isd:ied, jsd:jed, 1, liq_wat), q(isd:ied, jsd:jed, 1, &
&               rainwat), q(isd:ied, jsd:jed, 1, ice_wat), q(isd:ied, &
&               jsd:jed, 1, snowwat), q(isd:ied, jsd:jed, 1, graupel), &
&               check_negative=flagstruct%check_negative)
      END IF
      IF (flagstruct%fv_debug) THEN
        CALL PRT_MXM('T_dyn_a3', pt, is, ie, js, je, ng, npz, 1., &
&              gridstruct%area_64, domain)
        CALL PRT_MXM('SPHUM_dyn', q(isd:ied, jsd:jed, 1, sphum), is, ie&
&              , js, je, ng, npz, 1., gridstruct%area_64, domain)
        CALL PRT_MXM('liq_wat_dyn', q(isd:ied, jsd:jed, 1, liq_wat), is&
&              , ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
        CALL PRT_MXM('rainwat_dyn', q(isd:ied, jsd:jed, 1, rainwat), is&
&              , ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
        CALL PRT_MXM('ice_wat_dyn', q(isd:ied, jsd:jed, 1, ice_wat), is&
&              , ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
        CALL PRT_MXM('snowwat_dyn', q(isd:ied, jsd:jed, 1, snowwat), is&
&              , ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
        CALL PRT_MXM('graupel_dyn', q(isd:ied, jsd:jed, 1, graupel), is&
&              , ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
      END IF
    END IF
    IF (((flagstruct%consv_am .OR. idiag%id_amdt .GT. 0) .OR. idiag%&
&       id_aam .GT. 0) .AND. (.NOT.do_adiabatic_init)) THEN
      CALL COMPUTE_AAM(npz, is, ie, js, je, isd, ied, jsd, jed, &
&                gridstruct, bd, ptop, ua, va, u, v, delp, te_2d, ps, &
&                m_fac)
      IF (idiag%id_aam .GT. 0) THEN
        used = SEND_DATA(idiag%id_aam, te_2d, fv_time)
        IF (prt_minmax) gam = G_SUM(domain, te_2d, is, ie, js, je, ng, &
&           gridstruct%area_64, 0)
!if( is_master() ) write(6,*) 'Total AAM =', gam
      END IF
    END IF
    IF ((flagstruct%consv_am .OR. idiag%id_amdt .GT. 0) .AND. (.NOT.&
&       do_adiabatic_init)) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,te_2d,teq,dt2,ps2,ps,idiag)
      DO j=js,je
        DO i=is,ie
! Note: the mountain torque computation contains also numerical error
! The numerical error is mostly from the zonal gradient of the terrain (zxg)
          te_2d(i, j) = te_2d(i, j) - teq(i, j) + dt2*(ps2(i, j)+ps(i, j&
&           ))*idiag%zxg(i, j)
        END DO
      END DO
      IF (idiag%id_amdt .GT. 0) THEN
        arg12(:, :) = te_2d/bdt
        used = SEND_DATA(idiag%id_amdt, arg12(:, :), fv_time)
      END IF
      IF (flagstruct%consv_am .OR. prt_minmax) THEN
        amdt = G_SUM(domain, te_2d, is, ie, js, je, ng, gridstruct%&
&         area_64, 0, .true.)
        result1 = G_SUM(domain, m_fac, is, ie, js, je, ng, gridstruct%&
&         area_64, 0, .true.)
        u0 = -(radius*amdt/result1)
        IF (IS_MASTER() .AND. prt_minmax) WRITE(6, *) &
&                                       'Dynamic AM tendency (Hadleys)='&
&                                         , amdt/(bdt*1.e18), &
&                                         'del-u (per day)=', u0*86400./&
&                                         bdt
      END IF
!  consv_am
      IF (flagstruct%consv_am) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,m_fac,u0,gridstruct)
        DO j=js,je
          DO i=is,ie
            m_fac(i, j) = u0*COS(gridstruct%agrid(i, j, 2))
          END DO
        END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,hydrostatic,pt,m_fac,ua,cp_air, &
!$OMP                                  u,u0,gridstruct,v )
        DO k=1,npz
          DO j=js,je+1
            DO i=is,ie
              u(i, j, k) = u(i, j, k) + u0*gridstruct%l2c_u(i, j)
            END DO
          END DO
          DO j=js,je
            DO i=is,ie+1
              v(i, j, k) = v(i, j, k) + u0*gridstruct%l2c_v(i, j)
            END DO
          END DO
        END DO
      END IF
    END IF
!911  call cubed_to_latlon(u, v, ua, va, gridstruct, &
!          npx, npy, npz, 1, gridstruct%grid_type, domain, gridstruct%nested, flagstruct%c2l_ord, bd)
!deallocate(dp1)
!deallocate(cappa)
    IF (flagstruct%fv_debug) THEN
      CALL PRT_MXM('UA', ua, is, ie, js, je, ng, npz, 1., gridstruct%&
&            area_64, domain)
      CALL PRT_MXM('VA', va, is, ie, js, je, ng, npz, 1., gridstruct%&
&            area_64, domain)
      CALL PRT_MXM('TA', pt, is, ie, js, je, ng, npz, 1., gridstruct%&
&            area_64, domain)
      IF (.NOT.hydrostatic) CALL PRT_MXM('W ', w, is, ie, js, je, ng, &
&                                  npz, 1., gridstruct%area_64, domain)
    END IF
    IF (flagstruct%range_warn) THEN
      CALL RANGE_CHECK('UA_dyn', ua, is, ie, js, je, ng, npz, gridstruct&
&                %agrid, -280., 280., bad_range)
      CALL RANGE_CHECK('VA_dyn', ua, is, ie, js, je, ng, npz, gridstruct&
&                %agrid, -280., 280., bad_range)
      CALL RANGE_CHECK('TA_dyn', pt, is, ie, js, je, ng, npz, gridstruct&
&                %agrid, 150., 335., bad_range)
      IF (.NOT.hydrostatic) CALL RANGE_CHECK('W_dyn', w, is, ie, js, je&
&                                      , ng, npz, gridstruct%agrid, -50.&
&                                      , 100., bad_range)
    END IF
    IF (fpp%fpp_mapl_mode) dyn_timer = dyn_timer + (t2-t1)
!t2 = MPI_Wtime(status)
  END SUBROUTINE FV_DYNAMICS
!  Differentiation of rayleigh_super in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
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
!   gradient     of useful results: u v w ua va pt
!   with respect to varying inputs: u v w ua va pt
  SUBROUTINE RAYLEIGH_SUPER_FWD(dt, npx, npy, npz, ks, pm, phis, tau, u&
&   , v, w, pt, ua, va, delz, agrid, cp, rg, ptop, hydrostatic, conserve&
&   , rf_cutoff, rf, gridstruct, domain, bd)
    IMPLICIT NONE
!deallocate ( u2f )
    REAL, INTENT(IN) :: dt
! time scale (days)
    REAL, INTENT(IN) :: tau
    REAL, INTENT(IN) :: cp, rg, ptop, rf_cutoff
    INTEGER, INTENT(IN) :: npx, npy, npz, ks
    REAL, DIMENSION(npz), INTENT(IN) :: pm
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: conserve
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
! cell center vertical wind (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temp
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: rf(npz)
    REAL, INTENT(IN) :: agrid(bd%isd:bd%ied, bd%jsd:bd%jed, 2)
! Surface geopotential (g*Z_surf)
    REAL, INTENT(IN) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
!
    REAL :: u2f(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! scaling velocity
    REAL, PARAMETER :: u0=60.
    REAL, PARAMETER :: sday=86400.
    REAL :: rcv, tau0
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC LOG
    INTRINSIC SIN
    LOGICAL :: res
    INTEGER :: ad_count

    u2f = 0.0
    rcv = 0.0
    tau0 = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    ad_count = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    rcv = 1./(cp-rg)
    IF (.NOT.rf_initialized) THEN
      tau0 = tau*sday
!allocate( rf(npz) )
      rf(:) = 0.
      k = ks + 2
!if( is_master() ) write(6,*) k, 0.01*pm(k)
      res = IS_MASTER()
      IF (res) WRITE(6, *) 'Rayleigh friction E-folding time (days):'
      ad_count = 1
      DO k=1,npz
        IF (pm(k) .LT. rf_cutoff) THEN
          rf(k) = dt/tau0*SIN(0.5*pi*LOG(rf_cutoff/pm(k))/LOG(rf_cutoff/&
&           ptop))**2
!if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
          kmax = k
          ad_count = ad_count + 1
        ELSE
          GOTO 100
        END IF
      END DO
      CALL PUSHCONTROL(1,0)
      CALL PUSHINTEGER(ad_count)
      CALL PUSHCONTROL(1,0)
      GOTO 110
 100  CALL PUSHCONTROL(1,1)
      CALL PUSHINTEGER(ad_count)
      CALL PUSHCONTROL(1,0)
 110  rf_initialized = .true.
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    CALL C2L_ORD2_FWD(u, v, ua, va, gridstruct, npz, gridstruct%&
&               grid_type, bd, gridstruct%nested)
!allocate( u2f(isd:ied,jsd:jed,kmax) )
    u2f = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,pm,rf_cutoff,hydrostatic,ua,va,agrid, &
!$OMP                                  u2f,rf,w)
    DO k=1,kmax
      IF (pm(k) .LT. rf_cutoff) THEN
        CALL PUSHCONTROL(1,1)
        u2f(:, :, k) = 1./(1.+rf(k))
      ELSE
        CALL PUSHCONTROL(1,0)
        u2f(:, :, k) = 1.
      END IF
    END DO
    CALL MPP_UPDATE_DOMAINS(u2f, domain)
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,pm,rf_cutoff,w,rf,u,v, &
!$OMP                                  conserve,hydrostatic,pt,ua,va,u2f,cp,rg,ptop,rcv)
    DO k=1,kmax
      IF (pm(k) .LT. rf_cutoff) THEN
! Add heat so as to conserve TE
        IF (conserve) THEN
          IF (hydrostatic) THEN
            DO j=js,je
              DO i=is,ie
                CALL PUSHREALARRAY(pt(i, j, k))
                pt(i, j, k) = pt(i, j, k) + 0.5*(ua(i, j, k)**2+va(i, j&
&                 , k)**2)*(1.-u2f(i, j, k)**2)/(cp-rg*ptop/pm(k))
              END DO
            END DO
            CALL PUSHCONTROL(2,2)
          ELSE
            DO j=js,je
              DO i=is,ie
                CALL PUSHREALARRAY(pt(i, j, k))
                pt(i, j, k) = pt(i, j, k) + 0.5*(ua(i, j, k)**2+va(i, j&
&                 , k)**2+w(i, j, k)**2)*(1.-u2f(i, j, k)**2)*rcv
              END DO
            END DO
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE
          CALL PUSHCONTROL(2,0)
        END IF
        DO j=js,je+1
          DO i=is,ie
            CALL PUSHREALARRAY(u(i, j, k))
            u(i, j, k) = 0.5*(u2f(i, j-1, k)+u2f(i, j, k))*u(i, j, k)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            CALL PUSHREALARRAY(v(i, j, k))
            v(i, j, k) = 0.5*(u2f(i-1, j, k)+u2f(i, j, k))*v(i, j, k)
          END DO
        END DO
        IF (.NOT.hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(w(i, j, k))
              w(i, j, k) = u2f(i, j, k)*w(i, j, k)
            END DO
          END DO
          CALL PUSHCONTROL(2,2)
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,0)
      END IF
    END DO
    CALL PUSHREALARRAY(u2f, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(ie)
    CALL PUSHREALARRAY(rcv)
    CALL PUSHINTEGER(js)
  END SUBROUTINE RAYLEIGH_SUPER_FWD
!  Differentiation of rayleigh_super in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_e
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
!   gradient     of useful results: u v w ua va pt
!   with respect to varying inputs: u v w ua va pt
  SUBROUTINE RAYLEIGH_SUPER_BWD(dt, npx, npy, npz, ks, pm, phis, tau, u&
&   , u_ad, v, v_ad, w, w_ad, pt, pt_ad, ua, ua_ad, va, va_ad, delz, &
&   agrid, cp, rg, ptop, hydrostatic, conserve, rf_cutoff, rf, &
&   gridstruct, domain, bd)
    IMPLICIT NONE
!deallocate ( u2f )
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: tau
    REAL, INTENT(IN) :: cp, rg, ptop, rf_cutoff
    INTEGER, INTENT(IN) :: npx, npy, npz, ks
    REAL, DIMENSION(npz), INTENT(IN) :: pm
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: conserve
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: w_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ua_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: va_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: rf(npz)
    REAL, INTENT(IN) :: agrid(bd%isd:bd%ied, bd%jsd:bd%jed, 2)
    REAL, INTENT(IN) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL :: u2f(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, PARAMETER :: u0=60.
    REAL, PARAMETER :: sday=86400.
    REAL :: rcv, tau0
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC LOG
    INTRINSIC SIN
    REAL :: temp_ad
    REAL :: temp_ad0
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch

    u2f = 0.0
    rcv = 0.0
    tau0 = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    ad_count = 0
    branch = 0

    CALL POPINTEGER(js)
    CALL POPREALARRAY(rcv)
    CALL POPINTEGER(ie)
    CALL POPINTEGER(is)
    CALL POPINTEGER(je)
    CALL POPREALARRAY(u2f, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    DO k=kmax,1,-1
      CALL POPCONTROL(2,branch)
      IF (branch .NE. 0) THEN
        IF (branch .NE. 1) THEN
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(w(i, j, k))
              w_ad(i, j, k) = u2f(i, j, k)*w_ad(i, j, k)
            END DO
          END DO
        END IF
        DO j=je,js,-1
          DO i=ie+1,is,-1
            CALL POPREALARRAY(v(i, j, k))
            v_ad(i, j, k) = (u2f(i-1, j, k)+u2f(i, j, k))*0.5*v_ad(i, j&
&             , k)
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(u(i, j, k))
            u_ad(i, j, k) = (u2f(i, j-1, k)+u2f(i, j, k))*0.5*u_ad(i, j&
&             , k)
          END DO
        END DO
        CALL POPCONTROL(2,branch)
        IF (branch .NE. 0) THEN
          IF (branch .EQ. 1) THEN
            DO j=je,js,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(pt(i, j, k))
                temp_ad0 = (1.-u2f(i, j, k)**2)*rcv*0.5*pt_ad(i, j, k)
                ua_ad(i, j, k) = ua_ad(i, j, k) + 2*ua(i, j, k)*temp_ad0
                va_ad(i, j, k) = va_ad(i, j, k) + 2*va(i, j, k)*temp_ad0
                w_ad(i, j, k) = w_ad(i, j, k) + 2*w(i, j, k)*temp_ad0
              END DO
            END DO
          ELSE
            DO j=je,js,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(pt(i, j, k))
                temp_ad = (1.-u2f(i, j, k)**2)*0.5*pt_ad(i, j, k)/(cp-rg&
&                 *ptop/pm(k))
                ua_ad(i, j, k) = ua_ad(i, j, k) + 2*ua(i, j, k)*temp_ad
                va_ad(i, j, k) = va_ad(i, j, k) + 2*va(i, j, k)*temp_ad
              END DO
            END DO
          END IF
        END IF
      END IF
    END DO
    DO k=kmax,1,-1
      CALL POPCONTROL(1,branch)
    END DO
    CALL C2L_ORD2_BWD(u, u_ad, v, v_ad, ua, ua_ad, va, va_ad, gridstruct&
&               , npz, gridstruct%grid_type, bd, gridstruct%nested)
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(ad_count)
      DO i0=1,ad_count
        IF (i0 .EQ. 1) CALL POPCONTROL(1,branch)
      END DO
    END IF
  END SUBROUTINE RAYLEIGH_SUPER_BWD
  SUBROUTINE RAYLEIGH_SUPER(dt, npx, npy, npz, ks, pm, phis, tau, u, v, &
&   w, pt, ua, va, delz, agrid, cp, rg, ptop, hydrostatic, conserve, &
&   rf_cutoff, rf, gridstruct, domain, bd)
    IMPLICIT NONE
!deallocate ( u2f )
    REAL, INTENT(IN) :: dt
! time scale (days)
    REAL, INTENT(IN) :: tau
    REAL, INTENT(IN) :: cp, rg, ptop, rf_cutoff
    INTEGER, INTENT(IN) :: npx, npy, npz, ks
    REAL, DIMENSION(npz), INTENT(IN) :: pm
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: conserve
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
! cell center vertical wind (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temp
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: rf(npz)
    REAL, INTENT(IN) :: agrid(bd%isd:bd%ied, bd%jsd:bd%jed, 2)
! Surface geopotential (g*Z_surf)
    REAL, INTENT(IN) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
!
    REAL :: u2f(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! scaling velocity
    REAL, PARAMETER :: u0=60.
    REAL, PARAMETER :: sday=86400.
    REAL :: rcv, tau0
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC LOG
    INTRINSIC SIN
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    rcv = 1./(cp-rg)
    IF (.NOT.rf_initialized) THEN
      tau0 = tau*sday
!allocate( rf(npz) )
      rf(:) = 0.
      k = ks + 2
!if( is_master() ) write(6,*) k, 0.01*pm(k)
      IF (IS_MASTER()) WRITE(6, *) &
&                      'Rayleigh friction E-folding time (days):'
      DO k=1,npz
        IF (pm(k) .LT. rf_cutoff) THEN
          rf(k) = dt/tau0*SIN(0.5*pi*LOG(rf_cutoff/pm(k))/LOG(rf_cutoff/&
&           ptop))**2
!if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
          kmax = k
        ELSE
          GOTO 100
        END IF
      END DO
 100  rf_initialized = .true.
    END IF
    CALL C2L_ORD2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, &
&           bd, gridstruct%nested)
!allocate( u2f(isd:ied,jsd:jed,kmax) )
    u2f = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,pm,rf_cutoff,hydrostatic,ua,va,agrid, &
!$OMP                                  u2f,rf,w)
    DO k=1,kmax
      IF (pm(k) .LT. rf_cutoff) THEN
        u2f(:, :, k) = 1./(1.+rf(k))
      ELSE
        u2f(:, :, k) = 1.
      END IF
    END DO
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_DOMAINS(u2f, domain)
    CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,pm,rf_cutoff,w,rf,u,v, &
!$OMP                                  conserve,hydrostatic,pt,ua,va,u2f,cp,rg,ptop,rcv)
    DO k=1,kmax
      IF (pm(k) .LT. rf_cutoff) THEN
! Add heat so as to conserve TE
        IF (conserve) THEN
          IF (hydrostatic) THEN
            DO j=js,je
              DO i=is,ie
                pt(i, j, k) = pt(i, j, k) + 0.5*(ua(i, j, k)**2+va(i, j&
&                 , k)**2)*(1.-u2f(i, j, k)**2)/(cp-rg*ptop/pm(k))
              END DO
            END DO
          ELSE
            DO j=js,je
              DO i=is,ie
                pt(i, j, k) = pt(i, j, k) + 0.5*(ua(i, j, k)**2+va(i, j&
&                 , k)**2+w(i, j, k)**2)*(1.-u2f(i, j, k)**2)*rcv
              END DO
            END DO
          END IF
        END IF
        DO j=js,je+1
          DO i=is,ie
            u(i, j, k) = 0.5*(u2f(i, j-1, k)+u2f(i, j, k))*u(i, j, k)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            v(i, j, k) = 0.5*(u2f(i-1, j, k)+u2f(i, j, k))*v(i, j, k)
          END DO
        END DO
        IF (.NOT.hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              w(i, j, k) = u2f(i, j, k)*w(i, j, k)
            END DO
          END DO
        END IF
      END IF
    END DO
  END SUBROUTINE RAYLEIGH_SUPER
!  Differentiation of rayleigh_friction in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b
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
!   gradient     of useful results: u v w ua delz va pt
!   with respect to varying inputs: u v w ua delz va pt
  SUBROUTINE RAYLEIGH_FRICTION_FWD(dt, npx, npy, npz, ks, pm, tau, u, v&
&   , w, pt, ua, va, delz, cp, rg, ptop, hydrostatic, conserve, &
&   rf_cutoff, rf, gridstruct, domain, bd)
    IMPLICIT NONE
!deallocate ( u2f )
    REAL, INTENT(IN) :: dt
! time scale (days)
    REAL, INTENT(IN) :: tau
    REAL, INTENT(IN) :: cp, rg, ptop, rf_cutoff
    INTEGER, INTENT(IN) :: npx, npy, npz, ks
    REAL, DIMENSION(npz), INTENT(IN) :: pm
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: conserve
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
! cell center vertical wind (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temp
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: rf(npz)
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! local:
    REAL :: u2f(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, PARAMETER :: sday=86400.
! scaling velocity  **2
    REAL, PARAMETER :: u000=4900.
    REAL :: rcv
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC LOG
    INTRINSIC SIN
    INTRINSIC SQRT
    LOGICAL :: res
    INTEGER :: ad_count

    u2f = 0.0
    rcv = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    ad_count = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    rcv = 1./(cp-rg)
    IF (.NOT.rf_initialized) THEN
!allocate( rf(npz) )
      rf = 0.0
      res = IS_MASTER()
      IF (res) WRITE(6, *) 'Rayleigh friction E-folding time (days):'
      ad_count = 1
      DO k=1,npz
        IF (pm(k) .LT. rf_cutoff) THEN
          rf(k) = dt/(tau*sday)*SIN(0.5*pi*LOG(rf_cutoff/pm(k))/LOG(&
&           rf_cutoff/ptop))**2
!if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
          kmax = k
          ad_count = ad_count + 1
        ELSE
          GOTO 100
        END IF
      END DO
      CALL PUSHCONTROL(1,0)
      CALL PUSHINTEGER(ad_count)
      CALL PUSHCONTROL(1,0)
      GOTO 110
 100  CALL PUSHCONTROL(1,1)
      CALL PUSHINTEGER(ad_count)
      CALL PUSHCONTROL(1,0)
 110  rf_initialized = .true.
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
!allocate( u2f(isd:ied,jsd:jed,kmax) )
    CALL C2L_ORD2_FWD(u, v, ua, va, gridstruct, npz, gridstruct%&
&               grid_type, bd, gridstruct%nested)
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,u2f,hydrostatic,ua,va,w)
    DO k=1,kmax
      IF (hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        DO j=js,je
          DO i=is,ie
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2 + w(i, j, k)&
&             **2
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
    CALL MPP_UPDATE_DOMAINS(u2f, domain)
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,conserve,hydrostatic,pt,u2f,cp,rg, &
!$OMP                                  ptop,pm,rf,delz,rcv,u,v,w)
    DO k=1,kmax
      IF (conserve) THEN
        IF (hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(pt(i, j, k))
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)/(cp-rg*ptop/&
&               pm(k))*(1.-1./(1.+rf(k)*SQRT(u2f(i, j, k)/u000))**2)
            END DO
          END DO
          CALL PUSHCONTROL(2,2)
        ELSE
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(delz(i, j, k))
              delz(i, j, k) = delz(i, j, k)/pt(i, j, k)
              CALL PUSHREALARRAY(pt(i, j, k))
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)*rcv*(1.-1./(&
&               1.+rf(k)*SQRT(u2f(i, j, k)/u000))**2)
              CALL PUSHREALARRAY(delz(i, j, k))
              delz(i, j, k) = delz(i, j, k)*pt(i, j, k)
            END DO
          END DO
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,0)
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+1
          CALL PUSHREALARRAY(u2f(i, j, k))
          u2f(i, j, k) = rf(k)*SQRT(u2f(i, j, k)/u000)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(u(i, j, k))
          u(i, j, k) = u(i, j, k)/(1.+0.5*(u2f(i, j-1, k)+u2f(i, j, k)))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(v(i, j, k))
          v(i, j, k) = v(i, j, k)/(1.+0.5*(u2f(i-1, j, k)+u2f(i, j, k)))
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(w(i, j, k))
            w(i, j, k) = w(i, j, k)/(1.+u2f(i, j, k))
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
    CALL PUSHREALARRAY(u2f, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(ie)
    CALL PUSHREALARRAY(rcv)
    CALL PUSHINTEGER(js)
  END SUBROUTINE RAYLEIGH_FRICTION_FWD
!  Differentiation of rayleigh_friction in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2
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
!   gradient     of useful results: u v w ua delz va pt
!   with respect to varying inputs: u v w ua delz va pt
  SUBROUTINE RAYLEIGH_FRICTION_BWD(dt, npx, npy, npz, ks, pm, tau, u, &
&   u_ad, v, v_ad, w, w_ad, pt, pt_ad, ua, ua_ad, va, va_ad, delz, &
&   delz_ad, cp, rg, ptop, hydrostatic, conserve, rf_cutoff, rf, &
&   gridstruct, domain, bd)
    IMPLICIT NONE
!deallocate ( u2f )
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: tau
    REAL, INTENT(IN) :: cp, rg, ptop, rf_cutoff
    INTEGER, INTENT(IN) :: npx, npy, npz, ks
    REAL, DIMENSION(npz), INTENT(IN) :: pm
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: conserve
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: w_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ua_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: va_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delz_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: rf(npz)
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL :: u2f(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: u2f_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, PARAMETER :: sday=86400.
    REAL, PARAMETER :: u000=4900.
    REAL :: rcv
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC LOG
    INTRINSIC SIN
    INTRINSIC SQRT
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp2
    REAL :: temp3
    REAL :: temp4
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp8
    REAL :: temp_ad1
    REAL :: temp9
    REAL :: temp_ad2
    REAL :: temp_ad3
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch

    u2f = 0.0
    rcv = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    ad_count = 0
    branch = 0

    CALL POPINTEGER(js)
    CALL POPREALARRAY(rcv)
    CALL POPINTEGER(ie)
    CALL POPINTEGER(is)
    CALL POPINTEGER(je)
    CALL POPREALARRAY(u2f, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
    u2f_ad = 0.0
    DO k=kmax,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(w(i, j, k))
            temp_ad3 = w_ad(i, j, k)/(u2f(i, j, k)+1.)
            u2f_ad(i, j, k) = u2f_ad(i, j, k) - w(i, j, k)*temp_ad3/(u2f&
&             (i, j, k)+1.)
            w_ad(i, j, k) = temp_ad3
          END DO
        END DO
      END IF
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(v(i, j, k))
          temp9 = 0.5*(u2f(i-1, j, k)+u2f(i, j, k)) + 1.
          temp_ad2 = -(v(i, j, k)*0.5*v_ad(i, j, k)/temp9**2)
          u2f_ad(i-1, j, k) = u2f_ad(i-1, j, k) + temp_ad2
          u2f_ad(i, j, k) = u2f_ad(i, j, k) + temp_ad2
          v_ad(i, j, k) = v_ad(i, j, k)/temp9
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(u(i, j, k))
          temp8 = 0.5*(u2f(i, j-1, k)+u2f(i, j, k)) + 1.
          temp_ad1 = -(u(i, j, k)*0.5*u_ad(i, j, k)/temp8**2)
          u2f_ad(i, j-1, k) = u2f_ad(i, j-1, k) + temp_ad1
          u2f_ad(i, j, k) = u2f_ad(i, j, k) + temp_ad1
          u_ad(i, j, k) = u_ad(i, j, k)/temp8
        END DO
      END DO
      DO j=je+1,js-1,-1
        DO i=ie+1,is-1,-1
          CALL POPREALARRAY(u2f(i, j, k))
          IF (u2f(i, j, k)/u000 .EQ. 0.0) THEN
            u2f_ad(i, j, k) = 0.0
          ELSE
            u2f_ad(i, j, k) = rf(k)*u2f_ad(i, j, k)/(2.0*SQRT(u2f(i, j, &
&             k)/u000)*u000)
          END IF
        END DO
      END DO
      CALL POPCONTROL(2,branch)
      IF (branch .NE. 0) THEN
        IF (branch .EQ. 1) THEN
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(delz(i, j, k))
              pt_ad(i, j, k) = pt_ad(i, j, k) + delz(i, j, k)*delz_ad(i&
&               , j, k)
              delz_ad(i, j, k) = pt(i, j, k)*delz_ad(i, j, k)
              CALL POPREALARRAY(pt(i, j, k))
              temp7 = u2f(i, j, k)/u000
              temp6 = SQRT(temp7)
              temp5 = rf(k)*temp6 + 1.
              temp4 = temp5**2
              temp_ad = rcv*0.5*pt_ad(i, j, k)
              IF (temp7 .EQ. 0.0) THEN
                u2f_ad(i, j, k) = u2f_ad(i, j, k) + (1.-1.0/temp4)*&
&                 temp_ad
              ELSE
                u2f_ad(i, j, k) = u2f_ad(i, j, k) + (rf(k)*2*temp5*u2f(i&
&                 , j, k)/(2.0*temp6*temp4**2*u000)-1.0/temp4+1.)*&
&                 temp_ad
              END IF
              CALL POPREALARRAY(delz(i, j, k))
              temp_ad0 = delz_ad(i, j, k)/pt(i, j, k)
              pt_ad(i, j, k) = pt_ad(i, j, k) - delz(i, j, k)*temp_ad0/&
&               pt(i, j, k)
              delz_ad(i, j, k) = temp_ad0
            END DO
          END DO
        ELSE
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(pt(i, j, k))
              temp3 = cp - rg*ptop/pm(k)
              temp2 = u2f(i, j, k)/u000
              temp1 = SQRT(temp2)
              temp0 = rf(k)*temp1 + 1.
              temp = temp0**2
              IF (temp2 .EQ. 0.0) THEN
                u2f_ad(i, j, k) = u2f_ad(i, j, k) + (1.-1.0/temp)*0.5*&
&                 pt_ad(i, j, k)/temp3
              ELSE
                u2f_ad(i, j, k) = u2f_ad(i, j, k) + ((1.-1.0/temp)*0.5/&
&                 temp3+rf(k)*2*temp0*u2f(i, j, k)*0.5/(2.0*temp1*temp**&
&                 2*temp3*u000))*pt_ad(i, j, k)
              END IF
            END DO
          END DO
        END IF
      END IF
    END DO
    CALL MPP_UPDATE_DOMAINS_ADM(u2f, u2f_ad, domain)
    DO k=kmax,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            ua_ad(i, j, k) = ua_ad(i, j, k) + 2*ua(i, j, k)*u2f_ad(i, j&
&             , k)
            va_ad(i, j, k) = va_ad(i, j, k) + 2*va(i, j, k)*u2f_ad(i, j&
&             , k)
            w_ad(i, j, k) = w_ad(i, j, k) + 2*w(i, j, k)*u2f_ad(i, j, k)
            u2f_ad(i, j, k) = 0.0
          END DO
        END DO
      ELSE
        DO j=je,js,-1
          DO i=ie,is,-1
            ua_ad(i, j, k) = ua_ad(i, j, k) + 2*ua(i, j, k)*u2f_ad(i, j&
&             , k)
            va_ad(i, j, k) = va_ad(i, j, k) + 2*va(i, j, k)*u2f_ad(i, j&
&             , k)
            u2f_ad(i, j, k) = 0.0
          END DO
        END DO
      END IF
    END DO
    CALL C2L_ORD2_BWD(u, u_ad, v, v_ad, ua, ua_ad, va, va_ad, gridstruct&
&               , npz, gridstruct%grid_type, bd, gridstruct%nested)
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(ad_count)
      DO i0=1,ad_count
        IF (i0 .EQ. 1) CALL POPCONTROL(1,branch)
      END DO
    END IF
  END SUBROUTINE RAYLEIGH_FRICTION_BWD
  SUBROUTINE RAYLEIGH_FRICTION(dt, npx, npy, npz, ks, pm, tau, u, v, w, &
&   pt, ua, va, delz, cp, rg, ptop, hydrostatic, conserve, rf_cutoff, rf&
&   , gridstruct, domain, bd)
    IMPLICIT NONE
!deallocate ( u2f )
    REAL, INTENT(IN) :: dt
! time scale (days)
    REAL, INTENT(IN) :: tau
    REAL, INTENT(IN) :: cp, rg, ptop, rf_cutoff
    INTEGER, INTENT(IN) :: npx, npy, npz, ks
    REAL, DIMENSION(npz), INTENT(IN) :: pm
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: conserve
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
! cell center vertical wind (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temp
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: rf(npz)
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! local:
    REAL :: u2f(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, PARAMETER :: sday=86400.
! scaling velocity  **2
    REAL, PARAMETER :: u000=4900.
    REAL :: rcv
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC LOG
    INTRINSIC SIN
    INTRINSIC SQRT
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    rcv = 1./(cp-rg)
    IF (.NOT.rf_initialized) THEN
!allocate( rf(npz) )
      rf = 0.0
      IF (IS_MASTER()) WRITE(6, *) &
&                      'Rayleigh friction E-folding time (days):'
      DO k=1,npz
        IF (pm(k) .LT. rf_cutoff) THEN
          rf(k) = dt/(tau*sday)*SIN(0.5*pi*LOG(rf_cutoff/pm(k))/LOG(&
&           rf_cutoff/ptop))**2
!if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
          kmax = k
        ELSE
          GOTO 100
        END IF
      END DO
 100  rf_initialized = .true.
    END IF
!allocate( u2f(isd:ied,jsd:jed,kmax) )
    CALL C2L_ORD2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, &
&           bd, gridstruct%nested)
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,u2f,hydrostatic,ua,va,w)
    DO k=1,kmax
      IF (hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2 + w(i, j, k)&
&             **2
          END DO
        END DO
      END IF
    END DO
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_DOMAINS(u2f, domain)
    CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,conserve,hydrostatic,pt,u2f,cp,rg, &
!$OMP                                  ptop,pm,rf,delz,rcv,u,v,w)
    DO k=1,kmax
      IF (conserve) THEN
        IF (hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)/(cp-rg*ptop/&
&               pm(k))*(1.-1./(1.+rf(k)*SQRT(u2f(i, j, k)/u000))**2)
            END DO
          END DO
        ELSE
          DO j=js,je
            DO i=is,ie
              delz(i, j, k) = delz(i, j, k)/pt(i, j, k)
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)*rcv*(1.-1./(&
&               1.+rf(k)*SQRT(u2f(i, j, k)/u000))**2)
              delz(i, j, k) = delz(i, j, k)*pt(i, j, k)
            END DO
          END DO
        END IF
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+1
          u2f(i, j, k) = rf(k)*SQRT(u2f(i, j, k)/u000)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          u(i, j, k) = u(i, j, k)/(1.+0.5*(u2f(i, j-1, k)+u2f(i, j, k)))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v(i, j, k) = v(i, j, k)/(1.+0.5*(u2f(i-1, j, k)+u2f(i, j, k)))
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            w(i, j, k) = w(i, j, k)/(1.+u2f(i, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE RAYLEIGH_FRICTION
!  Differentiation of compute_aam in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: u v delp ua aam va m_fac ps
!   with respect to varying inputs: u v delp ua aam va m_fac ps
  SUBROUTINE COMPUTE_AAM_FWD(npz, is, ie, js, je, isd, ied, jsd, jed, &
&   gridstruct, bd, ptop, ua, va, u, v, delp, aam, ps, m_fac)
    IMPLICIT NONE
! Compute vertically (mass) integrated Atmospheric Angular Momentum
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: is, ie, js, je
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    REAL, INTENT(IN) :: ptop
! D grid zonal wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua, va
    REAL :: aam(is:ie, js:je)
    REAL :: m_fac(is:ie, js:je)
    REAL :: ps(isd:ied, jsd:jed)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
! local:
    REAL, DIMENSION(is:ie) :: r1, r2, dm
    INTEGER :: i, j, k
    INTRINSIC COS

    r1 = 0.0
    r2 = 0.0
    dm = 0.0

    CALL C2L_ORD2_FWD(u, v, ua, va, gridstruct, npz, gridstruct%&
&               grid_type, bd, gridstruct%nested)
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gridstruct,aam,m_fac,ps,ptop,delp,agrav,ua) &
!$OMP                          private(r1, r2, dm)
    DO j=js,je
      DO i=is,ie
        CALL PUSHREALARRAY(r1(i))
        r1(i) = radius*COS(gridstruct%agrid(i, j, 2))
        CALL PUSHREALARRAY(r2(i))
        r2(i) = r1(i)*r1(i)
        aam(i, j) = 0.
        m_fac(i, j) = 0.
        ps(i, j) = ptop
      END DO
      DO k=1,npz
        DO i=is,ie
          CALL PUSHREALARRAY(dm(i))
          dm(i) = delp(i, j, k)
          ps(i, j) = ps(i, j) + dm(i)
          dm(i) = dm(i)*agrav
          aam(i, j) = aam(i, j) + (r2(i)*omega+r1(i)*ua(i, j, k))*dm(i)
          m_fac(i, j) = m_fac(i, j) + dm(i)*r2(i)
        END DO
      END DO
    END DO
    CALL PUSHREALARRAY(r2, ie - is + 1)
    CALL PUSHREALARRAY(r1, ie - is + 1)
    CALL PUSHREALARRAY(dm, ie - is + 1)
  END SUBROUTINE COMPUTE_AAM_FWD
!  Differentiation of compute_aam in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: u v delp ua aam va m_fac ps
!   with respect to varying inputs: u v delp ua aam va m_fac ps
  SUBROUTINE COMPUTE_AAM_BWD(npz, is, ie, js, je, isd, ied, jsd, jed, &
&   gridstruct, bd, ptop, ua, ua_ad, va, va_ad, u, u_ad, v, v_ad, delp, &
&   delp_ad, aam, aam_ad, ps, ps_ad, m_fac, m_fac_ad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: is, ie, js, je
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    REAL, INTENT(IN) :: ptop
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: u_ad(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: v_ad(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_ad(isd:ied, jsd:jed, npz)
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua_ad, &
&   va_ad
    REAL :: aam(is:ie, js:je)
    REAL :: aam_ad(is:ie, js:je)
    REAL :: m_fac(is:ie, js:je)
    REAL :: m_fac_ad(is:ie, js:je)
    REAL :: ps(isd:ied, jsd:jed)
    REAL :: ps_ad(isd:ied, jsd:jed)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    REAL, DIMENSION(is:ie) :: r1, r2, dm
    REAL, DIMENSION(is:ie) :: dm_ad
    INTEGER :: i, j, k
    INTRINSIC COS

    r1 = 0.0
    r2 = 0.0
    dm = 0.0

    CALL POPREALARRAY(dm, ie - is + 1)
    CALL POPREALARRAY(r1, ie - is + 1)
    CALL POPREALARRAY(r2, ie - is + 1)
    dm_ad = 0.0
    DO j=je,js,-1
      DO k=npz,1,-1
        DO i=ie,is,-1
          dm_ad(i) = dm_ad(i) + (r2(i)*omega+r1(i)*ua(i, j, k))*aam_ad(i&
&           , j) + r2(i)*m_fac_ad(i, j)
          ua_ad(i, j, k) = ua_ad(i, j, k) + dm(i)*r1(i)*aam_ad(i, j)
          dm_ad(i) = ps_ad(i, j) + agrav*dm_ad(i)
          CALL POPREALARRAY(dm(i))
          delp_ad(i, j, k) = delp_ad(i, j, k) + dm_ad(i)
          dm_ad(i) = 0.0
        END DO
      END DO
      DO i=ie,is,-1
        ps_ad(i, j) = 0.0
        m_fac_ad(i, j) = 0.0
        aam_ad(i, j) = 0.0
        CALL POPREALARRAY(r2(i))
        CALL POPREALARRAY(r1(i))
      END DO
    END DO
    CALL C2L_ORD2_BWD(u, u_ad, v, v_ad, ua, ua_ad, va, va_ad, gridstruct&
&               , npz, gridstruct%grid_type, bd, gridstruct%nested)
  END SUBROUTINE COMPUTE_AAM_BWD
  SUBROUTINE COMPUTE_AAM(npz, is, ie, js, je, isd, ied, jsd, jed, &
&   gridstruct, bd, ptop, ua, va, u, v, delp, aam, ps, m_fac)
    IMPLICIT NONE
! Compute vertically (mass) integrated Atmospheric Angular Momentum
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: is, ie, js, je
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    REAL, INTENT(IN) :: ptop
! D grid zonal wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua, va
    REAL, INTENT(OUT) :: aam(is:ie, js:je)
    REAL, INTENT(OUT) :: m_fac(is:ie, js:je)
    REAL, INTENT(OUT) :: ps(isd:ied, jsd:jed)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
! local:
    REAL, DIMENSION(is:ie) :: r1, r2, dm
    INTEGER :: i, j, k
    INTRINSIC COS
    CALL C2L_ORD2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, &
&           bd, gridstruct%nested)
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gridstruct,aam,m_fac,ps,ptop,delp,agrav,ua) &
!$OMP                          private(r1, r2, dm)
    DO j=js,je
      DO i=is,ie
        r1(i) = radius*COS(gridstruct%agrid(i, j, 2))
        r2(i) = r1(i)*r1(i)
        aam(i, j) = 0.
        m_fac(i, j) = 0.
        ps(i, j) = ptop
      END DO
      DO k=1,npz
        DO i=is,ie
          dm(i) = delp(i, j, k)
          ps(i, j) = ps(i, j) + dm(i)
          dm(i) = dm(i)*agrav
          aam(i, j) = aam(i, j) + (r2(i)*omega+r1(i)*ua(i, j, k))*dm(i)
          m_fac(i, j) = m_fac(i, j) + dm(i)*r2(i)
        END DO
      END DO
    END DO
  END SUBROUTINE COMPUTE_AAM
end module fv_dynamics_adm_mod
