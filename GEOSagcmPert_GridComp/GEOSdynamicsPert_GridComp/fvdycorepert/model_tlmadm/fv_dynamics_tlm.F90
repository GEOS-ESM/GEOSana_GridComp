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
module fv_dynamics_tlm_mod
   use constants_mod,           only: grav, pi=>pi_8, radius, hlv, rdgas, omega, rvgas, cp_vapor
   use dyn_core_tlm_mod,        only: dyn_core, del2_cubed, init_ijk_mem
   use dyn_core_tlm_mod,        only: dyn_core_tlm, del2_cubed_tlm
   use fv_mapz_tlm_mod,         only: compute_total_energy, Lagrangian_to_Eulerian, moist_cv, moist_cp
   use fv_mapz_tlm_mod,         only: compute_total_energy_tlm, Lagrangian_to_Eulerian_tlm
   use fv_tracer2d_tlm_mod,     only: tracer_2d, tracer_2d_1L, tracer_2d_nested
   use fv_tracer2d_tlm_mod,     only: tracer_2d_tlm, tracer_2d_1L_tlm, tracer_2d_nested_tlm
   use fv_grid_utils_tlm_mod,   only: cubed_to_latlon, c2l_ord2, g_sum
   use fv_grid_utils_tlm_mod,   only: c2l_ord2_tlm, g_sum_tlm
   use fv_fill_mod,             only: fill2D
   use fv_mp_mod,               only: is_master
   use fv_mp_mod,               only: group_halo_update_type
   use fv_mp_tlm_mod,           only: start_group_halo_update, complete_group_halo_update
   use fv_mp_tlm_mod,           only: start_group_halo_update_tlm
   use fv_timing_mod,           only: timing_on, timing_off
   use diag_manager_mod,        only: send_data
   use fv_diagnostics_mod,      only: fv_time, prt_mxm, range_check, prt_minmax
   use mpp_domains_mod,         only: mpp_update_domains, DGRID_NE, CGRID_NE, domain2D
   use fv_mp_tlm_mod,           only: mpp_update_domains_tlm
   use mpp_mod,                 only: mpp_pe
   use field_manager_mod,       only: MODEL_ATMOS
   use tracer_manager_mod,      only: get_tracer_index
   use fv_sg_mod,               only: neg_adj3
   use fv_nesting_tlm_mod,      only: setup_nested_grid_BCs
   use fv_nesting_tlm_mod,      only: setup_nested_grid_BCs_tlm
   use boundary_tlm_mod,        only: nested_grid_BC_apply_intT
   use boundary_tlm_mod,        only: nested_grid_BC_apply_intT_tlm
   use fv_arrays_mod,           only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_grid_bounds_type
   use fv_nwp_nudge_mod,        only: do_adiabatic_init
#ifdef MAPL_MODE
   use fv_control_mod,          only: dyn_timer, comm_timer
#endif
   use fv_arrays_nlm_mod,       only: fv_flags_pert_type, fpp

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
public :: fv_dynamics, fv_dynamics_tlm

!---- version number -----
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of fv_dynamics in forward (tangent) mode:
!   variations   of useful results: q u v w delp delz pt
!   with respect to varying inputs: q u v w delp delz pt
!   RW status of diff variables: peln:(loc) q:in-out u:in-out v:in-out
!                w:in-out delp:in-out ua:(loc) uc:(loc) mfx:(loc)
!                delz:in-out mfy:(loc) omga:(loc) va:(loc) vc:(loc)
!                pkz:(loc) pe:(loc) pk:(loc) ps:(loc) pt:in-out
!                cx:(loc) cy:(loc)
!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
  SUBROUTINE FV_DYNAMICS_TLM(npx, npy, npz, nq_tot, ng, bdt, consv_te, &
&   fill, reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split, &
&   q_split, u, u_tl, v, v_tl, w, w_tl, delz, delz_tl, hydrostatic, pt, &
&   pt_tl, delp, delp_tl, q, q_tl, ps, ps_tl, pe, pe_tl, pk, pk_tl, peln&
&   , peln_tl, pkz, pkz_tl, phis, q_con, omga, omga_tl, ua, ua_tl, va, &
&   va_tl, uc, uc_tl, vc, vc_tl, ak, bk, mfx, mfx_tl, mfy, mfy_tl, cx, &
&   cx_tl, cy, cy_tl, ze0, hybrid_z, gridstruct, flagstruct, flagstructp&
&   , neststruct, idiag, bd, parent_grid, domain, time_total)
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
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u_tl
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v_tl
!  W (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: w_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! specific humidity and constituents
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst)
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst&
&   )
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delz_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
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
    REAL, INTENT(INOUT) :: ps_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
! edge pressure (pascal)
    REAL, INTENT(INOUT) :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(INOUT) :: pe_tl(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1&
&   )
! pe**kappa
    REAL, INTENT(INOUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL, INTENT(INOUT) :: pk_tl(bd%is:bd%ie, bd%js:bd%je, npz+1)
! ln(pe)
    REAL, INTENT(INOUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: peln_tl(bd%is:bd%ie, npz+1, bd%js:bd%je)
! finite-volume mean pk
    REAL, INTENT(INOUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: pkz_tl(bd%is:bd%ie, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: q_con(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: omga_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! (uc,vc) mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: uc_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: vc_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua_tl, va_tl
    REAL, DIMENSION(npz+1), INTENT(IN) :: ak, bk
! Accumulated Mass flux arrays: the "Flux Capacitor"
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_tl(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_tl(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(INOUT) :: flagstruct
    TYPE(FV_FLAGS_PERT_TYPE), INTENT(INOUT) :: flagstructp
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(FV_DIAG_TYPE), INTENT(IN) :: idiag
! Local Arrays
    REAL :: ws(bd%is:bd%ie, bd%js:bd%je)
    REAL :: ws_tl(bd%is:bd%ie, bd%js:bd%je)
    REAL :: te_2d(bd%is:bd%ie, bd%js:bd%je)
    REAL :: te_2d_tl(bd%is:bd%ie, bd%js:bd%je)
    REAL :: teq(bd%is:bd%ie, bd%js:bd%je)
    REAL :: teq_tl(bd%is:bd%ie, bd%js:bd%je)
    REAL :: ps2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: ps2_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: m_fac(bd%is:bd%ie, bd%js:bd%je)
    REAL :: m_fac_tl(bd%is:bd%ie, bd%js:bd%je)
    REAL :: pfull(npz)
    REAL, DIMENSION(bd%is:bd%ie) :: cvm
    REAL :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz), dtdt_m(bd%is:bd%ie, &
&   bd%js:bd%je, npz), cappa(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: dp1_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
#ifdef OVERLOAD_R4
    REAL(kind=4) :: psx(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL(kind=4) :: psx_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
#else
    REAL(kind=8) :: psx(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL(kind=8) :: psx_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
#endif
    REAL(kind=8) :: dpx(bd%is:bd%ie, bd%js:bd%je)
    REAL(kind=8) :: dpx_tl(bd%is:bd%ie, bd%js:bd%je)
    REAL :: akap, rdg, ph1, ph2, mdt, gam, amdt, u0
    REAL :: amdt_tl, u0_tl
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
    REAL :: gz_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pkc(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pkc_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: ptc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: ptc_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: crx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: cry(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cry_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: divgd(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL :: divgd_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL :: delpc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: delpc_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: ut(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: ut_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: vt_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: zh(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: zh_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pk3(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: pk3_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL :: du_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL :: dv_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    INTRINSIC ANY
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC ABS
    INTRINSIC REAL
    INTRINSIC COS
    REAL :: abs0
    REAL :: abs1
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: arg2
    REAL :: arg2_tl
    REAL :: result1
    REAL :: result1_tl
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
      CALL SETUP_NESTED_GRID_BCS_TLM(npx, npy, npz, zvir, ncnst, u, u_tl&
&                              , v, v_tl, w, pt, delp, delz, q, uc, &
&                              uc_tl, vc, vc_tl, pkz, neststruct%nested&
&                              , flagstruct%inline_q, flagstruct%make_nh&
&                              , ng, gridstruct, flagstruct, neststruct&
&                              , neststruct%nest_timestep, neststruct%&
&                              tracer_nest_timestep, domain, bd, nwat)
      IF (gridstruct%nested) CALL NESTED_GRID_BC_APPLY_INTT_TLM(pt, &
&                                                         pt_tl, 0, 0, &
&                                                         npx, npy, npz&
&                                                         , bd, 1., 1., &
&                                                         neststruct%&
&                                                         pt_bc, bctype=&
&                                                         neststruct%&
&                                                         nestbctype)
!Correct halo values have now been set up for BCs; we can go ahead and apply them too...
      CALL TIMING_OFF('NEST_BCs')
    ELSE
      uc_tl = 0.0
      vc_tl = 0.0
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
      dp1_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,zvir,nwat,q,q_con,sphum,liq_wat, &
!$OMP      rainwat,ice_wat,snowwat,graupel) private(cvm)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dp1_tl(i, j, k) = zvir*q_tl(i, j, k, sphum)
            dp1(i, j, k) = zvir*q(i, j, k, sphum)
          END DO
        END DO
      END DO
    ELSE
      dp1_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,zvir,q,q_con,sphum,liq_wat, &
!$OMP                                  rainwat,ice_wat,snowwat,graupel,pkz,flagstruct, &
!$OMP                                  cappa,kappa,rdg,delp,pt,delz,nwat)              &
!$OMP                          private(cvm)
      DO k=1,npz
        IF (flagstruct%moist_phys) THEN
          DO j=js,je
            DO i=is,ie
              dp1_tl(i, j, k) = zvir*q_tl(i, j, k, sphum)
              dp1(i, j, k) = zvir*q(i, j, k, sphum)
              arg1_tl = (rdg*((delp_tl(i, j, k)*pt(i, j, k)+delp(i, j, k&
&               )*pt_tl(i, j, k))*(1.+dp1(i, j, k))+delp(i, j, k)*pt(i, &
&               j, k)*dp1_tl(i, j, k))*delz(i, j, k)-rdg*delp(i, j, k)*&
&               pt(i, j, k)*(1.+dp1(i, j, k))*delz_tl(i, j, k))/delz(i, &
&               j, k)**2
              arg1 = rdg*delp(i, j, k)*pt(i, j, k)*(1.+dp1(i, j, k))/&
&               delz(i, j, k)
              arg2_tl = kappa*arg1_tl/arg1
              arg2 = kappa*LOG(arg1)
              pkz_tl(i, j, k) = arg2_tl*EXP(arg2)
              pkz(i, j, k) = EXP(arg2)
            END DO
          END DO
        ELSE
! Using dry pressure for the definition of the virtual potential temperature
!              pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!                                      (1.-q(i,j,k,sphum))/delz(i,j,k)) )
          DO j=js,je
            DO i=is,ie
              dp1_tl(i, j, k) = 0.0
              dp1(i, j, k) = 0.
              arg1_tl = (rdg*(delp_tl(i, j, k)*pt(i, j, k)+delp(i, j, k)&
&               *pt_tl(i, j, k))*delz(i, j, k)-rdg*delp(i, j, k)*pt(i, j&
&               , k)*delz_tl(i, j, k))/delz(i, j, k)**2
              arg1 = rdg*delp(i, j, k)*pt(i, j, k)/delz(i, j, k)
              arg2_tl = kappa*arg1_tl/arg1
              arg2 = kappa*LOG(arg1)
              pkz_tl(i, j, k) = arg2_tl*EXP(arg2)
              pkz(i, j, k) = EXP(arg2)
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
      CALL PRT_MXM('pk_b', pk, is, ie, js, je, 0, npz + 1, 1., &
&            gridstruct%area_64, domain)
      CALL PRT_MXM('pkz_b', pkz, is, ie, js, je, 0, npz, 1., gridstruct%&
&            area_64, domain)
    END IF
!---------------------
! Compute Total Energy
!---------------------
    IF (consv_te .GT. 0. .AND. (.NOT.do_adiabatic_init)) THEN
      CALL COMPUTE_TOTAL_ENERGY_TLM(is, ie, js, je, isd, ied, jsd, jed, &
&                             npz, u, u_tl, v, v_tl, w, w_tl, delz, &
&                             delz_tl, pt, pt_tl, delp, delp_tl, q, q_tl&
&                             , dp1, dp1_tl, pe, pe_tl, peln, peln_tl, &
&                             phis, gridstruct%rsin2, gridstruct%cosa_s&
&                             , zvir, cp_air, rdgas, hlv, te_2d, &
&                             te_2d_tl, ua, va, teq, teq_tl, flagstruct%&
&                             moist_phys, nwat, sphum, liq_wat, rainwat&
&                             , ice_wat, snowwat, graupel, hydrostatic, &
&                             idiag%id_te)
      IF (idiag%id_te .GT. 0) used = SEND_DATA(idiag%id_te, teq, fv_time&
&         )
!              te_den=1.E-9*g_sum(teq, is, ie, js, je, ng, area, 0)/(grav*4.*pi*radius**2)
!              if(is_master())  write(*,*) 'Total Energy Density (Giga J/m**2)=',te_den
    ELSE
      teq_tl = 0.0
      te_2d_tl = 0.0
    END IF
    IF ((flagstruct%consv_am .OR. idiag%id_amdt .GT. 0) .AND. (.NOT.&
&       do_adiabatic_init)) THEN
      m_fac_tl = 0.0
      ps2_tl = 0.0
      va_tl = 0.0
      ua_tl = 0.0
      CALL COMPUTE_AAM_TLM(npz, is, ie, js, je, isd, ied, jsd, jed, &
&                    gridstruct, bd, ptop, ua, ua_tl, va, va_tl, u, u_tl&
&                    , v, v_tl, delp, delp_tl, teq, teq_tl, ps2, ps2_tl&
&                    , m_fac, m_fac_tl)
    ELSE
      ua_tl = 0.0
      va_tl = 0.0
      ps2_tl = 0.0
      m_fac_tl = 0.0
    END IF
    IF (flagstruct%tau .GT. 0.) THEN
      IF (gridstruct%grid_type .LT. 4) THEN
        IF (bdt .GE. 0.) THEN
          abs0 = bdt
        ELSE
          abs0 = -bdt
        END IF
        CALL RAYLEIGH_SUPER_TLM(abs0, npx, npy, npz, ks, pfull, phis, &
&                         flagstruct%tau, u, u_tl, v, v_tl, w, w_tl, pt&
&                         , pt_tl, ua, ua_tl, va, va_tl, delz, &
&                         gridstruct%agrid, cp_air, rdgas, ptop, &
&                         hydrostatic, .NOT.neststruct%nested, &
&                         flagstruct%rf_cutoff, rf, gridstruct, domain, &
&                         bd)
      ELSE
        IF (bdt .GE. 0.) THEN
          abs1 = bdt
        ELSE
          abs1 = -bdt
        END IF
        CALL RAYLEIGH_FRICTION_TLM(abs1, npx, npy, npz, ks, pfull, &
&                            flagstruct%tau, u, u_tl, v, v_tl, w, w_tl, &
&                            pt, pt_tl, ua, ua_tl, va, va_tl, delz, &
&                            delz_tl, cp_air, rdgas, ptop, hydrostatic, &
&                            .true., flagstruct%rf_cutoff, rf, &
&                            gridstruct, domain, bd)
      END IF
    END IF
! Convert pt to virtual potential temperature on the first timestep
    IF (flagstruct%adiabatic) THEN
!$OMP parallel do default(none) shared(theta_d,is,ie,js,je,npz,pt,pkz,q)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            pt_tl(i, j, k) = (pt_tl(i, j, k)*pkz(i, j, k)-pt(i, j, k)*&
&             pkz_tl(i, j, k))/pkz(i, j, k)**2
            pt(i, j, k) = pt(i, j, k)/pkz(i, j, k)
          END DO
        END DO
        IF (theta_d .GT. 0) THEN
          DO j=js,je
            DO i=is,ie
              q_tl(i, j, k, theta_d) = pt_tl(i, j, k)
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
            pt_tl(i, j, k) = ((pt_tl(i, j, k)*(1.+dp1(i, j, k))+pt(i, j&
&             , k)*dp1_tl(i, j, k))*pkz(i, j, k)-pt(i, j, k)*(1.+dp1(i, &
&             j, k))*pkz_tl(i, j, k))/pkz(i, j, k)**2
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
      psx_tl = 0.0_8
      DO j=js,je
        DO i=is,ie
          psx_tl(i, j) = pe_tl(i, npz+1, j)
          psx(i, j) = pe(i, npz+1, j)
          dpx(i, j) = 0.0
        END DO
      END DO
    ELSE
      psx_tl = 0.0_8
    END IF
    CALL TIMING_ON('FV_DYN_LOOP')
    omga_tl = 0.0
    ps_tl = 0.0
    pk3_tl = 0.0
    xfx_tl = 0.0
    ws_tl = 0.0
    gz_tl = 0.0
    du_tl = 0.0
    dv_tl = 0.0
    ptc_tl = 0.0
    ut_tl = 0.0
    divgd_tl = 0.0
    pkc_tl = 0.0
    delpc_tl = 0.0
    yfx_tl = 0.0
    vt_tl = 0.0
    zh_tl = 0.0
    dpx_tl = 0.0_8
    crx_tl = 0.0
    cry_tl = 0.0
! first level of time-split
    DO n_map=1,k_split
      CALL TIMING_ON('COMM_TOTAL')
      CALL START_GROUP_HALO_UPDATE_TLM(i_pack(1), delp, delp_tl, domain&
&                                , complete=.true.)
      CALL START_GROUP_HALO_UPDATE_TLM(i_pack(1), pt, pt_tl, domain, &
&                                complete=.true.)
      CALL START_GROUP_HALO_UPDATE_TLM(i_pack(8), u, u_tl, v, v_tl, &
&                                domain, gridtype=dgrid_ne)
      CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,dp1,delp)
      DO k=1,npz
        DO j=jsd,jed
          DO i=isd,ied
            dp1_tl(i, j, k) = delp_tl(i, j, k)
            dp1(i, j, k) = delp(i, j, k)
          END DO
        END DO
      END DO
      IF (n_map .EQ. k_split) last_step = .true.
      CALL TIMING_ON('DYN_CORE')
      CALL DYN_CORE_TLM(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir&
&                 , cp_air, akap, cappa, grav, hydrostatic, u, u_tl, v, &
&                 v_tl, w, w_tl, delz, delz_tl, pt, pt_tl, q, q_tl, delp&
&                 , delp_tl, pe, pe_tl, pk, pk_tl, phis, ws, ws_tl, omga&
&                 , omga_tl, ptop, pfull, ua, ua_tl, va, va_tl, uc, &
&                 uc_tl, vc, vc_tl, mfx, mfx_tl, mfy, mfy_tl, cx, cx_tl&
&                 , cy, cy_tl, pkz, pkz_tl, peln, peln_tl, q_con, ak, bk&
&                 , dpx, dpx_tl, ks, gridstruct, flagstruct, flagstructp&
&                 , neststruct, idiag, bd, domain, n_map .EQ. 1, i_pack&
&                 , last_step, gz, gz_tl, pkc, pkc_tl, ptc, ptc_tl, crx&
&                 , crx_tl, xfx, xfx_tl, cry, cry_tl, yfx, yfx_tl, divgd&
&                 , divgd_tl, delpc, delpc_tl, ut, ut_tl, vt, vt_tl, zh&
&                 , zh_tl, pk3, pk3_tl, du, du_tl, dv, dv_tl, time_total&
&                )
      CALL TIMING_OFF('DYN_CORE')
!DryMassRoundoffControl
      IF (last_step) THEN
        IF (fpp%fpp_overload_r4) THEN
          DO j=js,je
            DO i=is,ie
              psx_tl(i, j) = psx_tl(i, j) + dpx_tl(i, j)
              psx(i, j) = psx(i, j) + dpx(i, j)
            END DO
          END DO
          CALL TIMING_ON('COMM_TOTAL')
          CALL MPP_UPDATE_DOMAINS_TLM(psx, psx_tl, domain)
          CALL TIMING_OFF('COMM_TOTAL')
          DO j=js-1,je+1
            DO i=is-1,ie+1
              pe_tl(i, npz+1, j) = psx_tl(i, j)
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
          CALL TRACER_2D_NESTED_TLM(q, q_tl, dp1, dp1_tl, mfx, mfx_tl, &
&                             mfy, mfy_tl, cx, cx_tl, cy, cy_tl, &
&                             gridstruct, bd, domain, npx, npy, npz, nq&
&                             , flagstruct%hord_tr, q_split, mdt, idiag%&
&                             id_divg, i_pack(10), flagstruct%nord_tr, &
&                             flagstruct%trdm2, k_split, neststruct, &
&                             parent_grid, flagstructp%hord_tr_pert, &
&                             flagstructp%nord_tr_pert, flagstructp%&
&                             trdm2_pert, flagstructp%split_damp_tr)
        ELSE IF (flagstruct%z_tracer) THEN
          CALL TRACER_2D_1L_TLM(q, q_tl, dp1, dp1_tl, mfx, mfx_tl, mfy, &
&                         mfy_tl, cx, cx_tl, cy, cy_tl, gridstruct, bd, &
&                         domain, npx, npy, npz, nq, flagstruct%hord_tr&
&                         , q_split, mdt, idiag%id_divg, i_pack(10), &
&                         flagstruct%nord_tr, flagstruct%trdm2, &
&                         flagstructp%hord_tr_pert, flagstructp%&
&                         nord_tr_pert, flagstructp%trdm2_pert, &
&                         flagstructp%split_damp_tr)
        ELSE
          CALL TRACER_2D_TLM(q, q_tl, dp1, dp1_tl, mfx, mfx_tl, mfy, &
&                      mfy_tl, cx, cx_tl, cy, cy_tl, gridstruct, bd, &
&                      domain, npx, npy, npz, nq, flagstruct%hord_tr, &
&                      q_split, mdt, idiag%id_divg, i_pack(10), &
&                      flagstruct%nord_tr, flagstruct%trdm2, flagstructp&
&                      %hord_tr_pert, flagstructp%nord_tr_pert, &
&                      flagstructp%trdm2_pert, flagstructp%split_damp_tr&
&                     )
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
        CALL LAGRANGIAN_TO_EULERIAN_TLM(last_step, consv_te, ps, ps_tl, &
&                                 pe, pe_tl, delp, delp_tl, pkz, pkz_tl&
&                                 , pk, pk_tl, mdt, bdt, npz, is, ie, js&
&                                 , je, isd, ied, jsd, jed, nq, nwat, &
&                                 sphum, q_con, u, u_tl, v, v_tl, w, &
&                                 w_tl, delz, delz_tl, pt, pt_tl, q, &
&                                 q_tl, phis, zvir, cp_air, akap, cappa&
&                                 , kord_mt, kord_wz, kord_tracer, &
&                                 kord_tm, peln, peln_tl, te_2d, &
&                                 te_2d_tl, ng, ua, ua_tl, va, omga, &
&                                 omga_tl, dp1, dp1_tl, ws, ws_tl, fill&
&                                 , reproduce_sum, idiag%id_mdt .GT. 0, &
&                                 dtdt_m, ptop, ak, bk, pfull, &
&                                 flagstruct, gridstruct, domain, &
&                                 flagstruct%do_sat_adj, hydrostatic, &
&                                 hybrid_z, do_omega, flagstruct%&
&                                 adiabatic, do_adiabatic_init, mfx, mfy&
&                                 , flagstruct%remap_option, &
&                                 kord_mt_pert, kord_wz_pert, &
&                                 kord_tracer_pert, kord_tm_pert)
        CALL TIMING_OFF('Remapping')
        IF (last_step) THEN
          IF (.NOT.hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga,delp,delz,w)
            DO k=1,npz
              DO j=js,je
                DO i=is,ie
                  omga_tl(i, j, k) = (delp_tl(i, j, k)*delz(i, j, k)-&
&                   delp(i, j, k)*delz_tl(i, j, k))*w(i, j, k)/delz(i, j&
&                   , k)**2 + delp(i, j, k)*w_tl(i, j, k)/delz(i, j, k)
                  omga(i, j, k) = delp(i, j, k)/delz(i, j, k)*w(i, j, k)
                END DO
              END DO
            END DO
          END IF
!--------------------------
! Filter omega for physics:
!--------------------------
          IF (flagstruct%nf_omega .GT. 0) CALL DEL2_CUBED_TLM(omga, &
&                                                       omga_tl, 0.18*&
&                                                       gridstruct%&
&                                                       da_min, &
&                                                       gridstruct, &
&                                                       domain, npx, npy&
&                                                       , npz, &
&                                                       flagstruct%&
&                                                       nf_omega, bd)
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
      CALL COMPUTE_AAM_TLM(npz, is, ie, js, je, isd, ied, jsd, jed, &
&                    gridstruct, bd, ptop, ua, ua_tl, va, va_tl, u, u_tl&
&                    , v, v_tl, delp, delp_tl, te_2d, te_2d_tl, ps, &
&                    ps_tl, m_fac, m_fac_tl)
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
          te_2d_tl(i, j) = te_2d_tl(i, j) - teq_tl(i, j) + dt2*idiag%zxg&
&           (i, j)*(ps2_tl(i, j)+ps_tl(i, j))
          te_2d(i, j) = te_2d(i, j) - teq(i, j) + dt2*(ps2(i, j)+ps(i, j&
&           ))*idiag%zxg(i, j)
        END DO
      END DO
      IF (idiag%id_amdt .GT. 0) used = SEND_DATA(idiag%id_amdt, te_2d/&
&         bdt, fv_time)
      IF (flagstruct%consv_am .OR. prt_minmax) THEN
        amdt_tl = G_SUM_TLM(domain, te_2d, te_2d_tl, is, ie, js, je, ng&
&         , gridstruct%area_64, 0, reproduce=.true., g_sum=amdt)
        result1_tl = G_SUM_TLM(domain, m_fac, m_fac_tl, is, ie, js, je, &
&         ng, gridstruct%area_64, 0, reproduce=.true., g_sum=result1)
        u0_tl = -((radius*amdt_tl*result1-radius*amdt*result1_tl)/&
&         result1**2)
        u0 = -(radius*amdt/result1)
        IF (IS_MASTER() .AND. prt_minmax) WRITE(6, *) &
&                                       'Dynamic AM tendency (Hadleys)='&
&                                         , amdt/(bdt*1.e18), &
&                                         'del-u (per day)=', u0*86400./&
&                                         bdt
      ELSE
        u0_tl = 0.0
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
              u_tl(i, j, k) = u_tl(i, j, k) + gridstruct%l2c_u(i, j)*&
&               u0_tl
              u(i, j, k) = u(i, j, k) + u0*gridstruct%l2c_u(i, j)
            END DO
          END DO
          DO j=js,je
            DO i=is,ie+1
              v_tl(i, j, k) = v_tl(i, j, k) + gridstruct%l2c_v(i, j)*&
&               u0_tl
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
  END SUBROUTINE FV_DYNAMICS_TLM
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
    REAL :: arg1
    REAL :: arg2
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
              arg1 = rdg*delp(i, j, k)*pt(i, j, k)*(1.+dp1(i, j, k))/&
&               delz(i, j, k)
              arg2 = kappa*LOG(arg1)
              pkz(i, j, k) = EXP(arg2)
            END DO
          END DO
        ELSE
! Using dry pressure for the definition of the virtual potential temperature
!              pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!                                      (1.-q(i,j,k,sphum))/delz(i,j,k)) )
          DO j=js,je
            DO i=is,ie
              dp1(i, j, k) = 0.
              arg1 = rdg*delp(i, j, k)*pt(i, j, k)/delz(i, j, k)
              arg2 = kappa*LOG(arg1)
              pkz(i, j, k) = EXP(arg2)
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
      CALL PRT_MXM('pk_b', pk, is, ie, js, je, 0, npz + 1, 1., &
&            gridstruct%area_64, domain)
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
        CALL RAYLEIGH_SUPER(abs0, npx, npy, npz, ks, pfull, phis, &
&                     flagstruct%tau, u, v, w, pt, ua, va, delz, &
&                     gridstruct%agrid, cp_air, rdgas, ptop, hydrostatic&
&                     , .NOT.neststruct%nested, flagstruct%rf_cutoff, rf&
&                     , gridstruct, domain, bd)
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
      CALL DYN_CORE(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir, &
&             cp_air, akap, cappa, grav, hydrostatic, u, v, w, delz, pt&
&             , q, delp, pe, pk, phis, ws, omga, ptop, pfull, ua, va, uc&
&             , vc, mfx, mfy, cx, cy, pkz, peln, q_con, ak, bk, dpx, ks&
&             , gridstruct, flagstruct, flagstructp, neststruct, idiag, &
&             bd, domain, n_map .EQ. 1, i_pack, last_step, gz, pkc, ptc&
&             , crx, xfx, cry, yfx, divgd, delpc, ut, vt, zh, pk3, du, &
&             dv, time_total)
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
        CALL LAGRANGIAN_TO_EULERIAN(last_step, consv_te, ps, pe, delp, &
&                             pkz, pk, mdt, bdt, npz, is, ie, js, je, &
&                             isd, ied, jsd, jed, nq, nwat, sphum, q_con&
&                             , u, v, w, delz, pt, q, phis, zvir, cp_air&
&                             , akap, cappa, kord_mt, kord_wz, &
&                             kord_tracer, kord_tm, peln, te_2d, ng, ua&
&                             , va, omga, dp1, ws, fill, reproduce_sum, &
&                             idiag%id_mdt .GT. 0, dtdt_m, ptop, ak, bk&
&                             , pfull, flagstruct, gridstruct, domain, &
&                             flagstruct%do_sat_adj, hydrostatic, &
&                             hybrid_z, do_omega, flagstruct%adiabatic, &
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
          IF (flagstruct%nf_omega .GT. 0) CALL DEL2_CUBED(omga, 0.18*&
&                                                   gridstruct%da_min, &
&                                                   gridstruct, domain, &
&                                                   npx, npy, npz, &
&                                                   flagstruct%nf_omega&
&                                                   , bd)
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
      IF (idiag%id_amdt .GT. 0) used = SEND_DATA(idiag%id_amdt, te_2d/&
&         bdt, fv_time)
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
!  Differentiation of rayleigh_super in forward (tangent) mode:
!   variations   of useful results: u v w ua va pt
!   with respect to varying inputs: u v w ua va pt
  SUBROUTINE RAYLEIGH_SUPER_TLM(dt, npx, npy, npz, ks, pm, phis, tau, u&
&   , u_tl, v, v_tl, w, w_tl, pt, pt_tl, ua, ua_tl, va, va_tl, delz, &
&   agrid, cp, rg, ptop, hydrostatic, conserve, rf_cutoff, rf, &
&   gridstruct, domain, bd)
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
    REAL, INTENT(INOUT) :: u_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
! cell center vertical wind (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: w_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temp
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ua_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: va_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
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
    REAL*8 :: arg1
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
          arg1 = 0.5*pi*LOG(rf_cutoff/pm(k))/LOG(rf_cutoff/ptop)
          rf(k) = dt/tau0*SIN(arg1)**2
!if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
          kmax = k
        ELSE
          EXIT
        END IF
      END DO
      rf_initialized = .true.
    END IF
    CALL C2L_ORD2_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, gridstruct&
&               , npz, gridstruct%grid_type, bd, gridstruct%nested)
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
                pt_tl(i, j, k) = pt_tl(i, j, k) + 0.5*(1.-u2f(i, j, k)**&
&                 2)*(2*ua(i, j, k)*ua_tl(i, j, k)+2*va(i, j, k)*va_tl(i&
&                 , j, k))/(cp-rg*ptop/pm(k))
                pt(i, j, k) = pt(i, j, k) + 0.5*(ua(i, j, k)**2+va(i, j&
&                 , k)**2)*(1.-u2f(i, j, k)**2)/(cp-rg*ptop/pm(k))
              END DO
            END DO
          ELSE
            DO j=js,je
              DO i=is,ie
                pt_tl(i, j, k) = pt_tl(i, j, k) + 0.5*(1.-u2f(i, j, k)**&
&                 2)*rcv*(2*ua(i, j, k)*ua_tl(i, j, k)+2*va(i, j, k)*&
&                 va_tl(i, j, k)+2*w(i, j, k)*w_tl(i, j, k))
                pt(i, j, k) = pt(i, j, k) + 0.5*(ua(i, j, k)**2+va(i, j&
&                 , k)**2+w(i, j, k)**2)*(1.-u2f(i, j, k)**2)*rcv
              END DO
            END DO
          END IF
        END IF
        DO j=js,je+1
          DO i=is,ie
            u_tl(i, j, k) = 0.5*(u2f(i, j-1, k)+u2f(i, j, k))*u_tl(i, j&
&             , k)
            u(i, j, k) = 0.5*(u2f(i, j-1, k)+u2f(i, j, k))*u(i, j, k)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            v_tl(i, j, k) = 0.5*(u2f(i-1, j, k)+u2f(i, j, k))*v_tl(i, j&
&             , k)
            v(i, j, k) = 0.5*(u2f(i-1, j, k)+u2f(i, j, k))*v(i, j, k)
          END DO
        END DO
        IF (.NOT.hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              w_tl(i, j, k) = u2f(i, j, k)*w_tl(i, j, k)
              w(i, j, k) = u2f(i, j, k)*w(i, j, k)
            END DO
          END DO
        END IF
      END IF
    END DO
  END SUBROUTINE RAYLEIGH_SUPER_TLM
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
    REAL*8 :: arg1
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
          arg1 = 0.5*pi*LOG(rf_cutoff/pm(k))/LOG(rf_cutoff/ptop)
          rf(k) = dt/tau0*SIN(arg1)**2
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
!  Differentiation of rayleigh_friction in forward (tangent) mode:
!   variations   of useful results: u v w ua delz va pt
!   with respect to varying inputs: u v w ua delz va pt
  SUBROUTINE RAYLEIGH_FRICTION_TLM(dt, npx, npy, npz, ks, pm, tau, u, &
&   u_tl, v, v_tl, w, w_tl, pt, pt_tl, ua, ua_tl, va, va_tl, delz, &
&   delz_tl, cp, rg, ptop, hydrostatic, conserve, rf_cutoff, rf, &
&   gridstruct, domain, bd)
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
    REAL, INTENT(INOUT) :: u_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
! cell center vertical wind (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: w_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temp
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ua_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: va_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delz_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: rf(npz)
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! local:
    REAL :: u2f(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: u2f_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
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
    REAL*8 :: arg1
    REAL :: arg10
    REAL :: arg10_tl
    REAL :: result1
    REAL :: result1_tl
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
          arg1 = 0.5*pi*LOG(rf_cutoff/pm(k))/LOG(rf_cutoff/ptop)
          rf(k) = dt/(tau*sday)*SIN(arg1)**2
!if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
          kmax = k
        ELSE
          EXIT
        END IF
      END DO
      rf_initialized = .true.
    END IF
!allocate( u2f(isd:ied,jsd:jed,kmax) )
    CALL C2L_ORD2_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, gridstruct&
&               , npz, gridstruct%grid_type, bd, gridstruct%nested)
    u2f_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,u2f,hydrostatic,ua,va,w)
    DO k=1,kmax
      IF (hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            u2f_tl(i, j, k) = 2*ua(i, j, k)*ua_tl(i, j, k) + 2*va(i, j, &
&             k)*va_tl(i, j, k)
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            u2f_tl(i, j, k) = 2*ua(i, j, k)*ua_tl(i, j, k) + 2*va(i, j, &
&             k)*va_tl(i, j, k) + 2*w(i, j, k)*w_tl(i, j, k)
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2 + w(i, j, k)&
&             **2
          END DO
        END DO
      END IF
    END DO
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_DOMAINS_TLM(u2f, u2f_tl, domain)
    CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,conserve,hydrostatic,pt,u2f,cp,rg, &
!$OMP                                  ptop,pm,rf,delz,rcv,u,v,w)
    DO k=1,kmax
      IF (conserve) THEN
        IF (hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              arg10_tl = u2f_tl(i, j, k)/u000
              arg10 = u2f(i, j, k)/u000
              IF (arg10 .EQ. 0.0) THEN
                result1_tl = 0.0
              ELSE
                result1_tl = arg10_tl/(2.0*SQRT(arg10))
              END IF
              result1 = SQRT(arg10)
              pt_tl(i, j, k) = pt_tl(i, j, k) + 0.5*u2f_tl(i, j, k)*(1.-&
&               1./(1.+rf(k)*result1)**2)/(cp-rg*ptop/pm(k)) + 0.5*u2f(i&
&               , j, k)*2*rf(k)*result1_tl/((cp-rg*ptop/pm(k))*(1.+rf(k)&
&               *result1)**3)
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)/(cp-rg*ptop/&
&               pm(k))*(1.-1./(1.+rf(k)*result1)**2)
            END DO
          END DO
        ELSE
          DO j=js,je
            DO i=is,ie
              delz_tl(i, j, k) = (delz_tl(i, j, k)*pt(i, j, k)-delz(i, j&
&               , k)*pt_tl(i, j, k))/pt(i, j, k)**2
              delz(i, j, k) = delz(i, j, k)/pt(i, j, k)
              arg10_tl = u2f_tl(i, j, k)/u000
              arg10 = u2f(i, j, k)/u000
              IF (arg10 .EQ. 0.0) THEN
                result1_tl = 0.0
              ELSE
                result1_tl = arg10_tl/(2.0*SQRT(arg10))
              END IF
              result1 = SQRT(arg10)
              pt_tl(i, j, k) = pt_tl(i, j, k) + 0.5*rcv*(u2f_tl(i, j, k)&
&               *(1.-1./(1.+rf(k)*result1)**2)+u2f(i, j, k)*2*rf(k)*&
&               result1_tl/(1.+rf(k)*result1)**3)
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)*rcv*(1.-1./(&
&               1.+rf(k)*result1)**2)
              delz_tl(i, j, k) = delz_tl(i, j, k)*pt(i, j, k) + delz(i, &
&               j, k)*pt_tl(i, j, k)
              delz(i, j, k) = delz(i, j, k)*pt(i, j, k)
            END DO
          END DO
        END IF
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+1
          arg10_tl = u2f_tl(i, j, k)/u000
          arg10 = u2f(i, j, k)/u000
          IF (arg10 .EQ. 0.0) THEN
            result1_tl = 0.0
          ELSE
            result1_tl = arg10_tl/(2.0*SQRT(arg10))
          END IF
          result1 = SQRT(arg10)
          u2f_tl(i, j, k) = rf(k)*result1_tl
          u2f(i, j, k) = rf(k)*result1
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          u_tl(i, j, k) = (u_tl(i, j, k)*(1.+0.5*(u2f(i, j-1, k)+u2f(i, &
&           j, k)))-u(i, j, k)*0.5*(u2f_tl(i, j-1, k)+u2f_tl(i, j, k)))/&
&           (1.+0.5*(u2f(i, j-1, k)+u2f(i, j, k)))**2
          u(i, j, k) = u(i, j, k)/(1.+0.5*(u2f(i, j-1, k)+u2f(i, j, k)))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v_tl(i, j, k) = (v_tl(i, j, k)*(1.+0.5*(u2f(i-1, j, k)+u2f(i, &
&           j, k)))-v(i, j, k)*0.5*(u2f_tl(i-1, j, k)+u2f_tl(i, j, k)))/&
&           (1.+0.5*(u2f(i-1, j, k)+u2f(i, j, k)))**2
          v(i, j, k) = v(i, j, k)/(1.+0.5*(u2f(i-1, j, k)+u2f(i, j, k)))
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            w_tl(i, j, k) = (w_tl(i, j, k)*(1.+u2f(i, j, k))-w(i, j, k)*&
&             u2f_tl(i, j, k))/(1.+u2f(i, j, k))**2
            w(i, j, k) = w(i, j, k)/(1.+u2f(i, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE RAYLEIGH_FRICTION_TLM
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
    REAL*8 :: arg1
    REAL :: arg10
    REAL :: result1
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
          arg1 = 0.5*pi*LOG(rf_cutoff/pm(k))/LOG(rf_cutoff/ptop)
          rf(k) = dt/(tau*sday)*SIN(arg1)**2
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
              arg10 = u2f(i, j, k)/u000
              result1 = SQRT(arg10)
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)/(cp-rg*ptop/&
&               pm(k))*(1.-1./(1.+rf(k)*result1)**2)
            END DO
          END DO
        ELSE
          DO j=js,je
            DO i=is,ie
              delz(i, j, k) = delz(i, j, k)/pt(i, j, k)
              arg10 = u2f(i, j, k)/u000
              result1 = SQRT(arg10)
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)*rcv*(1.-1./(&
&               1.+rf(k)*result1)**2)
              delz(i, j, k) = delz(i, j, k)*pt(i, j, k)
            END DO
          END DO
        END IF
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+1
          arg10 = u2f(i, j, k)/u000
          result1 = SQRT(arg10)
          u2f(i, j, k) = rf(k)*result1
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
!  Differentiation of compute_aam in forward (tangent) mode:
!   variations   of useful results: ua aam va m_fac ps
!   with respect to varying inputs: u v delp ua aam va m_fac ps
  SUBROUTINE COMPUTE_AAM_TLM(npz, is, ie, js, je, isd, ied, jsd, jed, &
&   gridstruct, bd, ptop, ua, ua_tl, va, va_tl, u, u_tl, v, v_tl, delp, &
&   delp_tl, aam, aam_tl, ps, ps_tl, m_fac, m_fac_tl)
    IMPLICIT NONE
! Compute vertically (mass) integrated Atmospheric Angular Momentum
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: is, ie, js, je
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    REAL, INTENT(IN) :: ptop
! D grid zonal wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: u_tl(isd:ied, jsd:jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: v_tl(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(isd:ied, jsd:jed, npz)
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua_tl, &
&   va_tl
    REAL, INTENT(OUT) :: aam(is:ie, js:je)
    REAL, INTENT(OUT) :: aam_tl(is:ie, js:je)
    REAL, INTENT(OUT) :: m_fac(is:ie, js:je)
    REAL, INTENT(OUT) :: m_fac_tl(is:ie, js:je)
    REAL, INTENT(OUT) :: ps(isd:ied, jsd:jed)
    REAL, INTENT(OUT) :: ps_tl(isd:ied, jsd:jed)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
! local:
    REAL, DIMENSION(is:ie) :: r1, r2, dm
    REAL, DIMENSION(is:ie) :: dm_tl
    INTEGER :: i, j, k
    INTRINSIC COS
    CALL C2L_ORD2_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, gridstruct&
&               , npz, gridstruct%grid_type, bd, gridstruct%nested)
    dm_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gridstruct,aam,m_fac,ps,ptop,delp,agrav,ua) &
!$OMP                          private(r1, r2, dm)
    DO j=js,je
      DO i=is,ie
        r1(i) = radius*COS(gridstruct%agrid(i, j, 2))
        r2(i) = r1(i)*r1(i)
        aam_tl(i, j) = 0.0
        aam(i, j) = 0.
        m_fac_tl(i, j) = 0.0
        m_fac(i, j) = 0.
        ps_tl(i, j) = 0.0
        ps(i, j) = ptop
      END DO
      DO k=1,npz
        DO i=is,ie
          dm_tl(i) = delp_tl(i, j, k)
          dm(i) = delp(i, j, k)
          ps_tl(i, j) = ps_tl(i, j) + dm_tl(i)
          ps(i, j) = ps(i, j) + dm(i)
          dm_tl(i) = agrav*dm_tl(i)
          dm(i) = dm(i)*agrav
          aam_tl(i, j) = aam_tl(i, j) + r1(i)*ua_tl(i, j, k)*dm(i) + (r2&
&           (i)*omega+r1(i)*ua(i, j, k))*dm_tl(i)
          aam(i, j) = aam(i, j) + (r2(i)*omega+r1(i)*ua(i, j, k))*dm(i)
          m_fac_tl(i, j) = m_fac_tl(i, j) + r2(i)*dm_tl(i)
          m_fac(i, j) = m_fac(i, j) + dm(i)*r2(i)
        END DO
      END DO
    END DO
  END SUBROUTINE COMPUTE_AAM_TLM
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
end module fv_dynamics_tlm_mod
