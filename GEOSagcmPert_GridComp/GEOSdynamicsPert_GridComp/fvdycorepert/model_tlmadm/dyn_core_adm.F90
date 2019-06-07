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
module dyn_core_adm_mod

  use constants_mod,      only: rdgas, radius, cp_air, pi
  use mpp_mod,            only: mpp_pe 
  use mpp_domains_mod,    only: CGRID_NE, DGRID_NE, mpp_get_boundary, mpp_update_domains, &
                                domain2d
  use fv_mp_adm_mod,      only: mpp_get_boundary_adm, mpp_update_domains_adm
  use mpp_parameter_mod,  only: CORNER
  use fv_mp_mod,          only: is_master
  use fv_mp_adm_mod,      only: start_group_halo_update, complete_group_halo_update
  use fv_mp_adm_mod,      only: start_group_halo_update_adm
  use fv_mp_mod,          only: group_halo_update_type
  use sw_core_adm_mod,    only: c_sw, d_sw
  use sw_core_adm_mod,    only: c_sw_fwd, c_sw_bwd, d_sw_fwd, d_sw_bwd
  use a2b_edge_adm_mod,   only: a2b_ord2, a2b_ord4
  use a2b_edge_adm_mod,   only: a2b_ord2_fwd, a2b_ord2_bwd, a2b_ord4_fwd, a2b_ord4_bwd
  use nh_core_adm_mod,    only: Riem_Solver3, Riem_Solver_C, update_dz_c, update_dz_d, nest_halo_nh
  use nh_core_adm_mod,    only: Riem_Solver3_fwd, Riem_Solver_C_fwd, update_dz_c_fwd, update_dz_d_fwd, nest_halo_nh_fwd
  use nh_core_adm_mod,    only: Riem_Solver3_bwd, Riem_Solver_C_bwd, update_dz_c_bwd, update_dz_d_bwd, nest_halo_nh_bwd
  use tp_core_adm_mod,    only: copy_corners
  use tp_core_adm_mod,    only: copy_corners_adm
  use fv_timing_mod,      only: timing_on, timing_off
  use fv_diagnostics_mod, only: prt_maxmin, fv_time, prt_mxm
#ifdef ROT3
  use fv_update_phys_mod, only: update_dwinds_phys
#endif
#if defined (ADA_NUDGE)
  use fv_ada_nudge_mod,   only: breed_slp_inline_ada
#else
  use fv_nwp_nudge_mod,   only: breed_slp_inline, do_adiabatic_init
#endif
  use diag_manager_mod,   only: send_data
  use fv_arrays_mod,      only: fv_grid_type, fv_flags_type, fv_nest_type, fv_diag_type, &
                                fv_grid_bounds_type, R_GRID

  use boundary_adm_mod,   only: extrapolation_BC,  nested_grid_BC_apply_intT
  use boundary_adm_mod,   only: nested_grid_BC_apply_intT_adm

#ifdef SW_DYNAMICS
  use test_cases_mod,     only: test_case, case9_forcing1, case9_forcing2
#endif

  use tapenade_iter,      only: pushcontrol, popcontrol, pushinteger, popinteger, &
                                pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

  use fv_arrays_nlm_mod,  only: fv_flags_pert_type, fpp

implicit none
private

public :: dyn_core, del2_cubed, init_ijk_mem
public :: dyn_core_fwd, del2_cubed_fwd
public :: dyn_core_bwd, del2_cubed_bwd

  real :: ptk, peln1, rgrav
  real :: d3_damp
!  real, allocatable, dimension(:,:,:) ::  ut, vt, crx, cry, xfx, yfx, divgd, &
!                                          zh, du, dv, pkc, delpc, pk3, ptc, gz
! real, parameter:: delt_max = 1.e-1   ! Max dissipative heating/cooling rate
                                       ! 6 deg per 10-min
  real(kind=R_GRID), parameter :: cnst_0p20=0.20d0

!  real, allocatable ::  rf(:)
  logical:: RFF_initialized = .false.
  integer :: kmax=1

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of dyn_core in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
!.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mo
!d.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix
!_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_S
!uper fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 
!fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rema
!p_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_
!mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv
!_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_res
!tart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_
!z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mo
!d.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.
!SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.n
!est_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2
!c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
! sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mo
!d.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_m
!od.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pk3 xfx ws peln q gz du u dv
!                v w delp ua uc ptc mfx delz mfy omga ut divgd
!                pkc delpc va vc yfx pkz pe vt pk zh pt cx cy dpx
!                crx cry
!   with respect to varying inputs: pk3 xfx ws peln q gz du u dv
!                v w delp ua uc ptc delz omga ut divgd pkc delpc
!                va vc yfx pkz pe vt pk zh pt dpx crx cry
!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
  SUBROUTINE DYN_CORE_FWD(npx, npy, npz, ng, sphum, nq, bdt, n_split, &
&   zvir, cp, akap, cappa, grav, hydrostatic, u, v, w, delz, pt, q, delp&
&   , pe, pk, phis, ws, omga, ptop, pfull, ua, va, uc, vc, mfx, mfy, cx&
&   , cy, pkz, peln, q_con, ak, bk, dpx, ks, gridstruct, flagstruct, &
&   flagstructp, neststruct, idiag, bd, domain, init_step, i_pack, &
&   end_step, gz, pkc, ptc, crx, xfx, cry, yfx, divgd, delpc, ut, vt, zh&
&   , pk3, du, dv, time_total)
    IMPLICIT NONE
! end init_step
! Start of the big dynamic time stepping
!allocate(    gz(isd:ied, jsd:jed ,npz+1) )
!  call init_ijk_mem(isd,ied, jsd,jed, npz+1, gz, huge_r)
!allocate(   pkc(isd:ied, jsd:jed ,npz+1) )
!allocate(   ptc(isd:ied, jsd:jed ,npz ) )
!allocate( crx(is :ie+1, jsd:jed,  npz) )
!allocate( xfx(is :ie+1, jsd:jed,  npz) )
!allocate( cry(isd:ied,  js :je+1, npz) )
!allocate( yfx(isd:ied,  js :je+1, npz) )
!allocate( divgd(isd:ied+1,jsd:jed+1,npz) )
!allocate( delpc(isd:ied, jsd:jed  ,npz  ) )
!                    call init_ijk_mem(isd,ied, jsd,jed, npz, delpc, 0.)
!allocate( ut(isd:ied, jsd:jed, npz) )
!                    call init_ijk_mem(isd,ied, jsd,jed, npz, ut, 0.)
!allocate( vt(isd:ied, jsd:jed, npz) )
!                    call init_ijk_mem(isd,ied, jsd,jed, npz, vt, 0.)
!allocate( zh(isd:ied, jsd:jed, npz+1) )
!              call init_ijk_mem(isd,ied, jsd,jed, npz+1, zh, huge_r )
!allocate ( pk3(isd:ied,jsd:jed,npz+1) )
!call init_ijk_mem(isd,ied, jsd,jed, npz+1, pk3, huge_r )
!if (allocated(heat_source)) deallocate( heat_source ) !If ncon == 0 but d_con > 1.e-5, this would not be deallocated in earlier 
!versions of the code
!deallocate(    gz )
!deallocate(   ptc )
!deallocate(   crx )
!deallocate(   xfx )
!deallocate(   cry )
!deallocate(   yfx )
!deallocate( divgd )
!deallocate(   pkc )
!deallocate( delpc )
!if( allocated(ut))   deallocate( ut )
!if( allocated(vt))   deallocate( vt )
!if ( allocated (du) ) deallocate( du )
!if ( allocated (dv) ) deallocate( dv )
!if ( .not. hydrostatic ) then
!     deallocate( zh )
!     if( allocated(pk3) )   deallocate ( pk3 )
!endif
!if( allocated(pem) )   deallocate ( pem )
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: ng, nq, sphum
    INTEGER, INTENT(IN) :: n_split
    REAL, INTENT(IN) :: bdt
    REAL, INTENT(IN) :: zvir, cp, akap, grav
    REAL, INTENT(IN) :: ptop
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: init_step, end_step
    REAL, INTENT(IN) :: pfull(npz)
    REAL, DIMENSION(npz+1), INTENT(IN) :: ak, bk
    INTEGER, INTENT(IN) :: ks
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: i_pack(*)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
! vertical vel. (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! delta-height (m, negative)
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! moist kappa
    REAL, INTENT(INOUT) :: cappa(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! total time (seconds) since start
    REAL, INTENT(IN), OPTIONAL :: time_total
!-----------------------------------------------------------------------
! Auxilliary pressure arrays:
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
! edge pressure (pascal)
    REAL, INTENT(INOUT) :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
! ln(pe)
    REAL, INTENT(INOUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
! pe**kappa
    REAL, INTENT(INOUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL(kind=8), INTENT(INOUT) :: dpx(bd%is:bd%ie, bd%js:bd%je)
!-----------------------------------------------------------------------
! Others:
    REAL, PARAMETER :: near0=1.e-8
    REAL, PARAMETER :: huge_r=1.e8
!-----------------------------------------------------------------------
! w at surface
    REAL :: ws(bd%is:bd%ie, bd%js:bd%je)
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! (uc, vc) are mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua, va
    REAL, INTENT(INOUT) :: q_con(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! The Flux capacitors: accumulated Mass flux arrays
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je, npz), INTENT(INOUT) :: pkz
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    TYPE(FV_FLAGS_PERT_TYPE), INTENT(IN), TARGET :: flagstructp
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(FV_DIAG_TYPE), INTENT(IN) :: idiag
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
!real, allocatable, dimension(:,:,:):: pem, heat_source
    REAL :: pem(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1), heat_source(bd&
&   %isd:bd%ied, bd%jsd:bd%jed, npz)
! Auto 1D & 2D arrays:
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ws3, z_rat
    REAL :: dp_ref(npz)
! surface height (m)
    REAL :: zs(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: p1d(bd%is:bd%ie)
    REAL :: om2d(bd%is:bd%ie, npz)
    REAL :: wbuffer(npy+2, npz)
    REAL :: ebuffer(npy+2, npz)
    REAL :: nbuffer(npx+2, npz)
    REAL :: sbuffer(npx+2, npz)
! ----   For external mode:
    REAL :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: fz(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: heat_s(bd%is:bd%ie, bd%js:bd%je)
    REAL :: damp_vt(npz+1)
    INTEGER :: nord_v(npz+1)
!-------------------------------------
    INTEGER :: hord_m, hord_v, hord_t, hord_p
    INTEGER :: nord_k, nord_w, nord_t
    INTEGER :: ms
!---------------------------------------
    INTEGER :: hord_m_pert, hord_v_pert, hord_t_pert, hord_p_pert
    INTEGER :: nord_k_pert, nord_w_pert, nord_t_pert, nord_v_pert(npz+1)
    REAL :: d2_divg_pert, damp_vt_pert(npz+1), damp_w_pert, damp_t_pert
!---------------------------------------
    INTEGER :: i, j, k, it, iq, n_con, nf_ke
    INTEGER :: iep1, jep1
    REAL :: beta, beta_d, d_con_k, damp_w, damp_t, kgb, cv_air
    REAL :: dt, dt2, rdt
    REAL :: d2_divg
    REAL :: k1k, rdg, dtmp, delt
    LOGICAL :: last_step, remap_step
    LOGICAL :: used
    REAL :: split_timestep_bc
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pkc(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: ptc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cry(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: divgd(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: delpc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ut(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: zh(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk3(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    INTRINSIC LOG
    INTRINSIC REAL
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC EXP
    INTRINSIC ABS
    INTRINSIC SIGN
    INTEGER :: max1
    INTEGER :: max2
    REAL :: min1
    REAL :: min2
    REAL :: abs0
    REAL :: arg1
    REAL :: arg2
    LOGICAL :: arg10
    REAL*8 :: arg11
    LOGICAL :: res
    LOGICAL :: res0
    REAL :: x1
    REAL :: y2
    REAL :: y1

    pem = 0.0
    heat_source = 0.0
    ws3 = 0.0
    z_rat = 0.0
    dp_ref = 0.0
    zs = 0.0
    p1d = 0.0
    om2d = 0.0
    wbuffer = 0.0
    ebuffer = 0.0
    nbuffer = 0.0
    sbuffer = 0.0
    divg2 = 0.0
    wk = 0.0
    fz = 0.0
    heat_s = 0.0
    damp_vt = 0.0
    d2_divg_pert = 0.0
    damp_vt_pert = 0.0
    damp_w_pert = 0.0
    damp_t_pert = 0.0
    beta = 0.0
    beta_d = 0.0
    d_con_k = 0.0
    damp_w = 0.0
    damp_t = 0.0
    kgb = 0.0
    cv_air = 0.0
    dt = 0.0
    dt2 = 0.0
    rdt = 0.0
    d2_divg = 0.0
    k1k = 0.0
    rdg = 0.0
    dtmp = 0.0
    delt = 0.0
    split_timestep_bc = 0.0
    min1 = 0.0
    min2 = 0.0
    abs0 = 0.0
    arg1 = 0.0
    arg2 = 0.0
    arg10 = 0.0
    x1 = 0.0
    y2 = 0.0
    y1 = 0.0
    nord_v = 0
    hord_m = 0
    hord_v = 0
    hord_t = 0
    hord_p = 0
    nord_k = 0
    nord_w = 0
    nord_t = 0
    ms = 0
    hord_m_pert = 0
    hord_v_pert = 0
    hord_t_pert = 0
    hord_p_pert = 0
    nord_k_pert = 0
    nord_w_pert = 0
    nord_t_pert = 0
    nord_v_pert = 0
    i = 0
    j = 0
    k = 0
    it = 0
    iq = 0
    n_con = 0
    nf_ke = 0
    iep1 = 0
    jep1 = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    peln1 = LOG(ptop)
    ptk = ptop**akap
    dt = bdt/REAL(n_split)
    dt2 = 0.5*dt
    rdt = 1.0/dt
    IF (1 .LT. flagstruct%m_split/2) THEN
      ms = flagstruct%m_split/2
    ELSE
      ms = 1
    END IF
    beta = flagstruct%beta
    rdg = -(rdgas/grav)
    cv_air = cp_air - rdgas
! Indexes:
    iep1 = ie + 1
    jep1 = je + 1
    IF (.NOT.hydrostatic) THEN
      rgrav = 1.0/grav
! rg/Cv=0.4
      k1k = akap/(1.-akap)
!$OMP parallel do default(none) shared(npz,dp_ref,ak,bk)
      DO k=1,npz
        dp_ref(k) = ak(k+1) - ak(k) + (bk(k+1)-bk(k))*1.e5
      END DO
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,zs,phis,rgrav)
      DO j=jsd,jed
        DO i=isd,ied
          zs(i, j) = phis(i, j)*rgrav
        END DO
      END DO
    END IF
!allocate( du(isd:ied,  jsd:jed+1,npz) )
!call init_ijk_mem(isd,ied,   jsd,jed+1, npz, du, 0.)
!allocate( dv(isd:ied+1,jsd:jed,  npz) )
!call init_ijk_mem(isd,ied+1, jsd,jed  , npz, dv, 0.)
! Empty the "flux capacitors"
!call init_ijk_mem(is, ie+1, js,  je,   npz, mfx, 0.)
    CALL PUSHREALARRAY(mfx, (bd%ie-bd%is+2)*(bd%je-bd%js+1)*npz)
    mfx = 0.0
!call init_ijk_mem(is, ie  , js,  je+1, npz, mfy, 0.)
    CALL PUSHREALARRAY(mfy, (bd%ie-bd%is+1)*(bd%je-bd%js+2)*npz)
    mfy = 0.0
!call init_ijk_mem(is, ie+1, jsd, jed,  npz, cx, 0.)
    CALL PUSHREALARRAY(cx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
    cx = 0.0
!call init_ijk_mem(isd, ied, js,  je+1, npz, cy, 0.)
    CALL PUSHREALARRAY(cy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    cy = 0.0
    IF (flagstruct%d_con .GT. 1.0e-5) heat_source = 0.0
!allocate( heat_source(isd:ied, jsd:jed, npz) )
!call init_ijk_mem(isd, ied, jsd, jed, npz, heat_source, 0.)
    IF (flagstruct%convert_ke .OR. flagstruct%vtdm4 .GT. 1.e-4) THEN
      CALL PUSHCONTROL(2,2)
      n_con = npz
    ELSE IF (flagstruct%d2_bg_k1 .LT. 1.e-3) THEN
      CALL PUSHCONTROL(2,1)
      n_con = 0
    ELSE IF (flagstruct%d2_bg_k2 .LT. 1.e-3) THEN
      CALL PUSHCONTROL(2,0)
      n_con = 1
    ELSE
      CALL PUSHCONTROL(2,0)
      n_con = 2
    END IF
!-----------------------------------------------------
    DO it=1,n_split
!-----------------------------------------------------
      IF (flagstruct%breed_vortex_inline .OR. it .EQ. n_split) THEN
        remap_step = .true.
      ELSE
        remap_step = .false.
      END IF
      IF (flagstruct%fv_debug) THEN
        res = IS_MASTER()
        IF (res) THEN
          CALL PUSHCONTROL(1,0)
          WRITE(*, *) 'n_split loop, it=', it
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%nested) split_timestep_bc = REAL(n_split*flagstruct&
&         %k_split + neststruct%nest_timestep)
!First split timestep has split_timestep_BC = n_split*k_split
!   to do time-extrapolation on BCs.
      IF (nq .GT. 0) THEN
        IF (flagstruct%inline_q) THEN
          CALL PUSHREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz&
&                       *nq)
          CALL START_GROUP_HALO_UPDATE(i_pack(10), q, domain)
          CALL PUSHCONTROL(2,0)
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,2)
      END IF
      IF (.NOT.hydrostatic) THEN
        CALL PUSHREALARRAY(w, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
        CALL START_GROUP_HALO_UPDATE(i_pack(7), w, domain)
        IF (it .EQ. 1) THEN
          IF (gridstruct%nested) THEN
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,gz,zs,delz)
            DO j=jsd,jed
              DO i=isd,ied
                CALL PUSHREALARRAY(gz(i, j, npz+1))
                gz(i, j, npz+1) = zs(i, j)
              END DO
              DO k=npz,1,-1
                DO i=isd,ied
                  CALL PUSHREALARRAY(gz(i, j, k))
                  gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)
                END DO
              END DO
            END DO
            CALL PUSHCONTROL(1,0)
          ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gz,zs,delz)
            DO j=js,je
              DO i=is,ie
                CALL PUSHREALARRAY(gz(i, j, npz+1))
                gz(i, j, npz+1) = zs(i, j)
              END DO
              DO k=npz,1,-1
                DO i=is,ie
                  CALL PUSHREALARRAY(gz(i, j, k))
                  gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)
                END DO
              END DO
            END DO
            CALL PUSHCONTROL(1,1)
          END IF
          CALL PUSHREALARRAY(gz, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(&
&                       npz+1))
          CALL START_GROUP_HALO_UPDATE(i_pack(5), gz, domain)
          CALL PUSHCONTROL(2,0)
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,2)
      END IF
      IF (it .EQ. 1) THEN
        CALL PUSHREALARRAY(beta_d)
        beta_d = 0.
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHREALARRAY(beta_d)
        beta_d = beta
        CALL PUSHCONTROL(1,1)
      END IF
      IF (it .EQ. n_split .AND. end_step) THEN
        IF (flagstruct%use_old_omega) THEN
          pem = 0.0
!allocate ( pem(is-1:ie+1,npz+1,js-1:je+1) )
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pem,delp,ptop)
          DO j=js-1,je+1
            DO i=is-1,ie+1
              pem(i, 1, j) = ptop
            END DO
            DO k=1,npz
              DO i=is-1,ie+1
                pem(i, k+1, j) = pem(i, k, j) + delp(i, j, k)
              END DO
            END DO
          END DO
          CALL PUSHCONTROL(2,0)
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
        last_step = .true.
      ELSE
        CALL PUSHCONTROL(2,2)
        last_step = .false.
      END IF
!$OMP parallel do default(none) shared(npz,isd,jsd,delpc,delp,ptc,pt,u,v,w,uc,vc,ua,va, &
!$OMP                                  omga,ut,vt,divgd,flagstruct,dt2,hydrostatic,bd,  &
!$OMP                                  gridstruct)
      DO k=1,npz
        CALL C_SW_FWD(delpc(isd:ied, jsd:jed, k), delp(isd:ied, jsd:jed&
&               , k), ptc(isd:ied, jsd:jed, k), pt(isd:ied, jsd:jed, k)&
&               , u(isd:ied, jsd:jed+1, k), v(isd:ied+1, jsd:jed, k), w(&
&               isd:ied, jsd:jed, k), uc(isd:ied+1, jsd:jed, k), vc(isd:&
&               ied, jsd:jed+1, k), ua(isd:ied, jsd:jed, k), va(isd:ied&
&               , jsd:jed, k), omga(isd:ied, jsd:jed, k), ut(isd:ied, &
&               jsd:jed, k), vt(isd:ied, jsd:jed, k), divgd(isd:ied+1, &
&               jsd:jed+1, k), flagstruct%nord, dt2, hydrostatic, .true.&
&               , bd, gridstruct, flagstruct)
      END DO
      IF (flagstruct%nord .GT. 0) THEN
        CALL START_GROUP_HALO_UPDATE(i_pack(3), divgd, domain, position=&
&                              corner)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%nested) THEN
        arg1 = split_timestep_bc + 0.5
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL PUSHREALARRAY(delpc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*&
&                     npz)
        CALL NESTED_GRID_BC_APPLY_INTT(delpc, 0, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 0.5, REAL(n_split*&
&                                flagstruct%k_split), neststruct%delp_bc&
&                                , bctype=neststruct%nestbctype)
        arg1 = split_timestep_bc + 0.5
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL PUSHREALARRAY(ptc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz&
&                    )
        CALL NESTED_GRID_BC_APPLY_INTT(ptc, 0, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 0.5, REAL(n_split*&
&                                flagstruct%k_split), neststruct%pt_bc, &
&                                bctype=neststruct%nestbctype)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
! end hydro check
      IF (hydrostatic) THEN
        CALL GEOPK_FWD(ptop, pe, peln, delpc, pkc, gz, phis, ptc, q_con&
&                , pkz, npz, akap, .true., gridstruct%nested, .false., &
&                npx, npy, flagstruct%a2b_ord, bd)
        CALL PUSHCONTROL(2,0)
      ELSE
        IF (it .EQ. 1) THEN
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,zh,gz)
          DO k=1,npz+1
            DO j=jsd,jed
              DO i=isd,ied
! Save edge heights for update_dz_d
                zh(i, j, k) = gz(i, j, k)
              END DO
            END DO
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,zh,gz)
          DO k=1,npz+1
            DO j=jsd,jed
              DO i=isd,ied
                CALL PUSHREALARRAY(gz(i, j, k))
                gz(i, j, k) = zh(i, j, k)
              END DO
            END DO
          END DO
          CALL PUSHCONTROL(1,1)
        END IF
        CALL UPDATE_DZ_C_FWD(is, ie, js, je, npz, ng, dt2, dp_ref, zs, &
&                      gridstruct%area, ut, vt, gz, ws3, npx, npy, &
&                      gridstruct%sw_corner, gridstruct%se_corner, &
&                      gridstruct%ne_corner, gridstruct%nw_corner, bd, &
&                      gridstruct%grid_type)
        CALL RIEM_SOLVER_C_FWD(ms, dt2, is, ie, js, je, npz, ng, akap, &
&                        cappa, cp, ptop, phis, omga, ptc, q_con, delpc&
&                        , gz, pkc, ws3, flagstruct%p_fac, flagstruct%&
&                        a_imp, flagstruct%scale_z)
        IF (gridstruct%nested) THEN
          arg1 = split_timestep_bc + 0.5
          arg2 = REAL(n_split*flagstruct%k_split)
          CALL PUSHREALARRAY(delz, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*&
&                       npz)
          CALL NESTED_GRID_BC_APPLY_INTT(delz, 0, 0, npx, npy, npz, bd, &
&                                  split_timestep_bc + 0.5, REAL(n_split&
&                                  *flagstruct%k_split), neststruct%&
&                                  delz_bc, bctype=neststruct%nestbctype&
&                                 )
!Compute gz/pkc
!NOTE: nominally only need to compute quantities one out in the halo for p_grad_c
!(instead of entire halo)
          CALL NEST_HALO_NH_FWD(ptop, grav, akap, cp, delpc, delz, ptc, &
&                         phis, pkc, gz, pk3, npx, npy, npz, gridstruct%&
&                         nested, .false., .false., .false., bd)
          CALL PUSHCONTROL(2,1)
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
      END IF
      CALL P_GRAD_C_FWD(dt2, npz, delpc, pkc, gz, uc, vc, bd, gridstruct&
&                 %rdxc, gridstruct%rdyc, hydrostatic)
      CALL START_GROUP_HALO_UPDATE(i_pack(9), uc, vc, domain, gridtype=&
&                            cgrid_ne)
      IF (gridstruct%nested) THEN
!On a nested grid we have to do SOMETHING with uc and vc in
! the boundary halo, particularly at the corners of the
! domain and of each processor element. We must either
! apply an interpolated BC, or extrapolate into the
! boundary halo
! NOTE:
!The update_domains calls for uc and vc need to go BEFORE the BCs to ensure cross-restart
!bitwise-consistent solutions when doing the spatial extrapolation; should not make a
!difference for interpolated BCs from the coarse grid.
        arg1 = split_timestep_bc + 0.5
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(vc, 0, 1, npx, npy, npz, bd, &
&                                split_timestep_bc + 0.5, REAL(n_split*&
&                                flagstruct%k_split), neststruct%vc_bc, &
&                                bctype=neststruct%nestbctype)
        arg1 = split_timestep_bc + 0.5
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(uc, 1, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 0.5, REAL(n_split*&
&                                flagstruct%k_split), neststruct%uc_bc, &
&                                bctype=neststruct%nestbctype)
!QUESTION: What to do with divgd in nested halo?
        arg1 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(divgd, 1, 1, npx, npy, npz, bd, &
&                                split_timestep_bc, REAL(n_split*&
&                                flagstruct%k_split), neststruct%divg_bc&
&                                , bctype=neststruct%nestbctype)
!!$            if (is == 1 .and. js == 1) then
!!$               do j=jsd,5
!!$                  write(mpp_pe()+2000,*) j, divg(isd:5,j,1)
!!$            endif
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%nested .AND. flagstruct%inline_q) THEN
        DO iq=1,nq
          arg1 = split_timestep_bc + 1
          arg2 = REAL(n_split*flagstruct%k_split)
          CALL PUSHREALARRAY(q(isd:ied, jsd:jed, :, iq), (ied-isd+1)*(&
&                       jed-jsd+1)*npz)
          CALL NESTED_GRID_BC_APPLY_INTT(q(isd:ied, jsd:jed, :, iq), 0, &
&                                  0, npx, npy, npz, bd, &
&                                  split_timestep_bc + 1, REAL(n_split*&
&                                  flagstruct%k_split), neststruct%q_bc(&
&                                  iq), bctype=neststruct%nestbctype)
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
!$OMP parallel do default(none) shared(npz,flagstruct,nord_v,pfull,damp_vt,hydrostatic,last_step, &
!$OMP                                  is,ie,js,je,isd,ied,jsd,jed,omga,delp,gridstruct,npx,npy,  &
!$OMP                                  ng,zh,vt,ptc,pt,u,v,w,uc,vc,ua,va,divgd,mfx,mfy,cx,cy,     &
!$OMP                                  crx,cry,xfx,yfx,q_con,zvir,sphum,nq,q,dt,bd,rdt,iep1,jep1, &
!$OMP                                  heat_source)                                               &
!$OMP                          private(nord_k, nord_w, nord_t, damp_w, damp_t, d2_divg,   &
!$OMP                          d_con_k,kgb, hord_m, hord_v, hord_t, hord_p, wk, heat_s, z_rat)
      DO k=1,npz
        hord_m = flagstruct%hord_mt
        hord_t = flagstruct%hord_tm
        hord_v = flagstruct%hord_vt
        hord_p = flagstruct%hord_dp
        CALL PUSHINTEGER(nord_k)
        nord_k = flagstruct%nord
!      if ( k==npz ) then
        kgb = flagstruct%ke_bg
        IF (2 .GT. flagstruct%nord) THEN
          CALL PUSHINTEGER(nord_v(k))
          nord_v(k) = flagstruct%nord
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHINTEGER(nord_v(k))
          nord_v(k) = 2
          CALL PUSHCONTROL(1,1)
        END IF
        IF (0.20 .GT. flagstruct%d2_bg) THEN
          CALL PUSHREALARRAY(d2_divg)
          d2_divg = flagstruct%d2_bg
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(d2_divg)
          d2_divg = 0.20
          CALL PUSHCONTROL(1,1)
        END IF
        IF (flagstruct%do_vort_damp) THEN
! for delp, delz, and vorticity
          CALL PUSHREALARRAY(damp_vt(k))
          damp_vt(k) = flagstruct%vtdm4
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(damp_vt(k))
          damp_vt(k) = 0.
          CALL PUSHCONTROL(1,1)
        END IF
        CALL PUSHINTEGER(nord_w)
        nord_w = nord_v(k)
        CALL PUSHINTEGER(nord_t)
        nord_t = nord_v(k)
        damp_w = damp_vt(k)
        damp_t = damp_vt(k)
        d_con_k = flagstruct%d_con
        IF (npz .EQ. 1 .OR. flagstruct%n_sponge .LT. 0) THEN
          CALL PUSHCONTROL(3,6)
          d2_divg = flagstruct%d2_bg
        ELSE IF (k .EQ. 1) THEN
! Sponge layers with del-2 damping on divergence, vorticity, w, z, and air mass (delp).
! no special damping of potential temperature in sponge layers
! Divergence damping:
          nord_k = 0
          IF (0.01 .LT. flagstruct%d2_bg) THEN
            IF (flagstruct%d2_bg .LT. flagstruct%d2_bg_k1) THEN
              CALL PUSHCONTROL(1,0)
              d2_divg = flagstruct%d2_bg_k1
            ELSE
              CALL PUSHCONTROL(1,0)
              d2_divg = flagstruct%d2_bg
            END IF
          ELSE IF (0.01 .LT. flagstruct%d2_bg_k1) THEN
            CALL PUSHCONTROL(1,1)
            d2_divg = flagstruct%d2_bg_k1
          ELSE
            CALL PUSHCONTROL(1,1)
            d2_divg = 0.01
          END IF
! Vertical velocity:
          nord_w = 0
          damp_w = d2_divg
          IF (flagstruct%do_vort_damp) THEN
! damping on delp and vorticity:
            CALL PUSHINTEGER(nord_v(k))
            nord_v(k) = 0
            CALL PUSHREALARRAY(damp_vt(k))
            damp_vt(k) = 0.5*d2_divg
            CALL PUSHCONTROL(3,4)
          ELSE
            CALL PUSHCONTROL(3,5)
          END IF
          d_con_k = 0.
        ELSE
          IF (2 .LT. flagstruct%n_sponge - 1) THEN
            max1 = flagstruct%n_sponge - 1
          ELSE
            max1 = 2
          END IF
          IF (k .EQ. max1 .AND. flagstruct%d2_bg_k2 .GT. 0.01) THEN
            nord_k = 0
            IF (flagstruct%d2_bg .LT. flagstruct%d2_bg_k2) THEN
              d2_divg = flagstruct%d2_bg_k2
            ELSE
              d2_divg = flagstruct%d2_bg
            END IF
            nord_w = 0
            damp_w = d2_divg
            IF (flagstruct%do_vort_damp) THEN
              CALL PUSHINTEGER(nord_v(k))
              nord_v(k) = 0
              CALL PUSHREALARRAY(damp_vt(k))
              damp_vt(k) = 0.5*d2_divg
              CALL PUSHCONTROL(3,2)
            ELSE
              CALL PUSHCONTROL(3,3)
            END IF
            d_con_k = 0.
          ELSE
            IF (3 .LT. flagstruct%n_sponge) THEN
              max2 = flagstruct%n_sponge
            ELSE
              max2 = 3
            END IF
            IF (k .EQ. max2 .AND. flagstruct%d2_bg_k2 .GT. 0.05) THEN
              nord_k = 0
              IF (flagstruct%d2_bg .LT. 0.2*flagstruct%d2_bg_k2) THEN
                CALL PUSHCONTROL(3,1)
                d2_divg = 0.2*flagstruct%d2_bg_k2
              ELSE
                CALL PUSHCONTROL(3,1)
                d2_divg = flagstruct%d2_bg
              END IF
              nord_w = 0
              damp_w = d2_divg
              d_con_k = 0.
            ELSE
              CALL PUSHCONTROL(3,0)
            END IF
          END IF
        END IF
        CALL PUSHINTEGER(hord_m_pert)
        hord_m_pert = flagstructp%hord_mt_pert
        CALL PUSHINTEGER(hord_t_pert)
        hord_t_pert = flagstructp%hord_tm_pert
        CALL PUSHINTEGER(hord_v_pert)
        hord_v_pert = flagstructp%hord_vt_pert
        CALL PUSHINTEGER(hord_p_pert)
        hord_p_pert = flagstructp%hord_dp_pert
        CALL PUSHINTEGER(nord_k_pert)
        nord_k_pert = flagstructp%nord_pert
        IF (2 .GT. flagstructp%nord_pert) THEN
          CALL PUSHINTEGER(nord_v_pert(k))
          nord_v_pert(k) = flagstructp%nord_pert
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHINTEGER(nord_v_pert(k))
          nord_v_pert(k) = 2
          CALL PUSHCONTROL(1,1)
        END IF
        IF (0.20 .GT. flagstructp%d2_bg_pert) THEN
          CALL PUSHREALARRAY(d2_divg_pert)
          d2_divg_pert = flagstructp%d2_bg_pert
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(d2_divg_pert)
          d2_divg_pert = 0.20
          CALL PUSHCONTROL(1,1)
        END IF
        IF (flagstructp%do_vort_damp_pert) THEN
! for delp, delz, and vorticity
          CALL PUSHREALARRAY(damp_vt_pert(k))
          damp_vt_pert(k) = flagstructp%vtdm4_pert
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(damp_vt_pert(k))
          damp_vt_pert(k) = 0.
          CALL PUSHCONTROL(1,1)
        END IF
        CALL PUSHINTEGER(nord_t_pert)
        nord_t_pert = nord_v_pert(k)
        CALL PUSHREALARRAY(damp_t_pert)
        damp_t_pert = damp_vt_pert(k)
!Sponge layers for the pertuabtiosn
        IF (k .LE. flagstructp%n_sponge_pert) THEN
          IF (k .LE. flagstructp%n_sponge_pert - 1) THEN
            IF (flagstructp%hord_ks_traj) THEN
              hord_m = flagstructp%hord_mt_ks_traj
              hord_t = flagstructp%hord_tm_ks_traj
              hord_v = flagstructp%hord_vt_ks_traj
              hord_p = flagstructp%hord_dp_ks_traj
            END IF
            IF (flagstructp%hord_ks_pert) THEN
              CALL PUSHCONTROL(1,0)
              hord_m_pert = flagstructp%hord_mt_ks_pert
              hord_t_pert = flagstructp%hord_tm_ks_pert
              hord_v_pert = flagstructp%hord_vt_ks_pert
              hord_p_pert = flagstructp%hord_dp_ks_pert
            ELSE
              CALL PUSHCONTROL(1,0)
            END IF
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
          nord_k_pert = 0
          IF (k .EQ. 1) THEN
            IF (0.01 .LT. flagstructp%d2_bg_pert) THEN
              IF (flagstructp%d2_bg_pert .LT. flagstructp%d2_bg_k1_pert&
&             ) THEN
                CALL PUSHCONTROL(3,1)
                d2_divg_pert = flagstructp%d2_bg_k1_pert
              ELSE
                CALL PUSHCONTROL(3,1)
                d2_divg_pert = flagstructp%d2_bg_pert
              END IF
            ELSE IF (0.01 .LT. flagstructp%d2_bg_k1_pert) THEN
              CALL PUSHCONTROL(3,0)
              d2_divg_pert = flagstructp%d2_bg_k1_pert
            ELSE
              CALL PUSHCONTROL(3,0)
              d2_divg_pert = 0.01
            END IF
          ELSE IF (k .EQ. 2) THEN
            IF (0.01 .LT. flagstructp%d2_bg_pert) THEN
              IF (flagstructp%d2_bg_pert .LT. flagstructp%d2_bg_k2_pert&
&             ) THEN
                CALL PUSHCONTROL(3,3)
                d2_divg_pert = flagstructp%d2_bg_k2_pert
              ELSE
                CALL PUSHCONTROL(3,3)
                d2_divg_pert = flagstructp%d2_bg_pert
              END IF
            ELSE IF (0.01 .LT. flagstructp%d2_bg_k2_pert) THEN
              CALL PUSHCONTROL(3,2)
              d2_divg_pert = flagstructp%d2_bg_k2_pert
            ELSE
              CALL PUSHCONTROL(3,2)
              d2_divg_pert = 0.01
            END IF
          ELSE IF (0.01 .LT. flagstructp%d2_bg_pert) THEN
            IF (flagstructp%d2_bg_pert .LT. flagstructp%d2_bg_ks_pert) &
&           THEN
              CALL PUSHCONTROL(3,5)
              d2_divg_pert = flagstructp%d2_bg_ks_pert
            ELSE
              CALL PUSHCONTROL(3,5)
              d2_divg_pert = flagstructp%d2_bg_pert
            END IF
          ELSE IF (0.01 .LT. flagstructp%d2_bg_ks_pert) THEN
            CALL PUSHCONTROL(3,4)
            d2_divg_pert = flagstructp%d2_bg_ks_pert
          ELSE
            CALL PUSHCONTROL(3,4)
            d2_divg_pert = 0.01
          END IF
          IF (flagstructp%do_vort_damp_pert) THEN
            CALL PUSHINTEGER(nord_v_pert(k))
            nord_v_pert(k) = 0
            CALL PUSHREALARRAY(damp_vt_pert(k))
            damp_vt_pert(k) = 0.5*d2_divg_pert
            CALL PUSHCONTROL(2,0)
          ELSE
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
!Tapenade issue if not defined at level npz+1
        CALL PUSHREALARRAY(damp_vt(npz+1))
        damp_vt(npz+1) = damp_vt(npz)
        CALL PUSHREALARRAY(damp_vt_pert(npz+1))
        damp_vt_pert(npz+1) = damp_vt_pert(npz)
        CALL PUSHINTEGER(nord_v(npz+1))
        nord_v(npz+1) = nord_v(npz)
        CALL PUSHINTEGER(nord_v_pert(npz+1))
        nord_v_pert(npz+1) = nord_v_pert(npz)
        IF (hydrostatic .AND. (.NOT.flagstruct%use_old_omega) .AND. &
&           last_step) THEN
! Average horizontal "convergence" to cell center
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(omga(i, j, k))
              omga(i, j, k) = delp(i, j, k)
            END DO
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
!--- external mode divergence damping ---
        IF (flagstruct%d_ext .GT. 0.) THEN
          CALL A2B_ORD2_FWD(delp(isd:ied, jsd:jed, k), wk, gridstruct, &
&                     npx, npy, is, ie, js, je, ng, .false.)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (.NOT.hydrostatic .AND. flagstruct%do_f3d) THEN
! Correction factor for 3D Coriolis force
          DO j=jsd,jed
            DO i=isd,ied
              z_rat(i, j) = 1. + (zh(i, j, k)+zh(i, j, k+1))/radius
            END DO
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        CALL D_SW_FWD(vt(isd:ied, jsd:jed, k), delp(isd:ied, jsd:jed, k)&
&               , ptc(isd:ied, jsd:jed, k), pt(isd:ied, jsd:jed, k), u(&
&               isd:ied, jsd:jed+1, k), v(isd:ied+1, jsd:jed, k), w(isd:&
&               ied, jsd:jed, k), uc(isd:ied+1, jsd:jed, k), vc(isd:ied&
&               , jsd:jed+1, k), ua(isd:ied, jsd:jed, k), va(isd:ied, &
&               jsd:jed, k), divgd(isd:ied+1, jsd:jed+1, k), mfx(is:ie+1&
&               , js:je, k), mfy(is:ie, js:je+1, k), cx(is:ie+1, jsd:jed&
&               , k), cy(isd:ied, js:je+1, k), crx(is:ie+1, jsd:jed, k)&
&               , cry(isd:ied, js:je+1, k), xfx(is:ie+1, jsd:jed, k), &
&               yfx(isd:ied, js:je+1, k), q_con(isd:ied, jsd:jed, 1), &
&               z_rat(isd:ied, jsd:jed), kgb, heat_s, dpx, zvir, sphum, &
&               nq, q, k, npz, flagstruct%inline_q, dt, flagstruct%&
&               hord_tr, hord_m, hord_v, hord_t, hord_p, nord_k, nord_v(&
&               k), nord_w, nord_t, flagstruct%dddmp, d2_divg, &
&               flagstruct%d4_bg, damp_vt(k), damp_w, damp_t, d_con_k, &
&               hydrostatic, gridstruct, flagstruct, bd, flagstructp%&
&               hord_tr_pert, hord_m_pert, hord_v_pert, hord_t_pert, &
&               hord_p_pert, flagstructp%split_damp, nord_k_pert, &
&               nord_v_pert(k), nord_w_pert, nord_t_pert, flagstructp%&
&               dddmp_pert, d2_divg_pert, flagstructp%d4_bg_pert, &
&               damp_vt_pert(k), damp_w_pert, damp_t_pert)
        IF (hydrostatic .AND. (.NOT.flagstruct%use_old_omega) .AND. &
&           last_step) THEN
! Average horizontal "convergence" to cell center
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(omga(i, j, k))
              omga(i, j, k) = omga(i, j, k)*(xfx(i, j, k)-xfx(i+1, j, k)&
&               +yfx(i, j, k)-yfx(i, j+1, k))*gridstruct%rarea(i, j)*rdt
            END DO
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (flagstruct%d_ext .GT. 0.) THEN
          DO j=js,jep1
            DO i=is,iep1
! delp at cell corners
              CALL PUSHREALARRAY(ptc(i, j, k))
              ptc(i, j, k) = wk(i, j)
            END DO
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (flagstruct%d_con .GT. 1.0e-5) THEN
! Average horizontal "convergence" to cell center
          DO j=js,je
            DO i=is,ie
              heat_source(i, j, k) = heat_source(i, j, k) + heat_s(i, j)
            END DO
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
! end openMP k-loop
      IF (flagstruct%fill_dp) THEN
        CALL MIX_DP_FWD(hydrostatic, w, delp, pt, npz, ak, bk, .false., &
&                 flagstruct%fv_debug, bd)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      CALL PUSHREALARRAY(delp, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE(i_pack(1), delp, domain, complete=&
&                            .true.)
      CALL PUSHREALARRAY(pt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE(i_pack(1), pt, domain, complete=&
&                            .true.)
      IF (flagstruct%d_ext .GT. 0.) THEN
        CALL PUSHREALARRAY(d2_divg)
        d2_divg = flagstruct%d_ext*gridstruct%da_min_c
!$OMP parallel do default(none) shared(is,iep1,js,jep1,npz,wk,ptc,divg2,vt,d2_divg)
        DO j=js,jep1
          DO i=is,iep1
            CALL PUSHREALARRAY(wk(i, j))
            wk(i, j) = ptc(i, j, 1)
            divg2(i, j) = wk(i, j)*vt(i, j, 1)
          END DO
          DO k=2,npz
            DO i=is,iep1
              CALL PUSHREALARRAY(wk(i, j))
              wk(i, j) = wk(i, j) + ptc(i, j, k)
              divg2(i, j) = divg2(i, j) + ptc(i, j, k)*vt(i, j, k)
            END DO
          END DO
          DO i=is,iep1
            CALL PUSHREALARRAY(divg2(i, j))
            divg2(i, j) = d2_divg*divg2(i, j)/wk(i, j)
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        divg2(:, :) = 0.
        CALL PUSHCONTROL(1,1)
      END IF
!Want to move this block into the hydro/nonhydro branch above and merge the two if structures
      IF (gridstruct%nested) THEN
        arg1 = split_timestep_bc + 1
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(delp, 0, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 1, REAL(n_split*&
&                                flagstruct%k_split), neststruct%delp_bc&
&                                , bctype=neststruct%nestbctype)
        arg1 = split_timestep_bc + 1
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(pt, 0, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 1, REAL(n_split*&
&                                flagstruct%k_split), neststruct%pt_bc, &
&                                bctype=neststruct%nestbctype)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
! end hydro check
      IF (hydrostatic) THEN
        CALL GEOPK_FWD(ptop, pe, peln, delp, pkc, gz, phis, pt, q_con, &
&                pkz, npz, akap, .false., gridstruct%nested, .true., npx&
&                , npy, flagstruct%a2b_ord, bd)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL UPDATE_DZ_D_FWD(nord_v, damp_vt, flagstruct%hord_tm, is, ie&
&                      , js, je, npz, ng, npx, npy, gridstruct%area, &
&                      gridstruct%rarea, dp_ref, zs, zh, crx, cry, xfx, &
&                      yfx, delz, ws, rdt, gridstruct, bd, flagstructp%&
&                      hord_tm_pert)
        arg10 = beta .LT. -0.1
        CALL RIEM_SOLVER3_FWD(flagstruct%m_split, dt, is, ie, js, je, &
&                       npz, ng, isd, ied, jsd, jed, akap, cappa, cp, &
&                       ptop, zs, q_con, w, delz, pt, delp, zh, pe, pkc&
&                       , pk3, pk, peln, ws, flagstruct%scale_z, &
&                       flagstruct%p_fac, flagstruct%a_imp, flagstruct%&
&                       use_logp, remap_step, arg10)
        IF (gridstruct%square_domain) THEN
          CALL PUSHREALARRAY(zh, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(&
&                       npz+1))
          CALL START_GROUP_HALO_UPDATE(i_pack(4), zh, domain)
          CALL PUSHREALARRAY(pkc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(&
&                       npz+1))
          CALL START_GROUP_HALO_UPDATE(i_pack(5), pkc, domain, whalo=2, &
&                                ehalo=2, shalo=2, nhalo=2)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(zh, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(&
&                       npz+1))
          CALL START_GROUP_HALO_UPDATE(i_pack(4), zh, domain, complete=&
&                                .true.)
          CALL PUSHREALARRAY(pkc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(&
&                       npz+1))
          CALL START_GROUP_HALO_UPDATE(i_pack(4), pkc, domain, complete=&
&                                .true.)
          CALL PUSHCONTROL(1,1)
        END IF
        IF (remap_step) THEN
          CALL PE_HALO_FWD(is, ie, js, je, isd, ied, jsd, jed, npz, ptop&
&                    , pe, delp)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (flagstruct%use_logp) THEN
          CALL PLN_HALO_FWD(is, ie, js, je, isd, ied, jsd, jed, npz, &
&                     ptop, pk3, delp)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PK3_HALO_FWD(is, ie, js, je, isd, ied, jsd, jed, npz, &
&                     ptop, akap, pk3, delp)
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%nested) THEN
          arg1 = split_timestep_bc + 1.
          arg2 = REAL(n_split*flagstruct%k_split)
          CALL PUSHREALARRAY(delz, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*&
&                       npz)
          CALL NESTED_GRID_BC_APPLY_INTT(delz, 0, 0, npx, npy, npz, bd, &
&                                  split_timestep_bc + 1., REAL(n_split*&
&                                  flagstruct%k_split), neststruct%&
&                                  delz_bc, bctype=neststruct%nestbctype&
&                                 )
!Compute gz/pkc/pk3; note that now pkc should be nonhydro pert'n pressure
          CALL NEST_HALO_NH_FWD(ptop, grav, akap, cp, delp, delz, pt, &
&                         phis, pkc, gz, pk3, npx, npy, npz, gridstruct%&
&                         nested, .true., .true., .true., bd)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gz,zh,grav)
        DO k=1,npz+1
          DO j=js-2,je+2
            DO i=is-2,ie+2
              CALL PUSHREALARRAY(gz(i, j, k))
              gz(i, j, k) = zh(i, j, k)*grav
            END DO
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      END IF
      IF (remap_step .AND. hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pk,pkc)
        DO k=1,npz+1
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(pk(i, j, k))
              pk(i, j, k) = pkc(i, j, k)
            END DO
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
!----------------------------
! Compute pressure gradient:
!----------------------------
      IF (hydrostatic) THEN
        IF (beta .GT. 0.) THEN
          CALL GRAD1_P_UPDATE_FWD(divg2, u, v, pkc, gz, du, dv, dt, ng, &
&                           gridstruct, bd, npx, npy, npz, ptop, beta_d&
&                           , flagstruct%a2b_ord)
          CALL PUSHCONTROL(3,0)
        ELSE
          CALL ONE_GRAD_P_FWD(u, v, pkc, gz, divg2, delp, dt, ng, &
&                       gridstruct, bd, npx, npy, npz, ptop, hydrostatic&
&                       , flagstruct%a2b_ord, flagstruct%d_ext)
          CALL PUSHCONTROL(3,1)
        END IF
      ELSE IF (beta .GT. 0.) THEN
        CALL SPLIT_P_GRAD_FWD(u, v, pkc, gz, du, dv, delp, pk3, beta_d, &
&                       dt, ng, gridstruct, bd, npx, npy, npz, &
&                       flagstruct%use_logp)
        CALL PUSHCONTROL(3,2)
      ELSE IF (beta .LT. -0.1) THEN
        CALL ONE_GRAD_P_FWD(u, v, pkc, gz, divg2, delp, dt, ng, &
&                     gridstruct, bd, npx, npy, npz, ptop, hydrostatic, &
&                     flagstruct%a2b_ord, flagstruct%d_ext)
        CALL PUSHCONTROL(3,3)
      ELSE
        CALL NH_P_GRAD_FWD(u, v, pkc, gz, delp, pk3, dt, ng, gridstruct&
&                    , bd, npx, npy, npz, flagstruct%use_logp)
        CALL PUSHCONTROL(3,4)
      END IF
! Inline Rayleigh friction here?
!-------------------------------------------------------------------------------------------------------
      IF (flagstruct%breed_vortex_inline) THEN
        IF (.NOT.hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pkz,cappa,rdg,delp,delz,pt,k1k)
          DO k=1,npz
            DO j=js,je
              DO i=is,ie
! Note: pt at this stage is Theta_m
                CALL PUSHREALARRAY(pkz(i, j, k))
                pkz(i, j, k) = EXP(k1k*LOG(rdg*delp(i, j, k)/delz(i, j, &
&                 k)*pt(i, j, k)))
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
!-------------------------------------------------------------------------------------------------------
      IF (it .EQ. n_split .AND. gridstruct%grid_type .LT. 4 .AND. (.NOT.&
&         gridstruct%nested)) THEN
! Prevent accumulation of rounding errors at overlapped domain edges:
        CALL PUSHREALARRAY(v, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1)*npz)
        CALL PUSHREALARRAY(u, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2)*npz)
        CALL MPP_GET_BOUNDARY(u, v, domain, ebuffery=ebuffer, nbufferx=&
&                       nbuffer, gridtype=dgrid_ne)
!$OMP parallel do default(none) shared(is,ie,js,je,npz,u,nbuffer,v,ebuffer)
        DO k=1,npz
          DO i=is,ie
            u(i, je+1, k) = nbuffer(i-is+1, k)
          END DO
          DO j=js,je
            v(ie+1, j, k) = ebuffer(j-js+1, k)
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (it .NE. n_split) THEN
        CALL PUSHREALARRAY(v, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1)*npz)
        CALL PUSHREALARRAY(u, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2)*npz)
        CALL START_GROUP_HALO_UPDATE(i_pack(8), u, v, domain, gridtype=&
&                              dgrid_ne)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%nested) neststruct%nest_timestep = neststruct%&
&         nest_timestep + 1
      IF (hydrostatic .AND. last_step) THEN
        IF (flagstruct%use_old_omega) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga,pe,pem,rdt)
          DO k=1,npz
            DO j=js,je
              DO i=is,ie
                CALL PUSHREALARRAY(omga(i, j, k))
                omga(i, j, k) = (pe(i, k+1, j)-pem(i, k+1, j))*rdt
              END DO
            END DO
          END DO
!------------------------------
! Compute the "advective term"
!------------------------------
          CALL ADV_PE_FWD(ua, va, pem, omga, gridstruct, bd, npx, npy, &
&                   npz, ng)
          CALL PUSHCONTROL(1,0)
        ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga) private(om2d)
          DO j=js,je
            DO k=1,npz
              DO i=is,ie
                om2d(i, k) = omga(i, j, k)
              END DO
            END DO
            DO k=2,npz
              DO i=is,ie
                om2d(i, k) = om2d(i, k-1) + omga(i, j, k)
              END DO
            END DO
            DO k=2,npz
              DO i=is,ie
                CALL PUSHREALARRAY(omga(i, j, k))
                omga(i, j, k) = om2d(i, k)
              END DO
            END DO
          END DO
          CALL PUSHCONTROL(1,1)
        END IF
        IF (idiag%id_ws .GT. 0 .AND. hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ws,delz,delp,omga)
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(ws(i, j))
              ws(i, j) = delz(i, j, npz)/delp(i, j, npz)*omga(i, j, npz)
            END DO
          END DO
          CALL PUSHCONTROL(2,0)
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,2)
      END IF
      IF (gridstruct%nested) THEN
        IF (.NOT.hydrostatic) THEN
          arg1 = split_timestep_bc + 1
          arg2 = REAL(n_split*flagstruct%k_split)
          CALL PUSHREALARRAY(w, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz&
&                      )
          CALL NESTED_GRID_BC_APPLY_INTT(w, 0, 0, npx, npy, npz, bd, &
&                                  split_timestep_bc + 1, REAL(n_split*&
&                                  flagstruct%k_split), neststruct%w_bc&
&                                  , bctype=neststruct%nestbctype)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        arg1 = split_timestep_bc + 1
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL PUSHREALARRAY(u, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2)*npz)
        CALL NESTED_GRID_BC_APPLY_INTT(u, 0, 1, npx, npy, npz, bd, &
&                                split_timestep_bc + 1, REAL(n_split*&
&                                flagstruct%k_split), neststruct%u_bc, &
&                                bctype=neststruct%nestbctype)
        arg1 = split_timestep_bc + 1
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL PUSHREALARRAY(v, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1)*npz)
        CALL NESTED_GRID_BC_APPLY_INTT(v, 1, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 1, REAL(n_split*&
&                                flagstruct%k_split), neststruct%v_bc, &
&                                bctype=neststruct%nestbctype)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
!-----------------------------------------------------
! time split loop
!-----------------------------------------------------
    IF (nq .GT. 0 .AND. (.NOT.flagstruct%inline_q)) THEN
      CALL PUSHREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz*nq)
      CALL START_GROUP_HALO_UPDATE(i_pack(10), q, domain)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (flagstruct%fv_debug) THEN
      res0 = IS_MASTER()
      IF (res0) THEN
        CALL PUSHCONTROL(1,0)
        WRITE(*, *) 'End of n_split loop'
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (n_con .NE. 0 .AND. flagstruct%d_con .GT. 1.e-5) THEN
      IF (3 .GT. flagstruct%nord + 1) THEN
        nf_ke = flagstruct%nord + 1
      ELSE
        nf_ke = 3
      END IF
      arg11 = cnst_0p20*gridstruct%da_min
      CALL DEL2_CUBED_FWD(heat_source, arg11, gridstruct, domain, npx, &
&                   npy, npz, nf_ke, bd)
! Note: pt here is cp*(Virtual_Temperature/pkz)
      IF (hydrostatic) THEN
!
! del(Cp*T) = - del(KE)
!
!$OMP parallel do default(none) shared(flagstruct,is,ie,js,je,n_con,pt,heat_source,delp,pkz,bdt) &
!$OMP                          private(dtmp)
        DO j=js,je
! n_con is usually less than 3;
          DO k=1,n_con
            IF (k .LT. 3) THEN
              DO i=is,ie
                CALL PUSHREALARRAY(pt(i, j, k))
                pt(i, j, k) = pt(i, j, k) + heat_source(i, j, k)/(cp_air&
&                 *delp(i, j, k)*pkz(i, j, k))
              END DO
              CALL PUSHCONTROL(1,1)
            ELSE
              DO i=is,ie
                dtmp = heat_source(i, j, k)/(cp_air*delp(i, j, k))
                IF (bdt .GE. 0.) THEN
                  abs0 = bdt
                ELSE
                  abs0 = -bdt
                END IF
                x1 = abs0*flagstruct%delt_max
                IF (dtmp .GE. 0.) THEN
                  y1 = dtmp
                  CALL PUSHCONTROL(1,0)
                ELSE
                  y1 = -dtmp
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
                CALL PUSHREALARRAY(pt(i, j, k))
                pt(i, j, k) = pt(i, j, k) + SIGN(min1, dtmp)/pkz(i, j, k&
&                 )
              END DO
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
        END DO
        CALL PUSHREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL PUSHREALARRAY(d2_divg)
        CALL PUSHINTEGER(ms)
        CALL PUSHINTEGER(hord_v_pert)
        CALL PUSHINTEGER(hord_t_pert)
        CALL PUSHREALARRAY(dp_ref, npz)
        CALL PUSHINTEGER(hord_p_pert)
        CALL PUSHINTEGER(je)
        CALL PUSHREALARRAY(heat_source, (bd%ied-bd%isd+1)*(bd%jed-bd%&
&                     jsd+1)*npz)
        CALL PUSHINTEGER(nord_w)
        CALL PUSHREALARRAY(min1)
        CALL PUSHINTEGER(nord_v, npz + 1)
        CALL PUSHREALARRAY(damp_vt, npz + 1)
        CALL PUSHINTEGER(nord_t)
        CALL PUSHREALARRAY(ws3, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL PUSHREALARRAY(damp_t_pert)
        CALL PUSHINTEGER(nord_k)
        CALL PUSHINTEGER(nord_k_pert)
        CALL PUSHREALARRAY(d2_divg_pert)
        CALL PUSHREALARRAY(rdt)
        CALL PUSHINTEGER(is)
        CALL PUSHREALARRAY(damp_vt_pert, npz + 1)
        CALL PUSHREALARRAY(rdg)
        CALL PUSHINTEGER(ie)
        CALL PUSHREALARRAY(k1k)
        CALL PUSHREALARRAY(dt2)
        CALL PUSHINTEGER(n_con)
        CALL PUSHREALARRAY(beta_d)
        CALL PUSHINTEGER(hord_m_pert)
        CALL PUSHINTEGER(nord_v_pert, npz + 1)
        CALL PUSHINTEGER(nord_t_pert)
        CALL PUSHREALARRAY(dt)
        CALL PUSHINTEGER(js)
        CALL PUSHCONTROL(2,1)
      ELSE
!$OMP parallel do default(none) shared(flagstruct,is,ie,js,je,n_con,pkz,cappa,rdg,delp,delz,pt, &
!$OMP                                  heat_source,k1k,cv_air,bdt) &
!$OMP                          private(dtmp, delt)
        DO k=1,n_con
          IF (bdt*flagstruct%delt_max .GE. 0.) THEN
            delt = bdt*flagstruct%delt_max
          ELSE
            delt = -(bdt*flagstruct%delt_max)
          END IF
! Sponge layers:
!         if ( k == 1 ) delt = 2.0*delt
!         if ( k == 2 ) delt = 1.5*delt
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(pkz(i, j, k))
              pkz(i, j, k) = EXP(k1k*LOG(rdg*delp(i, j, k)/delz(i, j, k)&
&               *pt(i, j, k)))
              dtmp = heat_source(i, j, k)/(cv_air*delp(i, j, k))
              IF (dtmp .GE. 0.) THEN
                y2 = dtmp
                CALL PUSHCONTROL(1,0)
              ELSE
                y2 = -dtmp
                CALL PUSHCONTROL(1,1)
              END IF
              IF (delt .GT. y2) THEN
                CALL PUSHREALARRAY(min2)
                min2 = y2
                CALL PUSHCONTROL(1,0)
              ELSE
                CALL PUSHREALARRAY(min2)
                min2 = delt
                CALL PUSHCONTROL(1,1)
              END IF
              CALL PUSHREALARRAY(pt(i, j, k))
              pt(i, j, k) = pt(i, j, k) + SIGN(min2, dtmp)/pkz(i, j, k)
            END DO
          END DO
        END DO
        CALL PUSHREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL PUSHREALARRAY(d2_divg)
        CALL PUSHINTEGER(ms)
        CALL PUSHINTEGER(hord_v_pert)
        CALL PUSHINTEGER(hord_t_pert)
        CALL PUSHREALARRAY(dp_ref, npz)
        CALL PUSHINTEGER(hord_p_pert)
        CALL PUSHINTEGER(je)
        CALL PUSHREALARRAY(heat_source, (bd%ied-bd%isd+1)*(bd%jed-bd%&
&                     jsd+1)*npz)
        CALL PUSHREALARRAY(min2)
        CALL PUSHINTEGER(nord_w)
        CALL PUSHINTEGER(nord_v, npz + 1)
        CALL PUSHREALARRAY(damp_vt, npz + 1)
        CALL PUSHINTEGER(nord_t)
        CALL PUSHREALARRAY(ws3, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL PUSHREALARRAY(damp_t_pert)
        CALL PUSHINTEGER(nord_k)
        CALL PUSHINTEGER(nord_k_pert)
        CALL PUSHREALARRAY(d2_divg_pert)
        CALL PUSHREALARRAY(rdt)
        CALL PUSHINTEGER(is)
        CALL PUSHREALARRAY(damp_vt_pert, npz + 1)
        CALL PUSHREALARRAY(cv_air)
        CALL PUSHREALARRAY(rdg)
        CALL PUSHINTEGER(ie)
        CALL PUSHREALARRAY(k1k)
        CALL PUSHREALARRAY(dt2)
        CALL PUSHINTEGER(n_con)
        CALL PUSHREALARRAY(beta_d)
        CALL PUSHINTEGER(hord_m_pert)
        CALL PUSHINTEGER(nord_v_pert, npz + 1)
        CALL PUSHINTEGER(nord_t_pert)
        CALL PUSHREALARRAY(dt)
        CALL PUSHINTEGER(js)
        CALL PUSHCONTROL(2,2)
      END IF
    ELSE
      CALL PUSHREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(d2_divg)
      CALL PUSHINTEGER(ms)
      CALL PUSHINTEGER(hord_v_pert)
      CALL PUSHINTEGER(hord_t_pert)
      CALL PUSHREALARRAY(dp_ref, npz)
      CALL PUSHINTEGER(hord_p_pert)
      CALL PUSHINTEGER(je)
      CALL PUSHINTEGER(nord_w)
      CALL PUSHINTEGER(nord_v, npz + 1)
      CALL PUSHREALARRAY(damp_vt, npz + 1)
      CALL PUSHINTEGER(nord_t)
      CALL PUSHREALARRAY(ws3, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(damp_t_pert)
      CALL PUSHINTEGER(nord_k)
      CALL PUSHINTEGER(nord_k_pert)
      CALL PUSHREALARRAY(d2_divg_pert)
      CALL PUSHREALARRAY(rdt)
      CALL PUSHINTEGER(is)
      CALL PUSHREALARRAY(damp_vt_pert, npz + 1)
      CALL PUSHREALARRAY(rdg)
      CALL PUSHINTEGER(ie)
      CALL PUSHREALARRAY(k1k)
      CALL PUSHREALARRAY(dt2)
      CALL PUSHREALARRAY(beta_d)
      CALL PUSHINTEGER(hord_m_pert)
      CALL PUSHINTEGER(nord_v_pert, npz + 1)
      CALL PUSHINTEGER(nord_t_pert)
      CALL PUSHREALARRAY(dt)
      CALL PUSHINTEGER(js)
      CALL PUSHCONTROL(2,0)
    END IF
  END SUBROUTINE DYN_CORE_FWD
!  Differentiation of dyn_core in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
!d.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_m
!od.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mi
!x_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_
!Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4
! fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rem
!ap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv
!_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters f
!v_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_re
!start_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid
!_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_m
!od.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod
!.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.
!nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a
!2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v_f
!b sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_m
!od.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_
!mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pk3 xfx ws peln q gz du u dv
!                v w delp ua uc ptc mfx delz mfy omga ut divgd
!                pkc delpc va vc yfx pkz pe vt pk zh pt cx cy dpx
!                crx cry
!   with respect to varying inputs: pk3 xfx ws peln q gz du u dv
!                v w delp ua uc ptc delz omga ut divgd pkc delpc
!                va vc yfx pkz pe vt pk zh pt dpx crx cry
!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
  SUBROUTINE DYN_CORE_BWD(npx, npy, npz, ng, sphum, nq, bdt, n_split, &
&   zvir, cp, akap, cappa, grav, hydrostatic, u, u_ad, v, v_ad, w, w_ad&
&   , delz, delz_ad, pt, pt_ad, q, q_ad, delp, delp_ad, pe, pe_ad, pk, &
&   pk_ad, phis, ws, ws_ad, omga, omga_ad, ptop, pfull, ua, ua_ad, va, &
&   va_ad, uc, uc_ad, vc, vc_ad, mfx, mfx_ad, mfy, mfy_ad, cx, cx_ad, cy&
&   , cy_ad, pkz, pkz_ad, peln, peln_ad, q_con, ak, bk, dpx, dpx_ad, ks&
&   , gridstruct, flagstruct, flagstructp, neststruct, idiag, bd, domain&
&   , init_step, i_pack, end_step, gz, gz_ad, pkc, pkc_ad, ptc, ptc_ad, &
&   crx, crx_ad, xfx, xfx_ad, cry, cry_ad, yfx, yfx_ad, divgd, divgd_ad&
&   , delpc, delpc_ad, ut, ut_ad, vt, vt_ad, zh, zh_ad, pk3, pk3_ad, du&
&   , du_ad, dv, dv_ad, time_total)
    IMPLICIT NONE
! end init_step
! Start of the big dynamic time stepping
!allocate(    gz(isd:ied, jsd:jed ,npz+1) )
!  call init_ijk_mem(isd,ied, jsd,jed, npz+1, gz, huge_r)
!allocate(   pkc(isd:ied, jsd:jed ,npz+1) )
!allocate(   ptc(isd:ied, jsd:jed ,npz ) )
!allocate( crx(is :ie+1, jsd:jed,  npz) )
!allocate( xfx(is :ie+1, jsd:jed,  npz) )
!allocate( cry(isd:ied,  js :je+1, npz) )
!allocate( yfx(isd:ied,  js :je+1, npz) )
!allocate( divgd(isd:ied+1,jsd:jed+1,npz) )
!allocate( delpc(isd:ied, jsd:jed  ,npz  ) )
!                    call init_ijk_mem(isd,ied, jsd,jed, npz, delpc, 0.)
!allocate( ut(isd:ied, jsd:jed, npz) )
!                    call init_ijk_mem(isd,ied, jsd,jed, npz, ut, 0.)
!allocate( vt(isd:ied, jsd:jed, npz) )
!                    call init_ijk_mem(isd,ied, jsd,jed, npz, vt, 0.)
!allocate( zh(isd:ied, jsd:jed, npz+1) )
!              call init_ijk_mem(isd,ied, jsd,jed, npz+1, zh, huge_r )
!allocate ( pk3(isd:ied,jsd:jed,npz+1) )
!call init_ijk_mem(isd,ied, jsd,jed, npz+1, pk3, huge_r )
!if (allocated(heat_source)) deallocate( heat_source ) !If ncon == 0 but d_con > 1.e-5, this would not be deallocated in earlier 
!versions of the code
!deallocate(    gz )
!deallocate(   ptc )
!deallocate(   crx )
!deallocate(   xfx )
!deallocate(   cry )
!deallocate(   yfx )
!deallocate( divgd )
!deallocate(   pkc )
!deallocate( delpc )
!if( allocated(ut))   deallocate( ut )
!if( allocated(vt))   deallocate( vt )
!if ( allocated (du) ) deallocate( du )
!if ( allocated (dv) ) deallocate( dv )
!if ( .not. hydrostatic ) then
!     deallocate( zh )
!     if( allocated(pk3) )   deallocate ( pk3 )
!endif
!if( allocated(pem) )   deallocate ( pem )
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: ng, nq, sphum
    INTEGER, INTENT(IN) :: n_split
    REAL, INTENT(IN) :: bdt
    REAL, INTENT(IN) :: zvir, cp, akap, grav
    REAL, INTENT(IN) :: ptop
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: init_step, end_step
    REAL, INTENT(IN) :: pfull(npz)
    REAL, DIMENSION(npz+1), INTENT(IN) :: ak, bk
    INTEGER, INTENT(IN) :: ks
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: i_pack(*)
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
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delz_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cappa(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
    REAL, INTENT(IN), OPTIONAL :: time_total
    REAL, INTENT(INOUT) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(INOUT) :: pe_ad(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1&
&   )
    REAL, INTENT(INOUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: peln_ad(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL, INTENT(INOUT) :: pk_ad(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL(kind=8), INTENT(INOUT) :: dpx(bd%is:bd%ie, bd%js:bd%je)
    REAL(kind=8), INTENT(INOUT) :: dpx_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL, PARAMETER :: near0=1.e-8
    REAL, PARAMETER :: huge_r=1.e8
    REAL :: ws(bd%is:bd%ie, bd%js:bd%je)
    REAL :: ws_ad(bd%is:bd%ie, bd%js:bd%je)
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
    REAL, INTENT(INOUT) :: q_con(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_ad(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_ad(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je, npz), INTENT(INOUT) :: pkz
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je, npz), INTENT(INOUT) :: &
&   pkz_ad
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    TYPE(FV_FLAGS_PERT_TYPE), INTENT(IN), TARGET :: flagstructp
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(FV_DIAG_TYPE), INTENT(IN) :: idiag
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL :: pem(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1), heat_source(bd&
&   %isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: pem_ad(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1), &
&   heat_source_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ws3, z_rat
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ws3_ad, z_rat_ad
    REAL :: dp_ref(npz)
    REAL :: zs(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: p1d(bd%is:bd%ie)
    REAL :: om2d(bd%is:bd%ie, npz)
    REAL :: om2d_ad(bd%is:bd%ie, npz)
    REAL :: wbuffer(npy+2, npz)
    REAL :: ebuffer(npy+2, npz)
    REAL :: ebuffer_ad(npy+2, npz)
    REAL :: nbuffer(npx+2, npz)
    REAL :: nbuffer_ad(npx+2, npz)
    REAL :: sbuffer(npx+2, npz)
    REAL :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: divg2_ad(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: fz(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: heat_s(bd%is:bd%ie, bd%js:bd%je)
    REAL :: heat_s_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL :: damp_vt(npz+1)
    INTEGER :: nord_v(npz+1)
    INTEGER :: hord_m, hord_v, hord_t, hord_p
    INTEGER :: nord_k, nord_w, nord_t
    INTEGER :: ms
    INTEGER :: hord_m_pert, hord_v_pert, hord_t_pert, hord_p_pert
    INTEGER :: nord_k_pert, nord_w_pert, nord_t_pert, nord_v_pert(npz+1)
    REAL :: d2_divg_pert, damp_vt_pert(npz+1), damp_w_pert, damp_t_pert
    INTEGER :: i, j, k, it, iq, n_con, nf_ke
    INTEGER :: iep1, jep1
    REAL :: beta, beta_d, d_con_k, damp_w, damp_t, kgb, cv_air
    REAL :: dt, dt2, rdt
    REAL :: d2_divg
    REAL :: k1k, rdg, dtmp, delt
    REAL :: dtmp_ad
    LOGICAL :: last_step, remap_step
    LOGICAL :: used
    REAL :: split_timestep_bc
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pkc(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pkc_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: ptc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ptc_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: crx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: xfx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cry(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cry_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: yfx_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: divgd(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: divgd_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, &
&   npz)
    REAL, INTENT(INOUT) :: delpc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delpc_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ut(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ut_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vt_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: zh(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: zh_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk3(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk3_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: du_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dv_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    INTRINSIC LOG
    INTRINSIC REAL
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC EXP
    INTRINSIC ABS
    INTRINSIC SIGN
    INTEGER :: max1
    INTEGER :: max2
    REAL :: min1
    REAL :: min1_ad
    REAL :: min2
    REAL :: min2_ad
    REAL :: abs0
    REAL :: arg1
    REAL :: arg2
    LOGICAL :: arg10
    REAL*8 :: arg11
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp2
    REAL :: temp3
    REAL :: temp_ad4
    REAL :: temp4
    REAL :: y1_ad
    REAL :: temp_ad5
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp8
    REAL :: temp_ad6
    REAL :: y2_ad
    REAL :: temp_ad7
    INTEGER :: branch
    REAL :: x1
    REAL :: y2
    REAL :: y1

    pem = 0.0
    heat_source = 0.0
    ws3 = 0.0
    z_rat = 0.0
    dp_ref = 0.0
    zs = 0.0
    p1d = 0.0
    om2d = 0.0
    wbuffer = 0.0
    ebuffer = 0.0
    nbuffer = 0.0
    sbuffer = 0.0
    divg2 = 0.0
    wk = 0.0
    fz = 0.0
    heat_s = 0.0
    damp_vt = 0.0
    d2_divg_pert = 0.0
    damp_vt_pert = 0.0
    damp_w_pert = 0.0
    damp_t_pert = 0.0
    beta = 0.0
    beta_d = 0.0
    d_con_k = 0.0
    damp_w = 0.0
    damp_t = 0.0
    kgb = 0.0
    cv_air = 0.0
    dt = 0.0
    dt2 = 0.0
    rdt = 0.0
    d2_divg = 0.0
    k1k = 0.0
    rdg = 0.0
    dtmp = 0.0
    delt = 0.0
    split_timestep_bc = 0.0
    min1 = 0.0
    min2 = 0.0
    abs0 = 0.0
    arg1 = 0.0
    arg2 = 0.0
    arg10 = 0.0
    x1 = 0.0
    y2 = 0.0
    y1 = 0.0
    nord_v = 0
    hord_m = 0
    hord_v = 0
    hord_t = 0
    hord_p = 0
    nord_k = 0
    nord_w = 0
    nord_t = 0
    ms = 0
    hord_m_pert = 0
    hord_v_pert = 0
    hord_t_pert = 0
    hord_p_pert = 0
    nord_k_pert = 0
    nord_w_pert = 0
    nord_t_pert = 0
    nord_v_pert = 0
    i = 0
    j = 0
    k = 0
    it = 0
    iq = 0
    n_con = 0
    nf_ke = 0
    iep1 = 0
    jep1 = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    branch = 0

    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(js)
      CALL POPREALARRAY(dt)
      CALL POPINTEGER(nord_t_pert)
      CALL POPINTEGER(nord_v_pert, npz + 1)
      CALL POPINTEGER(hord_m_pert)
      CALL POPREALARRAY(beta_d)
      CALL POPREALARRAY(dt2)
      CALL POPREALARRAY(k1k)
      CALL POPINTEGER(ie)
      CALL POPREALARRAY(rdg)
      CALL POPREALARRAY(damp_vt_pert, npz + 1)
      CALL POPINTEGER(is)
      CALL POPREALARRAY(rdt)
      CALL POPREALARRAY(d2_divg_pert)
      CALL POPINTEGER(nord_k_pert)
      CALL POPINTEGER(nord_k)
      CALL POPREALARRAY(damp_t_pert)
      CALL POPREALARRAY(ws3, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL POPINTEGER(nord_t)
      CALL POPREALARRAY(damp_vt, npz + 1)
      CALL POPINTEGER(nord_v, npz + 1)
      CALL POPINTEGER(nord_w)
      CALL POPINTEGER(je)
      CALL POPINTEGER(hord_p_pert)
      CALL POPREALARRAY(dp_ref, npz)
      CALL POPINTEGER(hord_t_pert)
      CALL POPINTEGER(hord_v_pert)
      CALL POPINTEGER(ms)
      CALL POPREALARRAY(d2_divg)
      CALL POPREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      heat_source_ad = 0.0
    ELSE
      IF (branch .EQ. 1) THEN
        CALL POPINTEGER(js)
        CALL POPREALARRAY(dt)
        CALL POPINTEGER(nord_t_pert)
        CALL POPINTEGER(nord_v_pert, npz + 1)
        CALL POPINTEGER(hord_m_pert)
        CALL POPREALARRAY(beta_d)
        CALL POPINTEGER(n_con)
        CALL POPREALARRAY(dt2)
        CALL POPREALARRAY(k1k)
        CALL POPINTEGER(ie)
        CALL POPREALARRAY(rdg)
        CALL POPREALARRAY(damp_vt_pert, npz + 1)
        CALL POPINTEGER(is)
        CALL POPREALARRAY(rdt)
        CALL POPREALARRAY(d2_divg_pert)
        CALL POPINTEGER(nord_k_pert)
        CALL POPINTEGER(nord_k)
        CALL POPREALARRAY(damp_t_pert)
        CALL POPREALARRAY(ws3, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL POPINTEGER(nord_t)
        CALL POPREALARRAY(damp_vt, npz + 1)
        CALL POPINTEGER(nord_v, npz + 1)
        CALL POPREALARRAY(min1)
        CALL POPINTEGER(nord_w)
        CALL POPREALARRAY(heat_source, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd&
&                    +1)*npz)
        CALL POPINTEGER(je)
        CALL POPINTEGER(hord_p_pert)
        CALL POPREALARRAY(dp_ref, npz)
        CALL POPINTEGER(hord_t_pert)
        CALL POPINTEGER(hord_v_pert)
        CALL POPINTEGER(ms)
        CALL POPREALARRAY(d2_divg)
        CALL POPREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        heat_source_ad = 0.0
        DO j=je,js,-1
          DO k=n_con,1,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              DO i=ie,is,-1
                dtmp = heat_source(i, j, k)/(cp_air*delp(i, j, k))
                CALL POPREALARRAY(pt(i, j, k))
                temp_ad5 = pt_ad(i, j, k)/pkz(i, j, k)
                min1_ad = SIGN(1.d0, min1*dtmp)*temp_ad5
                pkz_ad(i, j, k) = pkz_ad(i, j, k) - SIGN(min1, dtmp)*&
&                 temp_ad5/pkz(i, j, k)
                CALL POPCONTROL(1,branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY(min1)
                  y1_ad = min1_ad
                ELSE
                  CALL POPREALARRAY(min1)
                  y1_ad = 0.0
                END IF
                CALL POPCONTROL(1,branch)
                IF (branch .EQ. 0) THEN
                  dtmp_ad = y1_ad
                ELSE
                  dtmp_ad = -y1_ad
                END IF
                temp4 = cp_air*delp(i, j, k)
                heat_source_ad(i, j, k) = heat_source_ad(i, j, k) + &
&                 dtmp_ad/temp4
                delp_ad(i, j, k) = delp_ad(i, j, k) - heat_source(i, j, &
&                 k)*cp_air*dtmp_ad/temp4**2
              END DO
            ELSE
              DO i=ie,is,-1
                CALL POPREALARRAY(pt(i, j, k))
                temp3 = cp_air*delp(i, j, k)
                temp2 = temp3*pkz(i, j, k)
                temp_ad4 = -(heat_source(i, j, k)*pt_ad(i, j, k)/temp2**&
&                 2)
                heat_source_ad(i, j, k) = heat_source_ad(i, j, k) + &
&                 pt_ad(i, j, k)/temp2
                delp_ad(i, j, k) = delp_ad(i, j, k) + pkz(i, j, k)*&
&                 cp_air*temp_ad4
                pkz_ad(i, j, k) = pkz_ad(i, j, k) + temp3*temp_ad4
              END DO
            END IF
          END DO
        END DO
      ELSE
        CALL POPINTEGER(js)
        CALL POPREALARRAY(dt)
        CALL POPINTEGER(nord_t_pert)
        CALL POPINTEGER(nord_v_pert, npz + 1)
        CALL POPINTEGER(hord_m_pert)
        CALL POPREALARRAY(beta_d)
        CALL POPINTEGER(n_con)
        CALL POPREALARRAY(dt2)
        CALL POPREALARRAY(k1k)
        CALL POPINTEGER(ie)
        CALL POPREALARRAY(rdg)
        CALL POPREALARRAY(cv_air)
        CALL POPREALARRAY(damp_vt_pert, npz + 1)
        CALL POPINTEGER(is)
        CALL POPREALARRAY(rdt)
        CALL POPREALARRAY(d2_divg_pert)
        CALL POPINTEGER(nord_k_pert)
        CALL POPINTEGER(nord_k)
        CALL POPREALARRAY(damp_t_pert)
        CALL POPREALARRAY(ws3, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL POPINTEGER(nord_t)
        CALL POPREALARRAY(damp_vt, npz + 1)
        CALL POPINTEGER(nord_v, npz + 1)
        CALL POPINTEGER(nord_w)
        CALL POPREALARRAY(min2)
        CALL POPREALARRAY(heat_source, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd&
&                    +1)*npz)
        CALL POPINTEGER(je)
        CALL POPINTEGER(hord_p_pert)
        CALL POPREALARRAY(dp_ref, npz)
        CALL POPINTEGER(hord_t_pert)
        CALL POPINTEGER(hord_v_pert)
        CALL POPINTEGER(ms)
        CALL POPREALARRAY(d2_divg)
        CALL POPREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        heat_source_ad = 0.0
        DO k=n_con,1,-1
          DO j=je,js,-1
            DO i=ie,is,-1
              dtmp = heat_source(i, j, k)/(cv_air*delp(i, j, k))
              CALL POPREALARRAY(pt(i, j, k))
              temp_ad7 = pt_ad(i, j, k)/pkz(i, j, k)
              min2_ad = SIGN(1.d0, min2*dtmp)*temp_ad7
              pkz_ad(i, j, k) = pkz_ad(i, j, k) - SIGN(min2, dtmp)*&
&               temp_ad7/pkz(i, j, k)
              CALL POPCONTROL(1,branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY(min2)
                y2_ad = min2_ad
              ELSE
                CALL POPREALARRAY(min2)
                y2_ad = 0.0
              END IF
              CALL POPCONTROL(1,branch)
              IF (branch .EQ. 0) THEN
                dtmp_ad = y2_ad
              ELSE
                dtmp_ad = -y2_ad
              END IF
              temp7 = delz(i, j, k)
              temp6 = delp(i, j, k)*pt(i, j, k)
              temp5 = temp6/temp7
              temp_ad6 = k1k*EXP(k1k*LOG(rdg*temp5))*pkz_ad(i, j, k)/(&
&               temp5*temp7)
              temp8 = cv_air*delp(i, j, k)
              heat_source_ad(i, j, k) = heat_source_ad(i, j, k) + &
&               dtmp_ad/temp8
              delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*temp_ad6&
&               - heat_source(i, j, k)*cv_air*dtmp_ad/temp8**2
              CALL POPREALARRAY(pkz(i, j, k))
              pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*temp_ad6
              delz_ad(i, j, k) = delz_ad(i, j, k) - temp5*temp_ad6
              pkz_ad(i, j, k) = 0.0
            END DO
          END DO
        END DO
      END IF
      arg11 = cnst_0p20*gridstruct%da_min
      CALL DEL2_CUBED_BWD(heat_source, heat_source_ad, arg11, gridstruct&
&                   , domain, npx, npy, npz, nf_ke, bd)
    END IF
    CALL POPCONTROL(1,branch)
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz*nq)
      CALL START_GROUP_HALO_UPDATE_ADM(i_pack(10), q, q_ad, domain)
    END IF
    jep1 = je + 1
    jsd = bd%jsd
    ied = bd%ied
    iep1 = ie + 1
    isd = bd%isd
    jed = bd%jed
    om2d_ad = 0.0
    pem_ad = 0.0
    ws3_ad = 0.0
    z_rat_ad = 0.0
    heat_s_ad = 0.0
    wk_ad = 0.0
    divg2_ad = 0.0
    DO it=n_split,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        CALL POPREALARRAY(v, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1)*npz)
        CALL NESTED_GRID_BC_APPLY_INTT_ADM(v, v_ad, 1, 0, npx, npy, npz&
&                                    , bd, arg1, arg2, neststruct%v_bc, &
&                                    neststruct%nestbctype)
        CALL POPREALARRAY(u, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2)*npz)
        CALL NESTED_GRID_BC_APPLY_INTT_ADM(u, u_ad, 0, 1, npx, npy, npz&
&                                    , bd, arg1, arg2, neststruct%u_bc, &
&                                    neststruct%nestbctype)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(w, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
          CALL NESTED_GRID_BC_APPLY_INTT_ADM(w, w_ad, 0, 0, npx, npy, &
&                                      npz, bd, arg1, arg2, neststruct%&
&                                      w_bc, neststruct%nestbctype)
        END IF
      END IF
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(ws(i, j))
            temp_ad3 = ws_ad(i, j)/delp(i, j, npz)
            delz_ad(i, j, npz) = delz_ad(i, j, npz) + omga(i, j, npz)*&
&             temp_ad3
            omga_ad(i, j, npz) = omga_ad(i, j, npz) + delz(i, j, npz)*&
&             temp_ad3
            delp_ad(i, j, npz) = delp_ad(i, j, npz) - delz(i, j, npz)*&
&             omga(i, j, npz)*temp_ad3/delp(i, j, npz)
            ws_ad(i, j) = 0.0
          END DO
        END DO
      ELSE IF (branch .NE. 1) THEN
        GOTO 100
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL ADV_PE_BWD(ua, ua_ad, va, va_ad, pem, pem_ad, omga, omga_ad&
&                 , gridstruct, bd, npx, npy, npz, ng)
        DO k=npz,1,-1
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(omga(i, j, k))
              pe_ad(i, k+1, j) = pe_ad(i, k+1, j) + rdt*omga_ad(i, j, k)
              pem_ad(i, k+1, j) = pem_ad(i, k+1, j) - rdt*omga_ad(i, j, &
&               k)
              omga_ad(i, j, k) = 0.0
            END DO
          END DO
        END DO
      ELSE
        DO j=je,js,-1
          DO k=npz,2,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(omga(i, j, k))
              om2d_ad(i, k) = om2d_ad(i, k) + omga_ad(i, j, k)
              omga_ad(i, j, k) = 0.0
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ie,is,-1
              om2d_ad(i, k-1) = om2d_ad(i, k-1) + om2d_ad(i, k)
              omga_ad(i, j, k) = omga_ad(i, j, k) + om2d_ad(i, k)
              om2d_ad(i, k) = 0.0
            END DO
          END DO
          DO k=npz,1,-1
            DO i=ie,is,-1
              omga_ad(i, j, k) = omga_ad(i, j, k) + om2d_ad(i, k)
              om2d_ad(i, k) = 0.0
            END DO
          END DO
        END DO
      END IF
 100  CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(u, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2)*npz)
        CALL POPREALARRAY(v, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1)*npz)
        CALL START_GROUP_HALO_UPDATE_ADM(i_pack(8), u, u_ad, v, v_ad, &
&                                  domain, gridtype=dgrid_ne)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        nbuffer_ad = 0.0
        ebuffer_ad = 0.0
        DO k=npz,1,-1
          DO j=je,js,-1
            ebuffer_ad(j-js+1, k) = ebuffer_ad(j-js+1, k) + v_ad(ie+1, j&
&             , k)
            v_ad(ie+1, j, k) = 0.0
          END DO
          DO i=ie,is,-1
            nbuffer_ad(i-is+1, k) = nbuffer_ad(i-is+1, k) + u_ad(i, je+1&
&             , k)
            u_ad(i, je+1, k) = 0.0
          END DO
        END DO
        CALL POPREALARRAY(u, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2)*npz)
        CALL POPREALARRAY(v, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1)*npz)
        CALL MPP_GET_BOUNDARY_ADM(u, u_ad, v, v_ad, domain, ebuffery=&
&                           ebuffer, ebuffery_ad=ebuffer_ad, nbufferx=&
&                           nbuffer, nbufferx_ad=nbuffer_ad, gridtype=&
&                           dgrid_ne)
      END IF
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        DO k=npz,1,-1
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(pkz(i, j, k))
              temp1 = delz(i, j, k)
              temp0 = delp(i, j, k)*pt(i, j, k)
              temp = temp0/temp1
              temp_ad2 = k1k*EXP(k1k*LOG(rdg*temp))*pkz_ad(i, j, k)/(&
&               temp*temp1)
              delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*temp_ad2
              pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*temp_ad2
              delz_ad(i, j, k) = delz_ad(i, j, k) - temp*temp_ad2
              pkz_ad(i, j, k) = 0.0
            END DO
          END DO
        END DO
      END IF
      CALL POPCONTROL(3,branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          CALL GRAD1_P_UPDATE_BWD(divg2, divg2_ad, u, u_ad, v, v_ad, pkc&
&                           , pkc_ad, gz, gz_ad, du, du_ad, dv, dv_ad, &
&                           dt, ng, gridstruct, bd, npx, npy, npz, ptop&
&                           , beta_d, flagstruct%a2b_ord)
        ELSE
          CALL ONE_GRAD_P_BWD(u, u_ad, v, v_ad, pkc, pkc_ad, gz, gz_ad, &
&                       divg2, divg2_ad, delp, delp_ad, dt, ng, &
&                       gridstruct, bd, npx, npy, npz, ptop, hydrostatic&
&                       , flagstruct%a2b_ord, flagstruct%d_ext)
        END IF
      ELSE IF (branch .EQ. 2) THEN
        CALL SPLIT_P_GRAD_BWD(u, u_ad, v, v_ad, pkc, pkc_ad, gz, gz_ad, &
&                       du, du_ad, dv, dv_ad, delp, delp_ad, pk3, pk3_ad&
&                       , beta_d, dt, ng, gridstruct, bd, npx, npy, npz&
&                       , flagstruct%use_logp)
      ELSE IF (branch .EQ. 3) THEN
        CALL ONE_GRAD_P_BWD(u, u_ad, v, v_ad, pkc, pkc_ad, gz, gz_ad, &
&                     divg2, divg2_ad, delp, delp_ad, dt, ng, gridstruct&
&                     , bd, npx, npy, npz, ptop, hydrostatic, flagstruct&
&                     %a2b_ord, flagstruct%d_ext)
      ELSE
        CALL NH_P_GRAD_BWD(u, u_ad, v, v_ad, pkc, pkc_ad, gz, gz_ad, &
&                    delp, delp_ad, pk3, pk3_ad, dt, ng, gridstruct, bd&
&                    , npx, npy, npz, flagstruct%use_logp)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO k=npz+1,1,-1
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(pk(i, j, k))
              pkc_ad(i, j, k) = pkc_ad(i, j, k) + pk_ad(i, j, k)
              pk_ad(i, j, k) = 0.0
            END DO
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL GEOPK_BWD(ptop, pe, pe_ad, peln, peln_ad, delp, delp_ad, &
&                pkc, pkc_ad, gz, gz_ad, phis, pt, pt_ad, q_con, pkz, &
&                pkz_ad, npz, akap, .false., gridstruct%nested, .true., &
&                npx, npy, flagstruct%a2b_ord, bd)
      ELSE
        DO k=npz+1,1,-1
          DO j=je+2,js-2,-1
            DO i=ie+2,is-2,-1
              CALL POPREALARRAY(gz(i, j, k))
              zh_ad(i, j, k) = zh_ad(i, j, k) + grav*gz_ad(i, j, k)
              gz_ad(i, j, k) = 0.0
            END DO
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL NEST_HALO_NH_BWD(ptop, grav, akap, cp, delp, delp_ad, &
&                         delz, delz_ad, pt, pt_ad, phis, pkc, pkc_ad, &
&                         gz, gz_ad, pk3, pk3_ad, npx, npy, npz, &
&                         gridstruct%nested, .true., .true., .true., bd)
          CALL POPREALARRAY(delz, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*&
&                      npz)
          CALL NESTED_GRID_BC_APPLY_INTT_ADM(delz, delz_ad, 0, 0, npx, &
&                                      npy, npz, bd, arg1, arg2, &
&                                      neststruct%delz_bc, neststruct%&
&                                      nestbctype)
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL PLN_HALO_BWD(is, ie, js, je, isd, ied, jsd, jed, npz, &
&                     ptop, pk3, pk3_ad, delp, delp_ad)
        ELSE
          CALL PK3_HALO_BWD(is, ie, js, je, isd, ied, jsd, jed, npz, &
&                     ptop, akap, pk3, pk3_ad, delp, delp_ad)
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) CALL PE_HALO_BWD(is, ie, js, je, isd, ied, &
&                                     jsd, jed, npz, ptop, pe, pe_ad, &
&                                     delp, delp_ad)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(pkc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(&
&                      npz+1))
          CALL START_GROUP_HALO_UPDATE_ADM(i_pack(5), pkc, pkc_ad, &
&                                    domain, whalo=2, ehalo=2, shalo=2, &
&                                    nhalo=2)
          CALL POPREALARRAY(zh, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(&
&                      npz+1))
          CALL START_GROUP_HALO_UPDATE_ADM(i_pack(4), zh, zh_ad, domain)
        ELSE
          CALL POPREALARRAY(pkc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(&
&                      npz+1))
          CALL START_GROUP_HALO_UPDATE_ADM(i_pack(4), pkc, pkc_ad, &
&                                    domain, complete=.true.)
          CALL POPREALARRAY(zh, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(&
&                      npz+1))
          CALL START_GROUP_HALO_UPDATE_ADM(i_pack(4), zh, zh_ad, domain&
&                                    , complete=.true.)
        END IF
        CALL RIEM_SOLVER3_BWD(flagstruct%m_split, dt, is, ie, js, je, &
&                       npz, ng, isd, ied, jsd, jed, akap, cappa, cp, &
&                       ptop, zs, q_con, w, w_ad, delz, delz_ad, pt, &
&                       pt_ad, delp, delp_ad, zh, zh_ad, pe, pe_ad, pkc&
&                       , pkc_ad, pk3, pk3_ad, pk, pk_ad, peln, peln_ad&
&                       , ws, ws_ad, flagstruct%scale_z, flagstruct%&
&                       p_fac, flagstruct%a_imp, flagstruct%use_logp, &
&                       remap_step, arg10)
        CALL UPDATE_DZ_D_BWD(nord_v, damp_vt, flagstruct%hord_tm, is, ie&
&                      , js, je, npz, ng, npx, npy, gridstruct%area, &
&                      gridstruct%rarea, dp_ref, zs, zh, zh_ad, crx, &
&                      crx_ad, cry, cry_ad, xfx, xfx_ad, yfx, yfx_ad, &
&                      delz, ws, ws_ad, rdt, gridstruct, bd, flagstructp&
&                      %hord_tm_pert)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL NESTED_GRID_BC_APPLY_INTT_ADM(pt, pt_ad, 0, 0, npx, npy, &
&                                    npz, bd, arg1, arg2, neststruct%&
&                                    pt_bc, neststruct%nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT_ADM(delp, delp_ad, 0, 0, npx, npy&
&                                    , npz, bd, arg1, arg2, neststruct%&
&                                    delp_bc, neststruct%nestbctype)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=jep1,js,-1
          DO i=iep1,is,-1
            CALL POPREALARRAY(divg2(i, j))
            temp_ad1 = d2_divg*divg2_ad(i, j)/wk(i, j)
            wk_ad(i, j) = wk_ad(i, j) - divg2(i, j)*temp_ad1/wk(i, j)
            divg2_ad(i, j) = temp_ad1
          END DO
          DO k=npz,2,-1
            DO i=iep1,is,-1
              ptc_ad(i, j, k) = ptc_ad(i, j, k) + wk_ad(i, j) + vt(i, j&
&               , k)*divg2_ad(i, j)
              vt_ad(i, j, k) = vt_ad(i, j, k) + ptc(i, j, k)*divg2_ad(i&
&               , j)
              CALL POPREALARRAY(wk(i, j))
            END DO
          END DO
          DO i=iep1,is,-1
            wk_ad(i, j) = wk_ad(i, j) + vt(i, j, 1)*divg2_ad(i, j)
            vt_ad(i, j, 1) = vt_ad(i, j, 1) + wk(i, j)*divg2_ad(i, j)
            divg2_ad(i, j) = 0.0
            CALL POPREALARRAY(wk(i, j))
            ptc_ad(i, j, 1) = ptc_ad(i, j, 1) + wk_ad(i, j)
            wk_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPREALARRAY(d2_divg)
      ELSE
        divg2_ad = 0.0
      END IF
      CALL POPREALARRAY(pt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE_ADM(i_pack(1), pt, pt_ad, domain, &
&                                complete=.true.)
      CALL POPREALARRAY(delp, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE_ADM(i_pack(1), delp, delp_ad, domain&
&                                , complete=.true.)
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) CALL MIX_DP_BWD(hydrostatic, w, w_ad, delp, &
&                                  delp_ad, pt, pt_ad, npz, ak, bk, &
&                                  .false., flagstruct%fv_debug, bd)
      DO k=npz,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          DO j=je,js,-1
            DO i=ie,is,-1
              heat_s_ad(i, j) = heat_s_ad(i, j) + heat_source_ad(i, j, k&
&               )
            END DO
          END DO
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=jep1,js,-1
            DO i=iep1,is,-1
              CALL POPREALARRAY(ptc(i, j, k))
              wk_ad(i, j) = wk_ad(i, j) + ptc_ad(i, j, k)
              ptc_ad(i, j, k) = 0.0
            END DO
          END DO
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(omga(i, j, k))
              temp_ad = gridstruct%rarea(i, j)*rdt*omga_ad(i, j, k)
              temp_ad0 = omga(i, j, k)*temp_ad
              xfx_ad(i, j, k) = xfx_ad(i, j, k) + temp_ad0
              xfx_ad(i+1, j, k) = xfx_ad(i+1, j, k) - temp_ad0
              yfx_ad(i, j, k) = yfx_ad(i, j, k) + temp_ad0
              yfx_ad(i, j+1, k) = yfx_ad(i, j+1, k) - temp_ad0
              omga_ad(i, j, k) = (xfx(i, j, k)-xfx(i+1, j, k)+yfx(i, j, &
&               k)-yfx(i, j+1, k))*temp_ad
            END DO
          END DO
        END IF
        CALL D_SW_BWD(vt(isd:ied, jsd:jed, k), vt_ad(isd:ied, jsd:jed, k&
&               ), delp(isd:ied, jsd:jed, k), delp_ad(isd:ied, jsd:jed, &
&               k), ptc(isd:ied, jsd:jed, k), ptc_ad(isd:ied, jsd:jed, k&
&               ), pt(isd:ied, jsd:jed, k), pt_ad(isd:ied, jsd:jed, k), &
&               u(isd:ied, jsd:jed+1, k), u_ad(isd:ied, jsd:jed+1, k), v&
&               (isd:ied+1, jsd:jed, k), v_ad(isd:ied+1, jsd:jed, k), w(&
&               isd:ied, jsd:jed, k), w_ad(isd:ied, jsd:jed, k), uc(isd:&
&               ied+1, jsd:jed, k), uc_ad(isd:ied+1, jsd:jed, k), vc(isd&
&               :ied, jsd:jed+1, k), vc_ad(isd:ied, jsd:jed+1, k), ua(&
&               isd:ied, jsd:jed, k), ua_ad(isd:ied, jsd:jed, k), va(isd&
&               :ied, jsd:jed, k), va_ad(isd:ied, jsd:jed, k), divgd(isd&
&               :ied+1, jsd:jed+1, k), divgd_ad(isd:ied+1, jsd:jed+1, k)&
&               , mfx(is:ie+1, js:je, k), mfx_ad(is:ie+1, js:je, k), mfy&
&               (is:ie, js:je+1, k), mfy_ad(is:ie, js:je+1, k), cx(is:ie&
&               +1, jsd:jed, k), cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied&
&               , js:je+1, k), cy_ad(isd:ied, js:je+1, k), crx(is:ie+1, &
&               jsd:jed, k), crx_ad(is:ie+1, jsd:jed, k), cry(isd:ied, &
&               js:je+1, k), cry_ad(isd:ied, js:je+1, k), xfx(is:ie+1, &
&               jsd:jed, k), xfx_ad(is:ie+1, jsd:jed, k), yfx(isd:ied, &
&               js:je+1, k), yfx_ad(isd:ied, js:je+1, k), q_con(isd:ied&
&               , jsd:jed, 1), z_rat(isd:ied, jsd:jed), z_rat_ad(isd:ied&
&               , jsd:jed), kgb, heat_s, heat_s_ad, dpx, dpx_ad, zvir, &
&               sphum, nq, q, q_ad, k, npz, flagstruct%inline_q, dt, &
&               flagstruct%hord_tr, hord_m, hord_v, hord_t, hord_p, &
&               nord_k, nord_v(k), nord_w, nord_t, flagstruct%dddmp, &
&               d2_divg, flagstruct%d4_bg, damp_vt(k), damp_w, damp_t, &
&               d_con_k, hydrostatic, gridstruct, flagstruct, bd, &
&               flagstructp%hord_tr_pert, hord_m_pert, hord_v_pert, &
&               hord_t_pert, hord_p_pert, flagstructp%split_damp, &
&               nord_k_pert, nord_v_pert(k), nord_w_pert, nord_t_pert, &
&               flagstructp%dddmp_pert, d2_divg_pert, flagstructp%&
&               d4_bg_pert, damp_vt_pert(k), damp_w_pert, damp_t_pert)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=jed,jsd,-1
            DO i=ied,isd,-1
              zh_ad(i, j, k) = zh_ad(i, j, k) + z_rat_ad(i, j)/radius
              zh_ad(i, j, k+1) = zh_ad(i, j, k+1) + z_rat_ad(i, j)/&
&               radius
              z_rat_ad(i, j) = 0.0
            END DO
          END DO
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) CALL A2B_ORD2_BWD(delp(isd:ied, jsd:jed, k), &
&                                      delp_ad(isd:ied, jsd:jed, k), wk&
&                                      , wk_ad, gridstruct, npx, npy, is&
&                                      , ie, js, je, ng, .false.)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(omga(i, j, k))
              delp_ad(i, j, k) = delp_ad(i, j, k) + omga_ad(i, j, k)
              omga_ad(i, j, k) = 0.0
            END DO
          END DO
        END IF
        CALL POPINTEGER(nord_v_pert(npz+1))
        CALL POPINTEGER(nord_v(npz+1))
        CALL POPREALARRAY(damp_vt_pert(npz+1))
        CALL POPREALARRAY(damp_vt(npz+1))
        CALL POPCONTROL(2,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(damp_vt_pert(k))
          CALL POPINTEGER(nord_v_pert(k))
        ELSE IF (branch .NE. 1) THEN
          GOTO 110
        END IF
        CALL POPCONTROL(3,branch)
        CALL POPCONTROL(1,branch)
 110    CALL POPREALARRAY(damp_t_pert)
        CALL POPINTEGER(nord_t_pert)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(damp_vt_pert(k))
        ELSE
          CALL POPREALARRAY(damp_vt_pert(k))
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(d2_divg_pert)
        ELSE
          CALL POPREALARRAY(d2_divg_pert)
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPINTEGER(nord_v_pert(k))
        ELSE
          CALL POPINTEGER(nord_v_pert(k))
        END IF
        CALL POPINTEGER(nord_k_pert)
        CALL POPINTEGER(hord_p_pert)
        CALL POPINTEGER(hord_v_pert)
        CALL POPINTEGER(hord_t_pert)
        CALL POPINTEGER(hord_m_pert)
        CALL POPCONTROL(3,branch)
        IF (branch .LT. 3) THEN
          IF (branch .NE. 0) THEN
            IF (branch .NE. 1) THEN
              CALL POPREALARRAY(damp_vt(k))
              CALL POPINTEGER(nord_v(k))
            END IF
          END IF
        ELSE
          IF (branch .LT. 5) THEN
            IF (branch .EQ. 3) THEN
              GOTO 120
            ELSE
              CALL POPREALARRAY(damp_vt(k))
              CALL POPINTEGER(nord_v(k))
            END IF
          ELSE IF (branch .NE. 5) THEN
            GOTO 120
          END IF
          CALL POPCONTROL(1,branch)
        END IF
 120    CALL POPINTEGER(nord_t)
        CALL POPINTEGER(nord_w)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(damp_vt(k))
        ELSE
          CALL POPREALARRAY(damp_vt(k))
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(d2_divg)
        ELSE
          CALL POPREALARRAY(d2_divg)
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPINTEGER(nord_v(k))
        ELSE
          CALL POPINTEGER(nord_v(k))
        END IF
        CALL POPINTEGER(nord_k)
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO iq=nq,1,-1
          CALL POPREALARRAY(q(isd:ied, jsd:jed, :, iq), (ied-isd+1)*(&
&                      jed-jsd+1)*npz)
          CALL NESTED_GRID_BC_APPLY_INTT_ADM(q(isd:ied, jsd:jed, :, iq)&
&                                      , q_ad(isd:ied, jsd:jed, :, iq), &
&                                      0, 0, npx, npy, npz, bd, arg1, &
&                                      arg2, neststruct%q_bc(iq), &
&                                      neststruct%nestbctype)
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL NESTED_GRID_BC_APPLY_INTT_ADM(divgd, divgd_ad, 1, 1, npx, &
&                                    npy, npz, bd, split_timestep_bc, &
&                                    arg1, neststruct%divg_bc, &
&                                    neststruct%nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT_ADM(uc, uc_ad, 1, 0, npx, npy, &
&                                    npz, bd, arg1, arg2, neststruct%&
&                                    uc_bc, neststruct%nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT_ADM(vc, vc_ad, 0, 1, npx, npy, &
&                                    npz, bd, arg1, arg2, neststruct%&
&                                    vc_bc, neststruct%nestbctype)
      END IF
      CALL START_GROUP_HALO_UPDATE_ADM(i_pack(9), uc, uc_ad, vc, vc_ad, &
&                                domain, gridtype=cgrid_ne)
      CALL P_GRAD_C_BWD(dt2, npz, delpc, delpc_ad, pkc, pkc_ad, gz, &
&                 gz_ad, uc, uc_ad, vc, vc_ad, bd, gridstruct%rdxc, &
&                 gridstruct%rdyc, hydrostatic)
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        CALL GEOPK_BWD(ptop, pe, pe_ad, peln, peln_ad, delpc, delpc_ad, &
&                pkc, pkc_ad, gz, gz_ad, phis, ptc, ptc_ad, q_con, pkz, &
&                pkz_ad, npz, akap, .true., gridstruct%nested, .false., &
&                npx, npy, flagstruct%a2b_ord, bd)
      ELSE
        IF (branch .EQ. 1) THEN
          CALL NEST_HALO_NH_BWD(ptop, grav, akap, cp, delpc, delpc_ad, &
&                         delz, delz_ad, ptc, ptc_ad, phis, pkc, pkc_ad&
&                         , gz, gz_ad, pk3, pk3_ad, npx, npy, npz, &
&                         gridstruct%nested, .false., .false., .false., &
&                         bd)
          CALL POPREALARRAY(delz, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*&
&                      npz)
          CALL NESTED_GRID_BC_APPLY_INTT_ADM(delz, delz_ad, 0, 0, npx, &
&                                      npy, npz, bd, arg1, arg2, &
&                                      neststruct%delz_bc, neststruct%&
&                                      nestbctype)
        END IF
        CALL RIEM_SOLVER_C_BWD(ms, dt2, is, ie, js, je, npz, ng, akap, &
&                        cappa, cp, ptop, phis, omga, omga_ad, ptc, &
&                        ptc_ad, q_con, delpc, delpc_ad, gz, gz_ad, pkc&
&                        , pkc_ad, ws3, ws3_ad, flagstruct%p_fac, &
&                        flagstruct%a_imp, flagstruct%scale_z)
        CALL UPDATE_DZ_C_BWD(is, ie, js, je, npz, ng, dt2, dp_ref, zs, &
&                      gridstruct%area, ut, ut_ad, vt, vt_ad, gz, gz_ad&
&                      , ws3, ws3_ad, npx, npy, gridstruct%sw_corner, &
&                      gridstruct%se_corner, gridstruct%ne_corner, &
&                      gridstruct%nw_corner, bd, gridstruct%grid_type)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO k=npz+1,1,-1
            DO j=jed,jsd,-1
              DO i=ied,isd,-1
                gz_ad(i, j, k) = gz_ad(i, j, k) + zh_ad(i, j, k)
                zh_ad(i, j, k) = 0.0
              END DO
            END DO
          END DO
        ELSE
          DO k=npz+1,1,-1
            DO j=jed,jsd,-1
              DO i=ied,isd,-1
                CALL POPREALARRAY(gz(i, j, k))
                zh_ad(i, j, k) = zh_ad(i, j, k) + gz_ad(i, j, k)
                gz_ad(i, j, k) = 0.0
              END DO
            END DO
          END DO
        END IF
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(ptc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
        CALL NESTED_GRID_BC_APPLY_INTT_ADM(ptc, ptc_ad, 0, 0, npx, npy, &
&                                    npz, bd, arg1, arg2, neststruct%&
&                                    pt_bc, neststruct%nestbctype)
        CALL POPREALARRAY(delpc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*&
&                    npz)
        CALL NESTED_GRID_BC_APPLY_INTT_ADM(delpc, delpc_ad, 0, 0, npx, &
&                                    npy, npz, bd, arg1, arg2, &
&                                    neststruct%delp_bc, neststruct%&
&                                    nestbctype)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) CALL START_GROUP_HALO_UPDATE_ADM(i_pack(3), &
&                                                   divgd, divgd_ad, &
&                                                   domain, position=&
&                                                   corner)
      DO k=npz,1,-1
        CALL C_SW_BWD(delpc(isd:ied, jsd:jed, k), delpc_ad(isd:ied, jsd:&
&               jed, k), delp(isd:ied, jsd:jed, k), delp_ad(isd:ied, jsd&
&               :jed, k), ptc(isd:ied, jsd:jed, k), ptc_ad(isd:ied, jsd:&
&               jed, k), pt(isd:ied, jsd:jed, k), pt_ad(isd:ied, jsd:jed&
&               , k), u(isd:ied, jsd:jed+1, k), u_ad(isd:ied, jsd:jed+1&
&               , k), v(isd:ied+1, jsd:jed, k), v_ad(isd:ied+1, jsd:jed&
&               , k), w(isd:ied, jsd:jed, k), w_ad(isd:ied, jsd:jed, k)&
&               , uc(isd:ied+1, jsd:jed, k), uc_ad(isd:ied+1, jsd:jed, k&
&               ), vc(isd:ied, jsd:jed+1, k), vc_ad(isd:ied, jsd:jed+1, &
&               k), ua(isd:ied, jsd:jed, k), ua_ad(isd:ied, jsd:jed, k)&
&               , va(isd:ied, jsd:jed, k), va_ad(isd:ied, jsd:jed, k), &
&               omga(isd:ied, jsd:jed, k), omga_ad(isd:ied, jsd:jed, k)&
&               , ut(isd:ied, jsd:jed, k), ut_ad(isd:ied, jsd:jed, k), &
&               vt(isd:ied, jsd:jed, k), vt_ad(isd:ied, jsd:jed, k), &
&               divgd(isd:ied+1, jsd:jed+1, k), divgd_ad(isd:ied+1, jsd:&
&               jed+1, k), flagstruct%nord, dt2, hydrostatic, .true., bd&
&               , gridstruct, flagstruct)
      END DO
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        DO j=je+1,js-1,-1
          DO k=npz,1,-1
            DO i=ie+1,is-1,-1
              pem_ad(i, k, j) = pem_ad(i, k, j) + pem_ad(i, k+1, j)
              delp_ad(i, j, k) = delp_ad(i, j, k) + pem_ad(i, k+1, j)
              pem_ad(i, k+1, j) = 0.0
            END DO
          END DO
          DO i=ie+1,is-1,-1
            pem_ad(i, 1, j) = 0.0
          END DO
        END DO
        pem_ad = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(beta_d)
      ELSE
        CALL POPREALARRAY(beta_d)
      END IF
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(gz, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*(npz+&
&                    1))
        CALL START_GROUP_HALO_UPDATE_ADM(i_pack(5), gz, gz_ad, domain)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=jed,jsd,-1
            DO k=1,npz,1
              DO i=ied,isd,-1
                CALL POPREALARRAY(gz(i, j, k))
                gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + gz_ad(i, j, k)
                delz_ad(i, j, k) = delz_ad(i, j, k) - gz_ad(i, j, k)
                gz_ad(i, j, k) = 0.0
              END DO
            END DO
            DO i=ied,isd,-1
              CALL POPREALARRAY(gz(i, j, npz+1))
              gz_ad(i, j, npz+1) = 0.0
            END DO
          END DO
        ELSE
          DO j=je,js,-1
            DO k=1,npz,1
              DO i=ie,is,-1
                CALL POPREALARRAY(gz(i, j, k))
                gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + gz_ad(i, j, k)
                delz_ad(i, j, k) = delz_ad(i, j, k) - gz_ad(i, j, k)
                gz_ad(i, j, k) = 0.0
              END DO
            END DO
            DO i=ie,is,-1
              CALL POPREALARRAY(gz(i, j, npz+1))
              gz_ad(i, j, npz+1) = 0.0
            END DO
          END DO
        END IF
      ELSE IF (branch .NE. 1) THEN
        GOTO 130
      END IF
      CALL POPREALARRAY(w, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz)
      CALL START_GROUP_HALO_UPDATE_ADM(i_pack(7), w, w_ad, domain)
 130  CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz*nq&
&                   )
        CALL START_GROUP_HALO_UPDATE_ADM(i_pack(10), q, q_ad, domain)
      END IF
      CALL POPCONTROL(1,branch)
    END DO
    CALL POPCONTROL(2,branch)
    CALL POPREALARRAY(cy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    CALL POPREALARRAY(cx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
    CALL POPREALARRAY(mfy, (bd%ie-bd%is+1)*(bd%je-bd%js+2)*npz)
    CALL POPREALARRAY(mfx, (bd%ie-bd%is+2)*(bd%je-bd%js+1)*npz)
  END SUBROUTINE DYN_CORE_BWD
!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
  SUBROUTINE DYN_CORE(npx, npy, npz, ng, sphum, nq, bdt, n_split, zvir, &
&   cp, akap, cappa, grav, hydrostatic, u, v, w, delz, pt, q, delp, pe, &
&   pk, phis, ws, omga, ptop, pfull, ua, va, uc, vc, mfx, mfy, cx, cy, &
&   pkz, peln, q_con, ak, bk, dpx, ks, gridstruct, flagstruct, &
&   flagstructp, neststruct, idiag, bd, domain, init_step, i_pack, &
&   end_step, gz, pkc, ptc, crx, xfx, cry, yfx, divgd, delpc, ut, vt, zh&
&   , pk3, du, dv, time_total)
    IMPLICIT NONE
! end init_step
! Start of the big dynamic time stepping
!allocate(    gz(isd:ied, jsd:jed ,npz+1) )
!  call init_ijk_mem(isd,ied, jsd,jed, npz+1, gz, huge_r)
!allocate(   pkc(isd:ied, jsd:jed ,npz+1) )
!allocate(   ptc(isd:ied, jsd:jed ,npz ) )
!allocate( crx(is :ie+1, jsd:jed,  npz) )
!allocate( xfx(is :ie+1, jsd:jed,  npz) )
!allocate( cry(isd:ied,  js :je+1, npz) )
!allocate( yfx(isd:ied,  js :je+1, npz) )
!allocate( divgd(isd:ied+1,jsd:jed+1,npz) )
!allocate( delpc(isd:ied, jsd:jed  ,npz  ) )
!                    call init_ijk_mem(isd,ied, jsd,jed, npz, delpc, 0.)
!allocate( ut(isd:ied, jsd:jed, npz) )
!                    call init_ijk_mem(isd,ied, jsd,jed, npz, ut, 0.)
!allocate( vt(isd:ied, jsd:jed, npz) )
!                    call init_ijk_mem(isd,ied, jsd,jed, npz, vt, 0.)
!allocate( zh(isd:ied, jsd:jed, npz+1) )
!              call init_ijk_mem(isd,ied, jsd,jed, npz+1, zh, huge_r )
!allocate ( pk3(isd:ied,jsd:jed,npz+1) )
!call init_ijk_mem(isd,ied, jsd,jed, npz+1, pk3, huge_r )
!if (allocated(heat_source)) deallocate( heat_source ) !If ncon == 0 but d_con > 1.e-5, this would not be deallocated in earlier 
!versions of the code
!deallocate(    gz )
!deallocate(   ptc )
!deallocate(   crx )
!deallocate(   xfx )
!deallocate(   cry )
!deallocate(   yfx )
!deallocate( divgd )
!deallocate(   pkc )
!deallocate( delpc )
!if( allocated(ut))   deallocate( ut )
!if( allocated(vt))   deallocate( vt )
!if ( allocated (du) ) deallocate( du )
!if ( allocated (dv) ) deallocate( dv )
!if ( .not. hydrostatic ) then
!     deallocate( zh )
!     if( allocated(pk3) )   deallocate ( pk3 )
!endif
!if( allocated(pem) )   deallocate ( pem )
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: ng, nq, sphum
    INTEGER, INTENT(IN) :: n_split
    REAL, INTENT(IN) :: bdt
    REAL, INTENT(IN) :: zvir, cp, akap, grav
    REAL, INTENT(IN) :: ptop
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: init_step, end_step
    REAL, INTENT(IN) :: pfull(npz)
    REAL, DIMENSION(npz+1), INTENT(IN) :: ak, bk
    INTEGER, INTENT(IN) :: ks
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: i_pack(*)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
! vertical vel. (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! delta-height (m, negative)
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! moist kappa
    REAL, INTENT(INOUT) :: cappa(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! total time (seconds) since start
    REAL, INTENT(IN), OPTIONAL :: time_total
!-----------------------------------------------------------------------
! Auxilliary pressure arrays:
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
! edge pressure (pascal)
    REAL, INTENT(INOUT) :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
! ln(pe)
    REAL, INTENT(INOUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
! pe**kappa
    REAL, INTENT(INOUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL(kind=8), INTENT(INOUT) :: dpx(bd%is:bd%ie, bd%js:bd%je)
!-----------------------------------------------------------------------
! Others:
    REAL, PARAMETER :: near0=1.e-8
    REAL, PARAMETER :: huge_r=1.e8
!-----------------------------------------------------------------------
! w at surface
    REAL, INTENT(OUT) :: ws(bd%is:bd%ie, bd%js:bd%je)
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! (uc, vc) are mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua, va
    REAL, INTENT(INOUT) :: q_con(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! The Flux capacitors: accumulated Mass flux arrays
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je, npz), INTENT(INOUT) :: pkz
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    TYPE(FV_FLAGS_PERT_TYPE), INTENT(IN), TARGET :: flagstructp
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(FV_DIAG_TYPE), INTENT(IN) :: idiag
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
!real, allocatable, dimension(:,:,:):: pem, heat_source
    REAL :: pem(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1), heat_source(bd&
&   %isd:bd%ied, bd%jsd:bd%jed, npz)
! Auto 1D & 2D arrays:
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ws3, z_rat
    REAL :: dp_ref(npz)
! surface height (m)
    REAL :: zs(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: p1d(bd%is:bd%ie)
    REAL :: om2d(bd%is:bd%ie, npz)
    REAL :: wbuffer(npy+2, npz)
    REAL :: ebuffer(npy+2, npz)
    REAL :: nbuffer(npx+2, npz)
    REAL :: sbuffer(npx+2, npz)
! ----   For external mode:
    REAL :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: fz(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: heat_s(bd%is:bd%ie, bd%js:bd%je)
    REAL :: damp_vt(npz+1)
    INTEGER :: nord_v(npz+1)
!-------------------------------------
    INTEGER :: hord_m, hord_v, hord_t, hord_p
    INTEGER :: nord_k, nord_w, nord_t
    INTEGER :: ms
!---------------------------------------
    INTEGER :: hord_m_pert, hord_v_pert, hord_t_pert, hord_p_pert
    INTEGER :: nord_k_pert, nord_w_pert, nord_t_pert, nord_v_pert(npz+1)
    REAL :: d2_divg_pert, damp_vt_pert(npz+1), damp_w_pert, damp_t_pert
!---------------------------------------
    INTEGER :: i, j, k, it, iq, n_con, nf_ke
    INTEGER :: iep1, jep1
    REAL :: beta, beta_d, d_con_k, damp_w, damp_t, kgb, cv_air
    REAL :: dt, dt2, rdt
    REAL :: d2_divg
    REAL :: k1k, rdg, dtmp, delt
    LOGICAL :: last_step, remap_step
    LOGICAL :: used
    REAL :: split_timestep_bc
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pkc(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: ptc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cry(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: divgd(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: delpc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ut(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: zh(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk3(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    INTRINSIC LOG
    INTRINSIC REAL
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC EXP
    INTRINSIC ABS
    INTRINSIC SIGN
    INTEGER :: max1
    INTEGER :: max2
    REAL :: min1
    REAL :: min2
    REAL :: abs0
    REAL :: arg1
    REAL :: arg2
    LOGICAL :: arg10
    REAL*8 :: arg11
    REAL :: x1
    REAL :: y2
    REAL :: y1
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    peln1 = LOG(ptop)
    ptk = ptop**akap
    dt = bdt/REAL(n_split)
    dt2 = 0.5*dt
    rdt = 1.0/dt
    IF (1 .LT. flagstruct%m_split/2) THEN
      ms = flagstruct%m_split/2
    ELSE
      ms = 1
    END IF
    beta = flagstruct%beta
    rdg = -(rdgas/grav)
    cv_air = cp_air - rdgas
! Indexes:
    iep1 = ie + 1
    jep1 = je + 1
    IF (.NOT.hydrostatic) THEN
      rgrav = 1.0/grav
! rg/Cv=0.4
      k1k = akap/(1.-akap)
!$OMP parallel do default(none) shared(npz,dp_ref,ak,bk)
      DO k=1,npz
        dp_ref(k) = ak(k+1) - ak(k) + (bk(k+1)-bk(k))*1.e5
      END DO
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,zs,phis,rgrav)
      DO j=jsd,jed
        DO i=isd,ied
          zs(i, j) = phis(i, j)*rgrav
        END DO
      END DO
    END IF
!allocate( du(isd:ied,  jsd:jed+1,npz) )
!call init_ijk_mem(isd,ied,   jsd,jed+1, npz, du, 0.)
!allocate( dv(isd:ied+1,jsd:jed,  npz) )
!call init_ijk_mem(isd,ied+1, jsd,jed  , npz, dv, 0.)
! Empty the "flux capacitors"
!call init_ijk_mem(is, ie+1, js,  je,   npz, mfx, 0.)
    mfx = 0.0
!call init_ijk_mem(is, ie  , js,  je+1, npz, mfy, 0.)
    mfy = 0.0
!call init_ijk_mem(is, ie+1, jsd, jed,  npz, cx, 0.)
    cx = 0.0
!call init_ijk_mem(isd, ied, js,  je+1, npz, cy, 0.)
    cy = 0.0
    IF (flagstruct%d_con .GT. 1.0e-5) heat_source = 0.0
!allocate( heat_source(isd:ied, jsd:jed, npz) )
!call init_ijk_mem(isd, ied, jsd, jed, npz, heat_source, 0.)
    IF (flagstruct%convert_ke .OR. flagstruct%vtdm4 .GT. 1.e-4) THEN
      n_con = npz
    ELSE IF (flagstruct%d2_bg_k1 .LT. 1.e-3) THEN
      n_con = 0
    ELSE IF (flagstruct%d2_bg_k2 .LT. 1.e-3) THEN
      n_con = 1
    ELSE
      n_con = 2
    END IF
!-----------------------------------------------------
    DO it=1,n_split
!-----------------------------------------------------
      IF (flagstruct%breed_vortex_inline .OR. it .EQ. n_split) THEN
        remap_step = .true.
      ELSE
        remap_step = .false.
      END IF
      IF (flagstruct%fv_debug) THEN
        IF (IS_MASTER()) WRITE(*, *) 'n_split loop, it=', it
        IF (.NOT.flagstruct%hydrostatic) CALL PRT_MXM('delz', delz, is, &
&                                               ie, js, je, ng, npz, 1.&
&                                               , gridstruct%area_64, &
&                                               domain)
        CALL PRT_MXM('PT', pt, is, ie, js, je, ng, npz, 1., gridstruct%&
&              area_64, domain)
      END IF
      IF (gridstruct%nested) split_timestep_bc = REAL(n_split*flagstruct&
&         %k_split + neststruct%nest_timestep)
!First split timestep has split_timestep_BC = n_split*k_split
!   to do time-extrapolation on BCs.
      IF (nq .GT. 0) THEN
        CALL TIMING_ON('COMM_TOTAL')
        CALL TIMING_ON('COMM_TRACER')
        IF (flagstruct%inline_q) CALL START_GROUP_HALO_UPDATE(i_pack(10)&
&                                                       , q, domain)
        CALL TIMING_OFF('COMM_TRACER')
        CALL TIMING_OFF('COMM_TOTAL')
      END IF
      IF (.NOT.hydrostatic) THEN
        CALL TIMING_ON('COMM_TOTAL')
        CALL START_GROUP_HALO_UPDATE(i_pack(7), w, domain)
        CALL TIMING_OFF('COMM_TOTAL')
        IF (it .EQ. 1) THEN
          IF (gridstruct%nested) THEN
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,gz,zs,delz)
            DO j=jsd,jed
              DO i=isd,ied
                gz(i, j, npz+1) = zs(i, j)
              END DO
              DO k=npz,1,-1
                DO i=isd,ied
                  gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)
                END DO
              END DO
            END DO
          ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gz,zs,delz)
            DO j=js,je
              DO i=is,ie
                gz(i, j, npz+1) = zs(i, j)
              END DO
              DO k=npz,1,-1
                DO i=is,ie
                  gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)
                END DO
              END DO
            END DO
          END IF
          CALL TIMING_ON('COMM_TOTAL')
          CALL START_GROUP_HALO_UPDATE(i_pack(5), gz, domain)
          CALL TIMING_OFF('COMM_TOTAL')
        END IF
      END IF
      IF (it .EQ. 1) THEN
        CALL TIMING_ON('COMM_TOTAL')
        CALL COMPLETE_GROUP_HALO_UPDATE(i_pack(1), domain)
        CALL TIMING_OFF('COMM_TOTAL')
        beta_d = 0.
      ELSE
        beta_d = beta
      END IF
      IF (it .EQ. n_split .AND. end_step) THEN
        IF (flagstruct%use_old_omega) THEN
          pem = 0.0
!allocate ( pem(is-1:ie+1,npz+1,js-1:je+1) )
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pem,delp,ptop)
          DO j=js-1,je+1
            DO i=is-1,ie+1
              pem(i, 1, j) = ptop
            END DO
            DO k=1,npz
              DO i=is-1,ie+1
                pem(i, k+1, j) = pem(i, k, j) + delp(i, j, k)
              END DO
            END DO
          END DO
        END IF
        last_step = .true.
      ELSE
        last_step = .false.
      END IF
      CALL TIMING_ON('COMM_TOTAL')
      CALL COMPLETE_GROUP_HALO_UPDATE(i_pack(8), domain)
      IF (.NOT.hydrostatic) CALL COMPLETE_GROUP_HALO_UPDATE(i_pack(7), &
&                                                     domain)
      CALL TIMING_OFF('COMM_TOTAL')
      CALL TIMING_ON('c_sw')
!$OMP parallel do default(none) shared(npz,isd,jsd,delpc,delp,ptc,pt,u,v,w,uc,vc,ua,va, &
!$OMP                                  omga,ut,vt,divgd,flagstruct,dt2,hydrostatic,bd,  &
!$OMP                                  gridstruct)
      DO k=1,npz
        CALL C_SW(delpc(isd:ied, jsd:jed, k), delp(isd:ied, jsd:jed, k)&
&           , ptc(isd:ied, jsd:jed, k), pt(isd:ied, jsd:jed, k), u(isd:&
&           ied, jsd:jed+1, k), v(isd:ied+1, jsd:jed, k), w(isd:ied, jsd&
&           :jed, k), uc(isd:ied+1, jsd:jed, k), vc(isd:ied, jsd:jed+1, &
&           k), ua(isd:ied, jsd:jed, k), va(isd:ied, jsd:jed, k), omga(&
&           isd:ied, jsd:jed, k), ut(isd:ied, jsd:jed, k), vt(isd:ied, &
&           jsd:jed, k), divgd(isd:ied+1, jsd:jed+1, k), flagstruct%nord&
&           , dt2, hydrostatic, .true., bd, gridstruct, flagstruct)
      END DO
      CALL TIMING_OFF('c_sw')
      IF (flagstruct%nord .GT. 0) THEN
        CALL TIMING_ON('COMM_TOTAL')
        CALL START_GROUP_HALO_UPDATE(i_pack(3), divgd, domain, position=&
&                              corner)
        CALL TIMING_OFF('COMM_TOTAL')
      END IF
      IF (gridstruct%nested) THEN
        arg1 = split_timestep_bc + 0.5
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(delpc, 0, 0, npx, npy, npz, bd, &
&                                arg1, arg2, neststruct%delp_bc, &
&                                neststruct%nestbctype)
        arg1 = split_timestep_bc + 0.5
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(ptc, 0, 0, npx, npy, npz, bd, &
&                                arg1, arg2, neststruct%pt_bc, &
&                                neststruct%nestbctype)
      END IF
! end hydro check
      IF (hydrostatic) THEN
        CALL GEOPK(ptop, pe, peln, delpc, pkc, gz, phis, ptc, q_con, pkz&
&            , npz, akap, .true., gridstruct%nested, .false., npx, npy, &
&            flagstruct%a2b_ord, bd)
      ELSE
        IF (it .EQ. 1) THEN
          CALL TIMING_ON('COMM_TOTAL')
          CALL COMPLETE_GROUP_HALO_UPDATE(i_pack(5), domain)
          CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,zh,gz)
          DO k=1,npz+1
            DO j=jsd,jed
              DO i=isd,ied
! Save edge heights for update_dz_d
                zh(i, j, k) = gz(i, j, k)
              END DO
            END DO
          END DO
        ELSE
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,zh,gz)
          DO k=1,npz+1
            DO j=jsd,jed
              DO i=isd,ied
                gz(i, j, k) = zh(i, j, k)
              END DO
            END DO
          END DO
        END IF
        CALL TIMING_ON('UPDATE_DZ_C')
        CALL UPDATE_DZ_C(is, ie, js, je, npz, ng, dt2, dp_ref, zs, &
&                  gridstruct%area, ut, vt, gz, ws3, npx, npy, &
&                  gridstruct%sw_corner, gridstruct%se_corner, &
&                  gridstruct%ne_corner, gridstruct%nw_corner, bd, &
&                  gridstruct%grid_type)
        CALL TIMING_OFF('UPDATE_DZ_C')
        CALL TIMING_ON('Riem_Solver')
        CALL RIEM_SOLVER_C(ms, dt2, is, ie, js, je, npz, ng, akap, cappa&
&                    , cp, ptop, phis, omga, ptc, q_con, delpc, gz, pkc&
&                    , ws3, flagstruct%p_fac, flagstruct%a_imp, &
&                    flagstruct%scale_z)
        CALL TIMING_OFF('Riem_Solver')
        IF (gridstruct%nested) THEN
          arg1 = split_timestep_bc + 0.5
          arg2 = REAL(n_split*flagstruct%k_split)
          CALL NESTED_GRID_BC_APPLY_INTT(delz, 0, 0, npx, npy, npz, bd, &
&                                  arg1, arg2, neststruct%delz_bc, &
&                                  neststruct%nestbctype)
!Compute gz/pkc
!NOTE: nominally only need to compute quantities one out in the halo for p_grad_c
!(instead of entire halo)
          CALL NEST_HALO_NH(ptop, grav, akap, cp, delpc, delz, ptc, phis&
&                     , pkc, gz, pk3, npx, npy, npz, gridstruct%nested, &
&                     .false., .false., .false., bd)
        END IF
      END IF
      CALL P_GRAD_C(dt2, npz, delpc, pkc, gz, uc, vc, bd, gridstruct%&
&             rdxc, gridstruct%rdyc, hydrostatic)
      CALL TIMING_ON('COMM_TOTAL')
      CALL START_GROUP_HALO_UPDATE(i_pack(9), uc, vc, domain, gridtype=&
&                            cgrid_ne)
      CALL TIMING_OFF('COMM_TOTAL')
      CALL TIMING_ON('COMM_TOTAL')
      IF (flagstruct%inline_q .AND. nq .GT. 0) CALL &
&       COMPLETE_GROUP_HALO_UPDATE(i_pack(10), domain)
      IF (flagstruct%nord .GT. 0) CALL COMPLETE_GROUP_HALO_UPDATE(i_pack&
&                                                           (3), domain)
      CALL COMPLETE_GROUP_HALO_UPDATE(i_pack(9), domain)
      CALL TIMING_OFF('COMM_TOTAL')
      IF (gridstruct%nested) THEN
!On a nested grid we have to do SOMETHING with uc and vc in
! the boundary halo, particularly at the corners of the
! domain and of each processor element. We must either
! apply an interpolated BC, or extrapolate into the
! boundary halo
! NOTE:
!The update_domains calls for uc and vc need to go BEFORE the BCs to ensure cross-restart
!bitwise-consistent solutions when doing the spatial extrapolation; should not make a
!difference for interpolated BCs from the coarse grid.
        arg1 = split_timestep_bc + 0.5
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(vc, 0, 1, npx, npy, npz, bd, arg1&
&                                , arg2, neststruct%vc_bc, neststruct%&
&                                nestbctype)
        arg1 = split_timestep_bc + 0.5
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(uc, 1, 0, npx, npy, npz, bd, arg1&
&                                , arg2, neststruct%uc_bc, neststruct%&
&                                nestbctype)
!QUESTION: What to do with divgd in nested halo?
        arg1 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(divgd, 1, 1, npx, npy, npz, bd, &
&                                split_timestep_bc, arg1, neststruct%&
&                                divg_bc, neststruct%nestbctype)
!!$            if (is == 1 .and. js == 1) then
!!$               do j=jsd,5
!!$                  write(mpp_pe()+2000,*) j, divg(isd:5,j,1)
!!$            endif
      END IF
      IF (gridstruct%nested .AND. flagstruct%inline_q) THEN
        DO iq=1,nq
          arg1 = split_timestep_bc + 1
          arg2 = REAL(n_split*flagstruct%k_split)
          CALL NESTED_GRID_BC_APPLY_INTT(q(isd:ied, jsd:jed, :, iq), 0, &
&                                  0, npx, npy, npz, bd, arg1, arg2, &
&                                  neststruct%q_bc(iq), neststruct%&
&                                  nestbctype)
        END DO
      END IF
      CALL TIMING_ON('d_sw')
!$OMP parallel do default(none) shared(npz,flagstruct,nord_v,pfull,damp_vt,hydrostatic,last_step, &
!$OMP                                  is,ie,js,je,isd,ied,jsd,jed,omga,delp,gridstruct,npx,npy,  &
!$OMP                                  ng,zh,vt,ptc,pt,u,v,w,uc,vc,ua,va,divgd,mfx,mfy,cx,cy,     &
!$OMP                                  crx,cry,xfx,yfx,q_con,zvir,sphum,nq,q,dt,bd,rdt,iep1,jep1, &
!$OMP                                  heat_source)                                               &
!$OMP                          private(nord_k, nord_w, nord_t, damp_w, damp_t, d2_divg,   &
!$OMP                          d_con_k,kgb, hord_m, hord_v, hord_t, hord_p, wk, heat_s, z_rat)
      DO k=1,npz
        hord_m = flagstruct%hord_mt
        hord_t = flagstruct%hord_tm
        hord_v = flagstruct%hord_vt
        hord_p = flagstruct%hord_dp
        nord_k = flagstruct%nord
!      if ( k==npz ) then
        kgb = flagstruct%ke_bg
        IF (2 .GT. flagstruct%nord) THEN
          nord_v(k) = flagstruct%nord
        ELSE
          nord_v(k) = 2
        END IF
        IF (0.20 .GT. flagstruct%d2_bg) THEN
          d2_divg = flagstruct%d2_bg
        ELSE
          d2_divg = 0.20
        END IF
        IF (flagstruct%do_vort_damp) THEN
! for delp, delz, and vorticity
          damp_vt(k) = flagstruct%vtdm4
        ELSE
          damp_vt(k) = 0.
        END IF
        nord_w = nord_v(k)
        nord_t = nord_v(k)
        damp_w = damp_vt(k)
        damp_t = damp_vt(k)
        d_con_k = flagstruct%d_con
        IF (npz .EQ. 1 .OR. flagstruct%n_sponge .LT. 0) THEN
          d2_divg = flagstruct%d2_bg
        ELSE IF (k .EQ. 1) THEN
! Sponge layers with del-2 damping on divergence, vorticity, w, z, and air mass (delp).
! no special damping of potential temperature in sponge layers
! Divergence damping:
          nord_k = 0
          IF (0.01 .LT. flagstruct%d2_bg) THEN
            IF (flagstruct%d2_bg .LT. flagstruct%d2_bg_k1) THEN
              d2_divg = flagstruct%d2_bg_k1
            ELSE
              d2_divg = flagstruct%d2_bg
            END IF
          ELSE IF (0.01 .LT. flagstruct%d2_bg_k1) THEN
            d2_divg = flagstruct%d2_bg_k1
          ELSE
            d2_divg = 0.01
          END IF
! Vertical velocity:
          nord_w = 0
          damp_w = d2_divg
          IF (flagstruct%do_vort_damp) THEN
! damping on delp and vorticity:
            nord_v(k) = 0
            damp_vt(k) = 0.5*d2_divg
          END IF
          d_con_k = 0.
        ELSE
          IF (2 .LT. flagstruct%n_sponge - 1) THEN
            max1 = flagstruct%n_sponge - 1
          ELSE
            max1 = 2
          END IF
          IF (k .EQ. max1 .AND. flagstruct%d2_bg_k2 .GT. 0.01) THEN
            nord_k = 0
            IF (flagstruct%d2_bg .LT. flagstruct%d2_bg_k2) THEN
              d2_divg = flagstruct%d2_bg_k2
            ELSE
              d2_divg = flagstruct%d2_bg
            END IF
            nord_w = 0
            damp_w = d2_divg
            IF (flagstruct%do_vort_damp) THEN
              nord_v(k) = 0
              damp_vt(k) = 0.5*d2_divg
            END IF
            d_con_k = 0.
          ELSE
            IF (3 .LT. flagstruct%n_sponge) THEN
              max2 = flagstruct%n_sponge
            ELSE
              max2 = 3
            END IF
            IF (k .EQ. max2 .AND. flagstruct%d2_bg_k2 .GT. 0.05) THEN
              nord_k = 0
              IF (flagstruct%d2_bg .LT. 0.2*flagstruct%d2_bg_k2) THEN
                d2_divg = 0.2*flagstruct%d2_bg_k2
              ELSE
                d2_divg = flagstruct%d2_bg
              END IF
              nord_w = 0
              damp_w = d2_divg
              d_con_k = 0.
            END IF
          END IF
        END IF
        hord_m_pert = flagstructp%hord_mt_pert
        hord_t_pert = flagstructp%hord_tm_pert
        hord_v_pert = flagstructp%hord_vt_pert
        hord_p_pert = flagstructp%hord_dp_pert
        nord_k_pert = flagstructp%nord_pert
        IF (2 .GT. flagstructp%nord_pert) THEN
          nord_v_pert(k) = flagstructp%nord_pert
        ELSE
          nord_v_pert(k) = 2
        END IF
        IF (0.20 .GT. flagstructp%d2_bg_pert) THEN
          d2_divg_pert = flagstructp%d2_bg_pert
        ELSE
          d2_divg_pert = 0.20
        END IF
        IF (flagstructp%do_vort_damp_pert) THEN
! for delp, delz, and vorticity
          damp_vt_pert(k) = flagstructp%vtdm4_pert
        ELSE
          damp_vt_pert(k) = 0.
        END IF
        nord_w_pert = nord_v_pert(k)
        nord_t_pert = nord_v_pert(k)
        damp_w_pert = damp_vt_pert(k)
        damp_t_pert = damp_vt_pert(k)
!Sponge layers for the pertuabtiosn
        IF (k .LE. flagstructp%n_sponge_pert) THEN
          IF (k .LE. flagstructp%n_sponge_pert - 1) THEN
            IF (flagstructp%hord_ks_traj) THEN
              hord_m = flagstructp%hord_mt_ks_traj
              hord_t = flagstructp%hord_tm_ks_traj
              hord_v = flagstructp%hord_vt_ks_traj
              hord_p = flagstructp%hord_dp_ks_traj
            END IF
            IF (flagstructp%hord_ks_pert) THEN
              hord_m_pert = flagstructp%hord_mt_ks_pert
              hord_t_pert = flagstructp%hord_tm_ks_pert
              hord_v_pert = flagstructp%hord_vt_ks_pert
              hord_p_pert = flagstructp%hord_dp_ks_pert
            END IF
          END IF
          nord_k_pert = 0
          IF (k .EQ. 1) THEN
            IF (0.01 .LT. flagstructp%d2_bg_pert) THEN
              IF (flagstructp%d2_bg_pert .LT. flagstructp%d2_bg_k1_pert&
&             ) THEN
                d2_divg_pert = flagstructp%d2_bg_k1_pert
              ELSE
                d2_divg_pert = flagstructp%d2_bg_pert
              END IF
            ELSE IF (0.01 .LT. flagstructp%d2_bg_k1_pert) THEN
              d2_divg_pert = flagstructp%d2_bg_k1_pert
            ELSE
              d2_divg_pert = 0.01
            END IF
          ELSE IF (k .EQ. 2) THEN
            IF (0.01 .LT. flagstructp%d2_bg_pert) THEN
              IF (flagstructp%d2_bg_pert .LT. flagstructp%d2_bg_k2_pert&
&             ) THEN
                d2_divg_pert = flagstructp%d2_bg_k2_pert
              ELSE
                d2_divg_pert = flagstructp%d2_bg_pert
              END IF
            ELSE IF (0.01 .LT. flagstructp%d2_bg_k2_pert) THEN
              d2_divg_pert = flagstructp%d2_bg_k2_pert
            ELSE
              d2_divg_pert = 0.01
            END IF
          ELSE IF (0.01 .LT. flagstructp%d2_bg_pert) THEN
            IF (flagstructp%d2_bg_pert .LT. flagstructp%d2_bg_ks_pert) &
&           THEN
              d2_divg_pert = flagstructp%d2_bg_ks_pert
            ELSE
              d2_divg_pert = flagstructp%d2_bg_pert
            END IF
          ELSE IF (0.01 .LT. flagstructp%d2_bg_ks_pert) THEN
            d2_divg_pert = flagstructp%d2_bg_ks_pert
          ELSE
            d2_divg_pert = 0.01
          END IF
          nord_w_pert = 0
          damp_w_pert = d2_divg_pert
          IF (flagstructp%do_vort_damp_pert) THEN
            nord_v_pert(k) = 0
            damp_vt_pert(k) = 0.5*d2_divg_pert
          END IF
        END IF
!Tapenade issue if not defined at level npz+1
        damp_vt(npz+1) = damp_vt(npz)
        damp_vt_pert(npz+1) = damp_vt_pert(npz)
        nord_v(npz+1) = nord_v(npz)
        nord_v_pert(npz+1) = nord_v_pert(npz)
        IF (hydrostatic .AND. (.NOT.flagstruct%use_old_omega) .AND. &
&           last_step) THEN
! Average horizontal "convergence" to cell center
          DO j=js,je
            DO i=is,ie
              omga(i, j, k) = delp(i, j, k)
            END DO
          END DO
        END IF
!--- external mode divergence damping ---
        IF (flagstruct%d_ext .GT. 0.) CALL A2B_ORD2(delp(isd:ied, jsd:&
&                                             jed, k), wk, gridstruct, &
&                                             npx, npy, is, ie, js, je, &
&                                             ng, .false.)
        IF (.NOT.hydrostatic .AND. flagstruct%do_f3d) THEN
! Correction factor for 3D Coriolis force
          DO j=jsd,jed
            DO i=isd,ied
              z_rat(i, j) = 1. + (zh(i, j, k)+zh(i, j, k+1))/radius
            END DO
          END DO
        END IF
        CALL D_SW(vt(isd:ied, jsd:jed, k), delp(isd:ied, jsd:jed, k), &
&           ptc(isd:ied, jsd:jed, k), pt(isd:ied, jsd:jed, k), u(isd:ied&
&           , jsd:jed+1, k), v(isd:ied+1, jsd:jed, k), w(isd:ied, jsd:&
&           jed, k), uc(isd:ied+1, jsd:jed, k), vc(isd:ied, jsd:jed+1, k&
&           ), ua(isd:ied, jsd:jed, k), va(isd:ied, jsd:jed, k), divgd(&
&           isd:ied+1, jsd:jed+1, k), mfx(is:ie+1, js:je, k), mfy(is:ie&
&           , js:je+1, k), cx(is:ie+1, jsd:jed, k), cy(isd:ied, js:je+1&
&           , k), crx(is:ie+1, jsd:jed, k), cry(isd:ied, js:je+1, k), &
&           xfx(is:ie+1, jsd:jed, k), yfx(isd:ied, js:je+1, k), q_con(&
&           isd:ied, jsd:jed, 1), z_rat(isd:ied, jsd:jed), kgb, heat_s, &
&           dpx, zvir, sphum, nq, q, k, npz, flagstruct%inline_q, dt, &
&           flagstruct%hord_tr, hord_m, hord_v, hord_t, hord_p, nord_k, &
&           nord_v(k), nord_w, nord_t, flagstruct%dddmp, d2_divg, &
&           flagstruct%d4_bg, damp_vt(k), damp_w, damp_t, d_con_k, &
&           hydrostatic, gridstruct, flagstruct, bd, flagstructp%&
&           hord_tr_pert, hord_m_pert, hord_v_pert, hord_t_pert, &
&           hord_p_pert, flagstructp%split_damp, nord_k_pert, &
&           nord_v_pert(k), nord_w_pert, nord_t_pert, flagstructp%&
&           dddmp_pert, d2_divg_pert, flagstructp%d4_bg_pert, &
&           damp_vt_pert(k), damp_w_pert, damp_t_pert)
        IF (hydrostatic .AND. (.NOT.flagstruct%use_old_omega) .AND. &
&           last_step) THEN
! Average horizontal "convergence" to cell center
          DO j=js,je
            DO i=is,ie
              omga(i, j, k) = omga(i, j, k)*(xfx(i, j, k)-xfx(i+1, j, k)&
&               +yfx(i, j, k)-yfx(i, j+1, k))*gridstruct%rarea(i, j)*rdt
            END DO
          END DO
        END IF
        IF (flagstruct%d_ext .GT. 0.) THEN
          DO j=js,jep1
            DO i=is,iep1
! delp at cell corners
              ptc(i, j, k) = wk(i, j)
            END DO
          END DO
        END IF
        IF (flagstruct%d_con .GT. 1.0e-5) THEN
! Average horizontal "convergence" to cell center
          DO j=js,je
            DO i=is,ie
              heat_source(i, j, k) = heat_source(i, j, k) + heat_s(i, j)
            END DO
          END DO
        END IF
      END DO
! end openMP k-loop
      CALL TIMING_OFF('d_sw')
      IF (flagstruct%fill_dp) CALL MIX_DP(hydrostatic, w, delp, pt, npz&
&                                   , ak, bk, .false., flagstruct%&
&                                   fv_debug, bd)
      CALL TIMING_ON('COMM_TOTAL')
      CALL START_GROUP_HALO_UPDATE(i_pack(1), delp, domain, complete=&
&                            .true.)
      CALL START_GROUP_HALO_UPDATE(i_pack(1), pt, domain, complete=&
&                            .true.)
      CALL TIMING_OFF('COMM_TOTAL')
      IF (flagstruct%d_ext .GT. 0.) THEN
        d2_divg = flagstruct%d_ext*gridstruct%da_min_c
!$OMP parallel do default(none) shared(is,iep1,js,jep1,npz,wk,ptc,divg2,vt,d2_divg)
        DO j=js,jep1
          DO i=is,iep1
            wk(i, j) = ptc(i, j, 1)
            divg2(i, j) = wk(i, j)*vt(i, j, 1)
          END DO
          DO k=2,npz
            DO i=is,iep1
              wk(i, j) = wk(i, j) + ptc(i, j, k)
              divg2(i, j) = divg2(i, j) + ptc(i, j, k)*vt(i, j, k)
            END DO
          END DO
          DO i=is,iep1
            divg2(i, j) = d2_divg*divg2(i, j)/wk(i, j)
          END DO
        END DO
      ELSE
        divg2(:, :) = 0.
      END IF
      CALL TIMING_ON('COMM_TOTAL')
      CALL COMPLETE_GROUP_HALO_UPDATE(i_pack(1), domain)
      CALL TIMING_OFF('COMM_TOTAL')
      IF (flagstruct%fv_debug) THEN
        IF (.NOT.flagstruct%hydrostatic) CALL PRT_MXM('delz', delz, is, &
&                                               ie, js, je, ng, npz, 1.&
&                                               , gridstruct%area_64, &
&                                               domain)
      END IF
!Want to move this block into the hydro/nonhydro branch above and merge the two if structures
      IF (gridstruct%nested) THEN
        arg1 = split_timestep_bc + 1
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(delp, 0, 0, npx, npy, npz, bd, &
&                                arg1, arg2, neststruct%delp_bc, &
&                                neststruct%nestbctype)
        arg1 = split_timestep_bc + 1
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(pt, 0, 0, npx, npy, npz, bd, arg1&
&                                , arg2, neststruct%pt_bc, neststruct%&
&                                nestbctype)
      END IF
! end hydro check
      IF (hydrostatic) THEN
        CALL GEOPK(ptop, pe, peln, delp, pkc, gz, phis, pt, q_con, pkz, &
&            npz, akap, .false., gridstruct%nested, .true., npx, npy, &
&            flagstruct%a2b_ord, bd)
      ELSE
        CALL TIMING_ON('UPDATE_DZ')
        CALL UPDATE_DZ_D(nord_v, damp_vt, flagstruct%hord_tm, is, ie, js&
&                  , je, npz, ng, npx, npy, gridstruct%area, gridstruct%&
&                  rarea, dp_ref, zs, zh, crx, cry, xfx, yfx, delz, ws, &
&                  rdt, gridstruct, bd, flagstructp%hord_tm_pert)
        CALL TIMING_OFF('UPDATE_DZ')
        IF (flagstruct%fv_debug) THEN
          IF (.NOT.flagstruct%hydrostatic) CALL PRT_MXM('delz updated', &
&                                                 delz, is, ie, js, je, &
&                                                 ng, npz, 1., &
&                                                 gridstruct%area_64, &
&                                                 domain)
        END IF
        IF (idiag%id_ws .GT. 0 .AND. last_step) used = SEND_DATA(idiag%&
&           id_ws, ws, fv_time)
!           call prt_maxmin('WS', ws, is, ie, js, je, 0, 1, 1., master)
        CALL TIMING_ON('Riem_Solver')
        arg10 = beta .LT. -0.1
        CALL RIEM_SOLVER3(flagstruct%m_split, dt, is, ie, js, je, npz, &
&                   ng, isd, ied, jsd, jed, akap, cappa, cp, ptop, zs, &
&                   q_con, w, delz, pt, delp, zh, pe, pkc, pk3, pk, peln&
&                   , ws, flagstruct%scale_z, flagstruct%p_fac, &
&                   flagstruct%a_imp, flagstruct%use_logp, remap_step, &
&                   arg10)
        CALL TIMING_OFF('Riem_Solver')
        CALL TIMING_ON('COMM_TOTAL')
        IF (gridstruct%square_domain) THEN
          CALL START_GROUP_HALO_UPDATE(i_pack(4), zh, domain)
          CALL START_GROUP_HALO_UPDATE(i_pack(5), pkc, domain, whalo=2, &
&                                ehalo=2, shalo=2, nhalo=2)
        ELSE
          CALL START_GROUP_HALO_UPDATE(i_pack(4), zh, domain, complete=&
&                                .true.)
          CALL START_GROUP_HALO_UPDATE(i_pack(4), pkc, domain, complete=&
&                                .true.)
        END IF
        CALL TIMING_OFF('COMM_TOTAL')
        IF (remap_step) CALL PE_HALO(is, ie, js, je, isd, ied, jsd, jed&
&                              , npz, ptop, pe, delp)
        IF (flagstruct%use_logp) THEN
          CALL PLN_HALO(is, ie, js, je, isd, ied, jsd, jed, npz, ptop, &
&                 pk3, delp)
        ELSE
          CALL PK3_HALO(is, ie, js, je, isd, ied, jsd, jed, npz, ptop, &
&                 akap, pk3, delp)
        END IF
        IF (gridstruct%nested) THEN
          arg1 = split_timestep_bc + 1.
          arg2 = REAL(n_split*flagstruct%k_split)
          CALL NESTED_GRID_BC_APPLY_INTT(delz, 0, 0, npx, npy, npz, bd, &
&                                  arg1, arg2, neststruct%delz_bc, &
&                                  neststruct%nestbctype)
!Compute gz/pkc/pk3; note that now pkc should be nonhydro pert'n pressure
          CALL NEST_HALO_NH(ptop, grav, akap, cp, delp, delz, pt, phis, &
&                     pkc, gz, pk3, npx, npy, npz, gridstruct%nested, &
&                     .true., .true., .true., bd)
        END IF
        CALL TIMING_ON('COMM_TOTAL')
        CALL COMPLETE_GROUP_HALO_UPDATE(i_pack(4), domain)
        CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gz,zh,grav)
        DO k=1,npz+1
          DO j=js-2,je+2
            DO i=is-2,ie+2
              gz(i, j, k) = zh(i, j, k)*grav
            END DO
          END DO
        END DO
        IF (gridstruct%square_domain) THEN
          CALL TIMING_ON('COMM_TOTAL')
          CALL COMPLETE_GROUP_HALO_UPDATE(i_pack(5), domain)
          CALL TIMING_OFF('COMM_TOTAL')
        END IF
      END IF
      IF (remap_step .AND. hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pk,pkc)
        DO k=1,npz+1
          DO j=js,je
            DO i=is,ie
              pk(i, j, k) = pkc(i, j, k)
            END DO
          END DO
        END DO
      END IF
!----------------------------
! Compute pressure gradient:
!----------------------------
      CALL TIMING_ON('PG_D')
      IF (hydrostatic) THEN
        IF (beta .GT. 0.) THEN
          CALL GRAD1_P_UPDATE(divg2, u, v, pkc, gz, du, dv, dt, ng, &
&                       gridstruct, bd, npx, npy, npz, ptop, beta_d, &
&                       flagstruct%a2b_ord)
        ELSE
          CALL ONE_GRAD_P(u, v, pkc, gz, divg2, delp, dt, ng, gridstruct&
&                   , bd, npx, npy, npz, ptop, hydrostatic, flagstruct%&
&                   a2b_ord, flagstruct%d_ext)
        END IF
      ELSE IF (beta .GT. 0.) THEN
        CALL SPLIT_P_GRAD(u, v, pkc, gz, du, dv, delp, pk3, beta_d, dt, &
&                   ng, gridstruct, bd, npx, npy, npz, flagstruct%&
&                   use_logp)
      ELSE IF (beta .LT. -0.1) THEN
        CALL ONE_GRAD_P(u, v, pkc, gz, divg2, delp, dt, ng, gridstruct, &
&                 bd, npx, npy, npz, ptop, hydrostatic, flagstruct%&
&                 a2b_ord, flagstruct%d_ext)
      ELSE
        CALL NH_P_GRAD(u, v, pkc, gz, delp, pk3, dt, ng, gridstruct, bd&
&                , npx, npy, npz, flagstruct%use_logp)
      END IF
      CALL TIMING_OFF('PG_D')
! Inline Rayleigh friction here?
!-------------------------------------------------------------------------------------------------------
      IF (flagstruct%breed_vortex_inline) THEN
        IF (.NOT.hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pkz,cappa,rdg,delp,delz,pt,k1k)
          DO k=1,npz
            DO j=js,je
              DO i=is,ie
! Note: pt at this stage is Theta_m
                pkz(i, j, k) = EXP(k1k*LOG(rdg*delp(i, j, k)/delz(i, j, &
&                 k)*pt(i, j, k)))
              END DO
            END DO
          END DO
        END IF
        CALL BREED_SLP_INLINE(it, dt, npz, ak, bk, phis, pe, pk, peln, &
&                       pkz, delp, u, v, pt, q, flagstruct%nwat, zvir, &
&                       gridstruct, ks, domain, bd, hydrostatic)
      END IF
!-------------------------------------------------------------------------------------------------------
      CALL TIMING_ON('COMM_TOTAL')
      IF (it .EQ. n_split .AND. gridstruct%grid_type .LT. 4 .AND. (.NOT.&
&         gridstruct%nested)) THEN
! Prevent accumulation of rounding errors at overlapped domain edges:
        CALL MPP_GET_BOUNDARY(u, v, domain, ebuffery=ebuffer, nbufferx=&
&                       nbuffer, gridtype=dgrid_ne)
!$OMP parallel do default(none) shared(is,ie,js,je,npz,u,nbuffer,v,ebuffer)
        DO k=1,npz
          DO i=is,ie
            u(i, je+1, k) = nbuffer(i-is+1, k)
          END DO
          DO j=js,je
            v(ie+1, j, k) = ebuffer(j-js+1, k)
          END DO
        END DO
      END IF
      IF (it .NE. n_split) CALL START_GROUP_HALO_UPDATE(i_pack(8), u, v&
&                                                 , domain, gridtype=&
&                                                 dgrid_ne)
      CALL TIMING_OFF('COMM_TOTAL')
      IF (gridstruct%nested) neststruct%nest_timestep = neststruct%&
&         nest_timestep + 1
      IF (hydrostatic .AND. last_step) THEN
        IF (flagstruct%use_old_omega) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga,pe,pem,rdt)
          DO k=1,npz
            DO j=js,je
              DO i=is,ie
                omga(i, j, k) = (pe(i, k+1, j)-pem(i, k+1, j))*rdt
              END DO
            END DO
          END DO
!------------------------------
! Compute the "advective term"
!------------------------------
          CALL ADV_PE(ua, va, pem, omga, gridstruct, bd, npx, npy, npz, &
&               ng)
        ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga) private(om2d)
          DO j=js,je
            DO k=1,npz
              DO i=is,ie
                om2d(i, k) = omga(i, j, k)
              END DO
            END DO
            DO k=2,npz
              DO i=is,ie
                om2d(i, k) = om2d(i, k-1) + omga(i, j, k)
              END DO
            END DO
            DO k=2,npz
              DO i=is,ie
                omga(i, j, k) = om2d(i, k)
              END DO
            END DO
          END DO
        END IF
        IF (idiag%id_ws .GT. 0 .AND. hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ws,delz,delp,omga)
          DO j=js,je
            DO i=is,ie
              ws(i, j) = delz(i, j, npz)/delp(i, j, npz)*omga(i, j, npz)
            END DO
          END DO
          used = SEND_DATA(idiag%id_ws, ws, fv_time)
        END IF
      END IF
      IF (gridstruct%nested) THEN
        IF (.NOT.hydrostatic) THEN
          arg1 = split_timestep_bc + 1
          arg2 = REAL(n_split*flagstruct%k_split)
          CALL NESTED_GRID_BC_APPLY_INTT(w, 0, 0, npx, npy, npz, bd, &
&                                  arg1, arg2, neststruct%w_bc, &
&                                  neststruct%nestbctype)
        END IF
        arg1 = split_timestep_bc + 1
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(u, 0, 1, npx, npy, npz, bd, arg1&
&                                , arg2, neststruct%u_bc, neststruct%&
&                                nestbctype)
        arg1 = split_timestep_bc + 1
        arg2 = REAL(n_split*flagstruct%k_split)
        CALL NESTED_GRID_BC_APPLY_INTT(v, 1, 0, npx, npy, npz, bd, arg1&
&                                , arg2, neststruct%v_bc, neststruct%&
&                                nestbctype)
      END IF
    END DO
!-----------------------------------------------------
! time split loop
!-----------------------------------------------------
    IF (nq .GT. 0 .AND. (.NOT.flagstruct%inline_q)) THEN
      CALL TIMING_ON('COMM_TOTAL')
      CALL TIMING_ON('COMM_TRACER')
      CALL START_GROUP_HALO_UPDATE(i_pack(10), q, domain)
      CALL TIMING_OFF('COMM_TRACER')
      CALL TIMING_OFF('COMM_TOTAL')
    END IF
    IF (flagstruct%fv_debug) THEN
      IF (IS_MASTER()) WRITE(*, *) 'End of n_split loop'
    END IF
    IF (n_con .NE. 0 .AND. flagstruct%d_con .GT. 1.e-5) THEN
      IF (3 .GT. flagstruct%nord + 1) THEN
        nf_ke = flagstruct%nord + 1
      ELSE
        nf_ke = 3
      END IF
      arg11 = cnst_0p20*gridstruct%da_min
      CALL DEL2_CUBED(heat_source, arg11, gridstruct, domain, npx, npy, &
&               npz, nf_ke, bd)
! Note: pt here is cp*(Virtual_Temperature/pkz)
      IF (hydrostatic) THEN
!
! del(Cp*T) = - del(KE)
!
!$OMP parallel do default(none) shared(flagstruct,is,ie,js,je,n_con,pt,heat_source,delp,pkz,bdt) &
!$OMP                          private(dtmp)
        DO j=js,je
! n_con is usually less than 3;
          DO k=1,n_con
            IF (k .LT. 3) THEN
              DO i=is,ie
                pt(i, j, k) = pt(i, j, k) + heat_source(i, j, k)/(cp_air&
&                 *delp(i, j, k)*pkz(i, j, k))
              END DO
            ELSE
              DO i=is,ie
                dtmp = heat_source(i, j, k)/(cp_air*delp(i, j, k))
                IF (bdt .GE. 0.) THEN
                  abs0 = bdt
                ELSE
                  abs0 = -bdt
                END IF
                x1 = abs0*flagstruct%delt_max
                IF (dtmp .GE. 0.) THEN
                  y1 = dtmp
                ELSE
                  y1 = -dtmp
                END IF
                IF (x1 .GT. y1) THEN
                  min1 = y1
                ELSE
                  min1 = x1
                END IF
                pt(i, j, k) = pt(i, j, k) + SIGN(min1, dtmp)/pkz(i, j, k&
&                 )
              END DO
            END IF
          END DO
        END DO
      ELSE
!$OMP parallel do default(none) shared(flagstruct,is,ie,js,je,n_con,pkz,cappa,rdg,delp,delz,pt, &
!$OMP                                  heat_source,k1k,cv_air,bdt) &
!$OMP                          private(dtmp, delt)
        DO k=1,n_con
          IF (bdt*flagstruct%delt_max .GE. 0.) THEN
            delt = bdt*flagstruct%delt_max
          ELSE
            delt = -(bdt*flagstruct%delt_max)
          END IF
! Sponge layers:
!         if ( k == 1 ) delt = 2.0*delt
!         if ( k == 2 ) delt = 1.5*delt
          DO j=js,je
            DO i=is,ie
              pkz(i, j, k) = EXP(k1k*LOG(rdg*delp(i, j, k)/delz(i, j, k)&
&               *pt(i, j, k)))
              dtmp = heat_source(i, j, k)/(cv_air*delp(i, j, k))
              IF (dtmp .GE. 0.) THEN
                y2 = dtmp
              ELSE
                y2 = -dtmp
              END IF
              IF (delt .GT. y2) THEN
                min2 = y2
              ELSE
                min2 = delt
              END IF
              pt(i, j, k) = pt(i, j, k) + SIGN(min2, dtmp)/pkz(i, j, k)
            END DO
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE DYN_CORE
!  Differentiation of pk3_halo in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
!.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mo
!d.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix
!_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_S
!uper fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 
!fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rema
!p_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_
!mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv
!_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_res
!tart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_
!z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mo
!d.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.
!SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.n
!est_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2
!c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
! sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mo
!d.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_m
!od.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pk3 delp
!   with respect to varying inputs: pk3 delp
  SUBROUTINE PK3_HALO_FWD(is, ie, js, je, isd, ied, jsd, jed, npz, ptop&
&   , akap, pk3, delp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop, akap
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3
! Local:
    REAL :: pei(isd:ied)
    REAL :: pej(jsd:jed)
    INTEGER :: i, j, k
    INTRINSIC LOG
    INTRINSIC EXP

    pei = 0.0
    pej = 0.0

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,npz,ptop,delp,pk3,akap) &
!$OMP                          private(pei)
    DO j=js,je
      CALL PUSHREALARRAY(pei(is-2))
      pei(is-2) = ptop
      CALL PUSHREALARRAY(pei(is-1))
      pei(is-1) = ptop
      DO k=1,npz
        CALL PUSHREALARRAY(pei(is-2))
        pei(is-2) = pei(is-2) + delp(is-2, j, k)
        CALL PUSHREALARRAY(pei(is-1))
        pei(is-1) = pei(is-1) + delp(is-1, j, k)
        CALL PUSHREALARRAY(pk3(is-2, j, k+1))
        pk3(is-2, j, k+1) = EXP(akap*LOG(pei(is-2)))
        CALL PUSHREALARRAY(pk3(is-1, j, k+1))
        pk3(is-1, j, k+1) = EXP(akap*LOG(pei(is-1)))
      END DO
      CALL PUSHREALARRAY(pei(ie+1))
      pei(ie+1) = ptop
      CALL PUSHREALARRAY(pei(ie+2))
      pei(ie+2) = ptop
      DO k=1,npz
        CALL PUSHREALARRAY(pei(ie+1))
        pei(ie+1) = pei(ie+1) + delp(ie+1, j, k)
        CALL PUSHREALARRAY(pei(ie+2))
        pei(ie+2) = pei(ie+2) + delp(ie+2, j, k)
        CALL PUSHREALARRAY(pk3(ie+1, j, k+1))
        pk3(ie+1, j, k+1) = EXP(akap*LOG(pei(ie+1)))
        CALL PUSHREALARRAY(pk3(ie+2, j, k+1))
        pk3(ie+2, j, k+1) = EXP(akap*LOG(pei(ie+2)))
      END DO
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ptop,delp,pk3,akap) &
!$OMP                          private(pej)
    DO i=is-2,ie+2
      CALL PUSHREALARRAY(pej(js-2))
      pej(js-2) = ptop
      CALL PUSHREALARRAY(pej(js-1))
      pej(js-1) = ptop
      DO k=1,npz
        CALL PUSHREALARRAY(pej(js-2))
        pej(js-2) = pej(js-2) + delp(i, js-2, k)
        CALL PUSHREALARRAY(pej(js-1))
        pej(js-1) = pej(js-1) + delp(i, js-1, k)
        CALL PUSHREALARRAY(pk3(i, js-2, k+1))
        pk3(i, js-2, k+1) = EXP(akap*LOG(pej(js-2)))
        CALL PUSHREALARRAY(pk3(i, js-1, k+1))
        pk3(i, js-1, k+1) = EXP(akap*LOG(pej(js-1)))
      END DO
      CALL PUSHREALARRAY(pej(je+1))
      pej(je+1) = ptop
      CALL PUSHREALARRAY(pej(je+2))
      pej(je+2) = ptop
      DO k=1,npz
        CALL PUSHREALARRAY(pej(je+1))
        pej(je+1) = pej(je+1) + delp(i, je+1, k)
        CALL PUSHREALARRAY(pej(je+2))
        pej(je+2) = pej(je+2) + delp(i, je+2, k)
        CALL PUSHREALARRAY(pk3(i, je+1, k+1))
        pk3(i, je+1, k+1) = EXP(akap*LOG(pej(je+1)))
        CALL PUSHREALARRAY(pk3(i, je+2, k+1))
        pk3(i, je+2, k+1) = EXP(akap*LOG(pej(je+2)))
      END DO
    END DO
    CALL PUSHREALARRAY(pej, jed - jsd + 1)
    CALL PUSHREALARRAY(pei, ied - isd + 1)
  END SUBROUTINE PK3_HALO_FWD
!  Differentiation of pk3_halo in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
!d.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_m
!od.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mi
!x_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_
!Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4
! fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rem
!ap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv
!_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters f
!v_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_re
!start_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid
!_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_m
!od.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod
!.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.
!nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a
!2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v_f
!b sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_m
!od.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_
!mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pk3 delp
!   with respect to varying inputs: pk3 delp
  SUBROUTINE PK3_HALO_BWD(is, ie, js, je, isd, ied, jsd, jed, npz, ptop&
&   , akap, pk3, pk3_ad, delp, delp_ad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop, akap
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz) :: delp_ad
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3_ad
    REAL :: pei(isd:ied)
    REAL :: pei_ad(isd:ied)
    REAL :: pej(jsd:jed)
    REAL :: pej_ad(jsd:jed)
    INTEGER :: i, j, k
    INTRINSIC LOG
    INTRINSIC EXP

    pei = 0.0
    pej = 0.0

    CALL POPREALARRAY(pei, ied - isd + 1)
    CALL POPREALARRAY(pej, jed - jsd + 1)
    pej_ad = 0.0
    DO i=ie+2,is-2,-1
      DO k=npz,1,-1
        CALL POPREALARRAY(pk3(i, je+2, k+1))
        pej_ad(je+2) = pej_ad(je+2) + akap*EXP(akap*LOG(pej(je+2)))*&
&         pk3_ad(i, je+2, k+1)/pej(je+2)
        pk3_ad(i, je+2, k+1) = 0.0
        CALL POPREALARRAY(pk3(i, je+1, k+1))
        pej_ad(je+1) = pej_ad(je+1) + akap*EXP(akap*LOG(pej(je+1)))*&
&         pk3_ad(i, je+1, k+1)/pej(je+1)
        pk3_ad(i, je+1, k+1) = 0.0
        CALL POPREALARRAY(pej(je+2))
        delp_ad(i, je+2, k) = delp_ad(i, je+2, k) + pej_ad(je+2)
        CALL POPREALARRAY(pej(je+1))
        delp_ad(i, je+1, k) = delp_ad(i, je+1, k) + pej_ad(je+1)
      END DO
      CALL POPREALARRAY(pej(je+2))
      pej_ad(je+2) = 0.0
      CALL POPREALARRAY(pej(je+1))
      pej_ad(je+1) = 0.0
      DO k=npz,1,-1
        CALL POPREALARRAY(pk3(i, js-1, k+1))
        pej_ad(js-1) = pej_ad(js-1) + akap*EXP(akap*LOG(pej(js-1)))*&
&         pk3_ad(i, js-1, k+1)/pej(js-1)
        pk3_ad(i, js-1, k+1) = 0.0
        CALL POPREALARRAY(pk3(i, js-2, k+1))
        pej_ad(js-2) = pej_ad(js-2) + akap*EXP(akap*LOG(pej(js-2)))*&
&         pk3_ad(i, js-2, k+1)/pej(js-2)
        pk3_ad(i, js-2, k+1) = 0.0
        CALL POPREALARRAY(pej(js-1))
        delp_ad(i, js-1, k) = delp_ad(i, js-1, k) + pej_ad(js-1)
        CALL POPREALARRAY(pej(js-2))
        delp_ad(i, js-2, k) = delp_ad(i, js-2, k) + pej_ad(js-2)
      END DO
      CALL POPREALARRAY(pej(js-1))
      pej_ad(js-1) = 0.0
      CALL POPREALARRAY(pej(js-2))
      pej_ad(js-2) = 0.0
    END DO
    pei_ad = 0.0
    DO j=je,js,-1
      DO k=npz,1,-1
        CALL POPREALARRAY(pk3(ie+2, j, k+1))
        pei_ad(ie+2) = pei_ad(ie+2) + akap*EXP(akap*LOG(pei(ie+2)))*&
&         pk3_ad(ie+2, j, k+1)/pei(ie+2)
        pk3_ad(ie+2, j, k+1) = 0.0
        CALL POPREALARRAY(pk3(ie+1, j, k+1))
        pei_ad(ie+1) = pei_ad(ie+1) + akap*EXP(akap*LOG(pei(ie+1)))*&
&         pk3_ad(ie+1, j, k+1)/pei(ie+1)
        pk3_ad(ie+1, j, k+1) = 0.0
        CALL POPREALARRAY(pei(ie+2))
        delp_ad(ie+2, j, k) = delp_ad(ie+2, j, k) + pei_ad(ie+2)
        CALL POPREALARRAY(pei(ie+1))
        delp_ad(ie+1, j, k) = delp_ad(ie+1, j, k) + pei_ad(ie+1)
      END DO
      CALL POPREALARRAY(pei(ie+2))
      pei_ad(ie+2) = 0.0
      CALL POPREALARRAY(pei(ie+1))
      pei_ad(ie+1) = 0.0
      DO k=npz,1,-1
        CALL POPREALARRAY(pk3(is-1, j, k+1))
        pei_ad(is-1) = pei_ad(is-1) + akap*EXP(akap*LOG(pei(is-1)))*&
&         pk3_ad(is-1, j, k+1)/pei(is-1)
        pk3_ad(is-1, j, k+1) = 0.0
        CALL POPREALARRAY(pk3(is-2, j, k+1))
        pei_ad(is-2) = pei_ad(is-2) + akap*EXP(akap*LOG(pei(is-2)))*&
&         pk3_ad(is-2, j, k+1)/pei(is-2)
        pk3_ad(is-2, j, k+1) = 0.0
        CALL POPREALARRAY(pei(is-1))
        delp_ad(is-1, j, k) = delp_ad(is-1, j, k) + pei_ad(is-1)
        CALL POPREALARRAY(pei(is-2))
        delp_ad(is-2, j, k) = delp_ad(is-2, j, k) + pei_ad(is-2)
      END DO
      CALL POPREALARRAY(pei(is-1))
      pei_ad(is-1) = 0.0
      CALL POPREALARRAY(pei(is-2))
      pei_ad(is-2) = 0.0
    END DO
  END SUBROUTINE PK3_HALO_BWD
  SUBROUTINE PK3_HALO(is, ie, js, je, isd, ied, jsd, jed, npz, ptop, &
&   akap, pk3, delp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop, akap
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3
! Local:
    REAL :: pei(isd:ied)
    REAL :: pej(jsd:jed)
    INTEGER :: i, j, k
    INTRINSIC LOG
    INTRINSIC EXP
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,npz,ptop,delp,pk3,akap) &
!$OMP                          private(pei)
    DO j=js,je
      pei(is-2) = ptop
      pei(is-1) = ptop
      DO k=1,npz
        pei(is-2) = pei(is-2) + delp(is-2, j, k)
        pei(is-1) = pei(is-1) + delp(is-1, j, k)
        pk3(is-2, j, k+1) = EXP(akap*LOG(pei(is-2)))
        pk3(is-1, j, k+1) = EXP(akap*LOG(pei(is-1)))
      END DO
      pei(ie+1) = ptop
      pei(ie+2) = ptop
      DO k=1,npz
        pei(ie+1) = pei(ie+1) + delp(ie+1, j, k)
        pei(ie+2) = pei(ie+2) + delp(ie+2, j, k)
        pk3(ie+1, j, k+1) = EXP(akap*LOG(pei(ie+1)))
        pk3(ie+2, j, k+1) = EXP(akap*LOG(pei(ie+2)))
      END DO
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ptop,delp,pk3,akap) &
!$OMP                          private(pej)
    DO i=is-2,ie+2
      pej(js-2) = ptop
      pej(js-1) = ptop
      DO k=1,npz
        pej(js-2) = pej(js-2) + delp(i, js-2, k)
        pej(js-1) = pej(js-1) + delp(i, js-1, k)
        pk3(i, js-2, k+1) = EXP(akap*LOG(pej(js-2)))
        pk3(i, js-1, k+1) = EXP(akap*LOG(pej(js-1)))
      END DO
      pej(je+1) = ptop
      pej(je+2) = ptop
      DO k=1,npz
        pej(je+1) = pej(je+1) + delp(i, je+1, k)
        pej(je+2) = pej(je+2) + delp(i, je+2, k)
        pk3(i, je+1, k+1) = EXP(akap*LOG(pej(je+1)))
        pk3(i, je+2, k+1) = EXP(akap*LOG(pej(je+2)))
      END DO
    END DO
  END SUBROUTINE PK3_HALO
!  Differentiation of pln_halo in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
!.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mo
!d.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix
!_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_S
!uper fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 
!fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rema
!p_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_
!mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv
!_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_res
!tart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_
!z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mo
!d.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.
!SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.n
!est_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2
!c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
! sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mo
!d.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_m
!od.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pk3 delp
!   with respect to varying inputs: pk3 delp
  SUBROUTINE PLN_HALO_FWD(is, ie, js, je, isd, ied, jsd, jed, npz, ptop&
&   , pk3, delp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3
! Local:
    REAL :: pet
    INTEGER :: i, j, k
    INTRINSIC LOG

    pet = 0.0

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,npz,ptop,delp,pk3) &
!$OMP                          private(pet)
    DO j=js,je
      DO i=is-2,is-1
        CALL PUSHREALARRAY(pet)
        pet = ptop
        DO k=1,npz
          CALL PUSHREALARRAY(pet)
          pet = pet + delp(i, j, k)
          CALL PUSHREALARRAY(pk3(i, j, k+1))
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
      DO i=ie+1,ie+2
        CALL PUSHREALARRAY(pet)
        pet = ptop
        DO k=1,npz
          CALL PUSHREALARRAY(pet)
          pet = pet + delp(i, j, k)
          CALL PUSHREALARRAY(pk3(i, j, k+1))
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ptop,delp,pk3) &
!$OMP                          private(pet)
    DO i=is-2,ie+2
      DO j=js-2,js-1
        CALL PUSHREALARRAY(pet)
        pet = ptop
        DO k=1,npz
          CALL PUSHREALARRAY(pet)
          pet = pet + delp(i, j, k)
          CALL PUSHREALARRAY(pk3(i, j, k+1))
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
      DO j=je+1,je+2
        CALL PUSHREALARRAY(pet)
        pet = ptop
        DO k=1,npz
          CALL PUSHREALARRAY(pet)
          pet = pet + delp(i, j, k)
          CALL PUSHREALARRAY(pk3(i, j, k+1))
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
    END DO
    CALL PUSHREALARRAY(pet)
  END SUBROUTINE PLN_HALO_FWD
!  Differentiation of pln_halo in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
!d.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_m
!od.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mi
!x_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_
!Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4
! fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rem
!ap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv
!_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters f
!v_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_re
!start_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid
!_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_m
!od.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod
!.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.
!nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a
!2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v_f
!b sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_m
!od.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_
!mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: pk3 delp
!   with respect to varying inputs: pk3 delp
  SUBROUTINE PLN_HALO_BWD(is, ie, js, je, isd, ied, jsd, jed, npz, ptop&
&   , pk3, pk3_ad, delp, delp_ad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz) :: delp_ad
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3_ad
    REAL :: pet
    REAL :: pet_ad
    INTEGER :: i, j, k
    INTRINSIC LOG

    pet = 0.0

    CALL POPREALARRAY(pet)
    DO i=ie+2,is-2,-1
      DO j=je+2,je+1,-1
        pet_ad = 0.0
        DO k=npz,1,-1
          CALL POPREALARRAY(pk3(i, j, k+1))
          pet_ad = pet_ad + pk3_ad(i, j, k+1)/pet
          pk3_ad(i, j, k+1) = 0.0
          CALL POPREALARRAY(pet)
          delp_ad(i, j, k) = delp_ad(i, j, k) + pet_ad
        END DO
        CALL POPREALARRAY(pet)
      END DO
      DO j=js-1,js-2,-1
        pet_ad = 0.0
        DO k=npz,1,-1
          CALL POPREALARRAY(pk3(i, j, k+1))
          pet_ad = pet_ad + pk3_ad(i, j, k+1)/pet
          pk3_ad(i, j, k+1) = 0.0
          CALL POPREALARRAY(pet)
          delp_ad(i, j, k) = delp_ad(i, j, k) + pet_ad
        END DO
        CALL POPREALARRAY(pet)
      END DO
    END DO
    DO j=je,js,-1
      DO i=ie+2,ie+1,-1
        pet_ad = 0.0
        DO k=npz,1,-1
          CALL POPREALARRAY(pk3(i, j, k+1))
          pet_ad = pet_ad + pk3_ad(i, j, k+1)/pet
          pk3_ad(i, j, k+1) = 0.0
          CALL POPREALARRAY(pet)
          delp_ad(i, j, k) = delp_ad(i, j, k) + pet_ad
        END DO
        CALL POPREALARRAY(pet)
      END DO
      DO i=is-1,is-2,-1
        pet_ad = 0.0
        DO k=npz,1,-1
          CALL POPREALARRAY(pk3(i, j, k+1))
          pet_ad = pet_ad + pk3_ad(i, j, k+1)/pet
          pk3_ad(i, j, k+1) = 0.0
          CALL POPREALARRAY(pet)
          delp_ad(i, j, k) = delp_ad(i, j, k) + pet_ad
        END DO
        CALL POPREALARRAY(pet)
      END DO
    END DO
  END SUBROUTINE PLN_HALO_BWD
  SUBROUTINE PLN_HALO(is, ie, js, je, isd, ied, jsd, jed, npz, ptop, pk3&
&   , delp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3
! Local:
    REAL :: pet
    INTEGER :: i, j, k
    INTRINSIC LOG
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,npz,ptop,delp,pk3) &
!$OMP                          private(pet)
    DO j=js,je
      DO i=is-2,is-1
        pet = ptop
        DO k=1,npz
          pet = pet + delp(i, j, k)
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
      DO i=ie+1,ie+2
        pet = ptop
        DO k=1,npz
          pet = pet + delp(i, j, k)
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ptop,delp,pk3) &
!$OMP                          private(pet)
    DO i=is-2,ie+2
      DO j=js-2,js-1
        pet = ptop
        DO k=1,npz
          pet = pet + delp(i, j, k)
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
      DO j=je+1,je+2
        pet = ptop
        DO k=1,npz
          pet = pet + delp(i, j, k)
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
    END DO
  END SUBROUTINE PLN_HALO
!  Differentiation of pe_halo in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.
!a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod
!.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_
!dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Su
!per fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 f
!v_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap
!_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_m
!apz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_
!mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_rest
!art_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z
! main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod
!.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.S
!IM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.ne
!st_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c
!_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v 
!sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod
!.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mo
!d.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: delp pe
!   with respect to varying inputs: delp pe
  SUBROUTINE PE_HALO_FWD(is, ie, js, je, isd, ied, jsd, jed, npz, ptop, &
&   pe, delp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(is-1:ie+1, npz+1, js-1:je+1), INTENT(INOUT) :: pe
! Local:
    INTEGER :: i, j, k
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pe,delp,ptop)
    DO j=js,je
      CALL PUSHREALARRAY(pe(is-1, 1, j))
      pe(is-1, 1, j) = ptop
      CALL PUSHREALARRAY(pe(ie+1, 1, j))
      pe(ie+1, 1, j) = ptop
      DO k=1,npz
        CALL PUSHREALARRAY(pe(is-1, k+1, j))
        pe(is-1, k+1, j) = pe(is-1, k, j) + delp(is-1, j, k)
        CALL PUSHREALARRAY(pe(ie+1, k+1, j))
        pe(ie+1, k+1, j) = pe(ie+1, k, j) + delp(ie+1, j, k)
      END DO
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pe,delp,ptop)
    DO i=is-1,ie+1
      CALL PUSHREALARRAY(pe(i, 1, js-1))
      pe(i, 1, js-1) = ptop
      CALL PUSHREALARRAY(pe(i, 1, je+1))
      pe(i, 1, je+1) = ptop
      DO k=1,npz
        CALL PUSHREALARRAY(pe(i, k+1, js-1))
        pe(i, k+1, js-1) = pe(i, k, js-1) + delp(i, js-1, k)
        CALL PUSHREALARRAY(pe(i, k+1, je+1))
        pe(i, k+1, je+1) = pe(i, k, je+1) + delp(i, je+1, k)
      END DO
    END DO
  END SUBROUTINE PE_HALO_FWD
!  Differentiation of pe_halo in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
!.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mo
!d.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix
!_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_S
!uper fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 
!fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rema
!p_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_
!mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv
!_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_res
!tart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_
!z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mo
!d.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.
!SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.n
!est_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2
!c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
! sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mo
!d.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_m
!od.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: delp pe
!   with respect to varying inputs: delp pe
  SUBROUTINE PE_HALO_BWD(is, ie, js, je, isd, ied, jsd, jed, npz, ptop, &
&   pe, pe_ad, delp, delp_ad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz) :: delp_ad
    REAL, DIMENSION(is-1:ie+1, npz+1, js-1:je+1), INTENT(INOUT) :: pe
    REAL, DIMENSION(is-1:ie+1, npz+1, js-1:je+1), INTENT(INOUT) :: pe_ad
    INTEGER :: i, j, k
    DO i=ie+1,is-1,-1
      DO k=npz,1,-1
        CALL POPREALARRAY(pe(i, k+1, je+1))
        pe_ad(i, k, je+1) = pe_ad(i, k, je+1) + pe_ad(i, k+1, je+1)
        delp_ad(i, je+1, k) = delp_ad(i, je+1, k) + pe_ad(i, k+1, je+1)
        pe_ad(i, k+1, je+1) = 0.0
        CALL POPREALARRAY(pe(i, k+1, js-1))
        pe_ad(i, k, js-1) = pe_ad(i, k, js-1) + pe_ad(i, k+1, js-1)
        delp_ad(i, js-1, k) = delp_ad(i, js-1, k) + pe_ad(i, k+1, js-1)
        pe_ad(i, k+1, js-1) = 0.0
      END DO
      CALL POPREALARRAY(pe(i, 1, je+1))
      pe_ad(i, 1, je+1) = 0.0
      CALL POPREALARRAY(pe(i, 1, js-1))
      pe_ad(i, 1, js-1) = 0.0
    END DO
    DO j=je,js,-1
      DO k=npz,1,-1
        CALL POPREALARRAY(pe(ie+1, k+1, j))
        pe_ad(ie+1, k, j) = pe_ad(ie+1, k, j) + pe_ad(ie+1, k+1, j)
        delp_ad(ie+1, j, k) = delp_ad(ie+1, j, k) + pe_ad(ie+1, k+1, j)
        pe_ad(ie+1, k+1, j) = 0.0
        CALL POPREALARRAY(pe(is-1, k+1, j))
        pe_ad(is-1, k, j) = pe_ad(is-1, k, j) + pe_ad(is-1, k+1, j)
        delp_ad(is-1, j, k) = delp_ad(is-1, j, k) + pe_ad(is-1, k+1, j)
        pe_ad(is-1, k+1, j) = 0.0
      END DO
      CALL POPREALARRAY(pe(ie+1, 1, j))
      pe_ad(ie+1, 1, j) = 0.0
      CALL POPREALARRAY(pe(is-1, 1, j))
      pe_ad(is-1, 1, j) = 0.0
    END DO
  END SUBROUTINE PE_HALO_BWD
  SUBROUTINE PE_HALO(is, ie, js, je, isd, ied, jsd, jed, npz, ptop, pe, &
&   delp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(is-1:ie+1, npz+1, js-1:je+1), INTENT(INOUT) :: pe
! Local:
    INTEGER :: i, j, k
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pe,delp,ptop)
    DO j=js,je
      pe(is-1, 1, j) = ptop
      pe(ie+1, 1, j) = ptop
      DO k=1,npz
        pe(is-1, k+1, j) = pe(is-1, k, j) + delp(is-1, j, k)
        pe(ie+1, k+1, j) = pe(ie+1, k, j) + delp(ie+1, j, k)
      END DO
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pe,delp,ptop)
    DO i=is-1,ie+1
      pe(i, 1, js-1) = ptop
      pe(i, 1, je+1) = ptop
      DO k=1,npz
        pe(i, k+1, js-1) = pe(i, k, js-1) + delp(i, js-1, k)
        pe(i, k+1, je+1) = pe(i, k, je+1) + delp(i, je+1, k)
      END DO
    END DO
  END SUBROUTINE PE_HALO
!  Differentiation of adv_pe in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a
!2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.
!p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_d
!p dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Sup
!er fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv
!_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_
!z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_ma
!pz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_m
!apz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_resta
!rt_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z 
!main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.
!Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SI
!M3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nes
!t_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_
!vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v s
!w_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.
!copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod
!.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: ua om va pem
!   with respect to varying inputs: ua om va pem
  SUBROUTINE ADV_PE_FWD(ua, va, pem, om, gridstruct, bd, npx, npy, npz, &
&   ng)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, npz, ng
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! Contra-variant wind components:
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(IN) :: ua&
&   , va
! Pressure at edges:
    REAL, INTENT(IN) :: pem(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(INOUT) :: om(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: up, vp
    REAL :: v3(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: pin(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pb(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: grad(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: pdx(3, bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: pdy(3, bd%is:bd%ie+1, bd%js:bd%je)
    INTEGER :: i, j, k, n
    INTEGER :: is, ie, js, je

    up = 0.0
    vp = 0.0
    v3 = 0.0
    pin = 0.0
    pb = 0.0
    grad = 0.0
    pdx = 0.0
    pdy = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ua,va,gridstruct,pem,npx,npy,ng,om) &
!$OMP                          private(n, pdx, pdy, pin, pb, up, vp, grad, v3)
    DO k=1,npz
      IF (k .EQ. npz) THEN
        DO j=js,je
          DO i=is,ie
            up(i, j) = ua(i, j, npz)
            vp(i, j) = va(i, j, npz)
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        DO j=js,je
          DO i=is,ie
            up(i, j) = 0.5*(ua(i, j, k)+ua(i, j, k+1))
            vp(i, j) = 0.5*(va(i, j, k)+va(i, j, k+1))
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      END IF
! Compute Vect wind:
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            CALL PUSHREALARRAY(v3(n, i, j))
            v3(n, i, j) = up(i, j)*gridstruct%ec1(n, i, j) + vp(i, j)*&
&             gridstruct%ec2(n, i, j)
          END DO
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          pin(i, j) = pem(i, k+1, j)
        END DO
      END DO
! Compute pe at 4 cell corners:
      CALL A2B_ORD2_FWD(pin, pb, gridstruct, npx, npy, is, ie, js, je, &
&                 ng)
      DO j=js,je+1
        DO i=is,ie
          DO n=1,3
            pdx(n, i, j) = (pb(i, j)+pb(i+1, j))*gridstruct%dx(i, j)*&
&             gridstruct%en1(n, i, j)
          END DO
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          DO n=1,3
            pdy(n, i, j) = (pb(i, j)+pb(i, j+1))*gridstruct%dy(i, j)*&
&             gridstruct%en2(n, i, j)
          END DO
        END DO
      END DO
! Compute grad (pe) by Green's theorem
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            CALL PUSHREALARRAY(grad(n, i, j))
            grad(n, i, j) = pdx(n, i, j+1) - pdx(n, i, j) - pdy(n, i, j)&
&             + pdy(n, i+1, j)
          END DO
        END DO
      END DO
! Compute inner product: V3 * grad (pe)
      DO j=js,je
        DO i=is,ie
          CALL PUSHREALARRAY(om(i, j, k))
          om(i, j, k) = om(i, j, k) + 0.5*gridstruct%rarea(i, j)*(v3(1, &
&           i, j)*grad(1, i, j)+v3(2, i, j)*grad(2, i, j)+v3(3, i, j)*&
&           grad(3, i, j))
        END DO
      END DO
    END DO
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(ie)
    CALL PUSHREALARRAY(grad, 3*(bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL PUSHREALARRAY(v3, 3*(bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL PUSHINTEGER(js)
  END SUBROUTINE ADV_PE_FWD
!  Differentiation of adv_pe in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.
!a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod
!.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_
!dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Su
!per fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 f
!v_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap
!_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_m
!apz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_
!mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_rest
!art_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z
! main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod
!.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.S
!IM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.ne
!st_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c
!_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v 
!sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod
!.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mo
!d.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: ua om va pem
!   with respect to varying inputs: ua om va pem
  SUBROUTINE ADV_PE_BWD(ua, ua_ad, va, va_ad, pem, pem_ad, om, om_ad, &
&   gridstruct, bd, npx, npy, npz, ng)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, npz, ng
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(IN) :: ua&
&   , va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz) :: ua_ad, va_ad
    REAL, INTENT(IN) :: pem(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL :: pem_ad(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(INOUT) :: om(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: om_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: up, vp
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: up_ad, vp_ad
    REAL :: v3(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: v3_ad(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: pin(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pin_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pb(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pb_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: grad(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: grad_ad(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: pdx(3, bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: pdx_ad(3, bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: pdy(3, bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: pdy_ad(3, bd%is:bd%ie+1, bd%js:bd%je)
    INTEGER :: i, j, k, n
    INTEGER :: is, ie, js, je
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    INTEGER :: branch

    up = 0.0
    vp = 0.0
    v3 = 0.0
    pin = 0.0
    pb = 0.0
    grad = 0.0
    pdx = 0.0
    pdy = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    branch = 0

    CALL POPINTEGER(js)
    CALL POPREALARRAY(v3, 3*(bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL POPREALARRAY(grad, 3*(bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL POPINTEGER(ie)
    CALL POPINTEGER(is)
    CALL POPINTEGER(je)
    v3_ad = 0.0
    grad_ad = 0.0
    up_ad = 0.0
    pdx_ad = 0.0
    pdy_ad = 0.0
    pb_ad = 0.0
    vp_ad = 0.0
    pin_ad = 0.0
    DO k=npz,1,-1
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(om(i, j, k))
          temp_ad1 = gridstruct%rarea(i, j)*0.5*om_ad(i, j, k)
          v3_ad(1, i, j) = v3_ad(1, i, j) + grad(1, i, j)*temp_ad1
          grad_ad(1, i, j) = grad_ad(1, i, j) + v3(1, i, j)*temp_ad1
          v3_ad(2, i, j) = v3_ad(2, i, j) + grad(2, i, j)*temp_ad1
          grad_ad(2, i, j) = grad_ad(2, i, j) + v3(2, i, j)*temp_ad1
          v3_ad(3, i, j) = v3_ad(3, i, j) + grad(3, i, j)*temp_ad1
          grad_ad(3, i, j) = grad_ad(3, i, j) + v3(3, i, j)*temp_ad1
        END DO
      END DO
      DO j=je,js,-1
        DO i=ie,is,-1
          DO n=3,1,-1
            CALL POPREALARRAY(grad(n, i, j))
            pdx_ad(n, i, j+1) = pdx_ad(n, i, j+1) + grad_ad(n, i, j)
            pdx_ad(n, i, j) = pdx_ad(n, i, j) - grad_ad(n, i, j)
            pdy_ad(n, i+1, j) = pdy_ad(n, i+1, j) + grad_ad(n, i, j)
            pdy_ad(n, i, j) = pdy_ad(n, i, j) - grad_ad(n, i, j)
            grad_ad(n, i, j) = 0.0
          END DO
        END DO
      END DO
      DO j=je,js,-1
        DO i=ie+1,is,-1
          DO n=3,1,-1
            temp_ad0 = gridstruct%dy(i, j)*gridstruct%en2(n, i, j)*&
&             pdy_ad(n, i, j)
            pb_ad(i, j) = pb_ad(i, j) + temp_ad0
            pb_ad(i, j+1) = pb_ad(i, j+1) + temp_ad0
            pdy_ad(n, i, j) = 0.0
          END DO
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          DO n=3,1,-1
            temp_ad = gridstruct%dx(i, j)*gridstruct%en1(n, i, j)*pdx_ad&
&             (n, i, j)
            pb_ad(i, j) = pb_ad(i, j) + temp_ad
            pb_ad(i+1, j) = pb_ad(i+1, j) + temp_ad
            pdx_ad(n, i, j) = 0.0
          END DO
        END DO
      END DO
      CALL A2B_ORD2_BWD(pin, pin_ad, pb, pb_ad, gridstruct, npx, npy, is&
&                 , ie, js, je, ng)
      DO j=je+1,js-1,-1
        DO i=ie+1,is-1,-1
          pem_ad(i, k+1, j) = pem_ad(i, k+1, j) + pin_ad(i, j)
          pin_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je,js,-1
        DO i=ie,is,-1
          DO n=3,1,-1
            CALL POPREALARRAY(v3(n, i, j))
            up_ad(i, j) = up_ad(i, j) + gridstruct%ec1(n, i, j)*v3_ad(n&
&             , i, j)
            vp_ad(i, j) = vp_ad(i, j) + gridstruct%ec2(n, i, j)*v3_ad(n&
&             , i, j)
            v3_ad(n, i, j) = 0.0
          END DO
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            va_ad(i, j, k) = va_ad(i, j, k) + 0.5*vp_ad(i, j)
            va_ad(i, j, k+1) = va_ad(i, j, k+1) + 0.5*vp_ad(i, j)
            vp_ad(i, j) = 0.0
            ua_ad(i, j, k) = ua_ad(i, j, k) + 0.5*up_ad(i, j)
            ua_ad(i, j, k+1) = ua_ad(i, j, k+1) + 0.5*up_ad(i, j)
            up_ad(i, j) = 0.0
          END DO
        END DO
      ELSE
        DO j=je,js,-1
          DO i=ie,is,-1
            va_ad(i, j, npz) = va_ad(i, j, npz) + vp_ad(i, j)
            vp_ad(i, j) = 0.0
            ua_ad(i, j, npz) = ua_ad(i, j, npz) + up_ad(i, j)
            up_ad(i, j) = 0.0
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE ADV_PE_BWD
  SUBROUTINE ADV_PE(ua, va, pem, om, gridstruct, bd, npx, npy, npz, ng)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, npz, ng
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! Contra-variant wind components:
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(IN) :: ua&
&   , va
! Pressure at edges:
    REAL, INTENT(IN) :: pem(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(INOUT) :: om(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: up, vp
    REAL :: v3(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: pin(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pb(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: grad(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: pdx(3, bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: pdy(3, bd%is:bd%ie+1, bd%js:bd%je)
    INTEGER :: i, j, k, n
    INTEGER :: is, ie, js, je
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ua,va,gridstruct,pem,npx,npy,ng,om) &
!$OMP                          private(n, pdx, pdy, pin, pb, up, vp, grad, v3)
    DO k=1,npz
      IF (k .EQ. npz) THEN
        DO j=js,je
          DO i=is,ie
            up(i, j) = ua(i, j, npz)
            vp(i, j) = va(i, j, npz)
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            up(i, j) = 0.5*(ua(i, j, k)+ua(i, j, k+1))
            vp(i, j) = 0.5*(va(i, j, k)+va(i, j, k+1))
          END DO
        END DO
      END IF
! Compute Vect wind:
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            v3(n, i, j) = up(i, j)*gridstruct%ec1(n, i, j) + vp(i, j)*&
&             gridstruct%ec2(n, i, j)
          END DO
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          pin(i, j) = pem(i, k+1, j)
        END DO
      END DO
! Compute pe at 4 cell corners:
      CALL A2B_ORD2(pin, pb, gridstruct, npx, npy, is, ie, js, je, ng)
      DO j=js,je+1
        DO i=is,ie
          DO n=1,3
            pdx(n, i, j) = (pb(i, j)+pb(i+1, j))*gridstruct%dx(i, j)*&
&             gridstruct%en1(n, i, j)
          END DO
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          DO n=1,3
            pdy(n, i, j) = (pb(i, j)+pb(i, j+1))*gridstruct%dy(i, j)*&
&             gridstruct%en2(n, i, j)
          END DO
        END DO
      END DO
! Compute grad (pe) by Green's theorem
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            grad(n, i, j) = pdx(n, i, j+1) - pdx(n, i, j) - pdy(n, i, j)&
&             + pdy(n, i+1, j)
          END DO
        END DO
      END DO
! Compute inner product: V3 * grad (pe)
      DO j=js,je
        DO i=is,ie
          om(i, j, k) = om(i, j, k) + 0.5*gridstruct%rarea(i, j)*(v3(1, &
&           i, j)*grad(1, i, j)+v3(2, i, j)*grad(2, i, j)+v3(3, i, j)*&
&           grad(3, i, j))
        END DO
      END DO
    END DO
  END SUBROUTINE ADV_PE
!  Differentiation of p_grad_c in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
!.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mo
!d.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix
!_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_S
!uper fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 
!fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rema
!p_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_
!mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv
!_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_res
!tart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_
!z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mo
!d.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.
!SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.n
!est_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2
!c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v
! sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mo
!d.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_m
!od.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: gz uc pkc delpc vc
!   with respect to varying inputs: gz uc pkc delpc vc
  SUBROUTINE P_GRAD_C_FWD(dt2, npz, delpc, pkc, gz, uc, vc, bd, rdxc, &
&   rdyc, hydrostatic)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npz
    REAL, INTENT(IN) :: dt2
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:, bd%jsd:, :), INTENT(IN) :: delpc
! pkc is pe**cappa     if hydrostatic
! pkc is full pressure if non-hydrostatic
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(IN) :: &
&   pkc, gz
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(IN) :: rdxc(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL, INTENT(IN) :: rdyc(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: hydrostatic
! Local:
    REAL :: wk(bd%is-1:bd%ie+1, bd%js-1:bd%je+1)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je

    wk = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
!$OMP parallel do default(none) shared(is,ie,js,je,npz,hydrostatic,pkc,delpc,uc,dt2,rdxc,gz,vc,rdyc) &
!$OMP                          private(wk)
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js-1,je+1
          DO i=is-1,ie+1
            CALL PUSHREALARRAY(wk(i, j))
            wk(i, j) = pkc(i, j, k+1) - pkc(i, j, k)
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        DO j=js-1,je+1
          DO i=is-1,ie+1
            CALL PUSHREALARRAY(wk(i, j))
            wk(i, j) = delpc(i, j, k)
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      END IF
      DO j=js,je
        DO i=is,ie+1
          uc(i, j, k) = uc(i, j, k) + dt2*rdxc(i, j)/(wk(i-1, j)+wk(i, j&
&           ))*((gz(i-1, j, k+1)-gz(i, j, k))*(pkc(i, j, k+1)-pkc(i-1, j&
&           , k))+(gz(i-1, j, k)-gz(i, j, k+1))*(pkc(i-1, j, k+1)-pkc(i&
&           , j, k)))
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          vc(i, j, k) = vc(i, j, k) + dt2*rdyc(i, j)/(wk(i, j-1)+wk(i, j&
&           ))*((gz(i, j-1, k+1)-gz(i, j, k))*(pkc(i, j, k+1)-pkc(i, j-1&
&           , k))+(gz(i, j-1, k)-gz(i, j, k+1))*(pkc(i, j-1, k+1)-pkc(i&
&           , j, k)))
        END DO
      END DO
    END DO
    CALL PUSHREALARRAY(wk, (bd%ie-bd%is+3)*(bd%je-bd%js+3))
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(ie)
    CALL PUSHINTEGER(js)
  END SUBROUTINE P_GRAD_C_FWD
!  Differentiation of p_grad_c in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
!d.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_m
!od.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mi
!x_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_
!Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4
! fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rem
!ap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv
!_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters f
!v_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_re
!start_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid
!_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_m
!od.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod
!.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.
!nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a
!2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v_f
!b sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_m
!od.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_
!mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: gz uc pkc delpc vc
!   with respect to varying inputs: gz uc pkc delpc vc
  SUBROUTINE P_GRAD_C_BWD(dt2, npz, delpc, delpc_ad, pkc, pkc_ad, gz, &
&   gz_ad, uc, uc_ad, vc, vc_ad, bd, rdxc, rdyc, hydrostatic)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npz
    REAL, INTENT(IN) :: dt2
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:, bd%jsd:, :), INTENT(IN) :: delpc
    REAL, DIMENSION(bd%isd:, bd%jsd:, :) :: delpc_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(IN) :: &
&   pkc, gz
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1) :: pkc_ad, &
&   gz_ad
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: uc_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: vc_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(IN) :: rdxc(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL, INTENT(IN) :: rdyc(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: hydrostatic
    REAL :: wk(bd%is-1:bd%ie+1, bd%js-1:bd%je+1)
    REAL :: wk_ad(bd%is-1:bd%ie+1, bd%js-1:bd%je+1)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp2
    REAL :: temp3
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp4
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp8
    REAL :: temp_ad1
    REAL :: temp_ad2
    INTEGER :: branch

    wk = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    branch = 0

    CALL POPINTEGER(js)
    CALL POPINTEGER(ie)
    CALL POPINTEGER(is)
    CALL POPINTEGER(je)
    CALL POPREALARRAY(wk, (bd%ie-bd%is+3)*(bd%je-bd%js+3))
    wk_ad = 0.0
    DO k=npz,1,-1
      DO j=je+1,js,-1
        DO i=ie,is,-1
          temp4 = wk(i, j-1) + wk(i, j)
          temp8 = pkc(i, j-1, k+1) - pkc(i, j, k)
          temp7 = gz(i, j-1, k) - gz(i, j, k+1)
          temp6 = pkc(i, j, k+1) - pkc(i, j-1, k)
          temp5 = gz(i, j-1, k+1) - gz(i, j, k)
          temp_ad1 = dt2*rdyc(i, j)*vc_ad(i, j, k)/temp4
          temp_ad2 = -((temp5*temp6+temp7*temp8)*temp_ad1/temp4)
          gz_ad(i, j-1, k+1) = gz_ad(i, j-1, k+1) + temp6*temp_ad1
          gz_ad(i, j, k) = gz_ad(i, j, k) - temp6*temp_ad1
          pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) + temp5*temp_ad1
          pkc_ad(i, j-1, k) = pkc_ad(i, j-1, k) - temp5*temp_ad1
          gz_ad(i, j-1, k) = gz_ad(i, j-1, k) + temp8*temp_ad1
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) - temp8*temp_ad1
          pkc_ad(i, j-1, k+1) = pkc_ad(i, j-1, k+1) + temp7*temp_ad1
          pkc_ad(i, j, k) = pkc_ad(i, j, k) - temp7*temp_ad1
          wk_ad(i, j-1) = wk_ad(i, j-1) + temp_ad2
          wk_ad(i, j) = wk_ad(i, j) + temp_ad2
        END DO
      END DO
      DO j=je,js,-1
        DO i=ie+1,is,-1
          temp = wk(i-1, j) + wk(i, j)
          temp3 = pkc(i-1, j, k+1) - pkc(i, j, k)
          temp2 = gz(i-1, j, k) - gz(i, j, k+1)
          temp1 = pkc(i, j, k+1) - pkc(i-1, j, k)
          temp0 = gz(i-1, j, k+1) - gz(i, j, k)
          temp_ad = dt2*rdxc(i, j)*uc_ad(i, j, k)/temp
          temp_ad0 = -((temp0*temp1+temp2*temp3)*temp_ad/temp)
          gz_ad(i-1, j, k+1) = gz_ad(i-1, j, k+1) + temp1*temp_ad
          gz_ad(i, j, k) = gz_ad(i, j, k) - temp1*temp_ad
          pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) + temp0*temp_ad
          pkc_ad(i-1, j, k) = pkc_ad(i-1, j, k) - temp0*temp_ad
          gz_ad(i-1, j, k) = gz_ad(i-1, j, k) + temp3*temp_ad
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) - temp3*temp_ad
          pkc_ad(i-1, j, k+1) = pkc_ad(i-1, j, k+1) + temp2*temp_ad
          pkc_ad(i, j, k) = pkc_ad(i, j, k) - temp2*temp_ad
          wk_ad(i-1, j) = wk_ad(i-1, j) + temp_ad0
          wk_ad(i, j) = wk_ad(i, j) + temp_ad0
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=je+1,js-1,-1
          DO i=ie+1,is-1,-1
            CALL POPREALARRAY(wk(i, j))
            delpc_ad(i, j, k) = delpc_ad(i, j, k) + wk_ad(i, j)
            wk_ad(i, j) = 0.0
          END DO
        END DO
      ELSE
        DO j=je+1,js-1,-1
          DO i=ie+1,is-1,-1
            CALL POPREALARRAY(wk(i, j))
            pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) + wk_ad(i, j)
            pkc_ad(i, j, k) = pkc_ad(i, j, k) - wk_ad(i, j)
            wk_ad(i, j) = 0.0
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE P_GRAD_C_BWD
  SUBROUTINE P_GRAD_C(dt2, npz, delpc, pkc, gz, uc, vc, bd, rdxc, rdyc, &
&   hydrostatic)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npz
    REAL, INTENT(IN) :: dt2
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:, bd%jsd:, :), INTENT(IN) :: delpc
! pkc is pe**cappa     if hydrostatic
! pkc is full pressure if non-hydrostatic
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(IN) :: &
&   pkc, gz
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(IN) :: rdxc(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL, INTENT(IN) :: rdyc(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: hydrostatic
! Local:
    REAL :: wk(bd%is-1:bd%ie+1, bd%js-1:bd%je+1)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
!$OMP parallel do default(none) shared(is,ie,js,je,npz,hydrostatic,pkc,delpc,uc,dt2,rdxc,gz,vc,rdyc) &
!$OMP                          private(wk)
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js-1,je+1
          DO i=is-1,ie+1
            wk(i, j) = pkc(i, j, k+1) - pkc(i, j, k)
          END DO
        END DO
      ELSE
        DO j=js-1,je+1
          DO i=is-1,ie+1
            wk(i, j) = delpc(i, j, k)
          END DO
        END DO
      END IF
      DO j=js,je
        DO i=is,ie+1
          uc(i, j, k) = uc(i, j, k) + dt2*rdxc(i, j)/(wk(i-1, j)+wk(i, j&
&           ))*((gz(i-1, j, k+1)-gz(i, j, k))*(pkc(i, j, k+1)-pkc(i-1, j&
&           , k))+(gz(i-1, j, k)-gz(i, j, k+1))*(pkc(i-1, j, k+1)-pkc(i&
&           , j, k)))
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          vc(i, j, k) = vc(i, j, k) + dt2*rdyc(i, j)/(wk(i, j-1)+wk(i, j&
&           ))*((gz(i, j-1, k+1)-gz(i, j, k))*(pkc(i, j, k+1)-pkc(i, j-1&
&           , k))+(gz(i, j-1, k)-gz(i, j, k+1))*(pkc(i, j-1, k+1)-pkc(i&
&           , j, k)))
        END DO
      END DO
    END DO
  END SUBROUTINE P_GRAD_C
!  Differentiation of nh_p_grad in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
!d.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_m
!od.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mi
!x_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_
!Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4
! fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.rem
!ap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv
!_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters f
!v_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_re
!start_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid
!_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_m
!od.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod
!.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.
!nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a
!2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v_f
!b sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_m
!od.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_
!mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: gz u v delp pk pp
!   with respect to varying inputs: gz u v delp pk pp
  SUBROUTINE NH_P_GRAD_FWD(u, v, pp, gz, delp, pk, dt, ng, gridstruct, &
&   bd, npx, npy, npz, use_logp)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    LOGICAL, INTENT(IN) :: use_logp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! perturbation pressure
    REAL, INTENT(INOUT) :: pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! p**kappa
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! g * h
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL :: wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: du1, dv1, top_value
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed

    wk1 = 0.0
    wk = 0.0
    du1 = 0.0
    dv1 = 0.0
    top_value = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (use_logp) THEN
      top_value = peln1
    ELSE
      top_value = ptk
    END IF
!Remember that not all compilers set pp to zero by default
!$OMP parallel do default(none) shared(is,ie,js,je,pp,pk,top_value)
    DO j=js,je+1
      DO i=is,ie+1
        CALL PUSHREALARRAY(pp(i, j, 1))
        pp(i, j, 1) = 0.
        CALL PUSHREALARRAY(pk(i, j, 1))
        pk(i, j, 1) = top_value
      END DO
    END DO
!$OMP parallel do default(none) shared(isd,jsd,npz,pp,gridstruct,npx,npy,is,ie,js,je,ng,pk,gz) &
!$OMP                          private(wk1)
    DO k=1,npz+1
      IF (k .NE. 1) THEN
        CALL A2B_ORD4_FWD(pp(isd:ied, jsd:jed, k), wk1, gridstruct, &
&                      npx, npy, is, ie, js, je, ng, .true.)
        CALL A2B_ORD4_FWD(pk(isd:ied, jsd:jed, k), wk1, gridstruct, &
&                      npx, npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      CALL A2B_ORD4_FWD(gz(isd:ied, jsd:jed, k), wk1, gridstruct, npx&
&                    , npy, is, ie, js, je, ng, .true.)
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,delp,gridstruct,npx,npy,ng,isd,jsd, &
!$OMP                                  pk,dt,gz,u,pp,v) &
!$OMP                          private(wk1, wk, du1, dv1)
    DO k=1,npz
      CALL A2B_ORD4_FWD(delp(isd:ied, jsd:jed, k), wk1, gridstruct, &
&                    npx, npy, is, ie, js, je, ng)
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY(wk(i, j))
          wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
! hydrostatic contributions from past time-step already added in the "beta" part
! Current gradient from "hydrostatic" components:
          du1 = dt/(wk(i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*&
&           (pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*&
&           (pk(i, j, k+1)-pk(i+1, j, k)))
! Non-hydrostatic contribution
          CALL PUSHREALARRAY(u(i, j, k))
          u(i, j, k) = (u(i, j, k)+du1+dt/(wk1(i, j)+wk1(i+1, j))*((gz(i&
&           , j, k+1)-gz(i+1, j, k))*(pp(i+1, j, k+1)-pp(i, j, k))+(gz(i&
&           , j, k)-gz(i+1, j, k+1))*(pp(i, j, k+1)-pp(i+1, j, k))))*&
&           gridstruct%rdx(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
! Current gradient from "hydrostatic" components:
          dv1 = dt/(wk(i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*&
&           (pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*&
&           (pk(i, j, k+1)-pk(i, j+1, k)))
! Non-hydrostatic contribution
          CALL PUSHREALARRAY(v(i, j, k))
          v(i, j, k) = (v(i, j, k)+dv1+dt/(wk1(i, j)+wk1(i, j+1))*((gz(i&
&           , j, k+1)-gz(i, j+1, k))*(pp(i, j+1, k+1)-pp(i, j, k))+(gz(i&
&           , j, k)-gz(i, j+1, k+1))*(pp(i, j, k+1)-pp(i, j+1, k))))*&
&           gridstruct%rdy(i, j)
        END DO
      END DO
    END DO
    CALL PUSHINTEGER(jed)
    CALL PUSHREALARRAY(wk, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(isd)
    CALL PUSHINTEGER(ie)
    CALL PUSHINTEGER(ied)
    CALL PUSHINTEGER(jsd)
    CALL PUSHREALARRAY(wk1, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL PUSHINTEGER(js)
  END SUBROUTINE NH_P_GRAD_FWD
!  Differentiation of nh_p_grad in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_m
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
!   gradient     of useful results: gz u v delp pk pp
!   with respect to varying inputs: gz u v delp pk pp
  SUBROUTINE NH_P_GRAD_BWD(u, u_ad, v, v_ad, pp, pp_ad, gz, gz_ad, delp&
&   , delp_ad, pk, pk_ad, dt, ng, gridstruct, bd, npx, npy, npz, &
&   use_logp)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    LOGICAL, INTENT(IN) :: use_logp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pp_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    REAL :: wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk1_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk_ad(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: du1, dv1, top_value
    REAL :: du1_ad, dv1_ad
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp2
    REAL :: temp3
    REAL :: temp4
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp8
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp9
    REAL :: temp10
    REAL :: temp11
    REAL :: temp12
    REAL :: temp13
    REAL :: temp14
    REAL :: temp15
    REAL :: temp16
    REAL :: temp17
    REAL :: temp18
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    INTEGER :: branch

    wk1 = 0.0
    wk = 0.0
    du1 = 0.0
    dv1 = 0.0
    top_value = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    branch = 0

    CALL POPINTEGER(js)
    CALL POPREALARRAY(wk1, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL POPINTEGER(jsd)
    CALL POPINTEGER(ied)
    CALL POPINTEGER(ie)
    CALL POPINTEGER(isd)
    CALL POPINTEGER(is)
    CALL POPINTEGER(je)
    CALL POPREALARRAY(wk, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
    CALL POPINTEGER(jed)
    wk1_ad = 0.0
    wk_ad = 0.0
    DO k=npz,1,-1
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(v(i, j, k))
          temp14 = wk1(i, j) + wk1(i, j+1)
          temp18 = pp(i, j, k+1) - pp(i, j+1, k)
          temp17 = gz(i, j, k) - gz(i, j+1, k+1)
          temp16 = pp(i, j+1, k+1) - pp(i, j, k)
          temp15 = gz(i, j, k+1) - gz(i, j+1, k)
          temp_ad4 = gridstruct%rdy(i, j)*v_ad(i, j, k)
          temp_ad5 = dt*temp_ad4/temp14
          temp_ad6 = -((temp15*temp16+temp17*temp18)*temp_ad5/temp14)
          dv1_ad = temp_ad4
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp16*temp_ad5
          gz_ad(i, j+1, k) = gz_ad(i, j+1, k) - temp16*temp_ad5
          pp_ad(i, j+1, k+1) = pp_ad(i, j+1, k+1) + temp15*temp_ad5
          pp_ad(i, j, k) = pp_ad(i, j, k) - temp15*temp_ad5
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp18*temp_ad5
          gz_ad(i, j+1, k+1) = gz_ad(i, j+1, k+1) - temp18*temp_ad5
          pp_ad(i, j, k+1) = pp_ad(i, j, k+1) + temp17*temp_ad5
          pp_ad(i, j+1, k) = pp_ad(i, j+1, k) - temp17*temp_ad5
          wk1_ad(i, j) = wk1_ad(i, j) + temp_ad6
          wk1_ad(i, j+1) = wk1_ad(i, j+1) + temp_ad6
          v_ad(i, j, k) = temp_ad4
          temp9 = wk(i, j) + wk(i, j+1)
          temp13 = pk(i, j, k+1) - pk(i, j+1, k)
          temp12 = gz(i, j, k) - gz(i, j+1, k+1)
          temp11 = pk(i, j+1, k+1) - pk(i, j, k)
          temp10 = gz(i, j, k+1) - gz(i, j+1, k)
          temp_ad7 = dt*dv1_ad/temp9
          temp_ad8 = -((temp10*temp11+temp12*temp13)*temp_ad7/temp9)
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp11*temp_ad7
          gz_ad(i, j+1, k) = gz_ad(i, j+1, k) - temp11*temp_ad7
          pk_ad(i, j+1, k+1) = pk_ad(i, j+1, k+1) + temp10*temp_ad7
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp10*temp_ad7
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp13*temp_ad7
          gz_ad(i, j+1, k+1) = gz_ad(i, j+1, k+1) - temp13*temp_ad7
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp12*temp_ad7
          pk_ad(i, j+1, k) = pk_ad(i, j+1, k) - temp12*temp_ad7
          wk_ad(i, j) = wk_ad(i, j) + temp_ad8
          wk_ad(i, j+1) = wk_ad(i, j+1) + temp_ad8
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(u(i, j, k))
          temp4 = wk1(i, j) + wk1(i+1, j)
          temp8 = pp(i, j, k+1) - pp(i+1, j, k)
          temp7 = gz(i, j, k) - gz(i+1, j, k+1)
          temp6 = pp(i+1, j, k+1) - pp(i, j, k)
          temp5 = gz(i, j, k+1) - gz(i+1, j, k)
          temp_ad = gridstruct%rdx(i, j)*u_ad(i, j, k)
          temp_ad0 = dt*temp_ad/temp4
          temp_ad1 = -((temp5*temp6+temp7*temp8)*temp_ad0/temp4)
          du1_ad = temp_ad
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp6*temp_ad0
          gz_ad(i+1, j, k) = gz_ad(i+1, j, k) - temp6*temp_ad0
          pp_ad(i+1, j, k+1) = pp_ad(i+1, j, k+1) + temp5*temp_ad0
          pp_ad(i, j, k) = pp_ad(i, j, k) - temp5*temp_ad0
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp8*temp_ad0
          gz_ad(i+1, j, k+1) = gz_ad(i+1, j, k+1) - temp8*temp_ad0
          pp_ad(i, j, k+1) = pp_ad(i, j, k+1) + temp7*temp_ad0
          pp_ad(i+1, j, k) = pp_ad(i+1, j, k) - temp7*temp_ad0
          wk1_ad(i, j) = wk1_ad(i, j) + temp_ad1
          wk1_ad(i+1, j) = wk1_ad(i+1, j) + temp_ad1
          u_ad(i, j, k) = temp_ad
          temp = wk(i, j) + wk(i+1, j)
          temp3 = pk(i, j, k+1) - pk(i+1, j, k)
          temp2 = gz(i, j, k) - gz(i+1, j, k+1)
          temp1 = pk(i+1, j, k+1) - pk(i, j, k)
          temp0 = gz(i, j, k+1) - gz(i+1, j, k)
          temp_ad2 = dt*du1_ad/temp
          temp_ad3 = -((temp0*temp1+temp2*temp3)*temp_ad2/temp)
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp1*temp_ad2
          gz_ad(i+1, j, k) = gz_ad(i+1, j, k) - temp1*temp_ad2
          pk_ad(i+1, j, k+1) = pk_ad(i+1, j, k+1) + temp0*temp_ad2
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp0*temp_ad2
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp3*temp_ad2
          gz_ad(i+1, j, k+1) = gz_ad(i+1, j, k+1) - temp3*temp_ad2
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp2*temp_ad2
          pk_ad(i+1, j, k) = pk_ad(i+1, j, k) - temp2*temp_ad2
          wk_ad(i, j) = wk_ad(i, j) + temp_ad3
          wk_ad(i+1, j) = wk_ad(i+1, j) + temp_ad3
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(wk(i, j))
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + wk_ad(i, j)
          pk_ad(i, j, k) = pk_ad(i, j, k) - wk_ad(i, j)
          wk_ad(i, j) = 0.0
        END DO
      END DO
      CALL A2B_ORD4_BWD(delp(isd:ied, jsd:jed, k), delp_ad(isd:ied, &
&                    jsd:jed, k), wk1, wk1_ad, gridstruct, npx, npy, is&
&                    , ie, js, je, ng)
    END DO
    DO k=npz+1,1,-1
      CALL A2B_ORD4_BWD(gz(isd:ied, jsd:jed, k), gz_ad(isd:ied, jsd:&
&                    jed, k), wk1, wk1_ad, gridstruct, npx, npy, is, ie&
&                    , js, je, ng, .true.)
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD4_BWD(pk(isd:ied, jsd:jed, k), pk_ad(isd:ied, jsd&
&                      :jed, k), wk1, wk1_ad, gridstruct, npx, npy, is, &
&                      ie, js, je, ng, .true.)
        CALL A2B_ORD4_BWD(pp(isd:ied, jsd:jed, k), pp_ad(isd:ied, jsd&
&                      :jed, k), wk1, wk1_ad, gridstruct, npx, npy, is, &
&                      ie, js, je, ng, .true.)
      END IF
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(pk(i, j, 1))
        pk_ad(i, j, 1) = 0.0
        CALL POPREALARRAY(pp(i, j, 1))
        pp_ad(i, j, 1) = 0.0
      END DO
    END DO
  END SUBROUTINE NH_P_GRAD_BWD
  SUBROUTINE NH_P_GRAD(u, v, pp, gz, delp, pk, dt, ng, gridstruct, bd, &
&   npx, npy, npz, use_logp)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    LOGICAL, INTENT(IN) :: use_logp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! perturbation pressure
    REAL, INTENT(INOUT) :: pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! p**kappa
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! g * h
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL :: wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: du1, dv1, top_value
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (use_logp) THEN
      top_value = peln1
    ELSE
      top_value = ptk
    END IF
!Remember that not all compilers set pp to zero by default
!$OMP parallel do default(none) shared(is,ie,js,je,pp,pk,top_value)
    DO j=js,je+1
      DO i=is,ie+1
        pp(i, j, 1) = 0.
        pk(i, j, 1) = top_value
      END DO
    END DO
!$OMP parallel do default(none) shared(isd,jsd,npz,pp,gridstruct,npx,npy,is,ie,js,je,ng,pk,gz) &
!$OMP                          private(wk1)
    DO k=1,npz+1
      IF (k .NE. 1) THEN
        CALL A2B_ORD4(pp(isd:ied, jsd:jed, k), wk1, gridstruct, npx, &
&                  npy, is, ie, js, je, ng, .true.)
        CALL A2B_ORD4(pk(isd:ied, jsd:jed, k), wk1, gridstruct, npx, &
&                  npy, is, ie, js, je, ng, .true.)
      END IF
      CALL A2B_ORD4(gz(isd:ied, jsd:jed, k), wk1, gridstruct, npx, &
&                npy, is, ie, js, je, ng, .true.)
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,delp,gridstruct,npx,npy,ng,isd,jsd, &
!$OMP                                  pk,dt,gz,u,pp,v) &
!$OMP                          private(wk1, wk, du1, dv1)
    DO k=1,npz
      CALL A2B_ORD4(delp(isd:ied, jsd:jed, k), wk1, gridstruct, npx, &
&                npy, is, ie, js, je, ng)
      DO j=js,je+1
        DO i=is,ie+1
          wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
! hydrostatic contributions from past time-step already added in the "beta" part
! Current gradient from "hydrostatic" components:
          du1 = dt/(wk(i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*&
&           (pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*&
&           (pk(i, j, k+1)-pk(i+1, j, k)))
! Non-hydrostatic contribution
          u(i, j, k) = (u(i, j, k)+du1+dt/(wk1(i, j)+wk1(i+1, j))*((gz(i&
&           , j, k+1)-gz(i+1, j, k))*(pp(i+1, j, k+1)-pp(i, j, k))+(gz(i&
&           , j, k)-gz(i+1, j, k+1))*(pp(i, j, k+1)-pp(i+1, j, k))))*&
&           gridstruct%rdx(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
! Current gradient from "hydrostatic" components:
          dv1 = dt/(wk(i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*&
&           (pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*&
&           (pk(i, j, k+1)-pk(i, j+1, k)))
! Non-hydrostatic contribution
          v(i, j, k) = (v(i, j, k)+dv1+dt/(wk1(i, j)+wk1(i, j+1))*((gz(i&
&           , j, k+1)-gz(i, j+1, k))*(pp(i, j+1, k+1)-pp(i, j, k))+(gz(i&
&           , j, k)-gz(i, j+1, k+1))*(pp(i, j, k+1)-pp(i, j+1, k))))*&
&           gridstruct%rdy(i, j)
        END DO
      END DO
    END DO
  END SUBROUTINE NH_P_GRAD
!  Differentiation of split_p_grad in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: gz du u dv v delp pk pp
!   with respect to varying inputs: gz du u dv v delp pk pp
  SUBROUTINE SPLIT_P_GRAD_FWD(u, v, pp, gz, du, dv, delp, pk, beta, dt, &
&   ng, gridstruct, bd, npx, npy, npz, use_logp)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: beta, dt
    LOGICAL, INTENT(IN) :: use_logp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! perturbation pressure
    REAL, INTENT(INOUT) :: pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! p**kappa
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! g * h
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL :: wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: alpha, top_value
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed

    wk1 = 0.0
    wk = 0.0
    alpha = 0.0
    top_value = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (use_logp) THEN
      top_value = peln1
    ELSE
      top_value = ptk
    END IF
    alpha = 1. - beta
!$OMP parallel do default(none) shared(is,ie,js,je,pp,pk,top_value)
    DO j=js,je+1
      DO i=is,ie+1
        CALL PUSHREALARRAY(pp(i, j, 1))
        pp(i, j, 1) = 0.
        CALL PUSHREALARRAY(pk(i, j, 1))
        pk(i, j, 1) = top_value
      END DO
    END DO
!$OMP parallel do default(none) shared(isd,jsd,npz,pp,gridstruct,npx,npy,is,ie,js,je,ng,pk,gz) &
!$OMP                          private(wk1)
    DO k=1,npz+1
      IF (k .NE. 1) THEN
        CALL A2B_ORD4_FWD(pp(isd:ied, jsd:jed, k), wk1, gridstruct, &
&                      npx, npy, is, ie, js, je, ng, .true.)
        CALL A2B_ORD4_FWD(pk(isd:ied, jsd:jed, k), wk1, gridstruct, &
&                      npx, npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      CALL A2B_ORD4_FWD(gz(isd:ied, jsd:jed, k), wk1, gridstruct, npx&
&                    , npy, is, ie, js, je, ng, .true.)
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,isd,jsd,npz,delp,gridstruct,npx,npy,ng, &
!$OMP                                  pk,u,beta,du,dt,gz,alpha,pp,v,dv) &
!$OMP                          private(wk1, wk)
    DO k=1,npz
      CALL A2B_ORD4_FWD(delp(isd:ied, jsd:jed, k), wk1, gridstruct, &
&                    npx, npy, is, ie, js, je, ng)
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY(wk(i, j))
          wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(u(i, j, k))
          u(i, j, k) = u(i, j, k) + beta*du(i, j, k)
! hydrostatic contributions from past time-step already added in the "beta" part
! Current gradient from "hydrostatic" components:
!---------------------------------------------------------------------------------
          du(i, j, k) = dt/(wk(i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1&
&           , j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, &
&           j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
!---------------------------------------------------------------------------------
! Non-hydrostatic contribution
          u(i, j, k) = (u(i, j, k)+alpha*du(i, j, k)+dt/(wk1(i, j)+wk1(i&
&           +1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*(pp(i+1, j, k+1)-pp(i&
&           , j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*(pp(i, j, k+1)-pp(i+1&
&           , j, k))))*gridstruct%rdx(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(v(i, j, k))
          v(i, j, k) = v(i, j, k) + beta*dv(i, j, k)
! Current gradient from "hydrostatic" components:
!---------------------------------------------------------------------------------
          dv(i, j, k) = dt/(wk(i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j&
&           +1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1&
&           , k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
!---------------------------------------------------------------------------------
! Non-hydrostatic contribution
          v(i, j, k) = (v(i, j, k)+alpha*dv(i, j, k)+dt/(wk1(i, j)+wk1(i&
&           , j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*(pp(i, j+1, k+1)-pp(i&
&           , j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*(pp(i, j, k+1)-pp(i, &
&           j+1, k))))*gridstruct%rdy(i, j)
        END DO
      END DO
    END DO
    CALL PUSHINTEGER(jed)
    CALL PUSHREALARRAY(wk, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(isd)
    CALL PUSHINTEGER(ie)
    CALL PUSHINTEGER(ied)
    CALL PUSHINTEGER(jsd)
    CALL PUSHREALARRAY(wk1, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL PUSHREALARRAY(alpha)
    CALL PUSHINTEGER(js)
  END SUBROUTINE SPLIT_P_GRAD_FWD
!  Differentiation of split_p_grad in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: gz du u dv v delp pk pp
!   with respect to varying inputs: gz du u dv v delp pk pp
  SUBROUTINE SPLIT_P_GRAD_BWD(u, u_ad, v, v_ad, pp, pp_ad, gz, gz_ad, du&
&   , du_ad, dv, dv_ad, delp, delp_ad, pk, pk_ad, beta, dt, ng, &
&   gridstruct, bd, npx, npy, npz, use_logp)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: beta, dt
    LOGICAL, INTENT(IN) :: use_logp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pp_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: du_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dv_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    REAL :: wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk1_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk_ad(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: alpha, top_value
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp2
    REAL :: temp3
    REAL :: temp4
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp8
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp9
    REAL :: temp10
    REAL :: temp11
    REAL :: temp12
    REAL :: temp13
    REAL :: temp14
    REAL :: temp15
    REAL :: temp16
    REAL :: temp17
    REAL :: temp18
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    INTEGER :: branch

    wk1 = 0.0
    wk = 0.0
    alpha = 0.0
    top_value = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    branch = 0

    CALL POPINTEGER(js)
    CALL POPREALARRAY(alpha)
    CALL POPREALARRAY(wk1, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL POPINTEGER(jsd)
    CALL POPINTEGER(ied)
    CALL POPINTEGER(ie)
    CALL POPINTEGER(isd)
    CALL POPINTEGER(is)
    CALL POPINTEGER(je)
    CALL POPREALARRAY(wk, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
    CALL POPINTEGER(jed)
    wk1_ad = 0.0
    wk_ad = 0.0
    DO k=npz,1,-1
      DO j=je,js,-1
        DO i=ie+1,is,-1
          temp14 = wk1(i, j) + wk1(i, j+1)
          temp18 = pp(i, j, k+1) - pp(i, j+1, k)
          temp17 = gz(i, j, k) - gz(i, j+1, k+1)
          temp16 = pp(i, j+1, k+1) - pp(i, j, k)
          temp15 = gz(i, j, k+1) - gz(i, j+1, k)
          temp_ad4 = gridstruct%rdy(i, j)*v_ad(i, j, k)
          temp_ad5 = dt*temp_ad4/temp14
          temp_ad6 = -((temp15*temp16+temp17*temp18)*temp_ad5/temp14)
          dv_ad(i, j, k) = dv_ad(i, j, k) + alpha*temp_ad4
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp16*temp_ad5
          gz_ad(i, j+1, k) = gz_ad(i, j+1, k) - temp16*temp_ad5
          pp_ad(i, j+1, k+1) = pp_ad(i, j+1, k+1) + temp15*temp_ad5
          pp_ad(i, j, k) = pp_ad(i, j, k) - temp15*temp_ad5
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp18*temp_ad5
          gz_ad(i, j+1, k+1) = gz_ad(i, j+1, k+1) - temp18*temp_ad5
          pp_ad(i, j, k+1) = pp_ad(i, j, k+1) + temp17*temp_ad5
          pp_ad(i, j+1, k) = pp_ad(i, j+1, k) - temp17*temp_ad5
          wk1_ad(i, j) = wk1_ad(i, j) + temp_ad6
          wk1_ad(i, j+1) = wk1_ad(i, j+1) + temp_ad6
          v_ad(i, j, k) = temp_ad4
          temp9 = wk(i, j) + wk(i, j+1)
          temp13 = pk(i, j, k+1) - pk(i, j+1, k)
          temp12 = gz(i, j, k) - gz(i, j+1, k+1)
          temp11 = pk(i, j+1, k+1) - pk(i, j, k)
          temp10 = gz(i, j, k+1) - gz(i, j+1, k)
          temp_ad7 = dt*dv_ad(i, j, k)/temp9
          temp_ad8 = -((temp10*temp11+temp12*temp13)*temp_ad7/temp9)
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp11*temp_ad7
          gz_ad(i, j+1, k) = gz_ad(i, j+1, k) - temp11*temp_ad7
          pk_ad(i, j+1, k+1) = pk_ad(i, j+1, k+1) + temp10*temp_ad7
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp10*temp_ad7
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp13*temp_ad7
          gz_ad(i, j+1, k+1) = gz_ad(i, j+1, k+1) - temp13*temp_ad7
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp12*temp_ad7
          pk_ad(i, j+1, k) = pk_ad(i, j+1, k) - temp12*temp_ad7
          wk_ad(i, j) = wk_ad(i, j) + temp_ad8
          wk_ad(i, j+1) = wk_ad(i, j+1) + temp_ad8
          dv_ad(i, j, k) = beta*v_ad(i, j, k)
          CALL POPREALARRAY(v(i, j, k))
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          temp4 = wk1(i, j) + wk1(i+1, j)
          temp8 = pp(i, j, k+1) - pp(i+1, j, k)
          temp7 = gz(i, j, k) - gz(i+1, j, k+1)
          temp6 = pp(i+1, j, k+1) - pp(i, j, k)
          temp5 = gz(i, j, k+1) - gz(i+1, j, k)
          temp_ad = gridstruct%rdx(i, j)*u_ad(i, j, k)
          temp_ad0 = dt*temp_ad/temp4
          temp_ad1 = -((temp5*temp6+temp7*temp8)*temp_ad0/temp4)
          du_ad(i, j, k) = du_ad(i, j, k) + alpha*temp_ad
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp6*temp_ad0
          gz_ad(i+1, j, k) = gz_ad(i+1, j, k) - temp6*temp_ad0
          pp_ad(i+1, j, k+1) = pp_ad(i+1, j, k+1) + temp5*temp_ad0
          pp_ad(i, j, k) = pp_ad(i, j, k) - temp5*temp_ad0
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp8*temp_ad0
          gz_ad(i+1, j, k+1) = gz_ad(i+1, j, k+1) - temp8*temp_ad0
          pp_ad(i, j, k+1) = pp_ad(i, j, k+1) + temp7*temp_ad0
          pp_ad(i+1, j, k) = pp_ad(i+1, j, k) - temp7*temp_ad0
          wk1_ad(i, j) = wk1_ad(i, j) + temp_ad1
          wk1_ad(i+1, j) = wk1_ad(i+1, j) + temp_ad1
          u_ad(i, j, k) = temp_ad
          temp = wk(i, j) + wk(i+1, j)
          temp3 = pk(i, j, k+1) - pk(i+1, j, k)
          temp2 = gz(i, j, k) - gz(i+1, j, k+1)
          temp1 = pk(i+1, j, k+1) - pk(i, j, k)
          temp0 = gz(i, j, k+1) - gz(i+1, j, k)
          temp_ad2 = dt*du_ad(i, j, k)/temp
          temp_ad3 = -((temp0*temp1+temp2*temp3)*temp_ad2/temp)
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp1*temp_ad2
          gz_ad(i+1, j, k) = gz_ad(i+1, j, k) - temp1*temp_ad2
          pk_ad(i+1, j, k+1) = pk_ad(i+1, j, k+1) + temp0*temp_ad2
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp0*temp_ad2
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp3*temp_ad2
          gz_ad(i+1, j, k+1) = gz_ad(i+1, j, k+1) - temp3*temp_ad2
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp2*temp_ad2
          pk_ad(i+1, j, k) = pk_ad(i+1, j, k) - temp2*temp_ad2
          wk_ad(i, j) = wk_ad(i, j) + temp_ad3
          wk_ad(i+1, j) = wk_ad(i+1, j) + temp_ad3
          du_ad(i, j, k) = beta*u_ad(i, j, k)
          CALL POPREALARRAY(u(i, j, k))
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(wk(i, j))
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + wk_ad(i, j)
          pk_ad(i, j, k) = pk_ad(i, j, k) - wk_ad(i, j)
          wk_ad(i, j) = 0.0
        END DO
      END DO
      CALL A2B_ORD4_BWD(delp(isd:ied, jsd:jed, k), delp_ad(isd:ied, &
&                    jsd:jed, k), wk1, wk1_ad, gridstruct, npx, npy, is&
&                    , ie, js, je, ng)
    END DO
    DO k=npz+1,1,-1
      CALL A2B_ORD4_BWD(gz(isd:ied, jsd:jed, k), gz_ad(isd:ied, jsd:&
&                    jed, k), wk1, wk1_ad, gridstruct, npx, npy, is, ie&
&                    , js, je, ng, .true.)
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD4_BWD(pk(isd:ied, jsd:jed, k), pk_ad(isd:ied, jsd&
&                      :jed, k), wk1, wk1_ad, gridstruct, npx, npy, is, &
&                      ie, js, je, ng, .true.)
        CALL A2B_ORD4_BWD(pp(isd:ied, jsd:jed, k), pp_ad(isd:ied, jsd&
&                      :jed, k), wk1, wk1_ad, gridstruct, npx, npy, is, &
&                      ie, js, je, ng, .true.)
      END IF
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(pk(i, j, 1))
        pk_ad(i, j, 1) = 0.0
        CALL POPREALARRAY(pp(i, j, 1))
        pp_ad(i, j, 1) = 0.0
      END DO
    END DO
  END SUBROUTINE SPLIT_P_GRAD_BWD
  SUBROUTINE SPLIT_P_GRAD(u, v, pp, gz, du, dv, delp, pk, beta, dt, ng, &
&   gridstruct, bd, npx, npy, npz, use_logp)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: beta, dt
    LOGICAL, INTENT(IN) :: use_logp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! perturbation pressure
    REAL, INTENT(INOUT) :: pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! p**kappa
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! g * h
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL :: wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: alpha, top_value
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (use_logp) THEN
      top_value = peln1
    ELSE
      top_value = ptk
    END IF
    alpha = 1. - beta
!$OMP parallel do default(none) shared(is,ie,js,je,pp,pk,top_value)
    DO j=js,je+1
      DO i=is,ie+1
        pp(i, j, 1) = 0.
        pk(i, j, 1) = top_value
      END DO
    END DO
!$OMP parallel do default(none) shared(isd,jsd,npz,pp,gridstruct,npx,npy,is,ie,js,je,ng,pk,gz) &
!$OMP                          private(wk1)
    DO k=1,npz+1
      IF (k .NE. 1) THEN
        CALL A2B_ORD4(pp(isd:ied, jsd:jed, k), wk1, gridstruct, npx, &
&                  npy, is, ie, js, je, ng, .true.)
        CALL A2B_ORD4(pk(isd:ied, jsd:jed, k), wk1, gridstruct, npx, &
&                  npy, is, ie, js, je, ng, .true.)
      END IF
      CALL A2B_ORD4(gz(isd:ied, jsd:jed, k), wk1, gridstruct, npx, &
&                npy, is, ie, js, je, ng, .true.)
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,isd,jsd,npz,delp,gridstruct,npx,npy,ng, &
!$OMP                                  pk,u,beta,du,dt,gz,alpha,pp,v,dv) &
!$OMP                          private(wk1, wk)
    DO k=1,npz
      CALL A2B_ORD4(delp(isd:ied, jsd:jed, k), wk1, gridstruct, npx, &
&                npy, is, ie, js, je, ng)
      DO j=js,je+1
        DO i=is,ie+1
          wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          u(i, j, k) = u(i, j, k) + beta*du(i, j, k)
! hydrostatic contributions from past time-step already added in the "beta" part
! Current gradient from "hydrostatic" components:
!---------------------------------------------------------------------------------
          du(i, j, k) = dt/(wk(i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1&
&           , j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, &
&           j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
!---------------------------------------------------------------------------------
! Non-hydrostatic contribution
          u(i, j, k) = (u(i, j, k)+alpha*du(i, j, k)+dt/(wk1(i, j)+wk1(i&
&           +1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*(pp(i+1, j, k+1)-pp(i&
&           , j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*(pp(i, j, k+1)-pp(i+1&
&           , j, k))))*gridstruct%rdx(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v(i, j, k) = v(i, j, k) + beta*dv(i, j, k)
! Current gradient from "hydrostatic" components:
!---------------------------------------------------------------------------------
          dv(i, j, k) = dt/(wk(i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j&
&           +1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1&
&           , k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
!---------------------------------------------------------------------------------
! Non-hydrostatic contribution
          v(i, j, k) = (v(i, j, k)+alpha*dv(i, j, k)+dt/(wk1(i, j)+wk1(i&
&           , j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*(pp(i, j+1, k+1)-pp(i&
&           , j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*(pp(i, j, k+1)-pp(i, &
&           j+1, k))))*gridstruct%rdy(i, j)
        END DO
      END DO
    END DO
  END SUBROUTINE SPLIT_P_GRAD
!  Differentiation of one_grad_p in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_m
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
!   gradient     of useful results: gz u v delp pk divg2
!   with respect to varying inputs: gz u v delp pk divg2
  SUBROUTINE ONE_GRAD_P_FWD(u, v, pk, gz, divg2, delp, dt, ng, &
&   gridstruct, bd, npx, npy, npz, ptop, hydrostatic, a2b_ord, d_ext)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz, a2b_ord
    REAL, INTENT(IN) :: dt, ptop, d_ext
    LOGICAL, INTENT(IN) :: hydrostatic
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: wk
    REAL :: wk1(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk2(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: top_value
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed

    wk = 0.0
    wk1 = 0.0
    wk2 = 0.0
    top_value = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (hydrostatic) THEN
! pk is pe**kappa if hydrostatic
      top_value = ptk
    ELSE
! pk is full pressure if non-hydrostatic
      top_value = ptop
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,pk,top_value)
    DO j=js,je+1
      DO i=is,ie+1
        CALL PUSHREALARRAY(pk(i, j, 1))
        pk(i, j, 1) = top_value
      END DO
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,pk,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_FWD(pk(isd:ied, jsd:jed, k), wk, gridstruct, &
&                      npx, npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL A2B_ORD2_FWD(pk(isd:ied, jsd:jed, k), wk, gridstruct, npx, &
&                   npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,gz,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_FWD(gz(isd:ied, jsd:jed, k), wk, gridstruct, &
&                      npx, npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL A2B_ORD2_FWD(gz(isd:ied, jsd:jed, k), wk, gridstruct, npx, &
&                   npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
    IF (d_ext .GT. 0.) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,wk2,divg2)
      DO j=js,je+1
        DO i=is,ie
          wk2(i, j) = divg2(i, j) - divg2(i+1, j)
        END DO
      END DO
!$OMP parallel do default(none) shared(is,ie,js,je,wk1,divg2)
      DO j=js,je
        DO i=is,ie+1
          wk1(i, j) = divg2(i, j) - divg2(i, j+1)
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
!$OMP parallel do default(none) shared(is,ie,js,je,wk1,wk2)
      DO j=js,je+1
        DO i=is,ie
          wk2(i, j) = 0.
        END DO
        DO i=is,ie+1
          wk1(i, j) = 0.
        END DO
      END DO
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pk,delp,hydrostatic,a2b_ord,gridstruct, &
!$OMP                                  npx,npy,isd,jsd,ng,u,v,wk2,dt,gz,wk1) &
!$OMP                          private(wk)
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js,je+1
          DO i=is,ie+1
            CALL PUSHREALARRAY(wk(i, j))
            wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
          END DO
        END DO
        CALL PUSHCONTROL(2,2)
      ELSE IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_FWD(delp(isd:ied, jsd:jed, k), wk, gridstruct, &
&                      npx, npy, is, ie, js, je, ng)
        CALL PUSHCONTROL(2,1)
      ELSE
        CALL A2B_ORD2_FWD(delp(isd:ied, jsd:jed, k), wk, gridstruct, npx&
&                   , npy, is, ie, js, je, ng)
        CALL PUSHCONTROL(2,0)
      END IF
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(u(i, j, k))
          u(i, j, k) = gridstruct%rdx(i, j)*(wk2(i, j)+u(i, j, k)+dt/(wk&
&           (i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*(pk(i+1, j&
&           , k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*(pk(i, j, &
&           k+1)-pk(i+1, j, k))))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(v(i, j, k))
          v(i, j, k) = gridstruct%rdy(i, j)*(wk1(i, j)+v(i, j, k)+dt/(wk&
&           (i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*(pk(i, j+1&
&           , k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*(pk(i, j, &
&           k+1)-pk(i, j+1, k))))
        END DO
      END DO
    END DO
    CALL PUSHINTEGER(jed)
    CALL PUSHREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(isd)
    CALL PUSHINTEGER(ie)
    CALL PUSHINTEGER(ied)
    CALL PUSHINTEGER(jsd)
    CALL PUSHINTEGER(js)
  END SUBROUTINE ONE_GRAD_P_FWD
!  Differentiation of one_grad_p in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: gz u v delp pk divg2
!   with respect to varying inputs: gz u v delp pk divg2
  SUBROUTINE ONE_GRAD_P_BWD(u, u_ad, v, v_ad, pk, pk_ad, gz, gz_ad, &
&   divg2, divg2_ad, delp, delp_ad, dt, ng, gridstruct, bd, npx, npy, &
&   npz, ptop, hydrostatic, a2b_ord, d_ext)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz, a2b_ord
    REAL, INTENT(IN) :: dt, ptop, d_ext
    LOGICAL, INTENT(IN) :: hydrostatic
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: divg2_ad(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: wk
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: wk_ad
    REAL :: wk1(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk1_ad(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk2(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: wk2_ad(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: top_value
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp2
    REAL :: temp3
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp4
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp8
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    INTEGER :: branch

    wk = 0.0
    wk1 = 0.0
    wk2 = 0.0
    top_value = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    branch = 0

    CALL POPINTEGER(js)
    CALL POPINTEGER(jsd)
    CALL POPINTEGER(ied)
    CALL POPINTEGER(ie)
    CALL POPINTEGER(isd)
    CALL POPINTEGER(is)
    CALL POPINTEGER(je)
    CALL POPREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL POPINTEGER(jed)
    wk1_ad = 0.0
    wk2_ad = 0.0
    wk_ad = 0.0
    DO k=npz,1,-1
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(v(i, j, k))
          temp4 = wk(i, j) + wk(i, j+1)
          temp8 = pk(i, j, k+1) - pk(i, j+1, k)
          temp7 = gz(i, j, k) - gz(i, j+1, k+1)
          temp6 = pk(i, j+1, k+1) - pk(i, j, k)
          temp5 = gz(i, j, k+1) - gz(i, j+1, k)
          temp_ad2 = gridstruct%rdy(i, j)*v_ad(i, j, k)
          temp_ad3 = dt*temp_ad2/temp4
          temp_ad4 = -((temp5*temp6+temp7*temp8)*temp_ad3/temp4)
          wk1_ad(i, j) = wk1_ad(i, j) + temp_ad2
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp6*temp_ad3
          gz_ad(i, j+1, k) = gz_ad(i, j+1, k) - temp6*temp_ad3
          pk_ad(i, j+1, k+1) = pk_ad(i, j+1, k+1) + temp5*temp_ad3
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp5*temp_ad3
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp8*temp_ad3
          gz_ad(i, j+1, k+1) = gz_ad(i, j+1, k+1) - temp8*temp_ad3
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp7*temp_ad3
          pk_ad(i, j+1, k) = pk_ad(i, j+1, k) - temp7*temp_ad3
          wk_ad(i, j) = wk_ad(i, j) + temp_ad4
          wk_ad(i, j+1) = wk_ad(i, j+1) + temp_ad4
          v_ad(i, j, k) = temp_ad2
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(u(i, j, k))
          temp = wk(i, j) + wk(i+1, j)
          temp3 = pk(i, j, k+1) - pk(i+1, j, k)
          temp2 = gz(i, j, k) - gz(i+1, j, k+1)
          temp1 = pk(i+1, j, k+1) - pk(i, j, k)
          temp0 = gz(i, j, k+1) - gz(i+1, j, k)
          temp_ad = gridstruct%rdx(i, j)*u_ad(i, j, k)
          temp_ad0 = dt*temp_ad/temp
          temp_ad1 = -((temp0*temp1+temp2*temp3)*temp_ad0/temp)
          wk2_ad(i, j) = wk2_ad(i, j) + temp_ad
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp1*temp_ad0
          gz_ad(i+1, j, k) = gz_ad(i+1, j, k) - temp1*temp_ad0
          pk_ad(i+1, j, k+1) = pk_ad(i+1, j, k+1) + temp0*temp_ad0
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp0*temp_ad0
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp3*temp_ad0
          gz_ad(i+1, j, k+1) = gz_ad(i+1, j, k+1) - temp3*temp_ad0
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp2*temp_ad0
          pk_ad(i+1, j, k) = pk_ad(i+1, j, k) - temp2*temp_ad0
          wk_ad(i, j) = wk_ad(i, j) + temp_ad1
          wk_ad(i+1, j) = wk_ad(i+1, j) + temp_ad1
          u_ad(i, j, k) = temp_ad
        END DO
      END DO
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD2_BWD(delp(isd:ied, jsd:jed, k), delp_ad(isd:ied, &
&                   jsd:jed, k), wk, wk_ad, gridstruct, npx, npy, is, ie&
&                   , js, je, ng)
      ELSE IF (branch .EQ. 1) THEN
        CALL A2B_ORD4_BWD(delp(isd:ied, jsd:jed, k), delp_ad(isd:ied&
&                      , jsd:jed, k), wk, wk_ad, gridstruct, npx, npy, &
&                      is, ie, js, je, ng)
      ELSE
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPREALARRAY(wk(i, j))
            pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + wk_ad(i, j)
            pk_ad(i, j, k) = pk_ad(i, j, k) - wk_ad(i, j)
            wk_ad(i, j) = 0.0
          END DO
        END DO
      END IF
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      DO j=je,js,-1
        DO i=ie+1,is,-1
          divg2_ad(i, j) = divg2_ad(i, j) + wk1_ad(i, j)
          divg2_ad(i, j+1) = divg2_ad(i, j+1) - wk1_ad(i, j)
          wk1_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          divg2_ad(i, j) = divg2_ad(i, j) + wk2_ad(i, j)
          divg2_ad(i+1, j) = divg2_ad(i+1, j) - wk2_ad(i, j)
          wk2_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    DO k=npz+1,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD2_BWD(gz(isd:ied, jsd:jed, k), gz_ad(isd:ied, jsd:&
&                   jed, k), wk, wk_ad, gridstruct, npx, npy, is, ie, js&
&                   , je, ng, .true.)
      ELSE
        CALL A2B_ORD4_BWD(gz(isd:ied, jsd:jed, k), gz_ad(isd:ied, jsd&
&                      :jed, k), wk, wk_ad, gridstruct, npx, npy, is, ie&
&                      , js, je, ng, .true.)
      END IF
    END DO
    DO k=npz+1,2,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD2_BWD(pk(isd:ied, jsd:jed, k), pk_ad(isd:ied, jsd:&
&                   jed, k), wk, wk_ad, gridstruct, npx, npy, is, ie, js&
&                   , je, ng, .true.)
      ELSE
        CALL A2B_ORD4_BWD(pk(isd:ied, jsd:jed, k), pk_ad(isd:ied, jsd&
&                      :jed, k), wk, wk_ad, gridstruct, npx, npy, is, ie&
&                      , js, je, ng, .true.)
      END IF
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(pk(i, j, 1))
        pk_ad(i, j, 1) = 0.0
      END DO
    END DO
  END SUBROUTINE ONE_GRAD_P_BWD
  SUBROUTINE ONE_GRAD_P(u, v, pk, gz, divg2, delp, dt, ng, gridstruct, &
&   bd, npx, npy, npz, ptop, hydrostatic, a2b_ord, d_ext)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz, a2b_ord
    REAL, INTENT(IN) :: dt, ptop, d_ext
    LOGICAL, INTENT(IN) :: hydrostatic
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: wk
    REAL :: wk1(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk2(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: top_value
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (hydrostatic) THEN
! pk is pe**kappa if hydrostatic
      top_value = ptk
    ELSE
! pk is full pressure if non-hydrostatic
      top_value = ptop
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,pk,top_value)
    DO j=js,je+1
      DO i=is,ie+1
        pk(i, j, 1) = top_value
      END DO
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,pk,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(pk(isd:ied, jsd:jed, k), wk, gridstruct, npx, &
&                  npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2(pk(isd:ied, jsd:jed, k), wk, gridstruct, npx, npy&
&               , is, ie, js, je, ng, .true.)
      END IF
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,gz,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(gz(isd:ied, jsd:jed, k), wk, gridstruct, npx, &
&                  npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2(gz(isd:ied, jsd:jed, k), wk, gridstruct, npx, npy&
&               , is, ie, js, je, ng, .true.)
      END IF
    END DO
    IF (d_ext .GT. 0.) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,wk2,divg2)
      DO j=js,je+1
        DO i=is,ie
          wk2(i, j) = divg2(i, j) - divg2(i+1, j)
        END DO
      END DO
!$OMP parallel do default(none) shared(is,ie,js,je,wk1,divg2)
      DO j=js,je
        DO i=is,ie+1
          wk1(i, j) = divg2(i, j) - divg2(i, j+1)
        END DO
      END DO
    ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,wk1,wk2)
      DO j=js,je+1
        DO i=is,ie
          wk2(i, j) = 0.
        END DO
        DO i=is,ie+1
          wk1(i, j) = 0.
        END DO
      END DO
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pk,delp,hydrostatic,a2b_ord,gridstruct, &
!$OMP                                  npx,npy,isd,jsd,ng,u,v,wk2,dt,gz,wk1) &
!$OMP                          private(wk)
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js,je+1
          DO i=is,ie+1
            wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
          END DO
        END DO
      ELSE IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(delp(isd:ied, jsd:jed, k), wk, gridstruct, npx&
&                  , npy, is, ie, js, je, ng)
      ELSE
        CALL A2B_ORD2(delp(isd:ied, jsd:jed, k), wk, gridstruct, npx, &
&               npy, is, ie, js, je, ng)
      END IF
      DO j=js,je+1
        DO i=is,ie
          u(i, j, k) = gridstruct%rdx(i, j)*(wk2(i, j)+u(i, j, k)+dt/(wk&
&           (i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*(pk(i+1, j&
&           , k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*(pk(i, j, &
&           k+1)-pk(i+1, j, k))))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v(i, j, k) = gridstruct%rdy(i, j)*(wk1(i, j)+v(i, j, k)+dt/(wk&
&           (i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*(pk(i, j+1&
&           , k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*(pk(i, j, &
&           k+1)-pk(i, j+1, k))))
        END DO
      END DO
    END DO
  END SUBROUTINE ONE_GRAD_P
!  Differentiation of grad1_p_update in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
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
!   gradient     of useful results: gz du u dv v pk divg2
!   with respect to varying inputs: gz du u dv v pk divg2
  SUBROUTINE GRAD1_P_UPDATE_FWD(divg2, u, v, pk, gz, du, dv, dt, ng, &
&   gridstruct, bd, npx, npy, npz, ptop, beta, a2b_ord)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz, a2b_ord
    REAL, INTENT(IN) :: dt, ptop, beta
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: top_value, alpha
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed

    wk = 0.0
    top_value = 0.0
    alpha = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    alpha = 1. - beta
! pk is pe**kappa if hydrostatic
    top_value = ptk
!$OMP parallel do default(none) shared(is,ie,js,je,pk,top_value)
    DO j=js,je+1
      DO i=is,ie+1
        CALL PUSHREALARRAY(pk(i, j, 1))
        pk(i, j, 1) = top_value
      END DO
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,pk,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_FWD(pk(isd:ied, jsd:jed, k), wk, gridstruct, &
&                      npx, npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL A2B_ORD2_FWD(pk(isd:ied, jsd:jed, k), wk, gridstruct, npx, &
&                   npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,gz,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_FWD(gz(isd:ied, jsd:jed, k), wk, gridstruct, &
&                      npx, npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL A2B_ORD2_FWD(gz(isd:ied, jsd:jed, k), wk, gridstruct, npx, &
&                   npy, is, ie, js, je, ng, .true.)
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
!$OMP parallel do default(none) shared(npz,is,ie,js,je,pk,u,beta,gz,divg2,alpha, &
!$OMP                                  gridstruct,v,dt,du,dv) &
!$OMP                          private(wk)
    DO k=1,npz
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY(wk(i, j))
          wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(u(i, j, k))
          u(i, j, k) = u(i, j, k) + beta*du(i, j, k)
          du(i, j, k) = dt/(wk(i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1&
&           , j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, &
&           j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
          u(i, j, k) = (u(i, j, k)+divg2(i, j)-divg2(i+1, j)+alpha*du(i&
&           , j, k))*gridstruct%rdx(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(v(i, j, k))
          v(i, j, k) = v(i, j, k) + beta*dv(i, j, k)
          dv(i, j, k) = dt/(wk(i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j&
&           +1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1&
&           , k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
          v(i, j, k) = (v(i, j, k)+divg2(i, j)-divg2(i, j+1)+alpha*dv(i&
&           , j, k))*gridstruct%rdy(i, j)
        END DO
      END DO
    END DO
    CALL PUSHINTEGER(jed)
    CALL PUSHREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(isd)
    CALL PUSHINTEGER(ie)
    CALL PUSHINTEGER(ied)
    CALL PUSHINTEGER(jsd)
    CALL PUSHREALARRAY(alpha)
    CALL PUSHINTEGER(js)
  END SUBROUTINE GRAD1_P_UPDATE_FWD
!  Differentiation of grad1_p_update in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_e
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
!   gradient     of useful results: gz du u dv v pk divg2
!   with respect to varying inputs: gz du u dv v pk divg2
  SUBROUTINE GRAD1_P_UPDATE_BWD(divg2, divg2_ad, u, u_ad, v, v_ad, pk, &
&   pk_ad, gz, gz_ad, du, du_ad, dv, dv_ad, dt, ng, gridstruct, bd, npx&
&   , npy, npz, ptop, beta, a2b_ord)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz, a2b_ord
    REAL, INTENT(IN) :: dt, ptop, beta
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: divg2_ad(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: du_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dv_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: top_value, alpha
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp2
    REAL :: temp3
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp4
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp8
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    INTEGER :: branch

    wk = 0.0
    top_value = 0.0
    alpha = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    branch = 0

    CALL POPINTEGER(js)
    CALL POPREALARRAY(alpha)
    CALL POPINTEGER(jsd)
    CALL POPINTEGER(ied)
    CALL POPINTEGER(ie)
    CALL POPINTEGER(isd)
    CALL POPINTEGER(is)
    CALL POPINTEGER(je)
    CALL POPREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL POPINTEGER(jed)
    wk_ad = 0.0
    DO k=npz,1,-1
      DO j=je,js,-1
        DO i=ie+1,is,-1
          temp_ad2 = gridstruct%rdy(i, j)*v_ad(i, j, k)
          divg2_ad(i, j) = divg2_ad(i, j) + temp_ad2
          dv_ad(i, j, k) = dv_ad(i, j, k) + alpha*temp_ad2
          divg2_ad(i, j+1) = divg2_ad(i, j+1) - temp_ad2
          v_ad(i, j, k) = temp_ad2
          temp4 = wk(i, j) + wk(i, j+1)
          temp8 = pk(i, j, k+1) - pk(i, j+1, k)
          temp7 = gz(i, j, k) - gz(i, j+1, k+1)
          temp6 = pk(i, j+1, k+1) - pk(i, j, k)
          temp5 = gz(i, j, k+1) - gz(i, j+1, k)
          temp_ad3 = dt*dv_ad(i, j, k)/temp4
          temp_ad4 = -((temp5*temp6+temp7*temp8)*temp_ad3/temp4)
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp6*temp_ad3
          gz_ad(i, j+1, k) = gz_ad(i, j+1, k) - temp6*temp_ad3
          pk_ad(i, j+1, k+1) = pk_ad(i, j+1, k+1) + temp5*temp_ad3
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp5*temp_ad3
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp8*temp_ad3
          gz_ad(i, j+1, k+1) = gz_ad(i, j+1, k+1) - temp8*temp_ad3
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp7*temp_ad3
          pk_ad(i, j+1, k) = pk_ad(i, j+1, k) - temp7*temp_ad3
          wk_ad(i, j) = wk_ad(i, j) + temp_ad4
          wk_ad(i, j+1) = wk_ad(i, j+1) + temp_ad4
          dv_ad(i, j, k) = beta*v_ad(i, j, k)
          CALL POPREALARRAY(v(i, j, k))
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          temp_ad = gridstruct%rdx(i, j)*u_ad(i, j, k)
          divg2_ad(i, j) = divg2_ad(i, j) + temp_ad
          du_ad(i, j, k) = du_ad(i, j, k) + alpha*temp_ad
          divg2_ad(i+1, j) = divg2_ad(i+1, j) - temp_ad
          u_ad(i, j, k) = temp_ad
          temp = wk(i, j) + wk(i+1, j)
          temp3 = pk(i, j, k+1) - pk(i+1, j, k)
          temp2 = gz(i, j, k) - gz(i+1, j, k+1)
          temp1 = pk(i+1, j, k+1) - pk(i, j, k)
          temp0 = gz(i, j, k+1) - gz(i+1, j, k)
          temp_ad0 = dt*du_ad(i, j, k)/temp
          temp_ad1 = -((temp0*temp1+temp2*temp3)*temp_ad0/temp)
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + temp1*temp_ad0
          gz_ad(i+1, j, k) = gz_ad(i+1, j, k) - temp1*temp_ad0
          pk_ad(i+1, j, k+1) = pk_ad(i+1, j, k+1) + temp0*temp_ad0
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp0*temp_ad0
          gz_ad(i, j, k) = gz_ad(i, j, k) + temp3*temp_ad0
          gz_ad(i+1, j, k+1) = gz_ad(i+1, j, k+1) - temp3*temp_ad0
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp2*temp_ad0
          pk_ad(i+1, j, k) = pk_ad(i+1, j, k) - temp2*temp_ad0
          wk_ad(i, j) = wk_ad(i, j) + temp_ad1
          wk_ad(i+1, j) = wk_ad(i+1, j) + temp_ad1
          du_ad(i, j, k) = beta*u_ad(i, j, k)
          CALL POPREALARRAY(u(i, j, k))
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(wk(i, j))
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + wk_ad(i, j)
          pk_ad(i, j, k) = pk_ad(i, j, k) - wk_ad(i, j)
          wk_ad(i, j) = 0.0
        END DO
      END DO
    END DO
    DO k=npz+1,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD2_BWD(gz(isd:ied, jsd:jed, k), gz_ad(isd:ied, jsd:&
&                   jed, k), wk, wk_ad, gridstruct, npx, npy, is, ie, js&
&                   , je, ng, .true.)
      ELSE
        CALL A2B_ORD4_BWD(gz(isd:ied, jsd:jed, k), gz_ad(isd:ied, jsd&
&                      :jed, k), wk, wk_ad, gridstruct, npx, npy, is, ie&
&                      , js, je, ng, .true.)
      END IF
    END DO
    DO k=npz+1,2,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD2_BWD(pk(isd:ied, jsd:jed, k), pk_ad(isd:ied, jsd:&
&                   jed, k), wk, wk_ad, gridstruct, npx, npy, is, ie, js&
&                   , je, ng, .true.)
      ELSE
        CALL A2B_ORD4_BWD(pk(isd:ied, jsd:jed, k), pk_ad(isd:ied, jsd&
&                      :jed, k), wk, wk_ad, gridstruct, npx, npy, is, ie&
&                      , js, je, ng, .true.)
      END IF
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(pk(i, j, 1))
        pk_ad(i, j, 1) = 0.0
      END DO
    END DO
  END SUBROUTINE GRAD1_P_UPDATE_BWD
  SUBROUTINE GRAD1_P_UPDATE(divg2, u, v, pk, gz, du, dv, dt, ng, &
&   gridstruct, bd, npx, npy, npz, ptop, beta, a2b_ord)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz, a2b_ord
    REAL, INTENT(IN) :: dt, ptop, beta
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: top_value, alpha
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    alpha = 1. - beta
! pk is pe**kappa if hydrostatic
    top_value = ptk
!$OMP parallel do default(none) shared(is,ie,js,je,pk,top_value)
    DO j=js,je+1
      DO i=is,ie+1
        pk(i, j, 1) = top_value
      END DO
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,pk,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(pk(isd:ied, jsd:jed, k), wk, gridstruct, npx, &
&                  npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2(pk(isd:ied, jsd:jed, k), wk, gridstruct, npx, npy&
&               , is, ie, js, je, ng, .true.)
      END IF
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,gz,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(gz(isd:ied, jsd:jed, k), wk, gridstruct, npx, &
&                  npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2(gz(isd:ied, jsd:jed, k), wk, gridstruct, npx, npy&
&               , is, ie, js, je, ng, .true.)
      END IF
    END DO
!$OMP parallel do default(none) shared(npz,is,ie,js,je,pk,u,beta,gz,divg2,alpha, &
!$OMP                                  gridstruct,v,dt,du,dv) &
!$OMP                          private(wk)
    DO k=1,npz
      DO j=js,je+1
        DO i=is,ie+1
          wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          u(i, j, k) = u(i, j, k) + beta*du(i, j, k)
          du(i, j, k) = dt/(wk(i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1&
&           , j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, &
&           j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
          u(i, j, k) = (u(i, j, k)+divg2(i, j)-divg2(i+1, j)+alpha*du(i&
&           , j, k))*gridstruct%rdx(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v(i, j, k) = v(i, j, k) + beta*dv(i, j, k)
          dv(i, j, k) = dt/(wk(i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j&
&           +1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1&
&           , k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
          v(i, j, k) = (v(i, j, k)+divg2(i, j)-divg2(i, j+1)+alpha*dv(i&
&           , j, k))*gridstruct%rdy(i, j)
        END DO
      END DO
    END DO
  END SUBROUTINE GRAD1_P_UPDATE
!  Differentiation of mix_dp in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a
!2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.
!p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_d
!p dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Sup
!er fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv
!_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_
!z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_ma
!pz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_m
!apz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_resta
!rt_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z 
!main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.
!Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SI
!M3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nes
!t_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_
!vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v s
!w_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.
!copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod
!.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: w delp pt
!   with respect to varying inputs: w delp pt
  SUBROUTINE MIX_DP_FWD(hydrostatic, w, delp, pt, km, ak, bk, cg, &
&   fv_debug, bd)
    IMPLICIT NONE
!      if ( ip/=0 ) write(*,*) 'Warning: Mix_dp', mpp_pe(), j, ip
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: ak(km+1), bk(km+1)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   pt, delp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   w
    LOGICAL, INTENT(IN) :: hydrostatic, cg, fv_debug
! Local:
    REAL :: dp, dpmin
    INTEGER :: i, j, k, ip
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: res

    dp = 0.0
    dpmin = 0.0
    ifirst = 0
    ilast = 0
    jfirst = 0
    jlast = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    res = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    IF (cg) THEN
      ifirst = is - 1
      ilast = ie + 1
      jfirst = js - 1
      jlast = je + 1
    ELSE
      ifirst = is
      ilast = ie
      jfirst = js
      jlast = je
    END IF
!$OMP parallel do default(none) shared(jfirst,jlast,km,ifirst,ilast,delp,ak,bk,pt, &
!$OMP                                  hydrostatic,w,fv_debug) &
!$OMP                          private(ip, dpmin, dp)
    DO j=jfirst,jlast
      ip = 0
      DO k=1,km-1
        CALL PUSHREALARRAY(dpmin)
        dpmin = 0.01*(ak(k+1)-ak(k)+(bk(k+1)-bk(k))*1.e5)
        DO i=ifirst,ilast
          IF (delp(i, j, k) .LT. dpmin) THEN
!if (fv_debug) write(*,*) 'Mix_dp: ', i, j, k, mpp_pe(), delp(i,j,k), pt(i,j,k)
! Remap from below and mix pt
            CALL PUSHREALARRAY(dp)
            dp = dpmin - delp(i, j, k)
            CALL PUSHREALARRAY(pt(i, j, k))
            pt(i, j, k) = (pt(i, j, k)*delp(i, j, k)+pt(i, j, k+1)*dp)/&
&             dpmin
            IF (.NOT.hydrostatic) THEN
              CALL PUSHREALARRAY(w(i, j, k))
              w(i, j, k) = (w(i, j, k)*delp(i, j, k)+w(i, j, k+1)*dp)/&
&               dpmin
              CALL PUSHCONTROL(1,0)
            ELSE
              CALL PUSHCONTROL(1,1)
            END IF
            CALL PUSHREALARRAY(delp(i, j, k))
            delp(i, j, k) = dpmin
            CALL PUSHREALARRAY(delp(i, j, k+1))
            delp(i, j, k+1) = delp(i, j, k+1) - dp
            ip = ip + 1
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
! Bottom (k=km):
      CALL PUSHREALARRAY(dpmin)
      dpmin = 0.01*(ak(km+1)-ak(km)+(bk(km+1)-bk(km))*1.e5)
      DO i=ifirst,ilast
        IF (delp(i, j, km) .LT. dpmin) THEN
!if (fv_debug) write(*,*) 'Mix_dp: ', i, j, km, mpp_pe(), delp(i,j,km), pt(i,j,km)
! Remap from above and mix pt
          CALL PUSHREALARRAY(dp)
          dp = dpmin - delp(i, j, km)
          CALL PUSHREALARRAY(pt(i, j, km))
          pt(i, j, km) = (pt(i, j, km)*delp(i, j, km)+pt(i, j, km-1)*dp)&
&           /dpmin
          IF (.NOT.hydrostatic) THEN
            CALL PUSHREALARRAY(w(i, j, km))
            w(i, j, km) = (w(i, j, km)*delp(i, j, km)+w(i, j, km-1)*dp)/&
&             dpmin
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
          CALL PUSHREALARRAY(delp(i, j, km))
          delp(i, j, km) = dpmin
          CALL PUSHREALARRAY(delp(i, j, km-1))
          delp(i, j, km-1) = delp(i, j, km-1) - dp
          ip = ip + 1
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      IF (fv_debug .AND. ip .NE. 0) THEN
        CALL PUSHCONTROL(1,0)
        res = MPP_PE()
        WRITE(*, *) 'Warning: Mix_dp', res, j, ip
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
    END DO
    CALL PUSHINTEGER(ifirst)
    CALL PUSHINTEGER(jlast)
    CALL PUSHINTEGER(jfirst)
    CALL PUSHINTEGER(ilast)
    CALL PUSHREALARRAY(dp)
  END SUBROUTINE MIX_DP_FWD
!  Differentiation of mix_dp in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.
!a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod
!.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_
!dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Su
!per fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 f
!v_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap
!_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_m
!apz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_
!mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_rest
!art_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z
! main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod
!.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.S
!IM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.ne
!st_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c
!_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v 
!sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod
!.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mo
!d.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: w delp pt
!   with respect to varying inputs: w delp pt
  SUBROUTINE MIX_DP_BWD(hydrostatic, w, w_ad, delp, delp_ad, pt, pt_ad, &
&   km, ak, bk, cg, fv_debug, bd)
    IMPLICIT NONE
!      if ( ip/=0 ) write(*,*) 'Warning: Mix_dp', mpp_pe(), j, ip
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: ak(km+1), bk(km+1)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   pt, delp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   pt_ad, delp_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   w_ad
    LOGICAL, INTENT(IN) :: hydrostatic, cg, fv_debug
    REAL :: dp, dpmin
    REAL :: dp_ad
    INTEGER :: i, j, k, ip
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    INTEGER :: branch

    dp = 0.0
    dpmin = 0.0
    ifirst = 0
    ilast = 0
    jfirst = 0
    jlast = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    branch = 0

    CALL POPREALARRAY(dp)
    CALL POPINTEGER(ilast)
    CALL POPINTEGER(jfirst)
    CALL POPINTEGER(jlast)
    CALL POPINTEGER(ifirst)
    DO j=jlast,jfirst,-1
      CALL POPCONTROL(1,branch)
      dpmin = 0.01*(ak(km+1)-ak(km)+(bk(km+1)-bk(km))*1.e5)
      DO i=ilast,ifirst,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          CALL POPREALARRAY(delp(i, j, km-1))
          dp_ad = -delp_ad(i, j, km-1)
          CALL POPREALARRAY(delp(i, j, km))
          delp_ad(i, j, km) = 0.0
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(w(i, j, km))
            temp_ad2 = w_ad(i, j, km)/dpmin
            delp_ad(i, j, km) = delp_ad(i, j, km) + w(i, j, km)*temp_ad2
            w_ad(i, j, km-1) = w_ad(i, j, km-1) + dp*temp_ad2
            dp_ad = dp_ad + w(i, j, km-1)*temp_ad2
            w_ad(i, j, km) = delp(i, j, km)*temp_ad2
          END IF
          CALL POPREALARRAY(pt(i, j, km))
          temp_ad1 = pt_ad(i, j, km)/dpmin
          pt_ad(i, j, km-1) = pt_ad(i, j, km-1) + dp*temp_ad1
          dp_ad = dp_ad + pt(i, j, km-1)*temp_ad1
          delp_ad(i, j, km) = delp_ad(i, j, km) + pt(i, j, km)*temp_ad1 &
&           - dp_ad
          pt_ad(i, j, km) = delp(i, j, km)*temp_ad1
          CALL POPREALARRAY(dp)
        END IF
      END DO
      CALL POPREALARRAY(dpmin)
      DO k=km-1,1,-1
        DO i=ilast,ifirst,-1
          CALL POPCONTROL(1,branch)
          IF (branch .NE. 0) THEN
            CALL POPREALARRAY(delp(i, j, k+1))
            dp_ad = -delp_ad(i, j, k+1)
            CALL POPREALARRAY(delp(i, j, k))
            delp_ad(i, j, k) = 0.0
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY(w(i, j, k))
              temp_ad0 = w_ad(i, j, k)/dpmin
              delp_ad(i, j, k) = delp_ad(i, j, k) + w(i, j, k)*temp_ad0
              w_ad(i, j, k+1) = w_ad(i, j, k+1) + dp*temp_ad0
              dp_ad = dp_ad + w(i, j, k+1)*temp_ad0
              w_ad(i, j, k) = delp(i, j, k)*temp_ad0
            END IF
            CALL POPREALARRAY(pt(i, j, k))
            temp_ad = pt_ad(i, j, k)/dpmin
            pt_ad(i, j, k+1) = pt_ad(i, j, k+1) + dp*temp_ad
            dp_ad = dp_ad + pt(i, j, k+1)*temp_ad
            delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*temp_ad - &
&             dp_ad
            pt_ad(i, j, k) = delp(i, j, k)*temp_ad
            CALL POPREALARRAY(dp)
          END IF
        END DO
        CALL POPREALARRAY(dpmin)
      END DO
    END DO
  END SUBROUTINE MIX_DP_BWD
  SUBROUTINE MIX_DP(hydrostatic, w, delp, pt, km, ak, bk, cg, fv_debug, &
&   bd)
    IMPLICIT NONE
!      if ( ip/=0 ) write(*,*) 'Warning: Mix_dp', mpp_pe(), j, ip
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: ak(km+1), bk(km+1)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   pt, delp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   w
    LOGICAL, INTENT(IN) :: hydrostatic, cg, fv_debug
! Local:
    REAL :: dp, dpmin
    INTEGER :: i, j, k, ip
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (cg) THEN
      ifirst = is - 1
      ilast = ie + 1
      jfirst = js - 1
      jlast = je + 1
    ELSE
      ifirst = is
      ilast = ie
      jfirst = js
      jlast = je
    END IF
!$OMP parallel do default(none) shared(jfirst,jlast,km,ifirst,ilast,delp,ak,bk,pt, &
!$OMP                                  hydrostatic,w,fv_debug) &
!$OMP                          private(ip, dpmin, dp)
    DO j=jfirst,jlast
      ip = 0
      DO k=1,km-1
        dpmin = 0.01*(ak(k+1)-ak(k)+(bk(k+1)-bk(k))*1.e5)
        DO i=ifirst,ilast
          IF (delp(i, j, k) .LT. dpmin) THEN
!if (fv_debug) write(*,*) 'Mix_dp: ', i, j, k, mpp_pe(), delp(i,j,k), pt(i,j,k)
! Remap from below and mix pt
            dp = dpmin - delp(i, j, k)
            pt(i, j, k) = (pt(i, j, k)*delp(i, j, k)+pt(i, j, k+1)*dp)/&
&             dpmin
            IF (.NOT.hydrostatic) w(i, j, k) = (w(i, j, k)*delp(i, j, k)&
&               +w(i, j, k+1)*dp)/dpmin
            delp(i, j, k) = dpmin
            delp(i, j, k+1) = delp(i, j, k+1) - dp
            ip = ip + 1
          END IF
        END DO
      END DO
! Bottom (k=km):
      dpmin = 0.01*(ak(km+1)-ak(km)+(bk(km+1)-bk(km))*1.e5)
      DO i=ifirst,ilast
        IF (delp(i, j, km) .LT. dpmin) THEN
!if (fv_debug) write(*,*) 'Mix_dp: ', i, j, km, mpp_pe(), delp(i,j,km), pt(i,j,km)
! Remap from above and mix pt
          dp = dpmin - delp(i, j, km)
          pt(i, j, km) = (pt(i, j, km)*delp(i, j, km)+pt(i, j, km-1)*dp)&
&           /dpmin
          IF (.NOT.hydrostatic) w(i, j, km) = (w(i, j, km)*delp(i, j, km&
&             )+w(i, j, km-1)*dp)/dpmin
          delp(i, j, km) = dpmin
          delp(i, j, km-1) = delp(i, j, km-1) - dp
          ip = ip + 1
        END IF
      END DO
      IF (fv_debug .AND. ip .NE. 0) WRITE(*, *) 'Warning: Mix_dp', &
&                                   MPP_PE(), j, ip
    END DO
  END SUBROUTINE MIX_DP
!  Differentiation of geopk in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2
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
!   gradient     of useful results: peln gz delp pkz pe pk pt
!   with respect to varying inputs: peln gz delp pkz pe pk pt
  SUBROUTINE GEOPK_FWD(ptop, pe, peln, delp, pk, gz, hs, pt, q_con, pkz&
&   , km, akap, cg, nested, computehalo, npx, npy, a2b_ord, bd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km, npx, npy, a2b_ord
    REAL, INTENT(IN) :: akap, ptop
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: hs(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(IN) :: pt&
&   , delp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(IN) :: &
&   q_con
    LOGICAL, INTENT(IN) :: cg, nested, computehalo
! !OUTPUT PARAMETERS
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km+1) :: gz, pk
    REAL :: pe(bd%is-1:bd%ie+1, km+1, bd%js-1:bd%je+1)
! ln(pe)
    REAL :: peln(bd%is:bd%ie, km+1, bd%js:bd%je)
    REAL :: pkz(bd%is:bd%ie, bd%js:bd%je, km)
! !DESCRIPTION:
!    Calculates geopotential and pressure to the kappa.
! Local:
    REAL :: peg(bd%isd:bd%ied, km+1)
    REAL :: pkg(bd%isd:bd%ied, km+1)
    REAL(kind=8) :: p1d(bd%isd:bd%ied)
    REAL(kind=8) :: g1d(bd%isd:bd%ied)
    REAL :: logp(bd%isd:bd%ied)
    INTEGER :: i, j, k
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC LOG
    INTRINSIC EXP
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: ad_from
    INTEGER :: ad_from0

    peg = 0.0
    pkg = 0.0
    p1d = 0.0
    logp = 0.0
    ifirst = 0
    ilast = 0
    jfirst = 0
    jlast = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    max1 = 0
    max2 = 0
    min1 = 0
    min2 = 0
    ad_from = 0
    ad_from0 = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF ((.NOT.cg .AND. a2b_ord .EQ. 4) .OR. (nested .AND. (.NOT.cg))) &
&   THEN
! D-Grid
      ifirst = is - 2
      ilast = ie + 2
      jfirst = js - 2
      jlast = je + 2
    ELSE
      ifirst = is - 1
      ilast = ie + 1
      jfirst = js - 1
      jlast = je + 1
    END IF
    IF (nested .AND. computehalo) THEN
      IF (is .EQ. 1) ifirst = isd
      IF (ie .EQ. npx - 1) ilast = ied
      IF (js .EQ. 1) jfirst = jsd
      IF (je .EQ. npy - 1) THEN
        CALL PUSHCONTROL(1,1)
        jlast = jed
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
!$OMP parallel do default(none) shared(jfirst,jlast,ifirst,ilast,pk,km,gz,hs,ptop,ptk, &
!$OMP                                  js,je,is,ie,peln,peln1,pe,delp,akap,pt,CG,pkz,q_con) &
!$OMP                          private(peg, pkg, p1d, g1d, logp)
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        CALL PUSHREALARRAY(p1d(i))
        p1d(i) = ptop
        CALL PUSHREALARRAY(pk(i, j, 1))
        pk(i, j, 1) = ptk
        g1d(i) = hs(i, j)
        CALL PUSHREALARRAY(gz(i, j, km+1))
        gz(i, j, km+1) = hs(i, j)
      END DO
      IF (j .GE. js .AND. j .LE. je) THEN
        DO i=is,ie
          CALL PUSHREALARRAY(peln(i, 1, j))
          peln(i, 1, j) = peln1
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (j .GT. js - 2 .AND. j .LT. je + 2) THEN
        IF (ifirst .LT. is - 1) THEN
          max1 = is - 1
        ELSE
          max1 = ifirst
        END IF
        IF (ilast .GT. ie + 1) THEN
          min1 = ie + 1
        ELSE
          min1 = ilast
        END IF
        ad_from = max1
        DO i=ad_from,min1
          CALL PUSHREALARRAY(pe(i, 1, j))
          pe(i, 1, j) = ptop
        END DO
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
! Top down
      DO k=2,km+1
        DO i=ifirst,ilast
          CALL PUSHREALARRAY(p1d(i))
          p1d(i) = p1d(i) + delp(i, j, k-1)
          CALL PUSHREALARRAY(logp(i))
          logp(i) = LOG(p1d(i))
          CALL PUSHREALARRAY(pk(i, j, k))
          pk(i, j, k) = EXP(akap*logp(i))
        END DO
        IF (j .GT. js - 2 .AND. j .LT. je + 2) THEN
          IF (ifirst .LT. is - 1) THEN
            max2 = is - 1
          ELSE
            max2 = ifirst
          END IF
          IF (ilast .GT. ie + 1) THEN
            min2 = ie + 1
          ELSE
            min2 = ilast
          END IF
          ad_from0 = max2
          DO i=ad_from0,min2
            CALL PUSHREALARRAY(pe(i, k, j))
            pe(i, k, j) = p1d(i)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from0)
          IF (j .GE. js .AND. j .LE. je) THEN
            DO i=is,ie
              CALL PUSHREALARRAY(peln(i, k, j))
              peln(i, k, j) = logp(i)
            END DO
            CALL PUSHCONTROL(2,2)
          ELSE
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE
          CALL PUSHCONTROL(2,0)
        END IF
      END DO
! Bottom up
      DO k=km,1,-1
        DO i=ifirst,ilast
          g1d(i) = g1d(i) + cp_air*pt(i, j, k)*(pk(i, j, k+1)-pk(i, j, k&
&           ))
          CALL PUSHREALARRAY(gz(i, j, k))
          gz(i, j, k) = g1d(i)
        END DO
      END DO
      IF (.NOT.cg .AND. j .GE. js .AND. j .LE. je) THEN
        DO k=1,km
          DO i=is,ie
            CALL PUSHREALARRAY(pkz(i, j, k))
            pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(akap*(peln(i, k+&
&             1, j)-peln(i, k, j)))
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
    END DO
    CALL PUSHINTEGER(ifirst)
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(ie)
    CALL PUSHINTEGER(jlast)
    CALL PUSHREALARRAY(p1d, bd%ied - bd%isd + 1)
    CALL PUSHINTEGER(jfirst)
    CALL PUSHINTEGER(ilast)
    CALL PUSHREALARRAY(logp, bd%ied - bd%isd + 1)
  END SUBROUTINE GEOPK_FWD
!  Differentiation of geopk in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a
!2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.
!p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_d
!p dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Sup
!er fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv
!_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_
!z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_ma
!pz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_m
!apz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_resta
!rt_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z 
!main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.
!Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SI
!M3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nes
!t_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_
!vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v s
!w_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.
!copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod
!.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: peln gz delp pkz pe pk pt
!   with respect to varying inputs: peln gz delp pkz pe pk pt
  SUBROUTINE GEOPK_BWD(ptop, pe, pe_ad, peln, peln_ad, delp, delp_ad, pk&
&   , pk_ad, gz, gz_ad, hs, pt, pt_ad, q_con, pkz, pkz_ad, km, akap, cg&
&   , nested, computehalo, npx, npy, a2b_ord, bd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km, npx, npy, a2b_ord
    REAL, INTENT(IN) :: akap, ptop
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: hs(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(IN) :: pt&
&   , delp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km) :: pt_ad, delp_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(IN) :: &
&   q_con
    LOGICAL, INTENT(IN) :: cg, nested, computehalo
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km+1) :: gz, pk
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km+1) :: gz_ad, pk_ad
    REAL :: pe(bd%is-1:bd%ie+1, km+1, bd%js-1:bd%je+1)
    REAL :: pe_ad(bd%is-1:bd%ie+1, km+1, bd%js-1:bd%je+1)
    REAL :: peln(bd%is:bd%ie, km+1, bd%js:bd%je)
    REAL :: peln_ad(bd%is:bd%ie, km+1, bd%js:bd%je)
    REAL :: pkz(bd%is:bd%ie, bd%js:bd%je, km)
    REAL :: pkz_ad(bd%is:bd%ie, bd%js:bd%je, km)
    REAL :: peg(bd%isd:bd%ied, km+1)
    REAL :: pkg(bd%isd:bd%ied, km+1)
    REAL(kind=8) :: p1d(bd%isd:bd%ied)
    REAL(kind=8) :: p1d_ad(bd%isd:bd%ied)
    REAL(kind=8) :: g1d(bd%isd:bd%ied)
    REAL(kind=8) :: g1d_ad(bd%isd:bd%ied)
    REAL :: logp(bd%isd:bd%ied)
    REAL :: logp_ad(bd%isd:bd%ied)
    INTEGER :: i, j, k
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC LOG
    INTRINSIC EXP
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: min1
    INTEGER :: min2
    REAL :: temp_ad
    REAL :: temp
    REAL :: temp_ad0
    REAL :: temp_ad1
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: branch

    peg = 0.0
    pkg = 0.0
    p1d = 0.0
    logp = 0.0
    ifirst = 0
    ilast = 0
    jfirst = 0
    jlast = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    max1 = 0
    max2 = 0
    min1 = 0
    min2 = 0
    ad_from = 0
    ad_from0 = 0
    ad_to = 0
    ad_to0 = 0
    branch = 0

    CALL POPREALARRAY(logp, bd%ied - bd%isd + 1)
    CALL POPINTEGER(ilast)
    CALL POPINTEGER(jfirst)
    CALL POPREALARRAY(p1d, bd%ied - bd%isd + 1)
    CALL POPINTEGER(jlast)
    CALL POPINTEGER(ie)
    CALL POPINTEGER(is)
    CALL POPINTEGER(ifirst)
    g1d_ad = 0.0_8
    logp_ad = 0.0
    p1d_ad = 0.0_8
    DO j=jlast,jfirst,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO k=km,1,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(pkz(i, j, k))
            temp = akap*(peln(i, k+1, j)-peln(i, k, j))
            temp_ad0 = pkz_ad(i, j, k)/temp
            temp_ad1 = -((pk(i, j, k+1)-pk(i, j, k))*akap*temp_ad0/temp)
            pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp_ad0
            pk_ad(i, j, k) = pk_ad(i, j, k) - temp_ad0
            peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad1
            peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad1
            pkz_ad(i, j, k) = 0.0
          END DO
        END DO
      END IF
      DO k=1,km,1
        DO i=ilast,ifirst,-1
          CALL POPREALARRAY(gz(i, j, k))
          g1d_ad(i) = g1d_ad(i) + gz_ad(i, j, k)
          gz_ad(i, j, k) = 0.0
          temp_ad = cp_air*pt(i, j, k)*g1d_ad(i)
          pt_ad(i, j, k) = pt_ad(i, j, k) + cp_air*(pk(i, j, k+1)-pk(i, &
&           j, k))*g1d_ad(i)
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp_ad
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp_ad
        END DO
      END DO
      DO k=km+1,2,-1
        CALL POPCONTROL(2,branch)
        IF (branch .NE. 0) THEN
          IF (branch .NE. 1) THEN
            DO i=ie,is,-1
              CALL POPREALARRAY(peln(i, k, j))
              logp_ad(i) = logp_ad(i) + peln_ad(i, k, j)
              peln_ad(i, k, j) = 0.0
            END DO
          END IF
          CALL POPINTEGER(ad_from0)
          CALL POPINTEGER(ad_to0)
          DO i=ad_to0,ad_from0,-1
            CALL POPREALARRAY(pe(i, k, j))
            p1d_ad(i) = p1d_ad(i) + pe_ad(i, k, j)
            pe_ad(i, k, j) = 0.0
          END DO
        END IF
        DO i=ilast,ifirst,-1
          CALL POPREALARRAY(pk(i, j, k))
          logp_ad(i) = logp_ad(i) + EXP(akap*logp(i))*akap*pk_ad(i, j, k&
&           )
          pk_ad(i, j, k) = 0.0
          CALL POPREALARRAY(logp(i))
          p1d_ad(i) = p1d_ad(i) + logp_ad(i)/p1d(i)
          logp_ad(i) = 0.0
          CALL POPREALARRAY(p1d(i))
          delp_ad(i, j, k-1) = delp_ad(i, j, k-1) + p1d_ad(i)
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        CALL POPINTEGER(ad_from)
        CALL POPINTEGER(ad_to)
        DO i=ad_to,ad_from,-1
          CALL POPREALARRAY(pe(i, 1, j))
          pe_ad(i, 1, j) = 0.0
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO i=ie,is,-1
          CALL POPREALARRAY(peln(i, 1, j))
          peln_ad(i, 1, j) = 0.0
        END DO
      END IF
      DO i=ilast,ifirst,-1
        CALL POPREALARRAY(gz(i, j, km+1))
        gz_ad(i, j, km+1) = 0.0
        g1d_ad(i) = 0.0_8
        CALL POPREALARRAY(pk(i, j, 1))
        pk_ad(i, j, 1) = 0.0
        CALL POPREALARRAY(p1d(i))
        p1d_ad(i) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL(1,branch)
  END SUBROUTINE GEOPK_BWD
  SUBROUTINE GEOPK(ptop, pe, peln, delp, pk, gz, hs, pt, q_con, pkz, km&
&   , akap, cg, nested, computehalo, npx, npy, a2b_ord, bd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km, npx, npy, a2b_ord
    REAL, INTENT(IN) :: akap, ptop
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: hs(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(IN) :: pt&
&   , delp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(IN) :: &
&   q_con
    LOGICAL, INTENT(IN) :: cg, nested, computehalo
! !OUTPUT PARAMETERS
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km+1), INTENT(OUT) :: &
&   gz, pk
    REAL, INTENT(OUT) :: pe(bd%is-1:bd%ie+1, km+1, bd%js-1:bd%je+1)
! ln(pe)
    REAL, INTENT(OUT) :: peln(bd%is:bd%ie, km+1, bd%js:bd%je)
    REAL, INTENT(OUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, km)
! !DESCRIPTION:
!    Calculates geopotential and pressure to the kappa.
! Local:
    REAL :: peg(bd%isd:bd%ied, km+1)
    REAL :: pkg(bd%isd:bd%ied, km+1)
    REAL(kind=8) :: p1d(bd%isd:bd%ied)
    REAL(kind=8) :: g1d(bd%isd:bd%ied)
    REAL :: logp(bd%isd:bd%ied)
    INTEGER :: i, j, k
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC LOG
    INTRINSIC EXP
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: min1
    INTEGER :: min2
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF ((.NOT.cg .AND. a2b_ord .EQ. 4) .OR. (nested .AND. (.NOT.cg))) &
&   THEN
! D-Grid
      ifirst = is - 2
      ilast = ie + 2
      jfirst = js - 2
      jlast = je + 2
    ELSE
      ifirst = is - 1
      ilast = ie + 1
      jfirst = js - 1
      jlast = je + 1
    END IF
    IF (nested .AND. computehalo) THEN
      IF (is .EQ. 1) ifirst = isd
      IF (ie .EQ. npx - 1) ilast = ied
      IF (js .EQ. 1) jfirst = jsd
      IF (je .EQ. npy - 1) jlast = jed
    END IF
!$OMP parallel do default(none) shared(jfirst,jlast,ifirst,ilast,pk,km,gz,hs,ptop,ptk, &
!$OMP                                  js,je,is,ie,peln,peln1,pe,delp,akap,pt,CG,pkz,q_con) &
!$OMP                          private(peg, pkg, p1d, g1d, logp)
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        p1d(i) = ptop
        pk(i, j, 1) = ptk
        g1d(i) = hs(i, j)
        gz(i, j, km+1) = hs(i, j)
      END DO
      IF (j .GE. js .AND. j .LE. je) THEN
        DO i=is,ie
          peln(i, 1, j) = peln1
        END DO
      END IF
      IF (j .GT. js - 2 .AND. j .LT. je + 2) THEN
        IF (ifirst .LT. is - 1) THEN
          max1 = is - 1
        ELSE
          max1 = ifirst
        END IF
        IF (ilast .GT. ie + 1) THEN
          min1 = ie + 1
        ELSE
          min1 = ilast
        END IF
        DO i=max1,min1
          pe(i, 1, j) = ptop
        END DO
      END IF
! Top down
      DO k=2,km+1
        DO i=ifirst,ilast
          p1d(i) = p1d(i) + delp(i, j, k-1)
          logp(i) = LOG(p1d(i))
          pk(i, j, k) = EXP(akap*logp(i))
        END DO
        IF (j .GT. js - 2 .AND. j .LT. je + 2) THEN
          IF (ifirst .LT. is - 1) THEN
            max2 = is - 1
          ELSE
            max2 = ifirst
          END IF
          IF (ilast .GT. ie + 1) THEN
            min2 = ie + 1
          ELSE
            min2 = ilast
          END IF
          DO i=max2,min2
            pe(i, k, j) = p1d(i)
          END DO
          IF (j .GE. js .AND. j .LE. je) THEN
            DO i=is,ie
              peln(i, k, j) = logp(i)
            END DO
          END IF
        END IF
      END DO
! Bottom up
      DO k=km,1,-1
        DO i=ifirst,ilast
          g1d(i) = g1d(i) + cp_air*pt(i, j, k)*(pk(i, j, k+1)-pk(i, j, k&
&           ))
          gz(i, j, k) = g1d(i)
        END DO
      END DO
      IF (.NOT.cg .AND. j .GE. js .AND. j .LE. je) THEN
        DO k=1,km
          DO i=is,ie
            pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(akap*(peln(i, k+&
&             1, j)-peln(i, k, j)))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE GEOPK
!  Differentiation of del2_cubed in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_m
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
!   gradient     of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE DEL2_CUBED_FWD(q, cd, gridstruct, domain, npx, npy, km, &
&   nmax, bd)
    IMPLICIT NONE
!---------------------------------------------------------------
! This routine is for filtering the omega field for the physics
!---------------------------------------------------------------
    INTEGER, INTENT(IN) :: npx, npy, km, nmax
! cd = K * da_min;   0 < K < 0.25
    REAL(kind=r_grid), INTENT(IN) :: cd
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL, PARAMETER :: r3=1./3.
    REAL :: fx(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy(bd%isd:bd%ied, bd%jsd&
&   :bd%jed+1)
    REAL :: q2(bd%isd:bd%ied, bd%jsd:bd%jed)
    INTEGER :: i, j, k, n, nt, ntimes
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MIN
    INTEGER :: ad_from
    INTEGER :: ad_from0
    INTEGER :: ad_from1
    INTEGER :: ad_from2
    INTEGER :: ad_from3
    INTEGER :: ad_from4
    REAL :: tmp
    REAL :: tmp0
    REAL :: tmp1
    REAL :: tmp2
    REAL :: tmp3
    REAL :: tmp4
    REAL :: tmp5

    fx = 0.0
    fy = 0.0
    q2 = 0.0
    nt = 0
    ntimes = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    ad_from = 0
    ad_from0 = 0
    ad_from1 = 0
    ad_from2 = 0
    ad_from3 = 0
    ad_from4 = 0
    tmp = 0.0
    tmp0 = 0.0
    tmp1 = 0.0
    tmp2 = 0.0
    tmp3 = 0.0
    tmp4 = 0.0
    tmp5 = 0.0

!Local routine pointers
!     real, pointer, dimension(:,:) :: rarea
!     real, pointer, dimension(:,:) :: del6_u, del6_v
!     logical, pointer :: sw_corner, se_corner, ne_corner, nw_corner
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (3 .GT. nmax) THEN
      ntimes = nmax
    ELSE
      ntimes = 3
    END IF
    CALL PUSHREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*km)
    CALL MPP_UPDATE_DOMAINS(q, domain, complete=.true.)
    DO n=1,ntimes
      nt = ntimes - n
!$OMP parallel do default(none) shared(km,q,is,ie,js,je,npx,npy, &
!$OMP                                  nt,isd,jsd,gridstruct,bd, &
!$OMP                                  cd) &
!$OMP                          private(fx, fy)
      DO k=1,km
        IF (gridstruct%sw_corner) THEN
          q(1, 1, k) = (q(1, 1, k)+q(0, 1, k)+q(1, 0, k))*r3
          q(0, 1, k) = q(1, 1, k)
          q(1, 0, k) = q(1, 1, k)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%se_corner) THEN
          tmp0 = (q(ie, 1, k)+q(npx, 1, k)+q(ie, 0, k))*r3
          q(ie, 1, k) = tmp0
          tmp = q(ie, 1, k)
          q(npx, 1, k) = tmp
          q(ie, 0, k) = q(ie, 1, k)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%ne_corner) THEN
          tmp3 = (q(ie, je, k)+q(npx, je, k)+q(ie, npy, k))*r3
          q(ie, je, k) = tmp3
          tmp2 = q(ie, je, k)
          q(npx, je, k) = tmp2
          tmp1 = q(ie, je, k)
          q(ie, npy, k) = tmp1
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%nw_corner) THEN
          tmp5 = (q(1, je, k)+q(0, je, k)+q(1, npy, k))*r3
          q(1, je, k) = tmp5
          q(0, je, k) = q(1, je, k)
          tmp4 = q(1, je, k)
          q(1, npy, k) = tmp4
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (nt .GT. 0) THEN
          CALL COPY_CORNERS(q(isd:ied, jsd:jed, k), npx, npy, 1, &
&                     gridstruct%nested, bd, gridstruct%sw_corner, &
&                     gridstruct%se_corner, gridstruct%nw_corner, &
&                     gridstruct%ne_corner)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        ad_from0 = js - nt
        DO j=ad_from0,je+nt
          ad_from = is - nt
          DO i=ad_from,ie+1+nt
            fx(i, j) = gridstruct%del6_v(i, j)*(q(i-1, j, k)-q(i, j, k))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from)
        END DO
        CALL PUSHINTEGER(j - 1)
        CALL PUSHINTEGER(ad_from0)
        IF (nt .GT. 0) THEN
          CALL COPY_CORNERS(q(isd:ied, jsd:jed, k), npx, npy, 2, &
&                     gridstruct%nested, bd, gridstruct%sw_corner, &
&                     gridstruct%se_corner, gridstruct%nw_corner, &
&                     gridstruct%ne_corner)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        ad_from2 = js - nt
        DO j=ad_from2,je+1+nt
          ad_from1 = is - nt
          DO i=ad_from1,ie+nt
            fy(i, j) = gridstruct%del6_u(i, j)*(q(i, j-1, k)-q(i, j, k))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from1)
        END DO
        CALL PUSHINTEGER(j - 1)
        CALL PUSHINTEGER(ad_from2)
        ad_from4 = js - nt
        DO j=ad_from4,je+nt
          ad_from3 = is - nt
          DO i=ad_from3,ie+nt
            q(i, j, k) = q(i, j, k) + cd*gridstruct%rarea(i, j)*(fx(i, j&
&             )-fx(i+1, j)+fy(i, j)-fy(i, j+1))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from3)
        END DO
        CALL PUSHINTEGER(j - 1)
        CALL PUSHINTEGER(ad_from4)
      END DO
    END DO
    CALL PUSHINTEGER(jed)
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(isd)
    CALL PUSHINTEGER(ie)
    CALL PUSHINTEGER(ied)
    CALL PUSHINTEGER(jsd)
    CALL PUSHINTEGER(ntimes)
  END SUBROUTINE DEL2_CUBED_FWD
!  Differentiation of del2_cubed in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE DEL2_CUBED_BWD(q, q_ad, cd, gridstruct, domain, npx, npy, &
&   km, nmax, bd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, km, nmax
    REAL(kind=r_grid), INTENT(IN) :: cd
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL, PARAMETER :: r3=1./3.
    REAL :: fx(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy(bd%isd:bd%ied, bd%jsd&
&   :bd%jed+1)
    REAL :: fx_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy_ad(bd%isd:bd%ied, &
&   bd%jsd:bd%jed+1)
    REAL :: q2(bd%isd:bd%ied, bd%jsd:bd%jed)
    INTEGER :: i, j, k, n, nt, ntimes
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MIN
    REAL :: temp_ad
    REAL :: tmp_ad
    REAL :: tmp_ad0
    REAL :: temp_ad0
    REAL :: tmp_ad1
    REAL :: tmp_ad2
    REAL :: tmp_ad3
    REAL :: temp_ad1
    REAL :: tmp_ad4
    REAL :: tmp_ad5
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: ad_from1
    INTEGER :: ad_to1
    INTEGER :: ad_from2
    INTEGER :: ad_to2
    INTEGER :: ad_from3
    INTEGER :: ad_to3
    INTEGER :: ad_from4
    INTEGER :: ad_to4
    INTEGER :: branch

    fx = 0.0
    fy = 0.0
    q2 = 0.0
    nt = 0
    ntimes = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    ad_from = 0
    ad_from0 = 0
    ad_from1 = 0
    ad_from2 = 0
    ad_from3 = 0
    ad_from4 = 0
    ad_to = 0
    ad_to0 = 0
    ad_to1 = 0
    ad_to2 = 0
    ad_to3 = 0
    ad_to4 = 0
    branch = 0

    CALL POPINTEGER(ntimes)
    CALL POPINTEGER(jsd)
    CALL POPINTEGER(ied)
    CALL POPINTEGER(ie)
    CALL POPINTEGER(isd)
    CALL POPINTEGER(je)
    CALL POPINTEGER(jed)
    fx_ad = 0.0
    fy_ad = 0.0
    DO n=ntimes,1,-1
      DO k=km,1,-1
        CALL POPINTEGER(ad_from4)
        CALL POPINTEGER(ad_to4)
        DO j=ad_to4,ad_from4,-1
          CALL POPINTEGER(ad_from3)
          CALL POPINTEGER(ad_to3)
          DO i=ad_to3,ad_from3,-1
            temp_ad5 = cd*gridstruct%rarea(i, j)*q_ad(i, j, k)
            fx_ad(i, j) = fx_ad(i, j) + temp_ad5
            fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad5
            fy_ad(i, j) = fy_ad(i, j) + temp_ad5
            fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad5
          END DO
        END DO
        CALL POPINTEGER(ad_from2)
        CALL POPINTEGER(ad_to2)
        DO j=ad_to2,ad_from2,-1
          CALL POPINTEGER(ad_from1)
          CALL POPINTEGER(ad_to1)
          DO i=ad_to1,ad_from1,-1
            temp_ad4 = gridstruct%del6_u(i, j)*fy_ad(i, j)
            q_ad(i, j-1, k) = q_ad(i, j-1, k) + temp_ad4
            q_ad(i, j, k) = q_ad(i, j, k) - temp_ad4
            fy_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) CALL COPY_CORNERS_ADM(q(isd:ied, jsd:jed, k)&
&                                          , q_ad(isd:ied, jsd:jed, k), &
&                                          npx, npy, 2, gridstruct%&
&                                          nested, bd, gridstruct%&
&                                          sw_corner, gridstruct%&
&                                          se_corner, gridstruct%&
&                                          nw_corner, gridstruct%&
&                                          ne_corner)
        CALL POPINTEGER(ad_from0)
        CALL POPINTEGER(ad_to0)
        DO j=ad_to0,ad_from0,-1
          CALL POPINTEGER(ad_from)
          CALL POPINTEGER(ad_to)
          DO i=ad_to,ad_from,-1
            temp_ad3 = gridstruct%del6_v(i, j)*fx_ad(i, j)
            q_ad(i-1, j, k) = q_ad(i-1, j, k) + temp_ad3
            q_ad(i, j, k) = q_ad(i, j, k) - temp_ad3
            fx_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) CALL COPY_CORNERS_ADM(q(isd:ied, jsd:jed, k)&
&                                          , q_ad(isd:ied, jsd:jed, k), &
&                                          npx, npy, 1, gridstruct%&
&                                          nested, bd, gridstruct%&
&                                          sw_corner, gridstruct%&
&                                          se_corner, gridstruct%&
&                                          nw_corner, gridstruct%&
&                                          ne_corner)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          tmp_ad4 = q_ad(1, npy, k)
          q_ad(1, npy, k) = 0.0
          q_ad(1, je, k) = q_ad(1, je, k) + q_ad(0, je, k) + tmp_ad4
          q_ad(0, je, k) = 0.0
          tmp_ad5 = q_ad(1, je, k)
          temp_ad2 = r3*tmp_ad5
          q_ad(1, je, k) = temp_ad2
          q_ad(0, je, k) = q_ad(0, je, k) + temp_ad2
          q_ad(1, npy, k) = q_ad(1, npy, k) + temp_ad2
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          tmp_ad1 = q_ad(ie, npy, k)
          q_ad(ie, npy, k) = 0.0
          q_ad(ie, je, k) = q_ad(ie, je, k) + tmp_ad1
          tmp_ad2 = q_ad(npx, je, k)
          q_ad(npx, je, k) = 0.0
          q_ad(ie, je, k) = q_ad(ie, je, k) + tmp_ad2
          tmp_ad3 = q_ad(ie, je, k)
          temp_ad1 = r3*tmp_ad3
          q_ad(ie, je, k) = temp_ad1
          q_ad(npx, je, k) = q_ad(npx, je, k) + temp_ad1
          q_ad(ie, npy, k) = q_ad(ie, npy, k) + temp_ad1
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          q_ad(ie, 1, k) = q_ad(ie, 1, k) + q_ad(ie, 0, k)
          q_ad(ie, 0, k) = 0.0
          tmp_ad = q_ad(npx, 1, k)
          q_ad(npx, 1, k) = 0.0
          q_ad(ie, 1, k) = q_ad(ie, 1, k) + tmp_ad
          tmp_ad0 = q_ad(ie, 1, k)
          temp_ad0 = r3*tmp_ad0
          q_ad(ie, 1, k) = temp_ad0
          q_ad(npx, 1, k) = q_ad(npx, 1, k) + temp_ad0
          q_ad(ie, 0, k) = q_ad(ie, 0, k) + temp_ad0
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          q_ad(1, 1, k) = q_ad(1, 1, k) + q_ad(1, 0, k)
          q_ad(1, 0, k) = 0.0
          q_ad(1, 1, k) = q_ad(1, 1, k) + q_ad(0, 1, k)
          q_ad(0, 1, k) = 0.0
          temp_ad = r3*q_ad(1, 1, k)
          q_ad(0, 1, k) = q_ad(0, 1, k) + temp_ad
          q_ad(1, 0, k) = q_ad(1, 0, k) + temp_ad
          q_ad(1, 1, k) = temp_ad
        END IF
      END DO
    END DO
    CALL POPREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*km)
    CALL MPP_UPDATE_DOMAINS_ADM(q, q_ad, domain, complete=.true.)
  END SUBROUTINE DEL2_CUBED_BWD
  SUBROUTINE DEL2_CUBED(q, cd, gridstruct, domain, npx, npy, km, nmax, &
&   bd)
    IMPLICIT NONE
!---------------------------------------------------------------
! This routine is for filtering the omega field for the physics
!---------------------------------------------------------------
    INTEGER, INTENT(IN) :: npx, npy, km, nmax
! cd = K * da_min;   0 < K < 0.25
    REAL(kind=r_grid), INTENT(IN) :: cd
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL, PARAMETER :: r3=1./3.
    REAL :: fx(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy(bd%isd:bd%ied, bd%jsd&
&   :bd%jed+1)
    REAL :: q2(bd%isd:bd%ied, bd%jsd:bd%jed)
    INTEGER :: i, j, k, n, nt, ntimes
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MIN
!Local routine pointers
!     real, pointer, dimension(:,:) :: rarea
!     real, pointer, dimension(:,:) :: del6_u, del6_v
!     logical, pointer :: sw_corner, se_corner, ne_corner, nw_corner
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (3 .GT. nmax) THEN
      ntimes = nmax
    ELSE
      ntimes = 3
    END IF
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_DOMAINS(q, domain, complete=.true.)
    CALL TIMING_OFF('COMM_TOTAL')
    DO n=1,ntimes
      nt = ntimes - n
!$OMP parallel do default(none) shared(km,q,is,ie,js,je,npx,npy, &
!$OMP                                  nt,isd,jsd,gridstruct,bd, &
!$OMP                                  cd) &
!$OMP                          private(fx, fy)
      DO k=1,km
        IF (gridstruct%sw_corner) THEN
          q(1, 1, k) = (q(1, 1, k)+q(0, 1, k)+q(1, 0, k))*r3
          q(0, 1, k) = q(1, 1, k)
          q(1, 0, k) = q(1, 1, k)
        END IF
        IF (gridstruct%se_corner) THEN
          q(ie, 1, k) = (q(ie, 1, k)+q(npx, 1, k)+q(ie, 0, k))*r3
          q(npx, 1, k) = q(ie, 1, k)
          q(ie, 0, k) = q(ie, 1, k)
        END IF
        IF (gridstruct%ne_corner) THEN
          q(ie, je, k) = (q(ie, je, k)+q(npx, je, k)+q(ie, npy, k))*r3
          q(npx, je, k) = q(ie, je, k)
          q(ie, npy, k) = q(ie, je, k)
        END IF
        IF (gridstruct%nw_corner) THEN
          q(1, je, k) = (q(1, je, k)+q(0, je, k)+q(1, npy, k))*r3
          q(0, je, k) = q(1, je, k)
          q(1, npy, k) = q(1, je, k)
        END IF
        IF (nt .GT. 0) CALL COPY_CORNERS(q(isd:ied, jsd:jed, k), npx, &
&                                  npy, 1, gridstruct%nested, bd, &
&                                  gridstruct%sw_corner, gridstruct%&
&                                  se_corner, gridstruct%nw_corner, &
&                                  gridstruct%ne_corner)
        DO j=js-nt,je+nt
          DO i=is-nt,ie+1+nt
            fx(i, j) = gridstruct%del6_v(i, j)*(q(i-1, j, k)-q(i, j, k))
          END DO
        END DO
        IF (nt .GT. 0) CALL COPY_CORNERS(q(isd:ied, jsd:jed, k), npx, &
&                                  npy, 2, gridstruct%nested, bd, &
&                                  gridstruct%sw_corner, gridstruct%&
&                                  se_corner, gridstruct%nw_corner, &
&                                  gridstruct%ne_corner)
        DO j=js-nt,je+1+nt
          DO i=is-nt,ie+nt
            fy(i, j) = gridstruct%del6_u(i, j)*(q(i, j-1, k)-q(i, j, k))
          END DO
        END DO
        DO j=js-nt,je+nt
          DO i=is-nt,ie+nt
            q(i, j, k) = q(i, j, k) + cd*gridstruct%rarea(i, j)*(fx(i, j&
&             )-fx(i+1, j)+fy(i, j)-fy(i, j+1))
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE DEL2_CUBED
  SUBROUTINE INIT_IJK_MEM(i1, i2, j1, j2, km, array, var)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1, i2, j1, j2, km
    REAL, INTENT(INOUT) :: array(i1:i2, j1:j2, km)
    REAL, INTENT(IN) :: var
    INTEGER :: i, j, k
!$OMP parallel do default(none) shared(i1,i2,j1,j2,km,array,var)
    DO k=1,km
      DO j=j1,j2
        DO i=i1,i2
          array(i, j, k) = var
        END DO
      END DO
    END DO
  END SUBROUTINE INIT_IJK_MEM
  SUBROUTINE RAYLEIGH_FAST(dt, npx, npy, npz, pfull, tau, u, v, w, ptop&
&   , hydrostatic, rf_cutoff, bd)
    IMPLICIT NONE
! Simple "inline" version of the Rayleigh friction
    REAL, INTENT(IN) :: dt
! time scale (days)
    REAL, INTENT(IN) :: tau
    REAL, INTENT(IN) :: ptop, rf_cutoff
    INTEGER, INTENT(IN) :: npx, npy, npz
    REAL, DIMENSION(npz), INTENT(IN) :: pfull
    LOGICAL, INTENT(IN) :: hydrostatic
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
! cell center vertical wind (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:, bd%jsd:, :)
!
    REAL(kind=r_grid) :: rff(npz)
    REAL, PARAMETER :: sday=86400.
    REAL :: tau0
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL :: rf(npz)
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
    IF (.NOT.rff_initialized) THEN
      tau0 = tau*sday
!allocate( rf(npz) )
      rf(:) = 1.
      IF (IS_MASTER()) WRITE(6, *) &
&                      'Fast Rayleigh friction E-folding time (days):'
      DO k=1,npz
        IF (pfull(k) .LT. rf_cutoff) THEN
          rff(k) = dt/tau0*SIN(0.5*pi*LOG(rf_cutoff/pfull(k))/LOG(&
&           rf_cutoff/ptop))**2
! Re-FACTOR rf
!if( is_master() ) write(6,*) k, 0.01*pfull(k), dt/(rff(k)*sday)
          kmax = k
          rff(k) = 1.d0/(1.0d0+rff(k))
          rf(k) = rff(k)
        ELSE
          GOTO 100
        END IF
      END DO
 100  rff_initialized = .true.
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,pfull,rf_cutoff,w,rf,u,v,hydrostatic)
    DO k=1,kmax
      IF (pfull(k) .LT. rf_cutoff) THEN
        DO j=js,je+1
          DO i=is,ie
            u(i, j, k) = rf(k)*u(i, j, k)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            v(i, j, k) = rf(k)*v(i, j, k)
          END DO
        END DO
        IF (.NOT.hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              w(i, j, k) = rf(k)*w(i, j, k)
            END DO
          END DO
        END IF
      END IF
    END DO
  END SUBROUTINE RAYLEIGH_FAST
end module dyn_core_adm_mod
