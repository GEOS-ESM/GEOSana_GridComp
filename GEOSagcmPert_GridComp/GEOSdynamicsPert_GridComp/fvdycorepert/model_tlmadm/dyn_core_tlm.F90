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
module dyn_core_tlm_mod

  use constants_mod,      only: rdgas, radius, cp_air, pi
  use mpp_mod,            only: mpp_pe, mpp_root_pe
  use mpp_domains_mod,    only: CGRID_NE, DGRID_NE, mpp_get_boundary, mpp_update_domains, &
                                domain2d
  use fv_mp_tlm_mod,      only: mpp_get_boundary_tlm, mpp_update_domains_tlm
  use mpp_parameter_mod,  only: CORNER
  use fv_mp_mod,          only: is_master
  use fv_mp_tlm_mod,      only: start_group_halo_update, complete_group_halo_update
  use fv_mp_tlm_mod,      only: start_group_halo_update_tlm
  use fv_mp_mod,          only: group_halo_update_type
  use sw_core_tlm_mod,    only: c_sw, d_sw
  use sw_core_tlm_mod,    only: c_sw_tlm, d_sw_tlm
  use a2b_edge_tlm_mod,   only: a2b_ord2, a2b_ord4
  use a2b_edge_tlm_mod,   only: a2b_ord2_tlm, a2b_ord4_tlm
  use nh_core_tlm_mod,    only: Riem_Solver3, Riem_Solver_C, update_dz_c, update_dz_d, nest_halo_nh
  use nh_core_tlm_mod,    only: Riem_Solver3_tlm, Riem_Solver_C_tlm, update_dz_c_tlm, update_dz_d_tlm, nest_halo_nh_tlm
  use tp_core_tlm_mod,    only: copy_corners
  use tp_core_tlm_mod,    only: copy_corners_tlm
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

  use boundary_tlm_mod,   only: extrapolation_BC,  nested_grid_BC_apply_intT
  use boundary_tlm_mod,   only: nested_grid_BC_apply_intT_tlm

#ifdef SW_DYNAMICS
  use test_cases_mod,     only: test_case, case9_forcing1, case9_forcing2
#endif

  use fv_arrays_nlm_mod,  only: fv_flags_pert_type, fpp

implicit none
private

public :: dyn_core, del2_cubed, init_ijk_mem
public :: dyn_core_tlm, del2_cubed_tlm

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
!  Differentiation of dyn_core in forward (tangent) mode:
!   variations   of useful results: pk3 xfx ws peln q gz du u dv
!                v w delp ua uc ptc mfx delz mfy omga ut divgd
!                pkc delpc va vc yfx pkz pe vt pk zh pt cx cy crx
!                cry
!   with respect to varying inputs: pk3 xfx ws peln q gz du u dv
!                v w delp ua uc ptc delz omga ut divgd pkc delpc
!                va vc yfx pkz pe vt pk zh pt crx cry
!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
  SUBROUTINE DYN_CORE_TLM(npx, npy, npz, ng, sphum, nq, bdt, n_split, &
&   zvir, cp, akap, cappa, grav, hydrostatic, u, u_tl, v, v_tl, w, w_tl&
&   , delz, delz_tl, pt, pt_tl, q, q_tl, delp, delp_tl, pe, pe_tl, pk, &
&   pk_tl, phis, ws, ws_tl, omga, omga_tl, ptop, pfull, ua, ua_tl, va, &
&   va_tl, uc, uc_tl, vc, vc_tl, mfx, mfx_tl, mfy, mfy_tl, cx, cx_tl, cy&
&   , cy_tl, pkz, pkz_tl, peln, peln_tl, q_con, ak, bk, dpx, dpx_tl, ks&
&   , gridstruct, flagstruct, flagstructp, neststruct, idiag, bd, domain&
&   , init_step, i_pack, end_step, gz, gz_tl, pkc, pkc_tl, ptc, ptc_tl, &
&   crx, crx_tl, xfx, xfx_tl, cry, cry_tl, yfx, yfx_tl, divgd, divgd_tl&
&   , delpc, delpc_tl, ut, ut_tl, vt, vt_tl, zh, zh_tl, pk3, pk3_tl, du&
&   , du_tl, dv, dv_tl, time_total)
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
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u_tl
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v_tl
! vertical vel. (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: w_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! delta-height (m, negative)
    REAL, INTENT(INOUT) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delz_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! moist kappa
    REAL, INTENT(INOUT) :: cappa(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: pt_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
!
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
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
    REAL, INTENT(INOUT) :: pe_tl(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1&
&   )
! ln(pe)
    REAL, INTENT(INOUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: peln_tl(bd%is:bd%ie, npz+1, bd%js:bd%je)
! pe**kappa
    REAL, INTENT(INOUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL, INTENT(INOUT) :: pk_tl(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL(kind=8), INTENT(INOUT) :: dpx(bd%is:bd%ie, bd%js:bd%je)
    REAL(kind=8), INTENT(INOUT) :: dpx_tl(bd%is:bd%ie, bd%js:bd%je)
!-----------------------------------------------------------------------
! Others:
    REAL, PARAMETER :: near0=1.e-8
    REAL, PARAMETER :: huge_r=1.e8
!-----------------------------------------------------------------------
! w at surface
    REAL, INTENT(OUT) :: ws(bd%is:bd%ie, bd%js:bd%je)
    REAL, INTENT(OUT) :: ws_tl(bd%is:bd%ie, bd%js:bd%je)
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: omga_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! (uc, vc) are mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: uc_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: vc_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua_tl, va_tl
    REAL, INTENT(INOUT) :: q_con(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! The Flux capacitors: accumulated Mass flux arrays
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_tl(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_tl(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je, npz), INTENT(INOUT) :: pkz
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je, npz), INTENT(INOUT) :: &
&   pkz_tl
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    TYPE(FV_FLAGS_PERT_TYPE), INTENT(IN), TARGET :: flagstructp
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(FV_DIAG_TYPE), INTENT(IN) :: idiag
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
!real, allocatable, dimension(:,:,:):: pem, heat_source
    REAL :: pem(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1), heat_source(bd&
&   %isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: pem_tl(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1), &
&   heat_source_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Auto 1D & 2D arrays:
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ws3, z_rat
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ws3_tl, z_rat_tl
    REAL :: dp_ref(npz)
! surface height (m)
    REAL :: zs(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: p1d(bd%is:bd%ie)
    REAL :: om2d(bd%is:bd%ie, npz)
    REAL :: om2d_tl(bd%is:bd%ie, npz)
    REAL :: wbuffer(npy+2, npz)
    REAL :: ebuffer(npy+2, npz)
    REAL :: ebuffer_tl(npy+2, npz)
    REAL :: nbuffer(npx+2, npz)
    REAL :: nbuffer_tl(npx+2, npz)
    REAL :: sbuffer(npx+2, npz)
! ----   For external mode:
    REAL :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: divg2_tl(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: fz(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: heat_s(bd%is:bd%ie, bd%js:bd%je)
    REAL :: heat_s_tl(bd%is:bd%ie, bd%js:bd%je)
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
    REAL :: dtmp_tl
    LOGICAL :: last_step, remap_step
    LOGICAL :: used
    REAL :: split_timestep_bc
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pkc(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pkc_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: ptc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ptc_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: crx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: xfx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cry(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cry_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: yfx_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: divgd(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: divgd_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, &
&   npz)
    REAL, INTENT(INOUT) :: delpc(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delpc_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ut(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: ut_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vt_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: zh(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: zh_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk3(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk3_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: du_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dv_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
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
    REAL :: min1_tl
    REAL :: min2
    REAL :: min2_tl
    REAL :: abs0
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: arg2
    REAL :: arg2_tl
    REAL :: y1_tl
    REAL :: y2_tl
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
      mfx_tl = 0.0
      mfy_tl = 0.0
      cx_tl = 0.0
      cy_tl = 0.0
      om2d_tl = 0.0
      pem_tl = 0.0
      ws3_tl = 0.0
      z_rat_tl = 0.0
      heat_source_tl = 0.0
      heat_s_tl = 0.0
      wk_tl = 0.0
      divg2_tl = 0.0
    ELSE IF (flagstruct%d2_bg_k1 .LT. 1.e-3) THEN
      n_con = 0
      mfx_tl = 0.0
      mfy_tl = 0.0
      cx_tl = 0.0
      cy_tl = 0.0
      om2d_tl = 0.0
      pem_tl = 0.0
      ws3_tl = 0.0
      z_rat_tl = 0.0
      heat_source_tl = 0.0
      heat_s_tl = 0.0
      wk_tl = 0.0
      divg2_tl = 0.0
    ELSE IF (flagstruct%d2_bg_k2 .LT. 1.e-3) THEN
      n_con = 1
      mfx_tl = 0.0
      mfy_tl = 0.0
      cx_tl = 0.0
      cy_tl = 0.0
      om2d_tl = 0.0
      pem_tl = 0.0
      ws3_tl = 0.0
      z_rat_tl = 0.0
      heat_source_tl = 0.0
      heat_s_tl = 0.0
      wk_tl = 0.0
      divg2_tl = 0.0
    ELSE
      n_con = 2
      mfx_tl = 0.0
      mfy_tl = 0.0
      cx_tl = 0.0
      cy_tl = 0.0
      om2d_tl = 0.0
      pem_tl = 0.0
      ws3_tl = 0.0
      z_rat_tl = 0.0
      heat_source_tl = 0.0
      heat_s_tl = 0.0
      wk_tl = 0.0
      divg2_tl = 0.0
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
        IF (flagstruct%inline_q) CALL START_GROUP_HALO_UPDATE_TLM(i_pack&
&                                                           (10), q, &
&                                                           q_tl, domain&
&                                                          )
        CALL TIMING_OFF('COMM_TRACER')
        CALL TIMING_OFF('COMM_TOTAL')
      END IF
      IF (.NOT.hydrostatic) THEN
        CALL TIMING_ON('COMM_TOTAL')
        CALL START_GROUP_HALO_UPDATE_TLM(i_pack(7), w, w_tl, domain)
        CALL TIMING_OFF('COMM_TOTAL')
        IF (it .EQ. 1) THEN
          IF (gridstruct%nested) THEN
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,gz,zs,delz)
            DO j=jsd,jed
              DO i=isd,ied
                gz_tl(i, j, npz+1) = 0.0
                gz(i, j, npz+1) = zs(i, j)
              END DO
              DO k=npz,1,-1
                DO i=isd,ied
                  gz_tl(i, j, k) = gz_tl(i, j, k+1) - delz_tl(i, j, k)
                  gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)
                END DO
              END DO
            END DO
          ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gz,zs,delz)
            DO j=js,je
              DO i=is,ie
                gz_tl(i, j, npz+1) = 0.0
                gz(i, j, npz+1) = zs(i, j)
              END DO
              DO k=npz,1,-1
                DO i=is,ie
                  gz_tl(i, j, k) = gz_tl(i, j, k+1) - delz_tl(i, j, k)
                  gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)
                END DO
              END DO
            END DO
          END IF
          CALL TIMING_ON('COMM_TOTAL')
          CALL START_GROUP_HALO_UPDATE_TLM(i_pack(5), gz, gz_tl, domain)
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
          pem_tl = 0.0
!allocate ( pem(is-1:ie+1,npz+1,js-1:je+1) )
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pem,delp,ptop)
          DO j=js-1,je+1
            DO i=is-1,ie+1
              pem_tl(i, 1, j) = 0.0
              pem(i, 1, j) = ptop
            END DO
            DO k=1,npz
              DO i=is-1,ie+1
                pem_tl(i, k+1, j) = pem_tl(i, k, j) + delp_tl(i, j, k)
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
        CALL C_SW_TLM(delpc(isd:ied, jsd:jed, k), delpc_tl(isd:ied, jsd:&
&               jed, k), delp(isd:ied, jsd:jed, k), delp_tl(isd:ied, jsd&
&               :jed, k), ptc(isd:ied, jsd:jed, k), ptc_tl(isd:ied, jsd:&
&               jed, k), pt(isd:ied, jsd:jed, k), pt_tl(isd:ied, jsd:jed&
&               , k), u(isd:ied, jsd:jed+1, k), u_tl(isd:ied, jsd:jed+1&
&               , k), v(isd:ied+1, jsd:jed, k), v_tl(isd:ied+1, jsd:jed&
&               , k), w(isd:ied, jsd:jed, k), w_tl(isd:ied, jsd:jed, k)&
&               , uc(isd:ied+1, jsd:jed, k), uc_tl(isd:ied+1, jsd:jed, k&
&               ), vc(isd:ied, jsd:jed+1, k), vc_tl(isd:ied, jsd:jed+1, &
&               k), ua(isd:ied, jsd:jed, k), ua_tl(isd:ied, jsd:jed, k)&
&               , va(isd:ied, jsd:jed, k), va_tl(isd:ied, jsd:jed, k), &
&               omga(isd:ied, jsd:jed, k), omga_tl(isd:ied, jsd:jed, k)&
&               , ut(isd:ied, jsd:jed, k), ut_tl(isd:ied, jsd:jed, k), &
&               vt(isd:ied, jsd:jed, k), vt_tl(isd:ied, jsd:jed, k), &
&               divgd(isd:ied+1, jsd:jed+1, k), divgd_tl(isd:ied+1, jsd:&
&               jed+1, k), flagstruct%nord, dt2, hydrostatic, .true., bd&
&               , gridstruct, flagstruct)
      END DO
      CALL TIMING_OFF('c_sw')
      IF (flagstruct%nord .GT. 0) THEN
        CALL TIMING_ON('COMM_TOTAL')
        CALL START_GROUP_HALO_UPDATE_TLM(i_pack(3), divgd, divgd_tl, &
&                                  domain, position=corner)
        CALL TIMING_OFF('COMM_TOTAL')
      END IF
      IF (gridstruct%nested) THEN
        CALL NESTED_GRID_BC_APPLY_INTT_TLM(delpc, delpc_tl, 0, 0, npx, &
&                                    npy, npz, bd, split_timestep_bc + &
&                                    0.5, REAL(n_split*flagstruct%&
&                                    k_split), neststruct%delp_bc, &
&                                    bctype=neststruct%nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT_TLM(ptc, ptc_tl, 0, 0, npx, npy, &
&                                    npz, bd, split_timestep_bc + 0.5, &
&                                    REAL(n_split*flagstruct%k_split), &
&                                    neststruct%pt_bc, bctype=neststruct&
&                                    %nestbctype)
      END IF
! end hydro check
      IF (hydrostatic) THEN
        CALL GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, delpc, delpc_tl, &
&                pkc, pkc_tl, gz, gz_tl, phis, ptc, ptc_tl, q_con, pkz, &
&                pkz_tl, npz, akap, .true., gridstruct%nested, .false., &
&                npx, npy, flagstruct%a2b_ord, bd)
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
                zh_tl(i, j, k) = gz_tl(i, j, k)
                zh(i, j, k) = gz(i, j, k)
              END DO
            END DO
          END DO
        ELSE
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,zh,gz)
          DO k=1,npz+1
            DO j=jsd,jed
              DO i=isd,ied
                gz_tl(i, j, k) = zh_tl(i, j, k)
                gz(i, j, k) = zh(i, j, k)
              END DO
            END DO
          END DO
        END IF
        CALL TIMING_ON('UPDATE_DZ_C')
        CALL UPDATE_DZ_C_TLM(is, ie, js, je, npz, ng, dt2, dp_ref, zs, &
&                      gridstruct%area, ut, ut_tl, vt, vt_tl, gz, gz_tl&
&                      , ws3, ws3_tl, npx, npy, gridstruct%sw_corner, &
&                      gridstruct%se_corner, gridstruct%ne_corner, &
&                      gridstruct%nw_corner, bd, gridstruct%grid_type)
        CALL TIMING_OFF('UPDATE_DZ_C')
        CALL TIMING_ON('Riem_Solver')
        CALL RIEM_SOLVER_C_TLM(ms, dt2, is, ie, js, je, npz, ng, akap, &
&                        cappa, cp, ptop, phis, omga, omga_tl, ptc, &
&                        ptc_tl, q_con, delpc, delpc_tl, gz, gz_tl, pkc&
&                        , pkc_tl, ws3, ws3_tl, flagstruct%p_fac, &
&                        flagstruct%a_imp, flagstruct%scale_z)
        CALL TIMING_OFF('Riem_Solver')
        IF (gridstruct%nested) THEN
          CALL NESTED_GRID_BC_APPLY_INTT_TLM(delz, delz_tl, 0, 0, npx, &
&                                      npy, npz, bd, split_timestep_bc +&
&                                      0.5, REAL(n_split*flagstruct%&
&                                      k_split), neststruct%delz_bc, &
&                                      bctype=neststruct%nestbctype)
!Compute gz/pkc
!NOTE: nominally only need to compute quantities one out in the halo for p_grad_c
!(instead of entire halo)
          CALL NEST_HALO_NH_TLM(ptop, grav, akap, cp, delpc, delpc_tl, &
&                         delz, delz_tl, ptc, ptc_tl, phis, pkc, pkc_tl&
&                         , gz, gz_tl, pk3, pk3_tl, npx, npy, npz, &
&                         gridstruct%nested, .false., .false., .false., &
&                         bd)
        END IF
      END IF
      CALL P_GRAD_C_TLM(dt2, npz, delpc, delpc_tl, pkc, pkc_tl, gz, &
&                 gz_tl, uc, uc_tl, vc, vc_tl, bd, gridstruct%rdxc, &
&                 gridstruct%rdyc, hydrostatic)
      CALL TIMING_ON('COMM_TOTAL')
      CALL START_GROUP_HALO_UPDATE_TLM(i_pack(9), uc, uc_tl, vc, vc_tl, &
&                                domain, gridtype=cgrid_ne)
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
        CALL NESTED_GRID_BC_APPLY_INTT_TLM(vc, vc_tl, 0, 1, npx, npy, &
&                                    npz, bd, split_timestep_bc + 0.5, &
&                                    REAL(n_split*flagstruct%k_split), &
&                                    neststruct%vc_bc, bctype=neststruct&
&                                    %nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT_TLM(uc, uc_tl, 1, 0, npx, npy, &
&                                    npz, bd, split_timestep_bc + 0.5, &
&                                    REAL(n_split*flagstruct%k_split), &
&                                    neststruct%uc_bc, bctype=neststruct&
&                                    %nestbctype)
!QUESTION: What to do with divgd in nested halo?
        CALL NESTED_GRID_BC_APPLY_INTT_TLM(divgd, divgd_tl, 1, 1, npx, &
&                                    npy, npz, bd, split_timestep_bc, &
&                                    REAL(n_split*flagstruct%k_split), &
&                                    neststruct%divg_bc, bctype=&
&                                    neststruct%nestbctype)
!!$            if (is == 1 .and. js == 1) then
!!$               do j=jsd,5
!!$                  write(mpp_pe()+2000,*) j, divg(isd:5,j,1)
!!$            endif
      END IF
      IF (gridstruct%nested .AND. flagstruct%inline_q) THEN
        DO iq=1,nq
          CALL NESTED_GRID_BC_APPLY_INTT_TLM(q(isd:ied, jsd:jed, :, iq)&
&                                      , q_tl(isd:ied, jsd:jed, :, iq), &
&                                      0, 0, npx, npy, npz, bd, &
&                                      split_timestep_bc + 1, REAL(&
&                                      n_split*flagstruct%k_split), &
&                                      neststruct%q_bc(iq), bctype=&
&                                      neststruct%nestbctype)
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
              omga_tl(i, j, k) = delp_tl(i, j, k)
              omga(i, j, k) = delp(i, j, k)
            END DO
          END DO
        END IF
!--- external mode divergence damping ---
        IF (flagstruct%d_ext .GT. 0.) CALL A2B_ORD2_TLM(delp(isd:ied, &
&                                                 jsd:jed, k), delp_tl(&
&                                                 isd:ied, jsd:jed, k), &
&                                                 wk, wk_tl, gridstruct&
&                                                 , npx, npy, is, ie, js&
&                                                 , je, ng, .false.)
        IF (.NOT.hydrostatic .AND. flagstruct%do_f3d) THEN
! Correction factor for 3D Coriolis force
          DO j=jsd,jed
            DO i=isd,ied
              z_rat_tl(i, j) = (zh_tl(i, j, k)+zh_tl(i, j, k+1))/radius
              z_rat(i, j) = 1. + (zh(i, j, k)+zh(i, j, k+1))/radius
            END DO
          END DO
        END IF
        CALL D_SW_TLM(vt(isd:ied, jsd:jed, k), vt_tl(isd:ied, jsd:jed, k&
&               ), delp(isd:ied, jsd:jed, k), delp_tl(isd:ied, jsd:jed, &
&               k), ptc(isd:ied, jsd:jed, k), ptc_tl(isd:ied, jsd:jed, k&
&               ), pt(isd:ied, jsd:jed, k), pt_tl(isd:ied, jsd:jed, k), &
&               u(isd:ied, jsd:jed+1, k), u_tl(isd:ied, jsd:jed+1, k), v&
&               (isd:ied+1, jsd:jed, k), v_tl(isd:ied+1, jsd:jed, k), w(&
&               isd:ied, jsd:jed, k), w_tl(isd:ied, jsd:jed, k), uc(isd:&
&               ied+1, jsd:jed, k), uc_tl(isd:ied+1, jsd:jed, k), vc(isd&
&               :ied, jsd:jed+1, k), vc_tl(isd:ied, jsd:jed+1, k), ua(&
&               isd:ied, jsd:jed, k), ua_tl(isd:ied, jsd:jed, k), va(isd&
&               :ied, jsd:jed, k), va_tl(isd:ied, jsd:jed, k), divgd(isd&
&               :ied+1, jsd:jed+1, k), divgd_tl(isd:ied+1, jsd:jed+1, k)&
&               , mfx(is:ie+1, js:je, k), mfx_tl(is:ie+1, js:je, k), mfy&
&               (is:ie, js:je+1, k), mfy_tl(is:ie, js:je+1, k), cx(is:ie&
&               +1, jsd:jed, k), cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied&
&               , js:je+1, k), cy_tl(isd:ied, js:je+1, k), crx(is:ie+1, &
&               jsd:jed, k), crx_tl(is:ie+1, jsd:jed, k), cry(isd:ied, &
&               js:je+1, k), cry_tl(isd:ied, js:je+1, k), xfx(is:ie+1, &
&               jsd:jed, k), xfx_tl(is:ie+1, jsd:jed, k), yfx(isd:ied, &
&               js:je+1, k), yfx_tl(isd:ied, js:je+1, k), q_con(isd:ied&
&               , jsd:jed, 1), z_rat(isd:ied, jsd:jed), z_rat_tl(isd:ied&
&               , jsd:jed), kgb, heat_s, heat_s_tl, dpx, dpx_tl, zvir, &
&               sphum, nq, q, q_tl, k, npz, flagstruct%inline_q, dt, &
&               flagstruct%hord_tr, hord_m, hord_v, hord_t, hord_p, &
&               nord_k, nord_v(k), nord_w, nord_t, flagstruct%dddmp, &
&               d2_divg, flagstruct%d4_bg, damp_vt(k), damp_w, damp_t, &
&               d_con_k, hydrostatic, gridstruct, flagstruct, bd, &
&               flagstructp%hord_tr_pert, hord_m_pert, hord_v_pert, &
&               hord_t_pert, hord_p_pert, flagstructp%split_damp, &
&               nord_k_pert, nord_v_pert(k), nord_w_pert, nord_t_pert, &
&               flagstructp%dddmp_pert, d2_divg_pert, flagstructp%&
&               d4_bg_pert, damp_vt_pert(k), damp_w_pert, damp_t_pert)
        IF (hydrostatic .AND. (.NOT.flagstruct%use_old_omega) .AND. &
&           last_step) THEN
! Average horizontal "convergence" to cell center
          DO j=js,je
            DO i=is,ie
              omga_tl(i, j, k) = gridstruct%rarea(i, j)*rdt*(omga_tl(i, &
&               j, k)*(xfx(i, j, k)-xfx(i+1, j, k)+yfx(i, j, k)-yfx(i, j&
&               +1, k))+omga(i, j, k)*(xfx_tl(i, j, k)-xfx_tl(i+1, j, k)&
&               +yfx_tl(i, j, k)-yfx_tl(i, j+1, k)))
              omga(i, j, k) = omga(i, j, k)*(xfx(i, j, k)-xfx(i+1, j, k)&
&               +yfx(i, j, k)-yfx(i, j+1, k))*gridstruct%rarea(i, j)*rdt
            END DO
          END DO
        END IF
        IF (flagstruct%d_ext .GT. 0.) THEN
          DO j=js,jep1
            DO i=is,iep1
! delp at cell corners
              ptc_tl(i, j, k) = wk_tl(i, j)
              ptc(i, j, k) = wk(i, j)
            END DO
          END DO
        END IF
        IF (flagstruct%d_con .GT. 1.0e-5) THEN
! Average horizontal "convergence" to cell center
          DO j=js,je
            DO i=is,ie
              heat_source_tl(i, j, k) = heat_source_tl(i, j, k) + &
&               heat_s_tl(i, j)
              heat_source(i, j, k) = heat_source(i, j, k) + heat_s(i, j)
            END DO
          END DO
        END IF
      END DO
! end openMP k-loop
      CALL TIMING_OFF('d_sw')
      IF (flagstruct%fill_dp) CALL MIX_DP_TLM(hydrostatic, w, w_tl, delp&
&                                       , delp_tl, pt, pt_tl, npz, ak, &
&                                       bk, .false., flagstruct%fv_debug&
&                                       , bd)
      CALL TIMING_ON('COMM_TOTAL')
      CALL START_GROUP_HALO_UPDATE_TLM(i_pack(1), delp, delp_tl, domain&
&                                , complete=.true.)
      CALL START_GROUP_HALO_UPDATE_TLM(i_pack(1), pt, pt_tl, domain, &
&                                complete=.true.)
      CALL TIMING_OFF('COMM_TOTAL')
      IF (flagstruct%d_ext .GT. 0.) THEN
        d2_divg = flagstruct%d_ext*gridstruct%da_min_c
!$OMP parallel do default(none) shared(is,iep1,js,jep1,npz,wk,ptc,divg2,vt,d2_divg)
        DO j=js,jep1
          DO i=is,iep1
            wk_tl(i, j) = ptc_tl(i, j, 1)
            wk(i, j) = ptc(i, j, 1)
            divg2_tl(i, j) = wk_tl(i, j)*vt(i, j, 1) + wk(i, j)*vt_tl(i&
&             , j, 1)
            divg2(i, j) = wk(i, j)*vt(i, j, 1)
          END DO
          DO k=2,npz
            DO i=is,iep1
              wk_tl(i, j) = wk_tl(i, j) + ptc_tl(i, j, k)
              wk(i, j) = wk(i, j) + ptc(i, j, k)
              divg2_tl(i, j) = divg2_tl(i, j) + ptc_tl(i, j, k)*vt(i, j&
&               , k) + ptc(i, j, k)*vt_tl(i, j, k)
              divg2(i, j) = divg2(i, j) + ptc(i, j, k)*vt(i, j, k)
            END DO
          END DO
          DO i=is,iep1
            divg2_tl(i, j) = (d2_divg*divg2_tl(i, j)*wk(i, j)-d2_divg*&
&             divg2(i, j)*wk_tl(i, j))/wk(i, j)**2
            divg2(i, j) = d2_divg*divg2(i, j)/wk(i, j)
          END DO
        END DO
      ELSE
        divg2(:, :) = 0.
        divg2_tl = 0.0
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
        CALL NESTED_GRID_BC_APPLY_INTT_TLM(delp, delp_tl, 0, 0, npx, npy&
&                                    , npz, bd, split_timestep_bc + 1, &
&                                    REAL(n_split*flagstruct%k_split), &
&                                    neststruct%delp_bc, bctype=&
&                                    neststruct%nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT_TLM(pt, pt_tl, 0, 0, npx, npy, &
&                                    npz, bd, split_timestep_bc + 1, &
&                                    REAL(n_split*flagstruct%k_split), &
&                                    neststruct%pt_bc, bctype=neststruct&
&                                    %nestbctype)
      END IF
! end hydro check
      IF (hydrostatic) THEN
        CALL GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, delp, delp_tl, &
&                pkc, pkc_tl, gz, gz_tl, phis, pt, pt_tl, q_con, pkz, &
&                pkz_tl, npz, akap, .false., gridstruct%nested, .true., &
&                npx, npy, flagstruct%a2b_ord, bd)
      ELSE
        CALL TIMING_ON('UPDATE_DZ')
        CALL UPDATE_DZ_D_TLM(nord_v, damp_vt, flagstruct%hord_tm, is, ie&
&                      , js, je, npz, ng, npx, npy, gridstruct%area, &
&                      gridstruct%rarea, dp_ref, zs, zh, zh_tl, crx, &
&                      crx_tl, cry, cry_tl, xfx, xfx_tl, yfx, yfx_tl, &
&                      delz, ws, ws_tl, rdt, gridstruct, bd, flagstructp&
&                      %hord_tm_pert)
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
        CALL RIEM_SOLVER3_TLM(flagstruct%m_split, dt, is, ie, js, je, &
&                       npz, ng, isd, ied, jsd, jed, akap, cappa, cp, &
&                       ptop, zs, q_con, w, w_tl, delz, delz_tl, pt, &
&                       pt_tl, delp, delp_tl, zh, zh_tl, pe, pe_tl, pkc&
&                       , pkc_tl, pk3, pk3_tl, pk, pk_tl, peln, peln_tl&
&                       , ws, ws_tl, flagstruct%scale_z, flagstruct%&
&                       p_fac, flagstruct%a_imp, flagstruct%use_logp, &
&                       remap_step, beta .LT. -0.1)
        CALL TIMING_OFF('Riem_Solver')
        CALL TIMING_ON('COMM_TOTAL')
        IF (gridstruct%square_domain) THEN
          CALL START_GROUP_HALO_UPDATE_TLM(i_pack(4), zh, zh_tl, domain)
          CALL START_GROUP_HALO_UPDATE_TLM(i_pack(5), pkc, pkc_tl, &
&                                    domain, whalo=2, ehalo=2, shalo=2, &
&                                    nhalo=2)
        ELSE
          CALL START_GROUP_HALO_UPDATE_TLM(i_pack(4), zh, zh_tl, domain&
&                                    , complete=.true.)
          CALL START_GROUP_HALO_UPDATE_TLM(i_pack(4), pkc, pkc_tl, &
&                                    domain, complete=.true.)
        END IF
        CALL TIMING_OFF('COMM_TOTAL')
        IF (remap_step) CALL PE_HALO_TLM(is, ie, js, je, isd, ied, jsd, &
&                                  jed, npz, ptop, pe, pe_tl, delp, &
&                                  delp_tl)
        IF (flagstruct%use_logp) THEN
          CALL PLN_HALO_TLM(is, ie, js, je, isd, ied, jsd, jed, npz, &
&                     ptop, pk3, pk3_tl, delp, delp_tl)
        ELSE
          CALL PK3_HALO_TLM(is, ie, js, je, isd, ied, jsd, jed, npz, &
&                     ptop, akap, pk3, pk3_tl, delp, delp_tl)
        END IF
        IF (gridstruct%nested) THEN
          CALL NESTED_GRID_BC_APPLY_INTT_TLM(delz, delz_tl, 0, 0, npx, &
&                                      npy, npz, bd, split_timestep_bc +&
&                                      1., REAL(n_split*flagstruct%&
&                                      k_split), neststruct%delz_bc, &
&                                      bctype=neststruct%nestbctype)
!Compute gz/pkc/pk3; note that now pkc should be nonhydro pert'n pressure
          CALL NEST_HALO_NH_TLM(ptop, grav, akap, cp, delp, delp_tl, &
&                         delz, delz_tl, pt, pt_tl, phis, pkc, pkc_tl, &
&                         gz, gz_tl, pk3, pk3_tl, npx, npy, npz, &
&                         gridstruct%nested, .true., .true., .true., bd)
        END IF
        CALL TIMING_ON('COMM_TOTAL')
        CALL COMPLETE_GROUP_HALO_UPDATE(i_pack(4), domain)
        CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gz,zh,grav)
        DO k=1,npz+1
          DO j=js-2,je+2
            DO i=is-2,ie+2
              gz_tl(i, j, k) = grav*zh_tl(i, j, k)
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
              pk_tl(i, j, k) = pkc_tl(i, j, k)
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
          CALL GRAD1_P_UPDATE_TLM(divg2, divg2_tl, u, u_tl, v, v_tl, pkc&
&                           , pkc_tl, gz, gz_tl, du, du_tl, dv, dv_tl, &
&                           dt, ng, gridstruct, bd, npx, npy, npz, ptop&
&                           , beta_d, flagstruct%a2b_ord)
        ELSE
          CALL ONE_GRAD_P_TLM(u, u_tl, v, v_tl, pkc, pkc_tl, gz, gz_tl, &
&                       divg2, divg2_tl, delp, delp_tl, dt, ng, &
&                       gridstruct, bd, npx, npy, npz, ptop, hydrostatic&
&                       , flagstruct%a2b_ord, flagstruct%d_ext)
        END IF
      ELSE IF (beta .GT. 0.) THEN
        CALL SPLIT_P_GRAD_TLM(u, u_tl, v, v_tl, pkc, pkc_tl, gz, gz_tl, &
&                       du, du_tl, dv, dv_tl, delp, delp_tl, pk3, pk3_tl&
&                       , beta_d, dt, ng, gridstruct, bd, npx, npy, npz&
&                       , flagstruct%use_logp)
      ELSE IF (beta .LT. -0.1) THEN
        CALL ONE_GRAD_P_TLM(u, u_tl, v, v_tl, pkc, pkc_tl, gz, gz_tl, &
&                     divg2, divg2_tl, delp, delp_tl, dt, ng, gridstruct&
&                     , bd, npx, npy, npz, ptop, hydrostatic, flagstruct&
&                     %a2b_ord, flagstruct%d_ext)
      ELSE
        CALL NH_P_GRAD_TLM(u, u_tl, v, v_tl, pkc, pkc_tl, gz, gz_tl, &
&                    delp, delp_tl, pk3, pk3_tl, dt, ng, gridstruct, bd&
&                    , npx, npy, npz, flagstruct%use_logp)
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
                arg1_tl = (rdg*delp_tl(i, j, k)*delz(i, j, k)-rdg*delp(i&
&                 , j, k)*delz_tl(i, j, k))*pt(i, j, k)/delz(i, j, k)**2&
&                 + rdg*delp(i, j, k)*pt_tl(i, j, k)/delz(i, j, k)
                arg1 = rdg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
                arg2_tl = k1k*arg1_tl/arg1
                arg2 = k1k*LOG(arg1)
                pkz_tl(i, j, k) = arg2_tl*EXP(arg2)
                pkz(i, j, k) = EXP(arg2)
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
        CALL MPP_GET_BOUNDARY_TLM(u, u_tl, v, v_tl, domain, ebuffery=&
&                           ebuffer, ebuffery_tl=ebuffer_tl, nbufferx=&
&                           nbuffer, nbufferx_tl=nbuffer_tl, gridtype=&
&                           dgrid_ne)
!$OMP parallel do default(none) shared(is,ie,js,je,npz,u,nbuffer,v,ebuffer)
        DO k=1,npz
          DO i=is,ie
            u_tl(i, je+1, k) = nbuffer_tl(i-is+1, k)
            u(i, je+1, k) = nbuffer(i-is+1, k)
          END DO
          DO j=js,je
            v_tl(ie+1, j, k) = ebuffer_tl(j-js+1, k)
            v(ie+1, j, k) = ebuffer(j-js+1, k)
          END DO
        END DO
      END IF
      IF (it .NE. n_split) CALL START_GROUP_HALO_UPDATE_TLM(i_pack(8), u&
&                                                     , u_tl, v, v_tl, &
&                                                     domain, gridtype=&
&                                                     dgrid_ne)
      CALL TIMING_OFF('COMM_TOTAL')
      IF (gridstruct%nested) neststruct%nest_timestep = neststruct%&
&         nest_timestep + 1
      IF (hydrostatic .AND. last_step) THEN
        IF (flagstruct%use_old_omega) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga,pe,pem,rdt)
          DO k=1,npz
            DO j=js,je
              DO i=is,ie
                omga_tl(i, j, k) = rdt*(pe_tl(i, k+1, j)-pem_tl(i, k+1, &
&                 j))
                omga(i, j, k) = (pe(i, k+1, j)-pem(i, k+1, j))*rdt
              END DO
            END DO
          END DO
!------------------------------
! Compute the "advective term"
!------------------------------
          CALL ADV_PE_TLM(ua, ua_tl, va, va_tl, pem, pem_tl, omga, &
&                   omga_tl, gridstruct, bd, npx, npy, npz, ng)
        ELSE
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga) private(om2d)
          DO j=js,je
            DO k=1,npz
              DO i=is,ie
                om2d_tl(i, k) = omga_tl(i, j, k)
                om2d(i, k) = omga(i, j, k)
              END DO
            END DO
            DO k=2,npz
              DO i=is,ie
                om2d_tl(i, k) = om2d_tl(i, k-1) + omga_tl(i, j, k)
                om2d(i, k) = om2d(i, k-1) + omga(i, j, k)
              END DO
            END DO
            DO k=2,npz
              DO i=is,ie
                omga_tl(i, j, k) = om2d_tl(i, k)
                omga(i, j, k) = om2d(i, k)
              END DO
            END DO
          END DO
        END IF
        IF (idiag%id_ws .GT. 0 .AND. hydrostatic) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ws,delz,delp,omga)
          DO j=js,je
            DO i=is,ie
              ws_tl(i, j) = (delz_tl(i, j, npz)*delp(i, j, npz)-delz(i, &
&               j, npz)*delp_tl(i, j, npz))*omga(i, j, npz)/delp(i, j, &
&               npz)**2 + delz(i, j, npz)*omga_tl(i, j, npz)/delp(i, j, &
&               npz)
              ws(i, j) = delz(i, j, npz)/delp(i, j, npz)*omga(i, j, npz)
            END DO
          END DO
          used = SEND_DATA(idiag%id_ws, ws, fv_time)
        END IF
      END IF
      IF (gridstruct%nested) THEN
        IF (.NOT.hydrostatic) CALL NESTED_GRID_BC_APPLY_INTT_TLM(w, w_tl&
&                                                          , 0, 0, npx, &
&                                                          npy, npz, bd&
&                                                          , &
&                                                      split_timestep_bc&
&                                                          + 1, REAL(&
&                                                          n_split*&
&                                                          flagstruct%&
&                                                          k_split), &
&                                                          neststruct%&
&                                                          w_bc, bctype=&
&                                                          neststruct%&
&                                                          nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT_TLM(u, u_tl, 0, 1, npx, npy, npz&
&                                    , bd, split_timestep_bc + 1, REAL(&
&                                    n_split*flagstruct%k_split), &
&                                    neststruct%u_bc, bctype=neststruct%&
&                                    nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT_TLM(v, v_tl, 1, 0, npx, npy, npz&
&                                    , bd, split_timestep_bc + 1, REAL(&
&                                    n_split*flagstruct%k_split), &
&                                    neststruct%v_bc, bctype=neststruct%&
&                                    nestbctype)
      END IF
    END DO
!-----------------------------------------------------
! time split loop
!-----------------------------------------------------
    IF (nq .GT. 0 .AND. (.NOT.flagstruct%inline_q)) THEN
      CALL TIMING_ON('COMM_TOTAL')
      CALL TIMING_ON('COMM_TRACER')
      CALL START_GROUP_HALO_UPDATE_TLM(i_pack(10), q, q_tl, domain)
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
      CALL DEL2_CUBED_TLM(heat_source, heat_source_tl, cnst_0p20*&
&                   gridstruct%da_min, gridstruct, domain, npx, npy, npz&
&                   , nf_ke, bd)
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
                pt_tl(i, j, k) = pt_tl(i, j, k) + (heat_source_tl(i, j, &
&                 k)*cp_air*delp(i, j, k)*pkz(i, j, k)-heat_source(i, j&
&                 , k)*cp_air*(delp_tl(i, j, k)*pkz(i, j, k)+delp(i, j, &
&                 k)*pkz_tl(i, j, k)))/(cp_air*delp(i, j, k)*pkz(i, j, k&
&                 ))**2
                pt(i, j, k) = pt(i, j, k) + heat_source(i, j, k)/(cp_air&
&                 *delp(i, j, k)*pkz(i, j, k))
              END DO
            ELSE
              DO i=is,ie
                dtmp_tl = (heat_source_tl(i, j, k)*cp_air*delp(i, j, k)-&
&                 heat_source(i, j, k)*cp_air*delp_tl(i, j, k))/(cp_air*&
&                 delp(i, j, k))**2
                dtmp = heat_source(i, j, k)/(cp_air*delp(i, j, k))
                IF (bdt .GE. 0.) THEN
                  abs0 = bdt
                ELSE
                  abs0 = -bdt
                END IF
                x1 = abs0*flagstruct%delt_max
                IF (dtmp .GE. 0.) THEN
                  y1_tl = dtmp_tl
                  y1 = dtmp
                ELSE
                  y1_tl = -dtmp_tl
                  y1 = -dtmp
                END IF
                IF (x1 .GT. y1) THEN
                  min1_tl = y1_tl
                  min1 = y1
                ELSE
                  min1 = x1
                  min1_tl = 0.0
                END IF
                pt_tl(i, j, k) = pt_tl(i, j, k) + (min1_tl*SIGN(1.d0, &
&                 min1*dtmp)*pkz(i, j, k)-SIGN(min1, dtmp)*pkz_tl(i, j, &
&                 k))/pkz(i, j, k)**2
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
              arg1_tl = (rdg*delp_tl(i, j, k)*delz(i, j, k)-rdg*delp(i, &
&               j, k)*delz_tl(i, j, k))*pt(i, j, k)/delz(i, j, k)**2 + &
&               rdg*delp(i, j, k)*pt_tl(i, j, k)/delz(i, j, k)
              arg1 = rdg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
              arg2_tl = k1k*arg1_tl/arg1
              arg2 = k1k*LOG(arg1)
              pkz_tl(i, j, k) = arg2_tl*EXP(arg2)
              pkz(i, j, k) = EXP(arg2)
              dtmp_tl = (heat_source_tl(i, j, k)*cv_air*delp(i, j, k)-&
&               heat_source(i, j, k)*cv_air*delp_tl(i, j, k))/(cv_air*&
&               delp(i, j, k))**2
              dtmp = heat_source(i, j, k)/(cv_air*delp(i, j, k))
              IF (dtmp .GE. 0.) THEN
                y2_tl = dtmp_tl
                y2 = dtmp
              ELSE
                y2_tl = -dtmp_tl
                y2 = -dtmp
              END IF
              IF (delt .GT. y2) THEN
                min2_tl = y2_tl
                min2 = y2
              ELSE
                min2 = delt
                min2_tl = 0.0
              END IF
              pt_tl(i, j, k) = pt_tl(i, j, k) + (min2_tl*SIGN(1.d0, min2&
&               *dtmp)*pkz(i, j, k)-SIGN(min2, dtmp)*pkz_tl(i, j, k))/&
&               pkz(i, j, k)**2
              pt(i, j, k) = pt(i, j, k) + SIGN(min2, dtmp)/pkz(i, j, k)
            END DO
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE DYN_CORE_TLM
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
        CALL NESTED_GRID_BC_APPLY_INTT(delpc, 0, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 0.5, REAL(n_split*&
&                                flagstruct%k_split), neststruct%delp_bc&
&                                , neststruct%nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT(ptc, 0, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 0.5, REAL(n_split*&
&                                flagstruct%k_split), neststruct%pt_bc, &
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
          CALL NESTED_GRID_BC_APPLY_INTT(delz, 0, 0, npx, npy, npz, bd, &
&                                  split_timestep_bc + 0.5, REAL(n_split&
&                                  *flagstruct%k_split), neststruct%&
&                                  delz_bc, neststruct%nestbctype)
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
        CALL NESTED_GRID_BC_APPLY_INTT(vc, 0, 1, npx, npy, npz, bd, &
&                                split_timestep_bc + 0.5, REAL(n_split*&
&                                flagstruct%k_split), neststruct%vc_bc, &
&                                neststruct%nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT(uc, 1, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 0.5, REAL(n_split*&
&                                flagstruct%k_split), neststruct%uc_bc, &
&                                neststruct%nestbctype)
!QUESTION: What to do with divgd in nested halo?
        CALL NESTED_GRID_BC_APPLY_INTT(divgd, 1, 1, npx, npy, npz, bd, &
&                                split_timestep_bc, REAL(n_split*&
&                                flagstruct%k_split), neststruct%divg_bc&
&                                , neststruct%nestbctype)
!!$            if (is == 1 .and. js == 1) then
!!$               do j=jsd,5
!!$                  write(mpp_pe()+2000,*) j, divg(isd:5,j,1)
!!$            endif
      END IF
      IF (gridstruct%nested .AND. flagstruct%inline_q) THEN
        DO iq=1,nq
          CALL NESTED_GRID_BC_APPLY_INTT(q(isd:ied, jsd:jed, :, iq), 0, &
&                                  0, npx, npy, npz, bd, &
&                                  split_timestep_bc + 1, REAL(n_split*&
&                                  flagstruct%k_split), neststruct%q_bc(&
&                                  iq), neststruct%nestbctype)
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
        CALL NESTED_GRID_BC_APPLY_INTT(delp, 0, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 1, REAL(n_split*&
&                                flagstruct%k_split), neststruct%delp_bc&
&                                , neststruct%nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT(pt, 0, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 1, REAL(n_split*&
&                                flagstruct%k_split), neststruct%pt_bc, &
&                                neststruct%nestbctype)
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
        CALL RIEM_SOLVER3(flagstruct%m_split, dt, is, ie, js, je, npz, &
&                   ng, isd, ied, jsd, jed, akap, cappa, cp, ptop, zs, &
&                   q_con, w, delz, pt, delp, zh, pe, pkc, pk3, pk, peln&
&                   , ws, flagstruct%scale_z, flagstruct%p_fac, &
&                   flagstruct%a_imp, flagstruct%use_logp, remap_step, &
&                   beta .LT. -0.1)
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
          CALL NESTED_GRID_BC_APPLY_INTT(delz, 0, 0, npx, npy, npz, bd, &
&                                  split_timestep_bc + 1., REAL(n_split*&
&                                  flagstruct%k_split), neststruct%&
&                                  delz_bc, neststruct%nestbctype)
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
                arg1 = rdg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
                arg2 = k1k*LOG(arg1)
                pkz(i, j, k) = EXP(arg2)
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
        IF (.NOT.hydrostatic) CALL NESTED_GRID_BC_APPLY_INTT(w, 0, 0, &
&                                                      npx, npy, npz, bd&
&                                                      , &
&                                                      split_timestep_bc&
&                                                      + 1, REAL(n_split&
&                                                      *flagstruct%&
&                                                      k_split), &
&                                                      neststruct%w_bc, &
&                                                      neststruct%&
&                                                      nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT(u, 0, 1, npx, npy, npz, bd, &
&                                split_timestep_bc + 1, REAL(n_split*&
&                                flagstruct%k_split), neststruct%u_bc, &
&                                neststruct%nestbctype)
        CALL NESTED_GRID_BC_APPLY_INTT(v, 1, 0, npx, npy, npz, bd, &
&                                split_timestep_bc + 1, REAL(n_split*&
&                                flagstruct%k_split), neststruct%v_bc, &
&                                neststruct%nestbctype)
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
      CALL DEL2_CUBED(heat_source, cnst_0p20*gridstruct%da_min, &
&               gridstruct, domain, npx, npy, npz, nf_ke, bd)
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
              arg1 = rdg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
              arg2 = k1k*LOG(arg1)
              pkz(i, j, k) = EXP(arg2)
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
!  Differentiation of pk3_halo in forward (tangent) mode:
!   variations   of useful results: pk3
!   with respect to varying inputs: pk3 delp
  SUBROUTINE PK3_HALO_TLM(is, ie, js, je, isd, ied, jsd, jed, npz, ptop&
&   , akap, pk3, pk3_tl, delp, delp_tl)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop, akap
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp_tl
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3_tl
! Local:
    REAL :: pei(isd:ied)
    REAL :: pei_tl(isd:ied)
    REAL :: pej(jsd:jed)
    REAL :: pej_tl(jsd:jed)
    INTEGER :: i, j, k
    INTRINSIC LOG
    INTRINSIC EXP
    REAL :: arg1
    REAL :: arg1_tl
    pei_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,npz,ptop,delp,pk3,akap) &
!$OMP                          private(pei)
    DO j=js,je
      pei_tl(is-2) = 0.0
      pei(is-2) = ptop
      pei_tl(is-1) = 0.0
      pei(is-1) = ptop
      DO k=1,npz
        pei_tl(is-2) = pei_tl(is-2) + delp_tl(is-2, j, k)
        pei(is-2) = pei(is-2) + delp(is-2, j, k)
        pei_tl(is-1) = pei_tl(is-1) + delp_tl(is-1, j, k)
        pei(is-1) = pei(is-1) + delp(is-1, j, k)
        arg1_tl = akap*pei_tl(is-2)/pei(is-2)
        arg1 = akap*LOG(pei(is-2))
        pk3_tl(is-2, j, k+1) = arg1_tl*EXP(arg1)
        pk3(is-2, j, k+1) = EXP(arg1)
        arg1_tl = akap*pei_tl(is-1)/pei(is-1)
        arg1 = akap*LOG(pei(is-1))
        pk3_tl(is-1, j, k+1) = arg1_tl*EXP(arg1)
        pk3(is-1, j, k+1) = EXP(arg1)
      END DO
      pei_tl(ie+1) = 0.0
      pei(ie+1) = ptop
      pei_tl(ie+2) = 0.0
      pei(ie+2) = ptop
      DO k=1,npz
        pei_tl(ie+1) = pei_tl(ie+1) + delp_tl(ie+1, j, k)
        pei(ie+1) = pei(ie+1) + delp(ie+1, j, k)
        pei_tl(ie+2) = pei_tl(ie+2) + delp_tl(ie+2, j, k)
        pei(ie+2) = pei(ie+2) + delp(ie+2, j, k)
        arg1_tl = akap*pei_tl(ie+1)/pei(ie+1)
        arg1 = akap*LOG(pei(ie+1))
        pk3_tl(ie+1, j, k+1) = arg1_tl*EXP(arg1)
        pk3(ie+1, j, k+1) = EXP(arg1)
        arg1_tl = akap*pei_tl(ie+2)/pei(ie+2)
        arg1 = akap*LOG(pei(ie+2))
        pk3_tl(ie+2, j, k+1) = arg1_tl*EXP(arg1)
        pk3(ie+2, j, k+1) = EXP(arg1)
      END DO
    END DO
    pej_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ptop,delp,pk3,akap) &
!$OMP                          private(pej)
    DO i=is-2,ie+2
      pej_tl(js-2) = 0.0
      pej(js-2) = ptop
      pej_tl(js-1) = 0.0
      pej(js-1) = ptop
      DO k=1,npz
        pej_tl(js-2) = pej_tl(js-2) + delp_tl(i, js-2, k)
        pej(js-2) = pej(js-2) + delp(i, js-2, k)
        pej_tl(js-1) = pej_tl(js-1) + delp_tl(i, js-1, k)
        pej(js-1) = pej(js-1) + delp(i, js-1, k)
        arg1_tl = akap*pej_tl(js-2)/pej(js-2)
        arg1 = akap*LOG(pej(js-2))
        pk3_tl(i, js-2, k+1) = arg1_tl*EXP(arg1)
        pk3(i, js-2, k+1) = EXP(arg1)
        arg1_tl = akap*pej_tl(js-1)/pej(js-1)
        arg1 = akap*LOG(pej(js-1))
        pk3_tl(i, js-1, k+1) = arg1_tl*EXP(arg1)
        pk3(i, js-1, k+1) = EXP(arg1)
      END DO
      pej_tl(je+1) = 0.0
      pej(je+1) = ptop
      pej_tl(je+2) = 0.0
      pej(je+2) = ptop
      DO k=1,npz
        pej_tl(je+1) = pej_tl(je+1) + delp_tl(i, je+1, k)
        pej(je+1) = pej(je+1) + delp(i, je+1, k)
        pej_tl(je+2) = pej_tl(je+2) + delp_tl(i, je+2, k)
        pej(je+2) = pej(je+2) + delp(i, je+2, k)
        arg1_tl = akap*pej_tl(je+1)/pej(je+1)
        arg1 = akap*LOG(pej(je+1))
        pk3_tl(i, je+1, k+1) = arg1_tl*EXP(arg1)
        pk3(i, je+1, k+1) = EXP(arg1)
        arg1_tl = akap*pej_tl(je+2)/pej(je+2)
        arg1 = akap*LOG(pej(je+2))
        pk3_tl(i, je+2, k+1) = arg1_tl*EXP(arg1)
        pk3(i, je+2, k+1) = EXP(arg1)
      END DO
    END DO
  END SUBROUTINE PK3_HALO_TLM
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
    REAL :: arg1
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,npz,ptop,delp,pk3,akap) &
!$OMP                          private(pei)
    DO j=js,je
      pei(is-2) = ptop
      pei(is-1) = ptop
      DO k=1,npz
        pei(is-2) = pei(is-2) + delp(is-2, j, k)
        pei(is-1) = pei(is-1) + delp(is-1, j, k)
        arg1 = akap*LOG(pei(is-2))
        pk3(is-2, j, k+1) = EXP(arg1)
        arg1 = akap*LOG(pei(is-1))
        pk3(is-1, j, k+1) = EXP(arg1)
      END DO
      pei(ie+1) = ptop
      pei(ie+2) = ptop
      DO k=1,npz
        pei(ie+1) = pei(ie+1) + delp(ie+1, j, k)
        pei(ie+2) = pei(ie+2) + delp(ie+2, j, k)
        arg1 = akap*LOG(pei(ie+1))
        pk3(ie+1, j, k+1) = EXP(arg1)
        arg1 = akap*LOG(pei(ie+2))
        pk3(ie+2, j, k+1) = EXP(arg1)
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
        arg1 = akap*LOG(pej(js-2))
        pk3(i, js-2, k+1) = EXP(arg1)
        arg1 = akap*LOG(pej(js-1))
        pk3(i, js-1, k+1) = EXP(arg1)
      END DO
      pej(je+1) = ptop
      pej(je+2) = ptop
      DO k=1,npz
        pej(je+1) = pej(je+1) + delp(i, je+1, k)
        pej(je+2) = pej(je+2) + delp(i, je+2, k)
        arg1 = akap*LOG(pej(je+1))
        pk3(i, je+1, k+1) = EXP(arg1)
        arg1 = akap*LOG(pej(je+2))
        pk3(i, je+2, k+1) = EXP(arg1)
      END DO
    END DO
  END SUBROUTINE PK3_HALO
!  Differentiation of pln_halo in forward (tangent) mode:
!   variations   of useful results: pk3
!   with respect to varying inputs: pk3 delp
  SUBROUTINE PLN_HALO_TLM(is, ie, js, je, isd, ied, jsd, jed, npz, ptop&
&   , pk3, pk3_tl, delp, delp_tl)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp_tl
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3
    REAL, DIMENSION(isd:ied, jsd:jed, npz+1), INTENT(INOUT) :: pk3_tl
! Local:
    REAL :: pet
    REAL :: pet_tl
    INTEGER :: i, j, k
    INTRINSIC LOG
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,npz,ptop,delp,pk3) &
!$OMP                          private(pet)
    DO j=js,je
      DO i=is-2,is-1
        pet = ptop
        pet_tl = 0.0
        DO k=1,npz
          pet_tl = pet_tl + delp_tl(i, j, k)
          pet = pet + delp(i, j, k)
          pk3_tl(i, j, k+1) = pet_tl/pet
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
      DO i=ie+1,ie+2
        pet = ptop
        pet_tl = 0.0
        DO k=1,npz
          pet_tl = pet_tl + delp_tl(i, j, k)
          pet = pet + delp(i, j, k)
          pk3_tl(i, j, k+1) = pet_tl/pet
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ptop,delp,pk3) &
!$OMP                          private(pet)
    DO i=is-2,ie+2
      DO j=js-2,js-1
        pet = ptop
        pet_tl = 0.0
        DO k=1,npz
          pet_tl = pet_tl + delp_tl(i, j, k)
          pet = pet + delp(i, j, k)
          pk3_tl(i, j, k+1) = pet_tl/pet
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
      DO j=je+1,je+2
        pet = ptop
        pet_tl = 0.0
        DO k=1,npz
          pet_tl = pet_tl + delp_tl(i, j, k)
          pet = pet + delp(i, j, k)
          pk3_tl(i, j, k+1) = pet_tl/pet
          pk3(i, j, k+1) = LOG(pet)
        END DO
      END DO
    END DO
  END SUBROUTINE PLN_HALO_TLM
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
!  Differentiation of pe_halo in forward (tangent) mode:
!   variations   of useful results: pe
!   with respect to varying inputs: delp pe
  SUBROUTINE PE_HALO_TLM(is, ie, js, je, isd, ied, jsd, jed, npz, ptop, &
&   pe, pe_tl, delp, delp_tl)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed, npz
    REAL, INTENT(IN) :: ptop
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp_tl
    REAL, DIMENSION(is-1:ie+1, npz+1, js-1:je+1), INTENT(INOUT) :: pe
    REAL, DIMENSION(is-1:ie+1, npz+1, js-1:je+1), INTENT(INOUT) :: pe_tl
! Local:
    INTEGER :: i, j, k
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pe,delp,ptop)
    DO j=js,je
      pe_tl(is-1, 1, j) = 0.0
      pe(is-1, 1, j) = ptop
      pe_tl(ie+1, 1, j) = 0.0
      pe(ie+1, 1, j) = ptop
      DO k=1,npz
        pe_tl(is-1, k+1, j) = pe_tl(is-1, k, j) + delp_tl(is-1, j, k)
        pe(is-1, k+1, j) = pe(is-1, k, j) + delp(is-1, j, k)
        pe_tl(ie+1, k+1, j) = pe_tl(ie+1, k, j) + delp_tl(ie+1, j, k)
        pe(ie+1, k+1, j) = pe(ie+1, k, j) + delp(ie+1, j, k)
      END DO
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pe,delp,ptop)
    DO i=is-1,ie+1
      pe_tl(i, 1, js-1) = 0.0
      pe(i, 1, js-1) = ptop
      pe_tl(i, 1, je+1) = 0.0
      pe(i, 1, je+1) = ptop
      DO k=1,npz
        pe_tl(i, k+1, js-1) = pe_tl(i, k, js-1) + delp_tl(i, js-1, k)
        pe(i, k+1, js-1) = pe(i, k, js-1) + delp(i, js-1, k)
        pe_tl(i, k+1, je+1) = pe_tl(i, k, je+1) + delp_tl(i, je+1, k)
        pe(i, k+1, je+1) = pe(i, k, je+1) + delp(i, je+1, k)
      END DO
    END DO
  END SUBROUTINE PE_HALO_TLM
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
!  Differentiation of adv_pe in forward (tangent) mode:
!   variations   of useful results: om
!   with respect to varying inputs: ua om va pem
  SUBROUTINE ADV_PE_TLM(ua, ua_tl, va, va_tl, pem, pem_tl, om, om_tl, &
&   gridstruct, bd, npx, npy, npz, ng)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, npz, ng
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! Contra-variant wind components:
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(IN) :: ua&
&   , va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(IN) :: &
&   ua_tl, va_tl
! Pressure at edges:
    REAL, INTENT(IN) :: pem(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(IN) :: pem_tl(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(INOUT) :: om(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: om_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: up, vp
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: up_tl, vp_tl
    REAL :: v3(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: v3_tl(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: pin(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pin_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pb(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pb_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: grad(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: grad_tl(3, bd%is:bd%ie, bd%js:bd%je)
    REAL :: pdx(3, bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: pdx_tl(3, bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: pdy(3, bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: pdy_tl(3, bd%is:bd%ie+1, bd%js:bd%je)
    INTEGER :: i, j, k, n
    INTEGER :: is, ie, js, je
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    v3_tl = 0.0
    grad_tl = 0.0
    up_tl = 0.0
    pdx_tl = 0.0
    pdy_tl = 0.0
    pb_tl = 0.0
    vp_tl = 0.0
    pin_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ua,va,gridstruct,pem,npx,npy,ng,om) &
!$OMP                          private(n, pdx, pdy, pin, pb, up, vp, grad, v3)
    DO k=1,npz
      IF (k .EQ. npz) THEN
        DO j=js,je
          DO i=is,ie
            up_tl(i, j) = ua_tl(i, j, npz)
            up(i, j) = ua(i, j, npz)
            vp_tl(i, j) = va_tl(i, j, npz)
            vp(i, j) = va(i, j, npz)
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            up_tl(i, j) = 0.5*(ua_tl(i, j, k)+ua_tl(i, j, k+1))
            up(i, j) = 0.5*(ua(i, j, k)+ua(i, j, k+1))
            vp_tl(i, j) = 0.5*(va_tl(i, j, k)+va_tl(i, j, k+1))
            vp(i, j) = 0.5*(va(i, j, k)+va(i, j, k+1))
          END DO
        END DO
      END IF
! Compute Vect wind:
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            v3_tl(n, i, j) = gridstruct%ec1(n, i, j)*up_tl(i, j) + &
&             gridstruct%ec2(n, i, j)*vp_tl(i, j)
            v3(n, i, j) = up(i, j)*gridstruct%ec1(n, i, j) + vp(i, j)*&
&             gridstruct%ec2(n, i, j)
          END DO
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          pin_tl(i, j) = pem_tl(i, k+1, j)
          pin(i, j) = pem(i, k+1, j)
        END DO
      END DO
! Compute pe at 4 cell corners:
      CALL A2B_ORD2_TLM(pin, pin_tl, pb, pb_tl, gridstruct, npx, npy, is&
&                 , ie, js, je, ng)
      DO j=js,je+1
        DO i=is,ie
          DO n=1,3
            pdx_tl(n, i, j) = gridstruct%dx(i, j)*gridstruct%en1(n, i, j&
&             )*(pb_tl(i, j)+pb_tl(i+1, j))
            pdx(n, i, j) = (pb(i, j)+pb(i+1, j))*gridstruct%dx(i, j)*&
&             gridstruct%en1(n, i, j)
          END DO
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          DO n=1,3
            pdy_tl(n, i, j) = gridstruct%dy(i, j)*gridstruct%en2(n, i, j&
&             )*(pb_tl(i, j)+pb_tl(i, j+1))
            pdy(n, i, j) = (pb(i, j)+pb(i, j+1))*gridstruct%dy(i, j)*&
&             gridstruct%en2(n, i, j)
          END DO
        END DO
      END DO
! Compute grad (pe) by Green's theorem
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            grad_tl(n, i, j) = pdx_tl(n, i, j+1) - pdx_tl(n, i, j) - &
&             pdy_tl(n, i, j) + pdy_tl(n, i+1, j)
            grad(n, i, j) = pdx(n, i, j+1) - pdx(n, i, j) - pdy(n, i, j)&
&             + pdy(n, i+1, j)
          END DO
        END DO
      END DO
! Compute inner product: V3 * grad (pe)
      DO j=js,je
        DO i=is,ie
          om_tl(i, j, k) = om_tl(i, j, k) + 0.5*gridstruct%rarea(i, j)*(&
&           v3_tl(1, i, j)*grad(1, i, j)+v3(1, i, j)*grad_tl(1, i, j)+&
&           v3_tl(2, i, j)*grad(2, i, j)+v3(2, i, j)*grad_tl(2, i, j)+&
&           v3_tl(3, i, j)*grad(3, i, j)+v3(3, i, j)*grad_tl(3, i, j))
          om(i, j, k) = om(i, j, k) + 0.5*gridstruct%rarea(i, j)*(v3(1, &
&           i, j)*grad(1, i, j)+v3(2, i, j)*grad(2, i, j)+v3(3, i, j)*&
&           grad(3, i, j))
        END DO
      END DO
    END DO
  END SUBROUTINE ADV_PE_TLM
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
!  Differentiation of p_grad_c in forward (tangent) mode:
!   variations   of useful results: uc vc
!   with respect to varying inputs: gz uc pkc delpc vc
  SUBROUTINE P_GRAD_C_TLM(dt2, npz, delpc, delpc_tl, pkc, pkc_tl, gz, &
&   gz_tl, uc, uc_tl, vc, vc_tl, bd, rdxc, rdyc, hydrostatic)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npz
    REAL, INTENT(IN) :: dt2
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:, bd%jsd:, :), INTENT(IN) :: delpc
    REAL, DIMENSION(bd%isd:, bd%jsd:, :), INTENT(IN) :: delpc_tl
! pkc is pe**cappa     if hydrostatic
! pkc is full pressure if non-hydrostatic
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(IN) :: &
&   pkc, gz
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(IN) :: &
&   pkc_tl, gz_tl
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: uc_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: vc_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(IN) :: rdxc(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL, INTENT(IN) :: rdyc(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: hydrostatic
! Local:
    REAL :: wk(bd%is-1:bd%ie+1, bd%js-1:bd%je+1)
    REAL :: wk_tl(bd%is-1:bd%ie+1, bd%js-1:bd%je+1)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    wk_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,npz,hydrostatic,pkc,delpc,uc,dt2,rdxc,gz,vc,rdyc) &
!$OMP                          private(wk)
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js-1,je+1
          DO i=is-1,ie+1
            wk_tl(i, j) = pkc_tl(i, j, k+1) - pkc_tl(i, j, k)
            wk(i, j) = pkc(i, j, k+1) - pkc(i, j, k)
          END DO
        END DO
      ELSE
        DO j=js-1,je+1
          DO i=is-1,ie+1
            wk_tl(i, j) = delpc_tl(i, j, k)
            wk(i, j) = delpc(i, j, k)
          END DO
        END DO
      END IF
      DO j=js,je
        DO i=is,ie+1
          uc_tl(i, j, k) = uc_tl(i, j, k) + dt2*rdxc(i, j)*((gz_tl(i-1, &
&           j, k+1)-gz_tl(i, j, k))*(pkc(i, j, k+1)-pkc(i-1, j, k))+(gz(&
&           i-1, j, k+1)-gz(i, j, k))*(pkc_tl(i, j, k+1)-pkc_tl(i-1, j, &
&           k))+(gz_tl(i-1, j, k)-gz_tl(i, j, k+1))*(pkc(i-1, j, k+1)-&
&           pkc(i, j, k))+(gz(i-1, j, k)-gz(i, j, k+1))*(pkc_tl(i-1, j, &
&           k+1)-pkc_tl(i, j, k)))/(wk(i-1, j)+wk(i, j)) - dt2*rdxc(i, j&
&           )*(wk_tl(i-1, j)+wk_tl(i, j))*((gz(i-1, j, k+1)-gz(i, j, k))&
&           *(pkc(i, j, k+1)-pkc(i-1, j, k))+(gz(i-1, j, k)-gz(i, j, k+1&
&           ))*(pkc(i-1, j, k+1)-pkc(i, j, k)))/(wk(i-1, j)+wk(i, j))**2
          uc(i, j, k) = uc(i, j, k) + dt2*rdxc(i, j)/(wk(i-1, j)+wk(i, j&
&           ))*((gz(i-1, j, k+1)-gz(i, j, k))*(pkc(i, j, k+1)-pkc(i-1, j&
&           , k))+(gz(i-1, j, k)-gz(i, j, k+1))*(pkc(i-1, j, k+1)-pkc(i&
&           , j, k)))
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          vc_tl(i, j, k) = vc_tl(i, j, k) + dt2*rdyc(i, j)*((gz_tl(i, j-&
&           1, k+1)-gz_tl(i, j, k))*(pkc(i, j, k+1)-pkc(i, j-1, k))+(gz(&
&           i, j-1, k+1)-gz(i, j, k))*(pkc_tl(i, j, k+1)-pkc_tl(i, j-1, &
&           k))+(gz_tl(i, j-1, k)-gz_tl(i, j, k+1))*(pkc(i, j-1, k+1)-&
&           pkc(i, j, k))+(gz(i, j-1, k)-gz(i, j, k+1))*(pkc_tl(i, j-1, &
&           k+1)-pkc_tl(i, j, k)))/(wk(i, j-1)+wk(i, j)) - dt2*rdyc(i, j&
&           )*(wk_tl(i, j-1)+wk_tl(i, j))*((gz(i, j-1, k+1)-gz(i, j, k))&
&           *(pkc(i, j, k+1)-pkc(i, j-1, k))+(gz(i, j-1, k)-gz(i, j, k+1&
&           ))*(pkc(i, j-1, k+1)-pkc(i, j, k)))/(wk(i, j-1)+wk(i, j))**2
          vc(i, j, k) = vc(i, j, k) + dt2*rdyc(i, j)/(wk(i, j-1)+wk(i, j&
&           ))*((gz(i, j-1, k+1)-gz(i, j, k))*(pkc(i, j, k+1)-pkc(i, j-1&
&           , k))+(gz(i, j-1, k)-gz(i, j, k+1))*(pkc(i, j-1, k+1)-pkc(i&
&           , j, k)))
        END DO
      END DO
    END DO
  END SUBROUTINE P_GRAD_C_TLM
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
!  Differentiation of nh_p_grad in forward (tangent) mode:
!   variations   of useful results: gz u v delp pk pp
!   with respect to varying inputs: gz u v delp pk pp
  SUBROUTINE NH_P_GRAD_TLM(u, u_tl, v, v_tl, pp, pp_tl, gz, gz_tl, delp&
&   , delp_tl, pk, pk_tl, dt, ng, gridstruct, bd, npx, npy, npz, &
&   use_logp)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    LOGICAL, INTENT(IN) :: use_logp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! perturbation pressure
    REAL, INTENT(INOUT) :: pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pp_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! p**kappa
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! g * h
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL :: wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk1_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk_tl(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: du1, dv1, top_value
    REAL :: du1_tl, dv1_tl
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
        pp_tl(i, j, 1) = 0.0
        pp(i, j, 1) = 0.
        pk_tl(i, j, 1) = 0.0
        pk(i, j, 1) = top_value
      END DO
    END DO
    wk1_tl = 0.0
!$OMP parallel do default(none) shared(isd,jsd,npz,pp,gridstruct,npx,npy,is,ie,js,je,ng,pk,gz) &
!$OMP                          private(wk1)
    DO k=1,npz+1
      IF (k .NE. 1) THEN
        CALL A2B_ORD4_TLM(pp(isd:ied, jsd:jed, k), pp_tl(isd:ied, jsd&
&                      :jed, k), wk1, wk1_tl, gridstruct, npx, npy, is, &
&                      ie, js, je, ng, .true.)
        CALL A2B_ORD4_TLM(pk(isd:ied, jsd:jed, k), pk_tl(isd:ied, jsd&
&                      :jed, k), wk1, wk1_tl, gridstruct, npx, npy, is, &
&                      ie, js, je, ng, .true.)
      END IF
      CALL A2B_ORD4_TLM(gz(isd:ied, jsd:jed, k), gz_tl(isd:ied, jsd:&
&                    jed, k), wk1, wk1_tl, gridstruct, npx, npy, is, ie&
&                    , js, je, ng, .true.)
    END DO
    wk_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,npz,delp,gridstruct,npx,npy,ng,isd,jsd, &
!$OMP                                  pk,dt,gz,u,pp,v) &
!$OMP                          private(wk1, wk, du1, dv1)
    DO k=1,npz
      CALL A2B_ORD4_TLM(delp(isd:ied, jsd:jed, k), delp_tl(isd:ied, &
&                    jsd:jed, k), wk1, wk1_tl, gridstruct, npx, npy, is&
&                    , ie, js, je, ng)
      DO j=js,je+1
        DO i=is,ie+1
          wk_tl(i, j) = pk_tl(i, j, k+1) - pk_tl(i, j, k)
          wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
! hydrostatic contributions from past time-step already added in the "beta" part
! Current gradient from "hydrostatic" components:
          du1_tl = dt*((gz_tl(i, j, k+1)-gz_tl(i+1, j, k))*(pk(i+1, j, k&
&           +1)-pk(i, j, k))+(gz(i, j, k+1)-gz(i+1, j, k))*(pk_tl(i+1, j&
&           , k+1)-pk_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i+1, j, k+1))*(&
&           pk(i, j, k+1)-pk(i+1, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*(&
&           pk_tl(i, j, k+1)-pk_tl(i+1, j, k)))/(wk(i, j)+wk(i+1, j)) - &
&           dt*(wk_tl(i, j)+wk_tl(i+1, j))*((gz(i, j, k+1)-gz(i+1, j, k)&
&           )*(pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, j, k+1)&
&           )*(pk(i, j, k+1)-pk(i+1, j, k)))/(wk(i, j)+wk(i+1, j))**2
          du1 = dt/(wk(i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*&
&           (pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*&
&           (pk(i, j, k+1)-pk(i+1, j, k)))
! Non-hydrostatic contribution
          u_tl(i, j, k) = gridstruct%rdx(i, j)*(u_tl(i, j, k)+du1_tl+dt*&
&           ((gz_tl(i, j, k+1)-gz_tl(i+1, j, k))*(pp(i+1, j, k+1)-pp(i, &
&           j, k))+(gz(i, j, k+1)-gz(i+1, j, k))*(pp_tl(i+1, j, k+1)-&
&           pp_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i+1, j, k+1))*(pp(i, j&
&           , k+1)-pp(i+1, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*(pp_tl(i&
&           , j, k+1)-pp_tl(i+1, j, k)))/(wk1(i, j)+wk1(i+1, j))-dt*(&
&           wk1_tl(i, j)+wk1_tl(i+1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*&
&           (pp(i+1, j, k+1)-pp(i, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*&
&           (pp(i, j, k+1)-pp(i+1, j, k)))/(wk1(i, j)+wk1(i+1, j))**2)
          u(i, j, k) = (u(i, j, k)+du1+dt/(wk1(i, j)+wk1(i+1, j))*((gz(i&
&           , j, k+1)-gz(i+1, j, k))*(pp(i+1, j, k+1)-pp(i, j, k))+(gz(i&
&           , j, k)-gz(i+1, j, k+1))*(pp(i, j, k+1)-pp(i+1, j, k))))*&
&           gridstruct%rdx(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
! Current gradient from "hydrostatic" components:
          dv1_tl = dt*((gz_tl(i, j, k+1)-gz_tl(i, j+1, k))*(pk(i, j+1, k&
&           +1)-pk(i, j, k))+(gz(i, j, k+1)-gz(i, j+1, k))*(pk_tl(i, j+1&
&           , k+1)-pk_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i, j+1, k+1))*(&
&           pk(i, j, k+1)-pk(i, j+1, k))+(gz(i, j, k)-gz(i, j+1, k+1))*(&
&           pk_tl(i, j, k+1)-pk_tl(i, j+1, k)))/(wk(i, j)+wk(i, j+1)) - &
&           dt*(wk_tl(i, j)+wk_tl(i, j+1))*((gz(i, j, k+1)-gz(i, j+1, k)&
&           )*(pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1, k+1)&
&           )*(pk(i, j, k+1)-pk(i, j+1, k)))/(wk(i, j)+wk(i, j+1))**2
          dv1 = dt/(wk(i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*&
&           (pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*&
&           (pk(i, j, k+1)-pk(i, j+1, k)))
! Non-hydrostatic contribution
          v_tl(i, j, k) = gridstruct%rdy(i, j)*(v_tl(i, j, k)+dv1_tl+dt*&
&           ((gz_tl(i, j, k+1)-gz_tl(i, j+1, k))*(pp(i, j+1, k+1)-pp(i, &
&           j, k))+(gz(i, j, k+1)-gz(i, j+1, k))*(pp_tl(i, j+1, k+1)-&
&           pp_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i, j+1, k+1))*(pp(i, j&
&           , k+1)-pp(i, j+1, k))+(gz(i, j, k)-gz(i, j+1, k+1))*(pp_tl(i&
&           , j, k+1)-pp_tl(i, j+1, k)))/(wk1(i, j)+wk1(i, j+1))-dt*(&
&           wk1_tl(i, j)+wk1_tl(i, j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*&
&           (pp(i, j+1, k+1)-pp(i, j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*&
&           (pp(i, j, k+1)-pp(i, j+1, k)))/(wk1(i, j)+wk1(i, j+1))**2)
          v(i, j, k) = (v(i, j, k)+dv1+dt/(wk1(i, j)+wk1(i, j+1))*((gz(i&
&           , j, k+1)-gz(i, j+1, k))*(pp(i, j+1, k+1)-pp(i, j, k))+(gz(i&
&           , j, k)-gz(i, j+1, k+1))*(pp(i, j, k+1)-pp(i, j+1, k))))*&
&           gridstruct%rdy(i, j)
        END DO
      END DO
    END DO
  END SUBROUTINE NH_P_GRAD_TLM
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
!  Differentiation of split_p_grad in forward (tangent) mode:
!   variations   of useful results: gz du u dv v delp pk pp
!   with respect to varying inputs: gz du u dv v delp pk pp
  SUBROUTINE SPLIT_P_GRAD_TLM(u, u_tl, v, v_tl, pp, pp_tl, gz, gz_tl, du&
&   , du_tl, dv, dv_tl, delp, delp_tl, pk, pk_tl, beta, dt, ng, &
&   gridstruct, bd, npx, npy, npz, use_logp)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: beta, dt
    LOGICAL, INTENT(IN) :: use_logp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! perturbation pressure
    REAL, INTENT(INOUT) :: pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pp_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! p**kappa
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
! g * h
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: du_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dv_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL :: wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk1_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk_tl(bd%is:bd%ie+1, bd%js:bd%je+1)
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
        pp_tl(i, j, 1) = 0.0
        pp(i, j, 1) = 0.
        pk_tl(i, j, 1) = 0.0
        pk(i, j, 1) = top_value
      END DO
    END DO
    wk1_tl = 0.0
!$OMP parallel do default(none) shared(isd,jsd,npz,pp,gridstruct,npx,npy,is,ie,js,je,ng,pk,gz) &
!$OMP                          private(wk1)
    DO k=1,npz+1
      IF (k .NE. 1) THEN
        CALL A2B_ORD4_TLM(pp(isd:ied, jsd:jed, k), pp_tl(isd:ied, jsd&
&                      :jed, k), wk1, wk1_tl, gridstruct, npx, npy, is, &
&                      ie, js, je, ng, .true.)
        CALL A2B_ORD4_TLM(pk(isd:ied, jsd:jed, k), pk_tl(isd:ied, jsd&
&                      :jed, k), wk1, wk1_tl, gridstruct, npx, npy, is, &
&                      ie, js, je, ng, .true.)
      END IF
      CALL A2B_ORD4_TLM(gz(isd:ied, jsd:jed, k), gz_tl(isd:ied, jsd:&
&                    jed, k), wk1, wk1_tl, gridstruct, npx, npy, is, ie&
&                    , js, je, ng, .true.)
    END DO
    wk_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,isd,jsd,npz,delp,gridstruct,npx,npy,ng, &
!$OMP                                  pk,u,beta,du,dt,gz,alpha,pp,v,dv) &
!$OMP                          private(wk1, wk)
    DO k=1,npz
      CALL A2B_ORD4_TLM(delp(isd:ied, jsd:jed, k), delp_tl(isd:ied, &
&                    jsd:jed, k), wk1, wk1_tl, gridstruct, npx, npy, is&
&                    , ie, js, je, ng)
      DO j=js,je+1
        DO i=is,ie+1
          wk_tl(i, j) = pk_tl(i, j, k+1) - pk_tl(i, j, k)
          wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          u_tl(i, j, k) = u_tl(i, j, k) + beta*du_tl(i, j, k)
          u(i, j, k) = u(i, j, k) + beta*du(i, j, k)
! hydrostatic contributions from past time-step already added in the "beta" part
! Current gradient from "hydrostatic" components:
!---------------------------------------------------------------------------------
          du_tl(i, j, k) = dt*((gz_tl(i, j, k+1)-gz_tl(i+1, j, k))*(pk(i&
&           +1, j, k+1)-pk(i, j, k))+(gz(i, j, k+1)-gz(i+1, j, k))*(&
&           pk_tl(i+1, j, k+1)-pk_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i+1&
&           , j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k))+(gz(i, j, k)-gz(i+1&
&           , j, k+1))*(pk_tl(i, j, k+1)-pk_tl(i+1, j, k)))/(wk(i, j)+wk&
&           (i+1, j)) - dt*(wk_tl(i, j)+wk_tl(i+1, j))*((gz(i, j, k+1)-&
&           gz(i+1, j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz&
&           (i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))/(wk(i, j)+wk(i&
&           +1, j))**2
          du(i, j, k) = dt/(wk(i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1&
&           , j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, &
&           j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
!---------------------------------------------------------------------------------
! Non-hydrostatic contribution
          u_tl(i, j, k) = gridstruct%rdx(i, j)*(u_tl(i, j, k)+alpha*&
&           du_tl(i, j, k)+dt*((gz_tl(i, j, k+1)-gz_tl(i+1, j, k))*(pp(i&
&           +1, j, k+1)-pp(i, j, k))+(gz(i, j, k+1)-gz(i+1, j, k))*(&
&           pp_tl(i+1, j, k+1)-pp_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i+1&
&           , j, k+1))*(pp(i, j, k+1)-pp(i+1, j, k))+(gz(i, j, k)-gz(i+1&
&           , j, k+1))*(pp_tl(i, j, k+1)-pp_tl(i+1, j, k)))/(wk1(i, j)+&
&           wk1(i+1, j))-dt*(wk1_tl(i, j)+wk1_tl(i+1, j))*((gz(i, j, k+1&
&           )-gz(i+1, j, k))*(pp(i+1, j, k+1)-pp(i, j, k))+(gz(i, j, k)-&
&           gz(i+1, j, k+1))*(pp(i, j, k+1)-pp(i+1, j, k)))/(wk1(i, j)+&
&           wk1(i+1, j))**2)
          u(i, j, k) = (u(i, j, k)+alpha*du(i, j, k)+dt/(wk1(i, j)+wk1(i&
&           +1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*(pp(i+1, j, k+1)-pp(i&
&           , j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*(pp(i, j, k+1)-pp(i+1&
&           , j, k))))*gridstruct%rdx(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v_tl(i, j, k) = v_tl(i, j, k) + beta*dv_tl(i, j, k)
          v(i, j, k) = v(i, j, k) + beta*dv(i, j, k)
! Current gradient from "hydrostatic" components:
!---------------------------------------------------------------------------------
          dv_tl(i, j, k) = dt*((gz_tl(i, j, k+1)-gz_tl(i, j+1, k))*(pk(i&
&           , j+1, k+1)-pk(i, j, k))+(gz(i, j, k+1)-gz(i, j+1, k))*(&
&           pk_tl(i, j+1, k+1)-pk_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i, &
&           j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k))+(gz(i, j, k)-gz(i, &
&           j+1, k+1))*(pk_tl(i, j, k+1)-pk_tl(i, j+1, k)))/(wk(i, j)+wk&
&           (i, j+1)) - dt*(wk_tl(i, j)+wk_tl(i, j+1))*((gz(i, j, k+1)-&
&           gz(i, j+1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz&
&           (i, j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))/(wk(i, j)+wk(i&
&           , j+1))**2
          dv(i, j, k) = dt/(wk(i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j&
&           +1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1&
&           , k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
!---------------------------------------------------------------------------------
! Non-hydrostatic contribution
          v_tl(i, j, k) = gridstruct%rdy(i, j)*(v_tl(i, j, k)+alpha*&
&           dv_tl(i, j, k)+dt*((gz_tl(i, j, k+1)-gz_tl(i, j+1, k))*(pp(i&
&           , j+1, k+1)-pp(i, j, k))+(gz(i, j, k+1)-gz(i, j+1, k))*(&
&           pp_tl(i, j+1, k+1)-pp_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i, &
&           j+1, k+1))*(pp(i, j, k+1)-pp(i, j+1, k))+(gz(i, j, k)-gz(i, &
&           j+1, k+1))*(pp_tl(i, j, k+1)-pp_tl(i, j+1, k)))/(wk1(i, j)+&
&           wk1(i, j+1))-dt*(wk1_tl(i, j)+wk1_tl(i, j+1))*((gz(i, j, k+1&
&           )-gz(i, j+1, k))*(pp(i, j+1, k+1)-pp(i, j, k))+(gz(i, j, k)-&
&           gz(i, j+1, k+1))*(pp(i, j, k+1)-pp(i, j+1, k)))/(wk1(i, j)+&
&           wk1(i, j+1))**2)
          v(i, j, k) = (v(i, j, k)+alpha*dv(i, j, k)+dt/(wk1(i, j)+wk1(i&
&           , j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*(pp(i, j+1, k+1)-pp(i&
&           , j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*(pp(i, j, k+1)-pp(i, &
&           j+1, k))))*gridstruct%rdy(i, j)
        END DO
      END DO
    END DO
  END SUBROUTINE SPLIT_P_GRAD_TLM
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
!  Differentiation of one_grad_p in forward (tangent) mode:
!   variations   of useful results: gz u v delp pk
!   with respect to varying inputs: gz u v delp pk divg2
  SUBROUTINE ONE_GRAD_P_TLM(u, u_tl, v, v_tl, pk, pk_tl, gz, gz_tl, &
&   divg2, divg2_tl, delp, delp_tl, dt, ng, gridstruct, bd, npx, npy, &
&   npz, ptop, hydrostatic, a2b_ord, d_ext)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz, a2b_ord
    REAL, INTENT(IN) :: dt, ptop, d_ext
    LOGICAL, INTENT(IN) :: hydrostatic
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(IN) :: divg2_tl(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: wk
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: wk_tl
    REAL :: wk1(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk1_tl(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL :: wk2(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: wk2_tl(bd%is:bd%ie, bd%js:bd%je+1)
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
        pk_tl(i, j, 1) = 0.0
        pk(i, j, 1) = top_value
      END DO
    END DO
    wk_tl = 0.0
!$OMP parallel do default(none) shared(npz,isd,jsd,pk,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(pk(isd:ied, jsd:jed, k), pk_tl(isd:ied, jsd&
&                      :jed, k), wk, wk_tl, gridstruct, npx, npy, is, ie&
&                      , js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2_TLM(pk(isd:ied, jsd:jed, k), pk_tl(isd:ied, jsd:&
&                   jed, k), wk, wk_tl, gridstruct, npx, npy, is, ie, js&
&                   , je, ng, .true.)
      END IF
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,gz,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(gz(isd:ied, jsd:jed, k), gz_tl(isd:ied, jsd&
&                      :jed, k), wk, wk_tl, gridstruct, npx, npy, is, ie&
&                      , js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2_TLM(gz(isd:ied, jsd:jed, k), gz_tl(isd:ied, jsd:&
&                   jed, k), wk, wk_tl, gridstruct, npx, npy, is, ie, js&
&                   , je, ng, .true.)
      END IF
    END DO
    IF (d_ext .GT. 0.) THEN
      wk2_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,wk2,divg2)
      DO j=js,je+1
        DO i=is,ie
          wk2_tl(i, j) = divg2_tl(i, j) - divg2_tl(i+1, j)
          wk2(i, j) = divg2(i, j) - divg2(i+1, j)
        END DO
      END DO
      wk1_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,wk1,divg2)
      DO j=js,je
        DO i=is,ie+1
          wk1_tl(i, j) = divg2_tl(i, j) - divg2_tl(i, j+1)
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
      wk1_tl = 0.0
      wk2_tl = 0.0
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pk,delp,hydrostatic,a2b_ord,gridstruct, &
!$OMP                                  npx,npy,isd,jsd,ng,u,v,wk2,dt,gz,wk1) &
!$OMP                          private(wk)
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js,je+1
          DO i=is,ie+1
            wk_tl(i, j) = pk_tl(i, j, k+1) - pk_tl(i, j, k)
            wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
          END DO
        END DO
      ELSE IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(delp(isd:ied, jsd:jed, k), delp_tl(isd:ied&
&                      , jsd:jed, k), wk, wk_tl, gridstruct, npx, npy, &
&                      is, ie, js, je, ng)
      ELSE
        CALL A2B_ORD2_TLM(delp(isd:ied, jsd:jed, k), delp_tl(isd:ied, &
&                   jsd:jed, k), wk, wk_tl, gridstruct, npx, npy, is, ie&
&                   , js, je, ng)
      END IF
      DO j=js,je+1
        DO i=is,ie
          u_tl(i, j, k) = gridstruct%rdx(i, j)*(wk2_tl(i, j)+u_tl(i, j, &
&           k)+dt*((gz_tl(i, j, k+1)-gz_tl(i+1, j, k))*(pk(i+1, j, k+1)-&
&           pk(i, j, k))+(gz(i, j, k+1)-gz(i+1, j, k))*(pk_tl(i+1, j, k+&
&           1)-pk_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i+1, j, k+1))*(pk(i&
&           , j, k+1)-pk(i+1, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*(&
&           pk_tl(i, j, k+1)-pk_tl(i+1, j, k)))/(wk(i, j)+wk(i+1, j))-dt&
&           *(wk_tl(i, j)+wk_tl(i+1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*&
&           (pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*&
&           (pk(i, j, k+1)-pk(i+1, j, k)))/(wk(i, j)+wk(i+1, j))**2)
          u(i, j, k) = gridstruct%rdx(i, j)*(wk2(i, j)+u(i, j, k)+dt/(wk&
&           (i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1, j, k))*(pk(i+1, j&
&           , k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, j, k+1))*(pk(i, j, &
&           k+1)-pk(i+1, j, k))))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v_tl(i, j, k) = gridstruct%rdy(i, j)*(wk1_tl(i, j)+v_tl(i, j, &
&           k)+dt*((gz_tl(i, j, k+1)-gz_tl(i, j+1, k))*(pk(i, j+1, k+1)-&
&           pk(i, j, k))+(gz(i, j, k+1)-gz(i, j+1, k))*(pk_tl(i, j+1, k+&
&           1)-pk_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i, j+1, k+1))*(pk(i&
&           , j, k+1)-pk(i, j+1, k))+(gz(i, j, k)-gz(i, j+1, k+1))*(&
&           pk_tl(i, j, k+1)-pk_tl(i, j+1, k)))/(wk(i, j)+wk(i, j+1))-dt&
&           *(wk_tl(i, j)+wk_tl(i, j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*&
&           (pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*&
&           (pk(i, j, k+1)-pk(i, j+1, k)))/(wk(i, j)+wk(i, j+1))**2)
          v(i, j, k) = gridstruct%rdy(i, j)*(wk1(i, j)+v(i, j, k)+dt/(wk&
&           (i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j+1, k))*(pk(i, j+1&
&           , k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1, k+1))*(pk(i, j, &
&           k+1)-pk(i, j+1, k))))
        END DO
      END DO
    END DO
  END SUBROUTINE ONE_GRAD_P_TLM
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
!  Differentiation of grad1_p_update in forward (tangent) mode:
!   variations   of useful results: gz du u dv v pk
!   with respect to varying inputs: gz du u dv v pk divg2
  SUBROUTINE GRAD1_P_UPDATE_TLM(divg2, divg2_tl, u, u_tl, v, v_tl, pk, &
&   pk_tl, gz, gz_tl, du, du_tl, dv, dv_tl, dt, ng, gridstruct, bd, npx&
&   , npy, npz, ptop, beta, a2b_ord)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz, a2b_ord
    REAL, INTENT(IN) :: dt, ptop, beta
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(IN) :: divg2_tl(bd%is:bd%ie+1, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: pk_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: gz_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: u_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: v_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: du(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: du_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: dv(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dv_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    TYPE(FV_GRID_TYPE), INTENT(INOUT), TARGET :: gridstruct
! Local:
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
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
        pk_tl(i, j, 1) = 0.0
        pk(i, j, 1) = top_value
      END DO
    END DO
    wk_tl = 0.0
!$OMP parallel do default(none) shared(npz,isd,jsd,pk,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(pk(isd:ied, jsd:jed, k), pk_tl(isd:ied, jsd&
&                      :jed, k), wk, wk_tl, gridstruct, npx, npy, is, ie&
&                      , js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2_TLM(pk(isd:ied, jsd:jed, k), pk_tl(isd:ied, jsd:&
&                   jed, k), wk, wk_tl, gridstruct, npx, npy, is, ie, js&
&                   , je, ng, .true.)
      END IF
    END DO
!$OMP parallel do default(none) shared(npz,isd,jsd,gz,gridstruct,npx,npy,is,ie,js,je,ng,a2b_ord) &
!$OMP                          private(wk)
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(gz(isd:ied, jsd:jed, k), gz_tl(isd:ied, jsd&
&                      :jed, k), wk, wk_tl, gridstruct, npx, npy, is, ie&
&                      , js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2_TLM(gz(isd:ied, jsd:jed, k), gz_tl(isd:ied, jsd:&
&                   jed, k), wk, wk_tl, gridstruct, npx, npy, is, ie, js&
&                   , je, ng, .true.)
      END IF
    END DO
!$OMP parallel do default(none) shared(npz,is,ie,js,je,pk,u,beta,gz,divg2,alpha, &
!$OMP                                  gridstruct,v,dt,du,dv) &
!$OMP                          private(wk)
    DO k=1,npz
      DO j=js,je+1
        DO i=is,ie+1
          wk_tl(i, j) = pk_tl(i, j, k+1) - pk_tl(i, j, k)
          wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          u_tl(i, j, k) = u_tl(i, j, k) + beta*du_tl(i, j, k)
          u(i, j, k) = u(i, j, k) + beta*du(i, j, k)
          du_tl(i, j, k) = dt*((gz_tl(i, j, k+1)-gz_tl(i+1, j, k))*(pk(i&
&           +1, j, k+1)-pk(i, j, k))+(gz(i, j, k+1)-gz(i+1, j, k))*(&
&           pk_tl(i+1, j, k+1)-pk_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i+1&
&           , j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k))+(gz(i, j, k)-gz(i+1&
&           , j, k+1))*(pk_tl(i, j, k+1)-pk_tl(i+1, j, k)))/(wk(i, j)+wk&
&           (i+1, j)) - dt*(wk_tl(i, j)+wk_tl(i+1, j))*((gz(i, j, k+1)-&
&           gz(i+1, j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz&
&           (i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))/(wk(i, j)+wk(i&
&           +1, j))**2
          du(i, j, k) = dt/(wk(i, j)+wk(i+1, j))*((gz(i, j, k+1)-gz(i+1&
&           , j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i+1, &
&           j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
          u_tl(i, j, k) = gridstruct%rdx(i, j)*(u_tl(i, j, k)+divg2_tl(i&
&           , j)-divg2_tl(i+1, j)+alpha*du_tl(i, j, k))
          u(i, j, k) = (u(i, j, k)+divg2(i, j)-divg2(i+1, j)+alpha*du(i&
&           , j, k))*gridstruct%rdx(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v_tl(i, j, k) = v_tl(i, j, k) + beta*dv_tl(i, j, k)
          v(i, j, k) = v(i, j, k) + beta*dv(i, j, k)
          dv_tl(i, j, k) = dt*((gz_tl(i, j, k+1)-gz_tl(i, j+1, k))*(pk(i&
&           , j+1, k+1)-pk(i, j, k))+(gz(i, j, k+1)-gz(i, j+1, k))*(&
&           pk_tl(i, j+1, k+1)-pk_tl(i, j, k))+(gz_tl(i, j, k)-gz_tl(i, &
&           j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k))+(gz(i, j, k)-gz(i, &
&           j+1, k+1))*(pk_tl(i, j, k+1)-pk_tl(i, j+1, k)))/(wk(i, j)+wk&
&           (i, j+1)) - dt*(wk_tl(i, j)+wk_tl(i, j+1))*((gz(i, j, k+1)-&
&           gz(i, j+1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz&
&           (i, j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))/(wk(i, j)+wk(i&
&           , j+1))**2
          dv(i, j, k) = dt/(wk(i, j)+wk(i, j+1))*((gz(i, j, k+1)-gz(i, j&
&           +1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gz(i, j, k)-gz(i, j+1&
&           , k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
          v_tl(i, j, k) = gridstruct%rdy(i, j)*(v_tl(i, j, k)+divg2_tl(i&
&           , j)-divg2_tl(i, j+1)+alpha*dv_tl(i, j, k))
          v(i, j, k) = (v(i, j, k)+divg2(i, j)-divg2(i, j+1)+alpha*dv(i&
&           , j, k))*gridstruct%rdy(i, j)
        END DO
      END DO
    END DO
  END SUBROUTINE GRAD1_P_UPDATE_TLM
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
!  Differentiation of mix_dp in forward (tangent) mode:
!   variations   of useful results: w delp pt
!   with respect to varying inputs: w delp pt
  SUBROUTINE MIX_DP_TLM(hydrostatic, w, w_tl, delp, delp_tl, pt, pt_tl, &
&   km, ak, bk, cg, fv_debug, bd)
    IMPLICIT NONE
!      if ( ip/=0 ) write(*,*) 'Warning: Mix_dp', mpp_pe(), j, ip
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: ak(km+1), bk(km+1)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   pt, delp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   pt_tl, delp_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(INOUT) :: &
&   w_tl
    LOGICAL, INTENT(IN) :: hydrostatic, cg, fv_debug
! Local:
    REAL :: dp, dpmin
    REAL :: dp_tl
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
            dp_tl = -delp_tl(i, j, k)
            dp = dpmin - delp(i, j, k)
            pt_tl(i, j, k) = (pt_tl(i, j, k)*delp(i, j, k)+pt(i, j, k)*&
&             delp_tl(i, j, k)+pt_tl(i, j, k+1)*dp+pt(i, j, k+1)*dp_tl)/&
&             dpmin
            pt(i, j, k) = (pt(i, j, k)*delp(i, j, k)+pt(i, j, k+1)*dp)/&
&             dpmin
            IF (.NOT.hydrostatic) THEN
              w_tl(i, j, k) = (w_tl(i, j, k)*delp(i, j, k)+w(i, j, k)*&
&               delp_tl(i, j, k)+w_tl(i, j, k+1)*dp+w(i, j, k+1)*dp_tl)/&
&               dpmin
              w(i, j, k) = (w(i, j, k)*delp(i, j, k)+w(i, j, k+1)*dp)/&
&               dpmin
            END IF
            delp_tl(i, j, k) = 0.0
            delp(i, j, k) = dpmin
            delp_tl(i, j, k+1) = delp_tl(i, j, k+1) - dp_tl
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
          dp_tl = -delp_tl(i, j, km)
          dp = dpmin - delp(i, j, km)
          pt_tl(i, j, km) = (pt_tl(i, j, km)*delp(i, j, km)+pt(i, j, km)&
&           *delp_tl(i, j, km)+pt_tl(i, j, km-1)*dp+pt(i, j, km-1)*dp_tl&
&           )/dpmin
          pt(i, j, km) = (pt(i, j, km)*delp(i, j, km)+pt(i, j, km-1)*dp)&
&           /dpmin
          IF (.NOT.hydrostatic) THEN
            w_tl(i, j, km) = (w_tl(i, j, km)*delp(i, j, km)+w(i, j, km)*&
&             delp_tl(i, j, km)+w_tl(i, j, km-1)*dp+w(i, j, km-1)*dp_tl)&
&             /dpmin
            w(i, j, km) = (w(i, j, km)*delp(i, j, km)+w(i, j, km-1)*dp)/&
&             dpmin
          END IF
          delp_tl(i, j, km) = 0.0
          delp(i, j, km) = dpmin
          delp_tl(i, j, km-1) = delp_tl(i, j, km-1) - dp_tl
          delp(i, j, km-1) = delp(i, j, km-1) - dp
          ip = ip + 1
        END IF
      END DO
      IF (fv_debug .AND. ip .NE. 0) WRITE(*, *) 'Warning: Mix_dp', &
&                                   MPP_PE(), j, ip
    END DO
  END SUBROUTINE MIX_DP_TLM
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
!  Differentiation of geopk in forward (tangent) mode:
!   variations   of useful results: peln gz pkz pe pk
!   with respect to varying inputs: peln gz delp pkz pe pk pt
  SUBROUTINE GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, delp, delp_tl, pk&
&   , pk_tl, gz, gz_tl, hs, pt, pt_tl, q_con, pkz, pkz_tl, km, akap, cg&
&   , nested, computehalo, npx, npy, a2b_ord, bd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km, npx, npy, a2b_ord
    REAL, INTENT(IN) :: akap, ptop
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: hs(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(IN) :: pt&
&   , delp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(IN) :: &
&   pt_tl, delp_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km), INTENT(IN) :: &
&   q_con
    LOGICAL, INTENT(IN) :: cg, nested, computehalo
! !OUTPUT PARAMETERS
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km+1), INTENT(OUT) :: &
&   gz, pk
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, km+1), INTENT(OUT) :: &
&   gz_tl, pk_tl
    REAL, INTENT(OUT) :: pe(bd%is-1:bd%ie+1, km+1, bd%js-1:bd%je+1)
    REAL, INTENT(OUT) :: pe_tl(bd%is-1:bd%ie+1, km+1, bd%js-1:bd%je+1)
! ln(pe)
    REAL, INTENT(OUT) :: peln(bd%is:bd%ie, km+1, bd%js:bd%je)
    REAL, INTENT(OUT) :: peln_tl(bd%is:bd%ie, km+1, bd%js:bd%je)
    REAL, INTENT(OUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, km)
    REAL, INTENT(OUT) :: pkz_tl(bd%is:bd%ie, bd%js:bd%je, km)
! !DESCRIPTION:
!    Calculates geopotential and pressure to the kappa.
! Local:
    REAL :: peg(bd%isd:bd%ied, km+1)
    REAL :: pkg(bd%isd:bd%ied, km+1)
    REAL(kind=8) :: p1d(bd%isd:bd%ied)
    REAL(kind=8) :: p1d_tl(bd%isd:bd%ied)
    REAL(kind=8) :: g1d(bd%isd:bd%ied)
    REAL(kind=8) :: g1d_tl(bd%isd:bd%ied)
    REAL :: logp(bd%isd:bd%ied)
    REAL :: logp_tl(bd%isd:bd%ied)
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
      IF (je .EQ. npy - 1) THEN
        jlast = jed
        g1d_tl = 0.0_8
        logp_tl = 0.0
        p1d_tl = 0.0_8
      ELSE
        g1d_tl = 0.0_8
        logp_tl = 0.0
        p1d_tl = 0.0_8
      END IF
    ELSE
      g1d_tl = 0.0_8
      logp_tl = 0.0
      p1d_tl = 0.0_8
    END IF
!$OMP parallel do default(none) shared(jfirst,jlast,ifirst,ilast,pk,km,gz,hs,ptop,ptk, &
!$OMP                                  js,je,is,ie,peln,peln1,pe,delp,akap,pt,CG,pkz,q_con) &
!$OMP                          private(peg, pkg, p1d, g1d, logp)
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        p1d_tl(i) = 0.0_8
        p1d(i) = ptop
        pk_tl(i, j, 1) = 0.0
        pk(i, j, 1) = ptk
        g1d_tl(i) = 0.0_8
        g1d(i) = hs(i, j)
        gz_tl(i, j, km+1) = 0.0
        gz(i, j, km+1) = hs(i, j)
      END DO
      IF (j .GE. js .AND. j .LE. je) THEN
        DO i=is,ie
          peln_tl(i, 1, j) = 0.0
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
          pe_tl(i, 1, j) = 0.0
          pe(i, 1, j) = ptop
        END DO
      END IF
! Top down
      DO k=2,km+1
        DO i=ifirst,ilast
          p1d_tl(i) = p1d_tl(i) + delp_tl(i, j, k-1)
          p1d(i) = p1d(i) + delp(i, j, k-1)
          logp_tl(i) = p1d_tl(i)/p1d(i)
          logp(i) = LOG(p1d(i))
          pk_tl(i, j, k) = akap*logp_tl(i)*EXP(akap*logp(i))
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
            pe_tl(i, k, j) = p1d_tl(i)
            pe(i, k, j) = p1d(i)
          END DO
          IF (j .GE. js .AND. j .LE. je) THEN
            DO i=is,ie
              peln_tl(i, k, j) = logp_tl(i)
              peln(i, k, j) = logp(i)
            END DO
          END IF
        END IF
      END DO
! Bottom up
      DO k=km,1,-1
        DO i=ifirst,ilast
          g1d_tl(i) = g1d_tl(i) + cp_air*(pt_tl(i, j, k)*(pk(i, j, k+1)-&
&           pk(i, j, k))+pt(i, j, k)*(pk_tl(i, j, k+1)-pk_tl(i, j, k)))
          g1d(i) = g1d(i) + cp_air*pt(i, j, k)*(pk(i, j, k+1)-pk(i, j, k&
&           ))
          gz_tl(i, j, k) = g1d_tl(i)
          gz(i, j, k) = g1d(i)
        END DO
      END DO
      IF (.NOT.cg .AND. j .GE. js .AND. j .LE. je) THEN
        DO k=1,km
          DO i=is,ie
            pkz_tl(i, j, k) = ((pk_tl(i, j, k+1)-pk_tl(i, j, k))*akap*(&
&             peln(i, k+1, j)-peln(i, k, j))-(pk(i, j, k+1)-pk(i, j, k))&
&             *akap*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(akap*(peln(i&
&             , k+1, j)-peln(i, k, j)))**2
            pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(akap*(peln(i, k+&
&             1, j)-peln(i, k, j)))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE GEOPK_TLM
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
!  Differentiation of del2_cubed in forward (tangent) mode:
!   variations   of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE DEL2_CUBED_TLM(q, q_tl, cd, gridstruct, domain, npx, npy, &
&   km, nmax, bd)
    IMPLICIT NONE
!---------------------------------------------------------------
! This routine is for filtering the omega field for the physics
!---------------------------------------------------------------
    INTEGER, INTENT(IN) :: npx, npy, km, nmax
! cd = K * da_min;   0 < K < 0.25
    REAL(kind=r_grid), INTENT(IN) :: cd
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL, PARAMETER :: r3=1./3.
    REAL :: fx(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy(bd%isd:bd%ied, bd%jsd&
&   :bd%jed+1)
    REAL :: fx_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy_tl(bd%isd:bd%ied, &
&   bd%jsd:bd%jed+1)
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
    CALL MPP_UPDATE_DOMAINS_TLM(q, q_tl, domain, complete=.true.)
    CALL TIMING_OFF('COMM_TOTAL')
    fx_tl = 0.0
    fy_tl = 0.0
    DO n=1,ntimes
      nt = ntimes - n
!$OMP parallel do default(none) shared(km,q,is,ie,js,je,npx,npy, &
!$OMP                                  nt,isd,jsd,gridstruct,bd, &
!$OMP                                  cd) &
!$OMP                          private(fx, fy)
      DO k=1,km
        IF (gridstruct%sw_corner) THEN
          q_tl(1, 1, k) = r3*(q_tl(1, 1, k)+q_tl(0, 1, k)+q_tl(1, 0, k))
          q(1, 1, k) = (q(1, 1, k)+q(0, 1, k)+q(1, 0, k))*r3
          q_tl(0, 1, k) = q_tl(1, 1, k)
          q(0, 1, k) = q(1, 1, k)
          q_tl(1, 0, k) = q_tl(1, 1, k)
          q(1, 0, k) = q(1, 1, k)
        END IF
        IF (gridstruct%se_corner) THEN
          q_tl(ie, 1, k) = r3*(q_tl(ie, 1, k)+q_tl(npx, 1, k)+q_tl(ie, 0&
&           , k))
          q(ie, 1, k) = (q(ie, 1, k)+q(npx, 1, k)+q(ie, 0, k))*r3
          q_tl(npx, 1, k) = q_tl(ie, 1, k)
          q(npx, 1, k) = q(ie, 1, k)
          q_tl(ie, 0, k) = q_tl(ie, 1, k)
          q(ie, 0, k) = q(ie, 1, k)
        END IF
        IF (gridstruct%ne_corner) THEN
          q_tl(ie, je, k) = r3*(q_tl(ie, je, k)+q_tl(npx, je, k)+q_tl(ie&
&           , npy, k))
          q(ie, je, k) = (q(ie, je, k)+q(npx, je, k)+q(ie, npy, k))*r3
          q_tl(npx, je, k) = q_tl(ie, je, k)
          q(npx, je, k) = q(ie, je, k)
          q_tl(ie, npy, k) = q_tl(ie, je, k)
          q(ie, npy, k) = q(ie, je, k)
        END IF
        IF (gridstruct%nw_corner) THEN
          q_tl(1, je, k) = r3*(q_tl(1, je, k)+q_tl(0, je, k)+q_tl(1, npy&
&           , k))
          q(1, je, k) = (q(1, je, k)+q(0, je, k)+q(1, npy, k))*r3
          q_tl(0, je, k) = q_tl(1, je, k)
          q(0, je, k) = q(1, je, k)
          q_tl(1, npy, k) = q_tl(1, je, k)
          q(1, npy, k) = q(1, je, k)
        END IF
        IF (nt .GT. 0) CALL COPY_CORNERS_TLM(q(isd:ied, jsd:jed, k), &
&                                      q_tl(isd:ied, jsd:jed, k), npx, &
&                                      npy, 1, gridstruct%nested, bd, &
&                                      gridstruct%sw_corner, gridstruct%&
&                                      se_corner, gridstruct%nw_corner, &
&                                      gridstruct%ne_corner)
        DO j=js-nt,je+nt
          DO i=is-nt,ie+1+nt
            fx_tl(i, j) = gridstruct%del6_v(i, j)*(q_tl(i-1, j, k)-q_tl(&
&             i, j, k))
            fx(i, j) = gridstruct%del6_v(i, j)*(q(i-1, j, k)-q(i, j, k))
          END DO
        END DO
        IF (nt .GT. 0) CALL COPY_CORNERS_TLM(q(isd:ied, jsd:jed, k), &
&                                      q_tl(isd:ied, jsd:jed, k), npx, &
&                                      npy, 2, gridstruct%nested, bd, &
&                                      gridstruct%sw_corner, gridstruct%&
&                                      se_corner, gridstruct%nw_corner, &
&                                      gridstruct%ne_corner)
        DO j=js-nt,je+1+nt
          DO i=is-nt,ie+nt
            fy_tl(i, j) = gridstruct%del6_u(i, j)*(q_tl(i, j-1, k)-q_tl(&
&             i, j, k))
            fy(i, j) = gridstruct%del6_u(i, j)*(q(i, j-1, k)-q(i, j, k))
          END DO
        END DO
        DO j=js-nt,je+nt
          DO i=is-nt,ie+nt
            q_tl(i, j, k) = q_tl(i, j, k) + cd*gridstruct%rarea(i, j)*(&
&             fx_tl(i, j)-fx_tl(i+1, j)+fy_tl(i, j)-fy_tl(i, j+1))
            q(i, j, k) = q(i, j, k) + cd*gridstruct%rarea(i, j)*(fx(i, j&
&             )-fx(i+1, j)+fy(i, j)-fy(i, j+1))
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE DEL2_CUBED_TLM
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
    REAL :: arg1
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
          arg1 = 0.5*pi*LOG(rf_cutoff/pfull(k))/LOG(rf_cutoff/ptop)
          rff(k) = dt/tau0*SIN(arg1)**2
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

end module dyn_core_tlm_mod
