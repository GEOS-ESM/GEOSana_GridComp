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
module fv_nesting_adm_mod

   use mpp_domains_mod,     only: mpp_update_domains, mpp_global_field
   use fv_mp_adm_mod,       only: mpp_update_domains_adm
   use field_manager_mod,   only: MODEL_ATMOS
   use tracer_manager_mod,  only: get_tracer_index
   use fv_sg_mod,           only: neg_adj3
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
   use mpp_domains_mod,     only: DGRID_NE, domain2D
   use fv_restart_adm_mod,  only: d2a_setup, d2c_setup
   use fv_restart_adm_mod,  only: d2c_setup_fwd, d2c_setup_bwd
   use mpp_mod,             only: mpp_sync_self, mpp_sync, mpp_send, mpp_recv, mpp_error, FATAL
   use mpp_domains_mod,     only: mpp_global_sum, BITWISE_EFP_SUM, BITWISE_EXACT_SUM
   use boundary_mod,        only: update_coarse_grid
   use boundary_mod,        only: nested_grid_BC_send, nested_grid_BC_recv, nested_grid_BC_save_proc
   use fv_mp_mod,           only: is, ie, js, je, isd, ied, jsd, jed, isc, iec, jsc, jec
   use fv_arrays_mod,       only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_nest_BC_type_3D
   use fv_arrays_mod,       only: allocate_fv_nest_BC_type, fv_atmos_type, fv_grid_bounds_type
   use fv_grid_utils_mod,   only: ptop_min, g_sum, cubed_to_latlon, f_p
   use init_hydro_mod,      only: p_var
   use constants_mod,       only: grav, pi=>pi_8, radius, hlv, rdgas, cp_air, rvgas, cp_vapor, kappa
   use fv_mapz_mod,         only: mappm
   use fv_timing_mod,       only: timing_on, timing_off
   use fv_mp_mod,           only: is_master
   use fv_mp_mod,           only: mp_reduce_sum
   use fv_diagnostics_mod,  only: sphum_ll_fix, range_check
   use sw_core_mod,         only: divergence_corner, divergence_corner_nest

   use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                            pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

implicit none
   logical :: RF_initialized = .false.
   logical :: bad_range
   real, allocatable ::  rf(:), rw(:)
   integer :: kmax=1
   !Arrays for global grid total energy, used for grid nesting
   real, allocatable :: te_2d_coarse(:,:)
   real, allocatable :: dp1_coarse(:,:,:)

   !For nested grid buffers
	!Individual structures are allocated by nested_grid_BC_recv
   type(fv_nest_BC_type_3d) :: u_buf, v_buf, uc_buf, vc_buf, delp_buf, delz_buf, pt_buf, pkz_buf, w_buf, divg_buf
   type(fv_nest_BC_type_3d), allocatable:: q_buf(:)
!#ifdef USE_COND
   real, dimension(:,:,:), allocatable, target :: dum_West, dum_East, dum_North, dum_South
!#endif

private
public :: twoway_nesting, setup_nested_grid_BCs
public :: setup_nested_grid_BCs_adm

!---- version number -----
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of setup_nested_grid_bcs in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a
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
!   gradient     of useful results: u v uc vc
!   with respect to varying inputs: u v
!!!! NOTE: Many of the routines here and in boundary.F90 have a lot of
!!!!   redundant code, which could be cleaned up and simplified.
  SUBROUTINE SETUP_NESTED_GRID_BCS_ADM(npx, npy, npz, zvir, ncnst, u, &
&   u_ad, v, v_ad, w, pt, delp, delz, q, uc, uc_ad, vc, vc_ad, pkz, &
&   nested, inline_q, make_nh, ng, gridstruct, flagstruct, neststruct, &
&   nest_timestep, tracer_nest_timestep, domain, bd, nwat)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: zvir
    INTEGER, INTENT(IN) :: npx, npy, npz
    INTEGER, INTENT(IN) :: ncnst, ng, nwat
    LOGICAL, INTENT(IN) :: inline_q, make_nh, nested
! D grid zonal wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u_ad
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v_ad
!  W (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:, bd%jsd:, :)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! height thickness (m)
    REAL, INTENT(INOUT) :: delz(bd%isd:, bd%jsd:, :)
! specific humidity and constituents
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst)
! (uc,vc) mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: uc_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
    REAL, INTENT(INOUT) :: vc_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! finite-volume mean pk
    REAL, INTENT(INOUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
    INTEGER, INTENT(INOUT) :: nest_timestep, tracer_nest_timestep
    TYPE(FV_GRID_TYPE), INTENT(INOUT) :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(INOUT) :: flagstruct
    TYPE(FV_NEST_TYPE), INTENT(INOUT), TARGET :: neststruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL :: divg(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL :: ua(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: va(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pkz_coarse(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    INTEGER :: i, j, k, n, p, sphum
    LOGICAL :: do_pd
    TYPE(FV_NEST_BC_TYPE_3D) :: pkz_bc
!local pointers
    LOGICAL, POINTER :: child_grids(:)
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ANY
    INTRINSIC ALLOCATED
    INTRINSIC SIZE
    LOGICAL :: arg1
    INTEGER :: branch
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
!IF nested, set up nested grid BCs for time-interpolation
!(actually applying the BCs is done in dyn_core
!compute uc/vc for nested-grid BCs
!!! CLEANUP: if we compute uc/vc here we don't need to do on the first call of c_sw, right?
    IF (ANY(neststruct%child_grids)) THEN
!$OMP parallel do default(none) shared(isd,jsd,ied,jed,is,ie,js,je,npx,npy,npz, &
!$OMP       gridstruct,flagstruct,bd,u,v,uc,vc,nested,divg) &
!$OMP       private(ua,va)
      DO k=1,npz
        CALL D2C_SETUP_FWD(u(isd, jsd, k), v(isd, jsd, k), ua, va, uc(&
&                    isd, jsd, k), vc(isd, jsd, k), arg1, isd, ied, jsd&
&                    , jed, is, ie, js, je, npx, npy, gridstruct%&
&                    grid_type, gridstruct%nested, gridstruct%se_corner&
&                    , gridstruct%sw_corner, gridstruct%ne_corner, &
&                    gridstruct%nw_corner, gridstruct%rsin_u, gridstruct&
&                    %rsin_v, gridstruct%cosa_s, gridstruct%rsin2)
      END DO
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
!! Nested grid: receive from parent grid
    IF (neststruct%nested) THEN
      IF (.NOT.ALLOCATED(q_buf)) THEN
        ALLOCATE(q_buf(ncnst))
        DEALLOCATE(q_buf)
      END IF
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO k=npz,1,-1
        CALL D2C_SETUP_BWD(u(isd, jsd, k), u_ad(isd, jsd, k), v(isd, jsd&
&                    , k), v_ad(isd, jsd, k), ua, va, uc(isd, jsd, k), &
&                    uc_ad(isd, jsd, k), vc(isd, jsd, k), vc_ad(isd, jsd&
&                    , k), arg1, isd, ied, jsd, jed, is, ie, js, je, npx&
&                    , npy, gridstruct%grid_type, gridstruct%nested, &
&                    gridstruct%se_corner, gridstruct%sw_corner, &
&                    gridstruct%ne_corner, gridstruct%nw_corner, &
&                    gridstruct%rsin_u, gridstruct%rsin_v, gridstruct%&
&                    cosa_s, gridstruct%rsin2)
      END DO
      CALL MPP_UPDATE_DOMAINS_ADM(u, u_ad, v, v_ad, domain, gridtype=&
&                           dgrid_ne, complete=.true.)
    END IF
  END SUBROUTINE SETUP_NESTED_GRID_BCS_ADM
!!!! NOTE: Many of the routines here and in boundary.F90 have a lot of
!!!!   redundant code, which could be cleaned up and simplified.
  SUBROUTINE SETUP_NESTED_GRID_BCS(npx, npy, npz, zvir, ncnst, u, v, w, &
&   pt, delp, delz, q, uc, vc, pkz, nested, inline_q, make_nh, ng, &
&   gridstruct, flagstruct, neststruct, nest_timestep, &
&   tracer_nest_timestep, domain, bd, nwat)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: zvir
    INTEGER, INTENT(IN) :: npx, npy, npz
    INTEGER, INTENT(IN) :: ncnst, ng, nwat
    LOGICAL, INTENT(IN) :: inline_q, make_nh, nested
! D grid zonal wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
!  W (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:, bd%jsd:, :)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! height thickness (m)
    REAL, INTENT(INOUT) :: delz(bd%isd:, bd%jsd:, :)
! specific humidity and constituents
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst)
! (uc,vc) mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! finite-volume mean pk
    REAL, INTENT(INOUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
    INTEGER, INTENT(INOUT) :: nest_timestep, tracer_nest_timestep
    TYPE(FV_GRID_TYPE), INTENT(INOUT) :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(INOUT) :: flagstruct
    TYPE(FV_NEST_TYPE), INTENT(INOUT), TARGET :: neststruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL :: divg(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, npz)
    REAL :: ua(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: va(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pkz_coarse(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    INTEGER :: i, j, k, n, p, sphum
    LOGICAL :: do_pd
    TYPE(FV_NEST_BC_TYPE_3D) :: pkz_bc
!local pointers
    LOGICAL, POINTER :: child_grids(:)
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ANY
    INTRINSIC ALLOCATED
    INTRINSIC SIZE
    LOGICAL :: arg1
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    child_grids => neststruct%child_grids
!IF nested, set up nested grid BCs for time-interpolation
!(actually applying the BCs is done in dyn_core
    nest_timestep = 0
    IF (.NOT.inline_q) tracer_nest_timestep = 0
    IF (neststruct%nested .AND. ((.NOT.neststruct%first_step) .OR. &
&       make_nh)) THEN
      do_pd = .true.
      CALL SET_BCS_T0(ncnst, flagstruct%hydrostatic, neststruct)
    ELSE
!On first timestep the t0 BCs are not initialized and may contain garbage
      do_pd = .false.
    END IF
!compute uc/vc for nested-grid BCs
!!! CLEANUP: if we compute uc/vc here we don't need to do on the first call of c_sw, right?
    IF (ANY(neststruct%child_grids)) THEN
      CALL TIMING_ON('COMM_TOTAL')
!!! CLEANUP: could we make this a non-blocking operation?
!!! Is this needed? it is on the initialization step.
      CALL MPP_UPDATE_DOMAINS(u, v, domain, gridtype=dgrid_ne, complete=&
&                       .true.)
      CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(isd,jsd,ied,jed,is,ie,js,je,npx,npy,npz, &
!$OMP       gridstruct,flagstruct,bd,u,v,uc,vc,nested,divg) &
!$OMP       private(ua,va)
      DO k=1,npz
        arg1 = flagstruct%nord .GT. 0
        CALL D2C_SETUP(u(isd, jsd, k), v(isd, jsd, k), ua, va, uc(isd, &
&                jsd, k), vc(isd, jsd, k), arg1, isd, ied, jsd, jed, is&
&                , ie, js, je, npx, npy, gridstruct%grid_type, &
&                gridstruct%nested, gridstruct%se_corner, gridstruct%&
&                sw_corner, gridstruct%ne_corner, gridstruct%nw_corner, &
&                gridstruct%rsin_u, gridstruct%rsin_v, gridstruct%cosa_s&
&                , gridstruct%rsin2)
        IF (nested) THEN
          CALL DIVERGENCE_CORNER_NEST(u(isd, jsd, k), v(isd, jsd, k), ua&
&                               , va, divg(isd, jsd, k), gridstruct, &
&                               flagstruct, bd)
        ELSE
          CALL DIVERGENCE_CORNER(u(isd, jsd, k), v(isd, jsd, k), ua, va&
&                          , divg(isd, jsd, k), gridstruct, flagstruct, &
&                          bd)
        END IF
      END DO
    END IF
    IF (flagstruct%hydrostatic) THEN
!$OMP parallel do default(none) shared(npz,is,ie,js,je,pkz,pkz_coarse)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            pkz_coarse(i, j, k) = pkz(i, j, k)
          END DO
        END DO
      END DO
    END IF
!! Nested grid: receive from parent grid
    IF (neststruct%nested) THEN
      IF (.NOT.ALLOCATED(q_buf)) THEN
        ALLOCATE(q_buf(ncnst))
      END IF
      CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 0, 0, npz, bd, &
&                        delp_buf)
      DO n=1,ncnst
        CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 0, 0, npz, bd, &
&                          q_buf(n))
      END DO
      CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 0, 0, npz, bd, &
&                        pt_buf)
      IF (flagstruct%hydrostatic) THEN
        CALL ALLOCATE_FV_NEST_BC_TYPE(pkz_bc, is, ie, js, je, isd, ied, &
&                               jsd, jed, npx, npy, npz, ng, 0, 0, 0, &
&                               .false.)
        CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 0, 0, npz, bd, &
&                          pkz_buf)
      ELSE
        CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 0, 0, npz, bd, &
&                          w_buf)
        CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 0, 0, npz, bd, &
&                          delz_buf)
      END IF
      CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 0, 1, npz, bd, &
&                        u_buf)
      CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 0, 1, npz, bd, &
&                        vc_buf)
      CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 1, 0, npz, bd, &
&                        v_buf)
      CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 1, 0, npz, bd, &
&                        uc_buf)
      CALL NESTED_GRID_BC_RECV(neststruct%nest_domain, 1, 1, npz, bd, &
&                        divg_buf)
    END IF
!! Coarse grid: send to child grids
    DO p=1,SIZE(child_grids)
      IF (child_grids(p)) THEN
        CALL NESTED_GRID_BC_SEND(delp, neststruct%nest_domain_all(p), 0&
&                          , 0)
        DO n=1,ncnst
          CALL NESTED_GRID_BC_SEND(q(:, :, :, n), neststruct%&
&                            nest_domain_all(p), 0, 0)
        END DO
        CALL NESTED_GRID_BC_SEND(pt, neststruct%nest_domain_all(p), 0, 0&
&                         )
        IF (flagstruct%hydrostatic) THEN
!Working with PKZ is more complicated since it is only defined on the interior of the grid.
          CALL NESTED_GRID_BC_SEND(pkz_coarse, neststruct%&
&                            nest_domain_all(p), 0, 0)
        ELSE
          CALL NESTED_GRID_BC_SEND(w, neststruct%nest_domain_all(p), 0, &
&                            0)
          CALL NESTED_GRID_BC_SEND(delz, neststruct%nest_domain_all(p), &
&                            0, 0)
        END IF
        CALL NESTED_GRID_BC_SEND(u, neststruct%nest_domain_all(p), 0, 1)
        CALL NESTED_GRID_BC_SEND(vc, neststruct%nest_domain_all(p), 0, 1&
&                         )
        CALL NESTED_GRID_BC_SEND(v, neststruct%nest_domain_all(p), 1, 0)
        CALL NESTED_GRID_BC_SEND(uc, neststruct%nest_domain_all(p), 1, 0&
&                         )
        CALL NESTED_GRID_BC_SEND(divg, neststruct%nest_domain_all(p), 1&
&                          , 1)
      END IF
    END DO
!Nested grid: do computations
    IF (nested) THEN
      CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct%&
&                             ind_h, neststruct%wt_h, 0, 0, npx, npy, &
&                             npz, bd, neststruct%delp_bc, delp_buf, &
&                             do_pd)
      DO n=1,ncnst
        CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct&
&                               %ind_h, neststruct%wt_h, 0, 0, npx, npy&
&                               , npz, bd, neststruct%q_bc(n), q_buf(n)&
&                               , do_pd)
      END DO
      CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct%&
&                             ind_h, neststruct%wt_h, 0, 0, npx, npy, &
&                             npz, bd, neststruct%pt_bc, pt_buf)
      sphum = GET_TRACER_INDEX(model_atmos, 'sphum')
      IF (flagstruct%hydrostatic) THEN
        CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct&
&                               %ind_h, neststruct%wt_h, 0, 0, npx, npy&
&                               , npz, bd, pkz_bc, pkz_buf)
        CALL SETUP_PT_BC(neststruct%pt_bc, pkz_bc, neststruct%q_bc(sphum&
&                  ), npx, npy, npz, zvir, bd)
      ELSE
        CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct&
&                               %ind_h, neststruct%wt_h, 0, 0, npx, npy&
&                               , npz, bd, neststruct%w_bc, w_buf)
        CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct&
&                               %ind_h, neststruct%wt_h, 0, 0, npx, npy&
&                               , npz, bd, neststruct%delz_bc, delz_buf)
!Need a negative-definite method?
        CALL SETUP_PT_NH_BC(neststruct%pt_bc, neststruct%delp_bc, &
&                     neststruct%delz_bc, neststruct%q_bc(sphum), &
&                     neststruct%q_bc, ncnst, npx, npy, npz, zvir, bd)
      END IF
      CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct%&
&                             ind_u, neststruct%wt_u, 0, 1, npx, npy, &
&                             npz, bd, neststruct%u_bc, u_buf)
      CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct%&
&                             ind_u, neststruct%wt_u, 0, 1, npx, npy, &
&                             npz, bd, neststruct%vc_bc, vc_buf)
      CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct%&
&                             ind_v, neststruct%wt_v, 1, 0, npx, npy, &
&                             npz, bd, neststruct%v_bc, v_buf)
      CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct%&
&                             ind_v, neststruct%wt_v, 1, 0, npx, npy, &
&                             npz, bd, neststruct%uc_bc, uc_buf)
      CALL NESTED_GRID_BC_SAVE_PROC(neststruct%nest_domain, neststruct%&
&                             ind_b, neststruct%wt_b, 1, 1, npx, npy, &
&                             npz, bd, neststruct%divg_bc, divg_buf)
    END IF
    IF (neststruct%first_step) THEN
      IF (neststruct%nested) CALL SET_BCS_T0(ncnst, flagstruct%&
&                                      hydrostatic, neststruct)
      neststruct%first_step = .false.
      IF (.NOT.flagstruct%hydrostatic) flagstruct%make_nh = .false.
    ELSE IF (flagstruct%make_nh) THEN
      IF (neststruct%nested) CALL SET_NH_BCS_T0(neststruct)
      flagstruct%make_nh = .false.
    END IF
!Unnecessary?
!!$    if ( neststruct%nested .and. .not. neststruct%divg_BC%initialized) then
!!$       neststruct%divg_BC%east_t0  = neststruct%divg_BC%east_t1
!!$       neststruct%divg_BC%west_t0  = neststruct%divg_BC%west_t1
!!$       neststruct%divg_BC%north_t0 = neststruct%divg_BC%north_t1
!!$       neststruct%divg_BC%south_t0 = neststruct%divg_BC%south_t1
!!$       neststruct%divg_BC%initialized = .true.
!!$    endif
    CALL MPP_SYNC_SELF()
  END SUBROUTINE SETUP_NESTED_GRID_BCS
  SUBROUTINE SETUP_PT_BC(pt_bc, pkz_bc, sphum_bc, npx, npy, npz, zvir, &
&   bd)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_NEST_BC_TYPE_3D), INTENT(IN), TARGET :: pkz_bc, sphum_bc
    TYPE(FV_NEST_BC_TYPE_3D), INTENT(INOUT), TARGET :: pt_bc
    INTEGER, INTENT(IN) :: npx, npy, npz
    REAL, INTENT(IN) :: zvir
    REAL, DIMENSION(:, :, :), POINTER :: ptbc, pkzbc, sphumbc
    INTEGER :: i, j, k, istart, iend
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
    IF (is .EQ. 1) THEN
      ptbc => pt_bc%west_t1
      pkzbc => pkz_bc%west_t1
      sphumbc => sphum_bc%west_t1
!$OMP parallel do default(none) shared(npz,jsd,jed,isd,ptBC,pkzBC,zvir,sphumBC)
      DO k=1,npz
        DO j=jsd,jed
          DO i=isd,0
            ptbc(i, j, k) = ptbc(i, j, k)/pkzbc(i, j, k)*(1.+zvir*&
&             sphumbc(i, j, k))
          END DO
        END DO
      END DO
    END IF
    IF (js .EQ. 1) THEN
      ptbc => pt_bc%south_t1
      pkzbc => pkz_bc%south_t1
      sphumbc => sphum_bc%south_t1
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
!$OMP parallel do default(none) shared(npz,jsd,istart,iend,ptBC,pkzBC,zvir,sphumBC)
      DO k=1,npz
        DO j=jsd,0
          DO i=istart,iend
            ptbc(i, j, k) = ptbc(i, j, k)/pkzbc(i, j, k)*(1.+zvir*&
&             sphumbc(i, j, k))
          END DO
        END DO
      END DO
    END IF
    IF (ie .EQ. npx - 1) THEN
      ptbc => pt_bc%east_t1
      pkzbc => pkz_bc%east_t1
      sphumbc => sphum_bc%east_t1
!$OMP parallel do default(none) shared(npz,jsd,jed,npx,ied,ptBC,pkzBC,zvir,sphumBC)
      DO k=1,npz
        DO j=jsd,jed
          DO i=npx,ied
            ptbc(i, j, k) = ptbc(i, j, k)/pkzbc(i, j, k)*(1.+zvir*&
&             sphumbc(i, j, k))
          END DO
        END DO
      END DO
    END IF
    IF (je .EQ. npy - 1) THEN
      ptbc => pt_bc%north_t1
      pkzbc => pkz_bc%north_t1
      sphumbc => sphum_bc%north_t1
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
!$OMP parallel do default(none) shared(npz,npy,jed,npx,istart,iend,ptBC,pkzBC,zvir,sphumBC)
      DO k=1,npz
        DO j=npy,jed
          DO i=istart,iend
            ptbc(i, j, k) = ptbc(i, j, k)/pkzbc(i, j, k)*(1.+zvir*&
&             sphumbc(i, j, k))
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE SETUP_PT_BC
  SUBROUTINE SETUP_PT_NH_BC(pt_bc, delp_bc, delz_bc, sphum_bc, q_bc, nq&
&   , npx, npy, npz, zvir, bd)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_NEST_BC_TYPE_3D), INTENT(IN), TARGET :: delp_bc, delz_bc, &
&   sphum_bc
    TYPE(FV_NEST_BC_TYPE_3D), INTENT(INOUT), TARGET :: pt_bc
    INTEGER, INTENT(IN) :: nq
    TYPE(FV_NEST_BC_TYPE_3D), INTENT(IN), TARGET :: q_bc(nq)
    INTEGER, INTENT(IN) :: npx, npy, npz
    REAL, INTENT(IN) :: zvir
! heat capacity of water at 0C
    REAL, PARAMETER :: c_liq=4185.5
! heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice)
    REAL, PARAMETER :: c_ice=1972.
! 1384.5
    REAL, PARAMETER :: cv_vap=cp_vapor-rvgas
    REAL, DIMENSION(:, :, :), POINTER :: ptbc, sphumbc, qconbc, delpbc, &
&   delzbc, cappabc
    REAL, DIMENSION(:, :, :), POINTER :: liq_watbc_west, ice_watbc_west&
&   , rainwatbc_west, snowwatbc_west, graupelbc_west
    REAL, DIMENSION(:, :, :), POINTER :: liq_watbc_east, ice_watbc_east&
&   , rainwatbc_east, snowwatbc_east, graupelbc_east
    REAL, DIMENSION(:, :, :), POINTER :: liq_watbc_north, &
&   ice_watbc_north, rainwatbc_north, snowwatbc_north, graupelbc_north
    REAL, DIMENSION(:, :, :), POINTER :: liq_watbc_south, &
&   ice_watbc_south, rainwatbc_south, snowwatbc_south, graupelbc_south
    REAL :: dp1, q_liq, q_sol
    REAL, SAVE :: q_con=0.
    REAL :: cvm
    REAL :: pkz
    REAL :: rdg
    REAL :: cv_air
    INTEGER :: i, j, k, istart, iend
    INTEGER :: liq_wat, ice_wat, rainwat, snowwat, graupel
! For GFS Partitioning
    REAL, PARAMETER :: tice=273.16
    REAL, PARAMETER :: t_i0=15.
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ALLOCATED
    INTRINSIC LOG
    INTRINSIC EXP
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    rdg = -(rdgas/grav)
    cv_air = cp_air - rdgas
    liq_wat = GET_TRACER_INDEX(model_atmos, 'liq_wat')
    ice_wat = GET_TRACER_INDEX(model_atmos, 'ice_wat')
    rainwat = GET_TRACER_INDEX(model_atmos, 'rainwat')
    snowwat = GET_TRACER_INDEX(model_atmos, 'snowwat')
    graupel = GET_TRACER_INDEX(model_atmos, 'graupel')
    IF (is .EQ. 1) THEN
      IF (.NOT.ALLOCATED(dum_west)) THEN
        ALLOCATE(dum_west(isd:0, jsd:jed, npz))
!$OMP parallel do default(none) shared(npz,isd,jsd,jed,dum_West)
        DO k=1,npz
          DO j=jsd,jed
            DO i=isd,0
              dum_west(i, j, k) = 0.
            END DO
          END DO
        END DO
      END IF
    END IF
    IF (js .EQ. 1) THEN
      IF (.NOT.ALLOCATED(dum_south)) THEN
        ALLOCATE(dum_south(isd:ied, jsd:0, npz))
!$OMP parallel do default(none) shared(npz,isd,ied,jsd,dum_South)
        DO k=1,npz
          DO j=jsd,0
            DO i=isd,ied
              dum_south(i, j, k) = 0.
            END DO
          END DO
        END DO
      END IF
    END IF
    IF (ie .EQ. npx - 1) THEN
      IF (.NOT.ALLOCATED(dum_east)) THEN
        ALLOCATE(dum_east(npx:ied, jsd:jed, npz))
!$OMP parallel do default(none) shared(npx,npz,ied,jsd,jed,dum_East)
        DO k=1,npz
          DO j=jsd,jed
            DO i=npx,ied
              dum_east(i, j, k) = 0.
            END DO
          END DO
        END DO
      END IF
    END IF
    IF (je .EQ. npy - 1) THEN
      IF (.NOT.ALLOCATED(dum_north)) THEN
        ALLOCATE(dum_north(isd:ied, npy:jed, npz))
!$OMP parallel do default(none) shared(npy,npz,isd,ied,jed,dum_North)
        DO k=1,npz
          DO j=npy,jed
            DO i=isd,ied
              dum_north(i, j, k) = 0.
            END DO
          END DO
        END DO
      END IF
    END IF
    IF (liq_wat .GT. 0) THEN
      liq_watbc_west => q_bc(liq_wat)%west_t1
      liq_watbc_east => q_bc(liq_wat)%east_t1
      liq_watbc_north => q_bc(liq_wat)%north_t1
      liq_watbc_south => q_bc(liq_wat)%south_t1
    ELSE
      liq_watbc_west => dum_west
      liq_watbc_east => dum_east
      liq_watbc_north => dum_north
      liq_watbc_south => dum_south
    END IF
    IF (ice_wat .GT. 0) THEN
      ice_watbc_west => q_bc(ice_wat)%west_t1
      ice_watbc_east => q_bc(ice_wat)%east_t1
      ice_watbc_north => q_bc(ice_wat)%north_t1
      ice_watbc_south => q_bc(ice_wat)%south_t1
    ELSE
      ice_watbc_west => dum_west
      ice_watbc_east => dum_east
      ice_watbc_north => dum_north
      ice_watbc_south => dum_south
    END IF
    IF (rainwat .GT. 0) THEN
      rainwatbc_west => q_bc(rainwat)%west_t1
      rainwatbc_east => q_bc(rainwat)%east_t1
      rainwatbc_north => q_bc(rainwat)%north_t1
      rainwatbc_south => q_bc(rainwat)%south_t1
    ELSE
      rainwatbc_west => dum_west
      rainwatbc_east => dum_east
      rainwatbc_north => dum_north
      rainwatbc_south => dum_south
    END IF
    IF (snowwat .GT. 0) THEN
      snowwatbc_west => q_bc(snowwat)%west_t1
      snowwatbc_east => q_bc(snowwat)%east_t1
      snowwatbc_north => q_bc(snowwat)%north_t1
      snowwatbc_south => q_bc(snowwat)%south_t1
    ELSE
      snowwatbc_west => dum_west
      snowwatbc_east => dum_east
      snowwatbc_north => dum_north
      snowwatbc_south => dum_south
    END IF
    IF (graupel .GT. 0) THEN
      graupelbc_west => q_bc(graupel)%west_t1
      graupelbc_east => q_bc(graupel)%east_t1
      graupelbc_north => q_bc(graupel)%north_t1
      graupelbc_south => q_bc(graupel)%south_t1
    ELSE
      graupelbc_west => dum_west
      graupelbc_east => dum_east
      graupelbc_north => dum_north
      graupelbc_south => dum_south
    END IF
    IF (is .EQ. 1) THEN
      ptbc => pt_bc%west_t1
      sphumbc => sphum_bc%west_t1
      delpbc => delp_bc%west_t1
      delzbc => delz_bc%west_t1
!$OMP parallel do default(none) shared(npz,jsd,jed,isd,zvir,sphumBC,liq_watBC_west,rainwatBC_west,ice_watBC_west,snowwatBC_west,g
!raupelBC_west,qconBC,cappaBC, &
!$OMP      rdg,cv_air,delpBC,delzBC,ptBC) &
!$OMP      private(dp1,q_con,q_liq,q_sol,cvm,pkz)
      DO k=1,npz
        DO j=jsd,jed
          DO i=isd,0
            dp1 = zvir*sphumbc(i, j, k)
            pkz = EXP(kappa*LOG(rdg*delpbc(i, j, k)*ptbc(i, j, k)*(1.+&
&             dp1)/delzbc(i, j, k)))
            ptbc(i, j, k) = ptbc(i, j, k)*(1.+dp1)/pkz
          END DO
        END DO
      END DO
    END IF
    IF (js .EQ. 1) THEN
      ptbc => pt_bc%south_t1
      sphumbc => sphum_bc%south_t1
      delpbc => delp_bc%south_t1
      delzbc => delz_bc%south_t1
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
!$OMP parallel do default(none) shared(npz,jsd,istart,iend,zvir,sphumBC, &
!$OMP      liq_watBC_south,rainwatBC_south,ice_watBC_south,&
!$OMP      snowwatBC_south,graupelBC_south,qconBC,cappaBC, &
!$OMP      rdg,cv_air,delpBC,delzBC,ptBC) &
!$OMP      private(dp1,q_con,q_liq,q_sol,cvm,pkz)
      DO k=1,npz
        DO j=jsd,0
          DO i=istart,iend
            dp1 = zvir*sphumbc(i, j, k)
            pkz = EXP(kappa*LOG(rdg*delpbc(i, j, k)*ptbc(i, j, k)*(1.+&
&             dp1)/delzbc(i, j, k)))
            ptbc(i, j, k) = ptbc(i, j, k)*(1.+dp1)/pkz
          END DO
        END DO
      END DO
    END IF
    IF (ie .EQ. npx - 1) THEN
      ptbc => pt_bc%east_t1
      sphumbc => sphum_bc%east_t1
      delpbc => delp_bc%east_t1
      delzbc => delz_bc%east_t1
!$OMP parallel do default(none) shared(npz,jsd,jed,npx,ied,zvir,sphumBC, &
!$OMP      liq_watBC_east,rainwatBC_east,ice_watBC_east,snowwatBC_east,graupelBC_east,qconBC,cappaBC, &
!$OMP      rdg,cv_air,delpBC,delzBC,ptBC) &
!$OMP      private(dp1,q_con,q_liq,q_sol,cvm,pkz)
      DO k=1,npz
        DO j=jsd,jed
          DO i=npx,ied
            dp1 = zvir*sphumbc(i, j, k)
            pkz = EXP(kappa*LOG(rdg*delpbc(i, j, k)*ptbc(i, j, k)*(1.+&
&             dp1)/delzbc(i, j, k)))
            ptbc(i, j, k) = ptbc(i, j, k)*(1.+dp1)/pkz
          END DO
        END DO
      END DO
    END IF
    IF (je .EQ. npy - 1) THEN
      ptbc => pt_bc%north_t1
      sphumbc => sphum_bc%north_t1
      delpbc => delp_bc%north_t1
      delzbc => delz_bc%north_t1
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
!$OMP parallel do default(none) shared(npz,npy,jed,istart,iend,zvir, &
!$OMP      sphumBC,liq_watBC_north,rainwatBC_north,ice_watBC_north,snowwatBC_north,graupelBC_north,qconBC,cappaBC, &
!$OMP      rdg,cv_air,delpBC,delzBC,ptBC) &
!$OMP      private(dp1,q_con,q_liq,q_sol,cvm,pkz)
      DO k=1,npz
        DO j=npy,jed
          DO i=istart,iend
            dp1 = zvir*sphumbc(i, j, k)
            pkz = EXP(kappa*LOG(rdg*delpbc(i, j, k)*ptbc(i, j, k)*(1.+&
&             dp1)/delzbc(i, j, k)))
            ptbc(i, j, k) = ptbc(i, j, k)*(1.+dp1)/pkz
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE SETUP_PT_NH_BC
  SUBROUTINE SET_BCS_T0(ncnst, hydrostatic, neststruct)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ncnst
    LOGICAL, INTENT(IN) :: hydrostatic
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    INTEGER :: n
    neststruct%delp_bc%east_t0 = neststruct%delp_bc%east_t1
    neststruct%delp_bc%west_t0 = neststruct%delp_bc%west_t1
    neststruct%delp_bc%north_t0 = neststruct%delp_bc%north_t1
    neststruct%delp_bc%south_t0 = neststruct%delp_bc%south_t1
    DO n=1,ncnst
      neststruct%q_bc(n)%east_t0 = neststruct%q_bc(n)%east_t1
      neststruct%q_bc(n)%west_t0 = neststruct%q_bc(n)%west_t1
      neststruct%q_bc(n)%north_t0 = neststruct%q_bc(n)%north_t1
      neststruct%q_bc(n)%south_t0 = neststruct%q_bc(n)%south_t1
    END DO
    neststruct%pt_bc%east_t0 = neststruct%pt_bc%east_t1
    neststruct%pt_bc%west_t0 = neststruct%pt_bc%west_t1
    neststruct%pt_bc%north_t0 = neststruct%pt_bc%north_t1
    neststruct%pt_bc%south_t0 = neststruct%pt_bc%south_t1
    neststruct%pt_bc%east_t0 = neststruct%pt_bc%east_t1
    neststruct%pt_bc%west_t0 = neststruct%pt_bc%west_t1
    neststruct%pt_bc%north_t0 = neststruct%pt_bc%north_t1
    neststruct%pt_bc%south_t0 = neststruct%pt_bc%south_t1
    IF (.NOT.hydrostatic) CALL SET_NH_BCS_T0(neststruct)
    neststruct%u_bc%east_t0 = neststruct%u_bc%east_t1
    neststruct%u_bc%west_t0 = neststruct%u_bc%west_t1
    neststruct%u_bc%north_t0 = neststruct%u_bc%north_t1
    neststruct%u_bc%south_t0 = neststruct%u_bc%south_t1
    neststruct%v_bc%east_t0 = neststruct%v_bc%east_t1
    neststruct%v_bc%west_t0 = neststruct%v_bc%west_t1
    neststruct%v_bc%north_t0 = neststruct%v_bc%north_t1
    neststruct%v_bc%south_t0 = neststruct%v_bc%south_t1
    neststruct%vc_bc%east_t0 = neststruct%vc_bc%east_t1
    neststruct%vc_bc%west_t0 = neststruct%vc_bc%west_t1
    neststruct%vc_bc%north_t0 = neststruct%vc_bc%north_t1
    neststruct%vc_bc%south_t0 = neststruct%vc_bc%south_t1
    neststruct%uc_bc%east_t0 = neststruct%uc_bc%east_t1
    neststruct%uc_bc%west_t0 = neststruct%uc_bc%west_t1
    neststruct%uc_bc%north_t0 = neststruct%uc_bc%north_t1
    neststruct%uc_bc%south_t0 = neststruct%uc_bc%south_t1
    neststruct%divg_bc%east_t0 = neststruct%divg_bc%east_t1
    neststruct%divg_bc%west_t0 = neststruct%divg_bc%west_t1
    neststruct%divg_bc%north_t0 = neststruct%divg_bc%north_t1
    neststruct%divg_bc%south_t0 = neststruct%divg_bc%south_t1
  END SUBROUTINE SET_BCS_T0
  SUBROUTINE SET_NH_BCS_T0(neststruct)
    IMPLICIT NONE
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    neststruct%delz_bc%east_t0 = neststruct%delz_bc%east_t1
    neststruct%delz_bc%west_t0 = neststruct%delz_bc%west_t1
    neststruct%delz_bc%north_t0 = neststruct%delz_bc%north_t1
    neststruct%delz_bc%south_t0 = neststruct%delz_bc%south_t1
    neststruct%w_bc%east_t0 = neststruct%w_bc%east_t1
    neststruct%w_bc%west_t0 = neststruct%w_bc%west_t1
    neststruct%w_bc%north_t0 = neststruct%w_bc%north_t1
    neststruct%w_bc%south_t0 = neststruct%w_bc%south_t1
  END SUBROUTINE SET_NH_BCS_T0
!! nestupdate types
!! 1 - Interpolation update on all variables
!! 2 - Conserving update (over areas on cell-
!!     centered variables, over faces on winds) on all variables
!! 3 - Interpolation update on winds only
!! 4 - Interpolation update on all variables except delp (mass conserving)
!! 5 - Remap interpolating update, delp not updated
!! 6 - Remap conserving update, delp not updated
!! 7 - Remap conserving update, delp and q not updated
!! 8 - Remap conserving update, only winds updated
!! Note that nestupdate > 3 will not update delp.
!! "Remap update" remaps updated variables from the nested grid's
!!  vertical coordinate to that of the coarse grid. When delp is not
!!  updated (nestbctype >= 3) the vertical coordinates differ on
!!  the two grids, because the surface pressure will be different
!!  on the two grids.
!! Note: "conserving updates" do not guarantee global conservation
!!  unless flux nested grid BCs are specified, or if a quantity is
!!  not updated at all. This ability has not been implemented.
  SUBROUTINE TWOWAY_NESTING(atm, ngrids, grids_on_this_pe, zvir)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ngrids
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: atm(ngrids)
    LOGICAL, INTENT(IN) :: grids_on_this_pe(ngrids)
    REAL, INTENT(IN) :: zvir
    INTEGER :: n, p, sphum
! ngrids > 1
    IF (ngrids .GT. 1) THEN
!loop backwards to allow information to propagate from finest to coarsest grids
      DO n=ngrids,2,-1
!two-way updating
        IF (atm(n)%neststruct%twowaynest) THEN
          IF (grids_on_this_pe(n) .OR. grids_on_this_pe(atm(n)%&
&             parent_grid%grid_number)) THEN
            sphum = GET_TRACER_INDEX(model_atmos, 'sphum')
            CALL TWOWAY_NEST_UPDATE(atm(n)%npx, atm(n)%npy, atm(n)%npz, &
&                             zvir, atm(n)%ncnst, sphum, atm(n)%u, atm(n&
&                             )%v, atm(n)%w, atm(n)%omga, atm(n)%pt, atm&
&                             (n)%delp, atm(n)%q, atm(n)%uc, atm(n)%vc, &
&                             atm(n)%pkz, atm(n)%delz, atm(n)%ps, atm(n)&
&                             %ptop, atm(n)%gridstruct, atm(n)%&
&                             flagstruct, atm(n)%neststruct, atm(n)%&
&                             parent_grid, atm(n)%bd, .false.)
          END IF
        END IF
      END DO
!NOTE: these routines need to be used with any grid which has been updated to, not just the coarsest grid.
      DO n=1,ngrids
        IF (atm(n)%neststruct%parent_of_twoway .AND. grids_on_this_pe(n)&
&       ) CALL AFTER_TWOWAY_NEST_UPDATE(atm(n)%npx, atm(n)%npy, atm(n)%&
&                                 npz, atm(n)%ng, atm(n)%ncnst, atm(n)%u&
&                                 , atm(n)%v, atm(n)%w, atm(n)%delz, atm&
&                                 (n)%pt, atm(n)%delp, atm(n)%q, atm(n)%&
&                                 ps, atm(n)%pe, atm(n)%pk, atm(n)%peln&
&                                 , atm(n)%pkz, atm(n)%phis, atm(n)%ua, &
&                                 atm(n)%va, atm(n)%ptop, atm(n)%&
&                                 gridstruct, atm(n)%flagstruct, atm(n)%&
&                                 domain, atm(n)%bd)
      END DO
    END IF
  END SUBROUTINE TWOWAY_NESTING
!!!CLEANUP: this routine assumes that the PARENT GRID has pt = (regular) temperature,
!!!not potential temperature; which may cause problems when updating if this is not the case.
  SUBROUTINE TWOWAY_NEST_UPDATE(npx, npy, npz, zvir, ncnst, sphum, u, v&
&   , w, omga, pt, delp, q, uc, vc, pkz, delz, ps, ptop, gridstruct, &
&   flagstruct, neststruct, parent_grid, bd, conv_theta_in)
    IMPLICIT NONE
    REAL, INTENT(IN) :: zvir, ptop
    INTEGER, INTENT(IN) :: npx, npy, npz
    INTEGER, INTENT(IN) :: ncnst, sphum
    LOGICAL, INTENT(IN), OPTIONAL :: conv_theta_in
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! D grid zonal wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
!  W (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:, bd%jsd:, :)
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! specific humidity and constituents
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst)
! (uc,vc) C grid winds
    REAL, INTENT(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz)
! finite-volume mean pk
    REAL, INTENT(INOUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:, bd%jsd:, :)
! Surface pressure (pascal)
    REAL, INTENT(INOUT) :: ps(bd%isd:bd%ied, bd%jsd:bd%jed)
    TYPE(FV_GRID_TYPE), INTENT(INOUT) :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(INOUT) :: flagstruct
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    REAL, ALLOCATABLE :: t_nest(:, :, :), ps0(:, :)
    INTEGER :: i, j, k, n
    INTEGER :: isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p
    INTEGER :: isg, ieg, jsg, jeg, npx_p, npy_p
    INTEGER :: istart, iend
    REAL :: qmass_b, qmass_a
    REAL, SAVE :: fix=1.
    LOGICAL :: used
    LOGICAL, SAVE :: conv_theta=.true.
    REAL :: qdp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, ALLOCATABLE :: qdp_coarse(:, :, :)
    REAL(kind=f_p), ALLOCATABLE :: q_diff(:, :, :)
    REAL :: l_sum_b(npz), l_sum_a(npz)
    INTEGER :: upoff
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: isu, ieu, jsu, jeu
    INTRINSIC PRESENT
    INTRINSIC ALLOCATED
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    isu = neststruct%isu
    ieu = neststruct%ieu
    jsu = neststruct%jsu
    jeu = neststruct%jeu
    upoff = neststruct%upoff
!We update actual temperature, not theta.
!If pt is actual temperature, set conv_theta to .false.
    IF (PRESENT(conv_theta_in)) conv_theta = conv_theta_in
    IF (.NOT.neststruct%parent_proc .AND. (.NOT.neststruct%child_proc)) &
&   THEN
      RETURN
    ELSE
      CALL MPP_GET_DATA_DOMAIN(parent_grid%domain, isd_p, ied_p, jsd_p, &
&                        jed_p)
      CALL MPP_GET_COMPUTE_DOMAIN(parent_grid%domain, isc_p, iec_p, &
&                           jsc_p, jec_p)
!delp/ps
      IF (neststruct%nestupdate .LT. 3) THEN
        CALL UPDATE_COARSE_GRID(parent_grid%delp, delp, neststruct%&
&                         nest_domain, neststruct%ind_update_h, &
&                         gridstruct%dx, gridstruct%dy, gridstruct%area&
&                         , isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, &
&                         jed, neststruct%isu, neststruct%ieu, &
&                         neststruct%jsu, neststruct%jeu, npx, npy, npz&
&                         , 0, 0, neststruct%refinement, neststruct%&
&                         nestupdate, upoff, 0, neststruct%parent_proc, &
&                         neststruct%child_proc, parent_grid)
!self
        CALL MPP_SYNC()
      END IF
!if (neststruct%nestupdate /= 3 .and. neststruct%nestbctype /= 3) then
      IF (neststruct%nestupdate .NE. 3 .AND. neststruct%nestupdate .NE. &
&         7 .AND. neststruct%nestupdate .NE. 8) THEN
        ALLOCATE(qdp_coarse(isd_p:ied_p, jsd_p:jed_p, npz))
        IF (parent_grid%flagstruct%nwat .GT. 0) THEN
          ALLOCATE(q_diff(isd_p:ied_p, jsd_p:jed_p, npz))
          q_diff = 0.
        END IF
        DO n=1,parent_grid%flagstruct%nwat
          qdp_coarse = 0.
          IF (neststruct%child_proc) THEN
            DO k=1,npz
              DO j=jsd,jed
                DO i=isd,ied
                  qdp(i, j, k) = q(i, j, k, n)*delp(i, j, k)
                END DO
              END DO
            END DO
          ELSE
            qdp = 0.
          END IF
          IF (neststruct%parent_proc) THEN
!Add up ONLY region being replaced by nested grid
            DO k=1,npz
              DO j=jsu,jeu
                DO i=isu,ieu
                  qdp_coarse(i, j, k) = parent_grid%q(i, j, k, n)*&
&                   parent_grid%delp(i, j, k)
                END DO
              END DO
            END DO
            CALL LEVEL_SUM(qdp_coarse, parent_grid%gridstruct%area, &
&                    parent_grid%domain, parent_grid%bd, npz, l_sum_b)
          ELSE
            qdp_coarse = 0.
          END IF
          IF (neststruct%parent_proc) THEN
            IF (n .LE. parent_grid%flagstruct%nwat) THEN
              DO k=1,npz
                DO j=jsu,jeu
                  DO i=isu,ieu
                    q_diff(i, j, k) = q_diff(i, j, k) - qdp_coarse(i, j&
&                     , k)
                  END DO
                END DO
              END DO
            END IF
          END IF
          CALL UPDATE_COARSE_GRID(qdp_coarse, qdp, neststruct%&
&                           nest_domain, neststruct%ind_update_h, &
&                           gridstruct%dx, gridstruct%dy, gridstruct%&
&                           area, isd_p, ied_p, jsd_p, jed_p, isd, ied, &
&                           jsd, jed, neststruct%isu, neststruct%ieu, &
&                           neststruct%jsu, neststruct%jeu, npx, npy, &
&                           npz, 0, 0, neststruct%refinement, neststruct&
&                           %nestupdate, upoff, 0, neststruct%&
&                           parent_proc, neststruct%child_proc, &
&                           parent_grid)
!self
          CALL MPP_SYNC()
          IF (neststruct%parent_proc) THEN
            CALL LEVEL_SUM(qdp_coarse, parent_grid%gridstruct%area, &
&                    parent_grid%domain, parent_grid%bd, npz, l_sum_a)
            DO k=1,npz
              IF (l_sum_a(k) .GT. 0.) THEN
                fix = l_sum_b(k)/l_sum_a(k)
                DO j=jsu,jeu
                  DO i=isu,ieu
!Normalization mass fixer
                    parent_grid%q(i, j, k, n) = qdp_coarse(i, j, k)*fix
                  END DO
                END DO
              END IF
            END DO
            IF (n .EQ. 1) sphum_ll_fix = 1. - fix
          END IF
          IF (neststruct%parent_proc) THEN
            IF (n .LE. parent_grid%flagstruct%nwat) THEN
              DO k=1,npz
                DO j=jsu,jeu
                  DO i=isu,ieu
                    q_diff(i, j, k) = q_diff(i, j, k) + parent_grid%q(i&
&                     , j, k, n)
                  END DO
                END DO
              END DO
            END IF
          END IF
        END DO
        IF (neststruct%parent_proc) THEN
          IF (parent_grid%flagstruct%nwat .GT. 0) THEN
            DO k=1,npz
              DO j=jsu,jeu
                DO i=isu,ieu
                  parent_grid%delp(i, j, k) = parent_grid%delp(i, j, k) &
&                   + q_diff(i, j, k)
                END DO
              END DO
            END DO
          END IF
          DO n=1,parent_grid%flagstruct%nwat
            DO k=1,npz
              DO j=jsu,jeu
                DO i=isu,ieu
                  parent_grid%q(i, j, k, n) = parent_grid%q(i, j, k, n)/&
&                   parent_grid%delp(i, j, k)
                END DO
              END DO
            END DO
          END DO
        END IF
        DEALLOCATE(qdp_coarse)
        IF (ALLOCATED(q_diff)) THEN
          DEALLOCATE(q_diff)
        END IF
      END IF
!Neststruct%nestupdate /= 3
      IF (neststruct%nestupdate .NE. 3 .AND. neststruct%nestupdate .NE. &
&         8) THEN
!conv_theta
        IF (conv_theta) THEN
          IF (neststruct%child_proc) THEN
!pt is potential temperature on the nested grid, but actual
!temperature on the coarse grid. Compute actual temperature
!on the nested grid, then gather.
            ALLOCATE(t_nest(isd:ied, jsd:jed, 1:npz))
!$OMP parallel do default(none) shared(npz,js,je,is,ie,t_nest,pt,pkz,zvir,q,sphum)
            DO k=1,npz
              DO j=js,je
                DO i=is,ie
                  t_nest(i, j, k) = pt(i, j, k)*pkz(i, j, k)/(1.+zvir*q(&
&                   i, j, k, sphum))
                END DO
              END DO
            END DO
            DEALLOCATE(t_nest)
          END IF
          CALL UPDATE_COARSE_GRID(parent_grid%pt, t_nest, neststruct%&
&                           nest_domain, neststruct%ind_update_h, &
&                           gridstruct%dx, gridstruct%dy, gridstruct%&
&                           area, isd_p, ied_p, jsd_p, jed_p, isd, ied, &
&                           jsd, jed, neststruct%isu, neststruct%ieu, &
&                           neststruct%jsu, neststruct%jeu, npx, npy, &
&                           npz, 0, 0, neststruct%refinement, neststruct&
&                           %nestupdate, upoff, 0, neststruct%&
&                           parent_proc, neststruct%child_proc, &
&                           parent_grid)
        ELSE
          CALL UPDATE_COARSE_GRID(parent_grid%pt, pt, neststruct%&
&                           nest_domain, neststruct%ind_update_h, &
&                           gridstruct%dx, gridstruct%dy, gridstruct%&
&                           area, isd_p, ied_p, jsd_p, jed_p, isd, ied, &
&                           jsd, jed, neststruct%isu, neststruct%ieu, &
&                           neststruct%jsu, neststruct%jeu, npx, npy, &
&                           npz, 0, 0, neststruct%refinement, neststruct&
&                           %nestupdate, upoff, 0, neststruct%&
&                           parent_proc, neststruct%child_proc, &
&                           parent_grid)
        END IF
!self
        CALL MPP_SYNC()
        IF (.NOT.flagstruct%hydrostatic) THEN
          CALL UPDATE_COARSE_GRID(parent_grid%w, w, neststruct%&
&                           nest_domain, neststruct%ind_update_h, &
&                           gridstruct%dx, gridstruct%dy, gridstruct%&
&                           area, isd_p, ied_p, jsd_p, jed_p, isd, ied, &
&                           jsd, jed, neststruct%isu, neststruct%ieu, &
&                           neststruct%jsu, neststruct%jeu, npx, npy, &
&                           npz, 0, 0, neststruct%refinement, neststruct&
&                           %nestupdate, upoff, 0, neststruct%&
&                           parent_proc, neststruct%child_proc, &
&                           parent_grid)
!Updating for delz not yet implemented; may be problematic
!!$            call update_coarse_grid(parent_grid%delz, delz, neststruct%nest_domain, &
!!$                 neststruct%ind_update_h, &
!!$                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 0, &
!!$                 neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc)
!self
          CALL MPP_SYNC()
        END IF
      END IF
      CALL UPDATE_COARSE_GRID(parent_grid%u, u, neststruct%nest_domain, &
&                       neststruct%ind_update_h, gridstruct%dx, &
&                       gridstruct%dy, gridstruct%area, isd_p, ied_p, &
&                       jsd_p, jed_p, isd, ied, jsd, jed, neststruct%isu&
&                       , neststruct%ieu, neststruct%jsu, neststruct%jeu&
&                       , npx, npy, npz, 0, 1, neststruct%refinement, &
&                       neststruct%nestupdate, upoff, 0, neststruct%&
&                       parent_proc, neststruct%child_proc, parent_grid)
      CALL UPDATE_COARSE_GRID(parent_grid%v, v, neststruct%nest_domain, &
&                       neststruct%ind_update_h, gridstruct%dx, &
&                       gridstruct%dy, gridstruct%area, isd_p, ied_p, &
&                       jsd_p, jed_p, isd, ied, jsd, jed, neststruct%isu&
&                       , neststruct%ieu, neststruct%jsu, neststruct%jeu&
&                       , npx, npy, npz, 1, 0, neststruct%refinement, &
&                       neststruct%nestupdate, upoff, 0, neststruct%&
&                       parent_proc, neststruct%child_proc, parent_grid)
!self
      CALL MPP_SYNC()
      IF (neststruct%nestupdate .GE. 5 .AND. npz .GT. 4) THEN
!Use PS0 from nested grid, NOT the full delp. Also we assume the same number of levels on both grids.
!PS0 should be initially set to be ps so that this routine does NOTHING outside of the update region
!Re-compute nested (AND COARSE) grid ps
        ALLOCATE(ps0(isd_p:ied_p, jsd_p:jed_p))
        IF (neststruct%parent_proc) THEN
          parent_grid%ps = parent_grid%ptop
!This loop appears to cause problems with OMP
!$OMP parallel do default(none) shared(npz,jsd_p,jed_p,isd_p,ied_p,parent_grid)
          DO j=jsd_p,jed_p
            DO k=1,npz
              DO i=isd_p,ied_p
                parent_grid%ps(i, j) = parent_grid%ps(i, j) + &
&                 parent_grid%delp(i, j, k)
              END DO
            END DO
          END DO
          ps0 = parent_grid%ps
        END IF
        IF (neststruct%child_proc) THEN
          ps = ptop
!$OMP parallel do default(none) shared(npz,jsd,jed,isd,ied,ps,delp)
          DO j=jsd,jed
            DO k=1,npz
              DO i=isd,ied
                ps(i, j) = ps(i, j) + delp(i, j, k)
              END DO
            END DO
          END DO
        END IF
        CALL UPDATE_COARSE_GRID(ps0, ps, neststruct%nest_domain, &
&                         neststruct%ind_update_h, gridstruct%dx, &
&                         gridstruct%dy, gridstruct%area, isd_p, ied_p, &
&                         jsd_p, jed_p, isd, ied, jsd, jed, neststruct%&
&                         isu, neststruct%ieu, neststruct%jsu, &
&                         neststruct%jeu, npx, npy, 0, 0, neststruct%&
&                         refinement, neststruct%nestupdate, upoff, 0, &
&                         neststruct%parent_proc, neststruct%child_proc&
&                         , parent_grid)
!!! The mpp version of update_coarse_grid does not return a consistent value of ps
!!! across PEs, as it does not go into the haloes of a given coarse-grid PE. This
!!! update_domains call takes care of the problem.
        IF (neststruct%parent_proc) THEN
          CALL MPP_UPDATE_DOMAINS(parent_grid%ps, parent_grid%domain, &
&                           complete=.true.)
          CALL MPP_UPDATE_DOMAINS(ps0, parent_grid%domain, complete=&
&                           .true.)
        END IF
!self
        CALL MPP_SYNC()
        IF (parent_grid%tile .EQ. neststruct%parent_tile) THEN
!neststruct%parent_proc
          IF (neststruct%parent_proc) THEN
!comment out if statement to always remap theta instead of t in the remap-update.
!(In LtE typically we use remap_t = .true.: remapping t is better (except in
!idealized simulations with a background uniform theta) since near the top
!boundary theta is exponential, which is hard to accurately interpolate with a spline
            IF (parent_grid%flagstruct%remap_option .NE. 0) THEN
!$OMP parallel do default(none) shared(npz,jsc_p,jec_p,isc_p,iec_p,parent_grid,zvir,sphum)
              DO k=1,npz
                DO j=jsc_p,jec_p
                  DO i=isc_p,iec_p
                    parent_grid%pt(i, j, k) = parent_grid%pt(i, j, k)/&
&                     parent_grid%pkz(i, j, k)*(1.+zvir*parent_grid%q(i&
&                     , j, k, sphum))
                  END DO
                END DO
              END DO
            END IF
            CALL UPDATE_REMAP_TQW(npz, parent_grid%ak, parent_grid%bk, &
&                           parent_grid%ps, parent_grid%delp, &
&                           parent_grid%pt, parent_grid%q, parent_grid%w&
&                           , parent_grid%flagstruct%hydrostatic, npz, &
&                           ps0, zvir, parent_grid%ptop, ncnst, &
&                           parent_grid%flagstruct%kord_tm, parent_grid%&
&                           flagstruct%kord_tr, parent_grid%flagstruct%&
&                           kord_wz, isc_p, iec_p, jsc_p, jec_p, isd_p, &
&                           ied_p, jsd_p, jed_p, .false.)
!neststruct%nestupdate < 7)
            IF (parent_grid%flagstruct%remap_option .NE. 0) THEN
!$OMP parallel do default(none) shared(npz,jsc_p,jec_p,isc_p,iec_p,parent_grid,zvir,sphum)
              DO k=1,npz
                DO j=jsc_p,jec_p
                  DO i=isc_p,iec_p
                    parent_grid%pt(i, j, k) = parent_grid%pt(i, j, k)*&
&                     parent_grid%pkz(i, j, k)/(1.+zvir*parent_grid%q(i&
&                     , j, k, sphum))
                  END DO
                END DO
              END DO
            END IF
            CALL UPDATE_REMAP_UV(npz, parent_grid%ak, parent_grid%bk, &
&                          parent_grid%ps, parent_grid%u, parent_grid%v&
&                          , npz, ps0, parent_grid%flagstruct%kord_mt, &
&                          isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, &
&                          jsd_p, jed_p, parent_grid%ptop)
          END IF
        END IF
        IF (ALLOCATED(ps0)) THEN
          DEALLOCATE(ps0)
        END IF
      END IF
    END IF
  END SUBROUTINE TWOWAY_NEST_UPDATE
  SUBROUTINE LEVEL_SUM(q, area, domain, bd, npz, l_sum)
    IMPLICIT NONE
!       L_sum(k) = mpp_global_sum(domain, qA, flags=BITWISE_EXACT_SUM)
!       L_sum(k) = mpp_global_sum(domain, qA, flags=BITWISE_EFP_SUM) ! doesn't work??
    INTEGER, INTENT(IN) :: npz
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: area(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(OUT) :: l_sum(npz)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER :: i, j, k, n
!(bd%is:bd%ie, bd%js:bd%je)
    REAL :: qa
    DO k=1,npz
      qa = 0.
      DO j=bd%js,bd%je
        DO i=bd%is,bd%ie
!qA(i,j) = q(i,j,k)*area(i,j)
          qa = qa + q(i, j, k)*area(i, j)
        END DO
      END DO
      CALL MP_REDUCE_SUM(qa)
      l_sum(k) = qa
    END DO
  END SUBROUTINE LEVEL_SUM
  SUBROUTINE AFTER_TWOWAY_NEST_UPDATE(npx, npy, npz, ng, ncnst, u, v, w&
&   , delz, pt, delp, q, ps, pe, pk, peln, pkz, phis, ua, va, ptop, &
&   gridstruct, flagstruct, domain, bd)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: ptop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    INTEGER, INTENT(IN) :: ncnst
! D grid zonal wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), INTENT(INOUT) &
&   :: u
! D grid meridional wind (m/s)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), INTENT(INOUT) &
&   :: v
!  W (m/s)
    REAL, INTENT(INOUT) :: w(bd%isd:, bd%jsd:, :)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! specific humidity and constituents
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(bd%isd:, bd%jsd:, :)
!-----------------------------------------------------------------------
! Auxilliary pressure arrays:
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
! Surface pressure (pascal)
    REAL, INTENT(INOUT) :: ps(bd%isd:bd%ied, bd%jsd:bd%jed)
! edge pressure (pascal)
    REAL, INTENT(INOUT) :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
! pe**cappa
    REAL, INTENT(INOUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
! ln(pe)
    REAL, INTENT(INOUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
! finite-volume mean pk
    REAL, INTENT(INOUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(INOUT) ::&
&   ua, va
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN) :: flagstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    LOGICAL :: bad_range
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
    CALL CUBED_TO_LATLON(u, v, ua, va, gridstruct, npx, npy, npz, 1, &
&                  gridstruct%grid_type, domain, gridstruct%nested, &
&                  flagstruct%c2l_ord, bd)
!To get coarse grid pkz, etc right after a two-way update so
!that it is consistent across a restart:
!(should only be called after doing such an update)
!! CLEANUP: move to twoway_nest_update??
!mountain argument not used
    CALL P_VAR(npz, is, ie, js, je, ptop, ptop_min, delp, delz, pt, ps, &
&        pe, peln, pk, pkz, kappa, q, ng, flagstruct%ncnst, gridstruct%&
&        area_64, 0., .false., .false., flagstruct%moist_phys, &
&        flagstruct%hydrostatic, flagstruct%nwat, domain, .false.)
    IF (flagstruct%range_warn) THEN
      CALL RANGE_CHECK('TA update', pt, is, ie, js, je, ng, npz, &
&                gridstruct%agrid, 130., 350., bad_range)
      CALL RANGE_CHECK('UA update', ua, is, ie, js, je, ng, npz, &
&                gridstruct%agrid, -220., 250., bad_range)
      CALL RANGE_CHECK('VA update', va, is, ie, js, je, ng, npz, &
&                gridstruct%agrid, -220., 220., bad_range)
      IF (.NOT.flagstruct%hydrostatic) CALL RANGE_CHECK('W update', w, &
&                                                 is, ie, js, je, ng, &
&                                                 npz, gridstruct%agrid&
&                                                 , -50., 100., &
&                                                 bad_range)
    END IF
  END SUBROUTINE AFTER_TWOWAY_NEST_UPDATE
!Routines for remapping (interpolated) nested-grid data to the coarse-grid's vertical coordinate.
!This does not yet do anything for the tracers
  SUBROUTINE UPDATE_REMAP_TQW(npz, ak, bk, ps, delp, t, q, w, &
&   hydrostatic, kmd, ps0, zvir, ptop, nq, kord_tm, kord_tr, kord_wz, is&
&   , ie, js, je, isd, ied, jsd, jed, do_q)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npz, kmd, nq, kord_tm, kord_tr, kord_wz
    REAL, INTENT(IN) :: zvir, ptop
    REAL, INTENT(IN) :: ak(npz+1), bk(npz+1)
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(IN) :: ps0
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(IN) :: ps
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: delp
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: t, w
    REAL, DIMENSION(isd:ied, jsd:jed, npz, nq), INTENT(INOUT) :: q
    LOGICAL, INTENT(IN) :: hydrostatic, do_q
! local:
    REAL, DIMENSION(is:ie, kmd) :: tp, qp
    REAL, DIMENSION(is:ie, kmd+1) :: pe0, pn0
    REAL, DIMENSION(is:ie, npz) :: qn1
    REAL, DIMENSION(is:ie, npz+1) :: pe1, pn1
    INTEGER :: i, j, k, iq
    INTRINSIC LOG
    INTRINSIC ABS
    INTEGER :: abs0
!$OMP parallel do default(none) shared(js,je,kmd,is,ie,ak,bk,ps0,q,npz,ptop,do_q,&
!$OMP          t,w,ps,nq,hydrostatic,kord_tm,kord_tr,kord_wz) &
!$OMP          private(pe0,pn0,pe1,pn1,qp,tp,qn1)
    DO j=js,je
      DO k=1,kmd+1
        DO i=is,ie
          pe0(i, k) = ak(k) + bk(k)*ps0(i, j)
          pn0(i, k) = LOG(pe0(i, k))
        END DO
      END DO
      DO k=1,kmd+1
        DO i=is,ie
          pe1(i, k) = ak(k) + bk(k)*ps(i, j)
          pn1(i, k) = LOG(pe1(i, k))
        END DO
      END DO
      IF (do_q) THEN
        DO iq=1,nq
          DO k=1,kmd
            DO i=is,ie
              qp(i, k) = q(i, j, k, iq)
            END DO
          END DO
          CALL MAPPM(kmd, pe0, qp, npz, pe1, qn1, is, ie, 0, kord_tr, &
&              ptop)
          DO k=1,npz
            DO i=is,ie
              q(i, j, k, iq) = qn1(i, k)
            END DO
          END DO
        END DO
      END IF
      DO k=1,kmd
        DO i=is,ie
          tp(i, k) = t(i, j, k)
        END DO
      END DO
      IF (kord_tm .GE. 0.) THEN
        abs0 = kord_tm
      ELSE
        abs0 = -kord_tm
      END IF
!Remap T using logp
      CALL MAPPM(kmd, pn0, tp, npz, pn1, qn1, is, ie, 1, abs0, ptop)
      DO k=1,npz
        DO i=is,ie
          t(i, j, k) = qn1(i, k)
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        DO k=1,kmd
          DO i=is,ie
            tp(i, k) = w(i, j, k)
          END DO
        END DO
!Remap w using p
!Using iv == -1 instead of -2
        CALL MAPPM(kmd, pe0, tp, npz, pe1, qn1, is, ie, -1, kord_wz, &
&            ptop)
        DO k=1,npz
          DO i=is,ie
            w(i, j, k) = qn1(i, k)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE UPDATE_REMAP_TQW
!remap_uv as-is remaps only a-grid velocities. A new routine has been written to handle staggered grids.
  SUBROUTINE UPDATE_REMAP_UV(npz, ak, bk, ps, u, v, kmd, ps0, kord_mt, &
&   is, ie, js, je, isd, ied, jsd, jed, ptop)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npz
    REAL, INTENT(IN) :: ak(npz+1), bk(npz+1)
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    REAL, INTENT(IN) :: ps(isd:ied, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed+1, npz), INTENT(INOUT) :: u
    REAL, DIMENSION(isd:ied+1, jsd:jed, npz), INTENT(INOUT) :: v
!
    INTEGER, INTENT(IN) :: kmd, kord_mt
    REAL, INTENT(IN) :: ptop
    REAL, INTENT(IN) :: ps0(isd:ied, jsd:jed)
!
! local:
    REAL, DIMENSION(is:ie+1, kmd+1) :: pe0
    REAL, DIMENSION(is:ie+1, npz+1) :: pe1
    REAL, DIMENSION(is:ie+1, kmd) :: qt
    REAL, DIMENSION(is:ie+1, npz) :: qn1
    INTEGER :: i, j, k
    INTEGER :: arg1
!------
! map u
!------
!$OMP parallel do default(none) shared(js,je,kmd,is,ie,ak,bk,ps,ps0,npz,u,ptop,kord_mt) &
!$OMP          private(pe0,pe1,qt,qn1)
    DO j=js,je+1
!------
! Data
!------
      DO k=1,kmd+1
        DO i=is,ie
          pe0(i, k) = ak(k) + bk(k)*0.5*(ps0(i, j)+ps0(i, j-1))
        END DO
      END DO
!------
! Model
!------
      DO k=1,kmd+1
        DO i=is,ie
          pe1(i, k) = ak(k) + bk(k)*0.5*(ps(i, j)+ps(i, j-1))
        END DO
      END DO
!------
!Do map
!------
      qt = 0.
      DO k=1,kmd
        DO i=is,ie
          qt(i, k) = u(i, j, k)
        END DO
      END DO
      qn1 = 0.
      CALL MAPPM(kmd, pe0(is:ie, :), qt(is:ie, :), npz, pe1(is:ie, :), &
&          qn1(is:ie, :), is, ie, -1, kord_mt, ptop)
      DO k=1,npz
        DO i=is,ie
          u(i, j, k) = qn1(i, k)
        END DO
      END DO
    END DO
!------
! map v
!------
!$OMP parallel do default(none) shared(js,je,kmd,is,ie,ak,bk,ps,ps0,npz,v,ptop) &
!$OMP          private(pe0,pe1,qt,qn1)
    DO j=js,je
!------
! Data
!------
      DO k=1,kmd+1
        DO i=is,ie+1
          pe0(i, k) = ak(k) + bk(k)*0.5*(ps0(i, j)+ps0(i-1, j))
        END DO
      END DO
!------
! Model
!------
      DO k=1,kmd+1
        DO i=is,ie+1
          pe1(i, k) = ak(k) + bk(k)*0.5*(ps(i, j)+ps(i-1, j))
        END DO
      END DO
!------
!Do map
!------
      qt = 0.
      DO k=1,kmd
        DO i=is,ie+1
          qt(i, k) = v(i, j, k)
        END DO
      END DO
      qn1 = 0.
      arg1 = ie + 1
      CALL MAPPM(kmd, pe0(is:ie+1, :), qt(is:ie+1, :), npz, pe1(is:ie+1&
&          , :), qn1(is:ie+1, :), is, arg1, -1, 8, ptop)
      DO k=1,npz
        DO i=is,ie+1
          v(i, j, k) = qn1(i, k)
        END DO
      END DO
    END DO
  END SUBROUTINE UPDATE_REMAP_UV
end module fv_nesting_adm_mod
