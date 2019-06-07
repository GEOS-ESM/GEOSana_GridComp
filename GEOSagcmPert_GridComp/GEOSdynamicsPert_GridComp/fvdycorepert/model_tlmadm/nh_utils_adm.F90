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
module nh_utils_adm_mod
! Developer: S.-J. Lin, NOAA/GFDL
! To do list:
! include moisture effect in pt
!------------------------------
   use constants_mod,     only: rdgas, cp_air, grav
   use tp_core_adm_mod,   only: fv_tp_2d
   use tp_core_adm_mod,   only: fv_tp_2d_fwd, fv_tp_2d_bwd, fv_tp_2d_adm
   use sw_core_adm_mod,   only: fill_4corners, del6_vt_flux
   use sw_core_adm_mod,   only: fill_4corners_fwd, fill_4corners_bwd, del6_vt_flux_adm
   use fv_arrays_mod,     only: fv_grid_bounds_type, fv_grid_type

   use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                            pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

   implicit none
   private

   public update_dz_c, update_dz_d, nest_halo_nh
   public sim_solver, sim1_solver, sim3_solver
   public sim3p0_solver, rim_2d
   public Riem_Solver_c
   public update_dz_c_fwd, update_dz_d_fwd, nest_halo_nh_fwd
   public sim_solver_fwd, sim1_solver_fwd, sim3_solver_fwd
   public sim3p0_solver_fwd, rim_2d_fwd
   public Riem_Solver_c_fwd
   public update_dz_c_bwd, update_dz_d_bwd, nest_halo_nh_bwd
   public sim_solver_bwd, sim1_solver_bwd, sim3_solver_bwd
   public sim3p0_solver_bwd, rim_2d_bwd
   public Riem_Solver_c_bwd

   real, parameter:: dz_min = 2.
   real, parameter:: r3 = 1./3.

CONTAINS
!  Differentiation of update_dz_c in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: ws gz ut vt
!   with respect to varying inputs: ws gz ut vt
  SUBROUTINE UPDATE_DZ_C_FWD(is, ie, js, je, km, ng, dt, dp0, zs, area, &
&   ut, vt, gz, ws, npx, npy, sw_corner, se_corner, ne_corner, nw_corner&
&   , bd, grid_type)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km, npx, npy, grid_type
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: dp0(km)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: ut, vt
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng), INTENT(IN) :: area
    REAL, INTENT(INOUT) :: gz(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(IN) :: zs(is-ng:ie+ng, js-ng:je+ng)
    REAL :: ws(is-ng:ie+ng, js-ng:je+ng)
! Local Work array:
    REAL :: gz2(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-1:ie+2, js-1:je+1) :: xfx, fx
    REAL, DIMENSION(is-1:ie+1, js-1:je+2) :: yfx, fy
    REAL, PARAMETER :: r14=1./14.
    INTEGER :: i, j, k
    INTEGER :: is1, ie1, js1, je1
    INTEGER :: ie2, je2
    REAL :: rdt, top_ratio, bot_ratio, int_ratio
    INTRINSIC MAX

    gz2 = 0.0
    xfx = 0.0
    fx = 0.0
    yfx = 0.0
    fy = 0.0
    is1 = 0
    ie1 = 0
    js1 = 0
    je1 = 0
    ie2 = 0
    je2 = 0
    rdt = 0.0
    top_ratio = 0.0
    bot_ratio = 0.0
    int_ratio = 0.0

!--------------------------------------------------------------------
    rdt = 1./dt
    top_ratio = dp0(1)/(dp0(1)+dp0(2))
    bot_ratio = dp0(km)/(dp0(km-1)+dp0(km))
    is1 = is - 1
    js1 = js - 1
    ie1 = ie + 1
    je1 = je + 1
    ie2 = ie + 2
    je2 = je + 2
!$OMP parallel do default(none) shared(js1,je1,is1,ie2,km,je2,ie1,ut,top_ratio,vt, &
!$OMP                                  bot_ratio,dp0,js,je,ng,is,ie,gz,grid_type,  &
!$OMP                                  bd,npx,npy,sw_corner,se_corner,ne_corner,   &
!$OMP                                  nw_corner,area) &
!$OMP                          private(gz2, xfx, yfx, fx, fy, int_ratio)
    DO k=1,km+1
      IF (k .EQ. 1) THEN
        DO j=js1,je1
          DO i=is1,ie2
            CALL PUSHREALARRAY(xfx(i, j))
            xfx(i, j) = ut(i, j, 1) + (ut(i, j, 1)-ut(i, j, 2))*&
&             top_ratio
          END DO
        END DO
        DO j=js1,je2
          DO i=is1,ie1
            CALL PUSHREALARRAY(yfx(i, j))
            yfx(i, j) = vt(i, j, 1) + (vt(i, j, 1)-vt(i, j, 2))*&
&             top_ratio
          END DO
        END DO
        CALL PUSHCONTROL(2,2)
      ELSE IF (k .EQ. km + 1) THEN
! Bottom extrapolation
        DO j=js1,je1
          DO i=is1,ie2
            CALL PUSHREALARRAY(xfx(i, j))
            xfx(i, j) = ut(i, j, km) + (ut(i, j, km)-ut(i, j, km-1))*&
&             bot_ratio
          END DO
        END DO
!             xfx(i,j) = r14*(3.*ut(i,j,km-2)-13.*ut(i,j,km-1)+24.*ut(i,j,km))
!             if ( xfx(i,j)*ut(i,j,km)<0. ) xfx(i,j) = 0.
        DO j=js1,je2
          DO i=is1,ie1
            CALL PUSHREALARRAY(yfx(i, j))
            yfx(i, j) = vt(i, j, km) + (vt(i, j, km)-vt(i, j, km-1))*&
&             bot_ratio
          END DO
        END DO
        CALL PUSHCONTROL(2,1)
      ELSE
!             yfx(i,j) = r14*(3.*vt(i,j,km-2)-13.*vt(i,j,km-1)+24.*vt(i,j,km))
!             if ( yfx(i,j)*vt(i,j,km)<0. ) yfx(i,j) = 0.
        CALL PUSHREALARRAY(int_ratio)
        int_ratio = 1./(dp0(k-1)+dp0(k))
        DO j=js1,je1
          DO i=is1,ie2
            CALL PUSHREALARRAY(xfx(i, j))
            xfx(i, j) = (dp0(k)*ut(i, j, k-1)+dp0(k-1)*ut(i, j, k))*&
&             int_ratio
          END DO
        END DO
        DO j=js1,je2
          DO i=is1,ie1
            CALL PUSHREALARRAY(yfx(i, j))
            yfx(i, j) = (dp0(k)*vt(i, j, k-1)+dp0(k-1)*vt(i, j, k))*&
&             int_ratio
          END DO
        END DO
        CALL PUSHCONTROL(2,0)
      END IF
      DO j=js-ng,je+ng
        DO i=is-ng,ie+ng
          CALL PUSHREALARRAY(gz2(i, j))
          gz2(i, j) = gz(i, j, k)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_FWD(gz2, 1, bd, npx, npy, sw_corner, &
&                        se_corner, ne_corner, nw_corner)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
      DO j=js1,je1
        DO i=is1,ie2
          IF (xfx(i, j) .GT. 0.) THEN
            CALL PUSHREALARRAY(fx(i, j))
            fx(i, j) = gz2(i-1, j)
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHREALARRAY(fx(i, j))
            fx(i, j) = gz2(i, j)
            CALL PUSHCONTROL(1,1)
          END IF
          CALL PUSHREALARRAY(fx(i, j))
          fx(i, j) = xfx(i, j)*fx(i, j)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_FWD(gz2, 2, bd, npx, npy, sw_corner, &
&                        se_corner, ne_corner, nw_corner)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
      DO j=js1,je2
        DO i=is1,ie1
          IF (yfx(i, j) .GT. 0.) THEN
            CALL PUSHREALARRAY(fy(i, j))
            fy(i, j) = gz2(i, j-1)
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHREALARRAY(fy(i, j))
            fy(i, j) = gz2(i, j)
            CALL PUSHCONTROL(1,1)
          END IF
          CALL PUSHREALARRAY(fy(i, j))
          fy(i, j) = yfx(i, j)*fy(i, j)
        END DO
      END DO
      DO j=js1,je1
        DO i=is1,ie1
          CALL PUSHREALARRAY(gz(i, j, k))
          gz(i, j, k) = (gz2(i, j)*area(i, j)+(fx(i, j)-fx(i+1, j))+(fy(&
&           i, j)-fy(i, j+1)))/(area(i, j)+(xfx(i, j)-xfx(i+1, j))+(yfx(&
&           i, j)-yfx(i, j+1)))
        END DO
      END DO
    END DO
! Enforce monotonicity of height to prevent blowup
!$OMP parallel do default(none) shared(is1,ie1,js1,je1,ws,zs,gz,rdt,km)
    DO j=js1,je1
      DO i=is1,ie1
        CALL PUSHREALARRAY(ws(i, j))
        ws(i, j) = (zs(i, j)-gz(i, j, km+1))*rdt
      END DO
      DO k=km,1,-1
        DO i=is1,ie1
          IF (gz(i, j, k) .LT. gz(i, j, k+1) + dz_min) THEN
            CALL PUSHREALARRAY(gz(i, j, k))
            gz(i, j, k) = gz(i, j, k+1) + dz_min
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHREALARRAY(gz(i, j, k))
            gz(i, j, k) = gz(i, j, k)
            CALL PUSHCONTROL(1,1)
          END IF
        END DO
      END DO
    END DO
    CALL PUSHINTEGER(ie2)
    CALL PUSHINTEGER(ie1)
    CALL PUSHREALARRAY(fy, (ie-is+3)*(je-js+4))
    CALL PUSHREALARRAY(fx, (ie-is+4)*(je-js+3))
    CALL PUSHINTEGER(js1)
    CALL PUSHREALARRAY(rdt)
    CALL PUSHREALARRAY(yfx, (ie-is+3)*(je-js+4))
    CALL PUSHINTEGER(je2)
    CALL PUSHINTEGER(je1)
    CALL PUSHREALARRAY(int_ratio)
    CALL PUSHREALARRAY(top_ratio)
    CALL PUSHREALARRAY(bot_ratio)
    CALL PUSHINTEGER(is1)
    CALL PUSHREALARRAY(gz2, (ie+2*ng-is+1)*(je+2*ng-js+1))
    CALL PUSHREALARRAY(xfx, (ie-is+4)*(je-js+3))
  END SUBROUTINE UPDATE_DZ_C_FWD
!  Differentiation of update_dz_c in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: ws gz ut vt
!   with respect to varying inputs: ws gz ut vt
  SUBROUTINE UPDATE_DZ_C_BWD(is, ie, js, je, km, ng, dt, dp0, zs, area, &
&   ut, ut_ad, vt, vt_ad, gz, gz_ad, ws, ws_ad, npx, npy, sw_corner, &
&   se_corner, ne_corner, nw_corner, bd, grid_type)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km, npx, npy, grid_type
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: dp0(km)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: ut, vt
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km) :: ut_ad, vt_ad
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng), INTENT(IN) :: area
    REAL, INTENT(INOUT) :: gz(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(INOUT) :: gz_ad(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(IN) :: zs(is-ng:ie+ng, js-ng:je+ng)
    REAL :: ws(is-ng:ie+ng, js-ng:je+ng)
    REAL :: ws_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL :: gz2(is-ng:ie+ng, js-ng:je+ng)
    REAL :: gz2_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-1:ie+2, js-1:je+1) :: xfx, fx
    REAL, DIMENSION(is-1:ie+2, js-1:je+1) :: xfx_ad, fx_ad
    REAL, DIMENSION(is-1:ie+1, js-1:je+2) :: yfx, fy
    REAL, DIMENSION(is-1:ie+1, js-1:je+2) :: yfx_ad, fy_ad
    REAL, PARAMETER :: r14=1./14.
    INTEGER :: i, j, k
    INTEGER :: is1, ie1, js1, je1
    INTEGER :: ie2, je2
    REAL :: rdt, top_ratio, bot_ratio, int_ratio
    INTRINSIC MAX
    REAL :: temp
    REAL :: temp_ad
    REAL :: temp_ad0
    INTEGER :: branch

    gz2 = 0.0
    xfx = 0.0
    fx = 0.0
    yfx = 0.0
    fy = 0.0
    is1 = 0
    ie1 = 0
    js1 = 0
    je1 = 0
    ie2 = 0
    je2 = 0
    rdt = 0.0
    top_ratio = 0.0
    bot_ratio = 0.0
    int_ratio = 0.0
    branch = 0

    CALL POPREALARRAY(xfx, (ie-is+4)*(je-js+3))
    CALL POPREALARRAY(gz2, (ie+2*ng-is+1)*(je+2*ng-js+1))
    CALL POPINTEGER(is1)
    CALL POPREALARRAY(bot_ratio)
    CALL POPREALARRAY(top_ratio)
    CALL POPREALARRAY(int_ratio)
    CALL POPINTEGER(je1)
    CALL POPINTEGER(je2)
    CALL POPREALARRAY(yfx, (ie-is+3)*(je-js+4))
    CALL POPREALARRAY(rdt)
    CALL POPINTEGER(js1)
    CALL POPREALARRAY(fx, (ie-is+4)*(je-js+3))
    CALL POPREALARRAY(fy, (ie-is+3)*(je-js+4))
    CALL POPINTEGER(ie1)
    CALL POPINTEGER(ie2)
    DO j=je1,js1,-1
      DO k=1,km,1
        DO i=ie1,is1,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(gz(i, j, k))
            gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + gz_ad(i, j, k)
            gz_ad(i, j, k) = 0.0
          ELSE
            CALL POPREALARRAY(gz(i, j, k))
          END IF
        END DO
      END DO
      DO i=ie1,is1,-1
        CALL POPREALARRAY(ws(i, j))
        gz_ad(i, j, km+1) = gz_ad(i, j, km+1) - rdt*ws_ad(i, j)
        ws_ad(i, j) = 0.0
      END DO
    END DO
    xfx_ad = 0.0
    gz2_ad = 0.0
    yfx_ad = 0.0
    fx_ad = 0.0
    fy_ad = 0.0
    DO k=km+1,1,-1
      DO j=je1,js1,-1
        DO i=ie1,is1,-1
          CALL POPREALARRAY(gz(i, j, k))
          temp = area(i, j) + xfx(i, j) - xfx(i+1, j) + yfx(i, j) - yfx(&
&           i, j+1)
          temp_ad = gz_ad(i, j, k)/temp
          temp_ad0 = -((area(i, j)*gz2(i, j)+fx(i, j)-fx(i+1, j)+fy(i, j&
&           )-fy(i, j+1))*temp_ad/temp)
          gz2_ad(i, j) = gz2_ad(i, j) + area(i, j)*temp_ad
          fx_ad(i, j) = fx_ad(i, j) + temp_ad
          fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad
          fy_ad(i, j) = fy_ad(i, j) + temp_ad
          fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad
          xfx_ad(i, j) = xfx_ad(i, j) + temp_ad0
          xfx_ad(i+1, j) = xfx_ad(i+1, j) - temp_ad0
          yfx_ad(i, j) = yfx_ad(i, j) + temp_ad0
          yfx_ad(i, j+1) = yfx_ad(i, j+1) - temp_ad0
          gz_ad(i, j, k) = 0.0
        END DO
      END DO
      DO j=je2,js1,-1
        DO i=ie1,is1,-1
          CALL POPREALARRAY(fy(i, j))
          yfx_ad(i, j) = yfx_ad(i, j) + fy(i, j)*fy_ad(i, j)
          fy_ad(i, j) = yfx(i, j)*fy_ad(i, j)
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(fy(i, j))
            gz2_ad(i, j-1) = gz2_ad(i, j-1) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
          ELSE
            CALL POPREALARRAY(fy(i, j))
            gz2_ad(i, j) = gz2_ad(i, j) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
          END IF
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) CALL FILL_4CORNERS_BWD(gz2, gz2_ad, 2, bd, npx&
&                                         , npy, sw_corner, se_corner, &
&                                         ne_corner, nw_corner)
      DO j=je1,js1,-1
        DO i=ie2,is1,-1
          CALL POPREALARRAY(fx(i, j))
          xfx_ad(i, j) = xfx_ad(i, j) + fx(i, j)*fx_ad(i, j)
          fx_ad(i, j) = xfx(i, j)*fx_ad(i, j)
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(fx(i, j))
            gz2_ad(i-1, j) = gz2_ad(i-1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0
          ELSE
            CALL POPREALARRAY(fx(i, j))
            gz2_ad(i, j) = gz2_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0
          END IF
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) CALL FILL_4CORNERS_BWD(gz2, gz2_ad, 1, bd, npx&
&                                         , npy, sw_corner, se_corner, &
&                                         ne_corner, nw_corner)
      DO j=je+ng,js-ng,-1
        DO i=ie+ng,is-ng,-1
          CALL POPREALARRAY(gz2(i, j))
          gz_ad(i, j, k) = gz_ad(i, j, k) + gz2_ad(i, j)
          gz2_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        DO j=je2,js1,-1
          DO i=ie1,is1,-1
            CALL POPREALARRAY(yfx(i, j))
            vt_ad(i, j, k-1) = vt_ad(i, j, k-1) + int_ratio*dp0(k)*&
&             yfx_ad(i, j)
            vt_ad(i, j, k) = vt_ad(i, j, k) + int_ratio*dp0(k-1)*yfx_ad(&
&             i, j)
            yfx_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je1,js1,-1
          DO i=ie2,is1,-1
            CALL POPREALARRAY(xfx(i, j))
            ut_ad(i, j, k-1) = ut_ad(i, j, k-1) + int_ratio*dp0(k)*&
&             xfx_ad(i, j)
            ut_ad(i, j, k) = ut_ad(i, j, k) + int_ratio*dp0(k-1)*xfx_ad(&
&             i, j)
            xfx_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPREALARRAY(int_ratio)
      ELSE IF (branch .EQ. 1) THEN
        DO j=je2,js1,-1
          DO i=ie1,is1,-1
            CALL POPREALARRAY(yfx(i, j))
            vt_ad(i, j, km) = vt_ad(i, j, km) + (bot_ratio+1.0)*yfx_ad(i&
&             , j)
            vt_ad(i, j, km-1) = vt_ad(i, j, km-1) - bot_ratio*yfx_ad(i, &
&             j)
            yfx_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je1,js1,-1
          DO i=ie2,is1,-1
            CALL POPREALARRAY(xfx(i, j))
            ut_ad(i, j, km) = ut_ad(i, j, km) + (bot_ratio+1.0)*xfx_ad(i&
&             , j)
            ut_ad(i, j, km-1) = ut_ad(i, j, km-1) - bot_ratio*xfx_ad(i, &
&             j)
            xfx_ad(i, j) = 0.0
          END DO
        END DO
      ELSE
        DO j=je2,js1,-1
          DO i=ie1,is1,-1
            CALL POPREALARRAY(yfx(i, j))
            vt_ad(i, j, 1) = vt_ad(i, j, 1) + (top_ratio+1.0)*yfx_ad(i, &
&             j)
            vt_ad(i, j, 2) = vt_ad(i, j, 2) - top_ratio*yfx_ad(i, j)
            yfx_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je1,js1,-1
          DO i=ie2,is1,-1
            CALL POPREALARRAY(xfx(i, j))
            ut_ad(i, j, 1) = ut_ad(i, j, 1) + (top_ratio+1.0)*xfx_ad(i, &
&             j)
            ut_ad(i, j, 2) = ut_ad(i, j, 2) - top_ratio*xfx_ad(i, j)
            xfx_ad(i, j) = 0.0
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE UPDATE_DZ_C_BWD
  SUBROUTINE UPDATE_DZ_C(is, ie, js, je, km, ng, dt, dp0, zs, area, ut, &
&   vt, gz, ws, npx, npy, sw_corner, se_corner, ne_corner, nw_corner, bd&
&   , grid_type)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km, npx, npy, grid_type
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: dp0(km)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: ut, vt
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng), INTENT(IN) :: area
    REAL, INTENT(INOUT) :: gz(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(IN) :: zs(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(OUT) :: ws(is-ng:ie+ng, js-ng:je+ng)
! Local Work array:
    REAL :: gz2(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-1:ie+2, js-1:je+1) :: xfx, fx
    REAL, DIMENSION(is-1:ie+1, js-1:je+2) :: yfx, fy
    REAL, PARAMETER :: r14=1./14.
    INTEGER :: i, j, k
    INTEGER :: is1, ie1, js1, je1
    INTEGER :: ie2, je2
    REAL :: rdt, top_ratio, bot_ratio, int_ratio
    INTRINSIC MAX
!--------------------------------------------------------------------
    rdt = 1./dt
    top_ratio = dp0(1)/(dp0(1)+dp0(2))
    bot_ratio = dp0(km)/(dp0(km-1)+dp0(km))
    is1 = is - 1
    js1 = js - 1
    ie1 = ie + 1
    je1 = je + 1
    ie2 = ie + 2
    je2 = je + 2
!$OMP parallel do default(none) shared(js1,je1,is1,ie2,km,je2,ie1,ut,top_ratio,vt, &
!$OMP                                  bot_ratio,dp0,js,je,ng,is,ie,gz,grid_type,  &
!$OMP                                  bd,npx,npy,sw_corner,se_corner,ne_corner,   &
!$OMP                                  nw_corner,area) &
!$OMP                          private(gz2, xfx, yfx, fx, fy, int_ratio)
    DO k=1,km+1
      IF (k .EQ. 1) THEN
        DO j=js1,je1
          DO i=is1,ie2
            xfx(i, j) = ut(i, j, 1) + (ut(i, j, 1)-ut(i, j, 2))*&
&             top_ratio
          END DO
        END DO
        DO j=js1,je2
          DO i=is1,ie1
            yfx(i, j) = vt(i, j, 1) + (vt(i, j, 1)-vt(i, j, 2))*&
&             top_ratio
          END DO
        END DO
      ELSE IF (k .EQ. km + 1) THEN
! Bottom extrapolation
        DO j=js1,je1
          DO i=is1,ie2
            xfx(i, j) = ut(i, j, km) + (ut(i, j, km)-ut(i, j, km-1))*&
&             bot_ratio
          END DO
        END DO
!             xfx(i,j) = r14*(3.*ut(i,j,km-2)-13.*ut(i,j,km-1)+24.*ut(i,j,km))
!             if ( xfx(i,j)*ut(i,j,km)<0. ) xfx(i,j) = 0.
        DO j=js1,je2
          DO i=is1,ie1
            yfx(i, j) = vt(i, j, km) + (vt(i, j, km)-vt(i, j, km-1))*&
&             bot_ratio
          END DO
        END DO
      ELSE
!             yfx(i,j) = r14*(3.*vt(i,j,km-2)-13.*vt(i,j,km-1)+24.*vt(i,j,km))
!             if ( yfx(i,j)*vt(i,j,km)<0. ) yfx(i,j) = 0.
        int_ratio = 1./(dp0(k-1)+dp0(k))
        DO j=js1,je1
          DO i=is1,ie2
            xfx(i, j) = (dp0(k)*ut(i, j, k-1)+dp0(k-1)*ut(i, j, k))*&
&             int_ratio
          END DO
        END DO
        DO j=js1,je2
          DO i=is1,ie1
            yfx(i, j) = (dp0(k)*vt(i, j, k-1)+dp0(k-1)*vt(i, j, k))*&
&             int_ratio
          END DO
        END DO
      END IF
      DO j=js-ng,je+ng
        DO i=is-ng,ie+ng
          gz2(i, j) = gz(i, j, k)
        END DO
      END DO
      IF (grid_type .LT. 3) CALL FILL_4CORNERS(gz2, 1, bd, npx, npy, &
&                                        sw_corner, se_corner, ne_corner&
&                                        , nw_corner)
      DO j=js1,je1
        DO i=is1,ie2
          IF (xfx(i, j) .GT. 0.) THEN
            fx(i, j) = gz2(i-1, j)
          ELSE
            fx(i, j) = gz2(i, j)
          END IF
          fx(i, j) = xfx(i, j)*fx(i, j)
        END DO
      END DO
      IF (grid_type .LT. 3) CALL FILL_4CORNERS(gz2, 2, bd, npx, npy, &
&                                        sw_corner, se_corner, ne_corner&
&                                        , nw_corner)
      DO j=js1,je2
        DO i=is1,ie1
          IF (yfx(i, j) .GT. 0.) THEN
            fy(i, j) = gz2(i, j-1)
          ELSE
            fy(i, j) = gz2(i, j)
          END IF
          fy(i, j) = yfx(i, j)*fy(i, j)
        END DO
      END DO
      DO j=js1,je1
        DO i=is1,ie1
          gz(i, j, k) = (gz2(i, j)*area(i, j)+(fx(i, j)-fx(i+1, j))+(fy(&
&           i, j)-fy(i, j+1)))/(area(i, j)+(xfx(i, j)-xfx(i+1, j))+(yfx(&
&           i, j)-yfx(i, j+1)))
        END DO
      END DO
    END DO
! Enforce monotonicity of height to prevent blowup
!$OMP parallel do default(none) shared(is1,ie1,js1,je1,ws,zs,gz,rdt,km)
    DO j=js1,je1
      DO i=is1,ie1
        ws(i, j) = (zs(i, j)-gz(i, j, km+1))*rdt
      END DO
      DO k=km,1,-1
        DO i=is1,ie1
          IF (gz(i, j, k) .LT. gz(i, j, k+1) + dz_min) THEN
            gz(i, j, k) = gz(i, j, k+1) + dz_min
          ELSE
            gz(i, j, k) = gz(i, j, k)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE UPDATE_DZ_C
!  Differentiation of update_dz_d in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: xfx ws yfx zh crx cry
!   with respect to varying inputs: xfx ws yfx zh crx cry
  SUBROUTINE UPDATE_DZ_D_FWD(ndif, damp, hord, is, ie, js, je, km, ng, &
&   npx, npy, area, rarea, dp0, zs, zh, crx, cry, xfx, yfx, delz, ws, &
&   rdt, gridstruct, bd, hord_pert)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km, npx, npy
    INTEGER, INTENT(IN) :: hord, hord_pert
    REAL, INTENT(IN) :: rdt
    REAL, INTENT(IN) :: dp0(km)
    REAL, INTENT(IN) :: area(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(IN) :: rarea(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: damp(km+1)
    INTEGER, INTENT(INOUT) :: ndif(km+1)
    REAL, INTENT(IN) :: zs(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: zh(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(OUT) :: delz(is-ng:ie+ng, js-ng:je+ng, km)
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km), INTENT(INOUT) :: crx, xfx
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km), INTENT(INOUT) :: cry, yfx
    REAL :: ws(is:ie, js:je)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
!-----------------------------------------------------
! Local array:
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km+1) :: crx_adv, xfx_adv
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km+1) :: cry_adv, yfx_adv
    REAL, DIMENSION(is:ie+1, js:je) :: fx
    REAL, DIMENSION(is:ie, js:je+1) :: fy
    REAL, DIMENSION(is-ng:ie+ng+1, js-ng:je+ng) :: fx2
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng+1) :: fy2
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng) :: wk2, z2
    REAL :: ra_x(is:ie, js-ng:je+ng)
    REAL :: ra_y(is-ng:ie+ng, js:je)
!--------------------------------------------------------------------
    INTEGER :: i, j, k, isd, ied, jsd, jed
    LOGICAL :: uniform_grid
    INTRINSIC MAX
    INTEGER :: arg1

    crx_adv = 0.0
    xfx_adv = 0.0
    cry_adv = 0.0
    yfx_adv = 0.0
    fx = 0.0
    fy = 0.0
    fx2 = 0.0
    fy2 = 0.0
    wk2 = 0.0
    z2 = 0.0
    ra_x = 0.0
    ra_y = 0.0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    arg1 = 0

    uniform_grid = .false.
    CALL PUSHREALARRAY(damp(km+1))
    damp(km+1) = damp(km)
    CALL PUSHINTEGER(ndif(km+1))
    ndif(km+1) = ndif(km)
    isd = is - ng
    ied = ie + ng
    jsd = js - ng
    jed = je + ng
!$OMP parallel do default(none) shared(jsd,jed,crx,xfx,crx_adv,xfx_adv,is,ie,isd,ied, &
!$OMP                                  km,dp0,uniform_grid,js,je,cry,yfx,cry_adv,yfx_adv)
    DO j=jsd,jed
      arg1 = ie + 1
      CALL EDGE_PROFILE_FWD(crx, xfx, crx_adv, xfx_adv, is, arg1, jsd, &
&                     jed, j, km, dp0, uniform_grid, 0)
      IF (j .LE. je + 1 .AND. j .GE. js) THEN
        arg1 = je + 1
        CALL EDGE_PROFILE_FWD(cry, yfx, cry_adv, yfx_adv, isd, ied, js, &
&                       arg1, j, km, dp0, uniform_grid, 0)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,km,area,xfx_adv,yfx_adv, &
!$OMP                                  damp,zh,crx_adv,cry_adv,npx,npy,hord,gridstruct,bd,  &
!$OMP                                  ndif,rarea) &
!$OMP                          private(z2, fx2, fy2, ra_x, ra_y, fx, fy,wk2)
    DO k=1,km+1
      DO j=jsd,jed
        DO i=is,ie
          CALL PUSHREALARRAY(ra_x(i, j))
          ra_x(i, j) = area(i, j) + (xfx_adv(i, j, k)-xfx_adv(i+1, j, k)&
&           )
        END DO
      END DO
      DO j=js,je
        DO i=isd,ied
          CALL PUSHREALARRAY(ra_y(i, j))
          ra_y(i, j) = area(i, j) + (yfx_adv(i, j, k)-yfx_adv(i, j+1, k)&
&           )
        END DO
      END DO
      IF (damp(k) .GT. 1.e-5) THEN
        DO j=jsd,jed
          DO i=isd,ied
            CALL PUSHREALARRAY(z2(i, j))
            z2(i, j) = zh(i, j, k)
          END DO
        END DO
        IF (hord .EQ. hord_pert) THEN
          CALL FV_TP_2D_FWD(z2, crx_adv(is:ie+1, jsd:jed, k), cry_adv&
&                        (isd:ied, js:je+1, k), npx, npy, hord, fx, fy, &
&                        xfx_adv(is:ie+1, jsd:jed, k), yfx_adv(isd:ied, &
&                        js:je+1, k), gridstruct, bd, ra_x, ra_y)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(fy, (ie-is+1)*(je-js+2))
          CALL PUSHREALARRAY(fx, (ie-is+2)*(je-js+1))
          CALL PUSHREALARRAY(z2, (ie+2*ng-is+1)*(je+2*ng-js+1))
          CALL FV_TP_2D(z2, crx_adv(is:ie+1, jsd:jed, k), cry_adv(isd:&
&                 ied, js:je+1, k), npx, npy, hord_pert, fx, fy, xfx_adv&
&                 (is:ie+1, jsd:jed, k), yfx_adv(isd:ied, js:je+1, k), &
&                 gridstruct, bd, ra_x, ra_y)
          CALL PUSHCONTROL(1,1)
        END IF
        CALL DEL6_VT_FLUX(ndif(k), npx, npy, damp(k), z2, wk2, fx2, fy2&
&                   , gridstruct, bd)
        DO j=js,je
          DO i=is,ie
            zh(i, j, k) = (z2(i, j)*area(i, j)+(fx(i, j)-fx(i+1, j))+(fy&
&             (i, j)-fy(i, j+1)))/(ra_x(i, j)+ra_y(i, j)-area(i, j)) + (&
&             fx2(i, j)-fx2(i+1, j)+(fy2(i, j)-fy2(i, j+1)))*rarea(i, j)
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        IF (hord .EQ. hord_pert) THEN
          CALL FV_TP_2D_FWD(zh(isd:ied, jsd:jed, k), crx_adv(is:ie+1&
&                        , jsd:jed, k), cry_adv(isd:ied, js:je+1, k), &
&                        npx, npy, hord, fx, fy, xfx_adv(is:ie+1, jsd:&
&                        jed, k), yfx_adv(isd:ied, js:je+1, k), &
&                        gridstruct, bd, ra_x, ra_y)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHREALARRAY(fy, (ie-is+1)*(je-js+2))
          CALL PUSHREALARRAY(fx, (ie-is+2)*(je-js+1))
          CALL PUSHREALARRAY(zh(isd:ied, jsd:jed, k), (ied-isd+1)*(jed-&
&                       jsd+1))
          CALL FV_TP_2D(zh(isd:ied, jsd:jed, k), crx_adv(is:ie+1, jsd:&
&                 jed, k), cry_adv(isd:ied, js:je+1, k), npx, npy, &
&                 hord_pert, fx, fy, xfx_adv(is:ie+1, jsd:jed, k), &
&                 yfx_adv(isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                 ra_y)
          CALL PUSHCONTROL(1,0)
        END IF
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(zh(i, j, k))
            zh(i, j, k) = (zh(i, j, k)*area(i, j)+(fx(i, j)-fx(i+1, j))+&
&             (fy(i, j)-fy(i, j+1)))/(ra_x(i, j)+ra_y(i, j)-area(i, j))
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
!          zh(i,j,k) = rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))   &
!                    + zh(i,j,k)*(3.-rarea(i,j)*(ra_x(i,j) + ra_y(i,j)))
!$OMP parallel do default(none) shared(is,ie,js,je,km,ws,zs,zh,rdt)
    DO j=js,je
      DO i=is,ie
        CALL PUSHREALARRAY(ws(i, j))
        ws(i, j) = (zs(i, j)-zh(i, j, km+1))*rdt
      END DO
      DO k=km,1,-1
        DO i=is,ie
          IF (zh(i, j, k) .LT. zh(i, j, k+1) + dz_min) THEN
            CALL PUSHREALARRAY(zh(i, j, k))
            zh(i, j, k) = zh(i, j, k+1) + dz_min
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHREALARRAY(zh(i, j, k))
            zh(i, j, k) = zh(i, j, k)
            CALL PUSHCONTROL(1,1)
          END IF
        END DO
      END DO
    END DO
    CALL PUSHINTEGER(jed)
    CALL PUSHREALARRAY(fy, (ie-is+1)*(je-js+2))
    CALL PUSHREALARRAY(fx, (ie-is+2)*(je-js+1))
    CALL PUSHINTEGER(isd)
    CALL PUSHREALARRAY(ra_y, (ie+2*ng-is+1)*(je-js+1))
    CALL PUSHREALARRAY(ra_x, (ie-is+1)*(je+2*ng-js+1))
    CALL PUSHREALARRAY(cry_adv, (ie+2*ng-is+1)*(je-js+2)*(km+1))
    CALL PUSHINTEGER(ied)
    CALL PUSHREALARRAY(z2, (ie+2*ng-is+1)*(je+2*ng-js+1))
    CALL PUSHINTEGER(jsd)
    CALL PUSHREALARRAY(xfx_adv, (ie-is+2)*(je+2*ng-js+1)*(km+1))
    CALL PUSHREALARRAY(crx_adv, (ie-is+2)*(je+2*ng-js+1)*(km+1))
    CALL PUSHREALARRAY(yfx_adv, (ie+2*ng-is+1)*(je-js+2)*(km+1))
  END SUBROUTINE UPDATE_DZ_D_FWD
!  Differentiation of update_dz_d in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: xfx ws yfx zh crx cry
!   with respect to varying inputs: xfx ws yfx zh crx cry
  SUBROUTINE UPDATE_DZ_D_BWD(ndif, damp, hord, is, ie, js, je, km, ng, &
&   npx, npy, area, rarea, dp0, zs, zh, zh_ad, crx, crx_ad, cry, cry_ad&
&   , xfx, xfx_ad, yfx, yfx_ad, delz, ws, ws_ad, rdt, gridstruct, bd, &
&   hord_pert)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km, npx, npy
    INTEGER, INTENT(IN) :: hord, hord_pert
    REAL, INTENT(IN) :: rdt
    REAL, INTENT(IN) :: dp0(km)
    REAL, INTENT(IN) :: area(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(IN) :: rarea(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: damp(km+1)
    INTEGER, INTENT(INOUT) :: ndif(km+1)
    REAL, INTENT(IN) :: zs(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: zh(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(INOUT) :: zh_ad(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(OUT) :: delz(is-ng:ie+ng, js-ng:je+ng, km)
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km), INTENT(INOUT) :: crx, xfx
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km), INTENT(INOUT) :: crx_ad, &
&   xfx_ad
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km), INTENT(INOUT) :: cry, yfx
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km), INTENT(INOUT) :: cry_ad, &
&   yfx_ad
    REAL :: ws(is:ie, js:je)
    REAL :: ws_ad(is:ie, js:je)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km+1) :: crx_adv, xfx_adv
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km+1) :: crx_adv_ad, &
&   xfx_adv_ad
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km+1) :: cry_adv, yfx_adv
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km+1) :: cry_adv_ad, &
&   yfx_adv_ad
    REAL, DIMENSION(is:ie+1, js:je) :: fx
    REAL, DIMENSION(is:ie+1, js:je) :: fx_ad
    REAL, DIMENSION(is:ie, js:je+1) :: fy
    REAL, DIMENSION(is:ie, js:je+1) :: fy_ad
    REAL, DIMENSION(is-ng:ie+ng+1, js-ng:je+ng) :: fx2
    REAL, DIMENSION(is-ng:ie+ng+1, js-ng:je+ng) :: fx2_ad
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng+1) :: fy2
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng+1) :: fy2_ad
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng) :: wk2, z2
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng) :: wk2_ad, z2_ad
    REAL :: ra_x(is:ie, js-ng:je+ng)
    REAL :: ra_x_ad(is:ie, js-ng:je+ng)
    REAL :: ra_y(is-ng:ie+ng, js:je)
    REAL :: ra_y_ad(is-ng:ie+ng, js:je)
    INTEGER :: i, j, k, isd, ied, jsd, jed
    LOGICAL :: uniform_grid
    INTRINSIC MAX
    INTEGER :: arg1
    REAL :: temp
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp0
    REAL :: temp_ad2
    REAL :: temp_ad3
    INTEGER :: branch

    crx_adv = 0.0
    xfx_adv = 0.0
    cry_adv = 0.0
    yfx_adv = 0.0
    fx = 0.0
    fy = 0.0
    fx2 = 0.0
    fy2 = 0.0
    wk2 = 0.0
    z2 = 0.0
    ra_x = 0.0
    ra_y = 0.0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    arg1 = 0
    branch = 0

    CALL POPREALARRAY(yfx_adv, (ie+2*ng-is+1)*(je-js+2)*(km+1))
    CALL POPREALARRAY(crx_adv, (ie-is+2)*(je+2*ng-js+1)*(km+1))
    CALL POPREALARRAY(xfx_adv, (ie-is+2)*(je+2*ng-js+1)*(km+1))
    CALL POPINTEGER(jsd)
    CALL POPREALARRAY(z2, (ie+2*ng-is+1)*(je+2*ng-js+1))
    CALL POPINTEGER(ied)
    CALL POPREALARRAY(cry_adv, (ie+2*ng-is+1)*(je-js+2)*(km+1))
    CALL POPREALARRAY(ra_x, (ie-is+1)*(je+2*ng-js+1))
    CALL POPREALARRAY(ra_y, (ie+2*ng-is+1)*(je-js+1))
    CALL POPINTEGER(isd)
    CALL POPREALARRAY(fx, (ie-is+2)*(je-js+1))
    CALL POPREALARRAY(fy, (ie-is+1)*(je-js+2))
    CALL POPINTEGER(jed)
    DO j=je,js,-1
      DO k=1,km,1
        DO i=ie,is,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(zh(i, j, k))
            zh_ad(i, j, k+1) = zh_ad(i, j, k+1) + zh_ad(i, j, k)
            zh_ad(i, j, k) = 0.0
          ELSE
            CALL POPREALARRAY(zh(i, j, k))
          END IF
        END DO
      END DO
      DO i=ie,is,-1
        CALL POPREALARRAY(ws(i, j))
        zh_ad(i, j, km+1) = zh_ad(i, j, km+1) - rdt*ws_ad(i, j)
        ws_ad(i, j) = 0.0
      END DO
    END DO
    yfx_adv_ad = 0.0
    crx_adv_ad = 0.0
    fy2_ad = 0.0
    wk2_ad = 0.0
    xfx_adv_ad = 0.0
    z2_ad = 0.0
    cry_adv_ad = 0.0
    ra_x_ad = 0.0
    ra_y_ad = 0.0
    fx_ad = 0.0
    fy_ad = 0.0
    fx2_ad = 0.0
    DO k=km+1,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(zh(i, j, k))
            temp0 = ra_x(i, j) - area(i, j) + ra_y(i, j)
            temp_ad2 = zh_ad(i, j, k)/temp0
            temp_ad3 = -((area(i, j)*zh(i, j, k)+fx(i, j)-fx(i+1, j)+fy(&
&             i, j)-fy(i, j+1))*temp_ad2/temp0)
            fx_ad(i, j) = fx_ad(i, j) + temp_ad2
            fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad2
            fy_ad(i, j) = fy_ad(i, j) + temp_ad2
            fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad2
            ra_x_ad(i, j) = ra_x_ad(i, j) + temp_ad3
            ra_y_ad(i, j) = ra_y_ad(i, j) + temp_ad3
            zh_ad(i, j, k) = area(i, j)*temp_ad2
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(zh(isd:ied, jsd:jed, k), (ied-isd+1)*(jed-&
&                      jsd+1))
          CALL POPREALARRAY(fx, (ie-is+2)*(je-js+1))
          CALL POPREALARRAY(fy, (ie-is+1)*(je-js+2))
          CALL FV_TP_2D_ADM(zh(isd:ied, jsd:jed, k), zh_ad(isd:ied, jsd:&
&                     jed, k), crx_adv(is:ie+1, jsd:jed, k), crx_adv_ad(&
&                     is:ie+1, jsd:jed, k), cry_adv(isd:ied, js:je+1, k)&
&                     , cry_adv_ad(isd:ied, js:je+1, k), npx, npy, &
&                     hord_pert, fx, fx_ad, fy, fy_ad, xfx_adv(is:ie+1, &
&                     jsd:jed, k), xfx_adv_ad(is:ie+1, jsd:jed, k), &
&                     yfx_adv(isd:ied, js:je+1, k), yfx_adv_ad(isd:ied, &
&                     js:je+1, k), gridstruct, bd, ra_x, ra_x_ad, ra_y, &
&                     ra_y_ad)
        ELSE
          CALL FV_TP_2D_BWD(zh(isd:ied, jsd:jed, k), zh_ad(isd:ied, &
&                        jsd:jed, k), crx_adv(is:ie+1, jsd:jed, k), &
&                        crx_adv_ad(is:ie+1, jsd:jed, k), cry_adv(isd:&
&                        ied, js:je+1, k), cry_adv_ad(isd:ied, js:je+1, &
&                        k), npx, npy, hord, fx, fx_ad, fy, fy_ad, &
&                        xfx_adv(is:ie+1, jsd:jed, k), xfx_adv_ad(is:ie+&
&                        1, jsd:jed, k), yfx_adv(isd:ied, js:je+1, k), &
&                        yfx_adv_ad(isd:ied, js:je+1, k), gridstruct, bd&
&                        , ra_x, ra_x_ad, ra_y, ra_y_ad)
        END IF
      ELSE
        DO j=je,js,-1
          DO i=ie,is,-1
            temp = ra_x(i, j) - area(i, j) + ra_y(i, j)
            temp_ad = zh_ad(i, j, k)/temp
            temp_ad0 = -((area(i, j)*z2(i, j)+fx(i, j)-fx(i+1, j)+fy(i, &
&             j)-fy(i, j+1))*temp_ad/temp)
            temp_ad1 = rarea(i, j)*zh_ad(i, j, k)
            z2_ad(i, j) = z2_ad(i, j) + area(i, j)*temp_ad
            fx_ad(i, j) = fx_ad(i, j) + temp_ad
            fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad
            fy_ad(i, j) = fy_ad(i, j) + temp_ad
            fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad
            ra_x_ad(i, j) = ra_x_ad(i, j) + temp_ad0
            ra_y_ad(i, j) = ra_y_ad(i, j) + temp_ad0
            fx2_ad(i, j) = fx2_ad(i, j) + temp_ad1
            fx2_ad(i+1, j) = fx2_ad(i+1, j) - temp_ad1
            fy2_ad(i, j) = fy2_ad(i, j) + temp_ad1
            fy2_ad(i, j+1) = fy2_ad(i, j+1) - temp_ad1
            zh_ad(i, j, k) = 0.0
          END DO
        END DO
        CALL DEL6_VT_FLUX_ADM(ndif(k), npx, npy, damp(k), z2, z2_ad, wk2&
&                       , wk2_ad, fx2, fx2_ad, fy2, fy2_ad, gridstruct, &
&                       bd)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL FV_TP_2D_BWD(z2, z2_ad, crx_adv(is:ie+1, jsd:jed, k), &
&                        crx_adv_ad(is:ie+1, jsd:jed, k), cry_adv(isd:&
&                        ied, js:je+1, k), cry_adv_ad(isd:ied, js:je+1, &
&                        k), npx, npy, hord, fx, fx_ad, fy, fy_ad, &
&                        xfx_adv(is:ie+1, jsd:jed, k), xfx_adv_ad(is:ie+&
&                        1, jsd:jed, k), yfx_adv(isd:ied, js:je+1, k), &
&                        yfx_adv_ad(isd:ied, js:je+1, k), gridstruct, bd&
&                        , ra_x, ra_x_ad, ra_y, ra_y_ad)
        ELSE
          CALL POPREALARRAY(z2, (ie+2*ng-is+1)*(je+2*ng-js+1))
          CALL POPREALARRAY(fx, (ie-is+2)*(je-js+1))
          CALL POPREALARRAY(fy, (ie-is+1)*(je-js+2))
          CALL FV_TP_2D_ADM(z2, z2_ad, crx_adv(is:ie+1, jsd:jed, k), &
&                     crx_adv_ad(is:ie+1, jsd:jed, k), cry_adv(isd:ied, &
&                     js:je+1, k), cry_adv_ad(isd:ied, js:je+1, k), npx&
&                     , npy, hord_pert, fx, fx_ad, fy, fy_ad, xfx_adv(is&
&                     :ie+1, jsd:jed, k), xfx_adv_ad(is:ie+1, jsd:jed, k&
&                     ), yfx_adv(isd:ied, js:je+1, k), yfx_adv_ad(isd:&
&                     ied, js:je+1, k), gridstruct, bd, ra_x, ra_x_ad, &
&                     ra_y, ra_y_ad)
        END IF
        DO j=jed,jsd,-1
          DO i=ied,isd,-1
            CALL POPREALARRAY(z2(i, j))
            zh_ad(i, j, k) = zh_ad(i, j, k) + z2_ad(i, j)
            z2_ad(i, j) = 0.0
          END DO
        END DO
      END IF
      DO j=je,js,-1
        DO i=ied,isd,-1
          CALL POPREALARRAY(ra_y(i, j))
          yfx_adv_ad(i, j, k) = yfx_adv_ad(i, j, k) + ra_y_ad(i, j)
          yfx_adv_ad(i, j+1, k) = yfx_adv_ad(i, j+1, k) - ra_y_ad(i, j)
          ra_y_ad(i, j) = 0.0
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(ra_x(i, j))
          xfx_adv_ad(i, j, k) = xfx_adv_ad(i, j, k) + ra_x_ad(i, j)
          xfx_adv_ad(i+1, j, k) = xfx_adv_ad(i+1, j, k) - ra_x_ad(i, j)
          ra_x_ad(i, j) = 0.0
        END DO
      END DO
    END DO
    DO j=jed,jsd,-1
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) CALL EDGE_PROFILE_BWD(cry, cry_ad, yfx, yfx_ad&
&                                        , cry_adv, cry_adv_ad, yfx_adv&
&                                        , yfx_adv_ad, isd, ied, js, &
&                                        arg1, j, km, dp0, uniform_grid&
&                                        , 0)
      arg1 = ie + 1
      CALL EDGE_PROFILE_BWD(crx, crx_ad, xfx, xfx_ad, crx_adv, &
&                     crx_adv_ad, xfx_adv, xfx_adv_ad, is, arg1, jsd, &
&                     jed, j, km, dp0, uniform_grid, 0)
    END DO
    CALL POPINTEGER(ndif(km+1))
    CALL POPREALARRAY(damp(km+1))
  END SUBROUTINE UPDATE_DZ_D_BWD
  SUBROUTINE UPDATE_DZ_D(ndif, damp, hord, is, ie, js, je, km, ng, npx, &
&   npy, area, rarea, dp0, zs, zh, crx, cry, xfx, yfx, delz, ws, rdt, &
&   gridstruct, bd, hord_pert)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km, npx, npy
    INTEGER, INTENT(IN) :: hord, hord_pert
    REAL, INTENT(IN) :: rdt
    REAL, INTENT(IN) :: dp0(km)
    REAL, INTENT(IN) :: area(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(IN) :: rarea(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: damp(km+1)
    INTEGER, INTENT(INOUT) :: ndif(km+1)
    REAL, INTENT(IN) :: zs(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: zh(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(OUT) :: delz(is-ng:ie+ng, js-ng:je+ng, km)
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km), INTENT(INOUT) :: crx, xfx
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km), INTENT(INOUT) :: cry, yfx
    REAL, INTENT(OUT) :: ws(is:ie, js:je)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
!-----------------------------------------------------
! Local array:
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km+1) :: crx_adv, xfx_adv
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km+1) :: cry_adv, yfx_adv
    REAL, DIMENSION(is:ie+1, js:je) :: fx
    REAL, DIMENSION(is:ie, js:je+1) :: fy
    REAL, DIMENSION(is-ng:ie+ng+1, js-ng:je+ng) :: fx2
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng+1) :: fy2
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng) :: wk2, z2
    REAL :: ra_x(is:ie, js-ng:je+ng)
    REAL :: ra_y(is-ng:ie+ng, js:je)
!--------------------------------------------------------------------
    INTEGER :: i, j, k, isd, ied, jsd, jed
    LOGICAL :: uniform_grid
    INTRINSIC MAX
    INTEGER :: arg1
    uniform_grid = .false.
    damp(km+1) = damp(km)
    ndif(km+1) = ndif(km)
    isd = is - ng
    ied = ie + ng
    jsd = js - ng
    jed = je + ng
!$OMP parallel do default(none) shared(jsd,jed,crx,xfx,crx_adv,xfx_adv,is,ie,isd,ied, &
!$OMP                                  km,dp0,uniform_grid,js,je,cry,yfx,cry_adv,yfx_adv)
    DO j=jsd,jed
      arg1 = ie + 1
      CALL EDGE_PROFILE(crx, xfx, crx_adv, xfx_adv, is, arg1, jsd, jed, &
&                 j, km, dp0, uniform_grid, 0)
      IF (j .LE. je + 1 .AND. j .GE. js) THEN
        arg1 = je + 1
        CALL EDGE_PROFILE(cry, yfx, cry_adv, yfx_adv, isd, ied, js, arg1&
&                   , j, km, dp0, uniform_grid, 0)
      END IF
    END DO
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,km,area,xfx_adv,yfx_adv, &
!$OMP                                  damp,zh,crx_adv,cry_adv,npx,npy,hord,gridstruct,bd,  &
!$OMP                                  ndif,rarea) &
!$OMP                          private(z2, fx2, fy2, ra_x, ra_y, fx, fy,wk2)
    DO k=1,km+1
      DO j=jsd,jed
        DO i=is,ie
          ra_x(i, j) = area(i, j) + (xfx_adv(i, j, k)-xfx_adv(i+1, j, k)&
&           )
        END DO
      END DO
      DO j=js,je
        DO i=isd,ied
          ra_y(i, j) = area(i, j) + (yfx_adv(i, j, k)-yfx_adv(i, j+1, k)&
&           )
        END DO
      END DO
      IF (damp(k) .GT. 1.e-5) THEN
        DO j=jsd,jed
          DO i=isd,ied
            z2(i, j) = zh(i, j, k)
          END DO
        END DO
        IF (hord .EQ. hord_pert) THEN
          CALL FV_TP_2D(z2, crx_adv(is:ie+1, jsd:jed, k), cry_adv(isd&
&                    :ied, js:je+1, k), npx, npy, hord, fx, fy, xfx_adv(&
&                    is:ie+1, jsd:jed, k), yfx_adv(isd:ied, js:je+1, k)&
&                    , gridstruct, bd, ra_x, ra_y)
        ELSE
          CALL FV_TP_2D(z2, crx_adv(is:ie+1, jsd:jed, k), cry_adv(isd:&
&                 ied, js:je+1, k), npx, npy, hord_pert, fx, fy, xfx_adv&
&                 (is:ie+1, jsd:jed, k), yfx_adv(isd:ied, js:je+1, k), &
&                 gridstruct, bd, ra_x, ra_y)
        END IF
        CALL DEL6_VT_FLUX(ndif(k), npx, npy, damp(k), z2, wk2, fx2, fy2&
&                   , gridstruct, bd)
        DO j=js,je
          DO i=is,ie
            zh(i, j, k) = (z2(i, j)*area(i, j)+(fx(i, j)-fx(i+1, j))+(fy&
&             (i, j)-fy(i, j+1)))/(ra_x(i, j)+ra_y(i, j)-area(i, j)) + (&
&             fx2(i, j)-fx2(i+1, j)+(fy2(i, j)-fy2(i, j+1)))*rarea(i, j)
          END DO
        END DO
      ELSE
        IF (hord .EQ. hord_pert) THEN
          CALL FV_TP_2D(zh(isd:ied, jsd:jed, k), crx_adv(is:ie+1, jsd&
&                    :jed, k), cry_adv(isd:ied, js:je+1, k), npx, npy, &
&                    hord, fx, fy, xfx_adv(is:ie+1, jsd:jed, k), yfx_adv&
&                    (isd:ied, js:je+1, k), gridstruct, bd, ra_x, ra_y)
        ELSE
          CALL FV_TP_2D(zh(isd:ied, jsd:jed, k), crx_adv(is:ie+1, jsd:&
&                 jed, k), cry_adv(isd:ied, js:je+1, k), npx, npy, &
&                 hord_pert, fx, fy, xfx_adv(is:ie+1, jsd:jed, k), &
&                 yfx_adv(isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                 ra_y)
        END IF
        DO j=js,je
          DO i=is,ie
            zh(i, j, k) = (zh(i, j, k)*area(i, j)+(fx(i, j)-fx(i+1, j))+&
&             (fy(i, j)-fy(i, j+1)))/(ra_x(i, j)+ra_y(i, j)-area(i, j))
          END DO
        END DO
      END IF
    END DO
!          zh(i,j,k) = rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))   &
!                    + zh(i,j,k)*(3.-rarea(i,j)*(ra_x(i,j) + ra_y(i,j)))
!$OMP parallel do default(none) shared(is,ie,js,je,km,ws,zs,zh,rdt)
    DO j=js,je
      DO i=is,ie
        ws(i, j) = (zs(i, j)-zh(i, j, km+1))*rdt
      END DO
      DO k=km,1,-1
        DO i=is,ie
          IF (zh(i, j, k) .LT. zh(i, j, k+1) + dz_min) THEN
            zh(i, j, k) = zh(i, j, k+1) + dz_min
          ELSE
            zh(i, j, k) = zh(i, j, k)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE UPDATE_DZ_D
!  Differentiation of riem_solver_c in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: ws gz delp w3 pef pt
!   with respect to varying inputs: ws gz delp w3 pef pt
  SUBROUTINE RIEM_SOLVER_C_FWD(ms, dt, is, ie, js, je, km, ng, akap, &
&   cappa, cp, ptop, hs, w3, pt, q_con, delp, gz, pef, ws, p_fac, a_imp&
&   , scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km
    INTEGER, INTENT(IN) :: ms
    REAL, INTENT(IN) :: dt, akap, cp, ptop, p_fac, a_imp, scale_m
    REAL, INTENT(IN) :: ws(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: pt, &
&   delp
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: q_con, &
&   cappa
    REAL, INTENT(IN) :: hs(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: w3
! OUTPUT PARAMETERS
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1), INTENT(INOUT) :: gz
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1) :: pef
! Local:
    REAL, DIMENSION(is-1:ie+1, km) :: dm, dz2, w2, pm2, gm2, cp2
    REAL, DIMENSION(is-1:ie+1, km+1) :: pem, pe2, peg
    REAL :: gama, rgrav
    INTEGER :: i, j, k
    INTEGER :: is1, ie1
    INTRINSIC LOG

    dm = 0.0
    dz2 = 0.0
    w2 = 0.0
    pm2 = 0.0
    gm2 = 0.0
    cp2 = 0.0
    pem = 0.0
    pe2 = 0.0
    peg = 0.0
    gama = 0.0
    rgrav = 0.0
    is1 = 0
    ie1 = 0

    gama = 1./(1.-akap)
    rgrav = 1./grav
    is1 = is - 1
    ie1 = ie + 1
!$OMP parallel do default(none) shared(js,je,is1,ie1,km,delp,pef,ptop,gz,rgrav,w3,pt, &
!$OMP                                  a_imp,dt,gama,akap,ws,p_fac,scale_m,ms,hs,q_con,cappa) &
!$OMP                          private(cp2,gm2, dm, dz2, w2, pm2, pe2, pem, peg)
    DO j=js-1,je+1
      DO k=1,km
        DO i=is1,ie1
          CALL PUSHREALARRAY(dm(i, k))
          dm(i, k) = delp(i, j, k)
        END DO
      END DO
      DO i=is1,ie1
! full pressure at top
        CALL PUSHREALARRAY(pef(i, j, 1))
        pef(i, j, 1) = ptop
        CALL PUSHREALARRAY(pem(i, 1))
        pem(i, 1) = ptop
      END DO
      DO k=2,km+1
        DO i=is1,ie1
          CALL PUSHREALARRAY(pem(i, k))
          pem(i, k) = pem(i, k-1) + dm(i, k-1)
        END DO
      END DO
      DO k=1,km
        DO i=is1,ie1
          CALL PUSHREALARRAY(dz2(i, k))
          dz2(i, k) = gz(i, j, k+1) - gz(i, j, k)
          CALL PUSHREALARRAY(pm2(i, k))
          pm2(i, k) = dm(i, k)/LOG(pem(i, k+1)/pem(i, k))
          CALL PUSHREALARRAY(dm(i, k))
          dm(i, k) = dm(i, k)*rgrav
          CALL PUSHREALARRAY(w2(i, k))
          w2(i, k) = w3(i, j, k)
        END DO
      END DO
      IF (a_imp .LT. -0.01) THEN
        CALL SIM3P0_SOLVER_FWD(dt, is1, ie1, km, rdgas, gama, akap, pe2&
&                        , dm, pem, w2, dz2, pt(is1:ie1, j, 1:km), ws(&
&                        is1:ie1, j), p_fac, scale_m)
        CALL PUSHCONTROL(2,2)
      ELSE IF (a_imp .LE. 0.5) THEN
        CALL RIM_2D_FWD(ms, dt, is1, ie1, km, rdgas, gama, gm2, pe2, dm&
&                 , pm2, w2, dz2, pt(is1:ie1, j, 1:km), ws(is1:ie1, j), &
&                 .true.)
        CALL PUSHCONTROL(2,1)
      ELSE
        CALL SIM1_SOLVER_FWD(dt, is1, ie1, km, rdgas, gama, gm2, cp2, &
&                      akap, pe2, dm, pm2, pem, w2, dz2, pt(is1:ie1, j, &
&                      1:km), ws(is1:ie1, j), p_fac)
        CALL PUSHCONTROL(2,0)
      END IF
      DO k=2,km+1
        DO i=is1,ie1
! add hydrostatic full-component
          CALL PUSHREALARRAY(pef(i, j, k))
          pef(i, j, k) = pe2(i, k) + pem(i, k)
        END DO
      END DO
! Compute Height * grav (for p-gradient computation)
      DO i=is1,ie1
        CALL PUSHREALARRAY(gz(i, j, km+1))
        gz(i, j, km+1) = hs(i, j)
      END DO
      DO k=km,1,-1
        DO i=is1,ie1
          CALL PUSHREALARRAY(gz(i, j, k))
          gz(i, j, k) = gz(i, j, k+1) - dz2(i, k)*grav
        END DO
      END DO
    END DO
    CALL PUSHINTEGER(ie1)
    CALL PUSHREALARRAY(pem, (ie-is+3)*(km+1))
    CALL PUSHREALARRAY(gama)
    CALL PUSHREALARRAY(pm2, (ie-is+3)*km)
    CALL PUSHREALARRAY(rgrav)
    CALL PUSHREALARRAY(w2, (ie-is+3)*km)
    CALL PUSHREALARRAY(dz2, (ie-is+3)*km)
    CALL PUSHINTEGER(is1)
    CALL PUSHREALARRAY(pe2, (ie-is+3)*(km+1))
    CALL PUSHREALARRAY(dm, (ie-is+3)*km)
  END SUBROUTINE RIEM_SOLVER_C_FWD
!  Differentiation of riem_solver_c in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
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
!   gradient     of useful results: ws gz delp w3 pef pt
!   with respect to varying inputs: ws gz delp w3 pef pt
  SUBROUTINE RIEM_SOLVER_C_BWD(ms, dt, is, ie, js, je, km, ng, akap, &
&   cappa, cp, ptop, hs, w3, w3_ad, pt, pt_ad, q_con, delp, delp_ad, gz&
&   , gz_ad, pef, pef_ad, ws, ws_ad, p_fac, a_imp, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km
    INTEGER, INTENT(IN) :: ms
    REAL, INTENT(IN) :: dt, akap, cp, ptop, p_fac, a_imp, scale_m
    REAL, INTENT(IN) :: ws(is-ng:ie+ng, js-ng:je+ng)
    REAL :: ws_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: pt, &
&   delp
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km) :: pt_ad, delp_ad
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: q_con, &
&   cappa
    REAL, INTENT(IN) :: hs(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: w3
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km) :: w3_ad
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1), INTENT(INOUT) :: gz
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1), INTENT(INOUT) :: &
&   gz_ad
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1) :: pef
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1) :: pef_ad
    REAL, DIMENSION(is-1:ie+1, km) :: dm, dz2, w2, pm2, gm2, cp2
    REAL, DIMENSION(is-1:ie+1, km) :: dm_ad, dz2_ad, w2_ad, pm2_ad
    REAL, DIMENSION(is-1:ie+1, km+1) :: pem, pe2, peg
    REAL, DIMENSION(is-1:ie+1, km+1) :: pem_ad, pe2_ad
    REAL :: gama, rgrav
    INTEGER :: i, j, k
    INTEGER :: is1, ie1
    INTRINSIC LOG
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp_ad
    INTEGER :: branch

    dm = 0.0
    dz2 = 0.0
    w2 = 0.0
    pm2 = 0.0
    gm2 = 0.0
    cp2 = 0.0
    pem = 0.0
    pe2 = 0.0
    peg = 0.0
    gama = 0.0
    rgrav = 0.0
    is1 = 0
    ie1 = 0
    branch = 0

    CALL POPREALARRAY(dm, (ie-is+3)*km)
    CALL POPREALARRAY(pe2, (ie-is+3)*(km+1))
    CALL POPINTEGER(is1)
    CALL POPREALARRAY(dz2, (ie-is+3)*km)
    CALL POPREALARRAY(w2, (ie-is+3)*km)
    CALL POPREALARRAY(rgrav)
    CALL POPREALARRAY(pm2, (ie-is+3)*km)
    CALL POPREALARRAY(gama)
    CALL POPREALARRAY(pem, (ie-is+3)*(km+1))
    CALL POPINTEGER(ie1)
    dm_ad = 0.0
    pe2_ad = 0.0
    dz2_ad = 0.0
    w2_ad = 0.0
    pm2_ad = 0.0
    pem_ad = 0.0
    DO j=je+1,js-1,-1
      DO k=1,km,1
        DO i=ie1,is1,-1
          CALL POPREALARRAY(gz(i, j, k))
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + gz_ad(i, j, k)
          dz2_ad(i, k) = dz2_ad(i, k) - grav*gz_ad(i, j, k)
          gz_ad(i, j, k) = 0.0
        END DO
      END DO
      DO i=ie1,is1,-1
        CALL POPREALARRAY(gz(i, j, km+1))
        gz_ad(i, j, km+1) = 0.0
      END DO
      DO k=km+1,2,-1
        DO i=ie1,is1,-1
          CALL POPREALARRAY(pef(i, j, k))
          pe2_ad(i, k) = pe2_ad(i, k) + pef_ad(i, j, k)
          pem_ad(i, k) = pem_ad(i, k) + pef_ad(i, j, k)
          pef_ad(i, j, k) = 0.0
        END DO
      END DO
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        CALL SIM1_SOLVER_BWD(dt, is1, ie1, km, rdgas, gama, gm2, cp2, &
&                      akap, pe2, pe2_ad, dm, dm_ad, pm2, pm2_ad, pem, &
&                      pem_ad, w2, w2_ad, dz2, dz2_ad, pt(is1:ie1, j, 1:&
&                      km), pt_ad(is1:ie1, j, 1:km), ws(is1:ie1, j), &
&                      ws_ad(is1:ie1, j), p_fac)
      ELSE IF (branch .EQ. 1) THEN
        CALL RIM_2D_BWD(ms, dt, is1, ie1, km, rdgas, gama, gm2, pe2, &
&                 pe2_ad, dm, dm_ad, pm2, pm2_ad, w2, w2_ad, dz2, dz2_ad&
&                 , pt(is1:ie1, j, 1:km), pt_ad(is1:ie1, j, 1:km), ws(&
&                 is1:ie1, j), ws_ad(is1:ie1, j), .true.)
      ELSE
        CALL SIM3P0_SOLVER_BWD(dt, is1, ie1, km, rdgas, gama, akap, pe2&
&                        , pe2_ad, dm, dm_ad, pem, pem_ad, w2, w2_ad, &
&                        dz2, dz2_ad, pt(is1:ie1, j, 1:km), pt_ad(is1:&
&                        ie1, j, 1:km), ws(is1:ie1, j), ws_ad(is1:ie1, j&
&                        ), p_fac, scale_m)
      END IF
      DO k=km,1,-1
        DO i=ie1,is1,-1
          temp1 = pem(i, k)
          temp = pem(i, k+1)/temp1
          temp0 = LOG(temp)
          CALL POPREALARRAY(w2(i, k))
          w3_ad(i, j, k) = w3_ad(i, j, k) + w2_ad(i, k)
          w2_ad(i, k) = 0.0
          CALL POPREALARRAY(dm(i, k))
          dm_ad(i, k) = pm2_ad(i, k)/temp0 + rgrav*dm_ad(i, k)
          CALL POPREALARRAY(pm2(i, k))
          temp_ad = -(dm(i, k)*pm2_ad(i, k)/(temp*temp0**2*temp1))
          pem_ad(i, k+1) = pem_ad(i, k+1) + temp_ad
          pem_ad(i, k) = pem_ad(i, k) - temp*temp_ad
          pm2_ad(i, k) = 0.0
          CALL POPREALARRAY(dz2(i, k))
          gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + dz2_ad(i, k)
          gz_ad(i, j, k) = gz_ad(i, j, k) - dz2_ad(i, k)
          dz2_ad(i, k) = 0.0
        END DO
      END DO
      DO k=km+1,2,-1
        DO i=ie1,is1,-1
          CALL POPREALARRAY(pem(i, k))
          pem_ad(i, k-1) = pem_ad(i, k-1) + pem_ad(i, k)
          dm_ad(i, k-1) = dm_ad(i, k-1) + pem_ad(i, k)
          pem_ad(i, k) = 0.0
        END DO
      END DO
      DO i=ie1,is1,-1
        CALL POPREALARRAY(pem(i, 1))
        pem_ad(i, 1) = 0.0
        CALL POPREALARRAY(pef(i, j, 1))
        pef_ad(i, j, 1) = 0.0
      END DO
      DO k=km,1,-1
        DO i=ie1,is1,-1
          CALL POPREALARRAY(dm(i, k))
          delp_ad(i, j, k) = delp_ad(i, j, k) + dm_ad(i, k)
          dm_ad(i, k) = 0.0
        END DO
      END DO
    END DO
  END SUBROUTINE RIEM_SOLVER_C_BWD
  SUBROUTINE RIEM_SOLVER_C(ms, dt, is, ie, js, je, km, ng, akap, cappa, &
&   cp, ptop, hs, w3, pt, q_con, delp, gz, pef, ws, p_fac, a_imp, &
&   scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km
    INTEGER, INTENT(IN) :: ms
    REAL, INTENT(IN) :: dt, akap, cp, ptop, p_fac, a_imp, scale_m
    REAL, INTENT(IN) :: ws(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: pt, &
&   delp
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: q_con, &
&   cappa
    REAL, INTENT(IN) :: hs(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: w3
! OUTPUT PARAMETERS
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1), INTENT(INOUT) :: gz
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1), INTENT(OUT) :: pef
! Local:
    REAL, DIMENSION(is-1:ie+1, km) :: dm, dz2, w2, pm2, gm2, cp2
    REAL, DIMENSION(is-1:ie+1, km+1) :: pem, pe2, peg
    REAL :: gama, rgrav
    INTEGER :: i, j, k
    INTEGER :: is1, ie1
    INTRINSIC LOG
    gama = 1./(1.-akap)
    rgrav = 1./grav
    is1 = is - 1
    ie1 = ie + 1
!$OMP parallel do default(none) shared(js,je,is1,ie1,km,delp,pef,ptop,gz,rgrav,w3,pt, &
!$OMP                                  a_imp,dt,gama,akap,ws,p_fac,scale_m,ms,hs,q_con,cappa) &
!$OMP                          private(cp2,gm2, dm, dz2, w2, pm2, pe2, pem, peg)
    DO j=js-1,je+1
      DO k=1,km
        DO i=is1,ie1
          dm(i, k) = delp(i, j, k)
        END DO
      END DO
      DO i=is1,ie1
! full pressure at top
        pef(i, j, 1) = ptop
        pem(i, 1) = ptop
      END DO
      DO k=2,km+1
        DO i=is1,ie1
          pem(i, k) = pem(i, k-1) + dm(i, k-1)
        END DO
      END DO
      DO k=1,km
        DO i=is1,ie1
          dz2(i, k) = gz(i, j, k+1) - gz(i, j, k)
          pm2(i, k) = dm(i, k)/LOG(pem(i, k+1)/pem(i, k))
          dm(i, k) = dm(i, k)*rgrav
          w2(i, k) = w3(i, j, k)
        END DO
      END DO
      IF (a_imp .LT. -0.01) THEN
        CALL SIM3P0_SOLVER(dt, is1, ie1, km, rdgas, gama, akap, pe2, dm&
&                    , pem, w2, dz2, pt(is1:ie1, j, 1:km), ws(is1:ie1, j&
&                    ), p_fac, scale_m)
      ELSE IF (a_imp .LE. 0.5) THEN
        CALL RIM_2D(ms, dt, is1, ie1, km, rdgas, gama, gm2, pe2, dm, pm2&
&             , w2, dz2, pt(is1:ie1, j, 1:km), ws(is1:ie1, j), .true.)
      ELSE
        CALL SIM1_SOLVER(dt, is1, ie1, km, rdgas, gama, gm2, cp2, akap, &
&                  pe2, dm, pm2, pem, w2, dz2, pt(is1:ie1, j, 1:km), ws(&
&                  is1:ie1, j), p_fac)
      END IF
      DO k=2,km+1
        DO i=is1,ie1
! add hydrostatic full-component
          pef(i, j, k) = pe2(i, k) + pem(i, k)
        END DO
      END DO
! Compute Height * grav (for p-gradient computation)
      DO i=is1,ie1
        gz(i, j, km+1) = hs(i, j)
      END DO
      DO k=km,1,-1
        DO i=is1,ie1
          gz(i, j, k) = gz(i, j, k+1) - dz2(i, k)*grav
        END DO
      END DO
    END DO
  END SUBROUTINE RIEM_SOLVER_C
!GFDL - This routine will not give absoulte reproducibility when compiled with -fast-transcendentals.
!GFDL - It is now inside of nh_core.F90 and being compiled without -fast-transcendentals.
  SUBROUTINE RIEM_SOLVER3TEST(ms, dt, is, ie, js, je, km, ng, isd, ied, &
&   jsd, jed, akap, cappa, cp, ptop, zs, q_con, w, delz, pt, delp, zh, &
&   pe, ppe, pk3, pk, peln, ws, scale_m, p_fac, a_imp, use_logp, &
&   last_call, fp_out)
    IMPLICIT NONE
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: gz: grav*height at edges
!        pe: full     hydrostatic pressure
!       ppe: non-hydrostatic pressure perturbation
!--------------------------------------------
    INTEGER, INTENT(IN) :: ms, is, ie, js, je, km, ng
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
! the BIG horizontal Lagrangian time step
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: akap, cp, ptop, p_fac, a_imp, scale_m
    REAL, INTENT(IN) :: zs(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: last_call, use_logp, fp_out
    REAL, INTENT(IN) :: ws(is:ie, js:je)
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: q_con, cappa
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: delp, pt
    REAL, DIMENSION(isd:ied, jsd:jed, km+1), INTENT(INOUT) :: zh
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: w
    REAL, INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
! ln(pe)
    REAL, INTENT(OUT) :: peln(is:ie, km+1, js:je)
    REAL, DIMENSION(isd:ied, jsd:jed, km+1), INTENT(OUT) :: ppe
    REAL, INTENT(OUT) :: delz(is-ng:ie+ng, js-ng:je+ng, km)
    REAL, INTENT(OUT) :: pk(is:ie, js:je, km+1)
    REAL, INTENT(OUT) :: pk3(isd:ied, jsd:jed, km+1)
! Local:
    REAL, DIMENSION(is:ie, km) :: dm, dz2, pm2, w2, gm2, cp2
    REAL, DIMENSION(is:ie, km+1) :: pem, pe2, peln2, peg, pelng
    REAL :: gama, rgrav, ptk, peln1
    INTEGER :: i, j, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC ABS
    REAL :: abs0
    gama = 1./(1.-akap)
    rgrav = 1./grav
    peln1 = LOG(ptop)
    ptk = EXP(akap*peln1)
!$OMP parallel do default(none) shared(is,ie,js,je,km,delp,ptop,peln1,pk3,ptk,akap,rgrav,zh,pt, &
!$OMP                                  w,a_imp,dt,gama,ws,p_fac,scale_m,ms,delz,last_call,  &
!$OMP                                  peln,pk,fp_out,ppe,use_logp,zs,pe,cappa,q_con )          &
!$OMP                          private(cp2, gm2, dm, dz2, pm2, pem, peg, pelng, pe2, peln2, w2)
    DO j=js,je
      DO k=1,km
        DO i=is,ie
          dm(i, k) = delp(i, j, k)
        END DO
      END DO
      DO i=is,ie
        pem(i, 1) = ptop
        peln2(i, 1) = peln1
        pk3(i, j, 1) = ptk
      END DO
      DO k=2,km+1
        DO i=is,ie
          pem(i, k) = pem(i, k-1) + dm(i, k-1)
          peln2(i, k) = LOG(pem(i, k))
          pk3(i, j, k) = EXP(akap*peln2(i, k))
        END DO
      END DO
      DO k=1,km
        DO i=is,ie
          pm2(i, k) = dm(i, k)/(peln2(i, k+1)-peln2(i, k))
          dm(i, k) = dm(i, k)*rgrav
          dz2(i, k) = zh(i, j, k+1) - zh(i, j, k)
          w2(i, k) = w(i, j, k)
        END DO
      END DO
      IF (a_imp .LT. -0.999) THEN
        CALL SIM3P0_SOLVER(dt, is, ie, km, rdgas, gama, akap, pe2, dm, &
&                    pem, w2, dz2, pt(is:ie, j, 1:km), ws(is:ie, j), &
&                    p_fac, scale_m)
      ELSE IF (a_imp .LT. -0.5) THEN
        IF (a_imp .GE. 0.) THEN
          abs0 = a_imp
        ELSE
          abs0 = -a_imp
        END IF
        CALL SIM3_SOLVER(dt, is, ie, km, rdgas, gama, akap, pe2, dm, pem&
&                  , w2, dz2, pt(is:ie, j, 1:km), ws(is:ie, j), abs0, &
&                  p_fac, scale_m)
      ELSE IF (a_imp .LE. 0.5) THEN
        CALL RIM_2D(ms, dt, is, ie, km, rdgas, gama, gm2, pe2, dm, pm2, &
&             w2, dz2, pt(is:ie, j, 1:km), ws(is:ie, j), .false.)
      ELSE IF (a_imp .GT. 0.999) THEN
        CALL SIM1_SOLVER(dt, is, ie, km, rdgas, gama, gm2, cp2, akap, &
&                  pe2, dm, pm2, pem, w2, dz2, pt(is:ie, j, 1:km), ws(is&
&                  :ie, j), p_fac)
      ELSE
        CALL SIM_SOLVER(dt, is, ie, km, rdgas, gama, gm2, cp2, akap, pe2&
&                 , dm, pm2, pem, w2, dz2, pt(is:ie, j, 1:km), ws(is:ie&
&                 , j), a_imp, p_fac, scale_m)
      END IF
      DO k=1,km
        DO i=is,ie
          w(i, j, k) = w2(i, k)
          delz(i, j, k) = dz2(i, k)
        END DO
      END DO
      IF (last_call) THEN
        DO k=1,km+1
          DO i=is,ie
            peln(i, k, j) = peln2(i, k)
            pk(i, j, k) = pk3(i, j, k)
            pe(i, k, j) = pem(i, k)
          END DO
        END DO
      END IF
      IF (fp_out) THEN
        DO k=1,km+1
          DO i=is,ie
            ppe(i, j, k) = pe2(i, k) + pem(i, k)
          END DO
        END DO
      ELSE
        DO k=1,km+1
          DO i=is,ie
            ppe(i, j, k) = pe2(i, k)
          END DO
        END DO
      END IF
      IF (use_logp) THEN
        DO k=2,km+1
          DO i=is,ie
            pk3(i, j, k) = peln2(i, k)
          END DO
        END DO
      END IF
      DO i=is,ie
        zh(i, j, km+1) = zs(i, j)
      END DO
      DO k=km,1,-1
        DO i=is,ie
          zh(i, j, k) = zh(i, j, k+1) - dz2(i, k)
        END DO
      END DO
    END DO
  END SUBROUTINE RIEM_SOLVER3TEST
  SUBROUTINE IMP_DIFF_W(j, is, ie, js, je, ng, km, cd, delz, ws, w, w3)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j, is, ie, js, je, km, ng
    REAL, INTENT(IN) :: cd
! delta-height (m)
    REAL, INTENT(IN) :: delz(is-ng:ie+ng, km)
! vertical vel. (m/s)
    REAL, INTENT(IN) :: w(is:ie, km)
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, INTENT(OUT) :: w3(is-ng:ie+ng, js-ng:je+ng, km)
! Local:
    REAL, DIMENSION(is:ie, km) :: c, gam, dz, wt
    REAL :: bet(is:ie)
    REAL :: a
    INTEGER :: i, k
    DO k=2,km
      DO i=is,ie
        dz(i, k) = 0.5*(delz(i, k-1)+delz(i, k))
      END DO
    END DO
    DO k=1,km-1
      DO i=is,ie
        c(i, k) = -(cd/(dz(i, k+1)*delz(i, k)))
      END DO
    END DO
! model top:
    DO i=is,ie
! bet(i) = b
      bet(i) = 1. - c(i, 1)
      wt(i, 1) = w(i, 1)/bet(i)
    END DO
! Interior:
    DO k=2,km-1
      DO i=is,ie
        gam(i, k) = c(i, k-1)/bet(i)
        a = cd/(dz(i, k)*delz(i, k))
        bet(i) = 1. + a - c(i, k) + a*gam(i, k)
        wt(i, k) = (w(i, k)+a*wt(i, k-1))/bet(i)
      END DO
    END DO
! Bottom:
    DO i=is,ie
      gam(i, km) = c(i, km-1)/bet(i)
      a = cd/(dz(i, km)*delz(i, km))
      wt(i, km) = (w(i, km)+2.*ws(i)*cd/delz(i, km)**2+a*wt(i, km-1))/(&
&       1.+a+(cd+cd)/delz(i, km)**2+a*gam(i, km))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        wt(i, k) = wt(i, k) - gam(i, k+1)*wt(i, k+1)
      END DO
    END DO
    DO k=1,km
      DO i=is,ie
        w3(i, j, k) = wt(i, k)
      END DO
    END DO
  END SUBROUTINE IMP_DIFF_W
!  Differentiation of rim_2d in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a
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
!   gradient     of useful results: ws pe2 dm2 dz2 w2 pm2 pt2
!   with respect to varying inputs: ws pe2 dm2 dz2 w2 pm2 pt2
  SUBROUTINE RIM_2D_FWD(ms, bdt, is, ie, km, rgas, gama, gm2, pe2, dm2, &
&   pm2, w2, dz2, pt2, ws, c_core)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ms, is, ie, km
    REAL, INTENT(IN) :: bdt, gama, rgas
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pm2, gm2
    LOGICAL, INTENT(IN) :: c_core
    REAL, INTENT(IN) :: pt2(is:ie, km)
    REAL, INTENT(IN) :: ws(is:ie)
! IN/OUT:
    REAL, INTENT(INOUT) :: dz2(is:ie, km)
    REAL, INTENT(INOUT) :: w2(is:ie, km)
    REAL :: pe2(is:ie, km+1)
! Local:
    REAL :: ws2(is:ie)
    REAL, DIMENSION(km+1) :: m_bot, m_top, r_bot, r_top, pe1, pbar, wbar
    REAL, DIMENSION(km) :: r_hi, r_lo, dz, wm, dm, dts
    REAL, DIMENSION(km) :: pf1, wc, cm, pp, pt1
    REAL :: dt, rdt, grg, z_frac, ptmp1, rden, pf, time_left
    REAL :: m_surf
    INTEGER :: i, k, n, ke, kt1, ktop
    INTEGER :: ks0, ks1
    INTRINSIC REAL
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC SQRT
    INTRINSIC MAX
    INTEGER :: ad_count
    INTEGER :: ad_from
    INTEGER :: ad_count0
    INTEGER :: ad_from0
    INTEGER :: ad_from1
    INTEGER :: ad_count1
    INTEGER :: ad_from2
    INTEGER :: ad_count2
    INTEGER :: ad_count3
    INTEGER :: ad_from3
    INTEGER :: ad_from4
    INTEGER :: ad_from5
    INTEGER :: ad_from6
    INTEGER :: ad_from7
    INTEGER :: ad_from8

    ws2 = 0.0
    m_bot = 0.0
    m_top = 0.0
    r_bot = 0.0
    r_top = 0.0
    pe1 = 0.0
    pbar = 0.0
    wbar = 0.0
    r_hi = 0.0
    r_lo = 0.0
    dz = 0.0
    wm = 0.0
    dm = 0.0
    dts = 0.0
    pf1 = 0.0
    wc = 0.0
    cm = 0.0
    pp = 0.0
    pt1 = 0.0
    dt = 0.0
    rdt = 0.0
    grg = 0.0
    z_frac = 0.0
    ptmp1 = 0.0
    rden = 0.0
    pf = 0.0
    time_left = 0.0
    m_surf = 0.0
    n = 0.0
    ke = 0.0
    kt1 = 0.0
    ktop = 0.0
    ks0 = 0.0
    ks1 = 0.0
    ad_count = 0.0
    ad_from = 0.0
    ad_count0 = 0.0
    ad_from0 = 0.0
    ad_from1 = 0.0
    ad_count1 = 0.0
    ad_from2 = 0.0
    ad_count2 = 0.0
    ad_count3 = 0.0
    ad_from3 = 0.0
    ad_from4 = 0.0
    ad_from5 = 0.0
    ad_from6 = 0.0
    ad_from7 = 0.0
    ad_from8 = 0.0

    grg = gama*rgas
    rdt = 1./bdt
    dt = bdt/REAL(ms)
    pbar(:) = 0.
    wbar(:) = 0.
    DO i=is,ie
      ws2(i) = 2.*ws(i)
    END DO
! end i-loop
    DO 170 i=is,ie
      DO k=1,km
        CALL PUSHREALARRAY(dz(k))
        dz(k) = dz2(i, k)
        CALL PUSHREALARRAY(dm(k))
        dm(k) = dm2(i, k)
        CALL PUSHREALARRAY(wm(k))
        wm(k) = w2(i, k)*dm(k)
        CALL PUSHREALARRAY(pt1(k))
        pt1(k) = pt2(i, k)
      END DO
      pe1(:) = 0.
      CALL PUSHREALARRAY(wbar(km+1))
      wbar(km+1) = ws(i)
      CALL PUSHINTEGER(ks0)
      ks0 = 1
      IF (ms .GT. 1 .AND. ms .LT. 8) THEN
        CALL PUSHINTEGER(k)
        ad_count = 1
! Continuity of (pbar, wbar) is maintained
        DO k=1,km
          rden = -(rgas*dm(k)/dz(k))
          CALL PUSHREALARRAY(pf1(k))
          pf1(k) = EXP(gama*LOG(rden*pt1(k)))
          CALL PUSHREALARRAY(dts(k))
          dts(k) = -(dz(k)/SQRT(grg*pf1(k)/rden))
          IF (bdt .GT. dts(k)) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER(k)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
        CALL PUSHINTEGER(ad_count)
        ks0 = km
        GOTO 222
 100    CALL PUSHCONTROL(1,1)
        CALL PUSHINTEGER(ad_count)
        ks0 = k - 1
 222    IF (ks0 .EQ. 1) THEN
          CALL PUSHCONTROL(2,0)
        ELSE
          CALL PUSHINTEGER(k)
          DO k=1,ks0
            CALL PUSHREALARRAY(cm(k))
            cm(k) = dm(k)/dts(k)
            CALL PUSHREALARRAY(wc(k))
            wc(k) = wm(k)/dts(k)
            CALL PUSHREALARRAY(pp(k))
            pp(k) = pf1(k) - pm2(i, k)
          END DO
          CALL PUSHINTEGER(k - 1)
          CALL PUSHREALARRAY(wbar(1))
          wbar(1) = (wc(1)+pp(1))/cm(1)
          DO k=2,ks0
            CALL PUSHREALARRAY(wbar(k))
            wbar(k) = (wc(k-1)+wc(k)+pp(k)-pp(k-1))/(cm(k-1)+cm(k))
            CALL PUSHREALARRAY(pbar(k))
            pbar(k) = bdt*(cm(k-1)*wbar(k)-wc(k-1)+pp(k-1))
            pe1(k) = pbar(k)
          END DO
          CALL PUSHINTEGER(k - 1)
          IF (ks0 .EQ. km) THEN
            CALL PUSHREALARRAY(pbar(km+1))
            pbar(km+1) = bdt*(cm(km)*wbar(km+1)-wc(km)+pp(km))
            IF (c_core) THEN
              DO k=1,km
                CALL PUSHREALARRAY(dz2(i, k))
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
              END DO
              CALL PUSHCONTROL(1,0)
            ELSE
              DO k=1,km
                CALL PUSHREALARRAY(dz2(i, k))
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
                CALL PUSHREALARRAY(w2(i, k))
                w2(i, k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
              END DO
              CALL PUSHCONTROL(1,1)
            END IF
            CALL PUSHREALARRAY(pe2(i, 1))
            pe2(i, 1) = 0.
            DO k=2,km+1
              CALL PUSHREALARRAY(pe2(i, k))
              pe2(i, k) = pbar(k)*rdt
            END DO
            CALL PUSHCONTROL(1,0)
            GOTO 170
          ELSE
! next i
            IF (c_core) THEN
              DO k=1,ks0-1
                CALL PUSHREALARRAY(dz2(i, k))
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
              END DO
              CALL PUSHINTEGER(k - 1)
              CALL PUSHCONTROL(1,0)
            ELSE
              DO k=1,ks0-1
                CALL PUSHREALARRAY(dz2(i, k))
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
                CALL PUSHREALARRAY(w2(i, k))
                w2(i, k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
              END DO
              CALL PUSHINTEGER(k - 1)
              CALL PUSHCONTROL(1,1)
            END IF
            CALL PUSHREALARRAY(pbar(ks0))
            pbar(ks0) = pbar(ks0)/REAL(ms)
            CALL PUSHCONTROL(2,1)
          END IF
        END IF
      ELSE
        CALL PUSHCONTROL(2,2)
      END IF
      ks1 = ks0
      DO n=1,ms
        ad_from = ks1
        CALL PUSHINTEGER(k)
        DO k=ad_from,km
          rden = -(rgas*dm(k)/dz(k))
          CALL PUSHREALARRAY(pf)
          pf = EXP(gama*LOG(rden*pt1(k)))
          CALL PUSHREALARRAY(dts(k))
          dts(k) = -(dz(k)/SQRT(grg*pf/rden))
          ptmp1 = dts(k)*(pf-pm2(i, k))
          CALL PUSHREALARRAY(r_lo(k))
          r_lo(k) = wm(k) + ptmp1
          CALL PUSHREALARRAY(r_hi(k))
          r_hi(k) = wm(k) - ptmp1
        END DO
        CALL PUSHINTEGER(ad_from)
        ad_count0 = 1
        DO k=ks1,km
          IF (dt .GT. dts(k)) THEN
            GOTO 110
          ELSE
            ad_count0 = ad_count0 + 1
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
        CALL PUSHINTEGER(ad_count0)
        ktop = km
        GOTO 333
 110    CALL PUSHCONTROL(1,1)
        CALL PUSHINTEGER(ad_count0)
        ktop = k - 1
 333    IF (ktop .GE. ks1) THEN
          ad_from0 = ks1
          DO k=ad_from0,ktop
            z_frac = dt/dts(k)
            CALL PUSHREALARRAY(r_bot(k))
            r_bot(k) = z_frac*r_lo(k)
            CALL PUSHREALARRAY(r_top(k+1))
            r_top(k+1) = z_frac*r_hi(k)
            CALL PUSHREALARRAY(m_bot(k))
            m_bot(k) = z_frac*dm(k)
            CALL PUSHREALARRAY(m_top(k+1))
            m_top(k+1) = m_bot(k)
          END DO
          CALL PUSHINTEGER(k - 1)
          CALL PUSHINTEGER(ad_from0)
          IF (ktop .EQ. km) THEN
            CALL PUSHCONTROL(1,0)
            GOTO 666
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        ad_from1 = ktop + 2
        DO k=ad_from1,km+1
          CALL PUSHREALARRAY(m_top(k))
          m_top(k) = 0.
          CALL PUSHREALARRAY(r_top(k))
          r_top(k) = 0.
        END DO
        CALL PUSHINTEGER(ad_from1)
        IF (1 .LT. ktop) THEN
          kt1 = ktop
        ELSE
          kt1 = 1
        END IF
        DO 130 ke=km+1,ktop+2,-1
          CALL PUSHREALARRAY(time_left)
          time_left = dt
          CALL PUSHINTEGER(k)
          ad_count1 = 1
          DO k=ke-1,kt1,-1
            IF (time_left .GT. dts(k)) THEN
              time_left = time_left - dts(k)
              CALL PUSHREALARRAY(m_top(ke))
              m_top(ke) = m_top(ke) + dm(k)
              CALL PUSHREALARRAY(r_top(ke))
              r_top(ke) = r_top(ke) + r_hi(k)
              CALL PUSHINTEGER(k)
              ad_count1 = ad_count1 + 1
            ELSE
              GOTO 120
            END IF
          END DO
          CALL PUSHCONTROL(1,0)
          CALL PUSHINTEGER(ad_count1)
          CALL PUSHCONTROL(1,1)
          GOTO 130
 120      CALL PUSHCONTROL(1,1)
          CALL PUSHINTEGER(ad_count1)
          z_frac = time_left/dts(k)
          CALL PUSHREALARRAY(m_top(ke))
          m_top(ke) = m_top(ke) + z_frac*dm(k)
          CALL PUSHREALARRAY(r_top(ke))
          r_top(ke) = r_top(ke) + z_frac*r_hi(k)
          CALL PUSHCONTROL(1,0)
 130    CONTINUE
        CALL PUSHINTEGER(ke + 1)
        ad_from2 = ktop + 1
        CALL PUSHINTEGER(k)
! next level
        DO k=ad_from2,km
          CALL PUSHREALARRAY(m_bot(k))
          m_bot(k) = 0.
          CALL PUSHREALARRAY(r_bot(k))
          r_bot(k) = 0.
        END DO
        CALL PUSHINTEGER(ad_from2)
        ad_from3 = ktop + 1
        DO 160 ke=ad_from3,km
          CALL PUSHREALARRAY(time_left)
          time_left = dt
          CALL PUSHINTEGER(k)
          ad_count2 = 1
          DO k=ke,km
            IF (time_left .GT. dts(k)) THEN
              time_left = time_left - dts(k)
              CALL PUSHREALARRAY(m_bot(ke))
              m_bot(ke) = m_bot(ke) + dm(k)
              CALL PUSHREALARRAY(r_bot(ke))
              r_bot(ke) = r_bot(ke) + r_lo(k)
              CALL PUSHINTEGER(k)
              ad_count2 = ad_count2 + 1
            ELSE
              GOTO 150
            END IF
          END DO
          CALL PUSHCONTROL(1,0)
          CALL PUSHINTEGER(ad_count2)
! next interface
          CALL PUSHREALARRAY(m_surf)
          m_surf = m_bot(ke)
          CALL PUSHINTEGER(k)
          ad_count3 = 1
          DO k=km,kt1,-1
            IF (time_left .GT. dts(k)) THEN
              time_left = time_left - dts(k)
              CALL PUSHREALARRAY(m_bot(ke))
              m_bot(ke) = m_bot(ke) + dm(k)
              CALL PUSHREALARRAY(r_bot(ke))
              r_bot(ke) = r_bot(ke) - r_hi(k)
              CALL PUSHINTEGER(k)
              ad_count3 = ad_count3 + 1
            ELSE
              GOTO 140
            END IF
          END DO
          CALL PUSHCONTROL(1,0)
          CALL PUSHINTEGER(ad_count3)
          CALL PUSHCONTROL(2,2)
          GOTO 160
 140      CALL PUSHCONTROL(1,1)
          CALL PUSHINTEGER(ad_count3)
          z_frac = time_left/dts(k)
          CALL PUSHREALARRAY(m_bot(ke))
          m_bot(ke) = m_bot(ke) + z_frac*dm(k)
          CALL PUSHREALARRAY(r_bot(ke))
          r_bot(ke) = r_bot(ke) - z_frac*r_hi(k) + (m_bot(ke)-m_surf)*&
&           ws2(i)
          CALL PUSHCONTROL(2,1)
          GOTO 160
 150      CALL PUSHCONTROL(1,1)
          CALL PUSHINTEGER(ad_count2)
          z_frac = time_left/dts(k)
          CALL PUSHREALARRAY(m_bot(ke))
          m_bot(ke) = m_bot(ke) + z_frac*dm(k)
          CALL PUSHREALARRAY(r_bot(ke))
          r_bot(ke) = r_bot(ke) + z_frac*r_lo(k)
          CALL PUSHCONTROL(2,0)
 160    CONTINUE
        CALL PUSHINTEGER(ad_from3)
        CALL PUSHCONTROL(1,1)
! next interface
 666    IF (ks1 .EQ. 1) THEN
          CALL PUSHREALARRAY(wbar(1))
          wbar(1) = r_bot(1)/m_bot(1)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        ad_from4 = ks1 + 1
        CALL PUSHINTEGER(k)
        DO k=ad_from4,km
          CALL PUSHREALARRAY(wbar(k))
          wbar(k) = (r_bot(k)+r_top(k))/(m_top(k)+m_bot(k))
        END DO
        CALL PUSHINTEGER(ad_from4)
        ad_from5 = ks1 + 1
! pbar here is actually dt*pbar
        DO k=ad_from5,km+1
          CALL PUSHREALARRAY(pbar(k))
          pbar(k) = m_top(k)*wbar(k) - r_top(k)
          pe1(k) = pe1(k) + pbar(k)
        END DO
        CALL PUSHINTEGER(ad_from5)
        IF (n .EQ. ms) THEN
          IF (c_core) THEN
            ad_from6 = ks1
            DO k=ad_from6,km
              CALL PUSHREALARRAY(dz2(i, k))
              dz2(i, k) = dz(k) + dt*(wbar(k+1)-wbar(k))
            END DO
            CALL PUSHINTEGER(ad_from6)
            CALL PUSHCONTROL(2,0)
          ELSE
            ad_from7 = ks1
            DO k=ad_from7,km
              CALL PUSHREALARRAY(dz2(i, k))
              dz2(i, k) = dz(k) + dt*(wbar(k+1)-wbar(k))
              CALL PUSHREALARRAY(w2(i, k))
              w2(i, k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
            END DO
            CALL PUSHINTEGER(ad_from7)
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE
          ad_from8 = ks1
          DO k=ad_from8,km
            CALL PUSHREALARRAY(dz(k))
            dz(k) = dz(k) + dt*(wbar(k+1)-wbar(k))
            CALL PUSHREALARRAY(wm(k))
            wm(k) = wm(k) + pbar(k+1) - pbar(k)
          END DO
          CALL PUSHINTEGER(ad_from8)
          CALL PUSHCONTROL(2,2)
        END IF
      END DO
      CALL PUSHREALARRAY(pe2(i, 1))
      pe2(i, 1) = 0.
      CALL PUSHINTEGER(k)
      DO k=2,km+1
        CALL PUSHREALARRAY(pe2(i, k))
        pe2(i, k) = pe1(k)*rdt
      END DO
      CALL PUSHCONTROL(1,1)
 170 CONTINUE
    CALL PUSHREALARRAY(r_hi, km)
    CALL PUSHREALARRAY(wm, km)
    CALL PUSHREALARRAY(m_surf)
    CALL PUSHREALARRAY(pbar, km + 1)
    CALL PUSHREALARRAY(wc, km)
    CALL PUSHREALARRAY(pp, km)
    CALL PUSHREALARRAY(cm, km)
    CALL PUSHREALARRAY(pt1, km)
    CALL PUSHREALARRAY(time_left)
    CALL PUSHREALARRAY(ws2, ie - is + 1)
    CALL PUSHREALARRAY(r_bot, km + 1)
    CALL PUSHREALARRAY(pf)
    CALL PUSHREALARRAY(rdt)
    CALL PUSHREALARRAY(m_bot, km + 1)
    CALL PUSHREALARRAY(r_top, km + 1)
    CALL PUSHREALARRAY(m_top, km + 1)
    CALL PUSHREALARRAY(pf1, km)
    CALL PUSHINTEGER(ks0)
    CALL PUSHREALARRAY(wbar, km + 1)
    CALL PUSHREALARRAY(r_lo, km)
    CALL PUSHREALARRAY(dz, km)
    CALL PUSHREALARRAY(grg)
    CALL PUSHREALARRAY(dt)
    CALL PUSHREALARRAY(dts, km)
    CALL PUSHREALARRAY(dm, km)
  END SUBROUTINE RIM_2D_FWD
!  Differentiation of rim_2d in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.
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
!   gradient     of useful results: ws pe2 dm2 dz2 w2 pm2 pt2
!   with respect to varying inputs: ws pe2 dm2 dz2 w2 pm2 pt2
  SUBROUTINE RIM_2D_BWD(ms, bdt, is, ie, km, rgas, gama, gm2, pe2, &
&   pe2_ad, dm2, dm2_ad, pm2, pm2_ad, w2, w2_ad, dz2, dz2_ad, pt2, &
&   pt2_ad, ws, ws_ad, c_core)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ms, is, ie, km
    REAL, INTENT(IN) :: bdt, gama, rgas
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pm2, gm2
    REAL, DIMENSION(is:ie, km) :: dm2_ad, pm2_ad
    LOGICAL, INTENT(IN) :: c_core
    REAL, INTENT(IN) :: pt2(is:ie, km)
    REAL :: pt2_ad(is:ie, km)
    REAL, INTENT(IN) :: ws(is:ie)
    REAL :: ws_ad(is:ie)
    REAL, INTENT(INOUT) :: dz2(is:ie, km)
    REAL, INTENT(INOUT) :: dz2_ad(is:ie, km)
    REAL, INTENT(INOUT) :: w2(is:ie, km)
    REAL, INTENT(INOUT) :: w2_ad(is:ie, km)
    REAL :: pe2(is:ie, km+1)
    REAL :: pe2_ad(is:ie, km+1)
    REAL :: ws2(is:ie)
    REAL :: ws2_ad(is:ie)
    REAL, DIMENSION(km+1) :: m_bot, m_top, r_bot, r_top, pe1, pbar, wbar
    REAL, DIMENSION(km+1) :: m_bot_ad, m_top_ad, r_bot_ad, r_top_ad, &
&   pe1_ad, pbar_ad, wbar_ad
    REAL, DIMENSION(km) :: r_hi, r_lo, dz, wm, dm, dts
    REAL, DIMENSION(km) :: r_hi_ad, r_lo_ad, dz_ad, wm_ad, dm_ad, dts_ad
    REAL, DIMENSION(km) :: pf1, wc, cm, pp, pt1
    REAL, DIMENSION(km) :: pf1_ad, wc_ad, cm_ad, pp_ad, pt1_ad
    REAL :: dt, rdt, grg, z_frac, ptmp1, rden, pf, time_left
    REAL :: z_frac_ad, ptmp1_ad, rden_ad, pf_ad, time_left_ad
    REAL :: m_surf
    REAL :: m_surf_ad
    INTEGER :: i, k, n, ke, kt1, ktop
    INTEGER :: ks0, ks1
    INTRINSIC REAL
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC SQRT
    INTRINSIC MAX
    REAL :: temp
    REAL :: temp0
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
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
    REAL :: temp_ad15
    REAL :: temp_ad16
    REAL :: temp_ad17
    REAL :: temp_ad18
    REAL :: temp_ad19
    REAL :: temp_ad20
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_to
    INTEGER :: ad_to0
    INTEGER :: ad_to1
    INTEGER :: ad_to2
    INTEGER :: ad_from
    INTEGER :: ad_count0
    INTEGER :: i1
    INTEGER :: ad_from0
    INTEGER :: ad_to3
    INTEGER :: ad_from1
    INTEGER :: ad_count1
    INTEGER :: i2
    INTEGER :: ad_to4
    INTEGER :: ad_from2
    INTEGER :: ad_count2
    INTEGER :: i3
    INTEGER :: ad_count3
    INTEGER :: i4
    INTEGER :: ad_from3
    INTEGER :: ad_from4
    INTEGER :: ad_from5
    INTEGER :: ad_from6
    INTEGER :: ad_from7
    INTEGER :: ad_from8

    ws2 = 0.0
    m_bot = 0.0
    m_top = 0.0
    r_bot = 0.0
    r_top = 0.0
    pe1 = 0.0
    pbar = 0.0
    wbar = 0.0
    r_hi = 0.0
    r_lo = 0.0
    dz = 0.0
    wm = 0.0
    dm = 0.0
    dts = 0.0
    pf1 = 0.0
    wc = 0.0
    cm = 0.0
    pp = 0.0
    pt1 = 0.0
    dt = 0.0
    rdt = 0.0
    grg = 0.0
    z_frac = 0.0
    ptmp1 = 0.0
    rden = 0.0
    pf = 0.0
    time_left = 0.0
    m_surf = 0.0
    n = 0.0
    ke = 0.0
    kt1 = 0.0
    ktop = 0.0
    ks0 = 0.0
    ks1 = 0.0
    ad_count = 0.0
    ad_from = 0.0
    ad_count0 = 0.0
    ad_from0 = 0.0
    ad_from1 = 0.0
    ad_count1 = 0.0
    ad_from2 = 0.0
    ad_count2 = 0.0
    ad_count3 = 0.0
    ad_from3 = 0.0
    ad_from4 = 0.0
    ad_from5 = 0.0
    ad_from6 = 0.0
    ad_from7 = 0.0
    ad_from8 = 0.0
    ad_to = 0
    ad_to0 = 0
    ad_to1 = 0
    ad_to2 = 0
    ad_to3 = 0
    ad_to4 = 0
    branch = 0

    CALL POPREALARRAY(dm, km)
    CALL POPREALARRAY(dts, km)
    CALL POPREALARRAY(dt)
    CALL POPREALARRAY(grg)
    CALL POPREALARRAY(dz, km)
    CALL POPREALARRAY(r_lo, km)
    CALL POPREALARRAY(wbar, km + 1)
    CALL POPINTEGER(ks0)
    CALL POPREALARRAY(pf1, km)
    CALL POPREALARRAY(m_top, km + 1)
    CALL POPREALARRAY(r_top, km + 1)
    CALL POPREALARRAY(m_bot, km + 1)
    CALL POPREALARRAY(rdt)
    CALL POPREALARRAY(pf)
    CALL POPREALARRAY(r_bot, km + 1)
    CALL POPREALARRAY(ws2, ie - is + 1)
    CALL POPREALARRAY(time_left)
    CALL POPREALARRAY(pt1, km)
    CALL POPREALARRAY(cm, km)
    CALL POPREALARRAY(pp, km)
    CALL POPREALARRAY(wc, km)
    CALL POPREALARRAY(pbar, km + 1)
    CALL POPREALARRAY(m_surf)
    CALL POPREALARRAY(wm, km)
    CALL POPREALARRAY(r_hi, km)
    dm_ad = 0.0
    dts_ad = 0.0
    dz_ad = 0.0
    r_lo_ad = 0.0
    wbar_ad = 0.0
    pf1_ad = 0.0
    m_top_ad = 0.0
    r_top_ad = 0.0
    m_bot_ad = 0.0
    r_bot_ad = 0.0
    ws2_ad = 0.0
    pt1_ad = 0.0
    cm_ad = 0.0
    pp_ad = 0.0
    wc_ad = 0.0
    pbar_ad = 0.0
    wm_ad = 0.0
    r_hi_ad = 0.0
    DO i=ie,is,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO k=km+1,2,-1
          CALL POPREALARRAY(pe2(i, k))
          pbar_ad(k) = pbar_ad(k) + rdt*pe2_ad(i, k)
          pe2_ad(i, k) = 0.0
        END DO
        CALL POPREALARRAY(pe2(i, 1))
        pe2_ad(i, 1) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO k=km,1,-1
            CALL POPREALARRAY(dz2(i, k))
            dz_ad(k) = dz_ad(k) + dz2_ad(i, k)
            wbar_ad(k+1) = wbar_ad(k+1) + bdt*dz2_ad(i, k)
            wbar_ad(k) = wbar_ad(k) - bdt*dz2_ad(i, k)
            dz2_ad(i, k) = 0.0
          END DO
        ELSE
          DO k=km,1,-1
            CALL POPREALARRAY(w2(i, k))
            temp_ad9 = w2_ad(i, k)/dm(k)
            wm_ad(k) = wm_ad(k) + temp_ad9
            pbar_ad(k+1) = pbar_ad(k+1) + temp_ad9
            pbar_ad(k) = pbar_ad(k) - temp_ad9
            dm_ad(k) = dm_ad(k) - (wm(k)+pbar(k+1)-pbar(k))*temp_ad9/dm(&
&             k)
            w2_ad(i, k) = 0.0
            CALL POPREALARRAY(dz2(i, k))
            dz_ad(k) = dz_ad(k) + dz2_ad(i, k)
            wbar_ad(k+1) = wbar_ad(k+1) + bdt*dz2_ad(i, k)
            wbar_ad(k) = wbar_ad(k) - bdt*dz2_ad(i, k)
            dz2_ad(i, k) = 0.0
          END DO
        END IF
        CALL POPREALARRAY(pbar(km+1))
        temp_ad8 = bdt*pbar_ad(km+1)
        cm_ad(km) = cm_ad(km) + wbar(km+1)*temp_ad8
        wbar_ad(km+1) = wbar_ad(km+1) + cm(km)*temp_ad8
        wc_ad(km) = wc_ad(km) - temp_ad8
        pp_ad(km) = pp_ad(km) + temp_ad8
        pbar_ad(km+1) = 0.0
        pe1_ad = 0.0
      ELSE
        pe1_ad = 0.0
        DO k=km+1,2,-1
          CALL POPREALARRAY(pe2(i, k))
          pe1_ad(k) = pe1_ad(k) + rdt*pe2_ad(i, k)
          pe2_ad(i, k) = 0.0
        END DO
        CALL POPINTEGER(k)
        CALL POPREALARRAY(pe2(i, 1))
        pe2_ad(i, 1) = 0.0
        DO n=ms,1,-1
          CALL POPCONTROL(2,branch)
          IF (branch .EQ. 0) THEN
            CALL POPINTEGER(ad_from6)
            DO k=km,ad_from6,-1
              CALL POPREALARRAY(dz2(i, k))
              dz_ad(k) = dz_ad(k) + dz2_ad(i, k)
              wbar_ad(k+1) = wbar_ad(k+1) + dt*dz2_ad(i, k)
              wbar_ad(k) = wbar_ad(k) - dt*dz2_ad(i, k)
              dz2_ad(i, k) = 0.0
            END DO
          ELSE IF (branch .EQ. 1) THEN
            CALL POPINTEGER(ad_from7)
            DO k=km,ad_from7,-1
              CALL POPREALARRAY(w2(i, k))
              temp_ad20 = w2_ad(i, k)/dm(k)
              wm_ad(k) = wm_ad(k) + temp_ad20
              pbar_ad(k+1) = pbar_ad(k+1) + temp_ad20
              pbar_ad(k) = pbar_ad(k) - temp_ad20
              dm_ad(k) = dm_ad(k) - (wm(k)+pbar(k+1)-pbar(k))*temp_ad20/&
&               dm(k)
              w2_ad(i, k) = 0.0
              CALL POPREALARRAY(dz2(i, k))
              dz_ad(k) = dz_ad(k) + dz2_ad(i, k)
              wbar_ad(k+1) = wbar_ad(k+1) + dt*dz2_ad(i, k)
              wbar_ad(k) = wbar_ad(k) - dt*dz2_ad(i, k)
              dz2_ad(i, k) = 0.0
            END DO
          ELSE
            CALL POPINTEGER(ad_from8)
            DO k=km,ad_from8,-1
              CALL POPREALARRAY(wm(k))
              pbar_ad(k+1) = pbar_ad(k+1) + wm_ad(k)
              pbar_ad(k) = pbar_ad(k) - wm_ad(k)
              CALL POPREALARRAY(dz(k))
              wbar_ad(k+1) = wbar_ad(k+1) + dt*dz_ad(k)
              wbar_ad(k) = wbar_ad(k) - dt*dz_ad(k)
            END DO
          END IF
          CALL POPINTEGER(ad_from5)
          DO k=km+1,ad_from5,-1
            pbar_ad(k) = pbar_ad(k) + pe1_ad(k)
            CALL POPREALARRAY(pbar(k))
            m_top_ad(k) = m_top_ad(k) + wbar(k)*pbar_ad(k)
            wbar_ad(k) = wbar_ad(k) + m_top(k)*pbar_ad(k)
            r_top_ad(k) = r_top_ad(k) - pbar_ad(k)
            pbar_ad(k) = 0.0
          END DO
          CALL POPINTEGER(ad_from4)
          DO k=km,ad_from4,-1
            CALL POPREALARRAY(wbar(k))
            temp_ad18 = wbar_ad(k)/(m_top(k)+m_bot(k))
            temp_ad19 = -((r_bot(k)+r_top(k))*temp_ad18/(m_top(k)+m_bot(&
&             k)))
            r_bot_ad(k) = r_bot_ad(k) + temp_ad18
            r_top_ad(k) = r_top_ad(k) + temp_ad18
            m_top_ad(k) = m_top_ad(k) + temp_ad19
            m_bot_ad(k) = m_bot_ad(k) + temp_ad19
            wbar_ad(k) = 0.0
          END DO
          CALL POPINTEGER(k)
          CALL POPCONTROL(1,branch)
          IF (branch .NE. 0) THEN
            CALL POPREALARRAY(wbar(1))
            temp_ad17 = wbar_ad(1)/m_bot(1)
            r_bot_ad(1) = r_bot_ad(1) + temp_ad17
            m_bot_ad(1) = m_bot_ad(1) - r_bot(1)*temp_ad17/m_bot(1)
            wbar_ad(1) = 0.0
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .NE. 0) THEN
            CALL POPINTEGER(ad_from3)
            DO ke=km,ad_from3,-1
              CALL POPCONTROL(2,branch)
              IF (branch .EQ. 0) THEN
                z_frac = time_left/dts(k)
                CALL POPREALARRAY(r_bot(ke))
                z_frac_ad = dm(k)*m_bot_ad(ke) + r_lo(k)*r_bot_ad(ke)
                r_lo_ad(k) = r_lo_ad(k) + z_frac*r_bot_ad(ke)
                CALL POPREALARRAY(m_bot(ke))
                dm_ad(k) = dm_ad(k) + z_frac*m_bot_ad(ke)
                temp_ad16 = z_frac_ad/dts(k)
                time_left_ad = temp_ad16
                dts_ad(k) = dts_ad(k) - time_left*temp_ad16/dts(k)
              ELSE
                IF (branch .EQ. 1) THEN
                  m_bot_ad(ke) = m_bot_ad(ke) + ws2(i)*r_bot_ad(ke)
                  z_frac = time_left/dts(k)
                  CALL POPREALARRAY(r_bot(ke))
                  z_frac_ad = dm(k)*m_bot_ad(ke) - r_hi(k)*r_bot_ad(ke)
                  r_hi_ad(k) = r_hi_ad(k) - z_frac*r_bot_ad(ke)
                  m_surf_ad = -(ws2(i)*r_bot_ad(ke))
                  ws2_ad(i) = ws2_ad(i) + (m_bot(ke)-m_surf)*r_bot_ad(ke&
&                   )
                  CALL POPREALARRAY(m_bot(ke))
                  dm_ad(k) = dm_ad(k) + z_frac*m_bot_ad(ke)
                  temp_ad15 = z_frac_ad/dts(k)
                  time_left_ad = temp_ad15
                  dts_ad(k) = dts_ad(k) - time_left*temp_ad15/dts(k)
                ELSE
                  time_left_ad = 0.0
                  m_surf_ad = 0.0
                END IF
                CALL POPINTEGER(ad_count3)
                DO i4=1,ad_count3
                  IF (i4 .EQ. 1) THEN
                    CALL POPCONTROL(1,branch)
                    IF (branch .EQ. 0) THEN
                      time_left_ad = 0.0
                      m_surf_ad = 0.0
                    END IF
                  ELSE
                    CALL POPREALARRAY(r_bot(ke))
                    r_hi_ad(k) = r_hi_ad(k) - r_bot_ad(ke)
                    CALL POPREALARRAY(m_bot(ke))
                    dm_ad(k) = dm_ad(k) + m_bot_ad(ke)
                    dts_ad(k) = dts_ad(k) - time_left_ad
                  END IF
                  CALL POPINTEGER(k)
                END DO
                CALL POPREALARRAY(m_surf)
                m_bot_ad(ke) = m_bot_ad(ke) + m_surf_ad
              END IF
              CALL POPINTEGER(ad_count2)
              DO i3=1,ad_count2
                IF (i3 .EQ. 1) THEN
                  CALL POPCONTROL(1,branch)
                ELSE
                  CALL POPREALARRAY(r_bot(ke))
                  r_lo_ad(k) = r_lo_ad(k) + r_bot_ad(ke)
                  CALL POPREALARRAY(m_bot(ke))
                  dm_ad(k) = dm_ad(k) + m_bot_ad(ke)
                  dts_ad(k) = dts_ad(k) - time_left_ad
                END IF
                CALL POPINTEGER(k)
              END DO
              CALL POPREALARRAY(time_left)
            END DO
            CALL POPINTEGER(ad_from2)
            DO k=km,ad_from2,-1
              CALL POPREALARRAY(r_bot(k))
              r_bot_ad(k) = 0.0
              CALL POPREALARRAY(m_bot(k))
              m_bot_ad(k) = 0.0
            END DO
            CALL POPINTEGER(k)
            CALL POPINTEGER(ad_to4)
            DO ke=ad_to4,km+1,1
              CALL POPCONTROL(1,branch)
              IF (branch .EQ. 0) THEN
                z_frac = time_left/dts(k)
                CALL POPREALARRAY(r_top(ke))
                z_frac_ad = dm(k)*m_top_ad(ke) + r_hi(k)*r_top_ad(ke)
                r_hi_ad(k) = r_hi_ad(k) + z_frac*r_top_ad(ke)
                CALL POPREALARRAY(m_top(ke))
                dm_ad(k) = dm_ad(k) + z_frac*m_top_ad(ke)
                temp_ad14 = z_frac_ad/dts(k)
                time_left_ad = temp_ad14
                dts_ad(k) = dts_ad(k) - time_left*temp_ad14/dts(k)
              ELSE
                time_left_ad = 0.0
              END IF
              CALL POPINTEGER(ad_count1)
              DO i2=1,ad_count1
                IF (i2 .EQ. 1) THEN
                  CALL POPCONTROL(1,branch)
                  IF (branch .EQ. 0) time_left_ad = 0.0
                ELSE
                  CALL POPREALARRAY(r_top(ke))
                  r_hi_ad(k) = r_hi_ad(k) + r_top_ad(ke)
                  CALL POPREALARRAY(m_top(ke))
                  dm_ad(k) = dm_ad(k) + m_top_ad(ke)
                  dts_ad(k) = dts_ad(k) - time_left_ad
                END IF
                CALL POPINTEGER(k)
              END DO
              CALL POPREALARRAY(time_left)
            END DO
            CALL POPINTEGER(ad_from1)
            DO k=km+1,ad_from1,-1
              CALL POPREALARRAY(r_top(k))
              r_top_ad(k) = 0.0
              CALL POPREALARRAY(m_top(k))
              m_top_ad(k) = 0.0
            END DO
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) GOTO 100
          END IF
          CALL POPINTEGER(ad_from0)
          CALL POPINTEGER(ad_to3)
          DO k=ad_to3,ad_from0,-1
            CALL POPREALARRAY(m_top(k+1))
            m_bot_ad(k) = m_bot_ad(k) + m_top_ad(k+1)
            m_top_ad(k+1) = 0.0
            z_frac = dt/dts(k)
            CALL POPREALARRAY(m_bot(k))
            z_frac_ad = r_hi(k)*r_top_ad(k+1) + r_lo(k)*r_bot_ad(k) + dm&
&             (k)*m_bot_ad(k)
            dm_ad(k) = dm_ad(k) + z_frac*m_bot_ad(k)
            m_bot_ad(k) = 0.0
            CALL POPREALARRAY(r_top(k+1))
            r_hi_ad(k) = r_hi_ad(k) + z_frac*r_top_ad(k+1)
            r_top_ad(k+1) = 0.0
            CALL POPREALARRAY(r_bot(k))
            r_lo_ad(k) = r_lo_ad(k) + z_frac*r_bot_ad(k)
            r_bot_ad(k) = 0.0
            dts_ad(k) = dts_ad(k) - dt*z_frac_ad/dts(k)**2
          END DO
 100      CALL POPINTEGER(ad_count0)
          DO i1=1,ad_count0
            IF (i1 .EQ. 1) CALL POPCONTROL(1,branch)
          END DO
          CALL POPINTEGER(ad_from)
          DO k=km,ad_from,-1
            CALL POPREALARRAY(r_hi(k))
            wm_ad(k) = wm_ad(k) + r_lo_ad(k) + r_hi_ad(k)
            ptmp1_ad = r_lo_ad(k) - r_hi_ad(k)
            r_hi_ad(k) = 0.0
            CALL POPREALARRAY(r_lo(k))
            r_lo_ad(k) = 0.0
            dts_ad(k) = dts_ad(k) + (pf-pm2(i, k))*ptmp1_ad
            pm2_ad(i, k) = pm2_ad(i, k) - dts(k)*ptmp1_ad
            rden = -(rgas*dm(k)/dz(k))
            temp1 = pf/rden
            temp2 = SQRT(grg*temp1)
            IF (grg*temp1 .EQ. 0.0) THEN
              temp_ad11 = 0.0
            ELSE
              temp_ad11 = grg*dz(k)*dts_ad(k)/(2.0*temp2**3*rden)
            END IF
            pf_ad = temp_ad11 + dts(k)*ptmp1_ad
            CALL POPREALARRAY(dts(k))
            CALL POPREALARRAY(pf)
            temp_ad13 = gama*EXP(gama*LOG(rden*pt1(k)))*pf_ad/(rden*pt1(&
&             k))
            rden_ad = pt1(k)*temp_ad13 - temp1*temp_ad11
            pt1_ad(k) = pt1_ad(k) + rden*temp_ad13
            temp_ad12 = -(rgas*rden_ad/dz(k))
            dz_ad(k) = dz_ad(k) - dm(k)*temp_ad12/dz(k) - dts_ad(k)/&
&             temp2
            dts_ad(k) = 0.0
            dm_ad(k) = dm_ad(k) + temp_ad12
          END DO
          CALL POPINTEGER(k)
        END DO
        CALL POPCONTROL(2,branch)
        IF (branch .EQ. 0) THEN
          GOTO 110
        ELSE IF (branch .EQ. 1) THEN
          CALL POPREALARRAY(pbar(ks0))
          pbar_ad(ks0) = pbar_ad(ks0)/REAL(ms)
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPINTEGER(ad_to1)
            DO k=ad_to1,1,-1
              CALL POPREALARRAY(dz2(i, k))
              dz_ad(k) = dz_ad(k) + dz2_ad(i, k)
              wbar_ad(k+1) = wbar_ad(k+1) + bdt*dz2_ad(i, k)
              wbar_ad(k) = wbar_ad(k) - bdt*dz2_ad(i, k)
              dz2_ad(i, k) = 0.0
            END DO
          ELSE
            CALL POPINTEGER(ad_to2)
            DO k=ad_to2,1,-1
              CALL POPREALARRAY(w2(i, k))
              temp_ad10 = w2_ad(i, k)/dm(k)
              wm_ad(k) = wm_ad(k) + temp_ad10
              pbar_ad(k+1) = pbar_ad(k+1) + temp_ad10
              pbar_ad(k) = pbar_ad(k) - temp_ad10
              dm_ad(k) = dm_ad(k) - (wm(k)+pbar(k+1)-pbar(k))*temp_ad10/&
&               dm(k)
              w2_ad(i, k) = 0.0
              CALL POPREALARRAY(dz2(i, k))
              dz_ad(k) = dz_ad(k) + dz2_ad(i, k)
              wbar_ad(k+1) = wbar_ad(k+1) + bdt*dz2_ad(i, k)
              wbar_ad(k) = wbar_ad(k) - bdt*dz2_ad(i, k)
              dz2_ad(i, k) = 0.0
            END DO
          END IF
        ELSE
          GOTO 130
        END IF
      END IF
      CALL POPINTEGER(ad_to0)
      DO k=ad_to0,2,-1
        pbar_ad(k) = pbar_ad(k) + pe1_ad(k)
        pe1_ad(k) = 0.0
        CALL POPREALARRAY(pbar(k))
        temp_ad5 = bdt*pbar_ad(k)
        wbar_ad(k) = wbar_ad(k) + cm(k-1)*temp_ad5
        pp_ad(k-1) = pp_ad(k-1) + temp_ad5
        pbar_ad(k) = 0.0
        temp_ad7 = wbar_ad(k)/(cm(k-1)+cm(k))
        wc_ad(k-1) = wc_ad(k-1) + temp_ad7 - temp_ad5
        temp_ad6 = -((wc(k-1)+wc(k)+pp(k)-pp(k-1))*temp_ad7/(cm(k-1)+cm(&
&         k)))
        cm_ad(k-1) = cm_ad(k-1) + temp_ad6 + wbar(k)*temp_ad5
        CALL POPREALARRAY(wbar(k))
        wc_ad(k) = wc_ad(k) + temp_ad7
        pp_ad(k) = pp_ad(k) + temp_ad7
        pp_ad(k-1) = pp_ad(k-1) - temp_ad7
        cm_ad(k) = cm_ad(k) + temp_ad6
        wbar_ad(k) = 0.0
      END DO
      CALL POPREALARRAY(wbar(1))
      temp_ad4 = wbar_ad(1)/cm(1)
      wc_ad(1) = wc_ad(1) + temp_ad4
      pp_ad(1) = pp_ad(1) + temp_ad4
      cm_ad(1) = cm_ad(1) - (wc(1)+pp(1))*temp_ad4/cm(1)
      wbar_ad(1) = 0.0
      CALL POPINTEGER(ad_to)
      DO k=ad_to,1,-1
        temp_ad3 = cm_ad(k)/dts(k)
        CALL POPREALARRAY(pp(k))
        pf1_ad(k) = pf1_ad(k) + pp_ad(k)
        pm2_ad(i, k) = pm2_ad(i, k) - pp_ad(k)
        pp_ad(k) = 0.0
        CALL POPREALARRAY(wc(k))
        temp_ad2 = wc_ad(k)/dts(k)
        wm_ad(k) = wm_ad(k) + temp_ad2
        dts_ad(k) = dts_ad(k) - dm(k)*temp_ad3/dts(k) - wm(k)*temp_ad2/&
&         dts(k)
        wc_ad(k) = 0.0
        CALL POPREALARRAY(cm(k))
        dm_ad(k) = dm_ad(k) + temp_ad3
        cm_ad(k) = 0.0
      END DO
      CALL POPINTEGER(k)
 110  CALL POPINTEGER(ad_count)
      DO i0=1,ad_count
        IF (i0 .EQ. 1) THEN
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) GOTO 120
        END IF
        rden = -(rgas*dm(k)/dz(k))
        CALL POPREALARRAY(dts(k))
        temp = pf1(k)/rden
        temp0 = SQRT(grg*temp)
        IF (grg*temp .EQ. 0.0) THEN
          temp_ad = 0.0
        ELSE
          temp_ad = grg*dz(k)*dts_ad(k)/(2.0*temp0**3*rden)
        END IF
        pf1_ad(k) = pf1_ad(k) + temp_ad
        CALL POPREALARRAY(pf1(k))
        temp_ad1 = gama*EXP(gama*LOG(rden*pt1(k)))*pf1_ad(k)/(rden*pt1(k&
&         ))
        rden_ad = pt1(k)*temp_ad1 - temp*temp_ad
        pt1_ad(k) = pt1_ad(k) + rden*temp_ad1
        pf1_ad(k) = 0.0
        temp_ad0 = -(rgas*rden_ad/dz(k))
        dz_ad(k) = dz_ad(k) - dm(k)*temp_ad0/dz(k) - dts_ad(k)/temp0
        dts_ad(k) = 0.0
        dm_ad(k) = dm_ad(k) + temp_ad0
 120    CALL POPINTEGER(k)
      END DO
 130  CALL POPINTEGER(ks0)
      CALL POPREALARRAY(wbar(km+1))
      ws_ad(i) = ws_ad(i) + wbar_ad(km+1)
      wbar_ad(km+1) = 0.0
      DO k=km,1,-1
        CALL POPREALARRAY(pt1(k))
        pt2_ad(i, k) = pt2_ad(i, k) + pt1_ad(k)
        pt1_ad(k) = 0.0
        CALL POPREALARRAY(wm(k))
        w2_ad(i, k) = w2_ad(i, k) + dm(k)*wm_ad(k)
        dm_ad(k) = dm_ad(k) + w2(i, k)*wm_ad(k)
        wm_ad(k) = 0.0
        CALL POPREALARRAY(dm(k))
        dm2_ad(i, k) = dm2_ad(i, k) + dm_ad(k)
        dm_ad(k) = 0.0
        CALL POPREALARRAY(dz(k))
        dz2_ad(i, k) = dz2_ad(i, k) + dz_ad(k)
        dz_ad(k) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      ws_ad(i) = ws_ad(i) + 2.*ws2_ad(i)
      ws2_ad(i) = 0.0
    END DO
  END SUBROUTINE RIM_2D_BWD
  SUBROUTINE RIM_2D(ms, bdt, is, ie, km, rgas, gama, gm2, pe2, dm2, pm2&
&   , w2, dz2, pt2, ws, c_core)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ms, is, ie, km
    REAL, INTENT(IN) :: bdt, gama, rgas
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pm2, gm2
    LOGICAL, INTENT(IN) :: c_core
    REAL, INTENT(IN) :: pt2(is:ie, km)
    REAL, INTENT(IN) :: ws(is:ie)
! IN/OUT:
    REAL, INTENT(INOUT) :: dz2(is:ie, km)
    REAL, INTENT(INOUT) :: w2(is:ie, km)
    REAL, INTENT(OUT) :: pe2(is:ie, km+1)
! Local:
    REAL :: ws2(is:ie)
    REAL, DIMENSION(km+1) :: m_bot, m_top, r_bot, r_top, pe1, pbar, wbar
    REAL, DIMENSION(km) :: r_hi, r_lo, dz, wm, dm, dts
    REAL, DIMENSION(km) :: pf1, wc, cm, pp, pt1
    REAL :: dt, rdt, grg, z_frac, ptmp1, rden, pf, time_left
    REAL :: m_surf
    INTEGER :: i, k, n, ke, kt1, ktop
    INTEGER :: ks0, ks1
    INTRINSIC REAL
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC SQRT
    INTRINSIC MAX
    grg = gama*rgas
    rdt = 1./bdt
    dt = bdt/REAL(ms)
    pbar(:) = 0.
    wbar(:) = 0.
    DO i=is,ie
      ws2(i) = 2.*ws(i)
    END DO
! end i-loop
    DO i=is,ie
      DO k=1,km
        dz(k) = dz2(i, k)
        dm(k) = dm2(i, k)
        wm(k) = w2(i, k)*dm(k)
        pt1(k) = pt2(i, k)
      END DO
      pe1(:) = 0.
      wbar(km+1) = ws(i)
      ks0 = 1
      IF (ms .GT. 1 .AND. ms .LT. 8) THEN
! Continuity of (pbar, wbar) is maintained
        DO k=1,km
          rden = -(rgas*dm(k)/dz(k))
          pf1(k) = EXP(gama*LOG(rden*pt1(k)))
          dts(k) = -(dz(k)/SQRT(grg*pf1(k)/rden))
          IF (bdt .GT. dts(k)) THEN
            ks0 = k - 1
            GOTO 222
          END IF
        END DO
        ks0 = km
 222    IF (ks0 .NE. 1) THEN
          DO k=1,ks0
            cm(k) = dm(k)/dts(k)
            wc(k) = wm(k)/dts(k)
            pp(k) = pf1(k) - pm2(i, k)
          END DO
          wbar(1) = (wc(1)+pp(1))/cm(1)
          DO k=2,ks0
            wbar(k) = (wc(k-1)+wc(k)+pp(k)-pp(k-1))/(cm(k-1)+cm(k))
            pbar(k) = bdt*(cm(k-1)*wbar(k)-wc(k-1)+pp(k-1))
            pe1(k) = pbar(k)
          END DO
          IF (ks0 .EQ. km) THEN
            pbar(km+1) = bdt*(cm(km)*wbar(km+1)-wc(km)+pp(km))
            IF (c_core) THEN
              DO k=1,km
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
              END DO
            ELSE
              DO k=1,km
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
                w2(i, k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
              END DO
            END IF
            pe2(i, 1) = 0.
            DO k=2,km+1
              pe2(i, k) = pbar(k)*rdt
            END DO
            GOTO 6000
          ELSE
! next i
            IF (c_core) THEN
              DO k=1,ks0-1
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
              END DO
            ELSE
              DO k=1,ks0-1
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
                w2(i, k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
              END DO
            END IF
            pbar(ks0) = pbar(ks0)/REAL(ms)
          END IF
        END IF
      END IF
      ks1 = ks0
      DO n=1,ms
        DO k=ks1,km
          rden = -(rgas*dm(k)/dz(k))
          pf = EXP(gama*LOG(rden*pt1(k)))
          dts(k) = -(dz(k)/SQRT(grg*pf/rden))
          ptmp1 = dts(k)*(pf-pm2(i, k))
          r_lo(k) = wm(k) + ptmp1
          r_hi(k) = wm(k) - ptmp1
        END DO
        ktop = ks1
        DO k=ks1,km
          IF (dt .GT. dts(k)) THEN
            ktop = k - 1
            GOTO 333
          END IF
        END DO
        ktop = km
 333    IF (ktop .GE. ks1) THEN
          DO k=ks1,ktop
            z_frac = dt/dts(k)
            r_bot(k) = z_frac*r_lo(k)
            r_top(k+1) = z_frac*r_hi(k)
            m_bot(k) = z_frac*dm(k)
            m_top(k+1) = m_bot(k)
          END DO
          IF (ktop .EQ. km) GOTO 666
        END IF
        DO k=ktop+2,km+1
          m_top(k) = 0.
          r_top(k) = 0.
        END DO
        IF (1 .LT. ktop) THEN
          kt1 = ktop
        ELSE
          kt1 = 1
        END IF
        DO ke=km+1,ktop+2,-1
          time_left = dt
          DO k=ke-1,kt1,-1
            IF (time_left .GT. dts(k)) THEN
              time_left = time_left - dts(k)
              m_top(ke) = m_top(ke) + dm(k)
              r_top(ke) = r_top(ke) + r_hi(k)
            ELSE
              z_frac = time_left/dts(k)
              m_top(ke) = m_top(ke) + z_frac*dm(k)
              r_top(ke) = r_top(ke) + z_frac*r_hi(k)
              GOTO 444
            END IF
          END DO
 444      CONTINUE
        END DO
! next level
        DO k=ktop+1,km
          m_bot(k) = 0.
          r_bot(k) = 0.
        END DO
        DO ke=ktop+1,km
          time_left = dt
          DO k=ke,km
            IF (time_left .GT. dts(k)) THEN
              time_left = time_left - dts(k)
              m_bot(ke) = m_bot(ke) + dm(k)
              r_bot(ke) = r_bot(ke) + r_lo(k)
            ELSE
              z_frac = time_left/dts(k)
              m_bot(ke) = m_bot(ke) + z_frac*dm(k)
              r_bot(ke) = r_bot(ke) + z_frac*r_lo(k)
              GOTO 4000
            END IF
          END DO
! next interface
          m_surf = m_bot(ke)
          DO k=km,kt1,-1
            IF (time_left .GT. dts(k)) THEN
              time_left = time_left - dts(k)
              m_bot(ke) = m_bot(ke) + dm(k)
              r_bot(ke) = r_bot(ke) - r_hi(k)
            ELSE
              z_frac = time_left/dts(k)
              m_bot(ke) = m_bot(ke) + z_frac*dm(k)
              r_bot(ke) = r_bot(ke) - z_frac*r_hi(k) + (m_bot(ke)-m_surf&
&               )*ws2(i)
              GOTO 4000
            END IF
          END DO
 4000     CONTINUE
        END DO
! next interface
 666    IF (ks1 .EQ. 1) wbar(1) = r_bot(1)/m_bot(1)
        DO k=ks1+1,km
          wbar(k) = (r_bot(k)+r_top(k))/(m_top(k)+m_bot(k))
        END DO
! pbar here is actually dt*pbar
        DO k=ks1+1,km+1
          pbar(k) = m_top(k)*wbar(k) - r_top(k)
          pe1(k) = pe1(k) + pbar(k)
        END DO
        IF (n .EQ. ms) THEN
          IF (c_core) THEN
            DO k=ks1,km
              dz2(i, k) = dz(k) + dt*(wbar(k+1)-wbar(k))
            END DO
          ELSE
            DO k=ks1,km
              dz2(i, k) = dz(k) + dt*(wbar(k+1)-wbar(k))
              w2(i, k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
            END DO
          END IF
        ELSE
          DO k=ks1,km
            dz(k) = dz(k) + dt*(wbar(k+1)-wbar(k))
            wm(k) = wm(k) + pbar(k+1) - pbar(k)
          END DO
        END IF
      END DO
      pe2(i, 1) = 0.
      DO k=2,km+1
        pe2(i, k) = pe1(k)*rdt
      END DO
 6000 CONTINUE
    END DO
  END SUBROUTINE RIM_2D
!  Differentiation of sim3_solver in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: ws dm pe2 dz2 w2 pem pt2
!   with respect to varying inputs: ws dm pe2 dz2 w2 pem pt2
  SUBROUTINE SIM3_SOLVER_FWD(dt, is, ie, km, rgas, gama, kappa, pe2, dm&
&   , pem, w2, dz2, pt2, ws, alpha, p_fac, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, alpha, p_fac, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm, pt2
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL :: pe2(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, wk, g_rat, gam
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL :: beta, t2, t1g, rdt, ra, capa1, r2g, r6g
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX

    aa = 0.0
    bb = 0.0
    dd = 0.0
    w1 = 0.0
    wk = 0.0
    g_rat = 0.0
    gam = 0.0
    pp = 0.0
    p1 = 0.0
    wk1 = 0.0
    bet = 0.0
    beta = 0.0
    t2 = 0.0
    t1g = 0.0
    rdt = 0.0
    ra = 0.0
    capa1 = 0.0
    r2g = 0.0
    r6g = 0.0

    beta = 1. - alpha
    ra = 1./alpha
    t2 = beta/alpha
    t1g = gama*2.*(alpha*dt)**2
    rdt = 1./dt
    capa1 = kappa - 1.
    r2g = grav/2.
    r6g = grav/6.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
! Full pressure at center
        aa(i, k) = EXP(gama*LOG(-(dm(i, k)/dz2(i, k)*rgas*pt2(i, k))))
      END DO
    END DO
    DO k=1,km-1
      DO i=is,ie
! for profile reconstruction
        g_rat(i, k) = dm(i, k)/dm(i, k+1)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd(i, k) = 3.*(aa(i, k)+g_rat(i, k)*aa(i, k+1))
      END DO
    END DO
! pe2 is full p at edges
    DO i=is,ie
! Top:
      bet(i) = bb(i, 1)
      CALL PUSHREALARRAY(pe2(i, 1))
      pe2(i, 1) = pem(i, 1)
      CALL PUSHREALARRAY(pe2(i, 2))
      pe2(i, 2) = (dd(i, 1)-pem(i, 1))/bet(i)
! Bottom:
      bb(i, km) = 2.
      CALL PUSHREALARRAY(dd(i, km))
      dd(i, km) = 3.*aa(i, km) + r2g*dm(i, km)
    END DO
    DO k=2,km
      DO i=is,ie
        gam(i, k) = g_rat(i, k-1)/bet(i)
        CALL PUSHREALARRAY(bet(i))
        bet(i) = bb(i, k) - gam(i, k)
        CALL PUSHREALARRAY(pe2(i, k+1))
        pe2(i, k+1) = (dd(i, k)-pe2(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        CALL PUSHREALARRAY(pe2(i, k))
        pe2(i, k) = pe2(i, k) - gam(i, k)*pe2(i, k+1)
      END DO
    END DO
! done reconstruction of full:
! pp is pert. p at edges
    DO k=1,km+1
      DO i=is,ie
        pp(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
    DO k=2,km
      DO i=is,ie
        CALL PUSHREALARRAY(aa(i, k))
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*pe2(i, k)
        wk(i, k) = t2*aa(i, k)*(w1(i, k-1)-w1(i, k))
        CALL PUSHREALARRAY(aa(i, k))
        aa(i, k) = aa(i, k) - scale_m*dm(i, 1)
      END DO
    END DO
    DO i=is,ie
      CALL PUSHREALARRAY(bet(i))
      bet(i) = dm(i, 1) - aa(i, 2)
      CALL PUSHREALARRAY(w2(i, 1))
      w2(i, 1) = (dm(i, 1)*w1(i, 1)+dt*pp(i, 2)+wk(i, 2))/bet(i)
    END DO
    DO k=2,km-1
      DO i=is,ie
        CALL PUSHREALARRAY(gam(i, k))
        gam(i, k) = aa(i, k)/bet(i)
        CALL PUSHREALARRAY(bet(i))
        bet(i) = dm(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        CALL PUSHREALARRAY(w2(i, k))
        w2(i, k) = (dm(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))+wk(i, k+1&
&         )-wk(i, k)-aa(i, k)*w2(i, k-1))/bet(i)
      END DO
    END DO
    DO i=is,ie
      wk1(i) = t1g/dz2(i, km)*pe2(i, km+1)
      CALL PUSHREALARRAY(gam(i, km))
      gam(i, km) = aa(i, km)/bet(i)
      CALL PUSHREALARRAY(bet(i))
      bet(i) = dm(i, km) - (aa(i, km)+wk1(i)+aa(i, km)*gam(i, km))
      CALL PUSHREALARRAY(w2(i, km))
      w2(i, km) = (dm(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk(i, &
&       km)+wk1(i)*(t2*w1(i, km)-ra*ws(i))-aa(i, km)*w2(i, km-1))/bet(i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        CALL PUSHREALARRAY(w2(i, k))
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
! pe2 is updated perturbation p at edges
    DO i=is,ie
      CALL PUSHREALARRAY(pe2(i, 1))
      pe2(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        CALL PUSHREALARRAY(pe2(i, k+1))
        pe2(i, k+1) = pe2(i, k) + (dm(i, k)*(w2(i, k)-w1(i, k))*rdt-beta&
&         *(pp(i, k+1)-pp(i, k)))*ra
      END DO
    END DO
! Full non-hydro pressure at edges:
    DO i=is,ie
      CALL PUSHREALARRAY(pe2(i, 1))
      pe2(i, 1) = pem(i, 1)
    END DO
    DO k=2,km+1
      DO i=is,ie
        IF (p_fac*pem(i, k) .LT. pe2(i, k) + pem(i, k)) THEN
          CALL PUSHREALARRAY(pe2(i, k))
          pe2(i, k) = pe2(i, k) + pem(i, k)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(pe2(i, k))
          pe2(i, k) = p_fac*pem(i, k)
          CALL PUSHCONTROL(1,1)
        END IF
      END DO
    END DO
    DO i=is,ie
! Recover cell-averaged pressure
      p1(i) = (pe2(i, km)+2.*pe2(i, km+1))*r3 - r6g*dm(i, km)
      CALL PUSHREALARRAY(dz2(i, km))
      dz2(i, km) = -(dm(i, km)*rgas*pt2(i, km)*EXP(capa1*LOG(p1(i))))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        CALL PUSHREALARRAY(p1(i))
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        CALL PUSHREALARRAY(dz2(i, k))
        dz2(i, k) = -(dm(i, k)*rgas*pt2(i, k)*EXP(capa1*LOG(p1(i))))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        CALL PUSHREALARRAY(pe2(i, k))
        pe2(i, k) = pe2(i, k) - pem(i, k)
        pe2(i, k) = pe2(i, k) + beta*(pp(i, k)-pe2(i, k))
      END DO
    END DO
    CALL PUSHREALARRAY(gam, (ie-is+1)*km)
    CALL PUSHREALARRAY(t1g)
    CALL PUSHREALARRAY(wk, (ie-is+1)*km)
    CALL PUSHREALARRAY(capa1)
    CALL PUSHREALARRAY(g_rat, (ie-is+1)*km)
    CALL PUSHREALARRAY(pp, (ie-is+1)*(km+1))
    CALL PUSHREALARRAY(beta)
    CALL PUSHREALARRAY(rdt)
    CALL PUSHREALARRAY(t2)
    CALL PUSHREALARRAY(w1, (ie-is+1)*km)
    CALL PUSHREALARRAY(bb, (ie-is+1)*km)
    CALL PUSHREALARRAY(ra)
    CALL PUSHREALARRAY(r6g)
    CALL PUSHREALARRAY(wk1, ie - is + 1)
    CALL PUSHREALARRAY(bet, ie - is + 1)
    CALL PUSHREALARRAY(p1, ie - is + 1)
    CALL PUSHREALARRAY(r2g)
    CALL PUSHREALARRAY(aa, (ie-is+1)*km)
    CALL PUSHREALARRAY(dd, (ie-is+1)*km)
  END SUBROUTINE SIM3_SOLVER_FWD
!  Differentiation of sim3_solver in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: ws dm pe2 dz2 w2 pem pt2
!   with respect to varying inputs: ws dm pe2 dz2 w2 pem pt2
  SUBROUTINE SIM3_SOLVER_BWD(dt, is, ie, km, rgas, gama, kappa, pe2, &
&   pe2_ad, dm, dm_ad, pem, pem_ad, w2, w2_ad, dz2, dz2_ad, pt2, pt2_ad&
&   , ws, ws_ad, alpha, p_fac, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, alpha, p_fac, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm, pt2
    REAL, DIMENSION(is:ie, km) :: dm_ad, pt2_ad
    REAL, INTENT(IN) :: ws(is:ie)
    REAL :: ws_ad(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL, DIMENSION(is:ie, km+1) :: pem_ad
    REAL :: pe2(is:ie, km+1)
    REAL :: pe2_ad(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2_ad, w2_ad
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, wk, g_rat, gam
    REAL, DIMENSION(is:ie, km) :: aa_ad, bb_ad, dd_ad, w1_ad, wk_ad, &
&   g_rat_ad, gam_ad
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie, km+1) :: pp_ad
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL, DIMENSION(is:ie) :: p1_ad, wk1_ad, bet_ad
    REAL :: beta, t2, t1g, rdt, ra, capa1, r2g, r6g
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp2
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp3
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: temp4
    REAL :: temp_ad16
    REAL :: temp5
    REAL :: temp_ad17
    REAL :: temp_ad18
    INTEGER :: branch

    aa = 0.0
    bb = 0.0
    dd = 0.0
    w1 = 0.0
    wk = 0.0
    g_rat = 0.0
    gam = 0.0
    pp = 0.0
    p1 = 0.0
    wk1 = 0.0
    bet = 0.0
    beta = 0.0
    t2 = 0.0
    t1g = 0.0
    rdt = 0.0
    ra = 0.0
    capa1 = 0.0
    r2g = 0.0
    r6g = 0.0
    branch = 0

    CALL POPREALARRAY(dd, (ie-is+1)*km)
    CALL POPREALARRAY(aa, (ie-is+1)*km)
    CALL POPREALARRAY(r2g)
    CALL POPREALARRAY(p1, ie - is + 1)
    CALL POPREALARRAY(bet, ie - is + 1)
    CALL POPREALARRAY(wk1, ie - is + 1)
    CALL POPREALARRAY(r6g)
    CALL POPREALARRAY(ra)
    CALL POPREALARRAY(bb, (ie-is+1)*km)
    CALL POPREALARRAY(w1, (ie-is+1)*km)
    CALL POPREALARRAY(t2)
    CALL POPREALARRAY(rdt)
    CALL POPREALARRAY(beta)
    CALL POPREALARRAY(pp, (ie-is+1)*(km+1))
    CALL POPREALARRAY(g_rat, (ie-is+1)*km)
    CALL POPREALARRAY(capa1)
    CALL POPREALARRAY(wk, (ie-is+1)*km)
    CALL POPREALARRAY(t1g)
    CALL POPREALARRAY(gam, (ie-is+1)*km)
    pp_ad = 0.0
    DO k=km+1,1,-1
      DO i=ie,is,-1
        pp_ad(i, k) = pp_ad(i, k) + beta*pe2_ad(i, k)
        pe2_ad(i, k) = (1.0-beta)*pe2_ad(i, k)
        CALL POPREALARRAY(pe2(i, k))
        pem_ad(i, k) = pem_ad(i, k) - pe2_ad(i, k)
      END DO
    END DO
    p1_ad = 0.0
    bb_ad = 0.0
    g_rat_ad = 0.0
    DO k=1,km-1,1
      DO i=ie,is,-1
        CALL POPREALARRAY(dz2(i, k))
        temp5 = capa1*LOG(p1(i))
        temp_ad17 = -(rgas*EXP(temp5)*dz2_ad(i, k))
        dm_ad(i, k) = dm_ad(i, k) + pt2(i, k)*temp_ad17
        pt2_ad(i, k) = pt2_ad(i, k) + dm(i, k)*temp_ad17
        p1_ad(i) = p1_ad(i) - capa1*EXP(temp5)*dm(i, k)*pt2(i, k)*rgas*&
&         dz2_ad(i, k)/p1(i)
        dz2_ad(i, k) = 0.0
        CALL POPREALARRAY(p1(i))
        temp_ad18 = r3*p1_ad(i)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad18
        bb_ad(i, k) = bb_ad(i, k) + pe2(i, k+1)*temp_ad18
        pe2_ad(i, k+1) = pe2_ad(i, k+1) + bb(i, k)*temp_ad18
        g_rat_ad(i, k) = g_rat_ad(i, k) + pe2(i, k+2)*temp_ad18 - p1(i)*&
&         p1_ad(i)
        pe2_ad(i, k+2) = pe2_ad(i, k+2) + g_rat(i, k)*temp_ad18
        p1_ad(i) = -(g_rat(i, k)*p1_ad(i))
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(dz2(i, km))
      temp4 = capa1*LOG(p1(i))
      temp_ad16 = -(rgas*EXP(temp4)*dz2_ad(i, km))
      pt2_ad(i, km) = pt2_ad(i, km) + dm(i, km)*temp_ad16
      p1_ad(i) = p1_ad(i) - capa1*EXP(temp4)*dm(i, km)*pt2(i, km)*rgas*&
&       dz2_ad(i, km)/p1(i)
      dm_ad(i, km) = dm_ad(i, km) + pt2(i, km)*temp_ad16 - r6g*p1_ad(i)
      dz2_ad(i, km) = 0.0
      pe2_ad(i, km) = pe2_ad(i, km) + r3*p1_ad(i)
      pe2_ad(i, km+1) = pe2_ad(i, km+1) + r3*2.*p1_ad(i)
      p1_ad(i) = 0.0
    END DO
    DO k=km+1,2,-1
      DO i=ie,is,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(pe2(i, k))
          pem_ad(i, k) = pem_ad(i, k) + pe2_ad(i, k)
        ELSE
          CALL POPREALARRAY(pe2(i, k))
          pem_ad(i, k) = pem_ad(i, k) + p_fac*pe2_ad(i, k)
          pe2_ad(i, k) = 0.0
        END IF
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(pe2(i, 1))
      pem_ad(i, 1) = pem_ad(i, 1) + pe2_ad(i, 1)
      pe2_ad(i, 1) = 0.0
    END DO
    w1_ad = 0.0
    DO k=km,1,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k+1))
        temp_ad14 = ra*pe2_ad(i, k+1)
        temp_ad15 = rdt*dm(i, k)*temp_ad14
        pe2_ad(i, k) = pe2_ad(i, k) + pe2_ad(i, k+1)
        dm_ad(i, k) = dm_ad(i, k) + rdt*(w2(i, k)-w1(i, k))*temp_ad14
        w2_ad(i, k) = w2_ad(i, k) + temp_ad15
        w1_ad(i, k) = w1_ad(i, k) - temp_ad15
        pp_ad(i, k+1) = pp_ad(i, k+1) - beta*temp_ad14
        pp_ad(i, k) = pp_ad(i, k) + beta*temp_ad14
        pe2_ad(i, k+1) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(pe2(i, 1))
      pe2_ad(i, 1) = 0.0
    END DO
    gam_ad = 0.0
    DO k=1,km-1,1
      DO i=ie,is,-1
        CALL POPREALARRAY(w2(i, k))
        gam_ad(i, k+1) = gam_ad(i, k+1) - w2(i, k+1)*w2_ad(i, k)
        w2_ad(i, k+1) = w2_ad(i, k+1) - gam(i, k+1)*w2_ad(i, k)
      END DO
    END DO
    aa_ad = 0.0
    bet_ad = 0.0
    wk1_ad = 0.0
    wk_ad = 0.0
    DO i=ie,is,-1
      CALL POPREALARRAY(w2(i, km))
      temp_ad11 = w2_ad(i, km)/bet(i)
      temp3 = t2*w1(i, km) - ra*ws(i)
      w1_ad(i, km) = w1_ad(i, km) + (wk1(i)*t2+dm(i, km))*temp_ad11
      pp_ad(i, km+1) = pp_ad(i, km+1) + dt*temp_ad11
      pp_ad(i, km) = pp_ad(i, km) - dt*temp_ad11
      wk_ad(i, km) = wk_ad(i, km) - temp_ad11
      ws_ad(i) = ws_ad(i) - wk1(i)*ra*temp_ad11
      w2_ad(i, km-1) = w2_ad(i, km-1) - aa(i, km)*temp_ad11
      bet_ad(i) = bet_ad(i) - (dm(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i&
&       , km))-wk(i, km)+wk1(i)*temp3-aa(i, km)*w2(i, km-1))*temp_ad11/&
&       bet(i)
      dm_ad(i, km) = dm_ad(i, km) + bet_ad(i) + w1(i, km)*temp_ad11
      wk1_ad(i) = wk1_ad(i) + temp3*temp_ad11 - bet_ad(i)
      w2_ad(i, km) = 0.0
      CALL POPREALARRAY(bet(i))
      gam_ad(i, km) = gam_ad(i, km) - aa(i, km)*bet_ad(i)
      temp_ad12 = gam_ad(i, km)/bet(i)
      aa_ad(i, km) = aa_ad(i, km) + ((-1.0)-gam(i, km))*bet_ad(i) + &
&       temp_ad12 - w2(i, km-1)*temp_ad11
      bet_ad(i) = -(aa(i, km)*temp_ad12/bet(i))
      CALL POPREALARRAY(gam(i, km))
      gam_ad(i, km) = 0.0
      temp_ad13 = t1g*wk1_ad(i)/dz2(i, km)
      pe2_ad(i, km+1) = pe2_ad(i, km+1) + temp_ad13
      dz2_ad(i, km) = dz2_ad(i, km) - pe2(i, km+1)*temp_ad13/dz2(i, km)
      wk1_ad(i) = 0.0
    END DO
    DO k=km-1,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(w2(i, k))
        temp_ad9 = w2_ad(i, k)/bet(i)
        w1_ad(i, k) = w1_ad(i, k) + dm(i, k)*temp_ad9
        pp_ad(i, k+1) = pp_ad(i, k+1) + dt*temp_ad9
        pp_ad(i, k) = pp_ad(i, k) - dt*temp_ad9
        wk_ad(i, k+1) = wk_ad(i, k+1) + temp_ad9
        wk_ad(i, k) = wk_ad(i, k) - temp_ad9
        w2_ad(i, k-1) = w2_ad(i, k-1) - aa(i, k)*temp_ad9
        bet_ad(i) = bet_ad(i) - (dm(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, &
&         k))+wk(i, k+1)-wk(i, k)-aa(i, k)*w2(i, k-1))*temp_ad9/bet(i)
        dm_ad(i, k) = dm_ad(i, k) + bet_ad(i) + w1(i, k)*temp_ad9
        aa_ad(i, k) = aa_ad(i, k) + ((-1.0)-gam(i, k))*bet_ad(i) - w2(i&
&         , k-1)*temp_ad9
        w2_ad(i, k) = 0.0
        CALL POPREALARRAY(bet(i))
        aa_ad(i, k+1) = aa_ad(i, k+1) - bet_ad(i)
        gam_ad(i, k) = gam_ad(i, k) - aa(i, k)*bet_ad(i)
        CALL POPREALARRAY(gam(i, k))
        temp_ad10 = gam_ad(i, k)/bet(i)
        bet_ad(i) = -(aa(i, k)*temp_ad10/bet(i))
        aa_ad(i, k) = aa_ad(i, k) + temp_ad10
        gam_ad(i, k) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(w2(i, 1))
      temp_ad8 = w2_ad(i, 1)/bet(i)
      w1_ad(i, 1) = w1_ad(i, 1) + dm(i, 1)*temp_ad8
      pp_ad(i, 2) = pp_ad(i, 2) + dt*temp_ad8
      wk_ad(i, 2) = wk_ad(i, 2) + temp_ad8
      bet_ad(i) = bet_ad(i) - (dm(i, 1)*w1(i, 1)+dt*pp(i, 2)+wk(i, 2))*&
&       temp_ad8/bet(i)
      dm_ad(i, 1) = dm_ad(i, 1) + bet_ad(i) + w1(i, 1)*temp_ad8
      w2_ad(i, 1) = 0.0
      CALL POPREALARRAY(bet(i))
      aa_ad(i, 2) = aa_ad(i, 2) - bet_ad(i)
      bet_ad(i) = 0.0
    END DO
    DO k=km,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(aa(i, k))
        dm_ad(i, 1) = dm_ad(i, 1) - scale_m*aa_ad(i, k)
        temp_ad5 = t2*aa(i, k)*wk_ad(i, k)
        aa_ad(i, k) = aa_ad(i, k) + t2*(w1(i, k-1)-w1(i, k))*wk_ad(i, k)
        w1_ad(i, k-1) = w1_ad(i, k-1) + temp_ad5
        w1_ad(i, k) = w1_ad(i, k) - temp_ad5
        wk_ad(i, k) = 0.0
        CALL POPREALARRAY(aa(i, k))
        temp2 = dz2(i, k-1) + dz2(i, k)
        temp_ad6 = t1g*aa_ad(i, k)/temp2
        temp_ad7 = -(pe2(i, k)*temp_ad6/temp2)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad6
        dz2_ad(i, k-1) = dz2_ad(i, k-1) + temp_ad7
        dz2_ad(i, k) = dz2_ad(i, k) + temp_ad7
        aa_ad(i, k) = 0.0
      END DO
    END DO
    DO k=km+1,1,-1
      DO i=ie,is,-1
        pe2_ad(i, k) = pe2_ad(i, k) + pp_ad(i, k)
        pem_ad(i, k) = pem_ad(i, k) - pp_ad(i, k)
        pp_ad(i, k) = 0.0
      END DO
    END DO
    DO k=2,km,1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k))
        gam_ad(i, k) = gam_ad(i, k) - pe2(i, k+1)*pe2_ad(i, k)
        pe2_ad(i, k+1) = pe2_ad(i, k+1) - gam(i, k)*pe2_ad(i, k)
      END DO
    END DO
    dd_ad = 0.0
    DO k=km,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k+1))
        temp_ad3 = pe2_ad(i, k+1)/bet(i)
        dd_ad(i, k) = dd_ad(i, k) + temp_ad3
        pe2_ad(i, k) = pe2_ad(i, k) - temp_ad3
        bet_ad(i) = bet_ad(i) - (dd(i, k)-pe2(i, k))*temp_ad3/bet(i)
        pe2_ad(i, k+1) = 0.0
        CALL POPREALARRAY(bet(i))
        bb_ad(i, k) = bb_ad(i, k) + bet_ad(i)
        gam_ad(i, k) = gam_ad(i, k) - bet_ad(i)
        temp_ad4 = gam_ad(i, k)/bet(i)
        bet_ad(i) = -(g_rat(i, k-1)*temp_ad4/bet(i))
        g_rat_ad(i, k-1) = g_rat_ad(i, k-1) + temp_ad4
        gam_ad(i, k) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(dd(i, km))
      aa_ad(i, km) = aa_ad(i, km) + 3.*dd_ad(i, km)
      dm_ad(i, km) = dm_ad(i, km) + r2g*dd_ad(i, km)
      dd_ad(i, km) = 0.0
      bb_ad(i, km) = 0.0
      CALL POPREALARRAY(pe2(i, 2))
      temp_ad2 = pe2_ad(i, 2)/bet(i)
      dd_ad(i, 1) = dd_ad(i, 1) + temp_ad2
      bet_ad(i) = bet_ad(i) - (dd(i, 1)-pem(i, 1))*temp_ad2/bet(i)
      pe2_ad(i, 2) = 0.0
      pem_ad(i, 1) = pem_ad(i, 1) + pe2_ad(i, 1) - temp_ad2
      CALL POPREALARRAY(pe2(i, 1))
      pe2_ad(i, 1) = 0.0
      bb_ad(i, 1) = bb_ad(i, 1) + bet_ad(i)
      bet_ad(i) = 0.0
    END DO
    DO k=km-1,1,-1
      DO i=ie,is,-1
        temp_ad0 = 3.*dd_ad(i, k)
        aa_ad(i, k) = aa_ad(i, k) + temp_ad0
        g_rat_ad(i, k) = g_rat_ad(i, k) + 2.*bb_ad(i, k) + aa(i, k+1)*&
&         temp_ad0
        aa_ad(i, k+1) = aa_ad(i, k+1) + g_rat(i, k)*temp_ad0
        dd_ad(i, k) = 0.0
        bb_ad(i, k) = 0.0
        temp_ad1 = g_rat_ad(i, k)/dm(i, k+1)
        dm_ad(i, k) = dm_ad(i, k) + temp_ad1
        dm_ad(i, k+1) = dm_ad(i, k+1) - dm(i, k)*temp_ad1/dm(i, k+1)
        g_rat_ad(i, k) = 0.0
      END DO
    END DO
    DO k=km,1,-1
      DO i=ie,is,-1
        temp1 = dz2(i, k)
        temp0 = dm(i, k)*pt2(i, k)
        temp = temp0/temp1
        temp_ad = gama*EXP(gama*LOG(-(rgas*temp)))*aa_ad(i, k)/(temp*&
&         temp1)
        dm_ad(i, k) = dm_ad(i, k) + pt2(i, k)*temp_ad
        pt2_ad(i, k) = pt2_ad(i, k) + dm(i, k)*temp_ad
        dz2_ad(i, k) = dz2_ad(i, k) - temp*temp_ad
        aa_ad(i, k) = 0.0
        w2_ad(i, k) = w2_ad(i, k) + w1_ad(i, k)
        w1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE SIM3_SOLVER_BWD
  SUBROUTINE SIM3_SOLVER(dt, is, ie, km, rgas, gama, kappa, pe2, dm, pem&
&   , w2, dz2, pt2, ws, alpha, p_fac, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, alpha, p_fac, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm, pt2
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL, INTENT(OUT) :: pe2(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, wk, g_rat, gam
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL :: beta, t2, t1g, rdt, ra, capa1, r2g, r6g
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    beta = 1. - alpha
    ra = 1./alpha
    t2 = beta/alpha
    t1g = gama*2.*(alpha*dt)**2
    rdt = 1./dt
    capa1 = kappa - 1.
    r2g = grav/2.
    r6g = grav/6.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
! Full pressure at center
        aa(i, k) = EXP(gama*LOG(-(dm(i, k)/dz2(i, k)*rgas*pt2(i, k))))
      END DO
    END DO
    DO k=1,km-1
      DO i=is,ie
! for profile reconstruction
        g_rat(i, k) = dm(i, k)/dm(i, k+1)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd(i, k) = 3.*(aa(i, k)+g_rat(i, k)*aa(i, k+1))
      END DO
    END DO
! pe2 is full p at edges
    DO i=is,ie
! Top:
      bet(i) = bb(i, 1)
      pe2(i, 1) = pem(i, 1)
      pe2(i, 2) = (dd(i, 1)-pem(i, 1))/bet(i)
! Bottom:
      bb(i, km) = 2.
      dd(i, km) = 3.*aa(i, km) + r2g*dm(i, km)
    END DO
    DO k=2,km
      DO i=is,ie
        gam(i, k) = g_rat(i, k-1)/bet(i)
        bet(i) = bb(i, k) - gam(i, k)
        pe2(i, k+1) = (dd(i, k)-pe2(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        pe2(i, k) = pe2(i, k) - gam(i, k)*pe2(i, k+1)
      END DO
    END DO
! done reconstruction of full:
! pp is pert. p at edges
    DO k=1,km+1
      DO i=is,ie
        pp(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
    DO k=2,km
      DO i=is,ie
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*pe2(i, k)
        wk(i, k) = t2*aa(i, k)*(w1(i, k-1)-w1(i, k))
        aa(i, k) = aa(i, k) - scale_m*dm(i, 1)
      END DO
    END DO
    DO i=is,ie
      bet(i) = dm(i, 1) - aa(i, 2)
      w2(i, 1) = (dm(i, 1)*w1(i, 1)+dt*pp(i, 2)+wk(i, 2))/bet(i)
    END DO
    DO k=2,km-1
      DO i=is,ie
        gam(i, k) = aa(i, k)/bet(i)
        bet(i) = dm(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        w2(i, k) = (dm(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))+wk(i, k+1&
&         )-wk(i, k)-aa(i, k)*w2(i, k-1))/bet(i)
      END DO
    END DO
    DO i=is,ie
      wk1(i) = t1g/dz2(i, km)*pe2(i, km+1)
      gam(i, km) = aa(i, km)/bet(i)
      bet(i) = dm(i, km) - (aa(i, km)+wk1(i)+aa(i, km)*gam(i, km))
      w2(i, km) = (dm(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk(i, &
&       km)+wk1(i)*(t2*w1(i, km)-ra*ws(i))-aa(i, km)*w2(i, km-1))/bet(i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
! pe2 is updated perturbation p at edges
    DO i=is,ie
      pe2(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        pe2(i, k+1) = pe2(i, k) + (dm(i, k)*(w2(i, k)-w1(i, k))*rdt-beta&
&         *(pp(i, k+1)-pp(i, k)))*ra
      END DO
    END DO
! Full non-hydro pressure at edges:
    DO i=is,ie
      pe2(i, 1) = pem(i, 1)
    END DO
    DO k=2,km+1
      DO i=is,ie
        IF (p_fac*pem(i, k) .LT. pe2(i, k) + pem(i, k)) THEN
          pe2(i, k) = pe2(i, k) + pem(i, k)
        ELSE
          pe2(i, k) = p_fac*pem(i, k)
        END IF
      END DO
    END DO
    DO i=is,ie
! Recover cell-averaged pressure
      p1(i) = (pe2(i, km)+2.*pe2(i, km+1))*r3 - r6g*dm(i, km)
      dz2(i, km) = -(dm(i, km)*rgas*pt2(i, km)*EXP(capa1*LOG(p1(i))))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        dz2(i, k) = -(dm(i, k)*rgas*pt2(i, k)*EXP(capa1*LOG(p1(i))))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        pe2(i, k) = pe2(i, k) - pem(i, k)
        pe2(i, k) = pe2(i, k) + beta*(pp(i, k)-pe2(i, k))
      END DO
    END DO
  END SUBROUTINE SIM3_SOLVER
!  Differentiation of sim3p0_solver in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: ws dm pe2 dz2 w2 pem pt2
!   with respect to varying inputs: ws dm pe2 dz2 w2 pem pt2
  SUBROUTINE SIM3P0_SOLVER_FWD(dt, is, ie, km, rgas, gama, kappa, pe2, &
&   dm, pem, w2, dz2, pt2, ws, p_fac, scale_m)
    IMPLICIT NONE
! Sa SIM3, but for beta==0
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm, pt2
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, INTENT(IN) :: pem(is:ie, km+1)
    REAL :: pe2(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, g_rat, gam
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL :: t1g, rdt, capa1, r2g, r6g
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX

    aa = 0.0
    bb = 0.0
    dd = 0.0
    w1 = 0.0
    g_rat = 0.0
    gam = 0.0
    pp = 0.0
    p1 = 0.0
    wk1 = 0.0
    bet = 0.0
    t1g = 0.0
    rdt = 0.0
    capa1 = 0.0
    r2g = 0.0
    r6g = 0.0

    t1g = 2.*gama*dt**2
    rdt = 1./dt
    capa1 = kappa - 1.
    r2g = grav/2.
    r6g = grav/6.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
! Full pressure at center
        aa(i, k) = EXP(gama*LOG(-(dm(i, k)/dz2(i, k)*rgas*pt2(i, k))))
      END DO
    END DO
    DO k=1,km-1
      DO i=is,ie
! for profile reconstruction
        g_rat(i, k) = dm(i, k)/dm(i, k+1)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd(i, k) = 3.*(aa(i, k)+g_rat(i, k)*aa(i, k+1))
      END DO
    END DO
! pe2 is full p at edges
    DO i=is,ie
! Top:
      bet(i) = bb(i, 1)
      CALL PUSHREALARRAY(pe2(i, 1))
      pe2(i, 1) = pem(i, 1)
      CALL PUSHREALARRAY(pe2(i, 2))
      pe2(i, 2) = (dd(i, 1)-pem(i, 1))/bet(i)
! Bottom:
      bb(i, km) = 2.
      CALL PUSHREALARRAY(dd(i, km))
      dd(i, km) = 3.*aa(i, km) + r2g*dm(i, km)
    END DO
    DO k=2,km
      DO i=is,ie
        gam(i, k) = g_rat(i, k-1)/bet(i)
        CALL PUSHREALARRAY(bet(i))
        bet(i) = bb(i, k) - gam(i, k)
        CALL PUSHREALARRAY(pe2(i, k+1))
        pe2(i, k+1) = (dd(i, k)-pe2(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        CALL PUSHREALARRAY(pe2(i, k))
        pe2(i, k) = pe2(i, k) - gam(i, k)*pe2(i, k+1)
      END DO
    END DO
! done reconstruction of full:
! pp is pert. p at edges
    DO k=1,km+1
      DO i=is,ie
        pp(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
    DO k=2,km
      DO i=is,ie
        CALL PUSHREALARRAY(aa(i, k))
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*pe2(i, k) - scale_m*dm(i&
&         , 1)
      END DO
    END DO
    DO i=is,ie
      CALL PUSHREALARRAY(bet(i))
      bet(i) = dm(i, 1) - aa(i, 2)
      CALL PUSHREALARRAY(w2(i, 1))
      w2(i, 1) = (dm(i, 1)*w1(i, 1)+dt*pp(i, 2))/bet(i)
    END DO
    DO k=2,km-1
      DO i=is,ie
        CALL PUSHREALARRAY(gam(i, k))
        gam(i, k) = aa(i, k)/bet(i)
        CALL PUSHREALARRAY(bet(i))
        bet(i) = dm(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        CALL PUSHREALARRAY(w2(i, k))
        w2(i, k) = (dm(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))-aa(i, k)*&
&         w2(i, k-1))/bet(i)
      END DO
    END DO
    DO i=is,ie
      wk1(i) = t1g/dz2(i, km)*pe2(i, km+1)
      CALL PUSHREALARRAY(gam(i, km))
      gam(i, km) = aa(i, km)/bet(i)
      CALL PUSHREALARRAY(bet(i))
      bet(i) = dm(i, km) - (aa(i, km)+wk1(i)+aa(i, km)*gam(i, km))
      CALL PUSHREALARRAY(w2(i, km))
      w2(i, km) = (dm(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk1(i)&
&       *ws(i)-aa(i, km)*w2(i, km-1))/bet(i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        CALL PUSHREALARRAY(w2(i, k))
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
! pe2 is updated perturbation p at edges
    DO i=is,ie
      CALL PUSHREALARRAY(pe2(i, 1))
      pe2(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        CALL PUSHREALARRAY(pe2(i, k+1))
        pe2(i, k+1) = pe2(i, k) + dm(i, k)*(w2(i, k)-w1(i, k))*rdt
      END DO
    END DO
! Full non-hydro pressure at edges:
    DO i=is,ie
      CALL PUSHREALARRAY(pe2(i, 1))
      pe2(i, 1) = pem(i, 1)
    END DO
    DO k=2,km+1
      DO i=is,ie
        IF (p_fac*pem(i, k) .LT. pe2(i, k) + pem(i, k)) THEN
          CALL PUSHREALARRAY(pe2(i, k))
          pe2(i, k) = pe2(i, k) + pem(i, k)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(pe2(i, k))
          pe2(i, k) = p_fac*pem(i, k)
          CALL PUSHCONTROL(1,1)
        END IF
      END DO
    END DO
    DO i=is,ie
! Recover cell-averaged pressure
      p1(i) = (pe2(i, km)+2.*pe2(i, km+1))*r3 - r6g*dm(i, km)
      CALL PUSHREALARRAY(dz2(i, km))
      dz2(i, km) = -(dm(i, km)*rgas*pt2(i, km)*EXP(capa1*LOG(p1(i))))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        CALL PUSHREALARRAY(p1(i))
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        CALL PUSHREALARRAY(dz2(i, k))
        dz2(i, k) = -(dm(i, k)*rgas*pt2(i, k)*EXP(capa1*LOG(p1(i))))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        CALL PUSHREALARRAY(pe2(i, k))
        pe2(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
    CALL PUSHREALARRAY(gam, (ie-is+1)*km)
    CALL PUSHREALARRAY(t1g)
    CALL PUSHREALARRAY(capa1)
    CALL PUSHREALARRAY(g_rat, (ie-is+1)*km)
    CALL PUSHREALARRAY(pp, (ie-is+1)*(km+1))
    CALL PUSHREALARRAY(rdt)
    CALL PUSHREALARRAY(w1, (ie-is+1)*km)
    CALL PUSHREALARRAY(bb, (ie-is+1)*km)
    CALL PUSHREALARRAY(r6g)
    CALL PUSHREALARRAY(wk1, ie - is + 1)
    CALL PUSHREALARRAY(bet, ie - is + 1)
    CALL PUSHREALARRAY(p1, ie - is + 1)
    CALL PUSHREALARRAY(r2g)
    CALL PUSHREALARRAY(aa, (ie-is+1)*km)
    CALL PUSHREALARRAY(dd, (ie-is+1)*km)
  END SUBROUTINE SIM3P0_SOLVER_FWD
!  Differentiation of sim3p0_solver in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
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
!   gradient     of useful results: ws dm pe2 dz2 w2 pem pt2
!   with respect to varying inputs: ws dm pe2 dz2 w2 pem pt2
  SUBROUTINE SIM3P0_SOLVER_BWD(dt, is, ie, km, rgas, gama, kappa, pe2, &
&   pe2_ad, dm, dm_ad, pem, pem_ad, w2, w2_ad, dz2, dz2_ad, pt2, pt2_ad&
&   , ws, ws_ad, p_fac, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm, pt2
    REAL, DIMENSION(is:ie, km) :: dm_ad, pt2_ad
    REAL, INTENT(IN) :: ws(is:ie)
    REAL :: ws_ad(is:ie)
    REAL, INTENT(IN) :: pem(is:ie, km+1)
    REAL :: pem_ad(is:ie, km+1)
    REAL :: pe2(is:ie, km+1)
    REAL :: pe2_ad(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2_ad, w2_ad
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, g_rat, gam
    REAL, DIMENSION(is:ie, km) :: aa_ad, bb_ad, dd_ad, w1_ad, g_rat_ad, &
&   gam_ad
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie, km+1) :: pp_ad
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL, DIMENSION(is:ie) :: p1_ad, wk1_ad, bet_ad
    REAL :: t1g, rdt, capa1, r2g, r6g
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp2
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp3
    REAL :: temp_ad14
    REAL :: temp4
    REAL :: temp_ad15
    REAL :: temp_ad16
    INTEGER :: branch

    aa = 0.0
    bb = 0.0
    dd = 0.0
    w1 = 0.0
    g_rat = 0.0
    gam = 0.0
    pp = 0.0
    p1 = 0.0
    wk1 = 0.0
    bet = 0.0
    t1g = 0.0
    rdt = 0.0
    capa1 = 0.0
    r2g = 0.0
    r6g = 0.0
    branch = 0

    CALL POPREALARRAY(dd, (ie-is+1)*km)
    CALL POPREALARRAY(aa, (ie-is+1)*km)
    CALL POPREALARRAY(r2g)
    CALL POPREALARRAY(p1, ie - is + 1)
    CALL POPREALARRAY(bet, ie - is + 1)
    CALL POPREALARRAY(wk1, ie - is + 1)
    CALL POPREALARRAY(r6g)
    CALL POPREALARRAY(bb, (ie-is+1)*km)
    CALL POPREALARRAY(w1, (ie-is+1)*km)
    CALL POPREALARRAY(rdt)
    CALL POPREALARRAY(pp, (ie-is+1)*(km+1))
    CALL POPREALARRAY(g_rat, (ie-is+1)*km)
    CALL POPREALARRAY(capa1)
    CALL POPREALARRAY(t1g)
    CALL POPREALARRAY(gam, (ie-is+1)*km)
    DO k=km+1,1,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k))
        pem_ad(i, k) = pem_ad(i, k) - pe2_ad(i, k)
      END DO
    END DO
    p1_ad = 0.0
    bb_ad = 0.0
    g_rat_ad = 0.0
    DO k=1,km-1,1
      DO i=ie,is,-1
        CALL POPREALARRAY(dz2(i, k))
        temp4 = capa1*LOG(p1(i))
        temp_ad15 = -(rgas*EXP(temp4)*dz2_ad(i, k))
        dm_ad(i, k) = dm_ad(i, k) + pt2(i, k)*temp_ad15
        pt2_ad(i, k) = pt2_ad(i, k) + dm(i, k)*temp_ad15
        p1_ad(i) = p1_ad(i) - capa1*EXP(temp4)*dm(i, k)*pt2(i, k)*rgas*&
&         dz2_ad(i, k)/p1(i)
        dz2_ad(i, k) = 0.0
        CALL POPREALARRAY(p1(i))
        temp_ad16 = r3*p1_ad(i)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad16
        bb_ad(i, k) = bb_ad(i, k) + pe2(i, k+1)*temp_ad16
        pe2_ad(i, k+1) = pe2_ad(i, k+1) + bb(i, k)*temp_ad16
        g_rat_ad(i, k) = g_rat_ad(i, k) + pe2(i, k+2)*temp_ad16 - p1(i)*&
&         p1_ad(i)
        pe2_ad(i, k+2) = pe2_ad(i, k+2) + g_rat(i, k)*temp_ad16
        p1_ad(i) = -(g_rat(i, k)*p1_ad(i))
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(dz2(i, km))
      temp3 = capa1*LOG(p1(i))
      temp_ad14 = -(rgas*EXP(temp3)*dz2_ad(i, km))
      pt2_ad(i, km) = pt2_ad(i, km) + dm(i, km)*temp_ad14
      p1_ad(i) = p1_ad(i) - capa1*EXP(temp3)*dm(i, km)*pt2(i, km)*rgas*&
&       dz2_ad(i, km)/p1(i)
      dm_ad(i, km) = dm_ad(i, km) + pt2(i, km)*temp_ad14 - r6g*p1_ad(i)
      dz2_ad(i, km) = 0.0
      pe2_ad(i, km) = pe2_ad(i, km) + r3*p1_ad(i)
      pe2_ad(i, km+1) = pe2_ad(i, km+1) + r3*2.*p1_ad(i)
      p1_ad(i) = 0.0
    END DO
    DO k=km+1,2,-1
      DO i=ie,is,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(pe2(i, k))
          pem_ad(i, k) = pem_ad(i, k) + pe2_ad(i, k)
        ELSE
          CALL POPREALARRAY(pe2(i, k))
          pem_ad(i, k) = pem_ad(i, k) + p_fac*pe2_ad(i, k)
          pe2_ad(i, k) = 0.0
        END IF
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(pe2(i, 1))
      pem_ad(i, 1) = pem_ad(i, 1) + pe2_ad(i, 1)
      pe2_ad(i, 1) = 0.0
    END DO
    w1_ad = 0.0
    DO k=km,1,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k+1))
        temp_ad13 = rdt*dm(i, k)*pe2_ad(i, k+1)
        pe2_ad(i, k) = pe2_ad(i, k) + pe2_ad(i, k+1)
        dm_ad(i, k) = dm_ad(i, k) + rdt*(w2(i, k)-w1(i, k))*pe2_ad(i, k+&
&         1)
        w2_ad(i, k) = w2_ad(i, k) + temp_ad13
        w1_ad(i, k) = w1_ad(i, k) - temp_ad13
        pe2_ad(i, k+1) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(pe2(i, 1))
      pe2_ad(i, 1) = 0.0
    END DO
    gam_ad = 0.0
    DO k=1,km-1,1
      DO i=ie,is,-1
        CALL POPREALARRAY(w2(i, k))
        gam_ad(i, k+1) = gam_ad(i, k+1) - w2(i, k+1)*w2_ad(i, k)
        w2_ad(i, k+1) = w2_ad(i, k+1) - gam(i, k+1)*w2_ad(i, k)
      END DO
    END DO
    aa_ad = 0.0
    bet_ad = 0.0
    wk1_ad = 0.0
    pp_ad = 0.0
    DO i=ie,is,-1
      CALL POPREALARRAY(w2(i, km))
      temp_ad10 = w2_ad(i, km)/bet(i)
      w1_ad(i, km) = w1_ad(i, km) + dm(i, km)*temp_ad10
      pp_ad(i, km+1) = pp_ad(i, km+1) + dt*temp_ad10
      pp_ad(i, km) = pp_ad(i, km) - dt*temp_ad10
      ws_ad(i) = ws_ad(i) - wk1(i)*temp_ad10
      w2_ad(i, km-1) = w2_ad(i, km-1) - aa(i, km)*temp_ad10
      bet_ad(i) = bet_ad(i) - (dm(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i&
&       , km))-wk1(i)*ws(i)-aa(i, km)*w2(i, km-1))*temp_ad10/bet(i)
      dm_ad(i, km) = dm_ad(i, km) + bet_ad(i) + w1(i, km)*temp_ad10
      wk1_ad(i) = wk1_ad(i) - bet_ad(i) - ws(i)*temp_ad10
      w2_ad(i, km) = 0.0
      CALL POPREALARRAY(bet(i))
      gam_ad(i, km) = gam_ad(i, km) - aa(i, km)*bet_ad(i)
      temp_ad11 = gam_ad(i, km)/bet(i)
      aa_ad(i, km) = aa_ad(i, km) + ((-1.0)-gam(i, km))*bet_ad(i) + &
&       temp_ad11 - w2(i, km-1)*temp_ad10
      bet_ad(i) = -(aa(i, km)*temp_ad11/bet(i))
      CALL POPREALARRAY(gam(i, km))
      gam_ad(i, km) = 0.0
      temp_ad12 = t1g*wk1_ad(i)/dz2(i, km)
      pe2_ad(i, km+1) = pe2_ad(i, km+1) + temp_ad12
      dz2_ad(i, km) = dz2_ad(i, km) - pe2(i, km+1)*temp_ad12/dz2(i, km)
      wk1_ad(i) = 0.0
    END DO
    DO k=km-1,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(w2(i, k))
        temp_ad8 = w2_ad(i, k)/bet(i)
        w1_ad(i, k) = w1_ad(i, k) + dm(i, k)*temp_ad8
        pp_ad(i, k+1) = pp_ad(i, k+1) + dt*temp_ad8
        pp_ad(i, k) = pp_ad(i, k) - dt*temp_ad8
        w2_ad(i, k-1) = w2_ad(i, k-1) - aa(i, k)*temp_ad8
        bet_ad(i) = bet_ad(i) - (dm(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, &
&         k))-aa(i, k)*w2(i, k-1))*temp_ad8/bet(i)
        dm_ad(i, k) = dm_ad(i, k) + bet_ad(i) + w1(i, k)*temp_ad8
        aa_ad(i, k) = aa_ad(i, k) + ((-1.0)-gam(i, k))*bet_ad(i) - w2(i&
&         , k-1)*temp_ad8
        w2_ad(i, k) = 0.0
        CALL POPREALARRAY(bet(i))
        aa_ad(i, k+1) = aa_ad(i, k+1) - bet_ad(i)
        gam_ad(i, k) = gam_ad(i, k) - aa(i, k)*bet_ad(i)
        CALL POPREALARRAY(gam(i, k))
        temp_ad9 = gam_ad(i, k)/bet(i)
        bet_ad(i) = -(aa(i, k)*temp_ad9/bet(i))
        aa_ad(i, k) = aa_ad(i, k) + temp_ad9
        gam_ad(i, k) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(w2(i, 1))
      temp_ad7 = w2_ad(i, 1)/bet(i)
      w1_ad(i, 1) = w1_ad(i, 1) + dm(i, 1)*temp_ad7
      pp_ad(i, 2) = pp_ad(i, 2) + dt*temp_ad7
      bet_ad(i) = bet_ad(i) - (dm(i, 1)*w1(i, 1)+dt*pp(i, 2))*temp_ad7/&
&       bet(i)
      dm_ad(i, 1) = dm_ad(i, 1) + bet_ad(i) + w1(i, 1)*temp_ad7
      w2_ad(i, 1) = 0.0
      CALL POPREALARRAY(bet(i))
      aa_ad(i, 2) = aa_ad(i, 2) - bet_ad(i)
      bet_ad(i) = 0.0
    END DO
    DO k=km,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(aa(i, k))
        temp2 = dz2(i, k-1) + dz2(i, k)
        temp_ad5 = t1g*aa_ad(i, k)/temp2
        temp_ad6 = -(pe2(i, k)*temp_ad5/temp2)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad5
        dz2_ad(i, k-1) = dz2_ad(i, k-1) + temp_ad6
        dz2_ad(i, k) = dz2_ad(i, k) + temp_ad6
        dm_ad(i, 1) = dm_ad(i, 1) - scale_m*aa_ad(i, k)
        aa_ad(i, k) = 0.0
      END DO
    END DO
    DO k=km+1,1,-1
      DO i=ie,is,-1
        pe2_ad(i, k) = pe2_ad(i, k) + pp_ad(i, k)
        pem_ad(i, k) = pem_ad(i, k) - pp_ad(i, k)
        pp_ad(i, k) = 0.0
      END DO
    END DO
    DO k=2,km,1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k))
        gam_ad(i, k) = gam_ad(i, k) - pe2(i, k+1)*pe2_ad(i, k)
        pe2_ad(i, k+1) = pe2_ad(i, k+1) - gam(i, k)*pe2_ad(i, k)
      END DO
    END DO
    dd_ad = 0.0
    DO k=km,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k+1))
        temp_ad3 = pe2_ad(i, k+1)/bet(i)
        dd_ad(i, k) = dd_ad(i, k) + temp_ad3
        pe2_ad(i, k) = pe2_ad(i, k) - temp_ad3
        bet_ad(i) = bet_ad(i) - (dd(i, k)-pe2(i, k))*temp_ad3/bet(i)
        pe2_ad(i, k+1) = 0.0
        CALL POPREALARRAY(bet(i))
        bb_ad(i, k) = bb_ad(i, k) + bet_ad(i)
        gam_ad(i, k) = gam_ad(i, k) - bet_ad(i)
        temp_ad4 = gam_ad(i, k)/bet(i)
        bet_ad(i) = -(g_rat(i, k-1)*temp_ad4/bet(i))
        g_rat_ad(i, k-1) = g_rat_ad(i, k-1) + temp_ad4
        gam_ad(i, k) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(dd(i, km))
      aa_ad(i, km) = aa_ad(i, km) + 3.*dd_ad(i, km)
      dm_ad(i, km) = dm_ad(i, km) + r2g*dd_ad(i, km)
      dd_ad(i, km) = 0.0
      bb_ad(i, km) = 0.0
      CALL POPREALARRAY(pe2(i, 2))
      temp_ad2 = pe2_ad(i, 2)/bet(i)
      dd_ad(i, 1) = dd_ad(i, 1) + temp_ad2
      bet_ad(i) = bet_ad(i) - (dd(i, 1)-pem(i, 1))*temp_ad2/bet(i)
      pe2_ad(i, 2) = 0.0
      pem_ad(i, 1) = pem_ad(i, 1) + pe2_ad(i, 1) - temp_ad2
      CALL POPREALARRAY(pe2(i, 1))
      pe2_ad(i, 1) = 0.0
      bb_ad(i, 1) = bb_ad(i, 1) + bet_ad(i)
      bet_ad(i) = 0.0
    END DO
    DO k=km-1,1,-1
      DO i=ie,is,-1
        temp_ad0 = 3.*dd_ad(i, k)
        aa_ad(i, k) = aa_ad(i, k) + temp_ad0
        g_rat_ad(i, k) = g_rat_ad(i, k) + 2.*bb_ad(i, k) + aa(i, k+1)*&
&         temp_ad0
        aa_ad(i, k+1) = aa_ad(i, k+1) + g_rat(i, k)*temp_ad0
        dd_ad(i, k) = 0.0
        bb_ad(i, k) = 0.0
        temp_ad1 = g_rat_ad(i, k)/dm(i, k+1)
        dm_ad(i, k) = dm_ad(i, k) + temp_ad1
        dm_ad(i, k+1) = dm_ad(i, k+1) - dm(i, k)*temp_ad1/dm(i, k+1)
        g_rat_ad(i, k) = 0.0
      END DO
    END DO
    DO k=km,1,-1
      DO i=ie,is,-1
        temp1 = dz2(i, k)
        temp0 = dm(i, k)*pt2(i, k)
        temp = temp0/temp1
        temp_ad = gama*EXP(gama*LOG(-(rgas*temp)))*aa_ad(i, k)/(temp*&
&         temp1)
        dm_ad(i, k) = dm_ad(i, k) + pt2(i, k)*temp_ad
        pt2_ad(i, k) = pt2_ad(i, k) + dm(i, k)*temp_ad
        dz2_ad(i, k) = dz2_ad(i, k) - temp*temp_ad
        aa_ad(i, k) = 0.0
        w2_ad(i, k) = w2_ad(i, k) + w1_ad(i, k)
        w1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE SIM3P0_SOLVER_BWD
  SUBROUTINE SIM3P0_SOLVER(dt, is, ie, km, rgas, gama, kappa, pe2, dm, &
&   pem, w2, dz2, pt2, ws, p_fac, scale_m)
    IMPLICIT NONE
! Sa SIM3, but for beta==0
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm, pt2
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, INTENT(IN) :: pem(is:ie, km+1)
    REAL, INTENT(OUT) :: pe2(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, g_rat, gam
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL :: t1g, rdt, capa1, r2g, r6g
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    t1g = 2.*gama*dt**2
    rdt = 1./dt
    capa1 = kappa - 1.
    r2g = grav/2.
    r6g = grav/6.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
! Full pressure at center
        aa(i, k) = EXP(gama*LOG(-(dm(i, k)/dz2(i, k)*rgas*pt2(i, k))))
      END DO
    END DO
    DO k=1,km-1
      DO i=is,ie
! for profile reconstruction
        g_rat(i, k) = dm(i, k)/dm(i, k+1)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd(i, k) = 3.*(aa(i, k)+g_rat(i, k)*aa(i, k+1))
      END DO
    END DO
! pe2 is full p at edges
    DO i=is,ie
! Top:
      bet(i) = bb(i, 1)
      pe2(i, 1) = pem(i, 1)
      pe2(i, 2) = (dd(i, 1)-pem(i, 1))/bet(i)
! Bottom:
      bb(i, km) = 2.
      dd(i, km) = 3.*aa(i, km) + r2g*dm(i, km)
    END DO
    DO k=2,km
      DO i=is,ie
        gam(i, k) = g_rat(i, k-1)/bet(i)
        bet(i) = bb(i, k) - gam(i, k)
        pe2(i, k+1) = (dd(i, k)-pe2(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        pe2(i, k) = pe2(i, k) - gam(i, k)*pe2(i, k+1)
      END DO
    END DO
! done reconstruction of full:
! pp is pert. p at edges
    DO k=1,km+1
      DO i=is,ie
        pp(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
    DO k=2,km
      DO i=is,ie
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*pe2(i, k) - scale_m*dm(i&
&         , 1)
      END DO
    END DO
    DO i=is,ie
      bet(i) = dm(i, 1) - aa(i, 2)
      w2(i, 1) = (dm(i, 1)*w1(i, 1)+dt*pp(i, 2))/bet(i)
    END DO
    DO k=2,km-1
      DO i=is,ie
        gam(i, k) = aa(i, k)/bet(i)
        bet(i) = dm(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        w2(i, k) = (dm(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))-aa(i, k)*&
&         w2(i, k-1))/bet(i)
      END DO
    END DO
    DO i=is,ie
      wk1(i) = t1g/dz2(i, km)*pe2(i, km+1)
      gam(i, km) = aa(i, km)/bet(i)
      bet(i) = dm(i, km) - (aa(i, km)+wk1(i)+aa(i, km)*gam(i, km))
      w2(i, km) = (dm(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk1(i)&
&       *ws(i)-aa(i, km)*w2(i, km-1))/bet(i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
! pe2 is updated perturbation p at edges
    DO i=is,ie
      pe2(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        pe2(i, k+1) = pe2(i, k) + dm(i, k)*(w2(i, k)-w1(i, k))*rdt
      END DO
    END DO
! Full non-hydro pressure at edges:
    DO i=is,ie
      pe2(i, 1) = pem(i, 1)
    END DO
    DO k=2,km+1
      DO i=is,ie
        IF (p_fac*pem(i, k) .LT. pe2(i, k) + pem(i, k)) THEN
          pe2(i, k) = pe2(i, k) + pem(i, k)
        ELSE
          pe2(i, k) = p_fac*pem(i, k)
        END IF
      END DO
    END DO
    DO i=is,ie
! Recover cell-averaged pressure
      p1(i) = (pe2(i, km)+2.*pe2(i, km+1))*r3 - r6g*dm(i, km)
      dz2(i, km) = -(dm(i, km)*rgas*pt2(i, km)*EXP(capa1*LOG(p1(i))))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        dz2(i, k) = -(dm(i, k)*rgas*pt2(i, k)*EXP(capa1*LOG(p1(i))))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        pe2(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
  END SUBROUTINE SIM3P0_SOLVER
!  Differentiation of sim1_solver in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: ws dm2 dz2 w2 pm2 pem pe pt2
!   with respect to varying inputs: ws dm2 dz2 w2 pm2 pem pe pt2
  SUBROUTINE SIM1_SOLVER_FWD(dt, is, ie, km, rgas, gama, gm2, cp2, kappa&
&   , pe, dm2, pm2, pem, w2, dz2, pt2, ws, p_fac)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pt2, pm2, gm2, cp2
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL :: pe(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, g_rat, gam
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie) :: p1, bet
    REAL :: t1g, rdt, capa1
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: max1
    REAL :: max2

    aa = 0.0
    bb = 0.0
    dd = 0.0
    w1 = 0.0
    g_rat = 0.0
    gam = 0.0
    pp = 0.0
    p1 = 0.0
    bet = 0.0
    t1g = 0.0
    rdt = 0.0
    capa1 = 0.0
    max1 = 0.0
    max2 = 0.0

    t1g = gama*2.*dt*dt
    rdt = 1./dt
    capa1 = kappa - 1.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
        CALL PUSHREALARRAY(pe(i, k))
        pe(i, k) = EXP(gama*LOG(-(dm2(i, k)/dz2(i, k)*rgas*pt2(i, k)))) &
&         - pm2(i, k)
      END DO
    END DO
    DO k=1,km-1
      DO i=is,ie
        g_rat(i, k) = dm2(i, k)/dm2(i, k+1)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd(i, k) = 3.*(pe(i, k)+g_rat(i, k)*pe(i, k+1))
      END DO
    END DO
    DO i=is,ie
      bet(i) = bb(i, 1)
      pp(i, 1) = 0.
      pp(i, 2) = dd(i, 1)/bet(i)
      bb(i, km) = 2.
      CALL PUSHREALARRAY(dd(i, km))
      dd(i, km) = 3.*pe(i, km)
    END DO
    DO k=2,km
      DO i=is,ie
        gam(i, k) = g_rat(i, k-1)/bet(i)
        CALL PUSHREALARRAY(bet(i))
        bet(i) = bb(i, k) - gam(i, k)
        CALL PUSHREALARRAY(pp(i, k+1))
        pp(i, k+1) = (dd(i, k)-pp(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        CALL PUSHREALARRAY(pp(i, k))
        pp(i, k) = pp(i, k) - gam(i, k)*pp(i, k+1)
      END DO
    END DO
! Start the w-solver
    DO k=2,km
      DO i=is,ie
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*(pem(i, k)+pp(i, k))
      END DO
    END DO
    DO i=is,ie
      CALL PUSHREALARRAY(bet(i))
      bet(i) = dm2(i, 1) - aa(i, 2)
      CALL PUSHREALARRAY(w2(i, 1))
      w2(i, 1) = (dm2(i, 1)*w1(i, 1)+dt*pp(i, 2))/bet(i)
    END DO
    DO k=2,km-1
      DO i=is,ie
        CALL PUSHREALARRAY(gam(i, k))
        gam(i, k) = aa(i, k)/bet(i)
        CALL PUSHREALARRAY(bet(i))
        bet(i) = dm2(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        CALL PUSHREALARRAY(w2(i, k))
        w2(i, k) = (dm2(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))-aa(i, k)&
&         *w2(i, k-1))/bet(i)
      END DO
    END DO
    DO i=is,ie
      p1(i) = t1g/dz2(i, km)*(pem(i, km+1)+pp(i, km+1))
      CALL PUSHREALARRAY(gam(i, km))
      gam(i, km) = aa(i, km)/bet(i)
      CALL PUSHREALARRAY(bet(i))
      bet(i) = dm2(i, km) - (aa(i, km)+p1(i)+aa(i, km)*gam(i, km))
      CALL PUSHREALARRAY(w2(i, km))
      w2(i, km) = (dm2(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-p1(i)&
&       *ws(i)-aa(i, km)*w2(i, km-1))/bet(i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        CALL PUSHREALARRAY(w2(i, k))
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
    DO i=is,ie
      CALL PUSHREALARRAY(pe(i, 1))
      pe(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        CALL PUSHREALARRAY(pe(i, k+1))
        pe(i, k+1) = pe(i, k) + dm2(i, k)*(w2(i, k)-w1(i, k))*rdt
      END DO
    END DO
    DO i=is,ie
      CALL PUSHREALARRAY(p1(i))
      p1(i) = (pe(i, km)+2.*pe(i, km+1))*r3
      IF (p_fac*pm2(i, km) .LT. p1(i) + pm2(i, km)) THEN
        CALL PUSHREALARRAY(max1)
        max1 = p1(i) + pm2(i, km)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHREALARRAY(max1)
        max1 = p_fac*pm2(i, km)
        CALL PUSHCONTROL(1,1)
      END IF
      CALL PUSHREALARRAY(dz2(i, km))
      dz2(i, km) = -(dm2(i, km)*rgas*pt2(i, km)*EXP(capa1*LOG(max1)))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        CALL PUSHREALARRAY(p1(i))
        p1(i) = (pe(i, k)+bb(i, k)*pe(i, k+1)+g_rat(i, k)*pe(i, k+2))*r3&
&         - g_rat(i, k)*p1(i)
        IF (p_fac*pm2(i, k) .LT. p1(i) + pm2(i, k)) THEN
          CALL PUSHREALARRAY(max2)
          max2 = p1(i) + pm2(i, k)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(max2)
          max2 = p_fac*pm2(i, k)
          CALL PUSHCONTROL(1,1)
        END IF
        CALL PUSHREALARRAY(dz2(i, k))
        dz2(i, k) = -(dm2(i, k)*rgas*pt2(i, k)*EXP(capa1*LOG(max2)))
      END DO
    END DO
    CALL PUSHREALARRAY(gam, (ie-is+1)*km)
    CALL PUSHREALARRAY(t1g)
    CALL PUSHREALARRAY(capa1)
    CALL PUSHREALARRAY(g_rat, (ie-is+1)*km)
    CALL PUSHREALARRAY(pp, (ie-is+1)*(km+1))
    CALL PUSHREALARRAY(rdt)
    CALL PUSHREALARRAY(w1, (ie-is+1)*km)
    CALL PUSHREALARRAY(bb, (ie-is+1)*km)
    CALL PUSHREALARRAY(bet, ie - is + 1)
    CALL PUSHREALARRAY(p1, ie - is + 1)
    CALL PUSHREALARRAY(max2)
    CALL PUSHREALARRAY(max1)
    CALL PUSHREALARRAY(aa, (ie-is+1)*km)
    CALL PUSHREALARRAY(dd, (ie-is+1)*km)
  END SUBROUTINE SIM1_SOLVER_FWD
!  Differentiation of sim1_solver in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: ws dm2 dz2 w2 pm2 pem pe pt2
!   with respect to varying inputs: ws dm2 dz2 w2 pm2 pem pe pt2
  SUBROUTINE SIM1_SOLVER_BWD(dt, is, ie, km, rgas, gama, gm2, cp2, kappa&
&   , pe, pe_ad, dm2, dm2_ad, pm2, pm2_ad, pem, pem_ad, w2, w2_ad, dz2, &
&   dz2_ad, pt2, pt2_ad, ws, ws_ad, p_fac)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pt2, pm2, gm2, cp2
    REAL, DIMENSION(is:ie, km) :: dm2_ad, pt2_ad, pm2_ad
    REAL, INTENT(IN) :: ws(is:ie)
    REAL :: ws_ad(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL, DIMENSION(is:ie, km+1) :: pem_ad
    REAL :: pe(is:ie, km+1)
    REAL :: pe_ad(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2_ad, w2_ad
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, g_rat, gam
    REAL, DIMENSION(is:ie, km) :: aa_ad, bb_ad, dd_ad, w1_ad, g_rat_ad, &
&   gam_ad
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie, km+1) :: pp_ad
    REAL, DIMENSION(is:ie) :: p1, bet
    REAL, DIMENSION(is:ie) :: p1_ad, bet_ad
    REAL :: t1g, rdt, capa1
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: max1
    REAL :: max1_ad
    REAL :: max2
    REAL :: max2_ad
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp2
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp3
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: temp4
    REAL :: temp_ad16
    INTEGER :: branch

    aa = 0.0
    bb = 0.0
    dd = 0.0
    w1 = 0.0
    g_rat = 0.0
    gam = 0.0
    pp = 0.0
    p1 = 0.0
    bet = 0.0
    t1g = 0.0
    rdt = 0.0
    capa1 = 0.0
    branch = 0

    CALL POPREALARRAY(dd, (ie-is+1)*km)
    CALL POPREALARRAY(aa, (ie-is+1)*km)
    CALL POPREALARRAY(max1)
    CALL POPREALARRAY(max2)
    CALL POPREALARRAY(p1, ie - is + 1)
    CALL POPREALARRAY(bet, ie - is + 1)
    CALL POPREALARRAY(bb, (ie-is+1)*km)
    CALL POPREALARRAY(w1, (ie-is+1)*km)
    CALL POPREALARRAY(rdt)
    CALL POPREALARRAY(pp, (ie-is+1)*(km+1))
    CALL POPREALARRAY(g_rat, (ie-is+1)*km)
    CALL POPREALARRAY(capa1)
    CALL POPREALARRAY(t1g)
    CALL POPREALARRAY(gam, (ie-is+1)*km)
    p1_ad = 0.0
    bb_ad = 0.0
    g_rat_ad = 0.0
    DO k=1,km-1,1
      DO i=ie,is,-1
        CALL POPREALARRAY(dz2(i, k))
        temp4 = capa1*LOG(max2)
        temp_ad16 = -(rgas*EXP(temp4)*dz2_ad(i, k))
        dm2_ad(i, k) = dm2_ad(i, k) + pt2(i, k)*temp_ad16
        pt2_ad(i, k) = pt2_ad(i, k) + dm2(i, k)*temp_ad16
        max2_ad = -(capa1*EXP(temp4)*dm2(i, k)*pt2(i, k)*rgas*dz2_ad(i, &
&         k)/max2)
        dz2_ad(i, k) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(max2)
          p1_ad(i) = p1_ad(i) + max2_ad
          pm2_ad(i, k) = pm2_ad(i, k) + max2_ad
        ELSE
          CALL POPREALARRAY(max2)
          pm2_ad(i, k) = pm2_ad(i, k) + p_fac*max2_ad
        END IF
        CALL POPREALARRAY(p1(i))
        temp_ad15 = r3*p1_ad(i)
        pe_ad(i, k) = pe_ad(i, k) + temp_ad15
        bb_ad(i, k) = bb_ad(i, k) + pe(i, k+1)*temp_ad15
        pe_ad(i, k+1) = pe_ad(i, k+1) + bb(i, k)*temp_ad15
        g_rat_ad(i, k) = g_rat_ad(i, k) + pe(i, k+2)*temp_ad15 - p1(i)*&
&         p1_ad(i)
        pe_ad(i, k+2) = pe_ad(i, k+2) + g_rat(i, k)*temp_ad15
        p1_ad(i) = -(g_rat(i, k)*p1_ad(i))
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(dz2(i, km))
      temp3 = capa1*LOG(max1)
      temp_ad14 = -(rgas*EXP(temp3)*dz2_ad(i, km))
      dm2_ad(i, km) = dm2_ad(i, km) + pt2(i, km)*temp_ad14
      pt2_ad(i, km) = pt2_ad(i, km) + dm2(i, km)*temp_ad14
      max1_ad = -(capa1*EXP(temp3)*dm2(i, km)*pt2(i, km)*rgas*dz2_ad(i, &
&       km)/max1)
      dz2_ad(i, km) = 0.0
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(max1)
        p1_ad(i) = p1_ad(i) + max1_ad
        pm2_ad(i, km) = pm2_ad(i, km) + max1_ad
      ELSE
        CALL POPREALARRAY(max1)
        pm2_ad(i, km) = pm2_ad(i, km) + p_fac*max1_ad
      END IF
      CALL POPREALARRAY(p1(i))
      pe_ad(i, km) = pe_ad(i, km) + r3*p1_ad(i)
      pe_ad(i, km+1) = pe_ad(i, km+1) + r3*2.*p1_ad(i)
      p1_ad(i) = 0.0
    END DO
    w1_ad = 0.0
    DO k=km,1,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe(i, k+1))
        temp_ad13 = rdt*dm2(i, k)*pe_ad(i, k+1)
        pe_ad(i, k) = pe_ad(i, k) + pe_ad(i, k+1)
        dm2_ad(i, k) = dm2_ad(i, k) + rdt*(w2(i, k)-w1(i, k))*pe_ad(i, k&
&         +1)
        w2_ad(i, k) = w2_ad(i, k) + temp_ad13
        w1_ad(i, k) = w1_ad(i, k) - temp_ad13
        pe_ad(i, k+1) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(pe(i, 1))
      pe_ad(i, 1) = 0.0
    END DO
    gam_ad = 0.0
    DO k=1,km-1,1
      DO i=ie,is,-1
        CALL POPREALARRAY(w2(i, k))
        gam_ad(i, k+1) = gam_ad(i, k+1) - w2(i, k+1)*w2_ad(i, k)
        w2_ad(i, k+1) = w2_ad(i, k+1) - gam(i, k+1)*w2_ad(i, k)
      END DO
    END DO
    aa_ad = 0.0
    bet_ad = 0.0
    pp_ad = 0.0
    DO i=ie,is,-1
      CALL POPREALARRAY(w2(i, km))
      temp_ad10 = w2_ad(i, km)/bet(i)
      w1_ad(i, km) = w1_ad(i, km) + dm2(i, km)*temp_ad10
      pp_ad(i, km+1) = pp_ad(i, km+1) + dt*temp_ad10
      pp_ad(i, km) = pp_ad(i, km) - dt*temp_ad10
      ws_ad(i) = ws_ad(i) - p1(i)*temp_ad10
      w2_ad(i, km-1) = w2_ad(i, km-1) - aa(i, km)*temp_ad10
      bet_ad(i) = bet_ad(i) - (dm2(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i&
&       , km))-p1(i)*ws(i)-aa(i, km)*w2(i, km-1))*temp_ad10/bet(i)
      dm2_ad(i, km) = dm2_ad(i, km) + bet_ad(i) + w1(i, km)*temp_ad10
      p1_ad(i) = p1_ad(i) - bet_ad(i) - ws(i)*temp_ad10
      w2_ad(i, km) = 0.0
      CALL POPREALARRAY(bet(i))
      gam_ad(i, km) = gam_ad(i, km) - aa(i, km)*bet_ad(i)
      temp_ad11 = gam_ad(i, km)/bet(i)
      aa_ad(i, km) = aa_ad(i, km) + ((-1.0)-gam(i, km))*bet_ad(i) + &
&       temp_ad11 - w2(i, km-1)*temp_ad10
      bet_ad(i) = -(aa(i, km)*temp_ad11/bet(i))
      CALL POPREALARRAY(gam(i, km))
      gam_ad(i, km) = 0.0
      temp_ad12 = t1g*p1_ad(i)/dz2(i, km)
      pem_ad(i, km+1) = pem_ad(i, km+1) + temp_ad12
      pp_ad(i, km+1) = pp_ad(i, km+1) + temp_ad12
      dz2_ad(i, km) = dz2_ad(i, km) - (pem(i, km+1)+pp(i, km+1))*&
&       temp_ad12/dz2(i, km)
      p1_ad(i) = 0.0
    END DO
    DO k=km-1,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(w2(i, k))
        temp_ad8 = w2_ad(i, k)/bet(i)
        w1_ad(i, k) = w1_ad(i, k) + dm2(i, k)*temp_ad8
        pp_ad(i, k+1) = pp_ad(i, k+1) + dt*temp_ad8
        pp_ad(i, k) = pp_ad(i, k) - dt*temp_ad8
        w2_ad(i, k-1) = w2_ad(i, k-1) - aa(i, k)*temp_ad8
        bet_ad(i) = bet_ad(i) - (dm2(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i&
&         , k))-aa(i, k)*w2(i, k-1))*temp_ad8/bet(i)
        dm2_ad(i, k) = dm2_ad(i, k) + bet_ad(i) + w1(i, k)*temp_ad8
        aa_ad(i, k) = aa_ad(i, k) + ((-1.0)-gam(i, k))*bet_ad(i) - w2(i&
&         , k-1)*temp_ad8
        w2_ad(i, k) = 0.0
        CALL POPREALARRAY(bet(i))
        aa_ad(i, k+1) = aa_ad(i, k+1) - bet_ad(i)
        gam_ad(i, k) = gam_ad(i, k) - aa(i, k)*bet_ad(i)
        CALL POPREALARRAY(gam(i, k))
        temp_ad9 = gam_ad(i, k)/bet(i)
        bet_ad(i) = -(aa(i, k)*temp_ad9/bet(i))
        aa_ad(i, k) = aa_ad(i, k) + temp_ad9
        gam_ad(i, k) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(w2(i, 1))
      temp_ad7 = w2_ad(i, 1)/bet(i)
      w1_ad(i, 1) = w1_ad(i, 1) + dm2(i, 1)*temp_ad7
      pp_ad(i, 2) = pp_ad(i, 2) + dt*temp_ad7
      bet_ad(i) = bet_ad(i) - (dm2(i, 1)*w1(i, 1)+dt*pp(i, 2))*temp_ad7/&
&       bet(i)
      dm2_ad(i, 1) = dm2_ad(i, 1) + bet_ad(i) + w1(i, 1)*temp_ad7
      w2_ad(i, 1) = 0.0
      CALL POPREALARRAY(bet(i))
      aa_ad(i, 2) = aa_ad(i, 2) - bet_ad(i)
      bet_ad(i) = 0.0
    END DO
    DO k=km,2,-1
      DO i=ie,is,-1
        temp2 = dz2(i, k-1) + dz2(i, k)
        temp_ad5 = t1g*aa_ad(i, k)/temp2
        temp_ad6 = -((pem(i, k)+pp(i, k))*temp_ad5/temp2)
        pem_ad(i, k) = pem_ad(i, k) + temp_ad5
        pp_ad(i, k) = pp_ad(i, k) + temp_ad5
        dz2_ad(i, k-1) = dz2_ad(i, k-1) + temp_ad6
        dz2_ad(i, k) = dz2_ad(i, k) + temp_ad6
        aa_ad(i, k) = 0.0
      END DO
    END DO
    DO k=2,km,1
      DO i=ie,is,-1
        CALL POPREALARRAY(pp(i, k))
        gam_ad(i, k) = gam_ad(i, k) - pp(i, k+1)*pp_ad(i, k)
        pp_ad(i, k+1) = pp_ad(i, k+1) - gam(i, k)*pp_ad(i, k)
      END DO
    END DO
    dd_ad = 0.0
    DO k=km,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pp(i, k+1))
        temp_ad3 = pp_ad(i, k+1)/bet(i)
        dd_ad(i, k) = dd_ad(i, k) + temp_ad3
        pp_ad(i, k) = pp_ad(i, k) - temp_ad3
        bet_ad(i) = bet_ad(i) - (dd(i, k)-pp(i, k))*temp_ad3/bet(i)
        pp_ad(i, k+1) = 0.0
        CALL POPREALARRAY(bet(i))
        bb_ad(i, k) = bb_ad(i, k) + bet_ad(i)
        gam_ad(i, k) = gam_ad(i, k) - bet_ad(i)
        temp_ad4 = gam_ad(i, k)/bet(i)
        bet_ad(i) = -(g_rat(i, k-1)*temp_ad4/bet(i))
        g_rat_ad(i, k-1) = g_rat_ad(i, k-1) + temp_ad4
        gam_ad(i, k) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(dd(i, km))
      pe_ad(i, km) = pe_ad(i, km) + 3.*dd_ad(i, km)
      dd_ad(i, km) = 0.0
      bb_ad(i, km) = 0.0
      temp_ad2 = pp_ad(i, 2)/bet(i)
      dd_ad(i, 1) = dd_ad(i, 1) + temp_ad2
      bet_ad(i) = bet_ad(i) - dd(i, 1)*temp_ad2/bet(i)
      pp_ad(i, 2) = 0.0
      pp_ad(i, 1) = 0.0
      bb_ad(i, 1) = bb_ad(i, 1) + bet_ad(i)
      bet_ad(i) = 0.0
    END DO
    DO k=km-1,1,-1
      DO i=ie,is,-1
        temp_ad0 = 3.*dd_ad(i, k)
        pe_ad(i, k) = pe_ad(i, k) + temp_ad0
        g_rat_ad(i, k) = g_rat_ad(i, k) + 2.*bb_ad(i, k) + pe(i, k+1)*&
&         temp_ad0
        pe_ad(i, k+1) = pe_ad(i, k+1) + g_rat(i, k)*temp_ad0
        dd_ad(i, k) = 0.0
        bb_ad(i, k) = 0.0
        temp_ad1 = g_rat_ad(i, k)/dm2(i, k+1)
        dm2_ad(i, k) = dm2_ad(i, k) + temp_ad1
        dm2_ad(i, k+1) = dm2_ad(i, k+1) - dm2(i, k)*temp_ad1/dm2(i, k+1)
        g_rat_ad(i, k) = 0.0
      END DO
    END DO
    DO k=km,1,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe(i, k))
        temp1 = dz2(i, k)
        temp0 = dm2(i, k)*pt2(i, k)
        temp = temp0/temp1
        temp_ad = gama*EXP(gama*LOG(-(rgas*temp)))*pe_ad(i, k)/(temp*&
&         temp1)
        dm2_ad(i, k) = dm2_ad(i, k) + pt2(i, k)*temp_ad
        pt2_ad(i, k) = pt2_ad(i, k) + dm2(i, k)*temp_ad
        dz2_ad(i, k) = dz2_ad(i, k) - temp*temp_ad
        pm2_ad(i, k) = pm2_ad(i, k) - pe_ad(i, k)
        pe_ad(i, k) = 0.0
        w2_ad(i, k) = w2_ad(i, k) + w1_ad(i, k)
        w1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE SIM1_SOLVER_BWD
  SUBROUTINE SIM1_SOLVER(dt, is, ie, km, rgas, gama, gm2, cp2, kappa, pe&
&   , dm2, pm2, pem, w2, dz2, pt2, ws, p_fac)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pt2, pm2, gm2, cp2
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL, INTENT(OUT) :: pe(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, g_rat, gam
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie) :: p1, bet
    REAL :: t1g, rdt, capa1
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: max1
    REAL :: max2
    t1g = gama*2.*dt*dt
    rdt = 1./dt
    capa1 = kappa - 1.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
        pe(i, k) = EXP(gama*LOG(-(dm2(i, k)/dz2(i, k)*rgas*pt2(i, k)))) &
&         - pm2(i, k)
      END DO
    END DO
    DO k=1,km-1
      DO i=is,ie
        g_rat(i, k) = dm2(i, k)/dm2(i, k+1)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd(i, k) = 3.*(pe(i, k)+g_rat(i, k)*pe(i, k+1))
      END DO
    END DO
    DO i=is,ie
      bet(i) = bb(i, 1)
      pp(i, 1) = 0.
      pp(i, 2) = dd(i, 1)/bet(i)
      bb(i, km) = 2.
      dd(i, km) = 3.*pe(i, km)
    END DO
    DO k=2,km
      DO i=is,ie
        gam(i, k) = g_rat(i, k-1)/bet(i)
        bet(i) = bb(i, k) - gam(i, k)
        pp(i, k+1) = (dd(i, k)-pp(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        pp(i, k) = pp(i, k) - gam(i, k)*pp(i, k+1)
      END DO
    END DO
! Start the w-solver
    DO k=2,km
      DO i=is,ie
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*(pem(i, k)+pp(i, k))
      END DO
    END DO
    DO i=is,ie
      bet(i) = dm2(i, 1) - aa(i, 2)
      w2(i, 1) = (dm2(i, 1)*w1(i, 1)+dt*pp(i, 2))/bet(i)
    END DO
    DO k=2,km-1
      DO i=is,ie
        gam(i, k) = aa(i, k)/bet(i)
        bet(i) = dm2(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        w2(i, k) = (dm2(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))-aa(i, k)&
&         *w2(i, k-1))/bet(i)
      END DO
    END DO
    DO i=is,ie
      p1(i) = t1g/dz2(i, km)*(pem(i, km+1)+pp(i, km+1))
      gam(i, km) = aa(i, km)/bet(i)
      bet(i) = dm2(i, km) - (aa(i, km)+p1(i)+aa(i, km)*gam(i, km))
      w2(i, km) = (dm2(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-p1(i)&
&       *ws(i)-aa(i, km)*w2(i, km-1))/bet(i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
    DO i=is,ie
      pe(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        pe(i, k+1) = pe(i, k) + dm2(i, k)*(w2(i, k)-w1(i, k))*rdt
      END DO
    END DO
    DO i=is,ie
      p1(i) = (pe(i, km)+2.*pe(i, km+1))*r3
      IF (p_fac*pm2(i, km) .LT. p1(i) + pm2(i, km)) THEN
        max1 = p1(i) + pm2(i, km)
      ELSE
        max1 = p_fac*pm2(i, km)
      END IF
      dz2(i, km) = -(dm2(i, km)*rgas*pt2(i, km)*EXP(capa1*LOG(max1)))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1(i) = (pe(i, k)+bb(i, k)*pe(i, k+1)+g_rat(i, k)*pe(i, k+2))*r3&
&         - g_rat(i, k)*p1(i)
        IF (p_fac*pm2(i, k) .LT. p1(i) + pm2(i, k)) THEN
          max2 = p1(i) + pm2(i, k)
        ELSE
          max2 = p_fac*pm2(i, k)
        END IF
        dz2(i, k) = -(dm2(i, k)*rgas*pt2(i, k)*EXP(capa1*LOG(max2)))
      END DO
    END DO
  END SUBROUTINE SIM1_SOLVER
!  Differentiation of sim_solver in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_m
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
!   gradient     of useful results: ws pe2 dm2 dz2 w2 pm2 pem pt2
!   with respect to varying inputs: ws pe2 dm2 dz2 w2 pm2 pem pt2
  SUBROUTINE SIM_SOLVER_FWD(dt, is, ie, km, rgas, gama, gm2, cp2, kappa&
&   , pe2, dm2, pm2, pem, w2, dz2, pt2, ws, alpha, p_fac, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac, alpha, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pt2, pm2, gm2, cp2
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL :: pe2(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, wk, g_rat, gam
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL :: beta, t2, t1g, rdt, ra, capa1
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: max1
    REAL :: max2

    aa = 0.0
    bb = 0.0
    dd = 0.0
    w1 = 0.0
    wk = 0.0
    g_rat = 0.0
    gam = 0.0
    pp = 0.0
    p1 = 0.0
    wk1 = 0.0
    bet = 0.0
    t1g = 0.0
    rdt = 0.0
    capa1 = 0.0
    max1 = 0.0
    max2 = 0.0

    beta = 1. - alpha
    ra = 1./alpha
    t2 = beta/alpha
    t1g = 2.*gama*(alpha*dt)**2
    rdt = 1./dt
    capa1 = kappa - 1.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
! P_g perturbation
        CALL PUSHREALARRAY(pe2(i, k))
        pe2(i, k) = EXP(gama*LOG(-(dm2(i, k)/dz2(i, k)*rgas*pt2(i, k))))&
&         - pm2(i, k)
      END DO
    END DO
    DO k=1,km-1
      DO i=is,ie
        g_rat(i, k) = dm2(i, k)/dm2(i, k+1)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd(i, k) = 3.*(pe2(i, k)+g_rat(i, k)*pe2(i, k+1))
      END DO
    END DO
    DO i=is,ie
      bet(i) = bb(i, 1)
      pp(i, 1) = 0.
      pp(i, 2) = dd(i, 1)/bet(i)
      bb(i, km) = 2.
      CALL PUSHREALARRAY(dd(i, km))
      dd(i, km) = 3.*pe2(i, km)
    END DO
    DO k=2,km
      DO i=is,ie
        gam(i, k) = g_rat(i, k-1)/bet(i)
        CALL PUSHREALARRAY(bet(i))
        bet(i) = bb(i, k) - gam(i, k)
        CALL PUSHREALARRAY(pp(i, k+1))
        pp(i, k+1) = (dd(i, k)-pp(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        CALL PUSHREALARRAY(pp(i, k))
        pp(i, k) = pp(i, k) - gam(i, k)*pp(i, k+1)
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
! pe2 is Full p
        CALL PUSHREALARRAY(pe2(i, k))
        pe2(i, k) = pem(i, k) + pp(i, k)
      END DO
    END DO
    DO k=2,km
      DO i=is,ie
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*pe2(i, k)
        wk(i, k) = t2*aa(i, k)*(w1(i, k-1)-w1(i, k))
        CALL PUSHREALARRAY(aa(i, k))
        aa(i, k) = aa(i, k) - scale_m*dm2(i, 1)
      END DO
    END DO
! Top:
    DO i=is,ie
      CALL PUSHREALARRAY(bet(i))
      bet(i) = dm2(i, 1) - aa(i, 2)
      CALL PUSHREALARRAY(w2(i, 1))
      w2(i, 1) = (dm2(i, 1)*w1(i, 1)+dt*pp(i, 2)+wk(i, 2))/bet(i)
    END DO
! Interior:
    DO k=2,km-1
      DO i=is,ie
        CALL PUSHREALARRAY(gam(i, k))
        gam(i, k) = aa(i, k)/bet(i)
        CALL PUSHREALARRAY(bet(i))
        bet(i) = dm2(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        CALL PUSHREALARRAY(w2(i, k))
        w2(i, k) = (dm2(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))+wk(i, k+&
&         1)-wk(i, k)-aa(i, k)*w2(i, k-1))/bet(i)
      END DO
    END DO
! Bottom: k=km
    DO i=is,ie
      wk1(i) = t1g/dz2(i, km)*pe2(i, km+1)
      CALL PUSHREALARRAY(gam(i, km))
      gam(i, km) = aa(i, km)/bet(i)
      CALL PUSHREALARRAY(bet(i))
      bet(i) = dm2(i, km) - (aa(i, km)+wk1(i)+aa(i, km)*gam(i, km))
      CALL PUSHREALARRAY(w2(i, km))
      w2(i, km) = (dm2(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk(i&
&       , km)+wk1(i)*(t2*w1(i, km)-ra*ws(i))-aa(i, km)*w2(i, km-1))/bet(&
&       i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        CALL PUSHREALARRAY(w2(i, k))
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
    DO i=is,ie
      CALL PUSHREALARRAY(pe2(i, 1))
      pe2(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        CALL PUSHREALARRAY(pe2(i, k+1))
        pe2(i, k+1) = pe2(i, k) + (dm2(i, k)*(w2(i, k)-w1(i, k))*rdt-&
&         beta*(pp(i, k+1)-pp(i, k)))*ra
      END DO
    END DO
    DO i=is,ie
      p1(i) = (pe2(i, km)+2.*pe2(i, km+1))*r3
      IF (p_fac*pm2(i, km) .LT. p1(i) + pm2(i, km)) THEN
        CALL PUSHREALARRAY(max1)
        max1 = p1(i) + pm2(i, km)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHREALARRAY(max1)
        max1 = p_fac*pm2(i, km)
        CALL PUSHCONTROL(1,1)
      END IF
      CALL PUSHREALARRAY(dz2(i, km))
      dz2(i, km) = -(dm2(i, km)*rgas*pt2(i, km)*EXP(capa1*LOG(max1)))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        CALL PUSHREALARRAY(p1(i))
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        IF (p_fac*pm2(i, k) .LT. p1(i) + pm2(i, k)) THEN
          CALL PUSHREALARRAY(max2)
          max2 = p1(i) + pm2(i, k)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHREALARRAY(max2)
          max2 = p_fac*pm2(i, k)
          CALL PUSHCONTROL(1,1)
        END IF
! delz = -dm*R*T_m / p_gas
        CALL PUSHREALARRAY(dz2(i, k))
        dz2(i, k) = -(dm2(i, k)*rgas*pt2(i, k)*EXP(capa1*LOG(max2)))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        CALL PUSHREALARRAY(pe2(i, k))
        pe2(i, k) = pe2(i, k) + beta*(pp(i, k)-pe2(i, k))
      END DO
    END DO
    CALL PUSHREALARRAY(gam, (ie-is+1)*km)
    CALL PUSHREALARRAY(t1g)
    CALL PUSHREALARRAY(wk, (ie-is+1)*km)
    CALL PUSHREALARRAY(capa1)
    CALL PUSHREALARRAY(g_rat, (ie-is+1)*km)
    CALL PUSHREALARRAY(pp, (ie-is+1)*(km+1))
    CALL PUSHREALARRAY(beta)
    CALL PUSHREALARRAY(rdt)
    CALL PUSHREALARRAY(t2)
    CALL PUSHREALARRAY(w1, (ie-is+1)*km)
    CALL PUSHREALARRAY(bb, (ie-is+1)*km)
    CALL PUSHREALARRAY(ra)
    CALL PUSHREALARRAY(wk1, ie - is + 1)
    CALL PUSHREALARRAY(bet, ie - is + 1)
    CALL PUSHREALARRAY(max2)
    CALL PUSHREALARRAY(max1)
    CALL PUSHREALARRAY(aa, (ie-is+1)*km)
    CALL PUSHREALARRAY(dd, (ie-is+1)*km)
  END SUBROUTINE SIM_SOLVER_FWD
!  Differentiation of sim_solver in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: ws pe2 dm2 dz2 w2 pm2 pem pt2
!   with respect to varying inputs: ws pe2 dm2 dz2 w2 pm2 pem pt2
  SUBROUTINE SIM_SOLVER_BWD(dt, is, ie, km, rgas, gama, gm2, cp2, kappa&
&   , pe2, pe2_ad, dm2, dm2_ad, pm2, pm2_ad, pem, pem_ad, w2, w2_ad, dz2&
&   , dz2_ad, pt2, pt2_ad, ws, ws_ad, alpha, p_fac, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac, alpha, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pt2, pm2, gm2, cp2
    REAL, DIMENSION(is:ie, km) :: dm2_ad, pt2_ad, pm2_ad
    REAL, INTENT(IN) :: ws(is:ie)
    REAL :: ws_ad(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL, DIMENSION(is:ie, km+1) :: pem_ad
    REAL :: pe2(is:ie, km+1)
    REAL :: pe2_ad(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2_ad, w2_ad
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, wk, g_rat, gam
    REAL, DIMENSION(is:ie, km) :: aa_ad, bb_ad, dd_ad, w1_ad, wk_ad, &
&   g_rat_ad, gam_ad
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie, km+1) :: pp_ad
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL, DIMENSION(is:ie) :: p1_ad, wk1_ad, bet_ad
    REAL :: beta, t2, t1g, rdt, ra, capa1
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: max1
    REAL :: max1_ad
    REAL :: max2
    REAL :: max2_ad
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp2
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp3
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: temp4
    REAL :: temp_ad16
    REAL :: temp_ad17
    REAL :: temp5
    REAL :: temp_ad18
    INTEGER :: branch

    aa = 0.0
    bb = 0.0
    dd = 0.0
    w1 = 0.0
    wk = 0.0
    g_rat = 0.0
    gam = 0.0
    pp = 0.0
    p1 = 0.0
    wk1 = 0.0
    bet = 0.0
    t1g = 0.0
    rdt = 0.0
    capa1 = 0.0
    max1 = 0.0
    max2 = 0.0
    branch = 0

    CALL POPREALARRAY(dd, (ie-is+1)*km)
    CALL POPREALARRAY(aa, (ie-is+1)*km)
    CALL POPREALARRAY(max1)
    CALL POPREALARRAY(max2)
    CALL POPREALARRAY(bet, ie - is + 1)
    CALL POPREALARRAY(wk1, ie - is + 1)
    CALL POPREALARRAY(ra)
    CALL POPREALARRAY(bb, (ie-is+1)*km)
    CALL POPREALARRAY(w1, (ie-is+1)*km)
    CALL POPREALARRAY(t2)
    CALL POPREALARRAY(rdt)
    CALL POPREALARRAY(beta)
    CALL POPREALARRAY(pp, (ie-is+1)*(km+1))
    CALL POPREALARRAY(g_rat, (ie-is+1)*km)
    CALL POPREALARRAY(capa1)
    CALL POPREALARRAY(wk, (ie-is+1)*km)
    CALL POPREALARRAY(t1g)
    CALL POPREALARRAY(gam, (ie-is+1)*km)
    pp_ad = 0.0
    DO k=km+1,1,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k))
        pp_ad(i, k) = pp_ad(i, k) + beta*pe2_ad(i, k)
        pe2_ad(i, k) = (1.0-beta)*pe2_ad(i, k)
      END DO
    END DO
    p1_ad = 0.0
    bb_ad = 0.0
    g_rat_ad = 0.0
    DO k=1,km-1,1
      DO i=ie,is,-1
        CALL POPREALARRAY(dz2(i, k))
        temp5 = capa1*LOG(max2)
        temp_ad18 = -(rgas*EXP(temp5)*dz2_ad(i, k))
        dm2_ad(i, k) = dm2_ad(i, k) + pt2(i, k)*temp_ad18
        pt2_ad(i, k) = pt2_ad(i, k) + dm2(i, k)*temp_ad18
        max2_ad = -(capa1*EXP(temp5)*dm2(i, k)*pt2(i, k)*rgas*dz2_ad(i, &
&         k)/max2)
        dz2_ad(i, k) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(max2)
          p1_ad(i) = p1_ad(i) + max2_ad
          pm2_ad(i, k) = pm2_ad(i, k) + max2_ad
        ELSE
          CALL POPREALARRAY(max2)
          pm2_ad(i, k) = pm2_ad(i, k) + p_fac*max2_ad
        END IF
        CALL POPREALARRAY(p1(i))
        temp_ad17 = r3*p1_ad(i)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad17
        bb_ad(i, k) = bb_ad(i, k) + pe2(i, k+1)*temp_ad17
        pe2_ad(i, k+1) = pe2_ad(i, k+1) + bb(i, k)*temp_ad17
        g_rat_ad(i, k) = g_rat_ad(i, k) + pe2(i, k+2)*temp_ad17 - p1(i)*&
&         p1_ad(i)
        pe2_ad(i, k+2) = pe2_ad(i, k+2) + g_rat(i, k)*temp_ad17
        p1_ad(i) = -(g_rat(i, k)*p1_ad(i))
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(dz2(i, km))
      temp4 = capa1*LOG(max1)
      temp_ad16 = -(rgas*EXP(temp4)*dz2_ad(i, km))
      dm2_ad(i, km) = dm2_ad(i, km) + pt2(i, km)*temp_ad16
      pt2_ad(i, km) = pt2_ad(i, km) + dm2(i, km)*temp_ad16
      max1_ad = -(capa1*EXP(temp4)*dm2(i, km)*pt2(i, km)*rgas*dz2_ad(i, &
&       km)/max1)
      dz2_ad(i, km) = 0.0
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(max1)
        p1_ad(i) = p1_ad(i) + max1_ad
        pm2_ad(i, km) = pm2_ad(i, km) + max1_ad
      ELSE
        CALL POPREALARRAY(max1)
        pm2_ad(i, km) = pm2_ad(i, km) + p_fac*max1_ad
      END IF
      pe2_ad(i, km) = pe2_ad(i, km) + r3*p1_ad(i)
      pe2_ad(i, km+1) = pe2_ad(i, km+1) + r3*2.*p1_ad(i)
      p1_ad(i) = 0.0
    END DO
    w1_ad = 0.0
    DO k=km,1,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k+1))
        temp_ad14 = ra*pe2_ad(i, k+1)
        temp_ad15 = rdt*dm2(i, k)*temp_ad14
        pe2_ad(i, k) = pe2_ad(i, k) + pe2_ad(i, k+1)
        dm2_ad(i, k) = dm2_ad(i, k) + rdt*(w2(i, k)-w1(i, k))*temp_ad14
        w2_ad(i, k) = w2_ad(i, k) + temp_ad15
        w1_ad(i, k) = w1_ad(i, k) - temp_ad15
        pp_ad(i, k+1) = pp_ad(i, k+1) - beta*temp_ad14
        pp_ad(i, k) = pp_ad(i, k) + beta*temp_ad14
        pe2_ad(i, k+1) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(pe2(i, 1))
      pe2_ad(i, 1) = 0.0
    END DO
    gam_ad = 0.0
    DO k=1,km-1,1
      DO i=ie,is,-1
        CALL POPREALARRAY(w2(i, k))
        gam_ad(i, k+1) = gam_ad(i, k+1) - w2(i, k+1)*w2_ad(i, k)
        w2_ad(i, k+1) = w2_ad(i, k+1) - gam(i, k+1)*w2_ad(i, k)
      END DO
    END DO
    aa_ad = 0.0
    bet_ad = 0.0
    wk1_ad = 0.0
    wk_ad = 0.0
    DO i=ie,is,-1
      CALL POPREALARRAY(w2(i, km))
      temp_ad11 = w2_ad(i, km)/bet(i)
      temp3 = t2*w1(i, km) - ra*ws(i)
      w1_ad(i, km) = w1_ad(i, km) + (wk1(i)*t2+dm2(i, km))*temp_ad11
      pp_ad(i, km+1) = pp_ad(i, km+1) + dt*temp_ad11
      pp_ad(i, km) = pp_ad(i, km) - dt*temp_ad11
      wk_ad(i, km) = wk_ad(i, km) - temp_ad11
      ws_ad(i) = ws_ad(i) - wk1(i)*ra*temp_ad11
      w2_ad(i, km-1) = w2_ad(i, km-1) - aa(i, km)*temp_ad11
      bet_ad(i) = bet_ad(i) - (dm2(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i&
&       , km))-wk(i, km)+wk1(i)*temp3-aa(i, km)*w2(i, km-1))*temp_ad11/&
&       bet(i)
      dm2_ad(i, km) = dm2_ad(i, km) + bet_ad(i) + w1(i, km)*temp_ad11
      wk1_ad(i) = wk1_ad(i) + temp3*temp_ad11 - bet_ad(i)
      w2_ad(i, km) = 0.0
      CALL POPREALARRAY(bet(i))
      gam_ad(i, km) = gam_ad(i, km) - aa(i, km)*bet_ad(i)
      temp_ad12 = gam_ad(i, km)/bet(i)
      aa_ad(i, km) = aa_ad(i, km) + ((-1.0)-gam(i, km))*bet_ad(i) + &
&       temp_ad12 - w2(i, km-1)*temp_ad11
      bet_ad(i) = -(aa(i, km)*temp_ad12/bet(i))
      CALL POPREALARRAY(gam(i, km))
      gam_ad(i, km) = 0.0
      temp_ad13 = t1g*wk1_ad(i)/dz2(i, km)
      pe2_ad(i, km+1) = pe2_ad(i, km+1) + temp_ad13
      dz2_ad(i, km) = dz2_ad(i, km) - pe2(i, km+1)*temp_ad13/dz2(i, km)
      wk1_ad(i) = 0.0
    END DO
    DO k=km-1,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(w2(i, k))
        temp_ad9 = w2_ad(i, k)/bet(i)
        w1_ad(i, k) = w1_ad(i, k) + dm2(i, k)*temp_ad9
        pp_ad(i, k+1) = pp_ad(i, k+1) + dt*temp_ad9
        pp_ad(i, k) = pp_ad(i, k) - dt*temp_ad9
        wk_ad(i, k+1) = wk_ad(i, k+1) + temp_ad9
        wk_ad(i, k) = wk_ad(i, k) - temp_ad9
        w2_ad(i, k-1) = w2_ad(i, k-1) - aa(i, k)*temp_ad9
        bet_ad(i) = bet_ad(i) - (dm2(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i&
&         , k))+wk(i, k+1)-wk(i, k)-aa(i, k)*w2(i, k-1))*temp_ad9/bet(i)
        dm2_ad(i, k) = dm2_ad(i, k) + bet_ad(i) + w1(i, k)*temp_ad9
        aa_ad(i, k) = aa_ad(i, k) + ((-1.0)-gam(i, k))*bet_ad(i) - w2(i&
&         , k-1)*temp_ad9
        w2_ad(i, k) = 0.0
        CALL POPREALARRAY(bet(i))
        aa_ad(i, k+1) = aa_ad(i, k+1) - bet_ad(i)
        gam_ad(i, k) = gam_ad(i, k) - aa(i, k)*bet_ad(i)
        CALL POPREALARRAY(gam(i, k))
        temp_ad10 = gam_ad(i, k)/bet(i)
        bet_ad(i) = -(aa(i, k)*temp_ad10/bet(i))
        aa_ad(i, k) = aa_ad(i, k) + temp_ad10
        gam_ad(i, k) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(w2(i, 1))
      temp_ad8 = w2_ad(i, 1)/bet(i)
      w1_ad(i, 1) = w1_ad(i, 1) + dm2(i, 1)*temp_ad8
      pp_ad(i, 2) = pp_ad(i, 2) + dt*temp_ad8
      wk_ad(i, 2) = wk_ad(i, 2) + temp_ad8
      bet_ad(i) = bet_ad(i) - (dm2(i, 1)*w1(i, 1)+dt*pp(i, 2)+wk(i, 2))*&
&       temp_ad8/bet(i)
      dm2_ad(i, 1) = dm2_ad(i, 1) + bet_ad(i) + w1(i, 1)*temp_ad8
      w2_ad(i, 1) = 0.0
      CALL POPREALARRAY(bet(i))
      aa_ad(i, 2) = aa_ad(i, 2) - bet_ad(i)
      bet_ad(i) = 0.0
    END DO
    DO k=km,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(aa(i, k))
        dm2_ad(i, 1) = dm2_ad(i, 1) - scale_m*aa_ad(i, k)
        temp_ad5 = t2*aa(i, k)*wk_ad(i, k)
        aa_ad(i, k) = aa_ad(i, k) + t2*(w1(i, k-1)-w1(i, k))*wk_ad(i, k)
        w1_ad(i, k-1) = w1_ad(i, k-1) + temp_ad5
        w1_ad(i, k) = w1_ad(i, k) - temp_ad5
        wk_ad(i, k) = 0.0
        temp2 = dz2(i, k-1) + dz2(i, k)
        temp_ad6 = t1g*aa_ad(i, k)/temp2
        temp_ad7 = -(pe2(i, k)*temp_ad6/temp2)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad6
        dz2_ad(i, k-1) = dz2_ad(i, k-1) + temp_ad7
        dz2_ad(i, k) = dz2_ad(i, k) + temp_ad7
        aa_ad(i, k) = 0.0
      END DO
    END DO
    DO k=km+1,1,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k))
        pem_ad(i, k) = pem_ad(i, k) + pe2_ad(i, k)
        pp_ad(i, k) = pp_ad(i, k) + pe2_ad(i, k)
        pe2_ad(i, k) = 0.0
      END DO
    END DO
    DO k=2,km,1
      DO i=ie,is,-1
        CALL POPREALARRAY(pp(i, k))
        gam_ad(i, k) = gam_ad(i, k) - pp(i, k+1)*pp_ad(i, k)
        pp_ad(i, k+1) = pp_ad(i, k+1) - gam(i, k)*pp_ad(i, k)
      END DO
    END DO
    dd_ad = 0.0
    DO k=km,2,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pp(i, k+1))
        temp_ad3 = pp_ad(i, k+1)/bet(i)
        dd_ad(i, k) = dd_ad(i, k) + temp_ad3
        pp_ad(i, k) = pp_ad(i, k) - temp_ad3
        bet_ad(i) = bet_ad(i) - (dd(i, k)-pp(i, k))*temp_ad3/bet(i)
        pp_ad(i, k+1) = 0.0
        CALL POPREALARRAY(bet(i))
        bb_ad(i, k) = bb_ad(i, k) + bet_ad(i)
        gam_ad(i, k) = gam_ad(i, k) - bet_ad(i)
        temp_ad4 = gam_ad(i, k)/bet(i)
        bet_ad(i) = -(g_rat(i, k-1)*temp_ad4/bet(i))
        g_rat_ad(i, k-1) = g_rat_ad(i, k-1) + temp_ad4
        gam_ad(i, k) = 0.0
      END DO
    END DO
    DO i=ie,is,-1
      CALL POPREALARRAY(dd(i, km))
      pe2_ad(i, km) = pe2_ad(i, km) + 3.*dd_ad(i, km)
      dd_ad(i, km) = 0.0
      bb_ad(i, km) = 0.0
      temp_ad2 = pp_ad(i, 2)/bet(i)
      dd_ad(i, 1) = dd_ad(i, 1) + temp_ad2
      bet_ad(i) = bet_ad(i) - dd(i, 1)*temp_ad2/bet(i)
      pp_ad(i, 2) = 0.0
      pp_ad(i, 1) = 0.0
      bb_ad(i, 1) = bb_ad(i, 1) + bet_ad(i)
      bet_ad(i) = 0.0
    END DO
    DO k=km-1,1,-1
      DO i=ie,is,-1
        temp_ad0 = 3.*dd_ad(i, k)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad0
        g_rat_ad(i, k) = g_rat_ad(i, k) + 2.*bb_ad(i, k) + pe2(i, k+1)*&
&         temp_ad0
        pe2_ad(i, k+1) = pe2_ad(i, k+1) + g_rat(i, k)*temp_ad0
        dd_ad(i, k) = 0.0
        bb_ad(i, k) = 0.0
        temp_ad1 = g_rat_ad(i, k)/dm2(i, k+1)
        dm2_ad(i, k) = dm2_ad(i, k) + temp_ad1
        dm2_ad(i, k+1) = dm2_ad(i, k+1) - dm2(i, k)*temp_ad1/dm2(i, k+1)
        g_rat_ad(i, k) = 0.0
      END DO
    END DO
    DO k=km,1,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(pe2(i, k))
        temp1 = dz2(i, k)
        temp0 = dm2(i, k)*pt2(i, k)
        temp = temp0/temp1
        temp_ad = gama*EXP(gama*LOG(-(rgas*temp)))*pe2_ad(i, k)/(temp*&
&         temp1)
        dm2_ad(i, k) = dm2_ad(i, k) + pt2(i, k)*temp_ad
        pt2_ad(i, k) = pt2_ad(i, k) + dm2(i, k)*temp_ad
        dz2_ad(i, k) = dz2_ad(i, k) - temp*temp_ad
        pm2_ad(i, k) = pm2_ad(i, k) - pe2_ad(i, k)
        pe2_ad(i, k) = 0.0
        w2_ad(i, k) = w2_ad(i, k) + w1_ad(i, k)
        w1_ad(i, k) = 0.0
      END DO
    END DO
  END SUBROUTINE SIM_SOLVER_BWD
  SUBROUTINE SIM_SOLVER(dt, is, ie, km, rgas, gama, gm2, cp2, kappa, pe2&
&   , dm2, pm2, pem, w2, dz2, pt2, ws, alpha, p_fac, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac, alpha, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pt2, pm2, gm2, cp2
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL, INTENT(OUT) :: pe2(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, wk, g_rat, gam
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL :: beta, t2, t1g, rdt, ra, capa1
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: max1
    REAL :: max2
    beta = 1. - alpha
    ra = 1./alpha
    t2 = beta/alpha
    t1g = 2.*gama*(alpha*dt)**2
    rdt = 1./dt
    capa1 = kappa - 1.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
! P_g perturbation
        pe2(i, k) = EXP(gama*LOG(-(dm2(i, k)/dz2(i, k)*rgas*pt2(i, k))))&
&         - pm2(i, k)
      END DO
    END DO
    DO k=1,km-1
      DO i=is,ie
        g_rat(i, k) = dm2(i, k)/dm2(i, k+1)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd(i, k) = 3.*(pe2(i, k)+g_rat(i, k)*pe2(i, k+1))
      END DO
    END DO
    DO i=is,ie
      bet(i) = bb(i, 1)
      pp(i, 1) = 0.
      pp(i, 2) = dd(i, 1)/bet(i)
      bb(i, km) = 2.
      dd(i, km) = 3.*pe2(i, km)
    END DO
    DO k=2,km
      DO i=is,ie
        gam(i, k) = g_rat(i, k-1)/bet(i)
        bet(i) = bb(i, k) - gam(i, k)
        pp(i, k+1) = (dd(i, k)-pp(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        pp(i, k) = pp(i, k) - gam(i, k)*pp(i, k+1)
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
! pe2 is Full p
        pe2(i, k) = pem(i, k) + pp(i, k)
      END DO
    END DO
    DO k=2,km
      DO i=is,ie
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*pe2(i, k)
        wk(i, k) = t2*aa(i, k)*(w1(i, k-1)-w1(i, k))
        aa(i, k) = aa(i, k) - scale_m*dm2(i, 1)
      END DO
    END DO
! Top:
    DO i=is,ie
      bet(i) = dm2(i, 1) - aa(i, 2)
      w2(i, 1) = (dm2(i, 1)*w1(i, 1)+dt*pp(i, 2)+wk(i, 2))/bet(i)
    END DO
! Interior:
    DO k=2,km-1
      DO i=is,ie
        gam(i, k) = aa(i, k)/bet(i)
        bet(i) = dm2(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        w2(i, k) = (dm2(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))+wk(i, k+&
&         1)-wk(i, k)-aa(i, k)*w2(i, k-1))/bet(i)
      END DO
    END DO
! Bottom: k=km
    DO i=is,ie
      wk1(i) = t1g/dz2(i, km)*pe2(i, km+1)
      gam(i, km) = aa(i, km)/bet(i)
      bet(i) = dm2(i, km) - (aa(i, km)+wk1(i)+aa(i, km)*gam(i, km))
      w2(i, km) = (dm2(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk(i&
&       , km)+wk1(i)*(t2*w1(i, km)-ra*ws(i))-aa(i, km)*w2(i, km-1))/bet(&
&       i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
    DO i=is,ie
      pe2(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        pe2(i, k+1) = pe2(i, k) + (dm2(i, k)*(w2(i, k)-w1(i, k))*rdt-&
&         beta*(pp(i, k+1)-pp(i, k)))*ra
      END DO
    END DO
    DO i=is,ie
      p1(i) = (pe2(i, km)+2.*pe2(i, km+1))*r3
      IF (p_fac*pm2(i, km) .LT. p1(i) + pm2(i, km)) THEN
        max1 = p1(i) + pm2(i, km)
      ELSE
        max1 = p_fac*pm2(i, km)
      END IF
      dz2(i, km) = -(dm2(i, km)*rgas*pt2(i, km)*EXP(capa1*LOG(max1)))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        IF (p_fac*pm2(i, k) .LT. p1(i) + pm2(i, k)) THEN
          max2 = p1(i) + pm2(i, k)
        ELSE
          max2 = p_fac*pm2(i, k)
        END IF
! delz = -dm*R*T_m / p_gas
        dz2(i, k) = -(dm2(i, k)*rgas*pt2(i, k)*EXP(capa1*LOG(max2)))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        pe2(i, k) = pe2(i, k) + beta*(pp(i, k)-pe2(i, k))
      END DO
    END DO
  END SUBROUTINE SIM_SOLVER
  SUBROUTINE EDGE_SCALAR(q1, qe, i1, i2, km, id)
    IMPLICIT NONE
! Optimized for wind profile reconstruction:
    INTEGER, INTENT(IN) :: i1, i2, km
! 0: pp 1: wind
    INTEGER, INTENT(IN) :: id
    REAL, DIMENSION(i1:i2, km), INTENT(IN) :: q1
    REAL, DIMENSION(i1:i2, km+1), INTENT(OUT) :: qe
!-----------------------------------------------------------------------
    REAL, PARAMETER :: r2o3=2./3.
    REAL, PARAMETER :: r4o3=4./3.
    REAL :: gak(km)
    REAL :: bet
    INTEGER :: i, k
!------------------------------------------------
! Optimized coding for uniform grid: SJL Apr 2007
!------------------------------------------------
    IF (id .EQ. 1) THEN
      DO i=i1,i2
        qe(i, 1) = r4o3*q1(i, 1) + r2o3*q1(i, 2)
      END DO
    ELSE
      DO i=i1,i2
        qe(i, 1) = 1.e30
      END DO
    END IF
    gak(1) = 7./3.
    DO k=2,km
      gak(k) = 1./(4.-gak(k-1))
      DO i=i1,i2
        qe(i, k) = (3.*(q1(i, k-1)+q1(i, k))-qe(i, k-1))*gak(k)
      END DO
    END DO
    bet = 1./(1.5-3.5*gak(km))
    DO i=i1,i2
      qe(i, km+1) = (4.*q1(i, km)+q1(i, km-1)-3.5*qe(i, km))*bet
    END DO
    DO k=km,1,-1
      DO i=i1,i2
        qe(i, k) = qe(i, k) - gak(k)*qe(i, k+1)
      END DO
    END DO
  END SUBROUTINE EDGE_SCALAR
!  Differentiation of edge_profile in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: q1e q2e q1 q2
!   with respect to varying inputs: q1e q2e q1 q2
  SUBROUTINE EDGE_PROFILE_FWD(q1, q2, q1e, q2e, i1, i2, j1, j2, j, km, &
&   dp0, uniform_grid, limiter)
    IMPLICIT NONE
! Optimized for wind profile reconstruction:
    INTEGER, INTENT(IN) :: i1, i2, j1, j2
    INTEGER, INTENT(IN) :: j, km
    INTEGER, INTENT(IN) :: limiter
    LOGICAL, INTENT(IN) :: uniform_grid
    REAL, INTENT(IN) :: dp0(km)
    REAL, DIMENSION(i1:i2, j1:j2, km), INTENT(IN) :: q1, q2
    REAL, DIMENSION(i1:i2, j1:j2, km+1) :: q1e, q2e
!-----------------------------------------------------------------------
! edge values
    REAL, DIMENSION(i1:i2, km+1) :: qe1, qe2, gam
    REAL :: gak(km)
    REAL :: bet, r2o3, r4o3
    REAL :: g0, gk, xt1, xt2, a_bot
    INTEGER :: i, k

    qe1 = 0.0
    qe2 = 0.0
    gam = 0.0
    gak = 0.0
    bet = 0.0
    r2o3 = 0.0
    r4o3 = 0.0
    g0 = 0.0
    gk = 0.0
    xt1 = 0.0
    xt2 = 0.0
    a_bot = 0.0

    IF (uniform_grid) THEN
!------------------------------------------------
! Optimized coding for uniform grid: SJL Apr 2007
!------------------------------------------------
      r2o3 = 2./3.
      r4o3 = 4./3.
      DO i=i1,i2
        qe1(i, 1) = r4o3*q1(i, j, 1) + r2o3*q1(i, j, 2)
        qe2(i, 1) = r4o3*q2(i, j, 1) + r2o3*q2(i, j, 2)
      END DO
      gak(1) = 7./3.
      DO k=2,km
        CALL PUSHREALARRAY(gak(k))
        gak(k) = 1./(4.-gak(k-1))
        DO i=i1,i2
          qe1(i, k) = (3.*(q1(i, j, k-1)+q1(i, j, k))-qe1(i, k-1))*gak(k&
&           )
          qe2(i, k) = (3.*(q2(i, j, k-1)+q2(i, j, k))-qe2(i, k-1))*gak(k&
&           )
        END DO
      END DO
      bet = 1./(1.5-3.5*gak(km))
      DO i=i1,i2
        qe1(i, km+1) = (4.*q1(i, j, km)+q1(i, j, km-1)-3.5*qe1(i, km))*&
&         bet
        qe2(i, km+1) = (4.*q2(i, j, km)+q2(i, j, km-1)-3.5*qe2(i, km))*&
&         bet
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          qe1(i, k) = qe1(i, k) - gak(k)*qe1(i, k+1)
          qe2(i, k) = qe2(i, k) - gak(k)*qe2(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
! Assuming grid varying in vertical only
      g0 = dp0(2)/dp0(1)
      xt1 = 2.*g0*(g0+1.)
      bet = g0*(g0+0.5)
      DO i=i1,i2
        qe1(i, 1) = (xt1*q1(i, j, 1)+q1(i, j, 2))/bet
        qe2(i, 1) = (xt1*q2(i, j, 1)+q2(i, j, 2))/bet
        gam(i, 1) = (1.+g0*(g0+1.5))/bet
      END DO
      DO k=2,km
        CALL PUSHREALARRAY(gk)
        gk = dp0(k-1)/dp0(k)
        DO i=i1,i2
          CALL PUSHREALARRAY(bet)
          bet = 2. + 2.*gk - gam(i, k-1)
          qe1(i, k) = (3.*(q1(i, j, k-1)+gk*q1(i, j, k))-qe1(i, k-1))/&
&           bet
          qe2(i, k) = (3.*(q2(i, j, k-1)+gk*q2(i, j, k))-qe2(i, k-1))/&
&           bet
          gam(i, k) = gk/bet
        END DO
      END DO
      a_bot = 1. + gk*(gk+1.5)
      CALL PUSHREALARRAY(xt1)
      xt1 = 2.*gk*(gk+1.)
      DO i=i1,i2
        CALL PUSHREALARRAY(xt2)
        xt2 = gk*(gk+0.5) - a_bot*gam(i, km)
        qe1(i, km+1) = (xt1*q1(i, j, km)+q1(i, j, km-1)-a_bot*qe1(i, km)&
&         )/xt2
        qe2(i, km+1) = (xt1*q2(i, j, km)+q2(i, j, km-1)-a_bot*qe2(i, km)&
&         )/xt2
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          qe1(i, k) = qe1(i, k) - gam(i, k)*qe1(i, k+1)
          qe2(i, k) = qe2(i, k) - gam(i, k)*qe2(i, k+1)
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    END IF
!------------------
! Apply constraints
!------------------
    IF (limiter .NE. 0) THEN
! limit the top & bottom winds
      DO i=i1,i2
! Top
        IF (q1(i, j, 1)*qe1(i, 1) .LT. 0.) THEN
          qe1(i, 1) = 0.
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (q2(i, j, 1)*qe2(i, 1) .LT. 0.) THEN
          qe2(i, 1) = 0.
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! Surface:
        IF (q1(i, j, km)*qe1(i, km+1) .LT. 0.) THEN
          qe1(i, km+1) = 0.
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (q2(i, j, km)*qe2(i, km+1) .LT. 0.) THEN
          qe2(i, km+1) = 0.
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
    DO k=1,km+1
      DO i=i1,i2
        q1e(i, j, k) = qe1(i, k)
        q2e(i, j, k) = qe2(i, k)
      END DO
    END DO
    CALL PUSHREALARRAY(gam, (i2-i1+1)*(km+1))
    CALL PUSHREALARRAY(gak, km)
    CALL PUSHREALARRAY(xt2)
    CALL PUSHREALARRAY(xt1)
    CALL PUSHREALARRAY(bet)
    CALL PUSHREALARRAY(a_bot)
    CALL PUSHREALARRAY(gk)
  END SUBROUTINE EDGE_PROFILE_FWD
!  Differentiation of edge_profile in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: q1e q2e q1 q2
!   with respect to varying inputs: q1e q2e q1 q2
  SUBROUTINE EDGE_PROFILE_BWD(q1, q1_ad, q2, q2_ad, q1e, q1e_ad, q2e, &
&   q2e_ad, i1, i2, j1, j2, j, km, dp0, uniform_grid, limiter)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1, i2, j1, j2
    INTEGER, INTENT(IN) :: j, km
    INTEGER, INTENT(IN) :: limiter
    LOGICAL, INTENT(IN) :: uniform_grid
    REAL, INTENT(IN) :: dp0(km)
    REAL, DIMENSION(i1:i2, j1:j2, km), INTENT(IN) :: q1, q2
    REAL, DIMENSION(i1:i2, j1:j2, km) :: q1_ad, q2_ad
    REAL, DIMENSION(i1:i2, j1:j2, km+1) :: q1e, q2e
    REAL, DIMENSION(i1:i2, j1:j2, km+1) :: q1e_ad, q2e_ad
    REAL, DIMENSION(i1:i2, km+1) :: qe1, qe2, gam
    REAL, DIMENSION(i1:i2, km+1) :: qe1_ad, qe2_ad
    REAL :: gak(km)
    REAL :: bet, r2o3, r4o3
    REAL :: g0, gk, xt1, xt2, a_bot
    INTEGER :: i, k
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    INTEGER :: branch

    qe1 = 0.0
    qe2 = 0.0
    gam = 0.0
    gak = 0.0
    bet = 0.0
    r2o3 = 0.0
    r4o3 = 0.0
    g0 = 0.0
    gk = 0.0
    xt1 = 0.0
    xt2 = 0.0
    a_bot = 0.0
    branch = 0

    CALL POPREALARRAY(gk)
    CALL POPREALARRAY(a_bot)
    CALL POPREALARRAY(bet)
    CALL POPREALARRAY(xt1)
    CALL POPREALARRAY(xt2)
    CALL POPREALARRAY(gak, km)
    CALL POPREALARRAY(gam, (i2-i1+1)*(km+1))
    qe1_ad = 0.0
    qe2_ad = 0.0
    DO k=km+1,1,-1
      DO i=i2,i1,-1
        qe2_ad(i, k) = qe2_ad(i, k) + q2e_ad(i, j, k)
        q2e_ad(i, j, k) = 0.0
        qe1_ad(i, k) = qe1_ad(i, k) + q1e_ad(i, j, k)
        q1e_ad(i, j, k) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      DO i=i2,i1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) qe2_ad(i, km+1) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) qe1_ad(i, km+1) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) qe2_ad(i, 1) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) qe1_ad(i, 1) = 0.0
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO k=1,km,1
        DO i=i2,i1,-1
          qe2_ad(i, k+1) = qe2_ad(i, k+1) - gak(k)*qe2_ad(i, k)
          qe1_ad(i, k+1) = qe1_ad(i, k+1) - gak(k)*qe1_ad(i, k)
        END DO
      END DO
      DO i=i2,i1,-1
        temp_ad1 = bet*qe2_ad(i, km+1)
        q2_ad(i, j, km) = q2_ad(i, j, km) + 4.*temp_ad1
        q2_ad(i, j, km-1) = q2_ad(i, j, km-1) + temp_ad1
        qe2_ad(i, km) = qe2_ad(i, km) - 3.5*temp_ad1
        qe2_ad(i, km+1) = 0.0
        temp_ad2 = bet*qe1_ad(i, km+1)
        q1_ad(i, j, km) = q1_ad(i, j, km) + 4.*temp_ad2
        q1_ad(i, j, km-1) = q1_ad(i, j, km-1) + temp_ad2
        qe1_ad(i, km) = qe1_ad(i, km) - 3.5*temp_ad2
        qe1_ad(i, km+1) = 0.0
      END DO
      DO k=km,2,-1
        DO i=i2,i1,-1
          temp_ad = gak(k)*qe2_ad(i, k)
          q2_ad(i, j, k-1) = q2_ad(i, j, k-1) + 3.*temp_ad
          q2_ad(i, j, k) = q2_ad(i, j, k) + 3.*temp_ad
          qe2_ad(i, k-1) = qe2_ad(i, k-1) - temp_ad
          qe2_ad(i, k) = 0.0
          temp_ad0 = gak(k)*qe1_ad(i, k)
          q1_ad(i, j, k-1) = q1_ad(i, j, k-1) + 3.*temp_ad0
          q1_ad(i, j, k) = q1_ad(i, j, k) + 3.*temp_ad0
          qe1_ad(i, k-1) = qe1_ad(i, k-1) - temp_ad0
          qe1_ad(i, k) = 0.0
        END DO
        CALL POPREALARRAY(gak(k))
      END DO
      r2o3 = 2./3.
      r4o3 = 4./3.
      DO i=i2,i1,-1
        q2_ad(i, j, 1) = q2_ad(i, j, 1) + r4o3*qe2_ad(i, 1)
        q2_ad(i, j, 2) = q2_ad(i, j, 2) + r2o3*qe2_ad(i, 1)
        qe2_ad(i, 1) = 0.0
        q1_ad(i, j, 1) = q1_ad(i, j, 1) + r4o3*qe1_ad(i, 1)
        q1_ad(i, j, 2) = q1_ad(i, j, 2) + r2o3*qe1_ad(i, 1)
        qe1_ad(i, 1) = 0.0
      END DO
    ELSE
      DO k=1,km,1
        DO i=i2,i1,-1
          qe2_ad(i, k+1) = qe2_ad(i, k+1) - gam(i, k)*qe2_ad(i, k)
          qe1_ad(i, k+1) = qe1_ad(i, k+1) - gam(i, k)*qe1_ad(i, k)
        END DO
      END DO
      DO i=i2,i1,-1
        temp_ad5 = qe2_ad(i, km+1)/xt2
        q2_ad(i, j, km) = q2_ad(i, j, km) + xt1*temp_ad5
        q2_ad(i, j, km-1) = q2_ad(i, j, km-1) + temp_ad5
        qe2_ad(i, km) = qe2_ad(i, km) - a_bot*temp_ad5
        qe2_ad(i, km+1) = 0.0
        temp_ad6 = qe1_ad(i, km+1)/xt2
        q1_ad(i, j, km) = q1_ad(i, j, km) + xt1*temp_ad6
        q1_ad(i, j, km-1) = q1_ad(i, j, km-1) + temp_ad6
        qe1_ad(i, km) = qe1_ad(i, km) - a_bot*temp_ad6
        qe1_ad(i, km+1) = 0.0
        CALL POPREALARRAY(xt2)
      END DO
      CALL POPREALARRAY(xt1)
      DO k=km,2,-1
        DO i=i2,i1,-1
          temp_ad3 = qe2_ad(i, k)/bet
          q2_ad(i, j, k-1) = q2_ad(i, j, k-1) + 3.*temp_ad3
          q2_ad(i, j, k) = q2_ad(i, j, k) + 3.*gk*temp_ad3
          qe2_ad(i, k-1) = qe2_ad(i, k-1) - temp_ad3
          qe2_ad(i, k) = 0.0
          temp_ad4 = qe1_ad(i, k)/bet
          q1_ad(i, j, k-1) = q1_ad(i, j, k-1) + 3.*temp_ad4
          q1_ad(i, j, k) = q1_ad(i, j, k) + 3.*gk*temp_ad4
          qe1_ad(i, k-1) = qe1_ad(i, k-1) - temp_ad4
          qe1_ad(i, k) = 0.0
          CALL POPREALARRAY(bet)
        END DO
        CALL POPREALARRAY(gk)
      END DO
      DO i=i2,i1,-1
        q2_ad(i, j, 1) = q2_ad(i, j, 1) + xt1*qe2_ad(i, 1)/bet
        q2_ad(i, j, 2) = q2_ad(i, j, 2) + qe2_ad(i, 1)/bet
        qe2_ad(i, 1) = 0.0
        q1_ad(i, j, 1) = q1_ad(i, j, 1) + xt1*qe1_ad(i, 1)/bet
        q1_ad(i, j, 2) = q1_ad(i, j, 2) + qe1_ad(i, 1)/bet
        qe1_ad(i, 1) = 0.0
      END DO
    END IF
  END SUBROUTINE EDGE_PROFILE_BWD
  SUBROUTINE EDGE_PROFILE(q1, q2, q1e, q2e, i1, i2, j1, j2, j, km, dp0, &
&   uniform_grid, limiter)
    IMPLICIT NONE
! Optimized for wind profile reconstruction:
    INTEGER, INTENT(IN) :: i1, i2, j1, j2
    INTEGER, INTENT(IN) :: j, km
    INTEGER, INTENT(IN) :: limiter
    LOGICAL, INTENT(IN) :: uniform_grid
    REAL, INTENT(IN) :: dp0(km)
    REAL, DIMENSION(i1:i2, j1:j2, km), INTENT(IN) :: q1, q2
    REAL, DIMENSION(i1:i2, j1:j2, km+1), INTENT(OUT) :: q1e, q2e
!-----------------------------------------------------------------------
! edge values
    REAL, DIMENSION(i1:i2, km+1) :: qe1, qe2, gam
    REAL :: gak(km)
    REAL :: bet, r2o3, r4o3
    REAL :: g0, gk, xt1, xt2, a_bot
    INTEGER :: i, k
    IF (uniform_grid) THEN
!------------------------------------------------
! Optimized coding for uniform grid: SJL Apr 2007
!------------------------------------------------
      r2o3 = 2./3.
      r4o3 = 4./3.
      DO i=i1,i2
        qe1(i, 1) = r4o3*q1(i, j, 1) + r2o3*q1(i, j, 2)
        qe2(i, 1) = r4o3*q2(i, j, 1) + r2o3*q2(i, j, 2)
      END DO
      gak(1) = 7./3.
      DO k=2,km
        gak(k) = 1./(4.-gak(k-1))
        DO i=i1,i2
          qe1(i, k) = (3.*(q1(i, j, k-1)+q1(i, j, k))-qe1(i, k-1))*gak(k&
&           )
          qe2(i, k) = (3.*(q2(i, j, k-1)+q2(i, j, k))-qe2(i, k-1))*gak(k&
&           )
        END DO
      END DO
      bet = 1./(1.5-3.5*gak(km))
      DO i=i1,i2
        qe1(i, km+1) = (4.*q1(i, j, km)+q1(i, j, km-1)-3.5*qe1(i, km))*&
&         bet
        qe2(i, km+1) = (4.*q2(i, j, km)+q2(i, j, km-1)-3.5*qe2(i, km))*&
&         bet
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          qe1(i, k) = qe1(i, k) - gak(k)*qe1(i, k+1)
          qe2(i, k) = qe2(i, k) - gak(k)*qe2(i, k+1)
        END DO
      END DO
    ELSE
! Assuming grid varying in vertical only
      g0 = dp0(2)/dp0(1)
      xt1 = 2.*g0*(g0+1.)
      bet = g0*(g0+0.5)
      DO i=i1,i2
        qe1(i, 1) = (xt1*q1(i, j, 1)+q1(i, j, 2))/bet
        qe2(i, 1) = (xt1*q2(i, j, 1)+q2(i, j, 2))/bet
        gam(i, 1) = (1.+g0*(g0+1.5))/bet
      END DO
      DO k=2,km
        gk = dp0(k-1)/dp0(k)
        DO i=i1,i2
          bet = 2. + 2.*gk - gam(i, k-1)
          qe1(i, k) = (3.*(q1(i, j, k-1)+gk*q1(i, j, k))-qe1(i, k-1))/&
&           bet
          qe2(i, k) = (3.*(q2(i, j, k-1)+gk*q2(i, j, k))-qe2(i, k-1))/&
&           bet
          gam(i, k) = gk/bet
        END DO
      END DO
      a_bot = 1. + gk*(gk+1.5)
      xt1 = 2.*gk*(gk+1.)
      DO i=i1,i2
        xt2 = gk*(gk+0.5) - a_bot*gam(i, km)
        qe1(i, km+1) = (xt1*q1(i, j, km)+q1(i, j, km-1)-a_bot*qe1(i, km)&
&         )/xt2
        qe2(i, km+1) = (xt1*q2(i, j, km)+q2(i, j, km-1)-a_bot*qe2(i, km)&
&         )/xt2
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          qe1(i, k) = qe1(i, k) - gam(i, k)*qe1(i, k+1)
          qe2(i, k) = qe2(i, k) - gam(i, k)*qe2(i, k+1)
        END DO
      END DO
    END IF
!------------------
! Apply constraints
!------------------
    IF (limiter .NE. 0) THEN
! limit the top & bottom winds
      DO i=i1,i2
! Top
        IF (q1(i, j, 1)*qe1(i, 1) .LT. 0.) qe1(i, 1) = 0.
        IF (q2(i, j, 1)*qe2(i, 1) .LT. 0.) qe2(i, 1) = 0.
! Surface:
        IF (q1(i, j, km)*qe1(i, km+1) .LT. 0.) qe1(i, km+1) = 0.
        IF (q2(i, j, km)*qe2(i, km+1) .LT. 0.) qe2(i, km+1) = 0.
      END DO
    END IF
    DO k=1,km+1
      DO i=i1,i2
        q1e(i, j, k) = qe1(i, k)
        q2e(i, j, k) = qe2(i, k)
      END DO
    END DO
  END SUBROUTINE EDGE_PROFILE
!  Differentiation of nest_halo_nh in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: pk3 gz delp delz pkc pt
!   with respect to varying inputs: pk3 gz delp delz pkc pt
  SUBROUTINE NEST_HALO_NH_FWD(ptop, grav, kappa, cp, delp, delz, pt, &
&   phis, pkc, gz, pk3, npx, npy, npz, nested, pkc_pertn, computepk3, &
&   fullhalo, bd)
    IMPLICIT NONE
!INPUT: delp, delz, pt
!OUTPUT: gz, pkc, pk3 (optional)
    INTEGER, INTENT(IN) :: npx, npy, npz
    LOGICAL, INTENT(IN) :: pkc_pertn, computepk3, fullhalo, nested
    REAL, INTENT(IN) :: ptop, kappa, cp, grav
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(IN) :: pt&
&   , delp, delz
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(INOUT) &
&   :: gz, pkc, pk3
    INTEGER :: i, j, k
!'gamma'
    REAL :: gama
    REAL :: ptk, rgrav, rkap, peln1, rdg
    REAL, DIMENSION(bd%isd:bd%ied, npz+1, bd%jsd:bd%jed) :: pe, peln
    REAL, DIMENSION(bd%isd:bd%ied, npz) :: gam, bb, dd, pkz
    REAL, DIMENSION(bd%isd:bd%ied, npz-1) :: g_rat
    REAL, DIMENSION(bd%isd:bd%ied) :: bet
    REAL :: pm
    INTEGER :: ifirst, ilast, jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC LOG
    INTRINSIC EXP

    gama = 0.0
    ptk = 0.0
    rgrav = 0.0
    rkap = 0.0
    peln1 = 0.0
    rdg = 0.0
    pe = 0.0
    peln = 0.0
    gam = 0.0
    bb = 0.0
    dd = 0.0
    pkz = 0.0
    g_rat = 0.0
    bet = 0.0
    pm = 0.0
    ifirst = 0
    ilast = 0
    jfirst = 0
    jlast = 0
    is = 0
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
    IF (.NOT.nested) THEN
      CALL PUSHCONTROL(2,0)
    ELSE
      ifirst = isd
      jfirst = jsd
      ilast = ied
      jlast = jed
!Remember we want to compute these in the HALO. Note also this routine
!requires an appropriate
      rgrav = 1./grav
      gama = 1./(1.-kappa)
      ptk = ptop**kappa
      peln1 = LOG(ptop)
!NOTE: Compiler does NOT like this sort of nested-grid BC code. Is it trying to do some ugly optimization?
      IF (is .EQ. 1) THEN
        DO j=jfirst,jlast
!GZ
          DO i=ifirst,0
            CALL PUSHREALARRAY(gz(i, j, npz+1))
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=ifirst,0
              CALL PUSHREALARRAY(gz(i, j, k))
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=ifirst,0
            pe(i, 1, j) = ptop
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=ifirst,0
              CALL PUSHREALARRAY(pe(i, k, j))
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=ifirst,0
!Full p
              CALL PUSHREALARRAY(pkz(i, k))
              pkz(i, k) = EXP(gama*LOG(-(delp(i, j, k)*rgrav/delz(i, j, &
&               k)*rdgas*pt(i, j, k))))
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!pressure solver
          DO k=1,npz-1
            DO i=ifirst,0
              CALL PUSHREALARRAY(g_rat(i, k))
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              CALL PUSHREALARRAY(dd(i, k))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=ifirst,0
            CALL PUSHREALARRAY(bet(i))
            bet(i) = bb(i, 1)
            CALL PUSHREALARRAY(pkc(i, j, 1))
            pkc(i, j, 1) = 0.
            CALL PUSHREALARRAY(pkc(i, j, 2))
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb(i, npz) = 2.
            CALL PUSHREALARRAY(dd(i, npz))
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=ifirst,0
              CALL PUSHREALARRAY(gam(i, k))
              gam(i, k) = g_rat(i, k-1)/bet(i)
              CALL PUSHREALARRAY(bet(i))
              bet(i) = bb(i, k) - gam(i, k)
              CALL PUSHREALARRAY(pkc(i, j, k+1))
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ifirst,0
              CALL PUSHREALARRAY(pkc(i, j, k))
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=jfirst,jlast
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=ifirst,0
                CALL PUSHREALARRAY(pkc(i, j, k))
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
!pk3 if necessary; doesn't require condenstate loading calculation
          IF (computepk3) THEN
            DO i=ifirst,0
              CALL PUSHREALARRAY(pk3(i, j, 1))
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=ifirst,0
                CALL PUSHREALARRAY(pk3(i, j, k))
                pk3(i, j, k) = EXP(kappa*LOG(pe(i, k, j)))
              END DO
            END DO
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ie .EQ. npx - 1) THEN
        DO j=jfirst,jlast
!GZ
          DO i=npx,ilast
            CALL PUSHREALARRAY(gz(i, j, npz+1))
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=npx,ilast
              CALL PUSHREALARRAY(gz(i, j, k))
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=npx,ilast
            CALL PUSHREALARRAY(pe(i, 1, j))
            pe(i, 1, j) = ptop
            CALL PUSHREALARRAY(peln(i, 1, j))
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=npx,ilast
              CALL PUSHREALARRAY(pe(i, k, j))
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              CALL PUSHREALARRAY(peln(i, k, j))
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=npx,ilast
!Full p
              CALL PUSHREALARRAY(pkz(i, k))
              pkz(i, k) = EXP(gama*LOG(-(delp(i, j, k)*rgrav/delz(i, j, &
&               k)*rdgas*pt(i, j, k))))
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!pressure solver
          DO k=1,npz-1
            DO i=npx,ilast
              CALL PUSHREALARRAY(g_rat(i, k))
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              CALL PUSHREALARRAY(dd(i, k))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=npx,ilast
            CALL PUSHREALARRAY(bet(i))
            bet(i) = bb(i, 1)
            CALL PUSHREALARRAY(pkc(i, j, 1))
            pkc(i, j, 1) = 0.
            CALL PUSHREALARRAY(pkc(i, j, 2))
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb(i, npz) = 2.
            CALL PUSHREALARRAY(dd(i, npz))
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=npx,ilast
              CALL PUSHREALARRAY(gam(i, k))
              gam(i, k) = g_rat(i, k-1)/bet(i)
              CALL PUSHREALARRAY(bet(i))
              bet(i) = bb(i, k) - gam(i, k)
              CALL PUSHREALARRAY(pkc(i, j, k+1))
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=npx,ilast
              CALL PUSHREALARRAY(pkc(i, j, k))
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=jfirst,jlast
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=npx,ilast
                CALL PUSHREALARRAY(pkc(i, j, k))
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
!pk3 if necessary
          IF (computepk3) THEN
            DO i=npx,ilast
              CALL PUSHREALARRAY(pk3(i, j, 1))
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=npx,ilast
                CALL PUSHREALARRAY(pk3(i, j, k))
                pk3(i, j, k) = EXP(kappa*LOG(pe(i, k, j)))
              END DO
            END DO
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (js .EQ. 1) THEN
        DO j=jfirst,0
!GZ
          DO i=ifirst,ilast
            CALL PUSHREALARRAY(gz(i, j, npz+1))
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(gz(i, j, k))
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=ifirst,ilast
            CALL PUSHREALARRAY(pe(i, 1, j))
            pe(i, 1, j) = ptop
            CALL PUSHREALARRAY(peln(i, 1, j))
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(pe(i, k, j))
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              CALL PUSHREALARRAY(peln(i, k, j))
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=ifirst,ilast
!Full p
              CALL PUSHREALARRAY(pkz(i, k))
              pkz(i, k) = EXP(gama*LOG(-(delp(i, j, k)*rgrav/delz(i, j, &
&               k)*rdgas*pt(i, j, k))))
!hydro
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!pressure solver
          DO k=1,npz-1
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(g_rat(i, k))
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              CALL PUSHREALARRAY(dd(i, k))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=ifirst,ilast
            CALL PUSHREALARRAY(bet(i))
            bet(i) = bb(i, 1)
            CALL PUSHREALARRAY(pkc(i, j, 1))
            pkc(i, j, 1) = 0.
            CALL PUSHREALARRAY(pkc(i, j, 2))
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb(i, npz) = 2.
            CALL PUSHREALARRAY(dd(i, npz))
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(gam(i, k))
              gam(i, k) = g_rat(i, k-1)/bet(i)
              CALL PUSHREALARRAY(bet(i))
              bet(i) = bb(i, k) - gam(i, k)
              CALL PUSHREALARRAY(pkc(i, j, k+1))
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(pkc(i, j, k))
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=jfirst,0
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=ifirst,ilast
                CALL PUSHREALARRAY(pkc(i, j, k))
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
!pk3 if necessary
          IF (computepk3) THEN
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(pk3(i, j, 1))
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=ifirst,ilast
                CALL PUSHREALARRAY(pk3(i, j, k))
                pk3(i, j, k) = EXP(kappa*LOG(pe(i, k, j)))
              END DO
            END DO
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (je .EQ. npy - 1) THEN
        DO j=npy,jlast
!GZ
          DO i=ifirst,ilast
            CALL PUSHREALARRAY(gz(i, j, npz+1))
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(gz(i, j, k))
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=ifirst,ilast
            CALL PUSHREALARRAY(pe(i, 1, j))
            pe(i, 1, j) = ptop
            CALL PUSHREALARRAY(peln(i, 1, j))
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(pe(i, k, j))
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              CALL PUSHREALARRAY(peln(i, k, j))
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=ifirst,ilast
!Full p
              CALL PUSHREALARRAY(pkz(i, k))
              pkz(i, k) = EXP(gama*LOG(-(delp(i, j, k)*rgrav/delz(i, j, &
&               k)*rdgas*pt(i, j, k))))
!hydro
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!Reversible interpolation on layer NH pressure perturbation
!                 to recover  lastge NH pressure perturbation
          DO k=1,npz-1
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(g_rat(i, k))
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              CALL PUSHREALARRAY(dd(i, k))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=ifirst,ilast
            CALL PUSHREALARRAY(bet(i))
            bet(i) = bb(i, 1)
            CALL PUSHREALARRAY(pkc(i, j, 1))
            pkc(i, j, 1) = 0.
            CALL PUSHREALARRAY(pkc(i, j, 2))
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb(i, npz) = 2.
            CALL PUSHREALARRAY(dd(i, npz))
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(gam(i, k))
              gam(i, k) = g_rat(i, k-1)/bet(i)
              CALL PUSHREALARRAY(bet(i))
              bet(i) = bb(i, k) - gam(i, k)
              CALL PUSHREALARRAY(pkc(i, j, k+1))
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(pkc(i, j, k))
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=npy,jlast
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=ifirst,ilast
                CALL PUSHREALARRAY(pkc(i, j, k))
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
!pk3 if necessary
          IF (computepk3) THEN
            DO i=ifirst,ilast
              CALL PUSHREALARRAY(pk3(i, j, 1))
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=ifirst,ilast
                CALL PUSHREALARRAY(pk3(i, j, k))
                pk3(i, j, k) = EXP(kappa*LOG(pe(i, k, j)))
              END DO
            END DO
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        CALL PUSHREALARRAY(gam, (bd%ied-bd%isd+1)*npz)
        CALL PUSHINTEGER(ifirst)
        CALL PUSHREALARRAY(g_rat, (bd%ied-bd%isd+1)*(npz-1))
        CALL PUSHREALARRAY(pe, (bd%ied-bd%isd+1)*(npz+1)*(bd%jed-bd%jsd&
&                     +1))
        CALL PUSHREALARRAY(pkz, (bd%ied-bd%isd+1)*npz)
        CALL PUSHREALARRAY(gama)
        CALL PUSHINTEGER(jlast)
        CALL PUSHREALARRAY(rgrav)
        CALL PUSHINTEGER(ilast)
        CALL PUSHREALARRAY(bet, bd%ied - bd%isd + 1)
        CALL PUSHREALARRAY(peln, (bd%ied-bd%isd+1)*(npz+1)*(bd%jed-bd%&
&                     jsd+1))
        CALL PUSHREALARRAY(dd, (bd%ied-bd%isd+1)*npz)
        CALL PUSHCONTROL(2,2)
      ELSE
        CALL PUSHREALARRAY(gam, (bd%ied-bd%isd+1)*npz)
        CALL PUSHINTEGER(ifirst)
        CALL PUSHREALARRAY(g_rat, (bd%ied-bd%isd+1)*(npz-1))
        CALL PUSHREALARRAY(pe, (bd%ied-bd%isd+1)*(npz+1)*(bd%jed-bd%jsd&
&                     +1))
        CALL PUSHREALARRAY(pkz, (bd%ied-bd%isd+1)*npz)
        CALL PUSHREALARRAY(gama)
        CALL PUSHINTEGER(jlast)
        CALL PUSHREALARRAY(rgrav)
        CALL PUSHINTEGER(ilast)
        CALL PUSHREALARRAY(bet, bd%ied - bd%isd + 1)
        CALL PUSHREALARRAY(peln, (bd%ied-bd%isd+1)*(npz+1)*(bd%jed-bd%&
&                     jsd+1))
        CALL PUSHREALARRAY(dd, (bd%ied-bd%isd+1)*npz)
        CALL PUSHCONTROL(2,1)
      END IF
    END IF
  END SUBROUTINE NEST_HALO_NH_FWD
!  Differentiation of nest_halo_nh in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: pk3 gz delp delz pkc pt
!   with respect to varying inputs: pk3 gz delp delz pkc pt
  SUBROUTINE NEST_HALO_NH_BWD(ptop, grav, kappa, cp, delp, delp_ad, delz&
&   , delz_ad, pt, pt_ad, phis, pkc, pkc_ad, gz, gz_ad, pk3, pk3_ad, npx&
&   , npy, npz, nested, pkc_pertn, computepk3, fullhalo, bd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, npz
    LOGICAL, INTENT(IN) :: pkc_pertn, computepk3, fullhalo, nested
    REAL, INTENT(IN) :: ptop, kappa, cp, grav
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(IN) :: pt&
&   , delp, delz
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz) :: pt_ad, delp_ad&
&   , delz_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(INOUT) &
&   :: gz, pkc, pk3
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(INOUT) &
&   :: gz_ad, pkc_ad, pk3_ad
    INTEGER :: i, j, k
    REAL :: gama
    REAL :: ptk, rgrav, rkap, peln1, rdg
    REAL, DIMENSION(bd%isd:bd%ied, npz+1, bd%jsd:bd%jed) :: pe, peln
    REAL, DIMENSION(bd%isd:bd%ied, npz+1, bd%jsd:bd%jed) :: pe_ad, &
&   peln_ad
    REAL, DIMENSION(bd%isd:bd%ied, npz) :: gam, bb, dd, pkz
    REAL, DIMENSION(bd%isd:bd%ied, npz) :: gam_ad, bb_ad, dd_ad, pkz_ad
    REAL, DIMENSION(bd%isd:bd%ied, npz-1) :: g_rat
    REAL, DIMENSION(bd%isd:bd%ied, npz-1) :: g_rat_ad
    REAL, DIMENSION(bd%isd:bd%ied) :: bet
    REAL, DIMENSION(bd%isd:bd%ied) :: bet_ad
    REAL :: pm
    REAL :: pm_ad
    INTEGER :: ifirst, ilast, jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC LOG
    INTRINSIC EXP
    REAL :: temp
    REAL :: temp0
    REAL :: temp1
    REAL :: temp2
    REAL :: temp3
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp4
    REAL :: temp5
    REAL :: temp6
    REAL :: temp7
    REAL :: temp8
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp9
    REAL :: temp10
    REAL :: temp11
    REAL :: temp12
    REAL :: temp13
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: temp_ad16
    REAL :: temp_ad17
    REAL :: temp_ad18
    REAL :: temp_ad19
    REAL :: temp14
    REAL :: temp15
    REAL :: temp16
    REAL :: temp17
    REAL :: temp18
    REAL :: temp_ad20
    REAL :: temp_ad21
    REAL :: temp_ad22
    REAL :: temp_ad23
    REAL :: temp_ad24
    REAL :: temp_ad25
    REAL :: temp_ad26
    INTEGER :: branch

    gama = 0.0
    ptk = 0.0
    rgrav = 0.0
    rkap = 0.0
    peln1 = 0.0
    rdg = 0.0
    pe = 0.0
    peln = 0.0
    gam = 0.0
    bb = 0.0
    dd = 0.0
    pkz = 0.0
    g_rat = 0.0
    bet = 0.0
    pm = 0.0
    ifirst = 0
    ilast = 0
    jfirst = 0
    jlast = 0
    is = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    branch = 0

    CALL POPCONTROL(2,branch)
    IF (branch .NE. 0) THEN
      IF (branch .EQ. 1) THEN
        CALL POPREALARRAY(dd, (bd%ied-bd%isd+1)*npz)
        CALL POPREALARRAY(peln, (bd%ied-bd%isd+1)*(npz+1)*(bd%jed-bd%&
&                    jsd+1))
        CALL POPREALARRAY(bet, bd%ied - bd%isd + 1)
        CALL POPINTEGER(ilast)
        CALL POPREALARRAY(rgrav)
        CALL POPINTEGER(jlast)
        CALL POPREALARRAY(gama)
        CALL POPREALARRAY(pkz, (bd%ied-bd%isd+1)*npz)
        CALL POPREALARRAY(pe, (bd%ied-bd%isd+1)*(npz+1)*(bd%jed-bd%jsd+&
&                    1))
        CALL POPREALARRAY(g_rat, (bd%ied-bd%isd+1)*(npz-1))
        CALL POPINTEGER(ifirst)
        CALL POPREALARRAY(gam, (bd%ied-bd%isd+1)*npz)
        dd_ad = 0.0
        peln_ad = 0.0
        bet_ad = 0.0
        bb_ad = 0.0
        pkz_ad = 0.0
        pe_ad = 0.0
        g_rat_ad = 0.0
        gam_ad = 0.0
      ELSE
        CALL POPREALARRAY(dd, (bd%ied-bd%isd+1)*npz)
        CALL POPREALARRAY(peln, (bd%ied-bd%isd+1)*(npz+1)*(bd%jed-bd%&
&                    jsd+1))
        CALL POPREALARRAY(bet, bd%ied - bd%isd + 1)
        CALL POPINTEGER(ilast)
        CALL POPREALARRAY(rgrav)
        CALL POPINTEGER(jlast)
        CALL POPREALARRAY(gama)
        CALL POPREALARRAY(pkz, (bd%ied-bd%isd+1)*npz)
        CALL POPREALARRAY(pe, (bd%ied-bd%isd+1)*(npz+1)*(bd%jed-bd%jsd+&
&                    1))
        CALL POPREALARRAY(g_rat, (bd%ied-bd%isd+1)*(npz-1))
        CALL POPINTEGER(ifirst)
        CALL POPREALARRAY(gam, (bd%ied-bd%isd+1)*npz)
        pe_ad = 0.0
        DO j=jlast,npy,-1
          CALL POPCONTROL(1,branch)
          IF (branch .NE. 0) THEN
            DO k=npz+1,2,-1
              DO i=ilast,ifirst,-1
                CALL POPREALARRAY(pk3(i, j, k))
                pe_ad(i, k, j) = pe_ad(i, k, j) + kappa*EXP(kappa*LOG(pe&
&                 (i, k, j)))*pk3_ad(i, j, k)/pe(i, k, j)
                pk3_ad(i, j, k) = 0.0
              END DO
            END DO
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(pk3(i, j, 1))
              pk3_ad(i, j, 1) = 0.0
            END DO
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            DO k=1,npz+1,1
              DO i=ilast,ifirst,-1
                CALL POPREALARRAY(pkc(i, j, k))
                pe_ad(i, k, j) = pe_ad(i, k, j) + pkc_ad(i, j, k)
              END DO
            END DO
          END IF
        END DO
        dd_ad = 0.0
        peln_ad = 0.0
        bet_ad = 0.0
        bb_ad = 0.0
        pkz_ad = 0.0
        g_rat_ad = 0.0
        gam_ad = 0.0
        DO j=jlast,npy,-1
          DO k=2,npz,1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(pkc(i, j, k))
              gam_ad(i, k) = gam_ad(i, k) - pkc(i, j, k+1)*pkc_ad(i, j, &
&               k)
              pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) - gam(i, k)*pkc_ad(i&
&               , j, k)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(pkc(i, j, k+1))
              temp_ad25 = pkc_ad(i, j, k+1)/bet(i)
              dd_ad(i, k) = dd_ad(i, k) + temp_ad25
              pkc_ad(i, j, k) = pkc_ad(i, j, k) - temp_ad25
              bet_ad(i) = bet_ad(i) - (dd(i, k)-pkc(i, j, k))*temp_ad25/&
&               bet(i)
              pkc_ad(i, j, k+1) = 0.0
              CALL POPREALARRAY(bet(i))
              bb_ad(i, k) = bb_ad(i, k) + bet_ad(i)
              gam_ad(i, k) = gam_ad(i, k) - bet_ad(i)
              CALL POPREALARRAY(gam(i, k))
              temp_ad26 = gam_ad(i, k)/bet(i)
              bet_ad(i) = -(g_rat(i, k-1)*temp_ad26/bet(i))
              g_rat_ad(i, k-1) = g_rat_ad(i, k-1) + temp_ad26
              gam_ad(i, k) = 0.0
            END DO
          END DO
          DO i=ilast,ifirst,-1
            CALL POPREALARRAY(dd(i, npz))
            pkz_ad(i, npz) = pkz_ad(i, npz) + 3.*dd_ad(i, npz)
            dd_ad(i, npz) = 0.0
            bb_ad(i, npz) = 0.0
            CALL POPREALARRAY(pkc(i, j, 2))
            temp_ad24 = pkc_ad(i, j, 2)/bet(i)
            dd_ad(i, 1) = dd_ad(i, 1) + temp_ad24
            bet_ad(i) = bet_ad(i) - dd(i, 1)*temp_ad24/bet(i)
            pkc_ad(i, j, 2) = 0.0
            CALL POPREALARRAY(pkc(i, j, 1))
            pkc_ad(i, j, 1) = 0.0
            CALL POPREALARRAY(bet(i))
            bb_ad(i, 1) = bb_ad(i, 1) + bet_ad(i)
            bet_ad(i) = 0.0
          END DO
          DO k=npz-1,1,-1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(dd(i, k))
              temp_ad22 = 3.*dd_ad(i, k)
              pkz_ad(i, k) = pkz_ad(i, k) + temp_ad22
              g_rat_ad(i, k) = g_rat_ad(i, k) + 2.*bb_ad(i, k) + pkz(i, &
&               k+1)*temp_ad22
              pkz_ad(i, k+1) = pkz_ad(i, k+1) + g_rat(i, k)*temp_ad22
              dd_ad(i, k) = 0.0
              bb_ad(i, k) = 0.0
              CALL POPREALARRAY(g_rat(i, k))
              temp_ad23 = g_rat_ad(i, k)/delp(i, j, k+1)
              delp_ad(i, j, k) = delp_ad(i, j, k) + temp_ad23
              delp_ad(i, j, k+1) = delp_ad(i, j, k+1) - delp(i, j, k)*&
&               temp_ad23/delp(i, j, k+1)
              g_rat_ad(i, k) = 0.0
            END DO
          END DO
          DO k=npz,1,-1
            DO i=ilast,ifirst,-1
              temp17 = delz(i, j, k)
              temp16 = delp(i, j, k)*pt(i, j, k)
              temp14 = temp16/temp17
              temp15 = -(rgrav*rdgas*temp14)
              temp_ad21 = -(rgrav*rdgas*gama*EXP(gama*LOG(temp15))*&
&               pkz_ad(i, k)/(temp15*temp17))
              pm_ad = -pkz_ad(i, k)
              temp18 = peln(i, k+1, j) - peln(i, k, j)
              temp_ad20 = -(delp(i, j, k)*pm_ad/temp18**2)
              delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*&
&               temp_ad21 + pm_ad/temp18
              peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad20
              peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad20
              CALL POPREALARRAY(pkz(i, k))
              pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*temp_ad21
              delz_ad(i, j, k) = delz_ad(i, j, k) - temp14*temp_ad21
              pkz_ad(i, k) = 0.0
            END DO
          END DO
          DO k=npz+1,2,-1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(peln(i, k, j))
              pe_ad(i, k, j) = pe_ad(i, k, j) + peln_ad(i, k, j)/pe(i, k&
&               , j)
              peln_ad(i, k, j) = 0.0
              CALL POPREALARRAY(pe(i, k, j))
              pe_ad(i, k-1, j) = pe_ad(i, k-1, j) + pe_ad(i, k, j)
              delp_ad(i, j, k-1) = delp_ad(i, j, k-1) + pe_ad(i, k, j)
              pe_ad(i, k, j) = 0.0
            END DO
          END DO
          DO i=ilast,ifirst,-1
            CALL POPREALARRAY(peln(i, 1, j))
            peln_ad(i, 1, j) = 0.0
            CALL POPREALARRAY(pe(i, 1, j))
            pe_ad(i, 1, j) = 0.0
          END DO
          DO k=1,npz,1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(gz(i, j, k))
              gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + gz_ad(i, j, k)
              delz_ad(i, j, k) = delz_ad(i, j, k) - grav*gz_ad(i, j, k)
              gz_ad(i, j, k) = 0.0
            END DO
          END DO
          DO i=ilast,ifirst,-1
            CALL POPREALARRAY(gz(i, j, npz+1))
            gz_ad(i, j, npz+1) = 0.0
          END DO
        END DO
      END IF
      jsd = bd%jsd
      jfirst = jsd
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=0,jfirst,-1
          CALL POPCONTROL(1,branch)
          IF (branch .NE. 0) THEN
            DO k=npz+1,2,-1
              DO i=ilast,ifirst,-1
                CALL POPREALARRAY(pk3(i, j, k))
                pe_ad(i, k, j) = pe_ad(i, k, j) + kappa*EXP(kappa*LOG(pe&
&                 (i, k, j)))*pk3_ad(i, j, k)/pe(i, k, j)
                pk3_ad(i, j, k) = 0.0
              END DO
            END DO
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(pk3(i, j, 1))
              pk3_ad(i, j, 1) = 0.0
            END DO
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            DO k=1,npz+1,1
              DO i=ilast,ifirst,-1
                CALL POPREALARRAY(pkc(i, j, k))
                pe_ad(i, k, j) = pe_ad(i, k, j) + pkc_ad(i, j, k)
              END DO
            END DO
          END IF
        END DO
        DO j=0,jfirst,-1
          DO k=2,npz,1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(pkc(i, j, k))
              gam_ad(i, k) = gam_ad(i, k) - pkc(i, j, k+1)*pkc_ad(i, j, &
&               k)
              pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) - gam(i, k)*pkc_ad(i&
&               , j, k)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(pkc(i, j, k+1))
              temp_ad18 = pkc_ad(i, j, k+1)/bet(i)
              dd_ad(i, k) = dd_ad(i, k) + temp_ad18
              pkc_ad(i, j, k) = pkc_ad(i, j, k) - temp_ad18
              bet_ad(i) = bet_ad(i) - (dd(i, k)-pkc(i, j, k))*temp_ad18/&
&               bet(i)
              pkc_ad(i, j, k+1) = 0.0
              CALL POPREALARRAY(bet(i))
              bb_ad(i, k) = bb_ad(i, k) + bet_ad(i)
              gam_ad(i, k) = gam_ad(i, k) - bet_ad(i)
              CALL POPREALARRAY(gam(i, k))
              temp_ad19 = gam_ad(i, k)/bet(i)
              bet_ad(i) = -(g_rat(i, k-1)*temp_ad19/bet(i))
              g_rat_ad(i, k-1) = g_rat_ad(i, k-1) + temp_ad19
              gam_ad(i, k) = 0.0
            END DO
          END DO
          DO i=ilast,ifirst,-1
            CALL POPREALARRAY(dd(i, npz))
            pkz_ad(i, npz) = pkz_ad(i, npz) + 3.*dd_ad(i, npz)
            dd_ad(i, npz) = 0.0
            bb_ad(i, npz) = 0.0
            CALL POPREALARRAY(pkc(i, j, 2))
            temp_ad17 = pkc_ad(i, j, 2)/bet(i)
            dd_ad(i, 1) = dd_ad(i, 1) + temp_ad17
            bet_ad(i) = bet_ad(i) - dd(i, 1)*temp_ad17/bet(i)
            pkc_ad(i, j, 2) = 0.0
            CALL POPREALARRAY(pkc(i, j, 1))
            pkc_ad(i, j, 1) = 0.0
            CALL POPREALARRAY(bet(i))
            bb_ad(i, 1) = bb_ad(i, 1) + bet_ad(i)
            bet_ad(i) = 0.0
          END DO
          DO k=npz-1,1,-1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(dd(i, k))
              temp_ad15 = 3.*dd_ad(i, k)
              pkz_ad(i, k) = pkz_ad(i, k) + temp_ad15
              g_rat_ad(i, k) = g_rat_ad(i, k) + 2.*bb_ad(i, k) + pkz(i, &
&               k+1)*temp_ad15
              pkz_ad(i, k+1) = pkz_ad(i, k+1) + g_rat(i, k)*temp_ad15
              dd_ad(i, k) = 0.0
              bb_ad(i, k) = 0.0
              CALL POPREALARRAY(g_rat(i, k))
              temp_ad16 = g_rat_ad(i, k)/delp(i, j, k+1)
              delp_ad(i, j, k) = delp_ad(i, j, k) + temp_ad16
              delp_ad(i, j, k+1) = delp_ad(i, j, k+1) - delp(i, j, k)*&
&               temp_ad16/delp(i, j, k+1)
              g_rat_ad(i, k) = 0.0
            END DO
          END DO
          DO k=npz,1,-1
            DO i=ilast,ifirst,-1
              temp12 = delz(i, j, k)
              temp11 = delp(i, j, k)*pt(i, j, k)
              temp9 = temp11/temp12
              temp10 = -(rgrav*rdgas*temp9)
              temp_ad14 = -(rgrav*rdgas*gama*EXP(gama*LOG(temp10))*&
&               pkz_ad(i, k)/(temp10*temp12))
              pm_ad = -pkz_ad(i, k)
              temp13 = peln(i, k+1, j) - peln(i, k, j)
              temp_ad13 = -(delp(i, j, k)*pm_ad/temp13**2)
              delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*&
&               temp_ad14 + pm_ad/temp13
              peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad13
              peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad13
              CALL POPREALARRAY(pkz(i, k))
              pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*temp_ad14
              delz_ad(i, j, k) = delz_ad(i, j, k) - temp9*temp_ad14
              pkz_ad(i, k) = 0.0
            END DO
          END DO
          DO k=npz+1,2,-1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(peln(i, k, j))
              pe_ad(i, k, j) = pe_ad(i, k, j) + peln_ad(i, k, j)/pe(i, k&
&               , j)
              peln_ad(i, k, j) = 0.0
              CALL POPREALARRAY(pe(i, k, j))
              pe_ad(i, k-1, j) = pe_ad(i, k-1, j) + pe_ad(i, k, j)
              delp_ad(i, j, k-1) = delp_ad(i, j, k-1) + pe_ad(i, k, j)
              pe_ad(i, k, j) = 0.0
            END DO
          END DO
          DO i=ilast,ifirst,-1
            CALL POPREALARRAY(peln(i, 1, j))
            peln_ad(i, 1, j) = 0.0
            CALL POPREALARRAY(pe(i, 1, j))
            pe_ad(i, 1, j) = 0.0
          END DO
          DO k=1,npz,1
            DO i=ilast,ifirst,-1
              CALL POPREALARRAY(gz(i, j, k))
              gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + gz_ad(i, j, k)
              delz_ad(i, j, k) = delz_ad(i, j, k) - grav*gz_ad(i, j, k)
              gz_ad(i, j, k) = 0.0
            END DO
          END DO
          DO i=ilast,ifirst,-1
            CALL POPREALARRAY(gz(i, j, npz+1))
            gz_ad(i, j, npz+1) = 0.0
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=jlast,jfirst,-1
          CALL POPCONTROL(1,branch)
          IF (branch .NE. 0) THEN
            DO k=npz+1,2,-1
              DO i=ilast,npx,-1
                CALL POPREALARRAY(pk3(i, j, k))
                pe_ad(i, k, j) = pe_ad(i, k, j) + kappa*EXP(kappa*LOG(pe&
&                 (i, k, j)))*pk3_ad(i, j, k)/pe(i, k, j)
                pk3_ad(i, j, k) = 0.0
              END DO
            END DO
            DO i=ilast,npx,-1
              CALL POPREALARRAY(pk3(i, j, 1))
              pk3_ad(i, j, 1) = 0.0
            END DO
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            DO k=1,npz+1,1
              DO i=ilast,npx,-1
                CALL POPREALARRAY(pkc(i, j, k))
                pe_ad(i, k, j) = pe_ad(i, k, j) + pkc_ad(i, j, k)
              END DO
            END DO
          END IF
        END DO
        DO j=jlast,jfirst,-1
          DO k=2,npz,1
            DO i=ilast,npx,-1
              CALL POPREALARRAY(pkc(i, j, k))
              gam_ad(i, k) = gam_ad(i, k) - pkc(i, j, k+1)*pkc_ad(i, j, &
&               k)
              pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) - gam(i, k)*pkc_ad(i&
&               , j, k)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ilast,npx,-1
              CALL POPREALARRAY(pkc(i, j, k+1))
              temp_ad11 = pkc_ad(i, j, k+1)/bet(i)
              dd_ad(i, k) = dd_ad(i, k) + temp_ad11
              pkc_ad(i, j, k) = pkc_ad(i, j, k) - temp_ad11
              bet_ad(i) = bet_ad(i) - (dd(i, k)-pkc(i, j, k))*temp_ad11/&
&               bet(i)
              pkc_ad(i, j, k+1) = 0.0
              CALL POPREALARRAY(bet(i))
              bb_ad(i, k) = bb_ad(i, k) + bet_ad(i)
              gam_ad(i, k) = gam_ad(i, k) - bet_ad(i)
              CALL POPREALARRAY(gam(i, k))
              temp_ad12 = gam_ad(i, k)/bet(i)
              bet_ad(i) = -(g_rat(i, k-1)*temp_ad12/bet(i))
              g_rat_ad(i, k-1) = g_rat_ad(i, k-1) + temp_ad12
              gam_ad(i, k) = 0.0
            END DO
          END DO
          DO i=ilast,npx,-1
            CALL POPREALARRAY(dd(i, npz))
            pkz_ad(i, npz) = pkz_ad(i, npz) + 3.*dd_ad(i, npz)
            dd_ad(i, npz) = 0.0
            bb_ad(i, npz) = 0.0
            CALL POPREALARRAY(pkc(i, j, 2))
            temp_ad10 = pkc_ad(i, j, 2)/bet(i)
            dd_ad(i, 1) = dd_ad(i, 1) + temp_ad10
            bet_ad(i) = bet_ad(i) - dd(i, 1)*temp_ad10/bet(i)
            pkc_ad(i, j, 2) = 0.0
            CALL POPREALARRAY(pkc(i, j, 1))
            pkc_ad(i, j, 1) = 0.0
            CALL POPREALARRAY(bet(i))
            bb_ad(i, 1) = bb_ad(i, 1) + bet_ad(i)
            bet_ad(i) = 0.0
          END DO
          DO k=npz-1,1,-1
            DO i=ilast,npx,-1
              CALL POPREALARRAY(dd(i, k))
              temp_ad8 = 3.*dd_ad(i, k)
              pkz_ad(i, k) = pkz_ad(i, k) + temp_ad8
              g_rat_ad(i, k) = g_rat_ad(i, k) + 2.*bb_ad(i, k) + pkz(i, &
&               k+1)*temp_ad8
              pkz_ad(i, k+1) = pkz_ad(i, k+1) + g_rat(i, k)*temp_ad8
              dd_ad(i, k) = 0.0
              bb_ad(i, k) = 0.0
              CALL POPREALARRAY(g_rat(i, k))
              temp_ad9 = g_rat_ad(i, k)/delp(i, j, k+1)
              delp_ad(i, j, k) = delp_ad(i, j, k) + temp_ad9
              delp_ad(i, j, k+1) = delp_ad(i, j, k+1) - delp(i, j, k)*&
&               temp_ad9/delp(i, j, k+1)
              g_rat_ad(i, k) = 0.0
            END DO
          END DO
          DO k=npz,1,-1
            DO i=ilast,npx,-1
              temp7 = delz(i, j, k)
              temp6 = delp(i, j, k)*pt(i, j, k)
              temp4 = temp6/temp7
              temp5 = -(rgrav*rdgas*temp4)
              temp_ad7 = -(rgrav*rdgas*gama*EXP(gama*LOG(temp5))*pkz_ad(&
&               i, k)/(temp5*temp7))
              pm_ad = -pkz_ad(i, k)
              temp8 = peln(i, k+1, j) - peln(i, k, j)
              temp_ad6 = -(delp(i, j, k)*pm_ad/temp8**2)
              delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*temp_ad7&
&               + pm_ad/temp8
              peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad6
              peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad6
              CALL POPREALARRAY(pkz(i, k))
              pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*temp_ad7
              delz_ad(i, j, k) = delz_ad(i, j, k) - temp4*temp_ad7
              pkz_ad(i, k) = 0.0
            END DO
          END DO
          DO k=npz+1,2,-1
            DO i=ilast,npx,-1
              CALL POPREALARRAY(peln(i, k, j))
              pe_ad(i, k, j) = pe_ad(i, k, j) + peln_ad(i, k, j)/pe(i, k&
&               , j)
              peln_ad(i, k, j) = 0.0
              CALL POPREALARRAY(pe(i, k, j))
              pe_ad(i, k-1, j) = pe_ad(i, k-1, j) + pe_ad(i, k, j)
              delp_ad(i, j, k-1) = delp_ad(i, j, k-1) + pe_ad(i, k, j)
              pe_ad(i, k, j) = 0.0
            END DO
          END DO
          DO i=ilast,npx,-1
            CALL POPREALARRAY(peln(i, 1, j))
            peln_ad(i, 1, j) = 0.0
            CALL POPREALARRAY(pe(i, 1, j))
            pe_ad(i, 1, j) = 0.0
          END DO
          DO k=1,npz,1
            DO i=ilast,npx,-1
              CALL POPREALARRAY(gz(i, j, k))
              gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + gz_ad(i, j, k)
              delz_ad(i, j, k) = delz_ad(i, j, k) - grav*gz_ad(i, j, k)
              gz_ad(i, j, k) = 0.0
            END DO
          END DO
          DO i=ilast,npx,-1
            CALL POPREALARRAY(gz(i, j, npz+1))
            gz_ad(i, j, npz+1) = 0.0
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=jlast,jfirst,-1
          CALL POPCONTROL(1,branch)
          IF (branch .NE. 0) THEN
            DO k=npz+1,2,-1
              DO i=0,ifirst,-1
                CALL POPREALARRAY(pk3(i, j, k))
                pe_ad(i, k, j) = pe_ad(i, k, j) + kappa*EXP(kappa*LOG(pe&
&                 (i, k, j)))*pk3_ad(i, j, k)/pe(i, k, j)
                pk3_ad(i, j, k) = 0.0
              END DO
            END DO
            DO i=0,ifirst,-1
              CALL POPREALARRAY(pk3(i, j, 1))
              pk3_ad(i, j, 1) = 0.0
            END DO
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            DO k=1,npz+1,1
              DO i=0,ifirst,-1
                CALL POPREALARRAY(pkc(i, j, k))
                pe_ad(i, k, j) = pe_ad(i, k, j) + pkc_ad(i, j, k)
              END DO
            END DO
          END IF
        END DO
        DO j=jlast,jfirst,-1
          DO k=2,npz,1
            DO i=0,ifirst,-1
              CALL POPREALARRAY(pkc(i, j, k))
              gam_ad(i, k) = gam_ad(i, k) - pkc(i, j, k+1)*pkc_ad(i, j, &
&               k)
              pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) - gam(i, k)*pkc_ad(i&
&               , j, k)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=0,ifirst,-1
              CALL POPREALARRAY(pkc(i, j, k+1))
              temp_ad4 = pkc_ad(i, j, k+1)/bet(i)
              dd_ad(i, k) = dd_ad(i, k) + temp_ad4
              pkc_ad(i, j, k) = pkc_ad(i, j, k) - temp_ad4
              bet_ad(i) = bet_ad(i) - (dd(i, k)-pkc(i, j, k))*temp_ad4/&
&               bet(i)
              pkc_ad(i, j, k+1) = 0.0
              CALL POPREALARRAY(bet(i))
              bb_ad(i, k) = bb_ad(i, k) + bet_ad(i)
              gam_ad(i, k) = gam_ad(i, k) - bet_ad(i)
              CALL POPREALARRAY(gam(i, k))
              temp_ad5 = gam_ad(i, k)/bet(i)
              bet_ad(i) = -(g_rat(i, k-1)*temp_ad5/bet(i))
              g_rat_ad(i, k-1) = g_rat_ad(i, k-1) + temp_ad5
              gam_ad(i, k) = 0.0
            END DO
          END DO
          DO i=0,ifirst,-1
            CALL POPREALARRAY(dd(i, npz))
            pkz_ad(i, npz) = pkz_ad(i, npz) + 3.*dd_ad(i, npz)
            dd_ad(i, npz) = 0.0
            bb_ad(i, npz) = 0.0
            CALL POPREALARRAY(pkc(i, j, 2))
            temp_ad3 = pkc_ad(i, j, 2)/bet(i)
            dd_ad(i, 1) = dd_ad(i, 1) + temp_ad3
            bet_ad(i) = bet_ad(i) - dd(i, 1)*temp_ad3/bet(i)
            pkc_ad(i, j, 2) = 0.0
            CALL POPREALARRAY(pkc(i, j, 1))
            pkc_ad(i, j, 1) = 0.0
            CALL POPREALARRAY(bet(i))
            bb_ad(i, 1) = bb_ad(i, 1) + bet_ad(i)
            bet_ad(i) = 0.0
          END DO
          DO k=npz-1,1,-1
            DO i=0,ifirst,-1
              CALL POPREALARRAY(dd(i, k))
              temp_ad1 = 3.*dd_ad(i, k)
              pkz_ad(i, k) = pkz_ad(i, k) + temp_ad1
              g_rat_ad(i, k) = g_rat_ad(i, k) + 2.*bb_ad(i, k) + pkz(i, &
&               k+1)*temp_ad1
              pkz_ad(i, k+1) = pkz_ad(i, k+1) + g_rat(i, k)*temp_ad1
              dd_ad(i, k) = 0.0
              bb_ad(i, k) = 0.0
              CALL POPREALARRAY(g_rat(i, k))
              temp_ad2 = g_rat_ad(i, k)/delp(i, j, k+1)
              delp_ad(i, j, k) = delp_ad(i, j, k) + temp_ad2
              delp_ad(i, j, k+1) = delp_ad(i, j, k+1) - delp(i, j, k)*&
&               temp_ad2/delp(i, j, k+1)
              g_rat_ad(i, k) = 0.0
            END DO
          END DO
          DO k=npz,1,-1
            DO i=0,ifirst,-1
              temp2 = delz(i, j, k)
              temp1 = delp(i, j, k)*pt(i, j, k)
              temp = temp1/temp2
              temp0 = -(rgrav*rdgas*temp)
              temp_ad0 = -(rgrav*rdgas*gama*EXP(gama*LOG(temp0))*pkz_ad(&
&               i, k)/(temp0*temp2))
              pm_ad = -pkz_ad(i, k)
              temp3 = peln(i, k+1, j) - peln(i, k, j)
              temp_ad = -(delp(i, j, k)*pm_ad/temp3**2)
              delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*temp_ad0&
&               + pm_ad/temp3
              peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad
              peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad
              CALL POPREALARRAY(pkz(i, k))
              pt_ad(i, j, k) = pt_ad(i, j, k) + delp(i, j, k)*temp_ad0
              delz_ad(i, j, k) = delz_ad(i, j, k) - temp*temp_ad0
              pkz_ad(i, k) = 0.0
            END DO
          END DO
          DO k=npz+1,2,-1
            DO i=0,ifirst,-1
              pe_ad(i, k, j) = pe_ad(i, k, j) + peln_ad(i, k, j)/pe(i, k&
&               , j)
              peln_ad(i, k, j) = 0.0
              CALL POPREALARRAY(pe(i, k, j))
              pe_ad(i, k-1, j) = pe_ad(i, k-1, j) + pe_ad(i, k, j)
              delp_ad(i, j, k-1) = delp_ad(i, j, k-1) + pe_ad(i, k, j)
              pe_ad(i, k, j) = 0.0
            END DO
          END DO
          DO i=0,ifirst,-1
            peln_ad(i, 1, j) = 0.0
            pe_ad(i, 1, j) = 0.0
          END DO
          DO k=1,npz,1
            DO i=0,ifirst,-1
              CALL POPREALARRAY(gz(i, j, k))
              gz_ad(i, j, k+1) = gz_ad(i, j, k+1) + gz_ad(i, j, k)
              delz_ad(i, j, k) = delz_ad(i, j, k) - grav*gz_ad(i, j, k)
              gz_ad(i, j, k) = 0.0
            END DO
          END DO
          DO i=0,ifirst,-1
            CALL POPREALARRAY(gz(i, j, npz+1))
            gz_ad(i, j, npz+1) = 0.0
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE NEST_HALO_NH_BWD
  SUBROUTINE NEST_HALO_NH(ptop, grav, kappa, cp, delp, delz, pt, phis, &
&   pkc, gz, pk3, npx, npy, npz, nested, pkc_pertn, computepk3, fullhalo&
&   , bd)
    IMPLICIT NONE
!INPUT: delp, delz, pt
!OUTPUT: gz, pkc, pk3 (optional)
    INTEGER, INTENT(IN) :: npx, npy, npz
    LOGICAL, INTENT(IN) :: pkc_pertn, computepk3, fullhalo, nested
    REAL, INTENT(IN) :: ptop, kappa, cp, grav
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: phis(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(IN) :: pt&
&   , delp, delz
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(INOUT) &
&   :: gz, pkc, pk3
    INTEGER :: i, j, k
!'gamma'
    REAL :: gama
    REAL :: ptk, rgrav, rkap, peln1, rdg
    REAL, DIMENSION(bd%isd:bd%ied, npz+1, bd%jsd:bd%jed) :: pe, peln
    REAL, DIMENSION(bd%isd:bd%ied, npz) :: gam, bb, dd, pkz
    REAL, DIMENSION(bd%isd:bd%ied, npz-1) :: g_rat
    REAL, DIMENSION(bd%isd:bd%ied) :: bet
    REAL :: pm
    INTEGER :: ifirst, ilast, jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
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
    IF (.NOT.nested) THEN
      RETURN
    ELSE
      ifirst = isd
      jfirst = jsd
      ilast = ied
      jlast = jed
!Remember we want to compute these in the HALO. Note also this routine
!requires an appropriate
      rgrav = 1./grav
      gama = 1./(1.-kappa)
      ptk = ptop**kappa
      rkap = 1./kappa
      peln1 = LOG(ptop)
      rdg = -(rdgas*rgrav)
!NOTE: Compiler does NOT like this sort of nested-grid BC code. Is it trying to do some ugly optimization?
      IF (is .EQ. 1) THEN
        DO j=jfirst,jlast
!GZ
          DO i=ifirst,0
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=ifirst,0
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=ifirst,0
            pe(i, 1, j) = ptop
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=ifirst,0
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=ifirst,0
!Full p
              pkz(i, k) = EXP(gama*LOG(-(delp(i, j, k)*rgrav/delz(i, j, &
&               k)*rdgas*pt(i, j, k))))
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!pressure solver
          DO k=1,npz-1
            DO i=ifirst,0
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=ifirst,0
            bet(i) = bb(i, 1)
            pkc(i, j, 1) = 0.
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb(i, npz) = 2.
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=ifirst,0
              gam(i, k) = g_rat(i, k-1)/bet(i)
              bet(i) = bb(i, k) - gam(i, k)
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ifirst,0
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=jfirst,jlast
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=ifirst,0
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
          END IF
!pk3 if necessary; doesn't require condenstate loading calculation
          IF (computepk3) THEN
            DO i=ifirst,0
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=ifirst,0
                pk3(i, j, k) = EXP(kappa*LOG(pe(i, k, j)))
              END DO
            END DO
          END IF
        END DO
      END IF
      IF (ie .EQ. npx - 1) THEN
        DO j=jfirst,jlast
!GZ
          DO i=npx,ilast
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=npx,ilast
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=npx,ilast
            pe(i, 1, j) = ptop
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=npx,ilast
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=npx,ilast
!Full p
              pkz(i, k) = EXP(gama*LOG(-(delp(i, j, k)*rgrav/delz(i, j, &
&               k)*rdgas*pt(i, j, k))))
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!pressure solver
          DO k=1,npz-1
            DO i=npx,ilast
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=npx,ilast
            bet(i) = bb(i, 1)
            pkc(i, j, 1) = 0.
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb(i, npz) = 2.
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=npx,ilast
              gam(i, k) = g_rat(i, k-1)/bet(i)
              bet(i) = bb(i, k) - gam(i, k)
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=npx,ilast
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=jfirst,jlast
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=npx,ilast
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
          END IF
!pk3 if necessary
          IF (computepk3) THEN
            DO i=npx,ilast
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=npx,ilast
                pk3(i, j, k) = EXP(kappa*LOG(pe(i, k, j)))
              END DO
            END DO
          END IF
        END DO
      END IF
      IF (js .EQ. 1) THEN
        DO j=jfirst,0
!GZ
          DO i=ifirst,ilast
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=ifirst,ilast
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=ifirst,ilast
            pe(i, 1, j) = ptop
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=ifirst,ilast
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=ifirst,ilast
!Full p
              pkz(i, k) = EXP(gama*LOG(-(delp(i, j, k)*rgrav/delz(i, j, &
&               k)*rdgas*pt(i, j, k))))
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!pressure solver
          DO k=1,npz-1
            DO i=ifirst,ilast
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=ifirst,ilast
            bet(i) = bb(i, 1)
            pkc(i, j, 1) = 0.
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb(i, npz) = 2.
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=ifirst,ilast
              gam(i, k) = g_rat(i, k-1)/bet(i)
              bet(i) = bb(i, k) - gam(i, k)
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ifirst,ilast
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=jfirst,0
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=ifirst,ilast
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
          END IF
!pk3 if necessary
          IF (computepk3) THEN
            DO i=ifirst,ilast
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=ifirst,ilast
                pk3(i, j, k) = EXP(kappa*LOG(pe(i, k, j)))
              END DO
            END DO
          END IF
        END DO
      END IF
      IF (je .EQ. npy - 1) THEN
        DO j=npy,jlast
!GZ
          DO i=ifirst,ilast
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=ifirst,ilast
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=ifirst,ilast
            pe(i, 1, j) = ptop
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=ifirst,ilast
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=ifirst,ilast
!Full p
              pkz(i, k) = EXP(gama*LOG(-(delp(i, j, k)*rgrav/delz(i, j, &
&               k)*rdgas*pt(i, j, k))))
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!Reversible interpolation on layer NH pressure perturbation
!                 to recover  lastge NH pressure perturbation
          DO k=1,npz-1
            DO i=ifirst,ilast
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=ifirst,ilast
            bet(i) = bb(i, 1)
            pkc(i, j, 1) = 0.
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb(i, npz) = 2.
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=ifirst,ilast
              gam(i, k) = g_rat(i, k-1)/bet(i)
              bet(i) = bb(i, k) - gam(i, k)
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ifirst,ilast
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=npy,jlast
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=ifirst,ilast
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
          END IF
!pk3 if necessary
          IF (computepk3) THEN
            DO i=ifirst,ilast
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=ifirst,ilast
                pk3(i, j, k) = EXP(kappa*LOG(pe(i, k, j)))
              END DO
            END DO
          END IF
        END DO
      END IF
    END IF
  END SUBROUTINE NEST_HALO_NH
end module nh_utils_adm_mod
