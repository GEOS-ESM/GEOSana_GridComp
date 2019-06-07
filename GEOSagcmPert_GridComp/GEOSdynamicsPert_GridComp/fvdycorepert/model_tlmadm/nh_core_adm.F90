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
module nh_core_adm_mod
! Developer: S.-J. Lin, NOAA/GFDL
! To do list:
! include moisture effect in pt
!------------------------------
   use constants_mod,     only: rdgas, cp_air, grav
   use tp_core_adm_mod,   only: fv_tp_2d
   use tp_core_adm_mod,   only: fv_tp_2d_fwd, fv_tp_2d_bwd
   use nh_utils_adm_mod,  only: update_dz_c, update_dz_d, nest_halo_nh
   use nh_utils_adm_mod,  only: sim_solver, sim1_solver, sim3_solver
   use nh_utils_adm_mod,  only: sim3p0_solver, rim_2d
   use nh_utils_adm_mod,  only: Riem_Solver_c
   use nh_utils_adm_mod,  only: update_dz_c_fwd, update_dz_d_fwd, nest_halo_nh_fwd
   use nh_utils_adm_mod,  only: sim_solver_fwd, sim1_solver_fwd, sim3_solver_fwd
   use nh_utils_adm_mod,  only: sim3p0_solver_fwd, rim_2d_fwd
   use nh_utils_adm_mod,  only: Riem_Solver_c_fwd
   use nh_utils_adm_mod,  only: update_dz_c_bwd, update_dz_d_bwd, nest_halo_nh_bwd
   use nh_utils_adm_mod,  only: sim_solver_bwd, sim1_solver_bwd, sim3_solver_bwd
   use nh_utils_adm_mod,  only: sim3p0_solver_bwd, rim_2d_bwd
   use nh_utils_adm_mod,  only: Riem_Solver_c_bwd

   use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                            pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

   implicit none
   private

   public Riem_Solver3, Riem_Solver_c, update_dz_c, update_dz_d, nest_halo_nh
   public Riem_Solver3_fwd, Riem_Solver_c_fwd, update_dz_c_fwd, update_dz_d_fwd, nest_halo_nh_fwd
   public Riem_Solver3_bwd, Riem_Solver_c_bwd, update_dz_c_bwd, update_dz_d_bwd, nest_halo_nh_bwd
   real, parameter:: r3 = 1./3.

CONTAINS
!  Differentiation of riem_solver3 in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: pk3 ws ppe peln w delp delz
!                pe pk zh pt
!   with respect to varying inputs: pk3 ws ppe peln w delp delz
!                pe pk zh pt
  SUBROUTINE RIEM_SOLVER3_FWD(ms, dt, is, ie, js, je, km, ng, isd, ied, &
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
    REAL, DIMENSION(isd:, jsd:, :), INTENT(IN) :: q_con, cappa
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: delp, pt
    REAL, DIMENSION(isd:ied, jsd:jed, km+1), INTENT(INOUT) :: zh
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: w
    REAL, INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
! ln(pe)
    REAL :: peln(is:ie, km+1, js:je)
    REAL, DIMENSION(isd:ied, jsd:jed, km+1) :: ppe
    REAL :: delz(is-ng:ie+ng, js-ng:je+ng, km)
    REAL :: pk(is:ie, js:je, km+1)
    REAL :: pk3(isd:ied, jsd:jed, km+1)
! Local:
    REAL, DIMENSION(is:ie, km) :: dm, dz2, pm2, w2, gm2, cp2
    REAL, DIMENSION(is:ie, km+1) :: pem, pe2, peln2, peg, pelng
    REAL :: gama, rgrav, ptk, peln1
    INTEGER :: i, j, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC ABS
    REAL :: abs0

    dm = 0.0
    dz2 = 0.0
    pm2 = 0.0
    w2 = 0.0
    gm2 = 0.0
    cp2 = 0.0
    pem = 0.0
    pe2 = 0.0
    peln2 = 0.0
    peg = 0.0
    pelng = 0.0
    gama = 0.0
    rgrav = 0.0
    ptk = 0.0
    peln1 = 0.0
    abs0 = 0.0

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
          CALL PUSHREALARRAY(dm(i, k))
          dm(i, k) = delp(i, j, k)
        END DO
      END DO
      DO i=is,ie
        CALL PUSHREALARRAY(pem(i, 1))
        pem(i, 1) = ptop
        CALL PUSHREALARRAY(peln2(i, 1))
        peln2(i, 1) = peln1
        CALL PUSHREALARRAY(pk3(i, j, 1))
        pk3(i, j, 1) = ptk
      END DO
      DO k=2,km+1
        DO i=is,ie
          CALL PUSHREALARRAY(pem(i, k))
          pem(i, k) = pem(i, k-1) + dm(i, k-1)
          CALL PUSHREALARRAY(peln2(i, k))
          peln2(i, k) = LOG(pem(i, k))
          CALL PUSHREALARRAY(pk3(i, j, k))
          pk3(i, j, k) = EXP(akap*peln2(i, k))
        END DO
      END DO
      DO k=1,km
        DO i=is,ie
          CALL PUSHREALARRAY(pm2(i, k))
          pm2(i, k) = dm(i, k)/(peln2(i, k+1)-peln2(i, k))
          CALL PUSHREALARRAY(dm(i, k))
          dm(i, k) = dm(i, k)*rgrav
          CALL PUSHREALARRAY(dz2(i, k))
          dz2(i, k) = zh(i, j, k+1) - zh(i, j, k)
          CALL PUSHREALARRAY(w2(i, k))
          w2(i, k) = w(i, j, k)
        END DO
      END DO
      IF (a_imp .LT. -0.999) THEN
        CALL SIM3P0_SOLVER_FWD(dt, is, ie, km, rdgas, gama, akap, pe2, &
&                        dm, pem, w2, dz2, pt(is:ie, j, 1:km), ws(is:ie&
&                        , j), p_fac, scale_m)
        CALL PUSHCONTROL(3,4)
      ELSE IF (a_imp .LT. -0.5) THEN
        IF (a_imp .GE. 0.) THEN
          abs0 = a_imp
        ELSE
          abs0 = -a_imp
        END IF
        CALL SIM3_SOLVER_FWD(dt, is, ie, km, rdgas, gama, akap, pe2, dm&
&                      , pem, w2, dz2, pt(is:ie, j, 1:km), ws(is:ie, j)&
&                      , abs0, p_fac, scale_m)
        CALL PUSHCONTROL(3,3)
      ELSE IF (a_imp .LE. 0.5) THEN
        CALL RIM_2D_FWD(ms, dt, is, ie, km, rdgas, gama, gm2, pe2, dm, &
&                 pm2, w2, dz2, pt(is:ie, j, 1:km), ws(is:ie, j), &
&                 .false.)
        CALL PUSHCONTROL(3,2)
      ELSE IF (a_imp .GT. 0.999) THEN
        CALL SIM1_SOLVER_FWD(dt, is, ie, km, rdgas, gama, gm2, cp2, akap&
&                      , pe2, dm, pm2, pem, w2, dz2, pt(is:ie, j, 1:km)&
&                      , ws(is:ie, j), p_fac)
        CALL PUSHCONTROL(3,1)
      ELSE
        CALL SIM_SOLVER_FWD(dt, is, ie, km, rdgas, gama, gm2, cp2, akap&
&                     , pe2, dm, pm2, pem, w2, dz2, pt(is:ie, j, 1:km), &
&                     ws(is:ie, j), a_imp, p_fac, scale_m)
        CALL PUSHCONTROL(3,0)
      END IF
      DO k=1,km
        DO i=is,ie
          CALL PUSHREALARRAY(w(i, j, k))
          w(i, j, k) = w2(i, k)
          CALL PUSHREALARRAY(delz(i, j, k))
          delz(i, j, k) = dz2(i, k)
        END DO
      END DO
      IF (last_call) THEN
        DO k=1,km+1
          DO i=is,ie
            CALL PUSHREALARRAY(peln(i, k, j))
            peln(i, k, j) = peln2(i, k)
            CALL PUSHREALARRAY(pk(i, j, k))
            pk(i, j, k) = pk3(i, j, k)
            CALL PUSHREALARRAY(pe(i, k, j))
            pe(i, k, j) = pem(i, k)
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (fp_out) THEN
        DO k=1,km+1
          DO i=is,ie
            CALL PUSHREALARRAY(ppe(i, j, k))
            ppe(i, j, k) = pe2(i, k) + pem(i, k)
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        DO k=1,km+1
          DO i=is,ie
            CALL PUSHREALARRAY(ppe(i, j, k))
            ppe(i, j, k) = pe2(i, k)
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      END IF
      IF (use_logp) THEN
        DO k=2,km+1
          DO i=is,ie
            CALL PUSHREALARRAY(pk3(i, j, k))
            pk3(i, j, k) = peln2(i, k)
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
      DO i=is,ie
        CALL PUSHREALARRAY(zh(i, j, km+1))
        zh(i, j, km+1) = zs(i, j)
      END DO
      DO k=km,1,-1
        DO i=is,ie
          CALL PUSHREALARRAY(zh(i, j, k))
          zh(i, j, k) = zh(i, j, k+1) - dz2(i, k)
        END DO
      END DO
    END DO
    CALL PUSHREALARRAY(peln2, (ie-is+1)*(km+1))
    CALL PUSHREALARRAY(pem, (ie-is+1)*(km+1))
    CALL PUSHREALARRAY(gama)
    CALL PUSHREALARRAY(pm2, (ie-is+1)*km)
    CALL PUSHREALARRAY(rgrav)
    CALL PUSHREALARRAY(w2, (ie-is+1)*km)
    CALL PUSHREALARRAY(dz2, (ie-is+1)*km)
    CALL PUSHREALARRAY(pe2, (ie-is+1)*(km+1))
    CALL PUSHREALARRAY(dm, (ie-is+1)*km)
  END SUBROUTINE RIEM_SOLVER3_FWD
!  Differentiation of riem_solver3 in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: pk3 ws ppe peln w delp delz
!                pe pk zh pt
!   with respect to varying inputs: pk3 ws ppe peln w delp delz
!                pe pk zh pt
  SUBROUTINE RIEM_SOLVER3_BWD(ms, dt, is, ie, js, je, km, ng, isd, ied, &
&   jsd, jed, akap, cappa, cp, ptop, zs, q_con, w, w_ad, delz, delz_ad, &
&   pt, pt_ad, delp, delp_ad, zh, zh_ad, pe, pe_ad, ppe, ppe_ad, pk3, &
&   pk3_ad, pk, pk_ad, peln, peln_ad, ws, ws_ad, scale_m, p_fac, a_imp, &
&   use_logp, last_call, fp_out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ms, is, ie, js, je, km, ng
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: akap, cp, ptop, p_fac, a_imp, scale_m
    REAL, INTENT(IN) :: zs(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: last_call, use_logp, fp_out
    REAL, INTENT(IN) :: ws(is:ie, js:je)
    REAL :: ws_ad(is:ie, js:je)
    REAL, DIMENSION(isd:, jsd:, :), INTENT(IN) :: q_con, cappa
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: delp, pt
    REAL, DIMENSION(isd:ied, jsd:jed, km) :: delp_ad, pt_ad
    REAL, DIMENSION(isd:ied, jsd:jed, km+1), INTENT(INOUT) :: zh
    REAL, DIMENSION(isd:ied, jsd:jed, km+1), INTENT(INOUT) :: zh_ad
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: w
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: w_ad
    REAL, INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(INOUT) :: pe_ad(is-1:ie+1, km+1, js-1:je+1)
    REAL :: peln(is:ie, km+1, js:je)
    REAL :: peln_ad(is:ie, km+1, js:je)
    REAL, DIMENSION(isd:ied, jsd:jed, km+1) :: ppe
    REAL, DIMENSION(isd:ied, jsd:jed, km+1) :: ppe_ad
    REAL :: delz(is-ng:ie+ng, js-ng:je+ng, km)
    REAL :: delz_ad(is-ng:ie+ng, js-ng:je+ng, km)
    REAL :: pk(is:ie, js:je, km+1)
    REAL :: pk_ad(is:ie, js:je, km+1)
    REAL :: pk3(isd:ied, jsd:jed, km+1)
    REAL :: pk3_ad(isd:ied, jsd:jed, km+1)
    REAL, DIMENSION(is:ie, km) :: dm, dz2, pm2, w2, gm2, cp2
    REAL, DIMENSION(is:ie, km) :: dm_ad, dz2_ad, pm2_ad, w2_ad
    REAL, DIMENSION(is:ie, km+1) :: pem, pe2, peln2, peg, pelng
    REAL, DIMENSION(is:ie, km+1) :: pem_ad, pe2_ad, peln2_ad
    REAL :: gama, rgrav, ptk, peln1
    INTEGER :: i, j, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC ABS
    REAL :: abs0
    REAL :: temp
    REAL :: temp_ad
    INTEGER :: branch

    dm = 0.0
    dz2 = 0.0
    pm2 = 0.0
    w2 = 0.0
    gm2 = 0.0
    cp2 = 0.0
    pem = 0.0
    pe2 = 0.0
    peln2 = 0.0
    peg = 0.0
    pelng = 0.0
    gama = 0.0
    rgrav = 0.0
    ptk = 0.0
    peln1 = 0.0
    abs0 = 0.0
    branch = 0

    CALL POPREALARRAY(dm, (ie-is+1)*km)
    CALL POPREALARRAY(pe2, (ie-is+1)*(km+1))
    CALL POPREALARRAY(dz2, (ie-is+1)*km)
    CALL POPREALARRAY(w2, (ie-is+1)*km)
    CALL POPREALARRAY(rgrav)
    CALL POPREALARRAY(pm2, (ie-is+1)*km)
    CALL POPREALARRAY(gama)
    CALL POPREALARRAY(pem, (ie-is+1)*(km+1))
    CALL POPREALARRAY(peln2, (ie-is+1)*(km+1))
    dm_ad = 0.0
    pe2_ad = 0.0
    dz2_ad = 0.0
    w2_ad = 0.0
    pm2_ad = 0.0
    pem_ad = 0.0
    peln2_ad = 0.0
    DO j=je,js,-1
      DO k=1,km,1
        DO i=ie,is,-1
          CALL POPREALARRAY(zh(i, j, k))
          zh_ad(i, j, k+1) = zh_ad(i, j, k+1) + zh_ad(i, j, k)
          dz2_ad(i, k) = dz2_ad(i, k) - zh_ad(i, j, k)
          zh_ad(i, j, k) = 0.0
        END DO
      END DO
      DO i=ie,is,-1
        CALL POPREALARRAY(zh(i, j, km+1))
        zh_ad(i, j, km+1) = 0.0
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        DO k=km+1,2,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(pk3(i, j, k))
            peln2_ad(i, k) = peln2_ad(i, k) + pk3_ad(i, j, k)
            pk3_ad(i, j, k) = 0.0
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO k=km+1,1,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(ppe(i, j, k))
            pe2_ad(i, k) = pe2_ad(i, k) + ppe_ad(i, j, k)
            pem_ad(i, k) = pem_ad(i, k) + ppe_ad(i, j, k)
            ppe_ad(i, j, k) = 0.0
          END DO
        END DO
      ELSE
        DO k=km+1,1,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(ppe(i, j, k))
            pe2_ad(i, k) = pe2_ad(i, k) + ppe_ad(i, j, k)
            ppe_ad(i, j, k) = 0.0
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO k=km+1,1,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(pe(i, k, j))
            pem_ad(i, k) = pem_ad(i, k) + pe_ad(i, k, j)
            pe_ad(i, k, j) = 0.0
            CALL POPREALARRAY(pk(i, j, k))
            pk3_ad(i, j, k) = pk3_ad(i, j, k) + pk_ad(i, j, k)
            pk_ad(i, j, k) = 0.0
            CALL POPREALARRAY(peln(i, k, j))
            peln2_ad(i, k) = peln2_ad(i, k) + peln_ad(i, k, j)
            peln_ad(i, k, j) = 0.0
          END DO
        END DO
      END IF
      DO k=km,1,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(delz(i, j, k))
          dz2_ad(i, k) = dz2_ad(i, k) + delz_ad(i, j, k)
          delz_ad(i, j, k) = 0.0
          CALL POPREALARRAY(w(i, j, k))
          w2_ad(i, k) = w2_ad(i, k) + w_ad(i, j, k)
          w_ad(i, j, k) = 0.0
        END DO
      END DO
      CALL POPCONTROL(3,branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          CALL SIM_SOLVER_BWD(dt, is, ie, km, rdgas, gama, gm2, cp2, &
&                       akap, pe2, pe2_ad, dm, dm_ad, pm2, pm2_ad, pem, &
&                       pem_ad, w2, w2_ad, dz2, dz2_ad, pt(is:ie, j, 1:&
&                       km), pt_ad(is:ie, j, 1:km), ws(is:ie, j), ws_ad(&
&                       is:ie, j), a_imp, p_fac, scale_m)
        ELSE
          CALL SIM1_SOLVER_BWD(dt, is, ie, km, rdgas, gama, gm2, cp2, &
&                        akap, pe2, pe2_ad, dm, dm_ad, pm2, pm2_ad, pem&
&                        , pem_ad, w2, w2_ad, dz2, dz2_ad, pt(is:ie, j, &
&                        1:km), pt_ad(is:ie, j, 1:km), ws(is:ie, j), &
&                        ws_ad(is:ie, j), p_fac)
        END IF
      ELSE IF (branch .EQ. 2) THEN
        CALL RIM_2D_BWD(ms, dt, is, ie, km, rdgas, gama, gm2, pe2, &
&                 pe2_ad, dm, dm_ad, pm2, pm2_ad, w2, w2_ad, dz2, dz2_ad&
&                 , pt(is:ie, j, 1:km), pt_ad(is:ie, j, 1:km), ws(is:ie&
&                 , j), ws_ad(is:ie, j), .false.)
      ELSE IF (branch .EQ. 3) THEN
        CALL SIM3_SOLVER_BWD(dt, is, ie, km, rdgas, gama, akap, pe2, &
&                      pe2_ad, dm, dm_ad, pem, pem_ad, w2, w2_ad, dz2, &
&                      dz2_ad, pt(is:ie, j, 1:km), pt_ad(is:ie, j, 1:km)&
&                      , ws(is:ie, j), ws_ad(is:ie, j), abs0, p_fac, &
&                      scale_m)
      ELSE
        CALL SIM3P0_SOLVER_BWD(dt, is, ie, km, rdgas, gama, akap, pe2, &
&                        pe2_ad, dm, dm_ad, pem, pem_ad, w2, w2_ad, dz2&
&                        , dz2_ad, pt(is:ie, j, 1:km), pt_ad(is:ie, j, 1&
&                        :km), ws(is:ie, j), ws_ad(is:ie, j), p_fac, &
&                        scale_m)
      END IF
      DO k=km,1,-1
        DO i=ie,is,-1
          temp = peln2(i, k+1) - peln2(i, k)
          CALL POPREALARRAY(w2(i, k))
          w_ad(i, j, k) = w_ad(i, j, k) + w2_ad(i, k)
          w2_ad(i, k) = 0.0
          CALL POPREALARRAY(dz2(i, k))
          zh_ad(i, j, k+1) = zh_ad(i, j, k+1) + dz2_ad(i, k)
          zh_ad(i, j, k) = zh_ad(i, j, k) - dz2_ad(i, k)
          dz2_ad(i, k) = 0.0
          CALL POPREALARRAY(dm(i, k))
          dm_ad(i, k) = pm2_ad(i, k)/temp + rgrav*dm_ad(i, k)
          CALL POPREALARRAY(pm2(i, k))
          temp_ad = -(dm(i, k)*pm2_ad(i, k)/temp**2)
          peln2_ad(i, k+1) = peln2_ad(i, k+1) + temp_ad
          peln2_ad(i, k) = peln2_ad(i, k) - temp_ad
          pm2_ad(i, k) = 0.0
        END DO
      END DO
      DO k=km+1,2,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(pk3(i, j, k))
          peln2_ad(i, k) = peln2_ad(i, k) + EXP(akap*peln2(i, k))*akap*&
&           pk3_ad(i, j, k)
          pk3_ad(i, j, k) = 0.0
          CALL POPREALARRAY(peln2(i, k))
          pem_ad(i, k) = pem_ad(i, k) + peln2_ad(i, k)/pem(i, k)
          peln2_ad(i, k) = 0.0
          CALL POPREALARRAY(pem(i, k))
          pem_ad(i, k-1) = pem_ad(i, k-1) + pem_ad(i, k)
          dm_ad(i, k-1) = dm_ad(i, k-1) + pem_ad(i, k)
          pem_ad(i, k) = 0.0
        END DO
      END DO
      DO i=ie,is,-1
        CALL POPREALARRAY(pk3(i, j, 1))
        pk3_ad(i, j, 1) = 0.0
        CALL POPREALARRAY(peln2(i, 1))
        peln2_ad(i, 1) = 0.0
        CALL POPREALARRAY(pem(i, 1))
        pem_ad(i, 1) = 0.0
      END DO
      DO k=km,1,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(dm(i, k))
          delp_ad(i, j, k) = delp_ad(i, j, k) + dm_ad(i, k)
          dm_ad(i, k) = 0.0
        END DO
      END DO
    END DO
  END SUBROUTINE RIEM_SOLVER3_BWD
  SUBROUTINE RIEM_SOLVER3(ms, dt, is, ie, js, je, km, ng, isd, ied, jsd&
&   , jed, akap, cappa, cp, ptop, zs, q_con, w, delz, pt, delp, zh, pe, &
&   ppe, pk3, pk, peln, ws, scale_m, p_fac, a_imp, use_logp, last_call, &
&   fp_out)
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
    REAL, DIMENSION(isd:, jsd:, :), INTENT(IN) :: q_con, cappa
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
  END SUBROUTINE RIEM_SOLVER3
end module nh_core_adm_mod
