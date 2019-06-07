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
module fv_restart_adm_mod

use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                         pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

implicit none
private

public :: d2c_setup, d2a_setup
public :: d2c_setup_fwd, d2c_setup_bwd

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of d2c_setup in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
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
!   gradient     of useful results: u v uc vc
!   with respect to varying inputs: u v uc vc
  SUBROUTINE D2C_SETUP_FWD(u, v, ua, va, uc, vc, dord4, isd, ied, jsd, &
&   jed, is, ie, js, je, npx, npy, grid_type, nested, se_corner, &
&   sw_corner, ne_corner, nw_corner, rsin_u, rsin_v, cosa_s, rsin2)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: dord4
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, npx, npy&
&   , grid_type
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: va
    REAL, DIMENSION(isd:ied+1, jsd:jed) :: uc
    REAL, DIMENSION(isd:ied, jsd:jed+1) :: vc
    LOGICAL, INTENT(IN) :: nested, se_corner, sw_corner, ne_corner, &
&   nw_corner
    REAL, INTENT(IN) :: rsin_u(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rsin_v(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: cosa_s(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rsin2(isd:ied, jsd:jed)
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, PARAMETER :: t11=27./28., t12=-(13./28.), t13=3./7., t14=6./7.&
&   , t15=3./28.
    REAL, PARAMETER :: a1=0.5625
    REAL, PARAMETER :: a2=-0.0625
    REAL, PARAMETER :: c1=-(2./14.)
    REAL, PARAMETER :: c2=11./14.
    REAL, PARAMETER :: c3=5./14.
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    INTEGER :: ad_from
    INTEGER :: ad_from0

    utmp = 0.0
    vtmp = 0.0
    npt = 0
    ifirst = 0
    ilast = 0
    id = 0
    max1 = 0
    max2 = 0
    max3 = 0
    max4 = 0
    max5 = 0
    max6 = 0
    min1 = 0
    min2 = 0
    min3 = 0
    min4 = 0
    min5 = 0
    min6 = 0
    ad_from = 0
    ad_from0 = 0

    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
    IF (nested) THEN
      CALL PUSHCONTROL(2,0)
    ELSE
!----------
! Interior:
!----------
      IF (npt .LT. js - 1) THEN
        max1 = js - 1
      ELSE
        max1 = npt
      END IF
      IF (npy - npt .GT. je + 1) THEN
        min1 = je + 1
      ELSE
        min1 = npy - npt
      END IF
      DO j=max1,min1
        IF (npt .LT. isd) THEN
          max2 = isd
        ELSE
          max2 = npt
        END IF
        IF (npx - npt .GT. ied) THEN
          min2 = ied
        ELSE
          min2 = npx - npt
        END IF
        ad_from = max2
        i = min2 + 1
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from)
      END DO
      IF (npt .LT. jsd) THEN
        max3 = jsd
      ELSE
        max3 = npt
      END IF
      IF (npy - npt .GT. jed) THEN
        min3 = jed
      ELSE
        min3 = npy - npt
      END IF
      DO j=max3,min3
        IF (npt .LT. is - 1) THEN
          max4 = is - 1
        ELSE
          max4 = npt
        END IF
        IF (npx - npt .GT. ie + 1) THEN
          min4 = ie + 1
        ELSE
          min4 = npx - npt
        END IF
        ad_from0 = max4
        i = min4 + 1
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from0)
      END DO
!----------
! edges:
!----------
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (is .EQ. 1 .OR. isd .LT. npt) THEN
          IF (npt .LT. jsd) THEN
            max5 = jsd
          ELSE
            max5 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min5 = jed
          ELSE
            min5 = npy - npt
          END IF
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
          IF (npt .LT. jsd) THEN
            max6 = jsd
          ELSE
            max6 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min6 = jed
          ELSE
            min6 = npy - npt
          END IF
          CALL PUSHCONTROL(2,3)
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
      ELSE
        CALL PUSHCONTROL(2,1)
      END IF
    END IF
! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
    IF (sw_corner) THEN
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (se_corner) THEN
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (ne_corner) THEN
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (nw_corner) THEN
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      IF (3 .LT. is - 1) THEN
        ifirst = is - 1
      ELSE
        ifirst = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        CALL PUSHCONTROL(1,1)
        ilast = ie + 2
      ELSE
        CALL PUSHCONTROL(1,1)
        ilast = npx - 2
      END IF
    ELSE
      CALL PUSHCONTROL(1,0)
      ifirst = is - 1
      ilast = ie + 2
    END IF
    IF (grid_type .LT. 3) THEN
! Xdir:
      IF (is .EQ. 1 .AND. (.NOT.nested)) THEN
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ie + 1 .EQ. npx .AND. (.NOT.nested)) THEN
        CALL PUSHCONTROL(2,0)
      ELSE
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHCONTROL(2,2)
    END IF
!------
! Ydir:
!------
    IF (sw_corner) THEN
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (nw_corner) THEN
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (se_corner) THEN
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (ne_corner) THEN
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1 .AND. (.NOT.nested)) THEN
          CALL PUSHCONTROL(3,4)
        ELSE IF ((j .EQ. 0 .OR. j .EQ. npy - 1) .AND. (.NOT.nested)) &
&       THEN
          CALL PUSHCONTROL(3,3)
        ELSE IF ((j .EQ. 2 .OR. j .EQ. npy + 1) .AND. (.NOT.nested)) &
&       THEN
          CALL PUSHCONTROL(3,2)
        ELSE IF (j .EQ. npy .AND. (.NOT.nested)) THEN
          CALL PUSHCONTROL(3,1)
        ELSE
          CALL PUSHCONTROL(3,0)
        END IF
      END DO
      CALL PUSHINTEGER(npt)
      CALL PUSHINTEGER(ifirst)
      CALL PUSHINTEGER(min6)
      CALL PUSHINTEGER(min5)
      CALL PUSHINTEGER(min3)
      CALL PUSHINTEGER(min1)
      CALL PUSHINTEGER(id)
      CALL PUSHINTEGER(ilast)
      CALL PUSHINTEGER(max6)
      CALL PUSHINTEGER(max5)
      CALL PUSHINTEGER(max3)
      CALL PUSHINTEGER(max1)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHINTEGER(npt)
      CALL PUSHINTEGER(ifirst)
      CALL PUSHINTEGER(min6)
      CALL PUSHINTEGER(min5)
      CALL PUSHINTEGER(min3)
      CALL PUSHINTEGER(min1)
      CALL PUSHINTEGER(id)
      CALL PUSHINTEGER(ilast)
      CALL PUSHINTEGER(max6)
      CALL PUSHINTEGER(max5)
      CALL PUSHINTEGER(max3)
      CALL PUSHINTEGER(max1)
      CALL PUSHCONTROL(1,1)
    END IF
  END SUBROUTINE D2C_SETUP_FWD
!  Differentiation of d2c_setup in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_m
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
!   gradient     of useful results: u v uc vc
!   with respect to varying inputs: u v uc vc
  SUBROUTINE D2C_SETUP_BWD(u, u_ad, v, v_ad, ua, va, uc, uc_ad, vc, &
&   vc_ad, dord4, isd, ied, jsd, jed, is, ie, js, je, npx, npy, &
&   grid_type, nested, se_corner, sw_corner, ne_corner, nw_corner, &
&   rsin_u, rsin_v, cosa_s, rsin2)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: dord4
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, npx, npy&
&   , grid_type
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL :: u_ad(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL :: v_ad(isd:ied+1, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: va
    REAL, DIMENSION(isd:ied+1, jsd:jed) :: uc
    REAL, DIMENSION(isd:ied+1, jsd:jed) :: uc_ad
    REAL, DIMENSION(isd:ied, jsd:jed+1) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed+1) :: vc_ad
    LOGICAL, INTENT(IN) :: nested, se_corner, sw_corner, ne_corner, &
&   nw_corner
    REAL, INTENT(IN) :: rsin_u(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rsin_v(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: cosa_s(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rsin2(isd:ied, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp_ad, vtmp_ad
    REAL, PARAMETER :: t11=27./28., t12=-(13./28.), t13=3./7., t14=6./7.&
&   , t15=3./28.
    REAL, PARAMETER :: a1=0.5625
    REAL, PARAMETER :: a2=-0.0625
    REAL, PARAMETER :: c1=-(2./14.)
    REAL, PARAMETER :: c2=11./14.
    REAL, PARAMETER :: c3=5./14.
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: branch

    utmp = 0.0
    vtmp = 0.0
    npt = 0
    ifirst = 0
    ilast = 0
    id = 0
    max1 = 0
    max2 = 0
    max3 = 0
    max4 = 0
    max5 = 0
    max6 = 0
    min1 = 0
    min2 = 0
    min3 = 0
    min4 = 0
    min5 = 0
    min6 = 0
    ad_from = 0
    ad_to = 0
    ad_from0 = 0
    ad_to0 = 0
    branch = 0

    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(max1)
      CALL POPINTEGER(max3)
      CALL POPINTEGER(max5)
      CALL POPINTEGER(max6)
      CALL POPINTEGER(ilast)
      CALL POPINTEGER(id)
      CALL POPINTEGER(min1)
      CALL POPINTEGER(min3)
      CALL POPINTEGER(min5)
      CALL POPINTEGER(min6)
      CALL POPINTEGER(ifirst)
      CALL POPINTEGER(npt)
      vtmp_ad = 0.0
      DO j=je+2,js-1,-1
        CALL POPCONTROL(3,branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            DO i=ie+1,is-1,-1
              vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + a2*vc_ad(i, j)
              vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + a2*vc_ad(i, j)
              vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + a1*vc_ad(i, j)
              vtmp_ad(i, j) = vtmp_ad(i, j) + a1*vc_ad(i, j)
              vc_ad(i, j) = 0.0
            END DO
          ELSE
            DO i=ie+1,is-1,-1
              temp_ad2 = rsin_v(i, npy)*vc_ad(i, npy)
              vtmp_ad(i, npy-1) = vtmp_ad(i, npy-1) + t14*temp_ad2
              vtmp_ad(i, npy) = vtmp_ad(i, npy) + t14*temp_ad2
              vtmp_ad(i, npy-2) = vtmp_ad(i, npy-2) + t12*temp_ad2
              vtmp_ad(i, npy+1) = vtmp_ad(i, npy+1) + t12*temp_ad2
              vtmp_ad(i, npy-3) = vtmp_ad(i, npy-3) + t15*temp_ad2
              vtmp_ad(i, npy+2) = vtmp_ad(i, npy+2) + t15*temp_ad2
              vc_ad(i, npy) = 0.0
            END DO
          END IF
        ELSE IF (branch .EQ. 2) THEN
          DO i=ie+1,is-1,-1
            vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + c1*vc_ad(i, j)
            vtmp_ad(i, j) = vtmp_ad(i, j) + c2*vc_ad(i, j)
            vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + c3*vc_ad(i, j)
            vc_ad(i, j) = 0.0
          END DO
        ELSE IF (branch .EQ. 3) THEN
          DO i=ie+1,is-1,-1
            vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + c1*vc_ad(i, j)
            vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + c2*vc_ad(i, j)
            vtmp_ad(i, j) = vtmp_ad(i, j) + c3*vc_ad(i, j)
            vc_ad(i, j) = 0.0
          END DO
        ELSE
          DO i=ie+1,is-1,-1
            temp_ad1 = rsin_v(i, 1)*vc_ad(i, 1)
            vtmp_ad(i, 0) = vtmp_ad(i, 0) + t14*temp_ad1
            vtmp_ad(i, 1) = vtmp_ad(i, 1) + t14*temp_ad1
            vtmp_ad(i, -1) = vtmp_ad(i, -1) + t12*temp_ad1
            vtmp_ad(i, 2) = vtmp_ad(i, 2) + t12*temp_ad1
            vtmp_ad(i, -2) = vtmp_ad(i, -2) + t15*temp_ad1
            vtmp_ad(i, 3) = vtmp_ad(i, 3) + t15*temp_ad1
            vc_ad(i, 1) = 0.0
          END DO
        END IF
      END DO
    ELSE
      CALL POPINTEGER(max1)
      CALL POPINTEGER(max3)
      CALL POPINTEGER(max5)
      CALL POPINTEGER(max6)
      CALL POPINTEGER(ilast)
      CALL POPINTEGER(id)
      CALL POPINTEGER(min1)
      CALL POPINTEGER(min3)
      CALL POPINTEGER(min5)
      CALL POPINTEGER(min6)
      CALL POPINTEGER(ifirst)
      CALL POPINTEGER(npt)
      vtmp_ad = 0.0
      DO j=je+2,js-1,-1
        DO i=ie+1,is-1,-1
          vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + a2*vc_ad(i, j)
          vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + a2*vc_ad(i, j)
          vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + a1*vc_ad(i, j)
          vtmp_ad(i, j) = vtmp_ad(i, j) + a1*vc_ad(i, j)
          vc_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      utmp_ad = 0.0
      DO j=2,0,-1
        utmp_ad(ie-j, npy) = utmp_ad(ie-j, npy) - vtmp_ad(npx, npy+j)
        vtmp_ad(npx, npy+j) = 0.0
      END DO
    ELSE
      utmp_ad = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=0,-2,-1
        utmp_ad(ie+j, 0) = utmp_ad(ie+j, 0) + vtmp_ad(npx, j)
        vtmp_ad(npx, j) = 0.0
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=2,0,-1
        utmp_ad(j+1, npy) = utmp_ad(j+1, npy) + vtmp_ad(0, npy+j)
        vtmp_ad(0, npy+j) = 0.0
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=0,-2,-1
        utmp_ad(1-j, 0) = utmp_ad(1-j, 0) - vtmp_ad(0, j)
        vtmp_ad(0, j) = 0.0
      END DO
    END IF
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      DO j=je+1,js-1,-1
        utmp_ad(npx, j) = utmp_ad(npx, j) + c3*uc_ad(npx+1, j)
        utmp_ad(npx+1, j) = utmp_ad(npx+1, j) + c2*uc_ad(npx+1, j)
        utmp_ad(npx+2, j) = utmp_ad(npx+2, j) + c1*uc_ad(npx+1, j)
        uc_ad(npx+1, j) = 0.0
        temp_ad0 = rsin_u(npx, j)*uc_ad(npx, j)
        utmp_ad(npx-1, j) = utmp_ad(npx-1, j) + t14*temp_ad0
        utmp_ad(npx, j) = utmp_ad(npx, j) + t14*temp_ad0
        utmp_ad(npx-2, j) = utmp_ad(npx-2, j) + t12*temp_ad0
        utmp_ad(npx+1, j) = utmp_ad(npx+1, j) + t12*temp_ad0
        utmp_ad(npx-3, j) = utmp_ad(npx-3, j) + t15*temp_ad0
        utmp_ad(npx+2, j) = utmp_ad(npx+2, j) + t15*temp_ad0
        uc_ad(npx, j) = 0.0
        utmp_ad(npx-3, j) = utmp_ad(npx-3, j) + c1*uc_ad(npx-1, j)
        utmp_ad(npx-2, j) = utmp_ad(npx-2, j) + c2*uc_ad(npx-1, j)
        utmp_ad(npx-1, j) = utmp_ad(npx-1, j) + c3*uc_ad(npx-1, j)
        uc_ad(npx-1, j) = 0.0
      END DO
    ELSE IF (branch .NE. 1) THEN
      GOTO 100
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=je+1,js-1,-1
        utmp_ad(3, j) = utmp_ad(3, j) + c1*uc_ad(2, j)
        utmp_ad(2, j) = utmp_ad(2, j) + c2*uc_ad(2, j)
        utmp_ad(1, j) = utmp_ad(1, j) + c3*uc_ad(2, j)
        uc_ad(2, j) = 0.0
        temp_ad = rsin_u(1, j)*uc_ad(1, j)
        utmp_ad(0, j) = utmp_ad(0, j) + t14*temp_ad
        utmp_ad(1, j) = utmp_ad(1, j) + t14*temp_ad
        utmp_ad(-1, j) = utmp_ad(-1, j) + t12*temp_ad
        utmp_ad(2, j) = utmp_ad(2, j) + t12*temp_ad
        utmp_ad(-2, j) = utmp_ad(-2, j) + t15*temp_ad
        utmp_ad(3, j) = utmp_ad(3, j) + t15*temp_ad
        uc_ad(1, j) = 0.0
        utmp_ad(-2, j) = utmp_ad(-2, j) + c1*uc_ad(0, j)
        utmp_ad(-1, j) = utmp_ad(-1, j) + c2*uc_ad(0, j)
        utmp_ad(0, j) = utmp_ad(0, j) + c3*uc_ad(0, j)
        uc_ad(0, j) = 0.0
      END DO
    END IF
 100 DO j=je+1,js-1,-1
      DO i=ilast,ifirst,-1
        utmp_ad(i-1, j) = utmp_ad(i-1, j) + a1*uc_ad(i, j)
        utmp_ad(i, j) = utmp_ad(i, j) + a1*uc_ad(i, j)
        utmp_ad(i-2, j) = utmp_ad(i-2, j) + a2*uc_ad(i, j)
        utmp_ad(i+1, j) = utmp_ad(i+1, j) + a2*uc_ad(i, j)
        uc_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO i=0,-2,-1
        vtmp_ad(0, je+i) = vtmp_ad(0, je+i) + utmp_ad(i, npy)
        utmp_ad(i, npy) = 0.0
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO i=2,0,-1
        vtmp_ad(npx, je-i) = vtmp_ad(npx, je-i) - utmp_ad(npx+i, npy)
        utmp_ad(npx+i, npy) = 0.0
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO i=2,0,-1
        vtmp_ad(npx, i+1) = vtmp_ad(npx, i+1) + utmp_ad(npx+i, 0)
        utmp_ad(npx+i, 0) = 0.0
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO i=0,-2,-1
        vtmp_ad(0, 1-i) = vtmp_ad(0, 1-i) - utmp_ad(i, 0)
        utmp_ad(i, 0) = 0.0
      END DO
    END IF
    CALL POPCONTROL(2,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) THEN
        DO j=jed,jsd,-1
          v_ad(ied, j) = v_ad(ied, j) + 0.5*vtmp_ad(ied, j)
          v_ad(ied+1, j) = v_ad(ied+1, j) + 0.5*vtmp_ad(ied, j)
          vtmp_ad(ied, j) = 0.0
          v_ad(isd, j) = v_ad(isd, j) + 0.5*vtmp_ad(isd, j)
          v_ad(isd+1, j) = v_ad(isd+1, j) + 0.5*vtmp_ad(isd, j)
          vtmp_ad(isd, j) = 0.0
          DO i=ied-1,isd+1,-1
            v_ad(i-1, j) = v_ad(i-1, j) + a2*vtmp_ad(i, j)
            v_ad(i+2, j) = v_ad(i+2, j) + a2*vtmp_ad(i, j)
            v_ad(i, j) = v_ad(i, j) + a1*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + a1*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0
          END DO
        END DO
        DO i=ied,isd,-1
          u_ad(i, jed) = u_ad(i, jed) + 0.5*utmp_ad(i, jed)
          u_ad(i, jed+1) = u_ad(i, jed+1) + 0.5*utmp_ad(i, jed)
          utmp_ad(i, jed) = 0.0
          u_ad(i, jsd) = u_ad(i, jsd) + 0.5*utmp_ad(i, jsd)
          u_ad(i, jsd+1) = u_ad(i, jsd+1) + 0.5*utmp_ad(i, jsd)
          utmp_ad(i, jsd) = 0.0
        END DO
        DO j=jed-1,jsd+1,-1
          DO i=ied,isd,-1
            u_ad(i, j-1) = u_ad(i, j-1) + a2*utmp_ad(i, j)
            u_ad(i, j+2) = u_ad(i, j+2) + a2*utmp_ad(i, j)
            u_ad(i, j) = u_ad(i, j) + a1*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + a1*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0
          END DO
        END DO
        GOTO 110
      END IF
    ELSE
      IF (branch .NE. 2) THEN
        DO j=min6,max6,-1
          DO i=ied,npx-npt+1,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=min5,max5,-1
          DO i=npt-1,isd,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=jed,npy-npt+1,-1
          DO i=ied,isd,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=npt-1,jsd,-1
          DO i=ied,isd,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0
          END DO
        END DO
      END IF
    END IF
    DO j=min3,max3,-1
      CALL POPINTEGER(ad_from0)
      CALL POPINTEGER(ad_to0)
      DO i=ad_to0,ad_from0,-1
        v_ad(i-1, j) = v_ad(i-1, j) + a2*vtmp_ad(i, j)
        v_ad(i+2, j) = v_ad(i+2, j) + a2*vtmp_ad(i, j)
        v_ad(i, j) = v_ad(i, j) + a1*vtmp_ad(i, j)
        v_ad(i+1, j) = v_ad(i+1, j) + a1*vtmp_ad(i, j)
        vtmp_ad(i, j) = 0.0
      END DO
    END DO
    DO j=min1,max1,-1
      CALL POPINTEGER(ad_from)
      CALL POPINTEGER(ad_to)
      DO i=ad_to,ad_from,-1
        u_ad(i, j-1) = u_ad(i, j-1) + a2*utmp_ad(i, j)
        u_ad(i, j+2) = u_ad(i, j+2) + a2*utmp_ad(i, j)
        u_ad(i, j) = u_ad(i, j) + a1*utmp_ad(i, j)
        u_ad(i, j+1) = u_ad(i, j+1) + a1*utmp_ad(i, j)
        utmp_ad(i, j) = 0.0
      END DO
    END DO
 110 CONTINUE
  END SUBROUTINE D2C_SETUP_BWD
  SUBROUTINE D2C_SETUP(u, v, ua, va, uc, vc, dord4, isd, ied, jsd, jed, &
&   is, ie, js, je, npx, npy, grid_type, nested, se_corner, sw_corner, &
&   ne_corner, nw_corner, rsin_u, rsin_v, cosa_s, rsin2)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: dord4
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, npx, npy&
&   , grid_type
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: va
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(OUT) :: vc
    LOGICAL, INTENT(IN) :: nested, se_corner, sw_corner, ne_corner, &
&   nw_corner
    REAL, INTENT(IN) :: rsin_u(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rsin_v(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: cosa_s(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rsin2(isd:ied, jsd:jed)
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, PARAMETER :: t11=27./28., t12=-(13./28.), t13=3./7., t14=6./7.&
&   , t15=3./28.
    REAL, PARAMETER :: a1=0.5625
    REAL, PARAMETER :: a2=-0.0625
    REAL, PARAMETER :: c1=-(2./14.)
    REAL, PARAMETER :: c2=11./14.
    REAL, PARAMETER :: c3=5./14.
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    IF (dord4) THEN
      id = 1
    ELSE
      id = 0
    END IF
    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
    IF (nested) THEN
      DO j=jsd+1,jed-1
        DO i=isd,ied
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      DO i=isd,ied
!j = jsd
        utmp(i, jsd) = 0.5*(u(i, jsd)+u(i, jsd+1))
!j = jed
        utmp(i, jed) = 0.5*(u(i, jed)+u(i, jed+1))
      END DO
      DO j=jsd,jed
        DO i=isd+1,ied-1
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
!i = isd
        vtmp(isd, j) = 0.5*(v(isd, j)+v(isd+1, j))
!i = ied
        vtmp(ied, j) = 0.5*(v(ied, j)+v(ied+1, j))
      END DO
      DO j=jsd,jed
        DO i=isd,ied
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        END DO
      END DO
    ELSE
!----------
! Interior:
!----------
      utmp = 0.
      vtmp = 0.
      IF (npt .LT. js - 1) THEN
        max1 = js - 1
      ELSE
        max1 = npt
      END IF
      IF (npy - npt .GT. je + 1) THEN
        min1 = je + 1
      ELSE
        min1 = npy - npt
      END IF
      DO j=max1,min1
        IF (npt .LT. isd) THEN
          max2 = isd
        ELSE
          max2 = npt
        END IF
        IF (npx - npt .GT. ied) THEN
          min2 = ied
        ELSE
          min2 = npx - npt
        END IF
        DO i=max2,min2
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      IF (npt .LT. jsd) THEN
        max3 = jsd
      ELSE
        max3 = npt
      END IF
      IF (npy - npt .GT. jed) THEN
        min3 = jed
      ELSE
        min3 = npy - npt
      END IF
      DO j=max3,min3
        IF (npt .LT. is - 1) THEN
          max4 = is - 1
        ELSE
          max4 = npt
        END IF
        IF (npx - npt .GT. ie + 1) THEN
          min4 = ie + 1
        ELSE
          min4 = npx - npt
        END IF
        DO i=max4,min4
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
      END DO
!----------
! edges:
!----------
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
          DO j=jsd,npt-1
            DO i=isd,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
          DO j=npy-npt+1,jed
            DO i=isd,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (is .EQ. 1 .OR. isd .LT. npt) THEN
          IF (npt .LT. jsd) THEN
            max5 = jsd
          ELSE
            max5 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min5 = jed
          ELSE
            min5 = npy - npt
          END IF
          DO j=max5,min5
            DO i=isd,npt-1
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
          IF (npt .LT. jsd) THEN
            max6 = jsd
          ELSE
            max6 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min6 = jed
          ELSE
            min6 = npy - npt
          END IF
          DO j=max6,min6
            DO i=npx-npt+1,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
      END IF
      DO j=js-1-id,je+1+id
        DO i=is-1-id,ie+1+id
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        END DO
      END DO
    END IF
! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
    IF (sw_corner) THEN
      DO i=-2,0
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
    END IF
    IF (se_corner) THEN
      DO i=0,2
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
    END IF
    IF (ne_corner) THEN
      DO i=0,2
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
    END IF
    IF (nw_corner) THEN
      DO i=-2,0
        utmp(i, npy) = vtmp(0, je+i)
      END DO
    END IF
    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      IF (3 .LT. is - 1) THEN
        ifirst = is - 1
      ELSE
        ifirst = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        ilast = ie + 2
      ELSE
        ilast = npx - 2
      END IF
    ELSE
      ifirst = is - 1
      ilast = ie + 2
    END IF
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
    DO j=js-1,je+1
      DO i=ifirst,ilast
        uc(i, j) = a1*(utmp(i-1, j)+utmp(i, j)) + a2*(utmp(i-2, j)+utmp(&
&         i+1, j))
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
! Xdir:
      IF (is .EQ. 1 .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc(0, j) = c1*utmp(-2, j) + c2*utmp(-1, j) + c3*utmp(0, j)
          uc(1, j) = (t14*(utmp(0, j)+utmp(1, j))+t12*(utmp(-1, j)+utmp(&
&           2, j))+t15*(utmp(-2, j)+utmp(3, j)))*rsin_u(1, j)
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
          uc(npx, j) = (t14*(utmp(npx-1, j)+utmp(npx, j))+t12*(utmp(npx-&
&           2, j)+utmp(npx+1, j))+t15*(utmp(npx-3, j)+utmp(npx+2, j)))*&
&           rsin_u(npx, j)
          uc(npx+1, j) = c3*utmp(npx, j) + c2*utmp(npx+1, j) + c1*utmp(&
&           npx+2, j)
        END DO
      END IF
    END IF
!------
! Ydir:
!------
    IF (sw_corner) THEN
      DO j=-2,0
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
    END IF
    IF (nw_corner) THEN
      DO j=0,2
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
    END IF
    IF (se_corner) THEN
      DO j=-2,0
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
    END IF
    IF (ne_corner) THEN
      DO j=0,2
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1 .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vc(i, 1) = (t14*(vtmp(i, 0)+vtmp(i, 1))+t12*(vtmp(i, -1)+&
&             vtmp(i, 2))+t15*(vtmp(i, -2)+vtmp(i, 3)))*rsin_v(i, 1)
          END DO
        ELSE IF ((j .EQ. 0 .OR. j .EQ. npy - 1) .AND. (.NOT.nested)) &
&       THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
          END DO
        ELSE IF ((j .EQ. 2 .OR. j .EQ. npy + 1) .AND. (.NOT.nested)) &
&       THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j+1) + c2*vtmp(i, j) + c3*vtmp(i, j-1)
          END DO
        ELSE IF (j .EQ. npy .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vc(i, npy) = (t14*(vtmp(i, npy-1)+vtmp(i, npy))+t12*(vtmp(i&
&             , npy-2)+vtmp(i, npy+1))+t15*(vtmp(i, npy-3)+vtmp(i, npy+2&
&             )))*rsin_v(i, npy)
          END DO
        ELSE
! 4th order interpolation for interior points:
          DO i=is-1,ie+1
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
          END DO
        END IF
      END DO
    ELSE
! 4th order interpolation:
      DO j=js-1,je+2
        DO i=is-1,ie+1
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
        END DO
      END DO
    END IF
  END SUBROUTINE D2C_SETUP
  SUBROUTINE D2A_SETUP(u, v, ua, va, dord4, isd, ied, jsd, jed, is, ie, &
&   js, je, npx, npy, grid_type, nested, cosa_s, rsin2)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: dord4
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, npx, npy&
&   , grid_type
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: va
    REAL, INTENT(IN) :: cosa_s(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rsin2(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, PARAMETER :: t11=27./28., t12=-(13./28.), t13=3./7., t14=6./7.&
&   , t15=3./28.
    REAL, PARAMETER :: a1=0.5625
    REAL, PARAMETER :: a2=-0.0625
    REAL, PARAMETER :: c1=-(2./14.)
    REAL, PARAMETER :: c2=11./14.
    REAL, PARAMETER :: c3=5./14.
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    IF (dord4) THEN
      id = 1
    ELSE
      id = 0
    END IF
    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
    IF (nested) THEN
      DO j=jsd+1,jed-1
        DO i=isd,ied
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      DO i=isd,ied
!j = jsd
        utmp(i, jsd) = 0.5*(u(i, jsd)+u(i, jsd+1))
!j = jed
        utmp(i, jed) = 0.5*(u(i, jed)+u(i, jed+1))
      END DO
      DO j=jsd,jed
        DO i=isd+1,ied-1
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
!i = isd
        vtmp(isd, j) = 0.5*(v(isd, j)+v(isd+1, j))
!i = ied
        vtmp(ied, j) = 0.5*(v(ied, j)+v(ied+1, j))
      END DO
    ELSE
      IF (npt .LT. js - 1) THEN
        max1 = js - 1
      ELSE
        max1 = npt
      END IF
      IF (npy - npt .GT. je + 1) THEN
        min1 = je + 1
      ELSE
        min1 = npy - npt
      END IF
!----------
! Interior:
!----------
      DO j=max1,min1
        IF (npt .LT. isd) THEN
          max2 = isd
        ELSE
          max2 = npt
        END IF
        IF (npx - npt .GT. ied) THEN
          min2 = ied
        ELSE
          min2 = npx - npt
        END IF
        DO i=max2,min2
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      IF (npt .LT. jsd) THEN
        max3 = jsd
      ELSE
        max3 = npt
      END IF
      IF (npy - npt .GT. jed) THEN
        min3 = jed
      ELSE
        min3 = npy - npt
      END IF
      DO j=max3,min3
        IF (npt .LT. is - 1) THEN
          max4 = is - 1
        ELSE
          max4 = npt
        END IF
        IF (npx - npt .GT. ie + 1) THEN
          min4 = ie + 1
        ELSE
          min4 = npx - npt
        END IF
        DO i=max4,min4
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
      END DO
!----------
! edges:
!----------
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
          DO j=jsd,npt-1
            DO i=isd,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
          DO j=npy-npt+1,jed
            DO i=isd,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (is .EQ. 1 .OR. isd .LT. npt) THEN
          IF (npt .LT. jsd) THEN
            max5 = jsd
          ELSE
            max5 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min5 = jed
          ELSE
            min5 = npy - npt
          END IF
          DO j=max5,min5
            DO i=isd,npt-1
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
          IF (npt .LT. jsd) THEN
            max6 = jsd
          ELSE
            max6 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min6 = jed
          ELSE
            min6 = npy - npt
          END IF
          DO j=max6,min6
            DO i=npx-npt+1,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
      END IF
    END IF
    DO j=js-1-id,je+1+id
      DO i=is-1-id,ie+1+id
        ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
      END DO
    END DO
  END SUBROUTINE D2A_SETUP
end module fv_restart_adm_mod
