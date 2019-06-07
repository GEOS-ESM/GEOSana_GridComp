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
module a2b_edge_adm_mod

  use fv_grid_utils_mod, only: great_circle_dist
#ifdef VAN2
  use fv_grid_utils_mod, only: van2
#endif

  use fv_arrays_mod,     only: fv_grid_type, R_GRID

  use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                           pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

  implicit none

  real, parameter:: r3 = 1./3.
!----------------------------
! 4-pt Lagrange interpolation
!----------------------------
  real, parameter:: a1 =  0.5625  !  9/16
  real, parameter:: a2 = -0.0625  ! -1/16
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: b1 =  7./12.     ! 0.58333333
  real, parameter:: b2 = -1./12.

  private
  public :: a2b_ord2, a2b_ord4, &
            a2b_ord2_fwd, a2b_ord4_fwd, &
            a2b_ord2_bwd, a2b_ord4_bwd, a2b_ord4_adm

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of a2b_ord4 in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_c
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
!   gradient     of useful results: qin qout
!   with respect to varying inputs: qin
  SUBROUTINE A2B_ORD4_ADM(qin, qin_ad, qout, qout_ad, gridstruct, npx, &
&   npy, is, ie, js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qin_ad(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(INOUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qout_ad(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local: compact 4-pt cubic
    REAL, PARAMETER :: c1=2./3.
    REAL, PARAMETER :: c2=-(1./6.)
! Parabolic spline
! real, parameter:: c1 =  0.75
! real, parameter:: c2 = -0.25
    REAL :: qx(is:ie+1, js-ng:je+ng)
    REAL :: qx_ad(is:ie+1, js-ng:je+ng)
    REAL :: qy(is-ng:ie+ng, js:je+1)
    REAL :: qy_ad(is-ng:ie+ng, js:je+1)
    REAL :: qxx(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qxx_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL :: g_in, g_ou
    REAL :: p0(2)
    REAL :: q1(is-1:ie+1), q2(js-1:je+1)
    REAL :: q1_ad(is-1:ie+1), q2_ad(js-1:je+1)
    INTEGER :: i, j, is1, js1, is2, js2, ie1, je1
    REAL, DIMENSION(:, :, :), POINTER :: grid, agrid
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
    REAL(kind=r_grid), DIMENSION(:), POINTER :: edge_w, edge_e, edge_s, &
&   edge_n
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: max7
    INTEGER :: max8
    INTEGER :: max9
    INTEGER :: max10
    INTEGER :: max11
    INTEGER :: max12
    INTEGER :: max13
    INTEGER :: max14
    INTEGER :: max15
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    INTEGER :: min7
    INTEGER :: min8
    INTEGER :: min9
    INTEGER :: min10
    INTEGER :: min11
    INTEGER :: min12
    INTEGER :: min13
    INTEGER :: min14
    INTEGER :: min15
    REAL :: result1
    REAL :: result1_ad
    REAL :: result2
    REAL :: result2_ad
    REAL :: result3
    REAL :: result3_ad
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
    REAL :: temp_ad22
    REAL :: temp_ad23
    REAL :: temp_ad24
    REAL :: temp_ad25
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
    INTEGER :: branch
    edge_w => gridstruct%edge_w
    edge_e => gridstruct%edge_e
    edge_s => gridstruct%edge_s
    edge_n => gridstruct%edge_n
    agrid => gridstruct%agrid
    grid => gridstruct%grid
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    IF (gridstruct%grid_type .LT. 3) THEN
      IF (1 .LT. is - 1) THEN
        is1 = is - 1
      ELSE
        is1 = 1
      END IF
      IF (1 .LT. js - 1) THEN
        js1 = js - 1
      ELSE
        js1 = 1
      END IF
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (2 .LT. js) THEN
        js2 = js
      ELSE
        js2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        je1 = je + 1
      ELSE
        je1 = npy - 1
      END IF
! Corners:
! 3-way extrapolation
      IF (gridstruct%nested) THEN
        CALL PUSHCONTROL2B(0)
      ELSE
        IF (gridstruct%sw_corner) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (gridstruct%se_corner) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (gridstruct%ne_corner) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (gridstruct%nw_corner) THEN
          p0(1:2) = grid(1, npy, 1:2)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        IF (1 .LT. js - 2) THEN
          max1 = js - 2
        ELSE
          max1 = 1
        END IF
        IF (npy - 1 .GT. je + 2) THEN
          min1 = je + 2
        ELSE
          min1 = npy - 1
        END IF
!------------
! X-Interior:
!------------
        DO j=max1,min1
          IF (3 .LT. is) THEN
            max2 = is
          ELSE
            max2 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min2 = ie + 1
          ELSE
            min2 = npx - 2
          END IF
          ad_from = max2
          i = min2 + 1
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from)
        END DO
! *** West Edges:
        IF (is .EQ. 1) THEN
          IF (1 .LT. js - 2) THEN
            max3 = js - 2
          ELSE
            max3 = 1
          END IF
          IF (npy - 1 .GT. je + 2) THEN
            min3 = je + 2
          ELSE
            min3 = npy - 1
          END IF
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
! East Edges:
        IF (ie + 1 .EQ. npx) THEN
          IF (1 .LT. js - 2) THEN
            max4 = js - 2
          ELSE
            max4 = 1
          END IF
          IF (npy - 1 .GT. je + 2) THEN
            min4 = je + 2
          ELSE
            min4 = npy - 1
          END IF
          CALL PUSHCONTROL2B(1)
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
      END IF
!------------
! Y-Interior:
!------------
      IF (gridstruct%nested) THEN
        CALL PUSHCONTROL2B(0)
      ELSE
        IF (3 .LT. js) THEN
          max5 = js
        ELSE
          max5 = 3
        END IF
        IF (npy - 2 .GT. je + 1) THEN
          min5 = je + 1
        ELSE
          min5 = npy - 2
        END IF
        DO j=max5,min5
          IF (1 .LT. is - 2) THEN
            max6 = is - 2
          ELSE
            max6 = 1
          END IF
          IF (npx - 1 .GT. ie + 2) THEN
            min6 = ie + 2
          ELSE
            min6 = npx - 1
          END IF
          ad_from0 = max6
          i = min6 + 1
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from0)
        END DO
! South Edges:
        IF (js .EQ. 1) THEN
          IF (1 .LT. is - 2) THEN
            max7 = is - 2
          ELSE
            max7 = 1
          END IF
          IF (npx - 1 .GT. ie + 2) THEN
            min7 = ie + 2
          ELSE
            min7 = npx - 1
          END IF
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
! North Edges:
        IF (je + 1 .EQ. npy) THEN
          IF (1 .LT. is - 2) THEN
            max8 = is - 2
          ELSE
            max8 = 1
          END IF
          IF (npx - 1 .GT. ie + 2) THEN
            min8 = ie + 2
          ELSE
            min8 = npx - 1
          END IF
          CALL PUSHCONTROL2B(1)
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
      END IF
!--------------------------------------
      IF (gridstruct%nested) THEN
        CALL PUSHCONTROL2B(0)
      ELSE
        IF (3 .LT. js) THEN
          max9 = js
        ELSE
          max9 = 3
        END IF
        IF (npy - 2 .GT. je + 1) THEN
          min9 = je + 1
        ELSE
          min9 = npy - 2
        END IF
        DO j=max9,min9
          IF (2 .LT. is) THEN
            max10 = is
          ELSE
            max10 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min10 = ie + 1
          ELSE
            min10 = npx - 1
          END IF
          ad_from1 = max10
          i = min10 + 1
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from1)
        END DO
        IF (js .EQ. 1) THEN
          IF (2 .LT. is) THEN
            max11 = is
          ELSE
            max11 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min11 = ie + 1
          ELSE
            min11 = npx - 1
          END IF
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (je + 1 .EQ. npy) THEN
          IF (2 .LT. is) THEN
            max12 = is
          ELSE
            max12 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min12 = ie + 1
          ELSE
            min12 = npx - 1
          END IF
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        IF (2 .LT. js) THEN
          max13 = js
        ELSE
          max13 = 2
        END IF
        IF (npy - 1 .GT. je + 1) THEN
          min13 = je + 1
        ELSE
          min13 = npy - 1
        END IF
        DO j=max13,min13
          IF (3 .LT. is) THEN
            max14 = is
          ELSE
            max14 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min14 = ie + 1
          ELSE
            min14 = npx - 2
          END IF
          ad_from2 = max14
          i = min14 + 1
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from2)
          IF (is .EQ. 1) THEN
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
          IF (2 .LT. is) THEN
            max15 = is
          ELSE
            max15 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min15 = ie + 1
          ELSE
            min15 = npx - 1
          END IF
          ad_from3 = max15
          i = min15 + 1
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from3)
        END DO
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            qout_ad(i, j) = qout_ad(i, j) + qin_ad(i, j)
            qin_ad(i, j) = 0.0
          END DO
        END DO
      END IF
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      qy_ad = 0.0
      qxx_ad = 0.0
      qyy_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          qxx_ad(i, j) = qxx_ad(i, j) + 0.5*qout_ad(i, j)
          qyy_ad(i, j) = qyy_ad(i, j) + 0.5*qout_ad(i, j)
          qout_ad(i, j) = 0.0
        END DO
        DO i=ie+1,is,-1
          qy_ad(i-2, j) = qy_ad(i-2, j) + a2*qyy_ad(i, j)
          qy_ad(i+1, j) = qy_ad(i+1, j) + a2*qyy_ad(i, j)
          qy_ad(i-1, j) = qy_ad(i-1, j) + a1*qyy_ad(i, j)
          qy_ad(i, j) = qy_ad(i, j) + a1*qyy_ad(i, j)
          qyy_ad(i, j) = 0.0
        END DO
      END DO
      qx_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          qx_ad(i, j-2) = qx_ad(i, j-2) + a2*qxx_ad(i, j)
          qx_ad(i, j+1) = qx_ad(i, j+1) + a2*qxx_ad(i, j)
          qx_ad(i, j-1) = qx_ad(i, j-1) + a1*qxx_ad(i, j)
          qx_ad(i, j) = qx_ad(i, j) + a1*qxx_ad(i, j)
          qxx_ad(i, j) = 0.0
        END DO
      END DO
    ELSE IF (branch .EQ. 1) THEN
      qy_ad = 0.0
      qxx_ad = 0.0
      qyy_ad = 0.0
      DO j=min13,max13,-1
        CALL POPINTEGER4(ad_from3)
        CALL POPINTEGER4(ad_to3)
        DO i=ad_to3,ad_from3,-1
          qxx_ad(i, j) = qxx_ad(i, j) + 0.5*qout_ad(i, j)
          qyy_ad(i, j) = qyy_ad(i, j) + 0.5*qout_ad(i, j)
          qout_ad(i, j) = 0.0
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          qy_ad(npx-2, j) = qy_ad(npx-2, j) + c1*qyy_ad(npx-1, j)
          qy_ad(npx-1, j) = qy_ad(npx-1, j) + c1*qyy_ad(npx-1, j)
          qout_ad(npx, j) = qout_ad(npx, j) + c2*qyy_ad(npx-1, j)
          qyy_ad(npx-2, j) = qyy_ad(npx-2, j) + c2*qyy_ad(npx-1, j)
          qyy_ad(npx-1, j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          qy_ad(1, j) = qy_ad(1, j) + c1*qyy_ad(2, j)
          qy_ad(2, j) = qy_ad(2, j) + c1*qyy_ad(2, j)
          qout_ad(1, j) = qout_ad(1, j) + c2*qyy_ad(2, j)
          qyy_ad(3, j) = qyy_ad(3, j) + c2*qyy_ad(2, j)
          qyy_ad(2, j) = 0.0
        END IF
        CALL POPINTEGER4(ad_from2)
        CALL POPINTEGER4(ad_to2)
        DO i=ad_to2,ad_from2,-1
          qy_ad(i-2, j) = qy_ad(i-2, j) + a2*qyy_ad(i, j)
          qy_ad(i+1, j) = qy_ad(i+1, j) + a2*qyy_ad(i, j)
          qy_ad(i-1, j) = qy_ad(i-1, j) + a1*qyy_ad(i, j)
          qy_ad(i, j) = qy_ad(i, j) + a1*qyy_ad(i, j)
          qyy_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        qx_ad = 0.0
      ELSE
        qx_ad = 0.0
        DO i=min12,max12,-1
          qx_ad(i, npy-2) = qx_ad(i, npy-2) + c1*qxx_ad(i, npy-1)
          qx_ad(i, npy-1) = qx_ad(i, npy-1) + c1*qxx_ad(i, npy-1)
          qout_ad(i, npy) = qout_ad(i, npy) + c2*qxx_ad(i, npy-1)
          qxx_ad(i, npy-2) = qxx_ad(i, npy-2) + c2*qxx_ad(i, npy-1)
          qxx_ad(i, npy-1) = 0.0
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO i=min11,max11,-1
          qx_ad(i, 1) = qx_ad(i, 1) + c1*qxx_ad(i, 2)
          qx_ad(i, 2) = qx_ad(i, 2) + c1*qxx_ad(i, 2)
          qout_ad(i, 1) = qout_ad(i, 1) + c2*qxx_ad(i, 2)
          qxx_ad(i, 3) = qxx_ad(i, 3) + c2*qxx_ad(i, 2)
          qxx_ad(i, 2) = 0.0
        END DO
      END IF
      DO j=min9,max9,-1
        CALL POPINTEGER4(ad_from1)
        CALL POPINTEGER4(ad_to1)
        DO i=ad_to1,ad_from1,-1
          qx_ad(i, j-2) = qx_ad(i, j-2) + a2*qxx_ad(i, j)
          qx_ad(i, j+1) = qx_ad(i, j+1) + a2*qxx_ad(i, j)
          qx_ad(i, j-1) = qx_ad(i, j-1) + a1*qxx_ad(i, j)
          qx_ad(i, j) = qx_ad(i, j) + a1*qxx_ad(i, j)
          qxx_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      qx_ad = 0.0
      qy_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          temp_ad23 = 0.5*qout_ad(i, j)
          temp_ad24 = a1*temp_ad23
          temp_ad25 = a2*temp_ad23
          qx_ad(i, j-1) = qx_ad(i, j-1) + temp_ad24
          qx_ad(i, j) = qx_ad(i, j) + temp_ad24
          qy_ad(i-1, j) = qy_ad(i-1, j) + temp_ad24
          qy_ad(i, j) = qy_ad(i, j) + temp_ad24
          qx_ad(i, j-2) = qx_ad(i, j-2) + temp_ad25
          qx_ad(i, j+1) = qx_ad(i, j+1) + temp_ad25
          qy_ad(i-2, j) = qy_ad(i-2, j) + temp_ad25
          qy_ad(i+1, j) = qy_ad(i+1, j) + temp_ad25
          qout_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+2,is-2,-1
          qin_ad(i, j-1) = qin_ad(i, j-1) + b1*qy_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qy_ad(i, j)
          qin_ad(i, j-2) = qin_ad(i, j-2) + b2*qy_ad(i, j)
          qin_ad(i, j+1) = qin_ad(i, j+1) + b2*qy_ad(i, j)
          qy_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+2,js-2,-1
        DO i=ie+1,is,-1
          qin_ad(i-1, j) = qin_ad(i-1, j) + b1*qx_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qx_ad(i, j)
          qin_ad(i-2, j) = qin_ad(i-2, j) + b2*qx_ad(i, j)
          qin_ad(i+1, j) = qin_ad(i+1, j) + b2*qx_ad(i, j)
          qx_ad(i, j) = 0.0
        END DO
      END DO
      GOTO 100
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je+1,js,-1
        DO i=ie+2,is-2,-1
          qin_ad(i, j-2) = qin_ad(i, j-2) + b2*qy_ad(i, j)
          qin_ad(i, j+1) = qin_ad(i, j+1) + b2*qy_ad(i, j)
          qin_ad(i, j-1) = qin_ad(i, j-1) + b1*qy_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qy_ad(i, j)
          qy_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      IF (branch .EQ. 1) THEN
        DO i=min8,max8,-1
          g_in = dya(i, npy-2)/dya(i, npy-1)
          temp_ad19 = qy_ad(i, npy-1)/(g_in*2.+2.)
          qin_ad(i, npy-2) = qin_ad(i, npy-2) + 3.*temp_ad19
          qy_ad(i, npy) = qy_ad(i, npy) - g_in*temp_ad19
          qy_ad(i, npy-2) = qy_ad(i, npy-2) - temp_ad19
          qy_ad(i, npy-1) = 0.0
          g_ou = dya(i, npy+1)/dya(i, npy)
          temp_ad21 = 0.5*qy_ad(i, npy)
          temp_ad20 = temp_ad21/(g_in+1.)
          qin_ad(i, npy-1) = qin_ad(i, npy-1) + (g_in+2.)*temp_ad20 + 3.&
&           *g_in*temp_ad19
          temp_ad22 = temp_ad21/(g_ou+1.)
          qin_ad(i, npy-2) = qin_ad(i, npy-2) - temp_ad20
          qin_ad(i, npy) = qin_ad(i, npy) + (g_ou+2.)*temp_ad22
          qin_ad(i, npy+1) = qin_ad(i, npy+1) - temp_ad22
          qy_ad(i, npy) = 0.0
        END DO
        q1_ad = 0.0
        DO i=ie1,is2,-1
          q1_ad(i-1) = q1_ad(i-1) + edge_n(i)*qout_ad(i, npy)
          q1_ad(i) = q1_ad(i) + (1.-edge_n(i))*qout_ad(i, npy)
          qout_ad(i, npy) = 0.0
        END DO
        DO i=ie1,is1,-1
          temp_ad18 = q1_ad(i)/(dya(i, npy-1)+dya(i, npy))
          qin_ad(i, npy-1) = qin_ad(i, npy-1) + dya(i, npy)*temp_ad18
          qin_ad(i, npy) = qin_ad(i, npy) + dya(i, npy-1)*temp_ad18
          q1_ad(i) = 0.0
        END DO
      ELSE
        q1_ad = 0.0
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO i=min7,max7,-1
          g_in = dya(i, 2)/dya(i, 1)
          temp_ad14 = qy_ad(i, 2)/(g_in*2.+2.)
          qin_ad(i, 1) = qin_ad(i, 1) + 3.*g_in*temp_ad14
          qin_ad(i, 2) = qin_ad(i, 2) + 3.*temp_ad14
          qy_ad(i, 1) = qy_ad(i, 1) - g_in*temp_ad14
          qy_ad(i, 3) = qy_ad(i, 3) - temp_ad14
          qy_ad(i, 2) = 0.0
          g_ou = dya(i, -1)/dya(i, 0)
          temp_ad15 = 0.5*qy_ad(i, 1)
          temp_ad16 = temp_ad15/(g_in+1.)
          temp_ad17 = temp_ad15/(g_ou+1.)
          qin_ad(i, 1) = qin_ad(i, 1) + (g_in+2.)*temp_ad16
          qin_ad(i, 2) = qin_ad(i, 2) - temp_ad16
          qin_ad(i, 0) = qin_ad(i, 0) + (g_ou+2.)*temp_ad17
          qin_ad(i, -1) = qin_ad(i, -1) - temp_ad17
          qy_ad(i, 1) = 0.0
        END DO
        DO i=ie1,is2,-1
          q1_ad(i-1) = q1_ad(i-1) + edge_s(i)*qout_ad(i, 1)
          q1_ad(i) = q1_ad(i) + (1.-edge_s(i))*qout_ad(i, 1)
          qout_ad(i, 1) = 0.0
        END DO
        DO i=ie1,is1,-1
          temp_ad13 = q1_ad(i)/(dya(i, 0)+dya(i, 1))
          qin_ad(i, 0) = qin_ad(i, 0) + dya(i, 1)*temp_ad13
          qin_ad(i, 1) = qin_ad(i, 1) + dya(i, 0)*temp_ad13
          q1_ad(i) = 0.0
        END DO
      END IF
      DO j=min5,max5,-1
        CALL POPINTEGER4(ad_from0)
        CALL POPINTEGER4(ad_to0)
        DO i=ad_to0,ad_from0,-1
          qin_ad(i, j-2) = qin_ad(i, j-2) + b2*qy_ad(i, j)
          qin_ad(i, j+1) = qin_ad(i, j+1) + b2*qy_ad(i, j)
          qin_ad(i, j-1) = qin_ad(i, j-1) + b1*qy_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qy_ad(i, j)
          qy_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je+2,js-2,-1
        DO i=ie+1,is,-1
          qin_ad(i-2, j) = qin_ad(i-2, j) + b2*qx_ad(i, j)
          qin_ad(i+1, j) = qin_ad(i+1, j) + b2*qx_ad(i, j)
          qin_ad(i-1, j) = qin_ad(i-1, j) + b1*qx_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qx_ad(i, j)
          qx_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      IF (branch .EQ. 1) THEN
        DO j=min4,max4,-1
          g_in = dxa(npx-2, j)/dxa(npx-1, j)
          temp_ad9 = qx_ad(npx-1, j)/(g_in*2.+2.)
          qin_ad(npx-2, j) = qin_ad(npx-2, j) + 3.*temp_ad9
          qx_ad(npx, j) = qx_ad(npx, j) - g_in*temp_ad9
          qx_ad(npx-2, j) = qx_ad(npx-2, j) - temp_ad9
          qx_ad(npx-1, j) = 0.0
          g_ou = dxa(npx+1, j)/dxa(npx, j)
          temp_ad11 = 0.5*qx_ad(npx, j)
          temp_ad10 = temp_ad11/(g_in+1.)
          qin_ad(npx-1, j) = qin_ad(npx-1, j) + (g_in+2.)*temp_ad10 + 3.&
&           *g_in*temp_ad9
          temp_ad12 = temp_ad11/(g_ou+1.)
          qin_ad(npx-2, j) = qin_ad(npx-2, j) - temp_ad10
          qin_ad(npx, j) = qin_ad(npx, j) + (g_ou+2.)*temp_ad12
          qin_ad(npx+1, j) = qin_ad(npx+1, j) - temp_ad12
          qx_ad(npx, j) = 0.0
        END DO
        q2_ad = 0.0
        DO j=je1,js2,-1
          q2_ad(j-1) = q2_ad(j-1) + edge_e(j)*qout_ad(npx, j)
          q2_ad(j) = q2_ad(j) + (1.-edge_e(j))*qout_ad(npx, j)
          qout_ad(npx, j) = 0.0
        END DO
        DO j=je1,js1,-1
          temp_ad8 = q2_ad(j)/(dxa(npx-1, j)+dxa(npx, j))
          qin_ad(npx-1, j) = qin_ad(npx-1, j) + dxa(npx, j)*temp_ad8
          qin_ad(npx, j) = qin_ad(npx, j) + dxa(npx-1, j)*temp_ad8
          q2_ad(j) = 0.0
        END DO
      ELSE
        q2_ad = 0.0
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=min3,max3,-1
          g_in = dxa(2, j)/dxa(1, j)
          temp_ad4 = qx_ad(2, j)/(g_in*2.+2.)
          qin_ad(1, j) = qin_ad(1, j) + 3.*g_in*temp_ad4
          qin_ad(2, j) = qin_ad(2, j) + 3.*temp_ad4
          qx_ad(1, j) = qx_ad(1, j) - g_in*temp_ad4
          qx_ad(3, j) = qx_ad(3, j) - temp_ad4
          qx_ad(2, j) = 0.0
          g_ou = dxa(-1, j)/dxa(0, j)
          temp_ad5 = 0.5*qx_ad(1, j)
          temp_ad6 = temp_ad5/(g_in+1.)
          temp_ad7 = temp_ad5/(g_ou+1.)
          qin_ad(1, j) = qin_ad(1, j) + (g_in+2.)*temp_ad6
          qin_ad(2, j) = qin_ad(2, j) - temp_ad6
          qin_ad(0, j) = qin_ad(0, j) + (g_ou+2.)*temp_ad7
          qin_ad(-1, j) = qin_ad(-1, j) - temp_ad7
          qx_ad(1, j) = 0.0
        END DO
        DO j=je1,js2,-1
          q2_ad(j-1) = q2_ad(j-1) + edge_w(j)*qout_ad(1, j)
          q2_ad(j) = q2_ad(j) + (1.-edge_w(j))*qout_ad(1, j)
          qout_ad(1, j) = 0.0
        END DO
        DO j=je1,js1,-1
          temp_ad3 = q2_ad(j)/(dxa(0, j)+dxa(1, j))
          qin_ad(0, j) = qin_ad(0, j) + dxa(1, j)*temp_ad3
          qin_ad(1, j) = qin_ad(1, j) + dxa(0, j)*temp_ad3
          q2_ad(j) = 0.0
        END DO
      END IF
      DO j=min1,max1,-1
        CALL POPINTEGER4(ad_from)
        CALL POPINTEGER4(ad_to)
        DO i=ad_to,ad_from,-1
          qin_ad(i-2, j) = qin_ad(i-2, j) + b2*qx_ad(i, j)
          qin_ad(i+1, j) = qin_ad(i+1, j) + b2*qx_ad(i, j)
          qin_ad(i-1, j) = qin_ad(i-1, j) + b1*qx_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qx_ad(i, j)
          qx_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        temp_ad2 = r3*qout_ad(1, npy)
        result1_ad = temp_ad2
        result2_ad = temp_ad2
        result3_ad = temp_ad2
        qout_ad(1, npy) = 0.0
        CALL EXTRAP_CORNER_ADM(p0, agrid(1, npy, 1:2), agrid(2, npy+1, 1&
&                        :2), qin(1, npy), qin_ad(1, npy), qin(2, npy+1)&
&                        , qin_ad(2, npy+1), result3_ad)
        CALL EXTRAP_CORNER_ADM(p0, agrid(0, npy-1, 1:2), agrid(-1, npy-2&
&                        , 1:2), qin(0, npy-1), qin_ad(0, npy-1), qin(-1&
&                        , npy-2), qin_ad(-1, npy-2), result2_ad)
        CALL EXTRAP_CORNER_ADM(p0, agrid(1, npy-1, 1:2), agrid(2, npy-2&
&                        , 1:2), qin(1, npy-1), qin_ad(1, npy-1), qin(2&
&                        , npy-2), qin_ad(2, npy-2), result1_ad)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp_ad1 = r3*qout_ad(npx, npy)
        result1_ad = temp_ad1
        result2_ad = temp_ad1
        result3_ad = temp_ad1
        qout_ad(npx, npy) = 0.0
        p0(1:2) = grid(npx, npy, 1:2)
        CALL EXTRAP_CORNER_ADM(p0, agrid(npx-1, npy, 1:2), agrid(npx-2, &
&                        npy+1, 1:2), qin(npx-1, npy), qin_ad(npx-1, npy&
&                        ), qin(npx-2, npy+1), qin_ad(npx-2, npy+1), &
&                        result3_ad)
        CALL EXTRAP_CORNER_ADM(p0, agrid(npx, npy-1, 1:2), agrid(npx+1, &
&                        npy-2, 1:2), qin(npx, npy-1), qin_ad(npx, npy-1&
&                        ), qin(npx+1, npy-2), qin_ad(npx+1, npy-2), &
&                        result2_ad)
        CALL EXTRAP_CORNER_ADM(p0, agrid(npx-1, npy-1, 1:2), agrid(npx-2&
&                        , npy-2, 1:2), qin(npx-1, npy-1), qin_ad(npx-1&
&                        , npy-1), qin(npx-2, npy-2), qin_ad(npx-2, npy-&
&                        2), result1_ad)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp_ad0 = r3*qout_ad(npx, 1)
        result1_ad = temp_ad0
        result2_ad = temp_ad0
        result3_ad = temp_ad0
        qout_ad(npx, 1) = 0.0
        p0(1:2) = grid(npx, 1, 1:2)
        CALL EXTRAP_CORNER_ADM(p0, agrid(npx, 1, 1:2), agrid(npx+1, 2, 1&
&                        :2), qin(npx, 1), qin_ad(npx, 1), qin(npx+1, 2)&
&                        , qin_ad(npx+1, 2), result3_ad)
        CALL EXTRAP_CORNER_ADM(p0, agrid(npx-1, 0, 1:2), agrid(npx-2, -1&
&                        , 1:2), qin(npx-1, 0), qin_ad(npx-1, 0), qin(&
&                        npx-2, -1), qin_ad(npx-2, -1), result2_ad)
        CALL EXTRAP_CORNER_ADM(p0, agrid(npx-1, 1, 1:2), agrid(npx-2, 2&
&                        , 1:2), qin(npx-1, 1), qin_ad(npx-1, 1), qin(&
&                        npx-2, 2), qin_ad(npx-2, 2), result1_ad)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp_ad = r3*qout_ad(1, 1)
        result1_ad = temp_ad
        result2_ad = temp_ad
        result3_ad = temp_ad
        p0(1:2) = grid(1, 1, 1:2)
        CALL EXTRAP_CORNER_ADM(p0, agrid(1, 0, 1:2), agrid(2, -1, 1:2), &
&                        qin(1, 0), qin_ad(1, 0), qin(2, -1), qin_ad(2, &
&                        -1), result3_ad)
        CALL EXTRAP_CORNER_ADM(p0, agrid(0, 1, 1:2), agrid(-1, 2, 1:2), &
&                        qin(0, 1), qin_ad(0, 1), qin(-1, 2), qin_ad(-1&
&                        , 2), result2_ad)
        CALL EXTRAP_CORNER_ADM(p0, agrid(1, 1, 1:2), agrid(2, 2, 1:2), &
&                        qin(1, 1), qin_ad(1, 1), qin(2, 2), qin_ad(2, 2&
&                        ), result1_ad)
      END IF
    END IF
 100 CONTINUE
  END SUBROUTINE A2B_ORD4_ADM
  SUBROUTINE A2B_ORD4(qin, qout, gridstruct, npx, npy, is, ie, js, je, &
&   ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(INOUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local: compact 4-pt cubic
    REAL, PARAMETER :: c1=2./3.
    REAL, PARAMETER :: c2=-(1./6.)
! Parabolic spline
! real, parameter:: c1 =  0.75
! real, parameter:: c2 = -0.25
    REAL :: qx(is:ie+1, js-ng:je+ng)
    REAL :: qy(is-ng:ie+ng, js:je+1)
    REAL :: qxx(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy(is-ng:ie+ng, js-ng:je+ng)
    REAL :: g_in, g_ou
    REAL :: p0(2)
    REAL :: q1(is-1:ie+1), q2(js-1:je+1)
    INTEGER :: i, j, is1, js1, is2, js2, ie1, je1
    REAL, DIMENSION(:, :, :), POINTER :: grid, agrid
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
    REAL(kind=r_grid), DIMENSION(:), POINTER :: edge_w, edge_e, edge_s, &
&   edge_n
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: max7
    INTEGER :: max8
    INTEGER :: max9
    INTEGER :: max10
    INTEGER :: max11
    INTEGER :: max12
    INTEGER :: max13
    INTEGER :: max14
    INTEGER :: max15
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    INTEGER :: min7
    INTEGER :: min8
    INTEGER :: min9
    INTEGER :: min10
    INTEGER :: min11
    INTEGER :: min12
    INTEGER :: min13
    INTEGER :: min14
    INTEGER :: min15
    REAL :: result1
    REAL :: result2
    REAL :: result3
    edge_w => gridstruct%edge_w
    edge_e => gridstruct%edge_e
    edge_s => gridstruct%edge_s
    edge_n => gridstruct%edge_n
    grid => gridstruct%grid
    agrid => gridstruct%agrid
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    IF (gridstruct%grid_type .LT. 3) THEN
      IF (1 .LT. is - 1) THEN
        is1 = is - 1
      ELSE
        is1 = 1
      END IF
      IF (1 .LT. js - 1) THEN
        js1 = js - 1
      ELSE
        js1 = 1
      END IF
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (2 .LT. js) THEN
        js2 = js
      ELSE
        js2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        je1 = je + 1
      ELSE
        je1 = npy - 1
      END IF
! Corners:
! 3-way extrapolation
      IF (gridstruct%nested) THEN
        DO j=js-2,je+2
          DO i=is,ie+1
            qx(i, j) = b2*(qin(i-2, j)+qin(i+1, j)) + b1*(qin(i-1, j)+&
&             qin(i, j))
          END DO
        END DO
      ELSE
        IF (gridstruct%sw_corner) THEN
          p0(1:2) = grid(1, 1, 1:2)
          result1 = EXTRAP_CORNER(p0, agrid(1, 1, 1:2), agrid(2, 2, 1:2)&
&           , qin(1, 1), qin(2, 2))
          result2 = EXTRAP_CORNER(p0, agrid(0, 1, 1:2), agrid(-1, 2, 1:2&
&           ), qin(0, 1), qin(-1, 2))
          result3 = EXTRAP_CORNER(p0, agrid(1, 0, 1:2), agrid(2, -1, 1:2&
&           ), qin(1, 0), qin(2, -1))
          qout(1, 1) = (result1+result2+result3)*r3
        END IF
        IF (gridstruct%se_corner) THEN
          p0(1:2) = grid(npx, 1, 1:2)
          result1 = EXTRAP_CORNER(p0, agrid(npx-1, 1, 1:2), agrid(npx-2&
&           , 2, 1:2), qin(npx-1, 1), qin(npx-2, 2))
          result2 = EXTRAP_CORNER(p0, agrid(npx-1, 0, 1:2), agrid(npx-2&
&           , -1, 1:2), qin(npx-1, 0), qin(npx-2, -1))
          result3 = EXTRAP_CORNER(p0, agrid(npx, 1, 1:2), agrid(npx+1, 2&
&           , 1:2), qin(npx, 1), qin(npx+1, 2))
          qout(npx, 1) = (result1+result2+result3)*r3
        END IF
        IF (gridstruct%ne_corner) THEN
          p0(1:2) = grid(npx, npy, 1:2)
          result1 = EXTRAP_CORNER(p0, agrid(npx-1, npy-1, 1:2), agrid(&
&           npx-2, npy-2, 1:2), qin(npx-1, npy-1), qin(npx-2, npy-2))
          result2 = EXTRAP_CORNER(p0, agrid(npx, npy-1, 1:2), agrid(npx+&
&           1, npy-2, 1:2), qin(npx, npy-1), qin(npx+1, npy-2))
          result3 = EXTRAP_CORNER(p0, agrid(npx-1, npy, 1:2), agrid(npx-&
&           2, npy+1, 1:2), qin(npx-1, npy), qin(npx-2, npy+1))
          qout(npx, npy) = (result1+result2+result3)*r3
        END IF
        IF (gridstruct%nw_corner) THEN
          p0(1:2) = grid(1, npy, 1:2)
          result1 = EXTRAP_CORNER(p0, agrid(1, npy-1, 1:2), agrid(2, npy&
&           -2, 1:2), qin(1, npy-1), qin(2, npy-2))
          result2 = EXTRAP_CORNER(p0, agrid(0, npy-1, 1:2), agrid(-1, &
&           npy-2, 1:2), qin(0, npy-1), qin(-1, npy-2))
          result3 = EXTRAP_CORNER(p0, agrid(1, npy, 1:2), agrid(2, npy+1&
&           , 1:2), qin(1, npy), qin(2, npy+1))
          qout(1, npy) = (result1+result2+result3)*r3
        END IF
        IF (1 .LT. js - 2) THEN
          max1 = js - 2
        ELSE
          max1 = 1
        END IF
        IF (npy - 1 .GT. je + 2) THEN
          min1 = je + 2
        ELSE
          min1 = npy - 1
        END IF
!------------
! X-Interior:
!------------
        DO j=max1,min1
          IF (3 .LT. is) THEN
            max2 = is
          ELSE
            max2 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min2 = ie + 1
          ELSE
            min2 = npx - 2
          END IF
          DO i=max2,min2
            qx(i, j) = b2*(qin(i-2, j)+qin(i+1, j)) + b1*(qin(i-1, j)+&
&             qin(i, j))
          END DO
        END DO
! *** West Edges:
        IF (is .EQ. 1) THEN
          DO j=js1,je1
            q2(j) = (qin(0, j)*dxa(1, j)+qin(1, j)*dxa(0, j))/(dxa(0, j)&
&             +dxa(1, j))
          END DO
          DO j=js2,je1
            qout(1, j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
          END DO
          IF (1 .LT. js - 2) THEN
            max3 = js - 2
          ELSE
            max3 = 1
          END IF
          IF (npy - 1 .GT. je + 2) THEN
            min3 = je + 2
          ELSE
            min3 = npy - 1
          END IF
!
          DO j=max3,min3
            g_in = dxa(2, j)/dxa(1, j)
            g_ou = dxa(-1, j)/dxa(0, j)
            qx(1, j) = 0.5*(((2.+g_in)*qin(1, j)-qin(2, j))/(1.+g_in)+((&
&             2.+g_ou)*qin(0, j)-qin(-1, j))/(1.+g_ou))
            qx(2, j) = (3.*(g_in*qin(1, j)+qin(2, j))-(g_in*qx(1, j)+qx(&
&             3, j)))/(2.+2.*g_in)
          END DO
        END IF
! East Edges:
        IF (ie + 1 .EQ. npx) THEN
          DO j=js1,je1
            q2(j) = (qin(npx-1, j)*dxa(npx, j)+qin(npx, j)*dxa(npx-1, j)&
&             )/(dxa(npx-1, j)+dxa(npx, j))
          END DO
          DO j=js2,je1
            qout(npx, j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
          END DO
          IF (1 .LT. js - 2) THEN
            max4 = js - 2
          ELSE
            max4 = 1
          END IF
          IF (npy - 1 .GT. je + 2) THEN
            min4 = je + 2
          ELSE
            min4 = npy - 1
          END IF
!
          DO j=max4,min4
            g_in = dxa(npx-2, j)/dxa(npx-1, j)
            g_ou = dxa(npx+1, j)/dxa(npx, j)
            qx(npx, j) = 0.5*(((2.+g_in)*qin(npx-1, j)-qin(npx-2, j))/(&
&             1.+g_in)+((2.+g_ou)*qin(npx, j)-qin(npx+1, j))/(1.+g_ou))
            qx(npx-1, j) = (3.*(qin(npx-2, j)+g_in*qin(npx-1, j))-(g_in*&
&             qx(npx, j)+qx(npx-2, j)))/(2.+2.*g_in)
          END DO
        END IF
      END IF
!------------
! Y-Interior:
!------------
      IF (gridstruct%nested) THEN
        DO j=js,je+1
          DO i=is-2,ie+2
            qy(i, j) = b2*(qin(i, j-2)+qin(i, j+1)) + b1*(qin(i, j-1)+&
&             qin(i, j))
          END DO
        END DO
      ELSE
        IF (3 .LT. js) THEN
          max5 = js
        ELSE
          max5 = 3
        END IF
        IF (npy - 2 .GT. je + 1) THEN
          min5 = je + 1
        ELSE
          min5 = npy - 2
        END IF
        DO j=max5,min5
          IF (1 .LT. is - 2) THEN
            max6 = is - 2
          ELSE
            max6 = 1
          END IF
          IF (npx - 1 .GT. ie + 2) THEN
            min6 = ie + 2
          ELSE
            min6 = npx - 1
          END IF
          DO i=max6,min6
            qy(i, j) = b2*(qin(i, j-2)+qin(i, j+1)) + b1*(qin(i, j-1)+&
&             qin(i, j))
          END DO
        END DO
! South Edges:
        IF (js .EQ. 1) THEN
          DO i=is1,ie1
            q1(i) = (qin(i, 0)*dya(i, 1)+qin(i, 1)*dya(i, 0))/(dya(i, 0)&
&             +dya(i, 1))
          END DO
          DO i=is2,ie1
            qout(i, 1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
          END DO
          IF (1 .LT. is - 2) THEN
            max7 = is - 2
          ELSE
            max7 = 1
          END IF
          IF (npx - 1 .GT. ie + 2) THEN
            min7 = ie + 2
          ELSE
            min7 = npx - 1
          END IF
!
          DO i=max7,min7
            g_in = dya(i, 2)/dya(i, 1)
            g_ou = dya(i, -1)/dya(i, 0)
            qy(i, 1) = 0.5*(((2.+g_in)*qin(i, 1)-qin(i, 2))/(1.+g_in)+((&
&             2.+g_ou)*qin(i, 0)-qin(i, -1))/(1.+g_ou))
            qy(i, 2) = (3.*(g_in*qin(i, 1)+qin(i, 2))-(g_in*qy(i, 1)+qy(&
&             i, 3)))/(2.+2.*g_in)
          END DO
        END IF
! North Edges:
        IF (je + 1 .EQ. npy) THEN
          DO i=is1,ie1
            q1(i) = (qin(i, npy-1)*dya(i, npy)+qin(i, npy)*dya(i, npy-1)&
&             )/(dya(i, npy-1)+dya(i, npy))
          END DO
          DO i=is2,ie1
            qout(i, npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
          END DO
          IF (1 .LT. is - 2) THEN
            max8 = is - 2
          ELSE
            max8 = 1
          END IF
          IF (npx - 1 .GT. ie + 2) THEN
            min8 = ie + 2
          ELSE
            min8 = npx - 1
          END IF
!
          DO i=max8,min8
            g_in = dya(i, npy-2)/dya(i, npy-1)
            g_ou = dya(i, npy+1)/dya(i, npy)
            qy(i, npy) = 0.5*(((2.+g_in)*qin(i, npy-1)-qin(i, npy-2))/(&
&             1.+g_in)+((2.+g_ou)*qin(i, npy)-qin(i, npy+1))/(1.+g_ou))
            qy(i, npy-1) = (3.*(qin(i, npy-2)+g_in*qin(i, npy-1))-(g_in*&
&             qy(i, npy)+qy(i, npy-2)))/(2.+2.*g_in)
          END DO
        END IF
      END IF
!--------------------------------------
      IF (gridstruct%nested) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qxx(i, j) = a2*(qx(i, j-2)+qx(i, j+1)) + a1*(qx(i, j-1)+qx(i&
&             , j))
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            qyy(i, j) = a2*(qy(i-2, j)+qy(i+1, j)) + a1*(qy(i-1, j)+qy(i&
&             , j))
          END DO
          DO i=is,ie+1
! averaging
            qout(i, j) = 0.5*(qxx(i, j)+qyy(i, j))
          END DO
        END DO
      ELSE
        IF (3 .LT. js) THEN
          max9 = js
        ELSE
          max9 = 3
        END IF
        IF (npy - 2 .GT. je + 1) THEN
          min9 = je + 1
        ELSE
          min9 = npy - 2
        END IF
        DO j=max9,min9
          IF (2 .LT. is) THEN
            max10 = is
          ELSE
            max10 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min10 = ie + 1
          ELSE
            min10 = npx - 1
          END IF
          DO i=max10,min10
            qxx(i, j) = a2*(qx(i, j-2)+qx(i, j+1)) + a1*(qx(i, j-1)+qx(i&
&             , j))
          END DO
        END DO
        IF (js .EQ. 1) THEN
          IF (2 .LT. is) THEN
            max11 = is
          ELSE
            max11 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min11 = ie + 1
          ELSE
            min11 = npx - 1
          END IF
          DO i=max11,min11
            qxx(i, 2) = c1*(qx(i, 1)+qx(i, 2)) + c2*(qout(i, 1)+qxx(i, 3&
&             ))
          END DO
        END IF
        IF (je + 1 .EQ. npy) THEN
          IF (2 .LT. is) THEN
            max12 = is
          ELSE
            max12 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min12 = ie + 1
          ELSE
            min12 = npx - 1
          END IF
          DO i=max12,min12
            qxx(i, npy-1) = c1*(qx(i, npy-2)+qx(i, npy-1)) + c2*(qout(i&
&             , npy)+qxx(i, npy-2))
          END DO
        END IF
        IF (2 .LT. js) THEN
          max13 = js
        ELSE
          max13 = 2
        END IF
        IF (npy - 1 .GT. je + 1) THEN
          min13 = je + 1
        ELSE
          min13 = npy - 1
        END IF
        DO j=max13,min13
          IF (3 .LT. is) THEN
            max14 = is
          ELSE
            max14 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min14 = ie + 1
          ELSE
            min14 = npx - 2
          END IF
          DO i=max14,min14
            qyy(i, j) = a2*(qy(i-2, j)+qy(i+1, j)) + a1*(qy(i-1, j)+qy(i&
&             , j))
          END DO
          IF (is .EQ. 1) qyy(2, j) = c1*(qy(1, j)+qy(2, j)) + c2*(qout(1&
&             , j)+qyy(3, j))
          IF (ie + 1 .EQ. npx) qyy(npx-1, j) = c1*(qy(npx-2, j)+qy(npx-1&
&             , j)) + c2*(qout(npx, j)+qyy(npx-2, j))
          IF (2 .LT. is) THEN
            max15 = is
          ELSE
            max15 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min15 = ie + 1
          ELSE
            min15 = npx - 1
          END IF
          DO i=max15,min15
! averaging
            qout(i, j) = 0.5*(qxx(i, j)+qyy(i, j))
          END DO
        END DO
      END IF
    ELSE
! grid_type>=3
!------------------------
! Doubly periodic domain:
!------------------------
! X-sweep: PPM
      DO j=js-2,je+2
        DO i=is,ie+1
          qx(i, j) = b1*(qin(i-1, j)+qin(i, j)) + b2*(qin(i-2, j)+qin(i+&
&           1, j))
        END DO
      END DO
! Y-sweep: PPM
      DO j=js,je+1
        DO i=is-2,ie+2
          qy(i, j) = b1*(qin(i, j-1)+qin(i, j)) + b2*(qin(i, j-2)+qin(i&
&           , j+1))
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie+1
          qout(i, j) = 0.5*(a1*(qx(i, j-1)+qx(i, j)+qy(i-1, j)+qy(i, j))&
&           +a2*(qx(i, j-2)+qx(i, j+1)+qy(i-2, j)+qy(i+1, j)))
        END DO
      END DO
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qin(i, j) = qout(i, j)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE A2B_ORD4
!  Differentiation of a2b_ord4 in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: qin qout
!   with respect to varying inputs: qin qout
  SUBROUTINE A2B_ORD4_FWD(qin, qout, gridstruct, npx, npy, is, ie, js&
&   , je, ng, replace)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(INOUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local: compact 4-pt cubic
    REAL, PARAMETER :: c1=2./3.
    REAL, PARAMETER :: c2=-(1./6.)
! Parabolic spline
! real, parameter:: c1 =  0.75
! real, parameter:: c2 = -0.25
    REAL :: qx(is:ie+1, js-ng:je+ng)
    REAL :: qy(is-ng:ie+ng, js:je+1)
    REAL :: qxx(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy(is-ng:ie+ng, js-ng:je+ng)
    REAL :: g_in, g_ou
    REAL :: p0(2)
    REAL :: q1(is-1:ie+1), q2(js-1:je+1)
    INTEGER :: i, j, is1, js1, is2, js2, ie1, je1
    REAL, DIMENSION(:, :, :), POINTER :: grid, agrid
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
    REAL(kind=r_grid), DIMENSION(:), POINTER :: edge_w, edge_e, edge_s, &
&   edge_n
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: max7
    INTEGER :: max8
    INTEGER :: max9
    INTEGER :: max10
    INTEGER :: max11
    INTEGER :: max12
    INTEGER :: max13
    INTEGER :: max14
    INTEGER :: max15
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    INTEGER :: min7
    INTEGER :: min8
    INTEGER :: min9
    INTEGER :: min10
    INTEGER :: min11
    INTEGER :: min12
    INTEGER :: min13
    INTEGER :: min14
    INTEGER :: min15
    REAL :: result1
    REAL :: result2
    REAL :: result3
    INTEGER :: ad_from
    INTEGER :: ad_from0
    INTEGER :: ad_from1
    INTEGER :: ad_from2
    INTEGER :: ad_from3

    qx = 0.0
    qy = 0.0
    qxx = 0.0
    qyy = 0.0
    g_in = 0.0
    g_ou = 0.0
    p0 = 0.0
    q1 = 0.0
    q2 = 0.0
    result1 = 0.0
    result2 = 0.0
    result3 = 0.0
    i = 0
    j = 0
    is1 = 0
    js1 = 0
    is2 = 0
    js2 = 0
    ie1 = 0
    je1 = 0
    max1 = 0
    max2 = 0
    max3 = 0
    max4 = 0
    max5 = 0
    max6 = 0
    max7 = 0
    max8 = 0
    max9 = 0
    max10 = 0
    max11 = 0
    max12 = 0
    max13 = 0
    max14 = 0
    max15 = 0
    min1 = 0
    min2 = 0
    min3 = 0
    min4 = 0
    min5 = 0
    min6 = 0
    min7 = 0
    min8 = 0
    min9 = 0
    min10 = 0
    min11 = 0
    min12 = 0
    min13 = 0
    min14 = 0
    min15 = 0
    ad_from = 0
    ad_from0 = 0
    ad_from1 = 0
    ad_from2 = 0
    ad_from3 = 0

    edge_w => gridstruct%edge_w
    edge_e => gridstruct%edge_e
    edge_s => gridstruct%edge_s
    edge_n => gridstruct%edge_n
    grid => gridstruct%grid
    agrid => gridstruct%agrid
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    IF (gridstruct%grid_type .LT. 3) THEN
      IF (1 .LT. is - 1) THEN
        is1 = is - 1
      ELSE
        is1 = 1
      END IF
      IF (1 .LT. js - 1) THEN
        js1 = js - 1
      ELSE
        js1 = 1
      END IF
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (2 .LT. js) THEN
        js2 = js
      ELSE
        js2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        je1 = je + 1
      ELSE
        je1 = npy - 1
      END IF
! Corners:
! 3-way extrapolation
      IF (gridstruct%nested) THEN
        DO j=js-2,je+2
          DO i=is,ie+1
            qx(i, j) = b2*(qin(i-2, j)+qin(i+1, j)) + b1*(qin(i-1, j)+&
&             qin(i, j))
          END DO
        END DO
        CALL PUSHCONTROL(2,0)
      ELSE
        IF (gridstruct%sw_corner) THEN
          p0(1:2) = grid(1, 1, 1:2)
          result1 = EXTRAP_CORNER_FWD(p0, agrid(1, 1, 1:2), agrid(2, &
&           2, 1:2), qin(1, 1), qin(2, 2))
          result2 = EXTRAP_CORNER_FWD(p0, agrid(0, 1, 1:2), agrid(-1&
&           , 2, 1:2), qin(0, 1), qin(-1, 2))
          result3 = EXTRAP_CORNER_FWD(p0, agrid(1, 0, 1:2), agrid(2, &
&           -1, 1:2), qin(1, 0), qin(2, -1))
          CALL PUSHREALARRAY(qout(1, 1))
          qout(1, 1) = (result1+result2+result3)*r3
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%se_corner) THEN
          p0(1:2) = grid(npx, 1, 1:2)
          result1 = EXTRAP_CORNER_FWD(p0, agrid(npx-1, 1, 1:2), agrid&
&           (npx-2, 2, 1:2), qin(npx-1, 1), qin(npx-2, 2))
          result2 = EXTRAP_CORNER_FWD(p0, agrid(npx-1, 0, 1:2), agrid&
&           (npx-2, -1, 1:2), qin(npx-1, 0), qin(npx-2, -1))
          result3 = EXTRAP_CORNER_FWD(p0, agrid(npx, 1, 1:2), agrid(&
&           npx+1, 2, 1:2), qin(npx, 1), qin(npx+1, 2))
          CALL PUSHREALARRAY(qout(npx, 1))
          qout(npx, 1) = (result1+result2+result3)*r3
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%ne_corner) THEN
          p0(1:2) = grid(npx, npy, 1:2)
          result1 = EXTRAP_CORNER_FWD(p0, agrid(npx-1, npy-1, 1:2), &
&           agrid(npx-2, npy-2, 1:2), qin(npx-1, npy-1), qin(npx-2, npy-&
&           2))
          result2 = EXTRAP_CORNER_FWD(p0, agrid(npx, npy-1, 1:2), &
&           agrid(npx+1, npy-2, 1:2), qin(npx, npy-1), qin(npx+1, npy-2)&
&           )
          result3 = EXTRAP_CORNER_FWD(p0, agrid(npx-1, npy, 1:2), &
&           agrid(npx-2, npy+1, 1:2), qin(npx-1, npy), qin(npx-2, npy+1)&
&           )
          CALL PUSHREALARRAY(qout(npx, npy))
          qout(npx, npy) = (result1+result2+result3)*r3
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%nw_corner) THEN
          p0(1:2) = grid(1, npy, 1:2)
          result1 = EXTRAP_CORNER_FWD(p0, agrid(1, npy-1, 1:2), agrid&
&           (2, npy-2, 1:2), qin(1, npy-1), qin(2, npy-2))
          result2 = EXTRAP_CORNER_FWD(p0, agrid(0, npy-1, 1:2), agrid&
&           (-1, npy-2, 1:2), qin(0, npy-1), qin(-1, npy-2))
          result3 = EXTRAP_CORNER_FWD(p0, agrid(1, npy, 1:2), agrid(2&
&           , npy+1, 1:2), qin(1, npy), qin(2, npy+1))
          CALL PUSHREALARRAY(qout(1, npy))
          qout(1, npy) = (result1+result2+result3)*r3
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        IF (1 .LT. js - 2) THEN
          max1 = js - 2
        ELSE
          max1 = 1
        END IF
        IF (npy - 1 .GT. je + 2) THEN
          min1 = je + 2
        ELSE
          min1 = npy - 1
        END IF
!------------
! X-Interior:
!------------
        DO j=max1,min1
          IF (3 .LT. is) THEN
            max2 = is
          ELSE
            max2 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min2 = ie + 1
          ELSE
            min2 = npx - 2
          END IF
          ad_from = max2
          DO i=ad_from,min2
            qx(i, j) = b2*(qin(i-2, j)+qin(i+1, j)) + b1*(qin(i-1, j)+&
&             qin(i, j))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from)
        END DO
! *** West Edges:
        IF (is .EQ. 1) THEN
          DO j=js1,je1
            q2(j) = (qin(0, j)*dxa(1, j)+qin(1, j)*dxa(0, j))/(dxa(0, j)&
&             +dxa(1, j))
          END DO
          DO j=js2,je1
            CALL PUSHREALARRAY(qout(1, j))
            qout(1, j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
          END DO
          IF (1 .LT. js - 2) THEN
            max3 = js - 2
          ELSE
            max3 = 1
          END IF
          IF (npy - 1 .GT. je + 2) THEN
            min3 = je + 2
          ELSE
            min3 = npy - 1
          END IF
!
          DO j=max3,min3
            g_in = dxa(2, j)/dxa(1, j)
            g_ou = dxa(-1, j)/dxa(0, j)
            qx(1, j) = 0.5*(((2.+g_in)*qin(1, j)-qin(2, j))/(1.+g_in)+((&
&             2.+g_ou)*qin(0, j)-qin(-1, j))/(1.+g_ou))
            qx(2, j) = (3.*(g_in*qin(1, j)+qin(2, j))-(g_in*qx(1, j)+qx(&
&             3, j)))/(2.+2.*g_in)
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! East Edges:
        IF (ie + 1 .EQ. npx) THEN
          DO j=js1,je1
            q2(j) = (qin(npx-1, j)*dxa(npx, j)+qin(npx, j)*dxa(npx-1, j)&
&             )/(dxa(npx-1, j)+dxa(npx, j))
          END DO
          DO j=js2,je1
            CALL PUSHREALARRAY(qout(npx, j))
            qout(npx, j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
          END DO
          IF (1 .LT. js - 2) THEN
            max4 = js - 2
          ELSE
            max4 = 1
          END IF
          IF (npy - 1 .GT. je + 2) THEN
            min4 = je + 2
          ELSE
            min4 = npy - 1
          END IF
!
          DO j=max4,min4
            g_in = dxa(npx-2, j)/dxa(npx-1, j)
            g_ou = dxa(npx+1, j)/dxa(npx, j)
            qx(npx, j) = 0.5*(((2.+g_in)*qin(npx-1, j)-qin(npx-2, j))/(&
&             1.+g_in)+((2.+g_ou)*qin(npx, j)-qin(npx+1, j))/(1.+g_ou))
            qx(npx-1, j) = (3.*(qin(npx-2, j)+g_in*qin(npx-1, j))-(g_in*&
&             qx(npx, j)+qx(npx-2, j)))/(2.+2.*g_in)
          END DO
          CALL PUSHCONTROL(2,1)
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
      END IF
!------------
! Y-Interior:
!------------
      IF (gridstruct%nested) THEN
        DO j=js,je+1
          DO i=is-2,ie+2
            qy(i, j) = b2*(qin(i, j-2)+qin(i, j+1)) + b1*(qin(i, j-1)+&
&             qin(i, j))
          END DO
        END DO
        CALL PUSHCONTROL(2,0)
      ELSE
        IF (3 .LT. js) THEN
          max5 = js
        ELSE
          max5 = 3
        END IF
        IF (npy - 2 .GT. je + 1) THEN
          min5 = je + 1
        ELSE
          min5 = npy - 2
        END IF
        DO j=max5,min5
          IF (1 .LT. is - 2) THEN
            max6 = is - 2
          ELSE
            max6 = 1
          END IF
          IF (npx - 1 .GT. ie + 2) THEN
            min6 = ie + 2
          ELSE
            min6 = npx - 1
          END IF
          ad_from0 = max6
          DO i=ad_from0,min6
            qy(i, j) = b2*(qin(i, j-2)+qin(i, j+1)) + b1*(qin(i, j-1)+&
&             qin(i, j))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from0)
        END DO
! South Edges:
        IF (js .EQ. 1) THEN
          DO i=is1,ie1
            q1(i) = (qin(i, 0)*dya(i, 1)+qin(i, 1)*dya(i, 0))/(dya(i, 0)&
&             +dya(i, 1))
          END DO
          DO i=is2,ie1
            CALL PUSHREALARRAY(qout(i, 1))
            qout(i, 1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
          END DO
          IF (1 .LT. is - 2) THEN
            max7 = is - 2
          ELSE
            max7 = 1
          END IF
          IF (npx - 1 .GT. ie + 2) THEN
            min7 = ie + 2
          ELSE
            min7 = npx - 1
          END IF
!
          DO i=max7,min7
            g_in = dya(i, 2)/dya(i, 1)
            g_ou = dya(i, -1)/dya(i, 0)
            qy(i, 1) = 0.5*(((2.+g_in)*qin(i, 1)-qin(i, 2))/(1.+g_in)+((&
&             2.+g_ou)*qin(i, 0)-qin(i, -1))/(1.+g_ou))
            qy(i, 2) = (3.*(g_in*qin(i, 1)+qin(i, 2))-(g_in*qy(i, 1)+qy(&
&             i, 3)))/(2.+2.*g_in)
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! North Edges:
        IF (je + 1 .EQ. npy) THEN
          DO i=is1,ie1
            q1(i) = (qin(i, npy-1)*dya(i, npy)+qin(i, npy)*dya(i, npy-1)&
&             )/(dya(i, npy-1)+dya(i, npy))
          END DO
          DO i=is2,ie1
            CALL PUSHREALARRAY(qout(i, npy))
            qout(i, npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
          END DO
          IF (1 .LT. is - 2) THEN
            max8 = is - 2
          ELSE
            max8 = 1
          END IF
          IF (npx - 1 .GT. ie + 2) THEN
            min8 = ie + 2
          ELSE
            min8 = npx - 1
          END IF
!
          DO i=max8,min8
            g_in = dya(i, npy-2)/dya(i, npy-1)
            g_ou = dya(i, npy+1)/dya(i, npy)
            qy(i, npy) = 0.5*(((2.+g_in)*qin(i, npy-1)-qin(i, npy-2))/(&
&             1.+g_in)+((2.+g_ou)*qin(i, npy)-qin(i, npy+1))/(1.+g_ou))
            qy(i, npy-1) = (3.*(qin(i, npy-2)+g_in*qin(i, npy-1))-(g_in*&
&             qy(i, npy)+qy(i, npy-2)))/(2.+2.*g_in)
          END DO
          CALL PUSHCONTROL(2,1)
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
      END IF
!--------------------------------------
      IF (gridstruct%nested) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qxx(i, j) = a2*(qx(i, j-2)+qx(i, j+1)) + a1*(qx(i, j-1)+qx(i&
&             , j))
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            qyy(i, j) = a2*(qy(i-2, j)+qy(i+1, j)) + a1*(qy(i-1, j)+qy(i&
&             , j))
          END DO
          DO i=is,ie+1
! averaging
            CALL PUSHREALARRAY(qout(i, j))
            qout(i, j) = 0.5*(qxx(i, j)+qyy(i, j))
          END DO
        END DO
        CALL PUSHCONTROL(2,0)
      ELSE
        IF (3 .LT. js) THEN
          max9 = js
        ELSE
          max9 = 3
        END IF
        IF (npy - 2 .GT. je + 1) THEN
          min9 = je + 1
        ELSE
          min9 = npy - 2
        END IF
        DO j=max9,min9
          IF (2 .LT. is) THEN
            max10 = is
          ELSE
            max10 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min10 = ie + 1
          ELSE
            min10 = npx - 1
          END IF
          ad_from1 = max10
          DO i=ad_from1,min10
            qxx(i, j) = a2*(qx(i, j-2)+qx(i, j+1)) + a1*(qx(i, j-1)+qx(i&
&             , j))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from1)
        END DO
        IF (js .EQ. 1) THEN
          IF (2 .LT. is) THEN
            max11 = is
          ELSE
            max11 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min11 = ie + 1
          ELSE
            min11 = npx - 1
          END IF
          DO i=max11,min11
            qxx(i, 2) = c1*(qx(i, 1)+qx(i, 2)) + c2*(qout(i, 1)+qxx(i, 3&
&             ))
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (je + 1 .EQ. npy) THEN
          IF (2 .LT. is) THEN
            max12 = is
          ELSE
            max12 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min12 = ie + 1
          ELSE
            min12 = npx - 1
          END IF
          DO i=max12,min12
            qxx(i, npy-1) = c1*(qx(i, npy-2)+qx(i, npy-1)) + c2*(qout(i&
&             , npy)+qxx(i, npy-2))
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        IF (2 .LT. js) THEN
          max13 = js
        ELSE
          max13 = 2
        END IF
        IF (npy - 1 .GT. je + 1) THEN
          min13 = je + 1
        ELSE
          min13 = npy - 1
        END IF
        DO j=max13,min13
          IF (3 .LT. is) THEN
            max14 = is
          ELSE
            max14 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min14 = ie + 1
          ELSE
            min14 = npx - 2
          END IF
          ad_from2 = max14
          DO i=ad_from2,min14
            qyy(i, j) = a2*(qy(i-2, j)+qy(i+1, j)) + a1*(qy(i-1, j)+qy(i&
&             , j))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from2)
          IF (is .EQ. 1) THEN
            qyy(2, j) = c1*(qy(1, j)+qy(2, j)) + c2*(qout(1, j)+qyy(3, j&
&             ))
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            qyy(npx-1, j) = c1*(qy(npx-2, j)+qy(npx-1, j)) + c2*(qout(&
&             npx, j)+qyy(npx-2, j))
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
          IF (2 .LT. is) THEN
            max15 = is
          ELSE
            max15 = 2
          END IF
          IF (npx - 1 .GT. ie + 1) THEN
            min15 = ie + 1
          ELSE
            min15 = npx - 1
          END IF
          ad_from3 = max15
          DO i=ad_from3,min15
! averaging
            CALL PUSHREALARRAY(qout(i, j))
            qout(i, j) = 0.5*(qxx(i, j)+qyy(i, j))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from3)
        END DO
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
! grid_type>=3
!------------------------
! Doubly periodic domain:
!------------------------
! X-sweep: PPM
      DO j=js-2,je+2
        DO i=is,ie+1
          qx(i, j) = b1*(qin(i-1, j)+qin(i, j)) + b2*(qin(i-2, j)+qin(i+&
&           1, j))
        END DO
      END DO
! Y-sweep: PPM
      DO j=js,je+1
        DO i=is-2,ie+2
          qy(i, j) = b1*(qin(i, j-1)+qin(i, j)) + b2*(qin(i, j-2)+qin(i&
&           , j+1))
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY(qout(i, j))
          qout(i, j) = 0.5*(a1*(qx(i, j-1)+qx(i, j)+qy(i-1, j)+qy(i, j))&
&           +a2*(qx(i, j-2)+qx(i, j+1)+qy(i-2, j)+qy(i+1, j)))
        END DO
      END DO
      CALL PUSHCONTROL(2,2)
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            CALL PUSHREALARRAY(qin(i, j))
            qin(i, j) = qout(i, j)
          END DO
        END DO
        CALL PUSHINTEGER(min9)
        CALL PUSHINTEGER(min8)
        CALL PUSHINTEGER(min7)
        CALL PUSHINTEGER(min5)
        CALL PUSHINTEGER(min4)
        CALL PUSHINTEGER(ie1)
        CALL PUSHINTEGER(min3)
        CALL PUSHINTEGER(min1)
        CALL PUSHINTEGER(js2)
        CALL PUSHINTEGER(js1)
        CALL PUSHINTEGER(je1)
        CALL PUSHINTEGER(min13)
        CALL PUSHINTEGER(min12)
        CALL PUSHINTEGER(min11)
        CALL PUSHINTEGER(is2)
        CALL PUSHINTEGER(is1)
        !CALL PUSHPOINTER8(C_LOC(agrid))
        CALL PUSHINTEGER(max9)
        CALL PUSHINTEGER(max8)
        CALL PUSHINTEGER(max7)
        CALL PUSHINTEGER(max5)
        CALL PUSHINTEGER(max4)
        CALL PUSHINTEGER(max3)
        CALL PUSHINTEGER(max13)
        CALL PUSHINTEGER(max1)
        CALL PUSHINTEGER(max12)
        CALL PUSHINTEGER(max11)
        CALL PUSHCONTROL(2,2)
      ELSE
        CALL PUSHINTEGER(min9)
        CALL PUSHINTEGER(min8)
        CALL PUSHINTEGER(min7)
        CALL PUSHINTEGER(min5)
        CALL PUSHINTEGER(min4)
        CALL PUSHINTEGER(ie1)
        CALL PUSHINTEGER(min3)
        CALL PUSHINTEGER(min1)
        CALL PUSHINTEGER(js2)
        CALL PUSHINTEGER(js1)
        CALL PUSHINTEGER(je1)
        CALL PUSHINTEGER(min13)
        CALL PUSHINTEGER(min12)
        CALL PUSHINTEGER(min11)
        CALL PUSHINTEGER(is2)
        CALL PUSHINTEGER(is1)
        !CALL PUSHPOINTER8(C_LOC(agrid))
        CALL PUSHINTEGER(max9)
        CALL PUSHINTEGER(max8)
        CALL PUSHINTEGER(max7)
        CALL PUSHINTEGER(max5)
        CALL PUSHINTEGER(max4)
        CALL PUSHINTEGER(max3)
        CALL PUSHINTEGER(max13)
        CALL PUSHINTEGER(max1)
        CALL PUSHINTEGER(max12)
        CALL PUSHINTEGER(max11)
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHINTEGER(min9)
      CALL PUSHINTEGER(min8)
      CALL PUSHINTEGER(min7)
      CALL PUSHINTEGER(min5)
      CALL PUSHINTEGER(min4)
      CALL PUSHINTEGER(ie1)
      CALL PUSHINTEGER(min3)
      CALL PUSHINTEGER(min1)
      CALL PUSHINTEGER(js2)
      CALL PUSHINTEGER(js1)
      CALL PUSHINTEGER(je1)
      CALL PUSHINTEGER(min13)
      CALL PUSHINTEGER(min12)
      CALL PUSHINTEGER(min11)
      CALL PUSHINTEGER(is2)
      CALL PUSHINTEGER(is1)
      !CALL PUSHPOINTER8(C_LOC(agrid))
      CALL PUSHINTEGER(max9)
      CALL PUSHINTEGER(max8)
      CALL PUSHINTEGER(max7)
      CALL PUSHINTEGER(max5)
      CALL PUSHINTEGER(max4)
      CALL PUSHINTEGER(max3)
      CALL PUSHINTEGER(max13)
      CALL PUSHINTEGER(max1)
      CALL PUSHINTEGER(max12)
      CALL PUSHINTEGER(max11)
      CALL PUSHCONTROL(2,0)
    END IF
  END SUBROUTINE A2B_ORD4_FWD
!  Differentiation of a2b_ord4 in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: qin qout
!   with respect to varying inputs: qin qout
  SUBROUTINE A2B_ORD4_BWD(qin, qin_ad, qout, qout_ad, gridstruct, npx&
&   , npy, is, ie, js, je, ng, replace)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
    REAL, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qin_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qout_ad(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
    REAL, PARAMETER :: c1=2./3.
    REAL, PARAMETER :: c2=-(1./6.)
    REAL :: qx(is:ie+1, js-ng:je+ng)
    REAL :: qx_ad(is:ie+1, js-ng:je+ng)
    REAL :: qy(is-ng:ie+ng, js:je+1)
    REAL :: qy_ad(is-ng:ie+ng, js:je+1)
    REAL :: qxx(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qxx_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL :: g_in, g_ou
    REAL :: p0(2)
    REAL :: q1(is-1:ie+1), q2(js-1:je+1)
    REAL :: q1_ad(is-1:ie+1), q2_ad(js-1:je+1)
    INTEGER :: i, j, is1, js1, is2, js2, ie1, je1
    REAL, DIMENSION(:, :, :), POINTER :: grid, agrid
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
    REAL(kind=r_grid), DIMENSION(:), POINTER :: edge_w, edge_e, edge_s, &
&   edge_n
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: max7
    INTEGER :: max8
    INTEGER :: max9
    INTEGER :: max10
    INTEGER :: max11
    INTEGER :: max12
    INTEGER :: max13
    INTEGER :: max14
    INTEGER :: max15
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    INTEGER :: min7
    INTEGER :: min8
    INTEGER :: min9
    INTEGER :: min10
    INTEGER :: min11
    INTEGER :: min12
    INTEGER :: min13
    INTEGER :: min14
    INTEGER :: min15
    REAL :: result1
    REAL :: result1_ad
    REAL :: result2
    REAL :: result2_ad
    REAL :: result3
    REAL :: result3_ad
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
    REAL :: temp_ad22
    REAL :: temp_ad23
    REAL :: temp_ad24
    REAL :: temp_ad25
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
    INTEGER :: branch
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_a2b_ord4

    qx = 0.0
    qy = 0.0
    qxx = 0.0
    qyy = 0.0
    g_in = 0.0
    g_ou = 0.0
    p0 = 0.0
    q1 = 0.0
    q2 = 0.0
    result1 = 0.0
    result2 = 0.0
    result3 = 0.0
    i = 0
    j = 0
    is1 = 0
    js1 = 0
    is2 = 0
    js2 = 0
    ie1 = 0
    je1 = 0
    max1 = 0
    max2 = 0
    max3 = 0
    max4 = 0
    max5 = 0
    max6 = 0
    max7 = 0
    max8 = 0
    max9 = 0
    max10 = 0
    max11 = 0
    max12 = 0
    max13 = 0
    max14 = 0
    max15 = 0
    min1 = 0
    min2 = 0
    min3 = 0
    min4 = 0
    min5 = 0
    min6 = 0
    min7 = 0
    min8 = 0
    min9 = 0
    min10 = 0
    min11 = 0
    min12 = 0
    min13 = 0
    min14 = 0
    min15 = 0
    ad_from = 0
    ad_from0 = 0
    ad_from1 = 0
    ad_from2 = 0
    ad_from3 = 0
    ad_to = 0
    ad_to0 = 0
    ad_to1 = 0
    ad_to2 = 0
    ad_to3 = 0
    branch = 0

    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(max11)
      CALL POPINTEGER(max12)
      CALL POPINTEGER(max1)
      CALL POPINTEGER(max13)
      CALL POPINTEGER(max3)
      CALL POPINTEGER(max4)
      CALL POPINTEGER(max5)
      CALL POPINTEGER(max7)
      CALL POPINTEGER(max8)
      CALL POPINTEGER(max9)
      !CALL POPPOINTER8(cptr)
      agrid => gridstruct%agrid ! (/unknown_shape_in_a2b_ord4/))
      CALL POPINTEGER(is1)
      CALL POPINTEGER(is2)
      CALL POPINTEGER(min11)
      CALL POPINTEGER(min12)
      CALL POPINTEGER(min13)
      CALL POPINTEGER(je1)
      CALL POPINTEGER(js1)
      CALL POPINTEGER(js2)
      CALL POPINTEGER(min1)
      CALL POPINTEGER(min3)
      CALL POPINTEGER(ie1)
      CALL POPINTEGER(min4)
      CALL POPINTEGER(min5)
      CALL POPINTEGER(min7)
      CALL POPINTEGER(min8)
      CALL POPINTEGER(min9)
    ELSE IF (branch .EQ. 1) THEN
      CALL POPINTEGER(max11)
      CALL POPINTEGER(max12)
      CALL POPINTEGER(max1)
      CALL POPINTEGER(max13)
      CALL POPINTEGER(max3)
      CALL POPINTEGER(max4)
      CALL POPINTEGER(max5)
      CALL POPINTEGER(max7)
      CALL POPINTEGER(max8)
      CALL POPINTEGER(max9)
      !CALL POPPOINTER8(cptr)
      agrid => gridstruct%agrid ! (/unknown_shape_in_a2b_ord4/))
      CALL POPINTEGER(is1)
      CALL POPINTEGER(is2)
      CALL POPINTEGER(min11)
      CALL POPINTEGER(min12)
      CALL POPINTEGER(min13)
      CALL POPINTEGER(je1)
      CALL POPINTEGER(js1)
      CALL POPINTEGER(js2)
      CALL POPINTEGER(min1)
      CALL POPINTEGER(min3)
      CALL POPINTEGER(ie1)
      CALL POPINTEGER(min4)
      CALL POPINTEGER(min5)
      CALL POPINTEGER(min7)
      CALL POPINTEGER(min8)
      CALL POPINTEGER(min9)
    ELSE
      CALL POPINTEGER(max11)
      CALL POPINTEGER(max12)
      CALL POPINTEGER(max1)
      CALL POPINTEGER(max13)
      CALL POPINTEGER(max3)
      CALL POPINTEGER(max4)
      CALL POPINTEGER(max5)
      CALL POPINTEGER(max7)
      CALL POPINTEGER(max8)
      CALL POPINTEGER(max9)
      !CALL POPPOINTER8(cptr)
      agrid => gridstruct%agrid ! (/unknown_shape_in_a2b_ord4/))
      CALL POPINTEGER(is1)
      CALL POPINTEGER(is2)
      CALL POPINTEGER(min11)
      CALL POPINTEGER(min12)
      CALL POPINTEGER(min13)
      CALL POPINTEGER(je1)
      CALL POPINTEGER(js1)
      CALL POPINTEGER(js2)
      CALL POPINTEGER(min1)
      CALL POPINTEGER(min3)
      CALL POPINTEGER(ie1)
      CALL POPINTEGER(min4)
      CALL POPINTEGER(min5)
      CALL POPINTEGER(min7)
      CALL POPINTEGER(min8)
      CALL POPINTEGER(min9)
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(qin(i, j))
          qout_ad(i, j) = qout_ad(i, j) + qin_ad(i, j)
          qin_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      qy_ad = 0.0
      qxx_ad = 0.0
      qyy_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(qout(i, j))
          qxx_ad(i, j) = qxx_ad(i, j) + 0.5*qout_ad(i, j)
          qyy_ad(i, j) = qyy_ad(i, j) + 0.5*qout_ad(i, j)
          qout_ad(i, j) = 0.0
        END DO
        DO i=ie+1,is,-1
          qy_ad(i-2, j) = qy_ad(i-2, j) + a2*qyy_ad(i, j)
          qy_ad(i+1, j) = qy_ad(i+1, j) + a2*qyy_ad(i, j)
          qy_ad(i-1, j) = qy_ad(i-1, j) + a1*qyy_ad(i, j)
          qy_ad(i, j) = qy_ad(i, j) + a1*qyy_ad(i, j)
          qyy_ad(i, j) = 0.0
        END DO
      END DO
      qx_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          qx_ad(i, j-2) = qx_ad(i, j-2) + a2*qxx_ad(i, j)
          qx_ad(i, j+1) = qx_ad(i, j+1) + a2*qxx_ad(i, j)
          qx_ad(i, j-1) = qx_ad(i, j-1) + a1*qxx_ad(i, j)
          qx_ad(i, j) = qx_ad(i, j) + a1*qxx_ad(i, j)
          qxx_ad(i, j) = 0.0
        END DO
      END DO
    ELSE IF (branch .EQ. 1) THEN
      qy_ad = 0.0
      qxx_ad = 0.0
      qyy_ad = 0.0
      DO j=min13,max13,-1
        CALL POPINTEGER(ad_from3)
        CALL POPINTEGER(ad_to3)
        DO i=ad_to3,ad_from3,-1
          CALL POPREALARRAY(qout(i, j))
          qxx_ad(i, j) = qxx_ad(i, j) + 0.5*qout_ad(i, j)
          qyy_ad(i, j) = qyy_ad(i, j) + 0.5*qout_ad(i, j)
          qout_ad(i, j) = 0.0
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          qy_ad(npx-2, j) = qy_ad(npx-2, j) + c1*qyy_ad(npx-1, j)
          qy_ad(npx-1, j) = qy_ad(npx-1, j) + c1*qyy_ad(npx-1, j)
          qout_ad(npx, j) = qout_ad(npx, j) + c2*qyy_ad(npx-1, j)
          qyy_ad(npx-2, j) = qyy_ad(npx-2, j) + c2*qyy_ad(npx-1, j)
          qyy_ad(npx-1, j) = 0.0
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          qy_ad(1, j) = qy_ad(1, j) + c1*qyy_ad(2, j)
          qy_ad(2, j) = qy_ad(2, j) + c1*qyy_ad(2, j)
          qout_ad(1, j) = qout_ad(1, j) + c2*qyy_ad(2, j)
          qyy_ad(3, j) = qyy_ad(3, j) + c2*qyy_ad(2, j)
          qyy_ad(2, j) = 0.0
        END IF
        CALL POPINTEGER(ad_from2)
        CALL POPINTEGER(ad_to2)
        DO i=ad_to2,ad_from2,-1
          qy_ad(i-2, j) = qy_ad(i-2, j) + a2*qyy_ad(i, j)
          qy_ad(i+1, j) = qy_ad(i+1, j) + a2*qyy_ad(i, j)
          qy_ad(i-1, j) = qy_ad(i-1, j) + a1*qyy_ad(i, j)
          qy_ad(i, j) = qy_ad(i, j) + a1*qyy_ad(i, j)
          qyy_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        qx_ad = 0.0
      ELSE
        qx_ad = 0.0
        DO i=min12,max12,-1
          qx_ad(i, npy-2) = qx_ad(i, npy-2) + c1*qxx_ad(i, npy-1)
          qx_ad(i, npy-1) = qx_ad(i, npy-1) + c1*qxx_ad(i, npy-1)
          qout_ad(i, npy) = qout_ad(i, npy) + c2*qxx_ad(i, npy-1)
          qxx_ad(i, npy-2) = qxx_ad(i, npy-2) + c2*qxx_ad(i, npy-1)
          qxx_ad(i, npy-1) = 0.0
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO i=min11,max11,-1
          qx_ad(i, 1) = qx_ad(i, 1) + c1*qxx_ad(i, 2)
          qx_ad(i, 2) = qx_ad(i, 2) + c1*qxx_ad(i, 2)
          qout_ad(i, 1) = qout_ad(i, 1) + c2*qxx_ad(i, 2)
          qxx_ad(i, 3) = qxx_ad(i, 3) + c2*qxx_ad(i, 2)
          qxx_ad(i, 2) = 0.0
        END DO
      END IF
      DO j=min9,max9,-1
        CALL POPINTEGER(ad_from1)
        CALL POPINTEGER(ad_to1)
        DO i=ad_to1,ad_from1,-1
          qx_ad(i, j-2) = qx_ad(i, j-2) + a2*qxx_ad(i, j)
          qx_ad(i, j+1) = qx_ad(i, j+1) + a2*qxx_ad(i, j)
          qx_ad(i, j-1) = qx_ad(i, j-1) + a1*qxx_ad(i, j)
          qx_ad(i, j) = qx_ad(i, j) + a1*qxx_ad(i, j)
          qxx_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      qx_ad = 0.0
      qy_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(qout(i, j))
          temp_ad23 = 0.5*qout_ad(i, j)
          temp_ad24 = a1*temp_ad23
          temp_ad25 = a2*temp_ad23
          qx_ad(i, j-1) = qx_ad(i, j-1) + temp_ad24
          qx_ad(i, j) = qx_ad(i, j) + temp_ad24
          qy_ad(i-1, j) = qy_ad(i-1, j) + temp_ad24
          qy_ad(i, j) = qy_ad(i, j) + temp_ad24
          qx_ad(i, j-2) = qx_ad(i, j-2) + temp_ad25
          qx_ad(i, j+1) = qx_ad(i, j+1) + temp_ad25
          qy_ad(i-2, j) = qy_ad(i-2, j) + temp_ad25
          qy_ad(i+1, j) = qy_ad(i+1, j) + temp_ad25
          qout_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+2,is-2,-1
          qin_ad(i, j-1) = qin_ad(i, j-1) + b1*qy_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qy_ad(i, j)
          qin_ad(i, j-2) = qin_ad(i, j-2) + b2*qy_ad(i, j)
          qin_ad(i, j+1) = qin_ad(i, j+1) + b2*qy_ad(i, j)
          qy_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+2,js-2,-1
        DO i=ie+1,is,-1
          qin_ad(i-1, j) = qin_ad(i-1, j) + b1*qx_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qx_ad(i, j)
          qin_ad(i-2, j) = qin_ad(i-2, j) + b2*qx_ad(i, j)
          qin_ad(i+1, j) = qin_ad(i+1, j) + b2*qx_ad(i, j)
          qx_ad(i, j) = 0.0
        END DO
      END DO
      GOTO 100
    END IF
    dya => gridstruct%dya
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      DO j=je+1,js,-1
        DO i=ie+2,is-2,-1
          qin_ad(i, j-2) = qin_ad(i, j-2) + b2*qy_ad(i, j)
          qin_ad(i, j+1) = qin_ad(i, j+1) + b2*qy_ad(i, j)
          qin_ad(i, j-1) = qin_ad(i, j-1) + b1*qy_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qy_ad(i, j)
          qy_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      IF (branch .EQ. 1) THEN
        DO i=min8,max8,-1
          g_in = dya(i, npy-2)/dya(i, npy-1)
          temp_ad19 = qy_ad(i, npy-1)/(g_in*2.+2.)
          qin_ad(i, npy-2) = qin_ad(i, npy-2) + 3.*temp_ad19
          qy_ad(i, npy) = qy_ad(i, npy) - g_in*temp_ad19
          qy_ad(i, npy-2) = qy_ad(i, npy-2) - temp_ad19
          qy_ad(i, npy-1) = 0.0
          g_ou = dya(i, npy+1)/dya(i, npy)
          temp_ad21 = 0.5*qy_ad(i, npy)
          temp_ad20 = temp_ad21/(g_in+1.)
          qin_ad(i, npy-1) = qin_ad(i, npy-1) + (g_in+2.)*temp_ad20 + 3.&
&           *g_in*temp_ad19
          temp_ad22 = temp_ad21/(g_ou+1.)
          qin_ad(i, npy-2) = qin_ad(i, npy-2) - temp_ad20
          qin_ad(i, npy) = qin_ad(i, npy) + (g_ou+2.)*temp_ad22
          qin_ad(i, npy+1) = qin_ad(i, npy+1) - temp_ad22
          qy_ad(i, npy) = 0.0
        END DO
        edge_n => gridstruct%edge_n
        q1_ad = 0.0
        DO i=ie1,is2,-1
          CALL POPREALARRAY(qout(i, npy))
          q1_ad(i-1) = q1_ad(i-1) + edge_n(i)*qout_ad(i, npy)
          q1_ad(i) = q1_ad(i) + (1.-edge_n(i))*qout_ad(i, npy)
          qout_ad(i, npy) = 0.0
        END DO
        DO i=ie1,is1,-1
          temp_ad18 = q1_ad(i)/(dya(i, npy-1)+dya(i, npy))
          qin_ad(i, npy-1) = qin_ad(i, npy-1) + dya(i, npy)*temp_ad18
          qin_ad(i, npy) = qin_ad(i, npy) + dya(i, npy-1)*temp_ad18
          q1_ad(i) = 0.0
        END DO
      ELSE
        q1_ad = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO i=min7,max7,-1
          g_in = dya(i, 2)/dya(i, 1)
          temp_ad14 = qy_ad(i, 2)/(g_in*2.+2.)
          qin_ad(i, 1) = qin_ad(i, 1) + 3.*g_in*temp_ad14
          qin_ad(i, 2) = qin_ad(i, 2) + 3.*temp_ad14
          qy_ad(i, 1) = qy_ad(i, 1) - g_in*temp_ad14
          qy_ad(i, 3) = qy_ad(i, 3) - temp_ad14
          qy_ad(i, 2) = 0.0
          g_ou = dya(i, -1)/dya(i, 0)
          temp_ad15 = 0.5*qy_ad(i, 1)
          temp_ad16 = temp_ad15/(g_in+1.)
          temp_ad17 = temp_ad15/(g_ou+1.)
          qin_ad(i, 1) = qin_ad(i, 1) + (g_in+2.)*temp_ad16
          qin_ad(i, 2) = qin_ad(i, 2) - temp_ad16
          qin_ad(i, 0) = qin_ad(i, 0) + (g_ou+2.)*temp_ad17
          qin_ad(i, -1) = qin_ad(i, -1) - temp_ad17
          qy_ad(i, 1) = 0.0
        END DO
        edge_s => gridstruct%edge_s
        DO i=ie1,is2,-1
          CALL POPREALARRAY(qout(i, 1))
          q1_ad(i-1) = q1_ad(i-1) + edge_s(i)*qout_ad(i, 1)
          q1_ad(i) = q1_ad(i) + (1.-edge_s(i))*qout_ad(i, 1)
          qout_ad(i, 1) = 0.0
        END DO
        DO i=ie1,is1,-1
          temp_ad13 = q1_ad(i)/(dya(i, 0)+dya(i, 1))
          qin_ad(i, 0) = qin_ad(i, 0) + dya(i, 1)*temp_ad13
          qin_ad(i, 1) = qin_ad(i, 1) + dya(i, 0)*temp_ad13
          q1_ad(i) = 0.0
        END DO
      END IF
      DO j=min5,max5,-1
        CALL POPINTEGER(ad_from0)
        CALL POPINTEGER(ad_to0)
        DO i=ad_to0,ad_from0,-1
          qin_ad(i, j-2) = qin_ad(i, j-2) + b2*qy_ad(i, j)
          qin_ad(i, j+1) = qin_ad(i, j+1) + b2*qy_ad(i, j)
          qin_ad(i, j-1) = qin_ad(i, j-1) + b1*qy_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qy_ad(i, j)
          qy_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    dxa => gridstruct%dxa
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      DO j=je+2,js-2,-1
        DO i=ie+1,is,-1
          qin_ad(i-2, j) = qin_ad(i-2, j) + b2*qx_ad(i, j)
          qin_ad(i+1, j) = qin_ad(i+1, j) + b2*qx_ad(i, j)
          qin_ad(i-1, j) = qin_ad(i-1, j) + b1*qx_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qx_ad(i, j)
          qx_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      IF (branch .EQ. 1) THEN
        DO j=min4,max4,-1
          g_in = dxa(npx-2, j)/dxa(npx-1, j)
          temp_ad9 = qx_ad(npx-1, j)/(g_in*2.+2.)
          qin_ad(npx-2, j) = qin_ad(npx-2, j) + 3.*temp_ad9
          qx_ad(npx, j) = qx_ad(npx, j) - g_in*temp_ad9
          qx_ad(npx-2, j) = qx_ad(npx-2, j) - temp_ad9
          qx_ad(npx-1, j) = 0.0
          g_ou = dxa(npx+1, j)/dxa(npx, j)
          temp_ad11 = 0.5*qx_ad(npx, j)
          temp_ad10 = temp_ad11/(g_in+1.)
          qin_ad(npx-1, j) = qin_ad(npx-1, j) + (g_in+2.)*temp_ad10 + 3.&
&           *g_in*temp_ad9
          temp_ad12 = temp_ad11/(g_ou+1.)
          qin_ad(npx-2, j) = qin_ad(npx-2, j) - temp_ad10
          qin_ad(npx, j) = qin_ad(npx, j) + (g_ou+2.)*temp_ad12
          qin_ad(npx+1, j) = qin_ad(npx+1, j) - temp_ad12
          qx_ad(npx, j) = 0.0
        END DO
        edge_e => gridstruct%edge_e
        q2_ad = 0.0
        DO j=je1,js2,-1
          CALL POPREALARRAY(qout(npx, j))
          q2_ad(j-1) = q2_ad(j-1) + edge_e(j)*qout_ad(npx, j)
          q2_ad(j) = q2_ad(j) + (1.-edge_e(j))*qout_ad(npx, j)
          qout_ad(npx, j) = 0.0
        END DO
        DO j=je1,js1,-1
          temp_ad8 = q2_ad(j)/(dxa(npx-1, j)+dxa(npx, j))
          qin_ad(npx-1, j) = qin_ad(npx-1, j) + dxa(npx, j)*temp_ad8
          qin_ad(npx, j) = qin_ad(npx, j) + dxa(npx-1, j)*temp_ad8
          q2_ad(j) = 0.0
        END DO
      ELSE
        q2_ad = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=min3,max3,-1
          g_in = dxa(2, j)/dxa(1, j)
          temp_ad4 = qx_ad(2, j)/(g_in*2.+2.)
          qin_ad(1, j) = qin_ad(1, j) + 3.*g_in*temp_ad4
          qin_ad(2, j) = qin_ad(2, j) + 3.*temp_ad4
          qx_ad(1, j) = qx_ad(1, j) - g_in*temp_ad4
          qx_ad(3, j) = qx_ad(3, j) - temp_ad4
          qx_ad(2, j) = 0.0
          g_ou = dxa(-1, j)/dxa(0, j)
          temp_ad5 = 0.5*qx_ad(1, j)
          temp_ad6 = temp_ad5/(g_in+1.)
          temp_ad7 = temp_ad5/(g_ou+1.)
          qin_ad(1, j) = qin_ad(1, j) + (g_in+2.)*temp_ad6
          qin_ad(2, j) = qin_ad(2, j) - temp_ad6
          qin_ad(0, j) = qin_ad(0, j) + (g_ou+2.)*temp_ad7
          qin_ad(-1, j) = qin_ad(-1, j) - temp_ad7
          qx_ad(1, j) = 0.0
        END DO
        edge_w => gridstruct%edge_w
        DO j=je1,js2,-1
          CALL POPREALARRAY(qout(1, j))
          q2_ad(j-1) = q2_ad(j-1) + edge_w(j)*qout_ad(1, j)
          q2_ad(j) = q2_ad(j) + (1.-edge_w(j))*qout_ad(1, j)
          qout_ad(1, j) = 0.0
        END DO
        DO j=je1,js1,-1
          temp_ad3 = q2_ad(j)/(dxa(0, j)+dxa(1, j))
          qin_ad(0, j) = qin_ad(0, j) + dxa(1, j)*temp_ad3
          qin_ad(1, j) = qin_ad(1, j) + dxa(0, j)*temp_ad3
          q2_ad(j) = 0.0
        END DO
      END IF
      DO j=min1,max1,-1
        CALL POPINTEGER(ad_from)
        CALL POPINTEGER(ad_to)
        DO i=ad_to,ad_from,-1
          qin_ad(i-2, j) = qin_ad(i-2, j) + b2*qx_ad(i, j)
          qin_ad(i+1, j) = qin_ad(i+1, j) + b2*qx_ad(i, j)
          qin_ad(i-1, j) = qin_ad(i-1, j) + b1*qx_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qx_ad(i, j)
          qx_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        CALL POPREALARRAY(qout(1, npy))
        temp_ad2 = r3*qout_ad(1, npy)
        result1_ad = temp_ad2
        result2_ad = temp_ad2
        result3_ad = temp_ad2
        qout_ad(1, npy) = 0.0
        CALL EXTRAP_CORNER_BWD(p0, agrid(1, npy, 1:2), agrid(2, npy+1&
&                           , 1:2), qin(1, npy), qin_ad(1, npy), qin(2, &
&                           npy+1), qin_ad(2, npy+1), result3_ad)
        CALL EXTRAP_CORNER_BWD(p0, agrid(0, npy-1, 1:2), agrid(-1, &
&                           npy-2, 1:2), qin(0, npy-1), qin_ad(0, npy-1)&
&                           , qin(-1, npy-2), qin_ad(-1, npy-2), &
&                           result2_ad)
        CALL EXTRAP_CORNER_BWD(p0, agrid(1, npy-1, 1:2), agrid(2, npy&
&                           -2, 1:2), qin(1, npy-1), qin_ad(1, npy-1), &
&                           qin(2, npy-2), qin_ad(2, npy-2), result1_ad)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(qout(npx, npy))
        temp_ad1 = r3*qout_ad(npx, npy)
        result1_ad = temp_ad1
        result2_ad = temp_ad1
        result3_ad = temp_ad1
        qout_ad(npx, npy) = 0.0
        CALL EXTRAP_CORNER_BWD(p0, agrid(npx-1, npy, 1:2), agrid(npx-&
&                           2, npy+1, 1:2), qin(npx-1, npy), qin_ad(npx-&
&                           1, npy), qin(npx-2, npy+1), qin_ad(npx-2, &
&                           npy+1), result3_ad)
        CALL EXTRAP_CORNER_BWD(p0, agrid(npx, npy-1, 1:2), agrid(npx+&
&                           1, npy-2, 1:2), qin(npx, npy-1), qin_ad(npx&
&                           , npy-1), qin(npx+1, npy-2), qin_ad(npx+1, &
&                           npy-2), result2_ad)
        CALL EXTRAP_CORNER_BWD(p0, agrid(npx-1, npy-1, 1:2), agrid(&
&                           npx-2, npy-2, 1:2), qin(npx-1, npy-1), &
&                           qin_ad(npx-1, npy-1), qin(npx-2, npy-2), &
&                           qin_ad(npx-2, npy-2), result1_ad)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(qout(npx, 1))
        temp_ad0 = r3*qout_ad(npx, 1)
        result1_ad = temp_ad0
        result2_ad = temp_ad0
        result3_ad = temp_ad0
        qout_ad(npx, 1) = 0.0
        CALL EXTRAP_CORNER_BWD(p0, agrid(npx, 1, 1:2), agrid(npx+1, 2&
&                           , 1:2), qin(npx, 1), qin_ad(npx, 1), qin(npx&
&                           +1, 2), qin_ad(npx+1, 2), result3_ad)
        CALL EXTRAP_CORNER_BWD(p0, agrid(npx-1, 0, 1:2), agrid(npx-2&
&                           , -1, 1:2), qin(npx-1, 0), qin_ad(npx-1, 0)&
&                           , qin(npx-2, -1), qin_ad(npx-2, -1), &
&                           result2_ad)
        CALL EXTRAP_CORNER_BWD(p0, agrid(npx-1, 1, 1:2), agrid(npx-2&
&                           , 2, 1:2), qin(npx-1, 1), qin_ad(npx-1, 1), &
&                           qin(npx-2, 2), qin_ad(npx-2, 2), result1_ad)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(qout(1, 1))
        temp_ad = r3*qout_ad(1, 1)
        result1_ad = temp_ad
        result2_ad = temp_ad
        result3_ad = temp_ad
        qout_ad(1, 1) = 0.0
        CALL EXTRAP_CORNER_BWD(p0, agrid(1, 0, 1:2), agrid(2, -1, 1:2&
&                           ), qin(1, 0), qin_ad(1, 0), qin(2, -1), &
&                           qin_ad(2, -1), result3_ad)
        CALL EXTRAP_CORNER_BWD(p0, agrid(0, 1, 1:2), agrid(-1, 2, 1:2&
&                           ), qin(0, 1), qin_ad(0, 1), qin(-1, 2), &
&                           qin_ad(-1, 2), result2_ad)
        CALL EXTRAP_CORNER_BWD(p0, agrid(1, 1, 1:2), agrid(2, 2, 1:2)&
&                           , qin(1, 1), qin_ad(1, 1), qin(2, 2), qin_ad&
&                           (2, 2), result1_ad)
      END IF
    END IF
 100 CONTINUE
  END SUBROUTINE A2B_ORD4_BWD
!  Differentiation of a2b_ord2 in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
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
!   gradient     of useful results: qin qout
!   with respect to varying inputs: qin qout
  SUBROUTINE A2B_ORD2_FWD(qin, qout, gridstruct, npx, npy, is, ie, js, &
&   je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL :: qout(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local:
    REAL :: q1(npx), q2(npy)
    INTEGER :: i, j
    INTEGER :: is1, js1, is2, js2, ie1, je1
    REAL, DIMENSION(:, :, :), POINTER :: grid, agrid
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
    REAL(kind=r_grid), DIMENSION(:), POINTER :: edge_w, edge_e, edge_s, &
&   edge_n
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT

    q1 = 0.0
    q2 = 0.0
    i = 0
    j = 0
    is1 = 0
    js1 = 0
    is2 = 0
    js2 = 0
    ie1 = 0
    je1 = 0

    edge_w => gridstruct%edge_w
    edge_e => gridstruct%edge_e
    edge_s => gridstruct%edge_s
    edge_n => gridstruct%edge_n
    IF (gridstruct%grid_type .LT. 3) THEN
      IF (gridstruct%nested) THEN
        DO j=js-2,je+1+2
          DO i=is-2,ie+1+2
            CALL PUSHREALARRAY(qout(i, j))
            qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin&
&             (i, j))
          END DO
        END DO
        CALL PUSHCONTROL(2,0)
      ELSE
        IF (1 .LT. is - 1) THEN
          is1 = is - 1
        ELSE
          is1 = 1
        END IF
        IF (1 .LT. js - 1) THEN
          js1 = js - 1
        ELSE
          js1 = 1
        END IF
        IF (2 .LT. is) THEN
          is2 = is
        ELSE
          is2 = 2
        END IF
        IF (2 .LT. js) THEN
          js2 = js
        ELSE
          js2 = 2
        END IF
        IF (npx - 1 .GT. ie + 1) THEN
          ie1 = ie + 1
        ELSE
          ie1 = npx - 1
        END IF
        IF (npy - 1 .GT. je + 1) THEN
          je1 = je + 1
        ELSE
          je1 = npy - 1
        END IF
        DO j=js2,je1
          DO i=is2,ie1
            CALL PUSHREALARRAY(qout(i, j))
            qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin&
&             (i, j))
          END DO
        END DO
! Fix the 4 Corners:
        IF (gridstruct%sw_corner) THEN
          CALL PUSHREALARRAY(qout(1, 1))
          qout(1, 1) = r3*(qin(1, 1)+qin(1, 0)+qin(0, 1))
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%se_corner) THEN
          CALL PUSHREALARRAY(qout(npx, 1))
          qout(npx, 1) = r3*(qin(npx-1, 1)+qin(npx-1, 0)+qin(npx, 1))
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%ne_corner) THEN
          CALL PUSHREALARRAY(qout(npx, npy))
          qout(npx, npy) = r3*(qin(npx-1, npy-1)+qin(npx, npy-1)+qin(npx&
&           -1, npy))
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (gridstruct%nw_corner) THEN
          CALL PUSHREALARRAY(qout(1, npy))
          qout(1, npy) = r3*(qin(1, npy-1)+qin(0, npy-1)+qin(1, npy))
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! *** West Edges:
        IF (is .EQ. 1) THEN
          DO j=js1,je1
            q2(j) = 0.5*(qin(0, j)+qin(1, j))
          END DO
          DO j=js2,je1
            CALL PUSHREALARRAY(qout(1, j))
            qout(1, j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! East Edges:
        IF (ie + 1 .EQ. npx) THEN
          DO j=js1,je1
            q2(j) = 0.5*(qin(npx-1, j)+qin(npx, j))
          END DO
          DO j=js2,je1
            CALL PUSHREALARRAY(qout(npx, j))
            qout(npx, j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! South Edges:
        IF (js .EQ. 1) THEN
          DO i=is1,ie1
            q1(i) = 0.5*(qin(i, 0)+qin(i, 1))
          END DO
          DO i=is2,ie1
            CALL PUSHREALARRAY(qout(i, 1))
            qout(i, 1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! North Edges:
        IF (je + 1 .EQ. npy) THEN
          DO i=is1,ie1
            q1(i) = 0.5*(qin(i, npy-1)+qin(i, npy))
          END DO
          DO i=is2,ie1
            CALL PUSHREALARRAY(qout(i, npy))
            qout(i, npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
          END DO
          CALL PUSHCONTROL(2,1)
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY(qout(i, j))
          qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin(i&
&           , j))
        END DO
      END DO
      CALL PUSHCONTROL(2,3)
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            CALL PUSHREALARRAY(qin(i, j))
            qin(i, j) = qout(i, j)
          END DO
        END DO
        CALL PUSHINTEGER(ie1)
        CALL PUSHINTEGER(js2)
        CALL PUSHINTEGER(js1)
        CALL PUSHINTEGER(je1)
        CALL PUSHINTEGER(is2)
        CALL PUSHINTEGER(is1)
        CALL PUSHCONTROL(2,2)
      ELSE
        CALL PUSHINTEGER(ie1)
        CALL PUSHINTEGER(js2)
        CALL PUSHINTEGER(js1)
        CALL PUSHINTEGER(je1)
        CALL PUSHINTEGER(is2)
        CALL PUSHINTEGER(is1)
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHINTEGER(ie1)
      CALL PUSHINTEGER(js2)
      CALL PUSHINTEGER(js1)
      CALL PUSHINTEGER(je1)
      CALL PUSHINTEGER(is2)
      CALL PUSHINTEGER(is1)
      CALL PUSHCONTROL(2,0)
    END IF
  END SUBROUTINE A2B_ORD2_FWD
!  Differentiation of a2b_ord2 in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
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
!   gradient     of useful results: qin qout
!   with respect to varying inputs: qin qout
  SUBROUTINE A2B_ORD2_BWD(qin, qin_ad, qout, qout_ad, gridstruct, npx, &
&   npy, is, ie, js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
    REAL, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qin_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qout_ad(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
    REAL :: q1(npx), q2(npy)
    REAL :: q1_ad(npx), q2_ad(npy)
    INTEGER :: i, j
    INTEGER :: is1, js1, is2, js2, ie1, je1
    REAL, DIMENSION(:, :, :), POINTER :: grid, agrid
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
    REAL(kind=r_grid), DIMENSION(:), POINTER :: edge_w, edge_e, edge_s, &
&   edge_n
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    INTEGER :: branch

    q1 = 0.0
    q2 = 0.0
    i = 0
    j = 0
    is1 = 0
    js1 = 0
    is2 = 0
    js2 = 0
    ie1 = 0
    je1 = 0
    branch = 0

    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(is1)
      CALL POPINTEGER(is2)
      CALL POPINTEGER(je1)
      CALL POPINTEGER(js1)
      CALL POPINTEGER(js2)
      CALL POPINTEGER(ie1)
    ELSE IF (branch .EQ. 1) THEN
      CALL POPINTEGER(is1)
      CALL POPINTEGER(is2)
      CALL POPINTEGER(je1)
      CALL POPINTEGER(js1)
      CALL POPINTEGER(js2)
      CALL POPINTEGER(ie1)
    ELSE
      CALL POPINTEGER(is1)
      CALL POPINTEGER(is2)
      CALL POPINTEGER(je1)
      CALL POPINTEGER(js1)
      CALL POPINTEGER(js2)
      CALL POPINTEGER(ie1)
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(qin(i, j))
          qout_ad(i, j) = qout_ad(i, j) + qin_ad(i, j)
          qin_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    edge_n => gridstruct%edge_n
    CALL POPCONTROL(2,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) THEN
        DO j=je+3,js-2,-1
          DO i=ie+3,is-2,-1
            CALL POPREALARRAY(qout(i, j))
            temp_ad = 0.25*qout_ad(i, j)
            qin_ad(i-1, j-1) = qin_ad(i-1, j-1) + temp_ad
            qin_ad(i, j-1) = qin_ad(i, j-1) + temp_ad
            qin_ad(i-1, j) = qin_ad(i-1, j) + temp_ad
            qin_ad(i, j) = qin_ad(i, j) + temp_ad
            qout_ad(i, j) = 0.0
          END DO
        END DO
        GOTO 100
      ELSE
        q1_ad = 0.0
        DO i=ie1,is2,-1
          CALL POPREALARRAY(qout(i, npy))
          q1_ad(i-1) = q1_ad(i-1) + edge_n(i)*qout_ad(i, npy)
          q1_ad(i) = q1_ad(i) + (1.-edge_n(i))*qout_ad(i, npy)
          qout_ad(i, npy) = 0.0
        END DO
        DO i=ie1,is1,-1
          qin_ad(i, npy-1) = qin_ad(i, npy-1) + 0.5*q1_ad(i)
          qin_ad(i, npy) = qin_ad(i, npy) + 0.5*q1_ad(i)
          q1_ad(i) = 0.0
        END DO
      END IF
    ELSE IF (branch .EQ. 2) THEN
      q1_ad = 0.0
    ELSE
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(qout(i, j))
          temp_ad5 = 0.25*qout_ad(i, j)
          qin_ad(i-1, j-1) = qin_ad(i-1, j-1) + temp_ad5
          qin_ad(i, j-1) = qin_ad(i, j-1) + temp_ad5
          qin_ad(i-1, j) = qin_ad(i-1, j) + temp_ad5
          qin_ad(i, j) = qin_ad(i, j) + temp_ad5
          qout_ad(i, j) = 0.0
        END DO
      END DO
      GOTO 100
    END IF
    edge_s => gridstruct%edge_s
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO i=ie1,is2,-1
        CALL POPREALARRAY(qout(i, 1))
        q1_ad(i-1) = q1_ad(i-1) + edge_s(i)*qout_ad(i, 1)
        q1_ad(i) = q1_ad(i) + (1.-edge_s(i))*qout_ad(i, 1)
        qout_ad(i, 1) = 0.0
      END DO
      DO i=ie1,is1,-1
        qin_ad(i, 0) = qin_ad(i, 0) + 0.5*q1_ad(i)
        qin_ad(i, 1) = qin_ad(i, 1) + 0.5*q1_ad(i)
        q1_ad(i) = 0.0
      END DO
    END IF
    edge_e => gridstruct%edge_e
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      q2_ad = 0.0
      DO j=je1,js2,-1
        CALL POPREALARRAY(qout(npx, j))
        q2_ad(j-1) = q2_ad(j-1) + edge_e(j)*qout_ad(npx, j)
        q2_ad(j) = q2_ad(j) + (1.-edge_e(j))*qout_ad(npx, j)
        qout_ad(npx, j) = 0.0
      END DO
      DO j=je1,js1,-1
        qin_ad(npx-1, j) = qin_ad(npx-1, j) + 0.5*q2_ad(j)
        qin_ad(npx, j) = qin_ad(npx, j) + 0.5*q2_ad(j)
        q2_ad(j) = 0.0
      END DO
    ELSE
      q2_ad = 0.0
    END IF
    edge_w => gridstruct%edge_w
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=je1,js2,-1
        CALL POPREALARRAY(qout(1, j))
        q2_ad(j-1) = q2_ad(j-1) + edge_w(j)*qout_ad(1, j)
        q2_ad(j) = q2_ad(j) + (1.-edge_w(j))*qout_ad(1, j)
        qout_ad(1, j) = 0.0
      END DO
      DO j=je1,js1,-1
        qin_ad(0, j) = qin_ad(0, j) + 0.5*q2_ad(j)
        qin_ad(1, j) = qin_ad(1, j) + 0.5*q2_ad(j)
        q2_ad(j) = 0.0
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(qout(1, npy))
      temp_ad4 = r3*qout_ad(1, npy)
      qin_ad(1, npy-1) = qin_ad(1, npy-1) + temp_ad4
      qin_ad(0, npy-1) = qin_ad(0, npy-1) + temp_ad4
      qin_ad(1, npy) = qin_ad(1, npy) + temp_ad4
      qout_ad(1, npy) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(qout(npx, npy))
      temp_ad3 = r3*qout_ad(npx, npy)
      qin_ad(npx-1, npy-1) = qin_ad(npx-1, npy-1) + temp_ad3
      qin_ad(npx, npy-1) = qin_ad(npx, npy-1) + temp_ad3
      qin_ad(npx-1, npy) = qin_ad(npx-1, npy) + temp_ad3
      qout_ad(npx, npy) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(qout(npx, 1))
      temp_ad2 = r3*qout_ad(npx, 1)
      qin_ad(npx-1, 1) = qin_ad(npx-1, 1) + temp_ad2
      qin_ad(npx-1, 0) = qin_ad(npx-1, 0) + temp_ad2
      qin_ad(npx, 1) = qin_ad(npx, 1) + temp_ad2
      qout_ad(npx, 1) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(qout(1, 1))
      temp_ad1 = r3*qout_ad(1, 1)
      qin_ad(1, 1) = qin_ad(1, 1) + temp_ad1
      qin_ad(1, 0) = qin_ad(1, 0) + temp_ad1
      qin_ad(0, 1) = qin_ad(0, 1) + temp_ad1
      qout_ad(1, 1) = 0.0
    END IF
    DO j=je1,js2,-1
      DO i=ie1,is2,-1
        CALL POPREALARRAY(qout(i, j))
        temp_ad0 = 0.25*qout_ad(i, j)
        qin_ad(i-1, j-1) = qin_ad(i-1, j-1) + temp_ad0
        qin_ad(i, j-1) = qin_ad(i, j-1) + temp_ad0
        qin_ad(i-1, j) = qin_ad(i-1, j) + temp_ad0
        qin_ad(i, j) = qin_ad(i, j) + temp_ad0
        qout_ad(i, j) = 0.0
      END DO
    END DO
 100 CONTINUE
  END SUBROUTINE A2B_ORD2_BWD
  SUBROUTINE A2B_ORD2(qin, qout, gridstruct, npx, npy, is, ie, js, je, &
&   ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(OUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local:
    REAL :: q1(npx), q2(npy)
    INTEGER :: i, j
    INTEGER :: is1, js1, is2, js2, ie1, je1
    REAL, DIMENSION(:, :, :), POINTER :: grid, agrid
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
    REAL(kind=r_grid), DIMENSION(:), POINTER :: edge_w, edge_e, edge_s, &
&   edge_n
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    edge_w => gridstruct%edge_w
    edge_e => gridstruct%edge_e
    edge_s => gridstruct%edge_s
    edge_n => gridstruct%edge_n
    grid => gridstruct%grid
    agrid => gridstruct%agrid
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    IF (gridstruct%grid_type .LT. 3) THEN
      IF (gridstruct%nested) THEN
        DO j=js-2,je+1+2
          DO i=is-2,ie+1+2
            qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin&
&             (i, j))
          END DO
        END DO
      ELSE
        IF (1 .LT. is - 1) THEN
          is1 = is - 1
        ELSE
          is1 = 1
        END IF
        IF (1 .LT. js - 1) THEN
          js1 = js - 1
        ELSE
          js1 = 1
        END IF
        IF (2 .LT. is) THEN
          is2 = is
        ELSE
          is2 = 2
        END IF
        IF (2 .LT. js) THEN
          js2 = js
        ELSE
          js2 = 2
        END IF
        IF (npx - 1 .GT. ie + 1) THEN
          ie1 = ie + 1
        ELSE
          ie1 = npx - 1
        END IF
        IF (npy - 1 .GT. je + 1) THEN
          je1 = je + 1
        ELSE
          je1 = npy - 1
        END IF
        DO j=js2,je1
          DO i=is2,ie1
            qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin&
&             (i, j))
          END DO
        END DO
! Fix the 4 Corners:
        IF (gridstruct%sw_corner) qout(1, 1) = r3*(qin(1, 1)+qin(1, 0)+&
&           qin(0, 1))
        IF (gridstruct%se_corner) qout(npx, 1) = r3*(qin(npx-1, 1)+qin(&
&           npx-1, 0)+qin(npx, 1))
        IF (gridstruct%ne_corner) qout(npx, npy) = r3*(qin(npx-1, npy-1)&
&           +qin(npx, npy-1)+qin(npx-1, npy))
        IF (gridstruct%nw_corner) qout(1, npy) = r3*(qin(1, npy-1)+qin(0&
&           , npy-1)+qin(1, npy))
! *** West Edges:
        IF (is .EQ. 1) THEN
          DO j=js1,je1
            q2(j) = 0.5*(qin(0, j)+qin(1, j))
          END DO
          DO j=js2,je1
            qout(1, j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
          END DO
        END IF
! East Edges:
        IF (ie + 1 .EQ. npx) THEN
          DO j=js1,je1
            q2(j) = 0.5*(qin(npx-1, j)+qin(npx, j))
          END DO
          DO j=js2,je1
            qout(npx, j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
          END DO
        END IF
! South Edges:
        IF (js .EQ. 1) THEN
          DO i=is1,ie1
            q1(i) = 0.5*(qin(i, 0)+qin(i, 1))
          END DO
          DO i=is2,ie1
            qout(i, 1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
          END DO
        END IF
! North Edges:
        IF (je + 1 .EQ. npy) THEN
          DO i=is1,ie1
            q1(i) = 0.5*(qin(i, npy-1)+qin(i, npy))
          END DO
          DO i=is2,ie1
            qout(i, npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
          END DO
        END IF
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin(i&
&           , j))
        END DO
      END DO
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qin(i, j) = qout(i, j)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE A2B_ORD2
!  Differentiation of extrap_corner in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 
!dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c
! dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_co
!re_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dy
!namics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_ut
!ils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_map
!z_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.s
!calar_profile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.
!steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d
!2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod
!.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Sol
!ver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_sol
!ver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_n
!h sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_
!core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_m
!od.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_cor
!ners_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_c
!ircle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: extrap_corner q1 q2
!   with respect to varying inputs: q1 q2
  SUBROUTINE EXTRAP_CORNER_ADM(p0, p1, p2, q1, q1_ad, q2, q2_ad, &
&   extrap_corner_ad)
    IMPLICIT NONE
    REAL, DIMENSION(2), INTENT(IN) :: p0, p1, p2
    REAL, INTENT(IN) :: q1, q2
    REAL :: q1_ad, q2_ad
    REAL :: x1, x2
    INTRINSIC REAL
    REAL(r_grid), DIMENSION(2) :: arg1
    REAL(r_grid), DIMENSION(2) :: arg2
    REAL :: temp_ad
    REAL :: extrap_corner_ad
    REAL :: extrap_corner
    arg1(:) = REAL(p1, kind=r_grid)
    arg2(:) = REAL(p0, kind=r_grid)
    x1 = GREAT_CIRCLE_DIST(REAL(p1, kind=r_grid), REAL(p0, kind=r_grid))
    arg1(:) = REAL(p2, kind=r_grid)
    arg2(:) = REAL(p0, kind=r_grid)
    x2 = GREAT_CIRCLE_DIST(REAL(p2, kind=r_grid), REAL(p0, kind=r_grid))
    temp_ad = x1*extrap_corner_ad/(x2-x1)
    q1_ad = q1_ad + temp_ad + extrap_corner_ad
    q2_ad = q2_ad - temp_ad
  END SUBROUTINE EXTRAP_CORNER_ADM
  REAL FUNCTION EXTRAP_CORNER(p0, p1, p2, q1, q2)
    IMPLICIT NONE
    REAL, DIMENSION(2), INTENT(IN) :: p0, p1, p2
    REAL, INTENT(IN) :: q1, q2
    REAL :: x1, x2
    INTRINSIC REAL
    REAL(r_grid), DIMENSION(2) :: arg1
    REAL(r_grid), DIMENSION(2) :: arg2
    arg1(:) = REAL(p1, kind=r_grid)
    arg2(:) = REAL(p0, kind=r_grid)
    x1 = GREAT_CIRCLE_DIST(arg1(:), arg2(:))
    arg1(:) = REAL(p2, kind=r_grid)
    arg2(:) = REAL(p0, kind=r_grid)
    x2 = GREAT_CIRCLE_DIST(arg1(:), arg2(:))
    extrap_corner = q1 + x1/(x2-x1)*(q1-q2)
  END FUNCTION EXTRAP_CORNER
!  Differentiation of extrap_corner in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_
!edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn
!_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core
!_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Ra
!yleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c
!2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_
!mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.rema
!p_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_lim
!iters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubi
!c fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_
!subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_
!utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_ut
!ils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_uti
!ls_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_
!mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.
!ytp_v sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp
!_core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid
!_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: q1 q2 extrap_corner
!   with respect to varying inputs: q1 q2
  REAL FUNCTION EXTRAP_CORNER_FWD(p0, p1, p2, q1, q2)
    IMPLICIT NONE
    REAL, DIMENSION(2), INTENT(IN) :: p0, p1, p2
    REAL, INTENT(IN) :: q1, q2
    REAL :: x1, x2
    INTRINSIC REAL
    REAL(r_grid), DIMENSION(2) :: arg1
    REAL(r_grid), DIMENSION(2) :: arg2
    REAL :: extrap_corner
    arg1(:) = REAL(p1, kind=r_grid)
    arg2(:) = REAL(p0, kind=r_grid)
    x1 = GREAT_CIRCLE_DIST(REAL(p1, kind=r_grid), REAL(p0, kind=r_grid))
    arg1(:) = REAL(p2, kind=r_grid)
    arg2(:) = REAL(p0, kind=r_grid)
    x2 = GREAT_CIRCLE_DIST(REAL(p2, kind=r_grid), REAL(p0, kind=r_grid))
    extrap_corner = q1 + x1/(x2-x1)*(q1-q2)
    CALL PUSHREALARRAY(x2)
    CALL PUSHREALARRAY(x1)
    extrap_corner_fwd = extrap_corner
  END FUNCTION EXTRAP_CORNER_FWD
!  Differentiation of extrap_corner in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b
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
!   gradient     of useful results: q1 q2 extrap_corner
!   with respect to varying inputs: q1 q2
  SUBROUTINE EXTRAP_CORNER_BWD(p0, p1, p2, q1, q1_ad, q2, q2_ad, &
&   extrap_corner_ad)
    IMPLICIT NONE
    REAL, DIMENSION(2), INTENT(IN) :: p0, p1, p2
    REAL, INTENT(IN) :: q1, q2
    REAL :: q1_ad, q2_ad
    REAL :: x1, x2
    INTRINSIC REAL
    REAL(r_grid), DIMENSION(2) :: arg1
    REAL(r_grid), DIMENSION(2) :: arg2
    REAL :: temp_ad
    REAL :: extrap_corner
    REAL :: extrap_corner_ad
    CALL POPREALARRAY(x1)
    CALL POPREALARRAY(x2)
    temp_ad = x1*extrap_corner_ad/(x2-x1)
    q1_ad = q1_ad + temp_ad + extrap_corner_ad
    q2_ad = q2_ad - temp_ad
  END SUBROUTINE EXTRAP_CORNER_BWD
end module a2b_edge_adm_mod
