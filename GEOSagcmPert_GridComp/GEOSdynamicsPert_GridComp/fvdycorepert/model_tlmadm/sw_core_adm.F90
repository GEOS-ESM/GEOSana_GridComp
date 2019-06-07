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
 module sw_core_adm_mod

 use fv_mp_mod,        only: ng
 use tp_core_adm_mod,  only: fv_tp_2d, pert_ppm, copy_corners
 use tp_core_adm_mod,  only: fv_tp_2d_adm, fv_tp_2d_fwd, fv_tp_2d_bwd, copy_corners_adm, pert_ppm_adm
 use fv_mp_mod,        only: XDir, YDir
 use fv_mp_mod,        only: fill_corners
 use fv_mp_adm_mod,    only: fill_corners_adm
 use fv_arrays_mod,    only: fv_grid_type, fv_grid_bounds_type, fv_flags_type, r_grid
 use a2b_edge_adm_mod, only: a2b_ord4
 use a2b_edge_adm_mod, only: a2b_ord4_fwd, a2b_ord4_bwd, a2b_ord4_adm
 use fv_arrays_mod,    only: fvprc

 use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                          pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

#ifdef SW_DYNAMICS
 use test_cases_mod,   only: test_case
#endif

 use fv_arrays_nlm_mod, only: fpp

 implicit none

  real, parameter:: r3 =   1./3.
  real, parameter:: t11=27./28., t12=-13./28., t13=3./7., t14=6./7., t15=3./28.
  real, parameter:: s11=11./14., s13=-13./14., s14=4./7., s15=3./14.
  real, parameter:: near_zero = 1.E-9     ! for KE limiter
#ifdef OVERLOAD_R4
  real, parameter:: big_number = 1.E8
#else
  real, parameter:: big_number = 1.E30
#endif
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.
!----------------------------
! 4-pt Lagrange interpolation
!----------------------------
  real, parameter:: a1 =  0.5625
  real, parameter:: a2 = -0.0625
!----------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
! 3-pt off-center intp formular:
! real, parameter:: c1 = -0.125
! real, parameter:: c2 =  0.75
! real, parameter:: c3 =  0.375
!----------------------------------------------
! scheme 2.1: perturbation form
  REAL, PARAMETER :: b1=1./30.
  REAL, PARAMETER :: b2=-(13./60.)
  REAL, PARAMETER :: b3=-(13./60.)
  REAL, PARAMETER :: b4=0.45
  REAL, PARAMETER :: b5=-0.05
  PRIVATE 
  PUBLIC c_sw, d_sw, fill_4corners, del6_vt_flux, divergence_corner, &
& divergence_corner_nest
  PUBLIC c_sw_fwd, c_sw_bwd, d_sw_fwd, d_sw_bwd, fill_4corners_fwd, &
& fill_4corners_bwd, del6_vt_flux_adm, divergence_corner_fwd, &
& divergence_corner_bwd, divergence_corner_nest_fwd, &
& divergence_corner_nest_bwd
  PUBLIC d2a2c_vect
  PUBLIC d2a2c_vect_fwd, d2a2c_vect_bwd
  EXTERNAL D_SW_ADM
  EXTERNAL DIVERGENCE_CORNER_NEST_ADM
  EXTERNAL DIVERGENCE_CORNER_ADM
  EXTERNAL C_SW_ADM
  EXTERNAL FILL_4CORNERS_ADM

CONTAINS
!  Differentiation of c_sw in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b
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
!   gradient     of useful results: u v w delp ua uc ptc ut delpc
!                va vc vt divg_d wc pt
!   with respect to varying inputs: u v w delp ua uc ptc ut delpc
!                va vc vt divg_d wc pt
  SUBROUTINE C_SW_FWD(delpc, delp, ptc, pt, u, v, w, uc, vc, ua, va, wc&
&   , ut, vt, divg_d, nord, dt2, hydrostatic, dord4, bd, gridstruct, &
&   flagstruct)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: u&
&   , vc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: v&
&   , uc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: delp&
&   , pt, ua, va, ut, vt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: delpc, ptc, wc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1) :: divg_d
    INTEGER, INTENT(IN) :: nord
    REAL, INTENT(IN) :: dt2
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: dord4
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! Local:
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+1) :: vort, ke
    REAL, DIMENSION(bd%is-1:bd%ie+2, bd%js-1:bd%je+1) :: fx, fx1, fx2
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+2) :: fy, fy1, fy2
    REAL :: dt4
    INTEGER :: i, j, is2, ie1
    INTEGER :: iep1, jep1
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: npx, npy
    LOGICAL :: nested
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg, cos_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: dx, dy, dxc, dyc

    vort = 0.0
    ke = 0.0
    fx = 0.0
    fx1 = 0.0
    fx2 = 0.0
    fy = 0.0
    fy1 = 0.0
    fy2 = 0.0
    dt4 = 0.0
    is2 = 0
    ie1 = 0
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
    npx = 0
    npy = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    sin_sg => gridstruct%sin_sg
    cos_sg => gridstruct%cos_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    sina_u => gridstruct%sina_u
    sina_v => gridstruct%sina_v
    dx => gridstruct%dx
    dy => gridstruct%dy
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    sw_corner = gridstruct%sw_corner
    se_corner = gridstruct%se_corner
    nw_corner = gridstruct%nw_corner
    ne_corner = gridstruct%ne_corner
    iep1 = ie + 1
    jep1 = je + 1
    CALL D2A2C_VECT_FWD(u, v, ua, va, uc, vc, ut, vt, dord4, gridstruct&
&                 , bd, npx, npy, nested, flagstruct%grid_type)
    IF (nord .GT. 0) THEN
      IF (nested) THEN
        CALL DIVERGENCE_CORNER_NEST_FWD(u, v, ua, va, divg_d, gridstruct&
&                                 , flagstruct, bd)
        CALL PUSHCONTROL(2,2)
      ELSE
        CALL DIVERGENCE_CORNER_FWD(u, v, ua, va, divg_d, gridstruct, &
&                            flagstruct, bd)
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHCONTROL(2,0)
    END IF
    DO j=js-1,jep1
      DO i=is-1,iep1+1
        IF (ut(i, j) .GT. 0.) THEN
          CALL PUSHREALARRAY(ut(i, j))
          ut(i, j) = dt2*ut(i, j)*dy(i, j)*sin_sg(i-1, j, 3)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHREALARRAY(ut(i, j))
          ut(i, j) = dt2*ut(i, j)*dy(i, j)*sin_sg(i, j, 1)
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
    END DO
    DO j=js-1,je+2
      DO i=is-1,iep1
        IF (vt(i, j) .GT. 0.) THEN
          CALL PUSHREALARRAY(vt(i, j))
          vt(i, j) = dt2*vt(i, j)*dx(i, j)*sin_sg(i, j-1, 4)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHREALARRAY(vt(i, j))
          vt(i, j) = dt2*vt(i, j)*dx(i, j)*sin_sg(i, j, 2)
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
    END DO
!----------------
! Transport delp:
!----------------
! Xdir:
    IF (flagstruct%grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      CALL FILL2_4CORNERS_FWD(delp, pt, 1, bd, npx, npy, sw_corner, &
&                       se_corner, ne_corner, nw_corner)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (hydrostatic) THEN
      DO j=js-1,jep1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1(i, j) = delp(i-1, j)
            fx(i, j) = pt(i-1, j)
            CALL PUSHCONTROL(1,0)
          ELSE
            fx1(i, j) = delp(i, j)
            fx(i, j) = pt(i, j)
            CALL PUSHCONTROL(1,1)
          END IF
          CALL PUSHREALARRAY(fx1(i, j))
          fx1(i, j) = ut(i, j)*fx1(i, j)
          CALL PUSHREALARRAY(fx(i, j))
          fx(i, j) = fx1(i, j)*fx(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      IF (flagstruct%grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_FWD(w, 1, bd, npx, npy, sw_corner, se_corner&
&                        , ne_corner, nw_corner)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1(i, j) = delp(i-1, j)
            fx(i, j) = pt(i-1, j)
            fx2(i, j) = w(i-1, j)
            CALL PUSHCONTROL(1,0)
          ELSE
            fx1(i, j) = delp(i, j)
            fx(i, j) = pt(i, j)
            fx2(i, j) = w(i, j)
            CALL PUSHCONTROL(1,1)
          END IF
          CALL PUSHREALARRAY(fx1(i, j))
          fx1(i, j) = ut(i, j)*fx1(i, j)
          CALL PUSHREALARRAY(fx(i, j))
          fx(i, j) = fx1(i, j)*fx(i, j)
          CALL PUSHREALARRAY(fx2(i, j))
          fx2(i, j) = fx1(i, j)*fx2(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    END IF
! Ydir:
    IF (flagstruct%grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      CALL FILL2_4CORNERS_FWD(delp, pt, 2, bd, npx, npy, sw_corner, &
&                       se_corner, ne_corner, nw_corner)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (hydrostatic) THEN
      DO j=js-1,jep1+1
        DO i=is-1,iep1
          IF (vt(i, j) .GT. 0.) THEN
            fy1(i, j) = delp(i, j-1)
            fy(i, j) = pt(i, j-1)
            CALL PUSHCONTROL(1,0)
          ELSE
            fy1(i, j) = delp(i, j)
            fy(i, j) = pt(i, j)
            CALL PUSHCONTROL(1,1)
          END IF
          CALL PUSHREALARRAY(fy1(i, j))
          fy1(i, j) = vt(i, j)*fy1(i, j)
          CALL PUSHREALARRAY(fy(i, j))
          fy(i, j) = fy1(i, j)*fy(i, j)
        END DO
      END DO
      DO j=js-1,jep1
        DO i=is-1,iep1
          CALL PUSHREALARRAY(delpc(i, j))
          delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+(fy1(i, j)-&
&           fy1(i, j+1)))*gridstruct%rarea(i, j)
          CALL PUSHREALARRAY(ptc(i, j))
          ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+(fy(i, j&
&           )-fy(i, j+1)))*gridstruct%rarea(i, j))/delpc(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      IF (flagstruct%grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_FWD(w, 2, bd, npx, npy, sw_corner, se_corner&
&                        , ne_corner, nw_corner)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
      DO j=js-1,je+2
        DO i=is-1,ie+1
          IF (vt(i, j) .GT. 0.) THEN
            fy1(i, j) = delp(i, j-1)
            fy(i, j) = pt(i, j-1)
            fy2(i, j) = w(i, j-1)
            CALL PUSHCONTROL(1,0)
          ELSE
            fy1(i, j) = delp(i, j)
            fy(i, j) = pt(i, j)
            fy2(i, j) = w(i, j)
            CALL PUSHCONTROL(1,1)
          END IF
          CALL PUSHREALARRAY(fy1(i, j))
          fy1(i, j) = vt(i, j)*fy1(i, j)
          CALL PUSHREALARRAY(fy(i, j))
          fy(i, j) = fy1(i, j)*fy(i, j)
          CALL PUSHREALARRAY(fy2(i, j))
          fy2(i, j) = fy1(i, j)*fy2(i, j)
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          CALL PUSHREALARRAY(delpc(i, j))
          delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+(fy1(i, j)-&
&           fy1(i, j+1)))*gridstruct%rarea(i, j)
          CALL PUSHREALARRAY(ptc(i, j))
          ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+(fy(i, j&
&           )-fy(i, j+1)))*gridstruct%rarea(i, j))/delpc(i, j)
          CALL PUSHREALARRAY(wc(i, j))
          wc(i, j) = (w(i, j)*delp(i, j)+(fx2(i, j)-fx2(i+1, j)+(fy2(i, &
&           j)-fy2(i, j+1)))*gridstruct%rarea(i, j))/delpc(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    END IF
!------------
! Compute KE:
!------------
!Since uc = u*, i.e. the covariant wind perpendicular to the face edge, if we want to compute kinetic energy we will need the tru
!e coordinate-parallel covariant wind, computed through u = uc*sina + v*cosa.
!Use the alpha for the cell KE is being computed in.
!!! TO DO:
!!! Need separate versions for nesting/single-tile
!!!   and for cubed-sphere
    IF (nested .OR. flagstruct%grid_type .GE. 3) THEN
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (ua(i, j) .GT. 0.) THEN
            ke(i, j) = uc(i, j)
            CALL PUSHCONTROL(1,1)
          ELSE
            ke(i, j) = uc(i+1, j)
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (va(i, j) .GT. 0.) THEN
            vort(i, j) = vc(i, j)
            CALL PUSHCONTROL(1,1)
          ELSE
            vort(i, j) = vc(i, j+1)
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (ua(i, j) .GT. 0.) THEN
            IF (i .EQ. 1) THEN
              ke(1, j) = uc(1, j)*sin_sg(1, j, 1) + v(1, j)*cos_sg(1, j&
&               , 1)
              CALL PUSHCONTROL(3,5)
            ELSE IF (i .EQ. npx) THEN
              ke(i, j) = uc(npx, j)*sin_sg(npx, j, 1) + v(npx, j)*cos_sg&
&               (npx, j, 1)
              CALL PUSHCONTROL(3,4)
            ELSE
              ke(i, j) = uc(i, j)
              CALL PUSHCONTROL(3,3)
            END IF
          ELSE IF (i .EQ. 0) THEN
            ke(0, j) = uc(1, j)*sin_sg(0, j, 3) + v(1, j)*cos_sg(0, j, 3&
&             )
            CALL PUSHCONTROL(3,2)
          ELSE IF (i .EQ. npx - 1) THEN
            ke(i, j) = uc(npx, j)*sin_sg(npx-1, j, 3) + v(npx, j)*cos_sg&
&             (npx-1, j, 3)
            CALL PUSHCONTROL(3,1)
          ELSE
            ke(i, j) = uc(i+1, j)
            CALL PUSHCONTROL(3,0)
          END IF
        END DO
      END DO
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (va(i, j) .GT. 0.) THEN
            IF (j .EQ. 1) THEN
              vort(i, 1) = vc(i, 1)*sin_sg(i, 1, 2) + u(i, 1)*cos_sg(i, &
&               1, 2)
              CALL PUSHCONTROL(3,5)
            ELSE IF (j .EQ. npy) THEN
              vort(i, j) = vc(i, npy)*sin_sg(i, npy, 2) + u(i, npy)*&
&               cos_sg(i, npy, 2)
              CALL PUSHCONTROL(3,4)
            ELSE
              vort(i, j) = vc(i, j)
              CALL PUSHCONTROL(3,3)
            END IF
          ELSE IF (j .EQ. 0) THEN
            vort(i, 0) = vc(i, 1)*sin_sg(i, 0, 4) + u(i, 1)*cos_sg(i, 0&
&             , 4)
            CALL PUSHCONTROL(3,2)
          ELSE IF (j .EQ. npy - 1) THEN
            vort(i, j) = vc(i, npy)*sin_sg(i, npy-1, 4) + u(i, npy)*&
&             cos_sg(i, npy-1, 4)
            CALL PUSHCONTROL(3,1)
          ELSE
            vort(i, j) = vc(i, j+1)
            CALL PUSHCONTROL(3,0)
          END IF
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    END IF
    dt4 = 0.5*dt2
    DO j=js-1,jep1
      DO i=is-1,iep1
        CALL PUSHREALARRAY(ke(i, j))
        ke(i, j) = dt4*(ua(i, j)*ke(i, j)+va(i, j)*vort(i, j))
      END DO
    END DO
!------------------------------
! Compute circulation on C grid
!------------------------------
! To consider using true co-variant winds at face edges?
    DO j=js-1,je+1
      DO i=is,ie+1
        CALL PUSHREALARRAY(fx(i, j))
        fx(i, j) = uc(i, j)*dxc(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is-1,ie+1
        CALL PUSHREALARRAY(fy(i, j))
        fy(i, j) = vc(i, j)*dyc(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie+1
        CALL PUSHREALARRAY(vort(i, j))
        vort(i, j) = fx(i, j-1) - fx(i, j) + (fy(i, j)-fy(i-1, j))
      END DO
    END DO
! Remove the extra term at the corners:
    IF (sw_corner) THEN
      CALL PUSHREALARRAY(vort(1, 1))
      vort(1, 1) = vort(1, 1) + fy(0, 1)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (se_corner) THEN
      CALL PUSHREALARRAY(vort(npx, 1))
      vort(npx, 1) = vort(npx, 1) - fy(npx, 1)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (ne_corner) THEN
      CALL PUSHREALARRAY(vort(npx, npy))
      vort(npx, npy) = vort(npx, npy) - fy(npx, npy)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (nw_corner) THEN
      CALL PUSHREALARRAY(vort(1, npy))
      vort(1, npy) = vort(1, npy) + fy(0, npy)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
!----------------------------
! Compute absolute vorticity
!----------------------------
    DO j=js,je+1
      DO i=is,ie+1
        CALL PUSHREALARRAY(vort(i, j))
        vort(i, j) = gridstruct%fc(i, j) + gridstruct%rarea_c(i, j)*vort&
&         (i, j)
      END DO
    END DO
!----------------------------------
! Transport absolute vorticity:
!----------------------------------
!To go from v to contravariant v at the edges, we divide by sin_sg;
! but we then must multiply by sin_sg to get the proper flux.
! These cancel, leaving us with fy1 = dt2*v at the edges.
! (For the same reason we only divide by sin instead of sin**2 in the interior)
!! TO DO: separate versions for nesting/single-tile and cubed-sphere
    IF (nested .OR. flagstruct%grid_type .GE. 3) THEN
      DO j=js,je
        DO i=is,iep1
          CALL PUSHREALARRAY(fy1(i, j))
          fy1(i, j) = dt2*(v(i, j)-uc(i, j)*cosa_u(i, j))/sina_u(i, j)
          IF (fy1(i, j) .GT. 0.) THEN
            CALL PUSHREALARRAY(fy(i, j))
            fy(i, j) = vort(i, j)
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHREALARRAY(fy(i, j))
            fy(i, j) = vort(i, j+1)
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      DO j=js,jep1
        DO i=is,ie
          CALL PUSHREALARRAY(fx1(i, j))
          fx1(i, j) = dt2*(u(i, j)-vc(i, j)*cosa_v(i, j))/sina_v(i, j)
          IF (fx1(i, j) .GT. 0.) THEN
            CALL PUSHREALARRAY(fx(i, j))
            fx(i, j) = vort(i, j)
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHREALARRAY(fx(i, j))
            fx(i, j) = vort(i+1, j)
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    ELSE
      DO j=js,je
!DEC$ VECTOR ALWAYS
        DO i=is,iep1
          IF (i .EQ. 1 .OR. i .EQ. npx) THEN
            CALL PUSHREALARRAY(fy1(i, j))
            fy1(i, j) = dt2*v(i, j)
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHREALARRAY(fy1(i, j))
            fy1(i, j) = dt2*(v(i, j)-uc(i, j)*cosa_u(i, j))/sina_u(i, j)
            CALL PUSHCONTROL(1,1)
          END IF
          IF (fy1(i, j) .GT. 0.) THEN
            CALL PUSHREALARRAY(fy(i, j))
            fy(i, j) = vort(i, j)
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHREALARRAY(fy(i, j))
            fy(i, j) = vort(i, j+1)
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      DO j=js,jep1
        IF (j .EQ. 1 .OR. j .EQ. npy) THEN
!DEC$ VECTOR ALWAYS
          DO i=is,ie
            CALL PUSHREALARRAY(fx1(i, j))
            fx1(i, j) = dt2*u(i, j)
            IF (fx1(i, j) .GT. 0.) THEN
              CALL PUSHREALARRAY(fx(i, j))
              fx(i, j) = vort(i, j)
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHREALARRAY(fx(i, j))
              fx(i, j) = vort(i+1, j)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
!DEC$ VECTOR ALWAYS
          DO i=is,ie
            CALL PUSHREALARRAY(fx1(i, j))
            fx1(i, j) = dt2*(u(i, j)-vc(i, j)*cosa_v(i, j))/sina_v(i, j)
            IF (fx1(i, j) .GT. 0.) THEN
              CALL PUSHREALARRAY(fx(i, j))
              fx(i, j) = vort(i, j)
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHREALARRAY(fx(i, j))
              fx(i, j) = vort(i+1, j)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHCONTROL(1,0)
    END IF
! Update time-centered winds on the C-Grid
    DO j=js,je
      DO i=is,iep1
        uc(i, j) = uc(i, j) + fy1(i, j)*fy(i, j) + gridstruct%rdxc(i, j)&
&         *(ke(i-1, j)-ke(i, j))
      END DO
    END DO
    DO j=js,jep1
      DO i=is,ie
        vc(i, j) = vc(i, j) - fx1(i, j)*fx(i, j) + gridstruct%rdyc(i, j)&
&         *(ke(i, j-1)-ke(i, j))
      END DO
    END DO
    CALL PUSHREALARRAY(fx2, (bd%ie-bd%is+4)*(bd%je-bd%js+3))
    CALL PUSHREALARRAY(fx1, (bd%ie-bd%is+4)*(bd%je-bd%js+3))
    CALL PUSHINTEGER(je)
    CALL PUSHREALARRAY(fy, (bd%ie-bd%is+3)*(bd%je-bd%js+4))
    CALL PUSHREALARRAY(fx, (bd%ie-bd%is+4)*(bd%je-bd%js+3))
    CALL PUSHINTEGER(is)
    !CALL PUSHPOINTER8(C_LOC(sina_v))
    !CALL PUSHPOINTER8(C_LOC(sina_u))
    CALL PUSHREALARRAY(vort, (bd%ie-bd%is+3)*(bd%je-bd%js+3))
    CALL PUSHINTEGER(ie)
    !CALL PUSHPOINTER8(C_LOC(dyc))
    CALL PUSHREALARRAY(dt4)
    !CALL PUSHPOINTER8(C_LOC(sin_sg))
    CALL PUSHINTEGER(iep1)
    !CALL PUSHPOINTER8(C_LOC(cosa_v))
    CALL PUSHINTEGER(jep1)
    !CALL PUSHPOINTER8(C_LOC(cosa_u))
    CALL PUSHREALARRAY(fy2, (bd%ie-bd%is+3)*(bd%je-bd%js+4))
    CALL PUSHREALARRAY(fy1, (bd%ie-bd%is+3)*(bd%je-bd%js+4))
    !CALL PUSHPOINTER8(C_LOC(dy))
    !CALL PUSHPOINTER8(C_LOC(dxc))
    CALL PUSHINTEGER(npy)
    CALL PUSHINTEGER(npx)
    CALL PUSHINTEGER(js)
  END SUBROUTINE C_SW_FWD
!  Differentiation of c_sw in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2
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
!   gradient     of useful results: u v w delp ua uc ptc ut delpc
!                va vc vt divg_d wc pt
!   with respect to varying inputs: u v w delp ua uc ptc ut delpc
!                va vc vt divg_d wc pt
  SUBROUTINE C_SW_BWD(delpc, delpc_ad, delp, delp_ad, ptc, ptc_ad, pt, &
&   pt_ad, u, u_ad, v, v_ad, w, w_ad, uc, uc_ad, vc, vc_ad, ua, ua_ad, &
&   va, va_ad, wc, wc_ad, ut, ut_ad, vt, vt_ad, divg_d, divg_d_ad, nord&
&   , dt2, hydrostatic, dord4, bd, gridstruct, flagstruct)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: u&
&   , vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   u_ad, vc_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: v&
&   , uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   v_ad, uc_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: delp&
&   , pt, ua, va, ut, vt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delp_ad, pt_ad, ua_ad, va_ad, ut_ad, vt_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: delpc, ptc, wc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: delpc_ad, ptc_ad, &
&   wc_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1) :: divg_d
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1) :: divg_d_ad
    INTEGER, INTENT(IN) :: nord
    REAL, INTENT(IN) :: dt2
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: dord4
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+1) :: vort, ke
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+1) :: vort_ad, ke_ad
    REAL, DIMENSION(bd%is-1:bd%ie+2, bd%js-1:bd%je+1) :: fx, fx1, fx2
    REAL, DIMENSION(bd%is-1:bd%ie+2, bd%js-1:bd%je+1) :: fx_ad, fx1_ad, &
&   fx2_ad
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+2) :: fy, fy1, fy2
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+2) :: fy_ad, fy1_ad, &
&   fy2_ad
    REAL :: dt4
    INTEGER :: i, j, is2, ie1
    INTEGER :: iep1, jep1
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: npx, npy
    LOGICAL :: nested
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg, cos_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: dx, dy, dxc, dyc
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
    INTEGER :: branch
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_c_sw

    vort = 0.0
    ke = 0.0
    fx = 0.0
    fx1 = 0.0
    fx2 = 0.0
    fy = 0.0
    fy1 = 0.0
    fy2 = 0.0
    dt4 = 0.0
    is2 = 0
    ie1 = 0
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
    npx = 0
    npy = 0
    branch = 0

    CALL POPINTEGER(js)
    CALL POPINTEGER(npx)
    CALL POPINTEGER(npy)
    !CALL POPPOINTER8(cptr)
    dxc => gridstruct%dxc ! (/unknown_shape_in_c_sw/))
    !CALL POPPOINTER8(cptr)
    dy => gridstruct%dy ! (/unknown_shape_in_c_sw/))
    CALL POPREALARRAY(fy1, (bd%ie-bd%is+3)*(bd%je-bd%js+4))
    CALL POPREALARRAY(fy2, (bd%ie-bd%is+3)*(bd%je-bd%js+4))
    !CALL POPPOINTER8(cptr)
    cosa_u => gridstruct%cosa_u ! (/unknown_shape_in_c_sw/))
    CALL POPINTEGER(jep1)
    !CALL POPPOINTER8(cptr)
    cosa_v => gridstruct%cosa_v ! (/unknown_shape_in_c_sw/))
    CALL POPINTEGER(iep1)
    !CALL POPPOINTER8(cptr)
    sin_sg => gridstruct%sin_sg ! (/unknown_shape_in_c_sw/))
    CALL POPREALARRAY(dt4)
    !CALL POPPOINTER8(cptr)
    dyc => gridstruct%dyc ! (/unknown_shape_in_c_sw/))
    CALL POPINTEGER(ie)
    CALL POPREALARRAY(vort, (bd%ie-bd%is+3)*(bd%je-bd%js+3))
    !CALL POPPOINTER8(cptr)
    sina_u => gridstruct%sina_u ! (/unknown_shape_in_c_sw/))
    !CALL POPPOINTER8(cptr)
    sina_v => gridstruct%sina_v ! (/unknown_shape_in_c_sw/))
    CALL POPINTEGER(is)
    CALL POPREALARRAY(fx, (bd%ie-bd%is+4)*(bd%je-bd%js+3))
    CALL POPREALARRAY(fy, (bd%ie-bd%is+3)*(bd%je-bd%js+4))
    CALL POPINTEGER(je)
    CALL POPREALARRAY(fx1, (bd%ie-bd%is+4)*(bd%je-bd%js+3))
    CALL POPREALARRAY(fx2, (bd%ie-bd%is+4)*(bd%je-bd%js+3))
    ke_ad = 0.0
    fx_ad = 0.0
    fx1_ad = 0.0
    DO j=jep1,js,-1
      DO i=ie,is,-1
        temp_ad13 = gridstruct%rdyc(i, j)*vc_ad(i, j)
        fx1_ad(i, j) = fx1_ad(i, j) - fx(i, j)*vc_ad(i, j)
        fx_ad(i, j) = fx_ad(i, j) - fx1(i, j)*vc_ad(i, j)
        ke_ad(i, j-1) = ke_ad(i, j-1) + temp_ad13
        ke_ad(i, j) = ke_ad(i, j) - temp_ad13
      END DO
    END DO
    fy1_ad = 0.0
    fy_ad = 0.0
    DO j=je,js,-1
      DO i=iep1,is,-1
        temp_ad12 = gridstruct%rdxc(i, j)*uc_ad(i, j)
        fy1_ad(i, j) = fy1_ad(i, j) + fy(i, j)*uc_ad(i, j)
        fy_ad(i, j) = fy_ad(i, j) + fy1(i, j)*uc_ad(i, j)
        ke_ad(i-1, j) = ke_ad(i-1, j) + temp_ad12
        ke_ad(i, j) = ke_ad(i, j) - temp_ad12
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      vort_ad = 0.0
      DO j=jep1,js,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO i=ie,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY(fx(i, j))
              vort_ad(i+1, j) = vort_ad(i+1, j) + fx_ad(i, j)
              fx_ad(i, j) = 0.0
            ELSE
              CALL POPREALARRAY(fx(i, j))
              vort_ad(i, j) = vort_ad(i, j) + fx_ad(i, j)
              fx_ad(i, j) = 0.0
            END IF
            CALL POPREALARRAY(fx1(i, j))
            temp_ad11 = dt2*fx1_ad(i, j)/sina_v(i, j)
            u_ad(i, j) = u_ad(i, j) + temp_ad11
            vc_ad(i, j) = vc_ad(i, j) - cosa_v(i, j)*temp_ad11
            fx1_ad(i, j) = 0.0
          END DO
        ELSE
          DO i=ie,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY(fx(i, j))
              vort_ad(i+1, j) = vort_ad(i+1, j) + fx_ad(i, j)
              fx_ad(i, j) = 0.0
            ELSE
              CALL POPREALARRAY(fx(i, j))
              vort_ad(i, j) = vort_ad(i, j) + fx_ad(i, j)
              fx_ad(i, j) = 0.0
            END IF
            CALL POPREALARRAY(fx1(i, j))
            u_ad(i, j) = u_ad(i, j) + dt2*fx1_ad(i, j)
            fx1_ad(i, j) = 0.0
          END DO
        END IF
      END DO
      DO j=je,js,-1
        DO i=iep1,is,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(fy(i, j))
            vort_ad(i, j+1) = vort_ad(i, j+1) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
          ELSE
            CALL POPREALARRAY(fy(i, j))
            vort_ad(i, j) = vort_ad(i, j) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(fy1(i, j))
            v_ad(i, j) = v_ad(i, j) + dt2*fy1_ad(i, j)
            fy1_ad(i, j) = 0.0
          ELSE
            CALL POPREALARRAY(fy1(i, j))
            temp_ad10 = dt2*fy1_ad(i, j)/sina_u(i, j)
            v_ad(i, j) = v_ad(i, j) + temp_ad10
            uc_ad(i, j) = uc_ad(i, j) - cosa_u(i, j)*temp_ad10
            fy1_ad(i, j) = 0.0
          END IF
        END DO
      END DO
    ELSE
      vort_ad = 0.0
      DO j=jep1,js,-1
        DO i=ie,is,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(fx(i, j))
            vort_ad(i+1, j) = vort_ad(i+1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0
          ELSE
            CALL POPREALARRAY(fx(i, j))
            vort_ad(i, j) = vort_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0
          END IF
          CALL POPREALARRAY(fx1(i, j))
          temp_ad9 = dt2*fx1_ad(i, j)/sina_v(i, j)
          u_ad(i, j) = u_ad(i, j) + temp_ad9
          vc_ad(i, j) = vc_ad(i, j) - cosa_v(i, j)*temp_ad9
          fx1_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je,js,-1
        DO i=iep1,is,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(fy(i, j))
            vort_ad(i, j+1) = vort_ad(i, j+1) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
          ELSE
            CALL POPREALARRAY(fy(i, j))
            vort_ad(i, j) = vort_ad(i, j) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
          END IF
          CALL POPREALARRAY(fy1(i, j))
          temp_ad8 = dt2*fy1_ad(i, j)/sina_u(i, j)
          v_ad(i, j) = v_ad(i, j) + temp_ad8
          uc_ad(i, j) = uc_ad(i, j) - cosa_u(i, j)*temp_ad8
          fy1_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(vort(i, j))
        vort_ad(i, j) = gridstruct%rarea_c(i, j)*vort_ad(i, j)
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      CALL POPREALARRAY(vort(1, npy))
      fy_ad(0, npy) = fy_ad(0, npy) + vort_ad(1, npy)
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(vort(npx, npy))
      fy_ad(npx, npy) = fy_ad(npx, npy) - vort_ad(npx, npy)
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(vort(npx, 1))
      fy_ad(npx, 1) = fy_ad(npx, 1) - vort_ad(npx, 1)
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(vort(1, 1))
      fy_ad(0, 1) = fy_ad(0, 1) + vort_ad(1, 1)
    END IF
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(vort(i, j))
        fx_ad(i, j-1) = fx_ad(i, j-1) + vort_ad(i, j)
        fx_ad(i, j) = fx_ad(i, j) - vort_ad(i, j)
        fy_ad(i, j) = fy_ad(i, j) + vort_ad(i, j)
        fy_ad(i-1, j) = fy_ad(i-1, j) - vort_ad(i, j)
        vort_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is-1,-1
        CALL POPREALARRAY(fy(i, j))
        vc_ad(i, j) = vc_ad(i, j) + dyc(i, j)*fy_ad(i, j)
        fy_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js-1,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(fx(i, j))
        uc_ad(i, j) = uc_ad(i, j) + dxc(i, j)*fx_ad(i, j)
        fx_ad(i, j) = 0.0
      END DO
    END DO
    DO j=jep1,js-1,-1
      DO i=iep1,is-1,-1
        CALL POPREALARRAY(ke(i, j))
        temp_ad7 = dt4*ke_ad(i, j)
        ua_ad(i, j) = ua_ad(i, j) + ke(i, j)*temp_ad7
        va_ad(i, j) = va_ad(i, j) + vort(i, j)*temp_ad7
        vort_ad(i, j) = vort_ad(i, j) + va(i, j)*temp_ad7
        ke_ad(i, j) = ua(i, j)*temp_ad7
      END DO
    END DO
    cos_sg => gridstruct%cos_sg
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=jep1,js-1,-1
        DO i=iep1,is-1,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            vc_ad(i, j+1) = vc_ad(i, j+1) + vort_ad(i, j)
            vort_ad(i, j) = 0.0
          ELSE
            vc_ad(i, j) = vc_ad(i, j) + vort_ad(i, j)
            vort_ad(i, j) = 0.0
          END IF
        END DO
      END DO
      DO j=jep1,js-1,-1
        DO i=iep1,is-1,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            uc_ad(i+1, j) = uc_ad(i+1, j) + ke_ad(i, j)
            ke_ad(i, j) = 0.0
          ELSE
            uc_ad(i, j) = uc_ad(i, j) + ke_ad(i, j)
            ke_ad(i, j) = 0.0
          END IF
        END DO
      END DO
    ELSE
      DO j=jep1,js-1,-1
        DO i=iep1,is-1,-1
          CALL POPCONTROL(3,branch)
          IF (branch .LT. 3) THEN
            IF (branch .EQ. 0) THEN
              vc_ad(i, j+1) = vc_ad(i, j+1) + vort_ad(i, j)
              vort_ad(i, j) = 0.0
            ELSE IF (branch .EQ. 1) THEN
              vc_ad(i, npy) = vc_ad(i, npy) + sin_sg(i, npy-1, 4)*&
&               vort_ad(i, j)
              u_ad(i, npy) = u_ad(i, npy) + cos_sg(i, npy-1, 4)*vort_ad(&
&               i, j)
              vort_ad(i, j) = 0.0
            ELSE
              vc_ad(i, 1) = vc_ad(i, 1) + sin_sg(i, 0, 4)*vort_ad(i, 0)
              u_ad(i, 1) = u_ad(i, 1) + cos_sg(i, 0, 4)*vort_ad(i, 0)
              vort_ad(i, 0) = 0.0
            END IF
          ELSE IF (branch .EQ. 3) THEN
            vc_ad(i, j) = vc_ad(i, j) + vort_ad(i, j)
            vort_ad(i, j) = 0.0
          ELSE IF (branch .EQ. 4) THEN
            vc_ad(i, npy) = vc_ad(i, npy) + sin_sg(i, npy, 2)*vort_ad(i&
&             , j)
            u_ad(i, npy) = u_ad(i, npy) + cos_sg(i, npy, 2)*vort_ad(i, j&
&             )
            vort_ad(i, j) = 0.0
          ELSE
            vc_ad(i, 1) = vc_ad(i, 1) + sin_sg(i, 1, 2)*vort_ad(i, 1)
            u_ad(i, 1) = u_ad(i, 1) + cos_sg(i, 1, 2)*vort_ad(i, 1)
            vort_ad(i, 1) = 0.0
          END IF
        END DO
      END DO
      DO j=jep1,js-1,-1
        DO i=iep1,is-1,-1
          CALL POPCONTROL(3,branch)
          IF (branch .LT. 3) THEN
            IF (branch .EQ. 0) THEN
              uc_ad(i+1, j) = uc_ad(i+1, j) + ke_ad(i, j)
              ke_ad(i, j) = 0.0
            ELSE IF (branch .EQ. 1) THEN
              uc_ad(npx, j) = uc_ad(npx, j) + sin_sg(npx-1, j, 3)*ke_ad(&
&               i, j)
              v_ad(npx, j) = v_ad(npx, j) + cos_sg(npx-1, j, 3)*ke_ad(i&
&               , j)
              ke_ad(i, j) = 0.0
            ELSE
              uc_ad(1, j) = uc_ad(1, j) + sin_sg(0, j, 3)*ke_ad(0, j)
              v_ad(1, j) = v_ad(1, j) + cos_sg(0, j, 3)*ke_ad(0, j)
              ke_ad(0, j) = 0.0
            END IF
          ELSE IF (branch .EQ. 3) THEN
            uc_ad(i, j) = uc_ad(i, j) + ke_ad(i, j)
            ke_ad(i, j) = 0.0
          ELSE IF (branch .EQ. 4) THEN
            uc_ad(npx, j) = uc_ad(npx, j) + sin_sg(npx, j, 1)*ke_ad(i, j&
&             )
            v_ad(npx, j) = v_ad(npx, j) + cos_sg(npx, j, 1)*ke_ad(i, j)
            ke_ad(i, j) = 0.0
          ELSE
            uc_ad(1, j) = uc_ad(1, j) + sin_sg(1, j, 1)*ke_ad(1, j)
            v_ad(1, j) = v_ad(1, j) + cos_sg(1, j, 1)*ke_ad(1, j)
            ke_ad(1, j) = 0.0
          END IF
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=jep1,js-1,-1
        DO i=iep1,is-1,-1
          CALL POPREALARRAY(ptc(i, j))
          temp_ad = ptc_ad(i, j)/delpc(i, j)
          temp_ad0 = gridstruct%rarea(i, j)*temp_ad
          pt_ad(i, j) = pt_ad(i, j) + delp(i, j)*temp_ad
          fx_ad(i, j) = fx_ad(i, j) + temp_ad0
          fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad0
          fy_ad(i, j) = fy_ad(i, j) + temp_ad0
          fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad0
          delpc_ad(i, j) = delpc_ad(i, j) - (pt(i, j)*delp(i, j)+&
&           gridstruct%rarea(i, j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j&
&           +1)))*temp_ad/delpc(i, j)
          delp_ad(i, j) = delp_ad(i, j) + delpc_ad(i, j) + pt(i, j)*&
&           temp_ad
          ptc_ad(i, j) = 0.0
          CALL POPREALARRAY(delpc(i, j))
          temp_ad1 = gridstruct%rarea(i, j)*delpc_ad(i, j)
          fx1_ad(i, j) = fx1_ad(i, j) + temp_ad1
          fx1_ad(i+1, j) = fx1_ad(i+1, j) - temp_ad1
          fy1_ad(i, j) = fy1_ad(i, j) + temp_ad1
          fy1_ad(i, j+1) = fy1_ad(i, j+1) - temp_ad1
          delpc_ad(i, j) = 0.0
        END DO
      END DO
      DO j=jep1+1,js-1,-1
        DO i=iep1,is-1,-1
          CALL POPREALARRAY(fy(i, j))
          fy1_ad(i, j) = fy1_ad(i, j) + fy(i, j)*fy_ad(i, j)
          fy_ad(i, j) = fy1(i, j)*fy_ad(i, j)
          CALL POPREALARRAY(fy1(i, j))
          vt_ad(i, j) = vt_ad(i, j) + fy1(i, j)*fy1_ad(i, j)
          fy1_ad(i, j) = vt(i, j)*fy1_ad(i, j)
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            pt_ad(i, j-1) = pt_ad(i, j-1) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
            delp_ad(i, j-1) = delp_ad(i, j-1) + fy1_ad(i, j)
            fy1_ad(i, j) = 0.0
          ELSE
            pt_ad(i, j) = pt_ad(i, j) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
            delp_ad(i, j) = delp_ad(i, j) + fy1_ad(i, j)
            fy1_ad(i, j) = 0.0
          END IF
        END DO
      END DO
      fx2_ad = 0.0
    ELSE
      fy2_ad = 0.0
      fx2_ad = 0.0
      DO j=je+1,js-1,-1
        DO i=ie+1,is-1,-1
          temp_ad4 = ptc_ad(i, j)/delpc(i, j)
          CALL POPREALARRAY(wc(i, j))
          temp_ad2 = wc_ad(i, j)/delpc(i, j)
          temp_ad3 = gridstruct%rarea(i, j)*temp_ad2
          w_ad(i, j) = w_ad(i, j) + delp(i, j)*temp_ad2
          fx2_ad(i, j) = fx2_ad(i, j) + temp_ad3
          fx2_ad(i+1, j) = fx2_ad(i+1, j) - temp_ad3
          fy2_ad(i, j) = fy2_ad(i, j) + temp_ad3
          fy2_ad(i, j+1) = fy2_ad(i, j+1) - temp_ad3
          delpc_ad(i, j) = delpc_ad(i, j) - (pt(i, j)*delp(i, j)+&
&           gridstruct%rarea(i, j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j&
&           +1)))*temp_ad4/delpc(i, j) - (w(i, j)*delp(i, j)+gridstruct%&
&           rarea(i, j)*(fx2(i, j)-fx2(i+1, j)+fy2(i, j)-fy2(i, j+1)))*&
&           temp_ad2/delpc(i, j)
          delp_ad(i, j) = delp_ad(i, j) + pt(i, j)*temp_ad4 + delpc_ad(i&
&           , j) + w(i, j)*temp_ad2
          wc_ad(i, j) = 0.0
          CALL POPREALARRAY(ptc(i, j))
          temp_ad5 = gridstruct%rarea(i, j)*temp_ad4
          pt_ad(i, j) = pt_ad(i, j) + delp(i, j)*temp_ad4
          fx_ad(i, j) = fx_ad(i, j) + temp_ad5
          fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad5
          fy_ad(i, j) = fy_ad(i, j) + temp_ad5
          fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad5
          ptc_ad(i, j) = 0.0
          CALL POPREALARRAY(delpc(i, j))
          temp_ad6 = gridstruct%rarea(i, j)*delpc_ad(i, j)
          fx1_ad(i, j) = fx1_ad(i, j) + temp_ad6
          fx1_ad(i+1, j) = fx1_ad(i+1, j) - temp_ad6
          fy1_ad(i, j) = fy1_ad(i, j) + temp_ad6
          fy1_ad(i, j+1) = fy1_ad(i, j+1) - temp_ad6
          delpc_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+2,js-1,-1
        DO i=ie+1,is-1,-1
          CALL POPREALARRAY(fy2(i, j))
          CALL POPREALARRAY(fy(i, j))
          fy1_ad(i, j) = fy1_ad(i, j) + fy(i, j)*fy_ad(i, j) + fy2(i, j)&
&           *fy2_ad(i, j)
          fy2_ad(i, j) = fy1(i, j)*fy2_ad(i, j)
          fy_ad(i, j) = fy1(i, j)*fy_ad(i, j)
          CALL POPREALARRAY(fy1(i, j))
          vt_ad(i, j) = vt_ad(i, j) + fy1(i, j)*fy1_ad(i, j)
          fy1_ad(i, j) = vt(i, j)*fy1_ad(i, j)
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            w_ad(i, j-1) = w_ad(i, j-1) + fy2_ad(i, j)
            fy2_ad(i, j) = 0.0
            pt_ad(i, j-1) = pt_ad(i, j-1) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
            delp_ad(i, j-1) = delp_ad(i, j-1) + fy1_ad(i, j)
            fy1_ad(i, j) = 0.0
          ELSE
            w_ad(i, j) = w_ad(i, j) + fy2_ad(i, j)
            fy2_ad(i, j) = 0.0
            pt_ad(i, j) = pt_ad(i, j) + fy_ad(i, j)
            fy_ad(i, j) = 0.0
            delp_ad(i, j) = delp_ad(i, j) + fy1_ad(i, j)
            fy1_ad(i, j) = 0.0
          END IF
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) CALL FILL_4CORNERS_BWD(w, w_ad, 2, bd, npx, npy&
&                                         , sw_corner, se_corner, &
&                                         ne_corner, nw_corner)
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) CALL FILL2_4CORNERS_BWD(delp, delp_ad, pt, pt_ad&
&                                        , 2, bd, npx, npy, sw_corner, &
&                                        se_corner, ne_corner, nw_corner&
&                                       )
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=jep1,js-1,-1
        DO i=ie+2,is-1,-1
          CALL POPREALARRAY(fx(i, j))
          fx1_ad(i, j) = fx1_ad(i, j) + fx(i, j)*fx_ad(i, j)
          fx_ad(i, j) = fx1(i, j)*fx_ad(i, j)
          CALL POPREALARRAY(fx1(i, j))
          ut_ad(i, j) = ut_ad(i, j) + fx1(i, j)*fx1_ad(i, j)
          fx1_ad(i, j) = ut(i, j)*fx1_ad(i, j)
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            pt_ad(i-1, j) = pt_ad(i-1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0
            delp_ad(i-1, j) = delp_ad(i-1, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0
          ELSE
            pt_ad(i, j) = pt_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0
            delp_ad(i, j) = delp_ad(i, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0
          END IF
        END DO
      END DO
    ELSE
      DO j=je+1,js-1,-1
        DO i=ie+2,is-1,-1
          CALL POPREALARRAY(fx2(i, j))
          CALL POPREALARRAY(fx(i, j))
          fx1_ad(i, j) = fx1_ad(i, j) + fx(i, j)*fx_ad(i, j) + fx2(i, j)&
&           *fx2_ad(i, j)
          fx2_ad(i, j) = fx1(i, j)*fx2_ad(i, j)
          fx_ad(i, j) = fx1(i, j)*fx_ad(i, j)
          CALL POPREALARRAY(fx1(i, j))
          ut_ad(i, j) = ut_ad(i, j) + fx1(i, j)*fx1_ad(i, j)
          fx1_ad(i, j) = ut(i, j)*fx1_ad(i, j)
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            w_ad(i-1, j) = w_ad(i-1, j) + fx2_ad(i, j)
            fx2_ad(i, j) = 0.0
            pt_ad(i-1, j) = pt_ad(i-1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0
            delp_ad(i-1, j) = delp_ad(i-1, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0
          ELSE
            w_ad(i, j) = w_ad(i, j) + fx2_ad(i, j)
            fx2_ad(i, j) = 0.0
            pt_ad(i, j) = pt_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0
            delp_ad(i, j) = delp_ad(i, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0
          END IF
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) CALL FILL_4CORNERS_BWD(w, w_ad, 1, bd, npx, npy&
&                                         , sw_corner, se_corner, &
&                                         ne_corner, nw_corner)
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) CALL FILL2_4CORNERS_BWD(delp, delp_ad, pt, pt_ad&
&                                        , 1, bd, npx, npy, sw_corner, &
&                                        se_corner, ne_corner, nw_corner&
&                                       )
    dx => gridstruct%dx
    DO j=je+2,js-1,-1
      DO i=iep1,is-1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(vt(i, j))
          vt_ad(i, j) = dx(i, j)*sin_sg(i, j, 2)*dt2*vt_ad(i, j)
        ELSE
          CALL POPREALARRAY(vt(i, j))
          vt_ad(i, j) = dx(i, j)*sin_sg(i, j-1, 4)*dt2*vt_ad(i, j)
        END IF
      END DO
    END DO
    DO j=jep1,js-1,-1
      DO i=iep1+1,is-1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(ut(i, j))
          ut_ad(i, j) = dy(i, j)*sin_sg(i, j, 1)*dt2*ut_ad(i, j)
        ELSE
          CALL POPREALARRAY(ut(i, j))
          ut_ad(i, j) = dy(i, j)*sin_sg(i-1, j, 3)*dt2*ut_ad(i, j)
        END IF
      END DO
    END DO
    CALL POPCONTROL(2,branch)
    IF (branch .NE. 0) THEN
      IF (branch .EQ. 1) THEN
        CALL DIVERGENCE_CORNER_BWD(u, u_ad, v, v_ad, ua, ua_ad, va, &
&                            va_ad, divg_d, divg_d_ad, gridstruct, &
&                            flagstruct, bd)
      ELSE
        CALL DIVERGENCE_CORNER_NEST_BWD(u, u_ad, v, v_ad, ua, ua_ad, va&
&                                 , va_ad, divg_d, divg_d_ad, gridstruct&
&                                 , flagstruct, bd)
        divg_d_ad = 0.0
      END IF
    END IF
    CALL D2A2C_VECT_BWD(u, u_ad, v, v_ad, ua, ua_ad, va, va_ad, uc, &
&                 uc_ad, vc, vc_ad, ut, ut_ad, vt, vt_ad, dord4, &
&                 gridstruct, bd, npx, npy, nested, flagstruct%grid_type&
&                )
  END SUBROUTINE C_SW_BWD
  SUBROUTINE C_SW(delpc, delp, ptc, pt, u, v, w, uc, vc, ua, va, wc, ut&
&   , vt, divg_d, nord, dt2, hydrostatic, dord4, bd, gridstruct, &
&   flagstruct)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: u&
&   , vc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: v&
&   , uc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: delp&
&   , pt, ua, va, ut, vt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: delpc&
&   , ptc, wc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   divg_d
    INTEGER, INTENT(IN) :: nord
    REAL, INTENT(IN) :: dt2
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: dord4
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! Local:
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+1) :: vort, ke
    REAL, DIMENSION(bd%is-1:bd%ie+2, bd%js-1:bd%je+1) :: fx, fx1, fx2
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+2) :: fy, fy1, fy2
    REAL :: dt4
    INTEGER :: i, j, is2, ie1
    INTEGER :: iep1, jep1
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: npx, npy
    LOGICAL :: nested
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg, cos_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: dx, dy, dxc, dyc
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    sin_sg => gridstruct%sin_sg
    cos_sg => gridstruct%cos_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    sina_u => gridstruct%sina_u
    sina_v => gridstruct%sina_v
    dx => gridstruct%dx
    dy => gridstruct%dy
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    sw_corner = gridstruct%sw_corner
    se_corner = gridstruct%se_corner
    nw_corner = gridstruct%nw_corner
    ne_corner = gridstruct%ne_corner
    iep1 = ie + 1
    jep1 = je + 1
    CALL D2A2C_VECT(u, v, ua, va, uc, vc, ut, vt, dord4, gridstruct, bd&
&             , npx, npy, nested, flagstruct%grid_type)
    IF (nord .GT. 0) THEN
      IF (nested) THEN
        CALL DIVERGENCE_CORNER_NEST(u, v, ua, va, divg_d, gridstruct, &
&                             flagstruct, bd)
      ELSE
        CALL DIVERGENCE_CORNER(u, v, ua, va, divg_d, gridstruct, &
&                        flagstruct, bd)
      END IF
    END IF
    DO j=js-1,jep1
      DO i=is-1,iep1+1
        IF (ut(i, j) .GT. 0.) THEN
          ut(i, j) = dt2*ut(i, j)*dy(i, j)*sin_sg(i-1, j, 3)
        ELSE
          ut(i, j) = dt2*ut(i, j)*dy(i, j)*sin_sg(i, j, 1)
        END IF
      END DO
    END DO
    DO j=js-1,je+2
      DO i=is-1,iep1
        IF (vt(i, j) .GT. 0.) THEN
          vt(i, j) = dt2*vt(i, j)*dx(i, j)*sin_sg(i, j-1, 4)
        ELSE
          vt(i, j) = dt2*vt(i, j)*dx(i, j)*sin_sg(i, j, 2)
        END IF
      END DO
    END DO
!----------------
! Transport delp:
!----------------
! Xdir:
    IF (flagstruct%grid_type .LT. 3 .AND. (.NOT.nested)) CALL &
&     FILL2_4CORNERS(delp, pt, 1, bd, npx, npy, sw_corner, se_corner, &
&              ne_corner, nw_corner)
    IF (hydrostatic) THEN
      DO j=js-1,jep1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1(i, j) = delp(i-1, j)
            fx(i, j) = pt(i-1, j)
          ELSE
            fx1(i, j) = delp(i, j)
            fx(i, j) = pt(i, j)
          END IF
          fx1(i, j) = ut(i, j)*fx1(i, j)
          fx(i, j) = fx1(i, j)*fx(i, j)
        END DO
      END DO
    ELSE
      IF (flagstruct%grid_type .LT. 3) CALL FILL_4CORNERS(w, 1, bd, npx&
&                                                   , npy, sw_corner, &
&                                                   se_corner, ne_corner&
&                                                   , nw_corner)
      DO j=js-1,je+1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1(i, j) = delp(i-1, j)
            fx(i, j) = pt(i-1, j)
            fx2(i, j) = w(i-1, j)
          ELSE
            fx1(i, j) = delp(i, j)
            fx(i, j) = pt(i, j)
            fx2(i, j) = w(i, j)
          END IF
          fx1(i, j) = ut(i, j)*fx1(i, j)
          fx(i, j) = fx1(i, j)*fx(i, j)
          fx2(i, j) = fx1(i, j)*fx2(i, j)
        END DO
      END DO
    END IF
! Ydir:
    IF (flagstruct%grid_type .LT. 3 .AND. (.NOT.nested)) CALL &
&     FILL2_4CORNERS(delp, pt, 2, bd, npx, npy, sw_corner, se_corner, &
&              ne_corner, nw_corner)
    IF (hydrostatic) THEN
      DO j=js-1,jep1+1
        DO i=is-1,iep1
          IF (vt(i, j) .GT. 0.) THEN
            fy1(i, j) = delp(i, j-1)
            fy(i, j) = pt(i, j-1)
          ELSE
            fy1(i, j) = delp(i, j)
            fy(i, j) = pt(i, j)
          END IF
          fy1(i, j) = vt(i, j)*fy1(i, j)
          fy(i, j) = fy1(i, j)*fy(i, j)
        END DO
      END DO
      DO j=js-1,jep1
        DO i=is-1,iep1
          delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+(fy1(i, j)-&
&           fy1(i, j+1)))*gridstruct%rarea(i, j)
          ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+(fy(i, j&
&           )-fy(i, j+1)))*gridstruct%rarea(i, j))/delpc(i, j)
        END DO
      END DO
    ELSE
      IF (flagstruct%grid_type .LT. 3) CALL FILL_4CORNERS(w, 2, bd, npx&
&                                                   , npy, sw_corner, &
&                                                   se_corner, ne_corner&
&                                                   , nw_corner)
      DO j=js-1,je+2
        DO i=is-1,ie+1
          IF (vt(i, j) .GT. 0.) THEN
            fy1(i, j) = delp(i, j-1)
            fy(i, j) = pt(i, j-1)
            fy2(i, j) = w(i, j-1)
          ELSE
            fy1(i, j) = delp(i, j)
            fy(i, j) = pt(i, j)
            fy2(i, j) = w(i, j)
          END IF
          fy1(i, j) = vt(i, j)*fy1(i, j)
          fy(i, j) = fy1(i, j)*fy(i, j)
          fy2(i, j) = fy1(i, j)*fy2(i, j)
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+(fy1(i, j)-&
&           fy1(i, j+1)))*gridstruct%rarea(i, j)
          ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+(fy(i, j&
&           )-fy(i, j+1)))*gridstruct%rarea(i, j))/delpc(i, j)
          wc(i, j) = (w(i, j)*delp(i, j)+(fx2(i, j)-fx2(i+1, j)+(fy2(i, &
&           j)-fy2(i, j+1)))*gridstruct%rarea(i, j))/delpc(i, j)
        END DO
      END DO
    END IF
!------------
! Compute KE:
!------------
!Since uc = u*, i.e. the covariant wind perpendicular to the face edge, if we want to compute kinetic energy we will need the tru
!e coordinate-parallel covariant wind, computed through u = uc*sina + v*cosa.
!Use the alpha for the cell KE is being computed in.
!!! TO DO:
!!! Need separate versions for nesting/single-tile
!!!   and for cubed-sphere
    IF (nested .OR. flagstruct%grid_type .GE. 3) THEN
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (ua(i, j) .GT. 0.) THEN
            ke(i, j) = uc(i, j)
          ELSE
            ke(i, j) = uc(i+1, j)
          END IF
        END DO
      END DO
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (va(i, j) .GT. 0.) THEN
            vort(i, j) = vc(i, j)
          ELSE
            vort(i, j) = vc(i, j+1)
          END IF
        END DO
      END DO
    ELSE
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (ua(i, j) .GT. 0.) THEN
            IF (i .EQ. 1) THEN
              ke(1, j) = uc(1, j)*sin_sg(1, j, 1) + v(1, j)*cos_sg(1, j&
&               , 1)
            ELSE IF (i .EQ. npx) THEN
              ke(i, j) = uc(npx, j)*sin_sg(npx, j, 1) + v(npx, j)*cos_sg&
&               (npx, j, 1)
            ELSE
              ke(i, j) = uc(i, j)
            END IF
          ELSE IF (i .EQ. 0) THEN
            ke(0, j) = uc(1, j)*sin_sg(0, j, 3) + v(1, j)*cos_sg(0, j, 3&
&             )
          ELSE IF (i .EQ. npx - 1) THEN
            ke(i, j) = uc(npx, j)*sin_sg(npx-1, j, 3) + v(npx, j)*cos_sg&
&             (npx-1, j, 3)
          ELSE
            ke(i, j) = uc(i+1, j)
          END IF
        END DO
      END DO
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (va(i, j) .GT. 0.) THEN
            IF (j .EQ. 1) THEN
              vort(i, 1) = vc(i, 1)*sin_sg(i, 1, 2) + u(i, 1)*cos_sg(i, &
&               1, 2)
            ELSE IF (j .EQ. npy) THEN
              vort(i, j) = vc(i, npy)*sin_sg(i, npy, 2) + u(i, npy)*&
&               cos_sg(i, npy, 2)
            ELSE
              vort(i, j) = vc(i, j)
            END IF
          ELSE IF (j .EQ. 0) THEN
            vort(i, 0) = vc(i, 1)*sin_sg(i, 0, 4) + u(i, 1)*cos_sg(i, 0&
&             , 4)
          ELSE IF (j .EQ. npy - 1) THEN
            vort(i, j) = vc(i, npy)*sin_sg(i, npy-1, 4) + u(i, npy)*&
&             cos_sg(i, npy-1, 4)
          ELSE
            vort(i, j) = vc(i, j+1)
          END IF
        END DO
      END DO
    END IF
    dt4 = 0.5*dt2
    DO j=js-1,jep1
      DO i=is-1,iep1
        ke(i, j) = dt4*(ua(i, j)*ke(i, j)+va(i, j)*vort(i, j))
      END DO
    END DO
!------------------------------
! Compute circulation on C grid
!------------------------------
! To consider using true co-variant winds at face edges?
    DO j=js-1,je+1
      DO i=is,ie+1
        fx(i, j) = uc(i, j)*dxc(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is-1,ie+1
        fy(i, j) = vc(i, j)*dyc(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie+1
        vort(i, j) = fx(i, j-1) - fx(i, j) + (fy(i, j)-fy(i-1, j))
      END DO
    END DO
! Remove the extra term at the corners:
    IF (sw_corner) vort(1, 1) = vort(1, 1) + fy(0, 1)
    IF (se_corner) vort(npx, 1) = vort(npx, 1) - fy(npx, 1)
    IF (ne_corner) vort(npx, npy) = vort(npx, npy) - fy(npx, npy)
    IF (nw_corner) vort(1, npy) = vort(1, npy) + fy(0, npy)
!----------------------------
! Compute absolute vorticity
!----------------------------
    DO j=js,je+1
      DO i=is,ie+1
        vort(i, j) = gridstruct%fc(i, j) + gridstruct%rarea_c(i, j)*vort&
&         (i, j)
      END DO
    END DO
!----------------------------------
! Transport absolute vorticity:
!----------------------------------
!To go from v to contravariant v at the edges, we divide by sin_sg;
! but we then must multiply by sin_sg to get the proper flux.
! These cancel, leaving us with fy1 = dt2*v at the edges.
! (For the same reason we only divide by sin instead of sin**2 in the interior)
!! TO DO: separate versions for nesting/single-tile and cubed-sphere
    IF (nested .OR. flagstruct%grid_type .GE. 3) THEN
      DO j=js,je
        DO i=is,iep1
          fy1(i, j) = dt2*(v(i, j)-uc(i, j)*cosa_u(i, j))/sina_u(i, j)
          IF (fy1(i, j) .GT. 0.) THEN
            fy(i, j) = vort(i, j)
          ELSE
            fy(i, j) = vort(i, j+1)
          END IF
        END DO
      END DO
      DO j=js,jep1
        DO i=is,ie
          fx1(i, j) = dt2*(u(i, j)-vc(i, j)*cosa_v(i, j))/sina_v(i, j)
          IF (fx1(i, j) .GT. 0.) THEN
            fx(i, j) = vort(i, j)
          ELSE
            fx(i, j) = vort(i+1, j)
          END IF
        END DO
      END DO
    ELSE
      DO j=js,je
!DEC$ VECTOR ALWAYS
        DO i=is,iep1
          IF (i .EQ. 1 .OR. i .EQ. npx) THEN
            fy1(i, j) = dt2*v(i, j)
          ELSE
            fy1(i, j) = dt2*(v(i, j)-uc(i, j)*cosa_u(i, j))/sina_u(i, j)
          END IF
          IF (fy1(i, j) .GT. 0.) THEN
            fy(i, j) = vort(i, j)
          ELSE
            fy(i, j) = vort(i, j+1)
          END IF
        END DO
      END DO
      DO j=js,jep1
        IF (j .EQ. 1 .OR. j .EQ. npy) THEN
!DEC$ VECTOR ALWAYS
          DO i=is,ie
            fx1(i, j) = dt2*u(i, j)
            IF (fx1(i, j) .GT. 0.) THEN
              fx(i, j) = vort(i, j)
            ELSE
              fx(i, j) = vort(i+1, j)
            END IF
          END DO
        ELSE
!DEC$ VECTOR ALWAYS
          DO i=is,ie
            fx1(i, j) = dt2*(u(i, j)-vc(i, j)*cosa_v(i, j))/sina_v(i, j)
            IF (fx1(i, j) .GT. 0.) THEN
              fx(i, j) = vort(i, j)
            ELSE
              fx(i, j) = vort(i+1, j)
            END IF
          END DO
        END IF
      END DO
    END IF
! Update time-centered winds on the C-Grid
    DO j=js,je
      DO i=is,iep1
        uc(i, j) = uc(i, j) + fy1(i, j)*fy(i, j) + gridstruct%rdxc(i, j)&
&         *(ke(i-1, j)-ke(i, j))
      END DO
    END DO
    DO j=js,jep1
      DO i=is,ie
        vc(i, j) = vc(i, j) - fx1(i, j)*fx(i, j) + gridstruct%rdyc(i, j)&
&         *(ke(i, j-1)-ke(i, j))
      END DO
    END DO
  END SUBROUTINE C_SW
!  Differentiation of d_sw in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b
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
!   gradient     of useful results: yfx_adv q crx_adv u v w delp
!                ua xfx_adv uc ptc xflux cry_adv delpc va vc yflux
!                divg_d z_rat heat_source pt cx cy dpx
!   with respect to varying inputs: yfx_adv q crx_adv u v w delp
!                ua xfx_adv uc ptc xflux cry_adv delpc va vc yflux
!                divg_d z_rat heat_source pt cx cy dpx
!     d_sw :: D-Grid Shallow Water Routine
  SUBROUTINE D_SW_FWD(delpc, delp, ptc, pt, u, v, w, uc, vc, ua, va, &
&   divg_d, xflux, yflux, cx, cy, crx_adv, cry_adv, xfx_adv, yfx_adv, &
&   q_con, z_rat, kgb, heat_source, dpx, zvir, sphum, nq, q, k, km, &
&   inline_q, dt, hord_tr, hord_mt, hord_vt, hord_tm, hord_dp, nord, &
&   nord_v, nord_w, nord_t, dddmp, d2_bg, d4_bg, damp_v, damp_w, damp_t&
&   , d_con, hydrostatic, gridstruct, flagstruct, bd, hord_tr_pert, &
&   hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_dp_pert, split_damp, &
&   nord_pert, nord_v_pert, nord_w_pert, nord_t_pert, dddmp_pert, &
&   d2_bg_pert, d4_bg_pert, damp_v_pert, damp_w_pert, damp_t_pert)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: hord_tr, hord_mt, hord_vt, hord_tm, hord_dp
! nord=1 divergence damping; (del-4) or 3 (del-8)
    INTEGER, INTENT(IN) :: nord
! vorticity damping
    INTEGER, INTENT(IN) :: nord_v
! vertical velocity
    INTEGER, INTENT(IN) :: nord_w
! pt
    INTEGER, INTENT(IN) :: nord_t
    INTEGER, INTENT(IN) :: sphum, nq, k, km
    REAL, INTENT(IN) :: dt, dddmp, d2_bg, d4_bg, d_con
    REAL, INTENT(IN) :: zvir
    REAL, INTENT(IN) :: damp_v, damp_w, damp_t, kgb
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: hord_tr_pert, hord_mt_pert, hord_vt_pert, &
&   hord_tm_pert, hord_dp_pert, nord_pert, nord_v_pert, nord_w_pert, &
&   nord_t_pert
    LOGICAL, INTENT(IN) :: split_damp
    REAL, INTENT(IN) :: dddmp_pert, d2_bg_pert, d4_bg_pert, damp_v_pert&
&   , damp_w_pert, damp_t_pert
! divergence
    REAL, INTENT(INOUT) :: divg_d(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: z_rat
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: delp&
&   , pt, ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   q_con
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: u&
&   , vc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: v&
&   , uc
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, km, nq)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: delpc, ptc
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: heat_source
    REAL(kind=8), DIMENSION(bd%is:bd%ie, bd%js:bd%je), INTENT(INOUT) :: &
&   dpx
! The flux capacitors:
    REAL, INTENT(INOUT) :: xflux(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: yflux(bd%is:bd%ie, bd%js:bd%je+1)
!------------------------
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1)
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: inline_q
    REAL, DIMENSION(bd%is:bd%ie+1, bd%jsd:bd%jed) :: crx_adv, xfx_adv
    REAL, DIMENSION(bd%isd:bd%ied, bd%js:bd%je+1) :: cry_adv, yfx_adv
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! Local:
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL :: ut(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!---
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: fy2(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!  work array
    REAL :: dw(bd%is:bd%ie, bd%js:bd%je)
!---
    REAL, DIMENSION(bd%is:bd%ie+1, bd%js:bd%je+1) :: ub, vb
!  work array
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
!  needs this for corner_comm
    REAL :: ke(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
! Vorticity
    REAL :: vort(bd%isd:bd%ied, bd%jsd:bd%jed)
! 1-D X-direction Fluxes
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
! 1-D Y-direction Fluxes
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: gx(bd%is:bd%ie+1, bd%js:bd%je)
! work Y-dir flux array
    REAL :: gy(bd%is:bd%ie, bd%js:bd%je+1)
    LOGICAL :: fill_c
    REAL :: dt2, dt4, dt5, dt6
    REAL :: damp, damp2, damp4, dd8, u2, v2, du2, dv2
    REAL :: u_lon
    INTEGER :: i, j, is2, ie1, js2, je1, n, nt, n2, iq
    REAL, DIMENSION(:, :), POINTER :: area, area_c, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsina
    REAL, DIMENSION(:, :), POINTER :: f0, rsin2, divg_u, divg_v
    REAL, DIMENSION(:, :), POINTER :: cosa, dx, dy, dxc, dyc, rdxa, rdya&
&   , rdx, rdy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: npx, npy
    LOGICAL :: nested
    REAL :: delp_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pt_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: vort_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: delpc_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: ptc_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: ke_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL :: vc_tj(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: uc_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: divg_d_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL :: ut_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt_tj(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    REAL :: abs0
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4

    ut = 0.0
    vt = 0.0
    fx2 = 0.0
    fy2 = 0.0
    dw = 0.0
    ub = 0.0
    vb = 0.0
    wk = 0.0
    ke = 0.0
    vort = 0.0
    fx = 0.0
    fy = 0.0
    ra_x = 0.0
    ra_y = 0.0
    gx = 0.0
    gy = 0.0
    dt2 = 0.0
    dt4 = 0.0
    dt5 = 0.0
    dt6 = 0.0
    damp = 0.0
    damp2 = 0.0
    damp4 = 0.0
    dd8 = 0.0
    u2 = 0.0
    v2 = 0.0
    du2 = 0.0
    dv2 = 0.0
    u_lon = 0.0
    is2 = 0
    ie1 = 0
    js2 = 0
    je1 = 0
    n = 0
    nt = 0
    n2 = 0
    iq = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    npx = 0
    npy = 0
    max1 = 0
    max2 = 0
    max3 = 0
    max4 = 0
    abs0 = 0.0
    min1 = 0
    min2 = 0
    min3 = 0
    min4 = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    cosa_s => gridstruct%cosa_s
    rsin_u => gridstruct%rsin_u
    rsin_v => gridstruct%rsin_v
    rsina => gridstruct%rsina
    f0 => gridstruct%f0
    rsin2 => gridstruct%rsin2
    cosa => gridstruct%cosa
    dx => gridstruct%dx
    dy => gridstruct%dy
    rdxa => gridstruct%rdxa
    rdya => gridstruct%rdya
    rdx => gridstruct%rdx
    rdy => gridstruct%rdy
    sw_corner = gridstruct%sw_corner
    se_corner = gridstruct%se_corner
    nw_corner = gridstruct%nw_corner
    ne_corner = gridstruct%ne_corner
! end grid_type choices
    IF (flagstruct%grid_type .LT. 3) THEN
!!! TO DO: separate versions for nesting and for cubed-sphere
      IF (nested) THEN
        DO j=jsd,jed
          DO i=is-1,ie+2
            ut(i, j) = (uc(i, j)-0.25*cosa_u(i, j)*(vc(i-1, j)+vc(i, j)+&
&             vc(i-1, j+1)+vc(i, j+1)))*rsin_u(i, j)
          END DO
        END DO
        DO j=js-1,je+2
          DO i=isd,ied
            vt(i, j) = (vc(i, j)-0.25*cosa_v(i, j)*(uc(i, j-1)+uc(i+1, j&
&             -1)+uc(i, j)+uc(i+1, j)))*rsin_v(i, j)
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        DO j=jsd,jed
          IF (j .NE. 0 .AND. j .NE. 1 .AND. j .NE. npy - 1 .AND. j .NE. &
&             npy) THEN
            DO i=is-1,ie+2
              ut(i, j) = (uc(i, j)-0.25*cosa_u(i, j)*(vc(i-1, j)+vc(i, j&
&               )+vc(i-1, j+1)+vc(i, j+1)))*rsin_u(i, j)
            END DO
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        DO j=js-1,je+2
          IF (j .NE. 1 .AND. j .NE. npy) THEN
            DO i=isd,ied
              vt(i, j) = (vc(i, j)-0.25*cosa_v(i, j)*(uc(i, j-1)+uc(i+1&
&               , j-1)+uc(i, j)+uc(i+1, j)))*rsin_v(i, j)
            END DO
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        CALL PUSHCONTROL(1,1)
      END IF
!.not. nested
      IF (.NOT.nested) THEN
! West face
! West edge:
        IF (is .EQ. 1) THEN
          DO j=jsd,jed
            IF (uc(1, j)*dt .GT. 0.) THEN
              ut(1, j) = uc(1, j)/sin_sg(0, j, 3)
              CALL PUSHCONTROL(1,1)
            ELSE
              ut(1, j) = uc(1, j)/sin_sg(1, j, 1)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          IF (3 .LT. js) THEN
            max1 = js
          ELSE
            max1 = 3
          END IF
          IF (npy - 2 .GT. je + 1) THEN
            min1 = je + 1
          ELSE
            min1 = npy - 2
          END IF
          DO j=max1,min1
            vt(0, j) = vc(0, j) - 0.25*cosa_v(0, j)*(ut(0, j-1)+ut(1, j-&
&             1)+ut(0, j)+ut(1, j))
            vt(1, j) = vc(1, j) - 0.25*cosa_v(1, j)*(ut(1, j-1)+ut(2, j-&
&             1)+ut(1, j)+ut(2, j))
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! East edge:
        IF (ie + 1 .EQ. npx) THEN
          DO j=jsd,jed
            IF (uc(npx, j)*dt .GT. 0.) THEN
              ut(npx, j) = uc(npx, j)/sin_sg(npx-1, j, 3)
              CALL PUSHCONTROL(1,1)
            ELSE
              ut(npx, j) = uc(npx, j)/sin_sg(npx, j, 1)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          IF (3 .LT. js) THEN
            max2 = js
          ELSE
            max2 = 3
          END IF
          IF (npy - 2 .GT. je + 1) THEN
            min2 = je + 1
          ELSE
            min2 = npy - 2
          END IF
          DO j=max2,min2
            vt(npx-1, j) = vc(npx-1, j) - 0.25*cosa_v(npx-1, j)*(ut(npx-&
&             1, j-1)+ut(npx, j-1)+ut(npx-1, j)+ut(npx, j))
            vt(npx, j) = vc(npx, j) - 0.25*cosa_v(npx, j)*(ut(npx, j-1)+&
&             ut(npx+1, j-1)+ut(npx, j)+ut(npx+1, j))
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! South (Bottom) edge:
        IF (js .EQ. 1) THEN
          DO i=isd,ied
            IF (vc(i, 1)*dt .GT. 0.) THEN
              vt(i, 1) = vc(i, 1)/sin_sg(i, 0, 4)
              CALL PUSHCONTROL(1,1)
            ELSE
              vt(i, 1) = vc(i, 1)/sin_sg(i, 1, 2)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          IF (3 .LT. is) THEN
            max3 = is
          ELSE
            max3 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min3 = ie + 1
          ELSE
            min3 = npx - 2
          END IF
          DO i=max3,min3
            ut(i, 0) = uc(i, 0) - 0.25*cosa_u(i, 0)*(vt(i-1, 0)+vt(i, 0)&
&             +vt(i-1, 1)+vt(i, 1))
            ut(i, 1) = uc(i, 1) - 0.25*cosa_u(i, 1)*(vt(i-1, 1)+vt(i, 1)&
&             +vt(i-1, 2)+vt(i, 2))
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! North edge:
        IF (je + 1 .EQ. npy) THEN
          DO i=isd,ied
            IF (vc(i, npy)*dt .GT. 0.) THEN
              vt(i, npy) = vc(i, npy)/sin_sg(i, npy-1, 4)
              CALL PUSHCONTROL(1,1)
            ELSE
              vt(i, npy) = vc(i, npy)/sin_sg(i, npy, 2)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          IF (3 .LT. is) THEN
            max4 = is
          ELSE
            max4 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min4 = ie + 1
          ELSE
            min4 = npx - 2
          END IF
          DO i=max4,min4
            ut(i, npy-1) = uc(i, npy-1) - 0.25*cosa_u(i, npy-1)*(vt(i-1&
&             , npy-1)+vt(i, npy-1)+vt(i-1, npy)+vt(i, npy))
            ut(i, npy) = uc(i, npy) - 0.25*cosa_u(i, npy)*(vt(i-1, npy)+&
&             vt(i, npy)+vt(i-1, npy+1)+vt(i, npy+1))
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
! The following code solves a 2x2 system to get the interior parallel-to-edge uc,vc values
! near the corners (ex: for the sw corner ut(2,1) and vt(1,2) are solved for simultaneously).
! It then computes the halo uc, vc values so as to be consistent with the computations on
! the facing panel.
!The system solved is:
!  ut(2,1) = uc(2,1) - avg(vt)*cosa_u(2,1)
!  vt(1,2) = vc(1,2) - avg(ut)*cosa_v(1,2)
! in which avg(vt) includes vt(1,2) and avg(ut) includes ut(2,1)
        IF (sw_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(2, 0)*cosa_v(1, 0))
          ut(2, 0) = (uc(2, 0)-0.25*cosa_u(2, 0)*(vt(1, 1)+vt(2, 1)+vt(2&
&           , 0)+vc(1, 0)-0.25*cosa_v(1, 0)*(ut(1, 0)+ut(1, -1)+ut(2, -1&
&           ))))*damp
          damp = 1./(1.-0.0625*cosa_u(0, 1)*cosa_v(0, 2))
          vt(0, 2) = (vc(0, 2)-0.25*cosa_v(0, 2)*(ut(1, 1)+ut(1, 2)+ut(0&
&           , 2)+uc(0, 1)-0.25*cosa_u(0, 1)*(vt(0, 1)+vt(-1, 1)+vt(-1, 2&
&           ))))*damp
          damp = 1./(1.-0.0625*cosa_u(2, 1)*cosa_v(1, 2))
          ut(2, 1) = (uc(2, 1)-0.25*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(2&
&           , 2)+vc(1, 2)-0.25*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(2, 2))&
&           ))*damp
          vt(1, 2) = (vc(1, 2)-0.25*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(2&
&           , 2)+uc(2, 1)-0.25*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(2, 2))&
&           ))*damp
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (se_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(npx-1, 0)*cosa_v(npx-1, 0))
          ut(npx-1, 0) = (uc(npx-1, 0)-0.25*cosa_u(npx-1, 0)*(vt(npx-1, &
&           1)+vt(npx-2, 1)+vt(npx-2, 0)+vc(npx-1, 0)-0.25*cosa_v(npx-1&
&           , 0)*(ut(npx, 0)+ut(npx, -1)+ut(npx-1, -1))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx+1, 1)*cosa_v(npx, 2))
          vt(npx, 2) = (vc(npx, 2)-0.25*cosa_v(npx, 2)*(ut(npx, 1)+ut(&
&           npx, 2)+ut(npx+1, 2)+uc(npx+1, 1)-0.25*cosa_u(npx+1, 1)*(vt(&
&           npx, 1)+vt(npx+1, 1)+vt(npx+1, 2))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx-1, 1)*cosa_v(npx-1, 2))
          ut(npx-1, 1) = (uc(npx-1, 1)-0.25*cosa_u(npx-1, 1)*(vt(npx-1, &
&           1)+vt(npx-2, 1)+vt(npx-2, 2)+vc(npx-1, 2)-0.25*cosa_v(npx-1&
&           , 2)*(ut(npx, 1)+ut(npx, 2)+ut(npx-1, 2))))*damp
          vt(npx-1, 2) = (vc(npx-1, 2)-0.25*cosa_v(npx-1, 2)*(ut(npx, 1)&
&           +ut(npx, 2)+ut(npx-1, 2)+uc(npx-1, 1)-0.25*cosa_u(npx-1, 1)*&
&           (vt(npx-1, 1)+vt(npx-2, 1)+vt(npx-2, 2))))*damp
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (ne_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(npx-1, npy)*cosa_v(npx-1, npy+1))
          ut(npx-1, npy) = (uc(npx-1, npy)-0.25*cosa_u(npx-1, npy)*(vt(&
&           npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy+1)+vc(npx-1, npy+1)&
&           -0.25*cosa_v(npx-1, npy+1)*(ut(npx, npy)+ut(npx, npy+1)+ut(&
&           npx-1, npy+1))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx+1, npy-1)*cosa_v(npx, npy-1))
          vt(npx, npy-1) = (vc(npx, npy-1)-0.25*cosa_v(npx, npy-1)*(ut(&
&           npx, npy-1)+ut(npx, npy-2)+ut(npx+1, npy-2)+uc(npx+1, npy-1)&
&           -0.25*cosa_u(npx+1, npy-1)*(vt(npx, npy)+vt(npx+1, npy)+vt(&
&           npx+1, npy-1))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx-1, npy-1)*cosa_v(npx-1, npy-1)&
&           )
          ut(npx-1, npy-1) = (uc(npx-1, npy-1)-0.25*cosa_u(npx-1, npy-1)&
&           *(vt(npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy-1)+vc(npx-1, &
&           npy-1)-0.25*cosa_v(npx-1, npy-1)*(ut(npx, npy-1)+ut(npx, npy&
&           -2)+ut(npx-1, npy-2))))*damp
          vt(npx-1, npy-1) = (vc(npx-1, npy-1)-0.25*cosa_v(npx-1, npy-1)&
&           *(ut(npx, npy-1)+ut(npx, npy-2)+ut(npx-1, npy-2)+uc(npx-1, &
&           npy-1)-0.25*cosa_u(npx-1, npy-1)*(vt(npx-1, npy)+vt(npx-2, &
&           npy)+vt(npx-2, npy-1))))*damp
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (nw_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(2, npy)*cosa_v(1, npy+1))
          ut(2, npy) = (uc(2, npy)-0.25*cosa_u(2, npy)*(vt(1, npy)+vt(2&
&           , npy)+vt(2, npy+1)+vc(1, npy+1)-0.25*cosa_v(1, npy+1)*(ut(1&
&           , npy)+ut(1, npy+1)+ut(2, npy+1))))*damp
          damp = 1./(1.-0.0625*cosa_u(0, npy-1)*cosa_v(0, npy-1))
          vt(0, npy-1) = (vc(0, npy-1)-0.25*cosa_v(0, npy-1)*(ut(1, npy-&
&           1)+ut(1, npy-2)+ut(0, npy-2)+uc(0, npy-1)-0.25*cosa_u(0, npy&
&           -1)*(vt(0, npy)+vt(-1, npy)+vt(-1, npy-1))))*damp
          damp = 1./(1.-0.0625*cosa_u(2, npy-1)*cosa_v(1, npy-1))
          ut(2, npy-1) = (uc(2, npy-1)-0.25*cosa_u(2, npy-1)*(vt(1, npy)&
&           +vt(2, npy)+vt(2, npy-1)+vc(1, npy-1)-0.25*cosa_v(1, npy-1)*&
&           (ut(1, npy-1)+ut(1, npy-2)+ut(2, npy-2))))*damp
          vt(1, npy-1) = (vc(1, npy-1)-0.25*cosa_v(1, npy-1)*(ut(1, npy-&
&           1)+ut(1, npy-2)+ut(2, npy-2)+uc(2, npy-1)-0.25*cosa_u(2, npy&
&           -1)*(vt(1, npy)+vt(2, npy)+vt(2, npy-1))))*damp
          CALL PUSHCONTROL(2,3)
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
      ELSE
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
! flagstruct%grid_type >= 3
      DO j=jsd,jed
        DO i=is,ie+1
          ut(i, j) = uc(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          vt(i, j) = vc(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(2,0)
    END IF
    DO j=jsd,jed
      DO i=is,ie+1
        CALL PUSHREALARRAY(xfx_adv(i, j))
        xfx_adv(i, j) = dt*ut(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        CALL PUSHREALARRAY(yfx_adv(i, j))
        yfx_adv(i, j) = dt*vt(i, j)
      END DO
    END DO
! Explanation of the following code:
!    xfx_adv = dt*ut*dy
!    crx_adv = dt*ut/dx
    DO j=jsd,jed
!DEC$ VECTOR ALWAYS
      DO i=is,ie+1
        IF (xfx_adv(i, j) .GT. 0.) THEN
          CALL PUSHREALARRAY(crx_adv(i, j))
          crx_adv(i, j) = xfx_adv(i, j)*rdxa(i-1, j)
          CALL PUSHREALARRAY(xfx_adv(i, j))
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sin_sg(i-1, j, 3)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHREALARRAY(crx_adv(i, j))
          crx_adv(i, j) = xfx_adv(i, j)*rdxa(i, j)
          CALL PUSHREALARRAY(xfx_adv(i, j))
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sin_sg(i, j, 1)
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
    END DO
    DO j=js,je+1
!DEC$ VECTOR ALWAYS
      DO i=isd,ied
        IF (yfx_adv(i, j) .GT. 0.) THEN
          CALL PUSHREALARRAY(cry_adv(i, j))
          cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j-1)
          CALL PUSHREALARRAY(yfx_adv(i, j))
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sin_sg(i, j-1, 4)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHREALARRAY(cry_adv(i, j))
          cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j)
          CALL PUSHREALARRAY(yfx_adv(i, j))
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sin_sg(i, j, 2)
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
    END DO
    DO j=jsd,jed
      DO i=is,ie
        ra_x(i, j) = area(i, j) + (xfx_adv(i, j)-xfx_adv(i+1, j))
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        ra_y(i, j) = area(i, j) + (yfx_adv(i, j)-yfx_adv(i, j+1))
      END DO
    END DO
    IF (hord_dp .EQ. hord_dp_pert .AND. (.NOT.split_damp)) THEN
      CALL FV_TP_2D_FWD(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx&
&                    , fy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y&
&                    , nord=nord_v, damp_c=damp_v)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(delp, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL FV_TP_2D(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx, &
&             fy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, nord=&
&             nord_v, damp_c=damp_v)
      CALL PUSHCONTROL(1,0)
    END IF
! <<< Save the mass fluxes to the "Flux Capacitor" for tracer transport >>>
    DO j=jsd,jed
      DO i=is,ie+1
        cx(i, j) = cx(i, j) + crx_adv(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie+1
        xflux(i, j) = xflux(i, j) + fx(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        cy(i, j) = cy(i, j) + cry_adv(i, j)
      END DO
      DO i=is,ie
        yflux(i, j) = yflux(i, j) + fy(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie
        heat_source(i, j) = 0.
      END DO
    END DO
    IF (.NOT.hydrostatic) THEN
      IF (damp_w .GT. 1.e-5) THEN
        IF (dt .GE. 0.) THEN
          abs0 = dt
        ELSE
          abs0 = -dt
        END IF
        dd8 = kgb*abs0
        damp4 = (damp_w*gridstruct%da_min_c)**(nord_w+1)
        CALL DEL6_VT_FLUX(nord_w, npx, npy, damp4, w, wk, fx2, fy2, &
&                   gridstruct, bd)
        DO j=js,je
          DO i=is,ie
            dw(i, j) = (fx2(i, j)-fx2(i+1, j)+(fy2(i, j)-fy2(i, j+1)))*&
&             rarea(i, j)
! 0.5 * [ (w+dw)**2 - w**2 ] = w*dw + 0.5*dw*dw
!                   heat_source(i,j) = -d_con*dw(i,j)*(w(i,j)+0.5*dw(i,j))
            heat_source(i, j) = dd8 - dw(i, j)*(w(i, j)+0.5*dw(i, j))
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (hord_vt .EQ. hord_vt_pert) THEN
        CALL FV_TP_2D_FWD(w, crx_adv, cry_adv, npx, npy, hord_vt, gx&
&                      , gy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, &
&                      ra_y, mfx=fx, mfy=fy)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
        CALL PUSHREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
        CALL PUSHREALARRAY(w, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL FV_TP_2D(w, crx_adv, cry_adv, npx, npy, hord_vt, gx, &
&               gy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, mfx=fx&
&               , mfy=fy)
        CALL PUSHCONTROL(1,0)
      END IF
      DO j=js,je
        DO i=is,ie
          CALL PUSHREALARRAY(w(i, j))
          w(i, j) = delp(i, j)*w(i, j) + (gx(i, j)-gx(i+1, j)+(gy(i, j)-&
&           gy(i, j+1)))*rarea(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
!    if ( inline_q .and. zvir>0.01 ) then
!       do j=jsd,jed
!          do i=isd,ied
!             pt(i,j) = pt(i,j)/(1.+zvir*q(i,j,k,sphum))
!          enddo
!       enddo
!    endif
    IF (hord_tm .EQ. hord_tm_pert .AND. (.NOT.split_damp)) THEN
      CALL FV_TP_2D_FWD(pt, crx_adv, cry_adv, npx, npy, hord_tm, gx, &
&                    gy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, &
&                    mfx=fx, mfy=fy, mass=delp, nord=nord_t, damp_c=&
&                    damp_t)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL PUSHREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(pt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL FV_TP_2D(pt, crx_adv, cry_adv, npx, npy, hord_tm, gx, gy&
&             , xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, mfx=fx, &
&             mfy=fy, mass=delp, nord=nord_t, damp_c=damp_t)
      CALL PUSHCONTROL(1,1)
    END IF
    IF (inline_q) THEN
      DO j=js,je
        DO i=is,ie
          wk(i, j) = delp(i, j)
          CALL PUSHREALARRAY(delp(i, j))
          delp(i, j) = wk(i, j) + (fx(i, j)-fx(i+1, j)+(fy(i, j)-fy(i, j&
&           +1)))*rarea(i, j)
          CALL PUSHREALARRAY(pt(i, j))
          pt(i, j) = (pt(i, j)*wk(i, j)+(gx(i, j)-gx(i+1, j)+(gy(i, j)-&
&           gy(i, j+1)))*rarea(i, j))/delp(i, j)
        END DO
      END DO
      DO iq=1,nq
        IF (hord_tr .EQ. hord_tr_pert) THEN
          CALL FV_TP_2D_FWD(q(isd:ied, jsd:jed, k, iq), crx_adv, &
&                        cry_adv, npx, npy, hord_tr, gx, gy, xfx_adv, &
&                        yfx_adv, gridstruct, bd, ra_x, ra_y, mfx=fx, &
&                        mfy=fy, mass=delp, nord=nord_t, damp_c=damp_t)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
          CALL PUSHREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
          CALL PUSHREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1)*(&
&                       jed-jsd+1))
          CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), crx_adv, cry_adv, &
&                 npx, npy, hord_tr, gx, gy, xfx_adv, yfx_adv, &
&                 gridstruct, bd, ra_x, ra_y, mfx=fx, mfy=fy, mass=delp&
&                 , nord=nord_t, damp_c=damp_t)
          CALL PUSHCONTROL(1,0)
        END IF
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(q(i, j, k, iq))
            q(i, j, k, iq) = (q(i, j, k, iq)*wk(i, j)+(gx(i, j)-gx(i+1, &
&             j)+(gy(i, j)-gy(i, j+1)))*rarea(i, j))/delp(i, j)
          END DO
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
!     if ( zvir>0.01 ) then
!       do j=js,je
!          do i=is,ie
!             pt(i,j) = pt(i,j)*(1.+zvir*q(i,j,k,sphum))
!          enddo
!       enddo
!     endif
      DO j=js,je
        DO i=is,ie
          CALL PUSHREALARRAY(pt(i, j))
          pt(i, j) = pt(i, j)*delp(i, j) + (gx(i, j)-gx(i+1, j)+(gy(i, j&
&           )-gy(i, j+1)))*rarea(i, j)
          CALL PUSHREALARRAY(delp(i, j))
          delp(i, j) = delp(i, j) + (fx(i, j)-fx(i+1, j)+(fy(i, j)-fy(i&
&           , j+1)))*rarea(i, j)
          CALL PUSHREALARRAY(pt(i, j))
          pt(i, j) = pt(i, j)/delp(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    END IF
    IF (fpp%fpp_overload_r4) THEN
      DO j=js,je
        DO i=is,ie
          dpx(i, j) = dpx(i, j) + (fx(i, j)-fx(i+1, j)+(fy(i, j)-fy(i, j&
&           +1)))*rarea(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
!----------------------
! Kinetic Energy Fluxes
!----------------------
! Compute B grid contra-variant components for KE:
    dt5 = 0.5*dt
    dt4 = 0.25*dt
    IF (nested) THEN
      CALL PUSHCONTROL(1,0)
      is2 = is
      ie1 = ie + 1
      js2 = js
      je1 = je + 1
    ELSE
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
      IF (2 .LT. js) THEN
        js2 = js
      ELSE
        js2 = 2
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        CALL PUSHCONTROL(1,1)
        je1 = je + 1
      ELSE
        CALL PUSHCONTROL(1,1)
        je1 = npy - 1
      END IF
    END IF
!!! TO DO: separate versions for nested and for cubed-sphere
    IF (flagstruct%grid_type .LT. 3) THEN
      IF (nested) THEN
        DO j=js2,je1
          DO i=is2,ie1
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j)-(uc(i, j-1)+uc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
        END DO
        CALL PUSHCONTROL(2,0)
      ELSE
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb(i, 1) = dt5*(vt(i-1, 1)+vt(i, 1))
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        DO j=js2,je1
          DO i=is2,ie1
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j)-(uc(i, j-1)+uc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
          IF (is .EQ. 1) THEN
! 2-pt extrapolation from both sides:
            vb(1, j) = dt4*(-vt(-1, j)+3.*(vt(0, j)+vt(1, j))-vt(2, j))
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
! 2-pt extrapolation from both sides:
            vb(npx, j) = dt4*(-vt(npx-2, j)+3.*(vt(npx-1, j)+vt(npx, j))&
&             -vt(npx+1, j))
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb(i, npy) = dt5*(vt(i-1, npy)+vt(i, npy))
          END DO
          CALL PUSHCONTROL(2,1)
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          vb(i, j) = dt5*(vc(i-1, j)+vc(i, j))
        END DO
      END DO
      CALL PUSHCONTROL(2,3)
    END IF
    IF (hord_mt .EQ. hord_mt_pert) THEN
      CALL YTP_V_FWD(is, ie, js, je, isd, ied, jsd, jed, vb, u, v, ub&
&                 , hord_mt, gridstruct%dy, gridstruct%rdy, npx, npy, &
&                 flagstruct%grid_type, nested)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL YTP_V(is, ie, js, je, isd, ied, jsd, jed, vb, u, v, ub, &
&          hord_mt, gridstruct%dy, gridstruct%rdy, npx, npy, &
&          flagstruct%grid_type, nested)
      CALL PUSHCONTROL(1,0)
    END IF
    DO j=js,je+1
      DO i=is,ie+1
        ke(i, j) = vb(i, j)*ub(i, j)
      END DO
    END DO
    IF (flagstruct%grid_type .LT. 3) THEN
      IF (nested) THEN
        DO j=js,je+1
          DO i=is2,ie1
            CALL PUSHREALARRAY(ub(i, j))
            ub(i, j) = dt5*(uc(i, j-1)+uc(i, j)-(vc(i-1, j)+vc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
        END DO
        CALL PUSHCONTROL(2,0)
      ELSE
        IF (is .EQ. 1) THEN
          DO j=js,je+1
! corner values are incorrect
            CALL PUSHREALARRAY(ub(1, j))
            ub(1, j) = dt5*(ut(1, j-1)+ut(1, j))
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is2,ie1
! 2-pt extrapolation from both sides:
              CALL PUSHREALARRAY(ub(i, j))
              ub(i, j) = dt4*(-ut(i, j-2)+3.*(ut(i, j-1)+ut(i, j))-ut(i&
&               , j+1))
            END DO
            CALL PUSHCONTROL(1,1)
          ELSE
            DO i=is2,ie1
              CALL PUSHREALARRAY(ub(i, j))
              ub(i, j) = dt5*(uc(i, j-1)+uc(i, j)-(vc(i-1, j)+vc(i, j))*&
&               cosa(i, j))*rsina(i, j)
            END DO
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        IF (ie + 1 .EQ. npx) THEN
          DO j=js,je+1
! corner values are incorrect
            CALL PUSHREALARRAY(ub(npx, j))
            ub(npx, j) = dt5*(ut(npx, j-1)+ut(npx, j))
          END DO
          CALL PUSHCONTROL(2,1)
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY(ub(i, j))
          ub(i, j) = dt5*(uc(i, j-1)+uc(i, j))
        END DO
      END DO
      CALL PUSHCONTROL(2,3)
    END IF
    IF (hord_mt .EQ. hord_mt_pert) THEN
      CALL XTP_U_FWD(is, ie, js, je, isd, ied, jsd, jed, ub, u, v, vb&
&                 , hord_mt, gridstruct%dx, gridstruct%rdx, npx, npy, &
&                 flagstruct%grid_type, nested)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHREALARRAY(vb, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      CALL XTP_U(is, ie, js, je, isd, ied, jsd, jed, ub, u, v, vb, &
&          hord_mt, gridstruct%dx, gridstruct%rdx, npx, npy, &
&          flagstruct%grid_type, nested)
      CALL PUSHCONTROL(1,0)
    END IF
    DO j=js,je+1
      DO i=is,ie+1
        ke(i, j) = 0.5*(ke(i, j)+ub(i, j)*vb(i, j))
      END DO
    END DO
!-----------------------------------------
! Fix KE at the 4 corners of the face:
!-----------------------------------------
    IF (.NOT.nested) THEN
      dt6 = dt/6.
      IF (sw_corner) THEN
        ke(1, 1) = dt6*((ut(1, 1)+ut(1, 0))*u(1, 1)+(vt(1, 1)+vt(0, 1))*&
&         v(1, 1)+(ut(1, 1)+vt(1, 1))*u(0, 1))
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (se_corner) THEN
!i = npx
        ke(npx, 1) = dt6*((ut(npx, 1)+ut(npx, 0))*u(npx-1, 1)+(vt(npx, 1&
&         )+vt(npx-1, 1))*v(npx, 1)+(ut(npx, 1)-vt(npx-1, 1))*u(npx, 1))
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ne_corner) THEN
!i = npx;      j = npy
        ke(npx, npy) = dt6*((ut(npx, npy)+ut(npx, npy-1))*u(npx-1, npy)+&
&         (vt(npx, npy)+vt(npx-1, npy))*v(npx, npy-1)+(ut(npx, npy-1)+vt&
&         (npx-1, npy))*u(npx, npy))
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (nw_corner) THEN
!j = npy
        ke(1, npy) = dt6*((ut(1, npy)+ut(1, npy-1))*u(1, npy)+(vt(1, npy&
&         )+vt(0, npy))*v(1, npy-1)+(ut(1, npy-1)-vt(1, npy))*u(0, npy))
        CALL PUSHCONTROL(2,2)
      ELSE
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHCONTROL(2,0)
    END IF
! Compute vorticity:
    DO j=jsd,jed+1
      DO i=isd,ied
        CALL PUSHREALARRAY(vt(i, j))
        vt(i, j) = u(i, j)*dx(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied+1
        CALL PUSHREALARRAY(ut(i, j))
        ut(i, j) = v(i, j)*dy(i, j)
      END DO
    END DO
! wk is "volume-mean" relative vorticity
    DO j=jsd,jed
      DO i=isd,ied
        CALL PUSHREALARRAY(wk(i, j))
        wk(i, j) = rarea(i, j)*(vt(i, j)-vt(i, j+1)+(ut(i+1, j)-ut(i, j)&
&         ))
      END DO
    END DO
    IF (.NOT.hydrostatic) THEN
      IF (flagstruct%do_f3d) THEN
        CALL PUSHCONTROL(1,0)
      ELSE
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(w(i, j))
            w(i, j) = w(i, j)/delp(i, j)
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      END IF
      IF (damp_w .GT. 1.e-5) THEN
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(w(i, j))
            w(i, j) = w(i, j) + dw(i, j)
          END DO
        END DO
        CALL PUSHCONTROL(2,0)
      ELSE
        CALL PUSHCONTROL(2,1)
      END IF
    ELSE
      CALL PUSHCONTROL(2,2)
    END IF
!-----------------------------
! Compute divergence damping
!-----------------------------
!!  damp = dddmp * da_min_c
!
!   if ( nord==0 ) then
!!         area ~ dxb*dyb*sin(alpha)
!
!      if (nested) then
!
!         do j=js,je+1
!            do i=is-1,ie+1
!               ptc(i,j) = (u(i,j)-0.5*(va(i,j-1)+va(i,j))*cosa_v(i,j))   &
!                    *dyc(i,j)*sina_v(i,j)
!            enddo
!         enddo
!
!         do j=js-1,je+1
!            do i=is2,ie1
!               vort(i,j) = (v(i,j) - 0.5*(ua(i-1,j)+ua(i,j))*cosa_u(i,j))  &
!                    *dxc(i,j)*sina_u(i,j)
!            enddo
!         enddo
!
!      else
!         do j=js,je+1
!
!            if ( (j==1 .or. j==npy)  ) then
!               do i=is-1,ie+1
!                  if (vc(i,j) > 0) then
!                     ptc(i,j) = u(i,j)*dyc(i,j)*sin_sg(i,j-1,4)
!                  else
!                     ptc(i,j) = u(i,j)*dyc(i,j)*sin_sg(i,j,2)
!                  end if
!               enddo
!            else
!               do i=is-1,ie+1
!                  ptc(i,j) = (u(i,j)-0.5*(va(i,j-1)+va(i,j))*cosa_v(i,j))   &
!                       *dyc(i,j)*sina_v(i,j)
!               enddo
!            endif
!         enddo
!
!         do j=js-1,je+1
!            do i=is2,ie1
!               vort(i,j) = (v(i,j) - 0.5*(ua(i-1,j)+ua(i,j))*cosa_u(i,j))  &
!                    *dxc(i,j)*sina_u(i,j)
!            enddo
!            if ( is ==  1 ) then
!               if (uc(1,j) > 0) then
!                  vort(1,  j) = v(1,  j)*dxc(1,  j)*sin_sg(0,j,3)
!               else
!                  vort(1,  j) = v(1,  j)*dxc(1,  j)*sin_sg(1,j,1)
!               end if
!            end if
!            if ( (ie+1)==npx ) then 
!               if (uc(npx,j) > 0) then
!                  vort(npx,j) = v(npx,j)*dxc(npx,j)* & 
!                       sin_sg(npx-1,j,3)
!               else
!                  vort(npx,j) = v(npx,j)*dxc(npx,j)* &
!                       sin_sg(npx,j,1)
!               end if
!            end if
!         enddo
!      endif
!
!      do j=js,je+1
!         do i=is,ie+1
!            delpc(i,j) = vort(i,j-1) - vort(i,j) + ptc(i-1,j) - ptc(i,j)
!         enddo
!      enddo
!
!! Remove the extra term at the corners:
!      if (sw_corner) delpc(1,    1) = delpc(1,    1) - vort(1,    0)
!      if (se_corner) delpc(npx,  1) = delpc(npx,  1) - vort(npx,  0)
!      if (ne_corner) delpc(npx,npy) = delpc(npx,npy) + vort(npx,npy)
!      if (nw_corner) delpc(1,  npy) = delpc(1,  npy) + vort(1,  npy)
!
!      do j=js,je+1
!         do i=is,ie+1
!            delpc(i,j) = gridstruct%rarea_c(i,j)*delpc(i,j)
!                damp = gridstruct%da_min_c*max(d2_bg, min(0.20, dddmp*abs(delpc(i,j)*dt)))
!                vort(i,j) = damp*delpc(i,j)
!                ke(i,j) = ke(i,j) + vort(i,j)
!         enddo
!      enddo
!   else
!!--------------------------
!! Higher order divg damping
!!--------------------------
!     do j=js,je+1
!        do i=is,ie+1
!! Save divergence for external mode filter
!           delpc(i,j) = divg_d(i,j)
!        enddo
!     enddo
!
!     n2 = nord + 1    ! N > 1
!     do n=1,nord
!        nt = nord-n
!
!        fill_c = (nt/=0) .and. (flagstruct%grid_type<3) .and.               &
!                 ( sw_corner .or. se_corner .or. ne_corner .or. nw_corner ) &
!                  .and. .not. nested
!
!        if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=XDir, BGRID=.true.)
!        do j=js-nt,je+1+nt
!           do i=is-1-nt,ie+1+nt
!              vc(i,j) = (divg_d(i+1,j)-divg_d(i,j))*divg_u(i,j)
!           enddo
!        enddo
!
!        if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=YDir, BGRID=.true.)
!        do j=js-1-nt,je+1+nt
!           do i=is-nt,ie+1+nt
!              uc(i,j) = (divg_d(i,j+1)-divg_d(i,j))*divg_v(i,j)
!           enddo
!        enddo
!
!        if ( fill_c ) call fill_corners(vc, uc, npx, npy, VECTOR=.true., DGRID=.true.)
!        do j=js-nt,je+1+nt
!           do i=is-nt,ie+1+nt
!              divg_d(i,j) = uc(i,j-1) - uc(i,j) + vc(i-1,j) - vc(i,j)
!           enddo
!        enddo
!
!! Remove the extra term at the corners:
!        if (sw_corner) divg_d(1,    1) = divg_d(1,    1) - uc(1,    0)
!        if (se_corner) divg_d(npx,  1) = divg_d(npx,  1) - uc(npx,  0)
!        if (ne_corner) divg_d(npx,npy) = divg_d(npx,npy) + uc(npx,npy)
!        if (nw_corner) divg_d(1,  npy) = divg_d(1,  npy) + uc(1,  npy)
!
!     if ( .not. gridstruct%stretched_grid ) then
!        do j=js-nt,je+1+nt
!           do i=is-nt,ie+1+nt
!              divg_d(i,j) = divg_d(i,j)*gridstruct%rarea_c(i,j)
!           enddo
!        enddo
!     endif
!
!     enddo ! n-loop
!
!     if ( dddmp<1.E-5) then
!          vort(:,:) = 0.
!     else
!      if ( flagstruct%grid_type < 3 ) then
!! Interpolate relative vort to cell corners
!          call a2b_ord4(wk, vort, gridstruct, npx, npy, is, ie, js, je, ng, .false.)
!          do j=js,je+1
!             do i=is,ie+1
!! The following is an approxi form of Smagorinsky diffusion
!                vort(i,j) = abs(dt)*sqrt(delpc(i,j)**2 + vort(i,j)**2)
!             enddo
!          enddo
!      else  ! Correct form: works only for doubly preiodic domain
!          call smag_corner(abs(dt), u, v, ua, va, vort, bd, npx, npy, gridstruct, ng)
!      endif
!     endif
!
!     if (gridstruct%stretched_grid ) then
!! Stretched grid with variable damping ~ area
!         dd8 = gridstruct%da_min * d4_bg**n2
!     else
!         dd8 = ( gridstruct%da_min_c*d4_bg )**n2
!     endif
!
!     do j=js,je+1
!        do i=is,ie+1
!           damp2 =  gridstruct%da_min_c*max(d2_bg, min(0.20, dddmp*vort(i,j)))  ! del-2
!           vort(i,j) = damp2*delpc(i,j) + dd8*divg_d(i,j)
!             ke(i,j) = ke(i,j) + vort(i,j)
!        enddo
!     enddo
!
!   endif
    IF (.NOT.split_damp) THEN
      CALL COMPUTE_DIVERGENCE_DAMPING_FWD(nord, d2_bg, d4_bg, dddmp, &
&                                      dt, vort, ptc, delpc, ke, u, v, &
&                                      uc, vc, ua, va, divg_d, wk, &
&                                      gridstruct, flagstruct, bd)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(divg_d, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+2))
      CALL PUSHREALARRAY(vc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      CALL PUSHREALARRAY(uc, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(delpc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(ptc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL COMPUTE_DIVERGENCE_DAMPING(nord, d2_bg, d4_bg&
&                               , dddmp, dt, vort, ptc, delpc, ke, &
&                               u, v, uc, vc, ua, va, divg_d, wk, &
&                               gridstruct, flagstruct, bd)
      CALL PUSHCONTROL(1,1)
    END IF
    IF (d_con .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(ub(i, j))
          ub(i, j) = vort(i, j) - vort(i+1, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(vb(i, j))
          vb(i, j) = vort(i, j) - vort(i, j+1)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
! Vorticity transport
    IF (hydrostatic) THEN
      DO j=jsd,jed
        DO i=isd,ied
          vort(i, j) = wk(i, j) + f0(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(2,0)
    ELSE IF (flagstruct%do_f3d) THEN
      DO j=jsd,jed
        DO i=isd,ied
          vort(i, j) = wk(i, j) + f0(i, j)*z_rat(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(2,1)
    ELSE
      DO j=jsd,jed
        DO i=isd,ied
          vort(i, j) = wk(i, j) + f0(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(2,2)
    END IF
    IF (hord_vt .EQ. hord_vt_pert) THEN
      CALL FV_TP_2D_FWD(vort, crx_adv, cry_adv, npx, npy, hord_vt, fx&
&                    , fy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL FV_TP_2D(vort, crx_adv, cry_adv, npx, npy, hord_vt, fx, &
&             fy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y)
      CALL PUSHCONTROL(1,0)
    END IF
    DO j=js,je+1
      DO i=is,ie
        CALL PUSHREALARRAY(u(i, j))
        u(i, j) = vt(i, j) + (ke(i, j)-ke(i+1, j)) + fy(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie+1
        CALL PUSHREALARRAY(v(i, j))
        v(i, j) = ut(i, j) + (ke(i, j)-ke(i, j+1)) - fx(i, j)
      END DO
    END DO
!--------------------------------------------------------
! damping applied to relative vorticity (wk):
    IF (damp_v .GT. 1.e-5) THEN
      damp4 = (damp_v*gridstruct%da_min_c)**(nord_v+1)
      CALL PUSHREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      CALL PUSHREALARRAY(ut, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL DEL6_VT_FLUX(nord_v, npx, npy, damp4, wk, vort, ut, vt, &
&                 gridstruct, bd)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (damp_v_pert .GT. 1.e-5) THEN
      damp4 = (damp_v_pert*gridstruct%da_min_c)**(nord_v_pert+1)
      CALL PUSHREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      CALL PUSHREALARRAY(ut, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
!      CALL DEL6_VT_FLUX(nord_v_pert, npx, npy, damp4, wk, vort, ut, vt, &
!&                 gridstruct, bd)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (d_con .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(ub(i, j))
          ub(i, j) = (ub(i, j)+vt(i, j))*rdx(i, j)
          CALL PUSHREALARRAY(fy(i, j))
          fy(i, j) = u(i, j)*rdx(i, j)
          CALL PUSHREALARRAY(gy(i, j))
          gy(i, j) = fy(i, j)*ub(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(vb(i, j))
          vb(i, j) = (vb(i, j)-ut(i, j))*rdy(i, j)
          CALL PUSHREALARRAY(fx(i, j))
          fx(i, j) = v(i, j)*rdy(i, j)
          CALL PUSHREALARRAY(gx(i, j))
          gx(i, j) = fx(i, j)*vb(i, j)
        END DO
      END DO
!----------------------------------
! Heating due to damping:
!----------------------------------
      damp = 0.25*d_con
      DO j=js,je
        DO i=is,ie
          u2 = fy(i, j) + fy(i, j+1)
          du2 = ub(i, j) + ub(i, j+1)
          v2 = fx(i, j) + fx(i+1, j)
          dv2 = vb(i, j) + vb(i+1, j)
! Total energy conserving:
! Convert lost KE due to divergence damping to "heat"
          CALL PUSHREALARRAY(heat_source(i, j))
          heat_source(i, j) = delp(i, j)*(heat_source(i, j)-damp*rsin2(i&
&           , j)*(ub(i, j)**2+ub(i, j+1)**2+vb(i, j)**2+vb(i+1, j)**2+2.&
&           *(gy(i, j)+gy(i, j+1)+gx(i, j)+gx(i+1, j))-cosa_s(i, j)*(u2*&
&           dv2+v2*du2+du2*dv2)))
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
! Add diffusive fluxes to the momentum equation:
    IF (damp_v .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(u(i, j))
          u(i, j) = u(i, j) + vt(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(v(i, j))
          v(i, j) = v(i, j) - ut(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (damp_v_pert .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(u(i, j))
          !u(i, j) = u(i, j) + vt(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(v(i, j))
          !v(i, j) = v(i, j) - ut(i, j)
        END DO
      END DO
      CALL PUSHREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      !CALL PUSHPOINTER8(C_LOC(rdxa))
      CALL PUSHINTEGER(je)
      CALL PUSHINTEGER(min4)
      CALL PUSHINTEGER(ie1)
      CALL PUSHINTEGER(min3)
      CALL PUSHINTEGER(min2)
      CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL PUSHINTEGER(min1)
      CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL PUSHINTEGER(js2)
      CALL PUSHREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      !CALL PUSHPOINTER8(C_LOC(rdx))
      !CALL PUSHPOINTER8(C_LOC(rsina))
      CALL PUSHINTEGER(is)
      !CALL PUSHPOINTER8(C_LOC(rsin_u))
      CALL PUSHREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(vb, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      CALL PUSHREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
      CALL PUSHINTEGER(ie)
      CALL PUSHREALARRAY(dt6)
      CALL PUSHREALARRAY(dt4)
      !CALL PUSHPOINTER8(C_LOC(sin_sg))
      CALL PUSHINTEGER(je1)
      CALL PUSHREALARRAY(ut, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      !CALL PUSHPOINTER8(C_LOC(cosa))
      !CALL PUSHPOINTER8(C_LOC(cosa_v))
      !CALL PUSHPOINTER8(C_LOC(cosa_u))
      CALL PUSHREALARRAY(ub, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      !CALL PUSHPOINTER8(C_LOC(rdya))
      !CALL PUSHPOINTER8(C_LOC(dy))
      CALL PUSHINTEGER(is2)
      !CALL PUSHPOINTER8(C_LOC(dx))
      CALL PUSHREALARRAY(dw, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL PUSHREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(damp)
      CALL PUSHINTEGER(max4)
      CALL PUSHINTEGER(max3)
      CALL PUSHINTEGER(max2)
      CALL PUSHINTEGER(max1)
      CALL PUSHINTEGER(npy)
      CALL PUSHINTEGER(npx)
      CALL PUSHINTEGER(js)
      CALL PUSHREALARRAY(damp4)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      !CALL PUSHPOINTER8(C_LOC(rdxa))
      CALL PUSHINTEGER(je)
      CALL PUSHINTEGER(min4)
      CALL PUSHINTEGER(ie1)
      CALL PUSHINTEGER(min3)
      CALL PUSHINTEGER(min2)
      CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL PUSHINTEGER(min1)
      CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL PUSHINTEGER(js2)
      CALL PUSHREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      !CALL PUSHPOINTER8(C_LOC(rdx))
      !CALL PUSHPOINTER8(C_LOC(rsina))
      CALL PUSHINTEGER(is)
      !CALL PUSHPOINTER8(C_LOC(rsin_u))
      CALL PUSHREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(vb, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      CALL PUSHREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
      CALL PUSHINTEGER(ie)
      CALL PUSHREALARRAY(dt6)
      CALL PUSHREALARRAY(dt4)
      !CALL PUSHPOINTER8(C_LOC(sin_sg))
      CALL PUSHINTEGER(je1)
      CALL PUSHREALARRAY(ut, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      !CALL PUSHPOINTER8(C_LOC(cosa))
      !CALL PUSHPOINTER8(C_LOC(cosa_v))
      !CALL PUSHPOINTER8(C_LOC(cosa_u))
      CALL PUSHREALARRAY(ub, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      !CALL PUSHPOINTER8(C_LOC(rdya))
      !CALL PUSHPOINTER8(C_LOC(dy))
      CALL PUSHINTEGER(is2)
      !CALL PUSHPOINTER8(C_LOC(dx))
      CALL PUSHREALARRAY(dw, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL PUSHREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(damp)
      CALL PUSHINTEGER(max4)
      CALL PUSHINTEGER(max3)
      CALL PUSHINTEGER(max2)
      CALL PUSHINTEGER(max1)
      CALL PUSHINTEGER(npy)
      CALL PUSHINTEGER(npx)
      CALL PUSHINTEGER(js)
      CALL PUSHREALARRAY(damp4)
      CALL PUSHCONTROL(1,0)
    END IF
  END SUBROUTINE D_SW_FWD
!  Differentiation of d_sw in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2
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
!   gradient     of useful results: yfx_adv q crx_adv u v w delp
!                ua xfx_adv uc ptc xflux cry_adv delpc va vc yflux
!                divg_d z_rat heat_source pt cx cy dpx
!   with respect to varying inputs: yfx_adv q crx_adv u v w delp
!                ua xfx_adv uc ptc xflux cry_adv delpc va vc yflux
!                divg_d z_rat heat_source pt cx cy dpx
!     d_sw :: D-Grid Shallow Water Routine
  SUBROUTINE D_SW_BWD(delpc, delpc_ad, delp, delp_ad, ptc, ptc_ad, pt, &
&   pt_ad, u, u_ad, v, v_ad, w, w_ad, uc, uc_ad, vc, vc_ad, ua, ua_ad, &
&   va, va_ad, divg_d, divg_d_ad, xflux, xflux_ad, yflux, yflux_ad, cx, &
&   cx_ad, cy, cy_ad, crx_adv, crx_adv_ad, cry_adv, cry_adv_ad, xfx_adv&
&   , xfx_adv_ad, yfx_adv, yfx_adv_ad, q_con, z_rat, z_rat_ad, kgb, &
&   heat_source, heat_source_ad, dpx, dpx_ad, zvir, sphum, nq, q, q_ad, &
&   k, km, inline_q, dt, hord_tr, hord_mt, hord_vt, hord_tm, hord_dp, &
&   nord, nord_v, nord_w, nord_t, dddmp, d2_bg, d4_bg, damp_v, damp_w, &
&   damp_t, d_con, hydrostatic, gridstruct, flagstruct, bd, hord_tr_pert&
&   , hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_dp_pert, split_damp&
&   , nord_pert, nord_v_pert, nord_w_pert, nord_t_pert, dddmp_pert, &
&   d2_bg_pert, d4_bg_pert, damp_v_pert, damp_w_pert, damp_t_pert)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: hord_tr, hord_mt, hord_vt, hord_tm, hord_dp
    INTEGER, INTENT(IN) :: nord
    INTEGER, INTENT(IN) :: nord_v
    INTEGER, INTENT(IN) :: nord_w
    INTEGER, INTENT(IN) :: nord_t
    INTEGER, INTENT(IN) :: sphum, nq, k, km
    REAL, INTENT(IN) :: dt, dddmp, d2_bg, d4_bg, d_con
    REAL, INTENT(IN) :: zvir
    REAL, INTENT(IN) :: damp_v, damp_w, damp_t, kgb
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: hord_tr_pert, hord_mt_pert, hord_vt_pert, &
&   hord_tm_pert, hord_dp_pert, nord_pert, nord_v_pert, nord_w_pert, &
&   nord_t_pert
    LOGICAL, INTENT(IN) :: split_damp
    REAL, INTENT(IN) :: dddmp_pert, d2_bg_pert, d4_bg_pert, damp_v_pert&
&   , damp_w_pert, damp_t_pert
    REAL, INTENT(INOUT) :: divg_d(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL, INTENT(INOUT) :: divg_d_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: z_rat
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: z_rat_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: delp&
&   , pt, ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delp_ad, pt_ad, ua_ad, va_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   q_con
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: u&
&   , vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   u_ad, vc_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: v&
&   , uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   v_ad, uc_ad
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, km, nq)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed, km, nq)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: delpc, ptc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: delpc_ad, ptc_ad
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: heat_source
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je) :: heat_source_ad
    REAL(kind=8), DIMENSION(bd%is:bd%ie, bd%js:bd%je), INTENT(INOUT) :: &
&   dpx
    REAL(kind=8), DIMENSION(bd%is:bd%ie, bd%js:bd%je), INTENT(INOUT) :: &
&   dpx_ad
    REAL, INTENT(INOUT) :: xflux(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: xflux_ad(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: yflux(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: yflux_ad(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: cx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: cy_ad(bd%isd:bd%ied, bd%js:bd%je+1)
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: inline_q
    REAL, DIMENSION(bd%is:bd%ie+1, bd%jsd:bd%jed) :: crx_adv, xfx_adv
    REAL, DIMENSION(bd%is:bd%ie+1, bd%jsd:bd%jed) :: crx_adv_ad, &
&   xfx_adv_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%js:bd%je+1) :: cry_adv, yfx_adv
    REAL, DIMENSION(bd%isd:bd%ied, bd%js:bd%je+1) :: cry_adv_ad, &
&   yfx_adv_ad
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL :: ut(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: ut_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: vt_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: fx2_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: fy2(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: fy2_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: dw(bd%is:bd%ie, bd%js:bd%je)
    REAL :: dw_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL, DIMENSION(bd%is:bd%ie+1, bd%js:bd%je+1) :: ub, vb
    REAL, DIMENSION(bd%is:bd%ie+1, bd%js:bd%je+1) :: ub_ad, vb_ad
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: ke(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL :: ke_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL :: vort(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: vort_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_ad(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_ad(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_ad(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_ad(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: gx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: gx_ad(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: gy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: gy_ad(bd%is:bd%ie, bd%js:bd%je+1)
    LOGICAL :: fill_c
    REAL :: dt2, dt4, dt5, dt6
    REAL :: damp, damp2, damp4, dd8, u2, v2, du2, dv2
    REAL :: u2_ad, v2_ad, du2_ad, dv2_ad
    REAL :: u_lon
    INTEGER :: i, j, is2, ie1, js2, je1, n, nt, n2, iq
    REAL, DIMENSION(:, :), POINTER :: area, area_c, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsina
    REAL, DIMENSION(:, :), POINTER :: f0, rsin2, divg_u, divg_v
    REAL, DIMENSION(:, :), POINTER :: cosa, dx, dy, dxc, dyc, rdxa, rdya&
&   , rdx, rdy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: npx, npy
    LOGICAL :: nested
    REAL :: delp_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pt_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: vort_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: delpc_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: ptc_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: ke_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL :: vc_tj(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: uc_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: divg_d_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL :: ut_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt_tj(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    REAL :: abs0
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
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
    REAL :: temp_ad26
    REAL :: temp_ad27
    REAL :: temp_ad28
    REAL :: temp_ad29
    REAL :: temp_ad30
    REAL :: temp_ad31
    REAL :: temp_ad32
    REAL :: temp_ad33
    REAL :: temp_ad34
    REAL :: temp_ad35
    REAL :: temp_ad36
    REAL :: temp_ad37
    REAL :: temp_ad38
    REAL :: temp_ad39
    REAL :: temp_ad40
    REAL :: temp_ad41
    REAL :: temp_ad42
    REAL :: temp_ad43
    REAL :: temp_ad44
    REAL :: temp_ad45
    REAL :: temp_ad46
    REAL :: temp_ad47
    REAL :: temp_ad48
    REAL :: temp_ad49
    REAL :: temp_ad50
    REAL :: temp_ad51
    REAL :: temp_ad52
    REAL :: temp
    REAL :: temp_ad53
    REAL :: temp_ad54
    REAL :: temp_ad55
    REAL :: temp_ad56
    REAL :: temp_ad57
    REAL :: temp_ad58
    REAL :: temp_ad59
    REAL :: temp_ad60
    REAL :: temp_ad61
    REAL :: temp_ad62
    REAL :: temp_ad63
    REAL :: temp_ad64
    REAL :: temp_ad65
    REAL :: temp_ad66
    REAL :: temp_ad67
    REAL :: temp_ad68
    REAL :: temp_ad69
    REAL :: temp_ad70
    REAL :: temp_ad71
    REAL :: temp_ad72
    REAL :: temp_ad73
    REAL :: temp_ad74
    REAL :: temp_ad75
    REAL :: temp_ad76
    REAL :: temp_ad77
    REAL :: temp_ad78
    REAL :: temp_ad79
    REAL :: temp_ad80
    REAL :: temp_ad81
    REAL :: temp_ad82
    REAL :: temp_ad83
    REAL :: temp_ad84
    REAL :: temp_ad85
    REAL :: temp_ad86
    REAL :: temp_ad87
    REAL :: temp0
    REAL :: temp_ad88
    REAL :: temp_ad89
    REAL :: temp_ad90
    REAL :: temp_ad91
    INTEGER :: branch
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_d_sw

    ut = 0.0
    vt = 0.0
    fx2 = 0.0
    fy2 = 0.0
    dw = 0.0
    ub = 0.0
    vb = 0.0
    wk = 0.0
    ke = 0.0
    vort = 0.0
    fx = 0.0
    fy = 0.0
    ra_x = 0.0
    ra_y = 0.0
    gx = 0.0
    gy = 0.0
    dt2 = 0.0
    dt4 = 0.0
    dt5 = 0.0
    dt6 = 0.0
    damp = 0.0
    damp2 = 0.0
    damp4 = 0.0
    dd8 = 0.0
    u2 = 0.0
    v2 = 0.0
    du2 = 0.0
    dv2 = 0.0
    u_lon = 0.0
    is2 = 0
    ie1 = 0
    js2 = 0
    je1 = 0
    n = 0
    nt = 0
    n2 = 0
    iq = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    npx = 0
    npy = 0
    max1 = 0
    max2 = 0
    max3 = 0
    max4 = 0
    abs0 = 0.0
    min1 = 0
    min2 = 0
    min3 = 0
    min4 = 0
    branch = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    cosa_s => gridstruct%cosa_s
    rsin_u => gridstruct%rsin_u
    rsin_v => gridstruct%rsin_v
    rsina => gridstruct%rsina
    f0 => gridstruct%f0
    rsin2 => gridstruct%rsin2
    cosa => gridstruct%cosa
    dx => gridstruct%dx
    dy => gridstruct%dy
    rdxa => gridstruct%rdxa
    rdya => gridstruct%rdya
    rdx => gridstruct%rdx
    rdy => gridstruct%rdy
    sw_corner = gridstruct%sw_corner
    se_corner = gridstruct%se_corner
    nw_corner = gridstruct%nw_corner
    ne_corner = gridstruct%ne_corner
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(damp4)
      CALL POPINTEGER(js)
      CALL POPINTEGER(npx)
      CALL POPINTEGER(npy)
      CALL POPINTEGER(max1)
      CALL POPINTEGER(max2)
      CALL POPINTEGER(max3)
      CALL POPINTEGER(max4)
      CALL POPREALARRAY(damp)
      CALL POPREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL POPREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL POPREALARRAY(dw, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
      !CALL POPPOINTER8(cptr)
      dx => gridstruct%dx ! (/unknown_shape_in_d_sw/))
      CALL POPINTEGER(is2)
      !CALL POPPOINTER8(cptr)
      dy => gridstruct%dy ! (/unknown_shape_in_d_sw/))
      !CALL POPPOINTER8(cptr)
      rdya => gridstruct%rdya ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(ub, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      !CALL POPPOINTER8(cptr)
      cosa_u => gridstruct%cosa_u ! (/unknown_shape_in_d_sw/))
      !CALL POPPOINTER8(cptr)
      cosa_v => gridstruct%cosa_v ! (/unknown_shape_in_d_sw/))
      !CALL POPPOINTER8(cptr)
      cosa => gridstruct%cosa ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(ut, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL POPINTEGER(je1)
      !CALL POPPOINTER8(cptr)
      sin_sg => gridstruct%sin_sg ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(dt4)
      CALL POPREALARRAY(dt6)
      CALL POPINTEGER(ie)
      CALL POPREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(vb, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      CALL POPREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
      CALL POPREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      !CALL POPPOINTER8(cptr)
      rsin_u => gridstruct%rsin_u ! (/unknown_shape_in_d_sw/))
      CALL POPINTEGER(is)
      !CALL POPPOINTER8(cptr)
      rsina => gridstruct%rsina ! (/unknown_shape_in_d_sw/))
      !CALL POPPOINTER8(cptr)
      rdx => gridstruct%rdx ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      CALL POPINTEGER(js2)
      CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL POPINTEGER(min1)
      CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL POPINTEGER(min2)
      CALL POPINTEGER(min3)
      CALL POPINTEGER(ie1)
      CALL POPINTEGER(min4)
      CALL POPINTEGER(je)
      !CALL POPPOINTER8(cptr)
      rdxa => gridstruct%rdxa ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      ut_ad = 0.0
      vt_ad = 0.0
    ELSE
      CALL POPREALARRAY(damp4)
      CALL POPINTEGER(js)
      CALL POPINTEGER(npx)
      CALL POPINTEGER(npy)
      CALL POPINTEGER(max1)
      CALL POPINTEGER(max2)
      CALL POPINTEGER(max3)
      CALL POPINTEGER(max4)
      CALL POPREALARRAY(damp)
      CALL POPREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL POPREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL POPREALARRAY(dw, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
      !CALL POPPOINTER8(cptr)
      dx => gridstruct%dx ! (/unknown_shape_in_d_sw/))
      CALL POPINTEGER(is2)
      !CALL POPPOINTER8(cptr)
      dy => gridstruct%dy ! (/unknown_shape_in_d_sw/))
      !CALL POPPOINTER8(cptr)
      rdya => gridstruct%rdya ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(ub, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      !CALL POPPOINTER8(cptr)
      cosa_u => gridstruct%cosa_u ! (/unknown_shape_in_d_sw/))
      !CALL POPPOINTER8(cptr)
      cosa_v => gridstruct%cosa_v ! (/unknown_shape_in_d_sw/))
      !CALL POPPOINTER8(cptr)
      cosa => gridstruct%cosa ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(ut, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL POPINTEGER(je1)
      !CALL POPPOINTER8(cptr)
      sin_sg => gridstruct%sin_sg ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(dt4)
      CALL POPREALARRAY(dt6)
      CALL POPINTEGER(ie)
      CALL POPREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(vb, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      CALL POPREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
      CALL POPREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      !CALL POPPOINTER8(cptr)
      rsin_u => gridstruct%rsin_u ! (/unknown_shape_in_d_sw/))
      CALL POPINTEGER(is)
      !CALL POPPOINTER8(cptr)
      rsina => gridstruct%rsina ! (/unknown_shape_in_d_sw/))
      !CALL POPPOINTER8(cptr)
      rdx => gridstruct%rdx ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      CALL POPINTEGER(js2)
      CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL POPINTEGER(min1)
      CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL POPINTEGER(min2)
      CALL POPINTEGER(min3)
      CALL POPINTEGER(ie1)
      CALL POPINTEGER(min4)
      CALL POPINTEGER(je)
      !CALL POPPOINTER8(cptr)
      rdxa => gridstruct%rdxa ! (/unknown_shape_in_d_sw/))
      CALL POPREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      ut_ad = 0.0
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(v(i, j))
          ut_ad(i, j) = ut_ad(i, j) - v_ad(i, j)
        END DO
      END DO
      vt_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(u(i, j))
          vt_ad(i, j) = vt_ad(i, j) + u_ad(i, j)
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(v(i, j))
          !ut_ad(i, j) = ut_ad(i, j) - v_ad(i, j)
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(u(i, j))
          !vt_ad(i, j) = vt_ad(i, j) + u_ad(i, j)
        END DO
      END DO
    END IF
    rsin2 => gridstruct%rsin2
    cosa_s => gridstruct%cosa_s
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      gx_ad = 0.0
      gy_ad = 0.0
      ub_ad = 0.0
      vb_ad = 0.0
      fx_ad = 0.0
      fy_ad = 0.0
      DO j=je,js,-1
        DO i=ie,is,-1
          dv2 = vb(i, j) + vb(i+1, j)
          v2 = fx(i, j) + fx(i+1, j)
          du2 = ub(i, j) + ub(i, j+1)
          u2 = fy(i, j) + fy(i, j+1)
          CALL POPREALARRAY(heat_source(i, j))
          temp0 = damp*rsin2(i, j)
          temp_ad88 = delp(i, j)*heat_source_ad(i, j)
          temp_ad89 = -(temp0*temp_ad88)
          temp_ad90 = 2.*temp_ad89
          temp_ad91 = -(cosa_s(i, j)*temp_ad89)
          delp_ad(i, j) = delp_ad(i, j) + (heat_source(i, j)-temp0*(ub(i&
&           , j)**2+ub(i, j+1)**2+vb(i, j)**2+vb(i+1, j)**2+2.*(gy(i, j)&
&           +gy(i, j+1)+gx(i, j)+gx(i+1, j))-cosa_s(i, j)*(u2*dv2+v2*du2&
&           +du2*dv2)))*heat_source_ad(i, j)
          ub_ad(i, j) = ub_ad(i, j) + 2*ub(i, j)*temp_ad89
          ub_ad(i, j+1) = ub_ad(i, j+1) + 2*ub(i, j+1)*temp_ad89
          vb_ad(i, j) = vb_ad(i, j) + 2*vb(i, j)*temp_ad89
          vb_ad(i+1, j) = vb_ad(i+1, j) + 2*vb(i+1, j)*temp_ad89
          gy_ad(i, j) = gy_ad(i, j) + temp_ad90
          gy_ad(i, j+1) = gy_ad(i, j+1) + temp_ad90
          gx_ad(i, j) = gx_ad(i, j) + temp_ad90
          gx_ad(i+1, j) = gx_ad(i+1, j) + temp_ad90
          u2_ad = dv2*temp_ad91
          dv2_ad = (du2+u2)*temp_ad91
          v2_ad = du2*temp_ad91
          du2_ad = (dv2+v2)*temp_ad91
          heat_source_ad(i, j) = temp_ad88
          vb_ad(i, j) = vb_ad(i, j) + dv2_ad
          vb_ad(i+1, j) = vb_ad(i+1, j) + dv2_ad
          fx_ad(i, j) = fx_ad(i, j) + v2_ad
          fx_ad(i+1, j) = fx_ad(i+1, j) + v2_ad
          ub_ad(i, j) = ub_ad(i, j) + du2_ad
          ub_ad(i, j+1) = ub_ad(i, j+1) + du2_ad
          fy_ad(i, j) = fy_ad(i, j) + u2_ad
          fy_ad(i, j+1) = fy_ad(i, j+1) + u2_ad
        END DO
      END DO
      rdy => gridstruct%rdy
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(gx(i, j))
          fx_ad(i, j) = fx_ad(i, j) + vb(i, j)*gx_ad(i, j)
          vb_ad(i, j) = vb_ad(i, j) + fx(i, j)*gx_ad(i, j)
          gx_ad(i, j) = 0.0
          CALL POPREALARRAY(fx(i, j))
          v_ad(i, j) = v_ad(i, j) + rdy(i, j)*fx_ad(i, j)
          fx_ad(i, j) = 0.0
          CALL POPREALARRAY(vb(i, j))
          temp_ad87 = rdy(i, j)*vb_ad(i, j)
          ut_ad(i, j) = ut_ad(i, j) - temp_ad87
          vb_ad(i, j) = temp_ad87
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(gy(i, j))
          fy_ad(i, j) = fy_ad(i, j) + ub(i, j)*gy_ad(i, j)
          ub_ad(i, j) = ub_ad(i, j) + fy(i, j)*gy_ad(i, j)
          gy_ad(i, j) = 0.0
          CALL POPREALARRAY(fy(i, j))
          u_ad(i, j) = u_ad(i, j) + rdx(i, j)*fy_ad(i, j)
          fy_ad(i, j) = 0.0
          CALL POPREALARRAY(ub(i, j))
          temp_ad86 = rdx(i, j)*ub_ad(i, j)
          vt_ad(i, j) = vt_ad(i, j) + temp_ad86
          ub_ad(i, j) = temp_ad86
        END DO
      END DO
    ELSE
      gx_ad = 0.0
      gy_ad = 0.0
      ub_ad = 0.0
      vb_ad = 0.0
      fx_ad = 0.0
      fy_ad = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      damp4 = (damp_v_pert*gridstruct%da_min_c)**(nord_v_pert+1)
      npx = flagstruct%npx
      npy = flagstruct%npy
      CALL POPREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(ut, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      wk_ad = 0.0
      vort_ad = 0.0
      CALL DEL6_VT_FLUX_ADM(nord_v_pert, npx, npy, damp4, wk, wk_ad, &
&                     vort, vort_ad, ut, ut_ad, vt, vt_ad, gridstruct, &
&                     bd)
    ELSE
      vort_ad = 0.0
      wk_ad = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      damp4 = (damp_v*gridstruct%da_min_c)**(nord_v+1)
      CALL POPREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(ut, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(vt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
!      CALL DEL6_VT_FLUX_ADM(nord_v, npx, npy, damp4, wk, wk_ad, vort, &
!&                     vort_ad, ut, ut_ad, vt, vt_ad, gridstruct, bd)
    END IF
    ke_ad = 0.0
    DO j=je,js,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(v(i, j))
        ut_ad(i, j) = ut_ad(i, j) + v_ad(i, j)
        ke_ad(i, j) = ke_ad(i, j) + v_ad(i, j)
        ke_ad(i, j+1) = ke_ad(i, j+1) - v_ad(i, j)
        fx_ad(i, j) = fx_ad(i, j) - v_ad(i, j)
        v_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ie,is,-1
        CALL POPREALARRAY(u(i, j))
        vt_ad(i, j) = vt_ad(i, j) + u_ad(i, j)
        ke_ad(i, j) = ke_ad(i, j) + u_ad(i, j)
        fy_ad(i, j) = fy_ad(i, j) + u_ad(i, j)
        ke_ad(i+1, j) = ke_ad(i+1, j) - u_ad(i, j)
        u_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      ra_x_ad = 0.0
      ra_y_ad = 0.0
      CALL FV_TP_2D_ADM(vort, vort_ad, crx_adv, crx_adv_ad, cry_adv, &
&                 cry_adv_ad, npx, npy, hord_vt_pert, fx, fx_ad, fy, &
&                 fy_ad, xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, &
&                 gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad)
    ELSE
      ra_x_ad = 0.0
      ra_y_ad = 0.0
      CALL FV_TP_2D_BWD(vort, vort_ad, crx_adv, crx_adv_ad, cry_adv, &
&                    cry_adv_ad, npx, npy, hord_vt, fx, fx_ad, fy, fy_ad&
&                    , xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, &
&                    gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad)
    END IF
    f0 => gridstruct%f0
    jsd = bd%jsd
    ied = bd%ied
    isd = bd%isd
    jed = bd%jed
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      DO j=jed,jsd,-1
        DO i=ied,isd,-1
          wk_ad(i, j) = wk_ad(i, j) + vort_ad(i, j)
          vort_ad(i, j) = 0.0
        END DO
      END DO
    ELSE IF (branch .EQ. 1) THEN
      DO j=jed,jsd,-1
        DO i=ied,isd,-1
          wk_ad(i, j) = wk_ad(i, j) + vort_ad(i, j)
          z_rat_ad(i, j) = z_rat_ad(i, j) + f0(i, j)*vort_ad(i, j)
          vort_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      DO j=jed,jsd,-1
        DO i=ied,isd,-1
          wk_ad(i, j) = wk_ad(i, j) + vort_ad(i, j)
          vort_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(vb(i, j))
          vort_ad(i, j) = vort_ad(i, j) + vb_ad(i, j)
          vort_ad(i, j+1) = vort_ad(i, j+1) - vb_ad(i, j)
          vb_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(ub(i, j))
          vort_ad(i, j) = vort_ad(i, j) + ub_ad(i, j)
          vort_ad(i+1, j) = vort_ad(i+1, j) - ub_ad(i, j)
          ub_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL COMPUTE_DIVERGENCE_DAMPING_BWD(nord, d2_bg, d4_bg, dddmp, &
&                                      dt, vort, vort_ad, ptc, ptc_ad, &
&                                      delpc, delpc_ad, ke, ke_ad, u, &
&                                      u_ad, v, v_ad, uc, uc_ad, vc, &
&                                      vc_ad, ua, ua_ad, va, va_ad, &
&                                      divg_d, divg_d_ad, wk, wk_ad, &
&                                      gridstruct, flagstruct, bd)
    ELSE
      CALL POPREALARRAY(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(ptc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(delpc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(uc, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(vc, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      CALL POPREALARRAY(divg_d, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+2))
      CALL POPREALARRAY(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL COMPUTE_DIVERGENCE_DAMPING_ADM(nord_pert, d2_bg_pert, &
&                                   d4_bg_pert, dddmp_pert, dt, vort, &
&                                   vort_ad, ptc, ptc_ad, delpc, &
&                                   delpc_ad, ke, ke_ad, u, u_ad, v, &
&                                   v_ad, uc, uc_ad, vc, vc_ad, ua, &
&                                   ua_ad, va, va_ad, divg_d, divg_d_ad&
&                                   , wk, wk_ad, gridstruct, flagstruct&
&                                   , bd)
    END IF
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      dw_ad = 0.0
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(w(i, j))
          dw_ad(i, j) = dw_ad(i, j) + w_ad(i, j)
        END DO
      END DO
    ELSE IF (branch .EQ. 1) THEN
      dw_ad = 0.0
    ELSE
      dw_ad = 0.0
      GOTO 100
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(w(i, j))
          temp_ad85 = w_ad(i, j)/delp(i, j)
          delp_ad(i, j) = delp_ad(i, j) - w(i, j)*temp_ad85/delp(i, j)
          w_ad(i, j) = temp_ad85
        END DO
      END DO
    END IF
 100 rarea => gridstruct%rarea
    DO j=jed,jsd,-1
      DO i=ied,isd,-1
        CALL POPREALARRAY(wk(i, j))
        temp_ad84 = rarea(i, j)*wk_ad(i, j)
        vt_ad(i, j) = vt_ad(i, j) + temp_ad84
        vt_ad(i, j+1) = vt_ad(i, j+1) - temp_ad84
        ut_ad(i+1, j) = ut_ad(i+1, j) + temp_ad84
        ut_ad(i, j) = ut_ad(i, j) - temp_ad84
        wk_ad(i, j) = 0.0
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ied+1,isd,-1
        CALL POPREALARRAY(ut(i, j))
        v_ad(i, j) = v_ad(i, j) + dy(i, j)*ut_ad(i, j)
        ut_ad(i, j) = 0.0
      END DO
    END DO
    DO j=jed+1,jsd,-1
      DO i=ied,isd,-1
        CALL POPREALARRAY(vt(i, j))
        u_ad(i, j) = u_ad(i, j) + dx(i, j)*vt_ad(i, j)
        vt_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL(2,branch)
    IF (branch .NE. 0) THEN
      IF (branch .NE. 1) THEN
        dt6 = dt/6.
        temp_ad80 = dt6*ke_ad(1, npy)
        temp_ad81 = u(1, npy)*temp_ad80
        temp_ad82 = v(1, npy-1)*temp_ad80
        temp_ad83 = u(0, npy)*temp_ad80
        ut_ad(1, npy) = ut_ad(1, npy) + temp_ad81
        ut_ad(1, npy-1) = ut_ad(1, npy-1) + temp_ad83 + temp_ad81
        u_ad(1, npy) = u_ad(1, npy) + (ut(1, npy)+ut(1, npy-1))*&
&         temp_ad80
        vt_ad(1, npy) = vt_ad(1, npy) + temp_ad82 - temp_ad83
        vt_ad(0, npy) = vt_ad(0, npy) + temp_ad82
        v_ad(1, npy-1) = v_ad(1, npy-1) + (vt(1, npy)+vt(0, npy))*&
&         temp_ad80
        u_ad(0, npy) = u_ad(0, npy) + (ut(1, npy-1)-vt(1, npy))*&
&         temp_ad80
        ke_ad(1, npy) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        temp_ad76 = dt6*ke_ad(npx, npy)
        temp_ad77 = u(npx-1, npy)*temp_ad76
        temp_ad78 = v(npx, npy-1)*temp_ad76
        temp_ad79 = u(npx, npy)*temp_ad76
        ut_ad(npx, npy) = ut_ad(npx, npy) + temp_ad77
        ut_ad(npx, npy-1) = ut_ad(npx, npy-1) + temp_ad79 + temp_ad77
        u_ad(npx-1, npy) = u_ad(npx-1, npy) + (ut(npx, npy)+ut(npx, npy-&
&         1))*temp_ad76
        vt_ad(npx, npy) = vt_ad(npx, npy) + temp_ad78
        vt_ad(npx-1, npy) = vt_ad(npx-1, npy) + temp_ad79 + temp_ad78
        v_ad(npx, npy-1) = v_ad(npx, npy-1) + (vt(npx, npy)+vt(npx-1, &
&         npy))*temp_ad76
        u_ad(npx, npy) = u_ad(npx, npy) + (ut(npx, npy-1)+vt(npx-1, npy)&
&         )*temp_ad76
        ke_ad(npx, npy) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        temp_ad72 = dt6*ke_ad(npx, 1)
        temp_ad73 = u(npx-1, 1)*temp_ad72
        temp_ad74 = v(npx, 1)*temp_ad72
        temp_ad75 = u(npx, 1)*temp_ad72
        ut_ad(npx, 1) = ut_ad(npx, 1) + temp_ad75 + temp_ad73
        ut_ad(npx, 0) = ut_ad(npx, 0) + temp_ad73
        u_ad(npx-1, 1) = u_ad(npx-1, 1) + (ut(npx, 1)+ut(npx, 0))*&
&         temp_ad72
        vt_ad(npx, 1) = vt_ad(npx, 1) + temp_ad74
        vt_ad(npx-1, 1) = vt_ad(npx-1, 1) + temp_ad74 - temp_ad75
        v_ad(npx, 1) = v_ad(npx, 1) + (vt(npx, 1)+vt(npx-1, 1))*&
&         temp_ad72
        u_ad(npx, 1) = u_ad(npx, 1) + (ut(npx, 1)-vt(npx-1, 1))*&
&         temp_ad72
        ke_ad(npx, 1) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        temp_ad71 = dt6*ke_ad(1, 1)
        ut_ad(1, 1) = ut_ad(1, 1) + (u(0, 1)+u(1, 1))*temp_ad71
        ut_ad(1, 0) = ut_ad(1, 0) + u(1, 1)*temp_ad71
        u_ad(1, 1) = u_ad(1, 1) + (ut(1, 1)+ut(1, 0))*temp_ad71
        vt_ad(1, 1) = vt_ad(1, 1) + (u(0, 1)+v(1, 1))*temp_ad71
        vt_ad(0, 1) = vt_ad(0, 1) + v(1, 1)*temp_ad71
        v_ad(1, 1) = v_ad(1, 1) + (vt(1, 1)+vt(0, 1))*temp_ad71
        u_ad(0, 1) = u_ad(0, 1) + (ut(1, 1)+vt(1, 1))*temp_ad71
        ke_ad(1, 1) = 0.0
      END IF
    END IF
    nested = gridstruct%nested
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        temp_ad70 = 0.5*ke_ad(i, j)
        ub_ad(i, j) = ub_ad(i, j) + vb(i, j)*temp_ad70
        vb_ad(i, j) = vb_ad(i, j) + ub(i, j)*temp_ad70
        ke_ad(i, j) = temp_ad70
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(vb, (bd%ie-bd%is+2)*(bd%je-bd%js+2))
      CALL XTP_U_ADM(is, ie, js, je, isd, ied, jsd, jed, ub, ub_ad, u, &
&              u_ad, v, vb, vb_ad, hord_mt_pert, gridstruct%dx, &
&              gridstruct%rdx, npx, npy, flagstruct%grid_type, nested)
    ELSE
      CALL XTP_U_BWD(is, ie, js, je, isd, ied, jsd, jed, ub, ub_ad, u&
&                 , u_ad, v, vb, vb_ad, hord_mt, gridstruct%dx, &
&                 gridstruct%rdx, npx, npy, flagstruct%grid_type, nested&
&                )
    END IF
    dt5 = 0.5*dt
    CALL POPCONTROL(2,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) THEN
        DO j=je+1,js,-1
          DO i=ie1,is2,-1
            CALL POPREALARRAY(ub(i, j))
            temp_ad65 = dt5*rsina(i, j)*ub_ad(i, j)
            temp_ad66 = -(cosa(i, j)*temp_ad65)
            uc_ad(i, j-1) = uc_ad(i, j-1) + temp_ad65
            uc_ad(i, j) = uc_ad(i, j) + temp_ad65
            vc_ad(i-1, j) = vc_ad(i-1, j) + temp_ad66
            vc_ad(i, j) = vc_ad(i, j) + temp_ad66
            ub_ad(i, j) = 0.0
          END DO
        END DO
        GOTO 110
      ELSE
        DO j=je+1,js,-1
          CALL POPREALARRAY(ub(npx, j))
          ut_ad(npx, j-1) = ut_ad(npx, j-1) + dt5*ub_ad(npx, j)
          ut_ad(npx, j) = ut_ad(npx, j) + dt5*ub_ad(npx, j)
          ub_ad(npx, j) = 0.0
        END DO
      END IF
    ELSE IF (branch .NE. 2) THEN
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(ub(i, j))
          uc_ad(i, j-1) = uc_ad(i, j-1) + dt5*ub_ad(i, j)
          uc_ad(i, j) = uc_ad(i, j) + dt5*ub_ad(i, j)
          ub_ad(i, j) = 0.0
        END DO
      END DO
      GOTO 110
    END IF
    cosa => gridstruct%cosa
    dt4 = 0.25*dt
    rsina => gridstruct%rsina
    DO j=je+1,js,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO i=ie1,is2,-1
          CALL POPREALARRAY(ub(i, j))
          temp_ad68 = dt5*rsina(i, j)*ub_ad(i, j)
          temp_ad69 = -(cosa(i, j)*temp_ad68)
          uc_ad(i, j-1) = uc_ad(i, j-1) + temp_ad68
          uc_ad(i, j) = uc_ad(i, j) + temp_ad68
          vc_ad(i-1, j) = vc_ad(i-1, j) + temp_ad69
          vc_ad(i, j) = vc_ad(i, j) + temp_ad69
          ub_ad(i, j) = 0.0
        END DO
      ELSE
        DO i=ie1,is2,-1
          CALL POPREALARRAY(ub(i, j))
          temp_ad67 = dt4*ub_ad(i, j)
          ut_ad(i, j-1) = ut_ad(i, j-1) + 3.*temp_ad67
          ut_ad(i, j) = ut_ad(i, j) + 3.*temp_ad67
          ut_ad(i, j-2) = ut_ad(i, j-2) - temp_ad67
          ut_ad(i, j+1) = ut_ad(i, j+1) - temp_ad67
          ub_ad(i, j) = 0.0
        END DO
      END IF
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      DO j=je+1,js,-1
        CALL POPREALARRAY(ub(1, j))
        ut_ad(1, j-1) = ut_ad(1, j-1) + dt5*ub_ad(1, j)
        ut_ad(1, j) = ut_ad(1, j) + dt5*ub_ad(1, j)
        ub_ad(1, j) = 0.0
      END DO
    END IF
 110 DO j=je+1,js,-1
      DO i=ie+1,is,-1
        vb_ad(i, j) = vb_ad(i, j) + ub(i, j)*ke_ad(i, j)
        ub_ad(i, j) = ub_ad(i, j) + vb(i, j)*ke_ad(i, j)
        ke_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL YTP_V_ADM(is, ie, js, je, isd, ied, jsd, jed, vb, vb_ad, u, v&
&              , v_ad, ub, ub_ad, hord_mt_pert, gridstruct%dy, &
&              gridstruct%rdy, npx, npy, flagstruct%grid_type, nested)
    ELSE
      CALL YTP_V_BWD(is, ie, js, je, isd, ied, jsd, jed, vb, vb_ad, u&
&                 , v, v_ad, ub, ub_ad, hord_mt, gridstruct%dy, &
&                 gridstruct%rdy, npx, npy, flagstruct%grid_type, nested&
&                )
    END IF
    CALL POPCONTROL(2,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) THEN
        DO j=je1,js2,-1
          DO i=ie1,is2,-1
            temp_ad59 = dt5*rsina(i, j)*vb_ad(i, j)
            temp_ad60 = -(cosa(i, j)*temp_ad59)
            vc_ad(i-1, j) = vc_ad(i-1, j) + temp_ad59
            vc_ad(i, j) = vc_ad(i, j) + temp_ad59
            uc_ad(i, j-1) = uc_ad(i, j-1) + temp_ad60
            uc_ad(i, j) = uc_ad(i, j) + temp_ad60
            vb_ad(i, j) = 0.0
          END DO
        END DO
        GOTO 120
      ELSE
        DO i=ie+1,is,-1
          vt_ad(i-1, npy) = vt_ad(i-1, npy) + dt5*vb_ad(i, npy)
          vt_ad(i, npy) = vt_ad(i, npy) + dt5*vb_ad(i, npy)
          vb_ad(i, npy) = 0.0
        END DO
      END IF
    ELSE IF (branch .NE. 2) THEN
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          vc_ad(i-1, j) = vc_ad(i-1, j) + dt5*vb_ad(i, j)
          vc_ad(i, j) = vc_ad(i, j) + dt5*vb_ad(i, j)
          vb_ad(i, j) = 0.0
        END DO
      END DO
      GOTO 120
    END IF
    DO j=je1,js2,-1
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        temp_ad64 = dt4*vb_ad(npx, j)
        vt_ad(npx-1, j) = vt_ad(npx-1, j) + 3.*temp_ad64
        vt_ad(npx, j) = vt_ad(npx, j) + 3.*temp_ad64
        vt_ad(npx-2, j) = vt_ad(npx-2, j) - temp_ad64
        vt_ad(npx+1, j) = vt_ad(npx+1, j) - temp_ad64
        vb_ad(npx, j) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        temp_ad63 = dt4*vb_ad(1, j)
        vt_ad(0, j) = vt_ad(0, j) + 3.*temp_ad63
        vt_ad(1, j) = vt_ad(1, j) + 3.*temp_ad63
        vt_ad(-1, j) = vt_ad(-1, j) - temp_ad63
        vt_ad(2, j) = vt_ad(2, j) - temp_ad63
        vb_ad(1, j) = 0.0
      END IF
      DO i=ie1,is2,-1
        temp_ad61 = dt5*rsina(i, j)*vb_ad(i, j)
        temp_ad62 = -(cosa(i, j)*temp_ad61)
        vc_ad(i-1, j) = vc_ad(i-1, j) + temp_ad61
        vc_ad(i, j) = vc_ad(i, j) + temp_ad61
        uc_ad(i, j-1) = uc_ad(i, j-1) + temp_ad62
        uc_ad(i, j) = uc_ad(i, j) + temp_ad62
        vb_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      DO i=ie+1,is,-1
        vt_ad(i-1, 1) = vt_ad(i-1, 1) + dt5*vb_ad(i, 1)
        vt_ad(i, 1) = vt_ad(i, 1) + dt5*vb_ad(i, 1)
        vb_ad(i, 1) = 0.0
      END DO
    END IF
 120 CALL POPCONTROL(1,branch)
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        DO i=ie,is,-1
          temp_ad58 = rarea(i, j)*dpx_ad(i, j)
          fx_ad(i, j) = fx_ad(i, j) + temp_ad58
          fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad58
          fy_ad(i, j) = fy_ad(i, j) + temp_ad58
          fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad58
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO iq=nq,1,-1
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(q(i, j, k, iq))
            temp_ad53 = q_ad(i, j, k, iq)/delp(i, j)
            temp = q(i, j, k, iq)
            temp_ad54 = rarea(i, j)*temp_ad53
            wk_ad(i, j) = wk_ad(i, j) + temp*temp_ad53
            gx_ad(i, j) = gx_ad(i, j) + temp_ad54
            gx_ad(i+1, j) = gx_ad(i+1, j) - temp_ad54
            gy_ad(i, j) = gy_ad(i, j) + temp_ad54
            gy_ad(i, j+1) = gy_ad(i, j+1) - temp_ad54
            delp_ad(i, j) = delp_ad(i, j) - (temp*wk(i, j)+rarea(i, j)*(&
&             gx(i, j)-gx(i+1, j)+gy(i, j)-gy(i, j+1)))*temp_ad53/delp(i&
&             , j)
            q_ad(i, j, k, iq) = wk(i, j)*temp_ad53
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1)*(&
&                      jed-jsd+1))
          CALL POPREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
          CALL POPREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
          CALL FV_TP_2D_ADM(q(isd:ied, jsd:jed, k, iq), q_ad(isd:ied, &
&                     jsd:jed, k, iq), crx_adv, crx_adv_ad, cry_adv, &
&                     cry_adv_ad, npx, npy, hord_tr_pert, gx, gx_ad, gy&
&                     , gy_ad, xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad&
&                     , gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad, &
&                     mfx=fx, mfx_ad=fx_ad, mfy=fy, mfy_ad=fy_ad, mass=&
&                     delp, mass_ad=delp_ad, nord=nord_t_pert, damp_c=&
&                     damp_t_pert)
        ELSE
          CALL FV_TP_2D_BWD(q(isd:ied, jsd:jed, k, iq), q_ad(isd:ied&
&                        , jsd:jed, k, iq), crx_adv, crx_adv_ad, cry_adv&
&                        , cry_adv_ad, npx, npy, hord_tr, gx, gx_ad, gy&
&                        , gy_ad, xfx_adv, xfx_adv_ad, yfx_adv, &
&                        yfx_adv_ad, gridstruct, bd, ra_x, ra_x_ad, ra_y&
&                        , ra_y_ad, mfx=fx, mfx_ad=fx_ad, mfy=fy, mfy_ad=fy_ad&
&                        , mass=delp, mass_ad=delp_ad, nord=nord_t, damp_c=&
&                        damp_t)
        END IF
      END DO
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(pt(i, j))
          temp_ad50 = pt_ad(i, j)/delp(i, j)
          temp_ad51 = rarea(i, j)*temp_ad50
          gx_ad(i, j) = gx_ad(i, j) + temp_ad51
          gx_ad(i+1, j) = gx_ad(i+1, j) - temp_ad51
          gy_ad(i, j) = gy_ad(i, j) + temp_ad51
          gy_ad(i, j+1) = gy_ad(i, j+1) - temp_ad51
          delp_ad(i, j) = delp_ad(i, j) - (pt(i, j)*wk(i, j)+rarea(i, j)&
&           *(gx(i, j)-gx(i+1, j)+gy(i, j)-gy(i, j+1)))*temp_ad50/delp(i&
&           , j)
          wk_ad(i, j) = wk_ad(i, j) + delp_ad(i, j) + pt(i, j)*temp_ad50
          pt_ad(i, j) = wk(i, j)*temp_ad50
          CALL POPREALARRAY(delp(i, j))
          temp_ad52 = rarea(i, j)*delp_ad(i, j)
          fx_ad(i, j) = fx_ad(i, j) + temp_ad52
          fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad52
          fy_ad(i, j) = fy_ad(i, j) + temp_ad52
          fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad52
          delp_ad(i, j) = wk_ad(i, j)
          wk_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(pt(i, j))
          temp_ad55 = pt_ad(i, j)/delp(i, j)
          delp_ad(i, j) = delp_ad(i, j) - pt(i, j)*temp_ad55/delp(i, j)
          pt_ad(i, j) = temp_ad55
          CALL POPREALARRAY(delp(i, j))
          temp_ad56 = rarea(i, j)*delp_ad(i, j)
          fx_ad(i, j) = fx_ad(i, j) + temp_ad56
          fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad56
          fy_ad(i, j) = fy_ad(i, j) + temp_ad56
          fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad56
          CALL POPREALARRAY(pt(i, j))
          temp_ad57 = rarea(i, j)*pt_ad(i, j)
          delp_ad(i, j) = delp_ad(i, j) + pt(i, j)*pt_ad(i, j)
          gx_ad(i, j) = gx_ad(i, j) + temp_ad57
          gx_ad(i+1, j) = gx_ad(i+1, j) - temp_ad57
          gy_ad(i, j) = gy_ad(i, j) + temp_ad57
          gy_ad(i, j+1) = gy_ad(i, j+1) - temp_ad57
          pt_ad(i, j) = delp(i, j)*pt_ad(i, j)
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL FV_TP_2D_BWD(pt, pt_ad, crx_adv, crx_adv_ad, cry_adv, &
&                    cry_adv_ad, npx, npy, hord_tm, gx, gx_ad, gy, gy_ad&
&                    , xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, &
&                    gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad, mfx=&
&                    fx, mfx_ad=fx_ad, mfy=fy, mfy_ad=fy_ad, mass=delp, mass_ad=&
&                    delp_ad, nord=nord_t, damp_c=damp_t)
    ELSE
      CALL POPREALARRAY(pt, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL POPREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL FV_TP_2D_ADM(pt, pt_ad, crx_adv, crx_adv_ad, cry_adv, &
&                 cry_adv_ad, npx, npy, hord_tm_pert, gx, gx_ad, gy, &
&                 gy_ad, xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, &
&                 gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad, mfx=fx, &
&                 mfx_ad=fx_ad, mfy=fy, mfy_ad=fy_ad, mass=delp, mass_ad&
&                 =delp_ad, nord=nord_t_pert, damp_c=damp_t_pert)
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(w(i, j))
          temp_ad49 = rarea(i, j)*w_ad(i, j)
          delp_ad(i, j) = delp_ad(i, j) + w(i, j)*w_ad(i, j)
          gx_ad(i, j) = gx_ad(i, j) + temp_ad49
          gx_ad(i+1, j) = gx_ad(i+1, j) - temp_ad49
          gy_ad(i, j) = gy_ad(i, j) + temp_ad49
          gy_ad(i, j+1) = gy_ad(i, j+1) - temp_ad49
          w_ad(i, j) = delp(i, j)*w_ad(i, j)
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(w, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL POPREALARRAY(gx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
        CALL POPREALARRAY(gy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
        CALL FV_TP_2D_ADM(w, w_ad, crx_adv, crx_adv_ad, cry_adv, &
&                   cry_adv_ad, npx, npy, hord_vt_pert, gx, gx_ad, gy, &
&                   gy_ad, xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, &
&                   gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad, mfx=fx&
&                   , mfx_ad=fx_ad, mfy=fy, mfy_ad=fy_ad)
      ELSE
        CALL FV_TP_2D_BWD(w, w_ad, crx_adv, crx_adv_ad, cry_adv, &
&                      cry_adv_ad, npx, npy, hord_vt, gx, gx_ad, gy, &
&                      gy_ad, xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, &
&                      gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad, mfx&
&                      =fx, mfx_ad=fx_ad, mfy=fy, mfy_ad=fy_ad)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        fy2_ad = 0.0
        fx2_ad = 0.0
        DO j=je,js,-1
          DO i=ie,is,-1
            temp_ad47 = -(dw(i, j)*heat_source_ad(i, j))
            dw_ad(i, j) = dw_ad(i, j) + 0.5*temp_ad47 - (w(i, j)+0.5*dw(&
&             i, j))*heat_source_ad(i, j)
            w_ad(i, j) = w_ad(i, j) + temp_ad47
            heat_source_ad(i, j) = 0.0
            temp_ad48 = rarea(i, j)*dw_ad(i, j)
            fx2_ad(i, j) = fx2_ad(i, j) + temp_ad48
            fx2_ad(i+1, j) = fx2_ad(i+1, j) - temp_ad48
            fy2_ad(i, j) = fy2_ad(i, j) + temp_ad48
            fy2_ad(i, j+1) = fy2_ad(i, j+1) - temp_ad48
            dw_ad(i, j) = 0.0
          END DO
        END DO
        damp4 = (damp_w*gridstruct%da_min_c)**(nord_w+1)
        CALL DEL6_VT_FLUX_ADM(nord_w, npx, npy, damp4, w, w_ad, wk, &
&                       wk_ad, fx2, fx2_ad, fy2, fy2_ad, gridstruct, bd)
      END IF
    END IF
    DO j=je,js,-1
      DO i=ie,is,-1
        heat_source_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ie,is,-1
        fy_ad(i, j) = fy_ad(i, j) + yflux_ad(i, j)
      END DO
      DO i=ied,isd,-1
        cry_adv_ad(i, j) = cry_adv_ad(i, j) + cy_ad(i, j)
      END DO
    END DO
    DO j=je,js,-1
      DO i=ie+1,is,-1
        fx_ad(i, j) = fx_ad(i, j) + xflux_ad(i, j)
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ie+1,is,-1
        crx_adv_ad(i, j) = crx_adv_ad(i, j) + cx_ad(i, j)
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(delp, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL FV_TP_2D_ADM(delp, delp_ad, crx_adv, crx_adv_ad, cry_adv, &
&                 cry_adv_ad, npx, npy, hord_dp_pert, fx, fx_ad, fy, &
&                 fy_ad, xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, &
&                 gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad, nord=&
&                 nord_v_pert, damp_c=damp_v_pert)
    ELSE
      CALL FV_TP_2D_BWD(delp, delp_ad, crx_adv, crx_adv_ad, cry_adv, &
&                    cry_adv_ad, npx, npy, hord_dp, fx, fx_ad, fy, fy_ad&
&                    , xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, &
&                    gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad, nord=&
&                    nord_v, damp_c=damp_v)
    END IF
    DO j=je,js,-1
      DO i=ied,isd,-1
        yfx_adv_ad(i, j) = yfx_adv_ad(i, j) + ra_y_ad(i, j)
        yfx_adv_ad(i, j+1) = yfx_adv_ad(i, j+1) - ra_y_ad(i, j)
        ra_y_ad(i, j) = 0.0
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ie,is,-1
        xfx_adv_ad(i, j) = xfx_adv_ad(i, j) + ra_x_ad(i, j)
        xfx_adv_ad(i+1, j) = xfx_adv_ad(i+1, j) - ra_x_ad(i, j)
        ra_x_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ied,isd,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(yfx_adv(i, j))
          yfx_adv_ad(i, j) = rdya(i, j)*cry_adv_ad(i, j) + sin_sg(i, j, &
&           2)*dx(i, j)*yfx_adv_ad(i, j)
          CALL POPREALARRAY(cry_adv(i, j))
          cry_adv_ad(i, j) = 0.0
        ELSE
          CALL POPREALARRAY(yfx_adv(i, j))
          yfx_adv_ad(i, j) = rdya(i, j-1)*cry_adv_ad(i, j) + sin_sg(i, j&
&           -1, 4)*dx(i, j)*yfx_adv_ad(i, j)
          CALL POPREALARRAY(cry_adv(i, j))
          cry_adv_ad(i, j) = 0.0
        END IF
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ie+1,is,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(xfx_adv(i, j))
          xfx_adv_ad(i, j) = rdxa(i, j)*crx_adv_ad(i, j) + sin_sg(i, j, &
&           1)*dy(i, j)*xfx_adv_ad(i, j)
          CALL POPREALARRAY(crx_adv(i, j))
          crx_adv_ad(i, j) = 0.0
        ELSE
          CALL POPREALARRAY(xfx_adv(i, j))
          xfx_adv_ad(i, j) = rdxa(i-1, j)*crx_adv_ad(i, j) + sin_sg(i-1&
&           , j, 3)*dy(i, j)*xfx_adv_ad(i, j)
          CALL POPREALARRAY(crx_adv(i, j))
          crx_adv_ad(i, j) = 0.0
        END IF
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ied,isd,-1
        CALL POPREALARRAY(yfx_adv(i, j))
        vt_ad(i, j) = vt_ad(i, j) + dt*yfx_adv_ad(i, j)
        yfx_adv_ad(i, j) = 0.0
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(xfx_adv(i, j))
        ut_ad(i, j) = ut_ad(i, j) + dt*xfx_adv_ad(i, j)
        xfx_adv_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL(2,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) THEN
        DO j=je+1,js,-1
          DO i=ied,isd,-1
            vc_ad(i, j) = vc_ad(i, j) + vt_ad(i, j)
            vt_ad(i, j) = 0.0
          END DO
        END DO
        DO j=jed,jsd,-1
          DO i=ie+1,is,-1
            uc_ad(i, j) = uc_ad(i, j) + ut_ad(i, j)
            ut_ad(i, j) = 0.0
          END DO
        END DO
        GOTO 130
      END IF
    ELSE
      IF (branch .NE. 2) THEN
        cosa_u => gridstruct%cosa_u
        cosa_v => gridstruct%cosa_v
        damp = 1./(1.-0.0625*cosa_u(2, npy-1)*cosa_v(1, npy-1))
        temp_ad39 = -(damp*cosa_v(1, npy-1)*0.25*vt_ad(1, npy-1))
        temp_ad40 = -(cosa_u(2, npy-1)*0.25*temp_ad39)
        ut_ad(1, npy-1) = ut_ad(1, npy-1) + temp_ad39
        ut_ad(1, npy-2) = ut_ad(1, npy-2) + temp_ad39
        ut_ad(2, npy-2) = ut_ad(2, npy-2) + temp_ad39
        uc_ad(2, npy-1) = uc_ad(2, npy-1) + damp*ut_ad(2, npy-1) + &
&         temp_ad39
        temp_ad41 = -(damp*cosa_u(2, npy-1)*0.25*ut_ad(2, npy-1))
        vc_ad(1, npy-1) = vc_ad(1, npy-1) + temp_ad41 + damp*vt_ad(1, &
&         npy-1)
        vt_ad(1, npy) = vt_ad(1, npy) + temp_ad40
        vt_ad(2, npy) = vt_ad(2, npy) + temp_ad40
        vt_ad(2, npy-1) = vt_ad(2, npy-1) + temp_ad40
        vt_ad(1, npy-1) = 0.0
        temp_ad42 = -(cosa_v(1, npy-1)*0.25*temp_ad41)
        vt_ad(1, npy) = vt_ad(1, npy) + temp_ad41
        vt_ad(2, npy) = vt_ad(2, npy) + temp_ad41
        vt_ad(2, npy-1) = vt_ad(2, npy-1) + temp_ad41
        ut_ad(1, npy-1) = ut_ad(1, npy-1) + temp_ad42
        ut_ad(1, npy-2) = ut_ad(1, npy-2) + temp_ad42
        ut_ad(2, npy-2) = ut_ad(2, npy-2) + temp_ad42
        ut_ad(2, npy-1) = 0.0
        damp = 1./(1.-0.0625*cosa_u(0, npy-1)*cosa_v(0, npy-1))
        temp_ad43 = -(damp*cosa_v(0, npy-1)*0.25*vt_ad(0, npy-1))
        temp_ad44 = -(cosa_u(0, npy-1)*0.25*temp_ad43)
        vc_ad(0, npy-1) = vc_ad(0, npy-1) + damp*vt_ad(0, npy-1)
        ut_ad(1, npy-1) = ut_ad(1, npy-1) + temp_ad43
        ut_ad(1, npy-2) = ut_ad(1, npy-2) + temp_ad43
        ut_ad(0, npy-2) = ut_ad(0, npy-2) + temp_ad43
        uc_ad(0, npy-1) = uc_ad(0, npy-1) + temp_ad43
        vt_ad(0, npy) = vt_ad(0, npy) + temp_ad44
        vt_ad(-1, npy) = vt_ad(-1, npy) + temp_ad44
        vt_ad(-1, npy-1) = vt_ad(-1, npy-1) + temp_ad44
        vt_ad(0, npy-1) = 0.0
        damp = 1./(1.-0.0625*cosa_u(2, npy)*cosa_v(1, npy+1))
        temp_ad45 = -(damp*cosa_u(2, npy)*0.25*ut_ad(2, npy))
        temp_ad46 = -(cosa_v(1, npy+1)*0.25*temp_ad45)
        uc_ad(2, npy) = uc_ad(2, npy) + damp*ut_ad(2, npy)
        vt_ad(1, npy) = vt_ad(1, npy) + temp_ad45
        vt_ad(2, npy) = vt_ad(2, npy) + temp_ad45
        vt_ad(2, npy+1) = vt_ad(2, npy+1) + temp_ad45
        vc_ad(1, npy+1) = vc_ad(1, npy+1) + temp_ad45
        ut_ad(1, npy) = ut_ad(1, npy) + temp_ad46
        ut_ad(1, npy+1) = ut_ad(1, npy+1) + temp_ad46
        ut_ad(2, npy+1) = ut_ad(2, npy+1) + temp_ad46
        ut_ad(2, npy) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        damp = 1./(1.-0.0625*cosa_u(npx-1, npy-1)*cosa_v(npx-1, npy-1))
        temp_ad31 = -(damp*cosa_v(npx-1, npy-1)*0.25*vt_ad(npx-1, npy-1)&
&         )
        temp_ad32 = -(cosa_u(npx-1, npy-1)*0.25*temp_ad31)
        ut_ad(npx, npy-1) = ut_ad(npx, npy-1) + temp_ad31
        ut_ad(npx, npy-2) = ut_ad(npx, npy-2) + temp_ad31
        ut_ad(npx-1, npy-2) = ut_ad(npx-1, npy-2) + temp_ad31
        uc_ad(npx-1, npy-1) = uc_ad(npx-1, npy-1) + damp*ut_ad(npx-1, &
&         npy-1) + temp_ad31
        temp_ad33 = -(damp*cosa_u(npx-1, npy-1)*0.25*ut_ad(npx-1, npy-1)&
&         )
        vc_ad(npx-1, npy-1) = vc_ad(npx-1, npy-1) + temp_ad33 + damp*&
&         vt_ad(npx-1, npy-1)
        vt_ad(npx-1, npy) = vt_ad(npx-1, npy) + temp_ad32
        vt_ad(npx-2, npy) = vt_ad(npx-2, npy) + temp_ad32
        vt_ad(npx-2, npy-1) = vt_ad(npx-2, npy-1) + temp_ad32
        vt_ad(npx-1, npy-1) = 0.0
        temp_ad34 = -(cosa_v(npx-1, npy-1)*0.25*temp_ad33)
        vt_ad(npx-1, npy) = vt_ad(npx-1, npy) + temp_ad33
        vt_ad(npx-2, npy) = vt_ad(npx-2, npy) + temp_ad33
        vt_ad(npx-2, npy-1) = vt_ad(npx-2, npy-1) + temp_ad33
        ut_ad(npx, npy-1) = ut_ad(npx, npy-1) + temp_ad34
        ut_ad(npx, npy-2) = ut_ad(npx, npy-2) + temp_ad34
        ut_ad(npx-1, npy-2) = ut_ad(npx-1, npy-2) + temp_ad34
        ut_ad(npx-1, npy-1) = 0.0
        damp = 1./(1.-0.0625*cosa_u(npx+1, npy-1)*cosa_v(npx, npy-1))
        temp_ad35 = -(damp*cosa_v(npx, npy-1)*0.25*vt_ad(npx, npy-1))
        temp_ad36 = -(cosa_u(npx+1, npy-1)*0.25*temp_ad35)
        vc_ad(npx, npy-1) = vc_ad(npx, npy-1) + damp*vt_ad(npx, npy-1)
        ut_ad(npx, npy-1) = ut_ad(npx, npy-1) + temp_ad35
        ut_ad(npx, npy-2) = ut_ad(npx, npy-2) + temp_ad35
        ut_ad(npx+1, npy-2) = ut_ad(npx+1, npy-2) + temp_ad35
        uc_ad(npx+1, npy-1) = uc_ad(npx+1, npy-1) + temp_ad35
        vt_ad(npx, npy) = vt_ad(npx, npy) + temp_ad36
        vt_ad(npx+1, npy) = vt_ad(npx+1, npy) + temp_ad36
        vt_ad(npx+1, npy-1) = vt_ad(npx+1, npy-1) + temp_ad36
        vt_ad(npx, npy-1) = 0.0
        damp = 1./(1.-0.0625*cosa_u(npx-1, npy)*cosa_v(npx-1, npy+1))
        temp_ad37 = -(damp*cosa_u(npx-1, npy)*0.25*ut_ad(npx-1, npy))
        temp_ad38 = -(cosa_v(npx-1, npy+1)*0.25*temp_ad37)
        uc_ad(npx-1, npy) = uc_ad(npx-1, npy) + damp*ut_ad(npx-1, npy)
        vt_ad(npx-1, npy) = vt_ad(npx-1, npy) + temp_ad37
        vt_ad(npx-2, npy) = vt_ad(npx-2, npy) + temp_ad37
        vt_ad(npx-2, npy+1) = vt_ad(npx-2, npy+1) + temp_ad37
        vc_ad(npx-1, npy+1) = vc_ad(npx-1, npy+1) + temp_ad37
        ut_ad(npx, npy) = ut_ad(npx, npy) + temp_ad38
        ut_ad(npx, npy+1) = ut_ad(npx, npy+1) + temp_ad38
        ut_ad(npx-1, npy+1) = ut_ad(npx-1, npy+1) + temp_ad38
        ut_ad(npx-1, npy) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        damp = 1./(1.-0.0625*cosa_u(npx-1, 1)*cosa_v(npx-1, 2))
        temp_ad23 = -(damp*cosa_v(npx-1, 2)*0.25*vt_ad(npx-1, 2))
        temp_ad24 = -(cosa_u(npx-1, 1)*0.25*temp_ad23)
        ut_ad(npx, 1) = ut_ad(npx, 1) + temp_ad23
        ut_ad(npx, 2) = ut_ad(npx, 2) + temp_ad23
        ut_ad(npx-1, 2) = ut_ad(npx-1, 2) + temp_ad23
        uc_ad(npx-1, 1) = uc_ad(npx-1, 1) + damp*ut_ad(npx-1, 1) + &
&         temp_ad23
        temp_ad25 = -(damp*cosa_u(npx-1, 1)*0.25*ut_ad(npx-1, 1))
        vc_ad(npx-1, 2) = vc_ad(npx-1, 2) + temp_ad25 + damp*vt_ad(npx-1&
&         , 2)
        vt_ad(npx-1, 1) = vt_ad(npx-1, 1) + temp_ad24
        vt_ad(npx-2, 1) = vt_ad(npx-2, 1) + temp_ad24
        vt_ad(npx-2, 2) = vt_ad(npx-2, 2) + temp_ad24
        vt_ad(npx-1, 2) = 0.0
        temp_ad26 = -(cosa_v(npx-1, 2)*0.25*temp_ad25)
        vt_ad(npx-1, 1) = vt_ad(npx-1, 1) + temp_ad25
        vt_ad(npx-2, 1) = vt_ad(npx-2, 1) + temp_ad25
        vt_ad(npx-2, 2) = vt_ad(npx-2, 2) + temp_ad25
        ut_ad(npx, 1) = ut_ad(npx, 1) + temp_ad26
        ut_ad(npx, 2) = ut_ad(npx, 2) + temp_ad26
        ut_ad(npx-1, 2) = ut_ad(npx-1, 2) + temp_ad26
        ut_ad(npx-1, 1) = 0.0
        damp = 1./(1.-0.0625*cosa_u(npx+1, 1)*cosa_v(npx, 2))
        temp_ad27 = -(damp*cosa_v(npx, 2)*0.25*vt_ad(npx, 2))
        temp_ad28 = -(cosa_u(npx+1, 1)*0.25*temp_ad27)
        vc_ad(npx, 2) = vc_ad(npx, 2) + damp*vt_ad(npx, 2)
        ut_ad(npx, 1) = ut_ad(npx, 1) + temp_ad27
        ut_ad(npx, 2) = ut_ad(npx, 2) + temp_ad27
        ut_ad(npx+1, 2) = ut_ad(npx+1, 2) + temp_ad27
        uc_ad(npx+1, 1) = uc_ad(npx+1, 1) + temp_ad27
        vt_ad(npx, 1) = vt_ad(npx, 1) + temp_ad28
        vt_ad(npx+1, 1) = vt_ad(npx+1, 1) + temp_ad28
        vt_ad(npx+1, 2) = vt_ad(npx+1, 2) + temp_ad28
        vt_ad(npx, 2) = 0.0
        damp = 1./(1.-0.0625*cosa_u(npx-1, 0)*cosa_v(npx-1, 0))
        temp_ad29 = -(damp*cosa_u(npx-1, 0)*0.25*ut_ad(npx-1, 0))
        temp_ad30 = -(cosa_v(npx-1, 0)*0.25*temp_ad29)
        uc_ad(npx-1, 0) = uc_ad(npx-1, 0) + damp*ut_ad(npx-1, 0)
        vt_ad(npx-1, 1) = vt_ad(npx-1, 1) + temp_ad29
        vt_ad(npx-2, 1) = vt_ad(npx-2, 1) + temp_ad29
        vt_ad(npx-2, 0) = vt_ad(npx-2, 0) + temp_ad29
        vc_ad(npx-1, 0) = vc_ad(npx-1, 0) + temp_ad29
        ut_ad(npx, 0) = ut_ad(npx, 0) + temp_ad30
        ut_ad(npx, -1) = ut_ad(npx, -1) + temp_ad30
        ut_ad(npx-1, -1) = ut_ad(npx-1, -1) + temp_ad30
        ut_ad(npx-1, 0) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        damp = 1./(1.-0.0625*cosa_u(2, 1)*cosa_v(1, 2))
        temp_ad15 = -(damp*cosa_v(1, 2)*0.25*vt_ad(1, 2))
        temp_ad16 = -(cosa_u(2, 1)*0.25*temp_ad15)
        ut_ad(1, 1) = ut_ad(1, 1) + temp_ad15
        ut_ad(1, 2) = ut_ad(1, 2) + temp_ad15
        ut_ad(2, 2) = ut_ad(2, 2) + temp_ad15
        uc_ad(2, 1) = uc_ad(2, 1) + damp*ut_ad(2, 1) + temp_ad15
        temp_ad17 = -(damp*cosa_u(2, 1)*0.25*ut_ad(2, 1))
        vc_ad(1, 2) = vc_ad(1, 2) + temp_ad17 + damp*vt_ad(1, 2)
        vt_ad(1, 1) = vt_ad(1, 1) + temp_ad16
        vt_ad(2, 1) = vt_ad(2, 1) + temp_ad16
        vt_ad(2, 2) = vt_ad(2, 2) + temp_ad16
        vt_ad(1, 2) = 0.0
        temp_ad18 = -(cosa_v(1, 2)*0.25*temp_ad17)
        vt_ad(1, 1) = vt_ad(1, 1) + temp_ad17
        vt_ad(2, 1) = vt_ad(2, 1) + temp_ad17
        vt_ad(2, 2) = vt_ad(2, 2) + temp_ad17
        ut_ad(1, 1) = ut_ad(1, 1) + temp_ad18
        ut_ad(1, 2) = ut_ad(1, 2) + temp_ad18
        ut_ad(2, 2) = ut_ad(2, 2) + temp_ad18
        ut_ad(2, 1) = 0.0
        damp = 1./(1.-0.0625*cosa_u(0, 1)*cosa_v(0, 2))
        temp_ad19 = -(damp*cosa_v(0, 2)*0.25*vt_ad(0, 2))
        temp_ad20 = -(cosa_u(0, 1)*0.25*temp_ad19)
        vc_ad(0, 2) = vc_ad(0, 2) + damp*vt_ad(0, 2)
        ut_ad(1, 1) = ut_ad(1, 1) + temp_ad19
        ut_ad(1, 2) = ut_ad(1, 2) + temp_ad19
        ut_ad(0, 2) = ut_ad(0, 2) + temp_ad19
        uc_ad(0, 1) = uc_ad(0, 1) + temp_ad19
        vt_ad(0, 1) = vt_ad(0, 1) + temp_ad20
        vt_ad(-1, 1) = vt_ad(-1, 1) + temp_ad20
        vt_ad(-1, 2) = vt_ad(-1, 2) + temp_ad20
        vt_ad(0, 2) = 0.0
        damp = 1./(1.-0.0625*cosa_u(2, 0)*cosa_v(1, 0))
        temp_ad21 = -(damp*cosa_u(2, 0)*0.25*ut_ad(2, 0))
        temp_ad22 = -(cosa_v(1, 0)*0.25*temp_ad21)
        uc_ad(2, 0) = uc_ad(2, 0) + damp*ut_ad(2, 0)
        vt_ad(1, 1) = vt_ad(1, 1) + temp_ad21
        vt_ad(2, 1) = vt_ad(2, 1) + temp_ad21
        vt_ad(2, 0) = vt_ad(2, 0) + temp_ad21
        vc_ad(1, 0) = vc_ad(1, 0) + temp_ad21
        ut_ad(1, 0) = ut_ad(1, 0) + temp_ad22
        ut_ad(1, -1) = ut_ad(1, -1) + temp_ad22
        ut_ad(2, -1) = ut_ad(2, -1) + temp_ad22
        ut_ad(2, 0) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO i=min4,max4,-1
          temp_ad13 = -(cosa_u(i, npy)*0.25*ut_ad(i, npy))
          uc_ad(i, npy) = uc_ad(i, npy) + ut_ad(i, npy)
          vt_ad(i-1, npy) = vt_ad(i-1, npy) + temp_ad13
          vt_ad(i, npy) = vt_ad(i, npy) + temp_ad13
          vt_ad(i-1, npy+1) = vt_ad(i-1, npy+1) + temp_ad13
          vt_ad(i, npy+1) = vt_ad(i, npy+1) + temp_ad13
          ut_ad(i, npy) = 0.0
          temp_ad14 = -(cosa_u(i, npy-1)*0.25*ut_ad(i, npy-1))
          uc_ad(i, npy-1) = uc_ad(i, npy-1) + ut_ad(i, npy-1)
          vt_ad(i-1, npy-1) = vt_ad(i-1, npy-1) + temp_ad14
          vt_ad(i, npy-1) = vt_ad(i, npy-1) + temp_ad14
          vt_ad(i-1, npy) = vt_ad(i-1, npy) + temp_ad14
          vt_ad(i, npy) = vt_ad(i, npy) + temp_ad14
          ut_ad(i, npy-1) = 0.0
        END DO
        DO i=ied,isd,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            vc_ad(i, npy) = vc_ad(i, npy) + vt_ad(i, npy)/sin_sg(i, npy&
&             , 2)
            vt_ad(i, npy) = 0.0
          ELSE
            vc_ad(i, npy) = vc_ad(i, npy) + vt_ad(i, npy)/sin_sg(i, npy-&
&             1, 4)
            vt_ad(i, npy) = 0.0
          END IF
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO i=min3,max3,-1
          temp_ad11 = -(cosa_u(i, 1)*0.25*ut_ad(i, 1))
          uc_ad(i, 1) = uc_ad(i, 1) + ut_ad(i, 1)
          vt_ad(i-1, 1) = vt_ad(i-1, 1) + temp_ad11
          vt_ad(i, 1) = vt_ad(i, 1) + temp_ad11
          vt_ad(i-1, 2) = vt_ad(i-1, 2) + temp_ad11
          vt_ad(i, 2) = vt_ad(i, 2) + temp_ad11
          ut_ad(i, 1) = 0.0
          temp_ad12 = -(cosa_u(i, 0)*0.25*ut_ad(i, 0))
          uc_ad(i, 0) = uc_ad(i, 0) + ut_ad(i, 0)
          vt_ad(i-1, 0) = vt_ad(i-1, 0) + temp_ad12
          vt_ad(i, 0) = vt_ad(i, 0) + temp_ad12
          vt_ad(i-1, 1) = vt_ad(i-1, 1) + temp_ad12
          vt_ad(i, 1) = vt_ad(i, 1) + temp_ad12
          ut_ad(i, 0) = 0.0
        END DO
        DO i=ied,isd,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            vc_ad(i, 1) = vc_ad(i, 1) + vt_ad(i, 1)/sin_sg(i, 1, 2)
            vt_ad(i, 1) = 0.0
          ELSE
            vc_ad(i, 1) = vc_ad(i, 1) + vt_ad(i, 1)/sin_sg(i, 0, 4)
            vt_ad(i, 1) = 0.0
          END IF
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=min2,max2,-1
          temp_ad9 = -(cosa_v(npx, j)*0.25*vt_ad(npx, j))
          vc_ad(npx, j) = vc_ad(npx, j) + vt_ad(npx, j)
          ut_ad(npx, j-1) = ut_ad(npx, j-1) + temp_ad9
          ut_ad(npx+1, j-1) = ut_ad(npx+1, j-1) + temp_ad9
          ut_ad(npx, j) = ut_ad(npx, j) + temp_ad9
          ut_ad(npx+1, j) = ut_ad(npx+1, j) + temp_ad9
          vt_ad(npx, j) = 0.0
          temp_ad10 = -(cosa_v(npx-1, j)*0.25*vt_ad(npx-1, j))
          vc_ad(npx-1, j) = vc_ad(npx-1, j) + vt_ad(npx-1, j)
          ut_ad(npx-1, j-1) = ut_ad(npx-1, j-1) + temp_ad10
          ut_ad(npx, j-1) = ut_ad(npx, j-1) + temp_ad10
          ut_ad(npx-1, j) = ut_ad(npx-1, j) + temp_ad10
          ut_ad(npx, j) = ut_ad(npx, j) + temp_ad10
          vt_ad(npx-1, j) = 0.0
        END DO
        DO j=jed,jsd,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            uc_ad(npx, j) = uc_ad(npx, j) + ut_ad(npx, j)/sin_sg(npx, j&
&             , 1)
            ut_ad(npx, j) = 0.0
          ELSE
            uc_ad(npx, j) = uc_ad(npx, j) + ut_ad(npx, j)/sin_sg(npx-1, &
&             j, 3)
            ut_ad(npx, j) = 0.0
          END IF
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=min1,max1,-1
          temp_ad7 = -(cosa_v(1, j)*0.25*vt_ad(1, j))
          vc_ad(1, j) = vc_ad(1, j) + vt_ad(1, j)
          ut_ad(1, j-1) = ut_ad(1, j-1) + temp_ad7
          ut_ad(2, j-1) = ut_ad(2, j-1) + temp_ad7
          ut_ad(1, j) = ut_ad(1, j) + temp_ad7
          ut_ad(2, j) = ut_ad(2, j) + temp_ad7
          vt_ad(1, j) = 0.0
          temp_ad8 = -(cosa_v(0, j)*0.25*vt_ad(0, j))
          vc_ad(0, j) = vc_ad(0, j) + vt_ad(0, j)
          ut_ad(0, j-1) = ut_ad(0, j-1) + temp_ad8
          ut_ad(1, j-1) = ut_ad(1, j-1) + temp_ad8
          ut_ad(0, j) = ut_ad(0, j) + temp_ad8
          ut_ad(1, j) = ut_ad(1, j) + temp_ad8
          vt_ad(0, j) = 0.0
        END DO
        DO j=jed,jsd,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            uc_ad(1, j) = uc_ad(1, j) + ut_ad(1, j)/sin_sg(1, j, 1)
            ut_ad(1, j) = 0.0
          ELSE
            uc_ad(1, j) = uc_ad(1, j) + ut_ad(1, j)/sin_sg(0, j, 3)
            ut_ad(1, j) = 0.0
          END IF
        END DO
      END IF
    END IF
    rsin_v => gridstruct%rsin_v
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=je+2,js-1,-1
        DO i=ied,isd,-1
          temp_ad1 = rsin_v(i, j)*vt_ad(i, j)
          temp_ad2 = -(cosa_v(i, j)*0.25*temp_ad1)
          vc_ad(i, j) = vc_ad(i, j) + temp_ad1
          uc_ad(i, j-1) = uc_ad(i, j-1) + temp_ad2
          uc_ad(i+1, j-1) = uc_ad(i+1, j-1) + temp_ad2
          uc_ad(i, j) = uc_ad(i, j) + temp_ad2
          uc_ad(i+1, j) = uc_ad(i+1, j) + temp_ad2
          vt_ad(i, j) = 0.0
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+2,is-1,-1
          temp_ad = rsin_u(i, j)*ut_ad(i, j)
          temp_ad0 = -(cosa_u(i, j)*0.25*temp_ad)
          uc_ad(i, j) = uc_ad(i, j) + temp_ad
          vc_ad(i-1, j) = vc_ad(i-1, j) + temp_ad0
          vc_ad(i, j) = vc_ad(i, j) + temp_ad0
          vc_ad(i-1, j+1) = vc_ad(i-1, j+1) + temp_ad0
          vc_ad(i, j+1) = vc_ad(i, j+1) + temp_ad0
          ut_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      DO j=je+2,js-1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          DO i=ied,isd,-1
            temp_ad5 = rsin_v(i, j)*vt_ad(i, j)
            temp_ad6 = -(cosa_v(i, j)*0.25*temp_ad5)
            vc_ad(i, j) = vc_ad(i, j) + temp_ad5
            uc_ad(i, j-1) = uc_ad(i, j-1) + temp_ad6
            uc_ad(i+1, j-1) = uc_ad(i+1, j-1) + temp_ad6
            uc_ad(i, j) = uc_ad(i, j) + temp_ad6
            uc_ad(i+1, j) = uc_ad(i+1, j) + temp_ad6
            vt_ad(i, j) = 0.0
          END DO
        END IF
      END DO
      DO j=jed,jsd,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          DO i=ie+2,is-1,-1
            temp_ad3 = rsin_u(i, j)*ut_ad(i, j)
            temp_ad4 = -(cosa_u(i, j)*0.25*temp_ad3)
            uc_ad(i, j) = uc_ad(i, j) + temp_ad3
            vc_ad(i-1, j) = vc_ad(i-1, j) + temp_ad4
            vc_ad(i, j) = vc_ad(i, j) + temp_ad4
            vc_ad(i-1, j+1) = vc_ad(i-1, j+1) + temp_ad4
            vc_ad(i, j+1) = vc_ad(i, j+1) + temp_ad4
            ut_ad(i, j) = 0.0
          END DO
        END IF
      END DO
    END IF
 130 CONTINUE
  END SUBROUTINE D_SW_BWD
!     d_sw :: D-Grid Shallow Water Routine
  SUBROUTINE D_SW(delpc, delp, ptc, pt, u, v, w, uc, vc, ua, va, divg_d&
&   , xflux, yflux, cx, cy, crx_adv, cry_adv, xfx_adv, yfx_adv, q_con, &
&   z_rat, kgb, heat_source, dpx, zvir, sphum, nq, q, k, km, inline_q, &
&   dt, hord_tr, hord_mt, hord_vt, hord_tm, hord_dp, nord, nord_v, &
&   nord_w, nord_t, dddmp, d2_bg, d4_bg, damp_v, damp_w, damp_t, d_con, &
&   hydrostatic, gridstruct, flagstruct, bd, hord_tr_pert, hord_mt_pert&
&   , hord_vt_pert, hord_tm_pert, hord_dp_pert, split_damp, nord_pert, &
&   nord_v_pert, nord_w_pert, nord_t_pert, dddmp_pert, d2_bg_pert, &
&   d4_bg_pert, damp_v_pert, damp_w_pert, damp_t_pert)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: hord_tr, hord_mt, hord_vt, hord_tm, hord_dp
! nord=1 divergence damping; (del-4) or 3 (del-8)
    INTEGER, INTENT(IN) :: nord
! vorticity damping
    INTEGER, INTENT(IN) :: nord_v
! vertical velocity
    INTEGER, INTENT(IN) :: nord_w
! pt
    INTEGER, INTENT(IN) :: nord_t
    INTEGER, INTENT(IN) :: sphum, nq, k, km
    REAL, INTENT(IN) :: dt, dddmp, d2_bg, d4_bg, d_con
    REAL, INTENT(IN) :: zvir
    REAL, INTENT(IN) :: damp_v, damp_w, damp_t, kgb
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: hord_tr_pert, hord_mt_pert, hord_vt_pert, &
&   hord_tm_pert, hord_dp_pert, nord_pert, nord_v_pert, nord_w_pert, &
&   nord_t_pert
    LOGICAL, INTENT(IN) :: split_damp
    REAL, INTENT(IN) :: dddmp_pert, d2_bg_pert, d4_bg_pert, damp_v_pert&
&   , damp_w_pert, damp_t_pert
! divergence
    REAL, INTENT(INOUT) :: divg_d(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: z_rat
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: delp&
&   , pt, ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   q_con
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: u&
&   , vc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: v&
&   , uc
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, km, nq)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: delpc&
&   , ptc
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je), INTENT(OUT) :: &
&   heat_source
    REAL(kind=8), DIMENSION(bd%is:bd%ie, bd%js:bd%je), INTENT(INOUT) :: &
&   dpx
! The flux capacitors:
    REAL, INTENT(INOUT) :: xflux(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: yflux(bd%is:bd%ie, bd%js:bd%je+1)
!------------------------
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1)
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: inline_q
    REAL, DIMENSION(bd%is:bd%ie+1, bd%jsd:bd%jed), INTENT(OUT) :: &
&   crx_adv, xfx_adv
    REAL, DIMENSION(bd%isd:bd%ied, bd%js:bd%je+1), INTENT(OUT) :: &
&   cry_adv, yfx_adv
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! Local:
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL :: ut(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!---
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: fy2(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!  work array
    REAL :: dw(bd%is:bd%ie, bd%js:bd%je)
!---
    REAL, DIMENSION(bd%is:bd%ie+1, bd%js:bd%je+1) :: ub, vb
!  work array
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
!  needs this for corner_comm
    REAL :: ke(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
! Vorticity
    REAL :: vort(bd%isd:bd%ied, bd%jsd:bd%jed)
! 1-D X-direction Fluxes
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
! 1-D Y-direction Fluxes
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: gx(bd%is:bd%ie+1, bd%js:bd%je)
! work Y-dir flux array
    REAL :: gy(bd%is:bd%ie, bd%js:bd%je+1)
    LOGICAL :: fill_c
    REAL :: dt2, dt4, dt5, dt6
    REAL :: damp, damp2, damp4, dd8, u2, v2, du2, dv2
    REAL :: u_lon
    INTEGER :: i, j, is2, ie1, js2, je1, n, nt, n2, iq
    REAL, DIMENSION(:, :), POINTER :: area, area_c, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsina
    REAL, DIMENSION(:, :), POINTER :: f0, rsin2, divg_u, divg_v
    REAL, DIMENSION(:, :), POINTER :: cosa, dx, dy, dxc, dyc, rdxa, rdya&
&   , rdx, rdy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: npx, npy
    LOGICAL :: nested
    REAL :: delp_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: pt_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: vort_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: delpc_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: ptc_tj(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: ke_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL :: vc_tj(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: uc_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: divg_d_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL :: ut_tj(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt_tj(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    REAL :: abs0
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    cosa_s => gridstruct%cosa_s
    sina_u => gridstruct%sina_u
    sina_v => gridstruct%sina_v
    rsin_u => gridstruct%rsin_u
    rsin_v => gridstruct%rsin_v
    rsina => gridstruct%rsina
    f0 => gridstruct%f0
    rsin2 => gridstruct%rsin2
    divg_u => gridstruct%divg_u
    divg_v => gridstruct%divg_v
    cosa => gridstruct%cosa
    dx => gridstruct%dx
    dy => gridstruct%dy
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    rdxa => gridstruct%rdxa
    rdya => gridstruct%rdya
    rdx => gridstruct%rdx
    rdy => gridstruct%rdy
    sw_corner = gridstruct%sw_corner
    se_corner = gridstruct%se_corner
    nw_corner = gridstruct%nw_corner
    ne_corner = gridstruct%ne_corner
! end grid_type choices
    IF (flagstruct%grid_type .LT. 3) THEN
!!! TO DO: separate versions for nesting and for cubed-sphere
      IF (nested) THEN
        DO j=jsd,jed
          DO i=is-1,ie+2
            ut(i, j) = (uc(i, j)-0.25*cosa_u(i, j)*(vc(i-1, j)+vc(i, j)+&
&             vc(i-1, j+1)+vc(i, j+1)))*rsin_u(i, j)
          END DO
        END DO
        DO j=js-1,je+2
          DO i=isd,ied
            vt(i, j) = (vc(i, j)-0.25*cosa_v(i, j)*(uc(i, j-1)+uc(i+1, j&
&             -1)+uc(i, j)+uc(i+1, j)))*rsin_v(i, j)
          END DO
        END DO
      ELSE
        DO j=jsd,jed
          IF (j .NE. 0 .AND. j .NE. 1 .AND. j .NE. npy - 1 .AND. j .NE. &
&             npy) THEN
            DO i=is-1,ie+2
              ut(i, j) = (uc(i, j)-0.25*cosa_u(i, j)*(vc(i-1, j)+vc(i, j&
&               )+vc(i-1, j+1)+vc(i, j+1)))*rsin_u(i, j)
            END DO
          END IF
        END DO
        DO j=js-1,je+2
          IF (j .NE. 1 .AND. j .NE. npy) THEN
            DO i=isd,ied
              vt(i, j) = (vc(i, j)-0.25*cosa_v(i, j)*(uc(i, j-1)+uc(i+1&
&               , j-1)+uc(i, j)+uc(i+1, j)))*rsin_v(i, j)
            END DO
          END IF
        END DO
      END IF
!.not. nested
      IF (.NOT.nested) THEN
! West face
! West edge:
        IF (is .EQ. 1) THEN
          DO j=jsd,jed
            IF (uc(1, j)*dt .GT. 0.) THEN
              ut(1, j) = uc(1, j)/sin_sg(0, j, 3)
            ELSE
              ut(1, j) = uc(1, j)/sin_sg(1, j, 1)
            END IF
          END DO
          IF (3 .LT. js) THEN
            max1 = js
          ELSE
            max1 = 3
          END IF
          IF (npy - 2 .GT. je + 1) THEN
            min1 = je + 1
          ELSE
            min1 = npy - 2
          END IF
          DO j=max1,min1
            vt(0, j) = vc(0, j) - 0.25*cosa_v(0, j)*(ut(0, j-1)+ut(1, j-&
&             1)+ut(0, j)+ut(1, j))
            vt(1, j) = vc(1, j) - 0.25*cosa_v(1, j)*(ut(1, j-1)+ut(2, j-&
&             1)+ut(1, j)+ut(2, j))
          END DO
        END IF
! East edge:
        IF (ie + 1 .EQ. npx) THEN
          DO j=jsd,jed
            IF (uc(npx, j)*dt .GT. 0.) THEN
              ut(npx, j) = uc(npx, j)/sin_sg(npx-1, j, 3)
            ELSE
              ut(npx, j) = uc(npx, j)/sin_sg(npx, j, 1)
            END IF
          END DO
          IF (3 .LT. js) THEN
            max2 = js
          ELSE
            max2 = 3
          END IF
          IF (npy - 2 .GT. je + 1) THEN
            min2 = je + 1
          ELSE
            min2 = npy - 2
          END IF
          DO j=max2,min2
            vt(npx-1, j) = vc(npx-1, j) - 0.25*cosa_v(npx-1, j)*(ut(npx-&
&             1, j-1)+ut(npx, j-1)+ut(npx-1, j)+ut(npx, j))
            vt(npx, j) = vc(npx, j) - 0.25*cosa_v(npx, j)*(ut(npx, j-1)+&
&             ut(npx+1, j-1)+ut(npx, j)+ut(npx+1, j))
          END DO
        END IF
! South (Bottom) edge:
        IF (js .EQ. 1) THEN
          DO i=isd,ied
            IF (vc(i, 1)*dt .GT. 0.) THEN
              vt(i, 1) = vc(i, 1)/sin_sg(i, 0, 4)
            ELSE
              vt(i, 1) = vc(i, 1)/sin_sg(i, 1, 2)
            END IF
          END DO
          IF (3 .LT. is) THEN
            max3 = is
          ELSE
            max3 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min3 = ie + 1
          ELSE
            min3 = npx - 2
          END IF
          DO i=max3,min3
            ut(i, 0) = uc(i, 0) - 0.25*cosa_u(i, 0)*(vt(i-1, 0)+vt(i, 0)&
&             +vt(i-1, 1)+vt(i, 1))
            ut(i, 1) = uc(i, 1) - 0.25*cosa_u(i, 1)*(vt(i-1, 1)+vt(i, 1)&
&             +vt(i-1, 2)+vt(i, 2))
          END DO
        END IF
! North edge:
        IF (je + 1 .EQ. npy) THEN
          DO i=isd,ied
            IF (vc(i, npy)*dt .GT. 0.) THEN
              vt(i, npy) = vc(i, npy)/sin_sg(i, npy-1, 4)
            ELSE
              vt(i, npy) = vc(i, npy)/sin_sg(i, npy, 2)
            END IF
          END DO
          IF (3 .LT. is) THEN
            max4 = is
          ELSE
            max4 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min4 = ie + 1
          ELSE
            min4 = npx - 2
          END IF
          DO i=max4,min4
            ut(i, npy-1) = uc(i, npy-1) - 0.25*cosa_u(i, npy-1)*(vt(i-1&
&             , npy-1)+vt(i, npy-1)+vt(i-1, npy)+vt(i, npy))
            ut(i, npy) = uc(i, npy) - 0.25*cosa_u(i, npy)*(vt(i-1, npy)+&
&             vt(i, npy)+vt(i-1, npy+1)+vt(i, npy+1))
          END DO
        END IF
! The following code solves a 2x2 system to get the interior parallel-to-edge uc,vc values
! near the corners (ex: for the sw corner ut(2,1) and vt(1,2) are solved for simultaneously).
! It then computes the halo uc, vc values so as to be consistent with the computations on
! the facing panel.
!The system solved is:
!  ut(2,1) = uc(2,1) - avg(vt)*cosa_u(2,1)
!  vt(1,2) = vc(1,2) - avg(ut)*cosa_v(1,2)
! in which avg(vt) includes vt(1,2) and avg(ut) includes ut(2,1)
        IF (sw_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(2, 0)*cosa_v(1, 0))
          ut(2, 0) = (uc(2, 0)-0.25*cosa_u(2, 0)*(vt(1, 1)+vt(2, 1)+vt(2&
&           , 0)+vc(1, 0)-0.25*cosa_v(1, 0)*(ut(1, 0)+ut(1, -1)+ut(2, -1&
&           ))))*damp
          damp = 1./(1.-0.0625*cosa_u(0, 1)*cosa_v(0, 2))
          vt(0, 2) = (vc(0, 2)-0.25*cosa_v(0, 2)*(ut(1, 1)+ut(1, 2)+ut(0&
&           , 2)+uc(0, 1)-0.25*cosa_u(0, 1)*(vt(0, 1)+vt(-1, 1)+vt(-1, 2&
&           ))))*damp
          damp = 1./(1.-0.0625*cosa_u(2, 1)*cosa_v(1, 2))
          ut(2, 1) = (uc(2, 1)-0.25*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(2&
&           , 2)+vc(1, 2)-0.25*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(2, 2))&
&           ))*damp
          vt(1, 2) = (vc(1, 2)-0.25*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(2&
&           , 2)+uc(2, 1)-0.25*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(2, 2))&
&           ))*damp
        END IF
        IF (se_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(npx-1, 0)*cosa_v(npx-1, 0))
          ut(npx-1, 0) = (uc(npx-1, 0)-0.25*cosa_u(npx-1, 0)*(vt(npx-1, &
&           1)+vt(npx-2, 1)+vt(npx-2, 0)+vc(npx-1, 0)-0.25*cosa_v(npx-1&
&           , 0)*(ut(npx, 0)+ut(npx, -1)+ut(npx-1, -1))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx+1, 1)*cosa_v(npx, 2))
          vt(npx, 2) = (vc(npx, 2)-0.25*cosa_v(npx, 2)*(ut(npx, 1)+ut(&
&           npx, 2)+ut(npx+1, 2)+uc(npx+1, 1)-0.25*cosa_u(npx+1, 1)*(vt(&
&           npx, 1)+vt(npx+1, 1)+vt(npx+1, 2))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx-1, 1)*cosa_v(npx-1, 2))
          ut(npx-1, 1) = (uc(npx-1, 1)-0.25*cosa_u(npx-1, 1)*(vt(npx-1, &
&           1)+vt(npx-2, 1)+vt(npx-2, 2)+vc(npx-1, 2)-0.25*cosa_v(npx-1&
&           , 2)*(ut(npx, 1)+ut(npx, 2)+ut(npx-1, 2))))*damp
          vt(npx-1, 2) = (vc(npx-1, 2)-0.25*cosa_v(npx-1, 2)*(ut(npx, 1)&
&           +ut(npx, 2)+ut(npx-1, 2)+uc(npx-1, 1)-0.25*cosa_u(npx-1, 1)*&
&           (vt(npx-1, 1)+vt(npx-2, 1)+vt(npx-2, 2))))*damp
        END IF
        IF (ne_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(npx-1, npy)*cosa_v(npx-1, npy+1))
          ut(npx-1, npy) = (uc(npx-1, npy)-0.25*cosa_u(npx-1, npy)*(vt(&
&           npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy+1)+vc(npx-1, npy+1)&
&           -0.25*cosa_v(npx-1, npy+1)*(ut(npx, npy)+ut(npx, npy+1)+ut(&
&           npx-1, npy+1))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx+1, npy-1)*cosa_v(npx, npy-1))
          vt(npx, npy-1) = (vc(npx, npy-1)-0.25*cosa_v(npx, npy-1)*(ut(&
&           npx, npy-1)+ut(npx, npy-2)+ut(npx+1, npy-2)+uc(npx+1, npy-1)&
&           -0.25*cosa_u(npx+1, npy-1)*(vt(npx, npy)+vt(npx+1, npy)+vt(&
&           npx+1, npy-1))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx-1, npy-1)*cosa_v(npx-1, npy-1)&
&           )
          ut(npx-1, npy-1) = (uc(npx-1, npy-1)-0.25*cosa_u(npx-1, npy-1)&
&           *(vt(npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy-1)+vc(npx-1, &
&           npy-1)-0.25*cosa_v(npx-1, npy-1)*(ut(npx, npy-1)+ut(npx, npy&
&           -2)+ut(npx-1, npy-2))))*damp
          vt(npx-1, npy-1) = (vc(npx-1, npy-1)-0.25*cosa_v(npx-1, npy-1)&
&           *(ut(npx, npy-1)+ut(npx, npy-2)+ut(npx-1, npy-2)+uc(npx-1, &
&           npy-1)-0.25*cosa_u(npx-1, npy-1)*(vt(npx-1, npy)+vt(npx-2, &
&           npy)+vt(npx-2, npy-1))))*damp
        END IF
        IF (nw_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(2, npy)*cosa_v(1, npy+1))
          ut(2, npy) = (uc(2, npy)-0.25*cosa_u(2, npy)*(vt(1, npy)+vt(2&
&           , npy)+vt(2, npy+1)+vc(1, npy+1)-0.25*cosa_v(1, npy+1)*(ut(1&
&           , npy)+ut(1, npy+1)+ut(2, npy+1))))*damp
          damp = 1./(1.-0.0625*cosa_u(0, npy-1)*cosa_v(0, npy-1))
          vt(0, npy-1) = (vc(0, npy-1)-0.25*cosa_v(0, npy-1)*(ut(1, npy-&
&           1)+ut(1, npy-2)+ut(0, npy-2)+uc(0, npy-1)-0.25*cosa_u(0, npy&
&           -1)*(vt(0, npy)+vt(-1, npy)+vt(-1, npy-1))))*damp
          damp = 1./(1.-0.0625*cosa_u(2, npy-1)*cosa_v(1, npy-1))
          ut(2, npy-1) = (uc(2, npy-1)-0.25*cosa_u(2, npy-1)*(vt(1, npy)&
&           +vt(2, npy)+vt(2, npy-1)+vc(1, npy-1)-0.25*cosa_v(1, npy-1)*&
&           (ut(1, npy-1)+ut(1, npy-2)+ut(2, npy-2))))*damp
          vt(1, npy-1) = (vc(1, npy-1)-0.25*cosa_v(1, npy-1)*(ut(1, npy-&
&           1)+ut(1, npy-2)+ut(2, npy-2)+uc(2, npy-1)-0.25*cosa_u(2, npy&
&           -1)*(vt(1, npy)+vt(2, npy)+vt(2, npy-1))))*damp
        END IF
      END IF
    ELSE
! flagstruct%grid_type >= 3
      DO j=jsd,jed
        DO i=is,ie+1
          ut(i, j) = uc(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          vt(i, j) = vc(i, j)
        END DO
      END DO
    END IF
    DO j=jsd,jed
      DO i=is,ie+1
        xfx_adv(i, j) = dt*ut(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        yfx_adv(i, j) = dt*vt(i, j)
      END DO
    END DO
! Explanation of the following code:
!    xfx_adv = dt*ut*dy
!    crx_adv = dt*ut/dx
    DO j=jsd,jed
!DEC$ VECTOR ALWAYS
      DO i=is,ie+1
        IF (xfx_adv(i, j) .GT. 0.) THEN
          crx_adv(i, j) = xfx_adv(i, j)*rdxa(i-1, j)
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sin_sg(i-1, j, 3)
        ELSE
          crx_adv(i, j) = xfx_adv(i, j)*rdxa(i, j)
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sin_sg(i, j, 1)
        END IF
      END DO
    END DO
    DO j=js,je+1
!DEC$ VECTOR ALWAYS
      DO i=isd,ied
        IF (yfx_adv(i, j) .GT. 0.) THEN
          cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j-1)
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sin_sg(i, j-1, 4)
        ELSE
          cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j)
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sin_sg(i, j, 2)
        END IF
      END DO
    END DO
    DO j=jsd,jed
      DO i=is,ie
        ra_x(i, j) = area(i, j) + (xfx_adv(i, j)-xfx_adv(i+1, j))
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        ra_y(i, j) = area(i, j) + (yfx_adv(i, j)-yfx_adv(i, j+1))
      END DO
    END DO
    IF (hord_dp .EQ. hord_dp_pert .AND. (.NOT.split_damp)) THEN
      CALL FV_TP_2D(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx, fy&
&                , xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, nord=&
&                nord_v, damp_c=damp_v)
    ELSE
      CALL FV_TP_2D(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx, &
&             fy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, nord=&
&             nord_v, damp_c=damp_v)
    END IF
! <<< Save the mass fluxes to the "Flux Capacitor" for tracer transport >>>
    DO j=jsd,jed
      DO i=is,ie+1
        cx(i, j) = cx(i, j) + crx_adv(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie+1
        xflux(i, j) = xflux(i, j) + fx(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        cy(i, j) = cy(i, j) + cry_adv(i, j)
      END DO
      DO i=is,ie
        yflux(i, j) = yflux(i, j) + fy(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie
        heat_source(i, j) = 0.
      END DO
    END DO
    IF (.NOT.hydrostatic) THEN
      IF (damp_w .GT. 1.e-5) THEN
        IF (dt .GE. 0.) THEN
          abs0 = dt
        ELSE
          abs0 = -dt
        END IF
        dd8 = kgb*abs0
        damp4 = (damp_w*gridstruct%da_min_c)**(nord_w+1)
        CALL DEL6_VT_FLUX(nord_w, npx, npy, damp4, w, wk, fx2, fy2, &
&                   gridstruct, bd)
        DO j=js,je
          DO i=is,ie
            dw(i, j) = (fx2(i, j)-fx2(i+1, j)+(fy2(i, j)-fy2(i, j+1)))*&
&             rarea(i, j)
! 0.5 * [ (w+dw)**2 - w**2 ] = w*dw + 0.5*dw*dw
!                   heat_source(i,j) = -d_con*dw(i,j)*(w(i,j)+0.5*dw(i,j))
            heat_source(i, j) = dd8 - dw(i, j)*(w(i, j)+0.5*dw(i, j))
          END DO
        END DO
      END IF
      IF (hord_vt .EQ. hord_vt_pert) THEN
        CALL FV_TP_2D(w, crx_adv, cry_adv, npx, npy, hord_vt, gx, gy&
&                  , xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, mfx=&
&                  fx, mfy=fy)
      ELSE
        CALL FV_TP_2D(w, crx_adv, cry_adv, npx, npy, hord_vt, gx, &
&               gy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, mfx=fx&
&               , mfy=fy)
      END IF
      DO j=js,je
        DO i=is,ie
          w(i, j) = delp(i, j)*w(i, j) + (gx(i, j)-gx(i+1, j)+(gy(i, j)-&
&           gy(i, j+1)))*rarea(i, j)
        END DO
      END DO
    END IF
!    if ( inline_q .and. zvir>0.01 ) then
!       do j=jsd,jed
!          do i=isd,ied
!             pt(i,j) = pt(i,j)/(1.+zvir*q(i,j,k,sphum))
!          enddo
!       enddo
!    endif
    IF (hord_tm .EQ. hord_tm_pert .AND. (.NOT.split_damp)) THEN
      CALL FV_TP_2D(pt, crx_adv, cry_adv, npx, npy, hord_tm, gx, gy, &
&                xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, fx, fy, &
&                delp, nord_t, damp_t)
    ELSE
      CALL FV_TP_2D(pt, crx_adv, cry_adv, npx, npy, hord_tm, gx, gy&
&             , xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y, fx, fy, &
&             delp, nord_t, damp_t)
    END IF
    IF (inline_q) THEN
      DO j=js,je
        DO i=is,ie
          wk(i, j) = delp(i, j)
          delp(i, j) = wk(i, j) + (fx(i, j)-fx(i+1, j)+(fy(i, j)-fy(i, j&
&           +1)))*rarea(i, j)
          pt(i, j) = (pt(i, j)*wk(i, j)+(gx(i, j)-gx(i+1, j)+(gy(i, j)-&
&           gy(i, j+1)))*rarea(i, j))/delp(i, j)
        END DO
      END DO
      DO iq=1,nq
        IF (hord_tr .EQ. hord_tr_pert) THEN
          CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), crx_adv, cry_adv&
&                    , npx, npy, hord_tr, gx, gy, xfx_adv, yfx_adv, &
&                    gridstruct, bd, ra_x, ra_y, fx, fy, delp, nord_t, &
&                    damp_t)
        ELSE
          CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), crx_adv, cry_adv, &
&                 npx, npy, hord_tr, gx, gy, xfx_adv, yfx_adv, &
&                 gridstruct, bd, ra_x, ra_y, fx, fy, delp, nord_t&
&                 , damp_t)
        END IF
        DO j=js,je
          DO i=is,ie
            q(i, j, k, iq) = (q(i, j, k, iq)*wk(i, j)+(gx(i, j)-gx(i+1, &
&             j)+(gy(i, j)-gy(i, j+1)))*rarea(i, j))/delp(i, j)
          END DO
        END DO
      END DO
    ELSE
!     if ( zvir>0.01 ) then
!       do j=js,je
!          do i=is,ie
!             pt(i,j) = pt(i,j)*(1.+zvir*q(i,j,k,sphum))
!          enddo
!       enddo
!     endif
      DO j=js,je
        DO i=is,ie
          pt(i, j) = pt(i, j)*delp(i, j) + (gx(i, j)-gx(i+1, j)+(gy(i, j&
&           )-gy(i, j+1)))*rarea(i, j)
          delp(i, j) = delp(i, j) + (fx(i, j)-fx(i+1, j)+(fy(i, j)-fy(i&
&           , j+1)))*rarea(i, j)
          pt(i, j) = pt(i, j)/delp(i, j)
        END DO
      END DO
    END IF
    IF (fpp%fpp_overload_r4) THEN
      DO j=js,je
        DO i=is,ie
          dpx(i, j) = dpx(i, j) + (fx(i, j)-fx(i+1, j)+(fy(i, j)-fy(i, j&
&           +1)))*rarea(i, j)
        END DO
      END DO
    END IF
!----------------------
! Kinetic Energy Fluxes
!----------------------
! Compute B grid contra-variant components for KE:
    dt5 = 0.5*dt
    dt4 = 0.25*dt
    IF (nested) THEN
      is2 = is
      ie1 = ie + 1
      js2 = js
      je1 = je + 1
    ELSE
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
      IF (2 .LT. js) THEN
        js2 = js
      ELSE
        js2 = 2
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        je1 = je + 1
      ELSE
        je1 = npy - 1
      END IF
    END IF
!!! TO DO: separate versions for nested and for cubed-sphere
    IF (flagstruct%grid_type .LT. 3) THEN
      IF (nested) THEN
        DO j=js2,je1
          DO i=is2,ie1
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j)-(uc(i, j-1)+uc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
        END DO
      ELSE
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb(i, 1) = dt5*(vt(i-1, 1)+vt(i, 1))
          END DO
        END IF
        DO j=js2,je1
          DO i=is2,ie1
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j)-(uc(i, j-1)+uc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
          IF (is .EQ. 1) vb(1, j) = dt4*(-vt(-1, j)+3.*(vt(0, j)+vt(1, j&
&             ))-vt(2, j))
! 2-pt extrapolation from both sides:
          IF (ie + 1 .EQ. npx) vb(npx, j) = dt4*(-vt(npx-2, j)+3.*(vt(&
&             npx-1, j)+vt(npx, j))-vt(npx+1, j))
! 2-pt extrapolation from both sides:
        END DO
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb(i, npy) = dt5*(vt(i-1, npy)+vt(i, npy))
          END DO
        END IF
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          vb(i, j) = dt5*(vc(i-1, j)+vc(i, j))
        END DO
      END DO
    END IF
    IF (hord_mt .EQ. hord_mt_pert) THEN
      CALL YTP_V(is, ie, js, je, isd, ied, jsd, jed, vb, u, v, ub, &
&             hord_mt, gridstruct%dy, gridstruct%rdy, npx, npy, &
&             flagstruct%grid_type, nested)
    ELSE
      CALL YTP_V(is, ie, js, je, isd, ied, jsd, jed, vb, u, v, ub, &
&          hord_mt, gridstruct%dy, gridstruct%rdy, npx, npy, &
&          flagstruct%grid_type, nested)
    END IF
    DO j=js,je+1
      DO i=is,ie+1
        ke(i, j) = vb(i, j)*ub(i, j)
      END DO
    END DO
    IF (flagstruct%grid_type .LT. 3) THEN
      IF (nested) THEN
        DO j=js,je+1
          DO i=is2,ie1
            ub(i, j) = dt5*(uc(i, j-1)+uc(i, j)-(vc(i-1, j)+vc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
        END DO
      ELSE
        IF (is .EQ. 1) THEN
          DO j=js,je+1
! corner values are incorrect
            ub(1, j) = dt5*(ut(1, j-1)+ut(1, j))
          END DO
        END IF
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is2,ie1
! 2-pt extrapolation from both sides:
              ub(i, j) = dt4*(-ut(i, j-2)+3.*(ut(i, j-1)+ut(i, j))-ut(i&
&               , j+1))
            END DO
          ELSE
            DO i=is2,ie1
              ub(i, j) = dt5*(uc(i, j-1)+uc(i, j)-(vc(i-1, j)+vc(i, j))*&
&               cosa(i, j))*rsina(i, j)
            END DO
          END IF
        END DO
        IF (ie + 1 .EQ. npx) THEN
          DO j=js,je+1
! corner values are incorrect
            ub(npx, j) = dt5*(ut(npx, j-1)+ut(npx, j))
          END DO
        END IF
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          ub(i, j) = dt5*(uc(i, j-1)+uc(i, j))
        END DO
      END DO
    END IF
    IF (hord_mt .EQ. hord_mt_pert) THEN
      CALL XTP_U(is, ie, js, je, isd, ied, jsd, jed, ub, u, v, vb, &
&             hord_mt, gridstruct%dx, gridstruct%rdx, npx, npy, &
&             flagstruct%grid_type, nested)
    ELSE
      CALL XTP_U(is, ie, js, je, isd, ied, jsd, jed, ub, u, v, vb, &
&          hord_mt, gridstruct%dx, gridstruct%rdx, npx, npy, &
&          flagstruct%grid_type, nested)
    END IF
    DO j=js,je+1
      DO i=is,ie+1
        ke(i, j) = 0.5*(ke(i, j)+ub(i, j)*vb(i, j))
      END DO
    END DO
!-----------------------------------------
! Fix KE at the 4 corners of the face:
!-----------------------------------------
    IF (.NOT.nested) THEN
      dt6 = dt/6.
      IF (sw_corner) ke(1, 1) = dt6*((ut(1, 1)+ut(1, 0))*u(1, 1)+(vt(1, &
&         1)+vt(0, 1))*v(1, 1)+(ut(1, 1)+vt(1, 1))*u(0, 1))
      IF (se_corner) ke(npx, 1) = dt6*((ut(npx, 1)+ut(npx, 0))*u(npx-1, &
&         1)+(vt(npx, 1)+vt(npx-1, 1))*v(npx, 1)+(ut(npx, 1)-vt(npx-1, 1&
&         ))*u(npx, 1))
!i = npx
      IF (ne_corner) ke(npx, npy) = dt6*((ut(npx, npy)+ut(npx, npy-1))*u&
&         (npx-1, npy)+(vt(npx, npy)+vt(npx-1, npy))*v(npx, npy-1)+(ut(&
&         npx, npy-1)+vt(npx-1, npy))*u(npx, npy))
!i = npx;      j = npy
      IF (nw_corner) ke(1, npy) = dt6*((ut(1, npy)+ut(1, npy-1))*u(1, &
&         npy)+(vt(1, npy)+vt(0, npy))*v(1, npy-1)+(ut(1, npy-1)-vt(1, &
&         npy))*u(0, npy))
!j = npy
    END IF
! Compute vorticity:
    DO j=jsd,jed+1
      DO i=isd,ied
        vt(i, j) = u(i, j)*dx(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied+1
        ut(i, j) = v(i, j)*dy(i, j)
      END DO
    END DO
! wk is "volume-mean" relative vorticity
    DO j=jsd,jed
      DO i=isd,ied
        wk(i, j) = rarea(i, j)*(vt(i, j)-vt(i, j+1)+(ut(i+1, j)-ut(i, j)&
&         ))
      END DO
    END DO
    IF (.NOT.hydrostatic) THEN
      IF (.NOT.flagstruct%do_f3d) THEN
        DO j=js,je
          DO i=is,ie
            w(i, j) = w(i, j)/delp(i, j)
          END DO
        END DO
      END IF
      IF (damp_w .GT. 1.e-5) THEN
        DO j=js,je
          DO i=is,ie
            w(i, j) = w(i, j) + dw(i, j)
          END DO
        END DO
      END IF
    END IF
!-----------------------------
! Compute divergence damping
!-----------------------------
!!  damp = dddmp * da_min_c
!
!   if ( nord==0 ) then
!!         area ~ dxb*dyb*sin(alpha)
!
!      if (nested) then
!
!         do j=js,je+1
!            do i=is-1,ie+1
!               ptc(i,j) = (u(i,j)-0.5*(va(i,j-1)+va(i,j))*cosa_v(i,j))   &
!                    *dyc(i,j)*sina_v(i,j)
!            enddo
!         enddo
!
!         do j=js-1,je+1
!            do i=is2,ie1
!               vort(i,j) = (v(i,j) - 0.5*(ua(i-1,j)+ua(i,j))*cosa_u(i,j))  &
!                    *dxc(i,j)*sina_u(i,j)
!            enddo
!         enddo
!
!      else
!         do j=js,je+1
!
!            if ( (j==1 .or. j==npy)  ) then
!               do i=is-1,ie+1
!                  if (vc(i,j) > 0) then
!                     ptc(i,j) = u(i,j)*dyc(i,j)*sin_sg(i,j-1,4)
!                  else
!                     ptc(i,j) = u(i,j)*dyc(i,j)*sin_sg(i,j,2)
!                  end if
!               enddo
!            else
!               do i=is-1,ie+1
!                  ptc(i,j) = (u(i,j)-0.5*(va(i,j-1)+va(i,j))*cosa_v(i,j))   &
!                       *dyc(i,j)*sina_v(i,j)
!               enddo
!            endif
!         enddo
!
!         do j=js-1,je+1
!            do i=is2,ie1
!               vort(i,j) = (v(i,j) - 0.5*(ua(i-1,j)+ua(i,j))*cosa_u(i,j))  &
!                    *dxc(i,j)*sina_u(i,j)
!            enddo
!            if ( is ==  1 ) then
!               if (uc(1,j) > 0) then
!                  vort(1,  j) = v(1,  j)*dxc(1,  j)*sin_sg(0,j,3)
!               else
!                  vort(1,  j) = v(1,  j)*dxc(1,  j)*sin_sg(1,j,1)
!               end if
!            end if
!            if ( (ie+1)==npx ) then 
!               if (uc(npx,j) > 0) then
!                  vort(npx,j) = v(npx,j)*dxc(npx,j)* & 
!                       sin_sg(npx-1,j,3)
!               else
!                  vort(npx,j) = v(npx,j)*dxc(npx,j)* &
!                       sin_sg(npx,j,1)
!               end if
!            end if
!         enddo
!      endif
!
!      do j=js,je+1
!         do i=is,ie+1
!            delpc(i,j) = vort(i,j-1) - vort(i,j) + ptc(i-1,j) - ptc(i,j)
!         enddo
!      enddo
!
!! Remove the extra term at the corners:
!      if (sw_corner) delpc(1,    1) = delpc(1,    1) - vort(1,    0)
!      if (se_corner) delpc(npx,  1) = delpc(npx,  1) - vort(npx,  0)
!      if (ne_corner) delpc(npx,npy) = delpc(npx,npy) + vort(npx,npy)
!      if (nw_corner) delpc(1,  npy) = delpc(1,  npy) + vort(1,  npy)
!
!      do j=js,je+1
!         do i=is,ie+1
!            delpc(i,j) = gridstruct%rarea_c(i,j)*delpc(i,j)
!                damp = gridstruct%da_min_c*max(d2_bg, min(0.20, dddmp*abs(delpc(i,j)*dt)))
!                vort(i,j) = damp*delpc(i,j)
!                ke(i,j) = ke(i,j) + vort(i,j)
!         enddo
!      enddo
!   else
!!--------------------------
!! Higher order divg damping
!!--------------------------
!     do j=js,je+1
!        do i=is,ie+1
!! Save divergence for external mode filter
!           delpc(i,j) = divg_d(i,j)
!        enddo
!     enddo
!
!     n2 = nord + 1    ! N > 1
!     do n=1,nord
!        nt = nord-n
!
!        fill_c = (nt/=0) .and. (flagstruct%grid_type<3) .and.               &
!                 ( sw_corner .or. se_corner .or. ne_corner .or. nw_corner ) &
!                  .and. .not. nested
!
!        if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=XDir, BGRID=.true.)
!        do j=js-nt,je+1+nt
!           do i=is-1-nt,ie+1+nt
!              vc(i,j) = (divg_d(i+1,j)-divg_d(i,j))*divg_u(i,j)
!           enddo
!        enddo
!
!        if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=YDir, BGRID=.true.)
!        do j=js-1-nt,je+1+nt
!           do i=is-nt,ie+1+nt
!              uc(i,j) = (divg_d(i,j+1)-divg_d(i,j))*divg_v(i,j)
!           enddo
!        enddo
!
!        if ( fill_c ) call fill_corners(vc, uc, npx, npy, VECTOR=.true., DGRID=.true.)
!        do j=js-nt,je+1+nt
!           do i=is-nt,ie+1+nt
!              divg_d(i,j) = uc(i,j-1) - uc(i,j) + vc(i-1,j) - vc(i,j)
!           enddo
!        enddo
!
!! Remove the extra term at the corners:
!        if (sw_corner) divg_d(1,    1) = divg_d(1,    1) - uc(1,    0)
!        if (se_corner) divg_d(npx,  1) = divg_d(npx,  1) - uc(npx,  0)
!        if (ne_corner) divg_d(npx,npy) = divg_d(npx,npy) + uc(npx,npy)
!        if (nw_corner) divg_d(1,  npy) = divg_d(1,  npy) + uc(1,  npy)
!
!     if ( .not. gridstruct%stretched_grid ) then
!        do j=js-nt,je+1+nt
!           do i=is-nt,ie+1+nt
!              divg_d(i,j) = divg_d(i,j)*gridstruct%rarea_c(i,j)
!           enddo
!        enddo
!     endif
!
!     enddo ! n-loop
!
!     if ( dddmp<1.E-5) then
!          vort(:,:) = 0.
!     else
!      if ( flagstruct%grid_type < 3 ) then
!! Interpolate relative vort to cell corners
!          call a2b_ord4(wk, vort, gridstruct, npx, npy, is, ie, js, je, ng, .false.)
!          do j=js,je+1
!             do i=is,ie+1
!! The following is an approxi form of Smagorinsky diffusion
!                vort(i,j) = abs(dt)*sqrt(delpc(i,j)**2 + vort(i,j)**2)
!             enddo
!          enddo
!      else  ! Correct form: works only for doubly preiodic domain
!          call smag_corner(abs(dt), u, v, ua, va, vort, bd, npx, npy, gridstruct, ng)
!      endif
!     endif
!
!     if (gridstruct%stretched_grid ) then
!! Stretched grid with variable damping ~ area
!         dd8 = gridstruct%da_min * d4_bg**n2
!     else
!         dd8 = ( gridstruct%da_min_c*d4_bg )**n2
!     endif
!
!     do j=js,je+1
!        do i=is,ie+1
!           damp2 =  gridstruct%da_min_c*max(d2_bg, min(0.20, dddmp*vort(i,j)))  ! del-2
!           vort(i,j) = damp2*delpc(i,j) + dd8*divg_d(i,j)
!             ke(i,j) = ke(i,j) + vort(i,j)
!        enddo
!     enddo
!
!   endif
    IF (.NOT.split_damp) THEN
      CALL COMPUTE_DIVERGENCE_DAMPING(nord, d2_bg, d4_bg, dddmp, dt, &
&                                  vort, ptc, delpc, ke, u, v, uc, vc, &
&                                  ua, va, divg_d, wk, gridstruct, &
&                                  flagstruct, bd)
    ELSE
      CALL COMPUTE_DIVERGENCE_DAMPING(nord, d2_bg, d4_bg&
&                               , dddmp, dt, vort, ptc, delpc, ke, &
&                               u, v, uc, vc, ua, va, divg_d, wk, &
&                               gridstruct, flagstruct, bd)
    END IF
    IF (d_con .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          ub(i, j) = vort(i, j) - vort(i+1, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          vb(i, j) = vort(i, j) - vort(i, j+1)
        END DO
      END DO
    END IF
! Vorticity transport
    IF (hydrostatic) THEN
      DO j=jsd,jed
        DO i=isd,ied
          vort(i, j) = wk(i, j) + f0(i, j)
        END DO
      END DO
    ELSE IF (flagstruct%do_f3d) THEN
      DO j=jsd,jed
        DO i=isd,ied
          vort(i, j) = wk(i, j) + f0(i, j)*z_rat(i, j)
        END DO
      END DO
    ELSE
      DO j=jsd,jed
        DO i=isd,ied
          vort(i, j) = wk(i, j) + f0(i, j)
        END DO
      END DO
    END IF
    IF (hord_vt .EQ. hord_vt_pert) THEN
      CALL FV_TP_2D(vort, crx_adv, cry_adv, npx, npy, hord_vt, fx, fy&
&                , xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y)
    ELSE
      CALL FV_TP_2D(vort, crx_adv, cry_adv, npx, npy, hord_vt, fx, &
&             fy, xfx_adv, yfx_adv, gridstruct, bd, ra_x, ra_y)
    END IF
    DO j=js,je+1
      DO i=is,ie
        u(i, j) = vt(i, j) + (ke(i, j)-ke(i+1, j)) + fy(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie+1
        v(i, j) = ut(i, j) + (ke(i, j)-ke(i, j+1)) - fx(i, j)
      END DO
    END DO
!--------------------------------------------------------
! damping applied to relative vorticity (wk):
    IF (damp_v .GT. 1.e-5) THEN
      damp4 = (damp_v*gridstruct%da_min_c)**(nord_v+1)
      CALL DEL6_VT_FLUX(nord_v, npx, npy, damp4, wk, vort, ut, vt, &
&                 gridstruct, bd)
    END IF
    IF (d_con .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          ub(i, j) = (ub(i, j)+vt(i, j))*rdx(i, j)
          fy(i, j) = u(i, j)*rdx(i, j)
          gy(i, j) = fy(i, j)*ub(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          vb(i, j) = (vb(i, j)-ut(i, j))*rdy(i, j)
          fx(i, j) = v(i, j)*rdy(i, j)
          gx(i, j) = fx(i, j)*vb(i, j)
        END DO
      END DO
!----------------------------------
! Heating due to damping:
!----------------------------------
      damp = 0.25*d_con
      DO j=js,je
        DO i=is,ie
          u2 = fy(i, j) + fy(i, j+1)
          du2 = ub(i, j) + ub(i, j+1)
          v2 = fx(i, j) + fx(i+1, j)
          dv2 = vb(i, j) + vb(i+1, j)
! Total energy conserving:
! Convert lost KE due to divergence damping to "heat"
          heat_source(i, j) = delp(i, j)*(heat_source(i, j)-damp*rsin2(i&
&           , j)*(ub(i, j)**2+ub(i, j+1)**2+vb(i, j)**2+vb(i+1, j)**2+2.&
&           *(gy(i, j)+gy(i, j+1)+gx(i, j)+gx(i+1, j))-cosa_s(i, j)*(u2*&
&           dv2+v2*du2+du2*dv2)))
        END DO
      END DO
    END IF
! Add diffusive fluxes to the momentum equation:
    IF (damp_v .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          u(i, j) = u(i, j) + vt(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v(i, j) = v(i, j) - ut(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE D_SW
!  Differentiation of del6_vt_flux in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 d
!yn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c 
!dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_cor
!e_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dyn
!amics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_uti
!ls_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz
!_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.sc
!alar_profile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.s
!teepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2
!c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.
!compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solv
!er_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solv
!er nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh
! sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_c
!ore_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mo
!d.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corn
!ers_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_ci
!rcle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: q fy2 d2 fx2
!   with respect to varying inputs: q fy2 d2 fx2
  SUBROUTINE DEL6_VT_FLUX_ADM(nord, npx, npy, damp, q, q_ad, d2, d2_ad, &
&   fx2, fx2_ad, fy2, fy2_ad, gridstruct, bd)
    IMPLICIT NONE
! Del-nord damping for the relative vorticity
! nord must be <= 2
!------------------
! nord = 0:   del-2
! nord = 1:   del-4
! nord = 2:   del-6
!------------------
    INTEGER, INTENT(IN) :: nord, npx, npy
    REAL, INTENT(IN) :: damp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! rel. vorticity ghosted on input
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! Work arrays:
    REAL :: d2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: d2_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2(bd%isd:bd%ied, bd%&
&   jsd:bd%jed+1)
    REAL :: fx2_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2_ad(bd%isd:bd%ied&
&   , bd%jsd:bd%jed+1)
    INTEGER :: i, j, nt, n, i1, i2, j1, j2
    LOGICAL :: nested
    INTEGER :: is, ie, js, je
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
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
    nested = gridstruct%nested
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    i1 = is - 1 - nord
    i2 = ie + 1 + nord
    j1 = js - 1 - nord
    j2 = je + 1 + nord
    IF (nord .GT. 0) THEN
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    IF (nord .GT. 0) THEN
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    IF (nord .GT. 0) THEN
      DO n=1,nord
        nt = nord - n
        ad_from0 = js - nt - 1
        DO j=ad_from0,je+nt+1
          ad_from = is - nt - 1
          i = ie + nt + 2
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from0)
        ad_from2 = js - nt
        DO j=ad_from2,je+nt
          ad_from1 = is - nt
          i = ie + nt + 2
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from1)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from2)
        ad_from4 = js - nt
        DO j=ad_from4,je+nt+1
          ad_from3 = is - nt
          i = ie + nt + 1
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from3)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from4)
      END DO
      DO n=nord,1,-1
        CALL POPINTEGER4(ad_from4)
        CALL POPINTEGER4(ad_to4)
        DO j=ad_to4,ad_from4,-1
          CALL POPINTEGER4(ad_from3)
          CALL POPINTEGER4(ad_to3)
          DO i=ad_to3,ad_from3,-1
            temp_ad3 = gridstruct%del6_u(i, j)*fy2_ad(i, j)
            d2_ad(i, j) = d2_ad(i, j) + temp_ad3
            d2_ad(i, j-1) = d2_ad(i, j-1) - temp_ad3
            fy2_ad(i, j) = 0.0
          END DO
        END DO
        CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 2, nested, bd, &
&                       gridstruct%sw_corner, gridstruct%se_corner, &
&                       gridstruct%nw_corner, gridstruct%ne_corner)
        CALL POPINTEGER4(ad_from2)
        CALL POPINTEGER4(ad_to2)
        DO j=ad_to2,ad_from2,-1
          CALL POPINTEGER4(ad_from1)
          CALL POPINTEGER4(ad_to1)
          DO i=ad_to1,ad_from1,-1
            temp_ad2 = gridstruct%del6_v(i, j)*fx2_ad(i, j)
            d2_ad(i, j) = d2_ad(i, j) + temp_ad2
            d2_ad(i-1, j) = d2_ad(i-1, j) - temp_ad2
            fx2_ad(i, j) = 0.0
          END DO
        END DO
        CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 1, nested, bd, &
&                       gridstruct%sw_corner, gridstruct%se_corner, &
&                       gridstruct%nw_corner, gridstruct%ne_corner)
        CALL POPINTEGER4(ad_from0)
        CALL POPINTEGER4(ad_to0)
        DO j=ad_to0,ad_from0,-1
          CALL POPINTEGER4(ad_from)
          CALL POPINTEGER4(ad_to)
          DO i=ad_to,ad_from,-1
            temp_ad1 = gridstruct%rarea(i, j)*d2_ad(i, j)
            fx2_ad(i, j) = fx2_ad(i, j) + temp_ad1
            fx2_ad(i+1, j) = fx2_ad(i+1, j) - temp_ad1
            fy2_ad(i, j) = fy2_ad(i, j) + temp_ad1
            fy2_ad(i, j+1) = fy2_ad(i, j+1) - temp_ad1
            d2_ad(i, j) = 0.0
          END DO
        END DO
      END DO
    END IF
    DO j=je+nord+1,js-nord,-1
      DO i=ie+nord,is-nord,-1
        temp_ad0 = gridstruct%del6_u(i, j)*fy2_ad(i, j)
        d2_ad(i, j-1) = d2_ad(i, j-1) + temp_ad0
        d2_ad(i, j) = d2_ad(i, j) - temp_ad0
        fy2_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 2, &
&                                      nested, bd, gridstruct%sw_corner&
&                                      , gridstruct%se_corner, &
&                                      gridstruct%nw_corner, gridstruct%&
&                                      ne_corner)
    DO j=je+nord,js-nord,-1
      DO i=ie+nord+1,is-nord,-1
        temp_ad = gridstruct%del6_v(i, j)*fx2_ad(i, j)
        d2_ad(i-1, j) = d2_ad(i-1, j) + temp_ad
        d2_ad(i, j) = d2_ad(i, j) - temp_ad
        fx2_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 1, &
&                                      nested, bd, gridstruct%sw_corner&
&                                      , gridstruct%se_corner, &
&                                      gridstruct%nw_corner, gridstruct%&
&                                      ne_corner)
    DO j=j2,j1,-1
      DO i=i2,i1,-1
        q_ad(i, j) = q_ad(i, j) + damp*d2_ad(i, j)
        d2_ad(i, j) = 0.0
      END DO
    END DO
  END SUBROUTINE DEL6_VT_FLUX_ADM
  SUBROUTINE DEL6_VT_FLUX(nord, npx, npy, damp, q, d2, fx2, fy2, &
&   gridstruct, bd)
    IMPLICIT NONE
! Del-nord damping for the relative vorticity
! nord must be <= 2
!------------------
! nord = 0:   del-2
! nord = 1:   del-4
! nord = 2:   del-6
!------------------
    INTEGER, INTENT(IN) :: nord, npx, npy
    REAL, INTENT(IN) :: damp
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! rel. vorticity ghosted on input
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! Work arrays:
    REAL, INTENT(OUT) :: d2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(OUT) :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2(bd%isd&
&   :bd%ied, bd%jsd:bd%jed+1)
    INTEGER :: i, j, nt, n, i1, i2, j1, j2
    LOGICAL :: nested
    INTEGER :: is, ie, js, je
    nested = gridstruct%nested
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    i1 = is - 1 - nord
    i2 = ie + 1 + nord
    j1 = js - 1 - nord
    j2 = je + 1 + nord
    DO j=j1,j2
      DO i=i1,i2
        d2(i, j) = damp*q(i, j)
      END DO
    END DO
    IF (nord .GT. 0) CALL COPY_CORNERS(d2, npx, npy, 1, nested, bd, &
&                                gridstruct%sw_corner, gridstruct%&
&                                se_corner, gridstruct%nw_corner, &
&                                gridstruct%ne_corner)
    DO j=js-nord,je+nord
      DO i=is-nord,ie+nord+1
        fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i-1, j)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) CALL COPY_CORNERS(d2, npx, npy, 2, nested, bd, &
&                                gridstruct%sw_corner, gridstruct%&
&                                se_corner, gridstruct%nw_corner, &
&                                gridstruct%ne_corner)
    DO j=js-nord,je+nord+1
      DO i=is-nord,ie+nord
        fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j-1)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) THEN
      DO n=1,nord
        nt = nord - n
        DO j=js-nt-1,je+nt+1
          DO i=is-nt-1,ie+nt+1
            d2(i, j) = (fx2(i, j)-fx2(i+1, j)+(fy2(i, j)-fy2(i, j+1)))*&
&             gridstruct%rarea(i, j)
          END DO
        END DO
        CALL COPY_CORNERS(d2, npx, npy, 1, nested, bd, gridstruct%&
&                   sw_corner, gridstruct%se_corner, gridstruct%&
&                   nw_corner, gridstruct%ne_corner)
        DO j=js-nt,je+nt
          DO i=is-nt,ie+nt+1
            fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i, j)-d2(i-1, j))
          END DO
        END DO
        CALL COPY_CORNERS(d2, npx, npy, 2, nested, bd, gridstruct%&
&                   sw_corner, gridstruct%se_corner, gridstruct%&
&                   nw_corner, gridstruct%ne_corner)
        DO j=js-nt,je+nt+1
          DO i=is-nt,ie+nt
            fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j)-d2(i, j-1))
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE DEL6_VT_FLUX
!  Differentiation of divergence_corner in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b
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
!   gradient     of useful results: u v ua va divg_d
!   with respect to varying inputs: u v ua va divg_d
  SUBROUTINE DIVERGENCE_CORNER_FWD(u, v, ua, va, divg_d, gridstruct, &
&   flagstruct, bd)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1) :: divg_d
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! local
    REAL :: uf(bd%is-2:bd%ie+2, bd%js-1:bd%je+2)
    REAL :: vf(bd%is-1:bd%ie+2, bd%js-2:bd%je+2)
    INTEGER :: i, j
    INTEGER :: is2, ie1
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg, cos_sg
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc
    INTEGER :: is, ie, js, je
    INTEGER :: npx, npy
    LOGICAL :: nested
    INTRINSIC MAX
    INTRINSIC MIN

    uf = 0.0
    vf = 0.0
    is2 = 0
    ie1 = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    npx = 0
    npy = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    sin_sg => gridstruct%sin_sg
    cos_sg => gridstruct%cos_sg
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    IF (nested) THEN
      CALL PUSHCONTROL(1,0)
      is2 = is
      ie1 = ie + 1
    ELSE
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        CALL PUSHCONTROL(1,1)
        ie1 = ie + 1
      ELSE
        CALL PUSHCONTROL(1,1)
        ie1 = npx - 1
      END IF
    END IF
    IF (flagstruct%grid_type .EQ. 4) THEN
      DO j=js-1,je+2
        DO i=is-2,ie+2
          uf(i, j) = u(i, j)*dyc(i, j)
        END DO
      END DO
      DO j=js-2,je+2
        DO i=is-1,ie+2
          vf(i, j) = v(i, j)*dxc(i, j)
        END DO
      END DO
      DO j=js-1,je+2
        DO i=is-1,ie+2
          divg_d(i, j) = gridstruct%rarea_c(i, j)*(vf(i, j-1)-vf(i, j)+(&
&           uf(i-1, j)-uf(i, j)))
        END DO
      END DO
      CALL PUSHINTEGER(je)
      CALL PUSHINTEGER(is)
      CALL PUSHINTEGER(ie)
      !CALL PUSHPOINTER8(C_LOC(dyc))
      !CALL PUSHPOINTER8(C_LOC(dxc))
      CALL PUSHINTEGER(js)
      CALL PUSHCONTROL(1,0)
    ELSE
!     9---4---8
!     |       |
!     1   5   3
!     |       |
!     6---2---7
      DO j=js,je+1
        IF (j .EQ. 1 .OR. j .EQ. npy) THEN
          DO i=is-1,ie+1
            uf(i, j) = u(i, j)*dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+sin_sg(i&
&             , j, 2))
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          DO i=is-1,ie+1
            uf(i, j) = (u(i, j)-0.25*(va(i, j-1)+va(i, j))*(cos_sg(i, j-&
&             1, 4)+cos_sg(i, j, 2)))*dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+&
&             sin_sg(i, j, 2))
          END DO
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      DO j=js-1,je+1
        DO i=is2,ie1
          vf(i, j) = (v(i, j)-0.25*(ua(i-1, j)+ua(i, j))*(cos_sg(i-1, j&
&           , 3)+cos_sg(i, j, 1)))*dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+&
&           sin_sg(i, j, 1))
        END DO
        IF (is .EQ. 1) THEN
          vf(1, j) = v(1, j)*dxc(1, j)*0.5*(sin_sg(0, j, 3)+sin_sg(1, j&
&           , 1))
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (ie + 1 .EQ. npx) THEN
          vf(npx, j) = v(npx, j)*dxc(npx, j)*0.5*(sin_sg(npx-1, j, 3)+&
&           sin_sg(npx, j, 1))
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      DO j=js,je+1
        DO i=is,ie+1
          divg_d(i, j) = vf(i, j-1) - vf(i, j) + (uf(i-1, j)-uf(i, j))
        END DO
      END DO
! Remove the extra term at the corners:
      IF (gridstruct%sw_corner) THEN
        divg_d(1, 1) = divg_d(1, 1) - vf(1, 0)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%se_corner) THEN
        divg_d(npx, 1) = divg_d(npx, 1) - vf(npx, 0)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%ne_corner) THEN
        divg_d(npx, npy) = divg_d(npx, npy) + vf(npx, npy)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%nw_corner) THEN
        divg_d(1, npy) = divg_d(1, npy) + vf(1, npy)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          divg_d(i, j) = gridstruct%rarea_c(i, j)*divg_d(i, j)
        END DO
      END DO
      CALL PUSHINTEGER(je)
      CALL PUSHINTEGER(ie1)
      CALL PUSHINTEGER(is)
      CALL PUSHINTEGER(ie)
      !CALL PUSHPOINTER8(C_LOC(dyc))
      !CALL PUSHPOINTER8(C_LOC(sin_sg))
      CALL PUSHINTEGER(is2)
      !CALL PUSHPOINTER8(C_LOC(dxc))
      !CALL PUSHPOINTER8(C_LOC(cos_sg))
      CALL PUSHINTEGER(npy)
      CALL PUSHINTEGER(npx)
      CALL PUSHINTEGER(js)
      CALL PUSHCONTROL(1,1)
    END IF
  END SUBROUTINE DIVERGENCE_CORNER_FWD
!  Differentiation of divergence_corner in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2
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
!   gradient     of useful results: u v ua va divg_d
!   with respect to varying inputs: u v ua va divg_d
  SUBROUTINE DIVERGENCE_CORNER_BWD(u, u_ad, v, v_ad, ua, ua_ad, va, &
&   va_ad, divg_d, divg_d_ad, gridstruct, flagstruct, bd)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1) :: u_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed) :: v_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ua_ad, va_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1) :: divg_d
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1) :: divg_d_ad
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    REAL :: uf(bd%is-2:bd%ie+2, bd%js-1:bd%je+2)
    REAL :: uf_ad(bd%is-2:bd%ie+2, bd%js-1:bd%je+2)
    REAL :: vf(bd%is-1:bd%ie+2, bd%js-2:bd%je+2)
    REAL :: vf_ad(bd%is-1:bd%ie+2, bd%js-2:bd%je+2)
    INTEGER :: i, j
    INTEGER :: is2, ie1
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg, cos_sg
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc
    INTEGER :: is, ie, js, je
    INTEGER :: npx, npy
    LOGICAL :: nested
    INTRINSIC MAX
    INTRINSIC MIN
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    INTEGER :: branch
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_divergence_corner

    uf = 0.0
    vf = 0.0
    is2 = 0
    ie1 = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    npx = 0
    npy = 0
    branch = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    sin_sg => gridstruct%sin_sg
    cos_sg => gridstruct%cos_sg
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(js)
      !CALL POPPOINTER8(cptr)
      dxc => gridstruct%dxc ! (/unknown_shape_in_divergence_corner/)&
!&               )
      !CALL POPPOINTER8(cptr)
      dyc => gridstruct%dyc ! (/unknown_shape_in_divergence_corner/)&
!&               )
      CALL POPINTEGER(ie)
      CALL POPINTEGER(is)
      CALL POPINTEGER(je)
      uf_ad = 0.0
      vf_ad = 0.0
      DO j=je+2,js-1,-1
        DO i=ie+2,is-1,-1
          temp_ad = gridstruct%rarea_c(i, j)*divg_d_ad(i, j)
          vf_ad(i, j-1) = vf_ad(i, j-1) + temp_ad
          vf_ad(i, j) = vf_ad(i, j) - temp_ad
          uf_ad(i-1, j) = uf_ad(i-1, j) + temp_ad
          uf_ad(i, j) = uf_ad(i, j) - temp_ad
          divg_d_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+2,js-2,-1
        DO i=ie+2,is-1,-1
          v_ad(i, j) = v_ad(i, j) + dxc(i, j)*vf_ad(i, j)
          vf_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+2,js-1,-1
        DO i=ie+2,is-2,-1
          u_ad(i, j) = u_ad(i, j) + dyc(i, j)*uf_ad(i, j)
          uf_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      CALL POPINTEGER(js)
      CALL POPINTEGER(npx)
      CALL POPINTEGER(npy)
      !CALL POPPOINTER8(cptr)
      cos_sg => gridstruct%cos_sg ! (/&
!&                unknown_shape_in_divergence_corner/))
      !CALL POPPOINTER8(cptr)
      dxc => gridstruct%dxc ! (/unknown_shape_in_divergence_corner/)&
!&               )
      CALL POPINTEGER(is2)
      !CALL POPPOINTER8(cptr)
      sin_sg => gridstruct%sin_sg ! (/&
!&                unknown_shape_in_divergence_corner/))
      !CALL POPPOINTER8(cptr)
      dyc => gridstruct%dyc ! (/unknown_shape_in_divergence_corner/)&
!&               )
      CALL POPINTEGER(ie)
      CALL POPINTEGER(is)
      CALL POPINTEGER(ie1)
      CALL POPINTEGER(je)
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          divg_d_ad(i, j) = gridstruct%rarea_c(i, j)*divg_d_ad(i, j)
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        vf_ad = 0.0
      ELSE
        npy = flagstruct%npy
        vf_ad = 0.0
        vf_ad(1, npy) = vf_ad(1, npy) + divg_d_ad(1, npy)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        npx = flagstruct%npx
        vf_ad(npx, npy) = vf_ad(npx, npy) + divg_d_ad(npx, npy)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) vf_ad(npx, 0) = vf_ad(npx, 0) - divg_d_ad(npx, &
&         1)
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) vf_ad(1, 0) = vf_ad(1, 0) - divg_d_ad(1, 1)
      uf_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          vf_ad(i, j-1) = vf_ad(i, j-1) + divg_d_ad(i, j)
          vf_ad(i, j) = vf_ad(i, j) - divg_d_ad(i, j)
          uf_ad(i-1, j) = uf_ad(i-1, j) + divg_d_ad(i, j)
          uf_ad(i, j) = uf_ad(i, j) - divg_d_ad(i, j)
          divg_d_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+1,js-1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          v_ad(npx, j) = v_ad(npx, j) + dxc(npx, j)*(sin_sg(npx-1, j, 3)&
&           +sin_sg(npx, j, 1))*0.5*vf_ad(npx, j)
          vf_ad(npx, j) = 0.0
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          v_ad(1, j) = v_ad(1, j) + dxc(1, j)*(sin_sg(0, j, 3)+sin_sg(1&
&           , j, 1))*0.5*vf_ad(1, j)
          vf_ad(1, j) = 0.0
        END IF
        DO i=ie1,is2,-1
          temp_ad2 = dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+sin_sg(i, j, 1))*&
&           vf_ad(i, j)
          temp_ad3 = -((cos_sg(i-1, j, 3)+cos_sg(i, j, 1))*0.25*temp_ad2&
&           )
          v_ad(i, j) = v_ad(i, j) + temp_ad2
          ua_ad(i-1, j) = ua_ad(i-1, j) + temp_ad3
          ua_ad(i, j) = ua_ad(i, j) + temp_ad3
          vf_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+1,js,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO i=ie+1,is-1,-1
            temp_ad0 = dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+sin_sg(i, j, 2))&
&             *uf_ad(i, j)
            temp_ad1 = -((cos_sg(i, j-1, 4)+cos_sg(i, j, 2))*0.25*&
&             temp_ad0)
            u_ad(i, j) = u_ad(i, j) + temp_ad0
            va_ad(i, j-1) = va_ad(i, j-1) + temp_ad1
            va_ad(i, j) = va_ad(i, j) + temp_ad1
            uf_ad(i, j) = 0.0
          END DO
        ELSE
          DO i=ie+1,is-1,-1
            u_ad(i, j) = u_ad(i, j) + dyc(i, j)*(sin_sg(i, j-1, 4)+&
&             sin_sg(i, j, 2))*0.5*uf_ad(i, j)
            uf_ad(i, j) = 0.0
          END DO
        END IF
      END DO
    END IF
    CALL POPCONTROL(1,branch)
  END SUBROUTINE DIVERGENCE_CORNER_BWD
  SUBROUTINE DIVERGENCE_CORNER(u, v, ua, va, divg_d, gridstruct, &
&   flagstruct, bd)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   divg_d
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! local
    REAL :: uf(bd%is-2:bd%ie+2, bd%js-1:bd%je+2)
    REAL :: vf(bd%is-1:bd%ie+2, bd%js-2:bd%je+2)
    INTEGER :: i, j
    INTEGER :: is2, ie1
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg, cos_sg
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc
    INTEGER :: is, ie, js, je
    INTEGER :: npx, npy
    LOGICAL :: nested
    INTRINSIC MAX
    INTRINSIC MIN
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    sin_sg => gridstruct%sin_sg
    cos_sg => gridstruct%cos_sg
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    IF (nested) THEN
      is2 = is
      ie1 = ie + 1
    ELSE
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
    END IF
    IF (flagstruct%grid_type .EQ. 4) THEN
      DO j=js-1,je+2
        DO i=is-2,ie+2
          uf(i, j) = u(i, j)*dyc(i, j)
        END DO
      END DO
      DO j=js-2,je+2
        DO i=is-1,ie+2
          vf(i, j) = v(i, j)*dxc(i, j)
        END DO
      END DO
      DO j=js-1,je+2
        DO i=is-1,ie+2
          divg_d(i, j) = gridstruct%rarea_c(i, j)*(vf(i, j-1)-vf(i, j)+(&
&           uf(i-1, j)-uf(i, j)))
        END DO
      END DO
    ELSE
!     9---4---8
!     |       |
!     1   5   3
!     |       |
!     6---2---7
      DO j=js,je+1
        IF (j .EQ. 1 .OR. j .EQ. npy) THEN
          DO i=is-1,ie+1
            uf(i, j) = u(i, j)*dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+sin_sg(i&
&             , j, 2))
          END DO
        ELSE
          DO i=is-1,ie+1
            uf(i, j) = (u(i, j)-0.25*(va(i, j-1)+va(i, j))*(cos_sg(i, j-&
&             1, 4)+cos_sg(i, j, 2)))*dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+&
&             sin_sg(i, j, 2))
          END DO
        END IF
      END DO
      DO j=js-1,je+1
        DO i=is2,ie1
          vf(i, j) = (v(i, j)-0.25*(ua(i-1, j)+ua(i, j))*(cos_sg(i-1, j&
&           , 3)+cos_sg(i, j, 1)))*dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+&
&           sin_sg(i, j, 1))
        END DO
        IF (is .EQ. 1) vf(1, j) = v(1, j)*dxc(1, j)*0.5*(sin_sg(0, j, 3)&
&           +sin_sg(1, j, 1))
        IF (ie + 1 .EQ. npx) vf(npx, j) = v(npx, j)*dxc(npx, j)*0.5*(&
&           sin_sg(npx-1, j, 3)+sin_sg(npx, j, 1))
      END DO
      DO j=js,je+1
        DO i=is,ie+1
          divg_d(i, j) = vf(i, j-1) - vf(i, j) + (uf(i-1, j)-uf(i, j))
        END DO
      END DO
! Remove the extra term at the corners:
      IF (gridstruct%sw_corner) divg_d(1, 1) = divg_d(1, 1) - vf(1, 0)
      IF (gridstruct%se_corner) divg_d(npx, 1) = divg_d(npx, 1) - vf(npx&
&         , 0)
      IF (gridstruct%ne_corner) divg_d(npx, npy) = divg_d(npx, npy) + vf&
&         (npx, npy)
      IF (gridstruct%nw_corner) divg_d(1, npy) = divg_d(1, npy) + vf(1, &
&         npy)
      DO j=js,je+1
        DO i=is,ie+1
          divg_d(i, j) = gridstruct%rarea_c(i, j)*divg_d(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE DIVERGENCE_CORNER
!  Differentiation of divergence_corner_nest in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4_f
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
!   gradient     of useful results: u v ua va divg_d
!   with respect to varying inputs: u v ua va
  SUBROUTINE DIVERGENCE_CORNER_NEST_FWD(u, v, ua, va, divg_d, gridstruct&
&   , flagstruct, bd)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
!!$       !Edges
!!$
!!$       !West, East
!!$       do j=jsd+1,jed
!!$          divg_d(isd  ,j) = (vf(isd,j-1) - vf(isd,j) + uf(isd,j) - uf(isd+1,j))*rarea_c(isd,j)
!!$          divg_d(ied+1,j) = (vf(ied+1,j-1) - vf(ied+1,j) + uf(ied-1,j) - uf(ied,j))*rarea_c(ied,j)
!!$       end do
!!$
!!$       !North, South
!!$       do i=isd+1,ied
!!$          divg_d(i,jsd  ) = (vf(i,jsd) - vf(i,jsd+1) + uf(i-1,jsd) - uf(i,jsd))*rarea_c(i,jsd)
!!$          divg_d(i,jed+1) = (vf(i,jed-1) - vf(i,jed) + uf(i-1,jed+1) - uf(i,jed+1))*rarea_c(i,jed)
!!$       end do
!!$
!!$       !Corners (just use next corner value)
!!$       divg_d(isd,jsd)   = divg_d(isd+1,jsd+1)
!!$       divg_d(isd,jed+1) = divg_d(isd+1,jed)
!!$       divg_d(ied+1,jsd)   = divg_d(ied,jsd+1)
!!$       divg_d(ied+1,jed+1) = divg_d(ied,jed)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1) :: divg_d
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! local
    REAL :: uf(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: vf(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    INTEGER :: i, j
    REAL, DIMENSION(:, :), POINTER :: rarea_c
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg, cos_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: npx, npy
    LOGICAL :: nested

    uf = 0.0
    vf = 0.0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    npx = 0
    npy = 0

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    rarea_c => gridstruct%rarea_c
    sin_sg => gridstruct%sin_sg
    cos_sg => gridstruct%cos_sg
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    divg_d = 1.e25
    IF (flagstruct%grid_type .EQ. 4) THEN
      DO j=jsd,jed
        DO i=isd,ied
          uf(i, j) = u(i, j)*dyc(i, j)
        END DO
      END DO
      DO j=jsd,jed
        DO i=isd,ied
          vf(i, j) = v(i, j)*dxc(i, j)
        END DO
      END DO
      DO j=jsd+1,jed
        DO i=isd+1,ied
          divg_d(i, j) = rarea_c(i, j)*(vf(i, j-1)-vf(i, j)+(uf(i-1, j)-&
&           uf(i, j)))
        END DO
      END DO
      CALL PUSHINTEGER(jed)
      CALL PUSHINTEGER(isd)
      !CALL PUSHPOINTER8(C_LOC(dyc))
      CALL PUSHINTEGER(ied)
      CALL PUSHINTEGER(jsd)
      !CALL PUSHPOINTER8(C_LOC(rarea_c))
      !CALL PUSHPOINTER8(C_LOC(dxc))
      CALL PUSHCONTROL(1,0)
    ELSE
      DO j=jsd+1,jed
        DO i=isd,ied
          uf(i, j) = (u(i, j)-0.25*(va(i, j-1)+va(i, j))*(cos_sg(i, j-1&
&           , 4)+cos_sg(i, j, 2)))*dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+&
&           sin_sg(i, j, 2))
        END DO
      END DO
      DO j=jsd,jed
        DO i=isd+1,ied
          vf(i, j) = (v(i, j)-0.25*(ua(i-1, j)+ua(i, j))*(cos_sg(i-1, j&
&           , 3)+cos_sg(i, j, 1)))*dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+&
&           sin_sg(i, j, 1))
        END DO
      END DO
      DO j=jsd+1,jed
        DO i=isd+1,ied
          divg_d(i, j) = (vf(i, j-1)-vf(i, j)+(uf(i-1, j)-uf(i, j)))*&
&           rarea_c(i, j)
        END DO
      END DO
      CALL PUSHINTEGER(jed)
      CALL PUSHINTEGER(isd)
      !CALL PUSHPOINTER8(C_LOC(dyc))
      !CALL PUSHPOINTER8(C_LOC(sin_sg))
      CALL PUSHINTEGER(ied)
      CALL PUSHINTEGER(jsd)
      !CALL PUSHPOINTER8(C_LOC(rarea_c))
      !CALL PUSHPOINTER8(C_LOC(dxc))
      !CALL PUSHPOINTER8(C_LOC(cos_sg))
      CALL PUSHCONTROL(1,1)
    END IF
  END SUBROUTINE DIVERGENCE_CORNER_NEST_FWD
!  Differentiation of divergence_corner_nest in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4_
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
!   gradient     of useful results: u v ua va divg_d
!   with respect to varying inputs: u v ua va
  SUBROUTINE DIVERGENCE_CORNER_NEST_BWD(u, u_ad, v, v_ad, ua, ua_ad, va&
&   , va_ad, divg_d, divg_d_ad, gridstruct, flagstruct, bd)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
!!$       !Edges
!!$
!!$       !West, East
!!$       do j=jsd+1,jed
!!$          divg_d(isd  ,j) = (vf(isd,j-1) - vf(isd,j) + uf(isd,j) - uf(isd+1,j))*rarea_c(isd,j)
!!$          divg_d(ied+1,j) = (vf(ied+1,j-1) - vf(ied+1,j) + uf(ied-1,j) - uf(ied,j))*rarea_c(ied,j)
!!$       end do
!!$
!!$       !North, South
!!$       do i=isd+1,ied
!!$          divg_d(i,jsd  ) = (vf(i,jsd) - vf(i,jsd+1) + uf(i-1,jsd) - uf(i,jsd))*rarea_c(i,jsd)
!!$          divg_d(i,jed+1) = (vf(i,jed-1) - vf(i,jed) + uf(i-1,jed+1) - uf(i,jed+1))*rarea_c(i,jed)
!!$       end do
!!$
!!$       !Corners (just use next corner value)
!!$       divg_d(isd,jsd)   = divg_d(isd+1,jsd+1)
!!$       divg_d(isd,jed+1) = divg_d(isd+1,jed)
!!$       divg_d(ied+1,jsd)   = divg_d(ied,jsd+1)
!!$       divg_d(ied+1,jed+1) = divg_d(ied,jed)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1) :: u_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed) :: v_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ua_ad, va_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1) :: divg_d
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1) :: divg_d_ad
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    REAL :: uf(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: uf_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: vf(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vf_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    INTEGER :: i, j
    REAL, DIMENSION(:, :), POINTER :: rarea_c
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg, cos_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: npx, npy
    LOGICAL :: nested
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    INTEGER :: branch
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_divergence_corner_nest

    uf = 0.0
    vf = 0.0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    npx = 0
    npy = 0
    branch = 0

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    rarea_c => gridstruct%rarea_c
    sin_sg => gridstruct%sin_sg
    cos_sg => gridstruct%cos_sg
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      !CALL POPPOINTER8(cptr)
      dxc => gridstruct%dxc ! (/&
!&                unknown_shape_in_divergence_corner_nest/))
      !CALL POPPOINTER8(cptr)
      rarea_c => gridstruct%rarea_c ! (/&
!&                unknown_shape_in_divergence_corner_nest/))
      CALL POPINTEGER(jsd)
      CALL POPINTEGER(ied)
      !CALL POPPOINTER8(cptr)
      dyc => gridstruct%dyc ! (/&
!&                unknown_shape_in_divergence_corner_nest/))
      CALL POPINTEGER(isd)
      CALL POPINTEGER(jed)
      uf_ad = 0.0
      vf_ad = 0.0
      DO j=jed,jsd+1,-1
        DO i=ied,isd+1,-1
          temp_ad = rarea_c(i, j)*divg_d_ad(i, j)
          vf_ad(i, j-1) = vf_ad(i, j-1) + temp_ad
          vf_ad(i, j) = vf_ad(i, j) - temp_ad
          uf_ad(i-1, j) = uf_ad(i-1, j) + temp_ad
          uf_ad(i, j) = uf_ad(i, j) - temp_ad
          divg_d_ad(i, j) = 0.0
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ied,isd,-1
          v_ad(i, j) = v_ad(i, j) + dxc(i, j)*vf_ad(i, j)
          vf_ad(i, j) = 0.0
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ied,isd,-1
          u_ad(i, j) = u_ad(i, j) + dyc(i, j)*uf_ad(i, j)
          uf_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      !CALL POPPOINTER8(cptr)
      cos_sg => gridstruct%cos_sg ! (/&
!&                unknown_shape_in_divergence_corner_nest/))
      !CALL POPPOINTER8(cptr)
      dxc => gridstruct%dxc ! (/&
!&                unknown_shape_in_divergence_corner_nest/))
      !CALL POPPOINTER8(cptr)
      rarea_c => gridstruct%rarea_c ! (/&
!&                unknown_shape_in_divergence_corner_nest/))
      CALL POPINTEGER(jsd)
      CALL POPINTEGER(ied)
      !CALL POPPOINTER8(cptr)
      sin_sg => gridstruct%sin_sg ! (/&
!&                unknown_shape_in_divergence_corner_nest/))
      !CALL POPPOINTER8(cptr)
      dyc => gridstruct%dyc ! (/&
!&                unknown_shape_in_divergence_corner_nest/))
      CALL POPINTEGER(isd)
      CALL POPINTEGER(jed)
      uf_ad = 0.0
      vf_ad = 0.0
      DO j=jed,jsd+1,-1
        DO i=ied,isd+1,-1
          temp_ad4 = rarea_c(i, j)*divg_d_ad(i, j)
          vf_ad(i, j-1) = vf_ad(i, j-1) + temp_ad4
          vf_ad(i, j) = vf_ad(i, j) - temp_ad4
          uf_ad(i-1, j) = uf_ad(i-1, j) + temp_ad4
          uf_ad(i, j) = uf_ad(i, j) - temp_ad4
          divg_d_ad(i, j) = 0.0
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ied,isd+1,-1
          temp_ad2 = dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+sin_sg(i, j, 1))*&
&           vf_ad(i, j)
          temp_ad3 = -((cos_sg(i-1, j, 3)+cos_sg(i, j, 1))*0.25*temp_ad2&
&           )
          v_ad(i, j) = v_ad(i, j) + temp_ad2
          ua_ad(i-1, j) = ua_ad(i-1, j) + temp_ad3
          ua_ad(i, j) = ua_ad(i, j) + temp_ad3
          vf_ad(i, j) = 0.0
        END DO
      END DO
      DO j=jed,jsd+1,-1
        DO i=ied,isd,-1
          temp_ad0 = dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+sin_sg(i, j, 2))*&
&           uf_ad(i, j)
          temp_ad1 = -((cos_sg(i, j-1, 4)+cos_sg(i, j, 2))*0.25*temp_ad0&
&           )
          u_ad(i, j) = u_ad(i, j) + temp_ad0
          va_ad(i, j-1) = va_ad(i, j-1) + temp_ad1
          va_ad(i, j) = va_ad(i, j) + temp_ad1
          uf_ad(i, j) = 0.0
        END DO
      END DO
    END IF
  END SUBROUTINE DIVERGENCE_CORNER_NEST_BWD
  SUBROUTINE DIVERGENCE_CORNER_NEST(u, v, ua, va, divg_d, gridstruct, &
&   flagstruct, bd)
    IMPLICIT NONE
!!$       !Edges
!!$
!!$       !West, East
!!$       do j=jsd+1,jed
!!$          divg_d(isd  ,j) = (vf(isd,j-1) - vf(isd,j) + uf(isd,j) - uf(isd+1,j))*rarea_c(isd,j)
!!$          divg_d(ied+1,j) = (vf(ied+1,j-1) - vf(ied+1,j) + uf(ied-1,j) - uf(ied,j))*rarea_c(ied,j)
!!$       end do
!!$
!!$       !North, South
!!$       do i=isd+1,ied
!!$          divg_d(i,jsd  ) = (vf(i,jsd) - vf(i,jsd+1) + uf(i-1,jsd) - uf(i,jsd))*rarea_c(i,jsd)
!!$          divg_d(i,jed+1) = (vf(i,jed-1) - vf(i,jed) + uf(i-1,jed+1) - uf(i,jed+1))*rarea_c(i,jed)
!!$       end do
!!$
!!$       !Corners (just use next corner value)
!!$       divg_d(isd,jsd)   = divg_d(isd+1,jsd+1)
!!$       divg_d(isd,jed+1) = divg_d(isd+1,jed)
!!$       divg_d(ied+1,jsd)   = divg_d(ied,jsd+1)
!!$       divg_d(ied+1,jed+1) = divg_d(ied,jed)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   divg_d
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! local
    REAL :: uf(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: vf(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    INTEGER :: i, j
    REAL, DIMENSION(:, :), POINTER :: rarea_c
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg, cos_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: npx, npy
    LOGICAL :: nested
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    rarea_c => gridstruct%rarea_c
    sin_sg => gridstruct%sin_sg
    cos_sg => gridstruct%cos_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    sina_u => gridstruct%sina_u
    sina_v => gridstruct%sina_v
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    divg_d = 1.e25
    IF (flagstruct%grid_type .EQ. 4) THEN
      DO j=jsd,jed
        DO i=isd,ied
          uf(i, j) = u(i, j)*dyc(i, j)
        END DO
      END DO
      DO j=jsd,jed
        DO i=isd,ied
          vf(i, j) = v(i, j)*dxc(i, j)
        END DO
      END DO
      DO j=jsd+1,jed
        DO i=isd+1,ied
          divg_d(i, j) = rarea_c(i, j)*(vf(i, j-1)-vf(i, j)+(uf(i-1, j)-&
&           uf(i, j)))
        END DO
      END DO
    ELSE
      DO j=jsd+1,jed
        DO i=isd,ied
          uf(i, j) = (u(i, j)-0.25*(va(i, j-1)+va(i, j))*(cos_sg(i, j-1&
&           , 4)+cos_sg(i, j, 2)))*dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+&
&           sin_sg(i, j, 2))
        END DO
      END DO
      DO j=jsd,jed
        DO i=isd+1,ied
          vf(i, j) = (v(i, j)-0.25*(ua(i-1, j)+ua(i, j))*(cos_sg(i-1, j&
&           , 3)+cos_sg(i, j, 1)))*dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+&
&           sin_sg(i, j, 1))
        END DO
      END DO
      DO j=jsd+1,jed
        DO i=isd+1,ied
          divg_d(i, j) = (vf(i, j-1)-vf(i, j)+(uf(i-1, j)-uf(i, j)))*&
&           rarea_c(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE DIVERGENCE_CORNER_NEST
  SUBROUTINE SMAG_CORNER(dt, u, v, ua, va, smag_c, bd, npx, npy, &
&   gridstruct, ng)
    IMPLICIT NONE
! Compute the Tension_Shear strain at cell corners for Smagorinsky diffusion
!!!  work only if (grid_type==4)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: npx, npy, ng
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: smag_c
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! local
    REAL :: ut(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!  work array
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: sh(bd%isd:bd%ied, bd%jsd:bd%jed)
    INTEGER :: i, j
    INTEGER :: is2, ie1
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc, dx, dy, rarea, rarea_c
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SQRT
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    dx => gridstruct%dx
    dy => gridstruct%dy
    rarea => gridstruct%rarea
    rarea_c => gridstruct%rarea_c
    IF (2 .LT. is) THEN
      is2 = is
    ELSE
      is2 = 2
    END IF
    IF (npx - 1 .GT. ie + 1) THEN
      ie1 = ie + 1
    ELSE
      ie1 = npx - 1
    END IF
! Smag = sqrt [ T**2 + S**2 ]:  unit = 1/s
! where T = du/dx - dv/dy;   S = du/dy + dv/dx
! Compute tension strain at corners:
    DO j=js,je+1
      DO i=is-1,ie+1
        ut(i, j) = u(i, j)*dyc(i, j)
      END DO
    END DO
    DO j=js-1,je+1
      DO i=is,ie+1
        vt(i, j) = v(i, j)*dxc(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie+1
        smag_c(i, j) = rarea_c(i, j)*(vt(i, j-1)-vt(i, j)+(ut(i, j)-ut(i&
&         -1, j)))
      END DO
    END DO
! Fix the corners?? if grid_type /= 4
! Compute shear strain:
    DO j=jsd,jed+1
      DO i=isd,ied
        vt(i, j) = u(i, j)*dx(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied+1
        ut(i, j) = v(i, j)*dy(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied
        wk(i, j) = rarea(i, j)*(vt(i, j)-vt(i, j+1)+(ut(i, j)-ut(i+1, j)&
&         ))
      END DO
    END DO
    CALL A2B_ORD4(wk, sh, gridstruct, npx, npy, is, ie, js, je, ng, &
&           .false.)
    DO j=js,je+1
      DO i=is,ie+1
        smag_c(i, j) = dt*SQRT(sh(i, j)**2+smag_c(i, j)**2)
      END DO
    END DO
  END SUBROUTINE SMAG_CORNER
!  Differentiation of xtp_u in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_core
!_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn_cor
!e_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_mod.g
!eopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynamics_m
!od.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mod.
!c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.ma
!p_scalar_fb fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scalar_pr
!ofile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz f
!v_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_setup
! fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.compute
!_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_c nh
!_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver nh_u
!tils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw_cor
!e_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core_mod
!.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.compu
!te_divergence_damping_fb sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corners 
!tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circle_di
!st sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: flux u c
!   with respect to varying inputs: flux u c
  SUBROUTINE XTP_U_ADM(is, ie, js, je, isd, ied, jsd, jed, c, c_ad, u, &
&   u_ad, v, flux, flux_ad, iord, dx, rdx, npx, npy, grid_type, nested)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL :: u_ad(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL :: c_ad(is:ie+1, js:je+1)
    REAL :: flux(is:ie+1, js:je+1)
    REAL :: flux_ad(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: rdx(isd:ied, jsd:jed+1)
    INTEGER, INTENT(IN) :: iord, npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
! Local
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    REAL, DIMENSION(is-1:ie+1) :: bl_ad, br_ad, b0_ad
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: fx0_ad(is:ie+1)
    REAL :: al(is-1:ie+2), dm(is-2:ie+2)
    REAL :: al_ad(is-1:ie+2), dm_ad(is-2:ie+2)
    REAL :: dq(is-3:ie+2)
    REAL :: dq_ad(is-3:ie+2)
    REAL :: dl, dr, xt, pmp, lac, cfl
    REAL :: xt_ad, pmp_ad, lac_ad, cfl_ad
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: pmp_1_ad, lac_1_ad, pmp_2_ad, lac_2_ad
    REAL :: x0, x1, x0l, x0r
    REAL :: x0l_ad, x0r_ad
    INTEGER :: i, j
    INTEGER :: is3, ie3
    INTEGER :: is2, ie2
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min1_ad
    REAL :: min2
    REAL :: min2_ad
    REAL :: abs0
    REAL :: min3
    REAL :: min3_ad
    REAL :: min4
    REAL :: min4_ad
    REAL :: min5
    REAL :: min5_ad
    REAL :: abs1
    REAL :: abs2
    REAL :: abs3
    REAL :: abs4
    REAL :: max1
    REAL :: max1_ad
    REAL :: min6
    REAL :: min6_ad
    REAL :: abs5
    REAL :: abs6
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: x2_ad
    REAL :: y1_ad
    REAL :: x3_ad
    REAL :: y2_ad
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: x4_ad
    REAL :: y3_ad
    REAL :: z1_ad
    REAL :: x5_ad
    REAL :: y4_ad
    REAL :: x6_ad
    REAL :: y5_ad
    REAL :: x7_ad
    REAL :: y12_ad
    REAL :: y6_ad
    REAL :: x8_ad
    REAL :: y13_ad
    REAL :: y7_ad
    REAL :: x9_ad
    REAL :: y14_ad
    REAL :: y8_ad
    REAL :: x10_ad
    REAL :: y15_ad
    REAL :: y9_ad
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: x11_ad
    REAL :: y16_ad
    REAL :: y10_ad
    REAL :: x12_ad
    REAL :: y17_ad
    REAL :: y11_ad
    REAL :: temp
    REAL :: temp_ad13
    REAL :: temp_ad14
    INTEGER :: branch
    REAL :: x12
    REAL :: x11
    REAL :: x10
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: y17
    REAL :: y16
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
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
    IF (nested .OR. grid_type .GT. 3) THEN
      CALL PUSHCONTROL1B(0)
      is3 = is - 1
      ie3 = ie + 1
    ELSE
      IF (3 .LT. is - 1) THEN
        is3 = is - 1
      ELSE
        is3 = 3
      END IF
      IF (npx - 3 .GT. ie + 1) THEN
        CALL PUSHCONTROL1B(1)
        ie3 = ie + 1
      ELSE
        CALL PUSHCONTROL1B(1)
        ie3 = npx - 3
      END IF
    END IF
    IF (iord .EQ. 1) THEN
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
            flux_ad(i, j) = 0.0
          ELSE
            u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
            flux_ad(i, j) = 0.0
          END IF
        END DO
      END DO
    ELSE IF (iord .LT. 8) THEN
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
      DO j=js,je+1
        DO i=is3,ie3+1
          al(i) = p1*(u(i-1, j)+u(i, j)) + p2*(u(i-2, j)+u(i+1, j))
        END DO
        DO i=is3,ie3
          CALL PUSHREALARRAY_ADM(bl(i))
          bl(i) = al(i) - u(i, j)
          CALL PUSHREALARRAY_ADM(br(i))
          br(i) = al(i+1) - u(i, j)
        END DO
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            xt = c3*u(1, j) + c2*u(2, j) + c1*u(3, j)
            CALL PUSHREALARRAY_ADM(br(1))
            br(1) = xt - u(1, j)
            CALL PUSHREALARRAY_ADM(bl(2))
            bl(2) = xt - u(2, j)
            CALL PUSHREALARRAY_ADM(br(2))
            br(2) = al(3) - u(2, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! out
              CALL PUSHREALARRAY_ADM(bl(0))
              bl(0) = 0.
! edge
              CALL PUSHREALARRAY_ADM(br(0))
              br(0) = 0.
! edge
              CALL PUSHREALARRAY_ADM(bl(1))
              bl(1) = 0.
! in
              CALL PUSHREALARRAY_ADM(br(1))
              br(1) = 0.
              CALL PUSHCONTROL2B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(bl(0))
              bl(0) = c1*u(-2, j) + c2*u(-1, j) + c3*u(0, j) - u(0, j)
              xt = 0.5*(((2.*dx(0, j)+dx(-1, j))*u(0, j)-dx(0, j)*u(-1, &
&               j))/(dx(0, j)+dx(-1, j))+((2.*dx(1, j)+dx(2, j))*u(1, j)&
&               -dx(1, j)*u(2, j))/(dx(1, j)+dx(2, j)))
              CALL PUSHREALARRAY_ADM(br(0))
              br(0) = xt - u(0, j)
              CALL PUSHREALARRAY_ADM(bl(1))
              bl(1) = xt - u(1, j)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE
            CALL PUSHCONTROL2B(2)
          END IF
!       call pert_ppm(1, u(2,j), bl(2), br(2), -1)
          IF (ie + 1 .EQ. npx) THEN
            CALL PUSHREALARRAY_ADM(bl(npx-2))
            bl(npx-2) = al(npx-2) - u(npx-2, j)
            xt = c1*u(npx-3, j) + c2*u(npx-2, j) + c3*u(npx-1, j)
            CALL PUSHREALARRAY_ADM(br(npx-2))
            br(npx-2) = xt - u(npx-2, j)
            CALL PUSHREALARRAY_ADM(bl(npx-1))
            bl(npx-1) = xt - u(npx-1, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! in
              CALL PUSHREALARRAY_ADM(bl(npx-1))
              bl(npx-1) = 0.
! edge
              CALL PUSHREALARRAY_ADM(br(npx-1))
              br(npx-1) = 0.
! edge
              CALL PUSHREALARRAY_ADM(bl(npx))
              bl(npx) = 0.
! out
              CALL PUSHREALARRAY_ADM(br(npx))
              br(npx) = 0.
              CALL PUSHCONTROL2B(3)
            ELSE
              xt = 0.5*(((2.*dx(npx-1, j)+dx(npx-2, j))*u(npx-1, j)-dx(&
&               npx-1, j)*u(npx-2, j))/(dx(npx-1, j)+dx(npx-2, j))+((2.*&
&               dx(npx, j)+dx(npx+1, j))*u(npx, j)-dx(npx, j)*u(npx+1, j&
&               ))/(dx(npx, j)+dx(npx+1, j)))
              CALL PUSHREALARRAY_ADM(br(npx-1))
              br(npx-1) = xt - u(npx-1, j)
              CALL PUSHREALARRAY_ADM(bl(npx))
              bl(npx) = xt - u(npx, j)
              CALL PUSHREALARRAY_ADM(br(npx))
              br(npx) = c3*u(npx, j) + c2*u(npx+1, j) + c1*u(npx+2, j) -&
&               u(npx, j)
              CALL PUSHCONTROL2B(2)
            END IF
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(0)
        END IF
!       call pert_ppm(1, u(npx-2,j), bl(npx-2), br(npx-2), -1)
        DO i=is-1,ie+1
          CALL PUSHREALARRAY_ADM(b0(i))
          b0(i) = bl(i) + br(i)
        END DO
        IF (iord .EQ. 2) THEN
! Perfectly linear
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
          CALL PUSHCONTROL2B(3)
        ELSE IF (iord .EQ. 3) THEN
          DO i=is-1,ie+1
            IF (b0(i) .GE. 0.) THEN
              x0 = b0(i)
            ELSE
              x0 = -b0(i)
            END IF
            IF (bl(i) - br(i) .GE. 0.) THEN
              x1 = bl(i) - br(i)
            ELSE
              x1 = -(bl(i)-br(i))
            END IF
            smt5(i) = x0 .LT. x1
            smt6(i) = 3.*x0 .LT. x1
          END DO
          DO i=is,ie+1
            CALL PUSHREALARRAY_ADM(fx0(i))
            fx0(i) = 0.
          END DO
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdx(i-1, j)
              IF (smt6(i-1) .OR. smt5(i)) THEN
                CALL PUSHREALARRAY_ADM(fx0(i))
                fx0(i) = br(i-1) - cfl*b0(i-1)
                CALL PUSHCONTROL2B(0)
              ELSE IF (smt5(i-1)) THEN
                IF (bl(i-1) .GE. 0.) THEN
                  x2 = bl(i-1)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  x2 = -bl(i-1)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (br(i-1) .GE. 0.) THEN
                  y1 = br(i-1)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  y1 = -br(i-1)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (x2 .GT. y1) THEN
                  CALL PUSHREALARRAY_ADM(min1)
                  min1 = y1
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHREALARRAY_ADM(min1)
                  min1 = x2
                  CALL PUSHCONTROL1B(1)
                END IF
                CALL PUSHREALARRAY_ADM(fx0(i))
                fx0(i) = SIGN(min1, br(i-1))
                CALL PUSHCONTROL2B(1)
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
              CALL PUSHCONTROL1B(1)
            ELSE
              cfl = c(i, j)*rdx(i, j)
              IF (smt6(i) .OR. smt5(i-1)) THEN
                CALL PUSHREALARRAY_ADM(fx0(i))
                fx0(i) = bl(i) + cfl*b0(i)
                CALL PUSHCONTROL2B(0)
              ELSE IF (smt5(i)) THEN
                IF (bl(i) .GE. 0.) THEN
                  x3 = bl(i)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  x3 = -bl(i)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (br(i) .GE. 0.) THEN
                  y2 = br(i)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  y2 = -br(i)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (x3 .GT. y2) THEN
                  CALL PUSHREALARRAY_ADM(min2)
                  min2 = y2
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHREALARRAY_ADM(min2)
                  min2 = x3
                  CALL PUSHCONTROL1B(1)
                END IF
                CALL PUSHREALARRAY_ADM(fx0(i))
                fx0(i) = SIGN(min2, bl(i))
                CALL PUSHCONTROL2B(1)
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
          CALL PUSHCONTROL2B(2)
        ELSE IF (iord .EQ. 4) THEN
! more damp than ord5 but less damp than ord6
          DO i=is-1,ie+1
            IF (b0(i) .GE. 0.) THEN
              x0 = b0(i)
            ELSE
              x0 = -b0(i)
            END IF
            IF (bl(i) - br(i) .GE. 0.) THEN
              x1 = bl(i) - br(i)
            ELSE
              x1 = -(bl(i)-br(i))
            END IF
            smt5(i) = x0 .LT. x1
! if smt6 =.T. --> smt5=.T.
            smt6(i) = 3.*x0 .LT. x1
          END DO
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              IF (smt6(i-1) .OR. smt5(i)) THEN
                CALL PUSHCONTROL2B(3)
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
            ELSE IF (smt6(i) .OR. smt5(i-1)) THEN
              CALL PUSHCONTROL2B(1)
            ELSE
              CALL PUSHCONTROL2B(0)
            END IF
          END DO
          CALL PUSHCONTROL2B(1)
        ELSE
!  iord=5,6,7
          IF (iord .EQ. 5) THEN
            CALL PUSHCONTROL1B(1)
            DO i=is-1,ie+1
              smt5(i) = bl(i)*br(i) .LT. 0.
            END DO
          ELSE
            DO i=is-1,ie+1
              IF (3.*b0(i) .GE. 0.) THEN
                abs0 = 3.*b0(i)
              ELSE
                abs0 = -(3.*b0(i))
              END IF
              IF (bl(i) - br(i) .GE. 0.) THEN
                abs4 = bl(i) - br(i)
              ELSE
                abs4 = -(bl(i)-br(i))
              END IF
              smt5(i) = abs0 .LT. abs4
            END DO
            CALL PUSHCONTROL1B(0)
          END IF
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdx(i-1, j)
              CALL PUSHREALARRAY_ADM(fx0(i))
              fx0(i) = (1.-cfl)*(br(i-1)-cfl*b0(i-1))
              CALL PUSHCONTROL1B(0)
            ELSE
              cfl = c(i, j)*rdx(i, j)
              CALL PUSHREALARRAY_ADM(fx0(i))
              fx0(i) = (1.+cfl)*(bl(i)+cfl*b0(i))
              CALL PUSHCONTROL1B(1)
            END IF
            IF (smt5(i-1) .OR. smt5(i)) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
          CALL PUSHCONTROL2B(0)
        END IF
      END DO
      al_ad = 0.0
      bl_ad = 0.0
      br_ad = 0.0
      b0_ad = 0.0
      fx0_ad = 0.0
      DO j=je+1,js,-1
        CALL POPCONTROL2B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            DO i=ie+1,is,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) fx0_ad(i) = fx0_ad(i) + flux_ad(i, j)
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
                flux_ad(i, j) = 0.0
                cfl = c(i, j)*rdx(i-1, j)
                CALL POPREALARRAY_ADM(fx0(i))
                temp_ad7 = (1.-cfl)*fx0_ad(i)
                cfl_ad = -(b0(i-1)*temp_ad7) - (br(i-1)-cfl*b0(i-1))*&
&                 fx0_ad(i)
                br_ad(i-1) = br_ad(i-1) + temp_ad7
                b0_ad(i-1) = b0_ad(i-1) - cfl*temp_ad7
                fx0_ad(i) = 0.0
                c_ad(i, j) = c_ad(i, j) + rdx(i-1, j)*cfl_ad
              ELSE
                u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
                flux_ad(i, j) = 0.0
                cfl = c(i, j)*rdx(i, j)
                CALL POPREALARRAY_ADM(fx0(i))
                temp_ad8 = (cfl+1.)*fx0_ad(i)
                cfl_ad = b0(i)*temp_ad8 + (bl(i)+cfl*b0(i))*fx0_ad(i)
                bl_ad(i) = bl_ad(i) + temp_ad8
                b0_ad(i) = b0_ad(i) + cfl*temp_ad8
                fx0_ad(i) = 0.0
                c_ad(i, j) = c_ad(i, j) + rdx(i, j)*cfl_ad
              END IF
            END DO
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) i = is - 2
          ELSE
            DO i=ie+1,is,-1
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
                  flux_ad(i, j) = 0.0
                ELSE
                  cfl = c(i, j)*rdx(i, j)
                  temp_ad6 = (cfl+1.)*flux_ad(i, j)
                  u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
                  cfl_ad = b0(i)*temp_ad6 + (bl(i)+cfl*b0(i))*flux_ad(i&
&                   , j)
                  bl_ad(i) = bl_ad(i) + temp_ad6
                  b0_ad(i) = b0_ad(i) + cfl*temp_ad6
                  flux_ad(i, j) = 0.0
                  c_ad(i, j) = c_ad(i, j) + rdx(i, j)*cfl_ad
                END IF
              ELSE IF (branch .EQ. 2) THEN
                u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
                flux_ad(i, j) = 0.0
              ELSE
                cfl = c(i, j)*rdx(i-1, j)
                temp_ad5 = (1.-cfl)*flux_ad(i, j)
                u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
                cfl_ad = -(b0(i-1)*temp_ad5) - (br(i-1)-cfl*b0(i-1))*&
&                 flux_ad(i, j)
                br_ad(i-1) = br_ad(i-1) + temp_ad5
                b0_ad(i-1) = b0_ad(i-1) - cfl*temp_ad5
                flux_ad(i, j) = 0.0
                c_ad(i, j) = c_ad(i, j) + rdx(i-1, j)*cfl_ad
              END IF
            END DO
            i = is - 2
          END IF
        ELSE IF (branch .EQ. 2) THEN
          DO i=ie+1,is,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              cfl = c(i, j)*rdx(i, j)
              u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
              cfl_ad = fx0(i)*flux_ad(i, j)
              fx0_ad(i) = fx0_ad(i) + (cfl+1.)*flux_ad(i, j)
              flux_ad(i, j) = 0.0
              CALL POPCONTROL2B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(fx0(i))
                bl_ad(i) = bl_ad(i) + fx0_ad(i)
                cfl_ad = cfl_ad + b0(i)*fx0_ad(i)
                b0_ad(i) = b0_ad(i) + cfl*fx0_ad(i)
                fx0_ad(i) = 0.0
              ELSE IF (branch .EQ. 1) THEN
                CALL POPREALARRAY_ADM(fx0(i))
                min2_ad = SIGN(1.d0, min2*bl(i))*fx0_ad(i)
                fx0_ad(i) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(min2)
                  y2_ad = min2_ad
                  x3_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(min2)
                  x3_ad = min2_ad
                  y2_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  br_ad(i) = br_ad(i) + y2_ad
                ELSE
                  br_ad(i) = br_ad(i) - y2_ad
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  bl_ad(i) = bl_ad(i) + x3_ad
                ELSE
                  bl_ad(i) = bl_ad(i) - x3_ad
                END IF
              END IF
              c_ad(i, j) = c_ad(i, j) + rdx(i, j)*cfl_ad
            ELSE
              cfl = c(i, j)*rdx(i-1, j)
              u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
              cfl_ad = -(fx0(i)*flux_ad(i, j))
              fx0_ad(i) = fx0_ad(i) + (1.-cfl)*flux_ad(i, j)
              flux_ad(i, j) = 0.0
              CALL POPCONTROL2B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(fx0(i))
                br_ad(i-1) = br_ad(i-1) + fx0_ad(i)
                cfl_ad = cfl_ad - b0(i-1)*fx0_ad(i)
                b0_ad(i-1) = b0_ad(i-1) - cfl*fx0_ad(i)
                fx0_ad(i) = 0.0
              ELSE IF (branch .EQ. 1) THEN
                CALL POPREALARRAY_ADM(fx0(i))
                min1_ad = SIGN(1.d0, min1*br(i-1))*fx0_ad(i)
                fx0_ad(i) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(min1)
                  y1_ad = min1_ad
                  x2_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(min1)
                  x2_ad = min1_ad
                  y1_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  br_ad(i-1) = br_ad(i-1) + y1_ad
                ELSE
                  br_ad(i-1) = br_ad(i-1) - y1_ad
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  bl_ad(i-1) = bl_ad(i-1) + x2_ad
                ELSE
                  bl_ad(i-1) = bl_ad(i-1) - x2_ad
                END IF
              END IF
              c_ad(i, j) = c_ad(i, j) + rdx(i-1, j)*cfl_ad
            END IF
          END DO
          DO i=ie+1,is,-1
            CALL POPREALARRAY_ADM(fx0(i))
            fx0_ad(i) = 0.0
          END DO
          i = is - 2
        ELSE
          DO i=ie+1,is,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              cfl = c(i, j)*rdx(i, j)
              temp_ad4 = (cfl+1.)*flux_ad(i, j)
              u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
              cfl_ad = b0(i)*temp_ad4 + (bl(i)+cfl*b0(i))*flux_ad(i, j)
              bl_ad(i) = bl_ad(i) + temp_ad4
              b0_ad(i) = b0_ad(i) + cfl*temp_ad4
              flux_ad(i, j) = 0.0
              c_ad(i, j) = c_ad(i, j) + rdx(i, j)*cfl_ad
            ELSE
              cfl = c(i, j)*rdx(i-1, j)
              temp_ad3 = (1.-cfl)*flux_ad(i, j)
              u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
              cfl_ad = -(b0(i-1)*temp_ad3) - (br(i-1)-cfl*b0(i-1))*&
&               flux_ad(i, j)
              br_ad(i-1) = br_ad(i-1) + temp_ad3
              b0_ad(i-1) = b0_ad(i-1) - cfl*temp_ad3
              flux_ad(i, j) = 0.0
              c_ad(i, j) = c_ad(i, j) + rdx(i-1, j)*cfl_ad
            END IF
          END DO
        END IF
        DO i=ie+1,is-1,-1
          CALL POPREALARRAY_ADM(b0(i))
          bl_ad(i) = bl_ad(i) + b0_ad(i)
          br_ad(i) = br_ad(i) + b0_ad(i)
          b0_ad(i) = 0.0
        END DO
        CALL POPCONTROL2B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) GOTO 100
        ELSE
          IF (branch .EQ. 2) THEN
            CALL POPREALARRAY_ADM(br(npx))
            u_ad(npx, j) = u_ad(npx, j) + (c3-1.0)*br_ad(npx)
            u_ad(npx+1, j) = u_ad(npx+1, j) + c2*br_ad(npx)
            u_ad(npx+2, j) = u_ad(npx+2, j) + c1*br_ad(npx)
            br_ad(npx) = 0.0
            CALL POPREALARRAY_ADM(bl(npx))
            xt_ad = br_ad(npx-1) + bl_ad(npx)
            u_ad(npx, j) = u_ad(npx, j) - bl_ad(npx)
            bl_ad(npx) = 0.0
            CALL POPREALARRAY_ADM(br(npx-1))
            temp_ad1 = 0.5*xt_ad/(dx(npx-1, j)+dx(npx-2, j))
            u_ad(npx-1, j) = u_ad(npx-1, j) + (dx(npx-1, j)*2.+dx(npx-2&
&             , j))*temp_ad1 - br_ad(npx-1)
            br_ad(npx-1) = 0.0
            temp_ad2 = 0.5*xt_ad/(dx(npx, j)+dx(npx+1, j))
            u_ad(npx-2, j) = u_ad(npx-2, j) - dx(npx-1, j)*temp_ad1
            u_ad(npx, j) = u_ad(npx, j) + (dx(npx, j)*2.+dx(npx+1, j))*&
&             temp_ad2
            u_ad(npx+1, j) = u_ad(npx+1, j) - dx(npx, j)*temp_ad2
          ELSE
            CALL POPREALARRAY_ADM(br(npx))
            br_ad(npx) = 0.0
            CALL POPREALARRAY_ADM(bl(npx))
            bl_ad(npx) = 0.0
            CALL POPREALARRAY_ADM(br(npx-1))
            br_ad(npx-1) = 0.0
            CALL POPREALARRAY_ADM(bl(npx-1))
            bl_ad(npx-1) = 0.0
          END IF
          CALL POPREALARRAY_ADM(bl(npx-1))
          xt_ad = br_ad(npx-2) + bl_ad(npx-1)
          u_ad(npx-1, j) = u_ad(npx-1, j) - bl_ad(npx-1)
          bl_ad(npx-1) = 0.0
          CALL POPREALARRAY_ADM(br(npx-2))
          u_ad(npx-2, j) = u_ad(npx-2, j) - br_ad(npx-2)
          br_ad(npx-2) = 0.0
          u_ad(npx-3, j) = u_ad(npx-3, j) + c1*xt_ad
          u_ad(npx-2, j) = u_ad(npx-2, j) + c2*xt_ad
          u_ad(npx-1, j) = u_ad(npx-1, j) + c3*xt_ad
          CALL POPREALARRAY_ADM(bl(npx-2))
          al_ad(npx-2) = al_ad(npx-2) + bl_ad(npx-2)
          u_ad(npx-2, j) = u_ad(npx-2, j) - bl_ad(npx-2)
          bl_ad(npx-2) = 0.0
        END IF
        CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY_ADM(br(1))
          br_ad(1) = 0.0
          CALL POPREALARRAY_ADM(bl(1))
          bl_ad(1) = 0.0
          CALL POPREALARRAY_ADM(br(0))
          br_ad(0) = 0.0
          CALL POPREALARRAY_ADM(bl(0))
          bl_ad(0) = 0.0
        ELSE IF (branch .EQ. 1) THEN
          CALL POPREALARRAY_ADM(bl(1))
          xt_ad = br_ad(0) + bl_ad(1)
          u_ad(1, j) = u_ad(1, j) - bl_ad(1)
          bl_ad(1) = 0.0
          CALL POPREALARRAY_ADM(br(0))
          temp_ad = 0.5*xt_ad/(dx(0, j)+dx(-1, j))
          u_ad(0, j) = u_ad(0, j) + (dx(0, j)*2.+dx(-1, j))*temp_ad - &
&           br_ad(0)
          br_ad(0) = 0.0
          temp_ad0 = 0.5*xt_ad/(dx(1, j)+dx(2, j))
          u_ad(-1, j) = u_ad(-1, j) - dx(0, j)*temp_ad
          u_ad(1, j) = u_ad(1, j) + (dx(1, j)*2.+dx(2, j))*temp_ad0
          u_ad(2, j) = u_ad(2, j) - dx(1, j)*temp_ad0
          CALL POPREALARRAY_ADM(bl(0))
          u_ad(-2, j) = u_ad(-2, j) + c1*bl_ad(0)
          u_ad(-1, j) = u_ad(-1, j) + c2*bl_ad(0)
          u_ad(0, j) = u_ad(0, j) + (c3-1.0)*bl_ad(0)
          bl_ad(0) = 0.0
        ELSE
          GOTO 100
        END IF
        CALL POPREALARRAY_ADM(br(2))
        al_ad(3) = al_ad(3) + br_ad(2)
        u_ad(2, j) = u_ad(2, j) - bl_ad(2) - br_ad(2)
        br_ad(2) = 0.0
        CALL POPREALARRAY_ADM(bl(2))
        xt_ad = br_ad(1) + bl_ad(2)
        bl_ad(2) = 0.0
        CALL POPREALARRAY_ADM(br(1))
        u_ad(1, j) = u_ad(1, j) + c3*xt_ad - br_ad(1)
        br_ad(1) = 0.0
        u_ad(2, j) = u_ad(2, j) + c2*xt_ad
        u_ad(3, j) = u_ad(3, j) + c1*xt_ad
 100    DO i=ie3,is3,-1
          CALL POPREALARRAY_ADM(br(i))
          al_ad(i+1) = al_ad(i+1) + br_ad(i)
          u_ad(i, j) = u_ad(i, j) - bl_ad(i) - br_ad(i)
          br_ad(i) = 0.0
          CALL POPREALARRAY_ADM(bl(i))
          al_ad(i) = al_ad(i) + bl_ad(i)
          bl_ad(i) = 0.0
        END DO
        DO i=ie3+1,is3,-1
          u_ad(i-1, j) = u_ad(i-1, j) + p1*al_ad(i)
          u_ad(i, j) = u_ad(i, j) + p1*al_ad(i)
          u_ad(i-2, j) = u_ad(i-2, j) + p2*al_ad(i)
          u_ad(i+1, j) = u_ad(i+1, j) + p2*al_ad(i)
          al_ad(i) = 0.0
        END DO
      END DO
    ELSE
! iord = 8, 9, 10, 11
      DO j=js,je+1
        DO i=is-2,ie+2
          CALL PUSHREALARRAY_ADM(xt)
          xt = 0.25*(u(i+1, j)-u(i-1, j))
          IF (xt .GE. 0.) THEN
            x4 = xt
            CALL PUSHCONTROL1B(0)
          ELSE
            x4 = -xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (u(i-1, j) .LT. u(i, j)) THEN
            IF (u(i, j) .LT. u(i+1, j)) THEN
              max1 = u(i+1, j)
              CALL PUSHCONTROL2B(0)
            ELSE
              max1 = u(i, j)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (u(i-1, j) .LT. u(i+1, j)) THEN
            max1 = u(i+1, j)
            CALL PUSHCONTROL2B(2)
          ELSE
            max1 = u(i-1, j)
            CALL PUSHCONTROL2B(3)
          END IF
          y3 = max1 - u(i, j)
          IF (u(i-1, j) .GT. u(i, j)) THEN
            IF (u(i, j) .GT. u(i+1, j)) THEN
              min6 = u(i+1, j)
              CALL PUSHCONTROL2B(0)
            ELSE
              min6 = u(i, j)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (u(i-1, j) .GT. u(i+1, j)) THEN
            min6 = u(i+1, j)
            CALL PUSHCONTROL2B(2)
          ELSE
            min6 = u(i-1, j)
            CALL PUSHCONTROL2B(3)
          END IF
          z1 = u(i, j) - min6
          IF (x4 .GT. y3) THEN
            IF (y3 .GT. z1) THEN
              CALL PUSHREALARRAY_ADM(min3)
              min3 = z1
              CALL PUSHCONTROL2B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(min3)
              min3 = y3
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (x4 .GT. z1) THEN
            CALL PUSHREALARRAY_ADM(min3)
            min3 = z1
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHREALARRAY_ADM(min3)
            min3 = x4
            CALL PUSHCONTROL2B(3)
          END IF
          dm(i) = SIGN(min3, xt)
        END DO
        DO i=is-3,ie+2
          dq(i) = u(i+1, j) - u(i, j)
        END DO
        IF (grid_type .LT. 3) THEN
          DO i=is3,ie3+1
            al(i) = 0.5*(u(i-1, j)+u(i, j)) + r3*(dm(i-1)-dm(i))
          END DO
! Perturbation form:
          IF (iord .EQ. 8) THEN
            DO i=is3,ie3
              CALL PUSHREALARRAY_ADM(xt)
              xt = 2.*dm(i)
              IF (xt .GE. 0.) THEN
                x5 = xt
                CALL PUSHCONTROL1B(0)
              ELSE
                x5 = -xt
                CALL PUSHCONTROL1B(1)
              END IF
              IF (al(i) - u(i, j) .GE. 0.) THEN
                y4 = al(i) - u(i, j)
                CALL PUSHCONTROL1B(0)
              ELSE
                y4 = -(al(i)-u(i, j))
                CALL PUSHCONTROL1B(1)
              END IF
              IF (x5 .GT. y4) THEN
                CALL PUSHREALARRAY_ADM(min4)
                min4 = y4
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(min4)
                min4 = x5
                CALL PUSHCONTROL1B(1)
              END IF
              CALL PUSHREALARRAY_ADM(bl(i))
              bl(i) = -SIGN(min4, xt)
              IF (xt .GE. 0.) THEN
                x6 = xt
                CALL PUSHCONTROL1B(0)
              ELSE
                x6 = -xt
                CALL PUSHCONTROL1B(1)
              END IF
              IF (al(i+1) - u(i, j) .GE. 0.) THEN
                y5 = al(i+1) - u(i, j)
                CALL PUSHCONTROL1B(0)
              ELSE
                y5 = -(al(i+1)-u(i, j))
                CALL PUSHCONTROL1B(1)
              END IF
              IF (x6 .GT. y5) THEN
                CALL PUSHREALARRAY_ADM(min5)
                min5 = y5
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(min5)
                min5 = x6
                CALL PUSHCONTROL1B(1)
              END IF
              CALL PUSHREALARRAY_ADM(br(i))
              br(i) = SIGN(min5, xt)
            END DO
            CALL PUSHCONTROL2B(0)
          ELSE IF (iord .EQ. 9) THEN
            DO i=is3,ie3
              pmp_1 = -(2.*dq(i))
              lac_1 = pmp_1 + 1.5*dq(i+1)
              IF (0. .LT. pmp_1) THEN
                IF (pmp_1 .LT. lac_1) THEN
                  x7 = lac_1
                  CALL PUSHCONTROL2B(0)
                ELSE
                  x7 = pmp_1
                  CALL PUSHCONTROL2B(1)
                END IF
              ELSE IF (0. .LT. lac_1) THEN
                x7 = lac_1
                CALL PUSHCONTROL2B(2)
              ELSE
                CALL PUSHCONTROL2B(3)
                x7 = 0.
              END IF
              IF (0. .GT. pmp_1) THEN
                IF (pmp_1 .GT. lac_1) THEN
                  y12 = lac_1
                  CALL PUSHCONTROL2B(0)
                ELSE
                  y12 = pmp_1
                  CALL PUSHCONTROL2B(1)
                END IF
              ELSE IF (0. .GT. lac_1) THEN
                y12 = lac_1
                CALL PUSHCONTROL2B(2)
              ELSE
                y12 = 0.
                CALL PUSHCONTROL2B(3)
              END IF
              IF (al(i) - u(i, j) .LT. y12) THEN
                y6 = y12
                CALL PUSHCONTROL1B(0)
              ELSE
                y6 = al(i) - u(i, j)
                CALL PUSHCONTROL1B(1)
              END IF
              IF (x7 .GT. y6) THEN
                CALL PUSHREALARRAY_ADM(bl(i))
                bl(i) = y6
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(bl(i))
                bl(i) = x7
                CALL PUSHCONTROL1B(1)
              END IF
              pmp_2 = 2.*dq(i-1)
              lac_2 = pmp_2 - 1.5*dq(i-2)
              IF (0. .LT. pmp_2) THEN
                IF (pmp_2 .LT. lac_2) THEN
                  x8 = lac_2
                  CALL PUSHCONTROL2B(0)
                ELSE
                  x8 = pmp_2
                  CALL PUSHCONTROL2B(1)
                END IF
              ELSE IF (0. .LT. lac_2) THEN
                x8 = lac_2
                CALL PUSHCONTROL2B(2)
              ELSE
                CALL PUSHCONTROL2B(3)
                x8 = 0.
              END IF
              IF (0. .GT. pmp_2) THEN
                IF (pmp_2 .GT. lac_2) THEN
                  y13 = lac_2
                  CALL PUSHCONTROL2B(0)
                ELSE
                  y13 = pmp_2
                  CALL PUSHCONTROL2B(1)
                END IF
              ELSE IF (0. .GT. lac_2) THEN
                y13 = lac_2
                CALL PUSHCONTROL2B(2)
              ELSE
                y13 = 0.
                CALL PUSHCONTROL2B(3)
              END IF
              IF (al(i+1) - u(i, j) .LT. y13) THEN
                y7 = y13
                CALL PUSHCONTROL1B(0)
              ELSE
                y7 = al(i+1) - u(i, j)
                CALL PUSHCONTROL1B(1)
              END IF
              IF (x8 .GT. y7) THEN
                CALL PUSHREALARRAY_ADM(br(i))
                br(i) = y7
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(br(i))
                br(i) = x8
                CALL PUSHCONTROL1B(1)
              END IF
            END DO
            CALL PUSHCONTROL2B(1)
          ELSE IF (iord .EQ. 10) THEN
            DO i=is3,ie3
              CALL PUSHREALARRAY_ADM(bl(i))
              bl(i) = al(i) - u(i, j)
              CALL PUSHREALARRAY_ADM(br(i))
              br(i) = al(i+1) - u(i, j)
              IF (dm(i) .GE. 0.) THEN
                abs1 = dm(i)
              ELSE
                abs1 = -dm(i)
              END IF
!             if ( abs(dm(i-1))+abs(dm(i))+abs(dm(i+1)) < near_zero ) then
              IF (abs1 .LT. near_zero) THEN
                IF (dm(i-1) .GE. 0.) THEN
                  abs2 = dm(i-1)
                ELSE
                  abs2 = -dm(i-1)
                END IF
                IF (dm(i+1) .GE. 0.) THEN
                  abs5 = dm(i+1)
                ELSE
                  abs5 = -dm(i+1)
                END IF
                IF (abs2 + abs5 .LT. near_zero) THEN
! 2-delta-x structure detected within 3 cells
                  bl(i) = 0.
                  br(i) = 0.
                  CALL PUSHCONTROL3B(4)
                ELSE
                  CALL PUSHCONTROL3B(3)
                END IF
              ELSE
                IF (3.*(bl(i)+br(i)) .GE. 0.) THEN
                  abs3 = 3.*(bl(i)+br(i))
                ELSE
                  abs3 = -(3.*(bl(i)+br(i)))
                END IF
                IF (bl(i) - br(i) .GE. 0.) THEN
                  abs6 = bl(i) - br(i)
                ELSE
                  abs6 = -(bl(i)-br(i))
                END IF
                IF (abs3 .GT. abs6) THEN
                  pmp_1 = -(2.*dq(i))
                  lac_1 = pmp_1 + 1.5*dq(i+1)
                  IF (0. .LT. pmp_1) THEN
                    IF (pmp_1 .LT. lac_1) THEN
                      x9 = lac_1
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      x9 = pmp_1
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (0. .LT. lac_1) THEN
                    x9 = lac_1
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    CALL PUSHCONTROL2B(3)
                    x9 = 0.
                  END IF
                  IF (0. .GT. pmp_1) THEN
                    IF (pmp_1 .GT. lac_1) THEN
                      y14 = lac_1
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y14 = pmp_1
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (0. .GT. lac_1) THEN
                    y14 = lac_1
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y14 = 0.
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (bl(i) .LT. y14) THEN
                    y8 = y14
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    y8 = bl(i)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (x9 .GT. y8) THEN
                    bl(i) = y8
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    bl(i) = x9
                    CALL PUSHCONTROL1B(1)
                  END IF
                  pmp_2 = 2.*dq(i-1)
                  lac_2 = pmp_2 - 1.5*dq(i-2)
                  IF (0. .LT. pmp_2) THEN
                    IF (pmp_2 .LT. lac_2) THEN
                      x10 = lac_2
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      x10 = pmp_2
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (0. .LT. lac_2) THEN
                    x10 = lac_2
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    CALL PUSHCONTROL2B(3)
                    x10 = 0.
                  END IF
                  IF (0. .GT. pmp_2) THEN
                    IF (pmp_2 .GT. lac_2) THEN
                      y15 = lac_2
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y15 = pmp_2
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (0. .GT. lac_2) THEN
                    y15 = lac_2
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y15 = 0.
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (br(i) .LT. y15) THEN
                    y9 = y15
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    y9 = br(i)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (x10 .GT. y9) THEN
                    br(i) = y9
                    CALL PUSHCONTROL3B(1)
                  ELSE
                    br(i) = x10
                    CALL PUSHCONTROL3B(2)
                  END IF
                ELSE
                  CALL PUSHCONTROL3B(0)
                END IF
              END IF
            END DO
            CALL PUSHCONTROL2B(2)
          ELSE
! un-limited: 11
            DO i=is3,ie3
              CALL PUSHREALARRAY_ADM(bl(i))
              bl(i) = al(i) - u(i, j)
              CALL PUSHREALARRAY_ADM(br(i))
              br(i) = al(i+1) - u(i, j)
            END DO
            CALL PUSHCONTROL2B(3)
          END IF
!--------------
! fix the edges
!--------------
!!! TO DO: separate versions for nested and for cubed-sphere
          IF (is .EQ. 1 .AND. (.NOT.nested)) THEN
            CALL PUSHREALARRAY_ADM(br(2))
            br(2) = al(3) - u(2, j)
            CALL PUSHREALARRAY_ADM(xt)
            xt = s15*u(1, j) + s11*u(2, j) - s14*dm(2)
            CALL PUSHREALARRAY_ADM(bl(2))
            bl(2) = xt - u(2, j)
            CALL PUSHREALARRAY_ADM(br(1))
            br(1) = xt - u(1, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! out
              CALL PUSHREALARRAY_ADM(bl(0))
              bl(0) = 0.
! edge
              CALL PUSHREALARRAY_ADM(br(0))
              br(0) = 0.
! edge
              CALL PUSHREALARRAY_ADM(bl(1))
              bl(1) = 0.
! in
              CALL PUSHREALARRAY_ADM(br(1))
              br(1) = 0.
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(bl(0))
              bl(0) = s14*dm(-1) - s11*dq(-1)
              x0l = 0.5*((2.*dx(0, j)+dx(-1, j))*u(0, j)-dx(0, j)*u(-1, &
&               j))/(dx(0, j)+dx(-1, j))
              x0r = 0.5*((2.*dx(1, j)+dx(2, j))*u(1, j)-dx(1, j)*u(2, j)&
&               )/(dx(1, j)+dx(2, j))
              xt = x0l + x0r
              CALL PUSHREALARRAY_ADM(br(0))
              br(0) = xt - u(0, j)
              CALL PUSHREALARRAY_ADM(bl(1))
              bl(1) = xt - u(1, j)
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREALARRAY_ADM(br(2:2), 1)
            CALL PUSHREALARRAY_ADM(bl(2:2), 1)
            CALL PERT_PPM(1, u(2:2, j), bl(2:2), br(2:2), -1)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx .AND. (.NOT.nested)) THEN
            CALL PUSHREALARRAY_ADM(bl(npx-2))
            bl(npx-2) = al(npx-2) - u(npx-2, j)
            CALL PUSHREALARRAY_ADM(xt)
            xt = s15*u(npx-1, j) + s11*u(npx-2, j) + s14*dm(npx-2)
            CALL PUSHREALARRAY_ADM(br(npx-2))
            br(npx-2) = xt - u(npx-2, j)
            CALL PUSHREALARRAY_ADM(bl(npx-1))
            bl(npx-1) = xt - u(npx-1, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! in
              CALL PUSHREALARRAY_ADM(bl(npx-1))
              bl(npx-1) = 0.
! edge
              CALL PUSHREALARRAY_ADM(br(npx-1))
              br(npx-1) = 0.
! edge
              CALL PUSHREALARRAY_ADM(bl(npx))
              bl(npx) = 0.
! out
              CALL PUSHREALARRAY_ADM(br(npx))
              br(npx) = 0.
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(br(npx))
              br(npx) = s11*dq(npx) - s14*dm(npx+1)
              x0l = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*u(npx-1, j)-dx(&
&               npx-1, j)*u(npx-2, j))/(dx(npx-1, j)+dx(npx-2, j))
              x0r = 0.5*((2.*dx(npx, j)+dx(npx+1, j))*u(npx, j)-dx(npx, &
&               j)*u(npx+1, j))/(dx(npx, j)+dx(npx+1, j))
              xt = x0l + x0r
              CALL PUSHREALARRAY_ADM(br(npx-1))
              br(npx-1) = xt - u(npx-1, j)
              CALL PUSHREALARRAY_ADM(bl(npx))
              bl(npx) = xt - u(npx, j)
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREALARRAY_ADM(br(npx-2:npx-2), 1)
            CALL PUSHREALARRAY_ADM(bl(npx-2:npx-2), 1)
            CALL PERT_PPM(1, u(npx-2:npx-2, j), bl(npx-2:npx-2), br(npx-&
&                   2:npx-2), -1)
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
! Other grids:
          DO i=is-1,ie+2
            al(i) = 0.5*(u(i-1, j)+u(i, j)) + r3*(dm(i-1)-dm(i))
          END DO
          DO i=is-1,ie+1
            pmp = -(2.*dq(i))
            lac = pmp + 1.5*dq(i+1)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x11 = lac
                CALL PUSHCONTROL2B(0)
              ELSE
                x11 = pmp
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (0. .LT. lac) THEN
              x11 = lac
              CALL PUSHCONTROL2B(2)
            ELSE
              CALL PUSHCONTROL2B(3)
              x11 = 0.
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y16 = lac
                CALL PUSHCONTROL2B(0)
              ELSE
                y16 = pmp
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (0. .GT. lac) THEN
              y16 = lac
              CALL PUSHCONTROL2B(2)
            ELSE
              y16 = 0.
              CALL PUSHCONTROL2B(3)
            END IF
            IF (al(i) - u(i, j) .LT. y16) THEN
              y10 = y16
              CALL PUSHCONTROL1B(0)
            ELSE
              y10 = al(i) - u(i, j)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x11 .GT. y10) THEN
              CALL PUSHREALARRAY_ADM(bl(i))
              bl(i) = y10
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(bl(i))
              bl(i) = x11
              CALL PUSHCONTROL1B(1)
            END IF
            pmp = 2.*dq(i-1)
            lac = pmp - 1.5*dq(i-2)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x12 = lac
                CALL PUSHCONTROL2B(0)
              ELSE
                x12 = pmp
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (0. .LT. lac) THEN
              x12 = lac
              CALL PUSHCONTROL2B(2)
            ELSE
              CALL PUSHCONTROL2B(3)
              x12 = 0.
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y17 = lac
                CALL PUSHCONTROL2B(0)
              ELSE
                y17 = pmp
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (0. .GT. lac) THEN
              y17 = lac
              CALL PUSHCONTROL2B(2)
            ELSE
              y17 = 0.
              CALL PUSHCONTROL2B(3)
            END IF
            IF (al(i+1) - u(i, j) .LT. y17) THEN
              y11 = y17
              CALL PUSHCONTROL1B(0)
            ELSE
              y11 = al(i+1) - u(i, j)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x12 .GT. y11) THEN
              CALL PUSHREALARRAY_ADM(br(i))
              br(i) = y11
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(br(i))
              br(i) = x12
              CALL PUSHCONTROL1B(1)
            END IF
          END DO
          CALL PUSHCONTROL2B(0)
        END IF
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      dm_ad = 0.0
      dq_ad = 0.0
      al_ad = 0.0
      bl_ad = 0.0
      br_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            cfl = c(i, j)*rdx(i, j)
            temp_ad14 = (cfl+1.)*flux_ad(i, j)
            u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
            cfl_ad = (bl(i)+br(i))*temp_ad14 + (bl(i)+cfl*(bl(i)+br(i)))&
&             *flux_ad(i, j)
            bl_ad(i) = bl_ad(i) + (cfl+1.0)*temp_ad14
            br_ad(i) = br_ad(i) + cfl*temp_ad14
            flux_ad(i, j) = 0.0
            c_ad(i, j) = c_ad(i, j) + rdx(i, j)*cfl_ad
          ELSE
            cfl = c(i, j)*rdx(i-1, j)
            temp = bl(i-1) + br(i-1)
            temp_ad13 = (1.-cfl)*flux_ad(i, j)
            u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
            cfl_ad = -(temp*temp_ad13) - (br(i-1)-cfl*temp)*flux_ad(i, j&
&             )
            br_ad(i-1) = br_ad(i-1) + (1.0-cfl)*temp_ad13
            bl_ad(i-1) = bl_ad(i-1) - cfl*temp_ad13
            flux_ad(i, j) = 0.0
            c_ad(i, j) = c_ad(i, j) + rdx(i-1, j)*cfl_ad
          END IF
        END DO
        CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          DO i=ie+1,is-1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(br(i))
              y11_ad = br_ad(i)
              br_ad(i) = 0.0
              x12_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(br(i))
              x12_ad = br_ad(i)
              br_ad(i) = 0.0
              y11_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y17_ad = y11_ad
            ELSE
              al_ad(i+1) = al_ad(i+1) + y11_ad
              u_ad(i, j) = u_ad(i, j) - y11_ad
              y17_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_ad = y17_ad
                pmp_ad = 0.0
              ELSE
                pmp_ad = y17_ad
                lac_ad = 0.0
              END IF
            ELSE
              IF (branch .EQ. 2) THEN
                lac_ad = y17_ad
              ELSE
                lac_ad = 0.0
              END IF
              pmp_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_ad = lac_ad + x12_ad
              ELSE
                pmp_ad = pmp_ad + x12_ad
              END IF
            ELSE IF (branch .EQ. 2) THEN
              lac_ad = lac_ad + x12_ad
            END IF
            pmp_ad = pmp_ad + lac_ad
            dq_ad(i-2) = dq_ad(i-2) - 1.5*lac_ad
            dq_ad(i-1) = dq_ad(i-1) + 2.*pmp_ad
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(bl(i))
              y10_ad = bl_ad(i)
              bl_ad(i) = 0.0
              x11_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(bl(i))
              x11_ad = bl_ad(i)
              bl_ad(i) = 0.0
              y10_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y16_ad = y10_ad
            ELSE
              al_ad(i) = al_ad(i) + y10_ad
              u_ad(i, j) = u_ad(i, j) - y10_ad
              y16_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_ad = y16_ad
                pmp_ad = 0.0
              ELSE
                pmp_ad = y16_ad
                lac_ad = 0.0
              END IF
            ELSE
              IF (branch .EQ. 2) THEN
                lac_ad = y16_ad
              ELSE
                lac_ad = 0.0
              END IF
              pmp_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_ad = lac_ad + x11_ad
              ELSE
                pmp_ad = pmp_ad + x11_ad
              END IF
            ELSE IF (branch .EQ. 2) THEN
              lac_ad = lac_ad + x11_ad
            END IF
            pmp_ad = pmp_ad + lac_ad
            dq_ad(i+1) = dq_ad(i+1) + 1.5*lac_ad
            dq_ad(i) = dq_ad(i) - 2.*pmp_ad
          END DO
          DO i=ie+2,is-1,-1
            u_ad(i-1, j) = u_ad(i-1, j) + 0.5*al_ad(i)
            u_ad(i, j) = u_ad(i, j) + 0.5*al_ad(i)
            dm_ad(i-1) = dm_ad(i-1) + r3*al_ad(i)
            dm_ad(i) = dm_ad(i) - r3*al_ad(i)
            al_ad(i) = 0.0
          END DO
        ELSE
          IF (branch .NE. 1) THEN
            CALL POPREALARRAY_ADM(bl(npx-2:npx-2), 1)
            CALL POPREALARRAY_ADM(br(npx-2:npx-2), 1)
            CALL PERT_PPM_ADM(1, u(npx-2:npx-2, j), bl(npx-2:npx-2), &
&                       bl_ad(npx-2:npx-2), br(npx-2:npx-2), br_ad(npx-2&
&                       :npx-2), -1)
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(br(npx))
              br_ad(npx) = 0.0
              CALL POPREALARRAY_ADM(bl(npx))
              bl_ad(npx) = 0.0
              CALL POPREALARRAY_ADM(br(npx-1))
              br_ad(npx-1) = 0.0
              CALL POPREALARRAY_ADM(bl(npx-1))
              bl_ad(npx-1) = 0.0
            ELSE
              CALL POPREALARRAY_ADM(bl(npx))
              xt_ad = br_ad(npx-1) + bl_ad(npx)
              u_ad(npx, j) = u_ad(npx, j) - bl_ad(npx)
              bl_ad(npx) = 0.0
              CALL POPREALARRAY_ADM(br(npx-1))
              u_ad(npx-1, j) = u_ad(npx-1, j) - br_ad(npx-1)
              br_ad(npx-1) = 0.0
              x0l_ad = xt_ad
              x0r_ad = xt_ad
              temp_ad11 = 0.5*x0r_ad/(dx(npx, j)+dx(npx+1, j))
              u_ad(npx, j) = u_ad(npx, j) + (dx(npx, j)*2.+dx(npx+1, j))&
&               *temp_ad11
              u_ad(npx+1, j) = u_ad(npx+1, j) - dx(npx, j)*temp_ad11
              temp_ad12 = 0.5*x0l_ad/(dx(npx-1, j)+dx(npx-2, j))
              u_ad(npx-1, j) = u_ad(npx-1, j) + (dx(npx-1, j)*2.+dx(npx-&
&               2, j))*temp_ad12
              u_ad(npx-2, j) = u_ad(npx-2, j) - dx(npx-1, j)*temp_ad12
              CALL POPREALARRAY_ADM(br(npx))
              dq_ad(npx) = dq_ad(npx) + s11*br_ad(npx)
              dm_ad(npx+1) = dm_ad(npx+1) - s14*br_ad(npx)
              br_ad(npx) = 0.0
            END IF
            CALL POPREALARRAY_ADM(bl(npx-1))
            xt_ad = br_ad(npx-2) + bl_ad(npx-1)
            u_ad(npx-1, j) = u_ad(npx-1, j) - bl_ad(npx-1)
            bl_ad(npx-1) = 0.0
            CALL POPREALARRAY_ADM(br(npx-2))
            u_ad(npx-2, j) = u_ad(npx-2, j) - br_ad(npx-2)
            br_ad(npx-2) = 0.0
            CALL POPREALARRAY_ADM(xt)
            u_ad(npx-1, j) = u_ad(npx-1, j) + s15*xt_ad
            u_ad(npx-2, j) = u_ad(npx-2, j) + s11*xt_ad - bl_ad(npx-2)
            dm_ad(npx-2) = dm_ad(npx-2) + s14*xt_ad
            CALL POPREALARRAY_ADM(bl(npx-2))
            al_ad(npx-2) = al_ad(npx-2) + bl_ad(npx-2)
            bl_ad(npx-2) = 0.0
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY_ADM(bl(2:2), 1)
            CALL POPREALARRAY_ADM(br(2:2), 1)
            CALL PERT_PPM_ADM(1, u(2:2, j), bl(2:2), bl_ad(2:2), br(2:2)&
&                       , br_ad(2:2), -1)
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(br(1))
              br_ad(1) = 0.0
              CALL POPREALARRAY_ADM(bl(1))
              bl_ad(1) = 0.0
              CALL POPREALARRAY_ADM(br(0))
              br_ad(0) = 0.0
              CALL POPREALARRAY_ADM(bl(0))
              bl_ad(0) = 0.0
            ELSE
              CALL POPREALARRAY_ADM(bl(1))
              xt_ad = br_ad(0) + bl_ad(1)
              u_ad(1, j) = u_ad(1, j) - bl_ad(1)
              bl_ad(1) = 0.0
              CALL POPREALARRAY_ADM(br(0))
              u_ad(0, j) = u_ad(0, j) - br_ad(0)
              br_ad(0) = 0.0
              x0l_ad = xt_ad
              x0r_ad = xt_ad
              temp_ad9 = 0.5*x0r_ad/(dx(1, j)+dx(2, j))
              u_ad(1, j) = u_ad(1, j) + (dx(1, j)*2.+dx(2, j))*temp_ad9
              u_ad(2, j) = u_ad(2, j) - dx(1, j)*temp_ad9
              temp_ad10 = 0.5*x0l_ad/(dx(0, j)+dx(-1, j))
              u_ad(0, j) = u_ad(0, j) + (dx(0, j)*2.+dx(-1, j))*&
&               temp_ad10
              u_ad(-1, j) = u_ad(-1, j) - dx(0, j)*temp_ad10
              CALL POPREALARRAY_ADM(bl(0))
              dm_ad(-1) = dm_ad(-1) + s14*bl_ad(0)
              dq_ad(-1) = dq_ad(-1) - s11*bl_ad(0)
              bl_ad(0) = 0.0
            END IF
            CALL POPREALARRAY_ADM(br(1))
            xt_ad = bl_ad(2) + br_ad(1)
            u_ad(1, j) = u_ad(1, j) - br_ad(1)
            br_ad(1) = 0.0
            CALL POPREALARRAY_ADM(bl(2))
            u_ad(2, j) = u_ad(2, j) - bl_ad(2)
            bl_ad(2) = 0.0
            CALL POPREALARRAY_ADM(xt)
            u_ad(1, j) = u_ad(1, j) + s15*xt_ad
            u_ad(2, j) = u_ad(2, j) + s11*xt_ad - br_ad(2)
            dm_ad(2) = dm_ad(2) - s14*xt_ad
            CALL POPREALARRAY_ADM(br(2))
            al_ad(3) = al_ad(3) + br_ad(2)
            br_ad(2) = 0.0
          END IF
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              DO i=ie3,is3,-1
                CALL POPREALARRAY_ADM(br(i))
                min5_ad = SIGN(1.d0, min5*xt)*br_ad(i)
                br_ad(i) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(min5)
                  y5_ad = min5_ad
                  x6_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(min5)
                  x6_ad = min5_ad
                  y5_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  al_ad(i+1) = al_ad(i+1) + y5_ad
                  u_ad(i, j) = u_ad(i, j) - y5_ad
                ELSE
                  u_ad(i, j) = u_ad(i, j) + y5_ad
                  al_ad(i+1) = al_ad(i+1) - y5_ad
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  xt_ad = x6_ad
                ELSE
                  xt_ad = -x6_ad
                END IF
                CALL POPREALARRAY_ADM(bl(i))
                min4_ad = -(SIGN(1.d0, min4*xt)*bl_ad(i))
                bl_ad(i) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(min4)
                  y4_ad = min4_ad
                  x5_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(min4)
                  x5_ad = min4_ad
                  y4_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  al_ad(i) = al_ad(i) + y4_ad
                  u_ad(i, j) = u_ad(i, j) - y4_ad
                ELSE
                  u_ad(i, j) = u_ad(i, j) + y4_ad
                  al_ad(i) = al_ad(i) - y4_ad
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  xt_ad = xt_ad + x5_ad
                ELSE
                  xt_ad = xt_ad - x5_ad
                END IF
                CALL POPREALARRAY_ADM(xt)
                dm_ad(i) = dm_ad(i) + 2.*xt_ad
              END DO
            ELSE
              DO i=ie3,is3,-1
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(br(i))
                  y7_ad = br_ad(i)
                  br_ad(i) = 0.0
                  x8_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(br(i))
                  x8_ad = br_ad(i)
                  br_ad(i) = 0.0
                  y7_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y13_ad = y7_ad
                ELSE
                  al_ad(i+1) = al_ad(i+1) + y7_ad
                  u_ad(i, j) = u_ad(i, j) - y7_ad
                  y13_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = y13_ad
                    pmp_2_ad = 0.0
                  ELSE
                    pmp_2_ad = y13_ad
                    lac_2_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_2_ad = y13_ad
                  ELSE
                    lac_2_ad = 0.0
                  END IF
                  pmp_2_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = lac_2_ad + x8_ad
                  ELSE
                    pmp_2_ad = pmp_2_ad + x8_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_2_ad = lac_2_ad + x8_ad
                END IF
                pmp_2_ad = pmp_2_ad + lac_2_ad
                dq_ad(i-2) = dq_ad(i-2) - 1.5*lac_2_ad
                dq_ad(i-1) = dq_ad(i-1) + 2.*pmp_2_ad
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(bl(i))
                  y6_ad = bl_ad(i)
                  bl_ad(i) = 0.0
                  x7_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(bl(i))
                  x7_ad = bl_ad(i)
                  bl_ad(i) = 0.0
                  y6_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y12_ad = y6_ad
                ELSE
                  al_ad(i) = al_ad(i) + y6_ad
                  u_ad(i, j) = u_ad(i, j) - y6_ad
                  y12_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = y12_ad
                    pmp_1_ad = 0.0
                  ELSE
                    pmp_1_ad = y12_ad
                    lac_1_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_1_ad = y12_ad
                  ELSE
                    lac_1_ad = 0.0
                  END IF
                  pmp_1_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = lac_1_ad + x7_ad
                  ELSE
                    pmp_1_ad = pmp_1_ad + x7_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_1_ad = lac_1_ad + x7_ad
                END IF
                pmp_1_ad = pmp_1_ad + lac_1_ad
                dq_ad(i+1) = dq_ad(i+1) + 1.5*lac_1_ad
                dq_ad(i) = dq_ad(i) - 2.*pmp_1_ad
              END DO
            END IF
          ELSE IF (branch .EQ. 2) THEN
            DO i=ie3,is3,-1
              CALL POPCONTROL3B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  GOTO 110
                ELSE
                  y9_ad = br_ad(i)
                  br_ad(i) = 0.0
                  x10_ad = 0.0
                END IF
              ELSE IF (branch .EQ. 2) THEN
                x10_ad = br_ad(i)
                br_ad(i) = 0.0
                y9_ad = 0.0
              ELSE
                IF (branch .NE. 3) THEN
                  br_ad(i) = 0.0
                  bl_ad(i) = 0.0
                END IF
                GOTO 110
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                y15_ad = y9_ad
              ELSE
                br_ad(i) = br_ad(i) + y9_ad
                y15_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_2_ad = y15_ad
                  pmp_2_ad = 0.0
                ELSE
                  pmp_2_ad = y15_ad
                  lac_2_ad = 0.0
                END IF
              ELSE
                IF (branch .EQ. 2) THEN
                  lac_2_ad = y15_ad
                ELSE
                  lac_2_ad = 0.0
                END IF
                pmp_2_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_2_ad = lac_2_ad + x10_ad
                ELSE
                  pmp_2_ad = pmp_2_ad + x10_ad
                END IF
              ELSE IF (branch .EQ. 2) THEN
                lac_2_ad = lac_2_ad + x10_ad
              END IF
              pmp_2_ad = pmp_2_ad + lac_2_ad
              dq_ad(i-2) = dq_ad(i-2) - 1.5*lac_2_ad
              dq_ad(i-1) = dq_ad(i-1) + 2.*pmp_2_ad
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                y8_ad = bl_ad(i)
                bl_ad(i) = 0.0
                x9_ad = 0.0
              ELSE
                x9_ad = bl_ad(i)
                bl_ad(i) = 0.0
                y8_ad = 0.0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                y14_ad = y8_ad
              ELSE
                bl_ad(i) = bl_ad(i) + y8_ad
                y14_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_1_ad = y14_ad
                  pmp_1_ad = 0.0
                ELSE
                  pmp_1_ad = y14_ad
                  lac_1_ad = 0.0
                END IF
              ELSE
                IF (branch .EQ. 2) THEN
                  lac_1_ad = y14_ad
                ELSE
                  lac_1_ad = 0.0
                END IF
                pmp_1_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_1_ad = lac_1_ad + x9_ad
                ELSE
                  pmp_1_ad = pmp_1_ad + x9_ad
                END IF
              ELSE IF (branch .EQ. 2) THEN
                lac_1_ad = lac_1_ad + x9_ad
              END IF
              pmp_1_ad = pmp_1_ad + lac_1_ad
              dq_ad(i+1) = dq_ad(i+1) + 1.5*lac_1_ad
              dq_ad(i) = dq_ad(i) - 2.*pmp_1_ad
 110          CALL POPREALARRAY_ADM(br(i))
              al_ad(i+1) = al_ad(i+1) + br_ad(i)
              u_ad(i, j) = u_ad(i, j) - bl_ad(i) - br_ad(i)
              br_ad(i) = 0.0
              CALL POPREALARRAY_ADM(bl(i))
              al_ad(i) = al_ad(i) + bl_ad(i)
              bl_ad(i) = 0.0
            END DO
          ELSE
            DO i=ie3,is3,-1
              CALL POPREALARRAY_ADM(br(i))
              al_ad(i+1) = al_ad(i+1) + br_ad(i)
              u_ad(i, j) = u_ad(i, j) - bl_ad(i) - br_ad(i)
              br_ad(i) = 0.0
              CALL POPREALARRAY_ADM(bl(i))
              al_ad(i) = al_ad(i) + bl_ad(i)
              bl_ad(i) = 0.0
            END DO
          END IF
          DO i=ie3+1,is3,-1
            u_ad(i-1, j) = u_ad(i-1, j) + 0.5*al_ad(i)
            u_ad(i, j) = u_ad(i, j) + 0.5*al_ad(i)
            dm_ad(i-1) = dm_ad(i-1) + r3*al_ad(i)
            dm_ad(i) = dm_ad(i) - r3*al_ad(i)
            al_ad(i) = 0.0
          END DO
        END IF
        DO i=ie+2,is-3,-1
          u_ad(i+1, j) = u_ad(i+1, j) + dq_ad(i)
          u_ad(i, j) = u_ad(i, j) - dq_ad(i)
          dq_ad(i) = 0.0
        END DO
        DO i=ie+2,is-2,-1
          xt = 0.25*(u(i+1, j)-u(i-1, j))
          min3_ad = SIGN(1.d0, min3*xt)*dm_ad(i)
          dm_ad(i) = 0.0
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(min3)
              z1_ad = min3_ad
              y3_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min3)
              y3_ad = min3_ad
              z1_ad = 0.0
            END IF
            x4_ad = 0.0
          ELSE
            IF (branch .EQ. 2) THEN
              CALL POPREALARRAY_ADM(min3)
              z1_ad = min3_ad
              x4_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min3)
              x4_ad = min3_ad
              z1_ad = 0.0
            END IF
            y3_ad = 0.0
          END IF
          u_ad(i, j) = u_ad(i, j) + z1_ad
          min6_ad = -z1_ad
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              u_ad(i+1, j) = u_ad(i+1, j) + min6_ad
            ELSE
              u_ad(i, j) = u_ad(i, j) + min6_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            u_ad(i+1, j) = u_ad(i+1, j) + min6_ad
          ELSE
            u_ad(i-1, j) = u_ad(i-1, j) + min6_ad
          END IF
          max1_ad = y3_ad
          u_ad(i, j) = u_ad(i, j) - y3_ad
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              u_ad(i+1, j) = u_ad(i+1, j) + max1_ad
            ELSE
              u_ad(i, j) = u_ad(i, j) + max1_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            u_ad(i+1, j) = u_ad(i+1, j) + max1_ad
          ELSE
            u_ad(i-1, j) = u_ad(i-1, j) + max1_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            xt_ad = x4_ad
          ELSE
            xt_ad = -x4_ad
          END IF
          CALL POPREALARRAY_ADM(xt)
          u_ad(i+1, j) = u_ad(i+1, j) + 0.25*xt_ad
          u_ad(i-1, j) = u_ad(i-1, j) - 0.25*xt_ad
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
  END SUBROUTINE XTP_U_ADM
  SUBROUTINE XTP_U(is, ie, js, je, isd, ied, jsd, jed, c, u, v, flux, &
&   iord, dx, rdx, npx, npy, grid_type, nested)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL, INTENT(OUT) :: flux(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: rdx(isd:ied, jsd:jed+1)
    INTEGER, INTENT(IN) :: iord, npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
! Local
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: al(is-1:ie+2), dm(is-2:ie+2)
    REAL :: dq(is-3:ie+2)
    REAL :: dl, dr, xt, pmp, lac, cfl
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: x0, x1, x0l, x0r
    INTEGER :: i, j
    INTEGER :: is3, ie3
    INTEGER :: is2, ie2
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min2
    REAL :: abs0
    REAL :: min3
    REAL :: min4
    REAL :: min5
    REAL :: abs1
    REAL :: abs2
    REAL :: abs3
    REAL :: abs4
    REAL :: max1
    REAL :: min6
    REAL :: abs5
    REAL :: abs6
    REAL :: x12
    REAL :: x11
    REAL :: x10
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: y17
    REAL :: y16
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
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
    IF (nested .OR. grid_type .GT. 3) THEN
      is3 = is - 1
      ie3 = ie + 1
    ELSE
      IF (3 .LT. is - 1) THEN
        is3 = is - 1
      ELSE
        is3 = 3
      END IF
      IF (npx - 3 .GT. ie + 1) THEN
        ie3 = ie + 1
      ELSE
        ie3 = npx - 3
      END IF
    END IF
    IF (iord .EQ. 1) THEN
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux(i, j) = u(i-1, j)
          ELSE
            flux(i, j) = u(i, j)
          END IF
        END DO
      END DO
    ELSE IF (iord .LT. 8) THEN
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
      DO j=js,je+1
        DO i=is3,ie3+1
          al(i) = p1*(u(i-1, j)+u(i, j)) + p2*(u(i-2, j)+u(i+1, j))
        END DO
        DO i=is3,ie3
          bl(i) = al(i) - u(i, j)
          br(i) = al(i+1) - u(i, j)
        END DO
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            xt = c3*u(1, j) + c2*u(2, j) + c1*u(3, j)
            br(1) = xt - u(1, j)
            bl(2) = xt - u(2, j)
            br(2) = al(3) - u(2, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! out
              bl(0) = 0.
! edge
              br(0) = 0.
! edge
              bl(1) = 0.
! in
              br(1) = 0.
            ELSE
              bl(0) = c1*u(-2, j) + c2*u(-1, j) + c3*u(0, j) - u(0, j)
              xt = 0.5*(((2.*dx(0, j)+dx(-1, j))*u(0, j)-dx(0, j)*u(-1, &
&               j))/(dx(0, j)+dx(-1, j))+((2.*dx(1, j)+dx(2, j))*u(1, j)&
&               -dx(1, j)*u(2, j))/(dx(1, j)+dx(2, j)))
              br(0) = xt - u(0, j)
              bl(1) = xt - u(1, j)
            END IF
          END IF
!       call pert_ppm(1, u(2,j), bl(2), br(2), -1)
          IF (ie + 1 .EQ. npx) THEN
            bl(npx-2) = al(npx-2) - u(npx-2, j)
            xt = c1*u(npx-3, j) + c2*u(npx-2, j) + c3*u(npx-1, j)
            br(npx-2) = xt - u(npx-2, j)
            bl(npx-1) = xt - u(npx-1, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! in
              bl(npx-1) = 0.
! edge
              br(npx-1) = 0.
! edge
              bl(npx) = 0.
! out
              br(npx) = 0.
            ELSE
              xt = 0.5*(((2.*dx(npx-1, j)+dx(npx-2, j))*u(npx-1, j)-dx(&
&               npx-1, j)*u(npx-2, j))/(dx(npx-1, j)+dx(npx-2, j))+((2.*&
&               dx(npx, j)+dx(npx+1, j))*u(npx, j)-dx(npx, j)*u(npx+1, j&
&               ))/(dx(npx, j)+dx(npx+1, j)))
              br(npx-1) = xt - u(npx-1, j)
              bl(npx) = xt - u(npx, j)
              br(npx) = c3*u(npx, j) + c2*u(npx+1, j) + c1*u(npx+2, j) -&
&               u(npx, j)
            END IF
          END IF
        END IF
!       call pert_ppm(1, u(npx-2,j), bl(npx-2), br(npx-2), -1)
        DO i=is-1,ie+1
          b0(i) = bl(i) + br(i)
        END DO
        IF (iord .EQ. 2) THEN
! Perfectly linear
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdx(i-1, j)
              flux(i, j) = u(i-1, j) + (1.-cfl)*(br(i-1)-cfl*b0(i-1))
            ELSE
              cfl = c(i, j)*rdx(i, j)
              flux(i, j) = u(i, j) + (1.+cfl)*(bl(i)+cfl*b0(i))
            END IF
          END DO
        ELSE IF (iord .EQ. 3) THEN
          DO i=is-1,ie+1
            IF (b0(i) .GE. 0.) THEN
              x0 = b0(i)
            ELSE
              x0 = -b0(i)
            END IF
            IF (bl(i) - br(i) .GE. 0.) THEN
              x1 = bl(i) - br(i)
            ELSE
              x1 = -(bl(i)-br(i))
            END IF
            smt5(i) = x0 .LT. x1
            smt6(i) = 3.*x0 .LT. x1
          END DO
          DO i=is,ie+1
            fx0(i) = 0.
          END DO
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdx(i-1, j)
              IF (smt6(i-1) .OR. smt5(i)) THEN
                fx0(i) = br(i-1) - cfl*b0(i-1)
              ELSE IF (smt5(i-1)) THEN
                IF (bl(i-1) .GE. 0.) THEN
                  x2 = bl(i-1)
                ELSE
                  x2 = -bl(i-1)
                END IF
                IF (br(i-1) .GE. 0.) THEN
                  y1 = br(i-1)
                ELSE
                  y1 = -br(i-1)
                END IF
                IF (x2 .GT. y1) THEN
                  min1 = y1
                ELSE
                  min1 = x2
                END IF
                fx0(i) = SIGN(min1, br(i-1))
              END IF
              flux(i, j) = u(i-1, j) + (1.-cfl)*fx0(i)
            ELSE
              cfl = c(i, j)*rdx(i, j)
              IF (smt6(i) .OR. smt5(i-1)) THEN
                fx0(i) = bl(i) + cfl*b0(i)
              ELSE IF (smt5(i)) THEN
                IF (bl(i) .GE. 0.) THEN
                  x3 = bl(i)
                ELSE
                  x3 = -bl(i)
                END IF
                IF (br(i) .GE. 0.) THEN
                  y2 = br(i)
                ELSE
                  y2 = -br(i)
                END IF
                IF (x3 .GT. y2) THEN
                  min2 = y2
                ELSE
                  min2 = x3
                END IF
                fx0(i) = SIGN(min2, bl(i))
              END IF
              flux(i, j) = u(i, j) + (1.+cfl)*fx0(i)
            END IF
          END DO
        ELSE IF (iord .EQ. 4) THEN
! more damp than ord5 but less damp than ord6
          DO i=is-1,ie+1
            IF (b0(i) .GE. 0.) THEN
              x0 = b0(i)
            ELSE
              x0 = -b0(i)
            END IF
            IF (bl(i) - br(i) .GE. 0.) THEN
              x1 = bl(i) - br(i)
            ELSE
              x1 = -(bl(i)-br(i))
            END IF
            smt5(i) = x0 .LT. x1
! if smt6 =.T. --> smt5=.T.
            smt6(i) = 3.*x0 .LT. x1
          END DO
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              IF (smt6(i-1) .OR. smt5(i)) THEN
                cfl = c(i, j)*rdx(i-1, j)
                flux(i, j) = u(i-1, j) + (1.-cfl)*(br(i-1)-cfl*b0(i-1))
              ELSE
! 1st order ONLY_IF smt6(i-1)=.F.  .AND. smt5(i)=.F.
                flux(i, j) = u(i-1, j)
              END IF
            ELSE IF (smt6(i) .OR. smt5(i-1)) THEN
              cfl = c(i, j)*rdx(i, j)
              flux(i, j) = u(i, j) + (1.+cfl)*(bl(i)+cfl*b0(i))
            ELSE
              flux(i, j) = u(i, j)
            END IF
          END DO
        ELSE
!  iord=5,6,7
          IF (iord .EQ. 5) THEN
            DO i=is-1,ie+1
              smt5(i) = bl(i)*br(i) .LT. 0.
            END DO
          ELSE
            DO i=is-1,ie+1
              IF (3.*b0(i) .GE. 0.) THEN
                abs0 = 3.*b0(i)
              ELSE
                abs0 = -(3.*b0(i))
              END IF
              IF (bl(i) - br(i) .GE. 0.) THEN
                abs4 = bl(i) - br(i)
              ELSE
                abs4 = -(bl(i)-br(i))
              END IF
              smt5(i) = abs0 .LT. abs4
            END DO
          END IF
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdx(i-1, j)
              fx0(i) = (1.-cfl)*(br(i-1)-cfl*b0(i-1))
              flux(i, j) = u(i-1, j)
            ELSE
              cfl = c(i, j)*rdx(i, j)
              fx0(i) = (1.+cfl)*(bl(i)+cfl*b0(i))
              flux(i, j) = u(i, j)
            END IF
            IF (smt5(i-1) .OR. smt5(i)) flux(i, j) = flux(i, j) + fx0(i)
          END DO
        END IF
      END DO
    ELSE
! iord = 8, 9, 10, 11
      DO j=js,je+1
        DO i=is-2,ie+2
          xt = 0.25*(u(i+1, j)-u(i-1, j))
          IF (xt .GE. 0.) THEN
            x4 = xt
          ELSE
            x4 = -xt
          END IF
          IF (u(i-1, j) .LT. u(i, j)) THEN
            IF (u(i, j) .LT. u(i+1, j)) THEN
              max1 = u(i+1, j)
            ELSE
              max1 = u(i, j)
            END IF
          ELSE IF (u(i-1, j) .LT. u(i+1, j)) THEN
            max1 = u(i+1, j)
          ELSE
            max1 = u(i-1, j)
          END IF
          y3 = max1 - u(i, j)
          IF (u(i-1, j) .GT. u(i, j)) THEN
            IF (u(i, j) .GT. u(i+1, j)) THEN
              min6 = u(i+1, j)
            ELSE
              min6 = u(i, j)
            END IF
          ELSE IF (u(i-1, j) .GT. u(i+1, j)) THEN
            min6 = u(i+1, j)
          ELSE
            min6 = u(i-1, j)
          END IF
          z1 = u(i, j) - min6
          IF (x4 .GT. y3) THEN
            IF (y3 .GT. z1) THEN
              min3 = z1
            ELSE
              min3 = y3
            END IF
          ELSE IF (x4 .GT. z1) THEN
            min3 = z1
          ELSE
            min3 = x4
          END IF
          dm(i) = SIGN(min3, xt)
        END DO
        DO i=is-3,ie+2
          dq(i) = u(i+1, j) - u(i, j)
        END DO
        IF (grid_type .LT. 3) THEN
          DO i=is3,ie3+1
            al(i) = 0.5*(u(i-1, j)+u(i, j)) + r3*(dm(i-1)-dm(i))
          END DO
! Perturbation form:
          IF (iord .EQ. 8) THEN
            DO i=is3,ie3
              xt = 2.*dm(i)
              IF (xt .GE. 0.) THEN
                x5 = xt
              ELSE
                x5 = -xt
              END IF
              IF (al(i) - u(i, j) .GE. 0.) THEN
                y4 = al(i) - u(i, j)
              ELSE
                y4 = -(al(i)-u(i, j))
              END IF
              IF (x5 .GT. y4) THEN
                min4 = y4
              ELSE
                min4 = x5
              END IF
              bl(i) = -SIGN(min4, xt)
              IF (xt .GE. 0.) THEN
                x6 = xt
              ELSE
                x6 = -xt
              END IF
              IF (al(i+1) - u(i, j) .GE. 0.) THEN
                y5 = al(i+1) - u(i, j)
              ELSE
                y5 = -(al(i+1)-u(i, j))
              END IF
              IF (x6 .GT. y5) THEN
                min5 = y5
              ELSE
                min5 = x6
              END IF
              br(i) = SIGN(min5, xt)
            END DO
          ELSE IF (iord .EQ. 9) THEN
            DO i=is3,ie3
              pmp_1 = -(2.*dq(i))
              lac_1 = pmp_1 + 1.5*dq(i+1)
              IF (0. .LT. pmp_1) THEN
                IF (pmp_1 .LT. lac_1) THEN
                  x7 = lac_1
                ELSE
                  x7 = pmp_1
                END IF
              ELSE IF (0. .LT. lac_1) THEN
                x7 = lac_1
              ELSE
                x7 = 0.
              END IF
              IF (0. .GT. pmp_1) THEN
                IF (pmp_1 .GT. lac_1) THEN
                  y12 = lac_1
                ELSE
                  y12 = pmp_1
                END IF
              ELSE IF (0. .GT. lac_1) THEN
                y12 = lac_1
              ELSE
                y12 = 0.
              END IF
              IF (al(i) - u(i, j) .LT. y12) THEN
                y6 = y12
              ELSE
                y6 = al(i) - u(i, j)
              END IF
              IF (x7 .GT. y6) THEN
                bl(i) = y6
              ELSE
                bl(i) = x7
              END IF
              pmp_2 = 2.*dq(i-1)
              lac_2 = pmp_2 - 1.5*dq(i-2)
              IF (0. .LT. pmp_2) THEN
                IF (pmp_2 .LT. lac_2) THEN
                  x8 = lac_2
                ELSE
                  x8 = pmp_2
                END IF
              ELSE IF (0. .LT. lac_2) THEN
                x8 = lac_2
              ELSE
                x8 = 0.
              END IF
              IF (0. .GT. pmp_2) THEN
                IF (pmp_2 .GT. lac_2) THEN
                  y13 = lac_2
                ELSE
                  y13 = pmp_2
                END IF
              ELSE IF (0. .GT. lac_2) THEN
                y13 = lac_2
              ELSE
                y13 = 0.
              END IF
              IF (al(i+1) - u(i, j) .LT. y13) THEN
                y7 = y13
              ELSE
                y7 = al(i+1) - u(i, j)
              END IF
              IF (x8 .GT. y7) THEN
                br(i) = y7
              ELSE
                br(i) = x8
              END IF
            END DO
          ELSE IF (iord .EQ. 10) THEN
            DO i=is3,ie3
              bl(i) = al(i) - u(i, j)
              br(i) = al(i+1) - u(i, j)
              IF (dm(i) .GE. 0.) THEN
                abs1 = dm(i)
              ELSE
                abs1 = -dm(i)
              END IF
!             if ( abs(dm(i-1))+abs(dm(i))+abs(dm(i+1)) < near_zero ) then
              IF (abs1 .LT. near_zero) THEN
                IF (dm(i-1) .GE. 0.) THEN
                  abs2 = dm(i-1)
                ELSE
                  abs2 = -dm(i-1)
                END IF
                IF (dm(i+1) .GE. 0.) THEN
                  abs5 = dm(i+1)
                ELSE
                  abs5 = -dm(i+1)
                END IF
                IF (abs2 + abs5 .LT. near_zero) THEN
! 2-delta-x structure detected within 3 cells
                  bl(i) = 0.
                  br(i) = 0.
                END IF
              ELSE
                IF (3.*(bl(i)+br(i)) .GE. 0.) THEN
                  abs3 = 3.*(bl(i)+br(i))
                ELSE
                  abs3 = -(3.*(bl(i)+br(i)))
                END IF
                IF (bl(i) - br(i) .GE. 0.) THEN
                  abs6 = bl(i) - br(i)
                ELSE
                  abs6 = -(bl(i)-br(i))
                END IF
                IF (abs3 .GT. abs6) THEN
                  pmp_1 = -(2.*dq(i))
                  lac_1 = pmp_1 + 1.5*dq(i+1)
                  IF (0. .LT. pmp_1) THEN
                    IF (pmp_1 .LT. lac_1) THEN
                      x9 = lac_1
                    ELSE
                      x9 = pmp_1
                    END IF
                  ELSE IF (0. .LT. lac_1) THEN
                    x9 = lac_1
                  ELSE
                    x9 = 0.
                  END IF
                  IF (0. .GT. pmp_1) THEN
                    IF (pmp_1 .GT. lac_1) THEN
                      y14 = lac_1
                    ELSE
                      y14 = pmp_1
                    END IF
                  ELSE IF (0. .GT. lac_1) THEN
                    y14 = lac_1
                  ELSE
                    y14 = 0.
                  END IF
                  IF (bl(i) .LT. y14) THEN
                    y8 = y14
                  ELSE
                    y8 = bl(i)
                  END IF
                  IF (x9 .GT. y8) THEN
                    bl(i) = y8
                  ELSE
                    bl(i) = x9
                  END IF
                  pmp_2 = 2.*dq(i-1)
                  lac_2 = pmp_2 - 1.5*dq(i-2)
                  IF (0. .LT. pmp_2) THEN
                    IF (pmp_2 .LT. lac_2) THEN
                      x10 = lac_2
                    ELSE
                      x10 = pmp_2
                    END IF
                  ELSE IF (0. .LT. lac_2) THEN
                    x10 = lac_2
                  ELSE
                    x10 = 0.
                  END IF
                  IF (0. .GT. pmp_2) THEN
                    IF (pmp_2 .GT. lac_2) THEN
                      y15 = lac_2
                    ELSE
                      y15 = pmp_2
                    END IF
                  ELSE IF (0. .GT. lac_2) THEN
                    y15 = lac_2
                  ELSE
                    y15 = 0.
                  END IF
                  IF (br(i) .LT. y15) THEN
                    y9 = y15
                  ELSE
                    y9 = br(i)
                  END IF
                  IF (x10 .GT. y9) THEN
                    br(i) = y9
                  ELSE
                    br(i) = x10
                  END IF
                END IF
              END IF
            END DO
          ELSE
! un-limited: 11
            DO i=is3,ie3
              bl(i) = al(i) - u(i, j)
              br(i) = al(i+1) - u(i, j)
            END DO
          END IF
!--------------
! fix the edges
!--------------
!!! TO DO: separate versions for nested and for cubed-sphere
          IF (is .EQ. 1 .AND. (.NOT.nested)) THEN
            br(2) = al(3) - u(2, j)
            xt = s15*u(1, j) + s11*u(2, j) - s14*dm(2)
            bl(2) = xt - u(2, j)
            br(1) = xt - u(1, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! out
              bl(0) = 0.
! edge
              br(0) = 0.
! edge
              bl(1) = 0.
! in
              br(1) = 0.
            ELSE
              bl(0) = s14*dm(-1) - s11*dq(-1)
              x0l = 0.5*((2.*dx(0, j)+dx(-1, j))*u(0, j)-dx(0, j)*u(-1, &
&               j))/(dx(0, j)+dx(-1, j))
              x0r = 0.5*((2.*dx(1, j)+dx(2, j))*u(1, j)-dx(1, j)*u(2, j)&
&               )/(dx(1, j)+dx(2, j))
              xt = x0l + x0r
              br(0) = xt - u(0, j)
              bl(1) = xt - u(1, j)
            END IF
            CALL PERT_PPM(1, u(2:2, j), bl(2:2), br(2:2), -1)
          END IF
          IF (ie + 1 .EQ. npx .AND. (.NOT.nested)) THEN
            bl(npx-2) = al(npx-2) - u(npx-2, j)
            xt = s15*u(npx-1, j) + s11*u(npx-2, j) + s14*dm(npx-2)
            br(npx-2) = xt - u(npx-2, j)
            bl(npx-1) = xt - u(npx-1, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! in
              bl(npx-1) = 0.
! edge
              br(npx-1) = 0.
! edge
              bl(npx) = 0.
! out
              br(npx) = 0.
            ELSE
              br(npx) = s11*dq(npx) - s14*dm(npx+1)
              x0l = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*u(npx-1, j)-dx(&
&               npx-1, j)*u(npx-2, j))/(dx(npx-1, j)+dx(npx-2, j))
              x0r = 0.5*((2.*dx(npx, j)+dx(npx+1, j))*u(npx, j)-dx(npx, &
&               j)*u(npx+1, j))/(dx(npx, j)+dx(npx+1, j))
              xt = x0l + x0r
              br(npx-1) = xt - u(npx-1, j)
              bl(npx) = xt - u(npx, j)
            END IF
            CALL PERT_PPM(1, u(npx-2:npx-2, j), bl(npx-2:npx-2), br(npx-&
&                   2:npx-2), -1)
          END IF
        ELSE
! Other grids:
          DO i=is-1,ie+2
            al(i) = 0.5*(u(i-1, j)+u(i, j)) + r3*(dm(i-1)-dm(i))
          END DO
          DO i=is-1,ie+1
            pmp = -(2.*dq(i))
            lac = pmp + 1.5*dq(i+1)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x11 = lac
              ELSE
                x11 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x11 = lac
            ELSE
              x11 = 0.
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y16 = lac
              ELSE
                y16 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y16 = lac
            ELSE
              y16 = 0.
            END IF
            IF (al(i) - u(i, j) .LT. y16) THEN
              y10 = y16
            ELSE
              y10 = al(i) - u(i, j)
            END IF
            IF (x11 .GT. y10) THEN
              bl(i) = y10
            ELSE
              bl(i) = x11
            END IF
            pmp = 2.*dq(i-1)
            lac = pmp - 1.5*dq(i-2)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x12 = lac
              ELSE
                x12 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x12 = lac
            ELSE
              x12 = 0.
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y17 = lac
              ELSE
                y17 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y17 = lac
            ELSE
              y17 = 0.
            END IF
            IF (al(i+1) - u(i, j) .LT. y17) THEN
              y11 = y17
            ELSE
              y11 = al(i+1) - u(i, j)
            END IF
            IF (x12 .GT. y11) THEN
              br(i) = y11
            ELSE
              br(i) = x12
            END IF
          END DO
        END IF
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            cfl = c(i, j)*rdx(i-1, j)
            flux(i, j) = u(i-1, j) + (1.-cfl)*(br(i-1)-cfl*(bl(i-1)+br(i&
&             -1)))
          ELSE
            cfl = c(i, j)*rdx(i, j)
            flux(i, j) = u(i, j) + (1.+cfl)*(bl(i)+cfl*(bl(i)+br(i)))
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE XTP_U
!  Differentiation of ytp_v in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_core
!_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn_cor
!e_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_mod.g
!eopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynamics_m
!od.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mod.
!c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.ma
!p_scalar_fb fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scalar_pr
!ofile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz f
!v_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_setup
! fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.compute
!_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_c nh
!_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver nh_u
!tils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw_cor
!e_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core_mod
!.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.compu
!te_divergence_damping_fb sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corners 
!tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circle_di
!st sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: flux v c
!   with respect to varying inputs: v c
  SUBROUTINE YTP_V_ADM(is, ie, js, je, isd, ied, jsd, jed, c, c_ad, u, v&
&   , v_ad, flux, flux_ad, jord, dy, rdy, npx, npy, grid_type, nested)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL :: v_ad(isd:ied+1, jsd:jed)
!  Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL :: c_ad(is:ie+1, js:je+1)
    REAL :: flux(is:ie+1, js:je+1)
    REAL :: flux_ad(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rdy(isd:ied+1, jsd:jed)
    INTEGER, INTENT(IN) :: npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
! Local:
    LOGICAL, DIMENSION(is:ie+1, js-1:je+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: fx0_ad(is:ie+1)
    REAL :: dm(is:ie+1, js-2:je+2)
    REAL :: dm_ad(is:ie+1, js-2:je+2)
    REAL :: al(is:ie+1, js-1:je+2)
    REAL :: al_ad(is:ie+1, js-1:je+2)
    REAL, DIMENSION(is:ie+1, js-1:je+1) :: bl, br, b0
    REAL, DIMENSION(is:ie+1, js-1:je+1) :: bl_ad, br_ad, b0_ad
    REAL :: dq(is:ie+1, js-3:je+2)
    REAL :: dq_ad(is:ie+1, js-3:je+2)
    REAL :: xt, dl, dr, pmp, lac, cfl
    REAL :: xt_ad, pmp_ad, lac_ad, cfl_ad
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: pmp_1_ad, lac_1_ad, pmp_2_ad, lac_2_ad
    REAL :: x0, x1, x0r, x0l
    REAL :: x0r_ad, x0l_ad
    INTEGER :: i, j, is1, ie1, js3, je3
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min1_ad
    REAL :: min2
    REAL :: min2_ad
    REAL :: abs0
    REAL :: min3
    REAL :: min3_ad
    REAL :: min4
    REAL :: min4_ad
    REAL :: min5
    REAL :: min5_ad
    REAL :: abs1
    REAL :: abs2
    REAL :: abs3
    REAL :: abs4
    REAL :: max1
    REAL :: max1_ad
    REAL :: min6
    REAL :: min6_ad
    REAL :: abs5
    REAL :: abs6
    INTEGER :: arg1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: x2_ad
    REAL :: y1_ad
    REAL :: x3_ad
    REAL :: y2_ad
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: x4_ad
    REAL :: y3_ad
    REAL :: z1_ad
    REAL :: x5_ad
    REAL :: y4_ad
    REAL :: x6_ad
    REAL :: y5_ad
    REAL :: x7_ad
    REAL :: y12_ad
    REAL :: y6_ad
    REAL :: x8_ad
    REAL :: y13_ad
    REAL :: y7_ad
    REAL :: x9_ad
    REAL :: y14_ad
    REAL :: y8_ad
    REAL :: x10_ad
    REAL :: y15_ad
    REAL :: y9_ad
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: x11_ad
    REAL :: y16_ad
    REAL :: y10_ad
    REAL :: x12_ad
    REAL :: y17_ad
    REAL :: y11_ad
    REAL :: temp
    REAL :: temp_ad13
    REAL :: temp0
    REAL :: temp_ad14
    INTEGER :: branch
    REAL :: x12
    REAL :: x11
    REAL :: x10
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: y17
    REAL :: y16
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
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
    IF (nested .OR. grid_type .GT. 3) THEN
      CALL PUSHCONTROL1B(0)
      js3 = js - 1
      je3 = je + 1
    ELSE
      IF (3 .LT. js - 1) THEN
        js3 = js - 1
      ELSE
        js3 = 3
      END IF
      IF (npy - 3 .GT. je + 1) THEN
        CALL PUSHCONTROL1B(1)
        je3 = je + 1
      ELSE
        CALL PUSHCONTROL1B(1)
        je3 = npy - 3
      END IF
    END IF
    IF (jord .EQ. 1) THEN
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
            flux_ad(i, j) = 0.0
          ELSE
            v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
            flux_ad(i, j) = 0.0
          END IF
        END DO
      END DO
    ELSE IF (jord .LT. 8) THEN
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
      DO j=js3,je3+1
        DO i=is,ie+1
          al(i, j) = p1*(v(i, j-1)+v(i, j)) + p2*(v(i, j-2)+v(i, j+1))
        END DO
      END DO
      DO j=js3,je3
        DO i=is,ie+1
          bl(i, j) = al(i, j) - v(i, j)
          br(i, j) = al(i, j+1) - v(i, j)
        END DO
      END DO
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
            bl(i, 0) = c1*v(i, -2) + c2*v(i, -1) + c3*v(i, 0) - v(i, 0)
            xt = 0.5*(((2.*dy(i, 0)+dy(i, -1))*v(i, 0)-dy(i, 0)*v(i, -1)&
&             )/(dy(i, 0)+dy(i, -1))+((2.*dy(i, 1)+dy(i, 2))*v(i, 1)-dy(&
&             i, 1)*v(i, 2))/(dy(i, 1)+dy(i, 2)))
            br(i, 0) = xt - v(i, 0)
            bl(i, 1) = xt - v(i, 1)
            xt = c3*v(i, 1) + c2*v(i, 2) + c1*v(i, 3)
            br(i, 1) = xt - v(i, 1)
            bl(i, 2) = xt - v(i, 2)
            br(i, 2) = al(i, 3) - v(i, 2)
          END DO
          IF (is .EQ. 1) THEN
! out
            bl(1, 0) = 0.
! edge
            br(1, 0) = 0.
! edge
            bl(1, 1) = 0.
! in
            br(1, 1) = 0.
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
! out
            bl(npx, 0) = 0.
! edge
            br(npx, 0) = 0.
! edge
            bl(npx, 1) = 0.
! in
            br(npx, 1) = 0.
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
!      j=2
!      call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
            bl(i, npy-2) = al(i, npy-2) - v(i, npy-2)
            xt = c1*v(i, npy-3) + c2*v(i, npy-2) + c3*v(i, npy-1)
            br(i, npy-2) = xt - v(i, npy-2)
            bl(i, npy-1) = xt - v(i, npy-1)
            xt = 0.5*(((2.*dy(i, npy-1)+dy(i, npy-2))*v(i, npy-1)-dy(i, &
&             npy-1)*v(i, npy-2))/(dy(i, npy-1)+dy(i, npy-2))+((2.*dy(i&
&             , npy)+dy(i, npy+1))*v(i, npy)-dy(i, npy)*v(i, npy+1))/(dy&
&             (i, npy)+dy(i, npy+1)))
            br(i, npy-1) = xt - v(i, npy-1)
            bl(i, npy) = xt - v(i, npy)
            br(i, npy) = c3*v(i, npy) + c2*v(i, npy+1) + c1*v(i, npy+2) &
&             - v(i, npy)
          END DO
          IF (is .EQ. 1) THEN
! in
            bl(1, npy-1) = 0.
! edge
            br(1, npy-1) = 0.
! edge
            bl(1, npy) = 0.
! out
            br(1, npy) = 0.
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
! in
            bl(npx, npy-1) = 0.
! edge
            br(npx, npy-1) = 0.
! edge
            bl(npx, npy) = 0.
! out
            br(npx, npy) = 0.
            CALL PUSHCONTROL2B(3)
          ELSE
            CALL PUSHCONTROL2B(2)
          END IF
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
!      j=npy-2
!      call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
      DO j=js-1,je+1
        DO i=is,ie+1
          b0(i, j) = bl(i, j) + br(i, j)
        END DO
      END DO
      IF (jord .EQ. 2) THEN
! Perfectly linear
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
        END DO
        bl_ad = 0.0
        br_ad = 0.0
        b0_ad = 0.0
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              cfl = c(i, j)*rdy(i, j)
              temp_ad4 = (cfl+1.)*flux_ad(i, j)
              v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
              cfl_ad = b0(i, j)*temp_ad4 + (bl(i, j)+cfl*b0(i, j))*&
&               flux_ad(i, j)
              bl_ad(i, j) = bl_ad(i, j) + temp_ad4
              b0_ad(i, j) = b0_ad(i, j) + cfl*temp_ad4
              flux_ad(i, j) = 0.0
              c_ad(i, j) = c_ad(i, j) + rdy(i, j)*cfl_ad
            ELSE
              cfl = c(i, j)*rdy(i, j-1)
              temp_ad3 = (1.-cfl)*flux_ad(i, j)
              v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
              cfl_ad = -(b0(i, j-1)*temp_ad3) - (br(i, j-1)-cfl*b0(i, j-&
&               1))*flux_ad(i, j)
              br_ad(i, j-1) = br_ad(i, j-1) + temp_ad3
              b0_ad(i, j-1) = b0_ad(i, j-1) - cfl*temp_ad3
              flux_ad(i, j) = 0.0
              c_ad(i, j) = c_ad(i, j) + rdy(i, j-1)*cfl_ad
            END IF
          END DO
        END DO
      ELSE IF (jord .EQ. 3) THEN
        DO j=js-1,je+1
          DO i=is,ie+1
            IF (b0(i, j) .GE. 0.) THEN
              x0 = b0(i, j)
            ELSE
              x0 = -b0(i, j)
            END IF
            IF (bl(i, j) - br(i, j) .GE. 0.) THEN
              x1 = bl(i, j) - br(i, j)
            ELSE
              x1 = -(bl(i, j)-br(i, j))
            END IF
            smt5(i, j) = x0 .LT. x1
            smt6(i, j) = 3.*x0 .LT. x1
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            CALL PUSHREALARRAY_ADM(fx0(i))
            fx0(i) = 0.
          END DO
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdy(i, j-1)
              IF (smt6(i, j-1) .OR. smt5(i, j)) THEN
                CALL PUSHREALARRAY_ADM(fx0(i))
                fx0(i) = br(i, j-1) - cfl*b0(i, j-1)
                CALL PUSHCONTROL2B(0)
              ELSE IF (smt5(i, j-1)) THEN
                IF (bl(i, j-1) .GE. 0.) THEN
                  x2 = bl(i, j-1)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  x2 = -bl(i, j-1)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (br(i, j-1) .GE. 0.) THEN
                  y1 = br(i, j-1)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  y1 = -br(i, j-1)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (x2 .GT. y1) THEN
                  CALL PUSHREALARRAY_ADM(min1)
                  min1 = y1
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHREALARRAY_ADM(min1)
                  min1 = x2
                  CALL PUSHCONTROL1B(1)
                END IF
! piece-wise linear
                CALL PUSHREALARRAY_ADM(fx0(i))
                fx0(i) = SIGN(min1, br(i, j-1))
                CALL PUSHCONTROL2B(1)
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
              CALL PUSHCONTROL1B(1)
            ELSE
              cfl = c(i, j)*rdy(i, j)
              IF (smt6(i, j) .OR. smt5(i, j-1)) THEN
                CALL PUSHREALARRAY_ADM(fx0(i))
                fx0(i) = bl(i, j) + cfl*b0(i, j)
                CALL PUSHCONTROL2B(0)
              ELSE IF (smt5(i, j)) THEN
                IF (bl(i, j) .GE. 0.) THEN
                  x3 = bl(i, j)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  x3 = -bl(i, j)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (br(i, j) .GE. 0.) THEN
                  y2 = br(i, j)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  y2 = -br(i, j)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (x3 .GT. y2) THEN
                  CALL PUSHREALARRAY_ADM(min2)
                  min2 = y2
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHREALARRAY_ADM(min2)
                  min2 = x3
                  CALL PUSHCONTROL1B(1)
                END IF
                CALL PUSHREALARRAY_ADM(fx0(i))
                fx0(i) = SIGN(min2, bl(i, j))
                CALL PUSHCONTROL2B(1)
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
        END DO
        bl_ad = 0.0
        br_ad = 0.0
        b0_ad = 0.0
        fx0_ad = 0.0
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              cfl = c(i, j)*rdy(i, j)
              v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
              cfl_ad = fx0(i)*flux_ad(i, j)
              fx0_ad(i) = fx0_ad(i) + (cfl+1.)*flux_ad(i, j)
              flux_ad(i, j) = 0.0
              CALL POPCONTROL2B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(fx0(i))
                bl_ad(i, j) = bl_ad(i, j) + fx0_ad(i)
                cfl_ad = cfl_ad + b0(i, j)*fx0_ad(i)
                b0_ad(i, j) = b0_ad(i, j) + cfl*fx0_ad(i)
                fx0_ad(i) = 0.0
              ELSE IF (branch .EQ. 1) THEN
                CALL POPREALARRAY_ADM(fx0(i))
                min2_ad = SIGN(1.d0, min2*bl(i, j))*fx0_ad(i)
                fx0_ad(i) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(min2)
                  y2_ad = min2_ad
                  x3_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(min2)
                  x3_ad = min2_ad
                  y2_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  br_ad(i, j) = br_ad(i, j) + y2_ad
                ELSE
                  br_ad(i, j) = br_ad(i, j) - y2_ad
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  bl_ad(i, j) = bl_ad(i, j) + x3_ad
                ELSE
                  bl_ad(i, j) = bl_ad(i, j) - x3_ad
                END IF
              END IF
              c_ad(i, j) = c_ad(i, j) + rdy(i, j)*cfl_ad
            ELSE
              cfl = c(i, j)*rdy(i, j-1)
              v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
              cfl_ad = -(fx0(i)*flux_ad(i, j))
              fx0_ad(i) = fx0_ad(i) + (1.-cfl)*flux_ad(i, j)
              flux_ad(i, j) = 0.0
              CALL POPCONTROL2B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(fx0(i))
                br_ad(i, j-1) = br_ad(i, j-1) + fx0_ad(i)
                cfl_ad = cfl_ad - b0(i, j-1)*fx0_ad(i)
                b0_ad(i, j-1) = b0_ad(i, j-1) - cfl*fx0_ad(i)
                fx0_ad(i) = 0.0
              ELSE IF (branch .EQ. 1) THEN
                CALL POPREALARRAY_ADM(fx0(i))
                min1_ad = SIGN(1.d0, min1*br(i, j-1))*fx0_ad(i)
                fx0_ad(i) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(min1)
                  y1_ad = min1_ad
                  x2_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(min1)
                  x2_ad = min1_ad
                  y1_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  br_ad(i, j-1) = br_ad(i, j-1) + y1_ad
                ELSE
                  br_ad(i, j-1) = br_ad(i, j-1) - y1_ad
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  bl_ad(i, j-1) = bl_ad(i, j-1) + x2_ad
                ELSE
                  bl_ad(i, j-1) = bl_ad(i, j-1) - x2_ad
                END IF
              END IF
              c_ad(i, j) = c_ad(i, j) + rdy(i, j-1)*cfl_ad
            END IF
          END DO
          DO i=ie+1,is,-1
            CALL POPREALARRAY_ADM(fx0(i))
            fx0_ad(i) = 0.0
          END DO
        END DO
        DO j=je+1,js-1,-1
          i = is - 1
        END DO
      ELSE IF (jord .EQ. 4) THEN
        DO j=js-1,je+1
          DO i=is,ie+1
            IF (b0(i, j) .GE. 0.) THEN
              x0 = b0(i, j)
            ELSE
              x0 = -b0(i, j)
            END IF
            IF (bl(i, j) - br(i, j) .GE. 0.) THEN
              x1 = bl(i, j) - br(i, j)
            ELSE
              x1 = -(bl(i, j)-br(i, j))
            END IF
            smt5(i, j) = x0 .LT. x1
            smt6(i, j) = 3.*x0 .LT. x1
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              IF (smt6(i, j-1) .OR. smt5(i, j)) THEN
                CALL PUSHCONTROL2B(3)
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
            ELSE IF (smt6(i, j) .OR. smt5(i, j-1)) THEN
              CALL PUSHCONTROL2B(1)
            ELSE
              CALL PUSHCONTROL2B(0)
            END IF
          END DO
        END DO
        bl_ad = 0.0
        br_ad = 0.0
        b0_ad = 0.0
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
                flux_ad(i, j) = 0.0
              ELSE
                cfl = c(i, j)*rdy(i, j)
                temp_ad6 = (cfl+1.)*flux_ad(i, j)
                v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
                cfl_ad = b0(i, j)*temp_ad6 + (bl(i, j)+cfl*b0(i, j))*&
&                 flux_ad(i, j)
                bl_ad(i, j) = bl_ad(i, j) + temp_ad6
                b0_ad(i, j) = b0_ad(i, j) + cfl*temp_ad6
                flux_ad(i, j) = 0.0
                c_ad(i, j) = c_ad(i, j) + rdy(i, j)*cfl_ad
              END IF
            ELSE IF (branch .EQ. 2) THEN
              v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              cfl = c(i, j)*rdy(i, j-1)
              temp_ad5 = (1.-cfl)*flux_ad(i, j)
              v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
              cfl_ad = -(b0(i, j-1)*temp_ad5) - (br(i, j-1)-cfl*b0(i, j-&
&               1))*flux_ad(i, j)
              br_ad(i, j-1) = br_ad(i, j-1) + temp_ad5
              b0_ad(i, j-1) = b0_ad(i, j-1) - cfl*temp_ad5
              flux_ad(i, j) = 0.0
              c_ad(i, j) = c_ad(i, j) + rdy(i, j-1)*cfl_ad
            END IF
          END DO
        END DO
        DO j=je+1,js-1,-1
          i = is - 1
        END DO
      ELSE
! jord = 5,6,7
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7
        IF (jord .EQ. 5) THEN
          CALL PUSHCONTROL1B(1)
          DO j=js-1,je+1
            DO i=is,ie+1
              smt5(i, j) = bl(i, j)*br(i, j) .LT. 0.
            END DO
          END DO
        ELSE
! ord = 6, 7
          DO j=js-1,je+1
            DO i=is,ie+1
              IF (3.*b0(i, j) .GE. 0.) THEN
                abs0 = 3.*b0(i, j)
              ELSE
                abs0 = -(3.*b0(i, j))
              END IF
              IF (bl(i, j) - br(i, j) .GE. 0.) THEN
                abs4 = bl(i, j) - br(i, j)
              ELSE
                abs4 = -(bl(i, j)-br(i, j))
              END IF
              smt5(i, j) = abs0 .LT. abs4
            END DO
          END DO
          CALL PUSHCONTROL1B(0)
        END IF
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (smt5(i, j-1) .OR. smt5(i, j)) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
        END DO
        bl_ad = 0.0
        br_ad = 0.0
        b0_ad = 0.0
        fx0_ad = 0.0
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) fx0_ad(i) = fx0_ad(i) + flux_ad(i, j)
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
              cfl = c(i, j)*rdy(i, j-1)
              temp_ad7 = (1.-cfl)*fx0_ad(i)
              cfl_ad = -(b0(i, j-1)*temp_ad7) - (br(i, j-1)-cfl*b0(i, j-&
&               1))*fx0_ad(i)
              br_ad(i, j-1) = br_ad(i, j-1) + temp_ad7
              b0_ad(i, j-1) = b0_ad(i, j-1) - cfl*temp_ad7
              fx0_ad(i) = 0.0
              c_ad(i, j) = c_ad(i, j) + rdy(i, j-1)*cfl_ad
            ELSE
              v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
              cfl = c(i, j)*rdy(i, j)
              temp_ad8 = (cfl+1.)*fx0_ad(i)
              cfl_ad = b0(i, j)*temp_ad8 + (bl(i, j)+cfl*b0(i, j))*&
&               fx0_ad(i)
              bl_ad(i, j) = bl_ad(i, j) + temp_ad8
              b0_ad(i, j) = b0_ad(i, j) + cfl*temp_ad8
              fx0_ad(i) = 0.0
              c_ad(i, j) = c_ad(i, j) + rdy(i, j)*cfl_ad
            END IF
          END DO
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=je+1,js-1,-1
            i = is - 1
          END DO
        END IF
      END IF
      DO j=je+1,js-1,-1
        DO i=ie+1,is,-1
          bl_ad(i, j) = bl_ad(i, j) + b0_ad(i, j)
          br_ad(i, j) = br_ad(i, j) + b0_ad(i, j)
          b0_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          al_ad = 0.0
          GOTO 100
        ELSE
          al_ad = 0.0
        END IF
      ELSE
        IF (branch .NE. 2) THEN
          br_ad(npx, npy) = 0.0
          bl_ad(npx, npy) = 0.0
          br_ad(npx, npy-1) = 0.0
          bl_ad(npx, npy-1) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          br_ad(1, npy) = 0.0
          bl_ad(1, npy) = 0.0
          br_ad(1, npy-1) = 0.0
          bl_ad(1, npy-1) = 0.0
        END IF
        al_ad = 0.0
        DO i=ie+1,is,-1
          v_ad(i, npy) = v_ad(i, npy) + (c3-1.0)*br_ad(i, npy)
          v_ad(i, npy+1) = v_ad(i, npy+1) + c2*br_ad(i, npy)
          v_ad(i, npy+2) = v_ad(i, npy+2) + c1*br_ad(i, npy)
          br_ad(i, npy) = 0.0
          xt_ad = br_ad(i, npy-1) + bl_ad(i, npy)
          v_ad(i, npy) = v_ad(i, npy) - bl_ad(i, npy)
          bl_ad(i, npy) = 0.0
          temp_ad1 = 0.5*xt_ad/(dy(i, npy-1)+dy(i, npy-2))
          v_ad(i, npy-1) = v_ad(i, npy-1) + (dy(i, npy-1)*2.+dy(i, npy-2&
&           ))*temp_ad1 - br_ad(i, npy-1)
          br_ad(i, npy-1) = 0.0
          temp_ad2 = 0.5*xt_ad/(dy(i, npy)+dy(i, npy+1))
          v_ad(i, npy-2) = v_ad(i, npy-2) - dy(i, npy-1)*temp_ad1
          v_ad(i, npy) = v_ad(i, npy) + (dy(i, npy)*2.+dy(i, npy+1))*&
&           temp_ad2
          v_ad(i, npy+1) = v_ad(i, npy+1) - dy(i, npy)*temp_ad2
          xt_ad = br_ad(i, npy-2) + bl_ad(i, npy-1)
          v_ad(i, npy-1) = v_ad(i, npy-1) - bl_ad(i, npy-1)
          bl_ad(i, npy-1) = 0.0
          v_ad(i, npy-2) = v_ad(i, npy-2) - br_ad(i, npy-2)
          br_ad(i, npy-2) = 0.0
          v_ad(i, npy-3) = v_ad(i, npy-3) + c1*xt_ad
          v_ad(i, npy-2) = v_ad(i, npy-2) + c2*xt_ad
          v_ad(i, npy-1) = v_ad(i, npy-1) + c3*xt_ad
          al_ad(i, npy-2) = al_ad(i, npy-2) + bl_ad(i, npy-2)
          v_ad(i, npy-2) = v_ad(i, npy-2) - bl_ad(i, npy-2)
          bl_ad(i, npy-2) = 0.0
        END DO
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        br_ad(npx, 1) = 0.0
        bl_ad(npx, 1) = 0.0
        br_ad(npx, 0) = 0.0
        bl_ad(npx, 0) = 0.0
      ELSE IF (branch .NE. 1) THEN
        GOTO 100
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        br_ad(1, 1) = 0.0
        bl_ad(1, 1) = 0.0
        br_ad(1, 0) = 0.0
        bl_ad(1, 0) = 0.0
      END IF
      DO i=ie+1,is,-1
        al_ad(i, 3) = al_ad(i, 3) + br_ad(i, 2)
        v_ad(i, 2) = v_ad(i, 2) - bl_ad(i, 2) - br_ad(i, 2)
        br_ad(i, 2) = 0.0
        xt_ad = br_ad(i, 1) + bl_ad(i, 2)
        bl_ad(i, 2) = 0.0
        v_ad(i, 1) = v_ad(i, 1) + c3*xt_ad - br_ad(i, 1)
        br_ad(i, 1) = 0.0
        v_ad(i, 2) = v_ad(i, 2) + c2*xt_ad
        v_ad(i, 3) = v_ad(i, 3) + c1*xt_ad
        xt_ad = br_ad(i, 0) + bl_ad(i, 1)
        v_ad(i, 1) = v_ad(i, 1) - bl_ad(i, 1)
        bl_ad(i, 1) = 0.0
        temp_ad = 0.5*xt_ad/(dy(i, 0)+dy(i, -1))
        v_ad(i, 0) = v_ad(i, 0) + (dy(i, 0)*2.+dy(i, -1))*temp_ad - &
&         br_ad(i, 0)
        br_ad(i, 0) = 0.0
        temp_ad0 = 0.5*xt_ad/(dy(i, 1)+dy(i, 2))
        v_ad(i, -1) = v_ad(i, -1) - dy(i, 0)*temp_ad
        v_ad(i, 1) = v_ad(i, 1) + (dy(i, 1)*2.+dy(i, 2))*temp_ad0
        v_ad(i, 2) = v_ad(i, 2) - dy(i, 1)*temp_ad0
        v_ad(i, -2) = v_ad(i, -2) + c1*bl_ad(i, 0)
        v_ad(i, -1) = v_ad(i, -1) + c2*bl_ad(i, 0)
        v_ad(i, 0) = v_ad(i, 0) + (c3-1.0)*bl_ad(i, 0)
        bl_ad(i, 0) = 0.0
      END DO
 100  DO j=je3,js3,-1
        DO i=ie+1,is,-1
          al_ad(i, j+1) = al_ad(i, j+1) + br_ad(i, j)
          v_ad(i, j) = v_ad(i, j) - bl_ad(i, j) - br_ad(i, j)
          br_ad(i, j) = 0.0
          al_ad(i, j) = al_ad(i, j) + bl_ad(i, j)
          bl_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je3+1,js3,-1
        DO i=ie+1,is,-1
          v_ad(i, j-1) = v_ad(i, j-1) + p1*al_ad(i, j)
          v_ad(i, j) = v_ad(i, j) + p1*al_ad(i, j)
          v_ad(i, j-2) = v_ad(i, j-2) + p2*al_ad(i, j)
          v_ad(i, j+1) = v_ad(i, j+1) + p2*al_ad(i, j)
          al_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
! jord= 8, 9, 10
      DO j=js-2,je+2
        DO i=is,ie+1
          xt = 0.25*(v(i, j+1)-v(i, j-1))
          IF (xt .GE. 0.) THEN
            x4 = xt
            CALL PUSHCONTROL1B(0)
          ELSE
            x4 = -xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (v(i, j-1) .LT. v(i, j)) THEN
            IF (v(i, j) .LT. v(i, j+1)) THEN
              max1 = v(i, j+1)
              CALL PUSHCONTROL2B(0)
            ELSE
              max1 = v(i, j)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (v(i, j-1) .LT. v(i, j+1)) THEN
            max1 = v(i, j+1)
            CALL PUSHCONTROL2B(2)
          ELSE
            max1 = v(i, j-1)
            CALL PUSHCONTROL2B(3)
          END IF
          y3 = max1 - v(i, j)
          IF (v(i, j-1) .GT. v(i, j)) THEN
            IF (v(i, j) .GT. v(i, j+1)) THEN
              min6 = v(i, j+1)
              CALL PUSHCONTROL2B(0)
            ELSE
              min6 = v(i, j)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (v(i, j-1) .GT. v(i, j+1)) THEN
            min6 = v(i, j+1)
            CALL PUSHCONTROL2B(2)
          ELSE
            min6 = v(i, j-1)
            CALL PUSHCONTROL2B(3)
          END IF
          z1 = v(i, j) - min6
          IF (x4 .GT. y3) THEN
            IF (y3 .GT. z1) THEN
              CALL PUSHREALARRAY_ADM(min3)
              min3 = z1
              CALL PUSHCONTROL2B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(min3)
              min3 = y3
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (x4 .GT. z1) THEN
            CALL PUSHREALARRAY_ADM(min3)
            min3 = z1
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHREALARRAY_ADM(min3)
            min3 = x4
            CALL PUSHCONTROL2B(3)
          END IF
          dm(i, j) = SIGN(min3, xt)
        END DO
      END DO
      DO j=js-3,je+2
        DO i=is,ie+1
          dq(i, j) = v(i, j+1) - v(i, j)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        DO j=js3,je3+1
          DO i=is,ie+1
            al(i, j) = 0.5*(v(i, j-1)+v(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        IF (jord .EQ. 8) THEN
          DO j=js3,je3
            DO i=is,ie+1
              xt = 2.*dm(i, j)
              IF (xt .GE. 0.) THEN
                x5 = xt
                CALL PUSHCONTROL1B(0)
              ELSE
                x5 = -xt
                CALL PUSHCONTROL1B(1)
              END IF
              IF (al(i, j) - v(i, j) .GE. 0.) THEN
                y4 = al(i, j) - v(i, j)
                CALL PUSHCONTROL1B(0)
              ELSE
                y4 = -(al(i, j)-v(i, j))
                CALL PUSHCONTROL1B(1)
              END IF
              IF (x5 .GT. y4) THEN
                CALL PUSHREALARRAY_ADM(min4)
                min4 = y4
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(min4)
                min4 = x5
                CALL PUSHCONTROL1B(1)
              END IF
              bl(i, j) = -SIGN(min4, xt)
              IF (xt .GE. 0.) THEN
                x6 = xt
                CALL PUSHCONTROL1B(0)
              ELSE
                x6 = -xt
                CALL PUSHCONTROL1B(1)
              END IF
              IF (al(i, j+1) - v(i, j) .GE. 0.) THEN
                y5 = al(i, j+1) - v(i, j)
                CALL PUSHCONTROL1B(0)
              ELSE
                y5 = -(al(i, j+1)-v(i, j))
                CALL PUSHCONTROL1B(1)
              END IF
              IF (x6 .GT. y5) THEN
                CALL PUSHREALARRAY_ADM(min5)
                min5 = y5
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(min5)
                min5 = x6
                CALL PUSHCONTROL1B(1)
              END IF
              br(i, j) = SIGN(min5, xt)
            END DO
          END DO
          CALL PUSHCONTROL2B(0)
        ELSE IF (jord .EQ. 9) THEN
          DO j=js3,je3
            DO i=is,ie+1
              pmp_1 = -(2.*dq(i, j))
              lac_1 = pmp_1 + 1.5*dq(i, j+1)
              IF (0. .LT. pmp_1) THEN
                IF (pmp_1 .LT. lac_1) THEN
                  x7 = lac_1
                  CALL PUSHCONTROL2B(0)
                ELSE
                  x7 = pmp_1
                  CALL PUSHCONTROL2B(1)
                END IF
              ELSE IF (0. .LT. lac_1) THEN
                x7 = lac_1
                CALL PUSHCONTROL2B(2)
              ELSE
                CALL PUSHCONTROL2B(3)
                x7 = 0.
              END IF
              IF (0. .GT. pmp_1) THEN
                IF (pmp_1 .GT. lac_1) THEN
                  y12 = lac_1
                  CALL PUSHCONTROL2B(0)
                ELSE
                  y12 = pmp_1
                  CALL PUSHCONTROL2B(1)
                END IF
              ELSE IF (0. .GT. lac_1) THEN
                y12 = lac_1
                CALL PUSHCONTROL2B(2)
              ELSE
                y12 = 0.
                CALL PUSHCONTROL2B(3)
              END IF
              IF (al(i, j) - v(i, j) .LT. y12) THEN
                y6 = y12
                CALL PUSHCONTROL1B(0)
              ELSE
                y6 = al(i, j) - v(i, j)
                CALL PUSHCONTROL1B(1)
              END IF
              IF (x7 .GT. y6) THEN
                bl(i, j) = y6
                CALL PUSHCONTROL1B(0)
              ELSE
                bl(i, j) = x7
                CALL PUSHCONTROL1B(1)
              END IF
              pmp_2 = 2.*dq(i, j-1)
              lac_2 = pmp_2 - 1.5*dq(i, j-2)
              IF (0. .LT. pmp_2) THEN
                IF (pmp_2 .LT. lac_2) THEN
                  x8 = lac_2
                  CALL PUSHCONTROL2B(0)
                ELSE
                  x8 = pmp_2
                  CALL PUSHCONTROL2B(1)
                END IF
              ELSE IF (0. .LT. lac_2) THEN
                x8 = lac_2
                CALL PUSHCONTROL2B(2)
              ELSE
                CALL PUSHCONTROL2B(3)
                x8 = 0.
              END IF
              IF (0. .GT. pmp_2) THEN
                IF (pmp_2 .GT. lac_2) THEN
                  y13 = lac_2
                  CALL PUSHCONTROL2B(0)
                ELSE
                  y13 = pmp_2
                  CALL PUSHCONTROL2B(1)
                END IF
              ELSE IF (0. .GT. lac_2) THEN
                y13 = lac_2
                CALL PUSHCONTROL2B(2)
              ELSE
                y13 = 0.
                CALL PUSHCONTROL2B(3)
              END IF
              IF (al(i, j+1) - v(i, j) .LT. y13) THEN
                y7 = y13
                CALL PUSHCONTROL1B(0)
              ELSE
                y7 = al(i, j+1) - v(i, j)
                CALL PUSHCONTROL1B(1)
              END IF
              IF (x8 .GT. y7) THEN
                br(i, j) = y7
                CALL PUSHCONTROL1B(0)
              ELSE
                br(i, j) = x8
                CALL PUSHCONTROL1B(1)
              END IF
            END DO
          END DO
          CALL PUSHCONTROL2B(1)
        ELSE IF (jord .EQ. 10) THEN
          DO j=js3,je3
            DO i=is,ie+1
              bl(i, j) = al(i, j) - v(i, j)
              br(i, j) = al(i, j+1) - v(i, j)
              IF (dm(i, j) .GE. 0.) THEN
                abs1 = dm(i, j)
              ELSE
                abs1 = -dm(i, j)
              END IF
!           if ( abs(dm(i,j-1))+abs(dm(i,j))+abs(dm(i,j+1)) < near_zero ) then
              IF (abs1 .LT. near_zero) THEN
                IF (dm(i, j-1) .GE. 0.) THEN
                  abs2 = dm(i, j-1)
                ELSE
                  abs2 = -dm(i, j-1)
                END IF
                IF (dm(i, j+1) .GE. 0.) THEN
                  abs5 = dm(i, j+1)
                ELSE
                  abs5 = -dm(i, j+1)
                END IF
                IF (abs2 + abs5 .LT. near_zero) THEN
                  bl(i, j) = 0.
                  br(i, j) = 0.
                  CALL PUSHCONTROL3B(4)
                ELSE
                  CALL PUSHCONTROL3B(3)
                END IF
              ELSE
                IF (3.*(bl(i, j)+br(i, j)) .GE. 0.) THEN
                  abs3 = 3.*(bl(i, j)+br(i, j))
                ELSE
                  abs3 = -(3.*(bl(i, j)+br(i, j)))
                END IF
                IF (bl(i, j) - br(i, j) .GE. 0.) THEN
                  abs6 = bl(i, j) - br(i, j)
                ELSE
                  abs6 = -(bl(i, j)-br(i, j))
                END IF
                IF (abs3 .GT. abs6) THEN
                  pmp_1 = -(2.*dq(i, j))
                  lac_1 = pmp_1 + 1.5*dq(i, j+1)
                  IF (0. .LT. pmp_1) THEN
                    IF (pmp_1 .LT. lac_1) THEN
                      x9 = lac_1
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      x9 = pmp_1
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (0. .LT. lac_1) THEN
                    x9 = lac_1
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    CALL PUSHCONTROL2B(3)
                    x9 = 0.
                  END IF
                  IF (0. .GT. pmp_1) THEN
                    IF (pmp_1 .GT. lac_1) THEN
                      y14 = lac_1
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y14 = pmp_1
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (0. .GT. lac_1) THEN
                    y14 = lac_1
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y14 = 0.
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (bl(i, j) .LT. y14) THEN
                    y8 = y14
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    y8 = bl(i, j)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (x9 .GT. y8) THEN
                    bl(i, j) = y8
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    bl(i, j) = x9
                    CALL PUSHCONTROL1B(1)
                  END IF
                  pmp_2 = 2.*dq(i, j-1)
                  lac_2 = pmp_2 - 1.5*dq(i, j-2)
                  IF (0. .LT. pmp_2) THEN
                    IF (pmp_2 .LT. lac_2) THEN
                      x10 = lac_2
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      x10 = pmp_2
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (0. .LT. lac_2) THEN
                    x10 = lac_2
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    CALL PUSHCONTROL2B(3)
                    x10 = 0.
                  END IF
                  IF (0. .GT. pmp_2) THEN
                    IF (pmp_2 .GT. lac_2) THEN
                      y15 = lac_2
                      CALL PUSHCONTROL2B(0)
                    ELSE
                      y15 = pmp_2
                      CALL PUSHCONTROL2B(1)
                    END IF
                  ELSE IF (0. .GT. lac_2) THEN
                    y15 = lac_2
                    CALL PUSHCONTROL2B(2)
                  ELSE
                    y15 = 0.
                    CALL PUSHCONTROL2B(3)
                  END IF
                  IF (br(i, j) .LT. y15) THEN
                    y9 = y15
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    y9 = br(i, j)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (x10 .GT. y9) THEN
                    br(i, j) = y9
                    CALL PUSHCONTROL3B(1)
                  ELSE
                    br(i, j) = x10
                    CALL PUSHCONTROL3B(2)
                  END IF
                ELSE
                  CALL PUSHCONTROL3B(0)
                END IF
              END IF
            END DO
          END DO
          CALL PUSHCONTROL2B(2)
        ELSE
! Unlimited:
          DO j=js3,je3
            DO i=is,ie+1
              bl(i, j) = al(i, j) - v(i, j)
              br(i, j) = al(i, j+1) - v(i, j)
            END DO
          END DO
          CALL PUSHCONTROL2B(3)
        END IF
!--------------
! fix the edges
!--------------
        IF (js .EQ. 1 .AND. (.NOT.nested)) THEN
          DO i=is,ie+1
            br(i, 2) = al(i, 3) - v(i, 2)
            xt = s15*v(i, 1) + s11*v(i, 2) - s14*dm(i, 2)
            br(i, 1) = xt - v(i, 1)
            bl(i, 2) = xt - v(i, 2)
            bl(i, 0) = s14*dm(i, -1) - s11*dq(i, -1)
            x0l = 0.5*((2.*dy(i, 0)+dy(i, -1))*v(i, 0)-dy(i, 0)*v(i, -1)&
&             )/(dy(i, 0)+dy(i, -1))
            x0r = 0.5*((2.*dy(i, 1)+dy(i, 2))*v(i, 1)-dy(i, 1)*v(i, 2))/&
&             (dy(i, 1)+dy(i, 2))
            xt = x0l + x0r
            bl(i, 1) = xt - v(i, 1)
            br(i, 0) = xt - v(i, 0)
          END DO
          IF (is .EQ. 1) THEN
! out
            bl(1, 0) = 0.
! edge
            br(1, 0) = 0.
! edge
            bl(1, 1) = 0.
! in
            br(1, 1) = 0.
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
! out
            bl(npx, 0) = 0.
! edge
            br(npx, 0) = 0.
! edge
            bl(npx, 1) = 0.
! in
            br(npx, 1) = 0.
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          j = 2
          arg1 = ie - is + 2
          CALL PUSHREALARRAY_ADM(br(:, j), ie - is + 2)
          CALL PUSHREALARRAY_ADM(bl(:, j), ie - is + 2)
          CALL PERT_PPM(arg1, v(is:ie+1, j), bl(is:ie+1, j), br(is:ie+1&
&                 , j), -1)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (je + 1 .EQ. npy .AND. (.NOT.nested)) THEN
          DO i=is,ie+1
            bl(i, npy-2) = al(i, npy-2) - v(i, npy-2)
            xt = s15*v(i, npy-1) + s11*v(i, npy-2) + s14*dm(i, npy-2)
            br(i, npy-2) = xt - v(i, npy-2)
            bl(i, npy-1) = xt - v(i, npy-1)
            br(i, npy) = s11*dq(i, npy) - s14*dm(i, npy+1)
            x0l = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*v(i, npy-1)-dy(i, &
&             npy-1)*v(i, npy-2))/(dy(i, npy-1)+dy(i, npy-2))
            x0r = 0.5*((2.*dy(i, npy)+dy(i, npy+1))*v(i, npy)-dy(i, npy)&
&             *v(i, npy+1))/(dy(i, npy)+dy(i, npy+1))
            xt = x0l + x0r
            br(i, npy-1) = xt - v(i, npy-1)
            bl(i, npy) = xt - v(i, npy)
          END DO
          IF (is .EQ. 1) THEN
! in
            bl(1, npy-1) = 0.
! edge
            br(1, npy-1) = 0.
! edge
            bl(1, npy) = 0.
! out
            br(1, npy) = 0.
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
! in
            bl(npx, npy-1) = 0.
! edge
            br(npx, npy-1) = 0.
! edge
            bl(npx, npy) = 0.
! out
            br(npx, npy) = 0.
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          j = npy - 2
          arg1 = ie - is + 2
          CALL PUSHREALARRAY_ADM(br(:, j), ie - is + 2)
          CALL PUSHREALARRAY_ADM(bl(:, j), ie - is + 2)
          CALL PERT_PPM(arg1, v(is:ie+1, j), bl(is:ie+1, j), br(is:ie+1&
&                 , j), -1)
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        DO j=js-1,je+2
          DO i=is,ie+1
            al(i, j) = 0.5*(v(i, j-1)+v(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        DO j=js-1,je+1
          DO i=is,ie+1
            pmp = 2.*dq(i, j-1)
            lac = pmp - 1.5*dq(i, j-2)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x11 = lac
                CALL PUSHCONTROL2B(0)
              ELSE
                x11 = pmp
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (0. .LT. lac) THEN
              x11 = lac
              CALL PUSHCONTROL2B(2)
            ELSE
              CALL PUSHCONTROL2B(3)
              x11 = 0.
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y16 = lac
                CALL PUSHCONTROL2B(0)
              ELSE
                y16 = pmp
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (0. .GT. lac) THEN
              y16 = lac
              CALL PUSHCONTROL2B(2)
            ELSE
              y16 = 0.
              CALL PUSHCONTROL2B(3)
            END IF
            IF (al(i, j+1) - v(i, j) .LT. y16) THEN
              y10 = y16
              CALL PUSHCONTROL1B(0)
            ELSE
              y10 = al(i, j+1) - v(i, j)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x11 .GT. y10) THEN
              br(i, j) = y10
              CALL PUSHCONTROL1B(0)
            ELSE
              br(i, j) = x11
              CALL PUSHCONTROL1B(1)
            END IF
            pmp = -(2.*dq(i, j))
            lac = pmp + 1.5*dq(i, j+1)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x12 = lac
                CALL PUSHCONTROL2B(0)
              ELSE
                x12 = pmp
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (0. .LT. lac) THEN
              x12 = lac
              CALL PUSHCONTROL2B(2)
            ELSE
              CALL PUSHCONTROL2B(3)
              x12 = 0.
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y17 = lac
                CALL PUSHCONTROL2B(0)
              ELSE
                y17 = pmp
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (0. .GT. lac) THEN
              y17 = lac
              CALL PUSHCONTROL2B(2)
            ELSE
              y17 = 0.
              CALL PUSHCONTROL2B(3)
            END IF
            IF (al(i, j) - v(i, j) .LT. y17) THEN
              y11 = y17
              CALL PUSHCONTROL1B(0)
            ELSE
              y11 = al(i, j) - v(i, j)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x12 .GT. y11) THEN
              bl(i, j) = y11
              CALL PUSHCONTROL1B(0)
            ELSE
              bl(i, j) = x12
              CALL PUSHCONTROL1B(1)
            END IF
          END DO
        END DO
        CALL PUSHCONTROL2B(0)
      END IF
      CALL PUSHINTEGER4(j)
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      bl_ad = 0.0
      br_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            cfl = c(i, j)*rdy(i, j)
            temp0 = bl(i, j) + br(i, j)
            temp_ad14 = (cfl+1.)*flux_ad(i, j)
            v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
            cfl_ad = temp0*temp_ad14 + (bl(i, j)+cfl*temp0)*flux_ad(i, j&
&             )
            bl_ad(i, j) = bl_ad(i, j) + (cfl+1.0)*temp_ad14
            br_ad(i, j) = br_ad(i, j) + cfl*temp_ad14
            flux_ad(i, j) = 0.0
            c_ad(i, j) = c_ad(i, j) + rdy(i, j)*cfl_ad
          ELSE
            cfl = c(i, j)*rdy(i, j-1)
            temp = bl(i, j-1) + br(i, j-1)
            temp_ad13 = (1.-cfl)*flux_ad(i, j)
            v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
            cfl_ad = -(temp*temp_ad13) - (br(i, j-1)-cfl*temp)*flux_ad(i&
&             , j)
            br_ad(i, j-1) = br_ad(i, j-1) + (1.0-cfl)*temp_ad13
            bl_ad(i, j-1) = bl_ad(i, j-1) - cfl*temp_ad13
            flux_ad(i, j) = 0.0
            c_ad(i, j) = c_ad(i, j) + rdy(i, j-1)*cfl_ad
          END IF
        END DO
      END DO
      CALL POPINTEGER4(j)
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        dq_ad = 0.0
        al_ad = 0.0
        DO j=je+1,js-1,-1
          DO i=ie+1,is,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y11_ad = bl_ad(i, j)
              bl_ad(i, j) = 0.0
              x12_ad = 0.0
            ELSE
              x12_ad = bl_ad(i, j)
              bl_ad(i, j) = 0.0
              y11_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y17_ad = y11_ad
            ELSE
              al_ad(i, j) = al_ad(i, j) + y11_ad
              v_ad(i, j) = v_ad(i, j) - y11_ad
              y17_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_ad = y17_ad
                pmp_ad = 0.0
              ELSE
                pmp_ad = y17_ad
                lac_ad = 0.0
              END IF
            ELSE
              IF (branch .EQ. 2) THEN
                lac_ad = y17_ad
              ELSE
                lac_ad = 0.0
              END IF
              pmp_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_ad = lac_ad + x12_ad
              ELSE
                pmp_ad = pmp_ad + x12_ad
              END IF
            ELSE IF (branch .EQ. 2) THEN
              lac_ad = lac_ad + x12_ad
            END IF
            pmp_ad = pmp_ad + lac_ad
            dq_ad(i, j+1) = dq_ad(i, j+1) + 1.5*lac_ad
            dq_ad(i, j) = dq_ad(i, j) - 2.*pmp_ad
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y10_ad = br_ad(i, j)
              br_ad(i, j) = 0.0
              x11_ad = 0.0
            ELSE
              x11_ad = br_ad(i, j)
              br_ad(i, j) = 0.0
              y10_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y16_ad = y10_ad
            ELSE
              al_ad(i, j+1) = al_ad(i, j+1) + y10_ad
              v_ad(i, j) = v_ad(i, j) - y10_ad
              y16_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_ad = y16_ad
                pmp_ad = 0.0
              ELSE
                pmp_ad = y16_ad
                lac_ad = 0.0
              END IF
            ELSE
              IF (branch .EQ. 2) THEN
                lac_ad = y16_ad
              ELSE
                lac_ad = 0.0
              END IF
              pmp_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_ad = lac_ad + x11_ad
              ELSE
                pmp_ad = pmp_ad + x11_ad
              END IF
            ELSE IF (branch .EQ. 2) THEN
              lac_ad = lac_ad + x11_ad
            END IF
            pmp_ad = pmp_ad + lac_ad
            dq_ad(i, j-2) = dq_ad(i, j-2) - 1.5*lac_ad
            dq_ad(i, j-1) = dq_ad(i, j-1) + 2.*pmp_ad
          END DO
        END DO
        dm_ad = 0.0
        DO j=je+2,js-1,-1
          DO i=ie+1,is,-1
            v_ad(i, j-1) = v_ad(i, j-1) + 0.5*al_ad(i, j)
            v_ad(i, j) = v_ad(i, j) + 0.5*al_ad(i, j)
            dm_ad(i, j-1) = dm_ad(i, j-1) + r3*al_ad(i, j)
            dm_ad(i, j) = dm_ad(i, j) - r3*al_ad(i, j)
            al_ad(i, j) = 0.0
          END DO
        END DO
      ELSE
        IF (branch .EQ. 1) THEN
          dm_ad = 0.0
          dq_ad = 0.0
          al_ad = 0.0
        ELSE
          j = npy - 2
          CALL POPREALARRAY_ADM(bl(:, j), ie - is + 2)
          CALL POPREALARRAY_ADM(br(:, j), ie - is + 2)
          CALL PERT_PPM_ADM(arg1, v(is:ie+1, j), bl(is:ie+1, j), bl_ad(&
&                     is:ie+1, j), br(is:ie+1, j), br_ad(is:ie+1, j), -1&
&                    )
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            br_ad(npx, npy) = 0.0
            bl_ad(npx, npy) = 0.0
            br_ad(npx, npy-1) = 0.0
            bl_ad(npx, npy-1) = 0.0
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            br_ad(1, npy) = 0.0
            bl_ad(1, npy) = 0.0
            br_ad(1, npy-1) = 0.0
            bl_ad(1, npy-1) = 0.0
          END IF
          dm_ad = 0.0
          dq_ad = 0.0
          al_ad = 0.0
          DO i=ie+1,is,-1
            xt_ad = br_ad(i, npy-1) + bl_ad(i, npy)
            v_ad(i, npy) = v_ad(i, npy) - bl_ad(i, npy)
            bl_ad(i, npy) = 0.0
            v_ad(i, npy-1) = v_ad(i, npy-1) - br_ad(i, npy-1)
            br_ad(i, npy-1) = 0.0
            x0l_ad = xt_ad
            x0r_ad = xt_ad
            temp_ad11 = 0.5*x0r_ad/(dy(i, npy)+dy(i, npy+1))
            v_ad(i, npy) = v_ad(i, npy) + (dy(i, npy)*2.+dy(i, npy+1))*&
&             temp_ad11
            v_ad(i, npy+1) = v_ad(i, npy+1) - dy(i, npy)*temp_ad11
            temp_ad12 = 0.5*x0l_ad/(dy(i, npy-1)+dy(i, npy-2))
            v_ad(i, npy-1) = v_ad(i, npy-1) + (dy(i, npy-1)*2.+dy(i, npy&
&             -2))*temp_ad12
            v_ad(i, npy-2) = v_ad(i, npy-2) - dy(i, npy-1)*temp_ad12
            dq_ad(i, npy) = dq_ad(i, npy) + s11*br_ad(i, npy)
            dm_ad(i, npy+1) = dm_ad(i, npy+1) - s14*br_ad(i, npy)
            br_ad(i, npy) = 0.0
            xt_ad = br_ad(i, npy-2) + bl_ad(i, npy-1)
            v_ad(i, npy-1) = v_ad(i, npy-1) - bl_ad(i, npy-1)
            bl_ad(i, npy-1) = 0.0
            v_ad(i, npy-2) = v_ad(i, npy-2) - br_ad(i, npy-2)
            br_ad(i, npy-2) = 0.0
            v_ad(i, npy-1) = v_ad(i, npy-1) + s15*xt_ad
            v_ad(i, npy-2) = v_ad(i, npy-2) + s11*xt_ad - bl_ad(i, npy-2&
&             )
            dm_ad(i, npy-2) = dm_ad(i, npy-2) + s14*xt_ad
            al_ad(i, npy-2) = al_ad(i, npy-2) + bl_ad(i, npy-2)
            bl_ad(i, npy-2) = 0.0
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          j = 2
          arg1 = ie - is + 2
          CALL POPREALARRAY_ADM(bl(:, j), ie - is + 2)
          CALL POPREALARRAY_ADM(br(:, j), ie - is + 2)
          CALL PERT_PPM_ADM(arg1, v(is:ie+1, j), bl(is:ie+1, j), bl_ad(&
&                     is:ie+1, j), br(is:ie+1, j), br_ad(is:ie+1, j), -1&
&                    )
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            br_ad(npx, 1) = 0.0
            bl_ad(npx, 1) = 0.0
            br_ad(npx, 0) = 0.0
            bl_ad(npx, 0) = 0.0
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            br_ad(1, 1) = 0.0
            bl_ad(1, 1) = 0.0
            br_ad(1, 0) = 0.0
            bl_ad(1, 0) = 0.0
          END IF
          DO i=ie+1,is,-1
            xt_ad = bl_ad(i, 1) + br_ad(i, 0)
            v_ad(i, 0) = v_ad(i, 0) - br_ad(i, 0)
            br_ad(i, 0) = 0.0
            x0l_ad = xt_ad
            x0r_ad = xt_ad
            temp_ad9 = 0.5*x0r_ad/(dy(i, 1)+dy(i, 2))
            v_ad(i, 1) = v_ad(i, 1) + (dy(i, 1)*2.+dy(i, 2))*temp_ad9 - &
&             bl_ad(i, 1)
            bl_ad(i, 1) = 0.0
            v_ad(i, 2) = v_ad(i, 2) - dy(i, 1)*temp_ad9
            temp_ad10 = 0.5*x0l_ad/(dy(i, 0)+dy(i, -1))
            v_ad(i, 0) = v_ad(i, 0) + (dy(i, 0)*2.+dy(i, -1))*temp_ad10
            v_ad(i, -1) = v_ad(i, -1) - dy(i, 0)*temp_ad10
            dm_ad(i, -1) = dm_ad(i, -1) + s14*bl_ad(i, 0)
            dq_ad(i, -1) = dq_ad(i, -1) - s11*bl_ad(i, 0)
            bl_ad(i, 0) = 0.0
            xt_ad = br_ad(i, 1) + bl_ad(i, 2)
            v_ad(i, 2) = v_ad(i, 2) - bl_ad(i, 2)
            bl_ad(i, 2) = 0.0
            v_ad(i, 1) = v_ad(i, 1) + s15*xt_ad - br_ad(i, 1)
            br_ad(i, 1) = 0.0
            v_ad(i, 2) = v_ad(i, 2) + s11*xt_ad - br_ad(i, 2)
            dm_ad(i, 2) = dm_ad(i, 2) - s14*xt_ad
            al_ad(i, 3) = al_ad(i, 3) + br_ad(i, 2)
            br_ad(i, 2) = 0.0
          END DO
        END IF
        CALL POPCONTROL2B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            DO j=je3,js3,-1
              DO i=ie+1,is,-1
                xt = 2.*dm(i, j)
                min5_ad = SIGN(1.d0, min5*xt)*br_ad(i, j)
                br_ad(i, j) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(min5)
                  y5_ad = min5_ad
                  x6_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(min5)
                  x6_ad = min5_ad
                  y5_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  al_ad(i, j+1) = al_ad(i, j+1) + y5_ad
                  v_ad(i, j) = v_ad(i, j) - y5_ad
                ELSE
                  v_ad(i, j) = v_ad(i, j) + y5_ad
                  al_ad(i, j+1) = al_ad(i, j+1) - y5_ad
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  xt_ad = x6_ad
                ELSE
                  xt_ad = -x6_ad
                END IF
                min4_ad = -(SIGN(1.d0, min4*xt)*bl_ad(i, j))
                bl_ad(i, j) = 0.0
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY_ADM(min4)
                  y4_ad = min4_ad
                  x5_ad = 0.0
                ELSE
                  CALL POPREALARRAY_ADM(min4)
                  x5_ad = min4_ad
                  y4_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  al_ad(i, j) = al_ad(i, j) + y4_ad
                  v_ad(i, j) = v_ad(i, j) - y4_ad
                ELSE
                  v_ad(i, j) = v_ad(i, j) + y4_ad
                  al_ad(i, j) = al_ad(i, j) - y4_ad
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  xt_ad = xt_ad + x5_ad
                ELSE
                  xt_ad = xt_ad - x5_ad
                END IF
                dm_ad(i, j) = dm_ad(i, j) + 2.*xt_ad
              END DO
            END DO
          ELSE
            DO j=je3,js3,-1
              DO i=ie+1,is,-1
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y7_ad = br_ad(i, j)
                  br_ad(i, j) = 0.0
                  x8_ad = 0.0
                ELSE
                  x8_ad = br_ad(i, j)
                  br_ad(i, j) = 0.0
                  y7_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y13_ad = y7_ad
                ELSE
                  al_ad(i, j+1) = al_ad(i, j+1) + y7_ad
                  v_ad(i, j) = v_ad(i, j) - y7_ad
                  y13_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = y13_ad
                    pmp_2_ad = 0.0
                  ELSE
                    pmp_2_ad = y13_ad
                    lac_2_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_2_ad = y13_ad
                  ELSE
                    lac_2_ad = 0.0
                  END IF
                  pmp_2_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_2_ad = lac_2_ad + x8_ad
                  ELSE
                    pmp_2_ad = pmp_2_ad + x8_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_2_ad = lac_2_ad + x8_ad
                END IF
                pmp_2_ad = pmp_2_ad + lac_2_ad
                dq_ad(i, j-2) = dq_ad(i, j-2) - 1.5*lac_2_ad
                dq_ad(i, j-1) = dq_ad(i, j-1) + 2.*pmp_2_ad
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y6_ad = bl_ad(i, j)
                  bl_ad(i, j) = 0.0
                  x7_ad = 0.0
                ELSE
                  x7_ad = bl_ad(i, j)
                  bl_ad(i, j) = 0.0
                  y6_ad = 0.0
                END IF
                CALL POPCONTROL1B(branch)
                IF (branch .EQ. 0) THEN
                  y12_ad = y6_ad
                ELSE
                  al_ad(i, j) = al_ad(i, j) + y6_ad
                  v_ad(i, j) = v_ad(i, j) - y6_ad
                  y12_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = y12_ad
                    pmp_1_ad = 0.0
                  ELSE
                    pmp_1_ad = y12_ad
                    lac_1_ad = 0.0
                  END IF
                ELSE
                  IF (branch .EQ. 2) THEN
                    lac_1_ad = y12_ad
                  ELSE
                    lac_1_ad = 0.0
                  END IF
                  pmp_1_ad = 0.0
                END IF
                CALL POPCONTROL2B(branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    lac_1_ad = lac_1_ad + x7_ad
                  ELSE
                    pmp_1_ad = pmp_1_ad + x7_ad
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  lac_1_ad = lac_1_ad + x7_ad
                END IF
                pmp_1_ad = pmp_1_ad + lac_1_ad
                dq_ad(i, j+1) = dq_ad(i, j+1) + 1.5*lac_1_ad
                dq_ad(i, j) = dq_ad(i, j) - 2.*pmp_1_ad
              END DO
            END DO
          END IF
        ELSE IF (branch .EQ. 2) THEN
          DO j=je3,js3,-1
            DO i=ie+1,is,-1
              CALL POPCONTROL3B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  GOTO 110
                ELSE
                  y9_ad = br_ad(i, j)
                  br_ad(i, j) = 0.0
                  x10_ad = 0.0
                END IF
              ELSE IF (branch .EQ. 2) THEN
                x10_ad = br_ad(i, j)
                br_ad(i, j) = 0.0
                y9_ad = 0.0
              ELSE
                IF (branch .NE. 3) THEN
                  br_ad(i, j) = 0.0
                  bl_ad(i, j) = 0.0
                END IF
                GOTO 110
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                y15_ad = y9_ad
              ELSE
                br_ad(i, j) = br_ad(i, j) + y9_ad
                y15_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_2_ad = y15_ad
                  pmp_2_ad = 0.0
                ELSE
                  pmp_2_ad = y15_ad
                  lac_2_ad = 0.0
                END IF
              ELSE
                IF (branch .EQ. 2) THEN
                  lac_2_ad = y15_ad
                ELSE
                  lac_2_ad = 0.0
                END IF
                pmp_2_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_2_ad = lac_2_ad + x10_ad
                ELSE
                  pmp_2_ad = pmp_2_ad + x10_ad
                END IF
              ELSE IF (branch .EQ. 2) THEN
                lac_2_ad = lac_2_ad + x10_ad
              END IF
              pmp_2_ad = pmp_2_ad + lac_2_ad
              dq_ad(i, j-2) = dq_ad(i, j-2) - 1.5*lac_2_ad
              dq_ad(i, j-1) = dq_ad(i, j-1) + 2.*pmp_2_ad
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                y8_ad = bl_ad(i, j)
                bl_ad(i, j) = 0.0
                x9_ad = 0.0
              ELSE
                x9_ad = bl_ad(i, j)
                bl_ad(i, j) = 0.0
                y8_ad = 0.0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                y14_ad = y8_ad
              ELSE
                bl_ad(i, j) = bl_ad(i, j) + y8_ad
                y14_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_1_ad = y14_ad
                  pmp_1_ad = 0.0
                ELSE
                  pmp_1_ad = y14_ad
                  lac_1_ad = 0.0
                END IF
              ELSE
                IF (branch .EQ. 2) THEN
                  lac_1_ad = y14_ad
                ELSE
                  lac_1_ad = 0.0
                END IF
                pmp_1_ad = 0.0
              END IF
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  lac_1_ad = lac_1_ad + x9_ad
                ELSE
                  pmp_1_ad = pmp_1_ad + x9_ad
                END IF
              ELSE IF (branch .EQ. 2) THEN
                lac_1_ad = lac_1_ad + x9_ad
              END IF
              pmp_1_ad = pmp_1_ad + lac_1_ad
              dq_ad(i, j+1) = dq_ad(i, j+1) + 1.5*lac_1_ad
              dq_ad(i, j) = dq_ad(i, j) - 2.*pmp_1_ad
 110          al_ad(i, j+1) = al_ad(i, j+1) + br_ad(i, j)
              v_ad(i, j) = v_ad(i, j) - bl_ad(i, j) - br_ad(i, j)
              br_ad(i, j) = 0.0
              al_ad(i, j) = al_ad(i, j) + bl_ad(i, j)
              bl_ad(i, j) = 0.0
            END DO
          END DO
        ELSE
          DO j=je3,js3,-1
            DO i=ie+1,is,-1
              al_ad(i, j+1) = al_ad(i, j+1) + br_ad(i, j)
              v_ad(i, j) = v_ad(i, j) - bl_ad(i, j) - br_ad(i, j)
              br_ad(i, j) = 0.0
              al_ad(i, j) = al_ad(i, j) + bl_ad(i, j)
              bl_ad(i, j) = 0.0
            END DO
          END DO
        END IF
        DO j=je3+1,js3,-1
          DO i=ie+1,is,-1
            v_ad(i, j-1) = v_ad(i, j-1) + 0.5*al_ad(i, j)
            v_ad(i, j) = v_ad(i, j) + 0.5*al_ad(i, j)
            dm_ad(i, j-1) = dm_ad(i, j-1) + r3*al_ad(i, j)
            dm_ad(i, j) = dm_ad(i, j) - r3*al_ad(i, j)
            al_ad(i, j) = 0.0
          END DO
        END DO
      END IF
      DO j=je+2,js-3,-1
        DO i=ie+1,is,-1
          v_ad(i, j+1) = v_ad(i, j+1) + dq_ad(i, j)
          v_ad(i, j) = v_ad(i, j) - dq_ad(i, j)
          dq_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+2,js-2,-1
        DO i=ie+1,is,-1
          xt = 0.25*(v(i, j+1)-v(i, j-1))
          min3_ad = SIGN(1.d0, min3*xt)*dm_ad(i, j)
          dm_ad(i, j) = 0.0
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(min3)
              z1_ad = min3_ad
              y3_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min3)
              y3_ad = min3_ad
              z1_ad = 0.0
            END IF
            x4_ad = 0.0
          ELSE
            IF (branch .EQ. 2) THEN
              CALL POPREALARRAY_ADM(min3)
              z1_ad = min3_ad
              x4_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min3)
              x4_ad = min3_ad
              z1_ad = 0.0
            END IF
            y3_ad = 0.0
          END IF
          v_ad(i, j) = v_ad(i, j) + z1_ad
          min6_ad = -z1_ad
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              v_ad(i, j+1) = v_ad(i, j+1) + min6_ad
            ELSE
              v_ad(i, j) = v_ad(i, j) + min6_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            v_ad(i, j+1) = v_ad(i, j+1) + min6_ad
          ELSE
            v_ad(i, j-1) = v_ad(i, j-1) + min6_ad
          END IF
          max1_ad = y3_ad
          v_ad(i, j) = v_ad(i, j) - y3_ad
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              v_ad(i, j+1) = v_ad(i, j+1) + max1_ad
            ELSE
              v_ad(i, j) = v_ad(i, j) + max1_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            v_ad(i, j+1) = v_ad(i, j+1) + max1_ad
          ELSE
            v_ad(i, j-1) = v_ad(i, j-1) + max1_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            xt_ad = x4_ad
          ELSE
            xt_ad = -x4_ad
          END IF
          v_ad(i, j+1) = v_ad(i, j+1) + 0.25*xt_ad
          v_ad(i, j-1) = v_ad(i, j-1) - 0.25*xt_ad
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
  END SUBROUTINE YTP_V_ADM
  SUBROUTINE YTP_V(is, ie, js, je, isd, ied, jsd, jed, c, u, v, flux, &
&   jord, dy, rdy, npx, npy, grid_type, nested)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
!  Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL, INTENT(OUT) :: flux(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rdy(isd:ied+1, jsd:jed)
    INTEGER, INTENT(IN) :: npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
! Local:
    LOGICAL, DIMENSION(is:ie+1, js-1:je+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: dm(is:ie+1, js-2:je+2)
    REAL :: al(is:ie+1, js-1:je+2)
    REAL, DIMENSION(is:ie+1, js-1:je+1) :: bl, br, b0
    REAL :: dq(is:ie+1, js-3:je+2)
    REAL :: xt, dl, dr, pmp, lac, cfl
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: x0, x1, x0r, x0l
    INTEGER :: i, j, is1, ie1, js3, je3
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min2
    REAL :: abs0
    REAL :: min3
    REAL :: min4
    REAL :: min5
    REAL :: abs1
    REAL :: abs2
    REAL :: abs3
    REAL :: abs4
    REAL :: max1
    REAL :: min6
    REAL :: abs5
    REAL :: abs6
    INTEGER :: arg1
    REAL :: x12
    REAL :: x11
    REAL :: x10
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: y17
    REAL :: y16
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
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
    IF (nested .OR. grid_type .GT. 3) THEN
      js3 = js - 1
      je3 = je + 1
    ELSE
      IF (3 .LT. js - 1) THEN
        js3 = js - 1
      ELSE
        js3 = 3
      END IF
      IF (npy - 3 .GT. je + 1) THEN
        je3 = je + 1
      ELSE
        je3 = npy - 3
      END IF
    END IF
    IF (jord .EQ. 1) THEN
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux(i, j) = v(i, j-1)
          ELSE
            flux(i, j) = v(i, j)
          END IF
        END DO
      END DO
    ELSE IF (jord .LT. 8) THEN
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
      DO j=js3,je3+1
        DO i=is,ie+1
          al(i, j) = p1*(v(i, j-1)+v(i, j)) + p2*(v(i, j-2)+v(i, j+1))
        END DO
      END DO
      DO j=js3,je3
        DO i=is,ie+1
          bl(i, j) = al(i, j) - v(i, j)
          br(i, j) = al(i, j+1) - v(i, j)
        END DO
      END DO
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
            bl(i, 0) = c1*v(i, -2) + c2*v(i, -1) + c3*v(i, 0) - v(i, 0)
            xt = 0.5*(((2.*dy(i, 0)+dy(i, -1))*v(i, 0)-dy(i, 0)*v(i, -1)&
&             )/(dy(i, 0)+dy(i, -1))+((2.*dy(i, 1)+dy(i, 2))*v(i, 1)-dy(&
&             i, 1)*v(i, 2))/(dy(i, 1)+dy(i, 2)))
            br(i, 0) = xt - v(i, 0)
            bl(i, 1) = xt - v(i, 1)
            xt = c3*v(i, 1) + c2*v(i, 2) + c1*v(i, 3)
            br(i, 1) = xt - v(i, 1)
            bl(i, 2) = xt - v(i, 2)
            br(i, 2) = al(i, 3) - v(i, 2)
          END DO
          IF (is .EQ. 1) THEN
! out
            bl(1, 0) = 0.
! edge
            br(1, 0) = 0.
! edge
            bl(1, 1) = 0.
! in
            br(1, 1) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
! out
            bl(npx, 0) = 0.
! edge
            br(npx, 0) = 0.
! edge
            bl(npx, 1) = 0.
! in
            br(npx, 1) = 0.
          END IF
        END IF
!      j=2
!      call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
            bl(i, npy-2) = al(i, npy-2) - v(i, npy-2)
            xt = c1*v(i, npy-3) + c2*v(i, npy-2) + c3*v(i, npy-1)
            br(i, npy-2) = xt - v(i, npy-2)
            bl(i, npy-1) = xt - v(i, npy-1)
            xt = 0.5*(((2.*dy(i, npy-1)+dy(i, npy-2))*v(i, npy-1)-dy(i, &
&             npy-1)*v(i, npy-2))/(dy(i, npy-1)+dy(i, npy-2))+((2.*dy(i&
&             , npy)+dy(i, npy+1))*v(i, npy)-dy(i, npy)*v(i, npy+1))/(dy&
&             (i, npy)+dy(i, npy+1)))
            br(i, npy-1) = xt - v(i, npy-1)
            bl(i, npy) = xt - v(i, npy)
            br(i, npy) = c3*v(i, npy) + c2*v(i, npy+1) + c1*v(i, npy+2) &
&             - v(i, npy)
          END DO
          IF (is .EQ. 1) THEN
! in
            bl(1, npy-1) = 0.
! edge
            br(1, npy-1) = 0.
! edge
            bl(1, npy) = 0.
! out
            br(1, npy) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
! in
            bl(npx, npy-1) = 0.
! edge
            br(npx, npy-1) = 0.
! edge
            bl(npx, npy) = 0.
! out
            br(npx, npy) = 0.
          END IF
        END IF
      END IF
!      j=npy-2
!      call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
      DO j=js-1,je+1
        DO i=is,ie+1
          b0(i, j) = bl(i, j) + br(i, j)
        END DO
      END DO
      IF (jord .EQ. 2) THEN
! Perfectly linear
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdy(i, j-1)
              flux(i, j) = v(i, j-1) + (1.-cfl)*(br(i, j-1)-cfl*b0(i, j-&
&               1))
            ELSE
              cfl = c(i, j)*rdy(i, j)
              flux(i, j) = v(i, j) + (1.+cfl)*(bl(i, j)+cfl*b0(i, j))
            END IF
          END DO
        END DO
      ELSE IF (jord .EQ. 3) THEN
        DO j=js-1,je+1
          DO i=is,ie+1
            IF (b0(i, j) .GE. 0.) THEN
              x0 = b0(i, j)
            ELSE
              x0 = -b0(i, j)
            END IF
            IF (bl(i, j) - br(i, j) .GE. 0.) THEN
              x1 = bl(i, j) - br(i, j)
            ELSE
              x1 = -(bl(i, j)-br(i, j))
            END IF
            smt5(i, j) = x0 .LT. x1
            smt6(i, j) = 3.*x0 .LT. x1
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            fx0(i) = 0.
          END DO
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdy(i, j-1)
              IF (smt6(i, j-1) .OR. smt5(i, j)) THEN
                fx0(i) = br(i, j-1) - cfl*b0(i, j-1)
              ELSE IF (smt5(i, j-1)) THEN
                IF (bl(i, j-1) .GE. 0.) THEN
                  x2 = bl(i, j-1)
                ELSE
                  x2 = -bl(i, j-1)
                END IF
                IF (br(i, j-1) .GE. 0.) THEN
                  y1 = br(i, j-1)
                ELSE
                  y1 = -br(i, j-1)
                END IF
                IF (x2 .GT. y1) THEN
                  min1 = y1
                ELSE
                  min1 = x2
                END IF
! piece-wise linear
                fx0(i) = SIGN(min1, br(i, j-1))
              END IF
              flux(i, j) = v(i, j-1) + (1.-cfl)*fx0(i)
            ELSE
              cfl = c(i, j)*rdy(i, j)
              IF (smt6(i, j) .OR. smt5(i, j-1)) THEN
                fx0(i) = bl(i, j) + cfl*b0(i, j)
              ELSE IF (smt5(i, j)) THEN
                IF (bl(i, j) .GE. 0.) THEN
                  x3 = bl(i, j)
                ELSE
                  x3 = -bl(i, j)
                END IF
                IF (br(i, j) .GE. 0.) THEN
                  y2 = br(i, j)
                ELSE
                  y2 = -br(i, j)
                END IF
                IF (x3 .GT. y2) THEN
                  min2 = y2
                ELSE
                  min2 = x3
                END IF
                fx0(i) = SIGN(min2, bl(i, j))
              END IF
              flux(i, j) = v(i, j) + (1.+cfl)*fx0(i)
            END IF
          END DO
        END DO
      ELSE IF (jord .EQ. 4) THEN
        DO j=js-1,je+1
          DO i=is,ie+1
            IF (b0(i, j) .GE. 0.) THEN
              x0 = b0(i, j)
            ELSE
              x0 = -b0(i, j)
            END IF
            IF (bl(i, j) - br(i, j) .GE. 0.) THEN
              x1 = bl(i, j) - br(i, j)
            ELSE
              x1 = -(bl(i, j)-br(i, j))
            END IF
            smt5(i, j) = x0 .LT. x1
            smt6(i, j) = 3.*x0 .LT. x1
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              IF (smt6(i, j-1) .OR. smt5(i, j)) THEN
                cfl = c(i, j)*rdy(i, j-1)
                flux(i, j) = v(i, j-1) + (1.-cfl)*(br(i, j-1)-cfl*b0(i, &
&                 j-1))
              ELSE
                flux(i, j) = v(i, j-1)
              END IF
            ELSE IF (smt6(i, j) .OR. smt5(i, j-1)) THEN
              cfl = c(i, j)*rdy(i, j)
              flux(i, j) = v(i, j) + (1.+cfl)*(bl(i, j)+cfl*b0(i, j))
            ELSE
              flux(i, j) = v(i, j)
            END IF
          END DO
        END DO
      ELSE
! jord = 5,6,7
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7
        IF (jord .EQ. 5) THEN
          DO j=js-1,je+1
            DO i=is,ie+1
              smt5(i, j) = bl(i, j)*br(i, j) .LT. 0.
            END DO
          END DO
        ELSE
! ord = 6, 7
          DO j=js-1,je+1
            DO i=is,ie+1
              IF (3.*b0(i, j) .GE. 0.) THEN
                abs0 = 3.*b0(i, j)
              ELSE
                abs0 = -(3.*b0(i, j))
              END IF
              IF (bl(i, j) - br(i, j) .GE. 0.) THEN
                abs4 = bl(i, j) - br(i, j)
              ELSE
                abs4 = -(bl(i, j)-br(i, j))
              END IF
              smt5(i, j) = abs0 .LT. abs4
            END DO
          END DO
        END IF
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdy(i, j-1)
              fx0(i) = (1.-cfl)*(br(i, j-1)-cfl*b0(i, j-1))
              flux(i, j) = v(i, j-1)
            ELSE
              cfl = c(i, j)*rdy(i, j)
              fx0(i) = (1.+cfl)*(bl(i, j)+cfl*b0(i, j))
              flux(i, j) = v(i, j)
            END IF
            IF (smt5(i, j-1) .OR. smt5(i, j)) flux(i, j) = flux(i, j) + &
&               fx0(i)
          END DO
        END DO
      END IF
    ELSE
! jord= 8, 9, 10
      DO j=js-2,je+2
        DO i=is,ie+1
          xt = 0.25*(v(i, j+1)-v(i, j-1))
          IF (xt .GE. 0.) THEN
            x4 = xt
          ELSE
            x4 = -xt
          END IF
          IF (v(i, j-1) .LT. v(i, j)) THEN
            IF (v(i, j) .LT. v(i, j+1)) THEN
              max1 = v(i, j+1)
            ELSE
              max1 = v(i, j)
            END IF
          ELSE IF (v(i, j-1) .LT. v(i, j+1)) THEN
            max1 = v(i, j+1)
          ELSE
            max1 = v(i, j-1)
          END IF
          y3 = max1 - v(i, j)
          IF (v(i, j-1) .GT. v(i, j)) THEN
            IF (v(i, j) .GT. v(i, j+1)) THEN
              min6 = v(i, j+1)
            ELSE
              min6 = v(i, j)
            END IF
          ELSE IF (v(i, j-1) .GT. v(i, j+1)) THEN
            min6 = v(i, j+1)
          ELSE
            min6 = v(i, j-1)
          END IF
          z1 = v(i, j) - min6
          IF (x4 .GT. y3) THEN
            IF (y3 .GT. z1) THEN
              min3 = z1
            ELSE
              min3 = y3
            END IF
          ELSE IF (x4 .GT. z1) THEN
            min3 = z1
          ELSE
            min3 = x4
          END IF
          dm(i, j) = SIGN(min3, xt)
        END DO
      END DO
      DO j=js-3,je+2
        DO i=is,ie+1
          dq(i, j) = v(i, j+1) - v(i, j)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        DO j=js3,je3+1
          DO i=is,ie+1
            al(i, j) = 0.5*(v(i, j-1)+v(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        IF (jord .EQ. 8) THEN
          DO j=js3,je3
            DO i=is,ie+1
              xt = 2.*dm(i, j)
              IF (xt .GE. 0.) THEN
                x5 = xt
              ELSE
                x5 = -xt
              END IF
              IF (al(i, j) - v(i, j) .GE. 0.) THEN
                y4 = al(i, j) - v(i, j)
              ELSE
                y4 = -(al(i, j)-v(i, j))
              END IF
              IF (x5 .GT. y4) THEN
                min4 = y4
              ELSE
                min4 = x5
              END IF
              bl(i, j) = -SIGN(min4, xt)
              IF (xt .GE. 0.) THEN
                x6 = xt
              ELSE
                x6 = -xt
              END IF
              IF (al(i, j+1) - v(i, j) .GE. 0.) THEN
                y5 = al(i, j+1) - v(i, j)
              ELSE
                y5 = -(al(i, j+1)-v(i, j))
              END IF
              IF (x6 .GT. y5) THEN
                min5 = y5
              ELSE
                min5 = x6
              END IF
              br(i, j) = SIGN(min5, xt)
            END DO
          END DO
        ELSE IF (jord .EQ. 9) THEN
          DO j=js3,je3
            DO i=is,ie+1
              pmp_1 = -(2.*dq(i, j))
              lac_1 = pmp_1 + 1.5*dq(i, j+1)
              IF (0. .LT. pmp_1) THEN
                IF (pmp_1 .LT. lac_1) THEN
                  x7 = lac_1
                ELSE
                  x7 = pmp_1
                END IF
              ELSE IF (0. .LT. lac_1) THEN
                x7 = lac_1
              ELSE
                x7 = 0.
              END IF
              IF (0. .GT. pmp_1) THEN
                IF (pmp_1 .GT. lac_1) THEN
                  y12 = lac_1
                ELSE
                  y12 = pmp_1
                END IF
              ELSE IF (0. .GT. lac_1) THEN
                y12 = lac_1
              ELSE
                y12 = 0.
              END IF
              IF (al(i, j) - v(i, j) .LT. y12) THEN
                y6 = y12
              ELSE
                y6 = al(i, j) - v(i, j)
              END IF
              IF (x7 .GT. y6) THEN
                bl(i, j) = y6
              ELSE
                bl(i, j) = x7
              END IF
              pmp_2 = 2.*dq(i, j-1)
              lac_2 = pmp_2 - 1.5*dq(i, j-2)
              IF (0. .LT. pmp_2) THEN
                IF (pmp_2 .LT. lac_2) THEN
                  x8 = lac_2
                ELSE
                  x8 = pmp_2
                END IF
              ELSE IF (0. .LT. lac_2) THEN
                x8 = lac_2
              ELSE
                x8 = 0.
              END IF
              IF (0. .GT. pmp_2) THEN
                IF (pmp_2 .GT. lac_2) THEN
                  y13 = lac_2
                ELSE
                  y13 = pmp_2
                END IF
              ELSE IF (0. .GT. lac_2) THEN
                y13 = lac_2
              ELSE
                y13 = 0.
              END IF
              IF (al(i, j+1) - v(i, j) .LT. y13) THEN
                y7 = y13
              ELSE
                y7 = al(i, j+1) - v(i, j)
              END IF
              IF (x8 .GT. y7) THEN
                br(i, j) = y7
              ELSE
                br(i, j) = x8
              END IF
            END DO
          END DO
        ELSE IF (jord .EQ. 10) THEN
          DO j=js3,je3
            DO i=is,ie+1
              bl(i, j) = al(i, j) - v(i, j)
              br(i, j) = al(i, j+1) - v(i, j)
              IF (dm(i, j) .GE. 0.) THEN
                abs1 = dm(i, j)
              ELSE
                abs1 = -dm(i, j)
              END IF
!           if ( abs(dm(i,j-1))+abs(dm(i,j))+abs(dm(i,j+1)) < near_zero ) then
              IF (abs1 .LT. near_zero) THEN
                IF (dm(i, j-1) .GE. 0.) THEN
                  abs2 = dm(i, j-1)
                ELSE
                  abs2 = -dm(i, j-1)
                END IF
                IF (dm(i, j+1) .GE. 0.) THEN
                  abs5 = dm(i, j+1)
                ELSE
                  abs5 = -dm(i, j+1)
                END IF
                IF (abs2 + abs5 .LT. near_zero) THEN
                  bl(i, j) = 0.
                  br(i, j) = 0.
                END IF
              ELSE
                IF (3.*(bl(i, j)+br(i, j)) .GE. 0.) THEN
                  abs3 = 3.*(bl(i, j)+br(i, j))
                ELSE
                  abs3 = -(3.*(bl(i, j)+br(i, j)))
                END IF
                IF (bl(i, j) - br(i, j) .GE. 0.) THEN
                  abs6 = bl(i, j) - br(i, j)
                ELSE
                  abs6 = -(bl(i, j)-br(i, j))
                END IF
                IF (abs3 .GT. abs6) THEN
                  pmp_1 = -(2.*dq(i, j))
                  lac_1 = pmp_1 + 1.5*dq(i, j+1)
                  IF (0. .LT. pmp_1) THEN
                    IF (pmp_1 .LT. lac_1) THEN
                      x9 = lac_1
                    ELSE
                      x9 = pmp_1
                    END IF
                  ELSE IF (0. .LT. lac_1) THEN
                    x9 = lac_1
                  ELSE
                    x9 = 0.
                  END IF
                  IF (0. .GT. pmp_1) THEN
                    IF (pmp_1 .GT. lac_1) THEN
                      y14 = lac_1
                    ELSE
                      y14 = pmp_1
                    END IF
                  ELSE IF (0. .GT. lac_1) THEN
                    y14 = lac_1
                  ELSE
                    y14 = 0.
                  END IF
                  IF (bl(i, j) .LT. y14) THEN
                    y8 = y14
                  ELSE
                    y8 = bl(i, j)
                  END IF
                  IF (x9 .GT. y8) THEN
                    bl(i, j) = y8
                  ELSE
                    bl(i, j) = x9
                  END IF
                  pmp_2 = 2.*dq(i, j-1)
                  lac_2 = pmp_2 - 1.5*dq(i, j-2)
                  IF (0. .LT. pmp_2) THEN
                    IF (pmp_2 .LT. lac_2) THEN
                      x10 = lac_2
                    ELSE
                      x10 = pmp_2
                    END IF
                  ELSE IF (0. .LT. lac_2) THEN
                    x10 = lac_2
                  ELSE
                    x10 = 0.
                  END IF
                  IF (0. .GT. pmp_2) THEN
                    IF (pmp_2 .GT. lac_2) THEN
                      y15 = lac_2
                    ELSE
                      y15 = pmp_2
                    END IF
                  ELSE IF (0. .GT. lac_2) THEN
                    y15 = lac_2
                  ELSE
                    y15 = 0.
                  END IF
                  IF (br(i, j) .LT. y15) THEN
                    y9 = y15
                  ELSE
                    y9 = br(i, j)
                  END IF
                  IF (x10 .GT. y9) THEN
                    br(i, j) = y9
                  ELSE
                    br(i, j) = x10
                  END IF
                END IF
              END IF
            END DO
          END DO
        ELSE
! Unlimited:
          DO j=js3,je3
            DO i=is,ie+1
              bl(i, j) = al(i, j) - v(i, j)
              br(i, j) = al(i, j+1) - v(i, j)
            END DO
          END DO
        END IF
!--------------
! fix the edges
!--------------
        IF (js .EQ. 1 .AND. (.NOT.nested)) THEN
          DO i=is,ie+1
            br(i, 2) = al(i, 3) - v(i, 2)
            xt = s15*v(i, 1) + s11*v(i, 2) - s14*dm(i, 2)
            br(i, 1) = xt - v(i, 1)
            bl(i, 2) = xt - v(i, 2)
            bl(i, 0) = s14*dm(i, -1) - s11*dq(i, -1)
            x0l = 0.5*((2.*dy(i, 0)+dy(i, -1))*v(i, 0)-dy(i, 0)*v(i, -1)&
&             )/(dy(i, 0)+dy(i, -1))
            x0r = 0.5*((2.*dy(i, 1)+dy(i, 2))*v(i, 1)-dy(i, 1)*v(i, 2))/&
&             (dy(i, 1)+dy(i, 2))
            xt = x0l + x0r
            bl(i, 1) = xt - v(i, 1)
            br(i, 0) = xt - v(i, 0)
          END DO
          IF (is .EQ. 1) THEN
! out
            bl(1, 0) = 0.
! edge
            br(1, 0) = 0.
! edge
            bl(1, 1) = 0.
! in
            br(1, 1) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
! out
            bl(npx, 0) = 0.
! edge
            br(npx, 0) = 0.
! edge
            bl(npx, 1) = 0.
! in
            br(npx, 1) = 0.
          END IF
          j = 2
          arg1 = ie - is + 2
          CALL PERT_PPM(arg1, v(is:ie+1, j), bl(is:ie+1, j), br(is:ie+1&
&                 , j), -1)
        END IF
        IF (je + 1 .EQ. npy .AND. (.NOT.nested)) THEN
          DO i=is,ie+1
            bl(i, npy-2) = al(i, npy-2) - v(i, npy-2)
            xt = s15*v(i, npy-1) + s11*v(i, npy-2) + s14*dm(i, npy-2)
            br(i, npy-2) = xt - v(i, npy-2)
            bl(i, npy-1) = xt - v(i, npy-1)
            br(i, npy) = s11*dq(i, npy) - s14*dm(i, npy+1)
            x0l = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*v(i, npy-1)-dy(i, &
&             npy-1)*v(i, npy-2))/(dy(i, npy-1)+dy(i, npy-2))
            x0r = 0.5*((2.*dy(i, npy)+dy(i, npy+1))*v(i, npy)-dy(i, npy)&
&             *v(i, npy+1))/(dy(i, npy)+dy(i, npy+1))
            xt = x0l + x0r
            br(i, npy-1) = xt - v(i, npy-1)
            bl(i, npy) = xt - v(i, npy)
          END DO
          IF (is .EQ. 1) THEN
! in
            bl(1, npy-1) = 0.
! edge
            br(1, npy-1) = 0.
! edge
            bl(1, npy) = 0.
! out
            br(1, npy) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
! in
            bl(npx, npy-1) = 0.
! edge
            br(npx, npy-1) = 0.
! edge
            bl(npx, npy) = 0.
! out
            br(npx, npy) = 0.
          END IF
          j = npy - 2
          arg1 = ie - is + 2
          CALL PERT_PPM(arg1, v(is:ie+1, j), bl(is:ie+1, j), br(is:ie+1&
&                 , j), -1)
        END IF
      ELSE
        DO j=js-1,je+2
          DO i=is,ie+1
            al(i, j) = 0.5*(v(i, j-1)+v(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        DO j=js-1,je+1
          DO i=is,ie+1
            pmp = 2.*dq(i, j-1)
            lac = pmp - 1.5*dq(i, j-2)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x11 = lac
              ELSE
                x11 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x11 = lac
            ELSE
              x11 = 0.
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y16 = lac
              ELSE
                y16 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y16 = lac
            ELSE
              y16 = 0.
            END IF
            IF (al(i, j+1) - v(i, j) .LT. y16) THEN
              y10 = y16
            ELSE
              y10 = al(i, j+1) - v(i, j)
            END IF
            IF (x11 .GT. y10) THEN
              br(i, j) = y10
            ELSE
              br(i, j) = x11
            END IF
            pmp = -(2.*dq(i, j))
            lac = pmp + 1.5*dq(i, j+1)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x12 = lac
              ELSE
                x12 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x12 = lac
            ELSE
              x12 = 0.
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y17 = lac
              ELSE
                y17 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y17 = lac
            ELSE
              y17 = 0.
            END IF
            IF (al(i, j) - v(i, j) .LT. y17) THEN
              y11 = y17
            ELSE
              y11 = al(i, j) - v(i, j)
            END IF
            IF (x12 .GT. y11) THEN
              bl(i, j) = y11
            ELSE
              bl(i, j) = x12
            END IF
          END DO
        END DO
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            cfl = c(i, j)*rdy(i, j-1)
            flux(i, j) = v(i, j-1) + (1.-cfl)*(br(i, j-1)-cfl*(bl(i, j-1&
&             )+br(i, j-1)))
          ELSE
            cfl = c(i, j)*rdy(i, j)
            flux(i, j) = v(i, j) + (1.+cfl)*(bl(i, j)+cfl*(bl(i, j)+br(i&
&             , j)))
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE YTP_V
!  Differentiation of d2a2c_vect in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_m
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
!   gradient     of useful results: u v ua uc ut va vc vt
!   with respect to varying inputs: u v ua uc ut va vc vt
!There is a limit to how far this routine can fill uc and vc in the
! halo, and so either mpp_update_domains or some sort of boundary
!  routine (extrapolation, outflow, interpolation from a nested grid)
!   is needed after c_sw is completed if these variables are needed
!    in the halo
  SUBROUTINE D2A2C_VECT_FWD(u, v, ua, va, uc, vc, ut, vt, dord4, &
&   gridstruct, bd, npx, npy, nested, grid_type)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    LOGICAL, INTENT(IN) :: dord4
    REAL, INTENT(IN) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL, INTENT(IN) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed) :: uc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1) :: vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ua, va, ut, vt
    INTEGER, INTENT(IN) :: npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! Local
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: utmp, vtmp
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsin2
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
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

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    cosa_s => gridstruct%cosa_s
    rsin_u => gridstruct%rsin_u
    rsin_v => gridstruct%rsin_v
    rsin2 => gridstruct%rsin2
    dxa => gridstruct%dxa
    dya => gridstruct%dya
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
! Initialize the non-existing corner regions
    utmp(:, :) = big_number
    vtmp(:, :) = big_number
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
          CALL PUSHREALARRAY(ua(i, j))
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          CALL PUSHREALARRAY(va(i, j))
          va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
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
        ad_from = max2
        DO i=ad_from,min2
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
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
        DO i=ad_from0,min4
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
        CALL PUSHINTEGER(i - 1)
        CALL PUSHINTEGER(ad_from0)
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
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
          DO j=npy-npt+1,jed
            DO i=isd,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
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
          DO j=max5,min5
            DO i=isd,npt-1
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
          IF (npt .LT. jsd) THEN
            CALL PUSHCONTROL(1,1)
            max6 = jsd
          ELSE
            max6 = npt
            CALL PUSHCONTROL(1,0)
          END IF
          IF (npy - npt .GT. jed) THEN
            CALL PUSHCONTROL(1,1)
            min6 = jed
          ELSE
            min6 = npy - npt
            CALL PUSHCONTROL(1,0)
          END IF
          DO j=max6,min6
            DO i=npx-npt+1,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
          CALL PUSHCONTROL(2,2)
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,0)
      END IF
! Contra-variant components at cell center:
      DO j=js-1-id,je+1+id
        DO i=is-1-id,ie+1+id
          CALL PUSHREALARRAY(ua(i, j))
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          CALL PUSHREALARRAY(va(i, j))
          va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    END IF
! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
    IF (gridstruct%sw_corner) THEN
      DO i=-2,0
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%se_corner) THEN
      DO i=0,2
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%ne_corner) THEN
      DO i=0,2
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%nw_corner) THEN
      DO i=-2,0
        utmp(i, npy) = vtmp(0, je+i)
      END DO
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
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
    DO j=js-1,je+1
      DO i=ifirst,ilast
        uc(i, j) = a2*(utmp(i-2, j)+utmp(i+1, j)) + a1*(utmp(i-1, j)+&
&         utmp(i, j))
        CALL PUSHREALARRAY(ut(i, j))
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
! Xdir:
      IF (gridstruct%sw_corner) THEN
        CALL PUSHREALARRAY(ua(-1, 0))
        ua(-1, 0) = -va(0, 2)
        CALL PUSHREALARRAY(ua(0, 0))
        ua(0, 0) = -va(0, 1)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%se_corner) THEN
        CALL PUSHREALARRAY(ua(npx, 0))
        ua(npx, 0) = va(npx, 1)
        CALL PUSHREALARRAY(ua(npx+1, 0))
        ua(npx+1, 0) = va(npx, 2)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%ne_corner) THEN
        CALL PUSHREALARRAY(ua(npx, npy))
        ua(npx, npy) = -va(npx, npy-1)
        CALL PUSHREALARRAY(ua(npx+1, npy))
        ua(npx+1, npy) = -va(npx, npy-2)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (gridstruct%nw_corner) THEN
        CALL PUSHREALARRAY(ua(-1, npy))
        ua(-1, npy) = va(0, npy-2)
        CALL PUSHREALARRAY(ua(0, npy))
        ua(0, npy) = va(0, npy-1)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (is .EQ. 1 .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc(0, j) = c1*utmp(-2, j) + c2*utmp(-1, j) + c3*utmp(0, j)
          CALL PUSHREALARRAY(ut(1, j))
          ut(1, j) = EDGE_INTERPOLATE4_FWD(ua(-1:2, j), dxa(-1:2, j))
!Want to use the UPSTREAM value
          IF (ut(1, j) .GT. 0.) THEN
            uc(1, j) = ut(1, j)*sin_sg(0, j, 3)
            CALL PUSHCONTROL(1,0)
          ELSE
            uc(1, j) = ut(1, j)*sin_sg(1, j, 1)
            CALL PUSHCONTROL(1,1)
          END IF
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
          CALL PUSHREALARRAY(ut(0, j))
          ut(0, j) = (uc(0, j)-v(0, j)*cosa_u(0, j))*rsin_u(0, j)
          CALL PUSHREALARRAY(ut(2, j))
          ut(2, j) = (uc(2, j)-v(2, j)*cosa_u(2, j))*rsin_u(2, j)
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ie + 1 .EQ. npx .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
          CALL PUSHREALARRAY(ut(npx, j))
          ut(npx, j) = EDGE_INTERPOLATE4_FWD(ua(npx-2:npx+1, j), dxa(npx&
&           -2:npx+1, j))
          IF (ut(npx, j) .GT. 0.) THEN
            uc(npx, j) = ut(npx, j)*sin_sg(npx-1, j, 3)
            CALL PUSHCONTROL(1,0)
          ELSE
            uc(npx, j) = ut(npx, j)*sin_sg(npx, j, 1)
            CALL PUSHCONTROL(1,1)
          END IF
          uc(npx+1, j) = c3*utmp(npx, j) + c2*utmp(npx+1, j) + c1*utmp(&
&           npx+2, j)
          CALL PUSHREALARRAY(ut(npx-1, j))
          ut(npx-1, j) = (uc(npx-1, j)-v(npx-1, j)*cosa_u(npx-1, j))*&
&           rsin_u(npx-1, j)
          CALL PUSHREALARRAY(ut(npx+1, j))
          ut(npx+1, j) = (uc(npx+1, j)-v(npx+1, j)*cosa_u(npx+1, j))*&
&           rsin_u(npx+1, j)
        END DO
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
    IF (gridstruct%sw_corner) THEN
      DO j=-2,0
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%nw_corner) THEN
      DO j=0,2
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%se_corner) THEN
      DO j=-2,0
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%ne_corner) THEN
      DO j=0,2
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%sw_corner) THEN
      CALL PUSHREALARRAY(va(0, -1))
      va(0, -1) = -ua(2, 0)
      CALL PUSHREALARRAY(va(0, 0))
      va(0, 0) = -ua(1, 0)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%se_corner) THEN
      CALL PUSHREALARRAY(va(npx, 0))
      va(npx, 0) = ua(npx-1, 0)
      CALL PUSHREALARRAY(va(npx, -1))
      va(npx, -1) = ua(npx-2, 0)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%ne_corner) THEN
      CALL PUSHREALARRAY(va(npx, npy))
      va(npx, npy) = -ua(npx-1, npy)
      CALL PUSHREALARRAY(va(npx, npy+1))
      va(npx, npy+1) = -ua(npx-2, npy)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (gridstruct%nw_corner) THEN
      CALL PUSHREALARRAY(va(0, npy))
      va(0, npy) = ua(1, npy)
      CALL PUSHREALARRAY(va(0, npy+1))
      va(0, npy+1) = ua(2, npy)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1 .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            CALL PUSHREALARRAY(vt(i, j))
            vt(i, j) = EDGE_INTERPOLATE4_FWD(va(i, -1:2), dya(i, -1:2))
            IF (vt(i, j) .GT. 0.) THEN
              vc(i, j) = vt(i, j)*sin_sg(i, j-1, 4)
              CALL PUSHCONTROL(1,1)
            ELSE
              vc(i, j) = vt(i, j)*sin_sg(i, j, 2)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          CALL PUSHCONTROL(3,4)
        ELSE IF (j .EQ. 0 .OR. (j .EQ. npy - 1 .AND. (.NOT.nested))) &
&       THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
            CALL PUSHREALARRAY(vt(i, j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
          CALL PUSHCONTROL(3,3)
        ELSE IF (j .EQ. 2 .OR. (j .EQ. npy + 1 .AND. (.NOT.nested))) &
&       THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j+1) + c2*vtmp(i, j) + c3*vtmp(i, j-1)
            CALL PUSHREALARRAY(vt(i, j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
          CALL PUSHCONTROL(3,2)
        ELSE IF (j .EQ. npy .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            CALL PUSHREALARRAY(vt(i, j))
            vt(i, j) = EDGE_INTERPOLATE4_FWD(va(i, j-2:j+1), dya(i, j-2:&
&             j+1))
            IF (vt(i, j) .GT. 0.) THEN
              vc(i, j) = vt(i, j)*sin_sg(i, j-1, 4)
              CALL PUSHCONTROL(1,1)
            ELSE
              vc(i, j) = vt(i, j)*sin_sg(i, j, 2)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          CALL PUSHCONTROL(3,1)
        ELSE
! 4th order interpolation for interior points:
          DO i=is-1,ie+1
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
            CALL PUSHREALARRAY(vt(i, j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
          CALL PUSHCONTROL(3,0)
        END IF
      END DO
      CALL PUSHINTEGER(npt)
      CALL PUSHINTEGER(jed)
      CALL PUSHINTEGER(ifirst)
      CALL PUSHINTEGER(min6)
      CALL PUSHINTEGER(je)
      CALL PUSHINTEGER(min5)
      CALL PUSHINTEGER(min3)
      CALL PUSHINTEGER(min1)
      CALL PUSHINTEGER(is)
      CALL PUSHINTEGER(isd)
      !CALL PUSHPOINTER8(C_LOC(rsin_v))
      CALL PUSHINTEGER(ie)
      CALL PUSHINTEGER(id)
      !CALL PUSHPOINTER8(C_LOC(dya))
      !CALL PUSHPOINTER8(C_LOC(sin_sg))
      CALL PUSHINTEGER(ied)
      CALL PUSHINTEGER(ilast)
      CALL PUSHINTEGER(jsd)
      !CALL PUSHPOINTER8(C_LOC(cosa_v))
      !CALL PUSHPOINTER8(C_LOC(dxa))
      CALL PUSHINTEGER(max6)
      CALL PUSHINTEGER(max5)
      CALL PUSHINTEGER(max3)
      CALL PUSHINTEGER(max1)
      CALL PUSHINTEGER(js)
      CALL PUSHCONTROL(1,0)
    ELSE
! 4th order interpolation:
      DO j=js-1,je+2
        DO i=is-1,ie+1
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
          CALL PUSHREALARRAY(vt(i, j))
          vt(i, j) = vc(i, j)
        END DO
      END DO
      CALL PUSHINTEGER(npt)
      CALL PUSHINTEGER(jed)
      CALL PUSHINTEGER(ifirst)
      CALL PUSHINTEGER(min6)
      CALL PUSHINTEGER(je)
      CALL PUSHINTEGER(min5)
      CALL PUSHINTEGER(min3)
      CALL PUSHINTEGER(min1)
      CALL PUSHINTEGER(is)
      CALL PUSHINTEGER(isd)
      CALL PUSHINTEGER(ie)
      CALL PUSHINTEGER(id)
      !CALL PUSHPOINTER8(C_LOC(sin_sg))
      CALL PUSHINTEGER(ied)
      CALL PUSHINTEGER(ilast)
      CALL PUSHINTEGER(jsd)
      !CALL PUSHPOINTER8(C_LOC(dxa))
      CALL PUSHINTEGER(max6)
      CALL PUSHINTEGER(max5)
      CALL PUSHINTEGER(max3)
      CALL PUSHINTEGER(max1)
      CALL PUSHINTEGER(js)
      CALL PUSHCONTROL(1,1)
    END IF
  END SUBROUTINE D2A2C_VECT_FWD
!  Differentiation of d2a2c_vect in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: u v ua uc ut va vc vt
!   with respect to varying inputs: u v ua uc ut va vc vt
!There is a limit to how far this routine can fill uc and vc in the
! halo, and so either mpp_update_domains or some sort of boundary
!  routine (extrapolation, outflow, interpolation from a nested grid)
!   is needed after c_sw is completed if these variables are needed
!    in the halo
  SUBROUTINE D2A2C_VECT_BWD(u, u_ad, v, v_ad, ua, ua_ad, va, va_ad, uc, &
&   uc_ad, vc, vc_ad, ut, ut_ad, vt, vt_ad, dord4, gridstruct, bd, npx, &
&   npy, nested, grid_type)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    LOGICAL, INTENT(IN) :: dord4
    REAL, INTENT(IN) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: u_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL, INTENT(IN) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: v_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed) :: uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed) :: uc_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1) :: vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1) :: vc_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ua, va, ut, vt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ua_ad, va_ad, ut_ad&
&   , vt_ad
    INTEGER, INTENT(IN) :: npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: utmp, vtmp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: utmp_ad, vtmp_ad
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsin2
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
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
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: branch
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_d2a2c_vect

    utmp = 0.0
    vtmp = 0.0
    npt = 0
    ifirst = 0
    ilast = 0
    id = 0
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
    ad_to = 0
    ad_to0 = 0
    branch = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    cosa_s => gridstruct%cosa_s
    rsin_u => gridstruct%rsin_u
    rsin_v => gridstruct%rsin_v
    rsin2 => gridstruct%rsin2
    dxa => gridstruct%dxa
    dya => gridstruct%dya
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
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(js)
      CALL POPINTEGER(max1)
      CALL POPINTEGER(max3)
      CALL POPINTEGER(max5)
      CALL POPINTEGER(max6)
      !CALL POPPOINTER8(cptr)
      dxa => gridstruct%dxa ! (/unknown_shape_in_d2a2c_vect/))
      !CALL POPPOINTER8(cptr)
      cosa_v => gridstruct%cosa_v ! (/unknown_shape_in_d2a2c_vect/))
      CALL POPINTEGER(jsd)
      CALL POPINTEGER(ilast)
      CALL POPINTEGER(ied)
      !CALL POPPOINTER8(cptr)
      sin_sg => gridstruct%sin_sg ! (/unknown_shape_in_d2a2c_vect/))
      !CALL POPPOINTER8(cptr)
      dya => gridstruct%dya ! (/unknown_shape_in_d2a2c_vect/))
      CALL POPINTEGER(id)
      CALL POPINTEGER(ie)
      !CALL POPPOINTER8(cptr)
      rsin_v => gridstruct%rsin_v ! (/unknown_shape_in_d2a2c_vect/))
      CALL POPINTEGER(isd)
      CALL POPINTEGER(is)
      CALL POPINTEGER(min1)
      CALL POPINTEGER(min3)
      CALL POPINTEGER(min5)
      CALL POPINTEGER(je)
      CALL POPINTEGER(min6)
      CALL POPINTEGER(ifirst)
      CALL POPINTEGER(jed)
      CALL POPINTEGER(npt)
      vtmp_ad = 0.0
      DO j=je+2,js-1,-1
        CALL POPCONTROL(3,branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            DO i=ie+1,is-1,-1
              CALL POPREALARRAY(vt(i, j))
              temp_ad10 = rsin_v(i, j)*vt_ad(i, j)
              vc_ad(i, j) = vc_ad(i, j) + temp_ad10
              u_ad(i, j) = u_ad(i, j) - cosa_v(i, j)*temp_ad10
              vt_ad(i, j) = 0.0
              vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + a2*vc_ad(i, j)
              vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + a2*vc_ad(i, j)
              vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + a1*vc_ad(i, j)
              vtmp_ad(i, j) = vtmp_ad(i, j) + a1*vc_ad(i, j)
              vc_ad(i, j) = 0.0
            END DO
          ELSE
            DO i=ie+1,is-1,-1
              CALL POPCONTROL(1,branch)
              IF (branch .EQ. 0) THEN
                vt_ad(i, j) = vt_ad(i, j) + sin_sg(i, j, 2)*vc_ad(i, j)
                vc_ad(i, j) = 0.0
              ELSE
                vt_ad(i, j) = vt_ad(i, j) + sin_sg(i, j-1, 4)*vc_ad(i, j&
&                 )
                vc_ad(i, j) = 0.0
              END IF
              CALL EDGE_INTERPOLATE4_BWD(va(i, j-2:j+1), va_ad(i, j-2:j+&
&                                  1), dya(i, j-2:j+1), vt_ad(i, j))
              vt_ad(i, j) = 0.0
              CALL POPREALARRAY(vt(i, j))
            END DO
          END IF
        ELSE IF (branch .EQ. 2) THEN
          DO i=ie+1,is-1,-1
            CALL POPREALARRAY(vt(i, j))
            temp_ad9 = rsin_v(i, j)*vt_ad(i, j)
            vc_ad(i, j) = vc_ad(i, j) + temp_ad9
            u_ad(i, j) = u_ad(i, j) - cosa_v(i, j)*temp_ad9
            vt_ad(i, j) = 0.0
            vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + c1*vc_ad(i, j)
            vtmp_ad(i, j) = vtmp_ad(i, j) + c2*vc_ad(i, j)
            vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + c3*vc_ad(i, j)
            vc_ad(i, j) = 0.0
          END DO
        ELSE IF (branch .EQ. 3) THEN
          DO i=ie+1,is-1,-1
            CALL POPREALARRAY(vt(i, j))
            temp_ad8 = rsin_v(i, j)*vt_ad(i, j)
            vc_ad(i, j) = vc_ad(i, j) + temp_ad8
            u_ad(i, j) = u_ad(i, j) - cosa_v(i, j)*temp_ad8
            vt_ad(i, j) = 0.0
            vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + c1*vc_ad(i, j)
            vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + c2*vc_ad(i, j)
            vtmp_ad(i, j) = vtmp_ad(i, j) + c3*vc_ad(i, j)
            vc_ad(i, j) = 0.0
          END DO
        ELSE
          DO i=ie+1,is-1,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              vt_ad(i, j) = vt_ad(i, j) + sin_sg(i, j, 2)*vc_ad(i, j)
              vc_ad(i, j) = 0.0
            ELSE
              vt_ad(i, j) = vt_ad(i, j) + sin_sg(i, j-1, 4)*vc_ad(i, j)
              vc_ad(i, j) = 0.0
            END IF
            CALL EDGE_INTERPOLATE4_BWD(va(i, -1:2), va_ad(i, -1:2), dya(&
&                                i, -1:2), vt_ad(i, j))
            vt_ad(i, j) = 0.0
            CALL POPREALARRAY(vt(i, j))
          END DO
        END IF
      END DO
    ELSE
      CALL POPINTEGER(js)
      CALL POPINTEGER(max1)
      CALL POPINTEGER(max3)
      CALL POPINTEGER(max5)
      CALL POPINTEGER(max6)
      !CALL POPPOINTER8(cptr)
      dxa => gridstruct%dxa ! (/unknown_shape_in_d2a2c_vect/))
      CALL POPINTEGER(jsd)
      CALL POPINTEGER(ilast)
      CALL POPINTEGER(ied)
      !CALL POPPOINTER8(cptr)
      sin_sg => gridstruct%sin_sg ! (/unknown_shape_in_d2a2c_vect/))
      CALL POPINTEGER(id)
      CALL POPINTEGER(ie)
      CALL POPINTEGER(isd)
      CALL POPINTEGER(is)
      CALL POPINTEGER(min1)
      CALL POPINTEGER(min3)
      CALL POPINTEGER(min5)
      CALL POPINTEGER(je)
      CALL POPINTEGER(min6)
      CALL POPINTEGER(ifirst)
      CALL POPINTEGER(jed)
      CALL POPINTEGER(npt)
      vtmp_ad = 0.0
      DO j=je+2,js-1,-1
        DO i=ie+1,is-1,-1
          CALL POPREALARRAY(vt(i, j))
          vc_ad(i, j) = vc_ad(i, j) + vt_ad(i, j)
          vt_ad(i, j) = 0.0
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
      CALL POPREALARRAY(va(0, npy+1))
      ua_ad(2, npy) = ua_ad(2, npy) + va_ad(0, npy+1)
      va_ad(0, npy+1) = 0.0
      CALL POPREALARRAY(va(0, npy))
      ua_ad(1, npy) = ua_ad(1, npy) + va_ad(0, npy)
      va_ad(0, npy) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(va(npx, npy+1))
      ua_ad(npx-2, npy) = ua_ad(npx-2, npy) - va_ad(npx, npy+1)
      va_ad(npx, npy+1) = 0.0
      CALL POPREALARRAY(va(npx, npy))
      ua_ad(npx-1, npy) = ua_ad(npx-1, npy) - va_ad(npx, npy)
      va_ad(npx, npy) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(va(npx, -1))
      ua_ad(npx-2, 0) = ua_ad(npx-2, 0) + va_ad(npx, -1)
      va_ad(npx, -1) = 0.0
      CALL POPREALARRAY(va(npx, 0))
      ua_ad(npx-1, 0) = ua_ad(npx-1, 0) + va_ad(npx, 0)
      va_ad(npx, 0) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(va(0, 0))
      ua_ad(1, 0) = ua_ad(1, 0) - va_ad(0, 0)
      va_ad(0, 0) = 0.0
      CALL POPREALARRAY(va(0, -1))
      ua_ad(2, 0) = ua_ad(2, 0) - va_ad(0, -1)
      va_ad(0, -1) = 0.0
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
    cosa_u => gridstruct%cosa_u
    rsin_u => gridstruct%rsin_u
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      DO j=je+1,js-1,-1
        CALL POPREALARRAY(ut(npx+1, j))
        temp_ad6 = rsin_u(npx+1, j)*ut_ad(npx+1, j)
        uc_ad(npx+1, j) = uc_ad(npx+1, j) + temp_ad6
        v_ad(npx+1, j) = v_ad(npx+1, j) - cosa_u(npx+1, j)*temp_ad6
        ut_ad(npx+1, j) = 0.0
        CALL POPREALARRAY(ut(npx-1, j))
        temp_ad7 = rsin_u(npx-1, j)*ut_ad(npx-1, j)
        uc_ad(npx-1, j) = uc_ad(npx-1, j) + temp_ad7
        v_ad(npx-1, j) = v_ad(npx-1, j) - cosa_u(npx-1, j)*temp_ad7
        ut_ad(npx-1, j) = 0.0
        utmp_ad(npx, j) = utmp_ad(npx, j) + c3*uc_ad(npx+1, j)
        utmp_ad(npx+1, j) = utmp_ad(npx+1, j) + c2*uc_ad(npx+1, j)
        utmp_ad(npx+2, j) = utmp_ad(npx+2, j) + c1*uc_ad(npx+1, j)
        uc_ad(npx+1, j) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          ut_ad(npx, j) = ut_ad(npx, j) + sin_sg(npx-1, j, 3)*uc_ad(npx&
&           , j)
          uc_ad(npx, j) = 0.0
        ELSE
          ut_ad(npx, j) = ut_ad(npx, j) + sin_sg(npx, j, 1)*uc_ad(npx, j&
&           )
          uc_ad(npx, j) = 0.0
        END IF
        CALL EDGE_INTERPOLATE4_BWD(ua(npx-2:npx+1, j), ua_ad(npx-2:npx+1&
&                            , j), dxa(npx-2:npx+1, j), ut_ad(npx, j))
        ut_ad(npx, j) = 0.0
        CALL POPREALARRAY(ut(npx, j))
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
        CALL POPREALARRAY(ut(2, j))
        temp_ad4 = rsin_u(2, j)*ut_ad(2, j)
        uc_ad(2, j) = uc_ad(2, j) + temp_ad4
        v_ad(2, j) = v_ad(2, j) - cosa_u(2, j)*temp_ad4
        ut_ad(2, j) = 0.0
        CALL POPREALARRAY(ut(0, j))
        temp_ad5 = rsin_u(0, j)*ut_ad(0, j)
        uc_ad(0, j) = uc_ad(0, j) + temp_ad5
        v_ad(0, j) = v_ad(0, j) - cosa_u(0, j)*temp_ad5
        ut_ad(0, j) = 0.0
        utmp_ad(3, j) = utmp_ad(3, j) + c1*uc_ad(2, j)
        utmp_ad(2, j) = utmp_ad(2, j) + c2*uc_ad(2, j)
        utmp_ad(1, j) = utmp_ad(1, j) + c3*uc_ad(2, j)
        uc_ad(2, j) = 0.0
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          ut_ad(1, j) = ut_ad(1, j) + sin_sg(0, j, 3)*uc_ad(1, j)
          uc_ad(1, j) = 0.0
        ELSE
          ut_ad(1, j) = ut_ad(1, j) + sin_sg(1, j, 1)*uc_ad(1, j)
          uc_ad(1, j) = 0.0
        END IF
        CALL EDGE_INTERPOLATE4_BWD(ua(-1:2, j), ua_ad(-1:2, j), dxa(-1:2&
&                            , j), ut_ad(1, j))
        ut_ad(1, j) = 0.0
        CALL POPREALARRAY(ut(1, j))
        utmp_ad(-2, j) = utmp_ad(-2, j) + c1*uc_ad(0, j)
        utmp_ad(-1, j) = utmp_ad(-1, j) + c2*uc_ad(0, j)
        utmp_ad(0, j) = utmp_ad(0, j) + c3*uc_ad(0, j)
        uc_ad(0, j) = 0.0
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(ua(0, npy))
      va_ad(0, npy-1) = va_ad(0, npy-1) + ua_ad(0, npy)
      ua_ad(0, npy) = 0.0
      CALL POPREALARRAY(ua(-1, npy))
      va_ad(0, npy-2) = va_ad(0, npy-2) + ua_ad(-1, npy)
      ua_ad(-1, npy) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(ua(npx+1, npy))
      va_ad(npx, npy-2) = va_ad(npx, npy-2) - ua_ad(npx+1, npy)
      ua_ad(npx+1, npy) = 0.0
      CALL POPREALARRAY(ua(npx, npy))
      va_ad(npx, npy-1) = va_ad(npx, npy-1) - ua_ad(npx, npy)
      ua_ad(npx, npy) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(ua(npx+1, 0))
      va_ad(npx, 2) = va_ad(npx, 2) + ua_ad(npx+1, 0)
      ua_ad(npx+1, 0) = 0.0
      CALL POPREALARRAY(ua(npx, 0))
      va_ad(npx, 1) = va_ad(npx, 1) + ua_ad(npx, 0)
      ua_ad(npx, 0) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(ua(0, 0))
      va_ad(0, 1) = va_ad(0, 1) - ua_ad(0, 0)
      ua_ad(0, 0) = 0.0
      CALL POPREALARRAY(ua(-1, 0))
      va_ad(0, 2) = va_ad(0, 2) - ua_ad(-1, 0)
      ua_ad(-1, 0) = 0.0
    END IF
 100 DO j=je+1,js-1,-1
      DO i=ilast,ifirst,-1
        CALL POPREALARRAY(ut(i, j))
        temp_ad3 = rsin_u(i, j)*ut_ad(i, j)
        uc_ad(i, j) = uc_ad(i, j) + temp_ad3
        v_ad(i, j) = v_ad(i, j) - cosa_u(i, j)*temp_ad3
        ut_ad(i, j) = 0.0
        utmp_ad(i-2, j) = utmp_ad(i-2, j) + a2*uc_ad(i, j)
        utmp_ad(i+1, j) = utmp_ad(i+1, j) + a2*uc_ad(i, j)
        utmp_ad(i-1, j) = utmp_ad(i-1, j) + a1*uc_ad(i, j)
        utmp_ad(i, j) = utmp_ad(i, j) + a1*uc_ad(i, j)
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
    rsin2 => gridstruct%rsin2
    cosa_s => gridstruct%cosa_s
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=jed,jsd,-1
        DO i=ied,isd,-1
          temp_ad0 = rsin2(i, j)*ua_ad(i, j)
          CALL POPREALARRAY(va(i, j))
          temp_ad = rsin2(i, j)*va_ad(i, j)
          vtmp_ad(i, j) = vtmp_ad(i, j) + temp_ad - cosa_s(i, j)*&
&           temp_ad0
          utmp_ad(i, j) = utmp_ad(i, j) + temp_ad0 - cosa_s(i, j)*&
&           temp_ad
          va_ad(i, j) = 0.0
          CALL POPREALARRAY(ua(i, j))
          ua_ad(i, j) = 0.0
        END DO
      END DO
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
    ELSE
      DO j=je+id+1,js-1-id,-1
        DO i=ie+id+1,is-1-id,-1
          temp_ad2 = rsin2(i, j)*ua_ad(i, j)
          CALL POPREALARRAY(va(i, j))
          temp_ad1 = rsin2(i, j)*va_ad(i, j)
          vtmp_ad(i, j) = vtmp_ad(i, j) + temp_ad1 - cosa_s(i, j)*&
&           temp_ad2
          utmp_ad(i, j) = utmp_ad(i, j) + temp_ad2 - cosa_s(i, j)*&
&           temp_ad1
          va_ad(i, j) = 0.0
          CALL POPREALARRAY(ua(i, j))
          ua_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL(2,branch)
      IF (branch .NE. 0) THEN
        IF (branch .NE. 1) THEN
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
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) jed = bd%jed
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) jsd = bd%jsd
        END IF
        isd = bd%isd
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
    END IF
  END SUBROUTINE D2A2C_VECT_BWD
!There is a limit to how far this routine can fill uc and vc in the
! halo, and so either mpp_update_domains or some sort of boundary
!  routine (extrapolation, outflow, interpolation from a nested grid)
!   is needed after c_sw is completed if these variables are needed
!    in the halo
  SUBROUTINE D2A2C_VECT(u, v, ua, va, uc, vc, ut, vt, dord4, gridstruct&
&   , bd, npx, npy, nested, grid_type)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    LOGICAL, INTENT(IN) :: dord4
    REAL, INTENT(IN) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL, INTENT(IN) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(OUT) :: uc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(OUT) :: vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: ua, va&
&   , ut, vt
    INTEGER, INTENT(IN) :: npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! Local
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: utmp, vtmp
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsin2
    REAL, DIMENSION(:, :), POINTER :: dxa, dya
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
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    cosa_s => gridstruct%cosa_s
    rsin_u => gridstruct%rsin_u
    rsin_v => gridstruct%rsin_v
    rsin2 => gridstruct%rsin2
    dxa => gridstruct%dxa
    dya => gridstruct%dya
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
! Initialize the non-existing corner regions
    utmp(:, :) = big_number
    vtmp(:, :) = big_number
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
! Contra-variant components at cell center:
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
    IF (gridstruct%sw_corner) THEN
      DO i=-2,0
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
    END IF
    IF (gridstruct%se_corner) THEN
      DO i=0,2
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
    END IF
    IF (gridstruct%ne_corner) THEN
      DO i=0,2
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
    END IF
    IF (gridstruct%nw_corner) THEN
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
        uc(i, j) = a2*(utmp(i-2, j)+utmp(i+1, j)) + a1*(utmp(i-1, j)+&
&         utmp(i, j))
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
! Xdir:
      IF (gridstruct%sw_corner) THEN
        ua(-1, 0) = -va(0, 2)
        ua(0, 0) = -va(0, 1)
      END IF
      IF (gridstruct%se_corner) THEN
        ua(npx, 0) = va(npx, 1)
        ua(npx+1, 0) = va(npx, 2)
      END IF
      IF (gridstruct%ne_corner) THEN
        ua(npx, npy) = -va(npx, npy-1)
        ua(npx+1, npy) = -va(npx, npy-2)
      END IF
      IF (gridstruct%nw_corner) THEN
        ua(-1, npy) = va(0, npy-2)
        ua(0, npy) = va(0, npy-1)
      END IF
      IF (is .EQ. 1 .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc(0, j) = c1*utmp(-2, j) + c2*utmp(-1, j) + c3*utmp(0, j)
          ut(1, j) = EDGE_INTERPOLATE4(ua(-1:2, j), dxa(-1:2, j))
!Want to use the UPSTREAM value
          IF (ut(1, j) .GT. 0.) THEN
            uc(1, j) = ut(1, j)*sin_sg(0, j, 3)
          ELSE
            uc(1, j) = ut(1, j)*sin_sg(1, j, 1)
          END IF
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
          ut(0, j) = (uc(0, j)-v(0, j)*cosa_u(0, j))*rsin_u(0, j)
          ut(2, j) = (uc(2, j)-v(2, j)*cosa_u(2, j))*rsin_u(2, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
          ut(npx, j) = EDGE_INTERPOLATE4(ua(npx-2:npx+1, j), dxa(npx-2:&
&           npx+1, j))
          IF (ut(npx, j) .GT. 0.) THEN
            uc(npx, j) = ut(npx, j)*sin_sg(npx-1, j, 3)
          ELSE
            uc(npx, j) = ut(npx, j)*sin_sg(npx, j, 1)
          END IF
          uc(npx+1, j) = c3*utmp(npx, j) + c2*utmp(npx+1, j) + c1*utmp(&
&           npx+2, j)
          ut(npx-1, j) = (uc(npx-1, j)-v(npx-1, j)*cosa_u(npx-1, j))*&
&           rsin_u(npx-1, j)
          ut(npx+1, j) = (uc(npx+1, j)-v(npx+1, j)*cosa_u(npx+1, j))*&
&           rsin_u(npx+1, j)
        END DO
      END IF
    END IF
!------
! Ydir:
!------
    IF (gridstruct%sw_corner) THEN
      DO j=-2,0
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
    END IF
    IF (gridstruct%nw_corner) THEN
      DO j=0,2
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
    END IF
    IF (gridstruct%se_corner) THEN
      DO j=-2,0
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
    END IF
    IF (gridstruct%ne_corner) THEN
      DO j=0,2
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
    END IF
    IF (gridstruct%sw_corner) THEN
      va(0, -1) = -ua(2, 0)
      va(0, 0) = -ua(1, 0)
    END IF
    IF (gridstruct%se_corner) THEN
      va(npx, 0) = ua(npx-1, 0)
      va(npx, -1) = ua(npx-2, 0)
    END IF
    IF (gridstruct%ne_corner) THEN
      va(npx, npy) = -ua(npx-1, npy)
      va(npx, npy+1) = -ua(npx-2, npy)
    END IF
    IF (gridstruct%nw_corner) THEN
      va(0, npy) = ua(1, npy)
      va(0, npy+1) = ua(2, npy)
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1 .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vt(i, j) = EDGE_INTERPOLATE4(va(i, -1:2), dya(i, -1:2))
            IF (vt(i, j) .GT. 0.) THEN
              vc(i, j) = vt(i, j)*sin_sg(i, j-1, 4)
            ELSE
              vc(i, j) = vt(i, j)*sin_sg(i, j, 2)
            END IF
          END DO
        ELSE IF (j .EQ. 0 .OR. (j .EQ. npy - 1 .AND. (.NOT.nested))) &
&       THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        ELSE IF (j .EQ. 2 .OR. (j .EQ. npy + 1 .AND. (.NOT.nested))) &
&       THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j+1) + c2*vtmp(i, j) + c3*vtmp(i, j-1)
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        ELSE IF (j .EQ. npy .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vt(i, j) = EDGE_INTERPOLATE4(va(i, j-2:j+1), dya(i, j-2:j+1)&
&             )
            IF (vt(i, j) .GT. 0.) THEN
              vc(i, j) = vt(i, j)*sin_sg(i, j-1, 4)
            ELSE
              vc(i, j) = vt(i, j)*sin_sg(i, j, 2)
            END IF
          END DO
        ELSE
! 4th order interpolation for interior points:
          DO i=is-1,ie+1
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        END IF
      END DO
    ELSE
! 4th order interpolation:
      DO j=js-1,je+2
        DO i=is-1,ie+1
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
          vt(i, j) = vc(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE D2A2C_VECT
!  Differentiation of edge_interpolate4 in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b
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
!   gradient     of useful results: ua edge_interpolate4
!   with respect to varying inputs: ua
  REAL FUNCTION EDGE_INTERPOLATE4_FWD(ua, dxa)
    IMPLICIT NONE
    REAL, INTENT(IN) :: ua(4)
    REAL, INTENT(IN) :: dxa(4)
    REAL :: t1, t2
    REAL :: edge_interpolate4
    t1 = dxa(1) + dxa(2)
    t2 = dxa(3) + dxa(4)
    edge_interpolate4 = 0.5*(((t1+dxa(2))*ua(2)-dxa(2)*ua(1))/t1+((t2+&
&     dxa(3))*ua(3)-dxa(3)*ua(4))/t2)
    edge_interpolate4_fwd = edge_interpolate4
  END FUNCTION EDGE_INTERPOLATE4_FWD
!  Differentiation of edge_interpolate4 in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2
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
!   gradient     of useful results: ua edge_interpolate4
!   with respect to varying inputs: ua
  SUBROUTINE EDGE_INTERPOLATE4_BWD(ua, ua_ad, dxa, edge_interpolate4_ad)
    IMPLICIT NONE
    REAL, INTENT(IN) :: ua(4)
    REAL :: ua_ad(4)
    REAL, INTENT(IN) :: dxa(4)
    REAL :: t1, t2
    REAL :: temp_ad
    REAL :: edge_interpolate4_ad
    REAL :: edge_interpolate4
    t1 = dxa(1) + dxa(2)
    t2 = dxa(3) + dxa(4)
    temp_ad = 0.5*edge_interpolate4_ad
    ua_ad(2) = ua_ad(2) + (t1+dxa(2))*temp_ad/t1
    ua_ad(1) = ua_ad(1) - dxa(2)*temp_ad/t1
    ua_ad(3) = ua_ad(3) + (t2+dxa(3))*temp_ad/t2
    ua_ad(4) = ua_ad(4) - dxa(3)*temp_ad/t2
  END SUBROUTINE EDGE_INTERPOLATE4_BWD
  REAL FUNCTION EDGE_INTERPOLATE4(ua, dxa)
    IMPLICIT NONE
    REAL, INTENT(IN) :: ua(4)
    REAL, INTENT(IN) :: dxa(4)
    REAL :: t1, t2
    t1 = dxa(1) + dxa(2)
    t2 = dxa(3) + dxa(4)
    edge_interpolate4 = 0.5*(((t1+dxa(2))*ua(2)-dxa(2)*ua(1))/t1+((t2+&
&     dxa(3))*ua(3)-dxa(3)*ua(4))/t2)
  END FUNCTION EDGE_INTERPOLATE4
  SUBROUTINE FILL3_4CORNERS(q1, q2, q3, dir, bd, npx, npy, sw_corner, &
&   se_corner, ne_corner, nw_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q3(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER :: i, j
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
    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        q1(-1, 0) = q1(0, 2)
        q1(0, 0) = q1(0, 1)
        q1(0, -1) = q1(-1, 1)
        q2(-1, 0) = q2(0, 2)
        q2(0, 0) = q2(0, 1)
        q2(0, -1) = q2(-1, 1)
        q3(-1, 0) = q3(0, 2)
        q3(0, 0) = q3(0, 1)
        q3(0, -1) = q3(-1, 1)
      END IF
      IF (se_corner) THEN
        q1(npx+1, 0) = q1(npx, 2)
        q1(npx, 0) = q1(npx, 1)
        q1(npx, -1) = q1(npx+1, 1)
        q2(npx+1, 0) = q2(npx, 2)
        q2(npx, 0) = q2(npx, 1)
        q2(npx, -1) = q2(npx+1, 1)
        q3(npx+1, 0) = q3(npx, 2)
        q3(npx, 0) = q3(npx, 1)
        q3(npx, -1) = q3(npx+1, 1)
      END IF
      IF (ne_corner) THEN
        q1(npx, npy) = q1(npx, npy-1)
        q1(npx+1, npy) = q1(npx, npy-2)
        q1(npx, npy+1) = q1(npx+1, npy-1)
        q2(npx, npy) = q2(npx, npy-1)
        q2(npx+1, npy) = q2(npx, npy-2)
        q2(npx, npy+1) = q2(npx+1, npy-1)
        q3(npx, npy) = q3(npx, npy-1)
        q3(npx+1, npy) = q3(npx, npy-2)
        q3(npx, npy+1) = q3(npx+1, npy-1)
      END IF
      IF (nw_corner) THEN
        q1(0, npy) = q1(0, npy-1)
        q1(-1, npy) = q1(0, npy-2)
        q1(0, npy+1) = q1(-1, npy-1)
        q2(0, npy) = q2(0, npy-1)
        q2(-1, npy) = q2(0, npy-2)
        q2(0, npy+1) = q2(-1, npy-1)
        q3(0, npy) = q3(0, npy-1)
        q3(-1, npy) = q3(0, npy-2)
        q3(0, npy+1) = q3(-1, npy-1)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q1(0, 0) = q1(1, 0)
        q1(0, -1) = q1(2, 0)
        q1(-1, 0) = q1(1, -1)
        q2(0, 0) = q2(1, 0)
        q2(0, -1) = q2(2, 0)
        q2(-1, 0) = q2(1, -1)
        q3(0, 0) = q3(1, 0)
        q3(0, -1) = q3(2, 0)
        q3(-1, 0) = q3(1, -1)
      END IF
      IF (se_corner) THEN
        q1(npx, 0) = q1(npx-1, 0)
        q1(npx, -1) = q1(npx-2, 0)
        q1(npx+1, 0) = q1(npx-1, -1)
        q2(npx, 0) = q2(npx-1, 0)
        q2(npx, -1) = q2(npx-2, 0)
        q2(npx+1, 0) = q2(npx-1, -1)
        q3(npx, 0) = q3(npx-1, 0)
        q3(npx, -1) = q3(npx-2, 0)
        q3(npx+1, 0) = q3(npx-1, -1)
      END IF
      IF (ne_corner) THEN
        q1(npx, npy) = q1(npx-1, npy)
        q1(npx, npy+1) = q1(npx-2, npy)
        q1(npx+1, npy) = q1(npx-1, npy+1)
        q2(npx, npy) = q2(npx-1, npy)
        q2(npx, npy+1) = q2(npx-2, npy)
        q2(npx+1, npy) = q2(npx-1, npy+1)
        q3(npx, npy) = q3(npx-1, npy)
        q3(npx, npy+1) = q3(npx-2, npy)
        q3(npx+1, npy) = q3(npx-1, npy+1)
      END IF
      IF (nw_corner) THEN
        q1(0, npy) = q1(1, npy)
        q1(0, npy+1) = q1(2, npy)
        q1(-1, npy) = q1(1, npy+1)
        q2(0, npy) = q2(1, npy)
        q2(0, npy+1) = q2(2, npy)
        q2(-1, npy) = q2(1, npy+1)
        q3(0, npy) = q3(1, npy)
        q3(0, npy+1) = q3(2, npy)
        q3(-1, npy) = q3(1, npy+1)
      END IF
    END SELECT
  END SUBROUTINE FILL3_4CORNERS
!  Differentiation of fill2_4corners in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
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
!   gradient     of useful results: q1 q2
!   with respect to varying inputs: q1 q2
  SUBROUTINE FILL2_4CORNERS_FWD(q1, q2, dir, bd, npx, npy, sw_corner, &
&   se_corner, ne_corner, nw_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q2(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed

    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0

    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        CALL PUSHREALARRAY(q1(-1, 0))
        q1(-1, 0) = q1(0, 2)
        CALL PUSHREALARRAY(q1(0, 0))
        q1(0, 0) = q1(0, 1)
        CALL PUSHREALARRAY(q2(-1, 0))
        q2(-1, 0) = q2(0, 2)
        CALL PUSHREALARRAY(q2(0, 0))
        q2(0, 0) = q2(0, 1)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (se_corner) THEN
        CALL PUSHREALARRAY(q1(npx+1, 0))
        q1(npx+1, 0) = q1(npx, 2)
        CALL PUSHREALARRAY(q1(npx, 0))
        q1(npx, 0) = q1(npx, 1)
        CALL PUSHREALARRAY(q2(npx+1, 0))
        q2(npx+1, 0) = q2(npx, 2)
        CALL PUSHREALARRAY(q2(npx, 0))
        q2(npx, 0) = q2(npx, 1)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (nw_corner) THEN
        CALL PUSHREALARRAY(q1(0, npy))
        q1(0, npy) = q1(0, npy-1)
        CALL PUSHREALARRAY(q1(-1, npy))
        q1(-1, npy) = q1(0, npy-2)
        CALL PUSHREALARRAY(q2(0, npy))
        q2(0, npy) = q2(0, npy-1)
        CALL PUSHREALARRAY(q2(-1, npy))
        q2(-1, npy) = q2(0, npy-2)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHREALARRAY(q1(npx, npy))
        q1(npx, npy) = q1(npx, npy-1)
        CALL PUSHREALARRAY(q1(npx+1, npy))
        q1(npx+1, npy) = q1(npx, npy-2)
        CALL PUSHREALARRAY(q2(npx, npy))
        q2(npx, npy) = q2(npx, npy-1)
        CALL PUSHREALARRAY(q2(npx+1, npy))
        q2(npx+1, npy) = q2(npx, npy-2)
        CALL PUSHCONTROL(3,2)
      ELSE
        CALL PUSHCONTROL(3,1)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        CALL PUSHREALARRAY(q1(0, 0))
        q1(0, 0) = q1(1, 0)
        CALL PUSHREALARRAY(q1(0, -1))
        q1(0, -1) = q1(2, 0)
        CALL PUSHREALARRAY(q2(0, 0))
        q2(0, 0) = q2(1, 0)
        CALL PUSHREALARRAY(q2(0, -1))
        q2(0, -1) = q2(2, 0)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (se_corner) THEN
        CALL PUSHREALARRAY(q1(npx, 0))
        q1(npx, 0) = q1(npx-1, 0)
        CALL PUSHREALARRAY(q1(npx, -1))
        q1(npx, -1) = q1(npx-2, 0)
        CALL PUSHREALARRAY(q2(npx, 0))
        q2(npx, 0) = q2(npx-1, 0)
        CALL PUSHREALARRAY(q2(npx, -1))
        q2(npx, -1) = q2(npx-2, 0)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (nw_corner) THEN
        CALL PUSHREALARRAY(q1(0, npy))
        q1(0, npy) = q1(1, npy)
        CALL PUSHREALARRAY(q1(0, npy+1))
        q1(0, npy+1) = q1(2, npy)
        CALL PUSHREALARRAY(q2(0, npy))
        q2(0, npy) = q2(1, npy)
        CALL PUSHREALARRAY(q2(0, npy+1))
        q2(0, npy+1) = q2(2, npy)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHREALARRAY(q1(npx, npy))
        q1(npx, npy) = q1(npx-1, npy)
        CALL PUSHREALARRAY(q1(npx, npy+1))
        q1(npx, npy+1) = q1(npx-2, npy)
        CALL PUSHREALARRAY(q2(npx, npy))
        q2(npx, npy) = q2(npx-1, npy)
        CALL PUSHREALARRAY(q2(npx, npy+1))
        q2(npx, npy+1) = q2(npx-2, npy)
        CALL PUSHCONTROL(3,4)
      ELSE
        CALL PUSHCONTROL(3,3)
      END IF
    CASE DEFAULT
      CALL PUSHCONTROL(3,0)
    END SELECT
  END SUBROUTINE FILL2_4CORNERS_FWD
!  Differentiation of fill2_4corners in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_e
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
!   gradient     of useful results: q1 q2
!   with respect to varying inputs: q1 q2
  SUBROUTINE FILL2_4CORNERS_BWD(q1, q1_ad, q2, q2_ad, dir, bd, npx, npy&
&   , sw_corner, se_corner, ne_corner, nw_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q1_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q2_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: branch

    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    branch = 0

    CALL POPCONTROL(3,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) GOTO 100
    ELSE IF (branch .EQ. 2) THEN
      CALL POPREALARRAY(q2(npx+1, npy))
      q2_ad(npx, npy-2) = q2_ad(npx, npy-2) + q2_ad(npx+1, npy)
      q2_ad(npx+1, npy) = 0.0
      CALL POPREALARRAY(q2(npx, npy))
      q2_ad(npx, npy-1) = q2_ad(npx, npy-1) + q2_ad(npx, npy)
      q2_ad(npx, npy) = 0.0
      CALL POPREALARRAY(q1(npx+1, npy))
      q1_ad(npx, npy-2) = q1_ad(npx, npy-2) + q1_ad(npx+1, npy)
      q1_ad(npx+1, npy) = 0.0
      CALL POPREALARRAY(q1(npx, npy))
      q1_ad(npx, npy-1) = q1_ad(npx, npy-1) + q1_ad(npx, npy)
      q1_ad(npx, npy) = 0.0
    ELSE
      IF (branch .NE. 3) THEN
        CALL POPREALARRAY(q2(npx, npy+1))
        q2_ad(npx-2, npy) = q2_ad(npx-2, npy) + q2_ad(npx, npy+1)
        q2_ad(npx, npy+1) = 0.0
        CALL POPREALARRAY(q2(npx, npy))
        q2_ad(npx-1, npy) = q2_ad(npx-1, npy) + q2_ad(npx, npy)
        q2_ad(npx, npy) = 0.0
        CALL POPREALARRAY(q1(npx, npy+1))
        q1_ad(npx-2, npy) = q1_ad(npx-2, npy) + q1_ad(npx, npy+1)
        q1_ad(npx, npy+1) = 0.0
        CALL POPREALARRAY(q1(npx, npy))
        q1_ad(npx-1, npy) = q1_ad(npx-1, npy) + q1_ad(npx, npy)
        q1_ad(npx, npy) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(q2(0, npy+1))
        q2_ad(2, npy) = q2_ad(2, npy) + q2_ad(0, npy+1)
        q2_ad(0, npy+1) = 0.0
        CALL POPREALARRAY(q2(0, npy))
        q2_ad(1, npy) = q2_ad(1, npy) + q2_ad(0, npy)
        q2_ad(0, npy) = 0.0
        CALL POPREALARRAY(q1(0, npy+1))
        q1_ad(2, npy) = q1_ad(2, npy) + q1_ad(0, npy+1)
        q1_ad(0, npy+1) = 0.0
        CALL POPREALARRAY(q1(0, npy))
        q1_ad(1, npy) = q1_ad(1, npy) + q1_ad(0, npy)
        q1_ad(0, npy) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(q2(npx, -1))
        q2_ad(npx-2, 0) = q2_ad(npx-2, 0) + q2_ad(npx, -1)
        q2_ad(npx, -1) = 0.0
        CALL POPREALARRAY(q2(npx, 0))
        q2_ad(npx-1, 0) = q2_ad(npx-1, 0) + q2_ad(npx, 0)
        q2_ad(npx, 0) = 0.0
        CALL POPREALARRAY(q1(npx, -1))
        q1_ad(npx-2, 0) = q1_ad(npx-2, 0) + q1_ad(npx, -1)
        q1_ad(npx, -1) = 0.0
        CALL POPREALARRAY(q1(npx, 0))
        q1_ad(npx-1, 0) = q1_ad(npx-1, 0) + q1_ad(npx, 0)
        q1_ad(npx, 0) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(q2(0, -1))
        q2_ad(2, 0) = q2_ad(2, 0) + q2_ad(0, -1)
        q2_ad(0, -1) = 0.0
        CALL POPREALARRAY(q2(0, 0))
        q2_ad(1, 0) = q2_ad(1, 0) + q2_ad(0, 0)
        q2_ad(0, 0) = 0.0
        CALL POPREALARRAY(q1(0, -1))
        q1_ad(2, 0) = q1_ad(2, 0) + q1_ad(0, -1)
        q1_ad(0, -1) = 0.0
        CALL POPREALARRAY(q1(0, 0))
        q1_ad(1, 0) = q1_ad(1, 0) + q1_ad(0, 0)
        q1_ad(0, 0) = 0.0
      END IF
      GOTO 100
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(q2(-1, npy))
      q2_ad(0, npy-2) = q2_ad(0, npy-2) + q2_ad(-1, npy)
      q2_ad(-1, npy) = 0.0
      CALL POPREALARRAY(q2(0, npy))
      q2_ad(0, npy-1) = q2_ad(0, npy-1) + q2_ad(0, npy)
      q2_ad(0, npy) = 0.0
      CALL POPREALARRAY(q1(-1, npy))
      q1_ad(0, npy-2) = q1_ad(0, npy-2) + q1_ad(-1, npy)
      q1_ad(-1, npy) = 0.0
      CALL POPREALARRAY(q1(0, npy))
      q1_ad(0, npy-1) = q1_ad(0, npy-1) + q1_ad(0, npy)
      q1_ad(0, npy) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(q2(npx, 0))
      q2_ad(npx, 1) = q2_ad(npx, 1) + q2_ad(npx, 0)
      q2_ad(npx, 0) = 0.0
      CALL POPREALARRAY(q2(npx+1, 0))
      q2_ad(npx, 2) = q2_ad(npx, 2) + q2_ad(npx+1, 0)
      q2_ad(npx+1, 0) = 0.0
      CALL POPREALARRAY(q1(npx, 0))
      q1_ad(npx, 1) = q1_ad(npx, 1) + q1_ad(npx, 0)
      q1_ad(npx, 0) = 0.0
      CALL POPREALARRAY(q1(npx+1, 0))
      q1_ad(npx, 2) = q1_ad(npx, 2) + q1_ad(npx+1, 0)
      q1_ad(npx+1, 0) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(q2(0, 0))
      q2_ad(0, 1) = q2_ad(0, 1) + q2_ad(0, 0)
      q2_ad(0, 0) = 0.0
      CALL POPREALARRAY(q2(-1, 0))
      q2_ad(0, 2) = q2_ad(0, 2) + q2_ad(-1, 0)
      q2_ad(-1, 0) = 0.0
      CALL POPREALARRAY(q1(0, 0))
      q1_ad(0, 1) = q1_ad(0, 1) + q1_ad(0, 0)
      q1_ad(0, 0) = 0.0
      CALL POPREALARRAY(q1(-1, 0))
      q1_ad(0, 2) = q1_ad(0, 2) + q1_ad(-1, 0)
      q1_ad(-1, 0) = 0.0
    END IF
 100 CONTINUE
  END SUBROUTINE FILL2_4CORNERS_BWD
  SUBROUTINE FILL2_4CORNERS(q1, q2, dir, bd, npx, npy, sw_corner, &
&   se_corner, ne_corner, nw_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q2(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    INTEGER, INTENT(IN) :: npx, npy
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
    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        q1(-1, 0) = q1(0, 2)
        q1(0, 0) = q1(0, 1)
        q2(-1, 0) = q2(0, 2)
        q2(0, 0) = q2(0, 1)
      END IF
      IF (se_corner) THEN
        q1(npx+1, 0) = q1(npx, 2)
        q1(npx, 0) = q1(npx, 1)
        q2(npx+1, 0) = q2(npx, 2)
        q2(npx, 0) = q2(npx, 1)
      END IF
      IF (nw_corner) THEN
        q1(0, npy) = q1(0, npy-1)
        q1(-1, npy) = q1(0, npy-2)
        q2(0, npy) = q2(0, npy-1)
        q2(-1, npy) = q2(0, npy-2)
      END IF
      IF (ne_corner) THEN
        q1(npx, npy) = q1(npx, npy-1)
        q1(npx+1, npy) = q1(npx, npy-2)
        q2(npx, npy) = q2(npx, npy-1)
        q2(npx+1, npy) = q2(npx, npy-2)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q1(0, 0) = q1(1, 0)
        q1(0, -1) = q1(2, 0)
        q2(0, 0) = q2(1, 0)
        q2(0, -1) = q2(2, 0)
      END IF
      IF (se_corner) THEN
        q1(npx, 0) = q1(npx-1, 0)
        q1(npx, -1) = q1(npx-2, 0)
        q2(npx, 0) = q2(npx-1, 0)
        q2(npx, -1) = q2(npx-2, 0)
      END IF
      IF (nw_corner) THEN
        q1(0, npy) = q1(1, npy)
        q1(0, npy+1) = q1(2, npy)
        q2(0, npy) = q2(1, npy)
        q2(0, npy+1) = q2(2, npy)
      END IF
      IF (ne_corner) THEN
        q1(npx, npy) = q1(npx-1, npy)
        q1(npx, npy+1) = q1(npx-2, npy)
        q2(npx, npy) = q2(npx-1, npy)
        q2(npx, npy+1) = q2(npx-2, npy)
      END IF
    END SELECT
  END SUBROUTINE FILL2_4CORNERS
!  Differentiation of fill_4corners in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE FILL_4CORNERS_FWD(q, dir, bd, npx, npy, sw_corner, &
&   se_corner, ne_corner, nw_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed

    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0 

    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        CALL PUSHREALARRAY(q(-1, 0))
        q(-1, 0) = q(0, 2)
        CALL PUSHREALARRAY(q(0, 0))
        q(0, 0) = q(0, 1)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (se_corner) THEN
        CALL PUSHREALARRAY(q(npx+1, 0))
        q(npx+1, 0) = q(npx, 2)
        CALL PUSHREALARRAY(q(npx, 0))
        q(npx, 0) = q(npx, 1)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (nw_corner) THEN
        CALL PUSHREALARRAY(q(0, npy))
        q(0, npy) = q(0, npy-1)
        CALL PUSHREALARRAY(q(-1, npy))
        q(-1, npy) = q(0, npy-2)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHREALARRAY(q(npx, npy))
        q(npx, npy) = q(npx, npy-1)
        CALL PUSHREALARRAY(q(npx+1, npy))
        q(npx+1, npy) = q(npx, npy-2)
        CALL PUSHCONTROL(3,2)
      ELSE
        CALL PUSHCONTROL(3,1)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        CALL PUSHREALARRAY(q(0, 0))
        q(0, 0) = q(1, 0)
        CALL PUSHREALARRAY(q(0, -1))
        q(0, -1) = q(2, 0)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (se_corner) THEN
        CALL PUSHREALARRAY(q(npx, 0))
        q(npx, 0) = q(npx-1, 0)
        CALL PUSHREALARRAY(q(npx, -1))
        q(npx, -1) = q(npx-2, 0)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (nw_corner) THEN
        CALL PUSHREALARRAY(q(0, npy))
        q(0, npy) = q(1, npy)
        CALL PUSHREALARRAY(q(0, npy+1))
        q(0, npy+1) = q(2, npy)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHREALARRAY(q(npx, npy))
        q(npx, npy) = q(npx-1, npy)
        CALL PUSHREALARRAY(q(npx, npy+1))
        q(npx, npy+1) = q(npx-2, npy)
        CALL PUSHCONTROL(3,4)
      ELSE
        CALL PUSHCONTROL(3,3)
      END IF
    CASE DEFAULT
      CALL PUSHCONTROL(3,0)
    END SELECT
  END SUBROUTINE FILL_4CORNERS_FWD
!  Differentiation of fill_4corners in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
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
!   gradient     of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE FILL_4CORNERS_BWD(q, q_ad, dir, bd, npx, npy, sw_corner, &
&   se_corner, ne_corner, nw_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTEGER :: branch

    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    branch = 0

    CALL POPCONTROL(3,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) GOTO 100
    ELSE IF (branch .EQ. 2) THEN
      CALL POPREALARRAY(q(npx+1, npy))
      q_ad(npx, npy-2) = q_ad(npx, npy-2) + q_ad(npx+1, npy)
      q_ad(npx+1, npy) = 0.0
      CALL POPREALARRAY(q(npx, npy))
      q_ad(npx, npy-1) = q_ad(npx, npy-1) + q_ad(npx, npy)
      q_ad(npx, npy) = 0.0
    ELSE
      IF (branch .NE. 3) THEN
        CALL POPREALARRAY(q(npx, npy+1))
        q_ad(npx-2, npy) = q_ad(npx-2, npy) + q_ad(npx, npy+1)
        q_ad(npx, npy+1) = 0.0
        CALL POPREALARRAY(q(npx, npy))
        q_ad(npx-1, npy) = q_ad(npx-1, npy) + q_ad(npx, npy)
        q_ad(npx, npy) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(q(0, npy+1))
        q_ad(2, npy) = q_ad(2, npy) + q_ad(0, npy+1)
        q_ad(0, npy+1) = 0.0
        CALL POPREALARRAY(q(0, npy))
        q_ad(1, npy) = q_ad(1, npy) + q_ad(0, npy)
        q_ad(0, npy) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(q(npx, -1))
        q_ad(npx-2, 0) = q_ad(npx-2, 0) + q_ad(npx, -1)
        q_ad(npx, -1) = 0.0
        CALL POPREALARRAY(q(npx, 0))
        q_ad(npx-1, 0) = q_ad(npx-1, 0) + q_ad(npx, 0)
        q_ad(npx, 0) = 0.0
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(q(0, -1))
        q_ad(2, 0) = q_ad(2, 0) + q_ad(0, -1)
        q_ad(0, -1) = 0.0
        CALL POPREALARRAY(q(0, 0))
        q_ad(1, 0) = q_ad(1, 0) + q_ad(0, 0)
        q_ad(0, 0) = 0.0
      END IF
      GOTO 100
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(q(-1, npy))
      q_ad(0, npy-2) = q_ad(0, npy-2) + q_ad(-1, npy)
      q_ad(-1, npy) = 0.0
      CALL POPREALARRAY(q(0, npy))
      q_ad(0, npy-1) = q_ad(0, npy-1) + q_ad(0, npy)
      q_ad(0, npy) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(q(npx, 0))
      q_ad(npx, 1) = q_ad(npx, 1) + q_ad(npx, 0)
      q_ad(npx, 0) = 0.0
      CALL POPREALARRAY(q(npx+1, 0))
      q_ad(npx, 2) = q_ad(npx, 2) + q_ad(npx+1, 0)
      q_ad(npx+1, 0) = 0.0
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY(q(0, 0))
      q_ad(0, 1) = q_ad(0, 1) + q_ad(0, 0)
      q_ad(0, 0) = 0.0
      CALL POPREALARRAY(q(-1, 0))
      q_ad(0, 2) = q_ad(0, 2) + q_ad(-1, 0)
      q_ad(-1, 0) = 0.0
    END IF
 100 CONTINUE
  END SUBROUTINE FILL_4CORNERS_BWD
  SUBROUTINE FILL_4CORNERS(q, dir, bd, npx, npy, sw_corner, se_corner, &
&   ne_corner, nw_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    INTEGER, INTENT(IN) :: npx, npy
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
    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        q(-1, 0) = q(0, 2)
        q(0, 0) = q(0, 1)
      END IF
      IF (se_corner) THEN
        q(npx+1, 0) = q(npx, 2)
        q(npx, 0) = q(npx, 1)
      END IF
      IF (nw_corner) THEN
        q(0, npy) = q(0, npy-1)
        q(-1, npy) = q(0, npy-2)
      END IF
      IF (ne_corner) THEN
        q(npx, npy) = q(npx, npy-1)
        q(npx+1, npy) = q(npx, npy-2)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q(0, 0) = q(1, 0)
        q(0, -1) = q(2, 0)
      END IF
      IF (se_corner) THEN
        q(npx, 0) = q(npx-1, 0)
        q(npx, -1) = q(npx-2, 0)
      END IF
      IF (nw_corner) THEN
        q(0, npy) = q(1, npy)
        q(0, npy+1) = q(2, npy)
      END IF
      IF (ne_corner) THEN
        q(npx, npy) = q(npx-1, npy)
        q(npx, npy+1) = q(npx-2, npy)
      END IF
    END SELECT
  END SUBROUTINE FILL_4CORNERS
!  Differentiation of xtp_u in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
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
!   gradient     of useful results: flux u c
!   with respect to varying inputs: flux u c
  SUBROUTINE XTP_U_FWD(is, ie, js, je, isd, ied, jsd, jed, c, u, v, &
&   flux, iord, dx, rdx, npx, npy, grid_type, nested)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL :: flux(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: rdx(isd:ied, jsd:jed+1)
    INTEGER, INTENT(IN) :: iord, npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
! Local
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: al(is-1:ie+2), dm(is-2:ie+2)
    REAL :: dq(is-3:ie+2)
    REAL :: dl, dr, xt, pmp, lac, cfl
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: x0, x1, x0l, x0r
    INTEGER :: i, j
    INTEGER :: is3, ie3
    INTEGER :: is2, ie2
    INTRINSIC MAX
    INTRINSIC MIN

    bl = 0.0
    br = 0.0
    b0 = 0.0
    fx0 = 0.0
    al = 0.0
    dm = 0.0
    dq = 0.0
    dl = 0.0
    dr = 0.0
    xt = 0.0
    pmp = 0.0
    lac = 0.0
    cfl = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    x0 = 0.0
    x1 = 0.0
    x0l = 0.0
    x0r = 0.0
    is3 = 0
    ie3 = 0
    is2 = 0
    ie2 = 0

    IF (nested .OR. grid_type .GT. 3) THEN
      CALL PUSHCONTROL(1,0)
      is3 = is - 1
      ie3 = ie + 1
    ELSE
      IF (3 .LT. is - 1) THEN
        is3 = is - 1
      ELSE
        is3 = 3
      END IF
      IF (npx - 3 .GT. ie + 1) THEN
        CALL PUSHCONTROL(1,1)
        ie3 = ie + 1
      ELSE
        CALL PUSHCONTROL(1,1)
        ie3 = npx - 3
      END IF
    END IF
    IF (iord .EQ. 1) THEN
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHREALARRAY(flux(i, j))
            flux(i, j) = u(i-1, j)
            CALL PUSHCONTROL(1,1)
          ELSE
            CALL PUSHREALARRAY(flux(i, j))
            flux(i, j) = u(i, j)
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      CALL PUSHCONTROL(2,0)
    ELSE IF (iord .EQ. 333) THEN
      CALL PUSHREALARRAY(flux, (ie-is+2)*(je-js+2))
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux(i, j) = (2.0*u(i, j)+5.0*u(i-1, j)-u(i-2, j))/6.0 - 0.5&
&             *c(i, j)*rdx(i-1, j)*(u(i, j)-u(i-1, j)) + c(i, j)*rdx(i-1&
&             , j)*c(i, j)*rdx(i-1, j)/6.0*(u(i, j)-2.0*u(i-1, j)+u(i-2&
&             , j))
          ELSE
            flux(i, j) = (2.0*u(i-1, j)+5.0*u(i, j)-u(i+1, j))/6.0 - 0.5&
&             *c(i, j)*rdx(i, j)*(u(i, j)-u(i-1, j)) + c(i, j)*rdx(i, j)&
&             *c(i, j)*rdx(i, j)/6.0*(u(i+1, j)-2.0*u(i, j)+u(i-1, j))
          END IF
        END DO
      END DO
      CALL PUSHCONTROL(2,1)
    ELSE IF (iord .LT. 8) THEN
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
      DO j=js,je+1
        DO i=is3,ie3+1
          al(i) = p1*(u(i-1, j)+u(i, j)) + p2*(u(i-2, j)+u(i+1, j))
        END DO
        DO i=is3,ie3
          CALL PUSHREALARRAY(bl(i))
          bl(i) = al(i) - u(i, j)
          CALL PUSHREALARRAY(br(i))
          br(i) = al(i+1) - u(i, j)
        END DO
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            xt = c3*u(1, j) + c2*u(2, j) + c1*u(3, j)
            CALL PUSHREALARRAY(br(1))
            br(1) = xt - u(1, j)
            CALL PUSHREALARRAY(bl(2))
            bl(2) = xt - u(2, j)
            CALL PUSHREALARRAY(br(2))
            br(2) = al(3) - u(2, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! out
              CALL PUSHREALARRAY(bl(0))
              bl(0) = 0.
! edge
              CALL PUSHREALARRAY(br(0))
              br(0) = 0.
! edge
              CALL PUSHREALARRAY(bl(1))
              bl(1) = 0.
! in
              CALL PUSHREALARRAY(br(1))
              br(1) = 0.
              CALL PUSHCONTROL(2,0)
            ELSE
              CALL PUSHREALARRAY(bl(0))
              bl(0) = c1*u(-2, j) + c2*u(-1, j) + c3*u(0, j) - u(0, j)
              xt = 0.5*(((2.*dx(0, j)+dx(-1, j))*u(0, j)-dx(0, j)*u(-1, &
&               j))/(dx(0, j)+dx(-1, j))+((2.*dx(1, j)+dx(2, j))*u(1, j)&
&               -dx(1, j)*u(2, j))/(dx(1, j)+dx(2, j)))
              CALL PUSHREALARRAY(br(0))
              br(0) = xt - u(0, j)
              CALL PUSHREALARRAY(bl(1))
              bl(1) = xt - u(1, j)
              CALL PUSHCONTROL(2,1)
            END IF
          ELSE
            CALL PUSHCONTROL(2,2)
          END IF
!       call pert_ppm(1, u(2,j), bl(2), br(2), -1)
          IF (ie + 1 .EQ. npx) THEN
            CALL PUSHREALARRAY(bl(npx-2))
            bl(npx-2) = al(npx-2) - u(npx-2, j)
            xt = c1*u(npx-3, j) + c2*u(npx-2, j) + c3*u(npx-1, j)
            CALL PUSHREALARRAY(br(npx-2))
            br(npx-2) = xt - u(npx-2, j)
            CALL PUSHREALARRAY(bl(npx-1))
            bl(npx-1) = xt - u(npx-1, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! in
              CALL PUSHREALARRAY(bl(npx-1))
              bl(npx-1) = 0.
! edge
              CALL PUSHREALARRAY(br(npx-1))
              br(npx-1) = 0.
! edge
              CALL PUSHREALARRAY(bl(npx))
              bl(npx) = 0.
! out
              CALL PUSHREALARRAY(br(npx))
              br(npx) = 0.
              CALL PUSHCONTROL(2,3)
            ELSE
              xt = 0.5*(((2.*dx(npx-1, j)+dx(npx-2, j))*u(npx-1, j)-dx(&
&               npx-1, j)*u(npx-2, j))/(dx(npx-1, j)+dx(npx-2, j))+((2.*&
&               dx(npx, j)+dx(npx+1, j))*u(npx, j)-dx(npx, j)*u(npx+1, j&
&               ))/(dx(npx, j)+dx(npx+1, j)))
              CALL PUSHREALARRAY(br(npx-1))
              br(npx-1) = xt - u(npx-1, j)
              CALL PUSHREALARRAY(bl(npx))
              bl(npx) = xt - u(npx, j)
              CALL PUSHREALARRAY(br(npx))
              br(npx) = c3*u(npx, j) + c2*u(npx+1, j) + c1*u(npx+2, j) -&
&               u(npx, j)
              CALL PUSHCONTROL(2,2)
            END IF
          ELSE
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE
          CALL PUSHCONTROL(2,0)
        END IF
!       call pert_ppm(1, u(npx-2,j), bl(npx-2), br(npx-2), -1)
        DO i=is-1,ie+1
          CALL PUSHREALARRAY(b0(i))
          b0(i) = bl(i) + br(i)
        END DO
        IF (iord .EQ. 2) THEN
! Perfectly linear
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl = c(i, j)*rdx(i-1, j)
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = u(i-1, j) + (1.-cfl)*(br(i-1)-cfl*b0(i-1))
              CALL PUSHCONTROL(1,1)
            ELSE
              cfl = c(i, j)*rdx(i, j)
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = u(i, j) + (1.+cfl)*(bl(i)+cfl*b0(i))
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHINTEGER(ie3)
      CALL PUSHREALARRAY(b0, ie - is + 3)
      CALL PUSHREALARRAY(br, ie - is + 3)
      CALL PUSHREALARRAY(bl, ie - is + 3)
      CALL PUSHINTEGER(is3)
      CALL PUSHCONTROL(2,3)
    ELSE
      CALL PUSHCONTROL(2,2)
    END IF
  END SUBROUTINE XTP_U_FWD
!  Differentiation of xtp_u in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
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
!   gradient     of useful results: flux u c
!   with respect to varying inputs: flux u c
  SUBROUTINE XTP_U_BWD(is, ie, js, je, isd, ied, jsd, jed, c, c_ad, u&
&   , u_ad, v, flux, flux_ad, iord, dx, rdx, npx, npy, grid_type, nested&
& )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL :: u_ad(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL :: c_ad(is:ie+1, js:je+1)
    REAL :: flux(is:ie+1, js:je+1)
    REAL :: flux_ad(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: rdx(isd:ied, jsd:jed+1)
    INTEGER, INTENT(IN) :: iord, npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    REAL, DIMENSION(is-1:ie+1) :: bl_ad, br_ad, b0_ad
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: al(is-1:ie+2), dm(is-2:ie+2)
    REAL :: al_ad(is-1:ie+2)
    REAL :: dq(is-3:ie+2)
    REAL :: dl, dr, xt, pmp, lac, cfl
    REAL :: xt_ad, cfl_ad
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: x0, x1, x0l, x0r
    INTEGER :: i, j
    INTEGER :: is3, ie3
    INTEGER :: is2, ie2
    INTRINSIC MAX
    INTRINSIC MIN
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    INTEGER :: branch
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

    bl = 0.0
    br = 0.0
    b0 = 0.0
    fx0 = 0.0
    al = 0.0
    dm = 0.0
    dq = 0.0
    dl = 0.0
    dr = 0.0
    xt = 0.0
    pmp = 0.0
    lac = 0.0
    cfl = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    x0 = 0.0
    x1 = 0.0
    x0l = 0.0
    x0r = 0.0
    is3 = 0
    ie3 = 0
    is2 = 0
    ie2 = 0
    branch = 0

    CALL POPCONTROL(2,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY(flux(i, j))
              u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              CALL POPREALARRAY(flux(i, j))
              u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            END IF
          END DO
        END DO
      ELSE
        CALL POPREALARRAY(flux, (ie-is+2)*(je-js+2))
        DO j=js,je+1
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          DO i=ie+1,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              temp_ad10 = flux_ad(i, j)/6.0
              temp_ad11 = -(rdx(i, j)*0.5*flux_ad(i, j))
              temp_ad12 = c(i, j)*temp_ad11
              temp_ad13 = rdx(i, j)**2*flux_ad(i, j)
              temp_ad14 = c(i, j)**2*temp_ad13/6.0
              u_ad(i-1, j) = u_ad(i-1, j) + temp_ad14 - temp_ad12 + 2.0*&
&               temp_ad10
              u_ad(i, j) = u_ad(i, j) + temp_ad12 - 2.0*temp_ad14 + 5.0*&
&               temp_ad10
              u_ad(i+1, j) = u_ad(i+1, j) + temp_ad14 - temp_ad10
              c_ad(i, j) = c_ad(i, j) + (u(i+1, j)-2.0*u(i, j)+u(i-1, j)&
&               )*2*c(i, j)*temp_ad13/6.0 + (u(i, j)-u(i-1, j))*&
&               temp_ad11
              flux_ad(i, j) = 0.0
            ELSE
              temp_ad5 = flux_ad(i, j)/6.0
              temp_ad6 = -(rdx(i-1, j)*0.5*flux_ad(i, j))
              temp_ad7 = c(i, j)*temp_ad6
              temp_ad8 = rdx(i-1, j)**2*flux_ad(i, j)
              temp_ad9 = c(i, j)**2*temp_ad8/6.0
              u_ad(i, j) = u_ad(i, j) + temp_ad9 + temp_ad7 + 2.0*&
&               temp_ad5
              u_ad(i-1, j) = u_ad(i-1, j) + 5.0*temp_ad5 - temp_ad7 - &
&               2.0*temp_ad9
              u_ad(i-2, j) = u_ad(i-2, j) + temp_ad9 - temp_ad5
              c_ad(i, j) = c_ad(i, j) + (u(i, j)-2.0*u(i-1, j)+u(i-2, j)&
&               )*2*c(i, j)*temp_ad8/6.0 + (u(i, j)-u(i-1, j))*temp_ad6
              flux_ad(i, j) = 0.0
            END IF
          END DO
        END DO
      END IF
    ELSE IF (branch .NE. 2) THEN
      CALL POPINTEGER(is3)
      CALL POPREALARRAY(bl, ie - is + 3)
      CALL POPREALARRAY(br, ie - is + 3)
      CALL POPREALARRAY(b0, ie - is + 3)
      CALL POPINTEGER(ie3)
      al_ad = 0.0
      bl_ad = 0.0
      br_ad = 0.0
      b0_ad = 0.0
      DO j=je+1,js,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          DO i=ie+1,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              cfl = c(i, j)*rdx(i, j)
              CALL POPREALARRAY(flux(i, j))
              temp_ad4 = (cfl+1.)*flux_ad(i, j)
              u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
              cfl_ad = b0(i)*temp_ad4 + (bl(i)+cfl*b0(i))*flux_ad(i, j)
              bl_ad(i) = bl_ad(i) + temp_ad4
              b0_ad(i) = b0_ad(i) + cfl*temp_ad4
              flux_ad(i, j) = 0.0
              c_ad(i, j) = c_ad(i, j) + rdx(i, j)*cfl_ad
            ELSE
              cfl = c(i, j)*rdx(i-1, j)
              CALL POPREALARRAY(flux(i, j))
              temp_ad3 = (1.-cfl)*flux_ad(i, j)
              u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
              cfl_ad = -(b0(i-1)*temp_ad3) - (br(i-1)-cfl*b0(i-1))*&
&               flux_ad(i, j)
              br_ad(i-1) = br_ad(i-1) + temp_ad3
              b0_ad(i-1) = b0_ad(i-1) - cfl*temp_ad3
              flux_ad(i, j) = 0.0
              c_ad(i, j) = c_ad(i, j) + rdx(i-1, j)*cfl_ad
            END IF
          END DO
        END IF
        DO i=ie+1,is-1,-1
          CALL POPREALARRAY(b0(i))
          bl_ad(i) = bl_ad(i) + b0_ad(i)
          br_ad(i) = br_ad(i) + b0_ad(i)
          b0_ad(i) = 0.0
        END DO
        CALL POPCONTROL(2,branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) GOTO 100
        ELSE
          IF (branch .EQ. 2) THEN
            CALL POPREALARRAY(br(npx))
            u_ad(npx, j) = u_ad(npx, j) + (c3-1.0)*br_ad(npx)
            u_ad(npx+1, j) = u_ad(npx+1, j) + c2*br_ad(npx)
            u_ad(npx+2, j) = u_ad(npx+2, j) + c1*br_ad(npx)
            br_ad(npx) = 0.0
            CALL POPREALARRAY(bl(npx))
            xt_ad = br_ad(npx-1) + bl_ad(npx)
            u_ad(npx, j) = u_ad(npx, j) - bl_ad(npx)
            bl_ad(npx) = 0.0
            CALL POPREALARRAY(br(npx-1))
            temp_ad1 = 0.5*xt_ad/(dx(npx-1, j)+dx(npx-2, j))
            u_ad(npx-1, j) = u_ad(npx-1, j) + (dx(npx-1, j)*2.+dx(npx-2&
&             , j))*temp_ad1 - br_ad(npx-1)
            br_ad(npx-1) = 0.0
            temp_ad2 = 0.5*xt_ad/(dx(npx, j)+dx(npx+1, j))
            u_ad(npx-2, j) = u_ad(npx-2, j) - dx(npx-1, j)*temp_ad1
            u_ad(npx, j) = u_ad(npx, j) + (dx(npx, j)*2.+dx(npx+1, j))*&
&             temp_ad2
            u_ad(npx+1, j) = u_ad(npx+1, j) - dx(npx, j)*temp_ad2
          ELSE
            CALL POPREALARRAY(br(npx))
            br_ad(npx) = 0.0
            CALL POPREALARRAY(bl(npx))
            bl_ad(npx) = 0.0
            CALL POPREALARRAY(br(npx-1))
            br_ad(npx-1) = 0.0
            CALL POPREALARRAY(bl(npx-1))
            bl_ad(npx-1) = 0.0
          END IF
          CALL POPREALARRAY(bl(npx-1))
          xt_ad = br_ad(npx-2) + bl_ad(npx-1)
          u_ad(npx-1, j) = u_ad(npx-1, j) - bl_ad(npx-1)
          bl_ad(npx-1) = 0.0
          CALL POPREALARRAY(br(npx-2))
          u_ad(npx-2, j) = u_ad(npx-2, j) - br_ad(npx-2)
          br_ad(npx-2) = 0.0
          u_ad(npx-3, j) = u_ad(npx-3, j) + c1*xt_ad
          u_ad(npx-2, j) = u_ad(npx-2, j) + c2*xt_ad
          u_ad(npx-1, j) = u_ad(npx-1, j) + c3*xt_ad
          CALL POPREALARRAY(bl(npx-2))
          al_ad(npx-2) = al_ad(npx-2) + bl_ad(npx-2)
          u_ad(npx-2, j) = u_ad(npx-2, j) - bl_ad(npx-2)
          bl_ad(npx-2) = 0.0
        END IF
        CALL POPCONTROL(2,branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY(br(1))
          br_ad(1) = 0.0
          CALL POPREALARRAY(bl(1))
          bl_ad(1) = 0.0
          CALL POPREALARRAY(br(0))
          br_ad(0) = 0.0
          CALL POPREALARRAY(bl(0))
          bl_ad(0) = 0.0
        ELSE IF (branch .EQ. 1) THEN
          CALL POPREALARRAY(bl(1))
          xt_ad = br_ad(0) + bl_ad(1)
          u_ad(1, j) = u_ad(1, j) - bl_ad(1)
          bl_ad(1) = 0.0
          CALL POPREALARRAY(br(0))
          temp_ad = 0.5*xt_ad/(dx(0, j)+dx(-1, j))
          u_ad(0, j) = u_ad(0, j) + (dx(0, j)*2.+dx(-1, j))*temp_ad - &
&           br_ad(0)
          br_ad(0) = 0.0
          temp_ad0 = 0.5*xt_ad/(dx(1, j)+dx(2, j))
          u_ad(-1, j) = u_ad(-1, j) - dx(0, j)*temp_ad
          u_ad(1, j) = u_ad(1, j) + (dx(1, j)*2.+dx(2, j))*temp_ad0
          u_ad(2, j) = u_ad(2, j) - dx(1, j)*temp_ad0
          CALL POPREALARRAY(bl(0))
          u_ad(-2, j) = u_ad(-2, j) + c1*bl_ad(0)
          u_ad(-1, j) = u_ad(-1, j) + c2*bl_ad(0)
          u_ad(0, j) = u_ad(0, j) + (c3-1.0)*bl_ad(0)
          bl_ad(0) = 0.0
        ELSE
          GOTO 100
        END IF
        CALL POPREALARRAY(br(2))
        al_ad(3) = al_ad(3) + br_ad(2)
        u_ad(2, j) = u_ad(2, j) - bl_ad(2) - br_ad(2)
        br_ad(2) = 0.0
        CALL POPREALARRAY(bl(2))
        xt_ad = br_ad(1) + bl_ad(2)
        bl_ad(2) = 0.0
        CALL POPREALARRAY(br(1))
        u_ad(1, j) = u_ad(1, j) + c3*xt_ad - br_ad(1)
        br_ad(1) = 0.0
        u_ad(2, j) = u_ad(2, j) + c2*xt_ad
        u_ad(3, j) = u_ad(3, j) + c1*xt_ad
 100    DO i=ie3,is3,-1
          CALL POPREALARRAY(br(i))
          al_ad(i+1) = al_ad(i+1) + br_ad(i)
          u_ad(i, j) = u_ad(i, j) - bl_ad(i) - br_ad(i)
          br_ad(i) = 0.0
          CALL POPREALARRAY(bl(i))
          al_ad(i) = al_ad(i) + bl_ad(i)
          bl_ad(i) = 0.0
        END DO
        DO i=ie3+1,is3,-1
          u_ad(i-1, j) = u_ad(i-1, j) + p1*al_ad(i)
          u_ad(i, j) = u_ad(i, j) + p1*al_ad(i)
          u_ad(i-2, j) = u_ad(i-2, j) + p2*al_ad(i)
          u_ad(i+1, j) = u_ad(i+1, j) + p2*al_ad(i)
          al_ad(i) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
  END SUBROUTINE XTP_U_BWD
!  Differentiation of ytp_v in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
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
!   gradient     of useful results: flux v c
!   with respect to varying inputs: v c
  SUBROUTINE YTP_V_FWD(is, ie, js, je, isd, ied, jsd, jed, c, u, v, &
&   flux, jord, dy, rdy, npx, npy, grid_type, nested)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
!  Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL :: flux(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rdy(isd:ied+1, jsd:jed)
    INTEGER, INTENT(IN) :: npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
! Local:
    LOGICAL, DIMENSION(is:ie+1, js-1:je+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: dm(is:ie+1, js-2:je+2)
    REAL :: al(is:ie+1, js-1:je+2)
    REAL, DIMENSION(is:ie+1, js-1:je+1) :: bl, br, b0
    REAL :: dq(is:ie+1, js-3:je+2)
    REAL :: xt, dl, dr, pmp, lac, cfl
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: x0, x1, x0r, x0l
    INTEGER :: i, j, is1, ie1, js3, je3
    INTRINSIC MAX
    INTRINSIC MIN

    fx0 = 0.0
    dm = 0.0
    al = 0.0
    bl = 0.0
    br = 0.0
    b0 = 0.0
    dq = 0.0
    xt = 0.0
    dl = 0.0
    dr = 0.0
    pmp = 0.0
    lac = 0.0
    cfl = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    x0 = 0.0
    x1 = 0.0
    x0r = 0.0
    x0l = 0.0
    is1 = 0
    ie1 = 0
    js3 = 0
    je3 = 0

    IF (nested .OR. grid_type .GT. 3) THEN
      CALL PUSHCONTROL(1,0)
      js3 = js - 1
      je3 = je + 1
    ELSE
      IF (3 .LT. js - 1) THEN
        js3 = js - 1
      ELSE
        js3 = 3
      END IF
      IF (npy - 3 .GT. je + 1) THEN
        CALL PUSHCONTROL(1,1)
        je3 = je + 1
      ELSE
        CALL PUSHCONTROL(1,1)
        je3 = npy - 3
      END IF
    END IF
    IF (jord .EQ. 1) THEN
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux(i, j) = v(i, j-1)
            CALL PUSHCONTROL(1,1)
          ELSE
            flux(i, j) = v(i, j)
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      CALL PUSHCONTROL(3,0)
    ELSE IF (jord .EQ. 333) THEN
      CALL PUSHREALARRAY(c, (ie-is+2)*(je-js+2))
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux(i, j) = (2.0*v(i, j)+5.0*v(i, j-1)-v(i, j-2))/6.0 - 0.5&
&             *c(i, j)*rdy(i, j-1)*(v(i, j)-v(i, j-1)) + c(i, j)*rdy(i, &
&             j-1)*c(i, j)*rdy(i, j-1)/6.0*(v(i, j)-2.0*v(i, j-1)+v(i, j&
&             -2))
          ELSE
            flux(i, j) = (2.0*v(i, j-1)+5.0*v(i, j)-v(i, j+1))/6.0 - 0.5&
&             *c(i, j)*rdy(i, j)*(v(i, j)-v(i, j-1)) + c(i, j)*rdy(i, j)&
&             *c(i, j)*rdy(i, j)/6.0*(v(i, j+1)-2.0*v(i, j)+v(i, j-1))
          END IF
        END DO
      END DO
      CALL PUSHCONTROL(3,1)
    ELSE IF (jord .LT. 8) THEN
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
      DO j=js3,je3+1
        DO i=is,ie+1
          al(i, j) = p1*(v(i, j-1)+v(i, j)) + p2*(v(i, j-2)+v(i, j+1))
        END DO
      END DO
      DO j=js3,je3
        DO i=is,ie+1
          bl(i, j) = al(i, j) - v(i, j)
          br(i, j) = al(i, j+1) - v(i, j)
        END DO
      END DO
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
            bl(i, 0) = c1*v(i, -2) + c2*v(i, -1) + c3*v(i, 0) - v(i, 0)
            xt = 0.5*(((2.*dy(i, 0)+dy(i, -1))*v(i, 0)-dy(i, 0)*v(i, -1)&
&             )/(dy(i, 0)+dy(i, -1))+((2.*dy(i, 1)+dy(i, 2))*v(i, 1)-dy(&
&             i, 1)*v(i, 2))/(dy(i, 1)+dy(i, 2)))
            br(i, 0) = xt - v(i, 0)
            bl(i, 1) = xt - v(i, 1)
            xt = c3*v(i, 1) + c2*v(i, 2) + c1*v(i, 3)
            br(i, 1) = xt - v(i, 1)
            bl(i, 2) = xt - v(i, 2)
            br(i, 2) = al(i, 3) - v(i, 2)
          END DO
          IF (is .EQ. 1) THEN
! out
            bl(1, 0) = 0.
! edge
            br(1, 0) = 0.
! edge
            bl(1, 1) = 0.
! in
            br(1, 1) = 0.
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
! out
            bl(npx, 0) = 0.
! edge
            br(npx, 0) = 0.
! edge
            bl(npx, 1) = 0.
! in
            br(npx, 1) = 0.
            CALL PUSHCONTROL(2,0)
          ELSE
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
!      j=2
!      call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
            bl(i, npy-2) = al(i, npy-2) - v(i, npy-2)
            xt = c1*v(i, npy-3) + c2*v(i, npy-2) + c3*v(i, npy-1)
            br(i, npy-2) = xt - v(i, npy-2)
            bl(i, npy-1) = xt - v(i, npy-1)
            xt = 0.5*(((2.*dy(i, npy-1)+dy(i, npy-2))*v(i, npy-1)-dy(i, &
&             npy-1)*v(i, npy-2))/(dy(i, npy-1)+dy(i, npy-2))+((2.*dy(i&
&             , npy)+dy(i, npy+1))*v(i, npy)-dy(i, npy)*v(i, npy+1))/(dy&
&             (i, npy)+dy(i, npy+1)))
            br(i, npy-1) = xt - v(i, npy-1)
            bl(i, npy) = xt - v(i, npy)
            br(i, npy) = c3*v(i, npy) + c2*v(i, npy+1) + c1*v(i, npy+2) &
&             - v(i, npy)
          END DO
          IF (is .EQ. 1) THEN
! in
            bl(1, npy-1) = 0.
! edge
            br(1, npy-1) = 0.
! edge
            bl(1, npy) = 0.
! out
            br(1, npy) = 0.
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
! in
            bl(npx, npy-1) = 0.
! edge
            br(npx, npy-1) = 0.
! edge
            bl(npx, npy) = 0.
! out
            br(npx, npy) = 0.
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
!      j=npy-2
!      call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
      DO j=js-1,je+1
        DO i=is,ie+1
          b0(i, j) = bl(i, j) + br(i, j)
        END DO
      END DO
      IF (jord .EQ. 2) THEN
! Perfectly linear
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHREALARRAY(cfl)
              cfl = c(i, j)*rdy(i, j-1)
              flux(i, j) = v(i, j-1) + (1.-cfl)*(br(i, j-1)-cfl*b0(i, j-&
&               1))
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHREALARRAY(cfl)
              cfl = c(i, j)*rdy(i, j)
              flux(i, j) = v(i, j) + (1.+cfl)*(bl(i, j)+cfl*b0(i, j))
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
        END DO
        CALL PUSHREALARRAY(b0, (ie-is+2)*(je-js+3))
        CALL PUSHINTEGER(js3)
        CALL PUSHREALARRAY(br, (ie-is+2)*(je-js+3))
        CALL PUSHINTEGER(je3)
        CALL PUSHREALARRAY(bl, (ie-is+2)*(je-js+3))
        CALL PUSHREALARRAY(cfl)
        CALL PUSHCONTROL(3,4)
      ELSE
        CALL PUSHINTEGER(js3)
        CALL PUSHINTEGER(je3)
        CALL PUSHCONTROL(3,3)
      END IF
    ELSE
      CALL PUSHCONTROL(3,2)
    END IF
  END SUBROUTINE YTP_V_FWD
!  Differentiation of ytp_v in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
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
!   gradient     of useful results: flux v c
!   with respect to varying inputs: v c
  SUBROUTINE YTP_V_BWD(is, ie, js, je, isd, ied, jsd, jed, c, c_ad, u&
&   , v, v_ad, flux, flux_ad, jord, dy, rdy, npx, npy, grid_type, nested&
& )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL :: v_ad(isd:ied+1, jsd:jed)
    REAL, INTENT(INOUT) :: c(is:ie+1, js:je+1)
    REAL :: c_ad(is:ie+1, js:je+1)
    REAL :: flux(is:ie+1, js:je+1)
    REAL :: flux_ad(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rdy(isd:ied+1, jsd:jed)
    INTEGER, INTENT(IN) :: npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
    LOGICAL, DIMENSION(is:ie+1, js-1:je+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: dm(is:ie+1, js-2:je+2)
    REAL :: al(is:ie+1, js-1:je+2)
    REAL :: al_ad(is:ie+1, js-1:je+2)
    REAL, DIMENSION(is:ie+1, js-1:je+1) :: bl, br, b0
    REAL, DIMENSION(is:ie+1, js-1:je+1) :: bl_ad, br_ad, b0_ad
    REAL :: dq(is:ie+1, js-3:je+2)
    REAL :: xt, dl, dr, pmp, lac, cfl
    REAL :: xt_ad, cfl_ad
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: x0, x1, x0r, x0l
    INTEGER :: i, j, is1, ie1, js3, je3
    INTRINSIC MAX
    INTRINSIC MIN
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    INTEGER :: branch
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

    fx0 = 0.0
    dm = 0.0
    al = 0.0
    bl = 0.0
    br = 0.0
    b0 = 0.0
    dq = 0.0
    xt = 0.0
    dl = 0.0
    dr = 0.0
    pmp = 0.0
    lac = 0.0
    cfl = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    x0 = 0.0
    x1 = 0.0
    x0r = 0.0
    x0l = 0.0
    is1 = 0
    ie1 = 0
    js3 = 0
    je3 = 0
    branch = 0

    CALL POPCONTROL(3,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            END IF
          END DO
        END DO
      ELSE
        CALL POPREALARRAY(c, (ie-is+2)*(je-js+2))
        DO j=js,je+1
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          DO i=ie+1,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              temp_ad10 = flux_ad(i, j)/6.0
              temp_ad11 = -(rdy(i, j)*0.5*flux_ad(i, j))
              temp_ad12 = c(i, j)*temp_ad11
              temp_ad13 = rdy(i, j)**2*flux_ad(i, j)
              temp_ad14 = c(i, j)**2*temp_ad13/6.0
              v_ad(i, j-1) = v_ad(i, j-1) + temp_ad14 - temp_ad12 + 2.0*&
&               temp_ad10
              v_ad(i, j) = v_ad(i, j) + temp_ad12 - 2.0*temp_ad14 + 5.0*&
&               temp_ad10
              v_ad(i, j+1) = v_ad(i, j+1) + temp_ad14 - temp_ad10
              c_ad(i, j) = c_ad(i, j) + (v(i, j+1)-2.0*v(i, j)+v(i, j-1)&
&               )*2*c(i, j)*temp_ad13/6.0 + (v(i, j)-v(i, j-1))*&
&               temp_ad11
              flux_ad(i, j) = 0.0
            ELSE
              temp_ad5 = flux_ad(i, j)/6.0
              temp_ad6 = -(rdy(i, j-1)*0.5*flux_ad(i, j))
              temp_ad7 = c(i, j)*temp_ad6
              temp_ad8 = rdy(i, j-1)**2*flux_ad(i, j)
              temp_ad9 = c(i, j)**2*temp_ad8/6.0
              v_ad(i, j) = v_ad(i, j) + temp_ad9 + temp_ad7 + 2.0*&
&               temp_ad5
              v_ad(i, j-1) = v_ad(i, j-1) + 5.0*temp_ad5 - temp_ad7 - &
&               2.0*temp_ad9
              v_ad(i, j-2) = v_ad(i, j-2) + temp_ad9 - temp_ad5
              c_ad(i, j) = c_ad(i, j) + (v(i, j)-2.0*v(i, j-1)+v(i, j-2)&
&               )*2*c(i, j)*temp_ad8/6.0 + (v(i, j)-v(i, j-1))*temp_ad6
              flux_ad(i, j) = 0.0
            END IF
          END DO
        END DO
      END IF
    ELSE IF (branch .NE. 2) THEN
      IF (branch .EQ. 3) THEN
        CALL POPINTEGER(je3)
        CALL POPINTEGER(js3)
        bl_ad = 0.0
        br_ad = 0.0
        b0_ad = 0.0
      ELSE
        CALL POPREALARRAY(cfl)
        CALL POPREALARRAY(bl, (ie-is+2)*(je-js+3))
        CALL POPINTEGER(je3)
        CALL POPREALARRAY(br, (ie-is+2)*(je-js+3))
        CALL POPINTEGER(js3)
        CALL POPREALARRAY(b0, (ie-is+2)*(je-js+3))
        bl_ad = 0.0
        br_ad = 0.0
        b0_ad = 0.0
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              temp_ad4 = (cfl+1.)*flux_ad(i, j)
              v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
              cfl_ad = b0(i, j)*temp_ad4 + (bl(i, j)+cfl*b0(i, j))*&
&               flux_ad(i, j)
              bl_ad(i, j) = bl_ad(i, j) + temp_ad4
              b0_ad(i, j) = b0_ad(i, j) + cfl*temp_ad4
              flux_ad(i, j) = 0.0
              CALL POPREALARRAY(cfl)
              c_ad(i, j) = c_ad(i, j) + rdy(i, j)*cfl_ad
            ELSE
              temp_ad3 = (1.-cfl)*flux_ad(i, j)
              v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
              cfl_ad = -(b0(i, j-1)*temp_ad3) - (br(i, j-1)-cfl*b0(i, j-&
&               1))*flux_ad(i, j)
              br_ad(i, j-1) = br_ad(i, j-1) + temp_ad3
              b0_ad(i, j-1) = b0_ad(i, j-1) - cfl*temp_ad3
              flux_ad(i, j) = 0.0
              CALL POPREALARRAY(cfl)
              c_ad(i, j) = c_ad(i, j) + rdy(i, j-1)*cfl_ad
            END IF
          END DO
        END DO
      END IF
      DO j=je+1,js-1,-1
        DO i=ie+1,is,-1
          bl_ad(i, j) = bl_ad(i, j) + b0_ad(i, j)
          br_ad(i, j) = br_ad(i, j) + b0_ad(i, j)
          b0_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL(2,branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          al_ad = 0.0
          GOTO 100
        ELSE
          al_ad = 0.0
        END IF
      ELSE
        IF (branch .NE. 2) THEN
          br_ad(npx, npy) = 0.0
          bl_ad(npx, npy) = 0.0
          br_ad(npx, npy-1) = 0.0
          bl_ad(npx, npy-1) = 0.0
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          br_ad(1, npy) = 0.0
          bl_ad(1, npy) = 0.0
          br_ad(1, npy-1) = 0.0
          bl_ad(1, npy-1) = 0.0
        END IF
        al_ad = 0.0
        DO i=ie+1,is,-1
          v_ad(i, npy) = v_ad(i, npy) + (c3-1.0)*br_ad(i, npy)
          v_ad(i, npy+1) = v_ad(i, npy+1) + c2*br_ad(i, npy)
          v_ad(i, npy+2) = v_ad(i, npy+2) + c1*br_ad(i, npy)
          br_ad(i, npy) = 0.0
          xt_ad = br_ad(i, npy-1) + bl_ad(i, npy)
          v_ad(i, npy) = v_ad(i, npy) - bl_ad(i, npy)
          bl_ad(i, npy) = 0.0
          temp_ad1 = 0.5*xt_ad/(dy(i, npy-1)+dy(i, npy-2))
          v_ad(i, npy-1) = v_ad(i, npy-1) + (dy(i, npy-1)*2.+dy(i, npy-2&
&           ))*temp_ad1 - br_ad(i, npy-1)
          br_ad(i, npy-1) = 0.0
          temp_ad2 = 0.5*xt_ad/(dy(i, npy)+dy(i, npy+1))
          v_ad(i, npy-2) = v_ad(i, npy-2) - dy(i, npy-1)*temp_ad1
          v_ad(i, npy) = v_ad(i, npy) + (dy(i, npy)*2.+dy(i, npy+1))*&
&           temp_ad2
          v_ad(i, npy+1) = v_ad(i, npy+1) - dy(i, npy)*temp_ad2
          xt_ad = br_ad(i, npy-2) + bl_ad(i, npy-1)
          v_ad(i, npy-1) = v_ad(i, npy-1) - bl_ad(i, npy-1)
          bl_ad(i, npy-1) = 0.0
          v_ad(i, npy-2) = v_ad(i, npy-2) - br_ad(i, npy-2)
          br_ad(i, npy-2) = 0.0
          v_ad(i, npy-3) = v_ad(i, npy-3) + c1*xt_ad
          v_ad(i, npy-2) = v_ad(i, npy-2) + c2*xt_ad
          v_ad(i, npy-1) = v_ad(i, npy-1) + c3*xt_ad
          al_ad(i, npy-2) = al_ad(i, npy-2) + bl_ad(i, npy-2)
          v_ad(i, npy-2) = v_ad(i, npy-2) - bl_ad(i, npy-2)
          bl_ad(i, npy-2) = 0.0
        END DO
      END IF
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        br_ad(npx, 1) = 0.0
        bl_ad(npx, 1) = 0.0
        br_ad(npx, 0) = 0.0
        bl_ad(npx, 0) = 0.0
      ELSE IF (branch .NE. 1) THEN
        GOTO 100
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        br_ad(1, 1) = 0.0
        bl_ad(1, 1) = 0.0
        br_ad(1, 0) = 0.0
        bl_ad(1, 0) = 0.0
      END IF
      DO i=ie+1,is,-1
        al_ad(i, 3) = al_ad(i, 3) + br_ad(i, 2)
        v_ad(i, 2) = v_ad(i, 2) - bl_ad(i, 2) - br_ad(i, 2)
        br_ad(i, 2) = 0.0
        xt_ad = br_ad(i, 1) + bl_ad(i, 2)
        bl_ad(i, 2) = 0.0
        v_ad(i, 1) = v_ad(i, 1) + c3*xt_ad - br_ad(i, 1)
        br_ad(i, 1) = 0.0
        v_ad(i, 2) = v_ad(i, 2) + c2*xt_ad
        v_ad(i, 3) = v_ad(i, 3) + c1*xt_ad
        xt_ad = br_ad(i, 0) + bl_ad(i, 1)
        v_ad(i, 1) = v_ad(i, 1) - bl_ad(i, 1)
        bl_ad(i, 1) = 0.0
        temp_ad = 0.5*xt_ad/(dy(i, 0)+dy(i, -1))
        v_ad(i, 0) = v_ad(i, 0) + (dy(i, 0)*2.+dy(i, -1))*temp_ad - &
&         br_ad(i, 0)
        br_ad(i, 0) = 0.0
        temp_ad0 = 0.5*xt_ad/(dy(i, 1)+dy(i, 2))
        v_ad(i, -1) = v_ad(i, -1) - dy(i, 0)*temp_ad
        v_ad(i, 1) = v_ad(i, 1) + (dy(i, 1)*2.+dy(i, 2))*temp_ad0
        v_ad(i, 2) = v_ad(i, 2) - dy(i, 1)*temp_ad0
        v_ad(i, -2) = v_ad(i, -2) + c1*bl_ad(i, 0)
        v_ad(i, -1) = v_ad(i, -1) + c2*bl_ad(i, 0)
        v_ad(i, 0) = v_ad(i, 0) + (c3-1.0)*bl_ad(i, 0)
        bl_ad(i, 0) = 0.0
      END DO
 100  DO j=je3,js3,-1
        DO i=ie+1,is,-1
          al_ad(i, j+1) = al_ad(i, j+1) + br_ad(i, j)
          v_ad(i, j) = v_ad(i, j) - bl_ad(i, j) - br_ad(i, j)
          br_ad(i, j) = 0.0
          al_ad(i, j) = al_ad(i, j) + bl_ad(i, j)
          bl_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je3+1,js3,-1
        DO i=ie+1,is,-1
          v_ad(i, j-1) = v_ad(i, j-1) + p1*al_ad(i, j)
          v_ad(i, j) = v_ad(i, j) + p1*al_ad(i, j)
          v_ad(i, j-2) = v_ad(i, j-2) + p2*al_ad(i, j)
          v_ad(i, j+1) = v_ad(i, j+1) + p2*al_ad(i, j)
          al_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
  END SUBROUTINE YTP_V_BWD
!  Differentiation of compute_divergence_damping in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: u v ke ua uc ptc delpc va vc
!                vort divg_d wk
!   with respect to varying inputs: u v ke ua uc ptc delpc va vc
!                divg_d wk
  SUBROUTINE COMPUTE_DIVERGENCE_DAMPING_ADM(nord, d2_bg, d4_bg, dddmp, &
&   dt, vort, vort_ad, ptc, ptc_ad, delpc, delpc_ad, ke, ke_ad, u, u_ad&
&   , v, v_ad, uc, uc_ad, vc, vc_ad, ua, ua_ad, va, va_ad, divg_d, &
&   divg_d_ad, wk, wk_ad, gridstruct, flagstruct, bd)
    IMPLICIT NONE
!InOut Arguments
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    INTEGER, INTENT(IN) :: nord
    REAL, INTENT(IN) :: d2_bg, d4_bg, dddmp, dt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ua_ad, va_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1) :: u_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed) :: v_ad
!Intent is really in
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: wk
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   wk_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: vort
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   vort_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delpc, ptc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delpc_ad, ptc_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   ke
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   ke_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   vc_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   uc_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   divg_d
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   divg_d_ad
!Locals
    REAL :: damp, dd8, damp2, da_min, da_min_c, absdt
    REAL :: damp_ad, damp2_ad
    INTEGER :: is, ie, js, je, npx, npy, is2, ie1
    LOGICAL :: nested, fill_c
    INTEGER :: i, j, n, n2, nt
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, DIMENSION(:, :), POINTER :: area, area_c, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsina
    REAL, DIMENSION(:, :), POINTER :: f0, rsin2, divg_u, divg_v
    REAL, DIMENSION(:, :), POINTER :: cosa, dx, dy, dxc, dyc, rdxa, rdya&
&   , rdx, rdy
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SQRT
    REAL :: max1
    REAL :: max1_ad
    REAL :: abs0
    REAL :: abs1
    REAL :: max2
    REAL :: max2_ad
    REAL :: abs2
    REAL :: abs2_ad
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: y3_ad
    REAL :: y1_ad
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: y2_ad
    INTEGER :: branch
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
    INTEGER :: ad_from5
    INTEGER :: ad_to5
    INTEGER :: ad_from6
    INTEGER :: ad_to6
    REAL :: y3
    REAL :: y2
    REAL :: y1
    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    sina_u => gridstruct%sina_u
    sina_v => gridstruct%sina_v
    divg_u => gridstruct%divg_u
    divg_v => gridstruct%divg_v
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    sw_corner = gridstruct%sw_corner
    se_corner = gridstruct%se_corner
    nw_corner = gridstruct%nw_corner
    ne_corner = gridstruct%ne_corner
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    IF (nested) THEN
      CALL PUSHCONTROL1B(0)
      is2 = is
      ie1 = ie + 1
    ELSE
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        CALL PUSHCONTROL1B(1)
        ie1 = ie + 1
      ELSE
        CALL PUSHCONTROL1B(1)
        ie1 = npx - 1
      END IF
    END IF
!-----------------------------
! Compute divergence damping
!-----------------------------
!  damp = dddmp * da_min_c
    IF (nord .EQ. 0) THEN
!         area ~ dxb*dyb*sin(alpha)
      IF (nested) THEN
        DO j=js,je+1
          DO i=is-1,ie+1
            CALL PUSHREALARRAY_ADM(ptc(i, j))
            ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j))&
&             *dyc(i, j)*sina_v(i, j)
          END DO
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is-1,ie+1
              IF (vc(i, j) .GT. 0) THEN
                CALL PUSHREALARRAY_ADM(ptc(i, j))
                ptc(i, j) = u(i, j)*dyc(i, j)*sin_sg(i, j-1, 4)
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHREALARRAY_ADM(ptc(i, j))
                ptc(i, j) = u(i, j)*dyc(i, j)*sin_sg(i, j, 2)
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
            CALL PUSHCONTROL1B(1)
          ELSE
            DO i=is-1,ie+1
              CALL PUSHREALARRAY_ADM(ptc(i, j))
              ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j&
&               ))*dyc(i, j)*sina_v(i, j)
            END DO
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
          IF (is .EQ. 1) THEN
            IF (uc(1, j) .GT. 0) THEN
              vort(1, j) = v(1, j)*dxc(1, j)*sin_sg(0, j, 3)
              CALL PUSHCONTROL2B(0)
            ELSE
              vort(1, j) = v(1, j)*dxc(1, j)*sin_sg(1, j, 1)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE
            CALL PUSHCONTROL2B(2)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            IF (uc(npx, j) .GT. 0) THEN
              vort(npx, j) = v(npx, j)*dxc(npx, j)*sin_sg(npx-1, j, 3)
              CALL PUSHCONTROL2B(2)
            ELSE
              vort(npx, j) = v(npx, j)*dxc(npx, j)*sin_sg(npx, j, 1)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE
            CALL PUSHCONTROL2B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY_ADM(delpc(i, j))
          delpc(i, j) = vort(i, j-1) - vort(i, j) + ptc(i-1, j) - ptc(i&
&           , j)
        END DO
      END DO
! Remove the extra term at the corners:
      IF (sw_corner) THEN
        CALL PUSHREALARRAY_ADM(delpc(1, 1))
        delpc(1, 1) = delpc(1, 1) - vort(1, 0)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (se_corner) THEN
        CALL PUSHREALARRAY_ADM(delpc(npx, 1))
        delpc(npx, 1) = delpc(npx, 1) - vort(npx, 0)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHREALARRAY_ADM(delpc(npx, npy))
        delpc(npx, npy) = delpc(npx, npy) + vort(npx, npy)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (nw_corner) THEN
        CALL PUSHREALARRAY_ADM(delpc(1, npy))
        delpc(1, npy) = delpc(1, npy) + vort(1, npy)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY_ADM(delpc(i, j))
          delpc(i, j) = gridstruct%rarea_c(i, j)*delpc(i, j)
          IF (delpc(i, j)*dt .GE. 0.) THEN
            abs2 = delpc(i, j)*dt
            CALL PUSHCONTROL1B(0)
          ELSE
            abs2 = -(delpc(i, j)*dt)
            CALL PUSHCONTROL1B(1)
          END IF
          y3 = dddmp*abs2
          IF (0.20 .GT. y3) THEN
            y1 = y3
            CALL PUSHCONTROL1B(0)
          ELSE
            y1 = 0.20
            CALL PUSHCONTROL1B(1)
          END IF
          IF (d2_bg .LT. y1) THEN
            max1 = y1
            CALL PUSHCONTROL1B(0)
          ELSE
            max1 = d2_bg
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREALARRAY_ADM(damp)
          damp = gridstruct%da_min_c*max1
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          vort_ad(i, j) = vort_ad(i, j) + ke_ad(i, j)
          damp_ad = delpc(i, j)*vort_ad(i, j)
          delpc_ad(i, j) = delpc_ad(i, j) + damp*vort_ad(i, j)
          vort_ad(i, j) = 0.0
          CALL POPREALARRAY_ADM(damp)
          max1_ad = gridstruct%da_min_c*damp_ad
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y1_ad = max1_ad
          ELSE
            y1_ad = 0.0
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y3_ad = y1_ad
          ELSE
            y3_ad = 0.0
          END IF
          abs2_ad = dddmp*y3_ad
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            delpc_ad(i, j) = delpc_ad(i, j) + dt*abs2_ad
          ELSE
            delpc_ad(i, j) = delpc_ad(i, j) - dt*abs2_ad
          END IF
          CALL POPREALARRAY_ADM(delpc(i, j))
          delpc_ad(i, j) = gridstruct%rarea_c(i, j)*delpc_ad(i, j)
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        CALL POPREALARRAY_ADM(delpc(1, npy))
        vort_ad(1, npy) = vort_ad(1, npy) + delpc_ad(1, npy)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY_ADM(delpc(npx, npy))
        vort_ad(npx, npy) = vort_ad(npx, npy) + delpc_ad(npx, npy)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY_ADM(delpc(npx, 1))
        vort_ad(npx, 0) = vort_ad(npx, 0) - delpc_ad(npx, 1)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY_ADM(delpc(1, 1))
        vort_ad(1, 0) = vort_ad(1, 0) - delpc_ad(1, 1)
      END IF
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY_ADM(delpc(i, j))
          vort_ad(i, j-1) = vort_ad(i, j-1) + delpc_ad(i, j)
          vort_ad(i, j) = vort_ad(i, j) - delpc_ad(i, j)
          ptc_ad(i-1, j) = ptc_ad(i-1, j) + delpc_ad(i, j)
          ptc_ad(i, j) = ptc_ad(i, j) - delpc_ad(i, j)
          delpc_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je+1,js-1,-1
          CALL POPCONTROL2B(branch)
          IF (branch .NE. 0) THEN
            IF (branch .EQ. 1) THEN
              v_ad(npx, j) = v_ad(npx, j) + sin_sg(npx, j, 1)*dxc(npx, j&
&               )*vort_ad(npx, j)
              vort_ad(npx, j) = 0.0
            ELSE
              v_ad(npx, j) = v_ad(npx, j) + sin_sg(npx-1, j, 3)*dxc(npx&
&               , j)*vort_ad(npx, j)
              vort_ad(npx, j) = 0.0
            END IF
          END IF
          CALL POPCONTROL2B(branch)
          IF (branch .EQ. 0) THEN
            v_ad(1, j) = v_ad(1, j) + sin_sg(0, j, 3)*dxc(1, j)*vort_ad(&
&             1, j)
            vort_ad(1, j) = 0.0
          ELSE IF (branch .EQ. 1) THEN
            v_ad(1, j) = v_ad(1, j) + sin_sg(1, j, 1)*dxc(1, j)*vort_ad(&
&             1, j)
            vort_ad(1, j) = 0.0
          END IF
          DO i=ie1,is2,-1
            temp_ad5 = dxc(i, j)*sina_u(i, j)*vort_ad(i, j)
            temp_ad6 = -(cosa_u(i, j)*0.5*temp_ad5)
            v_ad(i, j) = v_ad(i, j) + temp_ad5
            ua_ad(i-1, j) = ua_ad(i-1, j) + temp_ad6
            ua_ad(i, j) = ua_ad(i, j) + temp_ad6
            vort_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je+1,js,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            DO i=ie+1,is-1,-1
              CALL POPREALARRAY_ADM(ptc(i, j))
              temp_ad3 = dyc(i, j)*sina_v(i, j)*ptc_ad(i, j)
              temp_ad4 = -(cosa_v(i, j)*0.5*temp_ad3)
              u_ad(i, j) = u_ad(i, j) + temp_ad3
              va_ad(i, j-1) = va_ad(i, j-1) + temp_ad4
              va_ad(i, j) = va_ad(i, j) + temp_ad4
              ptc_ad(i, j) = 0.0
            END DO
          ELSE
            DO i=ie+1,is-1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(ptc(i, j))
                u_ad(i, j) = u_ad(i, j) + sin_sg(i, j, 2)*dyc(i, j)*&
&                 ptc_ad(i, j)
                ptc_ad(i, j) = 0.0
              ELSE
                CALL POPREALARRAY_ADM(ptc(i, j))
                u_ad(i, j) = u_ad(i, j) + sin_sg(i, j-1, 4)*dyc(i, j)*&
&                 ptc_ad(i, j)
                ptc_ad(i, j) = 0.0
              END IF
            END DO
          END IF
        END DO
      ELSE
        DO j=je+1,js-1,-1
          DO i=ie1,is2,-1
            temp_ad1 = dxc(i, j)*sina_u(i, j)*vort_ad(i, j)
            temp_ad2 = -(cosa_u(i, j)*0.5*temp_ad1)
            v_ad(i, j) = v_ad(i, j) + temp_ad1
            ua_ad(i-1, j) = ua_ad(i-1, j) + temp_ad2
            ua_ad(i, j) = ua_ad(i, j) + temp_ad2
            vort_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ie+1,is-1,-1
            CALL POPREALARRAY_ADM(ptc(i, j))
            temp_ad = dyc(i, j)*sina_v(i, j)*ptc_ad(i, j)
            temp_ad0 = -(cosa_v(i, j)*0.5*temp_ad)
            u_ad(i, j) = u_ad(i, j) + temp_ad
            va_ad(i, j-1) = va_ad(i, j-1) + temp_ad0
            va_ad(i, j) = va_ad(i, j) + temp_ad0
            ptc_ad(i, j) = 0.0
          END DO
        END DO
      END IF
    ELSE
!--------------------------
! Higher order divg damping
!--------------------------
      DO j=js,je+1
        DO i=is,ie+1
! Save divergence for external mode filter
          CALL PUSHREALARRAY_ADM(delpc(i, j))
          delpc(i, j) = divg_d(i, j)
        END DO
      END DO
! N > 1
      n2 = nord + 1
      DO n=1,nord
        nt = nord - n
        fill_c = nt .NE. 0 .AND. flagstruct%grid_type .LT. 3 .AND. (((&
&         sw_corner .OR. se_corner) .OR. ne_corner) .OR. nw_corner) &
&         .AND. (.NOT.nested)
        IF (fill_c) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        ad_from0 = js - nt
        DO j=ad_from0,je+1+nt
          ad_from = is - 1 - nt
          i = ie + nt + 2
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from0)
        IF (fill_c) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        ad_from2 = js - 1 - nt
        DO j=ad_from2,je+1+nt
          ad_from1 = is - nt
          i = ie + nt + 2
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from1)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from2)
        IF (fill_c) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        ad_from4 = js - nt
        DO j=ad_from4,je+1+nt
          ad_from3 = is - nt
          i = ie + nt + 2
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from3)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from4)
! Remove the extra term at the corners:
        IF (sw_corner) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (se_corner) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (ne_corner) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (nw_corner) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (.NOT.gridstruct%stretched_grid) THEN
          ad_from6 = js - nt
          DO j=ad_from6,je+1+nt
            ad_from5 = is - nt
            i = ie + nt + 2
            CALL PUSHINTEGER4(i - 1)
            CALL PUSHINTEGER4(ad_from5)
          END DO
          CALL PUSHINTEGER4(j - 1)
          CALL PUSHINTEGER4(ad_from6)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
! n-loop
      IF (dddmp .LT. 1.e-5) THEN
        CALL PUSHCONTROL2B(0)
        vort(:, :) = 0.
      ELSE IF (flagstruct%grid_type .LT. 3) THEN
! Interpolate relative vort to cell corners
        CALL PUSHREALARRAY_ADM(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL A2B_ORD4(wk, vort, gridstruct, npx, npy, is, ie, js, je, ng&
&               , .false.)
        DO j=js,je+1
          DO i=is,ie+1
            IF (dt .GE. 0.) THEN
              CALL PUSHREALARRAY_ADM(abs0)
              abs0 = dt
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(abs0)
              abs0 = -dt
              CALL PUSHCONTROL1B(1)
            END IF
! The following is an approxi form of Smagorinsky diffusion
            CALL PUSHREALARRAY_ADM(vort(i, j))
            vort(i, j) = abs0*SQRT(delpc(i, j)**2+vort(i, j)**2)
          END DO
        END DO
        CALL PUSHCONTROL2B(1)
      ELSE
        IF (dt .GE. 0.) THEN
          abs1 = dt
        ELSE
          abs1 = -dt
        END IF
! Correct form: works only for doubly preiodic domain
        CALL PUSHREALARRAY_ADM(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
        CALL SMAG_CORNER(abs1, u, v, ua, va, vort, bd, npx, npy, &
&                  gridstruct, ng)
        CALL PUSHCONTROL2B(2)
      END IF
      IF (gridstruct%stretched_grid) THEN
! Stretched grid with variable damping ~ area
        dd8 = gridstruct%da_min*d4_bg**n2
      ELSE
        dd8 = (gridstruct%da_min_c*d4_bg)**n2
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          IF (0.20 .GT. dddmp*vort(i, j)) THEN
            y2 = dddmp*vort(i, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
            y2 = 0.20
          END IF
          IF (d2_bg .LT. y2) THEN
            max2 = y2
            CALL PUSHCONTROL1B(0)
          ELSE
            max2 = d2_bg
            CALL PUSHCONTROL1B(1)
          END IF
! del-2
          CALL PUSHREALARRAY_ADM(damp2)
          damp2 = gridstruct%da_min_c*max2
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          vort_ad(i, j) = vort_ad(i, j) + ke_ad(i, j)
          damp2_ad = delpc(i, j)*vort_ad(i, j)
          delpc_ad(i, j) = delpc_ad(i, j) + damp2*vort_ad(i, j)
          divg_d_ad(i, j) = divg_d_ad(i, j) + dd8*vort_ad(i, j)
          vort_ad(i, j) = 0.0
          CALL POPREALARRAY_ADM(damp2)
          max2_ad = gridstruct%da_min_c*damp2_ad
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y2_ad = max2_ad
          ELSE
            y2_ad = 0.0
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) vort_ad(i, j) = vort_ad(i, j) + dddmp*y2_ad
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .NE. 0) THEN
        IF (branch .EQ. 1) THEN
          DO j=je+1,js,-1
            DO i=ie+1,is,-1
              CALL POPREALARRAY_ADM(vort(i, j))
              IF (delpc(i, j)**2 + vort(i, j)**2 .EQ. 0.0) THEN
                temp_ad9 = 0.0
              ELSE
                temp_ad9 = abs0*vort_ad(i, j)/(2.0*SQRT(delpc(i, j)**2+&
&                 vort(i, j)**2))
              END IF
              delpc_ad(i, j) = delpc_ad(i, j) + 2*delpc(i, j)*temp_ad9
              vort_ad(i, j) = 2*vort(i, j)*temp_ad9
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(abs0)
              ELSE
                CALL POPREALARRAY_ADM(abs0)
              END IF
            END DO
          END DO
          CALL POPREALARRAY_ADM(wk, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
          CALL A2B_ORD4_ADM(wk, wk_ad, vort, vort_ad, gridstruct, npx, &
&                     npy, is, ie, js, je, ng, .false.)
        ELSE
          CALL POPREALARRAY_ADM(vort, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
          CALL SMAG_CORNER_ADM(abs1, u, u_ad, v, v_ad, ua, va, vort, &
&                        vort_ad, bd, npx, npy, gridstruct, ng)
        END IF
      END IF
      DO n=nord,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          CALL POPINTEGER4(ad_from6)
          CALL POPINTEGER4(ad_to6)
          DO j=ad_to6,ad_from6,-1
            CALL POPINTEGER4(ad_from5)
            CALL POPINTEGER4(ad_to5)
            DO i=ad_to5,ad_from5,-1
              divg_d_ad(i, j) = gridstruct%rarea_c(i, j)*divg_d_ad(i, j)
            END DO
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) uc_ad(1, npy) = uc_ad(1, npy) + divg_d_ad(1, &
&           npy)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) uc_ad(npx, npy) = uc_ad(npx, npy) + divg_d_ad&
&           (npx, npy)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) uc_ad(npx, 0) = uc_ad(npx, 0) - divg_d_ad(npx&
&           , 1)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) uc_ad(1, 0) = uc_ad(1, 0) - divg_d_ad(1, 1)
        CALL POPINTEGER4(ad_from4)
        CALL POPINTEGER4(ad_to4)
        DO j=ad_to4,ad_from4,-1
          CALL POPINTEGER4(ad_from3)
          CALL POPINTEGER4(ad_to3)
          DO i=ad_to3,ad_from3,-1
            uc_ad(i, j-1) = uc_ad(i, j-1) + divg_d_ad(i, j)
            uc_ad(i, j) = uc_ad(i, j) - divg_d_ad(i, j)
            vc_ad(i-1, j) = vc_ad(i-1, j) + divg_d_ad(i, j)
            vc_ad(i, j) = vc_ad(i, j) - divg_d_ad(i, j)
            divg_d_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) CALL FILL_CORNERS_ADM(vc, vc_ad, uc, uc_ad, &
&                                          npx, npy, dgrid=.true., &
&                                          vector=.true.)
        CALL POPINTEGER4(ad_from2)
        CALL POPINTEGER4(ad_to2)
        DO j=ad_to2,ad_from2,-1
          CALL POPINTEGER4(ad_from1)
          CALL POPINTEGER4(ad_to1)
          DO i=ad_to1,ad_from1,-1
            temp_ad8 = divg_v(i, j)*uc_ad(i, j)
            divg_d_ad(i, j+1) = divg_d_ad(i, j+1) + temp_ad8
            divg_d_ad(i, j) = divg_d_ad(i, j) - temp_ad8
            uc_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) CALL FILL_CORNERS_ADM(divg_d, divg_d_ad, npx&
&                                          , npy, fill=ydir, bgrid=&
&                                          .true.)
        CALL POPINTEGER4(ad_from0)
        CALL POPINTEGER4(ad_to0)
        DO j=ad_to0,ad_from0,-1
          CALL POPINTEGER4(ad_from)
          CALL POPINTEGER4(ad_to)
          DO i=ad_to,ad_from,-1
            temp_ad7 = divg_u(i, j)*vc_ad(i, j)
            divg_d_ad(i+1, j) = divg_d_ad(i+1, j) + temp_ad7
            divg_d_ad(i, j) = divg_d_ad(i, j) - temp_ad7
            vc_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) CALL FILL_CORNERS_ADM(divg_d, divg_d_ad, npx&
&                                          , npy, fill=xdir, bgrid=&
&                                          .true.)
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY_ADM(delpc(i, j))
          divg_d_ad(i, j) = divg_d_ad(i, j) + delpc_ad(i, j)
          delpc_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
  END SUBROUTINE COMPUTE_DIVERGENCE_DAMPING_ADM
!  Differentiation of smag_corner in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dy
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
!   gradient     of useful results: smag_c u v
!   with respect to varying inputs: u v
  SUBROUTINE SMAG_CORNER_ADM(dt, u, u_ad, v, v_ad, ua, va, smag_c, &
&   smag_c_ad, bd, npx, npy, gridstruct, ng)
    IMPLICIT NONE
! Compute the Tension_Shear strain at cell corners for Smagorinsky diffusion
!!!  work only if (grid_type==4)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: npx, npy, ng
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1) :: u_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed) :: v_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: smag_c
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: smag_c_ad
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! local
    REAL :: ut(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: ut_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: vt_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!  work array
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: sh(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: sh_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    INTEGER :: i, j
    INTEGER :: is2, ie1
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc, dx, dy, rarea, rarea_c
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SQRT
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    dx => gridstruct%dx
    dy => gridstruct%dy
    rarea => gridstruct%rarea
    rarea_c => gridstruct%rarea_c
! Smag = sqrt [ T**2 + S**2 ]:  unit = 1/s
! where T = du/dx - dv/dy;   S = du/dy + dv/dx
! Compute tension strain at corners:
    DO j=js,je+1
      DO i=is-1,ie+1
        ut(i, j) = u(i, j)*dyc(i, j)
      END DO
    END DO
    DO j=js-1,je+1
      DO i=is,ie+1
        vt(i, j) = v(i, j)*dxc(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie+1
        smag_c(i, j) = rarea_c(i, j)*(vt(i, j-1)-vt(i, j)+(ut(i, j)-ut(i&
&         -1, j)))
      END DO
    END DO
! Fix the corners?? if grid_type /= 4
! Compute shear strain:
    DO j=jsd,jed+1
      DO i=isd,ied
        vt(i, j) = u(i, j)*dx(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied+1
        ut(i, j) = v(i, j)*dy(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied
        wk(i, j) = rarea(i, j)*(vt(i, j)-vt(i, j+1)+(ut(i, j)-ut(i+1, j)&
&         ))
      END DO
    END DO
    CALL A2B_ORD4(wk, sh, gridstruct, npx, npy, is, ie, js, je, ng, &
&           .false.)
    sh_ad = 0.0
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        IF (sh(i, j)**2 + smag_c(i, j)**2 .EQ. 0.0) THEN
          temp_ad1 = 0.0
        ELSE
          temp_ad1 = dt*smag_c_ad(i, j)/(2.0*SQRT(sh(i, j)**2+smag_c(i, &
&           j)**2))
        END IF
        sh_ad(i, j) = sh_ad(i, j) + 2*sh(i, j)*temp_ad1
        smag_c_ad(i, j) = 2*smag_c(i, j)*temp_ad1
      END DO
    END DO
    wk_ad = 0.0
    CALL A2B_ORD4_ADM(wk, wk_ad, sh, sh_ad, gridstruct, npx, npy, is, ie&
&               , js, je, ng, .false.)
    ut_ad = 0.0
    vt_ad = 0.0
    DO j=jed,jsd,-1
      DO i=ied,isd,-1
        temp_ad0 = rarea(i, j)*wk_ad(i, j)
        vt_ad(i, j) = vt_ad(i, j) + temp_ad0
        vt_ad(i, j+1) = vt_ad(i, j+1) - temp_ad0
        ut_ad(i, j) = ut_ad(i, j) + temp_ad0
        ut_ad(i+1, j) = ut_ad(i+1, j) - temp_ad0
        wk_ad(i, j) = 0.0
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ied+1,isd,-1
        v_ad(i, j) = v_ad(i, j) + dy(i, j)*ut_ad(i, j)
        ut_ad(i, j) = 0.0
      END DO
    END DO
    DO j=jed+1,jsd,-1
      DO i=ied,isd,-1
        u_ad(i, j) = u_ad(i, j) + dx(i, j)*vt_ad(i, j)
        vt_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        temp_ad = rarea_c(i, j)*smag_c_ad(i, j)
        vt_ad(i, j-1) = vt_ad(i, j-1) + temp_ad
        vt_ad(i, j) = vt_ad(i, j) - temp_ad
        ut_ad(i, j) = ut_ad(i, j) + temp_ad
        ut_ad(i-1, j) = ut_ad(i-1, j) - temp_ad
        smag_c_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js-1,-1
      DO i=ie+1,is,-1
        v_ad(i, j) = v_ad(i, j) + dxc(i, j)*vt_ad(i, j)
        vt_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is-1,-1
        u_ad(i, j) = u_ad(i, j) + dyc(i, j)*ut_ad(i, j)
        ut_ad(i, j) = 0.0
      END DO
    END DO
  END SUBROUTINE SMAG_CORNER_ADM
  SUBROUTINE COMPUTE_DIVERGENCE_DAMPING(nord, d2_bg, d4_bg, dddmp, dt, &
&   vort, ptc, delpc, ke, u, v, uc, vc, ua, va, divg_d, wk, gridstruct, &
&   flagstruct, bd)
    IMPLICIT NONE
!InOut Arguments
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    INTEGER, INTENT(IN) :: nord
    REAL, INTENT(IN) :: d2_bg, d4_bg, dddmp, dt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
!Intent is really in
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: wk
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: vort
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delpc, ptc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   ke
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: vc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   divg_d
!Locals
    REAL :: damp, dd8, damp2, da_min, da_min_c, absdt
    INTEGER :: is, ie, js, je, npx, npy, is2, ie1
    LOGICAL :: nested, fill_c
    INTEGER :: i, j, n, n2, nt
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, DIMENSION(:, :), POINTER :: area, area_c, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsina
    REAL, DIMENSION(:, :), POINTER :: f0, rsin2, divg_u, divg_v
    REAL, DIMENSION(:, :), POINTER :: cosa, dx, dy, dxc, dyc, rdxa, rdya&
&   , rdx, rdy
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SQRT
    REAL :: max1
    REAL :: abs0
    REAL :: abs1
    REAL :: max2
    REAL :: abs2
    REAL :: y3
    REAL :: y2
    REAL :: y1
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    cosa_s => gridstruct%cosa_s
    sina_u => gridstruct%sina_u
    sina_v => gridstruct%sina_v
    rsin_u => gridstruct%rsin_u
    rsin_v => gridstruct%rsin_v
    rsina => gridstruct%rsina
    f0 => gridstruct%f0
    rsin2 => gridstruct%rsin2
    divg_u => gridstruct%divg_u
    divg_v => gridstruct%divg_v
    cosa => gridstruct%cosa
    dx => gridstruct%dx
    dy => gridstruct%dy
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    rdxa => gridstruct%rdxa
    rdya => gridstruct%rdya
    rdx => gridstruct%rdx
    rdy => gridstruct%rdy
    sw_corner = gridstruct%sw_corner
    se_corner = gridstruct%se_corner
    nw_corner = gridstruct%nw_corner
    ne_corner = gridstruct%ne_corner
    IF (dt .GE. 0.) THEN
      absdt = dt
    ELSE
      absdt = -dt
    END IF
    da_min = gridstruct%da_min
    da_min_c = gridstruct%da_min_c
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    IF (nested) THEN
      is2 = is
      ie1 = ie + 1
    ELSE
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
    END IF
!-----------------------------
! Compute divergence damping
!-----------------------------
!  damp = dddmp * da_min_c
    IF (nord .EQ. 0) THEN
!         area ~ dxb*dyb*sin(alpha)
      IF (nested) THEN
        DO j=js,je+1
          DO i=is-1,ie+1
            ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j))&
&             *dyc(i, j)*sina_v(i, j)
          END DO
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
        END DO
      ELSE
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is-1,ie+1
              IF (vc(i, j) .GT. 0) THEN
                ptc(i, j) = u(i, j)*dyc(i, j)*sin_sg(i, j-1, 4)
              ELSE
                ptc(i, j) = u(i, j)*dyc(i, j)*sin_sg(i, j, 2)
              END IF
            END DO
          ELSE
            DO i=is-1,ie+1
              ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j&
&               ))*dyc(i, j)*sina_v(i, j)
            END DO
          END IF
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
          IF (is .EQ. 1) THEN
            IF (uc(1, j) .GT. 0) THEN
              vort(1, j) = v(1, j)*dxc(1, j)*sin_sg(0, j, 3)
            ELSE
              vort(1, j) = v(1, j)*dxc(1, j)*sin_sg(1, j, 1)
            END IF
          END IF
          IF (ie + 1 .EQ. npx) THEN
            IF (uc(npx, j) .GT. 0) THEN
              vort(npx, j) = v(npx, j)*dxc(npx, j)*sin_sg(npx-1, j, 3)
            ELSE
              vort(npx, j) = v(npx, j)*dxc(npx, j)*sin_sg(npx, j, 1)
            END IF
          END IF
        END DO
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          delpc(i, j) = vort(i, j-1) - vort(i, j) + ptc(i-1, j) - ptc(i&
&           , j)
        END DO
      END DO
! Remove the extra term at the corners:
      IF (sw_corner) delpc(1, 1) = delpc(1, 1) - vort(1, 0)
      IF (se_corner) delpc(npx, 1) = delpc(npx, 1) - vort(npx, 0)
      IF (ne_corner) delpc(npx, npy) = delpc(npx, npy) + vort(npx, npy)
      IF (nw_corner) delpc(1, npy) = delpc(1, npy) + vort(1, npy)
      DO j=js,je+1
        DO i=is,ie+1
          delpc(i, j) = gridstruct%rarea_c(i, j)*delpc(i, j)
          IF (delpc(i, j)*dt .GE. 0.) THEN
            abs2 = delpc(i, j)*dt
          ELSE
            abs2 = -(delpc(i, j)*dt)
          END IF
          y3 = dddmp*abs2
          IF (0.20 .GT. y3) THEN
            y1 = y3
          ELSE
            y1 = 0.20
          END IF
          IF (d2_bg .LT. y1) THEN
            max1 = y1
          ELSE
            max1 = d2_bg
          END IF
          damp = gridstruct%da_min_c*max1
          vort(i, j) = damp*delpc(i, j)
          ke(i, j) = ke(i, j) + vort(i, j)
        END DO
      END DO
    ELSE
!--------------------------
! Higher order divg damping
!--------------------------
      DO j=js,je+1
        DO i=is,ie+1
! Save divergence for external mode filter
          delpc(i, j) = divg_d(i, j)
        END DO
      END DO
! N > 1
      n2 = nord + 1
      DO n=1,nord
        nt = nord - n
        fill_c = nt .NE. 0 .AND. flagstruct%grid_type .LT. 3 .AND. (((&
&         sw_corner .OR. se_corner) .OR. ne_corner) .OR. nw_corner) &
&         .AND. (.NOT.nested)
        IF (fill_c) CALL FILL_CORNERS(divg_d, npx, npy, fill=xdir, bgrid&
&                               =.true.)
        DO j=js-nt,je+1+nt
          DO i=is-1-nt,ie+1+nt
            vc(i, j) = (divg_d(i+1, j)-divg_d(i, j))*divg_u(i, j)
          END DO
        END DO
        IF (fill_c) CALL FILL_CORNERS(divg_d, npx, npy, fill=ydir, bgrid&
&                               =.true.)
        DO j=js-1-nt,je+1+nt
          DO i=is-nt,ie+1+nt
            uc(i, j) = (divg_d(i, j+1)-divg_d(i, j))*divg_v(i, j)
          END DO
        END DO
        IF (fill_c) CALL FILL_CORNERS(vc, uc, npx, npy, dgrid=.true., &
&                               vector=.true.)
        DO j=js-nt,je+1+nt
          DO i=is-nt,ie+1+nt
            divg_d(i, j) = uc(i, j-1) - uc(i, j) + vc(i-1, j) - vc(i, j)
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) divg_d(1, 1) = divg_d(1, 1) - uc(1, 0)
        IF (se_corner) divg_d(npx, 1) = divg_d(npx, 1) - uc(npx, 0)
        IF (ne_corner) divg_d(npx, npy) = divg_d(npx, npy) + uc(npx, npy&
&           )
        IF (nw_corner) divg_d(1, npy) = divg_d(1, npy) + uc(1, npy)
        IF (.NOT.gridstruct%stretched_grid) THEN
          DO j=js-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              divg_d(i, j) = divg_d(i, j)*gridstruct%rarea_c(i, j)
            END DO
          END DO
        END IF
      END DO
! n-loop
      IF (dddmp .LT. 1.e-5) THEN
        vort(:, :) = 0.
      ELSE IF (flagstruct%grid_type .LT. 3) THEN
! Interpolate relative vort to cell corners
        CALL A2B_ORD4(wk, vort, gridstruct, npx, npy, is, ie, js, je, ng&
&               , .false.)
        DO j=js,je+1
          DO i=is,ie+1
            IF (dt .GE. 0.) THEN
              abs0 = dt
            ELSE
              abs0 = -dt
            END IF
! The following is an approxi form of Smagorinsky diffusion
            vort(i, j) = abs0*SQRT(delpc(i, j)**2+vort(i, j)**2)
          END DO
        END DO
      ELSE
        IF (dt .GE. 0.) THEN
          abs1 = dt
        ELSE
          abs1 = -dt
        END IF
! Correct form: works only for doubly preiodic domain
        CALL SMAG_CORNER(abs1, u, v, ua, va, vort, bd, npx, npy, &
&                  gridstruct, ng)
      END IF
      IF (gridstruct%stretched_grid) THEN
! Stretched grid with variable damping ~ area
        dd8 = gridstruct%da_min*d4_bg**n2
      ELSE
        dd8 = (gridstruct%da_min_c*d4_bg)**n2
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          IF (0.20 .GT. dddmp*vort(i, j)) THEN
            y2 = dddmp*vort(i, j)
          ELSE
            y2 = 0.20
          END IF
          IF (d2_bg .LT. y2) THEN
            max2 = y2
          ELSE
            max2 = d2_bg
          END IF
! del-2
          damp2 = gridstruct%da_min_c*max2
          vort(i, j) = damp2*delpc(i, j) + dd8*divg_d(i, j)
          ke(i, j) = ke(i, j) + vort(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE COMPUTE_DIVERGENCE_DAMPING
!  Differentiation of compute_divergence_damping in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b
!_ord4_fb a2b_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_m
!od.adv_pe dyn_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_up
!date dyn_core_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dy
!namics_mod.Rayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_gri
!d_utils_mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.
!pkez fv_mapz_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_
!mapz_mod.remap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_map
!z_mod.ppm_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_
!mod.map1_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested 
!fv_sg_mod.fv_subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.up
!date_dz_d nh_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3
!_solver nh_utils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_p
!rofile nh_utils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_
!nest sw_core_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u
! sw_core_mod.ytp_v sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.f
!v_tp_2d_fb tp_core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corn
!er_fb fv_grid_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: u v ke ua uc ptc delpc va vc
!                vort divg_d wk
!   with respect to varying inputs: u v ke ua uc ptc delpc va vc
!                divg_d wk
  SUBROUTINE COMPUTE_DIVERGENCE_DAMPING_FWD(nord, d2_bg, d4_bg, dddmp&
&   , dt, vort, ptc, delpc, ke, u, v, uc, vc, ua, va, divg_d, wk, &
&   gridstruct, flagstruct, bd)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
!InOut Arguments
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    INTEGER, INTENT(IN) :: nord
    REAL, INTENT(IN) :: d2_bg, d4_bg, dddmp, dt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
!Intent is really in
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: wk
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: vort
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delpc, ptc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   ke
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: vc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   divg_d
!Locals
    REAL :: damp, dd8, damp2, da_min, da_min_c, absdt
    INTEGER :: is, ie, js, je, npx, npy, is2, ie1
    LOGICAL :: nested, fill_c
    INTEGER :: i, j, n, n2, nt
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, DIMENSION(:, :), POINTER :: area, area_c, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsina
    REAL, DIMENSION(:, :), POINTER :: f0, rsin2, divg_u, divg_v
    REAL, DIMENSION(:, :), POINTER :: cosa, dx, dy, dxc, dyc, rdxa, rdya&
&   , rdx, rdy
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SQRT
    REAL :: max1
    REAL :: abs0
    REAL :: abs1
    REAL :: max2
    REAL :: abs2
    INTEGER :: ad_from
    INTEGER :: ad_from0
    INTEGER :: ad_from1
    INTEGER :: ad_from2
    INTEGER :: ad_from3
    INTEGER :: ad_from4
    INTEGER :: ad_from5
    INTEGER :: ad_from6
    REAL :: y3
    REAL :: y2
    REAL :: y1

    damp = 0.0
    dd8 = 0.0
    damp2 = 0.0
    da_min = 0.0
    da_min_c = 0.0
    absdt = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    npx = 0
    npy = 0
    is2 = 0
    ie1 = 0
    n = 0
    n2 = 0
    nt = 0
    max1 = 0.0
    abs0 = 0.0
    abs1 = 0.0
    max2 = 0.0
    abs2 = 0.0
    ad_from = 0
    ad_from0 = 0
    ad_from1 = 0
    ad_from2 = 0
    ad_from3 = 0
    ad_from4 = 0
    ad_from5 = 0
    ad_from6 = 0
    y3 = 0.0
    y2 = 0.0
    y1 = 0.0

    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    sina_u => gridstruct%sina_u
    sina_v => gridstruct%sina_v
    divg_u => gridstruct%divg_u
    divg_v => gridstruct%divg_v
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    sw_corner = gridstruct%sw_corner
    se_corner = gridstruct%se_corner
    nw_corner = gridstruct%nw_corner
    ne_corner = gridstruct%ne_corner
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    IF (nested) THEN
      CALL PUSHCONTROL(1,0)
      is2 = is
      ie1 = ie + 1
    ELSE
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        CALL PUSHCONTROL(1,1)
        ie1 = ie + 1
      ELSE
        CALL PUSHCONTROL(1,1)
        ie1 = npx - 1
      END IF
    END IF
!-----------------------------
! Compute divergence damping
!-----------------------------
!  damp = dddmp * da_min_c
    IF (nord .EQ. 0) THEN
!         area ~ dxb*dyb*sin(alpha)
      IF (nested) THEN
        DO j=js,je+1
          DO i=is-1,ie+1
            CALL PUSHREALARRAY(ptc(i, j))
            ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j))&
&             *dyc(i, j)*sina_v(i, j)
          END DO
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is-1,ie+1
              IF (vc(i, j) .GT. 0) THEN
                CALL PUSHREALARRAY(ptc(i, j))
                ptc(i, j) = u(i, j)*dyc(i, j)*sin_sg(i, j-1, 4)
                CALL PUSHCONTROL(1,1)
              ELSE
                CALL PUSHREALARRAY(ptc(i, j))
                ptc(i, j) = u(i, j)*dyc(i, j)*sin_sg(i, j, 2)
                CALL PUSHCONTROL(1,0)
              END IF
            END DO
            CALL PUSHCONTROL(1,1)
          ELSE
            DO i=is-1,ie+1
              CALL PUSHREALARRAY(ptc(i, j))
              ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j&
&               ))*dyc(i, j)*sina_v(i, j)
            END DO
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
          IF (is .EQ. 1) THEN
            IF (uc(1, j) .GT. 0) THEN
              vort(1, j) = v(1, j)*dxc(1, j)*sin_sg(0, j, 3)
              CALL PUSHCONTROL(2,0)
            ELSE
              vort(1, j) = v(1, j)*dxc(1, j)*sin_sg(1, j, 1)
              CALL PUSHCONTROL(2,1)
            END IF
          ELSE
            CALL PUSHCONTROL(2,2)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            IF (uc(npx, j) .GT. 0) THEN
              vort(npx, j) = v(npx, j)*dxc(npx, j)*sin_sg(npx-1, j, 3)
              CALL PUSHCONTROL(2,2)
            ELSE
              vort(npx, j) = v(npx, j)*dxc(npx, j)*sin_sg(npx, j, 1)
              CALL PUSHCONTROL(2,1)
            END IF
          ELSE
            CALL PUSHCONTROL(2,0)
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY(delpc(i, j))
          delpc(i, j) = vort(i, j-1) - vort(i, j) + (ptc(i-1, j)-ptc(i, &
&           j))
        END DO
      END DO
! Remove the extra term at the corners:
      IF (sw_corner) THEN
        CALL PUSHREALARRAY(delpc(1, 1))
        delpc(1, 1) = delpc(1, 1) - vort(1, 0)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (se_corner) THEN
        CALL PUSHREALARRAY(delpc(npx, 1))
        delpc(npx, 1) = delpc(npx, 1) - vort(npx, 0)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHREALARRAY(delpc(npx, npy))
        delpc(npx, npy) = delpc(npx, npy) + vort(npx, npy)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (nw_corner) THEN
        CALL PUSHREALARRAY(delpc(1, npy))
        delpc(1, npy) = delpc(1, npy) + vort(1, npy)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          CALL PUSHREALARRAY(delpc(i, j))
          delpc(i, j) = gridstruct%rarea_c(i, j)*delpc(i, j)
          IF (delpc(i, j)*dt .GE. 0.) THEN
            abs2 = delpc(i, j)*dt
            CALL PUSHCONTROL(1,0)
          ELSE
            abs2 = -(delpc(i, j)*dt)
            CALL PUSHCONTROL(1,1)
          END IF
          y3 = dddmp*abs2
          IF (0.20 .GT. y3) THEN
            y1 = y3
            CALL PUSHCONTROL(1,0)
          ELSE
            y1 = 0.20
            CALL PUSHCONTROL(1,1)
          END IF
          IF (d2_bg .LT. y1) THEN
            max1 = y1
            CALL PUSHCONTROL(1,0)
          ELSE
            max1 = d2_bg
            CALL PUSHCONTROL(1,1)
          END IF
          CALL PUSHREALARRAY(damp)
          damp = gridstruct%da_min_c*max1
          vort(i, j) = damp*delpc(i, j)
          ke(i, j) = ke(i, j) + vort(i, j)
        END DO
      END DO
      CALL PUSHINTEGER(je)
      CALL PUSHINTEGER(ie1)
      CALL PUSHINTEGER(is)
      !CALL PUSHPOINTER8(C_LOC(sina_v))
      !CALL PUSHPOINTER8(C_LOC(sina_u))
      CALL PUSHINTEGER(ie)
      !CALL PUSHPOINTER8(C_LOC(dyc))
      !CALL PUSHPOINTER8(C_LOC(sin_sg))
      !CALL PUSHPOINTER8(C_LOC(cosa_v))
      !CALL PUSHPOINTER8(C_LOC(cosa_u))
      CALL PUSHINTEGER(is2)
      !CALL PUSHPOINTER8(C_LOC(dxc))
      CALL PUSHREALARRAY(damp)
      CALL PUSHINTEGER(npy)
      CALL PUSHINTEGER(npx)
      CALL PUSHINTEGER(js)
      CALL PUSHCONTROL(1,0)
    ELSE
!--------------------------
! Higher order divg damping
!--------------------------
      DO j=js,je+1
        DO i=is,ie+1
! Save divergence for external mode filter
          CALL PUSHREALARRAY(delpc(i, j))
          delpc(i, j) = divg_d(i, j)
        END DO
      END DO
! N > 1
      n2 = nord + 1
      DO n=1,nord
        nt = nord - n
        fill_c = nt .NE. 0 .AND. flagstruct%grid_type .LT. 3 .AND. (((&
&         sw_corner .OR. se_corner) .OR. ne_corner) .OR. nw_corner) &
&         .AND. (.NOT.nested)
        IF (fill_c) THEN
          CALL FILL_CORNERS(divg_d, npx, npy, fill=xdir, bgrid=.true.)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        ad_from0 = js - nt
        DO j=ad_from0,je+1+nt
          ad_from = is - 1 - nt
          DO i=ad_from,ie+1+nt
            vc(i, j) = (divg_d(i+1, j)-divg_d(i, j))*divg_u(i, j)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from)
        END DO
        CALL PUSHINTEGER(j - 1)
        CALL PUSHINTEGER(ad_from0)
        IF (fill_c) THEN
          CALL FILL_CORNERS(divg_d, npx, npy, fill=ydir, bgrid=.true.)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        ad_from2 = js - 1 - nt
        DO j=ad_from2,je+1+nt
          ad_from1 = is - nt
          DO i=ad_from1,ie+1+nt
            uc(i, j) = (divg_d(i, j+1)-divg_d(i, j))*divg_v(i, j)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from1)
        END DO
        CALL PUSHINTEGER(j - 1)
        CALL PUSHINTEGER(ad_from2)
        IF (fill_c) THEN
          CALL FILL_CORNERS(vc, uc, npx, npy, vector=.true., dgrid=&
&                     .true.)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        ad_from4 = js - nt
        DO j=ad_from4,je+1+nt
          ad_from3 = is - nt
          DO i=ad_from3,ie+1+nt
            divg_d(i, j) = uc(i, j-1) - uc(i, j) + (vc(i-1, j)-vc(i, j))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from3)
        END DO
        CALL PUSHINTEGER(j - 1)
        CALL PUSHINTEGER(ad_from4)
! Remove the extra term at the corners:
        IF (sw_corner) THEN
          divg_d(1, 1) = divg_d(1, 1) - uc(1, 0)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (se_corner) THEN
          divg_d(npx, 1) = divg_d(npx, 1) - uc(npx, 0)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (ne_corner) THEN
          divg_d(npx, npy) = divg_d(npx, npy) + uc(npx, npy)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (nw_corner) THEN
          divg_d(1, npy) = divg_d(1, npy) + uc(1, npy)
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (.NOT.gridstruct%stretched_grid) THEN
          ad_from6 = js - nt
          DO j=ad_from6,je+1+nt
            ad_from5 = is - nt
            DO i=ad_from5,ie+1+nt
              divg_d(i, j) = divg_d(i, j)*gridstruct%rarea_c(i, j)
            END DO
            CALL PUSHINTEGER(i - 1)
            CALL PUSHINTEGER(ad_from5)
          END DO
          CALL PUSHINTEGER(j - 1)
          CALL PUSHINTEGER(ad_from6)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
! n-loop
      IF (dddmp .LT. 1.e-5) THEN
        vort(:, :) = 0.
        CALL PUSHCONTROL(2,0)
      ELSE IF (flagstruct%grid_type .LT. 3) THEN
! Interpolate relative vort to cell corners
        CALL A2B_ORD4_FWD(wk, vort, gridstruct, npx, npy, is, ie, js&
&                      , je, ng, .false.)
        DO j=js,je+1
          DO i=is,ie+1
            IF (dt .GE. 0.) THEN
              CALL PUSHREALARRAY(abs0)
              abs0 = dt
              CALL PUSHCONTROL(1,0)
            ELSE
              CALL PUSHREALARRAY(abs0)
              abs0 = -dt
              CALL PUSHCONTROL(1,1)
            END IF
! The following is an approxi form of Smagorinsky diffusion
            CALL PUSHREALARRAY(vort(i, j))
            vort(i, j) = abs0*SQRT(delpc(i, j)**2+vort(i, j)**2)
          END DO
        END DO
        CALL PUSHCONTROL(2,1)
      ELSE
        IF (dt .GE. 0.) THEN
          abs1 = dt
        ELSE
          abs1 = -dt
        END IF
! Correct form: works only for doubly preiodic domain
        CALL SMAG_CORNER_FWD(abs1, u, v, ua, va, vort, bd, npx, npy, &
&                         gridstruct, ng)
        CALL PUSHCONTROL(2,2)
      END IF
      IF (gridstruct%stretched_grid) THEN
! Stretched grid with variable damping ~ area
        dd8 = gridstruct%da_min*d4_bg**n2
      ELSE
        dd8 = (gridstruct%da_min_c*d4_bg)**n2
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          IF (0.20 .GT. dddmp*vort(i, j)) THEN
            y2 = dddmp*vort(i, j)
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
            y2 = 0.20
          END IF
          IF (d2_bg .LT. y2) THEN
            max2 = y2
            CALL PUSHCONTROL(1,0)
          ELSE
            max2 = d2_bg
            CALL PUSHCONTROL(1,1)
          END IF
! del-2
          CALL PUSHREALARRAY(damp2)
          damp2 = gridstruct%da_min_c*max2
          vort(i, j) = damp2*delpc(i, j) + dd8*divg_d(i, j)
          ke(i, j) = ke(i, j) + vort(i, j)
        END DO
      END DO
      CALL PUSHREALARRAY(damp2)
      CALL PUSHINTEGER(je)
      CALL PUSHINTEGER(is)
      CALL PUSHREALARRAY(dd8)
      CALL PUSHINTEGER(ie)
      CALL PUSHREALARRAY(abs1)
      CALL PUSHREALARRAY(abs0)
      CALL PUSHINTEGER(js)
      CALL PUSHCONTROL(1,1)
    END IF
  END SUBROUTINE COMPUTE_DIVERGENCE_DAMPING_FWD
!  Differentiation of compute_divergence_damping in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2
!b_ord4_fb a2b_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_
!mod.adv_pe dyn_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_u
!pdate dyn_core_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_d
!ynamics_mod.Rayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_gr
!id_utils_mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod
!.pkez fv_mapz_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv
!_mapz_mod.remap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_ma
!pz_mod.ppm_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz
!_mod.map1_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested
! fv_sg_mod.fv_subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.u
!pdate_dz_d nh_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM
!3_solver nh_utils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_
!profile nh_utils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner
!_nest sw_core_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u_f
!b sw_core_mod.ytp_v sw_core_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.
!fv_tp_2d tp_core_mod.copy_corners tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_cor
!ner_fb fv_grid_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: u v ke ua uc ptc delpc va vc
!                vort divg_d wk
!   with respect to varying inputs: u v ke ua uc ptc delpc va vc
!                divg_d wk
  SUBROUTINE COMPUTE_DIVERGENCE_DAMPING_BWD(nord, d2_bg, d4_bg, dddmp&
&   , dt, vort, vort_ad, ptc, ptc_ad, delpc, delpc_ad, ke, ke_ad, u, &
&   u_ad, v, v_ad, uc, uc_ad, vc, vc_ad, ua, ua_ad, va, va_ad, divg_d, &
&   divg_d_ad, wk, wk_ad, gridstruct, flagstruct, bd)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    INTEGER, INTENT(IN) :: nord
    REAL, INTENT(IN) :: d2_bg, d4_bg, dddmp, dt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: ua_ad, va_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1) :: u_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed) :: v_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: wk
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   wk_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: vort
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   vort_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delpc, ptc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delpc_ad, ptc_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   ke
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   ke_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   vc_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   uc_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   divg_d
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   divg_d_ad
    REAL :: damp, dd8, damp2, da_min, da_min_c, absdt
    REAL :: damp_ad, damp2_ad
    INTEGER :: is, ie, js, je, npx, npy, is2, ie1
    LOGICAL :: nested, fill_c
    INTEGER :: i, j, n, n2, nt
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, DIMENSION(:, :), POINTER :: area, area_c, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: cosa_u, cosa_v, cosa_s
    REAL, DIMENSION(:, :), POINTER :: sina_u, sina_v
    REAL, DIMENSION(:, :), POINTER :: rsin_u, rsin_v, rsina
    REAL, DIMENSION(:, :), POINTER :: f0, rsin2, divg_u, divg_v
    REAL, DIMENSION(:, :), POINTER :: cosa, dx, dy, dxc, dyc, rdxa, rdya&
&   , rdx, rdy
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SQRT
    REAL :: max1
    REAL :: max1_ad
    REAL :: abs0
    REAL :: abs1
    REAL :: max2
    REAL :: max2_ad
    REAL :: abs2
    REAL :: abs2_ad
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: y3_ad
    REAL :: y1_ad
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: y2_ad
    INTEGER :: branch
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
    INTEGER :: ad_from5
    INTEGER :: ad_to5
    INTEGER :: ad_from6
    INTEGER :: ad_to6
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_compute_divergence_damping
    REAL :: y3
    REAL :: y2
    REAL :: y1

    damp = 0.0
    dd8 = 0.0
    damp2 = 0.0
    da_min = 0.0
    da_min_c = 0.0
    absdt = 0.0
    is = 0
    ie = 0
    js = 0
    je = 0
    npx = 0
    npy = 0
    is2 = 0
    ie1 = 0
    n = 0
    n2 = 0
    nt = 0
    max1 = 0.0
    abs0 = 0.0
    abs1 = 0.0
    max2 = 0.0
    abs2 = 0.0
    ad_from = 0
    ad_from0 = 0
    ad_from1 = 0
    ad_from2 = 0
    ad_from3 = 0
    ad_from4 = 0
    ad_from5 = 0
    ad_from6 = 0
    ad_to = 0
    ad_to0 = 0
    ad_to1 = 0
    ad_to2 = 0
    ad_to3 = 0
    ad_to4 = 0
    ad_to5 = 0
    ad_to6 = 0
    y3 = 0.0
    y2 = 0.0
    y1 = 0.0
    branch = 0

    sin_sg => gridstruct%sin_sg
    cosa_u => gridstruct%cosa_u
    cosa_v => gridstruct%cosa_v
    sina_u => gridstruct%sina_u
    sina_v => gridstruct%sina_v
    divg_u => gridstruct%divg_u
    divg_v => gridstruct%divg_v
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    sw_corner = gridstruct%sw_corner
    se_corner = gridstruct%se_corner
    nw_corner = gridstruct%nw_corner
    ne_corner = gridstruct%ne_corner
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    npx = flagstruct%npx
    npy = flagstruct%npy
    nested = gridstruct%nested
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(js)
      CALL POPINTEGER(npx)
      CALL POPINTEGER(npy)
      CALL POPREALARRAY(damp)
      !CALL POPPOINTER8(cptr)
      dxc => gridstruct%dxc ! (/&
!&                unknown_shape_in_compute_divergence_damping/))
      CALL POPINTEGER(is2)
      !CALL POPPOINTER8(cptr)
      cosa_u => gridstruct%cosa_u ! (/&
!&                unknown_shape_in_compute_divergence_damping/))
      !CALL POPPOINTER8(cptr)
      cosa_v => gridstruct%cosa_v ! (/&
!&                unknown_shape_in_compute_divergence_damping/))
      !CALL POPPOINTER8(cptr)
      sin_sg => gridstruct%sin_sg ! (/&
!&                unknown_shape_in_compute_divergence_damping/))
      !CALL POPPOINTER8(cptr)
      dyc => gridstruct%dyc ! (/&
!&                unknown_shape_in_compute_divergence_damping/))
      CALL POPINTEGER(ie)
      !CALL POPPOINTER8(cptr)
      sina_u => gridstruct%sina_u ! (/&
!&                unknown_shape_in_compute_divergence_damping/))
      !CALL POPPOINTER8(cptr)
      sina_v => gridstruct%sina_v ! (/&
!&                unknown_shape_in_compute_divergence_damping/))
      CALL POPINTEGER(is)
      CALL POPINTEGER(ie1)
      CALL POPINTEGER(je)
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          vort_ad(i, j) = vort_ad(i, j) + ke_ad(i, j)
          damp_ad = delpc(i, j)*vort_ad(i, j)
          delpc_ad(i, j) = delpc_ad(i, j) + damp*vort_ad(i, j)
          vort_ad(i, j) = 0.0
          CALL POPREALARRAY(damp)
          max1_ad = gridstruct%da_min_c*damp_ad
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            y1_ad = max1_ad
          ELSE
            y1_ad = 0.0
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            y3_ad = y1_ad
          ELSE
            y3_ad = 0.0
          END IF
          abs2_ad = dddmp*y3_ad
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            delpc_ad(i, j) = delpc_ad(i, j) + dt*abs2_ad
          ELSE
            delpc_ad(i, j) = delpc_ad(i, j) - dt*abs2_ad
          END IF
          CALL POPREALARRAY(delpc(i, j))
          delpc_ad(i, j) = gridstruct%rarea_c(i, j)*delpc_ad(i, j)
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        npy = flagstruct%npy
        CALL POPREALARRAY(delpc(1, npy))
        vort_ad(1, npy) = vort_ad(1, npy) + delpc_ad(1, npy)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        npx = flagstruct%npx
        CALL POPREALARRAY(delpc(npx, npy))
        vort_ad(npx, npy) = vort_ad(npx, npy) + delpc_ad(npx, npy)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(delpc(npx, 1))
        vort_ad(npx, 0) = vort_ad(npx, 0) - delpc_ad(npx, 1)
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(delpc(1, 1))
        vort_ad(1, 0) = vort_ad(1, 0) - delpc_ad(1, 1)
      END IF
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(delpc(i, j))
          vort_ad(i, j-1) = vort_ad(i, j-1) + delpc_ad(i, j)
          vort_ad(i, j) = vort_ad(i, j) - delpc_ad(i, j)
          ptc_ad(i-1, j) = ptc_ad(i-1, j) + delpc_ad(i, j)
          ptc_ad(i, j) = ptc_ad(i, j) - delpc_ad(i, j)
          delpc_ad(i, j) = 0.0
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=je+1,js-1,-1
          CALL POPCONTROL(2,branch)
          IF (branch .NE. 0) THEN
            IF (branch .EQ. 1) THEN
              v_ad(npx, j) = v_ad(npx, j) + sin_sg(npx, j, 1)*dxc(npx, j&
&               )*vort_ad(npx, j)
              vort_ad(npx, j) = 0.0
            ELSE
              v_ad(npx, j) = v_ad(npx, j) + sin_sg(npx-1, j, 3)*dxc(npx&
&               , j)*vort_ad(npx, j)
              vort_ad(npx, j) = 0.0
            END IF
          END IF
          CALL POPCONTROL(2,branch)
          IF (branch .EQ. 0) THEN
            v_ad(1, j) = v_ad(1, j) + sin_sg(0, j, 3)*dxc(1, j)*vort_ad(&
&             1, j)
            vort_ad(1, j) = 0.0
          ELSE IF (branch .EQ. 1) THEN
            v_ad(1, j) = v_ad(1, j) + sin_sg(1, j, 1)*dxc(1, j)*vort_ad(&
&             1, j)
            vort_ad(1, j) = 0.0
          END IF
          DO i=ie1,is2,-1
            temp_ad5 = dxc(i, j)*sina_u(i, j)*vort_ad(i, j)
            temp_ad6 = -(cosa_u(i, j)*0.5*temp_ad5)
            v_ad(i, j) = v_ad(i, j) + temp_ad5
            ua_ad(i-1, j) = ua_ad(i-1, j) + temp_ad6
            ua_ad(i, j) = ua_ad(i, j) + temp_ad6
            vort_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je+1,js,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            DO i=ie+1,is-1,-1
              CALL POPREALARRAY(ptc(i, j))
              temp_ad3 = dyc(i, j)*sina_v(i, j)*ptc_ad(i, j)
              temp_ad4 = -(cosa_v(i, j)*0.5*temp_ad3)
              u_ad(i, j) = u_ad(i, j) + temp_ad3
              va_ad(i, j-1) = va_ad(i, j-1) + temp_ad4
              va_ad(i, j) = va_ad(i, j) + temp_ad4
              ptc_ad(i, j) = 0.0
            END DO
          ELSE
            DO i=ie+1,is-1,-1
              CALL POPCONTROL(1,branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY(ptc(i, j))
                u_ad(i, j) = u_ad(i, j) + sin_sg(i, j, 2)*dyc(i, j)*&
&                 ptc_ad(i, j)
                ptc_ad(i, j) = 0.0
              ELSE
                CALL POPREALARRAY(ptc(i, j))
                u_ad(i, j) = u_ad(i, j) + sin_sg(i, j-1, 4)*dyc(i, j)*&
&                 ptc_ad(i, j)
                ptc_ad(i, j) = 0.0
              END IF
            END DO
          END IF
        END DO
      ELSE
        DO j=je+1,js-1,-1
          DO i=ie1,is2,-1
            temp_ad1 = dxc(i, j)*sina_u(i, j)*vort_ad(i, j)
            temp_ad2 = -(cosa_u(i, j)*0.5*temp_ad1)
            v_ad(i, j) = v_ad(i, j) + temp_ad1
            ua_ad(i-1, j) = ua_ad(i-1, j) + temp_ad2
            ua_ad(i, j) = ua_ad(i, j) + temp_ad2
            vort_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ie+1,is-1,-1
            CALL POPREALARRAY(ptc(i, j))
            temp_ad = dyc(i, j)*sina_v(i, j)*ptc_ad(i, j)
            temp_ad0 = -(cosa_v(i, j)*0.5*temp_ad)
            u_ad(i, j) = u_ad(i, j) + temp_ad
            va_ad(i, j-1) = va_ad(i, j-1) + temp_ad0
            va_ad(i, j) = va_ad(i, j) + temp_ad0
            ptc_ad(i, j) = 0.0
          END DO
        END DO
      END IF
    ELSE
      CALL POPINTEGER(js)
      CALL POPREALARRAY(abs0)
      CALL POPREALARRAY(abs1)
      CALL POPINTEGER(ie)
      CALL POPREALARRAY(dd8)
      CALL POPINTEGER(is)
      CALL POPINTEGER(je)
      CALL POPREALARRAY(damp2)
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          vort_ad(i, j) = vort_ad(i, j) + ke_ad(i, j)
          damp2_ad = delpc(i, j)*vort_ad(i, j)
          delpc_ad(i, j) = delpc_ad(i, j) + damp2*vort_ad(i, j)
          divg_d_ad(i, j) = divg_d_ad(i, j) + dd8*vort_ad(i, j)
          vort_ad(i, j) = 0.0
          CALL POPREALARRAY(damp2)
          max2_ad = gridstruct%da_min_c*damp2_ad
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            y2_ad = max2_ad
          ELSE
            y2_ad = 0.0
          END IF
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) vort_ad(i, j) = vort_ad(i, j) + dddmp*y2_ad
        END DO
      END DO
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        npx = flagstruct%npx
        npy = flagstruct%npy
      ELSE IF (branch .EQ. 1) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPREALARRAY(vort(i, j))
            IF (delpc(i, j)**2 + vort(i, j)**2 .EQ. 0.0) THEN
              temp_ad9 = 0.0
            ELSE
              temp_ad9 = abs0*vort_ad(i, j)/(2.0*SQRT(delpc(i, j)**2+&
&               vort(i, j)**2))
            END IF
            delpc_ad(i, j) = delpc_ad(i, j) + 2*delpc(i, j)*temp_ad9
            vort_ad(i, j) = 2*vort(i, j)*temp_ad9
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY(abs0)
            ELSE
              CALL POPREALARRAY(abs0)
            END IF
          END DO
        END DO
        npx = flagstruct%npx
        npy = flagstruct%npy
        CALL A2B_ORD4_BWD(wk, wk_ad, vort, vort_ad, gridstruct, npx, &
&                      npy, is, ie, js, je, ng, .false.)
      ELSE
        npx = flagstruct%npx
        npy = flagstruct%npy
        CALL SMAG_CORNER_BWD(abs1, u, u_ad, v, v_ad, ua, va, vort, &
&                         vort_ad, bd, npx, npy, gridstruct, ng)
      END IF
      divg_u => gridstruct%divg_u
      divg_v => gridstruct%divg_v
      DO n=nord,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          CALL POPINTEGER(ad_from6)
          CALL POPINTEGER(ad_to6)
          DO j=ad_to6,ad_from6,-1
            CALL POPINTEGER(ad_from5)
            CALL POPINTEGER(ad_to5)
            DO i=ad_to5,ad_from5,-1
              divg_d_ad(i, j) = gridstruct%rarea_c(i, j)*divg_d_ad(i, j)
            END DO
          END DO
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) uc_ad(1, npy) = uc_ad(1, npy) + divg_d_ad(1, &
&           npy)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) uc_ad(npx, npy) = uc_ad(npx, npy) + divg_d_ad&
&           (npx, npy)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) uc_ad(npx, 0) = uc_ad(npx, 0) - divg_d_ad(npx&
&           , 1)
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) uc_ad(1, 0) = uc_ad(1, 0) - divg_d_ad(1, 1)
        CALL POPINTEGER(ad_from4)
        CALL POPINTEGER(ad_to4)
        DO j=ad_to4,ad_from4,-1
          CALL POPINTEGER(ad_from3)
          CALL POPINTEGER(ad_to3)
          DO i=ad_to3,ad_from3,-1
            uc_ad(i, j-1) = uc_ad(i, j-1) + divg_d_ad(i, j)
            uc_ad(i, j) = uc_ad(i, j) - divg_d_ad(i, j)
            vc_ad(i-1, j) = vc_ad(i-1, j) + divg_d_ad(i, j)
            vc_ad(i, j) = vc_ad(i, j) - divg_d_ad(i, j)
            divg_d_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) CALL FILL_CORNERS_ADM(vc, vc_ad, uc, uc_ad, &
&                                          npx, npy, dgrid=.true., &
&                                          vector=.true.)
        CALL POPINTEGER(ad_from2)
        CALL POPINTEGER(ad_to2)
        DO j=ad_to2,ad_from2,-1
          CALL POPINTEGER(ad_from1)
          CALL POPINTEGER(ad_to1)
          DO i=ad_to1,ad_from1,-1
            temp_ad8 = divg_v(i, j)*uc_ad(i, j)
            divg_d_ad(i, j+1) = divg_d_ad(i, j+1) + temp_ad8
            divg_d_ad(i, j) = divg_d_ad(i, j) - temp_ad8
            uc_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) CALL FILL_CORNERS_ADM(divg_d, divg_d_ad, npx&
&                                          , npy, fill=ydir, bgrid=&
&                                          .true.)
        CALL POPINTEGER(ad_from0)
        CALL POPINTEGER(ad_to0)
        DO j=ad_to0,ad_from0,-1
          CALL POPINTEGER(ad_from)
          CALL POPINTEGER(ad_to)
          DO i=ad_to,ad_from,-1
            temp_ad7 = divg_u(i, j)*vc_ad(i, j)
            divg_d_ad(i+1, j) = divg_d_ad(i+1, j) + temp_ad7
            divg_d_ad(i, j) = divg_d_ad(i, j) - temp_ad7
            vc_ad(i, j) = 0.0
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) CALL FILL_CORNERS_ADM(divg_d, divg_d_ad, npx&
&                                          , npy, fill=xdir, bgrid=&
&                                          .true.)
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(delpc(i, j))
          divg_d_ad(i, j) = divg_d_ad(i, j) + delpc_ad(i, j)
          delpc_ad(i, j) = 0.0
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
  END SUBROUTINE COMPUTE_DIVERGENCE_DAMPING_BWD
!  Differentiation of smag_corner in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_ed
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
!   gradient     of useful results: smag_c u v
!   with respect to varying inputs: u v
  SUBROUTINE SMAG_CORNER_FWD(dt, u, v, ua, va, smag_c, bd, npx, npy, &
&   gridstruct, ng)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
! Compute the Tension_Shear strain at cell corners for Smagorinsky diffusion
!!!  work only if (grid_type==4)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: npx, npy, ng
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: smag_c
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! local
    REAL :: ut(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!  work array
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: sh(bd%isd:bd%ied, bd%jsd:bd%jed)
    INTEGER :: i, j
    INTEGER :: is2, ie1
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc, dx, dy, rarea, rarea_c
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SQRT

    ut = 0.0
    vt = 0.0
    wk = 0.0
    sh = 0.0
    is2 = 0
    ie1 = 0
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
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    dx => gridstruct%dx
    dy => gridstruct%dy
    rarea => gridstruct%rarea
    rarea_c => gridstruct%rarea_c
! Smag = sqrt [ T**2 + S**2 ]:  unit = 1/s
! where T = du/dx - dv/dy;   S = du/dy + dv/dx
! Compute tension strain at corners:
    DO j=js,je+1
      DO i=is-1,ie+1
        ut(i, j) = u(i, j)*dyc(i, j)
      END DO
    END DO
    DO j=js-1,je+1
      DO i=is,ie+1
        vt(i, j) = v(i, j)*dxc(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie+1
        smag_c(i, j) = rarea_c(i, j)*(vt(i, j-1)-vt(i, j)-ut(i-1, j)+ut(&
&         i, j))
      END DO
    END DO
! Fix the corners?? if grid_type /= 4
! Compute shear strain:
    DO j=jsd,jed+1
      DO i=isd,ied
        vt(i, j) = u(i, j)*dx(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied+1
        ut(i, j) = v(i, j)*dy(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied
        wk(i, j) = rarea(i, j)*(vt(i, j)-vt(i, j+1)+ut(i, j)-ut(i+1, j))
      END DO
    END DO
    CALL A2B_ORD4_FWD(wk, sh, gridstruct, npx, npy, is, ie, js, je, &
&                  ng, .false.)
    DO j=js,je+1
      DO i=is,ie+1
        CALL PUSHREALARRAY(smag_c(i, j))
        smag_c(i, j) = dt*SQRT(sh(i, j)**2+smag_c(i, j)**2)
      END DO
    END DO
    CALL PUSHINTEGER(jed)
    CALL PUSHINTEGER(je)
    CALL PUSHINTEGER(is)
    CALL PUSHREALARRAY(sh, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL PUSHINTEGER(ie)
    !CALL PUSHPOINTER8(C_LOC(dyc))
    !CALL PUSHPOINTER8(C_LOC(rarea_c))
    !CALL PUSHPOINTER8(C_LOC(dy))
    !CALL PUSHPOINTER8(C_LOC(dx))
    !CALL PUSHPOINTER8(C_LOC(dxc))
    CALL PUSHINTEGER(js)
  END SUBROUTINE SMAG_CORNER_FWD
!  Differentiation of smag_corner in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_e
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
!   gradient     of useful results: smag_c u v
!   with respect to varying inputs: u v
  SUBROUTINE SMAG_CORNER_BWD(dt, u, u_ad, v, v_ad, ua, va, smag_c, &
&   smag_c_ad, bd, npx, npy, gridstruct, ng)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: npx, npy, ng
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1) :: u_ad
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed) :: v_ad
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: smag_c
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: smag_c_ad
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    REAL :: ut(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: ut_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: vt_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: sh(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: sh_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    INTEGER :: i, j
    INTEGER :: is2, ie1
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc, dx, dy, rarea, rarea_c
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SQRT
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_smag_corner

    ut = 0.0
    vt = 0.0
    wk = 0.0
    sh = 0.0
    is2 = 0
    ie1 = 0
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
    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    dx => gridstruct%dx
    dy => gridstruct%dy
    rarea => gridstruct%rarea
    rarea_c => gridstruct%rarea_c
    CALL POPINTEGER(js)
    !CALL POPPOINTER8(cptr)
    dxc => gridstruct%dxc ! (/unknown_shape_in_smag_corner/))
    !CALL POPPOINTER8(cptr)
    dx => gridstruct%dx ! (/unknown_shape_in_smag_corner/))
    !CALL POPPOINTER8(cptr)
    dy => gridstruct%dy ! (/unknown_shape_in_smag_corner/))
    !CALL POPPOINTER8(cptr)
    rarea_c => gridstruct%rarea_c ! (/unknown_shape_in_smag_corner/))
    !CALL POPPOINTER8(cptr)
    dyc => gridstruct%dyc ! (/unknown_shape_in_smag_corner/))
    CALL POPINTEGER(ie)
    CALL POPREALARRAY(sh, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
    CALL POPINTEGER(is)
    CALL POPINTEGER(je)
    CALL POPINTEGER(jed)
    sh_ad = 0.0
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        CALL POPREALARRAY(smag_c(i, j))
        IF (sh(i, j)**2 + smag_c(i, j)**2 .EQ. 0.0) THEN
          temp_ad1 = 0.0
        ELSE
          temp_ad1 = dt*smag_c_ad(i, j)/(2.0*SQRT(sh(i, j)**2+smag_c(i, &
&           j)**2))
        END IF
        sh_ad(i, j) = sh_ad(i, j) + 2*sh(i, j)*temp_ad1
        smag_c_ad(i, j) = 2*smag_c(i, j)*temp_ad1
      END DO
    END DO
    wk_ad = 0.0
    CALL A2B_ORD4_BWD(wk, wk_ad, sh, sh_ad, gridstruct, npx, npy, is&
&                  , ie, js, je, ng, .false.)
    rarea => gridstruct%rarea
    jsd = bd%jsd
    ied = bd%ied
    isd = bd%isd
    ut_ad = 0.0
    vt_ad = 0.0
    DO j=jed,jsd,-1
      DO i=ied,isd,-1
        temp_ad0 = rarea(i, j)*wk_ad(i, j)
        vt_ad(i, j) = vt_ad(i, j) + temp_ad0
        vt_ad(i, j+1) = vt_ad(i, j+1) - temp_ad0
        ut_ad(i, j) = ut_ad(i, j) + temp_ad0
        ut_ad(i+1, j) = ut_ad(i+1, j) - temp_ad0
        wk_ad(i, j) = 0.0
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ied+1,isd,-1
        v_ad(i, j) = v_ad(i, j) + dy(i, j)*ut_ad(i, j)
        ut_ad(i, j) = 0.0
      END DO
    END DO
    DO j=jed+1,jsd,-1
      DO i=ied,isd,-1
        u_ad(i, j) = u_ad(i, j) + dx(i, j)*vt_ad(i, j)
        vt_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        temp_ad = rarea_c(i, j)*smag_c_ad(i, j)
        vt_ad(i, j-1) = vt_ad(i, j-1) + temp_ad
        vt_ad(i, j) = vt_ad(i, j) - temp_ad
        ut_ad(i, j) = ut_ad(i, j) + temp_ad
        ut_ad(i-1, j) = ut_ad(i-1, j) - temp_ad
        smag_c_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js-1,-1
      DO i=ie+1,is,-1
        v_ad(i, j) = v_ad(i, j) + dxc(i, j)*vt_ad(i, j)
        vt_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is-1,-1
        u_ad(i, j) = u_ad(i, j) + dyc(i, j)*ut_ad(i, j)
        ut_ad(i, j) = 0.0
      END DO
    END DO
  END SUBROUTINE SMAG_CORNER_BWD
 end module sw_core_adm_mod
