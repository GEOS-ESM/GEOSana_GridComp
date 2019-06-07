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
module tp_core_adm_mod
!BOP
!
! !MODULE: tp_core --- A collection of routines to support FV transport
!
 use fv_mp_mod,         only: ng 
 use fv_grid_utils_mod, only: big_number
 use fv_arrays_mod,     only: fv_grid_type, fv_grid_bounds_type, r_grid

 use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                          pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

 implicit none

 private
 public fv_tp_2d, pert_ppm, copy_corners, &
        fv_tp_2d_adm, copy_corners_adm, & 
        fv_tp_2d_fwd, fv_tp_2d_bwd, pert_ppm_adm

 real, parameter:: ppm_fac = 1.5   ! nonlinear scheme limiter: between 1 and 2
 real, parameter:: r3 = 1./3.
 real, parameter:: near_zero = 1.E-25
 real, parameter:: ppm_limiter = 2.0

#ifdef WAVE_FORM
! Suresh & Huynh scheme 2.2 (purtabation form)
! The wave-form is more diffusive than scheme 2.1
 real, parameter:: b1 =   0.0375
 real, parameter:: b2 =  -7./30.
 real, parameter:: b3 =  -23./120.
 real, parameter:: b4 =  13./30.
 real, parameter:: b5 = -11./240.
#else
! scheme 2.1: perturbation form
 real, parameter:: b1 =   1./30.
 real, parameter:: b2 = -13./60.
 real, parameter:: b3 = -13./60.
 real, parameter:: b4 =  0.45
 real, parameter:: b5 = -0.05
#endif
 real, parameter:: t11 = 27./28., t12 = -13./28., t13=3./7.
 real, parameter:: s11 = 11./14., s14 = 4./7.,    s15=3./14.
!----------------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
!----------------------------------------------------
! Non-monotonic
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.
!   q(i+0.5) = p1*(q(i-1)+q(i)) + p2*(q(i-2)+q(i+1))
! integer:: is, ie, js, je, isd, ied, jsd, jed

!---- version number -----
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

!
!EOP
!-----------------------------------------------------------------------

CONTAINS
!  Differentiation of fv_tp_2d in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_c
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
!   gradient     of useful results: xfx q mass mfx mfy ra_x ra_y
!                yfx fx fy crx cry
!   with respect to varying inputs: xfx q mass mfx mfy ra_x ra_y
!                yfx fx fy crx cry
!   q(i+0.5) = p1*(q(i-1)+q(i)) + p2*(q(i-2)+q(i+1))
! integer:: is, ie, js, je, isd, ied, jsd, jed
!
!EOP
!-----------------------------------------------------------------------
  SUBROUTINE FV_TP_2D_ADM(q, q_ad, crx, crx_ad, cry, cry_ad, npx, npy, &
&   hord, fx, fx_ad, fy, fy_ad, xfx, xfx_ad, yfx, yfx_ad, gridstruct, bd&
&   , ra_x, ra_x_ad, ra_y, ra_y_ad, mfx, mfx_ad, mfy, mfy_ad, mass, &
&   mass_ad, nord, damp_c)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL, INTENT(IN) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: crx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed)
!
    REAL, INTENT(IN) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: xfx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed)
!
    REAL, INTENT(IN) :: cry(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: cry_ad(bd%isd:bd%ied, bd%js:bd%je+1)
!
    REAL, INTENT(IN) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: yfx_ad(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(IN) :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_ad(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_ad(bd%isd:bd%ied, bd%js:bd%je)
! transported scalar
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
! Flux in x ( E )
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_ad(bd%is:bd%ie+1, bd%js:bd%je)
! Flux in y ( N )
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_ad(bd%is:bd%ie, bd%js:bd%je+1)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! optional Arguments:
! Mass Flux X-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, OPTIONAL :: mfx_ad(bd%is:bd%ie+1, bd%js:bd%je)
! Mass Flux Y-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, OPTIONAL :: mfy_ad(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL :: mass_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL, INTENT(IN) :: damp_c
    INTEGER, OPTIONAL, INTENT(IN) :: nord
! Local:
    INTEGER :: ord_ou, ord_in
    REAL :: q_i(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: q_i_ad(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: q_j(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: q_j_ad(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: fx2(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: fx2_ad(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: fy2(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fy2_ad(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fyy(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fyy_ad(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fx1(bd%is:bd%ie+1)
    REAL :: fx1_ad(bd%is:bd%ie+1)
    REAL :: damp
    INTEGER :: i, j
    INTEGER :: is, ie, js, je, isd, ied, jsd, jed
    INTRINSIC PRESENT
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    INTEGER :: branch
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (hord .EQ. 10) THEN
      ord_in = 8
    ELSE
      ord_in = hord
    END IF
    ord_ou = hord
    IF (.NOT.gridstruct%nested) THEN
      CALL PUSHREALARRAY_ADM(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL COPY_CORNERS(q, npx, npy, 2, gridstruct%nested, bd, &
&                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct&
&                 %nw_corner, gridstruct%ne_corner)
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL YPPM(fy2, q, cry, ord_in, isd, ied, isd, ied, js, je, jsd, jed&
&       , npx, npy, gridstruct%dya, gridstruct%nested, gridstruct%&
&       grid_type)
    DO j=js,je+1
      DO i=isd,ied
        fyy(i, j) = yfx(i, j)*fy2(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        q_i(i, j) = (q(i, j)*gridstruct%area(i, j)+fyy(i, j)-fyy(i, j+1)&
&         )/ra_y(i, j)
      END DO
    END DO
    CALL PUSHREALARRAY_ADM(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
    CALL XPPM(fx, q_i, crx(is:ie+1, js:je), ord_ou, is, ie, isd, ied, js&
&       , je, jsd, jed, npx, npy, gridstruct%dxa, gridstruct%nested, &
&       gridstruct%grid_type)
    IF (.NOT.gridstruct%nested) THEN
      CALL PUSHREALARRAY_ADM(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL COPY_CORNERS(q, npx, npy, 1, gridstruct%nested, bd, &
&                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct&
&                 %nw_corner, gridstruct%ne_corner)
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL XPPM(fx2, q, crx, ord_in, is, ie, isd, ied, jsd, jed, jsd, jed&
&       , npx, npy, gridstruct%dxa, gridstruct%nested, gridstruct%&
&       grid_type)
    DO j=jsd,jed
      DO i=is,ie+1
        CALL PUSHREALARRAY_ADM(fx1(i))
        fx1(i) = xfx(i, j)*fx2(i, j)
      END DO
      DO i=is,ie
        q_j(i, j) = (q(i, j)*gridstruct%area(i, j)+fx1(i)-fx1(i+1))/ra_x&
&         (i, j)
      END DO
    END DO
    CALL PUSHREALARRAY_ADM(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
    CALL YPPM(fy, q_j, cry, ord_ou, is, ie, isd, ied, js, je, jsd, jed, &
&       npx, npy, gridstruct%dya, gridstruct%nested, gridstruct%&
&       grid_type)
!----------------
! Flux averaging:
!----------------
    IF (PRESENT(mfx) .AND. PRESENT(mfy)) THEN
      IF (PRESENT(nord) .AND. PRESENT(damp_c) .AND. PRESENT(mass)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*gridstruct%da_min)**(nord+1)
          CALL DELN_FLUX_ADM(nord, is, ie, js, je, npx, npy, damp, q, &
&                      q_ad, fx, fx_ad, fy, fy_ad, gridstruct, bd, mass&
&                      , mass_ad)
        END IF
      END IF
      fy2_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie,is,-1
          temp_ad2 = 0.5*mfy(i, j)*fy_ad(i, j)
          fy2_ad(i, j) = fy2_ad(i, j) + temp_ad2
          mfy_ad(i, j) = mfy_ad(i, j) + 0.5*(fy(i, j)+fy2(i, j))*fy_ad(i&
&           , j)
          fy_ad(i, j) = temp_ad2
        END DO
      END DO
      fx2_ad = 0.0
      DO j=je,js,-1
        DO i=ie+1,is,-1
          temp_ad1 = 0.5*mfx(i, j)*fx_ad(i, j)
          fx2_ad(i, j) = fx2_ad(i, j) + temp_ad1
          mfx_ad(i, j) = mfx_ad(i, j) + 0.5*(fx(i, j)+fx2(i, j))*fx_ad(i&
&           , j)
          fx_ad(i, j) = temp_ad1
        END DO
      END DO
    ELSE
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*gridstruct%da_min)**(nord+1)
          CALL DELN_FLUX_ADM(nord, is, ie, js, je, npx, npy, damp, q, &
&                      q_ad, fx, fx_ad, fy, fy_ad, gridstruct, bd)
        END IF
      END IF
      fy2_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie,is,-1
          temp_ad4 = 0.5*yfx(i, j)*fy_ad(i, j)
          fy2_ad(i, j) = fy2_ad(i, j) + temp_ad4
          yfx_ad(i, j) = yfx_ad(i, j) + 0.5*(fy(i, j)+fy2(i, j))*fy_ad(i&
&           , j)
          fy_ad(i, j) = temp_ad4
        END DO
      END DO
      fx2_ad = 0.0
      DO j=je,js,-1
        DO i=ie+1,is,-1
          temp_ad3 = 0.5*xfx(i, j)*fx_ad(i, j)
          fx2_ad(i, j) = fx2_ad(i, j) + temp_ad3
          xfx_ad(i, j) = xfx_ad(i, j) + 0.5*(fx(i, j)+fx2(i, j))*fx_ad(i&
&           , j)
          fx_ad(i, j) = temp_ad3
        END DO
      END DO
    END IF
    CALL POPREALARRAY_ADM(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
    q_j_ad = 0.0
    CALL YPPM_ADM(fy, fy_ad, q_j, q_j_ad, cry, cry_ad, ord_ou, is, ie, &
&           isd, ied, js, je, jsd, jed, npx, npy, gridstruct%dya, &
&           gridstruct%nested, gridstruct%grid_type)
    fx1_ad = 0.0
    DO j=jed,jsd,-1
      DO i=ie,is,-1
        temp_ad0 = q_j_ad(i, j)/ra_x(i, j)
        q_ad(i, j) = q_ad(i, j) + gridstruct%area(i, j)*temp_ad0
        fx1_ad(i) = fx1_ad(i) + temp_ad0
        fx1_ad(i+1) = fx1_ad(i+1) - temp_ad0
        ra_x_ad(i, j) = ra_x_ad(i, j) - (gridstruct%area(i, j)*q(i, j)+&
&         fx1(i)-fx1(i+1))*temp_ad0/ra_x(i, j)
        q_j_ad(i, j) = 0.0
      END DO
      DO i=ie+1,is,-1
        CALL POPREALARRAY_ADM(fx1(i))
        xfx_ad(i, j) = xfx_ad(i, j) + fx2(i, j)*fx1_ad(i)
        fx2_ad(i, j) = fx2_ad(i, j) + xfx(i, j)*fx1_ad(i)
        fx1_ad(i) = 0.0
      END DO
    END DO
    CALL XPPM_ADM(fx2, fx2_ad, q, q_ad, crx, crx_ad, ord_in, is, ie, isd&
&           , ied, jsd, jed, jsd, jed, npx, npy, gridstruct%dxa, &
&           gridstruct%nested, gridstruct%grid_type)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY_ADM(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL COPY_CORNERS_ADM(q, q_ad, npx, npy, 1, gridstruct%nested, bd&
&                     , gridstruct%sw_corner, gridstruct%se_corner, &
&                     gridstruct%nw_corner, gridstruct%ne_corner)
    END IF
    CALL POPREALARRAY_ADM(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
    q_i_ad = 0.0
    CALL XPPM_ADM(fx, fx_ad, q_i, q_i_ad, crx(is:ie+1, js:je), crx_ad(is&
&           :ie+1, js:je), ord_ou, is, ie, isd, ied, js, je, jsd, jed, &
&           npx, npy, gridstruct%dxa, gridstruct%nested, gridstruct%&
&           grid_type)
    fyy_ad = 0.0
    DO j=je,js,-1
      DO i=ied,isd,-1
        temp_ad = q_i_ad(i, j)/ra_y(i, j)
        q_ad(i, j) = q_ad(i, j) + gridstruct%area(i, j)*temp_ad
        fyy_ad(i, j) = fyy_ad(i, j) + temp_ad
        fyy_ad(i, j+1) = fyy_ad(i, j+1) - temp_ad
        ra_y_ad(i, j) = ra_y_ad(i, j) - (gridstruct%area(i, j)*q(i, j)+&
&         fyy(i, j)-fyy(i, j+1))*temp_ad/ra_y(i, j)
        q_i_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ied,isd,-1
        yfx_ad(i, j) = yfx_ad(i, j) + fy2(i, j)*fyy_ad(i, j)
        fy2_ad(i, j) = fy2_ad(i, j) + yfx(i, j)*fyy_ad(i, j)
        fyy_ad(i, j) = 0.0
      END DO
    END DO
    CALL YPPM_ADM(fy2, fy2_ad, q, q_ad, cry, cry_ad, ord_in, isd, ied, &
&           isd, ied, js, je, jsd, jed, npx, npy, gridstruct%dya, &
&           gridstruct%nested, gridstruct%grid_type)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREALARRAY_ADM(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1))
      CALL COPY_CORNERS_ADM(q, q_ad, npx, npy, 2, gridstruct%nested, bd&
&                     , gridstruct%sw_corner, gridstruct%se_corner, &
&                     gridstruct%nw_corner, gridstruct%ne_corner)
    END IF
  END SUBROUTINE FV_TP_2D_ADM
!   q(i+0.5) = p1*(q(i-1)+q(i)) + p2*(q(i-2)+q(i+1))
! integer:: is, ie, js, je, isd, ied, jsd, jed
!
!EOP
!-----------------------------------------------------------------------
  SUBROUTINE FV_TP_2D(q, crx, cry, npx, npy, hord, fx, fy, xfx, yfx, &
&   gridstruct, bd, ra_x, ra_y, mfx, mfy, mass, nord, damp_c)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL, INTENT(IN) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed)
!
    REAL, INTENT(IN) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed)
!
    REAL, INTENT(IN) :: cry(bd%isd:bd%ied, bd%js:bd%je+1)
!
    REAL, INTENT(IN) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(IN) :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
! transported scalar
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
! Flux in x ( E )
    REAL, INTENT(OUT) :: fx(bd%is:bd%ie+1, bd%js:bd%je)
! Flux in y ( N )
    REAL, INTENT(OUT) :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! optional Arguments:
! Mass Flux X-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfx(bd%is:bd%ie+1, bd%js:bd%je)
! Mass Flux Y-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL, INTENT(IN) :: damp_c
    INTEGER, OPTIONAL, INTENT(IN) :: nord
! Local:
    INTEGER :: ord_ou, ord_in
    REAL :: q_i(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: q_j(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: fx2(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: fy2(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fyy(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fx1(bd%is:bd%ie+1)
    REAL :: damp
    INTEGER :: i, j
    INTEGER :: is, ie, js, je, isd, ied, jsd, jed
    INTRINSIC PRESENT
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (hord .EQ. 10) THEN
      ord_in = 8
    ELSE
      ord_in = hord
    END IF
    ord_ou = hord
    IF (.NOT.gridstruct%nested) CALL COPY_CORNERS(q, npx, npy, 2, &
&                                           gridstruct%nested, bd, &
&                                           gridstruct%sw_corner, &
&                                           gridstruct%se_corner, &
&                                           gridstruct%nw_corner, &
&                                           gridstruct%ne_corner)
    CALL YPPM(fy2, q, cry, ord_in, isd, ied, isd, ied, js, je, jsd, jed&
&       , npx, npy, gridstruct%dya, gridstruct%nested, gridstruct%&
&       grid_type)
    DO j=js,je+1
      DO i=isd,ied
        fyy(i, j) = yfx(i, j)*fy2(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        q_i(i, j) = (q(i, j)*gridstruct%area(i, j)+fyy(i, j)-fyy(i, j+1)&
&         )/ra_y(i, j)
      END DO
    END DO
    CALL XPPM(fx, q_i, crx(is:ie+1, js:je), ord_ou, is, ie, isd, ied, js&
&       , je, jsd, jed, npx, npy, gridstruct%dxa, gridstruct%nested, &
&       gridstruct%grid_type)
    IF (.NOT.gridstruct%nested) CALL COPY_CORNERS(q, npx, npy, 1, &
&                                           gridstruct%nested, bd, &
&                                           gridstruct%sw_corner, &
&                                           gridstruct%se_corner, &
&                                           gridstruct%nw_corner, &
&                                           gridstruct%ne_corner)
    CALL XPPM(fx2, q, crx, ord_in, is, ie, isd, ied, jsd, jed, jsd, jed&
&       , npx, npy, gridstruct%dxa, gridstruct%nested, gridstruct%&
&       grid_type)
    DO j=jsd,jed
      DO i=is,ie+1
        fx1(i) = xfx(i, j)*fx2(i, j)
      END DO
      DO i=is,ie
        q_j(i, j) = (q(i, j)*gridstruct%area(i, j)+fx1(i)-fx1(i+1))/ra_x&
&         (i, j)
      END DO
    END DO
    CALL YPPM(fy, q_j, cry, ord_ou, is, ie, isd, ied, js, je, jsd, jed, &
&       npx, npy, gridstruct%dya, gridstruct%nested, gridstruct%&
&       grid_type)
!----------------
! Flux averaging:
!----------------
    IF (PRESENT(mfx) .AND. PRESENT(mfy)) THEN
!---------------------------------
! For transport of pt and tracers
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*mfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*mfy(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c) .AND. PRESENT(mass)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*gridstruct%da_min)**(nord+1)
          CALL DELN_FLUX(nord, is, ie, js, je, npx, npy, damp, q, fx, fy&
&                  , gridstruct, bd, mass)
        END IF
      END IF
    ELSE
!---------------------------------
! For transport of delp, vorticity
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*xfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*yfx(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*gridstruct%da_min)**(nord+1)
          CALL DELN_FLUX(nord, is, ie, js, je, npx, npy, damp, q, fx, fy&
&                  , gridstruct, bd)
        END IF
      END IF
    END IF
  END SUBROUTINE FV_TP_2D
!Weird arguments are because this routine is called in a lot of
!places outside of tp_core, sometimes very deeply nested in the call tree.
  SUBROUTINE COPY_CORNERS(q, npx, npy, dir, nested, bd, sw_corner, &
&   se_corner, nw_corner, ne_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy, dir
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: nested, sw_corner, se_corner, nw_corner, &
&   ne_corner
    INTEGER :: i, j
    IF (nested) THEN
      RETURN
    ELSE IF (dir .EQ. 1) THEN
! XDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            q(i, j) = q(j, 1-i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q(i, j) = q(npy-j, i-npx+1)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q(i, j) = q(j, 2*npx-1-i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q(i, j) = q(npy-j, i-1+npx)
          END DO
        END DO
      END IF
    ELSE IF (dir .EQ. 2) THEN
! YDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            q(i, j) = q(1-j, i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q(i, j) = q(npy+j-1, npx-i)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q(i, j) = q(2*npy-1-j, i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q(i, j) = q(j+1-npx, npy-i)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE COPY_CORNERS
!  Differentiation of xppm in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_core_
!mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn_core
!_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_mod.ge
!opk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynamics_mo
!d.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mod.c
!2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.map
!_scalar_fb fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scalar_pro
!file_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz fv
!_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_setup 
!fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.compute_
!pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_c nh_
!utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver nh_ut
!ils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw_core
!_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core_mod.
!fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.comput
!e_divergence_damping_fb sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corners t
!p_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circle_dis
!t sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: q flux c
!   with respect to varying inputs: q flux c
  SUBROUTINE XPPM_ADM(flux, flux_ad, q, q_ad, c, c_ad, iord, is, ie, isd&
&   , ied, jfirst, jlast, jsd, jed, npx, npy, dxa, nested, grid_type)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, isd, ied, jsd, jed
! compute domain
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: iord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(isd:ied, jfirst:jlast)
    REAL :: q_ad(isd:ied, jfirst:jlast)
! Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(is:ie+1, jfirst:jlast)
    REAL :: c_ad(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: dxa(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
! !OUTPUT PARAMETERS:
!  Flux
    REAL :: flux(is:ie+1, jfirst:jlast)
    REAL :: flux_ad(is:ie+1, jfirst:jlast)
! Local
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    REAL, DIMENSION(is-1:ie+1) :: bl_ad, br_ad, b0_ad
    REAL :: q1(isd:ied)
    REAL :: q1_ad(isd:ied)
    REAL, DIMENSION(is:ie+1) :: fx0, fx1
    REAL, DIMENSION(is:ie+1) :: fx0_ad, fx1_ad
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: al(is-1:ie+2)
    REAL :: al_ad(is-1:ie+2)
    REAL :: dm(is-2:ie+2)
    REAL :: dm_ad(is-2:ie+2)
    REAL :: dq(is-3:ie+2)
    REAL :: dq_ad(is-3:ie+2)
    INTEGER :: i, j, ie3, is1, ie1
    REAL :: x0, x1, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2
    REAL :: xt_ad, qtmp_ad, pmp_1_ad, lac_1_ad, pmp_2_ad, lac_2_ad
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min1_ad
    REAL :: min2
    REAL :: min2_ad
    REAL :: abs0
    REAL :: abs0_ad
    REAL :: abs1
    REAL :: min3
    REAL :: min3_ad
    REAL :: min4
    REAL :: min4_ad
    REAL :: min5
    REAL :: min5_ad
    REAL :: min6
    REAL :: min6_ad
    REAL :: min7
    REAL :: min7_ad
    REAL :: abs2
    REAL :: abs3
    REAL :: abs4
    REAL :: max1
    REAL :: max1_ad
    REAL :: min8
    REAL :: min8_ad
    REAL :: abs5
    REAL :: abs6
    REAL :: abs7
    INTEGER :: arg1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp0
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: x2_ad
    REAL :: y1_ad
    REAL :: x3_ad
    REAL :: y2_ad
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: temp_ad16
    REAL :: x4_ad
    REAL :: y3_ad
    REAL :: z1_ad
    REAL :: x5_ad
    REAL :: y4_ad
    REAL :: x6_ad
    REAL :: y5_ad
    REAL :: x7_ad
    REAL :: y6_ad
    REAL :: x8_ad
    REAL :: y7_ad
    REAL :: x9_ad
    REAL :: y14_ad
    REAL :: y8_ad
    REAL :: x10_ad
    REAL :: y15_ad
    REAL :: y9_ad
    REAL :: temp_ad17
    REAL :: temp_ad18
    REAL :: z2_ad
    REAL :: y10_ad
    REAL :: z3_ad
    REAL :: y11_ad
    REAL :: temp_ad19
    REAL :: temp_ad20
    REAL :: z4_ad
    REAL :: y12_ad
    REAL :: z5_ad
    REAL :: y13_ad
    REAL :: temp1
    REAL :: temp_ad21
    REAL :: temp_ad22
    REAL :: temp_ad23
    REAL :: temp_ad24
    INTEGER :: branch
    REAL :: x10
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
    REAL :: z5
    REAL :: z4
    REAL :: z3
    REAL :: z2
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
    IF (.NOT.nested .AND. grid_type .LT. 3) THEN
      IF (3 .LT. is - 1) THEN
        is1 = is - 1
      ELSE
        is1 = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        ie3 = ie + 2
      ELSE
        ie3 = npx - 2
      END IF
      IF (npx - 3 .GT. ie + 1) THEN
        CALL PUSHCONTROL1B(1)
        ie1 = ie + 1
      ELSE
        CALL PUSHCONTROL1B(1)
        ie1 = npx - 3
      END IF
    ELSE
      CALL PUSHCONTROL1B(0)
      is1 = is - 1
      ie3 = ie + 2
      ie1 = ie + 1
    END IF
    DO j=jfirst,jlast
      DO i=isd,ied
        CALL PUSHREALARRAY_ADM(q1(i))
        q1(i) = q(i, j)
      END DO
      IF (iord .LT. 8 .OR. iord .EQ. 333) THEN
! ord = 2: perfectly linear ppm scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
        DO i=is1,ie3
          CALL PUSHREALARRAY_ADM(al(i))
          al(i) = p1*(q1(i-1)+q1(i)) + p2*(q1(i-2)+q1(i+1))
        END DO
        IF (iord .EQ. 7) THEN
          DO i=is1,ie3
            IF (al(i) .LT. 0.) THEN
              CALL PUSHREALARRAY_ADM(al(i))
              al(i) = 0.5*(q1(i-1)+q1(i))
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            CALL PUSHREALARRAY_ADM(al(0))
            al(0) = c1*q1(-2) + c2*q1(-1) + c3*q1(0)
            CALL PUSHREALARRAY_ADM(al(1))
            al(1) = 0.5*(((2.*dxa(0, j)+dxa(-1, j))*q1(0)-dxa(0, j)*q1(-&
&             1))/(dxa(-1, j)+dxa(0, j))+((2.*dxa(1, j)+dxa(2, j))*q1(1)&
&             -dxa(1, j)*q1(2))/(dxa(1, j)+dxa(2, j)))
            CALL PUSHREALARRAY_ADM(al(2))
            al(2) = c3*q1(1) + c2*q1(2) + c1*q1(3)
            IF (iord .EQ. 7) THEN
              IF (0. .LT. al(0)) THEN
                CALL PUSHREALARRAY_ADM(al(0))
                al(0) = al(0)
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(al(0))
                al(0) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
              IF (0. .LT. al(1)) THEN
                CALL PUSHREALARRAY_ADM(al(1))
                al(1) = al(1)
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(al(1))
                al(1) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
              IF (0. .LT. al(2)) THEN
                CALL PUSHREALARRAY_ADM(al(2))
                al(2) = al(2)
                CALL PUSHCONTROL2B(3)
              ELSE
                CALL PUSHREALARRAY_ADM(al(2))
                al(2) = 0.
                CALL PUSHCONTROL2B(2)
              END IF
            ELSE
              CALL PUSHCONTROL2B(0)
            END IF
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            CALL PUSHREALARRAY_ADM(al(npx-1))
            al(npx-1) = c1*q1(npx-3) + c2*q1(npx-2) + c3*q1(npx-1)
            CALL PUSHREALARRAY_ADM(al(npx))
            al(npx) = 0.5*(((2.*dxa(npx-1, j)+dxa(npx-2, j))*q1(npx-1)-&
&             dxa(npx-1, j)*q1(npx-2))/(dxa(npx-2, j)+dxa(npx-1, j))+((&
&             2.*dxa(npx, j)+dxa(npx+1, j))*q1(npx)-dxa(npx, j)*q1(npx+1&
&             ))/(dxa(npx, j)+dxa(npx+1, j)))
            CALL PUSHREALARRAY_ADM(al(npx+1))
            al(npx+1) = c3*q1(npx) + c2*q1(npx+1) + c1*q1(npx+2)
            IF (iord .EQ. 7) THEN
              IF (0. .LT. al(npx-1)) THEN
                CALL PUSHREALARRAY_ADM(al(npx-1))
                al(npx-1) = al(npx-1)
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(al(npx-1))
                al(npx-1) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
              IF (0. .LT. al(npx)) THEN
                CALL PUSHREALARRAY_ADM(al(npx))
                al(npx) = al(npx)
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(al(npx))
                al(npx) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
              IF (0. .LT. al(npx+1)) THEN
                CALL PUSHREALARRAY_ADM(al(npx+1))
                al(npx+1) = al(npx+1)
                CALL PUSHCONTROL3B(4)
              ELSE
                CALL PUSHREALARRAY_ADM(al(npx+1))
                al(npx+1) = 0.
                CALL PUSHCONTROL3B(3)
              END IF
            ELSE
              CALL PUSHCONTROL3B(0)
            END IF
          ELSE
            CALL PUSHCONTROL3B(1)
          END IF
        ELSE
          CALL PUSHCONTROL3B(2)
        END IF
        IF (iord .EQ. 1) THEN
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
          CALL PUSHCONTROL3B(0)
        ELSE IF (iord .EQ. 2) THEN
! perfectly linear scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            CALL PUSHREALARRAY_ADM(xt)
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              CALL PUSHREALARRAY_ADM(qtmp)
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHREALARRAY_ADM(qtmp)
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
          CALL PUSHCONTROL3B(1)
        ELSE IF (iord .EQ. 333) THEN
!        x0 = sign(dim(xt, 0.), 1.)
!        x1 = sign(dim(0., xt), 1.)
!        flux(i,j) = x0*(q1(i-1)+(1.-xt)*(al(i)-qtmp-xt*(al(i-1)+al(i)-(qtmp+qtmp))))     &
!                  + x1*(q1(i)  +(1.+xt)*(al(i)-qtmp+xt*(al(i)+al(i+1)-(qtmp+qtmp))))
! Perfectly linear scheme, more diffusive than ord=2 (HoldawayKent-2015-TellusA)
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            CALL PUSHREALARRAY_ADM(xt)
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
          CALL PUSHCONTROL3B(2)
        ELSE IF (iord .EQ. 3) THEN
          DO i=is-1,ie+1
            CALL PUSHREALARRAY_ADM(bl(i))
            bl(i) = al(i) - q1(i)
            CALL PUSHREALARRAY_ADM(br(i))
            br(i) = al(i+1) - q1(i)
            CALL PUSHREALARRAY_ADM(b0(i))
            b0(i) = bl(i) + br(i)
            IF (b0(i) .GE. 0.) THEN
              x0 = b0(i)
            ELSE
              x0 = -b0(i)
            END IF
            IF (bl(i) - br(i) .GE. 0.) THEN
              CALL PUSHREALARRAY_ADM(xt)
              xt = bl(i) - br(i)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(xt)
              xt = -(bl(i)-br(i))
              CALL PUSHCONTROL1B(1)
            END IF
            smt5(i) = x0 .LT. xt
            smt6(i) = 3.*x0 .LT. xt
          END DO
          DO i=is,ie+1
            CALL PUSHREALARRAY_ADM(fx1(i))
            fx1(i) = 0.
          END DO
          DO i=is,ie+1
            CALL PUSHREALARRAY_ADM(xt)
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              IF (smt6(i-1) .OR. smt5(i)) THEN
                CALL PUSHREALARRAY_ADM(fx1(i))
                fx1(i) = br(i-1) - xt*b0(i-1)
                CALL PUSHCONTROL3B(5)
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
! 2nd order, piece-wise linear
                CALL PUSHREALARRAY_ADM(fx1(i))
                fx1(i) = SIGN(min1, br(i-1))
                CALL PUSHCONTROL3B(4)
              ELSE
                CALL PUSHCONTROL3B(3)
              END IF
            ELSE IF (smt6(i) .OR. smt5(i-1)) THEN
              CALL PUSHREALARRAY_ADM(fx1(i))
              fx1(i) = bl(i) + xt*b0(i)
              CALL PUSHCONTROL3B(2)
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
              CALL PUSHREALARRAY_ADM(fx1(i))
              fx1(i) = SIGN(min2, bl(i))
              CALL PUSHCONTROL3B(1)
            ELSE
              CALL PUSHCONTROL3B(0)
            END IF
            IF (xt .GE. 0.) THEN
              CALL PUSHREALARRAY_ADM(abs0)
              abs0 = xt
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(abs0)
              abs0 = -xt
              CALL PUSHCONTROL1B(1)
            END IF
          END DO
          CALL PUSHCONTROL3B(3)
        ELSE IF (iord .EQ. 4) THEN
          DO i=is-1,ie+1
            CALL PUSHREALARRAY_ADM(bl(i))
            bl(i) = al(i) - q1(i)
            CALL PUSHREALARRAY_ADM(br(i))
            br(i) = al(i+1) - q1(i)
            CALL PUSHREALARRAY_ADM(b0(i))
            b0(i) = bl(i) + br(i)
            IF (b0(i) .GE. 0.) THEN
              x0 = b0(i)
            ELSE
              x0 = -b0(i)
            END IF
            IF (bl(i) - br(i) .GE. 0.) THEN
              CALL PUSHREALARRAY_ADM(xt)
              xt = bl(i) - br(i)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(xt)
              xt = -(bl(i)-br(i))
              CALL PUSHCONTROL1B(1)
            END IF
            smt5(i) = x0 .LT. xt
            smt6(i) = 3.*x0 .LT. xt
          END DO
          DO i=is,ie+1
            CALL PUSHREALARRAY_ADM(fx1(i))
            fx1(i) = 0.
          END DO
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              IF (smt6(i-1) .OR. smt5(i)) THEN
                CALL PUSHREALARRAY_ADM(fx1(i))
                fx1(i) = (1.-c(i, j))*(br(i-1)-c(i, j)*b0(i-1))
                CALL PUSHCONTROL2B(0)
              ELSE
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (smt6(i) .OR. smt5(i-1)) THEN
              CALL PUSHREALARRAY_ADM(fx1(i))
              fx1(i) = (1.+c(i, j))*(bl(i)+c(i, j)*b0(i))
              CALL PUSHCONTROL2B(2)
            ELSE
              CALL PUSHCONTROL2B(3)
            END IF
          END DO
          CALL PUSHCONTROL3B(4)
        ELSE
! iord = 5 & 6
          IF (iord .EQ. 5) THEN
            DO i=is-1,ie+1
              CALL PUSHREALARRAY_ADM(bl(i))
              bl(i) = al(i) - q1(i)
              CALL PUSHREALARRAY_ADM(br(i))
              br(i) = al(i+1) - q1(i)
              CALL PUSHREALARRAY_ADM(b0(i))
              b0(i) = bl(i) + br(i)
              smt5(i) = bl(i)*br(i) .LT. 0.
            END DO
            CALL PUSHCONTROL1B(1)
          ELSE
            DO i=is-1,ie+1
              CALL PUSHREALARRAY_ADM(bl(i))
              bl(i) = al(i) - q1(i)
              CALL PUSHREALARRAY_ADM(br(i))
              br(i) = al(i+1) - q1(i)
              CALL PUSHREALARRAY_ADM(b0(i))
              b0(i) = bl(i) + br(i)
              IF (3.*b0(i) .GE. 0.) THEN
                abs1 = 3.*b0(i)
              ELSE
                abs1 = -(3.*b0(i))
              END IF
              IF (bl(i) - br(i) .GE. 0.) THEN
                abs4 = bl(i) - br(i)
              ELSE
                abs4 = -(bl(i)-br(i))
              END IF
              smt5(i) = abs1 .LT. abs4
            END DO
            CALL PUSHCONTROL1B(0)
          END IF
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHREALARRAY_ADM(fx1(i))
              fx1(i) = (1.-c(i, j))*(br(i-1)-c(i, j)*b0(i-1))
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(fx1(i))
              fx1(i) = (1.+c(i, j))*(bl(i)+c(i, j)*b0(i))
              CALL PUSHCONTROL1B(1)
            END IF
            IF (smt5(i-1) .OR. smt5(i)) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
          CALL PUSHCONTROL3B(5)
        END IF
      ELSE
! Monotonic constraints:
! ord = 8: PPM with Lin's PPM fast monotone constraint
! ord = 10: PPM with Lin's modification of Huynh 2nd constraint
! ord = 13: 10 plus positive definite constraint
        DO i=is-2,ie+2
          CALL PUSHREALARRAY_ADM(xt)
          xt = 0.25*(q1(i+1)-q1(i-1))
          IF (xt .GE. 0.) THEN
            x4 = xt
            CALL PUSHCONTROL1B(0)
          ELSE
            x4 = -xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q1(i-1) .LT. q1(i)) THEN
            IF (q1(i) .LT. q1(i+1)) THEN
              max1 = q1(i+1)
              CALL PUSHCONTROL2B(0)
            ELSE
              max1 = q1(i)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (q1(i-1) .LT. q1(i+1)) THEN
            max1 = q1(i+1)
            CALL PUSHCONTROL2B(2)
          ELSE
            max1 = q1(i-1)
            CALL PUSHCONTROL2B(3)
          END IF
          y3 = max1 - q1(i)
          IF (q1(i-1) .GT. q1(i)) THEN
            IF (q1(i) .GT. q1(i+1)) THEN
              min8 = q1(i+1)
              CALL PUSHCONTROL2B(0)
            ELSE
              min8 = q1(i)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (q1(i-1) .GT. q1(i+1)) THEN
            min8 = q1(i+1)
            CALL PUSHCONTROL2B(2)
          ELSE
            min8 = q1(i-1)
            CALL PUSHCONTROL2B(3)
          END IF
          z1 = q1(i) - min8
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
        DO i=is1,ie1+1
          CALL PUSHREALARRAY_ADM(al(i))
          al(i) = 0.5*(q1(i-1)+q1(i)) + r3*(dm(i-1)-dm(i))
        END DO
        IF (iord .EQ. 8) THEN
          DO i=is1,ie1
            CALL PUSHREALARRAY_ADM(xt)
            xt = 2.*dm(i)
            IF (xt .GE. 0.) THEN
              x5 = xt
              CALL PUSHCONTROL1B(0)
            ELSE
              x5 = -xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (al(i) - q1(i) .GE. 0.) THEN
              y4 = al(i) - q1(i)
              CALL PUSHCONTROL1B(0)
            ELSE
              y4 = -(al(i)-q1(i))
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
            IF (al(i+1) - q1(i) .GE. 0.) THEN
              y5 = al(i+1) - q1(i)
              CALL PUSHCONTROL1B(0)
            ELSE
              y5 = -(al(i+1)-q1(i))
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
        ELSE IF (iord .EQ. 11) THEN
! This is emulation of 2nd van Leer scheme using PPM codes
          DO i=is1,ie1
            CALL PUSHREALARRAY_ADM(xt)
            xt = ppm_fac*dm(i)
            IF (xt .GE. 0.) THEN
              x7 = xt
              CALL PUSHCONTROL1B(0)
            ELSE
              x7 = -xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (al(i) - q1(i) .GE. 0.) THEN
              y6 = al(i) - q1(i)
              CALL PUSHCONTROL1B(0)
            ELSE
              y6 = -(al(i)-q1(i))
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x7 .GT. y6) THEN
              CALL PUSHREALARRAY_ADM(min6)
              min6 = y6
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(min6)
              min6 = x7
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREALARRAY_ADM(bl(i))
            bl(i) = -SIGN(min6, xt)
            IF (xt .GE. 0.) THEN
              x8 = xt
              CALL PUSHCONTROL1B(0)
            ELSE
              x8 = -xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (al(i+1) - q1(i) .GE. 0.) THEN
              y7 = al(i+1) - q1(i)
              CALL PUSHCONTROL1B(0)
            ELSE
              y7 = -(al(i+1)-q1(i))
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x8 .GT. y7) THEN
              CALL PUSHREALARRAY_ADM(min7)
              min7 = y7
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(min7)
              min7 = x8
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREALARRAY_ADM(br(i))
            br(i) = SIGN(min7, xt)
          END DO
          CALL PUSHCONTROL2B(1)
        ELSE
          DO i=is1-2,ie1+1
            dq(i) = 2.*(q1(i+1)-q1(i))
          END DO
          DO i=is1,ie1
            CALL PUSHREALARRAY_ADM(bl(i))
            bl(i) = al(i) - q1(i)
            CALL PUSHREALARRAY_ADM(br(i))
            br(i) = al(i+1) - q1(i)
            IF (dm(i-1) .GE. 0.) THEN
              abs2 = dm(i-1)
            ELSE
              abs2 = -dm(i-1)
            END IF
            IF (dm(i) .GE. 0.) THEN
              abs5 = dm(i)
            ELSE
              abs5 = -dm(i)
            END IF
            IF (dm(i+1) .GE. 0.) THEN
              abs7 = dm(i+1)
            ELSE
              abs7 = -dm(i+1)
            END IF
            IF (abs2 + abs5 + abs7 .LT. near_zero) THEN
              bl(i) = 0.
              br(i) = 0.
              CALL PUSHCONTROL2B(3)
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
                pmp_2 = dq(i-1)
                lac_2 = pmp_2 - 0.75*dq(i-2)
                IF (0. .LT. pmp_2) THEN
                  IF (pmp_2 .LT. lac_2) THEN
                    x9 = lac_2
                    CALL PUSHCONTROL2B(0)
                  ELSE
                    x9 = pmp_2
                    CALL PUSHCONTROL2B(1)
                  END IF
                ELSE IF (0. .LT. lac_2) THEN
                  x9 = lac_2
                  CALL PUSHCONTROL2B(2)
                ELSE
                  CALL PUSHCONTROL2B(3)
                  x9 = 0.
                END IF
                IF (0. .GT. pmp_2) THEN
                  IF (pmp_2 .GT. lac_2) THEN
                    y14 = lac_2
                    CALL PUSHCONTROL2B(0)
                  ELSE
                    y14 = pmp_2
                    CALL PUSHCONTROL2B(1)
                  END IF
                ELSE IF (0. .GT. lac_2) THEN
                  y14 = lac_2
                  CALL PUSHCONTROL2B(2)
                ELSE
                  y14 = 0.
                  CALL PUSHCONTROL2B(3)
                END IF
                IF (br(i) .LT. y14) THEN
                  y8 = y14
                  CALL PUSHCONTROL1B(0)
                ELSE
                  y8 = br(i)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (x9 .GT. y8) THEN
                  br(i) = y8
                  CALL PUSHCONTROL1B(0)
                ELSE
                  br(i) = x9
                  CALL PUSHCONTROL1B(1)
                END IF
                pmp_1 = -dq(i)
                lac_1 = pmp_1 + 0.75*dq(i+1)
                IF (0. .LT. pmp_1) THEN
                  IF (pmp_1 .LT. lac_1) THEN
                    x10 = lac_1
                    CALL PUSHCONTROL2B(0)
                  ELSE
                    x10 = pmp_1
                    CALL PUSHCONTROL2B(1)
                  END IF
                ELSE IF (0. .LT. lac_1) THEN
                  x10 = lac_1
                  CALL PUSHCONTROL2B(2)
                ELSE
                  CALL PUSHCONTROL2B(3)
                  x10 = 0.
                END IF
                IF (0. .GT. pmp_1) THEN
                  IF (pmp_1 .GT. lac_1) THEN
                    y15 = lac_1
                    CALL PUSHCONTROL2B(0)
                  ELSE
                    y15 = pmp_1
                    CALL PUSHCONTROL2B(1)
                  END IF
                ELSE IF (0. .GT. lac_1) THEN
                  y15 = lac_1
                  CALL PUSHCONTROL2B(2)
                ELSE
                  y15 = 0.
                  CALL PUSHCONTROL2B(3)
                END IF
                IF (bl(i) .LT. y15) THEN
                  y9 = y15
                  CALL PUSHCONTROL1B(0)
                ELSE
                  y9 = bl(i)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (x10 .GT. y9) THEN
                  bl(i) = y9
                  CALL PUSHCONTROL2B(1)
                ELSE
                  bl(i) = x10
                  CALL PUSHCONTROL2B(2)
                END IF
              ELSE
                CALL PUSHCONTROL2B(0)
              END IF
            END IF
          END DO
          CALL PUSHCONTROL2B(2)
        END IF
! Positive definite constraint:
        IF (iord .EQ. 9 .OR. iord .EQ. 13) THEN
          arg1 = ie1 - is1 + 1
          CALL PUSHREALARRAY_ADM(br(is1:ie1), ie1 - is1 + 1)
          CALL PUSHREALARRAY_ADM(bl(is1:ie1), ie1 - is1 + 1)
          CALL PERT_PPM(arg1, q1(is1:ie1), bl(is1:ie1), br(is1:ie1), 0)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            CALL PUSHREALARRAY_ADM(bl(0))
            bl(0) = s14*dm(-1) + s11*(q1(-1)-q1(0))
            CALL PUSHREALARRAY_ADM(xt)
            xt = 0.5*(((2.*dxa(0, j)+dxa(-1, j))*q1(0)-dxa(0, j)*q1(-1))&
&             /(dxa(-1, j)+dxa(0, j))+((2.*dxa(1, j)+dxa(2, j))*q1(1)-&
&             dxa(1, j)*q1(2))/(dxa(1, j)+dxa(2, j)))
            IF (q1(1) .GT. q1(2)) THEN
              z2 = q1(2)
              CALL PUSHCONTROL1B(0)
            ELSE
              z2 = q1(1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q1(-1) .GT. q1(0)) THEN
              IF (q1(0) .GT. z2) THEN
                y10 = z2
                CALL PUSHCONTROL2B(0)
              ELSE
                y10 = q1(0)
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (q1(-1) .GT. z2) THEN
              y10 = z2
              CALL PUSHCONTROL2B(2)
            ELSE
              y10 = q1(-1)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (xt .LT. y10) THEN
              xt = y10
              CALL PUSHCONTROL1B(0)
            ELSE
              xt = xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q1(1) .LT. q1(2)) THEN
              z3 = q1(2)
              CALL PUSHCONTROL1B(0)
            ELSE
              z3 = q1(1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q1(-1) .LT. q1(0)) THEN
              IF (q1(0) .LT. z3) THEN
                y11 = z3
                CALL PUSHCONTROL2B(0)
              ELSE
                y11 = q1(0)
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (q1(-1) .LT. z3) THEN
              y11 = z3
              CALL PUSHCONTROL2B(2)
            ELSE
              y11 = q1(-1)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (xt .GT. y11) THEN
              xt = y11
              CALL PUSHCONTROL1B(0)
            ELSE
              xt = xt
              CALL PUSHCONTROL1B(1)
            END IF
!        endif
            CALL PUSHREALARRAY_ADM(br(0))
            br(0) = xt - q1(0)
            CALL PUSHREALARRAY_ADM(bl(1))
            bl(1) = xt - q1(1)
            xt = s15*q1(1) + s11*q1(2) - s14*dm(2)
            CALL PUSHREALARRAY_ADM(br(1))
            br(1) = xt - q1(1)
            CALL PUSHREALARRAY_ADM(bl(2))
            bl(2) = xt - q1(2)
            CALL PUSHREALARRAY_ADM(br(2))
            br(2) = al(3) - q1(2)
            CALL PUSHREALARRAY_ADM(br(0:2), 3)
            CALL PUSHREALARRAY_ADM(bl(0:2), 3)
            CALL PERT_PPM(3, q1(0:2), bl(0:2), br(0:2), 1)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            CALL PUSHREALARRAY_ADM(bl(npx-2))
            bl(npx-2) = al(npx-2) - q1(npx-2)
            CALL PUSHREALARRAY_ADM(xt)
            xt = s15*q1(npx-1) + s11*q1(npx-2) + s14*dm(npx-2)
            CALL PUSHREALARRAY_ADM(br(npx-2))
            br(npx-2) = xt - q1(npx-2)
            CALL PUSHREALARRAY_ADM(bl(npx-1))
            bl(npx-1) = xt - q1(npx-1)
            xt = 0.5*(((2.*dxa(npx-1, j)+dxa(npx-2, j))*q1(npx-1)-dxa(&
&             npx-1, j)*q1(npx-2))/(dxa(npx-2, j)+dxa(npx-1, j))+((2.*&
&             dxa(npx, j)+dxa(npx+1, j))*q1(npx)-dxa(npx, j)*q1(npx+1))/&
&             (dxa(npx, j)+dxa(npx+1, j)))
            IF (q1(npx) .GT. q1(npx+1)) THEN
              z4 = q1(npx+1)
              CALL PUSHCONTROL1B(0)
            ELSE
              z4 = q1(npx)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q1(npx-2) .GT. q1(npx-1)) THEN
              IF (q1(npx-1) .GT. z4) THEN
                y12 = z4
                CALL PUSHCONTROL2B(0)
              ELSE
                y12 = q1(npx-1)
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (q1(npx-2) .GT. z4) THEN
              y12 = z4
              CALL PUSHCONTROL2B(2)
            ELSE
              y12 = q1(npx-2)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (xt .LT. y12) THEN
              xt = y12
              CALL PUSHCONTROL1B(0)
            ELSE
              xt = xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q1(npx) .LT. q1(npx+1)) THEN
              z5 = q1(npx+1)
              CALL PUSHCONTROL1B(0)
            ELSE
              z5 = q1(npx)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q1(npx-2) .LT. q1(npx-1)) THEN
              IF (q1(npx-1) .LT. z5) THEN
                y13 = z5
                CALL PUSHCONTROL2B(0)
              ELSE
                y13 = q1(npx-1)
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (q1(npx-2) .LT. z5) THEN
              y13 = z5
              CALL PUSHCONTROL2B(2)
            ELSE
              y13 = q1(npx-2)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (xt .GT. y13) THEN
              xt = y13
              CALL PUSHCONTROL1B(0)
            ELSE
              xt = xt
              CALL PUSHCONTROL1B(1)
            END IF
!        endif
            CALL PUSHREALARRAY_ADM(br(npx-1))
            br(npx-1) = xt - q1(npx-1)
            CALL PUSHREALARRAY_ADM(bl(npx))
            bl(npx) = xt - q1(npx)
            CALL PUSHREALARRAY_ADM(br(npx))
            br(npx) = s11*(q1(npx+1)-q1(npx)) - s14*dm(npx+1)
            CALL PUSHREALARRAY_ADM(br(npx-2:npx), 3)
            CALL PUSHREALARRAY_ADM(bl(npx-2:npx), 3)
            CALL PERT_PPM(3, q1(npx-2:npx), bl(npx-2:npx), br(npx-2:npx)&
&                   , 1)
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(0)
        END IF
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL3B(6)
      END IF
    END DO
    dm_ad = 0.0
    dq_ad = 0.0
    al_ad = 0.0
    bl_ad = 0.0
    q1_ad = 0.0
    br_ad = 0.0
    b0_ad = 0.0
    fx0_ad = 0.0
    fx1_ad = 0.0
    DO j=jlast,jfirst,-1
      CALL POPCONTROL3B(branch)
      IF (branch .LT. 3) THEN
        IF (branch .EQ. 0) THEN
          DO i=ie+1,is,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q1_ad(i) = q1_ad(i) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              q1_ad(i-1) = q1_ad(i-1) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            END IF
          END DO
        ELSE IF (branch .EQ. 1) THEN
          DO i=ie+1,is,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt = c(i, j)
              qtmp = q1(i)
              temp0 = al(i) + al(i+1) - 2*qtmp
              temp_ad5 = (xt+1.)*flux_ad(i, j)
              temp_ad6 = xt*temp_ad5
              qtmp_ad = flux_ad(i, j) - temp_ad5 - 2*temp_ad6
              xt_ad = temp0*temp_ad5 + (al(i)-qtmp+xt*temp0)*flux_ad(i, &
&               j)
              al_ad(i) = al_ad(i) + temp_ad6 + temp_ad5
              al_ad(i+1) = al_ad(i+1) + temp_ad6
              flux_ad(i, j) = 0.0
              CALL POPREALARRAY_ADM(qtmp)
              q1_ad(i) = q1_ad(i) + qtmp_ad
            ELSE
              xt = c(i, j)
              qtmp = q1(i-1)
              temp = al(i-1) + al(i) - 2*qtmp
              temp_ad3 = (1.-xt)*flux_ad(i, j)
              temp_ad4 = -(xt*temp_ad3)
              qtmp_ad = flux_ad(i, j) - temp_ad3 - 2*temp_ad4
              xt_ad = -(temp*temp_ad3) - (al(i)-qtmp-xt*temp)*flux_ad(i&
&               , j)
              al_ad(i) = al_ad(i) + temp_ad4 + temp_ad3
              al_ad(i-1) = al_ad(i-1) + temp_ad4
              flux_ad(i, j) = 0.0
              CALL POPREALARRAY_ADM(qtmp)
              q1_ad(i-1) = q1_ad(i-1) + qtmp_ad
            END IF
            CALL POPREALARRAY_ADM(xt)
            c_ad(i, j) = c_ad(i, j) + xt_ad
          END DO
        ELSE
          DO i=ie+1,is,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt = c(i, j)
              temp_ad10 = flux_ad(i, j)/6.0
              temp_ad11 = -(0.5*xt*flux_ad(i, j))
              temp_ad12 = xt**2*flux_ad(i, j)/6.0
              q1_ad(i-1) = q1_ad(i-1) + temp_ad12 - temp_ad11 + 2.0*&
&               temp_ad10
              q1_ad(i) = q1_ad(i) + temp_ad11 - 2.0*temp_ad12 + 5.0*&
&               temp_ad10
              q1_ad(i+1) = q1_ad(i+1) + temp_ad12 - temp_ad10
              xt_ad = ((q1(i+1)-2.0*q1(i)+q1(i-1))*2*xt/6.0-0.5*(q1(i)-&
&               q1(i-1)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              xt = c(i, j)
              temp_ad7 = flux_ad(i, j)/6.0
              temp_ad8 = -(0.5*xt*flux_ad(i, j))
              temp_ad9 = xt**2*flux_ad(i, j)/6.0
              q1_ad(i) = q1_ad(i) + temp_ad9 + temp_ad8 + 2.0*temp_ad7
              q1_ad(i-1) = q1_ad(i-1) + 5.0*temp_ad7 - temp_ad8 - 2.0*&
&               temp_ad9
              q1_ad(i-2) = q1_ad(i-2) + temp_ad9 - temp_ad7
              xt_ad = ((q1(i)-2.0*q1(i-1)+q1(i-2))*2*xt/6.0-0.5*(q1(i)-&
&               q1(i-1)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0
            END IF
            CALL POPREALARRAY_ADM(xt)
            c_ad(i, j) = c_ad(i, j) + xt_ad
          END DO
        END IF
      ELSE IF (branch .LT. 5) THEN
        IF (branch .EQ. 3) THEN
          DO i=ie+1,is,-1
            fx0_ad(i) = fx0_ad(i) + flux_ad(i, j)
            abs0_ad = -(fx1(i)*flux_ad(i, j))
            fx1_ad(i) = fx1_ad(i) + (1.-abs0)*flux_ad(i, j)
            flux_ad(i, j) = 0.0
            xt = c(i, j)
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(abs0)
              xt_ad = abs0_ad
            ELSE
              CALL POPREALARRAY_ADM(abs0)
              xt_ad = -abs0_ad
            END IF
            CALL POPCONTROL3B(branch)
            IF (branch .LT. 3) THEN
              IF (branch .NE. 0) THEN
                IF (branch .EQ. 1) THEN
                  CALL POPREALARRAY_ADM(fx1(i))
                  min2_ad = SIGN(1.d0, min2*bl(i))*fx1_ad(i)
                  fx1_ad(i) = 0.0
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
                ELSE
                  CALL POPREALARRAY_ADM(fx1(i))
                  bl_ad(i) = bl_ad(i) + fx1_ad(i)
                  xt_ad = xt_ad + b0(i)*fx1_ad(i)
                  b0_ad(i) = b0_ad(i) + xt*fx1_ad(i)
                  fx1_ad(i) = 0.0
                END IF
              END IF
              q1_ad(i) = q1_ad(i) + fx0_ad(i)
              fx0_ad(i) = 0.0
            ELSE
              IF (branch .NE. 3) THEN
                IF (branch .EQ. 4) THEN
                  CALL POPREALARRAY_ADM(fx1(i))
                  min1_ad = SIGN(1.d0, min1*br(i-1))*fx1_ad(i)
                  fx1_ad(i) = 0.0
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
                ELSE
                  CALL POPREALARRAY_ADM(fx1(i))
                  br_ad(i-1) = br_ad(i-1) + fx1_ad(i)
                  xt_ad = xt_ad - b0(i-1)*fx1_ad(i)
                  b0_ad(i-1) = b0_ad(i-1) - xt*fx1_ad(i)
                  fx1_ad(i) = 0.0
                END IF
              END IF
              q1_ad(i-1) = q1_ad(i-1) + fx0_ad(i)
              fx0_ad(i) = 0.0
            END IF
            CALL POPREALARRAY_ADM(xt)
            c_ad(i, j) = c_ad(i, j) + xt_ad
          END DO
          DO i=ie+1,is,-1
            CALL POPREALARRAY_ADM(fx1(i))
            fx1_ad(i) = 0.0
          END DO
          DO i=ie+1,is-1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(xt)
            ELSE
              CALL POPREALARRAY_ADM(xt)
            END IF
            CALL POPREALARRAY_ADM(b0(i))
            bl_ad(i) = bl_ad(i) + b0_ad(i)
            br_ad(i) = br_ad(i) + b0_ad(i)
            b0_ad(i) = 0.0
            CALL POPREALARRAY_ADM(br(i))
            al_ad(i+1) = al_ad(i+1) + br_ad(i)
            q1_ad(i) = q1_ad(i) - bl_ad(i) - br_ad(i)
            br_ad(i) = 0.0
            CALL POPREALARRAY_ADM(bl(i))
            al_ad(i) = al_ad(i) + bl_ad(i)
            bl_ad(i) = 0.0
          END DO
        ELSE
          DO i=ie+1,is,-1
            fx0_ad(i) = fx0_ad(i) + flux_ad(i, j)
            fx1_ad(i) = fx1_ad(i) + flux_ad(i, j)
            flux_ad(i, j) = 0.0
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY_ADM(fx1(i))
                temp_ad13 = (1.-c(i, j))*fx1_ad(i)
                c_ad(i, j) = c_ad(i, j) - b0(i-1)*temp_ad13 - (br(i-1)-c&
&                 (i, j)*b0(i-1))*fx1_ad(i)
                br_ad(i-1) = br_ad(i-1) + temp_ad13
                b0_ad(i-1) = b0_ad(i-1) - c(i, j)*temp_ad13
                fx1_ad(i) = 0.0
              END IF
              q1_ad(i-1) = q1_ad(i-1) + fx0_ad(i)
              fx0_ad(i) = 0.0
            ELSE
              IF (branch .EQ. 2) THEN
                CALL POPREALARRAY_ADM(fx1(i))
                temp_ad14 = (c(i, j)+1.)*fx1_ad(i)
                c_ad(i, j) = c_ad(i, j) + b0(i)*temp_ad14 + (bl(i)+c(i, &
&                 j)*b0(i))*fx1_ad(i)
                bl_ad(i) = bl_ad(i) + temp_ad14
                b0_ad(i) = b0_ad(i) + c(i, j)*temp_ad14
                fx1_ad(i) = 0.0
              END IF
              q1_ad(i) = q1_ad(i) + fx0_ad(i)
              fx0_ad(i) = 0.0
            END IF
          END DO
          DO i=ie+1,is,-1
            CALL POPREALARRAY_ADM(fx1(i))
            fx1_ad(i) = 0.0
          END DO
          DO i=ie+1,is-1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(xt)
            ELSE
              CALL POPREALARRAY_ADM(xt)
            END IF
            CALL POPREALARRAY_ADM(b0(i))
            bl_ad(i) = bl_ad(i) + b0_ad(i)
            br_ad(i) = br_ad(i) + b0_ad(i)
            b0_ad(i) = 0.0
            CALL POPREALARRAY_ADM(br(i))
            al_ad(i+1) = al_ad(i+1) + br_ad(i)
            q1_ad(i) = q1_ad(i) - bl_ad(i) - br_ad(i)
            br_ad(i) = 0.0
            CALL POPREALARRAY_ADM(bl(i))
            al_ad(i) = al_ad(i) + bl_ad(i)
            bl_ad(i) = 0.0
          END DO
        END IF
      ELSE IF (branch .EQ. 5) THEN
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) fx1_ad(i) = fx1_ad(i) + flux_ad(i, j)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q1_ad(i-1) = q1_ad(i-1) + flux_ad(i, j)
            flux_ad(i, j) = 0.0
            CALL POPREALARRAY_ADM(fx1(i))
            temp_ad15 = (1.-c(i, j))*fx1_ad(i)
            c_ad(i, j) = c_ad(i, j) - b0(i-1)*temp_ad15 - (br(i-1)-c(i, &
&             j)*b0(i-1))*fx1_ad(i)
            br_ad(i-1) = br_ad(i-1) + temp_ad15
            b0_ad(i-1) = b0_ad(i-1) - c(i, j)*temp_ad15
            fx1_ad(i) = 0.0
          ELSE
            q1_ad(i) = q1_ad(i) + flux_ad(i, j)
            flux_ad(i, j) = 0.0
            CALL POPREALARRAY_ADM(fx1(i))
            temp_ad16 = (c(i, j)+1.)*fx1_ad(i)
            c_ad(i, j) = c_ad(i, j) + b0(i)*temp_ad16 + (bl(i)+c(i, j)*&
&             b0(i))*fx1_ad(i)
            bl_ad(i) = bl_ad(i) + temp_ad16
            b0_ad(i) = b0_ad(i) + c(i, j)*temp_ad16
            fx1_ad(i) = 0.0
          END IF
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO i=ie+1,is-1,-1
            CALL POPREALARRAY_ADM(b0(i))
            bl_ad(i) = bl_ad(i) + b0_ad(i)
            br_ad(i) = br_ad(i) + b0_ad(i)
            b0_ad(i) = 0.0
            CALL POPREALARRAY_ADM(br(i))
            al_ad(i+1) = al_ad(i+1) + br_ad(i)
            q1_ad(i) = q1_ad(i) - bl_ad(i) - br_ad(i)
            br_ad(i) = 0.0
            CALL POPREALARRAY_ADM(bl(i))
            al_ad(i) = al_ad(i) + bl_ad(i)
            bl_ad(i) = 0.0
          END DO
        ELSE
          DO i=ie+1,is-1,-1
            CALL POPREALARRAY_ADM(b0(i))
            bl_ad(i) = bl_ad(i) + b0_ad(i)
            br_ad(i) = br_ad(i) + b0_ad(i)
            b0_ad(i) = 0.0
            CALL POPREALARRAY_ADM(br(i))
            al_ad(i+1) = al_ad(i+1) + br_ad(i)
            q1_ad(i) = q1_ad(i) - bl_ad(i) - br_ad(i)
            br_ad(i) = 0.0
            CALL POPREALARRAY_ADM(bl(i))
            al_ad(i) = al_ad(i) + bl_ad(i)
            bl_ad(i) = 0.0
          END DO
        END IF
      ELSE
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            temp_ad23 = (c(i, j)+1.)*flux_ad(i, j)
            temp_ad24 = c(i, j)*temp_ad23
            q1_ad(i) = q1_ad(i) + flux_ad(i, j)
            c_ad(i, j) = c_ad(i, j) + (bl(i)+br(i))*temp_ad23 + (bl(i)+c&
&             (i, j)*(bl(i)+br(i)))*flux_ad(i, j)
            bl_ad(i) = bl_ad(i) + temp_ad24 + temp_ad23
            br_ad(i) = br_ad(i) + temp_ad24
            flux_ad(i, j) = 0.0
          ELSE
            temp1 = bl(i-1) + br(i-1)
            temp_ad21 = (1.-c(i, j))*flux_ad(i, j)
            temp_ad22 = -(c(i, j)*temp_ad21)
            q1_ad(i-1) = q1_ad(i-1) + flux_ad(i, j)
            c_ad(i, j) = c_ad(i, j) - temp1*temp_ad21 - (br(i-1)-c(i, j)&
&             *temp1)*flux_ad(i, j)
            br_ad(i-1) = br_ad(i-1) + temp_ad22 + temp_ad21
            bl_ad(i-1) = bl_ad(i-1) + temp_ad22
            flux_ad(i, j) = 0.0
          END IF
        END DO
        CALL POPCONTROL2B(branch)
        IF (branch .NE. 0) THEN
          IF (branch .NE. 1) THEN
            CALL POPREALARRAY_ADM(bl(npx-2:npx), 3)
            CALL POPREALARRAY_ADM(br(npx-2:npx), 3)
            CALL PERT_PPM_ADM(3, q1(npx-2:npx), bl(npx-2:npx), bl_ad(npx&
&                       -2:npx), br(npx-2:npx), br_ad(npx-2:npx), 1)
            CALL POPREALARRAY_ADM(br(npx))
            q1_ad(npx+1) = q1_ad(npx+1) + s11*br_ad(npx)
            q1_ad(npx) = q1_ad(npx) - bl_ad(npx) - s11*br_ad(npx)
            dm_ad(npx+1) = dm_ad(npx+1) - s14*br_ad(npx)
            br_ad(npx) = 0.0
            CALL POPREALARRAY_ADM(bl(npx))
            xt_ad = br_ad(npx-1) + bl_ad(npx)
            bl_ad(npx) = 0.0
            CALL POPREALARRAY_ADM(br(npx-1))
            q1_ad(npx-1) = q1_ad(npx-1) - br_ad(npx-1)
            br_ad(npx-1) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y13_ad = xt_ad
              xt_ad = 0.0
            ELSE
              y13_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                z5_ad = y13_ad
              ELSE
                q1_ad(npx-1) = q1_ad(npx-1) + y13_ad
                z5_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              z5_ad = y13_ad
            ELSE
              q1_ad(npx-2) = q1_ad(npx-2) + y13_ad
              z5_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q1_ad(npx+1) = q1_ad(npx+1) + z5_ad
            ELSE
              q1_ad(npx) = q1_ad(npx) + z5_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y12_ad = xt_ad
              xt_ad = 0.0
            ELSE
              y12_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                z4_ad = y12_ad
              ELSE
                q1_ad(npx-1) = q1_ad(npx-1) + y12_ad
                z4_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              z4_ad = y12_ad
            ELSE
              q1_ad(npx-2) = q1_ad(npx-2) + y12_ad
              z4_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q1_ad(npx+1) = q1_ad(npx+1) + z4_ad
            ELSE
              q1_ad(npx) = q1_ad(npx) + z4_ad
            END IF
            temp_ad19 = 0.5*xt_ad/(dxa(npx-2, j)+dxa(npx-1, j))
            temp_ad20 = 0.5*xt_ad/(dxa(npx, j)+dxa(npx+1, j))
            q1_ad(npx-1) = q1_ad(npx-1) + (dxa(npx-1, j)*2.+dxa(npx-2, j&
&             ))*temp_ad19
            q1_ad(npx-2) = q1_ad(npx-2) - dxa(npx-1, j)*temp_ad19
            q1_ad(npx) = q1_ad(npx) + (dxa(npx, j)*2.+dxa(npx+1, j))*&
&             temp_ad20
            q1_ad(npx+1) = q1_ad(npx+1) - dxa(npx, j)*temp_ad20
            CALL POPREALARRAY_ADM(bl(npx-1))
            xt_ad = br_ad(npx-2) + bl_ad(npx-1)
            q1_ad(npx-1) = q1_ad(npx-1) - bl_ad(npx-1)
            bl_ad(npx-1) = 0.0
            CALL POPREALARRAY_ADM(br(npx-2))
            q1_ad(npx-2) = q1_ad(npx-2) - br_ad(npx-2)
            br_ad(npx-2) = 0.0
            CALL POPREALARRAY_ADM(xt)
            q1_ad(npx-1) = q1_ad(npx-1) + s15*xt_ad
            q1_ad(npx-2) = q1_ad(npx-2) + s11*xt_ad - bl_ad(npx-2)
            dm_ad(npx-2) = dm_ad(npx-2) + s14*xt_ad
            CALL POPREALARRAY_ADM(bl(npx-2))
            al_ad(npx-2) = al_ad(npx-2) + bl_ad(npx-2)
            bl_ad(npx-2) = 0.0
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY_ADM(bl(0:2), 3)
            CALL POPREALARRAY_ADM(br(0:2), 3)
            CALL PERT_PPM_ADM(3, q1(0:2), bl(0:2), bl_ad(0:2), br(0:2), &
&                       br_ad(0:2), 1)
            CALL POPREALARRAY_ADM(br(2))
            al_ad(3) = al_ad(3) + br_ad(2)
            q1_ad(2) = q1_ad(2) - bl_ad(2) - br_ad(2)
            br_ad(2) = 0.0
            CALL POPREALARRAY_ADM(bl(2))
            xt_ad = br_ad(1) + bl_ad(2)
            bl_ad(2) = 0.0
            CALL POPREALARRAY_ADM(br(1))
            q1_ad(1) = q1_ad(1) + s15*xt_ad - br_ad(1)
            br_ad(1) = 0.0
            q1_ad(2) = q1_ad(2) + s11*xt_ad
            dm_ad(2) = dm_ad(2) - s14*xt_ad
            CALL POPREALARRAY_ADM(bl(1))
            xt_ad = br_ad(0) + bl_ad(1)
            q1_ad(1) = q1_ad(1) - bl_ad(1)
            bl_ad(1) = 0.0
            CALL POPREALARRAY_ADM(br(0))
            q1_ad(0) = q1_ad(0) - br_ad(0)
            br_ad(0) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y11_ad = xt_ad
              xt_ad = 0.0
            ELSE
              y11_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                z3_ad = y11_ad
              ELSE
                q1_ad(0) = q1_ad(0) + y11_ad
                z3_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              z3_ad = y11_ad
            ELSE
              q1_ad(-1) = q1_ad(-1) + y11_ad
              z3_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q1_ad(2) = q1_ad(2) + z3_ad
            ELSE
              q1_ad(1) = q1_ad(1) + z3_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y10_ad = xt_ad
              xt_ad = 0.0
            ELSE
              y10_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                z2_ad = y10_ad
              ELSE
                q1_ad(0) = q1_ad(0) + y10_ad
                z2_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              z2_ad = y10_ad
            ELSE
              q1_ad(-1) = q1_ad(-1) + y10_ad
              z2_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q1_ad(2) = q1_ad(2) + z2_ad
            ELSE
              q1_ad(1) = q1_ad(1) + z2_ad
            END IF
            CALL POPREALARRAY_ADM(xt)
            temp_ad17 = 0.5*xt_ad/(dxa(-1, j)+dxa(0, j))
            temp_ad18 = 0.5*xt_ad/(dxa(1, j)+dxa(2, j))
            q1_ad(0) = q1_ad(0) + (dxa(0, j)*2.+dxa(-1, j))*temp_ad17
            q1_ad(-1) = q1_ad(-1) - dxa(0, j)*temp_ad17
            q1_ad(1) = q1_ad(1) + (dxa(1, j)*2.+dxa(2, j))*temp_ad18
            q1_ad(2) = q1_ad(2) - dxa(1, j)*temp_ad18
            CALL POPREALARRAY_ADM(bl(0))
            dm_ad(-1) = dm_ad(-1) + s14*bl_ad(0)
            q1_ad(-1) = q1_ad(-1) + s11*bl_ad(0)
            q1_ad(0) = q1_ad(0) - s11*bl_ad(0)
            bl_ad(0) = 0.0
          END IF
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          arg1 = ie1 - is1 + 1
          CALL POPREALARRAY_ADM(bl(is1:ie1), ie1 - is1 + 1)
          CALL POPREALARRAY_ADM(br(is1:ie1), ie1 - is1 + 1)
          CALL PERT_PPM_ADM(arg1, q1(is1:ie1), bl(is1:ie1), bl_ad(is1:&
&                     ie1), br(is1:ie1), br_ad(is1:ie1), 0)
        END IF
        CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          DO i=ie1,is1,-1
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
              q1_ad(i) = q1_ad(i) - y5_ad
            ELSE
              q1_ad(i) = q1_ad(i) + y5_ad
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
              q1_ad(i) = q1_ad(i) - y4_ad
            ELSE
              q1_ad(i) = q1_ad(i) + y4_ad
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
        ELSE IF (branch .EQ. 1) THEN
          DO i=ie1,is1,-1
            CALL POPREALARRAY_ADM(br(i))
            min7_ad = SIGN(1.d0, min7*xt)*br_ad(i)
            br_ad(i) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(min7)
              y7_ad = min7_ad
              x8_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min7)
              x8_ad = min7_ad
              y7_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              al_ad(i+1) = al_ad(i+1) + y7_ad
              q1_ad(i) = q1_ad(i) - y7_ad
            ELSE
              q1_ad(i) = q1_ad(i) + y7_ad
              al_ad(i+1) = al_ad(i+1) - y7_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt_ad = x8_ad
            ELSE
              xt_ad = -x8_ad
            END IF
            CALL POPREALARRAY_ADM(bl(i))
            min6_ad = -(SIGN(1.d0, min6*xt)*bl_ad(i))
            bl_ad(i) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(min6)
              y6_ad = min6_ad
              x7_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min6)
              x7_ad = min6_ad
              y6_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              al_ad(i) = al_ad(i) + y6_ad
              q1_ad(i) = q1_ad(i) - y6_ad
            ELSE
              q1_ad(i) = q1_ad(i) + y6_ad
              al_ad(i) = al_ad(i) - y6_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt_ad = xt_ad + x7_ad
            ELSE
              xt_ad = xt_ad - x7_ad
            END IF
            CALL POPREALARRAY_ADM(xt)
            dm_ad(i) = dm_ad(i) + ppm_fac*xt_ad
          END DO
        ELSE
          DO i=ie1,is1,-1
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                GOTO 100
              ELSE
                y9_ad = bl_ad(i)
                bl_ad(i) = 0.0
                x10_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              x10_ad = bl_ad(i)
              bl_ad(i) = 0.0
              y9_ad = 0.0
            ELSE
              br_ad(i) = 0.0
              bl_ad(i) = 0.0
              GOTO 100
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y15_ad = y9_ad
            ELSE
              bl_ad(i) = bl_ad(i) + y9_ad
              y15_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_1_ad = y15_ad
                pmp_1_ad = 0.0
              ELSE
                pmp_1_ad = y15_ad
                lac_1_ad = 0.0
              END IF
            ELSE
              IF (branch .EQ. 2) THEN
                lac_1_ad = y15_ad
              ELSE
                lac_1_ad = 0.0
              END IF
              pmp_1_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_1_ad = lac_1_ad + x10_ad
              ELSE
                pmp_1_ad = pmp_1_ad + x10_ad
              END IF
            ELSE IF (branch .EQ. 2) THEN
              lac_1_ad = lac_1_ad + x10_ad
            END IF
            pmp_1_ad = pmp_1_ad + lac_1_ad
            dq_ad(i+1) = dq_ad(i+1) + 0.75*lac_1_ad
            dq_ad(i) = dq_ad(i) - pmp_1_ad
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y8_ad = br_ad(i)
              br_ad(i) = 0.0
              x9_ad = 0.0
            ELSE
              x9_ad = br_ad(i)
              br_ad(i) = 0.0
              y8_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y14_ad = y8_ad
            ELSE
              br_ad(i) = br_ad(i) + y8_ad
              y14_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_2_ad = y14_ad
                pmp_2_ad = 0.0
              ELSE
                pmp_2_ad = y14_ad
                lac_2_ad = 0.0
              END IF
            ELSE
              IF (branch .EQ. 2) THEN
                lac_2_ad = y14_ad
              ELSE
                lac_2_ad = 0.0
              END IF
              pmp_2_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_2_ad = lac_2_ad + x9_ad
              ELSE
                pmp_2_ad = pmp_2_ad + x9_ad
              END IF
            ELSE IF (branch .EQ. 2) THEN
              lac_2_ad = lac_2_ad + x9_ad
            END IF
            pmp_2_ad = pmp_2_ad + lac_2_ad
            dq_ad(i-2) = dq_ad(i-2) - 0.75*lac_2_ad
            dq_ad(i-1) = dq_ad(i-1) + pmp_2_ad
 100        CALL POPREALARRAY_ADM(br(i))
            al_ad(i+1) = al_ad(i+1) + br_ad(i)
            q1_ad(i) = q1_ad(i) - bl_ad(i) - br_ad(i)
            br_ad(i) = 0.0
            CALL POPREALARRAY_ADM(bl(i))
            al_ad(i) = al_ad(i) + bl_ad(i)
            bl_ad(i) = 0.0
          END DO
          DO i=ie1+1,is1-2,-1
            q1_ad(i+1) = q1_ad(i+1) + 2.*dq_ad(i)
            q1_ad(i) = q1_ad(i) - 2.*dq_ad(i)
            dq_ad(i) = 0.0
          END DO
        END IF
        DO i=ie1+1,is1,-1
          CALL POPREALARRAY_ADM(al(i))
          q1_ad(i-1) = q1_ad(i-1) + 0.5*al_ad(i)
          q1_ad(i) = q1_ad(i) + 0.5*al_ad(i)
          dm_ad(i-1) = dm_ad(i-1) + r3*al_ad(i)
          dm_ad(i) = dm_ad(i) - r3*al_ad(i)
          al_ad(i) = 0.0
        END DO
        DO i=ie+2,is-2,-1
          xt = 0.25*(q1(i+1)-q1(i-1))
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
          q1_ad(i) = q1_ad(i) + z1_ad
          min8_ad = -z1_ad
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              q1_ad(i+1) = q1_ad(i+1) + min8_ad
            ELSE
              q1_ad(i) = q1_ad(i) + min8_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            q1_ad(i+1) = q1_ad(i+1) + min8_ad
          ELSE
            q1_ad(i-1) = q1_ad(i-1) + min8_ad
          END IF
          max1_ad = y3_ad
          q1_ad(i) = q1_ad(i) - y3_ad
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              q1_ad(i+1) = q1_ad(i+1) + max1_ad
            ELSE
              q1_ad(i) = q1_ad(i) + max1_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            q1_ad(i+1) = q1_ad(i+1) + max1_ad
          ELSE
            q1_ad(i-1) = q1_ad(i-1) + max1_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            xt_ad = x4_ad
          ELSE
            xt_ad = -x4_ad
          END IF
          CALL POPREALARRAY_ADM(xt)
          q1_ad(i+1) = q1_ad(i+1) + 0.25*xt_ad
          q1_ad(i-1) = q1_ad(i-1) - 0.25*xt_ad
        END DO
        GOTO 130
      END IF
      CALL POPCONTROL3B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .NE. 0) GOTO 110
      ELSE IF (branch .EQ. 2) THEN
        GOTO 120
      ELSE
        IF (branch .EQ. 3) THEN
          CALL POPREALARRAY_ADM(al(npx+1))
          al_ad(npx+1) = 0.0
        ELSE
          CALL POPREALARRAY_ADM(al(npx+1))
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY_ADM(al(npx))
        ELSE
          CALL POPREALARRAY_ADM(al(npx))
          al_ad(npx) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY_ADM(al(npx-1))
        ELSE
          CALL POPREALARRAY_ADM(al(npx-1))
          al_ad(npx-1) = 0.0
        END IF
      END IF
      CALL POPREALARRAY_ADM(al(npx+1))
      q1_ad(npx) = q1_ad(npx) + c3*al_ad(npx+1)
      q1_ad(npx+1) = q1_ad(npx+1) + c2*al_ad(npx+1)
      q1_ad(npx+2) = q1_ad(npx+2) + c1*al_ad(npx+1)
      al_ad(npx+1) = 0.0
      CALL POPREALARRAY_ADM(al(npx))
      temp_ad1 = 0.5*al_ad(npx)/(dxa(npx-2, j)+dxa(npx-1, j))
      temp_ad2 = 0.5*al_ad(npx)/(dxa(npx, j)+dxa(npx+1, j))
      q1_ad(npx-1) = q1_ad(npx-1) + (dxa(npx-1, j)*2.+dxa(npx-2, j))*&
&       temp_ad1
      q1_ad(npx-2) = q1_ad(npx-2) - dxa(npx-1, j)*temp_ad1
      q1_ad(npx) = q1_ad(npx) + (dxa(npx, j)*2.+dxa(npx+1, j))*temp_ad2
      q1_ad(npx+1) = q1_ad(npx+1) - dxa(npx, j)*temp_ad2
      al_ad(npx) = 0.0
      CALL POPREALARRAY_ADM(al(npx-1))
      q1_ad(npx-3) = q1_ad(npx-3) + c1*al_ad(npx-1)
      q1_ad(npx-2) = q1_ad(npx-2) + c2*al_ad(npx-1)
      q1_ad(npx-1) = q1_ad(npx-1) + c3*al_ad(npx-1)
      al_ad(npx-1) = 0.0
 110  CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .NE. 0) GOTO 120
      ELSE
        IF (branch .EQ. 2) THEN
          CALL POPREALARRAY_ADM(al(2))
          al_ad(2) = 0.0
        ELSE
          CALL POPREALARRAY_ADM(al(2))
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY_ADM(al(1))
        ELSE
          CALL POPREALARRAY_ADM(al(1))
          al_ad(1) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREALARRAY_ADM(al(0))
        ELSE
          CALL POPREALARRAY_ADM(al(0))
          al_ad(0) = 0.0
        END IF
      END IF
      CALL POPREALARRAY_ADM(al(2))
      q1_ad(1) = q1_ad(1) + c3*al_ad(2)
      q1_ad(2) = q1_ad(2) + c2*al_ad(2)
      q1_ad(3) = q1_ad(3) + c1*al_ad(2)
      al_ad(2) = 0.0
      CALL POPREALARRAY_ADM(al(1))
      temp_ad = 0.5*al_ad(1)/(dxa(-1, j)+dxa(0, j))
      temp_ad0 = 0.5*al_ad(1)/(dxa(1, j)+dxa(2, j))
      q1_ad(0) = q1_ad(0) + (dxa(0, j)*2.+dxa(-1, j))*temp_ad
      q1_ad(-1) = q1_ad(-1) - dxa(0, j)*temp_ad
      q1_ad(1) = q1_ad(1) + (dxa(1, j)*2.+dxa(2, j))*temp_ad0
      q1_ad(2) = q1_ad(2) - dxa(1, j)*temp_ad0
      al_ad(1) = 0.0
      CALL POPREALARRAY_ADM(al(0))
      q1_ad(-2) = q1_ad(-2) + c1*al_ad(0)
      q1_ad(-1) = q1_ad(-1) + c2*al_ad(0)
      q1_ad(0) = q1_ad(0) + c3*al_ad(0)
      al_ad(0) = 0.0
 120  CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO i=ie3,is1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            CALL POPREALARRAY_ADM(al(i))
            q1_ad(i-1) = q1_ad(i-1) + 0.5*al_ad(i)
            q1_ad(i) = q1_ad(i) + 0.5*al_ad(i)
            al_ad(i) = 0.0
          END IF
        END DO
      END IF
      DO i=ie3,is1,-1
        CALL POPREALARRAY_ADM(al(i))
        q1_ad(i-1) = q1_ad(i-1) + p1*al_ad(i)
        q1_ad(i) = q1_ad(i) + p1*al_ad(i)
        q1_ad(i-2) = q1_ad(i-2) + p2*al_ad(i)
        q1_ad(i+1) = q1_ad(i+1) + p2*al_ad(i)
        al_ad(i) = 0.0
      END DO
 130  DO i=ied,isd,-1
        CALL POPREALARRAY_ADM(q1(i))
        q_ad(i, j) = q_ad(i, j) + q1_ad(i)
        q1_ad(i) = 0.0
      END DO
    END DO
    CALL POPCONTROL1B(branch)
  END SUBROUTINE XPPM_ADM
  SUBROUTINE XPPM(flux, q, c, iord, is, ie, isd, ied, jfirst, jlast, jsd&
&   , jed, npx, npy, dxa, nested, grid_type)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, isd, ied, jsd, jed
! compute domain
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: iord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(isd:ied, jfirst:jlast)
! Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: dxa(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
! !OUTPUT PARAMETERS:
!  Flux
    REAL, INTENT(OUT) :: flux(is:ie+1, jfirst:jlast)
! Local
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    REAL :: q1(isd:ied)
    REAL, DIMENSION(is:ie+1) :: fx0, fx1
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: al(is-1:ie+2)
    REAL :: dm(is-2:ie+2)
    REAL :: dq(is-3:ie+2)
    INTEGER :: i, j, ie3, is1, ie1
    REAL :: x0, x1, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min2
    REAL :: abs0
    REAL :: abs1
    REAL :: min3
    REAL :: min4
    REAL :: min5
    REAL :: min6
    REAL :: min7
    REAL :: abs2
    REAL :: abs3
    REAL :: abs4
    REAL :: max1
    REAL :: min8
    REAL :: abs5
    REAL :: abs6
    REAL :: abs7
    INTEGER :: arg1
    REAL :: x10
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
    REAL :: z5
    REAL :: z4
    REAL :: z3
    REAL :: z2
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
    IF (.NOT.nested .AND. grid_type .LT. 3) THEN
      IF (3 .LT. is - 1) THEN
        is1 = is - 1
      ELSE
        is1 = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        ie3 = ie + 2
      ELSE
        ie3 = npx - 2
      END IF
      IF (npx - 3 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 3
      END IF
    ELSE
      is1 = is - 1
      ie3 = ie + 2
      ie1 = ie + 1
    END IF
    DO j=jfirst,jlast
      DO i=isd,ied
        q1(i) = q(i, j)
      END DO
      IF (iord .LT. 8 .OR. iord .EQ. 333) THEN
! ord = 2: perfectly linear ppm scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
        DO i=is1,ie3
          al(i) = p1*(q1(i-1)+q1(i)) + p2*(q1(i-2)+q1(i+1))
        END DO
        IF (iord .EQ. 7) THEN
          DO i=is1,ie3
            IF (al(i) .LT. 0.) al(i) = 0.5*(q1(i-1)+q1(i))
          END DO
        END IF
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            al(0) = c1*q1(-2) + c2*q1(-1) + c3*q1(0)
            al(1) = 0.5*(((2.*dxa(0, j)+dxa(-1, j))*q1(0)-dxa(0, j)*q1(-&
&             1))/(dxa(-1, j)+dxa(0, j))+((2.*dxa(1, j)+dxa(2, j))*q1(1)&
&             -dxa(1, j)*q1(2))/(dxa(1, j)+dxa(2, j)))
            al(2) = c3*q1(1) + c2*q1(2) + c1*q1(3)
            IF (iord .EQ. 7) THEN
              IF (0. .LT. al(0)) THEN
                al(0) = al(0)
              ELSE
                al(0) = 0.
              END IF
              IF (0. .LT. al(1)) THEN
                al(1) = al(1)
              ELSE
                al(1) = 0.
              END IF
              IF (0. .LT. al(2)) THEN
                al(2) = al(2)
              ELSE
                al(2) = 0.
              END IF
            END IF
          END IF
          IF (ie + 1 .EQ. npx) THEN
            al(npx-1) = c1*q1(npx-3) + c2*q1(npx-2) + c3*q1(npx-1)
            al(npx) = 0.5*(((2.*dxa(npx-1, j)+dxa(npx-2, j))*q1(npx-1)-&
&             dxa(npx-1, j)*q1(npx-2))/(dxa(npx-2, j)+dxa(npx-1, j))+((&
&             2.*dxa(npx, j)+dxa(npx+1, j))*q1(npx)-dxa(npx, j)*q1(npx+1&
&             ))/(dxa(npx, j)+dxa(npx+1, j)))
            al(npx+1) = c3*q1(npx) + c2*q1(npx+1) + c1*q1(npx+2)
            IF (iord .EQ. 7) THEN
              IF (0. .LT. al(npx-1)) THEN
                al(npx-1) = al(npx-1)
              ELSE
                al(npx-1) = 0.
              END IF
              IF (0. .LT. al(npx)) THEN
                al(npx) = al(npx)
              ELSE
                al(npx) = 0.
              END IF
              IF (0. .LT. al(npx+1)) THEN
                al(npx+1) = al(npx+1)
              ELSE
                al(npx+1) = 0.
              END IF
            END IF
          END IF
        END IF
        IF (iord .EQ. 1) THEN
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              flux(i, j) = q1(i-1)
            ELSE
              flux(i, j) = q1(i)
            END IF
          END DO
        ELSE IF (iord .EQ. 2) THEN
! perfectly linear scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              qtmp = q1(i-1)
              flux(i, j) = qtmp + (1.-xt)*(al(i)-qtmp-xt*(al(i-1)+al(i)-&
&               (qtmp+qtmp)))
            ELSE
              qtmp = q1(i)
              flux(i, j) = qtmp + (1.+xt)*(al(i)-qtmp+xt*(al(i)+al(i+1)-&
&               (qtmp+qtmp)))
            END IF
          END DO
        ELSE IF (iord .EQ. 333) THEN
!        x0 = sign(dim(xt, 0.), 1.)
!        x1 = sign(dim(0., xt), 1.)
!        flux(i,j) = x0*(q1(i-1)+(1.-xt)*(al(i)-qtmp-xt*(al(i-1)+al(i)-(qtmp+qtmp))))     &
!                  + x1*(q1(i)  +(1.+xt)*(al(i)-qtmp+xt*(al(i)+al(i+1)-(qtmp+qtmp))))
! Perfectly linear scheme, more diffusive than ord=2 (HoldawayKent-2015-TellusA)
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              flux(i, j) = (2.0*q1(i)+5.0*q1(i-1)-q1(i-2))/6.0 - 0.5*xt*&
&               (q1(i)-q1(i-1)) + xt*xt/6.0*(q1(i)-2.0*q1(i-1)+q1(i-2))
            ELSE
              flux(i, j) = (2.0*q1(i-1)+5.0*q1(i)-q1(i+1))/6.0 - 0.5*xt*&
&               (q1(i)-q1(i-1)) + xt*xt/6.0*(q1(i+1)-2.0*q1(i)+q1(i-1))
            END IF
          END DO
        ELSE IF (iord .EQ. 3) THEN
          DO i=is-1,ie+1
            bl(i) = al(i) - q1(i)
            br(i) = al(i+1) - q1(i)
            b0(i) = bl(i) + br(i)
            IF (b0(i) .GE. 0.) THEN
              x0 = b0(i)
            ELSE
              x0 = -b0(i)
            END IF
            IF (bl(i) - br(i) .GE. 0.) THEN
              xt = bl(i) - br(i)
            ELSE
              xt = -(bl(i)-br(i))
            END IF
            smt5(i) = x0 .LT. xt
            smt6(i) = 3.*x0 .LT. xt
          END DO
          DO i=is,ie+1
            fx1(i) = 0.
          END DO
          DO i=is,ie+1
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              fx0(i) = q1(i-1)
              IF (smt6(i-1) .OR. smt5(i)) THEN
                fx1(i) = br(i-1) - xt*b0(i-1)
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
! 2nd order, piece-wise linear
                fx1(i) = SIGN(min1, br(i-1))
              END IF
            ELSE
              fx0(i) = q1(i)
              IF (smt6(i) .OR. smt5(i-1)) THEN
                fx1(i) = bl(i) + xt*b0(i)
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
                fx1(i) = SIGN(min2, bl(i))
              END IF
            END IF
            IF (xt .GE. 0.) THEN
              abs0 = xt
            ELSE
              abs0 = -xt
            END IF
            flux(i, j) = fx0(i) + (1.-abs0)*fx1(i)
          END DO
        ELSE IF (iord .EQ. 4) THEN
          DO i=is-1,ie+1
            bl(i) = al(i) - q1(i)
            br(i) = al(i+1) - q1(i)
            b0(i) = bl(i) + br(i)
            IF (b0(i) .GE. 0.) THEN
              x0 = b0(i)
            ELSE
              x0 = -b0(i)
            END IF
            IF (bl(i) - br(i) .GE. 0.) THEN
              xt = bl(i) - br(i)
            ELSE
              xt = -(bl(i)-br(i))
            END IF
            smt5(i) = x0 .LT. xt
            smt6(i) = 3.*x0 .LT. xt
          END DO
          DO i=is,ie+1
            fx1(i) = 0.
          END DO
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              fx0(i) = q1(i-1)
              IF (smt6(i-1) .OR. smt5(i)) fx1(i) = (1.-c(i, j))*(br(i-1)&
&                 -c(i, j)*b0(i-1))
            ELSE
              fx0(i) = q1(i)
              IF (smt6(i) .OR. smt5(i-1)) fx1(i) = (1.+c(i, j))*(bl(i)+c&
&                 (i, j)*b0(i))
            END IF
            flux(i, j) = fx0(i) + fx1(i)
          END DO
        ELSE
! iord = 5 & 6
          IF (iord .EQ. 5) THEN
            DO i=is-1,ie+1
              bl(i) = al(i) - q1(i)
              br(i) = al(i+1) - q1(i)
              b0(i) = bl(i) + br(i)
              smt5(i) = bl(i)*br(i) .LT. 0.
            END DO
          ELSE
            DO i=is-1,ie+1
              bl(i) = al(i) - q1(i)
              br(i) = al(i+1) - q1(i)
              b0(i) = bl(i) + br(i)
              IF (3.*b0(i) .GE. 0.) THEN
                abs1 = 3.*b0(i)
              ELSE
                abs1 = -(3.*b0(i))
              END IF
              IF (bl(i) - br(i) .GE. 0.) THEN
                abs4 = bl(i) - br(i)
              ELSE
                abs4 = -(bl(i)-br(i))
              END IF
              smt5(i) = abs1 .LT. abs4
            END DO
          END IF
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              fx1(i) = (1.-c(i, j))*(br(i-1)-c(i, j)*b0(i-1))
              flux(i, j) = q1(i-1)
            ELSE
              fx1(i) = (1.+c(i, j))*(bl(i)+c(i, j)*b0(i))
              flux(i, j) = q1(i)
            END IF
            IF (smt5(i-1) .OR. smt5(i)) flux(i, j) = flux(i, j) + fx1(i)
          END DO
        END IF
      ELSE
! Monotonic constraints:
! ord = 8: PPM with Lin's PPM fast monotone constraint
! ord = 10: PPM with Lin's modification of Huynh 2nd constraint
! ord = 13: 10 plus positive definite constraint
        DO i=is-2,ie+2
          xt = 0.25*(q1(i+1)-q1(i-1))
          IF (xt .GE. 0.) THEN
            x4 = xt
          ELSE
            x4 = -xt
          END IF
          IF (q1(i-1) .LT. q1(i)) THEN
            IF (q1(i) .LT. q1(i+1)) THEN
              max1 = q1(i+1)
            ELSE
              max1 = q1(i)
            END IF
          ELSE IF (q1(i-1) .LT. q1(i+1)) THEN
            max1 = q1(i+1)
          ELSE
            max1 = q1(i-1)
          END IF
          y3 = max1 - q1(i)
          IF (q1(i-1) .GT. q1(i)) THEN
            IF (q1(i) .GT. q1(i+1)) THEN
              min8 = q1(i+1)
            ELSE
              min8 = q1(i)
            END IF
          ELSE IF (q1(i-1) .GT. q1(i+1)) THEN
            min8 = q1(i+1)
          ELSE
            min8 = q1(i-1)
          END IF
          z1 = q1(i) - min8
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
        DO i=is1,ie1+1
          al(i) = 0.5*(q1(i-1)+q1(i)) + r3*(dm(i-1)-dm(i))
        END DO
        IF (iord .EQ. 8) THEN
          DO i=is1,ie1
            xt = 2.*dm(i)
            IF (xt .GE. 0.) THEN
              x5 = xt
            ELSE
              x5 = -xt
            END IF
            IF (al(i) - q1(i) .GE. 0.) THEN
              y4 = al(i) - q1(i)
            ELSE
              y4 = -(al(i)-q1(i))
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
            IF (al(i+1) - q1(i) .GE. 0.) THEN
              y5 = al(i+1) - q1(i)
            ELSE
              y5 = -(al(i+1)-q1(i))
            END IF
            IF (x6 .GT. y5) THEN
              min5 = y5
            ELSE
              min5 = x6
            END IF
            br(i) = SIGN(min5, xt)
          END DO
        ELSE IF (iord .EQ. 11) THEN
! This is emulation of 2nd van Leer scheme using PPM codes
          DO i=is1,ie1
            xt = ppm_fac*dm(i)
            IF (xt .GE. 0.) THEN
              x7 = xt
            ELSE
              x7 = -xt
            END IF
            IF (al(i) - q1(i) .GE. 0.) THEN
              y6 = al(i) - q1(i)
            ELSE
              y6 = -(al(i)-q1(i))
            END IF
            IF (x7 .GT. y6) THEN
              min6 = y6
            ELSE
              min6 = x7
            END IF
            bl(i) = -SIGN(min6, xt)
            IF (xt .GE. 0.) THEN
              x8 = xt
            ELSE
              x8 = -xt
            END IF
            IF (al(i+1) - q1(i) .GE. 0.) THEN
              y7 = al(i+1) - q1(i)
            ELSE
              y7 = -(al(i+1)-q1(i))
            END IF
            IF (x8 .GT. y7) THEN
              min7 = y7
            ELSE
              min7 = x8
            END IF
            br(i) = SIGN(min7, xt)
          END DO
        ELSE
          DO i=is1-2,ie1+1
            dq(i) = 2.*(q1(i+1)-q1(i))
          END DO
          DO i=is1,ie1
            bl(i) = al(i) - q1(i)
            br(i) = al(i+1) - q1(i)
            IF (dm(i-1) .GE. 0.) THEN
              abs2 = dm(i-1)
            ELSE
              abs2 = -dm(i-1)
            END IF
            IF (dm(i) .GE. 0.) THEN
              abs5 = dm(i)
            ELSE
              abs5 = -dm(i)
            END IF
            IF (dm(i+1) .GE. 0.) THEN
              abs7 = dm(i+1)
            ELSE
              abs7 = -dm(i+1)
            END IF
            IF (abs2 + abs5 + abs7 .LT. near_zero) THEN
              bl(i) = 0.
              br(i) = 0.
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
                pmp_2 = dq(i-1)
                lac_2 = pmp_2 - 0.75*dq(i-2)
                IF (0. .LT. pmp_2) THEN
                  IF (pmp_2 .LT. lac_2) THEN
                    x9 = lac_2
                  ELSE
                    x9 = pmp_2
                  END IF
                ELSE IF (0. .LT. lac_2) THEN
                  x9 = lac_2
                ELSE
                  x9 = 0.
                END IF
                IF (0. .GT. pmp_2) THEN
                  IF (pmp_2 .GT. lac_2) THEN
                    y14 = lac_2
                  ELSE
                    y14 = pmp_2
                  END IF
                ELSE IF (0. .GT. lac_2) THEN
                  y14 = lac_2
                ELSE
                  y14 = 0.
                END IF
                IF (br(i) .LT. y14) THEN
                  y8 = y14
                ELSE
                  y8 = br(i)
                END IF
                IF (x9 .GT. y8) THEN
                  br(i) = y8
                ELSE
                  br(i) = x9
                END IF
                pmp_1 = -dq(i)
                lac_1 = pmp_1 + 0.75*dq(i+1)
                IF (0. .LT. pmp_1) THEN
                  IF (pmp_1 .LT. lac_1) THEN
                    x10 = lac_1
                  ELSE
                    x10 = pmp_1
                  END IF
                ELSE IF (0. .LT. lac_1) THEN
                  x10 = lac_1
                ELSE
                  x10 = 0.
                END IF
                IF (0. .GT. pmp_1) THEN
                  IF (pmp_1 .GT. lac_1) THEN
                    y15 = lac_1
                  ELSE
                    y15 = pmp_1
                  END IF
                ELSE IF (0. .GT. lac_1) THEN
                  y15 = lac_1
                ELSE
                  y15 = 0.
                END IF
                IF (bl(i) .LT. y15) THEN
                  y9 = y15
                ELSE
                  y9 = bl(i)
                END IF
                IF (x10 .GT. y9) THEN
                  bl(i) = y9
                ELSE
                  bl(i) = x10
                END IF
              END IF
            END IF
          END DO
        END IF
! Positive definite constraint:
        IF (iord .EQ. 9 .OR. iord .EQ. 13) THEN
          arg1 = ie1 - is1 + 1
          CALL PERT_PPM(arg1, q1(is1:ie1), bl(is1:ie1), br(is1:ie1), 0)
        END IF
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            bl(0) = s14*dm(-1) + s11*(q1(-1)-q1(0))
            xt = 0.5*(((2.*dxa(0, j)+dxa(-1, j))*q1(0)-dxa(0, j)*q1(-1))&
&             /(dxa(-1, j)+dxa(0, j))+((2.*dxa(1, j)+dxa(2, j))*q1(1)-&
&             dxa(1, j)*q1(2))/(dxa(1, j)+dxa(2, j)))
            IF (q1(1) .GT. q1(2)) THEN
              z2 = q1(2)
            ELSE
              z2 = q1(1)
            END IF
            IF (q1(-1) .GT. q1(0)) THEN
              IF (q1(0) .GT. z2) THEN
                y10 = z2
              ELSE
                y10 = q1(0)
              END IF
            ELSE IF (q1(-1) .GT. z2) THEN
              y10 = z2
            ELSE
              y10 = q1(-1)
            END IF
            IF (xt .LT. y10) THEN
              xt = y10
            ELSE
              xt = xt
            END IF
            IF (q1(1) .LT. q1(2)) THEN
              z3 = q1(2)
            ELSE
              z3 = q1(1)
            END IF
            IF (q1(-1) .LT. q1(0)) THEN
              IF (q1(0) .LT. z3) THEN
                y11 = z3
              ELSE
                y11 = q1(0)
              END IF
            ELSE IF (q1(-1) .LT. z3) THEN
              y11 = z3
            ELSE
              y11 = q1(-1)
            END IF
            IF (xt .GT. y11) THEN
              xt = y11
            ELSE
              xt = xt
            END IF
!        endif
            br(0) = xt - q1(0)
            bl(1) = xt - q1(1)
            xt = s15*q1(1) + s11*q1(2) - s14*dm(2)
            br(1) = xt - q1(1)
            bl(2) = xt - q1(2)
            br(2) = al(3) - q1(2)
            CALL PERT_PPM(3, q1(0:2), bl(0:2), br(0:2), 1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            bl(npx-2) = al(npx-2) - q1(npx-2)
            xt = s15*q1(npx-1) + s11*q1(npx-2) + s14*dm(npx-2)
            br(npx-2) = xt - q1(npx-2)
            bl(npx-1) = xt - q1(npx-1)
            xt = 0.5*(((2.*dxa(npx-1, j)+dxa(npx-2, j))*q1(npx-1)-dxa(&
&             npx-1, j)*q1(npx-2))/(dxa(npx-2, j)+dxa(npx-1, j))+((2.*&
&             dxa(npx, j)+dxa(npx+1, j))*q1(npx)-dxa(npx, j)*q1(npx+1))/&
&             (dxa(npx, j)+dxa(npx+1, j)))
            IF (q1(npx) .GT. q1(npx+1)) THEN
              z4 = q1(npx+1)
            ELSE
              z4 = q1(npx)
            END IF
            IF (q1(npx-2) .GT. q1(npx-1)) THEN
              IF (q1(npx-1) .GT. z4) THEN
                y12 = z4
              ELSE
                y12 = q1(npx-1)
              END IF
            ELSE IF (q1(npx-2) .GT. z4) THEN
              y12 = z4
            ELSE
              y12 = q1(npx-2)
            END IF
            IF (xt .LT. y12) THEN
              xt = y12
            ELSE
              xt = xt
            END IF
            IF (q1(npx) .LT. q1(npx+1)) THEN
              z5 = q1(npx+1)
            ELSE
              z5 = q1(npx)
            END IF
            IF (q1(npx-2) .LT. q1(npx-1)) THEN
              IF (q1(npx-1) .LT. z5) THEN
                y13 = z5
              ELSE
                y13 = q1(npx-1)
              END IF
            ELSE IF (q1(npx-2) .LT. z5) THEN
              y13 = z5
            ELSE
              y13 = q1(npx-2)
            END IF
            IF (xt .GT. y13) THEN
              xt = y13
            ELSE
              xt = xt
            END IF
!        endif
            br(npx-1) = xt - q1(npx-1)
            bl(npx) = xt - q1(npx)
            br(npx) = s11*(q1(npx+1)-q1(npx)) - s14*dm(npx+1)
            CALL PERT_PPM(3, q1(npx-2:npx), bl(npx-2:npx), br(npx-2:npx)&
&                   , 1)
          END IF
        END IF
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux(i, j) = q1(i-1) + (1.-c(i, j))*(br(i-1)-c(i, j)*(bl(i-1&
&             )+br(i-1)))
          ELSE
            flux(i, j) = q1(i) + (1.+c(i, j))*(bl(i)+c(i, j)*(bl(i)+br(i&
&             )))
          END IF
        END DO
      END IF
    END DO
  END SUBROUTINE XPPM
!  Differentiation of yppm in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_core_
!mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn_core
!_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_mod.ge
!opk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynamics_mo
!d.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mod.c
!2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.map
!_scalar_fb fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scalar_pro
!file_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz fv
!_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_setup 
!fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.compute_
!pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_c nh_
!utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver nh_ut
!ils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw_core
!_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core_mod.
!fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.comput
!e_divergence_damping_fb sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corners t
!p_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circle_dis
!t sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: q flux c
!   with respect to varying inputs: q flux c
  SUBROUTINE YPPM_ADM(flux, flux_ad, q, q_ad, c, c_ad, jord, ifirst, &
&   ilast, isd, ied, js, je, jsd, jed, npx, npy, dya, nested, grid_type)
    IMPLICIT NONE
! Compute domain
    INTEGER, INTENT(IN) :: ifirst, ilast
    INTEGER, INTENT(IN) :: isd, ied, js, je, jsd, jed
    INTEGER, INTENT(IN) :: jord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst:ilast, jsd:jed)
    REAL :: q_ad(ifirst:ilast, jsd:jed)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
    REAL :: c_ad(isd:ied, js:je+1)
!  Flux
    REAL :: flux(ifirst:ilast, js:je+1)
    REAL :: flux_ad(ifirst:ilast, js:je+1)
    REAL, INTENT(IN) :: dya(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
! Local:
    REAL :: dm(ifirst:ilast, js-2:je+2)
    REAL :: dm_ad(ifirst:ilast, js-2:je+2)
    REAL :: al(ifirst:ilast, js-1:je+2)
    REAL :: al_ad(ifirst:ilast, js-1:je+2)
    REAL, DIMENSION(ifirst:ilast, js-1:je+1) :: bl, br, b0
    REAL, DIMENSION(ifirst:ilast, js-1:je+1) :: bl_ad, br_ad, b0_ad
    REAL :: dq(ifirst:ilast, js-3:je+2)
    REAL :: dq_ad(ifirst:ilast, js-3:je+2)
    REAL, DIMENSION(ifirst:ilast) :: fx0, fx1
    REAL, DIMENSION(ifirst:ilast) :: fx0_ad, fx1_ad
    LOGICAL, DIMENSION(ifirst:ilast, js-1:je+1) :: smt5, smt6
    REAL :: x0, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2, r1
    REAL :: xt_ad, qtmp_ad, pmp_1_ad, lac_1_ad, pmp_2_ad, lac_2_ad
    INTEGER :: i, j, js1, je3, je1
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min1_ad
    REAL :: min2
    REAL :: min2_ad
    REAL :: abs0
    REAL :: abs0_ad
    REAL :: abs1
    REAL :: min3
    REAL :: min3_ad
    REAL :: min4
    REAL :: min4_ad
    REAL :: min5
    REAL :: min5_ad
    REAL :: min6
    REAL :: min6_ad
    REAL :: min7
    REAL :: min7_ad
    REAL :: abs2
    REAL :: abs3
    REAL :: abs4
    REAL :: max1
    REAL :: max1_ad
    REAL :: min8
    REAL :: min8_ad
    REAL :: abs5
    REAL :: abs6
    REAL :: abs7
    INTEGER :: arg1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp0
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: x1_ad
    REAL :: y1_ad
    REAL :: x2_ad
    REAL :: y2_ad
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: temp_ad16
    REAL :: x3_ad
    REAL :: y3_ad
    REAL :: z1_ad
    REAL :: x4_ad
    REAL :: y4_ad
    REAL :: x5_ad
    REAL :: y5_ad
    REAL :: x6_ad
    REAL :: y6_ad
    REAL :: x7_ad
    REAL :: y7_ad
    REAL :: x8_ad
    REAL :: y14_ad
    REAL :: y8_ad
    REAL :: x9_ad
    REAL :: y15_ad
    REAL :: y9_ad
    REAL :: temp_ad17
    REAL :: temp_ad18
    REAL :: z2_ad
    REAL :: y10_ad
    REAL :: z3_ad
    REAL :: y11_ad
    REAL :: temp_ad19
    REAL :: temp_ad20
    REAL :: z4_ad
    REAL :: y12_ad
    REAL :: z5_ad
    REAL :: y13_ad
    REAL :: temp1
    REAL :: temp_ad21
    REAL :: temp_ad22
    REAL :: temp2
    REAL :: temp_ad23
    REAL :: temp_ad24
    INTEGER :: branch
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: x1
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
    REAL :: z5
    REAL :: z4
    REAL :: z3
    REAL :: z2
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
    IF (.NOT.nested .AND. grid_type .LT. 3) THEN
      IF (3 .LT. js - 1) THEN
        js1 = js - 1
      ELSE
        js1 = 3
      END IF
      IF (npy - 2 .GT. je + 2) THEN
        je3 = je + 2
      ELSE
        je3 = npy - 2
      END IF
      IF (npy - 3 .GT. je + 1) THEN
        CALL PUSHCONTROL1B(1)
        je1 = je + 1
      ELSE
        CALL PUSHCONTROL1B(1)
        je1 = npy - 3
      END IF
    ELSE
      CALL PUSHCONTROL1B(0)
! Nested grid OR Doubly periodic domain:
      js1 = js - 1
      je3 = je + 2
      je1 = je + 1
    END IF
    IF (jord .LT. 8 .OR. jord .EQ. 333) THEN
      DO j=js1,je3
        DO i=ifirst,ilast
          al(i, j) = p1*(q(i, j-1)+q(i, j)) + p2*(q(i, j-2)+q(i, j+1))
        END DO
      END DO
      IF (jord .EQ. 7) THEN
        DO j=js1,je3
          DO i=ifirst,ilast
            IF (al(i, j) .LT. 0.) THEN
              al(i, j) = 0.5*(q(i, j)+q(i, j+1))
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            al(i, 0) = c1*q(i, -2) + c2*q(i, -1) + c3*q(i, 0)
            al(i, 1) = 0.5*(((2.*dya(i, 0)+dya(i, -1))*q(i, 0)-dya(i, 0)&
&             *q(i, -1))/(dya(i, -1)+dya(i, 0))+((2.*dya(i, 1)+dya(i, 2)&
&             )*q(i, 1)-dya(i, 1)*q(i, 2))/(dya(i, 1)+dya(i, 2)))
            al(i, 2) = c3*q(i, 1) + c2*q(i, 2) + c1*q(i, 3)
          END DO
          IF (jord .EQ. 7) THEN
            DO i=ifirst,ilast
              IF (0. .LT. al(i, 0)) THEN
                CALL PUSHCONTROL1B(0)
                al(i, 0) = al(i, 0)
              ELSE
                al(i, 0) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
              IF (0. .LT. al(i, 1)) THEN
                CALL PUSHCONTROL1B(0)
                al(i, 1) = al(i, 1)
              ELSE
                al(i, 1) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
              IF (0. .LT. al(i, 2)) THEN
                CALL PUSHCONTROL1B(0)
                al(i, 2) = al(i, 2)
              ELSE
                al(i, 2) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
            END DO
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            al(i, npy-1) = c1*q(i, npy-3) + c2*q(i, npy-2) + c3*q(i, npy&
&             -1)
            al(i, npy) = 0.5*(((2.*dya(i, npy-1)+dya(i, npy-2))*q(i, npy&
&             -1)-dya(i, npy-1)*q(i, npy-2))/(dya(i, npy-2)+dya(i, npy-1&
&             ))+((2.*dya(i, npy)+dya(i, npy+1))*q(i, npy)-dya(i, npy)*q&
&             (i, npy+1))/(dya(i, npy)+dya(i, npy+1)))
            al(i, npy+1) = c3*q(i, npy) + c2*q(i, npy+1) + c1*q(i, npy+2&
&             )
          END DO
          IF (jord .EQ. 7) THEN
            DO i=ifirst,ilast
              IF (0. .LT. al(i, npy-1)) THEN
                CALL PUSHCONTROL1B(0)
                al(i, npy-1) = al(i, npy-1)
              ELSE
                al(i, npy-1) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
              IF (0. .LT. al(i, npy)) THEN
                CALL PUSHCONTROL1B(0)
                al(i, npy) = al(i, npy)
              ELSE
                al(i, npy) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
              IF (0. .LT. al(i, npy+1)) THEN
                CALL PUSHCONTROL1B(0)
                al(i, npy+1) = al(i, npy+1)
              ELSE
                al(i, npy+1) = 0.
                CALL PUSHCONTROL1B(1)
              END IF
            END DO
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
      ELSE
        CALL PUSHCONTROL2B(3)
      END IF
      IF (jord .EQ. 1) THEN
        DO j=js,je+1
          DO i=ifirst,ilast
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ilast,ifirst,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q_ad(i, j) = q_ad(i, j) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              q_ad(i, j-1) = q_ad(i, j-1) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            END IF
          END DO
        END DO
        al_ad = 0.0
      ELSE IF (jord .EQ. 2) THEN
! Perfectly linear scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
        END DO
        al_ad = 0.0
        DO j=je+1,js,-1
          DO i=ilast,ifirst,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt = c(i, j)
              qtmp = q(i, j)
              temp0 = al(i, j) + al(i, j+1) - 2*qtmp
              temp_ad5 = (xt+1.)*flux_ad(i, j)
              temp_ad6 = xt*temp_ad5
              qtmp_ad = flux_ad(i, j) - temp_ad5 - 2*temp_ad6
              xt_ad = temp0*temp_ad5 + (al(i, j)-qtmp+xt*temp0)*flux_ad(&
&               i, j)
              al_ad(i, j) = al_ad(i, j) + temp_ad6 + temp_ad5
              al_ad(i, j+1) = al_ad(i, j+1) + temp_ad6
              flux_ad(i, j) = 0.0
              q_ad(i, j) = q_ad(i, j) + qtmp_ad
            ELSE
              xt = c(i, j)
              qtmp = q(i, j-1)
              temp = al(i, j-1) + al(i, j) - 2*qtmp
              temp_ad3 = (1.-xt)*flux_ad(i, j)
              temp_ad4 = -(xt*temp_ad3)
              qtmp_ad = flux_ad(i, j) - temp_ad3 - 2*temp_ad4
              xt_ad = -(temp*temp_ad3) - (al(i, j)-qtmp-xt*temp)*flux_ad&
&               (i, j)
              al_ad(i, j) = al_ad(i, j) + temp_ad4 + temp_ad3
              al_ad(i, j-1) = al_ad(i, j-1) + temp_ad4
              flux_ad(i, j) = 0.0
              q_ad(i, j-1) = q_ad(i, j-1) + qtmp_ad
            END IF
            c_ad(i, j) = c_ad(i, j) + xt_ad
          END DO
        END DO
      ELSE IF (jord .EQ. 333) THEN
! Perfectly linear scheme, more diffusive than ord=2 (HoldawayKent-2015-TellusA)
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ilast,ifirst,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt = c(i, j)
              temp_ad10 = flux_ad(i, j)/6.0
              temp_ad11 = -(0.5*xt*flux_ad(i, j))
              temp_ad12 = xt**2*flux_ad(i, j)/6.0
              q_ad(i, j-1) = q_ad(i, j-1) + temp_ad12 - temp_ad11 + 2.0*&
&               temp_ad10
              q_ad(i, j) = q_ad(i, j) + temp_ad11 - 2.0*temp_ad12 + 5.0*&
&               temp_ad10
              q_ad(i, j+1) = q_ad(i, j+1) + temp_ad12 - temp_ad10
              xt_ad = ((q(i, j+1)-2.0*q(i, j)+q(i, j-1))*2*xt/6.0-0.5*(q&
&               (i, j)-q(i, j-1)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              xt = c(i, j)
              temp_ad7 = flux_ad(i, j)/6.0
              temp_ad8 = -(0.5*xt*flux_ad(i, j))
              temp_ad9 = xt**2*flux_ad(i, j)/6.0
              q_ad(i, j) = q_ad(i, j) + temp_ad9 + temp_ad8 + 2.0*&
&               temp_ad7
              q_ad(i, j-1) = q_ad(i, j-1) + 5.0*temp_ad7 - temp_ad8 - &
&               2.0*temp_ad9
              q_ad(i, j-2) = q_ad(i, j-2) + temp_ad9 - temp_ad7
              xt_ad = ((q(i, j)-2.0*q(i, j-1)+q(i, j-2))*2*xt/6.0-0.5*(q&
&               (i, j)-q(i, j-1)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0
            END IF
            c_ad(i, j) = c_ad(i, j) + xt_ad
          END DO
        END DO
        al_ad = 0.0
      ELSE IF (jord .EQ. 3) THEN
        DO j=js-1,je+1
          DO i=ifirst,ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j+1) - q(i, j)
            b0(i, j) = bl(i, j) + br(i, j)
            IF (b0(i, j) .GE. 0.) THEN
              x0 = b0(i, j)
            ELSE
              x0 = -b0(i, j)
            END IF
            IF (bl(i, j) - br(i, j) .GE. 0.) THEN
              xt = bl(i, j) - br(i, j)
            ELSE
              xt = -(bl(i, j)-br(i, j))
            END IF
            smt5(i, j) = x0 .LT. xt
            smt6(i, j) = 3.*x0 .LT. xt
          END DO
        END DO
        DO j=js,je+1
          DO i=ifirst,ilast
            CALL PUSHREALARRAY_ADM(fx1(i))
            fx1(i) = 0.
          END DO
          DO i=ifirst,ilast
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              IF (smt6(i, j-1) .OR. smt5(i, j)) THEN
                CALL PUSHREALARRAY_ADM(fx1(i))
                fx1(i) = br(i, j-1) - xt*b0(i, j-1)
                CALL PUSHCONTROL3B(5)
              ELSE IF (smt5(i, j-1)) THEN
                IF (bl(i, j-1) .GE. 0.) THEN
                  x1 = bl(i, j-1)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  x1 = -bl(i, j-1)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (br(i, j-1) .GE. 0.) THEN
                  y1 = br(i, j-1)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  y1 = -br(i, j-1)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (x1 .GT. y1) THEN
                  CALL PUSHREALARRAY_ADM(min1)
                  min1 = y1
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHREALARRAY_ADM(min1)
                  min1 = x1
                  CALL PUSHCONTROL1B(1)
                END IF
! both up-downwind sides are noisy; 2nd order, piece-wise linear
                CALL PUSHREALARRAY_ADM(fx1(i))
                fx1(i) = SIGN(min1, br(i, j-1))
                CALL PUSHCONTROL3B(4)
              ELSE
                CALL PUSHCONTROL3B(3)
              END IF
            ELSE IF (smt6(i, j) .OR. smt5(i, j-1)) THEN
              CALL PUSHREALARRAY_ADM(fx1(i))
              fx1(i) = bl(i, j) + xt*b0(i, j)
              CALL PUSHCONTROL3B(2)
            ELSE IF (smt5(i, j)) THEN
              IF (bl(i, j) .GE. 0.) THEN
                x2 = bl(i, j)
                CALL PUSHCONTROL1B(0)
              ELSE
                x2 = -bl(i, j)
                CALL PUSHCONTROL1B(1)
              END IF
              IF (br(i, j) .GE. 0.) THEN
                y2 = br(i, j)
                CALL PUSHCONTROL1B(0)
              ELSE
                y2 = -br(i, j)
                CALL PUSHCONTROL1B(1)
              END IF
              IF (x2 .GT. y2) THEN
                CALL PUSHREALARRAY_ADM(min2)
                min2 = y2
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHREALARRAY_ADM(min2)
                min2 = x2
                CALL PUSHCONTROL1B(1)
              END IF
              CALL PUSHREALARRAY_ADM(fx1(i))
              fx1(i) = SIGN(min2, bl(i, j))
              CALL PUSHCONTROL3B(1)
            ELSE
              CALL PUSHCONTROL3B(0)
            END IF
            IF (xt .GE. 0.) THEN
              CALL PUSHREALARRAY_ADM(abs0)
              abs0 = xt
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(abs0)
              abs0 = -xt
              CALL PUSHCONTROL1B(1)
            END IF
          END DO
        END DO
        bl_ad = 0.0
        br_ad = 0.0
        b0_ad = 0.0
        fx0_ad = 0.0
        fx1_ad = 0.0
        DO j=je+1,js,-1
          DO i=ilast,ifirst,-1
            fx0_ad(i) = fx0_ad(i) + flux_ad(i, j)
            abs0_ad = -(fx1(i)*flux_ad(i, j))
            fx1_ad(i) = fx1_ad(i) + (1.-abs0)*flux_ad(i, j)
            flux_ad(i, j) = 0.0
            xt = c(i, j)
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(abs0)
              xt_ad = abs0_ad
            ELSE
              CALL POPREALARRAY_ADM(abs0)
              xt_ad = -abs0_ad
            END IF
            CALL POPCONTROL3B(branch)
            IF (branch .LT. 3) THEN
              IF (branch .NE. 0) THEN
                IF (branch .EQ. 1) THEN
                  CALL POPREALARRAY_ADM(fx1(i))
                  min2_ad = SIGN(1.d0, min2*bl(i, j))*fx1_ad(i)
                  fx1_ad(i) = 0.0
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(min2)
                    y2_ad = min2_ad
                    x2_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(min2)
                    x2_ad = min2_ad
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
                    bl_ad(i, j) = bl_ad(i, j) + x2_ad
                  ELSE
                    bl_ad(i, j) = bl_ad(i, j) - x2_ad
                  END IF
                ELSE
                  CALL POPREALARRAY_ADM(fx1(i))
                  bl_ad(i, j) = bl_ad(i, j) + fx1_ad(i)
                  xt_ad = xt_ad + b0(i, j)*fx1_ad(i)
                  b0_ad(i, j) = b0_ad(i, j) + xt*fx1_ad(i)
                  fx1_ad(i) = 0.0
                END IF
              END IF
              q_ad(i, j) = q_ad(i, j) + fx0_ad(i)
              fx0_ad(i) = 0.0
            ELSE
              IF (branch .NE. 3) THEN
                IF (branch .EQ. 4) THEN
                  CALL POPREALARRAY_ADM(fx1(i))
                  min1_ad = SIGN(1.d0, min1*br(i, j-1))*fx1_ad(i)
                  fx1_ad(i) = 0.0
                  CALL POPCONTROL1B(branch)
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY_ADM(min1)
                    y1_ad = min1_ad
                    x1_ad = 0.0
                  ELSE
                    CALL POPREALARRAY_ADM(min1)
                    x1_ad = min1_ad
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
                    bl_ad(i, j-1) = bl_ad(i, j-1) + x1_ad
                  ELSE
                    bl_ad(i, j-1) = bl_ad(i, j-1) - x1_ad
                  END IF
                ELSE
                  CALL POPREALARRAY_ADM(fx1(i))
                  br_ad(i, j-1) = br_ad(i, j-1) + fx1_ad(i)
                  xt_ad = xt_ad - b0(i, j-1)*fx1_ad(i)
                  b0_ad(i, j-1) = b0_ad(i, j-1) - xt*fx1_ad(i)
                  fx1_ad(i) = 0.0
                END IF
              END IF
              q_ad(i, j-1) = q_ad(i, j-1) + fx0_ad(i)
              fx0_ad(i) = 0.0
            END IF
            c_ad(i, j) = c_ad(i, j) + xt_ad
          END DO
          DO i=ilast,ifirst,-1
            CALL POPREALARRAY_ADM(fx1(i))
            fx1_ad(i) = 0.0
          END DO
        END DO
        al_ad = 0.0
        DO j=je+1,js-1,-1
          DO i=ilast,ifirst,-1
            bl_ad(i, j) = bl_ad(i, j) + b0_ad(i, j)
            br_ad(i, j) = br_ad(i, j) + b0_ad(i, j)
            b0_ad(i, j) = 0.0
            al_ad(i, j+1) = al_ad(i, j+1) + br_ad(i, j)
            q_ad(i, j) = q_ad(i, j) - bl_ad(i, j) - br_ad(i, j)
            br_ad(i, j) = 0.0
            al_ad(i, j) = al_ad(i, j) + bl_ad(i, j)
            bl_ad(i, j) = 0.0
          END DO
        END DO
      ELSE IF (jord .EQ. 4) THEN
        DO j=js-1,je+1
          DO i=ifirst,ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j+1) - q(i, j)
            b0(i, j) = bl(i, j) + br(i, j)
            IF (b0(i, j) .GE. 0.) THEN
              x0 = b0(i, j)
            ELSE
              x0 = -b0(i, j)
            END IF
            IF (bl(i, j) - br(i, j) .GE. 0.) THEN
              xt = bl(i, j) - br(i, j)
            ELSE
              xt = -(bl(i, j)-br(i, j))
            END IF
            smt5(i, j) = x0 .LT. xt
            smt6(i, j) = 3.*x0 .LT. xt
          END DO
        END DO
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
            IF (c(i, j) .GT. 0.) THEN
              IF (smt6(i, j-1) .OR. smt5(i, j)) THEN
                CALL PUSHCONTROL2B(0)
              ELSE
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (smt6(i, j) .OR. smt5(i, j-1)) THEN
              CALL PUSHCONTROL2B(2)
            ELSE
              CALL PUSHCONTROL2B(3)
            END IF
          END DO
        END DO
        bl_ad = 0.0
        br_ad = 0.0
        b0_ad = 0.0
        fx0_ad = 0.0
        fx1_ad = 0.0
        DO j=je+1,js,-1
          DO i=ilast,ifirst,-1
            fx0_ad(i) = fx0_ad(i) + flux_ad(i, j)
            fx1_ad(i) = fx1_ad(i) + flux_ad(i, j)
            flux_ad(i, j) = 0.0
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                temp_ad13 = (1.-c(i, j))*fx1_ad(i)
                c_ad(i, j) = c_ad(i, j) - b0(i, j-1)*temp_ad13 - (br(i, &
&                 j-1)-c(i, j)*b0(i, j-1))*fx1_ad(i)
                br_ad(i, j-1) = br_ad(i, j-1) + temp_ad13
                b0_ad(i, j-1) = b0_ad(i, j-1) - c(i, j)*temp_ad13
                fx1_ad(i) = 0.0
              END IF
              q_ad(i, j-1) = q_ad(i, j-1) + fx0_ad(i)
              fx0_ad(i) = 0.0
            ELSE
              IF (branch .EQ. 2) THEN
                temp_ad14 = (c(i, j)+1.)*fx1_ad(i)
                c_ad(i, j) = c_ad(i, j) + b0(i, j)*temp_ad14 + (bl(i, j)&
&                 +c(i, j)*b0(i, j))*fx1_ad(i)
                bl_ad(i, j) = bl_ad(i, j) + temp_ad14
                b0_ad(i, j) = b0_ad(i, j) + c(i, j)*temp_ad14
                fx1_ad(i) = 0.0
              END IF
              q_ad(i, j) = q_ad(i, j) + fx0_ad(i)
              fx0_ad(i) = 0.0
            END IF
          END DO
          DO i=ilast,ifirst,-1
            fx1_ad(i) = 0.0
          END DO
        END DO
        al_ad = 0.0
        DO j=je+1,js-1,-1
          DO i=ilast,ifirst,-1
            bl_ad(i, j) = bl_ad(i, j) + b0_ad(i, j)
            br_ad(i, j) = br_ad(i, j) + b0_ad(i, j)
            b0_ad(i, j) = 0.0
            al_ad(i, j+1) = al_ad(i, j+1) + br_ad(i, j)
            q_ad(i, j) = q_ad(i, j) - bl_ad(i, j) - br_ad(i, j)
            br_ad(i, j) = 0.0
            al_ad(i, j) = al_ad(i, j) + bl_ad(i, j)
            bl_ad(i, j) = 0.0
          END DO
        END DO
      ELSE
! jord=5,6,7
        IF (jord .EQ. 5) THEN
          DO j=js-1,je+1
            DO i=ifirst,ilast
              bl(i, j) = al(i, j) - q(i, j)
              br(i, j) = al(i, j+1) - q(i, j)
              b0(i, j) = bl(i, j) + br(i, j)
              smt5(i, j) = bl(i, j)*br(i, j) .LT. 0.
            END DO
          END DO
          CALL PUSHCONTROL1B(1)
        ELSE
          DO j=js-1,je+1
            DO i=ifirst,ilast
              bl(i, j) = al(i, j) - q(i, j)
              br(i, j) = al(i, j+1) - q(i, j)
              b0(i, j) = bl(i, j) + br(i, j)
              IF (3.*b0(i, j) .GE. 0.) THEN
                abs1 = 3.*b0(i, j)
              ELSE
                abs1 = -(3.*b0(i, j))
              END IF
              IF (bl(i, j) - br(i, j) .GE. 0.) THEN
                abs4 = bl(i, j) - br(i, j)
              ELSE
                abs4 = -(bl(i, j)-br(i, j))
              END IF
              smt5(i, j) = abs1 .LT. abs4
            END DO
          END DO
          CALL PUSHCONTROL1B(0)
        END IF
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
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
        fx1_ad = 0.0
        DO j=je+1,js,-1
          DO i=ilast,ifirst,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) fx1_ad(i) = fx1_ad(i) + flux_ad(i, j)
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q_ad(i, j-1) = q_ad(i, j-1) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
              temp_ad15 = (1.-c(i, j))*fx1_ad(i)
              c_ad(i, j) = c_ad(i, j) - b0(i, j-1)*temp_ad15 - (br(i, j-&
&               1)-c(i, j)*b0(i, j-1))*fx1_ad(i)
              br_ad(i, j-1) = br_ad(i, j-1) + temp_ad15
              b0_ad(i, j-1) = b0_ad(i, j-1) - c(i, j)*temp_ad15
              fx1_ad(i) = 0.0
            ELSE
              q_ad(i, j) = q_ad(i, j) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
              temp_ad16 = (c(i, j)+1.)*fx1_ad(i)
              c_ad(i, j) = c_ad(i, j) + b0(i, j)*temp_ad16 + (bl(i, j)+c&
&               (i, j)*b0(i, j))*fx1_ad(i)
              bl_ad(i, j) = bl_ad(i, j) + temp_ad16
              b0_ad(i, j) = b0_ad(i, j) + c(i, j)*temp_ad16
              fx1_ad(i) = 0.0
            END IF
          END DO
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          al_ad = 0.0
          DO j=je+1,js-1,-1
            DO i=ilast,ifirst,-1
              bl_ad(i, j) = bl_ad(i, j) + b0_ad(i, j)
              br_ad(i, j) = br_ad(i, j) + b0_ad(i, j)
              b0_ad(i, j) = 0.0
              al_ad(i, j+1) = al_ad(i, j+1) + br_ad(i, j)
              q_ad(i, j) = q_ad(i, j) - bl_ad(i, j) - br_ad(i, j)
              br_ad(i, j) = 0.0
              al_ad(i, j) = al_ad(i, j) + bl_ad(i, j)
              bl_ad(i, j) = 0.0
            END DO
          END DO
        ELSE
          al_ad = 0.0
          DO j=je+1,js-1,-1
            DO i=ilast,ifirst,-1
              bl_ad(i, j) = bl_ad(i, j) + b0_ad(i, j)
              br_ad(i, j) = br_ad(i, j) + b0_ad(i, j)
              b0_ad(i, j) = 0.0
              al_ad(i, j+1) = al_ad(i, j+1) + br_ad(i, j)
              q_ad(i, j) = q_ad(i, j) - bl_ad(i, j) - br_ad(i, j)
              br_ad(i, j) = 0.0
              al_ad(i, j) = al_ad(i, j) + bl_ad(i, j)
              bl_ad(i, j) = 0.0
            END DO
          END DO
        END IF
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          DO i=ilast,ifirst,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) al_ad(i, npy+1) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) al_ad(i, npy) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) al_ad(i, npy-1) = 0.0
          END DO
        END IF
        DO i=ilast,ifirst,-1
          q_ad(i, npy) = q_ad(i, npy) + c3*al_ad(i, npy+1)
          q_ad(i, npy+1) = q_ad(i, npy+1) + c2*al_ad(i, npy+1)
          q_ad(i, npy+2) = q_ad(i, npy+2) + c1*al_ad(i, npy+1)
          al_ad(i, npy+1) = 0.0
          temp_ad1 = 0.5*al_ad(i, npy)/(dya(i, npy-2)+dya(i, npy-1))
          temp_ad2 = 0.5*al_ad(i, npy)/(dya(i, npy)+dya(i, npy+1))
          q_ad(i, npy-1) = q_ad(i, npy-1) + (dya(i, npy-1)*2.+dya(i, npy&
&           -2))*temp_ad1
          q_ad(i, npy-2) = q_ad(i, npy-2) - dya(i, npy-1)*temp_ad1
          q_ad(i, npy) = q_ad(i, npy) + (dya(i, npy)*2.+dya(i, npy+1))*&
&           temp_ad2
          q_ad(i, npy+1) = q_ad(i, npy+1) - dya(i, npy)*temp_ad2
          al_ad(i, npy) = 0.0
          q_ad(i, npy-3) = q_ad(i, npy-3) + c1*al_ad(i, npy-1)
          q_ad(i, npy-2) = q_ad(i, npy-2) + c2*al_ad(i, npy-1)
          q_ad(i, npy-1) = q_ad(i, npy-1) + c3*al_ad(i, npy-1)
          al_ad(i, npy-1) = 0.0
        END DO
      ELSE IF (branch .NE. 2) THEN
        GOTO 100
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        DO i=ilast,ifirst,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) al_ad(i, 2) = 0.0
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) al_ad(i, 1) = 0.0
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) al_ad(i, 0) = 0.0
        END DO
      ELSE IF (branch .NE. 1) THEN
        GOTO 100
      END IF
      DO i=ilast,ifirst,-1
        q_ad(i, 1) = q_ad(i, 1) + c3*al_ad(i, 2)
        q_ad(i, 2) = q_ad(i, 2) + c2*al_ad(i, 2)
        q_ad(i, 3) = q_ad(i, 3) + c1*al_ad(i, 2)
        al_ad(i, 2) = 0.0
        temp_ad = 0.5*al_ad(i, 1)/(dya(i, -1)+dya(i, 0))
        temp_ad0 = 0.5*al_ad(i, 1)/(dya(i, 1)+dya(i, 2))
        q_ad(i, 0) = q_ad(i, 0) + (dya(i, 0)*2.+dya(i, -1))*temp_ad
        q_ad(i, -1) = q_ad(i, -1) - dya(i, 0)*temp_ad
        q_ad(i, 1) = q_ad(i, 1) + (dya(i, 1)*2.+dya(i, 2))*temp_ad0
        q_ad(i, 2) = q_ad(i, 2) - dya(i, 1)*temp_ad0
        al_ad(i, 1) = 0.0
        q_ad(i, -2) = q_ad(i, -2) + c1*al_ad(i, 0)
        q_ad(i, -1) = q_ad(i, -1) + c2*al_ad(i, 0)
        q_ad(i, 0) = q_ad(i, 0) + c3*al_ad(i, 0)
        al_ad(i, 0) = 0.0
      END DO
 100  CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je3,js1,-1
          DO i=ilast,ifirst,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) THEN
              q_ad(i, j) = q_ad(i, j) + 0.5*al_ad(i, j)
              q_ad(i, j+1) = q_ad(i, j+1) + 0.5*al_ad(i, j)
              al_ad(i, j) = 0.0
            END IF
          END DO
        END DO
      END IF
      DO j=je3,js1,-1
        DO i=ilast,ifirst,-1
          q_ad(i, j-1) = q_ad(i, j-1) + p1*al_ad(i, j)
          q_ad(i, j) = q_ad(i, j) + p1*al_ad(i, j)
          q_ad(i, j-2) = q_ad(i, j-2) + p2*al_ad(i, j)
          q_ad(i, j+1) = q_ad(i, j+1) + p2*al_ad(i, j)
          al_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
! Monotonic constraints:
! ord = 8: PPM with Lin's PPM fast monotone constraint
! ord > 8: PPM with Lin's modification of Huynh 2nd constraint
      DO j=js-2,je+2
        DO i=ifirst,ilast
          xt = 0.25*(q(i, j+1)-q(i, j-1))
          IF (xt .GE. 0.) THEN
            x3 = xt
            CALL PUSHCONTROL1B(0)
          ELSE
            x3 = -xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q(i, j-1) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i, j+1)) THEN
              max1 = q(i, j+1)
              CALL PUSHCONTROL2B(0)
            ELSE
              max1 = q(i, j)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (q(i, j-1) .LT. q(i, j+1)) THEN
            max1 = q(i, j+1)
            CALL PUSHCONTROL2B(2)
          ELSE
            max1 = q(i, j-1)
            CALL PUSHCONTROL2B(3)
          END IF
          y3 = max1 - q(i, j)
          IF (q(i, j-1) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i, j+1)) THEN
              min8 = q(i, j+1)
              CALL PUSHCONTROL2B(0)
            ELSE
              min8 = q(i, j)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (q(i, j-1) .GT. q(i, j+1)) THEN
            min8 = q(i, j+1)
            CALL PUSHCONTROL2B(2)
          ELSE
            min8 = q(i, j-1)
            CALL PUSHCONTROL2B(3)
          END IF
          z1 = q(i, j) - min8
          IF (x3 .GT. y3) THEN
            IF (y3 .GT. z1) THEN
              CALL PUSHREALARRAY_ADM(min3)
              min3 = z1
              CALL PUSHCONTROL2B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(min3)
              min3 = y3
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (x3 .GT. z1) THEN
            CALL PUSHREALARRAY_ADM(min3)
            min3 = z1
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHREALARRAY_ADM(min3)
            min3 = x3
            CALL PUSHCONTROL2B(3)
          END IF
          dm(i, j) = SIGN(min3, xt)
        END DO
      END DO
      DO j=js1,je1+1
        DO i=ifirst,ilast
          al(i, j) = 0.5*(q(i, j-1)+q(i, j)) + r3*(dm(i, j-1)-dm(i, j))
        END DO
      END DO
      IF (jord .EQ. 8) THEN
        DO j=js1,je1
          DO i=ifirst,ilast
            xt = 2.*dm(i, j)
            IF (xt .GE. 0.) THEN
              x4 = xt
              CALL PUSHCONTROL1B(0)
            ELSE
              x4 = -xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (al(i, j) - q(i, j) .GE. 0.) THEN
              y4 = al(i, j) - q(i, j)
              CALL PUSHCONTROL1B(0)
            ELSE
              y4 = -(al(i, j)-q(i, j))
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x4 .GT. y4) THEN
              CALL PUSHREALARRAY_ADM(min4)
              min4 = y4
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(min4)
              min4 = x4
              CALL PUSHCONTROL1B(1)
            END IF
            bl(i, j) = -SIGN(min4, xt)
            IF (xt .GE. 0.) THEN
              x5 = xt
              CALL PUSHCONTROL1B(0)
            ELSE
              x5 = -xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (al(i, j+1) - q(i, j) .GE. 0.) THEN
              y5 = al(i, j+1) - q(i, j)
              CALL PUSHCONTROL1B(0)
            ELSE
              y5 = -(al(i, j+1)-q(i, j))
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x5 .GT. y5) THEN
              CALL PUSHREALARRAY_ADM(min5)
              min5 = y5
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(min5)
              min5 = x5
              CALL PUSHCONTROL1B(1)
            END IF
            br(i, j) = SIGN(min5, xt)
          END DO
        END DO
        CALL PUSHCONTROL2B(0)
      ELSE IF (jord .EQ. 11) THEN
        DO j=js1,je1
          DO i=ifirst,ilast
            xt = ppm_fac*dm(i, j)
            IF (xt .GE. 0.) THEN
              x6 = xt
              CALL PUSHCONTROL1B(0)
            ELSE
              x6 = -xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (al(i, j) - q(i, j) .GE. 0.) THEN
              y6 = al(i, j) - q(i, j)
              CALL PUSHCONTROL1B(0)
            ELSE
              y6 = -(al(i, j)-q(i, j))
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x6 .GT. y6) THEN
              CALL PUSHREALARRAY_ADM(min6)
              min6 = y6
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(min6)
              min6 = x6
              CALL PUSHCONTROL1B(1)
            END IF
            bl(i, j) = -SIGN(min6, xt)
            IF (xt .GE. 0.) THEN
              x7 = xt
              CALL PUSHCONTROL1B(0)
            ELSE
              x7 = -xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (al(i, j+1) - q(i, j) .GE. 0.) THEN
              y7 = al(i, j+1) - q(i, j)
              CALL PUSHCONTROL1B(0)
            ELSE
              y7 = -(al(i, j+1)-q(i, j))
              CALL PUSHCONTROL1B(1)
            END IF
            IF (x7 .GT. y7) THEN
              CALL PUSHREALARRAY_ADM(min7)
              min7 = y7
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREALARRAY_ADM(min7)
              min7 = x7
              CALL PUSHCONTROL1B(1)
            END IF
            br(i, j) = SIGN(min7, xt)
          END DO
        END DO
        CALL PUSHCONTROL2B(1)
      ELSE
        DO j=js1-2,je1+1
          DO i=ifirst,ilast
            dq(i, j) = 2.*(q(i, j+1)-q(i, j))
          END DO
        END DO
        DO j=js1,je1
          DO i=ifirst,ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j+1) - q(i, j)
            IF (dm(i, j-1) .GE. 0.) THEN
              abs2 = dm(i, j-1)
            ELSE
              abs2 = -dm(i, j-1)
            END IF
            IF (dm(i, j) .GE. 0.) THEN
              abs5 = dm(i, j)
            ELSE
              abs5 = -dm(i, j)
            END IF
            IF (dm(i, j+1) .GE. 0.) THEN
              abs7 = dm(i, j+1)
            ELSE
              abs7 = -dm(i, j+1)
            END IF
            IF (abs2 + abs5 + abs7 .LT. near_zero) THEN
              bl(i, j) = 0.
              br(i, j) = 0.
              CALL PUSHCONTROL2B(3)
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
                pmp_2 = dq(i, j-1)
                lac_2 = pmp_2 - 0.75*dq(i, j-2)
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
                    y14 = lac_2
                    CALL PUSHCONTROL2B(0)
                  ELSE
                    y14 = pmp_2
                    CALL PUSHCONTROL2B(1)
                  END IF
                ELSE IF (0. .GT. lac_2) THEN
                  y14 = lac_2
                  CALL PUSHCONTROL2B(2)
                ELSE
                  y14 = 0.
                  CALL PUSHCONTROL2B(3)
                END IF
                IF (br(i, j) .LT. y14) THEN
                  y8 = y14
                  CALL PUSHCONTROL1B(0)
                ELSE
                  y8 = br(i, j)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (x8 .GT. y8) THEN
                  br(i, j) = y8
                  CALL PUSHCONTROL1B(0)
                ELSE
                  br(i, j) = x8
                  CALL PUSHCONTROL1B(1)
                END IF
                pmp_1 = -dq(i, j)
                lac_1 = pmp_1 + 0.75*dq(i, j+1)
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
                    y15 = lac_1
                    CALL PUSHCONTROL2B(0)
                  ELSE
                    y15 = pmp_1
                    CALL PUSHCONTROL2B(1)
                  END IF
                ELSE IF (0. .GT. lac_1) THEN
                  y15 = lac_1
                  CALL PUSHCONTROL2B(2)
                ELSE
                  y15 = 0.
                  CALL PUSHCONTROL2B(3)
                END IF
                IF (bl(i, j) .LT. y15) THEN
                  y9 = y15
                  CALL PUSHCONTROL1B(0)
                ELSE
                  y9 = bl(i, j)
                  CALL PUSHCONTROL1B(1)
                END IF
                IF (x9 .GT. y9) THEN
                  bl(i, j) = y9
                  CALL PUSHCONTROL2B(1)
                ELSE
                  bl(i, j) = x9
                  CALL PUSHCONTROL2B(2)
                END IF
              ELSE
                CALL PUSHCONTROL2B(0)
              END IF
            END IF
          END DO
        END DO
        CALL PUSHCONTROL2B(2)
      END IF
      IF (jord .EQ. 9 .OR. jord .EQ. 13) THEN
! Positive definite constraint:
        DO j=js1,je1
          arg1 = ilast - ifirst + 1
          CALL PUSHREALARRAY_ADM(br(:, j), ilast - ifirst + 1)
          CALL PUSHREALARRAY_ADM(bl(:, j), ilast - ifirst + 1)
          CALL PERT_PPM(arg1, q(ifirst:ilast, j), bl(ifirst:ilast, j), &
&                 br(ifirst:ilast, j), 0)
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            bl(i, 0) = s14*dm(i, -1) + s11*(q(i, -1)-q(i, 0))
            xt = 0.5*(((2.*dya(i, 0)+dya(i, -1))*q(i, 0)-dya(i, 0)*q(i, &
&             -1))/(dya(i, -1)+dya(i, 0))+((2.*dya(i, 1)+dya(i, 2))*q(i&
&             , 1)-dya(i, 1)*q(i, 2))/(dya(i, 1)+dya(i, 2)))
            IF (q(i, 1) .GT. q(i, 2)) THEN
              z2 = q(i, 2)
              CALL PUSHCONTROL1B(0)
            ELSE
              z2 = q(i, 1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, -1) .GT. q(i, 0)) THEN
              IF (q(i, 0) .GT. z2) THEN
                y10 = z2
                CALL PUSHCONTROL2B(0)
              ELSE
                y10 = q(i, 0)
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (q(i, -1) .GT. z2) THEN
              y10 = z2
              CALL PUSHCONTROL2B(2)
            ELSE
              y10 = q(i, -1)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (xt .LT. y10) THEN
              xt = y10
              CALL PUSHCONTROL1B(0)
            ELSE
              xt = xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, 1) .LT. q(i, 2)) THEN
              z3 = q(i, 2)
              CALL PUSHCONTROL1B(0)
            ELSE
              z3 = q(i, 1)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, -1) .LT. q(i, 0)) THEN
              IF (q(i, 0) .LT. z3) THEN
                y11 = z3
                CALL PUSHCONTROL2B(0)
              ELSE
                y11 = q(i, 0)
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (q(i, -1) .LT. z3) THEN
              y11 = z3
              CALL PUSHCONTROL2B(2)
            ELSE
              y11 = q(i, -1)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (xt .GT. y11) THEN
              xt = y11
              CALL PUSHCONTROL1B(0)
            ELSE
              xt = xt
              CALL PUSHCONTROL1B(1)
            END IF
!        endif
            br(i, 0) = xt - q(i, 0)
            bl(i, 1) = xt - q(i, 1)
            xt = s15*q(i, 1) + s11*q(i, 2) - s14*dm(i, 2)
            br(i, 1) = xt - q(i, 1)
            bl(i, 2) = xt - q(i, 2)
            br(i, 2) = al(i, 3) - q(i, 2)
          END DO
          arg1 = 3*(ilast-ifirst+1)
          CALL PUSHREALARRAY_ADM(br(:, 0:2), (ilast-ifirst+1)*3)
          CALL PUSHREALARRAY_ADM(bl(:, 0:2), (ilast-ifirst+1)*3)
          CALL PERT_PPM(arg1, q(ifirst:ilast, 0:2), bl(ifirst:ilast, 0:2&
&                 ), br(ifirst:ilast, 0:2), 1)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            bl(i, npy-2) = al(i, npy-2) - q(i, npy-2)
            xt = s15*q(i, npy-1) + s11*q(i, npy-2) + s14*dm(i, npy-2)
            br(i, npy-2) = xt - q(i, npy-2)
            bl(i, npy-1) = xt - q(i, npy-1)
            xt = 0.5*(((2.*dya(i, npy-1)+dya(i, npy-2))*q(i, npy-1)-dya(&
&             i, npy-1)*q(i, npy-2))/(dya(i, npy-2)+dya(i, npy-1))+((2.*&
&             dya(i, npy)+dya(i, npy+1))*q(i, npy)-dya(i, npy)*q(i, npy+&
&             1))/(dya(i, npy)+dya(i, npy+1)))
            IF (q(i, npy) .GT. q(i, npy+1)) THEN
              z4 = q(i, npy+1)
              CALL PUSHCONTROL1B(0)
            ELSE
              z4 = q(i, npy)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, npy-2) .GT. q(i, npy-1)) THEN
              IF (q(i, npy-1) .GT. z4) THEN
                y12 = z4
                CALL PUSHCONTROL2B(0)
              ELSE
                y12 = q(i, npy-1)
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (q(i, npy-2) .GT. z4) THEN
              y12 = z4
              CALL PUSHCONTROL2B(2)
            ELSE
              y12 = q(i, npy-2)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (xt .LT. y12) THEN
              xt = y12
              CALL PUSHCONTROL1B(0)
            ELSE
              xt = xt
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, npy) .LT. q(i, npy+1)) THEN
              z5 = q(i, npy+1)
              CALL PUSHCONTROL1B(0)
            ELSE
              z5 = q(i, npy)
              CALL PUSHCONTROL1B(1)
            END IF
            IF (q(i, npy-2) .LT. q(i, npy-1)) THEN
              IF (q(i, npy-1) .LT. z5) THEN
                y13 = z5
                CALL PUSHCONTROL2B(0)
              ELSE
                y13 = q(i, npy-1)
                CALL PUSHCONTROL2B(1)
              END IF
            ELSE IF (q(i, npy-2) .LT. z5) THEN
              y13 = z5
              CALL PUSHCONTROL2B(2)
            ELSE
              y13 = q(i, npy-2)
              CALL PUSHCONTROL2B(3)
            END IF
            IF (xt .GT. y13) THEN
              xt = y13
              CALL PUSHCONTROL1B(0)
            ELSE
              xt = xt
              CALL PUSHCONTROL1B(1)
            END IF
!        endif
            br(i, npy-1) = xt - q(i, npy-1)
            bl(i, npy) = xt - q(i, npy)
            br(i, npy) = s11*(q(i, npy+1)-q(i, npy)) - s14*dm(i, npy+1)
          END DO
          arg1 = 3*(ilast-ifirst+1)
          CALL PUSHREALARRAY_ADM(br(:, npy-2:npy), (ilast-ifirst+1)*3)
          CALL PUSHREALARRAY_ADM(bl(:, npy-2:npy), (ilast-ifirst+1)*3)
          CALL PERT_PPM(arg1, q(ifirst:ilast, npy-2:npy), bl(ifirst:&
&                 ilast, npy-2:npy), br(ifirst:ilast, npy-2:npy), 1)
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
      DO j=js,je+1
        DO i=ifirst,ilast
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
        DO i=ilast,ifirst,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            temp2 = bl(i, j) + br(i, j)
            temp_ad23 = (c(i, j)+1.)*flux_ad(i, j)
            temp_ad24 = c(i, j)*temp_ad23
            q_ad(i, j) = q_ad(i, j) + flux_ad(i, j)
            c_ad(i, j) = c_ad(i, j) + temp2*temp_ad23 + (bl(i, j)+c(i, j&
&             )*temp2)*flux_ad(i, j)
            bl_ad(i, j) = bl_ad(i, j) + temp_ad24 + temp_ad23
            br_ad(i, j) = br_ad(i, j) + temp_ad24
            flux_ad(i, j) = 0.0
          ELSE
            temp1 = bl(i, j-1) + br(i, j-1)
            temp_ad21 = (1.-c(i, j))*flux_ad(i, j)
            temp_ad22 = -(c(i, j)*temp_ad21)
            q_ad(i, j-1) = q_ad(i, j-1) + flux_ad(i, j)
            c_ad(i, j) = c_ad(i, j) - temp1*temp_ad21 - (br(i, j-1)-c(i&
&             , j)*temp1)*flux_ad(i, j)
            br_ad(i, j-1) = br_ad(i, j-1) + temp_ad22 + temp_ad21
            bl_ad(i, j-1) = bl_ad(i, j-1) + temp_ad22
            flux_ad(i, j) = 0.0
          END IF
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        dm_ad = 0.0
        al_ad = 0.0
      ELSE
        IF (branch .EQ. 1) THEN
          dm_ad = 0.0
          al_ad = 0.0
        ELSE
          CALL POPREALARRAY_ADM(bl(:, npy-2:npy), (ilast-ifirst+1)*3)
          CALL POPREALARRAY_ADM(br(:, npy-2:npy), (ilast-ifirst+1)*3)
          CALL PERT_PPM_ADM(arg1, q(ifirst:ilast, npy-2:npy), bl(ifirst:&
&                     ilast, npy-2:npy), bl_ad(ifirst:ilast, npy-2:npy)&
&                     , br(ifirst:ilast, npy-2:npy), br_ad(ifirst:ilast&
&                     , npy-2:npy), 1)
          dm_ad = 0.0
          al_ad = 0.0
          DO i=ilast,ifirst,-1
            q_ad(i, npy+1) = q_ad(i, npy+1) + s11*br_ad(i, npy)
            q_ad(i, npy) = q_ad(i, npy) - bl_ad(i, npy) - s11*br_ad(i, &
&             npy)
            dm_ad(i, npy+1) = dm_ad(i, npy+1) - s14*br_ad(i, npy)
            br_ad(i, npy) = 0.0
            xt_ad = br_ad(i, npy-1) + bl_ad(i, npy)
            bl_ad(i, npy) = 0.0
            q_ad(i, npy-1) = q_ad(i, npy-1) - br_ad(i, npy-1)
            br_ad(i, npy-1) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y13_ad = xt_ad
              xt_ad = 0.0
            ELSE
              y13_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                z5_ad = y13_ad
              ELSE
                q_ad(i, npy-1) = q_ad(i, npy-1) + y13_ad
                z5_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              z5_ad = y13_ad
            ELSE
              q_ad(i, npy-2) = q_ad(i, npy-2) + y13_ad
              z5_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q_ad(i, npy+1) = q_ad(i, npy+1) + z5_ad
            ELSE
              q_ad(i, npy) = q_ad(i, npy) + z5_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y12_ad = xt_ad
              xt_ad = 0.0
            ELSE
              y12_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                z4_ad = y12_ad
              ELSE
                q_ad(i, npy-1) = q_ad(i, npy-1) + y12_ad
                z4_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              z4_ad = y12_ad
            ELSE
              q_ad(i, npy-2) = q_ad(i, npy-2) + y12_ad
              z4_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q_ad(i, npy+1) = q_ad(i, npy+1) + z4_ad
            ELSE
              q_ad(i, npy) = q_ad(i, npy) + z4_ad
            END IF
            temp_ad19 = 0.5*xt_ad/(dya(i, npy-2)+dya(i, npy-1))
            temp_ad20 = 0.5*xt_ad/(dya(i, npy)+dya(i, npy+1))
            q_ad(i, npy-1) = q_ad(i, npy-1) + (dya(i, npy-1)*2.+dya(i, &
&             npy-2))*temp_ad19
            q_ad(i, npy-2) = q_ad(i, npy-2) - dya(i, npy-1)*temp_ad19
            q_ad(i, npy) = q_ad(i, npy) + (dya(i, npy)*2.+dya(i, npy+1))&
&             *temp_ad20
            q_ad(i, npy+1) = q_ad(i, npy+1) - dya(i, npy)*temp_ad20
            xt_ad = br_ad(i, npy-2) + bl_ad(i, npy-1)
            q_ad(i, npy-1) = q_ad(i, npy-1) - bl_ad(i, npy-1)
            bl_ad(i, npy-1) = 0.0
            q_ad(i, npy-2) = q_ad(i, npy-2) - br_ad(i, npy-2)
            br_ad(i, npy-2) = 0.0
            q_ad(i, npy-1) = q_ad(i, npy-1) + s15*xt_ad
            q_ad(i, npy-2) = q_ad(i, npy-2) + s11*xt_ad - bl_ad(i, npy-2&
&             )
            dm_ad(i, npy-2) = dm_ad(i, npy-2) + s14*xt_ad
            al_ad(i, npy-2) = al_ad(i, npy-2) + bl_ad(i, npy-2)
            bl_ad(i, npy-2) = 0.0
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          arg1 = 3*(ilast-ifirst+1)
          CALL POPREALARRAY_ADM(bl(:, 0:2), (ilast-ifirst+1)*3)
          CALL POPREALARRAY_ADM(br(:, 0:2), (ilast-ifirst+1)*3)
          CALL PERT_PPM_ADM(arg1, q(ifirst:ilast, 0:2), bl(ifirst:ilast&
&                     , 0:2), bl_ad(ifirst:ilast, 0:2), br(ifirst:ilast&
&                     , 0:2), br_ad(ifirst:ilast, 0:2), 1)
          DO i=ilast,ifirst,-1
            al_ad(i, 3) = al_ad(i, 3) + br_ad(i, 2)
            q_ad(i, 2) = q_ad(i, 2) - bl_ad(i, 2) - br_ad(i, 2)
            br_ad(i, 2) = 0.0
            xt_ad = br_ad(i, 1) + bl_ad(i, 2)
            bl_ad(i, 2) = 0.0
            q_ad(i, 1) = q_ad(i, 1) + s15*xt_ad - br_ad(i, 1)
            br_ad(i, 1) = 0.0
            q_ad(i, 2) = q_ad(i, 2) + s11*xt_ad
            dm_ad(i, 2) = dm_ad(i, 2) - s14*xt_ad
            xt_ad = br_ad(i, 0) + bl_ad(i, 1)
            q_ad(i, 1) = q_ad(i, 1) - bl_ad(i, 1)
            bl_ad(i, 1) = 0.0
            q_ad(i, 0) = q_ad(i, 0) - br_ad(i, 0)
            br_ad(i, 0) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y11_ad = xt_ad
              xt_ad = 0.0
            ELSE
              y11_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                z3_ad = y11_ad
              ELSE
                q_ad(i, 0) = q_ad(i, 0) + y11_ad
                z3_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              z3_ad = y11_ad
            ELSE
              q_ad(i, -1) = q_ad(i, -1) + y11_ad
              z3_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q_ad(i, 2) = q_ad(i, 2) + z3_ad
            ELSE
              q_ad(i, 1) = q_ad(i, 1) + z3_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y10_ad = xt_ad
              xt_ad = 0.0
            ELSE
              y10_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                z2_ad = y10_ad
              ELSE
                q_ad(i, 0) = q_ad(i, 0) + y10_ad
                z2_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              z2_ad = y10_ad
            ELSE
              q_ad(i, -1) = q_ad(i, -1) + y10_ad
              z2_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q_ad(i, 2) = q_ad(i, 2) + z2_ad
            ELSE
              q_ad(i, 1) = q_ad(i, 1) + z2_ad
            END IF
            temp_ad17 = 0.5*xt_ad/(dya(i, -1)+dya(i, 0))
            temp_ad18 = 0.5*xt_ad/(dya(i, 1)+dya(i, 2))
            q_ad(i, 0) = q_ad(i, 0) + (dya(i, 0)*2.+dya(i, -1))*&
&             temp_ad17
            q_ad(i, -1) = q_ad(i, -1) - dya(i, 0)*temp_ad17
            q_ad(i, 1) = q_ad(i, 1) + (dya(i, 1)*2.+dya(i, 2))*temp_ad18
            q_ad(i, 2) = q_ad(i, 2) - dya(i, 1)*temp_ad18
            dm_ad(i, -1) = dm_ad(i, -1) + s14*bl_ad(i, 0)
            q_ad(i, -1) = q_ad(i, -1) + s11*bl_ad(i, 0)
            q_ad(i, 0) = q_ad(i, 0) - s11*bl_ad(i, 0)
            bl_ad(i, 0) = 0.0
          END DO
        END IF
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je1,js1,-1
          arg1 = ilast - ifirst + 1
          CALL POPREALARRAY_ADM(bl(:, j), ilast - ifirst + 1)
          CALL POPREALARRAY_ADM(br(:, j), ilast - ifirst + 1)
          CALL PERT_PPM_ADM(arg1, q(ifirst:ilast, j), bl(ifirst:ilast, j&
&                     ), bl_ad(ifirst:ilast, j), br(ifirst:ilast, j), &
&                     br_ad(ifirst:ilast, j), 0)
        END DO
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je1,js1,-1
          DO i=ilast,ifirst,-1
            xt = 2.*dm(i, j)
            min5_ad = SIGN(1.d0, min5*xt)*br_ad(i, j)
            br_ad(i, j) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(min5)
              y5_ad = min5_ad
              x5_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min5)
              x5_ad = min5_ad
              y5_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              al_ad(i, j+1) = al_ad(i, j+1) + y5_ad
              q_ad(i, j) = q_ad(i, j) - y5_ad
            ELSE
              q_ad(i, j) = q_ad(i, j) + y5_ad
              al_ad(i, j+1) = al_ad(i, j+1) - y5_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt_ad = x5_ad
            ELSE
              xt_ad = -x5_ad
            END IF
            min4_ad = -(SIGN(1.d0, min4*xt)*bl_ad(i, j))
            bl_ad(i, j) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(min4)
              y4_ad = min4_ad
              x4_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min4)
              x4_ad = min4_ad
              y4_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              al_ad(i, j) = al_ad(i, j) + y4_ad
              q_ad(i, j) = q_ad(i, j) - y4_ad
            ELSE
              q_ad(i, j) = q_ad(i, j) + y4_ad
              al_ad(i, j) = al_ad(i, j) - y4_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt_ad = xt_ad + x4_ad
            ELSE
              xt_ad = xt_ad - x4_ad
            END IF
            dm_ad(i, j) = dm_ad(i, j) + 2.*xt_ad
          END DO
        END DO
      ELSE IF (branch .EQ. 1) THEN
        DO j=je1,js1,-1
          DO i=ilast,ifirst,-1
            xt = ppm_fac*dm(i, j)
            min7_ad = SIGN(1.d0, min7*xt)*br_ad(i, j)
            br_ad(i, j) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(min7)
              y7_ad = min7_ad
              x7_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min7)
              x7_ad = min7_ad
              y7_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              al_ad(i, j+1) = al_ad(i, j+1) + y7_ad
              q_ad(i, j) = q_ad(i, j) - y7_ad
            ELSE
              q_ad(i, j) = q_ad(i, j) + y7_ad
              al_ad(i, j+1) = al_ad(i, j+1) - y7_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt_ad = x7_ad
            ELSE
              xt_ad = -x7_ad
            END IF
            min6_ad = -(SIGN(1.d0, min6*xt)*bl_ad(i, j))
            bl_ad(i, j) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY_ADM(min6)
              y6_ad = min6_ad
              x6_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min6)
              x6_ad = min6_ad
              y6_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              al_ad(i, j) = al_ad(i, j) + y6_ad
              q_ad(i, j) = q_ad(i, j) - y6_ad
            ELSE
              q_ad(i, j) = q_ad(i, j) + y6_ad
              al_ad(i, j) = al_ad(i, j) - y6_ad
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              xt_ad = xt_ad + x6_ad
            ELSE
              xt_ad = xt_ad - x6_ad
            END IF
            dm_ad(i, j) = dm_ad(i, j) + ppm_fac*xt_ad
          END DO
        END DO
      ELSE
        dq_ad = 0.0
        DO j=je1,js1,-1
          DO i=ilast,ifirst,-1
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                GOTO 110
              ELSE
                y9_ad = bl_ad(i, j)
                bl_ad(i, j) = 0.0
                x9_ad = 0.0
              END IF
            ELSE IF (branch .EQ. 2) THEN
              x9_ad = bl_ad(i, j)
              bl_ad(i, j) = 0.0
              y9_ad = 0.0
            ELSE
              br_ad(i, j) = 0.0
              bl_ad(i, j) = 0.0
              GOTO 110
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y15_ad = y9_ad
            ELSE
              bl_ad(i, j) = bl_ad(i, j) + y9_ad
              y15_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_1_ad = y15_ad
                pmp_1_ad = 0.0
              ELSE
                pmp_1_ad = y15_ad
                lac_1_ad = 0.0
              END IF
            ELSE
              IF (branch .EQ. 2) THEN
                lac_1_ad = y15_ad
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
            dq_ad(i, j+1) = dq_ad(i, j+1) + 0.75*lac_1_ad
            dq_ad(i, j) = dq_ad(i, j) - pmp_1_ad
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y8_ad = br_ad(i, j)
              br_ad(i, j) = 0.0
              x8_ad = 0.0
            ELSE
              x8_ad = br_ad(i, j)
              br_ad(i, j) = 0.0
              y8_ad = 0.0
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y14_ad = y8_ad
            ELSE
              br_ad(i, j) = br_ad(i, j) + y8_ad
              y14_ad = 0.0
            END IF
            CALL POPCONTROL2B(branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                lac_2_ad = y14_ad
                pmp_2_ad = 0.0
              ELSE
                pmp_2_ad = y14_ad
                lac_2_ad = 0.0
              END IF
            ELSE
              IF (branch .EQ. 2) THEN
                lac_2_ad = y14_ad
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
            dq_ad(i, j-2) = dq_ad(i, j-2) - 0.75*lac_2_ad
            dq_ad(i, j-1) = dq_ad(i, j-1) + pmp_2_ad
 110        al_ad(i, j+1) = al_ad(i, j+1) + br_ad(i, j)
            q_ad(i, j) = q_ad(i, j) - bl_ad(i, j) - br_ad(i, j)
            br_ad(i, j) = 0.0
            al_ad(i, j) = al_ad(i, j) + bl_ad(i, j)
            bl_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je1+1,js1-2,-1
          DO i=ilast,ifirst,-1
            q_ad(i, j+1) = q_ad(i, j+1) + 2.*dq_ad(i, j)
            q_ad(i, j) = q_ad(i, j) - 2.*dq_ad(i, j)
            dq_ad(i, j) = 0.0
          END DO
        END DO
      END IF
      DO j=je1+1,js1,-1
        DO i=ilast,ifirst,-1
          q_ad(i, j-1) = q_ad(i, j-1) + 0.5*al_ad(i, j)
          q_ad(i, j) = q_ad(i, j) + 0.5*al_ad(i, j)
          dm_ad(i, j-1) = dm_ad(i, j-1) + r3*al_ad(i, j)
          dm_ad(i, j) = dm_ad(i, j) - r3*al_ad(i, j)
          al_ad(i, j) = 0.0
        END DO
      END DO
      DO j=je+2,js-2,-1
        DO i=ilast,ifirst,-1
          xt = 0.25*(q(i, j+1)-q(i, j-1))
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
            x3_ad = 0.0
          ELSE
            IF (branch .EQ. 2) THEN
              CALL POPREALARRAY_ADM(min3)
              z1_ad = min3_ad
              x3_ad = 0.0
            ELSE
              CALL POPREALARRAY_ADM(min3)
              x3_ad = min3_ad
              z1_ad = 0.0
            END IF
            y3_ad = 0.0
          END IF
          q_ad(i, j) = q_ad(i, j) + z1_ad
          min8_ad = -z1_ad
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              q_ad(i, j+1) = q_ad(i, j+1) + min8_ad
            ELSE
              q_ad(i, j) = q_ad(i, j) + min8_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            q_ad(i, j+1) = q_ad(i, j+1) + min8_ad
          ELSE
            q_ad(i, j-1) = q_ad(i, j-1) + min8_ad
          END IF
          max1_ad = y3_ad
          q_ad(i, j) = q_ad(i, j) - y3_ad
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              q_ad(i, j+1) = q_ad(i, j+1) + max1_ad
            ELSE
              q_ad(i, j) = q_ad(i, j) + max1_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            q_ad(i, j+1) = q_ad(i, j+1) + max1_ad
          ELSE
            q_ad(i, j-1) = q_ad(i, j-1) + max1_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            xt_ad = x3_ad
          ELSE
            xt_ad = -x3_ad
          END IF
          q_ad(i, j+1) = q_ad(i, j+1) + 0.25*xt_ad
          q_ad(i, j-1) = q_ad(i, j-1) - 0.25*xt_ad
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
  END SUBROUTINE YPPM_ADM
  SUBROUTINE YPPM(flux, q, c, jord, ifirst, ilast, isd, ied, js, je, jsd&
&   , jed, npx, npy, dya, nested, grid_type)
    IMPLICIT NONE
! Compute domain
    INTEGER, INTENT(IN) :: ifirst, ilast
    INTEGER, INTENT(IN) :: isd, ied, js, je, jsd, jed
    INTEGER, INTENT(IN) :: jord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst:ilast, jsd:jed)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
!  Flux
    REAL, INTENT(OUT) :: flux(ifirst:ilast, js:je+1)
    REAL, INTENT(IN) :: dya(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
! Local:
    REAL :: dm(ifirst:ilast, js-2:je+2)
    REAL :: al(ifirst:ilast, js-1:je+2)
    REAL, DIMENSION(ifirst:ilast, js-1:je+1) :: bl, br, b0
    REAL :: dq(ifirst:ilast, js-3:je+2)
    REAL, DIMENSION(ifirst:ilast) :: fx0, fx1
    LOGICAL, DIMENSION(ifirst:ilast, js-1:je+1) :: smt5, smt6
    REAL :: x0, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2, r1
    INTEGER :: i, j, js1, je3, je1
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: min1
    REAL :: min2
    REAL :: abs0
    REAL :: abs1
    REAL :: min3
    REAL :: min4
    REAL :: min5
    REAL :: min6
    REAL :: min7
    REAL :: abs2
    REAL :: abs3
    REAL :: abs4
    REAL :: max1
    REAL :: min8
    REAL :: abs5
    REAL :: abs6
    REAL :: abs7
    INTEGER :: arg1
    REAL :: x9
    REAL :: x8
    REAL :: x7
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: x3
    REAL :: x2
    REAL :: x1
    REAL :: y15
    REAL :: y14
    REAL :: y13
    REAL :: y12
    REAL :: y11
    REAL :: y10
    REAL :: z5
    REAL :: z4
    REAL :: z3
    REAL :: z2
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
    IF (.NOT.nested .AND. grid_type .LT. 3) THEN
      IF (3 .LT. js - 1) THEN
        js1 = js - 1
      ELSE
        js1 = 3
      END IF
      IF (npy - 2 .GT. je + 2) THEN
        je3 = je + 2
      ELSE
        je3 = npy - 2
      END IF
      IF (npy - 3 .GT. je + 1) THEN
        je1 = je + 1
      ELSE
        je1 = npy - 3
      END IF
    ELSE
! Nested grid OR Doubly periodic domain:
      js1 = js - 1
      je3 = je + 2
      je1 = je + 1
    END IF
    IF (jord .LT. 8 .OR. jord .EQ. 333) THEN
      DO j=js1,je3
        DO i=ifirst,ilast
          al(i, j) = p1*(q(i, j-1)+q(i, j)) + p2*(q(i, j-2)+q(i, j+1))
        END DO
      END DO
      IF (jord .EQ. 7) THEN
        DO j=js1,je3
          DO i=ifirst,ilast
            IF (al(i, j) .LT. 0.) al(i, j) = 0.5*(q(i, j)+q(i, j+1))
          END DO
        END DO
      END IF
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            al(i, 0) = c1*q(i, -2) + c2*q(i, -1) + c3*q(i, 0)
            al(i, 1) = 0.5*(((2.*dya(i, 0)+dya(i, -1))*q(i, 0)-dya(i, 0)&
&             *q(i, -1))/(dya(i, -1)+dya(i, 0))+((2.*dya(i, 1)+dya(i, 2)&
&             )*q(i, 1)-dya(i, 1)*q(i, 2))/(dya(i, 1)+dya(i, 2)))
            al(i, 2) = c3*q(i, 1) + c2*q(i, 2) + c1*q(i, 3)
          END DO
          IF (jord .EQ. 7) THEN
            DO i=ifirst,ilast
              IF (0. .LT. al(i, 0)) THEN
                al(i, 0) = al(i, 0)
              ELSE
                al(i, 0) = 0.
              END IF
              IF (0. .LT. al(i, 1)) THEN
                al(i, 1) = al(i, 1)
              ELSE
                al(i, 1) = 0.
              END IF
              IF (0. .LT. al(i, 2)) THEN
                al(i, 2) = al(i, 2)
              ELSE
                al(i, 2) = 0.
              END IF
            END DO
          END IF
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            al(i, npy-1) = c1*q(i, npy-3) + c2*q(i, npy-2) + c3*q(i, npy&
&             -1)
            al(i, npy) = 0.5*(((2.*dya(i, npy-1)+dya(i, npy-2))*q(i, npy&
&             -1)-dya(i, npy-1)*q(i, npy-2))/(dya(i, npy-2)+dya(i, npy-1&
&             ))+((2.*dya(i, npy)+dya(i, npy+1))*q(i, npy)-dya(i, npy)*q&
&             (i, npy+1))/(dya(i, npy)+dya(i, npy+1)))
            al(i, npy+1) = c3*q(i, npy) + c2*q(i, npy+1) + c1*q(i, npy+2&
&             )
          END DO
          IF (jord .EQ. 7) THEN
            DO i=ifirst,ilast
              IF (0. .LT. al(i, npy-1)) THEN
                al(i, npy-1) = al(i, npy-1)
              ELSE
                al(i, npy-1) = 0.
              END IF
              IF (0. .LT. al(i, npy)) THEN
                al(i, npy) = al(i, npy)
              ELSE
                al(i, npy) = 0.
              END IF
              IF (0. .LT. al(i, npy+1)) THEN
                al(i, npy+1) = al(i, npy+1)
              ELSE
                al(i, npy+1) = 0.
              END IF
            END DO
          END IF
        END IF
      END IF
      IF (jord .EQ. 1) THEN
        DO j=js,je+1
          DO i=ifirst,ilast
            IF (c(i, j) .GT. 0.) THEN
              flux(i, j) = q(i, j-1)
            ELSE
              flux(i, j) = q(i, j)
            END IF
          END DO
        END DO
      ELSE IF (jord .EQ. 2) THEN
! Perfectly linear scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              qtmp = q(i, j-1)
              flux(i, j) = qtmp + (1.-xt)*(al(i, j)-qtmp-xt*(al(i, j-1)+&
&               al(i, j)-(qtmp+qtmp)))
            ELSE
              qtmp = q(i, j)
              flux(i, j) = qtmp + (1.+xt)*(al(i, j)-qtmp+xt*(al(i, j)+al&
&               (i, j+1)-(qtmp+qtmp)))
            END IF
          END DO
        END DO
      ELSE IF (jord .EQ. 333) THEN
! Perfectly linear scheme, more diffusive than ord=2 (HoldawayKent-2015-TellusA)
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              flux(i, j) = (2.0*q(i, j)+5.0*q(i, j-1)-q(i, j-2))/6.0 - &
&               0.5*xt*(q(i, j)-q(i, j-1)) + xt*xt/6.0*(q(i, j)-2.0*q(i&
&               , j-1)+q(i, j-2))
            ELSE
              flux(i, j) = (2.0*q(i, j-1)+5.0*q(i, j)-q(i, j+1))/6.0 - &
&               0.5*xt*(q(i, j)-q(i, j-1)) + xt*xt/6.0*(q(i, j+1)-2.0*q(&
&               i, j)+q(i, j-1))
            END IF
          END DO
        END DO
      ELSE IF (jord .EQ. 3) THEN
        DO j=js-1,je+1
          DO i=ifirst,ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j+1) - q(i, j)
            b0(i, j) = bl(i, j) + br(i, j)
            IF (b0(i, j) .GE. 0.) THEN
              x0 = b0(i, j)
            ELSE
              x0 = -b0(i, j)
            END IF
            IF (bl(i, j) - br(i, j) .GE. 0.) THEN
              xt = bl(i, j) - br(i, j)
            ELSE
              xt = -(bl(i, j)-br(i, j))
            END IF
            smt5(i, j) = x0 .LT. xt
            smt6(i, j) = 3.*x0 .LT. xt
          END DO
        END DO
        DO j=js,je+1
          DO i=ifirst,ilast
            fx1(i) = 0.
          END DO
          DO i=ifirst,ilast
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              fx0(i) = q(i, j-1)
              IF (smt6(i, j-1) .OR. smt5(i, j)) THEN
                fx1(i) = br(i, j-1) - xt*b0(i, j-1)
              ELSE IF (smt5(i, j-1)) THEN
                IF (bl(i, j-1) .GE. 0.) THEN
                  x1 = bl(i, j-1)
                ELSE
                  x1 = -bl(i, j-1)
                END IF
                IF (br(i, j-1) .GE. 0.) THEN
                  y1 = br(i, j-1)
                ELSE
                  y1 = -br(i, j-1)
                END IF
                IF (x1 .GT. y1) THEN
                  min1 = y1
                ELSE
                  min1 = x1
                END IF
! both up-downwind sides are noisy; 2nd order, piece-wise linear
                fx1(i) = SIGN(min1, br(i, j-1))
              END IF
            ELSE
              fx0(i) = q(i, j)
              IF (smt6(i, j) .OR. smt5(i, j-1)) THEN
                fx1(i) = bl(i, j) + xt*b0(i, j)
              ELSE IF (smt5(i, j)) THEN
                IF (bl(i, j) .GE. 0.) THEN
                  x2 = bl(i, j)
                ELSE
                  x2 = -bl(i, j)
                END IF
                IF (br(i, j) .GE. 0.) THEN
                  y2 = br(i, j)
                ELSE
                  y2 = -br(i, j)
                END IF
                IF (x2 .GT. y2) THEN
                  min2 = y2
                ELSE
                  min2 = x2
                END IF
                fx1(i) = SIGN(min2, bl(i, j))
              END IF
            END IF
            IF (xt .GE. 0.) THEN
              abs0 = xt
            ELSE
              abs0 = -xt
            END IF
            flux(i, j) = fx0(i) + (1.-abs0)*fx1(i)
          END DO
        END DO
      ELSE IF (jord .EQ. 4) THEN
        DO j=js-1,je+1
          DO i=ifirst,ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j+1) - q(i, j)
            b0(i, j) = bl(i, j) + br(i, j)
            IF (b0(i, j) .GE. 0.) THEN
              x0 = b0(i, j)
            ELSE
              x0 = -b0(i, j)
            END IF
            IF (bl(i, j) - br(i, j) .GE. 0.) THEN
              xt = bl(i, j) - br(i, j)
            ELSE
              xt = -(bl(i, j)-br(i, j))
            END IF
            smt5(i, j) = x0 .LT. xt
            smt6(i, j) = 3.*x0 .LT. xt
          END DO
        END DO
        DO j=js,je+1
          DO i=ifirst,ilast
            fx1(i) = 0.
          END DO
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
            IF (c(i, j) .GT. 0.) THEN
              fx0(i) = q(i, j-1)
              IF (smt6(i, j-1) .OR. smt5(i, j)) fx1(i) = (1.-c(i, j))*(&
&                 br(i, j-1)-c(i, j)*b0(i, j-1))
            ELSE
              fx0(i) = q(i, j)
              IF (smt6(i, j) .OR. smt5(i, j-1)) fx1(i) = (1.+c(i, j))*(&
&                 bl(i, j)+c(i, j)*b0(i, j))
            END IF
            flux(i, j) = fx0(i) + fx1(i)
          END DO
        END DO
      ELSE
! jord=5,6,7
        IF (jord .EQ. 5) THEN
          DO j=js-1,je+1
            DO i=ifirst,ilast
              bl(i, j) = al(i, j) - q(i, j)
              br(i, j) = al(i, j+1) - q(i, j)
              b0(i, j) = bl(i, j) + br(i, j)
              smt5(i, j) = bl(i, j)*br(i, j) .LT. 0.
            END DO
          END DO
        ELSE
          DO j=js-1,je+1
            DO i=ifirst,ilast
              bl(i, j) = al(i, j) - q(i, j)
              br(i, j) = al(i, j+1) - q(i, j)
              b0(i, j) = bl(i, j) + br(i, j)
              IF (3.*b0(i, j) .GE. 0.) THEN
                abs1 = 3.*b0(i, j)
              ELSE
                abs1 = -(3.*b0(i, j))
              END IF
              IF (bl(i, j) - br(i, j) .GE. 0.) THEN
                abs4 = bl(i, j) - br(i, j)
              ELSE
                abs4 = -(bl(i, j)-br(i, j))
              END IF
              smt5(i, j) = abs1 .LT. abs4
            END DO
          END DO
        END IF
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
            IF (c(i, j) .GT. 0.) THEN
              fx1(i) = (1.-c(i, j))*(br(i, j-1)-c(i, j)*b0(i, j-1))
              flux(i, j) = q(i, j-1)
            ELSE
              fx1(i) = (1.+c(i, j))*(bl(i, j)+c(i, j)*b0(i, j))
              flux(i, j) = q(i, j)
            END IF
            IF (smt5(i, j-1) .OR. smt5(i, j)) flux(i, j) = flux(i, j) + &
&               fx1(i)
          END DO
        END DO
      END IF
      RETURN
    ELSE
! Monotonic constraints:
! ord = 8: PPM with Lin's PPM fast monotone constraint
! ord > 8: PPM with Lin's modification of Huynh 2nd constraint
      DO j=js-2,je+2
        DO i=ifirst,ilast
          xt = 0.25*(q(i, j+1)-q(i, j-1))
          IF (xt .GE. 0.) THEN
            x3 = xt
          ELSE
            x3 = -xt
          END IF
          IF (q(i, j-1) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i, j+1)) THEN
              max1 = q(i, j+1)
            ELSE
              max1 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .LT. q(i, j+1)) THEN
            max1 = q(i, j+1)
          ELSE
            max1 = q(i, j-1)
          END IF
          y3 = max1 - q(i, j)
          IF (q(i, j-1) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i, j+1)) THEN
              min8 = q(i, j+1)
            ELSE
              min8 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .GT. q(i, j+1)) THEN
            min8 = q(i, j+1)
          ELSE
            min8 = q(i, j-1)
          END IF
          z1 = q(i, j) - min8
          IF (x3 .GT. y3) THEN
            IF (y3 .GT. z1) THEN
              min3 = z1
            ELSE
              min3 = y3
            END IF
          ELSE IF (x3 .GT. z1) THEN
            min3 = z1
          ELSE
            min3 = x3
          END IF
          dm(i, j) = SIGN(min3, xt)
        END DO
      END DO
      DO j=js1,je1+1
        DO i=ifirst,ilast
          al(i, j) = 0.5*(q(i, j-1)+q(i, j)) + r3*(dm(i, j-1)-dm(i, j))
        END DO
      END DO
      IF (jord .EQ. 8) THEN
        DO j=js1,je1
          DO i=ifirst,ilast
            xt = 2.*dm(i, j)
            IF (xt .GE. 0.) THEN
              x4 = xt
            ELSE
              x4 = -xt
            END IF
            IF (al(i, j) - q(i, j) .GE. 0.) THEN
              y4 = al(i, j) - q(i, j)
            ELSE
              y4 = -(al(i, j)-q(i, j))
            END IF
            IF (x4 .GT. y4) THEN
              min4 = y4
            ELSE
              min4 = x4
            END IF
            bl(i, j) = -SIGN(min4, xt)
            IF (xt .GE. 0.) THEN
              x5 = xt
            ELSE
              x5 = -xt
            END IF
            IF (al(i, j+1) - q(i, j) .GE. 0.) THEN
              y5 = al(i, j+1) - q(i, j)
            ELSE
              y5 = -(al(i, j+1)-q(i, j))
            END IF
            IF (x5 .GT. y5) THEN
              min5 = y5
            ELSE
              min5 = x5
            END IF
            br(i, j) = SIGN(min5, xt)
          END DO
        END DO
      ELSE IF (jord .EQ. 11) THEN
        DO j=js1,je1
          DO i=ifirst,ilast
            xt = ppm_fac*dm(i, j)
            IF (xt .GE. 0.) THEN
              x6 = xt
            ELSE
              x6 = -xt
            END IF
            IF (al(i, j) - q(i, j) .GE. 0.) THEN
              y6 = al(i, j) - q(i, j)
            ELSE
              y6 = -(al(i, j)-q(i, j))
            END IF
            IF (x6 .GT. y6) THEN
              min6 = y6
            ELSE
              min6 = x6
            END IF
            bl(i, j) = -SIGN(min6, xt)
            IF (xt .GE. 0.) THEN
              x7 = xt
            ELSE
              x7 = -xt
            END IF
            IF (al(i, j+1) - q(i, j) .GE. 0.) THEN
              y7 = al(i, j+1) - q(i, j)
            ELSE
              y7 = -(al(i, j+1)-q(i, j))
            END IF
            IF (x7 .GT. y7) THEN
              min7 = y7
            ELSE
              min7 = x7
            END IF
            br(i, j) = SIGN(min7, xt)
          END DO
        END DO
      ELSE
        DO j=js1-2,je1+1
          DO i=ifirst,ilast
            dq(i, j) = 2.*(q(i, j+1)-q(i, j))
          END DO
        END DO
        DO j=js1,je1
          DO i=ifirst,ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j+1) - q(i, j)
            IF (dm(i, j-1) .GE. 0.) THEN
              abs2 = dm(i, j-1)
            ELSE
              abs2 = -dm(i, j-1)
            END IF
            IF (dm(i, j) .GE. 0.) THEN
              abs5 = dm(i, j)
            ELSE
              abs5 = -dm(i, j)
            END IF
            IF (dm(i, j+1) .GE. 0.) THEN
              abs7 = dm(i, j+1)
            ELSE
              abs7 = -dm(i, j+1)
            END IF
            IF (abs2 + abs5 + abs7 .LT. near_zero) THEN
              bl(i, j) = 0.
              br(i, j) = 0.
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
                pmp_2 = dq(i, j-1)
                lac_2 = pmp_2 - 0.75*dq(i, j-2)
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
                    y14 = lac_2
                  ELSE
                    y14 = pmp_2
                  END IF
                ELSE IF (0. .GT. lac_2) THEN
                  y14 = lac_2
                ELSE
                  y14 = 0.
                END IF
                IF (br(i, j) .LT. y14) THEN
                  y8 = y14
                ELSE
                  y8 = br(i, j)
                END IF
                IF (x8 .GT. y8) THEN
                  br(i, j) = y8
                ELSE
                  br(i, j) = x8
                END IF
                pmp_1 = -dq(i, j)
                lac_1 = pmp_1 + 0.75*dq(i, j+1)
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
                    y15 = lac_1
                  ELSE
                    y15 = pmp_1
                  END IF
                ELSE IF (0. .GT. lac_1) THEN
                  y15 = lac_1
                ELSE
                  y15 = 0.
                END IF
                IF (bl(i, j) .LT. y15) THEN
                  y9 = y15
                ELSE
                  y9 = bl(i, j)
                END IF
                IF (x9 .GT. y9) THEN
                  bl(i, j) = y9
                ELSE
                  bl(i, j) = x9
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
      IF (jord .EQ. 9 .OR. jord .EQ. 13) THEN
! Positive definite constraint:
        DO j=js1,je1
          arg1 = ilast - ifirst + 1
          CALL PERT_PPM(arg1, q(ifirst:ilast, j), bl(ifirst:ilast, j), &
&                 br(ifirst:ilast, j), 0)
        END DO
      END IF
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            bl(i, 0) = s14*dm(i, -1) + s11*(q(i, -1)-q(i, 0))
            xt = 0.5*(((2.*dya(i, 0)+dya(i, -1))*q(i, 0)-dya(i, 0)*q(i, &
&             -1))/(dya(i, -1)+dya(i, 0))+((2.*dya(i, 1)+dya(i, 2))*q(i&
&             , 1)-dya(i, 1)*q(i, 2))/(dya(i, 1)+dya(i, 2)))
            IF (q(i, 1) .GT. q(i, 2)) THEN
              z2 = q(i, 2)
            ELSE
              z2 = q(i, 1)
            END IF
            IF (q(i, -1) .GT. q(i, 0)) THEN
              IF (q(i, 0) .GT. z2) THEN
                y10 = z2
              ELSE
                y10 = q(i, 0)
              END IF
            ELSE IF (q(i, -1) .GT. z2) THEN
              y10 = z2
            ELSE
              y10 = q(i, -1)
            END IF
            IF (xt .LT. y10) THEN
              xt = y10
            ELSE
              xt = xt
            END IF
            IF (q(i, 1) .LT. q(i, 2)) THEN
              z3 = q(i, 2)
            ELSE
              z3 = q(i, 1)
            END IF
            IF (q(i, -1) .LT. q(i, 0)) THEN
              IF (q(i, 0) .LT. z3) THEN
                y11 = z3
              ELSE
                y11 = q(i, 0)
              END IF
            ELSE IF (q(i, -1) .LT. z3) THEN
              y11 = z3
            ELSE
              y11 = q(i, -1)
            END IF
            IF (xt .GT. y11) THEN
              xt = y11
            ELSE
              xt = xt
            END IF
!        endif
            br(i, 0) = xt - q(i, 0)
            bl(i, 1) = xt - q(i, 1)
            xt = s15*q(i, 1) + s11*q(i, 2) - s14*dm(i, 2)
            br(i, 1) = xt - q(i, 1)
            bl(i, 2) = xt - q(i, 2)
            br(i, 2) = al(i, 3) - q(i, 2)
          END DO
          arg1 = 3*(ilast-ifirst+1)
          CALL PERT_PPM(arg1, q(ifirst:ilast, 0:2), bl(ifirst:ilast, 0:2&
&                 ), br(ifirst:ilast, 0:2), 1)
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            bl(i, npy-2) = al(i, npy-2) - q(i, npy-2)
            xt = s15*q(i, npy-1) + s11*q(i, npy-2) + s14*dm(i, npy-2)
            br(i, npy-2) = xt - q(i, npy-2)
            bl(i, npy-1) = xt - q(i, npy-1)
            xt = 0.5*(((2.*dya(i, npy-1)+dya(i, npy-2))*q(i, npy-1)-dya(&
&             i, npy-1)*q(i, npy-2))/(dya(i, npy-2)+dya(i, npy-1))+((2.*&
&             dya(i, npy)+dya(i, npy+1))*q(i, npy)-dya(i, npy)*q(i, npy+&
&             1))/(dya(i, npy)+dya(i, npy+1)))
            IF (q(i, npy) .GT. q(i, npy+1)) THEN
              z4 = q(i, npy+1)
            ELSE
              z4 = q(i, npy)
            END IF
            IF (q(i, npy-2) .GT. q(i, npy-1)) THEN
              IF (q(i, npy-1) .GT. z4) THEN
                y12 = z4
              ELSE
                y12 = q(i, npy-1)
              END IF
            ELSE IF (q(i, npy-2) .GT. z4) THEN
              y12 = z4
            ELSE
              y12 = q(i, npy-2)
            END IF
            IF (xt .LT. y12) THEN
              xt = y12
            ELSE
              xt = xt
            END IF
            IF (q(i, npy) .LT. q(i, npy+1)) THEN
              z5 = q(i, npy+1)
            ELSE
              z5 = q(i, npy)
            END IF
            IF (q(i, npy-2) .LT. q(i, npy-1)) THEN
              IF (q(i, npy-1) .LT. z5) THEN
                y13 = z5
              ELSE
                y13 = q(i, npy-1)
              END IF
            ELSE IF (q(i, npy-2) .LT. z5) THEN
              y13 = z5
            ELSE
              y13 = q(i, npy-2)
            END IF
            IF (xt .GT. y13) THEN
              xt = y13
            ELSE
              xt = xt
            END IF
!        endif
            br(i, npy-1) = xt - q(i, npy-1)
            bl(i, npy) = xt - q(i, npy)
            br(i, npy) = s11*(q(i, npy+1)-q(i, npy)) - s14*dm(i, npy+1)
          END DO
          arg1 = 3*(ilast-ifirst+1)
          CALL PERT_PPM(arg1, q(ifirst:ilast, npy-2:npy), bl(ifirst:&
&                 ilast, npy-2:npy), br(ifirst:ilast, npy-2:npy), 1)
        END IF
      END IF
      DO j=js,je+1
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
            flux(i, j) = q(i, j-1) + (1.-c(i, j))*(br(i, j-1)-c(i, j)*(&
&             bl(i, j-1)+br(i, j-1)))
          ELSE
            flux(i, j) = q(i, j) + (1.+c(i, j))*(bl(i, j)+c(i, j)*(bl(i&
&             , j)+br(i, j)))
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE YPPM
  SUBROUTINE MP_GHOST_EW(im, jm, km, nq, ifirst, ilast, jfirst, jlast, &
&   kfirst, klast, ng_w, ng_e, ng_s, ng_n, q_ghst, q)
    IMPLICIT NONE
!
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: im, jm, km, nq
    INTEGER, INTENT(IN) :: ifirst, ilast
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: kfirst, klast
! eastern  zones to ghost
    INTEGER, INTENT(IN) :: ng_e
! western  zones to ghost
    INTEGER, INTENT(IN) :: ng_w
! southern zones to ghost
    INTEGER, INTENT(IN) :: ng_s
! northern zones to ghost
    INTEGER, INTENT(IN) :: ng_n
    REAL, INTENT(INOUT) :: q_ghst(ifirst-ng_w:ilast+ng_e, jfirst-ng_s:&
&   jlast+ng_n, kfirst:klast, nq)
    REAL, OPTIONAL, INTENT(IN) :: q(ifirst:ilast, jfirst:jlast, kfirst:&
&   klast, nq)
!
! !DESCRIPTION:
!
!     Ghost 4d east/west
!
! !REVISION HISTORY:
!    2005.08.22   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: i, j, k, n
    INTRINSIC PRESENT
    IF (PRESENT(q)) q_ghst(ifirst:ilast, jfirst:jlast, kfirst:klast, 1:&
&     nq) = q(ifirst:ilast, jfirst:jlast, kfirst:klast, 1:nq)
!      Assume Periodicity in X-dir and not overlapping
    DO n=1,nq
      DO k=kfirst,klast
        DO j=jfirst-ng_s,jlast+ng_n
          DO i=1,ng_w
            q_ghst(ifirst-i, j, k, n) = q_ghst(ilast-i+1, j, k, n)
          END DO
          DO i=1,ng_e
            q_ghst(ilast+i, j, k, n) = q_ghst(ifirst+i-1, j, k, n)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE MP_GHOST_EW
!  Differentiation of pert_ppm in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_c
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
!   gradient     of useful results: al ar
!   with respect to varying inputs: al ar
  SUBROUTINE PERT_PPM_ADM(im, a0, al, al_ad, ar, ar_ad, iv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    INTEGER, INTENT(IN) :: iv
    REAL, INTENT(IN) :: a0(im)
    REAL, INTENT(INOUT) :: al(im), ar(im)
    REAL, INTENT(INOUT) :: al_ad(im), ar_ad(im)
! Local:
    REAL :: a4, da1, da2, a6da, fmin
    INTEGER :: i
    REAL, PARAMETER :: r12=1./12.
    INTRINSIC ABS
    REAL :: abs0
    INTEGER :: branch
!-----------------------------------
! Optimized PPM in perturbation form:
!-----------------------------------
    IF (iv .EQ. 0) THEN
! Positive definite constraint
      DO i=1,im
        IF (a0(i) .LE. 0.) THEN
          CALL PUSHCONTROL3B(5)
        ELSE
          a4 = -(3.*(ar(i)+al(i)))
          da1 = ar(i) - al(i)
          IF (da1 .GE. 0.) THEN
            abs0 = da1
          ELSE
            abs0 = -da1
          END IF
          IF (abs0 .LT. -a4) THEN
            fmin = a0(i) + 0.25/a4*da1**2 + a4*r12
            IF (fmin .LT. 0.) THEN
              IF (ar(i) .GT. 0. .AND. al(i) .GT. 0.) THEN
                CALL PUSHCONTROL3B(4)
              ELSE IF (da1 .GT. 0.) THEN
                CALL PUSHCONTROL3B(3)
              ELSE
                CALL PUSHCONTROL3B(2)
              END IF
            ELSE
              CALL PUSHCONTROL3B(1)
            END IF
          ELSE
            CALL PUSHCONTROL3B(0)
          END IF
        END IF
      END DO
      DO i=im,1,-1
        CALL POPCONTROL3B(branch)
        IF (branch .LT. 3) THEN
          IF (branch .NE. 0) THEN
            IF (branch .NE. 1) THEN
              ar_ad(i) = ar_ad(i) - 2.*al_ad(i)
              al_ad(i) = 0.0
            END IF
          END IF
        ELSE IF (branch .EQ. 3) THEN
          al_ad(i) = al_ad(i) - 2.*ar_ad(i)
          ar_ad(i) = 0.0
        ELSE IF (branch .EQ. 4) THEN
          al_ad(i) = 0.0
          ar_ad(i) = 0.0
        ELSE
          ar_ad(i) = 0.0
          al_ad(i) = 0.0
        END IF
      END DO
    ELSE
! Standard PPM constraint
      DO i=1,im
        IF (al(i)*ar(i) .LT. 0.) THEN
          da1 = al(i) - ar(i)
          da2 = da1**2
          a6da = 3.*(al(i)+ar(i))*da1
! abs(a6da) > da2 --> 3.*abs(al+ar) > abs(al-ar)
          IF (a6da .LT. -da2) THEN
            CALL PUSHCONTROL2B(3)
          ELSE IF (a6da .GT. da2) THEN
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(0)
        END IF
      END DO
      DO i=im,1,-1
        CALL POPCONTROL2B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            ar_ad(i) = 0.0
            al_ad(i) = 0.0
          END IF
        ELSE IF (branch .EQ. 2) THEN
          ar_ad(i) = ar_ad(i) - 2.*al_ad(i)
          al_ad(i) = 0.0
        ELSE
          al_ad(i) = al_ad(i) - 2.*ar_ad(i)
          ar_ad(i) = 0.0
        END IF
      END DO
    END IF
  END SUBROUTINE PERT_PPM_ADM
  SUBROUTINE PERT_PPM(im, a0, al, ar, iv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    INTEGER, INTENT(IN) :: iv
    REAL, INTENT(IN) :: a0(im)
    REAL, INTENT(INOUT) :: al(im), ar(im)
! Local:
    REAL :: a4, da1, da2, a6da, fmin
    INTEGER :: i
    REAL, PARAMETER :: r12=1./12.
    INTRINSIC ABS
    REAL :: abs0
!-----------------------------------
! Optimized PPM in perturbation form:
!-----------------------------------
    IF (iv .EQ. 0) THEN
! Positive definite constraint
      DO i=1,im
        IF (a0(i) .LE. 0.) THEN
          al(i) = 0.
          ar(i) = 0.
        ELSE
          a4 = -(3.*(ar(i)+al(i)))
          da1 = ar(i) - al(i)
          IF (da1 .GE. 0.) THEN
            abs0 = da1
          ELSE
            abs0 = -da1
          END IF
          IF (abs0 .LT. -a4) THEN
            fmin = a0(i) + 0.25/a4*da1**2 + a4*r12
            IF (fmin .LT. 0.) THEN
              IF (ar(i) .GT. 0. .AND. al(i) .GT. 0.) THEN
                ar(i) = 0.
                al(i) = 0.
              ELSE IF (da1 .GT. 0.) THEN
                ar(i) = -(2.*al(i))
              ELSE
                al(i) = -(2.*ar(i))
              END IF
            END IF
          END IF
        END IF
      END DO
    ELSE
! Standard PPM constraint
      DO i=1,im
        IF (al(i)*ar(i) .LT. 0.) THEN
          da1 = al(i) - ar(i)
          da2 = da1**2
          a6da = 3.*(al(i)+ar(i))*da1
! abs(a6da) > da2 --> 3.*abs(al+ar) > abs(al-ar)
          IF (a6da .LT. -da2) THEN
            ar(i) = -(2.*al(i))
          ELSE IF (a6da .GT. da2) THEN
            al(i) = -(2.*ar(i))
          END IF
        ELSE
! effect of dm=0 included here
          al(i) = 0.
          ar(i) = 0.
        END IF
      END DO
    END IF
  END SUBROUTINE PERT_PPM
!  Differentiation of deln_flux in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_
!core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn
!_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_m
!od.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynami
!cs_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_
!mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mo
!d.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scala
!r_profile_fb fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.stee
!pz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_s
!etup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.com
!pute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_
!c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver 
!nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw
!_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core
!_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.c
!ompute_divergence_damping_fb sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corners
!_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circl
!e_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: q mass fx fy
!   with respect to varying inputs: q mass fx fy
  SUBROUTINE DELN_FLUX_ADM(nord, is, ie, js, je, npx, npy, damp, q, q_ad&
&   , fx, fx_ad, fy, fy_ad, gridstruct, bd, mass, mass_ad)
    IMPLICIT NONE
! Del-n damping for the cell-mean values (A grid)
!------------------
! nord = 0:   del-2
! nord = 1:   del-4
! nord = 2:   del-6
! nord = 3:   del-8 --> requires more ghosting than current
!------------------
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! del-n
    INTEGER, INTENT(IN) :: nord
    INTEGER, INTENT(IN) :: is, ie, js, je, npx, npy
    REAL, INTENT(IN) :: damp
! q ghosted on input
    REAL, INTENT(IN) :: q(bd%is-ng:bd%ie+ng, bd%js-ng:bd%je+ng)
    REAL :: q_ad(bd%is-ng:bd%ie+ng, bd%js-ng:bd%je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! q ghosted on input
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL :: mass_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
! diffusive fluxes:
    REAL, INTENT(INOUT) :: fx(bd%is:bd%ie+1, bd%js:bd%je), fy(bd%is:bd%&
&   ie, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: fx_ad(bd%is:bd%ie+1, bd%js:bd%je), fy_ad(bd%&
&   is:bd%ie, bd%js:bd%je+1)
! local:
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2(bd%isd:bd%ied, bd%&
&   jsd:bd%jed+1)
    REAL :: fx2_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2_ad(bd%isd:bd%ied&
&   , bd%jsd:bd%jed+1)
    REAL :: d2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: d2_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: damp2
    INTEGER :: i, j, n, nt, i1, i2, j1, j2
    INTRINSIC PRESENT
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
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
    i1 = is - 1 - nord
    i2 = ie + 1 + nord
    j1 = js - 1 - nord
    j2 = je + 1 + nord
    IF (.NOT.PRESENT(mass)) THEN
      DO j=j1,j2
        DO i=i1,i2
          d2(i, j) = damp*q(i, j)
        END DO
      END DO
      CALL PUSHCONTROL1B(0)
    ELSE
      DO j=j1,j2
        DO i=i1,i2
          d2(i, j) = q(i, j)
        END DO
      END DO
      CALL PUSHCONTROL1B(1)
    END IF
    IF (nord .GT. 0) THEN
      CALL COPY_CORNERS(d2, npx, npy, 1, gridstruct%nested, bd, &
&                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct&
&                 %nw_corner, gridstruct%ne_corner)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    DO j=js-nord,je+nord
      DO i=is-nord,ie+nord+1
        fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i-1, j)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) THEN
      CALL COPY_CORNERS(d2, npx, npy, 2, gridstruct%nested, bd, &
&                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct&
&                 %nw_corner, gridstruct%ne_corner)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    DO j=js-nord,je+nord+1
      DO i=is-nord,ie+nord
        fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j-1)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) THEN
!----------
! high-order
!----------
      DO n=1,nord
        nt = nord - n
        ad_from0 = js - nt - 1
        DO j=ad_from0,je+nt+1
          ad_from = is - nt - 1
          DO i=ad_from,ie+nt+1
            d2(i, j) = (fx2(i, j)-fx2(i+1, j)+fy2(i, j)-fy2(i, j+1))*&
&             gridstruct%rarea(i, j)
          END DO
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from0)
        CALL COPY_CORNERS(d2, npx, npy, 1, gridstruct%nested, bd, &
&                   gridstruct%sw_corner, gridstruct%se_corner, &
&                   gridstruct%nw_corner, gridstruct%ne_corner)
        ad_from2 = js - nt
        DO j=ad_from2,je+nt
          ad_from1 = is - nt
          DO i=ad_from1,ie+nt+1
            fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i, j)-d2(i-1, j))
          END DO
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from1)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from2)
        CALL COPY_CORNERS(d2, npx, npy, 2, gridstruct%nested, bd, &
&                   gridstruct%sw_corner, gridstruct%se_corner, &
&                   gridstruct%nw_corner, gridstruct%ne_corner)
        ad_from4 = js - nt
        DO j=ad_from4,je+nt+1
          ad_from3 = is - nt
          DO i=ad_from3,ie+nt
            fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j)-d2(i, j-1))
          END DO
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from3)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from4)
      END DO
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
!---------------------------------------------
! Add the diffusive fluxes to the flux arrays:
!---------------------------------------------
    IF (PRESENT(mass)) THEN
! Apply mass weighting to diffusive fluxes:
      damp2 = 0.5*damp
      fy2_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie,is,-1
          temp_ad5 = damp2*fy2(i, j)*fy_ad(i, j)
          mass_ad(i, j-1) = mass_ad(i, j-1) + temp_ad5
          mass_ad(i, j) = mass_ad(i, j) + temp_ad5
          fy2_ad(i, j) = fy2_ad(i, j) + damp2*(mass(i, j-1)+mass(i, j))*&
&           fy_ad(i, j)
        END DO
      END DO
      fx2_ad = 0.0
      DO j=je,js,-1
        DO i=ie+1,is,-1
          temp_ad4 = damp2*fx2(i, j)*fx_ad(i, j)
          mass_ad(i-1, j) = mass_ad(i-1, j) + temp_ad4
          mass_ad(i, j) = mass_ad(i, j) + temp_ad4
          fx2_ad(i, j) = fx2_ad(i, j) + damp2*(mass(i-1, j)+mass(i, j))*&
&           fx_ad(i, j)
        END DO
      END DO
    ELSE
      fy2_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie,is,-1
          fy2_ad(i, j) = fy2_ad(i, j) + fy_ad(i, j)
        END DO
      END DO
      fx2_ad = 0.0
      DO j=je,js,-1
        DO i=ie+1,is,-1
          fx2_ad(i, j) = fx2_ad(i, j) + fx_ad(i, j)
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      d2_ad = 0.0
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
        CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 2, gridstruct%nested&
&                       , bd, gridstruct%sw_corner, gridstruct%se_corner&
&                       , gridstruct%nw_corner, gridstruct%ne_corner)
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
        CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 1, gridstruct%nested&
&                       , bd, gridstruct%sw_corner, gridstruct%se_corner&
&                       , gridstruct%nw_corner, gridstruct%ne_corner)
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
    ELSE
      d2_ad = 0.0
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
&                                      gridstruct%nested, bd, gridstruct&
&                                      %sw_corner, gridstruct%se_corner&
&                                      , gridstruct%nw_corner, &
&                                      gridstruct%ne_corner)
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
&                                      gridstruct%nested, bd, gridstruct&
&                                      %sw_corner, gridstruct%se_corner&
&                                      , gridstruct%nw_corner, &
&                                      gridstruct%ne_corner)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=j2,j1,-1
        DO i=i2,i1,-1
          q_ad(i, j) = q_ad(i, j) + damp*d2_ad(i, j)
          d2_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      DO j=j2,j1,-1
        DO i=i2,i1,-1
          q_ad(i, j) = q_ad(i, j) + d2_ad(i, j)
          d2_ad(i, j) = 0.0
        END DO
      END DO
    END IF
  END SUBROUTINE DELN_FLUX_ADM
!  Differentiation of copy_corners in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 d
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
!   gradient     of useful results: q
!   with respect to varying inputs: q
!Weird arguments are because this routine is called in a lot of
!places outside of tp_core, sometimes very deeply nested in the call tree.
  SUBROUTINE COPY_CORNERS_ADM(q, q_ad, npx, npy, dir, nested, bd, &
&   sw_corner, se_corner, nw_corner, ne_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy, dir
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: nested, sw_corner, se_corner, nw_corner, &
&   ne_corner
    INTEGER :: i, j
    REAL :: tmp
    REAL :: tmp_ad
    REAL :: tmp0
    REAL :: tmp_ad0
    REAL :: tmp1
    REAL :: tmp_ad1
    REAL :: tmp2
    REAL :: tmp_ad2
    REAL :: tmp3
    REAL :: tmp_ad3
    REAL :: tmp4
    REAL :: tmp_ad4
    REAL :: tmp5
    REAL :: tmp_ad5
    REAL :: tmp6
    REAL :: tmp_ad6
    INTEGER :: branch
    IF (.NOT.nested) THEN
      IF (dir .EQ. 1) THEN
! XDir:
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
          DO j=npy+ng-1,npy,-1
            DO i=0,1-ng,-1
              tmp_ad2 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(npy-j, i-1+npx) = q_ad(npy-j, i-1+npx) + tmp_ad2
            END DO
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=npy+ng-1,npy,-1
            DO i=npx+ng-1,npx,-1
              tmp_ad1 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(j, 2*npx-1-i) = q_ad(j, 2*npx-1-i) + tmp_ad1
            END DO
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=0,1-ng,-1
            DO i=npx+ng-1,npx,-1
              tmp_ad0 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(npy-j, i-npx+1) = q_ad(npy-j, i-npx+1) + tmp_ad0
            END DO
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=0,1-ng,-1
            DO i=0,1-ng,-1
              tmp_ad = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(j, 1-i) = q_ad(j, 1-i) + tmp_ad
            END DO
          END DO
        END IF
      ELSE IF (dir .EQ. 2) THEN
! YDir:
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
          DO j=npy+ng-1,npy,-1
            DO i=0,1-ng,-1
              tmp_ad6 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(j+1-npx, npy-i) = q_ad(j+1-npx, npy-i) + tmp_ad6
            END DO
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=npy+ng-1,npy,-1
            DO i=npx+ng-1,npx,-1
              tmp_ad5 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(2*npy-1-j, i) = q_ad(2*npy-1-j, i) + tmp_ad5
            END DO
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=0,1-ng,-1
            DO i=npx+ng-1,npx,-1
              tmp_ad4 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(npy+j-1, npx-i) = q_ad(npy+j-1, npx-i) + tmp_ad4
            END DO
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=0,1-ng,-1
            DO i=0,1-ng,-1
              tmp_ad3 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(1-j, i) = q_ad(1-j, i) + tmp_ad3
            END DO
          END DO
        END IF
      END IF
    END IF
  END SUBROUTINE COPY_CORNERS_ADM
  SUBROUTINE DELN_FLUX(nord, is, ie, js, je, npx, npy, damp, q, fx, fy, &
&   gridstruct, bd, mass)
    IMPLICIT NONE
! Del-n damping for the cell-mean values (A grid)
!------------------
! nord = 0:   del-2
! nord = 1:   del-4
! nord = 2:   del-6
! nord = 3:   del-8 --> requires more ghosting than current
!------------------
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! del-n
    INTEGER, INTENT(IN) :: nord
    INTEGER, INTENT(IN) :: is, ie, js, je, npx, npy
    REAL, INTENT(IN) :: damp
! q ghosted on input
    REAL, INTENT(IN) :: q(bd%is-ng:bd%ie+ng, bd%js-ng:bd%je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! q ghosted on input
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
! diffusive fluxes:
    REAL, INTENT(INOUT) :: fx(bd%is:bd%ie+1, bd%js:bd%je), fy(bd%is:bd%&
&   ie, bd%js:bd%je+1)
! local:
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2(bd%isd:bd%ied, bd%&
&   jsd:bd%jed+1)
    REAL :: d2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: damp2
    INTEGER :: i, j, n, nt, i1, i2, j1, j2
    INTRINSIC PRESENT
    i1 = is - 1 - nord
    i2 = ie + 1 + nord
    j1 = js - 1 - nord
    j2 = je + 1 + nord
    IF (.NOT.PRESENT(mass)) THEN
      DO j=j1,j2
        DO i=i1,i2
          d2(i, j) = damp*q(i, j)
        END DO
      END DO
    ELSE
      DO j=j1,j2
        DO i=i1,i2
          d2(i, j) = q(i, j)
        END DO
      END DO
    END IF
    IF (nord .GT. 0) CALL COPY_CORNERS(d2, npx, npy, 1, gridstruct%&
&                                nested, bd, gridstruct%sw_corner, &
&                                gridstruct%se_corner, gridstruct%&
&                                nw_corner, gridstruct%ne_corner)
    DO j=js-nord,je+nord
      DO i=is-nord,ie+nord+1
        fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i-1, j)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) CALL COPY_CORNERS(d2, npx, npy, 2, gridstruct%&
&                                nested, bd, gridstruct%sw_corner, &
&                                gridstruct%se_corner, gridstruct%&
&                                nw_corner, gridstruct%ne_corner)
    DO j=js-nord,je+nord+1
      DO i=is-nord,ie+nord
        fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j-1)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) THEN
!----------
! high-order
!----------
      DO n=1,nord
        nt = nord - n
        DO j=js-nt-1,je+nt+1
          DO i=is-nt-1,ie+nt+1
            d2(i, j) = (fx2(i, j)-fx2(i+1, j)+fy2(i, j)-fy2(i, j+1))*&
&             gridstruct%rarea(i, j)
          END DO
        END DO
        CALL COPY_CORNERS(d2, npx, npy, 1, gridstruct%nested, bd, &
&                   gridstruct%sw_corner, gridstruct%se_corner, &
&                   gridstruct%nw_corner, gridstruct%ne_corner)
        DO j=js-nt,je+nt
          DO i=is-nt,ie+nt+1
            fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i, j)-d2(i-1, j))
          END DO
        END DO
        CALL COPY_CORNERS(d2, npx, npy, 2, gridstruct%nested, bd, &
&                   gridstruct%sw_corner, gridstruct%se_corner, &
&                   gridstruct%nw_corner, gridstruct%ne_corner)
        DO j=js-nt,je+nt+1
          DO i=is-nt,ie+nt
            fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j)-d2(i, j-1))
          END DO
        END DO
      END DO
    END IF
!---------------------------------------------
! Add the diffusive fluxes to the flux arrays:
!---------------------------------------------
    IF (PRESENT(mass)) THEN
! Apply mass weighting to diffusive fluxes:
      damp2 = 0.5*damp
      DO j=js,je
        DO i=is,ie+1
          fx(i, j) = fx(i, j) + damp2*(mass(i-1, j)+mass(i, j))*fx2(i, j&
&           )
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy(i, j) = fy(i, j) + damp2*(mass(i, j-1)+mass(i, j))*fy2(i, j&
&           )
        END DO
      END DO
    ELSE
      DO j=js,je
        DO i=is,ie+1
          fx(i, j) = fx(i, j) + fx2(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy(i, j) = fy(i, j) + fy2(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE DELN_FLUX
!  Differentiation of fv_tp_2d in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_
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
!   gradient     of useful results: xfx q mass mfx mfy ra_x ra_y
!                yfx fx fy crx cry
!   with respect to varying inputs: xfx q mass mfx mfy ra_x ra_y
!                yfx fx fy crx cry
  SUBROUTINE FV_TP_2D_FWD(q, crx, cry, npx, npy, hord, fx, fy, xfx, &
&   yfx, gridstruct, bd, ra_x, ra_y, mfx, mfy, mass, nord, damp_c)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL, INTENT(IN) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed)
!
    REAL, INTENT(IN) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed)
!
    REAL, INTENT(IN) :: cry(bd%isd:bd%ied, bd%js:bd%je+1)
!
    REAL, INTENT(IN) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(IN) :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
! transported scalar
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
! Flux in x ( E )
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
! Flux in y ( N )
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! optional Arguments:
! Mass Flux X-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfx(bd%is:bd%ie+1, bd%js:bd%je)
! Mass Flux Y-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL, INTENT(IN) :: damp_c
    INTEGER, OPTIONAL, INTENT(IN) :: nord
! Local:
    INTEGER :: ord_ou, ord_in
    REAL :: q_i(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: q_j(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: fx2(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: fy2(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fyy(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fx1(bd%is:bd%ie+1)
    REAL :: damp
    INTEGER :: i, j
    INTEGER :: is, ie, js, je, isd, ied, jsd, jed
    INTRINSIC PRESENT

    ord_ou = 0
    ord_in = 0
    q_i = 0.0
    q_j = 0.0
    fx2 = 0.0
    fy2 = 0.0
    fyy = 0.0
    fx1 = 0.0
    damp = 0.0
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
    IF (hord .EQ. 10) THEN
      ord_in = 8
    ELSE
      ord_in = hord
    END IF
    ord_ou = hord
    IF (.NOT.gridstruct%nested) THEN
      CALL COPY_CORNERS_FWD(q, npx, npy, 2, gridstruct%nested, bd, &
&                        gridstruct%sw_corner, gridstruct%se_corner, &
&                        gridstruct%nw_corner, gridstruct%ne_corner)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    CALL YPPM_FWD(fy2, q, cry, ord_in, isd, ied, isd, ied, js, je, &
&              jsd, jed, npx, npy, gridstruct%dya, gridstruct%nested, &
&              gridstruct%grid_type)
    DO j=js,je+1
      DO i=isd,ied
        fyy(i, j) = yfx(i, j)*fy2(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        q_i(i, j) = (q(i, j)*gridstruct%area(i, j)+fyy(i, j)-fyy(i, j+1)&
&         )/ra_y(i, j)
      END DO
    END DO
    CALL XPPM_FWD(fx, q_i, crx(is:ie+1, js:je), ord_ou, is, ie, isd, &
&              ied, js, je, jsd, jed, npx, npy, gridstruct%dxa, &
&              gridstruct%nested, gridstruct%grid_type)
    IF (.NOT.gridstruct%nested) THEN
      CALL COPY_CORNERS_FWD(q, npx, npy, 1, gridstruct%nested, bd, &
&                        gridstruct%sw_corner, gridstruct%se_corner, &
&                        gridstruct%nw_corner, gridstruct%ne_corner)
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
    CALL XPPM_FWD(fx2, q, crx, ord_in, is, ie, isd, ied, jsd, jed, &
&              jsd, jed, npx, npy, gridstruct%dxa, gridstruct%nested, &
&              gridstruct%grid_type)
    DO j=jsd,jed
      DO i=is,ie+1
        CALL PUSHREALARRAY(fx1(i))
        fx1(i) = xfx(i, j)*fx2(i, j)
      END DO
      DO i=is,ie
        q_j(i, j) = (q(i, j)*gridstruct%area(i, j)+fx1(i)-fx1(i+1))/ra_x&
&         (i, j)
      END DO
    END DO
    CALL YPPM_FWD(fy, q_j, cry, ord_ou, is, ie, isd, ied, js, je, jsd&
&              , jed, npx, npy, gridstruct%dya, gridstruct%nested, &
&              gridstruct%grid_type)
!----------------
! Flux averaging:
!----------------
    IF (PRESENT(mfx) .AND. PRESENT(mfy)) THEN
!---------------------------------
! For transport of pt and tracers
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(fx(i, j))
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*mfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(fy(i, j))
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*mfy(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c) .AND. PRESENT(mass)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*gridstruct%da_min)**(nord+1)
          CALL DELN_FLUX_FWD(nord, is, ie, js, je, npx, npy, damp, q&
&                         , fx, fy, gridstruct, bd, mass)
          CALL PUSHREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
          CALL PUSHREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
          CALL PUSHREALARRAY(fx1, bd%ie - bd%is + 2)
          CALL PUSHREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
          CALL PUSHREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
          CALL PUSHCONTROL(3,2)
        ELSE
          CALL PUSHREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
          CALL PUSHREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
          CALL PUSHREALARRAY(fx1, bd%ie - bd%is + 2)
          CALL PUSHINTEGER(je)
          CALL PUSHINTEGER(is)
          CALL PUSHINTEGER(ie)
          CALL PUSHREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
          CALL PUSHREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
          CALL PUSHINTEGER(js)
          CALL PUSHCONTROL(3,1)
        END IF
      ELSE
        CALL PUSHREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL PUSHREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
        CALL PUSHREALARRAY(fx1, bd%ie - bd%is + 2)
        CALL PUSHINTEGER(je)
        CALL PUSHINTEGER(is)
        CALL PUSHINTEGER(ie)
        CALL PUSHREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
        CALL PUSHREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL PUSHINTEGER(js)
        CALL PUSHCONTROL(3,0)
      END IF
    ELSE
!---------------------------------
! For transport of delp, vorticity
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(fx(i, j))
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*xfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(fy(i, j))
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*yfx(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*gridstruct%da_min)**(nord+1)
          CALL DELN_FLUX_FWD(nord, is, ie, js, je, npx, npy, damp, q&
&                         , fx, fy, gridstruct, bd)
          CALL PUSHREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
          CALL PUSHREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
          CALL PUSHREALARRAY(fx1, bd%ie - bd%is + 2)
          CALL PUSHREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
          CALL PUSHREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
          CALL PUSHCONTROL(3,5)
        ELSE
          CALL PUSHREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
          CALL PUSHREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
          CALL PUSHREALARRAY(fx1, bd%ie - bd%is + 2)
          CALL PUSHINTEGER(je)
          CALL PUSHINTEGER(is)
          CALL PUSHINTEGER(ie)
          CALL PUSHREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
          CALL PUSHREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
          CALL PUSHINTEGER(js)
          CALL PUSHCONTROL(3,4)
        END IF
      ELSE
        CALL PUSHREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL PUSHREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
        CALL PUSHREALARRAY(fx1, bd%ie - bd%is + 2)
        CALL PUSHINTEGER(je)
        CALL PUSHINTEGER(is)
        CALL PUSHINTEGER(ie)
        CALL PUSHREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
        CALL PUSHREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL PUSHINTEGER(js)
        CALL PUSHCONTROL(3,3)
      END IF
    END IF
  END SUBROUTINE FV_TP_2D_FWD
!  Differentiation of fv_tp_2d in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: xfx q mass mfx mfy ra_x ra_y
!                yfx fx fy crx cry
!   with respect to varying inputs: xfx q mass mfx mfy ra_x ra_y
!                yfx fx fy crx cry
  SUBROUTINE FV_TP_2D_BWD(q, q_ad, crx, crx_ad, cry, cry_ad, npx, npy&
&   , hord, fx, fx_ad, fy, fy_ad, xfx, xfx_ad, yfx, yfx_ad, gridstruct, &
&   bd, ra_x, ra_x_ad, ra_y, ra_y_ad, mfx, mfx_ad, mfy, mfy_ad, mass, &
&   mass_ad, nord, damp_c)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
    REAL, INTENT(IN) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: crx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: xfx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: cry(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: cry_ad(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(IN) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: yfx_ad(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(IN) :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_ad(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_ad(bd%isd:bd%ied, bd%js:bd%je)
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_ad(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_ad(bd%is:bd%ie, bd%js:bd%je+1)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    REAL, OPTIONAL, INTENT(IN) :: mfx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, OPTIONAL :: mfx_ad(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, OPTIONAL, INTENT(IN) :: mfy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, OPTIONAL :: mfy_ad(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL :: mass_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL, INTENT(IN) :: damp_c
    INTEGER, OPTIONAL, INTENT(IN) :: nord
    INTEGER :: ord_ou, ord_in
    REAL :: q_i(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: q_i_ad(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: q_j(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: q_j_ad(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: fx2(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: fx2_ad(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: fy2(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fy2_ad(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fyy(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fyy_ad(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fx1(bd%is:bd%ie+1)
    REAL :: fx1_ad(bd%is:bd%ie+1)
    REAL :: damp
    INTEGER :: i, j
    INTEGER :: is, ie, js, je, isd, ied, jsd, jed
    INTRINSIC PRESENT
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    INTEGER :: branch

    ord_ou = 0
    ord_in = 0
    q_i = 0.0
    q_j = 0.0
    fx2 = 0.0
    fy2 = 0.0
    fyy = 0.0
    fx1 = 0.0
    damp = 0.0
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
    IF (branch .LT. 3) THEN
      IF (branch .EQ. 0) THEN
        CALL POPINTEGER(js)
        CALL POPREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL POPREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
        CALL POPINTEGER(ie)
        CALL POPINTEGER(is)
        CALL POPINTEGER(je)
        CALL POPREALARRAY(fx1, bd%ie - bd%is + 2)
        CALL POPREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
        CALL POPREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
      ELSE IF (branch .EQ. 1) THEN
        CALL POPINTEGER(js)
        CALL POPREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL POPREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
        CALL POPINTEGER(ie)
        CALL POPINTEGER(is)
        CALL POPINTEGER(je)
        CALL POPREALARRAY(fx1, bd%ie - bd%is + 2)
        CALL POPREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
        CALL POPREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
      ELSE
        CALL POPREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL POPREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
        CALL POPREALARRAY(fx1, bd%ie - bd%is + 2)
        CALL POPREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
        CALL POPREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        js = bd%js
        damp = (damp_c*gridstruct%da_min)**(nord+1)
        ie = bd%ie
        is = bd%is
        je = bd%je
        CALL DELN_FLUX_BWD(nord, is, ie, js, je, npx, npy, damp, q, &
&                       q_ad, fx, fx_ad, fy, fy_ad, gridstruct, bd, mass&
&                       , mass_ad)
      END IF
      fy2_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(fy(i, j))
          temp_ad2 = 0.5*mfy(i, j)*fy_ad(i, j)
          fy2_ad(i, j) = fy2_ad(i, j) + temp_ad2
          mfy_ad(i, j) = mfy_ad(i, j) + 0.5*(fy(i, j)+fy2(i, j))*fy_ad(i&
&           , j)
          fy_ad(i, j) = temp_ad2
        END DO
      END DO
      fx2_ad = 0.0
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(fx(i, j))
          temp_ad1 = 0.5*mfx(i, j)*fx_ad(i, j)
          fx2_ad(i, j) = fx2_ad(i, j) + temp_ad1
          mfx_ad(i, j) = mfx_ad(i, j) + 0.5*(fx(i, j)+fx2(i, j))*fx_ad(i&
&           , j)
          fx_ad(i, j) = temp_ad1
        END DO
      END DO
    ELSE
      IF (branch .EQ. 3) THEN
        CALL POPINTEGER(js)
        CALL POPREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL POPREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
        CALL POPINTEGER(ie)
        CALL POPINTEGER(is)
        CALL POPINTEGER(je)
        CALL POPREALARRAY(fx1, bd%ie - bd%is + 2)
        CALL POPREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
        CALL POPREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
      ELSE IF (branch .EQ. 4) THEN
        CALL POPINTEGER(js)
        CALL POPREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL POPREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
        CALL POPINTEGER(ie)
        CALL POPINTEGER(is)
        CALL POPINTEGER(je)
        CALL POPREALARRAY(fx1, bd%ie - bd%is + 2)
        CALL POPREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
        CALL POPREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
      ELSE
        CALL POPREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        CALL POPREALARRAY(q_j, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
        CALL POPREALARRAY(fx1, bd%ie - bd%is + 2)
        CALL POPREALARRAY(fx2, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1))
        CALL POPREALARRAY(fyy, (bd%ied-bd%isd+1)*(bd%je-bd%js+2))
        js = bd%js
        damp = (damp_c*gridstruct%da_min)**(nord+1)
        ie = bd%ie
        is = bd%is
        je = bd%je
        CALL DELN_FLUX_BWD(nord, is, ie, js, je, npx, npy, damp, q, &
&                       q_ad, fx, fx_ad, fy, fy_ad, gridstruct, bd)
      END IF
      fy2_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(fy(i, j))
          temp_ad4 = 0.5*yfx(i, j)*fy_ad(i, j)
          fy2_ad(i, j) = fy2_ad(i, j) + temp_ad4
          yfx_ad(i, j) = yfx_ad(i, j) + 0.5*(fy(i, j)+fy2(i, j))*fy_ad(i&
&           , j)
          fy_ad(i, j) = temp_ad4
        END DO
      END DO
      fx2_ad = 0.0
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(fx(i, j))
          temp_ad3 = 0.5*xfx(i, j)*fx_ad(i, j)
          fx2_ad(i, j) = fx2_ad(i, j) + temp_ad3
          xfx_ad(i, j) = xfx_ad(i, j) + 0.5*(fx(i, j)+fx2(i, j))*fx_ad(i&
&           , j)
          fx_ad(i, j) = temp_ad3
        END DO
      END DO
    END IF
    jsd = bd%jsd
    ied = bd%ied
    isd = bd%isd
    jed = bd%jed
    q_j_ad = 0.0
    CALL YPPM_BWD(fy, fy_ad, q_j, q_j_ad, cry, cry_ad, ord_ou, is, ie&
&              , isd, ied, js, je, jsd, jed, npx, npy, gridstruct%dya, &
&              gridstruct%nested, gridstruct%grid_type)
    fx1_ad = 0.0
    DO j=jed,jsd,-1
      DO i=ie,is,-1
        temp_ad0 = q_j_ad(i, j)/ra_x(i, j)
        q_ad(i, j) = q_ad(i, j) + gridstruct%area(i, j)*temp_ad0
        fx1_ad(i) = fx1_ad(i) + temp_ad0
        fx1_ad(i+1) = fx1_ad(i+1) - temp_ad0
        ra_x_ad(i, j) = ra_x_ad(i, j) - (gridstruct%area(i, j)*q(i, j)+&
&         fx1(i)-fx1(i+1))*temp_ad0/ra_x(i, j)
        q_j_ad(i, j) = 0.0
      END DO
      DO i=ie+1,is,-1
        CALL POPREALARRAY(fx1(i))
        xfx_ad(i, j) = xfx_ad(i, j) + fx2(i, j)*fx1_ad(i)
        fx2_ad(i, j) = fx2_ad(i, j) + xfx(i, j)*fx1_ad(i)
        fx1_ad(i) = 0.0
      END DO
    END DO
    CALL XPPM_BWD(fx2, fx2_ad, q, q_ad, crx, crx_ad, ord_in, is, ie, &
&              isd, ied, jsd, jed, jsd, jed, npx, npy, gridstruct%dxa, &
&              gridstruct%nested, gridstruct%grid_type)
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) CALL COPY_CORNERS_BWD(q, q_ad, npx, npy, 1, &
&                                         gridstruct%nested, bd, &
&                                         gridstruct%sw_corner, &
&                                         gridstruct%se_corner, &
&                                         gridstruct%nw_corner, &
&                                         gridstruct%ne_corner)
    q_i_ad = 0.0
    CALL XPPM_BWD(fx, fx_ad, q_i, q_i_ad, crx(is:ie+1, js:je), crx_ad&
&              (is:ie+1, js:je), ord_ou, is, ie, isd, ied, js, je, jsd, &
&              jed, npx, npy, gridstruct%dxa, gridstruct%nested, &
&              gridstruct%grid_type)
    fyy_ad = 0.0
    DO j=je,js,-1
      DO i=ied,isd,-1
        temp_ad = q_i_ad(i, j)/ra_y(i, j)
        q_ad(i, j) = q_ad(i, j) + gridstruct%area(i, j)*temp_ad
        fyy_ad(i, j) = fyy_ad(i, j) + temp_ad
        fyy_ad(i, j+1) = fyy_ad(i, j+1) - temp_ad
        ra_y_ad(i, j) = ra_y_ad(i, j) - (gridstruct%area(i, j)*q(i, j)+&
&         fyy(i, j)-fyy(i, j+1))*temp_ad/ra_y(i, j)
        q_i_ad(i, j) = 0.0
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ied,isd,-1
        yfx_ad(i, j) = yfx_ad(i, j) + fy2(i, j)*fyy_ad(i, j)
        fy2_ad(i, j) = fy2_ad(i, j) + yfx(i, j)*fyy_ad(i, j)
        fyy_ad(i, j) = 0.0
      END DO
    END DO
    CALL YPPM_BWD(fy2, fy2_ad, q, q_ad, cry, cry_ad, ord_in, isd, ied&
&              , isd, ied, js, je, jsd, jed, npx, npy, gridstruct%dya, &
&              gridstruct%nested, gridstruct%grid_type)
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) CALL COPY_CORNERS_BWD(q, q_ad, npx, npy, 2, &
&                                         gridstruct%nested, bd, &
&                                         gridstruct%sw_corner, &
&                                         gridstruct%se_corner, &
&                                         gridstruct%nw_corner, &
&                                         gridstruct%ne_corner)
  END SUBROUTINE FV_TP_2D_BWD
!  Differentiation of xppm in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.
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
!   gradient     of useful results: q flux c
!   with respect to varying inputs: q flux c
  SUBROUTINE XPPM_FWD(flux, q, c, iord, is, ie, isd, ied, jfirst, &
&   jlast, jsd, jed, npx, npy, dxa, nested, grid_type)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, isd, ied, jsd, jed
! compute domain
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: iord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(isd:ied, jfirst:jlast)
! Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: dxa(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
! !OUTPUT PARAMETERS:
!  Flux
    REAL :: flux(is:ie+1, jfirst:jlast)
! Local
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    REAL :: q1(isd:ied)
    REAL, DIMENSION(is:ie+1) :: fx0, fx1
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: al(is-1:ie+2)
    REAL :: dm(is-2:ie+2)
    REAL :: dq(is-3:ie+2)
    INTEGER :: i, j, ie3, is1, ie1
    REAL :: x0, x1, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2
    INTRINSIC MAX
    INTRINSIC MIN

    bl = 0.0
    br = 0.0
    b0 = 0.0
    q1 = 0.0
    fx0 = 0.0
    fx1 = 0.0
    al = 0.0
    dm = 0.0
    dq = 0.0
    ie3 = 0
    is1 = 0
    ie1 = 0
    x0 = 0.0
    x1 = 0.0
    xt = 0.0
    qtmp = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0

    IF (.NOT.nested .AND. grid_type .LT. 3) THEN
      IF (3 .LT. is - 1) THEN
        is1 = is - 1
      ELSE
        is1 = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        CALL PUSHCONTROL(1,1)
        ie3 = ie + 2
      ELSE
        CALL PUSHCONTROL(1,1)
        ie3 = npx - 2
      END IF
    ELSE
      CALL PUSHCONTROL(1,0)
      is1 = is - 1
      ie3 = ie + 2
    END IF
    DO j=jfirst,jlast
      DO i=isd,ied
        CALL PUSHREALARRAY(q1(i))
        q1(i) = q(i, j)
      END DO
      IF (iord .LT. 8 .OR. iord .EQ. 333) THEN
! ord = 2: perfectly linear ppm scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
        DO i=is1,ie3
          CALL PUSHREALARRAY(al(i))
          al(i) = p1*(q1(i-1)+q1(i)) + p2*(q1(i-2)+q1(i+1))
        END DO
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            CALL PUSHREALARRAY(al(0))
            al(0) = c1*q1(-2) + c2*q1(-1) + c3*q1(0)
            CALL PUSHREALARRAY(al(1))
            al(1) = 0.5*(((2.*dxa(0, j)+dxa(-1, j))*q1(0)-dxa(0, j)*q1(-&
&             1))/(dxa(-1, j)+dxa(0, j))+((2.*dxa(1, j)+dxa(2, j))*q1(1)&
&             -dxa(1, j)*q1(2))/(dxa(1, j)+dxa(2, j)))
            CALL PUSHREALARRAY(al(2))
            al(2) = c3*q1(1) + c2*q1(2) + c1*q1(3)
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            CALL PUSHREALARRAY(al(npx-1))
            al(npx-1) = c1*q1(npx-3) + c2*q1(npx-2) + c3*q1(npx-1)
            CALL PUSHREALARRAY(al(npx))
            al(npx) = 0.5*(((2.*dxa(npx-1, j)+dxa(npx-2, j))*q1(npx-1)-&
&             dxa(npx-1, j)*q1(npx-2))/(dxa(npx-2, j)+dxa(npx-1, j))+((&
&             2.*dxa(npx, j)+dxa(npx+1, j))*q1(npx)-dxa(npx, j)*q1(npx+1&
&             ))/(dxa(npx, j)+dxa(npx+1, j)))
            CALL PUSHREALARRAY(al(npx+1))
            al(npx+1) = c3*q1(npx) + c2*q1(npx+1) + c1*q1(npx+2)
            CALL PUSHCONTROL(2,0)
          ELSE
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE
          CALL PUSHCONTROL(2,2)
        END IF
        IF (iord .EQ. 1) THEN
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = q1(i-1)
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = q1(i)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          CALL PUSHCONTROL(3,0)
        ELSE IF (iord .EQ. 2) THEN
! perfectly linear scheme
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              CALL PUSHREALARRAY(qtmp)
              qtmp = q1(i-1)
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = qtmp + (1.-xt)*(al(i)-qtmp-xt*(al(i-1)+al(i)-&
&               (qtmp+qtmp)))
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHREALARRAY(qtmp)
              qtmp = q1(i)
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = qtmp + (1.+xt)*(al(i)-qtmp+xt*(al(i)+al(i+1)-&
&               (qtmp+qtmp)))
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          CALL PUSHCONTROL(3,1)
        ELSE IF (iord .EQ. 333) THEN
! Perfectly linear scheme, more diffusive than ord=2 (HoldawayKent-2015-TellusA)
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = (2.0*q1(i)+5.0*q1(i-1)-q1(i-2))/6.0 - 0.5*xt*&
&               (q1(i)-q1(i-1)) + xt*xt/6.0*(q1(i)-2.0*q1(i-1)+q1(i-2))
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = (2.0*q1(i-1)+5.0*q1(i)-q1(i+1))/6.0 - 0.5*xt*&
&               (q1(i)-q1(i-1)) + xt*xt/6.0*(q1(i+1)-2.0*q1(i)+q1(i-1))
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
          CALL PUSHCONTROL(3,2)
        ELSE
          CALL PUSHCONTROL(3,3)
        END IF
      ELSE
        CALL PUSHCONTROL(3,4)
      END IF
    END DO
    CALL PUSHINTEGER(ie3)
    CALL PUSHREALARRAY(q1, ied - isd + 1)
    CALL PUSHREALARRAY(qtmp)
    CALL PUSHINTEGER(is1)
    CALL PUSHREALARRAY(al, ie - is + 4)
  END SUBROUTINE XPPM_FWD
!  Differentiation of xppm in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
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
!   gradient     of useful results: q flux c
!   with respect to varying inputs: q flux c
  SUBROUTINE XPPM_BWD(flux, flux_ad, q, q_ad, c, c_ad, iord, is, ie, &
&   isd, ied, jfirst, jlast, jsd, jed, npx, npy, dxa, nested, grid_type)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: iord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(isd:ied, jfirst:jlast)
    REAL :: q_ad(isd:ied, jfirst:jlast)
    REAL, INTENT(IN) :: c(is:ie+1, jfirst:jlast)
    REAL :: c_ad(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: dxa(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
    REAL :: flux(is:ie+1, jfirst:jlast)
    REAL :: flux_ad(is:ie+1, jfirst:jlast)
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    REAL :: q1(isd:ied)
    REAL :: q1_ad(isd:ied)
    REAL, DIMENSION(is:ie+1) :: fx0, fx1
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: al(is-1:ie+2)
    REAL :: al_ad(is-1:ie+2)
    REAL :: dm(is-2:ie+2)
    REAL :: dq(is-3:ie+2)
    INTEGER :: i, j, ie3, is1, ie1
    REAL :: x0, x1, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2
    REAL :: xt_ad, qtmp_ad
    INTRINSIC MAX
    INTRINSIC MIN
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp0
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    INTEGER :: branch

    bl = 0.0
    br = 0.0
    b0 = 0.0
    q1 = 0.0
    fx0 = 0.0
    fx1 = 0.0
    al = 0.0
    dm = 0.0
    dq = 0.0
    ie3 = 0
    is1 = 0
    ie1 = 0
    x0 = 0.0
    x1 = 0.0
    xt = 0.0
    qtmp = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    branch = 0

    CALL POPREALARRAY(al, ie - is + 4)
    CALL POPINTEGER(is1)
    CALL POPREALARRAY(qtmp)
    CALL POPREALARRAY(q1, ied - isd + 1)
    CALL POPINTEGER(ie3)
    al_ad = 0.0
    q1_ad = 0.0
    DO j=jlast,jfirst,-1
      CALL POPCONTROL(3,branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          DO i=ie+1,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY(flux(i, j))
              q1_ad(i) = q1_ad(i) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              CALL POPREALARRAY(flux(i, j))
              q1_ad(i-1) = q1_ad(i-1) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            END IF
          END DO
        ELSE
          DO i=ie+1,is,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              xt = c(i, j)
              qtmp = q1(i)
              CALL POPREALARRAY(flux(i, j))
              temp0 = al(i) + al(i+1) - 2*qtmp
              temp_ad5 = (xt+1.)*flux_ad(i, j)
              temp_ad6 = xt*temp_ad5
              qtmp_ad = flux_ad(i, j) - temp_ad5 - 2*temp_ad6
              xt_ad = temp0*temp_ad5 + (al(i)-qtmp+xt*temp0)*flux_ad(i, &
&               j)
              al_ad(i) = al_ad(i) + temp_ad6 + temp_ad5
              al_ad(i+1) = al_ad(i+1) + temp_ad6
              flux_ad(i, j) = 0.0
              CALL POPREALARRAY(qtmp)
              q1_ad(i) = q1_ad(i) + qtmp_ad
            ELSE
              xt = c(i, j)
              qtmp = q1(i-1)
              CALL POPREALARRAY(flux(i, j))
              temp = al(i-1) + al(i) - 2*qtmp
              temp_ad3 = (1.-xt)*flux_ad(i, j)
              temp_ad4 = -(xt*temp_ad3)
              qtmp_ad = flux_ad(i, j) - temp_ad3 - 2*temp_ad4
              xt_ad = -(temp*temp_ad3) - (al(i)-qtmp-xt*temp)*flux_ad(i&
&               , j)
              al_ad(i) = al_ad(i) + temp_ad4 + temp_ad3
              al_ad(i-1) = al_ad(i-1) + temp_ad4
              flux_ad(i, j) = 0.0
              CALL POPREALARRAY(qtmp)
              q1_ad(i-1) = q1_ad(i-1) + qtmp_ad
            END IF
            c_ad(i, j) = c_ad(i, j) + xt_ad
          END DO
        END IF
      ELSE IF (branch .EQ. 2) THEN
        DO i=ie+1,is,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            xt = c(i, j)
            CALL POPREALARRAY(flux(i, j))
            temp_ad10 = flux_ad(i, j)/6.0
            temp_ad11 = -(0.5*xt*flux_ad(i, j))
            temp_ad12 = xt**2*flux_ad(i, j)/6.0
            q1_ad(i-1) = q1_ad(i-1) + temp_ad12 - temp_ad11 + 2.0*&
&             temp_ad10
            q1_ad(i) = q1_ad(i) + temp_ad11 - 2.0*temp_ad12 + 5.0*&
&             temp_ad10
            q1_ad(i+1) = q1_ad(i+1) + temp_ad12 - temp_ad10
            xt_ad = ((q1(i+1)-2.0*q1(i)+q1(i-1))*2*xt/6.0-0.5*(q1(i)-q1(&
&             i-1)))*flux_ad(i, j)
            flux_ad(i, j) = 0.0
          ELSE
            xt = c(i, j)
            CALL POPREALARRAY(flux(i, j))
            temp_ad7 = flux_ad(i, j)/6.0
            temp_ad8 = -(0.5*xt*flux_ad(i, j))
            temp_ad9 = xt**2*flux_ad(i, j)/6.0
            q1_ad(i) = q1_ad(i) + temp_ad9 + temp_ad8 + 2.0*temp_ad7
            q1_ad(i-1) = q1_ad(i-1) + 5.0*temp_ad7 - temp_ad8 - 2.0*&
&             temp_ad9
            q1_ad(i-2) = q1_ad(i-2) + temp_ad9 - temp_ad7
            xt_ad = ((q1(i)-2.0*q1(i-1)+q1(i-2))*2*xt/6.0-0.5*(q1(i)-q1(&
&             i-1)))*flux_ad(i, j)
            flux_ad(i, j) = 0.0
          END IF
          c_ad(i, j) = c_ad(i, j) + xt_ad
        END DO
      ELSE IF (branch .NE. 3) THEN
        GOTO 110
      END IF
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(al(npx+1))
        q1_ad(npx) = q1_ad(npx) + c3*al_ad(npx+1)
        q1_ad(npx+1) = q1_ad(npx+1) + c2*al_ad(npx+1)
        q1_ad(npx+2) = q1_ad(npx+2) + c1*al_ad(npx+1)
        al_ad(npx+1) = 0.0
        CALL POPREALARRAY(al(npx))
        temp_ad1 = 0.5*al_ad(npx)/(dxa(npx-2, j)+dxa(npx-1, j))
        temp_ad2 = 0.5*al_ad(npx)/(dxa(npx, j)+dxa(npx+1, j))
        q1_ad(npx-1) = q1_ad(npx-1) + (dxa(npx-1, j)*2.+dxa(npx-2, j))*&
&         temp_ad1
        q1_ad(npx-2) = q1_ad(npx-2) - dxa(npx-1, j)*temp_ad1
        q1_ad(npx) = q1_ad(npx) + (dxa(npx, j)*2.+dxa(npx+1, j))*&
&         temp_ad2
        q1_ad(npx+1) = q1_ad(npx+1) - dxa(npx, j)*temp_ad2
        al_ad(npx) = 0.0
        CALL POPREALARRAY(al(npx-1))
        q1_ad(npx-3) = q1_ad(npx-3) + c1*al_ad(npx-1)
        q1_ad(npx-2) = q1_ad(npx-2) + c2*al_ad(npx-1)
        q1_ad(npx-1) = q1_ad(npx-1) + c3*al_ad(npx-1)
        al_ad(npx-1) = 0.0
      ELSE IF (branch .NE. 1) THEN
        GOTO 100
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        CALL POPREALARRAY(al(2))
        q1_ad(1) = q1_ad(1) + c3*al_ad(2)
        q1_ad(2) = q1_ad(2) + c2*al_ad(2)
        q1_ad(3) = q1_ad(3) + c1*al_ad(2)
        al_ad(2) = 0.0
        CALL POPREALARRAY(al(1))
        temp_ad = 0.5*al_ad(1)/(dxa(-1, j)+dxa(0, j))
        temp_ad0 = 0.5*al_ad(1)/(dxa(1, j)+dxa(2, j))
        q1_ad(0) = q1_ad(0) + (dxa(0, j)*2.+dxa(-1, j))*temp_ad
        q1_ad(-1) = q1_ad(-1) - dxa(0, j)*temp_ad
        q1_ad(1) = q1_ad(1) + (dxa(1, j)*2.+dxa(2, j))*temp_ad0
        q1_ad(2) = q1_ad(2) - dxa(1, j)*temp_ad0
        al_ad(1) = 0.0
        CALL POPREALARRAY(al(0))
        q1_ad(-2) = q1_ad(-2) + c1*al_ad(0)
        q1_ad(-1) = q1_ad(-1) + c2*al_ad(0)
        q1_ad(0) = q1_ad(0) + c3*al_ad(0)
        al_ad(0) = 0.0
      END IF
 100  DO i=ie3,is1,-1
        CALL POPREALARRAY(al(i))
        q1_ad(i-1) = q1_ad(i-1) + p1*al_ad(i)
        q1_ad(i) = q1_ad(i) + p1*al_ad(i)
        q1_ad(i-2) = q1_ad(i-2) + p2*al_ad(i)
        q1_ad(i+1) = q1_ad(i+1) + p2*al_ad(i)
        al_ad(i) = 0.0
      END DO
 110  DO i=ied,isd,-1
        CALL POPREALARRAY(q1(i))
        q_ad(i, j) = q_ad(i, j) + q1_ad(i)
        q1_ad(i) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
  END SUBROUTINE XPPM_BWD
!  Differentiation of yppm in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.
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
!   gradient     of useful results: q flux c
!   with respect to varying inputs: q flux c
  SUBROUTINE YPPM_FWD(flux, q, c, jord, ifirst, ilast, isd, ied, js, &
&   je, jsd, jed, npx, npy, dya, nested, grid_type)
    IMPLICIT NONE
! Compute domain
    INTEGER, INTENT(IN) :: ifirst, ilast
    INTEGER, INTENT(IN) :: isd, ied, js, je, jsd, jed
    INTEGER, INTENT(IN) :: jord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst:ilast, jsd:jed)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
!  Flux
    REAL :: flux(ifirst:ilast, js:je+1)
    REAL, INTENT(IN) :: dya(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
! Local:
    REAL :: dm(ifirst:ilast, js-2:je+2)
    REAL :: al(ifirst:ilast, js-1:je+2)
    REAL, DIMENSION(ifirst:ilast, js-1:je+1) :: bl, br, b0
    REAL :: dq(ifirst:ilast, js-3:je+2)
    REAL, DIMENSION(ifirst:ilast) :: fx0, fx1
    LOGICAL, DIMENSION(ifirst:ilast, js-1:je+1) :: smt5, smt6
    REAL :: x0, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2, r1
    INTEGER :: i, j, js1, je3, je1
    INTRINSIC MAX
    INTRINSIC MIN

    dm = 0.0
    al = 0.0
    bl = 0.0
    br = 0.0
    b0 = 0.0
    dq = 0.0
    fx0 = 0.0
    fx1 = 0.0
    x0 = 0.0
    xt = 0.0
    qtmp = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    r1 = 0.0
    js1 = 0
    je3 = 0
    je1 = 0

    IF (.NOT.nested .AND. grid_type .LT. 3) THEN
      IF (3 .LT. js - 1) THEN
        js1 = js - 1
      ELSE
        js1 = 3
      END IF
      IF (npy - 2 .GT. je + 2) THEN
        CALL PUSHCONTROL(1,1)
        je3 = je + 2
      ELSE
        CALL PUSHCONTROL(1,1)
        je3 = npy - 2
      END IF
    ELSE
      CALL PUSHCONTROL(1,0)
! Nested grid OR Doubly periodic domain:
      js1 = js - 1
      je3 = je + 2
    END IF
    IF (jord .LT. 8 .OR. jord .EQ. 333) THEN
      DO j=js1,je3
        DO i=ifirst,ilast
          al(i, j) = p1*(q(i, j-1)+q(i, j)) + p2*(q(i, j-2)+q(i, j+1))
        END DO
      END DO
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            al(i, 0) = c1*q(i, -2) + c2*q(i, -1) + c3*q(i, 0)
            al(i, 1) = 0.5*(((2.*dya(i, 0)+dya(i, -1))*q(i, 0)-dya(i, 0)&
&             *q(i, -1))/(dya(i, -1)+dya(i, 0))+((2.*dya(i, 1)+dya(i, 2)&
&             )*q(i, 1)-dya(i, 1)*q(i, 2))/(dya(i, 1)+dya(i, 2)))
            al(i, 2) = c3*q(i, 1) + c2*q(i, 2) + c1*q(i, 3)
          END DO
          CALL PUSHCONTROL(1,0)
        ELSE
          CALL PUSHCONTROL(1,1)
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            al(i, npy-1) = c1*q(i, npy-3) + c2*q(i, npy-2) + c3*q(i, npy&
&             -1)
            al(i, npy) = 0.5*(((2.*dya(i, npy-1)+dya(i, npy-2))*q(i, npy&
&             -1)-dya(i, npy-1)*q(i, npy-2))/(dya(i, npy-2)+dya(i, npy-1&
&             ))+((2.*dya(i, npy)+dya(i, npy+1))*q(i, npy)-dya(i, npy)*q&
&             (i, npy+1))/(dya(i, npy)+dya(i, npy+1)))
            al(i, npy+1) = c3*q(i, npy) + c2*q(i, npy+1) + c1*q(i, npy+2&
&             )
          END DO
          CALL PUSHCONTROL(2,0)
        ELSE
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,2)
      END IF
      IF (jord .EQ. 1) THEN
        DO j=js,je+1
          DO i=ifirst,ilast
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = q(i, j-1)
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = q(i, j)
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
        END DO
        CALL PUSHINTEGER(js1)
        CALL PUSHINTEGER(je3)
        CALL PUSHCONTROL(3,1)
      ELSE IF (jord .EQ. 2) THEN
! Perfectly linear scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              qtmp = q(i, j-1)
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = qtmp + (1.-xt)*(al(i, j)-qtmp-xt*(al(i, j-1)+&
&               al(i, j)-(qtmp+qtmp)))
              CALL PUSHCONTROL(1,1)
            ELSE
              qtmp = q(i, j)
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = qtmp + (1.+xt)*(al(i, j)-qtmp+xt*(al(i, j)+al&
&               (i, j+1)-(qtmp+qtmp)))
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
        END DO
        CALL PUSHINTEGER(js1)
        CALL PUSHINTEGER(je3)
        CALL PUSHREALARRAY(al, (ilast-ifirst+1)*(je-js+4))
        CALL PUSHCONTROL(3,2)
      ELSE IF (jord .EQ. 333) THEN
! Perfectly linear scheme, more diffusive than ord=2 (HoldawayKent-2015-TellusA)
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=ifirst,ilast
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = (2.0*q(i, j)+5.0*q(i, j-1)-q(i, j-2))/6.0 - &
&               0.5*xt*(q(i, j)-q(i, j-1)) + xt*xt/6.0*(q(i, j)-2.0*q(i&
&               , j-1)+q(i, j-2))
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHREALARRAY(flux(i, j))
              flux(i, j) = (2.0*q(i, j-1)+5.0*q(i, j)-q(i, j+1))/6.0 - &
&               0.5*xt*(q(i, j)-q(i, j-1)) + xt*xt/6.0*(q(i, j+1)-2.0*q(&
&               i, j)+q(i, j-1))
              CALL PUSHCONTROL(1,0)
            END IF
          END DO
        END DO
        CALL PUSHINTEGER(js1)
        CALL PUSHINTEGER(je3)
        CALL PUSHCONTROL(3,3)
      ELSE
        CALL PUSHINTEGER(js1)
        CALL PUSHINTEGER(je3)
        CALL PUSHCONTROL(3,4)
      END IF
    ELSE
      CALL PUSHCONTROL(3,0)
    END IF
  END SUBROUTINE YPPM_FWD
!  Differentiation of yppm in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
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
!   gradient     of useful results: q flux c
!   with respect to varying inputs: q flux c
  SUBROUTINE YPPM_BWD(flux, flux_ad, q, q_ad, c, c_ad, jord, ifirst, &
&   ilast, isd, ied, js, je, jsd, jed, npx, npy, dya, nested, grid_type)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ifirst, ilast
    INTEGER, INTENT(IN) :: isd, ied, js, je, jsd, jed
    INTEGER, INTENT(IN) :: jord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst:ilast, jsd:jed)
    REAL :: q_ad(ifirst:ilast, jsd:jed)
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
    REAL :: c_ad(isd:ied, js:je+1)
    REAL :: flux(ifirst:ilast, js:je+1)
    REAL :: flux_ad(ifirst:ilast, js:je+1)
    REAL, INTENT(IN) :: dya(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
    REAL :: dm(ifirst:ilast, js-2:je+2)
    REAL :: al(ifirst:ilast, js-1:je+2)
    REAL :: al_ad(ifirst:ilast, js-1:je+2)
    REAL, DIMENSION(ifirst:ilast, js-1:je+1) :: bl, br, b0
    REAL :: dq(ifirst:ilast, js-3:je+2)
    REAL, DIMENSION(ifirst:ilast) :: fx0, fx1
    LOGICAL, DIMENSION(ifirst:ilast, js-1:je+1) :: smt5, smt6
    REAL :: x0, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2, r1
    REAL :: xt_ad, qtmp_ad
    INTEGER :: i, j, js1, je3, je1
    INTRINSIC MAX
    INTRINSIC MIN
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp0
    REAL :: temp_ad5
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    INTEGER :: branch

    dm = 0.0
    al = 0.0
    bl = 0.0
    br = 0.0
    b0 = 0.0
    dq = 0.0
    fx0 = 0.0
    fx1 = 0.0
    x0 = 0.0
    xt = 0.0
    qtmp = 0.0
    pmp_1 = 0.0
    lac_1 = 0.0
    pmp_2 = 0.0
    lac_2 = 0.0
    r1 = 0.0
    js1 = 0
    je3 = 0
    je1 = 0
    branch = 0

    CALL POPCONTROL(3,branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) THEN
        GOTO 110
      ELSE
        CALL POPINTEGER(je3)
        CALL POPINTEGER(js1)
        DO j=je+1,js,-1
          DO i=ilast,ifirst,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY(flux(i, j))
              q_ad(i, j) = q_ad(i, j) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              CALL POPREALARRAY(flux(i, j))
              q_ad(i, j-1) = q_ad(i, j-1) + flux_ad(i, j)
              flux_ad(i, j) = 0.0
            END IF
          END DO
        END DO
        al_ad = 0.0
      END IF
    ELSE IF (branch .EQ. 2) THEN
      CALL POPREALARRAY(al, (ilast-ifirst+1)*(je-js+4))
      CALL POPINTEGER(je3)
      CALL POPINTEGER(js1)
      al_ad = 0.0
      DO j=je+1,js,-1
        DO i=ilast,ifirst,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            xt = c(i, j)
            qtmp = q(i, j)
            CALL POPREALARRAY(flux(i, j))
            temp0 = al(i, j) + al(i, j+1) - 2*qtmp
            temp_ad5 = (xt+1.)*flux_ad(i, j)
            temp_ad6 = xt*temp_ad5
            qtmp_ad = flux_ad(i, j) - temp_ad5 - 2*temp_ad6
            xt_ad = temp0*temp_ad5 + (al(i, j)-qtmp+xt*temp0)*flux_ad(i&
&             , j)
            al_ad(i, j) = al_ad(i, j) + temp_ad6 + temp_ad5
            al_ad(i, j+1) = al_ad(i, j+1) + temp_ad6
            flux_ad(i, j) = 0.0
            q_ad(i, j) = q_ad(i, j) + qtmp_ad
          ELSE
            xt = c(i, j)
            qtmp = q(i, j-1)
            CALL POPREALARRAY(flux(i, j))
            temp = al(i, j-1) + al(i, j) - 2*qtmp
            temp_ad3 = (1.-xt)*flux_ad(i, j)
            temp_ad4 = -(xt*temp_ad3)
            qtmp_ad = flux_ad(i, j) - temp_ad3 - 2*temp_ad4
            xt_ad = -(temp*temp_ad3) - (al(i, j)-qtmp-xt*temp)*flux_ad(i&
&             , j)
            al_ad(i, j) = al_ad(i, j) + temp_ad4 + temp_ad3
            al_ad(i, j-1) = al_ad(i, j-1) + temp_ad4
            flux_ad(i, j) = 0.0
            q_ad(i, j-1) = q_ad(i, j-1) + qtmp_ad
          END IF
          c_ad(i, j) = c_ad(i, j) + xt_ad
        END DO
      END DO
    ELSE
      IF (branch .EQ. 3) THEN
        CALL POPINTEGER(je3)
        CALL POPINTEGER(js1)
        DO j=je+1,js,-1
          DO i=ilast,ifirst,-1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              xt = c(i, j)
              CALL POPREALARRAY(flux(i, j))
              temp_ad10 = flux_ad(i, j)/6.0
              temp_ad11 = -(0.5*xt*flux_ad(i, j))
              temp_ad12 = xt**2*flux_ad(i, j)/6.0
              q_ad(i, j-1) = q_ad(i, j-1) + temp_ad12 - temp_ad11 + 2.0*&
&               temp_ad10
              q_ad(i, j) = q_ad(i, j) + temp_ad11 - 2.0*temp_ad12 + 5.0*&
&               temp_ad10
              q_ad(i, j+1) = q_ad(i, j+1) + temp_ad12 - temp_ad10
              xt_ad = ((q(i, j+1)-2.0*q(i, j)+q(i, j-1))*2*xt/6.0-0.5*(q&
&               (i, j)-q(i, j-1)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0
            ELSE
              xt = c(i, j)
              CALL POPREALARRAY(flux(i, j))
              temp_ad7 = flux_ad(i, j)/6.0
              temp_ad8 = -(0.5*xt*flux_ad(i, j))
              temp_ad9 = xt**2*flux_ad(i, j)/6.0
              q_ad(i, j) = q_ad(i, j) + temp_ad9 + temp_ad8 + 2.0*&
&               temp_ad7
              q_ad(i, j-1) = q_ad(i, j-1) + 5.0*temp_ad7 - temp_ad8 - &
&               2.0*temp_ad9
              q_ad(i, j-2) = q_ad(i, j-2) + temp_ad9 - temp_ad7
              xt_ad = ((q(i, j)-2.0*q(i, j-1)+q(i, j-2))*2*xt/6.0-0.5*(q&
&               (i, j)-q(i, j-1)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0
            END IF
            c_ad(i, j) = c_ad(i, j) + xt_ad
          END DO
        END DO
      ELSE
        CALL POPINTEGER(je3)
        CALL POPINTEGER(js1)
      END IF
      al_ad = 0.0
    END IF
    CALL POPCONTROL(2,branch)
    IF (branch .EQ. 0) THEN
      DO i=ilast,ifirst,-1
        q_ad(i, npy) = q_ad(i, npy) + c3*al_ad(i, npy+1)
        q_ad(i, npy+1) = q_ad(i, npy+1) + c2*al_ad(i, npy+1)
        q_ad(i, npy+2) = q_ad(i, npy+2) + c1*al_ad(i, npy+1)
        al_ad(i, npy+1) = 0.0
        temp_ad1 = 0.5*al_ad(i, npy)/(dya(i, npy-2)+dya(i, npy-1))
        temp_ad2 = 0.5*al_ad(i, npy)/(dya(i, npy)+dya(i, npy+1))
        q_ad(i, npy-1) = q_ad(i, npy-1) + (dya(i, npy-1)*2.+dya(i, npy-2&
&         ))*temp_ad1
        q_ad(i, npy-2) = q_ad(i, npy-2) - dya(i, npy-1)*temp_ad1
        q_ad(i, npy) = q_ad(i, npy) + (dya(i, npy)*2.+dya(i, npy+1))*&
&         temp_ad2
        q_ad(i, npy+1) = q_ad(i, npy+1) - dya(i, npy)*temp_ad2
        al_ad(i, npy) = 0.0
        q_ad(i, npy-3) = q_ad(i, npy-3) + c1*al_ad(i, npy-1)
        q_ad(i, npy-2) = q_ad(i, npy-2) + c2*al_ad(i, npy-1)
        q_ad(i, npy-1) = q_ad(i, npy-1) + c3*al_ad(i, npy-1)
        al_ad(i, npy-1) = 0.0
      END DO
    ELSE IF (branch .NE. 1) THEN
      GOTO 100
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO i=ilast,ifirst,-1
        q_ad(i, 1) = q_ad(i, 1) + c3*al_ad(i, 2)
        q_ad(i, 2) = q_ad(i, 2) + c2*al_ad(i, 2)
        q_ad(i, 3) = q_ad(i, 3) + c1*al_ad(i, 2)
        al_ad(i, 2) = 0.0
        temp_ad = 0.5*al_ad(i, 1)/(dya(i, -1)+dya(i, 0))
        temp_ad0 = 0.5*al_ad(i, 1)/(dya(i, 1)+dya(i, 2))
        q_ad(i, 0) = q_ad(i, 0) + (dya(i, 0)*2.+dya(i, -1))*temp_ad
        q_ad(i, -1) = q_ad(i, -1) - dya(i, 0)*temp_ad
        q_ad(i, 1) = q_ad(i, 1) + (dya(i, 1)*2.+dya(i, 2))*temp_ad0
        q_ad(i, 2) = q_ad(i, 2) - dya(i, 1)*temp_ad0
        al_ad(i, 1) = 0.0
        q_ad(i, -2) = q_ad(i, -2) + c1*al_ad(i, 0)
        q_ad(i, -1) = q_ad(i, -1) + c2*al_ad(i, 0)
        q_ad(i, 0) = q_ad(i, 0) + c3*al_ad(i, 0)
        al_ad(i, 0) = 0.0
      END DO
    END IF
 100 DO j=je3,js1,-1
      DO i=ilast,ifirst,-1
        q_ad(i, j-1) = q_ad(i, j-1) + p1*al_ad(i, j)
        q_ad(i, j) = q_ad(i, j) + p1*al_ad(i, j)
        q_ad(i, j-2) = q_ad(i, j-2) + p2*al_ad(i, j)
        q_ad(i, j+1) = q_ad(i, j+1) + p2*al_ad(i, j)
        al_ad(i, j) = 0.0
      END DO
    END DO
 110 CALL POPCONTROL(1,branch)
  END SUBROUTINE YPPM_BWD
!  Differentiation of deln_flux in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: q mass fx fy
!   with respect to varying inputs: q mass fx fy
  SUBROUTINE DELN_FLUX_FWD(nord, is, ie, js, je, npx, npy, damp, q, &
&   fx, fy, gridstruct, bd, mass)
    IMPLICIT NONE
! Del-n damping for the cell-mean values (A grid)
!------------------
! nord = 0:   del-2
! nord = 1:   del-4
! nord = 2:   del-6
! nord = 3:   del-8 --> requires more ghosting than current
!------------------
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! del-n
    INTEGER, INTENT(IN) :: nord
    INTEGER, INTENT(IN) :: is, ie, js, je, npx, npy
    REAL, INTENT(IN) :: damp
! q ghosted on input
    REAL, INTENT(IN) :: q(bd%is-ng:bd%ie+ng, bd%js-ng:bd%je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! q ghosted on input
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
! diffusive fluxes:
    REAL, INTENT(INOUT) :: fx(bd%is:bd%ie+1, bd%js:bd%je), fy(bd%is:bd%&
&   ie, bd%js:bd%je+1)
! local:
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2(bd%isd:bd%ied, bd%&
&   jsd:bd%jed+1)
    REAL :: d2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: damp2
    INTEGER :: i, j, n, nt, i1, i2, j1, j2
    INTRINSIC PRESENT
    INTEGER :: ad_from
    INTEGER :: ad_from0
    INTEGER :: ad_from1
    INTEGER :: ad_from2
    INTEGER :: ad_from3
    INTEGER :: ad_from4

    fx2 = 0.0
    fy2 = 0.0
    d2 = 0.0
    damp2 = 0.0
    nt = 0
    i1 = 0
    i2 = 0
    j1 = 0
    j2 = 0
    ad_from = 0
    ad_from0 = 0
    ad_from1 = 0
    ad_from2 = 0
    ad_from3 = 0
    ad_from4 = 0

    i1 = is - 1 - nord
    i2 = ie + 1 + nord
    j1 = js - 1 - nord
    j2 = je + 1 + nord
    IF (.NOT.PRESENT(mass)) THEN
      DO j=j1,j2
        DO i=i1,i2
          d2(i, j) = damp*q(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      DO j=j1,j2
        DO i=i1,i2
          d2(i, j) = q(i, j)
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    END IF
    IF (nord .GT. 0) THEN
      CALL COPY_CORNERS_FWD(d2, npx, npy, 1, gridstruct%nested, bd, &
&                        gridstruct%sw_corner, gridstruct%se_corner, &
&                        gridstruct%nw_corner, gridstruct%ne_corner)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
    DO j=js-nord,je+nord
      DO i=is-nord,ie+nord+1
        fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i-1, j)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) THEN
      CALL COPY_CORNERS_FWD(d2, npx, npy, 2, gridstruct%nested, bd, &
&                        gridstruct%sw_corner, gridstruct%se_corner, &
&                        gridstruct%nw_corner, gridstruct%ne_corner)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
    DO j=js-nord,je+nord+1
      DO i=is-nord,ie+nord
        fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j-1)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) THEN
!----------
! high-order
!----------
      DO n=1,nord
        nt = nord - n
        ad_from0 = js - nt - 1
        DO j=ad_from0,je+nt+1
          ad_from = is - nt - 1
          DO i=ad_from,ie+nt+1
            d2(i, j) = (fx2(i, j)-fx2(i+1, j)+fy2(i, j)-fy2(i, j+1))*&
&             gridstruct%rarea(i, j)
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from)
        END DO
        CALL PUSHINTEGER(j - 1)
        CALL PUSHINTEGER(ad_from0)
        CALL COPY_CORNERS_FWD(d2, npx, npy, 1, gridstruct%nested, bd&
&                          , gridstruct%sw_corner, gridstruct%se_corner&
&                          , gridstruct%nw_corner, gridstruct%ne_corner)
        ad_from2 = js - nt
        DO j=ad_from2,je+nt
          ad_from1 = is - nt
          DO i=ad_from1,ie+nt+1
            fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i, j)-d2(i-1, j))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from1)
        END DO
        CALL PUSHINTEGER(j - 1)
        CALL PUSHINTEGER(ad_from2)
        CALL COPY_CORNERS_FWD(d2, npx, npy, 2, gridstruct%nested, bd&
&                          , gridstruct%sw_corner, gridstruct%se_corner&
&                          , gridstruct%nw_corner, gridstruct%ne_corner)
        ad_from4 = js - nt
        DO j=ad_from4,je+nt+1
          ad_from3 = is - nt
          DO i=ad_from3,ie+nt
            fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j)-d2(i, j-1))
          END DO
          CALL PUSHINTEGER(i - 1)
          CALL PUSHINTEGER(ad_from3)
        END DO
        CALL PUSHINTEGER(j - 1)
        CALL PUSHINTEGER(ad_from4)
      END DO
      CALL PUSHCONTROL(1,0)
    ELSE
      CALL PUSHCONTROL(1,1)
    END IF
!---------------------------------------------
! Add the diffusive fluxes to the flux arrays:
!---------------------------------------------
    IF (PRESENT(mass)) THEN
! Apply mass weighting to diffusive fluxes:
      damp2 = 0.5*damp
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(fx(i, j))
          fx(i, j) = fx(i, j) + damp2*(mass(i-1, j)+mass(i, j))*fx2(i, j&
&           )
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(fy(i, j))
          fy(i, j) = fy(i, j) + damp2*(mass(i, j-1)+mass(i, j))*fy2(i, j&
&           )
        END DO
      END DO
      CALL PUSHREALARRAY(damp2)
      CALL PUSHREALARRAY(fx2, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL PUSHREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      CALL PUSHINTEGER(i1)
      CALL PUSHCONTROL(1,0)
    ELSE
      DO j=js,je
        DO i=is,ie+1
          CALL PUSHREALARRAY(fx(i, j))
          fx(i, j) = fx(i, j) + fx2(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          CALL PUSHREALARRAY(fy(i, j))
          fy(i, j) = fy(i, j) + fy2(i, j)
        END DO
      END DO
      CALL PUSHINTEGER(i1)
      CALL PUSHCONTROL(1,1)
    END IF
  END SUBROUTINE DELN_FLUX_FWD
!  Differentiation of deln_flux in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: q mass fx fy
!   with respect to varying inputs: q mass fx fy
  SUBROUTINE DELN_FLUX_BWD(nord, is, ie, js, je, npx, npy, damp, q, &
&   q_ad, fx, fx_ad, fy, fy_ad, gridstruct, bd, mass, mass_ad)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: nord
    INTEGER, INTENT(IN) :: is, ie, js, je, npx, npy
    REAL, INTENT(IN) :: damp
    REAL, INTENT(IN) :: q(bd%is-ng:bd%ie+ng, bd%js-ng:bd%je+ng)
    REAL :: q_ad(bd%is-ng:bd%ie+ng, bd%js-ng:bd%je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL :: mass_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: fx(bd%is:bd%ie+1, bd%js:bd%je), fy(bd%is:bd%&
&   ie, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: fx_ad(bd%is:bd%ie+1, bd%js:bd%je), fy_ad(bd%&
&   is:bd%ie, bd%js:bd%je+1)
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2(bd%isd:bd%ied, bd%&
&   jsd:bd%jed+1)
    REAL :: fx2_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2_ad(bd%isd:bd%ied&
&   , bd%jsd:bd%jed+1)
    REAL :: d2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: d2_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: damp2
    INTEGER :: i, j, n, nt, i1, i2, j1, j2
    INTRINSIC PRESENT
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
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

    fx2 = 0.0
    fy2 = 0.0
    d2 = 0.0
    damp2 = 0.0
    nt = 0
    i1 = 0
    i2 = 0
    j1 = 0
    j2 = 0
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

    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(i1)
      CALL POPREALARRAY(fy2, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+2))
      CALL POPREALARRAY(fx2, (bd%ied-bd%isd+2)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(damp2)
      fy2_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(fy(i, j))
          temp_ad5 = damp2*fy2(i, j)*fy_ad(i, j)
          mass_ad(i, j-1) = mass_ad(i, j-1) + temp_ad5
          mass_ad(i, j) = mass_ad(i, j) + temp_ad5
          fy2_ad(i, j) = fy2_ad(i, j) + damp2*(mass(i, j-1)+mass(i, j))*&
&           fy_ad(i, j)
        END DO
      END DO
      fx2_ad = 0.0
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(fx(i, j))
          temp_ad4 = damp2*fx2(i, j)*fx_ad(i, j)
          mass_ad(i-1, j) = mass_ad(i-1, j) + temp_ad4
          mass_ad(i, j) = mass_ad(i, j) + temp_ad4
          fx2_ad(i, j) = fx2_ad(i, j) + damp2*(mass(i-1, j)+mass(i, j))*&
&           fx_ad(i, j)
        END DO
      END DO
    ELSE
      CALL POPINTEGER(i1)
      fy2_ad = 0.0
      DO j=je+1,js,-1
        DO i=ie,is,-1
          CALL POPREALARRAY(fy(i, j))
          fy2_ad(i, j) = fy2_ad(i, j) + fy_ad(i, j)
        END DO
      END DO
      fx2_ad = 0.0
      DO j=je,js,-1
        DO i=ie+1,is,-1
          CALL POPREALARRAY(fx(i, j))
          fx2_ad(i, j) = fx2_ad(i, j) + fx_ad(i, j)
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      d2_ad = 0.0
      DO n=nord,1,-1
        CALL POPINTEGER(ad_from4)
        CALL POPINTEGER(ad_to4)
        DO j=ad_to4,ad_from4,-1
          CALL POPINTEGER(ad_from3)
          CALL POPINTEGER(ad_to3)
          DO i=ad_to3,ad_from3,-1
            temp_ad3 = gridstruct%del6_u(i, j)*fy2_ad(i, j)
            d2_ad(i, j) = d2_ad(i, j) + temp_ad3
            d2_ad(i, j-1) = d2_ad(i, j-1) - temp_ad3
            fy2_ad(i, j) = 0.0
          END DO
        END DO
        CALL COPY_CORNERS_BWD(d2, d2_ad, npx, npy, 2, gridstruct%&
&                          nested, bd, gridstruct%sw_corner, gridstruct%&
&                          se_corner, gridstruct%nw_corner, gridstruct%&
&                          ne_corner)
        CALL POPINTEGER(ad_from2)
        CALL POPINTEGER(ad_to2)
        DO j=ad_to2,ad_from2,-1
          CALL POPINTEGER(ad_from1)
          CALL POPINTEGER(ad_to1)
          DO i=ad_to1,ad_from1,-1
            temp_ad2 = gridstruct%del6_v(i, j)*fx2_ad(i, j)
            d2_ad(i, j) = d2_ad(i, j) + temp_ad2
            d2_ad(i-1, j) = d2_ad(i-1, j) - temp_ad2
            fx2_ad(i, j) = 0.0
          END DO
        END DO
        CALL COPY_CORNERS_BWD(d2, d2_ad, npx, npy, 1, gridstruct%&
&                          nested, bd, gridstruct%sw_corner, gridstruct%&
&                          se_corner, gridstruct%nw_corner, gridstruct%&
&                          ne_corner)
        CALL POPINTEGER(ad_from0)
        CALL POPINTEGER(ad_to0)
        DO j=ad_to0,ad_from0,-1
          CALL POPINTEGER(ad_from)
          CALL POPINTEGER(ad_to)
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
    ELSE
      d2_ad = 0.0
    END IF
    DO j=je+nord+1,js-nord,-1
      DO i=ie+nord,is-nord,-1
        temp_ad0 = gridstruct%del6_u(i, j)*fy2_ad(i, j)
        d2_ad(i, j-1) = d2_ad(i, j-1) + temp_ad0
        d2_ad(i, j) = d2_ad(i, j) - temp_ad0
        fy2_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) CALL COPY_CORNERS_BWD(d2, d2_ad, npx, npy, 2, &
&                                         gridstruct%nested, bd, &
&                                         gridstruct%sw_corner, &
&                                         gridstruct%se_corner, &
&                                         gridstruct%nw_corner, &
&                                         gridstruct%ne_corner)
    DO j=je+nord,js-nord,-1
      DO i=ie+nord+1,is-nord,-1
        temp_ad = gridstruct%del6_v(i, j)*fx2_ad(i, j)
        d2_ad(i-1, j) = d2_ad(i-1, j) + temp_ad
        d2_ad(i, j) = d2_ad(i, j) - temp_ad
        fx2_ad(i, j) = 0.0
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) CALL COPY_CORNERS_BWD(d2, d2_ad, npx, npy, 1, &
&                                         gridstruct%nested, bd, &
&                                         gridstruct%sw_corner, &
&                                         gridstruct%se_corner, &
&                                         gridstruct%nw_corner, &
&                                         gridstruct%ne_corner)
    i2 = ie + 1 + nord
    j1 = js - 1 - nord
    j2 = je + 1 + nord
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      DO j=j2,j1,-1
        DO i=i2,i1,-1
          q_ad(i, j) = q_ad(i, j) + damp*d2_ad(i, j)
          d2_ad(i, j) = 0.0
        END DO
      END DO
    ELSE
      DO j=j2,j1,-1
        DO i=i2,i1,-1
          q_ad(i, j) = q_ad(i, j) + d2_ad(i, j)
          d2_ad(i, j) = 0.0
        END DO
      END DO
    END IF
  END SUBROUTINE DELN_FLUX_BWD
!  Differentiation of copy_corners in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_e
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
!   gradient     of useful results: q
!   with respect to varying inputs: q
!Weird arguments are because this routine is called in a lot of
!places outside of tp_core, sometimes very deeply nested in the call tree.
  SUBROUTINE COPY_CORNERS_FWD(q, npx, npy, dir, nested, bd, sw_corner&
&   , se_corner, nw_corner, ne_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy, dir
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: nested, sw_corner, se_corner, nw_corner, &
&   ne_corner
    INTEGER :: i, j
    REAL :: tmp
    REAL :: tmp0
    REAL :: tmp1
    REAL :: tmp2
    REAL :: tmp3
    REAL :: tmp4
    REAL :: tmp5
    REAL :: tmp6
    IF (nested) THEN
      CALL PUSHCONTROL(3,0)
    ELSE IF (dir .EQ. 1) THEN
! XDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            tmp = q(j, 1-i)
            CALL PUSHREALARRAY(q(i, j))
            q(i, j) = tmp
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            tmp0 = q(npy-j, i-npx+1)
            CALL PUSHREALARRAY(q(i, j))
            q(i, j) = tmp0
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            tmp1 = q(j, 2*npx-1-i)
            CALL PUSHREALARRAY(q(i, j))
            q(i, j) = tmp1
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            tmp2 = q(npy-j, i-1+npx)
            CALL PUSHREALARRAY(q(i, j))
            q(i, j) = tmp2
          END DO
        END DO
        CALL PUSHCONTROL(3,2)
      ELSE
        CALL PUSHCONTROL(3,1)
      END IF
    ELSE IF (dir .EQ. 2) THEN
! YDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            tmp3 = q(1-j, i)
            CALL PUSHREALARRAY(q(i, j))
            q(i, j) = tmp3
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            tmp4 = q(npy+j-1, npx-i)
            CALL PUSHREALARRAY(q(i, j))
            q(i, j) = tmp4
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            tmp5 = q(2*npy-1-j, i)
            CALL PUSHREALARRAY(q(i, j))
            q(i, j) = tmp5
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            tmp6 = q(j+1-npx, npy-i)
            CALL PUSHREALARRAY(q(i, j))
            q(i, j) = tmp6
          END DO
        END DO
        CALL PUSHCONTROL(3,5)
      ELSE
        CALL PUSHCONTROL(3,4)
      END IF
    ELSE
      CALL PUSHCONTROL(3,3)
    END IF
  END SUBROUTINE COPY_CORNERS_FWD
!  Differentiation of copy_corners in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_
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
!   gradient     of useful results: q
!   with respect to varying inputs: q
!Weird arguments are because this routine is called in a lot of
!places outside of tp_core, sometimes very deeply nested in the call tree.
  SUBROUTINE COPY_CORNERS_BWD(q, q_ad, npx, npy, dir, nested, bd, &
&   sw_corner, se_corner, nw_corner, ne_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy, dir
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed)
    LOGICAL, INTENT(IN) :: nested, sw_corner, se_corner, nw_corner, &
&   ne_corner
    INTEGER :: i, j
    REAL :: tmp_ad
    REAL :: tmp_ad0
    REAL :: tmp_ad1
    REAL :: tmp_ad2
    REAL :: tmp_ad3
    REAL :: tmp_ad4
    REAL :: tmp_ad5
    REAL :: tmp_ad6
    INTEGER :: branch

    branch = 0

    CALL POPCONTROL(3,branch)
    IF (branch .LT. 3) THEN
      IF (branch .NE. 0) THEN
        IF (branch .NE. 1) THEN
          DO j=npy+ng-1,npy,-1
            DO i=0,1-ng,-1
              CALL POPREALARRAY(q(i, j))
              tmp_ad2 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(npy-j, i-1+npx) = q_ad(npy-j, i-1+npx) + tmp_ad2
            END DO
          END DO
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=npy+ng-1,npy,-1
            DO i=npx+ng-1,npx,-1
              CALL POPREALARRAY(q(i, j))
              tmp_ad1 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(j, 2*npx-1-i) = q_ad(j, 2*npx-1-i) + tmp_ad1
            END DO
          END DO
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=0,1-ng,-1
            DO i=npx+ng-1,npx,-1
              CALL POPREALARRAY(q(i, j))
              tmp_ad0 = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(npy-j, i-npx+1) = q_ad(npy-j, i-npx+1) + tmp_ad0
            END DO
          END DO
        END IF
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=0,1-ng,-1
            DO i=0,1-ng,-1
              CALL POPREALARRAY(q(i, j))
              tmp_ad = q_ad(i, j)
              q_ad(i, j) = 0.0
              q_ad(j, 1-i) = q_ad(j, 1-i) + tmp_ad
            END DO
          END DO
        END IF
      END IF
    ELSE IF (branch .NE. 3) THEN
      IF (branch .NE. 4) THEN
        DO j=npy+ng-1,npy,-1
          DO i=0,1-ng,-1
            CALL POPREALARRAY(q(i, j))
            tmp_ad6 = q_ad(i, j)
            q_ad(i, j) = 0.0
            q_ad(j+1-npx, npy-i) = q_ad(j+1-npx, npy-i) + tmp_ad6
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=npy+ng-1,npy,-1
          DO i=npx+ng-1,npx,-1
            CALL POPREALARRAY(q(i, j))
            tmp_ad5 = q_ad(i, j)
            q_ad(i, j) = 0.0
            q_ad(2*npy-1-j, i) = q_ad(2*npy-1-j, i) + tmp_ad5
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=0,1-ng,-1
          DO i=npx+ng-1,npx,-1
            CALL POPREALARRAY(q(i, j))
            tmp_ad4 = q_ad(i, j)
            q_ad(i, j) = 0.0
            q_ad(npy+j-1, npx-i) = q_ad(npy+j-1, npx-i) + tmp_ad4
          END DO
        END DO
      END IF
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=0,1-ng,-1
          DO i=0,1-ng,-1
            CALL POPREALARRAY(q(i, j))
            tmp_ad3 = q_ad(i, j)
            q_ad(i, j) = 0.0
            q_ad(1-j, i) = q_ad(1-j, i) + tmp_ad3
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE COPY_CORNERS_BWD
!Weird arguments are because this routine is called in a lot of
!places outside of tp_core, sometimes very deeply nested in the call tree.
end module tp_core_adm_mod
