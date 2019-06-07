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
module fv_tracer2d_adm_mod
   use tp_core_adm_mod,   only: fv_tp_2d
   use tp_core_adm_mod,   only: fv_tp_2d_fwd, fv_tp_2d_bwd, fv_tp_2d_adm
   use fv_mp_mod,         only: mp_reduce_max
   use fv_mp_mod,         only: ng, mp_gather, is_master
   use fv_mp_mod,         only: group_halo_update_type
   use fv_mp_adm_mod,     only: start_group_halo_update, complete_group_halo_update
   use fv_mp_adm_mod,     only: start_group_halo_update_adm
   use mpp_domains_mod,   only: mpp_update_domains, mpp_get_boundary, CGRID_NE, domain2d
   use fv_mp_adm_mod,     only: mpp_update_domains_adm
   use fv_timing_mod,     only: timing_on, timing_off
   use boundary_adm_mod,  only: nested_grid_BC_apply_intT
   use boundary_adm_mod,  only: nested_grid_BC_apply_intT_adm
   use fv_arrays_mod,     only: fv_grid_type, fv_flags_type, fv_nest_type, fv_atmos_type, fv_grid_bounds_type
   use mpp_mod,           only: mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum, mpp_max

   use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                            pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm
implicit none
private

public :: tracer_2d, tracer_2d_nested, tracer_2d_1L
public :: tracer_2d_fwd, tracer_2d_nested_fwd, tracer_2d_1L_fwd
public :: tracer_2d_bwd, tracer_2d_nested_bwd, tracer_2d_1L_bwd

real, allocatable, dimension(:,:,:) :: nest_fx_west_accum, nest_fx_east_accum, nest_fx_south_accum, nest_fx_north_accum

!---- version number -----
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of tracer_2d_1l in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: q dp1
!   with respect to varying inputs: q dp1 mfx mfy cx cy
!-----------------------------------------------------------------------
! !ROUTINE: Perform 2D horizontal-to-lagrangian transport
!-----------------------------------------------------------------------
  SUBROUTINE TRACER_2D_1L_FWD(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, &
&   domain, npx, npy, npz, nq, hord, q_split, dt, id_divg, q_pack, &
&   nord_tr, trdm, hord_pert, nord_tr_pert, trdm_pert, split_damp_tr, &
&   dpa)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord, nord_tr
    INTEGER, INTENT(IN) :: hord_pert, nord_tr_pert
    LOGICAL, INTENT(IN) :: split_damp_tr
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt, trdm, trdm_pert
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: q_pack
! Tracers
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
! DELP after advection
    REAL, OPTIONAL, INTENT(OUT) :: dpa(bd%is:bd%ie, bd%js:bd%je)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Local Arrays
! 3D tracers
    REAL :: qn2(bd%isd:bd%ied, bd%jsd:bd%jed, nq)
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cmax(npz)
    REAL :: frac
    INTEGER :: nsplt
    INTEGER :: i, j, k, it, iq
    REAL, DIMENSION(:, :), POINTER :: area, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: dxa, dya, dx, dy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    INTRINSIC PRESENT
    REAL :: max1
    REAL :: x1
    REAL :: z1
    REAL :: y3
    REAL :: y2
    REAL :: y1

    qn2 = 0.0
    dp2 = 0.0
    fx = 0.0
    fy = 0.0
    ra_x = 0.0
    ra_y = 0.0
    xfx = 0.0
    yfx = 0.0
    cmax = 0.0
    frac = 0.0
    nsplt = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    max1 = 0
    x1 = 0
    z1 = 0
    y3 = 0
    y2 = 0
    y1 = 0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    dx => gridstruct%dx
    dy => gridstruct%dy
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx,cmax)
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sin_sg(i-1, &
&             j, 3)
            CALL PUSHCONTROL(1,1)
          ELSE
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1&
&             )
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sin_sg(i, j-&
&             1, 4)
            CALL PUSHCONTROL(1,1)
          ELSE
            yfx(i, j, k) = cy(i, j, k)*dya(i, j)*dx(i, j)*sin_sg(i, j, 2&
&             )
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      cmax(k) = 0.
      IF (k .LT. npz/6) THEN
        DO j=js,je
          DO i=is,ie
            IF (cx(i, j, k) .GE. 0.) THEN
              y1 = cx(i, j, k)
            ELSE
              y1 = -cx(i, j, k)
            END IF
            IF (cy(i, j, k) .GE. 0.) THEN
              z1 = cy(i, j, k)
            ELSE
              z1 = -cy(i, j, k)
            END IF
            IF (cmax(k) .LT. y1) THEN
              IF (y1 .LT. z1) THEN
                CALL PUSHCONTROL(2,0)
                cmax(k) = z1
              ELSE
                CALL PUSHCONTROL(2,1)
                cmax(k) = y1
              END IF
            ELSE IF (cmax(k) .LT. z1) THEN
              CALL PUSHCONTROL(2,2)
              cmax(k) = z1
            ELSE
              CALL PUSHCONTROL(2,3)
              cmax(k) = cmax(k)
            END IF
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        DO j=js,je
          DO i=is,ie
            IF (cx(i, j, k) .GE. 0.) THEN
              x1 = cx(i, j, k)
            ELSE
              x1 = -cx(i, j, k)
            END IF
            IF (cy(i, j, k) .GE. 0.) THEN
              y3 = cy(i, j, k)
            ELSE
              y3 = -cy(i, j, k)
            END IF
            IF (x1 .LT. y3) THEN
              max1 = y3
            ELSE
              max1 = x1
            END IF
            y2 = max1 + 1. - sin_sg(i, j, 5)
            IF (cmax(k) .LT. y2) THEN
              CALL PUSHCONTROL(1,0)
              cmax(k) = y2
            ELSE
              CALL PUSHCONTROL(1,1)
              cmax(k) = cmax(k)
            END IF
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
! k-loop
    CALL MP_REDUCE_MAX(cmax, npz)
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx, &
!$OMP                                  cy,yfx,mfx,mfy,cmax)   &
!$OMP                          private(nsplt, frac)
    DO k=1,npz
      nsplt = INT(1. + cmax(k))
      IF (nsplt .GT. 1) THEN
        CALL PUSHREALARRAY(frac)
        frac = 1./REAL(nsplt)
        DO j=jsd,jed
          DO i=is,ie+1
            CALL PUSHREALARRAY(cx(i, j, k))
            cx(i, j, k) = cx(i, j, k)*frac
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            CALL PUSHREALARRAY(mfx(i, j, k))
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            CALL PUSHREALARRAY(cy(i, j, k))
            cy(i, j, k) = cy(i, j, k)*frac
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            CALL PUSHREALARRAY(mfy(i, j, k))
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
! Begin k-independent tracer transport; can not be OpenMPed because the mpp_update call.
    DO k=1,npz
!$OMP parallel do default(none) shared(k,is,ie,js,je,isd,ied,jsd,jed,xfx,area,yfx,ra_x,ra_y)
      DO j=jsd,jed
        DO i=is,ie
          CALL PUSHREALARRAY(ra_x(i, j))
          ra_x(i, j) = area(i, j) + (xfx(i, j, k)-xfx(i+1, j, k))
        END DO
        IF (j .GE. js .AND. j .LE. je) THEN
          DO i=isd,ied
            CALL PUSHREALARRAY(ra_y(i, j))
            ra_y(i, j) = area(i, j) + (yfx(i, j, k)-yfx(i, j+1, k))
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      nsplt = INT(1. + cmax(k))
      DO it=1,nsplt
!$OMP parallel do default(none) shared(k,is,ie,js,je,rarea,mfx,mfy,dp1,dp2)
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(dp2(i, j))
            dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+(mfy&
&             (i, j, k)-mfy(i, j+1, k)))*rarea(i, j)
          END DO
        END DO
!$OMP parallel do default(none) shared(k,nsplt,it,is,ie,js,je,isd,ied,jsd,jed,npx,npy,cx,xfx,hord,trdm, &
!$OMP                                  nord_tr,nq,gridstruct,bd,cy,yfx,mfx,mfy,qn2,q,ra_x,ra_y,dp1,dp2,rarea) &
!$OMP                          private(fx,fy)
        DO iq=1,nq
          IF (nsplt .NE. 1) THEN
            IF (it .EQ. 1) THEN
              DO j=jsd,jed
                DO i=isd,ied
                  CALL PUSHREALARRAY(qn2(i, j, iq))
                  qn2(i, j, iq) = q(i, j, k, iq)
                END DO
              END DO
              CALL PUSHCONTROL(1,0)
            ELSE
              CALL PUSHCONTROL(1,1)
            END IF
            IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D_FWD(qn2(isd:ied, jsd:jed, iq), cx(is:ie+1&
&                            , jsd:jed, k), cy(isd:ied, js:je+1, k), npx&
&                            , npy, hord, fx, fy, xfx(is:ie+1, jsd:jed, &
&                            k), yfx(isd:ied, js:je+1, k), gridstruct, &
&                            bd, ra_x, ra_y, mfx=mfx(is:ie+1, js:je, k)&
&                            , mfy=mfy(is:ie, js:je+1, k))
              CALL PUSHCONTROL(1,0)
            ELSE
              CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
              CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
              CALL PUSHREALARRAY(qn2(isd:ied, jsd:jed, iq), (ied-isd+1)&
&                           *(jed-jsd+1))
              CALL FV_TP_2D(qn2(isd:ied, jsd:jed, iq), cx(is:ie+1, jsd:&
&                     jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                     hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                     isd:ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, &
&                     mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1&
&                     , k))
              CALL PUSHCONTROL(1,1)
            END IF
            IF (it .LT. nsplt) THEN
! not last call
              DO j=js,je
                DO i=is,ie
                  CALL PUSHREALARRAY(qn2(i, j, iq))
                  qn2(i, j, iq) = (qn2(i, j, iq)*dp1(i, j, k)+(fx(i, j)-&
&                   fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i&
&                   , j)
                END DO
              END DO
              CALL PUSHCONTROL(2,2)
            ELSE
              DO j=js,je
                DO i=is,ie
                  CALL PUSHREALARRAY(q(i, j, k, iq))
                  q(i, j, k, iq) = (qn2(i, j, iq)*dp1(i, j, k)+(fx(i, j)&
&                   -fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(&
&                   i, j)
                END DO
              END DO
              CALL PUSHCONTROL(2,1)
            END IF
          ELSE
            IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D_FWD(q(isd:ied, jsd:jed, k, iq), cx(is:ie+&
&                            1, jsd:jed, k), cy(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fy, xfx(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), &
&                            gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1&
&                            , js:je, k), mfy=mfy(is:ie, js:je+1, k))
              CALL PUSHCONTROL(1,1)
            ELSE
              CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
              CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
              CALL PUSHREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1&
&                           )*(jed-jsd+1))
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd:&
&                     jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                     hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                     isd:ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, &
&                     mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1&
&                     , k))
              CALL PUSHCONTROL(1,0)
            END IF
            DO j=js,je
              DO i=is,ie
                CALL PUSHREALARRAY(q(i, j, k, iq))
                q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-&
&                 fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i, &
&                 j)
              END DO
            END DO
            CALL PUSHCONTROL(2,0)
          END IF
        END DO
!  tracer-loop
        IF (it .LT. nsplt) THEN
! not last call
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(dp1(i, j, k))
              dp1(i, j, k) = dp2(i, j)
            END DO
          END DO
          CALL PUSHREALARRAY(qn2, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*&
&                       nq)
          CALL MPP_UPDATE_DOMAINS(qn2, domain)
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL PUSHINTEGER(it - 1)
    END DO
    CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
    CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
    CALL PUSHREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    CALL PUSHREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
    CALL PUSHREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
    CALL PUSHREALARRAY(qn2, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*nq)
    CALL PUSHREALARRAY(frac)
    CALL PUSHREALARRAY(dp2, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL PUSHREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
  END SUBROUTINE TRACER_2D_1L_FWD
!  Differentiation of tracer_2d_1l in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: q dp1
!   with respect to varying inputs: q dp1 mfx mfy cx cy
!-----------------------------------------------------------------------
! !ROUTINE: Perform 2D horizontal-to-lagrangian transport
!-----------------------------------------------------------------------
  SUBROUTINE TRACER_2D_1L_BWD(q, q_ad, dp1, dp1_ad, mfx, mfx_ad, mfy, &
&   mfy_ad, cx, cx_ad, cy, cy_ad, gridstruct, bd, domain, npx, npy, npz&
&   , nq, hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, hord_pert, &
&   nord_tr_pert, trdm_pert, split_damp_tr, dpa)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord, nord_tr
    INTEGER, INTENT(IN) :: hord_pert, nord_tr_pert
    LOGICAL, INTENT(IN) :: split_damp_tr
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt, trdm, trdm_pert
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: q_pack
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dp1_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_ad(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_ad(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, OPTIONAL, INTENT(OUT) :: dpa(bd%is:bd%ie, bd%js:bd%je)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL :: qn2(bd%isd:bd%ied, bd%jsd:bd%jed, nq)
    REAL :: qn2_ad(bd%isd:bd%ied, bd%jsd:bd%jed, nq)
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: dp2_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_ad(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_ad(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_ad(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_ad(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cmax(npz)
    REAL :: frac
    INTEGER :: nsplt
    INTEGER :: i, j, k, it, iq
    REAL, DIMENSION(:, :), POINTER :: area, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: dxa, dya, dx, dy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    INTRINSIC PRESENT
    REAL :: max1
    REAL :: temp_ad
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp
    REAL :: temp_ad4
    REAL :: temp_ad5
    INTEGER :: branch
    INTEGER :: ad_to
    REAL :: x1
    REAL :: z1
    REAL :: y3
    REAL :: y2
    REAL :: y1

    qn2 = 0.0
    dp2 = 0.0
    fx = 0.0
    fy = 0.0
    ra_x = 0.0
    ra_y = 0.0
    xfx = 0.0
    yfx = 0.0
    cmax = 0.0
    frac = 0.0
    nsplt = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    max1 = 0
    x1 = 0
    z1 = 0
    y3 = 0
    y2 = 0
    y1 = 0
    ad_to = 0
    branch = 0

    CALL POPREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
    CALL POPREALARRAY(dp2, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL POPREALARRAY(frac)
    CALL POPREALARRAY(qn2, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*nq)
    CALL POPREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
    CALL POPREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
    CALL POPREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
    CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
    js = bd%js
    rarea => gridstruct%rarea
    jsd = bd%jsd
    ied = bd%ied
    ie = bd%ie
    isd = bd%isd
    is = bd%is
    je = bd%je
    jed = bd%jed
    mfx_ad = 0.0
    mfy_ad = 0.0
    cx_ad = 0.0
    cy_ad = 0.0
    xfx_ad = 0.0
    dp2_ad = 0.0
    qn2_ad = 0.0
    ra_x_ad = 0.0
    ra_y_ad = 0.0
    yfx_ad = 0.0
    fx_ad = 0.0
    fy_ad = 0.0
    DO k=npz,1,-1
      CALL POPINTEGER(ad_to)
      DO it=ad_to,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          CALL POPREALARRAY(qn2, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*nq&
&                     )
          CALL MPP_UPDATE_DOMAINS_ADM(qn2, qn2_ad, domain)
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(dp1(i, j, k))
              dp2_ad(i, j) = dp2_ad(i, j) + dp1_ad(i, j, k)
              dp1_ad(i, j, k) = 0.0
            END DO
          END DO
        END IF
        DO iq=nq,1,-1
          CALL POPCONTROL(2,branch)
          IF (branch .EQ. 0) THEN
            DO j=je,js,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(q(i, j, k, iq))
                temp_ad4 = q_ad(i, j, k, iq)/dp2(i, j)
                temp = q(i, j, k, iq)
                temp_ad5 = rarea(i, j)*temp_ad4
                dp1_ad(i, j, k) = dp1_ad(i, j, k) + temp*temp_ad4
                fx_ad(i, j) = fx_ad(i, j) + temp_ad5
                fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad5
                fy_ad(i, j) = fy_ad(i, j) + temp_ad5
                fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad5
                dp2_ad(i, j) = dp2_ad(i, j) - (temp*dp1(i, j, k)+rarea(i&
&                 , j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*&
&                 temp_ad4/dp2(i, j)
                q_ad(i, j, k, iq) = dp1(i, j, k)*temp_ad4
              END DO
            END DO
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1)&
&                          *(jed-jsd+1))
              CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
              CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
              CALL FV_TP_2D_ADM(q(isd:ied, jsd:jed, k, iq), q_ad(isd:ied&
&                         , jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k), &
&                         cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied, js:je+&
&                         1, k), cy_ad(isd:ied, js:je+1, k), npx, npy, &
&                         hord_pert, fx, fx_ad, fy, fy_ad, xfx(is:ie+1, &
&                         jsd:jed, k), xfx_ad(is:ie+1, jsd:jed, k), yfx(&
&                         isd:ied, js:je+1, k), yfx_ad(isd:ied, js:je+1&
&                         , k), gridstruct, bd, ra_x, ra_x_ad, ra_y, &
&                         ra_y_ad, mfx(is:ie+1, js:je, k), mfx_ad(is:ie+&
&                         1, js:je, k), mfy(is:ie, js:je+1, k), mfy_ad(&
&                         is:ie, js:je+1, k))
            ELSE
              CALL FV_TP_2D_BWD(q(isd:ied, jsd:jed, k, iq), q_ad(isd:&
&                            ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, &
&                            k), cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied&
&                            , js:je+1, k), cy_ad(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fx_ad, fy, fy_ad, xfx(&
&                            is:ie+1, jsd:jed, k), xfx_ad(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), yfx_ad(&
&                            isd:ied, js:je+1, k), gridstruct, bd, ra_x&
&                            , ra_x_ad, ra_y, ra_y_ad, mfx(is:ie+1, js:&
&                            je, k), mfx_ad(is:ie+1, js:je, k), mfy(is:&
&                            ie, js:je+1, k), mfy_ad(is:ie, js:je+1, k))
            END IF
          ELSE
            IF (branch .EQ. 1) THEN
              DO j=je,js,-1
                DO i=ie,is,-1
                  CALL POPREALARRAY(q(i, j, k, iq))
                  temp_ad2 = q_ad(i, j, k, iq)/dp2(i, j)
                  temp_ad3 = rarea(i, j)*temp_ad2
                  qn2_ad(i, j, iq) = qn2_ad(i, j, iq) + dp1(i, j, k)*&
&                   temp_ad2
                  dp1_ad(i, j, k) = dp1_ad(i, j, k) + qn2(i, j, iq)*&
&                   temp_ad2
                  fx_ad(i, j) = fx_ad(i, j) + temp_ad3
                  fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad3
                  fy_ad(i, j) = fy_ad(i, j) + temp_ad3
                  fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad3
                  dp2_ad(i, j) = dp2_ad(i, j) - (qn2(i, j, iq)*dp1(i, j&
&                   , k)+rarea(i, j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i&
&                   , j+1)))*temp_ad2/dp2(i, j)
                  q_ad(i, j, k, iq) = 0.0
                END DO
              END DO
            ELSE
              DO j=je,js,-1
                DO i=ie,is,-1
                  CALL POPREALARRAY(qn2(i, j, iq))
                  temp_ad0 = qn2_ad(i, j, iq)/dp2(i, j)
                  temp_ad1 = rarea(i, j)*temp_ad0
                  dp1_ad(i, j, k) = dp1_ad(i, j, k) + qn2(i, j, iq)*&
&                   temp_ad0
                  fx_ad(i, j) = fx_ad(i, j) + temp_ad1
                  fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad1
                  fy_ad(i, j) = fy_ad(i, j) + temp_ad1
                  fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad1
                  dp2_ad(i, j) = dp2_ad(i, j) - (qn2(i, j, iq)*dp1(i, j&
&                   , k)+rarea(i, j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i&
&                   , j+1)))*temp_ad0/dp2(i, j)
                  qn2_ad(i, j, iq) = dp1(i, j, k)*temp_ad0
                END DO
              END DO
            END IF
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL FV_TP_2D_BWD(qn2(isd:ied, jsd:jed, iq), qn2_ad(isd&
&                            :ied, jsd:jed, iq), cx(is:ie+1, jsd:jed, k)&
&                            , cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied, &
&                            js:je+1, k), cy_ad(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fx_ad, fy, fy_ad, xfx(&
&                            is:ie+1, jsd:jed, k), xfx_ad(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), yfx_ad(&
&                            isd:ied, js:je+1, k), gridstruct, bd, ra_x&
&                            , ra_x_ad, ra_y, ra_y_ad, mfx(is:ie+1, js:&
&                            je, k), mfx_ad(is:ie+1, js:je, k), mfy(is:&
&                            ie, js:je+1, k), mfy_ad(is:ie, js:je+1, k))
            ELSE
              CALL POPREALARRAY(qn2(isd:ied, jsd:jed, iq), (ied-isd+1)*&
&                          (jed-jsd+1))
              CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
              CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
              CALL FV_TP_2D_ADM(qn2(isd:ied, jsd:jed, iq), qn2_ad(isd:&
&                         ied, jsd:jed, iq), cx(is:ie+1, jsd:jed, k), &
&                         cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied, js:je+&
&                         1, k), cy_ad(isd:ied, js:je+1, k), npx, npy, &
&                         hord_pert, fx, fx_ad, fy, fy_ad, xfx(is:ie+1, &
&                         jsd:jed, k), xfx_ad(is:ie+1, jsd:jed, k), yfx(&
&                         isd:ied, js:je+1, k), yfx_ad(isd:ied, js:je+1&
&                         , k), gridstruct, bd, ra_x, ra_x_ad, ra_y, &
&                         ra_y_ad, mfx(is:ie+1, js:je, k), mfx_ad(is:ie+&
&                         1, js:je, k), mfy(is:ie, js:je+1, k), mfy_ad(&
&                         is:ie, js:je+1, k))
            END IF
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              DO j=jed,jsd,-1
                DO i=ied,isd,-1
                  CALL POPREALARRAY(qn2(i, j, iq))
                  q_ad(i, j, k, iq) = q_ad(i, j, k, iq) + qn2_ad(i, j, &
&                   iq)
                  qn2_ad(i, j, iq) = 0.0
                END DO
              END DO
            END IF
          END IF
        END DO
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(dp2(i, j))
            temp_ad = rarea(i, j)*dp2_ad(i, j)
            dp1_ad(i, j, k) = dp1_ad(i, j, k) + dp2_ad(i, j)
            mfx_ad(i, j, k) = mfx_ad(i, j, k) + temp_ad
            mfx_ad(i+1, j, k) = mfx_ad(i+1, j, k) - temp_ad
            mfy_ad(i, j, k) = mfy_ad(i, j, k) + temp_ad
            mfy_ad(i, j+1, k) = mfy_ad(i, j+1, k) - temp_ad
            dp2_ad(i, j) = 0.0
          END DO
        END DO
      END DO
      DO j=jed,jsd,-1
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          DO i=ied,isd,-1
            CALL POPREALARRAY(ra_y(i, j))
            yfx_ad(i, j, k) = yfx_ad(i, j, k) + ra_y_ad(i, j)
            yfx_ad(i, j+1, k) = yfx_ad(i, j+1, k) - ra_y_ad(i, j)
            ra_y_ad(i, j) = 0.0
          END DO
        END IF
        DO i=ie,is,-1
          CALL POPREALARRAY(ra_x(i, j))
          xfx_ad(i, j, k) = xfx_ad(i, j, k) + ra_x_ad(i, j)
          xfx_ad(i+1, j, k) = xfx_ad(i+1, j, k) - ra_x_ad(i, j)
          ra_x_ad(i, j) = 0.0
        END DO
      END DO
    END DO
    DO k=npz,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        DO j=je+1,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(mfy(i, j, k))
            mfy_ad(i, j, k) = frac*mfy_ad(i, j, k)
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ied,isd,-1
            yfx_ad(i, j, k) = frac*yfx_ad(i, j, k)
            CALL POPREALARRAY(cy(i, j, k))
            cy_ad(i, j, k) = frac*cy_ad(i, j, k)
          END DO
        END DO
        DO j=je,js,-1
          DO i=ie+1,is,-1
            CALL POPREALARRAY(mfx(i, j, k))
            mfx_ad(i, j, k) = frac*mfx_ad(i, j, k)
          END DO
        END DO
        DO j=jed,jsd,-1
          DO i=ie+1,is,-1
            xfx_ad(i, j, k) = frac*xfx_ad(i, j, k)
            CALL POPREALARRAY(cx(i, j, k))
            cx_ad(i, j, k) = frac*cx_ad(i, j, k)
          END DO
        END DO
        CALL POPREALARRAY(frac)
      END IF
    END DO
    dxa => gridstruct%dxa
    dx => gridstruct%dx
    dy => gridstruct%dy
    sin_sg => gridstruct%sin_sg
    dya => gridstruct%dya
    DO k=npz,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPCONTROL(1,branch)
          END DO
        END DO
      ELSE
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPCONTROL(2,branch)
          END DO
        END DO
      END IF
      DO j=je+1,js,-1
        DO i=ied,isd,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            cy_ad(i, j, k) = cy_ad(i, j, k) + dya(i, j)*dx(i, j)*sin_sg(&
&             i, j, 2)*yfx_ad(i, j, k)
            yfx_ad(i, j, k) = 0.0
          ELSE
            cy_ad(i, j, k) = cy_ad(i, j, k) + dya(i, j-1)*dx(i, j)*&
&             sin_sg(i, j-1, 4)*yfx_ad(i, j, k)
            yfx_ad(i, j, k) = 0.0
          END IF
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            cx_ad(i, j, k) = cx_ad(i, j, k) + dxa(i, j)*dy(i, j)*sin_sg(&
&             i, j, 1)*xfx_ad(i, j, k)
            xfx_ad(i, j, k) = 0.0
          ELSE
            cx_ad(i, j, k) = cx_ad(i, j, k) + dxa(i-1, j)*dy(i, j)*&
&             sin_sg(i-1, j, 3)*xfx_ad(i, j, k)
            xfx_ad(i, j, k) = 0.0
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE TRACER_2D_1L_BWD
!-----------------------------------------------------------------------
! !ROUTINE: Perform 2D horizontal-to-lagrangian transport
!-----------------------------------------------------------------------
  SUBROUTINE TRACER_2D_1L(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, &
&   domain, npx, npy, npz, nq, hord, q_split, dt, id_divg, q_pack, &
&   nord_tr, trdm, hord_pert, nord_tr_pert, trdm_pert, split_damp_tr, &
&   dpa)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord, nord_tr
    INTEGER, INTENT(IN) :: hord_pert, nord_tr_pert
    LOGICAL, INTENT(IN) :: split_damp_tr
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt, trdm, trdm_pert
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: q_pack
! Tracers
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
! DELP after advection
    REAL, OPTIONAL, INTENT(OUT) :: dpa(bd%is:bd%ie, bd%js:bd%je)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Local Arrays
! 3D tracers
    REAL :: qn2(bd%isd:bd%ied, bd%jsd:bd%jed, nq)
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cmax(npz)
    REAL :: frac
    INTEGER :: nsplt
    INTEGER :: i, j, k, it, iq
    REAL, DIMENSION(:, :), POINTER :: area, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: dxa, dya, dx, dy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    INTRINSIC PRESENT
    REAL :: max1
    REAL :: x1
    REAL :: z1
    REAL :: y3
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
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    dx => gridstruct%dx
    dy => gridstruct%dy
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx,cmax)
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sin_sg(i-1, &
&             j, 3)
          ELSE
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1&
&             )
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sin_sg(i, j-&
&             1, 4)
          ELSE
            yfx(i, j, k) = cy(i, j, k)*dya(i, j)*dx(i, j)*sin_sg(i, j, 2&
&             )
          END IF
        END DO
      END DO
      cmax(k) = 0.
      IF (k .LT. npz/6) THEN
        DO j=js,je
          DO i=is,ie
            IF (cx(i, j, k) .GE. 0.) THEN
              y1 = cx(i, j, k)
            ELSE
              y1 = -cx(i, j, k)
            END IF
            IF (cy(i, j, k) .GE. 0.) THEN
              z1 = cy(i, j, k)
            ELSE
              z1 = -cy(i, j, k)
            END IF
            IF (cmax(k) .LT. y1) THEN
              IF (y1 .LT. z1) THEN
                cmax(k) = z1
              ELSE
                cmax(k) = y1
              END IF
            ELSE IF (cmax(k) .LT. z1) THEN
              cmax(k) = z1
            ELSE
              cmax(k) = cmax(k)
            END IF
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            IF (cx(i, j, k) .GE. 0.) THEN
              x1 = cx(i, j, k)
            ELSE
              x1 = -cx(i, j, k)
            END IF
            IF (cy(i, j, k) .GE. 0.) THEN
              y3 = cy(i, j, k)
            ELSE
              y3 = -cy(i, j, k)
            END IF
            IF (x1 .LT. y3) THEN
              max1 = y3
            ELSE
              max1 = x1
            END IF
            y2 = max1 + 1. - sin_sg(i, j, 5)
            IF (cmax(k) .LT. y2) THEN
              cmax(k) = y2
            ELSE
              cmax(k) = cmax(k)
            END IF
          END DO
        END DO
      END IF
    END DO
! k-loop
    CALL MP_REDUCE_MAX(cmax, npz)
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx, &
!$OMP                                  cy,yfx,mfx,mfy,cmax)   &
!$OMP                          private(nsplt, frac)
    DO k=1,npz
      nsplt = INT(1. + cmax(k))
      IF (nsplt .GT. 1) THEN
        frac = 1./REAL(nsplt)
        DO j=jsd,jed
          DO i=is,ie+1
            cx(i, j, k) = cx(i, j, k)*frac
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            cy(i, j, k) = cy(i, j, k)*frac
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END IF
    END DO
    CALL TIMING_ON('COMM_TOTAL')
    CALL TIMING_ON('COMM_TRACER')
    CALL COMPLETE_GROUP_HALO_UPDATE(q_pack, domain)
    CALL TIMING_OFF('COMM_TRACER')
    CALL TIMING_OFF('COMM_TOTAL')
! Begin k-independent tracer transport; can not be OpenMPed because the mpp_update call.
    DO k=1,npz
!$OMP parallel do default(none) shared(k,is,ie,js,je,isd,ied,jsd,jed,xfx,area,yfx,ra_x,ra_y)
      DO j=jsd,jed
        DO i=is,ie
          ra_x(i, j) = area(i, j) + (xfx(i, j, k)-xfx(i+1, j, k))
        END DO
        IF (j .GE. js .AND. j .LE. je) THEN
          DO i=isd,ied
            ra_y(i, j) = area(i, j) + (yfx(i, j, k)-yfx(i, j+1, k))
          END DO
        END IF
      END DO
      nsplt = INT(1. + cmax(k))
      DO it=1,nsplt
!$OMP parallel do default(none) shared(k,is,ie,js,je,rarea,mfx,mfy,dp1,dp2)
        DO j=js,je
          DO i=is,ie
            dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+(mfy&
&             (i, j, k)-mfy(i, j+1, k)))*rarea(i, j)
          END DO
        END DO
!$OMP parallel do default(none) shared(k,nsplt,it,is,ie,js,je,isd,ied,jsd,jed,npx,npy,cx,xfx,hord,trdm, &
!$OMP                                  nord_tr,nq,gridstruct,bd,cy,yfx,mfx,mfy,qn2,q,ra_x,ra_y,dp1,dp2,rarea) &
!$OMP                          private(fx,fy)
        DO iq=1,nq
          IF (nsplt .NE. 1) THEN
            IF (it .EQ. 1) THEN
              DO j=jsd,jed
                DO i=isd,ied
                  qn2(i, j, iq) = q(i, j, k, iq)
                END DO
              END DO
            END IF
            IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D(qn2(isd:ied, jsd:jed, iq), cx(is:ie+1, &
&                        jsd:jed, k), cy(isd:ied, js:je+1, k), npx, npy&
&                        , hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                        isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                        ra_y, mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie&
&                        , js:je+1, k))
            ELSE
              CALL FV_TP_2D(qn2(isd:ied, jsd:jed, iq), cx(is:ie+1, jsd:&
&                     jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                     hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                     isd:ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, &
&                     mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1&
&                     , k))
            END IF
            IF (it .LT. nsplt) THEN
! not last call
              DO j=js,je
                DO i=is,ie
                  qn2(i, j, iq) = (qn2(i, j, iq)*dp1(i, j, k)+(fx(i, j)-&
&                   fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i&
&                   , j)
                END DO
              END DO
            ELSE
              DO j=js,je
                DO i=is,ie
                  q(i, j, k, iq) = (qn2(i, j, iq)*dp1(i, j, k)+(fx(i, j)&
&                   -fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(&
&                   i, j)
                END DO
              END DO
            END IF
          ELSE
            IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, &
&                        jsd:jed, k), cy(isd:ied, js:je+1, k), npx, npy&
&                        , hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                        isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                        ra_y, mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie&
&                        , js:je+1, k))
            ELSE
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd:&
&                     jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                     hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                     isd:ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, &
&                     mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1&
&                     , k))
            END IF
            DO j=js,je
              DO i=is,ie
                q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-&
&                 fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i, &
&                 j)
              END DO
            END DO
          END IF
        END DO
!  tracer-loop
        IF (it .LT. nsplt) THEN
! not last call
          DO j=js,je
            DO i=is,ie
              dp1(i, j, k) = dp2(i, j)
            END DO
          END DO
          CALL TIMING_ON('COMM_TOTAL')
          CALL TIMING_ON('COMM_TRACER')
          CALL MPP_UPDATE_DOMAINS(qn2, domain)
          CALL TIMING_OFF('COMM_TRACER')
          CALL TIMING_OFF('COMM_TOTAL')
        END IF
      END DO
    END DO
! time-split loop
! k-loop
    IF (PRESENT(dpa)) dpa = dp2
  END SUBROUTINE TRACER_2D_1L
!  Differentiation of tracer_2d in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
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
!   gradient     of useful results: q dp1
!   with respect to varying inputs: q dp1 mfx mfy cx cy
  SUBROUTINE TRACER_2D_FWD(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, &
&   domain, npx, npy, npz, nq, hord, q_split, dt, id_divg, q_pack, &
&   nord_tr, trdm, hord_pert, nord_tr_pert, trdm_pert, split_damp_tr, &
&   dpa)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord, nord_tr
    INTEGER, INTENT(IN) :: hord_pert, nord_tr_pert
    LOGICAL, INTENT(IN) :: split_damp_tr
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt, trdm, trdm_pert
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: q_pack
! Tracers
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
! DELP after advection
    REAL, OPTIONAL, INTENT(OUT) :: dpa(bd%is:bd%ie, bd%js:bd%je, npz)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Local Arrays
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cmax(npz)
    REAL :: c_global
    REAL :: frac, rdt
    INTEGER :: ksplt(npz)
    INTEGER :: nsplt
    INTEGER :: i, j, k, it, iq
    REAL, DIMENSION(:, :), POINTER :: area, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: dxa, dya, dx, dy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    INTRINSIC PRESENT
    REAL :: max1
    LOGICAL :: res
    REAL :: x1
    REAL :: z1
    REAL :: y3
    REAL :: y2
    REAL :: y1

    dp2 = 0.0
    fx = 0.0
    fy = 0.0
    ra_x = 0.0
    ra_y = 0.0
    xfx = 0.0
    yfx = 0.0
    cmax = 0.0
    c_global = 0.0
    frac = 0.0
    rdt = 0.0
    ksplt = 0
    nsplt = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    max1 = 0.0
    x1 = 0.0
    z1 = 0.0
    y3 = 0.0
    y2 = 0.0
    y1 = 0.0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    dx => gridstruct%dx
    dy => gridstruct%dy
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx,cmax,q_split,ksplt)
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sin_sg(i-1, &
&             j, 3)
            CALL PUSHCONTROL(1,1)
          ELSE
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1&
&             )
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sin_sg(i, j-&
&             1, 4)
            CALL PUSHCONTROL(1,1)
          ELSE
            yfx(i, j, k) = cy(i, j, k)*dya(i, j)*dx(i, j)*sin_sg(i, j, 2&
&             )
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      IF (q_split .EQ. 0) THEN
        cmax(k) = 0.
        IF (k .LT. npz/6) THEN
          DO j=js,je
            DO i=is,ie
              IF (cx(i, j, k) .GE. 0.) THEN
                y1 = cx(i, j, k)
              ELSE
                y1 = -cx(i, j, k)
              END IF
              IF (cy(i, j, k) .GE. 0.) THEN
                z1 = cy(i, j, k)
              ELSE
                z1 = -cy(i, j, k)
              END IF
              IF (cmax(k) .LT. y1) THEN
                IF (y1 .LT. z1) THEN
                  CALL PUSHCONTROL(2,0)
                  cmax(k) = z1
                ELSE
                  CALL PUSHCONTROL(2,1)
                  cmax(k) = y1
                END IF
              ELSE IF (cmax(k) .LT. z1) THEN
                CALL PUSHCONTROL(2,2)
                cmax(k) = z1
              ELSE
                CALL PUSHCONTROL(2,3)
                cmax(k) = cmax(k)
              END IF
            END DO
          END DO
          CALL PUSHCONTROL(2,0)
        ELSE
          DO j=js,je
            DO i=is,ie
              IF (cx(i, j, k) .GE. 0.) THEN
                x1 = cx(i, j, k)
              ELSE
                x1 = -cx(i, j, k)
              END IF
              IF (cy(i, j, k) .GE. 0.) THEN
                y3 = cy(i, j, k)
              ELSE
                y3 = -cy(i, j, k)
              END IF
              IF (x1 .LT. y3) THEN
                max1 = y3
              ELSE
                max1 = x1
              END IF
              y2 = max1 + 1. - sin_sg(i, j, 5)
              IF (cmax(k) .LT. y2) THEN
                CALL PUSHCONTROL(1,0)
                cmax(k) = y2
              ELSE
                CALL PUSHCONTROL(1,1)
                cmax(k) = cmax(k)
              END IF
            END DO
          END DO
          CALL PUSHCONTROL(2,1)
        END IF
      ELSE
        CALL PUSHCONTROL(2,2)
      END IF
      ksplt(k) = 1
    END DO
!--------------------------------------------------------------------------------
! Determine global nsplt:
    IF (q_split .EQ. 0) THEN
      CALL MP_REDUCE_MAX(cmax, npz)
! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      IF (npz .NE. 1) THEN
! if NOT shallow water test case
        DO k=2,npz
          IF (cmax(k) .LT. c_global) THEN
            CALL PUSHCONTROL(1,0)
            c_global = c_global
          ELSE
            CALL PUSHCONTROL(1,1)
            c_global = cmax(k)
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      nsplt = INT(1. + c_global)
      res = IS_MASTER()
      IF (res .AND. nsplt .GT. 4) THEN
        CALL PUSHCONTROL(1,0)
        WRITE(*, *) 'Tracer_2d_split=', nsplt, c_global
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    ELSE
      CALL PUSHCONTROL(1,1)
      nsplt = q_split
    END IF
!--------------------------------------------------------------------------------
    IF (nsplt .NE. 1) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,mfx,cy,yfx,mfy,cmax,nsplt,ksplt) &
!$OMP                          private( frac )
      DO k=1,npz
        ksplt(k) = INT(1. + cmax(k))
        CALL PUSHREALARRAY(frac)
        frac = 1./REAL(ksplt(k))
        DO j=jsd,jed
          DO i=is,ie+1
            CALL PUSHREALARRAY(cx(i, j, k))
            cx(i, j, k) = cx(i, j, k)*frac
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            CALL PUSHREALARRAY(mfx(i, j, k))
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            CALL PUSHREALARRAY(cy(i, j, k))
            cy(i, j, k) = cy(i, j, k)*frac
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            CALL PUSHREALARRAY(mfy(i, j, k))
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
    DO it=1,nsplt
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx,mfy,rarea,nq,ksplt,&
!$OMP                                  area,xfx,yfx,q,cx,cy,npx,npy,hord,gridstruct,bd,it,nsplt,nord_tr,trdm) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy)
      DO k=1,npz
! ksplt
        IF (it .LE. ksplt(k)) THEN
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(dp2(i, j))
              dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+(&
&               mfy(i, j, k)-mfy(i, j+1, k)))*rarea(i, j)
            END DO
          END DO
          DO j=jsd,jed
            DO i=is,ie
              CALL PUSHREALARRAY(ra_x(i, j))
              ra_x(i, j) = area(i, j) + (xfx(i, j, k)-xfx(i+1, j, k))
            END DO
          END DO
          DO j=js,je
            DO i=isd,ied
              CALL PUSHREALARRAY(ra_y(i, j))
              ra_y(i, j) = area(i, j) + (yfx(i, j, k)-yfx(i, j+1, k))
            END DO
          END DO
          DO iq=1,nq
            IF (it .EQ. 1 .AND. trdm .GT. 1.e-4) THEN
              IF (hord .EQ. hord_pert) THEN
                CALL FV_TP_2D_FWD(q(isd:ied, jsd:jed, k, iq), cx(is:&
&                              ie+1, jsd:jed, k), cy(isd:ied, js:je+1, k&
&                              ), npx, npy, hord, fx, fy, xfx(is:ie+1, &
&                              jsd:jed, k), yfx(isd:ied, js:je+1, k), &
&                              gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie&
&                              +1, js:je, k), mfy=mfy(is:ie, js:je+1, k)&
&                              , mass=dp1(isd:ied, jsd:jed, k), nord=&
&                              nord_tr, damp_c=trdm)
                CALL PUSHCONTROL(2,3)
              ELSE
                CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
                CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
                CALL PUSHREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd&
&                             +1)*(jed-jsd+1))
                CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, &
&                       jsd:jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                       hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx&
&                       (isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                       ra_y, mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie&
&                       , js:je+1, k), mass=dp1(isd:ied, jsd:jed, k), &
&                       nord=nord_tr, damp_c=trdm)
                CALL PUSHCONTROL(2,2)
              END IF
            ELSE IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D_FWD(q(isd:ied, jsd:jed, k, iq), cx(is:ie+&
&                            1, jsd:jed, k), cy(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fy, xfx(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), &
&                            gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1&
&                            , js:je, k), mfy=mfy(is:ie, js:je+1, k))
              CALL PUSHCONTROL(2,1)
            ELSE
              CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
              CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
              CALL PUSHREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1&
&                           )*(jed-jsd+1))
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd:&
&                     jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                     hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                     isd:ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, &
&                     mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1&
&                     , k))
              CALL PUSHCONTROL(2,0)
            END IF
            DO j=js,je
              DO i=is,ie
                CALL PUSHREALARRAY(q(i, j, k, iq))
                q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-&
&                 fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i, &
&                 j)
              END DO
            END DO
          END DO
          IF (it .NE. nsplt) THEN
            DO j=js,je
              DO i=is,ie
                CALL PUSHREALARRAY(dp1(i, j, k))
                dp1(i, j, k) = dp2(i, j)
              END DO
            END DO
            CALL PUSHCONTROL(2,2)
          ELSE
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE
          CALL PUSHCONTROL(2,0)
        END IF
      END DO
! npz
      IF (it .NE. nsplt) THEN
        CALL PUSHREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz*&
&                     nq)
        CALL START_GROUP_HALO_UPDATE(q_pack, q, domain)
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
    CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
    CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
    CALL PUSHINTEGER(nsplt)
    CALL PUSHREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    CALL PUSHREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
    CALL PUSHREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
    CALL PUSHREALARRAY(frac)
    CALL PUSHREALARRAY(dp2, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL PUSHREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
  END SUBROUTINE TRACER_2D_FWD
!  Differentiation of tracer_2d in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_m
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
!   gradient     of useful results: q dp1
!   with respect to varying inputs: q dp1 mfx mfy cx cy
  SUBROUTINE TRACER_2D_BWD(q, q_ad, dp1, dp1_ad, mfx, mfx_ad, mfy, &
&   mfy_ad, cx, cx_ad, cy, cy_ad, gridstruct, bd, domain, npx, npy, npz&
&   , nq, hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, hord_pert, &
&   nord_tr_pert, trdm_pert, split_damp_tr, dpa)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord, nord_tr
    INTEGER, INTENT(IN) :: hord_pert, nord_tr_pert
    LOGICAL, INTENT(IN) :: split_damp_tr
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt, trdm, trdm_pert
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: q_pack
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dp1_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_ad(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_ad(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, OPTIONAL, INTENT(OUT) :: dpa(bd%is:bd%ie, bd%js:bd%je, npz)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: dp2_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_ad(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_ad(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_ad(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_ad(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cmax(npz)
    REAL :: c_global
    REAL :: frac, rdt
    INTEGER :: ksplt(npz)
    INTEGER :: nsplt
    INTEGER :: i, j, k, it, iq
    REAL, DIMENSION(:, :), POINTER :: area, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: dxa, dya, dx, dy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    INTRINSIC PRESENT
    REAL :: max1
    REAL :: temp_ad
    REAL :: temp
    REAL :: temp_ad0
    REAL :: temp_ad1
    INTEGER :: branch
    REAL :: x1
    REAL :: z1
    REAL :: y3
    REAL :: y2
    REAL :: y1

    dp2 = 0.0
    fx = 0.0
    fy = 0.0
    ra_x = 0.0
    ra_y = 0.0
    xfx = 0.0
    yfx = 0.0
    cmax = 0.0
    c_global = 0.0
    frac = 0.0
    rdt = 0.0
    ksplt = 0
    nsplt = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    max1 = 0.0
    x1 = 0.0
    z1 = 0.0
    y3 = 0.0
    y2 = 0.0
    y1 = 0.0
    branch = 0

    CALL POPREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
    CALL POPREALARRAY(dp2, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
    CALL POPREALARRAY(frac)
    CALL POPREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
    CALL POPREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
    CALL POPREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
    CALL POPINTEGER(nsplt)
    CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
    CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
    js = bd%js
    rarea => gridstruct%rarea
    jsd = bd%jsd
    ied = bd%ied
    ie = bd%ie
    isd = bd%isd
    is = bd%is
    je = bd%je
    jed = bd%jed
    mfx_ad = 0.0
    mfy_ad = 0.0
    cx_ad = 0.0
    cy_ad = 0.0
    xfx_ad = 0.0
    dp2_ad = 0.0
    ra_x_ad = 0.0
    ra_y_ad = 0.0
    yfx_ad = 0.0
    fx_ad = 0.0
    fy_ad = 0.0
    DO it=nsplt,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
        CALL POPREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz*nq&
&                   )
        CALL START_GROUP_HALO_UPDATE_ADM(q_pack, q, q_ad, domain)
      END IF
      DO k=npz,1,-1
        CALL POPCONTROL(2,branch)
        IF (branch .NE. 0) THEN
          IF (branch .NE. 1) THEN
            DO j=je,js,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(dp1(i, j, k))
                dp2_ad(i, j) = dp2_ad(i, j) + dp1_ad(i, j, k)
                dp1_ad(i, j, k) = 0.0
              END DO
            END DO
          END IF
          DO iq=nq,1,-1
            DO j=je,js,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(q(i, j, k, iq))
                temp_ad0 = q_ad(i, j, k, iq)/dp2(i, j)
                temp = q(i, j, k, iq)
                temp_ad1 = rarea(i, j)*temp_ad0
                dp1_ad(i, j, k) = dp1_ad(i, j, k) + temp*temp_ad0
                fx_ad(i, j) = fx_ad(i, j) + temp_ad1
                fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad1
                fy_ad(i, j) = fy_ad(i, j) + temp_ad1
                fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad1
                dp2_ad(i, j) = dp2_ad(i, j) - (temp*dp1(i, j, k)+rarea(i&
&                 , j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*&
&                 temp_ad0/dp2(i, j)
                q_ad(i, j, k, iq) = dp1(i, j, k)*temp_ad0
              END DO
            END DO
            CALL POPCONTROL(2,branch)
            IF (branch .LT. 2) THEN
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+&
&                            1)*(jed-jsd+1))
                CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
                CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
                CALL FV_TP_2D_ADM(q(isd:ied, jsd:jed, k, iq), q_ad(isd:&
&                           ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k&
&                           ), cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied, &
&                           js:je+1, k), cy_ad(isd:ied, js:je+1, k), npx&
&                           , npy, hord_pert, fx, fx_ad, fy, fy_ad, xfx(&
&                           is:ie+1, jsd:jed, k), xfx_ad(is:ie+1, jsd:&
&                           jed, k), yfx(isd:ied, js:je+1, k), yfx_ad(&
&                           isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                           ra_x_ad, ra_y, ra_y_ad, mfx(is:ie+1, js:je, &
&                           k), mfx_ad(is:ie+1, js:je, k), mfy(is:ie, js&
&                           :je+1, k), mfy_ad(is:ie, js:je+1, k))
              ELSE
                CALL FV_TP_2D_BWD(q(isd:ied, jsd:jed, k, iq), q_ad(&
&                              isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd&
&                              :jed, k), cx_ad(is:ie+1, jsd:jed, k), cy(&
&                              isd:ied, js:je+1, k), cy_ad(isd:ied, js:&
&                              je+1, k), npx, npy, hord, fx, fx_ad, fy, &
&                              fy_ad, xfx(is:ie+1, jsd:jed, k), xfx_ad(&
&                              is:ie+1, jsd:jed, k), yfx(isd:ied, js:je+&
&                              1, k), yfx_ad(isd:ied, js:je+1, k), &
&                              gridstruct, bd, ra_x, ra_x_ad, ra_y, &
&                              ra_y_ad, mfx(is:ie+1, js:je, k), mfx_ad(&
&                              is:ie+1, js:je, k), mfy(is:ie, js:je+1, k&
&                              ), mfy_ad(is:ie, js:je+1, k))
              END IF
            ELSE IF (branch .EQ. 2) THEN
              CALL POPREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1)&
&                          *(jed-jsd+1))
              CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
              CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
              CALL FV_TP_2D_ADM(q(isd:ied, jsd:jed, k, iq), q_ad(isd:ied&
&                         , jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k), &
&                         cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied, js:je+&
&                         1, k), cy_ad(isd:ied, js:je+1, k), npx, npy, &
&                         hord_pert, fx, fx_ad, fy, fy_ad, xfx(is:ie+1, &
&                         jsd:jed, k), xfx_ad(is:ie+1, jsd:jed, k), yfx(&
&                         isd:ied, js:je+1, k), yfx_ad(isd:ied, js:je+1&
&                         , k), gridstruct, bd, ra_x, ra_x_ad, ra_y, &
&                         ra_y_ad, mfx(is:ie+1, js:je, k), mfx_ad(is:ie+&
&                         1, js:je, k), mfy(is:ie, js:je+1, k), mfy_ad(&
&                         is:ie, js:je+1, k), dp1(isd:ied, jsd:jed, k), &
&                         dp1_ad(isd:ied, jsd:jed, k), nord=nord_tr_pert&
&                         , damp_c=trdm_pert)
            ELSE
              CALL FV_TP_2D_BWD(q(isd:ied, jsd:jed, k, iq), q_ad(isd:&
&                            ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, &
&                            k), cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied&
&                            , js:je+1, k), cy_ad(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fx_ad, fy, fy_ad, xfx(&
&                            is:ie+1, jsd:jed, k), xfx_ad(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), yfx_ad(&
&                            isd:ied, js:je+1, k), gridstruct, bd, ra_x&
&                            , ra_x_ad, ra_y, ra_y_ad, mfx(is:ie+1, js:&
&                            je, k), mfx_ad(is:ie+1, js:je, k), mfy(is:&
&                            ie, js:je+1, k), mfy_ad(is:ie, js:je+1, k)&
&                            , dp1(isd:ied, jsd:jed, k), dp1_ad(isd:ied&
&                            , jsd:jed, k), nord=nord_tr, damp_c=trdm)
            END IF
          END DO
          DO j=je,js,-1
            DO i=ied,isd,-1
              CALL POPREALARRAY(ra_y(i, j))
              yfx_ad(i, j, k) = yfx_ad(i, j, k) + ra_y_ad(i, j)
              yfx_ad(i, j+1, k) = yfx_ad(i, j+1, k) - ra_y_ad(i, j)
              ra_y_ad(i, j) = 0.0
            END DO
          END DO
          DO j=jed,jsd,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(ra_x(i, j))
              xfx_ad(i, j, k) = xfx_ad(i, j, k) + ra_x_ad(i, j)
              xfx_ad(i+1, j, k) = xfx_ad(i+1, j, k) - ra_x_ad(i, j)
              ra_x_ad(i, j) = 0.0
            END DO
          END DO
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(dp2(i, j))
              temp_ad = rarea(i, j)*dp2_ad(i, j)
              dp1_ad(i, j, k) = dp1_ad(i, j, k) + dp2_ad(i, j)
              mfx_ad(i, j, k) = mfx_ad(i, j, k) + temp_ad
              mfx_ad(i+1, j, k) = mfx_ad(i+1, j, k) - temp_ad
              mfy_ad(i, j, k) = mfy_ad(i, j, k) + temp_ad
              mfy_ad(i, j+1, k) = mfy_ad(i, j+1, k) - temp_ad
              dp2_ad(i, j) = 0.0
            END DO
          END DO
        END IF
      END DO
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      DO k=npz,1,-1
        DO j=je+1,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(mfy(i, j, k))
            mfy_ad(i, j, k) = frac*mfy_ad(i, j, k)
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ied,isd,-1
            yfx_ad(i, j, k) = frac*yfx_ad(i, j, k)
            CALL POPREALARRAY(cy(i, j, k))
            cy_ad(i, j, k) = frac*cy_ad(i, j, k)
          END DO
        END DO
        DO j=je,js,-1
          DO i=ie+1,is,-1
            CALL POPREALARRAY(mfx(i, j, k))
            mfx_ad(i, j, k) = frac*mfx_ad(i, j, k)
          END DO
        END DO
        DO j=jed,jsd,-1
          DO i=ie+1,is,-1
            xfx_ad(i, j, k) = frac*xfx_ad(i, j, k)
            CALL POPREALARRAY(cx(i, j, k))
            cx_ad(i, j, k) = frac*cx_ad(i, j, k)
          END DO
        END DO
        CALL POPREALARRAY(frac)
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO k=npz,2,-1
          CALL POPCONTROL(1,branch)
        END DO
      END IF
    END IF
    dxa => gridstruct%dxa
    dx => gridstruct%dx
    dy => gridstruct%dy
    sin_sg => gridstruct%sin_sg
    dya => gridstruct%dya
    DO k=npz,1,-1
      CALL POPCONTROL(2,branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPCONTROL(2,branch)
          END DO
        END DO
      ELSE IF (branch .EQ. 1) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPCONTROL(1,branch)
          END DO
        END DO
      END IF
      DO j=je+1,js,-1
        DO i=ied,isd,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            cy_ad(i, j, k) = cy_ad(i, j, k) + dya(i, j)*dx(i, j)*sin_sg(&
&             i, j, 2)*yfx_ad(i, j, k)
            yfx_ad(i, j, k) = 0.0
          ELSE
            cy_ad(i, j, k) = cy_ad(i, j, k) + dya(i, j-1)*dx(i, j)*&
&             sin_sg(i, j-1, 4)*yfx_ad(i, j, k)
            yfx_ad(i, j, k) = 0.0
          END IF
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            cx_ad(i, j, k) = cx_ad(i, j, k) + dxa(i, j)*dy(i, j)*sin_sg(&
&             i, j, 1)*xfx_ad(i, j, k)
            xfx_ad(i, j, k) = 0.0
          ELSE
            cx_ad(i, j, k) = cx_ad(i, j, k) + dxa(i-1, j)*dy(i, j)*&
&             sin_sg(i-1, j, 3)*xfx_ad(i, j, k)
            xfx_ad(i, j, k) = 0.0
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE TRACER_2D_BWD
  SUBROUTINE TRACER_2D(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain&
&   , npx, npy, npz, nq, hord, q_split, dt, id_divg, q_pack, nord_tr, &
&   trdm, hord_pert, nord_tr_pert, trdm_pert, split_damp_tr, dpa)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord, nord_tr
    INTEGER, INTENT(IN) :: hord_pert, nord_tr_pert
    LOGICAL, INTENT(IN) :: split_damp_tr
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt, trdm, trdm_pert
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: q_pack
! Tracers
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
! DELP after advection
    REAL, OPTIONAL, INTENT(OUT) :: dpa(bd%is:bd%ie, bd%js:bd%je, npz)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Local Arrays
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cmax(npz)
    REAL :: c_global
    REAL :: frac, rdt
    INTEGER :: ksplt(npz)
    INTEGER :: nsplt
    INTEGER :: i, j, k, it, iq
    REAL, DIMENSION(:, :), POINTER :: area, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: dxa, dya, dx, dy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    INTRINSIC PRESENT
    REAL :: max1
    REAL :: x1
    REAL :: z1
    REAL :: y3
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
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    dx => gridstruct%dx
    dy => gridstruct%dy
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx,cmax,q_split,ksplt)
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sin_sg(i-1, &
&             j, 3)
          ELSE
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1&
&             )
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sin_sg(i, j-&
&             1, 4)
          ELSE
            yfx(i, j, k) = cy(i, j, k)*dya(i, j)*dx(i, j)*sin_sg(i, j, 2&
&             )
          END IF
        END DO
      END DO
      IF (q_split .EQ. 0) THEN
        cmax(k) = 0.
        IF (k .LT. npz/6) THEN
          DO j=js,je
            DO i=is,ie
              IF (cx(i, j, k) .GE. 0.) THEN
                y1 = cx(i, j, k)
              ELSE
                y1 = -cx(i, j, k)
              END IF
              IF (cy(i, j, k) .GE. 0.) THEN
                z1 = cy(i, j, k)
              ELSE
                z1 = -cy(i, j, k)
              END IF
              IF (cmax(k) .LT. y1) THEN
                IF (y1 .LT. z1) THEN
                  cmax(k) = z1
                ELSE
                  cmax(k) = y1
                END IF
              ELSE IF (cmax(k) .LT. z1) THEN
                cmax(k) = z1
              ELSE
                cmax(k) = cmax(k)
              END IF
            END DO
          END DO
        ELSE
          DO j=js,je
            DO i=is,ie
              IF (cx(i, j, k) .GE. 0.) THEN
                x1 = cx(i, j, k)
              ELSE
                x1 = -cx(i, j, k)
              END IF
              IF (cy(i, j, k) .GE. 0.) THEN
                y3 = cy(i, j, k)
              ELSE
                y3 = -cy(i, j, k)
              END IF
              IF (x1 .LT. y3) THEN
                max1 = y3
              ELSE
                max1 = x1
              END IF
              y2 = max1 + 1. - sin_sg(i, j, 5)
              IF (cmax(k) .LT. y2) THEN
                cmax(k) = y2
              ELSE
                cmax(k) = cmax(k)
              END IF
            END DO
          END DO
        END IF
      END IF
      ksplt(k) = 1
    END DO
!--------------------------------------------------------------------------------
! Determine global nsplt:
    IF (q_split .EQ. 0) THEN
      CALL MP_REDUCE_MAX(cmax, npz)
! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      IF (npz .NE. 1) THEN
! if NOT shallow water test case
        DO k=2,npz
          IF (cmax(k) .LT. c_global) THEN
            c_global = c_global
          ELSE
            c_global = cmax(k)
          END IF
        END DO
      END IF
      nsplt = INT(1. + c_global)
      IF (IS_MASTER() .AND. nsplt .GT. 4) WRITE(*, *) 'Tracer_2d_split='&
&                                         , nsplt, c_global
    ELSE
      nsplt = q_split
    END IF
!--------------------------------------------------------------------------------
    IF (nsplt .NE. 1) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,mfx,cy,yfx,mfy,cmax,nsplt,ksplt) &
!$OMP                          private( frac )
      DO k=1,npz
        ksplt(k) = INT(1. + cmax(k))
        frac = 1./REAL(ksplt(k))
        DO j=jsd,jed
          DO i=is,ie+1
            cx(i, j, k) = cx(i, j, k)*frac
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            cy(i, j, k) = cy(i, j, k)*frac
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END DO
    END IF
    DO it=1,nsplt
      CALL TIMING_ON('COMM_TOTAL')
      CALL TIMING_ON('COMM_TRACER')
      CALL COMPLETE_GROUP_HALO_UPDATE(q_pack, domain)
      CALL TIMING_OFF('COMM_TRACER')
      CALL TIMING_OFF('COMM_TOTAL')
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx,mfy,rarea,nq,ksplt,&
!$OMP                                  area,xfx,yfx,q,cx,cy,npx,npy,hord,gridstruct,bd,it,nsplt,nord_tr,trdm) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy)
      DO k=1,npz
! ksplt
        IF (it .LE. ksplt(k)) THEN
          DO j=js,je
            DO i=is,ie
              dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+(&
&               mfy(i, j, k)-mfy(i, j+1, k)))*rarea(i, j)
            END DO
          END DO
          DO j=jsd,jed
            DO i=is,ie
              ra_x(i, j) = area(i, j) + (xfx(i, j, k)-xfx(i+1, j, k))
            END DO
          END DO
          DO j=js,je
            DO i=isd,ied
              ra_y(i, j) = area(i, j) + (yfx(i, j, k)-yfx(i, j+1, k))
            END DO
          END DO
          DO iq=1,nq
            IF (it .EQ. 1 .AND. trdm .GT. 1.e-4) THEN
              IF (hord .EQ. hord_pert) THEN
                CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1&
&                          , jsd:jed, k), cy(isd:ied, js:je+1, k), npx, &
&                          npy, hord, fx, fy, xfx(is:ie+1, jsd:jed, k), &
&                          yfx(isd:ied, js:je+1, k), gridstruct, bd, &
&                          ra_x, ra_y, mfx(is:ie+1, js:je, k), mfy(is:ie&
&                          , js:je+1, k), dp1(isd:ied, jsd:jed, k), &
&                          nord_tr, trdm)
              ELSE
                CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, &
&                       jsd:jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                       hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx&
&                       (isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                       ra_y, mfx(is:ie+1, js:je, k), mfy(is:ie, js:je+1&
&                       , k), dp1(isd:ied, jsd:jed, k), nord_tr, trdm)
              END IF
            ELSE IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, &
&                        jsd:jed, k), cy(isd:ied, js:je+1, k), npx, npy&
&                        , hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                        isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                        ra_y, mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie&
&                        , js:je+1, k))
            ELSE
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd:&
&                     jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                     hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                     isd:ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, &
&                     mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1&
&                     , k))
            END IF
            DO j=js,je
              DO i=is,ie
                q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-&
&                 fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i, &
&                 j)
              END DO
            END DO
          END DO
          IF (it .NE. nsplt) THEN
            DO j=js,je
              DO i=is,ie
                dp1(i, j, k) = dp2(i, j)
              END DO
            END DO
          END IF
        END IF
      END DO
! npz
      IF (it .NE. nsplt) THEN
        CALL TIMING_ON('COMM_TOTAL')
        CALL TIMING_ON('COMM_TRACER')
        CALL START_GROUP_HALO_UPDATE(q_pack, q, domain)
        CALL TIMING_OFF('COMM_TRACER')
        CALL TIMING_OFF('COMM_TOTAL')
      END IF
    END DO
! nsplt
    IF (PRESENT(dpa)) dpa = dp1(bd%is:bd%ie, bd%js:bd%je, 1:npz)
  END SUBROUTINE TRACER_2D
!  Differentiation of tracer_2d_nested in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_
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
!   gradient     of useful results: q dp1
!   with respect to varying inputs: q dp1 mfx mfy cx cy
  SUBROUTINE TRACER_2D_NESTED_FWD(q, dp1, mfx, mfy, cx, cy, gridstruct, &
&   bd, domain, npx, npy, npz, nq, hord, q_split, dt, id_divg, q_pack, &
&   nord_tr, trdm, k_split, neststruct, parent_grid, hord_pert, &
&   nord_tr_pert, trdm_pert, split_damp_tr)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord, nord_tr
    INTEGER, INTENT(IN) :: hord_pert, nord_tr_pert
    LOGICAL, INTENT(IN) :: split_damp_tr
    INTEGER, INTENT(IN) :: q_split, k_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt, trdm, trdm_pert
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: q_pack
! Tracers
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Local Arrays
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cmax(npz)
    REAL :: cmax_t
    REAL :: c_global
    REAL :: frac, rdt
    INTEGER :: nsplt, nsplt_parent
    INTEGER, SAVE :: msg_split_steps=1
    INTEGER :: i, j, k, it, iq
    REAL, DIMENSION(:, :), POINTER :: area, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: dxa, dya, dx, dy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    REAL :: max1
    REAL :: arg1
    REAL :: arg2
    LOGICAL :: res
    REAL :: x2
    REAL :: x1
    REAL :: y2
    REAL :: y1

    dp2 = 0.0
    fx = 0.0
    fy = 0.0
    ra_x = 0.0
    ra_y = 0.0
    xfx = 0.0
    yfx = 0.0
    cmax = 0.0
    cmax_t = 0.0
    c_global = 0.0
    frac = 0.0
    rdt = 0.0
    nsplt = 0
    nsplt_parent = 0
    it = 0
    iq = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    max1 = 0.0
    arg1 = 0.0
    arg2 = 0.0
    x2 = 0.0
    x1 = 0.0
    y2 = 0.0
    y1 = 0.0

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    dx => gridstruct%dx
    dy => gridstruct%dy
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx)
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sin_sg(i-1, &
&             j, 3)
            CALL PUSHCONTROL(1,1)
          ELSE
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1&
&             )
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sin_sg(i, j-&
&             1, 4)
            CALL PUSHCONTROL(1,1)
          ELSE
            yfx(i, j, k) = cy(i, j, k)*dya(i, j)*dx(i, j)*sin_sg(i, j, 2&
&             )
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
    END DO
!--------------------------------------------------------------------------------
    IF (q_split .EQ. 0) THEN
! Determine nsplt
!$OMP parallel do default(none) shared(is,ie,js,je,npz,cmax,cx,cy,sin_sg) &
!$OMP                          private(cmax_t )
      DO k=1,npz
        cmax(k) = 0.
        IF (k .LT. 4) THEN
! Top layers: C < max( abs(c_x), abs(c_y) )
          DO j=js,je
            DO i=is,ie
              IF (cx(i, j, k) .GE. 0.) THEN
                x1 = cx(i, j, k)
              ELSE
                x1 = -cx(i, j, k)
              END IF
              IF (cy(i, j, k) .GE. 0.) THEN
                y1 = cy(i, j, k)
              ELSE
                y1 = -cy(i, j, k)
              END IF
              IF (x1 .LT. y1) THEN
                cmax_t = y1
              ELSE
                cmax_t = x1
              END IF
              IF (cmax_t .LT. cmax(k)) THEN
                CALL PUSHCONTROL(1,0)
                cmax(k) = cmax(k)
              ELSE
                CALL PUSHCONTROL(1,1)
                cmax(k) = cmax_t
              END IF
            END DO
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          DO j=js,je
            DO i=is,ie
              IF (cx(i, j, k) .GE. 0.) THEN
                x2 = cx(i, j, k)
              ELSE
                x2 = -cx(i, j, k)
              END IF
              IF (cy(i, j, k) .GE. 0.) THEN
                y2 = cy(i, j, k)
              ELSE
                y2 = -cy(i, j, k)
              END IF
              IF (x2 .LT. y2) THEN
                max1 = y2
              ELSE
                max1 = x2
              END IF
              cmax_t = max1 + 1. - sin_sg(i, j, 5)
              IF (cmax_t .LT. cmax(k)) THEN
                CALL PUSHCONTROL(1,0)
                cmax(k) = cmax(k)
              ELSE
                CALL PUSHCONTROL(1,1)
                cmax(k) = cmax_t
              END IF
            END DO
          END DO
          CALL PUSHCONTROL(1,0)
        END IF
      END DO
      CALL MP_REDUCE_MAX(cmax, npz)
! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      IF (npz .NE. 1) THEN
! if NOT shallow water test case
        DO k=2,npz
          IF (cmax(k) .LT. c_global) THEN
            CALL PUSHCONTROL(1,0)
            c_global = c_global
          ELSE
            CALL PUSHCONTROL(1,1)
            c_global = cmax(k)
          END IF
        END DO
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
      nsplt = INT(1. + c_global)
      res = IS_MASTER()
      IF (res .AND. nsplt .GT. 3) THEN
        CALL PUSHCONTROL(1,0)
        WRITE(*, *) 'Tracer_2d_split=', nsplt, c_global
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    ELSE
      nsplt = q_split
      CALL PUSHCONTROL(1,1)
    END IF
!--------------------------------------------------------------------------------
    frac = 1./REAL(nsplt)
    IF (nsplt .NE. 1) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,frac,xfx,mfx,cy,yfx,mfy)
      DO k=1,npz
        DO j=jsd,jed
          DO i=is,ie+1
            CALL PUSHREALARRAY(cx(i, j, k))
            cx(i, j, k) = cx(i, j, k)*frac
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            CALL PUSHREALARRAY(mfx(i, j, k))
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            CALL PUSHREALARRAY(cy(i, j, k))
            cy(i, j, k) = cy(i, j, k)*frac
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            CALL PUSHREALARRAY(mfy(i, j, k))
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END DO
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHCONTROL(1,0)
    END IF
    DO it=1,nsplt
      IF (gridstruct%nested) neststruct%tracer_nest_timestep = &
&         neststruct%tracer_nest_timestep + 1
      IF (gridstruct%nested) THEN
        DO iq=1,nq
          arg1 = REAL(neststruct%tracer_nest_timestep) + REAL(nsplt*&
&           k_split)
          arg2 = REAL(nsplt*k_split)
          CALL PUSHREALARRAY(q(isd:ied, jsd:jed, :, iq), (ied-isd+1)*(&
&                       jed-jsd+1)*npz)
          CALL NESTED_GRID_BC_APPLY_INTT(q(isd:ied, jsd:jed, :, iq), 0, &
&                                  0, npx, npy, npz, bd, REAL(neststruct&
&                                  %tracer_nest_timestep) + REAL(nsplt*&
&                                  k_split), REAL(nsplt*k_split), &
&                                  neststruct%q_bc(iq), bctype=&
&                                  neststruct%nestbctype)
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx,mfy,rarea,nq, &
!$OMP                                  area,xfx,yfx,q,cx,cy,npx,npy,hord,gridstruct,bd,it,nsplt,nord_tr,trdm) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(dp2(i, j))
            dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+(mfy&
&             (i, j, k)-mfy(i, j+1, k)))*rarea(i, j)
          END DO
        END DO
        DO j=jsd,jed
          DO i=is,ie
            CALL PUSHREALARRAY(ra_x(i, j))
            ra_x(i, j) = area(i, j) + (xfx(i, j, k)-xfx(i+1, j, k))
          END DO
        END DO
        DO j=js,je
          DO i=isd,ied
            CALL PUSHREALARRAY(ra_y(i, j))
            ra_y(i, j) = area(i, j) + (yfx(i, j, k)-yfx(i, j+1, k))
          END DO
        END DO
        DO iq=1,nq
          IF (it .EQ. 1 .AND. trdm .GT. 1.e-4) THEN
            IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D_FWD(q(isd:ied, jsd:jed, k, iq), cx(is:ie+&
&                            1, jsd:jed, k), cy(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fy, xfx(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), &
&                            gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1&
&                            , js:je, k), mfy=mfy(is:ie, js:je+1, k), &
&                            mass=dp1(isd:ied, jsd:jed, k), nord=nord_tr&
&                            , damp_c=trdm)
              CALL PUSHCONTROL(2,3)
            ELSE
              CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
              CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
              CALL PUSHREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1&
&                           )*(jed-jsd+1))
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd:&
&                     jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                     hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                     isd:ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, &
&                     mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1&
&                     , k), mass=dp1(isd:ied, jsd:jed, k), nord=nord_tr&
&                     , damp_c=trdm)
              CALL PUSHCONTROL(2,2)
            END IF
          ELSE IF (hord .EQ. hord_pert) THEN
            CALL FV_TP_2D_FWD(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1&
&                          , jsd:jed, k), cy(isd:ied, js:je+1, k), npx, &
&                          npy, hord, fx, fy, xfx(is:ie+1, jsd:jed, k), &
&                          yfx(isd:ied, js:je+1, k), gridstruct, bd, &
&                          ra_x, ra_y, mfx=mfx(is:ie+1, js:je, k), mfy=&
&                          mfy(is:ie, js:je+1, k))
            CALL PUSHCONTROL(2,1)
          ELSE
            CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
            CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
            CALL PUSHREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1)*&
&                         (jed-jsd+1))
            CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd:&
&                   jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                   hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(isd&
&                   :ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, mfx=&
&                   mfx(is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1, k))
            CALL PUSHCONTROL(2,0)
          END IF
          DO j=js,je
            DO i=is,ie
              CALL PUSHREALARRAY(q(i, j, k, iq))
              q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-fx&
&               (i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
        END DO
      END DO
! npz
      IF (it .NE. nsplt) THEN
        CALL PUSHREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz*&
&                     nq)
        CALL START_GROUP_HALO_UPDATE(q_pack, q, domain)
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
!Apply nested-grid BCs
      IF (gridstruct%nested) THEN
        DO iq=1,nq
          arg1 = REAL(neststruct%tracer_nest_timestep)
          arg2 = REAL(nsplt*k_split)
          CALL PUSHREALARRAY(q(isd:ied, jsd:jed, :, iq), (ied-isd+1)*(&
&                       jed-jsd+1)*npz)
          CALL NESTED_GRID_BC_APPLY_INTT(q(isd:ied, jsd:jed, :, iq), 0, &
&                                  0, npx, npy, npz, bd, REAL(neststruct&
&                                  %tracer_nest_timestep), REAL(nsplt*&
&                                  k_split), neststruct%q_bc(iq), bctype&
&                                  =neststruct%nestbctype)
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
! nsplt
    IF (id_divg .GT. 0) THEN
      rdt = 1./(frac*dt)
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dp1,xfx,yfx,rarea,rdt)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(dp1(i, j, k))
            dp1(i, j, k) = (xfx(i+1, j, k)-xfx(i, j, k)+(yfx(i, j+1, k)-&
&             yfx(i, j, k)))*rarea(i, j)*rdt
          END DO
        END DO
      END DO
      CALL PUSHINTEGER(je)
      CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL PUSHINTEGER(nsplt)
      CALL PUSHREALARRAY(rdt)
      CALL PUSHINTEGER(is)
      CALL PUSHREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
      CALL PUSHREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
      CALL PUSHINTEGER(ie)
      CALL PUSHREALARRAY(frac)
      CALL PUSHREALARRAY(dp2, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
      !CALL PUSHPOINTER8(C_LOC(rarea))
      CALL PUSHREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
      CALL PUSHINTEGER(js)
      CALL PUSHCONTROL(1,1)
    ELSE
      CALL PUSHINTEGER(je)
      CALL PUSHREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL PUSHREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL PUSHINTEGER(nsplt)
      CALL PUSHINTEGER(is)
      CALL PUSHREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
      CALL PUSHREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
      CALL PUSHREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
      CALL PUSHINTEGER(ie)
      CALL PUSHREALARRAY(frac)
      CALL PUSHREALARRAY(dp2, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
      !CALL PUSHPOINTER8(C_LOC(rarea))
      CALL PUSHREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
      CALL PUSHINTEGER(js)
      CALL PUSHCONTROL(1,0)
    END IF
  END SUBROUTINE TRACER_2D_NESTED_FWD
!  Differentiation of tracer_2d_nested in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b
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
!   gradient     of useful results: q dp1
!   with respect to varying inputs: q dp1 mfx mfy cx cy
  SUBROUTINE TRACER_2D_NESTED_BWD(q, q_ad, dp1, dp1_ad, mfx, mfx_ad, mfy&
&   , mfy_ad, cx, cx_ad, cy, cy_ad, gridstruct, bd, domain, npx, npy, &
&   npz, nq, hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, k_split&
&   , neststruct, parent_grid, hord_pert, nord_tr_pert, trdm_pert, &
&   split_damp_tr)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord, nord_tr
    INTEGER, INTENT(IN) :: hord_pert, nord_tr_pert
    LOGICAL, INTENT(IN) :: split_damp_tr
    INTEGER, INTENT(IN) :: q_split, k_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt, trdm, trdm_pert
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: q_pack
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
    REAL, INTENT(INOUT) :: q_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dp1_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_ad(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_ad(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: dp2_ad(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_ad(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_ad(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_ad(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_ad(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx_ad(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx_ad(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cmax(npz)
    REAL :: cmax_t
    REAL :: c_global
    REAL :: frac, rdt
    INTEGER :: nsplt, nsplt_parent
    INTEGER, SAVE :: msg_split_steps=1
    INTEGER :: i, j, k, it, iq
    REAL, DIMENSION(:, :), POINTER :: area, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: dxa, dya, dx, dy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    REAL :: max1
    REAL :: arg1
    REAL :: arg2
    REAL :: temp_ad
    REAL :: temp
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    INTEGER :: branch
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_tracer_2d_nested
    REAL :: x2
    REAL :: x1
    REAL :: y2
    REAL :: y1

    dp2 = 0.0
    fx = 0.0
    fy = 0.0
    ra_x = 0.0
    ra_y = 0.0
    xfx = 0.0
    yfx = 0.0
    cmax = 0.0
    cmax_t = 0.0
    c_global = 0.0
    frac = 0.0
    rdt = 0.0
    nsplt = 0
    nsplt_parent = 0
    it = 0
    iq = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    isd = 0
    ied = 0
    jsd = 0
    jed = 0
    max1 = 0.0
    arg1 = 0.0
    arg2 = 0.0
    x2 = 0.0
    x1 = 0.0
    y2 = 0.0
    y1 = 0.0
    branch = 0

    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER(js)
      CALL POPREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
      !CALL POPPOINTER8(cptr)
      rarea => gridstruct%rarea ! (/unknown_shape_in_tracer_2d_nested&
!&                /))
      CALL POPREALARRAY(dp2, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
      CALL POPREALARRAY(frac)
      CALL POPINTEGER(ie)
      CALL POPREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
      CALL POPREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
      CALL POPINTEGER(is)
      CALL POPINTEGER(nsplt)
      CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL POPINTEGER(je)
      xfx_ad = 0.0
      yfx_ad = 0.0
    ELSE
      CALL POPINTEGER(js)
      CALL POPREALARRAY(xfx, (bd%ie-bd%is+2)*(bd%jed-bd%jsd+1)*npz)
      !CALL POPPOINTER8(cptr)
      rarea => gridstruct%rarea ! (/unknown_shape_in_tracer_2d_nested&
!i&                /))
      CALL POPREALARRAY(dp2, (bd%ie-bd%is+1)*(bd%je-bd%js+1))
      CALL POPREALARRAY(frac)
      CALL POPINTEGER(ie)
      CALL POPREALARRAY(ra_x, (bd%ie-bd%is+1)*(bd%jed-bd%jsd+1))
      CALL POPREALARRAY(ra_y, (bd%ied-bd%isd+1)*(bd%je-bd%js+1))
      CALL POPREALARRAY(yfx, (bd%ied-bd%isd+1)*(bd%je-bd%js+2)*npz)
      CALL POPINTEGER(is)
      CALL POPREALARRAY(rdt)
      CALL POPINTEGER(nsplt)
      CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
      CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
      CALL POPINTEGER(je)
      xfx_ad = 0.0
      yfx_ad = 0.0
      DO k=npz,1,-1
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(dp1(i, j, k))
            temp_ad2 = rarea(i, j)*rdt*dp1_ad(i, j, k)
            xfx_ad(i+1, j, k) = xfx_ad(i+1, j, k) + temp_ad2
            xfx_ad(i, j, k) = xfx_ad(i, j, k) - temp_ad2
            yfx_ad(i, j+1, k) = yfx_ad(i, j+1, k) + temp_ad2
            yfx_ad(i, j, k) = yfx_ad(i, j, k) - temp_ad2
            dp1_ad(i, j, k) = 0.0
          END DO
        END DO
      END DO
    END IF
    jsd = bd%jsd
    ied = bd%ied
    isd = bd%isd
    jed = bd%jed
    mfx_ad = 0.0
    mfy_ad = 0.0
    cx_ad = 0.0
    cy_ad = 0.0
    dp2_ad = 0.0
    ra_x_ad = 0.0
    ra_y_ad = 0.0
    fx_ad = 0.0
    fy_ad = 0.0
    DO it=nsplt,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
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
        CALL POPREALARRAY(q, (bd%ied-bd%isd+1)*(bd%jed-bd%jsd+1)*npz*nq&
&                   )
        CALL START_GROUP_HALO_UPDATE_ADM(q_pack, q, q_ad, domain)
      END IF
      DO k=npz,1,-1
        DO iq=nq,1,-1
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(q(i, j, k, iq))
              temp_ad0 = q_ad(i, j, k, iq)/dp2(i, j)
              temp = q(i, j, k, iq)
              temp_ad1 = rarea(i, j)*temp_ad0
              dp1_ad(i, j, k) = dp1_ad(i, j, k) + temp*temp_ad0
              fx_ad(i, j) = fx_ad(i, j) + temp_ad1
              fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad1
              fy_ad(i, j) = fy_ad(i, j) + temp_ad1
              fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad1
              dp2_ad(i, j) = dp2_ad(i, j) - (temp*dp1(i, j, k)+rarea(i, &
&               j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*temp_ad0/&
&               dp2(i, j)
              q_ad(i, j, k, iq) = dp1(i, j, k)*temp_ad0
            END DO
          END DO
          CALL POPCONTROL(2,branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              CALL POPREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1)&
&                          *(jed-jsd+1))
              CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
              CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
              CALL FV_TP_2D_ADM(q(isd:ied, jsd:jed, k, iq), q_ad(isd:ied&
&                         , jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k), &
&                         cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied, js:je+&
&                         1, k), cy_ad(isd:ied, js:je+1, k), npx, npy, &
&                         hord_pert, fx, fx_ad, fy, fy_ad, xfx(is:ie+1, &
&                         jsd:jed, k), xfx_ad(is:ie+1, jsd:jed, k), yfx(&
&                         isd:ied, js:je+1, k), yfx_ad(isd:ied, js:je+1&
&                         , k), gridstruct, bd, ra_x, ra_x_ad, ra_y, &
&                         ra_y_ad, mfx(is:ie+1, js:je, k), mfx_ad(is:ie+&
&                         1, js:je, k), mfy(is:ie, js:je+1, k), mfy_ad(&
&                         is:ie, js:je+1, k))
            ELSE
              CALL FV_TP_2D_BWD(q(isd:ied, jsd:jed, k, iq), q_ad(isd:&
&                            ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, &
&                            k), cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied&
&                            , js:je+1, k), cy_ad(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fx_ad, fy, fy_ad, xfx(&
&                            is:ie+1, jsd:jed, k), xfx_ad(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), yfx_ad(&
&                            isd:ied, js:je+1, k), gridstruct, bd, ra_x&
&                            , ra_x_ad, ra_y, ra_y_ad, mfx(is:ie+1, js:&
&                            je, k), mfx_ad(is:ie+1, js:je, k), mfy(is:&
&                            ie, js:je+1, k), mfy_ad(is:ie, js:je+1, k))
            END IF
          ELSE IF (branch .EQ. 2) THEN
            CALL POPREALARRAY(q(isd:ied, jsd:jed, k, iq), (ied-isd+1)*(&
&                        jed-jsd+1))
            CALL POPREALARRAY(fx, (bd%ie-bd%is+2)*(bd%je-bd%js+1))
            CALL POPREALARRAY(fy, (bd%ie-bd%is+1)*(bd%je-bd%js+2))
            CALL FV_TP_2D_ADM(q(isd:ied, jsd:jed, k, iq), q_ad(isd:ied, &
&                       jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k), cx_ad(&
&                       is:ie+1, jsd:jed, k), cy(isd:ied, js:je+1, k), &
&                       cy_ad(isd:ied, js:je+1, k), npx, npy, hord_pert&
&                       , fx, fx_ad, fy, fy_ad, xfx(is:ie+1, jsd:jed, k)&
&                       , xfx_ad(is:ie+1, jsd:jed, k), yfx(isd:ied, js:&
&                       je+1, k), yfx_ad(isd:ied, js:je+1, k), &
&                       gridstruct, bd, ra_x, ra_x_ad, ra_y, ra_y_ad, &
&                       mfx(is:ie+1, js:je, k), mfx_ad(is:ie+1, js:je, k&
&                       ), mfy(is:ie, js:je+1, k), mfy_ad(is:ie, js:je+1&
&                       , k), dp1(isd:ied, jsd:jed, k), dp1_ad(isd:ied, &
&                       jsd:jed, k), nord=nord_tr_pert, damp_c=trdm_pert&
&                      )
          ELSE
            CALL FV_TP_2D_BWD(q(isd:ied, jsd:jed, k, iq), q_ad(isd:&
&                          ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k)&
&                          , cx_ad(is:ie+1, jsd:jed, k), cy(isd:ied, js:&
&                          je+1, k), cy_ad(isd:ied, js:je+1, k), npx, &
&                          npy, hord, fx, fx_ad, fy, fy_ad, xfx(is:ie+1&
&                          , jsd:jed, k), xfx_ad(is:ie+1, jsd:jed, k), &
&                          yfx(isd:ied, js:je+1, k), yfx_ad(isd:ied, js:&
&                          je+1, k), gridstruct, bd, ra_x, ra_x_ad, ra_y&
&                          , ra_y_ad, mfx(is:ie+1, js:je, k), mfx_ad(is:&
&                          ie+1, js:je, k), mfy(is:ie, js:je+1, k), &
&                          mfy_ad(is:ie, js:je+1, k), dp1(isd:ied, jsd:&
&                          jed, k), dp1_ad(isd:ied, jsd:jed, k), nord=&
&                          nord_tr, damp_c=trdm)
          END IF
        END DO
        DO j=je,js,-1
          DO i=ied,isd,-1
            CALL POPREALARRAY(ra_y(i, j))
            yfx_ad(i, j, k) = yfx_ad(i, j, k) + ra_y_ad(i, j)
            yfx_ad(i, j+1, k) = yfx_ad(i, j+1, k) - ra_y_ad(i, j)
            ra_y_ad(i, j) = 0.0
          END DO
        END DO
        DO j=jed,jsd,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(ra_x(i, j))
            xfx_ad(i, j, k) = xfx_ad(i, j, k) + ra_x_ad(i, j)
            xfx_ad(i+1, j, k) = xfx_ad(i+1, j, k) - ra_x_ad(i, j)
            ra_x_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(dp2(i, j))
            temp_ad = rarea(i, j)*dp2_ad(i, j)
            dp1_ad(i, j, k) = dp1_ad(i, j, k) + dp2_ad(i, j)
            mfx_ad(i, j, k) = mfx_ad(i, j, k) + temp_ad
            mfx_ad(i+1, j, k) = mfx_ad(i+1, j, k) - temp_ad
            mfy_ad(i, j, k) = mfy_ad(i, j, k) + temp_ad
            mfy_ad(i, j+1, k) = mfy_ad(i, j+1, k) - temp_ad
            dp2_ad(i, j) = 0.0
          END DO
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      IF (branch .NE. 0) THEN
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
    END DO
    CALL POPCONTROL(1,branch)
    IF (branch .NE. 0) THEN
      DO k=npz,1,-1
        DO j=je+1,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(mfy(i, j, k))
            mfy_ad(i, j, k) = frac*mfy_ad(i, j, k)
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ied,isd,-1
            yfx_ad(i, j, k) = frac*yfx_ad(i, j, k)
            CALL POPREALARRAY(cy(i, j, k))
            cy_ad(i, j, k) = frac*cy_ad(i, j, k)
          END DO
        END DO
        DO j=je,js,-1
          DO i=ie+1,is,-1
            CALL POPREALARRAY(mfx(i, j, k))
            mfx_ad(i, j, k) = frac*mfx_ad(i, j, k)
          END DO
        END DO
        DO j=jed,jsd,-1
          DO i=ie+1,is,-1
            xfx_ad(i, j, k) = frac*xfx_ad(i, j, k)
            CALL POPREALARRAY(cx(i, j, k))
            cx_ad(i, j, k) = frac*cx_ad(i, j, k)
          END DO
        END DO
      END DO
    END IF
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO k=npz,2,-1
          CALL POPCONTROL(1,branch)
        END DO
      END IF
      sin_sg => gridstruct%sin_sg
      DO k=npz,1,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPCONTROL(1,branch)
            END DO
          END DO
        ELSE
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPCONTROL(1,branch)
            END DO
          END DO
        END IF
      END DO
    ELSE
      sin_sg => gridstruct%sin_sg
    END IF
    dxa => gridstruct%dxa
    dx => gridstruct%dx
    dy => gridstruct%dy
    dya => gridstruct%dya
    DO k=npz,1,-1
      DO j=je+1,js,-1
        DO i=ied,isd,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            cy_ad(i, j, k) = cy_ad(i, j, k) + dya(i, j)*dx(i, j)*sin_sg(&
&             i, j, 2)*yfx_ad(i, j, k)
            yfx_ad(i, j, k) = 0.0
          ELSE
            cy_ad(i, j, k) = cy_ad(i, j, k) + dya(i, j-1)*dx(i, j)*&
&             sin_sg(i, j-1, 4)*yfx_ad(i, j, k)
            yfx_ad(i, j, k) = 0.0
          END IF
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL(1,branch)
          IF (branch .EQ. 0) THEN
            cx_ad(i, j, k) = cx_ad(i, j, k) + dxa(i, j)*dy(i, j)*sin_sg(&
&             i, j, 1)*xfx_ad(i, j, k)
            xfx_ad(i, j, k) = 0.0
          ELSE
            cx_ad(i, j, k) = cx_ad(i, j, k) + dxa(i-1, j)*dy(i, j)*&
&             sin_sg(i-1, j, 3)*xfx_ad(i, j, k)
            xfx_ad(i, j, k) = 0.0
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE TRACER_2D_NESTED_BWD
  SUBROUTINE TRACER_2D_NESTED(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, &
&   domain, npx, npy, npz, nq, hord, q_split, dt, id_divg, q_pack, &
&   nord_tr, trdm, k_split, neststruct, parent_grid, hord_pert, &
&   nord_tr_pert, trdm_pert, split_damp_tr)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord, nord_tr
    INTEGER, INTENT(IN) :: hord_pert, nord_tr_pert
    LOGICAL, INTENT(IN) :: split_damp_tr
    INTEGER, INTENT(IN) :: q_split, k_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt, trdm, trdm_pert
    TYPE(GROUP_HALO_UPDATE_TYPE), INTENT(INOUT) :: q_pack
! Tracers
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Local Arrays
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: cmax(npz)
    REAL :: cmax_t
    REAL :: c_global
    REAL :: frac, rdt
    INTEGER :: nsplt, nsplt_parent
    INTEGER, SAVE :: msg_split_steps=1
    INTEGER :: i, j, k, it, iq
    REAL, DIMENSION(:, :), POINTER :: area, rarea
    REAL, DIMENSION(:, :, :), POINTER :: sin_sg
    REAL, DIMENSION(:, :), POINTER :: dxa, dya, dx, dy
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    REAL :: max1
    REAL :: arg1
    REAL :: arg2
    REAL :: x2
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
    area => gridstruct%area
    rarea => gridstruct%rarea
    sin_sg => gridstruct%sin_sg
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    dx => gridstruct%dx
    dy => gridstruct%dy
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx)
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sin_sg(i-1, &
&             j, 3)
          ELSE
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1&
&             )
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sin_sg(i, j-&
&             1, 4)
          ELSE
            yfx(i, j, k) = cy(i, j, k)*dya(i, j)*dx(i, j)*sin_sg(i, j, 2&
&             )
          END IF
        END DO
      END DO
    END DO
!--------------------------------------------------------------------------------
    IF (q_split .EQ. 0) THEN
! Determine nsplt
!$OMP parallel do default(none) shared(is,ie,js,je,npz,cmax,cx,cy,sin_sg) &
!$OMP                          private(cmax_t )
      DO k=1,npz
        cmax(k) = 0.
        IF (k .LT. 4) THEN
! Top layers: C < max( abs(c_x), abs(c_y) )
          DO j=js,je
            DO i=is,ie
              IF (cx(i, j, k) .GE. 0.) THEN
                x1 = cx(i, j, k)
              ELSE
                x1 = -cx(i, j, k)
              END IF
              IF (cy(i, j, k) .GE. 0.) THEN
                y1 = cy(i, j, k)
              ELSE
                y1 = -cy(i, j, k)
              END IF
              IF (x1 .LT. y1) THEN
                cmax_t = y1
              ELSE
                cmax_t = x1
              END IF
              IF (cmax_t .LT. cmax(k)) THEN
                cmax(k) = cmax(k)
              ELSE
                cmax(k) = cmax_t
              END IF
            END DO
          END DO
        ELSE
          DO j=js,je
            DO i=is,ie
              IF (cx(i, j, k) .GE. 0.) THEN
                x2 = cx(i, j, k)
              ELSE
                x2 = -cx(i, j, k)
              END IF
              IF (cy(i, j, k) .GE. 0.) THEN
                y2 = cy(i, j, k)
              ELSE
                y2 = -cy(i, j, k)
              END IF
              IF (x2 .LT. y2) THEN
                max1 = y2
              ELSE
                max1 = x2
              END IF
              cmax_t = max1 + 1. - sin_sg(i, j, 5)
              IF (cmax_t .LT. cmax(k)) THEN
                cmax(k) = cmax(k)
              ELSE
                cmax(k) = cmax_t
              END IF
            END DO
          END DO
        END IF
      END DO
      CALL MP_REDUCE_MAX(cmax, npz)
! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      IF (npz .NE. 1) THEN
! if NOT shallow water test case
        DO k=2,npz
          IF (cmax(k) .LT. c_global) THEN
            c_global = c_global
          ELSE
            c_global = cmax(k)
          END IF
        END DO
      END IF
      nsplt = INT(1. + c_global)
      IF (IS_MASTER() .AND. nsplt .GT. 3) WRITE(*, *) 'Tracer_2d_split='&
&                                         , nsplt, c_global
    ELSE
      nsplt = q_split
      IF (gridstruct%nested .AND. neststruct%nestbctype .GT. 1) THEN
        IF (q_split/parent_grid%flagstruct%q_split .LT. 1) THEN
          msg_split_steps = 1
        ELSE
          msg_split_steps = q_split/parent_grid%flagstruct%q_split
        END IF
      END IF
    END IF
!--------------------------------------------------------------------------------
    frac = 1./REAL(nsplt)
    IF (nsplt .NE. 1) THEN
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,frac,xfx,mfx,cy,yfx,mfy)
      DO k=1,npz
        DO j=jsd,jed
          DO i=is,ie+1
            cx(i, j, k) = cx(i, j, k)*frac
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            cy(i, j, k) = cy(i, j, k)*frac
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END DO
    END IF
    DO it=1,nsplt
      IF (gridstruct%nested) neststruct%tracer_nest_timestep = &
&         neststruct%tracer_nest_timestep + 1
      CALL TIMING_ON('COMM_TOTAL')
      CALL TIMING_ON('COMM_TRACER')
      CALL COMPLETE_GROUP_HALO_UPDATE(q_pack, domain)
      CALL TIMING_OFF('COMM_TRACER')
      CALL TIMING_OFF('COMM_TOTAL')
      IF (gridstruct%nested) THEN
        DO iq=1,nq
          arg1 = REAL(neststruct%tracer_nest_timestep) + REAL(nsplt*&
&           k_split)
          arg2 = REAL(nsplt*k_split)
          CALL NESTED_GRID_BC_APPLY_INTT(q(isd:ied, jsd:jed, :, iq), 0, &
&                                  0, npx, npy, npz, bd, arg1, arg2, &
&                                  neststruct%q_bc(iq), neststruct%&
&                                  nestbctype)
        END DO
      END IF
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx,mfy,rarea,nq, &
!$OMP                                  area,xfx,yfx,q,cx,cy,npx,npy,hord,gridstruct,bd,it,nsplt,nord_tr,trdm) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+(mfy&
&             (i, j, k)-mfy(i, j+1, k)))*rarea(i, j)
          END DO
        END DO
        DO j=jsd,jed
          DO i=is,ie
            ra_x(i, j) = area(i, j) + (xfx(i, j, k)-xfx(i+1, j, k))
          END DO
        END DO
        DO j=js,je
          DO i=isd,ied
            ra_y(i, j) = area(i, j) + (yfx(i, j, k)-yfx(i, j+1, k))
          END DO
        END DO
        DO iq=1,nq
          IF (it .EQ. 1 .AND. trdm .GT. 1.e-4) THEN
            IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, &
&                        jsd:jed, k), cy(isd:ied, js:je+1, k), npx, npy&
&                        , hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                        isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                        ra_y, mfx(is:ie+1, js:je, k), mfy(is:ie, js:je+&
&                        1, k), dp1(isd:ied, jsd:jed, k), nord_tr, trdm)
            ELSE
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd:&
&                     jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                     hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                     isd:ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, &
&                     mfx(is:ie+1, js:je, k), mfy(is:ie, js:je+1, k), &
&                     dp1(isd:ied, jsd:jed, k), nord_tr, trdm)
            END IF
          ELSE IF (hord .EQ. hord_pert) THEN
            CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd&
&                      :jed, k), cy(isd:ied, js:je+1, k), npx, npy, hord&
&                      , fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(isd:ied, &
&                      js:je+1, k), gridstruct, bd, ra_x, ra_y, mfx=mfx(&
&                      is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1, k))
          ELSE
            CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd:&
&                   jed, k), cy(isd:ied, js:je+1, k), npx, npy, &
&                   hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(isd&
&                   :ied, js:je+1, k), gridstruct, bd, ra_x, ra_y, mfx=&
&                   mfx(is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1, k))
          END IF
          DO j=js,je
            DO i=is,ie
              q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-fx&
&               (i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
        END DO
      END DO
! npz
      IF (it .NE. nsplt) THEN
        CALL TIMING_ON('COMM_TOTAL')
        CALL TIMING_ON('COMM_TRACER')
        CALL START_GROUP_HALO_UPDATE(q_pack, q, domain)
        CALL TIMING_OFF('COMM_TRACER')
        CALL TIMING_OFF('COMM_TOTAL')
      END IF
!Apply nested-grid BCs
      IF (gridstruct%nested) THEN
        DO iq=1,nq
          arg1 = REAL(neststruct%tracer_nest_timestep)
          arg2 = REAL(nsplt*k_split)
          CALL NESTED_GRID_BC_APPLY_INTT(q(isd:ied, jsd:jed, :, iq), 0, &
&                                  0, npx, npy, npz, bd, arg1, arg2, &
&                                  neststruct%q_bc(iq), neststruct%&
&                                  nestbctype)
        END DO
      END IF
    END DO
! nsplt
    IF (id_divg .GT. 0) THEN
      rdt = 1./(frac*dt)
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dp1,xfx,yfx,rarea,rdt)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dp1(i, j, k) = (xfx(i+1, j, k)-xfx(i, j, k)+(yfx(i, j+1, k)-&
&             yfx(i, j, k)))*rarea(i, j)*rdt
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE TRACER_2D_NESTED
end module fv_tracer2d_adm_mod
