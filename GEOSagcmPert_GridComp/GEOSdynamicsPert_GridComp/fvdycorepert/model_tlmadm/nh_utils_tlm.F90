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
module nh_utils_tlm_mod
! Developer: S.-J. Lin, NOAA/GFDL
! To do list:
! include moisture effect in pt
!------------------------------
   use constants_mod,     only: rdgas, cp_air, grav
   use tp_core_tlm_mod,   only: fv_tp_2d
   use tp_core_tlm_mod,   only: fv_tp_2d_tlm
   use sw_core_tlm_mod,   only: fill_4corners, del6_vt_flux
   use sw_core_tlm_mod,   only: fill_4corners_tlm, del6_vt_flux_tlm
   use fv_arrays_mod,     only: fv_grid_bounds_type, fv_grid_type

   implicit none
   private

   public update_dz_c, update_dz_d, nest_halo_nh
   public sim_solver, sim1_solver, sim3_solver
   public sim3p0_solver, rim_2d
   public Riem_Solver_c
   public update_dz_c_tlm, update_dz_d_tlm, nest_halo_nh_tlm
   public sim_solver_tlm, sim1_solver_tlm, sim3_solver_tlm
   public sim3p0_solver_tlm, rim_2d_tlm
   public Riem_Solver_c_tlm

   real, parameter:: dz_min = 2.
   real, parameter:: r3 = 1./3.

CONTAINS
!  Differentiation of update_dz_c in forward (tangent) mode:
!   variations   of useful results: ws gz
!   with respect to varying inputs: ws gz ut vt
  SUBROUTINE UPDATE_DZ_C_TLM(is, ie, js, je, km, ng, dt, dp0, zs, area, &
&   ut, ut_tl, vt, vt_tl, gz, gz_tl, ws, ws_tl, npx, npy, sw_corner, &
&   se_corner, ne_corner, nw_corner, bd, grid_type)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km, npx, npy, grid_type
    LOGICAL, INTENT(IN) :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: dp0(km)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: ut, vt
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: ut_tl, &
&   vt_tl
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng), INTENT(IN) :: area
    REAL, INTENT(INOUT) :: gz(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(INOUT) :: gz_tl(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(IN) :: zs(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(OUT) :: ws(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(OUT) :: ws_tl(is-ng:ie+ng, js-ng:je+ng)
! Local Work array:
    REAL :: gz2(is-ng:ie+ng, js-ng:je+ng)
    REAL :: gz2_tl(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-1:ie+2, js-1:je+1) :: xfx, fx
    REAL, DIMENSION(is-1:ie+2, js-1:je+1) :: xfx_tl, fx_tl
    REAL, DIMENSION(is-1:ie+1, js-1:je+2) :: yfx, fy
    REAL, DIMENSION(is-1:ie+1, js-1:je+2) :: yfx_tl, fy_tl
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
    xfx_tl = 0.0
    gz2_tl = 0.0
    yfx_tl = 0.0
    fx_tl = 0.0
    fy_tl = 0.0
!$OMP parallel do default(none) shared(js1,je1,is1,ie2,km,je2,ie1,ut,top_ratio,vt, &
!$OMP                                  bot_ratio,dp0,js,je,ng,is,ie,gz,grid_type,  &
!$OMP                                  bd,npx,npy,sw_corner,se_corner,ne_corner,   &
!$OMP                                  nw_corner,area) &
!$OMP                          private(gz2, xfx, yfx, fx, fy, int_ratio)
    DO k=1,km+1
      IF (k .EQ. 1) THEN
        DO j=js1,je1
          DO i=is1,ie2
            xfx_tl(i, j) = ut_tl(i, j, 1) + top_ratio*(ut_tl(i, j, 1)-&
&             ut_tl(i, j, 2))
            xfx(i, j) = ut(i, j, 1) + (ut(i, j, 1)-ut(i, j, 2))*&
&             top_ratio
          END DO
        END DO
        DO j=js1,je2
          DO i=is1,ie1
            yfx_tl(i, j) = vt_tl(i, j, 1) + top_ratio*(vt_tl(i, j, 1)-&
&             vt_tl(i, j, 2))
            yfx(i, j) = vt(i, j, 1) + (vt(i, j, 1)-vt(i, j, 2))*&
&             top_ratio
          END DO
        END DO
      ELSE IF (k .EQ. km + 1) THEN
! Bottom extrapolation
        DO j=js1,je1
          DO i=is1,ie2
            xfx_tl(i, j) = ut_tl(i, j, km) + bot_ratio*(ut_tl(i, j, km)-&
&             ut_tl(i, j, km-1))
            xfx(i, j) = ut(i, j, km) + (ut(i, j, km)-ut(i, j, km-1))*&
&             bot_ratio
          END DO
        END DO
!             xfx(i,j) = r14*(3.*ut(i,j,km-2)-13.*ut(i,j,km-1)+24.*ut(i,j,km))
!             if ( xfx(i,j)*ut(i,j,km)<0. ) xfx(i,j) = 0.
        DO j=js1,je2
          DO i=is1,ie1
            yfx_tl(i, j) = vt_tl(i, j, km) + bot_ratio*(vt_tl(i, j, km)-&
&             vt_tl(i, j, km-1))
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
            xfx_tl(i, j) = int_ratio*(dp0(k)*ut_tl(i, j, k-1)+dp0(k-1)*&
&             ut_tl(i, j, k))
            xfx(i, j) = (dp0(k)*ut(i, j, k-1)+dp0(k-1)*ut(i, j, k))*&
&             int_ratio
          END DO
        END DO
        DO j=js1,je2
          DO i=is1,ie1
            yfx_tl(i, j) = int_ratio*(dp0(k)*vt_tl(i, j, k-1)+dp0(k-1)*&
&             vt_tl(i, j, k))
            yfx(i, j) = (dp0(k)*vt(i, j, k-1)+dp0(k-1)*vt(i, j, k))*&
&             int_ratio
          END DO
        END DO
      END IF
      DO j=js-ng,je+ng
        DO i=is-ng,ie+ng
          gz2_tl(i, j) = gz_tl(i, j, k)
          gz2(i, j) = gz(i, j, k)
        END DO
      END DO
      IF (grid_type .LT. 3) CALL FILL_4CORNERS_TLM(gz2, gz2_tl, 1, bd, &
&                                            npx, npy, sw_corner, &
&                                            se_corner, ne_corner, &
&                                            nw_corner)
      DO j=js1,je1
        DO i=is1,ie2
          IF (xfx(i, j) .GT. 0.) THEN
            fx_tl(i, j) = gz2_tl(i-1, j)
            fx(i, j) = gz2(i-1, j)
          ELSE
            fx_tl(i, j) = gz2_tl(i, j)
            fx(i, j) = gz2(i, j)
          END IF
          fx_tl(i, j) = xfx_tl(i, j)*fx(i, j) + xfx(i, j)*fx_tl(i, j)
          fx(i, j) = xfx(i, j)*fx(i, j)
        END DO
      END DO
      IF (grid_type .LT. 3) CALL FILL_4CORNERS_TLM(gz2, gz2_tl, 2, bd, &
&                                            npx, npy, sw_corner, &
&                                            se_corner, ne_corner, &
&                                            nw_corner)
      DO j=js1,je2
        DO i=is1,ie1
          IF (yfx(i, j) .GT. 0.) THEN
            fy_tl(i, j) = gz2_tl(i, j-1)
            fy(i, j) = gz2(i, j-1)
          ELSE
            fy_tl(i, j) = gz2_tl(i, j)
            fy(i, j) = gz2(i, j)
          END IF
          fy_tl(i, j) = yfx_tl(i, j)*fy(i, j) + yfx(i, j)*fy_tl(i, j)
          fy(i, j) = yfx(i, j)*fy(i, j)
        END DO
      END DO
      DO j=js1,je1
        DO i=is1,ie1
          gz_tl(i, j, k) = ((area(i, j)*gz2_tl(i, j)+fx_tl(i, j)-fx_tl(i&
&           +1, j)+fy_tl(i, j)-fy_tl(i, j+1))*(area(i, j)+(xfx(i, j)-xfx&
&           (i+1, j))+(yfx(i, j)-yfx(i, j+1)))-(gz2(i, j)*area(i, j)+(fx&
&           (i, j)-fx(i+1, j))+(fy(i, j)-fy(i, j+1)))*(xfx_tl(i, j)-&
&           xfx_tl(i+1, j)+yfx_tl(i, j)-yfx_tl(i, j+1)))/(area(i, j)+(&
&           xfx(i, j)-xfx(i+1, j))+(yfx(i, j)-yfx(i, j+1)))**2
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
        ws_tl(i, j) = -(rdt*gz_tl(i, j, km+1))
        ws(i, j) = (zs(i, j)-gz(i, j, km+1))*rdt
      END DO
      DO k=km,1,-1
        DO i=is1,ie1
          IF (gz(i, j, k) .LT. gz(i, j, k+1) + dz_min) THEN
            gz_tl(i, j, k) = gz_tl(i, j, k+1)
            gz(i, j, k) = gz(i, j, k+1) + dz_min
          ELSE
            gz(i, j, k) = gz(i, j, k)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE UPDATE_DZ_C_TLM
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
!  Differentiation of update_dz_d in forward (tangent) mode:
!   variations   of useful results: ws zh
!   with respect to varying inputs: xfx ws yfx zh crx cry
  SUBROUTINE UPDATE_DZ_D_TLM(ndif, damp, hord, is, ie, js, je, km, ng, &
&   npx, npy, area, rarea, dp0, zs, zh, zh_tl, crx, crx_tl, cry, cry_tl&
&   , xfx, xfx_tl, yfx, yfx_tl, delz, ws, ws_tl, rdt, gridstruct, bd, &
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
    REAL, INTENT(INOUT) :: zh_tl(is-ng:ie+ng, js-ng:je+ng, km+1)
    REAL, INTENT(OUT) :: delz(is-ng:ie+ng, js-ng:je+ng, km)
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km), INTENT(INOUT) :: crx, xfx
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km), INTENT(INOUT) :: crx_tl, &
&   xfx_tl
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km), INTENT(INOUT) :: cry, yfx
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km), INTENT(INOUT) :: cry_tl, &
&   yfx_tl
    REAL, INTENT(OUT) :: ws(is:ie, js:je)
    REAL, INTENT(OUT) :: ws_tl(is:ie, js:je)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
!-----------------------------------------------------
! Local array:
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km+1) :: crx_adv, xfx_adv
    REAL, DIMENSION(is:ie+1, js-ng:je+ng, km+1) :: crx_adv_tl, &
&   xfx_adv_tl
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km+1) :: cry_adv, yfx_adv
    REAL, DIMENSION(is-ng:ie+ng, js:je+1, km+1) :: cry_adv_tl, &
&   yfx_adv_tl
    REAL, DIMENSION(is:ie+1, js:je) :: fx
    REAL, DIMENSION(is:ie+1, js:je) :: fx_tl
    REAL, DIMENSION(is:ie, js:je+1) :: fy
    REAL, DIMENSION(is:ie, js:je+1) :: fy_tl
    REAL, DIMENSION(is-ng:ie+ng+1, js-ng:je+ng) :: fx2
    REAL, DIMENSION(is-ng:ie+ng+1, js-ng:je+ng) :: fx2_tl
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng+1) :: fy2
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng+1) :: fy2_tl
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng) :: wk2, z2
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng) :: wk2_tl, z2_tl
    REAL :: ra_x(is:ie, js-ng:je+ng)
    REAL :: ra_x_tl(is:ie, js-ng:je+ng)
    REAL :: ra_y(is-ng:ie+ng, js:je)
    REAL :: ra_y_tl(is-ng:ie+ng, js:je)
!--------------------------------------------------------------------
    INTEGER :: i, j, k, isd, ied, jsd, jed
    LOGICAL :: uniform_grid
    INTRINSIC MAX
    uniform_grid = .false.
    damp(km+1) = damp(km)
    ndif(km+1) = ndif(km)
    isd = is - ng
    ied = ie + ng
    jsd = js - ng
    jed = je + ng
    yfx_adv_tl = 0.0
    crx_adv_tl = 0.0
    xfx_adv_tl = 0.0
    cry_adv_tl = 0.0
!$OMP parallel do default(none) shared(jsd,jed,crx,xfx,crx_adv,xfx_adv,is,ie,isd,ied, &
!$OMP                                  km,dp0,uniform_grid,js,je,cry,yfx,cry_adv,yfx_adv)
    DO j=jsd,jed
      CALL EDGE_PROFILE_TLM(crx, crx_tl, xfx, xfx_tl, crx_adv, &
&                     crx_adv_tl, xfx_adv, xfx_adv_tl, is, ie + 1, jsd, &
&                     jed, j, km, dp0, uniform_grid, 0)
      IF (j .LE. je + 1 .AND. j .GE. js) CALL EDGE_PROFILE_TLM(cry, &
&                                                        cry_tl, yfx, &
&                                                        yfx_tl, cry_adv&
&                                                        , cry_adv_tl, &
&                                                        yfx_adv, &
&                                                        yfx_adv_tl, isd&
&                                                        , ied, js, je +&
&                                                        1, j, km, dp0, &
&                                                        uniform_grid, 0&
&                                                       )
    END DO
    fy2_tl = 0.0
    wk2_tl = 0.0
    z2_tl = 0.0
    ra_x_tl = 0.0
    ra_y_tl = 0.0
    fx_tl = 0.0
    fy_tl = 0.0
    fx2_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,km,area,xfx_adv,yfx_adv, &
!$OMP                                  damp,zh,crx_adv,cry_adv,npx,npy,hord,gridstruct,bd,  &
!$OMP                                  ndif,rarea) &
!$OMP                          private(z2, fx2, fy2, ra_x, ra_y, fx, fy,wk2)
    DO k=1,km+1
      DO j=jsd,jed
        DO i=is,ie
          ra_x_tl(i, j) = xfx_adv_tl(i, j, k) - xfx_adv_tl(i+1, j, k)
          ra_x(i, j) = area(i, j) + (xfx_adv(i, j, k)-xfx_adv(i+1, j, k)&
&           )
        END DO
      END DO
      DO j=js,je
        DO i=isd,ied
          ra_y_tl(i, j) = yfx_adv_tl(i, j, k) - yfx_adv_tl(i, j+1, k)
          ra_y(i, j) = area(i, j) + (yfx_adv(i, j, k)-yfx_adv(i, j+1, k)&
&           )
        END DO
      END DO
      IF (damp(k) .GT. 1.e-5) THEN
        DO j=jsd,jed
          DO i=isd,ied
            z2_tl(i, j) = zh_tl(i, j, k)
            z2(i, j) = zh(i, j, k)
          END DO
        END DO
        IF (hord .EQ. hord_pert) THEN
          CALL FV_TP_2D_TLM(z2, z2_tl, crx_adv(is:ie+1, jsd:jed, k), &
&                        crx_adv_tl(is:ie+1, jsd:jed, k), cry_adv(isd:&
&                        ied, js:je+1, k), cry_adv_tl(isd:ied, js:je+1, &
&                        k), npx, npy, hord, fx, fx_tl, fy, fy_tl, &
&                        xfx_adv(is:ie+1, jsd:jed, k), xfx_adv_tl(is:ie+&
&                        1, jsd:jed, k), yfx_adv(isd:ied, js:je+1, k), &
&                        yfx_adv_tl(isd:ied, js:je+1, k), gridstruct, bd&
&                        , ra_x, ra_x_tl, ra_y, ra_y_tl)
        ELSE
          CALL FV_TP_2D_TLM(z2, z2_tl, crx_adv(is:ie+1, jsd:jed, k), &
&                     crx_adv_tl(is:ie+1, jsd:jed, k), cry_adv(isd:ied, &
&                     js:je+1, k), cry_adv_tl(isd:ied, js:je+1, k), npx&
&                     , npy, hord_pert, fx, fx_tl, fy, fy_tl, xfx_adv(is&
&                     :ie+1, jsd:jed, k), xfx_adv_tl(is:ie+1, jsd:jed, k&
&                     ), yfx_adv(isd:ied, js:je+1, k), yfx_adv_tl(isd:&
&                     ied, js:je+1, k), gridstruct, bd, ra_x, ra_x_tl, &
&                     ra_y, ra_y_tl)
      call fv_tp_2d(z2, crx_adv(is:ie+1,jsd:jed,k), cry_adv(isd:ied,js:je+1,k), npx,  npy, hord, &
                   fx, fy, xfx_adv(is:ie+1,jsd:jed,k), yfx_adv(isd:ied,js:je+1,k), gridstruct, bd, ra_x, ra_y)
        END IF
        CALL DEL6_VT_FLUX_TLM(ndif(k), npx, npy, damp(k), z2, z2_tl, wk2&
&                       , wk2_tl, fx2, fx2_tl, fy2, fy2_tl, gridstruct, &
&                       bd)
        DO j=js,je
          DO i=is,ie
            zh_tl(i, j, k) = ((area(i, j)*z2_tl(i, j)+fx_tl(i, j)-fx_tl(&
&             i+1, j)+fy_tl(i, j)-fy_tl(i, j+1))*(ra_x(i, j)+ra_y(i, j)-&
&             area(i, j))-(z2(i, j)*area(i, j)+(fx(i, j)-fx(i+1, j))+(fy&
&             (i, j)-fy(i, j+1)))*(ra_x_tl(i, j)+ra_y_tl(i, j)))/(ra_x(i&
&             , j)+ra_y(i, j)-area(i, j))**2 + rarea(i, j)*(fx2_tl(i, j)&
&             -fx2_tl(i+1, j)+fy2_tl(i, j)-fy2_tl(i, j+1))
            zh(i, j, k) = (z2(i, j)*area(i, j)+(fx(i, j)-fx(i+1, j))+(fy&
&             (i, j)-fy(i, j+1)))/(ra_x(i, j)+ra_y(i, j)-area(i, j)) + (&
&             fx2(i, j)-fx2(i+1, j)+(fy2(i, j)-fy2(i, j+1)))*rarea(i, j)
          END DO
        END DO
      ELSE
        IF (hord .EQ. hord_pert) THEN
          CALL FV_TP_2D_TLM(zh(isd:ied, jsd:jed, k), zh_tl(isd:ied, &
&                        jsd:jed, k), crx_adv(is:ie+1, jsd:jed, k), &
&                        crx_adv_tl(is:ie+1, jsd:jed, k), cry_adv(isd:&
&                        ied, js:je+1, k), cry_adv_tl(isd:ied, js:je+1, &
&                        k), npx, npy, hord, fx, fx_tl, fy, fy_tl, &
&                        xfx_adv(is:ie+1, jsd:jed, k), xfx_adv_tl(is:ie+&
&                        1, jsd:jed, k), yfx_adv(isd:ied, js:je+1, k), &
&                        yfx_adv_tl(isd:ied, js:je+1, k), gridstruct, bd&
&                        , ra_x, ra_x_tl, ra_y, ra_y_tl)
        ELSE
          CALL FV_TP_2D_TLM(zh(isd:ied, jsd:jed, k), zh_tl(isd:ied, jsd:&
&                     jed, k), crx_adv(is:ie+1, jsd:jed, k), crx_adv_tl(&
&                     is:ie+1, jsd:jed, k), cry_adv(isd:ied, js:je+1, k)&
&                     , cry_adv_tl(isd:ied, js:je+1, k), npx, npy, &
&                     hord_pert, fx, fx_tl, fy, fy_tl, xfx_adv(is:ie+1, &
&                     jsd:jed, k), xfx_adv_tl(is:ie+1, jsd:jed, k), &
&                     yfx_adv(isd:ied, js:je+1, k), yfx_adv_tl(isd:ied, &
&                     js:je+1, k), gridstruct, bd, ra_x, ra_x_tl, ra_y, &
&                     ra_y_tl)
      call fv_tp_2d(zh(isd:ied,jsd:jed,k), crx_adv(is:ie+1,jsd:jed,k), cry_adv(isd:ied,js:je+1,k), npx,  npy, hord, &
                    fx, fy, xfx_adv(is:ie+1,jsd:jed,k), yfx_adv(isd:ied,js:je+1,k), gridstruct, bd, ra_x, ra_y)
        END IF
        DO j=js,je
          DO i=is,ie
            zh_tl(i, j, k) = ((area(i, j)*zh_tl(i, j, k)+fx_tl(i, j)-&
&             fx_tl(i+1, j)+fy_tl(i, j)-fy_tl(i, j+1))*(ra_x(i, j)+ra_y(&
&             i, j)-area(i, j))-(zh(i, j, k)*area(i, j)+(fx(i, j)-fx(i+1&
&             , j))+(fy(i, j)-fy(i, j+1)))*(ra_x_tl(i, j)+ra_y_tl(i, j))&
&             )/(ra_x(i, j)+ra_y(i, j)-area(i, j))**2
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
        ws_tl(i, j) = -(rdt*zh_tl(i, j, km+1))
        ws(i, j) = (zs(i, j)-zh(i, j, km+1))*rdt
      END DO
      DO k=km,1,-1
        DO i=is,ie
          IF (zh(i, j, k) .LT. zh(i, j, k+1) + dz_min) THEN
            zh_tl(i, j, k) = zh_tl(i, j, k+1)
            zh(i, j, k) = zh(i, j, k+1) + dz_min
          ELSE
            zh(i, j, k) = zh(i, j, k)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE UPDATE_DZ_D_TLM
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
      CALL EDGE_PROFILE(crx, xfx, crx_adv, xfx_adv, is, ie + 1, jsd, jed&
&                 , j, km, dp0, uniform_grid, 0)
      IF (j .LE. je + 1 .AND. j .GE. js) CALL EDGE_PROFILE(cry, yfx, &
&                                                    cry_adv, yfx_adv, &
&                                                    isd, ied, js, je + &
&                                                    1, j, km, dp0, &
&                                                    uniform_grid, 0)
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
      call fv_tp_2d(z2, crx_adv(is:ie+1,jsd:jed,k), cry_adv(isd:ied,js:je+1,k), npx,  npy, hord, &
                   fx, fy, xfx_adv(is:ie+1,jsd:jed,k), yfx_adv(isd:ied,js:je+1,k), gridstruct, bd, ra_x, ra_y)
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
      call fv_tp_2d(zh(isd:ied,jsd:jed,k), crx_adv(is:ie+1,jsd:jed,k), cry_adv(isd:ied,js:je+1,k), npx,  npy, hord, &
                    fx, fy, xfx_adv(is:ie+1,jsd:jed,k), yfx_adv(isd:ied,js:je+1,k), gridstruct, bd, ra_x, ra_y)
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
!  Differentiation of riem_solver_c in forward (tangent) mode:
!   variations   of useful results: gz pef
!   with respect to varying inputs: ws gz delp w3 pef pt
  SUBROUTINE RIEM_SOLVER_C_TLM(ms, dt, is, ie, js, je, km, ng, akap, &
&   cappa, cp, ptop, hs, w3, w3_tl, pt, pt_tl, q_con, delp, delp_tl, gz&
&   , gz_tl, pef, pef_tl, ws, ws_tl, p_fac, a_imp, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, km
    INTEGER, INTENT(IN) :: ms
    REAL, INTENT(IN) :: dt, akap, cp, ptop, p_fac, a_imp, scale_m
    REAL, INTENT(IN) :: ws(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(IN) :: ws_tl(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: pt, &
&   delp
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: pt_tl, &
&   delp_tl
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: q_con, &
&   cappa
    REAL, INTENT(IN) :: hs(is-ng:ie+ng, js-ng:je+ng)
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: w3
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km), INTENT(IN) :: w3_tl
! OUTPUT PARAMETERS
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1), INTENT(INOUT) :: gz
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1), INTENT(INOUT) :: &
&   gz_tl
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1), INTENT(OUT) :: pef
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, km+1), INTENT(OUT) :: &
&   pef_tl
! Local:
    REAL, DIMENSION(is-1:ie+1, km) :: dm, dz2, w2, pm2, gm2, cp2
    REAL, DIMENSION(is-1:ie+1, km) :: dm_tl, dz2_tl, w2_tl, pm2_tl
    REAL, DIMENSION(is-1:ie+1, km+1) :: pem, pe2, peg
    REAL, DIMENSION(is-1:ie+1, km+1) :: pem_tl, pe2_tl
    REAL :: gama, rgrav
    INTEGER :: i, j, k
    INTEGER :: is1, ie1
    INTRINSIC LOG
    REAL :: arg1
    REAL :: arg1_tl
    gama = 1./(1.-akap)
    rgrav = 1./grav
    is1 = is - 1
    ie1 = ie + 1
    dm_tl = 0.0
    pe2_tl = 0.0
    dz2_tl = 0.0
    w2_tl = 0.0
    pm2_tl = 0.0
    pem_tl = 0.0
!$OMP parallel do default(none) shared(js,je,is1,ie1,km,delp,pef,ptop,gz,rgrav,w3,pt, &
!$OMP                                  a_imp,dt,gama,akap,ws,p_fac,scale_m,ms,hs,q_con,cappa) &
!$OMP                          private(cp2,gm2, dm, dz2, w2, pm2, pe2, pem, peg)
    DO j=js-1,je+1
      DO k=1,km
        DO i=is1,ie1
          dm_tl(i, k) = delp_tl(i, j, k)
          dm(i, k) = delp(i, j, k)
        END DO
      END DO
      DO i=is1,ie1
! full pressure at top
        pef_tl(i, j, 1) = 0.0
        pef(i, j, 1) = ptop
        pem_tl(i, 1) = 0.0
        pem(i, 1) = ptop
      END DO
      DO k=2,km+1
        DO i=is1,ie1
          pem_tl(i, k) = pem_tl(i, k-1) + dm_tl(i, k-1)
          pem(i, k) = pem(i, k-1) + dm(i, k-1)
        END DO
      END DO
      DO k=1,km
        DO i=is1,ie1
          dz2_tl(i, k) = gz_tl(i, j, k+1) - gz_tl(i, j, k)
          dz2(i, k) = gz(i, j, k+1) - gz(i, j, k)
          arg1_tl = (pem_tl(i, k+1)*pem(i, k)-pem(i, k+1)*pem_tl(i, k))/&
&           pem(i, k)**2
          arg1 = pem(i, k+1)/pem(i, k)
          pm2_tl(i, k) = (dm_tl(i, k)*LOG(arg1)-dm(i, k)*arg1_tl/arg1)/&
&           LOG(arg1)**2
          pm2(i, k) = dm(i, k)/LOG(arg1)
          dm_tl(i, k) = rgrav*dm_tl(i, k)
          dm(i, k) = dm(i, k)*rgrav
          w2_tl(i, k) = w3_tl(i, j, k)
          w2(i, k) = w3(i, j, k)
        END DO
      END DO
      IF (a_imp .LT. -0.01) THEN
        CALL SIM3P0_SOLVER_TLM(dt, is1, ie1, km, rdgas, gama, akap, pe2&
&                        , pe2_tl, dm, dm_tl, pem, pem_tl, w2, w2_tl, &
&                        dz2, dz2_tl, pt(is1:ie1, j, 1:km), pt_tl(is1:&
&                        ie1, j, 1:km), ws(is1:ie1, j), ws_tl(is1:ie1, j&
&                        ), p_fac, scale_m)
      ELSE IF (a_imp .LE. 0.5) THEN
        CALL RIM_2D_TLM(ms, dt, is1, ie1, km, rdgas, gama, gm2, pe2, &
&                 pe2_tl, dm, dm_tl, pm2, pm2_tl, w2, w2_tl, dz2, dz2_tl&
&                 , pt(is1:ie1, j, 1:km), pt_tl(is1:ie1, j, 1:km), ws(&
&                 is1:ie1, j), ws_tl(is1:ie1, j), .true.)
      ELSE
        CALL SIM1_SOLVER_TLM(dt, is1, ie1, km, rdgas, gama, gm2, cp2, &
&                      akap, pe2, pe2_tl, dm, dm_tl, pm2, pm2_tl, pem, &
&                      pem_tl, w2, w2_tl, dz2, dz2_tl, pt(is1:ie1, j, 1:&
&                      km), pt_tl(is1:ie1, j, 1:km), ws(is1:ie1, j), &
&                      ws_tl(is1:ie1, j), p_fac)
      END IF
      DO k=2,km+1
        DO i=is1,ie1
! add hydrostatic full-component
          pef_tl(i, j, k) = pe2_tl(i, k) + pem_tl(i, k)
          pef(i, j, k) = pe2(i, k) + pem(i, k)
        END DO
      END DO
! Compute Height * grav (for p-gradient computation)
      DO i=is1,ie1
        gz_tl(i, j, km+1) = 0.0
        gz(i, j, km+1) = hs(i, j)
      END DO
      DO k=km,1,-1
        DO i=is1,ie1
          gz_tl(i, j, k) = gz_tl(i, j, k+1) - grav*dz2_tl(i, k)
          gz(i, j, k) = gz(i, j, k+1) - dz2(i, k)*grav
        END DO
      END DO
    END DO
  END SUBROUTINE RIEM_SOLVER_C_TLM
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
    REAL :: arg1
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
          arg1 = pem(i, k+1)/pem(i, k)
          pm2(i, k) = dm(i, k)/LOG(arg1)
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
!  Differentiation of rim_2d in forward (tangent) mode:
!   variations   of useful results: pe2 dz2 w2
!   with respect to varying inputs: ws pe2 dm2 dz2 w2 pm2 pt2
  SUBROUTINE RIM_2D_TLM(ms, bdt, is, ie, km, rgas, gama, gm2, pe2, &
&   pe2_tl, dm2, dm2_tl, pm2, pm2_tl, w2, w2_tl, dz2, dz2_tl, pt2, &
&   pt2_tl, ws, ws_tl, c_core)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ms, is, ie, km
    REAL, INTENT(IN) :: bdt, gama, rgas
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pm2, gm2
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2_tl, pm2_tl
    LOGICAL, INTENT(IN) :: c_core
    REAL, INTENT(IN) :: pt2(is:ie, km)
    REAL, INTENT(IN) :: pt2_tl(is:ie, km)
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, INTENT(IN) :: ws_tl(is:ie)
! IN/OUT:
    REAL, INTENT(INOUT) :: dz2(is:ie, km)
    REAL, INTENT(INOUT) :: dz2_tl(is:ie, km)
    REAL, INTENT(INOUT) :: w2(is:ie, km)
    REAL, INTENT(INOUT) :: w2_tl(is:ie, km)
    REAL, INTENT(OUT) :: pe2(is:ie, km+1)
    REAL, INTENT(OUT) :: pe2_tl(is:ie, km+1)
! Local:
    REAL :: ws2(is:ie)
    REAL :: ws2_tl(is:ie)
    REAL, DIMENSION(km+1) :: m_bot, m_top, r_bot, r_top, pe1, pbar, wbar
    REAL, DIMENSION(km+1) :: m_bot_tl, m_top_tl, r_bot_tl, r_top_tl, &
&   pe1_tl, pbar_tl, wbar_tl
    REAL, DIMENSION(km) :: r_hi, r_lo, dz, wm, dm, dts
    REAL, DIMENSION(km) :: r_hi_tl, r_lo_tl, dz_tl, wm_tl, dm_tl, dts_tl
    REAL, DIMENSION(km) :: pf1, wc, cm, pp, pt1
    REAL, DIMENSION(km) :: pf1_tl, wc_tl, cm_tl, pp_tl, pt1_tl
    REAL :: dt, rdt, grg, z_frac, ptmp1, rden, pf, time_left
    REAL :: z_frac_tl, ptmp1_tl, rden_tl, pf_tl, time_left_tl
    REAL :: m_surf
    REAL :: m_surf_tl
    INTEGER :: i, k, n, ke, kt1, ktop
    INTEGER :: ks0, ks1
    INTRINSIC REAL
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC SQRT
    INTRINSIC MAX
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: result1
    REAL :: result1_tl
    grg = gama*rgas
    rdt = 1./bdt
    dt = bdt/REAL(ms)
    pbar(:) = 0.
    wbar(:) = 0.
    ws2_tl = 0.0
    DO i=is,ie
      ws2_tl(i) = 2.*ws_tl(i)
      ws2(i) = 2.*ws(i)
    END DO
    dm_tl = 0.0
    dts_tl = 0.0
    dz_tl = 0.0
    r_lo_tl = 0.0
    wbar_tl = 0.0
    pf1_tl = 0.0
    m_top_tl = 0.0
    r_top_tl = 0.0
    m_bot_tl = 0.0
    r_bot_tl = 0.0
    pt1_tl = 0.0
    cm_tl = 0.0
    pp_tl = 0.0
    wc_tl = 0.0
    pbar_tl = 0.0
    wm_tl = 0.0
    r_hi_tl = 0.0
! end i-loop
    DO 6000 i=is,ie
      DO k=1,km
        dz_tl(k) = dz2_tl(i, k)
        dz(k) = dz2(i, k)
        dm_tl(k) = dm2_tl(i, k)
        dm(k) = dm2(i, k)
        wm_tl(k) = w2_tl(i, k)*dm(k) + w2(i, k)*dm_tl(k)
        wm(k) = w2(i, k)*dm(k)
        pt1_tl(k) = pt2_tl(i, k)
        pt1(k) = pt2(i, k)
      END DO
      pe1(:) = 0.
      wbar_tl(km+1) = ws_tl(i)
      wbar(km+1) = ws(i)
      ks0 = 1
      IF (ms .GT. 1 .AND. ms .LT. 8) THEN
! Continuity of (pbar, wbar) is maintained
        DO k=1,km
          rden_tl = -((rgas*dm_tl(k)*dz(k)-rgas*dm(k)*dz_tl(k))/dz(k)**2&
&           )
          rden = -(rgas*dm(k)/dz(k))
          arg1_tl = gama*(rden_tl*pt1(k)+rden*pt1_tl(k))/(rden*pt1(k))
          arg1 = gama*LOG(rden*pt1(k))
          pf1_tl(k) = arg1_tl*EXP(arg1)
          pf1(k) = EXP(arg1)
          arg1_tl = (grg*pf1_tl(k)*rden-grg*pf1(k)*rden_tl)/rden**2
          arg1 = grg*pf1(k)/rden
          IF (arg1 .EQ. 0.0) THEN
            result1_tl = 0.0
          ELSE
            result1_tl = arg1_tl/(2.0*SQRT(arg1))
          END IF
          result1 = SQRT(arg1)
          dts_tl(k) = -((dz_tl(k)*result1-dz(k)*result1_tl)/result1**2)
          dts(k) = -(dz(k)/result1)
          IF (bdt .GT. dts(k)) GOTO 100
        END DO
        ks0 = km
        GOTO 222
 100    ks0 = k - 1
 222    IF (ks0 .EQ. 1) THEN
          pe1_tl = 0.0
        ELSE
          DO k=1,ks0
            cm_tl(k) = (dm_tl(k)*dts(k)-dm(k)*dts_tl(k))/dts(k)**2
            cm(k) = dm(k)/dts(k)
            wc_tl(k) = (wm_tl(k)*dts(k)-wm(k)*dts_tl(k))/dts(k)**2
            wc(k) = wm(k)/dts(k)
            pp_tl(k) = pf1_tl(k) - pm2_tl(i, k)
            pp(k) = pf1(k) - pm2(i, k)
          END DO
          wbar_tl(1) = ((wc_tl(1)+pp_tl(1))*cm(1)-(wc(1)+pp(1))*cm_tl(1)&
&           )/cm(1)**2
          wbar(1) = (wc(1)+pp(1))/cm(1)
          pe1_tl = 0.0
          DO k=2,ks0
            wbar_tl(k) = ((wc_tl(k-1)+wc_tl(k)+pp_tl(k)-pp_tl(k-1))*(cm(&
&             k-1)+cm(k))-(wc(k-1)+wc(k)+pp(k)-pp(k-1))*(cm_tl(k-1)+&
&             cm_tl(k)))/(cm(k-1)+cm(k))**2
            wbar(k) = (wc(k-1)+wc(k)+pp(k)-pp(k-1))/(cm(k-1)+cm(k))
            pbar_tl(k) = bdt*(cm_tl(k-1)*wbar(k)+cm(k-1)*wbar_tl(k)-&
&             wc_tl(k-1)+pp_tl(k-1))
            pbar(k) = bdt*(cm(k-1)*wbar(k)-wc(k-1)+pp(k-1))
            pe1_tl(k) = pbar_tl(k)
            pe1(k) = pbar(k)
          END DO
          IF (ks0 .EQ. km) THEN
            pbar_tl(km+1) = bdt*(cm_tl(km)*wbar(km+1)+cm(km)*wbar_tl(km+&
&             1)-wc_tl(km)+pp_tl(km))
            pbar(km+1) = bdt*(cm(km)*wbar(km+1)-wc(km)+pp(km))
            IF (c_core) THEN
              DO k=1,km
                dz2_tl(i, k) = dz_tl(k) + bdt*(wbar_tl(k+1)-wbar_tl(k))
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
              END DO
            ELSE
              DO k=1,km
                dz2_tl(i, k) = dz_tl(k) + bdt*(wbar_tl(k+1)-wbar_tl(k))
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
                w2_tl(i, k) = ((wm_tl(k)+pbar_tl(k+1)-pbar_tl(k))*dm(k)-&
&                 (wm(k)+pbar(k+1)-pbar(k))*dm_tl(k))/dm(k)**2
                w2(i, k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
              END DO
            END IF
            pe2_tl(i, 1) = 0.0
            pe2(i, 1) = 0.
            DO k=2,km+1
              pe2_tl(i, k) = rdt*pbar_tl(k)
              pe2(i, k) = pbar(k)*rdt
            END DO
            GOTO 6000
          ELSE
! next i
            IF (c_core) THEN
              DO k=1,ks0-1
                dz2_tl(i, k) = dz_tl(k) + bdt*(wbar_tl(k+1)-wbar_tl(k))
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
              END DO
            ELSE
              DO k=1,ks0-1
                dz2_tl(i, k) = dz_tl(k) + bdt*(wbar_tl(k+1)-wbar_tl(k))
                dz2(i, k) = dz(k) + bdt*(wbar(k+1)-wbar(k))
                w2_tl(i, k) = ((wm_tl(k)+pbar_tl(k+1)-pbar_tl(k))*dm(k)-&
&                 (wm(k)+pbar(k+1)-pbar(k))*dm_tl(k))/dm(k)**2
                w2(i, k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
              END DO
            END IF
            pbar_tl(ks0) = pbar_tl(ks0)/REAL(ms)
            pbar(ks0) = pbar(ks0)/REAL(ms)
          END IF
        END IF
      ELSE
        pe1_tl = 0.0
      END IF
      ks1 = ks0
      DO n=1,ms
        DO k=ks1,km
          rden_tl = -((rgas*dm_tl(k)*dz(k)-rgas*dm(k)*dz_tl(k))/dz(k)**2&
&           )
          rden = -(rgas*dm(k)/dz(k))
          arg1_tl = gama*(rden_tl*pt1(k)+rden*pt1_tl(k))/(rden*pt1(k))
          arg1 = gama*LOG(rden*pt1(k))
          pf_tl = arg1_tl*EXP(arg1)
          pf = EXP(arg1)
          arg1_tl = (grg*pf_tl*rden-grg*pf*rden_tl)/rden**2
          arg1 = grg*pf/rden
          IF (arg1 .EQ. 0.0) THEN
            result1_tl = 0.0
          ELSE
            result1_tl = arg1_tl/(2.0*SQRT(arg1))
          END IF
          result1 = SQRT(arg1)
          dts_tl(k) = -((dz_tl(k)*result1-dz(k)*result1_tl)/result1**2)
          dts(k) = -(dz(k)/result1)
          ptmp1_tl = dts_tl(k)*(pf-pm2(i, k)) + dts(k)*(pf_tl-pm2_tl(i, &
&           k))
          ptmp1 = dts(k)*(pf-pm2(i, k))
          r_lo_tl(k) = wm_tl(k) + ptmp1_tl
          r_lo(k) = wm(k) + ptmp1
          r_hi_tl(k) = wm_tl(k) - ptmp1_tl
          r_hi(k) = wm(k) - ptmp1
        END DO
        ktop = ks1
        DO k=ks1,km
          IF (dt .GT. dts(k)) GOTO 110
        END DO
        ktop = km
        GOTO 333
 110    ktop = k - 1
 333    IF (ktop .GE. ks1) THEN
          DO k=ks1,ktop
            z_frac_tl = -(dt*dts_tl(k)/dts(k)**2)
            z_frac = dt/dts(k)
            r_bot_tl(k) = z_frac_tl*r_lo(k) + z_frac*r_lo_tl(k)
            r_bot(k) = z_frac*r_lo(k)
            r_top_tl(k+1) = z_frac_tl*r_hi(k) + z_frac*r_hi_tl(k)
            r_top(k+1) = z_frac*r_hi(k)
            m_bot_tl(k) = z_frac_tl*dm(k) + z_frac*dm_tl(k)
            m_bot(k) = z_frac*dm(k)
            m_top_tl(k+1) = m_bot_tl(k)
            m_top(k+1) = m_bot(k)
          END DO
          IF (ktop .EQ. km) GOTO 666
        END IF
        DO k=ktop+2,km+1
          m_top_tl(k) = 0.0
          m_top(k) = 0.
          r_top_tl(k) = 0.0
          r_top(k) = 0.
        END DO
        IF (1 .LT. ktop) THEN
          kt1 = ktop
        ELSE
          kt1 = 1
        END IF
        DO 444 ke=km+1,ktop+2,-1
          time_left = dt
          time_left_tl = 0.0
          DO k=ke-1,kt1,-1
            IF (time_left .GT. dts(k)) THEN
              time_left_tl = time_left_tl - dts_tl(k)
              time_left = time_left - dts(k)
              m_top_tl(ke) = m_top_tl(ke) + dm_tl(k)
              m_top(ke) = m_top(ke) + dm(k)
              r_top_tl(ke) = r_top_tl(ke) + r_hi_tl(k)
              r_top(ke) = r_top(ke) + r_hi(k)
            ELSE
              GOTO 120
            END IF
          END DO
          GOTO 444
 120      z_frac_tl = (time_left_tl*dts(k)-time_left*dts_tl(k))/dts(k)**&
&           2
          z_frac = time_left/dts(k)
          m_top_tl(ke) = m_top_tl(ke) + z_frac_tl*dm(k) + z_frac*dm_tl(k&
&           )
          m_top(ke) = m_top(ke) + z_frac*dm(k)
          r_top_tl(ke) = r_top_tl(ke) + z_frac_tl*r_hi(k) + z_frac*&
&           r_hi_tl(k)
          r_top(ke) = r_top(ke) + z_frac*r_hi(k)
 444    CONTINUE
! next level
        DO k=ktop+1,km
          m_bot_tl(k) = 0.0
          m_bot(k) = 0.
          r_bot_tl(k) = 0.0
          r_bot(k) = 0.
        END DO
        DO 4000 ke=ktop+1,km
          time_left = dt
          time_left_tl = 0.0
          DO k=ke,km
            IF (time_left .GT. dts(k)) THEN
              time_left_tl = time_left_tl - dts_tl(k)
              time_left = time_left - dts(k)
              m_bot_tl(ke) = m_bot_tl(ke) + dm_tl(k)
              m_bot(ke) = m_bot(ke) + dm(k)
              r_bot_tl(ke) = r_bot_tl(ke) + r_lo_tl(k)
              r_bot(ke) = r_bot(ke) + r_lo(k)
            ELSE
              GOTO 140
            END IF
          END DO
! next interface
          m_surf_tl = m_bot_tl(ke)
          m_surf = m_bot(ke)
          DO k=km,kt1,-1
            IF (time_left .GT. dts(k)) THEN
              time_left_tl = time_left_tl - dts_tl(k)
              time_left = time_left - dts(k)
              m_bot_tl(ke) = m_bot_tl(ke) + dm_tl(k)
              m_bot(ke) = m_bot(ke) + dm(k)
              r_bot_tl(ke) = r_bot_tl(ke) - r_hi_tl(k)
              r_bot(ke) = r_bot(ke) - r_hi(k)
            ELSE
              GOTO 130
            END IF
          END DO
          GOTO 4000
 130      z_frac_tl = (time_left_tl*dts(k)-time_left*dts_tl(k))/dts(k)**&
&           2
          z_frac = time_left/dts(k)
          m_bot_tl(ke) = m_bot_tl(ke) + z_frac_tl*dm(k) + z_frac*dm_tl(k&
&           )
          m_bot(ke) = m_bot(ke) + z_frac*dm(k)
          r_bot_tl(ke) = r_bot_tl(ke) - z_frac_tl*r_hi(k) - z_frac*&
&           r_hi_tl(k) + (m_bot_tl(ke)-m_surf_tl)*ws2(i) + (m_bot(ke)-&
&           m_surf)*ws2_tl(i)
          r_bot(ke) = r_bot(ke) - z_frac*r_hi(k) + (m_bot(ke)-m_surf)*&
&           ws2(i)
          GOTO 4000
 140      z_frac_tl = (time_left_tl*dts(k)-time_left*dts_tl(k))/dts(k)**&
&           2
          z_frac = time_left/dts(k)
          m_bot_tl(ke) = m_bot_tl(ke) + z_frac_tl*dm(k) + z_frac*dm_tl(k&
&           )
          m_bot(ke) = m_bot(ke) + z_frac*dm(k)
          r_bot_tl(ke) = r_bot_tl(ke) + z_frac_tl*r_lo(k) + z_frac*&
&           r_lo_tl(k)
          r_bot(ke) = r_bot(ke) + z_frac*r_lo(k)
 4000   CONTINUE
! next interface
 666    IF (ks1 .EQ. 1) THEN
          wbar_tl(1) = (r_bot_tl(1)*m_bot(1)-r_bot(1)*m_bot_tl(1))/m_bot&
&           (1)**2
          wbar(1) = r_bot(1)/m_bot(1)
        END IF
        DO k=ks1+1,km
          wbar_tl(k) = ((r_bot_tl(k)+r_top_tl(k))*(m_top(k)+m_bot(k))-(&
&           r_bot(k)+r_top(k))*(m_top_tl(k)+m_bot_tl(k)))/(m_top(k)+&
&           m_bot(k))**2
          wbar(k) = (r_bot(k)+r_top(k))/(m_top(k)+m_bot(k))
        END DO
! pbar here is actually dt*pbar
        DO k=ks1+1,km+1
          pbar_tl(k) = m_top_tl(k)*wbar(k) + m_top(k)*wbar_tl(k) - &
&           r_top_tl(k)
          pbar(k) = m_top(k)*wbar(k) - r_top(k)
          pe1_tl(k) = pe1_tl(k) + pbar_tl(k)
          pe1(k) = pe1(k) + pbar(k)
        END DO
        IF (n .EQ. ms) THEN
          IF (c_core) THEN
            DO k=ks1,km
              dz2_tl(i, k) = dz_tl(k) + dt*(wbar_tl(k+1)-wbar_tl(k))
              dz2(i, k) = dz(k) + dt*(wbar(k+1)-wbar(k))
            END DO
          ELSE
            DO k=ks1,km
              dz2_tl(i, k) = dz_tl(k) + dt*(wbar_tl(k+1)-wbar_tl(k))
              dz2(i, k) = dz(k) + dt*(wbar(k+1)-wbar(k))
              w2_tl(i, k) = ((wm_tl(k)+pbar_tl(k+1)-pbar_tl(k))*dm(k)-(&
&               wm(k)+pbar(k+1)-pbar(k))*dm_tl(k))/dm(k)**2
              w2(i, k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
            END DO
          END IF
        ELSE
          DO k=ks1,km
            dz_tl(k) = dz_tl(k) + dt*(wbar_tl(k+1)-wbar_tl(k))
            dz(k) = dz(k) + dt*(wbar(k+1)-wbar(k))
            wm_tl(k) = wm_tl(k) + pbar_tl(k+1) - pbar_tl(k)
            wm(k) = wm(k) + pbar(k+1) - pbar(k)
          END DO
        END IF
      END DO
      pe2_tl(i, 1) = 0.0
      pe2(i, 1) = 0.
      DO k=2,km+1
        pe2_tl(i, k) = rdt*pe1_tl(k)
        pe2(i, k) = pe1(k)*rdt
      END DO
 6000 CONTINUE
  END SUBROUTINE RIM_2D_TLM
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
    REAL :: arg1
    REAL :: result1
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
          arg1 = gama*LOG(rden*pt1(k))
          pf1(k) = EXP(arg1)
          arg1 = grg*pf1(k)/rden
          result1 = SQRT(arg1)
          dts(k) = -(dz(k)/result1)
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
          arg1 = gama*LOG(rden*pt1(k))
          pf = EXP(arg1)
          arg1 = grg*pf/rden
          result1 = SQRT(arg1)
          dts(k) = -(dz(k)/result1)
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
!  Differentiation of sim3_solver in forward (tangent) mode:
!   variations   of useful results: pe2 dz2 w2
!   with respect to varying inputs: ws dm pe2 dz2 w2 pem pt2
  SUBROUTINE SIM3_SOLVER_TLM(dt, is, ie, km, rgas, gama, kappa, pe2, &
&   pe2_tl, dm, dm_tl, pem, pem_tl, w2, w2_tl, dz2, dz2_tl, pt2, pt2_tl&
&   , ws, ws_tl, alpha, p_fac, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, alpha, p_fac, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm, pt2
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm_tl, pt2_tl
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, INTENT(IN) :: ws_tl(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem_tl
    REAL, INTENT(OUT) :: pe2(is:ie, km+1)
    REAL, INTENT(OUT) :: pe2_tl(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2_tl, w2_tl
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, wk, g_rat, gam
    REAL, DIMENSION(is:ie, km) :: aa_tl, bb_tl, dd_tl, w1_tl, wk_tl, &
&   g_rat_tl, gam_tl
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie, km+1) :: pp_tl
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL, DIMENSION(is:ie) :: p1_tl, wk1_tl, bet_tl
    REAL :: beta, t2, t1g, rdt, ra, capa1, r2g, r6g
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: arg2
    REAL :: arg2_tl
    beta = 1. - alpha
    ra = 1./alpha
    t2 = beta/alpha
    t1g = gama*2.*(alpha*dt)**2
    rdt = 1./dt
    capa1 = kappa - 1.
    r2g = grav/2.
    r6g = grav/6.
    aa_tl = 0.0
    w1_tl = 0.0
    DO k=1,km
      DO i=is,ie
        w1_tl(i, k) = w2_tl(i, k)
        w1(i, k) = w2(i, k)
! Full pressure at center
        arg1_tl = -(rgas*((dm_tl(i, k)*dz2(i, k)-dm(i, k)*dz2_tl(i, k))*&
&         pt2(i, k)/dz2(i, k)**2+dm(i, k)*pt2_tl(i, k)/dz2(i, k)))
        arg1 = -(dm(i, k)/dz2(i, k)*rgas*pt2(i, k))
        arg2_tl = gama*arg1_tl/arg1
        arg2 = gama*LOG(arg1)
        aa_tl(i, k) = arg2_tl*EXP(arg2)
        aa(i, k) = EXP(arg2)
      END DO
    END DO
    dd_tl = 0.0
    bb_tl = 0.0
    g_rat_tl = 0.0
    DO k=1,km-1
      DO i=is,ie
! for profile reconstruction
        g_rat_tl(i, k) = (dm_tl(i, k)*dm(i, k+1)-dm(i, k)*dm_tl(i, k+1))&
&         /dm(i, k+1)**2
        g_rat(i, k) = dm(i, k)/dm(i, k+1)
        bb_tl(i, k) = 2.*g_rat_tl(i, k)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd_tl(i, k) = 3.*(aa_tl(i, k)+g_rat_tl(i, k)*aa(i, k+1)+g_rat(i&
&         , k)*aa_tl(i, k+1))
        dd(i, k) = 3.*(aa(i, k)+g_rat(i, k)*aa(i, k+1))
      END DO
    END DO
    bet_tl = 0.0
! pe2 is full p at edges
    DO i=is,ie
! Top:
      bet_tl(i) = bb_tl(i, 1)
      bet(i) = bb(i, 1)
      pe2_tl(i, 1) = pem_tl(i, 1)
      pe2(i, 1) = pem(i, 1)
      pe2_tl(i, 2) = ((dd_tl(i, 1)-pem_tl(i, 1))*bet(i)-(dd(i, 1)-pem(i&
&       , 1))*bet_tl(i))/bet(i)**2
      pe2(i, 2) = (dd(i, 1)-pem(i, 1))/bet(i)
! Bottom:
      bb_tl(i, km) = 0.0
      bb(i, km) = 2.
      dd_tl(i, km) = 3.*aa_tl(i, km) + r2g*dm_tl(i, km)
      dd(i, km) = 3.*aa(i, km) + r2g*dm(i, km)
    END DO
    gam_tl = 0.0
    DO k=2,km
      DO i=is,ie
        gam_tl(i, k) = (g_rat_tl(i, k-1)*bet(i)-g_rat(i, k-1)*bet_tl(i))&
&         /bet(i)**2
        gam(i, k) = g_rat(i, k-1)/bet(i)
        bet_tl(i) = bb_tl(i, k) - gam_tl(i, k)
        bet(i) = bb(i, k) - gam(i, k)
        pe2_tl(i, k+1) = ((dd_tl(i, k)-pe2_tl(i, k))*bet(i)-(dd(i, k)-&
&         pe2(i, k))*bet_tl(i))/bet(i)**2
        pe2(i, k+1) = (dd(i, k)-pe2(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        pe2_tl(i, k) = pe2_tl(i, k) - gam_tl(i, k)*pe2(i, k+1) - gam(i, &
&         k)*pe2_tl(i, k+1)
        pe2(i, k) = pe2(i, k) - gam(i, k)*pe2(i, k+1)
      END DO
    END DO
    pp_tl = 0.0
! done reconstruction of full:
! pp is pert. p at edges
    DO k=1,km+1
      DO i=is,ie
        pp_tl(i, k) = pe2_tl(i, k) - pem_tl(i, k)
        pp(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
    wk_tl = 0.0
    DO k=2,km
      DO i=is,ie
        aa_tl(i, k) = t1g*pe2_tl(i, k)/(dz2(i, k-1)+dz2(i, k)) - t1g*(&
&         dz2_tl(i, k-1)+dz2_tl(i, k))*pe2(i, k)/(dz2(i, k-1)+dz2(i, k))&
&         **2
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*pe2(i, k)
        wk_tl(i, k) = t2*(aa_tl(i, k)*(w1(i, k-1)-w1(i, k))+aa(i, k)*(&
&         w1_tl(i, k-1)-w1_tl(i, k)))
        wk(i, k) = t2*aa(i, k)*(w1(i, k-1)-w1(i, k))
        aa_tl(i, k) = aa_tl(i, k) - scale_m*dm_tl(i, 1)
        aa(i, k) = aa(i, k) - scale_m*dm(i, 1)
      END DO
    END DO
    DO i=is,ie
      bet_tl(i) = dm_tl(i, 1) - aa_tl(i, 2)
      bet(i) = dm(i, 1) - aa(i, 2)
      w2_tl(i, 1) = ((dm_tl(i, 1)*w1(i, 1)+dm(i, 1)*w1_tl(i, 1)+dt*pp_tl&
&       (i, 2)+wk_tl(i, 2))*bet(i)-(dm(i, 1)*w1(i, 1)+dt*pp(i, 2)+wk(i, &
&       2))*bet_tl(i))/bet(i)**2
      w2(i, 1) = (dm(i, 1)*w1(i, 1)+dt*pp(i, 2)+wk(i, 2))/bet(i)
    END DO
    DO k=2,km-1
      DO i=is,ie
        gam_tl(i, k) = (aa_tl(i, k)*bet(i)-aa(i, k)*bet_tl(i))/bet(i)**2
        gam(i, k) = aa(i, k)/bet(i)
        bet_tl(i) = dm_tl(i, k) - aa_tl(i, k) - aa_tl(i, k+1) - aa_tl(i&
&         , k)*gam(i, k) - aa(i, k)*gam_tl(i, k)
        bet(i) = dm(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        w2_tl(i, k) = ((dm_tl(i, k)*w1(i, k)+dm(i, k)*w1_tl(i, k)+dt*(&
&         pp_tl(i, k+1)-pp_tl(i, k))+wk_tl(i, k+1)-wk_tl(i, k)-aa_tl(i, &
&         k)*w2(i, k-1)-aa(i, k)*w2_tl(i, k-1))*bet(i)-(dm(i, k)*w1(i, k&
&         )+dt*(pp(i, k+1)-pp(i, k))+wk(i, k+1)-wk(i, k)-aa(i, k)*w2(i, &
&         k-1))*bet_tl(i))/bet(i)**2
        w2(i, k) = (dm(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))+wk(i, k+1&
&         )-wk(i, k)-aa(i, k)*w2(i, k-1))/bet(i)
      END DO
    END DO
    wk1_tl = 0.0
    DO i=is,ie
      wk1_tl(i) = t1g*pe2_tl(i, km+1)/dz2(i, km) - t1g*dz2_tl(i, km)*pe2&
&       (i, km+1)/dz2(i, km)**2
      wk1(i) = t1g/dz2(i, km)*pe2(i, km+1)
      gam_tl(i, km) = (aa_tl(i, km)*bet(i)-aa(i, km)*bet_tl(i))/bet(i)**&
&       2
      gam(i, km) = aa(i, km)/bet(i)
      bet_tl(i) = dm_tl(i, km) - aa_tl(i, km) - wk1_tl(i) - aa_tl(i, km)&
&       *gam(i, km) - aa(i, km)*gam_tl(i, km)
      bet(i) = dm(i, km) - (aa(i, km)+wk1(i)+aa(i, km)*gam(i, km))
      w2_tl(i, km) = ((dm_tl(i, km)*w1(i, km)+dm(i, km)*w1_tl(i, km)+dt*&
&       (pp_tl(i, km+1)-pp_tl(i, km))-wk_tl(i, km)+wk1_tl(i)*(t2*w1(i, &
&       km)-ra*ws(i))+wk1(i)*(t2*w1_tl(i, km)-ra*ws_tl(i))-aa_tl(i, km)*&
&       w2(i, km-1)-aa(i, km)*w2_tl(i, km-1))*bet(i)-(dm(i, km)*w1(i, km&
&       )+dt*(pp(i, km+1)-pp(i, km))-wk(i, km)+wk1(i)*(t2*w1(i, km)-ra*&
&       ws(i))-aa(i, km)*w2(i, km-1))*bet_tl(i))/bet(i)**2
      w2(i, km) = (dm(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk(i, &
&       km)+wk1(i)*(t2*w1(i, km)-ra*ws(i))-aa(i, km)*w2(i, km-1))/bet(i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        w2_tl(i, k) = w2_tl(i, k) - gam_tl(i, k+1)*w2(i, k+1) - gam(i, k&
&         +1)*w2_tl(i, k+1)
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
! pe2 is updated perturbation p at edges
    DO i=is,ie
      pe2_tl(i, 1) = 0.0
      pe2(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        pe2_tl(i, k+1) = pe2_tl(i, k) + ra*(rdt*(dm_tl(i, k)*(w2(i, k)-&
&         w1(i, k))+dm(i, k)*(w2_tl(i, k)-w1_tl(i, k)))-beta*(pp_tl(i, k&
&         +1)-pp_tl(i, k)))
        pe2(i, k+1) = pe2(i, k) + (dm(i, k)*(w2(i, k)-w1(i, k))*rdt-beta&
&         *(pp(i, k+1)-pp(i, k)))*ra
      END DO
    END DO
! Full non-hydro pressure at edges:
    DO i=is,ie
      pe2_tl(i, 1) = pem_tl(i, 1)
      pe2(i, 1) = pem(i, 1)
    END DO
    DO k=2,km+1
      DO i=is,ie
        IF (p_fac*pem(i, k) .LT. pe2(i, k) + pem(i, k)) THEN
          pe2_tl(i, k) = pe2_tl(i, k) + pem_tl(i, k)
          pe2(i, k) = pe2(i, k) + pem(i, k)
        ELSE
          pe2_tl(i, k) = p_fac*pem_tl(i, k)
          pe2(i, k) = p_fac*pem(i, k)
        END IF
      END DO
    END DO
    p1_tl = 0.0
    DO i=is,ie
! Recover cell-averaged pressure
      p1_tl(i) = r3*(pe2_tl(i, km)+2.*pe2_tl(i, km+1)) - r6g*dm_tl(i, km&
&       )
      p1(i) = (pe2(i, km)+2.*pe2(i, km+1))*r3 - r6g*dm(i, km)
      arg1_tl = capa1*p1_tl(i)/p1(i)
      arg1 = capa1*LOG(p1(i))
      dz2_tl(i, km) = -(rgas*((dm_tl(i, km)*pt2(i, km)+dm(i, km)*pt2_tl(&
&       i, km))*EXP(arg1)+dm(i, km)*pt2(i, km)*arg1_tl*EXP(arg1)))
      dz2(i, km) = -(dm(i, km)*rgas*pt2(i, km)*EXP(arg1))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1_tl(i) = r3*(pe2_tl(i, k)+bb_tl(i, k)*pe2(i, k+1)+bb(i, k)*&
&         pe2_tl(i, k+1)+g_rat_tl(i, k)*pe2(i, k+2)+g_rat(i, k)*pe2_tl(i&
&         , k+2)) - g_rat_tl(i, k)*p1(i) - g_rat(i, k)*p1_tl(i)
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        arg1_tl = capa1*p1_tl(i)/p1(i)
        arg1 = capa1*LOG(p1(i))
        dz2_tl(i, k) = -(rgas*((dm_tl(i, k)*pt2(i, k)+dm(i, k)*pt2_tl(i&
&         , k))*EXP(arg1)+dm(i, k)*pt2(i, k)*arg1_tl*EXP(arg1)))
        dz2(i, k) = -(dm(i, k)*rgas*pt2(i, k)*EXP(arg1))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        pe2_tl(i, k) = pe2_tl(i, k) - pem_tl(i, k)
        pe2(i, k) = pe2(i, k) - pem(i, k)
        pe2_tl(i, k) = pe2_tl(i, k) + beta*(pp_tl(i, k)-pe2_tl(i, k))
        pe2(i, k) = pe2(i, k) + beta*(pp(i, k)-pe2(i, k))
      END DO
    END DO
  END SUBROUTINE SIM3_SOLVER_TLM
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
    REAL :: arg1
    REAL :: arg2
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
        arg1 = -(dm(i, k)/dz2(i, k)*rgas*pt2(i, k))
        arg2 = gama*LOG(arg1)
        aa(i, k) = EXP(arg2)
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
      arg1 = capa1*LOG(p1(i))
      dz2(i, km) = -(dm(i, km)*rgas*pt2(i, km)*EXP(arg1))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        arg1 = capa1*LOG(p1(i))
        dz2(i, k) = -(dm(i, k)*rgas*pt2(i, k)*EXP(arg1))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        pe2(i, k) = pe2(i, k) - pem(i, k)
        pe2(i, k) = pe2(i, k) + beta*(pp(i, k)-pe2(i, k))
      END DO
    END DO
  END SUBROUTINE SIM3_SOLVER
!  Differentiation of sim3p0_solver in forward (tangent) mode:
!   variations   of useful results: pe2 dz2 w2
!   with respect to varying inputs: ws dm pe2 dz2 w2 pem pt2
  SUBROUTINE SIM3P0_SOLVER_TLM(dt, is, ie, km, rgas, gama, kappa, pe2, &
&   pe2_tl, dm, dm_tl, pem, pem_tl, w2, w2_tl, dz2, dz2_tl, pt2, pt2_tl&
&   , ws, ws_tl, p_fac, scale_m)
    IMPLICIT NONE
! Sa SIM3, but for beta==0
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm, pt2
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm_tl, pt2_tl
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, INTENT(IN) :: ws_tl(is:ie)
    REAL, INTENT(IN) :: pem(is:ie, km+1)
    REAL, INTENT(IN) :: pem_tl(is:ie, km+1)
    REAL, INTENT(OUT) :: pe2(is:ie, km+1)
    REAL, INTENT(OUT) :: pe2_tl(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2_tl, w2_tl
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, g_rat, gam
    REAL, DIMENSION(is:ie, km) :: aa_tl, bb_tl, dd_tl, w1_tl, g_rat_tl, &
&   gam_tl
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie, km+1) :: pp_tl
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL, DIMENSION(is:ie) :: p1_tl, wk1_tl, bet_tl
    REAL :: t1g, rdt, capa1, r2g, r6g
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: arg2
    REAL :: arg2_tl
    t1g = 2.*gama*dt**2
    rdt = 1./dt
    capa1 = kappa - 1.
    r2g = grav/2.
    r6g = grav/6.
    aa_tl = 0.0
    w1_tl = 0.0
    DO k=1,km
      DO i=is,ie
        w1_tl(i, k) = w2_tl(i, k)
        w1(i, k) = w2(i, k)
! Full pressure at center
        arg1_tl = -(rgas*((dm_tl(i, k)*dz2(i, k)-dm(i, k)*dz2_tl(i, k))*&
&         pt2(i, k)/dz2(i, k)**2+dm(i, k)*pt2_tl(i, k)/dz2(i, k)))
        arg1 = -(dm(i, k)/dz2(i, k)*rgas*pt2(i, k))
        arg2_tl = gama*arg1_tl/arg1
        arg2 = gama*LOG(arg1)
        aa_tl(i, k) = arg2_tl*EXP(arg2)
        aa(i, k) = EXP(arg2)
      END DO
    END DO
    dd_tl = 0.0
    bb_tl = 0.0
    g_rat_tl = 0.0
    DO k=1,km-1
      DO i=is,ie
! for profile reconstruction
        g_rat_tl(i, k) = (dm_tl(i, k)*dm(i, k+1)-dm(i, k)*dm_tl(i, k+1))&
&         /dm(i, k+1)**2
        g_rat(i, k) = dm(i, k)/dm(i, k+1)
        bb_tl(i, k) = 2.*g_rat_tl(i, k)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd_tl(i, k) = 3.*(aa_tl(i, k)+g_rat_tl(i, k)*aa(i, k+1)+g_rat(i&
&         , k)*aa_tl(i, k+1))
        dd(i, k) = 3.*(aa(i, k)+g_rat(i, k)*aa(i, k+1))
      END DO
    END DO
    bet_tl = 0.0
! pe2 is full p at edges
    DO i=is,ie
! Top:
      bet_tl(i) = bb_tl(i, 1)
      bet(i) = bb(i, 1)
      pe2_tl(i, 1) = pem_tl(i, 1)
      pe2(i, 1) = pem(i, 1)
      pe2_tl(i, 2) = ((dd_tl(i, 1)-pem_tl(i, 1))*bet(i)-(dd(i, 1)-pem(i&
&       , 1))*bet_tl(i))/bet(i)**2
      pe2(i, 2) = (dd(i, 1)-pem(i, 1))/bet(i)
! Bottom:
      bb_tl(i, km) = 0.0
      bb(i, km) = 2.
      dd_tl(i, km) = 3.*aa_tl(i, km) + r2g*dm_tl(i, km)
      dd(i, km) = 3.*aa(i, km) + r2g*dm(i, km)
    END DO
    gam_tl = 0.0
    DO k=2,km
      DO i=is,ie
        gam_tl(i, k) = (g_rat_tl(i, k-1)*bet(i)-g_rat(i, k-1)*bet_tl(i))&
&         /bet(i)**2
        gam(i, k) = g_rat(i, k-1)/bet(i)
        bet_tl(i) = bb_tl(i, k) - gam_tl(i, k)
        bet(i) = bb(i, k) - gam(i, k)
        pe2_tl(i, k+1) = ((dd_tl(i, k)-pe2_tl(i, k))*bet(i)-(dd(i, k)-&
&         pe2(i, k))*bet_tl(i))/bet(i)**2
        pe2(i, k+1) = (dd(i, k)-pe2(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        pe2_tl(i, k) = pe2_tl(i, k) - gam_tl(i, k)*pe2(i, k+1) - gam(i, &
&         k)*pe2_tl(i, k+1)
        pe2(i, k) = pe2(i, k) - gam(i, k)*pe2(i, k+1)
      END DO
    END DO
    pp_tl = 0.0
! done reconstruction of full:
! pp is pert. p at edges
    DO k=1,km+1
      DO i=is,ie
        pp_tl(i, k) = pe2_tl(i, k) - pem_tl(i, k)
        pp(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
    DO k=2,km
      DO i=is,ie
        aa_tl(i, k) = t1g*pe2_tl(i, k)/(dz2(i, k-1)+dz2(i, k)) - t1g*(&
&         dz2_tl(i, k-1)+dz2_tl(i, k))*pe2(i, k)/(dz2(i, k-1)+dz2(i, k))&
&         **2 - scale_m*dm_tl(i, 1)
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*pe2(i, k) - scale_m*dm(i&
&         , 1)
      END DO
    END DO
    DO i=is,ie
      bet_tl(i) = dm_tl(i, 1) - aa_tl(i, 2)
      bet(i) = dm(i, 1) - aa(i, 2)
      w2_tl(i, 1) = ((dm_tl(i, 1)*w1(i, 1)+dm(i, 1)*w1_tl(i, 1)+dt*pp_tl&
&       (i, 2))*bet(i)-(dm(i, 1)*w1(i, 1)+dt*pp(i, 2))*bet_tl(i))/bet(i)&
&       **2
      w2(i, 1) = (dm(i, 1)*w1(i, 1)+dt*pp(i, 2))/bet(i)
    END DO
    DO k=2,km-1
      DO i=is,ie
        gam_tl(i, k) = (aa_tl(i, k)*bet(i)-aa(i, k)*bet_tl(i))/bet(i)**2
        gam(i, k) = aa(i, k)/bet(i)
        bet_tl(i) = dm_tl(i, k) - aa_tl(i, k) - aa_tl(i, k+1) - aa_tl(i&
&         , k)*gam(i, k) - aa(i, k)*gam_tl(i, k)
        bet(i) = dm(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        w2_tl(i, k) = ((dm_tl(i, k)*w1(i, k)+dm(i, k)*w1_tl(i, k)+dt*(&
&         pp_tl(i, k+1)-pp_tl(i, k))-aa_tl(i, k)*w2(i, k-1)-aa(i, k)*&
&         w2_tl(i, k-1))*bet(i)-(dm(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, &
&         k))-aa(i, k)*w2(i, k-1))*bet_tl(i))/bet(i)**2
        w2(i, k) = (dm(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))-aa(i, k)*&
&         w2(i, k-1))/bet(i)
      END DO
    END DO
    wk1_tl = 0.0
    DO i=is,ie
      wk1_tl(i) = t1g*pe2_tl(i, km+1)/dz2(i, km) - t1g*dz2_tl(i, km)*pe2&
&       (i, km+1)/dz2(i, km)**2
      wk1(i) = t1g/dz2(i, km)*pe2(i, km+1)
      gam_tl(i, km) = (aa_tl(i, km)*bet(i)-aa(i, km)*bet_tl(i))/bet(i)**&
&       2
      gam(i, km) = aa(i, km)/bet(i)
      bet_tl(i) = dm_tl(i, km) - aa_tl(i, km) - wk1_tl(i) - aa_tl(i, km)&
&       *gam(i, km) - aa(i, km)*gam_tl(i, km)
      bet(i) = dm(i, km) - (aa(i, km)+wk1(i)+aa(i, km)*gam(i, km))
      w2_tl(i, km) = ((dm_tl(i, km)*w1(i, km)+dm(i, km)*w1_tl(i, km)+dt*&
&       (pp_tl(i, km+1)-pp_tl(i, km))-wk1_tl(i)*ws(i)-wk1(i)*ws_tl(i)-&
&       aa_tl(i, km)*w2(i, km-1)-aa(i, km)*w2_tl(i, km-1))*bet(i)-(dm(i&
&       , km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk1(i)*ws(i)-aa(i, km&
&       )*w2(i, km-1))*bet_tl(i))/bet(i)**2
      w2(i, km) = (dm(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk1(i)&
&       *ws(i)-aa(i, km)*w2(i, km-1))/bet(i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        w2_tl(i, k) = w2_tl(i, k) - gam_tl(i, k+1)*w2(i, k+1) - gam(i, k&
&         +1)*w2_tl(i, k+1)
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
! pe2 is updated perturbation p at edges
    DO i=is,ie
      pe2_tl(i, 1) = 0.0
      pe2(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        pe2_tl(i, k+1) = pe2_tl(i, k) + rdt*(dm_tl(i, k)*(w2(i, k)-w1(i&
&         , k))+dm(i, k)*(w2_tl(i, k)-w1_tl(i, k)))
        pe2(i, k+1) = pe2(i, k) + dm(i, k)*(w2(i, k)-w1(i, k))*rdt
      END DO
    END DO
! Full non-hydro pressure at edges:
    DO i=is,ie
      pe2_tl(i, 1) = pem_tl(i, 1)
      pe2(i, 1) = pem(i, 1)
    END DO
    DO k=2,km+1
      DO i=is,ie
        IF (p_fac*pem(i, k) .LT. pe2(i, k) + pem(i, k)) THEN
          pe2_tl(i, k) = pe2_tl(i, k) + pem_tl(i, k)
          pe2(i, k) = pe2(i, k) + pem(i, k)
        ELSE
          pe2_tl(i, k) = p_fac*pem_tl(i, k)
          pe2(i, k) = p_fac*pem(i, k)
        END IF
      END DO
    END DO
    p1_tl = 0.0
    DO i=is,ie
! Recover cell-averaged pressure
      p1_tl(i) = r3*(pe2_tl(i, km)+2.*pe2_tl(i, km+1)) - r6g*dm_tl(i, km&
&       )
      p1(i) = (pe2(i, km)+2.*pe2(i, km+1))*r3 - r6g*dm(i, km)
      arg1_tl = capa1*p1_tl(i)/p1(i)
      arg1 = capa1*LOG(p1(i))
      dz2_tl(i, km) = -(rgas*((dm_tl(i, km)*pt2(i, km)+dm(i, km)*pt2_tl(&
&       i, km))*EXP(arg1)+dm(i, km)*pt2(i, km)*arg1_tl*EXP(arg1)))
      dz2(i, km) = -(dm(i, km)*rgas*pt2(i, km)*EXP(arg1))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1_tl(i) = r3*(pe2_tl(i, k)+bb_tl(i, k)*pe2(i, k+1)+bb(i, k)*&
&         pe2_tl(i, k+1)+g_rat_tl(i, k)*pe2(i, k+2)+g_rat(i, k)*pe2_tl(i&
&         , k+2)) - g_rat_tl(i, k)*p1(i) - g_rat(i, k)*p1_tl(i)
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        arg1_tl = capa1*p1_tl(i)/p1(i)
        arg1 = capa1*LOG(p1(i))
        dz2_tl(i, k) = -(rgas*((dm_tl(i, k)*pt2(i, k)+dm(i, k)*pt2_tl(i&
&         , k))*EXP(arg1)+dm(i, k)*pt2(i, k)*arg1_tl*EXP(arg1)))
        dz2(i, k) = -(dm(i, k)*rgas*pt2(i, k)*EXP(arg1))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        pe2_tl(i, k) = pe2_tl(i, k) - pem_tl(i, k)
        pe2(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
  END SUBROUTINE SIM3P0_SOLVER_TLM
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
    REAL :: arg1
    REAL :: arg2
    t1g = 2.*gama*dt**2
    rdt = 1./dt
    capa1 = kappa - 1.
    r2g = grav/2.
    r6g = grav/6.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
! Full pressure at center
        arg1 = -(dm(i, k)/dz2(i, k)*rgas*pt2(i, k))
        arg2 = gama*LOG(arg1)
        aa(i, k) = EXP(arg2)
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
      arg1 = capa1*LOG(p1(i))
      dz2(i, km) = -(dm(i, km)*rgas*pt2(i, km)*EXP(arg1))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        arg1 = capa1*LOG(p1(i))
        dz2(i, k) = -(dm(i, k)*rgas*pt2(i, k)*EXP(arg1))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        pe2(i, k) = pe2(i, k) - pem(i, k)
      END DO
    END DO
  END SUBROUTINE SIM3P0_SOLVER
!  Differentiation of sim1_solver in forward (tangent) mode:
!   variations   of useful results: dz2 w2 pe
!   with respect to varying inputs: ws dm2 dz2 w2 pm2 pem pe pt2
  SUBROUTINE SIM1_SOLVER_TLM(dt, is, ie, km, rgas, gama, gm2, cp2, kappa&
&   , pe, pe_tl, dm2, dm2_tl, pm2, pm2_tl, pem, pem_tl, w2, w2_tl, dz2, &
&   dz2_tl, pt2, pt2_tl, ws, ws_tl, p_fac)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pt2, pm2, gm2, cp2
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2_tl, pt2_tl, pm2_tl
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, INTENT(IN) :: ws_tl(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem_tl
    REAL, INTENT(OUT) :: pe(is:ie, km+1)
    REAL, INTENT(OUT) :: pe_tl(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2_tl, w2_tl
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, g_rat, gam
    REAL, DIMENSION(is:ie, km) :: aa_tl, bb_tl, dd_tl, w1_tl, g_rat_tl, &
&   gam_tl
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie, km+1) :: pp_tl
    REAL, DIMENSION(is:ie) :: p1, bet
    REAL, DIMENSION(is:ie) :: p1_tl, bet_tl
    REAL :: t1g, rdt, capa1
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: max1
    REAL :: max1_tl
    REAL :: max2
    REAL :: max2_tl
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: arg2
    REAL :: arg2_tl
    t1g = gama*2.*dt*dt
    rdt = 1./dt
    capa1 = kappa - 1.
    w1_tl = 0.0
    DO k=1,km
      DO i=is,ie
        w1_tl(i, k) = w2_tl(i, k)
        w1(i, k) = w2(i, k)
        arg1_tl = -(rgas*((dm2_tl(i, k)*dz2(i, k)-dm2(i, k)*dz2_tl(i, k)&
&         )*pt2(i, k)/dz2(i, k)**2+dm2(i, k)*pt2_tl(i, k)/dz2(i, k)))
        arg1 = -(dm2(i, k)/dz2(i, k)*rgas*pt2(i, k))
        arg2_tl = gama*arg1_tl/arg1
        arg2 = gama*LOG(arg1)
        pe_tl(i, k) = arg2_tl*EXP(arg2) - pm2_tl(i, k)
        pe(i, k) = EXP(arg2) - pm2(i, k)
      END DO
    END DO
    dd_tl = 0.0
    bb_tl = 0.0
    g_rat_tl = 0.0
    DO k=1,km-1
      DO i=is,ie
        g_rat_tl(i, k) = (dm2_tl(i, k)*dm2(i, k+1)-dm2(i, k)*dm2_tl(i, k&
&         +1))/dm2(i, k+1)**2
        g_rat(i, k) = dm2(i, k)/dm2(i, k+1)
        bb_tl(i, k) = 2.*g_rat_tl(i, k)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd_tl(i, k) = 3.*(pe_tl(i, k)+g_rat_tl(i, k)*pe(i, k+1)+g_rat(i&
&         , k)*pe_tl(i, k+1))
        dd(i, k) = 3.*(pe(i, k)+g_rat(i, k)*pe(i, k+1))
      END DO
    END DO
    bet_tl = 0.0
    pp_tl = 0.0
    DO i=is,ie
      bet_tl(i) = bb_tl(i, 1)
      bet(i) = bb(i, 1)
      pp_tl(i, 1) = 0.0
      pp(i, 1) = 0.
      pp_tl(i, 2) = (dd_tl(i, 1)*bet(i)-dd(i, 1)*bet_tl(i))/bet(i)**2
      pp(i, 2) = dd(i, 1)/bet(i)
      bb_tl(i, km) = 0.0
      bb(i, km) = 2.
      dd_tl(i, km) = 3.*pe_tl(i, km)
      dd(i, km) = 3.*pe(i, km)
    END DO
    gam_tl = 0.0
    DO k=2,km
      DO i=is,ie
        gam_tl(i, k) = (g_rat_tl(i, k-1)*bet(i)-g_rat(i, k-1)*bet_tl(i))&
&         /bet(i)**2
        gam(i, k) = g_rat(i, k-1)/bet(i)
        bet_tl(i) = bb_tl(i, k) - gam_tl(i, k)
        bet(i) = bb(i, k) - gam(i, k)
        pp_tl(i, k+1) = ((dd_tl(i, k)-pp_tl(i, k))*bet(i)-(dd(i, k)-pp(i&
&         , k))*bet_tl(i))/bet(i)**2
        pp(i, k+1) = (dd(i, k)-pp(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        pp_tl(i, k) = pp_tl(i, k) - gam_tl(i, k)*pp(i, k+1) - gam(i, k)*&
&         pp_tl(i, k+1)
        pp(i, k) = pp(i, k) - gam(i, k)*pp(i, k+1)
      END DO
    END DO
    aa_tl = 0.0
! Start the w-solver
    DO k=2,km
      DO i=is,ie
        aa_tl(i, k) = t1g*(pem_tl(i, k)+pp_tl(i, k))/(dz2(i, k-1)+dz2(i&
&         , k)) - t1g*(dz2_tl(i, k-1)+dz2_tl(i, k))*(pem(i, k)+pp(i, k))&
&         /(dz2(i, k-1)+dz2(i, k))**2
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*(pem(i, k)+pp(i, k))
      END DO
    END DO
    DO i=is,ie
      bet_tl(i) = dm2_tl(i, 1) - aa_tl(i, 2)
      bet(i) = dm2(i, 1) - aa(i, 2)
      w2_tl(i, 1) = ((dm2_tl(i, 1)*w1(i, 1)+dm2(i, 1)*w1_tl(i, 1)+dt*&
&       pp_tl(i, 2))*bet(i)-(dm2(i, 1)*w1(i, 1)+dt*pp(i, 2))*bet_tl(i))/&
&       bet(i)**2
      w2(i, 1) = (dm2(i, 1)*w1(i, 1)+dt*pp(i, 2))/bet(i)
    END DO
    DO k=2,km-1
      DO i=is,ie
        gam_tl(i, k) = (aa_tl(i, k)*bet(i)-aa(i, k)*bet_tl(i))/bet(i)**2
        gam(i, k) = aa(i, k)/bet(i)
        bet_tl(i) = dm2_tl(i, k) - aa_tl(i, k) - aa_tl(i, k+1) - aa_tl(i&
&         , k)*gam(i, k) - aa(i, k)*gam_tl(i, k)
        bet(i) = dm2(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        w2_tl(i, k) = ((dm2_tl(i, k)*w1(i, k)+dm2(i, k)*w1_tl(i, k)+dt*(&
&         pp_tl(i, k+1)-pp_tl(i, k))-aa_tl(i, k)*w2(i, k-1)-aa(i, k)*&
&         w2_tl(i, k-1))*bet(i)-(dm2(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i&
&         , k))-aa(i, k)*w2(i, k-1))*bet_tl(i))/bet(i)**2
        w2(i, k) = (dm2(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))-aa(i, k)&
&         *w2(i, k-1))/bet(i)
      END DO
    END DO
    p1_tl = 0.0
    DO i=is,ie
      p1_tl(i) = t1g*(pem_tl(i, km+1)+pp_tl(i, km+1))/dz2(i, km) - t1g*&
&       dz2_tl(i, km)*(pem(i, km+1)+pp(i, km+1))/dz2(i, km)**2
      p1(i) = t1g/dz2(i, km)*(pem(i, km+1)+pp(i, km+1))
      gam_tl(i, km) = (aa_tl(i, km)*bet(i)-aa(i, km)*bet_tl(i))/bet(i)**&
&       2
      gam(i, km) = aa(i, km)/bet(i)
      bet_tl(i) = dm2_tl(i, km) - aa_tl(i, km) - p1_tl(i) - aa_tl(i, km)&
&       *gam(i, km) - aa(i, km)*gam_tl(i, km)
      bet(i) = dm2(i, km) - (aa(i, km)+p1(i)+aa(i, km)*gam(i, km))
      w2_tl(i, km) = ((dm2_tl(i, km)*w1(i, km)+dm2(i, km)*w1_tl(i, km)+&
&       dt*(pp_tl(i, km+1)-pp_tl(i, km))-p1_tl(i)*ws(i)-p1(i)*ws_tl(i)-&
&       aa_tl(i, km)*w2(i, km-1)-aa(i, km)*w2_tl(i, km-1))*bet(i)-(dm2(i&
&       , km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-p1(i)*ws(i)-aa(i, km)&
&       *w2(i, km-1))*bet_tl(i))/bet(i)**2
      w2(i, km) = (dm2(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-p1(i)&
&       *ws(i)-aa(i, km)*w2(i, km-1))/bet(i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        w2_tl(i, k) = w2_tl(i, k) - gam_tl(i, k+1)*w2(i, k+1) - gam(i, k&
&         +1)*w2_tl(i, k+1)
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
    DO i=is,ie
      pe_tl(i, 1) = 0.0
      pe(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        pe_tl(i, k+1) = pe_tl(i, k) + rdt*(dm2_tl(i, k)*(w2(i, k)-w1(i, &
&         k))+dm2(i, k)*(w2_tl(i, k)-w1_tl(i, k)))
        pe(i, k+1) = pe(i, k) + dm2(i, k)*(w2(i, k)-w1(i, k))*rdt
      END DO
    END DO
    DO i=is,ie
      p1_tl(i) = r3*(pe_tl(i, km)+2.*pe_tl(i, km+1))
      p1(i) = (pe(i, km)+2.*pe(i, km+1))*r3
      IF (p_fac*pm2(i, km) .LT. p1(i) + pm2(i, km)) THEN
        max1_tl = p1_tl(i) + pm2_tl(i, km)
        max1 = p1(i) + pm2(i, km)
      ELSE
        max1_tl = p_fac*pm2_tl(i, km)
        max1 = p_fac*pm2(i, km)
      END IF
      arg1_tl = capa1*max1_tl/max1
      arg1 = capa1*LOG(max1)
      dz2_tl(i, km) = -(rgas*((dm2_tl(i, km)*pt2(i, km)+dm2(i, km)*&
&       pt2_tl(i, km))*EXP(arg1)+dm2(i, km)*pt2(i, km)*arg1_tl*EXP(arg1)&
&       ))
      dz2(i, km) = -(dm2(i, km)*rgas*pt2(i, km)*EXP(arg1))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1_tl(i) = r3*(pe_tl(i, k)+bb_tl(i, k)*pe(i, k+1)+bb(i, k)*pe_tl&
&         (i, k+1)+g_rat_tl(i, k)*pe(i, k+2)+g_rat(i, k)*pe_tl(i, k+2)) &
&         - g_rat_tl(i, k)*p1(i) - g_rat(i, k)*p1_tl(i)
        p1(i) = (pe(i, k)+bb(i, k)*pe(i, k+1)+g_rat(i, k)*pe(i, k+2))*r3&
&         - g_rat(i, k)*p1(i)
        IF (p_fac*pm2(i, k) .LT. p1(i) + pm2(i, k)) THEN
          max2_tl = p1_tl(i) + pm2_tl(i, k)
          max2 = p1(i) + pm2(i, k)
        ELSE
          max2_tl = p_fac*pm2_tl(i, k)
          max2 = p_fac*pm2(i, k)
        END IF
        arg1_tl = capa1*max2_tl/max2
        arg1 = capa1*LOG(max2)
        dz2_tl(i, k) = -(rgas*((dm2_tl(i, k)*pt2(i, k)+dm2(i, k)*pt2_tl(&
&         i, k))*EXP(arg1)+dm2(i, k)*pt2(i, k)*arg1_tl*EXP(arg1)))
        dz2(i, k) = -(dm2(i, k)*rgas*pt2(i, k)*EXP(arg1))
      END DO
    END DO
  END SUBROUTINE SIM1_SOLVER_TLM
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
    REAL :: arg1
    REAL :: arg2
    t1g = gama*2.*dt*dt
    rdt = 1./dt
    capa1 = kappa - 1.
    DO k=1,km
      DO i=is,ie
        w1(i, k) = w2(i, k)
        arg1 = -(dm2(i, k)/dz2(i, k)*rgas*pt2(i, k))
        arg2 = gama*LOG(arg1)
        pe(i, k) = EXP(arg2) - pm2(i, k)
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
      arg1 = capa1*LOG(max1)
      dz2(i, km) = -(dm2(i, km)*rgas*pt2(i, km)*EXP(arg1))
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
        arg1 = capa1*LOG(max2)
        dz2(i, k) = -(dm2(i, k)*rgas*pt2(i, k)*EXP(arg1))
      END DO
    END DO
  END SUBROUTINE SIM1_SOLVER
!  Differentiation of sim_solver in forward (tangent) mode:
!   variations   of useful results: pe2 dz2 w2
!   with respect to varying inputs: ws pe2 dm2 dz2 w2 pm2 pem pt2
  SUBROUTINE SIM_SOLVER_TLM(dt, is, ie, km, rgas, gama, gm2, cp2, kappa&
&   , pe2, pe2_tl, dm2, dm2_tl, pm2, pm2_tl, pem, pem_tl, w2, w2_tl, dz2&
&   , dz2_tl, pt2, pt2_tl, ws, ws_tl, alpha, p_fac, scale_m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, km
    REAL, INTENT(IN) :: dt, rgas, gama, kappa, p_fac, alpha, scale_m
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2, pt2, pm2, gm2, cp2
    REAL, DIMENSION(is:ie, km), INTENT(IN) :: dm2_tl, pt2_tl, pm2_tl
    REAL, INTENT(IN) :: ws(is:ie)
    REAL, INTENT(IN) :: ws_tl(is:ie)
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem
    REAL, DIMENSION(is:ie, km+1), INTENT(IN) :: pem_tl
    REAL, INTENT(OUT) :: pe2(is:ie, km+1)
    REAL, INTENT(OUT) :: pe2_tl(is:ie, km+1)
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2, w2
    REAL, DIMENSION(is:ie, km), INTENT(INOUT) :: dz2_tl, w2_tl
! Local
    REAL, DIMENSION(is:ie, km) :: aa, bb, dd, w1, wk, g_rat, gam
    REAL, DIMENSION(is:ie, km) :: aa_tl, bb_tl, dd_tl, w1_tl, wk_tl, &
&   g_rat_tl, gam_tl
    REAL, DIMENSION(is:ie, km+1) :: pp
    REAL, DIMENSION(is:ie, km+1) :: pp_tl
    REAL, DIMENSION(is:ie) :: p1, wk1, bet
    REAL, DIMENSION(is:ie) :: p1_tl, wk1_tl, bet_tl
    REAL :: beta, t2, t1g, rdt, ra, capa1
    INTEGER :: i, k
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC MAX
    REAL :: max1
    REAL :: max1_tl
    REAL :: max2
    REAL :: max2_tl
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: arg2
    REAL :: arg2_tl
    beta = 1. - alpha
    ra = 1./alpha
    t2 = beta/alpha
    t1g = 2.*gama*(alpha*dt)**2
    rdt = 1./dt
    capa1 = kappa - 1.
    w1_tl = 0.0
    DO k=1,km
      DO i=is,ie
        w1_tl(i, k) = w2_tl(i, k)
        w1(i, k) = w2(i, k)
! P_g perturbation
        arg1_tl = -(rgas*((dm2_tl(i, k)*dz2(i, k)-dm2(i, k)*dz2_tl(i, k)&
&         )*pt2(i, k)/dz2(i, k)**2+dm2(i, k)*pt2_tl(i, k)/dz2(i, k)))
        arg1 = -(dm2(i, k)/dz2(i, k)*rgas*pt2(i, k))
        arg2_tl = gama*arg1_tl/arg1
        arg2 = gama*LOG(arg1)
        pe2_tl(i, k) = arg2_tl*EXP(arg2) - pm2_tl(i, k)
        pe2(i, k) = EXP(arg2) - pm2(i, k)
      END DO
    END DO
    dd_tl = 0.0
    bb_tl = 0.0
    g_rat_tl = 0.0
    DO k=1,km-1
      DO i=is,ie
        g_rat_tl(i, k) = (dm2_tl(i, k)*dm2(i, k+1)-dm2(i, k)*dm2_tl(i, k&
&         +1))/dm2(i, k+1)**2
        g_rat(i, k) = dm2(i, k)/dm2(i, k+1)
        bb_tl(i, k) = 2.*g_rat_tl(i, k)
        bb(i, k) = 2.*(1.+g_rat(i, k))
        dd_tl(i, k) = 3.*(pe2_tl(i, k)+g_rat_tl(i, k)*pe2(i, k+1)+g_rat(&
&         i, k)*pe2_tl(i, k+1))
        dd(i, k) = 3.*(pe2(i, k)+g_rat(i, k)*pe2(i, k+1))
      END DO
    END DO
    bet_tl = 0.0
    pp_tl = 0.0
    DO i=is,ie
      bet_tl(i) = bb_tl(i, 1)
      bet(i) = bb(i, 1)
      pp_tl(i, 1) = 0.0
      pp(i, 1) = 0.
      pp_tl(i, 2) = (dd_tl(i, 1)*bet(i)-dd(i, 1)*bet_tl(i))/bet(i)**2
      pp(i, 2) = dd(i, 1)/bet(i)
      bb_tl(i, km) = 0.0
      bb(i, km) = 2.
      dd_tl(i, km) = 3.*pe2_tl(i, km)
      dd(i, km) = 3.*pe2(i, km)
    END DO
    gam_tl = 0.0
    DO k=2,km
      DO i=is,ie
        gam_tl(i, k) = (g_rat_tl(i, k-1)*bet(i)-g_rat(i, k-1)*bet_tl(i))&
&         /bet(i)**2
        gam(i, k) = g_rat(i, k-1)/bet(i)
        bet_tl(i) = bb_tl(i, k) - gam_tl(i, k)
        bet(i) = bb(i, k) - gam(i, k)
        pp_tl(i, k+1) = ((dd_tl(i, k)-pp_tl(i, k))*bet(i)-(dd(i, k)-pp(i&
&         , k))*bet_tl(i))/bet(i)**2
        pp(i, k+1) = (dd(i, k)-pp(i, k))/bet(i)
      END DO
    END DO
    DO k=km,2,-1
      DO i=is,ie
        pp_tl(i, k) = pp_tl(i, k) - gam_tl(i, k)*pp(i, k+1) - gam(i, k)*&
&         pp_tl(i, k+1)
        pp(i, k) = pp(i, k) - gam(i, k)*pp(i, k+1)
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
! pe2 is Full p
        pe2_tl(i, k) = pem_tl(i, k) + pp_tl(i, k)
        pe2(i, k) = pem(i, k) + pp(i, k)
      END DO
    END DO
    aa_tl = 0.0
    wk_tl = 0.0
    DO k=2,km
      DO i=is,ie
        aa_tl(i, k) = t1g*pe2_tl(i, k)/(dz2(i, k-1)+dz2(i, k)) - t1g*(&
&         dz2_tl(i, k-1)+dz2_tl(i, k))*pe2(i, k)/(dz2(i, k-1)+dz2(i, k))&
&         **2
        aa(i, k) = t1g/(dz2(i, k-1)+dz2(i, k))*pe2(i, k)
        wk_tl(i, k) = t2*(aa_tl(i, k)*(w1(i, k-1)-w1(i, k))+aa(i, k)*(&
&         w1_tl(i, k-1)-w1_tl(i, k)))
        wk(i, k) = t2*aa(i, k)*(w1(i, k-1)-w1(i, k))
        aa_tl(i, k) = aa_tl(i, k) - scale_m*dm2_tl(i, 1)
        aa(i, k) = aa(i, k) - scale_m*dm2(i, 1)
      END DO
    END DO
! Top:
    DO i=is,ie
      bet_tl(i) = dm2_tl(i, 1) - aa_tl(i, 2)
      bet(i) = dm2(i, 1) - aa(i, 2)
      w2_tl(i, 1) = ((dm2_tl(i, 1)*w1(i, 1)+dm2(i, 1)*w1_tl(i, 1)+dt*&
&       pp_tl(i, 2)+wk_tl(i, 2))*bet(i)-(dm2(i, 1)*w1(i, 1)+dt*pp(i, 2)+&
&       wk(i, 2))*bet_tl(i))/bet(i)**2
      w2(i, 1) = (dm2(i, 1)*w1(i, 1)+dt*pp(i, 2)+wk(i, 2))/bet(i)
    END DO
! Interior:
    DO k=2,km-1
      DO i=is,ie
        gam_tl(i, k) = (aa_tl(i, k)*bet(i)-aa(i, k)*bet_tl(i))/bet(i)**2
        gam(i, k) = aa(i, k)/bet(i)
        bet_tl(i) = dm2_tl(i, k) - aa_tl(i, k) - aa_tl(i, k+1) - aa_tl(i&
&         , k)*gam(i, k) - aa(i, k)*gam_tl(i, k)
        bet(i) = dm2(i, k) - (aa(i, k)+aa(i, k+1)+aa(i, k)*gam(i, k))
        w2_tl(i, k) = ((dm2_tl(i, k)*w1(i, k)+dm2(i, k)*w1_tl(i, k)+dt*(&
&         pp_tl(i, k+1)-pp_tl(i, k))+wk_tl(i, k+1)-wk_tl(i, k)-aa_tl(i, &
&         k)*w2(i, k-1)-aa(i, k)*w2_tl(i, k-1))*bet(i)-(dm2(i, k)*w1(i, &
&         k)+dt*(pp(i, k+1)-pp(i, k))+wk(i, k+1)-wk(i, k)-aa(i, k)*w2(i&
&         , k-1))*bet_tl(i))/bet(i)**2
        w2(i, k) = (dm2(i, k)*w1(i, k)+dt*(pp(i, k+1)-pp(i, k))+wk(i, k+&
&         1)-wk(i, k)-aa(i, k)*w2(i, k-1))/bet(i)
      END DO
    END DO
    wk1_tl = 0.0
! Bottom: k=km
    DO i=is,ie
      wk1_tl(i) = t1g*pe2_tl(i, km+1)/dz2(i, km) - t1g*dz2_tl(i, km)*pe2&
&       (i, km+1)/dz2(i, km)**2
      wk1(i) = t1g/dz2(i, km)*pe2(i, km+1)
      gam_tl(i, km) = (aa_tl(i, km)*bet(i)-aa(i, km)*bet_tl(i))/bet(i)**&
&       2
      gam(i, km) = aa(i, km)/bet(i)
      bet_tl(i) = dm2_tl(i, km) - aa_tl(i, km) - wk1_tl(i) - aa_tl(i, km&
&       )*gam(i, km) - aa(i, km)*gam_tl(i, km)
      bet(i) = dm2(i, km) - (aa(i, km)+wk1(i)+aa(i, km)*gam(i, km))
      w2_tl(i, km) = ((dm2_tl(i, km)*w1(i, km)+dm2(i, km)*w1_tl(i, km)+&
&       dt*(pp_tl(i, km+1)-pp_tl(i, km))-wk_tl(i, km)+wk1_tl(i)*(t2*w1(i&
&       , km)-ra*ws(i))+wk1(i)*(t2*w1_tl(i, km)-ra*ws_tl(i))-aa_tl(i, km&
&       )*w2(i, km-1)-aa(i, km)*w2_tl(i, km-1))*bet(i)-(dm2(i, km)*w1(i&
&       , km)+dt*(pp(i, km+1)-pp(i, km))-wk(i, km)+wk1(i)*(t2*w1(i, km)-&
&       ra*ws(i))-aa(i, km)*w2(i, km-1))*bet_tl(i))/bet(i)**2
      w2(i, km) = (dm2(i, km)*w1(i, km)+dt*(pp(i, km+1)-pp(i, km))-wk(i&
&       , km)+wk1(i)*(t2*w1(i, km)-ra*ws(i))-aa(i, km)*w2(i, km-1))/bet(&
&       i)
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        w2_tl(i, k) = w2_tl(i, k) - gam_tl(i, k+1)*w2(i, k+1) - gam(i, k&
&         +1)*w2_tl(i, k+1)
        w2(i, k) = w2(i, k) - gam(i, k+1)*w2(i, k+1)
      END DO
    END DO
    DO i=is,ie
      pe2_tl(i, 1) = 0.0
      pe2(i, 1) = 0.
    END DO
    DO k=1,km
      DO i=is,ie
        pe2_tl(i, k+1) = pe2_tl(i, k) + ra*(rdt*(dm2_tl(i, k)*(w2(i, k)-&
&         w1(i, k))+dm2(i, k)*(w2_tl(i, k)-w1_tl(i, k)))-beta*(pp_tl(i, &
&         k+1)-pp_tl(i, k)))
        pe2(i, k+1) = pe2(i, k) + (dm2(i, k)*(w2(i, k)-w1(i, k))*rdt-&
&         beta*(pp(i, k+1)-pp(i, k)))*ra
      END DO
    END DO
    p1_tl = 0.0
    DO i=is,ie
      p1_tl(i) = r3*(pe2_tl(i, km)+2.*pe2_tl(i, km+1))
      p1(i) = (pe2(i, km)+2.*pe2(i, km+1))*r3
      IF (p_fac*pm2(i, km) .LT. p1(i) + pm2(i, km)) THEN
        max1_tl = p1_tl(i) + pm2_tl(i, km)
        max1 = p1(i) + pm2(i, km)
      ELSE
        max1_tl = p_fac*pm2_tl(i, km)
        max1 = p_fac*pm2(i, km)
      END IF
      arg1_tl = capa1*max1_tl/max1
      arg1 = capa1*LOG(max1)
      dz2_tl(i, km) = -(rgas*((dm2_tl(i, km)*pt2(i, km)+dm2(i, km)*&
&       pt2_tl(i, km))*EXP(arg1)+dm2(i, km)*pt2(i, km)*arg1_tl*EXP(arg1)&
&       ))
      dz2(i, km) = -(dm2(i, km)*rgas*pt2(i, km)*EXP(arg1))
    END DO
    DO k=km-1,1,-1
      DO i=is,ie
        p1_tl(i) = r3*(pe2_tl(i, k)+bb_tl(i, k)*pe2(i, k+1)+bb(i, k)*&
&         pe2_tl(i, k+1)+g_rat_tl(i, k)*pe2(i, k+2)+g_rat(i, k)*pe2_tl(i&
&         , k+2)) - g_rat_tl(i, k)*p1(i) - g_rat(i, k)*p1_tl(i)
        p1(i) = (pe2(i, k)+bb(i, k)*pe2(i, k+1)+g_rat(i, k)*pe2(i, k+2))&
&         *r3 - g_rat(i, k)*p1(i)
        IF (p_fac*pm2(i, k) .LT. p1(i) + pm2(i, k)) THEN
          max2_tl = p1_tl(i) + pm2_tl(i, k)
          max2 = p1(i) + pm2(i, k)
        ELSE
          max2_tl = p_fac*pm2_tl(i, k)
          max2 = p_fac*pm2(i, k)
        END IF
! delz = -dm*R*T_m / p_gas
        arg1_tl = capa1*max2_tl/max2
        arg1 = capa1*LOG(max2)
        dz2_tl(i, k) = -(rgas*((dm2_tl(i, k)*pt2(i, k)+dm2(i, k)*pt2_tl(&
&         i, k))*EXP(arg1)+dm2(i, k)*pt2(i, k)*arg1_tl*EXP(arg1)))
        dz2(i, k) = -(dm2(i, k)*rgas*pt2(i, k)*EXP(arg1))
      END DO
    END DO
    DO k=1,km+1
      DO i=is,ie
        pe2_tl(i, k) = pe2_tl(i, k) + beta*(pp_tl(i, k)-pe2_tl(i, k))
        pe2(i, k) = pe2(i, k) + beta*(pp(i, k)-pe2(i, k))
      END DO
    END DO
  END SUBROUTINE SIM_SOLVER_TLM
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
    REAL :: arg1
    REAL :: arg2
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
        arg1 = -(dm2(i, k)/dz2(i, k)*rgas*pt2(i, k))
        arg2 = gama*LOG(arg1)
        pe2(i, k) = EXP(arg2) - pm2(i, k)
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
      arg1 = capa1*LOG(max1)
      dz2(i, km) = -(dm2(i, km)*rgas*pt2(i, km)*EXP(arg1))
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
        arg1 = capa1*LOG(max2)
        dz2(i, k) = -(dm2(i, k)*rgas*pt2(i, k)*EXP(arg1))
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
!  Differentiation of edge_profile in forward (tangent) mode:
!   variations   of useful results: q1e q2e
!   with respect to varying inputs: q1e q2e q1 q2
  SUBROUTINE EDGE_PROFILE_TLM(q1, q1_tl, q2, q2_tl, q1e, q1e_tl, q2e, &
&   q2e_tl, i1, i2, j1, j2, j, km, dp0, uniform_grid, limiter)
    IMPLICIT NONE
! Optimized for wind profile reconstruction:
    INTEGER, INTENT(IN) :: i1, i2, j1, j2
    INTEGER, INTENT(IN) :: j, km
    INTEGER, INTENT(IN) :: limiter
    LOGICAL, INTENT(IN) :: uniform_grid
    REAL, INTENT(IN) :: dp0(km)
    REAL, DIMENSION(i1:i2, j1:j2, km), INTENT(IN) :: q1, q2
    REAL, DIMENSION(i1:i2, j1:j2, km), INTENT(IN) :: q1_tl, q2_tl
    REAL, DIMENSION(i1:i2, j1:j2, km+1), INTENT(OUT) :: q1e, q2e
    REAL, DIMENSION(i1:i2, j1:j2, km+1), INTENT(OUT) :: q1e_tl, q2e_tl
!-----------------------------------------------------------------------
! edge values
    REAL, DIMENSION(i1:i2, km+1) :: qe1, qe2, gam
    REAL, DIMENSION(i1:i2, km+1) :: qe1_tl, qe2_tl
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
      qe1_tl = 0.0
      qe2_tl = 0.0
      DO i=i1,i2
        qe1_tl(i, 1) = r4o3*q1_tl(i, j, 1) + r2o3*q1_tl(i, j, 2)
        qe1(i, 1) = r4o3*q1(i, j, 1) + r2o3*q1(i, j, 2)
        qe2_tl(i, 1) = r4o3*q2_tl(i, j, 1) + r2o3*q2_tl(i, j, 2)
        qe2(i, 1) = r4o3*q2(i, j, 1) + r2o3*q2(i, j, 2)
      END DO
      gak(1) = 7./3.
      DO k=2,km
        gak(k) = 1./(4.-gak(k-1))
        DO i=i1,i2
          qe1_tl(i, k) = gak(k)*(3.*(q1_tl(i, j, k-1)+q1_tl(i, j, k))-&
&           qe1_tl(i, k-1))
          qe1(i, k) = (3.*(q1(i, j, k-1)+q1(i, j, k))-qe1(i, k-1))*gak(k&
&           )
          qe2_tl(i, k) = gak(k)*(3.*(q2_tl(i, j, k-1)+q2_tl(i, j, k))-&
&           qe2_tl(i, k-1))
          qe2(i, k) = (3.*(q2(i, j, k-1)+q2(i, j, k))-qe2(i, k-1))*gak(k&
&           )
        END DO
      END DO
      bet = 1./(1.5-3.5*gak(km))
      DO i=i1,i2
        qe1_tl(i, km+1) = bet*(4.*q1_tl(i, j, km)+q1_tl(i, j, km-1)-3.5*&
&         qe1_tl(i, km))
        qe1(i, km+1) = (4.*q1(i, j, km)+q1(i, j, km-1)-3.5*qe1(i, km))*&
&         bet
        qe2_tl(i, km+1) = bet*(4.*q2_tl(i, j, km)+q2_tl(i, j, km-1)-3.5*&
&         qe2_tl(i, km))
        qe2(i, km+1) = (4.*q2(i, j, km)+q2(i, j, km-1)-3.5*qe2(i, km))*&
&         bet
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          qe1_tl(i, k) = qe1_tl(i, k) - gak(k)*qe1_tl(i, k+1)
          qe1(i, k) = qe1(i, k) - gak(k)*qe1(i, k+1)
          qe2_tl(i, k) = qe2_tl(i, k) - gak(k)*qe2_tl(i, k+1)
          qe2(i, k) = qe2(i, k) - gak(k)*qe2(i, k+1)
        END DO
      END DO
    ELSE
! Assuming grid varying in vertical only
      g0 = dp0(2)/dp0(1)
      xt1 = 2.*g0*(g0+1.)
      bet = g0*(g0+0.5)
      qe1_tl = 0.0
      qe2_tl = 0.0
      DO i=i1,i2
        qe1_tl(i, 1) = (xt1*q1_tl(i, j, 1)+q1_tl(i, j, 2))/bet
        qe1(i, 1) = (xt1*q1(i, j, 1)+q1(i, j, 2))/bet
        qe2_tl(i, 1) = (xt1*q2_tl(i, j, 1)+q2_tl(i, j, 2))/bet
        qe2(i, 1) = (xt1*q2(i, j, 1)+q2(i, j, 2))/bet
        gam(i, 1) = (1.+g0*(g0+1.5))/bet
      END DO
      DO k=2,km
        gk = dp0(k-1)/dp0(k)
        DO i=i1,i2
          bet = 2. + 2.*gk - gam(i, k-1)
          qe1_tl(i, k) = (3.*(q1_tl(i, j, k-1)+gk*q1_tl(i, j, k))-qe1_tl&
&           (i, k-1))/bet
          qe1(i, k) = (3.*(q1(i, j, k-1)+gk*q1(i, j, k))-qe1(i, k-1))/&
&           bet
          qe2_tl(i, k) = (3.*(q2_tl(i, j, k-1)+gk*q2_tl(i, j, k))-qe2_tl&
&           (i, k-1))/bet
          qe2(i, k) = (3.*(q2(i, j, k-1)+gk*q2(i, j, k))-qe2(i, k-1))/&
&           bet
          gam(i, k) = gk/bet
        END DO
      END DO
      a_bot = 1. + gk*(gk+1.5)
      xt1 = 2.*gk*(gk+1.)
      DO i=i1,i2
        xt2 = gk*(gk+0.5) - a_bot*gam(i, km)
        qe1_tl(i, km+1) = (xt1*q1_tl(i, j, km)+q1_tl(i, j, km-1)-a_bot*&
&         qe1_tl(i, km))/xt2
        qe1(i, km+1) = (xt1*q1(i, j, km)+q1(i, j, km-1)-a_bot*qe1(i, km)&
&         )/xt2
        qe2_tl(i, km+1) = (xt1*q2_tl(i, j, km)+q2_tl(i, j, km-1)-a_bot*&
&         qe2_tl(i, km))/xt2
        qe2(i, km+1) = (xt1*q2(i, j, km)+q2(i, j, km-1)-a_bot*qe2(i, km)&
&         )/xt2
      END DO
      DO k=km,1,-1
        DO i=i1,i2
          qe1_tl(i, k) = qe1_tl(i, k) - gam(i, k)*qe1_tl(i, k+1)
          qe1(i, k) = qe1(i, k) - gam(i, k)*qe1(i, k+1)
          qe2_tl(i, k) = qe2_tl(i, k) - gam(i, k)*qe2_tl(i, k+1)
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
        IF (q1(i, j, 1)*qe1(i, 1) .LT. 0.) THEN
          qe1_tl(i, 1) = 0.0
          qe1(i, 1) = 0.
        END IF
        IF (q2(i, j, 1)*qe2(i, 1) .LT. 0.) THEN
          qe2_tl(i, 1) = 0.0
          qe2(i, 1) = 0.
        END IF
! Surface:
        IF (q1(i, j, km)*qe1(i, km+1) .LT. 0.) THEN
          qe1_tl(i, km+1) = 0.0
          qe1(i, km+1) = 0.
        END IF
        IF (q2(i, j, km)*qe2(i, km+1) .LT. 0.) THEN
          qe2_tl(i, km+1) = 0.0
          qe2(i, km+1) = 0.
        END IF
      END DO
    END IF
    DO k=1,km+1
      DO i=i1,i2
        q1e_tl(i, j, k) = qe1_tl(i, k)
        q1e(i, j, k) = qe1(i, k)
        q2e_tl(i, j, k) = qe2_tl(i, k)
        q2e(i, j, k) = qe2(i, k)
      END DO
    END DO
  END SUBROUTINE EDGE_PROFILE_TLM
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
!  Differentiation of nest_halo_nh in forward (tangent) mode:
!   variations   of useful results: pk3 gz pkc
!   with respect to varying inputs: pk3 gz delp delz pkc pt
  SUBROUTINE NEST_HALO_NH_TLM(ptop, grav, kappa, cp, delp, delp_tl, delz&
&   , delz_tl, pt, pt_tl, phis, pkc, pkc_tl, gz, gz_tl, pk3, pk3_tl, npx&
&   , npy, npz, nested, pkc_pertn, computepk3, fullhalo, bd)
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
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz), INTENT(IN) :: &
&   pt_tl, delp_tl, delz_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(INOUT) &
&   :: gz, pkc, pk3
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1), INTENT(INOUT) &
&   :: gz_tl, pkc_tl, pk3_tl
    INTEGER :: i, j, k
!'gamma'
    REAL :: gama
    REAL :: ptk, rgrav, rkap, peln1, rdg
    REAL, DIMENSION(bd%isd:bd%ied, npz+1, bd%jsd:bd%jed) :: pe, peln
    REAL, DIMENSION(bd%isd:bd%ied, npz+1, bd%jsd:bd%jed) :: pe_tl, &
&   peln_tl
    REAL, DIMENSION(bd%isd:bd%ied, npz) :: gam, bb, dd, pkz
    REAL, DIMENSION(bd%isd:bd%ied, npz) :: gam_tl, bb_tl, dd_tl, pkz_tl
    REAL, DIMENSION(bd%isd:bd%ied, npz-1) :: g_rat
    REAL, DIMENSION(bd%isd:bd%ied, npz-1) :: g_rat_tl
    REAL, DIMENSION(bd%isd:bd%ied) :: bet
    REAL, DIMENSION(bd%isd:bd%ied) :: bet_tl
    REAL :: pm
    REAL :: pm_tl
    INTEGER :: ifirst, ilast, jfirst, jlast
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC LOG
    INTRINSIC EXP
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: arg2
    REAL :: arg2_tl
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
        dd_tl = 0.0
        peln_tl = 0.0
        bet_tl = 0.0
        bb_tl = 0.0
        pkz_tl = 0.0
        pe_tl = 0.0
        g_rat_tl = 0.0
        gam_tl = 0.0
        DO j=jfirst,jlast
!GZ
          DO i=ifirst,0
            gz_tl(i, j, npz+1) = 0.0
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=ifirst,0
              gz_tl(i, j, k) = gz_tl(i, j, k+1) - grav*delz_tl(i, j, k)
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=ifirst,0
            pe_tl(i, 1, j) = 0.0
            pe(i, 1, j) = ptop
            peln_tl(i, 1, j) = 0.0
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=ifirst,0
              pe_tl(i, k, j) = pe_tl(i, k-1, j) + delp_tl(i, j, k-1)
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              peln_tl(i, k, j) = pe_tl(i, k, j)/pe(i, k, j)
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=ifirst,0
!Full p
              arg1_tl = -(rdgas*((rgrav*delp_tl(i, j, k)*delz(i, j, k)-&
&               delp(i, j, k)*rgrav*delz_tl(i, j, k))*pt(i, j, k)/delz(i&
&               , j, k)**2+delp(i, j, k)*rgrav*pt_tl(i, j, k)/delz(i, j&
&               , k)))
              arg1 = -(delp(i, j, k)*rgrav/delz(i, j, k)*rdgas*pt(i, j, &
&               k))
              arg2_tl = gama*arg1_tl/arg1
              arg2 = gama*LOG(arg1)
              pkz_tl(i, k) = arg2_tl*EXP(arg2)
              pkz(i, k) = EXP(arg2)
!hydro
              pm_tl = (delp_tl(i, j, k)*(peln(i, k+1, j)-peln(i, k, j))-&
&               delp(i, j, k)*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(&
&               peln(i, k+1, j)-peln(i, k, j))**2
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz_tl(i, k) = pkz_tl(i, k) - pm_tl
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!pressure solver
          DO k=1,npz-1
            DO i=ifirst,0
              g_rat_tl(i, k) = (delp_tl(i, j, k)*delp(i, j, k+1)-delp(i&
&               , j, k)*delp_tl(i, j, k+1))/delp(i, j, k+1)**2
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb_tl(i, k) = 2.*g_rat_tl(i, k)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              dd_tl(i, k) = 3.*(pkz_tl(i, k)+g_rat_tl(i, k)*pkz(i, k+1)+&
&               g_rat(i, k)*pkz_tl(i, k+1))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=ifirst,0
            bet_tl(i) = bb_tl(i, 1)
            bet(i) = bb(i, 1)
            pkc_tl(i, j, 1) = 0.0
            pkc(i, j, 1) = 0.
            pkc_tl(i, j, 2) = (dd_tl(i, 1)*bet(i)-dd(i, 1)*bet_tl(i))/&
&             bet(i)**2
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb_tl(i, npz) = 0.0
            bb(i, npz) = 2.
            dd_tl(i, npz) = 3.*pkz_tl(i, npz)
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=ifirst,0
              gam_tl(i, k) = (g_rat_tl(i, k-1)*bet(i)-g_rat(i, k-1)*&
&               bet_tl(i))/bet(i)**2
              gam(i, k) = g_rat(i, k-1)/bet(i)
              bet_tl(i) = bb_tl(i, k) - gam_tl(i, k)
              bet(i) = bb(i, k) - gam(i, k)
              pkc_tl(i, j, k+1) = ((dd_tl(i, k)-pkc_tl(i, j, k))*bet(i)-&
&               (dd(i, k)-pkc(i, j, k))*bet_tl(i))/bet(i)**2
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ifirst,0
              pkc_tl(i, j, k) = pkc_tl(i, j, k) - gam_tl(i, k)*pkc(i, j&
&               , k+1) - gam(i, k)*pkc_tl(i, j, k+1)
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=jfirst,jlast
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=ifirst,0
                pkc_tl(i, j, k) = pkc_tl(i, j, k) + pe_tl(i, k, j)
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
          END IF
!pk3 if necessary; doesn't require condenstate loading calculation
          IF (computepk3) THEN
            DO i=ifirst,0
              pk3_tl(i, j, 1) = 0.0
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=ifirst,0
                arg1_tl = kappa*pe_tl(i, k, j)/pe(i, k, j)
                arg1 = kappa*LOG(pe(i, k, j))
                pk3_tl(i, j, k) = arg1_tl*EXP(arg1)
                pk3(i, j, k) = EXP(arg1)
              END DO
            END DO
          END IF
        END DO
      ELSE
        dd_tl = 0.0
        peln_tl = 0.0
        bet_tl = 0.0
        bb_tl = 0.0
        pkz_tl = 0.0
        pe_tl = 0.0
        g_rat_tl = 0.0
        gam_tl = 0.0
      END IF
      IF (ie .EQ. npx - 1) THEN
        DO j=jfirst,jlast
!GZ
          DO i=npx,ilast
            gz_tl(i, j, npz+1) = 0.0
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=npx,ilast
              gz_tl(i, j, k) = gz_tl(i, j, k+1) - grav*delz_tl(i, j, k)
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=npx,ilast
            pe_tl(i, 1, j) = 0.0
            pe(i, 1, j) = ptop
            peln_tl(i, 1, j) = 0.0
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=npx,ilast
              pe_tl(i, k, j) = pe_tl(i, k-1, j) + delp_tl(i, j, k-1)
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              peln_tl(i, k, j) = pe_tl(i, k, j)/pe(i, k, j)
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=npx,ilast
!Full p
              arg1_tl = -(rdgas*((rgrav*delp_tl(i, j, k)*delz(i, j, k)-&
&               delp(i, j, k)*rgrav*delz_tl(i, j, k))*pt(i, j, k)/delz(i&
&               , j, k)**2+delp(i, j, k)*rgrav*pt_tl(i, j, k)/delz(i, j&
&               , k)))
              arg1 = -(delp(i, j, k)*rgrav/delz(i, j, k)*rdgas*pt(i, j, &
&               k))
              arg2_tl = gama*arg1_tl/arg1
              arg2 = gama*LOG(arg1)
              pkz_tl(i, k) = arg2_tl*EXP(arg2)
              pkz(i, k) = EXP(arg2)
!hydro
              pm_tl = (delp_tl(i, j, k)*(peln(i, k+1, j)-peln(i, k, j))-&
&               delp(i, j, k)*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(&
&               peln(i, k+1, j)-peln(i, k, j))**2
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz_tl(i, k) = pkz_tl(i, k) - pm_tl
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!pressure solver
          DO k=1,npz-1
            DO i=npx,ilast
              g_rat_tl(i, k) = (delp_tl(i, j, k)*delp(i, j, k+1)-delp(i&
&               , j, k)*delp_tl(i, j, k+1))/delp(i, j, k+1)**2
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb_tl(i, k) = 2.*g_rat_tl(i, k)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              dd_tl(i, k) = 3.*(pkz_tl(i, k)+g_rat_tl(i, k)*pkz(i, k+1)+&
&               g_rat(i, k)*pkz_tl(i, k+1))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=npx,ilast
            bet_tl(i) = bb_tl(i, 1)
            bet(i) = bb(i, 1)
            pkc_tl(i, j, 1) = 0.0
            pkc(i, j, 1) = 0.
            pkc_tl(i, j, 2) = (dd_tl(i, 1)*bet(i)-dd(i, 1)*bet_tl(i))/&
&             bet(i)**2
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb_tl(i, npz) = 0.0
            bb(i, npz) = 2.
            dd_tl(i, npz) = 3.*pkz_tl(i, npz)
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=npx,ilast
              gam_tl(i, k) = (g_rat_tl(i, k-1)*bet(i)-g_rat(i, k-1)*&
&               bet_tl(i))/bet(i)**2
              gam(i, k) = g_rat(i, k-1)/bet(i)
              bet_tl(i) = bb_tl(i, k) - gam_tl(i, k)
              bet(i) = bb(i, k) - gam(i, k)
              pkc_tl(i, j, k+1) = ((dd_tl(i, k)-pkc_tl(i, j, k))*bet(i)-&
&               (dd(i, k)-pkc(i, j, k))*bet_tl(i))/bet(i)**2
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=npx,ilast
              pkc_tl(i, j, k) = pkc_tl(i, j, k) - gam_tl(i, k)*pkc(i, j&
&               , k+1) - gam(i, k)*pkc_tl(i, j, k+1)
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=jfirst,jlast
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=npx,ilast
                pkc_tl(i, j, k) = pkc_tl(i, j, k) + pe_tl(i, k, j)
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
          END IF
!pk3 if necessary
          IF (computepk3) THEN
            DO i=npx,ilast
              pk3_tl(i, j, 1) = 0.0
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=npx,ilast
                arg1_tl = kappa*pe_tl(i, k, j)/pe(i, k, j)
                arg1 = kappa*LOG(pe(i, k, j))
                pk3_tl(i, j, k) = arg1_tl*EXP(arg1)
                pk3(i, j, k) = EXP(arg1)
              END DO
            END DO
          END IF
        END DO
      END IF
      IF (js .EQ. 1) THEN
        DO j=jfirst,0
!GZ
          DO i=ifirst,ilast
            gz_tl(i, j, npz+1) = 0.0
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=ifirst,ilast
              gz_tl(i, j, k) = gz_tl(i, j, k+1) - grav*delz_tl(i, j, k)
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=ifirst,ilast
            pe_tl(i, 1, j) = 0.0
            pe(i, 1, j) = ptop
            peln_tl(i, 1, j) = 0.0
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=ifirst,ilast
              pe_tl(i, k, j) = pe_tl(i, k-1, j) + delp_tl(i, j, k-1)
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              peln_tl(i, k, j) = pe_tl(i, k, j)/pe(i, k, j)
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=ifirst,ilast
!Full p
              arg1_tl = -(rdgas*((rgrav*delp_tl(i, j, k)*delz(i, j, k)-&
&               delp(i, j, k)*rgrav*delz_tl(i, j, k))*pt(i, j, k)/delz(i&
&               , j, k)**2+delp(i, j, k)*rgrav*pt_tl(i, j, k)/delz(i, j&
&               , k)))
              arg1 = -(delp(i, j, k)*rgrav/delz(i, j, k)*rdgas*pt(i, j, &
&               k))
              arg2_tl = gama*arg1_tl/arg1
              arg2 = gama*LOG(arg1)
              pkz_tl(i, k) = arg2_tl*EXP(arg2)
              pkz(i, k) = EXP(arg2)
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!hydro
              pm_tl = (delp_tl(i, j, k)*(peln(i, k+1, j)-peln(i, k, j))-&
&               delp(i, j, k)*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(&
&               peln(i, k+1, j)-peln(i, k, j))**2
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz_tl(i, k) = pkz_tl(i, k) - pm_tl
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!pressure solver
          DO k=1,npz-1
            DO i=ifirst,ilast
              g_rat_tl(i, k) = (delp_tl(i, j, k)*delp(i, j, k+1)-delp(i&
&               , j, k)*delp_tl(i, j, k+1))/delp(i, j, k+1)**2
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb_tl(i, k) = 2.*g_rat_tl(i, k)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              dd_tl(i, k) = 3.*(pkz_tl(i, k)+g_rat_tl(i, k)*pkz(i, k+1)+&
&               g_rat(i, k)*pkz_tl(i, k+1))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=ifirst,ilast
            bet_tl(i) = bb_tl(i, 1)
            bet(i) = bb(i, 1)
            pkc_tl(i, j, 1) = 0.0
            pkc(i, j, 1) = 0.
            pkc_tl(i, j, 2) = (dd_tl(i, 1)*bet(i)-dd(i, 1)*bet_tl(i))/&
&             bet(i)**2
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb_tl(i, npz) = 0.0
            bb(i, npz) = 2.
            dd_tl(i, npz) = 3.*pkz_tl(i, npz)
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=ifirst,ilast
              gam_tl(i, k) = (g_rat_tl(i, k-1)*bet(i)-g_rat(i, k-1)*&
&               bet_tl(i))/bet(i)**2
              gam(i, k) = g_rat(i, k-1)/bet(i)
              bet_tl(i) = bb_tl(i, k) - gam_tl(i, k)
              bet(i) = bb(i, k) - gam(i, k)
              pkc_tl(i, j, k+1) = ((dd_tl(i, k)-pkc_tl(i, j, k))*bet(i)-&
&               (dd(i, k)-pkc(i, j, k))*bet_tl(i))/bet(i)**2
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ifirst,ilast
              pkc_tl(i, j, k) = pkc_tl(i, j, k) - gam_tl(i, k)*pkc(i, j&
&               , k+1) - gam(i, k)*pkc_tl(i, j, k+1)
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=jfirst,0
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=ifirst,ilast
                pkc_tl(i, j, k) = pkc_tl(i, j, k) + pe_tl(i, k, j)
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
          END IF
!pk3 if necessary
          IF (computepk3) THEN
            DO i=ifirst,ilast
              pk3_tl(i, j, 1) = 0.0
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=ifirst,ilast
                arg1_tl = kappa*pe_tl(i, k, j)/pe(i, k, j)
                arg1 = kappa*LOG(pe(i, k, j))
                pk3_tl(i, j, k) = arg1_tl*EXP(arg1)
                pk3(i, j, k) = EXP(arg1)
              END DO
            END DO
          END IF
        END DO
      END IF
      IF (je .EQ. npy - 1) THEN
        DO j=npy,jlast
!GZ
          DO i=ifirst,ilast
            gz_tl(i, j, npz+1) = 0.0
            gz(i, j, npz+1) = phis(i, j)
          END DO
          DO k=npz,1,-1
            DO i=ifirst,ilast
              gz_tl(i, j, k) = gz_tl(i, j, k+1) - grav*delz_tl(i, j, k)
              gz(i, j, k) = gz(i, j, k+1) - delz(i, j, k)*grav
            END DO
          END DO
!Hydrostatic interface pressure
          DO i=ifirst,ilast
            pe_tl(i, 1, j) = 0.0
            pe(i, 1, j) = ptop
            peln_tl(i, 1, j) = 0.0
            peln(i, 1, j) = peln1
          END DO
          DO k=2,npz+1
            DO i=ifirst,ilast
              pe_tl(i, k, j) = pe_tl(i, k-1, j) + delp_tl(i, j, k-1)
              pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
              peln_tl(i, k, j) = pe_tl(i, k, j)/pe(i, k, j)
              peln(i, k, j) = LOG(pe(i, k, j))
            END DO
          END DO
!Perturbation nonhydro layer-mean pressure (NOT to the kappa)
          DO k=1,npz
            DO i=ifirst,ilast
!Full p
              arg1_tl = -(rdgas*((rgrav*delp_tl(i, j, k)*delz(i, j, k)-&
&               delp(i, j, k)*rgrav*delz_tl(i, j, k))*pt(i, j, k)/delz(i&
&               , j, k)**2+delp(i, j, k)*rgrav*pt_tl(i, j, k)/delz(i, j&
&               , k)))
              arg1 = -(delp(i, j, k)*rgrav/delz(i, j, k)*rdgas*pt(i, j, &
&               k))
              arg2_tl = gama*arg1_tl/arg1
              arg2 = gama*LOG(arg1)
              pkz_tl(i, k) = arg2_tl*EXP(arg2)
              pkz(i, k) = EXP(arg2)
!hydro
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!hydro
              pm_tl = (delp_tl(i, j, k)*(peln(i, k+1, j)-peln(i, k, j))-&
&               delp(i, j, k)*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(&
&               peln(i, k+1, j)-peln(i, k, j))**2
              pm = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
!Remove hydro cell-mean pressure
              pkz_tl(i, k) = pkz_tl(i, k) - pm_tl
              pkz(i, k) = pkz(i, k) - pm
            END DO
          END DO
!Reversible interpolation on layer NH pressure perturbation
!                 to recover  lastge NH pressure perturbation
          DO k=1,npz-1
            DO i=ifirst,ilast
              g_rat_tl(i, k) = (delp_tl(i, j, k)*delp(i, j, k+1)-delp(i&
&               , j, k)*delp_tl(i, j, k+1))/delp(i, j, k+1)**2
              g_rat(i, k) = delp(i, j, k)/delp(i, j, k+1)
              bb_tl(i, k) = 2.*g_rat_tl(i, k)
              bb(i, k) = 2.*(1.+g_rat(i, k))
              dd_tl(i, k) = 3.*(pkz_tl(i, k)+g_rat_tl(i, k)*pkz(i, k+1)+&
&               g_rat(i, k)*pkz_tl(i, k+1))
              dd(i, k) = 3.*(pkz(i, k)+g_rat(i, k)*pkz(i, k+1))
            END DO
          END DO
          DO i=ifirst,ilast
            bet_tl(i) = bb_tl(i, 1)
            bet(i) = bb(i, 1)
            pkc_tl(i, j, 1) = 0.0
            pkc(i, j, 1) = 0.
            pkc_tl(i, j, 2) = (dd_tl(i, 1)*bet(i)-dd(i, 1)*bet_tl(i))/&
&             bet(i)**2
            pkc(i, j, 2) = dd(i, 1)/bet(i)
            bb_tl(i, npz) = 0.0
            bb(i, npz) = 2.
            dd_tl(i, npz) = 3.*pkz_tl(i, npz)
            dd(i, npz) = 3.*pkz(i, npz)
          END DO
          DO k=2,npz
            DO i=ifirst,ilast
              gam_tl(i, k) = (g_rat_tl(i, k-1)*bet(i)-g_rat(i, k-1)*&
&               bet_tl(i))/bet(i)**2
              gam(i, k) = g_rat(i, k-1)/bet(i)
              bet_tl(i) = bb_tl(i, k) - gam_tl(i, k)
              bet(i) = bb(i, k) - gam(i, k)
              pkc_tl(i, j, k+1) = ((dd_tl(i, k)-pkc_tl(i, j, k))*bet(i)-&
&               (dd(i, k)-pkc(i, j, k))*bet_tl(i))/bet(i)**2
              pkc(i, j, k+1) = (dd(i, k)-pkc(i, j, k))/bet(i)
            END DO
          END DO
          DO k=npz,2,-1
            DO i=ifirst,ilast
              pkc_tl(i, j, k) = pkc_tl(i, j, k) - gam_tl(i, k)*pkc(i, j&
&               , k+1) - gam(i, k)*pkc_tl(i, j, k+1)
              pkc(i, j, k) = pkc(i, j, k) - gam(i, k)*pkc(i, j, k+1)
            END DO
          END DO
        END DO
        DO j=npy,jlast
          IF (.NOT.pkc_pertn) THEN
            DO k=npz+1,1,-1
              DO i=ifirst,ilast
                pkc_tl(i, j, k) = pkc_tl(i, j, k) + pe_tl(i, k, j)
                pkc(i, j, k) = pkc(i, j, k) + pe(i, k, j)
              END DO
            END DO
          END IF
!pk3 if necessary
          IF (computepk3) THEN
            DO i=ifirst,ilast
              pk3_tl(i, j, 1) = 0.0
              pk3(i, j, 1) = ptk
            END DO
            DO k=2,npz+1
              DO i=ifirst,ilast
                arg1_tl = kappa*pe_tl(i, k, j)/pe(i, k, j)
                arg1 = kappa*LOG(pe(i, k, j))
                pk3_tl(i, j, k) = arg1_tl*EXP(arg1)
                pk3(i, j, k) = EXP(arg1)
              END DO
            END DO
          END IF
        END DO
      END IF
    END IF
  END SUBROUTINE NEST_HALO_NH_TLM
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
    REAL :: arg1
    REAL :: arg2
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
              arg1 = -(delp(i, j, k)*rgrav/delz(i, j, k)*rdgas*pt(i, j, &
&               k))
              arg2 = gama*LOG(arg1)
              pkz(i, k) = EXP(arg2)
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
                arg1 = kappa*LOG(pe(i, k, j))
                pk3(i, j, k) = EXP(arg1)
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
              arg1 = -(delp(i, j, k)*rgrav/delz(i, j, k)*rdgas*pt(i, j, &
&               k))
              arg2 = gama*LOG(arg1)
              pkz(i, k) = EXP(arg2)
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
                arg1 = kappa*LOG(pe(i, k, j))
                pk3(i, j, k) = EXP(arg1)
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
              arg1 = -(delp(i, j, k)*rgrav/delz(i, j, k)*rdgas*pt(i, j, &
&               k))
              arg2 = gama*LOG(arg1)
              pkz(i, k) = EXP(arg2)
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
                arg1 = kappa*LOG(pe(i, k, j))
                pk3(i, j, k) = EXP(arg1)
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
              arg1 = -(delp(i, j, k)*rgrav/delz(i, j, k)*rdgas*pt(i, j, &
&               k))
              arg2 = gama*LOG(arg1)
              pkz(i, k) = EXP(arg2)
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
                arg1 = kappa*LOG(pe(i, k, j))
                pk3(i, j, k) = EXP(arg1)
              END DO
            END DO
          END IF
        END DO
      END IF
    END IF
  END SUBROUTINE NEST_HALO_NH

end module nh_utils_tlm_mod
