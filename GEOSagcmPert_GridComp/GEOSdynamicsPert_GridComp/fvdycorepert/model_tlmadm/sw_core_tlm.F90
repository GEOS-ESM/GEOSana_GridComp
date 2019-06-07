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
 module sw_core_tlm_mod

 use fv_mp_mod,        only: ng
 use tp_core_tlm_mod,  only: fv_tp_2d, pert_ppm, copy_corners
 use tp_core_tlm_mod,  only: fv_tp_2d_tlm, copy_corners_tlm
 use fv_mp_mod,        only: XDir, YDir
 use fv_mp_mod,        only: fill_corners
 use fv_mp_tlm_mod,    only: fill_corners_tlm
 use fv_arrays_mod,    only: fv_grid_type, fv_grid_bounds_type, fv_flags_type, r_grid
 use a2b_edge_tlm_mod, only: a2b_ord4
 use a2b_edge_tlm_mod, only: a2b_ord4_tlm
 use fv_arrays_nlm_mod,only: fpp

#ifdef SW_DYNAMICS
 use test_cases_mod,   only: test_case
#endif

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
  real, parameter:: b1 =   1./30.
  real, parameter:: b2 = -13./60.
  real, parameter:: b3 = -13./60.
  real, parameter:: b4 =  0.45
  real, parameter:: b5 = -0.05


!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

      private
      public :: c_sw, d_sw, fill_4corners, del6_vt_flux, divergence_corner, divergence_corner_nest
      public :: d2a2c_vect
      public :: c_sw_tlm, d_sw_tlm, fill_4corners_tlm, del6_vt_flux_tlm, divergence_corner_tlm, divergence_corner_nest_tlm
      public :: d2a2c_vect_tlm

CONTAINS
!  Differentiation of c_sw in forward (tangent) mode:
!   variations   of useful results: w delp ua uc ptc ut delpc va
!                vc vt divg_d wc pt
!   with respect to varying inputs: u v w delp ua uc ptc ut delpc
!                va vc vt divg_d wc pt
  SUBROUTINE C_SW_TLM(delpc, delpc_tl, delp, delp_tl, ptc, ptc_tl, pt, &
&   pt_tl, u, u_tl, v, v_tl, w, w_tl, uc, uc_tl, vc, vc_tl, ua, ua_tl, &
&   va, va_tl, wc, wc_tl, ut, ut_tl, vt, vt_tl, divg_d, divg_d_tl, nord&
&   , dt2, hydrostatic, dord4, bd, gridstruct, flagstruct)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: u&
&   , vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   u_tl, vc_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: v&
&   , uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   v_tl, uc_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: delp&
&   , pt, ua, va, ut, vt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delp_tl, pt_tl, ua_tl, va_tl, ut_tl, vt_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: delpc&
&   , ptc, wc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: &
&   delpc_tl, ptc_tl, wc_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   divg_d
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   divg_d_tl
    INTEGER, INTENT(IN) :: nord
    REAL, INTENT(IN) :: dt2
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: dord4
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! Local:
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+1) :: vort, ke
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+1) :: vort_tl, ke_tl
    REAL, DIMENSION(bd%is-1:bd%ie+2, bd%js-1:bd%je+1) :: fx, fx1, fx2
    REAL, DIMENSION(bd%is-1:bd%ie+2, bd%js-1:bd%je+1) :: fx_tl, fx1_tl, &
&   fx2_tl
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+2) :: fy, fy1, fy2
    REAL, DIMENSION(bd%is-1:bd%ie+1, bd%js-1:bd%je+2) :: fy_tl, fy1_tl, &
&   fy2_tl
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
    CALL D2A2C_VECT_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, uc, &
&                 uc_tl, vc, vc_tl, ut, ut_tl, vt, vt_tl, dord4, &
&                 gridstruct, bd, npx, npy, nested, flagstruct%grid_type&
&                )
    IF (nord .GT. 0) THEN
      IF (nested) THEN
        CALL DIVERGENCE_CORNER_NEST_TLM(u, u_tl, v, v_tl, ua, ua_tl, va&
&                                 , va_tl, divg_d, divg_d_tl, gridstruct&
&                                 , flagstruct, bd)
      ELSE
        CALL DIVERGENCE_CORNER_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, &
&                            va_tl, divg_d, divg_d_tl, gridstruct, &
&                            flagstruct, bd)
      END IF
    END IF
    DO j=js-1,jep1
      DO i=is-1,iep1+1
        IF (ut(i, j) .GT. 0.) THEN
          ut_tl(i, j) = dt2*dy(i, j)*sin_sg(i-1, j, 3)*ut_tl(i, j)
          ut(i, j) = dt2*ut(i, j)*dy(i, j)*sin_sg(i-1, j, 3)
        ELSE
          ut_tl(i, j) = dt2*dy(i, j)*sin_sg(i, j, 1)*ut_tl(i, j)
          ut(i, j) = dt2*ut(i, j)*dy(i, j)*sin_sg(i, j, 1)
        END IF
      END DO
    END DO
    DO j=js-1,je+2
      DO i=is-1,iep1
        IF (vt(i, j) .GT. 0.) THEN
          vt_tl(i, j) = dt2*dx(i, j)*sin_sg(i, j-1, 4)*vt_tl(i, j)
          vt(i, j) = dt2*vt(i, j)*dx(i, j)*sin_sg(i, j-1, 4)
        ELSE
          vt_tl(i, j) = dt2*dx(i, j)*sin_sg(i, j, 2)*vt_tl(i, j)
          vt(i, j) = dt2*vt(i, j)*dx(i, j)*sin_sg(i, j, 2)
        END IF
      END DO
    END DO
!----------------
! Transport delp:
!----------------
! Xdir:
    IF (flagstruct%grid_type .LT. 3 .AND. (.NOT.nested)) CALL &
&     FILL2_4CORNERS_TLM(delp, delp_tl, pt, pt_tl, 1, bd, npx, npy, &
&                  sw_corner, se_corner, ne_corner, nw_corner)
    IF (hydrostatic) THEN
      fx_tl = 0.0
      fx1_tl = 0.0
      DO j=js-1,jep1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1_tl(i, j) = delp_tl(i-1, j)
            fx1(i, j) = delp(i-1, j)
            fx_tl(i, j) = pt_tl(i-1, j)
            fx(i, j) = pt(i-1, j)
          ELSE
            fx1_tl(i, j) = delp_tl(i, j)
            fx1(i, j) = delp(i, j)
            fx_tl(i, j) = pt_tl(i, j)
            fx(i, j) = pt(i, j)
          END IF
          fx1_tl(i, j) = ut_tl(i, j)*fx1(i, j) + ut(i, j)*fx1_tl(i, j)
          fx1(i, j) = ut(i, j)*fx1(i, j)
          fx_tl(i, j) = fx1_tl(i, j)*fx(i, j) + fx1(i, j)*fx_tl(i, j)
          fx(i, j) = fx1(i, j)*fx(i, j)
        END DO
      END DO
      fx2_tl = 0.0
    ELSE
      IF (flagstruct%grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_TLM(w, w_tl, 1, bd, npx, npy, sw_corner, &
&                        se_corner, ne_corner, nw_corner)
        fx_tl = 0.0
        fx1_tl = 0.0
        fx2_tl = 0.0
      ELSE
        fx_tl = 0.0
        fx1_tl = 0.0
        fx2_tl = 0.0
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1_tl(i, j) = delp_tl(i-1, j)
            fx1(i, j) = delp(i-1, j)
            fx_tl(i, j) = pt_tl(i-1, j)
            fx(i, j) = pt(i-1, j)
            fx2_tl(i, j) = w_tl(i-1, j)
            fx2(i, j) = w(i-1, j)
          ELSE
            fx1_tl(i, j) = delp_tl(i, j)
            fx1(i, j) = delp(i, j)
            fx_tl(i, j) = pt_tl(i, j)
            fx(i, j) = pt(i, j)
            fx2_tl(i, j) = w_tl(i, j)
            fx2(i, j) = w(i, j)
          END IF
          fx1_tl(i, j) = ut_tl(i, j)*fx1(i, j) + ut(i, j)*fx1_tl(i, j)
          fx1(i, j) = ut(i, j)*fx1(i, j)
          fx_tl(i, j) = fx1_tl(i, j)*fx(i, j) + fx1(i, j)*fx_tl(i, j)
          fx(i, j) = fx1(i, j)*fx(i, j)
          fx2_tl(i, j) = fx1_tl(i, j)*fx2(i, j) + fx1(i, j)*fx2_tl(i, j)
          fx2(i, j) = fx1(i, j)*fx2(i, j)
        END DO
      END DO
    END IF
! Ydir:
    IF (flagstruct%grid_type .LT. 3 .AND. (.NOT.nested)) CALL &
&     FILL2_4CORNERS_TLM(delp, delp_tl, pt, pt_tl, 2, bd, npx, npy, &
&                  sw_corner, se_corner, ne_corner, nw_corner)
    IF (hydrostatic) THEN
      fy1_tl = 0.0
      fy_tl = 0.0
      DO j=js-1,jep1+1
        DO i=is-1,iep1
          IF (vt(i, j) .GT. 0.) THEN
            fy1_tl(i, j) = delp_tl(i, j-1)
            fy1(i, j) = delp(i, j-1)
            fy_tl(i, j) = pt_tl(i, j-1)
            fy(i, j) = pt(i, j-1)
          ELSE
            fy1_tl(i, j) = delp_tl(i, j)
            fy1(i, j) = delp(i, j)
            fy_tl(i, j) = pt_tl(i, j)
            fy(i, j) = pt(i, j)
          END IF
          fy1_tl(i, j) = vt_tl(i, j)*fy1(i, j) + vt(i, j)*fy1_tl(i, j)
          fy1(i, j) = vt(i, j)*fy1(i, j)
          fy_tl(i, j) = fy1_tl(i, j)*fy(i, j) + fy1(i, j)*fy_tl(i, j)
          fy(i, j) = fy1(i, j)*fy(i, j)
        END DO
      END DO
      DO j=js-1,jep1
        DO i=is-1,iep1
          delpc_tl(i, j) = delp_tl(i, j) + gridstruct%rarea(i, j)*(&
&           fx1_tl(i, j)-fx1_tl(i+1, j)+fy1_tl(i, j)-fy1_tl(i, j+1))
          delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+(fy1(i, j)-&
&           fy1(i, j+1)))*gridstruct%rarea(i, j)
          ptc_tl(i, j) = ((pt_tl(i, j)*delp(i, j)+pt(i, j)*delp_tl(i, j)&
&           +gridstruct%rarea(i, j)*(fx_tl(i, j)-fx_tl(i+1, j)+fy_tl(i, &
&           j)-fy_tl(i, j+1)))*delpc(i, j)-(pt(i, j)*delp(i, j)+(fx(i, j&
&           )-fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*gridstruct%rarea(i, j))*&
&           delpc_tl(i, j))/delpc(i, j)**2
          ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+(fy(i, j&
&           )-fy(i, j+1)))*gridstruct%rarea(i, j))/delpc(i, j)
        END DO
      END DO
    ELSE
      IF (flagstruct%grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_TLM(w, w_tl, 2, bd, npx, npy, sw_corner, &
&                        se_corner, ne_corner, nw_corner)
        fy1_tl = 0.0
        fy2_tl = 0.0
        fy_tl = 0.0
      ELSE
        fy1_tl = 0.0
        fy2_tl = 0.0
        fy_tl = 0.0
      END IF
      DO j=js-1,je+2
        DO i=is-1,ie+1
          IF (vt(i, j) .GT. 0.) THEN
            fy1_tl(i, j) = delp_tl(i, j-1)
            fy1(i, j) = delp(i, j-1)
            fy_tl(i, j) = pt_tl(i, j-1)
            fy(i, j) = pt(i, j-1)
            fy2_tl(i, j) = w_tl(i, j-1)
            fy2(i, j) = w(i, j-1)
          ELSE
            fy1_tl(i, j) = delp_tl(i, j)
            fy1(i, j) = delp(i, j)
            fy_tl(i, j) = pt_tl(i, j)
            fy(i, j) = pt(i, j)
            fy2_tl(i, j) = w_tl(i, j)
            fy2(i, j) = w(i, j)
          END IF
          fy1_tl(i, j) = vt_tl(i, j)*fy1(i, j) + vt(i, j)*fy1_tl(i, j)
          fy1(i, j) = vt(i, j)*fy1(i, j)
          fy_tl(i, j) = fy1_tl(i, j)*fy(i, j) + fy1(i, j)*fy_tl(i, j)
          fy(i, j) = fy1(i, j)*fy(i, j)
          fy2_tl(i, j) = fy1_tl(i, j)*fy2(i, j) + fy1(i, j)*fy2_tl(i, j)
          fy2(i, j) = fy1(i, j)*fy2(i, j)
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          delpc_tl(i, j) = delp_tl(i, j) + gridstruct%rarea(i, j)*(&
&           fx1_tl(i, j)-fx1_tl(i+1, j)+fy1_tl(i, j)-fy1_tl(i, j+1))
          delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+(fy1(i, j)-&
&           fy1(i, j+1)))*gridstruct%rarea(i, j)
          ptc_tl(i, j) = ((pt_tl(i, j)*delp(i, j)+pt(i, j)*delp_tl(i, j)&
&           +gridstruct%rarea(i, j)*(fx_tl(i, j)-fx_tl(i+1, j)+fy_tl(i, &
&           j)-fy_tl(i, j+1)))*delpc(i, j)-(pt(i, j)*delp(i, j)+(fx(i, j&
&           )-fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*gridstruct%rarea(i, j))*&
&           delpc_tl(i, j))/delpc(i, j)**2
          ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+(fy(i, j&
&           )-fy(i, j+1)))*gridstruct%rarea(i, j))/delpc(i, j)
          wc_tl(i, j) = ((w_tl(i, j)*delp(i, j)+w(i, j)*delp_tl(i, j)+&
&           gridstruct%rarea(i, j)*(fx2_tl(i, j)-fx2_tl(i+1, j)+fy2_tl(i&
&           , j)-fy2_tl(i, j+1)))*delpc(i, j)-(w(i, j)*delp(i, j)+(fx2(i&
&           , j)-fx2(i+1, j)+(fy2(i, j)-fy2(i, j+1)))*gridstruct%rarea(i&
&           , j))*delpc_tl(i, j))/delpc(i, j)**2
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
      ke_tl = 0.0
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (ua(i, j) .GT. 0.) THEN
            ke_tl(i, j) = uc_tl(i, j)
            ke(i, j) = uc(i, j)
          ELSE
            ke_tl(i, j) = uc_tl(i+1, j)
            ke(i, j) = uc(i+1, j)
          END IF
        END DO
      END DO
      vort_tl = 0.0
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (va(i, j) .GT. 0.) THEN
            vort_tl(i, j) = vc_tl(i, j)
            vort(i, j) = vc(i, j)
          ELSE
            vort_tl(i, j) = vc_tl(i, j+1)
            vort(i, j) = vc(i, j+1)
          END IF
        END DO
      END DO
    ELSE
      ke_tl = 0.0
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (ua(i, j) .GT. 0.) THEN
            IF (i .EQ. 1) THEN
              ke_tl(1, j) = sin_sg(1, j, 1)*uc_tl(1, j) + cos_sg(1, j, 1&
&               )*v_tl(1, j)
              ke(1, j) = uc(1, j)*sin_sg(1, j, 1) + v(1, j)*cos_sg(1, j&
&               , 1)
            ELSE IF (i .EQ. npx) THEN
              ke_tl(i, j) = sin_sg(npx, j, 1)*uc_tl(npx, j) + cos_sg(npx&
&               , j, 1)*v_tl(npx, j)
              ke(i, j) = uc(npx, j)*sin_sg(npx, j, 1) + v(npx, j)*cos_sg&
&               (npx, j, 1)
            ELSE
              ke_tl(i, j) = uc_tl(i, j)
              ke(i, j) = uc(i, j)
            END IF
          ELSE IF (i .EQ. 0) THEN
            ke_tl(0, j) = sin_sg(0, j, 3)*uc_tl(1, j) + cos_sg(0, j, 3)*&
&             v_tl(1, j)
            ke(0, j) = uc(1, j)*sin_sg(0, j, 3) + v(1, j)*cos_sg(0, j, 3&
&             )
          ELSE IF (i .EQ. npx - 1) THEN
            ke_tl(i, j) = sin_sg(npx-1, j, 3)*uc_tl(npx, j) + cos_sg(npx&
&             -1, j, 3)*v_tl(npx, j)
            ke(i, j) = uc(npx, j)*sin_sg(npx-1, j, 3) + v(npx, j)*cos_sg&
&             (npx-1, j, 3)
          ELSE
            ke_tl(i, j) = uc_tl(i+1, j)
            ke(i, j) = uc(i+1, j)
          END IF
        END DO
      END DO
      vort_tl = 0.0
      DO j=js-1,jep1
        DO i=is-1,iep1
          IF (va(i, j) .GT. 0.) THEN
            IF (j .EQ. 1) THEN
              vort_tl(i, 1) = sin_sg(i, 1, 2)*vc_tl(i, 1) + cos_sg(i, 1&
&               , 2)*u_tl(i, 1)
              vort(i, 1) = vc(i, 1)*sin_sg(i, 1, 2) + u(i, 1)*cos_sg(i, &
&               1, 2)
            ELSE IF (j .EQ. npy) THEN
              vort_tl(i, j) = sin_sg(i, npy, 2)*vc_tl(i, npy) + cos_sg(i&
&               , npy, 2)*u_tl(i, npy)
              vort(i, j) = vc(i, npy)*sin_sg(i, npy, 2) + u(i, npy)*&
&               cos_sg(i, npy, 2)
            ELSE
              vort_tl(i, j) = vc_tl(i, j)
              vort(i, j) = vc(i, j)
            END IF
          ELSE IF (j .EQ. 0) THEN
            vort_tl(i, 0) = sin_sg(i, 0, 4)*vc_tl(i, 1) + cos_sg(i, 0, 4&
&             )*u_tl(i, 1)
            vort(i, 0) = vc(i, 1)*sin_sg(i, 0, 4) + u(i, 1)*cos_sg(i, 0&
&             , 4)
          ELSE IF (j .EQ. npy - 1) THEN
            vort_tl(i, j) = sin_sg(i, npy-1, 4)*vc_tl(i, npy) + cos_sg(i&
&             , npy-1, 4)*u_tl(i, npy)
            vort(i, j) = vc(i, npy)*sin_sg(i, npy-1, 4) + u(i, npy)*&
&             cos_sg(i, npy-1, 4)
          ELSE
            vort_tl(i, j) = vc_tl(i, j+1)
            vort(i, j) = vc(i, j+1)
          END IF
        END DO
      END DO
    END IF
    dt4 = 0.5*dt2
    DO j=js-1,jep1
      DO i=is-1,iep1
        ke_tl(i, j) = dt4*(ua_tl(i, j)*ke(i, j)+ua(i, j)*ke_tl(i, j)+&
&         va_tl(i, j)*vort(i, j)+va(i, j)*vort_tl(i, j))
        ke(i, j) = dt4*(ua(i, j)*ke(i, j)+va(i, j)*vort(i, j))
      END DO
    END DO
!------------------------------
! Compute circulation on C grid
!------------------------------
! To consider using true co-variant winds at face edges?
    DO j=js-1,je+1
      DO i=is,ie+1
        fx_tl(i, j) = dxc(i, j)*uc_tl(i, j)
        fx(i, j) = uc(i, j)*dxc(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is-1,ie+1
        fy_tl(i, j) = dyc(i, j)*vc_tl(i, j)
        fy(i, j) = vc(i, j)*dyc(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie+1
        vort_tl(i, j) = fx_tl(i, j-1) - fx_tl(i, j) + fy_tl(i, j) - &
&         fy_tl(i-1, j)
        vort(i, j) = fx(i, j-1) - fx(i, j) + (fy(i, j)-fy(i-1, j))
      END DO
    END DO
! Remove the extra term at the corners:
    IF (sw_corner) THEN
      vort_tl(1, 1) = vort_tl(1, 1) + fy_tl(0, 1)
      vort(1, 1) = vort(1, 1) + fy(0, 1)
    END IF
    IF (se_corner) THEN
      vort_tl(npx, 1) = vort_tl(npx, 1) - fy_tl(npx, 1)
      vort(npx, 1) = vort(npx, 1) - fy(npx, 1)
    END IF
    IF (ne_corner) THEN
      vort_tl(npx, npy) = vort_tl(npx, npy) - fy_tl(npx, npy)
      vort(npx, npy) = vort(npx, npy) - fy(npx, npy)
    END IF
    IF (nw_corner) THEN
      vort_tl(1, npy) = vort_tl(1, npy) + fy_tl(0, npy)
      vort(1, npy) = vort(1, npy) + fy(0, npy)
    END IF
!----------------------------
! Compute absolute vorticity
!----------------------------
    DO j=js,je+1
      DO i=is,ie+1
        vort_tl(i, j) = gridstruct%rarea_c(i, j)*vort_tl(i, j)
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
          fy1_tl(i, j) = dt2*(v_tl(i, j)-cosa_u(i, j)*uc_tl(i, j))/&
&           sina_u(i, j)
          fy1(i, j) = dt2*(v(i, j)-uc(i, j)*cosa_u(i, j))/sina_u(i, j)
          IF (fy1(i, j) .GT. 0.) THEN
            fy_tl(i, j) = vort_tl(i, j)
            fy(i, j) = vort(i, j)
          ELSE
            fy_tl(i, j) = vort_tl(i, j+1)
            fy(i, j) = vort(i, j+1)
          END IF
        END DO
      END DO
      DO j=js,jep1
        DO i=is,ie
          fx1_tl(i, j) = dt2*(u_tl(i, j)-cosa_v(i, j)*vc_tl(i, j))/&
&           sina_v(i, j)
          fx1(i, j) = dt2*(u(i, j)-vc(i, j)*cosa_v(i, j))/sina_v(i, j)
          IF (fx1(i, j) .GT. 0.) THEN
            fx_tl(i, j) = vort_tl(i, j)
            fx(i, j) = vort(i, j)
          ELSE
            fx_tl(i, j) = vort_tl(i+1, j)
            fx(i, j) = vort(i+1, j)
          END IF
        END DO
      END DO
    ELSE
      DO j=js,je
!DEC$ VECTOR ALWAYS
        DO i=is,iep1
          IF (i .EQ. 1 .OR. i .EQ. npx) THEN
            fy1_tl(i, j) = dt2*v_tl(i, j)
            fy1(i, j) = dt2*v(i, j)
          ELSE
            fy1_tl(i, j) = dt2*(v_tl(i, j)-cosa_u(i, j)*uc_tl(i, j))/&
&             sina_u(i, j)
            fy1(i, j) = dt2*(v(i, j)-uc(i, j)*cosa_u(i, j))/sina_u(i, j)
          END IF
          IF (fy1(i, j) .GT. 0.) THEN
            fy_tl(i, j) = vort_tl(i, j)
            fy(i, j) = vort(i, j)
          ELSE
            fy_tl(i, j) = vort_tl(i, j+1)
            fy(i, j) = vort(i, j+1)
          END IF
        END DO
      END DO
      DO j=js,jep1
        IF (j .EQ. 1 .OR. j .EQ. npy) THEN
!DEC$ VECTOR ALWAYS
          DO i=is,ie
            fx1_tl(i, j) = dt2*u_tl(i, j)
            fx1(i, j) = dt2*u(i, j)
            IF (fx1(i, j) .GT. 0.) THEN
              fx_tl(i, j) = vort_tl(i, j)
              fx(i, j) = vort(i, j)
            ELSE
              fx_tl(i, j) = vort_tl(i+1, j)
              fx(i, j) = vort(i+1, j)
            END IF
          END DO
        ELSE
!DEC$ VECTOR ALWAYS
          DO i=is,ie
            fx1_tl(i, j) = dt2*(u_tl(i, j)-cosa_v(i, j)*vc_tl(i, j))/&
&             sina_v(i, j)
            fx1(i, j) = dt2*(u(i, j)-vc(i, j)*cosa_v(i, j))/sina_v(i, j)
            IF (fx1(i, j) .GT. 0.) THEN
              fx_tl(i, j) = vort_tl(i, j)
              fx(i, j) = vort(i, j)
            ELSE
              fx_tl(i, j) = vort_tl(i+1, j)
              fx(i, j) = vort(i+1, j)
            END IF
          END DO
        END IF
      END DO
    END IF
! Update time-centered winds on the C-Grid
    DO j=js,je
      DO i=is,iep1
        uc_tl(i, j) = uc_tl(i, j) + fy1_tl(i, j)*fy(i, j) + fy1(i, j)*&
&         fy_tl(i, j) + gridstruct%rdxc(i, j)*(ke_tl(i-1, j)-ke_tl(i, j)&
&         )
        uc(i, j) = uc(i, j) + fy1(i, j)*fy(i, j) + gridstruct%rdxc(i, j)&
&         *(ke(i-1, j)-ke(i, j))
      END DO
    END DO
    DO j=js,jep1
      DO i=is,ie
        vc_tl(i, j) = vc_tl(i, j) - fx1_tl(i, j)*fx(i, j) - fx1(i, j)*&
&         fx_tl(i, j) + gridstruct%rdyc(i, j)*(ke_tl(i, j-1)-ke_tl(i, j)&
&         )
        vc(i, j) = vc(i, j) - fx1(i, j)*fx(i, j) + gridstruct%rdyc(i, j)&
&         *(ke(i, j-1)-ke(i, j))
      END DO
    END DO
  END SUBROUTINE C_SW_TLM
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
!  Differentiation of d_sw in forward (tangent) mode:
!   variations   of useful results: yfx_adv q crx_adv u v w delp
!                xfx_adv uc ptc xflux cry_adv delpc vc yflux divg_d
!                heat_source pt cx cy dpx
!   with respect to varying inputs: yfx_adv q crx_adv u v w delp
!                ua xfx_adv uc ptc xflux cry_adv delpc va vc yflux
!                divg_d z_rat heat_source pt cx cy dpx
!     d_sw :: D-Grid Shallow Water Routine
  SUBROUTINE D_SW_TLM(delpc, delpc_tl, delp, delp_tl, ptc, ptc_tl, pt, &
&   pt_tl, u, u_tl, v, v_tl, w, w_tl, uc, uc_tl, vc, vc_tl, ua, ua_tl, &
&   va, va_tl, divg_d, divg_d_tl, xflux, xflux_tl, yflux, yflux_tl, cx, &
&   cx_tl, cy, cy_tl, crx_adv, crx_adv_tl, cry_adv, cry_adv_tl, xfx_adv&
&   , xfx_adv_tl, yfx_adv, yfx_adv_tl, q_con, z_rat, z_rat_tl, kgb, &
&   heat_source, heat_source_tl, dpx, dpx_tl, zvir, sphum, nq, q, q_tl, &
&   k, km, inline_q, dt, hord_tr, hord_mt, hord_vt, hord_tm, hord_dp, &
&   nord, nord_v, nord_w, nord_t, dddmp, d2_bg, d4_bg, damp_v, damp_w, &
&   damp_t, d_con, hydrostatic, gridstruct, flagstruct, bd, hord_tr_pert&
&   , hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_dp_pert, split_damp&
&   , nord_pert, nord_v_pert, nord_w_pert, nord_t_pert, dddmp_pert, &
&   d2_bg_pert, d4_bg_pert, damp_v_pert, damp_w_pert, damp_t_pert)
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
    REAL, INTENT(INOUT) :: divg_d_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: z_rat
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: &
&   z_rat_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: delp&
&   , pt, ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delp_tl, pt_tl, ua_tl, va_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: w_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   q_con
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: u&
&   , vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   u_tl, vc_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: v&
&   , uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   v_tl, uc_tl
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, km, nq)
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed, km, nq)
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: delpc&
&   , ptc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: &
&   delpc_tl, ptc_tl
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je), INTENT(OUT) :: &
&   heat_source
    REAL, DIMENSION(bd%is:bd%ie, bd%js:bd%je), INTENT(OUT) :: &
&   heat_source_tl
    REAL(kind=8), DIMENSION(bd%is:bd%ie, bd%js:bd%je), INTENT(INOUT) :: &
&   dpx
    REAL(kind=8), DIMENSION(bd%is:bd%ie, bd%js:bd%je), INTENT(INOUT) :: &
&   dpx_tl
! The flux capacitors:
    REAL, INTENT(INOUT) :: xflux(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: xflux_tl(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, INTENT(INOUT) :: yflux(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: yflux_tl(bd%is:bd%ie, bd%js:bd%je+1)
!------------------------
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: cx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: cy_tl(bd%isd:bd%ied, bd%js:bd%je+1)
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: inline_q
    REAL, DIMENSION(bd%is:bd%ie+1, bd%jsd:bd%jed), INTENT(OUT) :: &
&   crx_adv, xfx_adv
    REAL, DIMENSION(bd%is:bd%ie+1, bd%jsd:bd%jed), INTENT(OUT) :: &
&   crx_adv_tl, xfx_adv_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%js:bd%je+1), INTENT(OUT) :: &
&   cry_adv, yfx_adv
    REAL, DIMENSION(bd%isd:bd%ied, bd%js:bd%je+1), INTENT(OUT) :: &
&   cry_adv_tl, yfx_adv_tl
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! Local:
    LOGICAL :: sw_corner, se_corner, ne_corner, nw_corner
    REAL :: ut(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: ut_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: vt_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!---
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: fx2_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: fy2(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: fy2_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!  work array
    REAL :: dw(bd%is:bd%ie, bd%js:bd%je)
    REAL :: dw_tl(bd%is:bd%ie, bd%js:bd%je)
!---
    REAL, DIMENSION(bd%is:bd%ie+1, bd%js:bd%je+1) :: ub, vb
    REAL, DIMENSION(bd%is:bd%ie+1, bd%js:bd%je+1) :: ub_tl, vb_tl
!  work array
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
!  needs this for corner_comm
    REAL :: ke(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
    REAL :: ke_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed+1)
! Vorticity
    REAL :: vort(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: vort_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
! 1-D X-direction Fluxes
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_tl(bd%is:bd%ie+1, bd%js:bd%je)
! 1-D Y-direction Fluxes
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_tl(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_tl(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_tl(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: gx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: gx_tl(bd%is:bd%ie+1, bd%js:bd%je)
! work Y-dir flux array
    REAL :: gy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: gy_tl(bd%is:bd%ie, bd%js:bd%je+1)
    LOGICAL :: fill_c
    REAL :: dt2, dt4, dt5, dt6
    REAL :: damp, damp2, damp4, dd8, u2, v2, du2, dv2
    REAL :: u2_tl, v2_tl, du2_tl, dv2_tl
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
    REAL*8 :: pwx1
    INTEGER :: pwy1
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
        ut_tl = 0.0
        DO j=jsd,jed
          DO i=is-1,ie+2
            ut_tl(i, j) = rsin_u(i, j)*(uc_tl(i, j)-0.25*cosa_u(i, j)*(&
&             vc_tl(i-1, j)+vc_tl(i, j)+vc_tl(i-1, j+1)+vc_tl(i, j+1)))
            ut(i, j) = (uc(i, j)-0.25*cosa_u(i, j)*(vc(i-1, j)+vc(i, j)+&
&             vc(i-1, j+1)+vc(i, j+1)))*rsin_u(i, j)
          END DO
        END DO
        vt_tl = 0.0
        DO j=js-1,je+2
          DO i=isd,ied
            vt_tl(i, j) = rsin_v(i, j)*(vc_tl(i, j)-0.25*cosa_v(i, j)*(&
&             uc_tl(i, j-1)+uc_tl(i+1, j-1)+uc_tl(i, j)+uc_tl(i+1, j)))
            vt(i, j) = (vc(i, j)-0.25*cosa_v(i, j)*(uc(i, j-1)+uc(i+1, j&
&             -1)+uc(i, j)+uc(i+1, j)))*rsin_v(i, j)
          END DO
        END DO
      ELSE
        ut_tl = 0.0
        DO j=jsd,jed
          IF (j .NE. 0 .AND. j .NE. 1 .AND. j .NE. npy - 1 .AND. j .NE. &
&             npy) THEN
            DO i=is-1,ie+2
              ut_tl(i, j) = rsin_u(i, j)*(uc_tl(i, j)-0.25*cosa_u(i, j)*&
&               (vc_tl(i-1, j)+vc_tl(i, j)+vc_tl(i-1, j+1)+vc_tl(i, j+1)&
&               ))
              ut(i, j) = (uc(i, j)-0.25*cosa_u(i, j)*(vc(i-1, j)+vc(i, j&
&               )+vc(i-1, j+1)+vc(i, j+1)))*rsin_u(i, j)
            END DO
          END IF
        END DO
        vt_tl = 0.0
        DO j=js-1,je+2
          IF (j .NE. 1 .AND. j .NE. npy) THEN
            DO i=isd,ied
              vt_tl(i, j) = rsin_v(i, j)*(vc_tl(i, j)-0.25*cosa_v(i, j)*&
&               (uc_tl(i, j-1)+uc_tl(i+1, j-1)+uc_tl(i, j)+uc_tl(i+1, j)&
&               ))
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
              ut_tl(1, j) = uc_tl(1, j)/sin_sg(0, j, 3)
              ut(1, j) = uc(1, j)/sin_sg(0, j, 3)
            ELSE
              ut_tl(1, j) = uc_tl(1, j)/sin_sg(1, j, 1)
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
            vt_tl(0, j) = vc_tl(0, j) - 0.25*cosa_v(0, j)*(ut_tl(0, j-1)&
&             +ut_tl(1, j-1)+ut_tl(0, j)+ut_tl(1, j))
            vt(0, j) = vc(0, j) - 0.25*cosa_v(0, j)*(ut(0, j-1)+ut(1, j-&
&             1)+ut(0, j)+ut(1, j))
            vt_tl(1, j) = vc_tl(1, j) - 0.25*cosa_v(1, j)*(ut_tl(1, j-1)&
&             +ut_tl(2, j-1)+ut_tl(1, j)+ut_tl(2, j))
            vt(1, j) = vc(1, j) - 0.25*cosa_v(1, j)*(ut(1, j-1)+ut(2, j-&
&             1)+ut(1, j)+ut(2, j))
          END DO
        END IF
! East edge:
        IF (ie + 1 .EQ. npx) THEN
          DO j=jsd,jed
            IF (uc(npx, j)*dt .GT. 0.) THEN
              ut_tl(npx, j) = uc_tl(npx, j)/sin_sg(npx-1, j, 3)
              ut(npx, j) = uc(npx, j)/sin_sg(npx-1, j, 3)
            ELSE
              ut_tl(npx, j) = uc_tl(npx, j)/sin_sg(npx, j, 1)
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
            vt_tl(npx-1, j) = vc_tl(npx-1, j) - 0.25*cosa_v(npx-1, j)*(&
&             ut_tl(npx-1, j-1)+ut_tl(npx, j-1)+ut_tl(npx-1, j)+ut_tl(&
&             npx, j))
            vt(npx-1, j) = vc(npx-1, j) - 0.25*cosa_v(npx-1, j)*(ut(npx-&
&             1, j-1)+ut(npx, j-1)+ut(npx-1, j)+ut(npx, j))
            vt_tl(npx, j) = vc_tl(npx, j) - 0.25*cosa_v(npx, j)*(ut_tl(&
&             npx, j-1)+ut_tl(npx+1, j-1)+ut_tl(npx, j)+ut_tl(npx+1, j))
            vt(npx, j) = vc(npx, j) - 0.25*cosa_v(npx, j)*(ut(npx, j-1)+&
&             ut(npx+1, j-1)+ut(npx, j)+ut(npx+1, j))
          END DO
        END IF
! South (Bottom) edge:
        IF (js .EQ. 1) THEN
          DO i=isd,ied
            IF (vc(i, 1)*dt .GT. 0.) THEN
              vt_tl(i, 1) = vc_tl(i, 1)/sin_sg(i, 0, 4)
              vt(i, 1) = vc(i, 1)/sin_sg(i, 0, 4)
            ELSE
              vt_tl(i, 1) = vc_tl(i, 1)/sin_sg(i, 1, 2)
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
            ut_tl(i, 0) = uc_tl(i, 0) - 0.25*cosa_u(i, 0)*(vt_tl(i-1, 0)&
&             +vt_tl(i, 0)+vt_tl(i-1, 1)+vt_tl(i, 1))
            ut(i, 0) = uc(i, 0) - 0.25*cosa_u(i, 0)*(vt(i-1, 0)+vt(i, 0)&
&             +vt(i-1, 1)+vt(i, 1))
            ut_tl(i, 1) = uc_tl(i, 1) - 0.25*cosa_u(i, 1)*(vt_tl(i-1, 1)&
&             +vt_tl(i, 1)+vt_tl(i-1, 2)+vt_tl(i, 2))
            ut(i, 1) = uc(i, 1) - 0.25*cosa_u(i, 1)*(vt(i-1, 1)+vt(i, 1)&
&             +vt(i-1, 2)+vt(i, 2))
          END DO
        END IF
! North edge:
        IF (je + 1 .EQ. npy) THEN
          DO i=isd,ied
            IF (vc(i, npy)*dt .GT. 0.) THEN
              vt_tl(i, npy) = vc_tl(i, npy)/sin_sg(i, npy-1, 4)
              vt(i, npy) = vc(i, npy)/sin_sg(i, npy-1, 4)
            ELSE
              vt_tl(i, npy) = vc_tl(i, npy)/sin_sg(i, npy, 2)
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
            ut_tl(i, npy-1) = uc_tl(i, npy-1) - 0.25*cosa_u(i, npy-1)*(&
&             vt_tl(i-1, npy-1)+vt_tl(i, npy-1)+vt_tl(i-1, npy)+vt_tl(i&
&             , npy))
            ut(i, npy-1) = uc(i, npy-1) - 0.25*cosa_u(i, npy-1)*(vt(i-1&
&             , npy-1)+vt(i, npy-1)+vt(i-1, npy)+vt(i, npy))
            ut_tl(i, npy) = uc_tl(i, npy) - 0.25*cosa_u(i, npy)*(vt_tl(i&
&             -1, npy)+vt_tl(i, npy)+vt_tl(i-1, npy+1)+vt_tl(i, npy+1))
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
          ut_tl(2, 0) = damp*(uc_tl(2, 0)-0.25*cosa_u(2, 0)*(vt_tl(1, 1)&
&           +vt_tl(2, 1)+vt_tl(2, 0)+vc_tl(1, 0)-0.25*cosa_v(1, 0)*(&
&           ut_tl(1, 0)+ut_tl(1, -1)+ut_tl(2, -1))))
          ut(2, 0) = (uc(2, 0)-0.25*cosa_u(2, 0)*(vt(1, 1)+vt(2, 1)+vt(2&
&           , 0)+vc(1, 0)-0.25*cosa_v(1, 0)*(ut(1, 0)+ut(1, -1)+ut(2, -1&
&           ))))*damp
          damp = 1./(1.-0.0625*cosa_u(0, 1)*cosa_v(0, 2))
          vt_tl(0, 2) = damp*(vc_tl(0, 2)-0.25*cosa_v(0, 2)*(ut_tl(1, 1)&
&           +ut_tl(1, 2)+ut_tl(0, 2)+uc_tl(0, 1)-0.25*cosa_u(0, 1)*(&
&           vt_tl(0, 1)+vt_tl(-1, 1)+vt_tl(-1, 2))))
          vt(0, 2) = (vc(0, 2)-0.25*cosa_v(0, 2)*(ut(1, 1)+ut(1, 2)+ut(0&
&           , 2)+uc(0, 1)-0.25*cosa_u(0, 1)*(vt(0, 1)+vt(-1, 1)+vt(-1, 2&
&           ))))*damp
          damp = 1./(1.-0.0625*cosa_u(2, 1)*cosa_v(1, 2))
          ut_tl(2, 1) = damp*(uc_tl(2, 1)-0.25*cosa_u(2, 1)*(vt_tl(1, 1)&
&           +vt_tl(2, 1)+vt_tl(2, 2)+vc_tl(1, 2)-0.25*cosa_v(1, 2)*(&
&           ut_tl(1, 1)+ut_tl(1, 2)+ut_tl(2, 2))))
          ut(2, 1) = (uc(2, 1)-0.25*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(2&
&           , 2)+vc(1, 2)-0.25*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(2, 2))&
&           ))*damp
          vt_tl(1, 2) = damp*(vc_tl(1, 2)-0.25*cosa_v(1, 2)*(ut_tl(1, 1)&
&           +ut_tl(1, 2)+ut_tl(2, 2)+uc_tl(2, 1)-0.25*cosa_u(2, 1)*(&
&           vt_tl(1, 1)+vt_tl(2, 1)+vt_tl(2, 2))))
          vt(1, 2) = (vc(1, 2)-0.25*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(2&
&           , 2)+uc(2, 1)-0.25*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(2, 2))&
&           ))*damp
        END IF
        IF (se_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(npx-1, 0)*cosa_v(npx-1, 0))
          ut_tl(npx-1, 0) = damp*(uc_tl(npx-1, 0)-0.25*cosa_u(npx-1, 0)*&
&           (vt_tl(npx-1, 1)+vt_tl(npx-2, 1)+vt_tl(npx-2, 0)+vc_tl(npx-1&
&           , 0)-0.25*cosa_v(npx-1, 0)*(ut_tl(npx, 0)+ut_tl(npx, -1)+&
&           ut_tl(npx-1, -1))))
          ut(npx-1, 0) = (uc(npx-1, 0)-0.25*cosa_u(npx-1, 0)*(vt(npx-1, &
&           1)+vt(npx-2, 1)+vt(npx-2, 0)+vc(npx-1, 0)-0.25*cosa_v(npx-1&
&           , 0)*(ut(npx, 0)+ut(npx, -1)+ut(npx-1, -1))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx+1, 1)*cosa_v(npx, 2))
          vt_tl(npx, 2) = damp*(vc_tl(npx, 2)-0.25*cosa_v(npx, 2)*(ut_tl&
&           (npx, 1)+ut_tl(npx, 2)+ut_tl(npx+1, 2)+uc_tl(npx+1, 1)-0.25*&
&           cosa_u(npx+1, 1)*(vt_tl(npx, 1)+vt_tl(npx+1, 1)+vt_tl(npx+1&
&           , 2))))
          vt(npx, 2) = (vc(npx, 2)-0.25*cosa_v(npx, 2)*(ut(npx, 1)+ut(&
&           npx, 2)+ut(npx+1, 2)+uc(npx+1, 1)-0.25*cosa_u(npx+1, 1)*(vt(&
&           npx, 1)+vt(npx+1, 1)+vt(npx+1, 2))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx-1, 1)*cosa_v(npx-1, 2))
          ut_tl(npx-1, 1) = damp*(uc_tl(npx-1, 1)-0.25*cosa_u(npx-1, 1)*&
&           (vt_tl(npx-1, 1)+vt_tl(npx-2, 1)+vt_tl(npx-2, 2)+vc_tl(npx-1&
&           , 2)-0.25*cosa_v(npx-1, 2)*(ut_tl(npx, 1)+ut_tl(npx, 2)+&
&           ut_tl(npx-1, 2))))
          ut(npx-1, 1) = (uc(npx-1, 1)-0.25*cosa_u(npx-1, 1)*(vt(npx-1, &
&           1)+vt(npx-2, 1)+vt(npx-2, 2)+vc(npx-1, 2)-0.25*cosa_v(npx-1&
&           , 2)*(ut(npx, 1)+ut(npx, 2)+ut(npx-1, 2))))*damp
          vt_tl(npx-1, 2) = damp*(vc_tl(npx-1, 2)-0.25*cosa_v(npx-1, 2)*&
&           (ut_tl(npx, 1)+ut_tl(npx, 2)+ut_tl(npx-1, 2)+uc_tl(npx-1, 1)&
&           -0.25*cosa_u(npx-1, 1)*(vt_tl(npx-1, 1)+vt_tl(npx-2, 1)+&
&           vt_tl(npx-2, 2))))
          vt(npx-1, 2) = (vc(npx-1, 2)-0.25*cosa_v(npx-1, 2)*(ut(npx, 1)&
&           +ut(npx, 2)+ut(npx-1, 2)+uc(npx-1, 1)-0.25*cosa_u(npx-1, 1)*&
&           (vt(npx-1, 1)+vt(npx-2, 1)+vt(npx-2, 2))))*damp
        END IF
        IF (ne_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(npx-1, npy)*cosa_v(npx-1, npy+1))
          ut_tl(npx-1, npy) = damp*(uc_tl(npx-1, npy)-0.25*cosa_u(npx-1&
&           , npy)*(vt_tl(npx-1, npy)+vt_tl(npx-2, npy)+vt_tl(npx-2, npy&
&           +1)+vc_tl(npx-1, npy+1)-0.25*cosa_v(npx-1, npy+1)*(ut_tl(npx&
&           , npy)+ut_tl(npx, npy+1)+ut_tl(npx-1, npy+1))))
          ut(npx-1, npy) = (uc(npx-1, npy)-0.25*cosa_u(npx-1, npy)*(vt(&
&           npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy+1)+vc(npx-1, npy+1)&
&           -0.25*cosa_v(npx-1, npy+1)*(ut(npx, npy)+ut(npx, npy+1)+ut(&
&           npx-1, npy+1))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx+1, npy-1)*cosa_v(npx, npy-1))
          vt_tl(npx, npy-1) = damp*(vc_tl(npx, npy-1)-0.25*cosa_v(npx, &
&           npy-1)*(ut_tl(npx, npy-1)+ut_tl(npx, npy-2)+ut_tl(npx+1, npy&
&           -2)+uc_tl(npx+1, npy-1)-0.25*cosa_u(npx+1, npy-1)*(vt_tl(npx&
&           , npy)+vt_tl(npx+1, npy)+vt_tl(npx+1, npy-1))))
          vt(npx, npy-1) = (vc(npx, npy-1)-0.25*cosa_v(npx, npy-1)*(ut(&
&           npx, npy-1)+ut(npx, npy-2)+ut(npx+1, npy-2)+uc(npx+1, npy-1)&
&           -0.25*cosa_u(npx+1, npy-1)*(vt(npx, npy)+vt(npx+1, npy)+vt(&
&           npx+1, npy-1))))*damp
          damp = 1./(1.-0.0625*cosa_u(npx-1, npy-1)*cosa_v(npx-1, npy-1)&
&           )
          ut_tl(npx-1, npy-1) = damp*(uc_tl(npx-1, npy-1)-0.25*cosa_u(&
&           npx-1, npy-1)*(vt_tl(npx-1, npy)+vt_tl(npx-2, npy)+vt_tl(npx&
&           -2, npy-1)+vc_tl(npx-1, npy-1)-0.25*cosa_v(npx-1, npy-1)*(&
&           ut_tl(npx, npy-1)+ut_tl(npx, npy-2)+ut_tl(npx-1, npy-2))))
          ut(npx-1, npy-1) = (uc(npx-1, npy-1)-0.25*cosa_u(npx-1, npy-1)&
&           *(vt(npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy-1)+vc(npx-1, &
&           npy-1)-0.25*cosa_v(npx-1, npy-1)*(ut(npx, npy-1)+ut(npx, npy&
&           -2)+ut(npx-1, npy-2))))*damp
          vt_tl(npx-1, npy-1) = damp*(vc_tl(npx-1, npy-1)-0.25*cosa_v(&
&           npx-1, npy-1)*(ut_tl(npx, npy-1)+ut_tl(npx, npy-2)+ut_tl(npx&
&           -1, npy-2)+uc_tl(npx-1, npy-1)-0.25*cosa_u(npx-1, npy-1)*(&
&           vt_tl(npx-1, npy)+vt_tl(npx-2, npy)+vt_tl(npx-2, npy-1))))
          vt(npx-1, npy-1) = (vc(npx-1, npy-1)-0.25*cosa_v(npx-1, npy-1)&
&           *(ut(npx, npy-1)+ut(npx, npy-2)+ut(npx-1, npy-2)+uc(npx-1, &
&           npy-1)-0.25*cosa_u(npx-1, npy-1)*(vt(npx-1, npy)+vt(npx-2, &
&           npy)+vt(npx-2, npy-1))))*damp
        END IF
        IF (nw_corner) THEN
          damp = 1./(1.-0.0625*cosa_u(2, npy)*cosa_v(1, npy+1))
          ut_tl(2, npy) = damp*(uc_tl(2, npy)-0.25*cosa_u(2, npy)*(vt_tl&
&           (1, npy)+vt_tl(2, npy)+vt_tl(2, npy+1)+vc_tl(1, npy+1)-0.25*&
&           cosa_v(1, npy+1)*(ut_tl(1, npy)+ut_tl(1, npy+1)+ut_tl(2, npy&
&           +1))))
          ut(2, npy) = (uc(2, npy)-0.25*cosa_u(2, npy)*(vt(1, npy)+vt(2&
&           , npy)+vt(2, npy+1)+vc(1, npy+1)-0.25*cosa_v(1, npy+1)*(ut(1&
&           , npy)+ut(1, npy+1)+ut(2, npy+1))))*damp
          damp = 1./(1.-0.0625*cosa_u(0, npy-1)*cosa_v(0, npy-1))
          vt_tl(0, npy-1) = damp*(vc_tl(0, npy-1)-0.25*cosa_v(0, npy-1)*&
&           (ut_tl(1, npy-1)+ut_tl(1, npy-2)+ut_tl(0, npy-2)+uc_tl(0, &
&           npy-1)-0.25*cosa_u(0, npy-1)*(vt_tl(0, npy)+vt_tl(-1, npy)+&
&           vt_tl(-1, npy-1))))
          vt(0, npy-1) = (vc(0, npy-1)-0.25*cosa_v(0, npy-1)*(ut(1, npy-&
&           1)+ut(1, npy-2)+ut(0, npy-2)+uc(0, npy-1)-0.25*cosa_u(0, npy&
&           -1)*(vt(0, npy)+vt(-1, npy)+vt(-1, npy-1))))*damp
          damp = 1./(1.-0.0625*cosa_u(2, npy-1)*cosa_v(1, npy-1))
          ut_tl(2, npy-1) = damp*(uc_tl(2, npy-1)-0.25*cosa_u(2, npy-1)*&
&           (vt_tl(1, npy)+vt_tl(2, npy)+vt_tl(2, npy-1)+vc_tl(1, npy-1)&
&           -0.25*cosa_v(1, npy-1)*(ut_tl(1, npy-1)+ut_tl(1, npy-2)+&
&           ut_tl(2, npy-2))))
          ut(2, npy-1) = (uc(2, npy-1)-0.25*cosa_u(2, npy-1)*(vt(1, npy)&
&           +vt(2, npy)+vt(2, npy-1)+vc(1, npy-1)-0.25*cosa_v(1, npy-1)*&
&           (ut(1, npy-1)+ut(1, npy-2)+ut(2, npy-2))))*damp
          vt_tl(1, npy-1) = damp*(vc_tl(1, npy-1)-0.25*cosa_v(1, npy-1)*&
&           (ut_tl(1, npy-1)+ut_tl(1, npy-2)+ut_tl(2, npy-2)+uc_tl(2, &
&           npy-1)-0.25*cosa_u(2, npy-1)*(vt_tl(1, npy)+vt_tl(2, npy)+&
&           vt_tl(2, npy-1))))
          vt(1, npy-1) = (vc(1, npy-1)-0.25*cosa_v(1, npy-1)*(ut(1, npy-&
&           1)+ut(1, npy-2)+ut(2, npy-2)+uc(2, npy-1)-0.25*cosa_u(2, npy&
&           -1)*(vt(1, npy)+vt(2, npy)+vt(2, npy-1))))*damp
        END IF
      END IF
    ELSE
      ut_tl = 0.0
! flagstruct%grid_type >= 3
      DO j=jsd,jed
        DO i=is,ie+1
          ut_tl(i, j) = uc_tl(i, j)
          ut(i, j) = uc(i, j)
        END DO
      END DO
      vt_tl = 0.0
      DO j=js,je+1
        DO i=isd,ied
          vt_tl(i, j) = vc_tl(i, j)
          vt(i, j) = vc(i, j)
        END DO
      END DO
    END IF
    DO j=jsd,jed
      DO i=is,ie+1
        xfx_adv_tl(i, j) = dt*ut_tl(i, j)
        xfx_adv(i, j) = dt*ut(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        yfx_adv_tl(i, j) = dt*vt_tl(i, j)
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
          crx_adv_tl(i, j) = rdxa(i-1, j)*xfx_adv_tl(i, j)
          crx_adv(i, j) = xfx_adv(i, j)*rdxa(i-1, j)
          xfx_adv_tl(i, j) = dy(i, j)*sin_sg(i-1, j, 3)*xfx_adv_tl(i, j)
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sin_sg(i-1, j, 3)
        ELSE
          crx_adv_tl(i, j) = rdxa(i, j)*xfx_adv_tl(i, j)
          crx_adv(i, j) = xfx_adv(i, j)*rdxa(i, j)
          xfx_adv_tl(i, j) = dy(i, j)*sin_sg(i, j, 1)*xfx_adv_tl(i, j)
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sin_sg(i, j, 1)
        END IF
      END DO
    END DO
    DO j=js,je+1
!DEC$ VECTOR ALWAYS
      DO i=isd,ied
        IF (yfx_adv(i, j) .GT. 0.) THEN
          cry_adv_tl(i, j) = rdya(i, j-1)*yfx_adv_tl(i, j)
          cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j-1)
          yfx_adv_tl(i, j) = dx(i, j)*sin_sg(i, j-1, 4)*yfx_adv_tl(i, j)
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sin_sg(i, j-1, 4)
        ELSE
          cry_adv_tl(i, j) = rdya(i, j)*yfx_adv_tl(i, j)
          cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j)
          yfx_adv_tl(i, j) = dx(i, j)*sin_sg(i, j, 2)*yfx_adv_tl(i, j)
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sin_sg(i, j, 2)
        END IF
      END DO
    END DO
    ra_x_tl = 0.0
    DO j=jsd,jed
      DO i=is,ie
        ra_x_tl(i, j) = xfx_adv_tl(i, j) - xfx_adv_tl(i+1, j)
        ra_x(i, j) = area(i, j) + (xfx_adv(i, j)-xfx_adv(i+1, j))
      END DO
    END DO
    ra_y_tl = 0.0
    DO j=js,je
      DO i=isd,ied
        ra_y_tl(i, j) = yfx_adv_tl(i, j) - yfx_adv_tl(i, j+1)
        ra_y(i, j) = area(i, j) + (yfx_adv(i, j)-yfx_adv(i, j+1))
      END DO
    END DO
    IF (hord_dp .EQ. hord_dp_pert .AND. (.NOT.split_damp)) THEN
      fy_tl = 0.0
      fx_tl = 0.0
      CALL FV_TP_2D_TLM(delp, delp_tl, crx_adv, crx_adv_tl, cry_adv, &
&                    cry_adv_tl, npx, npy, hord_dp, fx, fx_tl, fy, fy_tl&
&                    , xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, &
&                    gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl, nord=&
&                    nord_v, damp_c=damp_v)
    ELSE
      fy_tl = 0.0
      fx_tl = 0.0
      CALL FV_TP_2D_TLM(delp, delp_tl, crx_adv, crx_adv_tl, cry_adv, &
&                 cry_adv_tl, npx, npy, hord_dp_pert, fx, fx_tl, fy, &
&                 fy_tl, xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, &
&                 gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl, nord=&
&                 nord_v_pert, damp_c=damp_v_pert)
      call fv_tp_2d(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx, fy,  &
                    xfx_adv,yfx_adv, gridstruct, bd, ra_x, ra_y, nord=nord_v, damp_c=damp_v)
    END IF
! <<< Save the mass fluxes to the "Flux Capacitor" for tracer transport >>>
    DO j=jsd,jed
      DO i=is,ie+1
        cx_tl(i, j) = cx_tl(i, j) + crx_adv_tl(i, j)
        cx(i, j) = cx(i, j) + crx_adv(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie+1
        xflux_tl(i, j) = xflux_tl(i, j) + fx_tl(i, j)
        xflux(i, j) = xflux(i, j) + fx(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        cy_tl(i, j) = cy_tl(i, j) + cry_adv_tl(i, j)
        cy(i, j) = cy(i, j) + cry_adv(i, j)
      END DO
      DO i=is,ie
        yflux_tl(i, j) = yflux_tl(i, j) + fy_tl(i, j)
        yflux(i, j) = yflux(i, j) + fy(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie
        heat_source_tl(i, j) = 0.0
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
        pwx1 = damp_w*gridstruct%da_min_c
        pwy1 = nord_w + 1
        damp4 = pwx1**pwy1
        fy2_tl = 0.0
        fx2_tl = 0.0
        wk_tl = 0.0
        CALL DEL6_VT_FLUX_TLM(nord_w, npx, npy, damp4, w, w_tl, wk, &
&                       wk_tl, fx2, fx2_tl, fy2, fy2_tl, gridstruct, bd)
        dw_tl = 0.0
        DO j=js,je
          DO i=is,ie
            dw_tl(i, j) = rarea(i, j)*(fx2_tl(i, j)-fx2_tl(i+1, j)+&
&             fy2_tl(i, j)-fy2_tl(i, j+1))
            dw(i, j) = (fx2(i, j)-fx2(i+1, j)+(fy2(i, j)-fy2(i, j+1)))*&
&             rarea(i, j)
! 0.5 * [ (w+dw)**2 - w**2 ] = w*dw + 0.5*dw*dw
!                   heat_source(i,j) = -d_con*dw(i,j)*(w(i,j)+0.5*dw(i,j))
            heat_source_tl(i, j) = -(dw_tl(i, j)*(w(i, j)+0.5*dw(i, j)))&
&             - dw(i, j)*(w_tl(i, j)+0.5*dw_tl(i, j))
            heat_source(i, j) = dd8 - dw(i, j)*(w(i, j)+0.5*dw(i, j))
          END DO
        END DO
      ELSE
        dw_tl = 0.0
        wk_tl = 0.0
      END IF
      IF (hord_vt .EQ. hord_vt_pert) THEN
        gy_tl = 0.0
        gx_tl = 0.0
        CALL FV_TP_2D_TLM(w, w_tl, crx_adv, crx_adv_tl, cry_adv, &
&                      cry_adv_tl, npx, npy, hord_vt, gx, gx_tl, gy, &
&                      gy_tl, xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, &
&                      gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl, mfx&
&                      =fx, mfx_tl=fx_tl, mfy=fy, mfy_tl=fy_tl)
      ELSE
        gy_tl = 0.0
        gx_tl = 0.0
        CALL FV_TP_2D_TLM(w, w_tl, crx_adv, crx_adv_tl, cry_adv, &
&                   cry_adv_tl, npx, npy, hord_vt_pert, gx, gx_tl, gy, &
&                   gy_tl, xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, &
&                   gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl, mfx=fx&
&                   , mfx_tl=fx_tl, mfy=fy, mfy_tl=fy_tl)
      call fv_tp_2d(w, crx_adv,cry_adv, npx, npy, hord_vt, gx, gy, xfx_adv, yfx_adv, &
                    gridstruct, bd, ra_x, ra_y, mfx=fx, mfy=fy)
      END IF
      DO j=js,je
        DO i=is,ie
          w_tl(i, j) = delp_tl(i, j)*w(i, j) + delp(i, j)*w_tl(i, j) + &
&           rarea(i, j)*(gx_tl(i, j)-gx_tl(i+1, j)+gy_tl(i, j)-gy_tl(i, &
&           j+1))
          w(i, j) = delp(i, j)*w(i, j) + (gx(i, j)-gx(i+1, j)+(gy(i, j)-&
&           gy(i, j+1)))*rarea(i, j)
        END DO
      END DO
    ELSE
      gx_tl = 0.0
      gy_tl = 0.0
      dw_tl = 0.0
      wk_tl = 0.0
    END IF
!    if ( inline_q .and. zvir>0.01 ) then
!       do j=jsd,jed
!          do i=isd,ied
!             pt(i,j) = pt(i,j)/(1.+zvir*q(i,j,k,sphum))
!          enddo
!       enddo
!    endif
    IF (hord_tm .EQ. hord_tm_pert .AND. (.NOT.split_damp)) THEN
      CALL FV_TP_2D_TLM(pt, pt_tl, crx_adv, crx_adv_tl, cry_adv, &
&                    cry_adv_tl, npx, npy, hord_tm, gx, gx_tl, gy, gy_tl&
&                    , xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, &
&                    gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl, mfx=&
&                    fx, mfx_tl=fx_tl, mfy=fy, mfy_tl=fy_tl, mass=delp, &
&                    mass_tl=delp_tl, nord=nord_t, damp_c=damp_t)
    ELSE
      CALL FV_TP_2D_TLM(pt, pt_tl, crx_adv, crx_adv_tl, cry_adv, &
&                 cry_adv_tl, npx, npy, hord_tm_pert, gx, gx_tl, gy, &
&                 gy_tl, xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, &
&                 gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl, mfx=fx, &
&                 mfx_tl=fx_tl, mfy=fy, mfy_tl=fy_tl, mass=delp, mass_tl&
&                 =delp_tl, nord=nord_t_pert, damp_c=damp_t_pert)
      call fv_tp_2d(pt, crx_adv,cry_adv, npx, npy, hord_tm, gx, gy,  &
                    xfx_adv,yfx_adv, gridstruct, bd, ra_x, ra_y,     &
                    mfx=fx, mfy=fy, mass=delp, nord=nord_t, damp_c=damp_t)
    END IF
    IF (inline_q) THEN
      DO j=js,je
        DO i=is,ie
          wk_tl(i, j) = delp_tl(i, j)
          wk(i, j) = delp(i, j)
          delp_tl(i, j) = wk_tl(i, j) + rarea(i, j)*(fx_tl(i, j)-fx_tl(i&
&           +1, j)+fy_tl(i, j)-fy_tl(i, j+1))
          delp(i, j) = wk(i, j) + (fx(i, j)-fx(i+1, j)+(fy(i, j)-fy(i, j&
&           +1)))*rarea(i, j)
          pt_tl(i, j) = ((pt_tl(i, j)*wk(i, j)+pt(i, j)*wk_tl(i, j)+&
&           rarea(i, j)*(gx_tl(i, j)-gx_tl(i+1, j)+gy_tl(i, j)-gy_tl(i, &
&           j+1)))*delp(i, j)-(pt(i, j)*wk(i, j)+(gx(i, j)-gx(i+1, j)+(&
&           gy(i, j)-gy(i, j+1)))*rarea(i, j))*delp_tl(i, j))/delp(i, j)&
&           **2
          pt(i, j) = (pt(i, j)*wk(i, j)+(gx(i, j)-gx(i+1, j)+(gy(i, j)-&
&           gy(i, j+1)))*rarea(i, j))/delp(i, j)
        END DO
      END DO
      DO iq=1,nq
        IF (hord_tr .EQ. hord_tr_pert) THEN
          CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:ied&
&                        , jsd:jed, k, iq), crx_adv, crx_adv_tl, cry_adv&
&                        , cry_adv_tl, npx, npy, hord_tr, gx, gx_tl, gy&
&                        , gy_tl, xfx_adv, xfx_adv_tl, yfx_adv, &
&                        yfx_adv_tl, gridstruct, bd, ra_x, ra_x_tl, ra_y&
&                        , ra_y_tl, mfx=fx, mfx_tl=fx_tl, mfy=fy, mfy_tl&
&                        =fy_tl, mass=delp, mass_tl=delp_tl, nord=nord_t&
&                        , damp_c=damp_t)
        ELSE
          CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:ied, &
&                     jsd:jed, k, iq), crx_adv, crx_adv_tl, cry_adv, &
&                     cry_adv_tl, npx, npy, hord_tr_pert, gx, gx_tl, gy&
&                     , gy_tl, xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl&
&                     , gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl, &
&                     mfx=fx, mfx_tl=fx_tl, mfy=fy, mfy_tl=fy_tl, mass=&
&                     delp, mass_tl=delp_tl, nord=nord_t_pert, damp_c=&
&                     damp_t_pert)
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), crx_adv,cry_adv, npx, npy, hord_tr, gx, gy,  &
                    xfx_adv,yfx_adv, gridstruct, bd, ra_x, ra_y,     &
                    mfx=fx, mfy=fy, mass=delp, nord=nord_t, damp_c=damp_t)
        END IF
        DO j=js,je
          DO i=is,ie
            q_tl(i, j, k, iq) = ((q_tl(i, j, k, iq)*wk(i, j)+q(i, j, k, &
&             iq)*wk_tl(i, j)+rarea(i, j)*(gx_tl(i, j)-gx_tl(i+1, j)+&
&             gy_tl(i, j)-gy_tl(i, j+1)))*delp(i, j)-(q(i, j, k, iq)*wk(&
&             i, j)+(gx(i, j)-gx(i+1, j)+(gy(i, j)-gy(i, j+1)))*rarea(i&
&             , j))*delp_tl(i, j))/delp(i, j)**2
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
          pt_tl(i, j) = pt_tl(i, j)*delp(i, j) + pt(i, j)*delp_tl(i, j) &
&           + rarea(i, j)*(gx_tl(i, j)-gx_tl(i+1, j)+gy_tl(i, j)-gy_tl(i&
&           , j+1))
          pt(i, j) = pt(i, j)*delp(i, j) + (gx(i, j)-gx(i+1, j)+(gy(i, j&
&           )-gy(i, j+1)))*rarea(i, j)
          delp_tl(i, j) = delp_tl(i, j) + rarea(i, j)*(fx_tl(i, j)-fx_tl&
&           (i+1, j)+fy_tl(i, j)-fy_tl(i, j+1))
          delp(i, j) = delp(i, j) + (fx(i, j)-fx(i+1, j)+(fy(i, j)-fy(i&
&           , j+1)))*rarea(i, j)
          pt_tl(i, j) = (pt_tl(i, j)*delp(i, j)-pt(i, j)*delp_tl(i, j))/&
&           delp(i, j)**2
          pt(i, j) = pt(i, j)/delp(i, j)
        END DO
      END DO
    END IF
    IF (fpp%fpp_overload_r4) THEN
      DO j=js,je
        DO i=is,ie
          dpx_tl(i, j) = dpx_tl(i, j) + rarea(i, j)*(fx_tl(i, j)-fx_tl(i&
&           +1, j)+fy_tl(i, j)-fy_tl(i, j+1))
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
        vb_tl = 0.0
        DO j=js2,je1
          DO i=is2,ie1
            vb_tl(i, j) = dt5*rsina(i, j)*(vc_tl(i-1, j)+vc_tl(i, j)-&
&             cosa(i, j)*(uc_tl(i, j-1)+uc_tl(i, j)))
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j)-(uc(i, j-1)+uc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
        END DO
      ELSE
        IF (js .EQ. 1) THEN
          vb_tl = 0.0
          DO i=is,ie+1
! corner values are incorrect
            vb_tl(i, 1) = dt5*(vt_tl(i-1, 1)+vt_tl(i, 1))
            vb(i, 1) = dt5*(vt(i-1, 1)+vt(i, 1))
          END DO
        ELSE
          vb_tl = 0.0
        END IF
        DO j=js2,je1
          DO i=is2,ie1
            vb_tl(i, j) = dt5*rsina(i, j)*(vc_tl(i-1, j)+vc_tl(i, j)-&
&             cosa(i, j)*(uc_tl(i, j-1)+uc_tl(i, j)))
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j)-(uc(i, j-1)+uc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
          IF (is .EQ. 1) THEN
! 2-pt extrapolation from both sides:
            vb_tl(1, j) = dt4*(3.*(vt_tl(0, j)+vt_tl(1, j))-vt_tl(-1, j)&
&             -vt_tl(2, j))
            vb(1, j) = dt4*(-vt(-1, j)+3.*(vt(0, j)+vt(1, j))-vt(2, j))
          END IF
          IF (ie + 1 .EQ. npx) THEN
! 2-pt extrapolation from both sides:
            vb_tl(npx, j) = dt4*(3.*(vt_tl(npx-1, j)+vt_tl(npx, j))-&
&             vt_tl(npx-2, j)-vt_tl(npx+1, j))
            vb(npx, j) = dt4*(-vt(npx-2, j)+3.*(vt(npx-1, j)+vt(npx, j))&
&             -vt(npx+1, j))
          END IF
        END DO
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb_tl(i, npy) = dt5*(vt_tl(i-1, npy)+vt_tl(i, npy))
            vb(i, npy) = dt5*(vt(i-1, npy)+vt(i, npy))
          END DO
        END IF
      END IF
    ELSE
      vb_tl = 0.0
      DO j=js,je+1
        DO i=is,ie+1
          vb_tl(i, j) = dt5*(vc_tl(i-1, j)+vc_tl(i, j))
          vb(i, j) = dt5*(vc(i-1, j)+vc(i, j))
        END DO
      END DO
    END IF
    IF (hord_mt .EQ. hord_mt_pert) THEN
      CALL YTP_V_TLM(is, ie, js, je, isd, ied, jsd, jed, vb, vb_tl, u&
&                 , v, v_tl, ub, ub_tl, hord_mt, gridstruct%dy, &
&                 gridstruct%rdy, npx, npy, flagstruct%grid_type, nested&
&                )
      ke_tl = 0.0
    ELSE
      CALL YTP_V_TLM(is, ie, js, je, isd, ied, jsd, jed, vb, vb_tl, u, v&
&              , v_tl, ub, ub_tl, hord_mt_pert, gridstruct%dy, &
&              gridstruct%rdy, npx, npy, flagstruct%grid_type, nested)
      call ytp_v(is,ie,js,je,isd,ied,jsd,jed, vb, u, v, ub, hord_mt, gridstruct%dy, gridstruct%rdy, &
                 npx, npy, flagstruct%grid_type, nested)
      ke_tl = 0.0
    END IF
    DO j=js,je+1
      DO i=is,ie+1
        ke_tl(i, j) = vb_tl(i, j)*ub(i, j) + vb(i, j)*ub_tl(i, j)
        ke(i, j) = vb(i, j)*ub(i, j)
      END DO
    END DO
    IF (flagstruct%grid_type .LT. 3) THEN
      IF (nested) THEN
        DO j=js,je+1
          DO i=is2,ie1
            ub_tl(i, j) = dt5*rsina(i, j)*(uc_tl(i, j-1)+uc_tl(i, j)-&
&             cosa(i, j)*(vc_tl(i-1, j)+vc_tl(i, j)))
            ub(i, j) = dt5*(uc(i, j-1)+uc(i, j)-(vc(i-1, j)+vc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
        END DO
      ELSE
        IF (is .EQ. 1) THEN
          DO j=js,je+1
! corner values are incorrect
            ub_tl(1, j) = dt5*(ut_tl(1, j-1)+ut_tl(1, j))
            ub(1, j) = dt5*(ut(1, j-1)+ut(1, j))
          END DO
        END IF
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is2,ie1
! 2-pt extrapolation from both sides:
              ub_tl(i, j) = dt4*(3.*(ut_tl(i, j-1)+ut_tl(i, j))-ut_tl(i&
&               , j-2)-ut_tl(i, j+1))
              ub(i, j) = dt4*(-ut(i, j-2)+3.*(ut(i, j-1)+ut(i, j))-ut(i&
&               , j+1))
            END DO
          ELSE
            DO i=is2,ie1
              ub_tl(i, j) = dt5*rsina(i, j)*(uc_tl(i, j-1)+uc_tl(i, j)-&
&               cosa(i, j)*(vc_tl(i-1, j)+vc_tl(i, j)))
              ub(i, j) = dt5*(uc(i, j-1)+uc(i, j)-(vc(i-1, j)+vc(i, j))*&
&               cosa(i, j))*rsina(i, j)
            END DO
          END IF
        END DO
        IF (ie + 1 .EQ. npx) THEN
          DO j=js,je+1
! corner values are incorrect
            ub_tl(npx, j) = dt5*(ut_tl(npx, j-1)+ut_tl(npx, j))
            ub(npx, j) = dt5*(ut(npx, j-1)+ut(npx, j))
          END DO
        END IF
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          ub_tl(i, j) = dt5*(uc_tl(i, j-1)+uc_tl(i, j))
          ub(i, j) = dt5*(uc(i, j-1)+uc(i, j))
        END DO
      END DO
    END IF
    IF (hord_mt .EQ. hord_mt_pert) THEN
      CALL XTP_U_TLM(is, ie, js, je, isd, ied, jsd, jed, ub, ub_tl, u&
&                 , u_tl, v, vb, vb_tl, hord_mt, gridstruct%dx, &
&                 gridstruct%rdx, npx, npy, flagstruct%grid_type, nested&
&                )
    ELSE
      CALL XTP_U_TLM(is, ie, js, je, isd, ied, jsd, jed, ub, ub_tl, u, &
&              u_tl, v, vb, vb_tl, hord_mt_pert, gridstruct%dx, &
&              gridstruct%rdx, npx, npy, flagstruct%grid_type, nested)
      call xtp_u(is,ie,js,je, isd,ied,jsd,jed, ub, u, v, vb, hord_mt, gridstruct%dx, gridstruct%rdx, &
                 npx, npy, flagstruct%grid_type, nested)
    END IF
    DO j=js,je+1
      DO i=is,ie+1
        ke_tl(i, j) = 0.5*(ke_tl(i, j)+ub_tl(i, j)*vb(i, j)+ub(i, j)*&
&         vb_tl(i, j))
        ke(i, j) = 0.5*(ke(i, j)+ub(i, j)*vb(i, j))
      END DO
    END DO
!-----------------------------------------
! Fix KE at the 4 corners of the face:
!-----------------------------------------
    IF (.NOT.nested) THEN
      dt6 = dt/6.
      IF (sw_corner) THEN
        ke_tl(1, 1) = dt6*((ut_tl(1, 1)+ut_tl(1, 0))*u(1, 1)+(ut(1, 1)+&
&         ut(1, 0))*u_tl(1, 1)+(vt_tl(1, 1)+vt_tl(0, 1))*v(1, 1)+(vt(1, &
&         1)+vt(0, 1))*v_tl(1, 1)+(ut_tl(1, 1)+vt_tl(1, 1))*u(0, 1)+(ut(&
&         1, 1)+vt(1, 1))*u_tl(0, 1))
        ke(1, 1) = dt6*((ut(1, 1)+ut(1, 0))*u(1, 1)+(vt(1, 1)+vt(0, 1))*&
&         v(1, 1)+(ut(1, 1)+vt(1, 1))*u(0, 1))
      END IF
      IF (se_corner) THEN
!i = npx
        ke_tl(npx, 1) = dt6*((ut_tl(npx, 1)+ut_tl(npx, 0))*u(npx-1, 1)+(&
&         ut(npx, 1)+ut(npx, 0))*u_tl(npx-1, 1)+(vt_tl(npx, 1)+vt_tl(npx&
&         -1, 1))*v(npx, 1)+(vt(npx, 1)+vt(npx-1, 1))*v_tl(npx, 1)+(&
&         ut_tl(npx, 1)-vt_tl(npx-1, 1))*u(npx, 1)+(ut(npx, 1)-vt(npx-1&
&         , 1))*u_tl(npx, 1))
        ke(npx, 1) = dt6*((ut(npx, 1)+ut(npx, 0))*u(npx-1, 1)+(vt(npx, 1&
&         )+vt(npx-1, 1))*v(npx, 1)+(ut(npx, 1)-vt(npx-1, 1))*u(npx, 1))
      END IF
      IF (ne_corner) THEN
!i = npx;      j = npy
        ke_tl(npx, npy) = dt6*((ut_tl(npx, npy)+ut_tl(npx, npy-1))*u(npx&
&         -1, npy)+(ut(npx, npy)+ut(npx, npy-1))*u_tl(npx-1, npy)+(vt_tl&
&         (npx, npy)+vt_tl(npx-1, npy))*v(npx, npy-1)+(vt(npx, npy)+vt(&
&         npx-1, npy))*v_tl(npx, npy-1)+(ut_tl(npx, npy-1)+vt_tl(npx-1, &
&         npy))*u(npx, npy)+(ut(npx, npy-1)+vt(npx-1, npy))*u_tl(npx, &
&         npy))
        ke(npx, npy) = dt6*((ut(npx, npy)+ut(npx, npy-1))*u(npx-1, npy)+&
&         (vt(npx, npy)+vt(npx-1, npy))*v(npx, npy-1)+(ut(npx, npy-1)+vt&
&         (npx-1, npy))*u(npx, npy))
      END IF
      IF (nw_corner) THEN
!j = npy
        ke_tl(1, npy) = dt6*((ut_tl(1, npy)+ut_tl(1, npy-1))*u(1, npy)+(&
&         ut(1, npy)+ut(1, npy-1))*u_tl(1, npy)+(vt_tl(1, npy)+vt_tl(0, &
&         npy))*v(1, npy-1)+(vt(1, npy)+vt(0, npy))*v_tl(1, npy-1)+(&
&         ut_tl(1, npy-1)-vt_tl(1, npy))*u(0, npy)+(ut(1, npy-1)-vt(1, &
&         npy))*u_tl(0, npy))
        ke(1, npy) = dt6*((ut(1, npy)+ut(1, npy-1))*u(1, npy)+(vt(1, npy&
&         )+vt(0, npy))*v(1, npy-1)+(ut(1, npy-1)-vt(1, npy))*u(0, npy))
      END IF
    END IF
! Compute vorticity:
    DO j=jsd,jed+1
      DO i=isd,ied
        vt_tl(i, j) = dx(i, j)*u_tl(i, j)
        vt(i, j) = u(i, j)*dx(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied+1
        ut_tl(i, j) = dy(i, j)*v_tl(i, j)
        ut(i, j) = v(i, j)*dy(i, j)
      END DO
    END DO
! wk is "volume-mean" relative vorticity
    DO j=jsd,jed
      DO i=isd,ied
        wk_tl(i, j) = rarea(i, j)*(vt_tl(i, j)-vt_tl(i, j+1)+ut_tl(i+1, &
&         j)-ut_tl(i, j))
        wk(i, j) = rarea(i, j)*(vt(i, j)-vt(i, j+1)+(ut(i+1, j)-ut(i, j)&
&         ))
      END DO
    END DO
    IF (.NOT.hydrostatic) THEN
      IF (.NOT.flagstruct%do_f3d) THEN
        DO j=js,je
          DO i=is,ie
            w_tl(i, j) = (w_tl(i, j)*delp(i, j)-w(i, j)*delp_tl(i, j))/&
&             delp(i, j)**2
            w(i, j) = w(i, j)/delp(i, j)
          END DO
        END DO
      END IF
      IF (damp_w .GT. 1.e-5) THEN
        DO j=js,je
          DO i=is,ie
            w_tl(i, j) = w_tl(i, j) + dw_tl(i, j)
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
      CALL COMPUTE_DIVERGENCE_DAMPING_TLM(nord, d2_bg, d4_bg, dddmp, &
&                                      dt, vort, vort_tl, ptc, ptc_tl, &
&                                      delpc, delpc_tl, ke, ke_tl, u, &
&                                      u_tl, v, v_tl, uc, uc_tl, vc, &
&                                      vc_tl, ua, ua_tl, va, va_tl, &
&                                      divg_d, divg_d_tl, wk, wk_tl, &
&                                      gridstruct, flagstruct, bd)
    ELSE
      wk_tj = wk
      vort_tj = vort
      delpc_tj = delpc
      ptc_tj = ptc
      ke_tj = ke
      vc_tj = vc
      uc_tj = uc
      divg_d_tj = divg_d
      CALL COMPUTE_DIVERGENCE_DAMPING_TLM(nord_pert, d2_bg_pert, &
&                                   d4_bg_pert, dddmp_pert, dt, vort_tj, &
&                                   vort_tl, ptc_tj, ptc_tl, delpc_tj, &
&                                   delpc_tl, ke_tj, ke_tl, u, u_tl, v, &
&                                   v_tl, uc_tj, uc_tl, vc_tj, vc_tl, ua, &
&                                   ua_tl, va, va_tl, divg_d_tj, divg_d_tl&
&                                   , wk_tj, wk_tl, gridstruct, flagstruct&
&                                   , bd)
      call compute_divergence_damping( nord,d2_bg,d4_bg,dddmp,dt, &
                              vort,ptc,delpc,ke,u,v,uc,vc,ua,va,divg_d,wk, &
                              gridstruct, flagstruct, bd)
    END IF
    IF (d_con .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          ub_tl(i, j) = vort_tl(i, j) - vort_tl(i+1, j)
          ub(i, j) = vort(i, j) - vort(i+1, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          vb_tl(i, j) = vort_tl(i, j) - vort_tl(i, j+1)
          vb(i, j) = vort(i, j) - vort(i, j+1)
        END DO
      END DO
    END IF
! Vorticity transport
    IF (hydrostatic) THEN
      DO j=jsd,jed
        DO i=isd,ied
          vort_tl(i, j) = wk_tl(i, j)
          vort(i, j) = wk(i, j) + f0(i, j)
        END DO
      END DO
    ELSE IF (flagstruct%do_f3d) THEN
      DO j=jsd,jed
        DO i=isd,ied
          vort_tl(i, j) = wk_tl(i, j) + f0(i, j)*z_rat_tl(i, j)
          vort(i, j) = wk(i, j) + f0(i, j)*z_rat(i, j)
        END DO
      END DO
    ELSE
      DO j=jsd,jed
        DO i=isd,ied
          vort_tl(i, j) = wk_tl(i, j)
          vort(i, j) = wk(i, j) + f0(i, j)
        END DO
      END DO
    END IF
    IF (hord_vt .EQ. hord_vt_pert) THEN
      CALL FV_TP_2D_TLM(vort, vort_tl, crx_adv, crx_adv_tl, cry_adv, &
&                    cry_adv_tl, npx, npy, hord_vt, fx, fx_tl, fy, fy_tl&
&                    , xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, &
&                    gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl)
    ELSE
      CALL FV_TP_2D_TLM(vort, vort_tl, crx_adv, crx_adv_tl, cry_adv, &
&                 cry_adv_tl, npx, npy, hord_vt_pert, fx, fx_tl, fy, &
&                 fy_tl, xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, &
&                 gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl)
      call fv_tp_2d(vort, crx_adv, cry_adv, npx, npy, hord_vt, fx, fy, &
                    xfx_adv,yfx_adv, gridstruct, bd, ra_x, ra_y)
    END IF
    DO j=js,je+1
      DO i=is,ie
        u_tl(i, j) = vt_tl(i, j) + ke_tl(i, j) - ke_tl(i+1, j) + fy_tl(i&
&         , j)
        u(i, j) = vt(i, j) + (ke(i, j)-ke(i+1, j)) + fy(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie+1
        v_tl(i, j) = ut_tl(i, j) + ke_tl(i, j) - ke_tl(i, j+1) - fx_tl(i&
&         , j)
        v(i, j) = ut(i, j) + (ke(i, j)-ke(i, j+1)) - fx(i, j)
      END DO
    END DO
!--------------------------------------------------------
! damping applied to relative vorticity (wk):
    IF (damp_v .GT. 1.e-5) THEN
      pwx1 = damp_v*gridstruct%da_min_c
      pwy1 = nord_v + 1
      damp4 = pwx1**pwy1
      call del6_vt_flux(nord_v, npx, npy, damp4, wk, vort, ut, vt, gridstruct, bd)
    END IF
    IF (damp_v_pert .GT. 1.e-5) THEN
      wk_tj = wk
      vort_tj = vort
      ut_tj = ut
      vt_tj = vt
      pwx1 = damp_v_pert*gridstruct%da_min_c
      pwy1 = nord_v_pert + 1
      damp4 = pwx1**pwy1
      CALL DEL6_VT_FLUX_TLM(nord_v_pert, npx, npy, damp4, wk_tj, wk_tl, vort_tj, &
&                     vort_tl, ut_tj, ut_tl, vt_tj, vt_tl, gridstruct, bd)
    END IF
    IF (d_con .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          ub_tl(i, j) = rdx(i, j)*(ub_tl(i, j)+vt_tl(i, j))
          ub(i, j) = (ub(i, j)+vt(i, j))*rdx(i, j)
          fy_tl(i, j) = rdx(i, j)*u_tl(i, j)
          fy(i, j) = u(i, j)*rdx(i, j)
          gy_tl(i, j) = fy_tl(i, j)*ub(i, j) + fy(i, j)*ub_tl(i, j)
          gy(i, j) = fy(i, j)*ub(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          vb_tl(i, j) = rdy(i, j)*(vb_tl(i, j)-ut_tl(i, j))
          vb(i, j) = (vb(i, j)-ut(i, j))*rdy(i, j)
          fx_tl(i, j) = rdy(i, j)*v_tl(i, j)
          fx(i, j) = v(i, j)*rdy(i, j)
          gx_tl(i, j) = fx_tl(i, j)*vb(i, j) + fx(i, j)*vb_tl(i, j)
          gx(i, j) = fx(i, j)*vb(i, j)
        END DO
      END DO
!----------------------------------
! Heating due to damping:
!----------------------------------
      damp = 0.25*d_con
      DO j=js,je
        DO i=is,ie
          u2_tl = fy_tl(i, j) + fy_tl(i, j+1)
          u2 = fy(i, j) + fy(i, j+1)
          du2_tl = ub_tl(i, j) + ub_tl(i, j+1)
          du2 = ub(i, j) + ub(i, j+1)
          v2_tl = fx_tl(i, j) + fx_tl(i+1, j)
          v2 = fx(i, j) + fx(i+1, j)
          dv2_tl = vb_tl(i, j) + vb_tl(i+1, j)
          dv2 = vb(i, j) + vb(i+1, j)
! Total energy conserving:
! Convert lost KE due to divergence damping to "heat"
          heat_source_tl(i, j) = delp_tl(i, j)*(heat_source(i, j)-damp*&
&           rsin2(i, j)*(ub(i, j)**2+ub(i, j+1)**2+vb(i, j)**2+vb(i+1, j&
&           )**2+2.*(gy(i, j)+gy(i, j+1)+gx(i, j)+gx(i+1, j))-cosa_s(i, &
&           j)*(u2*dv2+v2*du2+du2*dv2))) + delp(i, j)*(heat_source_tl(i&
&           , j)-damp*rsin2(i, j)*(2*ub(i, j)*ub_tl(i, j)+2*ub(i, j+1)*&
&           ub_tl(i, j+1)+2*vb(i, j)*vb_tl(i, j)+2*vb(i+1, j)*vb_tl(i+1&
&           , j)+2.*(gy_tl(i, j)+gy_tl(i, j+1)+gx_tl(i, j)+gx_tl(i+1, j)&
&           )-cosa_s(i, j)*(u2_tl*dv2+u2*dv2_tl+v2_tl*du2+v2*du2_tl+&
&           du2_tl*dv2+du2*dv2_tl)))
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
    IF (damp_v_pert .GT. 1.e-5) THEN
      DO j=js,je+1
        DO i=is,ie
          u_tl(i, j) = u_tl(i, j) + vt_tl(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v_tl(i, j) = v_tl(i, j) - ut_tl(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE D_SW_TLM
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
    REAL*8 :: pwx1
    INTEGER :: pwy1
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
      call fv_tp_2d(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx, fy,  &
                    xfx_adv,yfx_adv, gridstruct, bd, ra_x, ra_y, nord=nord_v, damp_c=damp_v)
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
        pwx1 = damp_w*gridstruct%da_min_c
        pwy1 = nord_w + 1
        damp4 = pwx1**pwy1
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
      call fv_tp_2d(w, crx_adv,cry_adv, npx, npy, hord_vt, gx, gy, xfx_adv, yfx_adv, &
                    gridstruct, bd, ra_x, ra_y, mfx=fx, mfy=fy)
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
      call fv_tp_2d(pt, crx_adv,cry_adv, npx, npy, hord_tm, gx, gy,  &
                    xfx_adv,yfx_adv, gridstruct, bd, ra_x, ra_y,     &
                    mfx=fx, mfy=fy, mass=delp, nord=nord_t, damp_c=damp_t)
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
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), crx_adv,cry_adv, npx, npy, hord_tr, gx, gy,  &
                    xfx_adv,yfx_adv, gridstruct, bd, ra_x, ra_y,     &
                    mfx=fx, mfy=fy, mass=delp, nord=nord_t, damp_c=damp_t)
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
      call ytp_v(is,ie,js,je,isd,ied,jsd,jed, vb, u, v, ub, hord_mt, gridstruct%dy, gridstruct%rdy, &
                 npx, npy, flagstruct%grid_type, nested)
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
      call xtp_u(is,ie,js,je, isd,ied,jsd,jed, ub, u, v, vb, hord_mt, gridstruct%dx, gridstruct%rdx, &
                 npx, npy, flagstruct%grid_type, nested)
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
      call compute_divergence_damping( nord,d2_bg,d4_bg,dddmp,dt, &
                              vort,ptc,delpc,ke,u,v,uc,vc,ua,va,divg_d,wk, &
                              gridstruct, flagstruct, bd)
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
      call fv_tp_2d(vort, crx_adv, cry_adv, npx, npy, hord_vt, fx, fy, &
                    xfx_adv,yfx_adv, gridstruct, bd, ra_x, ra_y)
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
      pwx1 = damp_v*gridstruct%da_min_c
      pwy1 = nord_v + 1
      damp4 = pwx1**pwy1
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
!  Differentiation of del6_vt_flux in forward (tangent) mode:
!   variations   of useful results: fy2 d2 fx2
!   with respect to varying inputs: q fy2 d2 fx2
  SUBROUTINE DEL6_VT_FLUX_TLM(nord, npx, npy, damp, q, q_tl, d2, d2_tl, &
&   fx2, fx2_tl, fy2, fy2_tl, gridstruct, bd)
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
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! Work arrays:
    REAL, INTENT(OUT) :: d2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(OUT) :: d2_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(OUT) :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2(bd%isd&
&   :bd%ied, bd%jsd:bd%jed+1)
    REAL, INTENT(OUT) :: fx2_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2_tl(&
&   bd%isd:bd%ied, bd%jsd:bd%jed+1)
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
        d2_tl(i, j) = damp*q_tl(i, j)
        d2(i, j) = damp*q(i, j)
      END DO
    END DO
    IF (nord .GT. 0) CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 1, &
&                                    nested, bd, gridstruct%sw_corner, &
&                                    gridstruct%se_corner, gridstruct%&
&                                    nw_corner, gridstruct%ne_corner)
    DO j=js-nord,je+nord
      DO i=is-nord,ie+nord+1
        fx2_tl(i, j) = gridstruct%del6_v(i, j)*(d2_tl(i-1, j)-d2_tl(i, j&
&         ))
        fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i-1, j)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 2, &
&                                    nested, bd, gridstruct%sw_corner, &
&                                    gridstruct%se_corner, gridstruct%&
&                                    nw_corner, gridstruct%ne_corner)
    DO j=js-nord,je+nord+1
      DO i=is-nord,ie+nord
        fy2_tl(i, j) = gridstruct%del6_u(i, j)*(d2_tl(i, j-1)-d2_tl(i, j&
&         ))
        fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j-1)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) THEN
      DO n=1,nord
        nt = nord - n
        DO j=js-nt-1,je+nt+1
          DO i=is-nt-1,ie+nt+1
            d2_tl(i, j) = gridstruct%rarea(i, j)*(fx2_tl(i, j)-fx2_tl(i+&
&             1, j)+fy2_tl(i, j)-fy2_tl(i, j+1))
            d2(i, j) = (fx2(i, j)-fx2(i+1, j)+(fy2(i, j)-fy2(i, j+1)))*&
&             gridstruct%rarea(i, j)
          END DO
        END DO
        CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 1, nested, bd, &
&                       gridstruct%sw_corner, gridstruct%se_corner, &
&                       gridstruct%nw_corner, gridstruct%ne_corner)
        DO j=js-nt,je+nt
          DO i=is-nt,ie+nt+1
            fx2_tl(i, j) = gridstruct%del6_v(i, j)*(d2_tl(i, j)-d2_tl(i-&
&             1, j))
            fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i, j)-d2(i-1, j))
          END DO
        END DO
        CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 2, nested, bd, &
&                       gridstruct%sw_corner, gridstruct%se_corner, &
&                       gridstruct%nw_corner, gridstruct%ne_corner)
        DO j=js-nt,je+nt+1
          DO i=is-nt,ie+nt
            fy2_tl(i, j) = gridstruct%del6_u(i, j)*(d2_tl(i, j)-d2_tl(i&
&             , j-1))
            fy2(i, j) = gridstruct%del6_u(i, j)*(d2(i, j)-d2(i, j-1))
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE DEL6_VT_FLUX_TLM
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
!  Differentiation of divergence_corner in forward (tangent) mode:
!   variations   of useful results: divg_d
!   with respect to varying inputs: u v ua va divg_d
  SUBROUTINE DIVERGENCE_CORNER_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, &
&   va_tl, divg_d, divg_d_tl, gridstruct, flagstruct, bd)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua_tl, &
&   va_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   divg_d
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   divg_d_tl
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! local
    REAL :: uf(bd%is-2:bd%ie+2, bd%js-1:bd%je+2)
    REAL :: uf_tl(bd%is-2:bd%ie+2, bd%js-1:bd%je+2)
    REAL :: vf(bd%is-1:bd%ie+2, bd%js-2:bd%je+2)
    REAL :: vf_tl(bd%is-1:bd%ie+2, bd%js-2:bd%je+2)
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
      uf_tl = 0.0
      DO j=js-1,je+2
        DO i=is-2,ie+2
          uf_tl(i, j) = dyc(i, j)*u_tl(i, j)
          uf(i, j) = u(i, j)*dyc(i, j)
        END DO
      END DO
      vf_tl = 0.0
      DO j=js-2,je+2
        DO i=is-1,ie+2
          vf_tl(i, j) = dxc(i, j)*v_tl(i, j)
          vf(i, j) = v(i, j)*dxc(i, j)
        END DO
      END DO
      DO j=js-1,je+2
        DO i=is-1,ie+2
          divg_d_tl(i, j) = gridstruct%rarea_c(i, j)*(vf_tl(i, j-1)-&
&           vf_tl(i, j)+uf_tl(i-1, j)-uf_tl(i, j))
          divg_d(i, j) = gridstruct%rarea_c(i, j)*(vf(i, j-1)-vf(i, j)+(&
&           uf(i-1, j)-uf(i, j)))
        END DO
      END DO
    ELSE
      uf_tl = 0.0
!     9---4---8
!     |       |
!     1   5   3
!     |       |
!     6---2---7
      DO j=js,je+1
        IF (j .EQ. 1 .OR. j .EQ. npy) THEN
          DO i=is-1,ie+1
            uf_tl(i, j) = dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+sin_sg(i, j, &
&             2))*u_tl(i, j)
            uf(i, j) = u(i, j)*dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+sin_sg(i&
&             , j, 2))
          END DO
        ELSE
          DO i=is-1,ie+1
            uf_tl(i, j) = dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+sin_sg(i, j, &
&             2))*(u_tl(i, j)-0.25*(cos_sg(i, j-1, 4)+cos_sg(i, j, 2))*(&
&             va_tl(i, j-1)+va_tl(i, j)))
            uf(i, j) = (u(i, j)-0.25*(va(i, j-1)+va(i, j))*(cos_sg(i, j-&
&             1, 4)+cos_sg(i, j, 2)))*dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+&
&             sin_sg(i, j, 2))
          END DO
        END IF
      END DO
      vf_tl = 0.0
      DO j=js-1,je+1
        DO i=is2,ie1
          vf_tl(i, j) = dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+sin_sg(i, j, 1)&
&           )*(v_tl(i, j)-0.25*(cos_sg(i-1, j, 3)+cos_sg(i, j, 1))*(&
&           ua_tl(i-1, j)+ua_tl(i, j)))
          vf(i, j) = (v(i, j)-0.25*(ua(i-1, j)+ua(i, j))*(cos_sg(i-1, j&
&           , 3)+cos_sg(i, j, 1)))*dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+&
&           sin_sg(i, j, 1))
        END DO
        IF (is .EQ. 1) THEN
          vf_tl(1, j) = dxc(1, j)*0.5*(sin_sg(0, j, 3)+sin_sg(1, j, 1))*&
&           v_tl(1, j)
          vf(1, j) = v(1, j)*dxc(1, j)*0.5*(sin_sg(0, j, 3)+sin_sg(1, j&
&           , 1))
        END IF
        IF (ie + 1 .EQ. npx) THEN
          vf_tl(npx, j) = dxc(npx, j)*0.5*(sin_sg(npx-1, j, 3)+sin_sg(&
&           npx, j, 1))*v_tl(npx, j)
          vf(npx, j) = v(npx, j)*dxc(npx, j)*0.5*(sin_sg(npx-1, j, 3)+&
&           sin_sg(npx, j, 1))
        END IF
      END DO
      DO j=js,je+1
        DO i=is,ie+1
          divg_d_tl(i, j) = vf_tl(i, j-1) - vf_tl(i, j) + uf_tl(i-1, j) &
&           - uf_tl(i, j)
          divg_d(i, j) = vf(i, j-1) - vf(i, j) + (uf(i-1, j)-uf(i, j))
        END DO
      END DO
! Remove the extra term at the corners:
      IF (gridstruct%sw_corner) THEN
        divg_d_tl(1, 1) = divg_d_tl(1, 1) - vf_tl(1, 0)
        divg_d(1, 1) = divg_d(1, 1) - vf(1, 0)
      END IF
      IF (gridstruct%se_corner) THEN
        divg_d_tl(npx, 1) = divg_d_tl(npx, 1) - vf_tl(npx, 0)
        divg_d(npx, 1) = divg_d(npx, 1) - vf(npx, 0)
      END IF
      IF (gridstruct%ne_corner) THEN
        divg_d_tl(npx, npy) = divg_d_tl(npx, npy) + vf_tl(npx, npy)
        divg_d(npx, npy) = divg_d(npx, npy) + vf(npx, npy)
      END IF
      IF (gridstruct%nw_corner) THEN
        divg_d_tl(1, npy) = divg_d_tl(1, npy) + vf_tl(1, npy)
        divg_d(1, npy) = divg_d(1, npy) + vf(1, npy)
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          divg_d_tl(i, j) = gridstruct%rarea_c(i, j)*divg_d_tl(i, j)
          divg_d(i, j) = gridstruct%rarea_c(i, j)*divg_d(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE DIVERGENCE_CORNER_TLM
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
!  Differentiation of divergence_corner_nest in forward (tangent) mode:
!   variations   of useful results: divg_d
!   with respect to varying inputs: u v ua va
  SUBROUTINE DIVERGENCE_CORNER_NEST_TLM(u, u_tl, v, v_tl, ua, ua_tl, va&
&   , va_tl, divg_d, divg_d_tl, gridstruct, flagstruct, bd)
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
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua_tl, &
&   va_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   divg_d
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   divg_d_tl
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
! local
    REAL :: uf(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: uf_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: vf(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vf_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed)
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
      uf_tl = 0.0
      DO j=jsd,jed
        DO i=isd,ied
          uf_tl(i, j) = dyc(i, j)*u_tl(i, j)
          uf(i, j) = u(i, j)*dyc(i, j)
        END DO
      END DO
      vf_tl = 0.0
      DO j=jsd,jed
        DO i=isd,ied
          vf_tl(i, j) = dxc(i, j)*v_tl(i, j)
          vf(i, j) = v(i, j)*dxc(i, j)
        END DO
      END DO
      divg_d_tl = 0.0
      DO j=jsd+1,jed
        DO i=isd+1,ied
          divg_d_tl(i, j) = rarea_c(i, j)*(vf_tl(i, j-1)-vf_tl(i, j)+&
&           uf_tl(i-1, j)-uf_tl(i, j))
          divg_d(i, j) = rarea_c(i, j)*(vf(i, j-1)-vf(i, j)+(uf(i-1, j)-&
&           uf(i, j)))
        END DO
      END DO
    ELSE
      uf_tl = 0.0
      DO j=jsd+1,jed
        DO i=isd,ied
          uf_tl(i, j) = dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+sin_sg(i, j, 2)&
&           )*(u_tl(i, j)-0.25*(cos_sg(i, j-1, 4)+cos_sg(i, j, 2))*(&
&           va_tl(i, j-1)+va_tl(i, j)))
          uf(i, j) = (u(i, j)-0.25*(va(i, j-1)+va(i, j))*(cos_sg(i, j-1&
&           , 4)+cos_sg(i, j, 2)))*dyc(i, j)*0.5*(sin_sg(i, j-1, 4)+&
&           sin_sg(i, j, 2))
        END DO
      END DO
      vf_tl = 0.0
      DO j=jsd,jed
        DO i=isd+1,ied
          vf_tl(i, j) = dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+sin_sg(i, j, 1)&
&           )*(v_tl(i, j)-0.25*(cos_sg(i-1, j, 3)+cos_sg(i, j, 1))*(&
&           ua_tl(i-1, j)+ua_tl(i, j)))
          vf(i, j) = (v(i, j)-0.25*(ua(i-1, j)+ua(i, j))*(cos_sg(i-1, j&
&           , 3)+cos_sg(i, j, 1)))*dxc(i, j)*0.5*(sin_sg(i-1, j, 3)+&
&           sin_sg(i, j, 1))
        END DO
      END DO
      divg_d_tl = 0.0
      DO j=jsd+1,jed
        DO i=isd+1,ied
          divg_d_tl(i, j) = rarea_c(i, j)*(vf_tl(i, j-1)-vf_tl(i, j)+&
&           uf_tl(i-1, j)-uf_tl(i, j))
          divg_d(i, j) = (vf(i, j-1)-vf(i, j)+(uf(i-1, j)-uf(i, j)))*&
&           rarea_c(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE DIVERGENCE_CORNER_NEST_TLM
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
          CALL PERT_PPM(ie - is + 2, v(is:ie+1, j), bl(is:ie+1, j), br(&
&                 is:ie+1, j), -1)
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
          CALL PERT_PPM(ie - is + 2, v(is:ie+1, j), bl(is:ie+1, j), br(&
&                 is:ie+1, j), -1)
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
!  Differentiation of d2a2c_vect in forward (tangent) mode:
!   variations   of useful results: ua uc ut va vc vt
!   with respect to varying inputs: u v ua uc ut va vc vt
!There is a limit to how far this routine can fill uc and vc in the
! halo, and so either mpp_update_domains or some sort of boundary
!  routine (extrapolation, outflow, interpolation from a nested grid)
!   is needed after c_sw is completed if these variables are needed
!    in the halo
  SUBROUTINE D2A2C_VECT_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, uc, &
&   uc_tl, vc, vc_tl, ut, ut_tl, vt, vt_tl, dord4, gridstruct, bd, npx, &
&   npy, nested, grid_type)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    LOGICAL, INTENT(IN) :: dord4
    REAL, INTENT(IN) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL, INTENT(IN) :: u_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL, INTENT(IN) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: v_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(OUT) :: uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(OUT) :: &
&   uc_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(OUT) :: vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(OUT) :: &
&   vc_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: ua, va&
&   , ut, vt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: ua_tl&
&   , va_tl, ut_tl, vt_tl
    INTEGER, INTENT(IN) :: npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! Local
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: utmp, vtmp
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed) :: utmp_tl, vtmp_tl
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
      utmp_tl = 0.0
      DO j=jsd+1,jed-1
        DO i=isd,ied
          utmp_tl(i, j) = a2*(u_tl(i, j-1)+u_tl(i, j+2)) + a1*(u_tl(i, j&
&           )+u_tl(i, j+1))
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      DO i=isd,ied
!j = jsd
        utmp_tl(i, jsd) = 0.5*(u_tl(i, jsd)+u_tl(i, jsd+1))
        utmp(i, jsd) = 0.5*(u(i, jsd)+u(i, jsd+1))
!j = jed
        utmp_tl(i, jed) = 0.5*(u_tl(i, jed)+u_tl(i, jed+1))
        utmp(i, jed) = 0.5*(u(i, jed)+u(i, jed+1))
      END DO
      vtmp_tl = 0.0
      DO j=jsd,jed
        DO i=isd+1,ied-1
          vtmp_tl(i, j) = a2*(v_tl(i-1, j)+v_tl(i+2, j)) + a1*(v_tl(i, j&
&           )+v_tl(i+1, j))
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
!i = isd
        vtmp_tl(isd, j) = 0.5*(v_tl(isd, j)+v_tl(isd+1, j))
        vtmp(isd, j) = 0.5*(v(isd, j)+v(isd+1, j))
!i = ied
        vtmp_tl(ied, j) = 0.5*(v_tl(ied, j)+v_tl(ied+1, j))
        vtmp(ied, j) = 0.5*(v(ied, j)+v(ied+1, j))
      END DO
      DO j=jsd,jed
        DO i=isd,ied
          ua_tl(i, j) = rsin2(i, j)*(utmp_tl(i, j)-cosa_s(i, j)*vtmp_tl(&
&           i, j))
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          va_tl(i, j) = rsin2(i, j)*(vtmp_tl(i, j)-cosa_s(i, j)*utmp_tl(&
&           i, j))
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
        utmp_tl = 0.0
      ELSE
        min1 = npy - npt
        utmp_tl = 0.0
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
          utmp_tl(i, j) = a2*(u_tl(i, j-1)+u_tl(i, j+2)) + a1*(u_tl(i, j&
&           )+u_tl(i, j+1))
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
        vtmp_tl = 0.0
      ELSE
        min3 = npy - npt
        vtmp_tl = 0.0
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
          vtmp_tl(i, j) = a2*(v_tl(i-1, j)+v_tl(i+2, j)) + a1*(v_tl(i, j&
&           )+v_tl(i+1, j))
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
              utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
          DO j=npy-npt+1,jed
            DO i=isd,ied
              utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
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
              utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
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
              utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
      END IF
! Contra-variant components at cell center:
      DO j=js-1-id,je+1+id
        DO i=is-1-id,ie+1+id
          ua_tl(i, j) = rsin2(i, j)*(utmp_tl(i, j)-cosa_s(i, j)*vtmp_tl(&
&           i, j))
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          va_tl(i, j) = rsin2(i, j)*(vtmp_tl(i, j)-cosa_s(i, j)*utmp_tl(&
&           i, j))
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
        utmp_tl(i, 0) = -vtmp_tl(0, 1-i)
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
    END IF
    IF (gridstruct%se_corner) THEN
      DO i=0,2
        utmp_tl(npx+i, 0) = vtmp_tl(npx, i+1)
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
    END IF
    IF (gridstruct%ne_corner) THEN
      DO i=0,2
        utmp_tl(npx+i, npy) = -vtmp_tl(npx, je-i)
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
    END IF
    IF (gridstruct%nw_corner) THEN
      DO i=-2,0
        utmp_tl(i, npy) = vtmp_tl(0, je+i)
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
        uc_tl(i, j) = a2*(utmp_tl(i-2, j)+utmp_tl(i+1, j)) + a1*(utmp_tl&
&         (i-1, j)+utmp_tl(i, j))
        uc(i, j) = a2*(utmp(i-2, j)+utmp(i+1, j)) + a1*(utmp(i-1, j)+&
&         utmp(i, j))
        ut_tl(i, j) = rsin_u(i, j)*(uc_tl(i, j)-cosa_u(i, j)*v_tl(i, j))
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
! Xdir:
      IF (gridstruct%sw_corner) THEN
        ua_tl(-1, 0) = -va_tl(0, 2)
        ua(-1, 0) = -va(0, 2)
        ua_tl(0, 0) = -va_tl(0, 1)
        ua(0, 0) = -va(0, 1)
      END IF
      IF (gridstruct%se_corner) THEN
        ua_tl(npx, 0) = va_tl(npx, 1)
        ua(npx, 0) = va(npx, 1)
        ua_tl(npx+1, 0) = va_tl(npx, 2)
        ua(npx+1, 0) = va(npx, 2)
      END IF
      IF (gridstruct%ne_corner) THEN
        ua_tl(npx, npy) = -va_tl(npx, npy-1)
        ua(npx, npy) = -va(npx, npy-1)
        ua_tl(npx+1, npy) = -va_tl(npx, npy-2)
        ua(npx+1, npy) = -va(npx, npy-2)
      END IF
      IF (gridstruct%nw_corner) THEN
        ua_tl(-1, npy) = va_tl(0, npy-2)
        ua(-1, npy) = va(0, npy-2)
        ua_tl(0, npy) = va_tl(0, npy-1)
        ua(0, npy) = va(0, npy-1)
      END IF
      IF (is .EQ. 1 .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc_tl(0, j) = c1*utmp_tl(-2, j) + c2*utmp_tl(-1, j) + c3*&
&           utmp_tl(0, j)
          uc(0, j) = c1*utmp(-2, j) + c2*utmp(-1, j) + c3*utmp(0, j)
          ut_tl(1, j) = EDGE_INTERPOLATE4_TLM(ua(-1:2, j), ua_tl(-1:2, j&
&           ), dxa(-1:2, j), ut(1, j))
!Want to use the UPSTREAM value
          IF (ut(1, j) .GT. 0.) THEN
            uc_tl(1, j) = sin_sg(0, j, 3)*ut_tl(1, j)
            uc(1, j) = ut(1, j)*sin_sg(0, j, 3)
          ELSE
            uc_tl(1, j) = sin_sg(1, j, 1)*ut_tl(1, j)
            uc(1, j) = ut(1, j)*sin_sg(1, j, 1)
          END IF
          uc_tl(2, j) = c1*utmp_tl(3, j) + c2*utmp_tl(2, j) + c3*utmp_tl&
&           (1, j)
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
          ut_tl(0, j) = rsin_u(0, j)*(uc_tl(0, j)-cosa_u(0, j)*v_tl(0, j&
&           ))
          ut(0, j) = (uc(0, j)-v(0, j)*cosa_u(0, j))*rsin_u(0, j)
          ut_tl(2, j) = rsin_u(2, j)*(uc_tl(2, j)-cosa_u(2, j)*v_tl(2, j&
&           ))
          ut(2, j) = (uc(2, j)-v(2, j)*cosa_u(2, j))*rsin_u(2, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc_tl(npx-1, j) = c1*utmp_tl(npx-3, j) + c2*utmp_tl(npx-2, j) &
&           + c3*utmp_tl(npx-1, j)
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
          ut_tl(npx, j) = EDGE_INTERPOLATE4_TLM(ua(npx-2:npx+1, j), &
&           ua_tl(npx-2:npx+1, j), dxa(npx-2:npx+1, j), ut(npx, j))
          IF (ut(npx, j) .GT. 0.) THEN
            uc_tl(npx, j) = sin_sg(npx-1, j, 3)*ut_tl(npx, j)
            uc(npx, j) = ut(npx, j)*sin_sg(npx-1, j, 3)
          ELSE
            uc_tl(npx, j) = sin_sg(npx, j, 1)*ut_tl(npx, j)
            uc(npx, j) = ut(npx, j)*sin_sg(npx, j, 1)
          END IF
          uc_tl(npx+1, j) = c3*utmp_tl(npx, j) + c2*utmp_tl(npx+1, j) + &
&           c1*utmp_tl(npx+2, j)
          uc(npx+1, j) = c3*utmp(npx, j) + c2*utmp(npx+1, j) + c1*utmp(&
&           npx+2, j)
          ut_tl(npx-1, j) = rsin_u(npx-1, j)*(uc_tl(npx-1, j)-cosa_u(npx&
&           -1, j)*v_tl(npx-1, j))
          ut(npx-1, j) = (uc(npx-1, j)-v(npx-1, j)*cosa_u(npx-1, j))*&
&           rsin_u(npx-1, j)
          ut_tl(npx+1, j) = rsin_u(npx+1, j)*(uc_tl(npx+1, j)-cosa_u(npx&
&           +1, j)*v_tl(npx+1, j))
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
        vtmp_tl(0, j) = -utmp_tl(1-j, 0)
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
    END IF
    IF (gridstruct%nw_corner) THEN
      DO j=0,2
        vtmp_tl(0, npy+j) = utmp_tl(j+1, npy)
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
    END IF
    IF (gridstruct%se_corner) THEN
      DO j=-2,0
        vtmp_tl(npx, j) = utmp_tl(ie+j, 0)
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
    END IF
    IF (gridstruct%ne_corner) THEN
      DO j=0,2
        vtmp_tl(npx, npy+j) = -utmp_tl(ie-j, npy)
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
    END IF
    IF (gridstruct%sw_corner) THEN
      va_tl(0, -1) = -ua_tl(2, 0)
      va(0, -1) = -ua(2, 0)
      va_tl(0, 0) = -ua_tl(1, 0)
      va(0, 0) = -ua(1, 0)
    END IF
    IF (gridstruct%se_corner) THEN
      va_tl(npx, 0) = ua_tl(npx-1, 0)
      va(npx, 0) = ua(npx-1, 0)
      va_tl(npx, -1) = ua_tl(npx-2, 0)
      va(npx, -1) = ua(npx-2, 0)
    END IF
    IF (gridstruct%ne_corner) THEN
      va_tl(npx, npy) = -ua_tl(npx-1, npy)
      va(npx, npy) = -ua(npx-1, npy)
      va_tl(npx, npy+1) = -ua_tl(npx-2, npy)
      va(npx, npy+1) = -ua(npx-2, npy)
    END IF
    IF (gridstruct%nw_corner) THEN
      va_tl(0, npy) = ua_tl(1, npy)
      va(0, npy) = ua(1, npy)
      va_tl(0, npy+1) = ua_tl(2, npy)
      va(0, npy+1) = ua(2, npy)
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1 .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vt_tl(i, j) = EDGE_INTERPOLATE4_TLM(va(i, -1:2), va_tl(i, -1&
&             :2), dya(i, -1:2), vt(i, j))
            IF (vt(i, j) .GT. 0.) THEN
              vc_tl(i, j) = sin_sg(i, j-1, 4)*vt_tl(i, j)
              vc(i, j) = vt(i, j)*sin_sg(i, j-1, 4)
            ELSE
              vc_tl(i, j) = sin_sg(i, j, 2)*vt_tl(i, j)
              vc(i, j) = vt(i, j)*sin_sg(i, j, 2)
            END IF
          END DO
        ELSE IF (j .EQ. 0 .OR. (j .EQ. npy - 1 .AND. (.NOT.nested))) &
&       THEN
          DO i=is-1,ie+1
            vc_tl(i, j) = c1*vtmp_tl(i, j-2) + c2*vtmp_tl(i, j-1) + c3*&
&             vtmp_tl(i, j)
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
            vt_tl(i, j) = rsin_v(i, j)*(vc_tl(i, j)-cosa_v(i, j)*u_tl(i&
&             , j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        ELSE IF (j .EQ. 2 .OR. (j .EQ. npy + 1 .AND. (.NOT.nested))) &
&       THEN
          DO i=is-1,ie+1
            vc_tl(i, j) = c1*vtmp_tl(i, j+1) + c2*vtmp_tl(i, j) + c3*&
&             vtmp_tl(i, j-1)
            vc(i, j) = c1*vtmp(i, j+1) + c2*vtmp(i, j) + c3*vtmp(i, j-1)
            vt_tl(i, j) = rsin_v(i, j)*(vc_tl(i, j)-cosa_v(i, j)*u_tl(i&
&             , j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        ELSE IF (j .EQ. npy .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vt_tl(i, j) = EDGE_INTERPOLATE4_TLM(va(i, j-2:j+1), va_tl(i&
&             , j-2:j+1), dya(i, j-2:j+1), vt(i, j))
            IF (vt(i, j) .GT. 0.) THEN
              vc_tl(i, j) = sin_sg(i, j-1, 4)*vt_tl(i, j)
              vc(i, j) = vt(i, j)*sin_sg(i, j-1, 4)
            ELSE
              vc_tl(i, j) = sin_sg(i, j, 2)*vt_tl(i, j)
              vc(i, j) = vt(i, j)*sin_sg(i, j, 2)
            END IF
          END DO
        ELSE
! 4th order interpolation for interior points:
          DO i=is-1,ie+1
            vc_tl(i, j) = a2*(vtmp_tl(i, j-2)+vtmp_tl(i, j+1)) + a1*(&
&             vtmp_tl(i, j-1)+vtmp_tl(i, j))
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
            vt_tl(i, j) = rsin_v(i, j)*(vc_tl(i, j)-cosa_v(i, j)*u_tl(i&
&             , j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        END IF
      END DO
    ELSE
! 4th order interpolation:
      DO j=js-1,je+2
        DO i=is-1,ie+1
          vc_tl(i, j) = a2*(vtmp_tl(i, j-2)+vtmp_tl(i, j+1)) + a1*(&
&           vtmp_tl(i, j-1)+vtmp_tl(i, j))
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
          vt_tl(i, j) = vc_tl(i, j)
          vt(i, j) = vc(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE D2A2C_VECT_TLM
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
!  Differentiation of edge_interpolate4 in forward (tangent) mode:
!   variations   of useful results: edge_interpolate4
!   with respect to varying inputs: ua
  REAL FUNCTION EDGE_INTERPOLATE4_TLM(ua, ua_tl, dxa, edge_interpolate4)
    IMPLICIT NONE
    REAL, INTENT(IN) :: ua(4)
    REAL, INTENT(IN) :: ua_tl(4)
    REAL, INTENT(IN) :: dxa(4)
    REAL :: t1, t2
    REAL :: edge_interpolate4
    t1 = dxa(1) + dxa(2)
    t2 = dxa(3) + dxa(4)
    edge_interpolate4_tlm = 0.5*(((t1+dxa(2))*ua_tl(2)-dxa(2)*ua_tl(1))/&
&     t1+((t2+dxa(3))*ua_tl(3)-dxa(3)*ua_tl(4))/t2)
    edge_interpolate4 = 0.5*(((t1+dxa(2))*ua(2)-dxa(2)*ua(1))/t1+((t2+&
&     dxa(3))*ua(3)-dxa(3)*ua(4))/t2)
  END FUNCTION EDGE_INTERPOLATE4_TLM
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
!  Differentiation of fill2_4corners in forward (tangent) mode:
!   variations   of useful results: q1 q2
!   with respect to varying inputs: q1 q2
  SUBROUTINE FILL2_4CORNERS_TLM(q1, q1_tl, q2, q2_tl, dir, bd, npx, npy&
&   , sw_corner, se_corner, ne_corner, nw_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q1(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q1_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q2_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
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
        q1_tl(-1, 0) = q1_tl(0, 2)
        q1(-1, 0) = q1(0, 2)
        q1_tl(0, 0) = q1_tl(0, 1)
        q1(0, 0) = q1(0, 1)
        q2_tl(-1, 0) = q2_tl(0, 2)
        q2(-1, 0) = q2(0, 2)
        q2_tl(0, 0) = q2_tl(0, 1)
        q2(0, 0) = q2(0, 1)
      END IF
      IF (se_corner) THEN
        q1_tl(npx+1, 0) = q1_tl(npx, 2)
        q1(npx+1, 0) = q1(npx, 2)
        q1_tl(npx, 0) = q1_tl(npx, 1)
        q1(npx, 0) = q1(npx, 1)
        q2_tl(npx+1, 0) = q2_tl(npx, 2)
        q2(npx+1, 0) = q2(npx, 2)
        q2_tl(npx, 0) = q2_tl(npx, 1)
        q2(npx, 0) = q2(npx, 1)
      END IF
      IF (nw_corner) THEN
        q1_tl(0, npy) = q1_tl(0, npy-1)
        q1(0, npy) = q1(0, npy-1)
        q1_tl(-1, npy) = q1_tl(0, npy-2)
        q1(-1, npy) = q1(0, npy-2)
        q2_tl(0, npy) = q2_tl(0, npy-1)
        q2(0, npy) = q2(0, npy-1)
        q2_tl(-1, npy) = q2_tl(0, npy-2)
        q2(-1, npy) = q2(0, npy-2)
      END IF
      IF (ne_corner) THEN
        q1_tl(npx, npy) = q1_tl(npx, npy-1)
        q1(npx, npy) = q1(npx, npy-1)
        q1_tl(npx+1, npy) = q1_tl(npx, npy-2)
        q1(npx+1, npy) = q1(npx, npy-2)
        q2_tl(npx, npy) = q2_tl(npx, npy-1)
        q2(npx, npy) = q2(npx, npy-1)
        q2_tl(npx+1, npy) = q2_tl(npx, npy-2)
        q2(npx+1, npy) = q2(npx, npy-2)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q1_tl(0, 0) = q1_tl(1, 0)
        q1(0, 0) = q1(1, 0)
        q1_tl(0, -1) = q1_tl(2, 0)
        q1(0, -1) = q1(2, 0)
        q2_tl(0, 0) = q2_tl(1, 0)
        q2(0, 0) = q2(1, 0)
        q2_tl(0, -1) = q2_tl(2, 0)
        q2(0, -1) = q2(2, 0)
      END IF
      IF (se_corner) THEN
        q1_tl(npx, 0) = q1_tl(npx-1, 0)
        q1(npx, 0) = q1(npx-1, 0)
        q1_tl(npx, -1) = q1_tl(npx-2, 0)
        q1(npx, -1) = q1(npx-2, 0)
        q2_tl(npx, 0) = q2_tl(npx-1, 0)
        q2(npx, 0) = q2(npx-1, 0)
        q2_tl(npx, -1) = q2_tl(npx-2, 0)
        q2(npx, -1) = q2(npx-2, 0)
      END IF
      IF (nw_corner) THEN
        q1_tl(0, npy) = q1_tl(1, npy)
        q1(0, npy) = q1(1, npy)
        q1_tl(0, npy+1) = q1_tl(2, npy)
        q1(0, npy+1) = q1(2, npy)
        q2_tl(0, npy) = q2_tl(1, npy)
        q2(0, npy) = q2(1, npy)
        q2_tl(0, npy+1) = q2_tl(2, npy)
        q2(0, npy+1) = q2(2, npy)
      END IF
      IF (ne_corner) THEN
        q1_tl(npx, npy) = q1_tl(npx-1, npy)
        q1(npx, npy) = q1(npx-1, npy)
        q1_tl(npx, npy+1) = q1_tl(npx-2, npy)
        q1(npx, npy+1) = q1(npx-2, npy)
        q2_tl(npx, npy) = q2_tl(npx-1, npy)
        q2(npx, npy) = q2(npx-1, npy)
        q2_tl(npx, npy+1) = q2_tl(npx-2, npy)
        q2(npx, npy+1) = q2(npx-2, npy)
      END IF
    END SELECT
  END SUBROUTINE FILL2_4CORNERS_TLM
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
!  Differentiation of fill_4corners in forward (tangent) mode:
!   variations   of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE FILL_4CORNERS_TLM(q, q_tl, dir, bd, npx, npy, sw_corner, &
&   se_corner, ne_corner, nw_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
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
        q_tl(-1, 0) = q_tl(0, 2)
        q(-1, 0) = q(0, 2)
        q_tl(0, 0) = q_tl(0, 1)
        q(0, 0) = q(0, 1)
      END IF
      IF (se_corner) THEN
        q_tl(npx+1, 0) = q_tl(npx, 2)
        q(npx+1, 0) = q(npx, 2)
        q_tl(npx, 0) = q_tl(npx, 1)
        q(npx, 0) = q(npx, 1)
      END IF
      IF (nw_corner) THEN
        q_tl(0, npy) = q_tl(0, npy-1)
        q(0, npy) = q(0, npy-1)
        q_tl(-1, npy) = q_tl(0, npy-2)
        q(-1, npy) = q(0, npy-2)
      END IF
      IF (ne_corner) THEN
        q_tl(npx, npy) = q_tl(npx, npy-1)
        q(npx, npy) = q(npx, npy-1)
        q_tl(npx+1, npy) = q_tl(npx, npy-2)
        q(npx+1, npy) = q(npx, npy-2)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q_tl(0, 0) = q_tl(1, 0)
        q(0, 0) = q(1, 0)
        q_tl(0, -1) = q_tl(2, 0)
        q(0, -1) = q(2, 0)
      END IF
      IF (se_corner) THEN
        q_tl(npx, 0) = q_tl(npx-1, 0)
        q(npx, 0) = q(npx-1, 0)
        q_tl(npx, -1) = q_tl(npx-2, 0)
        q(npx, -1) = q(npx-2, 0)
      END IF
      IF (nw_corner) THEN
        q_tl(0, npy) = q_tl(1, npy)
        q(0, npy) = q(1, npy)
        q_tl(0, npy+1) = q_tl(2, npy)
        q(0, npy+1) = q(2, npy)
      END IF
      IF (ne_corner) THEN
        q_tl(npx, npy) = q_tl(npx-1, npy)
        q(npx, npy) = q(npx-1, npy)
        q_tl(npx, npy+1) = q_tl(npx-2, npy)
        q(npx, npy+1) = q(npx-2, npy)
      END IF
    END SELECT
  END SUBROUTINE FILL_4CORNERS_TLM
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
!  Differentiation of xtp_u in forward (tangent) mode:
!   variations   of useful results: flux
!   with respect to varying inputs: flux u c
  SUBROUTINE XTP_U_TLM(is, ie, js, je, isd, ied, jsd, jed, c, c_tl, u&
&   , u_tl, v, flux, flux_tl, iord, dx, rdx, npx, npy, grid_type, nested&
& )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: u_tl(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: c_tl(is:ie+1, js:je+1)
    REAL, INTENT(OUT) :: flux(is:ie+1, js:je+1)
    REAL, INTENT(OUT) :: flux_tl(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: rdx(isd:ied, jsd:jed+1)
    INTEGER, INTENT(IN) :: iord, npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
! Local
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    REAL, DIMENSION(is-1:ie+1) :: bl_tl, br_tl, b0_tl
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: al(is-1:ie+2), dm(is-2:ie+2)
    REAL :: al_tl(is-1:ie+2)
    REAL :: dq(is-3:ie+2)
    REAL :: dl, dr, xt, pmp, lac, cfl
    REAL :: xt_tl, cfl_tl
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: x0, x1, x0l, x0r
    INTEGER :: i, j
    INTEGER :: is3, ie3
    INTEGER :: is2, ie2
    INTRINSIC MAX
    INTRINSIC MIN
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
            flux_tl(i, j) = u_tl(i-1, j)
            flux(i, j) = u(i-1, j)
          ELSE
            flux_tl(i, j) = u_tl(i, j)
            flux(i, j) = u(i, j)
          END IF
        END DO
      END DO
    ELSE IF (iord .EQ. 333) THEN
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux_tl(i, j) = (2.0*u_tl(i, j)+5.0*u_tl(i-1, j)-u_tl(i-2, j&
&             ))/6.0 - 0.5*rdx(i-1, j)*(c_tl(i, j)*(u(i, j)-u(i-1, j))+c&
&             (i, j)*(u_tl(i, j)-u_tl(i-1, j))) + rdx(i-1, j)**2*(c_tl(i&
&             , j)*c(i, j)+c(i, j)*c_tl(i, j))*(u(i, j)-2.0*u(i-1, j)+u(&
&             i-2, j))/6.0 + c(i, j)**2*rdx(i-1, j)**2*(u_tl(i, j)-2.0*&
&             u_tl(i-1, j)+u_tl(i-2, j))/6.0
            flux(i, j) = (2.0*u(i, j)+5.0*u(i-1, j)-u(i-2, j))/6.0 - 0.5&
&             *c(i, j)*rdx(i-1, j)*(u(i, j)-u(i-1, j)) + c(i, j)*rdx(i-1&
&             , j)*c(i, j)*rdx(i-1, j)/6.0*(u(i, j)-2.0*u(i-1, j)+u(i-2&
&             , j))
          ELSE
            flux_tl(i, j) = (2.0*u_tl(i-1, j)+5.0*u_tl(i, j)-u_tl(i+1, j&
&             ))/6.0 - 0.5*rdx(i, j)*(c_tl(i, j)*(u(i, j)-u(i-1, j))+c(i&
&             , j)*(u_tl(i, j)-u_tl(i-1, j))) + rdx(i, j)**2*(c_tl(i, j)&
&             *c(i, j)+c(i, j)*c_tl(i, j))*(u(i+1, j)-2.0*u(i, j)+u(i-1&
&             , j))/6.0 + c(i, j)**2*rdx(i, j)**2*(u_tl(i+1, j)-2.0*u_tl&
&             (i, j)+u_tl(i-1, j))/6.0
            flux(i, j) = (2.0*u(i-1, j)+5.0*u(i, j)-u(i+1, j))/6.0 - 0.5&
&             *c(i, j)*rdx(i, j)*(u(i, j)-u(i-1, j)) + c(i, j)*rdx(i, j)&
&             *c(i, j)*rdx(i, j)/6.0*(u(i+1, j)-2.0*u(i, j)+u(i-1, j))
          END IF
        END DO
      END DO
    ELSE IF (iord .LT. 8) THEN
      al_tl = 0.0
      bl_tl = 0.0
      br_tl = 0.0
      b0_tl = 0.0
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
      DO j=js,je+1
        DO i=is3,ie3+1
          al_tl(i) = p1*(u_tl(i-1, j)+u_tl(i, j)) + p2*(u_tl(i-2, j)+&
&           u_tl(i+1, j))
          al(i) = p1*(u(i-1, j)+u(i, j)) + p2*(u(i-2, j)+u(i+1, j))
        END DO
        DO i=is3,ie3
          bl_tl(i) = al_tl(i) - u_tl(i, j)
          bl(i) = al(i) - u(i, j)
          br_tl(i) = al_tl(i+1) - u_tl(i, j)
          br(i) = al(i+1) - u(i, j)
        END DO
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            xt_tl = c3*u_tl(1, j) + c2*u_tl(2, j) + c1*u_tl(3, j)
            xt = c3*u(1, j) + c2*u(2, j) + c1*u(3, j)
            br_tl(1) = xt_tl - u_tl(1, j)
            br(1) = xt - u(1, j)
            bl_tl(2) = xt_tl - u_tl(2, j)
            bl(2) = xt - u(2, j)
            br_tl(2) = al_tl(3) - u_tl(2, j)
            br(2) = al(3) - u(2, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! out
              bl_tl(0) = 0.0
              bl(0) = 0.
! edge
              br_tl(0) = 0.0
              br(0) = 0.
! edge
              bl_tl(1) = 0.0
              bl(1) = 0.
! in
              br_tl(1) = 0.0
              br(1) = 0.
            ELSE
              bl_tl(0) = c1*u_tl(-2, j) + c2*u_tl(-1, j) + c3*u_tl(0, j)&
&               - u_tl(0, j)
              bl(0) = c1*u(-2, j) + c2*u(-1, j) + c3*u(0, j) - u(0, j)
              xt_tl = 0.5*(((2.*dx(0, j)+dx(-1, j))*u_tl(0, j)-dx(0, j)*&
&               u_tl(-1, j))/(dx(0, j)+dx(-1, j))+((2.*dx(1, j)+dx(2, j)&
&               )*u_tl(1, j)-dx(1, j)*u_tl(2, j))/(dx(1, j)+dx(2, j)))
              xt = 0.5*(((2.*dx(0, j)+dx(-1, j))*u(0, j)-dx(0, j)*u(-1, &
&               j))/(dx(0, j)+dx(-1, j))+((2.*dx(1, j)+dx(2, j))*u(1, j)&
&               -dx(1, j)*u(2, j))/(dx(1, j)+dx(2, j)))
              br_tl(0) = xt_tl - u_tl(0, j)
              br(0) = xt - u(0, j)
              bl_tl(1) = xt_tl - u_tl(1, j)
              bl(1) = xt - u(1, j)
            END IF
          END IF
!       call pert_ppm(1, u(2,j), bl(2), br(2), -1)
          IF (ie + 1 .EQ. npx) THEN
            bl_tl(npx-2) = al_tl(npx-2) - u_tl(npx-2, j)
            bl(npx-2) = al(npx-2) - u(npx-2, j)
            xt_tl = c1*u_tl(npx-3, j) + c2*u_tl(npx-2, j) + c3*u_tl(npx-&
&             1, j)
            xt = c1*u(npx-3, j) + c2*u(npx-2, j) + c3*u(npx-1, j)
            br_tl(npx-2) = xt_tl - u_tl(npx-2, j)
            br(npx-2) = xt - u(npx-2, j)
            bl_tl(npx-1) = xt_tl - u_tl(npx-1, j)
            bl(npx-1) = xt - u(npx-1, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! in
              bl_tl(npx-1) = 0.0
              bl(npx-1) = 0.
! edge
              br_tl(npx-1) = 0.0
              br(npx-1) = 0.
! edge
              bl_tl(npx) = 0.0
              bl(npx) = 0.
! out
              br_tl(npx) = 0.0
              br(npx) = 0.
            ELSE
              xt_tl = 0.5*(((2.*dx(npx-1, j)+dx(npx-2, j))*u_tl(npx-1, j&
&               )-dx(npx-1, j)*u_tl(npx-2, j))/(dx(npx-1, j)+dx(npx-2, j&
&               ))+((2.*dx(npx, j)+dx(npx+1, j))*u_tl(npx, j)-dx(npx, j)&
&               *u_tl(npx+1, j))/(dx(npx, j)+dx(npx+1, j)))
              xt = 0.5*(((2.*dx(npx-1, j)+dx(npx-2, j))*u(npx-1, j)-dx(&
&               npx-1, j)*u(npx-2, j))/(dx(npx-1, j)+dx(npx-2, j))+((2.*&
&               dx(npx, j)+dx(npx+1, j))*u(npx, j)-dx(npx, j)*u(npx+1, j&
&               ))/(dx(npx, j)+dx(npx+1, j)))
              br_tl(npx-1) = xt_tl - u_tl(npx-1, j)
              br(npx-1) = xt - u(npx-1, j)
              bl_tl(npx) = xt_tl - u_tl(npx, j)
              bl(npx) = xt - u(npx, j)
              br_tl(npx) = c3*u_tl(npx, j) + c2*u_tl(npx+1, j) + c1*u_tl&
&               (npx+2, j) - u_tl(npx, j)
              br(npx) = c3*u(npx, j) + c2*u(npx+1, j) + c1*u(npx+2, j) -&
&               u(npx, j)
            END IF
          END IF
        END IF
!       call pert_ppm(1, u(npx-2,j), bl(npx-2), br(npx-2), -1)
        DO i=is-1,ie+1
          b0_tl(i) = bl_tl(i) + br_tl(i)
          b0(i) = bl(i) + br(i)
        END DO
        IF (iord .EQ. 2) THEN
! Perfectly linear
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl_tl = rdx(i-1, j)*c_tl(i, j)
              cfl = c(i, j)*rdx(i-1, j)
              flux_tl(i, j) = u_tl(i-1, j) + (1.-cfl)*(br_tl(i-1)-cfl_tl&
&               *b0(i-1)-cfl*b0_tl(i-1)) - cfl_tl*(br(i-1)-cfl*b0(i-1))
              flux(i, j) = u(i-1, j) + (1.-cfl)*(br(i-1)-cfl*b0(i-1))
            ELSE
              cfl_tl = rdx(i, j)*c_tl(i, j)
              cfl = c(i, j)*rdx(i, j)
              flux_tl(i, j) = u_tl(i, j) + cfl_tl*(bl(i)+cfl*b0(i)) + (&
&               1.+cfl)*(bl_tl(i)+cfl_tl*b0(i)+cfl*b0_tl(i))
              flux(i, j) = u(i, j) + (1.+cfl)*(bl(i)+cfl*b0(i))
            END IF
          END DO
        END IF
      END DO
    END IF
  END SUBROUTINE XTP_U_TLM
!  Differentiation of ytp_v in forward (tangent) mode:
!   variations   of useful results: flux
!   with respect to varying inputs: v c
  SUBROUTINE YTP_V_TLM(is, ie, js, je, isd, ied, jsd, jed, c, c_tl, u&
&   , v, v_tl, flux, flux_tl, jord, dy, rdy, npx, npy, grid_type, nested&
& )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: v_tl(isd:ied+1, jsd:jed)
!  Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: c_tl(is:ie+1, js:je+1)
    REAL, INTENT(OUT) :: flux(is:ie+1, js:je+1)
    REAL, INTENT(OUT) :: flux_tl(is:ie+1, js:je+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rdy(isd:ied+1, jsd:jed)
    INTEGER, INTENT(IN) :: npx, npy, grid_type
    LOGICAL, INTENT(IN) :: nested
! Local:
    LOGICAL, DIMENSION(is:ie+1, js-1:je+1) :: smt5, smt6
    REAL :: fx0(is:ie+1)
    REAL :: dm(is:ie+1, js-2:je+2)
    REAL :: al(is:ie+1, js-1:je+2)
    REAL :: al_tl(is:ie+1, js-1:je+2)
    REAL, DIMENSION(is:ie+1, js-1:je+1) :: bl, br, b0
    REAL, DIMENSION(is:ie+1, js-1:je+1) :: bl_tl, br_tl, b0_tl
    REAL :: dq(is:ie+1, js-3:je+2)
    REAL :: xt, dl, dr, pmp, lac, cfl
    REAL :: xt_tl, cfl_tl
    REAL :: pmp_1, lac_1, pmp_2, lac_2
    REAL :: x0, x1, x0r, x0l
    INTEGER :: i, j, is1, ie1, js3, je3
    INTRINSIC MAX
    INTRINSIC MIN
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
      flux_tl = 0.0
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux_tl(i, j) = v_tl(i, j-1)
            flux(i, j) = v(i, j-1)
          ELSE
            flux_tl(i, j) = v_tl(i, j)
            flux(i, j) = v(i, j)
          END IF
        END DO
      END DO
    ELSE IF (jord .EQ. 333) THEN
      flux_tl = 0.0
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux_tl(i, j) = (2.0*v_tl(i, j)+5.0*v_tl(i, j-1)-v_tl(i, j-2&
&             ))/6.0 - 0.5*rdy(i, j-1)*(c_tl(i, j)*(v(i, j)-v(i, j-1))+c&
&             (i, j)*(v_tl(i, j)-v_tl(i, j-1))) + rdy(i, j-1)**2*(c_tl(i&
&             , j)*c(i, j)+c(i, j)*c_tl(i, j))*(v(i, j)-2.0*v(i, j-1)+v(&
&             i, j-2))/6.0 + c(i, j)**2*rdy(i, j-1)**2*(v_tl(i, j)-2.0*&
&             v_tl(i, j-1)+v_tl(i, j-2))/6.0
            flux(i, j) = (2.0*v(i, j)+5.0*v(i, j-1)-v(i, j-2))/6.0 - 0.5&
&             *c(i, j)*rdy(i, j-1)*(v(i, j)-v(i, j-1)) + c(i, j)*rdy(i, &
&             j-1)*c(i, j)*rdy(i, j-1)/6.0*(v(i, j)-2.0*v(i, j-1)+v(i, j&
&             -2))
          ELSE
            flux_tl(i, j) = (2.0*v_tl(i, j-1)+5.0*v_tl(i, j)-v_tl(i, j+1&
&             ))/6.0 - 0.5*rdy(i, j)*(c_tl(i, j)*(v(i, j)-v(i, j-1))+c(i&
&             , j)*(v_tl(i, j)-v_tl(i, j-1))) + rdy(i, j)**2*(c_tl(i, j)&
&             *c(i, j)+c(i, j)*c_tl(i, j))*(v(i, j+1)-2.0*v(i, j)+v(i, j&
&             -1))/6.0 + c(i, j)**2*rdy(i, j)**2*(v_tl(i, j+1)-2.0*v_tl(&
&             i, j)+v_tl(i, j-1))/6.0
            flux(i, j) = (2.0*v(i, j-1)+5.0*v(i, j)-v(i, j+1))/6.0 - 0.5&
&             *c(i, j)*rdy(i, j)*(v(i, j)-v(i, j-1)) + c(i, j)*rdy(i, j)&
&             *c(i, j)*rdy(i, j)/6.0*(v(i, j+1)-2.0*v(i, j)+v(i, j-1))
          END IF
        END DO
      END DO
    ELSE IF (jord .LT. 8) THEN
      al_tl = 0.0
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
      DO j=js3,je3+1
        DO i=is,ie+1
          al_tl(i, j) = p1*(v_tl(i, j-1)+v_tl(i, j)) + p2*(v_tl(i, j-2)+&
&           v_tl(i, j+1))
          al(i, j) = p1*(v(i, j-1)+v(i, j)) + p2*(v(i, j-2)+v(i, j+1))
        END DO
      END DO
      bl_tl = 0.0
      br_tl = 0.0
      DO j=js3,je3
        DO i=is,ie+1
          bl_tl(i, j) = al_tl(i, j) - v_tl(i, j)
          bl(i, j) = al(i, j) - v(i, j)
          br_tl(i, j) = al_tl(i, j+1) - v_tl(i, j)
          br(i, j) = al(i, j+1) - v(i, j)
        END DO
      END DO
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
            bl_tl(i, 0) = c1*v_tl(i, -2) + c2*v_tl(i, -1) + c3*v_tl(i, 0&
&             ) - v_tl(i, 0)
            bl(i, 0) = c1*v(i, -2) + c2*v(i, -1) + c3*v(i, 0) - v(i, 0)
            xt_tl = 0.5*(((2.*dy(i, 0)+dy(i, -1))*v_tl(i, 0)-dy(i, 0)*&
&             v_tl(i, -1))/(dy(i, 0)+dy(i, -1))+((2.*dy(i, 1)+dy(i, 2))*&
&             v_tl(i, 1)-dy(i, 1)*v_tl(i, 2))/(dy(i, 1)+dy(i, 2)))
            xt = 0.5*(((2.*dy(i, 0)+dy(i, -1))*v(i, 0)-dy(i, 0)*v(i, -1)&
&             )/(dy(i, 0)+dy(i, -1))+((2.*dy(i, 1)+dy(i, 2))*v(i, 1)-dy(&
&             i, 1)*v(i, 2))/(dy(i, 1)+dy(i, 2)))
            br_tl(i, 0) = xt_tl - v_tl(i, 0)
            br(i, 0) = xt - v(i, 0)
            bl_tl(i, 1) = xt_tl - v_tl(i, 1)
            bl(i, 1) = xt - v(i, 1)
            xt_tl = c3*v_tl(i, 1) + c2*v_tl(i, 2) + c1*v_tl(i, 3)
            xt = c3*v(i, 1) + c2*v(i, 2) + c1*v(i, 3)
            br_tl(i, 1) = xt_tl - v_tl(i, 1)
            br(i, 1) = xt - v(i, 1)
            bl_tl(i, 2) = xt_tl - v_tl(i, 2)
            bl(i, 2) = xt - v(i, 2)
            br_tl(i, 2) = al_tl(i, 3) - v_tl(i, 2)
            br(i, 2) = al(i, 3) - v(i, 2)
          END DO
          IF (is .EQ. 1) THEN
! out
            bl_tl(1, 0) = 0.0
            bl(1, 0) = 0.
! edge
            br_tl(1, 0) = 0.0
            br(1, 0) = 0.
! edge
            bl_tl(1, 1) = 0.0
            bl(1, 1) = 0.
! in
            br_tl(1, 1) = 0.0
            br(1, 1) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
! out
            bl_tl(npx, 0) = 0.0
            bl(npx, 0) = 0.
! edge
            br_tl(npx, 0) = 0.0
            br(npx, 0) = 0.
! edge
            bl_tl(npx, 1) = 0.0
            bl(npx, 1) = 0.
! in
            br_tl(npx, 1) = 0.0
            br(npx, 1) = 0.
          END IF
        END IF
!      j=2
!      call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
            bl_tl(i, npy-2) = al_tl(i, npy-2) - v_tl(i, npy-2)
            bl(i, npy-2) = al(i, npy-2) - v(i, npy-2)
            xt_tl = c1*v_tl(i, npy-3) + c2*v_tl(i, npy-2) + c3*v_tl(i, &
&             npy-1)
            xt = c1*v(i, npy-3) + c2*v(i, npy-2) + c3*v(i, npy-1)
            br_tl(i, npy-2) = xt_tl - v_tl(i, npy-2)
            br(i, npy-2) = xt - v(i, npy-2)
            bl_tl(i, npy-1) = xt_tl - v_tl(i, npy-1)
            bl(i, npy-1) = xt - v(i, npy-1)
            xt_tl = 0.5*(((2.*dy(i, npy-1)+dy(i, npy-2))*v_tl(i, npy-1)-&
&             dy(i, npy-1)*v_tl(i, npy-2))/(dy(i, npy-1)+dy(i, npy-2))+(&
&             (2.*dy(i, npy)+dy(i, npy+1))*v_tl(i, npy)-dy(i, npy)*v_tl(&
&             i, npy+1))/(dy(i, npy)+dy(i, npy+1)))
            xt = 0.5*(((2.*dy(i, npy-1)+dy(i, npy-2))*v(i, npy-1)-dy(i, &
&             npy-1)*v(i, npy-2))/(dy(i, npy-1)+dy(i, npy-2))+((2.*dy(i&
&             , npy)+dy(i, npy+1))*v(i, npy)-dy(i, npy)*v(i, npy+1))/(dy&
&             (i, npy)+dy(i, npy+1)))
            br_tl(i, npy-1) = xt_tl - v_tl(i, npy-1)
            br(i, npy-1) = xt - v(i, npy-1)
            bl_tl(i, npy) = xt_tl - v_tl(i, npy)
            bl(i, npy) = xt - v(i, npy)
            br_tl(i, npy) = c3*v_tl(i, npy) + c2*v_tl(i, npy+1) + c1*&
&             v_tl(i, npy+2) - v_tl(i, npy)
            br(i, npy) = c3*v(i, npy) + c2*v(i, npy+1) + c1*v(i, npy+2) &
&             - v(i, npy)
          END DO
          IF (is .EQ. 1) THEN
! in
            bl_tl(1, npy-1) = 0.0
            bl(1, npy-1) = 0.
! edge
            br_tl(1, npy-1) = 0.0
            br(1, npy-1) = 0.
! edge
            bl_tl(1, npy) = 0.0
            bl(1, npy) = 0.
! out
            br_tl(1, npy) = 0.0
            br(1, npy) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
! in
            bl_tl(npx, npy-1) = 0.0
            bl(npx, npy-1) = 0.
! edge
            br_tl(npx, npy-1) = 0.0
            br(npx, npy-1) = 0.
! edge
            bl_tl(npx, npy) = 0.0
            bl(npx, npy) = 0.
! out
            br_tl(npx, npy) = 0.0
            br(npx, npy) = 0.
            b0_tl = 0.0
          ELSE
            b0_tl = 0.0
          END IF
        ELSE
          b0_tl = 0.0
        END IF
      ELSE
        b0_tl = 0.0
      END IF
!      j=npy-2
!      call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
      DO j=js-1,je+1
        DO i=is,ie+1
          b0_tl(i, j) = bl_tl(i, j) + br_tl(i, j)
          b0(i, j) = bl(i, j) + br(i, j)
        END DO
      END DO
      IF (jord .EQ. 2) THEN
        flux_tl = 0.0
! Perfectly linear
        DO j=js,je+1
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              cfl_tl = rdy(i, j-1)*c_tl(i, j)
              cfl = c(i, j)*rdy(i, j-1)
              flux_tl(i, j) = v_tl(i, j-1) + (1.-cfl)*(br_tl(i, j-1)-&
&               cfl_tl*b0(i, j-1)-cfl*b0_tl(i, j-1)) - cfl_tl*(br(i, j-1&
&               )-cfl*b0(i, j-1))
              flux(i, j) = v(i, j-1) + (1.-cfl)*(br(i, j-1)-cfl*b0(i, j-&
&               1))
            ELSE
              cfl_tl = rdy(i, j)*c_tl(i, j)
              cfl = c(i, j)*rdy(i, j)
              flux_tl(i, j) = v_tl(i, j) + cfl_tl*(bl(i, j)+cfl*b0(i, j)&
&               ) + (1.+cfl)*(bl_tl(i, j)+cfl_tl*b0(i, j)+cfl*b0_tl(i, j&
&               ))
              flux(i, j) = v(i, j) + (1.+cfl)*(bl(i, j)+cfl*b0(i, j))
            END IF
          END DO
        END DO
      ELSE
        flux_tl = 0.0
      END IF
    ELSE
      flux_tl = 0.0
    END IF
  END SUBROUTINE YTP_V_TLM
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
    REAL :: arg1
    REAL :: result1
    REAL :: pwr1
    REAL*8 :: pwx1
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
            arg1 = delpc(i, j)**2 + vort(i, j)**2
            result1 = SQRT(arg1)
            vort(i, j) = abs0*result1
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
        pwr1 = d4_bg**n2
        dd8 = gridstruct%da_min*pwr1
      ELSE
        pwx1 = gridstruct%da_min_c*d4_bg
        dd8 = pwx1**n2
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
    REAL :: arg1
    REAL :: result1
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
        arg1 = sh(i, j)**2 + smag_c(i, j)**2
        result1 = SQRT(arg1)
        smag_c(i, j) = dt*result1
      END DO
    END DO
  END SUBROUTINE SMAG_CORNER
!  Differentiation of compute_divergence_damping in forward (tangent) mode:
!   variations   of useful results: ke uc ptc delpc vc vort divg_d
!                wk
!   with respect to varying inputs: u v ke ua uc ptc delpc va vc
!                divg_d wk
  SUBROUTINE COMPUTE_DIVERGENCE_DAMPING_TLM(nord, d2_bg, d4_bg, dddmp&
&   , dt, vort, vort_tl, ptc, ptc_tl, delpc, delpc_tl, ke, ke_tl, u, &
&   u_tl, v, v_tl, uc, uc_tl, vc, vc_tl, ua, ua_tl, va, va_tl, divg_d, &
&   divg_d_tl, wk, wk_tl, gridstruct, flagstruct, bd)
    IMPLICIT NONE
!InOut Arguments
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_FLAGS_TYPE), INTENT(IN), TARGET :: flagstruct
    INTEGER, INTENT(IN) :: nord
    REAL, INTENT(IN) :: d2_bg, d4_bg, dddmp, dt
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua_tl, &
&   va_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v_tl
!Intent is really in
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: wk
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   wk_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: vort
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   vort_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delpc, ptc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   delpc_tl, ptc_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   ke
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   ke_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: vc
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   vc_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: uc
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(INOUT) :: &
&   uc_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   divg_d
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), INTENT(INOUT) :: &
&   divg_d_tl
!Locals
    REAL :: damp, dd8, damp2, da_min, da_min_c, absdt
    REAL :: damp_tl, damp2_tl
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
    REAL :: max1_tl
    REAL :: abs0
    REAL :: abs1
    REAL :: max2
    REAL :: max2_tl
    REAL :: abs2
    REAL :: abs2_tl
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: result1
    REAL :: result1_tl
    REAL :: pwr1
    REAL*8 :: pwx1
    REAL :: y3_tl
    REAL :: y1_tl
    REAL :: y2_tl
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
            ptc_tl(i, j) = dyc(i, j)*sina_v(i, j)*(u_tl(i, j)-0.5*cosa_v&
&             (i, j)*(va_tl(i, j-1)+va_tl(i, j)))
            ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j))&
&             *dyc(i, j)*sina_v(i, j)
          END DO
        END DO
        vort_tl = 0.0
        DO j=js-1,je+1
          DO i=is2,ie1
            vort_tl(i, j) = dxc(i, j)*sina_u(i, j)*(v_tl(i, j)-0.5*&
&             cosa_u(i, j)*(ua_tl(i-1, j)+ua_tl(i, j)))
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
        END DO
      ELSE
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is-1,ie+1
              IF (vc(i, j) .GT. 0) THEN
                ptc_tl(i, j) = dyc(i, j)*sin_sg(i, j-1, 4)*u_tl(i, j)
                ptc(i, j) = u(i, j)*dyc(i, j)*sin_sg(i, j-1, 4)
              ELSE
                ptc_tl(i, j) = dyc(i, j)*sin_sg(i, j, 2)*u_tl(i, j)
                ptc(i, j) = u(i, j)*dyc(i, j)*sin_sg(i, j, 2)
              END IF
            END DO
          ELSE
            DO i=is-1,ie+1
              ptc_tl(i, j) = dyc(i, j)*sina_v(i, j)*(u_tl(i, j)-0.5*&
&               cosa_v(i, j)*(va_tl(i, j-1)+va_tl(i, j)))
              ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j&
&               ))*dyc(i, j)*sina_v(i, j)
            END DO
          END IF
        END DO
        vort_tl = 0.0
        DO j=js-1,je+1
          DO i=is2,ie1
            vort_tl(i, j) = dxc(i, j)*sina_u(i, j)*(v_tl(i, j)-0.5*&
&             cosa_u(i, j)*(ua_tl(i-1, j)+ua_tl(i, j)))
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
          IF (is .EQ. 1) THEN
            IF (uc(1, j) .GT. 0) THEN
              vort_tl(1, j) = dxc(1, j)*sin_sg(0, j, 3)*v_tl(1, j)
              vort(1, j) = v(1, j)*dxc(1, j)*sin_sg(0, j, 3)
            ELSE
              vort_tl(1, j) = dxc(1, j)*sin_sg(1, j, 1)*v_tl(1, j)
              vort(1, j) = v(1, j)*dxc(1, j)*sin_sg(1, j, 1)
            END IF
          END IF
          IF (ie + 1 .EQ. npx) THEN
            IF (uc(npx, j) .GT. 0) THEN
              vort_tl(npx, j) = dxc(npx, j)*sin_sg(npx-1, j, 3)*v_tl(npx&
&               , j)
              vort(npx, j) = v(npx, j)*dxc(npx, j)*sin_sg(npx-1, j, 3)
            ELSE
              vort_tl(npx, j) = dxc(npx, j)*sin_sg(npx, j, 1)*v_tl(npx, &
&               j)
              vort(npx, j) = v(npx, j)*dxc(npx, j)*sin_sg(npx, j, 1)
            END IF
          END IF
        END DO
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          delpc_tl(i, j) = vort_tl(i, j-1) - vort_tl(i, j) + ptc_tl(i-1&
&           , j) - ptc_tl(i, j)
          delpc(i, j) = vort(i, j-1) - vort(i, j) + (ptc(i-1, j)-ptc(i, &
&           j))
        END DO
      END DO
! Remove the extra term at the corners:
      IF (sw_corner) THEN
        delpc_tl(1, 1) = delpc_tl(1, 1) - vort_tl(1, 0)
        delpc(1, 1) = delpc(1, 1) - vort(1, 0)
      END IF
      IF (se_corner) THEN
        delpc_tl(npx, 1) = delpc_tl(npx, 1) - vort_tl(npx, 0)
        delpc(npx, 1) = delpc(npx, 1) - vort(npx, 0)
      END IF
      IF (ne_corner) THEN
        delpc_tl(npx, npy) = delpc_tl(npx, npy) + vort_tl(npx, npy)
        delpc(npx, npy) = delpc(npx, npy) + vort(npx, npy)
      END IF
      IF (nw_corner) THEN
        delpc_tl(1, npy) = delpc_tl(1, npy) + vort_tl(1, npy)
        delpc(1, npy) = delpc(1, npy) + vort(1, npy)
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          delpc_tl(i, j) = gridstruct%rarea_c(i, j)*delpc_tl(i, j)
          delpc(i, j) = gridstruct%rarea_c(i, j)*delpc(i, j)
          IF (delpc(i, j)*dt .GE. 0.) THEN
            abs2_tl = dt*delpc_tl(i, j)
            abs2 = delpc(i, j)*dt
          ELSE
            abs2_tl = -(dt*delpc_tl(i, j))
            abs2 = -(delpc(i, j)*dt)
          END IF
          y3_tl = dddmp*abs2_tl
          y3 = dddmp*abs2
          IF (0.20 .GT. y3) THEN
            y1_tl = y3_tl
            y1 = y3
          ELSE
            y1 = 0.20
            y1_tl = 0.0
          END IF
          IF (d2_bg .LT. y1) THEN
            max1_tl = y1_tl
            max1 = y1
          ELSE
            max1 = d2_bg
            max1_tl = 0.0
          END IF
          damp_tl = gridstruct%da_min_c*max1_tl
          damp = gridstruct%da_min_c*max1
          vort_tl(i, j) = damp_tl*delpc(i, j) + damp*delpc_tl(i, j)
          vort(i, j) = damp*delpc(i, j)
          ke_tl(i, j) = ke_tl(i, j) + vort_tl(i, j)
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
          delpc_tl(i, j) = divg_d_tl(i, j)
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
        IF (fill_c) CALL FILL_CORNERS_TLM(divg_d, divg_d_tl, npx, npy, &
&                                   fill=xdir, bgrid=.true.)
        DO j=js-nt,je+1+nt
          DO i=is-1-nt,ie+1+nt
            vc_tl(i, j) = divg_u(i, j)*(divg_d_tl(i+1, j)-divg_d_tl(i, j&
&             ))
            vc(i, j) = (divg_d(i+1, j)-divg_d(i, j))*divg_u(i, j)
          END DO
        END DO
        IF (fill_c) CALL FILL_CORNERS_TLM(divg_d, divg_d_tl, npx, npy, &
&                                   fill=ydir, bgrid=.true.)
        DO j=js-1-nt,je+1+nt
          DO i=is-nt,ie+1+nt
            uc_tl(i, j) = divg_v(i, j)*(divg_d_tl(i, j+1)-divg_d_tl(i, j&
&             ))
            uc(i, j) = (divg_d(i, j+1)-divg_d(i, j))*divg_v(i, j)
          END DO
        END DO
        IF (fill_c) CALL FILL_CORNERS_TLM(vc, vc_tl, uc, uc_tl, npx, npy&
&                                   , dgrid=.true., vector=.true.)
        DO j=js-nt,je+1+nt
          DO i=is-nt,ie+1+nt
            divg_d_tl(i, j) = uc_tl(i, j-1) - uc_tl(i, j) + vc_tl(i-1, j&
&             ) - vc_tl(i, j)
            divg_d(i, j) = uc(i, j-1) - uc(i, j) + (vc(i-1, j)-vc(i, j))
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) THEN
          divg_d_tl(1, 1) = divg_d_tl(1, 1) - uc_tl(1, 0)
          divg_d(1, 1) = divg_d(1, 1) - uc(1, 0)
        END IF
        IF (se_corner) THEN
          divg_d_tl(npx, 1) = divg_d_tl(npx, 1) - uc_tl(npx, 0)
          divg_d(npx, 1) = divg_d(npx, 1) - uc(npx, 0)
        END IF
        IF (ne_corner) THEN
          divg_d_tl(npx, npy) = divg_d_tl(npx, npy) + uc_tl(npx, npy)
          divg_d(npx, npy) = divg_d(npx, npy) + uc(npx, npy)
        END IF
        IF (nw_corner) THEN
          divg_d_tl(1, npy) = divg_d_tl(1, npy) + uc_tl(1, npy)
          divg_d(1, npy) = divg_d(1, npy) + uc(1, npy)
        END IF
        IF (.NOT.gridstruct%stretched_grid) THEN
          DO j=js-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              divg_d_tl(i, j) = gridstruct%rarea_c(i, j)*divg_d_tl(i, j)
              divg_d(i, j) = divg_d(i, j)*gridstruct%rarea_c(i, j)
            END DO
          END DO
        END IF
      END DO
! n-loop
      IF (dddmp .LT. 1.e-5) THEN
        vort(:, :) = 0.
        vort_tl = 0.0
      ELSE IF (flagstruct%grid_type .LT. 3) THEN
! Interpolate relative vort to cell corners
        vort_tl = 0.0
        CALL A2B_ORD4_TLM(wk, wk_tl, vort, vort_tl, gridstruct, npx, &
&                      npy, is, ie, js, je, ng, .false.)
        DO j=js,je+1
          DO i=is,ie+1
            IF (dt .GE. 0.) THEN
              abs0 = dt
            ELSE
              abs0 = -dt
            END IF
! The following is an approxi form of Smagorinsky diffusion
            arg1_tl = 2*delpc(i, j)*delpc_tl(i, j) + 2*vort(i, j)*&
&             vort_tl(i, j)
            arg1 = delpc(i, j)**2 + vort(i, j)**2
            IF (arg1 .EQ. 0.0) THEN
              result1_tl = 0.0
            ELSE
              result1_tl = arg1_tl/(2.0*SQRT(arg1))
            END IF
            result1 = SQRT(arg1)
            vort_tl(i, j) = abs0*result1_tl
            vort(i, j) = abs0*result1
          END DO
        END DO
      ELSE
        IF (dt .GE. 0.) THEN
          abs1 = dt
        ELSE
          abs1 = -dt
        END IF
! Correct form: works only for doubly preiodic domain
        CALL SMAG_CORNER_TLM(abs1, u, u_tl, v, v_tl, ua, va, vort, &
&                         vort_tl, bd, npx, npy, gridstruct, ng)
      END IF
      IF (gridstruct%stretched_grid) THEN
! Stretched grid with variable damping ~ area
        pwr1 = d4_bg**n2
        dd8 = gridstruct%da_min*pwr1
      ELSE
        pwx1 = gridstruct%da_min_c*d4_bg
        dd8 = pwx1**n2
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          IF (0.20 .GT. dddmp*vort(i, j)) THEN
            y2_tl = dddmp*vort_tl(i, j)
            y2 = dddmp*vort(i, j)
          ELSE
            y2 = 0.20
            y2_tl = 0.0
          END IF
          IF (d2_bg .LT. y2) THEN
            max2_tl = y2_tl
            max2 = y2
          ELSE
            max2 = d2_bg
            max2_tl = 0.0
          END IF
! del-2
          damp2_tl = gridstruct%da_min_c*max2_tl
          damp2 = gridstruct%da_min_c*max2
          vort_tl(i, j) = damp2_tl*delpc(i, j) + damp2*delpc_tl(i, j) + &
&           dd8*divg_d_tl(i, j)
          vort(i, j) = damp2*delpc(i, j) + dd8*divg_d(i, j)
          ke_tl(i, j) = ke_tl(i, j) + vort_tl(i, j)
          ke(i, j) = ke(i, j) + vort(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE COMPUTE_DIVERGENCE_DAMPING_TLM
!  Differentiation of smag_corner in forward (tangent) mode:
!   variations   of useful results: smag_c
!   with respect to varying inputs: u v
  SUBROUTINE SMAG_CORNER_TLM(dt, u, u_tl, v, v_tl, ua, va, smag_c, &
&   smag_c_tl, bd, npx, npy, gridstruct, ng)
    IMPLICIT NONE
! Compute the Tension_Shear strain at cell corners for Smagorinsky diffusion
!!!  work only if (grid_type==4)
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    REAL, INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: npx, npy, ng
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed+1), INTENT(IN) :: u_tl
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v
    REAL, DIMENSION(bd%isd:bd%ied+1, bd%jsd:bd%jed), INTENT(IN) :: v_tl
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(IN) :: ua, va
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: smag_c
    REAL, DIMENSION(bd%isd:bd%ied, bd%jsd:bd%jed), INTENT(OUT) :: &
&   smag_c_tl
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! local
    REAL :: ut(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: ut_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed)
    REAL :: vt(bd%isd:bd%ied, bd%jsd:bd%jed+1)
    REAL :: vt_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1)
!  work array
    REAL :: wk(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: wk_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: sh(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: sh_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    INTEGER :: i, j
    INTEGER :: is2, ie1
    REAL, DIMENSION(:, :), POINTER :: dxc, dyc, dx, dy, rarea, rarea_c
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SQRT
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: result1
    REAL :: result1_tl
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
    ut_tl = 0.0
! Smag = sqrt [ T**2 + S**2 ]:  unit = 1/s
! where T = du/dx - dv/dy;   S = du/dy + dv/dx
! Compute tension strain at corners:
    DO j=js,je+1
      DO i=is-1,ie+1
        ut_tl(i, j) = dyc(i, j)*u_tl(i, j)
        ut(i, j) = u(i, j)*dyc(i, j)
      END DO
    END DO
    vt_tl = 0.0
    DO j=js-1,je+1
      DO i=is,ie+1
        vt_tl(i, j) = dxc(i, j)*v_tl(i, j)
        vt(i, j) = v(i, j)*dxc(i, j)
      END DO
    END DO
    smag_c_tl = 0.0
    DO j=js,je+1
      DO i=is,ie+1
        smag_c_tl(i, j) = rarea_c(i, j)*(vt_tl(i, j-1)-vt_tl(i, j)-ut_tl&
&         (i-1, j)+ut_tl(i, j))
        smag_c(i, j) = rarea_c(i, j)*(vt(i, j-1)-vt(i, j)-ut(i-1, j)+ut(&
&         i, j))
      END DO
    END DO
! Fix the corners?? if grid_type /= 4
! Compute shear strain:
    DO j=jsd,jed+1
      DO i=isd,ied
        vt_tl(i, j) = dx(i, j)*u_tl(i, j)
        vt(i, j) = u(i, j)*dx(i, j)
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied+1
        ut_tl(i, j) = dy(i, j)*v_tl(i, j)
        ut(i, j) = v(i, j)*dy(i, j)
      END DO
    END DO
    wk_tl = 0.0
    DO j=jsd,jed
      DO i=isd,ied
        wk_tl(i, j) = rarea(i, j)*(vt_tl(i, j)-vt_tl(i, j+1)+ut_tl(i, j)&
&         -ut_tl(i+1, j))
        wk(i, j) = rarea(i, j)*(vt(i, j)-vt(i, j+1)+ut(i, j)-ut(i+1, j))
      END DO
    END DO
    sh_tl = 0.0
    CALL A2B_ORD4_TLM(wk, wk_tl, sh, sh_tl, gridstruct, npx, npy, is&
&                  , ie, js, je, ng, .false.)
    DO j=js,je+1
      DO i=is,ie+1
        arg1_tl = 2*sh(i, j)*sh_tl(i, j) + 2*smag_c(i, j)*smag_c_tl(i, j&
&         )
        arg1 = sh(i, j)**2 + smag_c(i, j)**2
        IF (arg1 .EQ. 0.0) THEN
          result1_tl = 0.0
        ELSE
          result1_tl = arg1_tl/(2.0*SQRT(arg1))
        END IF
        result1 = SQRT(arg1)
        smag_c_tl(i, j) = dt*result1_tl
        smag_c(i, j) = dt*result1
      END DO
    END DO
  END SUBROUTINE SMAG_CORNER_TLM
 end module sw_core_tlm_mod
