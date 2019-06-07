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
 module fv_grid_utils_tlm_mod
 
#include <fms_platform.h>
 use constants_mod,   only: omega, pi=>pi_8, cnst_radius=>radius
 use mpp_mod,         only: FATAL, mpp_error, WARNING
 use external_sst_mod, only: i_sst, j_sst, sst_ncep, sst_anom
 use mpp_domains_mod, only: mpp_update_domains, DGRID_NE, mpp_global_sum
 use fv_mp_tlm_mod,   only: mpp_update_domains_tlm
 use mpp_domains_mod, only: BITWISE_EXACT_SUM, domain2d, BITWISE_EFP_SUM
 use mpp_parameter_mod, only: AGRID_PARAM=>AGRID, CGRID_NE_PARAM=>CGRID_NE
 use mpp_parameter_mod, only: CORNER, SCALAR_PAIR

 use fv_arrays_mod,   only: fv_atmos_type, fv_grid_type, fv_grid_bounds_type, &
                            R_GRID
 use fv_eta_mod,      only: set_eta
 use fv_mp_mod,       only: ng, is_master
 use fv_mp_mod,       only: mp_reduce_sum, mp_reduce_min, mp_reduce_max
 use fv_mp_mod,       only: fill_corners, XDir, YDir
 use fv_timing_mod,   only: timing_on, timing_off

 implicit none
 private
 logical:: symm_grid
#ifdef NO_QUAD_PRECISION
! 64-bit precision (kind=8)
 integer, parameter:: f_p = selected_real_kind(15)
#else
! Higher precision (kind=16) for grid geometrical factors:
 integer, parameter:: f_p = selected_real_kind(20)
#endif
 real, parameter::  big_number=1.d8
 real, parameter:: tiny_number=1.d-8

 real(kind=R_GRID) :: radius=cnst_radius

 real, parameter:: ptop_min=1.d-8

 public f_p 
 public ptop_min, big_number !CLEANUP: OK to keep since they are constants?
! public cos_angle
! public latlon2xyz, gnomonic_grids, &
!        global_mx, unit_vect_latlon,  &
!        cubed_to_latlon, c2l_ord2, g_sum, global_qsum, great_circle_dist,  &
!        v_prod, get_unit_vect2, project_sphere_v
! public mid_pt_sphere,  mid_pt_cart, vect_cross, grid_utils_init, grid_utils_end, &
!        spherical_angle, cell_center2, get_area, inner_prod, fill_ghost, direct_transform,  &
!        make_eta_level, expand_cell, cart_to_latlon, intp_great_circle, normalize_vect, &
!        dist2side_latlon, spherical_linear_interpolation, get_latlon_vector
! public symm_grid

public cubed_to_latlon, c2l_ord2, g_sum, great_circle_dist
public c2l_ord2_tlm, g_sum_tlm

! INTERFACE fill_ghost
!#ifdef OVERLOAD_R4
!   MODULE PROCEDURE fill_ghost_r4
!#endif
!   MODULE PROCEDURE fill_ghost_r8
! END INTERFACE

!---- version number -----
 character(len=128) :: version = '$Id$'
 character(len=128) :: tagname = '$Name$'

CONTAINS
  SUBROUTINE GRID_UTILS_INIT(atm, npx, npy, npz, non_ortho, grid_type, &
&   c2l_order)
    IMPLICIT NONE
! Initialize 2D memory and geometrical factors
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT), TARGET :: atm
    LOGICAL, INTENT(IN) :: non_ortho
    INTEGER, INTENT(IN) :: npx, npy, npz
    INTEGER, INTENT(IN) :: grid_type, c2l_order
  END SUBROUTINE GRID_UTILS_INIT
  SUBROUTINE GRID_UTILS_END()
    IMPLICIT NONE
  END SUBROUTINE GRID_UTILS_END
  REAL FUNCTION GREAT_CIRCLE_DIST(q1, q2, radius)
    IMPLICIT NONE
    REAL(kind=r_grid), INTENT(IN) :: q1(2), q2(2)
    REAL(kind=r_grid), INTENT(IN), OPTIONAL :: radius
    REAL(f_p) :: p1(2), p2(2)
    REAL(f_p) :: beta
    INTEGER :: n
    INTRINSIC SIN
    INTRINSIC COS
    INTRINSIC SQRT
    INTRINSIC ASIN
    INTRINSIC PRESENT
    REAL(f_p) :: arg1
    REAL(f_p) :: arg2
    REAL(f_p) :: arg3
    REAL(f_p) :: result1
    REAL(f_p) :: result2
    DO n=1,2
      p1(n) = q1(n)
      p2(n) = q2(n)
    END DO
    arg1 = (p1(2)-p2(2))/2.
    arg2 = (p1(1)-p2(1))/2.
    arg3 = SIN(arg1)**2 + COS(p1(2))*COS(p2(2))*SIN(arg2)**2
    result1 = SQRT(arg3)
    result2 = ASIN(result1)
    beta = result2*2.
    IF (PRESENT(radius)) THEN
      great_circle_dist = radius*beta
    ELSE
! Returns the angle
      great_circle_dist = beta
    END IF
  END FUNCTION GREAT_CIRCLE_DIST
  SUBROUTINE CUBED_TO_LATLON(u, v, ua, va, gridstruct, npx, npy, km, &
&   mode, grid_type, domain, nested, c2l_ord, bd)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: km, npx, npy, grid_type, c2l_ord
! update if present
    INTEGER, INTENT(IN) :: mode
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    LOGICAL, INTENT(IN) :: nested
    IF (c2l_ord .EQ. 2) THEN
      CALL C2L_ORD2(u, v, ua, va, gridstruct, km, grid_type, bd, .false.&
&            )
    ELSE
      CALL C2L_ORD4(u, v, ua, va, gridstruct, npx, npy, km, grid_type, &
&             domain, nested, mode, bd)
    END IF
  END SUBROUTINE CUBED_TO_LATLON
  SUBROUTINE C2L_ORD4(u, v, ua, va, gridstruct, npx, npy, km, grid_type&
&   , domain, nested, mode, bd)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: km, npx, npy, grid_type
! update if present
    INTEGER, INTENT(IN) :: mode
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    LOGICAL, INTENT(IN) :: nested
! Local
! 4-pt Lagrange interpolation
    REAL, SAVE :: a1=0.5625
    REAL, SAVE :: a2=-0.0625
    REAL, SAVE :: c1=1.125
    REAL, SAVE :: c2=-0.125
    REAL :: utmp(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: vtmp(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: wu(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: wv(bd%is:bd%ie+1, bd%js:bd%je)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    IF (mode .GT. 0) THEN
      CALL TIMING_ON('COMM_TOTAL')
      CALL MPP_UPDATE_DOMAINS(u, v, domain, gridtype=dgrid_ne)
      CALL TIMING_OFF('COMM_TOTAL')
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,km,npx,npy,grid_type,nested,c2,c1, &
!$OMP                                  u,v,gridstruct,ua,va,a1,a2)         &
!$OMP                          private(utmp, vtmp, wu, wv)
    DO k=1,km
      IF (grid_type .LT. 4) THEN
!nested
        IF (nested) THEN
          IF (1 .LT. js) THEN
            max1 = js
          ELSE
            max1 = 1
          END IF
          IF (npy - 1 .GT. je) THEN
            min1 = je
          ELSE
            min1 = npy - 1
          END IF
          DO j=max1,min1
            IF (1 .LT. is) THEN
              max2 = is
            ELSE
              max2 = 1
            END IF
            IF (npx - 1 .GT. ie) THEN
              min2 = ie
            ELSE
              min2 = npx - 1
            END IF
            DO i=max2,min2
              utmp(i, j) = c2*(u(i, j-1, k)+u(i, j+2, k)) + c1*(u(i, j, &
&               k)+u(i, j+1, k))
              vtmp(i, j) = c2*(v(i-1, j, k)+v(i+2, j, k)) + c1*(v(i, j, &
&               k)+v(i+1, j, k))
            END DO
          END DO
        ELSE
          IF (2 .LT. js) THEN
            max3 = js
          ELSE
            max3 = 2
          END IF
          IF (npy - 2 .GT. je) THEN
            min3 = je
          ELSE
            min3 = npy - 2
          END IF
          DO j=max3,min3
            IF (2 .LT. is) THEN
              max4 = is
            ELSE
              max4 = 2
            END IF
            IF (npx - 2 .GT. ie) THEN
              min4 = ie
            ELSE
              min4 = npx - 2
            END IF
            DO i=max4,min4
              utmp(i, j) = c2*(u(i, j-1, k)+u(i, j+2, k)) + c1*(u(i, j, &
&               k)+u(i, j+1, k))
              vtmp(i, j) = c2*(v(i-1, j, k)+v(i+2, j, k)) + c1*(v(i, j, &
&               k)+v(i+1, j, k))
            END DO
          END DO
          IF (js .EQ. 1) THEN
            DO i=is,ie+1
              wv(i, 1) = v(i, 1, k)*gridstruct%dy(i, 1)
            END DO
            DO i=is,ie
              vtmp(i, 1) = 2.*(wv(i, 1)+wv(i+1, 1))/(gridstruct%dy(i, 1)&
&               +gridstruct%dy(i+1, 1))
              utmp(i, 1) = 2.*(u(i, 1, k)*gridstruct%dx(i, 1)+u(i, 2, k)&
&               *gridstruct%dx(i, 2))/(gridstruct%dx(i, 1)+gridstruct%dx&
&               (i, 2))
            END DO
          END IF
!!!         vtmp(i,1) = (wv(i,1) + wv(i+1,1)) * gridstruct%rdya(i,1)
!!!         utmp(i,1) = (u(i,1,k)*gridstruct%dx(i,1) + u(i,2,k)*gridstruct%dx(i,2)) * gridstruct%rdxa(i,1)
          IF (je + 1 .EQ. npy) THEN
            j = npy - 1
            DO i=is,ie+1
              wv(i, j) = v(i, j, k)*gridstruct%dy(i, j)
            END DO
            DO i=is,ie
              vtmp(i, j) = 2.*(wv(i, j)+wv(i+1, j))/(gridstruct%dy(i, j)&
&               +gridstruct%dy(i+1, j))
              utmp(i, j) = 2.*(u(i, j, k)*gridstruct%dx(i, j)+u(i, j+1, &
&               k)*gridstruct%dx(i, j+1))/(gridstruct%dx(i, j)+&
&               gridstruct%dx(i, j+1))
            END DO
          END IF
!!!         vtmp(i,j) = (wv(i,j) + wv(i+1,j)) * gridstruct%rdya(i,j)
!!!         utmp(i,j) = (u(i,j,k)*gridstruct%dx(i,j) + u(i,j+1,k)*gridstruct%dx(i,j+1)) * gridstruct%rdxa(i,j)
          IF (is .EQ. 1) THEN
            i = 1
            DO j=js,je
              wv(1, j) = v(1, j, k)*gridstruct%dy(1, j)
              wv(2, j) = v(2, j, k)*gridstruct%dy(2, j)
            END DO
            DO j=js,je+1
              wu(i, j) = u(i, j, k)*gridstruct%dx(i, j)
            END DO
            DO j=js,je
              utmp(i, j) = 2.*(wu(i, j)+wu(i, j+1))/(gridstruct%dx(i, j)&
&               +gridstruct%dx(i, j+1))
              vtmp(i, j) = 2.*(wv(1, j)+wv(2, j))/(gridstruct%dy(1, j)+&
&               gridstruct%dy(2, j))
            END DO
          END IF
!!!      utmp(i,j) = (wu(i,j) + wu(i,  j+1)) * gridstruct%rdxa(i,j)
!!!      vtmp(i,j) = (wv(i,j) + wv(i+1,j  )) * gridstruct%rdya(i,j)
          IF (ie + 1 .EQ. npx) THEN
            i = npx - 1
            DO j=js,je
              wv(i, j) = v(i, j, k)*gridstruct%dy(i, j)
              wv(i+1, j) = v(i+1, j, k)*gridstruct%dy(i+1, j)
            END DO
            DO j=js,je+1
              wu(i, j) = u(i, j, k)*gridstruct%dx(i, j)
            END DO
            DO j=js,je
              utmp(i, j) = 2.*(wu(i, j)+wu(i, j+1))/(gridstruct%dx(i, j)&
&               +gridstruct%dx(i, j+1))
              vtmp(i, j) = 2.*(wv(i, j)+wv(i+1, j))/(gridstruct%dy(i, j)&
&               +gridstruct%dy(i+1, j))
            END DO
          END IF
        END IF
!!!      utmp(i,j) = (wu(i,j) + wu(i,  j+1)) * gridstruct%rdxa(i,j)
!!!      vtmp(i,j) = (wv(i,j) + wv(i+1,j  )) * gridstruct%rdya(i,j)
!Transform local a-grid winds into latitude-longitude coordinates
        DO j=js,je
          DO i=is,ie
            ua(i, j, k) = gridstruct%a11(i, j)*utmp(i, j) + gridstruct%&
&             a12(i, j)*vtmp(i, j)
            va(i, j, k) = gridstruct%a21(i, j)*utmp(i, j) + gridstruct%&
&             a22(i, j)*vtmp(i, j)
          END DO
        END DO
      ELSE
! Simple Cartesian Geometry:
        DO j=js,je
          DO i=is,ie
            ua(i, j, k) = a2*(u(i, j-1, k)+u(i, j+2, k)) + a1*(u(i, j, k&
&             )+u(i, j+1, k))
            va(i, j, k) = a2*(v(i-1, j, k)+v(i+2, j, k)) + a1*(v(i, j, k&
&             )+v(i+1, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE C2L_ORD4
!  Differentiation of c2l_ord2 in forward (tangent) mode:
!   variations   of useful results: ua va
!   with respect to varying inputs: u v ua va
  SUBROUTINE C2L_ORD2_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, &
&   gridstruct, km, grid_type, bd, do_halo)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: km, grid_type
    REAL, INTENT(IN) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(IN) :: u_tl(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(IN) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    REAL, INTENT(IN) :: v_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, INTENT(IN) :: do_halo
!
    REAL, INTENT(OUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: ua_tl(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: va_tl(bd%isd:bd%ied, bd%jsd:bd%jed, km)
!--------------------------------------------------------------
! Local
    REAL :: wu(bd%is-1:bd%ie+1, bd%js-1:bd%je+2)
    REAL :: wu_tl(bd%is-1:bd%ie+1, bd%js-1:bd%je+2)
    REAL :: wv(bd%is-1:bd%ie+2, bd%js-1:bd%je+1)
    REAL :: wv_tl(bd%is-1:bd%ie+2, bd%js-1:bd%je+1)
    REAL :: u1(bd%is-1:bd%ie+1), v1(bd%is-1:bd%ie+1)
    REAL :: u1_tl(bd%is-1:bd%ie+1), v1_tl(bd%is-1:bd%ie+1)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    REAL, DIMENSION(:, :), POINTER :: a11, a12, a21, a22
    REAL, DIMENSION(:, :), POINTER :: dx, dy, rdxa, rdya
    a11 => gridstruct%a11
    a12 => gridstruct%a12
    a21 => gridstruct%a21
    a22 => gridstruct%a22
    dx => gridstruct%dx
    dy => gridstruct%dy
    rdxa => gridstruct%rdxa
    rdya => gridstruct%rdya
    IF (do_halo) THEN
      is = bd%is - 1
      ie = bd%ie + 1
      js = bd%js - 1
      je = bd%je + 1
      wu_tl = 0.0
      v1_tl = 0.0
      wv_tl = 0.0
      u1_tl = 0.0
    ELSE
      is = bd%is
      ie = bd%ie
      js = bd%js
      je = bd%je
      wu_tl = 0.0
      v1_tl = 0.0
      wv_tl = 0.0
      u1_tl = 0.0
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,km,grid_type,u,dx,v,dy,ua,va,a11,a12,a21,a22) &
!$OMP                          private(u1, v1, wu, wv)
    DO k=1,km
      IF (grid_type .LT. 4) THEN
        DO j=js,je+1
          DO i=is,ie
            wu_tl(i, j) = dx(i, j)*u_tl(i, j, k)
            wu(i, j) = u(i, j, k)*dx(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            wv_tl(i, j) = dy(i, j)*v_tl(i, j, k)
            wv(i, j) = v(i, j, k)*dy(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
            u1_tl(i) = 2.*(wu_tl(i, j)+wu_tl(i, j+1))/(dx(i, j)+dx(i, j+&
&             1))
            u1(i) = 2.*(wu(i, j)+wu(i, j+1))/(dx(i, j)+dx(i, j+1))
            v1_tl(i) = 2.*(wv_tl(i, j)+wv_tl(i+1, j))/(dy(i, j)+dy(i+1, &
&             j))
            v1(i) = 2.*(wv(i, j)+wv(i+1, j))/(dy(i, j)+dy(i+1, j))
!!!          u1(i) = (wu(i,j) + wu(i,j+1)) * rdxa(i,j)
!!!          v1(i) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
! Cubed (cell center co-variant winds) to lat-lon:
            ua_tl(i, j, k) = a11(i, j)*u1_tl(i) + a12(i, j)*v1_tl(i)
            ua(i, j, k) = a11(i, j)*u1(i) + a12(i, j)*v1(i)
            va_tl(i, j, k) = a21(i, j)*u1_tl(i) + a22(i, j)*v1_tl(i)
            va(i, j, k) = a21(i, j)*u1(i) + a22(i, j)*v1(i)
          END DO
        END DO
      ELSE
! 2nd order:
        DO j=js,je
          DO i=is,ie
            ua_tl(i, j, k) = 0.5*(u_tl(i, j, k)+u_tl(i, j+1, k))
            ua(i, j, k) = 0.5*(u(i, j, k)+u(i, j+1, k))
            va_tl(i, j, k) = 0.5*(v_tl(i, j, k)+v_tl(i+1, j, k))
            va(i, j, k) = 0.5*(v(i, j, k)+v(i+1, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE C2L_ORD2_TLM
  SUBROUTINE C2L_ORD2(u, v, ua, va, gridstruct, km, grid_type, bd, &
&   do_halo)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: km, grid_type
    REAL, INTENT(IN) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(IN) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, INTENT(IN) :: do_halo
!
    REAL, INTENT(OUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
!--------------------------------------------------------------
! Local
    REAL :: wu(bd%is-1:bd%ie+1, bd%js-1:bd%je+2)
    REAL :: wv(bd%is-1:bd%ie+2, bd%js-1:bd%je+1)
    REAL :: u1(bd%is-1:bd%ie+1), v1(bd%is-1:bd%ie+1)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    REAL, DIMENSION(:, :), POINTER :: a11, a12, a21, a22
    REAL, DIMENSION(:, :), POINTER :: dx, dy, rdxa, rdya
    a11 => gridstruct%a11
    a12 => gridstruct%a12
    a21 => gridstruct%a21
    a22 => gridstruct%a22
    dx => gridstruct%dx
    dy => gridstruct%dy
    rdxa => gridstruct%rdxa
    rdya => gridstruct%rdya
    IF (do_halo) THEN
      is = bd%is - 1
      ie = bd%ie + 1
      js = bd%js - 1
      je = bd%je + 1
    ELSE
      is = bd%is
      ie = bd%ie
      js = bd%js
      je = bd%je
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,km,grid_type,u,dx,v,dy,ua,va,a11,a12,a21,a22) &
!$OMP                          private(u1, v1, wu, wv)
    DO k=1,km
      IF (grid_type .LT. 4) THEN
        DO j=js,je+1
          DO i=is,ie
            wu(i, j) = u(i, j, k)*dx(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            wv(i, j) = v(i, j, k)*dy(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
            u1(i) = 2.*(wu(i, j)+wu(i, j+1))/(dx(i, j)+dx(i, j+1))
            v1(i) = 2.*(wv(i, j)+wv(i+1, j))/(dy(i, j)+dy(i+1, j))
!!!          u1(i) = (wu(i,j) + wu(i,j+1)) * rdxa(i,j)
!!!          v1(i) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
! Cubed (cell center co-variant winds) to lat-lon:
            ua(i, j, k) = a11(i, j)*u1(i) + a12(i, j)*v1(i)
            va(i, j, k) = a21(i, j)*u1(i) + a22(i, j)*v1(i)
          END DO
        END DO
      ELSE
! 2nd order:
        DO j=js,je
          DO i=is,ie
            ua(i, j, k) = 0.5*(u(i, j, k)+u(i, j+1, k))
            va(i, j, k) = 0.5*(v(i, j, k)+v(i+1, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE C2L_ORD2
!  Differentiation of g_sum in forward (tangent) mode:
!   variations   of useful results: g_sum
!   with respect to varying inputs: p
 REAL FUNCTION G_SUM_TLM(domain, p, p_tl, ifirst, ilast, jfirst, jlast, ngc, area, mode, reproduce, g_sum)
! Fast version of globalsum 
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast, ngc
      integer, intent(IN) :: mode  ! if ==1 divided by area
      logical, intent(in), optional :: reproduce
      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
      real, intent(IN) :: p_tl(ifirst:ilast,jfirst:jlast)      ! field to be summed
      real(kind=R_GRID), intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
      type(domain2d), intent(IN) :: domain
      integer :: i,j
      real gsum, gsum_tl
      logical, SAVE :: g_sum_initialized = .false.
      real(kind=R_GRID), SAVE :: global_area
      real :: tmp(ifirst:ilast,jfirst:jlast) 
      real :: g_sum        

      if ( .not. g_sum_initialized ) then
         global_area = mpp_global_sum(domain, area, flags=BITWISE_EFP_SUM)
         if ( is_master() ) write(*,*) 'Global Area=',global_area
         g_sum_initialized = .true.
      end if
 
!-------------------------
! FMS global sum algorithm:
!-------------------------
      if ( present(reproduce) ) then
         if (reproduce) then
            gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast), &
                                  flags=BITWISE_EFP_SUM)
            gsum_tl = mpp_global_sum(domain, p_tl(:,:)*area(ifirst:ilast,jfirst:jlast), &
                                  flags=BITWISE_EFP_SUM)
         else
            gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast))
            gsum_tl = mpp_global_sum(domain, p_tl(:,:)*area(ifirst:ilast,jfirst:jlast))
         endif
      else
!-------------------------
! Quick local sum algorithm
!-------------------------
         gsum = 0.
         gsum_tl = 0.
         do j=jfirst,jlast
            do i=ifirst,ilast
               gsum = gsum + p(i,j)*area(i,j)
               gsum_tl = gsum_tl + p_tl(i,j)*area(i,j)
            enddo
         enddo
         call mp_reduce_sum(gsum)
         call mp_reduce_sum(gsum_tl)
      endif

      if ( mode==1 ) then
           g_sum = gsum / global_area
           g_sum_tlm = gsum_tl / global_area
      else
           g_sum = gsum
           g_sum_tlm = gsum_tl
      endif

      !Uncomment here to test without TLM/ADM
      g_sum_tlm = 0.0

 END FUNCTION G_SUM_TLM
 REAL FUNCTION G_SUM(domain, p, ifirst, ilast, jfirst, jlast, ngc, area, mode, reproduce)
! Fast version of globalsum 
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast, ngc
      integer, intent(IN) :: mode  ! if ==1 divided by area
      logical, intent(in), optional :: reproduce
      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
      real(kind=R_GRID), intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
      type(domain2d), intent(IN) :: domain
      integer :: i,j
      real gsum
      logical, SAVE :: g_sum_initialized = .false.
      real(kind=R_GRID), SAVE :: global_area
      real :: tmp(ifirst:ilast,jfirst:jlast) 
        
      if ( .not. g_sum_initialized ) then
         global_area = mpp_global_sum(domain, area, flags=BITWISE_EFP_SUM)
         if ( is_master() ) write(*,*) 'Global Area=',global_area
         g_sum_initialized = .true.
      end if
 
!-------------------------
! FMS global sum algorithm:
!-------------------------
      if ( present(reproduce) ) then
         if (reproduce) then
            gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast), &
                                  flags=BITWISE_EFP_SUM)
         else
            gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast))
         endif
      else
!-------------------------
! Quick local sum algorithm
!-------------------------
         gsum = 0.
         do j=jfirst,jlast
            do i=ifirst,ilast
               gsum = gsum + p(i,j)*area(i,j)
            enddo
         enddo
         call mp_reduce_sum(gsum)
      endif

      if ( mode==1 ) then
           g_sum = gsum / global_area
      else
           g_sum = gsum
      endif

 END FUNCTION G_SUM

 end module fv_grid_utils_tlm_mod
