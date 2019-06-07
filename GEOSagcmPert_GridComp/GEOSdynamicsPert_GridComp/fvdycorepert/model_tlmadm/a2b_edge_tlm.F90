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
module a2b_edge_tlm_mod

  use fv_grid_utils_mod, only: great_circle_dist
#ifdef VAN2
  use fv_grid_utils_mod, only: van2
#endif

  use fv_arrays_mod,     only: fv_grid_type, R_GRID

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
  public :: a2b_ord2, a2b_ord4
  public :: a2b_ord2_tlm, a2b_ord4_tlm

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

CONTAINS
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
!  Differentiation of a2b_ord4 in forward (tangent) mode:
!   variations   of useful results: qin qout
!   with respect to varying inputs: qin qout
  SUBROUTINE A2B_ORD4_TLM(qin, qin_tl, qout, qout_tl, gridstruct, npx&
&   , npy, is, ie, js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qin_tl(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(INOUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qout_tl(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local: compact 4-pt cubic
    REAL, PARAMETER :: c1=2./3.
    REAL, PARAMETER :: c2=-(1./6.)
! Parabolic spline
! real, parameter:: c1 =  0.75
! real, parameter:: c2 = -0.25
    REAL :: qx(is:ie+1, js-ng:je+ng)
    REAL :: qx_tl(is:ie+1, js-ng:je+ng)
    REAL :: qy(is-ng:ie+ng, js:je+1)
    REAL :: qy_tl(is-ng:ie+ng, js:je+1)
    REAL :: qxx(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qxx_tl(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy_tl(is-ng:ie+ng, js-ng:je+ng)
    REAL :: g_in, g_ou
    REAL :: p0(2)
    REAL :: q1(is-1:ie+1), q2(js-1:je+1)
    REAL :: q1_tl(is-1:ie+1), q2_tl(js-1:je+1)
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
    REAL :: result1_tl
    REAL :: result2
    REAL :: result2_tl
    REAL :: result3
    REAL :: result3_tl
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
        qx_tl = 0.0
        DO j=js-2,je+2
          DO i=is,ie+1
            qx_tl(i, j) = b2*(qin_tl(i-2, j)+qin_tl(i+1, j)) + b1*(&
&             qin_tl(i-1, j)+qin_tl(i, j))
            qx(i, j) = b2*(qin(i-2, j)+qin(i+1, j)) + b1*(qin(i-1, j)+&
&             qin(i, j))
          END DO
        END DO
      ELSE
        IF (gridstruct%sw_corner) THEN
          p0(1:2) = grid(1, 1, 1:2)
          result1_tl = EXTRAP_CORNER_TLM(p0, agrid(1, 1, 1:2), agrid(&
&           2, 2, 1:2), qin(1, 1), qin_tl(1, 1), qin(2, 2), qin_tl(2, 2)&
&           , result1)
          result2_tl = EXTRAP_CORNER_TLM(p0, agrid(0, 1, 1:2), agrid(&
&           -1, 2, 1:2), qin(0, 1), qin_tl(0, 1), qin(-1, 2), qin_tl(-1&
&           , 2), result2)
          result3_tl = EXTRAP_CORNER_TLM(p0, agrid(1, 0, 1:2), agrid(&
&           2, -1, 1:2), qin(1, 0), qin_tl(1, 0), qin(2, -1), qin_tl(2, &
&           -1), result3)
          qout_tl(1, 1) = r3*(result1_tl+result2_tl+result3_tl)
          qout(1, 1) = (result1+result2+result3)*r3
        END IF
        IF (gridstruct%se_corner) THEN
          p0(1:2) = grid(npx, 1, 1:2)
          result1_tl = EXTRAP_CORNER_TLM(p0, agrid(npx-1, 1, 1:2), &
&           agrid(npx-2, 2, 1:2), qin(npx-1, 1), qin_tl(npx-1, 1), qin(&
&           npx-2, 2), qin_tl(npx-2, 2), result1)
          result2_tl = EXTRAP_CORNER_TLM(p0, agrid(npx-1, 0, 1:2), &
&           agrid(npx-2, -1, 1:2), qin(npx-1, 0), qin_tl(npx-1, 0), qin(&
&           npx-2, -1), qin_tl(npx-2, -1), result2)
          result3_tl = EXTRAP_CORNER_TLM(p0, agrid(npx, 1, 1:2), &
&           agrid(npx+1, 2, 1:2), qin(npx, 1), qin_tl(npx, 1), qin(npx+1&
&           , 2), qin_tl(npx+1, 2), result3)
          qout_tl(npx, 1) = r3*(result1_tl+result2_tl+result3_tl)
          qout(npx, 1) = (result1+result2+result3)*r3
        END IF
        IF (gridstruct%ne_corner) THEN
          p0(1:2) = grid(npx, npy, 1:2)
          result1_tl = EXTRAP_CORNER_TLM(p0, agrid(npx-1, npy-1, 1:2)&
&           , agrid(npx-2, npy-2, 1:2), qin(npx-1, npy-1), qin_tl(npx-1&
&           , npy-1), qin(npx-2, npy-2), qin_tl(npx-2, npy-2), result1)
          result2_tl = EXTRAP_CORNER_TLM(p0, agrid(npx, npy-1, 1:2), &
&           agrid(npx+1, npy-2, 1:2), qin(npx, npy-1), qin_tl(npx, npy-1&
&           ), qin(npx+1, npy-2), qin_tl(npx+1, npy-2), result2)
          result3_tl = EXTRAP_CORNER_TLM(p0, agrid(npx-1, npy, 1:2), &
&           agrid(npx-2, npy+1, 1:2), qin(npx-1, npy), qin_tl(npx-1, npy&
&           ), qin(npx-2, npy+1), qin_tl(npx-2, npy+1), result3)
          qout_tl(npx, npy) = r3*(result1_tl+result2_tl+result3_tl)
          qout(npx, npy) = (result1+result2+result3)*r3
        END IF
        IF (gridstruct%nw_corner) THEN
          p0(1:2) = grid(1, npy, 1:2)
          result1_tl = EXTRAP_CORNER_TLM(p0, agrid(1, npy-1, 1:2), &
&           agrid(2, npy-2, 1:2), qin(1, npy-1), qin_tl(1, npy-1), qin(2&
&           , npy-2), qin_tl(2, npy-2), result1)
          result2_tl = EXTRAP_CORNER_TLM(p0, agrid(0, npy-1, 1:2), &
&           agrid(-1, npy-2, 1:2), qin(0, npy-1), qin_tl(0, npy-1), qin(&
&           -1, npy-2), qin_tl(-1, npy-2), result2)
          result3_tl = EXTRAP_CORNER_TLM(p0, agrid(1, npy, 1:2), &
&           agrid(2, npy+1, 1:2), qin(1, npy), qin_tl(1, npy), qin(2, &
&           npy+1), qin_tl(2, npy+1), result3)
          qout_tl(1, npy) = r3*(result1_tl+result2_tl+result3_tl)
          qout(1, npy) = (result1+result2+result3)*r3
        END IF
        IF (1 .LT. js - 2) THEN
          max1 = js - 2
        ELSE
          max1 = 1
        END IF
        IF (npy - 1 .GT. je + 2) THEN
          min1 = je + 2
          qx_tl = 0.0
        ELSE
          min1 = npy - 1
          qx_tl = 0.0
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
            qx_tl(i, j) = b2*(qin_tl(i-2, j)+qin_tl(i+1, j)) + b1*(&
&             qin_tl(i-1, j)+qin_tl(i, j))
            qx(i, j) = b2*(qin(i-2, j)+qin(i+1, j)) + b1*(qin(i-1, j)+&
&             qin(i, j))
          END DO
        END DO
! *** West Edges:
        IF (is .EQ. 1) THEN
          q2_tl = 0.0
          DO j=js1,je1
            q2_tl(j) = (dxa(1, j)*qin_tl(0, j)+dxa(0, j)*qin_tl(1, j))/(&
&             dxa(0, j)+dxa(1, j))
            q2(j) = (qin(0, j)*dxa(1, j)+qin(1, j)*dxa(0, j))/(dxa(0, j)&
&             +dxa(1, j))
          END DO
          DO j=js2,je1
            qout_tl(1, j) = edge_w(j)*q2_tl(j-1) + (1.-edge_w(j))*q2_tl(&
&             j)
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
            qx_tl(1, j) = 0.5*(((2.+g_in)*qin_tl(1, j)-qin_tl(2, j))/(1.&
&             +g_in)+((2.+g_ou)*qin_tl(0, j)-qin_tl(-1, j))/(1.+g_ou))
            qx(1, j) = 0.5*(((2.+g_in)*qin(1, j)-qin(2, j))/(1.+g_in)+((&
&             2.+g_ou)*qin(0, j)-qin(-1, j))/(1.+g_ou))
            qx_tl(2, j) = (3.*(g_in*qin_tl(1, j)+qin_tl(2, j))-g_in*&
&             qx_tl(1, j)-qx_tl(3, j))/(2.+2.*g_in)
            qx(2, j) = (3.*(g_in*qin(1, j)+qin(2, j))-(g_in*qx(1, j)+qx(&
&             3, j)))/(2.+2.*g_in)
          END DO
        ELSE
          q2_tl = 0.0
        END IF
! East Edges:
        IF (ie + 1 .EQ. npx) THEN
          DO j=js1,je1
            q2_tl(j) = (dxa(npx, j)*qin_tl(npx-1, j)+dxa(npx-1, j)*&
&             qin_tl(npx, j))/(dxa(npx-1, j)+dxa(npx, j))
            q2(j) = (qin(npx-1, j)*dxa(npx, j)+qin(npx, j)*dxa(npx-1, j)&
&             )/(dxa(npx-1, j)+dxa(npx, j))
          END DO
          DO j=js2,je1
            qout_tl(npx, j) = edge_e(j)*q2_tl(j-1) + (1.-edge_e(j))*&
&             q2_tl(j)
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
            qx_tl(npx, j) = 0.5*(((2.+g_in)*qin_tl(npx-1, j)-qin_tl(npx-&
&             2, j))/(1.+g_in)+((2.+g_ou)*qin_tl(npx, j)-qin_tl(npx+1, j&
&             ))/(1.+g_ou))
            qx(npx, j) = 0.5*(((2.+g_in)*qin(npx-1, j)-qin(npx-2, j))/(&
&             1.+g_in)+((2.+g_ou)*qin(npx, j)-qin(npx+1, j))/(1.+g_ou))
            qx_tl(npx-1, j) = (3.*(qin_tl(npx-2, j)+g_in*qin_tl(npx-1, j&
&             ))-g_in*qx_tl(npx, j)-qx_tl(npx-2, j))/(2.+2.*g_in)
            qx(npx-1, j) = (3.*(qin(npx-2, j)+g_in*qin(npx-1, j))-(g_in*&
&             qx(npx, j)+qx(npx-2, j)))/(2.+2.*g_in)
          END DO
        END IF
      END IF
!------------
! Y-Interior:
!------------
      IF (gridstruct%nested) THEN
        qy_tl = 0.0
        DO j=js,je+1
          DO i=is-2,ie+2
            qy_tl(i, j) = b2*(qin_tl(i, j-2)+qin_tl(i, j+1)) + b1*(&
&             qin_tl(i, j-1)+qin_tl(i, j))
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
          qy_tl = 0.0
        ELSE
          min5 = npy - 2
          qy_tl = 0.0
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
            qy_tl(i, j) = b2*(qin_tl(i, j-2)+qin_tl(i, j+1)) + b1*(&
&             qin_tl(i, j-1)+qin_tl(i, j))
            qy(i, j) = b2*(qin(i, j-2)+qin(i, j+1)) + b1*(qin(i, j-1)+&
&             qin(i, j))
          END DO
        END DO
! South Edges:
        IF (js .EQ. 1) THEN
          q1_tl = 0.0
          DO i=is1,ie1
            q1_tl(i) = (dya(i, 1)*qin_tl(i, 0)+dya(i, 0)*qin_tl(i, 1))/(&
&             dya(i, 0)+dya(i, 1))
            q1(i) = (qin(i, 0)*dya(i, 1)+qin(i, 1)*dya(i, 0))/(dya(i, 0)&
&             +dya(i, 1))
          END DO
          DO i=is2,ie1
            qout_tl(i, 1) = edge_s(i)*q1_tl(i-1) + (1.-edge_s(i))*q1_tl(&
&             i)
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
            qy_tl(i, 1) = 0.5*(((2.+g_in)*qin_tl(i, 1)-qin_tl(i, 2))/(1.&
&             +g_in)+((2.+g_ou)*qin_tl(i, 0)-qin_tl(i, -1))/(1.+g_ou))
            qy(i, 1) = 0.5*(((2.+g_in)*qin(i, 1)-qin(i, 2))/(1.+g_in)+((&
&             2.+g_ou)*qin(i, 0)-qin(i, -1))/(1.+g_ou))
            qy_tl(i, 2) = (3.*(g_in*qin_tl(i, 1)+qin_tl(i, 2))-g_in*&
&             qy_tl(i, 1)-qy_tl(i, 3))/(2.+2.*g_in)
            qy(i, 2) = (3.*(g_in*qin(i, 1)+qin(i, 2))-(g_in*qy(i, 1)+qy(&
&             i, 3)))/(2.+2.*g_in)
          END DO
        ELSE
          q1_tl = 0.0
        END IF
! North Edges:
        IF (je + 1 .EQ. npy) THEN
          DO i=is1,ie1
            q1_tl(i) = (dya(i, npy)*qin_tl(i, npy-1)+dya(i, npy-1)*&
&             qin_tl(i, npy))/(dya(i, npy-1)+dya(i, npy))
            q1(i) = (qin(i, npy-1)*dya(i, npy)+qin(i, npy)*dya(i, npy-1)&
&             )/(dya(i, npy-1)+dya(i, npy))
          END DO
          DO i=is2,ie1
            qout_tl(i, npy) = edge_n(i)*q1_tl(i-1) + (1.-edge_n(i))*&
&             q1_tl(i)
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
            qy_tl(i, npy) = 0.5*(((2.+g_in)*qin_tl(i, npy-1)-qin_tl(i, &
&             npy-2))/(1.+g_in)+((2.+g_ou)*qin_tl(i, npy)-qin_tl(i, npy+&
&             1))/(1.+g_ou))
            qy(i, npy) = 0.5*(((2.+g_in)*qin(i, npy-1)-qin(i, npy-2))/(&
&             1.+g_in)+((2.+g_ou)*qin(i, npy)-qin(i, npy+1))/(1.+g_ou))
            qy_tl(i, npy-1) = (3.*(qin_tl(i, npy-2)+g_in*qin_tl(i, npy-1&
&             ))-g_in*qy_tl(i, npy)-qy_tl(i, npy-2))/(2.+2.*g_in)
            qy(i, npy-1) = (3.*(qin(i, npy-2)+g_in*qin(i, npy-1))-(g_in*&
&             qy(i, npy)+qy(i, npy-2)))/(2.+2.*g_in)
          END DO
        END IF
      END IF
!--------------------------------------
      IF (gridstruct%nested) THEN
        qxx_tl = 0.0
        DO j=js,je+1
          DO i=is,ie+1
            qxx_tl(i, j) = a2*(qx_tl(i, j-2)+qx_tl(i, j+1)) + a1*(qx_tl(&
&             i, j-1)+qx_tl(i, j))
            qxx(i, j) = a2*(qx(i, j-2)+qx(i, j+1)) + a1*(qx(i, j-1)+qx(i&
&             , j))
          END DO
        END DO
        qyy_tl = 0.0
        DO j=js,je+1
          DO i=is,ie+1
            qyy_tl(i, j) = a2*(qy_tl(i-2, j)+qy_tl(i+1, j)) + a1*(qy_tl(&
&             i-1, j)+qy_tl(i, j))
            qyy(i, j) = a2*(qy(i-2, j)+qy(i+1, j)) + a1*(qy(i-1, j)+qy(i&
&             , j))
          END DO
          DO i=is,ie+1
! averaging
            qout_tl(i, j) = 0.5*(qxx_tl(i, j)+qyy_tl(i, j))
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
          qxx_tl = 0.0
        ELSE
          min9 = npy - 2
          qxx_tl = 0.0
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
            qxx_tl(i, j) = a2*(qx_tl(i, j-2)+qx_tl(i, j+1)) + a1*(qx_tl(&
&             i, j-1)+qx_tl(i, j))
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
            qxx_tl(i, 2) = c1*(qx_tl(i, 1)+qx_tl(i, 2)) + c2*(qout_tl(i&
&             , 1)+qxx_tl(i, 3))
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
            qxx_tl(i, npy-1) = c1*(qx_tl(i, npy-2)+qx_tl(i, npy-1)) + c2&
&             *(qout_tl(i, npy)+qxx_tl(i, npy-2))
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
          qyy_tl = 0.0
        ELSE
          min13 = npy - 1
          qyy_tl = 0.0
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
            qyy_tl(i, j) = a2*(qy_tl(i-2, j)+qy_tl(i+1, j)) + a1*(qy_tl(&
&             i-1, j)+qy_tl(i, j))
            qyy(i, j) = a2*(qy(i-2, j)+qy(i+1, j)) + a1*(qy(i-1, j)+qy(i&
&             , j))
          END DO
          IF (is .EQ. 1) THEN
            qyy_tl(2, j) = c1*(qy_tl(1, j)+qy_tl(2, j)) + c2*(qout_tl(1&
&             , j)+qyy_tl(3, j))
            qyy(2, j) = c1*(qy(1, j)+qy(2, j)) + c2*(qout(1, j)+qyy(3, j&
&             ))
          END IF
          IF (ie + 1 .EQ. npx) THEN
            qyy_tl(npx-1, j) = c1*(qy_tl(npx-2, j)+qy_tl(npx-1, j)) + c2&
&             *(qout_tl(npx, j)+qyy_tl(npx-2, j))
            qyy(npx-1, j) = c1*(qy(npx-2, j)+qy(npx-1, j)) + c2*(qout(&
&             npx, j)+qyy(npx-2, j))
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
          DO i=max15,min15
! averaging
            qout_tl(i, j) = 0.5*(qxx_tl(i, j)+qyy_tl(i, j))
            qout(i, j) = 0.5*(qxx(i, j)+qyy(i, j))
          END DO
        END DO
      END IF
    ELSE
      qx_tl = 0.0
! grid_type>=3
!------------------------
! Doubly periodic domain:
!------------------------
! X-sweep: PPM
      DO j=js-2,je+2
        DO i=is,ie+1
          qx_tl(i, j) = b1*(qin_tl(i-1, j)+qin_tl(i, j)) + b2*(qin_tl(i-&
&           2, j)+qin_tl(i+1, j))
          qx(i, j) = b1*(qin(i-1, j)+qin(i, j)) + b2*(qin(i-2, j)+qin(i+&
&           1, j))
        END DO
      END DO
      qy_tl = 0.0
! Y-sweep: PPM
      DO j=js,je+1
        DO i=is-2,ie+2
          qy_tl(i, j) = b1*(qin_tl(i, j-1)+qin_tl(i, j)) + b2*(qin_tl(i&
&           , j-2)+qin_tl(i, j+1))
          qy(i, j) = b1*(qin(i, j-1)+qin(i, j)) + b2*(qin(i, j-2)+qin(i&
&           , j+1))
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie+1
          qout_tl(i, j) = 0.5*(a1*(qx_tl(i, j-1)+qx_tl(i, j)+qy_tl(i-1, &
&           j)+qy_tl(i, j))+a2*(qx_tl(i, j-2)+qx_tl(i, j+1)+qy_tl(i-2, j&
&           )+qy_tl(i+1, j)))
          qout(i, j) = 0.5*(a1*(qx(i, j-1)+qx(i, j)+qy(i-1, j)+qy(i, j))&
&           +a2*(qx(i, j-2)+qx(i, j+1)+qy(i-2, j)+qy(i+1, j)))
        END DO
      END DO
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qin_tl(i, j) = qout_tl(i, j)
            qin(i, j) = qout(i, j)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE A2B_ORD4_TLM
!  Differentiation of a2b_ord2 in forward (tangent) mode:
!   variations   of useful results: qin qout
!   with respect to varying inputs: qin qout
  SUBROUTINE A2B_ORD2_TLM(qin, qin_tl, qout, qout_tl, gridstruct, npx, &
&   npy, is, ie, js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qin_tl(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(OUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(OUT) :: qout_tl(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local:
    REAL :: q1(npx), q2(npy)
    REAL :: q1_tl(npx), q2_tl(npy)
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
            qout_tl(i, j) = 0.25*(qin_tl(i-1, j-1)+qin_tl(i, j-1)+qin_tl&
&             (i-1, j)+qin_tl(i, j))
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
            qout_tl(i, j) = 0.25*(qin_tl(i-1, j-1)+qin_tl(i, j-1)+qin_tl&
&             (i-1, j)+qin_tl(i, j))
            qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin&
&             (i, j))
          END DO
        END DO
! Fix the 4 Corners:
        IF (gridstruct%sw_corner) THEN
          qout_tl(1, 1) = r3*(qin_tl(1, 1)+qin_tl(1, 0)+qin_tl(0, 1))
          qout(1, 1) = r3*(qin(1, 1)+qin(1, 0)+qin(0, 1))
        END IF
        IF (gridstruct%se_corner) THEN
          qout_tl(npx, 1) = r3*(qin_tl(npx-1, 1)+qin_tl(npx-1, 0)+qin_tl&
&           (npx, 1))
          qout(npx, 1) = r3*(qin(npx-1, 1)+qin(npx-1, 0)+qin(npx, 1))
        END IF
        IF (gridstruct%ne_corner) THEN
          qout_tl(npx, npy) = r3*(qin_tl(npx-1, npy-1)+qin_tl(npx, npy-1&
&           )+qin_tl(npx-1, npy))
          qout(npx, npy) = r3*(qin(npx-1, npy-1)+qin(npx, npy-1)+qin(npx&
&           -1, npy))
        END IF
        IF (gridstruct%nw_corner) THEN
          qout_tl(1, npy) = r3*(qin_tl(1, npy-1)+qin_tl(0, npy-1)+qin_tl&
&           (1, npy))
          qout(1, npy) = r3*(qin(1, npy-1)+qin(0, npy-1)+qin(1, npy))
        END IF
! *** West Edges:
        IF (is .EQ. 1) THEN
          q2_tl = 0.0
          DO j=js1,je1
            q2_tl(j) = 0.5*(qin_tl(0, j)+qin_tl(1, j))
            q2(j) = 0.5*(qin(0, j)+qin(1, j))
          END DO
          DO j=js2,je1
            qout_tl(1, j) = edge_w(j)*q2_tl(j-1) + (1.-edge_w(j))*q2_tl(&
&             j)
            qout(1, j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
          END DO
        ELSE
          q2_tl = 0.0
        END IF
! East Edges:
        IF (ie + 1 .EQ. npx) THEN
          DO j=js1,je1
            q2_tl(j) = 0.5*(qin_tl(npx-1, j)+qin_tl(npx, j))
            q2(j) = 0.5*(qin(npx-1, j)+qin(npx, j))
          END DO
          DO j=js2,je1
            qout_tl(npx, j) = edge_e(j)*q2_tl(j-1) + (1.-edge_e(j))*&
&             q2_tl(j)
            qout(npx, j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
          END DO
        END IF
! South Edges:
        IF (js .EQ. 1) THEN
          q1_tl = 0.0
          DO i=is1,ie1
            q1_tl(i) = 0.5*(qin_tl(i, 0)+qin_tl(i, 1))
            q1(i) = 0.5*(qin(i, 0)+qin(i, 1))
          END DO
          DO i=is2,ie1
            qout_tl(i, 1) = edge_s(i)*q1_tl(i-1) + (1.-edge_s(i))*q1_tl(&
&             i)
            qout(i, 1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
          END DO
        ELSE
          q1_tl = 0.0
        END IF
! North Edges:
        IF (je + 1 .EQ. npy) THEN
          DO i=is1,ie1
            q1_tl(i) = 0.5*(qin_tl(i, npy-1)+qin_tl(i, npy))
            q1(i) = 0.5*(qin(i, npy-1)+qin(i, npy))
          END DO
          DO i=is2,ie1
            qout_tl(i, npy) = edge_n(i)*q1_tl(i-1) + (1.-edge_n(i))*&
&             q1_tl(i)
            qout(i, npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
          END DO
        END IF
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          qout_tl(i, j) = 0.25*(qin_tl(i-1, j-1)+qin_tl(i, j-1)+qin_tl(i&
&           -1, j)+qin_tl(i, j))
          qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin(i&
&           , j))
        END DO
      END DO
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qin_tl(i, j) = qout_tl(i, j)
            qin(i, j) = qout(i, j)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE A2B_ORD2_TLM
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
  REAL FUNCTION EXTRAP_CORNER(p0, p1, p2, q1, q2)
    IMPLICIT NONE
    REAL, DIMENSION(2), INTENT(IN) :: p0, p1, p2
    REAL, INTENT(IN) :: q1, q2
    REAL :: x1, x2
    INTRINSIC REAL
    x1 = GREAT_CIRCLE_DIST(REAL(p1, kind=r_grid), REAL(p0, kind=r_grid))
    x2 = GREAT_CIRCLE_DIST(REAL(p2, kind=r_grid), REAL(p0, kind=r_grid))
    extrap_corner = q1 + x1/(x2-x1)*(q1-q2)
  END FUNCTION EXTRAP_CORNER
!  Differentiation of extrap_corner in forward (tangent) mode:
!   variations   of useful results: extrap_corner
!   with respect to varying inputs: q1 q2
  REAL FUNCTION EXTRAP_CORNER_TLM(p0, p1, p2, q1, q1_tl, q2, q2_tl, &
&   extrap_corner)
    IMPLICIT NONE
    REAL, DIMENSION(2), INTENT(IN) :: p0, p1, p2
    REAL, INTENT(IN) :: q1, q2
    REAL, INTENT(IN) :: q1_tl, q2_tl
    REAL :: x1, x2
    INTRINSIC REAL
    REAL :: extrap_corner
    x1 = GREAT_CIRCLE_DIST(REAL(p1, kind=r_grid), REAL(p0, kind=r_grid))
    x2 = GREAT_CIRCLE_DIST(REAL(p2, kind=r_grid), REAL(p0, kind=r_grid))
    extrap_corner_tlm = q1_tl + x1*(q1_tl-q2_tl)/(x2-x1)
    extrap_corner = q1 + x1/(x2-x1)*(q1-q2)
  END FUNCTION EXTRAP_CORNER_TLM
end module a2b_edge_tlm_mod
