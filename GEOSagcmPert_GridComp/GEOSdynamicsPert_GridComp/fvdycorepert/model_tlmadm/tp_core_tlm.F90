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
module tp_core_tlm_mod
!BOP
!
! !MODULE: tp_core --- A collection of routines to support FV transport
!
 use fv_mp_mod,         only: ng 
 use fv_grid_utils_mod, only: big_number
 use fv_arrays_mod,     only: fv_grid_type, fv_grid_bounds_type, r_grid

 implicit none

 private
 public fv_tp_2d, pert_ppm, copy_corners
 public fv_tp_2d_tlm, copy_corners_tlm

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
          pwx1 = damp_c*gridstruct%da_min
          pwy1 = nord + 1
          damp = pwx1**pwy1
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
          pwx1 = damp_c*gridstruct%da_min
          pwy1 = nord + 1
          damp = pwx1**pwy1
          CALL DELN_FLUX(nord, is, ie, js, je, npx, npy, damp, q, fx, fy&
&                  , gridstruct, bd)
        END IF
      END IF
    END IF
  END SUBROUTINE FV_TP_2D
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
        IF (iord .EQ. 9 .OR. iord .EQ. 13) CALL PERT_PPM(ie1 - is1 + 1, &
&                                                  q1(is1:ie1), bl(is1:&
&                                                  ie1), br(is1:ie1), 0)
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
          CALL PERT_PPM(ilast - ifirst + 1, q(ifirst:ilast, j), bl(&
&                 ifirst:ilast, j), br(ifirst:ilast, j), 0)
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
          CALL PERT_PPM(3*(ilast-ifirst+1), q(ifirst:ilast, 0:2), bl(&
&                 ifirst:ilast, 0:2), br(ifirst:ilast, 0:2), 1)
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
          CALL PERT_PPM(3*(ilast-ifirst+1), q(ifirst:ilast, npy-2:npy), &
&                 bl(ifirst:ilast, npy-2:npy), br(ifirst:ilast, npy-2:&
&                 npy), 1)
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
!  Differentiation of pert_ppm in forward (tangent) mode:
!   variations   of useful results: al ar
!   with respect to varying inputs: al ar
  SUBROUTINE PERT_PPM_TLM(im, a0, al, al_tl, ar, ar_tl, iv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    INTEGER, INTENT(IN) :: iv
    REAL, INTENT(IN) :: a0(im)
    REAL, INTENT(INOUT) :: al(im), ar(im)
    REAL, INTENT(INOUT) :: al_tl(im), ar_tl(im)
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
          al_tl(i) = 0.0
          al(i) = 0.
          ar_tl(i) = 0.0
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
                ar_tl(i) = 0.0
                ar(i) = 0.
                al_tl(i) = 0.0
                al(i) = 0.
              ELSE IF (da1 .GT. 0.) THEN
                ar_tl(i) = -(2.*al_tl(i))
                ar(i) = -(2.*al(i))
              ELSE
                al_tl(i) = -(2.*ar_tl(i))
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
            ar_tl(i) = -(2.*al_tl(i))
            ar(i) = -(2.*al(i))
          ELSE IF (a6da .GT. da2) THEN
            al_tl(i) = -(2.*ar_tl(i))
            al(i) = -(2.*ar(i))
          END IF
        ELSE
! effect of dm=0 included here
          al_tl(i) = 0.0
          al(i) = 0.
          ar_tl(i) = 0.0
          ar(i) = 0.
        END IF
      END DO
    END IF
  END SUBROUTINE PERT_PPM_TLM
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
!  Differentiation of fv_tp_2d in forward (tangent) mode:
!   variations   of useful results: q fx fy
!   with respect to varying inputs: xfx q mass mfx mfy ra_x ra_y
!                yfx fx fy crx cry
  SUBROUTINE FV_TP_2D_TLM(q, q_tl, crx, crx_tl, cry, cry_tl, npx, npy&
&   , hord, fx, fx_tl, fy, fy_tl, xfx, xfx_tl, yfx, yfx_tl, gridstruct, &
&   bd, ra_x, ra_x_tl, ra_y, ra_y_tl, mfx, mfx_tl, mfy, mfy_tl, mass, &
&   mass_tl, nord, damp_c)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL, INTENT(IN) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: crx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed)
!
    REAL, INTENT(IN) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: xfx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed)
!
    REAL, INTENT(IN) :: cry(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(IN) :: cry_tl(bd%isd:bd%ied, bd%js:bd%je+1)
!
    REAL, INTENT(IN) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(IN) :: yfx_tl(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL, INTENT(IN) :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: ra_x_tl(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL, INTENT(IN) :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL, INTENT(IN) :: ra_y_tl(bd%isd:bd%ied, bd%js:bd%je)
! transported scalar
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
! Flux in x ( E )
    REAL, INTENT(OUT) :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, INTENT(OUT) :: fx_tl(bd%is:bd%ie+1, bd%js:bd%je)
! Flux in y ( N )
    REAL, INTENT(OUT) :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, INTENT(OUT) :: fy_tl(bd%is:bd%ie, bd%js:bd%je+1)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! optional Arguments:
! Mass Flux X-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL, OPTIONAL, INTENT(IN) :: mfx_tl(bd%is:bd%ie+1, bd%js:bd%je)
! Mass Flux Y-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, OPTIONAL, INTENT(IN) :: mfy_tl(bd%is:bd%ie, bd%js:bd%je+1)
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL, INTENT(IN) :: mass_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL, INTENT(IN) :: damp_c
    INTEGER, OPTIONAL, INTENT(IN) :: nord
! Local:
    INTEGER :: ord_ou, ord_in
    REAL :: q_i(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: q_i_tl(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: q_j(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: q_j_tl(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: fx2(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: fx2_tl(bd%is:bd%ie+1, bd%jsd:bd%jed)
    REAL :: fy2(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fy2_tl(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fyy(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fyy_tl(bd%isd:bd%ied, bd%js:bd%je+1)
    REAL :: fx1(bd%is:bd%ie+1)
    REAL :: fx1_tl(bd%is:bd%ie+1)
    REAL :: damp
    INTEGER :: i, j
    INTEGER :: is, ie, js, je, isd, ied, jsd, jed
    INTRINSIC PRESENT
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
    IF (hord .EQ. 10) THEN
      ord_in = 8
    ELSE
      ord_in = hord
    END IF
    ord_ou = hord
    IF (.NOT.gridstruct%nested) CALL COPY_CORNERS_TLM(q, q_tl, npx, &
&                                                  npy, 2, gridstruct%&
&                                                  nested, bd, &
&                                                  gridstruct%sw_corner&
&                                                  , gridstruct%&
&                                                  se_corner, gridstruct&
&                                                  %nw_corner, &
&                                                  gridstruct%ne_corner)
    fy2_tl = 0.0
    CALL YPPM_TLM(fy2, fy2_tl, q, q_tl, cry, cry_tl, ord_in, isd, ied&
&              , isd, ied, js, je, jsd, jed, npx, npy, gridstruct%dya, &
&              gridstruct%nested, gridstruct%grid_type)
    fyy_tl = 0.0
    DO j=js,je+1
      DO i=isd,ied
        fyy_tl(i, j) = yfx_tl(i, j)*fy2(i, j) + yfx(i, j)*fy2_tl(i, j)
        fyy(i, j) = yfx(i, j)*fy2(i, j)
      END DO
    END DO
    q_i_tl = 0.0
    DO j=js,je
      DO i=isd,ied
        q_i_tl(i, j) = ((gridstruct%area(i, j)*q_tl(i, j)+fyy_tl(i, j)-&
&         fyy_tl(i, j+1))*ra_y(i, j)-(q(i, j)*gridstruct%area(i, j)+fyy(&
&         i, j)-fyy(i, j+1))*ra_y_tl(i, j))/ra_y(i, j)**2
        q_i(i, j) = (q(i, j)*gridstruct%area(i, j)+fyy(i, j)-fyy(i, j+1)&
&         )/ra_y(i, j)
      END DO
    END DO
    CALL XPPM_TLM(fx, fx_tl, q_i, q_i_tl, crx(is:ie+1, js:je), crx_tl&
&              (is:ie+1, js:je), ord_ou, is, ie, isd, ied, js, je, jsd, &
&              jed, npx, npy, gridstruct%dxa, gridstruct%nested, &
&              gridstruct%grid_type)
    IF (.NOT.gridstruct%nested) CALL COPY_CORNERS_TLM(q, q_tl, npx, &
&                                                  npy, 1, gridstruct%&
&                                                  nested, bd, &
&                                                  gridstruct%sw_corner&
&                                                  , gridstruct%&
&                                                  se_corner, gridstruct&
&                                                  %nw_corner, &
&                                                  gridstruct%ne_corner)
    fx2_tl = 0.0
    CALL XPPM_TLM(fx2, fx2_tl, q, q_tl, crx, crx_tl, ord_in, is, ie, &
&              isd, ied, jsd, jed, jsd, jed, npx, npy, gridstruct%dxa, &
&              gridstruct%nested, gridstruct%grid_type)
    q_j_tl = 0.0
    fx1_tl = 0.0
    DO j=jsd,jed
      DO i=is,ie+1
        fx1_tl(i) = xfx_tl(i, j)*fx2(i, j) + xfx(i, j)*fx2_tl(i, j)
        fx1(i) = xfx(i, j)*fx2(i, j)
      END DO
      DO i=is,ie
        q_j_tl(i, j) = ((gridstruct%area(i, j)*q_tl(i, j)+fx1_tl(i)-&
&         fx1_tl(i+1))*ra_x(i, j)-(q(i, j)*gridstruct%area(i, j)+fx1(i)-&
&         fx1(i+1))*ra_x_tl(i, j))/ra_x(i, j)**2
        q_j(i, j) = (q(i, j)*gridstruct%area(i, j)+fx1(i)-fx1(i+1))/ra_x&
&         (i, j)
      END DO
    END DO
    CALL YPPM_TLM(fy, fy_tl, q_j, q_j_tl, cry, cry_tl, ord_ou, is, ie&
&              , isd, ied, js, je, jsd, jed, npx, npy, gridstruct%dya, &
&              gridstruct%nested, gridstruct%grid_type)
!----------------
! Flux averaging:
!----------------
    IF (PRESENT(mfx) .AND. PRESENT(mfy)) THEN
!---------------------------------
! For transport of pt and tracers
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          fx_tl(i, j) = 0.5*((fx_tl(i, j)+fx2_tl(i, j))*mfx(i, j)+(fx(i&
&           , j)+fx2(i, j))*mfx_tl(i, j))
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*mfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy_tl(i, j) = 0.5*((fy_tl(i, j)+fy2_tl(i, j))*mfy(i, j)+(fy(i&
&           , j)+fy2(i, j))*mfy_tl(i, j))
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*mfy(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c) .AND. PRESENT(mass)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          pwx1 = damp_c*gridstruct%da_min
          pwy1 = nord + 1
          damp = pwx1**pwy1
          CALL DELN_FLUX_TLM(nord, is, ie, js, je, npx, npy, damp, q&
&                         , q_tl, fx, fx_tl, fy, fy_tl, gridstruct, bd, &
&                         mass, mass_tl)
        END IF
      END IF
    ELSE
!---------------------------------
! For transport of delp, vorticity
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          fx_tl(i, j) = 0.5*((fx_tl(i, j)+fx2_tl(i, j))*xfx(i, j)+(fx(i&
&           , j)+fx2(i, j))*xfx_tl(i, j))
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*xfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy_tl(i, j) = 0.5*((fy_tl(i, j)+fy2_tl(i, j))*yfx(i, j)+(fy(i&
&           , j)+fy2(i, j))*yfx_tl(i, j))
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*yfx(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          pwx1 = damp_c*gridstruct%da_min
          pwy1 = nord + 1
          damp = pwx1**pwy1
          CALL DELN_FLUX_TLM(nord, is, ie, js, je, npx, npy, damp, q&
&                         , q_tl, fx, fx_tl, fy, fy_tl, gridstruct, bd)
        END IF
      END IF
    END IF
  END SUBROUTINE FV_TP_2D_TLM
!  Differentiation of xppm in forward (tangent) mode:
!   variations   of useful results: flux
!   with respect to varying inputs: q flux c
  SUBROUTINE XPPM_TLM(flux, flux_tl, q, q_tl, c, c_tl, iord, is, ie, &
&   isd, ied, jfirst, jlast, jsd, jed, npx, npy, dxa, nested, grid_type)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, isd, ied, jsd, jed
! compute domain
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: iord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(isd:ied, jfirst:jlast)
    REAL, INTENT(IN) :: q_tl(isd:ied, jfirst:jlast)
! Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: c_tl(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: dxa(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
! !OUTPUT PARAMETERS:
!  Flux
    REAL, INTENT(OUT) :: flux(is:ie+1, jfirst:jlast)
    REAL, INTENT(OUT) :: flux_tl(is:ie+1, jfirst:jlast)
! Local
    REAL, DIMENSION(is-1:ie+1) :: bl, br, b0
    REAL :: q1(isd:ied)
    REAL :: q1_tl(isd:ied)
    REAL, DIMENSION(is:ie+1) :: fx0, fx1
    LOGICAL, DIMENSION(is-1:ie+1) :: smt5, smt6
    REAL :: al(is-1:ie+2)
    REAL :: al_tl(is-1:ie+2)
    REAL :: dm(is-2:ie+2)
    REAL :: dq(is-3:ie+2)
    INTEGER :: i, j, ie3, is1, ie1
    REAL :: x0, x1, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2
    REAL :: xt_tl, qtmp_tl
    INTRINSIC MAX
    INTRINSIC MIN
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
      al_tl = 0.0
      q1_tl = 0.0
    ELSE
      is1 = is - 1
      ie3 = ie + 2
      ie1 = ie + 1
      al_tl = 0.0
      q1_tl = 0.0
    END IF
    DO j=jfirst,jlast
      DO i=isd,ied
        q1_tl(i) = q_tl(i, j)
        q1(i) = q(i, j)
      END DO
      IF (iord .LT. 8 .OR. iord .EQ. 333) THEN
! ord = 2: perfectly linear ppm scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6
        DO i=is1,ie3
          al_tl(i) = p1*(q1_tl(i-1)+q1_tl(i)) + p2*(q1_tl(i-2)+q1_tl(i+1&
&           ))
          al(i) = p1*(q1(i-1)+q1(i)) + p2*(q1(i-2)+q1(i+1))
        END DO
        IF (.NOT.nested .AND. grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            al_tl(0) = c1*q1_tl(-2) + c2*q1_tl(-1) + c3*q1_tl(0)
            al(0) = c1*q1(-2) + c2*q1(-1) + c3*q1(0)
            al_tl(1) = 0.5*(((2.*dxa(0, j)+dxa(-1, j))*q1_tl(0)-dxa(0, j&
&             )*q1_tl(-1))/(dxa(-1, j)+dxa(0, j))+((2.*dxa(1, j)+dxa(2, &
&             j))*q1_tl(1)-dxa(1, j)*q1_tl(2))/(dxa(1, j)+dxa(2, j)))
            al(1) = 0.5*(((2.*dxa(0, j)+dxa(-1, j))*q1(0)-dxa(0, j)*q1(-&
&             1))/(dxa(-1, j)+dxa(0, j))+((2.*dxa(1, j)+dxa(2, j))*q1(1)&
&             -dxa(1, j)*q1(2))/(dxa(1, j)+dxa(2, j)))
            al_tl(2) = c3*q1_tl(1) + c2*q1_tl(2) + c1*q1_tl(3)
            al(2) = c3*q1(1) + c2*q1(2) + c1*q1(3)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            al_tl(npx-1) = c1*q1_tl(npx-3) + c2*q1_tl(npx-2) + c3*q1_tl(&
&             npx-1)
            al(npx-1) = c1*q1(npx-3) + c2*q1(npx-2) + c3*q1(npx-1)
            al_tl(npx) = 0.5*(((2.*dxa(npx-1, j)+dxa(npx-2, j))*q1_tl(&
&             npx-1)-dxa(npx-1, j)*q1_tl(npx-2))/(dxa(npx-2, j)+dxa(npx-&
&             1, j))+((2.*dxa(npx, j)+dxa(npx+1, j))*q1_tl(npx)-dxa(npx&
&             , j)*q1_tl(npx+1))/(dxa(npx, j)+dxa(npx+1, j)))
            al(npx) = 0.5*(((2.*dxa(npx-1, j)+dxa(npx-2, j))*q1(npx-1)-&
&             dxa(npx-1, j)*q1(npx-2))/(dxa(npx-2, j)+dxa(npx-1, j))+((&
&             2.*dxa(npx, j)+dxa(npx+1, j))*q1(npx)-dxa(npx, j)*q1(npx+1&
&             ))/(dxa(npx, j)+dxa(npx+1, j)))
            al_tl(npx+1) = c3*q1_tl(npx) + c2*q1_tl(npx+1) + c1*q1_tl(&
&             npx+2)
            al(npx+1) = c3*q1(npx) + c2*q1(npx+1) + c1*q1(npx+2)
          END IF
        END IF
        IF (iord .EQ. 1) THEN
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              flux_tl(i, j) = q1_tl(i-1)
              flux(i, j) = q1(i-1)
            ELSE
              flux_tl(i, j) = q1_tl(i)
              flux(i, j) = q1(i)
            END IF
          END DO
        ELSE IF (iord .EQ. 2) THEN
! perfectly linear scheme
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            xt_tl = c_tl(i, j)
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              qtmp_tl = q1_tl(i-1)
              qtmp = q1(i-1)
              flux_tl(i, j) = qtmp_tl + (1.-xt)*(al_tl(i)-qtmp_tl-xt_tl*&
&               (al(i-1)+al(i)-(qtmp+qtmp))-xt*(al_tl(i-1)+al_tl(i)-2*&
&               qtmp_tl)) - xt_tl*(al(i)-qtmp-xt*(al(i-1)+al(i)-(qtmp+&
&               qtmp)))
              flux(i, j) = qtmp + (1.-xt)*(al(i)-qtmp-xt*(al(i-1)+al(i)-&
&               (qtmp+qtmp)))
            ELSE
              qtmp_tl = q1_tl(i)
              qtmp = q1(i)
              flux_tl(i, j) = qtmp_tl + xt_tl*(al(i)-qtmp+xt*(al(i)+al(i&
&               +1)-(qtmp+qtmp))) + (1.+xt)*(al_tl(i)-qtmp_tl+xt_tl*(al(&
&               i)+al(i+1)-(qtmp+qtmp))+xt*(al_tl(i)+al_tl(i+1)-2*&
&               qtmp_tl))
              flux(i, j) = qtmp + (1.+xt)*(al(i)-qtmp+xt*(al(i)+al(i+1)-&
&               (qtmp+qtmp)))
            END IF
          END DO
        ELSE IF (iord .EQ. 333) THEN
! Perfectly linear scheme, more diffusive than ord=2 (HoldawayKent-2015-TellusA)
!DEC$ VECTOR ALWAYS
          DO i=is,ie+1
            xt_tl = c_tl(i, j)
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              flux_tl(i, j) = (2.0*q1_tl(i)+5.0*q1_tl(i-1)-q1_tl(i-2))/&
&               6.0 - 0.5*(xt_tl*(q1(i)-q1(i-1))+xt*(q1_tl(i)-q1_tl(i-1)&
&               )) + (xt_tl*xt+xt*xt_tl)*(q1(i)-2.0*q1(i-1)+q1(i-2))/6.0&
&               + xt**2*(q1_tl(i)-2.0*q1_tl(i-1)+q1_tl(i-2))/6.0
              flux(i, j) = (2.0*q1(i)+5.0*q1(i-1)-q1(i-2))/6.0 - 0.5*xt*&
&               (q1(i)-q1(i-1)) + xt*xt/6.0*(q1(i)-2.0*q1(i-1)+q1(i-2))
            ELSE
              flux_tl(i, j) = (2.0*q1_tl(i-1)+5.0*q1_tl(i)-q1_tl(i+1))/&
&               6.0 - 0.5*(xt_tl*(q1(i)-q1(i-1))+xt*(q1_tl(i)-q1_tl(i-1)&
&               )) + (xt_tl*xt+xt*xt_tl)*(q1(i+1)-2.0*q1(i)+q1(i-1))/6.0&
&               + xt**2*(q1_tl(i+1)-2.0*q1_tl(i)+q1_tl(i-1))/6.0
              flux(i, j) = (2.0*q1(i-1)+5.0*q1(i)-q1(i+1))/6.0 - 0.5*xt*&
&               (q1(i)-q1(i-1)) + xt*xt/6.0*(q1(i+1)-2.0*q1(i)+q1(i-1))
            END IF
          END DO
        END IF
      END IF
    END DO
  END SUBROUTINE XPPM_TLM
!  Differentiation of yppm in forward (tangent) mode:
!   variations   of useful results: flux
!   with respect to varying inputs: q flux c
  SUBROUTINE YPPM_TLM(flux, flux_tl, q, q_tl, c, c_tl, jord, ifirst, &
&   ilast, isd, ied, js, je, jsd, jed, npx, npy, dya, nested, grid_type)
    IMPLICIT NONE
! Compute domain
    INTEGER, INTENT(IN) :: ifirst, ilast
    INTEGER, INTENT(IN) :: isd, ied, js, je, jsd, jed
    INTEGER, INTENT(IN) :: jord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst:ilast, jsd:jed)
    REAL, INTENT(IN) :: q_tl(ifirst:ilast, jsd:jed)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
    REAL, INTENT(IN) :: c_tl(isd:ied, js:je+1)
!  Flux
    REAL, INTENT(OUT) :: flux(ifirst:ilast, js:je+1)
    REAL, INTENT(OUT) :: flux_tl(ifirst:ilast, js:je+1)
    REAL, INTENT(IN) :: dya(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
    INTEGER, INTENT(IN) :: grid_type
! Local:
    REAL :: dm(ifirst:ilast, js-2:je+2)
    REAL :: al(ifirst:ilast, js-1:je+2)
    REAL :: al_tl(ifirst:ilast, js-1:je+2)
    REAL, DIMENSION(ifirst:ilast, js-1:je+1) :: bl, br, b0
    REAL :: dq(ifirst:ilast, js-3:je+2)
    REAL, DIMENSION(ifirst:ilast) :: fx0, fx1
    LOGICAL, DIMENSION(ifirst:ilast, js-1:je+1) :: smt5, smt6
    REAL :: x0, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2, r1
    REAL :: xt_tl, qtmp_tl
    INTEGER :: i, j, js1, je3, je1
    INTRINSIC MAX
    INTRINSIC MIN
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
      al_tl = 0.0
      DO j=js1,je3
        DO i=ifirst,ilast
          al_tl(i, j) = p1*(q_tl(i, j-1)+q_tl(i, j)) + p2*(q_tl(i, j-2)+&
&           q_tl(i, j+1))
          al(i, j) = p1*(q(i, j-1)+q(i, j)) + p2*(q(i, j-2)+q(i, j+1))
        END DO
      END DO
      IF (.NOT.nested .AND. grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            al_tl(i, 0) = c1*q_tl(i, -2) + c2*q_tl(i, -1) + c3*q_tl(i, 0&
&             )
            al(i, 0) = c1*q(i, -2) + c2*q(i, -1) + c3*q(i, 0)
            al_tl(i, 1) = 0.5*(((2.*dya(i, 0)+dya(i, -1))*q_tl(i, 0)-dya&
&             (i, 0)*q_tl(i, -1))/(dya(i, -1)+dya(i, 0))+((2.*dya(i, 1)+&
&             dya(i, 2))*q_tl(i, 1)-dya(i, 1)*q_tl(i, 2))/(dya(i, 1)+dya&
&             (i, 2)))
            al(i, 1) = 0.5*(((2.*dya(i, 0)+dya(i, -1))*q(i, 0)-dya(i, 0)&
&             *q(i, -1))/(dya(i, -1)+dya(i, 0))+((2.*dya(i, 1)+dya(i, 2)&
&             )*q(i, 1)-dya(i, 1)*q(i, 2))/(dya(i, 1)+dya(i, 2)))
            al_tl(i, 2) = c3*q_tl(i, 1) + c2*q_tl(i, 2) + c1*q_tl(i, 3)
            al(i, 2) = c3*q(i, 1) + c2*q(i, 2) + c1*q(i, 3)
          END DO
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            al_tl(i, npy-1) = c1*q_tl(i, npy-3) + c2*q_tl(i, npy-2) + c3&
&             *q_tl(i, npy-1)
            al(i, npy-1) = c1*q(i, npy-3) + c2*q(i, npy-2) + c3*q(i, npy&
&             -1)
            al_tl(i, npy) = 0.5*(((2.*dya(i, npy-1)+dya(i, npy-2))*q_tl(&
&             i, npy-1)-dya(i, npy-1)*q_tl(i, npy-2))/(dya(i, npy-2)+dya&
&             (i, npy-1))+((2.*dya(i, npy)+dya(i, npy+1))*q_tl(i, npy)-&
&             dya(i, npy)*q_tl(i, npy+1))/(dya(i, npy)+dya(i, npy+1)))
            al(i, npy) = 0.5*(((2.*dya(i, npy-1)+dya(i, npy-2))*q(i, npy&
&             -1)-dya(i, npy-1)*q(i, npy-2))/(dya(i, npy-2)+dya(i, npy-1&
&             ))+((2.*dya(i, npy)+dya(i, npy+1))*q(i, npy)-dya(i, npy)*q&
&             (i, npy+1))/(dya(i, npy)+dya(i, npy+1)))
            al_tl(i, npy+1) = c3*q_tl(i, npy) + c2*q_tl(i, npy+1) + c1*&
&             q_tl(i, npy+2)
            al(i, npy+1) = c3*q(i, npy) + c2*q(i, npy+1) + c1*q(i, npy+2&
&             )
          END DO
        END IF
      END IF
      IF (jord .EQ. 1) THEN
        DO j=js,je+1
          DO i=ifirst,ilast
            IF (c(i, j) .GT. 0.) THEN
              flux_tl(i, j) = q_tl(i, j-1)
              flux(i, j) = q(i, j-1)
            ELSE
              flux_tl(i, j) = q_tl(i, j)
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
            xt_tl = c_tl(i, j)
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              qtmp_tl = q_tl(i, j-1)
              qtmp = q(i, j-1)
              flux_tl(i, j) = qtmp_tl + (1.-xt)*(al_tl(i, j)-qtmp_tl-&
&               xt_tl*(al(i, j-1)+al(i, j)-(qtmp+qtmp))-xt*(al_tl(i, j-1&
&               )+al_tl(i, j)-2*qtmp_tl)) - xt_tl*(al(i, j)-qtmp-xt*(al(&
&               i, j-1)+al(i, j)-(qtmp+qtmp)))
              flux(i, j) = qtmp + (1.-xt)*(al(i, j)-qtmp-xt*(al(i, j-1)+&
&               al(i, j)-(qtmp+qtmp)))
            ELSE
              qtmp_tl = q_tl(i, j)
              qtmp = q(i, j)
              flux_tl(i, j) = qtmp_tl + xt_tl*(al(i, j)-qtmp+xt*(al(i, j&
&               )+al(i, j+1)-(qtmp+qtmp))) + (1.+xt)*(al_tl(i, j)-&
&               qtmp_tl+xt_tl*(al(i, j)+al(i, j+1)-(qtmp+qtmp))+xt*(&
&               al_tl(i, j)+al_tl(i, j+1)-2*qtmp_tl))
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
            xt_tl = c_tl(i, j)
            xt = c(i, j)
            IF (xt .GT. 0.) THEN
              flux_tl(i, j) = (2.0*q_tl(i, j)+5.0*q_tl(i, j-1)-q_tl(i, j&
&               -2))/6.0 - 0.5*(xt_tl*(q(i, j)-q(i, j-1))+xt*(q_tl(i, j)&
&               -q_tl(i, j-1))) + (xt_tl*xt+xt*xt_tl)*(q(i, j)-2.0*q(i, &
&               j-1)+q(i, j-2))/6.0 + xt**2*(q_tl(i, j)-2.0*q_tl(i, j-1)&
&               +q_tl(i, j-2))/6.0
              flux(i, j) = (2.0*q(i, j)+5.0*q(i, j-1)-q(i, j-2))/6.0 - &
&               0.5*xt*(q(i, j)-q(i, j-1)) + xt*xt/6.0*(q(i, j)-2.0*q(i&
&               , j-1)+q(i, j-2))
            ELSE
              flux_tl(i, j) = (2.0*q_tl(i, j-1)+5.0*q_tl(i, j)-q_tl(i, j&
&               +1))/6.0 - 0.5*(xt_tl*(q(i, j)-q(i, j-1))+xt*(q_tl(i, j)&
&               -q_tl(i, j-1))) + (xt_tl*xt+xt*xt_tl)*(q(i, j+1)-2.0*q(i&
&               , j)+q(i, j-1))/6.0 + xt**2*(q_tl(i, j+1)-2.0*q_tl(i, j)&
&               +q_tl(i, j-1))/6.0
              flux(i, j) = (2.0*q(i, j-1)+5.0*q(i, j)-q(i, j+1))/6.0 - &
&               0.5*xt*(q(i, j)-q(i, j-1)) + xt*xt/6.0*(q(i, j+1)-2.0*q(&
&               i, j)+q(i, j-1))
            END IF
          END DO
        END DO
      END IF
      RETURN
    END IF
  END SUBROUTINE YPPM_TLM
!  Differentiation of deln_flux in forward (tangent) mode:
!   variations   of useful results: fx fy
!   with respect to varying inputs: q mass fx fy
  SUBROUTINE DELN_FLUX_TLM(nord, is, ie, js, je, npx, npy, damp, q, &
&   q_tl, fx, fx_tl, fy, fy_tl, gridstruct, bd, mass, mass_tl)
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
    REAL, INTENT(IN) :: q_tl(bd%is-ng:bd%ie+ng, bd%js-ng:bd%je+ng)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
! q ghosted on input
    REAL, OPTIONAL, INTENT(IN) :: mass(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, OPTIONAL, INTENT(IN) :: mass_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
! diffusive fluxes:
    REAL, INTENT(INOUT) :: fx(bd%is:bd%ie+1, bd%js:bd%je), fy(bd%is:bd%&
&   ie, bd%js:bd%je+1)
    REAL, INTENT(INOUT) :: fx_tl(bd%is:bd%ie+1, bd%js:bd%je), fy_tl(bd%&
&   is:bd%ie, bd%js:bd%je+1)
! local:
    REAL :: fx2(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2(bd%isd:bd%ied, bd%&
&   jsd:bd%jed+1)
    REAL :: fx2_tl(bd%isd:bd%ied+1, bd%jsd:bd%jed), fy2_tl(bd%isd:bd%ied&
&   , bd%jsd:bd%jed+1)
    REAL :: d2(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: d2_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL :: damp2
    INTEGER :: i, j, n, nt, i1, i2, j1, j2
    INTRINSIC PRESENT
    i1 = is - 1 - nord
    i2 = ie + 1 + nord
    j1 = js - 1 - nord
    j2 = je + 1 + nord
    IF (.NOT.PRESENT(mass)) THEN
      d2_tl = 0.0
      DO j=j1,j2
        DO i=i1,i2
          d2_tl(i, j) = damp*q_tl(i, j)
          d2(i, j) = damp*q(i, j)
        END DO
      END DO
    ELSE
      d2_tl = 0.0
      DO j=j1,j2
        DO i=i1,i2
          d2_tl(i, j) = q_tl(i, j)
          d2(i, j) = q(i, j)
        END DO
      END DO
    END IF
    IF (nord .GT. 0) THEN
      CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 1, gridstruct%nested&
&                        , bd, gridstruct%sw_corner, gridstruct%&
&                        se_corner, gridstruct%nw_corner, gridstruct%&
&                        ne_corner)
      fx2_tl = 0.0
    ELSE
      fx2_tl = 0.0
    END IF
    DO j=js-nord,je+nord
      DO i=is-nord,ie+nord+1
        fx2_tl(i, j) = gridstruct%del6_v(i, j)*(d2_tl(i-1, j)-d2_tl(i, j&
&         ))
        fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i-1, j)-d2(i, j))
      END DO
    END DO
    IF (nord .GT. 0) THEN
      CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 2, gridstruct%nested&
&                        , bd, gridstruct%sw_corner, gridstruct%&
&                        se_corner, gridstruct%nw_corner, gridstruct%&
&                        ne_corner)
      fy2_tl = 0.0
    ELSE
      fy2_tl = 0.0
    END IF
    DO j=js-nord,je+nord+1
      DO i=is-nord,ie+nord
        fy2_tl(i, j) = gridstruct%del6_u(i, j)*(d2_tl(i, j-1)-d2_tl(i, j&
&         ))
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
            d2_tl(i, j) = gridstruct%rarea(i, j)*(fx2_tl(i, j)-fx2_tl(i+&
&             1, j)+fy2_tl(i, j)-fy2_tl(i, j+1))
            d2(i, j) = (fx2(i, j)-fx2(i+1, j)+fy2(i, j)-fy2(i, j+1))*&
&             gridstruct%rarea(i, j)
          END DO
        END DO
        CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 1, gridstruct%&
&                          nested, bd, gridstruct%sw_corner, gridstruct%&
&                          se_corner, gridstruct%nw_corner, gridstruct%&
&                          ne_corner)
        DO j=js-nt,je+nt
          DO i=is-nt,ie+nt+1
            fx2_tl(i, j) = gridstruct%del6_v(i, j)*(d2_tl(i, j)-d2_tl(i-&
&             1, j))
            fx2(i, j) = gridstruct%del6_v(i, j)*(d2(i, j)-d2(i-1, j))
          END DO
        END DO
        CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 2, gridstruct%&
&                          nested, bd, gridstruct%sw_corner, gridstruct%&
&                          se_corner, gridstruct%nw_corner, gridstruct%&
&                          ne_corner)
        DO j=js-nt,je+nt+1
          DO i=is-nt,ie+nt
            fy2_tl(i, j) = gridstruct%del6_u(i, j)*(d2_tl(i, j)-d2_tl(i&
&             , j-1))
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
          fx_tl(i, j) = fx_tl(i, j) + damp2*((mass_tl(i-1, j)+mass_tl(i&
&           , j))*fx2(i, j)+(mass(i-1, j)+mass(i, j))*fx2_tl(i, j))
          fx(i, j) = fx(i, j) + damp2*(mass(i-1, j)+mass(i, j))*fx2(i, j&
&           )
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy_tl(i, j) = fy_tl(i, j) + damp2*((mass_tl(i, j-1)+mass_tl(i&
&           , j))*fy2(i, j)+(mass(i, j-1)+mass(i, j))*fy2_tl(i, j))
          fy(i, j) = fy(i, j) + damp2*(mass(i, j-1)+mass(i, j))*fy2(i, j&
&           )
        END DO
      END DO
    ELSE
      DO j=js,je
        DO i=is,ie+1
          fx_tl(i, j) = fx_tl(i, j) + fx2_tl(i, j)
          fx(i, j) = fx(i, j) + fx2(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy_tl(i, j) = fy_tl(i, j) + fy2_tl(i, j)
          fy(i, j) = fy(i, j) + fy2(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE DELN_FLUX_TLM
!  Differentiation of copy_corners in forward (tangent) mode:
!   variations   of useful results: q
!   with respect to varying inputs: q
!Weird arguments are because this routine is called in a lot of
!places outside of tp_core, sometimes very deeply nested in the call tree.
  SUBROUTINE COPY_CORNERS_TLM(q, q_tl, npx, npy, dir, nested, bd, &
&   sw_corner, se_corner, nw_corner, ne_corner)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npx, npy, dir
    REAL, INTENT(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed)
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed)
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
            q_tl(i, j) = q_tl(j, 1-i)
            q(i, j) = q(j, 1-i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q_tl(i, j) = q_tl(npy-j, i-npx+1)
            q(i, j) = q(npy-j, i-npx+1)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q_tl(i, j) = q_tl(j, 2*npx-1-i)
            q(i, j) = q(j, 2*npx-1-i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q_tl(i, j) = q_tl(npy-j, i-1+npx)
            q(i, j) = q(npy-j, i-1+npx)
          END DO
        END DO
      END IF
    ELSE IF (dir .EQ. 2) THEN
! YDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            q_tl(i, j) = q_tl(1-j, i)
            q(i, j) = q(1-j, i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q_tl(i, j) = q_tl(npy+j-1, npx-i)
            q(i, j) = q(npy+j-1, npx-i)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q_tl(i, j) = q_tl(2*npy-1-j, i)
            q(i, j) = q(2*npy-1-j, i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q_tl(i, j) = q_tl(j+1-npx, npy-i)
            q(i, j) = q(j+1-npx, npy-i)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE COPY_CORNERS_TLM
!Weird arguments are because this routine is called in a lot of
!places outside of tp_core, sometimes very deeply nested in the call tree.
end module tp_core_tlm_mod
