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
module fv_tracer2d_tlm_mod
   use tp_core_tlm_mod,   only: fv_tp_2d, copy_corners
   use tp_core_tlm_mod,   only: fv_tp_2d_tlm, copy_corners_tlm
   use fv_mp_mod,         only: mp_reduce_max
   use fv_mp_mod,         only: ng, mp_gather, is_master
   use fv_mp_mod,         only: group_halo_update_type
   use fv_mp_tlm_mod,     only: start_group_halo_update, complete_group_halo_update
   use fv_mp_tlm_mod,     only: start_group_halo_update_tlm
   use mpp_domains_mod,   only: mpp_update_domains, mpp_get_boundary, CGRID_NE, domain2d
   use fv_mp_tlm_mod,     only: mpp_update_domains_tlm, mpp_get_boundary_tlm
   use fv_timing_mod,     only: timing_on, timing_off
   use boundary_tlm_mod,  only: nested_grid_BC_apply_intT
   use boundary_tlm_mod,  only: nested_grid_BC_apply_intT_tlm
   use fv_arrays_mod,     only: fv_grid_type, fv_flags_type, fv_nest_type, fv_atmos_type, fv_grid_bounds_type
   use mpp_mod,           only: mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum, mpp_max

implicit none
private

public :: tracer_2d, tracer_2d_nested, tracer_2d_1L
public :: tracer_2d_tlm, tracer_2d_nested_tlm, tracer_2d_1L_tlm

real, allocatable, dimension(:,:,:) :: nest_fx_west_accum, nest_fx_east_accum, nest_fx_south_accum, nest_fx_north_accum

!---- version number -----
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of tracer_2d_1l in forward (tangent) mode:
!   variations   of useful results: q dp1
!   with respect to varying inputs: q dp1 mfx mfy cx cy
!-----------------------------------------------------------------------
! !ROUTINE: Perform 2D horizontal-to-lagrangian transport
!-----------------------------------------------------------------------
  SUBROUTINE TRACER_2D_1L_TLM(q, q_tl, dp1, dp1_tl, mfx, mfx_tl, mfy, &
&   mfy_tl, cx, cx_tl, cy, cy_tl, gridstruct, bd, domain, npx, npy, npz&
&   , nq, hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, hord_pert, &
&   nord_tr_pert, trdm_pert, split_damp_tr, dpa)
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
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dp1_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_tl(bd%is:bd%ie+1, bd%js:bd%je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_tl(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
! DELP after advection
    REAL, OPTIONAL, INTENT(OUT) :: dpa(bd%is:bd%ie, bd%js:bd%je)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Local Arrays
! 3D tracers
    REAL :: qn2(bd%isd:bd%ied, bd%jsd:bd%jed, nq)
    REAL :: qn2_tl(bd%isd:bd%ied, bd%jsd:bd%jed, nq)
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: dp2_tl(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_tl(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_tl(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_tl(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_tl(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
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
    xfx_tl = 0.0
    yfx_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx,cmax)
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx_tl(i, j, k) = dxa(i-1, j)*dy(i, j)*sin_sg(i-1, j, 3)*&
&             cx_tl(i, j, k)
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sin_sg(i-1, &
&             j, 3)
          ELSE
            xfx_tl(i, j, k) = dxa(i, j)*dy(i, j)*sin_sg(i, j, 1)*cx_tl(i&
&             , j, k)
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1&
&             )
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx_tl(i, j, k) = dya(i, j-1)*dx(i, j)*sin_sg(i, j-1, 4)*&
&             cy_tl(i, j, k)
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sin_sg(i, j-&
&             1, 4)
          ELSE
            yfx_tl(i, j, k) = dya(i, j)*dx(i, j)*sin_sg(i, j, 2)*cy_tl(i&
&             , j, k)
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
            cx_tl(i, j, k) = frac*cx_tl(i, j, k)
            cx(i, j, k) = cx(i, j, k)*frac
            xfx_tl(i, j, k) = frac*xfx_tl(i, j, k)
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            mfx_tl(i, j, k) = frac*mfx_tl(i, j, k)
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            cy_tl(i, j, k) = frac*cy_tl(i, j, k)
            cy(i, j, k) = cy(i, j, k)*frac
            yfx_tl(i, j, k) = frac*yfx_tl(i, j, k)
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy_tl(i, j, k) = frac*mfy_tl(i, j, k)
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
    dp2_tl = 0.0
    qn2_tl = 0.0
    ra_x_tl = 0.0
    ra_y_tl = 0.0
    fx_tl = 0.0
    fy_tl = 0.0
! Begin k-independent tracer transport; can not be OpenMPed because the mpp_update call.
    DO k=1,npz
!$OMP parallel do default(none) shared(k,is,ie,js,je,isd,ied,jsd,jed,xfx,area,yfx,ra_x,ra_y)
      DO j=jsd,jed
        DO i=is,ie
          ra_x_tl(i, j) = xfx_tl(i, j, k) - xfx_tl(i+1, j, k)
          ra_x(i, j) = area(i, j) + (xfx(i, j, k)-xfx(i+1, j, k))
        END DO
        IF (j .GE. js .AND. j .LE. je) THEN
          DO i=isd,ied
            ra_y_tl(i, j) = yfx_tl(i, j, k) - yfx_tl(i, j+1, k)
            ra_y(i, j) = area(i, j) + (yfx(i, j, k)-yfx(i, j+1, k))
          END DO
        END IF
      END DO
      nsplt = INT(1. + cmax(k))
      DO it=1,nsplt
!$OMP parallel do default(none) shared(k,is,ie,js,je,rarea,mfx,mfy,dp1,dp2)
        DO j=js,je
          DO i=is,ie
            dp2_tl(i, j) = dp1_tl(i, j, k) + rarea(i, j)*(mfx_tl(i, j, k&
&             )-mfx_tl(i+1, j, k)+mfy_tl(i, j, k)-mfy_tl(i, j+1, k))
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
                  qn2_tl(i, j, iq) = q_tl(i, j, k, iq)
                  qn2(i, j, iq) = q(i, j, k, iq)
                END DO
              END DO
            END IF
            IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D_TLM(qn2(isd:ied, jsd:jed, iq), qn2_tl(isd&
&                            :ied, jsd:jed, iq), cx(is:ie+1, jsd:jed, k)&
&                            , cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied, &
&                            js:je+1, k), cy_tl(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fx_tl, fy, fy_tl, xfx(&
&                            is:ie+1, jsd:jed, k), xfx_tl(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), yfx_tl(&
&                            isd:ied, js:je+1, k), gridstruct, bd, ra_x&
&                            , ra_x_tl, ra_y, ra_y_tl, mfx=mfx(is:ie+1, &
&                            js:je, k), mfx_tl=mfx_tl(is:ie+1, js:je, k)&
&                            , mfy=mfy(is:ie, js:je+1, k), mfy_tl=mfy_tl&
&                            (is:ie, js:je+1, k))
            ELSE
              CALL FV_TP_2D_TLM(qn2(isd:ied, jsd:jed, iq), qn2_tl(isd:&
&                         ied, jsd:jed, iq), cx(is:ie+1, jsd:jed, k), &
&                         cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied, js:je+&
&                         1, k), cy_tl(isd:ied, js:je+1, k), npx, npy, &
&                         hord_pert, fx, fx_tl, fy, fy_tl, xfx(is:ie+1, &
&                         jsd:jed, k), xfx_tl(is:ie+1, jsd:jed, k), yfx(&
&                         isd:ied, js:je+1, k), yfx_tl(isd:ied, js:je+1&
&                         , k), gridstruct, bd, ra_x, ra_x_tl, ra_y, &
&                         ra_y_tl, mfx=mfx(is:ie+1, js:je, k), mfx_tl=&
&                         mfx_tl(is:ie+1, js:je, k), mfy=mfy(is:ie, js:&
&                         je+1, k), mfy_tl=mfy_tl(is:ie, js:je+1, k))
      call fv_tp_2d(qn2(isd:ied,jsd:jed,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k))
            END IF
            IF (it .LT. nsplt) THEN
! not last call
              DO j=js,je
                DO i=is,ie
                  qn2_tl(i, j, iq) = ((qn2_tl(i, j, iq)*dp1(i, j, k)+qn2&
&                   (i, j, iq)*dp1_tl(i, j, k)+rarea(i, j)*(fx_tl(i, j)-&
&                   fx_tl(i+1, j)+fy_tl(i, j)-fy_tl(i, j+1)))*dp2(i, j)-&
&                   (qn2(i, j, iq)*dp1(i, j, k)+(fx(i, j)-fx(i+1, j)+(fy&
&                   (i, j)-fy(i, j+1)))*rarea(i, j))*dp2_tl(i, j))/dp2(i&
&                   , j)**2
                  qn2(i, j, iq) = (qn2(i, j, iq)*dp1(i, j, k)+(fx(i, j)-&
&                   fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i&
&                   , j)
                END DO
              END DO
            ELSE
              DO j=js,je
                DO i=is,ie
                  q_tl(i, j, k, iq) = ((qn2_tl(i, j, iq)*dp1(i, j, k)+&
&                   qn2(i, j, iq)*dp1_tl(i, j, k)+rarea(i, j)*(fx_tl(i, &
&                   j)-fx_tl(i+1, j)+fy_tl(i, j)-fy_tl(i, j+1)))*dp2(i, &
&                   j)-(qn2(i, j, iq)*dp1(i, j, k)+(fx(i, j)-fx(i+1, j)+&
&                   (fy(i, j)-fy(i, j+1)))*rarea(i, j))*dp2_tl(i, j))/&
&                   dp2(i, j)**2
                  q(i, j, k, iq) = (qn2(i, j, iq)*dp1(i, j, k)+(fx(i, j)&
&                   -fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(&
&                   i, j)
                END DO
              END DO
            END IF
          ELSE
            IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:&
&                            ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, &
&                            k), cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied&
&                            , js:je+1, k), cy_tl(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fx_tl, fy, fy_tl, xfx(&
&                            is:ie+1, jsd:jed, k), xfx_tl(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), yfx_tl(&
&                            isd:ied, js:je+1, k), gridstruct, bd, ra_x&
&                            , ra_x_tl, ra_y, ra_y_tl, mfx=mfx(is:ie+1, &
&                            js:je, k), mfx_tl=mfx_tl(is:ie+1, js:je, k)&
&                            , mfy=mfy(is:ie, js:je+1, k), mfy_tl=mfy_tl&
&                            (is:ie, js:je+1, k))
            ELSE
              CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:ied&
&                         , jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k), &
&                         cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied, js:je+&
&                         1, k), cy_tl(isd:ied, js:je+1, k), npx, npy, &
&                         hord_pert, fx, fx_tl, fy, fy_tl, xfx(is:ie+1, &
&                         jsd:jed, k), xfx_tl(is:ie+1, jsd:jed, k), yfx(&
&                         isd:ied, js:je+1, k), yfx_tl(isd:ied, js:je+1&
&                         , k), gridstruct, bd, ra_x, ra_x_tl, ra_y, &
&                         ra_y_tl, mfx=mfx(is:ie+1, js:je, k), mfx_tl=&
&                         mfx_tl(is:ie+1, js:je, k), mfy=mfy(is:ie, js:&
&                         je+1, k), mfy_tl=mfy_tl(is:ie, js:je+1, k))
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k))
            END IF
            DO j=js,je
              DO i=is,ie
                q_tl(i, j, k, iq) = ((q_tl(i, j, k, iq)*dp1(i, j, k)+q(i&
&                 , j, k, iq)*dp1_tl(i, j, k)+rarea(i, j)*(fx_tl(i, j)-&
&                 fx_tl(i+1, j)+fy_tl(i, j)-fy_tl(i, j+1)))*dp2(i, j)-(q&
&                 (i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-fx(i+1, j)+(fy(i&
&                 , j)-fy(i, j+1)))*rarea(i, j))*dp2_tl(i, j))/dp2(i, j)&
&                 **2
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
              dp1_tl(i, j, k) = dp2_tl(i, j)
              dp1(i, j, k) = dp2(i, j)
            END DO
          END DO
          CALL TIMING_ON('COMM_TOTAL')
          CALL TIMING_ON('COMM_TRACER')
          CALL MPP_UPDATE_DOMAINS_TLM(qn2, qn2_tl, domain)
          CALL TIMING_OFF('COMM_TRACER')
          CALL TIMING_OFF('COMM_TOTAL')
        END IF
      END DO
    END DO
! time-split loop
! k-loop
    IF (PRESENT(dpa)) dpa = dp2
  END SUBROUTINE TRACER_2D_1L_TLM
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
      call fv_tp_2d(qn2(isd:ied,jsd:jed,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k))
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
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k))
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
!  Differentiation of tracer_2d in forward (tangent) mode:
!   variations   of useful results: q dp1
!   with respect to varying inputs: q dp1 mfx mfy cx cy
  SUBROUTINE TRACER_2D_TLM(q, q_tl, dp1, dp1_tl, mfx, mfx_tl, mfy, &
&   mfy_tl, cx, cx_tl, cy, cy_tl, gridstruct, bd, domain, npx, npy, npz&
&   , nq, hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, hord_pert, &
&   nord_tr_pert, trdm_pert, split_damp_tr, dpa)
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
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dp1_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_tl(bd%is:bd%ie+1, bd%js:bd%je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_tl(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
! DELP after advection
    REAL, OPTIONAL, INTENT(OUT) :: dpa(bd%is:bd%ie, bd%js:bd%je, npz)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Local Arrays
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: dp2_tl(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_tl(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_tl(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_tl(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_tl(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
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
    xfx_tl = 0.0
    yfx_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx,cmax,q_split,ksplt)
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx_tl(i, j, k) = dxa(i-1, j)*dy(i, j)*sin_sg(i-1, j, 3)*&
&             cx_tl(i, j, k)
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sin_sg(i-1, &
&             j, 3)
          ELSE
            xfx_tl(i, j, k) = dxa(i, j)*dy(i, j)*sin_sg(i, j, 1)*cx_tl(i&
&             , j, k)
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1&
&             )
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx_tl(i, j, k) = dya(i, j-1)*dx(i, j)*sin_sg(i, j-1, 4)*&
&             cy_tl(i, j, k)
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sin_sg(i, j-&
&             1, 4)
          ELSE
            yfx_tl(i, j, k) = dya(i, j)*dx(i, j)*sin_sg(i, j, 2)*cy_tl(i&
&             , j, k)
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
            cx_tl(i, j, k) = frac*cx_tl(i, j, k)
            cx(i, j, k) = cx(i, j, k)*frac
            xfx_tl(i, j, k) = frac*xfx_tl(i, j, k)
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            mfx_tl(i, j, k) = frac*mfx_tl(i, j, k)
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            cy_tl(i, j, k) = frac*cy_tl(i, j, k)
            cy(i, j, k) = cy(i, j, k)*frac
            yfx_tl(i, j, k) = frac*yfx_tl(i, j, k)
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy_tl(i, j, k) = frac*mfy_tl(i, j, k)
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END DO
      dp2_tl = 0.0
      ra_x_tl = 0.0
      ra_y_tl = 0.0
      fx_tl = 0.0
      fy_tl = 0.0
    ELSE
      dp2_tl = 0.0
      ra_x_tl = 0.0
      ra_y_tl = 0.0
      fx_tl = 0.0
      fy_tl = 0.0
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
              dp2_tl(i, j) = dp1_tl(i, j, k) + rarea(i, j)*(mfx_tl(i, j&
&               , k)-mfx_tl(i+1, j, k)+mfy_tl(i, j, k)-mfy_tl(i, j+1, k)&
&               )
              dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+(&
&               mfy(i, j, k)-mfy(i, j+1, k)))*rarea(i, j)
            END DO
          END DO
          DO j=jsd,jed
            DO i=is,ie
              ra_x_tl(i, j) = xfx_tl(i, j, k) - xfx_tl(i+1, j, k)
              ra_x(i, j) = area(i, j) + (xfx(i, j, k)-xfx(i+1, j, k))
            END DO
          END DO
          DO j=js,je
            DO i=isd,ied
              ra_y_tl(i, j) = yfx_tl(i, j, k) - yfx_tl(i, j+1, k)
              ra_y(i, j) = area(i, j) + (yfx(i, j, k)-yfx(i, j+1, k))
            END DO
          END DO
          DO iq=1,nq
            IF (it .EQ. 1 .AND. trdm .GT. 1.e-4) THEN
              IF (hord .EQ. hord_pert) THEN
                CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(&
&                              isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd&
&                              :jed, k), cx_tl(is:ie+1, jsd:jed, k), cy(&
&                              isd:ied, js:je+1, k), cy_tl(isd:ied, js:&
&                              je+1, k), npx, npy, hord, fx, fx_tl, fy, &
&                              fy_tl, xfx(is:ie+1, jsd:jed, k), xfx_tl(&
&                              is:ie+1, jsd:jed, k), yfx(isd:ied, js:je+&
&                              1, k), yfx_tl(isd:ied, js:je+1, k), &
&                              gridstruct, bd, ra_x, ra_x_tl, ra_y, &
&                              ra_y_tl, mfx=mfx(is:ie+1, js:je, k), &
&                              mfx_tl=mfx_tl(is:ie+1, js:je, k), mfy=mfy&
&                              (is:ie, js:je+1, k), mfy_tl=mfy_tl(is:ie&
&                              , js:je+1, k), mass=dp1(isd:ied, jsd:jed&
&                              , k), mass_tl=dp1_tl(isd:ied, jsd:jed, k)&
&                              , nord=nord_tr, damp_c=trdm)
              ELSE
                CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:&
&                           ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k&
&                           ), cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied, &
&                           js:je+1, k), cy_tl(isd:ied, js:je+1, k), npx&
&                           , npy, hord_pert, fx, fx_tl, fy, fy_tl, xfx(&
&                           is:ie+1, jsd:jed, k), xfx_tl(is:ie+1, jsd:&
&                           jed, k), yfx(isd:ied, js:je+1, k), yfx_tl(&
&                           isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                           ra_x_tl, ra_y, ra_y_tl, mfx=mfx(is:ie+1, js:&
&                           je, k), mfx_tl=mfx_tl(is:ie+1, js:je, k), &
&                           mfy=mfy(is:ie, js:je+1, k), mfy_tl=mfy_tl(is&
&                           :ie, js:je+1, k), mass=dp1(isd:ied, jsd:jed&
&                           , k), mass_tl=dp1_tl(isd:ied, jsd:jed, k), &
&                           nord=nord_tr_pert, damp_c=trdm_pert)
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k),   &
                    mass=dp1(isd:ied,jsd:jed,k), nord=nord_tr, damp_c=trdm)
              END IF
            ELSE IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:&
&                            ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, &
&                            k), cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied&
&                            , js:je+1, k), cy_tl(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fx_tl, fy, fy_tl, xfx(&
&                            is:ie+1, jsd:jed, k), xfx_tl(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), yfx_tl(&
&                            isd:ied, js:je+1, k), gridstruct, bd, ra_x&
&                            , ra_x_tl, ra_y, ra_y_tl, mfx=mfx(is:ie+1, &
&                            js:je, k), mfx_tl=mfx_tl(is:ie+1, js:je, k)&
&                            , mfy=mfy(is:ie, js:je+1, k), mfy_tl=mfy_tl&
&                            (is:ie, js:je+1, k))
            ELSE
              CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:ied&
&                         , jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k), &
&                         cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied, js:je+&
&                         1, k), cy_tl(isd:ied, js:je+1, k), npx, npy, &
&                         hord_pert, fx, fx_tl, fy, fy_tl, xfx(is:ie+1, &
&                         jsd:jed, k), xfx_tl(is:ie+1, jsd:jed, k), yfx(&
&                         isd:ied, js:je+1, k), yfx_tl(isd:ied, js:je+1&
&                         , k), gridstruct, bd, ra_x, ra_x_tl, ra_y, &
&                         ra_y_tl, mfx=mfx(is:ie+1, js:je, k), mfx_tl=&
&                         mfx_tl(is:ie+1, js:je, k), mfy=mfy(is:ie, js:&
&                         je+1, k), mfy_tl=mfy_tl(is:ie, js:je+1, k))
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k))
            END IF
            DO j=js,je
              DO i=is,ie
                q_tl(i, j, k, iq) = ((q_tl(i, j, k, iq)*dp1(i, j, k)+q(i&
&                 , j, k, iq)*dp1_tl(i, j, k)+rarea(i, j)*(fx_tl(i, j)-&
&                 fx_tl(i+1, j)+fy_tl(i, j)-fy_tl(i, j+1)))*dp2(i, j)-(q&
&                 (i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-fx(i+1, j)+(fy(i&
&                 , j)-fy(i, j+1)))*rarea(i, j))*dp2_tl(i, j))/dp2(i, j)&
&                 **2
                q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-&
&                 fx(i+1, j)+(fy(i, j)-fy(i, j+1)))*rarea(i, j))/dp2(i, &
&                 j)
              END DO
            END DO
          END DO
          IF (it .NE. nsplt) THEN
            DO j=js,je
              DO i=is,ie
                dp1_tl(i, j, k) = dp2_tl(i, j)
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
        CALL START_GROUP_HALO_UPDATE_TLM(q_pack, q, q_tl, domain)
        CALL TIMING_OFF('COMM_TRACER')
        CALL TIMING_OFF('COMM_TOTAL')
      END IF
    END DO
! nsplt
    IF (PRESENT(dpa)) dpa = dp1(bd%is:bd%ie, bd%js:bd%je, 1:npz)
  END SUBROUTINE TRACER_2D_TLM
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
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k),   &
                    mass=dp1(isd:ied,jsd:jed,k), nord=nord_tr, damp_c=trdm)
              END IF
            ELSE IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, &
&                        jsd:jed, k), cy(isd:ied, js:je+1, k), npx, npy&
&                        , hord, fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(&
&                        isd:ied, js:je+1, k), gridstruct, bd, ra_x, &
&                        ra_y, mfx=mfx(is:ie+1, js:je, k), mfy=mfy(is:ie&
&                        , js:je+1, k))
            ELSE
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k))
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
!  Differentiation of tracer_2d_nested in forward (tangent) mode:
!   variations   of useful results: q dp1
!   with respect to varying inputs: q dp1 mfx mfy cx cy
  SUBROUTINE TRACER_2D_NESTED_TLM(q, q_tl, dp1, dp1_tl, mfx, mfx_tl, mfy&
&   , mfy_tl, cx, cx_tl, cy, cy_tl, gridstruct, bd, domain, npx, npy, &
&   npz, nq, hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, k_split&
&   , neststruct, parent_grid, hord_pert, nord_tr_pert, trdm_pert, &
&   split_damp_tr)
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
    REAL, INTENT(INOUT) :: q_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: dp1_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(bd%is:bd%ie+1, bd%js:bd%je, npz)
    REAL, INTENT(INOUT) :: mfx_tl(bd%is:bd%ie+1, bd%js:bd%je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: mfy_tl(bd%is:bd%ie, bd%js:bd%je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL, INTENT(INOUT) :: cx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL, INTENT(INOUT) :: cy_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    TYPE(FV_NEST_TYPE), INTENT(INOUT) :: neststruct
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
! Local Arrays
    REAL :: dp2(bd%is:bd%ie, bd%js:bd%je)
    REAL :: dp2_tl(bd%is:bd%ie, bd%js:bd%je)
    REAL :: fx(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fx_tl(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: fy(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: fy_tl(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_x_tl(bd%is:bd%ie, bd%jsd:bd%jed)
    REAL :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: ra_y_tl(bd%isd:bd%ied, bd%js:bd%je)
    REAL :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: xfx_tl(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    REAL :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, npz)
    REAL :: yfx_tl(bd%isd:bd%ied, bd%js:bd%je+1, npz)
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
    xfx_tl = 0.0
    yfx_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx)
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx_tl(i, j, k) = dxa(i-1, j)*dy(i, j)*sin_sg(i-1, j, 3)*&
&             cx_tl(i, j, k)
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sin_sg(i-1, &
&             j, 3)
          ELSE
            xfx_tl(i, j, k) = dxa(i, j)*dy(i, j)*sin_sg(i, j, 1)*cx_tl(i&
&             , j, k)
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1&
&             )
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx_tl(i, j, k) = dya(i, j-1)*dx(i, j)*sin_sg(i, j-1, 4)*&
&             cy_tl(i, j, k)
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sin_sg(i, j-&
&             1, 4)
          ELSE
            yfx_tl(i, j, k) = dya(i, j)*dx(i, j)*sin_sg(i, j, 2)*cy_tl(i&
&             , j, k)
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
            cx_tl(i, j, k) = frac*cx_tl(i, j, k)
            cx(i, j, k) = cx(i, j, k)*frac
            xfx_tl(i, j, k) = frac*xfx_tl(i, j, k)
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            mfx_tl(i, j, k) = frac*mfx_tl(i, j, k)
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            cy_tl(i, j, k) = frac*cy_tl(i, j, k)
            cy(i, j, k) = cy(i, j, k)*frac
            yfx_tl(i, j, k) = frac*yfx_tl(i, j, k)
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy_tl(i, j, k) = frac*mfy_tl(i, j, k)
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END DO
      dp2_tl = 0.0
      ra_x_tl = 0.0
      ra_y_tl = 0.0
      fx_tl = 0.0
      fy_tl = 0.0
    ELSE
      dp2_tl = 0.0
      ra_x_tl = 0.0
      ra_y_tl = 0.0
      fx_tl = 0.0
      fy_tl = 0.0
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
          CALL NESTED_GRID_BC_APPLY_INTT_TLM(q(isd:ied, jsd:jed, :, iq)&
&                                      , q_tl(isd:ied, jsd:jed, :, iq), &
&                                      0, 0, npx, npy, npz, bd, REAL(&
&                                      neststruct%tracer_nest_timestep) &
&                                      + REAL(nsplt*k_split), REAL(nsplt&
&                                      *k_split), neststruct%q_bc(iq), &
&                                      bctype=neststruct%nestbctype)
        END DO
      END IF
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx,mfy,rarea,nq, &
!$OMP                                  area,xfx,yfx,q,cx,cy,npx,npy,hord,gridstruct,bd,it,nsplt,nord_tr,trdm) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dp2_tl(i, j) = dp1_tl(i, j, k) + rarea(i, j)*(mfx_tl(i, j, k&
&             )-mfx_tl(i+1, j, k)+mfy_tl(i, j, k)-mfy_tl(i, j+1, k))
            dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+(mfy&
&             (i, j, k)-mfy(i, j+1, k)))*rarea(i, j)
          END DO
        END DO
        DO j=jsd,jed
          DO i=is,ie
            ra_x_tl(i, j) = xfx_tl(i, j, k) - xfx_tl(i+1, j, k)
            ra_x(i, j) = area(i, j) + (xfx(i, j, k)-xfx(i+1, j, k))
          END DO
        END DO
        DO j=js,je
          DO i=isd,ied
            ra_y_tl(i, j) = yfx_tl(i, j, k) - yfx_tl(i, j+1, k)
            ra_y(i, j) = area(i, j) + (yfx(i, j, k)-yfx(i, j+1, k))
          END DO
        END DO
        DO iq=1,nq
          IF (it .EQ. 1 .AND. trdm .GT. 1.e-4) THEN
            IF (hord .EQ. hord_pert) THEN
              CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:&
&                            ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, &
&                            k), cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied&
&                            , js:je+1, k), cy_tl(isd:ied, js:je+1, k), &
&                            npx, npy, hord, fx, fx_tl, fy, fy_tl, xfx(&
&                            is:ie+1, jsd:jed, k), xfx_tl(is:ie+1, jsd:&
&                            jed, k), yfx(isd:ied, js:je+1, k), yfx_tl(&
&                            isd:ied, js:je+1, k), gridstruct, bd, ra_x&
&                            , ra_x_tl, ra_y, ra_y_tl, mfx=mfx(is:ie+1, &
&                            js:je, k), mfx_tl=mfx_tl(is:ie+1, js:je, k)&
&                            , mfy=mfy(is:ie, js:je+1, k), mfy_tl=mfy_tl&
&                            (is:ie, js:je+1, k), mass=dp1(isd:ied, jsd:&
&                            jed, k), mass_tl=dp1_tl(isd:ied, jsd:jed, k&
&                            ), nord=nord_tr, damp_c=trdm)
            ELSE
              CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:ied&
&                         , jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k), &
&                         cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied, js:je+&
&                         1, k), cy_tl(isd:ied, js:je+1, k), npx, npy, &
&                         hord_pert, fx, fx_tl, fy, fy_tl, xfx(is:ie+1, &
&                         jsd:jed, k), xfx_tl(is:ie+1, jsd:jed, k), yfx(&
&                         isd:ied, js:je+1, k), yfx_tl(isd:ied, js:je+1&
&                         , k), gridstruct, bd, ra_x, ra_x_tl, ra_y, &
&                         ra_y_tl, mfx=mfx(is:ie+1, js:je, k), mfx_tl=&
&                         mfx_tl(is:ie+1, js:je, k), mfy=mfy(is:ie, js:&
&                         je+1, k), mfy_tl=mfy_tl(is:ie, js:je+1, k), &
&                         mass=dp1(isd:ied, jsd:jed, k), mass_tl=dp1_tl(&
&                         isd:ied, jsd:jed, k), nord=nord_tr_pert, &
&                         damp_c=trdm_pert)
       call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k),   &
                    mass=dp1(isd:ied,jsd:jed,k), nord=nord_tr, damp_c=trdm) 
            END IF
          ELSE IF (hord .EQ. hord_pert) THEN
            CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:&
&                          ied, jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k)&
&                          , cx_tl(is:ie+1, jsd:jed, k), cy(isd:ied, js:&
&                          je+1, k), cy_tl(isd:ied, js:je+1, k), npx, &
&                          npy, hord, fx, fx_tl, fy, fy_tl, xfx(is:ie+1&
&                          , jsd:jed, k), xfx_tl(is:ie+1, jsd:jed, k), &
&                          yfx(isd:ied, js:je+1, k), yfx_tl(isd:ied, js:&
&                          je+1, k), gridstruct, bd, ra_x, ra_x_tl, ra_y&
&                          , ra_y_tl, mfx=mfx(is:ie+1, js:je, k), mfx_tl&
&                          =mfx_tl(is:ie+1, js:je, k), mfy=mfy(is:ie, js&
&                          :je+1, k), mfy_tl=mfy_tl(is:ie, js:je+1, k))
          ELSE
            CALL FV_TP_2D_TLM(q(isd:ied, jsd:jed, k, iq), q_tl(isd:ied, &
&                       jsd:jed, k, iq), cx(is:ie+1, jsd:jed, k), cx_tl(&
&                       is:ie+1, jsd:jed, k), cy(isd:ied, js:je+1, k), &
&                       cy_tl(isd:ied, js:je+1, k), npx, npy, hord_pert&
&                       , fx, fx_tl, fy, fy_tl, xfx(is:ie+1, jsd:jed, k)&
&                       , xfx_tl(is:ie+1, jsd:jed, k), yfx(isd:ied, js:&
&                       je+1, k), yfx_tl(isd:ied, js:je+1, k), &
&                       gridstruct, bd, ra_x, ra_x_tl, ra_y, ra_y_tl, &
&                       mfx=mfx(is:ie+1, js:je, k), mfx_tl=mfx_tl(is:ie+&
&                       1, js:je, k), mfy=mfy(is:ie, js:je+1, k), mfy_tl&
&                       =mfy_tl(is:ie, js:je+1, k))
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                   npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                   gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k)) 
          END IF
          DO j=js,je
            DO i=is,ie
              q_tl(i, j, k, iq) = ((q_tl(i, j, k, iq)*dp1(i, j, k)+q(i, &
&               j, k, iq)*dp1_tl(i, j, k)+rarea(i, j)*(fx_tl(i, j)-fx_tl&
&               (i+1, j)+fy_tl(i, j)-fy_tl(i, j+1)))*dp2(i, j)-(q(i, j, &
&               k, iq)*dp1(i, j, k)+(fx(i, j)-fx(i+1, j)+(fy(i, j)-fy(i&
&               , j+1)))*rarea(i, j))*dp2_tl(i, j))/dp2(i, j)**2
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
        CALL START_GROUP_HALO_UPDATE_TLM(q_pack, q, q_tl, domain)
        CALL TIMING_OFF('COMM_TRACER')
        CALL TIMING_OFF('COMM_TOTAL')
      END IF
!Apply nested-grid BCs
      IF (gridstruct%nested) THEN
        DO iq=1,nq
          CALL NESTED_GRID_BC_APPLY_INTT_TLM(q(isd:ied, jsd:jed, :, iq)&
&                                      , q_tl(isd:ied, jsd:jed, :, iq), &
&                                      0, 0, npx, npy, npz, bd, REAL(&
&                                      neststruct%tracer_nest_timestep)&
&                                      , REAL(nsplt*k_split), neststruct&
&                                      %q_bc(iq), bctype=neststruct%&
&                                      nestbctype)
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
            dp1_tl(i, j, k) = rarea(i, j)*rdt*(xfx_tl(i+1, j, k)-xfx_tl(&
&             i, j, k)+yfx_tl(i, j+1, k)-yfx_tl(i, j, k))
            dp1(i, j, k) = (xfx(i+1, j, k)-xfx(i, j, k)+(yfx(i, j+1, k)-&
&             yfx(i, j, k)))*rarea(i, j)*rdt
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE TRACER_2D_NESTED_TLM
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
          CALL NESTED_GRID_BC_APPLY_INTT(q(isd:ied, jsd:jed, :, iq), 0, &
&                                  0, npx, npy, npz, bd, REAL(neststruct&
&                                  %tracer_nest_timestep) + REAL(nsplt*&
&                                  k_split), REAL(nsplt*k_split), &
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
       call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                    npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                    gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k),   &
                    mass=dp1(isd:ied,jsd:jed,k), nord=nord_tr, damp_c=trdm) 
            END IF
          ELSE IF (hord .EQ. hord_pert) THEN
            CALL FV_TP_2D(q(isd:ied, jsd:jed, k, iq), cx(is:ie+1, jsd&
&                      :jed, k), cy(isd:ied, js:je+1, k), npx, npy, hord&
&                      , fx, fy, xfx(is:ie+1, jsd:jed, k), yfx(isd:ied, &
&                      js:je+1, k), gridstruct, bd, ra_x, ra_y, mfx=mfx(&
&                      is:ie+1, js:je, k), mfy=mfy(is:ie, js:je+1, k))
          ELSE
      call fv_tp_2d(q(isd:ied,jsd:jed,k,iq), cx(is:ie+1,jsd:jed,k), cy(isd:ied,js:je+1,k), &
                   npx, npy, hord, fx, fy, xfx(is:ie+1,jsd:jed,k), yfx(isd:ied,js:je+1,k), &
                   gridstruct, bd, ra_x, ra_y, mfx=mfx(is:ie+1,js:je,k), mfy=mfy(is:ie,js:je+1,k)) 
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
          CALL NESTED_GRID_BC_APPLY_INTT(q(isd:ied, jsd:jed, :, iq), 0, &
&                                  0, npx, npy, npz, bd, REAL(neststruct&
&                                  %tracer_nest_timestep), REAL(nsplt*&
&                                  k_split), neststruct%q_bc(iq), &
&                                  neststruct%nestbctype)
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

end module fv_tracer2d_tlm_mod
