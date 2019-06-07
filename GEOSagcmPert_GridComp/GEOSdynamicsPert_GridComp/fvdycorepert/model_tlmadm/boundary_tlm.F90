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
module boundary_tlm_mod

  use fv_mp_mod,         only: ng, isc,jsc,iec,jec, isd,jsd,ied,jed, is,js,ie,je, is_master
  use constants_mod,     only: grav

  use mpp_domains_mod,    only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,    only: CENTER, CORNER, NORTH, EAST
  use mpp_domains_mod,    only: mpp_global_field, mpp_get_pelist
  use mpp_mod,            only: mpp_error, FATAL, mpp_sum, mpp_sync, mpp_npes, mpp_broadcast, WARNING, mpp_pe

  use fv_mp_mod,          only: mp_bcst
  use fv_arrays_mod,      only: fv_atmos_type, fv_nest_BC_type_3D, fv_grid_bounds_type
  use mpp_mod,            only: mpp_send, mpp_recv
  use fv_timing_mod,      only: timing_on, timing_off
  use mpp_domains_mod, only : nest_domain_type, WEST, SOUTH
  use mpp_domains_mod, only : mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod, only : mpp_get_F2C_index, mpp_update_nest_coarse
  !use mpp_domains_mod, only : mpp_get_domain_shift

  implicit none
  public extrapolation_BC
  public nested_grid_bc, update_coarse_grid
  public fill_nested_grid, nested_grid_BC_apply_intT
  public nested_grid_BC_send, nested_grid_BC_recv, nested_grid_BC_save_proc

  public nested_grid_bc_apply_intt_tlm

  interface nested_grid_BC
     module procedure nested_grid_BC_2d
     module procedure nested_grid_BC_mpp
     module procedure nested_grid_BC_mpp_send
     module procedure nested_grid_BC_2D_mpp
     module procedure nested_grid_BC_3d
  end interface


  interface fill_nested_grid
     module procedure fill_nested_grid_2d
     module procedure fill_nested_grid_3d
  end interface

  interface update_coarse_grid
     module procedure update_coarse_grid_mpp
     module procedure update_coarse_grid_mpp_2d
  end interface

CONTAINS
!Linear extrapolation into halo region
!Not to be confused with extrapolated-in-time nested BCs
  SUBROUTINE EXTRAPOLATION_BC(q, istag, jstag, npx, npy, bd, pd_in, &
&   debug_in)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: istag, jstag, npx, npy
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag), INTENT(&
&   INOUT) :: q
    LOGICAL, INTENT(IN), OPTIONAL :: pd_in, debug_in
    INTEGER :: i, j, istart, iend, jstart, jend
    LOGICAL :: pd, debug
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    INTRINSIC REAL
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (isd .LT. 1) THEN
      istart = 1
    ELSE
      istart = isd
    END IF
    IF (ied .GT. npx - 1) THEN
      iend = npx - 1
    ELSE
      iend = ied
    END IF
    IF (jsd .LT. 1) THEN
      jstart = 1
    ELSE
      jstart = jsd
    END IF
    IF (jed .GT. npy - 1) THEN
      jend = npy - 1
    ELSE
      jend = jed
    END IF
!Positive-definite extrapolation: shift from linear extrapolation to zero-gradient when the extrapolated value turns negative.
    IF (PRESENT(pd_in)) THEN
      pd = pd_in
    ELSE
      pd = .false.
    END IF
    IF (PRESENT(debug_in)) THEN
      debug = debug_in
    ELSE
      debug = .false.
    END IF
    IF (is .EQ. 1) THEN
      IF (pd) THEN
        DO j=jstart,jend+jstag
          DO i=0,isd,-1
            IF (REAL(i) .LE. 1. - q(1, j)/(q(2, j)-q(1, j)+1.e-12) .AND.&
&               q(1, j) .LT. q(2, j)) THEN
              q(i, j) = q(i+1, j)
            ELSE
              q(i, j) = REAL(2-i)*q(1, j) - REAL(1-i)*q(2, j)
            END IF
          END DO
        END DO
      ELSE
        DO j=jstart,jend+jstag
          DO i=0,isd,-1
            q(i, j) = REAL(2-i)*q(1, j) - REAL(1-i)*q(2, j)
          END DO
        END DO
      END IF
    END IF
    IF (js .EQ. 1) THEN
      IF (pd) THEN
        DO j=0,jsd,-1
          DO i=istart,iend+istag
            IF (REAL(j) .LE. 1. - q(i, 1)/(q(i, 2)-q(i, 1)+1.e-12) .AND.&
&               q(i, 1) .LT. q(i, 2)) THEN
              q(i, j) = q(i, j+1)
            ELSE
              q(i, j) = REAL(2-j)*q(i, 1) - REAL(1-j)*q(i, 2)
            END IF
          END DO
        END DO
      ELSE
        DO j=0,jsd,-1
          DO i=istart,iend+istag
            q(i, j) = REAL(2-j)*q(i, 1) - REAL(1-j)*q(i, 2)
          END DO
        END DO
      END IF
    END IF
    IF (ie .EQ. npx - 1) THEN
      IF (pd) THEN
        DO j=jstart,jend+jstag
          DO i=ie+1+istag,ied+istag
            IF (REAL(i) .GE. ie + istag + q(ie+istag, j)/(q(ie+istag-1, &
&               j)-q(ie+istag, j)+1.e-12) .AND. q(ie+istag, j) .LT. q(ie&
&               +istag-1, j)) THEN
              q(i, j) = q(i-1, j)
            ELSE
              q(i, j) = REAL(i-(ie+istag-1))*q(ie+istag, j) + REAL(ie+&
&               istag-i)*q(ie+istag-1, j)
            END IF
          END DO
        END DO
      ELSE
        DO j=jstart,jend+jstag
          DO i=ie+1+istag,ied+istag
            q(i, j) = REAL(i-(ie+istag-1))*q(ie+istag, j) + REAL(ie+&
&             istag-i)*q(ie+istag-1, j)
          END DO
        END DO
      END IF
    END IF
    IF (je .EQ. npy - 1) THEN
      IF (pd) THEN
        DO j=je+1+jstag,jed+jstag
          DO i=istart,iend+istag
            IF (REAL(j) .GE. je + jstag + q(i, je+jstag)/(q(i, je+jstag-&
&               1)-q(i, je+jstag)+1.e-12) .AND. q(i, je+jstag-1) .GT. q(&
&               i, je+jstag)) THEN
              q(i, j) = q(i, j-1)
            ELSE
              q(i, j) = REAL(j-(je+jstag-1))*q(i, je+jstag) + REAL(je+&
&               jstag-j)*q(i, je+jstag-1)
            END IF
          END DO
        END DO
      ELSE
        DO j=je+1+jstag,jed+jstag
          DO i=istart,iend+istag
            q(i, j) = REAL(j-(je+jstag-1))*q(i, je+jstag) + REAL(je+&
&             jstag-j)*q(i, je+jstag-1)
          END DO
        END DO
      END IF
    END IF
!CORNERS: Average of extrapolations
    IF (is .EQ. 1 .AND. js .EQ. 1) THEN
      IF (pd) THEN
        DO j=0,jsd,-1
          DO i=0,isd,-1
            IF (REAL(i) .LE. 1. - q(1, j)/(q(2, j)-q(1, j)+1.e-12) .AND.&
&               q(2, j) .GT. q(1, j)) THEN
              q(i, j) = 0.5*q(i+1, j)
            ELSE
              q(i, j) = 0.5*(REAL(2-i)*q(1, j)-REAL(1-i)*q(2, j))
            END IF
            IF (REAL(j) .LE. 1. - q(i, 1)/(q(i, 2)-q(i, 1)+1.e-12) .AND.&
&               q(i, 2) .GT. q(i, 1)) THEN
              q(i, j) = q(i, j) + 0.5*q(i, j+1)
            ELSE
              q(i, j) = q(i, j) + 0.5*(REAL(2-j)*q(i, 1)-REAL(1-j)*q(i, &
&               2))
            END IF
          END DO
        END DO
      ELSE
        DO j=jsd,0
          DO i=isd,0
            q(i, j) = 0.5*(REAL(2-i)*q(1, j)-REAL(1-i)*q(2, j)) + 0.5*(&
&             REAL(2-j)*q(i, 1)-REAL(1-j)*q(i, 2))
          END DO
        END DO
      END IF
    END IF
    IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
      IF (pd) THEN
        DO j=je+1+jstag,jed+jstag
          DO i=0,isd,-1
            IF (REAL(i) .LE. 1. - q(1, j)/(q(2, j)-q(1, j)+1.e-12) .AND.&
&               q(2, j) .GT. q(1, j)) THEN
              q(i, j) = 0.5*q(i+1, j)
            ELSE
              q(i, j) = 0.5*(REAL(2-i)*q(1, j)-REAL(1-i)*q(2, j))
            END IF
!'Unary plus' removed to appease IBM compiler
!if (real(j) >= je+jstag - q(i,je+jstag)/(q(i,je+jstag-1)-q(i,je+jstag)+1.e-12) .and. &
            IF (REAL(j) .GE. je + jstag - q(i, je+jstag)/(q(i, je+jstag-&
&               1)-q(i, je+jstag)+1.e-12) .AND. q(i, je+jstag-1) .GT. q(&
&               i, je+jstag)) THEN
              q(i, j) = q(i, j) + 0.5*q(i, j-1)
            ELSE
              q(i, j) = q(i, j) + 0.5*(REAL(j-(je+jstag-1))*q(i, je+&
&               jstag)+REAL(je+jstag-j)*q(i, je+jstag-1))
            END IF
          END DO
        END DO
      ELSE
        DO j=je+1+jstag,jed+jstag
          DO i=isd,0
            q(i, j) = 0.5*(REAL(2-i)*q(1, j)-REAL(1-i)*q(2, j)) + 0.5*(&
&             REAL(j-(je+jstag-1))*q(i, je+jstag)+REAL(je+jstag-j)*q(i, &
&             je+jstag-1))
          END DO
        END DO
      END IF
    END IF
    IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
      IF (pd) THEN
        DO j=je+1+jstag,jed+jstag
          DO i=ie+1+istag,ied+istag
            IF (REAL(i) .GE. ie + istag + q(ie+istag, j)/(q(ie+istag-1, &
&               j)-q(ie+istag, j)+1.e-12) .AND. q(ie+istag-1, j) .GT. q(&
&               ie+istag, j)) THEN
              q(i, j) = 0.5*q(i-1, j)
            ELSE
              q(i, j) = 0.5*(REAL(i-(ie+istag-1))*q(ie+istag, j)+REAL(ie&
&               +istag-i)*q(ie+istag-1, j))
            END IF
            IF (REAL(j) .GE. je + jstag + q(i, je+jstag)/(q(i, je+jstag-&
&               1)-q(i, je+jstag)+1.e-12) .AND. q(i, je+jstag-1) .GT. q(&
&               i, je+jstag)) THEN
              q(i, j) = q(i, j) + 0.5*q(i, j-1)
            ELSE
              q(i, j) = q(i, j) + 0.5*(REAL(j-(je+jstag-1))*q(i, je+&
&               jstag)+REAL(je+jstag-j)*q(i, je+jstag-1))
            END IF
          END DO
        END DO
      ELSE
        DO j=je+1+jstag,jed+jstag
          DO i=ie+1+istag,ied+istag
            q(i, j) = 0.5*(REAL(i-(ie+istag-1))*q(ie+istag, j)+REAL(ie+&
&             istag-i)*q(ie+istag-1, j)) + 0.5*(REAL(j-(je+jstag-1))*q(i&
&             , je+jstag)+REAL(je+jstag-j)*q(i, je+jstag-1))
          END DO
        END DO
      END IF
    END IF
    IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
      IF (pd) THEN
        DO j=0,jsd,-1
          DO i=ie+1+istag,ied+istag
            IF (REAL(i) .GE. ie + istag + q(ie+istag, j)/(q(ie+istag-1, &
&               j)-q(ie+istag, j)+1.e-12) .AND. q(ie+istag-1, j) .GT. q(&
&               ie+istag, j)) THEN
              q(i, j) = 0.5*q(i-1, j)
            ELSE
              q(i, j) = 0.5*(REAL(i-(ie+istag-1))*q(ie+istag, j)+REAL(ie&
&               +istag-i)*q(ie+istag-1, j))
            END IF
            IF (REAL(j) .LE. 1. - q(i, 1)/(q(i, 2)-q(i, 1)+1.e-12) .AND.&
&               q(i, 2) .GT. q(i, 1)) THEN
              q(i, j) = q(i, j) + 0.5*q(i, j+1)
            ELSE
              q(i, j) = q(i, j) + 0.5*(REAL(2-j)*q(i, 1)-REAL(1-j)*q(i, &
&               2))
            END IF
          END DO
        END DO
      ELSE
        DO j=jsd,0
          DO i=ie+1+istag,ied+istag
            q(i, j) = 0.5*(REAL(i-(ie+istag-1))*q(ie+istag, j)+REAL(ie+&
&             istag-i)*q(ie+istag-1, j)) + 0.5*(REAL(2-j)*q(i, 1)-REAL(1&
&             -j)*q(i, 2))
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE EXTRAPOLATION_BC
  SUBROUTINE FILL_NESTED_GRID_2D(var_nest, var_coarse, ind, wt, istag, &
&   jstag, isg, ieg, jsg, jeg, bd, istart_in, iend_in, jstart_in, &
&   jend_in)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: istag, jstag, isg, ieg, jsg, jeg
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag), INTENT(&
&   INOUT) :: var_nest
    REAL, DIMENSION(isg:ieg+istag, jsg:jeg+jstag), INTENT(IN) :: &
&   var_coarse
    INTEGER, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 2), &
&   INTENT(IN) :: ind
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 4), &
&   INTENT(IN) :: wt
    INTEGER, INTENT(IN), OPTIONAL :: istart_in, iend_in, jstart_in, &
&   jend_in
    INTEGER :: i, j, ic, jc
    INTEGER :: istart, iend, jstart, jend
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC PRESENT
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (PRESENT(istart_in)) THEN
      istart = istart_in
    ELSE
      istart = isd
    END IF
    IF (PRESENT(iend_in)) THEN
      iend = iend_in + istag
    ELSE
      iend = ied + istag
    END IF
    IF (PRESENT(jstart_in)) THEN
      jstart = jstart_in
    ELSE
      jstart = jsd
    END IF
    IF (PRESENT(jend_in)) THEN
      jend = jend_in + jstag
    ELSE
      jend = jed + jstag
    END IF
    DO j=jstart,jend
      DO i=istart,iend
        ic = ind(i, j, 1)
        jc = ind(i, j, 2)
        var_nest(i, j) = wt(i, j, 1)*var_coarse(ic, jc) + wt(i, j, 2)*&
&         var_coarse(ic, jc+1) + wt(i, j, 3)*var_coarse(ic+1, jc+1) + wt&
&         (i, j, 4)*var_coarse(ic+1, jc)
      END DO
    END DO
  END SUBROUTINE FILL_NESTED_GRID_2D
  SUBROUTINE FILL_NESTED_GRID_3D(var_nest, var_coarse, ind, wt, istag, &
&   jstag, isg, ieg, jsg, jeg, npz, bd, istart_in, iend_in, jstart_in, &
&   jend_in)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: istag, jstag, isg, ieg, jsg, jeg, npz
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, npz), &
&   INTENT(INOUT) :: var_nest
    REAL, DIMENSION(isg:ieg+istag, jsg:jeg+jstag, npz), INTENT(IN) :: &
&   var_coarse
    INTEGER, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 2), &
&   INTENT(IN) :: ind
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 4), &
&   INTENT(IN) :: wt
    INTEGER, INTENT(IN), OPTIONAL :: istart_in, iend_in, jstart_in, &
&   jend_in
    INTEGER :: i, j, ic, jc, k
    INTEGER :: istart, iend, jstart, jend
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC PRESENT
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (PRESENT(istart_in)) THEN
      istart = istart_in
    ELSE
      istart = isd
    END IF
    IF (PRESENT(iend_in)) THEN
      iend = iend_in + istag
    ELSE
      iend = ied + istag
    END IF
    IF (PRESENT(jstart_in)) THEN
      jstart = jstart_in
    ELSE
      jstart = jsd
    END IF
    IF (PRESENT(jend_in)) THEN
      jend = jend_in + jstag
    ELSE
      jend = jed + jstag
    END IF
    DO k=1,npz
      DO j=jstart,jend
        DO i=istart,iend
          ic = ind(i, j, 1)
          jc = ind(i, j, 2)
          var_nest(i, j, k) = wt(i, j, 1)*var_coarse(ic, jc, k) + wt(i, &
&           j, 2)*var_coarse(ic, jc+1, k) + wt(i, j, 3)*var_coarse(ic+1&
&           , jc+1, k) + wt(i, j, 4)*var_coarse(ic+1, jc, k)
        END DO
      END DO
    END DO
  END SUBROUTINE FILL_NESTED_GRID_3D
  SUBROUTINE NESTED_GRID_BC_MPP(var_nest, var_coarse, nest_domain, ind, &
&   wt, istag, jstag, npx, npy, npz, bd, isg, ieg, jsg, jeg, nstep_in, &
&   nsplit_in, proc_in)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: istag, jstag, npx, npy, npz, isg, ieg, jsg, &
&   jeg
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, npz), &
&   INTENT(INOUT) :: var_nest
    REAL, DIMENSION(isg:ieg+istag, jsg:jeg+jstag, npz), INTENT(IN) :: &
&   var_coarse
    TYPE(NEST_DOMAIN_TYPE), INTENT(INOUT) :: nest_domain
    INTEGER, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 2), &
&   INTENT(IN) :: ind
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 4), &
&   INTENT(IN) :: wt
    INTEGER, INTENT(IN), OPTIONAL :: nstep_in, nsplit_in
    LOGICAL, INTENT(IN), OPTIONAL :: proc_in
    INTEGER :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
    INTEGER :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
    INTEGER :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
    INTEGER :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c
    REAL, ALLOCATABLE :: wbuffer(:, :, :)
    REAL, ALLOCATABLE :: ebuffer(:, :, :)
    REAL, ALLOCATABLE :: sbuffer(:, :, :)
    REAL, ALLOCATABLE :: nbuffer(:, :, :)
    INTEGER :: i, j, ic, jc, istart, iend, k
    INTEGER :: position
    LOGICAL :: process
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC PRESENT
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (PRESENT(proc_in)) THEN
      process = proc_in
    ELSE
      process = .true.
    END IF
    IF (istag .EQ. 1 .AND. jstag .EQ. 1) THEN
      position = corner
    ELSE IF (istag .EQ. 0 .AND. jstag .EQ. 1) THEN
      position = north
    ELSE IF (istag .EQ. 1 .AND. jstag .EQ. 0) THEN
      position = east
    ELSE
      position = center
    END IF
    CALL MPP_GET_C2F_INDEX(nest_domain, isw_f, iew_f, jsw_f, jew_f, &
&                    isw_c, iew_c, jsw_c, jew_c, west, position)
    CALL MPP_GET_C2F_INDEX(nest_domain, ise_f, iee_f, jse_f, jee_f, &
&                    ise_c, iee_c, jse_c, jee_c, east, position)
    CALL MPP_GET_C2F_INDEX(nest_domain, iss_f, ies_f, jss_f, jes_f, &
&                    iss_c, ies_c, jss_c, jes_c, south, position)
    CALL MPP_GET_C2F_INDEX(nest_domain, isn_f, ien_f, jsn_f, jen_f, &
&                    isn_c, ien_c, jsn_c, jen_c, north, position)
    IF (iew_c .GE. isw_c .AND. jew_c .GE. jsw_c) THEN
      ALLOCATE(wbuffer(isw_c:iew_c, jsw_c:jew_c, npz))
    ELSE
      ALLOCATE(wbuffer(1, 1, 1))
    END IF
    wbuffer = 0.0
    IF (iee_c .GE. ise_c .AND. jee_c .GE. jse_c) THEN
      ALLOCATE(ebuffer(ise_c:iee_c, jse_c:jee_c, npz))
    ELSE
      ALLOCATE(ebuffer(1, 1, 1))
    END IF
    ebuffer = 0.0
    IF (ies_c .GE. iss_c .AND. jes_c .GE. jss_c) THEN
      ALLOCATE(sbuffer(iss_c:ies_c, jss_c:jes_c, npz))
    ELSE
      ALLOCATE(sbuffer(1, 1, 1))
    END IF
    sbuffer = 0.0
    IF (ien_c .GE. isn_c .AND. jen_c .GE. jsn_c) THEN
      ALLOCATE(nbuffer(isn_c:ien_c, jsn_c:jen_c, npz))
    ELSE
      ALLOCATE(nbuffer(1, 1, 1))
    END IF
    nbuffer = 0.0
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_NEST_FINE(var_coarse, nest_domain, wbuffer, sbuffer&
&                       , ebuffer, nbuffer, position)
    CALL TIMING_OFF('COMM_TOTAL')
!process
    IF (process) THEN
      IF (is .EQ. 1) THEN
        DO k=1,npz
          DO j=jsd,jed+jstag
            DO i=isd,0
              ic = ind(i, j, 1)
              jc = ind(i, j, 2)
              var_nest(i, j, k) = wt(i, j, 1)*wbuffer(ic, jc, k) + wt(i&
&               , j, 2)*wbuffer(ic, jc+1, k) + wt(i, j, 3)*wbuffer(ic+1&
&               , jc+1, k) + wt(i, j, 4)*wbuffer(ic+1, jc, k)
            END DO
          END DO
        END DO
      END IF
      IF (js .EQ. 1) THEN
        IF (is .EQ. 1) THEN
          istart = is
        ELSE
          istart = isd
        END IF
        IF (ie .EQ. npx - 1) THEN
          iend = ie
        ELSE
          iend = ied
        END IF
        DO k=1,npz
          DO j=jsd,0
            DO i=istart,iend+istag
              ic = ind(i, j, 1)
              jc = ind(i, j, 2)
              var_nest(i, j, k) = wt(i, j, 1)*sbuffer(ic, jc, k) + wt(i&
&               , j, 2)*sbuffer(ic, jc+1, k) + wt(i, j, 3)*sbuffer(ic+1&
&               , jc+1, k) + wt(i, j, 4)*sbuffer(ic+1, jc, k)
            END DO
          END DO
        END DO
      END IF
      IF (ie .EQ. npx - 1) THEN
        DO k=1,npz
          DO j=jsd,jed+jstag
            DO i=npx+istag,ied+istag
              ic = ind(i, j, 1)
              jc = ind(i, j, 2)
              var_nest(i, j, k) = wt(i, j, 1)*ebuffer(ic, jc, k) + wt(i&
&               , j, 2)*ebuffer(ic, jc+1, k) + wt(i, j, 3)*ebuffer(ic+1&
&               , jc+1, k) + wt(i, j, 4)*ebuffer(ic+1, jc, k)
            END DO
          END DO
        END DO
      END IF
      IF (je .EQ. npy - 1) THEN
        IF (is .EQ. 1) THEN
          istart = is
        ELSE
          istart = isd
        END IF
        IF (ie .EQ. npx - 1) THEN
          iend = ie
        ELSE
          iend = ied
        END IF
        DO k=1,npz
          DO j=npy+jstag,jed+jstag
            DO i=istart,iend+istag
              ic = ind(i, j, 1)
              jc = ind(i, j, 2)
              var_nest(i, j, k) = wt(i, j, 1)*nbuffer(ic, jc, k) + wt(i&
&               , j, 2)*nbuffer(ic, jc+1, k) + wt(i, j, 3)*nbuffer(ic+1&
&               , jc+1, k) + wt(i, j, 4)*nbuffer(ic+1, jc, k)
            END DO
          END DO
        END DO
      END IF
    END IF
    DEALLOCATE(wbuffer)
    DEALLOCATE(ebuffer)
    DEALLOCATE(sbuffer)
    DEALLOCATE(nbuffer)
  END SUBROUTINE NESTED_GRID_BC_MPP
  SUBROUTINE NESTED_GRID_BC_MPP_SEND(var_coarse, nest_domain, istag, &
&   jstag)
    IMPLICIT NONE
    REAL, DIMENSION(:, :, :), INTENT(IN) :: var_coarse
    TYPE(NEST_DOMAIN_TYPE), INTENT(INOUT) :: nest_domain
    INTEGER, INTENT(IN) :: istag, jstag
    REAL, ALLOCATABLE :: wbuffer(:, :, :)
    REAL, ALLOCATABLE :: ebuffer(:, :, :)
    REAL, ALLOCATABLE :: sbuffer(:, :, :)
    REAL, ALLOCATABLE :: nbuffer(:, :, :)
    INTEGER :: i, j, ic, jc, istart, iend, k
    INTEGER :: position
    IF (istag .EQ. 1 .AND. jstag .EQ. 1) THEN
      position = corner
    ELSE IF (istag .EQ. 0 .AND. jstag .EQ. 1) THEN
      position = north
    ELSE IF (istag .EQ. 1 .AND. jstag .EQ. 0) THEN
      position = east
    ELSE
      position = center
    END IF
    ALLOCATE(wbuffer(1, 1, 1))
    ALLOCATE(ebuffer(1, 1, 1))
    ALLOCATE(sbuffer(1, 1, 1))
    ALLOCATE(nbuffer(1, 1, 1))
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_NEST_FINE(var_coarse, nest_domain, wbuffer, sbuffer&
&                       , ebuffer, nbuffer, position)
    CALL TIMING_OFF('COMM_TOTAL')
    DEALLOCATE(wbuffer)
    DEALLOCATE(ebuffer)
    DEALLOCATE(sbuffer)
    DEALLOCATE(nbuffer)
  END SUBROUTINE NESTED_GRID_BC_MPP_SEND
  SUBROUTINE NESTED_GRID_BC_2D_MPP(var_nest, var_coarse, nest_domain, &
&   ind, wt, istag, jstag, npx, npy, bd, isg, ieg, jsg, jeg, nstep_in, &
&   nsplit_in, proc_in)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag), INTENT(&
&   INOUT) :: var_nest
    REAL, DIMENSION(isg:ieg+istag, jsg:jeg+jstag), INTENT(IN) :: &
&   var_coarse
    TYPE(NEST_DOMAIN_TYPE), INTENT(INOUT) :: nest_domain
    INTEGER, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 2), &
&   INTENT(IN) :: ind
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 4), &
&   INTENT(IN) :: wt
    INTEGER, INTENT(IN), OPTIONAL :: nstep_in, nsplit_in
    LOGICAL, INTENT(IN), OPTIONAL :: proc_in
    INTEGER :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
    INTEGER :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
    INTEGER :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
    INTEGER :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c
    REAL, ALLOCATABLE :: wbuffer(:, :)
    REAL, ALLOCATABLE :: ebuffer(:, :)
    REAL, ALLOCATABLE :: sbuffer(:, :)
    REAL, ALLOCATABLE :: nbuffer(:, :)
    INTEGER :: i, j, ic, jc, istart, iend, k
    INTEGER :: position
    LOGICAL :: process
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC PRESENT
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (PRESENT(proc_in)) THEN
      process = proc_in
    ELSE
      process = .true.
    END IF
    IF (istag .EQ. 1 .AND. jstag .EQ. 1) THEN
      position = corner
    ELSE IF (istag .EQ. 0 .AND. jstag .EQ. 1) THEN
      position = north
    ELSE IF (istag .EQ. 1 .AND. jstag .EQ. 0) THEN
      position = east
    ELSE
      position = center
    END IF
    CALL MPP_GET_C2F_INDEX(nest_domain, isw_f, iew_f, jsw_f, jew_f, &
&                    isw_c, iew_c, jsw_c, jew_c, west, position)
    CALL MPP_GET_C2F_INDEX(nest_domain, ise_f, iee_f, jse_f, jee_f, &
&                    ise_c, iee_c, jse_c, jee_c, east, position)
    CALL MPP_GET_C2F_INDEX(nest_domain, iss_f, ies_f, jss_f, jes_f, &
&                    iss_c, ies_c, jss_c, jes_c, south, position)
    CALL MPP_GET_C2F_INDEX(nest_domain, isn_f, ien_f, jsn_f, jen_f, &
&                    isn_c, ien_c, jsn_c, jen_c, north, position)
    IF (iew_c .GE. isw_c .AND. jew_c .GE. jsw_c) THEN
      ALLOCATE(wbuffer(isw_c:iew_c, jsw_c:jew_c))
    ELSE
      ALLOCATE(wbuffer(1, 1))
    END IF
    wbuffer = 0.0
    IF (iee_c .GE. ise_c .AND. jee_c .GE. jse_c) THEN
      ALLOCATE(ebuffer(ise_c:iee_c, jse_c:jee_c))
    ELSE
      ALLOCATE(ebuffer(1, 1))
    END IF
    ebuffer = 0.0
    IF (ies_c .GE. iss_c .AND. jes_c .GE. jss_c) THEN
      ALLOCATE(sbuffer(iss_c:ies_c, jss_c:jes_c))
    ELSE
      ALLOCATE(sbuffer(1, 1))
    END IF
    sbuffer = 0.0
    IF (ien_c .GE. isn_c .AND. jen_c .GE. jsn_c) THEN
      ALLOCATE(nbuffer(isn_c:ien_c, jsn_c:jen_c))
    ELSE
      ALLOCATE(nbuffer(1, 1))
    END IF
    nbuffer = 0.0
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_NEST_FINE(var_coarse, nest_domain, wbuffer, sbuffer&
&                       , ebuffer, nbuffer, position)
    CALL TIMING_OFF('COMM_TOTAL')
!process
    IF (process) THEN
      IF (is .EQ. 1) THEN
        DO j=jsd,jed+jstag
          DO i=isd,0
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_nest(i, j) = wt(i, j, 1)*wbuffer(ic, jc) + wt(i, j, 2)*&
&             wbuffer(ic, jc+1) + wt(i, j, 3)*wbuffer(ic+1, jc+1) + wt(i&
&             , j, 4)*wbuffer(ic+1, jc)
          END DO
        END DO
      END IF
      IF (js .EQ. 1) THEN
        IF (is .EQ. 1) THEN
          istart = is
        ELSE
          istart = isd
        END IF
        IF (ie .EQ. npx - 1) THEN
          iend = ie
        ELSE
          iend = ied
        END IF
        DO j=jsd,0
          DO i=istart,iend+istag
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_nest(i, j) = wt(i, j, 1)*sbuffer(ic, jc) + wt(i, j, 2)*&
&             sbuffer(ic, jc+1) + wt(i, j, 3)*sbuffer(ic+1, jc+1) + wt(i&
&             , j, 4)*sbuffer(ic+1, jc)
          END DO
        END DO
      END IF
      IF (ie .EQ. npx - 1) THEN
        DO j=jsd,jed+jstag
          DO i=npx+istag,ied+istag
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_nest(i, j) = wt(i, j, 1)*ebuffer(ic, jc) + wt(i, j, 2)*&
&             ebuffer(ic, jc+1) + wt(i, j, 3)*ebuffer(ic+1, jc+1) + wt(i&
&             , j, 4)*ebuffer(ic+1, jc)
          END DO
        END DO
      END IF
      IF (je .EQ. npy - 1) THEN
        IF (is .EQ. 1) THEN
          istart = is
        ELSE
          istart = isd
        END IF
        IF (ie .EQ. npx - 1) THEN
          iend = ie
        ELSE
          iend = ied
        END IF
        DO j=npy+jstag,jed+jstag
          DO i=istart,iend+istag
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_nest(i, j) = wt(i, j, 1)*nbuffer(ic, jc) + wt(i, j, 2)*&
&             nbuffer(ic, jc+1) + wt(i, j, 3)*nbuffer(ic+1, jc+1) + wt(i&
&             , j, 4)*nbuffer(ic+1, jc)
          END DO
        END DO
      END IF
    END IF
    DEALLOCATE(wbuffer)
    DEALLOCATE(ebuffer)
    DEALLOCATE(sbuffer)
    DEALLOCATE(nbuffer)
  END SUBROUTINE NESTED_GRID_BC_2D_MPP
  SUBROUTINE NESTED_GRID_BC_2D(var_nest, var_coarse, ind, wt, istag, &
&   jstag, npx, npy, bd, isg, ieg, jsg, jeg, nstep_in, nsplit_in)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag), INTENT(&
&   INOUT) :: var_nest
    REAL, DIMENSION(isg:ieg+istag, jsg:jeg+jstag), INTENT(IN) :: &
&   var_coarse
    INTEGER, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 2), &
&   INTENT(IN) :: ind
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 4), &
&   INTENT(IN) :: wt
    INTEGER, INTENT(IN), OPTIONAL :: nstep_in, nsplit_in
    INTEGER :: nstep, nsplit
    INTEGER :: i, j, ic, jc, istart, iend
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC PRESENT
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF ((.NOT.PRESENT(nstep_in)) .OR. (.NOT.PRESENT(nsplit_in))) THEN
      nstep = 1
      nsplit = 2
    ELSE
      nstep = nstep_in
      nsplit = nsplit_in
    END IF
    IF (is .EQ. 1) THEN
      DO j=jsd,jed+jstag
        DO i=isd,0
          ic = ind(i, j, 1)
          jc = ind(i, j, 2)
          var_nest(i, j) = wt(i, j, 1)*var_coarse(ic, jc) + wt(i, j, 2)*&
&           var_coarse(ic, jc+1) + wt(i, j, 3)*var_coarse(ic+1, jc+1) + &
&           wt(i, j, 4)*var_coarse(ic+1, jc)
        END DO
      END DO
    END IF
    IF (js .EQ. 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
      DO j=jsd,0
        DO i=istart,iend+istag
          ic = ind(i, j, 1)
          jc = ind(i, j, 2)
          var_nest(i, j) = wt(i, j, 1)*var_coarse(ic, jc) + wt(i, j, 2)*&
&           var_coarse(ic, jc+1) + wt(i, j, 3)*var_coarse(ic+1, jc+1) + &
&           wt(i, j, 4)*var_coarse(ic+1, jc)
        END DO
      END DO
    END IF
    IF (ie .EQ. npx - 1) THEN
      DO j=jsd,jed+jstag
        DO i=npx+istag,ied+istag
          ic = ind(i, j, 1)
          jc = ind(i, j, 2)
          var_nest(i, j) = wt(i, j, 1)*var_coarse(ic, jc) + wt(i, j, 2)*&
&           var_coarse(ic, jc+1) + wt(i, j, 3)*var_coarse(ic+1, jc+1) + &
&           wt(i, j, 4)*var_coarse(ic+1, jc)
        END DO
      END DO
    END IF
    IF (je .EQ. npy - 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
      DO j=npy+jstag,jed+jstag
        DO i=istart,iend+istag
          ic = ind(i, j, 1)
          jc = ind(i, j, 2)
          var_nest(i, j) = wt(i, j, 1)*var_coarse(ic, jc) + wt(i, j, 2)*&
&           var_coarse(ic, jc+1) + wt(i, j, 3)*var_coarse(ic+1, jc+1) + &
&           wt(i, j, 4)*var_coarse(ic+1, jc)
        END DO
      END DO
    END IF
  END SUBROUTINE NESTED_GRID_BC_2D
  SUBROUTINE NESTED_GRID_BC_3D(var_nest, var_coarse, ind, wt, istag, &
&   jstag, npx, npy, npz, bd, isg, ieg, jsg, jeg, nstep_in, nsplit_in)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg, &
&   npz
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, npz), &
&   INTENT(INOUT) :: var_nest
    REAL, DIMENSION(isg:ieg+istag, jsg:jeg+jstag, npz), INTENT(IN) :: &
&   var_coarse
    INTEGER, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 2), &
&   INTENT(IN) :: ind
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 4), &
&   INTENT(IN) :: wt
    INTEGER, INTENT(IN), OPTIONAL :: nstep_in, nsplit_in
    INTEGER :: nstep, nsplit
    INTEGER :: i, j, ic, jc, istart, iend, k
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC PRESENT
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF ((.NOT.PRESENT(nstep_in)) .OR. (.NOT.PRESENT(nsplit_in))) THEN
      nstep = 1
      nsplit = 2
    ELSE
      nstep = nstep_in
      nsplit = nsplit_in
    END IF
    IF (is .EQ. 1) THEN
      DO k=1,npz
        DO j=jsd,jed+jstag
          DO i=isd,0
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_nest(i, j, k) = wt(i, j, 1)*var_coarse(ic, jc, k) + wt(i&
&             , j, 2)*var_coarse(ic, jc+1, k) + wt(i, j, 3)*var_coarse(&
&             ic+1, jc+1, k) + wt(i, j, 4)*var_coarse(ic+1, jc, k)
          END DO
        END DO
      END DO
    END IF
    IF (js .EQ. 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
      DO k=1,npz
        DO j=jsd,0
          DO i=istart,iend+istag
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_nest(i, j, k) = wt(i, j, 1)*var_coarse(ic, jc, k) + wt(i&
&             , j, 2)*var_coarse(ic, jc+1, k) + wt(i, j, 3)*var_coarse(&
&             ic+1, jc+1, k) + wt(i, j, 4)*var_coarse(ic+1, jc, k)
          END DO
        END DO
      END DO
    END IF
    IF (ie .EQ. npx - 1) THEN
      DO k=1,npz
        DO j=jsd,jed+jstag
          DO i=npx+istag,ied+istag
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_nest(i, j, k) = wt(i, j, 1)*var_coarse(ic, jc, k) + wt(i&
&             , j, 2)*var_coarse(ic, jc+1, k) + wt(i, j, 3)*var_coarse(&
&             ic+1, jc+1, k) + wt(i, j, 4)*var_coarse(ic+1, jc, k)
          END DO
        END DO
      END DO
    END IF
    IF (je .EQ. npy - 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
      DO k=1,npz
        DO j=npy+jstag,jed+jstag
          DO i=istart,iend+istag
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_nest(i, j, k) = wt(i, j, 1)*var_coarse(ic, jc, k) + wt(i&
&             , j, 2)*var_coarse(ic, jc+1, k) + wt(i, j, 3)*var_coarse(&
&             ic+1, jc+1, k) + wt(i, j, 4)*var_coarse(ic+1, jc, k)
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE NESTED_GRID_BC_3D
  SUBROUTINE NESTED_GRID_BC_SEND(var_coarse, nest_domain, istag, jstag)
    IMPLICIT NONE
    REAL, DIMENSION(:, :, :), INTENT(IN) :: var_coarse
    TYPE(NEST_DOMAIN_TYPE), INTENT(INOUT) :: nest_domain
    INTEGER, INTENT(IN) :: istag, jstag
    INTEGER :: position
    REAL :: wbuffer(1, 1, 1)
    REAL :: ebuffer(1, 1, 1)
    REAL :: sbuffer(1, 1, 1)
    REAL :: nbuffer(1, 1, 1)
    IF (istag .EQ. 1 .AND. jstag .EQ. 1) THEN
      position = corner
    ELSE IF (istag .EQ. 0 .AND. jstag .EQ. 1) THEN
      position = north
    ELSE IF (istag .EQ. 1 .AND. jstag .EQ. 0) THEN
      position = east
    ELSE
      position = center
    END IF
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_NEST_FINE(var_coarse, nest_domain, wbuffer, sbuffer&
&                       , ebuffer, nbuffer, position)
    CALL TIMING_OFF('COMM_TOTAL')
  END SUBROUTINE NESTED_GRID_BC_SEND
  SUBROUTINE NESTED_GRID_BC_RECV(nest_domain, istag, jstag, npz, bd, &
&   nest_bc_buffers)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(NEST_DOMAIN_TYPE), INTENT(INOUT) :: nest_domain
    INTEGER, INTENT(IN) :: istag, jstag, npz
    TYPE(FV_NEST_BC_TYPE_3D), INTENT(INOUT), TARGET :: nest_bc_buffers
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, npz) :: &
&   var_coarse_dummy
    INTEGER :: position
    INTEGER :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
    INTEGER :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
    INTEGER :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
    INTEGER :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c
    INTEGER :: i, j, k
    INTRINSIC ALLOCATED
    var_coarse_dummy = 0.0
    IF (istag .EQ. 1 .AND. jstag .EQ. 1) THEN
      position = corner
    ELSE IF (istag .EQ. 0 .AND. jstag .EQ. 1) THEN
      position = north
    ELSE IF (istag .EQ. 1 .AND. jstag .EQ. 0) THEN
      position = east
    ELSE
      position = center
    END IF
    IF (.NOT.ALLOCATED(nest_bc_buffers%west_t1)) THEN
      CALL MPP_GET_C2F_INDEX(nest_domain, isw_f, iew_f, jsw_f, jew_f, &
&                      isw_c, iew_c, jsw_c, jew_c, west, position)
      CALL MPP_GET_C2F_INDEX(nest_domain, ise_f, iee_f, jse_f, jee_f, &
&                      ise_c, iee_c, jse_c, jee_c, east, position)
      CALL MPP_GET_C2F_INDEX(nest_domain, iss_f, ies_f, jss_f, jes_f, &
&                      iss_c, ies_c, jss_c, jes_c, south, position)
      CALL MPP_GET_C2F_INDEX(nest_domain, isn_f, ien_f, jsn_f, jen_f, &
&                      isn_c, ien_c, jsn_c, jen_c, north, position)
      IF (iew_c .GE. isw_c .AND. jew_c .GE. jsw_c) THEN
        IF (.NOT.ALLOCATED(nest_bc_buffers%west_t1)) THEN
          ALLOCATE(nest_bc_buffers%west_t1(isw_c:iew_c, jsw_c:jew_c, npz&
&         ))
        END IF
!compatible with first touch principle
        DO k=1,npz
          DO j=jsw_c,jew_c
            DO i=isw_c,iew_c
              nest_bc_buffers%west_t1(i, j, k) = 0.
            END DO
          END DO
        END DO
      ELSE
        ALLOCATE(nest_bc_buffers%west_t1(1, 1, 1))
        nest_bc_buffers%west_t1(1, 1, 1) = 0.
      END IF
      IF (iee_c .GE. ise_c .AND. jee_c .GE. jse_c) THEN
        IF (.NOT.ALLOCATED(nest_bc_buffers%east_t1)) THEN
          ALLOCATE(nest_bc_buffers%east_t1(ise_c:iee_c, jse_c:jee_c, npz&
&         ))
        END IF
        DO k=1,npz
          DO j=jse_c,jee_c
            DO i=ise_c,iee_c
              nest_bc_buffers%east_t1(i, j, k) = 0.
            END DO
          END DO
        END DO
      ELSE
        ALLOCATE(nest_bc_buffers%east_t1(1, 1, 1))
        nest_bc_buffers%east_t1(1, 1, 1) = 0.
      END IF
      IF (ies_c .GE. iss_c .AND. jes_c .GE. jss_c) THEN
        IF (.NOT.ALLOCATED(nest_bc_buffers%south_t1)) THEN
          ALLOCATE(nest_bc_buffers%south_t1(iss_c:ies_c, jss_c:jes_c, &
&         npz))
        END IF
        DO k=1,npz
          DO j=jss_c,jes_c
            DO i=iss_c,ies_c
              nest_bc_buffers%south_t1(i, j, k) = 0.
            END DO
          END DO
        END DO
      ELSE
        ALLOCATE(nest_bc_buffers%south_t1(1, 1, 1))
        nest_bc_buffers%south_t1(1, 1, 1) = 0.
      END IF
      IF (ien_c .GE. isn_c .AND. jen_c .GE. jsn_c) THEN
        IF (.NOT.ALLOCATED(nest_bc_buffers%north_t1)) THEN
          ALLOCATE(nest_bc_buffers%north_t1(isn_c:ien_c, jsn_c:jen_c, &
&         npz))
        END IF
        DO k=1,npz
          DO j=jsn_c,jen_c
            DO i=isn_c,ien_c
              nest_bc_buffers%north_t1(i, j, k) = 0.
            END DO
          END DO
        END DO
      ELSE
        ALLOCATE(nest_bc_buffers%north_t1(1, 1, 1))
        nest_bc_buffers%north_t1(1, 1, 1) = 0
      END IF
    END IF
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_NEST_FINE(var_coarse_dummy, nest_domain, &
&                       nest_bc_buffers%west_t1, nest_bc_buffers%&
&                       south_t1, nest_bc_buffers%east_t1, &
&                       nest_bc_buffers%north_t1, position)
    CALL TIMING_OFF('COMM_TOTAL')
  END SUBROUTINE NESTED_GRID_BC_RECV
  SUBROUTINE NESTED_GRID_BC_SAVE_PROC(nest_domain, ind, wt, istag, jstag&
&   , npx, npy, npz, bd, nest_bc, nest_bc_buffers, pd_in)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    TYPE(NEST_DOMAIN_TYPE), INTENT(INOUT) :: nest_domain
    INTEGER, INTENT(IN) :: istag, jstag, npx, npy, npz
    INTEGER, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 2), &
&   INTENT(IN) :: ind
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, 4), &
&   INTENT(IN) :: wt
    LOGICAL, INTENT(IN), OPTIONAL :: pd_in
!!NOTE: if declaring an ALLOCATABLE array with intent(OUT), the resulting dummy array
!!      will NOT be allocated! This goes for allocatable members of derived types as well.
    TYPE(FV_NEST_BC_TYPE_3D), INTENT(INOUT), TARGET :: nest_bc, &
&   nest_bc_buffers
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, npz) :: &
&   var_coarse_dummy
    REAL, DIMENSION(:, :, :), POINTER :: var_east, var_west, var_south, &
&   var_north
    REAL, DIMENSION(:, :, :), POINTER :: buf_east, buf_west, buf_south, &
&   buf_north
    INTEGER :: position
    INTEGER :: i, j, k, ic, jc, istart, iend
    LOGICAL :: process
    LOGICAL, SAVE :: pd=.false.
    INTEGER :: is, ie, js, je
    INTEGER :: isd, ied, jsd, jed
    INTRINSIC PRESENT
    INTRINSIC MAX
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    IF (PRESENT(pd_in)) THEN
      pd = pd_in
    ELSE
      pd = .false.
    END IF
    var_east => nest_bc%east_t1
    var_west => nest_bc%west_t1
    var_north => nest_bc%north_t1
    var_south => nest_bc%south_t1
    buf_east => nest_bc_buffers%east_t1
    buf_west => nest_bc_buffers%west_t1
    buf_north => nest_bc_buffers%north_t1
    buf_south => nest_bc_buffers%south_t1
! ?buffer has uninterpolated coarse-grid data; need to perform interpolation ourselves
!To do this more securely, instead of using is/etc we could use the fine-grid indices defined above
    IF (is .EQ. 1) THEN
!$NO-MP parallel do default(none) shared(npz,isd,ied,jsd,jed,jstag,ind,var_west,wt,buf_west) private(ic,jc)
      DO k=1,npz
        DO j=jsd,jed+jstag
          DO i=isd,0
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_west(i, j, k) = wt(i, j, 1)*buf_west(ic, jc, k) + wt(i, &
&             j, 2)*buf_west(ic, jc+1, k) + wt(i, j, 3)*buf_west(ic+1, &
&             jc+1, k) + wt(i, j, 4)*buf_west(ic+1, jc, k)
          END DO
        END DO
      END DO
      IF (pd) THEN
!$NO-MP parallel do default(none) shared(npz,jsd,jed,jstag,isd,var_west,nest_BC)
        DO k=1,npz
          DO j=jsd,jed+jstag
            DO i=isd,0
              IF (var_west(i, j, k) .LT. 0.5*nest_bc%west_t0(i, j, k)) &
&             THEN
                var_west(i, j, k) = 0.5*nest_bc%west_t0(i, j, k)
              ELSE
                var_west(i, j, k) = var_west(i, j, k)
              END IF
            END DO
          END DO
        END DO
      END IF
    END IF
    IF (js .EQ. 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
!$NO-MP parallel do default(none) shared(npz,istart,iend,jsd,jed,istag,ind,var_south,wt,buf_south) private(ic,jc)
      DO k=1,npz
        DO j=jsd,0
          DO i=istart,iend+istag
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_south(i, j, k) = wt(i, j, 1)*buf_south(ic, jc, k) + wt(i&
&             , j, 2)*buf_south(ic, jc+1, k) + wt(i, j, 3)*buf_south(ic+&
&             1, jc+1, k) + wt(i, j, 4)*buf_south(ic+1, jc, k)
          END DO
        END DO
      END DO
      IF (pd) THEN
!$NO-MP parallel do default(none) shared(npz,jsd,jed,istart,iend,istag,var_south,nest_BC)
        DO k=1,npz
          DO j=jsd,0
            DO i=istart,iend+istag
              IF (var_south(i, j, k) .LT. 0.5*nest_bc%south_t0(i, j, k)&
&             ) THEN
                var_south(i, j, k) = 0.5*nest_bc%south_t0(i, j, k)
              ELSE
                var_south(i, j, k) = var_south(i, j, k)
              END IF
            END DO
          END DO
        END DO
      END IF
    END IF
    IF (ie .EQ. npx - 1) THEN
!$NO-MP parallel do default(none) shared(npx,npz,isd,ied,jsd,jed,istag,jstag,ind,var_east,wt,buf_east) private(ic,jc)
      DO k=1,npz
        DO j=jsd,jed+jstag
          DO i=npx+istag,ied+istag
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_east(i, j, k) = wt(i, j, 1)*buf_east(ic, jc, k) + wt(i, &
&             j, 2)*buf_east(ic, jc+1, k) + wt(i, j, 3)*buf_east(ic+1, &
&             jc+1, k) + wt(i, j, 4)*buf_east(ic+1, jc, k)
          END DO
        END DO
      END DO
      IF (pd) THEN
!$NO-MP parallel do default(none) shared(npx,npz,jsd,jed,istag,jstag,ied,var_east,nest_BC)
        DO k=1,npz
          DO j=jsd,jed+jstag
            DO i=npx+istag,ied+istag
              IF (var_east(i, j, k) .LT. 0.5*nest_bc%east_t0(i, j, k)) &
&             THEN
                var_east(i, j, k) = 0.5*nest_bc%east_t0(i, j, k)
              ELSE
                var_east(i, j, k) = var_east(i, j, k)
              END IF
            END DO
          END DO
        END DO
      END IF
    END IF
    IF (je .EQ. npy - 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
!$NO-MP parallel do default(none) shared(npy,npz,istart,iend,jsd,jed,istag,jstag,ind,var_north,wt,buf_north) private(ic,jc)
      DO k=1,npz
        DO j=npy+jstag,jed+jstag
          DO i=istart,iend+istag
            ic = ind(i, j, 1)
            jc = ind(i, j, 2)
            var_north(i, j, k) = wt(i, j, 1)*buf_north(ic, jc, k) + wt(i&
&             , j, 2)*buf_north(ic, jc+1, k) + wt(i, j, 3)*buf_north(ic+&
&             1, jc+1, k) + wt(i, j, 4)*buf_north(ic+1, jc, k)
          END DO
        END DO
      END DO
      IF (pd) THEN
!$NO-MP parallel do default(none) shared(npy,npz,jsd,jed,istart,iend,istag,jstag,ied,var_north,nest_BC)
        DO k=1,npz
          DO j=npy+jstag,jed+jstag
            DO i=istart,iend+istag
              IF (var_north(i, j, k) .LT. 0.5*nest_bc%north_t0(i, j, k)&
&             ) THEN
                var_north(i, j, k) = 0.5*nest_bc%north_t0(i, j, k)
              ELSE
                var_north(i, j, k) = var_north(i, j, k)
              END IF
            END DO
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE NESTED_GRID_BC_SAVE_PROC
!  Differentiation of nested_grid_bc_apply_intt in forward (tangent) mode:
!   variations   of useful results: var_nest
!   with respect to varying inputs: var_nest
! A NOTE ON BCTYPE: currently only an interpolation BC is implemented,
! bctype >= 2 currently correspond
! to a flux BC on the tracers ONLY, which is implemented in fv_tracer.
  SUBROUTINE NESTED_GRID_BC_APPLY_INTT_TLM(var_nest, var_nest_tl, istag&
&   , jstag, npx, npy, npz, bd, step, split, bc, bctype)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: istag, jstag, npx, npy, npz
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, npz), &
&   INTENT(INOUT) :: var_nest
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, npz), &
&   INTENT(INOUT) :: var_nest_tl
    REAL, INTENT(IN) :: split, step
    INTEGER, INTENT(IN) :: bctype
    TYPE(FV_NEST_BC_TYPE_3D), INTENT(IN), TARGET :: bc
    REAL, DIMENSION(:, :, :), POINTER :: var_t0, var_t1
    INTEGER :: i, j, istart, iend, k
    REAL :: denom
    LOGICAL, SAVE :: printdiag=.true.
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
    denom = 1./split
    IF (is .EQ. 1) THEN
      var_t0 => bc%west_t0
      var_t1 => bc%west_t1
      DO k=1,npz
        DO j=jsd,jed+jstag
          DO i=isd,0
            var_nest_tl(i, j, k) = 0.0
            var_nest(i, j, k) = (var_t0(i, j, k)*(split-step)+step*&
&             var_t1(i, j, k))*denom
          END DO
        END DO
      END DO
    END IF
    IF (js .EQ. 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
      var_t0 => bc%south_t0
      var_t1 => bc%south_t1
      DO k=1,npz
        DO j=jsd,0
          DO i=istart,iend+istag
            var_nest_tl(i, j, k) = 0.0
            var_nest(i, j, k) = (var_t0(i, j, k)*(split-step)+step*&
&             var_t1(i, j, k))*denom
          END DO
        END DO
      END DO
    END IF
    IF (ie .EQ. npx - 1) THEN
      var_t0 => bc%east_t0
      var_t1 => bc%east_t1
      DO k=1,npz
        DO j=jsd,jed+jstag
          DO i=npx+istag,ied+istag
            var_nest_tl(i, j, k) = 0.0
            var_nest(i, j, k) = (var_t0(i, j, k)*(split-step)+step*&
&             var_t1(i, j, k))*denom
          END DO
        END DO
      END DO
    END IF
    IF (je .EQ. npy - 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
      var_t0 => bc%north_t0
      var_t1 => bc%north_t1
      DO k=1,npz
        DO j=npy+jstag,jed+jstag
          DO i=istart,iend+istag
            var_nest_tl(i, j, k) = 0.0
            var_nest(i, j, k) = (var_t0(i, j, k)*(split-step)+step*&
&             var_t1(i, j, k))*denom
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE NESTED_GRID_BC_APPLY_INTT_TLM
! A NOTE ON BCTYPE: currently only an interpolation BC is implemented,
! bctype >= 2 currently correspond
! to a flux BC on the tracers ONLY, which is implemented in fv_tracer.
  SUBROUTINE NESTED_GRID_BC_APPLY_INTT(var_nest, istag, jstag, npx, npy&
&   , npz, bd, step, split, bc, bctype)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: istag, jstag, npx, npy, npz
    REAL, DIMENSION(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, npz), &
&   INTENT(INOUT) :: var_nest
    REAL, INTENT(IN) :: split, step
    INTEGER, INTENT(IN) :: bctype
    TYPE(FV_NEST_BC_TYPE_3D), INTENT(IN), TARGET :: bc
    REAL, DIMENSION(:, :, :), POINTER :: var_t0, var_t1
    INTEGER :: i, j, istart, iend, k
    REAL :: denom
    LOGICAL, SAVE :: printdiag=.true.
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
    denom = 1./split
    IF (is .EQ. 1) THEN
      var_t0 => bc%west_t0
      var_t1 => bc%west_t1
      DO k=1,npz
        DO j=jsd,jed+jstag
          DO i=isd,0
            var_nest(i, j, k) = (var_t0(i, j, k)*(split-step)+step*&
&             var_t1(i, j, k))*denom
          END DO
        END DO
      END DO
    END IF
    IF (js .EQ. 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
      var_t0 => bc%south_t0
      var_t1 => bc%south_t1
      DO k=1,npz
        DO j=jsd,0
          DO i=istart,iend+istag
            var_nest(i, j, k) = (var_t0(i, j, k)*(split-step)+step*&
&             var_t1(i, j, k))*denom
          END DO
        END DO
      END DO
    END IF
    IF (ie .EQ. npx - 1) THEN
      var_t0 => bc%east_t0
      var_t1 => bc%east_t1
      DO k=1,npz
        DO j=jsd,jed+jstag
          DO i=npx+istag,ied+istag
            var_nest(i, j, k) = (var_t0(i, j, k)*(split-step)+step*&
&             var_t1(i, j, k))*denom
          END DO
        END DO
      END DO
    END IF
    IF (je .EQ. npy - 1) THEN
      IF (is .EQ. 1) THEN
        istart = is
      ELSE
        istart = isd
      END IF
      IF (ie .EQ. npx - 1) THEN
        iend = ie
      ELSE
        iend = ied
      END IF
      var_t0 => bc%north_t0
      var_t1 => bc%north_t1
      DO k=1,npz
        DO j=npy+jstag,jed+jstag
          DO i=istart,iend+istag
            var_nest(i, j, k) = (var_t0(i, j, k)*(split-step)+step*&
&             var_t1(i, j, k))*denom
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE NESTED_GRID_BC_APPLY_INTT
  SUBROUTINE UPDATE_COARSE_GRID_MPP_2D(var_coarse, var_nest, nest_domain&
&   , ind_update, dx, dy, area, isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, &
&   js_n, je_n, isu, ieu, jsu, jeu, npx, npy, istag, jstag, r, &
&   nestupdate, upoff, nsponge, parent_proc, child_proc, parent_grid)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n&
&   , je_n
    INTEGER, INTENT(IN) :: isu, ieu, jsu, jeu
    INTEGER, INTENT(IN) :: istag, jstag, r, nestupdate, upoff, nsponge
    INTEGER, INTENT(IN) :: ind_update(isd_p:ied_p+1, jsd_p:jed_p+1, 2)
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: var_nest(is_n:ie_n+istag, js_n:je_n+jstag)
    REAL, INTENT(INOUT) :: var_coarse(isd_p:ied_p+istag, jsd_p:jed_p+&
&   jstag)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: area(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: parent_proc, child_proc
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(NEST_DOMAIN_TYPE), INTENT(INOUT) :: nest_domain
    REAL :: var_nest_3d(is_n:ie_n+istag, js_n:je_n+jstag, 1)
    REAL :: var_coarse_3d(isd_p:ied_p+istag, jsd_p:jed_p+jstag, 1)
    INTRINSIC SIZE
    IF (child_proc .AND. SIZE(var_nest) .GT. 1) var_nest_3d(is_n:ie_n+&
&     istag, js_n:je_n+jstag, 1) = var_nest(is_n:ie_n+istag, js_n:je_n+&
&       jstag)
    IF (parent_proc .AND. SIZE(var_coarse) .GT. 1) var_coarse_3d(isd_p:&
&     ied_p+istag, jsd_p:jed_p, 1) = var_coarse(isd_p:ied_p+istag, jsd_p&
&       :jed_p+jstag)
    CALL UPDATE_COARSE_GRID_MPP(var_coarse_3d, var_nest_3d, nest_domain&
&                         , ind_update, dx, dy, area, isd_p, ied_p, &
&                         jsd_p, jed_p, is_n, ie_n, js_n, je_n, isu, ieu&
&                         , jsu, jeu, npx, npy, 1, istag, jstag, r, &
&                         nestupdate, upoff, nsponge, parent_proc, &
&                         child_proc, parent_grid)
    IF (SIZE(var_coarse) .GT. 1 .AND. parent_proc) var_coarse(isd_p:&
&     ied_p+istag, jsd_p:jed_p+jstag) = var_coarse_3d(isd_p:ied_p+istag&
&       , jsd_p:jed_p, 1)
  END SUBROUTINE UPDATE_COARSE_GRID_MPP_2D
  SUBROUTINE UPDATE_COARSE_GRID_MPP(var_coarse, var_nest, nest_domain, &
&   ind_update, dx, dy, area, isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, &
&   js_n, je_n, isu, ieu, jsu, jeu, npx, npy, npz, istag, jstag, r, &
&   nestupdate, upoff, nsponge, parent_proc, child_proc, parent_grid)
    IMPLICIT NONE
!This routine assumes the coarse and nested grids are properly
! aligned, and that in particular for odd refinement ratios all
! coarse-grid points coincide with nested-grid points
    INTEGER, INTENT(IN) :: isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n&
&   , je_n
    INTEGER, INTENT(IN) :: isu, ieu, jsu, jeu
    INTEGER, INTENT(IN) :: istag, jstag, npx, npy, npz, r, nestupdate, &
&   upoff, nsponge
    INTEGER, INTENT(IN) :: ind_update(isd_p:ied_p+1, jsd_p:jed_p+1, 2)
    REAL, INTENT(IN) :: var_nest(is_n:ie_n+istag, js_n:je_n+jstag, npz)
    REAL, INTENT(INOUT) :: var_coarse(isd_p:ied_p+istag, jsd_p:jed_p+&
&   jstag, npz)
    REAL, INTENT(IN) :: area(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    LOGICAL, INTENT(IN) :: parent_proc, child_proc
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT) :: parent_grid
    TYPE(NEST_DOMAIN_TYPE), INTENT(INOUT) :: nest_domain
    INTEGER :: in, jn, ini, jnj, s, qr
    INTEGER :: is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f
    INTEGER :: istart, istop, jstart, jstop, ishift, jshift, j, i, k
    REAL :: val
    REAL, DIMENSION(:, :, :), ALLOCATABLE :: nest_dat
    REAL :: var_nest_send(is_n:ie_n+istag, js_n:je_n+jstag, npz)
    INTEGER :: position
    IF (istag .EQ. 1 .AND. jstag .EQ. 1) THEN
      position = corner
    ELSE IF (istag .EQ. 0 .AND. jstag .EQ. 1) THEN
      position = north
    ELSE IF (istag .EQ. 1 .AND. jstag .EQ. 0) THEN
      position = east
    ELSE
      position = center
    END IF
    CALL MPP_GET_F2C_INDEX(nest_domain, is_c, ie_c, js_c, je_c, is_f, &
&                    ie_f, js_f, je_f, position)
    IF (ie_f .GT. is_f .AND. je_f .GT. js_f) THEN
      ALLOCATE(nest_dat(is_f:ie_f, js_f:je_f, npz))
    ELSE
      ALLOCATE(nest_dat(1, 1, 1))
    END IF
    nest_dat = -600.0
    IF (child_proc) THEN
!! IF an area average (for istag == jstag == 0) or a linear average then multiply in the areas before sending data
      IF (istag .EQ. 0 .AND. jstag .EQ. 0) THEN
        SELECT CASE  (nestupdate) 
        CASE (1, 2, 6, 7, 8) 
!$NO-MP parallel do default(none) shared(npz,js_n,je_n,is_n,ie_n,var_nest_send,var_nest,area)
          DO k=1,npz
            DO j=js_n,je_n
              DO i=is_n,ie_n
                var_nest_send(i, j, k) = var_nest(i, j, k)*area(i, j)
              END DO
            END DO
          END DO
        END SELECT
      ELSE IF (istag .EQ. 0 .AND. jstag .GT. 0) THEN
        SELECT CASE  (nestupdate) 
        CASE (1, 6, 7, 8) 
!$NO-MP parallel do default(none) shared(npz,js_n,je_n,is_n,ie_n,var_nest_send,var_nest,dx)
          DO k=1,npz
            DO j=js_n,je_n+1
              DO i=is_n,ie_n
                var_nest_send(i, j, k) = var_nest(i, j, k)*dx(i, j)
              END DO
            END DO
          END DO
        CASE DEFAULT
          CALL MPP_ERROR(fatal, 'nestupdate type not implemented')
        END SELECT
      ELSE IF (istag .GT. 0 .AND. jstag .EQ. 0) THEN
        SELECT CASE  (nestupdate) 
        CASE (1, 6, 7, 8) 
!averaging update; in-line average for face-averaged values instead of areal average
!$NO-MP parallel do default(none) shared(npz,js_n,je_n,is_n,ie_n,var_nest_send,var_nest,dy)
          DO k=1,npz
            DO j=js_n,je_n
              DO i=is_n,ie_n+1
                var_nest_send(i, j, k) = var_nest(i, j, k)*dy(i, j)
              END DO
            END DO
          END DO
        CASE DEFAULT
          CALL MPP_ERROR(fatal, 'nestupdate type not implemented')
        END SELECT
      ELSE
        CALL MPP_ERROR(fatal, &
&                'Cannot have both nonzero istag and jstag.')
      END IF
    END IF
    CALL TIMING_ON('COMM_TOTAL')
    CALL MPP_UPDATE_NEST_COARSE(var_nest_send, nest_domain, nest_dat, &
&                         position=position)
    CALL TIMING_OFF('COMM_TOTAL')
!rounds down (since r > 0)
    s = r/2
    qr = r*upoff + nsponge - s
    IF (parent_proc .AND. (.NOT.(ieu .LT. isu .OR. jeu .LT. jsu))) THEN
      IF (istag .EQ. 0 .AND. jstag .EQ. 0) THEN
        SELECT CASE  (nestupdate) 
        CASE (1, 2, 6, 7, 8) 
! 1 = Conserving update on all variables; 2 = conserving update for cell-centered values; 6 = conserving remap-update
!$NO-MP parallel do default(none) shared(npz,jsu,jeu,isu,ieu,ind_update,nest_dat,parent_grid,var_coarse,r) &
!$NO-MP          private(in,jn,val)
          DO k=1,npz
            DO j=jsu,jeu
              DO i=isu,ieu
                in = ind_update(i, j, 1)
                jn = ind_update(i, j, 2)
!!$            if (in < max(1+qr,is_f) .or. in > min(npx-1-qr-r+1,ie_f) .or. &
!!$                 jn < max(1+qr,js_f) .or. jn > min(npy-1-qr-r+1,je_f)) then
!!$               write(mpp_pe()+3000,'(A, 14I6)') 'SKIP: ', i, j, in, jn, 1+qr, is_f, ie_f, js_f, je_f, npy-1-qr-r+1, isu, ieu, 
!jsu, jeu
!!$               cycle
!!$            endif
                val = 0.
                DO jnj=jn,jn+r-1
                  DO ini=in,in+r-1
                    val = val + nest_dat(ini, jnj, k)
                  END DO
                END DO
!var_coarse(i,j,k) = val/r**2.
!!! CLEANUP: Couldn't rarea and rdx and rdy be built into the weight arrays?
!!!    Two-way updates do not yet have weights, tho
                var_coarse(i, j, k) = val*parent_grid%gridstruct%rarea(i&
&                 , j)
              END DO
            END DO
          END DO
        CASE DEFAULT
          CALL MPP_ERROR(fatal, 'nestupdate type not implemented')
        END SELECT
      ELSE IF (istag .EQ. 0 .AND. jstag .GT. 0) THEN
        SELECT CASE  (nestupdate) 
        CASE (1, 6, 7, 8) 
!$NO-MP parallel do default(none) shared(npz,jsu,jeu,isu,ieu,ind_update,nest_dat,parent_grid,var_coarse,r) &
!$NO-MP          private(in,jn,val)
          DO k=1,npz
            DO j=jsu,jeu+1
              DO i=isu,ieu
                in = ind_update(i, j, 1)
                jn = ind_update(i, j, 2)
!!$            if (in < max(1+qr,is_f) .or. in > min(npx-1-qr-r+1,ie_f) .or. &
!!$                 jn < max(1+qr+s,js_f) .or. jn > min(npy-1-qr-s+1,je_f)) then
!!$               write(mpp_pe()+3000,'(A, 14I)') 'SKIP u: ', i, j, in, jn, 1+qr, is_f, ie_f, js_f, je_f, npy-1-qr-s+1, isu, ieu,
! jsu, jeu
!!$               cycle
!!$            endif
                val = 0.
                DO ini=in,in+r-1
                  val = val + nest_dat(ini, jn, k)
                END DO
!            var_coarse(i,j,k) = val/r
                var_coarse(i, j, k) = val*parent_grid%gridstruct%rdx(i, &
&                 j)
              END DO
            END DO
          END DO
        CASE DEFAULT
          CALL MPP_ERROR(fatal, 'nestupdate type not implemented')
        END SELECT
      ELSE IF (istag .GT. 0 .AND. jstag .EQ. 0) THEN
        SELECT CASE  (nestupdate) 
        CASE (1, 6, 7, 8) 
!averaging update; in-line average for face-averaged values instead of areal average
!$NO-MP parallel do default(none) shared(npz,jsu,jeu,isu,ieu,ind_update,nest_dat,parent_grid,var_coarse,r) &
!$NO-MP          private(in,jn,val)
          DO k=1,npz
            DO j=jsu,jeu
              DO i=isu,ieu+1
                in = ind_update(i, j, 1)
                jn = ind_update(i, j, 2)
!!$            if (in < max(1+qr+s,is_f) .or. in > min(npx-1-qr-s+1,ie_f) .or. &
!!$                 jn < max(1+qr,js_f) .or. jn > min(npy-1-qr-r+1,je_f)) then
!!$               write(mpp_pe()+3000,'(A, 14I6)') 'SKIP v: ', i, j, in, jn, 1+qr, is_f, ie_f, js_f, je_f, npx-1-qr-s+1, isu, ieu
!, jsu, jeu
!!$               cycle
!!$            endif
                val = 0.
                DO jnj=jn,jn+r-1
                  val = val + nest_dat(in, jnj, k)
                END DO
!            var_coarse(i,j,k) = val/r
                var_coarse(i, j, k) = val*parent_grid%gridstruct%rdy(i, &
&                 j)
              END DO
            END DO
          END DO
        CASE DEFAULT
          CALL MPP_ERROR(fatal, 'nestupdate type not implemented')
        END SELECT
      END IF
    END IF
    DEALLOCATE(nest_dat)
  END SUBROUTINE UPDATE_COARSE_GRID_MPP
   
end module boundary_tlm_mod
