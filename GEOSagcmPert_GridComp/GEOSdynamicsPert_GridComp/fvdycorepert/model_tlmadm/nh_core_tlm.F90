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
module nh_core_tlm_mod
! Developer: S.-J. Lin, NOAA/GFDL
! To do list:
! include moisture effect in pt
!------------------------------
   use constants_mod,     only: rdgas, cp_air, grav
   use tp_core_tlm_mod,   only: fv_tp_2d
   use tp_core_tlm_mod,   only: fv_tp_2d_tlm
   use nh_utils_tlm_mod,  only: update_dz_c, update_dz_d, nest_halo_nh
   use nh_utils_tlm_mod,  only: sim_solver, sim1_solver, sim3_solver
   use nh_utils_tlm_mod,  only: sim3p0_solver, rim_2d
   use nh_utils_tlm_mod,  only: Riem_Solver_c
   use nh_utils_tlm_mod,  only: update_dz_c_tlm, update_dz_d_tlm, nest_halo_nh_tlm
   use nh_utils_tlm_mod,  only: sim_solver_tlm, sim1_solver_tlm, sim3_solver_tlm
   use nh_utils_tlm_mod,  only: sim3p0_solver_tlm, rim_2d_tlm
   use nh_utils_tlm_mod,  only: Riem_Solver_c_tlm

   implicit none
   private

   public Riem_Solver3, Riem_Solver_c, update_dz_c, update_dz_d, nest_halo_nh
   public Riem_Solver3_tlm, Riem_Solver_c_tlm, update_dz_c_tlm, update_dz_d_tlm, nest_halo_nh_tlm
   real, parameter:: r3 = 1./3.

CONTAINS
!  Differentiation of riem_solver3 in forward (tangent) mode:
!   variations   of useful results: pk3 ppe peln w delz pe pk zh
!   with respect to varying inputs: pk3 ws ppe peln w delp delz
!                pe pk zh pt
  SUBROUTINE RIEM_SOLVER3_TLM(ms, dt, is, ie, js, je, km, ng, isd, ied, &
&   jsd, jed, akap, cappa, cp, ptop, zs, q_con, w, w_tl, delz, delz_tl, &
&   pt, pt_tl, delp, delp_tl, zh, zh_tl, pe, pe_tl, ppe, ppe_tl, pk3, &
&   pk3_tl, pk, pk_tl, peln, peln_tl, ws, ws_tl, scale_m, p_fac, a_imp, &
&   use_logp, last_call, fp_out)
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
    REAL, INTENT(IN) :: ws_tl(is:ie, js:je)
    REAL, DIMENSION(isd:, jsd:, :), INTENT(IN) :: q_con, cappa
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: delp, pt
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: delp_tl, pt_tl
    REAL, DIMENSION(isd:ied, jsd:jed, km+1), INTENT(INOUT) :: zh
    REAL, DIMENSION(isd:ied, jsd:jed, km+1), INTENT(INOUT) :: zh_tl
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: w
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: w_tl
    REAL, INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(INOUT) :: pe_tl(is-1:ie+1, km+1, js-1:je+1)
! ln(pe)
    REAL, INTENT(OUT) :: peln(is:ie, km+1, js:je)
    REAL, INTENT(OUT) :: peln_tl(is:ie, km+1, js:je)
    REAL, DIMENSION(isd:ied, jsd:jed, km+1), INTENT(OUT) :: ppe
    REAL, DIMENSION(isd:ied, jsd:jed, km+1), INTENT(OUT) :: ppe_tl
    REAL, INTENT(OUT) :: delz(is-ng:ie+ng, js-ng:je+ng, km)
    REAL, INTENT(OUT) :: delz_tl(is-ng:ie+ng, js-ng:je+ng, km)
    REAL, INTENT(OUT) :: pk(is:ie, js:je, km+1)
    REAL, INTENT(OUT) :: pk_tl(is:ie, js:je, km+1)
    REAL, INTENT(OUT) :: pk3(isd:ied, jsd:jed, km+1)
    REAL, INTENT(OUT) :: pk3_tl(isd:ied, jsd:jed, km+1)
! Local:
    REAL, DIMENSION(is:ie, km) :: dm, dz2, pm2, w2, gm2, cp2
    REAL, DIMENSION(is:ie, km) :: dm_tl, dz2_tl, pm2_tl, w2_tl
    REAL, DIMENSION(is:ie, km+1) :: pem, pe2, peln2, peg, pelng
    REAL, DIMENSION(is:ie, km+1) :: pem_tl, pe2_tl, peln2_tl
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
    dm_tl = 0.0
    pe2_tl = 0.0
    dz2_tl = 0.0
    w2_tl = 0.0
    pm2_tl = 0.0
    pem_tl = 0.0
    peln2_tl = 0.0
!$OMP parallel do default(none) shared(is,ie,js,je,km,delp,ptop,peln1,pk3,ptk,akap,rgrav,zh,pt, &
!$OMP                                  w,a_imp,dt,gama,ws,p_fac,scale_m,ms,delz,last_call,  &
!$OMP                                  peln,pk,fp_out,ppe,use_logp,zs,pe,cappa,q_con )          &
!$OMP                          private(cp2, gm2, dm, dz2, pm2, pem, peg, pelng, pe2, peln2, w2)
    DO j=js,je
      DO k=1,km
        DO i=is,ie
          dm_tl(i, k) = delp_tl(i, j, k)
          dm(i, k) = delp(i, j, k)
        END DO
      END DO
      DO i=is,ie
        pem_tl(i, 1) = 0.0
        pem(i, 1) = ptop
        peln2_tl(i, 1) = 0.0
        peln2(i, 1) = peln1
        pk3_tl(i, j, 1) = 0.0
        pk3(i, j, 1) = ptk
      END DO
      DO k=2,km+1
        DO i=is,ie
          pem_tl(i, k) = pem_tl(i, k-1) + dm_tl(i, k-1)
          pem(i, k) = pem(i, k-1) + dm(i, k-1)
          peln2_tl(i, k) = pem_tl(i, k)/pem(i, k)
          peln2(i, k) = LOG(pem(i, k))
          pk3_tl(i, j, k) = akap*peln2_tl(i, k)*EXP(akap*peln2(i, k))
          pk3(i, j, k) = EXP(akap*peln2(i, k))
        END DO
      END DO
      DO k=1,km
        DO i=is,ie
          pm2_tl(i, k) = (dm_tl(i, k)*(peln2(i, k+1)-peln2(i, k))-dm(i, &
&           k)*(peln2_tl(i, k+1)-peln2_tl(i, k)))/(peln2(i, k+1)-peln2(i&
&           , k))**2
          pm2(i, k) = dm(i, k)/(peln2(i, k+1)-peln2(i, k))
          dm_tl(i, k) = rgrav*dm_tl(i, k)
          dm(i, k) = dm(i, k)*rgrav
          dz2_tl(i, k) = zh_tl(i, j, k+1) - zh_tl(i, j, k)
          dz2(i, k) = zh(i, j, k+1) - zh(i, j, k)
          w2_tl(i, k) = w_tl(i, j, k)
          w2(i, k) = w(i, j, k)
        END DO
      END DO
      IF (a_imp .LT. -0.999) THEN
        CALL SIM3P0_SOLVER_TLM(dt, is, ie, km, rdgas, gama, akap, pe2, &
&                        pe2_tl, dm, dm_tl, pem, pem_tl, w2, w2_tl, dz2&
&                        , dz2_tl, pt(is:ie, j, 1:km), pt_tl(is:ie, j, 1&
&                        :km), ws(is:ie, j), ws_tl(is:ie, j), p_fac, &
&                        scale_m)
      ELSE IF (a_imp .LT. -0.5) THEN
        IF (a_imp .GE. 0.) THEN
          abs0 = a_imp
        ELSE
          abs0 = -a_imp
        END IF
        CALL SIM3_SOLVER_TLM(dt, is, ie, km, rdgas, gama, akap, pe2, &
&                      pe2_tl, dm, dm_tl, pem, pem_tl, w2, w2_tl, dz2, &
&                      dz2_tl, pt(is:ie, j, 1:km), pt_tl(is:ie, j, 1:km)&
&                      , ws(is:ie, j), ws_tl(is:ie, j), abs0, p_fac, &
&                      scale_m)
      ELSE IF (a_imp .LE. 0.5) THEN
        CALL RIM_2D_TLM(ms, dt, is, ie, km, rdgas, gama, gm2, pe2, &
&                 pe2_tl, dm, dm_tl, pm2, pm2_tl, w2, w2_tl, dz2, dz2_tl&
&                 , pt(is:ie, j, 1:km), pt_tl(is:ie, j, 1:km), ws(is:ie&
&                 , j), ws_tl(is:ie, j), .false.)
      ELSE IF (a_imp .GT. 0.999) THEN
        CALL SIM1_SOLVER_TLM(dt, is, ie, km, rdgas, gama, gm2, cp2, akap&
&                      , pe2, pe2_tl, dm, dm_tl, pm2, pm2_tl, pem, &
&                      pem_tl, w2, w2_tl, dz2, dz2_tl, pt(is:ie, j, 1:km&
&                      ), pt_tl(is:ie, j, 1:km), ws(is:ie, j), ws_tl(is:&
&                      ie, j), p_fac)
      ELSE
        CALL SIM_SOLVER_TLM(dt, is, ie, km, rdgas, gama, gm2, cp2, akap&
&                     , pe2, pe2_tl, dm, dm_tl, pm2, pm2_tl, pem, pem_tl&
&                     , w2, w2_tl, dz2, dz2_tl, pt(is:ie, j, 1:km), &
&                     pt_tl(is:ie, j, 1:km), ws(is:ie, j), ws_tl(is:ie, &
&                     j), a_imp, p_fac, scale_m)
      END IF
      DO k=1,km
        DO i=is,ie
          w_tl(i, j, k) = w2_tl(i, k)
          w(i, j, k) = w2(i, k)
          delz_tl(i, j, k) = dz2_tl(i, k)
          delz(i, j, k) = dz2(i, k)
        END DO
      END DO
      IF (last_call) THEN
        DO k=1,km+1
          DO i=is,ie
            peln_tl(i, k, j) = peln2_tl(i, k)
            peln(i, k, j) = peln2(i, k)
            pk_tl(i, j, k) = pk3_tl(i, j, k)
            pk(i, j, k) = pk3(i, j, k)
            pe_tl(i, k, j) = pem_tl(i, k)
            pe(i, k, j) = pem(i, k)
          END DO
        END DO
      END IF
      IF (fp_out) THEN
        DO k=1,km+1
          DO i=is,ie
            ppe_tl(i, j, k) = pe2_tl(i, k) + pem_tl(i, k)
            ppe(i, j, k) = pe2(i, k) + pem(i, k)
          END DO
        END DO
      ELSE
        DO k=1,km+1
          DO i=is,ie
            ppe_tl(i, j, k) = pe2_tl(i, k)
            ppe(i, j, k) = pe2(i, k)
          END DO
        END DO
      END IF
      IF (use_logp) THEN
        DO k=2,km+1
          DO i=is,ie
            pk3_tl(i, j, k) = peln2_tl(i, k)
            pk3(i, j, k) = peln2(i, k)
          END DO
        END DO
      END IF
      DO i=is,ie
        zh_tl(i, j, km+1) = 0.0
        zh(i, j, km+1) = zs(i, j)
      END DO
      DO k=km,1,-1
        DO i=is,ie
          zh_tl(i, j, k) = zh_tl(i, j, k+1) - dz2_tl(i, k)
          zh(i, j, k) = zh(i, j, k+1) - dz2(i, k)
        END DO
      END DO
    END DO
  END SUBROUTINE RIEM_SOLVER3_TLM
  SUBROUTINE RIEM_SOLVER3(ms, dt, is, ie, js, je, km, ng, isd, ied, jsd&
&   , jed, akap, cappa, cp, ptop, zs, q_con, w, delz, pt, delp, zh, pe, &
&   ppe, pk3, pk, peln, ws, scale_m, p_fac, a_imp, use_logp, last_call, &
&   fp_out)
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
    REAL, DIMENSION(isd:, jsd:, :), INTENT(IN) :: q_con, cappa
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
  END SUBROUTINE RIEM_SOLVER3

end module nh_core_tlm_mod
