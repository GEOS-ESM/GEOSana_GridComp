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
module fv_sg_tlm_mod

!-----------------------------------------------------------------------
! FV sub-grid mixing
!-----------------------------------------------------------------------
  use constants_mod,      only: rdgas, rvgas, cp_air, cp_vapor, hlv, hlf, kappa, grav
  use tracer_manager_mod, only: get_tracer_index
  use field_manager_mod,  only: MODEL_ATMOS
  use lin_cld_microphys_mod, only: wqs2, wqsat2_moist
  use fv_mp_mod,          only: mp_reduce_min, is_master

implicit none
private

public  fv_subgrid_z_tlm
public  fv_subgrid_z, neg_adj3

  real, parameter:: esl = 0.621971831
  real, parameter:: tice = 273.16
! real, parameter:: c_ice = 2106.  ! Emanuel table, page 566
  real, parameter:: c_ice = 1972.  !  -15 C
  real, parameter:: c_liq = 4.1855e+3    ! GFS
! real, parameter:: c_liq = 4218.        ! ECMWF-IFS
  real, parameter:: cv_vap = cp_vapor - rvgas  ! 1384.5
  real, parameter:: c_con = c_ice

! real, parameter:: dc_vap =  cp_vapor - c_liq   ! = -2368.
  real, parameter:: dc_vap =  cv_vap - c_liq   ! = -2368.
  real, parameter:: dc_ice =  c_liq - c_ice      ! = 2112.
! Values at 0 Deg C
  real, parameter:: hlv0 = 2.5e6
  real, parameter:: hlf0 = 3.3358e5
! real, parameter:: hlv0 = 2.501e6   ! Emanual Appendix-2
! real, parameter:: hlf0 = 3.337e5   ! Emanual
  real, parameter:: t_ice = 273.16
  real, parameter:: ri_max = 1.
  real, parameter:: ri_min = 0.25
  real, parameter:: t1_min = 160.
  real, parameter:: t2_min = 165.
  real, parameter:: t2_max = 315.
  real, parameter:: t3_max = 325.
  real, parameter:: Lv0 =  hlv0 - dc_vap*t_ice   ! = 3.147782e6
  real, parameter:: Li0 =  hlf0 - dc_ice*t_ice   ! = -2.431928e5 

  real, parameter:: zvir =  rvgas/rdgas - 1.     ! = 0.607789855
  real, allocatable:: table(:),des(:)
  real:: lv00, d0_vap

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of fv_subgrid_z in forward (tangent) mode:
!   variations   of useful results: qa w ta
!   with respect to varying inputs: qa peln w delp ua delz va pkz
!                pe ta
  SUBROUTINE FV_SUBGRID_Z_TLM(isd, ied, jsd, jed, is, ie, js, je, km, nq&
&   , dt, tau, nwat, delp, delp_tl, pe, pe_tl, peln, peln_tl, pkz, &
&   pkz_tl, ta, ta_tl, qa, qa_tl, ua, ua_tl, va, va_tl, hydrostatic, w, &
&   w_tl, delz, delz_tl, u_dt, v_dt, t_dt, k_bot)
    IMPLICIT NONE
! Dry convective adjustment-mixing
!-------------------------------------------
    INTEGER, INTENT(IN) :: is, ie, js, je, km, nq, nwat
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
! Relaxation time scale
    INTEGER, INTENT(IN) :: tau
! model time step
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(IN) :: pe_tl(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(IN) :: peln(is:ie, km+1, js:je)
    REAL, INTENT(IN) :: peln_tl(is:ie, km+1, js:je)
! Delta p at each model level
    REAL, INTENT(IN) :: delp(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delp_tl(isd:ied, jsd:jed, km)
! Delta z at each model level
    REAL, INTENT(IN) :: delz(isd:, jsd:, :)
    REAL, INTENT(IN) :: delz_tl(isd:, jsd:, :)
    REAL, INTENT(IN) :: pkz(is:ie, js:je, km)
    REAL, INTENT(IN) :: pkz_tl(is:ie, js:je, km)
    LOGICAL, INTENT(IN) :: hydrostatic
    INTEGER, INTENT(IN), OPTIONAL :: k_bot
!
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: ua_tl(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: va(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: va_tl(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: w(isd:, jsd:, :)
    REAL, INTENT(INOUT) :: w_tl(isd:, jsd:, :)
! Temperature
    REAL, INTENT(INOUT) :: ta(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: ta_tl(isd:ied, jsd:jed, km)
! Specific humidity & tracers
    REAL, INTENT(INOUT) :: qa(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: qa_tl(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: u_dt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_dt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: t_dt(is:ie, js:je, km)
!---------------------------Local variables-----------------------------
    REAL, DIMENSION(is:ie, km) :: u0, v0, w0, t0, hd, te, gz, tvm, pm, &
&   den
    REAL, DIMENSION(is:ie, km) :: u0_tl, v0_tl, w0_tl, t0_tl, hd_tl, &
&   te_tl, gz_tl, tvm_tl, pm_tl
    REAL :: q0(is:ie, km, nq), qcon(is:ie, km)
    REAL :: q0_tl(is:ie, km, nq), qcon_tl(is:ie, km)
    REAL, DIMENSION(is:ie) :: gzh, lcp2, icp2, cvm, cpm, qs
    REAL, DIMENSION(is:ie) :: gzh_tl, cvm_tl, cpm_tl
    REAL :: ri_ref, ri, pt1, pt2, ratio, tv, cv, tmp, q_liq, q_sol
    REAL :: ri_ref_tl, ri_tl, pt1_tl, pt2_tl, tv_tl, tmp_tl, q_liq_tl, &
&   q_sol_tl
    REAL :: tv1, tv2, g2, h0, mc, fra, rk, rz, rdt, tvd, tv_surf
    REAL :: tv1_tl, tv2_tl, h0_tl, mc_tl
    REAL :: dh, dq, qsw, dqsdt, tcp3, t_max, t_min
    INTEGER :: i, j, k, kk, n, m, iq, km1, im, kbot
    REAL, PARAMETER :: ustar2=1.e-4
    REAL :: cv_air, xvir
    INTEGER :: sphum, liq_wat, rainwat, snowwat, graupel, ice_wat, &
&   cld_amt
    INTRINSIC PRESENT
    INTRINSIC MIN
    INTRINSIC REAL
    INTRINSIC MAX
    INTEGER :: min1
    REAL :: max1
    REAL :: max1_tl
    REAL :: max2
    REAL :: max2_tl
    REAL :: y1_tl
    REAL :: y1
! = rdgas * (7/2-1) = 2.5*rdgas=717.68
    cv_air = cp_air - rdgas
    rk = cp_air/rdgas + 1.
    cv = cp_air - rdgas
    g2 = 0.5*grav
    rdt = 1./dt
    im = ie - is + 1
    IF (PRESENT(k_bot)) THEN
      IF (k_bot .LT. 3) THEN
        RETURN
      ELSE
        kbot = k_bot
      END IF
    ELSE
      kbot = km
    END IF
    IF (pe(is, 1, js) .LT. 2.) THEN
      t_min = t1_min
    ELSE
      t_min = t2_min
    END IF
    IF (km .GT. 24) THEN
      min1 = 24
    ELSE
      min1 = km
    END IF
    IF (k_bot .LT. min1) THEN
      t_max = t2_max
    ELSE
      t_max = t3_max
    END IF
    sphum = 1
    rainwat = -1
    snowwat = -1
    graupel = -1
    IF (nwat .EQ. 0) THEN
      xvir = 0.
      rz = 0.
    ELSE
      xvir = zvir
! rz = zvir * rdgas
      rz = rvgas - rdgas
      IF (nwat .EQ. 3) THEN
        liq_wat = 2
        ice_wat = 3
      END IF
    END IF
!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
    m = 3
    fra = dt/REAL(tau)
    cvm_tl = 0.0
    v0_tl = 0.0
    gz_tl = 0.0
    cpm_tl = 0.0
    hd_tl = 0.0
    w0_tl = 0.0
    tvm_tl = 0.0
    t0_tl = 0.0
    q0_tl = 0.0
    gzh_tl = 0.0
    qcon_tl = 0.0
    u0_tl = 0.0
    pm_tl = 0.0
    te_tl = 0.0
!$OMP parallel do default(none) shared(im,is,ie,js,je,nq,kbot,qa,ta,sphum,ua,va,delp,peln,   &
!$OMP                                  hydrostatic,pe,delz,g2,w,liq_wat,rainwat,ice_wat,     &
!$OMP                                  snowwat,cv_air,m,graupel,pkz,rk,rz,fra, t_max, t_min, &
!$OMP                                  u_dt,rdt,v_dt,xvir,nwat)                              &
!$OMP                          private(kk,lcp2,icp2,tcp3,dh,dq,den,qs,qsw,dqsdt,qcon,q0,     &
!$OMP                                  t0,u0,v0,w0,h0,pm,gzh,tvm,tmp,cpm,cvm,q_liq,q_sol,    &
!$OMP                                  tv,gz,hd,te,ratio,pt1,pt2,tv1,tv2,ri_ref, ri,mc,km1)
    DO j=js,je
      DO iq=1,nq
        DO k=1,kbot
          DO i=is,ie
            q0_tl(i, k, iq) = qa_tl(i, j, k, iq)
            q0(i, k, iq) = qa(i, j, k, iq)
          END DO
        END DO
      END DO
      DO k=1,kbot
        DO i=is,ie
          t0_tl(i, k) = ta_tl(i, j, k)
          t0(i, k) = ta(i, j, k)
          tvm_tl(i, k) = t0_tl(i, k)*(1.+xvir*q0(i, k, sphum)) + t0(i, k&
&           )*xvir*q0_tl(i, k, sphum)
          tvm(i, k) = t0(i, k)*(1.+xvir*q0(i, k, sphum))
          u0_tl(i, k) = ua_tl(i, j, k)
          u0(i, k) = ua(i, j, k)
          v0_tl(i, k) = va_tl(i, j, k)
          v0(i, k) = va(i, j, k)
          pm_tl(i, k) = (delp_tl(i, j, k)*(peln(i, k+1, j)-peln(i, k, j)&
&           )-delp(i, j, k)*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(peln&
&           (i, k+1, j)-peln(i, k, j))**2
          pm(i, k) = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
        END DO
      END DO
      DO i=is,ie
        gzh_tl(i) = 0.0
        gzh(i) = 0.
      END DO
      IF (hydrostatic) THEN
        DO k=kbot,1,-1
          DO i=is,ie
            tv_tl = rdgas*tvm_tl(i, k)
            tv = rdgas*tvm(i, k)
            den(i, k) = pm(i, k)/tv
            gz_tl(i, k) = gzh_tl(i) + tv_tl*(1.-pe(i, k, j)/pm(i, k)) - &
&             tv*(pe_tl(i, k, j)*pm(i, k)-pe(i, k, j)*pm_tl(i, k))/pm(i&
&             , k)**2
            gz(i, k) = gzh(i) + tv*(1.-pe(i, k, j)/pm(i, k))
            hd_tl(i, k) = cp_air*tvm_tl(i, k) + gz_tl(i, k) + 0.5*(2*u0(&
&             i, k)*u0_tl(i, k)+2*v0(i, k)*v0_tl(i, k))
            hd(i, k) = cp_air*tvm(i, k) + gz(i, k) + 0.5*(u0(i, k)**2+v0&
&             (i, k)**2)
            gzh_tl(i) = gzh_tl(i) + tv_tl*(peln(i, k+1, j)-peln(i, k, j)&
&             ) + tv*(peln_tl(i, k+1, j)-peln_tl(i, k, j))
            gzh(i) = gzh(i) + tv*(peln(i, k+1, j)-peln(i, k, j))
          END DO
        END DO
      ELSE
        DO k=kbot,1,-1
          IF (nwat .EQ. 0) THEN
            DO i=is,ie
              cpm_tl(i) = 0.0
              cpm(i) = cp_air
              cvm_tl(i) = 0.0
              cvm(i) = cv_air
            END DO
          ELSE IF (nwat .EQ. 1) THEN
            DO i=is,ie
              cpm_tl(i) = cp_vapor*q0_tl(i, k, sphum) - cp_air*q0_tl(i, &
&               k, sphum)
              cpm(i) = (1.-q0(i, k, sphum))*cp_air + q0(i, k, sphum)*&
&               cp_vapor
              cvm_tl(i) = cv_vap*q0_tl(i, k, sphum) - cv_air*q0_tl(i, k&
&               , sphum)
              cvm(i) = (1.-q0(i, k, sphum))*cv_air + q0(i, k, sphum)*&
&               cv_vap
            END DO
          ELSE IF (nwat .EQ. 2) THEN
! GFS
            DO i=is,ie
              cpm_tl(i) = cp_vapor*q0_tl(i, k, sphum) - cp_air*q0_tl(i, &
&               k, sphum)
              cpm(i) = (1.-q0(i, k, sphum))*cp_air + q0(i, k, sphum)*&
&               cp_vapor
              cvm_tl(i) = cv_vap*q0_tl(i, k, sphum) - cv_air*q0_tl(i, k&
&               , sphum)
              cvm(i) = (1.-q0(i, k, sphum))*cv_air + q0(i, k, sphum)*&
&               cv_vap
            END DO
          ELSE IF (nwat .EQ. 3) THEN
            DO i=is,ie
              q_liq_tl = q0_tl(i, k, liq_wat)
              q_liq = q0(i, k, liq_wat)
              q_sol_tl = q0_tl(i, k, ice_wat)
              q_sol = q0(i, k, ice_wat)
              cpm_tl(i) = cp_air*(-q0_tl(i, k, sphum)-q_liq_tl-q_sol_tl)&
&               + cp_vapor*q0_tl(i, k, sphum) + c_liq*q_liq_tl + c_ice*&
&               q_sol_tl
              cpm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cp_air + q0(i&
&               , k, sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
              cvm_tl(i) = cv_air*(-q0_tl(i, k, sphum)-q_liq_tl-q_sol_tl)&
&               + cv_vap*q0_tl(i, k, sphum) + c_liq*q_liq_tl + c_ice*&
&               q_sol_tl
              cvm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cv_air + q0(i&
&               , k, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
            END DO
          ELSE IF (nwat .EQ. 4) THEN
            DO i=is,ie
              q_liq_tl = q0_tl(i, k, liq_wat) + q0_tl(i, k, rainwat)
              q_liq = q0(i, k, liq_wat) + q0(i, k, rainwat)
              cpm_tl(i) = cp_air*(-q0_tl(i, k, sphum)-q_liq_tl) + &
&               cp_vapor*q0_tl(i, k, sphum) + c_liq*q_liq_tl
              cpm(i) = (1.-(q0(i, k, sphum)+q_liq))*cp_air + q0(i, k, &
&               sphum)*cp_vapor + q_liq*c_liq
              cvm_tl(i) = cv_air*(-q0_tl(i, k, sphum)-q_liq_tl) + cv_vap&
&               *q0_tl(i, k, sphum) + c_liq*q_liq_tl
              cvm(i) = (1.-(q0(i, k, sphum)+q_liq))*cv_air + q0(i, k, &
&               sphum)*cv_vap + q_liq*c_liq
            END DO
          ELSE
            DO i=is,ie
              q_liq_tl = q0_tl(i, k, liq_wat) + q0_tl(i, k, rainwat)
              q_liq = q0(i, k, liq_wat) + q0(i, k, rainwat)
              q_sol_tl = q0_tl(i, k, ice_wat) + q0_tl(i, k, snowwat) + &
&               q0_tl(i, k, graupel)
              q_sol = q0(i, k, ice_wat) + q0(i, k, snowwat) + q0(i, k, &
&               graupel)
              cpm_tl(i) = cp_air*(-q0_tl(i, k, sphum)-q_liq_tl-q_sol_tl)&
&               + cp_vapor*q0_tl(i, k, sphum) + c_liq*q_liq_tl + c_ice*&
&               q_sol_tl
              cpm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cp_air + q0(i&
&               , k, sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
              cvm_tl(i) = cv_air*(-q0_tl(i, k, sphum)-q_liq_tl-q_sol_tl)&
&               + cv_vap*q0_tl(i, k, sphum) + c_liq*q_liq_tl + c_ice*&
&               q_sol_tl
              cvm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cv_air + q0(i&
&               , k, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
            END DO
          END IF
          DO i=is,ie
            den(i, k) = -(delp(i, j, k)/(grav*delz(i, j, k)))
            w0_tl(i, k) = w_tl(i, j, k)
            w0(i, k) = w(i, j, k)
            gz_tl(i, k) = gzh_tl(i) - g2*delz_tl(i, j, k)
            gz(i, k) = gzh(i) - g2*delz(i, j, k)
            tmp_tl = gz_tl(i, k) + 0.5*(2*u0(i, k)*u0_tl(i, k)+2*v0(i, k&
&             )*v0_tl(i, k)+2*w0(i, k)*w0_tl(i, k))
            tmp = gz(i, k) + 0.5*(u0(i, k)**2+v0(i, k)**2+w0(i, k)**2)
            hd_tl(i, k) = cpm_tl(i)*t0(i, k) + cpm(i)*t0_tl(i, k) + &
&             tmp_tl
            hd(i, k) = cpm(i)*t0(i, k) + tmp
            te_tl(i, k) = cvm_tl(i)*t0(i, k) + cvm(i)*t0_tl(i, k) + &
&             tmp_tl
            te(i, k) = cvm(i)*t0(i, k) + tmp
            gzh_tl(i) = gzh_tl(i) - grav*delz_tl(i, j, k)
            gzh(i) = gzh(i) - grav*delz(i, j, k)
          END DO
        END DO
      END IF
      DO n=1,m
        IF (m .EQ. 3) THEN
          IF (n .EQ. 1) ratio = 0.25
          IF (n .EQ. 2) ratio = 0.5
          IF (n .EQ. 3) ratio = 0.999
        ELSE
          ratio = REAL(n)/REAL(m)
        END IF
        DO i=is,ie
          gzh_tl(i) = 0.0
          gzh(i) = 0.
        END DO
! Compute total condensate
        IF (nwat .LT. 2) THEN
          DO k=1,kbot
            DO i=is,ie
              qcon_tl(i, k) = 0.0
              qcon(i, k) = 0.
            END DO
          END DO
        ELSE IF (nwat .EQ. 2) THEN
! GFS_2015
          DO k=1,kbot
            DO i=is,ie
              qcon_tl(i, k) = q0_tl(i, k, liq_wat)
              qcon(i, k) = q0(i, k, liq_wat)
            END DO
          END DO
        ELSE IF (nwat .EQ. 3) THEN
          DO k=1,kbot
            DO i=is,ie
              qcon_tl(i, k) = q0_tl(i, k, liq_wat) + q0_tl(i, k, ice_wat&
&               )
              qcon(i, k) = q0(i, k, liq_wat) + q0(i, k, ice_wat)
            END DO
          END DO
        ELSE IF (nwat .EQ. 4) THEN
          DO k=1,kbot
            DO i=is,ie
              qcon_tl(i, k) = q0_tl(i, k, liq_wat) + q0_tl(i, k, rainwat&
&               )
              qcon(i, k) = q0(i, k, liq_wat) + q0(i, k, rainwat)
            END DO
          END DO
        ELSE
          DO k=1,kbot
            DO i=is,ie
              qcon_tl(i, k) = q0_tl(i, k, liq_wat) + q0_tl(i, k, ice_wat&
&               ) + q0_tl(i, k, snowwat) + q0_tl(i, k, rainwat) + q0_tl(&
&               i, k, graupel)
              qcon(i, k) = q0(i, k, liq_wat) + q0(i, k, ice_wat) + q0(i&
&               , k, snowwat) + q0(i, k, rainwat) + q0(i, k, graupel)
            END DO
          END DO
        END IF
        DO k=kbot,2,-1
          km1 = k - 1
          DO i=is,ie
! Richardson number = g*delz * del_theta/theta / (del_u**2 + del_v**2)
! Use exact form for "density temperature"
            tv1_tl = t0_tl(i, km1)*(1.+xvir*q0(i, km1, sphum)-qcon(i, &
&             km1)) + t0(i, km1)*(xvir*q0_tl(i, km1, sphum)-qcon_tl(i, &
&             km1))
            tv1 = t0(i, km1)*(1.+xvir*q0(i, km1, sphum)-qcon(i, km1))
            tv2_tl = t0_tl(i, k)*(1.+xvir*q0(i, k, sphum)-qcon(i, k)) + &
&             t0(i, k)*(xvir*q0_tl(i, k, sphum)-qcon_tl(i, k))
            tv2 = t0(i, k)*(1.+xvir*q0(i, k, sphum)-qcon(i, k))
            pt1_tl = (tv1_tl*pkz(i, j, km1)-tv1*pkz_tl(i, j, km1))/pkz(i&
&             , j, km1)**2
            pt1 = tv1/pkz(i, j, km1)
            pt2_tl = (tv2_tl*pkz(i, j, k)-tv2*pkz_tl(i, j, k))/pkz(i, j&
&             , k)**2
            pt2 = tv2/pkz(i, j, k)
!
            ri_tl = (((gz_tl(i, km1)-gz_tl(i, k))*(pt1-pt2)+(gz(i, km1)-&
&             gz(i, k))*(pt1_tl-pt2_tl))*0.5*(pt1+pt2)*((u0(i, km1)-u0(i&
&             , k))**2+(v0(i, km1)-v0(i, k))**2+ustar2)-(gz(i, km1)-gz(i&
&             , k))*(pt1-pt2)*0.5*((pt1_tl+pt2_tl)*((u0(i, km1)-u0(i, k)&
&             )**2+(v0(i, km1)-v0(i, k))**2+ustar2)+(pt1+pt2)*(2*(u0(i, &
&             km1)-u0(i, k))*(u0_tl(i, km1)-u0_tl(i, k))+2*(v0(i, km1)-&
&             v0(i, k))*(v0_tl(i, km1)-v0_tl(i, k)))))/(0.5*(pt1+pt2)*((&
&             u0(i, km1)-u0(i, k))**2+(v0(i, km1)-v0(i, k))**2+ustar2))&
&             **2
            ri = (gz(i, km1)-gz(i, k))*(pt1-pt2)/(0.5*(pt1+pt2)*((u0(i, &
&             km1)-u0(i, k))**2+(v0(i, km1)-v0(i, k))**2+ustar2))
            IF (tv1 .GT. t_max .AND. tv1 .GT. tv2) THEN
! top layer unphysically warm
              ri = 0.
              ri_tl = 0.0
            ELSE IF (tv2 .LT. t_min) THEN
              IF (ri .GT. 0.2) THEN
                ri = 0.2
                ri_tl = 0.0
              ELSE
                ri = ri
              END IF
            END IF
            IF (400.e2 - pm(i, k) .LT. 0.) THEN
              max2 = 0.
              max2_tl = 0.0
            ELSE
              max2_tl = -pm_tl(i, k)
              max2 = 400.e2 - pm(i, k)
            END IF
            y1_tl = (ri_max-ri_min)*max2_tl/200.e2
            y1 = ri_min + (ri_max-ri_min)*max2/200.e2
            IF (ri_max .GT. y1) THEN
              ri_ref_tl = y1_tl
              ri_ref = y1
            ELSE
              ri_ref = ri_max
              ri_ref_tl = 0.0
            END IF
! Enhancing mixing at the model top
            IF (k .EQ. 2) THEN
              ri_ref_tl = 4.*ri_ref_tl
              ri_ref = 4.*ri_ref
            ELSE IF (k .EQ. 3) THEN
              ri_ref_tl = 2.*ri_ref_tl
              ri_ref = 2.*ri_ref
            ELSE IF (k .EQ. 4) THEN
              ri_ref_tl = 1.5*ri_ref_tl
              ri_ref = 1.5*ri_ref
            END IF
            IF (ri .LT. ri_ref) THEN
              IF (0.0 .LT. ri/ri_ref) THEN
                max1_tl = (ri_tl*ri_ref-ri*ri_ref_tl)/ri_ref**2
                max1 = ri/ri_ref
              ELSE
                max1 = 0.0
                max1_tl = 0.0
              END IF
              mc_tl = (ratio*(delp_tl(i, j, km1)*delp(i, j, k)+delp(i, j&
&               , km1)*delp_tl(i, j, k))*(delp(i, j, km1)+delp(i, j, k))&
&               -ratio*delp(i, j, km1)*delp(i, j, k)*(delp_tl(i, j, km1)&
&               +delp_tl(i, j, k)))*(1.-max1)**2/(delp(i, j, km1)+delp(i&
&               , j, k))**2 - ratio*delp(i, j, km1)*delp(i, j, k)*2*(1.-&
&               max1)*max1_tl/(delp(i, j, km1)+delp(i, j, k))
              mc = ratio*delp(i, j, km1)*delp(i, j, k)/(delp(i, j, km1)+&
&               delp(i, j, k))*(1.-max1)**2
              DO iq=1,nq
                h0_tl = mc_tl*(q0(i, k, iq)-q0(i, km1, iq)) + mc*(q0_tl(&
&                 i, k, iq)-q0_tl(i, km1, iq))
                h0 = mc*(q0(i, k, iq)-q0(i, km1, iq))
                q0_tl(i, km1, iq) = q0_tl(i, km1, iq) + (h0_tl*delp(i, j&
&                 , km1)-h0*delp_tl(i, j, km1))/delp(i, j, km1)**2
                q0(i, km1, iq) = q0(i, km1, iq) + h0/delp(i, j, km1)
                q0_tl(i, k, iq) = q0_tl(i, k, iq) - (h0_tl*delp(i, j, k)&
&                 -h0*delp_tl(i, j, k))/delp(i, j, k)**2
                q0(i, k, iq) = q0(i, k, iq) - h0/delp(i, j, k)
              END DO
! Recompute qcon
              IF (nwat .LT. 2) THEN
                qcon_tl(i, km1) = 0.0
                qcon(i, km1) = 0.
              ELSE IF (nwat .EQ. 2) THEN
! GFS_2015
                qcon_tl(i, km1) = q0_tl(i, km1, liq_wat)
                qcon(i, km1) = q0(i, km1, liq_wat)
              ELSE IF (nwat .EQ. 3) THEN
! AM3/AM4
                qcon_tl(i, km1) = q0_tl(i, km1, liq_wat) + q0_tl(i, km1&
&                 , ice_wat)
                qcon(i, km1) = q0(i, km1, liq_wat) + q0(i, km1, ice_wat)
              ELSE IF (nwat .EQ. 4) THEN
! K_warm_rain scheme with fake ice
                qcon_tl(i, km1) = q0_tl(i, km1, liq_wat) + q0_tl(i, km1&
&                 , rainwat)
                qcon(i, km1) = q0(i, km1, liq_wat) + q0(i, km1, rainwat)
              ELSE
                qcon_tl(i, km1) = q0_tl(i, km1, liq_wat) + q0_tl(i, km1&
&                 , ice_wat) + q0_tl(i, km1, snowwat) + q0_tl(i, km1, &
&                 rainwat) + q0_tl(i, km1, graupel)
                qcon(i, km1) = q0(i, km1, liq_wat) + q0(i, km1, ice_wat)&
&                 + q0(i, km1, snowwat) + q0(i, km1, rainwat) + q0(i, &
&                 km1, graupel)
              END IF
! u:
              h0_tl = mc_tl*(u0(i, k)-u0(i, k-1)) + mc*(u0_tl(i, k)-&
&               u0_tl(i, k-1))
              h0 = mc*(u0(i, k)-u0(i, k-1))
              u0_tl(i, k-1) = u0_tl(i, k-1) + (h0_tl*delp(i, j, k-1)-h0*&
&               delp_tl(i, j, k-1))/delp(i, j, k-1)**2
              u0(i, k-1) = u0(i, k-1) + h0/delp(i, j, k-1)
              u0_tl(i, k) = u0_tl(i, k) - (h0_tl*delp(i, j, k)-h0*&
&               delp_tl(i, j, k))/delp(i, j, k)**2
              u0(i, k) = u0(i, k) - h0/delp(i, j, k)
! v:
              h0_tl = mc_tl*(v0(i, k)-v0(i, k-1)) + mc*(v0_tl(i, k)-&
&               v0_tl(i, k-1))
              h0 = mc*(v0(i, k)-v0(i, k-1))
              v0_tl(i, k-1) = v0_tl(i, k-1) + (h0_tl*delp(i, j, k-1)-h0*&
&               delp_tl(i, j, k-1))/delp(i, j, k-1)**2
              v0(i, k-1) = v0(i, k-1) + h0/delp(i, j, k-1)
              v0_tl(i, k) = v0_tl(i, k) - (h0_tl*delp(i, j, k)-h0*&
&               delp_tl(i, j, k))/delp(i, j, k)**2
              v0(i, k) = v0(i, k) - h0/delp(i, j, k)
              IF (hydrostatic) THEN
! Static energy
                h0_tl = mc_tl*(hd(i, k)-hd(i, k-1)) + mc*(hd_tl(i, k)-&
&                 hd_tl(i, k-1))
                h0 = mc*(hd(i, k)-hd(i, k-1))
                hd_tl(i, k-1) = hd_tl(i, k-1) + (h0_tl*delp(i, j, k-1)-&
&                 h0*delp_tl(i, j, k-1))/delp(i, j, k-1)**2
                hd(i, k-1) = hd(i, k-1) + h0/delp(i, j, k-1)
                hd_tl(i, k) = hd_tl(i, k) - (h0_tl*delp(i, j, k)-h0*&
&                 delp_tl(i, j, k))/delp(i, j, k)**2
                hd(i, k) = hd(i, k) - h0/delp(i, j, k)
              ELSE
! Total energy
                h0_tl = mc_tl*(hd(i, k)-hd(i, k-1)) + mc*(hd_tl(i, k)-&
&                 hd_tl(i, k-1))
                h0 = mc*(hd(i, k)-hd(i, k-1))
                te_tl(i, k-1) = te_tl(i, k-1) + (h0_tl*delp(i, j, k-1)-&
&                 h0*delp_tl(i, j, k-1))/delp(i, j, k-1)**2
                te(i, k-1) = te(i, k-1) + h0/delp(i, j, k-1)
                te_tl(i, k) = te_tl(i, k) - (h0_tl*delp(i, j, k)-h0*&
&                 delp_tl(i, j, k))/delp(i, j, k)**2
                te(i, k) = te(i, k) - h0/delp(i, j, k)
! w:
                h0_tl = mc_tl*(w0(i, k)-w0(i, k-1)) + mc*(w0_tl(i, k)-&
&                 w0_tl(i, k-1))
                h0 = mc*(w0(i, k)-w0(i, k-1))
                w0_tl(i, k-1) = w0_tl(i, k-1) + (h0_tl*delp(i, j, k-1)-&
&                 h0*delp_tl(i, j, k-1))/delp(i, j, k-1)**2
                w0(i, k-1) = w0(i, k-1) + h0/delp(i, j, k-1)
                w0_tl(i, k) = w0_tl(i, k) - (h0_tl*delp(i, j, k)-h0*&
&                 delp_tl(i, j, k))/delp(i, j, k)**2
                w0(i, k) = w0(i, k) - h0/delp(i, j, k)
              END IF
            END IF
          END DO
!--------------
! Retrive Temp:
!--------------
          IF (hydrostatic) THEN
            kk = k
            DO i=is,ie
              t0_tl(i, kk) = ((hd_tl(i, kk)-gzh_tl(i)-0.5*(2*u0(i, kk)*&
&               u0_tl(i, kk)+2*v0(i, kk)*v0_tl(i, kk)))*(rk-pe(i, kk, j)&
&               /pm(i, kk))+(hd(i, kk)-gzh(i)-0.5*(u0(i, kk)**2+v0(i, kk&
&               )**2))*(pe_tl(i, kk, j)*pm(i, kk)-pe(i, kk, j)*pm_tl(i, &
&               kk))/pm(i, kk)**2)/(rk-pe(i, kk, j)/pm(i, kk))**2
              t0(i, kk) = (hd(i, kk)-gzh(i)-0.5*(u0(i, kk)**2+v0(i, kk)&
&               **2))/(rk-pe(i, kk, j)/pm(i, kk))
              gzh_tl(i) = gzh_tl(i) + t0_tl(i, kk)*(peln(i, kk+1, j)-&
&               peln(i, kk, j)) + t0(i, kk)*(peln_tl(i, kk+1, j)-peln_tl&
&               (i, kk, j))
              gzh(i) = gzh(i) + t0(i, kk)*(peln(i, kk+1, j)-peln(i, kk, &
&               j))
              t0_tl(i, kk) = (t0_tl(i, kk)*(rdgas+rz*q0(i, kk, sphum))-&
&               t0(i, kk)*rz*q0_tl(i, kk, sphum))/(rdgas+rz*q0(i, kk, &
&               sphum))**2
              t0(i, kk) = t0(i, kk)/(rdgas+rz*q0(i, kk, sphum))
            END DO
            kk = k - 1
            DO i=is,ie
              t0_tl(i, kk) = ((hd_tl(i, kk)-gzh_tl(i)-0.5*(2*u0(i, kk)*&
&               u0_tl(i, kk)+2*v0(i, kk)*v0_tl(i, kk)))*(rk-pe(i, kk, j)&
&               /pm(i, kk))*(rdgas+rz*q0(i, kk, sphum))-(hd(i, kk)-gzh(i&
&               )-0.5*(u0(i, kk)**2+v0(i, kk)**2))*((rk-pe(i, kk, j)/pm(&
&               i, kk))*rz*q0_tl(i, kk, sphum)-(pe_tl(i, kk, j)*pm(i, kk&
&               )-pe(i, kk, j)*pm_tl(i, kk))*(rdgas+rz*q0(i, kk, sphum))&
&               /pm(i, kk)**2))/((rk-pe(i, kk, j)/pm(i, kk))*(rdgas+rz*&
&               q0(i, kk, sphum)))**2
              t0(i, kk) = (hd(i, kk)-gzh(i)-0.5*(u0(i, kk)**2+v0(i, kk)&
&               **2))/((rk-pe(i, kk, j)/pm(i, kk))*(rdgas+rz*q0(i, kk, &
&               sphum)))
            END DO
          ELSE
! Non-hydrostatic under constant volume heating/cooling
            DO kk=k-1,k
              IF (nwat .EQ. 0) THEN
                DO i=is,ie
                  cpm_tl(i) = 0.0
                  cpm(i) = cp_air
                  cvm_tl(i) = 0.0
                  cvm(i) = cv_air
                END DO
              ELSE IF (nwat .EQ. 1) THEN
                DO i=is,ie
                  cpm_tl(i) = cp_vapor*q0_tl(i, kk, sphum) - cp_air*&
&                   q0_tl(i, kk, sphum)
                  cpm(i) = (1.-q0(i, kk, sphum))*cp_air + q0(i, kk, &
&                   sphum)*cp_vapor
                  cvm_tl(i) = cv_vap*q0_tl(i, kk, sphum) - cv_air*q0_tl(&
&                   i, kk, sphum)
                  cvm(i) = (1.-q0(i, kk, sphum))*cv_air + q0(i, kk, &
&                   sphum)*cv_vap
                END DO
              ELSE IF (nwat .EQ. 2) THEN
                DO i=is,ie
                  cpm_tl(i) = cp_vapor*q0_tl(i, kk, sphum) - cp_air*&
&                   q0_tl(i, kk, sphum)
                  cpm(i) = (1.-q0(i, kk, sphum))*cp_air + q0(i, kk, &
&                   sphum)*cp_vapor
                  cvm_tl(i) = cv_vap*q0_tl(i, kk, sphum) - cv_air*q0_tl(&
&                   i, kk, sphum)
                  cvm(i) = (1.-q0(i, kk, sphum))*cv_air + q0(i, kk, &
&                   sphum)*cv_vap
                END DO
              ELSE IF (nwat .EQ. 3) THEN
                DO i=is,ie
                  q_liq_tl = q0_tl(i, kk, liq_wat)
                  q_liq = q0(i, kk, liq_wat)
                  q_sol_tl = q0_tl(i, kk, ice_wat)
                  q_sol = q0(i, kk, ice_wat)
                  cpm_tl(i) = cp_air*(-q0_tl(i, kk, sphum)-q_liq_tl-&
&                   q_sol_tl) + cp_vapor*q0_tl(i, kk, sphum) + c_liq*&
&                   q_liq_tl + c_ice*q_sol_tl
                  cpm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cp_air + &
&                   q0(i, kk, sphum)*cp_vapor + q_liq*c_liq + q_sol*&
&                   c_ice
                  cvm_tl(i) = cv_air*(-q0_tl(i, kk, sphum)-q_liq_tl-&
&                   q_sol_tl) + cv_vap*q0_tl(i, kk, sphum) + c_liq*&
&                   q_liq_tl + c_ice*q_sol_tl
                  cvm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cv_air + &
&                   q0(i, kk, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
                END DO
              ELSE IF (nwat .EQ. 4) THEN
                DO i=is,ie
                  q_liq_tl = q0_tl(i, kk, liq_wat) + q0_tl(i, kk, &
&                   rainwat)
                  q_liq = q0(i, kk, liq_wat) + q0(i, kk, rainwat)
                  cpm_tl(i) = cp_air*(-q0_tl(i, kk, sphum)-q_liq_tl) + &
&                   cp_vapor*q0_tl(i, kk, sphum) + c_liq*q_liq_tl
                  cpm(i) = (1.-(q0(i, kk, sphum)+q_liq))*cp_air + q0(i, &
&                   kk, sphum)*cp_vapor + q_liq*c_liq
                  cvm_tl(i) = cv_air*(-q0_tl(i, kk, sphum)-q_liq_tl) + &
&                   cv_vap*q0_tl(i, kk, sphum) + c_liq*q_liq_tl
                  cvm(i) = (1.-(q0(i, kk, sphum)+q_liq))*cv_air + q0(i, &
&                   kk, sphum)*cv_vap + q_liq*c_liq
                END DO
              ELSE
                DO i=is,ie
                  q_liq_tl = q0_tl(i, kk, liq_wat) + q0_tl(i, kk, &
&                   rainwat)
                  q_liq = q0(i, kk, liq_wat) + q0(i, kk, rainwat)
                  q_sol_tl = q0_tl(i, kk, ice_wat) + q0_tl(i, kk, &
&                   snowwat) + q0_tl(i, kk, graupel)
                  q_sol = q0(i, kk, ice_wat) + q0(i, kk, snowwat) + q0(i&
&                   , kk, graupel)
                  cpm_tl(i) = cp_air*(-q0_tl(i, kk, sphum)-q_liq_tl-&
&                   q_sol_tl) + cp_vapor*q0_tl(i, kk, sphum) + c_liq*&
&                   q_liq_tl + c_ice*q_sol_tl
                  cpm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cp_air + &
&                   q0(i, kk, sphum)*cp_vapor + q_liq*c_liq + q_sol*&
&                   c_ice
                  cvm_tl(i) = cv_air*(-q0_tl(i, kk, sphum)-q_liq_tl-&
&                   q_sol_tl) + cv_vap*q0_tl(i, kk, sphum) + c_liq*&
&                   q_liq_tl + c_ice*q_sol_tl
                  cvm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cv_air + &
&                   q0(i, kk, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
                END DO
              END IF
              DO i=is,ie
                tv_tl = gz_tl(i, kk) + 0.5*(2*u0(i, kk)*u0_tl(i, kk)+2*&
&                 v0(i, kk)*v0_tl(i, kk)+2*w0(i, kk)*w0_tl(i, kk))
                tv = gz(i, kk) + 0.5*(u0(i, kk)**2+v0(i, kk)**2+w0(i, kk&
&                 )**2)
                t0_tl(i, kk) = ((te_tl(i, kk)-tv_tl)*cvm(i)-(te(i, kk)-&
&                 tv)*cvm_tl(i))/cvm(i)**2
                t0(i, kk) = (te(i, kk)-tv)/cvm(i)
                hd_tl(i, kk) = cpm_tl(i)*t0(i, kk) + cpm(i)*t0_tl(i, kk)&
&                 + tv_tl
                hd(i, kk) = cpm(i)*t0(i, kk) + tv
              END DO
            END DO
          END IF
        END DO
      END DO
! k-loop
! n-loop
!--------------------
      IF (fra .LT. 1.) THEN
        DO k=1,kbot
          DO i=is,ie
            t0_tl(i, k) = ta_tl(i, j, k) + fra*(t0_tl(i, k)-ta_tl(i, j, &
&             k))
            t0(i, k) = ta(i, j, k) + (t0(i, k)-ta(i, j, k))*fra
            u0_tl(i, k) = ua_tl(i, j, k) + fra*(u0_tl(i, k)-ua_tl(i, j, &
&             k))
            u0(i, k) = ua(i, j, k) + (u0(i, k)-ua(i, j, k))*fra
            v0_tl(i, k) = va_tl(i, j, k) + fra*(v0_tl(i, k)-va_tl(i, j, &
&             k))
            v0(i, k) = va(i, j, k) + (v0(i, k)-va(i, j, k))*fra
          END DO
        END DO
        IF (.NOT.hydrostatic) THEN
          DO k=1,kbot
            DO i=is,ie
              w0_tl(i, k) = w_tl(i, j, k) + fra*(w0_tl(i, k)-w_tl(i, j, &
&               k))
              w0(i, k) = w(i, j, k) + (w0(i, k)-w(i, j, k))*fra
            END DO
          END DO
        END IF
        DO iq=1,nq
          DO k=1,kbot
            DO i=is,ie
              q0_tl(i, k, iq) = qa_tl(i, j, k, iq) + fra*(q0_tl(i, k, iq&
&               )-qa_tl(i, j, k, iq))
              q0(i, k, iq) = qa(i, j, k, iq) + (q0(i, k, iq)-qa(i, j, k&
&               , iq))*fra
            END DO
          END DO
        END DO
      END IF
      DO k=1,kbot
        DO i=is,ie
          u_dt(i, j, k) = rdt*(u0(i, k)-ua(i, j, k))
          v_dt(i, j, k) = rdt*(v0(i, k)-va(i, j, k))
! *** temperature updated ***
          ta_tl(i, j, k) = t0_tl(i, k)
          ta(i, j, k) = t0(i, k)
          ua_tl(i, j, k) = u0_tl(i, k)
          ua(i, j, k) = u0(i, k)
          va_tl(i, j, k) = v0_tl(i, k)
          va(i, j, k) = v0(i, k)
        END DO
        DO iq=1,nq
          DO i=is,ie
            qa_tl(i, j, k, iq) = q0_tl(i, k, iq)
            qa(i, j, k, iq) = q0(i, k, iq)
          END DO
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        DO k=1,kbot
          DO i=is,ie
! w updated
            w_tl(i, j, k) = w0_tl(i, k)
            w(i, j, k) = w0(i, k)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE FV_SUBGRID_Z_TLM
  SUBROUTINE FV_SUBGRID_Z(isd, ied, jsd, jed, is, ie, js, je, km, nq, dt&
&   , tau, nwat, delp, pe, peln, pkz, ta, qa, ua, va, hydrostatic, w, &
&   delz, u_dt, v_dt, t_dt, k_bot)
    IMPLICIT NONE
! Dry convective adjustment-mixing
!-------------------------------------------
    INTEGER, INTENT(IN) :: is, ie, js, je, km, nq, nwat
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
! Relaxation time scale
    INTEGER, INTENT(IN) :: tau
! model time step
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(IN) :: peln(is:ie, km+1, js:je)
! Delta p at each model level
    REAL, INTENT(IN) :: delp(isd:ied, jsd:jed, km)
! Delta z at each model level
    REAL, INTENT(IN) :: delz(isd:, jsd:, :)
    REAL, INTENT(IN) :: pkz(is:ie, js:je, km)
    LOGICAL, INTENT(IN) :: hydrostatic
    INTEGER, INTENT(IN), OPTIONAL :: k_bot
!
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: va(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: w(isd:, jsd:, :)
! Temperature
    REAL, INTENT(INOUT) :: ta(isd:ied, jsd:jed, km)
! Specific humidity & tracers
    REAL, INTENT(INOUT) :: qa(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: u_dt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_dt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: t_dt(is:ie, js:je, km)
!---------------------------Local variables-----------------------------
    REAL, DIMENSION(is:ie, km) :: u0, v0, w0, t0, hd, te, gz, tvm, pm, &
&   den
    REAL :: q0(is:ie, km, nq), qcon(is:ie, km)
    REAL, DIMENSION(is:ie) :: gzh, lcp2, icp2, cvm, cpm, qs
    REAL :: ri_ref, ri, pt1, pt2, ratio, tv, cv, tmp, q_liq, q_sol
    REAL :: tv1, tv2, g2, h0, mc, fra, rk, rz, rdt, tvd, tv_surf
    REAL :: dh, dq, qsw, dqsdt, tcp3, t_max, t_min
    INTEGER :: i, j, k, kk, n, m, iq, km1, im, kbot
    REAL, PARAMETER :: ustar2=1.e-4
    REAL :: cv_air, xvir
    INTEGER :: sphum, liq_wat, rainwat, snowwat, graupel, ice_wat, &
&   cld_amt
    INTRINSIC PRESENT
    INTRINSIC MIN
    INTRINSIC REAL
    INTRINSIC MAX
    INTEGER :: min1
    REAL :: max1
    REAL :: max2
    REAL :: y1
! = rdgas * (7/2-1) = 2.5*rdgas=717.68
    cv_air = cp_air - rdgas
    rk = cp_air/rdgas + 1.
    cv = cp_air - rdgas
    g2 = 0.5*grav
    rdt = 1./dt
    im = ie - is + 1
    IF (PRESENT(k_bot)) THEN
      IF (k_bot .LT. 3) THEN
        RETURN
      ELSE
        kbot = k_bot
      END IF
    ELSE
      kbot = km
    END IF
    IF (pe(is, 1, js) .LT. 2.) THEN
      t_min = t1_min
    ELSE
      t_min = t2_min
    END IF
    IF (km .GT. 24) THEN
      min1 = 24
    ELSE
      min1 = km
    END IF
    IF (k_bot .LT. min1) THEN
      t_max = t2_max
    ELSE
      t_max = t3_max
    END IF
    sphum = 1
    rainwat = -1
    snowwat = -1
    graupel = -1
    IF (nwat .EQ. 0) THEN
      xvir = 0.
      rz = 0.
    ELSE
      xvir = zvir
! rz = zvir * rdgas
      rz = rvgas - rdgas
      IF (nwat .EQ. 3) THEN
        liq_wat = 2
        ice_wat = 3
      END IF
    END IF
!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
    m = 3
    fra = dt/REAL(tau)
!$OMP parallel do default(none) shared(im,is,ie,js,je,nq,kbot,qa,ta,sphum,ua,va,delp,peln,   &
!$OMP                                  hydrostatic,pe,delz,g2,w,liq_wat,rainwat,ice_wat,     &
!$OMP                                  snowwat,cv_air,m,graupel,pkz,rk,rz,fra, t_max, t_min, &
!$OMP                                  u_dt,rdt,v_dt,xvir,nwat)                              &
!$OMP                          private(kk,lcp2,icp2,tcp3,dh,dq,den,qs,qsw,dqsdt,qcon,q0,     &
!$OMP                                  t0,u0,v0,w0,h0,pm,gzh,tvm,tmp,cpm,cvm,q_liq,q_sol,    &
!$OMP                                  tv,gz,hd,te,ratio,pt1,pt2,tv1,tv2,ri_ref, ri,mc,km1)
    DO j=js,je
      DO iq=1,nq
        DO k=1,kbot
          DO i=is,ie
            q0(i, k, iq) = qa(i, j, k, iq)
          END DO
        END DO
      END DO
      DO k=1,kbot
        DO i=is,ie
          t0(i, k) = ta(i, j, k)
          tvm(i, k) = t0(i, k)*(1.+xvir*q0(i, k, sphum))
          u0(i, k) = ua(i, j, k)
          v0(i, k) = va(i, j, k)
          pm(i, k) = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
        END DO
      END DO
      DO i=is,ie
        gzh(i) = 0.
      END DO
      IF (hydrostatic) THEN
        DO k=kbot,1,-1
          DO i=is,ie
            tv = rdgas*tvm(i, k)
            den(i, k) = pm(i, k)/tv
            gz(i, k) = gzh(i) + tv*(1.-pe(i, k, j)/pm(i, k))
            hd(i, k) = cp_air*tvm(i, k) + gz(i, k) + 0.5*(u0(i, k)**2+v0&
&             (i, k)**2)
            gzh(i) = gzh(i) + tv*(peln(i, k+1, j)-peln(i, k, j))
          END DO
        END DO
      ELSE
        DO k=kbot,1,-1
          IF (nwat .EQ. 0) THEN
            DO i=is,ie
              cpm(i) = cp_air
              cvm(i) = cv_air
            END DO
          ELSE IF (nwat .EQ. 1) THEN
            DO i=is,ie
              cpm(i) = (1.-q0(i, k, sphum))*cp_air + q0(i, k, sphum)*&
&               cp_vapor
              cvm(i) = (1.-q0(i, k, sphum))*cv_air + q0(i, k, sphum)*&
&               cv_vap
            END DO
          ELSE IF (nwat .EQ. 2) THEN
! GFS
            DO i=is,ie
              cpm(i) = (1.-q0(i, k, sphum))*cp_air + q0(i, k, sphum)*&
&               cp_vapor
              cvm(i) = (1.-q0(i, k, sphum))*cv_air + q0(i, k, sphum)*&
&               cv_vap
            END DO
          ELSE IF (nwat .EQ. 3) THEN
            DO i=is,ie
              q_liq = q0(i, k, liq_wat)
              q_sol = q0(i, k, ice_wat)
              cpm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cp_air + q0(i&
&               , k, sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
              cvm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cv_air + q0(i&
&               , k, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
            END DO
          ELSE IF (nwat .EQ. 4) THEN
            DO i=is,ie
              q_liq = q0(i, k, liq_wat) + q0(i, k, rainwat)
              cpm(i) = (1.-(q0(i, k, sphum)+q_liq))*cp_air + q0(i, k, &
&               sphum)*cp_vapor + q_liq*c_liq
              cvm(i) = (1.-(q0(i, k, sphum)+q_liq))*cv_air + q0(i, k, &
&               sphum)*cv_vap + q_liq*c_liq
            END DO
          ELSE
            DO i=is,ie
              q_liq = q0(i, k, liq_wat) + q0(i, k, rainwat)
              q_sol = q0(i, k, ice_wat) + q0(i, k, snowwat) + q0(i, k, &
&               graupel)
              cpm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cp_air + q0(i&
&               , k, sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
              cvm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cv_air + q0(i&
&               , k, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
            END DO
          END IF
          DO i=is,ie
            den(i, k) = -(delp(i, j, k)/(grav*delz(i, j, k)))
            w0(i, k) = w(i, j, k)
            gz(i, k) = gzh(i) - g2*delz(i, j, k)
            tmp = gz(i, k) + 0.5*(u0(i, k)**2+v0(i, k)**2+w0(i, k)**2)
            hd(i, k) = cpm(i)*t0(i, k) + tmp
            te(i, k) = cvm(i)*t0(i, k) + tmp
            gzh(i) = gzh(i) - grav*delz(i, j, k)
          END DO
        END DO
      END IF
      DO n=1,m
        IF (m .EQ. 3) THEN
          IF (n .EQ. 1) ratio = 0.25
          IF (n .EQ. 2) ratio = 0.5
          IF (n .EQ. 3) ratio = 0.999
        ELSE
          ratio = REAL(n)/REAL(m)
        END IF
        DO i=is,ie
          gzh(i) = 0.
        END DO
! Compute total condensate
        IF (nwat .LT. 2) THEN
          DO k=1,kbot
            DO i=is,ie
              qcon(i, k) = 0.
            END DO
          END DO
        ELSE IF (nwat .EQ. 2) THEN
! GFS_2015
          DO k=1,kbot
            DO i=is,ie
              qcon(i, k) = q0(i, k, liq_wat)
            END DO
          END DO
        ELSE IF (nwat .EQ. 3) THEN
          DO k=1,kbot
            DO i=is,ie
              qcon(i, k) = q0(i, k, liq_wat) + q0(i, k, ice_wat)
            END DO
          END DO
        ELSE IF (nwat .EQ. 4) THEN
          DO k=1,kbot
            DO i=is,ie
              qcon(i, k) = q0(i, k, liq_wat) + q0(i, k, rainwat)
            END DO
          END DO
        ELSE
          DO k=1,kbot
            DO i=is,ie
              qcon(i, k) = q0(i, k, liq_wat) + q0(i, k, ice_wat) + q0(i&
&               , k, snowwat) + q0(i, k, rainwat) + q0(i, k, graupel)
            END DO
          END DO
        END IF
        DO k=kbot,2,-1
          km1 = k - 1
          DO i=is,ie
! Richardson number = g*delz * del_theta/theta / (del_u**2 + del_v**2)
! Use exact form for "density temperature"
            tv1 = t0(i, km1)*(1.+xvir*q0(i, km1, sphum)-qcon(i, km1))
            tv2 = t0(i, k)*(1.+xvir*q0(i, k, sphum)-qcon(i, k))
            pt1 = tv1/pkz(i, j, km1)
            pt2 = tv2/pkz(i, j, k)
!
            ri = (gz(i, km1)-gz(i, k))*(pt1-pt2)/(0.5*(pt1+pt2)*((u0(i, &
&             km1)-u0(i, k))**2+(v0(i, km1)-v0(i, k))**2+ustar2))
            IF (tv1 .GT. t_max .AND. tv1 .GT. tv2) THEN
! top layer unphysically warm
              ri = 0.
            ELSE IF (tv2 .LT. t_min) THEN
              IF (ri .GT. 0.2) THEN
                ri = 0.2
              ELSE
                ri = ri
              END IF
            END IF
            IF (400.e2 - pm(i, k) .LT. 0.) THEN
              max2 = 0.
            ELSE
              max2 = 400.e2 - pm(i, k)
            END IF
            y1 = ri_min + (ri_max-ri_min)*max2/200.e2
            IF (ri_max .GT. y1) THEN
              ri_ref = y1
            ELSE
              ri_ref = ri_max
            END IF
! Enhancing mixing at the model top
            IF (k .EQ. 2) THEN
              ri_ref = 4.*ri_ref
            ELSE IF (k .EQ. 3) THEN
              ri_ref = 2.*ri_ref
            ELSE IF (k .EQ. 4) THEN
              ri_ref = 1.5*ri_ref
            END IF
            IF (ri .LT. ri_ref) THEN
              IF (0.0 .LT. ri/ri_ref) THEN
                max1 = ri/ri_ref
              ELSE
                max1 = 0.0
              END IF
              mc = ratio*delp(i, j, km1)*delp(i, j, k)/(delp(i, j, km1)+&
&               delp(i, j, k))*(1.-max1)**2
              DO iq=1,nq
                h0 = mc*(q0(i, k, iq)-q0(i, km1, iq))
                q0(i, km1, iq) = q0(i, km1, iq) + h0/delp(i, j, km1)
                q0(i, k, iq) = q0(i, k, iq) - h0/delp(i, j, k)
              END DO
! Recompute qcon
              IF (nwat .LT. 2) THEN
                qcon(i, km1) = 0.
              ELSE IF (nwat .EQ. 2) THEN
! GFS_2015
                qcon(i, km1) = q0(i, km1, liq_wat)
              ELSE IF (nwat .EQ. 3) THEN
! AM3/AM4
                qcon(i, km1) = q0(i, km1, liq_wat) + q0(i, km1, ice_wat)
              ELSE IF (nwat .EQ. 4) THEN
! K_warm_rain scheme with fake ice
                qcon(i, km1) = q0(i, km1, liq_wat) + q0(i, km1, rainwat)
              ELSE
                qcon(i, km1) = q0(i, km1, liq_wat) + q0(i, km1, ice_wat)&
&                 + q0(i, km1, snowwat) + q0(i, km1, rainwat) + q0(i, &
&                 km1, graupel)
              END IF
! u:
              h0 = mc*(u0(i, k)-u0(i, k-1))
              u0(i, k-1) = u0(i, k-1) + h0/delp(i, j, k-1)
              u0(i, k) = u0(i, k) - h0/delp(i, j, k)
! v:
              h0 = mc*(v0(i, k)-v0(i, k-1))
              v0(i, k-1) = v0(i, k-1) + h0/delp(i, j, k-1)
              v0(i, k) = v0(i, k) - h0/delp(i, j, k)
              IF (hydrostatic) THEN
! Static energy
                h0 = mc*(hd(i, k)-hd(i, k-1))
                hd(i, k-1) = hd(i, k-1) + h0/delp(i, j, k-1)
                hd(i, k) = hd(i, k) - h0/delp(i, j, k)
              ELSE
! Total energy
                h0 = mc*(hd(i, k)-hd(i, k-1))
                te(i, k-1) = te(i, k-1) + h0/delp(i, j, k-1)
                te(i, k) = te(i, k) - h0/delp(i, j, k)
! w:
                h0 = mc*(w0(i, k)-w0(i, k-1))
                w0(i, k-1) = w0(i, k-1) + h0/delp(i, j, k-1)
                w0(i, k) = w0(i, k) - h0/delp(i, j, k)
              END IF
            END IF
          END DO
!--------------
! Retrive Temp:
!--------------
          IF (hydrostatic) THEN
            kk = k
            DO i=is,ie
              t0(i, kk) = (hd(i, kk)-gzh(i)-0.5*(u0(i, kk)**2+v0(i, kk)&
&               **2))/(rk-pe(i, kk, j)/pm(i, kk))
              gzh(i) = gzh(i) + t0(i, kk)*(peln(i, kk+1, j)-peln(i, kk, &
&               j))
              t0(i, kk) = t0(i, kk)/(rdgas+rz*q0(i, kk, sphum))
            END DO
            kk = k - 1
            DO i=is,ie
              t0(i, kk) = (hd(i, kk)-gzh(i)-0.5*(u0(i, kk)**2+v0(i, kk)&
&               **2))/((rk-pe(i, kk, j)/pm(i, kk))*(rdgas+rz*q0(i, kk, &
&               sphum)))
            END DO
          ELSE
! Non-hydrostatic under constant volume heating/cooling
            DO kk=k-1,k
              IF (nwat .EQ. 0) THEN
                DO i=is,ie
                  cpm(i) = cp_air
                  cvm(i) = cv_air
                END DO
              ELSE IF (nwat .EQ. 1) THEN
                DO i=is,ie
                  cpm(i) = (1.-q0(i, kk, sphum))*cp_air + q0(i, kk, &
&                   sphum)*cp_vapor
                  cvm(i) = (1.-q0(i, kk, sphum))*cv_air + q0(i, kk, &
&                   sphum)*cv_vap
                END DO
              ELSE IF (nwat .EQ. 2) THEN
                DO i=is,ie
                  cpm(i) = (1.-q0(i, kk, sphum))*cp_air + q0(i, kk, &
&                   sphum)*cp_vapor
                  cvm(i) = (1.-q0(i, kk, sphum))*cv_air + q0(i, kk, &
&                   sphum)*cv_vap
                END DO
              ELSE IF (nwat .EQ. 3) THEN
                DO i=is,ie
                  q_liq = q0(i, kk, liq_wat)
                  q_sol = q0(i, kk, ice_wat)
                  cpm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cp_air + &
&                   q0(i, kk, sphum)*cp_vapor + q_liq*c_liq + q_sol*&
&                   c_ice
                  cvm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cv_air + &
&                   q0(i, kk, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
                END DO
              ELSE IF (nwat .EQ. 4) THEN
                DO i=is,ie
                  q_liq = q0(i, kk, liq_wat) + q0(i, kk, rainwat)
                  cpm(i) = (1.-(q0(i, kk, sphum)+q_liq))*cp_air + q0(i, &
&                   kk, sphum)*cp_vapor + q_liq*c_liq
                  cvm(i) = (1.-(q0(i, kk, sphum)+q_liq))*cv_air + q0(i, &
&                   kk, sphum)*cv_vap + q_liq*c_liq
                END DO
              ELSE
                DO i=is,ie
                  q_liq = q0(i, kk, liq_wat) + q0(i, kk, rainwat)
                  q_sol = q0(i, kk, ice_wat) + q0(i, kk, snowwat) + q0(i&
&                   , kk, graupel)
                  cpm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cp_air + &
&                   q0(i, kk, sphum)*cp_vapor + q_liq*c_liq + q_sol*&
&                   c_ice
                  cvm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cv_air + &
&                   q0(i, kk, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
                END DO
              END IF
              DO i=is,ie
                tv = gz(i, kk) + 0.5*(u0(i, kk)**2+v0(i, kk)**2+w0(i, kk&
&                 )**2)
                t0(i, kk) = (te(i, kk)-tv)/cvm(i)
                hd(i, kk) = cpm(i)*t0(i, kk) + tv
              END DO
            END DO
          END IF
        END DO
      END DO
! k-loop
! n-loop
!--------------------
      IF (fra .LT. 1.) THEN
        DO k=1,kbot
          DO i=is,ie
            t0(i, k) = ta(i, j, k) + (t0(i, k)-ta(i, j, k))*fra
            u0(i, k) = ua(i, j, k) + (u0(i, k)-ua(i, j, k))*fra
            v0(i, k) = va(i, j, k) + (v0(i, k)-va(i, j, k))*fra
          END DO
        END DO
        IF (.NOT.hydrostatic) THEN
          DO k=1,kbot
            DO i=is,ie
              w0(i, k) = w(i, j, k) + (w0(i, k)-w(i, j, k))*fra
            END DO
          END DO
        END IF
        DO iq=1,nq
          DO k=1,kbot
            DO i=is,ie
              q0(i, k, iq) = qa(i, j, k, iq) + (q0(i, k, iq)-qa(i, j, k&
&               , iq))*fra
            END DO
          END DO
        END DO
      END IF
      DO k=1,kbot
        DO i=is,ie
          u_dt(i, j, k) = rdt*(u0(i, k)-ua(i, j, k))
          v_dt(i, j, k) = rdt*(v0(i, k)-va(i, j, k))
! *** temperature updated ***
          ta(i, j, k) = t0(i, k)
          ua(i, j, k) = u0(i, k)
          va(i, j, k) = v0(i, k)
        END DO
        DO iq=1,nq
          DO i=is,ie
            qa(i, j, k, iq) = q0(i, k, iq)
          END DO
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        DO k=1,kbot
          DO i=is,ie
! w updated
            w(i, j, k) = w0(i, k)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE FV_SUBGRID_Z
  SUBROUTINE NEG_ADJ3(is, ie, js, je, ng, kbot, hydrostatic, peln, delz&
&   , pt, dp, qv, ql, qr, qi, qs, qg, qa, check_negative)
    IMPLICIT NONE
!Stubbed, would require careful checking of nonlinearity
! This is designed for 6-class micro-physics schemes
    INTEGER, INTENT(IN) :: is, ie, js, je, ng, kbot
    LOGICAL, INTENT(IN) :: hydrostatic
! total delp-p
    REAL, INTENT(IN) :: dp(is-ng:ie+ng, js-ng:je+ng, kbot)
    REAL, INTENT(IN) :: delz(is-ng:, js-ng:, :)
! ln(pe)
    REAL, INTENT(IN) :: peln(is:ie, kbot+1, js:je)
    LOGICAL, INTENT(IN), OPTIONAL :: check_negative
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, kbot), INTENT(INOUT) :: pt&
&   , qv, ql, qr, qi, qs, qg
    REAL, DIMENSION(is-ng:ie+ng, js-ng:je+ng, kbot), INTENT(INOUT), &
&   OPTIONAL :: qa
  END SUBROUTINE NEG_ADJ3

end module fv_sg_tlm_mod
