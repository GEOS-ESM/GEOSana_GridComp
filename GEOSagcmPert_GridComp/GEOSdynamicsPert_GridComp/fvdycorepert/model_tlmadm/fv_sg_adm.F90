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
module fv_sg_adm_mod

!-----------------------------------------------------------------------
! FV sub-grid mixing
!-----------------------------------------------------------------------
  use constants_mod,      only: rdgas, rvgas, cp_air, cp_vapor, hlv, hlf, kappa, grav
  use tracer_manager_mod, only: get_tracer_index
  use field_manager_mod,  only: MODEL_ATMOS
  use lin_cld_microphys_mod, only: wqs2, wqsat2_moist
  use fv_mp_mod,          only: mp_reduce_min, is_master
  use tapenade_iter,      only: pushcontrol, popcontrol, pushinteger, popinteger, &
                                pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

implicit none
private

public  fv_subgrid_z_fwd, fv_subgrid_z_bwd
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
!  Differentiation of fv_subgrid_z in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge
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
!   gradient     of useful results: qa w delp delz pkz ta
!   with respect to varying inputs: qa peln w delp ua delz va pkz
!                pe ta
  SUBROUTINE FV_SUBGRID_Z_FWD(isd, ied, jsd, jed, is, ie, js, je, km, nq&
&   , dt, tau, nwat, delp, pe, peln, pkz, ta, qa, ua, va, hydrostatic, w&
&   , delz, u_dt, v_dt, t_dt, k_bot)
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
    INTEGER :: ad_from
    REAL :: y1
! = rdgas * (7/2-1) = 2.5*rdgas=717.68
    cv_air = cp_air - rdgas
    rk = cp_air/rdgas + 1.
    g2 = 0.5*grav
    IF (PRESENT(k_bot)) THEN
      IF (k_bot .LT. 3) THEN
        CALL PUSHCONTROL(1,0)
        GOTO 100
      ELSE
        CALL PUSHCONTROL(1,0)
        kbot = k_bot
      END IF
    ELSE
      CALL PUSHCONTROL(1,1)
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
      CALL PUSHCONTROL(1,0)
      xvir = 0.
      rz = 0.
    ELSE
      xvir = zvir
! rz = zvir * rdgas
      rz = rvgas - rdgas
      IF (nwat .EQ. 3) THEN
        CALL PUSHCONTROL(1,1)
        liq_wat = 2
        ice_wat = 3
      ELSE
        CALL PUSHCONTROL(1,1)
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
            CALL PUSHREALARRAY(q0(i, k, iq))
            q0(i, k, iq) = qa(i, j, k, iq)
          END DO
        END DO
      END DO
      DO k=1,kbot
        DO i=is,ie
          CALL PUSHREALARRAY(t0(i, k))
          t0(i, k) = ta(i, j, k)
          tvm(i, k) = t0(i, k)*(1.+xvir*q0(i, k, sphum))
          CALL PUSHREALARRAY(u0(i, k))
          u0(i, k) = ua(i, j, k)
          CALL PUSHREALARRAY(v0(i, k))
          v0(i, k) = va(i, j, k)
          CALL PUSHREALARRAY(pm(i, k))
          pm(i, k) = delp(i, j, k)/(peln(i, k+1, j)-peln(i, k, j))
        END DO
      END DO
      DO i=is,ie
        CALL PUSHREALARRAY(gzh(i))
        gzh(i) = 0.
      END DO
      IF (hydrostatic) THEN
        DO k=kbot,1,-1
          DO i=is,ie
            CALL PUSHREALARRAY(tv)
            tv = rdgas*tvm(i, k)
            CALL PUSHREALARRAY(gz(i, k))
            gz(i, k) = gzh(i) + tv*(1.-pe(i, k, j)/pm(i, k))
            CALL PUSHREALARRAY(hd(i, k))
            hd(i, k) = cp_air*tvm(i, k) + gz(i, k) + 0.5*(u0(i, k)**2+v0&
&             (i, k)**2)
            CALL PUSHREALARRAY(gzh(i))
            gzh(i) = gzh(i) + tv*(peln(i, k+1, j)-peln(i, k, j))
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        DO k=kbot,1,-1
          IF (nwat .EQ. 0) THEN
            DO i=is,ie
              CALL PUSHREALARRAY(cpm(i))
              cpm(i) = cp_air
              CALL PUSHREALARRAY(cvm(i))
              cvm(i) = cv_air
            END DO
            CALL PUSHCONTROL(3,5)
          ELSE IF (nwat .EQ. 1) THEN
            DO i=is,ie
              CALL PUSHREALARRAY(cpm(i))
              cpm(i) = (1.-q0(i, k, sphum))*cp_air + q0(i, k, sphum)*&
&               cp_vapor
              CALL PUSHREALARRAY(cvm(i))
              cvm(i) = (1.-q0(i, k, sphum))*cv_air + q0(i, k, sphum)*&
&               cv_vap
            END DO
            CALL PUSHCONTROL(3,4)
          ELSE IF (nwat .EQ. 2) THEN
! GFS
            DO i=is,ie
              CALL PUSHREALARRAY(cpm(i))
              cpm(i) = (1.-q0(i, k, sphum))*cp_air + q0(i, k, sphum)*&
&               cp_vapor
              CALL PUSHREALARRAY(cvm(i))
              cvm(i) = (1.-q0(i, k, sphum))*cv_air + q0(i, k, sphum)*&
&               cv_vap
            END DO
            CALL PUSHCONTROL(3,3)
          ELSE IF (nwat .EQ. 3) THEN
            DO i=is,ie
              q_liq = q0(i, k, liq_wat)
              q_sol = q0(i, k, ice_wat)
              CALL PUSHREALARRAY(cpm(i))
              cpm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cp_air + q0(i&
&               , k, sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
              CALL PUSHREALARRAY(cvm(i))
              cvm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cv_air + q0(i&
&               , k, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
            END DO
            CALL PUSHCONTROL(3,2)
          ELSE IF (nwat .EQ. 4) THEN
            DO i=is,ie
              q_liq = q0(i, k, liq_wat) + q0(i, k, rainwat)
              CALL PUSHREALARRAY(cpm(i))
              cpm(i) = (1.-(q0(i, k, sphum)+q_liq))*cp_air + q0(i, k, &
&               sphum)*cp_vapor + q_liq*c_liq
              CALL PUSHREALARRAY(cvm(i))
              cvm(i) = (1.-(q0(i, k, sphum)+q_liq))*cv_air + q0(i, k, &
&               sphum)*cv_vap + q_liq*c_liq
            END DO
            CALL PUSHCONTROL(3,1)
          ELSE
            DO i=is,ie
              q_liq = q0(i, k, liq_wat) + q0(i, k, rainwat)
              q_sol = q0(i, k, ice_wat) + q0(i, k, snowwat) + q0(i, k, &
&               graupel)
              CALL PUSHREALARRAY(cpm(i))
              cpm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cp_air + q0(i&
&               , k, sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
              CALL PUSHREALARRAY(cvm(i))
              cvm(i) = (1.-(q0(i, k, sphum)+q_liq+q_sol))*cv_air + q0(i&
&               , k, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
            END DO
            CALL PUSHCONTROL(3,0)
          END IF
          DO i=is,ie
            CALL PUSHREALARRAY(w0(i, k))
            w0(i, k) = w(i, j, k)
            CALL PUSHREALARRAY(gz(i, k))
            gz(i, k) = gzh(i) - g2*delz(i, j, k)
            tmp = gz(i, k) + 0.5*(u0(i, k)**2+v0(i, k)**2+w0(i, k)**2)
            CALL PUSHREALARRAY(hd(i, k))
            hd(i, k) = cpm(i)*t0(i, k) + tmp
            CALL PUSHREALARRAY(te(i, k))
            te(i, k) = cvm(i)*t0(i, k) + tmp
            CALL PUSHREALARRAY(gzh(i))
            gzh(i) = gzh(i) - grav*delz(i, j, k)
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      END IF
      DO n=1,m
        IF (m .EQ. 3) THEN
          IF (n .EQ. 1) THEN
            CALL PUSHREALARRAY(ratio)
            ratio = 0.25
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
          IF (n .EQ. 2) THEN
            CALL PUSHREALARRAY(ratio)
            ratio = 0.5
            CALL PUSHCONTROL(1,0)
          ELSE
            CALL PUSHCONTROL(1,1)
          END IF
          IF (n .EQ. 3) THEN
            CALL PUSHREALARRAY(ratio)
            ratio = 0.999
            CALL PUSHCONTROL(2,2)
          ELSE
            CALL PUSHCONTROL(2,1)
          END IF
        ELSE
          CALL PUSHREALARRAY(ratio)
          ratio = REAL(n)/REAL(m)
          CALL PUSHCONTROL(2,0)
        END IF
        DO i=is,ie
          CALL PUSHREALARRAY(gzh(i))
          gzh(i) = 0.
        END DO
! Compute total condensate
        IF (nwat .LT. 2) THEN
          DO k=1,kbot
            DO i=is,ie
              CALL PUSHREALARRAY(qcon(i, k))
              qcon(i, k) = 0.
            END DO
          END DO
          CALL PUSHCONTROL(3,4)
        ELSE IF (nwat .EQ. 2) THEN
! GFS_2015
          DO k=1,kbot
            DO i=is,ie
              CALL PUSHREALARRAY(qcon(i, k))
              qcon(i, k) = q0(i, k, liq_wat)
            END DO
          END DO
          CALL PUSHCONTROL(3,3)
        ELSE IF (nwat .EQ. 3) THEN
          DO k=1,kbot
            DO i=is,ie
              CALL PUSHREALARRAY(qcon(i, k))
              qcon(i, k) = q0(i, k, liq_wat) + q0(i, k, ice_wat)
            END DO
          END DO
          CALL PUSHCONTROL(3,2)
        ELSE IF (nwat .EQ. 4) THEN
          DO k=1,kbot
            DO i=is,ie
              CALL PUSHREALARRAY(qcon(i, k))
              qcon(i, k) = q0(i, k, liq_wat) + q0(i, k, rainwat)
            END DO
          END DO
          CALL PUSHCONTROL(3,1)
        ELSE
          DO k=1,kbot
            DO i=is,ie
              CALL PUSHREALARRAY(qcon(i, k))
              qcon(i, k) = q0(i, k, liq_wat) + q0(i, k, ice_wat) + q0(i&
&               , k, snowwat) + q0(i, k, rainwat) + q0(i, k, graupel)
            END DO
          END DO
          CALL PUSHCONTROL(3,0)
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
            CALL PUSHREALARRAY(ri)
            ri = (gz(i, km1)-gz(i, k))*(pt1-pt2)/(0.5*(pt1+pt2)*((u0(i, &
&             km1)-u0(i, k))**2+(v0(i, km1)-v0(i, k))**2+ustar2))
            IF (tv1 .GT. t_max .AND. tv1 .GT. tv2) THEN
! top layer unphysically warm
              ri = 0.
              CALL PUSHCONTROL(2,0)
            ELSE IF (tv2 .LT. t_min) THEN
              IF (ri .GT. 0.2) THEN
                ri = 0.2
                CALL PUSHCONTROL(2,3)
              ELSE
                CALL PUSHCONTROL(2,2)
                ri = ri
              END IF
            ELSE
              CALL PUSHCONTROL(2,1)
            END IF
            IF (400.e2 - pm(i, k) .LT. 0.) THEN
              CALL PUSHCONTROL(1,0)
              max2 = 0.
            ELSE
              max2 = 400.e2 - pm(i, k)
              CALL PUSHCONTROL(1,1)
            END IF
            y1 = ri_min + (ri_max-ri_min)*max2/200.e2
            IF (ri_max .GT. y1) THEN
              CALL PUSHREALARRAY(ri_ref)
              ri_ref = y1
              CALL PUSHCONTROL(1,0)
            ELSE
              CALL PUSHREALARRAY(ri_ref)
              ri_ref = ri_max
              CALL PUSHCONTROL(1,1)
            END IF
! Enhancing mixing at the model top
            IF (k .EQ. 2) THEN
              ri_ref = 4.*ri_ref
              CALL PUSHCONTROL(2,0)
            ELSE IF (k .EQ. 3) THEN
              ri_ref = 2.*ri_ref
              CALL PUSHCONTROL(2,1)
            ELSE IF (k .EQ. 4) THEN
              ri_ref = 1.5*ri_ref
              CALL PUSHCONTROL(2,2)
            ELSE
              CALL PUSHCONTROL(2,3)
            END IF
            IF (ri .LT. ri_ref) THEN
              IF (0.0 .LT. ri/ri_ref) THEN
                CALL PUSHREALARRAY(max1)
                max1 = ri/ri_ref
                CALL PUSHCONTROL(1,0)
              ELSE
                CALL PUSHREALARRAY(max1)
                max1 = 0.0
                CALL PUSHCONTROL(1,1)
              END IF
              CALL PUSHREALARRAY(mc)
              mc = ratio*delp(i, j, km1)*delp(i, j, k)/(delp(i, j, km1)+&
&               delp(i, j, k))*(1.-max1)**2
              DO iq=1,nq
                CALL PUSHREALARRAY(h0)
                h0 = mc*(q0(i, k, iq)-q0(i, km1, iq))
                CALL PUSHREALARRAY(q0(i, km1, iq))
                q0(i, km1, iq) = q0(i, km1, iq) + h0/delp(i, j, km1)
                CALL PUSHREALARRAY(q0(i, k, iq))
                q0(i, k, iq) = q0(i, k, iq) - h0/delp(i, j, k)
              END DO
! Recompute qcon
              IF (nwat .LT. 2) THEN
                CALL PUSHREALARRAY(qcon(i, km1))
                qcon(i, km1) = 0.
                CALL PUSHCONTROL(3,0)
              ELSE IF (nwat .EQ. 2) THEN
! GFS_2015
                CALL PUSHREALARRAY(qcon(i, km1))
                qcon(i, km1) = q0(i, km1, liq_wat)
                CALL PUSHCONTROL(3,1)
              ELSE IF (nwat .EQ. 3) THEN
! AM3/AM4
                CALL PUSHREALARRAY(qcon(i, km1))
                qcon(i, km1) = q0(i, km1, liq_wat) + q0(i, km1, ice_wat)
                CALL PUSHCONTROL(3,2)
              ELSE IF (nwat .EQ. 4) THEN
! K_warm_rain scheme with fake ice
                CALL PUSHREALARRAY(qcon(i, km1))
                qcon(i, km1) = q0(i, km1, liq_wat) + q0(i, km1, rainwat)
                CALL PUSHCONTROL(3,3)
              ELSE
                CALL PUSHREALARRAY(qcon(i, km1))
                qcon(i, km1) = q0(i, km1, liq_wat) + q0(i, km1, ice_wat)&
&                 + q0(i, km1, snowwat) + q0(i, km1, rainwat) + q0(i, &
&                 km1, graupel)
                CALL PUSHCONTROL(3,4)
              END IF
! u:
              CALL PUSHREALARRAY(h0)
              h0 = mc*(u0(i, k)-u0(i, k-1))
              CALL PUSHREALARRAY(u0(i, k-1))
              u0(i, k-1) = u0(i, k-1) + h0/delp(i, j, k-1)
              CALL PUSHREALARRAY(u0(i, k))
              u0(i, k) = u0(i, k) - h0/delp(i, j, k)
! v:
              CALL PUSHREALARRAY(h0)
              h0 = mc*(v0(i, k)-v0(i, k-1))
              CALL PUSHREALARRAY(v0(i, k-1))
              v0(i, k-1) = v0(i, k-1) + h0/delp(i, j, k-1)
              CALL PUSHREALARRAY(v0(i, k))
              v0(i, k) = v0(i, k) - h0/delp(i, j, k)
              IF (hydrostatic) THEN
! Static energy
                CALL PUSHREALARRAY(h0)
                h0 = mc*(hd(i, k)-hd(i, k-1))
                CALL PUSHREALARRAY(hd(i, k-1))
                hd(i, k-1) = hd(i, k-1) + h0/delp(i, j, k-1)
                CALL PUSHREALARRAY(hd(i, k))
                hd(i, k) = hd(i, k) - h0/delp(i, j, k)
                CALL PUSHCONTROL(2,2)
              ELSE
! Total energy
                CALL PUSHREALARRAY(h0)
                h0 = mc*(hd(i, k)-hd(i, k-1))
                CALL PUSHREALARRAY(te(i, k-1))
                te(i, k-1) = te(i, k-1) + h0/delp(i, j, k-1)
                CALL PUSHREALARRAY(te(i, k))
                te(i, k) = te(i, k) - h0/delp(i, j, k)
! w:
                h0 = mc*(w0(i, k)-w0(i, k-1))
                CALL PUSHREALARRAY(w0(i, k-1))
                w0(i, k-1) = w0(i, k-1) + h0/delp(i, j, k-1)
                CALL PUSHREALARRAY(w0(i, k))
                w0(i, k) = w0(i, k) - h0/delp(i, j, k)
                CALL PUSHCONTROL(2,1)
              END IF
            ELSE
              CALL PUSHCONTROL(2,0)
            END IF
          END DO
!--------------
! Retrive Temp:
!--------------
          IF (hydrostatic) THEN
            CALL PUSHINTEGER(kk)
            kk = k
            DO i=is,ie
              CALL PUSHREALARRAY(t0(i, kk))
              t0(i, kk) = (hd(i, kk)-gzh(i)-0.5*(u0(i, kk)**2+v0(i, kk)&
&               **2))/(rk-pe(i, kk, j)/pm(i, kk))
              CALL PUSHREALARRAY(gzh(i))
              gzh(i) = gzh(i) + t0(i, kk)*(peln(i, kk+1, j)-peln(i, kk, &
&               j))
              CALL PUSHREALARRAY(t0(i, kk))
              t0(i, kk) = t0(i, kk)/(rdgas+rz*q0(i, kk, sphum))
            END DO
            kk = k - 1
            DO i=is,ie
              CALL PUSHREALARRAY(t0(i, kk))
              t0(i, kk) = (hd(i, kk)-gzh(i)-0.5*(u0(i, kk)**2+v0(i, kk)&
&               **2))/((rk-pe(i, kk, j)/pm(i, kk))*(rdgas+rz*q0(i, kk, &
&               sphum)))
            END DO
            CALL PUSHCONTROL(1,1)
          ELSE
            ad_from = k - 1
            CALL PUSHINTEGER(kk)
! Non-hydrostatic under constant volume heating/cooling
            DO kk=ad_from,k
              IF (nwat .EQ. 0) THEN
                DO i=is,ie
                  CALL PUSHREALARRAY(cpm(i))
                  cpm(i) = cp_air
                  CALL PUSHREALARRAY(cvm(i))
                  cvm(i) = cv_air
                END DO
                CALL PUSHCONTROL(3,5)
              ELSE IF (nwat .EQ. 1) THEN
                DO i=is,ie
                  CALL PUSHREALARRAY(cpm(i))
                  cpm(i) = (1.-q0(i, kk, sphum))*cp_air + q0(i, kk, &
&                   sphum)*cp_vapor
                  CALL PUSHREALARRAY(cvm(i))
                  cvm(i) = (1.-q0(i, kk, sphum))*cv_air + q0(i, kk, &
&                   sphum)*cv_vap
                END DO
                CALL PUSHCONTROL(3,4)
              ELSE IF (nwat .EQ. 2) THEN
                DO i=is,ie
                  CALL PUSHREALARRAY(cpm(i))
                  cpm(i) = (1.-q0(i, kk, sphum))*cp_air + q0(i, kk, &
&                   sphum)*cp_vapor
                  CALL PUSHREALARRAY(cvm(i))
                  cvm(i) = (1.-q0(i, kk, sphum))*cv_air + q0(i, kk, &
&                   sphum)*cv_vap
                END DO
                CALL PUSHCONTROL(3,3)
              ELSE IF (nwat .EQ. 3) THEN
                DO i=is,ie
                  q_liq = q0(i, kk, liq_wat)
                  q_sol = q0(i, kk, ice_wat)
                  CALL PUSHREALARRAY(cpm(i))
                  cpm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cp_air + &
&                   q0(i, kk, sphum)*cp_vapor + q_liq*c_liq + q_sol*&
&                   c_ice
                  CALL PUSHREALARRAY(cvm(i))
                  cvm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cv_air + &
&                   q0(i, kk, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
                END DO
                CALL PUSHCONTROL(3,2)
              ELSE IF (nwat .EQ. 4) THEN
                DO i=is,ie
                  q_liq = q0(i, kk, liq_wat) + q0(i, kk, rainwat)
                  CALL PUSHREALARRAY(cpm(i))
                  cpm(i) = (1.-(q0(i, kk, sphum)+q_liq))*cp_air + q0(i, &
&                   kk, sphum)*cp_vapor + q_liq*c_liq
                  CALL PUSHREALARRAY(cvm(i))
                  cvm(i) = (1.-(q0(i, kk, sphum)+q_liq))*cv_air + q0(i, &
&                   kk, sphum)*cv_vap + q_liq*c_liq
                END DO
                CALL PUSHCONTROL(3,1)
              ELSE
                DO i=is,ie
                  q_liq = q0(i, kk, liq_wat) + q0(i, kk, rainwat)
                  q_sol = q0(i, kk, ice_wat) + q0(i, kk, snowwat) + q0(i&
&                   , kk, graupel)
                  CALL PUSHREALARRAY(cpm(i))
                  cpm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cp_air + &
&                   q0(i, kk, sphum)*cp_vapor + q_liq*c_liq + q_sol*&
&                   c_ice
                  CALL PUSHREALARRAY(cvm(i))
                  cvm(i) = (1.-(q0(i, kk, sphum)+q_liq+q_sol))*cv_air + &
&                   q0(i, kk, sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
                END DO
                CALL PUSHCONTROL(3,0)
              END IF
              DO i=is,ie
                CALL PUSHREALARRAY(tv)
                tv = gz(i, kk) + 0.5*(u0(i, kk)**2+v0(i, kk)**2+w0(i, kk&
&                 )**2)
                CALL PUSHREALARRAY(t0(i, kk))
                t0(i, kk) = (te(i, kk)-tv)/cvm(i)
                CALL PUSHREALARRAY(hd(i, kk))
                hd(i, kk) = cpm(i)*t0(i, kk) + tv
              END DO
            END DO
            CALL PUSHINTEGER(kk - 1)
            CALL PUSHINTEGER(ad_from)
            CALL PUSHCONTROL(1,0)
          END IF
        END DO
      END DO
! k-loop
! n-loop
!--------------------
      IF (fra .LT. 1.) THEN
        DO k=1,kbot
          DO i=is,ie
            CALL PUSHREALARRAY(t0(i, k))
            t0(i, k) = ta(i, j, k) + (t0(i, k)-ta(i, j, k))*fra
            CALL PUSHREALARRAY(u0(i, k))
            u0(i, k) = ua(i, j, k) + (u0(i, k)-ua(i, j, k))*fra
            CALL PUSHREALARRAY(v0(i, k))
            v0(i, k) = va(i, j, k) + (v0(i, k)-va(i, j, k))*fra
          END DO
        END DO
        IF (.NOT.hydrostatic) THEN
          DO k=1,kbot
            DO i=is,ie
              CALL PUSHREALARRAY(w0(i, k))
              w0(i, k) = w(i, j, k) + (w0(i, k)-w(i, j, k))*fra
            END DO
          END DO
          CALL PUSHCONTROL(1,1)
        ELSE
          CALL PUSHCONTROL(1,0)
        END IF
        DO iq=1,nq
          DO k=1,kbot
            DO i=is,ie
              CALL PUSHREALARRAY(q0(i, k, iq))
              q0(i, k, iq) = qa(i, j, k, iq) + (q0(i, k, iq)-qa(i, j, k&
&               , iq))*fra
            END DO
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
        CALL PUSHCONTROL(1,0)
      END IF
      DO k=1,kbot
        DO i=is,ie
! *** temperature updated ***
          CALL PUSHREALARRAY(ta(i, j, k))
          ta(i, j, k) = t0(i, k)
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        CALL PUSHCONTROL(1,0)
      ELSE
        CALL PUSHCONTROL(1,1)
      END IF
    END DO
    CALL PUSHREALARRAY(te, (ie-is+1)*km)
    CALL PUSHINTEGER(liq_wat)
    CALL PUSHINTEGER(ice_wat)
    CALL PUSHREALARRAY(xvir)
    CALL PUSHINTEGER(sphum)
    CALL PUSHREALARRAY(pm, (ie-is+1)*km)
    CALL PUSHREALARRAY(mc)
    CALL PUSHREALARRAY(u0, (ie-is+1)*km)
    CALL PUSHINTEGER(rainwat)
    CALL PUSHREALARRAY(qcon, (ie-is+1)*km)
    CALL PUSHREALARRAY(h0)
    CALL PUSHINTEGER(kbot)
    CALL PUSHREALARRAY(gzh, ie - is + 1)
    CALL PUSHINTEGER(snowwat)
    CALL PUSHREALARRAY(ratio)
    CALL PUSHREALARRAY(cv_air)
    CALL PUSHREALARRAY(rz)
    CALL PUSHREALARRAY(ri_ref)
    CALL PUSHREALARRAY(fra)
    CALL PUSHREALARRAY(q0, (ie-is+1)*km*nq)
    CALL PUSHREALARRAY(t0, (ie-is+1)*km)
    CALL PUSHREALARRAY(g2)
    CALL PUSHREALARRAY(rk)
    CALL PUSHREALARRAY(ri)
    CALL PUSHREALARRAY(w0, (ie-is+1)*km)
    CALL PUSHREALARRAY(hd, (ie-is+1)*km)
    CALL PUSHINTEGER(kk)
    CALL PUSHREALARRAY(cpm, ie - is + 1)
    CALL PUSHREALARRAY(gz, (ie-is+1)*km)
    CALL PUSHREALARRAY(tv)
    CALL PUSHINTEGER(m)
    CALL PUSHREALARRAY(v0, (ie-is+1)*km)
    CALL PUSHREALARRAY(cvm, ie - is + 1)
    CALL PUSHREALARRAY(max1)
    CALL PUSHINTEGER(graupel)
    CALL PUSHCONTROL(1,1)
 100 CONTINUE
  END SUBROUTINE FV_SUBGRID_Z_FWD
!  Differentiation of fv_subgrid_z in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edg
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
!   gradient     of useful results: qa w delp delz pkz ta
!   with respect to varying inputs: qa peln w delp ua delz va pkz
!                pe ta
  SUBROUTINE FV_SUBGRID_Z_BWD(isd, ied, jsd, jed, is, ie, js, je, km, nq&
&   , dt, tau, nwat, delp, delp_ad, pe, pe_ad, peln, peln_ad, pkz, &
&   pkz_ad, ta, ta_ad, qa, qa_ad, ua, ua_ad, va, va_ad, hydrostatic, w, &
&   w_ad, delz, delz_ad, u_dt, v_dt, t_dt, k_bot)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, ie, js, je, km, nq, nwat
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed
    INTEGER, INTENT(IN) :: tau
    REAL, INTENT(IN) :: dt
    REAL, INTENT(IN) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL :: pe_ad(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(IN) :: peln(is:ie, km+1, js:je)
    REAL :: peln_ad(is:ie, km+1, js:je)
    REAL, INTENT(IN) :: delp(isd:ied, jsd:jed, km)
    REAL :: delp_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delz(isd:, jsd:, :)
    REAL :: delz_ad(isd:, jsd:, :)
    REAL, INTENT(IN) :: pkz(is:ie, js:je, km)
    REAL :: pkz_ad(is:ie, js:je, km)
    LOGICAL, INTENT(IN) :: hydrostatic
    INTEGER, INTENT(IN), OPTIONAL :: k_bot
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: ua_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: va(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: va_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: w(isd:, jsd:, :)
    REAL, INTENT(INOUT) :: w_ad(isd:, jsd:, :)
    REAL, INTENT(INOUT) :: ta(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: ta_ad(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: qa(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: qa_ad(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: u_dt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_dt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: t_dt(is:ie, js:je, km)
    REAL, DIMENSION(is:ie, km) :: u0, v0, w0, t0, hd, te, gz, tvm, pm, &
&   den
    REAL, DIMENSION(is:ie, km) :: u0_ad, v0_ad, w0_ad, t0_ad, hd_ad, &
&   te_ad, gz_ad, tvm_ad, pm_ad
    REAL :: q0(is:ie, km, nq), qcon(is:ie, km)
    REAL :: q0_ad(is:ie, km, nq), qcon_ad(is:ie, km)
    REAL, DIMENSION(is:ie) :: gzh, lcp2, icp2, cvm, cpm, qs
    REAL, DIMENSION(is:ie) :: gzh_ad, cvm_ad, cpm_ad
    REAL :: ri_ref, ri, pt1, pt2, ratio, tv, cv, tmp, q_liq, q_sol
    REAL :: ri_ref_ad, ri_ad, pt1_ad, pt2_ad, tv_ad, tmp_ad, q_liq_ad, &
&   q_sol_ad
    REAL :: tv1, tv2, g2, h0, mc, fra, rk, rz, rdt, tvd, tv_surf
    REAL :: tv1_ad, tv2_ad, h0_ad, mc_ad
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
    REAL :: max1_ad
    REAL :: max2
    REAL :: max2_ad
    REAL :: temp
    REAL :: temp_ad
    REAL :: temp0
    REAL :: temp1
    REAL :: temp_ad0
    REAL :: temp_ad1
    REAL :: temp_ad2
    REAL :: temp_ad3
    REAL :: temp_ad4
    REAL :: temp_ad5
    REAL :: temp2
    REAL :: temp3
    REAL :: temp4
    REAL :: temp5
    REAL :: temp6
    REAL :: temp_ad6
    REAL :: temp_ad7
    REAL :: temp_ad8
    REAL :: temp_ad9
    REAL :: temp_ad10
    REAL :: temp_ad11
    REAL :: temp_ad12
    REAL :: temp_ad13
    REAL :: temp_ad14
    REAL :: temp_ad15
    REAL :: y1_ad
    REAL :: temp7
    REAL :: temp8
    REAL :: temp_ad16
    REAL :: temp_ad17
    REAL :: temp_ad18
    REAL :: temp_ad19
    REAL :: temp_ad20
    REAL :: temp_ad21
    REAL :: temp_ad22
    REAL :: temp_ad23
    REAL :: temp_ad24
    REAL :: temp_ad25
    REAL :: temp_ad26
    REAL :: temp_ad27
    REAL :: temp_ad28
    REAL :: temp_ad29
    REAL :: temp_ad30
    REAL :: temp9
    REAL :: temp10
    REAL :: temp11
    REAL :: temp_ad31
    REAL :: temp_ad32
    REAL :: temp_ad33
    REAL :: temp12
    REAL :: temp13
    REAL :: temp14
    REAL :: temp15
    REAL :: temp_ad34
    REAL :: temp_ad35
    REAL :: temp_ad36
    REAL :: temp_ad37
    REAL :: temp_ad38
    REAL :: temp_ad39
    REAL :: temp_ad40
    REAL :: temp_ad41
    REAL :: temp_ad42
    INTEGER :: branch
    INTEGER :: ad_from
    INTEGER :: ad_to
    REAL :: y1
    CALL POPCONTROL(1,branch)
    IF (branch .EQ. 0) THEN
      peln_ad = 0.0
      ua_ad = 0.0
      va_ad = 0.0
      pe_ad = 0.0
    ELSE
      CALL POPINTEGER(graupel)
      CALL POPREALARRAY(max1)
      CALL POPREALARRAY(cvm, ie - is + 1)
      CALL POPREALARRAY(v0, (ie-is+1)*km)
      CALL POPINTEGER(m)
      CALL POPREALARRAY(tv)
      CALL POPREALARRAY(gz, (ie-is+1)*km)
      CALL POPREALARRAY(cpm, ie - is + 1)
      CALL POPINTEGER(kk)
      CALL POPREALARRAY(hd, (ie-is+1)*km)
      CALL POPREALARRAY(w0, (ie-is+1)*km)
      CALL POPREALARRAY(ri)
      CALL POPREALARRAY(rk)
      CALL POPREALARRAY(g2)
      CALL POPREALARRAY(t0, (ie-is+1)*km)
      CALL POPREALARRAY(q0, (ie-is+1)*km*nq)
      CALL POPREALARRAY(fra)
      CALL POPREALARRAY(ri_ref)
      CALL POPREALARRAY(rz)
      CALL POPREALARRAY(cv_air)
      CALL POPREALARRAY(ratio)
      CALL POPINTEGER(snowwat)
      CALL POPREALARRAY(gzh, ie - is + 1)
      CALL POPINTEGER(kbot)
      CALL POPREALARRAY(h0)
      CALL POPREALARRAY(qcon, (ie-is+1)*km)
      CALL POPINTEGER(rainwat)
      CALL POPREALARRAY(u0, (ie-is+1)*km)
      CALL POPREALARRAY(mc)
      CALL POPREALARRAY(pm, (ie-is+1)*km)
      CALL POPINTEGER(sphum)
      CALL POPREALARRAY(xvir)
      CALL POPINTEGER(ice_wat)
      CALL POPINTEGER(liq_wat)
      CALL POPREALARRAY(te, (ie-is+1)*km)
      peln_ad = 0.0
      ua_ad = 0.0
      va_ad = 0.0
      pe_ad = 0.0
      cvm_ad = 0.0
      v0_ad = 0.0
      gz_ad = 0.0
      cpm_ad = 0.0
      hd_ad = 0.0
      w0_ad = 0.0
      tvm_ad = 0.0
      t0_ad = 0.0
      q0_ad = 0.0
      gzh_ad = 0.0
      qcon_ad = 0.0
      u0_ad = 0.0
      pm_ad = 0.0
      te_ad = 0.0
      DO j=je,js,-1
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO k=kbot,1,-1
            DO i=ie,is,-1
              w0_ad(i, k) = w0_ad(i, k) + w_ad(i, j, k)
              w_ad(i, j, k) = 0.0
            END DO
          END DO
        END IF
        DO k=kbot,1,-1
          DO iq=nq,1,-1
            DO i=ie,is,-1
              q0_ad(i, k, iq) = q0_ad(i, k, iq) + qa_ad(i, j, k, iq)
              qa_ad(i, j, k, iq) = 0.0
            END DO
          END DO
          DO i=ie,is,-1
            v0_ad(i, k) = v0_ad(i, k) + va_ad(i, j, k)
            va_ad(i, j, k) = 0.0
            u0_ad(i, k) = u0_ad(i, k) + ua_ad(i, j, k)
            ua_ad(i, j, k) = 0.0
            CALL POPREALARRAY(ta(i, j, k))
            t0_ad(i, k) = t0_ad(i, k) + ta_ad(i, j, k)
            ta_ad(i, j, k) = 0.0
          END DO
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .NE. 0) THEN
          DO iq=nq,1,-1
            DO k=kbot,1,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(q0(i, k, iq))
                qa_ad(i, j, k, iq) = qa_ad(i, j, k, iq) + (1.0-fra)*&
&                 q0_ad(i, k, iq)
                q0_ad(i, k, iq) = fra*q0_ad(i, k, iq)
              END DO
            END DO
          END DO
          CALL POPCONTROL(1,branch)
          IF (branch .NE. 0) THEN
            DO k=kbot,1,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(w0(i, k))
                w_ad(i, j, k) = w_ad(i, j, k) + (1.0-fra)*w0_ad(i, k)
                w0_ad(i, k) = fra*w0_ad(i, k)
              END DO
            END DO
          END IF
          DO k=kbot,1,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(v0(i, k))
              va_ad(i, j, k) = va_ad(i, j, k) + (1.0-fra)*v0_ad(i, k)
              v0_ad(i, k) = fra*v0_ad(i, k)
              CALL POPREALARRAY(u0(i, k))
              ua_ad(i, j, k) = ua_ad(i, j, k) + (1.0-fra)*u0_ad(i, k)
              u0_ad(i, k) = fra*u0_ad(i, k)
              CALL POPREALARRAY(t0(i, k))
              ta_ad(i, j, k) = ta_ad(i, j, k) + (1.0-fra)*t0_ad(i, k)
              t0_ad(i, k) = fra*t0_ad(i, k)
            END DO
          END DO
        END IF
        DO n=m,1,-1
          DO k=2,kbot,1
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) THEN
              CALL POPINTEGER(ad_from)
              CALL POPINTEGER(ad_to)
              DO kk=ad_to,ad_from,-1
                DO i=ie,is,-1
                  CALL POPREALARRAY(hd(i, kk))
                  cpm_ad(i) = cpm_ad(i) + t0(i, kk)*hd_ad(i, kk)
                  t0_ad(i, kk) = t0_ad(i, kk) + cpm(i)*hd_ad(i, kk)
                  CALL POPREALARRAY(t0(i, kk))
                  temp_ad41 = t0_ad(i, kk)/cvm(i)
                  tv_ad = hd_ad(i, kk) - temp_ad41
                  hd_ad(i, kk) = 0.0
                  te_ad(i, kk) = te_ad(i, kk) + temp_ad41
                  cvm_ad(i) = cvm_ad(i) - (te(i, kk)-tv)*temp_ad41/cvm(i&
&                   )
                  t0_ad(i, kk) = 0.0
                  CALL POPREALARRAY(tv)
                  temp_ad42 = 0.5*tv_ad
                  gz_ad(i, kk) = gz_ad(i, kk) + tv_ad
                  u0_ad(i, kk) = u0_ad(i, kk) + 2*u0(i, kk)*temp_ad42
                  v0_ad(i, kk) = v0_ad(i, kk) + 2*v0(i, kk)*temp_ad42
                  w0_ad(i, kk) = w0_ad(i, kk) + 2*w0(i, kk)*temp_ad42
                END DO
                CALL POPCONTROL(3,branch)
                IF (branch .LT. 3) THEN
                  IF (branch .EQ. 0) THEN
                    DO i=ie,is,-1
                      temp_ad40 = cp_air*cpm_ad(i)
                      CALL POPREALARRAY(cvm(i))
                      temp_ad39 = cv_air*cvm_ad(i)
                      q_liq_ad = c_liq*cpm_ad(i) - temp_ad40 + c_liq*&
&                       cvm_ad(i) - temp_ad39
                      q_sol_ad = c_ice*cpm_ad(i) - temp_ad40 + c_ice*&
&                       cvm_ad(i) - temp_ad39
                      q0_ad(i, kk, sphum) = q0_ad(i, kk, sphum) + &
&                       cp_vapor*cpm_ad(i) - temp_ad40 + cv_vap*cvm_ad(i&
&                       ) - temp_ad39
                      cvm_ad(i) = 0.0
                      CALL POPREALARRAY(cpm(i))
                      cpm_ad(i) = 0.0
                      q0_ad(i, kk, ice_wat) = q0_ad(i, kk, ice_wat) + &
&                       q_sol_ad
                      q0_ad(i, kk, snowwat) = q0_ad(i, kk, snowwat) + &
&                       q_sol_ad
                      q0_ad(i, kk, graupel) = q0_ad(i, kk, graupel) + &
&                       q_sol_ad
                      q0_ad(i, kk, liq_wat) = q0_ad(i, kk, liq_wat) + &
&                       q_liq_ad
                      q0_ad(i, kk, rainwat) = q0_ad(i, kk, rainwat) + &
&                       q_liq_ad
                    END DO
                  ELSE IF (branch .EQ. 1) THEN
                    DO i=ie,is,-1
                      CALL POPREALARRAY(cvm(i))
                      q_liq_ad = (c_liq-cp_air)*cpm_ad(i) + (c_liq-&
&                       cv_air)*cvm_ad(i)
                      q0_ad(i, kk, sphum) = q0_ad(i, kk, sphum) + (&
&                       cp_vapor-cp_air)*cpm_ad(i) + (cv_vap-cv_air)*&
&                       cvm_ad(i)
                      cvm_ad(i) = 0.0
                      CALL POPREALARRAY(cpm(i))
                      cpm_ad(i) = 0.0
                      q0_ad(i, kk, liq_wat) = q0_ad(i, kk, liq_wat) + &
&                       q_liq_ad
                      q0_ad(i, kk, rainwat) = q0_ad(i, kk, rainwat) + &
&                       q_liq_ad
                    END DO
                  ELSE
                    DO i=ie,is,-1
                      temp_ad38 = cp_air*cpm_ad(i)
                      CALL POPREALARRAY(cvm(i))
                      temp_ad37 = cv_air*cvm_ad(i)
                      q_liq_ad = c_liq*cpm_ad(i) - temp_ad38 + c_liq*&
&                       cvm_ad(i) - temp_ad37
                      q_sol_ad = c_ice*cpm_ad(i) - temp_ad38 + c_ice*&
&                       cvm_ad(i) - temp_ad37
                      q0_ad(i, kk, sphum) = q0_ad(i, kk, sphum) + &
&                       cp_vapor*cpm_ad(i) - temp_ad38 + cv_vap*cvm_ad(i&
&                       ) - temp_ad37
                      cvm_ad(i) = 0.0
                      CALL POPREALARRAY(cpm(i))
                      cpm_ad(i) = 0.0
                      q0_ad(i, kk, ice_wat) = q0_ad(i, kk, ice_wat) + &
&                       q_sol_ad
                      q0_ad(i, kk, liq_wat) = q0_ad(i, kk, liq_wat) + &
&                       q_liq_ad
                    END DO
                  END IF
                ELSE IF (branch .EQ. 3) THEN
                  DO i=ie,is,-1
                    CALL POPREALARRAY(cvm(i))
                    q0_ad(i, kk, sphum) = q0_ad(i, kk, sphum) + (&
&                     cp_vapor-cp_air)*cpm_ad(i) + (cv_vap-cv_air)*&
&                     cvm_ad(i)
                    cvm_ad(i) = 0.0
                    CALL POPREALARRAY(cpm(i))
                    cpm_ad(i) = 0.0
                  END DO
                ELSE IF (branch .EQ. 4) THEN
                  DO i=ie,is,-1
                    CALL POPREALARRAY(cvm(i))
                    q0_ad(i, kk, sphum) = q0_ad(i, kk, sphum) + (&
&                     cp_vapor-cp_air)*cpm_ad(i) + (cv_vap-cv_air)*&
&                     cvm_ad(i)
                    cvm_ad(i) = 0.0
                    CALL POPREALARRAY(cpm(i))
                    cpm_ad(i) = 0.0
                  END DO
                ELSE
                  DO i=ie,is,-1
                    CALL POPREALARRAY(cvm(i))
                    cvm_ad(i) = 0.0
                    CALL POPREALARRAY(cpm(i))
                    cpm_ad(i) = 0.0
                  END DO
                END IF
              END DO
              CALL POPINTEGER(kk)
            ELSE
              DO i=ie,is,-1
                CALL POPREALARRAY(t0(i, kk))
                temp15 = rdgas + rz*q0(i, kk, sphum)
                temp14 = pm(i, kk)
                temp13 = pe(i, kk, j)/temp14
                temp12 = (rk-temp13)*temp15
                temp_ad34 = t0_ad(i, kk)/temp12
                temp_ad35 = -((hd(i, kk)-gzh(i)-0.5*(u0(i, kk)**2+v0(i, &
&                 kk)**2))*temp_ad34/temp12)
                temp_ad36 = -(temp15*temp_ad35/temp14)
                hd_ad(i, kk) = hd_ad(i, kk) + temp_ad34
                gzh_ad(i) = gzh_ad(i) - temp_ad34
                u0_ad(i, kk) = u0_ad(i, kk) - 0.5*2*u0(i, kk)*temp_ad34
                v0_ad(i, kk) = v0_ad(i, kk) - 0.5*2*v0(i, kk)*temp_ad34
                pe_ad(i, kk, j) = pe_ad(i, kk, j) + temp_ad36
                pm_ad(i, kk) = pm_ad(i, kk) - temp13*temp_ad36
                q0_ad(i, kk, sphum) = q0_ad(i, kk, sphum) + (rk-temp13)*&
&                 rz*temp_ad35
                t0_ad(i, kk) = 0.0
              END DO
              kk = k
              DO i=ie,is,-1
                CALL POPREALARRAY(t0(i, kk))
                temp11 = rdgas + rz*q0(i, kk, sphum)
                q0_ad(i, kk, sphum) = q0_ad(i, kk, sphum) - t0(i, kk)*rz&
&                 *t0_ad(i, kk)/temp11**2
                t0_ad(i, kk) = (peln(i, kk+1, j)-peln(i, kk, j))*gzh_ad(&
&                 i) + t0_ad(i, kk)/temp11
                CALL POPREALARRAY(gzh(i))
                temp_ad31 = t0(i, kk)*gzh_ad(i)
                peln_ad(i, kk+1, j) = peln_ad(i, kk+1, j) + temp_ad31
                peln_ad(i, kk, j) = peln_ad(i, kk, j) - temp_ad31
                CALL POPREALARRAY(t0(i, kk))
                temp10 = pm(i, kk)
                temp9 = pe(i, kk, j)/temp10
                temp_ad32 = t0_ad(i, kk)/(rk-temp9)
                temp_ad33 = (hd(i, kk)-gzh(i)-0.5*(u0(i, kk)**2+v0(i, kk&
&                 )**2))*temp_ad32/((rk-temp9)*temp10)
                hd_ad(i, kk) = hd_ad(i, kk) + temp_ad32
                gzh_ad(i) = gzh_ad(i) - temp_ad32
                u0_ad(i, kk) = u0_ad(i, kk) - 0.5*2*u0(i, kk)*temp_ad32
                v0_ad(i, kk) = v0_ad(i, kk) - 0.5*2*v0(i, kk)*temp_ad32
                pe_ad(i, kk, j) = pe_ad(i, kk, j) + temp_ad33
                pm_ad(i, kk) = pm_ad(i, kk) - temp9*temp_ad33
                t0_ad(i, kk) = 0.0
              END DO
              CALL POPINTEGER(kk)
            END IF
            km1 = k - 1
            DO i=ie,is,-1
              CALL POPCONTROL(2,branch)
              IF (branch .EQ. 0) THEN
                ri_ad = 0.0
                ri_ref_ad = 0.0
              ELSE
                IF (branch .EQ. 1) THEN
                  temp_ad30 = te_ad(i, k-1)/delp(i, j, k-1)
                  temp_ad28 = w0_ad(i, k-1)/delp(i, j, k-1)
                  CALL POPREALARRAY(w0(i, k))
                  temp_ad27 = -(w0_ad(i, k)/delp(i, j, k))
                  h0_ad = temp_ad28 + temp_ad27
                  delp_ad(i, j, k) = delp_ad(i, j, k) - h0*temp_ad27/&
&                   delp(i, j, k)
                  CALL POPREALARRAY(w0(i, k-1))
                  delp_ad(i, j, k-1) = delp_ad(i, j, k-1) - h0*temp_ad28&
&                   /delp(i, j, k-1)
                  mc_ad = (w0(i, k)-w0(i, k-1))*h0_ad
                  w0_ad(i, k) = w0_ad(i, k) + mc*h0_ad
                  w0_ad(i, k-1) = w0_ad(i, k-1) - mc*h0_ad
                  h0 = mc*(hd(i, k)-hd(i, k-1))
                  CALL POPREALARRAY(te(i, k))
                  temp_ad29 = -(te_ad(i, k)/delp(i, j, k))
                  h0_ad = temp_ad30 + temp_ad29
                  delp_ad(i, j, k) = delp_ad(i, j, k) - h0*temp_ad29/&
&                   delp(i, j, k)
                  CALL POPREALARRAY(te(i, k-1))
                  delp_ad(i, j, k-1) = delp_ad(i, j, k-1) - h0*temp_ad30&
&                   /delp(i, j, k-1)
                  CALL POPREALARRAY(h0)
                  mc_ad = mc_ad + (hd(i, k)-hd(i, k-1))*h0_ad
                  hd_ad(i, k) = hd_ad(i, k) + mc*h0_ad
                  hd_ad(i, k-1) = hd_ad(i, k-1) - mc*h0_ad
                ELSE
                  temp_ad26 = hd_ad(i, k-1)/delp(i, j, k-1)
                  CALL POPREALARRAY(hd(i, k))
                  temp_ad25 = -(hd_ad(i, k)/delp(i, j, k))
                  h0_ad = temp_ad26 + temp_ad25
                  delp_ad(i, j, k) = delp_ad(i, j, k) - h0*temp_ad25/&
&                   delp(i, j, k)
                  CALL POPREALARRAY(hd(i, k-1))
                  delp_ad(i, j, k-1) = delp_ad(i, j, k-1) - h0*temp_ad26&
&                   /delp(i, j, k-1)
                  CALL POPREALARRAY(h0)
                  mc_ad = (hd(i, k)-hd(i, k-1))*h0_ad
                  hd_ad(i, k) = hd_ad(i, k) + mc*h0_ad
                  hd_ad(i, k-1) = hd_ad(i, k-1) - mc*h0_ad
                END IF
                temp_ad24 = u0_ad(i, k-1)/delp(i, j, k-1)
                temp_ad22 = v0_ad(i, k-1)/delp(i, j, k-1)
                CALL POPREALARRAY(v0(i, k))
                temp_ad21 = -(v0_ad(i, k)/delp(i, j, k))
                h0_ad = temp_ad22 + temp_ad21
                delp_ad(i, j, k) = delp_ad(i, j, k) - h0*temp_ad21/delp(&
&                 i, j, k)
                CALL POPREALARRAY(v0(i, k-1))
                delp_ad(i, j, k-1) = delp_ad(i, j, k-1) - h0*temp_ad22/&
&                 delp(i, j, k-1)
                CALL POPREALARRAY(h0)
                mc_ad = mc_ad + (v0(i, k)-v0(i, k-1))*h0_ad
                v0_ad(i, k) = v0_ad(i, k) + mc*h0_ad
                v0_ad(i, k-1) = v0_ad(i, k-1) - mc*h0_ad
                CALL POPREALARRAY(u0(i, k))
                temp_ad23 = -(u0_ad(i, k)/delp(i, j, k))
                h0_ad = temp_ad24 + temp_ad23
                delp_ad(i, j, k) = delp_ad(i, j, k) - h0*temp_ad23/delp(&
&                 i, j, k)
                CALL POPREALARRAY(u0(i, k-1))
                delp_ad(i, j, k-1) = delp_ad(i, j, k-1) - h0*temp_ad24/&
&                 delp(i, j, k-1)
                CALL POPREALARRAY(h0)
                mc_ad = mc_ad + (u0(i, k)-u0(i, k-1))*h0_ad
                u0_ad(i, k) = u0_ad(i, k) + mc*h0_ad
                u0_ad(i, k-1) = u0_ad(i, k-1) - mc*h0_ad
                CALL POPCONTROL(3,branch)
                IF (branch .LT. 2) THEN
                  IF (branch .EQ. 0) THEN
                    CALL POPREALARRAY(qcon(i, km1))
                    qcon_ad(i, km1) = 0.0
                  ELSE
                    CALL POPREALARRAY(qcon(i, km1))
                    q0_ad(i, km1, liq_wat) = q0_ad(i, km1, liq_wat) + &
&                     qcon_ad(i, km1)
                    qcon_ad(i, km1) = 0.0
                  END IF
                ELSE IF (branch .EQ. 2) THEN
                  CALL POPREALARRAY(qcon(i, km1))
                  q0_ad(i, km1, liq_wat) = q0_ad(i, km1, liq_wat) + &
&                   qcon_ad(i, km1)
                  q0_ad(i, km1, ice_wat) = q0_ad(i, km1, ice_wat) + &
&                   qcon_ad(i, km1)
                  qcon_ad(i, km1) = 0.0
                ELSE IF (branch .EQ. 3) THEN
                  CALL POPREALARRAY(qcon(i, km1))
                  q0_ad(i, km1, liq_wat) = q0_ad(i, km1, liq_wat) + &
&                   qcon_ad(i, km1)
                  q0_ad(i, km1, rainwat) = q0_ad(i, km1, rainwat) + &
&                   qcon_ad(i, km1)
                  qcon_ad(i, km1) = 0.0
                ELSE
                  CALL POPREALARRAY(qcon(i, km1))
                  q0_ad(i, km1, liq_wat) = q0_ad(i, km1, liq_wat) + &
&                   qcon_ad(i, km1)
                  q0_ad(i, km1, ice_wat) = q0_ad(i, km1, ice_wat) + &
&                   qcon_ad(i, km1)
                  q0_ad(i, km1, snowwat) = q0_ad(i, km1, snowwat) + &
&                   qcon_ad(i, km1)
                  q0_ad(i, km1, rainwat) = q0_ad(i, km1, rainwat) + &
&                   qcon_ad(i, km1)
                  q0_ad(i, km1, graupel) = q0_ad(i, km1, graupel) + &
&                   qcon_ad(i, km1)
                  qcon_ad(i, km1) = 0.0
                END IF
                DO iq=nq,1,-1
                  temp_ad20 = q0_ad(i, km1, iq)/delp(i, j, km1)
                  CALL POPREALARRAY(q0(i, k, iq))
                  temp_ad19 = -(q0_ad(i, k, iq)/delp(i, j, k))
                  h0_ad = temp_ad20 + temp_ad19
                  delp_ad(i, j, k) = delp_ad(i, j, k) - h0*temp_ad19/&
&                   delp(i, j, k)
                  CALL POPREALARRAY(q0(i, km1, iq))
                  delp_ad(i, j, km1) = delp_ad(i, j, km1) - h0*temp_ad20&
&                   /delp(i, j, km1)
                  CALL POPREALARRAY(h0)
                  mc_ad = mc_ad + (q0(i, k, iq)-q0(i, km1, iq))*h0_ad
                  q0_ad(i, k, iq) = q0_ad(i, k, iq) + mc*h0_ad
                  q0_ad(i, km1, iq) = q0_ad(i, km1, iq) - mc*h0_ad
                END DO
                CALL POPREALARRAY(mc)
                temp8 = delp(i, j, km1) + delp(i, j, k)
                temp7 = delp(i, j, k)/temp8
                temp_ad16 = ratio*mc_ad
                temp_ad17 = delp(i, j, km1)*(1.-max1)**2*temp_ad16/temp8
                temp_ad18 = -(temp7*temp_ad17)
                delp_ad(i, j, km1) = delp_ad(i, j, km1) + temp_ad18 + &
&                 temp7*(1.-max1)**2*temp_ad16
                max1_ad = -(2*(1.-max1)*delp(i, j, km1)*temp7*temp_ad16)
                delp_ad(i, j, k) = delp_ad(i, j, k) + temp_ad18 + &
&                 temp_ad17
                CALL POPCONTROL(1,branch)
                IF (branch .EQ. 0) THEN
                  CALL POPREALARRAY(max1)
                  ri_ad = max1_ad/ri_ref
                  ri_ref_ad = -(ri*max1_ad/ri_ref**2)
                ELSE
                  CALL POPREALARRAY(max1)
                  ri_ad = 0.0
                  ri_ref_ad = 0.0
                END IF
              END IF
              CALL POPCONTROL(2,branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  ri_ref_ad = 4.*ri_ref_ad
                ELSE
                  ri_ref_ad = 2.*ri_ref_ad
                END IF
              ELSE IF (branch .EQ. 2) THEN
                ri_ref_ad = 1.5*ri_ref_ad
              END IF
              CALL POPCONTROL(1,branch)
              IF (branch .EQ. 0) THEN
                CALL POPREALARRAY(ri_ref)
                y1_ad = ri_ref_ad
              ELSE
                CALL POPREALARRAY(ri_ref)
                y1_ad = 0.0
              END IF
              max2_ad = (ri_max-ri_min)*y1_ad/200.e2
              CALL POPCONTROL(1,branch)
              IF (branch .NE. 0) pm_ad(i, k) = pm_ad(i, k) - max2_ad
              CALL POPCONTROL(2,branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) THEN
                  tv2 = t0(i, k)*(1.+xvir*q0(i, k, sphum)-qcon(i, k))
                  ri_ad = 0.0
                  GOTO 100
                END IF
              ELSE IF (branch .NE. 2) THEN
                ri_ad = 0.0
              END IF
              tv2 = t0(i, k)*(1.+xvir*q0(i, k, sphum)-qcon(i, k))
 100          tv1 = t0(i, km1)*(1.+xvir*q0(i, km1, sphum)-qcon(i, km1))
              pt1 = tv1/pkz(i, j, km1)
              pt2 = tv2/pkz(i, j, k)
              CALL POPREALARRAY(ri)
              temp5 = v0(i, km1) - v0(i, k)
              temp4 = u0(i, km1) - u0(i, k)
              temp3 = ustar2 + temp4**2 + temp5**2
              temp2 = 0.5*(pt1+pt2)*temp3
              temp_ad6 = ri_ad/temp2
              temp6 = gz(i, km1) - gz(i, k)
              temp_ad7 = -(temp6*(pt1-pt2)*temp_ad6/temp2)
              temp_ad8 = temp3*0.5*temp_ad7
              temp_ad9 = 0.5*(pt1+pt2)*temp_ad7
              temp_ad10 = 2*temp4*temp_ad9
              temp_ad11 = 2*temp5*temp_ad9
              gz_ad(i, km1) = gz_ad(i, km1) + (pt1-pt2)*temp_ad6
              gz_ad(i, k) = gz_ad(i, k) - (pt1-pt2)*temp_ad6
              pt1_ad = temp_ad8 + temp6*temp_ad6
              pt2_ad = temp_ad8 - temp6*temp_ad6
              u0_ad(i, km1) = u0_ad(i, km1) + temp_ad10
              u0_ad(i, k) = u0_ad(i, k) - temp_ad10
              v0_ad(i, km1) = v0_ad(i, km1) + temp_ad11
              v0_ad(i, k) = v0_ad(i, k) - temp_ad11
              temp_ad12 = pt2_ad/pkz(i, j, k)
              tv2_ad = temp_ad12
              pkz_ad(i, j, k) = pkz_ad(i, j, k) - tv2*temp_ad12/pkz(i, j&
&               , k)
              temp_ad13 = pt1_ad/pkz(i, j, km1)
              tv1_ad = temp_ad13
              pkz_ad(i, j, km1) = pkz_ad(i, j, km1) - tv1*temp_ad13/pkz(&
&               i, j, km1)
              temp_ad14 = t0(i, k)*tv2_ad
              t0_ad(i, k) = t0_ad(i, k) + (xvir*q0(i, k, sphum)-qcon(i, &
&               k)+1.)*tv2_ad
              q0_ad(i, k, sphum) = q0_ad(i, k, sphum) + xvir*temp_ad14
              qcon_ad(i, k) = qcon_ad(i, k) - temp_ad14
              temp_ad15 = t0(i, km1)*tv1_ad
              t0_ad(i, km1) = t0_ad(i, km1) + (xvir*q0(i, km1, sphum)-&
&               qcon(i, km1)+1.)*tv1_ad
              q0_ad(i, km1, sphum) = q0_ad(i, km1, sphum) + xvir*&
&               temp_ad15
              qcon_ad(i, km1) = qcon_ad(i, km1) - temp_ad15
            END DO
          END DO
          CALL POPCONTROL(3,branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              DO k=kbot,1,-1
                DO i=ie,is,-1
                  CALL POPREALARRAY(qcon(i, k))
                  q0_ad(i, k, liq_wat) = q0_ad(i, k, liq_wat) + qcon_ad(&
&                   i, k)
                  q0_ad(i, k, ice_wat) = q0_ad(i, k, ice_wat) + qcon_ad(&
&                   i, k)
                  q0_ad(i, k, snowwat) = q0_ad(i, k, snowwat) + qcon_ad(&
&                   i, k)
                  q0_ad(i, k, rainwat) = q0_ad(i, k, rainwat) + qcon_ad(&
&                   i, k)
                  q0_ad(i, k, graupel) = q0_ad(i, k, graupel) + qcon_ad(&
&                   i, k)
                  qcon_ad(i, k) = 0.0
                END DO
              END DO
            ELSE
              DO k=kbot,1,-1
                DO i=ie,is,-1
                  CALL POPREALARRAY(qcon(i, k))
                  q0_ad(i, k, liq_wat) = q0_ad(i, k, liq_wat) + qcon_ad(&
&                   i, k)
                  q0_ad(i, k, rainwat) = q0_ad(i, k, rainwat) + qcon_ad(&
&                   i, k)
                  qcon_ad(i, k) = 0.0
                END DO
              END DO
            END IF
          ELSE IF (branch .EQ. 2) THEN
            DO k=kbot,1,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(qcon(i, k))
                q0_ad(i, k, liq_wat) = q0_ad(i, k, liq_wat) + qcon_ad(i&
&                 , k)
                q0_ad(i, k, ice_wat) = q0_ad(i, k, ice_wat) + qcon_ad(i&
&                 , k)
                qcon_ad(i, k) = 0.0
              END DO
            END DO
          ELSE IF (branch .EQ. 3) THEN
            DO k=kbot,1,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(qcon(i, k))
                q0_ad(i, k, liq_wat) = q0_ad(i, k, liq_wat) + qcon_ad(i&
&                 , k)
                qcon_ad(i, k) = 0.0
              END DO
            END DO
          ELSE
            DO k=kbot,1,-1
              DO i=ie,is,-1
                CALL POPREALARRAY(qcon(i, k))
                qcon_ad(i, k) = 0.0
              END DO
            END DO
          END IF
          DO i=ie,is,-1
            CALL POPREALARRAY(gzh(i))
            gzh_ad(i) = 0.0
          END DO
          CALL POPCONTROL(2,branch)
          IF (branch .EQ. 0) THEN
            CALL POPREALARRAY(ratio)
          ELSE
            IF (branch .NE. 1) CALL POPREALARRAY(ratio)
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) CALL POPREALARRAY(ratio)
            CALL POPCONTROL(1,branch)
            IF (branch .EQ. 0) CALL POPREALARRAY(ratio)
          END IF
        END DO
        CALL POPCONTROL(1,branch)
        IF (branch .EQ. 0) THEN
          DO k=1,kbot,1
            DO i=ie,is,-1
              tmp_ad = hd_ad(i, k) + te_ad(i, k)
              gz_ad(i, k) = gz_ad(i, k) + tmp_ad
              CALL POPREALARRAY(gzh(i))
              delz_ad(i, j, k) = delz_ad(i, j, k) - g2*gz_ad(i, k) - &
&               grav*gzh_ad(i)
              CALL POPREALARRAY(te(i, k))
              cvm_ad(i) = cvm_ad(i) + t0(i, k)*te_ad(i, k)
              t0_ad(i, k) = t0_ad(i, k) + cpm(i)*hd_ad(i, k) + cvm(i)*&
&               te_ad(i, k)
              te_ad(i, k) = 0.0
              CALL POPREALARRAY(hd(i, k))
              cpm_ad(i) = cpm_ad(i) + t0(i, k)*hd_ad(i, k)
              hd_ad(i, k) = 0.0
              temp_ad5 = 0.5*tmp_ad
              u0_ad(i, k) = u0_ad(i, k) + 2*u0(i, k)*temp_ad5
              v0_ad(i, k) = v0_ad(i, k) + 2*v0(i, k)*temp_ad5
              w0_ad(i, k) = w0_ad(i, k) + 2*w0(i, k)*temp_ad5
              CALL POPREALARRAY(gz(i, k))
              gzh_ad(i) = gzh_ad(i) + gz_ad(i, k)
              gz_ad(i, k) = 0.0
              CALL POPREALARRAY(w0(i, k))
              w_ad(i, j, k) = w_ad(i, j, k) + w0_ad(i, k)
              w0_ad(i, k) = 0.0
            END DO
            CALL POPCONTROL(3,branch)
            IF (branch .LT. 3) THEN
              IF (branch .EQ. 0) THEN
                DO i=ie,is,-1
                  temp_ad4 = cp_air*cpm_ad(i)
                  CALL POPREALARRAY(cvm(i))
                  temp_ad3 = cv_air*cvm_ad(i)
                  q_liq_ad = c_liq*cpm_ad(i) - temp_ad4 + c_liq*cvm_ad(i&
&                   ) - temp_ad3
                  q_sol_ad = c_ice*cpm_ad(i) - temp_ad4 + c_ice*cvm_ad(i&
&                   ) - temp_ad3
                  q0_ad(i, k, sphum) = q0_ad(i, k, sphum) + cp_vapor*&
&                   cpm_ad(i) - temp_ad4 + cv_vap*cvm_ad(i) - temp_ad3
                  cvm_ad(i) = 0.0
                  CALL POPREALARRAY(cpm(i))
                  cpm_ad(i) = 0.0
                  q0_ad(i, k, ice_wat) = q0_ad(i, k, ice_wat) + q_sol_ad
                  q0_ad(i, k, snowwat) = q0_ad(i, k, snowwat) + q_sol_ad
                  q0_ad(i, k, graupel) = q0_ad(i, k, graupel) + q_sol_ad
                  q0_ad(i, k, liq_wat) = q0_ad(i, k, liq_wat) + q_liq_ad
                  q0_ad(i, k, rainwat) = q0_ad(i, k, rainwat) + q_liq_ad
                END DO
              ELSE IF (branch .EQ. 1) THEN
                DO i=ie,is,-1
                  CALL POPREALARRAY(cvm(i))
                  q_liq_ad = (c_liq-cp_air)*cpm_ad(i) + (c_liq-cv_air)*&
&                   cvm_ad(i)
                  q0_ad(i, k, sphum) = q0_ad(i, k, sphum) + (cp_vapor-&
&                   cp_air)*cpm_ad(i) + (cv_vap-cv_air)*cvm_ad(i)
                  cvm_ad(i) = 0.0
                  CALL POPREALARRAY(cpm(i))
                  cpm_ad(i) = 0.0
                  q0_ad(i, k, liq_wat) = q0_ad(i, k, liq_wat) + q_liq_ad
                  q0_ad(i, k, rainwat) = q0_ad(i, k, rainwat) + q_liq_ad
                END DO
              ELSE
                DO i=ie,is,-1
                  temp_ad2 = cp_air*cpm_ad(i)
                  CALL POPREALARRAY(cvm(i))
                  temp_ad1 = cv_air*cvm_ad(i)
                  q_liq_ad = c_liq*cpm_ad(i) - temp_ad2 + c_liq*cvm_ad(i&
&                   ) - temp_ad1
                  q_sol_ad = c_ice*cpm_ad(i) - temp_ad2 + c_ice*cvm_ad(i&
&                   ) - temp_ad1
                  q0_ad(i, k, sphum) = q0_ad(i, k, sphum) + cp_vapor*&
&                   cpm_ad(i) - temp_ad2 + cv_vap*cvm_ad(i) - temp_ad1
                  cvm_ad(i) = 0.0
                  CALL POPREALARRAY(cpm(i))
                  cpm_ad(i) = 0.0
                  q0_ad(i, k, ice_wat) = q0_ad(i, k, ice_wat) + q_sol_ad
                  q0_ad(i, k, liq_wat) = q0_ad(i, k, liq_wat) + q_liq_ad
                END DO
              END IF
            ELSE IF (branch .EQ. 3) THEN
              DO i=ie,is,-1
                CALL POPREALARRAY(cvm(i))
                q0_ad(i, k, sphum) = q0_ad(i, k, sphum) + (cp_vapor-&
&                 cp_air)*cpm_ad(i) + (cv_vap-cv_air)*cvm_ad(i)
                cvm_ad(i) = 0.0
                CALL POPREALARRAY(cpm(i))
                cpm_ad(i) = 0.0
              END DO
            ELSE IF (branch .EQ. 4) THEN
              DO i=ie,is,-1
                CALL POPREALARRAY(cvm(i))
                q0_ad(i, k, sphum) = q0_ad(i, k, sphum) + (cp_vapor-&
&                 cp_air)*cpm_ad(i) + (cv_vap-cv_air)*cvm_ad(i)
                cvm_ad(i) = 0.0
                CALL POPREALARRAY(cpm(i))
                cpm_ad(i) = 0.0
              END DO
            ELSE
              DO i=ie,is,-1
                CALL POPREALARRAY(cvm(i))
                cvm_ad(i) = 0.0
                CALL POPREALARRAY(cpm(i))
                cpm_ad(i) = 0.0
              END DO
            END IF
          END DO
        ELSE
          DO k=1,kbot,1
            DO i=ie,is,-1
              temp1 = pm(i, k)
              temp0 = pe(i, k, j)/temp1
              gz_ad(i, k) = gz_ad(i, k) + hd_ad(i, k)
              CALL POPREALARRAY(gzh(i))
              tv_ad = (1.-temp0)*gz_ad(i, k) + (peln(i, k+1, j)-peln(i, &
&               k, j))*gzh_ad(i)
              peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + tv*gzh_ad(i)
              peln_ad(i, k, j) = peln_ad(i, k, j) - tv*gzh_ad(i)
              CALL POPREALARRAY(hd(i, k))
              tvm_ad(i, k) = tvm_ad(i, k) + rdgas*tv_ad + cp_air*hd_ad(i&
&               , k)
              u0_ad(i, k) = u0_ad(i, k) + 0.5*2*u0(i, k)*hd_ad(i, k)
              v0_ad(i, k) = v0_ad(i, k) + 0.5*2*v0(i, k)*hd_ad(i, k)
              hd_ad(i, k) = 0.0
              CALL POPREALARRAY(gz(i, k))
              temp_ad0 = -(tv*gz_ad(i, k)/temp1)
              gzh_ad(i) = gzh_ad(i) + gz_ad(i, k)
              pe_ad(i, k, j) = pe_ad(i, k, j) + temp_ad0
              pm_ad(i, k) = pm_ad(i, k) - temp0*temp_ad0
              gz_ad(i, k) = 0.0
              CALL POPREALARRAY(tv)
            END DO
          END DO
        END IF
        DO i=ie,is,-1
          CALL POPREALARRAY(gzh(i))
          gzh_ad(i) = 0.0
        END DO
        DO k=kbot,1,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(pm(i, k))
            temp = peln(i, k+1, j) - peln(i, k, j)
            temp_ad = -(delp(i, j, k)*pm_ad(i, k)/temp**2)
            delp_ad(i, j, k) = delp_ad(i, j, k) + pm_ad(i, k)/temp
            peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad
            peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad
            pm_ad(i, k) = 0.0
            CALL POPREALARRAY(v0(i, k))
            va_ad(i, j, k) = va_ad(i, j, k) + v0_ad(i, k)
            v0_ad(i, k) = 0.0
            CALL POPREALARRAY(u0(i, k))
            ua_ad(i, j, k) = ua_ad(i, j, k) + u0_ad(i, k)
            u0_ad(i, k) = 0.0
            t0_ad(i, k) = t0_ad(i, k) + (xvir*q0(i, k, sphum)+1.)*tvm_ad&
&             (i, k)
            q0_ad(i, k, sphum) = q0_ad(i, k, sphum) + t0(i, k)*xvir*&
&             tvm_ad(i, k)
            tvm_ad(i, k) = 0.0
            CALL POPREALARRAY(t0(i, k))
            ta_ad(i, j, k) = ta_ad(i, j, k) + t0_ad(i, k)
            t0_ad(i, k) = 0.0
          END DO
        END DO
        DO iq=nq,1,-1
          DO k=kbot,1,-1
            DO i=ie,is,-1
              CALL POPREALARRAY(q0(i, k, iq))
              qa_ad(i, j, k, iq) = qa_ad(i, j, k, iq) + q0_ad(i, k, iq)
              q0_ad(i, k, iq) = 0.0
            END DO
          END DO
        END DO
      END DO
      CALL POPCONTROL(1,branch)
      CALL POPCONTROL(1,branch)
    END IF
  END SUBROUTINE FV_SUBGRID_Z_BWD
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

end module fv_sg_adm_mod
