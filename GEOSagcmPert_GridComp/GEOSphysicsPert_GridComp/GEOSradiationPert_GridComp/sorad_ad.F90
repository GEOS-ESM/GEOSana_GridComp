module SORAD_AD

USE SORADMOD

IMPLICIT NONE

PRIVATE
PUBLIC :: sorad_b

contains

SUBROUTINE SORAD_B(m, np, nb, cosz_dev, pl_dev, ta_dev, ta_devb, wa_dev&
& , wa_devb, oa_dev, oa_devb, co2, cwc_dev, cwc_devb, fcld_dev, &
& fcld_devb, ict, icb, reff_dev, reff_devb, hk_uv, hk_ir, taua_dev, &
& taua_devb, ssaa_dev, ssaa_devb, asya_dev, asya_devb, rsuvbm_dev, &
& rsuvdf_dev, rsirbm_dev, rsirdf_dev, flx_dev, flx_devb, cons_grav, &
& wk_uv, zk_uv, ry_uv, xk_ir, ry_ir, cah, coa, aig_uv, awg_uv, arg_uv, &
& aib_uv, awb_uv, arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, awa_nir, &
& ara_nir, aig_nir, awg_nir, arg_nir, caib, caif)
  IMPLICIT NONE
!end do RUN_LOOP
! Parameters
! ----------
  INTEGER, PARAMETER :: nu=43
  INTEGER, PARAMETER :: nw=37
  INTEGER, PARAMETER :: nx=62
  INTEGER, PARAMETER :: ny=101
  INTEGER, PARAMETER :: nband_uv=5
  INTEGER, PARAMETER :: nk_ir=10
  INTEGER, PARAMETER :: nband_ir=3
  INTEGER, PARAMETER :: nband=nband_uv+nband_ir
  REAL*8, PARAMETER :: dsm=0.602
!-----input values
!-----input parameters
  INTEGER :: m, np, ict, icb, nb
  REAL*8 :: cosz_dev(m), pl_dev(m, np+1), ta_dev(m, np), wa_dev(m, np), &
& oa_dev(m, np), co2
  REAL*8 :: ta_devb(m, np), wa_devb(m, np), oa_devb(m, np)
  REAL*8 :: cwc_dev(m, np, 4), fcld_dev(m, np), reff_dev(m, np, 4), &
& hk_uv(5), hk_ir(3, 10)
  REAL*8 :: cwc_devb(m, np, 4), fcld_devb(m, np), reff_devb(m, np, 4)
  REAL*8 :: rsuvbm_dev, rsuvdf_dev, rsirbm_dev, rsirdf_dev
  REAL*8 :: taua_dev(m, np, nb)
  REAL*8 :: taua_devb(m, np, nb)
  REAL*8 :: ssaa_dev(m, np, nb)
  REAL*8 :: ssaa_devb(m, np, nb)
  REAL*8 :: asya_dev(m, np, nb)
  REAL*8 :: asya_devb(m, np, nb)
  LOGICAL :: overcast
! Constants
  REAL*8, INTENT(IN) :: wk_uv(5), zk_uv(5), ry_uv(5)
  REAL*8, INTENT(IN) :: xk_ir(10), ry_ir(3)
  REAL*8, INTENT(IN) :: cah(43, 37), coa(62, 101)
  REAL*8, INTENT(IN) :: aig_uv(3), awg_uv(3), arg_uv(3)
  REAL*8, INTENT(IN) :: aib_uv, awb_uv(2), arb_uv(2)
  REAL*8, INTENT(IN) :: aib_nir, awb_nir(3, 2), arb_nir(3, 2)
  REAL*8, INTENT(IN) :: aia_nir(3, 3), awa_nir(3, 3), ara_nir(3, 3)
  REAL*8, INTENT(IN) :: aig_nir(3, 3), awg_nir(3, 3), arg_nir(3, 3)
  REAL*8, INTENT(IN) :: caib(11, 9, 11), caif(9, 11)
  REAL*8, INTENT(IN) :: cons_grav
!-----output parameters
  REAL*8 :: flx_dev(m, np+1), flc_dev(m, np+1)
  REAL*8 :: flx_devb(m, np+1)
  REAL*8 :: flxu_dev(m, np+1), flcu_dev(m, np+1)
  REAL*8 :: fdiruv_dev(m), fdifuv_dev(m)
  REAL*8 :: fdirpar_dev(m), fdifpar_dev(m)
  REAL*8 :: fdirir_dev(m), fdifir_dev(m)
  REAL*8 :: flx_sfc_band_dev(m, nband)
!-----temporary arrays
  INTEGER :: i, j, k, l, in, ntop
  REAL*8 :: dp(np), wh(np), oh(np)
  REAL*8 :: whb(np), ohb(np)
  REAL*8 :: scal(np)
  REAL*8 :: swh(np+1), so2(np+1), df(0:np+1)
  REAL*8 :: swhb(np+1), dfb(0:np+1)
  REAL*8 :: scal0, wvtoa, o3toa, pa
  REAL*8 :: wvtoab, o3toab
  REAL*8 :: snt, cnt, x, xx4, xtoa
  REAL*8 :: xx4b
  REAL*8 :: dp_pa(np)
!-----parameters for co2 transmission tables
  REAL*8 :: w1, dw, u1, du
  INTEGER :: ib, rc
  REAL*8 :: tauclb(np), tauclf(np), asycl(np)
  REAL*8 :: tauclbb(np), tauclfb(np), asyclb(np)
  REAL*8 :: taubeam(np, 4), taudiff(np, 4)
  REAL*8 :: taubeamb(np, 4), taudiffb(np, 4)
  REAL*8 :: fcld_col(np)
  REAL*8 :: fcld_colb(np)
  REAL*8 :: cwc_col(np, 4)
  REAL*8 :: cwc_colb(np, 4)
  REAL*8 :: reff_col(np, 4)
  REAL*8 :: reff_colb(np, 4)
  REAL*8 :: taurs, tauoz, tauwv
  REAL*8 :: tauozb, tauwvb
  REAL*8 :: tausto, ssatau, asysto
  REAL*8 :: taustob, ssataub, asystob
  REAL*8 :: tautob, ssatob, asytob
  REAL*8 :: tautobb, ssatobb, asytobb
  REAL*8 :: tautof, ssatof, asytof
  REAL*8 :: tautofb, ssatofb, asytofb
  REAL*8 :: rr(0:np+1, 2), tt(0:np+1, 2), td(0:np+1, 2)
  REAL*8 :: rrb(0:np+1, 2), ttb(0:np+1, 2), tdb(0:np+1, 2)
  REAL*8 :: rs(0:np+1, 2), ts(0:np+1, 2)
  REAL*8 :: rsb(0:np+1, 2), tsb(0:np+1, 2)
  REAL*8 :: fall(np+1), fclr(np+1), fsdir, fsdif
  REAL*8 :: fallb(np+1)
  REAL*8 :: fupa(np+1), fupc(np+1)
  REAL*8 :: cc1, cc2, cc3
  REAL*8 :: cc1b, cc2b, cc3b
  REAL*8 :: rrt, ttt, tdt, rst, tst
  REAL*8 :: rrtb, tttb, tdtb, rstb, tstb
  INTEGER :: iv, ik
  REAL*8 :: ssacl(np)
  REAL*8 :: ssaclb(np)
  INTEGER :: im
  INTEGER :: ic, iw
  REAL*8 :: ulog, wlog, dc, dd, x0, x1, x2, y0, y1, y2, du2, dw2
  REAL*8 :: wlogb, ddb, x2b, y2b
  INTEGER :: ih
!if (overcast == true) then
!real(8) :: rra(0:np+1),rxa(0:np+1)
!real(8) :: ttaold,tdaold,rsaold
!real(8) :: ttanew,tdanew,rsanew 
!else
  REAL*8 :: rra(0:np+1, 2, 2), tta(0:np, 2, 2)
  REAL*8 :: rrab(0:np+1, 2, 2), ttab(0:np, 2, 2)
  REAL*8 :: tda(0:np, 2, 2)
  REAL*8 :: tdab(0:np, 2, 2)
  REAL*8 :: rsa(0:np, 2, 2), rxa(0:np+1, 2, 2)
  REAL*8 :: rsab(0:np, 2, 2), rxab(0:np+1, 2, 2)
!endif
  REAL*8 :: flxdn
  REAL*8 :: flxdnb
  REAL*8 :: fdndir, fdndif, fupdif
  REAL*8 :: fdndirb, fdndifb, fupdifb
  REAL*8 :: denm, yy
  REAL*8 :: denmb, yyb
  INTEGER :: is
  REAL*8 :: ch, cm, ct
  REAL*8 :: chb, cmb, ctb
  INTEGER :: foundtop
  REAL*8 :: dftop
  REAL*8 :: dftopb
!-----Variables for aerosols
  INTEGER :: ii, jj, irhp1, an
  REAL*8 :: dum
  REAL*8 :: dumb
  INTRINSIC MAX
  INTRINSIC EXP
  INTRINSIC MIN
  INTRINSIC SQRT
  INTRINSIC REAL
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC ABS
  INTRINSIC EPSILON
  REAL :: result1
  INTEGER :: branch
  REAL*8 :: temp3
  REAL*8 :: temp29
  REAL*8 :: tempb52
  REAL*8 :: temp2
  REAL*8 :: temp28
  REAL*8 :: tempb51
  REAL*8 :: temp1
  REAL*8 :: temp27
  REAL*8 :: tempb50
  REAL*8 :: temp0
  REAL*8 :: temp26
  REAL*8 :: temp25
  REAL*8 :: temp24
  REAL*8 :: temp23
  REAL*8 :: temp22
  REAL*8 :: temp21
  REAL*8 :: temp20
  REAL*8 :: tempb9
  REAL*8 :: tempb8
  REAL*8 :: tempb7
  REAL*8 :: tempb6
  REAL*8 :: tempb5
  REAL*8 :: tempb4
  REAL*8 :: tempb19
  REAL*8 :: tempb3
  REAL*8 :: tempb18
  REAL*8 :: tempb2
  REAL*8 :: tempb17
  REAL*8 :: tempb1
  REAL*8 :: tempb16
  REAL*8 :: tempb0
  REAL*8 :: tempb15
  REAL*8 :: tempb14
  REAL*8 :: tempb13
  REAL*8 :: tempb12
  REAL*8 :: tempb49
  REAL*8 :: x6
  REAL*8 :: tempb11
  REAL*8 :: tempb48
  REAL*8 :: x5
  REAL*8 :: tempb10
  REAL*8 :: tempb47
  REAL*8 :: x4
  REAL*8 :: tempb46
  REAL*8 :: x3
  REAL*8 :: tempb45
  REAL*8 :: tempb44
  REAL*8 :: tempb43
  REAL*8 :: temp19
  REAL*8 :: tempb42
  REAL*8 :: temp18
  REAL*8 :: tempb41
  REAL*8 :: temp17
  REAL*8 :: tempb40
  REAL*8 :: temp16
  REAL*8 :: temp15
  REAL*8 :: temp14
  REAL*8 :: temp13
  REAL*8 :: temp12
  REAL*8 :: temp11
  REAL*8 :: temp10
  REAL*8 :: tempb
  REAL*8 :: tempb39
  REAL*8 :: tempb38
  REAL*8 :: tempb37
  REAL*8 :: tempb36
  REAL*8 :: tempb35
  REAL*8 :: tempb34
  REAL*8 :: tempb33
  REAL*8 :: tempb32
  REAL*8 :: tempb31
  REAL*8 :: tempb30
  REAL*8 :: x4b
  REAL*8 :: temp38
  REAL*8 :: temp37
  REAL*8 :: temp36
  REAL*8 :: temp35
  REAL*8 :: temp34
  REAL*8 :: temp33
  REAL*8 :: temp32
  REAL*8 :: temp31
  REAL*8 :: temp30
  REAL*8 :: abs0
  REAL*8 :: tempb29
  REAL*8 :: tempb28
  REAL*8 :: tempb27
  REAL*8 :: tempb26
  REAL*8 :: tempb25
  REAL*8 :: temp
  REAL*8 :: tempb24
  REAL*8 :: tempb23
  REAL*8 :: tempb22
  REAL*8 :: temp9
  REAL*8 :: tempb21
  REAL*8 :: temp8
  REAL*8 :: tempb20
  REAL*8 :: temp7
  REAL*8 :: tempb56
  REAL*8 :: temp6
  REAL*8 :: tempb55
  REAL*8 :: temp5
  REAL*8 :: tempb54
  REAL*8 :: temp4
  REAL*8 :: tempb53
  i = 1
!RUN_LOOP: do i=1,m
  ntop = 0
!-----Beginning of sorad code
!-----wvtoa and o3toa are the water vapor and o3 amounts of the region 
!     above the pl(1) level.
!     snt is the secant of the solar zenith angle
  snt = 1.0/cosz_dev(i)
  IF (pl_dev(i, 1) .LT. 1.e-3) THEN
    xtoa = 1.e-3
  ELSE
    xtoa = pl_dev(i, 1)
  END IF
  scal0 = xtoa*(0.5*xtoa/300.)**.8
  o3toa = 1.02*oa_dev(i, 1)*xtoa*466.7 + 1.0e-8
  wvtoa = 1.02*wa_dev(i, 1)*scal0*(1.0+0.00135*(ta_dev(i, 1)-240.)) + &
&   1.0e-9
  swh(1) = wvtoa
  DO k=1,np
!-----compute layer thickness. indices for the surface level and
!     surface layer are np+1 and np, respectively.
    dp(k) = pl_dev(i, k+1) - pl_dev(i, k)
! dp in pascals
    dp_pa(k) = dp(k)*100.
!-----compute scaled water vapor amount following Eqs. (3.3) and (3.5) 
!     unit is g/cm**2
!
    pa = 0.5*(pl_dev(i, k)+pl_dev(i, k+1))
    scal(k) = dp(k)*(pa/300.)**.8
    wh(k) = 1.02*wa_dev(i, k)*scal(k)*(1.+0.00135*(ta_dev(i, k)-240.)) +&
&     1.e-9
    swh(k+1) = swh(k) + wh(k)
!-----compute ozone amount, unit is (cm-atm)stp
!     the number 466.7 is the unit conversion factor
!     from g/cm**2 to (cm-atm)stp
    oh(k) = 1.02*oa_dev(i, k)*dp(k)*466.7 + 1.e-8
!-----Fill the reff, cwc, and fcld for the column
    fcld_col(k) = fcld_dev(i, k)
    DO l=1,4
      reff_col(k, l) = reff_dev(i, k, l)
      cwc_col(k, l) = cwc_dev(i, k, l)
    END DO
  END DO
!-----Initialize temporary arrays to zero to avoid UMR
  rr = 0.0
  tt = 0.0
  td = 0.0
  rs = 0.0
  ts = 0.0
  rra = 0.0
  rxa = 0.0
!if( OVERCAST == .false. ) then
  tta = 0.0
  tda = 0.0
  rsa = 0.0
!endif
!-----initialize fluxes for all-sky (flx), clear-sky (flc), and
!     flux reduction (df)
!
  DO k=1,np+1
    flx_dev(i, k) = 0.
  END DO
!-----Begin inline of SOLUV
!-----compute solar uv and par fluxes
!-----initialize fdiruv, fdifuv, surface reflectances and transmittances.
!     the reflectance and transmittance of the clear and cloudy portions
!     of a layer are denoted by 1 and 2, respectively.
!     cc is the maximum cloud cover in each of the high, middle, and low
!     cloud groups.
!     1/dsm=1/cos(53) = 1.66
  rr(np+1, 1) = rsuvbm_dev
  rr(np+1, 2) = rsuvbm_dev
  rs(np+1, 1) = rsuvdf_dev
  rs(np+1, 2) = rsuvdf_dev
  td(np+1, 1) = 0.0
  td(np+1, 2) = 0.0
  tt(np+1, 1) = 0.0
  tt(np+1, 2) = 0.0
  ts(np+1, 1) = 0.0
  ts(np+1, 2) = 0.0
  rr(0, 1) = 0.0
  rr(0, 2) = 0.0
  rs(0, 1) = 0.0
  rs(0, 2) = 0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
  tt(0, 1) = 1.0
  tt(0, 2) = 1.0
  ts(0, 1) = 1.0
  ts(0, 2) = 1.0
  cc1 = 0.0
  cc2 = 0.0
  cc3 = 0.0
!-----options for scaling cloud optical thickness
!if ( OVERCAST == .true. ) then
!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.
!         call getvistau1(np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
!                        taubeam,taudiff,asycl,                                  &
!                         aig_uv, awg_uv, arg_uv,                                 &
!                         aib_uv, awb_uv, arb_uv,                                 &
!                         aib_nir, awb_nir, arb_nir,                              &
!                         aia_nir, awa_nir, ara_nir,                              &
!                         aig_nir, awg_nir, arg_nir,                              &
!                         caib, caif,                                             &
!                         CONS_GRAV                                               )
!else
!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.
  CALL GETVISTAU1(np, cosz_dev(i), dp_pa, fcld_col, reff_col, cwc_col, &
&           ict, icb, taubeam, taudiff, asycl, aig_uv, awg_uv, arg_uv, &
&           aib_uv, awb_uv, arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, &
&           awa_nir, ara_nir, aig_nir, awg_nir, arg_nir, caib, caif, &
&           cons_grav)
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
!     The cc1,2,3 are still needed in the flux calculations below
!MAT---DO NOT FUSE THIS LOOP
!MAT---Loop must run to completion so that cc[1,2,3] are correct.
  DO k=1,np
    IF (k .LT. ict) THEN
      IF (cc1 .LT. fcld_dev(i, k)) THEN
        cc1 = fcld_dev(i, k)
        CALL PUSHCONTROL3B(4)
      ELSE
        CALL PUSHCONTROL3B(5)
        cc1 = cc1
      END IF
    ELSE IF (k .LT. icb) THEN
      IF (cc2 .LT. fcld_dev(i, k)) THEN
        cc2 = fcld_dev(i, k)
        CALL PUSHCONTROL3B(2)
      ELSE
        CALL PUSHCONTROL3B(3)
        cc2 = cc2
      END IF
    ELSE IF (cc3 .LT. fcld_dev(i, k)) THEN
      cc3 = fcld_dev(i, k)
      CALL PUSHCONTROL3B(0)
    ELSE
      CALL PUSHCONTROL3B(1)
      cc3 = cc3
    END IF
  END DO
!MAT---DO NOT FUSE THIS LOOP
!endif !overcast
  DO k=1,np
    tauclb(k) = taubeam(k, 1) + taubeam(k, 2) + taubeam(k, 3) + taubeam(&
&     k, 4)
    tauclf(k) = taudiff(k, 1) + taudiff(k, 2) + taudiff(k, 3) + taudiff(&
&     k, 4)
  END DO
!-----integration over spectral bands
!-----Compute optical thickness, single-scattering albedo and asymmetry
!     factor for a mixture of "na" aerosol types. [Eqs. (4.16)-(4.18)]
  DO ib=1,nband_uv
!-----compute direct beam transmittances of the layer above pl(1)
    CALL PUSHREAL8(td(0, 1))
    td(0, 1) = EXP(-((wvtoa*wk_uv(ib)+o3toa*zk_uv(ib))/cosz_dev(i)))
    CALL PUSHREAL8(td(0, 2))
    td(0, 2) = td(0, 1)
    DO k=1,np
!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor (Eqs. 6.2-6.4)
      taurs = ry_uv(ib)*dp(k)
      tauoz = zk_uv(ib)*oh(k)
      tauwv = wk_uv(ib)*wh(k)
      tausto = taurs + tauoz + tauwv + taua_dev(i, k, ib) + 1.0e-7
      ssatau = ssaa_dev(i, k, ib) + taurs
      asysto = asya_dev(i, k, ib)
      CALL PUSHREAL8(tautob)
      tautob = tausto
      asytob = asysto/ssatau
      CALL PUSHREAL8(ssatob)
      ssatob = ssatau/tautob + 1.0e-8
      IF (ssatob .GT. 0.999999) THEN
        ssatob = 0.999999
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        ssatob = ssatob
      END IF
!-----for direct incident radiation
      CALL DELEDD(tautob, ssatob, asytob, cosz_dev(i), rrt, ttt, tdt)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
      CALL DELEDD(tautob, ssatob, asytob, dsm, rst, tst, dum)
      CALL PUSHREAL8(rr(k, 1))
      rr(k, 1) = rrt
      CALL PUSHREAL8(tt(k, 1))
      tt(k, 1) = ttt
      CALL PUSHREAL8(td(k, 1))
      td(k, 1) = tdt
      CALL PUSHREAL8(rs(k, 1))
      rs(k, 1) = rst
      CALL PUSHREAL8(ts(k, 1))
      ts(k, 1) = tst
!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer
!-----for direct incident radiation
!     The effective layer optical properties. Eqs. (6.2)-(6.4)
      tautob = tausto + tauclb(k)
      CALL PUSHREAL8(ssatob)
      ssatob = (ssatau+tauclb(k))/tautob + 1.0e-8
      IF (ssatob .GT. 0.999999) THEN
        ssatob = 0.999999
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        ssatob = ssatob
      END IF
      asytob = (asysto+asycl(k)*tauclb(k))/(ssatob*tautob)
!-----for diffuse incident radiation
      CALL PUSHREAL8(tautof)
      tautof = tausto + tauclf(k)
      CALL PUSHREAL8(ssatof)
      ssatof = (ssatau+tauclf(k))/tautof + 1.0e-8
      IF (ssatof .GT. 0.999999) THEN
        ssatof = 0.999999
        CALL PUSHCONTROL1B(0)
      ELSE
        ssatof = ssatof
        CALL PUSHCONTROL1B(1)
      END IF
      asytof = (asysto+asycl(k)*tauclf(k))/(ssatof*tautof)
!-----for direct incident radiation
!     note that the cloud optical thickness is scaled differently 
!     for direct and diffuse insolation, Eqs. (7.3) and (7.4).
      CALL DELEDD(tautob, ssatob, asytob, cosz_dev(i), rrt, ttt, tdt)
!-----diffuse incident radiation is approximated by beam radiation 
!     with an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
      CALL DELEDD(tautof, ssatof, asytof, dsm, rst, tst, dum)
      CALL PUSHREAL8(rr(k, 2))
      rr(k, 2) = rrt
      CALL PUSHREAL8(tt(k, 2))
      tt(k, 2) = ttt
      CALL PUSHREAL8(td(k, 2))
      td(k, 2) = tdt
      CALL PUSHREAL8(rs(k, 2))
      rs(k, 2) = rst
      CALL PUSHREAL8(ts(k, 2))
      ts(k, 2) = tst
    END DO
!-----flux calculations
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)
    DO k=1,np+1
      fall(k) = 0.0
    END DO
!if ( OVERCAST == .true. ) then
!-----Inline CLDFLXY
!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.
!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)
!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)
!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux
!-----ih=1 for clear sky; ih=2 for cloudy sky.
!-----First set is ih = 1
!            rra(np+1)=rr(np+1,1)
!            rxa(np+1)=rs(np+1,1)
!
!            do k=np,0,-1
!               denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
!               rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
!               rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
!            end do
!
!            do k=1,np+1
!               if (k <= np) then
!                  if (k == 1) then
!                     tdaold = td(0,1)
!                     ttaold = tt(0,1)
!                     rsaold = rs(0,1)
!
!                     tdanew = 0.0
!                     ttanew = 0.0
!                     rsanew = 0.0
!                  end if
!                  denm=ts(k,1)/(1.-rsaold*rs(k,1))
!                  tdanew=tdaold*td(k,1)
!                  ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
!                  rsanew=rs(k,1)+ts(k,1)*rsaold*denm
!               end if
!
!               denm=1./(1.-rsaold*rxa(k))
!               fdndir=tdaold
!               xx4=tdaold*rra(k)
!               yy=ttaold-tdaold
!               fdndif=(xx4*rsaold+yy)*denm
!               fupdif=(xx4+yy*rxa(k))*denm
!               flxdn=fdndir+fdndif-fupdif
!               fupc(k)=fupdif 
!               fclr(k)=flxdn
!
!               tdaold = tdanew
!               ttaold = ttanew
!               rsaold = rsanew
!
!               tdanew = 0.0
!               ttanew = 0.0
!               rsanew = 0.0
!            end do
!
!!-----Second set is ih = 2
!
!            rra(np+1)=rr(np+1,2)
!            rxa(np+1)=rs(np+1,2)
!
!            do k=np,0,-1
!               denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
!               rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
!               rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
!            end do
!
!            do k=1,np+1
!               if (k <= np) then
!                  if (k == 1) then
!                     tdaold = td(0,2)
!                     ttaold = tt(0,2)
!                     rsaold = rs(0,2)
!                     tdanew = 0.0
!                     ttanew = 0.0
!                     rsanew = 0.0
!                  end if
!                  denm=ts(k,2)/(1.-rsaold*rs(k,2))
!                  tdanew=tdaold*td(k,2)
!                  ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
!                  rsanew=rs(k,2)+ts(k,2)*rsaold*denm
!               end if
!
!               denm=1./(1.-rsaold*rxa(k))
!               fdndir=tdaold
!               xx4=tdaold*rra(k)
!               yy=ttaold-tdaold
!               fdndif=(xx4*rsaold+yy)*denm
!               fupdif=(xx4+yy*rxa(k))*denm
!               flxdn=fdndir+fdndif-fupdif
!
!               fupa(k)=fupdif
!               fall(k)=flxdn
!
!               tdaold = tdanew
!               ttaold = ttanew
!               rsaold = rsanew
!
!               tdanew = 0.0
!               ttanew = 0.0
!               rsanew = 0.0
!            end do
!
!            fsdir=fdndir
!            fsdif=fdndif
!
!!-----End CLDFLXY inline
!
!else
!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)
!-----Inline CLDFLX
!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)
!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.
!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition
    DO ih=1,2
      CALL PUSHREAL8(tda(0, ih, 1))
      tda(0, ih, 1) = td(0, ih)
      CALL PUSHREAL8(tta(0, ih, 1))
      tta(0, ih, 1) = tt(0, ih)
      CALL PUSHREAL8(rsa(0, ih, 1))
      rsa(0, ih, 1) = rs(0, ih)
      CALL PUSHREAL8(tda(0, ih, 2))
      tda(0, ih, 2) = td(0, ih)
      CALL PUSHREAL8(tta(0, ih, 2))
      tta(0, ih, 2) = tt(0, ih)
      CALL PUSHREAL8(rsa(0, ih, 2))
      rsa(0, ih, 2) = rs(0, ih)
      DO k=1,ict-1
        CALL PUSHREAL8(denm)
        denm = ts(k, ih)/(1.-rsa(k-1, ih, 1)*rs(k, ih))
        CALL PUSHREAL8(tda(k, ih, 1))
        tda(k, ih, 1) = tda(k-1, ih, 1)*td(k, ih)
        CALL PUSHREAL8(tta(k, ih, 1))
        tta(k, ih, 1) = tda(k-1, ih, 1)*tt(k, ih) + (tda(k-1, ih, 1)*rsa&
&         (k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1, ih, 1))*denm
        CALL PUSHREAL8(rsa(k, ih, 1))
        rsa(k, ih, 1) = rs(k, ih) + ts(k, ih)*rsa(k-1, ih, 1)*denm
        CALL PUSHREAL8(tda(k, ih, 2))
        tda(k, ih, 2) = tda(k, ih, 1)
        CALL PUSHREAL8(tta(k, ih, 2))
        tta(k, ih, 2) = tta(k, ih, 1)
        CALL PUSHREAL8(rsa(k, ih, 2))
        rsa(k, ih, 2) = rsa(k, ih, 1)
      END DO
! k loop
!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition
      DO k=ict,icb-1
        DO im=1,2
          CALL PUSHREAL8(denm)
          denm = ts(k, im)/(1.-rsa(k-1, ih, im)*rs(k, im))
          CALL PUSHREAL8(tda(k, ih, im))
          tda(k, ih, im) = tda(k-1, ih, im)*td(k, im)
          CALL PUSHREAL8(tta(k, ih, im))
          tta(k, ih, im) = tda(k-1, ih, im)*tt(k, im) + (tda(k-1, ih, im&
&           )*rsa(k-1, ih, im)*rr(k, im)+tta(k-1, ih, im)-tda(k-1, ih, &
&           im))*denm
          CALL PUSHREAL8(rsa(k, ih, im))
          rsa(k, ih, im) = rs(k, im) + ts(k, im)*rsa(k-1, ih, im)*denm
        END DO
      END DO
    END DO
! im loop
! k loop
! ih loop
!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)
!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.
!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition
    DO is=1,2
      CALL PUSHREAL8(rra(np+1, 1, is))
      rra(np+1, 1, is) = rr(np+1, is)
      CALL PUSHREAL8(rxa(np+1, 1, is))
      rxa(np+1, 1, is) = rs(np+1, is)
      CALL PUSHREAL8(rra(np+1, 2, is))
      rra(np+1, 2, is) = rr(np+1, is)
      CALL PUSHREAL8(rxa(np+1, 2, is))
      rxa(np+1, 2, is) = rs(np+1, is)
      DO k=np,icb,-1
        CALL PUSHREAL8(denm)
        denm = ts(k, is)/(1.-rs(k, is)*rxa(k+1, 1, is))
        CALL PUSHREAL8(rra(k, 1, is))
        rra(k, 1, is) = rr(k, is) + (td(k, is)*rra(k+1, 1, is)+(tt(k, is&
&         )-td(k, is))*rxa(k+1, 1, is))*denm
        CALL PUSHREAL8(rxa(k, 1, is))
        rxa(k, 1, is) = rs(k, is) + ts(k, is)*rxa(k+1, 1, is)*denm
        CALL PUSHREAL8(rra(k, 2, is))
        rra(k, 2, is) = rra(k, 1, is)
        CALL PUSHREAL8(rxa(k, 2, is))
        rxa(k, 2, is) = rxa(k, 1, is)
      END DO
! k loop
!-----for middle clouds
      DO k=icb-1,ict,-1
        DO im=1,2
          CALL PUSHREAL8(denm)
          denm = ts(k, im)/(1.-rs(k, im)*rxa(k+1, im, is))
          CALL PUSHREAL8(rra(k, im, is))
          rra(k, im, is) = rr(k, im) + (td(k, im)*rra(k+1, im, is)+(tt(k&
&           , im)-td(k, im))*rxa(k+1, im, is))*denm
          CALL PUSHREAL8(rxa(k, im, is))
          rxa(k, im, is) = rs(k, im) + ts(k, im)*rxa(k+1, im, is)*denm
        END DO
      END DO
    END DO
! im loop
! k loop
! is loop
!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.
    DO ih=1,2
!-----clear portion 
      IF (ih .EQ. 1) THEN
        CALL PUSHREAL8(ch)
        ch = 1.0 - cc1
!-----cloudy portion
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHREAL8(ch)
        ch = cc1
        CALL PUSHCONTROL1B(0)
      END IF
      DO im=1,2
!-----clear portion 
        IF (im .EQ. 1) THEN
          CALL PUSHREAL8(cm)
          cm = ch*(1.0-cc2)
!-----cloudy portion
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHREAL8(cm)
          cm = ch*cc2
          CALL PUSHCONTROL1B(0)
        END IF
        DO is=1,2
!-----clear portion 
          IF (is .EQ. 1) THEN
            CALL PUSHREAL8(ct)
            ct = cm*(1.0-cc3)
!-----cloudy portion
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHREAL8(ct)
            ct = cm*cc3
            CALL PUSHCONTROL1B(0)
          END IF
!-----add one layer at a time, going down.
          DO k=icb,np
            CALL PUSHREAL8(denm)
            denm = ts(k, is)/(1.-rsa(k-1, ih, im)*rs(k, is))
            CALL PUSHREAL8(tda(k, ih, im))
            tda(k, ih, im) = tda(k-1, ih, im)*td(k, is)
            CALL PUSHREAL8(tta(k, ih, im))
            tta(k, ih, im) = tda(k-1, ih, im)*tt(k, is) + (tda(k-1, ih, &
&             im)*rr(k, is)*rsa(k-1, ih, im)+tta(k-1, ih, im)-tda(k-1, &
&             ih, im))*denm
            CALL PUSHREAL8(rsa(k, ih, im))
            rsa(k, ih, im) = rs(k, is) + ts(k, is)*rsa(k-1, ih, im)*denm
          END DO
! k loop
!-----add one layer at a time, going up.
          DO k=ict-1,0,-1
            CALL PUSHREAL8(denm)
            denm = ts(k, ih)/(1.-rs(k, ih)*rxa(k+1, im, is))
            CALL PUSHREAL8(rra(k, im, is))
            rra(k, im, is) = rr(k, ih) + (td(k, ih)*rra(k+1, im, is)+(tt&
&             (k, ih)-td(k, ih))*rxa(k+1, im, is))*denm
            CALL PUSHREAL8(rxa(k, im, is))
            rxa(k, im, is) = rs(k, ih) + ts(k, ih)*rxa(k+1, im, is)*denm
          END DO
! k loop
!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux
          DO k=1,np+1
            CALL PUSHREAL8(denm)
            denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
            fdndir = tda(k-1, ih, im)
            xx4 = tda(k-1, ih, im)*rra(k, im, is)
            yy = tta(k-1, ih, im) - tda(k-1, ih, im)
            fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
            fupdif = (xx4+yy*rxa(k, im, is))*denm
            flxdn = fdndir + fdndif - fupdif
!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)
            fall(k) = fall(k) + flxdn*ct
          END DO
        END DO
      END DO
    END DO
! is loop
! im loop
! ih loop
!-----End CLDFLX inline
!endif !overcast
!-----flux integration, Eq. (6.1)
    DO k=1,np+1
      flx_dev(i, k) = flx_dev(i, k) + fall(k)*hk_uv(ib)
    END DO
  END DO
!-----Inline SOLIR
!-----compute and update solar ir fluxes
  CALL PUSHREAL8(rr(np+1, 1))
  rr(np+1, 1) = rsirbm_dev
  CALL PUSHREAL8(rr(np+1, 2))
  rr(np+1, 2) = rsirbm_dev
  CALL PUSHREAL8(rs(np+1, 1))
  rs(np+1, 1) = rsirdf_dev
  CALL PUSHREAL8(rs(np+1, 2))
  rs(np+1, 2) = rsirdf_dev
  CALL PUSHREAL8(td(np+1, 1))
  td(np+1, 1) = 0.0
  CALL PUSHREAL8(td(np+1, 2))
  td(np+1, 2) = 0.0
  CALL PUSHREAL8(tt(np+1, 1))
  tt(np+1, 1) = 0.0
  CALL PUSHREAL8(tt(np+1, 2))
  tt(np+1, 2) = 0.0
  CALL PUSHREAL8(ts(np+1, 1))
  ts(np+1, 1) = 0.0
  CALL PUSHREAL8(ts(np+1, 2))
  ts(np+1, 2) = 0.0
  CALL PUSHREAL8(rr(0, 1))
  rr(0, 1) = 0.0
  CALL PUSHREAL8(rr(0, 2))
  rr(0, 2) = 0.0
  CALL PUSHREAL8(rs(0, 1))
  rs(0, 1) = 0.0
  CALL PUSHREAL8(rs(0, 2))
  rs(0, 2) = 0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
  CALL PUSHREAL8(tt(0, 1))
  tt(0, 1) = 1.0
  CALL PUSHREAL8(tt(0, 2))
  tt(0, 2) = 1.0
  CALL PUSHREAL8(ts(0, 1))
  ts(0, 1) = 1.0
  CALL PUSHREAL8(ts(0, 2))
  ts(0, 2) = 1.0
  cc1 = 0.0
  CALL PUSHREAL8(cc2)
  cc2 = 0.0
  CALL PUSHREAL8(cc3)
  cc3 = 0.0
!-----integration over spectral bands
!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10)
!     The indices 1, 2, 3 are for ice, water, rain particles,
!     respectively.
  DO ib=1,nband_ir
    CALL PUSHINTEGER4(iv)
    iv = ib + 5
!-----options for scaling cloud optical thickness
!if ( OVERCAST == .true. ) then
!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.
!            call getnirtau1(ib,np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
!                           taubeam,taudiff,asycl,ssacl,                               &
!                            aig_uv, awg_uv, arg_uv,                                    &
!                            aib_uv, awb_uv, arb_uv,                                    &
!                            aib_nir, awb_nir, arb_nir,                                 &
!                            aia_nir, awa_nir, ara_nir,                                 &
!                            aig_nir, awg_nir, arg_nir,                                 &
!                            caib, caif,                                                &
!                            CONS_GRAV                                                  )
!else
!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.
    CALL PUSHREAL8ARRAY(ssacl, np)
    CALL PUSHREAL8ARRAY(asycl, np)
    CALL GETNIRTAU1(ib, np, cosz_dev(i), dp_pa, fcld_col, reff_col, &
&             cwc_col, ict, icb, taubeam, taudiff, asycl, ssacl, aig_uv&
&             , awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, aib_nir, awb_nir&
&             , arb_nir, aia_nir, awa_nir, ara_nir, aig_nir, awg_nir, &
&             arg_nir, caib, caif, cons_grav)
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
!MAT--DO NOT FUSE THIS LOOP
!MAT  Loop must run to completion so that cc[1,2,3] are correct.
    DO k=1,np
      IF (k .LT. ict) THEN
        IF (cc1 .LT. fcld_dev(i, k)) THEN
          cc1 = fcld_dev(i, k)
          CALL PUSHCONTROL3B(4)
        ELSE
          CALL PUSHCONTROL3B(5)
          cc1 = cc1
        END IF
      ELSE IF (k .LT. icb) THEN
        IF (cc2 .LT. fcld_dev(i, k)) THEN
          CALL PUSHREAL8(cc2)
          cc2 = fcld_dev(i, k)
          CALL PUSHCONTROL3B(2)
        ELSE
          CALL PUSHREAL8(cc2)
          cc2 = cc2
          CALL PUSHCONTROL3B(3)
        END IF
      ELSE IF (cc3 .LT. fcld_dev(i, k)) THEN
        CALL PUSHREAL8(cc3)
        cc3 = fcld_dev(i, k)
        CALL PUSHCONTROL3B(0)
      ELSE
        CALL PUSHREAL8(cc3)
        cc3 = cc3
        CALL PUSHCONTROL3B(1)
      END IF
    END DO
!MAT--DO NOT FUSE THIS LOOP
!endif !overcast
    DO k=1,np
      CALL PUSHREAL8(tauclb(k))
      tauclb(k) = taubeam(k, 1) + taubeam(k, 2) + taubeam(k, 3) + &
&       taubeam(k, 4)
      CALL PUSHREAL8(tauclf(k))
      tauclf(k) = taudiff(k, 1) + taudiff(k, 2) + taudiff(k, 3) + &
&       taudiff(k, 4)
    END DO
!-----integration over the k-distribution function
    DO ik=1,nk_ir
!-----compute direct beam transmittances of the layer above pl(1)
      CALL PUSHREAL8(td(0, 1))
      td(0, 1) = EXP(-(wvtoa*xk_ir(ik)/cosz_dev(i)))
      CALL PUSHREAL8(td(0, 2))
      td(0, 2) = td(0, 1)
      DO k=1,np
        taurs = ry_ir(ib)*dp(k)
        tauwv = xk_ir(ik)*wh(k)
!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor. Eqs.(6.2)-(6.4)
        tausto = taurs + tauwv + taua_dev(i, k, iv) + 1.0e-7
        ssatau = ssaa_dev(i, k, iv) + taurs + 1.0e-8
        asysto = asya_dev(i, k, iv)
        CALL PUSHREAL8(tautob)
        tautob = tausto
        asytob = asysto/ssatau
        CALL PUSHREAL8(ssatob)
        ssatob = ssatau/tautob + 1.0e-8
        IF (ssatob .GT. 0.999999) THEN
          ssatob = 0.999999
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          ssatob = ssatob
        END IF
!-----Compute reflectance and transmittance of the clear portion 
!     of a layer
!-----for direct incident radiation
        CALL DELEDD(tautob, ssatob, asytob, cosz_dev(i), rrt, ttt, tdt)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
        CALL DELEDD(tautob, ssatob, asytob, dsm, rst, tst, dum)
        CALL PUSHREAL8(rr(k, 1))
        rr(k, 1) = rrt
        CALL PUSHREAL8(tt(k, 1))
        tt(k, 1) = ttt
        CALL PUSHREAL8(td(k, 1))
        td(k, 1) = tdt
        CALL PUSHREAL8(rs(k, 1))
        rs(k, 1) = rst
        CALL PUSHREAL8(ts(k, 1))
        ts(k, 1) = tst
!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer
!-----for direct incident radiation. Eqs.(6.2)-(6.4)
        tautob = tausto + tauclb(k)
        CALL PUSHREAL8(ssatob)
        ssatob = (ssatau+ssacl(k)*tauclb(k))/tautob + 1.0e-8
        IF (ssatob .GT. 0.999999) THEN
          ssatob = 0.999999
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          ssatob = ssatob
        END IF
        asytob = (asysto+asycl(k)*ssacl(k)*tauclb(k))/(ssatob*tautob)
!-----for diffuse incident radiation
        CALL PUSHREAL8(tautof)
        tautof = tausto + tauclf(k)
        CALL PUSHREAL8(ssatof)
        ssatof = (ssatau+ssacl(k)*tauclf(k))/tautof + 1.0e-8
        IF (ssatof .GT. 0.999999) THEN
          ssatof = 0.999999
          CALL PUSHCONTROL1B(0)
        ELSE
          ssatof = ssatof
          CALL PUSHCONTROL1B(1)
        END IF
        asytof = (asysto+asycl(k)*ssacl(k)*tauclf(k))/(ssatof*tautof)
!-----for direct incident radiation
        CALL DELEDD(tautob, ssatob, asytob, cosz_dev(i), rrt, ttt, tdt)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs.(6.5) and (6.6)
        CALL DELEDD(tautof, ssatof, asytof, dsm, rst, tst, dum)
        CALL PUSHREAL8(rr(k, 2))
        rr(k, 2) = rrt
        CALL PUSHREAL8(tt(k, 2))
        tt(k, 2) = ttt
        CALL PUSHREAL8(td(k, 2))
        td(k, 2) = tdt
        CALL PUSHREAL8(rs(k, 2))
        rs(k, 2) = rst
        CALL PUSHREAL8(ts(k, 2))
        ts(k, 2) = tst
      END DO
!-----FLUX CALCULATIONS
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)
      DO k=1,np+1
        fall(k) = 0.0
      END DO
!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.
!if ( OVERCAST == .true. ) then
!-----Inline CLDFLXY
!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)
!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)
!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux
!-----ih=1 for clear sky; ih=2 for cloudy sky.
!-----First set is ih = 1
!               rra(np+1)=rr(np+1,1)
!               rxa(np+1)=rs(np+1,1)
!
!               do k=np,0,-1
!                  denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
!                  rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
!                  rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
!               end do
!
!               do k=1,np+1
!                  if (k <= np) then
!                     if (k == 1) then
!                        tdaold = td(0,1)
!                        ttaold = tt(0,1)
!                        rsaold = rs(0,1)
!
!                        tdanew = 0.0
!                        ttanew = 0.0
!                        rsanew = 0.0
!                     end if
!                     denm=ts(k,1)/(1.-rsaold*rs(k,1))
!                     tdanew=tdaold*td(k,1)
!                     ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
!                     rsanew=rs(k,1)+ts(k,1)*rsaold*denm
!                  end if
!
!                  denm=1./(1.-rsaold*rxa(k))
!                  fdndir=tdaold
!                  xx4=tdaold*rra(k)
!                  yy=ttaold-tdaold
!                  fdndif=(xx4*rsaold+yy)*denm
!                  fupdif=(xx4+yy*rxa(k))*denm
!                  flxdn=fdndir+fdndif-fupdif
!
!                  fupc(k)=fupdif
!                  fclr(k)=flxdn
!
!                  tdaold = tdanew
!                  ttaold = ttanew
!                  rsaold = rsanew
!
!                  tdanew = 0.0
!                  ttanew = 0.0
!                  rsanew = 0.0
!               end do
!
!!-----Second set is ih = 2
!
!               rra(np+1)=rr(np+1,2)
!               rxa(np+1)=rs(np+1,2)
!
!               do k=np,0,-1
!                  denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
!                  rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
!                  rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
!               end do
!
!               do k=1,np+1
!                  if (k <= np) then
!                     if (k == 1) then
!                        tdaold = td(0,2)
!                        ttaold = tt(0,2)
!                        rsaold = rs(0,2)
!
!                        tdanew = 0.0
!                        ttanew = 0.0
!                        rsanew = 0.0
!                     end if
!                     denm=ts(k,2)/(1.-rsaold*rs(k,2))
!                     tdanew=tdaold*td(k,2)
!                     ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
!                     rsanew=rs(k,2)+ts(k,2)*rsaold*denm
!                  end if
!
!                  denm=1./(1.-rsaold*rxa(k))
!                  fdndir=tdaold
!                  xx4=tdaold*rra(k)
!                  yy=ttaold-tdaold
!                  fdndif=(xx4*rsaold+yy)*denm
!                  fupdif=(xx4+yy*rxa(k))*denm
!                  flxdn=fdndir+fdndif-fupdif
!
!                  fupa(k)=fupdif
!                  fall(k)=flxdn
!
!                  tdaold = tdanew
!                  ttaold = ttanew
!                  rsaold = rsanew
!
!                  tdanew = 0.0
!                  ttanew = 0.0
!                  rsanew = 0.0
!               end do
!
!               fsdir=fdndir
!               fsdif=fdndif
!
!!-----End CLDFLXY inline
!
!else
!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)
!-----Inline CLDFLX
!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)
!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.
!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition
      DO ih=1,2
        CALL PUSHREAL8(tda(0, ih, 1))
        tda(0, ih, 1) = td(0, ih)
        CALL PUSHREAL8(tta(0, ih, 1))
        tta(0, ih, 1) = tt(0, ih)
        CALL PUSHREAL8(rsa(0, ih, 1))
        rsa(0, ih, 1) = rs(0, ih)
        CALL PUSHREAL8(tda(0, ih, 2))
        tda(0, ih, 2) = td(0, ih)
        CALL PUSHREAL8(tta(0, ih, 2))
        tta(0, ih, 2) = tt(0, ih)
        CALL PUSHREAL8(rsa(0, ih, 2))
        rsa(0, ih, 2) = rs(0, ih)
        DO k=1,ict-1
          CALL PUSHREAL8(denm)
          denm = ts(k, ih)/(1.-rsa(k-1, ih, 1)*rs(k, ih))
          CALL PUSHREAL8(tda(k, ih, 1))
          tda(k, ih, 1) = tda(k-1, ih, 1)*td(k, ih)
          CALL PUSHREAL8(tta(k, ih, 1))
          tta(k, ih, 1) = tda(k-1, ih, 1)*tt(k, ih) + (tda(k-1, ih, 1)*&
&           rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1, ih, 1))*&
&           denm
          CALL PUSHREAL8(rsa(k, ih, 1))
          rsa(k, ih, 1) = rs(k, ih) + ts(k, ih)*rsa(k-1, ih, 1)*denm
          CALL PUSHREAL8(tda(k, ih, 2))
          tda(k, ih, 2) = tda(k, ih, 1)
          CALL PUSHREAL8(tta(k, ih, 2))
          tta(k, ih, 2) = tta(k, ih, 1)
          CALL PUSHREAL8(rsa(k, ih, 2))
          rsa(k, ih, 2) = rsa(k, ih, 1)
        END DO
! k loop
!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition
        DO k=ict,icb-1
          DO im=1,2
            CALL PUSHREAL8(denm)
            denm = ts(k, im)/(1.-rsa(k-1, ih, im)*rs(k, im))
            CALL PUSHREAL8(tda(k, ih, im))
            tda(k, ih, im) = tda(k-1, ih, im)*td(k, im)
            CALL PUSHREAL8(tta(k, ih, im))
            tta(k, ih, im) = tda(k-1, ih, im)*tt(k, im) + (tda(k-1, ih, &
&             im)*rsa(k-1, ih, im)*rr(k, im)+tta(k-1, ih, im)-tda(k-1, &
&             ih, im))*denm
            CALL PUSHREAL8(rsa(k, ih, im))
            rsa(k, ih, im) = rs(k, im) + ts(k, im)*rsa(k-1, ih, im)*denm
          END DO
        END DO
      END DO
! im loop
! k loop
! ih loop
!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)
!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.
!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition
      DO is=1,2
        CALL PUSHREAL8(rra(np+1, 1, is))
        rra(np+1, 1, is) = rr(np+1, is)
        CALL PUSHREAL8(rxa(np+1, 1, is))
        rxa(np+1, 1, is) = rs(np+1, is)
        CALL PUSHREAL8(rra(np+1, 2, is))
        rra(np+1, 2, is) = rr(np+1, is)
        CALL PUSHREAL8(rxa(np+1, 2, is))
        rxa(np+1, 2, is) = rs(np+1, is)
        DO k=np,icb,-1
          CALL PUSHREAL8(denm)
          denm = ts(k, is)/(1.-rs(k, is)*rxa(k+1, 1, is))
          CALL PUSHREAL8(rra(k, 1, is))
          rra(k, 1, is) = rr(k, is) + (td(k, is)*rra(k+1, 1, is)+(tt(k, &
&           is)-td(k, is))*rxa(k+1, 1, is))*denm
          CALL PUSHREAL8(rxa(k, 1, is))
          rxa(k, 1, is) = rs(k, is) + ts(k, is)*rxa(k+1, 1, is)*denm
          CALL PUSHREAL8(rra(k, 2, is))
          rra(k, 2, is) = rra(k, 1, is)
          CALL PUSHREAL8(rxa(k, 2, is))
          rxa(k, 2, is) = rxa(k, 1, is)
        END DO
! k loop
!-----for middle clouds
        DO k=icb-1,ict,-1
          DO im=1,2
            CALL PUSHREAL8(denm)
            denm = ts(k, im)/(1.-rs(k, im)*rxa(k+1, im, is))
            CALL PUSHREAL8(rra(k, im, is))
            rra(k, im, is) = rr(k, im) + (td(k, im)*rra(k+1, im, is)+(tt&
&             (k, im)-td(k, im))*rxa(k+1, im, is))*denm
            CALL PUSHREAL8(rxa(k, im, is))
            rxa(k, im, is) = rs(k, im) + ts(k, im)*rxa(k+1, im, is)*denm
          END DO
        END DO
      END DO
! im loop
! k loop
! is loop
!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.
      DO ih=1,2
!-----clear portion 
        IF (ih .EQ. 1) THEN
          CALL PUSHREAL8(ch)
          ch = 1.0 - cc1
!-----cloudy portion
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHREAL8(ch)
          ch = cc1
          CALL PUSHCONTROL1B(0)
        END IF
        DO im=1,2
!-----clear portion 
          IF (im .EQ. 1) THEN
            CALL PUSHREAL8(cm)
            cm = ch*(1.0-cc2)
!-----cloudy portion
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHREAL8(cm)
            cm = ch*cc2
            CALL PUSHCONTROL1B(0)
          END IF
          DO is=1,2
!-----clear portion 
            IF (is .EQ. 1) THEN
              CALL PUSHREAL8(ct)
              ct = cm*(1.0-cc3)
!-----cloudy portion
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHREAL8(ct)
              ct = cm*cc3
              CALL PUSHCONTROL1B(0)
            END IF
!-----add one layer at a time, going down.
            DO k=icb,np
              CALL PUSHREAL8(denm)
              denm = ts(k, is)/(1.-rsa(k-1, ih, im)*rs(k, is))
              CALL PUSHREAL8(tda(k, ih, im))
              tda(k, ih, im) = tda(k-1, ih, im)*td(k, is)
              CALL PUSHREAL8(tta(k, ih, im))
              tta(k, ih, im) = tda(k-1, ih, im)*tt(k, is) + (tda(k-1, ih&
&               , im)*rr(k, is)*rsa(k-1, ih, im)+tta(k-1, ih, im)-tda(k-&
&               1, ih, im))*denm
              CALL PUSHREAL8(rsa(k, ih, im))
              rsa(k, ih, im) = rs(k, is) + ts(k, is)*rsa(k-1, ih, im)*&
&               denm
            END DO
! k loop
!-----add one layer at a time, going up.
            DO k=ict-1,0,-1
              CALL PUSHREAL8(denm)
              denm = ts(k, ih)/(1.-rs(k, ih)*rxa(k+1, im, is))
              CALL PUSHREAL8(rra(k, im, is))
              rra(k, im, is) = rr(k, ih) + (td(k, ih)*rra(k+1, im, is)+(&
&               tt(k, ih)-td(k, ih))*rxa(k+1, im, is))*denm
              CALL PUSHREAL8(rxa(k, im, is))
              rxa(k, im, is) = rs(k, ih) + ts(k, ih)*rxa(k+1, im, is)*&
&               denm
            END DO
! k loop
!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux
            DO k=1,np+1
              CALL PUSHREAL8(denm)
              denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
              fdndir = tda(k-1, ih, im)
              xx4 = tda(k-1, ih, im)*rra(k, im, is)
              yy = tta(k-1, ih, im) - tda(k-1, ih, im)
              fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
              fupdif = (xx4+yy*rxa(k, im, is))*denm
              flxdn = fdndir + fdndif - fupdif
!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)
              fall(k) = fall(k) + flxdn*ct
            END DO
          END DO
        END DO
      END DO
! is loop
! im loop
! ih loop
!-----End CLDFLX inline
!endif !overcast
!-----flux integration following Eq. (6.1)
      DO k=1,np+1
        flx_dev(i, k) = flx_dev(i, k) + fall(k)*hk_ir(ib, ik)
      END DO
    END DO
  END DO
! ik loop
!-----compute pressure-scaled o2 amount following Eq. (3.5) with f=1.
!     unit is (cm-atm)stp. 165.22 = (1000/980)*23.14%*(22400/32)
!     compute flux reduction due to oxygen following Eq. (3.18). 0.0633 is the
!     fraction of insolation contained in the oxygen bands
  df(0) = 0.0
  cnt = 165.22*snt
  so2(1) = scal0*cnt
! LLT increased parameter 145 to 155 to enhance effect
  df(1) = 0.0633*(1.-EXP(-(0.000155*SQRT(so2(1)))))
  DO k=1,np
    so2(k+1) = so2(k) + scal(k)*cnt
! LLT increased parameter 145 to 155 to enhance effect
    df(k+1) = 0.0633*(1.0-EXP(-(0.000155*SQRT(so2(k+1)))))
  END DO
!-----for solar heating due to co2 scaling follows Eq(3.5) with f=1.
!     unit is (cm-atm)stp. 789 = (1000/980)*(44/28.97)*(22400/44)
  so2(1) = 789.*co2*scal0
  DO k=1,np
    so2(k+1) = so2(k) + 789.*co2*scal(k)
  END DO
!-----The updated flux reduction for co2 absorption in Band 7 where absorption due to
!     water vapor and co2 are both moderate. df is given by the second term on the
!     right-hand-side of Eq. (3.24) divided by So. so2 and swh are the co2 and
!     water vapor amounts integrated from the top of the atmosphere
  u1 = -3.0
  du = 0.15
  w1 = -4.0
  dw = 0.15
!-----Inline RFLX
  x0 = u1 + REAL(nu)*du
  y0 = w1 + REAL(nw)*dw
  x1 = u1 - 0.5*du
  y1 = w1 - 0.5*dw
  DO k=1,np+1
    x3 = LOG10(so2(k)*snt)
    IF (x3 .GT. x0) THEN
      ulog = x0
    ELSE
      ulog = x3
    END IF
    x4 = LOG10(swh(k)*snt)
    IF (x4 .GT. y0) THEN
      wlog = y0
      CALL PUSHCONTROL1B(0)
    ELSE
      wlog = x4
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHINTEGER4(ic)
    ic = INT((ulog-x1)/du + 1.)
    CALL PUSHINTEGER4(iw)
    iw = INT((wlog-y1)/dw + 1.)
    IF (ic .LT. 2) ic = 2
    IF (iw .LT. 2) iw = 2
    IF (ic .GT. nu) ic = nu
    IF (iw .GT. nw) iw = nw
    dc = ulog - REAL(ic-2)*du - u1
    dd = wlog - REAL(iw-2)*dw - w1
    x2 = cah(ic-1, iw-1) + (cah(ic-1, iw)-cah(ic-1, iw-1))/dw*dd
    y2 = x2 + (cah(ic, iw-1)-cah(ic-1, iw-1))/du*dc
    IF (y2 .LT. 0.0) THEN
      y2 = 0.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
      y2 = y2
    END IF
! LLT increase CO2 effect to help reduce cold tropopause bias
    df(k) = df(k) + 1.5*y2
  END DO
!-----df is the updated flux reduction for co2 absorption
!     in Band 8 where the co2 absorption has a large impact
!     on the heating of middle atmosphere. From the table
!     given by Eq. (3.19)
  u1 = 0.000250
  du = 0.000050
  w1 = -2.0
  dw = 0.05
!-----Inline RFLX
  x0 = u1 + REAL(nx)*du
  y0 = w1 + REAL(ny)*dw
  x1 = u1 - 0.5*du
  y1 = w1 - 0.5*dw
  DO k=1,np+1
    IF (co2*snt .GT. x0) THEN
      ulog = x0
    ELSE
      ulog = co2*snt
    END IF
    x5 = LOG10(pl_dev(i, k))
    IF (x5 .GT. y0) THEN
      wlog = y0
    ELSE
      wlog = x5
    END IF
    CALL PUSHINTEGER4(ic)
    ic = INT((ulog-x1)/du + 1.)
    CALL PUSHINTEGER4(iw)
    iw = INT((wlog-y1)/dw + 1.)
    IF (ic .LT. 2) ic = 2
    IF (iw .LT. 2) iw = 2
    IF (ic .GT. nx) ic = nx
    IF (iw .GT. ny) iw = ny
    dc = ulog - REAL(ic-2)*du - u1
    dd = wlog - REAL(iw-2)*dw - w1
    x2 = coa(ic-1, iw-1) + (coa(ic-1, iw)-coa(ic-1, iw-1))/dw*dd
    y2 = x2 + (coa(ic, iw-1)-coa(ic-1, iw-1))/du*dc
    IF (y2 .LT. 0.0) THEN
      y2 = 0.0
    ELSE
      y2 = y2
    END IF
! LLT increase CO2 effect to help reduce cold tropopause bias
    df(k) = df(k) + 1.5*y2
  END DO
!-----adjust the o2-co2 reduction below cloud top following Eq. (6.18)
  foundtop = 0
  DO k=1,np
    IF (fcld_dev(i, k) .GT. 0.02 .AND. foundtop .EQ. 0) THEN
      foundtop = 1
      ntop = k
    END IF
  END DO
  IF (foundtop .EQ. 0) ntop = np + 1
  dftop = df(ntop)
  DO k=1,np+1
    IF (k .GT. ntop) THEN
      xx4 = flx_dev(i, k)/flx_dev(i, ntop)
      CALL PUSHREAL8(df(k))
      df(k) = dftop + xx4*(df(k)-dftop)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
!-----update the net fluxes
  DO k=1,np+1
    IF (df(k) .GT. flx_dev(i, k) - 1.0e-8) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
  END DO
  dfb = 0.0_8
  DO k=np+1,1,-1
    dfb(k) = dfb(k) - flx_devb(i, k)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      flx_devb(i, k) = flx_devb(i, k) + dfb(k)
      dfb(k) = 0.0_8
    END IF
  END DO
  dftopb = 0.0_8
  DO k=np+1,1,-1
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      xx4 = flx_dev(i, k)/flx_dev(i, ntop)
      CALL POPREAL8(df(k))
      dftopb = dftopb + (1.0_8-xx4)*dfb(k)
      xx4b = (df(k)-dftop)*dfb(k)
      dfb(k) = xx4*dfb(k)
      tempb56 = xx4b/flx_dev(i, ntop)
      flx_devb(i, k) = flx_devb(i, k) + tempb56
      flx_devb(i, ntop) = flx_devb(i, ntop) - flx_dev(i, k)*tempb56/&
&       flx_dev(i, ntop)
    END IF
  END DO
  dfb(ntop) = dfb(ntop) + dftopb
  DO k=np+1,1,-1
    CALL POPINTEGER4(iw)
    CALL POPINTEGER4(ic)
  END DO
  dw = 0.15
  swhb = 0.0_8
  DO k=np+1,1,-1
    y2b = 1.5*dfb(k)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) y2b = 0.0_8
    x2b = y2b
    ddb = (cah(ic-1, iw)-cah(ic-1, iw-1))*x2b/dw
    wlogb = ddb
    CALL POPINTEGER4(iw)
    CALL POPINTEGER4(ic)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      x4b = 0.0_8
    ELSE
      x4b = wlogb
    END IF
    swhb(k) = swhb(k) + x4b/(swh(k)*LOG(10.0))
  END DO
  tsb = 0.0_8
  ttb = 0.0_8
  rsab = 0.0_8
  reff_colb = 0.0_8
  rrb = 0.0_8
  rsb = 0.0_8
  wvtoab = 0.0_8
  fallb = 0.0_8
  cc1b = 0.0_8
  cc2b = 0.0_8
  cc3b = 0.0_8
  asyclb = 0.0_8
  ttab = 0.0_8
  rxab = 0.0_8
  tdab = 0.0_8
  ssaclb = 0.0_8
  fcld_colb = 0.0_8
  cwc_colb = 0.0_8
  rrab = 0.0_8
  tauclbb = 0.0_8
  tauclfb = 0.0_8
  whb = 0.0_8
  tdb = 0.0_8
  DO ib=nband_ir,1,-1
    DO ik=nk_ir,1,-1
      DO k=np+1,1,-1
        fallb(k) = fallb(k) + hk_ir(ib, ik)*flx_devb(i, k)
      END DO
      DO ih=2,1,-1
        chb = 0.0_8
        DO im=2,1,-1
          cmb = 0.0_8
          DO is=2,1,-1
            ctb = 0.0_8
            DO k=np+1,1,-1
              yy = tta(k-1, ih, im) - tda(k-1, ih, im)
              denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
              xx4 = tda(k-1, ih, im)*rra(k, im, is)
              fupdif = (xx4+yy*rxa(k, im, is))*denm
              fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
              fdndir = tda(k-1, ih, im)
              flxdn = fdndir + fdndif - fupdif
              flxdnb = ct*fallb(k)
              ctb = ctb + flxdn*fallb(k)
              fdndirb = flxdnb
              fdndifb = flxdnb
              fupdifb = -flxdnb
              tempb53 = denm*fupdifb
              denmb = (xx4*rsa(k-1, ih, im)+yy)*fdndifb + (xx4+yy*rxa(k&
&               , im, is))*fupdifb
              tempb54 = denm*fdndifb
              xx4b = rsa(k-1, ih, im)*tempb54 + tempb53
              yyb = tempb54 + rxa(k, im, is)*tempb53
              ttab(k-1, ih, im) = ttab(k-1, ih, im) + yyb
              tdab(k-1, ih, im) = tdab(k-1, ih, im) + rra(k, im, is)*&
&               xx4b + fdndirb - yyb
              rrab(k, im, is) = rrab(k, im, is) + tda(k-1, ih, im)*xx4b
              CALL POPREAL8(denm)
              temp38 = -(rsa(k-1, ih, im)*rxa(k, im, is)) + 1.
              tempb55 = -(denmb/temp38**2)
              rxab(k, im, is) = rxab(k, im, is) + yy*tempb53 - rsa(k-1, &
&               ih, im)*tempb55
              rsab(k-1, ih, im) = rsab(k-1, ih, im) + xx4*tempb54 - rxa(&
&               k, im, is)*tempb55
            END DO
            DO k=0,ict-1,1
              CALL POPREAL8(rxa(k, im, is))
              tempb50 = rxa(k+1, im, is)*rxab(k, im, is)
              CALL POPREAL8(rra(k, im, is))
              tempb52 = denm*rrab(k, im, is)
              temp37 = rxa(k+1, im, is)
              temp36 = tt(k, ih) - td(k, ih)
              rrb(k, ih) = rrb(k, ih) + rrab(k, im, is)
              tdb(k, ih) = tdb(k, ih) + (rra(k+1, im, is)-temp37)*&
&               tempb52
              rrab(k+1, im, is) = rrab(k+1, im, is) + td(k, ih)*tempb52
              denmb = (td(k, ih)*rra(k+1, im, is)+temp36*temp37)*rrab(k&
&               , im, is) + ts(k, ih)*tempb50
              ttb(k, ih) = ttb(k, ih) + temp37*tempb52
              rrab(k, im, is) = 0.0_8
              temp35 = -(rs(k, ih)*rxa(k+1, im, is)) + 1.
              tsb(k, ih) = tsb(k, ih) + denmb/temp35 + denm*tempb50
              tempb51 = -(ts(k, ih)*denmb/temp35**2)
              rsb(k, ih) = rsb(k, ih) + rxab(k, im, is) - rxa(k+1, im, &
&               is)*tempb51
              rxab(k+1, im, is) = rxab(k+1, im, is) + ts(k, ih)*denm*&
&               rxab(k, im, is)
              rxab(k, im, is) = 0.0_8
              rxab(k+1, im, is) = rxab(k+1, im, is) + temp36*tempb52 - &
&               rs(k, ih)*tempb51
              CALL POPREAL8(denm)
            END DO
            DO k=np,icb,-1
              CALL POPREAL8(rsa(k, ih, im))
              tempb47 = rsa(k-1, ih, im)*rsab(k, ih, im)
              CALL POPREAL8(tta(k, ih, im))
              tempb49 = denm*ttab(k, ih, im)
              temp34 = rsa(k-1, ih, im)
              temp33 = tda(k-1, ih, im)*rr(k, is)
              tdab(k-1, ih, im) = tdab(k-1, ih, im) + td(k, is)*tdab(k, &
&               ih, im) + (temp34*rr(k, is)-1.0)*tempb49 + tt(k, is)*&
&               ttab(k, ih, im)
              ttb(k, is) = ttb(k, is) + tda(k-1, ih, im)*ttab(k, ih, im)
              rrb(k, is) = rrb(k, is) + temp34*tda(k-1, ih, im)*tempb49
              ttab(k-1, ih, im) = ttab(k-1, ih, im) + tempb49
              denmb = (temp33*temp34+tta(k-1, ih, im)-tda(k-1, ih, im))*&
&               ttab(k, ih, im) + ts(k, is)*tempb47
              ttab(k, ih, im) = 0.0_8
              CALL POPREAL8(tda(k, ih, im))
              tdb(k, is) = tdb(k, is) + tda(k-1, ih, im)*tdab(k, ih, im)
              tdab(k, ih, im) = 0.0_8
              temp32 = -(rsa(k-1, ih, im)*rs(k, is)) + 1.
              tsb(k, is) = tsb(k, is) + denmb/temp32 + denm*tempb47
              tempb48 = -(ts(k, is)*denmb/temp32**2)
              rsb(k, is) = rsb(k, is) + rsab(k, ih, im) - rsa(k-1, ih, &
&               im)*tempb48
              rsab(k-1, ih, im) = rsab(k-1, ih, im) + ts(k, is)*denm*&
&               rsab(k, ih, im)
              rsab(k, ih, im) = 0.0_8
              rsab(k-1, ih, im) = rsab(k-1, ih, im) + temp33*tempb49 - &
&               rs(k, is)*tempb48
              CALL POPREAL8(denm)
            END DO
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL8(ct)
              cmb = cmb + cc3*ctb
              cc3b = cc3b + cm*ctb
            ELSE
              CALL POPREAL8(ct)
              cmb = cmb + (1.0-cc3)*ctb
              cc3b = cc3b - cm*ctb
            END IF
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(cm)
            chb = chb + cc2*cmb
            cc2b = cc2b + ch*cmb
          ELSE
            CALL POPREAL8(cm)
            chb = chb + (1.0-cc2)*cmb
            cc2b = cc2b - ch*cmb
          END IF
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(ch)
          cc1b = cc1b + chb
        ELSE
          CALL POPREAL8(ch)
          cc1b = cc1b - chb
        END IF
      END DO
      DO is=2,1,-1
        DO k=ict,icb-1,1
          DO im=2,1,-1
            CALL POPREAL8(rxa(k, im, is))
            tempb44 = rxa(k+1, im, is)*rxab(k, im, is)
            CALL POPREAL8(rra(k, im, is))
            tempb46 = denm*rrab(k, im, is)
            temp31 = rxa(k+1, im, is)
            temp30 = tt(k, im) - td(k, im)
            rrb(k, im) = rrb(k, im) + rrab(k, im, is)
            tdb(k, im) = tdb(k, im) + (rra(k+1, im, is)-temp31)*tempb46
            rrab(k+1, im, is) = rrab(k+1, im, is) + td(k, im)*tempb46
            denmb = (td(k, im)*rra(k+1, im, is)+temp30*temp31)*rrab(k, &
&             im, is) + ts(k, im)*tempb44
            ttb(k, im) = ttb(k, im) + temp31*tempb46
            rrab(k, im, is) = 0.0_8
            temp29 = -(rs(k, im)*rxa(k+1, im, is)) + 1.
            tsb(k, im) = tsb(k, im) + denmb/temp29 + denm*tempb44
            tempb45 = -(ts(k, im)*denmb/temp29**2)
            rsb(k, im) = rsb(k, im) + rxab(k, im, is) - rxa(k+1, im, is)&
&             *tempb45
            rxab(k+1, im, is) = rxab(k+1, im, is) + ts(k, im)*denm*rxab(&
&             k, im, is)
            rxab(k, im, is) = 0.0_8
            rxab(k+1, im, is) = rxab(k+1, im, is) + temp30*tempb46 - rs(&
&             k, im)*tempb45
            CALL POPREAL8(denm)
          END DO
        END DO
        DO k=icb,np,1
          CALL POPREAL8(rxa(k, 2, is))
          rxab(k, 1, is) = rxab(k, 1, is) + rxab(k, 2, is)
          rxab(k, 2, is) = 0.0_8
          CALL POPREAL8(rra(k, 2, is))
          rrab(k, 1, is) = rrab(k, 1, is) + rrab(k, 2, is)
          rrab(k, 2, is) = 0.0_8
          CALL POPREAL8(rxa(k, 1, is))
          tempb41 = rxa(k+1, 1, is)*rxab(k, 1, is)
          CALL POPREAL8(rra(k, 1, is))
          tempb43 = denm*rrab(k, 1, is)
          temp28 = rxa(k+1, 1, is)
          temp27 = tt(k, is) - td(k, is)
          rrb(k, is) = rrb(k, is) + rrab(k, 1, is)
          tdb(k, is) = tdb(k, is) + (rra(k+1, 1, is)-temp28)*tempb43
          rrab(k+1, 1, is) = rrab(k+1, 1, is) + td(k, is)*tempb43
          denmb = (td(k, is)*rra(k+1, 1, is)+temp27*temp28)*rrab(k, 1, &
&           is) + ts(k, is)*tempb41
          ttb(k, is) = ttb(k, is) + temp28*tempb43
          rrab(k, 1, is) = 0.0_8
          temp26 = -(rs(k, is)*rxa(k+1, 1, is)) + 1.
          tsb(k, is) = tsb(k, is) + denmb/temp26 + denm*tempb41
          tempb42 = -(ts(k, is)*denmb/temp26**2)
          rsb(k, is) = rsb(k, is) + rxab(k, 1, is) - rxa(k+1, 1, is)*&
&           tempb42
          rxab(k+1, 1, is) = rxab(k+1, 1, is) + ts(k, is)*denm*rxab(k, 1&
&           , is)
          rxab(k, 1, is) = 0.0_8
          rxab(k+1, 1, is) = rxab(k+1, 1, is) + temp27*tempb43 - rs(k, &
&           is)*tempb42
          CALL POPREAL8(denm)
        END DO
        CALL POPREAL8(rxa(np+1, 2, is))
        rsb(np+1, is) = rsb(np+1, is) + rxab(np+1, 2, is)
        rxab(np+1, 2, is) = 0.0_8
        CALL POPREAL8(rra(np+1, 2, is))
        rrb(np+1, is) = rrb(np+1, is) + rrab(np+1, 2, is)
        rrab(np+1, 2, is) = 0.0_8
        CALL POPREAL8(rxa(np+1, 1, is))
        rsb(np+1, is) = rsb(np+1, is) + rxab(np+1, 1, is)
        rxab(np+1, 1, is) = 0.0_8
        CALL POPREAL8(rra(np+1, 1, is))
        rrb(np+1, is) = rrb(np+1, is) + rrab(np+1, 1, is)
        rrab(np+1, 1, is) = 0.0_8
      END DO
      DO ih=2,1,-1
        DO k=icb-1,ict,-1
          DO im=2,1,-1
            CALL POPREAL8(rsa(k, ih, im))
            tempb38 = rsa(k-1, ih, im)*rsab(k, ih, im)
            CALL POPREAL8(tta(k, ih, im))
            tempb40 = denm*ttab(k, ih, im)
            temp25 = rsa(k-1, ih, im)
            temp24 = tda(k-1, ih, im)*rr(k, im)
            tdab(k-1, ih, im) = tdab(k-1, ih, im) + td(k, im)*tdab(k, ih&
&             , im) + (temp25*rr(k, im)-1.0)*tempb40 + tt(k, im)*ttab(k&
&             , ih, im)
            ttb(k, im) = ttb(k, im) + tda(k-1, ih, im)*ttab(k, ih, im)
            rrb(k, im) = rrb(k, im) + temp25*tda(k-1, ih, im)*tempb40
            ttab(k-1, ih, im) = ttab(k-1, ih, im) + tempb40
            denmb = (temp24*temp25+tta(k-1, ih, im)-tda(k-1, ih, im))*&
&             ttab(k, ih, im) + ts(k, im)*tempb38
            ttab(k, ih, im) = 0.0_8
            CALL POPREAL8(tda(k, ih, im))
            tdb(k, im) = tdb(k, im) + tda(k-1, ih, im)*tdab(k, ih, im)
            tdab(k, ih, im) = 0.0_8
            temp23 = -(rsa(k-1, ih, im)*rs(k, im)) + 1.
            tsb(k, im) = tsb(k, im) + denmb/temp23 + denm*tempb38
            tempb39 = -(ts(k, im)*denmb/temp23**2)
            rsb(k, im) = rsb(k, im) + rsab(k, ih, im) - rsa(k-1, ih, im)&
&             *tempb39
            rsab(k-1, ih, im) = rsab(k-1, ih, im) + ts(k, im)*denm*rsab(&
&             k, ih, im)
            rsab(k, ih, im) = 0.0_8
            rsab(k-1, ih, im) = rsab(k-1, ih, im) + temp24*tempb40 - rs(&
&             k, im)*tempb39
            CALL POPREAL8(denm)
          END DO
        END DO
        DO k=ict-1,1,-1
          CALL POPREAL8(rsa(k, ih, 2))
          rsab(k, ih, 1) = rsab(k, ih, 1) + rsab(k, ih, 2)
          rsab(k, ih, 2) = 0.0_8
          CALL POPREAL8(tta(k, ih, 2))
          ttab(k, ih, 1) = ttab(k, ih, 1) + ttab(k, ih, 2)
          ttab(k, ih, 2) = 0.0_8
          CALL POPREAL8(tda(k, ih, 2))
          tdab(k, ih, 1) = tdab(k, ih, 1) + tdab(k, ih, 2)
          tdab(k, ih, 2) = 0.0_8
          CALL POPREAL8(rsa(k, ih, 1))
          tempb35 = rsa(k-1, ih, 1)*rsab(k, ih, 1)
          CALL POPREAL8(tta(k, ih, 1))
          tempb37 = denm*ttab(k, ih, 1)
          temp22 = rsa(k-1, ih, 1)
          temp21 = tda(k-1, ih, 1)*rr(k, ih)
          tdab(k-1, ih, 1) = tdab(k-1, ih, 1) + td(k, ih)*tdab(k, ih, 1)&
&           + (temp22*rr(k, ih)-1.0)*tempb37 + tt(k, ih)*ttab(k, ih, 1)
          ttb(k, ih) = ttb(k, ih) + tda(k-1, ih, 1)*ttab(k, ih, 1)
          rrb(k, ih) = rrb(k, ih) + temp22*tda(k-1, ih, 1)*tempb37
          ttab(k-1, ih, 1) = ttab(k-1, ih, 1) + tempb37
          denmb = (temp21*temp22+tta(k-1, ih, 1)-tda(k-1, ih, 1))*ttab(k&
&           , ih, 1) + ts(k, ih)*tempb35
          ttab(k, ih, 1) = 0.0_8
          CALL POPREAL8(tda(k, ih, 1))
          tdb(k, ih) = tdb(k, ih) + tda(k-1, ih, 1)*tdab(k, ih, 1)
          tdab(k, ih, 1) = 0.0_8
          temp20 = -(rsa(k-1, ih, 1)*rs(k, ih)) + 1.
          tsb(k, ih) = tsb(k, ih) + denmb/temp20 + denm*tempb35
          tempb36 = -(ts(k, ih)*denmb/temp20**2)
          rsb(k, ih) = rsb(k, ih) + rsab(k, ih, 1) - rsa(k-1, ih, 1)*&
&           tempb36
          rsab(k-1, ih, 1) = rsab(k-1, ih, 1) + ts(k, ih)*denm*rsab(k, &
&           ih, 1)
          rsab(k, ih, 1) = 0.0_8
          rsab(k-1, ih, 1) = rsab(k-1, ih, 1) + temp21*tempb37 - rs(k, &
&           ih)*tempb36
          CALL POPREAL8(denm)
        END DO
        CALL POPREAL8(rsa(0, ih, 2))
        rsb(0, ih) = rsb(0, ih) + rsab(0, ih, 2)
        rsab(0, ih, 2) = 0.0_8
        CALL POPREAL8(tta(0, ih, 2))
        ttb(0, ih) = ttb(0, ih) + ttab(0, ih, 2)
        ttab(0, ih, 2) = 0.0_8
        CALL POPREAL8(tda(0, ih, 2))
        tdb(0, ih) = tdb(0, ih) + tdab(0, ih, 2)
        tdab(0, ih, 2) = 0.0_8
        CALL POPREAL8(rsa(0, ih, 1))
        rsb(0, ih) = rsb(0, ih) + rsab(0, ih, 1)
        rsab(0, ih, 1) = 0.0_8
        CALL POPREAL8(tta(0, ih, 1))
        ttb(0, ih) = ttb(0, ih) + ttab(0, ih, 1)
        ttab(0, ih, 1) = 0.0_8
        CALL POPREAL8(tda(0, ih, 1))
        tdb(0, ih) = tdb(0, ih) + tdab(0, ih, 1)
        tdab(0, ih, 1) = 0.0_8
      END DO
      DO k=np+1,1,-1
        fallb(k) = 0.0_8
      END DO
      DO k=np,1,-1
        CALL POPREAL8(ts(k, 2))
        tstb = tsb(k, 2)
        tsb(k, 2) = 0.0_8
        CALL POPREAL8(rs(k, 2))
        rstb = rsb(k, 2)
        rsb(k, 2) = 0.0_8
        CALL POPREAL8(td(k, 2))
        tdtb = tdb(k, 2)
        tdb(k, 2) = 0.0_8
        CALL POPREAL8(tt(k, 2))
        tttb = ttb(k, 2)
        ttb(k, 2) = 0.0_8
        CALL POPREAL8(rr(k, 2))
        rrtb = rrb(k, 2)
        rrb(k, 2) = 0.0_8
        tauwv = xk_ir(ik)*wh(k)
        taurs = ry_ir(ib)*dp(k)
        tausto = taurs + tauwv + taua_dev(i, k, iv) + 1.0e-7
        tautof = tausto + tauclf(k)
        asysto = asya_dev(i, k, iv)
        asytof = (asysto+asycl(k)*ssacl(k)*tauclf(k))/(ssatof*tautof)
        tautofb = 0.0_8
        ssatofb = 0.0_8
        asytofb = 0.0_8
        dumb = 0.0_8
        CALL DELEDD_B(tautof, tautofb, ssatof, ssatofb, asytof, asytofb&
&               , dsm, rst, rstb, tst, tstb, dum, dumb)
        tautob = tausto + tauclb(k)
        asytob = (asysto+asycl(k)*ssacl(k)*tauclb(k))/(ssatob*tautob)
        tautobb = 0.0_8
        ssatobb = 0.0_8
        asytobb = 0.0_8
        CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, asytobb&
&               , cosz_dev(i), rrt, rrtb, ttt, tttb, tdt, tdtb)
        tempb33 = asytofb/(ssatof*tautof)
        temp19 = asycl(k)*ssacl(k)
        tempb34 = -((asysto+temp19*tauclf(k))*tempb33/(ssatof*tautof))
        asystob = tempb33
        asyclb(k) = asyclb(k) + tauclf(k)*ssacl(k)*tempb33
        ssaclb(k) = ssaclb(k) + tauclf(k)*asycl(k)*tempb33
        tauclfb(k) = tauclfb(k) + temp19*tempb33
        ssatofb = ssatofb + tautof*tempb34
        tautofb = tautofb + ssatof*tempb34
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          ssatau = ssaa_dev(i, k, iv) + taurs + 1.0e-8
          ssatofb = 0.0_8
        ELSE
          ssatau = ssaa_dev(i, k, iv) + taurs + 1.0e-8
        END IF
        tempb31 = asytobb/(ssatob*tautob)
        CALL POPREAL8(ssatof)
        tempb30 = ssatofb/tautof
        ssataub = tempb30
        ssaclb(k) = ssaclb(k) + tauclb(k)*asycl(k)*tempb31 + tauclf(k)*&
&         tempb30
        tautofb = tautofb - (ssatau+ssacl(k)*tauclf(k))*tempb30/tautof
        tauclfb(k) = tauclfb(k) + tautofb + ssacl(k)*tempb30
        CALL POPREAL8(tautof)
        taustob = tautofb
        temp18 = asycl(k)*ssacl(k)
        tempb32 = -((asysto+temp18*tauclb(k))*tempb31/(ssatob*tautob))
        asystob = asystob + tempb31
        asyclb(k) = asyclb(k) + tauclb(k)*ssacl(k)*tempb31
        tauclbb(k) = tauclbb(k) + temp18*tempb31
        ssatobb = ssatobb + tautob*tempb32
        tautobb = tautobb + ssatob*tempb32
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) ssatobb = 0.0_8
        CALL POPREAL8(ssatob)
        tempb29 = ssatobb/tautob
        ssataub = ssataub + tempb29
        ssaclb(k) = ssaclb(k) + tauclb(k)*tempb29
        tautobb = tautobb - (ssatau+ssacl(k)*tauclb(k))*tempb29/tautob
        tauclbb(k) = tauclbb(k) + tautobb + ssacl(k)*tempb29
        taustob = taustob + tautobb
        CALL POPREAL8(ts(k, 1))
        tstb = tsb(k, 1)
        tsb(k, 1) = 0.0_8
        CALL POPREAL8(rs(k, 1))
        rstb = rsb(k, 1)
        rsb(k, 1) = 0.0_8
        CALL POPREAL8(td(k, 1))
        tdtb = tdb(k, 1)
        tdb(k, 1) = 0.0_8
        CALL POPREAL8(tt(k, 1))
        tttb = ttb(k, 1)
        ttb(k, 1) = 0.0_8
        CALL POPREAL8(rr(k, 1))
        rrtb = rrb(k, 1)
        rrb(k, 1) = 0.0_8
        tautob = tausto
        asytob = asysto/ssatau
        tautobb = 0.0_8
        ssatobb = 0.0_8
        asytobb = 0.0_8
        dumb = 0.0_8
        CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, asytobb&
&               , dsm, rst, rstb, tst, tstb, dum, dumb)
        CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, asytobb&
&               , cosz_dev(i), rrt, rrtb, ttt, tttb, tdt, tdtb)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) ssatobb = 0.0_8
        CALL POPREAL8(ssatob)
        ssataub = ssataub + ssatobb/tautob - asysto*asytobb/ssatau**2
        tautobb = tautobb - ssatau*ssatobb/tautob**2
        asystob = asystob + asytobb/ssatau
        CALL POPREAL8(tautob)
        taustob = taustob + tautobb
        asya_devb(i, k, iv) = asya_devb(i, k, iv) + asystob
        ssaa_devb(i, k, iv) = ssaa_devb(i, k, iv) + ssataub
        tauwvb = taustob
        taua_devb(i, k, iv) = taua_devb(i, k, iv) + taustob
        whb(k) = whb(k) + xk_ir(ik)*tauwvb
      END DO
      CALL POPREAL8(td(0, 2))
      tdb(0, 1) = tdb(0, 1) + tdb(0, 2)
      tdb(0, 2) = 0.0_8
      CALL POPREAL8(td(0, 1))
      wvtoab = wvtoab - xk_ir(ik)*EXP(-(xk_ir(ik)*(wvtoa/cosz_dev(i))))*&
&       tdb(0, 1)/cosz_dev(i)
      tdb(0, 1) = 0.0_8
    END DO
    taudiffb = 0.0_8
    taubeamb = 0.0_8
    DO k=np,1,-1
      CALL POPREAL8(tauclf(k))
      taudiffb(k, 1) = taudiffb(k, 1) + tauclfb(k)
      taudiffb(k, 2) = taudiffb(k, 2) + tauclfb(k)
      taudiffb(k, 3) = taudiffb(k, 3) + tauclfb(k)
      taudiffb(k, 4) = taudiffb(k, 4) + tauclfb(k)
      tauclfb(k) = 0.0_8
      CALL POPREAL8(tauclb(k))
      taubeamb(k, 1) = taubeamb(k, 1) + tauclbb(k)
      taubeamb(k, 2) = taubeamb(k, 2) + tauclbb(k)
      taubeamb(k, 3) = taubeamb(k, 3) + tauclbb(k)
      taubeamb(k, 4) = taubeamb(k, 4) + tauclbb(k)
      tauclbb(k) = 0.0_8
    END DO
    DO k=np,1,-1
      CALL POPCONTROL3B(branch)
      IF (branch .LT. 3) THEN
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(cc3)
          fcld_devb(i, k) = fcld_devb(i, k) + cc3b
          cc3b = 0.0_8
        ELSE IF (branch .EQ. 1) THEN
          CALL POPREAL8(cc3)
        ELSE
          CALL POPREAL8(cc2)
          fcld_devb(i, k) = fcld_devb(i, k) + cc2b
          cc2b = 0.0_8
        END IF
      ELSE IF (branch .EQ. 3) THEN
        CALL POPREAL8(cc2)
      ELSE IF (branch .EQ. 4) THEN
        fcld_devb(i, k) = fcld_devb(i, k) + cc1b
        cc1b = 0.0_8
      END IF
    END DO
    CALL POPREAL8ARRAY(asycl, np)
    CALL POPREAL8ARRAY(ssacl, np)
    CALL GETNIRTAU1_B(ib, np, cosz_dev(i), dp_pa, fcld_col, fcld_colb, &
&               reff_col, reff_colb, cwc_col, cwc_colb, ict, icb, &
&               taubeam, taubeamb, taudiff, taudiffb, asycl, asyclb, &
&               ssacl, ssaclb, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, &
&               arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, awa_nir, &
&               ara_nir, aig_nir, awg_nir, arg_nir, caib, caif, &
&               cons_grav)
    CALL POPINTEGER4(iv)
  END DO
  CALL POPREAL8(cc3)
  CALL POPREAL8(cc2)
  CALL POPREAL8(ts(0, 2))
  tsb(0, 2) = 0.0_8
  CALL POPREAL8(ts(0, 1))
  tsb(0, 1) = 0.0_8
  CALL POPREAL8(tt(0, 2))
  ttb(0, 2) = 0.0_8
  CALL POPREAL8(tt(0, 1))
  ttb(0, 1) = 0.0_8
  CALL POPREAL8(rs(0, 2))
  rsb(0, 2) = 0.0_8
  CALL POPREAL8(rs(0, 1))
  rsb(0, 1) = 0.0_8
  CALL POPREAL8(rr(0, 2))
  rrb(0, 2) = 0.0_8
  CALL POPREAL8(rr(0, 1))
  rrb(0, 1) = 0.0_8
  CALL POPREAL8(ts(np+1, 2))
  tsb(np+1, 2) = 0.0_8
  CALL POPREAL8(ts(np+1, 1))
  tsb(np+1, 1) = 0.0_8
  CALL POPREAL8(tt(np+1, 2))
  ttb(np+1, 2) = 0.0_8
  CALL POPREAL8(tt(np+1, 1))
  ttb(np+1, 1) = 0.0_8
  CALL POPREAL8(td(np+1, 2))
  tdb(np+1, 2) = 0.0_8
  CALL POPREAL8(td(np+1, 1))
  tdb(np+1, 1) = 0.0_8
  CALL POPREAL8(rs(np+1, 2))
  rsb(np+1, 2) = 0.0_8
  CALL POPREAL8(rs(np+1, 1))
  rsb(np+1, 1) = 0.0_8
  CALL POPREAL8(rr(np+1, 2))
  rrb(np+1, 2) = 0.0_8
  CALL POPREAL8(rr(np+1, 1))
  rrb(np+1, 1) = 0.0_8
  ohb = 0.0_8
  cc1b = 0.0_8
  cc2b = 0.0_8
  cc3b = 0.0_8
  o3toab = 0.0_8
  DO ib=nband_uv,1,-1
    DO k=np+1,1,-1
      fallb(k) = fallb(k) + hk_uv(ib)*flx_devb(i, k)
    END DO
    DO ih=2,1,-1
      chb = 0.0_8
      DO im=2,1,-1
        cmb = 0.0_8
        DO is=2,1,-1
          ctb = 0.0_8
          DO k=np+1,1,-1
            yy = tta(k-1, ih, im) - tda(k-1, ih, im)
            denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
            xx4 = tda(k-1, ih, im)*rra(k, im, is)
            fupdif = (xx4+yy*rxa(k, im, is))*denm
            fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
            fdndir = tda(k-1, ih, im)
            flxdn = fdndir + fdndif - fupdif
            flxdnb = ct*fallb(k)
            ctb = ctb + flxdn*fallb(k)
            fdndirb = flxdnb
            fdndifb = flxdnb
            fupdifb = -flxdnb
            tempb26 = denm*fupdifb
            denmb = (xx4*rsa(k-1, ih, im)+yy)*fdndifb + (xx4+yy*rxa(k, &
&             im, is))*fupdifb
            tempb27 = denm*fdndifb
            xx4b = rsa(k-1, ih, im)*tempb27 + tempb26
            yyb = tempb27 + rxa(k, im, is)*tempb26
            ttab(k-1, ih, im) = ttab(k-1, ih, im) + yyb
            tdab(k-1, ih, im) = tdab(k-1, ih, im) + rra(k, im, is)*xx4b &
&             + fdndirb - yyb
            rrab(k, im, is) = rrab(k, im, is) + tda(k-1, ih, im)*xx4b
            CALL POPREAL8(denm)
            temp17 = -(rsa(k-1, ih, im)*rxa(k, im, is)) + 1.
            tempb28 = -(denmb/temp17**2)
            rxab(k, im, is) = rxab(k, im, is) + yy*tempb26 - rsa(k-1, ih&
&             , im)*tempb28
            rsab(k-1, ih, im) = rsab(k-1, ih, im) + xx4*tempb27 - rxa(k&
&             , im, is)*tempb28
          END DO
          DO k=0,ict-1,1
            CALL POPREAL8(rxa(k, im, is))
            tempb23 = rxa(k+1, im, is)*rxab(k, im, is)
            CALL POPREAL8(rra(k, im, is))
            tempb25 = denm*rrab(k, im, is)
            temp16 = rxa(k+1, im, is)
            temp15 = tt(k, ih) - td(k, ih)
            rrb(k, ih) = rrb(k, ih) + rrab(k, im, is)
            tdb(k, ih) = tdb(k, ih) + (rra(k+1, im, is)-temp16)*tempb25
            rrab(k+1, im, is) = rrab(k+1, im, is) + td(k, ih)*tempb25
            denmb = (td(k, ih)*rra(k+1, im, is)+temp15*temp16)*rrab(k, &
&             im, is) + ts(k, ih)*tempb23
            ttb(k, ih) = ttb(k, ih) + temp16*tempb25
            rrab(k, im, is) = 0.0_8
            temp14 = -(rs(k, ih)*rxa(k+1, im, is)) + 1.
            tsb(k, ih) = tsb(k, ih) + denmb/temp14 + denm*tempb23
            tempb24 = -(ts(k, ih)*denmb/temp14**2)
            rsb(k, ih) = rsb(k, ih) + rxab(k, im, is) - rxa(k+1, im, is)&
&             *tempb24
            rxab(k+1, im, is) = rxab(k+1, im, is) + ts(k, ih)*denm*rxab(&
&             k, im, is)
            rxab(k, im, is) = 0.0_8
            rxab(k+1, im, is) = rxab(k+1, im, is) + temp15*tempb25 - rs(&
&             k, ih)*tempb24
            CALL POPREAL8(denm)
          END DO
          DO k=np,icb,-1
            CALL POPREAL8(rsa(k, ih, im))
            tempb20 = rsa(k-1, ih, im)*rsab(k, ih, im)
            CALL POPREAL8(tta(k, ih, im))
            tempb22 = denm*ttab(k, ih, im)
            temp13 = rsa(k-1, ih, im)
            temp12 = tda(k-1, ih, im)*rr(k, is)
            tdab(k-1, ih, im) = tdab(k-1, ih, im) + td(k, is)*tdab(k, ih&
&             , im) + (temp13*rr(k, is)-1.0)*tempb22 + tt(k, is)*ttab(k&
&             , ih, im)
            ttb(k, is) = ttb(k, is) + tda(k-1, ih, im)*ttab(k, ih, im)
            rrb(k, is) = rrb(k, is) + temp13*tda(k-1, ih, im)*tempb22
            ttab(k-1, ih, im) = ttab(k-1, ih, im) + tempb22
            denmb = (temp12*temp13+tta(k-1, ih, im)-tda(k-1, ih, im))*&
&             ttab(k, ih, im) + ts(k, is)*tempb20
            ttab(k, ih, im) = 0.0_8
            CALL POPREAL8(tda(k, ih, im))
            tdb(k, is) = tdb(k, is) + tda(k-1, ih, im)*tdab(k, ih, im)
            tdab(k, ih, im) = 0.0_8
            temp11 = -(rsa(k-1, ih, im)*rs(k, is)) + 1.
            tsb(k, is) = tsb(k, is) + denmb/temp11 + denm*tempb20
            tempb21 = -(ts(k, is)*denmb/temp11**2)
            rsb(k, is) = rsb(k, is) + rsab(k, ih, im) - rsa(k-1, ih, im)&
&             *tempb21
            rsab(k-1, ih, im) = rsab(k-1, ih, im) + ts(k, is)*denm*rsab(&
&             k, ih, im)
            rsab(k, ih, im) = 0.0_8
            rsab(k-1, ih, im) = rsab(k-1, ih, im) + temp12*tempb22 - rs(&
&             k, is)*tempb21
            CALL POPREAL8(denm)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(ct)
            cmb = cmb + cc3*ctb
            cc3b = cc3b + cm*ctb
          ELSE
            CALL POPREAL8(ct)
            cmb = cmb + (1.0-cc3)*ctb
            cc3b = cc3b - cm*ctb
          END IF
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(cm)
          chb = chb + cc2*cmb
          cc2b = cc2b + ch*cmb
        ELSE
          CALL POPREAL8(cm)
          chb = chb + (1.0-cc2)*cmb
          cc2b = cc2b - ch*cmb
        END IF
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8(ch)
        cc1b = cc1b + chb
      ELSE
        CALL POPREAL8(ch)
        cc1b = cc1b - chb
      END IF
    END DO
    DO is=2,1,-1
      DO k=ict,icb-1,1
        DO im=2,1,-1
          CALL POPREAL8(rxa(k, im, is))
          tempb17 = rxa(k+1, im, is)*rxab(k, im, is)
          CALL POPREAL8(rra(k, im, is))
          tempb19 = denm*rrab(k, im, is)
          temp10 = rxa(k+1, im, is)
          temp9 = tt(k, im) - td(k, im)
          rrb(k, im) = rrb(k, im) + rrab(k, im, is)
          tdb(k, im) = tdb(k, im) + (rra(k+1, im, is)-temp10)*tempb19
          rrab(k+1, im, is) = rrab(k+1, im, is) + td(k, im)*tempb19
          denmb = (td(k, im)*rra(k+1, im, is)+temp9*temp10)*rrab(k, im, &
&           is) + ts(k, im)*tempb17
          ttb(k, im) = ttb(k, im) + temp10*tempb19
          rrab(k, im, is) = 0.0_8
          temp8 = -(rs(k, im)*rxa(k+1, im, is)) + 1.
          tsb(k, im) = tsb(k, im) + denmb/temp8 + denm*tempb17
          tempb18 = -(ts(k, im)*denmb/temp8**2)
          rsb(k, im) = rsb(k, im) + rxab(k, im, is) - rxa(k+1, im, is)*&
&           tempb18
          rxab(k+1, im, is) = rxab(k+1, im, is) + ts(k, im)*denm*rxab(k&
&           , im, is)
          rxab(k, im, is) = 0.0_8
          rxab(k+1, im, is) = rxab(k+1, im, is) + temp9*tempb19 - rs(k, &
&           im)*tempb18
          CALL POPREAL8(denm)
        END DO
      END DO
      DO k=icb,np,1
        CALL POPREAL8(rxa(k, 2, is))
        rxab(k, 1, is) = rxab(k, 1, is) + rxab(k, 2, is)
        rxab(k, 2, is) = 0.0_8
        CALL POPREAL8(rra(k, 2, is))
        rrab(k, 1, is) = rrab(k, 1, is) + rrab(k, 2, is)
        rrab(k, 2, is) = 0.0_8
        CALL POPREAL8(rxa(k, 1, is))
        tempb14 = rxa(k+1, 1, is)*rxab(k, 1, is)
        CALL POPREAL8(rra(k, 1, is))
        tempb16 = denm*rrab(k, 1, is)
        temp7 = rxa(k+1, 1, is)
        temp6 = tt(k, is) - td(k, is)
        rrb(k, is) = rrb(k, is) + rrab(k, 1, is)
        tdb(k, is) = tdb(k, is) + (rra(k+1, 1, is)-temp7)*tempb16
        rrab(k+1, 1, is) = rrab(k+1, 1, is) + td(k, is)*tempb16
        denmb = (td(k, is)*rra(k+1, 1, is)+temp6*temp7)*rrab(k, 1, is) +&
&         ts(k, is)*tempb14
        ttb(k, is) = ttb(k, is) + temp7*tempb16
        rrab(k, 1, is) = 0.0_8
        temp5 = -(rs(k, is)*rxa(k+1, 1, is)) + 1.
        tsb(k, is) = tsb(k, is) + denmb/temp5 + denm*tempb14
        tempb15 = -(ts(k, is)*denmb/temp5**2)
        rsb(k, is) = rsb(k, is) + rxab(k, 1, is) - rxa(k+1, 1, is)*&
&         tempb15
        rxab(k+1, 1, is) = rxab(k+1, 1, is) + ts(k, is)*denm*rxab(k, 1, &
&         is)
        rxab(k, 1, is) = 0.0_8
        rxab(k+1, 1, is) = rxab(k+1, 1, is) + temp6*tempb16 - rs(k, is)*&
&         tempb15
        CALL POPREAL8(denm)
      END DO
      CALL POPREAL8(rxa(np+1, 2, is))
      rsb(np+1, is) = rsb(np+1, is) + rxab(np+1, 2, is)
      rxab(np+1, 2, is) = 0.0_8
      CALL POPREAL8(rra(np+1, 2, is))
      rrb(np+1, is) = rrb(np+1, is) + rrab(np+1, 2, is)
      rrab(np+1, 2, is) = 0.0_8
      CALL POPREAL8(rxa(np+1, 1, is))
      rsb(np+1, is) = rsb(np+1, is) + rxab(np+1, 1, is)
      rxab(np+1, 1, is) = 0.0_8
      CALL POPREAL8(rra(np+1, 1, is))
      rrb(np+1, is) = rrb(np+1, is) + rrab(np+1, 1, is)
      rrab(np+1, 1, is) = 0.0_8
    END DO
    DO ih=2,1,-1
      DO k=icb-1,ict,-1
        DO im=2,1,-1
          CALL POPREAL8(rsa(k, ih, im))
          tempb11 = rsa(k-1, ih, im)*rsab(k, ih, im)
          CALL POPREAL8(tta(k, ih, im))
          tempb13 = denm*ttab(k, ih, im)
          temp4 = rsa(k-1, ih, im)
          temp3 = tda(k-1, ih, im)*rr(k, im)
          tdab(k-1, ih, im) = tdab(k-1, ih, im) + td(k, im)*tdab(k, ih, &
&           im) + (temp4*rr(k, im)-1.0)*tempb13 + tt(k, im)*ttab(k, ih, &
&           im)
          ttb(k, im) = ttb(k, im) + tda(k-1, ih, im)*ttab(k, ih, im)
          rrb(k, im) = rrb(k, im) + temp4*tda(k-1, ih, im)*tempb13
          ttab(k-1, ih, im) = ttab(k-1, ih, im) + tempb13
          denmb = (temp3*temp4+tta(k-1, ih, im)-tda(k-1, ih, im))*ttab(k&
&           , ih, im) + ts(k, im)*tempb11
          ttab(k, ih, im) = 0.0_8
          CALL POPREAL8(tda(k, ih, im))
          tdb(k, im) = tdb(k, im) + tda(k-1, ih, im)*tdab(k, ih, im)
          tdab(k, ih, im) = 0.0_8
          temp2 = -(rsa(k-1, ih, im)*rs(k, im)) + 1.
          tsb(k, im) = tsb(k, im) + denmb/temp2 + denm*tempb11
          tempb12 = -(ts(k, im)*denmb/temp2**2)
          rsb(k, im) = rsb(k, im) + rsab(k, ih, im) - rsa(k-1, ih, im)*&
&           tempb12
          rsab(k-1, ih, im) = rsab(k-1, ih, im) + ts(k, im)*denm*rsab(k&
&           , ih, im)
          rsab(k, ih, im) = 0.0_8
          rsab(k-1, ih, im) = rsab(k-1, ih, im) + temp3*tempb13 - rs(k, &
&           im)*tempb12
          CALL POPREAL8(denm)
        END DO
      END DO
      DO k=ict-1,1,-1
        CALL POPREAL8(rsa(k, ih, 2))
        rsab(k, ih, 1) = rsab(k, ih, 1) + rsab(k, ih, 2)
        rsab(k, ih, 2) = 0.0_8
        CALL POPREAL8(tta(k, ih, 2))
        ttab(k, ih, 1) = ttab(k, ih, 1) + ttab(k, ih, 2)
        ttab(k, ih, 2) = 0.0_8
        CALL POPREAL8(tda(k, ih, 2))
        tdab(k, ih, 1) = tdab(k, ih, 1) + tdab(k, ih, 2)
        tdab(k, ih, 2) = 0.0_8
        CALL POPREAL8(rsa(k, ih, 1))
        tempb8 = rsa(k-1, ih, 1)*rsab(k, ih, 1)
        CALL POPREAL8(tta(k, ih, 1))
        tempb10 = denm*ttab(k, ih, 1)
        temp1 = rsa(k-1, ih, 1)
        temp0 = tda(k-1, ih, 1)*rr(k, ih)
        tdab(k-1, ih, 1) = tdab(k-1, ih, 1) + td(k, ih)*tdab(k, ih, 1) +&
&         (temp1*rr(k, ih)-1.0)*tempb10 + tt(k, ih)*ttab(k, ih, 1)
        ttb(k, ih) = ttb(k, ih) + tda(k-1, ih, 1)*ttab(k, ih, 1)
        rrb(k, ih) = rrb(k, ih) + temp1*tda(k-1, ih, 1)*tempb10
        ttab(k-1, ih, 1) = ttab(k-1, ih, 1) + tempb10
        denmb = (temp0*temp1+tta(k-1, ih, 1)-tda(k-1, ih, 1))*ttab(k, ih&
&         , 1) + ts(k, ih)*tempb8
        ttab(k, ih, 1) = 0.0_8
        CALL POPREAL8(tda(k, ih, 1))
        tdb(k, ih) = tdb(k, ih) + tda(k-1, ih, 1)*tdab(k, ih, 1)
        tdab(k, ih, 1) = 0.0_8
        temp = -(rsa(k-1, ih, 1)*rs(k, ih)) + 1.
        tsb(k, ih) = tsb(k, ih) + denmb/temp + denm*tempb8
        tempb9 = -(ts(k, ih)*denmb/temp**2)
        rsb(k, ih) = rsb(k, ih) + rsab(k, ih, 1) - rsa(k-1, ih, 1)*&
&         tempb9
        rsab(k-1, ih, 1) = rsab(k-1, ih, 1) + ts(k, ih)*denm*rsab(k, ih&
&         , 1)
        rsab(k, ih, 1) = 0.0_8
        rsab(k-1, ih, 1) = rsab(k-1, ih, 1) + temp0*tempb10 - rs(k, ih)*&
&         tempb9
        CALL POPREAL8(denm)
      END DO
      CALL POPREAL8(rsa(0, ih, 2))
      rsb(0, ih) = rsb(0, ih) + rsab(0, ih, 2)
      rsab(0, ih, 2) = 0.0_8
      CALL POPREAL8(tta(0, ih, 2))
      ttb(0, ih) = ttb(0, ih) + ttab(0, ih, 2)
      ttab(0, ih, 2) = 0.0_8
      CALL POPREAL8(tda(0, ih, 2))
      tdb(0, ih) = tdb(0, ih) + tdab(0, ih, 2)
      tdab(0, ih, 2) = 0.0_8
      CALL POPREAL8(rsa(0, ih, 1))
      rsb(0, ih) = rsb(0, ih) + rsab(0, ih, 1)
      rsab(0, ih, 1) = 0.0_8
      CALL POPREAL8(tta(0, ih, 1))
      ttb(0, ih) = ttb(0, ih) + ttab(0, ih, 1)
      ttab(0, ih, 1) = 0.0_8
      CALL POPREAL8(tda(0, ih, 1))
      tdb(0, ih) = tdb(0, ih) + tdab(0, ih, 1)
      tdab(0, ih, 1) = 0.0_8
    END DO
    DO k=np+1,1,-1
      fallb(k) = 0.0_8
    END DO
    DO k=np,1,-1
      CALL POPREAL8(ts(k, 2))
      tstb = tsb(k, 2)
      tsb(k, 2) = 0.0_8
      CALL POPREAL8(rs(k, 2))
      rstb = rsb(k, 2)
      rsb(k, 2) = 0.0_8
      CALL POPREAL8(td(k, 2))
      tdtb = tdb(k, 2)
      tdb(k, 2) = 0.0_8
      CALL POPREAL8(tt(k, 2))
      tttb = ttb(k, 2)
      ttb(k, 2) = 0.0_8
      CALL POPREAL8(rr(k, 2))
      rrtb = rrb(k, 2)
      rrb(k, 2) = 0.0_8
      tauwv = wk_uv(ib)*wh(k)
      taurs = ry_uv(ib)*dp(k)
      tauoz = zk_uv(ib)*oh(k)
      tausto = taurs + tauoz + tauwv + taua_dev(i, k, ib) + 1.0e-7
      tautof = tausto + tauclf(k)
      asysto = asya_dev(i, k, ib)
      asytof = (asysto+asycl(k)*tauclf(k))/(ssatof*tautof)
      tautofb = 0.0_8
      ssatofb = 0.0_8
      asytofb = 0.0_8
      dumb = 0.0_8
      CALL DELEDD_B(tautof, tautofb, ssatof, ssatofb, asytof, asytofb, &
&             dsm, rst, rstb, tst, tstb, dum, dumb)
      tautob = tausto + tauclb(k)
      asytob = (asysto+asycl(k)*tauclb(k))/(ssatob*tautob)
      tautobb = 0.0_8
      ssatobb = 0.0_8
      asytobb = 0.0_8
      CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, asytobb, &
&             cosz_dev(i), rrt, rrtb, ttt, tttb, tdt, tdtb)
      tempb6 = asytofb/(ssatof*tautof)
      tempb7 = -((asysto+asycl(k)*tauclf(k))*tempb6/(ssatof*tautof))
      asystob = tempb6
      asyclb(k) = asyclb(k) + tauclf(k)*tempb6
      tauclfb(k) = tauclfb(k) + asycl(k)*tempb6
      ssatofb = ssatofb + tautof*tempb7
      tautofb = tautofb + ssatof*tempb7
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        ssatau = ssaa_dev(i, k, ib) + taurs
        ssatofb = 0.0_8
      ELSE
        ssatau = ssaa_dev(i, k, ib) + taurs
      END IF
      CALL POPREAL8(ssatof)
      tempb3 = ssatofb/tautof
      ssataub = tempb3
      tautofb = tautofb - (ssatau+tauclf(k))*tempb3/tautof
      tauclfb(k) = tauclfb(k) + tautofb + tempb3
      CALL POPREAL8(tautof)
      taustob = tautofb
      tempb4 = asytobb/(ssatob*tautob)
      tempb5 = -((asysto+asycl(k)*tauclb(k))*tempb4/(ssatob*tautob))
      asystob = asystob + tempb4
      asyclb(k) = asyclb(k) + tauclb(k)*tempb4
      tauclbb(k) = tauclbb(k) + asycl(k)*tempb4
      ssatobb = ssatobb + tautob*tempb5
      tautobb = tautobb + ssatob*tempb5
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) ssatobb = 0.0_8
      CALL POPREAL8(ssatob)
      tempb2 = ssatobb/tautob
      ssataub = ssataub + tempb2
      tautobb = tautobb - (ssatau+tauclb(k))*tempb2/tautob
      tauclbb(k) = tauclbb(k) + tautobb + tempb2
      taustob = taustob + tautobb
      CALL POPREAL8(ts(k, 1))
      tstb = tsb(k, 1)
      tsb(k, 1) = 0.0_8
      CALL POPREAL8(rs(k, 1))
      rstb = rsb(k, 1)
      rsb(k, 1) = 0.0_8
      CALL POPREAL8(td(k, 1))
      tdtb = tdb(k, 1)
      tdb(k, 1) = 0.0_8
      CALL POPREAL8(tt(k, 1))
      tttb = ttb(k, 1)
      ttb(k, 1) = 0.0_8
      CALL POPREAL8(rr(k, 1))
      rrtb = rrb(k, 1)
      rrb(k, 1) = 0.0_8
      tautob = tausto
      asytob = asysto/ssatau
      tautobb = 0.0_8
      ssatobb = 0.0_8
      asytobb = 0.0_8
      dumb = 0.0_8
      CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, asytobb, &
&             dsm, rst, rstb, tst, tstb, dum, dumb)
      CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, asytobb, &
&             cosz_dev(i), rrt, rrtb, ttt, tttb, tdt, tdtb)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) ssatobb = 0.0_8
      CALL POPREAL8(ssatob)
      ssataub = ssataub + ssatobb/tautob - asysto*asytobb/ssatau**2
      tautobb = tautobb - ssatau*ssatobb/tautob**2
      asystob = asystob + asytobb/ssatau
      CALL POPREAL8(tautob)
      taustob = taustob + tautobb
      asya_devb(i, k, ib) = asya_devb(i, k, ib) + asystob
      ssaa_devb(i, k, ib) = ssaa_devb(i, k, ib) + ssataub
      tauozb = taustob
      tauwvb = taustob
      taua_devb(i, k, ib) = taua_devb(i, k, ib) + taustob
      whb(k) = whb(k) + wk_uv(ib)*tauwvb
      ohb(k) = ohb(k) + zk_uv(ib)*tauozb
    END DO
    CALL POPREAL8(td(0, 2))
    tdb(0, 1) = tdb(0, 1) + tdb(0, 2)
    tdb(0, 2) = 0.0_8
    CALL POPREAL8(td(0, 1))
    tempb1 = -(EXP(-((wk_uv(ib)*wvtoa+zk_uv(ib)*o3toa)/cosz_dev(i)))*tdb&
&     (0, 1)/cosz_dev(i))
    wvtoab = wvtoab + wk_uv(ib)*tempb1
    o3toab = o3toab + zk_uv(ib)*tempb1
    tdb(0, 1) = 0.0_8
  END DO
  taudiffb = 0.0_8
  taubeamb = 0.0_8
  DO k=np,1,-1
    taudiffb(k, 1) = taudiffb(k, 1) + tauclfb(k)
    taudiffb(k, 2) = taudiffb(k, 2) + tauclfb(k)
    taudiffb(k, 3) = taudiffb(k, 3) + tauclfb(k)
    taudiffb(k, 4) = taudiffb(k, 4) + tauclfb(k)
    tauclfb(k) = 0.0_8
    taubeamb(k, 1) = taubeamb(k, 1) + tauclbb(k)
    taubeamb(k, 2) = taubeamb(k, 2) + tauclbb(k)
    taubeamb(k, 3) = taubeamb(k, 3) + tauclbb(k)
    taubeamb(k, 4) = taubeamb(k, 4) + tauclbb(k)
    tauclbb(k) = 0.0_8
  END DO
  DO k=np,1,-1
    CALL POPCONTROL3B(branch)
    IF (branch .LT. 3) THEN
      IF (branch .EQ. 0) THEN
        fcld_devb(i, k) = fcld_devb(i, k) + cc3b
        cc3b = 0.0_8
      ELSE IF (branch .NE. 1) THEN
        fcld_devb(i, k) = fcld_devb(i, k) + cc2b
        cc2b = 0.0_8
      END IF
    ELSE IF (branch .NE. 3) THEN
      IF (branch .EQ. 4) THEN
        fcld_devb(i, k) = fcld_devb(i, k) + cc1b
        cc1b = 0.0_8
      END IF
    END IF
  END DO
  CALL GETVISTAU1_B(np, cosz_dev(i), dp_pa, fcld_col, fcld_colb, &
&             reff_col, reff_colb, cwc_col, cwc_colb, ict, icb, taubeam&
&             , taubeamb, taudiff, taudiffb, asycl, asyclb, aig_uv, &
&             awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, aib_nir, awb_nir, &
&             arb_nir, aia_nir, awa_nir, ara_nir, aig_nir, awg_nir, &
&             arg_nir, caib, caif, cons_grav)
  DO k=np+1,1,-1
    flx_devb(i, k) = 0.0_8
  END DO
  DO k=np,1,-1
    DO l=4,1,-1
      cwc_devb(i, k, l) = cwc_devb(i, k, l) + cwc_colb(k, l)
      cwc_colb(k, l) = 0.0_8
      reff_devb(i, k, l) = reff_devb(i, k, l) + reff_colb(k, l)
      reff_colb(k, l) = 0.0_8
    END DO
    fcld_devb(i, k) = fcld_devb(i, k) + fcld_colb(k)
    fcld_colb(k) = 0.0_8
    oa_devb(i, k) = oa_devb(i, k) + dp(k)*466.7*1.02*ohb(k)
    ohb(k) = 0.0_8
    swhb(k) = swhb(k) + swhb(k+1)
    whb(k) = whb(k) + swhb(k+1)
    swhb(k+1) = 0.0_8
    tempb0 = scal(k)*1.02*whb(k)
    wa_devb(i, k) = wa_devb(i, k) + (0.00135*(ta_dev(i, k)-240.)+1.)*&
&     tempb0
    ta_devb(i, k) = ta_devb(i, k) + wa_dev(i, k)*0.00135*tempb0
    whb(k) = 0.0_8
  END DO
  wvtoab = wvtoab + swhb(1)
  tempb = scal0*1.02*wvtoab
  wa_devb(i, 1) = wa_devb(i, 1) + (0.00135*(ta_dev(i, 1)-240.)+1.0)*&
&   tempb
  ta_devb(i, 1) = ta_devb(i, 1) + wa_dev(i, 1)*0.00135*tempb
  oa_devb(i, 1) = oa_devb(i, 1) + xtoa*466.7*1.02*o3toab
END SUBROUTINE SORAD_B

!  Differentiation of deledd in reverse (adjoint) mode:
!   gradient     of useful results: g01 tt1 td1 rr1 tau1 ssc1
!   with respect to varying inputs: g01 tau1 ssc1
!*********************************************************************
SUBROUTINE DELEDD_B(tau1, tau1b, ssc1, ssc1b, g01, g01b, cza1, rr1, rr1b&
& , tt1, tt1b, td1, td1b)
  IMPLICIT NONE
! 8 byte real
  INTEGER, PARAMETER :: real_de=8
!integer,parameter :: REAL_SP = 4 ! 4 byte real
!-----input parameters
  REAL*8, INTENT(IN) :: tau1, ssc1, g01, cza1
  REAL*8 :: tau1b, ssc1b, g01b
!-----output parameters
  REAL*8 :: rr1, tt1, td1
  REAL*8 :: rr1b, tt1b, td1b
!-----temporary parameters
  REAL*8, PARAMETER :: zero=0.0_REAL_DE
  REAL*8, PARAMETER :: one=1.0_REAL_DE
  REAL*8, PARAMETER :: two=2.0_REAL_DE
  REAL*8, PARAMETER :: three=3.0_REAL_DE
  REAL*8, PARAMETER :: four=4.0_REAL_DE
  REAL*8, PARAMETER :: fourth=0.25_REAL_DE
  REAL*8, PARAMETER :: seven=7.0_REAL_DE
  REAL*8, PARAMETER :: thresh=1.e-8_REAL_DE
  REAL*8 :: tau, ssc, g0, rr, tt, td
  REAL*8 :: taub, sscb, g0b, rrb, ttb, tdb
  REAL*8 :: zth, ff, xx, taup, sscp, gp, gm1, gm2, gm3, akk, alf1, alf2
  REAL*8 :: ffb, xxb, taupb, sscpb, gpb, gm1b, gm2b, gm3b, akkb, alf1b, &
& alf2b
  REAL*8 :: all, bll, st7, st8, cll, dll, fll, ell, st1, st2, st3, st4
  REAL*8 :: allb, bllb, st7b, st8b, cllb, dllb, fllb, ellb, st1b, st2b, &
& st3b, st4b
  INTRINSIC SQRT
  INTRINSIC ABS
  INTRINSIC EXP
  INTRINSIC MAX
  INTRINSIC REAL
  INTEGER :: branch
  REAL*8 :: tempb9
  REAL*8 :: tempb8
  REAL*8 :: tempb7
  REAL*8 :: tempb6
  REAL*8 :: tempb5
  REAL*8 :: tempb4
  REAL*8 :: tempb3
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  REAL*8 :: tempb10
  REAL*8 :: tempb
  REAL*8 :: abs0
  REAL*8 :: temp
!zth = real(cza1,kind=REAL_DE)
!g0  = real(g01 ,kind=REAL_DE)
!tau = real(tau1,kind=REAL_DE)
!ssc = real(ssc1,kind=REAL_DE)
!zth = dble(cza1)
!g0  = dble(g01)
!tau = dble(tau1)
!ssc = dble(ssc1)
  zth = cza1
  g0 = g01
  tau = tau1
  ssc = ssc1
  ff = g0*g0
  xx = one - ff*ssc
  taup = tau*xx
  sscp = ssc*(one-ff)/xx
  gp = g0/(one+g0)
  xx = three*gp
  gm1 = (seven-sscp*(four+xx))*fourth
  gm2 = -((one-sscp*(four-xx))*fourth)
  akk = SQRT((gm1+gm2)*(gm1-gm2))
  xx = akk*zth
  st7 = one - xx
  st8 = one + xx
  st3 = st7*st8
  IF (st3 .GE. 0.) THEN
    abs0 = st3
  ELSE
    abs0 = -st3
  END IF
  IF (abs0 .LT. thresh) THEN
    CALL PUSHREAL8(zth)
    zth = zth + 0.0010
    IF (zth .GT. 1.0) zth = zth - 0.0020
    xx = akk*zth
    CALL PUSHREAL8(st7)
    st7 = one - xx
    CALL PUSHREAL8(st8)
    st8 = one + xx
    st3 = st7*st8
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
  END IF
  td = EXP(-(taup/zth))
  gm3 = (two-zth*three*gp)*fourth
  xx = gm1 - gm2
  alf1 = gm1 - gm3*xx
  alf2 = gm2 + gm3*xx
  xx = akk*two
  all = (gm3-alf2*zth)*xx*td
  bll = (one-gm3+alf1*zth)*xx
  xx = akk*gm3
  cll = (alf2+xx)*st7
  dll = (alf2-xx)*st8
  xx = akk*(one-gm3)
  fll = (alf1+xx)*st8
  ell = (alf1-xx)*st7
  st2 = EXP(-(akk*taup))
  st4 = st2*st2
  st1 = sscp/((akk+gm1+(akk-gm1)*st4)*st3)
  rr = (cll-dll*st4-all*st2)*st1
  tt = -(((fll-ell*st4)*td-bll*st2)*st1)
  IF (rr .LT. zero) THEN
    rr = zero
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    rr = rr
  END IF
  IF (tt .LT. zero) THEN
    tt = zero
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    tt = tt
  END IF
  tt = tt + td
!td1 = real(td,kind=REAL_SP)
!rr1 = real(rr,kind=REAL_SP)
!tt1 = real(tt,kind=REAL_SP)
  ttb = tt1b
  rrb = rr1b
  tdb = ttb + td1b
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) ttb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rrb = 0.0_8
  st1b = (cll-dll*st4-all*st2)*rrb - ((fll-ell*st4)*td-bll*st2)*ttb
  tempb4 = st1*rrb
  temp = akk + gm1 + (akk-gm1)*st4
  tempb7 = st1b/(temp*st3)
  tempb8 = -(sscp*tempb7/(temp*st3))
  tempb5 = st3*tempb8
  tempb2 = -(st1*ttb)
  tempb3 = td*tempb2
  fllb = tempb3
  ellb = -(st4*tempb3)
  st4b = (akk-gm1)*tempb5 - dll*tempb4 - ell*tempb3
  bllb = -(st2*tempb2)
  st2b = 2*st2*st4b - all*tempb4 - bll*tempb2
  cllb = tempb4
  dllb = -(st4*tempb4)
  allb = -(st2*tempb4)
  sscpb = tempb7
  st3b = temp*tempb8
  tempb9 = EXP(-(akk*taup))*st2b
  xxb = st8*fllb - st7*ellb
  akkb = (one-gm3)*xxb - taup*tempb9 + (st4+1.0_8)*tempb5
  st7b = (alf1-xx)*ellb
  st8b = (alf1+xx)*fllb
  gm3b = -(akk*xxb)
  xx = akk*gm3
  xxb = st7*cllb - st8*dllb
  st8b = st8b + (alf2-xx)*dllb
  st7b = st7b + (alf2+xx)*cllb
  akkb = akkb + gm3*xxb
  xx = akk*two
  alf1b = st8*fllb + xx*zth*bllb + st7*ellb
  gm3b = gm3b + akk*xxb - xx*bllb
  tempb10 = xx*td*allb
  alf2b = st7*cllb - zth*tempb10 + st8*dllb
  tempb6 = (gm3-zth*alf2)*allb
  tdb = tdb + xx*tempb6 + (fll-ell*st4)*tempb2
  taupb = -(EXP(-(taup/zth))*tdb/zth) - akk*tempb9
  xxb = td*tempb6 + (one+zth*alf1-gm3)*bllb
  akkb = akkb + two*xxb
  xx = gm1 - gm2
  gm3b = gm3b + xx*alf2b - xx*alf1b + tempb10
  xxb = gm3*alf2b - gm3*alf1b
  gm1b = alf1b + xxb + (1.0_8-st4)*tempb5
  gm2b = alf2b - xxb
  gpb = -(fourth*zth*three*gm3b)
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    st7b = st7b + st8*st3b
    st8b = st8b + st7*st3b
    CALL POPREAL8(st8)
    xxb = st8b - st7b
    CALL POPREAL8(st7)
    akkb = akkb + zth*xxb
    CALL POPREAL8(zth)
    st3b = 0.0_8
    st7b = 0.0_8
    st8b = 0.0_8
  END IF
  st7b = st7b + st8*st3b
  st8b = st8b + st7*st3b
  xxb = st8b - st7b
  akkb = akkb + zth*xxb
  xx = three*gp
  IF ((gm1+gm2)*(gm1-gm2) .EQ. 0.0) THEN
    tempb = 0.0
  ELSE
    tempb = akkb/(2.0*SQRT((gm1+gm2)*(gm1-gm2)))
  END IF
  gm1b = gm1b + 2*gm1*tempb
  gm2b = gm2b - 2*gm2*tempb
  sscpb = sscpb + fourth*(four-xx)*gm2b - fourth*(four+xx)*gm1b
  xxb = -(fourth*sscp*gm1b) - sscp*fourth*gm2b
  gpb = gpb + three*xxb
  xx = one - ff*ssc
  tempb0 = gpb/(one+g0)
  tempb1 = (one-ff)*sscpb/xx
  xxb = tau*taupb - ssc*tempb1/xx
  ffb = -(ssc*xxb) - ssc*sscpb/xx
  g0b = 2*g0*ffb + (1.0_8-g0/(one+g0))*tempb0
  sscb = tempb1 - ff*xxb
  taub = xx*taupb
  ssc1b = ssc1b + sscb
  tau1b = tau1b + taub
  g01b = g01b + g0b
END SUBROUTINE DELEDD_B

!  Differentiation of getvistau1 in reverse (adjoint) mode:
!   gradient     of useful results: hydromets asycl taudiff fcld
!                taubeam reff
!   with respect to varying inputs: hydromets fcld reff
SUBROUTINE GETVISTAU1_B(nlevs, cosz, dp, fcld, fcldb, reff, reffb, &
& hydromets, hydrometsb, ict, icb, taubeam, taubeamb, taudiff, taudiffb&
& , asycl, asyclb, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, &
& aib_nir, awb_nir, arb_nir, aia_nir, awa_nir, ara_nir, aig_nir, awg_nir&
& , arg_nir, caib, caif, cons_grav)
  IMPLICIT NONE
!EOC
! !INPUT PARAMETERS:
!  Number of levels
  INTEGER, INTENT(IN) :: nlevs
!  Cosine of solar zenith angle
  REAL*8, INTENT(IN) :: cosz
!  Delta pressure (Pa)
  REAL*8, INTENT(IN) :: dp(nlevs)
!  Cloud fraction (used sometimes)
  REAL*8, INTENT(IN) :: fcld(nlevs)
  REAL*8 :: fcldb(nlevs)
!  Effective radius (microns)
  REAL*8, INTENT(IN) :: reff(nlevs, 4)
  REAL*8 :: reffb(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL*8, INTENT(IN) :: hydromets(nlevs, 4)
  REAL*8 :: hydrometsb(nlevs, 4)
!  Flags for various uses 
  INTEGER, INTENT(IN) :: ict, icb
!                 ict  = 0   Indicates that in-cloud values have been given
!                            and are expected
!                 ict != 0   Indicates that overlap computation is needed, and:
!                               ict is the level of the mid-high boundary
!                               icb is the level of the low-mid  boundary
!                
! !OUTPUT PARAMETERS:
!  Optical Depth for Beam Radiation
  REAL*8 :: taubeam(nlevs, 4)
  REAL*8 :: taubeamb(nlevs, 4)
!  Optical Depth for Diffuse Radiation
  REAL*8 :: taudiff(nlevs, 4)
  REAL*8 :: taudiffb(nlevs, 4)
!  Cloud Asymmetry Factor
  REAL*8 :: asycl(nlevs)
  REAL*8 :: asyclb(nlevs)
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for visible wavelengths
!  In general will compute in-cloud - will do grid mean when called
!  for diagnostic use only. ict flag will indicate which to do.
!  Slots for reff, hydrometeors, taubeam, taudiff, and asycl are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.10.27   Molod moved to Radiation_Shared and revised arg list, units
!    2011.11.16   MAT: Generalized to a call that is per-column
!
!EOP
!------------------------------------------------------------------------------
!BOC
  INTEGER :: k, in, im, it, ia, kk
  REAL*8 :: fm, ft, fa, xai, tauc, asyclt
  REAL*8 :: ftb, fab, xaib, taucb, asycltb
  REAL*8 :: cc(3)
  REAL*8 :: ccb(3)
  REAL*8 :: taucld1, taucld2, taucld3, taucld4
  REAL*8 :: taucld1b, taucld2b, taucld3b, taucld4b
  REAL*8 :: g1, g2, g3, g4
  REAL*8 :: g1b, g2b, g3b, g4b
  REAL*8 :: reff_snow
  REAL*8 :: reff_snowb
  INTEGER, PARAMETER :: nm=11, nt=9, na=11
  REAL*8, PARAMETER :: dm=0.1, dt=0.30103, da=0.1, t1=-0.9031
  REAL*8, INTENT(IN) :: aig_uv(3), awg_uv(3), arg_uv(3)
  REAL*8, INTENT(IN) :: aib_uv, awb_uv(2), arb_uv(2)
  REAL*8, INTENT(IN) :: aib_nir, awb_nir(3, 2), arb_nir(3, 2)
  REAL*8, INTENT(IN) :: aia_nir(3, 3), awa_nir(3, 3), ara_nir(3, 3)
  REAL*8, INTENT(IN) :: aig_nir(3, 3), awg_nir(3, 3), arg_nir(3, 3)
  REAL*8, INTENT(IN) :: caib(11, 9, 11), caif(9, 11)
  REAL*8, INTENT(IN) :: cons_grav
  INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC REAL
  INTEGER :: branch
  REAL*8 :: temp3
  REAL*8 :: temp2
  REAL*8 :: temp1
  REAL*8 :: temp0
  REAL*8 :: tempb7
  REAL*8 :: tempb6
  REAL*8 :: tempb5
  REAL*8 :: tempb4
  REAL*8 :: tempb3
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  REAL*8 :: tempb
  REAL*8 :: temp
  IF (ict .NE. 0) THEN
!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
    cc = 0.0
    DO k=1,ict-1
      IF (cc(1) .LT. fcld(k)) THEN
        cc(1) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(1) = cc(1)
      END IF
    END DO
    DO k=ict,icb-1
      IF (cc(2) .LT. fcld(k)) THEN
        cc(2) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(2) = cc(2)
      END IF
    END DO
    DO k=icb,nlevs
      IF (cc(3) .LT. fcld(k)) THEN
        cc(3) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(3) = cc(3)
      END IF
    END DO
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
  END IF
!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.) THEN
      CALL PUSHREAL8(taucld1)
      taucld1 = 0.
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL8(taucld1)
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*aib_uv/reff(k, 1)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff(k, 2) .LE. 0.) THEN
      CALL PUSHREAL8(taucld2)
      taucld2 = 0.
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL8(taucld2)
      taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_uv(1)+awb_uv(&
&       2)/reff(k, 2))
      CALL PUSHCONTROL1B(0)
    END IF
    CALL PUSHREAL8(taucld3)
    taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_uv(1)
    IF (reff(k, 4) .GT. 112.0) THEN
      CALL PUSHREAL8(reff_snow)
      reff_snow = 112.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL8(reff_snow)
      reff_snow = reff(k, 4)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff_snow .LE. 0.) THEN
      CALL PUSHREAL8(taucld4)
      taucld4 = 0.
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL8(taucld4)
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*aib_uv/reff_snow
      CALL PUSHCONTROL1B(1)
    END IF
    IF (ict .NE. 0) THEN
!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
      IF (k .LT. ict) THEN
        CALL PUSHINTEGER4(kk)
        kk = 1
        CALL PUSHCONTROL2B(0)
      ELSE IF (k .GE. ict .AND. k .LT. icb) THEN
        CALL PUSHINTEGER4(kk)
        kk = 2
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHINTEGER4(kk)
        kk = 3
        CALL PUSHCONTROL2B(2)
      END IF
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
!-----normalize cloud cover following Eq. (7.8)
        CALL PUSHREAL8(fa)
        fa = fcld(k)/cc(kk)
        IF (tauc .GT. 32.) THEN
          tauc = 32.
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          tauc = tauc
        END IF
        fm = cosz/dm
        CALL PUSHREAL8(ft)
        ft = (LOG10(tauc)-t1)/dt
        fa = fa/da
        CALL PUSHINTEGER4(im)
        im = INT(fm + 1.5)
        CALL PUSHINTEGER4(it)
        it = INT(ft + 1.5)
        CALL PUSHINTEGER4(ia)
        ia = INT(fa + 1.5)
        IF (im .LT. 2) THEN
          im = 2
        ELSE
          im = im
        END IF
        IF (it .LT. 2) THEN
          it = 2
        ELSE
          it = it
        END IF
        IF (ia .LT. 2) THEN
          ia = 2
        ELSE
          ia = ia
        END IF
        IF (im .GT. nm - 1) THEN
          im = nm - 1
        ELSE
          im = im
        END IF
        IF (it .GT. nt - 1) THEN
          it = nt - 1
        ELSE
          it = it
        END IF
        IF (ia .GT. na - 1) THEN
          ia = na - 1
        ELSE
          ia = ia
        END IF
        fm = fm - REAL(im - 1)
        ft = ft - REAL(it - 1)
        fa = fa - REAL(ia - 1)
!-----scale cloud optical thickness for beam radiation following 
!     Eq. (7.3).
!     the scaling factor, xai, is a function of the solar zenith
!     angle, optical thickness, and cloud cover.
        CALL PUSHREAL8(xai)
        xai = (-(caib(im-1, it, ia)*(1.-fm))+caib(im+1, it, ia)*(1.+fm))&
&         *fm*.5 + caib(im, it, ia)*(1.-fm*fm)
        xai = xai + (-(caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(&
&         1.+ft))*ft*.5 + caib(im, it, ia)*(1.-ft*ft)
        xai = xai + (-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(&
&         1.+fa))*fa*.5 + caib(im, it, ia)*(1.-fa*fa)
        xai = xai - 2.*caib(im, it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
        CALL PUSHREAL8(xai)
        xai = (-(caif(it-1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ft*.5 +&
&         caif(it, ia)*(1.-ft*ft)
        xai = xai + (-(caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*&
&         fa*.5 + caif(it, ia)*(1.-fa*fa)
        xai = xai - caif(it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
!-----cloud asymmetry factor for a mixture of liquid and ice particles.
!     unit of reff is micrometers. Eqs. (4.8) and (6.4)
    CALL PUSHREAL8(tauc)
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
  END DO
  ccb = 0.0_8
  DO k=nlevs,1,-1
    asycltb = asyclb(k)
    asyclb(k) = 0.0_8
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      g1 = (aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff(k, 1))*reff(k, 1))*&
&       taucld1
      g2 = (awg_uv(1)+(awg_uv(2)+awg_uv(3)*reff(k, 2))*reff(k, 2))*&
&       taucld2
      taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_uv(1)
      g3 = arg_uv(1)*taucld3
      g4 = (aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff_snow)*reff_snow)*taucld4
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      tempb7 = asycltb/tauc
      g1b = tempb7
      g2b = tempb7
      g3b = tempb7
      g4b = tempb7
      taucb = -((g1+g2+g3+g4)*tempb7/tauc)
      temp3 = aig_uv(2) + aig_uv(3)*reff_snow
      reff_snowb = (taucld4*temp3+reff_snow*taucld4*aig_uv(3))*g4b
      taucld4b = (aig_uv(1)+temp3*reff_snow)*g4b
      taucld3b = arg_uv(1)*g3b
      temp2 = awg_uv(2) + awg_uv(3)*reff(k, 2)
      reffb(k, 2) = reffb(k, 2) + (taucld2*temp2+reff(k, 2)*taucld2*&
&       awg_uv(3))*g2b
      taucld2b = (awg_uv(1)+temp2*reff(k, 2))*g2b
      temp1 = aig_uv(2) + aig_uv(3)*reff(k, 1)
      reffb(k, 1) = reffb(k, 1) + (taucld1*temp1+reff(k, 1)*taucld1*&
&       aig_uv(3))*g1b
      taucld1b = (aig_uv(1)+temp1*reff(k, 1))*g1b
    ELSE
      reff_snowb = 0.0_8
      taucld1b = 0.0_8
      taucld2b = 0.0_8
      taucld3b = 0.0_8
      taucld4b = 0.0_8
      taucb = 0.0_8
    END IF
    CALL POPREAL8(tauc)
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      taucld4b = taucld4b + xai*taudiffb(k, 4)
      xaib = taucld4*taudiffb(k, 4)
      taudiffb(k, 4) = 0.0_8
      taucld3b = taucld3b + xai*taudiffb(k, 3)
      xaib = xaib + taucld3*taudiffb(k, 3)
      taudiffb(k, 3) = 0.0_8
      taucld2b = taucld2b + xai*taudiffb(k, 2)
      xaib = xaib + taucld2*taudiffb(k, 2)
      taudiffb(k, 2) = 0.0_8
      taucld1b = taucld1b + xai*taudiffb(k, 1)
      xaib = xaib + taucld1*taudiffb(k, 1)
      taudiffb(k, 1) = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0_8
      tempb5 = .5*fa*xaib
      fab = (.5*(caif(it, ia+1)*(fa+1.)-caif(it, ia-1)*(1.-fa))-caif(it&
&       , ia)*2*fa)*xaib + (caif(it, ia-1)+caif(it, ia+1))*tempb5
      CALL POPREAL8(xai)
      tempb6 = .5*ft*xaib
      ftb = (.5*(caif(it+1, ia)*(ft+1.)-caif(it-1, ia)*(1.-ft))-caif(it&
&       , ia)*2*ft)*xaib + (caif(it-1, ia)+caif(it+1, ia))*tempb6
      taucld4b = taucld4b + xai*taubeamb(k, 4)
      xaib = taucld4*taubeamb(k, 4)
      taubeamb(k, 4) = 0.0_8
      taucld3b = taucld3b + xai*taubeamb(k, 3)
      xaib = xaib + taucld3*taubeamb(k, 3)
      taubeamb(k, 3) = 0.0_8
      taucld2b = taucld2b + xai*taubeamb(k, 2)
      xaib = xaib + taucld2*taubeamb(k, 2)
      taubeamb(k, 2) = 0.0_8
      taucld1b = taucld1b + xai*taubeamb(k, 1)
      xaib = xaib + taucld1*taubeamb(k, 1)
      taubeamb(k, 1) = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0_8
      tempb3 = .5*fa*xaib
      fab = fab + (.5*(caib(im, it, ia+1)*(fa+1.)-caib(im, it, ia-1)*(1.&
&       -fa))-caib(im, it, ia)*2*fa)*xaib + (caib(im, it, ia-1)+caib(im&
&       , it, ia+1))*tempb3
      tempb4 = .5*ft*xaib
      ftb = ftb + (.5*(caib(im, it+1, ia)*(ft+1.)-caib(im, it-1, ia)*(1.&
&       -ft))-caib(im, it, ia)*2*ft)*xaib + (caib(im, it-1, ia)+caib(im&
&       , it+1, ia))*tempb4
      CALL POPREAL8(xai)
      CALL POPINTEGER4(ia)
      CALL POPINTEGER4(it)
      CALL POPINTEGER4(im)
      fab = fab/da
      CALL POPREAL8(ft)
      taucb = ftb/(dt*tauc*LOG(10.0))
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) taucb = 0.0_8
      CALL POPREAL8(fa)
      tempb2 = fab/cc(kk)
      fcldb(k) = fcldb(k) + tempb2
      ccb(kk) = ccb(kk) - fcld(k)*tempb2/cc(kk)
    ELSE IF (branch .EQ. 1) THEN
      taucb = 0.0_8
    ELSE
      taucld4b = taucld4b + taubeamb(k, 4) + taudiffb(k, 4)
      taudiffb(k, 4) = 0.0_8
      taubeamb(k, 4) = 0.0_8
      taucld3b = taucld3b + taubeamb(k, 3) + taudiffb(k, 3)
      taudiffb(k, 3) = 0.0_8
      taubeamb(k, 3) = 0.0_8
      taucld2b = taucld2b + taubeamb(k, 2) + taudiffb(k, 2)
      taudiffb(k, 2) = 0.0_8
      taubeamb(k, 2) = 0.0_8
      taucld1b = taucld1b + taubeamb(k, 1) + taudiffb(k, 1)
      taudiffb(k, 1) = 0.0_8
      taubeamb(k, 1) = 0.0_8
      GOTO 100
    END IF
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER4(kk)
    ELSE IF (branch .EQ. 1) THEN
      CALL POPINTEGER4(kk)
    ELSE
      CALL POPINTEGER4(kk)
    END IF
 100 CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(taucld4)
    ELSE
      CALL POPREAL8(taucld4)
      tempb1 = dp(k)*aib_uv*1.0e3*taucld4b/(cons_grav*reff_snow)
      hydrometsb(k, 4) = hydrometsb(k, 4) + tempb1
      reff_snowb = reff_snowb - hydromets(k, 4)*tempb1/reff_snow
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(reff_snow)
    ELSE
      CALL POPREAL8(reff_snow)
      reffb(k, 4) = reffb(k, 4) + reff_snowb
    END IF
    CALL POPREAL8(taucld3)
    hydrometsb(k, 3) = hydrometsb(k, 3) + dp(k)*arb_uv(1)*1.0e3*taucld3b&
&     /cons_grav
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(taucld2)
      temp0 = awb_uv(2)/reff(k, 2)
      tempb0 = dp(k)*1.0e3*taucld2b
      hydrometsb(k, 2) = hydrometsb(k, 2) + (awb_uv(1)+temp0)*tempb0/&
&       cons_grav
      reffb(k, 2) = reffb(k, 2) - hydromets(k, 2)*temp0*tempb0/(reff(k, &
&       2)*cons_grav)
    ELSE
      CALL POPREAL8(taucld2)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(taucld1)
    ELSE
      CALL POPREAL8(taucld1)
      temp = cons_grav*reff(k, 1)
      tempb = dp(k)*aib_uv*1.0e3*taucld1b/temp
      hydrometsb(k, 1) = hydrometsb(k, 1) + tempb
      reffb(k, 1) = reffb(k, 1) - hydromets(k, 1)*cons_grav*tempb/temp
    END IF
  END DO
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    DO k=nlevs,icb,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(3)
        ccb(3) = 0.0_8
      END IF
    END DO
    DO k=icb-1,ict,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(2)
        ccb(2) = 0.0_8
      END IF
    END DO
    DO k=ict-1,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(1)
        ccb(1) = 0.0_8
      END IF
    END DO
  END IF
END SUBROUTINE GETVISTAU1_B

!  Differentiation of getnirtau1 in reverse (adjoint) mode:
!   gradient     of useful results: hydromets asycl taudiff fcld
!                ssacl taubeam reff
!   with respect to varying inputs: hydromets asycl fcld ssacl
!                reff
SUBROUTINE GETNIRTAU1_B(ib, nlevs, cosz, dp, fcld, fcldb, reff, reffb, &
& hydromets, hydrometsb, ict, icb, taubeam, taubeamb, taudiff, taudiffb&
& , asycl, asyclb, ssacl, ssaclb, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv&
& , arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, awa_nir, ara_nir, &
& aig_nir, awg_nir, arg_nir, caib, caif, cons_grav)
  IMPLICIT NONE
! !INPUT PARAMETERS:
!  Band number
  INTEGER, INTENT(IN) :: ib
!  Number of levels
  INTEGER, INTENT(IN) :: nlevs
!  Cosine of solar zenith angle
  REAL*8, INTENT(IN) :: cosz
!  Delta pressure in Pa
  REAL*8, INTENT(IN) :: dp(nlevs)
!  Cloud fraction (used sometimes)
  REAL*8, INTENT(IN) :: fcld(nlevs)
  REAL*8 :: fcldb(nlevs)
!  Effective radius (microns)
  REAL*8, INTENT(IN) :: reff(nlevs, 4)
  REAL*8 :: reffb(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL*8, INTENT(IN) :: hydromets(nlevs, 4)
  REAL*8 :: hydrometsb(nlevs, 4)
!  Flags for various uses 
  INTEGER, INTENT(IN) :: ict, icb
  REAL*8, INTENT(IN) :: aig_uv(3), awg_uv(3), arg_uv(3)
  REAL*8, INTENT(IN) :: aib_uv, awb_uv(2), arb_uv(2)
  REAL*8, INTENT(IN) :: aib_nir, awb_nir(3, 2), arb_nir(3, 2)
  REAL*8, INTENT(IN) :: aia_nir(3, 3), awa_nir(3, 3), ara_nir(3, 3)
  REAL*8, INTENT(IN) :: aig_nir(3, 3), awg_nir(3, 3), arg_nir(3, 3)
  REAL*8, INTENT(IN) :: caib(11, 9, 11), caif(9, 11)
  REAL*8, INTENT(IN) :: cons_grav
! !OUTPUT PARAMETERS:
!  Optical depth for beam radiation
  REAL*8 :: taubeam(nlevs, 4)
  REAL*8 :: taubeamb(nlevs, 4)
!  Optical depth for diffuse radiation
  REAL*8 :: taudiff(nlevs, 4)
  REAL*8 :: taudiffb(nlevs, 4)
!  Cloud single scattering albedo
  REAL*8 :: ssacl(nlevs)
  REAL*8 :: ssaclb(nlevs)
!  Cloud asymmetry factor
  REAL*8 :: asycl(nlevs)
  REAL*8 :: asyclb(nlevs)
  INTEGER :: k, in, im, it, ia, kk
  REAL*8 :: fm, ft, fa, xai, tauc, asyclt, ssaclt
  REAL*8 :: ftb, fab, xaib, taucb, asycltb, ssacltb
  REAL*8 :: cc(3)
  REAL*8 :: ccb(3)
  REAL*8 :: taucld1, taucld2, taucld3, taucld4
  REAL*8 :: taucld1b, taucld2b, taucld3b, taucld4b
  REAL*8 :: g1, g2, g3, g4
  REAL*8 :: g1b, g2b, g3b, g4b
  REAL*8 :: w1, w2, w3, w4
  REAL*8 :: w1b, w2b, w3b, w4b
  REAL*8 :: reff_snow
  REAL*8 :: reff_snowb
  INTEGER, PARAMETER :: nm=11, nt=9, na=11
  REAL*8, PARAMETER :: dm=0.1, dt=0.30103, da=0.1, t1=-0.9031
  INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC REAL
  INTEGER :: branch
  REAL*8 :: temp3
  REAL*8 :: temp2
  REAL*8 :: temp1
  REAL*8 :: temp0
  REAL*8 :: tempb9
  REAL*8 :: tempb8
  REAL*8 :: tempb7
  REAL*8 :: tempb6
  REAL*8 :: tempb5
  REAL*8 :: tempb4
  REAL*8 :: tempb3
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  REAL*8 :: tempb
  REAL*8 :: temp
  REAL*8 :: temp6
  REAL*8 :: temp5
  REAL*8 :: temp4
  IF (ict .NE. 0) THEN
!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
    cc = 0.0
    DO k=1,ict-1
      IF (cc(1) .LT. fcld(k)) THEN
        cc(1) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(1) = cc(1)
      END IF
    END DO
    DO k=ict,icb-1
      IF (cc(2) .LT. fcld(k)) THEN
        cc(2) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(2) = cc(2)
      END IF
    END DO
    DO k=icb,nlevs
      IF (cc(3) .LT. fcld(k)) THEN
        cc(3) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(3) = cc(3)
      END IF
    END DO
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
  END IF
!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.) THEN
      CALL PUSHREAL8(taucld1)
      taucld1 = 0.
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL8(taucld1)
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*aib_nir/reff(k, 1)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff(k, 2) .LE. 0.) THEN
      CALL PUSHREAL8(taucld2)
      taucld2 = 0.
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL8(taucld2)
      taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_nir(ib, 1)+&
&       awb_nir(ib, 2)/reff(k, 2))
      CALL PUSHCONTROL1B(0)
    END IF
    CALL PUSHREAL8(taucld3)
    taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_nir(ib, 1)
    IF (reff(k, 4) .GT. 112.0) THEN
      CALL PUSHREAL8(reff_snow)
      reff_snow = 112.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL8(reff_snow)
      reff_snow = reff(k, 4)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff_snow .LE. 0.) THEN
      CALL PUSHREAL8(taucld4)
      taucld4 = 0.
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL8(taucld4)
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*aib_nir/reff_snow
      CALL PUSHCONTROL1B(1)
    END IF
    IF (ict .NE. 0) THEN
!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
      IF (k .LT. ict) THEN
        CALL PUSHINTEGER4(kk)
        kk = 1
        CALL PUSHCONTROL2B(0)
      ELSE IF (k .GE. ict .AND. k .LT. icb) THEN
        CALL PUSHINTEGER4(kk)
        kk = 2
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHINTEGER4(kk)
        kk = 3
        CALL PUSHCONTROL2B(2)
      END IF
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
!-----normalize cloud cover following Eq. (7.8)
        IF (cc(kk) .NE. 0.0) THEN
          CALL PUSHREAL8(fa)
          fa = fcld(k)/cc(kk)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHREAL8(fa)
          fa = 0.0
          CALL PUSHCONTROL1B(0)
        END IF
        IF (tauc .GT. 32.) THEN
          tauc = 32.
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          tauc = tauc
        END IF
        fm = cosz/dm
        CALL PUSHREAL8(ft)
        ft = (LOG10(tauc)-t1)/dt
        fa = fa/da
        CALL PUSHINTEGER4(im)
        im = INT(fm + 1.5)
        CALL PUSHINTEGER4(it)
        it = INT(ft + 1.5)
        CALL PUSHINTEGER4(ia)
        ia = INT(fa + 1.5)
        IF (im .LT. 2) THEN
          im = 2
        ELSE
          im = im
        END IF
        IF (it .LT. 2) THEN
          it = 2
        ELSE
          it = it
        END IF
        IF (ia .LT. 2) THEN
          ia = 2
        ELSE
          ia = ia
        END IF
        IF (im .GT. nm - 1) THEN
          im = nm - 1
        ELSE
          im = im
        END IF
        IF (it .GT. nt - 1) THEN
          it = nt - 1
        ELSE
          it = it
        END IF
        IF (ia .GT. na - 1) THEN
          ia = na - 1
        ELSE
          ia = ia
        END IF
        fm = fm - REAL(im - 1)
        ft = ft - REAL(it - 1)
        fa = fa - REAL(ia - 1)
!-----scale cloud optical thickness for beam radiation following 
!     Eq. (7.3).
!     the scaling factor, xai, is a function of the solar zenith
!     angle, optical thickness, and cloud cover.
        CALL PUSHREAL8(xai)
        xai = (-(caib(im-1, it, ia)*(1.-fm))+caib(im+1, it, ia)*(1.+fm))&
&         *fm*.5 + caib(im, it, ia)*(1.-fm*fm)
        xai = xai + (-(caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(&
&         1.+ft))*ft*.5 + caib(im, it, ia)*(1.-ft*ft)
        xai = xai + (-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(&
&         1.+fa))*fa*.5 + caib(im, it, ia)*(1.-fa*fa)
        xai = xai - 2.*caib(im, it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
        CALL PUSHREAL8(xai)
        xai = (-(caif(it-1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ft*.5 +&
&         caif(it, ia)*(1.-ft*ft)
        xai = xai + (-(caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*&
&         fa*.5 + caif(it, ia)*(1.-fa*fa)
        xai = xai - caif(it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
!-----compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
    CALL PUSHREAL8(tauc)
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      CALL PUSHREAL8(w1)
      w1 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff(k, 1)&
&       )*reff(k, 1)))*taucld1
      CALL PUSHREAL8(w2)
      w2 = (1.-(awa_nir(ib, 1)+(awa_nir(ib, 2)+awa_nir(ib, 3)*reff(k, 2)&
&       )*reff(k, 2)))*taucld2
      CALL PUSHREAL8(w3)
      w3 = (1.-ara_nir(ib, 1))*taucld3
      CALL PUSHREAL8(w4)
      w4 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff_snow)&
&       *reff_snow))*taucld4
      g1 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 1))*&
&       reff(k, 1))*w1
      g2 = (awg_nir(ib, 1)+(awg_nir(ib, 2)+awg_nir(ib, 3)*reff(k, 2))*&
&       reff(k, 2))*w2
      g3 = arg_nir(ib, 1)*w3
      g4 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 4))*&
&       reff(k, 4))*w4
      IF (w1 + w2 + w3 + w4 .NE. 0.0) THEN
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
  END DO
  ccb = 0.0_8
  DO k=nlevs,1,-1
    asycltb = asyclb(k)
    asyclb(k) = 0.0_8
    ssacltb = ssaclb(k)
    ssaclb(k) = 0.0_8
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      w1 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff(k, 1)&
&       )*reff(k, 1)))*taucld1
      w2 = (1.-(awa_nir(ib, 1)+(awa_nir(ib, 2)+awa_nir(ib, 3)*reff(k, 2)&
&       )*reff(k, 2)))*taucld2
      taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_nir(ib, 1)
      w3 = (1.-ara_nir(ib, 1))*taucld3
      w4 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff_snow)&
&       *reff_snow))*taucld4
      g1 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 1))*&
&       reff(k, 1))*w1
      g2 = (awg_nir(ib, 1)+(awg_nir(ib, 2)+awg_nir(ib, 3)*reff(k, 2))*&
&       reff(k, 2))*w2
      g3 = arg_nir(ib, 1)*w3
      g4 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 4))*&
&       reff(k, 4))*w4
      tempb8 = asycltb/(w1+w2+w3+w4)
      tempb9 = -((g1+g2+g3+g4)*tempb8/(w1+w2+w3+w4))
      g1b = tempb8
      g2b = tempb8
      g3b = tempb8
      g4b = tempb8
      w1b = tempb9
      w2b = tempb9
      w3b = tempb9
      w4b = tempb9
    ELSE IF (branch .EQ. 1) THEN
      w1b = 0.0_8
      w2b = 0.0_8
      w3b = 0.0_8
      w4b = 0.0_8
      g1b = 0.0_8
      g2b = 0.0_8
      g3b = 0.0_8
      g4b = 0.0_8
    ELSE
      reff_snowb = 0.0_8
      taucld1b = 0.0_8
      taucld2b = 0.0_8
      taucld3b = 0.0_8
      taucld4b = 0.0_8
      taucb = 0.0_8
      GOTO 100
    END IF
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    tempb7 = ssacltb/tauc
    temp6 = aig_nir(ib, 2) + aig_nir(ib, 3)*reff(k, 4)
    reffb(k, 4) = reffb(k, 4) + (w4*temp6+reff(k, 4)*w4*aig_nir(ib, 3))*&
&     g4b
    w4b = w4b + tempb7 + (aig_nir(ib, 1)+temp6*reff(k, 4))*g4b
    w3b = w3b + tempb7 + arg_nir(ib, 1)*g3b
    temp5 = awg_nir(ib, 2) + awg_nir(ib, 3)*reff(k, 2)
    reffb(k, 2) = reffb(k, 2) + (w2*temp5+reff(k, 2)*w2*awg_nir(ib, 3))*&
&     g2b
    w2b = w2b + tempb7 + (awg_nir(ib, 1)+temp5*reff(k, 2))*g2b
    temp4 = aig_nir(ib, 2) + aig_nir(ib, 3)*reff(k, 1)
    reffb(k, 1) = reffb(k, 1) + (w1*temp4+reff(k, 1)*w1*aig_nir(ib, 3))*&
&     g1b
    w1b = w1b + tempb7 + (aig_nir(ib, 1)+temp4*reff(k, 1))*g1b
    taucb = -((w1+w2+w3+w4)*tempb7/tauc)
    CALL POPREAL8(w4)
    temp3 = aia_nir(ib, 2) + aia_nir(ib, 3)*reff_snow
    reff_snowb = (-(taucld4*temp3)-reff_snow*taucld4*aia_nir(ib, 3))*w4b
    taucld4b = (1.-temp3*reff_snow-aia_nir(ib, 1))*w4b
    CALL POPREAL8(w3)
    taucld3b = (1.-ara_nir(ib, 1))*w3b
    CALL POPREAL8(w2)
    temp2 = awa_nir(ib, 2) + awa_nir(ib, 3)*reff(k, 2)
    reffb(k, 2) = reffb(k, 2) + (-(taucld2*temp2)-reff(k, 2)*taucld2*&
&     awa_nir(ib, 3))*w2b
    taucld2b = (1.-temp2*reff(k, 2)-awa_nir(ib, 1))*w2b
    CALL POPREAL8(w1)
    temp1 = aia_nir(ib, 2) + aia_nir(ib, 3)*reff(k, 1)
    reffb(k, 1) = reffb(k, 1) + (-(taucld1*temp1)-reff(k, 1)*taucld1*&
&     aia_nir(ib, 3))*w1b
    taucld1b = (1.-temp1*reff(k, 1)-aia_nir(ib, 1))*w1b
 100 CALL POPREAL8(tauc)
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      taucld4b = taucld4b + xai*taudiffb(k, 4)
      xaib = taucld4*taudiffb(k, 4)
      taudiffb(k, 4) = 0.0_8
      taucld3b = taucld3b + xai*taudiffb(k, 3)
      xaib = xaib + taucld3*taudiffb(k, 3)
      taudiffb(k, 3) = 0.0_8
      taucld2b = taucld2b + xai*taudiffb(k, 2)
      xaib = xaib + taucld2*taudiffb(k, 2)
      taudiffb(k, 2) = 0.0_8
      taucld1b = taucld1b + xai*taudiffb(k, 1)
      xaib = xaib + taucld1*taudiffb(k, 1)
      taudiffb(k, 1) = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0_8
      tempb5 = .5*fa*xaib
      fab = (.5*(caif(it, ia+1)*(fa+1.)-caif(it, ia-1)*(1.-fa))-caif(it&
&       , ia)*2*fa)*xaib + (caif(it, ia-1)+caif(it, ia+1))*tempb5
      CALL POPREAL8(xai)
      tempb6 = .5*ft*xaib
      ftb = (.5*(caif(it+1, ia)*(ft+1.)-caif(it-1, ia)*(1.-ft))-caif(it&
&       , ia)*2*ft)*xaib + (caif(it-1, ia)+caif(it+1, ia))*tempb6
      taucld4b = taucld4b + xai*taubeamb(k, 4)
      xaib = taucld4*taubeamb(k, 4)
      taubeamb(k, 4) = 0.0_8
      taucld3b = taucld3b + xai*taubeamb(k, 3)
      xaib = xaib + taucld3*taubeamb(k, 3)
      taubeamb(k, 3) = 0.0_8
      taucld2b = taucld2b + xai*taubeamb(k, 2)
      xaib = xaib + taucld2*taubeamb(k, 2)
      taubeamb(k, 2) = 0.0_8
      taucld1b = taucld1b + xai*taubeamb(k, 1)
      xaib = xaib + taucld1*taubeamb(k, 1)
      taubeamb(k, 1) = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0_8
      tempb3 = .5*fa*xaib
      fab = fab + (.5*(caib(im, it, ia+1)*(fa+1.)-caib(im, it, ia-1)*(1.&
&       -fa))-caib(im, it, ia)*2*fa)*xaib + (caib(im, it, ia-1)+caib(im&
&       , it, ia+1))*tempb3
      tempb4 = .5*ft*xaib
      ftb = ftb + (.5*(caib(im, it+1, ia)*(ft+1.)-caib(im, it-1, ia)*(1.&
&       -ft))-caib(im, it, ia)*2*ft)*xaib + (caib(im, it-1, ia)+caib(im&
&       , it+1, ia))*tempb4
      CALL POPREAL8(xai)
      CALL POPINTEGER4(ia)
      CALL POPINTEGER4(it)
      CALL POPINTEGER4(im)
      fab = fab/da
      CALL POPREAL8(ft)
      taucb = ftb/(dt*tauc*LOG(10.0))
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) taucb = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8(fa)
      ELSE
        CALL POPREAL8(fa)
        tempb2 = fab/cc(kk)
        fcldb(k) = fcldb(k) + tempb2
        ccb(kk) = ccb(kk) - fcld(k)*tempb2/cc(kk)
      END IF
    ELSE IF (branch .EQ. 1) THEN
      taucb = 0.0_8
    ELSE
      taucld4b = taucld4b + taubeamb(k, 4) + taudiffb(k, 4)
      taudiffb(k, 4) = 0.0_8
      taubeamb(k, 4) = 0.0_8
      taucld3b = taucld3b + taubeamb(k, 3) + taudiffb(k, 3)
      taudiffb(k, 3) = 0.0_8
      taubeamb(k, 3) = 0.0_8
      taucld2b = taucld2b + taubeamb(k, 2) + taudiffb(k, 2)
      taudiffb(k, 2) = 0.0_8
      taubeamb(k, 2) = 0.0_8
      taucld1b = taucld1b + taubeamb(k, 1) + taudiffb(k, 1)
      taudiffb(k, 1) = 0.0_8
      taubeamb(k, 1) = 0.0_8
      GOTO 110
    END IF
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER4(kk)
    ELSE IF (branch .EQ. 1) THEN
      CALL POPINTEGER4(kk)
    ELSE
      CALL POPINTEGER4(kk)
    END IF
 110 CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(taucld4)
    ELSE
      CALL POPREAL8(taucld4)
      tempb1 = dp(k)*aib_nir*1.0e3*taucld4b/(cons_grav*reff_snow)
      hydrometsb(k, 4) = hydrometsb(k, 4) + tempb1
      reff_snowb = reff_snowb - hydromets(k, 4)*tempb1/reff_snow
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(reff_snow)
    ELSE
      CALL POPREAL8(reff_snow)
      reffb(k, 4) = reffb(k, 4) + reff_snowb
    END IF
    CALL POPREAL8(taucld3)
    hydrometsb(k, 3) = hydrometsb(k, 3) + dp(k)*1.0e3*arb_nir(ib, 1)*&
&     taucld3b/cons_grav
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(taucld2)
      temp0 = awb_nir(ib, 2)/reff(k, 2)
      tempb0 = dp(k)*1.0e3*taucld2b
      hydrometsb(k, 2) = hydrometsb(k, 2) + (awb_nir(ib, 1)+temp0)*&
&       tempb0/cons_grav
      reffb(k, 2) = reffb(k, 2) - hydromets(k, 2)*temp0*tempb0/(reff(k, &
&       2)*cons_grav)
    ELSE
      CALL POPREAL8(taucld2)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(taucld1)
    ELSE
      CALL POPREAL8(taucld1)
      temp = cons_grav*reff(k, 1)
      tempb = dp(k)*aib_nir*1.0e3*taucld1b/temp
      hydrometsb(k, 1) = hydrometsb(k, 1) + tempb
      reffb(k, 1) = reffb(k, 1) - hydromets(k, 1)*cons_grav*tempb/temp
    END IF
  END DO
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    DO k=nlevs,icb,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(3)
        ccb(3) = 0.0_8
      END IF
    END DO
    DO k=icb-1,ict,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(2)
        ccb(2) = 0.0_8
      END IF
    END DO
    DO k=ict-1,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(1)
        ccb(1) = 0.0_8
      END IF
    END DO
  END IF
END SUBROUTINE GETNIRTAU1_B

end module SORAD_AD
