module SORAD_TL

use SORADMOD

IMPLICIT NONE

PRIVATE
PUBLIC :: sorad_d

contains

SUBROUTINE SORAD_D(m, np, nb, cosz_dev, pl_dev, ta_dev, ta_devd, wa_dev&
& , wa_devd, oa_dev, oa_devd, co2, cwc_dev, cwc_devd, fcld_dev, &
& fcld_devd, ict, icb, reff_dev, reff_devd, hk_uv, hk_ir, taua_dev, &
& taua_devd, ssaa_dev, ssaa_devd, asya_dev, asya_devd, rsuvbm_dev, &
& rsuvdf_dev, rsirbm_dev, rsirdf_dev, flx_dev, flx_devd, cons_grav, &
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
  REAL*8 :: ta_devd(m, np), wa_devd(m, np), oa_devd(m, np)
  REAL*8 :: cwc_dev(m, np, 4), fcld_dev(m, np), reff_dev(m, np, 4), &
& hk_uv(5), hk_ir(3, 10)
  REAL*8 :: cwc_devd(m, np, 4), fcld_devd(m, np), reff_devd(m, np, 4)
  REAL*8 :: rsuvbm_dev, rsuvdf_dev, rsirbm_dev, rsirdf_dev
  REAL*8 :: taua_dev(m, np, nb)
  REAL*8 :: taua_devd(m, np, nb)
  REAL*8 :: ssaa_dev(m, np, nb)
  REAL*8 :: ssaa_devd(m, np, nb)
  REAL*8 :: asya_dev(m, np, nb)
  REAL*8 :: asya_devd(m, np, nb)
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
  REAL*8 :: flx_devd(m, np+1)
  REAL*8 :: flxu_dev(m, np+1), flcu_dev(m, np+1)
  REAL*8 :: fdiruv_dev(m), fdifuv_dev(m)
  REAL*8 :: fdirpar_dev(m), fdifpar_dev(m)
  REAL*8 :: fdirir_dev(m), fdifir_dev(m)
  REAL*8 :: flx_sfc_band_dev(m, nband)
!-----temporary arrays
  INTEGER :: i, j, k, l, in, ntop
  REAL*8 :: dp(np), wh(np), oh(np)
  REAL*8 :: dpd(np), whd(np), ohd(np)
  REAL*8 :: scal(np)
  REAL*8 :: scald(np)
  REAL*8 :: swh(np+1), so2(np+1), df(0:np+1)
  REAL*8 :: swhd(np+1), so2d(np+1), dfd(0:np+1)
  REAL*8 :: scal0, wvtoa, o3toa, pa
  REAL*8 :: wvtoad, o3toad
  REAL*8 :: snt, cnt, x, xx4, xtoa
  REAL*8 :: xx4d
  REAL*8 :: dp_pa(np)
  REAL*8 :: dp_pad(np)
!-----parameters for co2 transmission tables
  REAL*8 :: w1, dw, u1, du
  INTEGER :: ib, rc
  REAL*8 :: tauclb(np), tauclf(np), asycl(np)
  REAL*8 :: tauclbd(np), tauclfd(np), asycld(np)
  REAL*8 :: taubeam(np, 4), taudiff(np, 4)
  REAL*8 :: taubeamd(np, 4), taudiffd(np, 4)
  REAL*8 :: fcld_col(np)
  REAL*8 :: fcld_cold(np)
  REAL*8 :: cwc_col(np, 4)
  REAL*8 :: cwc_cold(np, 4)
  REAL*8 :: reff_col(np, 4)
  REAL*8 :: reff_cold(np, 4)
  REAL*8 :: taurs, tauoz, tauwv
  REAL*8 :: tauozd, tauwvd
  REAL*8 :: tausto, ssatau, asysto
  REAL*8 :: taustod, ssataud, asystod
  REAL*8 :: tautob, ssatob, asytob
  REAL*8 :: tautobd, ssatobd, asytobd
  REAL*8 :: tautof, ssatof, asytof
  REAL*8 :: tautofd, ssatofd, asytofd
  REAL*8 :: rr(0:np+1, 2), tt(0:np+1, 2), td(0:np+1, 2)
  REAL*8 :: rrd(0:np+1, 2), ttd(0:np+1, 2), tdd(0:np+1, 2)
  REAL*8 :: rs(0:np+1, 2), ts(0:np+1, 2)
  REAL*8 :: rsd(0:np+1, 2), tsd(0:np+1, 2)
  REAL*8 :: fall(np+1), fclr(np+1), fsdir, fsdif
  REAL*8 :: falld(np+1)
  REAL*8 :: fupa(np+1), fupc(np+1)
  REAL*8 :: cc1, cc2, cc3
  REAL*8 :: cc1d, cc2d, cc3d
  REAL*8 :: rrt, ttt, tdt, rst, tst
  REAL*8 :: rrtd, tttd, tdtd, rstd, tstd
  INTEGER :: iv, ik
  REAL*8 :: ssacl(np)
  REAL*8 :: ssacld(np)
  INTEGER :: im
  INTEGER :: ic, iw
  REAL*8 :: ulog, wlog, dc, dd, x0, x1, x2, y0, y1, y2, du2, dw2
  REAL*8 :: wlogd, ddd, x2d, y2d
  INTEGER :: ih
!if (overcast == true) then
!real(8) :: rra(0:np+1),rxa(0:np+1)
!real(8) :: ttaold,tdaold,rsaold
!real(8) :: ttanew,tdanew,rsanew 
!else
  REAL*8 :: rra(0:np+1, 2, 2), tta(0:np, 2, 2)
  REAL*8 :: rrad(0:np+1, 2, 2), ttad(0:np, 2, 2)
  REAL*8 :: tda(0:np, 2, 2)
  REAL*8 :: tdad(0:np, 2, 2)
  REAL*8 :: rsa(0:np, 2, 2), rxa(0:np+1, 2, 2)
  REAL*8 :: rsad(0:np, 2, 2), rxad(0:np+1, 2, 2)
!endif
  REAL*8 :: flxdn
  REAL*8 :: flxdnd
  REAL*8 :: fdndir, fdndif, fupdif
  REAL*8 :: fdndird, fdndifd, fupdifd
  REAL*8 :: denm, yy
  REAL*8 :: denmd, yyd
  INTEGER :: is
  REAL*8 :: ch, cm, ct
  REAL*8 :: chd, cmd, ctd
  INTEGER :: foundtop
  REAL*8 :: dftop
  REAL*8 :: dftopd
!-----Variables for aerosols
  INTEGER :: ii, jj, irhp1, an
  REAL*8 :: dum
  REAL*8 :: dumd
  INTRINSIC MAX
  INTRINSIC EXP
  INTRINSIC MIN
  INTRINSIC SQRT
  INTRINSIC REAL
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC ABS
  INTRINSIC EPSILON
  REAL*8 :: arg1
  REAL*8 :: arg1d
  REAL*8 :: result1
  REAL :: result10
  REAL*8 :: x6
  REAL*8 :: x5
  REAL*8 :: x4
  REAL*8 :: x3
  REAL*8 :: x4d
  REAL*8 :: abs0
  i = 1
!RUN_LOOP: do i=1,m
  ntop = 0
  fdndir = 0.0
  fdndif = 0.0
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
  o3toad = 1.02*xtoa*466.7*oa_devd(i, 1)
  o3toa = 1.02*oa_dev(i, 1)*xtoa*466.7 + 1.0e-8
  wvtoad = 1.02*scal0*(wa_devd(i, 1)*(1.0+0.00135*(ta_dev(i, 1)-240.))+&
&   wa_dev(i, 1)*0.00135*ta_devd(i, 1))
  wvtoa = 1.02*wa_dev(i, 1)*scal0*(1.0+0.00135*(ta_dev(i, 1)-240.)) + &
&   1.0e-9
  swhd = 0.0_8
  swhd(1) = wvtoad
  swh(1) = wvtoa
  reff_cold = 0.0_8
  ohd = 0.0_8
  fcld_cold = 0.0_8
  cwc_cold = 0.0_8
  whd = 0.0_8
  DO k=1,np
!-----compute layer thickness. indices for the surface level and
!     surface layer are np+1 and np, respectively.
    dpd(k) = 0.0_8
    dp(k) = pl_dev(i, k+1) - pl_dev(i, k)
! dp in pascals
    dp_pad(k) = 0.0_8
    dp_pa(k) = dp(k)*100.
!-----compute scaled water vapor amount following Eqs. (3.3) and (3.5) 
!     unit is g/cm**2
!
    pa = 0.5*(pl_dev(i, k)+pl_dev(i, k+1))
    scald(k) = 0.0_8
    scal(k) = dp(k)*(pa/300.)**.8
    whd(k) = 1.02*scal(k)*(wa_devd(i, k)*(1.+0.00135*(ta_dev(i, k)-240.)&
&     )+wa_dev(i, k)*0.00135*ta_devd(i, k))
    wh(k) = 1.02*wa_dev(i, k)*scal(k)*(1.+0.00135*(ta_dev(i, k)-240.)) +&
&     1.e-9
    swhd(k+1) = swhd(k) + whd(k)
    swh(k+1) = swh(k) + wh(k)
!-----compute ozone amount, unit is (cm-atm)stp
!     the number 466.7 is the unit conversion factor
!     from g/cm**2 to (cm-atm)stp
    ohd(k) = 1.02*dp(k)*466.7*oa_devd(i, k)
    oh(k) = 1.02*oa_dev(i, k)*dp(k)*466.7 + 1.e-8
!-----Fill the reff, cwc, and fcld for the column
    fcld_cold(k) = fcld_devd(i, k)
    fcld_col(k) = fcld_dev(i, k)
    DO l=1,4
      reff_cold(k, l) = reff_devd(i, k, l)
      reff_col(k, l) = reff_dev(i, k, l)
      cwc_cold(k, l) = cwc_devd(i, k, l)
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
    flx_devd(i, k) = 0.0_8
    flx_dev(i, k) = 0.
    flc_dev(i, k) = 0.
    flxu_dev(i, k) = 0.
    flcu_dev(i, k) = 0.
  END DO
!-----Initialize new per-band surface fluxes
  DO ib=1,nband
    flx_sfc_band_dev(i, ib) = 0.
  END DO
!-----Begin inline of SOLUV
!-----compute solar uv and par fluxes
!-----initialize fdiruv, fdifuv, surface reflectances and transmittances.
!     the reflectance and transmittance of the clear and cloudy portions
!     of a layer are denoted by 1 and 2, respectively.
!     cc is the maximum cloud cover in each of the high, middle, and low
!     cloud groups.
!     1/dsm=1/cos(53) = 1.66
  fdiruv_dev(i) = 0.0
  fdifuv_dev(i) = 0.0
  rrd(np+1, 1) = 0.0_8
  rr(np+1, 1) = rsuvbm_dev
  rrd(np+1, 2) = 0.0_8
  rr(np+1, 2) = rsuvbm_dev
  rsd(np+1, 1) = 0.0_8
  rs(np+1, 1) = rsuvdf_dev
  rsd(np+1, 2) = 0.0_8
  rs(np+1, 2) = rsuvdf_dev
  tdd(np+1, 1) = 0.0_8
  td(np+1, 1) = 0.0
  tdd(np+1, 2) = 0.0_8
  td(np+1, 2) = 0.0
  ttd(np+1, 1) = 0.0_8
  tt(np+1, 1) = 0.0
  ttd(np+1, 2) = 0.0_8
  tt(np+1, 2) = 0.0
  tsd(np+1, 1) = 0.0_8
  ts(np+1, 1) = 0.0
  tsd(np+1, 2) = 0.0_8
  ts(np+1, 2) = 0.0
  rrd(0, 1) = 0.0_8
  rr(0, 1) = 0.0
  rrd(0, 2) = 0.0_8
  rr(0, 2) = 0.0
  rsd(0, 1) = 0.0_8
  rs(0, 1) = 0.0
  rsd(0, 2) = 0.0_8
  rs(0, 2) = 0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
  ttd(0, 1) = 0.0_8
  tt(0, 1) = 1.0
  ttd(0, 2) = 0.0_8
  tt(0, 2) = 1.0
  tsd(0, 1) = 0.0_8
  ts(0, 1) = 1.0
  tsd(0, 2) = 0.0_8
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
  CALL GETVISTAU1_D(np, cosz_dev(i), dp_pa, fcld_col, fcld_cold, &
&             reff_col, reff_cold, cwc_col, cwc_cold, ict, icb, taubeam&
&             , taubeamd, taudiff, taudiffd, asycl, asycld, aig_uv, &
&             awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, aib_nir, awb_nir, &
&             arb_nir, aia_nir, awa_nir, ara_nir, aig_nir, awg_nir, &
&             arg_nir, caib, caif, cons_grav)
  cc1d = 0.0_8
  cc2d = 0.0_8
  cc3d = 0.0_8
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
        cc1d = fcld_devd(i, k)
        cc1 = fcld_dev(i, k)
      ELSE
        cc1 = cc1
      END IF
    ELSE IF (k .LT. icb) THEN
      IF (cc2 .LT. fcld_dev(i, k)) THEN
        cc2d = fcld_devd(i, k)
        cc2 = fcld_dev(i, k)
      ELSE
        cc2 = cc2
      END IF
    ELSE IF (cc3 .LT. fcld_dev(i, k)) THEN
      cc3d = fcld_devd(i, k)
      cc3 = fcld_dev(i, k)
    ELSE
      cc3 = cc3
    END IF
  END DO
  tauclbd = 0.0_8
  tauclfd = 0.0_8
!MAT---DO NOT FUSE THIS LOOP
!endif !overcast
  DO k=1,np
    tauclbd(k) = taubeamd(k, 1) + taubeamd(k, 2) + taubeamd(k, 3) + &
&     taubeamd(k, 4)
    tauclb(k) = taubeam(k, 1) + taubeam(k, 2) + taubeam(k, 3) + taubeam(&
&     k, 4)
    tauclfd(k) = taudiffd(k, 1) + taudiffd(k, 2) + taudiffd(k, 3) + &
&     taudiffd(k, 4)
    tauclf(k) = taudiff(k, 1) + taudiff(k, 2) + taudiff(k, 3) + taudiff(&
&     k, 4)
  END DO
  tsd = 0.0_8
  ttd = 0.0_8
  rsad = 0.0_8
  rrd = 0.0_8
  rsd = 0.0_8
  falld = 0.0_8
  ttad = 0.0_8
  rxad = 0.0_8
  tdad = 0.0_8
  rrad = 0.0_8
  tdd = 0.0_8
!-----integration over spectral bands
!-----Compute optical thickness, single-scattering albedo and asymmetry
!     factor for a mixture of "na" aerosol types. [Eqs. (4.16)-(4.18)]
  DO ib=1,nband_uv
!-----compute direct beam transmittances of the layer above pl(1)
    arg1d = -((wk_uv(ib)*wvtoad+zk_uv(ib)*o3toad)/cosz_dev(i))
    arg1 = -((wvtoa*wk_uv(ib)+o3toa*zk_uv(ib))/cosz_dev(i))
    tdd(0, 1) = arg1d*EXP(arg1)
    td(0, 1) = EXP(arg1)
    tdd(0, 2) = tdd(0, 1)
    td(0, 2) = td(0, 1)
    DO k=1,np
!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor (Eqs. 6.2-6.4)
      taurs = ry_uv(ib)*dp(k)
      tauozd = zk_uv(ib)*ohd(k)
      tauoz = zk_uv(ib)*oh(k)
      tauwvd = wk_uv(ib)*whd(k)
      tauwv = wk_uv(ib)*wh(k)
      taustod = tauozd + tauwvd + taua_devd(i, k, ib)
      tausto = taurs + tauoz + tauwv + taua_dev(i, k, ib) + 1.0e-7
      ssataud = ssaa_devd(i, k, ib)
      ssatau = ssaa_dev(i, k, ib) + taurs
      asystod = asya_devd(i, k, ib)
      asysto = asya_dev(i, k, ib)
      tautobd = taustod
      tautob = tausto
      asytobd = (asystod*ssatau-asysto*ssataud)/ssatau**2
      asytob = asysto/ssatau
      ssatobd = (ssataud*tautob-ssatau*tautobd)/tautob**2
      ssatob = ssatau/tautob + 1.0e-8
      IF (ssatob .GT. 0.999999) THEN
        ssatob = 0.999999
        ssatobd = 0.0_8
      ELSE
        ssatob = ssatob
      END IF
!-----for direct incident radiation
      CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, asytobd, &
&             cosz_dev(i), rrt, rrtd, ttt, tttd, tdt, tdtd)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
      CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, asytobd, &
&             dsm, rst, rstd, tst, tstd, dum, dumd)
      rrd(k, 1) = rrtd
      rr(k, 1) = rrt
      ttd(k, 1) = tttd
      tt(k, 1) = ttt
      tdd(k, 1) = tdtd
      td(k, 1) = tdt
      rsd(k, 1) = rstd
      rs(k, 1) = rst
      tsd(k, 1) = tstd
      ts(k, 1) = tst
!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer
!-----for direct incident radiation
!     The effective layer optical properties. Eqs. (6.2)-(6.4)
      tautobd = taustod + tauclbd(k)
      tautob = tausto + tauclb(k)
      ssatobd = ((ssataud+tauclbd(k))*tautob-(ssatau+tauclb(k))*tautobd)&
&       /tautob**2
      ssatob = (ssatau+tauclb(k))/tautob + 1.0e-8
      IF (ssatob .GT. 0.999999) THEN
        ssatob = 0.999999
        ssatobd = 0.0_8
      ELSE
        ssatob = ssatob
      END IF
      asytobd = ((asystod+asycld(k)*tauclb(k)+asycl(k)*tauclbd(k))*&
&       ssatob*tautob-(asysto+asycl(k)*tauclb(k))*(ssatobd*tautob+ssatob&
&       *tautobd))/(ssatob*tautob)**2
      asytob = (asysto+asycl(k)*tauclb(k))/(ssatob*tautob)
!-----for diffuse incident radiation
      tautofd = taustod + tauclfd(k)
      tautof = tausto + tauclf(k)
      ssatofd = ((ssataud+tauclfd(k))*tautof-(ssatau+tauclf(k))*tautofd)&
&       /tautof**2
      ssatof = (ssatau+tauclf(k))/tautof + 1.0e-8
      IF (ssatof .GT. 0.999999) THEN
        ssatof = 0.999999
        ssatofd = 0.0_8
      ELSE
        ssatof = ssatof
      END IF
      asytofd = ((asystod+asycld(k)*tauclf(k)+asycl(k)*tauclfd(k))*&
&       ssatof*tautof-(asysto+asycl(k)*tauclf(k))*(ssatofd*tautof+ssatof&
&       *tautofd))/(ssatof*tautof)**2
      asytof = (asysto+asycl(k)*tauclf(k))/(ssatof*tautof)
!-----for direct incident radiation
!     note that the cloud optical thickness is scaled differently 
!     for direct and diffuse insolation, Eqs. (7.3) and (7.4).
      CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, asytobd, &
&             cosz_dev(i), rrt, rrtd, ttt, tttd, tdt, tdtd)
!-----diffuse incident radiation is approximated by beam radiation 
!     with an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
      CALL DELEDD_D(tautof, tautofd, ssatof, ssatofd, asytof, asytofd, &
&             dsm, rst, rstd, tst, tstd, dum, dumd)
      rrd(k, 2) = rrtd
      rr(k, 2) = rrt
      ttd(k, 2) = tttd
      tt(k, 2) = ttt
      tdd(k, 2) = tdtd
      td(k, 2) = tdt
      rsd(k, 2) = rstd
      rs(k, 2) = rst
      tsd(k, 2) = tstd
      ts(k, 2) = tst
    END DO
!-----flux calculations
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)
    DO k=1,np+1
      fclr(k) = 0.0
      falld(k) = 0.0_8
      fall(k) = 0.0
      fupa(k) = 0.0
      fupc(k) = 0.0
    END DO
    fsdir = 0.0
    fsdif = 0.0
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
      tdad(0, ih, 1) = tdd(0, ih)
      tda(0, ih, 1) = td(0, ih)
      ttad(0, ih, 1) = ttd(0, ih)
      tta(0, ih, 1) = tt(0, ih)
      rsad(0, ih, 1) = rsd(0, ih)
      rsa(0, ih, 1) = rs(0, ih)
      tdad(0, ih, 2) = tdd(0, ih)
      tda(0, ih, 2) = td(0, ih)
      ttad(0, ih, 2) = ttd(0, ih)
      tta(0, ih, 2) = tt(0, ih)
      rsad(0, ih, 2) = rsd(0, ih)
      rsa(0, ih, 2) = rs(0, ih)
      DO k=1,ict-1
        denmd = (tsd(k, ih)*(1.-rsa(k-1, ih, 1)*rs(k, ih))-ts(k, ih)*(-(&
&         rsad(k-1, ih, 1)*rs(k, ih))-rsa(k-1, ih, 1)*rsd(k, ih)))/(1.-&
&         rsa(k-1, ih, 1)*rs(k, ih))**2
        denm = ts(k, ih)/(1.-rsa(k-1, ih, 1)*rs(k, ih))
        tdad(k, ih, 1) = tdad(k-1, ih, 1)*td(k, ih) + tda(k-1, ih, 1)*&
&         tdd(k, ih)
        tda(k, ih, 1) = tda(k-1, ih, 1)*td(k, ih)
        ttad(k, ih, 1) = tdad(k-1, ih, 1)*tt(k, ih) + tda(k-1, ih, 1)*&
&         ttd(k, ih) + ((tdad(k-1, ih, 1)*rr(k, ih)+tda(k-1, ih, 1)*rrd(&
&         k, ih))*rsa(k-1, ih, 1)+tda(k-1, ih, 1)*rr(k, ih)*rsad(k-1, ih&
&         , 1)+ttad(k-1, ih, 1)-tdad(k-1, ih, 1))*denm + (tda(k-1, ih, 1&
&         )*rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1, ih, 1))*&
&         denmd
        tta(k, ih, 1) = tda(k-1, ih, 1)*tt(k, ih) + (tda(k-1, ih, 1)*rsa&
&         (k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1, ih, 1))*denm
        rsad(k, ih, 1) = rsd(k, ih) + (tsd(k, ih)*denm+ts(k, ih)*denmd)*&
&         rsa(k-1, ih, 1) + ts(k, ih)*denm*rsad(k-1, ih, 1)
        rsa(k, ih, 1) = rs(k, ih) + ts(k, ih)*rsa(k-1, ih, 1)*denm
        tdad(k, ih, 2) = tdad(k, ih, 1)
        tda(k, ih, 2) = tda(k, ih, 1)
        ttad(k, ih, 2) = ttad(k, ih, 1)
        tta(k, ih, 2) = tta(k, ih, 1)
        rsad(k, ih, 2) = rsad(k, ih, 1)
        rsa(k, ih, 2) = rsa(k, ih, 1)
      END DO
! k loop
!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition
      DO k=ict,icb-1
        DO im=1,2
          denmd = (tsd(k, im)*(1.-rsa(k-1, ih, im)*rs(k, im))-ts(k, im)*&
&           (-(rsad(k-1, ih, im)*rs(k, im))-rsa(k-1, ih, im)*rsd(k, im))&
&           )/(1.-rsa(k-1, ih, im)*rs(k, im))**2
          denm = ts(k, im)/(1.-rsa(k-1, ih, im)*rs(k, im))
          tdad(k, ih, im) = tdad(k-1, ih, im)*td(k, im) + tda(k-1, ih, &
&           im)*tdd(k, im)
          tda(k, ih, im) = tda(k-1, ih, im)*td(k, im)
          ttad(k, ih, im) = tdad(k-1, ih, im)*tt(k, im) + tda(k-1, ih, &
&           im)*ttd(k, im) + ((tdad(k-1, ih, im)*rr(k, im)+tda(k-1, ih, &
&           im)*rrd(k, im))*rsa(k-1, ih, im)+tda(k-1, ih, im)*rr(k, im)*&
&           rsad(k-1, ih, im)+ttad(k-1, ih, im)-tdad(k-1, ih, im))*denm &
&           + (tda(k-1, ih, im)*rsa(k-1, ih, im)*rr(k, im)+tta(k-1, ih, &
&           im)-tda(k-1, ih, im))*denmd
          tta(k, ih, im) = tda(k-1, ih, im)*tt(k, im) + (tda(k-1, ih, im&
&           )*rsa(k-1, ih, im)*rr(k, im)+tta(k-1, ih, im)-tda(k-1, ih, &
&           im))*denm
          rsad(k, ih, im) = rsd(k, im) + (tsd(k, im)*denm+ts(k, im)*&
&           denmd)*rsa(k-1, ih, im) + ts(k, im)*denm*rsad(k-1, ih, im)
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
      rrad(np+1, 1, is) = rrd(np+1, is)
      rra(np+1, 1, is) = rr(np+1, is)
      rxad(np+1, 1, is) = rsd(np+1, is)
      rxa(np+1, 1, is) = rs(np+1, is)
      rrad(np+1, 2, is) = rrd(np+1, is)
      rra(np+1, 2, is) = rr(np+1, is)
      rxad(np+1, 2, is) = rsd(np+1, is)
      rxa(np+1, 2, is) = rs(np+1, is)
      DO k=np,icb,-1
        denmd = (tsd(k, is)*(1.-rs(k, is)*rxa(k+1, 1, is))-ts(k, is)*(-(&
&         rsd(k, is)*rxa(k+1, 1, is))-rs(k, is)*rxad(k+1, 1, is)))/(1.-&
&         rs(k, is)*rxa(k+1, 1, is))**2
        denm = ts(k, is)/(1.-rs(k, is)*rxa(k+1, 1, is))
        rrad(k, 1, is) = rrd(k, is) + (tdd(k, is)*rra(k+1, 1, is)+td(k, &
&         is)*rrad(k+1, 1, is)+(ttd(k, is)-tdd(k, is))*rxa(k+1, 1, is)+(&
&         tt(k, is)-td(k, is))*rxad(k+1, 1, is))*denm + (td(k, is)*rra(k&
&         +1, 1, is)+(tt(k, is)-td(k, is))*rxa(k+1, 1, is))*denmd
        rra(k, 1, is) = rr(k, is) + (td(k, is)*rra(k+1, 1, is)+(tt(k, is&
&         )-td(k, is))*rxa(k+1, 1, is))*denm
        rxad(k, 1, is) = rsd(k, is) + (tsd(k, is)*denm+ts(k, is)*denmd)*&
&         rxa(k+1, 1, is) + ts(k, is)*denm*rxad(k+1, 1, is)
        rxa(k, 1, is) = rs(k, is) + ts(k, is)*rxa(k+1, 1, is)*denm
        rrad(k, 2, is) = rrad(k, 1, is)
        rra(k, 2, is) = rra(k, 1, is)
        rxad(k, 2, is) = rxad(k, 1, is)
        rxa(k, 2, is) = rxa(k, 1, is)
      END DO
! k loop
!-----for middle clouds
      DO k=icb-1,ict,-1
        DO im=1,2
          denmd = (tsd(k, im)*(1.-rs(k, im)*rxa(k+1, im, is))-ts(k, im)*&
&           (-(rsd(k, im)*rxa(k+1, im, is))-rs(k, im)*rxad(k+1, im, is))&
&           )/(1.-rs(k, im)*rxa(k+1, im, is))**2
          denm = ts(k, im)/(1.-rs(k, im)*rxa(k+1, im, is))
          rrad(k, im, is) = rrd(k, im) + (tdd(k, im)*rra(k+1, im, is)+td&
&           (k, im)*rrad(k+1, im, is)+(ttd(k, im)-tdd(k, im))*rxa(k+1, &
&           im, is)+(tt(k, im)-td(k, im))*rxad(k+1, im, is))*denm + (td(&
&           k, im)*rra(k+1, im, is)+(tt(k, im)-td(k, im))*rxa(k+1, im, &
&           is))*denmd
          rra(k, im, is) = rr(k, im) + (td(k, im)*rra(k+1, im, is)+(tt(k&
&           , im)-td(k, im))*rxa(k+1, im, is))*denm
          rxad(k, im, is) = rsd(k, im) + (tsd(k, im)*denm+ts(k, im)*&
&           denmd)*rxa(k+1, im, is) + ts(k, im)*denm*rxad(k+1, im, is)
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
        chd = -cc1d
        ch = 1.0 - cc1
!-----cloudy portion
      ELSE
        chd = cc1d
        ch = cc1
      END IF
      DO im=1,2
!-----clear portion 
        IF (im .EQ. 1) THEN
          cmd = chd*(1.0-cc2) - ch*cc2d
          cm = ch*(1.0-cc2)
!-----cloudy portion
        ELSE
          cmd = chd*cc2 + ch*cc2d
          cm = ch*cc2
        END IF
        DO is=1,2
!-----clear portion 
          IF (is .EQ. 1) THEN
            ctd = cmd*(1.0-cc3) - cm*cc3d
            ct = cm*(1.0-cc3)
!-----cloudy portion
          ELSE
            ctd = cmd*cc3 + cm*cc3d
            ct = cm*cc3
          END IF
!-----add one layer at a time, going down.
          DO k=icb,np
            denmd = (tsd(k, is)*(1.-rsa(k-1, ih, im)*rs(k, is))-ts(k, is&
&             )*(-(rsad(k-1, ih, im)*rs(k, is))-rsa(k-1, ih, im)*rsd(k, &
&             is)))/(1.-rsa(k-1, ih, im)*rs(k, is))**2
            denm = ts(k, is)/(1.-rsa(k-1, ih, im)*rs(k, is))
            tdad(k, ih, im) = tdad(k-1, ih, im)*td(k, is) + tda(k-1, ih&
&             , im)*tdd(k, is)
            tda(k, ih, im) = tda(k-1, ih, im)*td(k, is)
            ttad(k, ih, im) = tdad(k-1, ih, im)*tt(k, is) + tda(k-1, ih&
&             , im)*ttd(k, is) + ((tdad(k-1, ih, im)*rr(k, is)+tda(k-1, &
&             ih, im)*rrd(k, is))*rsa(k-1, ih, im)+tda(k-1, ih, im)*rr(k&
&             , is)*rsad(k-1, ih, im)+ttad(k-1, ih, im)-tdad(k-1, ih, im&
&             ))*denm + (tda(k-1, ih, im)*rr(k, is)*rsa(k-1, ih, im)+tta&
&             (k-1, ih, im)-tda(k-1, ih, im))*denmd
            tta(k, ih, im) = tda(k-1, ih, im)*tt(k, is) + (tda(k-1, ih, &
&             im)*rr(k, is)*rsa(k-1, ih, im)+tta(k-1, ih, im)-tda(k-1, &
&             ih, im))*denm
            rsad(k, ih, im) = rsd(k, is) + (tsd(k, is)*denm+ts(k, is)*&
&             denmd)*rsa(k-1, ih, im) + ts(k, is)*denm*rsad(k-1, ih, im)
            rsa(k, ih, im) = rs(k, is) + ts(k, is)*rsa(k-1, ih, im)*denm
          END DO
! k loop
!-----add one layer at a time, going up.
          DO k=ict-1,0,-1
            denmd = (tsd(k, ih)*(1.-rs(k, ih)*rxa(k+1, im, is))-ts(k, ih&
&             )*(-(rsd(k, ih)*rxa(k+1, im, is))-rs(k, ih)*rxad(k+1, im, &
&             is)))/(1.-rs(k, ih)*rxa(k+1, im, is))**2
            denm = ts(k, ih)/(1.-rs(k, ih)*rxa(k+1, im, is))
            rrad(k, im, is) = rrd(k, ih) + (tdd(k, ih)*rra(k+1, im, is)+&
&             td(k, ih)*rrad(k+1, im, is)+(ttd(k, ih)-tdd(k, ih))*rxa(k+&
&             1, im, is)+(tt(k, ih)-td(k, ih))*rxad(k+1, im, is))*denm +&
&             (td(k, ih)*rra(k+1, im, is)+(tt(k, ih)-td(k, ih))*rxa(k+1&
&             , im, is))*denmd
            rra(k, im, is) = rr(k, ih) + (td(k, ih)*rra(k+1, im, is)+(tt&
&             (k, ih)-td(k, ih))*rxa(k+1, im, is))*denm
            rxad(k, im, is) = rsd(k, ih) + (tsd(k, ih)*denm+ts(k, ih)*&
&             denmd)*rxa(k+1, im, is) + ts(k, ih)*denm*rxad(k+1, im, is)
            rxa(k, im, is) = rs(k, ih) + ts(k, ih)*rxa(k+1, im, is)*denm
          END DO
! k loop
!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux
          DO k=1,np+1
            denmd = -((-(rsad(k-1, ih, im)*rxa(k, im, is))-rsa(k-1, ih, &
&             im)*rxad(k, im, is))/(1.-rsa(k-1, ih, im)*rxa(k, im, is))&
&             **2)
            denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
            fdndird = tdad(k-1, ih, im)
            fdndir = tda(k-1, ih, im)
            xx4d = tdad(k-1, ih, im)*rra(k, im, is) + tda(k-1, ih, im)*&
&             rrad(k, im, is)
            xx4 = tda(k-1, ih, im)*rra(k, im, is)
            yyd = ttad(k-1, ih, im) - tdad(k-1, ih, im)
            yy = tta(k-1, ih, im) - tda(k-1, ih, im)
            fdndifd = (xx4d*rsa(k-1, ih, im)+xx4*rsad(k-1, ih, im)+yyd)*&
&             denm + (xx4*rsa(k-1, ih, im)+yy)*denmd
            fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
            fupdifd = (xx4d+yyd*rxa(k, im, is)+yy*rxad(k, im, is))*denm &
&             + (xx4+yy*rxa(k, im, is))*denmd
            fupdif = (xx4+yy*rxa(k, im, is))*denm
            flxdnd = fdndird + fdndifd - fupdifd
            flxdn = fdndir + fdndif - fupdif
!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)
            IF (ih .EQ. 1 .AND. im .EQ. 1 .AND. is .EQ. 1) THEN
              fupc(k) = fupdif
              fclr(k) = flxdn
            END IF
            fupa(k) = fupa(k) + fupdif*ct
            falld(k) = falld(k) + flxdnd*ct + flxdn*ctd
            fall(k) = fall(k) + flxdn*ct
          END DO
! k loop
          fsdir = fsdir + fdndir*ct
          fsdif = fsdif + fdndif*ct
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
      flx_devd(i, k) = flx_devd(i, k) + hk_uv(ib)*falld(k)
      flx_dev(i, k) = flx_dev(i, k) + fall(k)*hk_uv(ib)
      flc_dev(i, k) = flc_dev(i, k) + fclr(k)*hk_uv(ib)
      flxu_dev(i, k) = flxu_dev(i, k) + fupa(k)*hk_uv(ib)
      flcu_dev(i, k) = flcu_dev(i, k) + fupc(k)*hk_uv(ib)
    END DO
!-----get surface flux for each band
    flx_sfc_band_dev(i, ib) = flx_sfc_band_dev(i, ib) + fall(np+1)*hk_uv&
&     (ib)
!-----compute direct and diffuse downward surface fluxes in the UV
!     and par regions
    IF (ib .LT. 5) THEN
      fdiruv_dev(i) = fdiruv_dev(i) + fsdir*hk_uv(ib)
      fdifuv_dev(i) = fdifuv_dev(i) + fsdif*hk_uv(ib)
    ELSE
      fdirpar_dev(i) = fsdir*hk_uv(ib)
      fdifpar_dev(i) = fsdif*hk_uv(ib)
    END IF
  END DO
!-----Inline SOLIR
!-----compute and update solar ir fluxes
  fdirir_dev(i) = 0.0
  fdifir_dev(i) = 0.0
  rrd(np+1, 1) = 0.0_8
  rr(np+1, 1) = rsirbm_dev
  rrd(np+1, 2) = 0.0_8
  rr(np+1, 2) = rsirbm_dev
  rsd(np+1, 1) = 0.0_8
  rs(np+1, 1) = rsirdf_dev
  rsd(np+1, 2) = 0.0_8
  rs(np+1, 2) = rsirdf_dev
  tdd(np+1, 1) = 0.0_8
  td(np+1, 1) = 0.0
  tdd(np+1, 2) = 0.0_8
  td(np+1, 2) = 0.0
  ttd(np+1, 1) = 0.0_8
  tt(np+1, 1) = 0.0
  ttd(np+1, 2) = 0.0_8
  tt(np+1, 2) = 0.0
  tsd(np+1, 1) = 0.0_8
  ts(np+1, 1) = 0.0
  tsd(np+1, 2) = 0.0_8
  ts(np+1, 2) = 0.0
  rrd(0, 1) = 0.0_8
  rr(0, 1) = 0.0
  rrd(0, 2) = 0.0_8
  rr(0, 2) = 0.0
  rsd(0, 1) = 0.0_8
  rs(0, 1) = 0.0
  rsd(0, 2) = 0.0_8
  rs(0, 2) = 0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
  ttd(0, 1) = 0.0_8
  tt(0, 1) = 1.0
  ttd(0, 2) = 0.0_8
  tt(0, 2) = 1.0
  tsd(0, 1) = 0.0_8
  ts(0, 1) = 1.0
  tsd(0, 2) = 0.0_8
  ts(0, 2) = 1.0
  cc1 = 0.0
  cc2 = 0.0
  cc3 = 0.0
  cc1d = 0.0_8
  cc2d = 0.0_8
  cc3d = 0.0_8
  ssacld = 0.0_8
!-----integration over spectral bands
!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10)
!     The indices 1, 2, 3 are for ice, water, rain particles,
!     respectively.
  DO ib=1,nband_ir
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
    CALL GETNIRTAU1_D(ib, np, cosz_dev(i), dp_pa, fcld_col, fcld_cold, &
&               reff_col, reff_cold, cwc_col, cwc_cold, ict, icb, &
&               taubeam, taubeamd, taudiff, taudiffd, asycl, asycld, &
&               ssacl, ssacld, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, &
&               arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, awa_nir, &
&               ara_nir, aig_nir, awg_nir, arg_nir, caib, caif, &
&               cons_grav)
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
!MAT--DO NOT FUSE THIS LOOP
!MAT  Loop must run to completion so that cc[1,2,3] are correct.
    DO k=1,np
      IF (k .LT. ict) THEN
        IF (cc1 .LT. fcld_dev(i, k)) THEN
          cc1d = fcld_devd(i, k)
          cc1 = fcld_dev(i, k)
        ELSE
          cc1 = cc1
        END IF
      ELSE IF (k .LT. icb) THEN
        IF (cc2 .LT. fcld_dev(i, k)) THEN
          cc2d = fcld_devd(i, k)
          cc2 = fcld_dev(i, k)
        ELSE
          cc2 = cc2
        END IF
      ELSE IF (cc3 .LT. fcld_dev(i, k)) THEN
        cc3d = fcld_devd(i, k)
        cc3 = fcld_dev(i, k)
      ELSE
        cc3 = cc3
      END IF
    END DO
!MAT--DO NOT FUSE THIS LOOP
!endif !overcast
    DO k=1,np
      tauclbd(k) = taubeamd(k, 1) + taubeamd(k, 2) + taubeamd(k, 3) + &
&       taubeamd(k, 4)
      tauclb(k) = taubeam(k, 1) + taubeam(k, 2) + taubeam(k, 3) + &
&       taubeam(k, 4)
      tauclfd(k) = taudiffd(k, 1) + taudiffd(k, 2) + taudiffd(k, 3) + &
&       taudiffd(k, 4)
      tauclf(k) = taudiff(k, 1) + taudiff(k, 2) + taudiff(k, 3) + &
&       taudiff(k, 4)
    END DO
!-----integration over the k-distribution function
    DO ik=1,nk_ir
!-----compute direct beam transmittances of the layer above pl(1)
      tdd(0, 1) = -(xk_ir(ik)*wvtoad*EXP(-(wvtoa*xk_ir(ik)/cosz_dev(i)))&
&       /cosz_dev(i))
      td(0, 1) = EXP(-(wvtoa*xk_ir(ik)/cosz_dev(i)))
      tdd(0, 2) = tdd(0, 1)
      td(0, 2) = td(0, 1)
      DO k=1,np
        taurs = ry_ir(ib)*dp(k)
        tauwvd = xk_ir(ik)*whd(k)
        tauwv = xk_ir(ik)*wh(k)
!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor. Eqs.(6.2)-(6.4)
        taustod = tauwvd + taua_devd(i, k, iv)
        tausto = taurs + tauwv + taua_dev(i, k, iv) + 1.0e-7
        ssataud = ssaa_devd(i, k, iv)
        ssatau = ssaa_dev(i, k, iv) + taurs + 1.0e-8
        asystod = asya_devd(i, k, iv)
        asysto = asya_dev(i, k, iv)
        tautobd = taustod
        tautob = tausto
        asytobd = (asystod*ssatau-asysto*ssataud)/ssatau**2
        asytob = asysto/ssatau
        ssatobd = (ssataud*tautob-ssatau*tautobd)/tautob**2
        ssatob = ssatau/tautob + 1.0e-8
        IF (ssatob .GT. 0.999999) THEN
          ssatob = 0.999999
          ssatobd = 0.0_8
        ELSE
          ssatob = ssatob
        END IF
!-----Compute reflectance and transmittance of the clear portion 
!     of a layer
!-----for direct incident radiation
        CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, asytobd&
&               , cosz_dev(i), rrt, rrtd, ttt, tttd, tdt, tdtd)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
        CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, asytobd&
&               , dsm, rst, rstd, tst, tstd, dum, dumd)
        rrd(k, 1) = rrtd
        rr(k, 1) = rrt
        ttd(k, 1) = tttd
        tt(k, 1) = ttt
        tdd(k, 1) = tdtd
        td(k, 1) = tdt
        rsd(k, 1) = rstd
        rs(k, 1) = rst
        tsd(k, 1) = tstd
        ts(k, 1) = tst
!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer
!-----for direct incident radiation. Eqs.(6.2)-(6.4)
        tautobd = taustod + tauclbd(k)
        tautob = tausto + tauclb(k)
        ssatobd = ((ssataud+ssacld(k)*tauclb(k)+ssacl(k)*tauclbd(k))*&
&         tautob-(ssatau+ssacl(k)*tauclb(k))*tautobd)/tautob**2
        ssatob = (ssatau+ssacl(k)*tauclb(k))/tautob + 1.0e-8
        IF (ssatob .GT. 0.999999) THEN
          ssatob = 0.999999
          ssatobd = 0.0_8
        ELSE
          ssatob = ssatob
        END IF
        asytobd = ((asystod+(asycld(k)*ssacl(k)+asycl(k)*ssacld(k))*&
&         tauclb(k)+asycl(k)*ssacl(k)*tauclbd(k))*ssatob*tautob-(asysto+&
&         asycl(k)*ssacl(k)*tauclb(k))*(ssatobd*tautob+ssatob*tautobd))/&
&         (ssatob*tautob)**2
        asytob = (asysto+asycl(k)*ssacl(k)*tauclb(k))/(ssatob*tautob)
!-----for diffuse incident radiation
        tautofd = taustod + tauclfd(k)
        tautof = tausto + tauclf(k)
        ssatofd = ((ssataud+ssacld(k)*tauclf(k)+ssacl(k)*tauclfd(k))*&
&         tautof-(ssatau+ssacl(k)*tauclf(k))*tautofd)/tautof**2
        ssatof = (ssatau+ssacl(k)*tauclf(k))/tautof + 1.0e-8
        IF (ssatof .GT. 0.999999) THEN
          ssatof = 0.999999
          ssatofd = 0.0_8
        ELSE
          ssatof = ssatof
        END IF
        asytofd = ((asystod+(asycld(k)*ssacl(k)+asycl(k)*ssacld(k))*&
&         tauclf(k)+asycl(k)*ssacl(k)*tauclfd(k))*ssatof*tautof-(asysto+&
&         asycl(k)*ssacl(k)*tauclf(k))*(ssatofd*tautof+ssatof*tautofd))/&
&         (ssatof*tautof)**2
        asytof = (asysto+asycl(k)*ssacl(k)*tauclf(k))/(ssatof*tautof)
!-----for direct incident radiation
        CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, asytobd&
&               , cosz_dev(i), rrt, rrtd, ttt, tttd, tdt, tdtd)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs.(6.5) and (6.6)
        CALL DELEDD_D(tautof, tautofd, ssatof, ssatofd, asytof, asytofd&
&               , dsm, rst, rstd, tst, tstd, dum, dumd)
        rrd(k, 2) = rrtd
        rr(k, 2) = rrt
        ttd(k, 2) = tttd
        tt(k, 2) = ttt
        tdd(k, 2) = tdtd
        td(k, 2) = tdt
        rsd(k, 2) = rstd
        rs(k, 2) = rst
        tsd(k, 2) = tstd
        ts(k, 2) = tst
      END DO
!-----FLUX CALCULATIONS
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)
      DO k=1,np+1
        fclr(k) = 0.0
        falld(k) = 0.0_8
        fall(k) = 0.0
        fupc(k) = 0.0
        fupa(k) = 0.0
      END DO
      fsdir = 0.0
      fsdif = 0.0
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
        tdad(0, ih, 1) = tdd(0, ih)
        tda(0, ih, 1) = td(0, ih)
        ttad(0, ih, 1) = ttd(0, ih)
        tta(0, ih, 1) = tt(0, ih)
        rsad(0, ih, 1) = rsd(0, ih)
        rsa(0, ih, 1) = rs(0, ih)
        tdad(0, ih, 2) = tdd(0, ih)
        tda(0, ih, 2) = td(0, ih)
        ttad(0, ih, 2) = ttd(0, ih)
        tta(0, ih, 2) = tt(0, ih)
        rsad(0, ih, 2) = rsd(0, ih)
        rsa(0, ih, 2) = rs(0, ih)
        DO k=1,ict-1
          denmd = (tsd(k, ih)*(1.-rsa(k-1, ih, 1)*rs(k, ih))-ts(k, ih)*(&
&           -(rsad(k-1, ih, 1)*rs(k, ih))-rsa(k-1, ih, 1)*rsd(k, ih)))/(&
&           1.-rsa(k-1, ih, 1)*rs(k, ih))**2
          denm = ts(k, ih)/(1.-rsa(k-1, ih, 1)*rs(k, ih))
          tdad(k, ih, 1) = tdad(k-1, ih, 1)*td(k, ih) + tda(k-1, ih, 1)*&
&           tdd(k, ih)
          tda(k, ih, 1) = tda(k-1, ih, 1)*td(k, ih)
          ttad(k, ih, 1) = tdad(k-1, ih, 1)*tt(k, ih) + tda(k-1, ih, 1)*&
&           ttd(k, ih) + ((tdad(k-1, ih, 1)*rr(k, ih)+tda(k-1, ih, 1)*&
&           rrd(k, ih))*rsa(k-1, ih, 1)+tda(k-1, ih, 1)*rr(k, ih)*rsad(k&
&           -1, ih, 1)+ttad(k-1, ih, 1)-tdad(k-1, ih, 1))*denm + (tda(k-&
&           1, ih, 1)*rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1&
&           , ih, 1))*denmd
          tta(k, ih, 1) = tda(k-1, ih, 1)*tt(k, ih) + (tda(k-1, ih, 1)*&
&           rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1, ih, 1))*&
&           denm
          rsad(k, ih, 1) = rsd(k, ih) + (tsd(k, ih)*denm+ts(k, ih)*denmd&
&           )*rsa(k-1, ih, 1) + ts(k, ih)*denm*rsad(k-1, ih, 1)
          rsa(k, ih, 1) = rs(k, ih) + ts(k, ih)*rsa(k-1, ih, 1)*denm
          tdad(k, ih, 2) = tdad(k, ih, 1)
          tda(k, ih, 2) = tda(k, ih, 1)
          ttad(k, ih, 2) = ttad(k, ih, 1)
          tta(k, ih, 2) = tta(k, ih, 1)
          rsad(k, ih, 2) = rsad(k, ih, 1)
          rsa(k, ih, 2) = rsa(k, ih, 1)
        END DO
! k loop
!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition
        DO k=ict,icb-1
          DO im=1,2
            denmd = (tsd(k, im)*(1.-rsa(k-1, ih, im)*rs(k, im))-ts(k, im&
&             )*(-(rsad(k-1, ih, im)*rs(k, im))-rsa(k-1, ih, im)*rsd(k, &
&             im)))/(1.-rsa(k-1, ih, im)*rs(k, im))**2
            denm = ts(k, im)/(1.-rsa(k-1, ih, im)*rs(k, im))
            tdad(k, ih, im) = tdad(k-1, ih, im)*td(k, im) + tda(k-1, ih&
&             , im)*tdd(k, im)
            tda(k, ih, im) = tda(k-1, ih, im)*td(k, im)
            ttad(k, ih, im) = tdad(k-1, ih, im)*tt(k, im) + tda(k-1, ih&
&             , im)*ttd(k, im) + ((tdad(k-1, ih, im)*rr(k, im)+tda(k-1, &
&             ih, im)*rrd(k, im))*rsa(k-1, ih, im)+tda(k-1, ih, im)*rr(k&
&             , im)*rsad(k-1, ih, im)+ttad(k-1, ih, im)-tdad(k-1, ih, im&
&             ))*denm + (tda(k-1, ih, im)*rsa(k-1, ih, im)*rr(k, im)+tta&
&             (k-1, ih, im)-tda(k-1, ih, im))*denmd
            tta(k, ih, im) = tda(k-1, ih, im)*tt(k, im) + (tda(k-1, ih, &
&             im)*rsa(k-1, ih, im)*rr(k, im)+tta(k-1, ih, im)-tda(k-1, &
&             ih, im))*denm
            rsad(k, ih, im) = rsd(k, im) + (tsd(k, im)*denm+ts(k, im)*&
&             denmd)*rsa(k-1, ih, im) + ts(k, im)*denm*rsad(k-1, ih, im)
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
        rrad(np+1, 1, is) = rrd(np+1, is)
        rra(np+1, 1, is) = rr(np+1, is)
        rxad(np+1, 1, is) = rsd(np+1, is)
        rxa(np+1, 1, is) = rs(np+1, is)
        rrad(np+1, 2, is) = rrd(np+1, is)
        rra(np+1, 2, is) = rr(np+1, is)
        rxad(np+1, 2, is) = rsd(np+1, is)
        rxa(np+1, 2, is) = rs(np+1, is)
        DO k=np,icb,-1
          denmd = (tsd(k, is)*(1.-rs(k, is)*rxa(k+1, 1, is))-ts(k, is)*(&
&           -(rsd(k, is)*rxa(k+1, 1, is))-rs(k, is)*rxad(k+1, 1, is)))/(&
&           1.-rs(k, is)*rxa(k+1, 1, is))**2
          denm = ts(k, is)/(1.-rs(k, is)*rxa(k+1, 1, is))
          rrad(k, 1, is) = rrd(k, is) + (tdd(k, is)*rra(k+1, 1, is)+td(k&
&           , is)*rrad(k+1, 1, is)+(ttd(k, is)-tdd(k, is))*rxa(k+1, 1, &
&           is)+(tt(k, is)-td(k, is))*rxad(k+1, 1, is))*denm + (td(k, is&
&           )*rra(k+1, 1, is)+(tt(k, is)-td(k, is))*rxa(k+1, 1, is))*&
&           denmd
          rra(k, 1, is) = rr(k, is) + (td(k, is)*rra(k+1, 1, is)+(tt(k, &
&           is)-td(k, is))*rxa(k+1, 1, is))*denm
          rxad(k, 1, is) = rsd(k, is) + (tsd(k, is)*denm+ts(k, is)*denmd&
&           )*rxa(k+1, 1, is) + ts(k, is)*denm*rxad(k+1, 1, is)
          rxa(k, 1, is) = rs(k, is) + ts(k, is)*rxa(k+1, 1, is)*denm
          rrad(k, 2, is) = rrad(k, 1, is)
          rra(k, 2, is) = rra(k, 1, is)
          rxad(k, 2, is) = rxad(k, 1, is)
          rxa(k, 2, is) = rxa(k, 1, is)
        END DO
! k loop
!-----for middle clouds
        DO k=icb-1,ict,-1
          DO im=1,2
            denmd = (tsd(k, im)*(1.-rs(k, im)*rxa(k+1, im, is))-ts(k, im&
&             )*(-(rsd(k, im)*rxa(k+1, im, is))-rs(k, im)*rxad(k+1, im, &
&             is)))/(1.-rs(k, im)*rxa(k+1, im, is))**2
            denm = ts(k, im)/(1.-rs(k, im)*rxa(k+1, im, is))
            rrad(k, im, is) = rrd(k, im) + (tdd(k, im)*rra(k+1, im, is)+&
&             td(k, im)*rrad(k+1, im, is)+(ttd(k, im)-tdd(k, im))*rxa(k+&
&             1, im, is)+(tt(k, im)-td(k, im))*rxad(k+1, im, is))*denm +&
&             (td(k, im)*rra(k+1, im, is)+(tt(k, im)-td(k, im))*rxa(k+1&
&             , im, is))*denmd
            rra(k, im, is) = rr(k, im) + (td(k, im)*rra(k+1, im, is)+(tt&
&             (k, im)-td(k, im))*rxa(k+1, im, is))*denm
            rxad(k, im, is) = rsd(k, im) + (tsd(k, im)*denm+ts(k, im)*&
&             denmd)*rxa(k+1, im, is) + ts(k, im)*denm*rxad(k+1, im, is)
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
          chd = -cc1d
          ch = 1.0 - cc1
!-----cloudy portion
        ELSE
          chd = cc1d
          ch = cc1
        END IF
        DO im=1,2
!-----clear portion 
          IF (im .EQ. 1) THEN
            cmd = chd*(1.0-cc2) - ch*cc2d
            cm = ch*(1.0-cc2)
!-----cloudy portion
          ELSE
            cmd = chd*cc2 + ch*cc2d
            cm = ch*cc2
          END IF
          DO is=1,2
!-----clear portion 
            IF (is .EQ. 1) THEN
              ctd = cmd*(1.0-cc3) - cm*cc3d
              ct = cm*(1.0-cc3)
!-----cloudy portion
            ELSE
              ctd = cmd*cc3 + cm*cc3d
              ct = cm*cc3
            END IF
!-----add one layer at a time, going down.
            DO k=icb,np
              denmd = (tsd(k, is)*(1.-rsa(k-1, ih, im)*rs(k, is))-ts(k, &
&               is)*(-(rsad(k-1, ih, im)*rs(k, is))-rsa(k-1, ih, im)*rsd&
&               (k, is)))/(1.-rsa(k-1, ih, im)*rs(k, is))**2
              denm = ts(k, is)/(1.-rsa(k-1, ih, im)*rs(k, is))
              tdad(k, ih, im) = tdad(k-1, ih, im)*td(k, is) + tda(k-1, &
&               ih, im)*tdd(k, is)
              tda(k, ih, im) = tda(k-1, ih, im)*td(k, is)
              ttad(k, ih, im) = tdad(k-1, ih, im)*tt(k, is) + tda(k-1, &
&               ih, im)*ttd(k, is) + ((tdad(k-1, ih, im)*rr(k, is)+tda(k&
&               -1, ih, im)*rrd(k, is))*rsa(k-1, ih, im)+tda(k-1, ih, im&
&               )*rr(k, is)*rsad(k-1, ih, im)+ttad(k-1, ih, im)-tdad(k-1&
&               , ih, im))*denm + (tda(k-1, ih, im)*rr(k, is)*rsa(k-1, &
&               ih, im)+tta(k-1, ih, im)-tda(k-1, ih, im))*denmd
              tta(k, ih, im) = tda(k-1, ih, im)*tt(k, is) + (tda(k-1, ih&
&               , im)*rr(k, is)*rsa(k-1, ih, im)+tta(k-1, ih, im)-tda(k-&
&               1, ih, im))*denm
              rsad(k, ih, im) = rsd(k, is) + (tsd(k, is)*denm+ts(k, is)*&
&               denmd)*rsa(k-1, ih, im) + ts(k, is)*denm*rsad(k-1, ih, &
&               im)
              rsa(k, ih, im) = rs(k, is) + ts(k, is)*rsa(k-1, ih, im)*&
&               denm
            END DO
! k loop
!-----add one layer at a time, going up.
            DO k=ict-1,0,-1
              denmd = (tsd(k, ih)*(1.-rs(k, ih)*rxa(k+1, im, is))-ts(k, &
&               ih)*(-(rsd(k, ih)*rxa(k+1, im, is))-rs(k, ih)*rxad(k+1, &
&               im, is)))/(1.-rs(k, ih)*rxa(k+1, im, is))**2
              denm = ts(k, ih)/(1.-rs(k, ih)*rxa(k+1, im, is))
              rrad(k, im, is) = rrd(k, ih) + (tdd(k, ih)*rra(k+1, im, is&
&               )+td(k, ih)*rrad(k+1, im, is)+(ttd(k, ih)-tdd(k, ih))*&
&               rxa(k+1, im, is)+(tt(k, ih)-td(k, ih))*rxad(k+1, im, is)&
&               )*denm + (td(k, ih)*rra(k+1, im, is)+(tt(k, ih)-td(k, ih&
&               ))*rxa(k+1, im, is))*denmd
              rra(k, im, is) = rr(k, ih) + (td(k, ih)*rra(k+1, im, is)+(&
&               tt(k, ih)-td(k, ih))*rxa(k+1, im, is))*denm
              rxad(k, im, is) = rsd(k, ih) + (tsd(k, ih)*denm+ts(k, ih)*&
&               denmd)*rxa(k+1, im, is) + ts(k, ih)*denm*rxad(k+1, im, &
&               is)
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
              denmd = -((-(rsad(k-1, ih, im)*rxa(k, im, is))-rsa(k-1, ih&
&               , im)*rxad(k, im, is))/(1.-rsa(k-1, ih, im)*rxa(k, im, &
&               is))**2)
              denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
              fdndird = tdad(k-1, ih, im)
              fdndir = tda(k-1, ih, im)
              xx4d = tdad(k-1, ih, im)*rra(k, im, is) + tda(k-1, ih, im)&
&               *rrad(k, im, is)
              xx4 = tda(k-1, ih, im)*rra(k, im, is)
              yyd = ttad(k-1, ih, im) - tdad(k-1, ih, im)
              yy = tta(k-1, ih, im) - tda(k-1, ih, im)
              fdndifd = (xx4d*rsa(k-1, ih, im)+xx4*rsad(k-1, ih, im)+yyd&
&               )*denm + (xx4*rsa(k-1, ih, im)+yy)*denmd
              fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
              fupdifd = (xx4d+yyd*rxa(k, im, is)+yy*rxad(k, im, is))*&
&               denm + (xx4+yy*rxa(k, im, is))*denmd
              fupdif = (xx4+yy*rxa(k, im, is))*denm
              flxdnd = fdndird + fdndifd - fupdifd
              flxdn = fdndir + fdndif - fupdif
!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)
              IF (ih .EQ. 1 .AND. im .EQ. 1 .AND. is .EQ. 1) THEN
                fupc(k) = fupdif
                fclr(k) = flxdn
              END IF
              fupa(k) = fupa(k) + fupdif*ct
              falld(k) = falld(k) + flxdnd*ct + flxdn*ctd
              fall(k) = fall(k) + flxdn*ct
            END DO
! k loop
            fsdir = fsdir + fdndir*ct
            fsdif = fsdif + fdndif*ct
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
        flx_devd(i, k) = flx_devd(i, k) + hk_ir(ib, ik)*falld(k)
        flx_dev(i, k) = flx_dev(i, k) + fall(k)*hk_ir(ib, ik)
        flc_dev(i, k) = flc_dev(i, k) + fclr(k)*hk_ir(ib, ik)
        flxu_dev(i, k) = flxu_dev(i, k) + fupa(k)*hk_ir(ib, ik)
        flcu_dev(i, k) = flcu_dev(i, k) + fupc(k)*hk_ir(ib, ik)
      END DO
!-----compute downward surface fluxes in the ir region
      fdirir_dev(i) = fdirir_dev(i) + fsdir*hk_ir(ib, ik)
      fdifir_dev(i) = fdifir_dev(i) + fsdif*hk_ir(ib, ik)
!-----tabulate surface flux at ir bands
      flx_sfc_band_dev(i, iv) = flx_sfc_band_dev(i, iv) + fall(np+1)*&
&       hk_ir(ib, ik)
    END DO
  END DO
! ik loop
!-----compute pressure-scaled o2 amount following Eq. (3.5) with f=1.
!     unit is (cm-atm)stp. 165.22 = (1000/980)*23.14%*(22400/32)
!     compute flux reduction due to oxygen following Eq. (3.18). 0.0633 is the
!     fraction of insolation contained in the oxygen bands
  dfd(0) = 0.0_8
  df(0) = 0.0
  cnt = 165.22*snt
  so2d(1) = 0.0_8
  so2(1) = scal0*cnt
! LLT increased parameter 145 to 155 to enhance effect
  result1 = SQRT(so2(1))
  dfd(1) = 0.0_8
  df(1) = 0.0633*(1.-EXP(-(0.000155*result1)))
  DO k=1,np
    so2d(k+1) = 0.0_8
    so2(k+1) = so2(k) + scal(k)*cnt
! LLT increased parameter 145 to 155 to enhance effect
    result1 = SQRT(so2(k+1))
    dfd(k+1) = 0.0_8
    df(k+1) = 0.0633*(1.0-EXP(-(0.000155*result1)))
  END DO
!-----for solar heating due to co2 scaling follows Eq(3.5) with f=1.
!     unit is (cm-atm)stp. 789 = (1000/980)*(44/28.97)*(22400/44)
  so2d(1) = 0.0_8
  so2(1) = 789.*co2*scal0
  DO k=1,np
    so2d(k+1) = 0.0_8
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
  du2 = du*du
  dw2 = dw*dw
  x0 = u1 + REAL(nu)*du
  y0 = w1 + REAL(nw)*dw
  x1 = u1 - 0.5*du
  y1 = w1 - 0.5*dw
  dfd = 0.0_8
  DO k=1,np+1
    x3 = LOG10(so2(k)*snt)
    IF (x3 .GT. x0) THEN
      ulog = x0
    ELSE
      ulog = x3
    END IF
    x4d = swhd(k)/(swh(k)*LOG(10.0))
    x4 = LOG10(swh(k)*snt)
    IF (x4 .GT. y0) THEN
      wlog = y0
      wlogd = 0.0_8
    ELSE
      wlogd = x4d
      wlog = x4
    END IF
    ic = INT((ulog-x1)/du + 1.)
    iw = INT((wlog-y1)/dw + 1.)
    IF (ic .LT. 2) ic = 2
    IF (iw .LT. 2) iw = 2
    IF (ic .GT. nu) ic = nu
    IF (iw .GT. nw) iw = nw
    dc = ulog - REAL(ic-2)*du - u1
    ddd = wlogd
    dd = wlog - REAL(iw-2)*dw - w1
    x2d = (cah(ic-1, iw)-cah(ic-1, iw-1))*ddd/dw
    x2 = cah(ic-1, iw-1) + (cah(ic-1, iw)-cah(ic-1, iw-1))/dw*dd
    y2d = x2d
    y2 = x2 + (cah(ic, iw-1)-cah(ic-1, iw-1))/du*dc
    IF (y2 .LT. 0.0) THEN
      y2 = 0.0
      y2d = 0.0_8
    ELSE
      y2 = y2
    END IF
! LLT increase CO2 effect to help reduce cold tropopause bias
    dfd(k) = dfd(k) + 1.5*y2d
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
  du2 = du*du
  dw2 = dw*dw
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
    ic = INT((ulog-x1)/du + 1.)
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
  dftopd = dfd(ntop)
  dftop = df(ntop)
  DO k=1,np+1
    IF (k .GT. ntop) THEN
      xx4d = (flx_devd(i, k)*flx_dev(i, ntop)-flx_dev(i, k)*flx_devd(i, &
&       ntop))/flx_dev(i, ntop)**2
      xx4 = flx_dev(i, k)/flx_dev(i, ntop)
      dfd(k) = dftopd + xx4d*(df(k)-dftop) + xx4*(dfd(k)-dftopd)
      df(k) = dftop + xx4*(df(k)-dftop)
    END IF
  END DO
!-----update the net fluxes
  DO k=1,np+1
    IF (df(k) .GT. flx_dev(i, k) - 1.0e-8) THEN
      dfd(k) = flx_devd(i, k)
      df(k) = flx_dev(i, k) - 1.0e-8
    ELSE
      df(k) = df(k)
    END IF
!           df(k) = 0.0
    flx_devd(i, k) = flx_devd(i, k) - dfd(k)
    flx_dev(i, k) = flx_dev(i, k) - df(k)
    flc_dev(i, k) = flc_dev(i, k) - df(k)
  END DO
!-----update the downward surface fluxes 
!        xx4 = fdirir (i) + fdifir (i) +&
!              fdiruv (i) + fdifuv (i) +&
!              fdirpar(i) + fdifpar(i)
  xx4 = flx_dev(i, np+1) + df(np+1)
  IF (xx4 .GE. 0.) THEN
    abs0 = xx4
  ELSE
    abs0 = -xx4
  END IF
  result10 = EPSILON(1.0)
  IF (abs0 .GT. result10) THEN
    IF (1.0 - df(np+1)/xx4 .GT. 1.) THEN
      x6 = 1.
    ELSE
      x6 = 1.0 - df(np+1)/xx4
    END IF
    IF (x6 .LT. 0.) THEN
      xx4 = 0.
    ELSE
      xx4 = x6
    END IF
  ELSE
    xx4 = 0.0
  END IF
  fdirir_dev(i) = xx4*fdirir_dev(i)
  fdifir_dev(i) = xx4*fdifir_dev(i)
  fdiruv_dev(i) = xx4*fdiruv_dev(i)
  fdifuv_dev(i) = xx4*fdifuv_dev(i)
  fdirpar_dev(i) = xx4*fdirpar_dev(i)
  fdifpar_dev(i) = xx4*fdifpar_dev(i)
  DO ib=1,nband
    flx_sfc_band_dev(i, ib) = xx4*flx_sfc_band_dev(i, ib)
  END DO
END SUBROUTINE SORAD_D

!  Differentiation of deledd in forward (tangent) mode:
!   variations   of useful results: tt1 td1 rr1
!   with respect to varying inputs: g01 tau1 ssc1
!*********************************************************************
SUBROUTINE DELEDD_D(tau1, tau1d, ssc1, ssc1d, g01, g01d, cza1, rr1, rr1d&
& , tt1, tt1d, td1, td1d)
  IMPLICIT NONE
! 8 byte real
  INTEGER, PARAMETER :: real_de=8
!integer,parameter :: REAL_SP = 4 ! 4 byte real
!-----input parameters
  REAL*8, INTENT(IN) :: tau1, ssc1, g01, cza1
  REAL*8, INTENT(IN) :: tau1d, ssc1d, g01d
!-----output parameters
  REAL*8, INTENT(OUT) :: rr1, tt1, td1
  REAL*8, INTENT(OUT) :: rr1d, tt1d, td1d
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
  REAL*8 :: taud, sscd, g0d, rrd, ttd, tdd
  REAL*8 :: zth, ff, xx, taup, sscp, gp, gm1, gm2, gm3, akk, alf1, alf2
  REAL*8 :: ffd, xxd, taupd, sscpd, gpd, gm1d, gm2d, gm3d, akkd, alf1d, &
& alf2d
  REAL*8 :: all, bll, st7, st8, cll, dll, fll, ell, st1, st2, st3, st4
  REAL*8 :: alld, blld, st7d, st8d, clld, dlld, flld, elld, st1d, st2d, &
& st3d, st4d
  INTRINSIC SQRT
  INTRINSIC ABS
  INTRINSIC EXP
  INTRINSIC MAX
  INTRINSIC REAL
  REAL*8 :: arg1
  REAL*8 :: arg1d
  REAL*8 :: abs0
!zth = real(cza1,kind=REAL_DE)
!g0  = real(g01 ,kind=REAL_DE)
!tau = real(tau1,kind=REAL_DE)
!ssc = real(ssc1,kind=REAL_DE)
!zth = dble(cza1)
!g0  = dble(g01)
!tau = dble(tau1)
!ssc = dble(ssc1)
  zth = cza1
  g0d = g01d
  g0 = g01
  taud = tau1d
  tau = tau1
  sscd = ssc1d
  ssc = ssc1
  ffd = g0d*g0 + g0*g0d
  ff = g0*g0
  xxd = -(ffd*ssc) - ff*sscd
  xx = one - ff*ssc
  taupd = taud*xx + tau*xxd
  taup = tau*xx
  sscpd = ((sscd*(one-ff)-ssc*ffd)*xx-ssc*(one-ff)*xxd)/xx**2
  sscp = ssc*(one-ff)/xx
  gpd = (g0d*(one+g0)-g0*g0d)/(one+g0)**2
  gp = g0/(one+g0)
  xxd = three*gpd
  xx = three*gp
  gm1d = fourth*(-(sscpd*(four+xx))-sscp*xxd)
  gm1 = (seven-sscp*(four+xx))*fourth
  gm2d = -(fourth*(sscp*xxd-sscpd*(four-xx)))
  gm2 = -((one-sscp*(four-xx))*fourth)
  arg1d = (gm1d+gm2d)*(gm1-gm2) + (gm1+gm2)*(gm1d-gm2d)
  arg1 = (gm1+gm2)*(gm1-gm2)
  IF (arg1 .EQ. 0.0) THEN
    akkd = 0.0_8
  ELSE
    akkd = arg1d/(2.0*SQRT(arg1))
  END IF
  akk = SQRT(arg1)
  xxd = zth*akkd
  xx = akk*zth
  st7d = -xxd
  st7 = one - xx
  st8d = xxd
  st8 = one + xx
  st3d = st7d*st8 + st7*st8d
  st3 = st7*st8
  IF (st3 .GE. 0.) THEN
    abs0 = st3
  ELSE
    abs0 = -st3
  END IF
  IF (abs0 .LT. thresh) THEN
    zth = zth + 0.0010
    IF (zth .GT. 1.0) zth = zth - 0.0020
    xxd = zth*akkd
    xx = akk*zth
    st7d = -xxd
    st7 = one - xx
    st8d = xxd
    st8 = one + xx
    st3d = st7d*st8 + st7*st8d
    st3 = st7*st8
  END IF
  tdd = -(taupd*EXP(-(taup/zth))/zth)
  td = EXP(-(taup/zth))
  gm3d = -(fourth*zth*three*gpd)
  gm3 = (two-zth*three*gp)*fourth
  xxd = gm1d - gm2d
  xx = gm1 - gm2
  alf1d = gm1d - gm3d*xx - gm3*xxd
  alf1 = gm1 - gm3*xx
  alf2d = gm2d + gm3d*xx + gm3*xxd
  alf2 = gm2 + gm3*xx
  xxd = two*akkd
  xx = akk*two
  alld = (gm3d-zth*alf2d)*xx*td + (gm3-alf2*zth)*(xxd*td+xx*tdd)
  all = (gm3-alf2*zth)*xx*td
  blld = (zth*alf1d-gm3d)*xx + (one-gm3+alf1*zth)*xxd
  bll = (one-gm3+alf1*zth)*xx
  xxd = akkd*gm3 + akk*gm3d
  xx = akk*gm3
  clld = (alf2d+xxd)*st7 + (alf2+xx)*st7d
  cll = (alf2+xx)*st7
  dlld = (alf2d-xxd)*st8 + (alf2-xx)*st8d
  dll = (alf2-xx)*st8
  xxd = akkd*(one-gm3) - akk*gm3d
  xx = akk*(one-gm3)
  flld = (alf1d+xxd)*st8 + (alf1+xx)*st8d
  fll = (alf1+xx)*st8
  elld = (alf1d-xxd)*st7 + (alf1-xx)*st7d
  ell = (alf1-xx)*st7
  st2d = -((akkd*taup+akk*taupd)*EXP(-(akk*taup)))
  st2 = EXP(-(akk*taup))
  st4d = st2d*st2 + st2*st2d
  st4 = st2*st2
  st1d = (sscpd*(akk+gm1+(akk-gm1)*st4)*st3-sscp*((akkd+gm1d+(akkd-gm1d)&
&   *st4+(akk-gm1)*st4d)*st3+(akk+gm1+(akk-gm1)*st4)*st3d))/((akk+gm1+(&
&   akk-gm1)*st4)*st3)**2
  st1 = sscp/((akk+gm1+(akk-gm1)*st4)*st3)
  rrd = (clld-dlld*st4-dll*st4d-alld*st2-all*st2d)*st1 + (cll-dll*st4-&
&   all*st2)*st1d
  rr = (cll-dll*st4-all*st2)*st1
  ttd = -(((flld-elld*st4-ell*st4d)*td+(fll-ell*st4)*tdd-blld*st2-bll*&
&   st2d)*st1+((fll-ell*st4)*td-bll*st2)*st1d)
  tt = -(((fll-ell*st4)*td-bll*st2)*st1)
  IF (rr .LT. zero) THEN
    rr = zero
    rrd = 0.0_8
  ELSE
    rr = rr
  END IF
  IF (tt .LT. zero) THEN
    tt = zero
    ttd = 0.0_8
  ELSE
    tt = tt
  END IF
  ttd = ttd + tdd
  tt = tt + td
!td1 = real(td,kind=REAL_SP)
!rr1 = real(rr,kind=REAL_SP)
!tt1 = real(tt,kind=REAL_SP)
  td1d = tdd
  td1 = REAL(td)
  rr1d = rrd
  rr1 = REAL(rr)
  tt1d = ttd
  tt1 = REAL(tt)
END SUBROUTINE DELEDD_D

!  Differentiation of getvistau1 in forward (tangent) mode:
!   variations   of useful results: asycl taudiff taubeam
!   with respect to varying inputs: hydromets fcld reff
SUBROUTINE GETVISTAU1_D(nlevs, cosz, dp, fcld, fcldd, reff, reffd, &
& hydromets, hydrometsd, ict, icb, taubeam, taubeamd, taudiff, taudiffd&
& , asycl, asycld, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, &
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
  REAL*8, INTENT(IN) :: fcldd(nlevs)
!  Effective radius (microns)
  REAL*8, INTENT(IN) :: reff(nlevs, 4)
  REAL*8, INTENT(IN) :: reffd(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL*8, INTENT(IN) :: hydromets(nlevs, 4)
  REAL*8, INTENT(IN) :: hydrometsd(nlevs, 4)
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
  REAL*8, INTENT(OUT) :: taubeam(nlevs, 4)
  REAL*8, INTENT(OUT) :: taubeamd(nlevs, 4)
!  Optical Depth for Diffuse Radiation
  REAL*8, INTENT(OUT) :: taudiff(nlevs, 4)
  REAL*8, INTENT(OUT) :: taudiffd(nlevs, 4)
!  Cloud Asymmetry Factor
  REAL*8, INTENT(OUT) :: asycl(nlevs)
  REAL*8, INTENT(OUT) :: asycld(nlevs)
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
  REAL*8 :: ftd, fad, xaid, taucd, asycltd
  REAL*8 :: cc(3)
  REAL*8 :: ccd(3)
  REAL*8 :: taucld1, taucld2, taucld3, taucld4
  REAL*8 :: taucld1d, taucld2d, taucld3d, taucld4d
  REAL*8 :: g1, g2, g3, g4
  REAL*8 :: g1d, g2d, g3d, g4d
  REAL*8 :: reff_snow
  REAL*8 :: reff_snowd
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
  taubeam = 0.0
  taudiff = 0.0
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
    ccd = 0.0_8
    DO k=1,ict-1
      IF (cc(1) .LT. fcld(k)) THEN
        ccd(1) = fcldd(k)
        cc(1) = fcld(k)
      ELSE
        cc(1) = cc(1)
      END IF
    END DO
    DO k=ict,icb-1
      IF (cc(2) .LT. fcld(k)) THEN
        ccd(2) = fcldd(k)
        cc(2) = fcld(k)
      ELSE
        cc(2) = cc(2)
      END IF
    END DO
    DO k=icb,nlevs
      IF (cc(3) .LT. fcld(k)) THEN
        ccd(3) = fcldd(k)
        cc(3) = fcld(k)
      ELSE
        cc(3) = cc(3)
      END IF
    END DO
    asycld = 0.0_8
    taudiffd = 0.0_8
    taubeamd = 0.0_8
  ELSE
    asycld = 0.0_8
    taudiffd = 0.0_8
    taubeamd = 0.0_8
    ccd = 0.0_8
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
      taucld1 = 0.
      taucld1d = 0.0_8
    ELSE
      taucld1d = (dp(k)*1.0e3*aib_uv*hydrometsd(k, 1)*reff(k, 1)/&
&       cons_grav-dp(k)*1.0e3*hydromets(k, 1)*aib_uv*reffd(k, 1)/&
&       cons_grav)/reff(k, 1)**2
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*aib_uv/reff(k, 1)
    END IF
    IF (reff(k, 2) .LE. 0.) THEN
      taucld2 = 0.
      taucld2d = 0.0_8
    ELSE
      taucld2d = dp(k)*1.0e3*(hydrometsd(k, 2)*(awb_uv(1)+awb_uv(2)/reff&
&       (k, 2))-hydromets(k, 2)*awb_uv(2)*reffd(k, 2)/reff(k, 2)**2)/&
&       cons_grav
      taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_uv(1)+awb_uv(&
&       2)/reff(k, 2))
    END IF
    taucld3d = dp(k)*1.0e3*arb_uv(1)*hydrometsd(k, 3)/cons_grav
    taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_uv(1)
    IF (reff(k, 4) .GT. 112.0) THEN
      reff_snow = 112.0
      reff_snowd = 0.0_8
    ELSE
      reff_snowd = reffd(k, 4)
      reff_snow = reff(k, 4)
    END IF
    IF (reff_snow .LE. 0.) THEN
      taucld4 = 0.
      taucld4d = 0.0_8
    ELSE
      taucld4d = (dp(k)*1.0e3*aib_uv*hydrometsd(k, 4)*reff_snow/&
&       cons_grav-dp(k)*1.0e3*hydromets(k, 4)*aib_uv*reff_snowd/&
&       cons_grav)/reff_snow**2
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*aib_uv/reff_snow
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
        kk = 1
      ELSE IF (k .GE. ict .AND. k .LT. icb) THEN
        kk = 2
      ELSE
        kk = 3
      END IF
      taucd = taucld1d + taucld2d + taucld3d + taucld4d
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
!-----normalize cloud cover following Eq. (7.8)
        fad = (fcldd(k)*cc(kk)-fcld(k)*ccd(kk))/cc(kk)**2
        fa = fcld(k)/cc(kk)
        IF (tauc .GT. 32.) THEN
          tauc = 32.
          taucd = 0.0_8
        ELSE
          tauc = tauc
        END IF
        fm = cosz/dm
        ftd = taucd/(tauc*LOG(10.0))/dt
        ft = (LOG10(tauc)-t1)/dt
        fad = fad/da
        fa = fa/da
        im = INT(fm + 1.5)
        it = INT(ft + 1.5)
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
        xai = (-(caib(im-1, it, ia)*(1.-fm))+caib(im+1, it, ia)*(1.+fm))&
&         *fm*.5 + caib(im, it, ia)*(1.-fm*fm)
        xaid = .5*((caib(im, it-1, ia)*ftd+caib(im, it+1, ia)*ftd)*ft+(-&
&         (caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(1.+ft))*ftd) &
&         + caib(im, it, ia)*(-(ftd*ft)-ft*ftd)
        xai = xai + (-(caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(&
&         1.+ft))*ft*.5 + caib(im, it, ia)*(1.-ft*ft)
        xaid = xaid + .5*((caib(im, it, ia-1)*fad+caib(im, it, ia+1)*fad&
&         )*fa+(-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(1.+fa)&
&         )*fad) + caib(im, it, ia)*(-(fad*fa)-fa*fad)
        xai = xai + (-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(&
&         1.+fa))*fa*.5 + caib(im, it, ia)*(1.-fa*fa)
        xai = xai - 2.*caib(im, it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          xaid = 0.0_8
        ELSE
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          xaid = 0.0_8
        ELSE
          xai = xai
        END IF
        taubeamd(k, 1) = taucld1d*xai + taucld1*xaid
        taubeam(k, 1) = taucld1*xai
        taubeamd(k, 2) = taucld2d*xai + taucld2*xaid
        taubeam(k, 2) = taucld2*xai
        taubeamd(k, 3) = taucld3d*xai + taucld3*xaid
        taubeam(k, 3) = taucld3*xai
        taubeamd(k, 4) = taucld4d*xai + taucld4*xaid
        taubeam(k, 4) = taucld4*xai
!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
        xaid = .5*((caif(it-1, ia)*ftd+caif(it+1, ia)*ftd)*ft+(-(caif(it&
&         -1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ftd) + caif(it, ia)*(&
&         -(ftd*ft)-ft*ftd)
        xai = (-(caif(it-1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ft*.5 +&
&         caif(it, ia)*(1.-ft*ft)
        xaid = xaid + .5*((caif(it, ia-1)*fad+caif(it, ia+1)*fad)*fa+(-(&
&         caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*fad) + caif(it&
&         , ia)*(-(fad*fa)-fa*fad)
        xai = xai + (-(caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*&
&         fa*.5 + caif(it, ia)*(1.-fa*fa)
        xai = xai - caif(it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          xaid = 0.0_8
        ELSE
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          xaid = 0.0_8
        ELSE
          xai = xai
        END IF
        taudiffd(k, 1) = taucld1d*xai + taucld1*xaid
        taudiff(k, 1) = taucld1*xai
        taudiffd(k, 2) = taucld2d*xai + taucld2*xaid
        taudiff(k, 2) = taucld2*xai
        taudiffd(k, 3) = taucld3d*xai + taucld3*xaid
        taudiff(k, 3) = taucld3*xai
        taudiffd(k, 4) = taucld4d*xai + taucld4*xaid
        taudiff(k, 4) = taucld4*xai
      END IF
    ELSE
! Overlap calculation scaling not needed
      taubeamd(k, 1) = taucld1d
      taubeam(k, 1) = taucld1
      taubeamd(k, 2) = taucld2d
      taubeam(k, 2) = taucld2
      taubeamd(k, 3) = taucld3d
      taubeam(k, 3) = taucld3
      taubeamd(k, 4) = taucld4d
      taubeam(k, 4) = taucld4
      taudiffd(k, 1) = taucld1d
      taudiff(k, 1) = taucld1
      taudiffd(k, 2) = taucld2d
      taudiff(k, 2) = taucld2
      taudiffd(k, 3) = taucld3d
      taudiff(k, 3) = taucld3
      taudiffd(k, 4) = taucld4d
      taudiff(k, 4) = taucld4
    END IF
!-----cloud asymmetry factor for a mixture of liquid and ice particles.
!     unit of reff is micrometers. Eqs. (4.8) and (6.4)
    asyclt = 1.0
    taucd = taucld1d + taucld2d + taucld3d + taucld4d
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      g1d = (aig_uv(3)*reffd(k, 1)*reff(k, 1)+(aig_uv(2)+aig_uv(3)*reff(&
&       k, 1))*reffd(k, 1))*taucld1 + (aig_uv(1)+(aig_uv(2)+aig_uv(3)*&
&       reff(k, 1))*reff(k, 1))*taucld1d
      g1 = (aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff(k, 1))*reff(k, 1))*&
&       taucld1
      g2d = (awg_uv(3)*reffd(k, 2)*reff(k, 2)+(awg_uv(2)+awg_uv(3)*reff(&
&       k, 2))*reffd(k, 2))*taucld2 + (awg_uv(1)+(awg_uv(2)+awg_uv(3)*&
&       reff(k, 2))*reff(k, 2))*taucld2d
      g2 = (awg_uv(1)+(awg_uv(2)+awg_uv(3)*reff(k, 2))*reff(k, 2))*&
&       taucld2
      g3d = arg_uv(1)*taucld3d
      g3 = arg_uv(1)*taucld3
      g4d = (aig_uv(3)*reff_snowd*reff_snow+(aig_uv(2)+aig_uv(3)*&
&       reff_snow)*reff_snowd)*taucld4 + (aig_uv(1)+(aig_uv(2)+aig_uv(3)&
&       *reff_snow)*reff_snow)*taucld4d
      g4 = (aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff_snow)*reff_snow)*taucld4
      asycltd = ((g1d+g2d+g3d+g4d)*tauc-(g1+g2+g3+g4)*taucd)/tauc**2
      asyclt = (g1+g2+g3+g4)/tauc
    ELSE
      asycltd = 0.0_8
    END IF
    asycld(k) = asycltd
    asycl(k) = asyclt
  END DO
  RETURN
END SUBROUTINE GETVISTAU1_D

!  Differentiation of getnirtau1 in forward (tangent) mode:
!   variations   of useful results: asycl taudiff ssacl taubeam
!   with respect to varying inputs: hydromets asycl fcld ssacl
!                reff
SUBROUTINE GETNIRTAU1_D(ib, nlevs, cosz, dp, fcld, fcldd, reff, reffd, &
& hydromets, hydrometsd, ict, icb, taubeam, taubeamd, taudiff, taudiffd&
& , asycl, asycld, ssacl, ssacld, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv&
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
  REAL*8, INTENT(IN) :: fcldd(nlevs)
!  Effective radius (microns)
  REAL*8, INTENT(IN) :: reff(nlevs, 4)
  REAL*8, INTENT(IN) :: reffd(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL*8, INTENT(IN) :: hydromets(nlevs, 4)
  REAL*8, INTENT(IN) :: hydrometsd(nlevs, 4)
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
  REAL*8, INTENT(OUT) :: taubeam(nlevs, 4)
  REAL*8, INTENT(OUT) :: taubeamd(nlevs, 4)
!  Optical depth for diffuse radiation
  REAL*8, INTENT(OUT) :: taudiff(nlevs, 4)
  REAL*8, INTENT(OUT) :: taudiffd(nlevs, 4)
!  Cloud single scattering albedo
  REAL*8, INTENT(OUT) :: ssacl(nlevs)
  REAL*8, INTENT(OUT) :: ssacld(nlevs)
!  Cloud asymmetry factor
  REAL*8, INTENT(OUT) :: asycl(nlevs)
  REAL*8, INTENT(OUT) :: asycld(nlevs)
  INTEGER :: k, in, im, it, ia, kk
  REAL*8 :: fm, ft, fa, xai, tauc, asyclt, ssaclt
  REAL*8 :: ftd, fad, xaid, taucd, asycltd, ssacltd
  REAL*8 :: cc(3)
  REAL*8 :: ccd(3)
  REAL*8 :: taucld1, taucld2, taucld3, taucld4
  REAL*8 :: taucld1d, taucld2d, taucld3d, taucld4d
  REAL*8 :: g1, g2, g3, g4
  REAL*8 :: g1d, g2d, g3d, g4d
  REAL*8 :: w1, w2, w3, w4
  REAL*8 :: w1d, w2d, w3d, w4d
  REAL*8 :: reff_snow
  REAL*8 :: reff_snowd
  INTEGER, PARAMETER :: nm=11, nt=9, na=11
  REAL*8, PARAMETER :: dm=0.1, dt=0.30103, da=0.1, t1=-0.9031
  INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC REAL
  taubeam = 0.0
  taudiff = 0.0
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
    ccd = 0.0_8
    DO k=1,ict-1
      IF (cc(1) .LT. fcld(k)) THEN
        ccd(1) = fcldd(k)
        cc(1) = fcld(k)
      ELSE
        cc(1) = cc(1)
      END IF
    END DO
    DO k=ict,icb-1
      IF (cc(2) .LT. fcld(k)) THEN
        ccd(2) = fcldd(k)
        cc(2) = fcld(k)
      ELSE
        cc(2) = cc(2)
      END IF
    END DO
    DO k=icb,nlevs
      IF (cc(3) .LT. fcld(k)) THEN
        ccd(3) = fcldd(k)
        cc(3) = fcld(k)
      ELSE
        cc(3) = cc(3)
      END IF
    END DO
    taudiffd = 0.0_8
    taubeamd = 0.0_8
  ELSE
    taudiffd = 0.0_8
    taubeamd = 0.0_8
    ccd = 0.0_8
  END IF
!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.) THEN
      taucld1 = 0.
      taucld1d = 0.0_8
    ELSE
      taucld1d = (dp(k)*1.0e3*aib_nir*hydrometsd(k, 1)*reff(k, 1)/&
&       cons_grav-dp(k)*1.0e3*hydromets(k, 1)*aib_nir*reffd(k, 1)/&
&       cons_grav)/reff(k, 1)**2
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*aib_nir/reff(k, 1)
    END IF
    IF (reff(k, 2) .LE. 0.) THEN
      taucld2 = 0.
      taucld2d = 0.0_8
    ELSE
      taucld2d = dp(k)*1.0e3*(hydrometsd(k, 2)*(awb_nir(ib, 1)+awb_nir(&
&       ib, 2)/reff(k, 2))-hydromets(k, 2)*awb_nir(ib, 2)*reffd(k, 2)/&
&       reff(k, 2)**2)/cons_grav
      taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_nir(ib, 1)+&
&       awb_nir(ib, 2)/reff(k, 2))
    END IF
    taucld3d = dp(k)*1.0e3*arb_nir(ib, 1)*hydrometsd(k, 3)/cons_grav
    taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_nir(ib, 1)
    IF (reff(k, 4) .GT. 112.0) THEN
      reff_snow = 112.0
      reff_snowd = 0.0_8
    ELSE
      reff_snowd = reffd(k, 4)
      reff_snow = reff(k, 4)
    END IF
    IF (reff_snow .LE. 0.) THEN
      taucld4 = 0.
      taucld4d = 0.0_8
    ELSE
      taucld4d = (dp(k)*1.0e3*aib_nir*hydrometsd(k, 4)*reff_snow/&
&       cons_grav-dp(k)*1.0e3*hydromets(k, 4)*aib_nir*reff_snowd/&
&       cons_grav)/reff_snow**2
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*aib_nir/reff_snow
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
        kk = 1
      ELSE IF (k .GE. ict .AND. k .LT. icb) THEN
        kk = 2
      ELSE
        kk = 3
      END IF
      taucd = taucld1d + taucld2d + taucld3d + taucld4d
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
!-----normalize cloud cover following Eq. (7.8)
        IF (cc(kk) .NE. 0.0) THEN
          fad = (fcldd(k)*cc(kk)-fcld(k)*ccd(kk))/cc(kk)**2
          fa = fcld(k)/cc(kk)
        ELSE
          fa = 0.0
          fad = 0.0_8
        END IF
        IF (tauc .GT. 32.) THEN
          tauc = 32.
          taucd = 0.0_8
        ELSE
          tauc = tauc
        END IF
        fm = cosz/dm
        ftd = taucd/(tauc*LOG(10.0))/dt
        ft = (LOG10(tauc)-t1)/dt
        fad = fad/da
        fa = fa/da
        im = INT(fm + 1.5)
        it = INT(ft + 1.5)
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
        xai = (-(caib(im-1, it, ia)*(1.-fm))+caib(im+1, it, ia)*(1.+fm))&
&         *fm*.5 + caib(im, it, ia)*(1.-fm*fm)
        xaid = .5*((caib(im, it-1, ia)*ftd+caib(im, it+1, ia)*ftd)*ft+(-&
&         (caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(1.+ft))*ftd) &
&         + caib(im, it, ia)*(-(ftd*ft)-ft*ftd)
        xai = xai + (-(caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(&
&         1.+ft))*ft*.5 + caib(im, it, ia)*(1.-ft*ft)
        xaid = xaid + .5*((caib(im, it, ia-1)*fad+caib(im, it, ia+1)*fad&
&         )*fa+(-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(1.+fa)&
&         )*fad) + caib(im, it, ia)*(-(fad*fa)-fa*fad)
        xai = xai + (-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(&
&         1.+fa))*fa*.5 + caib(im, it, ia)*(1.-fa*fa)
        xai = xai - 2.*caib(im, it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          xaid = 0.0_8
        ELSE
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          xaid = 0.0_8
        ELSE
          xai = xai
        END IF
        taubeamd(k, 1) = taucld1d*xai + taucld1*xaid
        taubeam(k, 1) = taucld1*xai
        taubeamd(k, 2) = taucld2d*xai + taucld2*xaid
        taubeam(k, 2) = taucld2*xai
        taubeamd(k, 3) = taucld3d*xai + taucld3*xaid
        taubeam(k, 3) = taucld3*xai
        taubeamd(k, 4) = taucld4d*xai + taucld4*xaid
        taubeam(k, 4) = taucld4*xai
!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
        xaid = .5*((caif(it-1, ia)*ftd+caif(it+1, ia)*ftd)*ft+(-(caif(it&
&         -1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ftd) + caif(it, ia)*(&
&         -(ftd*ft)-ft*ftd)
        xai = (-(caif(it-1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ft*.5 +&
&         caif(it, ia)*(1.-ft*ft)
        xaid = xaid + .5*((caif(it, ia-1)*fad+caif(it, ia+1)*fad)*fa+(-(&
&         caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*fad) + caif(it&
&         , ia)*(-(fad*fa)-fa*fad)
        xai = xai + (-(caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*&
&         fa*.5 + caif(it, ia)*(1.-fa*fa)
        xai = xai - caif(it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          xaid = 0.0_8
        ELSE
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          xaid = 0.0_8
        ELSE
          xai = xai
        END IF
        taudiffd(k, 1) = taucld1d*xai + taucld1*xaid
        taudiff(k, 1) = taucld1*xai
        taudiffd(k, 2) = taucld2d*xai + taucld2*xaid
        taudiff(k, 2) = taucld2*xai
        taudiffd(k, 3) = taucld3d*xai + taucld3*xaid
        taudiff(k, 3) = taucld3*xai
        taudiffd(k, 4) = taucld4d*xai + taucld4*xaid
        taudiff(k, 4) = taucld4*xai
      END IF
    ELSE
! Overlap calculation scaling not needed
      taubeamd(k, 1) = taucld1d
      taubeam(k, 1) = taucld1
      taubeamd(k, 2) = taucld2d
      taubeam(k, 2) = taucld2
      taubeamd(k, 3) = taucld3d
      taubeam(k, 3) = taucld3
      taubeamd(k, 4) = taucld4d
      taubeam(k, 4) = taucld4
      taudiffd(k, 1) = taucld1d
      taudiff(k, 1) = taucld1
      taudiffd(k, 2) = taucld2d
      taudiff(k, 2) = taucld2
      taudiffd(k, 3) = taucld3d
      taudiff(k, 3) = taucld3
      taudiffd(k, 4) = taucld4d
      taudiff(k, 4) = taucld4
    END IF
!-----compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
    ssaclt = 0.99999
    asyclt = 1.0
    taucd = taucld1d + taucld2d + taucld3d + taucld4d
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      w1d = (-(aia_nir(ib, 3)*reffd(k, 1)*reff(k, 1))-(aia_nir(ib, 2)+&
&       aia_nir(ib, 3)*reff(k, 1))*reffd(k, 1))*taucld1 + (1.-(aia_nir(&
&       ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff(k, 1))*reff(k, 1)))*&
&       taucld1d
      w1 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff(k, 1)&
&       )*reff(k, 1)))*taucld1
      w2d = (-(awa_nir(ib, 3)*reffd(k, 2)*reff(k, 2))-(awa_nir(ib, 2)+&
&       awa_nir(ib, 3)*reff(k, 2))*reffd(k, 2))*taucld2 + (1.-(awa_nir(&
&       ib, 1)+(awa_nir(ib, 2)+awa_nir(ib, 3)*reff(k, 2))*reff(k, 2)))*&
&       taucld2d
      w2 = (1.-(awa_nir(ib, 1)+(awa_nir(ib, 2)+awa_nir(ib, 3)*reff(k, 2)&
&       )*reff(k, 2)))*taucld2
      w3d = (1.-ara_nir(ib, 1))*taucld3d
      w3 = (1.-ara_nir(ib, 1))*taucld3
      w4d = (-(aia_nir(ib, 3)*reff_snowd*reff_snow)-(aia_nir(ib, 2)+&
&       aia_nir(ib, 3)*reff_snow)*reff_snowd)*taucld4 + (1.-(aia_nir(ib&
&       , 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff_snow)*reff_snow))*&
&       taucld4d
      w4 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff_snow)&
&       *reff_snow))*taucld4
      ssacltd = ((w1d+w2d+w3d+w4d)*tauc-(w1+w2+w3+w4)*taucd)/tauc**2
      ssaclt = (w1+w2+w3+w4)/tauc
      g1d = (aig_nir(ib, 3)*reffd(k, 1)*reff(k, 1)+(aig_nir(ib, 2)+&
&       aig_nir(ib, 3)*reff(k, 1))*reffd(k, 1))*w1 + (aig_nir(ib, 1)+(&
&       aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 1))*reff(k, 1))*w1d
      g1 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 1))*&
&       reff(k, 1))*w1
      g2d = (awg_nir(ib, 3)*reffd(k, 2)*reff(k, 2)+(awg_nir(ib, 2)+&
&       awg_nir(ib, 3)*reff(k, 2))*reffd(k, 2))*w2 + (awg_nir(ib, 1)+(&
&       awg_nir(ib, 2)+awg_nir(ib, 3)*reff(k, 2))*reff(k, 2))*w2d
      g2 = (awg_nir(ib, 1)+(awg_nir(ib, 2)+awg_nir(ib, 3)*reff(k, 2))*&
&       reff(k, 2))*w2
      g3d = arg_nir(ib, 1)*w3d
      g3 = arg_nir(ib, 1)*w3
      g4d = (aig_nir(ib, 3)*reffd(k, 4)*reff(k, 4)+(aig_nir(ib, 2)+&
&       aig_nir(ib, 3)*reff(k, 4))*reffd(k, 4))*w4 + (aig_nir(ib, 1)+(&
&       aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 4))*reff(k, 4))*w4d
      g4 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 4))*&
&       reff(k, 4))*w4
      IF (w1 + w2 + w3 + w4 .NE. 0.0) THEN
        asycltd = ((g1d+g2d+g3d+g4d)*(w1+w2+w3+w4)-(g1+g2+g3+g4)*(w1d+&
&         w2d+w3d+w4d))/(w1+w2+w3+w4)**2
        asyclt = (g1+g2+g3+g4)/(w1+w2+w3+w4)
      ELSE
        asycltd = 0.0_8
      END IF
    ELSE
      ssacltd = 0.0_8
      asycltd = 0.0_8
    END IF
    ssacld(k) = ssacltd
    ssacl(k) = ssaclt
    asycld(k) = asycltd
    asycl(k) = asyclt
  END DO
  RETURN
END SUBROUTINE GETNIRTAU1_D

end module SORAD_TL
