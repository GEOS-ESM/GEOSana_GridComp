module IRRAD_TL

use IRRADMOD, only: sfcflux, mkicx

IMPLICIT NONE

PRIVATE
PUBLIC :: irrad_d

contains

SUBROUTINE IRRAD_D(np, ple_dev, ta_dev, ta_devd, wa_dev, wa_devd, oa_dev&
& , oa_devd, tb_dev, tb_devd, co2, trace, n2o_dev, ch4_dev, cfc11_dev, &
& cfc12_dev, cfc22_dev, cwc_dev, cwc_devd, fcld_dev, fcld_devd, ict, icb&
& , reff_dev, reff_devd, ns, fs_dev, tg_dev, eg_dev, tv_dev, ev_dev, &
& rv_dev, na, nb, taua_dev, taua_devd, ssaa_dev, ssaa_devd, asya_dev, &
& asya_devd, flxu_dev, flxu_devd, flxd_dev, flxd_devd, dfdts_dev, &
& dfdts_devd, aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, awg_ir, xkw, xke, &
& mw, aw, bw, pm, fkw, gkw, cb, dcb, w11, w12, w13, p11, p12, p13, dwe, &
& dpe, c1, c2, c3, oo1, oo2, oo3, h11, h12, h13, h21, h22, h23, h81, h82&
& , h83)
  IMPLICIT NONE
!Radiation constants, these need to be inputs for the autodiff tool
  REAL*8, INTENT(IN) :: aib_ir(3, 10), awb_ir(4, 10), aiw_ir(4, 10)
  REAL*8, INTENT(IN) :: aww_ir(4, 10), aig_ir(4, 10), awg_ir(4, 10)
  INTEGER, INTENT(IN) :: mw(9)
  REAL*8, INTENT(IN) :: xkw(9), xke(9), aw(9), bw(9), pm(9)
  REAL*8, INTENT(IN) :: fkw(6, 9), gkw(6, 3), cb(6, 10), dcb(5, 10)
  REAL*8, INTENT(IN) :: w11, w12, w13, p11, p12
  REAL*8, INTENT(IN) :: p13, dwe, dpe
  REAL*8, INTENT(IN) :: c1(26, 30), c2(26, 30), c3(26, 30)
  REAL*8, INTENT(IN) :: oo1(26, 21), oo2(26, 21), oo3(26, 21)
  REAL*8, INTENT(IN) :: h11(26, 31), h12(26, 31), h13(26, 31)
  REAL*8, INTENT(IN) :: h21(26, 31), h22(26, 31), h23(26, 31)
  REAL*8, INTENT(IN) :: h81(26, 31), h82(26, 31), h83(26, 31)
!----- INPUTS -----
  INTEGER, INTENT(IN) :: np, ict, icb, ns, na, nb
  LOGICAL, INTENT(IN) :: trace
  REAL*8, INTENT(IN) :: co2
!Rank 2 inputs
  REAL*8, INTENT(IN) :: tb_dev
  REAL*8, INTENT(IN) :: tb_devd
!Rank 3 (Prognostic variables and tracers)
  REAL*8, DIMENSION(np), INTENT(IN) :: ta_dev, wa_dev, oa_dev, fcld_dev
  REAL*8, DIMENSION(np), INTENT(IN) :: ta_devd, wa_devd, oa_devd, &
& fcld_devd
  REAL*8, DIMENSION(np), INTENT(IN) :: n2o_dev, ch4_dev, cfc11_dev, &
& cfc12_dev, cfc22_dev
  REAL*8, DIMENSION(np + 1), INTENT(IN) :: ple_dev
!Rank 3 (surface types)
  REAL*8, DIMENSION(ns), INTENT(IN) :: fs_dev, tg_dev, tv_dev
  REAL*8, DIMENSION(ns, 10), INTENT(IN) :: eg_dev, ev_dev, rv_dev
!Rank 3 (diagnostic cloud parts)
  REAL*8, DIMENSION(np, 4), INTENT(IN) :: cwc_dev, reff_dev
  REAL*8, DIMENSION(np, 4), INTENT(IN) :: cwc_devd, reff_devd
!Rank 3 (aerosols)
  REAL*8, DIMENSION(np, nb), INTENT(INOUT) :: taua_dev, ssaa_dev, &
& asya_dev
  REAL*8, DIMENSION(np, nb), INTENT(INOUT) :: taua_devd, ssaa_devd, &
& asya_devd
!----- OUPUTS -----
  REAL*8, DIMENSION(np + 1), INTENT(OUT) :: flxu_dev
  REAL*8, DIMENSION(np+1), INTENT(OUT) :: flxu_devd
  REAL*8, DIMENSION(np + 1), INTENT(OUT) :: flxd_dev
  REAL*8, DIMENSION(np+1), INTENT(OUT) :: flxd_devd
  REAL*8, DIMENSION(np + 1), INTENT(OUT) :: dfdts_dev
  REAL*8, DIMENSION(np+1), INTENT(OUT) :: dfdts_devd
!----- LOCALS -----
  REAL*8, PARAMETER :: cons_grav=9.80665
  INTEGER, PARAMETER :: nx1=26
  INTEGER, PARAMETER :: no1=21
  INTEGER, PARAMETER :: nc1=30
  INTEGER, PARAMETER :: nh1=31
!Temporary arrays
  REAL*8 :: pa(0:np), dt(0:np)
  REAL*8 :: pad(0:np), dtd(0:np)
  REAL*8 :: x1, x2, x3
  REAL*8 :: x1d, x2d, x3d
  REAL*8 :: dh2o(0:np), dcont(0:np), dco2(0:np), do3(0:np)
  REAL*8 :: dh2od(0:np), dcontd(0:np), dco2d(0:np), do3d(0:np)
  REAL*8 :: dn2o(0:np), dch4(0:np)
  REAL*8 :: dn2od(0:np), dch4d(0:np)
  REAL*8 :: df11(0:np), df12(0:np), df22(0:np)
  REAL*8 :: df11d(0:np), df12d(0:np), df22d(0:np)
  REAL*8 :: th2o(6), tcon(3), tco2(6)
  REAL*8 :: th2od(6), tcond(3), tco2d(6)
  REAL*8 :: tn2o(4), tch4(4), tcom(6)
  REAL*8 :: tn2od(4), tch4d(4), tcomd(6)
  REAL*8 :: tf11, tf12, tf22
  REAL*8 :: tf11d, tf12d, tf22d
  REAL*8 :: blayer(0:np+1), blevel(0:np+1)
  REAL*8 :: blayerd(0:np+1), bleveld(0:np+1)
  REAL*8 :: bd(0:np+1), bu(0:np+1)
  REAL*8 :: bdd(0:np+1), bud(0:np+1)
  REAL*8 :: bs, dbs, rflxs
  REAL*8 :: dp(0:np)
  REAL*8 :: dpd(0:np)
  REAL*8 :: trant, tranal
  REAL*8 :: trantd, tranald
  REAL*8 :: transfc(0:np+1)
  REAL*8 :: transfcd(0:np+1)
  REAL*8 :: flxu(0:np+1), flxd(0:np+1)
  REAL*8 :: flxud(0:np+1), flxdd(0:np+1)
  REAL*8 :: taerlyr(0:np)
  REAL*8 :: taerlyrd(0:np)
!OVERCAST
  INTEGER :: ncld(3)
  INTEGER :: icx(0:np)
!OVERCAST
  INTEGER :: idx, rc
  INTEGER :: k, l, ip, iw, ibn, ik, iq, isb, k1, k2, ne
  REAL*8 :: enn(0:np)
  REAL*8 :: ennd(0:np)
  REAL*8 :: cldhi, cldmd, cldlw, tcldlyr(0:np), fclr
  REAL*8 :: cldhid, cldmdd, cldlwd, tcldlyrd(0:np), fclrd
  REAL*8 :: x, xx, yy, p1, a1, b1, fk1, a2, b2, fk2
  REAL*8 :: xxd, yyd
  REAL*8 :: w1, ff
  REAL*8 :: ffd
  LOGICAL :: oznbnd, co2bnd, h2otable, conbnd, n2obnd
  LOGICAL :: ch4bnd, combnd, f11bnd, f12bnd, f22bnd, b10bnd
  LOGICAL :: do_aerosol
!Temp arrays and variables for consolidation of tables
  INTEGER, PARAMETER :: max_num_tables=17
  REAL*8 :: exptbl(0:np, max_num_tables)
  REAL*8 :: exptbld(0:np, max_num_tables)
  TYPE BAND_TABLE
      INTEGER :: start
      INTEGER :: end
  END TYPE BAND_TABLE
  TYPE(BAND_TABLE) :: h2oexp
  TYPE(BAND_TABLE) :: conexp
  TYPE(BAND_TABLE) :: co2exp
  TYPE(BAND_TABLE) :: n2oexp
  TYPE(BAND_TABLE) :: ch4exp
  TYPE(BAND_TABLE) :: comexp
  TYPE(BAND_TABLE) :: f11exp
  TYPE(BAND_TABLE) :: f12exp
  TYPE(BAND_TABLE) :: f22exp
!Variables for new getirtau routine
  REAL*8 :: dp_pa(np)
  REAL*8 :: dp_pad(np)
  REAL*8 :: fcld_col(np)
  REAL*8 :: fcld_cold(np)
  REAL*8 :: reff_col(np, 4)
  REAL*8 :: reff_cold(np, 4)
  REAL*8 :: cwc_col(np, 4)
  REAL*8 :: cwc_cold(np, 4)
  REAL*8 :: h2oexp_tmp(0:np, 5), conexp_tmp(0:np), co2exp_tmp(0:np, 6), &
& n2oexp_tmp(0:np, 2)
  REAL*8 :: h2oexp_tmpd(0:np, 5), co2exp_tmpd(0:np, 6), n2oexp_tmpd(0:np&
& , 2)
  INTRINSIC MAX
  INTRINSIC EXP
  INTRINSIC MIN
  INTRINSIC LOG
  dh2od = 0.0_8
  dcontd = 0.0_8
  dtd = 0.0_8
  reff_cold = 0.0_8
  fcld_cold = 0.0_8
  cwc_cold = 0.0_8
  do3d = 0.0_8
!BEGIN CALCULATIONS ...
!compute layer pressure (pa) and layer temperature minus 250K (dt) 
  DO k=1,np
    pad(k) = 0.0_8
    pa(k) = 0.5*(ple_dev(k+1)+ple_dev(k))*0.01
    dpd(k) = 0.0_8
    dp(k) = (ple_dev(k+1)-ple_dev(k))*0.01
! dp in Pascals for getirtau
    dp_pad(k) = 0.0_8
    dp_pa(k) = ple_dev(k+1) - ple_dev(k)
    dtd(k) = ta_devd(k)
    dt(k) = ta_dev(k) - 250.0
!compute layer absorber amount
!dh2o : water vapor amount (g/cm**2)
!dcont: scaled water vapor amount for continuum absorption
!       (g/cm**2)
!dco2 : co2 amount (cm-atm)stp
!do3  : o3 amount (cm-atm)stp
!dn2o : n2o amount (cm-atm)stp
!dch4 : ch4 amount (cm-atm)stp
!df11 : cfc11 amount (cm-atm)stp
!df12 : cfc12 amount (cm-atm)stp
!df22 : cfc22 amount (cm-atm)stp
!the factor 1.02 is equal to 1000/980
!factors 789 and 476 are for unit conversion
!the factor 0.001618 is equal to 1.02/(.622*1013.25) 
!the factor 6.081 is equal to 1800/296
    dh2od(k) = 1.02*dp(k)*wa_devd(k)
    dh2o(k) = 1.02*wa_dev(k)*dp(k)
    do3d(k) = 476.*dp(k)*oa_devd(k)
    do3(k) = 476.*oa_dev(k)*dp(k)
    dco2d(k) = 0.0_8
    dco2(k) = 789.*co2*dp(k)
    dch4d(k) = 0.0_8
    dch4(k) = 789.*ch4_dev(k)*dp(k)
    dn2od(k) = 0.0_8
    dn2o(k) = 789.*n2o_dev(k)*dp(k)
    df11d(k) = 0.0_8
    df11(k) = 789.*cfc11_dev(k)*dp(k)
    df12d(k) = 0.0_8
    df12(k) = 789.*cfc12_dev(k)*dp(k)
    df22d(k) = 0.0_8
    df22(k) = 789.*cfc22_dev(k)*dp(k)
    IF (dh2o(k) .LT. 1.e-10) THEN
      dh2od(k) = 0.0_8
      dh2o(k) = 1.e-10
    ELSE
      dh2o(k) = dh2o(k)
    END IF
    IF (do3(k) .LT. 1.e-6) THEN
      do3d(k) = 0.0_8
      do3(k) = 1.e-6
    ELSE
      do3(k) = do3(k)
    END IF
    IF (dco2(k) .LT. 1.e-4) THEN
      dco2d(k) = 0.0_8
      dco2(k) = 1.e-4
    ELSE
      dco2d(k) = 0.0_8
      dco2(k) = dco2(k)
    END IF
!Compute scaled water vapor amount for h2o continuum absorption
!following eq. (4.21).
    xxd = pa(k)*0.001618*dp(k)*(wa_devd(k)*wa_dev(k)+wa_dev(k)*wa_devd(k&
&     ))
    xx = pa(k)*0.001618*wa_dev(k)*wa_dev(k)*dp(k)
    dcontd(k) = xxd*EXP(1800./ta_dev(k)-6.081) - xx*1800.*ta_devd(k)*EXP&
&     (1800./ta_dev(k)-6.081)/ta_dev(k)**2
    dcont(k) = xx*EXP(1800./ta_dev(k)-6.081)
!Fill the reff, cwc, and fcld for the column   
    fcld_cold(k) = fcld_devd(k)
    fcld_col(k) = fcld_dev(k)
    DO l=1,4
      reff_cold(k, l) = reff_devd(k, l)
      reff_col(k, l) = reff_dev(k, l)
      cwc_cold(k, l) = cwc_devd(k, l)
      cwc_col(k, l) = cwc_dev(k, l)
    END DO
  END DO
  IF (ple_dev(1)*0.01 .LT. 0.005) THEN
    dpd(0) = 0.0_8
    dp(0) = 0.005
  ELSE
    dpd(0) = 0.0_8
    dp(0) = ple_dev(1)*0.01
  END IF
  pad(0) = 0.0_8
  pa(0) = 0.5*dp(0)
  dtd(0) = ta_devd(1)
  dt(0) = ta_dev(1) - 250.0
  dh2od(0) = 1.02*dp(0)*wa_devd(1)
  dh2o(0) = 1.02*wa_dev(1)*dp(0)
  do3d(0) = 476.*dp(0)*oa_devd(1)
  do3(0) = 476.*oa_dev(1)*dp(0)
  dco2d(0) = 0.0_8
  dco2(0) = 789.*co2*dp(0)
  dch4d(0) = 0.0_8
  dch4(0) = 789.*ch4_dev(1)*dp(0)
  dn2od(0) = 0.0_8
  dn2o(0) = 789.*n2o_dev(1)*dp(0)
  df11d(0) = 0.0_8
  df11(0) = 789.*cfc11_dev(1)*dp(0)
  df12d(0) = 0.0_8
  df12(0) = 789.*cfc12_dev(1)*dp(0)
  df22d(0) = 0.0_8
  df22(0) = 789.*cfc22_dev(1)*dp(0)
  IF (dh2o(0) .LT. 1.e-10) THEN
    dh2od(0) = 0.0_8
    dh2o(0) = 1.e-10
  ELSE
    dh2o(0) = dh2o(0)
  END IF
  IF (do3(0) .LT. 1.e-6) THEN
    do3d(0) = 0.0_8
    do3(0) = 1.e-6
  ELSE
    do3(0) = do3(0)
  END IF
  IF (dco2(0) .LT. 1.e-4) THEN
    dco2d(0) = 0.0_8
    dco2(0) = 1.e-4
  ELSE
    dco2d(0) = 0.0_8
    dco2(0) = dco2(0)
  END IF
  xxd = pa(0)*0.001618*dp(0)*(wa_devd(1)*wa_dev(1)+wa_dev(1)*wa_devd(1))
  xx = pa(0)*0.001618*wa_dev(1)*wa_dev(1)*dp(0)
  dcontd(0) = xxd*EXP(1800./ta_dev(1)-6.081) - xx*1800.*ta_devd(1)*EXP(&
&   1800./ta_dev(1)-6.081)/ta_dev(1)**2
  dcont(0) = xx*EXP(1800./ta_dev(1)-6.081)
!The surface (np+1) is treated as a layer filled with black clouds.
!transfc is the transmittance between the surface and a pressure
!level.
  transfcd(np+1) = 0.0_8
  transfc(np+1) = 1.0
!Initialize fluxes
  DO k=1,np+1
    flxu_devd(k) = 0.0_8
    flxu_dev(k) = 0.0
    flxd_devd(k) = 0.0_8
    flxd_dev(k) = 0.0
    dfdts_devd(k) = 0.0_8
    dfdts_dev(k) = 0.0
  END DO
  n2oexp_tmpd = 0.0_8
  blayerd = 0.0_8
  tn2od = 0.0_8
  transfcd = 0.0_8
  bdd = 0.0_8
  h2oexp_tmpd = 0.0_8
  tcomd = 0.0_8
  tf11d = 0.0_8
  tf12d = 0.0_8
  trantd = 0.0_8
  bud = 0.0_8
  th2od = 0.0_8
  tcldlyrd = 0.0_8
  tch4d = 0.0_8
  tf22d = 0.0_8
  ennd = 0.0_8
  fclrd = 0.0_8
  tco2d = 0.0_8
  co2exp_tmpd = 0.0_8
  taerlyrd = 0.0_8
  bleveld = 0.0_8
!Integration over spectral bands
  DO ibn=1,10
    IF (ibn .EQ. 10 .AND. (.NOT.trace)) THEN
      GOTO 100
    ELSE
!if h2otable, compute h2o (line) transmittance using table look-up.
!if conbnd,   compute h2o (continuum) transmittance in bands 2-7.
!if co2bnd,   compute co2 transmittance in band 3.
!if oznbnd,   compute  o3 transmittance in band 5.
!if n2obnd,   compute n2o transmittance in bands 6 and 7.
!if ch4bnd,   compute ch4 transmittance in bands 6 and 7.
!if combnd,   compute co2-minor transmittance in bands 4 and 5.
!if f11bnd,   compute cfc11 transmittance in bands 4 and 5.
!if f12bnd,   compute cfc12 transmittance in bands 4 and 6.
!if f22bnd,   compute cfc22 transmittance in bands 4 and 6.
!if b10bnd,   compute flux reduction due to n2o in band 10.
      h2otable = (ibn .EQ. 1 .OR. ibn .EQ. 2) .OR. ibn .EQ. 8
      conbnd = ibn .GE. 2 .AND. ibn .LE. 7
      co2bnd = ibn .EQ. 3
      oznbnd = ibn .EQ. 5
      n2obnd = ibn .EQ. 6 .OR. ibn .EQ. 7
      ch4bnd = ibn .EQ. 6 .OR. ibn .EQ. 7
      combnd = ibn .EQ. 4 .OR. ibn .EQ. 5
      f11bnd = ibn .EQ. 4 .OR. ibn .EQ. 5
      f12bnd = ibn .EQ. 4 .OR. ibn .EQ. 6
      f22bnd = ibn .EQ. 4 .OR. ibn .EQ. 6
      b10bnd = ibn .EQ. 10
      do_aerosol = na .GT. 0
      exptbl = 0.0
!Control packing of the new exponential tables by band
      SELECT CASE  (ibn) 
      CASE (2) 
        conexp%start = 1
        conexp%end = 1
      CASE (3) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 9
      CASE (4) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 7
        comexp%start = 8
        comexp%end = 13
        f11exp%start = 14
        f11exp%end = 14
        f12exp%start = 15
        f12exp%end = 15
        f22exp%start = 16
        f22exp%end = 16
      CASE (5) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 7
        comexp%start = 8
        comexp%end = 13
        f11exp%start = 14
        f11exp%end = 14
      CASE (6) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 7
        n2oexp%start = 8
        n2oexp%end = 11
        ch4exp%start = 12
        ch4exp%end = 15
        f12exp%start = 16
        f12exp%end = 16
        f22exp%start = 17
        f22exp%end = 17
      CASE (7) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 7
        n2oexp%start = 8
        n2oexp%end = 11
        ch4exp%start = 12
        ch4exp%end = 15
      CASE (9) 
        h2oexp%start = 1
        h2oexp%end = 6
      CASE (10) 
        h2oexp%start = 1
        h2oexp%end = 5
        conexp%start = 6
        conexp%end = 6
        co2exp%start = 7
        co2exp%end = 12
        n2oexp%start = 13
        n2oexp%end = 14
      END SELECT
!blayer is the spectrally integrated planck flux of the mean layer
!temperature derived from eq. (3.11)
!The fitting for the planck flux is valid for the range 160-345 K.
      DO k=1,np
        CALL PLANCK_D(ibn, cb, ta_dev(k), ta_devd(k), blayer(k), blayerd&
&               (k))
      END DO
!Index "0" is the layer above the top of the atmosphere.
      blayerd(0) = blayerd(1)
      blayer(0) = blayer(1)
      bleveld(0) = blayerd(1)
      blevel(0) = blayer(1)
!Surface emission and reflectivity. See Section 9.
!bs and dbs include the effect of surface emissivity.
      CALL SFCFLUX(ibn, cb, dcb, ns, fs_dev, tg_dev, eg_dev, tv_dev, &
&            ev_dev, rv_dev, bs, dbs, rflxs)
      blayerd(np+1) = 0.0_8
      blayer(np+1) = bs
!interpolate Planck function at model levels (linear in p)
      DO k=2,np
        bleveld(k) = (dp(k)*blayerd(k-1)+dp(k-1)*blayerd(k))/(dp(k-1)+dp&
&         (k))
        blevel(k) = (blayer(k-1)*dp(k)+blayer(k)*dp(k-1))/(dp(k-1)+dp(k)&
&         )
      END DO
!Extrapolate blevel(1) from blayer(2) and blayer(1)
      bleveld(1) = blayerd(1) + dp(1)*(blayerd(1)-blayerd(2))/(dp(1)+dp(&
&       2))
      blevel(1) = blayer(1) + (blayer(1)-blayer(2))*dp(1)/(dp(1)+dp(2))
      bleveld(0) = bleveld(1)
      blevel(0) = blevel(1)
!If the surface air temperature tb is known, compute blevel(np+1)
      CALL PLANCK_D(ibn, cb, tb_dev, tb_devd, blevel(np+1), bleveld(np+1&
&             ))
!Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!     NOTE: dp_pa is only dims(1:np) as the 0'th level isn't needed in getirtau.
!           Plus, the pressures in getirtau *MUST* be in Pascals.
!     Slots for reff, hydrometeors and tauall are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
      CALL GETIRTAU1_D(ibn, np, dp_pa, fcld_col, fcld_cold, reff_col, &
&                reff_cold, cwc_col, cwc_cold, tcldlyr, tcldlyrd, enn, &
&                ennd, aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, awg_ir, &
&                cons_grav)
      DO k=0,np
        icx(k) = k
      END DO
      CALL MKICX(np, ict, icb, enn, icx, ncld)
!Compute optical thickness, single-scattering albedo and asymmetry
!factor for a mixture of "na" aerosol types. Eqs. (7.1)-(7.3)
      IF (do_aerosol) THEN
        taerlyrd(0) = 0.0_8
        taerlyr(0) = 1.0
        DO k=1,np
!-----taerlyr is the aerosol diffuse transmittance
          taerlyrd(k) = 0.0_8
          taerlyr(k) = 1.0
          IF (taua_dev(k, ibn) .GT. 0.001) THEN
            IF (ssaa_dev(k, ibn) .GT. 0.001) THEN
              asya_devd(k, ibn) = (asya_devd(k, ibn)*ssaa_dev(k, ibn)-&
&               asya_dev(k, ibn)*ssaa_devd(k, ibn))/ssaa_dev(k, ibn)**2
              asya_dev(k, ibn) = asya_dev(k, ibn)/ssaa_dev(k, ibn)
              ssaa_devd(k, ibn) = (ssaa_devd(k, ibn)*taua_dev(k, ibn)-&
&               ssaa_dev(k, ibn)*taua_devd(k, ibn))/taua_dev(k, ibn)**2
              ssaa_dev(k, ibn) = ssaa_dev(k, ibn)/taua_dev(k, ibn)
!Parameterization of aerosol scattering following
              ffd = (0.1185*asya_devd(k, ibn)*asya_dev(k, ibn)+(0.0076+&
&               0.1185*asya_dev(k, ibn))*asya_devd(k, ibn))*asya_dev(k, &
&               ibn) + (.3739+(0.0076+0.1185*asya_dev(k, ibn))*asya_dev(&
&               k, ibn))*asya_devd(k, ibn)
              ff = .5 + (.3739+(0.0076+0.1185*asya_dev(k, ibn))*asya_dev&
&               (k, ibn))*asya_dev(k, ibn)
              taua_devd(k, ibn) = taua_devd(k, ibn)*(1.-ssaa_dev(k, ibn)&
&               *ff) + taua_dev(k, ibn)*(-(ssaa_devd(k, ibn)*ff)-&
&               ssaa_dev(k, ibn)*ffd)
              taua_dev(k, ibn) = taua_dev(k, ibn)*(1.-ssaa_dev(k, ibn)*&
&               ff)
            END IF
            taerlyrd(k) = -(1.66*taua_devd(k, ibn)*EXP(-(1.66*taua_dev(k&
&             , ibn))))
            taerlyr(k) = EXP(-(1.66*taua_dev(k, ibn)))
          END IF
        END DO
      END IF
!Compute the exponential terms (Eq. 8.21) at each layer due to
!water vapor line absorption when k-distribution is used
      IF (.NOT.h2otable .AND. (.NOT.b10bnd)) THEN
        CALL H2OEXPS_D(ibn, np, dh2o, dh2od, pa, dt, dtd, xkw, aw, bw, &
&                pm, mw, exptbl(:, h2oexp%start:h2oexp%end), exptbld(:, &
&                h2oexp%start:h2oexp%end))
      ELSE
        exptbld = 0.0_8
      END IF
!compute the exponential terms (Eq. 4.24) at each layer due to
!water vapor continuum absorption.
      ne = 0
      IF (conbnd) THEN
        ne = 1
        IF (ibn .EQ. 3) ne = 3
        CALL CONEXPS_D(ibn, np, dcont, dcontd, xke, exptbl(:, conexp%&
&                start:conexp%end), exptbld(:, conexp%start:conexp%end))
      END IF
      IF (trace) THEN
!compute the exponential terms at each layer due to n2o absorption
        IF (n2obnd) CALL N2OEXPS_D(ibn, np, dn2o, pa, dt, dtd, exptbl(:&
&                            , n2oexp%start:n2oexp%end), exptbld(:, &
&                            n2oexp%start:n2oexp%end))
!compute the exponential terms at each layer due to ch4 absorption
        IF (ch4bnd) CALL CH4EXPS_D(ibn, np, dch4, pa, dt, dtd, exptbl(:&
&                            , ch4exp%start:ch4exp%end), exptbld(:, &
&                            ch4exp%start:ch4exp%end))
!Compute the exponential terms due to co2 minor absorption
        IF (combnd) CALL COMEXPS_D(ibn, np, dco2, dt, dtd, exptbl(:, &
&                            comexp%start:comexp%end), exptbld(:, comexp&
&                            %start:comexp%end))
!Compute the exponential terms due to cfc11 absorption.
!The values of the parameters are given in Table 7.
        IF (f11bnd) THEN
          a1 = 1.26610e-3
          b1 = 3.55940e-6
          fk1 = 1.89736e+1
          a2 = 8.19370e-4
          b2 = 4.67810e-6
          fk2 = 1.01487e+1
          CALL CFCEXPS_D(ibn, np, a1, b1, fk1, a2, b2, fk2, df11, dt, &
&                  dtd, exptbl(:, f11exp%start:f11exp%end), exptbld(:, &
&                  f11exp%start:f11exp%end))
        END IF
!Compute the exponential terms due to cfc12 absorption.
        IF (f12bnd) THEN
          a1 = 8.77370e-4
          b1 = -5.88440e-6
          fk1 = 1.58104e+1
          a2 = 8.62000e-4
          b2 = -4.22500e-6
          fk2 = 3.70107e+1
          CALL CFCEXPS_D(ibn, np, a1, b1, fk1, a2, b2, fk2, df12, dt, &
&                  dtd, exptbl(:, f12exp%start:f12exp%end), exptbld(:, &
&                  f12exp%start:f12exp%end))
        END IF
!Compute the exponential terms due to cfc22 absorption.
        IF (f22bnd) THEN
          a1 = 9.65130e-4
          b1 = 1.31280e-5
          fk1 = 6.18536e+0
          a2 = -3.00010e-5
          b2 = 5.25010e-7
          fk2 = 3.27912e+1
          CALL CFCEXPS_D(ibn, np, a1, b1, fk1, a2, b2, fk2, df22, dt, &
&                  dtd, exptbl(:, f22exp%start:f22exp%end), exptbld(:, &
&                  f22exp%start:f22exp%end))
        END IF
!Compute the exponential terms at each layer in band 10 due to
!h2o line and continuum, co2, and n2o absorption
        IF (b10bnd) THEN
          CALL B10EXPS_D(np, dh2o, dh2od, dcont, dcontd, dco2, dn2o, pa&
&                  , dt, dtd, h2oexp_tmp, h2oexp_tmpd, exptbl(:, conexp%&
&                  start:conexp%end), exptbld(:, conexp%start:conexp%end&
&                  ), co2exp_tmp, co2exp_tmpd, n2oexp_tmp, n2oexp_tmpd)
          exptbld(:, h2oexp%start:h2oexp%end) = h2oexp_tmpd
          exptbl(:, h2oexp%start:h2oexp%end) = h2oexp_tmp
          exptbld(:, co2exp%start:co2exp%end) = co2exp_tmpd
          exptbl(:, co2exp%start:co2exp%end) = co2exp_tmp
          exptbld(:, n2oexp%start:n2oexp%end) = n2oexp_tmpd
          exptbl(:, n2oexp%start:n2oexp%end) = n2oexp_tmp
        END IF
      END IF
!blayer(np+1) includes the effect of surface emissivity.
! ALT: this was undefined, check with Max if 0.0 is good value
      bud(0) = 0.0_8
      bu(0) = 0.0
      bdd(0) = blayerd(1)
      bd(0) = blayer(1)
      bud(np+1) = blayerd(np+1)
      bu(np+1) = blayer(np+1)
!do-loop 1500 is for computing upward (bu) and downward (bd)
!Here, trant is the transmittance of the layer k2-1.
      DO k2=1,np+1
!for h2o line transmission
        IF (.NOT.h2otable) THEN
          th2o = 1.0
          th2od = 0.0_8
        END IF
!for h2o continuum transmission
        tcon = 1.0
        x1 = 0.0
        x2 = 0.0
        x3 = 0.0
        trant = 1.0
        IF (h2otable) THEN
!Compute water vapor transmittance using table look-up.
          IF (ibn .EQ. 1) THEN
            trantd = 0.0_8
            x3d = 0.0_8
            x2d = 0.0_8
            x1d = 0.0_8
            CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2-1), pa(k2-1), &
&                   dt(k2-1), dtd(k2-1), x1, x1d, x2, x2d, x3, x3d, w11&
&                   , p11, dwe, dpe, h11, h12, h13, trant, trantd)
          ELSE
            trantd = 0.0_8
            x1d = 0.0_8
            x2d = 0.0_8
            x3d = 0.0_8
          END IF
          IF (ibn .EQ. 2) CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2-1&
&                                 ), pa(k2-1), dt(k2-1), dtd(k2-1), x1, &
&                                 x1d, x2, x2d, x3, x3d, w11, p11, dwe, &
&                                 dpe, h21, h22, h23, trant, trantd)
          IF (ibn .EQ. 8) CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2-1&
&                                 ), pa(k2-1), dt(k2-1), dtd(k2-1), x1, &
&                                 x1d, x2, x2d, x3, x3d, w11, p11, dwe, &
&                                 dpe, h81, h82, h83, trant, trantd)
!for water vapor continuum absorption
          IF (conbnd) THEN
! Only the first exp
            tcond = 0.0_8
            tcond(1) = tcon(1)*exptbld(k2-1, conexp%start)
            tcon(1) = tcon(1)*exptbl(k2-1, conexp%start)
            trantd = trantd*tcon(1) + trant*tcond(1)
            trant = trant*tcon(1)
          END IF
        ELSE IF (.NOT.b10bnd) THEN
!compute water vapor transmittance using k-distribution
          tcond = 0.0_8
          CALL H2OKDIS_D(ibn, np, k2 - 1, fkw, gkw, ne, exptbl(:, h2oexp&
&                  %start:h2oexp%end), exptbld(:, h2oexp%start:h2oexp%&
&                  end), exptbl(:, conexp%start:conexp%end), exptbld(:, &
&                  conexp%start:conexp%end), th2o, th2od, tcon, tcond, &
&                  trant, trantd)
          x1d = 0.0_8
          x2d = 0.0_8
          x3d = 0.0_8
        ELSE
          trantd = 0.0_8
          x1d = 0.0_8
          x2d = 0.0_8
          x3d = 0.0_8
        END IF
        IF (co2bnd) THEN
!Compute co2 transmittance using table look-up method
          dco2d = 0.0_8
          CALL TABLUP_D(nx1, nc1, dco2(k2-1), dco2d(k2-1), pa(k2-1), dt(&
&                 k2-1), dtd(k2-1), x1, x1d, x2, x2d, x3, x3d, w12, p12&
&                 , dwe, dpe, c1, c2, c3, trant, trantd)
        END IF
!Always use table look-up to compute o3 transmittance.
        IF (oznbnd) CALL TABLUP_D(nx1, no1, do3(k2-1), do3d(k2-1), pa(k2&
&                           -1), dt(k2-1), dtd(k2-1), x1, x1d, x2, x2d, &
&                           x3, x3d, w13, p13, dwe, dpe, oo1, oo2, oo3, &
&                           trant, trantd)
!include aerosol effect
        IF (do_aerosol) THEN
          trantd = trantd*taerlyr(k2-1) + trant*taerlyrd(k2-1)
          trant = trant*taerlyr(k2-1)
        END IF
!Compute upward (bu) and downward (bd) emission of the layer k2-1
        xxd = (1.-enn(k2-1))*trantd - ennd(k2-1)*trant
        xx = (1.-enn(k2-1))*trant
        IF (0.9999 .GT. xx) THEN
          yyd = xxd
          yy = xx
        ELSE
          yy = 0.9999
          yyd = 0.0_8
        END IF
        IF (0.00001 .LT. yy) THEN
          yy = yy
        ELSE
          yy = 0.00001
          yyd = 0.0_8
        END IF
        xxd = ((bleveld(k2-1)-bleveld(k2))*LOG(yy)-(blevel(k2-1)-blevel(&
&         k2))*yyd/yy)/LOG(yy)**2
        xx = (blevel(k2-1)-blevel(k2))/LOG(yy)
        bdd(k2-1) = ((bleveld(k2)-bleveld(k2-1)*yy-blevel(k2-1)*yyd)*(&
&         1.0-yy)+(blevel(k2)-blevel(k2-1)*yy)*yyd)/(1.0-yy)**2 - xxd
        bd(k2-1) = (blevel(k2)-blevel(k2-1)*yy)/(1.0-yy) - xx
        bud(k2-1) = bleveld(k2-1) + bleveld(k2) - bdd(k2-1)
        bu(k2-1) = blevel(k2-1) + blevel(k2) - bd(k2-1)
      END DO
!Initialize fluxes
      flxu = 0.0
      flxd = 0.0
      flxud = 0.0_8
      flxdd = 0.0_8
!Compute upward and downward fluxes for each spectral band, ibn.
      DO k1=0,np
!initialization
        cldlw = 0.0
        cldmd = 0.0
        cldhi = 0.0
        tranal = 1.0
!for h2o line transmission
        IF (.NOT.h2otable) THEN
          th2o = 1.0
          th2od = 0.0_8
        END IF
!for h2o continuum transmission
        tcon = 1.0
        IF (trace) THEN
!Add trace gases contribution
          IF (n2obnd) THEN
!n2o
            tn2o = 1.0
            tn2od = 0.0_8
          END IF
          IF (ch4bnd) THEN
!ch4
            tch4 = 1.0
            tch4d = 0.0_8
          END IF
          IF (combnd) THEN
!co2-minor
            tcom = 1.0
            tcomd = 0.0_8
          END IF
          IF (f11bnd) THEN
!cfc-11
            tf11 = 1.0
            tf11d = 0.0_8
          END IF
          IF (f12bnd) THEN
!cfc-12
            tf12 = 1.0
            tf12d = 0.0_8
          END IF
          IF (f22bnd) THEN
!cfc-22
            tf22 = 1.0
            tf22d = 0.0_8
          END IF
          IF (b10bnd) THEN
!
            th2o = 1.0
            tco2 = 1.0
            tcond(1) = 0.0_8
            tcon(1) = 1.0
            tn2o = 1.0
            tn2od = 0.0_8
            th2od = 0.0_8
            tco2d = 0.0_8
          END IF
        END IF
        x1 = 0.0
        x2 = 0.0
        x3 = 0.0
        cldhid = 0.0_8
        tcond = 0.0_8
        tranald = 0.0_8
        x1d = 0.0_8
        x2d = 0.0_8
        x3d = 0.0_8
        cldlwd = 0.0_8
        cldmdd = 0.0_8
        DO k2=k1+1,np+1
          trant = 1.0
          fclr = 1.0
          IF (h2otable) THEN
!Compute water vapor transmittance using table look-up.
            IF (ibn .EQ. 1) THEN
              trantd = 0.0_8
              CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2-1), pa(k2-1)&
&                     , dt(k2-1), dtd(k2-1), x1, x1d, x2, x2d, x3, x3d, &
&                     w11, p11, dwe, dpe, h11, h12, h13, trant, trantd)
            ELSE
              trantd = 0.0_8
            END IF
            IF (ibn .EQ. 2) CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2&
&                                   -1), pa(k2-1), dt(k2-1), dtd(k2-1), &
&                                   x1, x1d, x2, x2d, x3, x3d, w11, p11&
&                                   , dwe, dpe, h21, h22, h23, trant, &
&                                   trantd)
            IF (ibn .EQ. 8) CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2&
&                                   -1), pa(k2-1), dt(k2-1), dtd(k2-1), &
&                                   x1, x1d, x2, x2d, x3, x3d, w11, p11&
&                                   , dwe, dpe, h81, h82, h83, trant, &
&                                   trantd)
            IF (conbnd) THEN
! Only the first exp
              tcond(1) = tcond(1)*exptbl(k2-1, conexp%start) + tcon(1)*&
&               exptbld(k2-1, conexp%start)
              tcon(1) = tcon(1)*exptbl(k2-1, conexp%start)
              trantd = trantd*tcon(1) + trant*tcond(1)
              trant = trant*tcon(1)
            END IF
          ELSE IF (.NOT.b10bnd) THEN
!Compute water vapor transmittance using k-distribution
            CALL H2OKDIS_D(ibn, np, k2 - 1, fkw, gkw, ne, exptbl(:, &
&                    h2oexp%start:h2oexp%end), exptbld(:, h2oexp%start:&
&                    h2oexp%end), exptbl(:, conexp%start:conexp%end), &
&                    exptbld(:, conexp%start:conexp%end), th2o, th2od, &
&                    tcon, tcond, trant, trantd)
          ELSE
            trantd = 0.0_8
          END IF
          IF (co2bnd) THEN
!Compute co2 transmittance using table look-up method.
            dco2d = 0.0_8
            CALL TABLUP_D(nx1, nc1, dco2(k2-1), dco2d(k2-1), pa(k2-1), &
&                   dt(k2-1), dtd(k2-1), x1, x1d, x2, x2d, x3, x3d, w12&
&                   , p12, dwe, dpe, c1, c2, c3, trant, trantd)
          END IF
          IF (oznbnd) CALL TABLUP_D(nx1, no1, do3(k2-1), do3d(k2-1), pa(&
&                             k2-1), dt(k2-1), dtd(k2-1), x1, x1d, x2, &
&                             x2d, x3, x3d, w13, p13, dwe, dpe, oo1, oo2&
&                             , oo3, trant, trantd)
!Always use table look-up to compute o3 transmittanc
          IF (trace) THEN
!Add trace gas effects
            IF (n2obnd) CALL N2OKDIS_D(ibn, np, k2 - 1, exptbl(:, n2oexp&
&                                %start:n2oexp%end), exptbld(:, n2oexp%&
&                                start:n2oexp%end), tn2o, tn2od, trant, &
&                                trantd)
!n2o
            IF (ch4bnd) CALL CH4KDIS_D(ibn, np, k2 - 1, exptbl(:, ch4exp&
&                                %start:ch4exp%end), exptbld(:, ch4exp%&
&                                start:ch4exp%end), tch4, tch4d, trant, &
&                                trantd)
!ch4
            IF (combnd) CALL COMKDIS_D(ibn, np, k2 - 1, exptbl(:, comexp&
&                                %start:comexp%end), exptbld(:, comexp%&
&                                start:comexp%end), tcom, tcomd, trant, &
&                                trantd)
!co2-minor
            IF (f11bnd) CALL CFCKDIS_D(np, k2 - 1, exptbl(:, f11exp%&
&                                start:f11exp%end), exptbld(:, f11exp%&
&                                start:f11exp%end), tf11, tf11d, trant, &
&                                trantd)
!cfc11
            IF (f12bnd) CALL CFCKDIS_D(np, k2 - 1, exptbl(:, f12exp%&
&                                start:f12exp%end), exptbld(:, f12exp%&
&                                start:f12exp%end), tf12, tf12d, trant, &
&                                trantd)
!cfc12
            IF (f22bnd) CALL CFCKDIS_D(np, k2 - 1, exptbl(:, f22exp%&
&                                start:f22exp%end), exptbld(:, f22exp%&
&                                start:f22exp%end), tf22, tf22d, trant, &
&                                trantd)
!CFC22
            IF (b10bnd) CALL B10KDIS_D(np, k2 - 1, exptbl(:, h2oexp%&
&                                start:h2oexp%end), exptbld(:, h2oexp%&
&                                start:h2oexp%end), exptbl(:, conexp%&
&                                start:conexp%end), exptbld(:, conexp%&
&                                start:conexp%end), exptbl(:, co2exp%&
&                                start:co2exp%end), exptbld(:, co2exp%&
&                                start:co2exp%end), exptbl(:, n2oexp%&
&                                start:n2oexp%end), exptbld(:, n2oexp%&
&                                start:n2oexp%end), th2o, th2od, tcon, &
&                                tcond, tco2, tco2d, tn2o, tn2od, trant&
&                                , trantd)
          END IF
          IF (do_aerosol) THEN
            tranald = tranald*taerlyr(k2-1) + tranal*taerlyrd(k2-1)
            tranal = tranal*taerlyr(k2-1)
            trantd = trantd*tranal + trant*tranald
            trant = trant*tranal
          END IF
          IF (enn(k2-1) .GE. 0.001) CALL CLDOVLP_D(np, k1, k2, ict, icb&
&                                            , icx, ncld, enn, ennd, &
&                                            tcldlyr, tcldlyrd, cldhi, &
&                                            cldhid, cldmd, cldmdd, &
&                                            cldlw, cldlwd)
          fclrd = (-(cldhid*(1.0-cldmd))-(1.0-cldhi)*cldmdd)*(1.0-cldlw)&
&           - (1.0-cldhi)*(1.0-cldmd)*cldlwd
          fclr = (1.0-cldhi)*(1.0-cldmd)*(1.0-cldlw)
          IF (k2 .EQ. k1 + 1 .AND. ibn .NE. 10) THEN
            flxud(k1) = flxud(k1) - bud(k1)
            flxu(k1) = flxu(k1) - bu(k1)
            flxdd(k2) = flxdd(k2) + bdd(k1)
            flxd(k2) = flxd(k2) + bd(k1)
          END IF
          xxd = trantd*(bu(k2-1)-bu(k2)) + trant*(bud(k2-1)-bud(k2))
          xx = trant*(bu(k2-1)-bu(k2))
          flxud(k1) = flxud(k1) + xxd*fclr + xx*fclrd
          flxu(k1) = flxu(k1) + xx*fclr
          IF (k1 .EQ. 0) THEN
!mjs  bd(-1) is not defined
            xxd = -(trantd*bd(k1)+trant*bdd(k1))
            xx = -(trant*bd(k1))
          ELSE
            xxd = trantd*(bd(k1-1)-bd(k1)) + trant*(bdd(k1-1)-bdd(k1))
            xx = trant*(bd(k1-1)-bd(k1))
          END IF
          flxdd(k2) = flxdd(k2) + xxd*fclr + xx*fclrd
          flxd(k2) = flxd(k2) + xx*fclr
        END DO
        transfcd(k1) = trantd*fclr + trant*fclrd
        transfc(k1) = trant*fclr
        IF (k1 .GT. 0) THEN
          dfdts_devd(k1) = dfdts_devd(k1) - dbs*transfcd(k1)
          dfdts_dev(k1) = dfdts_dev(k1) - dbs*transfc(k1)
        END IF
      END DO
      IF (.NOT.b10bnd) THEN
!surface emission
        flxud(np+1) = -blayerd(np+1)
        flxu(np+1) = -blayer(np+1)
        dfdts_dev(np+1) = dfdts_dev(np+1) - dbs
        DO k=1,np+1
          flxud(k) = flxud(k) - rflxs*(flxdd(np+1)*transfc(k)+flxd(np+1)&
&           *transfcd(k))
          flxu(k) = flxu(k) - flxd(np+1)*transfc(k)*rflxs
        END DO
      END IF
!Summation of fluxes over spectral bands
      DO k=1,np+1
        flxu_devd(k) = flxu_devd(k) + flxud(k)
        flxu_dev(k) = flxu_dev(k) + flxu(k)
        flxd_devd(k) = flxd_devd(k) + flxdd(k)
        flxd_dev(k) = flxd_dev(k) + flxd(k)
      END DO
    END IF
  END DO
  GOTO 110
 100 RETURN
 110 CONTINUE
END SUBROUTINE IRRAD_D

!  Differentiation of planck in forward (tangent) mode:
!   variations   of useful results: xlayer
!   with respect to varying inputs: t
!***********************************************************************
SUBROUTINE PLANCK_D(ibn, cb, t, td, xlayer, xlayerd)
  IMPLICIT NONE
! spectral band index
  INTEGER, INTENT(IN) :: ibn
! Planck table coefficients
  REAL*8, INTENT(IN) :: cb(6, 10)
! temperature (K)
  REAL*8, INTENT(IN) :: t
  REAL*8, INTENT(IN) :: td
! planck flux (w/m2)
  REAL*8, INTENT(OUT) :: xlayer
  REAL*8, INTENT(OUT) :: xlayerd
  xlayerd = td*(t*(t*(t*(t*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+cb(3, ibn)&
&   )+cb(2, ibn)) + t*(td*(t*(t*(t*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+cb&
&   (3, ibn))+t*(td*(t*(t*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+t*(td*(t*cb&
&   (6, ibn)+cb(5, ibn))+t*cb(6, ibn)*td)))
  xlayer = t*(t*(t*(t*(t*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+cb(3, ibn))+&
&   cb(2, ibn)) + cb(1, ibn)
END SUBROUTINE PLANCK_D

!  Differentiation of h2oexps in forward (tangent) mode:
!   variations   of useful results: h2oexp
!   with respect to varying inputs: dh2o dt
!**********************************************************************
SUBROUTINE H2OEXPS_D(ib, np, dh2o, dh2od, pa, dt, dtd, xkw, aw, bw, pm, &
& mw, h2oexp, h2oexpd)
  IMPLICIT NONE
  INTEGER :: ib, np, ik, k
!---- input parameters ------
  REAL*8 :: dh2o(0:np), pa(0:np), dt(0:np)
  REAL*8 :: dh2od(0:np), dtd(0:np)
!---- output parameters -----
  REAL*8 :: h2oexp(0:np, 6)
  REAL*8 :: h2oexpd(0:np, 6)
!---- static data -----
  INTEGER :: mw(9)
  REAL*8 :: xkw(9), aw(9), bw(9), pm(9)
!---- temporary arrays -----
  REAL*8 :: xh
  REAL*8 :: xhd
  INTRINSIC EXP
  REAL*8 :: pwx1
  REAL*8 :: pwr1
  h2oexpd = 0.0_8
!**********************************************************************
!    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
!    and bw,  therefore, h2oexp for these sub-bands are identical.
!**********************************************************************
  DO k=0,np
!-----xh is the scaled water vapor amount for line absorption
!     computed from Eq. (4.4).
    pwx1 = pa(k)/500.
    pwr1 = pwx1**pm(ib)
    xhd = pwr1*(dh2od(k)*(1.+(aw(ib)+bw(ib)*dt(k))*dt(k))+dh2o(k)*(bw(ib&
&     )*dtd(k)*dt(k)+(aw(ib)+bw(ib)*dt(k))*dtd(k)))
    xh = dh2o(k)*pwr1*(1.+(aw(ib)+bw(ib)*dt(k))*dt(k))
!-----h2oexp is the water vapor transmittance of the layer k
!     due to line absorption
    h2oexpd(k, 1) = -(xkw(ib)*xhd*EXP(-(xh*xkw(ib))))
    h2oexp(k, 1) = EXP(-(xh*xkw(ib)))
!-----compute transmittances from Eq. (8.22)
    DO ik=2,6
      IF (mw(ib) .EQ. 6) THEN
        xhd = h2oexpd(k, ik-1)*h2oexp(k, ik-1) + h2oexp(k, ik-1)*h2oexpd&
&         (k, ik-1)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        h2oexpd(k, ik) = (xhd*xh+xh*xhd)*xh + xh**2*xhd
        h2oexp(k, ik) = xh*xh*xh
      ELSE IF (mw(ib) .EQ. 8) THEN
        xhd = h2oexpd(k, ik-1)*h2oexp(k, ik-1) + h2oexp(k, ik-1)*h2oexpd&
&         (k, ik-1)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        xhd = xhd*xh + xh*xhd
        xh = xh*xh
        h2oexpd(k, ik) = xhd*xh + xh*xhd
        h2oexp(k, ik) = xh*xh
      ELSE IF (mw(ib) .EQ. 9) THEN
        xhd = (h2oexpd(k, ik-1)*h2oexp(k, ik-1)+h2oexp(k, ik-1)*h2oexpd(&
&         k, ik-1))*h2oexp(k, ik-1) + h2oexp(k, ik-1)**2*h2oexpd(k, ik-1&
&         )
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)*h2oexp(k, ik-1)
        h2oexpd(k, ik) = (xhd*xh+xh*xhd)*xh + xh**2*xhd
        h2oexp(k, ik) = xh*xh*xh
      ELSE
        xhd = h2oexpd(k, ik-1)*h2oexp(k, ik-1) + h2oexp(k, ik-1)*h2oexpd&
&         (k, ik-1)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        xhd = xhd*xh + xh*xhd
        xh = xh*xh
        xhd = xhd*xh + xh*xhd
        xh = xh*xh
        h2oexpd(k, ik) = xhd*xh + xh*xhd
        h2oexp(k, ik) = xh*xh
      END IF
    END DO
  END DO
END SUBROUTINE H2OEXPS_D

!  Differentiation of conexps in forward (tangent) mode:
!   variations   of useful results: conexp
!   with respect to varying inputs: dcont conexp
!**********************************************************************
SUBROUTINE CONEXPS_D(ib, np, dcont, dcontd, xke, conexp, conexpd)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters ------
  REAL*8 :: dcont(0:np)
  REAL*8 :: dcontd(0:np)
!---- updated parameters -----
  REAL*8 :: conexp(0:np, 3)
  REAL*8 :: conexpd(0:np, 3)
!---- static data -----
  REAL*8 :: xke(9)
  INTRINSIC EXP
!****************************************************************
  DO k=0,np
    conexpd(k, 1) = -(xke(ib)*dcontd(k)*EXP(-(dcont(k)*xke(ib))))
    conexp(k, 1) = EXP(-(dcont(k)*xke(ib)))
!-----The absorption coefficients for sub-bands 3b and 3a are, respectively,
!     two and four times the absorption coefficient for sub-band 3c (Table 9).
!     Note that conexp(3) is for sub-band 3a. 
    IF (ib .EQ. 3) THEN
      conexpd(k, 2) = conexpd(k, 1)*conexp(k, 1) + conexp(k, 1)*conexpd(&
&       k, 1)
      conexp(k, 2) = conexp(k, 1)*conexp(k, 1)
      conexpd(k, 3) = conexpd(k, 2)*conexp(k, 2) + conexp(k, 2)*conexpd(&
&       k, 2)
      conexp(k, 3) = conexp(k, 2)*conexp(k, 2)
    END IF
  END DO
END SUBROUTINE CONEXPS_D

!  Differentiation of n2oexps in forward (tangent) mode:
!   variations   of useful results: n2oexp
!   with respect to varying inputs: dt n2oexp
!**********************************************************************
SUBROUTINE N2OEXPS_D(ib, np, dn2o, pa, dt, dtd, n2oexp, n2oexpd)
  IMPLICIT NONE
  INTEGER :: ib, k, np
!---- input parameters -----
  REAL*8 :: dn2o(0:np), pa(0:np), dt(0:np)
  REAL*8 :: dtd(0:np)
!---- output parameters -----
  REAL*8 :: n2oexp(0:np, 4)
  REAL*8 :: n2oexpd(0:np, 4)
!---- temporary arrays -----
  REAL*8 :: xc, xc1, xc2
  REAL*8 :: xcd, xc1d, xc2d
  INTRINSIC EXP
!-----Scaling and absorption data are given in Table 5.
!     Transmittances are computed using Eqs. (8.21) and (8.22).
  DO k=0,np
!-----four exponential by powers of 21 for band 6.
    IF (ib .EQ. 6) THEN
      xcd = dn2o(k)*(4.3750e-6*dtd(k)*dt(k)+(1.9297e-3+4.3750e-6*dt(k))*&
&       dtd(k))
      xc = dn2o(k)*(1.+(1.9297e-3+4.3750e-6*dt(k))*dt(k))
      n2oexpd(k, 1) = -(6.31582e-2*xcd*EXP(-(xc*6.31582e-2)))
      n2oexp(k, 1) = EXP(-(xc*6.31582e-2))
      xcd = (n2oexpd(k, 1)*n2oexp(k, 1)+n2oexp(k, 1)*n2oexpd(k, 1))*&
&       n2oexp(k, 1) + n2oexp(k, 1)**2*n2oexpd(k, 1)
      xc = n2oexp(k, 1)*n2oexp(k, 1)*n2oexp(k, 1)
      xc1d = xcd*xc + xc*xcd
      xc1 = xc*xc
      xc2d = xc1d*xc1 + xc1*xc1d
      xc2 = xc1*xc1
      n2oexpd(k, 2) = (xcd*xc1+xc*xc1d)*xc2 + xc*xc1*xc2d
      n2oexp(k, 2) = xc*xc1*xc2
!-----four exponential by powers of 8 for band 7
    ELSE
      xcd = dn2o(k)*(pa(k)/500.0)**0.48*(7.4838e-6*dtd(k)*dt(k)+(&
&       1.3804e-3+7.4838e-6*dt(k))*dtd(k))
      xc = dn2o(k)*(pa(k)/500.0)**0.48*(1.+(1.3804e-3+7.4838e-6*dt(k))*&
&       dt(k))
      n2oexpd(k, 1) = -(5.35779e-2*xcd*EXP(-(xc*5.35779e-2)))
      n2oexp(k, 1) = EXP(-(xc*5.35779e-2))
      xcd = n2oexpd(k, 1)*n2oexp(k, 1) + n2oexp(k, 1)*n2oexpd(k, 1)
      xc = n2oexp(k, 1)*n2oexp(k, 1)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      n2oexpd(k, 2) = xcd*xc + xc*xcd
      n2oexp(k, 2) = xc*xc
      xcd = n2oexpd(k, 2)*n2oexp(k, 2) + n2oexp(k, 2)*n2oexpd(k, 2)
      xc = n2oexp(k, 2)*n2oexp(k, 2)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      n2oexpd(k, 3) = xcd*xc + xc*xcd
      n2oexp(k, 3) = xc*xc
      xcd = n2oexpd(k, 3)*n2oexp(k, 3) + n2oexp(k, 3)*n2oexpd(k, 3)
      xc = n2oexp(k, 3)*n2oexp(k, 3)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      n2oexpd(k, 4) = xcd*xc + xc*xcd
      n2oexp(k, 4) = xc*xc
    END IF
  END DO
END SUBROUTINE N2OEXPS_D

!  Differentiation of ch4exps in forward (tangent) mode:
!   variations   of useful results: ch4exp
!   with respect to varying inputs: dt ch4exp
!**********************************************************************
SUBROUTINE CH4EXPS_D(ib, np, dch4, pa, dt, dtd, ch4exp, ch4expd)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: dch4(0:np), pa(0:np), dt(0:np)
  REAL*8 :: dtd(0:np)
!---- output parameters -----
  REAL*8 :: ch4exp(0:np, 4)
  REAL*8 :: ch4expd(0:np, 4)
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcd
  INTRINSIC EXP
!*****  Scaling and absorption data are given in Table 5  *****
  DO k=0,np
!-----four exponentials for band 6
    IF (ib .EQ. 6) THEN
      xcd = dch4(k)*(1.5826e-4*dtd(k)*dt(k)+(1.7007e-2+1.5826e-4*dt(k))*&
&       dtd(k))
      xc = dch4(k)*(1.+(1.7007e-2+1.5826e-4*dt(k))*dt(k))
      ch4expd(k, 1) = -(5.80708e-3*xcd*EXP(-(xc*5.80708e-3)))
      ch4exp(k, 1) = EXP(-(xc*5.80708e-3))
!-----four exponentials by powers of 12 for band 7
    ELSE
      xcd = dch4(k)*(pa(k)/500.0)**0.65*((5.9590e-4-2.2931e-6*dt(k))*dtd&
&       (k)-2.2931e-6*dtd(k)*dt(k))
      xc = dch4(k)*(pa(k)/500.0)**0.65*(1.+(5.9590e-4-2.2931e-6*dt(k))*&
&       dt(k))
      ch4expd(k, 1) = -(6.29247e-2*xcd*EXP(-(xc*6.29247e-2)))
      ch4exp(k, 1) = EXP(-(xc*6.29247e-2))
      xcd = (ch4expd(k, 1)*ch4exp(k, 1)+ch4exp(k, 1)*ch4expd(k, 1))*&
&       ch4exp(k, 1) + ch4exp(k, 1)**2*ch4expd(k, 1)
      xc = ch4exp(k, 1)*ch4exp(k, 1)*ch4exp(k, 1)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      ch4expd(k, 2) = xcd*xc + xc*xcd
      ch4exp(k, 2) = xc*xc
      xcd = (ch4expd(k, 2)*ch4exp(k, 2)+ch4exp(k, 2)*ch4expd(k, 2))*&
&       ch4exp(k, 2) + ch4exp(k, 2)**2*ch4expd(k, 2)
      xc = ch4exp(k, 2)*ch4exp(k, 2)*ch4exp(k, 2)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      ch4expd(k, 3) = xcd*xc + xc*xcd
      ch4exp(k, 3) = xc*xc
      xcd = (ch4expd(k, 3)*ch4exp(k, 3)+ch4exp(k, 3)*ch4expd(k, 3))*&
&       ch4exp(k, 3) + ch4exp(k, 3)**2*ch4expd(k, 3)
      xc = ch4exp(k, 3)*ch4exp(k, 3)*ch4exp(k, 3)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      ch4expd(k, 4) = xcd*xc + xc*xcd
      ch4exp(k, 4) = xc*xc
    END IF
  END DO
END SUBROUTINE CH4EXPS_D

!  Differentiation of comexps in forward (tangent) mode:
!   variations   of useful results: comexp
!   with respect to varying inputs: dt comexp
!**********************************************************************
SUBROUTINE COMEXPS_D(ib, np, dcom, dt, dtd, comexp, comexpd)
  IMPLICIT NONE
  INTEGER :: ib, ik, np, k
!---- input parameters -----
  REAL*8 :: dcom(0:np), dt(0:np)
  REAL*8 :: dtd(0:np)
!---- output parameters -----
  REAL*8 :: comexp(0:np, 6)
  REAL*8 :: comexpd(0:np, 6)
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcd
  INTRINSIC EXP
  xcd = 0.0_8
!*****  Scaling and absorpton data are given in Table 6  *****
  DO k=0,np
    IF (ib .EQ. 4) THEN
      xcd = dcom(k)*(4.0447e-4*dtd(k)*dt(k)+(3.5775e-2+4.0447e-4*dt(k))*&
&       dtd(k))
      xc = dcom(k)*(1.+(3.5775e-2+4.0447e-4*dt(k))*dt(k))
    END IF
    IF (ib .EQ. 5) THEN
      xcd = dcom(k)*(3.7401e-4*dtd(k)*dt(k)+(3.4268e-2+3.7401e-4*dt(k))*&
&       dtd(k))
      xc = dcom(k)*(1.+(3.4268e-2+3.7401e-4*dt(k))*dt(k))
    END IF
    comexpd(k, 1) = -(1.922e-7*xcd*EXP(-(xc*1.922e-7)))
    comexp(k, 1) = EXP(-(xc*1.922e-7))
    DO ik=2,6
      xcd = comexpd(k, ik-1)*comexp(k, ik-1) + comexp(k, ik-1)*comexpd(k&
&       , ik-1)
      xc = comexp(k, ik-1)*comexp(k, ik-1)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      comexpd(k, ik) = xcd*comexp(k, ik-1) + xc*comexpd(k, ik-1)
      comexp(k, ik) = xc*comexp(k, ik-1)
    END DO
  END DO
END SUBROUTINE COMEXPS_D

!  Differentiation of cfcexps in forward (tangent) mode:
!   variations   of useful results: cfcexp
!   with respect to varying inputs: dt cfcexp
!**********************************************************************
SUBROUTINE CFCEXPS_D(ib, np, a1, b1, fk1, a2, b2, fk2, dcfc, dt, dtd, &
& cfcexp, cfcexpd)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: dcfc(0:np), dt(0:np)
  REAL*8 :: dtd(0:np)
!---- output parameters -----
  REAL*8 :: cfcexp(0:np)
  REAL*8 :: cfcexpd(0:np)
!---- static data -----
  REAL*8 :: a1, b1, fk1, a2, b2, fk2
!---- temporary arrays -----
  REAL*8 :: xf
  REAL*8 :: xfd
  INTRINSIC EXP
!**********************************************************************
  DO k=0,np
!-----compute the scaled cfc amount (xf) and exponential (cfcexp)
    IF (ib .EQ. 4) THEN
      xfd = dcfc(k)*(b1*dtd(k)*dt(k)+(a1+b1*dt(k))*dtd(k))
      xf = dcfc(k)*(1.+(a1+b1*dt(k))*dt(k))
      cfcexpd(k) = -(fk1*xfd*EXP(-(xf*fk1)))
      cfcexp(k) = EXP(-(xf*fk1))
    ELSE
      xfd = dcfc(k)*(b2*dtd(k)*dt(k)+(a2+b2*dt(k))*dtd(k))
      xf = dcfc(k)*(1.+(a2+b2*dt(k))*dt(k))
      cfcexpd(k) = -(fk2*xfd*EXP(-(xf*fk2)))
      cfcexp(k) = EXP(-(xf*fk2))
    END IF
  END DO
END SUBROUTINE CFCEXPS_D

!  Differentiation of b10exps in forward (tangent) mode:
!   variations   of useful results: co2exp h2oexp n2oexp conexp
!   with respect to varying inputs: dh2o dcont dt co2exp h2oexp
!                n2oexp conexp
!**********************************************************************
SUBROUTINE B10EXPS_D(np, dh2o, dh2od, dcont, dcontd, dco2, dn2o, pa, dt&
& , dtd, h2oexp, h2oexpd, conexp, conexpd, co2exp, co2expd, n2oexp, &
& n2oexpd)
  IMPLICIT NONE
  INTEGER :: np, k
!---- input parameters -----
  REAL*8 :: dh2o(0:np), dcont(0:np), dn2o(0:np)
  REAL*8 :: dh2od(0:np), dcontd(0:np)
  REAL*8 :: dco2(0:np), pa(0:np), dt(0:np)
  REAL*8 :: dtd(0:np)
!---- output parameters -----
  REAL*8 :: h2oexp(0:np, 5), conexp(0:np), co2exp(0:np, 6), n2oexp(0:np&
& , 2)
  REAL*8 :: h2oexpd(0:np, 5), conexpd(0:np), co2expd(0:np, 6), n2oexpd(0&
& :np, 2)
!---- temporary arrays -----
  REAL*8 :: xx, xx1, xx2, xx3
  REAL*8 :: xxd, xx1d, xx2d, xx3d
  INTRINSIC EXP
!**********************************************************************
  DO k=0,np
!-----Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3).
    xxd = pa(k)*(dh2od(k)*(1.+(0.0149+6.20e-5*dt(k))*dt(k))+dh2o(k)*(&
&     6.20e-5*dtd(k)*dt(k)+(0.0149+6.20e-5*dt(k))*dtd(k)))/500.0
    xx = dh2o(k)*(pa(k)/500.0)*(1.+(0.0149+6.20e-5*dt(k))*dt(k))
!-----six exponentials by powers of 8
    h2oexpd(k, 1) = -(0.10624*xxd*EXP(-(xx*0.10624)))
    h2oexp(k, 1) = EXP(-(xx*0.10624))
    xxd = h2oexpd(k, 1)*h2oexp(k, 1) + h2oexp(k, 1)*h2oexpd(k, 1)
    xx = h2oexp(k, 1)*h2oexp(k, 1)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    h2oexpd(k, 2) = xxd*xx + xx*xxd
    h2oexp(k, 2) = xx*xx
    xxd = h2oexpd(k, 2)*h2oexp(k, 2) + h2oexp(k, 2)*h2oexpd(k, 2)
    xx = h2oexp(k, 2)*h2oexp(k, 2)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    h2oexpd(k, 3) = xxd*xx + xx*xxd
    h2oexp(k, 3) = xx*xx
    xxd = h2oexpd(k, 3)*h2oexp(k, 3) + h2oexp(k, 3)*h2oexpd(k, 3)
    xx = h2oexp(k, 3)*h2oexp(k, 3)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    h2oexpd(k, 4) = xxd*xx + xx*xxd
    h2oexp(k, 4) = xx*xx
    xxd = h2oexpd(k, 4)*h2oexp(k, 4) + h2oexp(k, 4)*h2oexpd(k, 4)
    xx = h2oexp(k, 4)*h2oexp(k, 4)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    h2oexpd(k, 5) = xxd*xx + xx*xxd
    h2oexp(k, 5) = xx*xx
!-----one exponential of h2o continuum for sub-band 3a (Table 9).
    conexpd(k) = -(109.0*dcontd(k)*EXP(-(dcont(k)*109.0)))
    conexp(k) = EXP(-(dcont(k)*109.0))
!-----Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).
    xxd = dco2(k)*(pa(k)/300.0)**0.5*(1.02e-4*dtd(k)*dt(k)+(0.0179+&
&     1.02e-4*dt(k))*dtd(k))
    xx = dco2(k)*(pa(k)/300.0)**0.5*(1.+(0.0179+1.02e-4*dt(k))*dt(k))
!-----six exponentials by powers of 8
    co2expd(k, 1) = -(2.656e-5*xxd*EXP(-(xx*2.656e-5)))
    co2exp(k, 1) = EXP(-(xx*2.656e-5))
    xxd = co2expd(k, 1)*co2exp(k, 1) + co2exp(k, 1)*co2expd(k, 1)
    xx = co2exp(k, 1)*co2exp(k, 1)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 2) = xxd*xx + xx*xxd
    co2exp(k, 2) = xx*xx
    xxd = co2expd(k, 2)*co2exp(k, 2) + co2exp(k, 2)*co2expd(k, 2)
    xx = co2exp(k, 2)*co2exp(k, 2)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 3) = xxd*xx + xx*xxd
    co2exp(k, 3) = xx*xx
    xxd = co2expd(k, 3)*co2exp(k, 3) + co2exp(k, 3)*co2expd(k, 3)
    xx = co2exp(k, 3)*co2exp(k, 3)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 4) = xxd*xx + xx*xxd
    co2exp(k, 4) = xx*xx
    xxd = co2expd(k, 4)*co2exp(k, 4) + co2exp(k, 4)*co2expd(k, 4)
    xx = co2exp(k, 4)*co2exp(k, 4)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 5) = xxd*xx + xx*xxd
    co2exp(k, 5) = xx*xx
    xxd = co2expd(k, 5)*co2exp(k, 5) + co2exp(k, 5)*co2expd(k, 5)
    xx = co2exp(k, 5)*co2exp(k, 5)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 6) = xxd*xx + xx*xxd
    co2exp(k, 6) = xx*xx
!-----Compute the scaled n2o amount for Band 10 (Table 5).
    xxd = dn2o(k)*(3.6656e-6*dtd(k)*dt(k)+(1.4476e-3+3.6656e-6*dt(k))*&
&     dtd(k))
    xx = dn2o(k)*(1.+(1.4476e-3+3.6656e-6*dt(k))*dt(k))
!-----Two exponentials by powers of 58
    n2oexpd(k, 1) = -(0.25238*xxd*EXP(-(xx*0.25238)))
    n2oexp(k, 1) = EXP(-(xx*0.25238))
    xxd = n2oexpd(k, 1)*n2oexp(k, 1) + n2oexp(k, 1)*n2oexpd(k, 1)
    xx = n2oexp(k, 1)*n2oexp(k, 1)
    xx1d = xxd*xx + xx*xxd
    xx1 = xx*xx
    xx1d = xx1d*xx1 + xx1*xx1d
    xx1 = xx1*xx1
    xx2d = xx1d*xx1 + xx1*xx1d
    xx2 = xx1*xx1
    xx3d = xx2d*xx2 + xx2*xx2d
    xx3 = xx2*xx2
    n2oexpd(k, 2) = (xxd*xx1+xx*xx1d)*xx2*xx3 + xx*xx1*(xx2d*xx3+xx2*&
&     xx3d)
    n2oexp(k, 2) = xx*xx1*xx2*xx3
  END DO
END SUBROUTINE B10EXPS_D

!  Differentiation of tablup in forward (tangent) mode:
!   variations   of useful results: s1 s2 s3 tran
!   with respect to varying inputs: s1 s2 s3 dt dw tran
!**********************************************************************
SUBROUTINE TABLUP_D(nx1, nh1, dw, dwd, p, dt, dtd, s1, s1d, s2, s2d, s3&
& , s3d, w1, p1, dwe, dpe, coef1, coef2, coef3, tran, trand)
  IMPLICIT NONE
  INTEGER :: nx1, nh1
!---- input parameters -----
  REAL*8 :: w1, p1, dwe, dpe
  REAL*8 :: dw, p, dt
  REAL*8 :: dwd, dtd
  REAL*8 :: coef1(nx1, nh1), coef2(nx1, nh1), coef3(nx1, nh1)
!---- update parameter -----
  REAL*8 :: s1, s2, s3, tran
  REAL*8 :: s1d, s2d, s3d, trand
!---- temporary variables -----
  REAL*8 :: we, pe, fw, fp, pa, pb, pc, ax, ba, bb, t1, ca, cb, t2
  REAL*8 :: wed, ped, fwd, fpd, pad, pbd, pcd, axd, bad, bbd, t1d, cad, &
& cbd, t2d
  REAL*8 :: x1, x2, x3, xx, x1c
  REAL*8 :: x1d, x2d, x3d, xxd, x1cd
  INTEGER :: iw, ip
  INTRINSIC LOG10
  INTRINSIC REAL
  INTRINSIC MIN
  INTRINSIC INT
  INTRINSIC MAX
  REAL*8 :: y2
  REAL*8 :: y1
!-----Compute effective pressure (x2) and temperature (x3) following 
!     Eqs. (8.28) and (8.29)
  s1d = s1d + dwd
  s1 = s1 + dw
  s2d = s2d + p*dwd
  s2 = s2 + p*dw
  s3d = s3d + dtd*dw + dt*dwd
  s3 = s3 + dt*dw
  x1d = s1d
  x1 = s1
  x1cd = -(s1d/s1**2)
  x1c = 1.0/s1
  x2d = s2d*x1c + s2*x1cd
  x2 = s2*x1c
  x3d = s3d*x1c + s3*x1cd
  x3 = s3*x1c
!-----normalize we and pe
!       we=(log10(x1)-w1)/dwe
!       pe=(log10(x2)-p1)/dpe
  wed = dwe*x1d/(x1*LOG(10.0))
  we = (LOG10(x1)-w1)*dwe
  ped = dpe*x2d/(x2*LOG(10.0))
  pe = (LOG10(x2)-p1)*dpe
  y1 = REAL(nh1 - 1)
  IF (we .GT. y1) THEN
    we = y1
    wed = 0.0_8
  ELSE
    we = we
  END IF
  y2 = REAL(nx1 - 1)
  IF (pe .GT. y2) THEN
    pe = y2
    ped = 0.0_8
  ELSE
    pe = pe
  END IF
!-----assign iw and ip and compute the distance of we and pe 
!     from iw and ip.
  iw = INT(we + 1.0)
  IF (iw .GT. nh1 - 1) THEN
    iw = nh1 - 1
  ELSE
    iw = iw
  END IF
  IF (iw .LT. 2) THEN
    iw = 2
  ELSE
    iw = iw
  END IF
  fwd = wed
  fw = we - REAL(iw - 1)
  ip = INT(pe + 1.0)
  IF (ip .GT. nx1 - 1) THEN
    ip = nx1 - 1
  ELSE
    ip = ip
  END IF
  IF (ip .LT. 1) THEN
    ip = 1
  ELSE
    ip = ip
  END IF
  fpd = ped
  fp = pe - REAL(ip - 1)
!-----linear interpolation in pressure
  pad = (coef1(ip+1, iw-1)-coef1(ip, iw-1))*fpd
  pa = coef1(ip, iw-1) + (coef1(ip+1, iw-1)-coef1(ip, iw-1))*fp
  pbd = (coef1(ip+1, iw)-coef1(ip, iw))*fpd
  pb = coef1(ip, iw) + (coef1(ip+1, iw)-coef1(ip, iw))*fp
  pcd = (coef1(ip+1, iw+1)-coef1(ip, iw+1))*fpd
  pc = coef1(ip, iw+1) + (coef1(ip+1, iw+1)-coef1(ip, iw+1))*fp
!-----quadratic interpolation in absorber amount for coef1
  axd = 0.5*(((pcd+pad)*fw+(pc+pa)*fwd+pcd-pad)*fw+((pc+pa)*fw+(pc-pa))*&
&   fwd) + pbd*(1.-fw*fw) + pb*(-(fwd*fw)-fw*fwd)
  ax = ((pc+pa)*fw+(pc-pa))*fw*0.5 + pb*(1.-fw*fw)
!-----linear interpolation in absorber amount for coef2 and coef3
  bad = (coef2(ip+1, iw)-coef2(ip, iw))*fpd
  ba = coef2(ip, iw) + (coef2(ip+1, iw)-coef2(ip, iw))*fp
  bbd = (coef2(ip+1, iw+1)-coef2(ip, iw+1))*fpd
  bb = coef2(ip, iw+1) + (coef2(ip+1, iw+1)-coef2(ip, iw+1))*fp
  t1d = bad + (bbd-bad)*fw + (bb-ba)*fwd
  t1 = ba + (bb-ba)*fw
  cad = (coef3(ip+1, iw)-coef3(ip, iw))*fpd
  ca = coef3(ip, iw) + (coef3(ip+1, iw)-coef3(ip, iw))*fp
  cbd = (coef3(ip+1, iw+1)-coef3(ip, iw+1))*fpd
  cb = coef3(ip, iw+1) + (coef3(ip+1, iw+1)-coef3(ip, iw+1))*fp
  t2d = cad + (cbd-cad)*fw + (cb-ca)*fwd
  t2 = ca + (cb-ca)*fw
!-----update the total transmittance between levels k1 and k2
  xxd = axd + (t1d+t2d*x3+t2*x3d)*x3 + (t1+t2*x3)*x3d
  xx = ax + (t1+t2*x3)*x3
  IF (xx .GT. 0.9999999) THEN
    xx = 0.9999999
    xxd = 0.0_8
  ELSE
    xx = xx
  END IF
  IF (xx .LT. 0.0000001) THEN
    xx = 0.0000001
    xxd = 0.0_8
  ELSE
    xx = xx
  END IF
  trand = trand*xx + tran*xxd
  tran = tran*xx
END SUBROUTINE TABLUP_D

!  Differentiation of h2okdis in forward (tangent) mode:
!   variations   of useful results: tran tcon th2o
!   with respect to varying inputs: tcon h2oexp th2o conexp
!**********************************************************************
SUBROUTINE H2OKDIS_D(ib, np, k, fkw, gkw, ne, h2oexp, h2oexpd, conexp, &
& conexpd, th2o, th2od, tcon, tcond, tran, trand)
  IMPLICIT NONE
!---- input parameters ------
  INTEGER :: ib, ne, np, k
  REAL*8 :: h2oexp(0:np, 6), conexp(0:np, 3)
  REAL*8 :: h2oexpd(0:np, 6), conexpd(0:np, 3)
  REAL*8 :: fkw(6, 9), gkw(6, 3)
!---- updated parameters -----
  REAL*8 :: th2o(6), tcon(3), tran
  REAL*8 :: th2od(6), tcond(3), trand
!---- temporary arrays -----
  REAL*8 :: trnth2o
  REAL*8 :: trnth2od
!-----tco2 are the six exp factors between levels k1 and k2 
!     tran is the updated total transmittance between levels k1 and k2
!-----th2o is the 6 exp factors between levels k1 and k2 due to
!     h2o line absorption. 
!-----tcon is the 3 exp factors between levels k1 and k2 due to
!     h2o continuum absorption.
!-----trnth2o is the total transmittance between levels k1 and k2 due
!     to both line and continuum absorption.
!-----Compute th2o following Eq. (8.23).
  th2od(1) = th2od(1)*h2oexp(k, 1) + th2o(1)*h2oexpd(k, 1)
  th2o(1) = th2o(1)*h2oexp(k, 1)
  th2od(2) = th2od(2)*h2oexp(k, 2) + th2o(2)*h2oexpd(k, 2)
  th2o(2) = th2o(2)*h2oexp(k, 2)
  th2od(3) = th2od(3)*h2oexp(k, 3) + th2o(3)*h2oexpd(k, 3)
  th2o(3) = th2o(3)*h2oexp(k, 3)
  th2od(4) = th2od(4)*h2oexp(k, 4) + th2o(4)*h2oexpd(k, 4)
  th2o(4) = th2o(4)*h2oexp(k, 4)
  th2od(5) = th2od(5)*h2oexp(k, 5) + th2o(5)*h2oexpd(k, 5)
  th2o(5) = th2o(5)*h2oexp(k, 5)
  th2od(6) = th2od(6)*h2oexp(k, 6) + th2o(6)*h2oexpd(k, 6)
  th2o(6) = th2o(6)*h2oexp(k, 6)
  IF (ne .EQ. 0) THEN
!-----Compute trnh2o following Eq. (8.25). fkw is given in Table 4.
    trnth2od = fkw(1, ib)*th2od(1) + fkw(2, ib)*th2od(2) + fkw(3, ib)*&
&     th2od(3) + fkw(4, ib)*th2od(4) + fkw(5, ib)*th2od(5) + fkw(6, ib)*&
&     th2od(6)
    trnth2o = fkw(1, ib)*th2o(1) + fkw(2, ib)*th2o(2) + fkw(3, ib)*th2o(&
&     3) + fkw(4, ib)*th2o(4) + fkw(5, ib)*th2o(5) + fkw(6, ib)*th2o(6)
    trand = tran*trnth2od
    tran = tran*trnth2o
  ELSE IF (ne .EQ. 1) THEN
!-----Compute trnh2o following Eqs. (8.25) and (4.27).
    tcond(1) = tcond(1)*conexp(k, 1) + tcon(1)*conexpd(k, 1)
    tcon(1) = tcon(1)*conexp(k, 1)
    trnth2od = (fkw(1, ib)*th2od(1)+fkw(2, ib)*th2od(2)+fkw(3, ib)*th2od&
&     (3)+fkw(4, ib)*th2od(4)+fkw(5, ib)*th2od(5)+fkw(6, ib)*th2od(6))*&
&     tcon(1) + (fkw(1, ib)*th2o(1)+fkw(2, ib)*th2o(2)+fkw(3, ib)*th2o(3&
&     )+fkw(4, ib)*th2o(4)+fkw(5, ib)*th2o(5)+fkw(6, ib)*th2o(6))*tcond(&
&     1)
    trnth2o = (fkw(1, ib)*th2o(1)+fkw(2, ib)*th2o(2)+fkw(3, ib)*th2o(3)+&
&     fkw(4, ib)*th2o(4)+fkw(5, ib)*th2o(5)+fkw(6, ib)*th2o(6))*tcon(1)
    trand = tran*trnth2od
    tran = tran*trnth2o
  ELSE
!-----For band 3. This band is divided into 3 subbands.
    tcond(1) = tcond(1)*conexp(k, 1) + tcon(1)*conexpd(k, 1)
    tcon(1) = tcon(1)*conexp(k, 1)
    tcond(2) = tcond(2)*conexp(k, 2) + tcon(2)*conexpd(k, 2)
    tcon(2) = tcon(2)*conexp(k, 2)
    tcond(3) = tcond(3)*conexp(k, 3) + tcon(3)*conexpd(k, 3)
    tcon(3) = tcon(3)*conexp(k, 3)
!-----Compute trnh2o following Eqs. (4.29) and (8.25).
    trnth2od = (gkw(1, 1)*th2od(1)+gkw(2, 1)*th2od(2)+gkw(3, 1)*th2od(3)&
&     +gkw(4, 1)*th2od(4)+gkw(5, 1)*th2od(5)+gkw(6, 1)*th2od(6))*tcon(1)&
&     + (gkw(1, 1)*th2o(1)+gkw(2, 1)*th2o(2)+gkw(3, 1)*th2o(3)+gkw(4, 1)&
&     *th2o(4)+gkw(5, 1)*th2o(5)+gkw(6, 1)*th2o(6))*tcond(1) + (gkw(1, 2&
&     )*th2od(1)+gkw(2, 2)*th2od(2)+gkw(3, 2)*th2od(3)+gkw(4, 2)*th2od(4&
&     )+gkw(5, 2)*th2od(5)+gkw(6, 2)*th2od(6))*tcon(2) + (gkw(1, 2)*th2o&
&     (1)+gkw(2, 2)*th2o(2)+gkw(3, 2)*th2o(3)+gkw(4, 2)*th2o(4)+gkw(5, 2&
&     )*th2o(5)+gkw(6, 2)*th2o(6))*tcond(2) + (gkw(1, 3)*th2od(1)+gkw(2&
&     , 3)*th2od(2)+gkw(3, 3)*th2od(3)+gkw(4, 3)*th2od(4)+gkw(5, 3)*&
&     th2od(5)+gkw(6, 3)*th2od(6))*tcon(3) + (gkw(1, 3)*th2o(1)+gkw(2, 3&
&     )*th2o(2)+gkw(3, 3)*th2o(3)+gkw(4, 3)*th2o(4)+gkw(5, 3)*th2o(5)+&
&     gkw(6, 3)*th2o(6))*tcond(3)
    trnth2o = (gkw(1, 1)*th2o(1)+gkw(2, 1)*th2o(2)+gkw(3, 1)*th2o(3)+gkw&
&     (4, 1)*th2o(4)+gkw(5, 1)*th2o(5)+gkw(6, 1)*th2o(6))*tcon(1) + (gkw&
&     (1, 2)*th2o(1)+gkw(2, 2)*th2o(2)+gkw(3, 2)*th2o(3)+gkw(4, 2)*th2o(&
&     4)+gkw(5, 2)*th2o(5)+gkw(6, 2)*th2o(6))*tcon(2) + (gkw(1, 3)*th2o(&
&     1)+gkw(2, 3)*th2o(2)+gkw(3, 3)*th2o(3)+gkw(4, 3)*th2o(4)+gkw(5, 3)&
&     *th2o(5)+gkw(6, 3)*th2o(6))*tcon(3)
    trand = tran*trnth2od
    tran = tran*trnth2o
  END IF
END SUBROUTINE H2OKDIS_D

!  Differentiation of n2okdis in forward (tangent) mode:
!   variations   of useful results: tran tn2o
!   with respect to varying inputs: tran tn2o n2oexp
!**********************************************************************
SUBROUTINE N2OKDIS_D(ib, np, k, n2oexp, n2oexpd, tn2o, tn2od, tran, &
& trand)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: n2oexp(0:np, 4)
  REAL*8 :: n2oexpd(0:np, 4)
!---- updated parameters -----
  REAL*8 :: tn2o(4), tran
  REAL*8 :: tn2od(4), trand
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcd
!-----tn2o is computed from Eq. (8.23). 
!     xc is the total n2o transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.
!-----band 6
  IF (ib .EQ. 6) THEN
    tn2od(1) = tn2od(1)*n2oexp(k, 1) + tn2o(1)*n2oexpd(k, 1)
    tn2o(1) = tn2o(1)*n2oexp(k, 1)
    xcd = 0.940414*tn2od(1)
    xc = 0.940414*tn2o(1)
    tn2od(2) = tn2od(2)*n2oexp(k, 2) + tn2o(2)*n2oexpd(k, 2)
    tn2o(2) = tn2o(2)*n2oexp(k, 2)
    xcd = xcd + 0.059586*tn2od(2)
    xc = xc + 0.059586*tn2o(2)
!-----band 7
  ELSE
    tn2od(1) = tn2od(1)*n2oexp(k, 1) + tn2o(1)*n2oexpd(k, 1)
    tn2o(1) = tn2o(1)*n2oexp(k, 1)
    xcd = 0.561961*tn2od(1)
    xc = 0.561961*tn2o(1)
    tn2od(2) = tn2od(2)*n2oexp(k, 2) + tn2o(2)*n2oexpd(k, 2)
    tn2o(2) = tn2o(2)*n2oexp(k, 2)
    xcd = xcd + 0.138707*tn2od(2)
    xc = xc + 0.138707*tn2o(2)
    tn2od(3) = tn2od(3)*n2oexp(k, 3) + tn2o(3)*n2oexpd(k, 3)
    tn2o(3) = tn2o(3)*n2oexp(k, 3)
    xcd = xcd + 0.240670*tn2od(3)
    xc = xc + 0.240670*tn2o(3)
    tn2od(4) = tn2od(4)*n2oexp(k, 4) + tn2o(4)*n2oexpd(k, 4)
    tn2o(4) = tn2o(4)*n2oexp(k, 4)
    xcd = xcd + 0.058662*tn2od(4)
    xc = xc + 0.058662*tn2o(4)
  END IF
  trand = trand*xc + tran*xcd
  tran = tran*xc
END SUBROUTINE N2OKDIS_D

!  Differentiation of ch4kdis in forward (tangent) mode:
!   variations   of useful results: tran tch4
!   with respect to varying inputs: tran ch4exp tch4
!**********************************************************************
SUBROUTINE CH4KDIS_D(ib, np, k, ch4exp, ch4expd, tch4, tch4d, tran, &
& trand)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: ch4exp(0:np, 4)
  REAL*8 :: ch4expd(0:np, 4)
!---- updated parameters -----
  REAL*8 :: tch4(4), tran
  REAL*8 :: tch4d(4), trand
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcd
!-----tch4 is computed from Eq. (8.23). 
!     xc is the total ch4 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.
!-----band 6
  IF (ib .EQ. 6) THEN
    tch4d(1) = tch4d(1)*ch4exp(k, 1) + tch4(1)*ch4expd(k, 1)
    tch4(1) = tch4(1)*ch4exp(k, 1)
    xcd = tch4d(1)
    xc = tch4(1)
!-----band 7
  ELSE
    tch4d(1) = tch4d(1)*ch4exp(k, 1) + tch4(1)*ch4expd(k, 1)
    tch4(1) = tch4(1)*ch4exp(k, 1)
    xcd = 0.610650*tch4d(1)
    xc = 0.610650*tch4(1)
    tch4d(2) = tch4d(2)*ch4exp(k, 2) + tch4(2)*ch4expd(k, 2)
    tch4(2) = tch4(2)*ch4exp(k, 2)
    xcd = xcd + 0.280212*tch4d(2)
    xc = xc + 0.280212*tch4(2)
    tch4d(3) = tch4d(3)*ch4exp(k, 3) + tch4(3)*ch4expd(k, 3)
    tch4(3) = tch4(3)*ch4exp(k, 3)
    xcd = xcd + 0.107349*tch4d(3)
    xc = xc + 0.107349*tch4(3)
    tch4d(4) = tch4d(4)*ch4exp(k, 4) + tch4(4)*ch4expd(k, 4)
    tch4(4) = tch4(4)*ch4exp(k, 4)
    xcd = xcd + 0.001789*tch4d(4)
    xc = xc + 0.001789*tch4(4)
  END IF
  trand = trand*xc + tran*xcd
  tran = tran*xc
END SUBROUTINE CH4KDIS_D

!  Differentiation of comkdis in forward (tangent) mode:
!   variations   of useful results: tran tcom
!   with respect to varying inputs: tran tcom comexp
!**********************************************************************
SUBROUTINE COMKDIS_D(ib, np, k, comexp, comexpd, tcom, tcomd, tran, &
& trand)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: comexp(0:np, 6)
  REAL*8 :: comexpd(0:np, 6)
!---- updated parameters -----
  REAL*8 :: tcom(6), tran
  REAL*8 :: tcomd(6), trand
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcd
!-----tcom is computed from Eq. (8.23). 
!     xc is the total co2 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 6.
!-----band 4
  IF (ib .EQ. 4) THEN
    tcomd(1) = tcomd(1)*comexp(k, 1) + tcom(1)*comexpd(k, 1)
    tcom(1) = tcom(1)*comexp(k, 1)
    xcd = 0.12159*tcomd(1)
    xc = 0.12159*tcom(1)
    tcomd(2) = tcomd(2)*comexp(k, 2) + tcom(2)*comexpd(k, 2)
    tcom(2) = tcom(2)*comexp(k, 2)
    xcd = xcd + 0.24359*tcomd(2)
    xc = xc + 0.24359*tcom(2)
    tcomd(3) = tcomd(3)*comexp(k, 3) + tcom(3)*comexpd(k, 3)
    tcom(3) = tcom(3)*comexp(k, 3)
    xcd = xcd + 0.24981*tcomd(3)
    xc = xc + 0.24981*tcom(3)
    tcomd(4) = tcomd(4)*comexp(k, 4) + tcom(4)*comexpd(k, 4)
    tcom(4) = tcom(4)*comexp(k, 4)
    xcd = xcd + 0.26427*tcomd(4)
    xc = xc + 0.26427*tcom(4)
    tcomd(5) = tcomd(5)*comexp(k, 5) + tcom(5)*comexpd(k, 5)
    tcom(5) = tcom(5)*comexp(k, 5)
    xcd = xcd + 0.07807*tcomd(5)
    xc = xc + 0.07807*tcom(5)
    tcomd(6) = tcomd(6)*comexp(k, 6) + tcom(6)*comexpd(k, 6)
    tcom(6) = tcom(6)*comexp(k, 6)
    xcd = xcd + 0.04267*tcomd(6)
    xc = xc + 0.04267*tcom(6)
!-----band 5
  ELSE
    tcomd(1) = tcomd(1)*comexp(k, 1) + tcom(1)*comexpd(k, 1)
    tcom(1) = tcom(1)*comexp(k, 1)
    xcd = 0.06869*tcomd(1)
    xc = 0.06869*tcom(1)
    tcomd(2) = tcomd(2)*comexp(k, 2) + tcom(2)*comexpd(k, 2)
    tcom(2) = tcom(2)*comexp(k, 2)
    xcd = xcd + 0.14795*tcomd(2)
    xc = xc + 0.14795*tcom(2)
    tcomd(3) = tcomd(3)*comexp(k, 3) + tcom(3)*comexpd(k, 3)
    tcom(3) = tcom(3)*comexp(k, 3)
    xcd = xcd + 0.19512*tcomd(3)
    xc = xc + 0.19512*tcom(3)
    tcomd(4) = tcomd(4)*comexp(k, 4) + tcom(4)*comexpd(k, 4)
    tcom(4) = tcom(4)*comexp(k, 4)
    xcd = xcd + 0.33446*tcomd(4)
    xc = xc + 0.33446*tcom(4)
    tcomd(5) = tcomd(5)*comexp(k, 5) + tcom(5)*comexpd(k, 5)
    tcom(5) = tcom(5)*comexp(k, 5)
    xcd = xcd + 0.17199*tcomd(5)
    xc = xc + 0.17199*tcom(5)
    tcomd(6) = tcomd(6)*comexp(k, 6) + tcom(6)*comexpd(k, 6)
    tcom(6) = tcom(6)*comexp(k, 6)
    xcd = xcd + 0.08179*tcomd(6)
    xc = xc + 0.08179*tcom(6)
  END IF
  trand = trand*xc + tran*xcd
  tran = tran*xc
END SUBROUTINE COMKDIS_D

!  Differentiation of cfckdis in forward (tangent) mode:
!   variations   of useful results: tran tcfc
!   with respect to varying inputs: tran tcfc cfcexp
!**********************************************************************
SUBROUTINE CFCKDIS_D(np, k, cfcexp, cfcexpd, tcfc, tcfcd, tran, trand)
  IMPLICIT NONE
!---- input parameters -----
  INTEGER :: k, np
  REAL*8 :: cfcexp(0:np)
  REAL*8 :: cfcexpd(0:np)
!---- updated parameters -----
  REAL*8 :: tcfc, tran
  REAL*8 :: tcfcd, trand
!-----tcfc is the exp factors between levels k1 and k2. 
  tcfcd = tcfcd*cfcexp(k) + tcfc*cfcexpd(k)
  tcfc = tcfc*cfcexp(k)
  trand = trand*tcfc + tran*tcfcd
  tran = tran*tcfc
END SUBROUTINE CFCKDIS_D

!  Differentiation of b10kdis in forward (tangent) mode:
!   variations   of useful results: tran tn2o tcon th2o tco2
!   with respect to varying inputs: tn2o co2exp tcon h2oexp n2oexp
!                th2o tco2 conexp
!**********************************************************************
SUBROUTINE B10KDIS_D(np, k, h2oexp, h2oexpd, conexp, conexpd, co2exp, &
& co2expd, n2oexp, n2oexpd, th2o, th2od, tcon, tcond, tco2, tco2d, tn2o&
& , tn2od, tran, trand)
  IMPLICIT NONE
  INTEGER :: np, k
!---- input parameters -----
  REAL*8 :: h2oexp(0:np, 5), conexp(0:np), co2exp(0:np, 6), n2oexp(0:np&
& , 2)
  REAL*8 :: h2oexpd(0:np, 5), conexpd(0:np), co2expd(0:np, 6), n2oexpd(0&
& :np, 2)
!---- updated parameters -----
  REAL*8 :: th2o(6), tcon(3), tco2(6), tn2o(4), tran
  REAL*8 :: th2od(6), tcond(3), tco2d(6), tn2od(4), trand
!---- temporary arrays -----
  REAL*8 :: xx
  REAL*8 :: xxd
!-----For h2o line. The k-distribution functions are given in Table 4.
  th2od(1) = th2od(1)*h2oexp(k, 1) + th2o(1)*h2oexpd(k, 1)
  th2o(1) = th2o(1)*h2oexp(k, 1)
  xxd = 0.3153*th2od(1)
  xx = 0.3153*th2o(1)
  th2od(2) = th2od(2)*h2oexp(k, 2) + th2o(2)*h2oexpd(k, 2)
  th2o(2) = th2o(2)*h2oexp(k, 2)
  xxd = xxd + 0.4604*th2od(2)
  xx = xx + 0.4604*th2o(2)
  th2od(3) = th2od(3)*h2oexp(k, 3) + th2o(3)*h2oexpd(k, 3)
  th2o(3) = th2o(3)*h2oexp(k, 3)
  xxd = xxd + 0.1326*th2od(3)
  xx = xx + 0.1326*th2o(3)
  th2od(4) = th2od(4)*h2oexp(k, 4) + th2o(4)*h2oexpd(k, 4)
  th2o(4) = th2o(4)*h2oexp(k, 4)
  xxd = xxd + 0.0798*th2od(4)
  xx = xx + 0.0798*th2o(4)
  th2od(5) = th2od(5)*h2oexp(k, 5) + th2o(5)*h2oexpd(k, 5)
  th2o(5) = th2o(5)*h2oexp(k, 5)
  xxd = xxd + 0.0119*th2od(5)
  xx = xx + 0.0119*th2o(5)
  trand = xxd
  tran = xx
!-----For h2o continuum. Note that conexp(k,3) is for subband 3a.
  tcond(1) = tcond(1)*conexp(k) + tcon(1)*conexpd(k)
  tcon(1) = tcon(1)*conexp(k)
  trand = trand*tcon(1) + tran*tcond(1)
  tran = tran*tcon(1)
!-----For co2 (Table 6)
  tco2d(1) = tco2d(1)*co2exp(k, 1) + tco2(1)*co2expd(k, 1)
  tco2(1) = tco2(1)*co2exp(k, 1)
  xxd = 0.2673*tco2d(1)
  xx = 0.2673*tco2(1)
  tco2d(2) = tco2d(2)*co2exp(k, 2) + tco2(2)*co2expd(k, 2)
  tco2(2) = tco2(2)*co2exp(k, 2)
  xxd = xxd + 0.2201*tco2d(2)
  xx = xx + 0.2201*tco2(2)
  tco2d(3) = tco2d(3)*co2exp(k, 3) + tco2(3)*co2expd(k, 3)
  tco2(3) = tco2(3)*co2exp(k, 3)
  xxd = xxd + 0.2106*tco2d(3)
  xx = xx + 0.2106*tco2(3)
  tco2d(4) = tco2d(4)*co2exp(k, 4) + tco2(4)*co2expd(k, 4)
  tco2(4) = tco2(4)*co2exp(k, 4)
  xxd = xxd + 0.2409*tco2d(4)
  xx = xx + 0.2409*tco2(4)
  tco2d(5) = tco2d(5)*co2exp(k, 5) + tco2(5)*co2expd(k, 5)
  tco2(5) = tco2(5)*co2exp(k, 5)
  xxd = xxd + 0.0196*tco2d(5)
  xx = xx + 0.0196*tco2(5)
  tco2d(6) = tco2d(6)*co2exp(k, 6) + tco2(6)*co2expd(k, 6)
  tco2(6) = tco2(6)*co2exp(k, 6)
  xxd = xxd + 0.0415*tco2d(6)
  xx = xx + 0.0415*tco2(6)
  trand = trand*xx + tran*xxd
  tran = tran*xx
!-----For n2o (Table 5)
  tn2od(1) = tn2od(1)*n2oexp(k, 1) + tn2o(1)*n2oexpd(k, 1)
  tn2o(1) = tn2o(1)*n2oexp(k, 1)
  xxd = 0.970831*tn2od(1)
  xx = 0.970831*tn2o(1)
  tn2od(2) = tn2od(2)*n2oexp(k, 2) + tn2o(2)*n2oexpd(k, 2)
  tn2o(2) = tn2o(2)*n2oexp(k, 2)
  xxd = xxd + 0.029169*tn2od(2)
  xx = xx + 0.029169*tn2o(2)
  trand = trand*(xx-1.0) + tran*xxd
  tran = tran*(xx-1.0)
END SUBROUTINE B10KDIS_D

!  Differentiation of cldovlp in forward (tangent) mode:
!   variations   of useful results: cldhi cldlw cldmd
!   with respect to varying inputs: cldhi ett cldlw enn cldmd
!mjs
!***********************************************************************
SUBROUTINE CLDOVLP_D(np, k1, k2, ict, icb, icx, ncld, enn, ennd, ett, &
& ettd, cldhi, cldhid, cldmd, cldmdd, cldlw, cldlwd)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: np, k1, k2, ict, icb, icx(0:np), ncld(3)
  REAL*8, INTENT(IN) :: enn(0:np), ett(0:np)
  REAL*8, INTENT(IN) :: ennd(0:np), ettd(0:np)
  REAL*8, INTENT(INOUT) :: cldhi, cldmd, cldlw
  REAL*8, INTENT(INOUT) :: cldhid, cldmdd, cldlwd
  INTEGER :: j, k, km, kx
  km = k2 - 1
  IF (km .LT. ict) THEN
! do high clouds
    kx = ncld(1)
    IF (kx .EQ. 1 .OR. cldhi .EQ. 0.) THEN
      cldhid = ennd(km)
      cldhi = enn(km)
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      cldhi = 0.0
      IF (kx .NE. 0) THEN
        cldhid = 0.0_8
        DO k=ict-kx,ict-1
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            cldhid = ennd(j) + ettd(j)*cldhi + ett(j)*cldhid
            cldhi = enn(j) + ett(j)*cldhi
          END IF
        END DO
      ELSE
        cldhid = 0.0_8
      END IF
    END IF
  ELSE IF (km .GE. ict .AND. km .LT. icb) THEN
! do middle clouds
    kx = ncld(2)
    IF (kx .EQ. 1 .OR. cldmd .EQ. 0.) THEN
      cldmdd = ennd(km)
      cldmd = enn(km)
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      cldmd = 0.0
      IF (kx .NE. 0) THEN
        cldmdd = 0.0_8
        DO k=icb-kx,icb-1
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            cldmdd = ennd(j) + ettd(j)*cldmd + ett(j)*cldmdd
            cldmd = enn(j) + ett(j)*cldmd
          END IF
        END DO
      ELSE
        cldmdd = 0.0_8
      END IF
    END IF
  ELSE
! do low clouds
    kx = ncld(3)
    IF (kx .EQ. 1 .OR. cldlw .EQ. 0.) THEN
      cldlwd = ennd(km)
      cldlw = enn(km)
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      cldlw = 0.0
      IF (kx .NE. 0) THEN
        cldlwd = 0.0_8
        DO k=np+1-kx,np
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            cldlwd = ennd(j) + ettd(j)*cldlw + ett(j)*cldlwd
            cldlw = enn(j) + ett(j)*cldlw
          END IF
        END DO
      ELSE
        cldlwd = 0.0_8
      END IF
    END IF
  END IF
END SUBROUTINE CLDOVLP_D

!  Differentiation of getirtau1 in forward (tangent) mode:
!   variations   of useful results: tcldlyr enn
!   with respect to varying inputs: hydromets tcldlyr fcld enn
!                reff
SUBROUTINE GETIRTAU1_D(ib, nlevs, dp, fcld, fcldd, reff, reffd, &
& hydromets, hydrometsd, tcldlyr, tcldlyrd, enn, ennd, aib_ir1, awb_ir1&
& , aiw_ir1, aww_ir1, aig_ir1, awg_ir1, cons_grav)
  IMPLICIT NONE
! !INPUT PARAMETERS:
!  Band number
  INTEGER, INTENT(IN) :: ib
!  Number of levels
  INTEGER, INTENT(IN) :: nlevs
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
  REAL*8, INTENT(IN) :: aib_ir1(3, 10), awb_ir1(4, 10), aiw_ir1(4, 10)
  REAL*8, INTENT(IN) :: aww_ir1(4, 10), aig_ir1(4, 10), awg_ir1(4, 10)
  REAL*8, INTENT(IN) :: cons_grav
! !OUTPUT PARAMETERS:
!  Flux transmissivity?
  REAL*8, INTENT(OUT) :: tcldlyr(0:nlevs)
  REAL*8, INTENT(OUT) :: tcldlyrd(0:nlevs)
!  Flux transmissivity of a cloud layer?
  REAL*8, INTENT(OUT) :: enn(0:nlevs)
  REAL*8, INTENT(OUT) :: ennd(0:nlevs)
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for infrared wavelengths
!  Slots for reff, hydrometeors and tauall are as follows:
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
!    2011.11.18   MAT moved to Radiation_Shared and revised arg list, units
!
!EOP
!------------------------------------------------------------------------------
!BOC
  INTEGER :: k
  REAL*8 :: taucld1, taucld2, taucld3, taucld4
  REAL*8 :: taucld1d, taucld2d, taucld3d, taucld4d
  REAL*8 :: g1, g2, g3, g4, gg
  REAL*8 :: g1d, g2d, g3d, g4d, ggd
  REAL*8 :: w1, w2, w3, w4, ww
  REAL*8 :: w1d, w2d, w3d, w4d, wwd
  REAL*8 :: ff, tauc
  REAL*8 :: ffd, taucd
  REAL*8 :: reff_snow
  REAL*8 :: reff_snowd
  INTRINSIC MIN
  INTRINSIC ABS
  INTRINSIC MAX
  INTRINSIC EXP
  REAL*8 :: pwr1
  REAL*8 :: pwr1d
  REAL*8 :: max1d
  REAL*8 :: abs0
  REAL*8 :: max1
!-----Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!     Rain optical thickness is set to 0.00307 /(gm/m**2).
!     It is for a specific drop size distribution provided by Q. Fu.
  tcldlyrd(0) = 0.0_8
  tcldlyr(0) = 1.0
  ennd(0) = 0.0_8
  enn(0) = 0.0
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.0) THEN
      taucld1 = 0.0
      taucld1d = 0.0_8
    ELSE
      IF (reff(k, 1) .GT. 0.0 .OR. (reff(k, 1) .LT. 0.0 .AND. aib_ir1(3&
&         , ib) .EQ. INT(aib_ir1(3, ib)))) THEN
        pwr1d = aib_ir1(3, ib)*reff(k, 1)**(aib_ir1(3, ib)-1)*reffd(k, 1&
&         )
      ELSE IF (reff(k, 1) .EQ. 0.0 .AND. aib_ir1(3, ib) .EQ. 1.0) THEN
        pwr1d = reffd(k, 1)
      ELSE
        pwr1d = 0.0
      END IF
      pwr1 = reff(k, 1)**aib_ir1(3, ib)
      taucld1d = dp(k)*1.0e3*(hydrometsd(k, 1)*(aib_ir1(1, ib)+aib_ir1(2&
&       , ib)/pwr1)-hydromets(k, 1)*aib_ir1(2, ib)*pwr1d/pwr1**2)/&
&       cons_grav
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*(aib_ir1(1, ib)+&
&       aib_ir1(2, ib)/pwr1)
    END IF
    taucld2d = dp(k)*1.0e3*(hydrometsd(k, 2)*(awb_ir1(1, ib)+(awb_ir1(2&
&     , ib)+(awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reff(&
&     k, 2))+hydromets(k, 2)*((awb_ir1(4, ib)*reffd(k, 2)*reff(k, 2)+(&
&     awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reffd(k, 2))*reff(k, 2)+&
&     (awb_ir1(2, ib)+(awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reff(k&
&     , 2))*reffd(k, 2)))/cons_grav
    taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_ir1(1, ib)+(&
&     awb_ir1(2, ib)+(awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reff(k, &
&     2))*reff(k, 2))
    taucld3d = 0.00307*dp(k)*1.0e3*hydrometsd(k, 3)/cons_grav
    taucld3 = 0.00307*(dp(k)*1.0e3/cons_grav*hydromets(k, 3))
    IF (reff(k, 4) .GT. 112.0) THEN
      reff_snow = 112.0
      reff_snowd = 0.0_8
    ELSE
      reff_snowd = reffd(k, 4)
      reff_snow = reff(k, 4)
    END IF
    IF (reff_snow .LE. 0.0) THEN
      taucld4 = 0.0
      taucld4d = 0.0_8
    ELSE
      IF (reff_snow .GT. 0.0 .OR. (reff_snow .LT. 0.0 .AND. aib_ir1(3, &
&         ib) .EQ. INT(aib_ir1(3, ib)))) THEN
        pwr1d = aib_ir1(3, ib)*reff_snow**(aib_ir1(3, ib)-1)*reff_snowd
      ELSE IF (reff_snow .EQ. 0.0 .AND. aib_ir1(3, ib) .EQ. 1.0) THEN
        pwr1d = reff_snowd
      ELSE
        pwr1d = 0.0
      END IF
      pwr1 = reff_snow**aib_ir1(3, ib)
      taucld4d = dp(k)*1.0e3*(hydrometsd(k, 4)*(aib_ir1(1, ib)+aib_ir1(2&
&       , ib)/pwr1)-hydromets(k, 4)*aib_ir1(2, ib)*pwr1d/pwr1**2)/&
&       cons_grav
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*(aib_ir1(1, ib)+&
&       aib_ir1(2, ib)/pwr1)
    END IF
!-----Compute cloud single-scattering albedo and asymmetry factor for
!     a mixture of ice particles and liquid drops following 
!     Eqs. (6.5), (6.6), (6.15) and (6.16).
!     Single-scattering albedo and asymmetry factor of rain are set
!     to 0.54 and 0.95, respectively, based on the information provided
!     by Prof. Qiang Fu.
    taucd = taucld1d + taucld2d + taucld3d + taucld4d
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      w1d = taucld1d*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reff(k, 1)) + taucld1*((&
&       aiw_ir1(4, ib)*reffd(k, 1)*reff(k, 1)+(aiw_ir1(3, ib)+aiw_ir1(4&
&       , ib)*reff(k, 1))*reffd(k, 1))*reff(k, 1)+(aiw_ir1(2, ib)+(&
&       aiw_ir1(3, ib)+aiw_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reffd(k, 1&
&       ))
      w1 = taucld1*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reff(k, 1))
      w2d = taucld2d*(aww_ir1(1, ib)+(aww_ir1(2, ib)+(aww_ir1(3, ib)+&
&       aww_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reff(k, 2)) + taucld2*((&
&       aww_ir1(4, ib)*reffd(k, 2)*reff(k, 2)+(aww_ir1(3, ib)+aww_ir1(4&
&       , ib)*reff(k, 2))*reffd(k, 2))*reff(k, 2)+(aww_ir1(2, ib)+(&
&       aww_ir1(3, ib)+aww_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reffd(k, 2&
&       ))
      w2 = taucld2*(aww_ir1(1, ib)+(aww_ir1(2, ib)+(aww_ir1(3, ib)+&
&       aww_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reff(k, 2))
      w3d = 0.54*taucld3d
      w3 = taucld3*0.54
      w4d = taucld4d*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff_snow)*reff_snow)*reff_snow) + taucld4*((&
&       aiw_ir1(4, ib)*reff_snowd*reff_snow+(aiw_ir1(3, ib)+aiw_ir1(4, &
&       ib)*reff_snow)*reff_snowd)*reff_snow+(aiw_ir1(2, ib)+(aiw_ir1(3&
&       , ib)+aiw_ir1(4, ib)*reff_snow)*reff_snow)*reff_snowd)
      w4 = taucld4*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff_snow)*reff_snow)*reff_snow)
      wwd = ((w1d+w2d+w3d+w4d)*tauc-(w1+w2+w3+w4)*taucd)/tauc**2
      ww = (w1+w2+w3+w4)/tauc
      g1d = w1d*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(&
&       4, ib)*reff(k, 1))*reff(k, 1))*reff(k, 1)) + w1*((aig_ir1(4, ib)&
&       *reffd(k, 1)*reff(k, 1)+(aig_ir1(3, ib)+aig_ir1(4, ib)*reff(k, 1&
&       ))*reffd(k, 1))*reff(k, 1)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+&
&       aig_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reffd(k, 1))
      g1 = w1*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff(k, 1))*reff(k, 1))*reff(k, 1))
      g2d = w2d*(awg_ir1(1, ib)+(awg_ir1(2, ib)+(awg_ir1(3, ib)+awg_ir1(&
&       4, ib)*reff(k, 2))*reff(k, 2))*reff(k, 2)) + w2*((awg_ir1(4, ib)&
&       *reffd(k, 2)*reff(k, 2)+(awg_ir1(3, ib)+awg_ir1(4, ib)*reff(k, 2&
&       ))*reffd(k, 2))*reff(k, 2)+(awg_ir1(2, ib)+(awg_ir1(3, ib)+&
&       awg_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reffd(k, 2))
      g2 = w2*(awg_ir1(1, ib)+(awg_ir1(2, ib)+(awg_ir1(3, ib)+awg_ir1(4&
&       , ib)*reff(k, 2))*reff(k, 2))*reff(k, 2))
      g3d = 0.95*w3d
      g3 = w3*0.95
      g4d = w4d*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(&
&       4, ib)*reff_snow)*reff_snow)*reff_snow) + w4*((aig_ir1(4, ib)*&
&       reff_snowd*reff_snow+(aig_ir1(3, ib)+aig_ir1(4, ib)*reff_snow)*&
&       reff_snowd)*reff_snow+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff_snow)*reff_snow)*reff_snowd)
      g4 = w4*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff_snow)*reff_snow)*reff_snow)
      IF (w1 + w2 + w3 + w4 .GE. 0.) THEN
        abs0 = w1 + w2 + w3 + w4
      ELSE
        abs0 = -(w1+w2+w3+w4)
      END IF
      IF (abs0 .GT. 0.0) THEN
        ggd = ((g1d+g2d+g3d+g4d)*(w1+w2+w3+w4)-(g1+g2+g3+g4)*(w1d+w2d+&
&         w3d+w4d))/(w1+w2+w3+w4)**2
        gg = (g1+g2+g3+g4)/(w1+w2+w3+w4)
      ELSE
        gg = 0.5
        ggd = 0.0_8
      END IF
!-----Parameterization of LW scattering following Eqs. (6.11)
!     and (6.12). 
      ffd = (0.1185*ggd*gg+(0.0076+0.1185*gg)*ggd)*gg + (0.3739+(0.0076+&
&       0.1185*gg)*gg)*ggd
      ff = 0.5 + (0.3739+(0.0076+0.1185*gg)*gg)*gg
      IF (1. - ww*ff .LT. 0.0) THEN
        max1 = 0.0
        max1d = 0.0_8
      ELSE
        max1d = -(wwd*ff) - ww*ffd
        max1 = 1. - ww*ff
      END IF
!ALT: temporary protection against negative cloud optical thickness
      taucd = max1d*tauc + max1*taucd
      tauc = max1*tauc
!-----compute cloud diffuse transmittance. It is approximated by using 
!     a diffusivity factor of 1.66.
      tcldlyrd(k) = -(1.66*taucd*EXP(-(1.66*tauc)))
      tcldlyr(k) = EXP(-(1.66*tauc))
! N in the documentation (6.13)
      ennd(k) = fcldd(k)*(1.0-tcldlyr(k)) - fcld(k)*tcldlyrd(k)
      enn(k) = fcld(k)*(1.0-tcldlyr(k))
    ELSE
      tcldlyrd(k) = 0.0_8
      tcldlyr(k) = 1.0
      ennd(k) = 0.0_8
      enn(k) = 0.0
    END IF
  END DO
END SUBROUTINE GETIRTAU1_D

end module IRRAD_TL
