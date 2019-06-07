module IRRAD_AD

USE IRRADMOD

IMPLICIT NONE

PRIVATE
PUBLIC :: irrad_b

contains

SUBROUTINE IRRAD_B(np, ple_dev, ta_dev, ta_devb, wa_dev, wa_devb, oa_dev&
& , oa_devb, tb_dev, tb_devb, co2, trace, n2o_dev, ch4_dev, cfc11_dev, &
& cfc12_dev, cfc22_dev, cwc_dev, cwc_devb, fcld_dev, fcld_devb, ict, icb&
& , reff_dev, reff_devb, ns, fs_dev, tg_dev, eg_dev, tv_dev, ev_dev, &
& rv_dev, na, nb, taua_dev, taua_devb, ssaa_dev, ssaa_devb, asya_dev, &
& asya_devb, flxu_dev, flxu_devb, flxd_dev, flxd_devb, dfdts_dev, &
& dfdts_devb, aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, awg_ir, xkw, xke, &
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
  REAL*8 :: tb_devb
!Rank 3 (Prognostic variables and tracers)
  REAL*8, DIMENSION(np), INTENT(IN) :: ta_dev, wa_dev, oa_dev, fcld_dev
  REAL*8, DIMENSION(np) :: ta_devb, wa_devb, oa_devb, fcld_devb
  REAL*8, DIMENSION(np), INTENT(IN) :: n2o_dev, ch4_dev, cfc11_dev, &
& cfc12_dev, cfc22_dev
  REAL*8, DIMENSION(np + 1), INTENT(IN) :: ple_dev
!Rank 3 (surface types)
  REAL*8, DIMENSION(ns), INTENT(IN) :: fs_dev, tg_dev, tv_dev
  REAL*8, DIMENSION(ns, 10), INTENT(IN) :: eg_dev, ev_dev, rv_dev
!Rank 3 (diagnostic cloud parts)
  REAL*8, DIMENSION(np, 4), INTENT(IN) :: cwc_dev, reff_dev
  REAL*8, DIMENSION(np, 4) :: cwc_devb, reff_devb
!Rank 3 (aerosols)
  REAL*8, DIMENSION(np, nb), INTENT(INOUT) :: taua_dev, ssaa_dev, &
& asya_dev
  REAL*8, DIMENSION(np, nb), INTENT(INOUT) :: taua_devb
!----- OUPUTS -----
  REAL*8, DIMENSION(np + 1) :: flxu_dev
  REAL*8, DIMENSION(np+1) :: flxu_devb
  REAL*8, DIMENSION(np + 1) :: flxd_dev
  REAL*8, DIMENSION(np+1) :: flxd_devb
  REAL*8, DIMENSION(np + 1) :: dfdts_dev
  REAL*8, DIMENSION(np+1) :: dfdts_devb
!----- LOCALS -----
  REAL*8, PARAMETER :: cons_grav=9.80665
  INTEGER, PARAMETER :: nx1=26
  INTEGER, PARAMETER :: no1=21
  INTEGER, PARAMETER :: nc1=30
  INTEGER, PARAMETER :: nh1=31
!Temporary arrays
  REAL*8 :: pa(0:np), dt(0:np)
  REAL*8 :: dtb(0:np)
  REAL*8 :: x1, x2, x3
  REAL*8 :: x1b, x2b, x3b
  REAL*8 :: dh2o(0:np), dcont(0:np), dco2(0:np), do3(0:np)
  REAL*8 :: dh2ob(0:np), dcontb(0:np), dco2b(0:np), do3b(0:np)
  REAL*8 :: dn2o(0:np), dch4(0:np)
  REAL*8 :: df11(0:np), df12(0:np), df22(0:np)
  REAL*8 :: th2o(6), tcon(3), tco2(6)
  REAL*8 :: th2ob(6), tconb(3), tco2b(6)
  REAL*8 :: tn2o(4), tch4(4), tcom(6)
  REAL*8 :: tn2ob(4), tch4b(4), tcomb(6)
  REAL*8 :: tf11, tf12, tf22
  REAL*8 :: tf11b, tf12b, tf22b
  REAL*8 :: blayer(0:np+1), blevel(0:np+1)
  REAL*8 :: blayerb(0:np+1), blevelb(0:np+1)
  REAL*8 :: bd(0:np+1), bu(0:np+1)
  REAL*8 :: bdb(0:np+1), bub(0:np+1)
  REAL*8 :: bs, dbs, rflxs
  REAL*8 :: dp(0:np)
  REAL*8 :: trant, tranal
  REAL*8 :: trantb, tranalb
  REAL*8 :: transfc(0:np+1)
  REAL*8 :: transfcb(0:np+1)
  REAL*8 :: flxu(0:np+1), flxd(0:np+1)
  REAL*8 :: flxub(0:np+1), flxdb(0:np+1)
  REAL*8 :: taerlyr(0:np)
  REAL*8 :: taerlyrb(0:np)
!OVERCAST
  INTEGER :: ncld(3)
  INTEGER :: icx(0:np)
!OVERCAST
  INTEGER :: idx, rc
  INTEGER :: k, l, ip, iw, ibn, ik, iq, isb, k1, k2, ne
  REAL*8 :: enn(0:np)
  REAL*8 :: ennb(0:np)
  REAL*8 :: cldhi, cldmd, cldlw, tcldlyr(0:np), fclr
  REAL*8 :: cldhib, cldmdb, cldlwb, tcldlyrb(0:np), fclrb
  REAL*8 :: x, xx, yy, p1, a1, b1, fk1, a2, b2, fk2
  REAL*8 :: xxb, yyb
  REAL*8 :: w1, ff
  REAL*8 :: ffb
  LOGICAL :: oznbnd, co2bnd, h2otable, conbnd, n2obnd
  LOGICAL :: ch4bnd, combnd, f11bnd, f12bnd, f22bnd, b10bnd
  LOGICAL :: do_aerosol
!Temp arrays and variables for consolidation of tables
  INTEGER, PARAMETER :: max_num_tables=17
  REAL*8 :: exptbl(0:np, max_num_tables)
  REAL*8 :: exptblb(0:np, max_num_tables)
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
  REAL*8 :: fcld_col(np)
  REAL*8 :: fcld_colb(np)
  REAL*8 :: reff_col(np, 4)
  REAL*8 :: reff_colb(np, 4)
  REAL*8 :: cwc_col(np, 4)
  REAL*8 :: cwc_colb(np, 4)
  REAL*8 :: h2oexp_tmp(0:np, 5), conexp_tmp(0:np), co2exp_tmp(0:np, 6), &
& n2oexp_tmp(0:np, 2)
  REAL*8 :: h2oexp_tmpb(0:np, 5), co2exp_tmpb(0:np, 6), n2oexp_tmpb(0:np&
& , 2)
  INTRINSIC MAX
  INTRINSIC EXP
  INTRINSIC MIN
  INTRINSIC LOG
  INTEGER :: arg1
  INTEGER :: branch
  INTEGER :: ad_from
  INTEGER :: ad_count
  INTEGER :: i
  REAL*8, DIMENSION(np, nb), INTENT(INOUT) :: asya_devb
  REAL*8, DIMENSION(np, nb), INTENT(INOUT) :: ssaa_devb
  REAL*8 :: temp2
  REAL*8 :: temp1
  REAL*8 :: temp0
  REAL*8 :: tempb6
  REAL*8 :: tempb5
  REAL*8 :: tempb4
  REAL*8 :: tempb3
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  REAL*8 :: tempb
  REAL*8 :: temp
!BEGIN CALCULATIONS ...
!compute layer pressure (pa) and layer temperature minus 250K (dt) 
  DO k=1,np
    pa(k) = 0.5*(ple_dev(k+1)+ple_dev(k))*0.01
    dp(k) = (ple_dev(k+1)-ple_dev(k))*0.01
! dp in Pascals for getirtau
    dp_pa(k) = ple_dev(k+1) - ple_dev(k)
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
    dh2o(k) = 1.02*wa_dev(k)*dp(k)
    do3(k) = 476.*oa_dev(k)*dp(k)
    dco2(k) = 789.*co2*dp(k)
    dch4(k) = 789.*ch4_dev(k)*dp(k)
    dn2o(k) = 789.*n2o_dev(k)*dp(k)
    df11(k) = 789.*cfc11_dev(k)*dp(k)
    df12(k) = 789.*cfc12_dev(k)*dp(k)
    df22(k) = 789.*cfc22_dev(k)*dp(k)
    IF (dh2o(k) .LT. 1.e-10) THEN
      dh2o(k) = 1.e-10
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
      dh2o(k) = dh2o(k)
    END IF
    IF (do3(k) .LT. 1.e-6) THEN
      do3(k) = 1.e-6
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
      do3(k) = do3(k)
    END IF
    IF (dco2(k) .LT. 1.e-4) THEN
      dco2(k) = 1.e-4
    ELSE
      dco2(k) = dco2(k)
    END IF
!Compute scaled water vapor amount for h2o continuum absorption
!following eq. (4.21).
    xx = pa(k)*0.001618*wa_dev(k)*wa_dev(k)*dp(k)
    dcont(k) = xx*EXP(1800./ta_dev(k)-6.081)
!Fill the reff, cwc, and fcld for the column   
    fcld_col(k) = fcld_dev(k)
    DO l=1,4
      reff_col(k, l) = reff_dev(k, l)
      cwc_col(k, l) = cwc_dev(k, l)
    END DO
  END DO
  IF (ple_dev(1)*0.01 .LT. 0.005) THEN
    CALL PUSHREAL8(dp(0))
    dp(0) = 0.005
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL8(dp(0))
    dp(0) = ple_dev(1)*0.01
    CALL PUSHCONTROL1B(1)
  END IF
  CALL PUSHREAL8(pa(0))
  pa(0) = 0.5*dp(0)
  dt(0) = ta_dev(1) - 250.0
  dh2o(0) = 1.02*wa_dev(1)*dp(0)
  do3(0) = 476.*oa_dev(1)*dp(0)
  dco2(0) = 789.*co2*dp(0)
  dch4(0) = 789.*ch4_dev(1)*dp(0)
  dn2o(0) = 789.*n2o_dev(1)*dp(0)
  df11(0) = 789.*cfc11_dev(1)*dp(0)
  df12(0) = 789.*cfc12_dev(1)*dp(0)
  df22(0) = 789.*cfc22_dev(1)*dp(0)
  IF (dh2o(0) .LT. 1.e-10) THEN
    dh2o(0) = 1.e-10
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    dh2o(0) = dh2o(0)
  END IF
  IF (do3(0) .LT. 1.e-6) THEN
    do3(0) = 1.e-6
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    do3(0) = do3(0)
  END IF
  IF (dco2(0) .LT. 1.e-4) THEN
    dco2(0) = 1.e-4
  ELSE
    dco2(0) = dco2(0)
  END IF
  xx = pa(0)*0.001618*wa_dev(1)*wa_dev(1)*dp(0)
  dcont(0) = xx*EXP(1800./ta_dev(1)-6.081)
!The surface (np+1) is treated as a layer filled with black clouds.
!transfc is the transmittance between the surface and a pressure
!level.
  transfc(np+1) = 1.0
  CALL PUSHINTEGER4(ibn)
  ad_count = 1
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
      CALL PUSHREAL8ARRAY(exptbl, (np+1)*17)
      exptbl = 0.0
!Control packing of the new exponential tables by band
      SELECT CASE  (ibn) 
      CASE (2) 
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 1
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 1
        CALL PUSHCONTROL4B(7)
      CASE (3) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 9
        CALL PUSHCONTROL4B(6)
      CASE (4) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 7
        CALL PUSHINTEGER4(comexp%start)
        comexp%start = 8
        CALL PUSHINTEGER4(comexp%end)
        comexp%end = 13
        CALL PUSHINTEGER4(f11exp%start)
        f11exp%start = 14
        CALL PUSHINTEGER4(f11exp%end)
        f11exp%end = 14
        CALL PUSHINTEGER4(f12exp%start)
        f12exp%start = 15
        CALL PUSHINTEGER4(f12exp%end)
        f12exp%end = 15
        CALL PUSHINTEGER4(f22exp%start)
        f22exp%start = 16
        CALL PUSHINTEGER4(f22exp%end)
        f22exp%end = 16
        CALL PUSHCONTROL4B(5)
      CASE (5) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 7
        CALL PUSHINTEGER4(comexp%start)
        comexp%start = 8
        CALL PUSHINTEGER4(comexp%end)
        comexp%end = 13
        CALL PUSHINTEGER4(f11exp%start)
        f11exp%start = 14
        CALL PUSHINTEGER4(f11exp%end)
        f11exp%end = 14
        CALL PUSHCONTROL4B(4)
      CASE (6) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 7
        CALL PUSHINTEGER4(n2oexp%start)
        n2oexp%start = 8
        CALL PUSHINTEGER4(n2oexp%end)
        n2oexp%end = 11
        CALL PUSHINTEGER4(ch4exp%start)
        ch4exp%start = 12
        CALL PUSHINTEGER4(ch4exp%end)
        ch4exp%end = 15
        CALL PUSHINTEGER4(f12exp%start)
        f12exp%start = 16
        CALL PUSHINTEGER4(f12exp%end)
        f12exp%end = 16
        CALL PUSHINTEGER4(f22exp%start)
        f22exp%start = 17
        CALL PUSHINTEGER4(f22exp%end)
        f22exp%end = 17
        CALL PUSHCONTROL4B(3)
      CASE (7) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 7
        CALL PUSHINTEGER4(n2oexp%start)
        n2oexp%start = 8
        CALL PUSHINTEGER4(n2oexp%end)
        n2oexp%end = 11
        CALL PUSHINTEGER4(ch4exp%start)
        ch4exp%start = 12
        CALL PUSHINTEGER4(ch4exp%end)
        ch4exp%end = 15
        CALL PUSHCONTROL4B(2)
      CASE (9) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHCONTROL4B(1)
      CASE (10) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 5
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 6
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 6
        CALL PUSHINTEGER4(co2exp%start)
        co2exp%start = 7
        CALL PUSHINTEGER4(co2exp%end)
        co2exp%end = 12
        CALL PUSHINTEGER4(n2oexp%start)
        n2oexp%start = 13
        CALL PUSHINTEGER4(n2oexp%end)
        n2oexp%end = 14
        CALL PUSHCONTROL4B(0)
      CASE DEFAULT
        CALL PUSHCONTROL4B(8)
      END SELECT
!blayer is the spectrally integrated planck flux of the mean layer
!temperature derived from eq. (3.11)
!The fitting for the planck flux is valid for the range 160-345 K.
      DO k=1,np
        CALL PLANCK(ibn, cb, ta_dev(k), blayer(k))
      END DO
!Index "0" is the layer above the top of the atmosphere.
      blayer(0) = blayer(1)
      CALL PUSHREAL8(blevel(0))
      blevel(0) = blayer(1)
!Surface emission and reflectivity. See Section 9.
!bs and dbs include the effect of surface emissivity.
      CALL PUSHREAL8(rflxs)
      CALL PUSHREAL8(dbs)
      CALL SFCFLUX(ibn, cb, dcb, ns, fs_dev, tg_dev, eg_dev, tv_dev, &
&            ev_dev, rv_dev, bs, dbs, rflxs)
      blayer(np+1) = bs
!interpolate Planck function at model levels (linear in p)
      DO k=2,np
        CALL PUSHREAL8(blevel(k))
        blevel(k) = (blayer(k-1)*dp(k)+blayer(k)*dp(k-1))/(dp(k-1)+dp(k)&
&         )
      END DO
!Extrapolate blevel(1) from blayer(2) and blayer(1)
      CALL PUSHREAL8(blevel(1))
      blevel(1) = blayer(1) + (blayer(1)-blayer(2))*dp(1)/(dp(1)+dp(2))
      CALL PUSHREAL8(blevel(0))
      blevel(0) = blevel(1)
!If the surface air temperature tb is known, compute blevel(np+1)
      CALL PUSHREAL8(blevel(np+1))
      CALL PLANCK(ibn, cb, tb_dev, blevel(np+1))
!Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!     NOTE: dp_pa is only dims(1:np) as the 0'th level isn't needed in getirtau.
!           Plus, the pressures in getirtau *MUST* be in Pascals.
!     Slots for reff, hydrometeors and tauall are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
      CALL PUSHREAL8ARRAY(enn, np + 1)
      CALL PUSHREAL8ARRAY(tcldlyr, np + 1)
      CALL GETIRTAU1(ibn, np, dp_pa, fcld_col, reff_col, cwc_col, &
&              tcldlyr, enn, aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, &
&              awg_ir, cons_grav)
      DO k=0,np
        CALL PUSHINTEGER4(icx(k))
        icx(k) = k
      END DO
      CALL PUSHINTEGER4ARRAY(ncld, 3)
      CALL PUSHINTEGER4ARRAY(icx, np + 1)
      CALL MKICX(np, ict, icb, enn, icx, ncld)
!Compute optical thickness, single-scattering albedo and asymmetry
!factor for a mixture of "na" aerosol types. Eqs. (7.1)-(7.3)
      IF (do_aerosol) THEN
        CALL PUSHREAL8(taerlyr(0))
        taerlyr(0) = 1.0
        DO k=1,np
!-----taerlyr is the aerosol diffuse transmittance
          CALL PUSHREAL8(taerlyr(k))
          taerlyr(k) = 1.0
          IF (taua_dev(k, ibn) .GT. 0.001) THEN
            IF (ssaa_dev(k, ibn) .GT. 0.001) THEN
              CALL PUSHREAL8(asya_dev(k, ibn))
              asya_dev(k, ibn) = asya_dev(k, ibn)/ssaa_dev(k, ibn)
              CALL PUSHREAL8(ssaa_dev(k, ibn))
              ssaa_dev(k, ibn) = ssaa_dev(k, ibn)/taua_dev(k, ibn)
!Parameterization of aerosol scattering following
              ff = .5 + (.3739+(0.0076+0.1185*asya_dev(k, ibn))*asya_dev&
&               (k, ibn))*asya_dev(k, ibn)
              CALL PUSHREAL8(taua_dev(k, ibn))
              taua_dev(k, ibn) = taua_dev(k, ibn)*(1.-ssaa_dev(k, ibn)*&
&               ff)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            taerlyr(k) = EXP(-(1.66*taua_dev(k, ibn)))
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!Compute the exponential terms (Eq. 8.21) at each layer due to
!water vapor line absorption when k-distribution is used
      IF (.NOT.h2otable .AND. (.NOT.b10bnd)) THEN
        CALL PUSHREAL8ARRAY(exptbl(:, h2oexp%start:h2oexp%end), (np+1)*(&
&                     h2oexp%end-h2oexp%start+1))
        CALL H2OEXPS(ibn, np, dh2o, pa, dt, xkw, aw, bw, pm, mw, exptbl(&
&              :, h2oexp%start:h2oexp%end))
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!compute the exponential terms (Eq. 4.24) at each layer due to
!water vapor continuum absorption.
      CALL PUSHINTEGER4(ne)
      ne = 0
      IF (conbnd) THEN
        ne = 1
        IF (ibn .EQ. 3) ne = 3
        CALL PUSHREAL8ARRAY(exptbl(:, conexp%start:conexp%end), (np+1)*(&
&                     conexp%end-conexp%start+1))
        CALL CONEXPS(ibn, np, dcont, xke, exptbl(:, conexp%start:conexp%&
&              end))
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (trace) THEN
!compute the exponential terms at each layer due to n2o absorption
        IF (n2obnd) THEN
          CALL PUSHREAL8ARRAY(exptbl(:, n2oexp%start:n2oexp%end), (np+1)&
&                       *(n2oexp%end-n2oexp%start+1))
          CALL N2OEXPS(ibn, np, dn2o, pa, dt, exptbl(:, n2oexp%start:&
&                n2oexp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!compute the exponential terms at each layer due to ch4 absorption
        IF (ch4bnd) THEN
          CALL PUSHREAL8ARRAY(exptbl(:, ch4exp%start:ch4exp%end), (np+1)&
&                       *(ch4exp%end-ch4exp%start+1))
          CALL CH4EXPS(ibn, np, dch4, pa, dt, exptbl(:, ch4exp%start:&
&                ch4exp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!Compute the exponential terms due to co2 minor absorption
        IF (combnd) THEN
          CALL PUSHREAL8ARRAY(exptbl(:, comexp%start:comexp%end), (np+1)&
&                       *(comexp%end-comexp%start+1))
          CALL COMEXPS(ibn, np, dco2, dt, exptbl(:, comexp%start:comexp%&
&                end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!Compute the exponential terms due to cfc11 absorption.
!The values of the parameters are given in Table 7.
        IF (f11bnd) THEN
          a1 = 1.26610e-3
          b1 = 3.55940e-6
          fk1 = 1.89736e+1
          a2 = 8.19370e-4
          b2 = 4.67810e-6
          fk2 = 1.01487e+1
          CALL CFCEXPS(ibn, np, a1, b1, fk1, a2, b2, fk2, df11, dt, &
&                exptbl(:, f11exp%start:f11exp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!Compute the exponential terms due to cfc12 absorption.
        IF (f12bnd) THEN
          a1 = 8.77370e-4
          b1 = -5.88440e-6
          fk1 = 1.58104e+1
          a2 = 8.62000e-4
          b2 = -4.22500e-6
          fk2 = 3.70107e+1
          CALL CFCEXPS(ibn, np, a1, b1, fk1, a2, b2, fk2, df12, dt, &
&                exptbl(:, f12exp%start:f12exp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!Compute the exponential terms due to cfc22 absorption.
        IF (f22bnd) THEN
          a1 = 9.65130e-4
          b1 = 1.31280e-5
          fk1 = 6.18536e+0
          a2 = -3.00010e-5
          b2 = 5.25010e-7
          fk2 = 3.27912e+1
          CALL CFCEXPS(ibn, np, a1, b1, fk1, a2, b2, fk2, df22, dt, &
&                exptbl(:, f22exp%start:f22exp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!Compute the exponential terms at each layer in band 10 due to
!h2o line and continuum, co2, and n2o absorption
        IF (b10bnd) THEN
          CALL PUSHREAL8ARRAY(n2oexp_tmp, (np+1)*2)
          CALL PUSHREAL8ARRAY(co2exp_tmp, (np+1)*6)
          CALL PUSHREAL8ARRAY(h2oexp_tmp, (np+1)*5)
          CALL B10EXPS(np, dh2o, dcont, dco2, dn2o, pa, dt, h2oexp_tmp, &
&                exptbl(:, conexp%start:conexp%end), co2exp_tmp, &
&                n2oexp_tmp)
          exptbl(:, h2oexp%start:h2oexp%end) = h2oexp_tmp
          exptbl(:, co2exp%start:co2exp%end) = co2exp_tmp
          exptbl(:, n2oexp%start:n2oexp%end) = n2oexp_tmp
          CALL PUSHCONTROL2B(0)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(2)
      END IF
!blayer(np+1) includes the effect of surface emissivity.
! ALT: this was undefined, check with Max if 0.0 is good value
      CALL PUSHREAL8(bu(0))
      bu(0) = 0.0
      CALL PUSHREAL8(bd(0))
      bd(0) = blayer(1)
      CALL PUSHREAL8(bu(np+1))
      bu(np+1) = blayer(np+1)
!do-loop 1500 is for computing upward (bu) and downward (bd)
!Here, trant is the transmittance of the layer k2-1.
      DO k2=1,np+1
!for h2o line transmission
        IF (.NOT.h2otable) THEN
          th2o = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!for h2o continuum transmission
        CALL PUSHREAL8ARRAY(tcon, 3)
        tcon = 1.0
        x1 = 0.0
        x2 = 0.0
        x3 = 0.0
        CALL PUSHREAL8(trant)
        trant = 1.0
        IF (h2otable) THEN
!Compute water vapor transmittance using table look-up.
          IF (ibn .EQ. 1) THEN
            CALL PUSHREAL8(trant)
            CALL PUSHREAL8(x3)
            CALL PUSHREAL8(x2)
            CALL PUSHREAL8(x1)
            CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w11, p11, dwe, dpe, h11, h12, h13, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ibn .EQ. 2) THEN
            CALL PUSHREAL8(trant)
            CALL PUSHREAL8(x3)
            CALL PUSHREAL8(x2)
            CALL PUSHREAL8(x1)
            CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w11, p11, dwe, dpe, h21, h22, h23, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ibn .EQ. 8) THEN
            CALL PUSHREAL8(trant)
            CALL PUSHREAL8(x3)
            CALL PUSHREAL8(x2)
            CALL PUSHREAL8(x1)
            CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w11, p11, dwe, dpe, h81, h82, h83, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!for water vapor continuum absorption
          IF (conbnd) THEN
! Only the first exp
            CALL PUSHREAL8(tcon(1))
            tcon(1) = tcon(1)*exptbl(k2-1, conexp%start)
            CALL PUSHREAL8(trant)
            trant = trant*tcon(1)
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE IF (.NOT.b10bnd) THEN
!compute water vapor transmittance using k-distribution
          arg1 = k2 - 1
          CALL PUSHREAL8(trant)
          CALL PUSHREAL8ARRAY(tcon, 3)
          CALL PUSHREAL8ARRAY(th2o, 6)
          CALL H2OKDIS(ibn, np, arg1, fkw, gkw, ne, exptbl(:, h2oexp%&
&                start:h2oexp%end), exptbl(:, conexp%start:conexp%end), &
&                th2o, tcon, trant)
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(3)
        END IF
        IF (co2bnd) THEN
!Compute co2 transmittance using table look-up method
          CALL PUSHREAL8(trant)
          CALL PUSHREAL8(x3)
          CALL PUSHREAL8(x2)
          CALL PUSHREAL8(x1)
          CALL TABLUP(nx1, nc1, dco2(k2-1), pa(k2-1), dt(k2-1), x1, x2, &
&               x3, w12, p12, dwe, dpe, c1, c2, c3, trant)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!Always use table look-up to compute o3 transmittance.
        IF (oznbnd) THEN
          CALL PUSHREAL8(trant)
          CALL PUSHREAL8(x3)
          CALL PUSHREAL8(x2)
          CALL PUSHREAL8(x1)
          CALL TABLUP(nx1, no1, do3(k2-1), pa(k2-1), dt(k2-1), x1, x2, &
&               x3, w13, p13, dwe, dpe, oo1, oo2, oo3, trant)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!include aerosol effect
        IF (do_aerosol) THEN
          CALL PUSHREAL8(trant)
          trant = trant*taerlyr(k2-1)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
!Compute upward (bu) and downward (bd) emission of the layer k2-1
        CALL PUSHREAL8(xx)
        xx = (1.-enn(k2-1))*trant
        IF (0.9999 .GT. xx) THEN
          CALL PUSHREAL8(yy)
          yy = xx
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREAL8(yy)
          yy = 0.9999
          CALL PUSHCONTROL1B(1)
        END IF
        IF (0.00001 .LT. yy) THEN
          CALL PUSHCONTROL1B(0)
          yy = yy
        ELSE
          yy = 0.00001
          CALL PUSHCONTROL1B(1)
        END IF
        xx = (blevel(k2-1)-blevel(k2))/LOG(yy)
        CALL PUSHREAL8(bd(k2-1))
        bd(k2-1) = (blevel(k2)-blevel(k2-1)*yy)/(1.0-yy) - xx
        CALL PUSHREAL8(bu(k2-1))
        bu(k2-1) = blevel(k2-1) + blevel(k2) - bd(k2-1)
      END DO
!Initialize fluxes
      CALL PUSHREAL8ARRAY(flxd, np + 2)
      flxd = 0.0
!Compute upward and downward fluxes for each spectral band, ibn.
      DO k1=0,np
!initialization
        CALL PUSHREAL8(cldlw)
        cldlw = 0.0
        CALL PUSHREAL8(cldmd)
        cldmd = 0.0
        CALL PUSHREAL8(cldhi)
        cldhi = 0.0
        CALL PUSHREAL8(tranal)
        tranal = 1.0
!for h2o line transmission
        IF (.NOT.h2otable) THEN
          th2o = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!for h2o continuum transmission
        CALL PUSHREAL8ARRAY(tcon, 3)
        tcon = 1.0
        IF (trace) THEN
!Add trace gases contribution
          IF (n2obnd) THEN
!n2o
            tn2o = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ch4bnd) THEN
!ch4
            tch4 = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (combnd) THEN
!co2-minor
            tcom = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (f11bnd) THEN
!cfc-11
            tf11 = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (f12bnd) THEN
!cfc-12
            tf12 = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (f22bnd) THEN
!cfc-22
            tf22 = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (b10bnd) THEN
!
            th2o = 1.0
            tco2 = 1.0
            tcon(1) = 1.0
            tn2o = 1.0
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
        x1 = 0.0
        x2 = 0.0
        x3 = 0.0
        ad_from = k1 + 1
        DO k2=ad_from,np+1
          CALL PUSHREAL8(trant)
          trant = 1.0
          IF (h2otable) THEN
!Compute water vapor transmittance using table look-up.
            IF (ibn .EQ. 1) THEN
              CALL PUSHREAL8(trant)
              CALL PUSHREAL8(x3)
              CALL PUSHREAL8(x2)
              CALL PUSHREAL8(x1)
              CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, &
&                   x2, x3, w11, p11, dwe, dpe, h11, h12, h13, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (ibn .EQ. 2) THEN
              CALL PUSHREAL8(trant)
              CALL PUSHREAL8(x3)
              CALL PUSHREAL8(x2)
              CALL PUSHREAL8(x1)
              CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, &
&                   x2, x3, w11, p11, dwe, dpe, h21, h22, h23, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (ibn .EQ. 8) THEN
              CALL PUSHREAL8(trant)
              CALL PUSHREAL8(x3)
              CALL PUSHREAL8(x2)
              CALL PUSHREAL8(x1)
              CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, &
&                   x2, x3, w11, p11, dwe, dpe, h81, h82, h83, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (conbnd) THEN
! Only the first exp
              CALL PUSHREAL8(tcon(1))
              tcon(1) = tcon(1)*exptbl(k2-1, conexp%start)
              CALL PUSHREAL8(trant)
              trant = trant*tcon(1)
              CALL PUSHCONTROL2B(0)
            ELSE
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (.NOT.b10bnd) THEN
!Compute water vapor transmittance using k-distribution
            arg1 = k2 - 1
            CALL PUSHREAL8(trant)
            CALL PUSHREAL8ARRAY(tcon, 3)
            CALL PUSHREAL8ARRAY(th2o, 6)
            CALL H2OKDIS(ibn, np, arg1, fkw, gkw, ne, exptbl(:, h2oexp%&
&                  start:h2oexp%end), exptbl(:, conexp%start:conexp%end)&
&                  , th2o, tcon, trant)
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHCONTROL2B(3)
          END IF
          IF (co2bnd) THEN
!Compute co2 transmittance using table look-up method.
            CALL PUSHREAL8(trant)
            CALL PUSHREAL8(x3)
            CALL PUSHREAL8(x2)
            CALL PUSHREAL8(x1)
            CALL TABLUP(nx1, nc1, dco2(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w12, p12, dwe, dpe, c1, c2, c3, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (oznbnd) THEN
!Always use table look-up to compute o3 transmittanc
            CALL PUSHREAL8(trant)
            CALL PUSHREAL8(x3)
            CALL PUSHREAL8(x2)
            CALL PUSHREAL8(x1)
            CALL TABLUP(nx1, no1, do3(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w13, p13, dwe, dpe, oo1, oo2, oo3, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (trace) THEN
!Add trace gas effects
            IF (n2obnd) THEN
!n2o
              arg1 = k2 - 1
              CALL PUSHREAL8(trant)
              CALL PUSHREAL8ARRAY(tn2o, 4)
              CALL N2OKDIS(ibn, np, arg1, exptbl(:, n2oexp%start:n2oexp%&
&                    end), tn2o, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (ch4bnd) THEN
!ch4
              arg1 = k2 - 1
              CALL PUSHREAL8(trant)
              CALL PUSHREAL8ARRAY(tch4, 4)
              CALL CH4KDIS(ibn, np, arg1, exptbl(:, ch4exp%start:ch4exp%&
&                    end), tch4, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (combnd) THEN
!co2-minor
              arg1 = k2 - 1
              CALL PUSHREAL8(trant)
              CALL PUSHREAL8ARRAY(tcom, 6)
              CALL COMKDIS(ibn, np, arg1, exptbl(:, comexp%start:comexp%&
&                    end), tcom, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (f11bnd) THEN
!cfc11
              arg1 = k2 - 1
              CALL PUSHREAL8(trant)
              CALL PUSHREAL8(tf11)
              CALL CFCKDIS(np, arg1, exptbl(:, f11exp%start:f11exp%end)&
&                    , tf11, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (f12bnd) THEN
!cfc12
              arg1 = k2 - 1
              CALL PUSHREAL8(trant)
              CALL PUSHREAL8(tf12)
              CALL CFCKDIS(np, arg1, exptbl(:, f12exp%start:f12exp%end)&
&                    , tf12, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (f22bnd) THEN
!CFC22
              arg1 = k2 - 1
              CALL PUSHREAL8(trant)
              CALL PUSHREAL8(tf22)
              CALL CFCKDIS(np, arg1, exptbl(:, f22exp%start:f22exp%end)&
&                    , tf22, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (b10bnd) THEN
              arg1 = k2 - 1
              CALL PUSHREAL8ARRAY(tn2o, 4)
              CALL PUSHREAL8ARRAY(tco2, 6)
              CALL PUSHREAL8ARRAY(tcon, 3)
              CALL PUSHREAL8ARRAY(th2o, 6)
              CALL B10KDIS(np, arg1, exptbl(:, h2oexp%start:h2oexp%end)&
&                    , exptbl(:, conexp%start:conexp%end), exptbl(:, &
&                    co2exp%start:co2exp%end), exptbl(:, n2oexp%start:&
&                    n2oexp%end), th2o, tcon, tco2, tn2o, trant)
              CALL PUSHCONTROL2B(0)
            ELSE
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE
            CALL PUSHCONTROL2B(2)
          END IF
          IF (do_aerosol) THEN
            CALL PUSHREAL8(tranal)
            tranal = tranal*taerlyr(k2-1)
            CALL PUSHREAL8(trant)
            trant = trant*tranal
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (enn(k2-1) .GE. 0.001) THEN
            CALL PUSHREAL8(cldlw)
            CALL PUSHREAL8(cldmd)
            CALL PUSHREAL8(cldhi)
            CALL CLDOVLP(np, k1, k2, ict, icb, icx, ncld, enn, tcldlyr, &
&                  cldhi, cldmd, cldlw)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(fclr)
          fclr = (1.0-cldhi)*(1.0-cldmd)*(1.0-cldlw)
          IF (k2 .EQ. k1 + 1 .AND. ibn .NE. 10) THEN
            flxd(k2) = flxd(k2) + bd(k1)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(xx)
          IF (k1 .EQ. 0) THEN
!mjs  bd(-1) is not defined
            xx = -(trant*bd(k1))
            CALL PUSHCONTROL1B(0)
          ELSE
            xx = trant*(bd(k1-1)-bd(k1))
            CALL PUSHCONTROL1B(1)
          END IF
          flxd(k2) = flxd(k2) + xx*fclr
        END DO
        CALL PUSHINTEGER4(ad_from)
        CALL PUSHREAL8(transfc(k1))
        transfc(k1) = trant*fclr
        IF (k1 .GT. 0) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
      IF (.NOT.b10bnd) THEN
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      CALL PUSHINTEGER4(ibn)
      ad_count = ad_count + 1
    END IF
  END DO
  CALL PUSHCONTROL1B(0)
  CALL PUSHINTEGER4(ad_count)
  GOTO 110
 100 CALL PUSHCONTROL1B(1)
  CALL PUSHINTEGER4(ad_count)
 110 CALL POPINTEGER4(ad_count)
  DO i=1,ad_count
    IF (i .EQ. 1) THEN
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        dh2ob = 0.0_8
        dcontb = 0.0_8
        n2oexp_tmpb = 0.0_8
        dtb = 0.0_8
        reff_colb = 0.0_8
        blayerb = 0.0_8
        tn2ob = 0.0_8
        transfcb = 0.0_8
        bdb = 0.0_8
        h2oexp_tmpb = 0.0_8
        tcomb = 0.0_8
        tf11b = 0.0_8
        tf12b = 0.0_8
        trantb = 0.0_8
        bub = 0.0_8
        th2ob = 0.0_8
        tcldlyrb = 0.0_8
        tch4b = 0.0_8
        tf22b = 0.0_8
        ennb = 0.0_8
        fcld_colb = 0.0_8
        fclrb = 0.0_8
        cwc_colb = 0.0_8
        tco2b = 0.0_8
        co2exp_tmpb = 0.0_8
        taerlyrb = 0.0_8
        do3b = 0.0_8
        blevelb = 0.0_8
      ELSE
        dh2ob = 0.0_8
        dcontb = 0.0_8
        n2oexp_tmpb = 0.0_8
        dtb = 0.0_8
        reff_colb = 0.0_8
        blayerb = 0.0_8
        tn2ob = 0.0_8
        transfcb = 0.0_8
        bdb = 0.0_8
        h2oexp_tmpb = 0.0_8
        tcomb = 0.0_8
        tf11b = 0.0_8
        tf12b = 0.0_8
        trantb = 0.0_8
        bub = 0.0_8
        th2ob = 0.0_8
        tcldlyrb = 0.0_8
        tch4b = 0.0_8
        tf22b = 0.0_8
        ennb = 0.0_8
        fcld_colb = 0.0_8
        fclrb = 0.0_8
        cwc_colb = 0.0_8
        tco2b = 0.0_8
        co2exp_tmpb = 0.0_8
        taerlyrb = 0.0_8
        do3b = 0.0_8
        blevelb = 0.0_8
      END IF
    ELSE
      flxub = 0.0_8
      flxdb = 0.0_8
      DO k=np+1,1,-1
        flxdb(k) = flxdb(k) + flxd_devb(k)
        flxub(k) = flxub(k) + flxu_devb(k)
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO k=np+1,1,-1
          flxdb(np+1) = flxdb(np+1) - rflxs*transfc(k)*flxub(k)
          transfcb(k) = transfcb(k) - rflxs*flxd(np+1)*flxub(k)
        END DO
        blayerb(np+1) = blayerb(np+1) - flxub(np+1)
        flxub(np+1) = 0.0_8
      END IF
      exptblb = 0.0_8
      DO k1=np,0,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) transfcb(k1) = transfcb(k1) - dbs*dfdts_devb(&
&           k1)
        CALL POPREAL8(transfc(k1))
        trantb = trantb + fclr*transfcb(k1)
        fclrb = fclrb + trant*transfcb(k1)
        transfcb(k1) = 0.0_8
        cldhib = 0.0_8
        tconb = 0.0_8
        tranalb = 0.0_8
        x1b = 0.0_8
        x2b = 0.0_8
        x3b = 0.0_8
        cldlwb = 0.0_8
        cldmdb = 0.0_8
        CALL POPINTEGER4(ad_from)
        DO k2=np+1,ad_from,-1
          xxb = fclr*flxdb(k2)
          fclrb = fclrb + xx*flxdb(k2)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            trantb = trantb - bd(k1)*xxb
            bdb(k1) = bdb(k1) - trant*xxb
          ELSE
            trantb = trantb + (bd(k1-1)-bd(k1))*xxb
            bdb(k1-1) = bdb(k1-1) + trant*xxb
            bdb(k1) = bdb(k1) - trant*xxb
          END IF
          xx = trant*(bu(k2-1)-bu(k2))
          xxb = fclr*flxub(k1)
          fclrb = fclrb + xx*flxub(k1)
          CALL POPREAL8(xx)
          trantb = trantb + (bu(k2-1)-bu(k2))*xxb
          bub(k2-1) = bub(k2-1) + trant*xxb
          bub(k2) = bub(k2) - trant*xxb
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            bdb(k1) = bdb(k1) + flxdb(k2)
            bub(k1) = bub(k1) - flxub(k1)
          END IF
          CALL POPREAL8(fclr)
          cldhib = cldhib - (1.0-cldmd)*(1.0-cldlw)*fclrb
          cldmdb = cldmdb - (1.0-cldhi)*(1.0-cldlw)*fclrb
          cldlwb = cldlwb - (1.0-cldhi)*(1.0-cldmd)*fclrb
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(cldhi)
            CALL POPREAL8(cldmd)
            CALL POPREAL8(cldlw)
            CALL CLDOVLP_B(np, k1, k2, ict, icb, icx, ncld, enn, ennb, &
&                    tcldlyr, tcldlyrb, cldhi, cldhib, cldmd, cldmdb, &
&                    cldlw, cldlwb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(trant)
            tranalb = tranalb + trant*trantb
            trantb = tranal*trantb
            CALL POPREAL8(tranal)
            taerlyrb(k2-1) = taerlyrb(k2-1) + tranal*tranalb
            tranalb = taerlyr(k2-1)*tranalb
          END IF
          CALL POPCONTROL2B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL8ARRAY(th2o, 6)
            CALL POPREAL8ARRAY(tcon, 3)
            CALL POPREAL8ARRAY(tco2, 6)
            CALL POPREAL8ARRAY(tn2o, 4)
            CALL B10KDIS_B(np, arg1, exptbl(:, h2oexp%start:h2oexp%end)&
&                    , exptblb(:, h2oexp%start:h2oexp%end), exptbl(:, &
&                    conexp%start:conexp%end), exptblb(:, conexp%start:&
&                    conexp%end), exptbl(:, co2exp%start:co2exp%end), &
&                    exptblb(:, co2exp%start:co2exp%end), exptbl(:, &
&                    n2oexp%start:n2oexp%end), exptblb(:, n2oexp%start:&
&                    n2oexp%end), th2o, th2ob, tcon, tconb, tco2, tco2b&
&                    , tn2o, tn2ob, trant, trantb)
            trantb = 0.0_8
          ELSE IF (branch .NE. 1) THEN
            GOTO 120
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL8(tf22)
            CALL POPREAL8(trant)
            CALL CFCKDIS_B(np, arg1, exptbl(:, f22exp%start:f22exp%end)&
&                    , exptblb(:, f22exp%start:f22exp%end), tf22, tf22b&
&                    , trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL8(tf12)
            CALL POPREAL8(trant)
            CALL CFCKDIS_B(np, arg1, exptbl(:, f12exp%start:f12exp%end)&
&                    , exptblb(:, f12exp%start:f12exp%end), tf12, tf12b&
&                    , trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL8(tf11)
            CALL POPREAL8(trant)
            CALL CFCKDIS_B(np, arg1, exptbl(:, f11exp%start:f11exp%end)&
&                    , exptblb(:, f11exp%start:f11exp%end), tf11, tf11b&
&                    , trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL8ARRAY(tcom, 6)
            CALL POPREAL8(trant)
            CALL COMKDIS_B(ibn, np, arg1, exptbl(:, comexp%start:comexp%&
&                    end), exptblb(:, comexp%start:comexp%end), tcom, &
&                    tcomb, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL8ARRAY(tch4, 4)
            CALL POPREAL8(trant)
            CALL CH4KDIS_B(ibn, np, arg1, exptbl(:, ch4exp%start:ch4exp%&
&                    end), exptblb(:, ch4exp%start:ch4exp%end), tch4, &
&                    tch4b, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL8ARRAY(tn2o, 4)
            CALL POPREAL8(trant)
            CALL N2OKDIS_B(ibn, np, arg1, exptbl(:, n2oexp%start:n2oexp%&
&                    end), exptblb(:, n2oexp%start:n2oexp%end), tn2o, &
&                    tn2ob, trant, trantb)
          END IF
 120      CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(x1)
            CALL POPREAL8(x2)
            CALL POPREAL8(x3)
            CALL POPREAL8(trant)
            CALL TABLUP_B(nx1, no1, do3(k2-1), do3b(k2-1), pa(k2-1), dt(&
&                   k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w13, &
&                   p13, dwe, dpe, oo1, oo2, oo3, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(x1)
            CALL POPREAL8(x2)
            CALL POPREAL8(x3)
            CALL POPREAL8(trant)
            dco2b = 0.0_8
            CALL TABLUP_B(nx1, nc1, dco2(k2-1), dco2b(k2-1), pa(k2-1), &
&                   dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w12&
&                   , p12, dwe, dpe, c1, c2, c3, trant, trantb)
          END IF
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              CALL POPREAL8(trant)
              tconb(1) = tconb(1) + trant*trantb
              trantb = tcon(1)*trantb
              CALL POPREAL8(tcon(1))
              exptblb(k2-1, conexp%start) = exptblb(k2-1, conexp%start) &
&               + tcon(1)*tconb(1)
              tconb(1) = exptbl(k2-1, conexp%start)*tconb(1)
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL8(x1)
              CALL POPREAL8(x2)
              CALL POPREAL8(x3)
              CALL POPREAL8(trant)
              CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1)&
&                     , dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, &
&                     w11, p11, dwe, dpe, h81, h82, h83, trant, trantb)
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL8(x1)
              CALL POPREAL8(x2)
              CALL POPREAL8(x3)
              CALL POPREAL8(trant)
              CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1)&
&                     , dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, &
&                     w11, p11, dwe, dpe, h21, h22, h23, trant, trantb)
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL8(x1)
              CALL POPREAL8(x2)
              CALL POPREAL8(x3)
              CALL POPREAL8(trant)
              CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1)&
&                     , dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, &
&                     w11, p11, dwe, dpe, h11, h12, h13, trant, trantb)
            END IF
          ELSE IF (branch .EQ. 2) THEN
            arg1 = k2 - 1
            CALL POPREAL8ARRAY(th2o, 6)
            CALL POPREAL8ARRAY(tcon, 3)
            CALL POPREAL8(trant)
            CALL H2OKDIS_B(ibn, np, arg1, fkw, gkw, ne, exptbl(:, h2oexp&
&                    %start:h2oexp%end), exptblb(:, h2oexp%start:h2oexp%&
&                    end), exptbl(:, conexp%start:conexp%end), exptblb(:&
&                    , conexp%start:conexp%end), th2o, th2ob, tcon, &
&                    tconb, trant, trantb)
          END IF
          CALL POPREAL8(trant)
          trantb = 0.0_8
          fclrb = 0.0_8
        END DO
        CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          tn2ob = 0.0_8
          th2ob = 0.0_8
          tco2b = 0.0_8
        ELSE IF (branch .NE. 1) THEN
          GOTO 130
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tf22b = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tf12b = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tf11b = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tcomb = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tch4b = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tn2ob = 0.0_8
 130    CALL POPREAL8ARRAY(tcon, 3)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) th2ob = 0.0_8
        CALL POPREAL8(tranal)
        CALL POPREAL8(cldhi)
        CALL POPREAL8(cldmd)
        CALL POPREAL8(cldlw)
      END DO
      CALL POPREAL8ARRAY(flxd, np + 2)
      DO k2=np+1,1,-1
        bdb(k2-1) = bdb(k2-1) - bub(k2-1)
        tempb5 = bdb(k2-1)/(1.0-yy)
        xxb = -bdb(k2-1)
        temp2 = LOG(yy)
        tempb6 = xxb/temp2
        CALL POPREAL8(bu(k2-1))
        blevelb(k2-1) = blevelb(k2-1) + bub(k2-1)
        blevelb(k2) = blevelb(k2) + tempb5 + bub(k2-1)
        bub(k2-1) = 0.0_8
        CALL POPREAL8(bd(k2-1))
        blevelb(k2-1) = blevelb(k2-1) + tempb6 - yy*tempb5
        yyb = ((blevel(k2)-blevel(k2-1)*yy)/(1.0-yy)-blevel(k2-1))*&
&         tempb5 - (blevel(k2-1)-blevel(k2))*tempb6/(temp2*yy)
        bdb(k2-1) = 0.0_8
        blevelb(k2) = blevelb(k2) - tempb6
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) yyb = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(yy)
          xxb = yyb
        ELSE
          CALL POPREAL8(yy)
          xxb = 0.0_8
        END IF
        CALL POPREAL8(xx)
        ennb(k2-1) = ennb(k2-1) - trant*xxb
        trantb = trantb + (1.-enn(k2-1))*xxb
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          CALL POPREAL8(trant)
          taerlyrb(k2-1) = taerlyrb(k2-1) + trant*trantb
          trantb = taerlyr(k2-1)*trantb
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(x1)
          CALL POPREAL8(x2)
          CALL POPREAL8(x3)
          CALL POPREAL8(trant)
          x1b = 0.0_8
          x2b = 0.0_8
          x3b = 0.0_8
          CALL TABLUP_B(nx1, no1, do3(k2-1), do3b(k2-1), pa(k2-1), dt(k2&
&                 -1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w13, p13, &
&                 dwe, dpe, oo1, oo2, oo3, trant, trantb)
        ELSE
          x1b = 0.0_8
          x2b = 0.0_8
          x3b = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(x1)
          CALL POPREAL8(x2)
          CALL POPREAL8(x3)
          CALL POPREAL8(trant)
          dco2b = 0.0_8
          CALL TABLUP_B(nx1, nc1, dco2(k2-1), dco2b(k2-1), pa(k2-1), dt(&
&                 k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w12, p12&
&                 , dwe, dpe, c1, c2, c3, trant, trantb)
        END IF
        CALL POPCONTROL2B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            tconb = 0.0_8
            CALL POPREAL8(trant)
            tconb(1) = tconb(1) + trant*trantb
            trantb = tcon(1)*trantb
            CALL POPREAL8(tcon(1))
            exptblb(k2-1, conexp%start) = exptblb(k2-1, conexp%start) + &
&             tcon(1)*tconb(1)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(x1)
            CALL POPREAL8(x2)
            CALL POPREAL8(x3)
            CALL POPREAL8(trant)
            CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1), &
&                   dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w11&
&                   , p11, dwe, dpe, h81, h82, h83, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(x1)
            CALL POPREAL8(x2)
            CALL POPREAL8(x3)
            CALL POPREAL8(trant)
            CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1), &
&                   dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w11&
&                   , p11, dwe, dpe, h21, h22, h23, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(x1)
            CALL POPREAL8(x2)
            CALL POPREAL8(x3)
            CALL POPREAL8(trant)
            CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1), &
&                   dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w11&
&                   , p11, dwe, dpe, h11, h12, h13, trant, trantb)
          END IF
        ELSE IF (branch .EQ. 2) THEN
          arg1 = k2 - 1
          CALL POPREAL8ARRAY(th2o, 6)
          CALL POPREAL8ARRAY(tcon, 3)
          CALL POPREAL8(trant)
          tconb = 0.0_8
          CALL H2OKDIS_B(ibn, np, arg1, fkw, gkw, ne, exptbl(:, h2oexp%&
&                  start:h2oexp%end), exptblb(:, h2oexp%start:h2oexp%end&
&                  ), exptbl(:, conexp%start:conexp%end), exptblb(:, &
&                  conexp%start:conexp%end), th2o, th2ob, tcon, tconb, &
&                  trant, trantb)
        END IF
        CALL POPREAL8(trant)
        CALL POPREAL8ARRAY(tcon, 3)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) th2ob = 0.0_8
        trantb = 0.0_8
      END DO
      CALL POPREAL8(bu(np+1))
      blayerb(np+1) = blayerb(np+1) + bub(np+1)
      bub(np+1) = 0.0_8
      CALL POPREAL8(bd(0))
      blayerb(1) = blayerb(1) + bdb(0)
      bdb(0) = 0.0_8
      CALL POPREAL8(bu(0))
      bub(0) = 0.0_8
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        n2oexp_tmpb = n2oexp_tmpb + exptblb(:, n2oexp%start:n2oexp%end)
        exptblb(:, n2oexp%start:n2oexp%end) = 0.0_8
        co2exp_tmpb = co2exp_tmpb + exptblb(:, co2exp%start:co2exp%end)
        exptblb(:, co2exp%start:co2exp%end) = 0.0_8
        h2oexp_tmpb = h2oexp_tmpb + exptblb(:, h2oexp%start:h2oexp%end)
        exptblb(:, h2oexp%start:h2oexp%end) = 0.0_8
        CALL POPREAL8ARRAY(h2oexp_tmp, (np+1)*5)
        CALL POPREAL8ARRAY(co2exp_tmp, (np+1)*6)
        CALL POPREAL8ARRAY(n2oexp_tmp, (np+1)*2)
        CALL B10EXPS_B(np, dh2o, dh2ob, dcont, dcontb, dco2, dn2o, pa, &
&                dt, dtb, h2oexp_tmp, h2oexp_tmpb, exptbl(:, conexp%&
&                start:conexp%end), exptblb(:, conexp%start:conexp%end)&
&                , co2exp_tmp, co2exp_tmpb, n2oexp_tmp, n2oexp_tmpb)
      ELSE IF (branch .NE. 1) THEN
        GOTO 140
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        a1 = 9.65130e-4
        a2 = -3.00010e-5
        b1 = 1.31280e-5
        b2 = 5.25010e-7
        fk1 = 6.18536e+0
        fk2 = 3.27912e+1
        CALL CFCEXPS_B(ibn, np, a1, b1, fk1, a2, b2, fk2, df22, dt, dtb&
&                , exptbl(:, f22exp%start:f22exp%end), exptblb(:, f22exp&
&                %start:f22exp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        a1 = 8.77370e-4
        a2 = 8.62000e-4
        b1 = -5.88440e-6
        b2 = -4.22500e-6
        fk1 = 1.58104e+1
        fk2 = 3.70107e+1
        CALL CFCEXPS_B(ibn, np, a1, b1, fk1, a2, b2, fk2, df12, dt, dtb&
&                , exptbl(:, f12exp%start:f12exp%end), exptblb(:, f12exp&
&                %start:f12exp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        a1 = 1.26610e-3
        a2 = 8.19370e-4
        b1 = 3.55940e-6
        b2 = 4.67810e-6
        fk1 = 1.89736e+1
        fk2 = 1.01487e+1
        CALL CFCEXPS_B(ibn, np, a1, b1, fk1, a2, b2, fk2, df11, dt, dtb&
&                , exptbl(:, f11exp%start:f11exp%end), exptblb(:, f11exp&
&                %start:f11exp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(exptbl(:, comexp%start:comexp%end), (np+1)*(&
&                    comexp%end-comexp%start+1))
        CALL COMEXPS_B(ibn, np, dco2, dt, dtb, exptbl(:, comexp%start:&
&                comexp%end), exptblb(:, comexp%start:comexp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(exptbl(:, ch4exp%start:ch4exp%end), (np+1)*(&
&                    ch4exp%end-ch4exp%start+1))
        CALL CH4EXPS_B(ibn, np, dch4, pa, dt, dtb, exptbl(:, ch4exp%&
&                start:ch4exp%end), exptblb(:, ch4exp%start:ch4exp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(exptbl(:, n2oexp%start:n2oexp%end), (np+1)*(&
&                    n2oexp%end-n2oexp%start+1))
        CALL N2OEXPS_B(ibn, np, dn2o, pa, dt, dtb, exptbl(:, n2oexp%&
&                start:n2oexp%end), exptblb(:, n2oexp%start:n2oexp%end))
      END IF
 140  CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(exptbl(:, conexp%start:conexp%end), (np+1)*(&
&                    conexp%end-conexp%start+1))
        CALL CONEXPS_B(ibn, np, dcont, dcontb, xke, exptbl(:, conexp%&
&                start:conexp%end), exptblb(:, conexp%start:conexp%end))
      END IF
      CALL POPINTEGER4(ne)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(exptbl(:, h2oexp%start:h2oexp%end), (np+1)*(&
&                    h2oexp%end-h2oexp%start+1))
        CALL H2OEXPS_B(ibn, np, dh2o, dh2ob, pa, dt, dtb, xkw, aw, bw, &
&                pm, mw, exptbl(:, h2oexp%start:h2oexp%end), exptblb(:, &
&                h2oexp%start:h2oexp%end))
        exptblb(:, h2oexp%start:h2oexp%end) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO k=np,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            taua_devb(k, ibn) = taua_devb(k, ibn) - EXP(-(1.66*taua_dev(&
&             k, ibn)))*1.66*taerlyrb(k)
            taerlyrb(k) = 0.0_8
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              ff = .5 + (.3739+(0.0076+0.1185*asya_dev(k, ibn))*asya_dev&
&               (k, ibn))*asya_dev(k, ibn)
              CALL POPREAL8(taua_dev(k, ibn))
              tempb1 = taua_dev(k, ibn)*taua_devb(k, ibn)
              ssaa_devb(k, ibn) = ssaa_devb(k, ibn) - ff*tempb1
              ffb = -(ssaa_dev(k, ibn)*tempb1)
              taua_devb(k, ibn) = (1.-ssaa_dev(k, ibn)*ff)*taua_devb(k, &
&               ibn)
              tempb2 = asya_dev(k, ibn)*ffb
              temp1 = 0.1185*asya_dev(k, ibn) + 0.0076
              asya_devb(k, ibn) = asya_devb(k, ibn) + (temp1*asya_dev(k&
&               , ibn)+.3739)*ffb + (temp1+asya_dev(k, ibn)*0.1185)*&
&               tempb2
              CALL POPREAL8(ssaa_dev(k, ibn))
              tempb3 = ssaa_devb(k, ibn)/taua_dev(k, ibn)
              taua_devb(k, ibn) = taua_devb(k, ibn) - ssaa_dev(k, ibn)*&
&               tempb3/taua_dev(k, ibn)
              CALL POPREAL8(asya_dev(k, ibn))
              tempb4 = asya_devb(k, ibn)/ssaa_dev(k, ibn)
              ssaa_devb(k, ibn) = tempb3 - asya_dev(k, ibn)*tempb4/&
&               ssaa_dev(k, ibn)
              asya_devb(k, ibn) = tempb4
            END IF
          END IF
          CALL POPREAL8(taerlyr(k))
          taerlyrb(k) = 0.0_8
        END DO
        CALL POPREAL8(taerlyr(0))
        taerlyrb(0) = 0.0_8
      END IF
      CALL POPINTEGER4ARRAY(icx, np + 1)
      CALL POPINTEGER4ARRAY(ncld, 3)
      DO k=np,0,-1
        CALL POPINTEGER4(icx(k))
      END DO
      CALL POPREAL8ARRAY(tcldlyr, np + 1)
      CALL POPREAL8ARRAY(enn, np + 1)
      CALL GETIRTAU1_B(ibn, np, dp_pa, fcld_col, fcld_colb, reff_col, &
&                reff_colb, cwc_col, cwc_colb, tcldlyr, tcldlyrb, enn, &
&                ennb, aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, awg_ir, &
&                cons_grav)
      CALL POPREAL8(blevel(np+1))
      CALL PLANCK_B(ibn, cb, tb_dev, tb_devb, blevel(np+1), blevelb(np+1&
&             ))
      blevelb(np+1) = 0.0_8
      CALL POPREAL8(blevel(0))
      blevelb(1) = blevelb(1) + blevelb(0)
      blevelb(0) = 0.0_8
      CALL POPREAL8(blevel(1))
      tempb0 = dp(1)*blevelb(1)/(dp(1)+dp(2))
      blayerb(1) = blayerb(1) + tempb0 + blevelb(1)
      blayerb(2) = blayerb(2) - tempb0
      blevelb(1) = 0.0_8
      DO k=np,2,-1
        CALL POPREAL8(blevel(k))
        tempb = blevelb(k)/(dp(k-1)+dp(k))
        blayerb(k-1) = blayerb(k-1) + dp(k)*tempb
        blayerb(k) = blayerb(k) + dp(k-1)*tempb
        blevelb(k) = 0.0_8
      END DO
      blayerb(np+1) = 0.0_8
      CALL POPREAL8(dbs)
      CALL POPREAL8(rflxs)
      CALL POPREAL8(blevel(0))
      blayerb(1) = blayerb(1) + blayerb(0) + blevelb(0)
      blevelb(0) = 0.0_8
      blayerb(0) = 0.0_8
      DO k=np,1,-1
        CALL PLANCK_B(ibn, cb, ta_dev(k), ta_devb(k), blayer(k), blayerb&
&               (k))
        blayerb(k) = 0.0_8
      END DO
      CALL POPCONTROL4B(branch)
      IF (branch .LT. 4) THEN
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            CALL POPINTEGER4(n2oexp%end)
            CALL POPINTEGER4(n2oexp%start)
            CALL POPINTEGER4(co2exp%end)
            CALL POPINTEGER4(co2exp%start)
            CALL POPINTEGER4(conexp%end)
            CALL POPINTEGER4(conexp%start)
            CALL POPINTEGER4(h2oexp%end)
            CALL POPINTEGER4(h2oexp%start)
          ELSE
            CALL POPINTEGER4(h2oexp%end)
            CALL POPINTEGER4(h2oexp%start)
          END IF
        ELSE IF (branch .EQ. 2) THEN
          CALL POPINTEGER4(ch4exp%end)
          CALL POPINTEGER4(ch4exp%start)
          CALL POPINTEGER4(n2oexp%end)
          CALL POPINTEGER4(n2oexp%start)
          CALL POPINTEGER4(conexp%end)
          CALL POPINTEGER4(conexp%start)
          CALL POPINTEGER4(h2oexp%end)
          CALL POPINTEGER4(h2oexp%start)
        ELSE
          CALL POPINTEGER4(f22exp%end)
          CALL POPINTEGER4(f22exp%start)
          CALL POPINTEGER4(f12exp%end)
          CALL POPINTEGER4(f12exp%start)
          CALL POPINTEGER4(ch4exp%end)
          CALL POPINTEGER4(ch4exp%start)
          CALL POPINTEGER4(n2oexp%end)
          CALL POPINTEGER4(n2oexp%start)
          CALL POPINTEGER4(conexp%end)
          CALL POPINTEGER4(conexp%start)
          CALL POPINTEGER4(h2oexp%end)
          CALL POPINTEGER4(h2oexp%start)
        END IF
      ELSE IF (branch .LT. 6) THEN
        IF (branch .EQ. 4) THEN
          CALL POPINTEGER4(f11exp%end)
          CALL POPINTEGER4(f11exp%start)
          CALL POPINTEGER4(comexp%end)
          CALL POPINTEGER4(comexp%start)
          CALL POPINTEGER4(conexp%end)
          CALL POPINTEGER4(conexp%start)
          CALL POPINTEGER4(h2oexp%end)
          CALL POPINTEGER4(h2oexp%start)
        ELSE
          CALL POPINTEGER4(f22exp%end)
          CALL POPINTEGER4(f22exp%start)
          CALL POPINTEGER4(f12exp%end)
          CALL POPINTEGER4(f12exp%start)
          CALL POPINTEGER4(f11exp%end)
          CALL POPINTEGER4(f11exp%start)
          CALL POPINTEGER4(comexp%end)
          CALL POPINTEGER4(comexp%start)
          CALL POPINTEGER4(conexp%end)
          CALL POPINTEGER4(conexp%start)
          CALL POPINTEGER4(h2oexp%end)
          CALL POPINTEGER4(h2oexp%start)
        END IF
      ELSE IF (branch .EQ. 6) THEN
        CALL POPINTEGER4(conexp%end)
        CALL POPINTEGER4(conexp%start)
        CALL POPINTEGER4(h2oexp%end)
        CALL POPINTEGER4(h2oexp%start)
      ELSE IF (branch .EQ. 7) THEN
        CALL POPINTEGER4(conexp%end)
        CALL POPINTEGER4(conexp%start)
      END IF
      CALL POPREAL8ARRAY(exptbl, (np+1)*17)
    END IF
    CALL POPINTEGER4(ibn)
  END DO
  DO k=np+1,1,-1
    dfdts_devb(k) = 0.0_8
    flxd_devb(k) = 0.0_8
    flxu_devb(k) = 0.0_8
  END DO
  xx = pa(0)*0.001618*wa_dev(1)*wa_dev(1)*dp(0)
  temp0 = 1800./ta_dev(1)
  xxb = EXP(temp0-6.081)*dcontb(0)
  ta_devb(1) = ta_devb(1) - EXP(temp0-6.081)*xx*temp0*dcontb(0)/ta_dev(1&
&   )
  dcontb(0) = 0.0_8
  wa_devb(1) = wa_devb(1) + pa(0)*0.001618*dp(0)*2*wa_dev(1)*xxb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) do3b(0) = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) dh2ob(0) = 0.0_8
  oa_devb(1) = oa_devb(1) + dp(0)*476.*do3b(0)
  do3b(0) = 0.0_8
  wa_devb(1) = wa_devb(1) + dp(0)*1.02*dh2ob(0)
  dh2ob(0) = 0.0_8
  ta_devb(1) = ta_devb(1) + dtb(0)
  dtb(0) = 0.0_8
  CALL POPREAL8(pa(0))
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    CALL POPREAL8(dp(0))
  ELSE
    CALL POPREAL8(dp(0))
  END IF
  DO k=np,1,-1
    DO l=4,1,-1
      cwc_devb(k, l) = cwc_devb(k, l) + cwc_colb(k, l)
      cwc_colb(k, l) = 0.0_8
      reff_devb(k, l) = reff_devb(k, l) + reff_colb(k, l)
      reff_colb(k, l) = 0.0_8
    END DO
    fcld_devb(k) = fcld_devb(k) + fcld_colb(k)
    fcld_colb(k) = 0.0_8
    xx = pa(k)*0.001618*wa_dev(k)*wa_dev(k)*dp(k)
    temp = 1800./ta_dev(k)
    xxb = EXP(temp-6.081)*dcontb(k)
    ta_devb(k) = ta_devb(k) - EXP(temp-6.081)*xx*temp*dcontb(k)/ta_dev(k&
&     )
    dcontb(k) = 0.0_8
    wa_devb(k) = wa_devb(k) + pa(k)*0.001618*dp(k)*2*wa_dev(k)*xxb
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) do3b(k) = 0.0_8
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) dh2ob(k) = 0.0_8
    oa_devb(k) = oa_devb(k) + dp(k)*476.*do3b(k)
    do3b(k) = 0.0_8
    wa_devb(k) = wa_devb(k) + dp(k)*1.02*dh2ob(k)
    dh2ob(k) = 0.0_8
    ta_devb(k) = ta_devb(k) + dtb(k)
    dtb(k) = 0.0_8
  END DO
END SUBROUTINE IRRAD_B

!  Differentiation of planck in reverse (adjoint) mode:
!   gradient     of useful results: t xlayer
!   with respect to varying inputs: t
!***********************************************************************
SUBROUTINE PLANCK_B(ibn, cb, t, tb, xlayer, xlayerb)
  IMPLICIT NONE
! spectral band index
  INTEGER, INTENT(IN) :: ibn
! Planck table coefficients
  REAL*8, INTENT(IN) :: cb(6, 10)
! temperature (K)
  REAL*8, INTENT(IN) :: t
  REAL*8 :: tb
! planck flux (w/m2)
  REAL*8 :: xlayer
  REAL*8 :: xlayerb
  REAL*8 :: temp1
  REAL*8 :: temp0
  REAL*8 :: tempb
  REAL*8 :: temp
  temp1 = cb(5, ibn) + cb(6, ibn)*t
  temp0 = cb(4, ibn) + t*temp1
  temp = cb(3, ibn) + t*temp0
  tempb = t**2*xlayerb
  tb = tb + (t**2*cb(6, ibn)+t*temp1+temp0)*tempb + (2*(t*temp)+cb(2, &
&   ibn))*xlayerb
END SUBROUTINE PLANCK_B

!  Differentiation of h2oexps in reverse (adjoint) mode:
!   gradient     of useful results: dh2o dt h2oexp
!   with respect to varying inputs: dh2o dt
!**********************************************************************
SUBROUTINE H2OEXPS_B(ib, np, dh2o, dh2ob, pa, dt, dtb, xkw, aw, bw, pm, &
& mw, h2oexp, h2oexpb)
  IMPLICIT NONE
  INTEGER :: ib, np, ik, k
!---- input parameters ------
  REAL*8 :: dh2o(0:np), pa(0:np), dt(0:np)
  REAL*8 :: dh2ob(0:np), dtb(0:np)
!---- output parameters -----
  REAL*8 :: h2oexp(0:np, 6)
  REAL*8 :: h2oexpb(0:np, 6)
!---- static data -----
  INTEGER :: mw(9)
  REAL*8 :: xkw(9), aw(9), bw(9), pm(9)
!---- temporary arrays -----
  REAL*8 :: xh
  REAL*8 :: xhb
  INTRINSIC EXP
  INTEGER :: branch
  REAL*8 :: tempb
  REAL*8 :: temp
!**********************************************************************
!    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
!    and bw,  therefore, h2oexp for these sub-bands are identical.
!**********************************************************************
  DO k=0,np
!-----xh is the scaled water vapor amount for line absorption
!     computed from Eq. (4.4).
    CALL PUSHREAL8(xh)
    xh = dh2o(k)*(pa(k)/500.)**pm(ib)*(1.+(aw(ib)+bw(ib)*dt(k))*dt(k))
!-----h2oexp is the water vapor transmittance of the layer k
!     due to line absorption
    h2oexp(k, 1) = EXP(-(xh*xkw(ib)))
!-----compute transmittances from Eq. (8.22)
    DO ik=2,6
      IF (mw(ib) .EQ. 6) THEN
        CALL PUSHREAL8(xh)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        CALL PUSHREAL8(h2oexp(k, ik))
        h2oexp(k, ik) = xh*xh*xh
        CALL PUSHCONTROL2B(3)
      ELSE IF (mw(ib) .EQ. 8) THEN
        CALL PUSHREAL8(xh)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        CALL PUSHREAL8(xh)
        xh = xh*xh
        CALL PUSHREAL8(h2oexp(k, ik))
        h2oexp(k, ik) = xh*xh
        CALL PUSHCONTROL2B(2)
      ELSE IF (mw(ib) .EQ. 9) THEN
        CALL PUSHREAL8(xh)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)*h2oexp(k, ik-1)
        CALL PUSHREAL8(h2oexp(k, ik))
        h2oexp(k, ik) = xh*xh*xh
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHREAL8(xh)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        CALL PUSHREAL8(xh)
        xh = xh*xh
        CALL PUSHREAL8(xh)
        xh = xh*xh
        CALL PUSHREAL8(h2oexp(k, ik))
        h2oexp(k, ik) = xh*xh
        CALL PUSHCONTROL2B(0)
      END IF
    END DO
  END DO
  DO k=np,0,-1
    DO ik=6,2,-1
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(h2oexp(k, ik))
          xhb = 2*xh*h2oexpb(k, ik)
          h2oexpb(k, ik) = 0.0_8
          CALL POPREAL8(xh)
          xhb = 2*xh*xhb
          CALL POPREAL8(xh)
          xhb = 2*xh*xhb
          CALL POPREAL8(xh)
          h2oexpb(k, ik-1) = h2oexpb(k, ik-1) + 2*h2oexp(k, ik-1)*xhb
        ELSE
          CALL POPREAL8(h2oexp(k, ik))
          xhb = 3*xh**2*h2oexpb(k, ik)
          h2oexpb(k, ik) = 0.0_8
          CALL POPREAL8(xh)
          h2oexpb(k, ik-1) = h2oexpb(k, ik-1) + 3*h2oexp(k, ik-1)**2*xhb
        END IF
      ELSE IF (branch .EQ. 2) THEN
        CALL POPREAL8(h2oexp(k, ik))
        xhb = 2*xh*h2oexpb(k, ik)
        h2oexpb(k, ik) = 0.0_8
        CALL POPREAL8(xh)
        xhb = 2*xh*xhb
        CALL POPREAL8(xh)
        h2oexpb(k, ik-1) = h2oexpb(k, ik-1) + 2*h2oexp(k, ik-1)*xhb
      ELSE
        CALL POPREAL8(h2oexp(k, ik))
        xhb = 3*xh**2*h2oexpb(k, ik)
        h2oexpb(k, ik) = 0.0_8
        CALL POPREAL8(xh)
        h2oexpb(k, ik-1) = h2oexpb(k, ik-1) + 2*h2oexp(k, ik-1)*xhb
      END IF
    END DO
    xhb = -(EXP(-(xkw(ib)*xh))*xkw(ib)*h2oexpb(k, 1))
    h2oexpb(k, 1) = 0.0_8
    CALL POPREAL8(xh)
    temp = aw(ib) + bw(ib)*dt(k)
    tempb = (pa(k)/500.)**pm(ib)*xhb
    dh2ob(k) = dh2ob(k) + (temp*dt(k)+1.)*tempb
    dtb(k) = dtb(k) + (dh2o(k)*temp+dt(k)*dh2o(k)*bw(ib))*tempb
  END DO
END SUBROUTINE H2OEXPS_B

!  Differentiation of conexps in reverse (adjoint) mode:
!   gradient     of useful results: dcont conexp
!   with respect to varying inputs: dcont conexp
!**********************************************************************
SUBROUTINE CONEXPS_B(ib, np, dcont, dcontb, xke, conexp, conexpb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters ------
  REAL*8 :: dcont(0:np)
  REAL*8 :: dcontb(0:np)
!---- updated parameters -----
  REAL*8 :: conexp(0:np, 3)
  REAL*8 :: conexpb(0:np, 3)
!---- static data -----
  REAL*8 :: xke(9)
  INTRINSIC EXP
  INTEGER :: branch
!****************************************************************
  DO k=0,np
    conexp(k, 1) = EXP(-(dcont(k)*xke(ib)))
!-----The absorption coefficients for sub-bands 3b and 3a are, respectively,
!     two and four times the absorption coefficient for sub-band 3c (Table 9).
!     Note that conexp(3) is for sub-band 3a. 
    IF (ib .EQ. 3) THEN
      CALL PUSHREAL8(conexp(k, 2))
      conexp(k, 2) = conexp(k, 1)*conexp(k, 1)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=np,0,-1
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      conexpb(k, 2) = conexpb(k, 2) + 2*conexp(k, 2)*conexpb(k, 3)
      conexpb(k, 3) = 0.0_8
      CALL POPREAL8(conexp(k, 2))
      conexpb(k, 1) = conexpb(k, 1) + 2*conexp(k, 1)*conexpb(k, 2)
      conexpb(k, 2) = 0.0_8
    END IF
    dcontb(k) = dcontb(k) - EXP(-(xke(ib)*dcont(k)))*xke(ib)*conexpb(k, &
&     1)
    conexpb(k, 1) = 0.0_8
  END DO
END SUBROUTINE CONEXPS_B

!  Differentiation of n2oexps in reverse (adjoint) mode:
!   gradient     of useful results: dt n2oexp
!   with respect to varying inputs: dt n2oexp
!**********************************************************************
SUBROUTINE N2OEXPS_B(ib, np, dn2o, pa, dt, dtb, n2oexp, n2oexpb)
  IMPLICIT NONE
  INTEGER :: ib, k, np
!---- input parameters -----
  REAL*8 :: dn2o(0:np), pa(0:np), dt(0:np)
  REAL*8 :: dtb(0:np)
!---- output parameters -----
  REAL*8 :: n2oexp(0:np, 4)
  REAL*8 :: n2oexpb(0:np, 4)
!---- temporary arrays -----
  REAL*8 :: xc, xc1, xc2
  REAL*8 :: xcb, xc1b, xc2b
  INTRINSIC EXP
  INTEGER :: branch
  REAL*8 :: tempb
!-----Scaling and absorption data are given in Table 5.
!     Transmittances are computed using Eqs. (8.21) and (8.22).
  DO k=0,np
!-----four exponential by powers of 21 for band 6.
    IF (ib .EQ. 6) THEN
      CALL PUSHREAL8(xc)
      xc = dn2o(k)*(1.+(1.9297e-3+4.3750e-6*dt(k))*dt(k))
      n2oexp(k, 1) = EXP(-(xc*6.31582e-2))
      xc = n2oexp(k, 1)*n2oexp(k, 1)*n2oexp(k, 1)
!-----four exponential by powers of 8 for band 7
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL8(xc)
      xc = dn2o(k)*(pa(k)/500.0)**0.48*(1.+(1.3804e-3+7.4838e-6*dt(k))*&
&       dt(k))
      n2oexp(k, 1) = EXP(-(xc*5.35779e-2))
      xc = n2oexp(k, 1)*n2oexp(k, 1)
      CALL PUSHREAL8(xc)
      xc = xc*xc
      CALL PUSHREAL8(n2oexp(k, 2))
      n2oexp(k, 2) = xc*xc
      CALL PUSHREAL8(xc)
      xc = n2oexp(k, 2)*n2oexp(k, 2)
      CALL PUSHREAL8(xc)
      xc = xc*xc
      CALL PUSHREAL8(n2oexp(k, 3))
      n2oexp(k, 3) = xc*xc
      CALL PUSHREAL8(xc)
      xc = n2oexp(k, 3)*n2oexp(k, 3)
      CALL PUSHREAL8(xc)
      xc = xc*xc
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=np,0,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xcb = 2*xc*n2oexpb(k, 4)
      n2oexpb(k, 4) = 0.0_8
      CALL POPREAL8(xc)
      xcb = 2*xc*xcb
      CALL POPREAL8(xc)
      n2oexpb(k, 3) = n2oexpb(k, 3) + 2*n2oexp(k, 3)*xcb
      CALL POPREAL8(n2oexp(k, 3))
      xcb = 2*xc*n2oexpb(k, 3)
      n2oexpb(k, 3) = 0.0_8
      CALL POPREAL8(xc)
      xcb = 2*xc*xcb
      CALL POPREAL8(xc)
      n2oexpb(k, 2) = n2oexpb(k, 2) + 2*n2oexp(k, 2)*xcb
      CALL POPREAL8(n2oexp(k, 2))
      xcb = 2*xc*n2oexpb(k, 2)
      n2oexpb(k, 2) = 0.0_8
      CALL POPREAL8(xc)
      xcb = 2*xc*xcb
      n2oexpb(k, 1) = n2oexpb(k, 1) + 2*n2oexp(k, 1)*xcb
      xc = dn2o(k)*(pa(k)/500.0)**0.48*(1.+(1.3804e-3+7.4838e-6*dt(k))*&
&       dt(k))
      xcb = -(EXP(-(5.35779e-2*xc))*5.35779e-2*n2oexpb(k, 1))
      n2oexpb(k, 1) = 0.0_8
      CALL POPREAL8(xc)
      tempb = (pa(k)/500.0)**0.48*dn2o(k)*xcb
      dtb(k) = dtb(k) + (7.4838e-6*dt(k)+dt(k)*7.4838e-6+1.3804e-3)*&
&       tempb
    ELSE
      xc1 = xc*xc
      xc2 = xc1*xc1
      xc2b = xc*xc1*n2oexpb(k, 2)
      xc1b = 2*xc1*xc2b + xc2*xc*n2oexpb(k, 2)
      xcb = 2*xc*xc1b + xc2*xc1*n2oexpb(k, 2)
      n2oexpb(k, 2) = 0.0_8
      n2oexpb(k, 1) = n2oexpb(k, 1) + 3*n2oexp(k, 1)**2*xcb
      xc = dn2o(k)*(1.+(1.9297e-3+4.3750e-6*dt(k))*dt(k))
      xcb = -(EXP(-(6.31582e-2*xc))*6.31582e-2*n2oexpb(k, 1))
      n2oexpb(k, 1) = 0.0_8
      CALL POPREAL8(xc)
      dtb(k) = dtb(k) + (dn2o(k)*(4.3750e-6*dt(k)+1.9297e-3)+dt(k)*dn2o(&
&       k)*4.3750e-6)*xcb
    END IF
  END DO
END SUBROUTINE N2OEXPS_B

!  Differentiation of ch4exps in reverse (adjoint) mode:
!   gradient     of useful results: dt ch4exp
!   with respect to varying inputs: dt ch4exp
!**********************************************************************
SUBROUTINE CH4EXPS_B(ib, np, dch4, pa, dt, dtb, ch4exp, ch4expb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: dch4(0:np), pa(0:np), dt(0:np)
  REAL*8 :: dtb(0:np)
!---- output parameters -----
  REAL*8 :: ch4exp(0:np, 4)
  REAL*8 :: ch4expb(0:np, 4)
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcb
  INTRINSIC EXP
  INTEGER :: branch
  REAL*8 :: tempb
!*****  Scaling and absorption data are given in Table 5  *****
  DO k=0,np
!-----four exponentials for band 6
    IF (ib .EQ. 6) THEN
      CALL PUSHREAL8(xc)
!-----four exponentials by powers of 12 for band 7
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL8(xc)
      xc = dch4(k)*(pa(k)/500.0)**0.65*(1.+(5.9590e-4-2.2931e-6*dt(k))*&
&       dt(k))
      ch4exp(k, 1) = EXP(-(xc*6.29247e-2))
      xc = ch4exp(k, 1)*ch4exp(k, 1)*ch4exp(k, 1)
      CALL PUSHREAL8(xc)
      xc = xc*xc
      CALL PUSHREAL8(ch4exp(k, 2))
      ch4exp(k, 2) = xc*xc
      CALL PUSHREAL8(xc)
      xc = ch4exp(k, 2)*ch4exp(k, 2)*ch4exp(k, 2)
      CALL PUSHREAL8(xc)
      xc = xc*xc
      CALL PUSHREAL8(ch4exp(k, 3))
      ch4exp(k, 3) = xc*xc
      CALL PUSHREAL8(xc)
      xc = ch4exp(k, 3)*ch4exp(k, 3)*ch4exp(k, 3)
      CALL PUSHREAL8(xc)
      xc = xc*xc
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=np,0,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xcb = 2*xc*ch4expb(k, 4)
      ch4expb(k, 4) = 0.0_8
      CALL POPREAL8(xc)
      xcb = 2*xc*xcb
      CALL POPREAL8(xc)
      ch4expb(k, 3) = ch4expb(k, 3) + 3*ch4exp(k, 3)**2*xcb
      CALL POPREAL8(ch4exp(k, 3))
      xcb = 2*xc*ch4expb(k, 3)
      ch4expb(k, 3) = 0.0_8
      CALL POPREAL8(xc)
      xcb = 2*xc*xcb
      CALL POPREAL8(xc)
      ch4expb(k, 2) = ch4expb(k, 2) + 3*ch4exp(k, 2)**2*xcb
      CALL POPREAL8(ch4exp(k, 2))
      xcb = 2*xc*ch4expb(k, 2)
      ch4expb(k, 2) = 0.0_8
      CALL POPREAL8(xc)
      xcb = 2*xc*xcb
      ch4expb(k, 1) = ch4expb(k, 1) + 3*ch4exp(k, 1)**2*xcb
      xc = dch4(k)*(pa(k)/500.0)**0.65*(1.+(5.9590e-4-2.2931e-6*dt(k))*&
&       dt(k))
      xcb = -(EXP(-(6.29247e-2*xc))*6.29247e-2*ch4expb(k, 1))
      ch4expb(k, 1) = 0.0_8
      CALL POPREAL8(xc)
      tempb = (pa(k)/500.0)**0.65*dch4(k)*xcb
      dtb(k) = dtb(k) + (5.9590e-4-dt(k)*2.2931e-6-2.2931e-6*dt(k))*&
&       tempb
    ELSE
      xc = dch4(k)*(1.+(1.7007e-2+1.5826e-4*dt(k))*dt(k))
      xcb = -(EXP(-(5.80708e-3*xc))*5.80708e-3*ch4expb(k, 1))
      ch4expb(k, 1) = 0.0_8
      CALL POPREAL8(xc)
      dtb(k) = dtb(k) + (dch4(k)*(1.5826e-4*dt(k)+1.7007e-2)+dt(k)*dch4(&
&       k)*1.5826e-4)*xcb
    END IF
  END DO
END SUBROUTINE CH4EXPS_B

!  Differentiation of comexps in reverse (adjoint) mode:
!   gradient     of useful results: dt comexp
!   with respect to varying inputs: dt comexp
!**********************************************************************
SUBROUTINE COMEXPS_B(ib, np, dcom, dt, dtb, comexp, comexpb)
  IMPLICIT NONE
  INTEGER :: ib, ik, np, k
!---- input parameters -----
  REAL*8 :: dcom(0:np), dt(0:np)
  REAL*8 :: dtb(0:np)
!---- output parameters -----
  REAL*8 :: comexp(0:np, 6)
  REAL*8 :: comexpb(0:np, 6)
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcb
  INTRINSIC EXP
  INTEGER :: branch
!*****  Scaling and absorpton data are given in Table 6  *****
  DO k=0,np
    IF (ib .EQ. 4) THEN
      CALL PUSHREAL8(xc)
      xc = dcom(k)*(1.+(3.5775e-2+4.0447e-4*dt(k))*dt(k))
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (ib .EQ. 5) THEN
      CALL PUSHREAL8(xc)
      xc = dcom(k)*(1.+(3.4268e-2+3.7401e-4*dt(k))*dt(k))
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    comexp(k, 1) = EXP(-(xc*1.922e-7))
    DO ik=2,6
      CALL PUSHREAL8(xc)
      xc = comexp(k, ik-1)*comexp(k, ik-1)
      CALL PUSHREAL8(xc)
      xc = xc*xc
      CALL PUSHREAL8(comexp(k, ik))
      comexp(k, ik) = xc*comexp(k, ik-1)
    END DO
  END DO
  xcb = 0.0_8
  DO k=np,0,-1
    DO ik=6,2,-1
      CALL POPREAL8(comexp(k, ik))
      xcb = xcb + comexp(k, ik-1)*comexpb(k, ik)
      comexpb(k, ik-1) = comexpb(k, ik-1) + xc*comexpb(k, ik)
      comexpb(k, ik) = 0.0_8
      CALL POPREAL8(xc)
      xcb = 2*xc*xcb
      CALL POPREAL8(xc)
      comexpb(k, ik-1) = comexpb(k, ik-1) + 2*comexp(k, ik-1)*xcb
      xcb = 0.0_8
    END DO
    xcb = xcb - EXP(-(1.922e-7*xc))*1.922e-7*comexpb(k, 1)
    comexpb(k, 1) = 0.0_8
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(xc)
      dtb(k) = dtb(k) + (dcom(k)*(3.7401e-4*dt(k)+3.4268e-2)+dt(k)*dcom(&
&       k)*3.7401e-4)*xcb
      xcb = 0.0_8
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(xc)
      dtb(k) = dtb(k) + (dcom(k)*(4.0447e-4*dt(k)+3.5775e-2)+dt(k)*dcom(&
&       k)*4.0447e-4)*xcb
      xcb = 0.0_8
    END IF
  END DO
END SUBROUTINE COMEXPS_B

!  Differentiation of cfcexps in reverse (adjoint) mode:
!   gradient     of useful results: dt cfcexp
!   with respect to varying inputs: dt cfcexp
!**********************************************************************
SUBROUTINE CFCEXPS_B(ib, np, a1, b1, fk1, a2, b2, fk2, dcfc, dt, dtb, &
& cfcexp, cfcexpb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: dcfc(0:np), dt(0:np)
  REAL*8 :: dtb(0:np)
!---- output parameters -----
  REAL*8 :: cfcexp(0:np)
  REAL*8 :: cfcexpb(0:np)
!---- static data -----
  REAL*8 :: a1, b1, fk1, a2, b2, fk2
!---- temporary arrays -----
  REAL*8 :: xf
  REAL*8 :: xfb
  INTRINSIC EXP
  INTEGER :: branch
!**********************************************************************
  DO k=0,np
!-----compute the scaled cfc amount (xf) and exponential (cfcexp)
    IF (ib .EQ. 4) THEN
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=np,0,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xf = dcfc(k)*(1.+(a2+b2*dt(k))*dt(k))
      xfb = -(EXP(-(fk2*xf))*fk2*cfcexpb(k))
      cfcexpb(k) = 0.0_8
      dtb(k) = dtb(k) + (dcfc(k)*(a2+b2*dt(k))+dt(k)*dcfc(k)*b2)*xfb
    ELSE
      xf = dcfc(k)*(1.+(a1+b1*dt(k))*dt(k))
      xfb = -(EXP(-(fk1*xf))*fk1*cfcexpb(k))
      cfcexpb(k) = 0.0_8
      dtb(k) = dtb(k) + (dcfc(k)*(a1+b1*dt(k))+dt(k)*dcfc(k)*b1)*xfb
    END IF
  END DO
END SUBROUTINE CFCEXPS_B

!  Differentiation of b10exps in reverse (adjoint) mode:
!   gradient     of useful results: dh2o dcont dt co2exp h2oexp
!                n2oexp conexp
!   with respect to varying inputs: dh2o dcont dt co2exp h2oexp
!                n2oexp conexp
!**********************************************************************
SUBROUTINE B10EXPS_B(np, dh2o, dh2ob, dcont, dcontb, dco2, dn2o, pa, dt&
& , dtb, h2oexp, h2oexpb, conexp, conexpb, co2exp, co2expb, n2oexp, &
& n2oexpb)
  IMPLICIT NONE
  INTEGER :: np, k
!---- input parameters -----
  REAL*8 :: dh2o(0:np), dcont(0:np), dn2o(0:np)
  REAL*8 :: dh2ob(0:np), dcontb(0:np)
  REAL*8 :: dco2(0:np), pa(0:np), dt(0:np)
  REAL*8 :: dtb(0:np)
!---- output parameters -----
  REAL*8 :: h2oexp(0:np, 5), conexp(0:np), co2exp(0:np, 6), n2oexp(0:np&
& , 2)
  REAL*8 :: h2oexpb(0:np, 5), conexpb(0:np), co2expb(0:np, 6), n2oexpb(0&
& :np, 2)
!---- temporary arrays -----
  REAL*8 :: xx, xx1, xx2, xx3
  REAL*8 :: xxb, xx1b, xx2b, xx3b
  INTRINSIC EXP
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  REAL*8 :: tempb
!**********************************************************************
  DO k=0,np
!-----Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3).
    CALL PUSHREAL8(xx)
    xx = dh2o(k)*(pa(k)/500.0)*(1.+(0.0149+6.20e-5*dt(k))*dt(k))
!-----six exponentials by powers of 8
    h2oexp(k, 1) = EXP(-(xx*0.10624))
    xx = h2oexp(k, 1)*h2oexp(k, 1)
    CALL PUSHREAL8(xx)
    xx = xx*xx
    CALL PUSHREAL8(h2oexp(k, 2))
    h2oexp(k, 2) = xx*xx
    CALL PUSHREAL8(xx)
    xx = h2oexp(k, 2)*h2oexp(k, 2)
    CALL PUSHREAL8(xx)
    xx = xx*xx
    CALL PUSHREAL8(h2oexp(k, 3))
    h2oexp(k, 3) = xx*xx
    CALL PUSHREAL8(xx)
    xx = h2oexp(k, 3)*h2oexp(k, 3)
    CALL PUSHREAL8(xx)
    xx = xx*xx
    CALL PUSHREAL8(h2oexp(k, 4))
    h2oexp(k, 4) = xx*xx
    CALL PUSHREAL8(xx)
    xx = h2oexp(k, 4)*h2oexp(k, 4)
    CALL PUSHREAL8(xx)
    xx = xx*xx
!-----one exponential of h2o continuum for sub-band 3a (Table 9).
!-----Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).
    CALL PUSHREAL8(xx)
    xx = dco2(k)*(pa(k)/300.0)**0.5*(1.+(0.0179+1.02e-4*dt(k))*dt(k))
!-----six exponentials by powers of 8
    co2exp(k, 1) = EXP(-(xx*2.656e-5))
    xx = co2exp(k, 1)*co2exp(k, 1)
    CALL PUSHREAL8(xx)
    xx = xx*xx
    CALL PUSHREAL8(co2exp(k, 2))
    co2exp(k, 2) = xx*xx
    CALL PUSHREAL8(xx)
    xx = co2exp(k, 2)*co2exp(k, 2)
    CALL PUSHREAL8(xx)
    xx = xx*xx
    CALL PUSHREAL8(co2exp(k, 3))
    co2exp(k, 3) = xx*xx
    CALL PUSHREAL8(xx)
    xx = co2exp(k, 3)*co2exp(k, 3)
    CALL PUSHREAL8(xx)
    xx = xx*xx
    CALL PUSHREAL8(co2exp(k, 4))
    co2exp(k, 4) = xx*xx
    CALL PUSHREAL8(xx)
    xx = co2exp(k, 4)*co2exp(k, 4)
    CALL PUSHREAL8(xx)
    xx = xx*xx
    CALL PUSHREAL8(co2exp(k, 5))
    co2exp(k, 5) = xx*xx
    CALL PUSHREAL8(xx)
    xx = co2exp(k, 5)*co2exp(k, 5)
    CALL PUSHREAL8(xx)
    xx = xx*xx
!-----Compute the scaled n2o amount for Band 10 (Table 5).
    CALL PUSHREAL8(xx)
    xx = dn2o(k)*(1.+(1.4476e-3+3.6656e-6*dt(k))*dt(k))
!-----Two exponentials by powers of 58
    n2oexp(k, 1) = EXP(-(xx*0.25238))
    xx = n2oexp(k, 1)*n2oexp(k, 1)
    CALL PUSHREAL8(xx1)
    xx1 = xx*xx
    CALL PUSHREAL8(xx1)
    xx1 = xx1*xx1
  END DO
  DO k=np,0,-1
    xx2 = xx1*xx1
    xx3 = xx2*xx2
    tempb = xx2*xx3*n2oexpb(k, 2)
    tempb0 = xx*xx1*n2oexpb(k, 2)
    xxb = xx1*tempb
    xx3b = xx2*tempb0
    xx2b = 2*xx2*xx3b + xx3*tempb0
    xx1b = 2*xx1*xx2b + xx*tempb
    n2oexpb(k, 2) = 0.0_8
    CALL POPREAL8(xx1)
    xx1b = 2*xx1*xx1b
    CALL POPREAL8(xx1)
    xxb = xxb + 2*xx*xx1b
    n2oexpb(k, 1) = n2oexpb(k, 1) + 2*n2oexp(k, 1)*xxb
    xx = dn2o(k)*(1.+(1.4476e-3+3.6656e-6*dt(k))*dt(k))
    xxb = -(EXP(-(0.25238*xx))*0.25238*n2oexpb(k, 1))
    n2oexpb(k, 1) = 0.0_8
    CALL POPREAL8(xx)
    dtb(k) = dtb(k) + (dn2o(k)*(3.6656e-6*dt(k)+1.4476e-3)+dt(k)*dn2o(k)&
&     *3.6656e-6)*xxb
    xxb = 2*xx*co2expb(k, 6)
    co2expb(k, 6) = 0.0_8
    CALL POPREAL8(xx)
    xxb = 2*xx*xxb
    CALL POPREAL8(xx)
    co2expb(k, 5) = co2expb(k, 5) + 2*co2exp(k, 5)*xxb
    CALL POPREAL8(co2exp(k, 5))
    xxb = 2*xx*co2expb(k, 5)
    co2expb(k, 5) = 0.0_8
    CALL POPREAL8(xx)
    xxb = 2*xx*xxb
    CALL POPREAL8(xx)
    co2expb(k, 4) = co2expb(k, 4) + 2*co2exp(k, 4)*xxb
    CALL POPREAL8(co2exp(k, 4))
    xxb = 2*xx*co2expb(k, 4)
    co2expb(k, 4) = 0.0_8
    CALL POPREAL8(xx)
    xxb = 2*xx*xxb
    CALL POPREAL8(xx)
    co2expb(k, 3) = co2expb(k, 3) + 2*co2exp(k, 3)*xxb
    CALL POPREAL8(co2exp(k, 3))
    xxb = 2*xx*co2expb(k, 3)
    co2expb(k, 3) = 0.0_8
    CALL POPREAL8(xx)
    xxb = 2*xx*xxb
    CALL POPREAL8(xx)
    co2expb(k, 2) = co2expb(k, 2) + 2*co2exp(k, 2)*xxb
    CALL POPREAL8(co2exp(k, 2))
    xxb = 2*xx*co2expb(k, 2)
    co2expb(k, 2) = 0.0_8
    CALL POPREAL8(xx)
    xxb = 2*xx*xxb
    co2expb(k, 1) = co2expb(k, 1) + 2*co2exp(k, 1)*xxb
    xx = dco2(k)*(pa(k)/300.0)**0.5*(1.+(0.0179+1.02e-4*dt(k))*dt(k))
    xxb = -(EXP(-(2.656e-5*xx))*2.656e-5*co2expb(k, 1))
    co2expb(k, 1) = 0.0_8
    CALL POPREAL8(xx)
    tempb1 = (pa(k)/300.0)**0.5*dco2(k)*xxb
    dcontb(k) = dcontb(k) - EXP(-(109.0*dcont(k)))*109.0*conexpb(k)
    conexpb(k) = 0.0_8
    xxb = 2*xx*h2oexpb(k, 5)
    h2oexpb(k, 5) = 0.0_8
    CALL POPREAL8(xx)
    xxb = 2*xx*xxb
    CALL POPREAL8(xx)
    h2oexpb(k, 4) = h2oexpb(k, 4) + 2*h2oexp(k, 4)*xxb
    CALL POPREAL8(h2oexp(k, 4))
    xxb = 2*xx*h2oexpb(k, 4)
    h2oexpb(k, 4) = 0.0_8
    CALL POPREAL8(xx)
    xxb = 2*xx*xxb
    CALL POPREAL8(xx)
    h2oexpb(k, 3) = h2oexpb(k, 3) + 2*h2oexp(k, 3)*xxb
    CALL POPREAL8(h2oexp(k, 3))
    xxb = 2*xx*h2oexpb(k, 3)
    h2oexpb(k, 3) = 0.0_8
    CALL POPREAL8(xx)
    xxb = 2*xx*xxb
    CALL POPREAL8(xx)
    h2oexpb(k, 2) = h2oexpb(k, 2) + 2*h2oexp(k, 2)*xxb
    CALL POPREAL8(h2oexp(k, 2))
    xxb = 2*xx*h2oexpb(k, 2)
    h2oexpb(k, 2) = 0.0_8
    CALL POPREAL8(xx)
    xxb = 2*xx*xxb
    h2oexpb(k, 1) = h2oexpb(k, 1) + 2*h2oexp(k, 1)*xxb
    xx = dh2o(k)*(pa(k)/500.0)*(1.+(0.0149+6.20e-5*dt(k))*dt(k))
    xxb = -(EXP(-(0.10624*xx))*0.10624*h2oexpb(k, 1))
    h2oexpb(k, 1) = 0.0_8
    CALL POPREAL8(xx)
    tempb2 = pa(k)*dh2o(k)*xxb/500.0
    dtb(k) = dtb(k) + (6.20e-5*dt(k)+dt(k)*6.20e-5+0.0149)*tempb2 + (&
&     1.02e-4*dt(k)+dt(k)*1.02e-4+0.0179)*tempb1
    dh2ob(k) = dh2ob(k) + ((6.20e-5*dt(k)+0.0149)*dt(k)+1.)*pa(k)*xxb/&
&     500.0
  END DO
END SUBROUTINE B10EXPS_B

!  Differentiation of tablup in reverse (adjoint) mode:
!   gradient     of useful results: s1 s2 s3 dt dw tran
!   with respect to varying inputs: s1 s2 s3 dt dw tran
!**********************************************************************
SUBROUTINE TABLUP_B(nx1, nh1, dw, dwb, p, dt, dtb, s1, s1b, s2, s2b, s3&
& , s3b, w1, p1, dwe, dpe, coef1, coef2, coef3, tran, tranb)
  IMPLICIT NONE
  INTEGER :: nx1, nh1
!---- input parameters -----
  REAL*8 :: w1, p1, dwe, dpe
  REAL*8 :: dw, p, dt
  REAL*8 :: dwb, dtb
  REAL*8 :: coef1(nx1, nh1), coef2(nx1, nh1), coef3(nx1, nh1)
!---- update parameter -----
  REAL*8 :: s1, s2, s3, tran
  REAL*8 :: s1b, s2b, s3b, tranb
!---- temporary variables -----
  REAL*8 :: we, pe, fw, fp, pa, pb, pc, ax, ba, bb, t1, ca, cb, t2
  REAL*8 :: web, peb, fwb, fpb, pab, pbb, pcb, axb, bab, bbb, t1b, cab, &
& cbb, t2b
  REAL*8 :: x1, x2, x3, xx, x1c
  REAL*8 :: x1b, x2b, x3b, xxb, x1cb
  INTEGER :: iw, ip
  INTRINSIC LOG10
  INTRINSIC REAL
  INTRINSIC MIN
  INTRINSIC INT
  INTRINSIC MAX
  INTEGER :: branch
  REAL*8 :: tempb0
  REAL*8 :: tempb
  REAL*8 :: y2
  REAL*8 :: y1
!-----Compute effective pressure (x2) and temperature (x3) following 
!     Eqs. (8.28) and (8.29)
  s1 = s1 + dw
  s2 = s2 + p*dw
  s3 = s3 + dt*dw
  x1 = s1
  x1c = 1.0/s1
  x2 = s2*x1c
  x3 = s3*x1c
!-----normalize we and pe
!       we=(log10(x1)-w1)/dwe
!       pe=(log10(x2)-p1)/dpe
  we = (LOG10(x1)-w1)*dwe
  pe = (LOG10(x2)-p1)*dpe
  y1 = REAL(nh1 - 1)
  IF (we .GT. y1) THEN
    we = y1
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    we = we
  END IF
  y2 = REAL(nx1 - 1)
  IF (pe .GT. y2) THEN
    pe = y2
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
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
  fp = pe - REAL(ip - 1)
!-----linear interpolation in pressure
  pa = coef1(ip, iw-1) + (coef1(ip+1, iw-1)-coef1(ip, iw-1))*fp
  pb = coef1(ip, iw) + (coef1(ip+1, iw)-coef1(ip, iw))*fp
  pc = coef1(ip, iw+1) + (coef1(ip+1, iw+1)-coef1(ip, iw+1))*fp
!-----quadratic interpolation in absorber amount for coef1
  ax = ((pc+pa)*fw+(pc-pa))*fw*0.5 + pb*(1.-fw*fw)
!-----linear interpolation in absorber amount for coef2 and coef3
  ba = coef2(ip, iw) + (coef2(ip+1, iw)-coef2(ip, iw))*fp
  bb = coef2(ip, iw+1) + (coef2(ip+1, iw+1)-coef2(ip, iw+1))*fp
  t1 = ba + (bb-ba)*fw
  ca = coef3(ip, iw) + (coef3(ip+1, iw)-coef3(ip, iw))*fp
  cb = coef3(ip, iw+1) + (coef3(ip+1, iw+1)-coef3(ip, iw+1))*fp
  t2 = ca + (cb-ca)*fw
!-----update the total transmittance between levels k1 and k2
  xx = ax + (t1+t2*x3)*x3
  IF (xx .GT. 0.9999999) THEN
    xx = 0.9999999
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    xx = xx
  END IF
  IF (xx .LT. 0.0000001) THEN
    xx = 0.0000001
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    xx = xx
  END IF
  xxb = tran*tranb
  tranb = xx*tranb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) xxb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) xxb = 0.0_8
  tempb = x3*xxb
  axb = xxb
  t1b = tempb
  t2b = x3*tempb
  x3b = (t1+t2*x3)*xxb + t2*tempb
  cab = (1.0_8-fw)*t2b
  cbb = fw*t2b
  bab = (1.0_8-fw)*t1b
  bbb = fw*t1b
  tempb0 = 0.5*fw*axb
  fwb = (bb-ba)*t1b + (0.5*((pc+pa)*fw+pc-pa)-pb*2*fw)*axb + (pc+pa)*&
&   tempb0 + (cb-ca)*t2b
  pcb = (fw+1.0_8)*tempb0
  pab = (fw-1.0)*tempb0
  pbb = (1.-fw**2)*axb
  fpb = (coef3(ip+1, iw)-coef3(ip, iw))*cab + (coef2(ip+1, iw)-coef2(ip&
&   , iw))*bab + (coef1(ip+1, iw)-coef1(ip, iw))*pbb + (coef1(ip+1, iw-1&
&   )-coef1(ip, iw-1))*pab + (coef1(ip+1, iw+1)-coef1(ip, iw+1))*pcb + (&
&   coef2(ip+1, iw+1)-coef2(ip, iw+1))*bbb + (coef3(ip+1, iw+1)-coef3(ip&
&   , iw+1))*cbb
  peb = fpb
  web = fwb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) peb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) web = 0.0_8
  x2b = dpe*peb/(x2*LOG(10.0))
  x1b = dwe*web/(x1*LOG(10.0))
  s3b = s3b + x1c*x3b
  x1cb = s2*x2b + s3*x3b
  s2b = s2b + x1c*x2b
  s1b = s1b + x1b - x1cb/s1**2
  dtb = dtb + dw*s3b
  dwb = dwb + p*s2b + s1b + dt*s3b
END SUBROUTINE TABLUP_B

!  Differentiation of h2okdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tcon h2oexp th2o conexp
!   with respect to varying inputs: tcon h2oexp th2o conexp
!**********************************************************************
SUBROUTINE H2OKDIS_B(ib, np, k, fkw, gkw, ne, h2oexp, h2oexpb, conexp, &
& conexpb, th2o, th2ob, tcon, tconb, tran, tranb)
  IMPLICIT NONE
!---- input parameters ------
  INTEGER :: ib, ne, np, k
  REAL*8 :: h2oexp(0:np, 6), conexp(0:np, 3)
  REAL*8 :: h2oexpb(0:np, 6), conexpb(0:np, 3)
  REAL*8 :: fkw(6, 9), gkw(6, 3)
!---- updated parameters -----
  REAL*8 :: th2o(6), tcon(3), tran
  REAL*8 :: th2ob(6), tconb(3), tranb
!---- temporary arrays -----
  REAL*8 :: trnth2o
  REAL*8 :: trnth2ob
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  REAL*8 :: tempb
!-----tco2 are the six exp factors between levels k1 and k2 
!     tran is the updated total transmittance between levels k1 and k2
!-----th2o is the 6 exp factors between levels k1 and k2 due to
!     h2o line absorption. 
!-----tcon is the 3 exp factors between levels k1 and k2 due to
!     h2o continuum absorption.
!-----trnth2o is the total transmittance between levels k1 and k2 due
!     to both line and continuum absorption.
!-----Compute th2o following Eq. (8.23).
  CALL PUSHREAL8(th2o(1))
  th2o(1) = th2o(1)*h2oexp(k, 1)
  CALL PUSHREAL8(th2o(2))
  th2o(2) = th2o(2)*h2oexp(k, 2)
  CALL PUSHREAL8(th2o(3))
  th2o(3) = th2o(3)*h2oexp(k, 3)
  CALL PUSHREAL8(th2o(4))
  th2o(4) = th2o(4)*h2oexp(k, 4)
  CALL PUSHREAL8(th2o(5))
  th2o(5) = th2o(5)*h2oexp(k, 5)
  CALL PUSHREAL8(th2o(6))
  th2o(6) = th2o(6)*h2oexp(k, 6)
  IF (ne .EQ. 0) THEN
    trnth2ob = tran*tranb
    th2ob(1) = th2ob(1) + fkw(1, ib)*trnth2ob
    th2ob(2) = th2ob(2) + fkw(2, ib)*trnth2ob
    th2ob(3) = th2ob(3) + fkw(3, ib)*trnth2ob
    th2ob(4) = th2ob(4) + fkw(4, ib)*trnth2ob
    th2ob(5) = th2ob(5) + fkw(5, ib)*trnth2ob
    th2ob(6) = th2ob(6) + fkw(6, ib)*trnth2ob
  ELSE IF (ne .EQ. 1) THEN
!-----Compute trnh2o following Eqs. (8.25) and (4.27).
    CALL PUSHREAL8(tcon(1))
    tcon(1) = tcon(1)*conexp(k, 1)
    trnth2ob = tran*tranb
    tempb = tcon(1)*trnth2ob
    th2ob(1) = th2ob(1) + fkw(1, ib)*tempb
    th2ob(2) = th2ob(2) + fkw(2, ib)*tempb
    th2ob(3) = th2ob(3) + fkw(3, ib)*tempb
    th2ob(4) = th2ob(4) + fkw(4, ib)*tempb
    th2ob(5) = th2ob(5) + fkw(5, ib)*tempb
    th2ob(6) = th2ob(6) + fkw(6, ib)*tempb
    tconb(1) = tconb(1) + (fkw(1, ib)*th2o(1)+fkw(2, ib)*th2o(2)+fkw(3, &
&     ib)*th2o(3)+fkw(4, ib)*th2o(4)+fkw(5, ib)*th2o(5)+fkw(6, ib)*th2o(&
&     6))*trnth2ob
    CALL POPREAL8(tcon(1))
    conexpb(k, 1) = conexpb(k, 1) + tcon(1)*tconb(1)
    tconb(1) = conexp(k, 1)*tconb(1)
  ELSE
!-----For band 3. This band is divided into 3 subbands.
    CALL PUSHREAL8(tcon(1))
    tcon(1) = tcon(1)*conexp(k, 1)
    CALL PUSHREAL8(tcon(2))
    tcon(2) = tcon(2)*conexp(k, 2)
    CALL PUSHREAL8(tcon(3))
    tcon(3) = tcon(3)*conexp(k, 3)
!-----Compute trnh2o following Eqs. (4.29) and (8.25).
    trnth2ob = tran*tranb
    tempb0 = tcon(1)*trnth2ob
    tempb1 = tcon(2)*trnth2ob
    tempb2 = tcon(3)*trnth2ob
    th2ob(1) = th2ob(1) + gkw(1, 3)*tempb2 + gkw(1, 2)*tempb1 + gkw(1, 1&
&     )*tempb0
    th2ob(2) = th2ob(2) + gkw(2, 3)*tempb2 + gkw(2, 2)*tempb1 + gkw(2, 1&
&     )*tempb0
    th2ob(3) = th2ob(3) + gkw(3, 3)*tempb2 + gkw(3, 2)*tempb1 + gkw(3, 1&
&     )*tempb0
    th2ob(4) = th2ob(4) + gkw(4, 3)*tempb2 + gkw(4, 2)*tempb1 + gkw(4, 1&
&     )*tempb0
    th2ob(5) = th2ob(5) + gkw(5, 3)*tempb2 + gkw(5, 2)*tempb1 + gkw(5, 1&
&     )*tempb0
    th2ob(6) = th2ob(6) + gkw(6, 3)*tempb2 + gkw(6, 2)*tempb1 + gkw(6, 1&
&     )*tempb0
    tconb(1) = tconb(1) + (gkw(1, 1)*th2o(1)+gkw(2, 1)*th2o(2)+gkw(3, 1)&
&     *th2o(3)+gkw(4, 1)*th2o(4)+gkw(5, 1)*th2o(5)+gkw(6, 1)*th2o(6))*&
&     trnth2ob
    tconb(2) = tconb(2) + (gkw(1, 2)*th2o(1)+gkw(2, 2)*th2o(2)+gkw(3, 2)&
&     *th2o(3)+gkw(4, 2)*th2o(4)+gkw(5, 2)*th2o(5)+gkw(6, 2)*th2o(6))*&
&     trnth2ob
    tconb(3) = tconb(3) + (gkw(1, 3)*th2o(1)+gkw(2, 3)*th2o(2)+gkw(3, 3)&
&     *th2o(3)+gkw(4, 3)*th2o(4)+gkw(5, 3)*th2o(5)+gkw(6, 3)*th2o(6))*&
&     trnth2ob
    CALL POPREAL8(tcon(3))
    conexpb(k, 3) = conexpb(k, 3) + tcon(3)*tconb(3)
    tconb(3) = conexp(k, 3)*tconb(3)
    CALL POPREAL8(tcon(2))
    conexpb(k, 2) = conexpb(k, 2) + tcon(2)*tconb(2)
    tconb(2) = conexp(k, 2)*tconb(2)
    CALL POPREAL8(tcon(1))
    conexpb(k, 1) = conexpb(k, 1) + tcon(1)*tconb(1)
    tconb(1) = conexp(k, 1)*tconb(1)
  END IF
  CALL POPREAL8(th2o(6))
  h2oexpb(k, 6) = h2oexpb(k, 6) + th2o(6)*th2ob(6)
  th2ob(6) = h2oexp(k, 6)*th2ob(6)
  CALL POPREAL8(th2o(5))
  h2oexpb(k, 5) = h2oexpb(k, 5) + th2o(5)*th2ob(5)
  th2ob(5) = h2oexp(k, 5)*th2ob(5)
  CALL POPREAL8(th2o(4))
  h2oexpb(k, 4) = h2oexpb(k, 4) + th2o(4)*th2ob(4)
  th2ob(4) = h2oexp(k, 4)*th2ob(4)
  CALL POPREAL8(th2o(3))
  h2oexpb(k, 3) = h2oexpb(k, 3) + th2o(3)*th2ob(3)
  th2ob(3) = h2oexp(k, 3)*th2ob(3)
  CALL POPREAL8(th2o(2))
  h2oexpb(k, 2) = h2oexpb(k, 2) + th2o(2)*th2ob(2)
  th2ob(2) = h2oexp(k, 2)*th2ob(2)
  CALL POPREAL8(th2o(1))
  h2oexpb(k, 1) = h2oexpb(k, 1) + th2o(1)*th2ob(1)
  th2ob(1) = h2oexp(k, 1)*th2ob(1)
END SUBROUTINE H2OKDIS_B

!  Differentiation of n2okdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tn2o n2oexp
!   with respect to varying inputs: tran tn2o n2oexp
!**********************************************************************
SUBROUTINE N2OKDIS_B(ib, np, k, n2oexp, n2oexpb, tn2o, tn2ob, tran, &
& tranb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: n2oexp(0:np, 4)
  REAL*8 :: n2oexpb(0:np, 4)
!---- updated parameters -----
  REAL*8 :: tn2o(4), tran
  REAL*8 :: tn2ob(4), tranb
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcb
  INTEGER :: branch
!-----tn2o is computed from Eq. (8.23). 
!     xc is the total n2o transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.
!-----band 6
  IF (ib .EQ. 6) THEN
    CALL PUSHREAL8(tn2o(1))
    tn2o(1) = tn2o(1)*n2oexp(k, 1)
    xc = 0.940414*tn2o(1)
    CALL PUSHREAL8(tn2o(2))
    tn2o(2) = tn2o(2)*n2oexp(k, 2)
    xc = xc + 0.059586*tn2o(2)
!-----band 7
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL8(tn2o(1))
    tn2o(1) = tn2o(1)*n2oexp(k, 1)
    xc = 0.561961*tn2o(1)
    CALL PUSHREAL8(tn2o(2))
    tn2o(2) = tn2o(2)*n2oexp(k, 2)
    xc = xc + 0.138707*tn2o(2)
    CALL PUSHREAL8(tn2o(3))
    tn2o(3) = tn2o(3)*n2oexp(k, 3)
    xc = xc + 0.240670*tn2o(3)
    CALL PUSHREAL8(tn2o(4))
    tn2o(4) = tn2o(4)*n2oexp(k, 4)
    xc = xc + 0.058662*tn2o(4)
    CALL PUSHCONTROL1B(1)
  END IF
  xcb = tran*tranb
  tranb = xc*tranb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    tn2ob(2) = tn2ob(2) + 0.059586*xcb
    CALL POPREAL8(tn2o(2))
    n2oexpb(k, 2) = n2oexpb(k, 2) + tn2o(2)*tn2ob(2)
    tn2ob(2) = n2oexp(k, 2)*tn2ob(2)
    tn2ob(1) = tn2ob(1) + 0.940414*xcb
    CALL POPREAL8(tn2o(1))
    n2oexpb(k, 1) = n2oexpb(k, 1) + tn2o(1)*tn2ob(1)
    tn2ob(1) = n2oexp(k, 1)*tn2ob(1)
  ELSE
    tn2ob(4) = tn2ob(4) + 0.058662*xcb
    CALL POPREAL8(tn2o(4))
    n2oexpb(k, 4) = n2oexpb(k, 4) + tn2o(4)*tn2ob(4)
    tn2ob(4) = n2oexp(k, 4)*tn2ob(4)
    tn2ob(3) = tn2ob(3) + 0.240670*xcb
    CALL POPREAL8(tn2o(3))
    n2oexpb(k, 3) = n2oexpb(k, 3) + tn2o(3)*tn2ob(3)
    tn2ob(3) = n2oexp(k, 3)*tn2ob(3)
    tn2ob(2) = tn2ob(2) + 0.138707*xcb
    CALL POPREAL8(tn2o(2))
    n2oexpb(k, 2) = n2oexpb(k, 2) + tn2o(2)*tn2ob(2)
    tn2ob(2) = n2oexp(k, 2)*tn2ob(2)
    tn2ob(1) = tn2ob(1) + 0.561961*xcb
    CALL POPREAL8(tn2o(1))
    n2oexpb(k, 1) = n2oexpb(k, 1) + tn2o(1)*tn2ob(1)
    tn2ob(1) = n2oexp(k, 1)*tn2ob(1)
  END IF
END SUBROUTINE N2OKDIS_B

!  Differentiation of ch4kdis in reverse (adjoint) mode:
!   gradient     of useful results: tran ch4exp tch4
!   with respect to varying inputs: tran ch4exp tch4
!**********************************************************************
SUBROUTINE CH4KDIS_B(ib, np, k, ch4exp, ch4expb, tch4, tch4b, tran, &
& tranb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: ch4exp(0:np, 4)
  REAL*8 :: ch4expb(0:np, 4)
!---- updated parameters -----
  REAL*8 :: tch4(4), tran
  REAL*8 :: tch4b(4), tranb
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcb
  INTEGER :: branch
!-----tch4 is computed from Eq. (8.23). 
!     xc is the total ch4 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.
!-----band 6
  IF (ib .EQ. 6) THEN
    CALL PUSHREAL8(tch4(1))
    tch4(1) = tch4(1)*ch4exp(k, 1)
    xc = tch4(1)
!-----band 7
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL8(tch4(1))
    tch4(1) = tch4(1)*ch4exp(k, 1)
    xc = 0.610650*tch4(1)
    CALL PUSHREAL8(tch4(2))
    tch4(2) = tch4(2)*ch4exp(k, 2)
    xc = xc + 0.280212*tch4(2)
    CALL PUSHREAL8(tch4(3))
    tch4(3) = tch4(3)*ch4exp(k, 3)
    xc = xc + 0.107349*tch4(3)
    CALL PUSHREAL8(tch4(4))
    tch4(4) = tch4(4)*ch4exp(k, 4)
    xc = xc + 0.001789*tch4(4)
    CALL PUSHCONTROL1B(1)
  END IF
  xcb = tran*tranb
  tranb = xc*tranb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    tch4b(1) = tch4b(1) + xcb
    CALL POPREAL8(tch4(1))
    ch4expb(k, 1) = ch4expb(k, 1) + tch4(1)*tch4b(1)
    tch4b(1) = ch4exp(k, 1)*tch4b(1)
  ELSE
    tch4b(4) = tch4b(4) + 0.001789*xcb
    CALL POPREAL8(tch4(4))
    ch4expb(k, 4) = ch4expb(k, 4) + tch4(4)*tch4b(4)
    tch4b(4) = ch4exp(k, 4)*tch4b(4)
    tch4b(3) = tch4b(3) + 0.107349*xcb
    CALL POPREAL8(tch4(3))
    ch4expb(k, 3) = ch4expb(k, 3) + tch4(3)*tch4b(3)
    tch4b(3) = ch4exp(k, 3)*tch4b(3)
    tch4b(2) = tch4b(2) + 0.280212*xcb
    CALL POPREAL8(tch4(2))
    ch4expb(k, 2) = ch4expb(k, 2) + tch4(2)*tch4b(2)
    tch4b(2) = ch4exp(k, 2)*tch4b(2)
    tch4b(1) = tch4b(1) + 0.610650*xcb
    CALL POPREAL8(tch4(1))
    ch4expb(k, 1) = ch4expb(k, 1) + tch4(1)*tch4b(1)
    tch4b(1) = ch4exp(k, 1)*tch4b(1)
  END IF
END SUBROUTINE CH4KDIS_B

!  Differentiation of comkdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tcom comexp
!   with respect to varying inputs: tran tcom comexp
!**********************************************************************
SUBROUTINE COMKDIS_B(ib, np, k, comexp, comexpb, tcom, tcomb, tran, &
& tranb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL*8 :: comexp(0:np, 6)
  REAL*8 :: comexpb(0:np, 6)
!---- updated parameters -----
  REAL*8 :: tcom(6), tran
  REAL*8 :: tcomb(6), tranb
!---- temporary arrays -----
  REAL*8 :: xc
  REAL*8 :: xcb
  INTEGER :: branch
!-----tcom is computed from Eq. (8.23). 
!     xc is the total co2 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 6.
!-----band 4
  IF (ib .EQ. 4) THEN
    CALL PUSHREAL8(tcom(1))
    tcom(1) = tcom(1)*comexp(k, 1)
    xc = 0.12159*tcom(1)
    CALL PUSHREAL8(tcom(2))
    tcom(2) = tcom(2)*comexp(k, 2)
    xc = xc + 0.24359*tcom(2)
    CALL PUSHREAL8(tcom(3))
    tcom(3) = tcom(3)*comexp(k, 3)
    xc = xc + 0.24981*tcom(3)
    CALL PUSHREAL8(tcom(4))
    tcom(4) = tcom(4)*comexp(k, 4)
    xc = xc + 0.26427*tcom(4)
    CALL PUSHREAL8(tcom(5))
    tcom(5) = tcom(5)*comexp(k, 5)
    xc = xc + 0.07807*tcom(5)
    CALL PUSHREAL8(tcom(6))
    tcom(6) = tcom(6)*comexp(k, 6)
    xc = xc + 0.04267*tcom(6)
!-----band 5
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL8(tcom(1))
    tcom(1) = tcom(1)*comexp(k, 1)
    xc = 0.06869*tcom(1)
    CALL PUSHREAL8(tcom(2))
    tcom(2) = tcom(2)*comexp(k, 2)
    xc = xc + 0.14795*tcom(2)
    CALL PUSHREAL8(tcom(3))
    tcom(3) = tcom(3)*comexp(k, 3)
    xc = xc + 0.19512*tcom(3)
    CALL PUSHREAL8(tcom(4))
    tcom(4) = tcom(4)*comexp(k, 4)
    xc = xc + 0.33446*tcom(4)
    CALL PUSHREAL8(tcom(5))
    tcom(5) = tcom(5)*comexp(k, 5)
    xc = xc + 0.17199*tcom(5)
    CALL PUSHREAL8(tcom(6))
    tcom(6) = tcom(6)*comexp(k, 6)
    xc = xc + 0.08179*tcom(6)
    CALL PUSHCONTROL1B(1)
  END IF
  xcb = tran*tranb
  tranb = xc*tranb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    tcomb(6) = tcomb(6) + 0.04267*xcb
    CALL POPREAL8(tcom(6))
    comexpb(k, 6) = comexpb(k, 6) + tcom(6)*tcomb(6)
    tcomb(6) = comexp(k, 6)*tcomb(6)
    tcomb(5) = tcomb(5) + 0.07807*xcb
    CALL POPREAL8(tcom(5))
    comexpb(k, 5) = comexpb(k, 5) + tcom(5)*tcomb(5)
    tcomb(5) = comexp(k, 5)*tcomb(5)
    tcomb(4) = tcomb(4) + 0.26427*xcb
    CALL POPREAL8(tcom(4))
    comexpb(k, 4) = comexpb(k, 4) + tcom(4)*tcomb(4)
    tcomb(4) = comexp(k, 4)*tcomb(4)
    tcomb(3) = tcomb(3) + 0.24981*xcb
    CALL POPREAL8(tcom(3))
    comexpb(k, 3) = comexpb(k, 3) + tcom(3)*tcomb(3)
    tcomb(3) = comexp(k, 3)*tcomb(3)
    tcomb(2) = tcomb(2) + 0.24359*xcb
    CALL POPREAL8(tcom(2))
    comexpb(k, 2) = comexpb(k, 2) + tcom(2)*tcomb(2)
    tcomb(2) = comexp(k, 2)*tcomb(2)
    tcomb(1) = tcomb(1) + 0.12159*xcb
    CALL POPREAL8(tcom(1))
    comexpb(k, 1) = comexpb(k, 1) + tcom(1)*tcomb(1)
    tcomb(1) = comexp(k, 1)*tcomb(1)
  ELSE
    tcomb(6) = tcomb(6) + 0.08179*xcb
    CALL POPREAL8(tcom(6))
    comexpb(k, 6) = comexpb(k, 6) + tcom(6)*tcomb(6)
    tcomb(6) = comexp(k, 6)*tcomb(6)
    tcomb(5) = tcomb(5) + 0.17199*xcb
    CALL POPREAL8(tcom(5))
    comexpb(k, 5) = comexpb(k, 5) + tcom(5)*tcomb(5)
    tcomb(5) = comexp(k, 5)*tcomb(5)
    tcomb(4) = tcomb(4) + 0.33446*xcb
    CALL POPREAL8(tcom(4))
    comexpb(k, 4) = comexpb(k, 4) + tcom(4)*tcomb(4)
    tcomb(4) = comexp(k, 4)*tcomb(4)
    tcomb(3) = tcomb(3) + 0.19512*xcb
    CALL POPREAL8(tcom(3))
    comexpb(k, 3) = comexpb(k, 3) + tcom(3)*tcomb(3)
    tcomb(3) = comexp(k, 3)*tcomb(3)
    tcomb(2) = tcomb(2) + 0.14795*xcb
    CALL POPREAL8(tcom(2))
    comexpb(k, 2) = comexpb(k, 2) + tcom(2)*tcomb(2)
    tcomb(2) = comexp(k, 2)*tcomb(2)
    tcomb(1) = tcomb(1) + 0.06869*xcb
    CALL POPREAL8(tcom(1))
    comexpb(k, 1) = comexpb(k, 1) + tcom(1)*tcomb(1)
    tcomb(1) = comexp(k, 1)*tcomb(1)
  END IF
END SUBROUTINE COMKDIS_B

!  Differentiation of cfckdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tcfc cfcexp
!   with respect to varying inputs: tran tcfc cfcexp
!**********************************************************************
SUBROUTINE CFCKDIS_B(np, k, cfcexp, cfcexpb, tcfc, tcfcb, tran, tranb)
  IMPLICIT NONE
!---- input parameters -----
  INTEGER :: k, np
  REAL*8 :: cfcexp(0:np)
  REAL*8 :: cfcexpb(0:np)
!---- updated parameters -----
  REAL*8 :: tcfc, tran
  REAL*8 :: tcfcb, tranb
!-----tcfc is the exp factors between levels k1 and k2. 
  CALL PUSHREAL8(tcfc)
  tcfc = tcfc*cfcexp(k)
  tcfcb = tcfcb + tran*tranb
  tranb = tcfc*tranb
  CALL POPREAL8(tcfc)
  cfcexpb(k) = cfcexpb(k) + tcfc*tcfcb
  tcfcb = cfcexp(k)*tcfcb
END SUBROUTINE CFCKDIS_B

!  Differentiation of b10kdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tn2o co2exp tcon h2oexp
!                n2oexp th2o tco2 conexp
!   with respect to varying inputs: tn2o co2exp tcon h2oexp n2oexp
!                th2o tco2 conexp
!**********************************************************************
SUBROUTINE B10KDIS_B(np, k, h2oexp, h2oexpb, conexp, conexpb, co2exp, &
& co2expb, n2oexp, n2oexpb, th2o, th2ob, tcon, tconb, tco2, tco2b, tn2o&
& , tn2ob, tran, tranb)
  IMPLICIT NONE
  INTEGER :: np, k
!---- input parameters -----
  REAL*8 :: h2oexp(0:np, 5), conexp(0:np), co2exp(0:np, 6), n2oexp(0:np&
& , 2)
  REAL*8 :: h2oexpb(0:np, 5), conexpb(0:np), co2expb(0:np, 6), n2oexpb(0&
& :np, 2)
!---- updated parameters -----
  REAL*8 :: th2o(6), tcon(3), tco2(6), tn2o(4), tran
  REAL*8 :: th2ob(6), tconb(3), tco2b(6), tn2ob(4), tranb
!---- temporary arrays -----
  REAL*8 :: xx
  REAL*8 :: xxb
!-----For h2o line. The k-distribution functions are given in Table 4.
  CALL PUSHREAL8(th2o(1))
  th2o(1) = th2o(1)*h2oexp(k, 1)
  xx = 0.3153*th2o(1)
  CALL PUSHREAL8(th2o(2))
  th2o(2) = th2o(2)*h2oexp(k, 2)
  xx = xx + 0.4604*th2o(2)
  CALL PUSHREAL8(th2o(3))
  th2o(3) = th2o(3)*h2oexp(k, 3)
  xx = xx + 0.1326*th2o(3)
  CALL PUSHREAL8(th2o(4))
  th2o(4) = th2o(4)*h2oexp(k, 4)
  xx = xx + 0.0798*th2o(4)
  CALL PUSHREAL8(th2o(5))
  th2o(5) = th2o(5)*h2oexp(k, 5)
  xx = xx + 0.0119*th2o(5)
  tran = xx
!-----For h2o continuum. Note that conexp(k,3) is for subband 3a.
  CALL PUSHREAL8(tcon(1))
  tcon(1) = tcon(1)*conexp(k)
  CALL PUSHREAL8(tran)
  tran = tran*tcon(1)
!-----For co2 (Table 6)
  CALL PUSHREAL8(tco2(1))
  tco2(1) = tco2(1)*co2exp(k, 1)
  xx = 0.2673*tco2(1)
  CALL PUSHREAL8(tco2(2))
  tco2(2) = tco2(2)*co2exp(k, 2)
  xx = xx + 0.2201*tco2(2)
  CALL PUSHREAL8(tco2(3))
  tco2(3) = tco2(3)*co2exp(k, 3)
  xx = xx + 0.2106*tco2(3)
  CALL PUSHREAL8(tco2(4))
  tco2(4) = tco2(4)*co2exp(k, 4)
  xx = xx + 0.2409*tco2(4)
  CALL PUSHREAL8(tco2(5))
  tco2(5) = tco2(5)*co2exp(k, 5)
  xx = xx + 0.0196*tco2(5)
  CALL PUSHREAL8(tco2(6))
  tco2(6) = tco2(6)*co2exp(k, 6)
  xx = xx + 0.0415*tco2(6)
  CALL PUSHREAL8(tran)
  tran = tran*xx
!-----For n2o (Table 5)
  CALL PUSHREAL8(tn2o(1))
  tn2o(1) = tn2o(1)*n2oexp(k, 1)
  CALL PUSHREAL8(xx)
  xx = 0.970831*tn2o(1)
  CALL PUSHREAL8(tn2o(2))
  tn2o(2) = tn2o(2)*n2oexp(k, 2)
  xx = xx + 0.029169*tn2o(2)
  xxb = tran*tranb
  tranb = (xx-1.0)*tranb
  tn2ob(2) = tn2ob(2) + 0.029169*xxb
  CALL POPREAL8(tn2o(2))
  n2oexpb(k, 2) = n2oexpb(k, 2) + tn2o(2)*tn2ob(2)
  tn2ob(2) = n2oexp(k, 2)*tn2ob(2)
  CALL POPREAL8(xx)
  tn2ob(1) = tn2ob(1) + 0.970831*xxb
  CALL POPREAL8(tn2o(1))
  n2oexpb(k, 1) = n2oexpb(k, 1) + tn2o(1)*tn2ob(1)
  tn2ob(1) = n2oexp(k, 1)*tn2ob(1)
  CALL POPREAL8(tran)
  xxb = tran*tranb
  tranb = xx*tranb
  tco2b(6) = tco2b(6) + 0.0415*xxb
  CALL POPREAL8(tco2(6))
  co2expb(k, 6) = co2expb(k, 6) + tco2(6)*tco2b(6)
  tco2b(6) = co2exp(k, 6)*tco2b(6)
  tco2b(5) = tco2b(5) + 0.0196*xxb
  CALL POPREAL8(tco2(5))
  co2expb(k, 5) = co2expb(k, 5) + tco2(5)*tco2b(5)
  tco2b(5) = co2exp(k, 5)*tco2b(5)
  tco2b(4) = tco2b(4) + 0.2409*xxb
  CALL POPREAL8(tco2(4))
  co2expb(k, 4) = co2expb(k, 4) + tco2(4)*tco2b(4)
  tco2b(4) = co2exp(k, 4)*tco2b(4)
  tco2b(3) = tco2b(3) + 0.2106*xxb
  CALL POPREAL8(tco2(3))
  co2expb(k, 3) = co2expb(k, 3) + tco2(3)*tco2b(3)
  tco2b(3) = co2exp(k, 3)*tco2b(3)
  tco2b(2) = tco2b(2) + 0.2201*xxb
  CALL POPREAL8(tco2(2))
  co2expb(k, 2) = co2expb(k, 2) + tco2(2)*tco2b(2)
  tco2b(2) = co2exp(k, 2)*tco2b(2)
  tco2b(1) = tco2b(1) + 0.2673*xxb
  CALL POPREAL8(tco2(1))
  co2expb(k, 1) = co2expb(k, 1) + tco2(1)*tco2b(1)
  tco2b(1) = co2exp(k, 1)*tco2b(1)
  CALL POPREAL8(tran)
  tconb(1) = tconb(1) + tran*tranb
  tranb = tcon(1)*tranb
  CALL POPREAL8(tcon(1))
  conexpb(k) = conexpb(k) + tcon(1)*tconb(1)
  tconb(1) = conexp(k)*tconb(1)
  xxb = tranb
  th2ob(5) = th2ob(5) + 0.0119*xxb
  CALL POPREAL8(th2o(5))
  h2oexpb(k, 5) = h2oexpb(k, 5) + th2o(5)*th2ob(5)
  th2ob(5) = h2oexp(k, 5)*th2ob(5)
  th2ob(4) = th2ob(4) + 0.0798*xxb
  CALL POPREAL8(th2o(4))
  h2oexpb(k, 4) = h2oexpb(k, 4) + th2o(4)*th2ob(4)
  th2ob(4) = h2oexp(k, 4)*th2ob(4)
  th2ob(3) = th2ob(3) + 0.1326*xxb
  CALL POPREAL8(th2o(3))
  h2oexpb(k, 3) = h2oexpb(k, 3) + th2o(3)*th2ob(3)
  th2ob(3) = h2oexp(k, 3)*th2ob(3)
  th2ob(2) = th2ob(2) + 0.4604*xxb
  CALL POPREAL8(th2o(2))
  h2oexpb(k, 2) = h2oexpb(k, 2) + th2o(2)*th2ob(2)
  th2ob(2) = h2oexp(k, 2)*th2ob(2)
  th2ob(1) = th2ob(1) + 0.3153*xxb
  CALL POPREAL8(th2o(1))
  h2oexpb(k, 1) = h2oexpb(k, 1) + th2o(1)*th2ob(1)
  th2ob(1) = h2oexp(k, 1)*th2ob(1)
END SUBROUTINE B10KDIS_B

!  Differentiation of cldovlp in reverse (adjoint) mode:
!   gradient     of useful results: cldhi ett cldlw enn cldmd
!   with respect to varying inputs: cldhi ett cldlw enn cldmd
!mjs
!***********************************************************************
SUBROUTINE CLDOVLP_B(np, k1, k2, ict, icb, icx, ncld, enn, ennb, ett, &
& ettb, cldhi, cldhib, cldmd, cldmdb, cldlw, cldlwb)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: np, k1, k2, ict, icb, icx(0:np), ncld(3)
  REAL*8, INTENT(IN) :: enn(0:np), ett(0:np)
  REAL*8 :: ennb(0:np), ettb(0:np)
  REAL*8, INTENT(INOUT) :: cldhi, cldmd, cldlw
  REAL*8 :: cldhib, cldmdb, cldlwb
  INTEGER :: j, k, km, kx
  INTEGER :: branch
  km = k2 - 1
  IF (km .LT. ict) THEN
! do high clouds
    kx = ncld(1)
    IF (kx .EQ. 1 .OR. cldhi .EQ. 0.) THEN
      ennb(km) = ennb(km) + cldhib
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      CALL PUSHREAL8(cldhi)
      cldhi = 0.0
      IF (kx .NE. 0) THEN
        DO k=ict-kx,ict-1
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            CALL PUSHREAL8(cldhi)
            cldhi = enn(j) + ett(j)*cldhi
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO k=ict-1,ict-kx,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            j = icx(k)
            CALL POPREAL8(cldhi)
            ennb(j) = ennb(j) + cldhib
            ettb(j) = ettb(j) + cldhi*cldhib
            cldhib = ett(j)*cldhib
          END IF
        END DO
      END IF
      CALL POPREAL8(cldhi)
    END IF
    cldhib = 0.0_8
  ELSE IF (km .GE. ict .AND. km .LT. icb) THEN
! do middle clouds
    kx = ncld(2)
    IF (kx .EQ. 1 .OR. cldmd .EQ. 0.) THEN
      ennb(km) = ennb(km) + cldmdb
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      CALL PUSHREAL8(cldmd)
      cldmd = 0.0
      IF (kx .NE. 0) THEN
        DO k=icb-kx,icb-1
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            CALL PUSHREAL8(cldmd)
            cldmd = enn(j) + ett(j)*cldmd
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO k=icb-1,icb-kx,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            j = icx(k)
            CALL POPREAL8(cldmd)
            ennb(j) = ennb(j) + cldmdb
            ettb(j) = ettb(j) + cldmd*cldmdb
            cldmdb = ett(j)*cldmdb
          END IF
        END DO
      END IF
      CALL POPREAL8(cldmd)
    END IF
    cldmdb = 0.0_8
  ELSE
! do low clouds
    kx = ncld(3)
    IF (kx .EQ. 1 .OR. cldlw .EQ. 0.) THEN
      ennb(km) = ennb(km) + cldlwb
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      CALL PUSHREAL8(cldlw)
      cldlw = 0.0
      IF (kx .NE. 0) THEN
        DO k=np+1-kx,np
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            CALL PUSHREAL8(cldlw)
            cldlw = enn(j) + ett(j)*cldlw
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO k=np,np+1-kx,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            j = icx(k)
            CALL POPREAL8(cldlw)
            ennb(j) = ennb(j) + cldlwb
            ettb(j) = ettb(j) + cldlw*cldlwb
            cldlwb = ett(j)*cldlwb
          END IF
        END DO
      END IF
      CALL POPREAL8(cldlw)
    END IF
    cldlwb = 0.0_8
  END IF
END SUBROUTINE CLDOVLP_B

!  Differentiation of getirtau1 in reverse (adjoint) mode:
!   gradient     of useful results: hydromets tcldlyr fcld enn
!                reff
!   with respect to varying inputs: hydromets tcldlyr fcld enn
!                reff
SUBROUTINE GETIRTAU1_B(ib, nlevs, dp, fcld, fcldb, reff, reffb, &
& hydromets, hydrometsb, tcldlyr, tcldlyrb, enn, ennb, aib_ir1, awb_ir1&
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
  REAL*8 :: fcldb(nlevs)
!  Effective radius (microns)
  REAL*8, INTENT(IN) :: reff(nlevs, 4)
  REAL*8 :: reffb(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL*8, INTENT(IN) :: hydromets(nlevs, 4)
  REAL*8 :: hydrometsb(nlevs, 4)
  REAL*8, INTENT(IN) :: aib_ir1(3, 10), awb_ir1(4, 10), aiw_ir1(4, 10)
  REAL*8, INTENT(IN) :: aww_ir1(4, 10), aig_ir1(4, 10), awg_ir1(4, 10)
  REAL*8, INTENT(IN) :: cons_grav
! !OUTPUT PARAMETERS:
!  Flux transmissivity?
  REAL*8 :: tcldlyr(0:nlevs)
  REAL*8 :: tcldlyrb(0:nlevs)
!  Flux transmissivity of a cloud layer?
  REAL*8 :: enn(0:nlevs)
  REAL*8 :: ennb(0:nlevs)
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
  REAL*8 :: taucld1b, taucld2b, taucld3b, taucld4b
  REAL*8 :: g1, g2, g3, g4, gg
  REAL*8 :: g1b, g2b, g3b, g4b, ggb
  REAL*8 :: w1, w2, w3, w4, ww
  REAL*8 :: w1b, w2b, w3b, w4b, wwb
  REAL*8 :: ff, tauc
  REAL*8 :: ffb, taucb
  REAL*8 :: reff_snow
  REAL*8 :: reff_snowb
  INTRINSIC MIN
  INTRINSIC ABS
  INTRINSIC MAX
  INTRINSIC EXP
!-----Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!     Rain optical thickness is set to 0.00307 /(gm/m**2).
!     It is for a specific drop size distribution provided by Q. Fu.
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
  REAL*8 :: tempb12
  REAL*8 :: tempb11
  REAL*8 :: tempb10
  REAL*8 :: temp16
  REAL*8 :: temp15
  REAL*8 :: temp14
  REAL*8 :: temp13
  REAL*8 :: temp12
  REAL*8 :: temp11
  REAL*8 :: temp10
  REAL*8 :: max1b
  REAL*8 :: tempb
  REAL*8 :: abs0
  REAL*8 :: temp
  REAL*8 :: max1
  REAL*8 :: temp9
  REAL*8 :: temp8
  REAL*8 :: temp7
  REAL*8 :: temp6
  REAL*8 :: temp5
  REAL*8 :: temp4
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.0) THEN
      CALL PUSHREAL8(taucld1)
      taucld1 = 0.0
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL8(taucld1)
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*(aib_ir1(1, ib)+&
&       aib_ir1(2, ib)/reff(k, 1)**aib_ir1(3, ib))
      CALL PUSHCONTROL1B(0)
    END IF
    CALL PUSHREAL8(taucld2)
    taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_ir1(1, ib)+(&
&     awb_ir1(2, ib)+(awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reff(k, &
&     2))*reff(k, 2))
    taucld3 = 0.00307*(dp(k)*1.0e3/cons_grav*hydromets(k, 3))
    IF (reff(k, 4) .GT. 112.0) THEN
      CALL PUSHREAL8(reff_snow)
      reff_snow = 112.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL8(reff_snow)
      reff_snow = reff(k, 4)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff_snow .LE. 0.0) THEN
      CALL PUSHREAL8(taucld4)
      taucld4 = 0.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL8(taucld4)
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*(aib_ir1(1, ib)+&
&       aib_ir1(2, ib)/reff_snow**aib_ir1(3, ib))
      CALL PUSHCONTROL1B(1)
    END IF
!-----Compute cloud single-scattering albedo and asymmetry factor for
!     a mixture of ice particles and liquid drops following 
!     Eqs. (6.5), (6.6), (6.15) and (6.16).
!     Single-scattering albedo and asymmetry factor of rain are set
!     to 0.54 and 0.95, respectively, based on the information provided
!     by Prof. Qiang Fu.
    CALL PUSHREAL8(tauc)
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      CALL PUSHREAL8(w1)
      w1 = taucld1*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reff(k, 1))
      CALL PUSHREAL8(w2)
      w2 = taucld2*(aww_ir1(1, ib)+(aww_ir1(2, ib)+(aww_ir1(3, ib)+&
&       aww_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reff(k, 2))
      w3 = taucld3*0.54
      CALL PUSHREAL8(w4)
      w4 = taucld4*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff_snow)*reff_snow)*reff_snow)
      ww = (w1+w2+w3+w4)/tauc
      CALL PUSHREAL8(g1)
      g1 = w1*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff(k, 1))*reff(k, 1))*reff(k, 1))
      CALL PUSHREAL8(g2)
      g2 = w2*(awg_ir1(1, ib)+(awg_ir1(2, ib)+(awg_ir1(3, ib)+awg_ir1(4&
&       , ib)*reff(k, 2))*reff(k, 2))*reff(k, 2))
      g3 = w3*0.95
      CALL PUSHREAL8(g4)
      g4 = w4*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff_snow)*reff_snow)*reff_snow)
      IF (w1 + w2 + w3 + w4 .GE. 0.) THEN
        abs0 = w1 + w2 + w3 + w4
      ELSE
        abs0 = -(w1+w2+w3+w4)
      END IF
      IF (abs0 .GT. 0.0) THEN
        CALL PUSHREAL8(gg)
        gg = (g1+g2+g3+g4)/(w1+w2+w3+w4)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHREAL8(gg)
        gg = 0.5
        CALL PUSHCONTROL1B(0)
      END IF
!-----Parameterization of LW scattering following Eqs. (6.11)
!     and (6.12). 
      ff = 0.5 + (0.3739+(0.0076+0.1185*gg)*gg)*gg
      IF (1. - ww*ff .LT. 0.0) THEN
        CALL PUSHREAL8(max1)
        max1 = 0.0
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHREAL8(max1)
        max1 = 1. - ww*ff
        CALL PUSHCONTROL1B(1)
      END IF
!ALT: temporary protection against negative cloud optical thickness
      CALL PUSHREAL8(tauc)
      tauc = max1*tauc
!-----compute cloud diffuse transmittance. It is approximated by using 
!     a diffusivity factor of 1.66.
      CALL PUSHREAL8(tcldlyr(k))
      tcldlyr(k) = EXP(-(1.66*tauc))
! N in the documentation (6.13)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=nlevs,1,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      ennb(k) = 0.0_8
      tcldlyrb(k) = 0.0_8
      reff_snowb = 0.0_8
      taucld1b = 0.0_8
      taucld2b = 0.0_8
      taucld3b = 0.0_8
      taucld4b = 0.0_8
      taucb = 0.0_8
    ELSE
      fcldb(k) = fcldb(k) + (1.0-tcldlyr(k))*ennb(k)
      tcldlyrb(k) = tcldlyrb(k) - fcld(k)*ennb(k)
      ennb(k) = 0.0_8
      CALL POPREAL8(tcldlyr(k))
      taucb = -(EXP(-(1.66*tauc))*1.66*tcldlyrb(k))
      tcldlyrb(k) = 0.0_8
      CALL POPREAL8(tauc)
      max1b = tauc*taucb
      taucb = max1*taucb
      taucld3 = 0.00307*(dp(k)*1.0e3/cons_grav*hydromets(k, 3))
      w3 = taucld3*0.54
      ww = (w1+w2+w3+w4)/tauc
      ff = 0.5 + (0.3739+(0.0076+0.1185*gg)*gg)*gg
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8(max1)
        wwb = 0.0_8
        ffb = 0.0_8
      ELSE
        CALL POPREAL8(max1)
        wwb = -(ff*max1b)
        ffb = -(ww*max1b)
      END IF
      ggb = ((0.1185*gg+0.0076)*gg+gg*(0.1185*gg+0.0076)+gg**2*0.1185+&
&       0.3739)*ffb
      g3 = w3*0.95
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8(gg)
        w1b = 0.0_8
        w2b = 0.0_8
        w3b = 0.0_8
        w4b = 0.0_8
        g1b = 0.0_8
        g2b = 0.0_8
        g3b = 0.0_8
        g4b = 0.0_8
      ELSE
        CALL POPREAL8(gg)
        tempb11 = ggb/(w1+w2+w3+w4)
        tempb12 = -((g1+g2+g3+g4)*tempb11/(w1+w2+w3+w4))
        g1b = tempb11
        g2b = tempb11
        g3b = tempb11
        g4b = tempb11
        w1b = tempb12
        w2b = tempb12
        w3b = tempb12
        w4b = tempb12
      END IF
      tempb5 = wwb/tauc
      CALL POPREAL8(g4)
      temp16 = aig_ir1(3, ib) + aig_ir1(4, ib)*reff_snow
      temp15 = aig_ir1(2, ib) + temp16*reff_snow
      tempb4 = w4*g4b
      w4b = w4b + tempb5 + (aig_ir1(1, ib)+temp15*reff_snow)*g4b
      w3b = w3b + tempb5 + 0.95*g3b
      CALL POPREAL8(g2)
      temp14 = awg_ir1(3, ib) + awg_ir1(4, ib)*reff(k, 2)
      temp13 = awg_ir1(2, ib) + temp14*reff(k, 2)
      tempb7 = w2*reff(k, 2)*g2b
      w2b = w2b + tempb5 + (awg_ir1(1, ib)+temp13*reff(k, 2))*g2b
      reffb(k, 2) = reffb(k, 2) + w2*temp13*g2b + (temp14+reff(k, 2)*&
&       awg_ir1(4, ib))*tempb7
      CALL POPREAL8(g1)
      temp12 = aig_ir1(3, ib) + aig_ir1(4, ib)*reff(k, 1)
      temp11 = aig_ir1(2, ib) + temp12*reff(k, 1)
      tempb8 = w1*reff(k, 1)*g1b
      w1b = w1b + tempb5 + (aig_ir1(1, ib)+temp11*reff(k, 1))*g1b
      reffb(k, 1) = reffb(k, 1) + w1*temp11*g1b + (temp12+reff(k, 1)*&
&       aig_ir1(4, ib))*tempb8
      taucb = taucb - (w1+w2+w3+w4)*tempb5/tauc
      CALL POPREAL8(w4)
      temp10 = aiw_ir1(3, ib) + aiw_ir1(4, ib)*reff_snow
      temp9 = aiw_ir1(2, ib) + temp10*reff_snow
      tempb6 = taucld4*w4b
      reff_snowb = (temp9+reff_snow*temp10+reff_snow**2*aiw_ir1(4, ib))*&
&       tempb6 + (temp15+reff_snow*temp16+reff_snow**2*aig_ir1(4, ib))*&
&       tempb4
      taucld4b = (aiw_ir1(1, ib)+temp9*reff_snow)*w4b
      taucld3b = 0.54*w3b
      CALL POPREAL8(w2)
      temp8 = aww_ir1(3, ib) + aww_ir1(4, ib)*reff(k, 2)
      temp7 = aww_ir1(2, ib) + temp8*reff(k, 2)
      tempb9 = taucld2*reff(k, 2)*w2b
      taucld2b = (aww_ir1(1, ib)+temp7*reff(k, 2))*w2b
      reffb(k, 2) = reffb(k, 2) + taucld2*temp7*w2b + (temp8+reff(k, 2)*&
&       aww_ir1(4, ib))*tempb9
      CALL POPREAL8(w1)
      temp6 = aiw_ir1(3, ib) + aiw_ir1(4, ib)*reff(k, 1)
      temp5 = aiw_ir1(2, ib) + temp6*reff(k, 1)
      tempb10 = taucld1*reff(k, 1)*w1b
      taucld1b = (aiw_ir1(1, ib)+temp5*reff(k, 1))*w1b
      reffb(k, 1) = reffb(k, 1) + taucld1*temp5*w1b + (temp6+reff(k, 1)*&
&       aiw_ir1(4, ib))*tempb10
    END IF
    CALL POPREAL8(tauc)
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(taucld4)
    ELSE
      CALL POPREAL8(taucld4)
      temp4 = reff_snow**aib_ir1(3, ib)
      temp3 = aib_ir1(2, ib)/temp4
      tempb3 = dp(k)*1.0e3*taucld4b
      hydrometsb(k, 4) = hydrometsb(k, 4) + (aib_ir1(1, ib)+temp3)*&
&       tempb3/cons_grav
      IF (.NOT.(reff_snow .LE. 0.0 .AND. (aib_ir1(3, ib) .EQ. 0.0 .OR. &
&         aib_ir1(3, ib) .NE. INT(aib_ir1(3, ib))))) reff_snowb = &
&         reff_snowb - temp3*hydromets(k, 4)*aib_ir1(3, ib)*reff_snow**(&
&         aib_ir1(3, ib)-1)*tempb3/(temp4*cons_grav)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(reff_snow)
    ELSE
      CALL POPREAL8(reff_snow)
      reffb(k, 4) = reffb(k, 4) + reff_snowb
    END IF
    hydrometsb(k, 3) = hydrometsb(k, 3) + dp(k)*1.0e3*0.00307*taucld3b/&
&     cons_grav
    CALL POPREAL8(taucld2)
    temp2 = awb_ir1(3, ib) + awb_ir1(4, ib)*reff(k, 2)
    temp1 = awb_ir1(2, ib) + temp2*reff(k, 2)
    tempb0 = dp(k)*1.0e3*taucld2b
    tempb1 = hydromets(k, 2)*tempb0/cons_grav
    tempb2 = reff(k, 2)*tempb1
    hydrometsb(k, 2) = hydrometsb(k, 2) + (awb_ir1(1, ib)+temp1*reff(k, &
&     2))*tempb0/cons_grav
    reffb(k, 2) = reffb(k, 2) + temp1*tempb1 + (temp2+reff(k, 2)*awb_ir1&
&     (4, ib))*tempb2
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(taucld1)
      temp0 = reff(k, 1)**aib_ir1(3, ib)
      temp = aib_ir1(2, ib)/temp0
      tempb = dp(k)*1.0e3*taucld1b
      hydrometsb(k, 1) = hydrometsb(k, 1) + (aib_ir1(1, ib)+temp)*tempb/&
&       cons_grav
      IF (.NOT.(reff(k, 1) .LE. 0.0 .AND. (aib_ir1(3, ib) .EQ. 0.0 .OR. &
&         aib_ir1(3, ib) .NE. INT(aib_ir1(3, ib))))) reffb(k, 1) = reffb&
&         (k, 1) - temp*hydromets(k, 1)*aib_ir1(3, ib)*reff(k, 1)**(&
&         aib_ir1(3, ib)-1)*tempb/(temp0*cons_grav)
    ELSE
      CALL POPREAL8(taucld1)
    END IF
  END DO
  ennb(0) = 0.0_8
  tcldlyrb(0) = 0.0_8
END SUBROUTINE GETIRTAU1_B

end module IRRAD_AD
