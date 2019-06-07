module CLOUDRADCOUPLE

IMPLICIT NONE

PRIVATE
PUBLIC :: RADCOUPLE, RADCOUPLE_d, RADCOUPLE_b

contains

subroutine RADCOUPLE(  TE,              & 
                       PL,              & 
                       CF,              & 
                       AF,              & 
                       QClLS,           & 
                       QCiLS,           & 
                       QClAN,           & 
                       QCiAN,           & 
                       RAD_QL,          &  
                       RAD_QI,          & 
                       RAD_CF,          & 
                       RAD_RL,          & 
                       RAD_RI,          & 
                       TEMPOR           )

 IMPLICIT NONE

 !Inputs
 real(8), intent(in ) :: TE, PL, TEMPOR
 real(8), intent(in ) :: AF, CF, QClAN, QCiAN, QClLS, QCiLS
! real(8), intent(in ) :: QRN_ALL, QSN_ALL

 !Outputs
 real(8), intent(out) :: RAD_QL,RAD_QI,RAD_CF,RAD_RL,RAD_RI
! real(8), intent(out) :: RAD_QR,RAD_QS

 !Locals
 real(8) :: ss, RAD_RI_AN, AFx, ALPH

 real(8), parameter :: MIN_RI = 20.e-6, MAX_RI = 40.e-6, RI_ANV = 30.e-6

 !Initialize outputs
 RAD_QL = 0.0
 RAD_QI = 0.0
 RAD_CF = 0.0
 RAD_RL = 0.0
 RAD_RI = 0.0
 !RAD_QR = 0.0
 !RAD_QS = 0.0

 ! Adjust Anvil fractions for warm clouds
 ALPH =  0.1
 SS   =  (280.-TE)/20.
 SS   =  MIN( 1.0 , SS )
 SS   =  MAX( 0.0 , SS )

 SS   =  ALPH + (SS**3) * ( 1.0 - ALPH )

 AFx  =  AF * SS * 0.5

 !Total cloud fraction
 RAD_CF = MIN( CF + AFx, 1.00 )

 !Total In-cloud liquid
 if ( RAD_CF > 10.0e-8 ) then  !0 -> 10e-8 FOR LINEARIZATION PROTECTION
    RAD_QL = ( QClLS + QClAN ) / RAD_CF
 else
    RAD_QL = 0.0
 end if
 RAD_QL = MIN( RAD_QL, 0.01 )

 ! Total In-cloud ice
 if (  RAD_CF > 10.0e-8 ) then !0 -> 10e-8 FOR LINEARIZATION PROTECTION
    RAD_QI = ( QCiLS + QCiAN ) / RAD_CF
 else
    RAD_QI = 0.0
 end if
 RAD_QI = MIN( RAD_QI, 0.01 )

 ! Total In-cloud precipitation
! if (  RAD_CF >0. ) then
!    RAD_QR = ( QRN_ALL ) / RAD_CF
!    RAD_QS = ( QSN_ALL ) / RAD_CF
! else
!    RAD_QR = 0.0
!    RAD_QS = 0.0
! end if
! RAD_QR = MIN( RAD_QR, 0.01 )
! RAD_QS = MIN( RAD_QS, 0.01 )

 if (PL < 150. ) then
    RAD_RI = MAX_RI
 end if
 if (PL >= 150. ) then
    RAD_RI = MAX_RI*150./PL
 end if

 ! Weigh in a separate R_ice for Anvil Ice according to
 RAD_RI_AN  =  RAD_RI  

 if ( ( QCiLS + QCiAN ) > 0.0 ) then
    if (qcils/rad_ri+qcian/ri_anv .gt. 10e-8) then !LINEARIZATION PROTECTION
       RAD_RI_AN  = ( QCiLS + QCiAN ) / ( (QCiLS/RAD_RI) + (QCiAN/RI_ANV) )
    endif
 end if

 RAD_RI = MIN( RAD_RI, RAD_RI_AN )
 RAD_RI = MAX( RAD_RI, MIN_RI )

 ! Implement ramps for gradual change in effective radius
 if (PL < 300. ) then
    RAD_RL = 21.e-6
 end if
 if (PL >= 300. ) then
    RAD_RL = 21.e-6*300./PL
 end if
 RAD_RL = MAX( RAD_RL, 10.e-6 )

 ! Thicken low high lat clouds
 if ( PL .GE. 775.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
    RAD_RL = max(min(-0.1 * PL + 87.5, 10.),5.)*1.e-6
 end if
 if ( PL .GE. 825.  .AND. TE .LE.  282. .AND. (tempor.eq.1.) ) then
    RAD_RL = max(0.71 * TE - 190.25, 5.)*1.e-6
 end if
 if ( PL .GE. 775.  .AND. PL .LT. 825. .AND. TE .LE.  282. .AND. TE .GT. 275. .AND. (tempor.eq.1.) ) then
    RAD_RL = min(-0.1*PL + 0.71 * TE - 107.75, 10.)*1.e-6
 end if
 if ( PL .GE. 825.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
    RAD_RL = 5.*1.e-6
 end if

 ! Thin low tropical clouds
 if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
    RAD_RL = min(2.2 * TE - 617., 21.)*1.e-6
 end if
 if ( PL .GE. 925.  .AND. TE .GE.  290. ) then
    RAD_RL = min(0.44 * PL - 397., 21.)*1.e-6
 end if
 if ( PL .GE. 925.  .AND. PL .LT. 950. .AND. TE .GT.  285. .AND. TE .LT. 290.) then
    RAD_RL = max(min(0.44*PL + 2.2 * TE - 1035., 21.),10.)*1.e-6
 end if
 if ( PL .GE. 950.  .AND. TE .GE.  290. ) then
    RAD_RL = 21.*1.e-6
 end if

 if ( RAD_CF < 1.e-5 ) then
    RAD_QL = 0.
    RAD_QI = 0.
    RAD_CF = 0.
    !RAD_QR = 0.
    !RAD_QS = 0.
 end if

end subroutine RADCOUPLE

SUBROUTINE RADCOUPLE_D(te, ted, pl, cf, cfd, af, afd, qclls, qcllsd, &
& qcils, qcilsd, qclan, qcland, qcian, qciand, rad_ql, rad_qld, rad_qi, &
& rad_qid, rad_cf, rad_cfd, rad_rl, rad_rld, rad_ri, rad_rid, tempor)
  IMPLICIT NONE
!Inputs
  REAL*8, INTENT(IN) :: te, pl, tempor
  REAL*8, INTENT(IN) :: ted
  REAL*8, INTENT(IN) :: af, cf, qclan, qcian, qclls, qcils
  REAL*8, INTENT(IN) :: afd, cfd, qcland, qciand, qcllsd, qcilsd
! real(8), intent(in ) :: QRN_ALL, QSN_ALL
!Outputs
  REAL*8, INTENT(OUT) :: rad_ql, rad_qi, rad_cf, rad_rl, rad_ri
  REAL*8, INTENT(OUT) :: rad_qld, rad_qid, rad_cfd, rad_rld, rad_rid
! real(8), intent(out) :: RAD_QR,RAD_QS
!Locals
  REAL*8 :: ss, rad_ri_an, afx, alph
  REAL*8 :: ssd, rad_ri_and, afxd
  REAL*8, PARAMETER :: min_ri=20.e-6, max_ri=40.e-6, ri_anv=30.e-6
  INTRINSIC MIN
  INTRINSIC MAX
  REAL*8 :: max2d
  REAL*8 :: min3
  REAL*8 :: min2
  REAL*8 :: min1
  REAL*8 :: min1d
  REAL*8 :: x2
  REAL*8 :: x2d
  REAL*8 :: x1
  REAL*8 :: max3d
  REAL*8 :: max3
  REAL*8 :: max2
  REAL*8 :: max1
  REAL*8 :: min2d
!Initialize outputs
  rad_ql = 0.0
  rad_qi = 0.0
  rad_cf = 0.0
  rad_rl = 0.0
  rad_ri = 0.0
!RAD_QR = 0.0
!RAD_QS = 0.0
! Adjust Anvil fractions for warm clouds
  alph = 0.1
  ssd = (-ted)/20.
  ss = (280.-te)/20.
  IF (1.0 .GT. ss) THEN
    ss = ss
  ELSE
    ss = 1.0
    ssd = 0.0_8
  END IF
  IF (0.0 .LT. ss) THEN
    ss = ss
  ELSE
    ss = 0.0
    ssd = 0.0_8
  END IF
  ssd = (1.0-alph)*3*ss**2*ssd
  ss = alph + ss**3*(1.0-alph)
  afxd = 0.5*(afd*ss+af*ssd)
  afx = af*ss*0.5
  IF (cf + afx .GT. 1.00) THEN
    rad_cf = 1.00
    rad_cfd = 0.0_8
  ELSE
    rad_cfd = cfd + afxd
    rad_cf = cf + afx
  END IF
!Total In-cloud liquid
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_qld = ((qcllsd+qcland)*rad_cf-(qclls+qclan)*rad_cfd)/rad_cf**2
    rad_ql = (qclls+qclan)/rad_cf
  ELSE
    rad_ql = 0.0
    rad_qld = 0.0_8
  END IF
  IF (rad_ql .GT. 0.01) THEN
    rad_ql = 0.01
    rad_qld = 0.0_8
  ELSE
    rad_ql = rad_ql
  END IF
! Total In-cloud ice
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_qid = ((qcilsd+qciand)*rad_cf-(qcils+qcian)*rad_cfd)/rad_cf**2
    rad_qi = (qcils+qcian)/rad_cf
  ELSE
    rad_qi = 0.0
    rad_qid = 0.0_8
  END IF
  IF (rad_qi .GT. 0.01) THEN
    rad_qi = 0.01
    rad_qid = 0.0_8
  ELSE
    rad_qi = rad_qi
  END IF
! Total In-cloud precipitation
! if (  RAD_CF >0. ) then
!    RAD_QR = ( QRN_ALL ) / RAD_CF
!    RAD_QS = ( QSN_ALL ) / RAD_CF
! else
!    RAD_QR = 0.0
!    RAD_QS = 0.0
! end if
! RAD_QR = MIN( RAD_QR, 0.01 )
! RAD_QS = MIN( RAD_QS, 0.01 )
  IF (pl .LT. 150.) rad_ri = max_ri
  IF (pl .GE. 150.) rad_ri = max_ri*150./pl
! Weigh in a separate R_ice for Anvil Ice according to
  rad_ri_an = rad_ri
  IF (qcils + qcian .GT. 0.0) THEN
    IF (qcils/rad_ri + qcian/ri_anv .GT. 10e-8) THEN
!LINEARIZATION PROTECTION
      rad_ri_and = ((qcilsd+qciand)*(qcils/rad_ri+qcian/ri_anv)-(qcils+&
&       qcian)*(qcilsd/rad_ri+qciand/ri_anv))/(qcils/rad_ri+qcian/ri_anv&
&       )**2
      rad_ri_an = (qcils+qcian)/(qcils/rad_ri+qcian/ri_anv)
    ELSE
      rad_ri_and = 0.0_8
    END IF
  ELSE
    rad_ri_and = 0.0_8
  END IF
  IF (rad_ri .GT. rad_ri_an) THEN
    rad_rid = rad_ri_and
    rad_ri = rad_ri_an
  ELSE
    rad_ri = rad_ri
    rad_rid = 0.0_8
  END IF
  IF (rad_ri .LT. min_ri) THEN
    rad_ri = min_ri
    rad_rid = 0.0_8
  ELSE
    rad_ri = rad_ri
  END IF
! Implement ramps for gradual change in effective radius
  IF (pl .LT. 300.) rad_rl = 21.e-6
  IF (pl .GE. 300.) rad_rl = 21.e-6*300./pl
  IF (rad_rl .LT. 10.e-6) THEN
    rad_rl = 10.e-6
  ELSE
    rad_rl = rad_rl
  END IF
! Thicken low high lat clouds
  IF (pl .GE. 775. .AND. te .LE. 275. .AND. tempor .EQ. 1.) THEN
    IF (-(0.1*pl) + 87.5 .GT. 10.) THEN
      x1 = 10.
    ELSE
      x1 = -(0.1*pl) + 87.5
    END IF
    IF (x1 .LT. 5.) THEN
      max1 = 5.
    ELSE
      max1 = x1
    END IF
    rad_rl = max1*1.e-6
  END IF
  IF (pl .GE. 825. .AND. te .LE. 282. .AND. tempor .EQ. 1.) THEN
    IF (0.71*te - 190.25 .LT. 5.) THEN
      max2 = 5.
      max2d = 0.0_8
    ELSE
      max2d = 0.71*ted
      max2 = 0.71*te - 190.25
    END IF
    rad_rld = 1.e-6*max2d
    rad_rl = max2*1.e-6
  ELSE
    rad_rld = 0.0_8
  END IF
  IF (pl .GE. 775. .AND. pl .LT. 825. .AND. te .LE. 282. .AND. te .GT. &
&     275. .AND. tempor .EQ. 1.) THEN
    IF (-(0.1*pl) + 0.71*te - 107.75 .GT. 10.) THEN
      min1 = 10.
      min1d = 0.0_8
    ELSE
      min1d = 0.71*ted
      min1 = -(0.1*pl) + 0.71*te - 107.75
    END IF
    rad_rld = 1.e-6*min1d
    rad_rl = min1*1.e-6
  END IF
  IF (pl .GE. 825. .AND. te .LE. 275. .AND. tempor .EQ. 1.) THEN
    rad_rld = 0.0_8
    rad_rl = 5.*1.e-6
  END IF
! Thin low tropical clouds
  IF (pl .GE. 950. .AND. te .GE. 285.) THEN
    IF (2.2*te - 617. .GT. 21.) THEN
      min2 = 21.
      min2d = 0.0_8
    ELSE
      min2d = 2.2*ted
      min2 = 2.2*te - 617.
    END IF
    rad_rld = 1.e-6*min2d
    rad_rl = min2*1.e-6
  END IF
  IF (pl .GE. 925. .AND. te .GE. 290.) THEN
    IF (0.44*pl - 397. .GT. 21.) THEN
      min3 = 21.
    ELSE
      min3 = 0.44*pl - 397.
    END IF
    rad_rl = min3*1.e-6
    rad_rld = 0.0_8
  END IF
  IF (pl .GE. 925. .AND. pl .LT. 950. .AND. te .GT. 285. .AND. te .LT. &
&     290.) THEN
    IF (0.44*pl + 2.2*te - 1035. .GT. 21.) THEN
      x2 = 21.
      x2d = 0.0_8
    ELSE
      x2d = 2.2*ted
      x2 = 0.44*pl + 2.2*te - 1035.
    END IF
    IF (x2 .LT. 10.) THEN
      max3 = 10.
      max3d = 0.0_8
    ELSE
      max3d = x2d
      max3 = x2
    END IF
    rad_rld = 1.e-6*max3d
    rad_rl = max3*1.e-6
  END IF
  IF (pl .GE. 950. .AND. te .GE. 290.) THEN
    rad_rld = 0.0_8
    rad_rl = 21.*1.e-6
  END IF
  IF (rad_cf .LT. 1.e-5) THEN
    rad_ql = 0.
    rad_qi = 0.
    rad_cf = 0.
!RAD_QR = 0.
!RAD_QS = 0.
    rad_cfd = 0.0_8
    rad_qid = 0.0_8
    rad_qld = 0.0_8
  END IF
END SUBROUTINE RADCOUPLE_D

SUBROUTINE RADCOUPLE_B(te, teb, pl, cf, cfb, af, afb, qclls, qcllsb, &
& qcils, qcilsb, qclan, qclanb, qcian, qcianb, rad_ql, rad_qlb, rad_qi, &
& rad_qib, rad_cf, rad_cfb, rad_rl, rad_rlb, rad_ri, rad_rib, tempor)
  IMPLICIT NONE
!Inputs
  REAL*8, INTENT(IN) :: te, pl, tempor
  REAL*8 :: teb
  REAL*8, INTENT(IN) :: af, cf, qclan, qcian, qclls, qcils
  REAL*8 :: afb, cfb, qclanb, qcianb, qcllsb, qcilsb
! real(8), intent(in ) :: QRN_ALL, QSN_ALL
!Outputs
  REAL*8 :: rad_ql, rad_qi, rad_cf, rad_rl, rad_ri
  REAL*8 :: rad_qlb, rad_qib, rad_cfb, rad_rlb, rad_rib
! real(8), intent(out) :: RAD_QR,RAD_QS
!Locals
  REAL*8 :: ss, rad_ri_an, afx, alph
  REAL*8 :: ssb, rad_ri_anb, afxb
  REAL*8, PARAMETER :: min_ri=20.e-6, max_ri=40.e-6, ri_anv=30.e-6
  INTRINSIC MIN
  INTRINSIC MAX
  INTEGER :: branch
  REAL*8 :: min3
  REAL*8 :: min2
  REAL*8 :: max2b
  REAL*8 :: min1
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  REAL*8 :: min1b
  REAL*8 :: x2
  REAL*8 :: x1
  REAL*8 :: x2b
  REAL*8 :: tempb
  REAL*8 :: max3b
  REAL*8 :: temp
  REAL*8 :: max3
  REAL*8 :: max2
  REAL*8 :: max1
  REAL*8 :: min2b
!Initialize outputs
  rad_ri = 0.0
!RAD_QR = 0.0
!RAD_QS = 0.0
! Adjust Anvil fractions for warm clouds
  alph = 0.1
  ss = (280.-te)/20.
  IF (1.0 .GT. ss) THEN
    CALL PUSHCONTROL1B(0)
    ss = ss
  ELSE
    ss = 1.0
    CALL PUSHCONTROL1B(1)
  END IF
  IF (0.0 .LT. ss) THEN
    CALL PUSHCONTROL1B(0)
    ss = ss
  ELSE
    ss = 0.0
    CALL PUSHCONTROL1B(1)
  END IF
  CALL PUSHREAL8(ss)
  ss = alph + ss**3*(1.0-alph)
  afx = af*ss*0.5
  IF (cf + afx .GT. 1.00) THEN
    rad_cf = 1.00
    CALL PUSHCONTROL1B(0)
  ELSE
    rad_cf = cf + afx
    CALL PUSHCONTROL1B(1)
  END IF
!Total In-cloud liquid
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_ql = (qclls+qclan)/rad_cf
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
    rad_ql = 0.0
  END IF
  IF (rad_ql .GT. 0.01) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
! Total In-cloud ice
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_qi = (qcils+qcian)/rad_cf
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
    rad_qi = 0.0
  END IF
  IF (rad_qi .GT. 0.01) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
! Total In-cloud precipitation
! if (  RAD_CF >0. ) then
!    RAD_QR = ( QRN_ALL ) / RAD_CF
!    RAD_QS = ( QSN_ALL ) / RAD_CF
! else
!    RAD_QR = 0.0
!    RAD_QS = 0.0
! end if
! RAD_QR = MIN( RAD_QR, 0.01 )
! RAD_QS = MIN( RAD_QS, 0.01 )
  IF (pl .LT. 150.) rad_ri = max_ri
  IF (pl .GE. 150.) rad_ri = max_ri*150./pl
! Weigh in a separate R_ice for Anvil Ice according to
  rad_ri_an = rad_ri
  IF (qcils + qcian .GT. 0.0) THEN
    IF (qcils/rad_ri + qcian/ri_anv .GT. 10e-8) THEN
!LINEARIZATION PROTECTION
      rad_ri_an = (qcils+qcian)/(qcils/rad_ri+qcian/ri_anv)
      CALL PUSHCONTROL2B(2)
    ELSE
      CALL PUSHCONTROL2B(1)
    END IF
  ELSE
    CALL PUSHCONTROL2B(0)
  END IF
  IF (rad_ri .GT. rad_ri_an) THEN
    CALL PUSHREAL8(rad_ri)
    rad_ri = rad_ri_an
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL8(rad_ri)
    rad_ri = rad_ri
    CALL PUSHCONTROL1B(1)
  END IF
  IF (rad_ri .LT. min_ri) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 825. .AND. te .LE. 282. .AND. tempor .EQ. 1.) THEN
    IF (0.71*te - 190.25 .LT. 5.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 775. .AND. pl .LT. 825. .AND. te .LE. 282. .AND. te .GT. &
&     275. .AND. tempor .EQ. 1.) THEN
    IF (-(0.1*pl) + 0.71*te - 107.75 .GT. 10.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 825. .AND. te .LE. 275. .AND. tempor .EQ. 1.) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
! Thin low tropical clouds
  IF (pl .GE. 950. .AND. te .GE. 285.) THEN
    IF (2.2*te - 617. .GT. 21.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 925. .AND. te .GE. 290.) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 925. .AND. pl .LT. 950. .AND. te .GT. 285. .AND. te .LT. &
&     290.) THEN
    IF (0.44*pl + 2.2*te - 1035. .GT. 21.) THEN
      CALL PUSHCONTROL1B(0)
      x2 = 21.
    ELSE
      x2 = 0.44*pl + 2.2*te - 1035.
      CALL PUSHCONTROL1B(1)
    END IF
    IF (x2 .LT. 10.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 950. .AND. te .GE. 290.) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (rad_cf .LT. 1.e-5) THEN
    rad_cfb = 0.0_8
    rad_qib = 0.0_8
    rad_qlb = 0.0_8
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rlb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    max3b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      x2b = 0.0_8
    ELSE
      x2b = max3b
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 2.2*x2b
    rad_rlb = 0.0_8
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rlb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    min2b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 2.2*min2b
    rad_rlb = 0.0_8
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rlb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    min1b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 0.71*min1b
    rad_rlb = 0.0_8
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    max2b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 0.71*max2b
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rib = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    CALL POPREAL8(rad_ri)
    rad_ri_anb = rad_rib
  ELSE
    CALL POPREAL8(rad_ri)
    rad_ri_anb = 0.0_8
  END IF
  CALL POPCONTROL2B(branch)
  IF (branch .NE. 0) THEN
    IF (branch .NE. 1) THEN
      temp = qcils/rad_ri + qcian/ri_anv
      tempb1 = rad_ri_anb/temp
      tempb2 = -((qcils+qcian)*tempb1/temp)
      qcilsb = qcilsb + tempb2/rad_ri + tempb1
      qcianb = qcianb + tempb2/ri_anv + tempb1
    END IF
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_qib = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    tempb0 = rad_qib/rad_cf
    qcilsb = qcilsb + tempb0
    qcianb = qcianb + tempb0
    rad_cfb = rad_cfb - (qcils+qcian)*tempb0/rad_cf
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_qlb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    tempb = rad_qlb/rad_cf
    qcllsb = qcllsb + tempb
    qclanb = qclanb + tempb
    rad_cfb = rad_cfb - (qclls+qclan)*tempb/rad_cf
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    afxb = 0.0_8
  ELSE
    cfb = cfb + rad_cfb
    afxb = rad_cfb
  END IF
  afb = afb + 0.5*ss*afxb
  ssb = 0.5*af*afxb
  CALL POPREAL8(ss)
  ssb = (1.0-alph)*3*ss**2*ssb
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) ssb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) ssb = 0.0_8
  teb = teb - ssb/20.
END SUBROUTINE RADCOUPLE_B

end module CLOUDRADCOUPLE
