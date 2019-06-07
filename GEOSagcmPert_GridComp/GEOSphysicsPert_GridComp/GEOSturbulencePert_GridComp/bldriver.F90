MODULE BLDRIVER

!This subroutine takes the model prognostic variables as its inputs and produces
!the three components of the tri-diagonal matrix used to compute BL diffusion. 
!Diffusion coefficient is computed using Louis-Lock formulation.

USE MAPL_ConstantsMod

IMPLICIT NONE

PRIVATE
PUBLIC BL_DRIVER

!Saturation vapor pressure table ESTBLX initialization.
!To be defined for entire module:
integer, parameter :: DEGSUBS    =  100
real(4), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

CONTAINS

subroutine BL_DRIVER ( IM, JM, LM, DT, U, V, TH, Q, P, QIT, QLT, FRLAND, FROCEAN, VARFLT, &
                       ZPBL, CM, CT, CQ, TURBPARAMS, TURBPARAMSI, &
                       USTAR, BSTAR, &
                       AKS, BKS, CKS, AKQ, BKQ, CKQ, AKV, BKV, CKV, EKV, FKV &
                     )

IMPLICIT NONE

!INPUTS
INTEGER, INTENT(IN) :: IM, JM, LM
REAL, INTENT(IN) :: DT

REAL, INTENT(IN), DIMENSION(:) :: TURBPARAMS
INTEGER, INTENT(IN), DIMENSION(:) :: TURBPARAMSI

REAL, INTENT(IN), DIMENSION(IM,JM,LM) :: U, V, TH, Q, QIT, QLT
REAL, INTENT(IN), DIMENSION(IM,JM,0:LM) :: P
REAL, INTENT(IN), DIMENSION(IM,JM) :: FRLAND, VARFLT, FROCEAN

REAL, INTENT(IN), DIMENSION(IM,JM) :: CM, CQ

REAL, INTENT(IN), DIMENSION(IM,JM)    :: USTAR, BSTAR

!INOUTS
REAL, INTENT(INOUT), DIMENSION(IM,JM) :: ZPBL, CT

!OUTPUTS
REAL, INTENT(OUT), DIMENSION(IM,JM,LM) :: AKS, BKS, CKS
REAL, INTENT(OUT), DIMENSION(IM,JM,LM) :: AKQ, BKQ, CKQ
REAL, INTENT(OUT), DIMENSION(IM,JM,LM) :: AKV, BKV, CKV, EKV, FKV

!LOCALS
INTEGER :: I,J
REAL, DIMENSION(IM,JM,LM) :: Pf, PIf, QI, QL, TV, THV, Zf, DMI
REAL, DIMENSION(IM,JM,0:LM) :: PI, Z
REAL, DIMENSION(IM,JM,1:LM-1) :: RDZ
REAL, DIMENSION(IM,JM,LM) :: T
REAL, DIMENSION(IM,JM,0:LM) :: Km, Kh

REAL, DIMENSION(IM,JM,0:LM) :: Ri, DU
REAL, DIMENSION(IM,JM,0:LM) :: ALH_X, KMLS_X, KHLS_X

REAL, DIMENSION(IM,JM,LM) :: RADLW !dh Note that this is really an input but we would rather not put it in the traj
                                   !so for now do not use. Later on add approximate numbers.

!Turb constants
REAL :: LOUIS,LAMBDAM,LAMBDAM2,LAMBDAH,LAMBDAH2,ZKMENV,ZKHENV,MINTHICK,MINSHEAR,C_B,LAMBDA_B,AKHMMAX
REAL :: PRANDTLSFC,PRANDTLRAD,BETA_RAD,BETA_SURF,KHRADFAC,KHSFCFAC,TPFAC_SURF,ENTRATE_SURF,PCEFF_SURF,LOUIS_MEMORY
INTEGER :: KPBLMIN, LOCK_ON, PBLHT_OPTION, RADLW_DEP

!QSATVP TABLE
integer, parameter :: DEGSUBS    =  100
real(4), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1
real(4), dimension(TABLESIZE) :: ESTBLX

 !Compute saturation vapour pressure lookup table
 call ESINIT(ESTBLX)

 !Initialize outputs
 AKS = 0.0
 BKS = 0.0
 CKS = 0.0
 AKQ = 0.0
 BKQ = 0.0
 CKQ = 0.0
 AKV = 0.0
 BKV = 0.0
 CKV = 0.0
 EKV = 0.0
 FKV = 0.0

 !Caclulate temperature from Exner pressure.
 Pf = 0.5*(P(:,:,0:LM-1) + P(:,:,1:LM))
 PIf = (Pf/MAPL_P00)**(MAPL_RGAS/MAPL_CP)
 T = PIf*TH

 !TurbParams, this goes at the level above.
 !TurbParams(1)  = 5.0          !LOUIS
 !TurbParams(2)  = 160.0        !LAMBDAM
 !TurbParams(3)  = 1.0          !LAMBDAM2
 !TurbParams(4)  = 160.0        !LAMBDAH
 !TurbParams(5)  = 1.0          !LAMBDAH2
 !TurbParams(6)  = 3000.        !ZKMENV
 !TurbParams(7)  = 3000.        !ZKHENV
 !TurbParams(8)  = 0.1          !MINTHICK
 !TurbParams(9)  = 0.0030       !MINSHEAR
 !TurbParams(10) = 2.5101471e-8 !C_B
 !TurbParams(11) = 1500.        !LAMBDA_B 
 !TurbParams(12) = 500.         !AKHMMAX
 !TurbParams(13) = 1.0          !PRANDTLSFC
 !TurbParams(14) = 0.75         !PRANDTLRAD
 !TurbParams(15) = 0.50         !BETA_RAD
 !TurbParams(16) = 0.25         !BETA_SURF
 !TurbParams(17) = 0.85         !KHRADFAC
 !TurbParams(18) = 0.45         !KHSFCFAC
 !TurbParams(19) = 20.0         !TPFAC_SURF
 !TurbParams(20) = 1.5e-3       !ENTRATE_SURF
 !TurbParams(21) = 0.5          !PCEFF_SURF
 !TurbParams(22) = -999.        !LOUIS_MEMORY

 !TurbParamsI(1) = count(PREF < 50000.)   !KPBLMIN
 !TurbParamsI(2) = 1                      !LOCK_ON
 !TurbParamsI(3) = 1                      !PBLHT_OPTION 
 !TurbParamsI(4) = 0                      !RADLW_DEP

 !Turbulence constants
 LOUIS        = TurbParams(1)
 LAMBDAM      = TurbParams(2)
 LAMBDAM2     = TurbParams(3) 
 LAMBDAH      = TurbParams(4)
 LAMBDAH2     = TurbParams(5)
 ZKMENV       = TurbParams(6)
 ZKHENV       = TurbParams(7)
 MINTHICK     = TurbParams(8)
 MINSHEAR     = TurbParams(9)
 C_B          = TurbParams(10)
 LAMBDA_B     = TurbParams(11)
 AKHMMAX      = TurbParams(12)
 PRANDTLSFC   = TurbParams(13)
 PRANDTLRAD   = TurbParams(14)
 BETA_RAD     = TurbParams(15)
 BETA_SURF    = TurbParams(16)
 KHRADFAC     = TurbParams(17)
 KHSFCFAC     = TurbParams(18)
 TPFAC_SURF   = TurbParams(19)
 ENTRATE_SURF = TurbParams(20)
 PCEFF_SURF   = TurbParams(21)
 LOUIS_MEMORY = TurbParams(22)

 KPBLMIN      = TurbParamsI(1)
 LOCK_ON      = TurbParamsI(2)
 PBLHT_OPTION = TurbParamsI(3)
 RADLW_DEP    = TurbParamsI(4)

 !Initialise arrays
 KH    = 0.0
 KM    = 0.0
 RI    = 0.0
 DU    = 0.0

 !If over the ocean then Ts' = 0. That is the perturbation of surface temperature
 !is negligible. If over land then it is not. But temperature transfer is quick so
 !T1'-Ts' = 0. So set CT to 0 over land and 1 over ocean. 
 !In between use fraction to take percentage of increment.
 DO I = 1,IM
    DO J = 1,JM
       IF (FROCEAN(I,J).eq.1.0) then
          CT(I,J) = FROCEAN(I,J) * CT(I,J)
       endif
    enddo
 enddo

 !Call the subroutine to get levels, RDZ and DMI
 call PRELIMINARY( IM*JM      , & !IN
                   LM         , & !IN
                   DT         , & !IN
                   T          , & !IN
                   Q          , & !IN
                   P          , & !IN
                   TH         , & !IN
                   QIT        , & !IN
                   QLT        , & !IN
                   Zf         , & !OUT
                   Z          , & !OUT
                   TV         , & !OUT
                   THV        , & !OUT
                   RDZ        , & !OUT
                   DMI        , & !OUT
                   Pf           & !OUT
                 )

 !GET LOUIS DIFFUSIVITIES
 call LOUIS_DIFF( IM*JM       , & !IN  
                  LM          , & !IN
                  KH          , & !OUT
                  KM          , & !OUT
                  RI          , & !OUT
                  DU          , & !OUT
                  ZPBL        , & !IN
                  Zf          , & !IN
                  Z           , & !IN
                  THV         , & !IN
                  U           , & !IN
                  V           , & !IN
                  LOUIS       , & !IN
                  MINSHEAR    , & !IN
                  MINTHICK    , & !IN
                  LAMBDAM     , & !IN
                  LAMBDAM2    , & !IN
                  LAMBDAH     , & !IN
                  LAMBDAH2    , & !IN
                  ZKMENV      , & !IN
                  ZKHENV      , & !IN
                  AKHMMAX     , & !IN
                  ALH_X       , & !OUT
                  KMLS_X      , & !OUT
                  KHLS_X        & !OUT
                )

 !Optional - overwrite diffusivities using the Lock et al scheme.
 CALL LOCK_DIFF( IM*JM        , & !IN
                 LM           , & !IN
                 RADLW        , & !IN
                 USTAR        , & !IN
                 BSTAR        , & !IN
                 FRLAND       , & !IN
                 T            , & !IN
                 Q            , & !IN
                 QIT          , & !IN
                 QLT          , & !IN
                 U            , & !IN
                 V            , & !IN
                 Zf           , & !IN
                 Pf           , & !IN
                 Z            , & !IN
                 P            , & !IN
                 KM           , & !INOUT
                 KH           , & !INOUT
                 PRANDTLSFC   , & !IN
                 PRANDTLRAD   , & !IN
                 BETA_SURF    , & !IN
                 BETA_RAD     , & !IN
                 TPFAC_SURF   , & !IN
                 ENTRATE_SURF , & !IN
                 PCEFF_SURF   , & !IN
                 KHRADFAC     , & !IN
                 KHSFCFAC     , & !IN
                 RADLW_DEP    , & !IN
                 ESTBLX         & !IN
               )
 
 call TRIDIAG_SETUP( IM*JM          , & !IN
                     LM             , & !IN
                     DT             , & !IN
                     KPBLMIN        , & !IN
                     Zf             , & !IN
                     Pf             , & !IN
                     RDZ            , & !IN
                     DMI            , & !IN
                     P              , & !IN
                     TV             , & !IN
                     CT             , & !IN
                     CQ             , & !IN
                     CM             , & !IN
                     T              , & !IN
                     Q              , & !IN
                     TH             , & !IN
                     U              , & !IN
                     V              , & !IN
                     KH             , & !INOUT
                     KM             , & !INOUT
                     AKQ, AKS, AKV  , & !OUT
                     BKQ, BKS, BKV  , & !OUT
                     CKQ, CKS, CKV  , & !OUT
                               EKV  , & !OUT
                     ZPBL             & !OUT
                   )

 call ORODRAG( IM*JM          , & !IN
               LM             , & !IN
               DT             , & !IN
               LAMBDA_B       , & !IN
               C_B            , & !IN
               FKV            , & !OUT
               BKV            , & !INOUT
               U              , & !IN
               V              , & !IN
               Zf             , & !IN
               VARFLT         , & !IN
               P                & !IN
             )


end subroutine BL_DRIVER


subroutine PRELIMINARY( IRUN, LM, DT, T, QV, PHALF, TH, QIT, QLT, &
                        ZFULL, ZHALF, TV, PV, RDZ, DMI, PFULL &
                      )

 !Subroutine to set up priliminary

 IMPLICIT NONE

 !INTPUTS
 integer, intent(IN)                    :: IRUN, LM
 real, intent(IN)                       :: DT
 real, intent(IN), dimension(IRUN,LM)   :: T, QV, TH, QIT, QLT
 real, intent(IN), dimension(IRUN,LM+1) :: PHALF

 !OUTPUTS
 real, intent(OUT), dimension(IRUN,  LM  ) :: ZFULL, TV, PV, DMI, PFULL
 real, intent(OUT), dimension(IRUN,  LM+1) :: ZHALF
 real, intent(OUT), dimension(IRUN,1:LM-1) :: RDZ

 !LOCALS
 integer :: I,L
 real :: PKE_atL, PKE_atLp1, TVE
 real :: QL_tot, QI_tot 

 do I = 1, IRUN

    !Compute the edge heights using Arakawa-Suarez hydrostatic equation
    ZHALF(I,LM+1) = 0.0
    do L = LM, 1, -1
       PKE_atLp1  = (PHALF(I,L+1)/MAPL_P00)**MAPL_KAPPA
       PKE_atL    = (PHALF(I,L  )/MAPL_P00)**MAPL_KAPPA

       ZHALF(I,L) = ZHALF(I,L+1) + (MAPL_CP/MAPL_GRAV)*TH(I,L)*(PKE_atLp1-PKE_atL)
    end do

    !Layer height, pressure, and virtual temperatures
    do L = 1, LM
       QL_tot = QIT(I,L)
       QI_tot = QLT(I,L)

       ZFULL(I,L) = 0.5*(ZHALF(I,L)+ZHALF(I,L+1))
       PFULL(I,L) = 0.5*(PHALF(I,L)+PHALF(I,L+1))

       TV(I,L)  = T(I,L) *( 1.0 + MAPL_VIREPS *QV(I,L) - QL_tot - QI_tot ) 
       PV(I,L) = TV(I,L)*(TH(I,L)/T(I,L))
    end do

    do L = 1, LM-1
       TVE = (TV(I,L) + TV(I,L+1))*0.5
  
       ! Miscellaneous factors
       RDZ(I,L) = PHALF(I,L+1) / ( MAPL_RGAS * TVE )
       RDZ(I,L) = RDZ(I,L) / (ZFULL(I,L)-ZFULL(I,L+1))
    end do

    do L = 1, LM
       DMI(I,L) = (MAPL_GRAV*DT)/(PHALF(I,L+1)-PHALF(I,L))
    end do

    !Running 1-2-1 smooth of bottom 5 levels of Virtual Pot. Temp.
    PV(I,LM  ) = PV(I,LM-1)*0.25 + PV(I,LM  )*0.75
    PV(I,LM-1) = PV(I,LM-2)*0.25 + PV(I,LM-1)*0.50 + PV(I,LM  )*0.25 
    PV(I,LM-2) = PV(I,LM-3)*0.25 + PV(I,LM-2)*0.50 + PV(I,LM-1)*0.25 
    PV(I,LM-3) = PV(I,LM-4)*0.25 + PV(I,LM-3)*0.50 + PV(I,LM-2)*0.25 
    PV(I,LM-4) = PV(I,LM-5)*0.25 + PV(I,LM-4)*0.50 + PV(I,LM-3)*0.25 
    PV(I,LM-5) = PV(I,LM-6)*0.25 + PV(I,LM-5)*0.50 + PV(I,LM-4)*0.25 

 end do 

end subroutine PRELIMINARY



subroutine LOUIS_DIFF( IRUN,LM,KH,KM,RI,DU,ZPBL,ZZ,ZE,PV,UU,VV,LOUIS,MINSHEAR,MINTHICK,&
                       LAMBDAM,LAMBDAM2,LAMBDAH,LAMBDAH2,ZKMENV,ZKHENV,AKHMMAX,ALH_DIAG,KMLS_DIAG,KHLS_DIAG)

 !Subroutine to compute eddy diffussivity using the Louis (1979) method.
 !Simplifications are made compared with the NL stability functions f(Ri). 
 !Gradient is very steep for near neutral boundary layer conditions.

 IMPLICIT NONE

 !INPUTS
 integer, intent(IN   ) :: IRUN,LM
 real, intent(IN   ) :: LOUIS         ! Louis scheme parameters (usually 5).
 real, intent(IN   ) :: MINSHEAR      ! Min shear allowed in Ri calculation (s-1).
 real, intent(IN   ) :: MINTHICK      ! Min layer thickness (m).
 real, intent(IN   ) :: AKHMMAX       ! Maximum allowe diffusivity (m+2 s-1).
 real, intent(IN   ) :: LAMBDAM       ! Blackadar(1962) length scale parameter for momentum (m).
 real, intent(IN   ) :: LAMBDAM2      ! Second Blackadar parameter for momentum (m).
 real, intent(IN   ) :: LAMBDAH       ! Blackadar(1962) length scale parameter for heat (m).
 real, intent(IN   ) :: LAMBDAH2      ! Second Blackadar parameter for heat (m).
 real, intent(IN   ) :: ZKMENV        ! Transition height for Blackadar param for momentum (m)
 real, intent(IN   ) :: ZKHENV        ! Transition height for Blackadar param for heat     (m)
 real, intent(IN   ) :: ZZ(IRUN,LM)   ! Height of layer center above the surface (m).
 real, intent(IN   ) :: PV(IRUN,LM)   ! Virtual potential temperature at layer center (K).
 real, intent(IN   ) :: UU(IRUN,LM)   ! Eastward velocity at layer center (m s-1).
 real, intent(IN   ) :: VV(IRUN,LM)   ! Northward velocity at layer center (m s-1).
 real, intent(IN   ) :: ZE(IRUN,LM+1) ! Height of layer base above the surface (m).
 real, intent(IN   ) :: ZPBL(IRUN)    ! PBL Depth (m)

 !OUTPUTS
 !These are 1:LM+1 here but 0:LM in the GC - MAYBE JUST SORT THIS OUT?
 !Old code only passed in 1:LM-1 from GC which is 2:LM here.
 real, intent(  OUT) :: KM(IRUN,LM+1) ! Momentum diffusivity at base of each layer (m+2 s-1).
 real, intent(  OUT) :: KH(IRUN,LM+1) ! Heat diffusivity at base of each layer  (m+2 s-1).
 real, intent(  OUT) :: RI(IRUN,LM+1) ! Richardson number
 real, intent(  OUT) :: DU(IRUN,LM+1) ! Magnitude of wind shear (s-1).
 real, intent(  OUT) :: ALH_DIAG(IRUN,LM+1)  ! Blackadar Length Scale diagnostic (m) [Optional] 
 real, intent(  OUT) :: KMLS_DIAG(IRUN,LM+1) ! Momentum diffusivity at base of each layer (m+2 s-1).
 real, intent(  OUT) :: KHLS_DIAG(IRUN,LM+1) ! Heat diffusivity at base of each layer  (m+2 s-1).

 !LOCALS
 integer :: L,Levs,i,j,k
 real    :: ALH, ALM, DZ, DT, TM, PS, LAMBDAM_X, LAMBDAH_X
 real    :: pbllocal, alhfac,almfac

 almfac = 1.2
 alhfac = 1.2

 do i = 1, irun !Loop of columns

    ALH_DIAG(i,   1) = 0.0
    ALH_DIAG(i,LM+1) = 0.0

    KH(i,   1) = 0.0
    KH(i,LM+1) = 0.0
    KM(i,   1) = 0.0
    KM(i,LM+1) = 0.0
    DU(i,   1) = 0.0
    DU(i,LM+1) = 0.0
    RI(i,   1) = 0.0
    RI(i,LM+1) = 0.0
    KMLS_DIAG(i,   1) = 0.0
    KMLS_DIAG(i,LM+1) = 0.0
    KHLS_DIAG(i,   1) = 0.0
    KHLS_DIAG(i,LM+1) = 0.0

    pbllocal = ZPBL(i)

    if ( pbllocal .LE. ZZ(i,LM) ) then
       pbllocal = ZZ(i,LM)
    endif

    do k = 2, lm !This is equivalent to 1->LM-1 but we move all 0:72 elements
       KH(i,k) = 0.0
       KM(i,k) = 0.0
       DU(i,k) = 0.0
       RI(i,k) = 0.0

       DZ      = (ZZ(i,k-1) - ZZ(i,k))
       TM      = (PV(i,k-1) + PV(i,k))*0.5
       DT      = (PV(i,k-1) - PV(i,k))
       DU(i,k) = (UU(i,k-1) - UU(i,k))**2 + (VV(i,k-1) - VV(i,k))**2

       !Limits distance between layer centers and vertical shear at edges.
       DZ      = max(DZ, MINTHICK)
       DU(i,k) = sqrt(DU(i,k))/DZ

       !Calc Richardson
       RI(i,k) = MAPL_GRAV*(DT/DZ)/(TM*( max(DU(i,k), MINSHEAR)**2))

       !Blackadar (1962) length scale
       LAMBDAM_X = MAX( 0.1 * PBLLOCAL * EXP( -(ZE(i,k) / ZKMENV )**2 ) , LAMBDAM2 )
       LAMBDAH_X = MAX( 0.1 * PBLLOCAL * EXP( -(ZE(i,k) / ZKHENV )**2 ) , LAMBDAH2 )

       ALM = almfac * ( MAPL_KARMAN*ZE(i,k)/( 1.0 + MAPL_KARMAN*(ZE(i,k)/LAMBDAM_X) ) )**2
       ALH = alhfac * ( MAPL_KARMAN*ZE(i,k)/( 1.0 + MAPL_KARMAN*(ZE(i,k)/LAMBDAH_X) ) )**2

       ALH_DIAG(i,k) = SQRT( ALH )

       !Unstable convectve boundary layer
       if ( RI(i,k)  < 0.0 ) then
          PS = ( (ZZ(i,k-1)/ZZ(i,k))**(1./3.) - 1.0 ) ** 3
          PS = ALH*sqrt( PS/(ZE(i,k)*(DZ**3)) )
          PS = RI(i,k)/(1.0 + (3.0*LOUIS*LOUIS)*PS*sqrt(-RI(i,k)))
      
          KH(i,k) = 1.0 - (LOUIS*3.0)*PS
          KM(i,k) = 1.0 - (LOUIS*2.0)*PS
       end if

       !Stable boundary layer
       if ( RI(i,k) >= 0.0 ) then
          PS      = sqrt  (1.0 +  LOUIS     *RI(i,k)   )
         
          KH(i,k) = 1.0 / (1.0 + (LOUIS*3.0)*RI(i,k)*PS)
          KM(i,k) = PS  / (PS  + (LOUIS*2.0)*RI(i,k)   )
       end if

       !Dimensionalize Kz and limit diffusvity
       ALM = DU(i,k)*ALM
       ALH = DU(i,k)*ALH

       KM(i,k) = min(KM(i,k)*ALM, AKHMMAX)
       KH(i,k) = min(KH(i,k)*ALH, AKHMMAX)
       KMLS_DIAG(i,k) = KM(i,k)
       KHLS_DIAG(i,k) = KH(i,k)

    end do !Loop over levels

 end do !Loop over columns.

end subroutine LOUIS_DIFF

subroutine TRIDIAG_SETUP(IRUN,LM,DT,KPBLMIN,ZFULL,PFULL,RDZ,DMI,PHALF,&
                         TV,CT,CQ,CU,T,Q,TH,U,V,DIFF_T,DIFF_M,&
                         AKQ,AKS,AKV,BKQ,BKS,BKV,CKQ,CKS,CKV,EKV,ZPBL)

 IMPLICIT NONE

 !INPUTS
 integer, intent(IN)                          :: IRUN,LM,KPBLMIN
 real,    intent(IN)                          :: DT
 real,    intent(IN),    dimension(IRUN,LM)   :: T,Q,TH,U,V,ZFULL,PFULL,DMI,TV
 real,    intent(IN),    dimension(IRUN,LM-1) :: RDZ
 real,    intent(IN),    dimension(IRUN,LM+1) :: PHALF
 real,    intent(IN),    dimension(IRUN  )    :: CT, CQ, CU

 !INOUTS
 real,    intent(INOUT), dimension(IRUN,LM+1) :: DIFF_T, DIFF_M

 !OUTPUTS
 real,    intent(OUT),   dimension(IRUN,LM)   :: AKQ, AKS, AKV
 real,    intent(OUT),   dimension(IRUN,LM)   :: BKQ, BKS, BKV
 real,    intent(OUT),   dimension(IRUN,LM)   :: CKQ, CKS, CKV
 real,    intent(OUT),   dimension(IRUN,LM)   ::           EKV
 real,    intent(OUT),   dimension(IRUN)      :: ZPBL

 !LOCALS
 real, dimension(LM+1) :: temparray
 integer :: I,L,locmax
 real    :: maxkh,minlval,thetavs,thetavh,uv2h
 real, parameter :: tcri_crit = 0.25 !PBL-top diagnostic
 
 real,             dimension(LM ) :: tcrib

 do I = 1, IRUN

    !Update zpbl
    ZPBL(I) = 10.e15
    do L=LM,2,-1
       if ( (DIFF_T(I,L) < 2.) .and. (DIFF_T(I,L+1) >= 2.) .and. (ZPBL(I) == 10.e15 ) ) then
          ZPBL(I) = ZFULL(I,L) 
       end if 
    end do
    if (  ZPBL(I) .eq. 10.e15 ) then
       ZPBL(I) = ZFULL(I,LM)
    endif
    ZPBL(I) = MIN(ZPBL(I),ZFULL(I,KPBLMIN))

    !Second difference coefficients for scalars; RDZ is RHO/DZ, DMI is (G DT)/DP
    do L = 1, LM
       if (L <= LM-1) then
          CKS(I,L) = -DIFF_T(I,L+1) * RDZ(I,L)
          if (L == 1) then
             AKS(I,L) = 0.0
          endif
          AKS(I,L+1) = CKS(I,L) * DMI(I,L+1)
          CKS(I,L) = CKS(I,L) * DMI(I,L)
       end if

       if (L == LM) then
          CKS(I,L) = -CT(I) * DMI(I,L)
       endif

       !Fill DIFF_T at level LM+1 with CT * RDZ for diagnostic output
       if (L == LM) then
          DIFF_T(I,L+1) = CT(I) * (PHALF(I,L+1)/(MAPL_RGAS * TV(I,L))) / ZFULL(I,L)
       endif

       !Water vapor can differ at the surface
       CKQ(I,L) = CKS(I,L)
       AKQ(I,L) = AKS(I,L)
       if (L == LM) then
          CKQ(I,L) = -CQ(I) * DMI(I,L)
       endif

       !Second difference coefficients for winds
       !EKV is saved to use in the frictional heating calc.
       if (L <= LM-1) then
          EKV(I,L) = -DIFF_M(I,L+1) * RDZ(I,L)
          if (L == 1) then
             AKV(I,L) = 0.0
          endif
          AKV(I,L+1) = EKV(I,L) * DMI(I,L+1)
          CKV(I,L) = EKV(I,L) * DMI(I,L)
          EKV(I,L) = -MAPL_GRAV * EKV(I,L)
       end if

       if (L == LM) then
          CKV(I,L) = -  CU(I) * DMI(I,L)
          EKV(I,L) =  MAPL_GRAV *  CU(I)
       end if

       !Fill DIFF_M at level LM+1 with CU * RDZ for diagnostic output
       if (L == LM) then
          DIFF_M(I,L+1) = CU(I) * (PHALF(I,L+1)/(MAPL_RGAS * TV(I,L))) / ZFULL(I,L)
       endif

       !Setup the tridiagonal matrix
       BKS(I,L) = 1.00 - (AKS(I,L)+CKS(I,L))
       BKQ(I,L) = 1.00 - (AKQ(I,L)+CKQ(I,L))
       BKV(I,L) = 1.00 - (AKV(I,L)+CKV(I,L))

    end do

 end do 
 

end subroutine TRIDIAG_SETUP


subroutine ORODRAG(IRUN,LM,DT,LAMBDA_B,C_B,FKV,BKV,U,V,ZFULL,VARFLT,PHALF)

 !Subroutine to compute orograpic drag, Beljaars (2003).

 IMPLICIT NONE

 !INPUTS
 integer, intent(IN)                        :: IRUN, LM
 real,    intent(IN)                        :: DT, LAMBDA_B, C_B
 real,    intent(IN), dimension(IRUN)       :: VARFLT 
 real,    intent(IN), dimension(IRUN,LM)    :: U, V, ZFULL
 real,    intent(IN), dimension(IRUN,LM+1)  :: PHALF

 !INOUTS
 real,    intent(INOUT), dimension(IRUN,LM) :: BKV

 !OUTPUTS
 real,    intent(OUT), dimension(IRUN,LM) :: FKV
      
 !LOCALS
 integer :: I, L
 real    :: FKV_temp

 do I = 1, IRUN
    do L=LM,1,-1

       FKV(I,L) = 0.0

       if (ZFULL(I,L) < 4.0*LAMBDA_B) then
          FKV_temp = ZFULL(I,L)*(1.0/LAMBDA_B)
          FKV_temp = VARFLT(I) * exp(-FKV_temp*sqrt(FKV_temp))*(FKV_temp**(-1.2))
          FKV_temp = (C_B/LAMBDA_B)*min( sqrt(U(I,L)**2+V(I,L)**2),5.0 )*FKV_temp

          BKV(I,L) = BKV(I,L) + DT*FKV_temp
          FKV(I,L) = FKV_temp * (PHALF(I,L+1)-PHALF(I,L))
       end if
    end do

 end do 

end subroutine ORODRAG


subroutine LOCK_DIFF( ncol,nlev,tdtlw_in_dev,u_star_dev,b_star_dev,frland_dev,t_dev,qv_dev, &
                      qit_dev,qlt_dev,u_dev,v_dev,zfull_dev,pfull_dev,zhalf_dev,phalf_dev, &
                      !INOUTS
                      diff_m_dev,diff_t_dev, &
                      !Constants
                      prandtlsfc_const,prandtlrad_const,beta_surf_const,beta_rad_const, &
                      tpfac_sfc_const,entrate_sfc_const,pceff_sfc_const,khradfac_const,khsfcfac_const,radlw_dep,&
                      ESTBLX )

 IMPLICIT NONE

 !INPUTS
 integer, intent(in)                               :: ncol,nlev    
 real,    intent(in)                               :: ESTBLX(:) 
 real,    intent(in),    dimension(ncol)           :: u_star_dev,b_star_dev,frland_dev
 real,    intent(in),    dimension(ncol,nlev)      :: tdtlw_in_dev,t_dev,qv_dev,qit_dev,qlt_dev
 real,    intent(in),    dimension(ncol,nlev)      :: u_dev,v_dev,zfull_dev,pfull_dev
 real,    intent(in),    dimension(ncol,1:nlev+1)  :: zhalf_dev, phalf_dev ! 0:72 in GC, 1:73 here.
 real,    intent(in)                               :: prandtlsfc_const,prandtlrad_const,beta_surf_const,beta_rad_const
 real,    intent(in)                               :: khradfac_const,tpfac_sfc_const,entrate_sfc_const
 real,    intent(in)                               :: pceff_sfc_const,khsfcfac_const

 !INOUTS
 real,    intent(inout), dimension(ncol,1:nlev+1)  :: diff_m_dev,diff_t_dev

 !LOCALS
 real, parameter :: akmax = 1.e4, zcldtopmax = 3.e3, ramp = 20.
 integer :: i,j,k,ibot,itmp,ipbl,kmax,kcldtop,kcldbot,kcldtop2,radlw_dep
 logical :: used,keeplook,do_jump_exit
 real :: maxradf, tmpradf,stab
 real :: maxqdot, tmpqdot,wentrmax
 real :: maxinv, qlcrit,ql00,qlm1,Abuoy,Ashear
 real :: wentr_tmp,hlf,vbulkshr,vbulk_scale
 real :: k_entr_tmp,tmpjump,critjump,radperturb,buoypert
 real :: tmp1, tmp2, slmixture
 real :: vsurf3, vshear3,vrad3, vbr3,dsiems
 real :: dslvcptmp,ztmp
 real :: zradtop
 real :: vrad,svpcp,chis,vbrv
 real :: vsurf, qs_dummy, ptmp
 real :: wentr_rad,wentr_pbl,wentr_brv
 real :: convpbl, radpbl
 real, dimension(nlev) :: slv,density,qc,slvcp,hleff,radfq,pblfq,k_m_troen,k_t_troen,k_rad,dqs
 
 !Moved to local as not needed in outputs, but still used in calculations
 real, dimension(ncol,1:nlev+1)  :: k_m_entr_dev,k_t_entr_dev
 real, dimension(ncol)  :: zsml_dev,zradml_dev,zcloud_dev,zradbase_dev

 qlcrit     = 1.0e-6
 Abuoy      = 0.23
 Ashear     = 25.0
 wentrmax   = 0.05

 ibot = nlev

 do i = 1, ncol

    zradml_dev(i)     = 0.0
    zcloud_dev(i)     = 0.0
    zradbase_dev(i)   = 0.0
    zsml_dev(i)       = 0.0 ! note that this must be zero, indicates stable surface layer is used in gravity wave drag scheme

    zradtop    = 0.0
    convpbl    = 0.0
    wentr_pbl  = 0.0
    vsurf      = 0.0
    radpbl     = 0.0
    svpcp      = 0.0
    vrad       = 0.0
    wentr_rad  = 0.0

    do k = 1, nlev
       pblfq(k)      = 0.0
       k_t_troen(k)  = 0.0
       k_m_troen(k)  = 0.0
       radfq(k)      = 0.0
       k_rad(k)      = 0.0
       k_t_entr_dev(i,k)   = 0.0
       k_m_entr_dev(i,k)   = 0.0

    end do

    !For GPU we must zero out even the LM+1'th position

    k_t_entr_dev(i,nlev+1)   = 0.0
    k_m_entr_dev(i,nlev+1)   = 0.0

    !set up specific humidities and static energies and compute airdensity
    do k = 1, nlev
       if ( t_dev(i,k) <= MAPL_TICE-ramp ) then
          HLEFF(k) = MAPL_ALHS
       else if ( (t_dev(i,k) > MAPL_TICE-ramp) .and. (t_dev(i,k) < MAPL_TICE) ) then
          HLEFF(k) =  ( (t_dev(i,k)-MAPL_TICE+ramp)*MAPL_ALHL +  (MAPL_TICE -t_dev(i,k)    )*MAPL_ALHS   ) / ramp
       else
          HLEFF(k) = MAPL_ALHL
       end if

       !Compute sat specific humidity and deriv, Qs not actually needed
       !dqs(k) = dqsat(t_dev(i,k), pfull_dev(i,k), qsat=qs(k), pascals=.true. )
       !dqs(k) = dqsat(t_dev(i,k), pfull_dev(i,k), qsat=qs_dummy, pascals=.true. )

       ptmp = pfull_dev(i,k)*0.01
       call DQSAT_sub_sca(dqs(k),qs_dummy,t_dev(i,k),ptmp,ESTBLX)

       !Compute total cloud condensate - qc. These are grid box mean values.
       !qc(k) = 0.0 !(qlls_dev(i,k) + qils_dev(i,k))
       qc(k) = (qit_dev(i,k) + qlt_dev(i,k))

       !Compute liquid static energy.
       slv(k) = MAPL_CP*t_dev(i,k)*(1+MAPL_VIREPS*qv_dev(i,k)-qc(k)) + MAPL_GRAV*zfull_dev(i,k) - hleff(k)*qc(k)

       density(k) = pfull_dev(i,k)/(MAPL_RGAS*(t_dev(i,k) *(1.+MAPL_VIREPS*qv_dev(i,k)-qc(k))))              

    end do

    !Reset indices
    ipbl    = -1
    kcldtop = -1

    !Surface driven convective layers (if b_star_dev > 0., that is,upward surface buoyancy flux)
    if (b_star_dev(i) .gt. 0.) then

       !Find depth of surface driven mixing by raising a parcel from the surface with some excess buoyancy to its level 
       !of neutral buoyancy. Note the use slv as the density variable permits one to goes through phase changes to find parcel top
       call mpbl_depth(i,ncol,nlev,tpfac_sfc_const,entrate_sfc_const,pceff_sfc_const, &
                       t_dev,qv_dev,u_dev,v_dev,zfull_dev,pfull_dev,b_star_dev,u_star_dev,ipbl,zsml_dev,&
                       ESTBLX)
   
       !Define velocity scales vsurf and vshear
       vsurf3   = u_star_dev(i)*b_star_dev(i)*zsml_dev(i) !cubed
       vshear3  = Ashear*u_star_dev(i)*u_star_dev(i)*u_star_dev(i)
   
       vsurf = vsurf3**(1./3.)
   
       !Define velocity scale vbulkshr      
       vbulkshr   =  sqrt ( ( u_dev(i,ipbl)-u_dev(i,ibot) )**2 + ( v_dev(i,ipbl)-v_dev(i,ibot) )**2  ) /  zsml_dev(i)
       vbulk_scale = 3.0 / 1000.   ! Non-local (PBL-deep) shear scale

       !Following Lock et al. 2000, limit height of surface well mixed layer if interior stable interface is
       !found. An interior stable interface is diagnosed if the slope between 2 full levels is greater than critjump
       critjump = 2.0
       if (ipbl .lt. ibot) then 
          do k = ibot, ipbl+1, -1
             tmpjump =(slv(k-1)-slv(k))/MAPL_CP 
             if (tmpjump .gt. critjump) then
                ipbl = k
                zsml_dev(i) = zhalf_dev(i,ipbl)
                exit
             end if
          enddo
       end if

       !Compute entrainment rate
       tmp1 = MAPL_GRAV*max(0.1,(slv(ipbl-1)-slv(ipbl))/ MAPL_CP)/(slv(ipbl)/MAPL_CP)
       tmp2 = ((vsurf3+vshear3)**(2./3.)) / zsml_dev(i)

       wentr_tmp = min( wentrmax,  max(0., (beta_surf_const * (vsurf3 + vshear3)/zsml_dev(i))/(tmp1+tmp2) ) )

       !Fudgey adjustment of entrainment to reduce it for shallow boundary layers, increase for deep ones
       if ( zsml_dev(i) .lt. 1600. ) then 
          wentr_tmp = wentr_tmp * ( zsml_dev(i) / 800. )
       else
          wentr_tmp = 2.*wentr_tmp
       endif

       !More fudgey adjustment of entrainment.
       k_entr_tmp = wentr_tmp*(zfull_dev(i,ipbl-1)-zfull_dev(i,ipbl))  
       k_entr_tmp = min ( k_entr_tmp, akmax )

       do k = ipbl, ibot
          pblfq(k)  = 1.
       end do
       convpbl           = 1.
       wentr_pbl         = wentr_tmp
       k_t_troen(ipbl)   = k_entr_tmp
       k_m_troen(ipbl)   = k_entr_tmp
       k_t_entr_dev (i,ipbl) = k_t_entr_dev(i,ipbl) + k_entr_tmp
       k_m_entr_dev (i,ipbl) = k_m_entr_dev(i,ipbl) + k_entr_tmp

       !Compute diffusion coefficients in the interior of the PBL

       if (ipbl .lt. ibot) then
   
          call diffusivity_pbl2(i, ncol, nlev,zsml_dev,khsfcfac_const, k_entr_tmp,vsurf,&
                                frland_dev,zhalf_dev,k_m_troen,k_t_troen)

          do k = ipbl+1, ibot
             k_t_entr_dev(i,k) =  k_t_entr_dev(i,k) + k_t_troen(k)
             k_m_entr_dev(i,k) =  k_m_entr_dev(i,k) + k_m_troen(k)*prandtlsfc_const
          end do

       end if

    end if 


    !NEGATIVELY BUOYANT PLUMES DRIVEN BY LW RADIATIVE COOLING AND/OR BUOYANCY REVERSAL
    !This part is done only if a level kcldtop can be found below zcldtopmax

    kmax = ibot+1
    do k = 1, ibot
       if ( zhalf_dev(i,k) < zcldtopmax) then
          kmax = k
          exit
       end if
    end do

    !Find cloud top and bottom using GRID BOX MEAN or IN-CLOUD value of qc. Decision occurs where qc is calculated
    kcldtop  = ibot+1
    do k = ibot,kmax,-1
       qlm1 = qc(k-1)  ! qc one level UP
       ql00 = qc(k)
       stab = slv(k-1) - slv(k) 
       if ( ( ql00  .ge. qlcrit ) .and. ( qlm1 .lt. qlcrit) .and. (stab .gt. 0.) ) then
           kcldtop  = k   
           exit
       end if
    end do

    !if (kcldtop .ge. ibot+1) then 
    if (kcldtop .lt. ibot+1) then 
       !Only do this is above is true, else skip to the end

       kcldtop2=min( kcldtop+1,nlev)
       ! Look one level further down in case first guess is a thin diffusive veil
       if ( (qc(kcldtop) .lt. 10*qlcrit ) .and. (qc(kcldtop2) .ge. 10*qc(kcldtop) ) ) then
          kcldtop=kcldtop2
       endif

       kcldbot  = ibot+1
       do k = ibot,kcldtop,-1
          qlm1 = qc(k-1)  ! qc one level UP
          ql00 = qc(k)
          if ( ( ql00  .lt. qlcrit ) .and. ( qlm1 .ge. qlcrit) ) then
             kcldbot  = k   
             exit
          end if
       end do

       if (radlw_dep .eq. 1) then

          !if (kcldtop .eq. kcldbot) then 
          if (kcldtop .ne. kcldbot) then 
             !Only do this is above is true, else skip to the end

             !With diffusion of ql, qi "cloud top" found via these quantities may be above radiation max
             kcldtop2=min( kcldtop+2,nlev)
             maxradf = maxval( -1.*tdtlw_in_dev(i,kcldtop:kcldtop2) )
             maxradf = maxradf*MAPL_CP*( (phalf_dev(i,kcldtop+1)-phalf_dev(i,kcldtop)) / MAPL_GRAV )
             maxradf = max( maxradf , 0. ) ! do not consider cloud tops that are heating

             !Calculate optimal mixing fraction (chis) for buoyancy reversal. Use effective heat of evap/subl. Ignore diffs across cldtop
             hlf = hleff(kcldtop)

             tmp1 = ( slv(kcldtop-1)  -  hlf*qc(kcldtop-1) ) - (slv(kcldtop) - hlf*qc(kcldtop) )
             tmp1 = dqs(kcldtop)*tmp1/MAPL_CP

             tmp2 = ( qv_dev(i,kcldtop-1)   +  qc(kcldtop-1) ) - ( qv_dev(i,kcldtop)  +  qc(kcldtop)   )  
   
             chis = -qc(kcldtop)*( 1 + hlf * dqs(kcldtop) / MAPL_CP )

             if ( ( tmp2 - tmp1 ) >= 0.0 ) then
                chis = 0.
             else
                chis = chis / ( tmp2 - tmp1 ) 
             endif

             if ( chis .gt. 1.0 ) then
                chis=1.0
             endif
            
             slmixture = ( 1.0-chis ) * ( slv(kcldtop) - hlf*qc(kcldtop)   )   +  chis  * ( slv(kcldtop-1)  -  hlf*qc(kcldtop-1) )

             !Compute temperature of parcel at cloud top, svpcp.
             svpcp = slmixture /MAPL_CP

             buoypert   = ( slmixture - slv(kcldtop) )/MAPL_CP

             !Calculate my best guess at the LCs' D parameter attributed to Siems et al.
             stab       = slv(kcldtop-1) - slv(kcldtop)
             if (stab .eq. 0.) then 
                dsiems     =  ( slv(kcldtop) - slmixture ) ! / 1.0  ! arbitrary, needs to be re-thought 
             else
                dsiems     =  ( slv(kcldtop) - slmixture ) / stab
             endif
             dsiems     =  min( dsiems, 10. )
             dsiems     =  max( dsiems,  0. )
             zradtop = zhalf_dev(i,kcldtop)

             !Find depth of radiatively driven convection 
             !Expose radperturb and other funny business outside of radml_depth
             radperturb   = min( maxradf/100. , 0.3 ) ! dim. argument based on 100m deep cloud over 1000s
             do_jump_exit = .true.
             critjump     = 0.3
             svpcp = svpcp - radperturb

             slvcp = slv/MAPL_CP

             if (kcldtop .lt. ibot) then 
                call radml_depth(i,ncol,nlev,kcldtop,ibot,svpcp,zradtop,critjump,do_jump_exit,&
                                 slvcp,zfull_dev,zhalf_dev,zradbase_dev,zradml_dev )      
             else
                zradbase_dev(i) = 0.0
                zradml_dev(i)   = zradtop  
             end if

             zcloud_dev(i) = zhalf_dev(i,kcldtop) - zhalf_dev(i,kcldbot)

             !if (zradml_dev(i) .le. 0.0 ) then 
             if (zradml_dev(i) .gt. 0.0 ) then 
                !Only do this is above is true, else skip to the end
   
                !Compute radiation driven scale
                vrad3 = MAPL_GRAV*zradml_dev(i)*maxradf/density(kcldtop)/slv(kcldtop)   

                !tmp1 here should be the buoyancy jump at cloud top SAK has it w/ resp to parcel property - svpcp. Im not sure about that.
                tmp1 = MAPL_GRAV*max(0.1,((slv(kcldtop-1)/MAPL_CP)-svpcp))/(slv(kcldtop)/MAPL_CP)

                !Straightforward buoyancy jump across cloud top
                tmp1 = MAPL_GRAV*max( 0.1, ( slv(kcldtop-1)-slv(kcldtop) )/MAPL_CP ) / ( slv(kcldtop) /MAPL_CP )

                !Compute buoy rev driven scale
                vbr3  = ( max( tmp1*zcloud_dev(i), 0.)**3 )
                vbr3  = Abuoy*(chis**2)*max(dsiems,0.)*SQRT( vbr3 ) 

                !Adjust velocity scales to prevent jumps near zradtop=zcldtopmax
                if ( zradtop .gt. zcldtopmax-500. ) then 
                   vrad3 = vrad3*(zcldtopmax - zradtop)/500.  
                   vbr3  = vbr3 *(zcldtopmax - zradtop)/500.  
                endif
                vrad3=max( vrad3, 0. ) ! these really should not be needed
                vbr3 =max( vbr3,  0. )

                vrad = vrad3 ** (1./3.)    
                vbrv = vbr3**(1./3.)
      
                tmp2 = (  vrad**2 + vbrv**2  ) / zradml_dev(i)
                wentr_rad = min(wentrmax,beta_rad_const*(vrad3+vbr3)/zradml_dev(i)/(tmp1+tmp2))
   
                wentr_brv =  beta_rad_const*vbr3/zradml_dev(i)/(tmp1+tmp2)

                !Fudgey adjustment of entrainment to reduce it for shallow boundary layers, increase for deep ones
                if ( zradtop .lt. 500. ) then
                   wentr_rad = 0.00
                endif
                if (( zradtop .gt. 500.) .and. (zradtop .le. 800. )) then
                   wentr_rad = wentr_rad * ( zradtop-500.) / 300.
                endif

                if ( zradtop .lt. 2400. ) then 
                   wentr_rad = wentr_rad * ( zradtop / 800. )
                else
                   wentr_rad = 3.*wentr_rad
                endif

                k_entr_tmp = min ( akmax, wentr_rad*(zfull_dev(i,kcldtop-1)-zfull_dev(i,kcldtop)) )

                radfq(kcldtop)          = 1.
                radpbl                  = 1.
                k_rad(kcldtop)          = k_entr_tmp
                k_t_entr_dev(i,kcldtop) = k_t_entr_dev(i,kcldtop) + k_entr_tmp
                k_m_entr_dev(i,kcldtop) = k_m_entr_dev(i,kcldtop) + k_entr_tmp
   
                !Handle case of radiatively driven top being the same top as surface driven top
                if (ipbl .eq. kcldtop .and. ipbl .gt. 0) then

                   tmp2 = ((vbr3+vrad3+vsurf3+vshear3)**(2./3.)) / zradml_dev(i)
      
                   wentr_rad = min( wentrmax,  max(0., ((beta_surf_const *(vsurf3 + vshear3) + &
                                             beta_rad_const*(vrad3+vbr3) )/ zradml_dev(i))/(tmp1+tmp2) ) )
                   wentr_pbl       = wentr_rad

                   k_entr_tmp = min ( akmax, wentr_rad*(zfull_dev(i,kcldtop-1)-zfull_dev(i,kcldtop)) )

                   pblfq(ipbl)             = 1.
                   radfq(kcldtop)          = 1.
                   radpbl                  = 1.
                   k_rad(kcldtop)          = k_entr_tmp
                   k_t_troen(ipbl)         = k_entr_tmp
                   k_m_troen(ipbl)         = k_entr_tmp
                   k_t_entr_dev(i,kcldtop) = k_entr_tmp
                   k_m_entr_dev(i,kcldtop) = k_entr_tmp
   
                end if

                !If there are any interior layers to calculate diffusivity
                if ( kcldtop .lt. ibot ) then   
   
                   do k = kcldtop+1,ibot
   
                      ztmp = max(0.,(zhalf_dev(i,k)-zradbase_dev(i))/zradml_dev(i) )
   
                      if (ztmp.gt.0.) then
   
                         radfq(k) = 1.
                         k_entr_tmp = khradfac_const*MAPL_KARMAN*( vrad+vbrv )*ztmp*zradml_dev(i)*ztmp*((1.-ztmp)**0.5)
                         k_entr_tmp = min ( k_entr_tmp, akmax )
                         k_rad(k) = k_entr_tmp
                         k_t_entr_dev (i,k) = k_t_entr_dev(i,k) + k_entr_tmp
                         k_m_entr_dev (i,k) = k_m_entr_dev(i,k) + k_entr_tmp*prandtlrad_const
   
                      end if
                   enddo
      
                end if

                !Handle special case of zradbase_dev < zsml_dev, in this case there should be no entrainment from the surface.
                if (zradbase_dev(i) .lt. zsml_dev(i) .and. convpbl .eq. 1. .and. ipbl .gt. kcldtop) then
                   wentr_pbl            = 0.
                   pblfq(ipbl)          = 0.
                   k_t_entr_dev(i,ipbl) = k_t_entr_dev (i,ipbl) - k_t_troen(ipbl)
                   k_m_entr_dev(i,ipbl) = k_m_entr_dev (i,ipbl) - k_m_troen(ipbl)          
                   k_t_troen(ipbl)      = 0.
                   k_m_troen(ipbl)      = 0. 
                end if

             endif

          endif !kcldtop .ne. kcldbot
 
       endif !radlw_dep    
 
    endif

    !Modify diffusivity coefficients using MAX( A , B )        
    do k = 2, ibot     
       diff_t_dev(i,k)   = max( k_t_entr_dev(i,k), diff_t_dev(i,k) )
       diff_m_dev(i,k)   = max( k_m_entr_dev(i,k), diff_m_dev(i,k) )
       k_t_entr_dev(i,k) = max( k_t_entr_dev(i,k),             0.0 )
       k_m_entr_dev(i,k) = max( k_m_entr_dev(i,k),             0.0 )
    enddo


 end do

end subroutine LOCK_DIFF


!======================================================================= 
subroutine mpbl_depth(i,ncol,nlev,tpfac, entrate, pceff, t, q, u, v, z, p, b_star, u_star , ipbl, ztop, ESTBLX)

 !Subroutine to calculate pbl depth

 IMPLICIT NONE

 !INPUTS
 integer, intent(in)                       :: i, nlev, ncol
 real,    intent(in)                       :: ESTBLX(:)
 real,    intent(in), dimension(ncol,nlev) :: t, z, q, p, u, v
 real,    intent(in), dimension(ncol)      :: b_star, u_star
 real,    intent(in)                       :: tpfac, entrate, pceff

 !OUTPUTS
 integer, intent(out)                       :: ipbl
 real,    intent(out),dimension(ncol)       :: ztop

 !LOCALS
 real :: tep,z1,z2,t1,t2,qp,pp,qsp,dqp,dqsp,u1,v1,u2,v2,du
 real :: ws,k_t_ref,entfr, tpfac_x, entrate_x,vscale,ptmp
 integer :: k

 !Calculate surface parcel properties
 tep  = t(i,nlev)
 qp   = q(i,nlev)

 !Wind dependence of plume character. 
 vscale    = 0.25 / 100. ! change of .25 m s-1 in 100 m

 !tpfac scales up bstar by inv. ratio of heat-bubble area to stagnant area
 tep  = tep * (1.+ tpfac * b_star(i)/MAPL_GRAV)

 !Search for level where this is exceeded              
 t1   = t(i,nlev)
 v1   = v(i,nlev)
 u1   = u(i,nlev)
 z1   = z(i,nlev)
 ztop(i) = z1

 do k = nlev-1 , 2, -1
    z2 = z(i,k)
    t2 = t(i,k)
    u2 = u(i,k)
    v2 = v(i,k)
    pp = p(i,k)

    du = sqrt ( ( u2 - u1 )**2 + ( v2 - v1 )**2 ) / (z2-z1)
    du = min(du,1.0e-8)

    entrate_x = entrate * ( 1.0 + du / vscale )

    entfr= min( entrate_x*(z2-z1), 0.99 )

    qp  = qp  + entfr*(q(i,k)-qp)

    !Dry adiabatic ascent through one layer. Static energy conserved. 
    tep = tep - MAPL_GRAV*( z2-z1 )/MAPL_CP

    !Environmental air entrained
    tep = tep + entfr*(t(i,k)-tep)

    !dqsp = dqsat(tep , pp , qsat=qsp,  pascals=.true. )
    ptmp = pp*0.01
    call DQSAT_sub_sca(dqsp,qsp,tep,ptmp,ESTBLX)

    dqp = max( qp - qsp, 0. )/(1.+(MAPL_ALHL/MAPL_CP)*dqsp )
    qp  = qp - dqp
    !"Precipitation efficiency" basically means fraction of condensation heating that gets applied to parcel
    tep = tep  + pceff * MAPL_ALHL * dqp/MAPL_CP  

    !If parcel temperature (tep) colder than env (t2) OR if entrainment too big, declare this the PBL top
    if ( (t2 .ge. tep) .or. ( entfr .ge. 0.9899 ) ) then
       ztop(i) = 0.5*(z2+z1)
       ipbl = k+1
       exit
    end if

    z1 = z2
    t1 = t2
    u1 = u2
    v1 = v2

 enddo

end subroutine mpbl_depth

subroutine radml_depth(i, ncol, nlev, toplev, botlev, svp, zt, critjump, do_jump_exit, t, zf, zh, zb, zml)

 !Subroutine to calculate bottom and depth of radiatively driven mixed layer

 IMPLICIT NONE

 !INPUTS
 integer, intent(in)                         :: i, toplev, botlev, ncol, nlev
 real,    intent(in)                         :: svp, zt, critjump
 real,    intent(in), dimension(nlev)        :: t
 real,    intent(in), dimension(ncol,nlev)   :: zf
 real,    intent(in), dimension(ncol,nlev+1) :: zh
 logical, intent(in)                         :: do_jump_exit

 !OUTPUTS
 real,    intent(out), dimension(ncol)        :: zb, zml

 !LOCALS
 real    :: svpar,h1,h2,t1,t2,entrate,entfr
 integer :: k

 !Initialize zml
 zml(i) = 0.

 !Calculate cloud top parcel properties
 svpar  = svp
 h1   = zf(i,toplev)
 t1   = t(toplev)
 entrate = 0.2/200.

 !Search for level where parcel is warmer than env. First cut out if parcel is already warmer than cloudtop. 
 if (t1.lt.svpar) then
    zb(i)  = h1
    zml(i) = 0.00
    return !udh - get rid of this if linearising
 endif

 !If above is false keep looking
 do k = 1,botlev
    if (k > toplev) then
       h2 = zf(i,k)
       t2 = t(k)

       if (t2.lt.svpar) then
          if ( abs(t1 - t2 ) .gt. 0.2 ) then 
             zb(i) = h2 + (h1 - h2)*(svpar - t2)/(t1 - t2)
             zb(i) = MAX( zb(i) , 0. )  ! final protection against interp problems
          else
             zb(i) = h2
          endif
          zml(i) = zt - zb(i)
          return !dh get rid of this if linearising
       end if

       if (do_jump_exit .and. (t1-t2) .gt. critjump .and. k .gt. toplev+1) then
          zb(i) = zh(i,k)
          zml(i) = zt - zb(i)
          return !dh get rid of this if linearising
       end if

       entfr = min( entrate*(h1-h2), 1.0 )
       svpar = svpar + entfr*(t2-svpar)

       h1 = h2
       t1 = t2

    end if

 enddo

 zb(i) = 0.
 zml(i) = zt

end subroutine radml_depth

subroutine diffusivity_pbl2(i, ncol, lm, h, kfac, k_ent, vsurf, frland, zm, k_m, k_t)

 !Subroutine to return the vertical K-profile of diffusion coefficients for the surface driven convective mixed layer.

 IMPLICIT NONE

 !INPUTS
 integer, intent(in)                         :: i, lm, ncol
 real,    intent(in), dimension(ncol)        :: h, frland
 real,    intent(in)                         :: k_ent, vsurf, kfac 
 real,    intent(in), dimension(ncol,1:lm+1) :: zm

 !OUTPUTS
 real,    intent(out), dimension(lm)         :: k_m, k_t

 !LOCALS
 real :: EE, hin, kfacx 
 integer :: k, kk

 !Kluge. Raise KHs over land
 if ( frland(i) < 0.5 ) then 
    kfacx = kfac 
 else
    kfacx = kfac*2.0
 endif

 hin    = 0.0 ! 200.  ! "Organization" scale for plume (m).

 do k = 1, lm
    k_m(k) = 0.0
    k_t(k) = 0.0
 end do

 if ( vsurf*h(i) .gt. 0. ) then
    EE  = 1.0 - sqrt( k_ent / ( kfacx * MAPL_KARMAN * vsurf * h(i) ) )
    EE  = max( EE , 0.7 )  ! If EE is too small, then punt, as LCs
    do k = 1, lm
       if ( ( zm(i,k) .le. h(i) ) .and.  ( zm(i,k) .gt. hin )  ) then
          k_t(k) = kfacx * MAPL_KARMAN * vsurf * ( zm(i,k)-hin ) * ( 1. - EE*( (zm(i,k)-hin)/(h(i)-hin) ))**2
       end if
    end do
 endif

 do k = 1, lm
    k_m(k) = k_t(k)
 end do

end subroutine diffusivity_pbl2

subroutine ESINIT(ESTBLX)
!Initialise look-up table

 IMPLICIT NONE

 !OUTPUT
 real(4), dimension(TABLESIZE) :: ESTBLX

 !LOCALS
 real(4), parameter :: ZEROC = 273.16, TMIX = -20.0

 real(4), dimension(TABLESIZE) :: ESTBLE, ESTBLW

 integer :: I
 real(4)    :: T, DELTA_T

 DELTA_T = 1.0/DEGSUBS

 do I=1,TABLESIZE

    T = (I-1)*DELTA_T + TMINTBL

    if(T>ZEROC) then
       call QSATLQU0(ESTBLE(I),T)
    else
       call QSATICE0(ESTBLE(I),T)
    end if

    call QSATLQU0(ESTBLW(I),T)

    T = T-ZEROC
    if(T>=TMIX .and. T<0.0) then
       ESTBLX(I) = ( T/TMIX )*( ESTBLE(I) - ESTBLW(I) ) + ESTBLW(I)
    else
       ESTBLX(I) = ESTBLE(I)
    end if

 end do

 end subroutine ESINIT


subroutine QSATLQU0(QS,TL)
!SUPERSATURATED AS LIQUID

 IMPLICIT NONE

 !INPUTS
 real(4) :: TL!, TMAXTBL

 !OUTPUTS
 real(4) :: QS

 !LOCALS
 real(4), parameter :: ZEROC   = 273.16
 real(4), parameter :: TMINLQU = ZEROC - 40.0

 real*8,  parameter :: B6 = 6.136820929E-11*100.0
 real*8,  parameter :: B5 = 2.034080948E-8 *100.0
 real*8,  parameter :: B4 = 3.031240396E-6 *100.0
 real*8,  parameter :: B3 = 2.650648471E-4 *100.0
 real*8,  parameter :: B2 = 1.428945805E-2 *100.0
 real*8,  parameter :: B1 = 4.436518521E-1 *100.0
 real*8,  parameter :: B0 = 6.107799961E+0 *100.0

 real(8) :: TX, EX, TI, TT

 TX = TL

 if    (TX<TMINLQU) then
    TI = TMINLQU
 elseif(TX>TMAXTBL) then
    TI = TMAXTBL
 else
    TI = TX
 end if

 TT = TI-ZEROC  !Starr polynomial fit
 EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)

 TL = TX
 QS = EX

 return

end subroutine QSATLQU0


subroutine QSATICE0(QS,TL)
!SUPERSATURATED AS ICE

 IMPLICIT NONE   
      
 !INPUTS
 real(4) :: TL

 !OUTPUTS
 real(4) :: QS

 !LOCALS
 real(4), parameter :: ZEROC = 273.16, TMINSTR = -95.0
 real(4), parameter :: TMINICE = ZEROC + TMINSTR

 real(4), parameter :: TSTARR1 = -75.0, TSTARR2 = -65.0, TSTARR3 = -50.0,  TSTARR4 = -40.0

 real*8,  parameter :: BI6= 1.838826904E-10*100.0
 real*8,  parameter :: BI5= 4.838803174E-8 *100.0
 real*8,  parameter :: BI4= 5.824720280E-6 *100.0
 real*8,  parameter :: BI3= 4.176223716E-4 *100.0
 real*8,  parameter :: BI2= 1.886013408E-2 *100.0
 real*8,  parameter :: BI1= 5.034698970E-1 *100.0
 real*8,  parameter :: BI0= 6.109177956E+0 *100.0
 real*8,  parameter :: S16= 0.516000335E-11*100.0
 real*8,  parameter :: S15= 0.276961083E-8 *100.0
 real*8,  parameter :: S14= 0.623439266E-6 *100.0
 real*8,  parameter :: S13= 0.754129933E-4 *100.0
 real*8,  parameter :: S12= 0.517609116E-2 *100.0
 real*8,  parameter :: S11= 0.191372282E+0 *100.0
 real*8,  parameter :: S10= 0.298152339E+1 *100.0
 real*8,  parameter :: S26= 0.314296723E-10*100.0
 real*8,  parameter :: S25= 0.132243858E-7 *100.0
 real*8,  parameter :: S24= 0.236279781E-5 *100.0
 real*8,  parameter :: S23= 0.230325039E-3 *100.0
 real*8,  parameter :: S22= 0.129690326E-1 *100.0
 real*8,  parameter :: S21= 0.401390832E+0 *100.0
 real*8,  parameter :: S20= 0.535098336E+1 *100.0

 real(4) :: TX, TI, TT, W, EX

 TX = TL

 if (TX<TMINICE) then
    TI = TMINICE
 elseif(TX>ZEROC  ) then
    TI = ZEROC
 else
    TI = TX
 end if

 TT = TI - ZEROC
 if (TT < TSTARR1 ) then
     EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
 elseif(TT >= TSTARR1 .and. TT < TSTARR2) then
     W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
     EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
              + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
 elseif(TT >= TSTARR2 .and. TT < TSTARR3) then
     EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
 elseif(TT >= TSTARR3 .and. TT < TSTARR4) then
     W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
     EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
              + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
 else
     EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
 endif

 QS = EX

 return
 
end subroutine QSATICE0

subroutine DQSAT_sub_sca(DQSi,QSSi,TEMP,PLO,ESTBLX)
!COMPUTES SATURATION VAPOUR PRESSURE QSSi AND GRADIENT w.r.t TEMPERATURE DQSi.
!INPUTS ARE TEMPERATURE AND PLO (PRESSURE AT T-LEVELS)
!VALES ARE COMPUTED FROM LOOK-UP TALBE (PIECEWISE LINEAR)

 IMPLICIT NONE

 !Inputs
 real(4) :: TEMP, PLO
 real(4) :: ESTBLX(:)

 !Outputs
 real(4) :: DQSi, QSSi

 !Locals
 real(4), parameter :: MAX_MIXING_RATIO = 1.0
 real(4),    parameter :: ESFAC = MAPL_H2OMW/MAPL_AIRMW

 real(4) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 integer :: IT

 TL = TEMP
 PL = PLO

 PP = PL*100.0
 
 if (TL<=TMINTBL) then
    TI = TMINTBL
 elseif(TL>=TMAXTBL-.001) then
    TI = TMAXTBL-.001
 else
    TI = TL
 end if

 TT = (TI - TMINTBL)*DEGSUBS+1
 IT = int(TT)

 DQQ =  ESTBLX(IT+1) - ESTBLX(IT)
 QQ  =  (TT-IT)*DQQ + ESTBLX(IT)

 if (PP <= QQ) then
    QSAT = MAX_MIXING_RATIO
    DQSAT = 0.0
 else
    DD = 1.0/(PP - (1.0-ESFAC)*QQ)
    QSAT = ESFAC*QQ*DD
    DQSAT = (ESFAC*DEGSUBS)*DQQ*PP*(DD*DD)
 end if

 DQSi = DQSAT
 QSSi = QSAT

end subroutine DQSAT_sub_sca



END MODULE BLDRIVER
