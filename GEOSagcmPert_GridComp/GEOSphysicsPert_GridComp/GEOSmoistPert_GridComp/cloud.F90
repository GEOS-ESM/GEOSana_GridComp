MODULE CLOUD

IMPLICIT NONE

PRIVATE

!Make subroutines public for adjoint checkpointing
PUBLIC CLOUD_DRIVER, CLOUD_TIDY, MELTFREEZE, GET_ICE_FRACTION, CONVEC_SRC, PDF_WIDTH, LS_CLOUD, &
       AUTOCONVERSION_LS, AUTOCONVERSION_CNV, ICE_SETTLEFALL_LS, ICE_SETTLEFALL_CNV, PRECIPANDEVAP, &
       PDFFRAC, PDFCONDENSATE, CONS_SUNDQ3, CONS_ALHX, CONS_MICROPHYS, MARSHPALM, EVAP_CNV, SUBL_CNV, LDRADIUS, &
       DQSAT_BAC, DQSATs_BAC

CONTAINS

subroutine cloud_driver ( DT, IM, JM, LM, th, q, ple,                  &
                          CNV_DQLDT, CNV_MFD, CNV_PRC3, CNV_UPDF,      &
                          QI_ls, QL_ls, QI_con, QL_con, CF_LS, CF_con, &
                          FRLAND, PHYSPARAMS, ESTBLX, KHu, KHl,        &
                          CONS_RUNIV, CONS_KAPPA,  CONS_AIRMW,         &
                          CONS_H2OMW, CONS_GRAV,   CONS_ALHL,          &
                          CONS_ALHF,  CONS_PI,     CONS_RGAS,          &
                          CONS_CP,    CONS_VIREPS, CONS_ALHS,          &
                          CONS_TICE,  CONS_RVAP,   CONS_P00, do_moist_physics            )

!NAME LIST
!DT              - Physics timestep
!IM              - Number of points in x-dir (per processor)
!JM              - Number of points in y-dir 
!LM              - Number of vertical levels
!th           (d)- Potential temperature (K), defined with p0 = 1000hPa (th=T at the surface)
!q            (d)- Specific humidity (kg kg^-1)
!ple             - Pressure at level edge (Pa)
!CNV_DQLDT    (d)- Convective condensate source from RAS (kg m^-2 s^-1)
!CNV_MFD      (d)- Convective detraining mass flux from RAS (kg m^-2 s^-1)
!CNV_PRC3     (d)- Convective precipitation from RAS (kg m^-2 s^-1)
!CNV_UPDF     (d)- Convective updraft areal fraction from RAS (1)
!QI_ls        (d)- Mass fraction of large scale cloud ice water (kg kg^-1)
!QL_ls        (d)- Mass fraction of large scale cloud liquid water (kg kg^-1)
!QI_con       (d)- Mass fraction of convective cloud ice water (kg kg^-1)
!QL_con       (d)- Mass fraction of convective cloud liquid water (kg kg^-1)
!CF_LS        (d)- Large scale cloud area fraction (1)
!CF_con       (d)- Convective cloud area fraction (1)
!FRLAND          - Fraction of land, 0 to 1
!PHYSPARAMS      - Set of constants
!ESTBLX          - Table of values for calculating qsat and dqsatdT

!(d) = TLM Dependent variable

IMPLICIT NONE

!INPUTS
integer, intent(in) :: IM, JM, LM, do_moist_physics
real(8), intent(in) :: DT, FRLAND(IM,JM), PHYSPARAMS(:)
real(8), intent(in), dimension(IM,JM,LM) :: CNV_DQLDT, CNV_MFD, CNV_UPDF, CNV_PRC3

real(8), intent(in) :: ESTBLX(:)

real(8), intent(in), dimension(IM,JM,0:LM) :: ple
integer, intent(in), dimension(IM,JM) :: KHu, KHl

!MAPL_CONSTANTS REDEFINED FOR USE IN AUTODIFF TOOL
real(8), intent(in) :: CONS_RUNIV, CONS_KAPPA,  CONS_AIRMW
real(8), intent(in) :: CONS_H2OMW, CONS_GRAV,   CONS_ALHL
real(8), intent(in) :: CONS_ALHF,  CONS_PI,     CONS_RGAS
real(8), intent(in) :: CONS_CP,    CONS_VIREPS, CONS_ALHS
real(8), intent(in) :: CONS_TICE,  CONS_RVAP,   CONS_P00 

!PROGNOSTICS
real(8), intent(inout), dimension(IM,JM,LM) :: th, q
real(8), intent(inout), dimension(IM,JM,LM) :: QI_ls, QL_ls, QI_con, QL_con, CF_con, CF_ls

!OUTPUTS (DIAGNOSTICS)

!LOCALS
integer :: i, j, k, l, ktop

real(8), dimension (IM,JM,LM) :: ph, pih, MASS, iMASS, T, DZET, QDDF3, RH, DP, DM
real(8), dimension (IM,JM,LM+1) :: ZET
real(8), dimension (IM,JM,0:LM) :: p, pi
real(8), dimension (IM,JM,LM) :: QS, DQSDT, DQS

real(8), dimension(IM,JM) :: VMIP

real(8) :: CF_tot !Precip amounts and fall rate

real(8) :: ALPHA, ALHX3, RHCRIT
real(8) :: AA, BB !Microphyiscal constants

real(8), parameter :: PI_0  = 4.*atan(1.)
real(8), parameter :: RHO_W = 1.0e3  !Density of liquid water in kg/m^3
real(8)            :: T_ICE_MAX

!PHYSPARAMS constants
real(8) :: CNV_BETA,ANV_BETA,LS_BETA,RH00,C_00,LWCRIT,C_ACC,C_EV_R,C_EV_S,CLDVOL2FRC
real(8) :: RHSUP_ICE,SHR_EVAP_FAC,MIN_CLD_WATER,CLD_EVP_EFF,LS_SDQV2,LS_SDQV3,LS_SDQVT1
real(8) :: ANV_SDQV2,ANV_SDQV3,ANV_SDQVT1,ANV_TO_LS,N_WARM,N_ICE,N_ANVIL,N_PBL
real(8) :: ANV_ICEFALL_C,LS_ICEFALL_C,REVAP_OFF_P,CNVENVFC,WRHODEP,T_ICE_ALL,CNVICEPARAM
real(8) :: CNVDDRFC,ANVDDRFC,LSDDRFC,minrhcrit,maxrhcrit,turnrhcrit,maxrhcritland
real(8) :: MIN_RL,MIN_RI,MAX_RL,MAX_RI,RI_ANV
integer :: NSMAX,DISABLE_RAD,ICEFRPWR,tanhrhcrit,fr_ls_wat,fr_ls_ice,fr_an_wat,fr_an_ice,pdfflag

real(8) :: LSENVFC, ANVENVFC
real(8) :: QRN_cu, QSN_cu, QRN_an, QSN_an, QRN_ls, QSN_ls, QRN_cu_1D
real(8) :: QT_tmpi_1, QT_tmpi_2, QLT_tmp, QIT_tmp

real(8) :: PRN_above_cu_new, PRN_above_an_new, PRN_above_ls_new
real(8) :: PRN_above_cu_old, PRN_above_an_old, PRN_above_ls_old
real(8) :: PSN_above_cu_new, PSN_above_an_new, PSN_above_ls_new
real(8) :: PSN_above_cu_old, PSN_above_an_old, PSN_above_ls_old
real(8) :: EVAP_DD_CU_above_new, EVAP_DD_AN_above_new, EVAP_DD_LS_above_new
real(8) :: EVAP_DD_CU_above_old, EVAP_DD_AN_above_old, EVAP_DD_LS_above_old
real(8) :: SUBL_DD_CU_above_new, SUBL_DD_AN_above_new, SUBL_DD_LS_above_new
real(8) :: SUBL_DD_CU_above_old, SUBL_DD_AN_above_old, SUBL_DD_LS_above_old

real(8) :: AREA_LS_PRC1, AREA_UPD_PRC1, AREA_ANV_PRC1
real(8) :: TOT_PREC_UPD, TOT_PREC_ANV, TOT_PREC_LS, AREA_UPD_PRC, AREA_ANV_PRC, AREA_LS_PRC
real(8) :: QTMP2

real(8) :: RHEXCESS, TPW, NEGTPW

  !LS_CLOUD FILTERING
  integer :: ii, cloud_pertmod
  real(8) :: xx(8) 
  real(8) :: ttraj, qtraj, qi_lstraj, qi_contraj, ql_lstraj, ql_contraj, cf_lstraj, cf_contraj, phtraj 
  real(8) :: tpert, qpert, qi_lspert, qi_conpert, ql_lspert, ql_conpert, cf_lspert, cf_conpert
  real(8) :: Jacobian(8,8), A(8,8)



  !Total Filtering
  real(8) :: TOTfilt_T, TOTfilt_ql, TOTfilt_qi
  real(8) :: t_p_preall, ql_ls_p_preall, ql_con_p_preall, qi_ls_p_preall, qi_con_p_preall

  !Sink Filtering
  real(8) :: SINKfilt_ql, SINKfilt_qi, SINKfilt_CF
  real(8) :: t_p_presink, q_p_presink
  real(8) ::  ql_ls_p_presink, ql_con_p_presink
  real(8) ::  qi_ls_p_presink, qi_con_p_presink
  real(8) ::  cf_con_p_presink


!Highest level of calculations
KTOP = 30

!Get Constants from CLOUDPARAMS
CNV_BETA      = PHYSPARAMS(1)  ! Area factor for convective rain showers (non-dim)
ANV_BETA      = PHYSPARAMS(2)  ! Area factor for anvil rain showers (non-dim)
LS_BETA       = PHYSPARAMS(3)  ! Area factor for Large Scale rain showers (non-dim)
RH00          = PHYSPARAMS(4)  ! Critical relative humidity
C_00          = PHYSPARAMS(5)
LWCRIT        = PHYSPARAMS(6)
C_ACC         = PHYSPARAMS(7)
C_EV_R        = PHYSPARAMS(8)
C_EV_S        = PHYSPARAMS(56)
CLDVOL2FRC    = PHYSPARAMS(9)
RHSUP_ICE     = PHYSPARAMS(10)
SHR_EVAP_FAC  = PHYSPARAMS(11)
MIN_CLD_WATER = PHYSPARAMS(12)
CLD_EVP_EFF   = PHYSPARAMS(13)
NSMAX         = INT( PHYSPARAMS(14)  )
LS_SDQV2      = PHYSPARAMS(15)
LS_SDQV3      = PHYSPARAMS(16)
LS_SDQVT1     = PHYSPARAMS(17)
ANV_SDQV2     = PHYSPARAMS(18)
ANV_SDQV3     = PHYSPARAMS(19)
ANV_SDQVT1    = PHYSPARAMS(20)
ANV_TO_LS     = PHYSPARAMS(21)
N_WARM        = PHYSPARAMS(22)
N_ICE         = PHYSPARAMS(23)
N_ANVIL       = PHYSPARAMS(24)
N_PBL         = PHYSPARAMS(25)
DISABLE_RAD   = INT( PHYSPARAMS(26) )
ANV_ICEFALL_C = PHYSPARAMS(28)
LS_ICEFALL_C  = PHYSPARAMS(29)
REVAP_OFF_P   = PHYSPARAMS(30)
CNVENVFC      = PHYSPARAMS(31)
WRHODEP       = PHYSPARAMS(32)
T_ICE_ALL     = PHYSPARAMS(33) + CONS_TICE
CNVICEPARAM   = PHYSPARAMS(34)
ICEFRPWR      = INT( PHYSPARAMS(35) + .001 )
CNVDDRFC      = PHYSPARAMS(36)
ANVDDRFC      = PHYSPARAMS(37)
LSDDRFC       = PHYSPARAMS(38)
tanhrhcrit    = INT( PHYSPARAMS(41) )
minrhcrit     = PHYSPARAMS(42)
maxrhcrit     = PHYSPARAMS(43)
turnrhcrit    = PHYSPARAMS(45)
maxrhcritland = PHYSPARAMS(46)
fr_ls_wat     = INT( PHYSPARAMS(47) )
fr_ls_ice     = INT( PHYSPARAMS(48) )
fr_an_wat     = INT( PHYSPARAMS(49) )
fr_an_ice     = INT( PHYSPARAMS(50) )
MIN_RL        = PHYSPARAMS(51)
MIN_RI        = PHYSPARAMS(52)
MAX_RL        = PHYSPARAMS(53)
MAX_RI        = PHYSPARAMS(54)
RI_ANV        = PHYSPARAMS(55)
pdfflag       = INT(PHYSPARAMS(57))

T_ICE_MAX = CONS_TICE

!Initialize the saving of downdraft values.
PRN_above_cu_new = 0.
PRN_above_an_new = 0.
PRN_above_ls_new = 0.
PSN_above_cu_new = 0.
PSN_above_an_new = 0.
PSN_above_ls_new = 0.
EVAP_DD_CU_above_new = 0.
EVAP_DD_AN_above_new = 0.
EVAP_DD_LS_above_new = 0.
SUBL_DD_CU_above_new = 0.
SUBL_DD_AN_above_new = 0.
SUBL_DD_LS_above_new = 0.

!Convert to hPa and average pressure to temperature levels
p  = ple*0.01
ph = 0.5*(p(:,:,0:LM-1) +  p(:,:,1:LM  ) ) 

!Calculate Exner Pressure at temperature levels
pi  = (p /1000.)**(CONS_RGAS/CONS_CP)
pih = (ph/1000.)**(CONS_RGAS/CONS_CP)

!Calculate temperature
T = th*pih

!Compute QS and DQSDT
 call DQSAT_BAC(DQSDT,QS,T,Ph,im,jm,lm,ESTBLX,CONS_H2OMW,CONS_AIRMW)

!Relative humidity
RH = Q/QS

!Compute layer mass and 1/mass
MASS =  ( p(:,:,1:LM) - p(:,:,0:LM-1) )*100./CONS_GRAV
iMASS = 1 / MASS

!Level thickness
DZET(:,:,1:LM) = th(:,:,1:LM) * (pi(:,:,1:LM) - pi(:,:,0:LM-1)) * CONS_CP/CONS_GRAV

!Level heights
ZET(:,:,LM+1) = 0.0
DO K = LM, 1, -1
   ZET(:,:,K) = ZET(:,:,K+1)+DZET(:,:,K)
END DO

WHERE ( ZET(:,:,1:LM) < 3000. )
   QDDF3 = -( ZET(:,:,1:LM)-3000. ) * ZET(:,:,1:LM) * MASS
ELSEWHERE
   QDDF3 = 0.
END WHERE   

DO I = 1,IM
   DO J = 1,JM
      VMIP(I,J) = SUM( QDDF3(I,J,:) )
   END DO
END DO
DO K = 1,LM
   QDDF3(:,:,K) = QDDF3(:,:,K) / VMIP
END DO

!Pressure and mass thickness for use in cleanup.
DP = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
DM = DP*(1./CONS_GRAV)

!Begin loop over all grid boxes.
DO i = 1,IM
   DO j = 1,JM
      DO k = KTOP,LM

        !Save the inputs to the scheme for filtering
        t_p_preall = t(i, j, k)
        ql_ls_p_preall = ql_ls(i, j, k)
        ql_con_p_preall = ql_con(i, j, k)
        qi_ls_p_preall = qi_ls(i, j, k)
        qi_con_p_preall = qi_con(i, j, k)

         if (K == KTOP) then
            TOT_PREC_UPD = 0.
            TOT_PREC_ANV = 0.
            TOT_PREC_LS  = 0.

            AREA_UPD_PRC = 0.
            AREA_ANV_PRC = 0.
            AREA_LS_PRC  = 0.
         end if

         !Initialize precips, except QRN_CU which comes from RAS 
         QRN_LS    = 0.
         QRN_AN    = 0.
         QRN_CU_1D = 0.
         QSN_LS    = 0.
         QSN_AN    = 0.
         QSN_CU    = 0.

         QRN_CU_1D = CNV_PRC3(i,j,k) !Ras Rain         

         !Tidy up where fractions or cloud is too low
         call cloud_tidy( q(i,j,k),      &
                          T(i,j,k),      &
                          QL_ls(i,j,k),  &
                          QI_ls(i,j,k),  &
                          CF_ls(i,j,k),  &
                          QL_con(i,j,k), &
                          QI_con(i,j,k), &
                          CF_con(i,j,k), &
                          CONS_ALHL,     &
                          CONS_ALHS,     &
                          CONS_CP        )

         !Phase changes for large scale cloud.
         call meltfreeze ( DT,             &
                           T(i,j,k),       &
                           QL_ls(i,j,k),   & 
                           QI_ls(i,j,k),   &
                           T_ICE_ALL,      &
                           T_ICE_MAX,      &
                           ICEFRPWR,       &
                           CONS_ALHL,      &
                           CONS_ALHS,      &
                           CONS_CP         )

         !Phase changes for convective cloud.
         call meltfreeze ( DT,             &
                           T(i,j,k),       &
                           QL_con(i,j,k),  & 
                           QI_con(i,j,k),  &
                           T_ICE_ALL,      &
                           T_ICE_MAX,      &
                           ICEFRPWR,       &
                           CONS_ALHL,      &
                           CONS_ALHS,      &
                           CONS_CP         )




         !STAGE 1 - Compute convective clouds from RAS diagnostics
         call convec_src( DT,               &
                          MASS(i,j,k),      &
                          iMASS(i,j,k),     &
                          T(i,j,k),         &
                          q(i,j,k),         &
                          CNV_DQLDT(i,j,k), &
                          CNV_MFD(i,j,k),   &
                          QL_con(i,j,k),    &
                          QI_con(i,j,k),    &
                          CF_con(i,j,k),    &
                          QS(i,j,k),        &
                          CONS_ALHS,        &
                          CONS_ALHL,        &
                          CONS_CP,          &
                          T_ICE_ALL,        &
                          T_ICE_MAX,        &
                          ICEFRPWR          )
       

         !STAGE 2a - Get PDF attributes
         call pdf_width ( ph(i,j,k),        &
                          FRLAND(i,j),      &
                          maxrhcrit,        &
                          maxrhcritland,    &
                          turnrhcrit,       &
                          minrhcrit,        &
                          pi_0,             &
                          ALPHA             )

         ALPHA  = MAX(  ALPHA , 1.0 - RH00 )
         RHCRIT = 1.0 - ALPHA


         cloud_pertmod = 1

         !STAGE 2b - Use PDF to compute large scale cloud effects and diagnostics,
         !           also update convection clouds
         call ls_cloud( DT,                 &
                        ALPHA,              &
                        PDFFLAG,            &                        
                        ph(i,j,k),          &
                        T(i,j,k),           &
                        q(i,j,k),           &
                        QL_ls(i,j,k),       &
                        QL_con(i,j,k),      &
                        QI_ls(i,j,k),       &
                        QI_con(i,j,k),      &
                        CF_ls(i,j,k),       &
                        CF_con(i,j,k),      &
                        CONS_ALHL,          &
                        CONS_ALHF,          &
                        CONS_ALHS,          &
                        CONS_CP,            &
                        CONS_H2OMW,         &
                        CONS_AIRMW,         &
                        T_ICE_ALL,          &
                        T_ICE_MAX,          &
                        ICEFRPWR,           &
                        ESTBLX,             &
                        cloud_pertmod,      &
                        do_moist_physics    )

        !SAVE PRESINKS INPUTS
        t_p_presink = t(i,j,k)
        q_p_presink = q(i,j,k)
        qi_ls_p_presink = qi_ls(i,j,k)
        qi_con_p_presink = qi_con(i,j,k)
        ql_ls_p_presink = ql_ls(i,j,k)
        ql_con_p_presink = ql_con(i,j,k)
        cf_con_p_presink = cf_con(i,j,k)


         !Clean up where too much overall cloud.
         CF_tot = CF_ls(i,j,k) + CF_con(i,j,k)
         IF ( CF_tot > 1.00 ) THEN
            CF_ls(i,j,k)  = CF_ls (i,j,k)*(1.00 / CF_tot )
            CF_con(i,j,k) = CF_con(i,j,k)*(1.00 / CF_tot )
         END IF
         CF_tot = CF_ls(i,j,k) + CF_con(i,j,k)


         !STAGE 3 - Evap, Sublimation and Autoconversion

         !Evaporation and sublimation of anvil cloud
         call evap_cnv( DT,            &
                        RHCRIT,        &
                        ph(i,j,k),     &
                        T(i,j,k),      &
                        Q(i,j,k),      &
                        QL_con(i,j,k), &
                        QI_con(i,j,k), &
                        CF_con(i,j,k), &
                        CF_ls(i,j,k),  &
                        QS(i,j,k),     &  
                        RHO_W,         &
                        CLD_EVP_EFF,   &
                        CONS_H2OMW,    &
                        CONS_AIRMW,    &
                        CONS_ALHL,     &
                        CONS_RVAP,     &
                        CONS_RGAS,     &
                        CONS_PI,       &
                        CONS_CP        )

         call subl_cnv( DT,            & 
                        RHCRIT,        &
                        ph(i,j,k),     &
                        T(i,j,k),      &
                        Q(i,j,k),      &                           
                        QL_con(i,j,k), &
                        QI_con(i,j,k), &
                        CF_con(i,j,k), &
                        CF_ls(i,j,k),  &
                        QS(i,j,k),     &  
                        RHO_W,         &
                        CLD_EVP_EFF,   &
                        CONS_H2OMW,    &
                        CONS_AIRMW,    &
                        CONS_ALHL,     &
                        CONS_RVAP,     &
                        CONS_RGAS,     &
                        CONS_PI,       &
                        CONS_CP,       &
                        CONS_ALHS      )

         !Autoconversion
         call autoconversion_ls ( DT,            &
                                  QL_ls(i,j,k),  &
                                  QRN_LS,        & 
                                  T(i,j,k),      & 
                                  ph(i,j,k),     & 
                                  CF_ls(i,j,k),  &
                                  LS_SDQV2,      &
                                  LS_SDQV3,      &   
                                  LS_SDQVT1,     &
                                  C_00,          &
                                  LWCRIT,        &
                                  DZET(i,j,k)    )

         call autoconversion_cnv ( DT,            &
                                   QL_con(i,j,k), &
                                   QRN_AN,        & 
                                   T(i,j,k),      & 
                                   ph(i,j,k),     & 
                                   CF_con(i,j,k), &
                                   ANV_SDQV2,     &
                                   ANV_SDQV3,     &   
                                   ANV_SDQVT1,    &
                                   C_00,          &
                                   LWCRIT,        &
                                   DZET(i,j,k)    )


         !STAGE 4 - Fall and Re-evaporation of precip
         call ice_settlefall_cnv( WRHODEP,       &
                                  QI_con(i,j,k), &
                                  ph(i,j,k),     & 
                                  T(i,j,k),      & 
                                  CF_con(i,j,k), & 
                                  CONS_RGAS,     &
                                  KHu(i,j),      &
                                  KHl(i,j),      &
                                  k,             &
                                  DT,            &
                                  DZET(i,j,k),   &
                                  QSN_an,        &
                                  ANV_ICEFALL_C  )

         call ice_settlefall_ls( WRHODEP,       &
                                 QI_ls(i,j,k),  &
                                 ph(i,j,k),     & 
                                 T(i,j,k),      & 
                                 CF_ls(i,j,k),  & 
                                 CONS_RGAS,     &
                                 KHu(i,j),      &
                                 KHl(i,j),      &
                                 k,             &
                                 DT,            &
                                 DZET(i,j,k),   &
                                 QSN_ls,        &
                                 LS_ICEFALL_C   )

         !"Freeze" out any conv. precip, not done in RAS. This is
         ! precip w/ large particles, so freezing is strict.
         QTMP2 = 0.
         if ( T(i,j,k) < CONS_TICE ) then
            QTMP2     = QRN_CU_1D
            QSN_CU    = QRN_CU_1D
            QRN_CU_1D = 0.
            T(i,j,k)  = T(i,j,k) + QSN_CU*(CONS_ALHS-CONS_ALHL)/CONS_CP
         end if

         !Area
         AREA_LS_PRC1  = 0.0
         AREA_UPD_PRC1 = 0.0
         AREA_ANV_PRC1 = 0.0

         TOT_PREC_UPD  = TOT_PREC_UPD + ( ( QRN_CU_1D + QSN_CU ) * MASS(i,j,k) )
         AREA_UPD_PRC  = AREA_UPD_PRC + ( CNV_UPDF(i,j,k)* ( QRN_CU_1D + QSN_CU )* MASS(i,j,k) )

         TOT_PREC_ANV  = TOT_PREC_ANV + ( ( QRN_AN + QSN_AN ) * MASS(i,j,k) )
         AREA_ANV_PRC  = AREA_ANV_PRC + ( CF_con(i,j,k)* ( QRN_AN + QSN_AN )* MASS(i,j,k) )

         TOT_PREC_LS   = TOT_PREC_LS  + ( ( QRN_LS + QSN_LS ) * MASS(i,j,k) )
         AREA_LS_PRC   = AREA_LS_PRC  + ( CF_ls(i,j,k)*  ( QRN_LS + QSN_LS )* MASS(i,j,k) )

         if ( TOT_PREC_ANV > 0.0 ) AREA_ANV_PRC1 = MAX( AREA_ANV_PRC/TOT_PREC_ANV, 1.E-6 )
         if ( TOT_PREC_UPD > 0.0 ) AREA_UPD_PRC1 = MAX( AREA_UPD_PRC/TOT_PREC_UPD, 1.E-6 )
         if ( TOT_PREC_LS  > 0.0 ) AREA_LS_PRC1  = MAX( AREA_LS_PRC/TOT_PREC_LS,   1.E-6 )

         AREA_LS_PRC1  = LS_BETA  * AREA_LS_PRC1
         AREA_UPD_PRC1 = CNV_BETA * AREA_UPD_PRC1
         AREA_ANV_PRC1 = ANV_BETA * AREA_ANV_PRC1

         IF (K == LM) THEN ! We have accumulated over the whole column
            if ( TOT_PREC_ANV > 0.0 ) AREA_ANV_PRC = MAX( AREA_ANV_PRC/TOT_PREC_ANV, 1.E-6 )
            if ( TOT_PREC_UPD > 0.0 ) AREA_UPD_PRC = MAX( AREA_UPD_PRC/TOT_PREC_UPD, 1.E-6 )
            if ( TOT_PREC_LS  > 0.0 ) AREA_LS_PRC  = MAX( AREA_LS_PRC/TOT_PREC_LS,   1.E-6 )

            AREA_LS_PRC  = LS_BETA  * AREA_LS_PRC
            AREA_UPD_PRC = CNV_BETA * AREA_UPD_PRC
            AREA_ANV_PRC = ANV_BETA * AREA_ANV_PRC
         END IF

         !Get micro-physical constants
         call cons_alhx( T(i,j,k),          &
                         ALHX3,             &
                         T_ICE_MAX,         &
                         T_ICE_ALL,         &
                         CONS_ALHS,         &
                         CONS_ALHL          )

         call cons_microphys( T(i,j,k),          &
                              ph(i,j,k),         &
                              Qs(i,j,k),         &
                              AA,                &
                              BB,                &
                              CONS_H2OMW,        &
                              CONS_AIRMW,        &
                              CONS_RVAP,         &
                              ALHX3              )


         !Precip Scheme Expects Total Cloud Liquid
         QLT_tmp = QL_ls(i,j,k) + QL_con(i,j,k)
         QIT_tmp = QI_ls(i,j,k) + QI_con(i,j,k)

         PRN_above_CU_old = PRN_above_CU_new
         PSN_above_CU_old = PSN_above_CU_new
         EVAP_DD_CU_above_old = EVAP_DD_CU_above_new
         SUBL_DD_CU_above_old = SUBL_DD_CU_above_new

         !Precip and Evap for Convection
         call  precipandevap( K                         , &
                              KTOP                      , &
                              LM                        , &
                              DT                        , &
                              FRLAND(i,j)               , &
                              RHCRIT                    , &
                              QRN_CU_1D                 , &
                              QSN_CU                    , &
                              QLT_tmp                   , &
                              QIT_tmp                   , &
                              T(i,j,k)                  , &
                              Q(i,j,k)                  , &
                              mass(i,j,k)               , &
                              imass(i,j,k)              , &
                              ph(i,j,k)                 , &
                              DZET(i,j,k)               , &
                              QDDF3(i,j,k)              , &
                              AA                        , &
                              BB                        , &
                              AREA_UPD_PRC1             , &
                              PRN_above_CU_old          , &
                              PRN_above_CU_new          , &
                              PSN_above_CU_old          , &
                              PSN_above_CU_new          , &
                              EVAP_DD_CU_above_old      , &
                              EVAP_DD_CU_above_new      , &
                              SUBL_DD_CU_above_old      , &
                              SUBL_DD_CU_above_new      , &
                              CNVENVFC                  , &
                              CNVDDRFC                  , & 
                              CONS_ALHF                 , &
                              CONS_ALHS                 , &
                              CONS_ALHL                 , &
                              CONS_CP                   , &
                              CONS_TICE                 , &
                              CONS_H2OMW                , & 
                              CONS_AIRMW                , &
                              REVAP_OFF_P               , &
                              C_ACC                     , &
                              C_EV_R                    , &
                              C_EV_S                    , &
                              RHO_W                     , &
                              ESTBLX                      ) 

         PRN_above_AN_old = PRN_above_AN_new
         PSN_above_AN_old = PSN_above_AN_new
         EVAP_DD_AN_above_old = EVAP_DD_AN_above_new
         SUBL_DD_AN_above_old = SUBL_DD_AN_above_new

         !Precip and Evap for Anvil
         ANVENVFC = 1.0
         call  precipandevap( K                         , &
                              KTOP                      , &
                              LM                        , &
                              DT                        , &
                              FRLAND(i,j)               , &
                              RHCRIT                    , &
                              QRN_AN                    , &
                              QSN_AN                    , &
                              QLT_tmp                   , &
                              QIT_tmp                   , &
                              T(i,j,k)                  , &
                              Q(i,j,k)                  , &
                              mass(i,j,k)               , &
                              imass(i,j,k)              , &
                              ph(i,j,k)                 , &
                              DZET(i,j,k)               , &
                              QDDF3(i,j,k)              , &
                              AA                        , &
                              BB                        , &
                              AREA_ANV_PRC1             , &
                              PRN_above_AN_old          , &
                              PRN_above_AN_new          , &
                              PSN_above_AN_old          , &
                              PSN_above_AN_new          , &
                              EVAP_DD_AN_above_old      , &
                              EVAP_DD_AN_above_new      , &
                              SUBL_DD_AN_above_old      , &
                              SUBL_DD_AN_above_new      , &
                              ANVENVFC                  , &
                              ANVDDRFC                  , & 
                              CONS_ALHF                 , &
                              CONS_ALHS                 , &
                              CONS_ALHL                 , &
                              CONS_CP                   , &
                              CONS_TICE                 , &
                              CONS_H2OMW                , & 
                              CONS_AIRMW                , &
                              REVAP_OFF_P               , &
                              C_ACC                     , &
                              C_EV_R                    , &
                              C_EV_S                    , &
                              RHO_W                     , &
                              ESTBLX                      )
 
         PRN_above_LS_old = PRN_above_LS_new
         PSN_above_LS_old = PSN_above_LS_new
         EVAP_DD_LS_above_old = EVAP_DD_LS_above_new
         SUBL_DD_LS_above_old = SUBL_DD_LS_above_new

         !Precip and Evap for Large Scale
         LSENVFC = 1.0
         call  precipandevap( K                         , &
                              KTOP                      , &
                              LM                        , &
                              DT                        , &
                              FRLAND(i,j)               , &
                              RHCRIT                    , &
                              QRN_LS                    , &
                              QSN_LS                    , &
                              QLT_tmp                   , &
                              QIT_tmp                   , &
                              T(i,j,k)                  , &
                              Q(i,j,k)                  , &
                              mass(i,j,k)               , &
                              imass(i,j,k)              , &
                              ph(i,j,k)                 , &
                              DZET(i,j,k)               , &
                              QDDF3(i,j,k)              , &
                              AA                        , &
                              BB                        , &
                              AREA_LS_PRC1              , &
                              PRN_above_LS_old          , &
                              PRN_above_LS_new          , &
                              PSN_above_LS_old          , &
                              PSN_above_LS_new          , &
                              EVAP_DD_LS_above_old      , &
                              EVAP_DD_LS_above_new      , &
                              SUBL_DD_LS_above_old      , &
                              SUBL_DD_LS_above_new      , &
                              LSENVFC                   , &
                              LSDDRFC                   , & 
                              CONS_ALHF                 , &
                              CONS_ALHS                 , &
                              CONS_ALHL                 , &
                              CONS_CP                   , &
                              CONS_TICE                 , &
                              CONS_H2OMW                , & 
                              CONS_AIRMW                , &
                              REVAP_OFF_P               , &
                              C_ACC                     , &
                              C_EV_R                    , &
                              C_EV_S                    , &
                              RHO_W                     , &
                              ESTBLX                      )


         if ( (QL_ls(i,j,k) + QL_con(i,j,k)) > 0.00 ) then
            QT_tmpi_1 = 1./(QL_ls(i,j,k) + QL_con(i,j,k))
         else
            QT_tmpi_1 = 0.0
         endif
         QL_ls(i,j,k)  = QL_ls(i,j,k)  * QLT_tmp * QT_tmpi_1
         QL_con(i,j,k) = QL_con(i,j,k) * QLT_tmp * QT_tmpi_1
  
         if ( (QI_LS(i,j,k) + QI_con(i,j,k)) > 0.00 ) then
            QT_tmpi_2 = 1./(QI_ls(i,j,k) + QI_con(i,j,k))
         else
            QT_tmpi_2 = 0.0
         endif
         QI_ls(i,j,k)  = QI_ls(i,j,k)  * QIT_tmp * QT_tmpi_2
         QI_con(i,j,k) = QI_con(i,j,k) * QIT_tmp * QT_tmpi_2
                


        !SINK FILTERING
!        if (do_moist_physics == 1) then
!           SINKfilt_ql  = 0.65
!           SINKfilt_qi  = 0.65
!           SINKfilt_CF  = 1.0
!        elseif (do_moist_physics == 2) then
!           SINKfilt_ql  = 0.9
!           SINKfilt_qi  = 0.9
!           SINKfilt_CF  = 1.0
!        endif
!
!        if (k < 50) then
!            qi_ls(i,j,k)  = SINKfilt_qi * qi_ls(i,j,k)  +  (1.0-SINKfilt_qi) * qi_ls_p_presink
!            qi_con(i,j,k) = SINKfilt_qi * qi_con(i,j,k) +  (1.0-SINKfilt_qi) * qi_con_p_presink
!            q(i,j,k)      = SINKfilt_qi * q(i,j,k)      +  (1.0-SINKfilt_qi) * q_p_presink
!        endif
!
!        if ( abs(k-62) .le. 2) then
!           ql_ls(i,j,k)  = SINKfilt_ql * ql_ls(i,j,k)  +  (1.0-SINKfilt_ql) * ql_ls_p_presink
!           ql_con(i,j,k) = SINKfilt_ql * ql_con(i,j,k) +  (1.0-SINKfilt_ql) * ql_con_p_presink
!        endif
!
!        cf_con(i,j,k) = SINKfilt_CF * cf_con(i,j,k) +  (1.0-SINKfilt_CF) * cf_con_p_presink
!
!
!
!        !TOTAL FILTERING
!        TOTfilt_T  = 0.25
!        t(i,j,k)      = TOTfilt_T  * t(i,j,k)      +  (1.0-TOTfilt_T) * t_p_preall
!
!        if (do_moist_physics == 1) then
!           TOTfilt_ql = 0.75
!           TOTfilt_qi = 1.0
!        elseif (do_moist_physics == 2) then
!           TOTfilt_ql = 0.5
!           TOTfilt_qi = 0.5
!        endif
!
!        ql_ls(i,j,k)  = TOTfilt_ql * ql_ls(i,j,k)  +  (1.0-TOTfilt_ql) * ql_ls_p_preall
!        ql_con(i,j,k) = TOTfilt_ql * ql_con(i,j,k) +  (1.0-TOTfilt_ql) * ql_con_p_preall
!
!        qi_ls(i,j,k)  = TOTfilt_qi * qi_ls(i,j,k)  +  (1.0-TOTfilt_qi) * qi_ls_p_preall
!        qi_con(i,j,k) = TOTfilt_qi * qi_con(i,j,k) +  (1.0-TOTfilt_qi) * qi_con_p_preall


   
      endDO
   endDO
endDO

!Clean up of excess relative humidity
RHEXCESS = 1.1
call DQSAT_BAC(DQSDT,QS,T,Ph,im,jm,lm,ESTBLX,CONS_H2OMW,CONS_AIRMW)

where ( Q > RHEXCESS*QS )
   DQS = (Q - RHEXCESS*QS)/( 1.0 + RHEXCESS*DQSDT*CONS_ALHL/CONS_CP )
elsewhere
   DQS = 0.0
endwhere
           
Q = Q - DQS
T = T + (CONS_ALHL/CONS_CP)*DQS

!Clean up Q<0
do j=1,JM
  do i=1,IM

     !Total precipitable water
     TPW = SUM( Q(i,j,:)*DM(i,j,:) )

     NEGTPW = 0.
     do l=1,LM
        if ( Q(i,j,l) < 0.0 ) then 
           NEGTPW   = NEGTPW + ( Q(i,j,l)*DM( i,j,l ) )
           Q(i,j,l) = 0.0
        endif
     enddo

     do l=1,LM
        if ( Q(i,j,l) >= 0.0 ) then 
           Q(i,j,l) = Q(i,j,l)*( 1.0+NEGTPW/(TPW-NEGTPW) )
        endif
     enddo

   end do
end do

!Convert temperature back to potential temperature
th = T/pih

end subroutine cloud_driver



!SUBROUTINES
subroutine cloud_tidy( QV, TE, QLC, QIC, CF, QLA, QIA, AF, CONS_ALHL, CONS_ALHS, CONS_CP )

Implicit None

real(8), intent(inout) :: TE,QV,QLC,CF,QLA,AF,QIC,QIA
real(8), intent(in) :: CONS_ALHL, CONS_ALHS, CONS_CP

 !Fix if Anvil cloud fraction too small
 if (AF < 1.E-5) then
         QV  = QV + QLA + QIA
         TE  = TE - (CONS_ALHL/CONS_CP)*QLA - (CONS_ALHS/CONS_CP)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
 end if

 !Fix if LS cloud fraction too small
! if ( CF < 1.E-5 ) then
!         QV = QV + QLC + QIC
!         TE = TE - (CONS_ALHL/CONS_CP)*QLC - (CONS_ALHS/CONS_CP)*QIC
!         CF  = 0.
!         QLC = 0.
!         QIC = 0.
! end if
      
 !LS LIQUID too small
 if ( QLC  < 1.E-8 ) then
    QV = QV + QLC 
    TE = TE - (CONS_ALHL/CONS_CP)*QLC
    QLC = 0.
 end if
 !LS ICE too small
 if ( QIC  < 1.E-8 ) then
    QV = QV + QIC 
    TE = TE - (CONS_ALHS/CONS_CP)*QIC
    QIC = 0.
 end if

 !Anvil LIQUID too small
 if ( QLA  < 1.E-8 ) then
    QV = QV + QLA 
    TE = TE - (CONS_ALHL/CONS_CP)*QLA
    QLA = 0.
 end if
 !Anvil ICE too small
 if ( QIA  < 1.E-8 ) then
    QV = QV + QIA 
    TE = TE - (CONS_ALHS/CONS_CP)*QIA
    QIA = 0.
 end if

 !Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
 if ( ( QLA + QIA ) < 1.E-8 ) then
    QV = QV + QLA + QIA
    TE = TE - (CONS_ALHL/CONS_CP)*QLA - (CONS_ALHS/CONS_CP)*QIA
    AF  = 0.
    QLA = 0.
    QIA = 0.
 end if
 !Ditto if LS cloud LIQUID+ICE too small
 if ( ( QLC + QIC ) < 1.E-8 ) then
    QV = QV + QLC + QIC
    TE = TE - (CONS_ALHL/CONS_CP)*QLC - (CONS_ALHS/CONS_CP)*QIC
    CF  = 0.
    QLC = 0.
    QIC = 0.
 end if

end subroutine cloud_tidy

subroutine meltfreeze( DT, TE, QL, QI, T_ICE_ALL, T_ICE_MAX, ICEFRPWR, &
                       CONS_ALHL, CONS_ALHS, CONS_CP )

Implicit None

!Inputs
real(8), intent(in) :: DT, T_ICE_ALL, T_ICE_MAX
integer, intent(in) :: ICEFRPWR
real(8), intent(in) :: CONS_ALHL, CONS_ALHS, CONS_CP

!Prognostic
real(8), intent(inout) :: TE,QL,QI

!Locals
real(8)  :: fQi, dQil
real(8), parameter :: taufrz = 1000.

 fQi = 0.0
 dQil = 0.0

 call get_ice_fraction( TE, T_ICE_ALL, T_ICE_MAX, ICEFRPWR, fQi )

 !Freeze liquid
 if ( TE <= T_ICE_MAX ) then
    dQil = Ql *(1.0 - EXP( -Dt * fQi / taufrz ) )
 end if

 dQil = max(  0., dQil )
 Qi   = Qi + dQil
 Ql   = Ql - dQil
 TE   = TE + (CONS_ALHS-CONS_ALHL)*dQil/CONS_CP

 dQil = 0.
 !Melt ice instantly above 0^C
 if ( TE > T_ICE_MAX ) then
    dQil = -Qi 
 end if

 dQil = min(  0., dQil )
 Qi   = Qi + dQil
 Ql   = Ql - dQil
 TE   = TE + (CONS_ALHS-CONS_ALHL)*dQil/CONS_CP

end subroutine meltfreeze


subroutine convec_src( DT, MASS, iMASS, TE, QV, DCF, DMF, QLA, QIA, AF, QS, CONS_ALHS, &
                       CONS_ALHL, CONS_CP, T_ICE_ALL, T_ICE_MAX, ICEFRPWR )

IMPLICIT NONE

!Inputs
real(8), intent(IN) :: DT, T_ICE_ALL, T_ICE_MAX
integer, intent(in) :: ICEFRPWR
real(8), intent(in) :: MASS, iMASS, QS
real(8), intent(in) :: DMF, DCF
real(8), intent(in) :: CONS_ALHS, CONS_ALHL, CONS_CP

!Prognostic
real(8), intent(inout) :: TE, QV
real(8), intent(inout) :: QLA, QIA, AF

!Locals
real(8), parameter :: minrhx = 0.001 !Minimum allowed env RH
real(8) :: TEND, QVx, fQi

 !Namelist
 !DT         - Timestep
 !MASS       - Level Mass
 !iMASS      - 1/Level Mass
 !TE         - Temperature
 !QV         - Specific Humidity
 !DCF        - CNV_DQL from RAS
 !DMF        - CNV_MFD from RAS
 !QLA        - Convective cloud liquid water
 !QIA        - Convective cloud liquid ice
 !AF         - Convective cloud fraction
 !QS         - Qsat

 !Zero out locals
 TEND = 0.0
 QVx = 0.0
 fQi = 0.0

 !Addition of condensate from RAS
 TEND = DCF*iMASS
 call get_ice_fraction( TE, T_ICE_ALL, T_ICE_MAX, ICEFRPWR, fQi )
 QLA  = QLA + (1.0-fQi)* TEND*DT
 QIA  = QIA +      fQi * TEND*DT

 !Convective condensation has never frozen so latent heat of fusion
 TE   = TE +   (CONS_ALHS-CONS_ALHL) * fQi * TEND*DT/CONS_CP

 !Compute Tiedtke-style anvil fraction
 TEND=DMF*iMASS
 AF = AF + TEND*DT
 AF = MIN( AF , 0.99 )

 ! Check for funny (tiny, negative) external QV, resulting from assumed QV=QSAT within anvil.
 ! Simply constrain AF assume condensate just gets "packed" in     
      
 if ( AF < 1.0 ) then
    QVx  = ( QV - QS * AF )/(1.-AF)
 else
    QVx  = QS
 end if

 !If saturated over critial value and there is Anvil
 if ( (( QVx - minrhx*QS ) < 0.0 ) .and. (AF > 0.) ) then
    AF  = (QV  - minrhx*QS )/( QS*(1.0-minrhx) )
 end if

 if ( AF < 0. ) then  ! If still cant make suitable env RH then destroy anvil
    AF  = 0.0
    QV  = QV + QLA + QIA
    TE  = TE - (CONS_ALHL*QLA + CONS_ALHS*QIA)/CONS_CP        
    QLA = 0.0
    QIA = 0.0
 end if

end subroutine convec_src



subroutine pdf_width (PP,FRLAND,maxrhcrit,maxrhcritland,turnrhcrit,minrhcrit,pi_0,ALPHA )  

 Implicit none

 !Inputs
 real(8), intent(in) :: PP, FRLAND

 real(8), intent(in) :: maxrhcrit, maxrhcritland, turnrhcrit, minrhcrit, pi_0

 !Prognostic
 real(8), intent(inout) :: ALPHA

 !Locals
 real(8) :: A1, tempmaxrh

 !Namelist
 !PP             - Pressure (hPa)
 !FRLAND         - Fraction of land
 !maxrhcrit      - Constant
 !maxrhcritland  - Constant
 !turnrhcrit     - Constant
 !minrhcrit      - Constant
 !pi_0           - Constant
 !ALPHA          - 1/2*PDFwidth so RH_crit=1.0-alpha

 ! Use Slingo-Ritter (1985) formulation for critical relative humidity
 ! array a1 holds the critical rh, ranges from 0.8 to 1

 tempmaxrh = maxrhcrit
 if (frland > 0.05) then
    tempmaxrh = maxrhcritland
 endif
 
 a1 = 1.0
 if (pp .le. turnrhcrit) then

    a1 = minrhcrit

 else
    
    a1 = minrhcrit + (tempmaxrh-minrhcrit)/(19.) * &
         ((atan( (2.*(pp- turnrhcrit)/(1020.-turnrhcrit)-1.) * &
         tan(20.*pi_0/21.-0.5*pi_0) ) + 0.5*pi_0) * 21./pi_0 - 1.)

 end if

 a1 = min(a1,1.)
 alpha = 1. - a1

 ALPHA = MIN( ALPHA , 0.25 )  ! restrict RHcrit to > 75% 

end subroutine pdf_width







subroutine ls_cloud( DT, ALPHA, PDFSHAPE, PL, TE, QV, QCl, QAl, QCi, QAi, CF, AF, &
                     CONS_ALHL, CONS_ALHF, CONS_ALHS, CONS_CP, CONS_H2OMW, CONS_AIRMW, T_ICE_ALL, &
                     T_ICE_MAX, ICEFRPWR, ESTBLX, cloud_pertmod, dmp )

IMPLICIT NONE

!Inputs
real(8), intent(in) :: DT, ALPHA, PL, T_ICE_ALL, T_ICE_MAX
integer, intent(in) :: PDFSHAPE, cloud_pertmod, dmp
integer, intent(in) :: ICEFRPWR
real(8), intent(in) :: CONS_ALHL, CONS_ALHF, CONS_ALHS, CONS_CP, CONS_H2OMW, CONS_AIRMW
real(8), intent(in) :: ESTBLX(:)

!Prognostic
real(8), intent(inout) :: TE, QV, QCl, QCi, QAl, QAi, CF, AF

!Locals
integer :: n

real(8) :: QCO, CFO, QAO, QT, QMX, QMN, DQ
real(8) :: TEO, QSx, DQsx, QS, DQs, tmpARR
real(8) :: QCx, QVx, CFx, QAx, QC, QA, fQi, fQi_A, dQAi, dQAl, dQCi, dQCl

real(8) :: TEn, QSp, CFp, QVp, QCp 
real(8) :: TEp, QSn, CFn, QVn, QCn
real(8) :: ALHX, SIGMAQT1, SIGMAQT2

real(8), dimension(1) :: DQSx1,QSx1,TEo1,PL1

 !Namelist
 !DT      - Timestep 
 !ALPHA   - PDF half width
 !PL      - Pressure (hPa)
 !TE      - Temperature
 !QV      - Specific humidity
 !QCl     - Convective cloud liquid water
 !QAl     - Large scale cloud liquid water
 !QCi     - Convective cloud liquid ice
 !QAi     - Large scale cloud liquid ice
 !CF      - Large scale cloud fraction
 !AF      - Convective cloud fraction

 QC = QCl + QCi
 QA = QAl + QAi

 if ( QA > 0.0 ) then
    fQi_A = QAi / QA 
 else
    fQi_A = 0.0
 end if

 TEo = TE

 call DQSATs_BAC(DQSx, QSx,TEo,PL,ESTBLX,CONS_H2OMW, CONS_AIRMW)

 if ( AF < 1.0 ) then
    if (dmp == 1) then
       if ( (1.-af) .gt. 0.02) then
           tmpARR = 1./(1.-AF)
       else
          tmpARR = 1./(1.-AF)
       endif
    elseif (dmp == 2) then
       tmpARR = 1./(1.-AF)
    endif
 else
    tmpARR = 0.0
 end if

 CFx = CF*tmpARR
 QCx = QC*tmpARR
 QVx = ( QV - QSx*AF )*tmpARR

 if ( AF >= 1.0 ) then
    QVx = QSx*1.e-4 
 end if

 if ( AF > 0. ) then
    QAx = QA/AF
 else
    QAx = 0.
 end if

 QT  = QCx + QVx

 TEp = TEo
 QSn = QSx
 TEn = TEo
 CFn = CFx
 QVn = QVx
 QCn = QCx

 DQS = DQSx

 !Begin iteration
 !do n=1,4
 n = 1

    QSp = QSn
    QVp = QVn
    QCp = QCn
    CFp = CFn

    !Dont call again as not looping
    DQS = DQSx
    QSn = QSx
    !call DQSATs_BAC(DQS, QSn, TEn, PL, ESTBLX, CONS_H2OMW, CONS_AIRMW)

    TEp = TEn
    call get_ice_fraction( TEp, T_ICE_ALL, T_ICE_MAX, ICEFRPWR, fQi )

    sigmaqt1  = ALPHA*QSn
    sigmaqt2  = ALPHA*QSn

    if (PDFSHAPE .eq. 2) then
      !For triangular, symmetric: sigmaqt1 = sigmaqt2 = alpha*qsn (alpha is half width)
      !For triangular, skewed r : sigmaqt1 < sigmaqt2
      sigmaqt1  = ALPHA*QSn
      sigmaqt2  = ALPHA*QSn
    endif

    !Compute cloud fraction
    if (cloud_pertmod == 0) then
       call pdffrac(1,qt,sigmaqt1,sigmaqt2,qsn,CFn)
    elseif (cloud_pertmod == 1) then
       call pdffrac(4,qt,sigmaqt1,sigmaqt2,qsn,CFn)
    endif

    !Compute cloud condensate
    call pdfcondensate(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,QCn)

    !Adjustments to anvil condensate due to the assumption of a stationary TOTAL 
    !water PDF subject to a varying QSAT value during the iteration
    if ( AF > 0. ) then
       QAo = QAx  ! + QSx - QS 
    else
       QAo = 0.
    end if

    ALHX = (1.0-fQi)*CONS_ALHL + fQi*CONS_ALHS

    if (PDFSHAPE .eq. 1) then
       QCn = QCp + ( QCn - QCp ) / ( 1. - (CFn * (ALPHA-1.) - (QCn/QSn))*DQS*ALHX/CONS_CP)
    elseif (PDFSHAPE .eq. 2) then
       !This next line needs correcting - need proper d(del qc)/dT derivative for triangular
       !for now, just use relaxation of 1/2.
       if (n.ne.4) then
          QCn = QCp + ( QCn - QCp ) *0.5
       endif
    endif
    
    QVn = QVp - (QCn - QCp)
    TEn = TEp + (1.0-fQi)*(CONS_ALHL/CONS_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF ) &
              +      fQi* (CONS_ALHS/CONS_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF )

 !enddo ! qsat iteration

 CFo = CFn
 CF  = CFn
 QCo = QCn
 TEo = TEn

 ! Update prognostic variables. QCo, QAo become updated grid means.
 if ( AF < 1.0 ) then

    CF  = CFo * ( 1.-AF)
    QCo = QCo * ( 1.-AF)
    QAo = QAo *   AF  

 else !Grid box filled with anvil

    ! Special case AF=1, i.e., box filled with anvil. Note: no guarantee QV_box > QS_box

    CF  = 0.          ! Remove any other cloud
    QAo = QA  + QC    ! Add any LS condensate to anvil type
    QCo = 0.          ! Remove same from LS   

    QT  = QAo + QV    ! Total water

    ! Now set anvil condensate to any excess of total water over QSx (saturation value at top)
    QAo = MAX( QT - QSx, 0. )

 end if

 !Partition new condensate into ice and liquid taking care to keep both >=0 separately. 
 !New condensate can be less than old, so Delta can be < 0.

 QCx  = QCo - QC
 dQCl = (1.0-fQi)*QCx
 dQCi =      fQi *QCx

 !Large Scale Partition
 if ((QCl+dQCl)<0.) then
    dQCi = dQCi + (QCl+dQCl)
    dQCl = -QCl !== dQCl - (QCl+dQCl)
 end if
 if ((QCi+dQCi)<0.) then
    dQCl = dQCl + (QCi+dQCi)
    dQCi = -QCi !== dQCi - (QCi+dQCi)
 end if

 QAx   = QAo - QA
 dQAl  = QAx ! (1.0-fQi)*QAx
 dQAi  = 0.  !  fQi  * QAx

 !Convective partition
 if ((QAl+dQAl)<0.) then
    dQAi = dQAi + (QAl+dQAl)
    dQAl = -QAl
 end if
 if ((QAi+dQAi)<0.) then
    dQAl = dQAl + (QAi+dQAi)
    dQAi = -QAi 
 end if

 ! Clean-up cloud if fractions are too small
 if ( AF < 1.e-5 ) then
    dQAi = -QAi
    dQAl = -QAl
 end if
 if ( CF < 1.e-5 ) then
    dQCi = -QCi
    dQCl = -QCl
 end if

 QAi = QAi + dQAi
 QAl = QAl + dQAl
 QCi = QCi + dQCi
 QCl = QCl + dQCl

 !Update specific humidity
 QV  = QV  - ( dQAi+dQCi+dQAl+dQCl) 

 !Update temperature
 TE = TE + (CONS_ALHL*( dQAi+dQCi+dQAl+dQCl)+CONS_ALHF*(dQAi+dQCi))/ CONS_CP

 !Take care of situations where QS moves past QA during QSAT iteration (QAo negative). 
 !"Evaporate" offending QA
 if ( QAo <= 0. ) then
    QV  = QV + QAi + QAl
    TE  = TE - (CONS_ALHS/CONS_CP)*QAi - (CONS_ALHL/CONS_CP)*QAl
    QAi = 0.
    QAl = 0.
    AF  = 0.  
 end if

end subroutine ls_cloud

subroutine pdffrac (flag,qtmean,sigmaqt1,sigmaqt2,qstar,clfrac)
 
 IMPLICIT NONE

 !Inputs
 INTEGER, INTENT(IN) :: flag            
 real(8), INTENT(IN) :: qtmean, sigmaqt1, sigmaqt2, qstar

 !Prognostic
 real(8), INTENT(INOUT) :: clfrac

 !LOCALS
 REAL(8) :: qtmode, qtmin, qtmax
 REAL(8) :: RH, RHd, q1, q2


 if (flag.eq.1) then !Tophat PDF
 
    if ((qtmean+sigmaqt1).lt.qstar) then
       clfrac = 0.
    else
       if (sigmaqt1.gt.0.) then
          clfrac = min((qtmean + sigmaqt1 - qstar),2.*sigmaqt1)/(2.*sigmaqt1)
       else
          clfrac = 1.
       endif
    endif

 elseif(flag.eq.2) then !Triangular PDF
 
    qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.
    qtmin = min(qtmode-sigmaqt1,0.)
    qtmax = qtmode + sigmaqt2
 
    if (qtmax.lt.qstar) then
       clfrac = 0.
    elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
       clfrac = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
    elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
       clfrac = 1. - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
    elseif ( qstar.le.qtmin ) then
       clfrac = 1.
    endif

 elseif(flag.eq.3) then !TANH function for the perturabtions

    !Tophat PDF for the reference part
    if ((qtmean+sigmaqt1).lt.qstar) then
       clfrac = 0.
    else
       if (sigmaqt1.gt.0.) then
          clfrac = min((qtmean + sigmaqt1 - qstar),2.*sigmaqt1)/(2.*sigmaqt1)
       else
          clfrac = 1.
       endif
    endif

 elseif(flag.eq.4) then !Linear function for the perturabtions

    !Tophat PDF for the reference part
    if ((qtmean+sigmaqt1).lt.qstar) then
       clfrac = 0.
    else
       if (sigmaqt1.gt.0.) then
          clfrac = min((qtmean + sigmaqt1 - qstar),2.*sigmaqt1)/(2.*sigmaqt1)
       else
          clfrac = 1.
       endif
    endif

 endif

end subroutine pdffrac


subroutine pdfcondensate (flag,qtmean4,sigmaqt14,sigmaqt24,qstar4,condensate4)
      
 IMPLICIT NONE

 !Inputs
 INTEGER, INTENT(IN) :: flag            
 real(8), INTENT(IN) :: qtmean4, sigmaqt14, sigmaqt24, qstar4

 !Prognostic
 real(8), INTENT(INOUT) :: condensate4

 !Locals
 real(8) :: qtmode, qtmin, qtmax, constA, constB, cloudf
 real(8) :: term1, term2, term3
 real(8) :: qtmean, sigmaqt1, sigmaqt2, qstar, condensate

 qtmean = qtmean4
 sigmaqt1 = sigmaqt14
 sigmaqt2 = sigmaqt24
 qstar = qstar4

 if (flag.eq.1) then

    if (qtmean+sigmaqt1.lt.qstar) then
       condensate = 0.d0
    elseif (qstar.gt.qtmean-sigmaqt1) then
       if (sigmaqt1.gt.0.d0) then
          condensate = (min(qtmean + sigmaqt1 - qstar,2.d0*sigmaqt1)**2)/ (4.d0*sigmaqt1)
       else
          condensate = qtmean-qstar
       endif
    else
       condensate = qtmean-qstar
    endif

 elseif (flag.eq.2) then

    qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.d0
    qtmin = min(qtmode-sigmaqt1,0.d0)
    qtmax = qtmode + sigmaqt2
    
    if ( qtmax.lt.qstar ) then
       condensate = 0.d0
    elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
       constB = 2.d0 / ( (qtmax - qtmin)*(qtmax-qtmode) )
       cloudf = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
       term1 = (qstar*qstar*qstar)/3.d0
       term2 = (qtmax*qstar*qstar)/2.d0
       term3 = (qtmax*qtmax*qtmax)/6.d0
       condensate = constB * (term1-term2+term3) - qstar*cloudf
    elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
       constA = 2.d0 / ( (qtmax - qtmin)*(qtmode-qtmin) )
       cloudf = 1.d0 - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
       term1 = (qstar*qstar*qstar)/3.d0
       term2 = (qtmin*qstar*qstar)/2.d0
       term3 = (qtmin*qtmin*qtmin)/6.d0
       condensate = qtmean - ( constA * (term1-term2+term3) ) - qstar*cloudf
    elseif ( qstar.le.qtmin ) then
       condensate = qtmean-qstar
    endif

 elseif (flag.eq.3) then

    !Reference part done normally
    !IF (qtmean + sigmaqt1 .LT. qstar) THEN
    !  condensate = 0.d0
    !ELSE IF (qstar .GT. qtmean - sigmaqt1) THEN
    !  IF (sigmaqt1 .GT. 0.d0) THEN
    !    IF (qtmean + sigmaqt1 - qstar .GT. 2.d0*sigmaqt1) THEN
    !      min1 = 2.d0*sigmaqt1
    !    ELSE
    !      min1 = qtmean + sigmaqt1 - qstar
    !    END IF
    !    condensate = min1**2/(4.d0*sigmaqt1)
    !  ELSE
    !    condensate = qtmean - qstar
    !  END IF
    !ELSE
    !  condensate = qtmean - qstar
    !END IF

    !Perturbation part from linear
    if (qtmean - qstar > -0.5e-3) then
       condensate = qtmean - qstar
    else
       condensate = 0.0
    endif

 endif
 
 condensate4 = condensate

end subroutine pdfcondensate



subroutine evap_cnv( DT, RHCR, PL, TE, QV, QL, QI, F, XF, QS, RHO_W, CLD_EVP_EFF, &
                     CONS_H2OMW, CONS_AIRMW, CONS_ALHL, CONS_RVAP, CONS_RGAS, CONS_PI, CONS_CP)

IMPLICIT NONE

!Inputs
real(8), intent(in) :: DT, RHCR, PL, XF, QS, RHO_W, CLD_EVP_EFF
real(8), intent(in) :: CONS_H2OMW, CONS_AIRMW, CONS_ALHL, CONS_RVAP, CONS_RGAS, CONS_PI, CONS_CP

!Prognostics
real(8), intent(inout) :: TE, QV, QL, QI, F

!Locals
real(8) :: ES, RADIUS, K1, K2, TEFF, QCm, EVAP, RHx, QC, A_eff, EPSILON

real(8), parameter :: K_COND  =  2.4e-2
real(8), parameter :: DIFFU   =  2.2e-5
real(8), parameter :: NN = 50.*1.0e6

 EPSILON = CONS_H2OMW/CONS_AIRMW
 A_EFF = CLD_EVP_EFF

 !EVAPORATION OF CLOUD WATER.
 ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100 <-^ convert from mbar to Pa)

 RHx = MIN( QV/QS , 1.00 )

 K1 = (CONS_ALHL**2) * RHO_W / ( K_COND*CONS_RVAP*(TE**2))
 K2 = CONS_RVAP * TE * RHO_W / ( DIFFU * (1000./PL) * ES )
 !Here DIFFU is given for 1000 mb so 1000./PR accounts for increased diffusivity at lower pressure. 

 if ( ( F > 0.) .and. ( QL > 0. ) ) then
    QCm=QL/F
 else
    QCm=0.
 end if

 call LDRADIUS(PL,TE,QCm,NN,RHO_W,RADIUS,CONS_RGAS,CONS_PI)

 if ( (RHx < RHCR ) .and.(RADIUS > 0.0) ) then
    TEFF   =   (RHCR - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
 else
    TEFF   = 0.0 ! -999.
 end if

 EVAP = a_eff*QL*DT*TEFF
 EVAP = MIN( EVAP , QL  )

 QC=QL+QI
 if (QC > 0.) then
    F = F * ( QC - EVAP ) / QC
 end if

 QV   = QV + EVAP
 QL   = QL - EVAP
 TE   = TE - (CONS_ALHL/CONS_CP)*EVAP

end subroutine evap_cnv

subroutine subl_cnv( DT, RHCR, PL, TE, QV, QL, QI, F, XF, QS, RHO_W, CLD_EVP_EFF, &
                     CONS_H2OMW, CONS_AIRMW, CONS_ALHL, CONS_RVAP, CONS_RGAS, CONS_PI, CONS_CP, CONS_ALHS)

IMPLICIT NONE

!INPUTS
real(8), intent(in) :: DT, RHCR, PL, XF, QS, RHO_W, CLD_EVP_EFF
real(8), intent(in) :: CONS_H2OMW, CONS_AIRMW, CONS_ALHL, CONS_RVAP, CONS_RGAS, CONS_PI, CONS_CP, CONS_ALHS

!PROGNOSTIC
real(8), intent(inout) :: TE, QV, QL, QI, F

!LOCALS
real(8) :: ES, RADIUS, K1, K2, TEFF, QCm, SUBL, RHx, QC, A_eff, NN, EPSILON

real(8), parameter :: K_COND  =  2.4e-2
real(8), parameter :: DIFFU   =  2.2e-5

 EPSILON =  CONS_H2OMW/CONS_AIRMW

 A_EFF = CLD_EVP_EFF

 NN = 5.*1.0e6

 ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100 s <-^ convert from mbar to Pa)

 RHx = MIN( QV/QS , 1.00 )

 K1 = (CONS_ALHL**2) * RHO_W / ( K_COND*CONS_RVAP*(TE**2))
 K2 = CONS_RVAP * TE * RHO_W / ( DIFFU * (1000./PL) * ES )

 !Here DIFFU is given for 1000 mb so 1000./PR accounts for increased diffusivity at lower pressure.
 if ( ( F > 0.) .and. ( QI > 0. ) ) then
    QCm=QI/F
 else
    QCm=0.
 end if

 call LDRADIUS(PL,TE,QCm,NN,RHO_W,RADIUS,CONS_RGAS,CONS_PI)

 if ( (RHx < RHCR ) .and.(RADIUS > 0.0) ) then
    TEFF   =   ( RHCR - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
 else
    TEFF   = 0.0 ! -999.
 end if

 SUBL = a_eff*QI*DT*TEFF
 SUBL = MIN( SUBL , QI  )

 QC=QL+QI
 if (QC > 0.) then
    F    = F * ( QC - SUBL ) / QC
 end if

 QV   = QV   + SUBL
 QI   = QI   - SUBL
 TE   = TE   - (CONS_ALHS/CONS_CP)*SUBL

end subroutine subl_cnv

subroutine LDRADIUS(PL,TE,QCL,NN,RHO_W,RADIUS,CONS_RGAS,CONS_PI)

IMPLICIT NONE

!Inputs
real(8), intent(in) :: TE,PL,NN,QCL,RHO_W
real(8), intent(in) :: CONS_RGAS,CONS_PI

!Outputs      
real(8), intent(out) :: RADIUS

!Equiv. Spherical Cloud Particle Radius in m
RADIUS = ((QCL * (100.*PL / (CONS_RGAS*TE )))/(NN*RHO_W*(4./3.)*CONS_PI))**(1./3.)

end subroutine LDRADIUS

subroutine autoconversion_ls( DT, QC, QP, TE, PL, F, SUNDQV2, SUNDQV3, SUNDQT1, &
                           C_00, LWCRIT, DZET)

IMPLICIT NONE

!Inputs
real(8), intent(in) :: DT, TE, PL, DZET, SUNDQV2, SUNDQV3, SUNDQT1, C_00, LWCRIT

!Prognostic
real(8), intent(inout) :: QC, QP, F

!Locals
real(8) :: ACF0, ACF, C00x, iQCcrx, F2, F3, RATE, dQP, QCm, dqfac

 !Zero Locals
 ACF0 = 0.0
 ACF = 0.0
 C00x = 0.0
 iQCcrx = 0.0
 F2 = 0.0
 F3 = 0.0
 RATE = 0.0
 dQP = 0.0
 QCm = 0.0
 dqfac = 0.0

 CALL cons_sundq3(TE, SUNDQV2, SUNDQV3, SUNDQT1, F2, F3 )

 C00x  = C_00 * F2 * F3
 iQCcrx = F2 * F3 / LWCRIT

 if ( ( F > 0.) .and. ( QC > 0. ) )then
    QCm = QC/F
 else
    QCm = 0.
 end if

 RATE = C00x * ( 1.0 - EXP( - ( QCm * iQCcrx )**2 ) )

 !Temporary kluge until we can figure a better to make thicker low clouds.
 F2 = 1.0
 F3 = 1.0 

 !Implement ramps for gradual change in autoconv

 !Thicken low high lat clouds
 if ( PL .GE. 775.  .AND. TE .LE.  275. ) then
    !F3 = max(-0.016 * PL + 13.4, 0.2)
    F3 = 0.2
 end if
 if ( PL .GE. 825.  .AND. TE .LE.  282. ) then
    !F3 = max(0.11 * TE - 30.02, 0.2)
    F3 = 0.2
 end if
 if ( PL .GE. 775.  .AND. PL .LT. 825. .AND. TE .LE.  282. .AND. TE .GT. 275.) then
    !F3 = min(max(-0.016*PL + 0.11 * TE - 16.85, 0.2),1.)
    F3 = 0.2
 end if
 if ( PL .GE. 825.  .AND. TE .LE.  275. ) then
    F3 = 0.2
 end if
 if ( PL .LE. 775.  .OR. TE .GT.  282. ) then
    F3 = 1.
 end if

 !Thin-out low tropical clouds
 if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
    F3 = min(0.2 * TE - 56, 2.)
 end if
 if ( PL .GE. 925.  .AND. TE .GE.  290. ) then
    F3 = min(0.04 * PL - 36., 2.)
 end if
 if ( PL .GE. 925.  .AND. PL .LT. 950. .AND. TE .GT.  285. .AND. TE .LT. 290.) then
    F3 = max(min(0.04*PL + 0.2 * TE - 94., 2.),1.)
 end if
 if ( PL .GE. 950.  .AND. TE .GE.  290. ) then
    F3 = 2.
 end if

 F3   = MAX( F3, 0.1 )
 RATE = F3 * RATE

 dQP  =  QC*( 1.0 - EXP( -RATE * DT ) )
 dQP  =  MAX( dQP , 0.0 )  ! Protects against floating point problems for tiny RATE


 !Wipe-out warm fogs
 dqfac = 0.
 if ( PL .GE. 975.  .AND. TE .GE.  280. ) then
    dqfac = max(min(0.2 * TE - 56., 1.),0.)
 end if
 if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
    dqfac = max(min(0.04 * PL - 38., 1.),0.)
 end if
 if ( PL .GE. 950.  .AND. PL .LT. 975. .AND. TE .GT.  280. .AND. TE .LT. 285.) then
    dqfac = max(min(0.04*PL + 0.2 * TE - 95., 1.),0.)
 end if
 if ( ( PL >= 975. ) .AND. (TE >= 285. ) ) then
    dqfac = 1.
 end if

 dQP = max(dQP, dqfac*QC)

 QC = QC - dQP  
 QP = QP + dQP

 !IF LARGE SCALE THEN
 if ( ((QC + dQP) > 0.) ) then
    F = QC * F / (QC + dQP )
 end if


END SUBROUTINE autoconversion_ls

subroutine autoconversion_cnv( DT, QC, QP, TE, PL, F, SUNDQV2, SUNDQV3, SUNDQT1, &
                           C_00, LWCRIT, DZET)

IMPLICIT NONE

!Inputs
real(8), intent(in) :: DT, TE, PL, DZET, SUNDQV2, SUNDQV3, SUNDQT1, C_00, LWCRIT

!Prognostic
real(8), intent(inout) :: QC, QP, F

!Locals
real(8) :: ACF0, ACF, C00x, iQCcrx, F2, F3, RATE, dQP, QCm, dqfac

 !Zero Locals
 ACF0 = 0.0
 ACF = 0.0
 C00x = 0.0
 iQCcrx = 0.0
 F2 = 0.0
 F3 = 0.0
 RATE = 0.0
 dQP = 0.0
 QCm = 0.0
 dqfac = 0.0

 CALL cons_sundq3(TE, SUNDQV2, SUNDQV3, SUNDQT1, F2, F3 )

 C00x  = C_00 * F2 * F3
 iQCcrx = F2 * F3 / LWCRIT

 if ( ( F > 0.) .and. ( QC > 0. ) )then
    QCm = QC/F
 else
    QCm = 0.
 end if

 RATE = C00x * ( 1.0 - EXP( - ( QCm * iQCcrx )**2 ) )

 !Temporary kluge until we can figure a better to make thicker low clouds.
 F2 = 1.0
 F3 = 1.0 

 !Implement ramps for gradual change in autoconv

 !Thicken low high lat clouds
 if ( PL .GE. 775.  .AND. TE .LE.  275. ) then
    !F3 = max(-0.016 * PL + 13.4, 0.2)
    F3 = 0.2
 end if
 if ( PL .GE. 825.  .AND. TE .LE.  282. ) then
    !F3 = max(0.11 * TE - 30.02, 0.2)
    F3 = 0.2
 end if
 if ( PL .GE. 775.  .AND. PL .LT. 825. .AND. TE .LE.  282. .AND. TE .GT. 275.) then
    !F3 = min(max(-0.016*PL + 0.11 * TE - 16.85, 0.2),1.)
    F3 = 0.2
 end if
 if ( PL .GE. 825.  .AND. TE .LE.  275. ) then
    F3 = 0.2
 end if
 if ( PL .LE. 775.  .OR. TE .GT.  282. ) then
    F3 = 1.
 end if

 !Thin-out low tropical clouds
 if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
    F3 = min(0.2 * TE - 56, 2.)
 end if
 if ( PL .GE. 925.  .AND. TE .GE.  290. ) then
    F3 = min(0.04 * PL - 36., 2.)
 end if
 if ( PL .GE. 925.  .AND. PL .LT. 950. .AND. TE .GT.  285. .AND. TE .LT. 290.) then
    F3 = max(min(0.04*PL + 0.2 * TE - 94., 2.),1.)
 end if
 if ( PL .GE. 950.  .AND. TE .GE.  290. ) then
    F3 = 2.
 end if

 F3   = MAX( F3, 0.1 )
 RATE = F3 * RATE

 dQP  =  QC*( 1.0 - EXP( -RATE * DT ) )
 dQP  =  MAX( dQP , 0.0 )  ! Protects against floating point problems for tiny RATE


 !Wipe-out warm fogs
 dqfac = 0.
 if ( PL .GE. 975.  .AND. TE .GE.  280. ) then
    dqfac = max(min(0.2 * TE - 56., 1.),0.)
 end if
 if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
    dqfac = max(min(0.04 * PL - 38., 1.),0.)
 end if
 if ( PL .GE. 950.  .AND. PL .LT. 975. .AND. TE .GT.  280. .AND. TE .LT. 285.) then
    dqfac = max(min(0.04*PL + 0.2 * TE - 95., 1.),0.)
 end if
 if ( ( PL >= 975. ) .AND. (TE >= 285. ) ) then
    dqfac = 1.
 end if

 dQP = max(dQP, dqfac*QC)

 QC = QC - dQP  
 QP = QP + dQP

END SUBROUTINE autoconversion_cnv

subroutine get_ice_fraction(TEMP, T_ICE_ALL, T_ICE_MAX, ICEFRPWR, ICEFRCT)

IMPLICIT NONE

!Inputs
real(8), intent(in) :: TEMP, T_ICE_ALL, T_ICE_MAX
integer, intent(in) :: ICEFRPWR

!Outputs
real(8), intent(out) :: ICEFRCT


 ICEFRCT  = 0.00
 if ( TEMP <= T_ICE_ALL ) then
    ICEFRCT = 1.000
 else if ( (TEMP > T_ICE_ALL) .AND. (TEMP <= T_ICE_MAX) ) then
    ICEFRCT = 1.00 -  ( TEMP - T_ICE_ALL ) / ( T_ICE_MAX - T_ICE_ALL ) 
 end if
 
 ICEFRCT = MIN(ICEFRCT,1.00)
 ICEFRCT = MAX(ICEFRCT,0.00)

 ICEFRCT = ICEFRCT**ICEFRPWR


end subroutine get_ice_fraction


SUBROUTINE cons_sundq3(TEMP,RATE2,RATE3,TE1, F2, F3)

IMPLICIT NONE

!Inputs
real(8), intent(in) :: RATE2, RATE3, TE1, TEMP

!Outputs
real(8), intent(out) :: F2, F3

!Locals
real(8), parameter :: TE0 = 273.
real(8), parameter :: TE2 = 200.
real(8) :: JUMP1


 JUMP1=  (RATE2-1.0) / ( ( TE0-TE1 )**0.333 ) 

 !Ice - phase treatment 
 IF ( TEMP .GE. TE0 ) THEN
    F2   = 1.0
    F3   = 1.0
 END IF
 IF ( ( TEMP .GE. TE1 ) .AND. ( TEMP .LT. TE0 ) ) THEN
    if (abs(TE0 - TEMP) .gt. 0.0) then  !Linearisation security
       F2   = 1.0 + JUMP1 * (( TE0 - TEMP )**0.3333)
    else
       F2 = 1.0
    endif
    F3   = 1.0
 END IF
 IF ( TEMP .LT. TE1 ) THEN
    F2   = RATE2 + (RATE3-RATE2)*(TE1-TEMP)/(TE1-TE2)
    F3   = 1.0
 END IF

 F2 = MIN(F2,27.0)


end subroutine cons_sundq3


subroutine cons_microphys(TEMP,PR,Q_SAT,AA,BB,CONS_H2OMW,CONS_AIRMW,CONS_RVAP,ALHX3)

IMPLICIT NONE

!Inputs
real(8), intent(in) :: TEMP, Q_sat, PR, ALHX3
real(8), intent(in) :: CONS_H2OMW, CONS_AIRMW, CONS_RVAP

!Outputs
real(8), intent(out) :: AA, BB

!Locals
real(8), parameter :: K_COND  =  2.4e-2
real(8), parameter :: DIFFU   =  2.2e-5
real(8) :: E_SAT, EPSI

 EPSI = CONS_H2OMW/CONS_AIRMW

 E_SAT = 100.* PR * Q_SAT /( (EPSI) + (1.0-(EPSI))*Q_SAT )  ! (100 converts from mbar to Pa)
   
 AA  = ( ALHX3**2 ) / ( K_COND*CONS_RVAP*(TEMP**2) )
 BB  = CONS_RVAP*TEMP / ( DIFFU*(1000./PR)*E_SAT )

end subroutine cons_microphys



subroutine cons_alhx(T,ALHX3,T_ICE_MAX,T_ICE_ALL,CONS_ALHS,CONS_ALHL)

IMPLICIT NONE

!Inputs
real(8), intent(in) :: T, T_ICE_MAX, T_ICE_ALL 
real(8),  intent(in) :: CONS_ALHS, CONS_ALHL

!Outputs
real(8), intent(out) :: ALHX3

 if ( T < T_ICE_ALL ) then
    ALHX3 = CONS_ALHS
 end if

 if ( T > T_ICE_MAX ) then
    ALHX3 = CONS_ALHL
 end if

 if ( (T <= T_ICE_MAX) .and. (T >= T_ICE_ALL) ) then
    ALHX3 = CONS_ALHS + (CONS_ALHL - CONS_ALHS)*( T - T_ICE_ALL ) /( T_ICE_MAX - T_ICE_ALL )
 end if

end subroutine cons_alhx

subroutine marshpalm(RAIN,PR,DIAM3,NTOTAL,W,VE)

Implicit None

!Inputs
real(8), intent(in ) :: RAIN, PR     ! in kg m^-2 s^-1, mbar

!Outputs
real(8), intent(out) :: DIAM3, NTOTAL, W, VE

!Locals
integer :: IQD
real(8), parameter  :: N0 = 0.08  !cm^-3
real(8) :: RAIN_DAY,SLOPR,DIAM1

real(8) :: RX(8) , D3X(8)

 RAIN_DAY = 0.0
 SLOPR = 0.0
 DIAM1 = 0.0

 !Marshall-Palmer sizes at different rain-rates: avg(D^3)
 !RX = (/ 0.   , 5.   , 20.  , 80.  , 320. , 1280., 5120., 20480. /)  ! rain per in mm/day
 RX(1) = 0.
 RX(2) = 5.
 RX(3) = 20.
 RX(4) = 80.
 RX(5) = 320.
 RX(6) = 1280.
 RX(7) = 5120.
 RX(8) = 20480.

 !D3X= (/ 0.019, 0.032, 0.043, 0.057, 0.076, 0.102, 0.137, 0.183  /)
 D3X(1) = 0.019
 D3X(2) = 0.032
 D3X(3) = 0.043
 D3X(4) = 0.057
 D3X(5) = 0.076
 D3X(6) = 0.102
 D3X(7) = 0.137
 D3X(8) = 0.183

 RAIN_DAY = RAIN * 3600. *24.

 IF ( RAIN_DAY <= 0.00 ) THEN
    DIAM1 = 0.00
    DIAM3 = 0.00
    NTOTAL= 0.00
    W     = 0.00
 END IF

 DO IQD = 1,7
    IF ( (RAIN_DAY <= RX(IQD+1)) .AND. (RAIN_DAY > RX(IQD) ) ) THEN
       SLOPR =( D3X(IQD+1)-D3X(IQD) ) / ( RX(IQD+1)-RX(IQD) )
       DIAM3 = D3X(IQD) + (RAIN_DAY-RX(IQD))*SLOPR
    END IF
 END DO

 IF ( RAIN_DAY >= RX(8) ) THEN
    DIAM3=D3X(8)
 END IF

 NTOTAL = 0.019*DIAM3

 DIAM3  = 0.664 * DIAM3  

 W      = (2483.8 * DIAM3 + 80.)*SQRT(1000./PR)

 VE     = MAX( 0.99*W/100. , 1.000 )

 DIAM1  = 3.0*DIAM3

 DIAM1  = DIAM1/100.
 DIAM3  = DIAM3/100.
 W      = W/100.
 NTOTAL = NTOTAL*1.0e6

end subroutine marshpalm

subroutine ice_settlefall_cnv( WXR, QI, PL, TE, F, CONS_RGAS, KHu, KHl, k, DT, DZ, QP, ANV_ICEFALL_C )

IMPLICIT NONE

!Inputs
real(8), intent(in) :: WXR, PL, TE, DZ, DT, ANV_ICEFALL_C
real(8), intent(in) :: CONS_RGAS
integer, intent(in) :: KHu, KHl, k

real(8), intent(inout) :: QI, F, QP

!Locals
real(8) :: RHO, XIm, LXIm, QIxP, VF


 RHO = 1000.*100.*PL/(CONS_RGAS*TE)  ! 1000 TAKES TO g m^-3 ; 100 takes mb TO Pa

 if ( ( F > 0.) .and. ( QI > 0. ) ) then
    XIm = (QI/F)*RHO
 else
    XIm = 0.
 end if

 if ( XIm > 0.) then
    LXIm = LOG10(XIm)
 else
    LXIm = 0.0
 end if

 VF = 128.6 + 53.2*LXIm + 5.5*LXIm**2

 !VF = VF*100./MAX(PL,10.) ! Reduce/increase fall speeds for high/low pressure (NOT in LC98!!! ) 
 ! Assume unmodified they represent situation at 100 mb

 if (WXR > 0.) then
    VF = VF * ( 100./MAX(PL,10.) )**WXR
 endif

 VF = VF/100.

 if (KHu .gt. 0 .and. KHl .gt. 0) then
    if ( (k-1 .ge. KHu) .and. (k-1 .le. KHl)  ) then
       VF = 0.01 * VF
    end if
 endif

 VF = ANV_ICEFALL_C * VF

 QIxP = 0.0

 QIxP = QI * ( VF * DT / DZ ) 
 QIxP = MIN( QIxP , QI )

 QIxP = MAX( QIxP, 0.0 ) ! protects against precision problem

 QP = QP + QIxP
 QI = QI - QIxP



end SUBROUTINE ice_settlefall_cnv


subroutine ice_settlefall_ls( WXR, QI, PL, TE, F, CONS_RGAS, KHu, KHl, k, DT, DZ, QP, LS_ICEFALL_C )

IMPLICIT NONE

!Inputs
real(8), intent(in) :: WXR, PL, TE, DZ, DT, LS_ICEFALL_C
real(8), intent(in) :: CONS_RGAS
integer, intent(in) :: KHu, KHl, k

real(8), intent(inout) :: QI, F, QP

!Locals
real(8) :: RHO, XIm, LXIm, QIxP, VF

 RHO = 1000.*100.*PL/(CONS_RGAS*TE)  ! 1000 TAKES TO g m^-3 ; 100 takes mb TO Pa

 if ( ( F > 0.) .and. ( QI > 0. ) ) then
    XIm = (QI/F)*RHO
 else
    XIm = 0.
 end if

 if ( XIm > 0.) then
    LXIm = LOG10(XIm)
 else
    LXIm = 0.0
 end if

 if (abs(XIm) .gt. 0.0) then !Linearisation security
    VF = 109.0*(XIm**0.16)
 else
    VF = 0.0
 endif

 !VF = VF*100./MAX(PL,10.) ! Reduce/increase fall speeds for high/low pressure (NOT in LC98!!! ) 
 ! Assume unmodified they represent situation at 100 mb

 if (WXR > 0.) then
    VF = VF * ( 100./MAX(PL,10.) )**WXR
 endif

 VF = VF/100.

 if (KHu .gt. 0 .and. KHl .gt. 0) then
    if ( (k-1 .ge. KHu) .and. (k-1 .le. KHl)  ) then
       VF = 0.01 * VF
    end if
 endif

 VF = LS_ICEFALL_C * VF

 QIxP = 0.0

 QIxP = QI * ( VF * DT / DZ ) 
 QIxP = MIN( QIxP , QI )

 QIxP = MAX( QIxP, 0.0 ) ! protects against precision problem

 QP = QP + QIxP
 QI = QI - QIxP

 if ( ((QI + QIxP) > 0.) ) then
    F = QI * F / (QI + QIxP )
 end if


end SUBROUTINE ice_settlefall_ls


subroutine precipandevap(K, KTOP, LM, DT, FRLAND, RHCR3, QPl, QPi, QCl, QCi, TE, QV, mass, imass, &
                            PL, dZE, QDDF3, AA, BB, AREA, PFl_above_in, PFl_above_out, PFi_above_in, PFi_above_out, &
                            EVAP_DD_above_in, EVAP_DD_above_out, SUBL_DD_above_in, SUBL_DD_above_out, &
                            ENVFC, DDRFC,  &
                            CONS_ALHF, CONS_ALHS, CONS_ALHL, CONS_CP, CONS_TICE,CONS_H2OMW,CONS_AIRMW, REVAP_OFF_P, &
                            C_ACC, C_EV_R, C_EV_S, RHO_W, ESTBLX )

IMPLICIT NONE

!Inputs
integer, intent(in) :: K, LM, KTOP
real(8), intent(in) :: DT, mass, imass, PL, AA, BB, RHCR3, dZE, QDDF3, AREA, FRLAND, ENVFC, DDRFC
real(8), intent(in) :: CONS_ALHF, CONS_ALHS, CONS_ALHL, CONS_CP, CONS_TICE, CONS_H2OMW, CONS_AIRMW
real(8), intent(in) :: REVAP_OFF_P
real(8), intent(in) :: C_ACC, C_EV_R, C_EV_S, RHO_W
real(8), intent(in) :: ESTBLX(:)

!Prognostics
real(8), intent(inout) :: QV, QPl, QPi, QCl, QCi, TE
real(8), intent(inout) :: PFl_above_in, Pfl_above_out, PFi_above_in, PFi_above_out
real(8), intent(inout) :: EVAP_DD_above_in, EVAP_DD_above_out, SUBL_DD_above_in, SUBL_DD_above_out

!Locals
integer :: NS, NSMX, itr,L

real(8) :: PFi, PFl, QS, dQS, ENVFRAC, TKo, QKo, QSTKo, DQSTKo, RH_BOX, T_ED, QPlKo, QPiKo
real(8) :: Ifactor, RAINRAT0, SNOWRAT0, FALLRN, FALLSN, VEsn, VErn, NRAIN, NSNOW, Efactor
real(8) :: TinLAYERrn, DIAMrn, DROPRAD, TinLAYERsn, DIAMsn, FLAKRAD
real(8) :: EVAP, SUBL, ACCR, MLTFRZ, EVAPx, SUBLx, EVAP_DD,SUBL_DD,DDFRACT, LANDSEAF
real(8) :: TAU_FRZ, TAU_MLT

real(8), parameter :: TRMV_L = 1.0   !m/s
logical, parameter :: taneff = .false.

!Fraction of precip falling through "environment" vs through cloud
real(8), parameter :: B_SUB = 1.00

 ENVFRAC = ENVFC

 IF ( AREA > 0. ) THEN
    Ifactor = 1./ ( AREA )
 ELSE
    Ifactor = 1.00
 END if

 Ifactor = MAX( Ifactor, 1.) !

 !Start at top of precip column:
 !
 !   a) Accrete                   
 !   b) Evaporate/Sublimate  
 !   c) Rain/Snow-out to next level down 
 !   d) return to (a)

 !Update saturated humidity
 call DQSATs_BAC(dQS,QS,TE,PL,ESTBLX,CONS_H2OMW,CONS_AIRMW)

 DDFRACT = DDRFC 

 IF (K == KTOP) THEN

    PFl=QPl*MASS
    PFi=QPi*MASS

    EVAP_DD = 0.
    SUBL_DD = 0.

 ELSE 

    QPl = QPl + PFl_above_in * iMASS
    PFl = 0.00

    QPi = QPi + PFi_above_in * iMASS
    PFi = 0.00

    ACCR = B_SUB * C_ACC * ( QPl*MASS ) *QCl   

    ACCR = MIN(  ACCR , QCl  )

    QPl     = QPl + ACCR
    QCl     = QCl - ACCR

    !Accretion of liquid condensate by falling ice/snow
    ACCR = B_SUB * C_ACC * ( QPi*MASS ) *QCl   

    ACCR = MIN(  ACCR , QCl  )

    QPi = QPi + ACCR
    QCl = QCl - ACCR

    !! Liquid freezes when accreted by snow
    TE  = TE + CONS_ALHF*ACCR/CONS_CP

    RAINRAT0 = Ifactor*QPl*MASS/DT
    SNOWRAT0 = Ifactor*QPi*MASS/DT

    call MARSHPALM(RAINRAT0,PL,DIAMrn,NRAIN,FALLrn,VErn)
    call MARSHPALM(SNOWRAT0,PL,DIAMsn,NSNOW,FALLsn,VEsn)

    TinLAYERrn = dZE / ( FALLrn+0.01 )
    TinLAYERsn = dZE / ( FALLsn+0.01 )

    !Melting of Frozen precipitation      
    TAU_FRZ = 5000.  ! time scale for freezing (s). 

    MLTFRZ = 0.0
    IF ( (TE > CONS_TICE ) .and.(TE <= CONS_TICE+5. ) ) THEN
       MLTFRZ = TinLAYERsn * QPi *( TE - CONS_TICE ) / TAU_FRZ 
       MLTFRZ = MIN( QPi , MLTFRZ )
       TE  = TE  - CONS_ALHF*MLTFRZ/CONS_CP
       QPl = QPl + MLTFRZ
       QPi = QPi - MLTFRZ
    END IF

    MLTFRZ = 0.0
    IF ( TE > CONS_TICE+5.  ) THEN  ! Go Ahead and melt any snow/hail left above 5 C 
       MLTFRZ= QPi 
       TE  = TE  - CONS_ALHF*MLTFRZ/CONS_CP
       QPl = QPl + MLTFRZ
       QPi = QPi - MLTFRZ
    END IF

    MLTFRZ = 0.0
    if ( K >= LM-1 ) THEN
       IF ( TE > CONS_TICE+0.  ) THEN ! Go Ahead and melt any snow/hail left above 0 C in lowest layers 
          MLTFRZ= QPi 
          TE  = TE  - CONS_ALHF*MLTFRZ/CONS_CP
          QPl = QPl + MLTFRZ
          QPi = QPi - MLTFRZ
       END IF
    endif

    !Freezing of liquid precipitation      
    MLTFRZ = 0.0
    IF ( TE <= CONS_TICE ) THEN
       TE  = TE + CONS_ALHF*QPl/CONS_CP
       QPi = QPl + QPi
       MLTFRZ = QPl
       QPl = 0.
    END IF

    !In the exp below, evaporation time scale is determined "microphysically" from temp, 
    !press, and drop size. In this context C_EV becomes a dimensionless fudge-fraction. 
    !Also remember that these microphysics are still only for liquid.

    QKo   = QV
    TKo   = TE
    QPlKo = QPl
    QPiKo = QPi

    !do itr = 1,1
    itr = 1

       DQSTKo = dQS
       QSTKo  = QS  + DQSTKo * ( TKo - TE )
       QSTKo  = MAX( QSTKo , 1.0e-7 )

       RH_BOX = QKo/QSTKo

       QKo   = QV
       TKo   = TE

       IF ( RH_BOX < RHCR3 ) THEN
          Efactor =  RHO_W * ( AA + BB )    / (RHCR3 - RH_BOX )
       else
          Efactor = 9.99e9
       end if

       if ( FRLAND < 0.1 ) then
          LANDSEAF = 0.5  ! Over Ocean
       else
          LANDSEAF = 0.5  ! Over Land
       end if

       LANDSEAF = 1.00

       !Rain falling
       if ( (RH_BOX < RHCR3) .AND. (DIAMrn > 0.00) .AND. (PL > 100.) .AND. (PL < REVAP_OFF_P) ) then
          DROPRAD=0.5*DIAMrn
          T_ED =  Efactor * DROPRAD**2 
          T_ED =  T_ED * ( 1.0 + DQSTKo*CONS_ALHL/CONS_CP )
          EVAP =  QPl*(1.0 - EXP( -C_EV_R * VErn * LANDSEAF * ENVFRAC * TinLAYERrn / T_ED ) )
       ELSE
          EVAP = 0.0
       END if

       !Snow falling
       if ( (RH_BOX < RHCR3) .AND. (DIAMsn > 0.00) .AND. (PL > 100.) .AND. (PL < REVAP_OFF_P) ) then
          FLAKRAD=0.5*DIAMsn
          T_ED =  Efactor * FLAKRAD**2   
          T_ED =  T_ED * ( 1.0 + DQSTKo*CONS_ALHS/CONS_CP )
          SUBL =  QPi*(1.0 - EXP( -C_EV_S * VEsn * LANDSEAF * ENVFRAC * TinLAYERsn / T_ED ) )
       ELSE
          SUBL = 0.0
       END IF

       !if (itr == 1) then 
       !   EVAPx  = EVAP
       !   SUBLx  = SUBL
       !else
       !   EVAP   = (EVAP+EVAPx) /2.0
       !   SUBL   = (SUBL+SUBLx) /2.0
       !endif

       QKo = QV + EVAP + SUBL
       TKo = TE - EVAP * CONS_ALHL / CONS_CP - SUBL * CONS_ALHS / CONS_CP

    !enddo
      
    QPi  = QPi - SUBL
    QPl  = QPl - EVAP

    !Put some re-evap/re-subl precip in to a \quote{downdraft} to be applied later
    EVAP_DD = EVAP_DD_above_in + DDFRACT*EVAP*MASS 
    EVAP    = EVAP          - DDFRACT*EVAP
    SUBL_DD = SUBL_DD_above_in + DDFRACT*SUBL*MASS 
    SUBL    = SUBL          - DDFRACT*SUBL

    QV   = QV  + EVAP + SUBL
    TE   = TE  - EVAP * CONS_ALHL / CONS_CP - SUBL * CONS_ALHS / CONS_CP

    PFl  = QPl*MASS
    PFi  = QPi*MASS

 end if

 EVAP = QDDF3*EVAP_DD/MASS
 SUBL = QDDF3*SUBL_DD/MASS
 QV   = QV  + EVAP + SUBL
 TE   = TE  - EVAP * CONS_ALHL / CONS_CP - SUBL * CONS_ALHS / CONS_CP

 QPi = 0.
 QPl = 0.

 PFl_above_out = PFl
 PFi_above_out = Pfi

 EVAP_DD_above_out = EVAP_DD
 SUBL_DD_above_out = SUBL_DD

end subroutine precipandevap


subroutine DQSAT_BAC(DQSi,QSSi,TEMP,PLO,im,jm,lm,ESTBLX,CONS_H2OMW,CONS_AIRMW)
!COMPUTES SATURATION VAPOUR PRESSURE QSSi AND GRADIENT w.r.t TEMPERATURE DQSi.
!INPUTS ARE TEMPERATURE AND PLO (PRESSURE AT T-LEVELS)
!VALES ARE COMPUTED FROM LOOK-UP TALBE (PIECEWISE LINEAR)

 IMPLICIT NONE

 !Inputs
 integer :: im,jm,lm
 real(8), dimension(im,jm,lm) :: TEMP, PLO
 real(8) :: ESTBLX(:)
 real(8)  :: CONS_H2OMW, CONS_AIRMW

 !Outputs
 real(8), dimension(im,jm,lm) :: DQSi, QSSi

 !Locals
 real(8), parameter :: MAX_MIXING_RATIO = 1.0
 real(8) :: ESFAC

 integer :: i, j, k

 real(8) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 integer :: IT

 integer, parameter :: DEGSUBS    =  100
 real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
 integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

 ESFAC = CONS_H2OMW/CONS_AIRMW

 do K=1,LM
    do J=1,JM
       do I=1,IM

          TL = TEMP(I,J,K)
          PL = PLO(I,J,K)

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

          DQSi(I,J,K) = DQSAT
          QSSi(I,J,K) = QSAT

       end do
    end do
 end do

end subroutine DQSAT_BAC

subroutine DQSATs_BAC(DQSi,QSSi,TEMP,PLO,ESTBLX,CONS_H2OMW,CONS_AIRMW)
!COMPUTES SATURATION VAPOUR PRESSURE QSSi AND GRADIENT w.r.t TEMPERATURE DQSi.
!INPUTS ARE TEMPERATURE AND PLO (PRESSURE AT T-LEVELS)
!VALES ARE COMPUTED FROM LOOK-UP TALBE (PIECEWISE LINEAR)

 IMPLICIT NONE

 !Inputs
 real(8) :: TEMP, PLO
 real(8) :: ESTBLX(:)
 real(8)  :: CONS_H2OMW,CONS_AIRMW

 !Outputs
 real(8) :: DQSi, QSSi

 !Locals
 real(8), parameter :: MAX_MIXING_RATIO = 1.0
 real(8) :: ESFAC

 real(8) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 integer :: IT

 integer, parameter :: DEGSUBS    =  100
 real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
 integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

 ESFAC = CONS_H2OMW/CONS_AIRMW

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

end subroutine DQSATs_BAC

END MODULE CLOUD
