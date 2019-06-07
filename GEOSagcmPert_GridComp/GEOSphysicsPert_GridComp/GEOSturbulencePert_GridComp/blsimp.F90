MODULE BLSIMP

!This subroutine takes the model prognostic variables as its inputs and produces
!the three components of the tri-diagonal matrix used to compute BL diffusion. 
!The same diffussivity Kh is applied to all model variables above the surface. 
!Diffusion coefficient is computed using a simplified Louis type approximation.

IMPLICIT NONE

PRIVATE
PUBLIC BL_simp

CONTAINS

subroutine BL_simp( IM, JM, LM, DT, UT, VT, PTT, QVT, QLT, QIT, PET, PKT, FROCEAN, &
                    MAPL_GRAV, MAPL_VIREPS, MAPL_KAPPA, MAPL_CP, MAPL_RGAS, &
                    AKQ, AKS, AKV, BKQ, BKS, BKV, CKQ, CKS, CKV &
                  )

IMPLICIT NONE


!INPUTS
INTEGER, INTENT(IN)                        :: IM, JM, LM
REAL,    INTENT(IN)                        :: DT
REAL,    INTENT(IN), DIMENSION(IM,JM,LM)   :: UT, VT, PTT, QVT, QLT, QIT
REAL,    INTENT(IN), DIMENSION(IM,JM,0:LM) :: PET
REAL,    INTENT(IN), DIMENSION(IM,JM,LM)   :: PKT
REAL,    INTENT(IN), DIMENSION(IM,JM)      :: FROCEAN

REAL,    INTENT(IN)                        :: MAPL_GRAV, MAPL_VIREPS, MAPL_KAPPA, MAPL_CP, MAPL_RGAS

!OUTPUTS
REAL, INTENT(OUT), DIMENSION(IM,JM,LM)     :: AKQ, AKS, AKV, BKQ, BKS, BKV, CKQ, CKS, CKV

!LOCALS
INTEGER :: I, J, L
REAL :: WS, RI, RIN
REAL :: ELL2, CDRAG_LAND, CDRAG_SEA, CDRAG, TCOEF
REAL :: DMI, PKH, DZ, CKX, TVB, TVE, TVT

REAL, DIMENSION(IM,JM) :: BKSQ, BKSV, BKST
REAL, DIMENSION(LM)    :: Kh

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

 !These two are not used
 AKS(:,:,1)  = 0.0
 CKS(:,:,LM) = 0.0

 ELL2       = 30.*30.
 CDRAG_LAND = 0.002
 CDRAG_SEA  = 0.0015

 ! the following lines solves  ------------------------------------
 ! dQ/dT=-G/dP*Rho*CDrag*WS*Q, Q_(t+dt)=Q_t+dt/(Rho_t*dZ_t)*dTAU
 !                             dTAU=TAU_(z) - TAU_(z-1),   
 !                             TAU_(z)= Rho_t*CDrag*WS*Q_(t+dt), 
 !                             TAU_(z-1)=Rho_t*K_(z-1)/dZ*dQ_(t+dt)
 !                             dZ=RT*dP/(P*G), DMI=(G*dt)/dP
 ! ----------------------------------------------------------------
 
 do j=1,jm
    do i=1,im

       DMI    = (MAPL_GRAV*DT)/(PET(i,j,1)-PET(i,j,0))
       TVT    = PTT(I,J,1)*PKT(I,J,1) &
                 * (1.0 + MAPL_VIREPS *QVT(I,J,1) - QLT(I,J,1) - QIT(I,J,1) )
       do L=2,LM
          PKH      = PET(I,J,L)**MAPL_KAPPA
          DZ       = MAPL_CP*(PTT(I,J,L-1)*(PKH-PKT(I,J,L-1)) + PTT(I,J,L  )*(PKT(I,J,L )-PKH))
          WS       = (UT(I,J,L-1)-UT(I,J,L))**2+(VT(I,J,L-1)-VT(I,J,L))**2 + 0.01
          RI       = MAPL_GRAV*((PTT(I,J,L-1)-PTT(I,J,L)) / (0.5*(PTT(I,J,L-1)+PTT(I,J,L))))*DZ/WS
          RIN      = ELL2*SQRT(WS)/DZ
          IF (RI < 0.) THEN
             KH(L) = MAX(0.01, RIN*SQRT(1.-18.*RI))
          ELSE
             KH(L) = MAX(0.01, RIN/(1.+10.*RI*(1.+8.*RI)))
          ENDIF
     
          TVB      = PTT(I,J,L)*PKT(I,J,L)*(1.0 + MAPL_VIREPS *QVT(I,J,L) - QLT(I,J,L) - QIT(I,J,L))
          TVE      = 0.5*(TVT+TVB)

          TVT      = TVB 
          CKX      = -KH(L)*PET(i,j,L)/( MAPL_RGAS * TVE )/DZ
          CKS(I,J,L-1) = CKX * DMI
          DMI      = (MAPL_GRAV*DT)/(PET(i,j,L)-PET(i,j,L-1))
          AKS(I,J,L)   = CKX * DMI
          BKS(I,J,L-1) = 1.0 - (AKS(I,J,L-1)+CKS(I,J,L-1))

          if (L==LM) then
             BKS(I,J,L) = 1.0 - (AKS(I,J,L)+CKS(I,J,L))
             BKSQ(I,J)   = BKS(I,J,L)
             WS     = SQRT(UT(I,J,L)**2 + VT(I,J,L)**2 + 1.0)
             IF (FROCEAN(I,J).eq.1.0) then
                CDRAG=CDRAG_SEA
                TCOEF=1.! CDRAG       ! assume sea surface T unperturbed
             else
                CDRAG=CDRAG_LAND
                TCOEF=0.          ! assume land surface T same as air T
             endif
             KH(L)  = -CDRAG*DMI*WS*PET(I,J,L)/(MAPL_RGAS*TVB)
             BKSV(I,J)   = 1.0 - (AKS(I,J,L)+CKS(I,J,L)+KH(L))
             BKST(I,J)   = 1.0 - (AKS(I,J,L)+CKS(I,J,L)+KH(L)*TCOEF)
          endif

       end do

    end do
 end do

 !Make copy for wind and tracers as not doing differently for this
 !simplified boundary layer.
 AKV = AKS
 BKV = BKS
 CKV = CKS
 AKQ = AKS
 BKQ = BKS
 CKQ = CKS

 !Apply surface layer diffusivities
 BKV(:,:,LM) = BKSV(:,:)
 BKS(:,:,LM) = BKST(:,:)
 BKQ(:,:,LM) = BKSQ(:,:)

endsubroutine BL_simp


ENDMODULE BLSIMP
