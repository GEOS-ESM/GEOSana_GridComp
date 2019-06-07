#include <params.h>
      subroutine surfvar(jlat,taux,tauy,shflx,cflx1,hs,vk,virtcon,tnlev,
     +                   us,vs,thvs,ths,shs,plks,pls,rgas,grav,cp,
     +                   pzero,qbustr,qs,qqs,
     +                   dus,dvs,dthvs,strts,dw2s,ris,irun)

c-----This subroutine calculates the Brunt Vaisala frequencies, and the
c     shears at the top of the surface layer, and sets the value of
c     the turbulent kinetic energy and the tk velocity at the top of the SL
c     
c-----input
       real taux(irun)   ! x-surface stress
       real tauy(irun)   ! y-surface stress
       real shflx(irun)  ! surface sensible heat flux
       real cflx1(irun)  ! surface constuent flux
       real hs(irun)     ! height of surface layer
       real us(irun)     ! u wind at nlev
       real vs(irun)     ! v wind at nlev
       real tnlev(irun)  ! temperature at nlev
       real thvs(irun)   ! model virtual potential temperature at nlev
       real ths(irun)    ! model potential temperature at nlev
       real shs(irun)    ! specific humidity at nelv
       real plks(irun),pls(irun) ! pressure  and pressure **k at last mid level
       real vk           ! Von Karman constant
       real grav,  cp, qbustr,rgas
       integer irun
       integer jlat      ! index for the lat looping
c------output
       real qs(irun)     ! tk velocity at nlev
       real qqs(irun)    ! TKE at nlev
       real strts(irun)  ! Brunt-Vaisala frequency at nlev
       real dw2s(irun)   ! Shear at nlev
       real ris(irun)    ! Richardson number at nlev
c-----local workspace
      real tvs, tmp1,temp1
c      real tvs, gvk, tmp1,tmp2
      real part1,  part2, part3
      real ustar(irun)      !surface friction velocity 
      real ustar2(irun)      !ustar**2 
      real rrhosur(irun)    !1/rho on the last midlevel before ground
      real thstar(irun)     ! temperature scale
      real qstar(irun)      ! humidity scale
c      real temp1(irun)      ! temporany space
c      real temp2(irun)
      real virtfac(irun)   ! 1+0.61sh
      real akhs(irun)      ! 1./(k*hs)
      real zeta(irun)      ! hs/L, with L MOnin Ob. length
      real phim(irun)      ! similarity functions for momentum
      real phih(irun)      ! similarity functions for heat and moisture
      real dus(irun)       ! du/dz at top of surface layer
      real dvs(irun)       ! dv/dz at top of surface layer
      real dthvs(irun)     ! dthv/dz at top of surface layer

      real mol             ! Monin Obukov length
      real hvir            ! Virtual sensible heat flux
      real mpe             !prevents overflow error if division by zero
      real pzero           ! pzero ** kappa



c-----mpe
      mpe = 1.e-06
c-----some needed computation
      do i=1,irun
         akhs(i)=1./(vk*hs(i))
         virtfac(i) = 1. + virtcon*shs(i)
      enddo

c-----First calculates surface friction velocity
c     ustar**2=sqrt(taux**2+tauy**2)/rhosur(last midlevel before ground)
c     and sets the TKE at the ground

      do i=1,irun
      rrhosur(i)=rgas*thvs(i)*plks(i)/pls(i)
      ustar2(i)=sqrt(taux(i)**2 + tauy(i)**2)*rrhosur(i)
      ustar(i)=sqrt(ustar2(i))
      qs(i) = qbustr*ustar(i)
      qqs(i)=0.5*qs(i)*qs(i)
      enddo

c***************CMT momentanamente cambiato per farlo come lsm***********************************
c-----zeta = hs/L, with L monin Obukov length
c      gvk=grav*vk
c      do i = 1,irun
c      tvs = tnlev(i)*virtfac(i)
c      tmp1 = (gvk * hs(i))/(tvs*ustar2(i)*ustar(i))
c      tmp2 = temp1(i) + (0.61 * tnlev(i) * temp2(i))
c      zeta(i) = - tmp1 * tmp2
c      enddo
c      write(6,*) 'zeta',zeta(1)
c***************CMT momentanamente cambiato per farlo come lsm***********************************
c-----Temperature and moisture scales
     
      do i=1,irun
      temp1=(shflx(i)*rrhosur(i))/cp 
      thstar(i) = -1.* temp1/ustar(i)
      qstar(i) = -1.*cflx1(i)*rrhosur(i)/ustar(i)


c-----zeta=hs/L calculated as in the LSM model (for consistency)

            hvir = shflx(i) + virtcon*cp*tnlev(i) *cflx1(i)
            tvs = tnlev(i)*virtfac(i)
            tmp1 = vk * (grav/tvs) * hvir*rrhosur(i)/cp
            if (abs(tmp1) .le. mpe) tmp1 = mpe
            mol = -1. * ustar2(i)*ustar(i)/ tmp1
               zeta(i) = min(hs(i)/mol, 1.)

c-----similarity function
      if (zeta(i). lt. 0.) then
         argo = (1.- 16.*zeta(i))
         phim(i) = 1./(argo**0.25)
         phih(i) = 1./(argo**0.5)
      else if ((zeta(i). ge.0.). and.(zeta(i). le.1.)) then
         phim(i) = 1. + 5.*zeta(i)
         phih(i) = 1. + 5.*zeta(i)
      else if (zeta(i).gt.1.) then
         write(*,*) 'zeta greater than one in surfval'
         stop
      endif
      enddo


c-----gradients
      do i = 1, irun
      dus(i) = ustar(i) * phim(i)*akhs(i)
      dvs(i) = vs(i)*dus(i)
      dus(i) = us(i)*dus(i)
      dw2s(i) = dus(i)*dus(i)+dvs(i)*dvs(i)
      if (dw2s(i).le.1.e-4) dw2s(i)=1.e-4
       part1  = phih(i)*akhs(i)
       part2 = (thstar(i)*virtfac(i))/pzero
       part3 = part2 + virtcon*ths(i)*qstar(i) 
       dthvs(i)=part1*part3
       strts(i) = (grav/thvs(i))*dthvs(i)
       ris(i) = strts(i)/dw2s(i)
      enddo

      return
      end
     
    

      subroutine initurb(kh,km,qqe,qe,p3,qq,q,ibitstb,intq,lminq,
     +lminq1,irun,nlev)

c-----Subroutine to initialize some of the variables needed in trbflx 
#include<implicit.h>
      real km(irun,nlev)       ! Diffusion coefficients for momentum
      real kh(irun,nlev)       ! Diffusion coefficient for heat
      real qqe(irun,nlev-1)    ! TKE from level 2.0 model 
      real qe(irun,nlev-1)     ! TK velocity from  level 2.0 model 
      real p3(irun,nlev)       ! Production term for tke prognostic equation
      real qq(irun,nlev)       ! TKE
      real q(irun,nlev)        ! TK velocity
      
      integer ibitstb(irun,nlev) ! flag for growing or decaying turbulence
                                   ! initialize to zero
      integer intq(irun,nlev)    ! =1 for qq gt e-8, 0 for qq le e-8
      
      integer lminq             !first level from top where qq gt e-8
      integer lminq1            !lminq+1

      integer ibit              ! local workspace

      integer irun, nlev, l, i  ! grid 

c-----Initialize diffusion coefficients kh,km; equilibrium value of
c     turbulent kinetic energy qqe, production terms.
c     Eliminate small qq, select levels where tke is greater than threshold
c     find the first level for which tke is different than zero.

c-----Initialize the diffusion coefficients kh,km
c      do l=1,nlev-1
      do l=1,nlev
         do i = 1,irun
            kh(i,l)=0.
            km(i,l)=0.
         enddo
       enddo

c-----Initialize equilibrium value  of tke(qqe) and production term
      do l=1,nlev-1
         do i=1,irun
            qqe(i,l)=0.
            qe(i,l)=0.
            p3(i,l)=0.
            ibitstb(i,l)=0
          enddo
      enddo

c-----Eliminate small qq and looks for the levels where tke is greater
c     than e-8
      do l=1,nlev-1
         do i=1,irun
            if(qq(i,l).gt.1.e-8) then
              intq(i,l) = 1
            else
              intq(i,l)=0
              qq(i,l)=0.
              q(i,l)=0.
            endif
          enddo
      enddo

c------Determine first level (from top) where qq is different than zero, lminq
C
       DO 300 LMINQ = 1,nlev-1
        IBIT = 0
        DO 9116 I = 1,irun
         IF ( QQ(I,LMINQ).GT.1.e-8 ) IBIT = IBIT + 1
 9116   CONTINUE
        IF(IBIT.GE.1)GO TO 310
 300   CONTINUE
       LMINQ = NLEV-1
 310   CONTINUE
       LMINQ1 = 1
       IF(LMINQ.GT.1)LMINQ1 = LMINQ - 1

       return
       end

       Subroutine trbgrow(irun,nlevml,ibitstb,q,qe,qq,qqe,kh,km,p3)
c-------------------------------------------------------------------------
c
c       Thise subroutine performs three tasks:
c       First it looks the level where qe > 0 and q < qe,(growing turbulence)
c       and for these levels sets ibitstb = 1. For the levels where 
c       growing turbulence is found, it computes the dimensionless
c       diffusion coefficients kh=(q/qe)*Sh and km=(q/qe)*Sm,
c       where Sh and Sm are the dimensionless coefficients calculated 
c       in the level 2 model. In case the TKE qq is less than
c       1 % of the TKE qqe from the 2.0 level model,it forces qq 
c       to 0.01 qqe. The last task is to compute the production term
c       p3. For growing turbulence the production term is 
c       p3=(qe-q)*2*qq/(b1*l). For same strange reasons, here p3
c       is not yet divided by l (length scale), and is not yet
c       multiplied by qq. It will be divided by l later in the computations
c       while the term qq is directly included in the diffusion equation.
c       Therefore in this subroutine p3=2*(qe-q)/b1   

#include <implicit.h>
     
       real b1,b3
       parameter (b1=16.6)
       parameter (b3=1/b1)   

       integer ibitstb(irun,nlevml)  ! flag for occurence of growing
                                     !turbulence.=1 for growing turbulence
                                     !=0 otherwise
       real qq(irun,nlevml)  ! TKE
       real q(irun,nlevml)   ! TK velocity
       real qqe(irun,nlevml) ! TKE from equilibrium model
       real qe(irun,nlevml)  ! TK velocity  from equilibrium model
       real kh(irun,nlevml)  ! dimensionless diffusion coefficients for
                             ! the case of growing turbulence:
       real km(irun,nlevml)  ! Same as above but for momenetum

       real p3(irun,nlevml)  ! Part of the production term expression
       integer irun,nlevml, l, i

       do l=1,nlevml
          do i=1,irun
            if ( (ibitstb(i,l).eq.1) .and. (q(i,l).le.qe(i,l))) then
               ibitstb(i,l)=1
c***************************************************************
                 if ( qe(i,l).eq.0.) then
                  write(6,*) 'qe eq zero'
                  stop
                 endif	
c***************************************************************
               kh(i,l)=(q(i,l)/qe(i,l))*kh(i,l)
               km(i,l)=(q(i,l)/qe(i,l))*km(i,l)
                      if (qq(i,l).le.(0.01*qqe(i,l))) then
                          qq(i,l)=0.01*qqe(i,l)
                          q(i,l)=0.1*qe(i,l)
                      endif
               p3(i,l)=2.*b3*(qe(i,l)-q(i,l))
            else
               ibitstb(i,l)=0
            endif
c
       enddo
       enddo

       return
       end
              

      SUBROUTINE TRBLEN (STRT,DW2,DZ3,Q,VKZE,VKZM,DTHV,DPK,DU,DV,XL,
     1                   QXLM,NLEV,INIT,LMIN,LMINQ,LMINQ1,CP,INT1,INT2,
     2  DZITRP,STBFCN,XL0,Q1,WRKIT1,WRKIT2,WRKIT3,WRKIT4,irun)
C**********************************************************************
C
C  SUBROUTINE TRBLEN - COMPUTES TURBULENT LENGTH SCALE
C                    - CALLED FROM TRBFLX
C   ARGUMENTS ::
C
C     INPUT:
C     ------
C    STRT          -         BRUNT VAISALA FREQUENCY
C    DW2           -         SQUARED SHEAR
C    DZ3           -         LAYER THICKNESS FOR LENGTH SCALE CALC.
C    Q             -         TURBULENCE VELOCITY
C    VKZE          -         VK * Z AT LAYER EDGES
C    VKZM          -         VK * Z AT LAYER CENTERS
C    DTHV          -         VERTICAL GRADIENT OF THV
C    DPK           -         VERTICAL GRADIENT OF PK
C    DU            -         VERTICAL GRADIENT OF U
C    DV            -         VERTICAL GRADIENT OF V
C    NLEV          -         NUMBER OF ATMOSPHERIC LEVELS
C    INIT          -         INPUT FLAG : 1 = INITIAL START
C                                         2 = 2ND CALL FOR INITIAL STAR
C                                         0 = NON-INITIAL START
C
C     OUTPUT:
C     -------
C    XL            -         TURBULENT LENGTH SCALE
C    QXLM          -         TURBULENT LENGTH SCALE * Q AT LAYER CENTER
C    LMIN          -         HIGHEST LAYER WHERE INSTABILITY OCCURS
C    LMINQ         -         HIGHEST LAYER WHERE TURBULENCE OCCURS
C
C   SUBPROGRAMS NEEDED ::
C
C     TRBITP -  INTERPOLATES TO HEIGHT WHERE RI = RICR
C
C**********************************************************************
C
C
      PARAMETER ( RF1     = 0.2340678 )
      PARAMETER ( RF2     = 0.2231172 )
      PARAMETER ( E5      = 49.66     )  
      PARAMETER ( D4      = 2.6532122E-2 )
      PARAMETER ( D1      = D4 * E5   )
      PARAMETER ( RFC     = 0.1912323 )
      PARAMETER ( RICR    = ( (RF1-RFC)*RFC ) / ( (RF2-RFC)*D1 )  )
      PARAMETER ( ALPHA   = 0.1       )
      PARAMETER ( DZCNV   = 100.      )
      PARAMETER ( XL0CNV  = DZCNV * ALPHA )
      PARAMETER ( XL0MIN  = 1.        )
      PARAMETER ( CLMT    = 0.23      )  
      PARAMETER ( CLMT53  = 5. * CLMT / 3. )
 
 
      DIMENSION STRT(irun,NLEV),   DW2(irun,NLEV), DZ3(irun,NLEV)
      DIMENSION DTHV(irun,NLEV),   DPK(irun,NLEV),  DU(irun,NLEV)
      DIMENSION   DV(irun,NLEV),     Q(irun,NLEV)
      DIMENSION VKZM(irun,NLEV-1),VKZE(irun,NLEV-1)
      DIMENSION QXLM(irun,NLEV-1),  XL(irun,NLEV-1)
      DIMENSION DZITRP(irun,nlev-1), STBFCN(irun,nlev)
      DIMENSION    XL0(irun,nlev),       Q1(irun,nlev-1)
      DIMENSION WRKIT1(irun,nlev-1)
      DIMENSION WRKIT2(irun,nlev-1)
      DIMENSION WRKIT3(irun,nlev-1)
      DIMENSION WRKIT4(irun,nlev-1)
C
      INTEGER INT1(irun,nlev), INT2(irun,nlev-1)
C
      NLEVM1 = NLEV - 1
      NLEVP1 = NLEV + 1
      ISTNLV = irun * NLEV
      ISTNM1 = irun * NLEVM1
C
      IF ( INIT.EQ.2 ) GO TO 1200
C
C COMPUTE DEPTHS OF UNSTABLE LAYERS
C =================================
      DO 10 I=1,ISTNLV
          STBFCN(I,1) = STRT(I,1) - RICR * DW2(I,1)
                                INT1(I,1) = 0
      IF( STBFCN(I,1).LE.0. ) INT1(I,1) = 1
   10 CONTINUE
      DO 20 I=1,ISTNM1
                                                INT2(I,1) = 0
!      IF( INT1(I,1).EQ.1 .XOR. INT1(I,2).EQ.1 ) INT2(I,1) = 1
      IF( INT1(I,1).eq.1 .and. INT1(I,2).ne.1 ) INT2(I,1) = 1
      IF( INT1(I,1).ne.1 .and. INT1(I,2).eq.1 ) INT2(I,1) = 1
   20 CONTINUE
C
      DO 40 LMIN = 1,NLEV
      IBIT = 0
      DO 30 I=1,irun
      IBIT = IBIT + INT1(I,LMIN)
   30 CONTINUE
      IF(IBIT.GE.1) GO TO 50
   40 CONTINUE
      LMIN = NLEVP1
   50 CONTINUE
C
      DO 60 I=1,ISTNM1
      XL0(I,1) = 0.
   60 CONTINUE
      DO 70 I=1,irun
      XL0(I,NLEV) = DZ3(I,NLEV)
   70 CONTINUE
C
      IF(LMIN.GE.NLEVP1) GOTO 1100
                         LMIN1 = LMIN - 1
      IF(LMIN1.EQ.0)     LMIN1 = 1
         NLEVML = NLEV - LMIN1
         ISTNML = irun*NLEVML
      CALL TRBITP ( STBFCN(1,LMIN1),INT2(1,LMIN1),DTHV(1,LMIN1),
     .                 DPK(1,LMIN1),  DU(1,LMIN1),  DV(1,LMIN1),
     .              DZITRP(1,LMIN1),  NLEVML,
     .              WRKIT1,WRKIT2,WRKIT3,WRKIT4,CP,irun )
      LP1 = LMIN1 + 1
C
      DO 80 I=1,ISTNML
                                                     INT2(I,LMIN1) = 0
      IF( INT1(I,LMIN1).EQ.1 .OR. INT1(I,LP1).EQ.1 ) INT2(I,LMIN1) = 1
      IF( INT2(I,LMIN1).EQ.1 )
     .     XL0(I,LMIN1) = (0.5+DZITRP(I,LMIN1)) * DZ3(I,LP1)
   80 CONTINUE
      DO 90 I=1,irun
      INT2(I,NLEVM1) = INT1(I,NLEV)
   90 CONTINUE
C
      DO 100 I=1,ISTNML
      IF( INT2(I,LMIN1).EQ.1 ) THEN
      XL0(I,LP1) = XL0(I,LP1) + ( (0.5-DZITRP(I,LMIN1)) * DZ3(I,LP1) )
      ENDIF
  100 CONTINUE
      IF (LMIN.GT.1) GOTO 400
      DO 110 I=1,irun
      IF( INT1(I,1).EQ.1 ) XL0(I,1) = XL0(I,1) + DZ3(I,1)
  110 CONTINUE
  400 CONTINUE
C
      LMINP = LMIN + 1
      IF(LMINP.GT.NLEVM1) GOTO 550
      DO 500 L = LMINP,NLEVM1
           LM1 = L-1
      DO 120 I = 1,irun
      IF( INT1(I,LM1).EQ.1 ) XL0(I,L) = XL0(I,L) + XL0(I,LM1)
  120 CONTINUE
  500 CONTINUE
  550 CONTINUE
      IF(LMIN.GT.NLEVM1) GOTO 600
      DO 130 I = 1,irun
      IF( INT1(I,NLEVM1).EQ.1 .AND. INT1(I,NLEV).EQ.1 ) THEN
      XL0(I,NLEV) = XL0(I,NLEV) + XL0(I,NLEVM1)
      ENDIF
  130 CONTINUE
      IF(LMIN.GT.NLEV) GOTO 1100
  600 CONTINUE
      DO 1000 LL  = LMIN,NLEV
              L   = NLEVM1 + LMIN - LL
              LP1 = L+1
      DO 140 I = 1,irun
      IF( INT1(I,LP1).EQ.1 ) THEN
      IF( INT1(I,L)  .EQ.1 ) THEN
      XL0(I,L) = XL0(I,LP1)
      ELSE
      XL0(I,L) = XL0(I,L) + XL0(I,LP1)
      ENDIF
      ENDIF
  140 CONTINUE
 1000 CONTINUE
 1100 CONTINUE
C
      DO 150 I = 1,ISTNLV
      IF( XL0(I,1).LT.XL0CNV ) XL0(I,1) = XL0CNV
  150 CONTINUE
C
C *********************************************************************
C ****          DETERMINE MIXING LENGTHS FOR STABLE LAYERS          ***
C *********************************************************************
C
      IF(INIT.EQ.1) GOTO 1400
C
      IF(LMINQ.GT.1) THEN
                    ISTLMQ = irun * LMINQ1
      DO 160 I  = 1,ISTLMQ
      INT2(I,1) = 1 - INT1(I,1)
  160 CONTINUE
      ENDIF
      IF(LMINQ.LT.NLEV) THEN
        ISTNMQ = irun * (NLEV-LMINQ)
      DO 170 I = 1,ISTNMQ
      IF( INT1(I,LMINQ).EQ.0 ) THEN
           XL0(I,LMINQ) =         Q(I,LMINQ) / XL0(I,LMINQ)
           XL0(I,LMINQ) =       XL0(I,LMINQ) * XL0(I,LMINQ) + 1.0E-20
           XL0(I,LMINQ) =    STBFCN(I,LMINQ) + XL0(I,LMINQ)
           XL0(I,LMINQ) = SQRT( XL0(I,LMINQ) )
           XL0(I,LMINQ) =         Q(I,LMINQ) / XL0(I,LMINQ)
      ENDIF
                                   INT2(I,LMINQ)  = 0
      IF( XL0(I,LMINQ).LT.XL0MIN ) INT2(I,LMINQ)  = 1
  170 CONTINUE
      ENDIF
C
 1200 CONTINUE
C
      IF(INIT.EQ.2) THEN
      DO 180 I = 1,ISTNM1
      INT2(I,1) = 1 - INT1(I,1)
  180 CONTINUE
      ENDIF
      DO 190 I = 1,ISTNM1
      IF( INT2(I,1).EQ.1 ) XL0(I,1) = XL0MIN
  190 CONTINUE
C
C *********************************************************************
C ****             LENGTH SCALE XL FROM XL0 AND VKZE              ****
C *********************************************************************
C
 1400 CONTINUE
C
      DO 200 I = 1,ISTNM1
      XL(I,1) = XL0(I,1) * VKZE(I,1) / ( XL0(I,1)+VKZE(I,1) )
  200 CONTINUE
C
C *********************************************************************
C ****       CLMT53 TIMES Q TIMES LENGTH SCALE AT MID LEVELS        ***
C *********************************************************************
C
      IF(INIT.EQ.1) GOTO 1700
                   ISTNMQ = irun * (NLEV-LMINQ1)
      DO 210 I = 1,ISTNMQ
                                              Q1(I,LMINQ1) = Q(I,LMINQ1)
                                            INT1(I,LMINQ1) = 0
      IF(    Q(I,LMINQ1).LE.Q(I,LMINQ1+1) ) INT1(I,LMINQ1) = 1
      IF( INT1(I,LMINQ1).EQ.1 )  THEN
           XL0(I,LMINQ1) = XL0(I,LMINQ1+1)
            Q1(I,LMINQ1) =   Q(I,LMINQ1+1)
      ENDIF
  210 CONTINUE
C
      DO 240 I = 1,ISTNMQ
      QXLM(I,LMINQ1) =   XL0(I,LMINQ1)*VKZM(I,LMINQ1)
     .               / ( XL0(I,LMINQ1)+VKZM(I,LMINQ1) )
      QXLM(I,LMINQ1) = CLMT53 * Q1(I,LMINQ1)*QXLM(I,LMINQ1)
  240 CONTINUE
C
 1700 CONTINUE
C
      RETURN
      END
      SUBROUTINE TRBITP ( STBFCN,INTCHG,DTHV,DPK,DU,DV,DZITRP,NLEV,
     .                    AAA,BBB,CCC,DDD,CP,irun )
C**********************************************************************
C
C  SUBROUTINE TRBITP - INTERPOLATES TO THE HEIGHT WHERE RI EQUALS RICR
C                    - CALLED FROM TRBLEN
C   ARGUMENTS ::
C
C     INPUT:
C     ------
C    STBFCN        -         DTHV * DPK - RICR*( DU*DU + DV*DV)
C    INTCHG        -         INT '1' AT LEVELS WHERE STBFCN CHANGES SIG
C    DTHV          -         VERTICAL GRADIENT OF THV
C    DPK           -         VERTICAL GRADIENT OF PK
C    DU            -         VERTICAL GRADIENT OF U
C    DV            -         VERTICAL GRADIENT OF V
C    NLEV          -         NUMBER OF LEVELS TO BE PROCESSED
C
C     OUTPUT:
C     -------
C    DZITRP        -         INTERPOLATION COEFFICIENT
C
C**********************************************************************
C
C
      PARAMETER ( RF1     = 0.2340678 )
      PARAMETER ( RF2     = 0.2231172 )
      PARAMETER ( E5      = 49.66     )  
      PARAMETER ( D4      = 2.6532122E-2 )
      PARAMETER ( D1      = D4 * E5   )
      PARAMETER ( RFC     = 0.1912323 )
      PARAMETER ( RICR    = ( (RF1-RFC)*RFC ) / ( (RF2-RFC)*D1 )  )

      DIMENSION STBFCN(irun,NLEV+1), INTCHG(irun,NLEV)
      DIMENSION DTHV  (irun,NLEV+1), DPK   (irun,NLEV+1)
      DIMENSION DU    (irun,NLEV+1), DV    (irun,NLEV+1)
      DIMENSION DZITRP(irun,NLEV+1)
      DIMENSION AAA   (irun,NLEV), BBB   (irun,NLEV)
      DIMENSION CCC   (irun,NLEV), DDD   (irun,NLEV)
C
C *********************************************************************
C ****           QUADRATIC INTERPOLATION OF RI TO RICR VIA          ***
C ****           LINEAR INTERPOLATION OF DTHV, DPK, DU & DV         ***
C *********************************************************************
C
      ISTNLV = irun*NLEV
      DO 10 I=1,ISTNLV
      DZITRP(I,1) = 0.
   10 CONTINUE
      DO 20 I=1,ISTNLV
      IF( INTCHG(I,1).EQ.1 ) THEN
      DDD(I,1) = (         CP  *(DTHV(I,2)*DPK(I,1)
     .                         + DTHV(I,1)*DPK(I,2)) )
     .         - ( (2.*RICR)  * ( DU(I,2)* DU(I,1)
     .                         +   DV(I,2)* DV(I,1)) )
      AAA(I,1) = STBFCN(I,1) + STBFCN(I,2)
      BBB(I,1) = STBFCN(I,1) - STBFCN(I,2)
      CCC(I,1) =            1. / BBB(I,1)
      DZITRP(I,1) = AAA(I,1) * CCC(I,1)
         AAA(I,1) = AAA(I,1) - DDD(I,1)
         DDD(I,1) =        (   DDD(I,1) *    DDD(I,1) )
     .            - 4. * (STBFCN(I,2) * STBFCN(I,1) )
         DDD(I,1) = DDD(I,1)*CCC(I,1)*CCC(I,1)
         DDD(I,1) =    SQRT( DDD(I,1) )
      ENDIF
C
      IF( INTCHG(I,1).EQ.1 .AND. AAA(I,1).NE.0. ) THEN
      DZITRP(I,1) = ( BBB(I,1)*(1.-DDD(I,1)) ) / AAA(I,1)
      ENDIF
C
      DZITRP(I,1) = 0.5 * DZITRP(I,1)
   20 CONTINUE
C
      RETURN
      END

      SUBROUTINE TRBL20 (RI,STRT,DW2,XL,ZKM,ZKH,QE,QQE,INTSTB,NLEV,
     1                        nlay,irun)
C**********************************************************************
C
C  SUBROUTINE TRBL20 - COMPUTES QE AND DIMLESS COEFS FROM
C                       MELLOR-YAMADA LEVEL 2 MODEL
C                    - CALLED FROM AND FROM TRBFLX
C   ARGUMENTS ::
C
C     INPUT:
C     ------
C    RI            -         RICHARDSON NUMBER
C    STRT          -         BRUNT VAISALA FREQUENCY
C    DW2           -         SQUARED SHEAR
C    XL            -         TURBULENT LENGTH SCALE
C    NLEV          -         NUMBER OF LEVELS TO BE PROCESSED
C
C     OUTPUT:
C     -------
C    ZKM           -         MOMENTUM TRANSPORT COEFFICIENT
C    ZKH           -         HEAT TRANSPORT COEFFICIENT
C    QE            -         EQUILIBRIUM TURBULENT VELOCITY SCALE
C    QQE           -         EQUILIBRIUM TURBULENT KINETIC ENERGY
C    BITSTB        -         BIT '1' WHERE QE GREATER THAN ZERO
C
C**********************************************************************
C
C
 
      PARAMETER ( B1      =   16.6    )  
      PARAMETER ( B2      =   10.1    )
      PARAMETER ( D3      = 0.29397643 )
      PARAMETER ( RF1     = 0.2340678 )
      PARAMETER ( RF2     = 0.2231172 )
      PARAMETER ( D3B2    = D3 / RF1  )
      PARAMETER ( D2      = RF1       )
      PARAMETER ( E5      = 49.66     )  
      PARAMETER ( D4      = 2.6532122E-2 )
      PARAMETER ( D1      = D4 * E5   )
      PARAMETER ( D1HALF = 0.5 * D1 )
      PARAMETER ( D2HALF = 0.5 * D2 )
      PARAMETER ( RFC     = 0.1912323 )
      PARAMETER ( RICR    = ( (RF1-RFC)*RFC ) / ( (RF2-RFC)*D1 )  )
      PARAMETER ( CH      = 2.5828674 )
      PARAMETER ( CM      = CH / D1   )

 
      DIMENSION RI(irun,NLEV), STRT(irun,NLEV), DW2(irun,NLEV)
      DIMENSION XL(irun,NLEV),  ZKM(irun,NLEV), ZKH(irun,NLEV)
      DIMENSION QE(irun,NLEV),  QQE(irun,NLEV)
      INTEGER INTSTB(irun,nlev)
      DIMENSION     EE(irun,nlay-1),  RF(irun,nlay-1)
C
      ISTNLV = irun * NLEV
C
C *********************************************************************
C ****               COMPUTE FLUX RICHARDSON NUMBER                 ***
C *********************************************************************
C
      DO 10 I=1,ISTNLV
      EE(I,1) = D1HALF * RI(I,1) + D2HALF
      RF(I,1) = EE(I,1)* EE(I,1)
      RF(I,1) = RF(I,1)- D3*RI(I,1)
      RF(I,1) = SQRT( RF(I,1) )
      RF(I,1) = EE(I,1) - RF(I,1)
C
      IF( RI(I,1).LE.1.e-4 .AND. RI(I,1).GE.-1.e-4 ) THEN
          RF(I,1) = D3B2*RI(I,1)
      ENDIF
C
C *********************************************************************
C ****           QE AND DIMENSIONLESS DIFFUSION COEFICIENTS         ***
C ****                   FROM LEVEL 2 CLOSURE MODEL                 ***
C *********************************************************************
C
      IF( RI(I,1).LT.RICR .AND. RF(I,1).LT.RFC ) THEN
      ZKH(I,1) =   ( RFC-RF(I,1) ) / (1.-RF(I,1))
      ZKM(I,1) =                CM * (RF1-RF(I,1))
      ZKM(I,1) = ZKH(I,1)*ZKM(I,1) / (RF2-RF(I,1))
      ZKH(I,1) = CH      *ZKH(I,1)
       QE(I,1) = ZKM(I,1)*DW2(I,1) - ZKH(I,1)*STRT(I,1)
      ENDIF
C
      IF( QE(I,1).LT.1.e-14 ) THEN
      INTSTB(I,1) = 0
          QE(I,1) = 0.
      ELSE
      INTSTB(I,1) = 1
          QE(I,1) = B1*QE(I,1)
          QE(I,1) = SQRT( QE(I,1) )
          QE(I,1) = XL(I,1)*QE(I,1)
      ENDIF
      QQE(I,1) = 0.5 * QE(I,1) * QE(I,1)
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE TRBL25(Q,XL,STRT,DW2,INTSTB,INTQ,ZKM,ZKH,P3,NLEV,
     1                      nlay,irun)
C**********************************************************************
C
C  SUBROUTINE TRBL25 - COMPUTES P3 AND DIMLESS COEFS FROM
C                       MELLOR-YAMADA LEVEL 2.5 MODEL
C                    - CALLED FROM TRBFLX
C
C   ARGUMENTS ::
C
C     INPUT:
C     ------
C    Q             -         TURBULENCE VELOCITY
C    XL            -         TURBULENT LENGTH SCALE
C    STRT          -         BRUNT VAISALA FREQUENCY
C    DW2           -         SQUARED SHEAR
C    BITSTB        -         BIT '1' WHERE QE GREATER THAN ZERO
C    NLEV          -         NUMBER OF LEVELS TO BE PROCESSED
C
C     OUTPUT:
C     -------
C    ZKM           -         MOMENTUM TRANSPORT COEFFICIENT
C    ZKH           -         HEAT TRANSPORT COEFFICIENT
C    P3            -         PRODUCTION RATE OF TURBULENT KINETIC ENERG
C
C**********************************************************************
C
C
 
      PARAMETER ( A1      =   0.92    )
      PARAMETER ( A2      =   0.74    )
      PARAMETER ( A4      = 6. * A1 * A1)
      PARAMETER ( C1      =   0.08    )
      PARAMETER ( A5      = 3.*C1*(-1.) )
      PARAMETER ( A3      = A4 * A5*(-1.) )
      PARAMETER ( B1      =   16.6    )  
      PARAMETER ( B2      =   10.1    )
      PARAMETER ( B3      = 1. / B1  )  
      PARAMETER ( FF2     = 9. * A1 * A2 )
      PARAMETER ( FF3     = (3.*A2*B2) - (9.*A2*A2 )  )
      PARAMETER ( FF4     = (3.*A2*B2) + (12.*A1*A2 )  )

 
      DIMENSION   Q(irun,NLEV),  XL(irun,NLEV), STRT(irun,NLEV)
      DIMENSION DW2(irun,NLEV)
      DIMENSION ZKM(irun,NLEV), ZKH(irun,NLEV),   P3(irun,NLEV)
C
      DIMENSION F2(irun,nlay-1),  F3(irun,nlay-1)
      DIMENSION F4(irun,nlay-1),  XQ(irun,nlay-1)
      INTEGER INTSTB(irun,nlay), INTQ(irun,nlay)
C
      ISTNLV = irun * NLEV
C
C *********************************************************************
C ****           P3 AND DIMENSIONLESS DIFFUSION COEFICIENTS         ***
C ****                  FROM LEVEL 2.5 CLOSURE MODEL                ***
C *********************************************************************
C
      DO 10 I=1,ISTNLV
      IF( INTQ(I,1).EQ.1 .AND. INTSTB(I,1).EQ.0 ) THEN
        XQ(I,1) = XL(I,1)  /      Q(I,1)
        XQ(I,1) = XQ(I,1)  *     XQ(I,1)
      STRT(I,1) = XQ(I,1)  *   STRT(I,1)
       DW2(I,1) = XQ(I,1)  *    DW2(I,1)
        F2(I,1) = 1.+FF2  *   STRT(I,1)
        F3(I,1) = 1.+FF3  *   STRT(I,1)
        F4(I,1) = 1.+FF4  *   STRT(I,1)
       ZKH(I,1) = (F4(I,1) *     F2(I,1))
     .     + A4 * (F3(I,1) *    DW2(I,1))
       ZKH(I,1) = (F2(I,1) + A3*DW2(I,1))
     .          / ZKH(I,1)
       ZKM(I,1) = A1 * (F3(I,1)*ZKH(I,1)+A5)
     .               /  F2(I,1)
       ZKH(I,1) = A2 * ZKH(I,1)
        P3(I,1) = ZKH(I,1)*STRT(I,1) + B3
        P3(I,1) = 2. * ( ZKM(I,1)*DW2(I,1) - P3(I,1) )
        P3(I,1) = P3(I,1)*Q(I,1)
C
      ENDIF
   10 CONTINUE
C
      RETURN
      END
      
      subroutine trbdifq (qq,p3,rhokdz,flxfce,nlevmq,irun)

c     Solves the equation for TKE: dqq=(rho/dz)*(5/3*qxl)(gdt/dp)qq+p3*qq
c     The equation is solved with a backward implicit scheme, using
c     a tridiagonal matrix with elements aa,bb,cc.

 
c-----input
      real p3(irun,nlevmq)       ! production term
      real rhokdz(irun,nlevmq)   ! (rho/dz)(xqlm), triple correlation term
      real flxfce(irun,nlevmq)   ! (gdt)/dp

      integer nlevmq                 ! number of atm level to be processed
      integer irun                   ! number of point in the horizontal

c-----input/output
      real qq(irun,nlevmq+1)    !TKE, updated
 
c-----local space
      real aa(irun,nlevmq)      ! term for  solution of tridiagonal system
!                                      (1+bb+cc+forcing)
      real bb(irun,nlevmq)      ! term for solution of tridiagonal system
!                                      (upper term)
      real cc(irun,nlevmq)      ! term  for solution of tridiagonal system
!                                      (lower term)
c
c-----define tridiagonal matrix
c     (includes implicit backward forcing for qq ( production term) added to aa)
      do i=1,irun
        cc(i,1)=0.
      enddo

      do l = 2,nlevmq
          do i = 1,irun
             cc(i,l) = rhokdz(i,l-1)*flxfce(i,l)
          enddo
       enddo

      do l = 1,nlevmq
         do i = 1, irun
            bb(i,l) = rhokdz(i,l)*flxfce(i,l)
            aa(i,l) = (1. + cc(i,l) + bb(i,l)) - p3(i,l)
         enddo
      enddo
c
c-----solve the tridiagonal system for qq
      call vtri0(aa,bb,cc,qq,qq,nlevmq,irun)

c-----eliminate underflow
      do l = 1,nlevmq
         do i = 1, irun
            if (qq(i,l).lt.1.0e-20) qq(i,l) = 0.
         enddo
      enddo

      return
      end
C
      subroutine trbdif1 (xx,rhokdz,flxfac,sfcflx,nlev,irun)

c     Solves the diffusion equation for temperature and humidity.
c     The equation is solved with a backward implicit scheme, using
c     a tridiagonal matrix with elements aa,bb,cc.

 
c-----input
      real rhokdz(irun,nlev-1)    ! (rho/dz)
      real flxfac(irun,nlev)      ! (gdt)/dp
      real sfcflx(irun)           ! g dt /(dp*p0**k)*(H/cp) 
!                                        sfc heat flux into the lowest model
!                                        g dt /dp * E sfc humidity flux into 
!                                        the lowest model       

      integer nlev                    ! n. of atmospheric level to be processed
      integer irun                    ! number of point in the horizontal

c-----input/output
      real xx(irun,nlev)      ! potential temperature or specific humidity
 
c-----local space
      real aa(irun,nlev)       ! term for  solution of tridiagonal system
!                                     (1+bb+cc)
      real bb(irun,nlev)        ! term for solution of tridiagonal system
!                                      (upper term)
      real cc(irun,nlev)        ! term  for solution of tridiagonal system
!                                      (lower term)
c
c-----define tridiagonal matrix
c     NB: bb(irun,nlev) is set = to zero as kh(nlev) is zero. 

      do i=1,irun
        cc(i,1)=0.
        bb(i,nlev)=0.
      enddo

      do l = 2,nlev
          do i = 1,irun
             cc(i,l) = rhokdz(i,l-1)*flxfac(i,l)
          enddo
       enddo

      do l = 1,nlev-1
         do i = 1, irun
            bb(i,l) = rhokdz(i,l)*flxfac(i,l)
         enddo
      enddo
      do l = 1,nlev
         do i = 1, irun
            aa(i,l) = 1. + cc(i,l) + bb(i,l)
         enddo
      enddo
c
c-----solve the tridiagonal system for xx
      call vtri1(aa,bb,cc,xx,xx,sfcflx,nlev,irun)

      return
      end

      subroutine trbdif2 (xx1,xx2,rhokdz,flxfac,sfcflx1,
     +sfcflx2,nlev,irun)

c     Solves the diffusion equation for the u and v wind component
c     The equation is solved with a backward implicit scheme, using
c     a tridiagonal matrix with elements aa,bb,cc.

 
c-----input
      real rhokdz(irun,nlev-1)    ! (rho/dz)
      real flxfac(irun,nlev)      ! (gdt)/dp
      real sfcflx1(irun)           ! g dt /dp* taux 
      real sfcflx2(irun)           ! g dt /dp* tauy 

      integer nlev                    ! n. of atmospheric level to be processed
      integer irun                    ! number of point in the horizontal

c-----input/output
      real xx1(irun,nlev)      ! u-wind
      real xx2(irun,nlev)      ! v-wind
 
c-----local space
      real aa(irun,nlev)       ! term for  solution of tridiagonal system
!                                     (1+bb+cc)
      real bb(irun,nlev)        ! term for solution of tridiagonal system
!                                      (upper term)
      real cc(irun,nlev)        ! term  for solution of tridiagonal system
!                                      (lower term)
c
c-----define tridiagonal matrix
c     NB: bb(irun,nlev) is set = to zero as kh(nlev) is zero. 

      do i=1,irun
        cc(i,1)=0.
        bb(i,nlev)=0.
      enddo

      do l = 2,nlev
          do i = 1,irun
             cc(i,l) = rhokdz(i,l-1)*flxfac(i,l)
          enddo
       enddo

      do l = 1,nlev-1
         do i = 1, irun
            bb(i,l) = rhokdz(i,l)*flxfac(i,l)
         enddo
      enddo
      do l = 1,nlev
         do i = 1, irun
            aa(i,l) = 1. + cc(i,l) + bb(i,l)
         enddo
      enddo
c
c-----solve the tridiagonal system for  the u-component of wind
      call vtri1(aa,bb,cc,xx1,xx1,sfcflx1,nlev,irun)

c-----solve the tridiagonal system for  the u-component of wind
      call vtri2(aa,bb,cc,xx2,xx2,sfcflx2,nlev,irun)
      return
      end
C
      SUBROUTINE VTRI0 ( A,B,C,F,Y,K,irun)
      DIMENSION A(irun,K),B(irun,K),C(irun,K),Y(irun,K+1)
      DIMENSION F(irun,K)
C
      DO 9000 I = 1,irun
       A(I,1) = 1. / A(I,1)
 9000 CONTINUE
C
      DO 100 L = 2,K
       LM1 = L - 1
       DO 9002 I = 1,irun
        C(I,L) = C(I,L) * A(I,LM1)
        A(I,L) = 1. / ( A(I,L) - B(I,LM1) * C(I,L) )
        F(I,L) = F(I,L) + F(I,LM1) * C(I,L)
 9002  CONTINUE
 100  CONTINUE
C
      DO 200 L = K,1,-1
       DO 9004 I = 1,irun
       Y(I,L) = (F(I,L) + B(I,L) * Y(I,L+1)) * A(I,L)
 9004  CONTINUE
 200  CONTINUE
C
      RETURN
      END
C
C
      subroutine vtri1(a,b,c,f,y,sfcflx,nlev,irun)
C
c-----solves the diffusion equation with explicit surface flux

      real a(irun,nlev)
      real b(irun,nlev)
      real c(irun,nlev)
      real y(irun,nlev),f(irun,nlev)
      real sfcflx(irun)    
      integer irun, nlev,lm1,i,l
      
C
      do i=1,irun
      if (a(i,1).eq.0.) then
      write(6,*) 'in vtri1 a = 0'
      stop
      endif
      a(i,1) = 1./a(i,1)
      enddo

C
      do l=2,nlev
          lm1 = l-1
         do i=1,irun
            c(i,l)=c(i,l)*a(i,lm1)
            a(i,l)=1./(a(i,l)-b(i,lm1)*c(i,l))
            f(i,l) = f(i,l) + f(i,lm1)*c(i,l)
         enddo
      enddo
C
c-----Bottom level: includes surface fluxes
     
      do i=1,irun
         y(i,nlev) = (f(i,nlev) + sfcflx(i))*a(i,nlev)
      enddo
c-----perform back substitution
      do l=nlev-1,1,-1
           do i=1,irun
               Y(I,L) = (F(I,L) + B(I,L) * Y(I,L+1)) * A(I,L)
           enddo
      enddo

      return
      end

      subroutine vtri2 (a,b,c,f,y,sfcflx,nlev,irun)

c-----Solves tridiagpnaal system for the v-component using
c     the coefficients already computed in vtri1 (a,b,c)

      real a(irun,nlev),b(irun,nlev), c(irun,nlev)
      real F(irun,nlev),y(irun,nlev)
      real sfcflx(irun)
C
      DO 100 L = 2,nlev
       DO 9000 I = 1,irun
        F(I,L) = F(I,L) + F(I,L-1) * C(I,L)
 9000  CONTINUE
 100  CONTINUE
c-----Bottom level: includes surface fluxes
      do i = 1,irun
         y(i,nlev) = (f(i,nlev)+sfcflx(i))*a(i,nlev)
      enddo
C
c----performs back substitution
      DO 200 L = nlev-1,1,-1
       DO 9002 I = 1,irun
        Y(I,L) = (F(I,L) + B(I,L) * Y(I,L+1)) *  A(I,L)
 9002  CONTINUE
 200  CONTINUE
C
      RETURN
      END
c
      subroutine moist(th,sh,adz2,pl,plk,plke,khsave,xl,
     +                 dthv,irun,nlev,cp,rgas,virtcon,latvap,epsilo)

c-------Heat generated by condensation generates turbulence.The effect
c-------of condensation is included by recomputing the virtual
c-------potential temperature gradient, dthv.
c-------This subroutine calculates the gradient of the virtual temperature
c-------for the moist-turbulence scheme. In this calculation thv
c-------is defined as:thv=th*(1+0.61sh-1.61shl), where sh is the specific
c-------humidity of total water, and shl the one for liquid water.
c-------Equations used here are derived from : 'A Numerical simulation 
c-------of BOMEX data using a turbulence closure model coupled with ensemble
c-------cloud relations, Yamada and Mellor, Quart.J.Met.Soc, 1979,105.
c-------Calculations are performed only for levels above the surface
c       layer, i.e. surface layer is considered to be a dry layer.

#include<implicit.h>

        real b2,z1pem25,argmax,pi,rsq2pi,p5sr

        parameter(b2=10.1)
        parameter(z1pem25=1.e-25)
        parameter(argmax=30)
c-------erf function----------------------------------------------
#if (!defined LINUX)
        external function errf
#endif
        real errf

c-------input
        real th(irun,nlev)         !  model potential temperature
        real sh(irun,nlev)         !  specific humidity
        real adz2(irun,nlev-1)     !  distance between midpoints
        real pl(irun,nlev)         ! pressure at midpoints
        real plk(irun,nlev)        !  (pressure**k) at mid points 
        real plke(irun,nlev+1)     !  (pressure**k) at edge
        real khsave(irun,nlev-1)   ! non dimensional kh(kh/l*q)
        real xl(irun,nlev-1)       ! turbulent length scale
        integer irun,nlev,i,l
        real rgas,cp,virtcon,epsilo,latvap,vrt1con


c-------output
        real dthv(irun,nlev)    ! gradient for virtual potential temperature

c------local working space
       real dum
       real faceps,clh            ! constants
       real dth(irun,nlev-1)      ! potential temperature gradient
       real dsh(irun,nlev-1)      ! specific humidity gradient
       real tl(irun,nlev-1)       ! liquid water temperature
       real shsat(irun,nlev-1)    ! saturation specific humidity
       real qslt(irun,nlev-1)     ! (dshsat/dt)at liquid water temp tl(eq.26)
       real aa(irun,nlev-1)       ! eq 25.a
       real bb(irun,nlev-1)       ! b as in eq 25.b / (T/TH)
       real adq(irun,nlev-1)      ! adq=a*(qw-shsat)(eq 25 c * a)
       real ssdev(irun,nlev-1)    ! standard deviation for joint Gaussian R
       real svar(irun,nlev-1)     ! sqrt(ssdev)
       real fcc(irun,nlev-1)      ! R: cloud volume(eq. 30)
       real shl(irun,nlev-1)      ! mean liquid water specific humidity(eq.31)
       real q1m(irun,nlev-1)      ! adq/svar(eq.32)
       real temp(irun,nlev-1)     ! temporany space
       real funct(irun,nlev-1)    ! temporany space
       real arg(irun,nlev-1)      ! 
       real betat(irun,nlev-1)    ! parameter for calc. of dthv (eq.40a)
       real betaw(irun,nlev-1)    !     "                       (eq.40c)
       real betal(irun,nlev-1)    !     "                       (eq.40b)
       real betat1(irun,nlev-1)   !     "                       (eq.42a)
       real betaw1(irun,nlev-1)   !     "                       (eq.42b)
       

c------constants computation
       faceps=epsilo*latvap/rgas
       clh=latvap/cp
       pi = 4.*atan(1.)
       rsq2pi = 1./((2.*pi)**0.5)
       p5sr = 0.5**0.5
       vrt1con = 1.+virtcon

c------Gradients of potential temperature and  specifiy humidity 
c------Liquid water temperature

       do l=1,nlev-1
       do i=1,irun
          dth(i,l)=(th(i,l)-th(i,l+1))*adz2(i,l)
          dsh(i,l)=(sh(i,l)-sh(i,l+1))*adz2(i,l)
          tl(i,l)=th(i,l)*plk(i,l)
       enddo
       enddo
  
c------calculates shsat(Q_sl:mean saturation mixing ration of water vapour
c      at liquid water temperature tl)

       call aqsat(tl,pl,temp,shsat,irun,irun,nlev-1,1,nlev-1)


c------mean cloud fraction and mean liquid water specific humidity:
c      Computes a,b,adq,qslt,svar to be used in the expression for
c      fcc(mean cloud fraction) and shl(mean liquid water specific
c      humidity. 
c
       do l=1,nlev-1
       do i=1,irun
          qslt(i,l)=faceps*shsat(i,l)/(tl(i,l)*tl(i,l))
          aa(i,l)=1./(1.+clh*qslt(i,l))
          bb(i,l)=qslt(i,l)*aa(i,l)
          adq(i,l)=aa(i,l)*(sh(i,l)-shsat(i,l))
       enddo
       enddo
c
       do i = 1, irun
       ssdev(i,1) = xl(i,1)*(aa(i,1)*dsh(i,1)- 
     +              bb(i,1)*plke(i,2)*dth(i,1))
       ssdev(i,1) = b2*khsave(i,1)*ssdev(i,1)*ssdev(i,1)
       svar(i,1) = sqrt(ssdev(i,1))
       if ( svar(i,1).lt.z1pem25) svar(i,1) = z1pem25
       enddo

       do l=2,nlev-1
       do i = 1,irun
           ssdev(i,l) = xl(i,l-1)*(aa(i,l)*dsh(i,l-1)-bb(i,l)*
     +     plke(i,l)*dth(i,l-1))
           ssdev(i,l) = b2*khsave(i,l-1)*ssdev(i,l-1)*ssdev(i,l-1)
           svar(i,l) = sqrt(ssdev(i,l))
           ssdev(i,l) = xl(i,l)*(aa(i,l)*dsh(i,l)-bb(i,l)*
     +     plke(i,l+1)*dth(i,l))
           ssdev(i,l) = b2*khsave(i,l)*ssdev(i,l)*ssdev(i,l)
           temp(i,l)=sqrt(ssdev(i,l))
           svar(i,l)=(1./2.)*(svar(i,l)+temp(i,l))
           if ( svar(i,l).lt.z1pem25) svar(i,l) = z1pem25
       enddo
       enddo

       do l=1,nlev-1
       do i=1,irun
           q1m(i,l)= adq(i,l)/svar(i,l)
           funct(i,l)=errf(p5sr*q1m(i,l))
           fcc(i,l)=(1./2.)*(1.+funct(i,l))
           shl(i,l)=fcc(i,l)*adq(i,l)
           arg(i,l)=(1./2.)*q1m(i,l)*q1m(i,l)
       if(arg(i,l).le.argmax) then
            shl(i,l)=shl(i,l)+rsq2pi*svar(i,l)*exp(-arg(i,l))
       endif
c------Once the cloud fraction and the mean liquid water specific humidity
c      are computed, virtual potential temperature gradient
c      can be computed
c
       betat(i,l) = 1.+virtcon*sh(i,l)-vrt1con*shl(i,l)
       betaw(i,l) = virtcon*(th(i,l)+(clh/plk(i,l))*shl(i,l))
       betal(i,l) = (1. + virtcon*sh(i,l)-2.0*vrt1con*shl(i,l))*
     1              (clh/plke(i,l+1)) - vrt1con*th(i,l)

        betat1(i,l) = betat(i,l) - bb(i,l)*plk(i,l)*fcc(i,l)*betal(i,l)
        betaw1(i,l)  = betaw(i,l) + aa(i,l)*fcc(i,l)*betal(i,l)
        
        enddo
        enddo
   
c-------Compute dthv (edge) by averaging betat1 and betaw1 on the
c       adiacent midlevels. This is done for nlev -2 levels.
        do l=1,nlev-2
        do i =1,irun
            dthv(i,l)  = (1./2.)*((betat1(i,l)+betat1(i,l+1))*dth(i,l)+
     +                 (betaw1(i,l)+betaw1(i,l+1))*dsh(i,l))
        enddo
        enddo

c-------Compute dthv(nlev-1) skipping the averaging as we
c       don't have beta1(nlev)

        do i = 1,irun
             dthv(i,nlev-1) = betat1(i,nlev-1)*dth(i,nlev-1) +
     +                        betaw1(i,nlev-1)*dsh(i,nlev-1)
        enddo

        return

        end

        real function errf(arg)
C***********************************************************************
C  FUNCTION ERRF
C  PURPOSE
C     COMPUTES ERROR FUNCTION OF ARGUMENT
C  USAGE
C     CALLED BY TRBFLX
C  DESCRIPTION OF PARAMETERS
C     ARG   -  INPUTED ARGUMENT
C  REMARKS:
C        USED TO COMPUTE FRACTIONAL CLOUD COVER AND LIQUID WATER CONTENT
C          FROM TURBULENCE STATISTICS
C **********************************************************************

      PARAMETER ( AA1  =   0.254829592 )
      PARAMETER ( AA2  =  -0.284496736 )
      PARAMETER ( AA3  =   1.421413741 )
      PARAMETER ( AA4  =  -1.453152027 )
      PARAMETER ( AA5  =   1.061405429 )
      PARAMETER ( PP   =   0.3275911   )
      PARAMETER ( X2   =   AA2 / AA1   )
      PARAMETER ( X3   =   AA3 / AA2   )
      PARAMETER ( X4   =   AA5 / AA3   )
      PARAMETER ( X5   =   AA5 / AA4   )

      ERRF = 1.
      AARG=ABS(ARG)

      IF ( AARG .LT. 4.0 )  THEN
       TT = 1./(1.+PP*AARG)
       ERRF = 1. -
     1  (AA1*TT*(1.+X2*TT*(1.+X3*TT*(1.+X4*TT*(1.+X5*TT)))))
     2   * EXP(-AARG*AARG)
      ENDIF

      IF ( ARG .LT. 0.0 )  ERRF = -ERRF

      RETURN
      END
