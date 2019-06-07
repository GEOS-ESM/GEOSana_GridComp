      SUBROUTINE TRBFLX (tnlev,TH,THV,SH,U,V,QQ,PL,PLE,PLK,plke,DELTAP,
     &    taux,tauy,shflx,cflx,DTAU,ITRTRB,KMBG,KHBG, qbeg,tprof,
     &    jlat,irun,nlev,pcnst,grav,cp,rgas,cappa,
     &    WU,WV,WT,WSH,et,eu,pbldpth,trbfast)



c--------------------------------------------------------------------
c  Subroutine trbflx-computes turbulent adjustments to atmospheric
c                    profiles.
c                   - called from turb1
c
c  needed subroutines for trbflx:
c
c  ---->  initurb :       initialize variables for calculations
c  ---->  trblen  :       computes turbulent length scale
c       -----> trbitp:    interpolation needed in trblen
c  ---->  trbl20  :       Mellor-Yamda level 2 model
c  ---->  trbgrow :       Dimensionless diffusion coeff. and
c                         production term for growing turbulence
c  ---->  trbl25  :       Dimensionless diffusion coeff. and production
c                         term for decaying turbulence 
c  ---->  trbdifq :       solves diffusion equation for tke
c  ---->  trbdif1 :       solves diffusion equation for th and sh
c  ---->  trbdif2 :       solves diffusion equation for u and v
c---------------------------------------------------------------------
CC    version vc7.0: turbulence, no lsm, no moist turbulence
C     February/March 01 : modification of the trbflx to use it
C     in the fvcc3 model. 
C     The lsm will supply surface fluxes as input the turbulence
C     package.
C     Trbflx updates TKE, calculates the diffusion coefficients
C     and thereafter updates u,v,T and sh
c     Internal timestep 
C
#include <implicit.h>

      real B1, B3, ALPHA, HALPHA, QQMIN, QBUSTR, AITR
      PARAMETER ( B1      =   16.6    )  
      PARAMETER ( B3      = 1. / B1  )  
      PARAMETER ( ALPHA   = 0.1       )
      PARAMETER ( HALPHA = ALPHA * 0.5 )
      PARAMETER ( QQMIN  = 0.005      )
      PARAMETER ( QBUSTR = 2.550952   )

c-----------------------------------------------------------------------
c-----input

      real KMBG,KHBG         ! background diffusion coeff
      LOGICAL qbeg           ! qbeg=true for first run of the model
      logical tprof          ! =true to calculate diagnostics
      integer init            ! init=0 for initial start
      real rgas,cp,virtcon,grav,rcp
      real cappa
      integer jlat,m
      integer irun,nlev,itrtrb
      integer pcnst          ! number of constituent
      real dtau              !time step
      real sumz,temppp       !temporany storage for pbldpth calc.
      real akap              ! rgas/cp
      real pzero             !p_0**akap

      real tnlev(irun)  ! temperature at nlev
      real pl(irun,nlev)   ! pressure at mid-points
      real ple(irun,nlev+1)  ! pressure at edge levels
      real plk(irun,nlev)  ! (pressure**k)  at mid-points
      real plke(irun,nlev+1) ! (pressure**k)  at edge
      real deltap(irun,nlev)    ! pressure intervals (edge)
      real rpdele(irun,nlev)    ! reciprocal of pressure (edge)

      real taux(irun)        ! surface ustress
      real tauy(irun)        ! surface vstress
      real shflx(irun)       ! surface heat flux
      real cflx(irun,pcnst)        ! surface constituent flux


      

c-----------------------------------------------------------------------
c-----input/output
      real th(irun,nlev)   ! model potential temperature(theta/po**k)
      real thv(irun,nlev)  ! model virtual potential temperature
      real sh(irun,nlev)   ! specific humidity
      real u(irun,nlev)    ! u-wind
      real v(irun,nlev)    ! v-wind
      real qq(irun,nlev)   ! TKE, defined on edges


C----------------------------------------------------------------------
c-------------------------local work space-----------------------------
      integer istnm1, nlevm1,l,i,ll,lp1,iter,nlevml
      integer lmin, lminq, lminq1,nlevmq,lp
      real gbycp, fac1, fac2,vk
      real hs(irun)          ! depth of the surface layer

c-----for pbl height
      real qqcolmin(irun) ! variable for calculation of pbl depth
      real qqcolmax(irun) ! variable for calculation of pbl depth
      real levpbl(irun)   ! variable for calculation of pbl depth

c-----vertical grid
      real adz1(irun,nlev)   ! distance between edges
      real adz2(irun,nlev-1) ! distance between mid-points
      real dz3(irun,nlev)    ! alpha*adz1
      real vkze(irun,nlev-1) ! k*height on edge level
      real vkzm(irun,nlev-1) ! k*height on mid level

c-----gradients, shears, B-V frquence and Ri number
      real du(irun,nlev)   ! u-gradient
      real dv(irun,nlev)   ! v-gradient
      real dthv(irun,nlev) ! virtual potential temperature gradient
      real strt(irun,nlev) ! Brunt-vaisal frequence
      real dw2(irun,nlev)  ! shears:du**2 + dv**2
      real ri(irun,nlev)   ! Richardson number

c-----variables for the diffusion equation
      real rhodz1(irun,nlev-1)    !rho/dz at midpoints
      real rhodz2(irun,nlev-1)    !rho/dz at edges
      real rhozpk(irun,nlev-1)    !rhodz2*plke 
      real rhokdz(irun,nlev-1)    !rhodz(1,2)*km(kh)
      real flxfac(irun,nlev)      !gdt/dpe for the u,v,sh diffusion eq. 
      real flxfpk(irun,nlev)      !gdtpk/dpe for the th diffusion eq. 
      real flxfce(irun,nlev-1)    !gdt/dp for the qq diffusion eq. 
      real dpk(irun,nlev)         ! gradients of p**k on midlevels

      real km(irun,nlev) ! diff coe for momentum inside the internal timestep
      real kh(irun,nlev) ! difn coeff for heat  inside the internal time step

      real tmp1(irun)      ! local space
      real dusurf(irun)    ! lowest layer u change from stress
      real dvsurf(irun)    ! lowest layer v change from stress
      real dthsurf(irun)   ! lowest layer model th  change from heat flux
      real dqsurf(irun)    ! lowest layer sh  change from constituent flux

c-----TKE, turbulent length scale, production term for solution of
c     tke prognostic equation
      real q(irun,nlev)      ! turbulent kinetic velocity
      real qqe(irun,nlev-1)  ! TKE from level 2 model
      real qe(irun,nlev-1)   ! TK velocity  from level 2 model
      real p3(irun,nlev-1)   ! production term for tke prognostic equation
      real xq(irun,nlev-1)   ! tke velocity / turbulent length scale
      real xl(irun,nlev-1)   ! turbulent length scale
      real qxlm(irun,nlev-1) ! (5/3*a*xl)*q(triple coerr term in tke eq. 

      real temp(irun,nlev)
 
      integer ibitstb(irun,nlev) ! =1 for growing turbulence, 0 for dec.
      integer intq(irun,nlev)    ! =1 for points with qq.e-8

c-----Variables used in trblen
      real DZITRP(irun,nlev-1), STBFCN(irun,nlev)
      real    XL0(irun,nlev),       Q1(irun,nlev-1)
      real WRKIT1(irun,nlev-1)
      real WRKIT2(irun,nlev-1)
      real WRKIT3(irun,nlev-1)
      real WRKIT4(irun,nlev-1)
      INTEGER INT1(irun,nlev), INT2(irun,nlev-1)
c------------------------------------------------------------------------------
c-----output
      real pbldpth(irun) ! PBL depth(in meters)
      real pbldpthpa(irun) ! PBL depth(in pa)

      real wt(irun,nlev-1)  !heat flux
      real wsh(irun,nlev-1) !flux for sp humidity
      real wu(irun,nlev-1)  !flux for u-wind
      real wv(irun,nlev-1)  !flux for v-wind
C----------------------------------------------------------------
c-----Added variables for using an internal timestep
      real et(irun,nlev)    !final diffusion coeff
      real eu(irun,nlev)    !final diffusion coeff
      real trbfast(irun,nlev)   !flag for growing to fast cond.
      real check(irun,nlev)    


      LOGICAL FIRST,LAST         !first=true for first iteration of trbflx
!                                 (in this case first is always true as
!                                  we are not going to have internal timestep


      AITR   = 1. / FLOAT(ITRTRB)
      NLEVM1 = NLEV - 1
      ISTNM1 = irun * NLEVM1
      GBYCP  = GRAV / CP
      vk = 0.4
      virtcon=0.61
c-------------------------------------------------------------------------
c-----determine if initial start
      INIT = 0
      IF(qbeg) INIT = 1
c-------------------------------------------------------------------------
c-----if diagnostics are needed, initialize diagnostics variables
      do l=1,nlev-1
          do i=1,irun
             wu(i,l)=0.
             wv(i,l)=0.
             wt(i,l)=0.
             wsh(i,l)=0.
           enddo
      enddo
      do l=1,nlev
          do i=1,irun
             eu(i,l)=0.
             et(i,l)=0.
             trbfast(i,l)=0.
           enddo
      enddo
           do i = 1,irun
              pbldpthpa(i)=0.
              pbldpth(i)=0.
           enddo

C-----Starts calculation of  several variables to be used after
C
c-----Calculates p_0**k, pzero
      pzero = (100000.)**cappa
c-----Compute vertical grid
c-----Distance between edges
c
c*************************************************************************
      do  l=1,nlev
        do  i=1,irun
          adz1(i,l)=(cp/grav)*thv(i,l)*(plke(i,l+1)-plke(i,l))
        enddo
      enddo
c
c-----Distance between midpoints
c
      do  l=1,nlev-1
        do  i=1,irun
          adz2(i,l)=0.5*(adz1(i,l)+adz1(i,l+1))
        enddo
      enddo

c-----depth of surface layer
      DO 9042 I =1,irun
      hs(I) = 0.5 * ADZ1(I,NLEV)
9042  CONTINUE
c
c-----alpha*layer depths for trblen
c
      do  l=2,nlev-1
        do i=1,irun
           dz3(i,l)=alpha*adz1(i,l)
        enddo
      enddo
      do i=1,irun
        dz3(i,1)=halpha*adz1(i,1)
        dz3(i,nlev)=alpha*hs(i)
      enddo
c
c-----vk*heights at edge levels (vkze) and at mid levels(vkzm)
c     vkze(1) is defined on the first edge-level below the first midlevel
c     (from the top), vkzm(1) is defined on the second mid-level
c
      do l=2,nlev
         do i=1,irun
            temp(i,l)=vk*adz1(i,l)
         enddo
      enddo

      do  i =1,irun
      vkze(i,nlevm1)=temp(i,nlev)
      enddo

      do 100 ll = 2,nlevm1
      l = nlev-ll
      lp1 = l + 1
        do i =1,irun
           vkze(i,l) = vkze(i,lp1) + temp(i,lp1)
        enddo
 100  CONTINUE
c     
      do l=1,nlev-1
         do i=1,irun
            vkzm(i,l)=vkze(i,l)-0.5*temp(i,l+1)
         enddo
      enddo
C --------------------------------------------------------------------
c    compute rho/dz to be used later in the diffusion equation
c    rhodz2: rho/dz at edges.rhodz2(1) is defined on the first
c            edge level below the first mid level(from top) 
c    rhodz1: rho/dz at midpoints.rhodz1(1) is defined on the second mid-level
c    rhozpk: rhodz2*plke: includes additional term plke for the
c            diffusion equation of temperature


c-----computes thv on edges by weighting the midlevels values
       do  l = 1,nlevm1
         lp1 = l + 1
            do i =1,irun
               fac1=deltap(i,l)/(deltap(i,l)+deltap(i,lp1))
               fac2=1.-fac1
               rhodz2(i,l)=fac1*thv(i,lp1)+fac2*thv(i,l)
            enddo
       enddo


      do l=1,nlev-1
         do i=1,irun
            rhodz2(i,l)=(rgas*0.01)*rhodz2(i,l)*plke(i,l+1)*adz2(i,l)
            rhodz2(i,l)=ple(i,l+1)/rhodz2(i,l)
            rhozpk(i,l)=plke(i,l+1)*rhodz2(i,l)
            rhodz1(i,l)=(rgas*0.01)*thv(i,l+1)*plk(i,l+1)*adz1(i,l+1)
            rhodz1(i,l)=pl(i,l+1)/rhodz1(i,l)
         enddo
      enddo
c
c-------------------------------------------------------------------------
c     compute the other term in the diffusion equation,i.e.
c     flxfac=gdt/dpe for the u,v and sh diffusion equation
c            flxfac(1) is defined on the first mid level(as p(1))
c     flxfpk=gdt(pk)/dpe for the temperature diffusion equation
c            flxfpk(1) is defined on the first mid level(as p(1))
c     flxfce=gdt/dp for the tke equation
c            flxfce(1)is defined on the first edge level below the first
c            midlevel(as for rhodz1,vkzm)

      do l=1,nlev
         do i=1,irun
            flxfpk(i,l)=plk(i,l)*(ple(i,l+1)-ple(i,l))
c           flxfpk(i,l)=plk(i,l)*deltap(i,l)
            flxfpk(i,l)=(grav*dtau*0.01)/flxfpk(i,l)
         enddo
      enddo

      do l=1,nlev
         do i=1,irun
            flxfac(i,l)=flxfpk(i,l)*plk(i,l)
         enddo
      enddo


      do l=1,nlev-1
         do i=1,irun
            flxfce(i,l)=(grav*dtau*0.01)/(pl(i,l+1)-pl(i,l))
         enddo
      enddo

c-------------------------------------------------------------------------
c-------Compute reciprocals of dz1,dz2,hs
c      do l=1,nlev
c         do i=1,irun
c            adz1(i,l)=1./adz1(i,l)
c         enddo
c      enddo
      do l=1,nlev-1
         do i=1,irun
            adz2(i,l)=1./adz2(i,l)
         enddo
      enddo

c-----compute gradients of p**kappa
      do l=1,nlev-1
         do i=1,irun
            dpk(i,l)=(plk(i,l+1)-plk(i,l))*adz2(i,l)
         enddo
      enddo
         do i=1,irun
            dpk(i,nlev)=gbycp/thv(i,nlev)
         enddo
c------compute reciprocal of pressure depth on layers
       do l=1,nlev
         do i=1,irun
          rpdele(i,l) = 1./deltap(i,l)
         enddo
       enddo

c---------------------------------------------------------------------
c-----initialize Q array
      do l=1,nlev-1
         do i=1,irun
            q(i,l)=sqrt(2.*qq(i,l))
         enddo
      enddo

      FIRST = .TRUE.
      LAST = .FALSE.
c---------------------------------------------------------------------
C                              MAIN LOOP
C
c*********************************************************************
      DO 2000 ITER = 1, ITRTRB
C
      IF ( ITER .GE. ITRTRB ) LAST = .TRUE.
C
c-----calculates Brunt-Vaisala frequencies and Shears for suface layer
c     calculates qq(nlev) and q(nlev)
      call surfvar(jlat,taux,tauy,shflx,cflx(1,1),hs,vk,virtcon,tnlev,
     +             u(1,nlev),v(1,nlev), thv(1,nlev),
     +             th(1,nlev),sh(1,nlev),plk(1,nlev),pl(1,nlev),
     +             rgas,grav,cp,
     +             pzero,qbustr, q(1,nlev),qq(1,nlev),
     +             du(1,nlev),dv(1,nlev),dthv(1,nlev),
     +             strt(1,nlev),dw2(1,nlev),ri(1,nlev),irun)

c-----gradients
      do l=1,nlev-1
         do i=1,irun
            du(i,l)=(u(i,l)-u(i,l+1))*adz2(i,l)
            dv(i,l)=(v(i,l)-v(i,l+1))*adz2(i,l)
            dthv(i,l)=(thv(i,l)-thv(i,l+1))*adz2(i,l)
         enddo
      enddo

c------Brunt-Vaisala frequencies,shears,Ri numbers
      do l=1,nlev-1
         do i=1,irun
           strt(i,l)=cp*dthv(i,l)*dpk(i,l)
           dw2(i,l)=du(i,l)*du(i,l)+dv(i,l)*dv(i,l)
           if(dw2(i,l).le.1.e-4) dw2(i,l)=1.e-4
           ri(i,l)=strt(i,l)/dw2(i,l)
         enddo
      enddo
c------Added for checking if turbulence grows too fast
      do l=1,nlev
       do i=1,irun
          check(i,l)=0.0
       enddo
      enddo
c---------------------------------------------------------
c-----Call the subroutine to initialize the calculations

      call initurb(kh,km,qqe,qe,p3,qq,q,ibitstb,intq,lminq,lminq1,
     +irun,nlev)

c ----Determine the length scale      
       CALL TRBLEN(STRT,DW2,DZ3,Q,VKZE,VKZM,DTHV,DPK,DU,DV,XL,QXLM,
     1  NLEV,INIT,lmin,LMINQ,LMINQ1,CP,INT1,INT2,
     2  DZITRP,STBFCN,XL0,Q1,WRKIT1,WRKIT2,WRKIT3,WRKIT4,irun)
c---------------------------------------------------------------------------
c ---Determine QE and dimensionless coeff from level 2 model
c    This is done only if instability was found above the 
c    surface layer (nlev)

c------Computes first nlevml=nlev-lmin :number of level to be processed,
c      which are the levels where instability has been found.
c      NO calculations of qe or diffusion coefficients are done on 
c      the levels where no instability has been found.
c      Therefore these calucaltions are skipped if the first
c      level where instabilty occurs(lmin) is below the last model level(nlev)

      if ( lmin .lt. nlev) then
          nlevml = nlev-lmin

          CALL TRBL20(RI(1,LMIN),STRT(1,LMIN),DW2(1,LMIN),XL(1,LMIN),
     1     KM(1,LMIN),KH(1,LMIN),QE(1,LMIN),QQE(1,LMIN),IBITSTB(1,LMIN),
     2     NLEVML,nlev,irun)
       endif
c-------------------------------------------------------------------------------
c-----for initial start only: use equilibrium model(50)
      IF ( INIT .EQ. 1 ) THEN
         do l=1,nlev-1
            do i=1,irun
               qq(i,l)=qqe(i,l)
               q(i,l)=qe(i,l)
            enddo
         enddo
        INIT = 2
        CALL TRBLEN(STRT,DW2,DZ3,Q,VKZE,VKZM,DTHV,DPK,DU,DV,XL,QXLM,
     1  NLEV,INIT,LMIN,LMINQ,LMINQ1,CP,INT1,INT2,
     2  DZITRP,STBFCN,XL0,Q1,WRKIT1,WRKIT2,WRKIT3,WRKIT4,irun)

       INIT = 0
       GO TO 550
      ENDIF

c------------------------------------------------------------------------------
c-----calculations of dimensionless diffusion coefficients (kh,km)
c     and of  a part of the production term P3 for *growing* 
c     turbulence(q le qe)
c     Calculations are performed only for levels where instability was
c     found

      if(lmin.lt.nlev) then
          call trbgrow(irun,nlevml,ibitstb(1,lmin),q(1,lmin),
     +                 qe(1,lmin),qq(1,lmin),qqe(1,lmin),
     +                 kh(1,lmin),km(1,lmin),p3(1,lmin))
       endif

              
c------------------------------------------------------------------------------
c-----calculations of dimensionless diffusion ocefficients (kh,km)
c     and of part of the production term p3 for *decaying* 
c     turbulence(q gt qe)
c     Calculations are performed only for levels where
c     TKE is different form zero(qq gt e-8) and where q gt qe (ibitstb=0)

       NLEVML = NLEV - LMINQ
       CALL TRBL25(Q(1,LMINQ),XL(1,LMINQ),STRT(1,LMINQ),DW2(1,LMINQ),
     1    IBITSTB(1,LMINQ),INTQ(1,LMINQ),KM(1,LMINQ),KH(1,LMINQ),
     2    P3(1,LMINQ),NLEVML,nlev,irun)

c-----------------------------------------------------------------------------
c
c-----calculate source term P3
c-----If turbulence occurs at a higher level of where instability occurs
c     (i.e. the first level where qq gt e-8 (lminq) is  closer to the top
c     of the first level where instability occurs (lmin)),
c     then lminis set equal to lminq.
      
      IF ( LMINQ .LT. LMIN ) THEN
       LMIN = LMINQ
      ENDIF

c-----The final form of the production term is calculated for both
c     growing and decaying turbulence. 

      IF( LMIN .LT. NLEV ) THEN

c       nlmn=nlev-lmin 
c       do l=lmin,lmin+(nlmn-1)
        do l=lmin,nlev-1
          do i=1,irun
             p3(i,l)=p3(i,l)*dtau/xl(i,l)
          enddo
       enddo
c-----If turbulence grows too fast, allowing q from q<qe
c     to get q>qe in only one time step, so going
c     from a condition of growing turbulence to a condition
c     of decaying turbulence in just one timr step, or
c     if the opposite occurs, i.e. q>qe gets q<qe in one
c     time step, turbulence decays too fast, we prevent this
c     to happen by redefing an appropriate production term
c
c       do l=lmin,lmin+(nlmn-1)
        do l=lmin,nlev-1
          do i=1,irun
           temp(i,l)=qqe(i,l)*p3(i,l)
           xq(i,l)=qqe(i,l)-qq(i,l)
              if ( ((ibitstb(i,l).eq.1).and.(xq(i,l).lt.temp(i,l)) )
     +                   .or.
     +        ((ibitstb(i,l).eq.0).and.(xq(i,l).gt.temp(i,l)))) then
              p3(i,l)=xq(i,l)/qqe(i,l)
              check(i,l)=1. 
              endif 
          enddo
       enddo
     
      ENDIF
c-----------------------------------------------------------------------

 550  CONTINUE
c
c--------Prepares to solve the TKE equation
c        The diffusion equation to be solved is:
c        D(QQ)=d/dz(aq d/dz(qq))dt +p3*qq
c        where aq=qxlm was computed in trblen
c        
      nlevmq=nlev-lminq1

      do l=lminq1,nlev-1
          do i=1,irun
             rhokdz(i,l)= rhodz1(i,l)*qxlm(i,l)
          enddo
      enddo

       call trbdifq (qq(1,lminq1),p3(1,lminq1),rhokdz(1,lminq1),
     +               flxfce(1,lminq1),nlevmq,irun)

c-------------------------------/---------------------------------------
c
C   DIMENSIONAL DIFFUSION COEFFS INCLUDING BACKGROUND AMOUNTS;
c   Kh and km are defined on edges, kh(irun,nlev)  defines
c   the turbulent coefficient at the ground, (which in this case
c   is not defined as calucalations in the sl
c   are performed in the LSM), kh(irun,1) defines
c   it at the first edge level below the first mid-level(from the top)
c   In this version, surface calculations are performed by the LSM, therefore
c   kh,km(irun,nlev)=0.
c   For the levels going from 1 to lminq1, where therefore there is
c   no turbulence (lminq1 is the first level (from the top) where qq gt e-8
c   plus 1), the turbulence coefficients Kh and Km are set equal to
c   a background amount (khbg,kmbg),which for this version is zero.
c   For the levels from lminq1 to nlev-1, kh are expressed by
c   kh=q*l*Sh for decaying turbulence, kh=(q/qe)*q*l*Sh for growing
c   turbulence (same expression for km with Sm). The dimensionless 
c   turbulent coefficients (Sh for decaying and (q/qe)Sh for
c   growing turbulence) have been already computed, so the final
c   diffusion coefficients(dimensional) to be used for the
c   diffusion equation for u,v,T sh are simple computed
c   by multiplying the dimensionless coeff. by q and l.

c-----Kh and Km for levels with no turbulence

      IF(LMINQ1.GT.1)THEN
       do l=1,lminq1-1
          do i=1,irun
             km(i,l)=kmbg
             kh(i,l)=khbg
          enddo
       enddo
      ENDIF
C
c-----Dimensional diffusion coefficients for levels from lminq1-1 to nlev-1

c      do l=lminq1,lminq1+(nlevmq-1)
       do l=lminq1,nlev-1
         do i=1,irun
            q(i,l)=sqrt(2.*qq(i,l))
            km(i,l)=xl(i,l)*q(i,l)*km(i,l)+kmbg
            kh(i,l)=xl(i,l)*q(i,l)*kh(i,l)+khbg
         enddo
      enddo

c-----convert the surface fluxes in tendencies on the lowest level
      
       rcp = 1./cp
       do i = 1,irun
       tmp1(i) = dtau*grav*rpdele(i,nlev)
       dusurf(i) = taux(i)*tmp1(i)
       dvsurf(i) = tauy(i)*tmp1(i)
c       dthsurf(i) = (shflx(i)*tmp1(i)*rcp)*plke(i,nlev+1)
c       dthsurf(i) = (shflx(i)*tmp1(i)*rcp)/pzero
       dthsurf(i) = (shflx(i)*tmp1(i)*rcp)/plk(i,nlev)
       dqsurf(i) = cflx(i,1)*tmp1(i)
       enddo

c---------------------------------------------------------------------------
C   CALCULATE INTERNAL FLUXES AND UPDATE PROGNOSTIC VARIABLES: TH AND S
c
c-----prepares coefficients for diffusion equation for th
      do l=1,nlev-1
         do i=1,irun
            temp(i,l)=rhozpk(i,l)*kh(i,l)
         enddo
      enddo
      

       call trbdif1(th,temp,flxfpk,dthsurf,nlev,irun)


c-----preprares coefficient for diffusion equation for sh
      do l=1,nlev-1
         do i=1,irun
            rhokdz(i,l)=rhodz2(i,l)*kh(i,l)
         enddo
      enddo

       call trbdif1(sh,rhokdz,flxfac,dqsurf,nlev,irun)
C
c-----if we are using the internal time step, we need to know the
c     virtual potential temperature obtained with the updated
c     potential temperature and specific humidity

      IF(.NOT.LAST) THEN

      do l=1,nlev
         do i=1,irun
            thv(i,l) = th(i,l)*(1.+virtcon*sh(i,l))
         enddo
      enddo

      ENDIF
C
c---------------------------------------------------------------------------
C   CALCULATE INTERNAL FLUXES AND UPDATE PROGNOSTIC VARIABLES: U AND V
C
      do l=1,nlev-1
         do i=1,irun
            rhokdz(i,l)=rhodz2(i,l)*km(i,l)
         enddo
      enddo

       call trbdif2(u,v,rhokdz,flxfac,dusurf,dvsurf,nlev,irun)
c---------------------------------------------------------------------------
c      In the following part the pbl height is computed
c      We don/t need that if pbldif_simpl is called in tphysac
c      after turb, as pbl height is then computed in pbldif_simpl
c      together with tpert and qpert using original NCAR formulation.
c      In case we want pbl computed with MH scheme, we need
c      to 'uncomment' the following lines, and call
c      pbldif_sompl_tq  in tphysac
c
c---------------------------------------------------------------------------
C-----Calculates depth of PBL
c       do i = 1,irun
c        qqcolmin(i) = qq(i,nlev)*0.1
c        qqcolmax(i) = qq(i,nlev)
c        levpbl(i) = nlev
c       enddo
c       DO L = nlev-1,1,-1
c        DO I = 1,irun
c         IF ( (qq(i,l).gt.qqcolmax(I)).and.(levpbl(i).eq.nlev))then
c          qqcolmax(i) = qq(i,l)
c          qqcolmin(i) = 0.1*qqcolmax(I)
c         endif
c         if((qq(i,l).lt.qqcolmin(i)).and.(levpbl(i).eq.nlev))
c     1                                                levpbl(i)=L
c        enddo
c       enddo
cc**************************************************************************
cc-------We are trying to calculate the pbl in meters rather than in
cc        pascal. To be able to compare we are going to leave the 
cc        computation of the pbl i pascals, calling in pbldpthpa
cc*******old version***************************
cc       do i = 1,irun
cc        lp = levpbl(i)
cc        if(lp.lt.nlev)then
cc         pbldpth(I) = pbldpth(I) + ( (PLE(I,nlev+1)-PLE(I,Lp+2)) +
cc     1    ( (ple(i,lp+2)-ple(i,lp+1))*(qq(i,lp+1)-qqcolmin(i))
cc     2          / (qq(i,lp+1)-qq(i,lp)) ) ) * aitr
cc        else
cc         pbldpth(I) = pbldpth(I) + ( (PLE(I,nlev+1)-PLE(I,2)) +
cc     1    ( (ple(i,2)-ple(i,1))*(qq(i,1)-qqcolmin(i))
cc     2          / qq(i,1) ) ) * aitr
cc        endif
cc       enddo
cc*******old version with new name for pbldpth***************************
c       do i = 1,irun
c        lp = levpbl(i)
c        if(lp.lt.nlev)then
c         pbldpthpa(I) = pbldpthpa(I) + ( (PLE(I,nlev+1)-PLE(I,Lp+2)) +
c     1    ( (ple(i,lp+2)-ple(i,lp+1))*(qq(i,lp+1)-qqcolmin(i))
c     2          / (qq(i,lp+1)-qq(i,lp)) ) ) * aitr
c        else
c         pbldpthpa(I) = pbldpthpa(I) + ( (PLE(I,nlev+1)-PLE(I,2)) +
c     1    ( (ple(i,2)-ple(i,1))*(qq(i,1)-qqcolmin(i))
c     2          / qq(i,1) ) ) * aitr
c        endif
c       enddo
cc*******new  version with pbldpth in meters***************************
c       do i = 1,irun
c        lp = levpbl(i)
c        sumz=0.0
c        if(lp.lt.nlev)then
c
c              if (lp.eq.(nlev-1)) then
c                sumz=0.0
c                else
c                do l = nlev,lp+2,-1
c                sumz=sumz+adz1(i,l)
c                enddo
c              endif
c         pbldpth(I) = pbldpth(I) + (sumz +
c     1    ( adz1(i,lp+1)*(qq(i,lp+1)-qqcolmin(i))
c     2          / (qq(i,lp+1)-qq(i,lp)) ) ) * aitr
c
c        else
c
c        do l = nlev,2,-1
c        sumz = sumz + adz1(i,l)
c        enddo
c        pbldpth(I) = pbldpth(I) + ( sumz +
c     1    ( adz1(i,1)*(qq(i,1)-qqcolmin(i))
c     2          / qq(i,1) ) ) * aitr
c        endif
c       enddo
cC**************************************************************************
C      Diagnostics are saved here. This part of the code has to be
c      changed according to the kind of diagnostics required from
c      fvgcm 
c      NB: boundaries need to be changed, as variables are not defined
c     on nlev+1. This is just an example!!
c     NB: the diagnostic variables have to be initialized by 
c     setting them to zero out of the time loop!

c      if (tprof) then
c      do l=1,nlev-1
c         do i=1,irun
c            wt(i,l)=wt(i,l)+temp(i,l)*(th(i,l+1)-th(i,l))
c            wsh(i,l)=wsh(i,l)+rhokdz(i,l)*(sh(i,l+1)-sh(i,l))
c            wu(i,l)=wu(i,l)+rhokdz(i,l)*(u(i,l+1)-u(i,l))
c            wv(i,l)=wv(i,l)+rhokdz(i,l)*(v(i,l+1)-v(i,l))
c         enddo
c      enddo

c-----Diffusion coefficients
      do l=1,nlev
         do i=1,irun
            et(i,l)=et(i,l)+aitr*kh(i,l)
            eu(i,l)=eu(i,l)+aitr*km(i,l)
            trbfast(i,l)=trbfast(i,l)+aitr*check(i,l)
         enddo
      enddo
c      endif

c*********************************end of diagnostics**************************

      FIRST = .FALSE.
C
C   END OF MAIN LOOP
C
 2000 CONTINUE
C*****************************more diagnostics*****************************
c      if (tprof) then
c      do l=1,nlev-1
c         do i=1,irun
c            wu(i,l)=wu(i,l)*aitr
c            wv(i,l)=wv(i,l)*aitr
c            wt(i,l)=wt(i,l)*aitr
c            wsh(i,l)=wsh(i,l)*aitr
c         enddo
c      enddo
c      endif
C*****************************more diagnostics*****************************
C
C
      return
      END
