#include <params.h>
       subroutine turb(lat,hbuf,ztodt,
     .                 p,pe,pdelm1,um1,vm1,tm1,qm1,
     .                 taux,tauy,shflx,cflx,fu,fv,t2,
     .                 duv,dvv,dtv,up1,vp1, 
     .                 qbeg,prsk,prske,qq,
     .                 hattr,pblh) 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC      Version vc7.0, no lsm, no moisture turbulence.
C         Modification to use it with the fvcc3.f:
C         Change of vertical coordinate from sigma to p
C         Stripping of variables is done out of turb.f
C         So variables in turb.f are 2 dimensional
C         The main change consists in the different boundary conditions
c         i.e. in the different way the turbulnce code is linked
C         with the LSM. To be able to use the geos turbulnce
C         code majour modification have to be done:
C      1. remove of all the surface calucaltions. Surface
C         calculations are done by the LSM which outputs
C         the surface fluxes (momentum stress, sensible heat and
C         latent heat).
C      2. There are no values at the lower boundary, instead
C         we need to use the surface fluxes directly.
C      3. Includes the internal time stepping in trbflx. 
C         However the surface fluxes don't change over 30 min period

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         Turb.f is called where  vdintr.F is usually called, i.e.
C         directly fron tphysac.F

#include <implicit.h>


#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comtrcnm.h>
C-----------------------------------------------------------------------
#include <comcon.h>
C------------------------------Diagnostics-------------
#include <diag.h>
c------------------------------------------------------------------------------
       logical qbeg                      !true for cold start
       logical tprof                     !true if diagnostic needed
c********************check what to use here*********************************
c----input 
      real ztodt                          ! delta-t
      integer hattr(5,pdiag)              ! History fiels attribute
      real*4 hbuf(*)                      ! History buffer
      integer lat

      real    pe(plond,plevp)   ! pressure at edge level
      real    p(plond,plev)      ! pressure at even level
      real    prske(plond,plevp)  ! pressure at edge level ** k
      real    prsk(plond,plev)     ! pressure at even level **k
      real    pdelm1(plond,plev) ! press inter in the vertical(on edges)

      real    um1(plond,plev)     ! u-wind(m/s) input
      real    vm1(plond,plev)     ! v-wind(m/s) input
      real    tm1(plond,plev)     ! temperatura(K) input
      real    taux(plond)         ! x surface stress
      real    tauy(plond)         ! y surface stress
      real    shflx(plond)         ! surface sensible heat flux (w/m2)
      real    cflx(plond,pcnst)    ! surface constituent flux (kg/m2s)



C-----updated

      real  qq(plond,plev)  ! Turbulent kinteic energy
      real fu(plond,plev)   ! Cumulative u tendency
      real fv(plond,plev)   ! Cumulative v tendency
      real t2(plond,plev)   ! Cumulative T tendency
      real qm1(plond,plev,pcnst) ! init/fin moisture const 

C-----output
C     Turbulence Tendency: Change in u,v,t or q  due to turbulent
c      processes per unit time in m/sec**2
       real  duv(plond,plev)  ! u-wind Turbulence Tendency
       real  dvv(plond,plev)  ! v-wind Turbulence Tendency
       real  dtv(plond,plev)  ! temperature Turbulence Tendency
       real up1(plond,plev)      ! u-wind after vertical diffusion
       real vp1(plond,plev)      ! v-wind after vertical diffusion
c
c-----------------------------------------------------------------                                         
c-------------variable needed for version exp2, where pblh is
c             calculated according to MH scheme
       real pblh(plond)
c-----------------------------------------------------------------                                         
c             needed for checking on condition of fast g. turbulnece
        real trbfast(plond,plev)





c Local Variables
c ---------------

c------------------------CMT-------------------------------------
c
      real rztodt                      ! 1./ztodt
      real tmp
      real dqv(plond,plev,pcnst)       ! constituent diffusion tendency
      real denom

      real edle,rmu           !background heat/momentum diffusion
      real ed,pref,tref,fmu   ! parameters for computing edle and rmu
      real dttrb              ! internal time step for trbflx 
      real const
      integer itrtrb          ! number of internal time steps
      integer i,j,L,nn        ! grid index

c-----parameters for defining the background coefficients
      parameter ( fmu = 0.00067)
      PARAMETER ( TREF   =   258.     )  
      PARAMETER ( PREF   =   500.     )
      PARAMETER ( ED     =   0.5      )

 
      real  U(plond,plev)  ! u wind
      real  V(plond,plev)  ! v wind  
      real  TH(plond,plev)  !  model theta (theta/p0**k=tm1/p**k)
      real  SH(plond,plev)  !  specific humidity (KG/KG)
      real  THV(plond,plev) ! virtual potential temperature 


C Diagnostics Variables(controlla dopo quando si controlla trbflx)
C------------------------
      real kvh(plond,plev)   ! Final diffusion coefficient for heat
      real kvm(plond,plev)   ! Final diffusion coefficient for momentum
      real wt(plond,plev-1)  !heat flux
      real wsh(plond,plev-1) !flux for sp humidity
      real wu(plond,plev-1)  !flux for u-wind
      real wv(plond,plev-1)  !flux for v-wind
C----------------------------------------------------------------
c-----if calculation of diagnostics are needed
c      tprof=.true.
        tprof=.false.
C----------------------------------------------------------------
c      ITRTRB  is the number of time the time step in going to be split
c      If no splitting of dt is needed, set itrtrb=1     

       ITRTRB   =   6
       DTTRB    =   ZTODT / FLOAT(ITRTRB)
C
C     BACKGROUND HEAT/ MOMENTUM DIFFUSION
C
c************************************************************************
c      BACKGROUND COEFF ARE SET TO ZERO
c      EDLE    = ED * 0.2
c      RMU     = fmu * TREF * RAIR / PREF
       EDLE    = 0.0
       RMU     = 0.0
c*************************************************************************
c     Temperatures are transformed in model temperature,
c     theta _model=theta/p0**k=T/P**k
      
      do l=1,plev
      do i=1,plond
         th(i,l) = tm1(i,l)/prsk(i,l)
      enddo
      enddo   

C     We initialize a new wind and humidity field, so 
C     that all calculations are performed on the new
C     variables. We need to do that because at the end
C     we calculate tendencies and we need the 'old'
C     values of variables.
C
C     sh is the specific humidity, and it corresponds
C     to qm1(i,l,1). In this version we diffuse sh, in the
C     future we can diffuse different constituents.

      DO L=1,plev
      DO I =1,plond
       U(I,L)  = um1(i,l)
       V(I,L)  = vm1(i,l)
       SH(I,L) = qm1(i,l,1)
      ENDDO   
      ENDDO   



C COMPUTE VIRTUAL POTENTIAL TEMPERATURES
C --------------------------------------
c******************************************************************************
      do l=1,plev
      DO I =1,plond
      THV(I,L) = TH(I,L)*(1. + 0.61 * SH(I,L))
      enddo
      enddo
c

C ZERO OUT QQ ARRAY FOR BEGIN CASE
C --------------------------------
      IF(qbeg) then
        do l=1,plev
           do i=1,plond
              qq(i,l)=0.
            enddo
        enddo
      endif

C *********************************************************************
!-----Vertical diffusion with the 2.5 level DRY scheme
      CALL trbflx(tm1(1,plev),TH,THV,SH,U,V,QQ,P,PE,prsk,prske,PDELM1,
     & taux,tauy,shflx,cflx,DTTRB,ITRTRB,RMU,EDLE,qbeg,tprof,
     & lat,plond,plev,pcnst,gravit,cpair,rair,cappa,
     & wu,wv,wt,wsh,kvh,kvm,pblh,trbfast)
c
c****************************************************************************
C     Set the output winds to up1 and vp1

      do l=1,plev
      do i=1,plond
      up1(i,l)=u(i,l)
      vp1(i,l)=v(i,l)
      enddo
      enddo

C TENDENCY UPDATES
C ----------------
C
      rztodt=1./ztodt
      do  L=1,plev
        do I =1,plond
           duv(i,l)=(UP1(I,L)-UM1(I,L))*rztodt
           dvv(i,l)=(VP1(I,L)-VM1(I,L))*rztodt

           fu(i,l)=fu(i,l)+duv(i,l)
           fv(i,l)=fv(i,l)+dvv(i,l)

            denom=cpair*(1.+cpvir*qm1(i,l,1))

            tmp=tm1(i,l)    ! save original value

c     Calculates the form of temperature increment due to frictional heating
c     tm1 will be the final adjusted temperature     

      tm1(i,l)=th(i,l)*prsk(i,l)+
     &          ( (up1(i,l)+um1(i,l))*(um1(i,l)-up1(i,l)) +
     &                (vp1(i,l)+vm1(i,l))*(vm1(i,l)-vp1(i,l))
     &                 ) / (2.*denom)

c     Total heating rate due to pbl diffusion

      dtv(I,L)=(tm1(I,L)-tmp)*rztodt
c     Total (cumulative) heating rate :
      t2(i,l)=t2(i,l)+dtv(i,l)
c     
C     Calculates the humidity tendency dqv(i,l,1). In the future
C     this will just be a part of dqv(i,l,pcnst), for pcnst=1

      dqv(i,l,1) = (sh(i,l)-qm1(i,l,1))*rztodt
c**********************************momentaneo fino a quando non si implementa
c                                  il calcolo per i tracers sotto*****************
c------Sets value for new specific humidity
       qm1(i,l,1)=sh(i,l)
      enddo
      enddo
c********************************************************************************

C     Still to be implemented (calculation of tendency for constituents)
c      do m=1,pcnst
c        do l=1,plev
c          do i=1,plon
c              dqv(I,L,m)=(QP1(I,L,m)-QMM1(I,L,m))*rztodt
c            enddo
c           enddo
c          enddo



CSave the pbl output variables on the history file. The other
C vertical diffusion variables are saved later, after Rayleigh
C friction tendencies have been added.
C
c      call outfld(hattr(1,iPBLH), pblh ,plond,lat,hbuf)
c      call outfld(hattr(1,iUSTAR),ustar,plond,lat,hbuf)
c
C If you want to uncomment the following outfld calls you must modify
C the argument list to vdiff
C
      call outfld(hattr(1,iKVH),kvh,plond,lat,hbuf)
      call outfld(hattr(1,iKVM),kvm,plond,lat,hbuf)
c      call outfld(hattr(1,iTRBFAST),trbfast,plond,lat,hbuf)
c     call outfld(vdiffnam(1),dqv(1,1,1),plond,lat,hbuf)
c      call outfld(hattr(1,iVD01),dqv(1,1,1),plond,lat,hbuf)
c      do m=2,pcnst
c        call outfld(vdiffnam(m),dqv(1,1,m),plond,lat,hbuf)
c        call outfld(hattr(1,iVD01+m-1),dqv(1,1,m),plond,lat,hbuf)
c      end do
      return
      end
