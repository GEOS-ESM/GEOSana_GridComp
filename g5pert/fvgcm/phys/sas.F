      subroutine sas(imx2,imx22,kmx,delx,delt,
     &               q1,t1,psfc,pres, spd,
     &               rdp, pe, pln, slimsk,dot,
     &               g,rd,rv,cpa,hvap,t0c,
     &               cldwrk,rain,kbot,ktop,kuo,qci,qrs)
 
c subprogram documentation block
c                .      .    .                                       .
c subprogram:    sas      computes convective heating and moisng
c   prgmmr: hua-lu pan       org: w/nmc23    date: 92-03-01
!   Modified: SJ Lin, June 10, 2001
c
c abstract: computes convective heating and moistening using a one
c   cloud type arakawa-schubert convection scheme originally developed
c   by georg grell. the scheme includes updraft and downdraft effects.
c   the closure is the cloud work function. both updraft and downdraft
c   are assumed to be saturated and the heating and moistening are
c   accomplished by the compensating environment. the name comes from
c   "simplified arakawa-schubert convection parameterization".
c
! program history log:
!   92-03-01  hua-lu pan
!   99-02-19  kanamitsu & hong          plug-compatible version
c
c usage:    call sas(imx2,imx22,kmx,delx,delt,
c    &               q1,t1,psfc,pres,spd,slimsk,dot,
c    &               g,rd,rv,cpa,hvap,
c    &               cldwrk,rain,kbot,ktop,kuo,qci,qrs)
c
c   input argument list:
c     imx2     - integer number of points
c     imx22    - integer leading dimension of q1,t1, normally =imx2
c     kmx      - integer number of levels (index increases upward)
c     delt     - real time step in seconds
c     q1       - real (ix,km) current specific humidity in kg/kg
c     t1       - real (ix,km) current temperature in kelvin
c     psfc     - real (ix) surface pressure in pascal
c     pres     - real (ix,km) pressure where T,q are defined (in pascal).
c     spd      - real (ix,km) wind speed in m/sec.
c     slimsk   - real (im) land(1),sea(0),ice(2) flag
c     dot      - real (ix,km) omega in Pascal/sec
c     delx     - real model grid distance in meter
c     g        - real gravitataional acceleration in meter/sec**2
c     rd       - real gas constant of dry air (J/kg/K)
c     rv       - real gas constant of water vapor (J/kg/K)
c     cpa      - real specific heat of dry air at constant pressure (J/kg/K)
c     hvap     - real latent heat of condensation (J/kg)
c
c   output argument list:
c     q1       - real (ix,km) adjusted specific humidity in kg/kg
c     t1       - real (ix,km) adjusted temperature in kelvin
c     cldwrk   - real (imx2) cloud work function (J/kg)
c     rain     - real (imx2) convective rain in meters
c     kbot     - integer (imx2) cloud bottom level (model level index)
c     ktop     - integer (imx2) cloud top level (model level index)
c     kuo      - integer (imx2) bit flag indicating deep convection
c     qci      - real (imx2,kmx) cloud ice mixing ratio (gm/kg)
c     qrs      - real (imx2,kmx) cloud liquid water mixing ratio (gm/kg)
 
      implicit none

      integer, intent(in):: imx2,imx22,kmx

      real, intent(in):: delx          ! Model horizontal resoution
      real, intent(in):: delt          ! Model time step: leap-frog
      real, intent(in):: g,rd,rv,cpa,hvap,t0c
      real, intent(in):: psfc(imx2)
      real, intent(in):: pres(imx2,kmx)
      real, intent(in):: spd(imx2,kmx),slimsk(imx2),dot(imx2,kmx)
      real, intent(in):: pe(imx2,kmx+1)       ! pressure at layer edges
      real, intent(in):: pln(imx2,kmx+1)      ! log-pressure at layer edges
      real, intent(in):: rdp(imx2,kmx)        ! inverse of pressure thickness

      real, intent(inout):: t1(imx22,kmx)     ! temperature (K)
      real, intent(inout):: q1(imx22,kmx)     ! specific humidity

      integer, intent(out):: kbot(imx2),ktop(imx2),kuo(imx2)
      real, intent(out):: rain(imx2)
      real, intent(out):: cldwrk(imx2)
      real, intent(out):: qci(imx2,kmx)
      real, intent(out):: qrs(imx2,kmx)
 
!  local variables and arrays

!  physical parameters
 
      real qmin
      parameter(qmin=1.e-30)
      real cpoel,elocp,el2orc,eps,epsm1,fv
      real terr,c0
      parameter(terr=0.,c0=.002)
      real fact1,fact2
      real edtk1
      parameter(edtk1=3.e4)               ! 30 km
 
      real pdot(imx2),acrtfct(imx2)
      real to(imx2,kmx),qo(imx2,kmx)
      real qeso(imx2,kmx)
      real zo(imx2,kmx)            ! mid-layer height
      real tvo(imx2,kmx) 
      real dbyo(imx2,kmx)
      real heo(imx2,kmx),heso(imx2,kmx)
      real qrcd(imx2,kmx),dellah(imx2,kmx),dellaq(imx2,kmx)
      real hcko(imx2,kmx),qcko(imx2,kmx),eta(imx2,kmx)
      real etad(imx2,kmx),qrcdo(imx2,kmx),pwo(imx2,kmx)
      real pwdo(imx2,kmx),deltv(imx2),acrt(imx2)
      real dtconv(imx2)
 
      integer kb(imx2)
      integer kbcon(imx2)
      integer lmin(imx2), jmin(imx2)
      integer kds(imx2)
      integer ktcon(imx2)
      integer kbdtr(imx2)

      real hmax(imx2), delq(imx2)
      real hkbo(imx2),qkbo(imx2),pbcdif(imx2)
      real vmax(imx2),hmin(imx2)
      real pwavo(imx2),aa1(imx2),vshear(imx2)
      real edt(imx2),edto(imx2),pwevo(imx2),qcond(imx2)
      real hcdo(imx2),qcdo(imx2),ddp(imx2),pp2(imx2)
      real adet(imx2),aatmp(imx2),xhkb(imx2),xqkb(imx2)
      real xpwav(imx2),xpwev(imx2),xhcd(imx2),xaa0(imx2)
      real f(imx2),xk(imx2),xmb(imx2),edtx(imx2)
      real xqcd(imx2),hsbar(imx2),xmbmax(imx2),xlamb(imx2)
      real xlamd(imx2)
 
      real delhbar(imx2),delqbar(imx2),deltbar(imx2)
      real dpres(imx2)
 
      logical totflg
      logical cnvflg(imx2),dwnflg(imx2)
      logical dwnflg2(imx2),flg(imx2) 
 
      real pcrit(15), acritt(15), acrit(15)
      data pcrit/85000.,80000.,75000.,70000.,65000.,60000.,
     &           55000.,50000.,45000.,40000.,
     &           35000.,30000.,25000.,20000.,15000./
      data acritt/.0633,.0445,.0553,.0664,.075 , .1082, .1521,.2216,
     &            .3151,.3677,.41  ,.5255,.7663,1.1686,1.6851/
 
      integer im,km,i,k
      integer kbmax,kbm,kmax
      integer indx
      integer knumb 
      integer jmn

      real mbdt
      real dt2,dtmin,dtmax
      real edtmaxl,edtmaxs,alphal,alphas,betal,betas
      real evfact,pdpdwn,pdetrn
      real evef       ! evaporation efficiency (SJL)
      real t62delx,wratio,w1l,w2l,w3l,w4l,w1s,w2s,w3s,w4s
      real sl           ! psudo-sigma level
      real dlnsig
      real dz,dp,es,pprime,qs,dqsdp,desdt,dqsdt,gamma,dt,dq,po
      real alpha
      real factor,onemf
      real dz1,qrch,etah,qlk,qc
      real rfact
      real shear,e1
      real beta
      real dh,detad
      real edtmax
      real dhh,dg
      real aup,adw,dv1,dv2,dv3,dv1q,dv2q,dv3q,deta,dvq1,xqrch
      real xdby,dqrch,xqc,xpw,xpwd
      real w1,w2,w3,w4
      real dellat,qevap

#include <eslookup.h>

      real rdg
      real hkl
      real zi(imx2,kmx+1)
 
      im = imx2
      km = kmx

!  initialize arrays

      cpoel=cpa/hvap
      elocp=hvap/cpa
      el2orc=hvap*hvap/(rv*cpa)
      eps = rd/rv
      epsm1=rd/rv-1.
      fv=1./eps-1.
      fact1=(1.8460e+3-4.1855e+3)/rv
      fact2=hvap/rv-fact1*t0c
 
      rdg = rd/g

      do i=1,im
         zi(i,1) = 0.
      end do

      do k=1,kmx
        do i=1,im
           hkl = pln(i,k) - pln(i,k+1)
           tvo(i,k)  = t1(i,k) * (1.+fv*q1(i,k))
           zo(i,k)   = zi(i,k) + rdg * tvo(i,k) * (1.-hkl*pe(i,k+1)*rdp(i,k))
           zi(i,k+1) = zi(i,k) + rdg * tvo(i,k) * hkl
         end do
      end do

      do i=1,im
        rain(i) = 0.
        kbot(i) = km+1
        ktop(i) = 0
         kuo(i) = 0
        cnvflg(i) = .true.
        dtconv(i) = 3600.
        cldwrk(i) = 0.
        pdot(i) = 0.
        xmbmax(i) = .3
      enddo

      do k = 1, 15
        acrit(k) = acritt(k) * (97500. - pcrit(k))
      enddo

      dt2 = 2. * delt
      dtmin = max(dt2,1200.)
      dtmax = max(dt2,3600.)
!
!  model tunable parameters are all here
!
      mbdt = 10.
      edtmaxl = min(edtk1/delx,.9)
      edtmaxs = .3         ! original
!     edtmaxs = .2         ! over sea         
      alphal = .5
      alphas = .5
      betal = .05
      betas = .05
      evef = 0.07          ! original value (over land and ice)
!     evef = 0.08
      evfact = 1. - edtmaxs  ! evap factor over sea
      pdpdwn = 0.
      pdetrn = 20000.
 
      t62delx=40000000./192.
      wratio=t62delx/delx
      w1l = -2.e-1 * wratio
      w2l = -1.e-0 * wratio
      w3l = -2.e-1 * wratio
      w4l = -2.e-2 * wratio
      w1s = -2.e-1 * wratio
      w2s = -1.e-0 * wratio
      w3s = -2.e-1 * wratio
      w4s = -2.e-2 * wratio
 
!  define top layer for search of the downdraft originating layer
!  and the maximum thetae for updraft
 
      kbmax=1
      kbm=1
      kmax=1
      do k = 1, km
        do i = 1, im
         sl = pres(i,k)/psfc(i)
         if(sl > .45) kbmax = max(kbmax,k)
         if(sl > .70) kbm = max(kbm,k)
         if(sl > .05) kmax = max(kmax,k)
        enddo
      enddo
 
      do k = 1, km
        do i = 1, im
          qci(i,k) = 0.0
          qrs(i,k) = 0.0
          pwo(i,k) = 0.
          pwdo(i,k) = 0.
        enddo
      enddo
      do k = 1, kmax
        do i = 1, im
          pwo(i,k) = 0.
          pwdo(i,k) = 0.
          to(i,k) = t1(i,k)
          qo(i,k) = q1(i,k)
          dbyo(i,k) = 0.
        enddo
      enddo
 
c  column variables
c  pres is pressure of the layer (mb) changed to Pascal
c  to is temperature at t+dt (k)... this is after advection and turbulan
c  qo is mixing ratio at t+dt (kg/kg)..q1
 
      do k = 1, kmax
        do i = 1, im
          qeso(i,k) = estblf(t1(i,k))
          qeso(i,k) = eps*qeso(i,k)/(pres(i,k)+epsm1*qeso(i,k))
          qeso(i,k) = max(qeso(i,k),qmin)
        enddo
      enddo
 
!  compute moist static energy
 
      do k = 1, kmax
        do i = 1, im
           heo(i,k) = g*zo(i,k) + cpa*to(i,k) + hvap*qo(i,k)
          heso(i,k) = g*zo(i,k) + cpa*to(i,k) + hvap*qeso(i,k)
        enddo
      enddo
 
c  determine level with largest moist static energy
c  this is the level where updraft starts
 
      do i = 1, im
        hmax(i) = heo(i,1)
        kb(i) = 1
      enddo
      do k = 2, kbm
        do i = 1, im
          if(heo(i,k).gt.hmax(i).and.cnvflg(i)) then
            kb(i) = k
            hmax(i) = heo(i,k)
          endif
        enddo
      enddo
      do i = 1, im
        if(qo(i,kb(i)).lt.qmin) cnvflg(i) = .false.
      enddo
c
c  search for downdraft originating level above theta-e minimum
c
      do i = 1, im
        hmin(i) = heso(i,1)
        lmin(i) = kbmax
        jmin(i) = kbmax
      enddo
      do k = 2, kbmax
        do i = 1, im
          if(heso(i,k).lt.hmin(i).and.cnvflg(i)) then
             lmin(i) = k + 1
             hmin(i) = heso(i,k)
          endif
        enddo
      enddo
 
      do k = 1, kmax - 1
        do i = 1, im
        if(cnvflg(i)) then
          dz = .5 * (zo(i,k+1) - zo(i,k))
          dp = .5 * (pres(i,k+1) - pres(i,k))
          es = estblf(to(i,k+1))
          pprime = pres(i,k+1) + epsm1 * es
          qs = eps * es / pprime
          dqsdp = - qs / pprime
          desdt = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
          dqsdt = qs * pres(i,k+1) * desdt / (es * pprime)
          gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
          dt = (g * dz + hvap * dqsdp * dp) / (cpa * (1. + gamma))
          dq = dqsdt * dt + dqsdp * dp
          to(i,k) = to(i,k+1) + dt
          qo(i,k) = qo(i,k+1) + dq
! SJL         po = .5 * (pres(i,k) + pres(i,k+1))
          po = pe(i,k+1)
          qeso(i,k) = estblf(to(i,k))
          qeso(i,k) = eps * qeso(i,k) / (po + epsm1 * qeso(i,k))
          qeso(i,k) = max(qeso(i,k),qmin)
          heo(i,k) = g * zi(i,k+1) +
     &                cpa * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zi(i,k+1) +
     &                cpa * to(i,k) + hvap * qeso(i,k)
        endif
        enddo
      enddo
c
c  look for convective cloud base as the level of free convection
c
      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          hkbo(i) = heo(i,indx)
          qkbo(i) = qo(i,indx)
        endif
      enddo
      do i = 1, im
         flg(i) = cnvflg(i)
         kbcon(i) = kmax
      enddo
c
      do k = 1, kbmax
        do i = 1, im
          if(flg(i).and.k.gt.kb(i)) then
            hsbar(i) = heso(i,k)
            if(hkbo(i).gt.hsbar(i)) then
              flg(i) = .false.
              kbcon(i) = k
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          pbcdif(i) = -pres(i,kbcon(i)) + pres(i,kb(i))
          pdot(i) = dot(i,kbcon(i))
          if(pbcdif(i).gt.15000.) cnvflg(i) = .false.
          if(kbcon(i).eq.kmax) cnvflg(i) = .false.
        endif
      enddo
 
      do i = 1,im
         if(cnvflg(i)) then
            if(kbcon(i).ge.kbmax) cnvflg(i) = .false.
            pdot(i) = dot(i,kbcon(i))
         endif
      enddo
 
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
c
c  found lfc, can define rest of variables
c
c  determine entrainment rate between kb and kbcon
c
      do i = 1, im
        alpha = alphas
        if(slimsk(i) == 1.) alpha = alphal                ! over land
        if(cnvflg(i)) then
          if(kb(i).eq.1) then
            dz = zi(i,kbcon(i)) - zo(i,1)
          else
            dz = zi(i,kbcon(i)) - zi(i,kb(i))
          endif
          if(kbcon(i).ne.kb(i)) then
            xlamb(i) = -log(alpha) / dz
          else
            xlamb(i) = 0.
          endif
        endif
      enddo
c
c  determine updraft mass flux
c
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i)) then
            eta(i,k) = 1.
          endif
        enddo
      enddo
      do k = kbmax, 2, -1
        do i = 1, im
          if(cnvflg(i).and.k.lt.kbcon(i).and.k.ge.kb(i)) then
!           dz = .5 * (zo(i,k+1) - zo(i,k-1))
            dz = zi(i,k+1) - zi(i,k)
            eta(i,k) = eta(i,k+1) * exp(-xlamb(i) * dz)
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.kb(i).eq.1.and.kbcon(i).gt.1) then
          dz = .5 * (zo(i,2) - zo(i,1))
          eta(i,1) = eta(i,2) * exp(-xlamb(i) * dz)
        endif
      enddo
c
c  work up updraft cloud properties
c
      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          hcko(i,indx) = hkbo(i)
          qcko(i,indx) = qkbo(i)
          pwavo(i) = 0.
        endif
      enddo
c
c  cloud property below cloud base is modified by the entrainment proces
c
      do k = 2, kmax - 1
        do i = 1, im
          if(cnvflg(i).and.k.gt.kb(i).and.k.le.kbcon(i)) then
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            hcko(i,k) = factor * hcko(i,k-1) + onemf *
     &                  .5 * (heo(i,k) + heo(i,k+1))
            dbyo(i,k) = hcko(i,k) - heso(i,k)
          endif
          if(cnvflg(i).and.k.gt.kbcon(i)) then
            hcko(i,k) = hcko(i,k-1)
            dbyo(i,k) = hcko(i,k) - heso(i,k)
          endif
        enddo
      enddo
c
c  determine cloud top
c
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon(i) = 1
      enddo
c
c
c
      do k = 2, kmax
        do i = 1, im
          if(dbyo(i,k).lt.0..and.flg(i).and.k.gt.kbcon(i)) then
            ktcon(i) = min(kmax-1,k)
            flg(i) = .false.
          endif
        enddo
      enddo
c
      do i = 1, im
        dpres(i)=pres(i,kbcon(i))-pres(i,ktcon(i))
      enddo
c
c    check cloud depth
c
      do i = 1, im
        if(cnvflg(i).and.dpres(i).lt.15000.) cnvflg(i) = .false.
      enddo
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
 
!  detraining cloud
 
      do i = 1, im
        if(cnvflg(i)) then
          dz = zi(i,ktcon(i)) - zi(i,kbcon(i))
          xlamb(i) = -1. / dz
        endif
        dwnflg(i) = cnvflg(i)
        if(cnvflg(i).and.dpres(i) > pdetrn) dwnflg(i)=.false.
        dwnflg2(i) = cnvflg(i)
        if(cnvflg(i).and.dpres(i) < pdpdwn) dwnflg2(i)=.false.
      enddo

      do k = 2, kmax - 1
        do i = 1, im
          if(dwnflg(i).and.k.gt.kbcon(i).and.k.le.ktcon(i)) then
!           dz = .5 * (zo(i,k+1) - zo(i,k-1))
            dz = zi(i,k+1) - zi(i,k)
            eta(i,k) = eta(i,k-1) * exp( xlamb(i) * dz)
          endif
        enddo
      enddo
c
c  cloud property above cloud top is modified by the detrainment process
c
      do k = 2, kmax - 1
        do i = 1, im
          if(dwnflg(i).and.k.gt.kbcon(i).and.k.le.ktcon(i)) then
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            hcko(i,k) = factor * hcko(i,k-1) + onemf *
     &                  .5 * (heo(i,k) + heo(i,k+1))
            dbyo(i,k) = hcko(i,k) - heso(i,k)
          endif
        enddo
      enddo
c
c  make sure that jmin is within the cloud
c
      do i = 1, im
        if(cnvflg(i)) then
          jmin(i) = min(lmin(i),ktcon(i)-1)
          jmin(i) = max(jmin(i),kbcon(i)+1)
        endif
      enddo
      do i = 1, im
        if(cnvflg(i).and.dwnflg2(i).and.jmin(i).le.kbcon(i))
     &     then
          cnvflg(i) = .false.
          dwnflg(i) = .false.
          dwnflg2(i) = .false.
        endif
      enddo
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
c
c  compute cloud moisture property and precipitation
c
      do i = 1, im
          aa1(i) = 0.
      enddo
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.k.gt.kb(i).and.k.lt.ktcon(i)) then
!           dz = .5 * (zo(i,k+1) - zo(i,k-1))
            dz = zi(i,k+1) - zi(i,k)
            dz1 = zo(i,k) - zo(i,k-1)
            gamma = el2orc * qeso(i,k) / (to(i,k)**2)
            qrch = qeso(i,k)
     &           + gamma * dbyo(i,k) / (hvap * (1. + gamma))
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            qcko(i,k) = factor * qcko(i,k-1) + onemf *
     &                  .5 * (qo(i,k) + qo(i,k+1))
            dq = eta(i,k) * qcko(i,k) - eta(i,k) * qrch
c
c  below lfc check if there is excess moisture to release latent heat
c
            if(dq.gt.0.) then
              etah = .5 * (eta(i,k) + eta(i,k-1))
              qlk = dq / (eta(i,k) + etah * c0 * dz)
              aa1(i) = aa1(i) - dz1 * g * qlk
              qc = qlk + qrch
              pwo(i,k) = etah * c0 * dz * qlk
              qcko(i,k) = qc
              pwavo(i) = pwavo(i) + pwo(i,k)
              qci(i,k) = qlk
            endif
          endif
        enddo
      enddo
c
c  calculate cloud work function at t+dt
c
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.k.gt.kbcon(i).and.k.lt.ktcon(i)) then
            dz1 = zo(i,k) - zo(i,k-1)
            gamma = el2orc * qeso(i,k-1) / (to(i,k-1)**2)
            rfact =  1. + fv * cpa * gamma
     &               * to(i,k-1) / hvap
            aa1(i) = aa1(i) +
     &               dz1 * (g / (cpa * to(i,k-1)))
     &               * dbyo(i,k-1) / (1. + gamma)
     &               * rfact
            aa1(i)=aa1(i)+
     &               dz1 * g * fv *
     &               max(0.,(qeso(i,k-1) - qo(i,k-1)))
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) dwnflg(i) = .false.
        if(cnvflg(i).and.aa1(i).le.0.) dwnflg2(i) = .false.
        if(cnvflg(i).and.aa1(i).le.0.) cnvflg(i) = .false.
      enddo
c
c------- downdraft calculations
c
c  determine level with largest wind speed
      do i = 1, im
          vmax(i) = spd(i,kb(i))
          kds(i) = kb(i)
      enddo
      do k = 2, kmax - 1
        do i = 1, im
          if(k.ge.kb(i).and.spd(i,k).gt.vmax(i).and.cnvflg(i)
     &       .and.k.le.ktcon(i)) then
            vmax(i) = spd(i,k)
            kds(i) = k
          endif
        enddo
      enddo
c
c--- determine downdraft strength in terms of windshear
c
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 1, kmax
        do i = 1, im
          if(k.ge.kb(i).and.k.le.ktcon(i).and.cnvflg(i)) then
            shear = abs((spd(i,k+1)-spd(i,k))/(zo(i,k+1)-zo(i,k)))
            vshear(i) = vshear(i) + shear
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          knumb = ktcon(i) - kb(i) + 1
          knumb = max(knumb,1)
          vshear(i) = 1.e3 * vshear(i) / float(knumb)
          if(vshear(i).ge.1.35) then
          e1 = 1.591-.639*vshear(i)
     &       +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          else
          e1 = 0.9
          endif
          edt(i)=1.-e1
          edt(i) = min(edt(i),1.0)
          edt(i) = max(edt(i),.1)
          edto(i)=edt(i)
          edtx(i)=edt(i)
        endif
      enddo
c
c  determine detrainment rate between 1 and kbdtr
c
      do i = 1, im
        kbdtr(i) = kbcon(i)
        beta = betas
        if(slimsk(i) == 1.) beta = betal
        if(cnvflg(i)) then
          kbdtr(i) = kbcon(i)
          kbdtr(i) = max(kbdtr(i),1)
          xlamd(i) = 0.
          if(kbdtr(i).gt.1) then
            dz = zi(i,kbdtr(i)) - zo(i,1)
            xlamd(i) = log(beta) / dz
          endif
        endif
      enddo
c
c  determine downdraft mass flux
c
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i)) then
            etad(i,k) = 1.
          endif
          qrcdo(i,k) = 0.
        enddo
      enddo

      do k = kbmax, 2, -1
        do i = 1, im
          if(cnvflg(i).and.k.lt.kbdtr(i)) then
!           dz = .5 * (zo(i,k+1) - zo(i,k-1))
            dz = zi(i,k+1) - zi(i,k)
            etad(i,k) = etad(i,k+1) * exp(xlamd(i) * dz)
          endif
        enddo
      enddo

      k = 1
      do i = 1, im
        if(cnvflg(i).and.kbdtr(i).gt.1) then
          dz = .5 * (zo(i,2) - zo(i,1))
          etad(i,k) = etad(i,k+1) * exp(xlamd(i) * dz)
        endif
      enddo
 
!--- downdraft moisture properties
 
      do i = 1, im
        pwevo(i) = 0.
        flg(i) = cnvflg(i)
      enddo

      do i = 1, im
        if(cnvflg(i)) then
          jmn = jmin(i)
          hcdo(i) = heo(i,jmn)
          qcdo(i) = qo(i,jmn)
          qrcdo(i,jmn) = qeso(i,jmn)
        endif
      enddo

      do k = kmax-1, 1, -1
        do i = 1, im
          if(cnvflg(i).and.k.lt.jmin(i)) then
            dq = qeso(i,k)
            dt = to(i,k)
            gamma = el2orc * dq / dt**2
            dh = hcdo(i) - heso(i,k)
            qrcdo(i,k)=dq+(1./hvap)*(gamma/(1.+gamma))*dh
            detad = etad(i,k+1) - etad(i,k)
            pwdo(i,k) = etad(i,k+1) * qcdo(i) -
     &                 etad(i,k) * qrcdo(i,k)
              pwdo(i,k) = pwdo(i,k) - detad *
     &                   .5 * (qrcdo(i,k) + qrcdo(i,k+1))
            qcdo(i) = qrcdo(i,k)
            pwevo(i) = pwevo(i) + pwdo(i,k)
          endif
        enddo
      enddo
 
!--- final downdraft strength dependent on precip
!--- efficiency (edt), normalized condensate (pwav), and
!--- evaporate (pwev)
 
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i) == 0.) edtmax = edtmaxs     ! over sea
        if(dwnflg2(i)) then
          if(pwevo(i).lt.0.) then
            edto(i) = -edto(i) * pwavo(i) / pwevo(i)
            if(edto(i).ge..90) then
               cnvflg(i) = .false.
               dwnflg2(i) = .false.
            endif
            edto(i) = min(edto(i),edtmax)
          else
            edto(i) = 0.
          endif
        else
          edto(i) = 0.
        endif
      enddo
 
!--- downdraft cloudwork functions
 
      do k = kmax-1, 1, -1
        do i = 1, im
          if(dwnflg2(i).and.k.lt.jmin(i)) then
            gamma = el2orc * qeso(i,k) / to(i,k)**2
            dhh=hcdo(i)
            dt=to(i,k)
            dg=gamma
            dh=heso(i,k)
            dz=-1.*(zo(i,k+1)-zo(i,k))
            aa1(i)=aa1(i)+edto(i)*dz*(g/(cpa*dt))*((dhh-dh)/(1.+dg))
     &             *(1.+fv*cpa*dg*dt/hvap)
            aa1(i)=aa1(i)+edto(i)*
     &      dz*g*fv*max(0.,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
      do i = 1, im
        if(aa1(i).le.0.) cnvflg(i) = .false.
        if(aa1(i).le.0.) dwnflg(i) = .false.
        if(aa1(i).le.0.) dwnflg2(i) = .false.
      enddo
 
c--- what would the change be, that a cloud with unit mass
c--- will do to the environment?
 
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
          endif
        enddo
      enddo

      do i = 1, im
        if(cnvflg(i)) then
! SJL          dp =psfc(i)-(pres(i,1)+pres(i,2))*0.5
          dp = psfc(i) - pe(i,2)
          dellah(i,1) = edto(i) * etad(i,1) * (hcdo(i)
     &                - heo(i,1)) * g / dp
          dellaq(i,1) = edto(i) * etad(i,1) * (qcdo(i)
     &                - qo(i,1)) * g / dp
        endif
      enddo
c
c--- changed due to subsidence and entrainment
c
      do k = 2, kmax-1
        do i = 1, im
          if(cnvflg(i).and.k.lt.ktcon(i)) then
            aup = 1.
            if(k.le.kb(i)) aup = 0.
            adw = 1.
            if(k.gt.jmin(i)) adw = 0.
            dv1= heo(i,k)
            dv2 = .5 * (heo(i,k) + heo(i,k+1))
            dv3= heo(i,k-1)
            dv1q= qo(i,k)
            dv2q = .5 * (qo(i,k) + qo(i,k+1))
            dv3q= qo(i,k-1)
! SJL           dp = (pres(i,k-1)-pres(i,k+1))*0.5
            dp = pe(i,k) - pe(i,k+1)
!           dz = .5 * (zo(i,k+1) - zo(i,k-1))
            dz = zi(i,k+1) - zi(i,k)
            deta = eta(i,k) - eta(i,k-1)
            detad = etad(i,k) - etad(i,k-1)
            dellah(i,k) = dellah(i,k) +
     &          ((aup * eta(i,k) - adw * edto(i) * etad(i,k)) * dv1
     &      - (aup * eta(i,k-1) - adw * edto(i) * etad(i,k-1))* dv3
     &                   - aup * deta * dv2
     &                  + adw * edto(i) * detad * hcdo(i)) * g / dp
            dellaq(i,k) = dellaq(i,k) +
     &          ((aup * eta(i,k) - adw * edto(i) * etad(i,k)) * dv1q
     &      - (aup * eta(i,k-1) - adw * edto(i) * etad(i,k-1))* dv3q
     &                  - aup * deta * dv2q
     &     +adw*edto(i)*detad*.5*(qrcdo(i,k)+qrcdo(i,k-1))) * g / dp
          endif
        enddo
      enddo
 
c------- cloud top
 
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
!         dp = (pres(i,indx-1)-pres(i,indx+1))*0.5
          dp = pe(i,indx) - pe(i,indx+1)
          dv1 = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *
     &                     (hcko(i,indx-1) - dv1) * g / dp
          dvq1 = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *
     &                     (qcko(i,indx-1) - dvq1) * g / dp
        endif
      enddo
 
c------- final changed variable per unit mass flux
 
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.k.gt.ktcon(i)) then
            qo(i,k) = q1(i,k)
            to(i,k) = t1(i,k)
          endif
          if(cnvflg(i).and.k.le.ktcon(i)) then
            qo(i,k) = dellaq(i,k) * mbdt + q1(i,k)
            dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cpa
            to(i,k) = dellat * mbdt + t1(i,k)
            qo(i,k) = max(qo(i,k),qmin)
          endif
        enddo
      enddo
 
c--- the above changed environment is now used to calulate the
c--- effect the arbitrary cloud (with unit mass flux)
c--- would have on the stability,
c--- which then is used to calculate the real mass flux,
c--- necessary to keep this change in balance with the large-scale
c--- destabilization.
c--- environmental conditions again, first heights
 
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i)) then
            qeso(i,k) = estblf(to(i,k))
            qeso(i,k) = eps*qeso(i,k)/(pres(i,k)+epsm1*qeso(i,k))
            qeso(i,k) = max(qeso(i,k),qmin)
            tvo(i,k) = to(i,k) + fv * to(i,k) * max(qo(i,k),qmin)
          endif
        enddo
      enddo

      do i = 1, im
        if(cnvflg(i)) then
          xaa0(i) = 0.
          xpwav(i) = 0.
        endif
      enddo
 
!  hydrostatic height assume zero terr

      do k=1,kmx
        do i=1,im
        if(cnvflg(i)) then
           hkl = pln(i,k) - pln(i,k+1)
           zo(i,k)   = zi(i,k) + rdg * tvo(i,k) * (1.-hkl*pe(i,k+1)*rdp(i,k))
           zi(i,k+1) = zi(i,k) + rdg * tvo(i,k) * hkl
        endif
        end do
      end do
 
!--- moist static energy
 
      do k = 1, kmax - 1
        do i = 1, im
          if(cnvflg(i)) then
            dz = .5 * (zo(i,k+1) - zo(i,k))
            dp = .5 * (pres(i,k+1) - pres(i,k))
            es = estblf(to(i,k+1))
            pprime = pres(i,k+1) + epsm1 * es
            qs = eps * es / pprime
            dqsdp = - qs / pprime
            desdt = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt = qs * pres(i,k+1) * desdt / (es * pprime)
            gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt = (g * dz + hvap * dqsdp * dp) / (cpa * (1. + gamma))
            dq = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
!           po = .5 * (pres(i,k) + pres(i,k+1))
            po = pe(i,k+1)
            qeso(i,k) = estblf(to(i,k))
            qeso(i,k) = eps * qeso(i,k) / (po + epsm1 * qeso(i,k))
            qeso(i,k) = max(qeso(i,k),qmin)
            heo(i,k) = g * zi(i,k+1) +
     &                  cpa * to(i,k) + hvap * qo(i,k)
            heso(i,k) = g * zi(i,k+1) +
     &                  cpa * to(i,k) + hvap * qeso(i,k)
          endif
        enddo
      enddo

      k = kmax
      do i = 1, im
        if(cnvflg(i)) then
           heo(i,k) = g * zo(i,k) + cpa * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zo(i,k) + cpa * to(i,k) + hvap * qeso(i,k)
        endif
      enddo

      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          xhkb(i) = heo(i,indx)
          xqkb(i) = qo(i,indx)
          hcko(i,indx) = xhkb(i)
          qcko(i,indx) = xqkb(i)
        endif
      enddo
 
!   static control
 
!   moisture and cloud work functions
 
      do k = 2, kmax - 1
        do i = 1, im
          if(cnvflg(i).and.k.gt.kb(i).and.k.le.ktcon(i)) then
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            hcko(i,k) = factor * hcko(i,k-1) + onemf *
     &                  .5 * (heo(i,k) + heo(i,k+1))
          endif
        enddo
      enddo

      do k = 2, kmax - 1
        do i = 1, im
          if(cnvflg(i).and.k.gt.kb(i).and.k.lt.ktcon(i)) then
!           dz = .5 * (zo(i,k+1) - zo(i,k-1))
            dz = zi(i,k+1) - zi(i,k)
            gamma = el2orc * qeso(i,k) / (to(i,k)**2)
            xdby = hcko(i,k) - heso(i,k)
            xdby = max(xdby,0.)
            xqrch = qeso(i,k)
     &           + gamma * xdby / (hvap * (1. + gamma))
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            qcko(i,k) = factor * qcko(i,k-1) + onemf *
     &                  .5 * (qo(i,k) + qo(i,k+1))
            dq = eta(i,k) * qcko(i,k) - eta(i,k) * xqrch
            if(dq.gt.0.) then
              etah = .5 * (eta(i,k) + eta(i,k-1))
              qlk = dq / (eta(i,k) + etah * c0 * dz)
              xaa0(i) = xaa0(i) - (zo(i,k) - zo(i,k-1)) * g * qlk
              xqc = qlk + xqrch
              xpw = etah * c0 * dz * qlk
              qcko(i,k) = xqc
              xpwav(i) = xpwav(i) + xpw
            endif
          endif
          if(cnvflg(i).and.k.gt.kbcon(i).and.k.lt.ktcon(i)) then
            dz1 = zo(i,k) - zo(i,k-1)
            gamma = el2orc * qeso(i,k-1) / (to(i,k-1)**2)
            rfact =  1. + fv * cpa * gamma
     &               * to(i,k-1) / hvap
            xdby = hcko(i,k-1) - heso(i,k-1)
            xaa0(i) = xaa0(i)
     &              + dz1 * (g / (cpa * to(i,k-1)))
     &              * xdby / (1. + gamma)
     &              * rfact
            xaa0(i)=xaa0(i)+
     &               dz1 * g * fv *
     &               max(0.,(qeso(i,k-1) - qo(i,k-1)))
          endif
        enddo
      enddo
c
c------- downdraft calculations
c
c--- downdraft moisture properties
c
      do i = 1, im
        xpwev(i) = 0.
      enddo

      do i = 1, im
        if(dwnflg2(i)) then
          jmn = jmin(i)
          xhcd(i) = heo(i,jmn)
          xqcd(i) = qo(i,jmn)
          qrcd(i,jmn) = qeso(i,jmn)
        endif
      enddo

      do k = kmax-1, 1, -1
        do i = 1, im
          if(dwnflg2(i).and.k.lt.jmin(i)) then
            dq = qeso(i,k)
            dt = to(i,k)
            gamma = el2orc * dq / dt**2
            dh = xhcd(i) - heso(i,k)
            qrcd(i,k)=dq+(1./hvap)*(gamma/(1.+gamma))*dh
            detad = etad(i,k+1) - etad(i,k)
            xpwd = etad(i,k+1) * qrcd(i,k+1) -
     &                 etad(i,k) * qrcd(i,k)
            xpwd = xpwd - detad *
     &             .5 * (qrcd(i,k) + qrcd(i,k+1))
            xpwev(i) = xpwev(i) + xpwd
          endif
        enddo
      enddo
c
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i) == 0.) edtmax = edtmaxs
        if(dwnflg2(i)) then
          if(xpwev(i).ge.0.) then
            edtx(i) = 0.
          else
            edtx(i) = -edtx(i) * xpwav(i) / xpwev(i)
            if(edtx(i).ge..90) then
               cnvflg(i) = .false.
               dwnflg2(i) = .false.
            endif
            edtx(i) = min(edtx(i),edtmax)
          endif
        else
          edtx(i) = 0.
        endif
      enddo
c
c--- downdraft cloudwork functions
c
      do k = kmax-1, 1, -1
        do i = 1, im
          if(dwnflg2(i).and.k.lt.jmin(i)) then
            gamma = el2orc * qeso(i,k+1) / to(i,k+1)**2
            dhh=xhcd(i)
            dt= to(i,k)
            dg= gamma
            dh= heso(i,k)
            dz=-1.*(zo(i,k+1)-zo(i,k))
            xaa0(i)=xaa0(i)+edtx(i)*dz*(g/(cpa*dt))*((dhh-dh)/(1.+dg))
     &              *(1.+fv*cpa*dg*dt/hvap)
            xaa0(i)=xaa0(i)+edtx(i)*
     &      dz*g*fv*max(0.,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
c
c  calculate critical cloud work function
c
      do i = 1, im
        acrt(i) = 0.
        if(cnvflg(i)) then
          if(pres(i,ktcon(i)) < pcrit(15))then
            acrt(i)=acrit(15)*(97500.-pres(i,ktcon(i)))
     &              /(97500.-pcrit(15))
          else if(pres(i,ktcon(i)).gt.pcrit(1))then
            acrt(i)=acrit(1)
          else
            k = int((85000. - pres(i,ktcon(i)))/5000.) + 2
            k = min(k,15)
            k = max(k,2)
            acrt(i)=acrit(k)+(acrit(k-1)-acrit(k))*
     *           (pres(i,ktcon(i))-pcrit(k))/(pcrit(k-1)-pcrit(k))
           endif
         endif
      enddo

      do i = 1, im
        acrtfct(i) = 1.
        w1 = w1s
        w2 = w2s
        w3 = w3s
        w4 = w4s
        if(slimsk(i) == 1.) then
          w1 = w1l
          w2 = w2l
          w3 = w3l
          w4 = w4l
        endif
        if(cnvflg(i)) then
          if(pdot(i).le.w4) then
            acrtfct(i) = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            acrtfct(i) = (pdot(i) + w4) / (w4 - w3)
          else
            acrtfct(i) = 0.
          endif
          acrtfct(i) = max(acrtfct(i),-1.)
          acrtfct(i) = min(acrtfct(i),1.)
          acrtfct(i) = 1. - acrtfct(i)
          dtconv(i) = 1800. * (pdot(i) - w2) / (w1 - w2)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = min(dtconv(i),dtmax)
        endif
      enddo
c
c--- large scale forcing
c
      do i= 1, im
        flg(i) = cnvflg(i)
        if(cnvflg(i)) then
          f(i) = (aa1(i) - acrt(i)/100. * acrtfct(i)) / dtconv(i)
          if(f(i).le.0.) flg(i) = .false.
        endif
        cnvflg(i) = flg(i)
        if(cnvflg(i)) then
          xk(i) = (xaa0(i) - aa1(i)) / mbdt
          if(xk(i).ge.0.) flg(i) = .false.
        endif
c
c--- kernel, cloud base mass flux
c
        cnvflg(i) = flg(i)
        if(cnvflg(i)) then
          xmb(i) = -f(i) / xk(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
c
c  restore t0 and qo to t1 and q1 in case convection stops
c
      do k = 1, kmax
        do i = 1, im
          to(i,k) = t1(i,k)
          qo(i,k) = q1(i,k)
          qeso(i,k) = estblf(t1(i,k))
          qeso(i,k) = eps * qeso(i,k) / (pres(i,k) + epsm1 * qeso(i,k))
          qeso(i,k) = max(qeso(i,k),qmin)
        enddo
      enddo
c
c--- feedback: simply the changes from the cloud with unit mass flux
c---           multiplied by  the mass flux necessary to keep the
c---           equilibrium with the larger-scale.
c
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
          qcond(i) = 0.
      enddo

      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.k.le.ktcon(i)) then
            aup = 1.
            if(k.le.kb(i)) aup = 0.
            adw = 1.
            if(k.gt.jmin(i)) adw = 0.
            dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cpa
            t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
            q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
!
!           if(k == 1) then
!             dp = psfc(i)-.5*(pres(i,1)+pres(i,2))
!           else
!             dp = (pres(i,k-1)-pres(i,k+1))*0.5
!           endif
!
              dp = pe(i,k) - pe(i,k+1)

            delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
            delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
            deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
          endif
        enddo
      enddo

      do i = 1, im
        delqbar(i) = 0.
        deltbar(i) = 0.
      enddo

      do k = kmax, 1, -1
        do i = 1, im
          dellat = 0.
          dellaq(i,k) = 0.
          if(cnvflg(i).and.k.le.ktcon(i)) then
            if(slimsk(i) == 0.) evef = evfact * edt(i)
            aup = 1.
            if(k.le.kb(i)) aup = 0.
            adw = 1.
            if(k.gt.jmin(i)) adw = 0.
            rain(i) = rain(i)
     &            + (aup * pwo(i,k) + adw * edto(i) * pwdo(i,k))
     &            * xmb(i) * .001 * dt2
            qrs(i,k) = aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
            qcond(i) = evef * (qo(i,k) - qeso(i,k)) / (1. + el2orc *
     &               qeso(i,k) / to(i,k)**2)

!           if(k.eq.1) then
!             dp = psfc(i)-.5*(pres(i,1)+pres(i,2))
!           else
!             dp = (pres(i,k-1)-pres(i,k+1))*0.5
!           endif
!
              dp = pe(i,k) - pe(i,k+1)
!
            if(rain(i).gt.0..and.qcond(i).le.0.) then
              qevap = -qcond(i) * (1. - exp(-.32 * sqrt(dt2 * rain(i))))
              qevap = min(qevap, rain(i)*1000.*g/dp)
              q1(i,k) = q1(i,k) + qevap
              t1(i,k) = t1(i,k) - elocp * qevap
              rain(i) = rain(i) - .001 * qevap * dp / g
              dellat = - elocp*qevap/xmb(i)/dt2
              dellaq(i,k) =  + qevap/xmb(i)/dt2
            endif
            delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
            deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
          endif
        enddo
      enddo
 
c  precipitation rate converted to actual precip
c  in unit of m instead of kg
 
      do i = 1, im
        if(cnvflg(i)) then
 
c  in the event of upper level rain evaporation and lower level downdraf
c  moistening, rain can become negative, in this case, we back out of the
c  heating and the moistening
 
          if(rain(i).le.0.) then
            rain(i) = 0.
          else
            ktop(i) = ktcon(i)
            kbot(i) = kbcon(i)
            kuo(i) = 1
            cldwrk(i) = aa1(i)
          endif
        endif
      enddo

      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.rain(i).le.0.) then
             t1(i,k) = to(i,k)
             q1(i,k) = qo(i,k)
          endif
        enddo
      enddo

      return
      end
