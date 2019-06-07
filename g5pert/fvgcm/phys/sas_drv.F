      subroutine sas_drv(im,    jm,  km,  k1,  dt,  coslat,
     &                   qstar, tstar,    ps,  pm,   wspd,
     &                   pe,    peln,    rpdel,
     &                   oro, omga, ae,  g,    rair,
     &                   rh2o, cpair, latvap, tmelt, 
     &                   precc, cnt, cnb)

! Driver for NCEP's SAS cumulus parameterization
! Programmer: S.-J. Lin, June 11, 2001

      implicit none
! Input
      integer im, jm, km, k1
      real qstar(im,km)
      real tstar(im,km)
      real  wspd(im,km)
      real  omga(im,km)
      real    pm(im,km)
      real rpdel(im,km)
      real    pe(im,km+1)
      real  peln(im,km+1)
      real    ps(im)
      real   oro(im)
      real ae               ! radius of earth
      real g
      real rair, rh2o, cpair, latvap, tmelt
      real dt
      real coslat

! Output
      real cnt(im)               ! Top level of convective activity
      real cnb(im)               ! Bottom level of convective activity
      real precc(im)             ! 

! Local
      real delh                  ! model's grid resolution in meters
      real pi
      real rdt
      real dts
      integer i, k, it

      real tn(im,k1)
      real qn(im,k1)
      real wn(im,k1)
      real sn(im,k1)
      real pn(im,k1)
      real qsat(im,k1)
      real est(im,k1)

      real pen(im,k1+1)
      real pln(im,k1+1)
      real rdp(im,k1)

      real cldi(im,k1)           ! cloud ice mixing ratio (gm/kg)
      real cldw(im,k1)           ! cloud liquid water mixing ratio (gm/kg)
      real cldwrk(im)            ! cloud work function
      real rain(im)
      integer kbot(im), kb(im)
      integer ktop(im), kt(im)
      integer kuo (im)
      integer k2
      integer itmax
      real evf, pflx, dq, dpog, rh
      real hice                   ! latent heat ice phase
      real rsum(im)
      real isum(im)

      evf  = 5.e-6
      hice = 3.34E5

! Determine subcycling counts:
      itmax = nint(dt/1800.)    ! 30-min base time step
      itmax = max(1, itmax)
      dts = 0.5*dt/float(itmax)    ! NCEP assumed leap-frog is THE scheme

      pi = 4.*atan(1.0)
! Equivalent resolution: delh = sqrt(dx*dy)
      delh = ae*pi*sqrt(2.*coslat/float((jm-1)*im))

! The NCEP model's k increases from bottom to top

      do k=1,k1
            k2 = km-k+1
         do i=1,im
            tn(i,k) = tstar(i,k2)
            qn(i,k) = qstar(i,k2)
            wn(i,k) =  omga(i,k2)
            sn(i,k) =  wspd(i,k2)
            pn(i,k) =    pm(i,k2)
           rdp(i,k) = rpdel(i,k2)
         enddo
      enddo

      do k=1,k1+1
            k2 = km-k+2
         do i=1,im
           pen(i,k) =   pe(i,k2)
           pln(i,k) = peln(i,k2)
         enddo
      enddo

      do i=1,im
         precc(i) = 0.
         ktop(i) = 1
         kbot(i) = k1+1
      enddo

      do it=1, itmax
!
! Note: NCEP folks think everyone should be using leap-frog scheme
! therefore the input dts is the actual time step * 0.5
!
      call sas(im, im, k1, delh, dts,
     &         qn, tn, ps, pn, sn,
     &         rdp, pen, pln, oro, wn,
     &         g, rair, rh2o, cpair, latvap, tmelt,
     &         cldwrk, rain, kb, kt, kuo, cldi, cldw)

! Due to the lack of a cloud water-ice scheme we re-evaporate cloud H2O
! using very simple Sundqvist-like scheme and cloud ice is melted if
! t > tmelt

      call aqsat (tn,  pn,  est,   qsat,   im,
     &            im,  k1,  k1,    k1   )

      do i=1,im
         rsum(i) = 0.
         isum(i) = 0.
      enddo

      do k=k1,1,-1
        do i=1,im
          dpog = (pen(i,k)-pen(i,k+1)) / g

! Convert to kg / kg
	  cldi(i,k) = 1.e-3*cldi(i,k)
          cldw(i,k) = 1.e-3*cldw(i,k)

          isum(i) = isum(i) + cldi(i,k)*dpog

          if ( tn(i,k) >  tmelt ) then
! Melting of cloud ice
             dq = min( isum(i)/dpog,
     &                (tn(i,k)-tmelt)*cpair/hice )
             cldw(i,k) = cldw(i,k) + dq
             tn(i,k) = tn(i,k) - dq*hice/cpair
             isum(i) = isum(i) - dq*dpog
          endif

! pflx unit: kg / (s*m^2)
          rsum(i) = rsum(i) + cldw(i,k)*dpog
          pflx = rsum(i) / (2.*dts)
          rh = qn(i,k) / qsat(i,k)
          dq = max(0., evf*(1.-rh)*sqrt(pflx)*2.*dts)
          dq = min(dq,rsum(i)/dpog,max(0.,0.8*(qsat(i,k)-qn(i,k))))
          qn(i,k) = qn(i,k) + dq
          tn(i,k) = tn(i,k) - dq * (latvap/cpair)
          rsum(i) = rsum(i) - dq*dpog
         enddo
      enddo

      do i=1,im
         precc(i) = precc(i) + rain(i) + 1.e-3*(rsum(i)+isum(i))
         ktop(i) = max(ktop(i), kt(i))
         kbot(i) = min(kbot(i), kb(i))
      enddo
 
      enddo  ! end iteraction

! Flip back to original
      do k=km-k1+1,km
            k2 = km-k+1
         do i=1,im
            tstar(i,k) = tn(i,k2) 
            qstar(i,k) = qn(i,k2) 
         enddo
      enddo

! Convert to m/sec
      rdt = 1./ dt

      do i=1,im
         precc(i) = precc(i) * rdt
         cnt(i) = km - (ktop(i) - 1)
         cnb(i) = km - (kbot(i) - 1)
      enddo

      return
      end
