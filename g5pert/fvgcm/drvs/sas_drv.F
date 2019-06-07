      subroutine sas_drv(im,    jm,  km,  k1,  dt,  coslat,
     &                   qstar, tstar,    ps,  pm,   wspd,
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
      real cldi(im,k1)           ! cloud ice mixing ratio
      real cldw(im,k1)           ! cloud liquid water mixing ratio
      real cldwrk(im)            ! cloud work function
      real rain(im)
      integer kbot(im), kb(im)
      integer ktop(im), kt(im)
      integer kuo (im)
      integer k2
      integer itmax

! Determine subcycling counts:
      itmax = nint ( dt / 1800. )    ! 30-min base time step
      itmax = max(1, itmax)
      dts = 0.5*dt / float(itmax)    ! NCEP assumed leap-frog is THE scheme

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
         enddo
      enddo

      do i=1,im
         precc(i) = 0.
         ktop(i) = 1
         kbot(i) = k1
      enddo

      do it=1, itmax
!
! Note: NCEP folks think everyone should be using leap-frog scheme
! therefore the input sdt is the actual time step * 0.5
!
      call sas(im, im, k1, delh, dts,
     &         qn, tn, ps, pn, sn, oro, wn,
     &         g, rair, rh2o, cpair, latvap, tmelt,
     &         cldwrk, rain, kb, kt, kuo, cldi, cldw)

        do i=1,im
           precc(i) = precc(i) + rain(i)
           ktop(i) = max(ktop(i), kt(i))
           kbot(i) = min(kbot(i), kb(i))
        enddo

      enddo

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
