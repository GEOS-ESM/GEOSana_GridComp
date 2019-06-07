#include <params.h>
      subroutine physdrv(abstot, absnxt, emstot, qrs, qrl,
     &                   diagattr, diagbuf,   tdt,  bdt,
     &                   nq,  u,  v,   pt,  q,
     &                   doalb, pk, pe, delp, ps, precp,
     &                   rayf, co2f, omga, peln, phis,  sgh,
     &                   pblht, tpert, qpert,  pkz, 
     &                   coslon, sinlon, nymd, nhms,
     &                   ih2o, h2odata, dcaf, gwdf, qbeg, tke,
     &                   jfirst, jlast, ng_d, ng_s, nslice)


#if defined( TIMING )
      use timingModule
#endif

#if defined( CLM2 )
! JDR Added clm2 module
      use atm_lndMod
#endif

      implicit none

      integer    imr, jnp, nl, nc
      parameter (imr = PLON)
      parameter (jnp = PLAT)
      parameter (nl  = PLEV)
      parameter (nc  = FVGCM_TRACER + PNATS)

      integer jfirst, jlast
      integer ng_d, ng_s

! To-dos: trim 2D arrays in comsrf.h to [jfirst:jlast]

#include <pmgrid.h>
#include <comcon.h>
#include <comsrf.h>
#include <commap.h>
#include <comctl.h>
#include <comtim.h>
#include <physbuf.h>

#if defined( CLM2 )
! JDR Added clm2 module for orbital parameters
#include <comsol.h>
#endif

#if defined( R4BYTE)
      real*4 abstot(imr,nl+1,nl+1,jfirst:jlast)  !Non-adjacent layer absorptivites
      real*4 absnxt(imr,nl,4,jfirst:jlast)       !Nearest layer absorptivities
      real*4 emstot(imr,nl+1,jfirst:jlast)       !Total emissivity
#else
      real   abstot(imr,nl+1,nl+1,jfirst:jlast)  !Non-adjacent layer absorptivites
      real   absnxt(imr,nl,4,jfirst:jlast)       !Nearest layer absorptivities
      real   emstot(imr,nl+1,jfirst:jlast)       !Total emissivity
#endif
      real   qrs(imr,nl,jfirst:jlast)    ! Shortwave heating rate
      real   qrl(imr,nl,jfirst:jlast)    ! Longwave  heating rate

#include <diag.h>

      real    tdt
      integer nq
      logical rayf
      logical doalb        ! do albedos in lsm
      logical dcaf         ! flag for dry convection
      logical gwdf         ! flag for gwd

      real   sinlon(imr)
      real   coslon(imr)
      integer nymd
      integer nhms
      integer ih2o
      character*(*) h2odata

      logical qbeg

! Ghosted prog arrays:
      real, intent(inout):: u(imr,jfirst-ng_d:jlast+ng_s,nl)  ! u-wind (m/s)
      real, intent(inout):: v(imr,jfirst-ng_s:jlast+ng_d,nl)  ! v-wind (m/s)
      real, intent(inout):: pt(imr,jfirst-ng_d:jlast+ng_d,nl)
      real delp(imr,jfirst:jlast,nl)
      real q(imr,jfirst-ng_d:jlast+ng_d,nl,nq)

      real pk(imr,jfirst:jlast,nl+1)
      real pe(imr,nl+1,jfirst:jlast)
      real peln(imr,nl+1,jfirst:jlast)
      real pkz(imr,jfirst:jlast,nl)
      real ps(imr,jfirst:jlast)

      real wspd(imr,nl)        ! wind speed needed by NCEP's SAS
      real tke(imr,jfirst:jlast,nl)

! 2D arrays to be trimmed.

      real phis(imr,jfirst:jlast)         !Sfc geopotential
      real  sgh(imr,jfirst:jlast)         !std. of topg
! CCM Phys:
      real pblht(imr,jnp)        !Planetary boundary layer height
      real tpert(imr,jnp)        !Thermal temperature excess
      real qpert(imr,pcnst,jnp)  ! Thermal humidity & constituent excess
      real precp(imr,jnp)

      integer    maxslice
      integer    nslice
!      parameter (maxslice    = pd2d + pd3d * nl)
      real*4     diagbuf(imr,nslice,jfirst:jlast) ! diagnostics buffer
      integer    diagattr(5,pdiag)           ! field attributes

! Local permanent
      real rfac(nl)        ! Rayleigh friction coefficient
      real wfac(nl)        ! Relexation coeff. for water vapor
      integer ks
      integer kw_top, kw_bot

! Local temporary
! JDC modify begin
      real wet1(imr,jnp)   ! top soil layer wetness (fraction)
      real wet(imr,jnp)    ! root zone soil wetness (fraction)
      real trad(imr,jnp)   !surface brightness temperature (K)
      real emsfc(imr,jnp)  !bulk surface emissivity (fraction)
! JDC modify end
      real q_ref(nl,jfirst:jlast)  ! UARS Ref H2O at model levels
      real ak(nl+1)
      real bk(nl+1)
      real ptop, pint
      integer i, j, k, ic
      real tq(imr)                 ! Column precipitable water (kg/m**2)
      real dt5

! Ghosted:
      real ua(imr,jfirst:jlast,nl)
      real va(imr,jfirst:jlast,nl)
      real psomc(imr,nl,jfirst:jlast)       ! 
      real  zm(imr,nl,jfirst:jlast)         ! geopotential height at mid level
      real thm(imr,nl,jfirst:jlast)         ! column potential temperature
      real tstar(imr, nl, jfirst:jlast)     ! 3D temperature
      real qstar(imr, nl,nc,jfirst:jlast)   ! adjusted q fields
      real omga(imr,nl,jfirst:jlast)         ! Vertical velocity (dp/dt)
      real pm(imr,nl)
      real pmln(imr,nl)

      real tvm(imr,nl)
      real zi(imr,nl+1)
      real t2(imr,nl)
      real fu(imr,nl)
      real fv(imr,nl)
      real ustar(imr,nl)
      real vstar(imr,nl)
      real pdel(imr,nl)
      real press(nl)

      real fdq(imr,nl)
      real c1
      real c2
      real pc
      real bdt

      real precsl(imr,jfirst:jlast) ! L.S. precipitation rate (snow)
      real precsc(imr,jfirst:jlast) ! C.S. precipitation rate (snow)

      real zm2(imr,nl)              ! Height of midpoints (above surface)
      real rayfdt(plond,plev)       ! dT/dt due to Rayleigh Friction
      real dcafdt(plond,plev)       ! dT/dt due to dry convective adjustment
      real diabdt(plond,plev)       ! dT/dt due to all physics parameterizations
! local
      real wk1(imr)
      real wk2(imr)
      real q2c(imr,nl,nc)           !2D  q adjusted
      real q2 (imr,nl)              !2D  specific humidity
      real pt2(imr,nl)

! work arrays for pft2d:
      real tmp1((imr+2)*(jlast-jfirst+1))
      real tmp2((imr+1)*(jlast-jfirst+1))

      real slp(imr)         ! sea-level-pressure
      real hyam,  hyam1,  hybm,  hybm1
      real prsfac
      real rgrav
      real rrg

      real r1, durf, dvrf
      integer ndt
      real rdt                 ! inverse of the time step
      real rcp                 ! inverse of cpair

      real qmin, qmax, vmax
      real pek
      real ::  undef = 1.e25	! undefined for pressure level diagnostics
      real ::  pout        	! output level log(pressure) for diagnostics

      integer nsecf, n

! Chern 01/2001
! diagnostics at 2 and 10 meter
! local
      real t10m(imr)           ! 10 meter temperature
      real q10m(imr)           ! 10 meter specify humidity
      real u10m(imr)           ! 10 meter U wind
      real v10m(imr)           ! 10 meter V wind
      real t2m(imr)            ! 2 meter temperature
      real q2m(imr)            ! 2 meter specify humidity
      real u2m(imr)            ! 2 meter U wind
      real v2m(imr)            ! 2 meter V wind

! Chern 01/2002
! diagnostics of bulk moisture availability
! local
      real bma(imr)            ! bulk moisture availability (fraction)

! SJL 00-07-11
! Polar filter related:
      real tmm(imr,jfirst:jlast,nl)     ! pt tendency
      logical filp             ! polar filter for phys tendencies
      integer  ifax(13)                      !ECMWF fft
      real, allocatable, save :: trigs(:)
      real, allocatable, save :: dc(:,:), sc(:)
      real, allocatable       :: de(:,:), se(:)

#if defined ( SPMD )
      real, allocatable, save :: uasouth(:,:)
#endif

      real  sinp(jnp)
      real  sine(jnp)
      real  cosp(jnp)
      real  cose(jnp)
      real dl, dp
      real qz

      logical h2o
      logical co2f
      logical filz
      logical first

      data h2o    /.true./
      data filz   /.true./
      data filp   /.true./
      data first  /.true./

      save press
      save rfac
      save wfac
      save ks
      save kw_top, kw_bot
      save prsfac
      save ifax
      integer js2g0, jn2g0, jn1g1

      nsecf(n)  = n/10000*3600 + mod(n,10000)/100* 60 + mod(n,100)

      js2g0 = max(2,jfirst)
      jn2g0 = min(jnp-1,jlast)
      jn1g1 = min(jnp,jlast+1)

      if ( imr <= 144 ) filp = .false.

! Set CO2 volume mixing ratio

      if ( co2f ) call set_co2(nymd)

      if(first) then

#if defined ( SPMD )
        allocate( uasouth(imr,nl) )
#endif

! Setup reference coordinate
         call set_eta(nl, ks, ptop, pint, ak, bk)

         if ( rayf ) then
! Initialization for Rayleigh friction
         if(masterproc) write(6,*) 'time step (seconds) for ray_f=', tdt

! e-folding time
         if( ak(1) <= 50. ) then
!           c1 = 1. / (5.*24*3600.)
! SJL 06/08/2001
            c1 = 1. / (6.*24*3600.)
         else
            c1 = 1. / (13.*24.*3600.)
         endif

         pc = max(10., ak(1))

      do k=1,ks
         press(k) = 0.5*(ak(k) + ak(k+1))
         rfac(k) = tdt*c1*(1.+tanh(1.5*log10(pc/press(k))))
         if(masterproc) write(6,*) k, press(k), rfac(k)
      enddo

      if(masterproc) write(6,*) ' '
      endif

      if( h2o ) then
         if(masterproc) write(6,*) 'Relaxation coef for H2o (Days) '
         if(masterproc) write(6,*) ' '
! Compute relexation coef for water vapor.
         c2 = tdt / (7.*24.*3600.)
         pc = 100.                   ! 1 mb

         kw_top = 1

      do k=1,ks
! Cutoff below 30 mb
         press(k) = 0.5*(ak(k) + ak(k+1))
         if(press(k) .lt. 3000. ) then
            wfac(k) = c2*(1.+tanh(log10(pc/press(k))))
            if(masterproc) write(6,*) k, press(k)/100., tdt/(wfac(k)*24.*3600.)
            kw_bot = k
         else
            wfac(k) = 0.
         endif
      enddo
      endif              ! end h2o check

! Constants for CCM SLP algorithm:
        hyam  = 0.5*(ak(nl) + ak(nl+1)) * 1.e-5
        hyam1 = 0.5*(ak(nl) + ak(nl-1)) * 1.e-5
        hybm  = 0.5*(bk(nl) + bk(nl+1))
        hybm1 = 0.5*(bk(nl) + bk(nl-1))
        prsfac = log(hyam + hybm) / log( (hyam  + hybm ) /
     &                                 (hyam1 + hybm1)   )

! polar filter initialization:
        if ( filp ) then
          allocate( trigs(3*imr/2+1) )
          allocate( sc(js2g0:jn2g0),    se(js2g0:jn1g1)    )
          allocate( dc(imr,js2g0:jn2g0), de(imr,js2g0:jn1g1) )

          call setrig(imr, jnp, dp, dl, cosp, cose, sinp, sine)

          call fftfax(imr, ifax, trigs)
          call pft_cf(imr, jnp, js2g0, jn2g0, jn1g1, sc, se, dc, de,
     &                cosp, cose, 45.)
          deallocate( se )
          deallocate( de )
        endif

        first = .false.
      endif

! Read in reference water vapor in the stratosphere
! monthly and zonal mean (data prepared by F. Verter from Randell)

      if( h2o ) then
         call haloe_h2o(ih2o, h2odata, nymd, nhms, kw_top, kw_bot,
     &                  press, q_ref, jfirst, jlast)
      endif

      call d2a3d(u, v, ua, va, imr, jnp, nl, jfirst, jlast,
     &           ng_d, ng_s, coslon, sinlon)

      call timing_on('TPHYSBC')

!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j,k,ic,c1,c2, wspd)
!$omp& private(t2, q2c,wk1,wk2,pdel)
!$omp& private(pmln, pm, tvm, zm2, qz)

      do 1000 j=jfirst,jlast

      do ic=1,nq
         do k=1,nl
            do i=1,imr
               q2c(i,k,ic) = q(i,j,k,ic)
            enddo
         enddo
      enddo

      do k=1,nl
         do i=1,imr
            pmln(i,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))
             tvm(i,k)   =  pt(i,j,k)
              t2(i,k)   =  tvm(i,k) / (1.+zvir*q2c(i,k,1))
            pdel(i,k)   = delp(i,j,k)
              pm(i,k) = 0.5*(pe(i,k,j) + pe(i,k+1,j))
#if defined (ALT_CPS)
            wspd(i,k) = sqrt(ua(i,j,k)**2 + va(i,j,k)**2)
#endif
          enddo
      enddo

!   Relaxation of stratosphere water vapor to UARS (HALOE/MSL) Reference
!   atmosphere: q_ref(k,j) defined at model pressure levels

      if(h2o) then
         if( latdeg(j) .ge. -82. .and. latdeg(j) .le. 82. ) then
         do k=kw_top, kw_bot
            if(q_ref(k,j) .lt. 1.e-5) then
               c1 = 1. / (1. + wfac(k))
               c2 = wfac(k)*q_ref(k,j)
! Compute zonal mean
!               qz = SUM(q2c(1:imr,k,1)) / imr
               qz = 0.0
               do i=1, imr
                 qz = qz + q2c(i,k,1)
               enddo
               qz = qz / imr
               do i=1, imr
!                 q2c(i,k,1) = c1 * (q2c(i,k,1) + c2)
! Apply h20 source/sink to zonal mean
                  q2c(i,k,1) = max(1.e-12, q2c(i,k,1) + wfac(k)*(q_ref(k,j)-qz))
               enddo
            endif
          enddo
          endif
      endif
 
      if(filz) call fillz(imr, nl, nq, q2c, pdel)


! Call NCAR CCM tendency physics BEFORE coupling.

      call tphysbc(
     &     j              ,clat(j)      ,bdt          ,diagbuf(1,1,j) ,
     &     oro(1,j)       ,phis(1,j)    ,pm           ,pe(1,1,j)      ,
     &     pdel           ,pmln         ,peln(1,1,j)  ,
     &     omga(1,1,j)    ,pblht(1,j)   ,pe(1,nl+1,j) ,tpert(1,j)     ,
     &     tvm            ,t2           ,q2c          ,ts(1,j)        ,
     &     qpert(1,1,j)   ,precl(1,j)   ,precc(1,j)   ,precsl(1,j)    ,
     &     precsc(1,j)    ,zm2          ,asdir(1,j)   ,asdif(1,j)     ,
     &     aldir(1,j)     ,aldif(1,j)   ,snowh(1,j)   ,qrs(1,1,j)     ,
     &     qrl(1,1,j)     ,flwds(1,j)   ,fsns(1,j)    ,lwup(1,j)      ,
     &     qstar(1,1,1,j) ,tstar(1,1,j) ,srfrad(1,j)  ,sols(1,j)      ,
     &     soll(1,j)      ,solsd(1,j)   ,solld(1,j)   ,precst(1,j)    ,
     &     absnxt(1,1,1,j),abstot(1,1,1,j) ,emstot(1,1,j),
     $     diagattr       ,wspd )

! Convert precip to mm/day
      do i=1, imr
          wk1(i) = precl(i,j) * 86400000.
          wk2(i) = precc(i,j) * 86400000.
          precp(i,j) = wk1(i) + wk2(i)
      enddo

      call outfld(diagattr(1,iOMEGA), omga(1,1,j), imr, j, 
     &            diagbuf(1,1,j))

1000  continue

      call timing_off('TPHYSBC')

!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j,k,tvm,zi, rrg, slp, pmln)

      do 3000 j=jfirst,jlast

      call virtem(tstar(1,1,j), qstar(1,1,1,j), zvir, tvm)

      call slp_das(imr, nl, pe(1,nl+1,j), phis(1,j), slp, pe(1,1,j),
     &             tvm, rair, gravit)

      call outfld(diagattr(1,iSLP), slp, imr, j, diagbuf(1,1,j))

      call zint(peln(1,1,j) ,pmln  ,rair    ,gravit  ,tvm,
     $          zm(1,1,j)   ,zi )

      rrg = rair / gravit

      do k=1, nl
        do i=1,imr
! SJL: finite-volume algorithm:
            zm(i,k,j) = zi(i,k+1) + rrg*tvm(i,k)*
     &          (1. - pe(i,k,j)*(peln(i,k+1,j)-peln(i,k,j))
     &            /(pe(i,k+1,j)-pe(i,k,j)) )
         psomc(i,k,j) = pk(i,j,nl+1) / pkz(i,j,k)
           thm(i,k,j) = tstar(i,k,j) * psomc(i,k,j)
        end do
      end do

! Stuff global fluxes and state variables into common for CCM

      do i=1,imr
         tbot(i,j) = tstar(i,nl,j)
         thbot(i,j) = thm(i,nl,j)
         zbot(i,j) = zm(i,nl,j)
         ubot(i,j) = ua(i,j,nl)
         vbot(i,j) = va(i,j,nl)
         qbot(i,j) = qstar(i,nl,1,j)
         pbot(i,j) = 0.5*(pe(i,nl+1,j) + pe(i,nl,j))
         wet1(i,j) = 1.0
         wet(i,j)  = 1.0
      enddo
3000  continue


! need to init nlend
!	if(nlend) return
! Call "surface models"
! 30 min. step for land surface model.

#if ( !defined AQUA )
#if defined( CLM2 )

! JDR Added call to clm2 driver

      call timing_on('atmlnd_drv')
      call atmlnd_drv(iradsw, eccen, obliqr,
     &                lambm0, mvelpp)

      call timing_off('atmlnd_drv')
!
! Over land set tssub to ts
!
      do j=jfirst,jlast
         do k=1,plevmx
            do i=1,imr
               if (nint(oro(i,j)).eq.1) tssub(i,k,j) = ts(i,j)
            end do
         end do
      end do

#else

      mcsec = nsecf(nhms)
      call timing_on('LSMMAIN')
      call lsmmain( imr,     jnp,     begj,    endj,    pcnst,
     $              iradsw,  nstep,   mbdate,  mdbase,  msbase,
     $              mbsec,   tdt,     pbot(1,begj),     thbot(1,begj),
     $              ubot(1,begj),     vbot(1,begj),      tbot(1,begj),
     $              qbot(1,begj),     zbot(1,begj),     flwds(1,begj),
     $             precc(1,begj),    precl(1,begj),
     $              sols(1,begj),     soll(1,begj),     solsd(1,begj),
     $             solld(1,begj),      oro(1,begj), 
     $               shf(1,begj),      lhf(1,begj),     cflx(1,1,begj),
     $               wsx(1,begj),      wsy(1,begj),         ts(1,begj), 
     $              tref(1,begj),    asdir(1,begj),      asdif(1,begj), 
     $             aldir(1,begj),    aldif(1,begj),      snowh(1,begj),
     $              lwup(1,begj),  tssub(1,1,begj),    plevmx,
     $              doalb,           nymd,              mcsec ,
     $              z0m(1,begj),       z0h(1,begj),        zpd(1,begj),
     $              wet1(1,begj),    wet(1,begj)  ,       trad(1,begj),
     $              emsfc(1,begj)    )
      call timing_off('LSMMAIN')

#endif
#endif

! Call ocean model to get SST/sea_ice
      ndt = nint ( tdt )
      call ccmoce (nymd, nhms, ndt)

      call timing_on('TPHYSAC')

!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j,k, pdel,ustar,vstar,dt5,fu,fv,t2, zi, rgrav)
!$omp& private(zm2,slp,r1,durf,dvrf,rcp,rdt, fdq, q2,tq,ic)
!$omp& private(pek, pt2, rayfdt, dcafdt, pm, pmln, pout)
!$omp& private(t10m,q10m,u10m,v10m,t2m,q2m,u2m,v2m,bma)

      do 4000 j=jfirst,jlast

      call outfld(diagattr(1,iLHFX), lhf(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iSHFX), shf(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iORO),  oro(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iQFLX), cflx(1,1,j), imr, j, 
     &            diagbuf(1,1,j))
      call outfld(diagattr(1,iSNOWH),snowh(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iTAUX), wsx(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iTAUY), wsy(1,j), imr, j, diagbuf(1,1,j))
! JDC modify begin
      call outfld(diagattr(1,iZ0M),  z0m(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iZ0H),  z0h(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iZPD),  zpd(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iGWET1),wet1(1,j),imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iGWET), wet(1,j), imr, j, diagbuf(1,1,j))
        do i = 1, imr
          if( oro(i,j) .ne. 1.0 ) then
            trad(i,j) = ts(i,j)
            emsfc(i,j) = 1.0
          end if
        end do
      call outfld(diagattr(1,iTSKIN),    ts(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iTRAD),   trad(1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iEMSFC), emsfc(1,j), imr, j, diagbuf(1,1,j))
!      call outfld(diagattr(1,iTREFHT),tref(1,j), imr, j, diagbuf(1,1,j))
! JDC modify end
!
! JDC: 2 and 10 meter diagnostics
!
      call diagref(                                                
     &     j,imr,                                                   
     &     pbot(1,j),   ubot(1,j),   vbot(1,j),   tbot(1,j),     
     &     qbot(1,j),   thbot(1,j),  zbot(1,j),                 
     &     wsx(1,j),    wsy(1,j),    shf(1,j),    lhf(1,j),    
     &     ts(1,j),     oro(1,j),
     &     t10m,        q10m,        u10m,        v10m,       
     &     t2m,         q2m,         u2m,         v2m )

      call outfld(diagattr(1,iT10M), t10m, imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iQ10M), q10m, imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iU10M), u10m, imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iV10M), v10m, imr, j, diagbuf(1,1,j))

      call outfld(diagattr(1,iT2M),  t2m,  imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iQ2M),  q2m,  imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iU2M),  u2m,  imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iV2M),  v2m,  imr, j, diagbuf(1,1,j))
!
! JDC: bulk moisture avalibility diagnostics
!
      call diagbma(
     &     j, imr,
     &     thbot(1,j),  qbot(1,j),  pe(1,nl+1,j),  ts(1,j),  oro(1,j),
     &     shf(1,j),    lhf(1,j),   wet1(1,j),      bma   )

      call outfld(diagattr(1,iBMA), bma, imr, j, diagbuf(1,1,j))

      do k=1,nl
         do i=1,imr
            pdel(i,k)  =   pe(i,k+1,j) - pe(i,k,j)
              pm(i,k)  = 0.5*(pe(i,k,j) + pe(i,k+1,j))
            ustar(i,k) =  ua(i,j,k)
            vstar(i,k) =  va(i,j,k)
               fu(i,k) =  0.
               fv(i,k) =  0.
               t2(i,k) =  0.
          enddo
      enddo

#if defined (ALT_PBL)
      do k=1,nl
         do i=1,imr
            zm2(i,k) = pkz(i,j,k)
             q2(i,k) = tke(i,j,k)
          enddo
      enddo

      do k=1,nl+1
         do i=1,imr
            zi(i,k) = pk(i,j,k)
          enddo
      enddo
#endif

      call tphysac(j  ,diagbuf(1,1,j),  tdt   ,pm ,pe(1,1,j)  ,
     $             pdel ,peln(1,1,j),   zm2,  zi,  q2,  t2,
     $             fu,  fv,  ustar   ,vstar   ,
     &             tstar(1,1,j)  ,qstar(1,1,1,j)  ,psomc(1,1,j),
     $             thm(1,1,j)   ,zm(1,1,j) ,pblht(1,j) ,qpert(1,1,j) ,
     $             tpert(1,j) ,clat(j)  ,shf(1,j) ,wsx(1,j),
     $             wsy(1,j)  ,cflx(1,1,j),sgh(1,j) ,lhf(1,j) ,
     $             oro(1,j), ts(1,j)   ,tref(1,j), precl(1,j),
     &             precc(1,j),  diagattr, gwdf, qbeg)


#if defined (ALT_PBL)
      do k=1,nl
         do i=1,imr
            tke(i,j,k) = q2(i,k)
          enddo
       enddo
#endif

! Rayleigh friction (implicit form)
         rdt = 1. / tdt
      if(rayf) then
         do k=1, nl
            do i=1,imr
              rayfdt(i,k) = tstar(i,k,j)
            enddo
         enddo

         do k=1, ks
         if(pm(1,k) .le. 3000.) then           ! only above 30 mb
            rcp = 1. / ( cpair * ( 1. - cappa*pe(1,k,j)*
     &                  (peln(1,k+1,j)-peln(1,k,j))/pdel(1,k) ) )
            do i=1,imr
               r1   =  rfac(k)
               r1   = -r1/(1.+r1)
              durf  =  ustar(i,k)*r1
              dvrf  =  vstar(i,k)*r1
         tstar(i,k,j) =  tstar(i,k,j) - (durf*(ustar(i,k)+0.5*durf)
     &                  + dvrf*(vstar(i,k)+0.5*dvrf) ) * rcp              
         ustar(i,k)  =  ustar(i,k) + durf
         vstar(i,k)  =  vstar(i,k) + dvrf
            fu(i,k)  =  fu(i,k) + durf * rdt
            fv(i,k)  =  fv(i,k) + dvrf * rdt
         enddo
         endif    ! P-check
         enddo    ! k-loop
         do k=1, nl
            do i=1,imr
              rayfdt(i,k) = (tstar(i,k,j) - rayfdt(i,k)) * rdt
            enddo
         enddo
      else if (diagattr(1,iRAYFDT) .eq. 1) then
         do k=1, nl
            do i=1,imr
              rayfdt(i,k) = 0.
            enddo
         enddo

      endif    ! Rayf check

      call outfld(diagattr(1,iRAYFDT), rayfdt, imr, j, diagbuf(1,1,j))



      do k=1,nl

         if( j == 1 .or. j == jnp ) then
                call xpavg( tstar(1,k,j), imr)
             do ic=1, nc
                call xpavg( qstar(1,k,ic,j), imr)
             enddo
         endif

!-----------------------------------------------------
! delp to be modified by gain/loss in water vapor.
!-----------------------------------------------------

         do i=1,imr
           fdq(i,k) = 1. + qstar(i,k,1,j) - q(i,j,k,1)
         enddo

         if( j == 1 .or. j == jnp ) then
             call xpavg( fdq(1,k),  imr)
         endif

      enddo

! Precipitable water
      do i=1,imr
         tq(i) = 0.
      enddo

      do k=1,nl
         do i=1,imr
           pdel(i,k)   = pdel(i,k) * fdq(i,k)
             q2(i,k)   = qstar(i,k,1,j) / fdq(i,k)
           tq(i)       = tq(i) + pdel(i,k)*q2(i,k)
           pe(i,k+1,j) = pe(i,k,j) + pdel(i,k)
         enddo
      enddo

! Convert to kg/m**2
      do i=1,imr
         tq(i) = tq(i) * rga
      enddo

      call outfld(diagattr(1,iTQ), tq, imr, j, diagbuf(1,1,j))

! zi is 2D temp storage for pk
! t2 is 2D temp storage for pkz

      pek = pk(1,j,1)
      do i=1,imr
          zi(i,1) = pek
      enddo

! update peln and pk
      do k=2,nl+1
         do i=1,imr
            peln(i,k,j) =  log(pe(i,k,j))
              zi(i,k)   =  pe(i,k,j)**cappa
         enddo
      enddo

! update pkz
      do k=1,nl
         do i=1,imr
            t2(i,k) = (zi(i,k+1) - zi(i,k) )  /
     &                (cappa*(peln(i,k+1,j) - peln(i,k,j)) )
           pt2(i,k) = tstar(i,k,j) * (1. + zvir*q2(i,k)) / t2(i,k)
           pmln(i,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))     ! needed for 2d pressure diag
         enddo
      enddo

!-----------------------------------------------------
! Perform dry convective adjustment
!-----------------------------------------------------

      if(dcaf) then

        call dry_adj(imr,  nl, rdt, pt2, fu, fv, 
     &               ustar, vstar, pdel, j)

        do k=1,nl
           do i=1,imr
             dcafdt(i,k) = (pt2(i,k)*t2(i,k)/(1.+zvir*q2(i,k))
     &                     - tstar(i,k,j)) * rdt
          end do
        end do
      else if (diagattr(1,iDCAFDT) .eq. 1) then
        do k=1,nl
           do i=1,imr
             dcafdt(i,k) = 0.
          end do
        end do
      endif

      call outfld(diagattr(1,iDCAFDT),dcafdt ,imr,j,diagbuf(1,1,j))


      do k=1,nl
         do i=1,imr
! Copy (u, v) tendency to 3D array
            ua(i,j,k)   = fu(i,k)
            va(i,j,k)   = fv(i,k)
            q(i,j,k,1) = q2(i,k)
            pt(i,j,k)   = pt(i,j,k) / pkz(i,j,k)
            tmm(i,j,k)   = pt2(i,k) - pt(i,j,k)
            pkz(i,j,k)   = t2(i,k)
            pk(i,j,k+1) = zi(i,k+1)
            delp(i,j,k)   = pdel(i,k)
        enddo
      enddo

      if( nq .gt. 1 ) then
        do ic=2,nq
           do k=1,nl
              do i=1,imr
                 q(i,j,k,ic) = qstar(i,k,ic,j) / fdq(i,k)
              enddo
           enddo
        enddo
      endif

      do i=1,imr
         ps(i,j) = pe(i,nl+1,j)
      enddo

      dt5 = 0.5*tdt

! Sea-Level-Pressure (SLP): CCM3 Algorithm:
!     call ccmslp(prsfac, phis(1,j) ,ps(1,j)   ,tstar(1,nl-1,j)
!    &           ,tstar(1,nl,j)    ,slp   )
!     call outfld(diagattr(1,iSLP), slp, imr, j, diagbuf(1,1,j))

      call outfld(diagattr(1,iSURFP), ps(1,j), imr, j, diagbuf(1,1,j))

      call outfld(diagattr(1,iU), ustar, imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iV), vstar, imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iT), tstar(1,1,j), imr, j, diagbuf(1,1,j))
      call outfld(diagattr(1,iQ), q2,    imr, j, diagbuf(1,1,j))

! Compute height by adding surface height to zm
! Note: zm is not really up to date.

      rgrav = 1./ gravit

         do i=1,imr
            slp(i) = phis(i,j) * rgrav
         enddo

      do k=1, nl
         do i=1,imr
            zm2(i,k) = slp(i) + zm(i,k,j)
         enddo
      enddo
      call outfld(diagattr(1,iH), zm2  ,   imr, j, diagbuf(1,1,j))

!  Compute 2d pressure level diagnostics, use tq for temporary variable
      if (diagattr(1,iU200) .eq. 1) then
        pout = log(20000.)
        call linint(imr,nl,pmln,ustar,undef,pout,tq)
        call outfld(diagattr(1,iU200), tq, imr, j, diagbuf(1,1,j))
      endif
      if (diagattr(1,iU850) .eq. 1) then
        pout = log(85000.)
        call linint(imr,nl,pmln,ustar,undef,pout,tq)
        call outfld(diagattr(1,iU850), tq, imr, j, diagbuf(1,1,j))
      endif
      if (diagattr(1,iV200) .eq. 1) then
        pout = log(20000.)
        call linint(imr,nl,pmln,vstar,undef,pout,tq)
        call outfld(diagattr(1,iV200), tq, imr, j, diagbuf(1,1,j))
      endif
      if (diagattr(1,iV850) .eq. 1) then
        pout = log(85000.)
        call linint(imr,nl,pmln,vstar,undef,pout,tq)
        call outfld(diagattr(1,iV850), tq, imr, j, diagbuf(1,1,j))
      endif
      if (diagattr(1,iT200) .eq. 1) then
        pout = log(20000.)
        call linint(imr,nl,pmln,tstar(1,1,j),undef,pout,tq)
        call outfld(diagattr(1,iT200), tq, imr, j, diagbuf(1,1,j))
      endif
      if (diagattr(1,iT850) .eq. 1) then
        pout = log(85000.)
        call linint(imr,nl,pmln,tstar(1,1,j),undef,pout,tq)
        call outfld(diagattr(1,iT850), tq, imr, j, diagbuf(1,1,j))
      endif

4000  continue
      call timing_off('TPHYSAC')

	dt5 = 0.5*tdt

!$omp  parallel do
!$omp& default(shared)
!$omp& private(i, j, k, tmp1, tmp2)

! Adjust D-grid  u-winds
      do 5000 k=1,nl

! Polar filter diabatic heating tendency: tmm
       if(filp) call pft2d(tmm(1,js2g0,k), sc(js2g0), dc(1,js2g0),
     &                     imr, jn2g0-js2g0+1, ifax, trigs, tmp1, tmp2)

      do j=jfirst,jlast
         do i=1,imr
            pt(i,j,k) = pt(i,j,k) + tmm(i,j,k)
         enddo
      enddo

!-----------------------------------------------------
! Pole averaging of delp and pt is needed to prevent PSAS from creating 
! rounding differences at the poles.
!-----------------------------------------------------

         call polavg(delp(1,jfirst,k), imr, jnp, jfirst, jlast)
         call polavg(  pt(1,jfirst,k), imr, jnp, jfirst, jlast)

#if ( !defined ALT_PFT )
! Filter physics tendencies
       if(filp) then
          call pft2d(ua(1,js2g0,k), sc(js2g0), dc(1,js2g0),
     &               imr, jn2g0-js2g0+1, ifax, trigs, tmp1, tmp2)
          call pft2d(va(1,js2g0,k), sc(js2g0), dc(1,js2g0),
     &               imr, jn2g0-js2g0+1, ifax, trigs, tmp1, tmp2)
       endif
#endif


! Compute u at jfirst later; need ua(j=jfirst-1)
      do j=jfirst+1,jlast
         do i=1,imr
            u(i,j,k) = u(i,j,k) + dt5*(ua(i,j,k)+ua(i,j-1,k))
         enddo
      enddo

! Adjust D-grid v-winds
      do j=js2g0,jn2g0
            v(1,j,k) = v(1,j,k) + dt5*(va(1,j,k)+va(imr,j,k))
         do i=2,imr
            v(i,j,k) = v(i,j,k) + dt5*(va(i,j,k)+va(i-1,j,k))
         enddo
      enddo
5000  continue

#if defined ( SPMD )
      call mp_send_ua(imr, jnp, jfirst, jlast, 1, nl, ua)
#endif

      rdt = 1. / tdt

      if ((diagattr(1,iH300) .eq. 1 .or. diagattr(1,iH500) .eq. 1) .or.
     &    (diagattr(1,iDIABDT) .eq. 1)) then

!$omp parallel do private(i,j,k, zi, t2, diabdt) 
!$omp& private(tvm, tq, pmln, rgrav, pout)

      do j=jfirst,jlast

        do k=1,nl
          do i=1,imr
            t2(i,k) = (pk(i,j,k+1) - pk(i,j,k) )  /
     &                (cappa*(peln(i,k+1,j) - peln(i,k,j)) )
            diabdt(i,k) = tmm(i,j,k) * t2(i,k) /
     &                   (1. + zvir*q(i,j,k,1)) * rdt
            tvm(i,k) = pt(i,j,k) * t2(i,k) 
            pmln(i,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))  
          enddo
        enddo
        call outfld(diagattr(1,iDIABDT), diabdt, imr, j, diagbuf(1,1,j))
        do k=1, nl+1
           do i=1,imr
              zi(i,k) = 0.
           enddo
        enddo
        call zint(peln(1,1,j) ,pmln  ,rair    ,gravit  ,tvm,
     $            zm(1,1,j)   ,zi )
        rgrav = 1./ gravit
        do k=1, nl+1
           do i=1,imr
              zi(i,k) = zi(i,k) + phis(i,j) * rgrav
           enddo
        enddo
        call outfld(diagattr(1,iDIABDT), diabdt, imr, j, diagbuf(1,1,j))
        pout = log(50000.)
        call linint(imr,nl+1,peln(1,1,j),zi,undef,pout,tq)
        call outfld(diagattr(1,iH500), tq, imr, j, diagbuf(1,1,j))
        pout = log(30000.)
        call linint(imr,nl+1,peln(1,1,j),zi,undef,pout,tq)
        call outfld(diagattr(1,iH300), tq, imr, j, diagbuf(1,1,j))
      enddo
      endif

#if defined ( SPMD )
      call mp_barrier
      call mp_recv_ua(imr, jnp, jfirst, jlast, 1, nl, uasouth)

      if ( jfirst > 1 ) then
!$omp parallel do private(i,k)
      do k=1,nl
         do i=1,imr
            u(i,jfirst,k) = u(i,jfirst,k) + 
     &                     dt5*(ua(i,jfirst,k)+uasouth(i,k))
         enddo
      enddo
      endif
#endif

#if defined (ALT_PBL)
      qbeg = .false.
#endif

      return
      end

      subroutine set_co2(nymd)
      implicit none
      integer  nymd

#include <crdcon.h>

      integer  year

       integer NYR

!---------------------- GISS data from 1950 to 2000 -----------
!      parameter ( NYR = 11 )
!      data co2 / 310.7,  313.1,  316.5,  319.6,  324.7,  330.2,
!    &            337.9,  344.9,  353.0,  359.5,  368.4 /
!---------------------- GISS data -----------------------------

! CAM1 data          ! BAU ppmv from 1870 to 2100
       parameter ( NYR = 231 )
       real co2(NYR)
       data co2 /
     & 289.263, 289.416, 289.577, 289.745, 289.919,
     & 290.102, 290.293, 290.491, 290.696, 290.909,
     & 291.129, 291.355, 291.587, 291.824, 292.066,
     & 292.313, 292.563, 292.815, 293.071, 293.328,
     & 293.586, 293.843, 294.098, 294.350, 294.598,
     & 294.842, 295.082, 295.320, 295.558, 295.797,
     & 296.038, 296.284, 296.535, 296.794, 297.062,
     & 297.338, 297.620, 297.910, 298.204, 298.504,
     & 298.806, 299.111, 299.419, 299.729, 300.040,
     & 300.352, 300.666, 300.980, 301.294, 301.608,
     & 301.923, 302.237, 302.551, 302.863, 303.172,
     & 303.478, 303.779, 304.075, 304.366, 304.651,
     & 304.930, 305.206, 305.478, 305.746, 306.013,
     & 306.280, 306.546, 306.815, 307.087, 307.365,
     & 307.650, 307.943, 308.246, 308.560, 308.887,
     & 309.228, 309.584, 309.956, 310.344, 310.749,
     & 311.172, 311.614, 312.077, 312.561, 313.068,
     & 313.599, 314.154, 314.737, 315.347, 315.984,
     & 316.646, 317.328, 318.026, 318.742, 319.489,
     & 320.282, 321.133, 322.045, 323.021, 324.060,
     & 325.155, 326.299, 327.484, 328.698, 329.933,
     & 331.194, 332.499, 333.854, 335.254, 336.690,
     & 338.150, 339.628, 341.125, 342.650, 344.206,
     & 345.797, 347.397, 348.980, 350.551, 352.100,
     & 353.636, 355.197, 356.755, 358.286, 359.850,
     & 361.499, 363.231, 365.026, 366.880, 368.792,
     & 370.763, 372.759, 374.753, 376.750, 378.751,
     & 380.759, 382.788, 384.852, 386.947, 389.076,
     & 391.237, 393.441, 395.696, 398.001, 400.357,
     & 402.761, 406.186, 408.635, 411.112, 413.616,
     & 416.147, 418.703, 421.281, 423.883, 426.507,
     & 429.154, 431.823, 434.514, 437.226, 439.961,
     & 442.717, 445.496, 448.297, 451.122, 453.969,
     & 456.839, 459.745, 462.697, 465.695, 468.739,
     & 471.827, 474.960, 478.138, 481.361, 484.630,
     & 487.946, 491.308, 494.716, 498.172, 501.675,
     & 505.226, 508.812, 512.425, 516.065, 519.735,
     & 523.434, 527.163, 530.922, 534.712, 538.534,
     & 542.388, 546.273, 550.191, 554.142, 558.126,
     & 562.144, 566.179, 570.217, 574.260, 578.309,
     & 582.365, 586.429, 590.502, 594.584, 598.676,
     & 602.779, 606.896, 611.033, 615.189, 619.364,
     & 623.558, 627.769, 631.993, 636.231, 640.484,
     & 644.753, 649.036, 653.336, 657.653, 661.986,
     & 666.337, 670.705, 675.091, 679.496, 683.918,
     & 688.360, 692.820, 697.300, 701.800, 706.320,
     & 710.859 /

       year = nymd / 10000 - 1869
       year = max(year, 1)
       year = min(year, NYR)

      co2vmr = co2(year) * 1.e-6

!     write(*,*) 'Setting CO2 volume mixing ratio for year', nymd/10000,' to ', co2vmr

      return
      end
