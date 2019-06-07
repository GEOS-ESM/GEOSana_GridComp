!-----------------------------------------------------------------------
!BOP
! !ROUTINE: refout --- Reference output
!
! !INTERFACE:
      subroutine refout(im,    jm,     km,     jfirst, jlast, ccm_on,    &
                        iavg,  pe,     delp,   u,      v,     pt,        &
                        q,     hs,     pkz,    peln,   precp, nq,        &
                        rair,  coslon, sinlon, ae,     iout,  nymd,      &
                        nhms,  omga,   psx,    ux,     vx,    tx,        &
                        qx,    wx,     pvx,    usus,   usvs,  tsts,      &
                        vsts,  vspv,   zstat,  tg,     zvir,  grav,      &
                        tke,   ng_s,   ng_d,   oflnm,  undef, idt,       &
                        vcoord)

#if defined( SPMD )
      use mod_comm, only : gid, mp_gather4d
#define  CPP_SLP   slptmp
#define  CPP_HS    hstmp
#define  CPP_PS    pstmp
#define  CPP_TG    tgtmp
#define  CPP_PRECP precptmp
#define  CPP_EPV   epvtmp
#define  CPP_WZ    wztmp
#define  CPP_UA    uatmp
#define  CPP_VA    vatmp
#define  CPP_TA    tatmp
#define  CPP_TKE   tketmp
#else
#define  CPP_SLP   slp
#define  CPP_HS    hs
#define  CPP_PS    ps
#define  CPP_TG    tg
#define  CPP_PRECP precp
#define  CPP_EPV   epv
#define  CPP_WZ    wz
#define  CPP_UA    ua
#define  CPP_VA    va
#define  CPP_TA    ta
#define  CPP_TKE   tke
#endif

#if defined (GFIO)
      use m_die, only: die
#endif
      use precision
      use gmean, only: gmean_of

      implicit none

! !INPUT PARAMETERS:
      integer im                ! dimension in east-west
      integer jm                ! dimension in North-South
      integer km                ! number of Lagrangian layers
      integer jfirst            ! starting latitude index for MPI
      integer jlast             ! ending latitude index for MPI
      integer ng_s, ng_d
      logical ccm_on            ! running with CCM physics?
      integer nq                ! number of tracers
      integer iout              ! output unit number
      integer nymd, nhms        ! date and time
      integer iavg
      real(r8) rair
      real(r8) ae
      real(r8) zvir
      real(r8) grav
      logical zstat

      real(r8) u(im,jfirst-ng_d:jlast+ng_s,km)   ! u-wind (m/s)
      real(r8) v(im,jfirst-ng_s:jlast+ng_d,km)   ! v-wind (m/s)
      real(r8) hs(im,jfirst:jlast)               ! surface geopotential (grav*zs)

      real(r8) pt(im,jfirst-ng_d:jlast+ng_d,km)    ! virtual  potential temperature = T*/pkz
      real(r8) pkz(im,jfirst:jlast,km)   ! layer-mean value of pk
      real(r8) q(im,jfirst-ng_d:jlast+ng_d,km,*)   ! scalar tracers (e.g., specific humidity)
      real(r8) peln(im,km+1,jfirst:jlast)! log pressure (pe) at layer edges
      real(r8) delp(im,jfirst:jlast,km)  ! pressure thickness (pascal)
      real(r8) pe(im,km+1,jfirst:jlast)  ! pressure (pascal unit) at layer edges
      real(r8) precp(im,jm)
      real(r8) omga(im,km,jfirst:jlast)
      real(r8) psx(jfirst:jlast)
      real(r8) ux(jfirst:jlast,km)
      real(r8) vx(jfirst:jlast,km)
      real(r8) qx(jfirst:jlast,km)
      real(r8) tx(jfirst:jlast,km)
      real(r8) wx(jfirst:jlast,km)
      real(r8) pvx(jfirst:jlast,km)
      real(r8) usus(jfirst:jlast,km)
      real(r8) usvs(jfirst:jlast,km)
      real(r8) vsts(jfirst:jlast,km)
      real(r8) vspv(jfirst:jlast,km)
      real(r8) tsts(jfirst:jlast,km)
      real(r8) tg(im,jm)

! FastOpt begin
#if defined (ALT_PBL)
      real(r8) tke(im,jfirst:jlast,km)   ! PBL
#else
      real(r8) tke(1,1,1)   ! PBL
#endif
! FastOpt end

      real(r8) coslon(im)
      real(r8) sinlon(im)

! !DESCRIPTION:
!    Author: Shian-Jiann Lin NASA/GSFC
!
! !REVISION HISTORY:
!   SJL 99.04.13:  Delivery
!   WS  99.11.03:  Documentation; indentation
!   SJL 99.12.31:  Implementing parallel routine pmaxmin
!   WS  00.07.08:  Use precision module
!   WS  00.07.10:  Incorporate simplifications to PILGRIM
!
!EOP
!-----------------------------------------------------------------------
!BOC

      integer  i, j, k, ic, ks, n
      real(r8) epv(im,jfirst:jlast,km)
      real(r8)  ua(im,jfirst:jlast,km)
      real(r8)  va(im,jfirst:jlast,km)
      real(r8)  ps(im,jfirst:jlast)
      real(r8)  slp(im,jfirst:jlast)
      real(r8) ta(im,jfirst:jlast,km)
      real(r8) wz(im,jfirst:jlast,km+1)
      real(r8) ax(jfirst:jlast)

#if defined( SPMD )
      real(r8) slptmp(im,jm)
      real (r8) hstmp(im,jm)
      real (r8) pstmp(im,jm)
      real (r8) tgtmp(im,jm)
      real (r8) precptmp(im,jm)
      real (r8) epvtmp(im,jm,km)
      real (r8) wztmp(im,jm,km)		! do not need km+1 gather only for km
      real (r8) uatmp(im,jm,km)
      real (r8) vatmp(im,jm,km)
      real (r8) tatmp(im,jm,km)
      real (r8) tketmp(im,jm,km)
#endif

      real(r8) us(im,jfirst:jlast),vs(im,jfirst:jlast)
      real(r8) ts(im,jfirst:jlast),pvs(im,jfirst:jlast)
      real(r8) uu(im,jfirst:jlast),uv(im,jfirst:jlast)
      real(r8) tt(im,jfirst:jlast),vt(im,jfirst:jlast)
      real(r8) vpv(im,jfirst:jlast)
      real(r8) qs(im,jfirst:jlast)

      real(r8) tv(im,km,jfirst:jlast)
      real(r8) dl, dp
      real(r8) sinp(jm)
      real(r8) spec(im/2+1)
      real(r8) qmax, qmin, fac
      real(r8) wtmp

      real(r8), allocatable, save :: cosp(:), cose(:), acosp(:), sine(:)

! Polar filter related

! Work arrays for pft2d:
      real(r8) wk1( (im+2)*(jlast-jfirst+1) )
      real(r8) wk2( (im+1)*(jlast-jfirst+1) )
      integer  ifax(13)                      !ECMWF fft
      real(r8), allocatable, save :: trigs(:)
      real(r8), allocatable, save :: dc(:,:), sc(:)
      real(r8) se(jm)
      real(r8) de(im,jm)
      real(r8) ycrit
      real(r8) ptmp
      real(r8) rgf
      real(r8) ztop
      integer js2g0
      integer jn2g0
      integer jn1g1
      integer nrec
      logical filter

      logical first
      data first /.true./
      data filter /.true./
      save dl, dp
      save ifax, nrec

! HDF related variables
      integer         tVars
      parameter    ( tVars = 12 )
      character*128   :: source = 'Data Assimilation Office, NASA/GSFC'
      character*128   ::   contact = 'data@dao.gsfc.nasa.gov'
      character*128   ::   oTitle = 'FVGCM Dynamics State Vector'
      character(len=*), parameter :: myname = 'refout'
      integer         ::   out_fmode = 0        ! 0 for READ-WRITE
      integer              out_fid              ! output file ID
      integer              rc                   ! return error code
      integer         ::   outPrec = 0          ! Output file precision:
                                                ! 0 = 32 bits,  1 = 64bits
      character*(*)      oflnm                  ! output file name
      real          vcoord(*)                   ! vertical levels
      real          undef
      integer       idt, inc                    ! time step and time increment
      character*128 lName(tVars)                ! long name
      character*128 oName(tVars)                ! output variable name
      character*128 oUnit(tVars)                ! output unit
      integer       oVdim(tVars)                ! number of vertical levels
      real          valid_range(2, tVars)       ! range values for HDF output
      real          packing_range(2, tVars)     ! packing range for HDF output
      real lon(im), lat(jm)
      real          wk3d(im,jm,km)		! 3d space to gather if spmd
      real          wk2d(im,jm)			! 2d space to gather if spmd

#if ( !defined SPMD )
      integer gid
      gid = 0
#endif

      do j = 1, jm
         lat(j) = -90 + 180. / real(jm - 1) * (j-1)
      enddo

      do i = 1, im
         lon(i) = 360. / real(im) * (i-1)
      enddo
      inc  = 10000 * (idt/60/60) + 100 * mod(idt/60, 60) + mod(idt, 60)

! no TKE in hdf output yet
      lName(1) = 'Geopotential height at top edge of model layer'
      lName(2) = 'Surface geopotential'
      lName(3) = 'Sea-level pressure'
      lName(4) = 'Surface pressure'
      lName(5) = 'Zonal wind'    
      lName(6) = 'Meridional wind'
      lName(7) = 'Temperature'
      lName(8) = 'Ertels Potential Vorticity'
      lName(9) = 'Specific humidity'
      lName(10) = 'Surface temperature'
      lName(11) = 'Instaneous total precipitation'
      lName(12) = 'Vertical pressure velocity'

      oName(1) = 'ze'
      oName(2) = 'zs'
      oName(3) = 'slp'
      oName(4) = 'ps'
      oName(5) = 'u'
      oName(6) = 'v'
      oName(7) = 't'
      oName(8) = 'pv'
      oName(9) = 'q'
      oName(10) = 'tg'
      oName(11) = 'precp'
      oName(12) = 'omega'

      oUnit(1) = '(m/s)2'
      oUnit(2) = '(m/s)2'
      oUnit(3) = 'Pa'
      oUnit(4) = 'Pa'
      oUnit(5) = 'm/s'
      oUnit(6) = 'm/s'
      oUnit(7) = 'K'
      oUnit(8) = 'm2/(kg*sec)'
      oUnit(9) = 'kg/kg'
      oUnit(10) = 'K'
      oUnit(11) = 'mm/day'
      oUnit(12) = 'Pa/sec'

      oVdim(1) = km
      oVdim(2) = 0
      oVdim(3) = 0
      oVdim(4) = 0
      oVdim(5) = km
      oVdim(6) = km
      oVdim(7) = km
      oVdim(8) = km
      oVdim(9) = km
      oVdim(10) = 0
      oVdim(11) = 0
      oVdim(12) = km

      do i = 1, tVars
           valid_range(1,i) = undef
           valid_range(2,i) = undef
           packing_range(1,i) = undef
           packing_range(2,i) = undef
       enddo

        js2g0  = max(2,jfirst)
        jn2g0  = min(jm-1,jlast)
        jn1g1  = min(jm,jlast+1)

         rgf = 0.001 / grav

      if(first) then
        allocate( cosp(jm), cose(jm), acosp(jm), sine(jm) )
        call setrig(im,jm,dp,dl,cosp,cose,sinp,sine)

        do j=2,jm-1
          acosp(j) = 1./ cosp(j)
        enddo

! Initialization of polar filter

        allocate( trigs(3*im/2+1) )
        allocate( sc(js2g0:jn2g0) )
        allocate( dc(im,js2g0:jn2g0) )

        call fftfax(im, ifax, trigs)

        ycrit = 60.                      ! critical latitude for starting polar filter

        call pft_cf(im, jm, js2g0, jn2g0, jn1g1, sc, se, dc, de,  &
                    cosp, cose, ycrit)

!-----------------------------------
! Open direct access file for refout
!-----------------------------------
       if ( gid == 0 ) then
       open (iout, file='fort.82', form='unformatted',           &
             status='unknown', access='direct',recl=im*jm*4)
       nrec = 0
       endif
       first = .false.
      endif

      call d2a3d(u,v,ua,va,im,jm,km,jfirst,jlast,ng_d,ng_s,coslon,sinlon)

      if(zstat) then
        iavg = iavg + 1

! Compute EPV directly on the D-grid
        call epvd(im, jm, km, jfirst, jlast, u, v, pt, delp, epv, ng_s, ng_d)

      endif

!$omp  parallel do default(shared)     &
!$omp  private(i,j,k,ax,us,vs,ts,pvs,uu,uv,vt,tt,vpv, wtmp, wk1, wk2)

      do 1000 k=1,km

        if ( ccm_on ) then
          do j=jfirst,jlast
            do i=1,im
              tv(i,k,j) = pt(i,j,k)*pkz(i,j,k)
              ta(i,j,k) = tv(i,k,j)/(1.+zvir*q(i,j,k,1))
            enddo
          enddo
        else
          do j=jfirst,jlast
            do i=1,im
              tv(i,k,j) = pt(i,j,k)*pkz(i,j,k)
              ta(i,j,k) = tv(i,k,j)
            enddo
          enddo
        endif

        if(zstat) then

           if(filter) call pft2d(epv(1,js2g0,k), sc(js2g0), dc(1,js2g0),  &
                                 im, jn2g0-js2g0+1, ifax, trigs, wk1, wk2)
! Compute zonal means, zonal mean eddies and variances
          call zsmean(ua(1,jfirst,k),im,jfirst,jlast,us,ax)
          do j=jfirst,jlast
            ux(j,k) = ux(j,k) + ax(j)
          enddo

          call zsmean(va(1,jfirst,k),im,jfirst,jlast,vs,ax)
          do j=jfirst,jlast
            vx(j,k) = vx(j,k) + ax(j)
          enddo

          call zsmean(ta(1,jfirst,k),im,jfirst,jlast,ts,ax)
          do j=jfirst,jlast
            tx(j,k) = tx(j,k) + ax(j)
          enddo

          do j=jfirst,jlast
               wtmp = omga(1,k,j)
             do i=2,im
                wtmp = wtmp + omga(i,k,j)
             enddo
             wx(j,k) = wx(j,k) + wtmp/im
          enddo

          call zsmean(q(1,jfirst,k,1),im,jfirst,jlast,qs,ax)

          do j=jfirst,jlast
            qx(j,k) = qx(j,k) + ax(j)
          enddo

!         call avgp2(epv(1,jfirst,k),im,jm,jfirst,jlast)
          call zsmean(epv(1,jfirst,k),im,jfirst,jlast,pvs,ax)
          do j=jfirst,jlast
            pvx(j,k) = pvx(j,k) + ax(j)
          enddo

          do j=jfirst,jlast
	    do i=1,im
! U
	      uu(i,j) = us(i,j)**2
              uv(i,j) = us(i,j)*vs(i,j)
! T
	      tt(i,j) = ts(i,j)**2
              vt(i,j) = vs(i,j)*ts(i,j)
! meridional eddy PV transport
              vpv(i,j) = vs(i,j)*pvs(i,j)
            enddo
          enddo

          call zmean(uu,im,jm,jfirst,jlast,ax)
          do j=jfirst,jlast
            usus(j,k) = usus(j,k) + ax(j)
          enddo

          call zmean(uv,im,jm,jfirst,jlast,ax)
          do j=jfirst,jlast
            usvs(j,k) = usvs(j,k) + ax(j)
          enddo

          call zmean(tt,im,jm,jfirst,jlast,ax)
          do j=jfirst,jlast
            tsts(j,k) = tsts(j,k) + ax(j)
          enddo

          call zmean(vt,im,jm,jfirst,jlast,ax)
          do j=jfirst,jlast
            vsts(j,k) = vsts(j,k) + ax(j)
          enddo

          call zmean(vpv,im,jm,jfirst,jlast,ax)
          do j=jfirst,jlast
            vspv(j,k) = vspv(j,k) + ax(j)
          enddo
        endif
 1000  continue

!$omp  parallel do          &
!$omp  default(shared)      &
!$omp  private(i,j,k,ptmp)

      do 2000 j=jfirst,jlast
        do i=1,im
          ps(i,j) = pe(i,km+1,j)
        enddo

          ptmp = 0.
        do i=1,im
          ptmp = ptmp + ps(i,j)
        enddo
        psx(j) = psx(j) + ptmp/im

        call slp_das(im, km, ps(1,j), hs(1,j), slp(1,j), pe(1,1,j), tv(1,1,j),  &
                     rair,  grav)

! Compute geop. height
        do i=1,im
          wz(i,j,km+1) = hs(i,j) 
        enddo

        do k=km,1,-1
          do i=1,im
            wz(i,j,k) = wz(i,j,k+1) +          &
                        rair*tv(i,k,j)*(peln(i,k+1,j)-peln(i,k,j))
          enddo
        enddo
2000  continue

#if defined (GFIO)
!     If output file is no prog.bin, try to open a HDF output file.
!     If the HDF output file doesn't exist, create a new one.
      if ( index(oflnm, 'prog.bin') .le. 0 ) then
         call GFIO_Open(oflnm, out_fmode, out_fid, rc)
         if ( rc /= 0 )  then
            write(6,*) ' calling gfio create for ',trim(oflnm)
            call GFIO_Create ( oflnm, oTitle, source, contact, undef,          &
                        im, jm, oVdim(1), lon, lat, vcoord, "hPa",             &
                          nymd, nhms, inc, tVars, oName, lName, oUnit,         &
                          oVdim, valid_range, packing_range, outPrec,          &
                          out_fid, rc )
           if (rc /= 0) call die (myname,'wrong in GFIO_Create or GFIO_Open')
         endif
      end if
#endif

! Binary output
      if ( index(oflnm, 'prog.bin') .gt. 0 ) then
        call wrt3d(iout,nrec,im,jm,km+1,wz, jfirst, jlast)
        call wrt3d(iout,nrec,im,jm,1,slp, jfirst, jlast)
        call wrt3d(iout,nrec,im,jm,1,ps, jfirst, jlast)
        call wrt3d(iout,nrec,im,jm,km,ua, jfirst, jlast)
        call wrt3d(iout,nrec,im,jm,km,va, jfirst, jlast)
        call wrt3d(iout,nrec,im,jm,km,ta, jfirst, jlast)
#if defined (ALT_PBL)
        call wrt3d(iout,nrec,im,jm,km,tke, jfirst, jlast)
#endif

#if defined (GFIO)
! HDF output
      else

#if defined( SPMD )
         call mp_gather4d(wz, wztmp, im, jm, km, 1, jfirst, jlast,           &
                          1,  km,    0,  0,  0)
#endif
         if (gid == 0) then
         call GFIO_PutVar (out_fid, oName(1), nymd, nhms,                    &
                        im, jm, 1, km, CPP_WZ, rc )
         if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put ze')
         endif
#if defined( SPMD )
         call mp_gather4d(hs, hstmp, im, jm, 1, 1, jfirst, jlast,            &
                          1,  1,     0,  0,  0)
#endif
         if (gid == 0) then
         call GFIO_PutVar (out_fid, oName(2), nymd, nhms,                    &
                        im, jm, 0, 1, CPP_HS, rc )
         if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put zs')
         end if
#if defined( SPMD )
         call mp_gather4d(slp, slptmp, im, jm, 1, 1, jfirst, jlast,          &
                          1,   1,      0,  0,  0)
#endif
         if (gid == 0) then
         call GFIO_PutVar (out_fid, oName(3), nymd, nhms,                    &
                        im, jm, 0, 1, CPP_SLP, rc )
         if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put slp')
         end if
#if defined( SPMD )
         call mp_gather4d(ps, pstmp, im, jm, 1, 1, jfirst, jlast,            &
                          1,  1,     0,  0,  0)
#endif
         if (gid == 0) then
         call GFIO_PutVar (out_fid, oName(4), nymd, nhms,                    &
                        im, jm, 0, 1, CPP_PS, rc )
         if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put ps')
         end if
#if defined( SPMD )
         call mp_gather4d(ua, uatmp, im, jm, km, 1, jfirst, jlast,           &
                          1,  km,    0,  0,  0)
#endif
         if (gid == 0) then
         call GFIO_PutVar (out_fid, oName(5), nymd, nhms,                    &
                        im, jm, 1, km, CPP_UA, rc )
         if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put ua')
         end if
#if defined( SPMD )
         call mp_gather4d(va, vatmp, im, jm, km, 1, jfirst, jlast,           &
                          1,  km,    0,  0,  0)
#endif
         if (gid == 0) then
         call GFIO_PutVar (out_fid, oName(6), nymd, nhms,                    &
                        im, jm, 1, km, CPP_VA, rc )
         if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put va')
         end if
#if defined( SPMD )
         call mp_gather4d(ta, tatmp, im, jm, km, 1, jfirst, jlast,           &
                          1,  km,    0,  0,  0)
#endif
         if (gid == 0) then
         call GFIO_PutVar (out_fid, oName(7), nymd, nhms,                    &
                        im, jm, 1, km, CPP_TA, rc )
         if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put ta')
         end if
#endif
      endif


      if ( index(oflnm, 'prog.bin') .gt. 0 ) then
        if(zstat) call wrt3d(iout,nrec,im,jm,km,epv, jfirst, jlast)
#if defined (GFIO)
      else
#if defined( SPMD )
         call mp_gather4d(epv, epvtmp, im, jm, km, 1, jfirst, jlast,         &
                          1,   km,     0,  0,  0)
#endif
         if (gid == 0) then
         if(zstat) call GFIO_PutVar (out_fid, oName(8), nymd, nhms,           &
                         im, jm, 1, km, CPP_EPV, rc )                         
         if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put pv')
         end if
#endif
      end if

 
! Do not remove the following statement; it is for fvPSAS internal usage 
      if(gid == 0) then
         write(*,*) ' '
         write(*,'(a,I10.8,i10.6)') 'refout_time:', nymd, nhms
      endif
 
! Check global maximum/minimum

      call pmaxmin('PS', ps, qmin, qmax, im, (jlast-jfirst+1), 0.01)
      call pmaxmin('Ztop', wz(1,jfirst,1), qmin, qmax, im, (jlast-jfirst+1), rgf)
      ztop = gmean_of(im, jm, jfirst, jlast, wz(1,jfirst,1)) * rgf
      if( gid == 0) write(6,*) 'Mean model top height (km) =', ztop 

      if(nq.ne.0) then
        do ic=1,nq

!$omp parallel do private(i,j,k)
          do k=1,km
            do j=jfirst,jlast
              do i=1,im
                epv(i,j,k) = q(i,j,k,ic)
              enddo
            enddo
          enddo

          call pmaxmin('Q ', epv,qmin, qmax, im*(jlast-jfirst+1), km, 1.)
          if ( index(oflnm, 'prog.bin') .gt. 0 ) then
            call wrt3d(iout,nrec,im,jm,km,epv, jfirst, jlast)
#if defined (GFIO)
          else
#if defined( SPMD )
         call mp_gather4d(epv, epvtmp, im, jm, km, 1, jfirst, jlast,            &
                          1,   km,     0,  0,  0)
#endif
         if (gid == 0) then
             call GFIO_PutVar (out_fid, oName(9), nymd, nhms,                   &
                         im, jm, 1, km, CPP_EPV, rc )
             if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put Q')
         end if
#endif
          end if
        enddo

      endif
 
      call pmaxmin('U ', ua, qmin, qmax, im*(jlast-jfirst+1),km, 1.)
      call pmaxmin('V ', va, qmin, qmax, im*(jlast-jfirst+1),km, 1.)
      call pmaxmin('T ', ta, qmin, qmax, im*(jlast-jfirst+1),km, 1.)

      if ( ccm_on ) then 
! Writing ground temperature from lsm

        call pmaxmin('TG', tg(1,jfirst), qmin, qmax, im,(jlast-jfirst+1),1.)
        if ( index(oflnm, 'prog.bin') .gt. 0 ) then
          call wrt3d(iout,nrec,im,jm,1,tg(1,jfirst), jfirst, jlast)
#if defined (GFIO)
        else
#if defined( SPMD )
         call mp_gather4d(tg, tgtmp, im, jm, 1, 1, jfirst, jlast,              &
                          1,  1,     0,  0,  0)
#endif
         if (gid == 0) then
           call GFIO_PutVar (out_fid, oName(10), nymd, nhms,                   &
                         im, jm, 0, 1, CPP_TG, rc )
           if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put TG')
         end if
#endif
        end if

! SJL: Total precip (mm/dy)
        call pmaxmin('PRECP', precp(1,jfirst), qmin, qmax, im,(jlast-jfirst+1),1.)
        if ( index(oflnm, 'prog.bin') .gt. 0 ) then
          call wrt3d(iout,nrec,im,jm,1,precp(1,jfirst), jfirst, jlast)
#if defined (GFIO)
        else
#if defined( SPMD )
         call mp_gather4d(precp, precptmp, im, jm, 1, 1, jfirst, jlast,        &
                          1,     1,        0,  0,  0)
#endif
         if (gid == 0) then
          call GFIO_PutVar (out_fid, oName(11), nymd, nhms,                    &
                         im, jm, 0, 1, CPP_PRECP, rc )
          if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put precp')
         end if
#endif
        end if

      endif		! ccm_on

!$omp parallel do private(i,j,k)
       do k=1,km
          do j=jfirst,jlast
            do i=1,im
              epv(i,j,k) = omga(i,k,j)
            enddo
          enddo
       enddo
       if ( index(oflnm, 'prog.bin') .gt. 0 ) then
         call wrt3d(iout,nrec,im,jm,km,epv, jfirst, jlast)
#if defined (GFIO)
       else
#if defined( SPMD )
         call mp_gather4d(epv, epvtmp, im, jm, km, 1, jfirst, jlast,           &
                          1,   km,     0,  0,  0)
#endif
         if (gid == 0) then
          call GFIO_PutVar (out_fid, oName(12), nymd, nhms,                    &
                         im, jm, 1, km, CPP_EPV, rc )
           if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Put TG')
         end if
#endif
       end if


! dpdt unit is pascal per second; convert to mb/day
      fac = 24 * 3600  / 100
      call pmaxmin('W(mb/day)', epv,qmin,qmax,im*(jlast-jfirst+1),km, fac)

#if defined (GFIO)
      if (gid == 0) then
      if ( index(oflnm, 'prog.bin') .le. 0 ) then
         call GFIO_Close ( out_fid, rc )
      end if
      end if
#endif


      return
!EOC
      end

