      subroutine setup(iuhsx,  iuicx,  ioutx,  iphysx, izoux ,  ih2ox, 
     &                 zstatx, ccm3x,  rayfx,  co2fx,  mounx,  ndoutx, ngout2dx,
     &                 ngoutx, nymdbx, nhmsbx, nymdex, nhmsex,  pdtx,
     &                 nsplitx, mpdt,  umaxx,  iordx,  jordx,  kordx,
     &                 ictpx,  jobx,   mstep,  isw, ilw, iae,
     &                 lcroot, bndtvs, sstcycx, bndtvo, ozncycx,
     &                 h2odatax, srfdatax, diagx, dcaf, gwdf,
     &                 psasdt, nymd_p, nhms_p,
     &                 fcstbegx, fcstfrqx, fcsthrsx, nhms_dfx, cutoff_dfx,
     &                 diag2dnmx, diag3dnmx, prog3dnmx, prec_dasx,
     &                 trajfrqx, trajtmplx, ptrjfrqx, ptrjtmplx, rmaptrjx)

#if defined( SPMD )
#define CPP_PRT_PREFIX  if(gid.eq.0)
      use mod_comm,  only : gid
#else
#define CPP_PRT_PREFIX
#endif

      implicit none
#include <comadj.h>

      character*80 job, jobx
      character*80  lcroot,  bndtvs,  bndtvo
      character*80  ncpath, sstdata, ozndata
      character*80  h2odata, h2odatax
      character*80  srfdata, srfdatax
      integer  pdt, pdtx
      integer  nsplit, nsplitx
      integer  mdt, mpdt
      logical    moun,  rayf,  co2f
      logical    mounx, rayfx, co2fx
      logical    zstat, zstatx
      logical    ccm3, ccm3x
      logical    diag, diagx
      logical    sstcyc, sstcycx
      logical    ozncyc, ozncycx
      logical    dcaf, gwdf

      integer    iuhs, iuic, iout
      integer    iphys,  iphysx
      integer    iuhsx, iuicx, ioutx, ictpx
      integer    izou,  izoux
      integer    ih2o,  ih2ox

      integer       nymdb, nhmsb
      integer       nymdbx, nhmsbx
      integer       nymde, nhmse
      integer       nymdex, nhmsex
      integer       ndout, ndoutx
      integer       ngout, ngoutx
      integer       ngout2d, ngout2dx
      integer       ngout3d
      character*256 diag2dnmx, diag3dnmx, prog3dnmx
      character*256 trajtmplx, ptrjtmplx
      integer       nday
      integer       mstep
      integer       iord, jord, kord
      integer       iordx, jordx, kordx
      integer       fcstbeg, fcstbegx
      integer       fcstfrq, fcstfrqx
      integer       trajfrq, trajfrqx
      integer       ptrjfrq, ptrjfrqx
      integer       rmaptrj, rmaptrjx
      integer       fcsthrs(24), fcsthrsx(24)
      integer       nhms_dfx
      integer       cutoff_dfx
      integer       prec_dasx
      real          umax, umaxx
      integer iradsw, isw
      integer iradlw, ilw
      integer iradae, iae
      integer lenchr
      integer nsrest
      integer n, nn
      integer nhmse24

! For PSAS
! time step for PSAS_RUC
      integer psasdt, fvpsasdt
! Ending Date/time for PSAS
      integer nymd_p, nymd_psas
      integer nhms_p, nhms_psas
! Precision for DAS output    
! prec_das = 0 : 32 bits
! prec_das = 1 : 64 bits
      integer prec_das

! Digital filter
      integer nhms_df
      integer cutoff_df

      character*256 trajtmpl, ptrjtmpl
      character*256 diag2dnm, diag3dnm, prog3dnm ! output file names for diag2d,
                                                 ! diag3d and prog.


      namelist /input/  job     , nymdb   , nhmsb   , nymde   ,
     .                  nhmse   , ndout   , ngout   , ngout2d , 
     .                  pdt     ,
     .                  fvpsasdt, nymd_psas, nhms_psas,
     .                  nsplit  , mdt     , umax    , iord    ,
     .                  jord    , kord    , nsrest  , nday    ,
     .                  iuhs    , iuic    , iout    , iphys   ,
     .                  izou    , ih2o    , zstat   , ccm3    ,
     .                  rayf    , co2f    , moun    , iradsw  ,
     .                  iradlw  , iradae  , ncpath  , sstdata ,
     .                  sstcyc  , ozndata , ozncyc  , h2odata ,
     .                  srfdata , diag    , fcstbeg , fcstfrq ,
     .                  fcsthrs, nhms_df  , cutoff_df, diag2dnm,
     .                  diag3dnm, prog3dnm, prec_das,
     .                  trajfrq , trajtmpl, ptrjfrq, ptrjtmpl,
     .                  rmaptrj


      namelist /flags/ dcaf, nlvdry, gwdf

#if ( !defined SPMD )
      integer gid
      gid = 0
#endif

! Defaults:
      job = 'fv_gcm00'
      lcroot = './'
      bndtvs = 'sst.nc'
      sstcyc = .true.
      bndtvo = 'ozone.nc'
      ozncyc = .true.
      h2odata = 'h2o.data'
      srfdata = 'surf.data'

      iuhs  = 80         ! UNIT to read topography data set.
      iuic  = 81         ! UNIT to restart file
      iout  = 82         ! UNIT to write A-grid output.
      iphys = 83         ! UNIT to physics restart file
      izou  = 99         ! UNIT to write zonal mean stat
      ih2o  = 13         ! UNIT to physics restart file

      zstat = .true.
      ccm3  = .false. 
      rayf  = .false.
      co2f  = .false.
      moun  = .true.

! OUTPUT INTERVAL:
      ndout = 240000       ! refout
      ngout = 240000       ! diag
      ngout2d = -1         ! diag
      prog3dnm = 'prog.bin'
      diag2dnm = 'diag.sfc.bin'
      diag3dnm = 'diag.bin'
      trajtmpl = 'traj.bin'  ! state vector trajectory default filename
      ptrjtmpl = 'ptrj.bin'  ! partial physics trjectory default filename



! Physics time step is PDT
! Remapping interval is MDT

      pdt = 1800
      nsplit = 0
      mdt = 3600
      umax = 300.

      iradsw = -1
      iradlw = -1
      iradae = -12

      iord = 4
      jord = 4
      kord = 4

! restart type:

      nsrest = 1          ! Read from a restart file
      nday = 365
! NCAR phys.
      nlvdry = 6          ! comadj.h

! fvPSAS realted defaults:
      fvpsasdt = 0       ! Default for Safety 
      trajfrq  = -99     ! Default: do not write out trajectory in hhmmss
      ptrjfrq  = -99     ! Default: do not write out physics trajectory in hhmmss
      rmaptrj  = -99     ! Default: freq of remap of input trajectory (hhmmss)
      nymd_psas = 0
      nhms_psas = 0
      fcstbeg = 0
      fcstfrq = 0
      fcsthrs(1) = 0
      do n = 2, 24
        fcsthrs(n) = -1
      enddo
      nhms_df = 0   ! no digital filter by default
      cutoff_df = 8 ! 8 hours cutoff period digital filter by default
! Precision for DAS output    
! prec_das = 0 : 32 bits
! prec_das = 1 : 64 bits
      prec_das  = 1    ! 64 bits by default

! Multiple read/close in MPI mode? (works on the SGI) 

      read (811,input) ! RT: unit number to avoid conflict w/ other codes
      close(811)

      read(813,flags)  ! RT: unit number to avoid conflict w/ other codes
      close(813)

      lcroot = ncpath
      sstcycx = sstcyc
      ozncycx = ozncyc

      if(gid==0) then
         bndtvo = lcroot(:lenchr(lcroot))//'/'//ozndata
         bndtvs = lcroot(:lenchr(lcroot))//'/'//sstdata
         h2odatax = lcroot(:lenchr(lcroot))//'/'//h2odata
         srfdatax = lcroot(:lenchr(lcroot))//'/'//srfdata
      endif

! Copy fields
      jobx = job
      nymdbx = nymdb
      nhmsbx = nhmsb

      nymdex = nymde
      nhmsex = nhmse

      iuhsx  = iuhs
      iuicx  = iuic
      ioutx  = iout
      iphysx = iphys
      izoux  = izou
      ih2ox  = ih2o

      zstatx = zstat
      ccm3x  = ccm3
      rayfx  = rayf
       co2fx = co2f       ! flag to turn on BAU co2 
      mounx  = moun

      diagx  = diag

! OUTPUT INTERVAL:

      ndoutx  = ndout
      pdtx    = pdt

      ngoutx  = ngout
      ngout2dx = ngout2d

      prog3dnmx = prog3dnm
      diag2dnmx = diag2dnm
      diag3dnmx = diag3dnm

      if (ngout2dx .eq. -1 .and. trim(diag2dnmx) .ne. 'diag.sfc.bin') then
        ngout2dx = ngoutx
      endif

      rmaptrjx  = rmaptrj
      trajfrqx  = trajfrq
      ptrjfrqx  = ptrjfrq
      trajtmplx = trajtmpl
      ptrjtmplx = ptrjtmpl

! PSAS related:
      psasdt = fvpsasdt

      if (nymd_psas .eq. 0) then
          nymd_p = nymde
      else
          nymd_p = nymd_psas 
      endif

      if (nhms_psas .eq. 0) then
          nhms_p = nhmse
      else
          nhms_p = nhms_psas 
      endif
      fcstbegx = fcstbeg
      fcstfrqx = fcstfrq
      nhms_dfx = nhms_df
      cutoff_dfx = cutoff_df
      prec_dasx  = prec_das

!  check for duplicate hours
      do n = 1, 24
        if (fcsthrs(n) .gt. 23) then
          fcsthrs(n) = -1
        endif
      enddo
      do n = 1, 23
        do nn = n+1, 24
          if (fcsthrs(nn) .eq. fcsthrs(n)) then
            fcsthrs(nn) = -1
          endif
        enddo
      enddo
      do n = 1, 24
        fcsthrsx(n) = fcsthrs(n)
      enddo
      if (fcstfrqx .eq. 0) then
        fcstbegx = 0
      endif

      nsplitx = nsplit
      mpdt    = mdt
      umaxx   = umax

      iordx = iord
      jordx = jord
      kordx = kord

      ictpx = nsrest

      mstep =  nday * 24 * 3600 / pdt
      if (nhmse.ne.nhmsb) then  ! RG: run for other than whole-day periods
          nhmse24 = nhmse
          if (nhmse24.lt.nhmsb) nhmse24 = nhmse24 + 240000
          mstep = mstep + nint((nhmse24 - nhmsb)/10000.) * 3600 / pdt
      end if  !  RG
        
      CPP_PRT_PREFIX write(6,*) 'max. # of steps for this run =', mstep
      CPP_PRT_PREFIX write(6,*) 'size of time step in seconds = ', pdt

      if ( .not. ccm3 ) return
! Convert iradsw and iradlw from hours to timesteps if necessary
 
      if (iradsw.lt.0) iradsw = nint((-iradsw*3600.)/float(pdt))
      if (iradlw.lt.0) iradlw = nint((-iradlw*3600.)/float(pdt))
 
! Convert iradae from hours to timesteps if necessary
 
      if (iradae.lt.0) iradae = nint((-iradae*3600.)/float(pdt))
 
! iradae must be an even multiple of iradlw
 
      if (mod(iradae,iradlw).ne.0) then
        write(6,*)'iradae must be an even multiple of iradlw.'
        write(6,*)'iradae = ',iradae,', iradlw = ',iradlw
        stop
      end if

      CPP_PRT_PREFIX write(6,*) 'job name:', job
!     CPP_PRT_PREFIX write(6,*) 'iradsw = ', iradsw
!     CPP_PRT_PREFIX write(6,*) 'iradlw = ', iradlw
!     CPP_PRT_PREFIX write(6,*) 'iradae = ', iradae

      isw = iradsw
      ilw = iradlw
      iae = iradae

      return
      end
