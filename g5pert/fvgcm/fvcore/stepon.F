      module stepon

!
! Driver for the finite-volume GCM and DAS using NCAR CCM 3 physics
!
! Original developer:                  S.-J. Lin 1998
! Initial software engineering:        G. Lai 1998/1999
! Message-passing implementation:      W. Sawyer 1999
! Hybrid MPI-OpenMP performance tuning: SJL Jan 2000
! fvPSAS interface:                    da Silva, Feb 2000

! !REVISION HISTORY:
!
!  22Feb2000  da Silva  Now calling fvPSAS at bottom of time loop.
!  18May2000  da Silva  Using file resolver for SST file names (fvPSAS only)
!  18Jun2000  sawyer    Changed to precision module
!  19Jun2000  mirin     Initial changes for 2D decomposition
!  13Jul2000  sawyer    Incorporated PILGRIM simplifications
!  23Jan2001  g. lou    Digital filter implementation
!  20Apr2001  j. Chern  Moist physics adjustment after fvPSAS
!  22Nov2006  Oloso     Ifdef'd timer for SPMD
!  14May2007  Todling   Introduced dyn_prog; global change.
!  25Oct2007  Todling   Add opt to skip rst's input
!  21Nov2007  Todling   Resolution-free executable
!

      use precision
      use prognostics
      use m_pmgrid, only: pmgrid_set
      use m_comsrf, only: comsrf_init, oro, ts

#if defined( SPMD )
        use mod_comm, only : gid, y_decomp
#if defined (PIO)
        use mod_pio, only  : pio_init, pio_exit,
     &                       allocate_shm_4gout, deallocate_shm_4gout
#endif
#endif

#if defined ( fvPSAS_ON )
      use m_dfilter        ! digital filter parameters
      use m_FileResolv     ! File resolver 
      use m_fvpsas
      use  mapz_module, only: te_map
#endif
#if defined (GFIO)
      use m_StrTemplate    ! grads style template
#endif                  

      use timingModule, only: timing_init,timing_on,timing_off

      implicit      none

      real(r8)     undef
      parameter   (undef = 1.e25)

      logical masterproc
! =================
! Dynamical fields:
! =================

! Geometric arrays
      real(r8),allocatable::   sine(:)
      real(r8),allocatable::   cosp(:)
      real(r8),allocatable::   sinp(:)
      real(r8),allocatable::   cose(:)
      real(r8),allocatable::   sinlon(:)
      real(r8),allocatable::   coslon(:)

! Eulerian vertical coordinate
      real(r8),allocatable:: ak(:)
      real(r8),allocatable:: bk(:)
      real(r8) ptop                     ! pressure at model top
      real(r8) ptopk                    ! ptop**cappa
      real(r8) pint                     ! pressure at the interface
      integer ks                    ! number of layers in pure p region

      real(r8), allocatable :: omga(:,:,:)  ! dp/dt; vertical pressure velocity
      real(r8), allocatable :: pref(:)       !reference pressure @ eta layer center (mb)
      real(r8), allocatable :: pref_rev(:)   ! pref bottom-up
      real(r8)       psref
      parameter (psref = 100000.)   !reference surface pressure (Pascal)

! GSFC level 2.5 turbulence clousure
      real(r8), allocatable :: tke(:,:,:) !turbulent kinetic energy
      logical qbeg

! Include files from CCM
#     include <comcon.h>
#     include <comctl.h>
#     include <comtim.h>
#     include <commss.h>

! Heating rates:
      real(r8), allocatable:: qrs(:,:,:)    ! Shortwave heating rate
      real(r8), allocatable:: qrl(:,:,:)    ! longwave heating rate
      real(r4), allocatable:: diagbuf(:,:,:)! diag buffer

#     include <diag.h>

      integer, save :: rstpick_ = .false.
      integer, save :: rstskip_ = .false.
      integer      cdiag(pdiag)      ! counter
      integer      inx(pdiag)        ! sorted field table indices
      integer      qvdim(pdiag)      ! vertical dimension of field
      integer      fldloc(pdiag)     ! location of field in diagbuf
      integer      pick(pdiag)       ! flag for selected fields (1=true)
      integer      tavg(pdiag)       ! flag for time-averaged flds (1=true)
      integer      diagattr(5,pdiag) ! field attributes
      integer      nrec2d            ! diagnostics output record number for 2d
      integer      nrec              ! diagnostics output record number
      integer      nslice            ! tot. records of selected diag. flds
      character*8  qname(pdiag)          ! prescribed field table
      character*8  fldtbl(pdiag)         ! sorted field table
      character*16 qunit(pdiag)          ! unit for diagnosis fields
      character*80 qdesc(pdiag)          ! description for diagnosis fields
      real(r8), allocatable ::  pblht(:, :)
      real(r8), allocatable ::  tpert(:, :)
      real(r8), allocatable ::  qpert(:,:,:)
#if defined( CCM )
#     include <physbuf.h>
#endif
      real(r8),allocatable :: precp(:,:)

      character*80 job
      character*80 h2odata
      character*80 srfdata
      character*256 diag2dnm, diag3dnm, prog3dnm ! file names from ccmrun.namelist
      character*256 progflnm                     ! prog file name
      character*256 trajtmpl, ptrjtmpl  ! templates for traj/phy-traj filenames


! Zonal mean stat.
      real(r8), allocatable :: psx(:), pvx(:,:)
      real(r8), allocatable :: ux(:,:), vx(:,:), tx(:,:), qx(:,:)
      real(r8), allocatable :: wx(:,:)
      real(r8), allocatable :: usus(:,:), usvs(:,:)
      real(r8), allocatable :: vsts(:,:), vspv(:,:), tsts(:,:)
!
!JDC save pt for moist physics adjustment after fvPSAS
!
#if defined ( fvPSAS_ON )
#if defined ( DIAPER )
      real(r8), allocatable :: ptold(:,:,:)   ! Virture potential temperature before fvPSAS
#endif
#endif

! Time steps for various components of the model
      integer pdt         !
      integer mpdt        ! mapping time step
      integer nsplit      ! sub_steps for the dynamics
      real(r8) tdt        ! tme step for tphysac and surface models
      real(r8) bdt        ! time step for tphysbc

      integer fvpsasdt    ! analysis time step (in NHMS format)
      integer nymd_psas   ! YYMMDD to end psas
      integer nhms_psas   ! HHMMSS to end psas
      integer dasdt       ! analysis time step (in seconds)
      integer fcstbeg     ! YYYYMMDD of first restart write for forecast (0=today)
      integer fcstfrq     ! frequency in days for restart write for forecast (1=every day)
      integer fcsthrs(24) ! hours of restart write for forecast (0-23 GMT)
      integer fcstnymd    ! next YYYYMMDD to write out forecasts
      integer trajfrq     ! frequency in HHMMSS of trajectory output
      integer ptrjfrq     ! frequency in HHMMSS of physics trajectory output
      integer rmaptrj     ! frequency in HHMMSS of remap of external trajectory
      integer hour        ! current model hour
      logical dofvpsas        ! determine if it is analysis time 
      logical forecast        ! Forecast mode under fvPSAS

      logical moun
      logical rayf
      logical co2f            ! flag to turn on "BAU" co2
      logical dcaf            ! flag to to do dry convection
      logical gwdf            ! flag to do GWD

      logical zstat           ! Compute zonal mean stats
      logical ncar
      logical diag
      logical touch           ! first touch

      logical aoa             ! flag to do age of air computation
      logical force_dycore    ! flag to do dycore forcing
      logical convt           ! flag to control pt output from fvcore()
                              ! Set convt =.T. to output (virtual) temperature
      logical consv           ! flag to force conservation of Total
                              ! Energy
      logical ccm_on          ! running with CCM physics? 

      real(r8)  model_time    ! model time in seconds since init

! I/O units:
      integer  iuhs       ! unit to terrian data
      integer  iout       ! unit to A-grid dynamics output
      integer  jout2d     ! unit to surface diagnostics  
      integer  jout       ! unit to diag
      integer  iuic       ! Unit to dynamics IC
      integer  iphy       ! unit to physics IC
      integer  ilsm       ! unit to LSM restarts (forecast)
      integer  izou       ! unit to zonal mean stat
      integer  ih2o       ! unit to h2o data
      integer  ictl2d     ! unit to GrADS control file
      integer  ictl       ! unit to GrADS control file
      integer  aout       ! unit to Age_of_Air output

      integer  iavg
      integer  nymd
      integer  nhms

      integer  nymdb
      integer  nhmsb

      integer  nymde
      integer  nhmse

      integer  nymdd
      integer  nhmsd

      integer  gdate2d
      integer  gtime2d
      integer  gdate
      integer  gtime          
      integer  ndout
      integer  ngout2d  
      integer  ngout
      integer  idtout
      integer  igtout2d
      integer  igtout
      integer  n2ddiag              ! number of user activated 2d diags
      integer  n3ddiag              ! number of user activated 3d diags
      integer  diagstart      ! position of diag to begin output
      integer  diagend        ! position of diag to end output
      logical  diagfsfc             ! true if surface diag in separate file (diag.sfc.bin)
      logical  diagf                ! true if diag file needed (diag.bin)   
      integer  istep          ! time loop index
      integer  mstep          ! integration steps for this run
      integer  i, j, k, ic
      integer  nt
      integer  ndt            ! temporary value of time in seconds
      integer  idummy 
      integer  snhms, strajfrq, sptrjfrq

      integer(i8)  flsz

      real(r8)      umax
      real(r8)      dp, dl

! Remapping after analysis related variables
#if defined ( fvPSAS_ON )
      real(r8)    delpmin     ! minimum delp at the lowest level after analysis
      real(r8)    te0
      integer     nx
      integer     nremap      ! Total number of points need to remap
#endif

!  Digital filter related variables

#if defined ( fvPSAS_ON )
! nhms_df and cutoff_df are fed in through namelist in setup routine. 
! All other digital filter related parameters are derived.
      integer       ncount_df       ! filter window steps
      integer       mstep_hf        ! half of filter window steps
      integer       nstep_hf        ! half way through filter window steps
      integer       nhms_hf         ! half way through filter window hours
      integer       nymd_hf         ! half way through filter window ymd
#endif
      integer       nhms_df         ! filter window length in hhmmss 
      integer       cutoff_df       ! digital filter cutoff period in hh
      integer       prec_das        ! precision for DAS output
                                       ! 0: 32 bits
                                       ! 1: 64 bits (default)

! End of digital filter variables

      character*80 drstflnm, prstflnm, lrstflnm
      character*80 diagflnm2d, ctlflnm2d, diagflnm, ctlflnm     

      integer   ng_d, ng_s

! NCAR CCM block data

#if defined( CCM )
      external lakconi
      external phyconi
      external radconi
      external snoconi
      external soiconi
      external vegconi
      external vegtypi
#endif
      integer ierror, ier

      integer stepon_tape_rec                   ! FastOpt
      integer nouter                            ! FastOpt
      integer ninner                            ! FastOpt
      integer ninn                              ! FastOpt

!_RT  integer im                                ! FastOpt allow allocation use im
!_RT  parameter( im = imr )                     ! FastOpt allow allocation use im

      real(r8), pointer ::    u(:,:,:)   ! zonal wind on D-grid
      real(r8), pointer ::    v(:,:,:)   ! meridional wind
      real(r8), pointer ::   pt(:,:,:)   ! virtual potential temperature
      real(r8), pointer :: delp(:,:,:)   ! pressure thickness (pascal)
      real(r8), pointer ::    q(:,:,:,:) ! specific humidity & tracer mixing ratios

!     real(r8), pointer ::   pe(:,:,:)   ! pressure at layer edges
!     real(r8), pointer ::   pk(:,:,:)   ! pe**kappa
!     real(r8), pointer ::  pkz(:,:,:)   ! layer-mean value of pk
!     real(r8), pointer ::   ps(:,:)     ! surface pressure (pascal)
!     real(r8), pointer :: peln(:,:,:)   ! pressure edges

      interface stepon_set; module procedure stepon_set_; end interface stepon_set

      contains


      integer function nsecf( n )
      implicit none
      integer n
      nsecf = n/10000*3600 + mod(n,10000)/100* 60 + mod(n,100)
      end function nsecf


      integer function nhmsf( n )
      implicit none
      integer n
      nhmsf = n/3600*10000 + mod(n,3600 )/ 60*100 + mod(n, 60)
      end function nhmsf


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_initialize --- Initialization for stepon
!
! !INTERFACE:
      subroutine stepon_initialize ( prog )

      use fvcore, only: fvcore_initialize
      use gmean, only: gmean_initialize

      implicit none

      type(dyn_prog), TARGET :: prog

      progflnm = 'prog.bin'

      kfirst = 1
      klast = nl

#if defined( SPMD )
#if defined (PIO)
      call pio_init()
#endif
      call y_decomp( jnp, nl, jfirst, jlast, kfirst, klast, gid)
      if ( gid == 0 ) then
          masterproc = .true.
      else
          masterproc = .false.
      endif
#if (defined PIO) && (!defined NO_PIO_BUFFER )
      call allocate_shm_4gout(jfirst, jlast)
#endif
#else
      jfirst = 1
      jlast  = jnp
      masterproc = .true.
#endif
      call pmgrid_set ( imr, jnp, nl, jfirst, jlast )
      call comsrf_init()

      call timing_on('TOTAL')

       klastp = klast
      if (klast .eq. nl) klastp = klast + 1

      nq = nc
      aout = 98                 ! Output unit for AOA computation
      aoa = .false. 
      touch = .true.
      ilsm = 83			! Reuse phys restart unit number

! Get control parameters from namelist
! Note: need to broadcast within setup for portability

      call setup(iuhs, iuic, iout, iphy, izou, ih2o,
     &           zstat, ncar, rayf, co2f, moun, ndout, ngout2d,
     &           ngout, nymdb, nhmsb, nymde, nhmse, pdt,
     &           nsplit, mpdt, umax, iord, jord, kord,
     &           nsrest, job,  mstep, iradsw, iradlw, iradae,
     &           lcroot, bndtvs, sstcyc, bndtvo, ozncyc,
     &           h2odata, srfdata, diag, dcaf, gwdf,
     &           fvpsasdt, nymd_psas, nhms_psas,
     &           fcstbeg, fcstfrq, fcsthrs, nhms_df, cutoff_df,
     &           diag2dnm, diag3dnm, prog3dnm, prec_das,
     &           trajfrq, trajtmpl, ptrjfrq, ptrjtmpl, rmaptrj )

#if defined( SPMD )
        ng_d = max(2,min(abs(jord),3)) 
        ng_s = 3
#else
        ng_d = 0
        ng_s = 0
#endif

!
! Allocate variables
!
      call prognostics_allocate ( prog )

! Set pointers
! ------------
       u   => prog%u
       v   => prog%v
      pt   => prog%pt
      delp => prog%delp
      q    => prog%q

!     pe   => prog%pe
!     pk   => prog%pk
!     pkz  => prog%pkz
!     ps   => prog%ps
!     peln => prog%peln

      allocate( omga(imr,kfirst:klast,jfirst:jlast) )

#if defined ( CCM )
      allocate( qrs(imr,nl,jfirst:jlast) )
      allocate( qrl(imr,nl,jfirst:jlast) )
#endif
                                                           !mixing ratios
#if defined (ALT_PBL)
      allocate(  tke(imr,jfirst:jlast,nl) )      ! turbulent kinetic energy
! FastOpt begin make sure tke is allocated for refout
#else
      allocate(  tke(1,1,1) )                    ! turbulent kinetic energy
! FastOpt end
#endif

! Zonal mean stat.
      allocate( psx(jfirst:jlast) )
      allocate( ux(jfirst:jlast,kfirst:klast) )
      allocate( vx(jfirst:jlast,kfirst:klast) )
      allocate( tx(jfirst:jlast,kfirst:klast) )
      allocate( qx(jfirst:jlast,kfirst:klast) )
      allocate( wx(jfirst:jlast,kfirst:klast) )
      allocate( pvx(jfirst:jlast,kfirst:klast) )
      allocate( usus(jfirst:jlast,kfirst:klast) )
      allocate( usvs(jfirst:jlast,kfirst:klast) )
      allocate( vsts(jfirst:jlast,kfirst:klast) )
      allocate( vspv(jfirst:jlast,kfirst:klast) )
      allocate( tsts(jfirst:jlast,kfirst:klast) )

      allocate( sine(jnp) )
      allocate( cosp(jnp) )
      allocate( sinp(jnp) )
      allocate( cose(jnp) )
      allocate( sinlon(imr) )
      allocate( coslon(imr) )
      allocate( ak(nl+1) )
      allocate( bk(nl+1) )

      allocate ( pref(nl) )
      allocate ( pref_rev(nl) )

      allocate ( precp(imr,jnp) )
      allocate ( pblht(imr, jnp) )
      allocate ( tpert(imr, jnp) )
      allocate ( qpert(imr,nc,jnp) )

      idtout = nsecf(ndout)
      if (ngout2d .ne. -1) then
        igtout2d = nsecf(ngout2d)
      else
        igtout2d = -1
      endif                   
      igtout = nsecf(ngout)

      dtime = pdt
      bdt   = pdt
      tdt   = pdt

      if(masterproc) print *
      if(masterproc) print *, 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
 
!
! WS: By default assume CCM physics is *off*
!
      ccm_on = .false.
      force_dycore = .false. 

#if defined( CCM )
      ccm_on = .true.
      call blk_dat              !This replaced CCM's block data
      convt = .true.

#if defined ( fvPSAS_ON )
      if(masterproc) print *, 'NASA/DAO Finite-Volume/Physical-Space DAS'
      consv = .false.
      aoa = .false. 
#else
      if(masterproc) print *, '     NASa/nCAR Finite-Volume GCM '
      consv = .true.
#endif

#if defined (ALT_GWD)
      if(masterproc) print *, 'Alternative GWD scheme will be used'
#endif
#if defined (ALT_CPS)
      if(masterproc) print *, 'Alternative cumulus scheme will be used'
#endif

#else
      if(masterproc) print *, 'Finite-Volume Dynamical Core (MLP+OpenMP)'
      convt = .false.
      force_dycore = .false.
      consv = .true. 
#endif

! set consv = F so as to pass the restart test with -O3
      consv = .false.

      if(masterproc) print *, '          Version ' // 'fvgcm'

      if(masterproc) print *,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      if(masterproc) print *

! GEOMETRIC FACTORS
      call setrig(imr, jnp, dp, dl, cosp, cose, sinp, sine)

! Set up reference vertical coordinate
      call set_eta(nl, ks, ptop, pint, ak, bk)

      do k = 1, nl
        pref(nl-k+1) = .01 * .5 * (ak(k) + ak(k+1) +
     .                  (bk(k) + bk(k+1)) * psref)
      end do
      do k = 1, nl
         pref_rev(k) = pref(nl-k+1)
      end do


! Define logitude at the center of the finite-volume
! i=1, Zamda = -pi

      do i=1,imr/2
         coslon(i)       =  cos((i-1)*dl)
         coslon(i+imr/2) = -coslon(i)
         sinlon(i)       =  sin((i-1)*dl)
         sinlon(i+imr/2) = -sinlon(i)
      enddo

! Initialize  dynamical fields in "parallel" (SMP)
      call initmem(imr, jfirst, jlast, kfirst, klast, klastp,
     &             u, v, pt, q, nq, pk, pe, delp, undef, ng_s, ng_d)
      call initp3d(imr, jfirst, jlast, kfirst, klast, omga, 0.0)

      do j=jfirst,jlast
         do i=1, imr
            precp(i,j) = undef
         enddo
      enddo

#if defined( CCM )
      call initp3d(imr, jfirst, jlast, kfirst, klast, qrs, undef)
      call initp3d(imr, jfirst, jlast, kfirst, klast, qrl, undef)

      if( touch ) then
          if(masterproc)
     &       print *, 'Phys arrays will be initialized in parallel'
      endif
#endif

! Call initialization routine for CCM physics
      call initcon(ncar)


! FastOpt
      call fvcore_initialize( imr, jnp, nl, jfirst, jlast, ng_d, ng_s,
     &                        nsplit, jord, pdt )

      call gmean_initialize( imr, jnp )

! initialise time stepping splitting into out and inner loop
! currently number of timesteps are restricted to multiple of ninner
#ifdef CHECKPOINTING
#  ifdef    EXTERNAL_TRAJECTORY
      ninn   = 1
      ninner = ninn
      nouter = mstep / ninn
#  else
      ninn   = 1
      ninner = ninn
      nouter = mstep / ninn
#  endif /* EXTERNAL_TRAJECTORY */
#else   /*  CHECKPOINTING */
      ninner = mstep
      nouter = 1
#endif  /*  CHECKPOINTING */
!RG test test test
!      nouter = 1
!      ninner = 2
!RG test test test

      if(masterproc) print *,' mstep  = ', mstep
      if(masterproc) print *,' nouter = ', nouter
      if(masterproc) print *,' ninner = ', ninner

!EOC
      end subroutine stepon_initialize
!-----------------------------------------------------------------------


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_set_ --- Initialization for stepon
!
! !INTERFACE:
      subroutine stepon_set_ ( prog, rstpick, rstskip )
      use fvcore, only: fvcore_initialize
      use m_die, only  : die

      implicit none

      character(len=*), parameter :: myname_ = 'stepon_set_'
      type(dyn_prog) :: prog
      logical, optional, intent(in) :: rstpick
      logical, optional, intent(in) :: rstskip

      integer g_s, ng_d, ierr

#if defined( SPMD )
        ng_d = max(2,min(abs(jord),3))
        ng_s = 3
#else
        ng_d = 0
        ng_s = 0
#endif


      allocate(  u(imr,jfirst-ng_d:jlast+ng_s,kfirst:klast), stat=ierr )
         if(ierr/=0) call die(myname_,'Alloc(u)',ierr)
      allocate(  v(imr,jfirst-ng_s:jlast+ng_d,kfirst:klast), stat=ierr )
         if(ierr/=0) call die(myname_,'Alloc(v)',ierr)
      allocate( pt(imr,jfirst-ng_d:jlast+ng_d,kfirst:klast), stat=ierr )  ! virtual potential temperature
         if(ierr/=0) call die(myname_,'Alloc(pt)',ierr)
      allocate( delp(imr,jfirst:jlast,kfirst:klast), stat=ierr )          ! pressure thickness (pascal)
         if(ierr/=0) call die(myname_,'Alloc(delp)',ierr)
      allocate( q(imr,jfirst-ng_d:jlast+ng_d,kfirst:klast,nc), stat=ierr )! specific humidity & tracer
         if(ierr/=0) call die(myname_,'Alloc(q)',ierr)

      if (present(rstpick) ) then
          rstpick_ = rstpick
      endif
      if (present(rstskip) ) then
          rstskip_ = rstskip
      endif
      if (rstskip_) return

! Get initial condition for the dynamics.
      call timing_on('rst_dyn')
      call rst_dyn( nsrest, iuic, trajtmpl, job, rstpick_, imr, jnp, nl, kfirst, klast,
     &              jfirst, jlast, ccm_on, delp, u, v, pt, 
     &              nq, q, ps, ptop, nymd, nhms, nstep, ng_s, ng_d)
      call timing_off('rst_dyn')
 

      ptopk = ptop ** cappa

!$omp parallel do private(i,j,k)
      do j=jfirst,jlast
            do i=1,imr
               pe(i,1,j) = ptop
            enddo

            do k=2,nl+1
            do i=1,imr
               pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            enddo
            enddo
      enddo


      if(  nstep .eq. 0 .and. force_dycore ) then
         nymd = 101
         nhms = 0
         if ( aoa ) then
!$omp parallel do private(i,j,k)
             do k=1,nl
                do j=jfirst,jlast
                   do i=1,imr
                      q(i,j,k,1) = 0.
                   enddo
                enddo
             enddo
         endif
      endif

! In fvPSAS mode, use file resolver for SST filename templates
#if defined ( fvPSAS_ON )
      if (masterproc) then
        call FileResolv ( trim(job), nymd, nhms, bndtvs, bndtvs )
      endif
#endif

! Check if end of time marching
      if( nymd .ge. nymde .and. nhms .ge. nhmse ) then
         if(masterproc) then
            write(6,*) 'End of Simulation --'
            write(6,*) 'This restart file contains the ending time step of'
            write(6,*) 'the desired simulation, nothing needs to be done.'
         endif
         stop
      endif

#if defined( CCM )

      call mcalday(nymd, nhms, calday)

         mcdate = nymd
         mcsec = nsecf( nhms )
         mscur = mcsec
         mdbase = 0
         msbase = 0

      if(nstep .eq. 0) then
         mbdate  =  nymd
         mbsec   =  nsecf(nhms)
      else
         mbdate  =  nymdb
         mbsec   =  nsecf(nhmsb)
      endif

! Compute nestep for lsm
! Needs to account for shorter than MSTEP runs
 
         nestep = nstep + mstep

         nymdd = nymd
         nhmsd = nhms

      do istep = nstep, nstep + mstep
         if( nymdd .ge. nymde .and. nhmsd .ge. nhmse ) then
             nestep = istep
             goto 666
         endif
         call tick(nymdd, nhmsd, pdt)
      enddo
666   continue

      if (masterproc) then
          write(6,*) ' '
          write(6,*) 'mstep=', mstep
          write(6,*) 'nstep=', nstep
          write(6,*) 'lsm: ending time step= ', nestep
          write(6,*) ' '
      endif

      nlend = .false.

!
      if (.not. diag ) then
         nslice=1
         do k=1, pdiag
            diagattr(1, k) = 0     !not picked in case diag=.false.
         enddo
         allocate( diagbuf(1,nslice,jfirst:jlast) )
      endif

        if (diag) then 
          call diaginit (qname, fldtbl, qunit, qdesc, qvdim, cdiag,
     &                   inx, fldloc, pick, tavg, nslice, nrec2d,
     &                   nrec, n2ddiag, n3ddiag, diagattr, 
     &                   jfirst, jlast)
           allocate( diagbuf(imr,nslice,jfirst:jlast) )
           call init_diagbuf(diagbuf, imr, nslice, jfirst, jlast)
!          gdate2d = nymd
!          gtime2d = nhms      
!  Determine which diagnostic files are needed
          if (ngout2d .ne. -1 .and. n2ddiag .gt. 0) then
            diagfsfc = .true.
          else
            diagfsfc = .false.
          endif
          if (ngout .ne. -1)  then
             if (ngout2d .eq. -1 .and. n2ddiag .gt. 0) then
               diagf = .true.
             else if (n3ddiag .gt. 0) then
               diagf = .true.
             else
               diagf = .false.
             endif
          else
             diagf = .false.
          endif


          if (diagfsfc) then
            call time_shift(nymd, nhms, 0, ngout2d/2, gdate2d, gtime2d)   
            diagflnm2d = 'diag.sfc.bin'
            ctlflnm2d  = 'diag.sfc.ctl'
            ictl2d = 72
            jout2d = 85
          endif
          if (diagf) then
            call time_shift(nymd, nhms, 0, ngout/2, gdate, gtime)
            diagflnm = 'diag.bin'
            ctlflnm  = 'diag.ctl'
            ictl   = 71
            jout   = 84
          endif

          if (diagfsfc) then
            if ( masterproc ) then
            if ( index(diag2dnm, 'diag.sfc.bin') .gt. 0 ) then
            open (jout2d, file=diagflnm2d, form='unformatted',
     .            status='unknown', access='direct',    
     .            recl=imr*jnp*4)
            endif 
            endif 
          endif                                 

          if (diagf) then 
            if ( masterproc ) then
            if ( index(diag3dnm, 'diag.bin') .gt. 0 ) then
            open (jout, file=diagflnm, form='unformatted',
     .            status='unknown', access='direct',
     .            recl=imr*jnp*4)
            endif 
            endif 
          endif

!     ninavg(1) = ' '  		!tphysbc.F

      else
        diagf    = .false.
        diagfsfc = .false.
      endif

      call timing_on('rst_phys')
      call rst_phys(0, iphy, 'p_rst', imr, jnp, nl, icemask, snowh,
     &              lwup, ts, tssub, plevmx, sicthk, 
     &              asdir, asdif, aldir, aldif,
     &              pblht, tpert, qpert,  srfrad,
     &              fsns,  precst, flwds, precc, precl,
     &              sols,  soll,  solsd,  solld, qbeg, tke, 
     &              nc, nstep,  undef, touch, jfirst, jlast )
      call timing_off('rst_phys')

      if(nstep .eq. 0) then
         call drymadj(imr, jnp, nl, jfirst, jlast, ng_d,
     &                moun, ptop, ps, delp, pe, nq, q, 1)
      endif

! oro from (or modified) lsm
! icemask from ic

#if defined( AQUA )
! oro = 0  --> over ocean
! oro = 1  --> over land
! oro = 2  --> over sea ice
       if(nstep .eq. 0) then
       do j=jfirst,jlast
          do i=1,imr
             icemask(i,j)  = 0.
             oro(i,j) =  icemask(i,j)
          enddo
       enddo
       endif
#else
       do j=jfirst,jlast
          do i=1,imr
             oro(i,j) = icemask(i,j)
          enddo
       enddo
#endif
      call initext
#endif

! Get surface geopotential and standard deviation
      call get_srf(imr, jnp, jfirst, jlast, iuhs, srfdata,
     &             moun, phis, sgh, gravit)

        iavg = 0
      if(zstat) then
! Initialize zonal mean statistics
        call init2dz(jnp,  1, jfirst, jlast,psx,   0.)
        call init2dz(jnp, nl, jfirst, jlast, ux,   0.)
        call init2dz(jnp, nl, jfirst, jlast, vx,   0.)
        call init2dz(jnp, nl, jfirst, jlast, tx,   0.)
        call init2dz(jnp, nl, jfirst, jlast, qx,   0.)
        call init2dz(jnp, nl, jfirst, jlast, wx,   0.)
        call init2dz(jnp, nl, jfirst, jlast, pvx,  0.)
        call init2dz(jnp, nl, jfirst, jlast, usus, 0.)
        call init2dz(jnp, nl, jfirst, jlast, usvs, 0.)
        call init2dz(jnp, nl, jfirst, jlast, tsts, 0.)
        call init2dz(jnp, nl, jfirst, jlast, vsts, 0.)
        call init2dz(jnp, nl, jfirst, jlast, vspv, 0.)
      endif

! Cant do this in setup because nymd not known until restart is read
      if (fcstbeg .gt. 0) then
        fcstnymd = fcstbeg
        idummy = 0
        do while (fcstnymd .lt. nymd)
          ndt = fcstfrq * 24 * 3600
          call tick(fcstnymd, idummy, ndt)
        enddo
      else if (fcstbeg .eq. 0) then
        fcstnymd = nymd
      else
        fcstnymd = 0            ! no forecast restarts written
      endif

#if defined ( fvPSAS_ON )
      dasdt = nsecf(fvpsasdt)                ! DAS time step in seconds
#else
      dasdt = 0
#endif

      if ( consv ) then

!$omp parallel do default(shared) private(i,j,k)
      do j=jfirst,jlast
         do i=1,imr
            pk(i,j,1) = ptopk
         enddo

         do k=2,nl+1
            do i=1,imr
               pk(i,j,k) = pe(i,k,j)**cappa
            enddo
         enddo
      enddo

        call pkez(1, imr, nl, jfirst, jlast,
     &            pe, pk, cappa, ks, peln, pkz, .false.)

      endif

      if (masterproc) then
          write(6,*) ' '
          write(6,*) ' Begin time marching:', iord, jord, kord
          write(6,*) ' '
      endif

#if defined ( fvPSAS_ON ) 

! initialize the digital filter related arrays and parameters
       if ( nhms_df .gt. 0 ) then

          mstep = nhms_df/10000 * 3600 / pdt + 1  ! redefine m_step

          mstep_hf = ( mstep + 1 ) / 2
          ncount_df = 1
          call df_Init ( imr, jnp, nl, pdt, cutoff_df, mstep, ierror )
          if ( ierror .ne. 0 ) then
             print *, 'fvpsas: error from df_Init(),  ', ierror
             call exit(1)
          endif

!  accumulate initial conditions
          call df_accum ( imr, jnp, nl, mstep, 
     &                    u, v, pt, q(1:imr,1:jnp,1:nl,1:1), delp,
     &                    ps, ncount_df, ierror )
          if ( ierror .ne. 0 ) then
             print *, 'fvpsas: error from df_accum(),  ', ierror
             call exit(1)
          endif

       endif

#endif

      prog%u    = u
      prog%v    = v
      prog%pt   = pt
      prog%delp = delp
      prog%q    = q

      deallocate( u    )
      deallocate( v    )
      deallocate( pt   )
      deallocate( delp )
      deallocate( q    )

!EOC
      end subroutine stepon_set_
!-----------------------------------------------------------------------


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_finalize --- Finalization
!
! !INTERFACE:
      subroutine stepon_finalize ( prog )
      use fvcore, only: fvcore_finalize
      use gmean, only: gmean_finalize

      implicit none

      type(dyn_prog) :: prog
      integer jout

      deallocate( tsts )
      deallocate( vspv )
      deallocate( vsts )
      deallocate( usvs )
      deallocate( usus )

      deallocate( wx   )
      deallocate( pvx   )
      deallocate( qx   )
      deallocate( tx   )
      deallocate( vx   )
      deallocate( ux   )
      deallocate( psx   )

      deallocate( omga )

#if defined( CCM )
      deallocate( qrs )
      deallocate( qrl )
#endif

      deallocate(  tke )

      call fvcore_finalize

      nullify(u)
      nullify(v)
      nullify(pt)
      nullify(delp)
      nullify(q)

      call prognostics_deallocate ( prog )
      call gmean_finalize

#if defined( CCM )
      if(masterproc) print *, 'FVGCM: Normal Execution.'
#else
      if(masterproc) print *, 'DYCORE: Normal Execution.'
#endif

      jout = 14                             ! SJL: 00.01.02
#if defined( SPMD )
      if (Gid .eq. 0 ) then
#endif
        open(jout, file='status.log')
        write(jout,*) 'Normal Execution.'
        close(jout)

#if defined( SPMD )
      endif
#if (defined PIO)
#if (!defined NO_PIO_BUFFER )
      call deallocate_shm_4gout()
#endif
      call pio_exit()
#endif
#endif


      end subroutine stepon_finalize
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! IROUTINE: stepon_do

! !INTERFACE:
      subroutine stepon_do ( prog )

c$$$#if defined( _TRAJECTORY_ )
c$$$      use m_iostate, only : putstate
c$$$#endif
c$$$
c$$$#if defined( EXTERNAL_TRAJECTORY )
c$$$       use m_iostate, only : getstate
c$$$#endif

      use fvcore, only: fvcore_do, ns
      use fvcore, only: fvcore_tape_rec            ! FastOpt
      use fvcore, only: dynpkg_n2, dynpkg_nsplit   ! FastOpt
      use fvcore, only: ng_c                       ! FastOpt
      use hswf,   only: hswf_do                    ! FastOpt

      implicit none

#if defined( CCM )
#if defined( R4BYTE )
      real(r4), allocatable:: abstot(:,:,:,:) !Non-adjacent layer absorptivites
      real(r4), allocatable:: absnxt(:,:,:,:) !Nearest layer absorptivities
      real(r4), allocatable:: emstot(:,:,:)   !Total emissivity
#else
      real(r8), allocatable:: abstot(:,:,:,:) !Non-adjacent layer absorptivites
      real(r8), allocatable:: absnxt(:,:,:,:) !Nearest layer absorptivities
      real(r8), allocatable:: emstot(:,:,:)   !Total emissivity
#endif
#endif

#if defined ( _TRAJECTORY_ )
      character*255 fnametrj
#endif
#if defined( EXTERNAL_TRAJECTORY )
      character*255 fnametrj
#endif

      logical outd
      logical outg
      logical outg2d
      logical doalb           ! determine if time for albedo computation

      integer iouter, iinner  ! FastOpt
      integer istep           ! FastOpt

      real(r8)  fc
      real(r8)  om            ! FastOpt fixed bug: wrong type kind
      real(r8)  cp            ! FastOpt fixed bug: wrong type kind
      real(r8)  rg            ! FastOpt fixed bug: wrong type kind
      real(r8)  ae            ! FastOpt fixed bug: wrong type kind
      real(r8)  grav          ! FastOpt fixed bug: wrong type kind
      real(r8)  akap          ! FastOpt fixed bug: wrong type kind

      integer, parameter :: nsplt_max = 5    ! FastOpt

#ifdef CHECKPOINTING
#ifdef    EXTERNAL_TRAJECTORY
      integer stat
#endif /* EXTERNAL_TRAJECTORY */
#endif
      integer nymd_trj, nhms_trj
      
#if defined ( CCM )
      allocate( abstot(imr,nl+1,nl+1,jfirst:jlast) )
      allocate( absnxt(imr,nl+1,4,   jfirst:jlast) )
      allocate( emstot(imr,nl+1,jfirst:jlast) )
#endif

      type(dyn_prog), TARGET :: prog

#ifdef CHECKPOINTING
#ifdef    EXTERNAL_TRAJECTORY
!$taf init outer_time = static, nouter
#endif /* EXTERNAL_TRAJECTORY */

! Set pointers
! ------------
       u   => prog%u
       v   => prog%v
      pt   => prog%pt
      delp => prog%delp
      q    => prog%q

!----------------------------------------------------------
! Beginning of basic time step loop
!----------------------------------------------------------
      do iouter = 1, nouter

#ifdef    EXTERNAL_TRAJECTORY
! this store directive breaks the data flow dependency in the adjoint
! it is not intended to actually store something
!$taf store delp,pe,pk,pkz,ps,pt,q,u,v = dummytape, rec=1

! this call provides all prognostic variables
! its data flow is given by directives

!$taf store nymd,nhms             = outer_time, rec=iouter

       nymd_trj = nymd
       nhms_trj = nhms
       call tick(nymd_trj, nhms_trj, pdt)
!$taf store nymd_trj,nhms_trj = outer_time, rec=iouter
       if (trajfrq>0) then
         strajfrq = nsecf(trajfrq)
         snhms    = nsecf(nhms)
         if ( mod(snhms,strajfrq) == 0 ) then
           call strTemplate( fnametrj, trajtmpl,
     &                       xid=trim(job), nymd=nymd, nhms=nhms )
           call getstate( .true., fnametrj, job, nymd_trj, nhms_trj, ptop, ks,
     &                     stat=ierror ) ! intentionally ignore ierror
          endif
       endif

#else

! store all prognostic variables i.e. a full restart state
!$taf store delp,pe,pk,pkz,ps,pt,q,u,v = outer, rec=iouter

#endif /* EXTERNAL_TRAJECTORY */


!$taf init inner = static, ninner

#ifndef  TAPING_IN_DYNPKG
!$taf init dynpkg_1     = static, ninner
!$taf init dynpkg_n2    = static, ninner*dynpkg_n2
!$taf init cd_core_tape = static, ninner*dynpkg_nsplit*dynpkg_n2
!$taf init te_map_tape  = static, ninner
#endif

#ifndef  D_SW_TAPING
!$taf init d_sw_tape1   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*1
!$taf init d_sw_tapej   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*(jnp-1)
#endif

#ifndef  SW_CORE_STORING
!$taf init c_sw_tape1   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*1
!$taf init c_sw_tape2   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*(jnp-1)
#endif

#ifndef  TP_CORE_STORING
!$taf init tpcc1_tape   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*(jnp-1)
!$taf init tpcc2_tape   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*(jnp-1)
#endif

#define  TRAC2D_STORING
#ifndef  TRAC2D_STORING
!$taf init trac2d_tape    = static, ninner*dynpkg_n2
!$taf init trac2d_tape_ns = static, ninner*dynpkg_n2*nsplt_max
#endif

#else    /* CHECKPOINTING */

       iouter = 1

#ifdef    EXTERNAL_TRAJECTORY
! this store directive breaks the data flow dependency in the adjoint
! it is not intended to actually store something
!$taf store delp,pe,pk,pkz,ps,pt,q,u,v = dummytape, rec=1

! the getstate call provides all prognostic variables
! its data flow is given by directives
       nymd_trj = nymd
       nhms_trj = nhms
       call tick(nymd_trj, nhms_trj, pdt)
       if (trajfrq>0) then
         strajfrq = nsecf(trajfrq)
         snhms    = nsecf(nhms)
         if ( mod(snhms,strajfrq) == 0 ) then
           call strTemplate( fnametrj, trajtmpl,
     &                       xid=trim(job), nymd=nymd, nhms=nhms )
           call getstate( .true., fnametrj, job, nymd_trj, nhms_trj, ptop, ks,
     &                     stat=ierror ) ! intentionally ignore ierror
          endif
       endif

#else

! store all prognostic variables i.e. a full restart state
!$taf store delp,pe,pk,pkz,ps,pt,q,u,v = outer, rec=iouter

#endif /* EXTERNAL_TRAJECTORY */

#endif   /* CHECKPOINTING */

#if defined ( _TRAJECTORY_ ) && defined ( GFIO )
       if ( trajfrq .gt. 0 ) then
         strajfrq = nsecf(trajfrq)
         snhms    = nsecf(nhms)
         if ( mod(snhms,strajfrq) == 0 ) then
           call strTemplate( fnametrj, trajtmpl,
     &                       xid=trim(job), nymd=nymd, nhms=nhms )
           call putstate ( fnametrj, 
     &                     job, nymd, nhms, trajfrq, nstep,
     &                     ptop, ks, ak, bk, Ts, oro,
     &                     stat=ierror )
                if ( ierror .ne. 0 ) then
                     print *, 'stepon_do: error from putstate(),  ', ierror
                    call exit(1)
                endif
         endif
       endif
#endif

      do iinner = 1, ninner

       istep =  iinner    + (iouter-1)*ninner
       nstep = (iinner-1) + (iouter-1)*ninner

       fvcore_tape_rec  = iinner-1 + stepon_tape_rec*mstep

!$taf store delp,pe,pk,pkz,pt,u,v = inner, rec=fvcore_tape_rec+1
!$taf store omega                 = inner, rec=fvcore_tape_rec+1

! ************************************
! Calling Finite-Volume Dynamical Core
! ************************************
      call timing_on('FVCORE')
      om = omega                            ! FastOpt convert type
      cp = cpair                            ! FastOpt convert type
      rg = rair                             ! FastOpt convert type
      ae = rearth                           ! FastOpt convert type
      call fvcore_do(imr, jnp, nl, nc, jfirst, jlast,
     &            ng_d, ng_s, nq, ps, pe,
     &            delp, u, v, pt, q, pk, pkz, phis, nsplit, pdt,
     &            ptop, om, cp, rg, ae, iord, jord, kord, 
     &            umax,  omga, peln, consv, convt)
      call timing_off('FVCORE')

      nstep = nstep + 1
      call    tick(nymd, nhms, pdt)
      call mcalday(nymd, nhms, calday)
      mcdate = nymd
      mcsec = nsecf( nhms )
      mscur = mcsec

#if defined ( _TRAJECTORY_ )
       if ( trajfrq .gt. 0 ) then
       if ( mod(nsecf(nhms),nsecf(trajfrq)) .eq. 0 ) then
            call strTemplate( fnametrj, trajtmpl,
     &                        xid=trim(job), nymd=nymd, nhms=nhms )
            call putstate ( fnametrj, 
     &                      job, nymd, nhms, trajfrq, nstep,
     &                      ptop, ks, ak, bk, Ts, oro,
     &                      stat=ierror )
                 if ( ierror .ne. 0 ) then
                      print *, 'stepon_do: error from putstate(),  ', ierror
                      call exit(1)
                 endif
       endif
       endif
#endif

#if defined( CCM )

      if(istep .eq. mstep) nlend = .true.

! Determine whether it is time for a shortwave or longwave radiation
! calculation

      dosw =  nstep.eq.1    .or.   iradsw.eq.1
     $                      .or. mod(istep-1,iradsw).eq.0

      dolw =  nstep.eq.1    .or.   iradlw.eq.1
     $                      .or. mod(istep-1,iradlw).eq.0

! For lsm:
! Albedos over ocean computed every time step.
      doalb = iradsw.eq.1   .or.    mod(istep,iradsw).eq.0

! Determine whether it is time for an absorptivity/emissivity calculation
      doabsems = nstep.eq.1   .or.   iradae.eq.1 .or.
     $           mod(istep-1,iradae).eq.0

! If aeres true, a/e data will be stored on restart file; not implemented
      aeres = mod(istep-1,iradae) .ne. 0
 
! Always do annual cycle
      if( .not. anncyc ) then
          if(masterproc) write(6,*) 'perpetual run is not implemented.'
          stop
      endif

! Update ozone data on shortwave or longwave time step.
! Note that the ozone data is not needed on a longwave time step unless the
! absorptivities are being updated ("doabsems").
 
      if (dosw .or. dolw) then
        call oznint
      end if

      call sstint

! *********************************
! Calling NCAR CCM physics
! *********************************
! On input pt is virtual temp. becomes virtual pt on output.

       call timing_on('PHYSDRV')
       call physdrv(abstot, absnxt, emstot, qrs, qrl, diagattr, diagbuf,
     &            tdt, bdt,  nq, u,  v,  pt, q,  doalb,
     &            pk, pe, delp, ps, precp, rayf, co2f, omga, peln,
     &            phis , sgh,  pblht, tpert, qpert,
     &            pkz,  coslon, sinlon , nymd, nhms,
     &            ih2o, h2odata, dcaf, gwdf, qbeg, tke,
     &            jfirst, jlast, ng_d, ng_s, nslice)

       call timing_off('PHYSDRV')

#else

! Adiabatic case
       if( force_dycore ) then
!$taf store u,v,pt,pe,pkz = inner, rec=fvcore_tape_rec+1
         rg   = rair
         grav = gravit
         akap = cappa
         call hswf_do(imr, jnp, nl, jfirst, jlast,
     &            u, v, pt, pe, pkz, pdt, akap, grav, rg, dcaf,
     &           .true., rayf, sinp,cosp,sine,cose, coslon, sinlon,
     &            ng_s, ng_d)
       endif

#endif

      if( aoa ) then
!$taf store pe                     = inner, rec=fvcore_tape_rec+1
          model_time = nstep*pdt
          call age_of_air(imr, jnp, nl, jfirst, jlast, ng_d, model_time,
     &                    pe,  q(1,jfirst-ng_d,1,nq) )
      endif


      if( idtout .eq. 0) then
        outd = .false.
      else
        outd = mod(mcsec,idtout) .eq. 0
      endif

      if( igtout2d .le. 0) then
        outg2d = .false.
      else
        outg2d = mod(mcsec,igtout2d) .eq. 0
      endif                                         
      if( igtout .eq. 0) then
        outg = .false.
      else
        outg = mod(mcsec,igtout) .eq. 0
      endif

      if ( (nymd .gt. nymd_psas) .OR.
     &     (nymd .eq. nymd_psas .AND. nhms .gt. nhms_psas) ) then
           forecast = .true.
      else
           forecast = .false.
      end if

      if ( dasdt .eq. 0 .or. forecast) then
         dofvpsas = .false.
      else
         dofvpsas = mod(mcsec,dasdt) .eq. 0
      end if

#if defined ( fvPSAS_ON ) 

! Subject dynamical fields to the digital filter and accumulat them
       if ( nhms_df .gt. 0 ) then  !and accumulat them

          ncount_df = ncount_df + 1
          call df_accum ( imr, jnp, nl, mstep, 
     &                    u, v, pt, q(1:imr,1:jnp,1:nl,1:1), delp,
     &                    ps, ncount_df, ierror )

! Absolutely no output or analysis during Digital Filter accumulation phase

          dofvpsas = .false.
          outd     = .false.
          outg2d   = .false.
          outg     = .false.
          forecast = .false.
          fcstfrq  = 0

       endif

#endif

#if defined( CCM )

        call timing_on('DIAGOUT')
      if (outg2d) then

        if (ncar .and. diagfsfc) then
          ! NOTE: the fvPSAS script relies on the next statement
          if(masterproc) write(*,'(a,I10.8,i10.6)') 'diagout1_time:', nymd, nhms
#if defined (GFIO)
          call strTemplate( diagflnm2d, diag2dnm, 
     &                      xid=trim(job), nymd=nymd, nhms=nhms)
#endif
          call diagout(jout2d, qname, qvdim, cdiag, fldloc, pick,
     &         tavg, 1, pd2d, nrec2d, diagattr, diagbuf, jfirst, jlast,
     &         nslice,diagflnm2d, undef, nymd, nhms, igtout2d, pref_rev,
     &         qunit, qdesc)
        endif

      endif

      if (outg) then

        if (ncar .and. diagf) then
          if (.not.diagfsfc .and. n2ddiag .gt. 0) then
            diagstart = 1
          else
            diagstart = pd2d+1
          endif
          diagend = pd2d + pd3d          

          ! NOTE: the fvPSAS script relies on the next statement
          if(masterproc) write(*,'(a,I10.8,i10.6)') 'diagout2_time:', nymd, nhms
#if defined (GFIO)
          call strTemplate( diagflnm, diag3dnm, 
     &                      xid=trim(job), nymd=nymd, nhms=nhms )
#endif
          call diagout(jout, qname, qvdim, cdiag, fldloc, pick, tavg,
     &                 diagstart, diagend, nrec, diagattr,
     &                 diagbuf, jfirst, jlast, nslice,    
     &                 diagflnm, undef, nymd, nhms, igtout, pref_rev,
     &                 qunit, qdesc)

        endif
      endif

      call timing_off('DIAGOUT')
#endif

#ifdef REFOUT_FASTOPT
      if (outd) then

        call timing_on('REFOUT')

#if defined (GFIO)
        call strTemplate( progflnm, prog3dnm, 
     &                    xid=trim(job), nymd=nymd, nhms=nhms )
#endif
        call refout(imr, jnp, nl, jfirst, jlast, ccm_on,
     &              iavg, pe,    delp , u,  v, pt, q, phis, pkz,
     &              peln, precp, nq ,  rair, coslon, sinlon, rearth,
     &              iout, nymd,  nhms,  omga,  psx,
     &              ux, vx, tx, qx, wx, pvx, usus, usvs, tsts, vsts,
     &              vspv, zstat, ts, zvir, gravit, tke, ng_s, ng_d,
     &              progflnm, undef, idtout, pref_rev )

         call drymadj(imr, jnp, nl, jfirst, jlast, ng_d,
     &                moun, ptop, ps, delp, pe, nq, q, 0)

        call timing_off('REFOUT')
      endif
#endif

#if defined ( fvPSAS_ON )

! ******************************************************
! Calling NASA/DAO PSAS in Rapid Update Cycle (RUC) mode 
! ******************************************************

! Note: by choosing not to correct TE consevation the interface is greatly
!       simplified (no update of pe, pk, and pkz in the lowest three 
!       layers are needed). Otherwise, update pe, pk, and pkz here.

        if ( dofvpsas ) then

! SJL:
! Interface to pSAS needs to be rewritten for MLP/MPI
! Passing in local arrays for now
!

#if defined ( DIAPER )
! JDC: save the first guess pt
!
           allocate( ptold(imr,jfirst-ng_d:jlast+ng_d,kfirst:klast) )  !virture potential temperature

!$omp parallel do private(i,j,k)

           do k = 1, nl
             do j = jfirst, jlast
               do i = 1, imr
                 ptold(i,j,k)   = pt(i,j,k)
               enddo
             enddo
           enddo
#endif

           call fvPSAS ( job,     nymd,   nhms,     fvpsasdt, nstep,
     &                   imr,     jnp,    nl, 
     &                   jfirst,  jlast,  ng_s,     ng_d,
     &                   ptop,    pint,   ks,       ak,       bk,    
     &                   phis,    sgh,    ts,       oro,      ps,
     &                   delp,    u,      v,        pt,       q, 
     &                   delpmin, nremap, prec_das, ierror ) 
     
           if ( ierror .ne. 0 ) then
               print *, 'fvpsas: error from psasdrv, rc ', ierror
               call exit(1)
           else
             if(masterproc) write(6,*) 'psasdrv successfully called at '
     &                                  ,nymd,'/',nhms
           endif

#if defined ( DIAPER )
!
! JDC: call moist physics to avoid excessive precipitation  after fvPSAS
!
           call moistdrv(
     &          nstep,   bdt,     coslon,  sinlon,
     &          ptop,    phis,    pblht,   ps,      ptold,
     &          tpert,   qpert,   delp,    u,       v,
     &          pt,      q,
     &          jfirst,  jlast,   ng_d,    ng_s)

           deallocate( ptold )
! End of moist physics adjustment

#endif
!
! Remapping after analysis
!

           if( delpmin < 10.0 ) then
             if(masterproc) then
               write(6,*) ' '
               write(6,*) 'WARNING: Remmaping is needed after analysis'
               write(6,*) 'Check restart file and observations to see ',
     &                    'why this happen.'
               write(6,*) 'Total number need to remap is: ', nremap
               write(6,*) 'The minimum lowest layer thickness is: ', delpmin
               write(6,*) ' '
             end if

!$omp  parallel do default(shared) private(i,j,k)

             do j = jfirst, jlast
               do i = 1, imr
                 pe(i,1,j) = ptop
                 pk(i,j,1) = ptopk
               enddo

               do k = 2, nl+1
                 do i = 1, imr
                   pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
                   pk(i,j,k) = pe(i,k,j) ** cappa
                 enddo
               enddo
             enddo

             nx = 4
! te0 is not used because consv = .false.
             te0 = 1.0
             call te_map(.false., .false., ps, omga, pe, delp, pkz, pk, mpdt,
     &                  imr, jnp, nl, nx, jfirst, jlast, nq,  u,  v,
     &                  pt, q, phis,  cpair, cappa, kord, peln, te0,
     &                  ng_d, ng_s, 0 )
           endif

! End of remapping

        endif

        if (fcstfrq .gt. 0) then
          if (fcstnymd .lt. nymd) then
            ndt = fcstfrq * 24 * 3600
            idummy = 0
            call tick(fcstnymd, idummy, ndt)
          endif
          if (nymd .eq. fcstnymd) then
            if ((nymd .eq. nymde  .and. nhms .eq. nhmse ) .or.
     &         (istep .eq. mstep)) then
            else

              hour = int(real(nhms)/10000.)
              do n = 1, 24              ! 24 hours in the day
                if (nhms .eq. fcsthrs(n)*10000) then

! Check global mean dry air mass; non-conservation due to
! machine rounding.
                  call drymadj(imr, jnp, nl, jfirst, jlast, ng_d,
     &                  moun, ptop, ps, delp, pe, nq, q, 0)
 
                  call strTemplate ( drstflnm,
     &              '%s.rst.lcv.%y4%m2%d2_%h2z.bin', 'GRADS',
     &                               trim(job), nymd, nhms, ier )
#ifdef BINRST
                  call timing_on('writerst')
                  call writerst(iuic, drstflnm, imr, jnp, nl, nq, 
     &                 jfirst, jlast, nstep, nymd, nhms, ps, delp,
     &                 u, v, pt, q, ng_s, ng_d)
                  call timing_off('writerst')
#endif

#if defined( CCM )
                  call strTemplate ( prstflnm,
     &             '%s.rst.phys.%y4%m2%d2_%h2z.bin', 'GRADS',
     &                             trim(job), nymd, nhms, ier )
                  call rst_phys(1, iphy, prstflnm, imr, jnp, nl, oro,
     &                        snowh, lwup, ts, tssub, plevmx, sicthk,
     &                        asdir, asdif, aldir, aldif,
     &                        pblht, tpert, qpert,  srfrad,
     &                        fsns,  precst, flwds, precc, precl,
     &                        sols,  soll,  solsd,  solld, qbeg, tke,
     &                        nc, nstep , undef,  touch, jfirst, jlast )

                  call strTemplate ( lrstflnm,
     &                 '%s.rst.lsm.%y4%m2%d2_%h2z.bin', 'GRADS',
     &                             trim(job), nymd, nhms, ier )
#if (! defined CLM2)
                  call rst_lsm(lrstflnm, ilsm, nstep)
#endif

#endif

                endif
              enddo
            endif
          endif
        endif


#endif

#if defined ( fvPSAS_ON ) 

! write out rst_lsm and rst_phys at the half of the filter window steps.
        if ( nhms_df .gt. 0 .and. mstep_hf .eq. ncount_df) then

           nymd_hf = nymd
           nhms_hf = nhms
           nstep_hf = nstep
           call rst_phys(1, iphy, 'p_rst.df', imr, jnp, nl, oro,
     &          snowh, lwup, ts, tssub, plevmx, sicthk,
     &          asdir, asdif, aldir, aldif,
     &          pblht, tpert, qpert,  srfrad,
     &          fsns,  precst, flwds, precc, precl,
     &          sols,  soll,  solsd,  solld, qbeg, tke,
     &          nc, nstep , undef,  touch, jfirst, jlast )

#if (! defined CLM2)
           call rst_lsm('lsm_rst.df', ilsm, nstep)
#endif

        endif

#endif
       
crg      if(nymd .ge. nymde  .and. nhms .ge. nhmse )  goto 123

!----------------------------------------------------------
! Ending time step logic
!----------------------------------------------------------
#ifdef CHECKPOINTING
      end do
#endif
      end do

123   if(masterproc) write(6,*) 'End of time marching'

       if( aoa ) then
!$taf store q(:,:,:,nq),nstep                 = onetape, rec=stepon_tape_rec+1
          model_time = nstep*pdt
          call aoa_wrt(imr, jnp, nl, jfirst, jlast, ng_d, model_time,
     &                 q(1,jfirst-ng_d,1,nq), aout )
      endif

#ifdef REFOUT_FASTOPT

#if defined( CCM )

        call timing_on('DIAGOUT')

      if (igtout2d .eq. 0 .and. ngout2d .ge. 0) then

        if (ncar .and. diagfsfc) then
          ! NOTE: the fvPSAS script relies on the next statement
          if(masterproc) write(*,'(a,I10.8,i10.6)') 'diagout1_time:', nymd, nhms
#if defined (GFIO)
          call strTemplate( diagflnm2d, diag2dnm, 
     &                      xid=trim(job), nymd=nymd, nhms=nhms )
#endif
          call diagout(jout2d, qname, qvdim, cdiag, fldloc, pick,
     .         tavg, 1, pd2d, nrec2d, diagattr, diagbuf, jfirst, jlast,
     .         nslice,diagflnm2d, undef, nymd, nhms, igtout2d, pref_rev,
     .         qunit, qdesc)
        endif

      endif

      if (igtout .eq. 0 .and. ngout .ge. 0) then

        if (ncar .and. diagf) then
          if (.not.diagfsfc .and. n2ddiag .gt. 0) then
            diagstart = 1
          else
            diagstart = pd2d+1      
          endif
          diagend = pd2d + pd3d    
          ! NOTE: the fvPSAS script relies on the next statement
          if(masterproc) write(*,'(a,I10.8,i10.6)') 'diagout2_time:', nymd, nhms
#if defined (GFIO)
          call strTemplate( diagflnm, diag3dnm, 
     &                      xid=trim(job), nymd=nymd, nhms=nhms )
#endif
          call diagout(jout, qname, qvdim, cdiag, fldloc, pick, tavg,
     &                 diagstart, diagend, nrec, diagattr,
     &                 diagbuf, jfirst, jlast, nslice,       
     &                 diagflnm, undef, nymd, nhms, igtout, pref_rev,
     &                 qunit, qdesc)

        endif  
      endif  

        call timing_off('DIAGOUT')

#endif

      if (idtout .gt. 0 .and. ndout .ge. 0) then

        call timing_on('REFOUT')

#if defined (GFIO)
        call strTemplate( progflnm, prog3dnm, 
     &                    xid=trim(job), nymd=nymd, nhms=nhms )
#endif
        call refout(imr, jnp, nl, jfirst, jlast, ccm_on,
     &              iavg, pe,    delp , u,  v, pt, q, phis, pkz,
     &              peln, precp, nq ,  rair, coslon, sinlon, rearth,
     &              iout, nymd,  nhms,  omga, psx,
     &              ux, vx, tx, qx, wx, pvx, usus, usvs, tsts, vsts,
     &              vspv, zstat, ts, zvir, gravit, tke, ng_s, ng_d,
     &              progflnm, undef, idtout, pref_rev )


        call timing_off('REFOUT')

      endif

! Write restart file for the dynamics.

#if ( fvPSAS_ON ) 

!    output (filtered) dynamics restart

       if ( nhms_df .gt. 0 ) then  

          nymd = nymd_hf
          nhms = nhms_hf
          nstep = nstep_hf
          call df_update (imr, jnp, nl,
     &         u, v, pt, q(1:imr,1:jnp,1:nl,1:1), delp, ps, ierror )
          call df_clean()

#ifdef BINRST
          call timing_on('writerst')
          call writerst(iuic, 'd_rst.df', imr, jnp, nl, nq, jfirst, jlast,
     &                  nstep, nymd, nhms, ps, delp, u, v, pt, q,
     &                  ng_s, ng_d)
          call timing_off('writerst')
#endif

       else

#endif

! Check global mean dry air mass; non-conservation due to
! machine rounding.

#ifdef BINRST
          call timing_on('writerst')
          call writerst(iuic, 'd_rst', imr, jnp, nl, nq, jfirst, jlast,
     &                  nstep, nymd, nhms, ps, delp, u, v, pt, q,
     &                  ng_s, ng_d)
          call timing_off('writerst')
#endif

#if defined( CCM )
! Write restart file CCM phys.
         call timing_on('rst_phys')
         call rst_phys(1, iphy, 'p_rst', imr, jnp, nl, oro, snowh, lwup,
     &                 ts, tssub, plevmx, sicthk,
     &                 asdir, asdif, aldir, aldif,
     &                 pblht, tpert, qpert,  srfrad,
     &                 fsns,  precst, flwds, precc, precl,
     &                 sols,  soll,  solsd,  solld, qbeg, tke,
     &                 nc, nstep , undef,  touch, jfirst, jlast )
         call timing_off('rst_phys')

        if (diagfsfc) then
           flsz = real(nrec2d-1)*imr*jnp*4
           if(masterproc) print *, 
     &                   'Expected surface diagnostic file size = ', flsz
           nt = (nrec2d - 1) / n2ddiag
           do k = 1, nl
             pref(nl-k+1) = .01 * .5 * (ak(k) + ak(k+1) +
     .                      (bk(k) + bk(k+1)) * psref)
           end do

         call grads_ctl (ictl2d, ctlflnm2d, diagflnm2d, job, imr, jnp,
     .                  1, nt, undef, gdate2d, gtime2d, igtout2d,
     .                  1, pd2d,
     .                  pdiag, qname, pick, pref, qvdim, qunit,
     .                  qdesc)
       endif

       if (diagf) then               
           flsz = real(nrec-1)*imr*jnp*4
           if(masterproc) print *, 
     &                   'Expected eta diagnostic file size = ', flsz
           if (diagfsfc) then
             nt = (nrec - 1) / (n3ddiag*nl)
           else
             nt = (nrec - 1) / nslice
           endif               
           do k = 1, nl
             pref(nl-k+1) = .01 * .5 * (ak(k) + ak(k+1) +
     .                      (bk(k) + bk(k+1)) * psref)
           end do
           if (.not.diagfsfc .and. n2ddiag .gt. 0) then
             diagstart = 1
           else
             diagstart = pd2d+1
           endif
           diagend = pd2d + pd3d
           call grads_ctl (ictl, ctlflnm, diagflnm, job, imr, jnp,
     .                    nl, nt, undef, gdate, gtime, igtout,
     .                    diagstart, diagend,
     .                    pdiag, qname, pick, pref, qvdim, qunit,
     .                    qdesc)
        endif    

#endif

      if(zstat .and. iavg .gt. 0) then
        call zout(izou, iavg, jnp, nl, jfirst, jlast, psx, ux, vx,
     &            tx, qx, wx, pvx, usus, usvs, tsts, vsts, vspv)
      endif


#if defined ( fvPSAS_ON ) 
      end if    ! doing digital filter
#endif

#endif  /* REFOUT_FASTOPT */

      call timing_off('TOTAL')

#if defined( CCM )
      deallocate( abstot )
      deallocate( absnxt )
      deallocate( emstot )
#endif

!EOC
      end subroutine stepon_do
!-----------------------------------------------------------------------

      end module stepon
!-----------------------------------------------------------------------
