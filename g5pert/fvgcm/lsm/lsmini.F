#include <preproc.h>
      subroutine lsmini (londim  ,latdim  ,beglatspmd,endlatspmd,
     &                   nlonatm ,oroatm  ,xmdbase ,xmsbase ,
     &                   xmbdate ,xmbsec  ,xnstep  ,dtatm   ,
     &                   xnsrest ,xnestep ,xnhtfrq, xmfilt , 
     &                   xninavg, xirt    ,xirad   ,lwupxy  ,tsxy  ,
     &                   asdirxy ,asdifxy ,aldirxy ,aldifxy ,
     &                   snowxy, xmcdate, xmcsec)

#if defined (SPMD)
      use mod_comm, only: gid, mp_bcst_int, mp_scat3d, mp_scat3d_int
#endif

#include <vartyp.h>
#include <lsmpar.h>
#include <basinpar.h>

* ------------------------ code history --------------------------------
* source file:       lsmini.F
* purpose:           lsm initialization. return land surface variables to
*                    atm model if initial run
* date last revised: August 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      
* reviewed:         
* ----------------------------------------------------------------------

* ------------------------ notes ---------------------------------------
* initialization routine for land surface model. 
*
* the land surface model works by gathering all the land points on a
* [lsmlon] x [lsmlat] grid into a vector of [lpt] land points. this is
* then expanded into a "big" vector of [kpt] subgrid points, allowing for
* up to [msub] subgrid points per land point. the [kpt] subgrid points
* are processed as [numlv] "little" vectors of [numkpt] points for [numlv] 
* calls to the vector-based land surface model. 
*
* model constants are set in block data subprograms. this subroutine: 
*   o determines time-invariant surface data on [lsmlon] x [lsmlat] grid 
*     if initial run with arbitrary initialization this dataset is
*     created at run time, otherwise it is obtained from the initial 
*     or restart dataset
*   o initializes run control variables for the [lsmctl], [lsmhis] and [lsmhisdyn]
*     common blocks via the [lsmexp] namelist 
*   o defines the multiple plant types and fraction areas for each surface type 
*   o builds the appropriate subgrid <-> grid mapping indices  
*   o assigns subgrid points the appropriate time-invariant data. i.e., 
*     initializes the time constant [lsmtc] common block
*
* for an initial run only:
*   o reads initial data and initializes the time variant [lsmtv] and [lsmtvdyn] 
*     common blocks
*   o returns required land surface variables to atm model
*   o writes variables in the [lsmtc] and [lsmtcdyn] common blocks to a 
*     history file (branch also)
*   o initializes history file variables
*   o reads restart data if a continuation run 
*
* ----------------------------------------------------------------------

* ------------------------ input/output variables ----------------------
* input from atmospheric model

      integer londim            !longitude points
      integer latdim            !latitude points
      integer beglatspmd    !beg lat on proc, degenerates to 1 for nonspmd
      integer endlatspmd    !end lat on proc, degenerates to atmlat for nonspmd
#ifdef COUP_CCM
      integer nlonatm(latdim)   !number of longitudes per latitude band from atm
      real oroatm(londim,latdim)!oro flag from atm, 0 = non-land. 1 = land 
      integer xmdbase           !base day of run (e.g., 0)
      integer xmsbase           !base seconds of base day (e.g., 0)
      integer xmbdate           !base date of run (yymmdd format) (e.g., 000901)
      integer xmbsec            !base seconds of base date (e.g., 0)
      integer xnstep            !time step index
      integer xnestep           !ending time step index
      integer xnhtfrq           !history interval(iterations)
      integer xmfilt            !max number of time samples per history file
      integer xirad             !solar radiation frequency (iterations)
      character*1 xninavg       !equals "q" or "Q" if using monthly a
      integer xirt              !mss retention time
      integer xnsrest           !0 = initial run, > 0 = continuation run
      real dtatm                !time step (s)
#endif

* output to atmospheric model
 
      real asdirxy(londim,latdim)   !albedo - visible waveband, direct
      real asdifxy(londim,latdim)   !albedo - visible waveband, diffuse
      real aldirxy(londim,latdim)   !albedo - near infrared waveband, direct
      real aldifxy(londim,latdim)   !albedo - near infrared waveband, diffuse
      real tsxy   (londim,latdim)   !surface (skin) temperature
      real lwupxy (londim,latdim)   !emitted longwave radiation (w/m**2)
      real snowxy (londim,latdim)   !water equivalent snow (m)
* ----------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <dynmem.h>
#include <lsmtc.h>
#include <lsmio.h>
#include <vegcon.h>
#include <soicon.h>
#include <snocon.h>
#include <vegtyp.h>
#include <phycon.h>
#include <lsmctl.h>
#include <lsmhis.h>
#ifndef COUP_CCM
#include <cpl.h>
#endif
#include <lsmspmd.h>
#include <basin.h>

! SJL
      integer xmcdate, xmcsec
* ----------------------------------------------------------------------
 
* ------------------------ local variables -----------------------------

      integer i,j,k,l,m,n             !indices
      integer getavu                  !next available unit number
      integer nflnd                   !fractional land data data unit number
      integer nfini                   !initial data unit number
      integer nfpnt                   !pointer to restart data unit number
      integer nfrst                   !restart data unit number
      integer nflon                   !number of lons per lat unit number
      integer nio                     !generic unit number
      integer minlon                  !minimum longitude points
      integer maxlon                  !maximum longitude points
      integer veg(lsmlon,lsmlat,msub) !vegetation type
      real wt(lsmlon,lsmlat,msub)     !subgrid weights
      real pi                         !pi
      logical crtlsf                  !true if need to create land surface data

      integer kyr                     !current year (0, ...)
      integer kmo                     !current month (1 -> 12)
      integer kda                     !current day (1 -> 31)
      real calday                     !current gmt calendar day (1.00 -> 365.00)

      integer lencs                   !integer function: length of character string
      character*80 fname              !full pathname of lsm restart pointer file
      character*80 fnamer             !full pathname of lsm restart file
      character*80 fnameh             !full pathname of lsm history file
      character*80 locfn              !local file name

      real wedge(lsmlon+1,lsmlat)     !grid cell longitude, western edge  (degrees) 
      real sedge(lsmlat+1)            !grid cell latitude , southern edge (degrees) 
      real clon (lsmlon  ,lsmlat)     !grid cell longitude, center (degrees) 
      real clat (lsmlat)              !grid cell latitude , center (degrees) 

#ifndef COUP_CCM
      real clon_flxc(lsmlon)          !grid cell longitude, center, sent to coupler 
      real wedge_flxc(lsmlon+1)       !grid cell longitude, western edge, sent to coupler  
      integer mland (lsmlon,lsmlat)   !land mask sent to coupler
#endif

#if ( !defined SPMD )
      integer gid
#endif
      integer numcpu
      character*80 evalue
      integer nhms, nhmsf, nn
      nhmsf(nn)  = nn/3600*10000 + mod(nn,3600 )/ 60*100 + mod(nn, 60)


* ----------------------------------------------------------------------

#if defined ( SPMD )
      if( gid == 0 ) then
          masterproc = .true.
      else
          masterproc = .false.
      endif
#else
          masterproc = .true.
          gid = 0
#endif

* begin executable code: echo initialization to standard output

      if (masterproc) then
         write (6,*)
         write (6,*) '***************************************************'
         call lsmhdr
         write (6,*) 'attempting to initialize LSM'
         write (6,*)
      end if		

* ----------------------------------------------------------------------
* initialize run control variables
* ----------------------------------------------------------------------

#ifdef COUP_CCM
      call lsmctli (xmdbase, xmsbase, xmbdate, xmbsec,
     &              xnsrest, xnstep , xnestep, dtatm ,
     &              xnhtfrq, xmfilt , xninavg, xirt  ,
     &              xirad )
#else
      call lsmctli
#endif

* ----------------------------------------------------------------------
* initialize fortran unit numbers 1 to 99 to inactive: except standard 
* input (5) and standard output (6) 
* ----------------------------------------------------------------------

      do i = 1, 99
         lsmiou(i) = .false.
      end do
      lsmiou(5) = .true.
      lsmiou(6) = .true.

      if (masterproc) then
         write (6,*) 'preset fortran unit numbers:'
         write (6,*) '   unit  5 = standard input'
         write (6,*) '   unit  6 = standard output'
         write (6,*)
      end if             ! (if masterproc) 

* ----------------------------------------------------------------------
* set mass store pathnames for history and restart files
* ----------------------------------------------------------------------

      call msspni

* ----------------------------------------------------------------------
* determine lsm surface data on [lsmlon] x [lsmlat] grid.
* 
* if run is coupled to ccm, then this grid (size, lat and lon indexing) 
* must be the same as the atmospheric grid. this constraint no longer 
* holds if the land model runs as a separate executable through the 
* flux coupler. a brief description follows:
* if initial run:
*  o if arbitrary initialization the surface dataset is calculated
*    and is written out as diagnostic output (see below)
*  o if the run is initiated from an initial dataset, the dataset is
*    opened and the surface data part of the dataset is read here.
*    the remaining part of the dataset will be read in subroutine lsmtvi
* if restart or branch run:
*  o get appropriate restart and history file names if continuation run. 
*  o restart: the file lsmr_rstrt contains the full mass store or 
*    local names of the current restart file lsmr_xxxx and history 
*    file lsmh_xxxx
*  o branch: the variable nrevsn is the full mass store or local name
*    of the appropriate restart file lsmr_xxxx. do not need to have 
*    history file name because will open a new history file starting 
*    with lsmh_0001
*  o read surface data part of the dataset, subroutine lsmres will then 
*    read the rest of the restart file and initialize the lsmtvdyn 
*    common block
* ----------------------------------------------------------------------

* determine if lsm surface data needs to be created

      if (nsrest.eq.0  .and. finidat.eq.'arbitrary initialization') then
         crtlsf = .true.
      else
         crtlsf = .false.
      end if

* if lsm surface data already exists then obtain it and read it

      if (.not. crtlsf) then
         if (masterproc) then

* obtain file containing lsm surface data

            if (nsrest.eq.0 ) then  ! initial run
               nfini = getavu()
               call getfil (finidat, locfn)
               call opnfil (locfn, nfini, 'u')
               nio = nfini
            else
               if (nsrest.eq.1) then ! restart run
                  fname = rpnt_dir(1:lencs(rpnt_dir)) //'/'
     &                    //rpnt_fil(1:lencs(rpnt_fil))
                  nfpnt = getavu()
                  call opnfil (fname, nfpnt, 'f')
                  read (nfpnt,'(a80)') fnamer
                  read (nfpnt,'(a80)') fnameh
                  call relavu (nfpnt)
               else                 ! branch run
                  fnamer = nrevsn
                  fnameh = ' '
               end if
               nfrst = getavu()
               call getfil (fnamer, locfn)
               call opnfil (locfn, nfrst, 'u')
               nio = nfrst
            endif

*     read lsm surface data

            read(nio) ( numlon(j)                      ,j=1,lsmlat)
            read(nio) ((latixy(i,j)     ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) ((longxy(i,j)     ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) ((surf2d(i,j)     ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) ((soic2d(i,j)     ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) ((sand2d(i,j)     ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) ((silt2d(i,j)     ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) ((clay2d(i,j)     ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) ((pctlak(i,j)     ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) ((pctwet(i,j)     ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) ((drngarea(i,j)   ,i=1,numlon(j)),j=1,lsmlat)
            read(nio) (((drngfrac(n,i,j),n=1,ndrnmax),i=1,numlon(j)),
     $           j=1,lsmlat)
            read(nio) (((drngbasn(n,i,j),n=1,ndrnmax),i=1,numlon(j)),
     $           j=1,lsmlat)

         endif                  !if masterproc
#if defined(SPMD)
         call mp_scat3d_int(numlon, 1, lsmlat, 1, 1,1, 1,lsmlat, 1,1, 0)
         call mp_scat3d(latixy, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d(longxy, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d_int(surf2d, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d_int(soic2d, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d(sand2d, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d(silt2d, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d(clay2d, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d(pctlak, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d(pctwet, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d(drngarea, lsmlon, lsmlat, 1, 1,lsmlon, 1,lsmlat, 1,1, 0)
         call mp_scat3d(drngfrac, ndrnmax, lsmlon, lsmlat, 1,ndrnmax, 1,lsmlon, 1,lsmlat, 0)
         call mp_scat3d_int(drngbasn, ndrnmax, lsmlon, lsmlat, 1,ndrnmax, 1,lsmlon, 1,lsmlat, 0)
#endif
      endif

* ----------------------------------------------------------------------
* determine lsm grid only if need to great lsm surface data 
* ----------------------------------------------------------------------

      if (crtlsf) then

#ifndef COUP_CCM
* determine number of longitudes for each latitude strip from input data
* only if running through flux coupler 
* ***NOTE *** for now we are HARD-WIRING in a rectangular grid. If at
* a future time a non-uniform grid will be implemented the following 
* commented code can be used.
         
          do j=1,lsmlat
             numlon(j) = lsmlon
          end do

!         if (masterproc) then
!           nflon = getavu()
!           call getfil (flondat, locfn)
!           call opnfil (locfn, nflon, 'f')
!           do j = 1, lsmlat
!              read (nflon,*) numlon(j)
!           end do
!           call relavu (nflon)
!         end if
#else
* set up number of longitudes per latitude band for lsm to be the 
* same as that passed in by the atm

         do j=1,lsmlat
            numlon(j) = nlonatm(j)
         end do
#endif
         if (masterproc) then

* determine lsm cell in degrees

            call cell_box(lsmgeo, lsmlat, lsmlon, numlon, wedge, sedge,
     &                    clon  , clat)

* determine latitude and longitude in radians for lsm routines

            pi = 4.*atan(1.) 
            do j = 1, lsmlat
               do i = 1, numlon(j)
                  latixy(i,j) = clat(  j)*pi/180.
                  longxy(i,j) = clon(i,j)*pi/180.
               end do
            end do

         end if                 ! end of maseterproc block
      endif                     ! end of crtlsf block

* ----------------------------------------------------------------------
* create lsm surface data if initial run with arbitrary initialization
* ----------------------------------------------------------------------

      if (crtlsf) then

#ifdef COUP_CCM
* obtain fractional land data if need to create surface dataset
*  o from oro field passed from ccm3 (as default) or
*  o from specified fractional land dataset (as specified in namelist input)

         if (flnddat .eq.
     &      'fractional land data obtained from ccm model') then
            do j = 1, lsmlat
               do i = 1, numlon(j)
                  if (nint(oroatm(i,j)) .eq. 1) then
                     fland(i,j) = 1.
                  else
                     fland(i,j) = 0.
                  endif
               end do
            end do
         else
            if (masterproc) then
               nflnd = getavu()
               call getfil (flnddat, locfn)
               call opnfil (locfn, nflnd, 'f')
               do j = 1, lsmlat
                  do i = 1, numlon(j)
                     read (nflnd,*) fland(i,j)
                  end do
               end do
               call relavu (nflnd) 
            endif
         endif

#else
* obtain fractional land from flux coupler ONLY if surface data 
* needs to be determined.
*  o arget(i,j,1) => grid land mask values (real)
*  o arget(i,j,2) => grid fractional land values

         do i=1,nibuff
            ibuff(i) = 0
         end do
         msg_id = 70
         call msg_recv_i(ibuff, nibuff, tid_cpl, msg_id)
         call msg_recv_r(arget, nsizr , tid_cpl, msg_id)
         if(masterproc) write(6,*) '(lsmini) recd d->l land fraction, msg_id= ',msg_id

         do j=1,lsmlat
            do i=1,numlon(j)
               mland(i,j) = int(arget(i,j,1))
               fland(i,j) =     arget(i,j,2)
               if (mland(i,j) .eq. 0) fland(i,j) = 0.
            end do
         end do
#endif

* create lsm surface data 

         if (masterproc) then
            call fsurdat(lsmgeo  , numlon  , fland   , surf2d, soic2d,
     &                   sand2d  , silt2d  , clay2d  , pctlak, pctwet,
     &                   drngbasn, drngarea, drngfrac)
         endif

      endif                    ! end of crtlsf block

* ----------------------------------------------------------------------
* set up vegetation type [veg] and weight [wt] arrays for [msub] subgrid 
* points on [lsmlon] x [lsmlat] grid. [wt] is used to build mapping indices
* in subroutine [lsmmap]. [veg] and [wt] are used to build subgrid vectors
* in subroutine [lsmtci].
* ----------------------------------------------------------------------
 
* m=1 to m=3 are "vegetated" subgrid points. fractional areas for
* these points [cover] pertain to "vegetated" area not to total grid area.
* so need to adjust them for fraction of grid covered by inland water.

      do m = 1, 3
!$omp parallel do default(shared) private(i,j)
         do j = 1, lsmlat
            do i = 1, numlon(j)
               wt(i,j,m) = cover(surf2d(i,j)+1,m) *
     &                     (100.-pctlak(i,j)-pctwet(i,j))/100.
               veg(i,j,m) = plant(surf2d(i,j)+1,m)
            end do
         end do
      end do

* m=4 is lake    subgrid point. set veg to no veg. set weight to lake fraction
* m=5 is wetland subgrid point. set veg to no veg. set weight to wetland fract

!$omp parallel do default(shared) private(i,j)
      do j = 1, lsmlat
         do i = 1, numlon(j)
            if (surf2d(i,j) .gt. 0) then            !land
               veg(i,j,4) = noveg                   
               wt(i,j,4)  = pctlak(i,j)/100.
               veg(i,j,5) = noveg                   
               wt(i,j,5)  = pctwet(i,j)/100.
            else                                    !ocean
               veg(i,j,4) = 0                   
               wt(i,j,4)  = 0.
               veg(i,j,5) = 0                   
               wt(i,j,5)  = 0.
            end if
         end do
      end do

* ----------------------------------------------------------------------
* determine lpt, kpt and numlv
* ----------------------------------------------------------------------

      k = 0
      l = 0
      do j = 1, lsmlat
         do i = 1, numlon(j)
            if (surf2d(i,j) .gt. 0) then           !land point
               l = l+1                             !land index
               do m = 1, msub
                  if (wt(i,j,m) .gt. 0.) then      !valid subgrid point
                     k = k+1                       !subgrid index
                  end if
               end do
            end if
         end do
      end do
      lpt = l
      kpt = k

!     numlv = max0(kpt/128,1)

      if ( masterproc ) then
        call getenv('NUMBER_CPUS_PER_MLP_PROCESS',evalue)
        read(evalue,*) numcpu             
      endif
#if defined(SPMD)
      call mp_bcst_int(lpt)
      call mp_bcst_int(kpt)
      call mp_bcst_int(numcpu)
#endif
      numlv = numcpu               

      if( masterproc ) write(*,*) 'lsmlat=',lsmlat
      if( masterproc ) write(*,*) 'lpt= ',lpt, ' kpt= ',kpt,' numlv= ', numlv

* ----------------------------------------------------------------------
* dynamically allocate memory for common block variables that
* depend on kpt and lpt
* ----------------------------------------------------------------------

      call alloctc
      call alloctv
      call allocht
      call allocbas

* ----------------------------------------------------------------------
* use [wt] to build mapping indices and weights: 
* [lsmlon] x [lsmlat] grid <-> 
* [lpt] vector of land points <-> 
* [kpt] vector of subgrid points. 
* ----------------------------------------------------------------------

      call lsmmap (wt)
      numks = endkptspmd(endlatspmd) - begkptspmd(beglatspmd) + 1

* echo lsm grid parameters

      minlon =  999999
      maxlon = -999999
      do j = 1, lsmlat
         minlon = min0(minlon, numlon(j))
         maxlon = max0(maxlon, numlon(j))
      end do

      if (masterproc) then
         write (6,*) 'grid characteristics:'
         write (6,*) '   latitude points = ',lsmlat
         write (6,*) '   maximun longitude points = ',maxlon
         write (6,*) '   minimun longitude points = ',minlon
         write (6,*) '   total grid land points   = ',lpt
         write (6,*) '   total grid land points incl subgrid points = ',
     &        kpt
      endif
* ----------------------------------------------------------------------
* initialize time invariant variables as subgrid vectors of length 
* [kpt] (for [lsmtcdyn] common block) and set up multiprocessing.
* ----------------------------------------------------------------------

      call lsmtci (veg, wt, beglatspmd)
      call basinini
 
* ----------------------------------------------------------------------
* if initial or branch run, write history file with time invariant 
* variables (from lsmtc common block)
* ----------------------------------------------------------------------

#if defined(LSMH_off)
      if(masterproc) write(6,*) 'LSM history file output is turned off'
#else
      if (nsrest.eq.0 .or. nsrest.eq.3) then
         ncbasin = .false.
         call histtc(beglatspmd)
         ncbasin = .true.
      endif
#endif

* ----------------------------------------------------------------------
* if continuation run or branch run, read restart files 
* ----------------------------------------------------------------------

      if (nsrest .gt. 0) call lsmres(nfrst)

* ----------------------------------------------------------------------
* initialize variables for history files
* ----------------------------------------------------------------------

#if defined(LSMH_off)
#else
      call histini (nsrest, fnameh)
#endif

* ----------------------------------------------------------------------
* determine calendar information for current time step
* ----------------------------------------------------------------------

! SJL
!     call calendr(nbstep  ,dtlsm  ,mdbase  ,msbase  ,mbdate  ,
!    &             mbsec   ,mdcur  ,mscur   ,mcdate  ,mcsec   ,
!    &             calday  )

      mcdate = xmcdate
      mcsec = xmcsec
      mscur = mcsec
      nhms = nhmsf( mcsec )

      call mcalday(mcdate, nhms, calday)
      mdcur = mcdate - mdbase     ! quick & dirty
                                  ! Does not give the correct value
                                  ! But it is not used if lsm history file
                                  ! is turned off
! SJL

      kyr = mcdate/10000
      kmo = mod(mcdate,10000)/100
      kda = mod(mcdate,100)

* for monthly averaged history tapes, the following variables
* must be initialized here if the initial run is not started from 
* the beginning of the month. if the run is started at the beginning 
* of the month, then the values are overwritten in routine histend.

      if (ninavg .eq. 'Q' .and. nsrest .eq. 0) then 
         nmon     = kmo           
         nyr      = kyr            
         mcdate_f = mcdate
         mcsec_f  = mcsec
         mdcur_f  = mdcur
         mscur_f  = mscur
      endif

* ----------------------------------------------------------------------
* if continuation or branch run: end of lsm initialization because 
* time-varying data in [lsmtvdyn] common block has been read in by lsm
* in subroutine lsmres. (also note that atmospheric model has albedos, 
* upward longwave flux, etc. in its own restart file).
* ----------------------------------------------------------------------

      if (nsrest .gt. 0) go to 1000
 
* ----------------------------------------------------------------------
* if initial run, initialize time-varying data in [lsmtvdyn] common block 
* ----------------------------------------------------------------------

      call lsmtvi (finidat, nfini, calday, nbstep, kmo, kda)

* ----------------------------------------------------------------------
* if initial run, return surface variables for atmospheric model
* ----------------------------------------------------------------------
      
#ifdef COUP_CCM
      call atmini (londim ,latdim ,beglatspmd,endlatspmd,lwupxy ,
     &             tsxy   ,asdirxy,asdifxy   ,aldirxy   ,aldifxy,
     &             snowxy, oroatm ) 
#else
      call atmini (londim ,latdim ,beglatspmd,endlatspmd,lwupxy ,
     $             tsxy   ,asdirxy,asdifxy   ,aldirxy   ,aldifxy,
     $             snowxy ) 
#endif

* ----------------------------------------------------------------------
* end lsm initialization
* ----------------------------------------------------------------------

1000  continue
      if (masterproc) then
         write (6,*) ('successful initialization of LSM')
         if (nsrest .eq. 0) then
            write (6,*) 'begin initial run at:'
         else
            write (6,*) 'begin continuation run at:'
         end if
         write (6,*) '   nstep= ',nbstep,' year= ',kyr,' month= ',kmo,
     &        ' day= ',kda,' seconds= ',mscur
         write (6,*)
         write(6,*)'**************************************************'
         write (6,*)
      end if                    !masterproc
      return
      end
 
