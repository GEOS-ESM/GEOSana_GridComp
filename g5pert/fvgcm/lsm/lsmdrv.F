#include <preproc.h>
      subroutine lsmdrv (londim ,latdim  ,beglatspmd,endlatspmd,pcnst  ,
     &                   pgcmxy ,thgcmxy ,ugcmxy    ,vgcmxy    ,tgcmxy ,
     &                   qgcmxy ,flwdsxy ,zgcmxy    ,prcxy     ,prlxy  ,
     &                   solsxy ,sollxy  ,solsdxy   ,solldxy   ,doalb  ,
     &                   nstep  ,tsxy    ,shxy      ,lhxy      ,cfxy   ,
     &                   tauxxy ,tauyxy  ,trefxy    ,asdirxy   ,asdifxy,
     &                   aldirxy,aldifxy ,lwupxy    ,snowxy    , 
     &                   xmcdate,xmcsec  ,z0mxy     ,z0hxy     ,zpdxy  ,
     &                   wet1xy ,wetxy   ,tradxy   ,emsfcxy   )

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>
* ------------------------ code history ---------------------------
* source file:       lsmdrv.F
* purpose:           2-d model driver using [lsmlon] by [lsmlat] grid
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      
* reviewed:         
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* the land surface model works by gathering all the land points on a
* [lsmlon] by [lsmlat] grid into a vector of [lpt] land points. this is
* then expanded into a "big" vector of [kpt] subgrid points, allowing for
* up to [msub] subgrid points per land point. the [kpt] subgrid points
* are processed as [numlv] "little" vectors of [numkpt] points for [numlv] 
* calls to subroutine lsm. this subroutine:

*    o generates the calendar day (1.00 -> 365.99), month (1 -> 12),
*      and day (1 -> 31), which are used to calculate the surface 
*      albedos and leaf and stem areas for the next time step

*    o maps atmospheric fields from the [lsmlon] by [lsmlat] grid to 
*      subgrid vectors of length [kpt]

*    o calls the vector land surface model code [numlv] times in 
*      strips of [numkpt] points 

*    o maps fields from the subgrid vectors of length [kpt] to the 
*      [lsmlon] by [lsmlat] grid

*    o determines if end of history interval and writes history 
*      and restart files if appropriate

* -----------------------------------------------------------------

* ------------------- input variables -----------------------------
      integer londim           !atm number of longitudes
      integer latdim           !atm number of latitudes
      integer beglatspmd    !beg lat on proc, degenerates to 1 for nonspmd
      integer endlatspmd    !end lat on proc, degenerates to atmlat for nonspmd
      integer pcnst            !atm number of constituents (water is first)
      integer nstep            !atm time step index
      logical doalb            !true if surface albedo calculation time step

      integer xmcdate     
      integer xmcsec     

      real  pgcmxy(londim,beglatspmd:endlatspmd) !atm bottom level pressure (pa)
      real thgcmxy(londim,beglatspmd:endlatspmd) !atm btm lvl pot temp (kelvin)
      real  ugcmxy(londim,beglatspmd:endlatspmd) !atm btm level zonal wind (m/s)
      real  vgcmxy(londim,beglatspmd:endlatspmd) !atm btm lev merid wind (m/s)
      real  tgcmxy(londim,beglatspmd:endlatspmd) !atm bottom level temp (kelvin)
      real  qgcmxy(londim,beglatspmd:endlatspmd) !atm btm lvl spc humidity(kg/kg)
      real  zgcmxy(londim,beglatspmd:endlatspmd) !atm btm lev hght above srf (m)
      real   prcxy(londim,beglatspmd:endlatspmd) !conv precip rate (mm h2o/s)
      real   prlxy(londim,beglatspmd:endlatspmd) !large-scale prec rate(mm h2o/s)
      real flwdsxy(londim,beglatspmd:endlatspmd) !downward longwave rad onto surface (w/m**2)
      real solsxy (londim,beglatspmd:endlatspmd) !vis direct beam solar rad onto srf (w/m**2)
      real sollxy (londim,beglatspmd:endlatspmd) !nir direct beam solar rad onto srf (w/m**2)
      real solsdxy(londim,beglatspmd:endlatspmd) !vis diffuse solar rad onto srf (w/m**2)
      real solldxy(londim,beglatspmd:endlatspmd) !nir diffuse solar rad onto srf(w/m**2)

! ------------------- output variables ----------------------------
      real shxy(londim,beglatspmd:endlatspmd)   !sens heat flux (w/m**2) 
!                                               ![+ to atm]
      real lhxy(londim,beglatspmd:endlatspmd)   !latent heat flux (w/m**2)
!                                               ![+ to atm]
      real cfxy(londim,pcnst,beglatspmd:endlatspmd)!constituent fluxes 
!	                                        !(kg/m**2/s) [+ to atm]
      real tauxxy(londim,beglatspmd:endlatspmd) !zonal surf stress (kg/m/s**2)
      real tauyxy(londim,beglatspmd:endlatspmd) !merid surf stress (kg/m/s**2)
!
      real tsxy(londim,beglatspmd:endlatspmd)          !surface skin temperature (kelvin)
      real lwupxy(londim,beglatspmd:endlatspmd)        !emitted longwave radiation (w/m**2)
      real trefxy(londim,beglatspmd:endlatspmd)        !2 m height air temp (kelvin)
!
      real asdirxy(londim,beglatspmd:endlatspmd)       !albedo -- visible waveband, direct
      real asdifxy(londim,beglatspmd:endlatspmd)       !albedo -- visible waveband, diffuse
      real aldirxy(londim,beglatspmd:endlatspmd)       !albedo -- near infrared waveband, direct
      real aldifxy(londim,beglatspmd:endlatspmd)       !albedo -- near infrared waveband, diffuse
      real snowxy (londim,beglatspmd:endlatspmd)       !snow, liquid water equivalent (mm)
! JDC modify begin
      real z0mxy(londim,beglatspmd:endlatspmd)         !roughness length, momentum (m)
      real z0hxy(londim,beglatspmd:endlatspmd)         !roughness length, sensible heat (m)
      real zpdxy(londim,beglatspmd:endlatspmd)         !displacement height (m)
      real wet1xy(londim,beglatspmd:endlatspmd)        !top soil layer wetness (fraction)
      real wetxy(londim,beglatspmd:endlatspmd)         !root zone soil wetness (fraction)
      real tradxy(londim,beglatspmd:endlatspmd)        !surface brightness temperature (k)
      real emsfcxy(londim,beglatspmd:endlatspmd)       !bulk surface emissivity (fraction)
! JDC modify end

! ------------------- common block variables ----------------------
#include <dynmem.h>
#include <lsmtc.h>
#include <lsmtcdyn.h>
#include <lsmtvdyn.h>
#include <phycon.h>
#include <lsmctl.h>
! -----------------------------------------------------------------

! ---------------------- local variables --------------------------
      integer i,j,k,l,m        !loop/array indices
      real wt                  !subgrid weight

! calendar info
! -----------------------------------------------------------------
      integer mdcur            !current day (0, ...)
      integer mscur            !current seconds of current day (0, ..., 86400)
      integer mcdate           !current date (yymmdd format) (e.g., 030131)
      integer mcsec            !current seconds of current date (0, ..., 86400)
      integer kyr              !year (0, ...)
      integer kmo              !month (1, ..., 12)
      integer kda              !day of month (1, ..., 31)
      real calday              !calendar day at greenwich (1.00, ..., 365.99)

! constant atmospheric co2 and o2
! -----------------------------------------------------------------
      real po2                 !partial pressure  o2 (mol/mol)
      real pco2                !partial pressure co2 (mol/mol)
      data po2,pco2 /0.209,355.e-06/
      save po2,pco2

! subgrid vector version of atm model input
! -----------------------------------------------------------------

      real, allocatable :: pgcm(:)           !atmospheric pressure (pa)
      real, allocatable :: tgcm(:)           !atmospheric temperature (kelvin)
      real, allocatable :: qgcm(:)           !atmospheric specific humidity (kg/kg)
      real, allocatable ::  firgcm(:)        !downward infrared (longwave) radiation (w/m**2)
      real, allocatable ::  solad(:,:)       !direct beam solar radiation (w/m**2) (vis, nir)
      real, allocatable ::  solai(:,:)       !diffuse solar radiation (w/m**2) (vis, nir)
      real, allocatable ::  qprecc(:)        !convective precipitation (mm h2o /s)
      real, allocatable ::  qprecl(:)        !large-scale precipitation (mm h2o /s)
      real, allocatable ::  ugcm(:)          !atmospheric wind speed in east direction (m/s)
      real, allocatable ::  vgcm(:)          !atmospheric wind speed in north direction (m/s)
      real, allocatable ::  hgcm(:)          !atmospheric reference height (m)
      real, allocatable ::  thgcm(:)         !atmospheric potential temperature (kelvin)
      real, allocatable ::  egcm(:)          !atmospheric vapor pressure (pa)
      real, allocatable ::  rhogcm(:)        !density (kg/m**3)
      real, allocatable ::  co2gcm(:)        !atmospheric co2 concentration (pa)
      real, allocatable ::  o2gcm(:)         !atmospheric o2 concentration (pa)

! subgrid vector version of lsm output
! -----------------------------------------------------------------

      real, allocatable ::  taux(:)           !wind stress: e-w (kg/m/s**2)
      real, allocatable ::  tauy(:)           !wind stress: n-s (kg/m/s**2)
      real, allocatable ::  fire(:)           !emitted infrared (longwave) radiation (w/m**2) 
      real, allocatable ::  fsh(:)            !sensible heat flux (w/m**2) [+ to atm]
      real, allocatable ::  fcev(:)           !veg evaporation heat flux (w/m**2) [+ to atm]
      real, allocatable ::  fgev(:)           !ground evaporation heat flux (w/m**2) [+ to atm]
      real, allocatable ::  fctr(:)           !veg transpiration heat flux (w/m**2) [+ to atm]
      real, allocatable ::  trad(:)           !radiative temperature (kelvin)
      real, allocatable ::  tsa(:)            !2 m height surface air temperature (kelvin)
      real, allocatable ::  qsoil(:)          !soil evaporation (mm h2o/s) (+ = to atm)
      real, allocatable ::  qvege(:)          !vegetation evaporation (mm h2o/s) (+ = to atm)
      real, allocatable ::  qvegt(:)          !vegetation transpiration (mm h2o/s) (+ = to atm)
      real, allocatable ::  fpsn(:)           !photosynthesis (umol co2 /m**2 /s)
      real, allocatable ::  frm(:)            !total maintenance respiration (umol co2 /m**2/s)
      real, allocatable ::  frg(:)            !growth respiration (umol co2 /m**2 /s)
      real, allocatable ::  fmicr(:)          !microbial respiration (umol co2 /m**2 /s)
      real, allocatable ::  fco2(:)           !net co2 flux (umol co2 /m**2 /s) [+ = to atm]
      real, allocatable ::  qover(:)          !surface runoff (mm h2o /s) 
      real, allocatable ::  qdrai(:)          !sub-surface runoff (mm h2o /s) 
! JDC modify begin
      real, allocatable ::  z0m(:)            !roughness length, momentum (m)
      real, allocatable ::  z0h(:)            !roughness length, sensible heat (m)
      real, allocatable ::  zpd(:)            !displacement height (m)
      real, allocatable ::  wet1(:)           !top soil layer wetness (fraction)
      real, allocatable ::  wet(:)            !root zone soil wetness (fraction)
      real, allocatable ::  tskin4(:)         !4th power of skin temperature
      real, allocatable ::  trad4(:)          !4th power of radiative temperature
! JDC modify end 

      integer nhms, nhmsf, n
      nhmsf(n)  = n/3600*10000 + mod(n,3600 )/ 60*100 + mod(n, 60)
! -----------------------------------------------------------------

      allocate (pgcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (tgcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (qgcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (firgcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate 
     &    (solad(mband,begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate 
     &    (solai(mband,begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (qprecc(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (qprecl(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (ugcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (vgcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (hgcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (thgcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (egcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (rhogcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (co2gcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (o2gcm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))

      allocate (taux(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (tauy(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (fire(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (fsh(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (fcev(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (fgev(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (fctr(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (trad(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (tsa(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (qsoil(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (qvege(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (qvegt(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (fpsn(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (frm(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (frg(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (fmicr(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (fco2(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (qover(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (qdrai(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
! JDC modify begin
      allocate (z0m(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (z0h(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (zpd(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (wet1(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (wet(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (tskin4(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
      allocate (trad4(begkptspmd(beglatspmd):endkptspmd(endlatspmd)))
! JDC modify end

* ----------------------------------------------------------------------
* determine if end of history interval
* ----------------------------------------------------------------------

      call histend (nstep, xmcdate, xmcsec)          ! SJL

* ----------------------------------------------------------------------
* calendar information for next time step
*    o calday = calendar day (1.00 -> 365.99) for cosine solar zenith angle
*               calday is based on greenwich time
*    o kmo    = month (1 -> 12) for leaf area index and stem area index
*    o kda    = day (1 -> 31)   for leaf area index and stem area index
* ----------------------------------------------------------------------

! SJL
!     call calendr(nstep+1,dtlsm  ,mdbase  ,msbase  ,mbdate  ,
!    &             mbsec  ,mdcur  ,mscur   ,mcdate  ,mcsec   ,
!    &             calday )

      mcdate = xmcdate
      mcsec = xmcsec
      nhms = nhmsf( mcsec )

      call tick ( mcdate, nhms, nint(dtlsm) )
      call mcalday(mcdate, nhms, calday)
! SJL

      kyr = mcdate/10000
      kmo = mod(mcdate,10000)/100
      kda = mod(mcdate,100)

* ----------------------------------------------------------------------
* map atmospheric fields to force lsm: [lsmlon] x [lsmlat] grid ->
* [lpt] vector of land points -> [kpt] vector of subgrid points
* ----------------------------------------------------------------------

      do k = begkptspmd(beglatspmd),endkptspmd(endlatspmd) 
         l = klnd(k)                            !land point index
         i = ixy(l)                             !longitude index
         j = jxy(l)                             !latitude index 
         tgcm(k)    = tgcmxy(i,j)      
         hgcm(k)    = zgcmxy(i,j)     
         ugcm(k)    = ugcmxy(i,j)    
         vgcm(k)    = vgcmxy(i,j)   
         qgcm(k)    = qgcmxy(i,j)  
         pgcm(k)    = pgcmxy(i,j)      
         thgcm(k)   = thgcmxy(i,j)     
         qprecc(k)  = prcxy(i,j)   
         qprecl(k)  = prlxy(i,j)  
         solad(1,k) = solsxy(i,j)     
         solad(2,k) = sollxy(i,j)    
         solai(1,k) = solsdxy(i,j)      
         solai(2,k) = solldxy(i,j)     
         firgcm(k)  = flwdsxy(i,j)     
      end do

* potential temperature, vapor pressure, air density at atm reference 
* height. atmospheric co2 and o2 are currently constants

!$omp  parallel do default(shared) private(k)
      do k = begkptspmd(beglatspmd), endkptspmd(endlatspmd)
         egcm(k)   = qgcm(k)*pgcm(k) / (0.622+0.378*qgcm(k))
         rhogcm(k) = (pgcm(k)-0.378*egcm(k)) / (rair*tgcm(k))
         co2gcm(k) = pco2*pgcm(k)
         o2gcm(k)  = po2*pgcm(k) 
      end do

* ----------------------------------------------------------------------
* process the "big" vector of [kpt] points as [numlv] "little" 
* vectors of [numkpt] points. [begkpt] is the starting location of the 
* [numkpt] points for the "little" vector in the "big" [kpt] vector.
* multitask [numlv] processes
* ----------------------------------------------------------------------

!$omp parallel do default(shared) private(i,k)

      do i = 1, numlv
         k = begkpt(i)
      if( k /= 0 ) then  
         call lsm (nstep      ,numkpt(i)  ,dtlsm      ,dtsoi      ,
     &  doalb     ,pgcm (k)   ,tgcm(k)    ,qgcm(k)    ,firgcm(k)  ,
     &  lati(k)   ,long(k)    ,solad(1,k) ,solai(1,k) ,co2gcm(k)  ,
     &  o2gcm(k)  ,qprecc(k)  ,qprecl(k)  ,thgcm(k)   ,rhogcm(k)  ,
     &  egcm(k)   ,igs(k)     ,hsno(k)    ,fsno(k)    ,fwet(k)    ,
     &  ugcm(k)   ,vgcm(k)    ,hgcm(k)    ,ivt(k)     ,ist(k)     ,
     &  isc(k)    ,dzsoi(1,k) ,zsoi(1,k)  ,tsoi(1,k)  ,foln(k)    ,
     &  tg(k)     ,tv(k)      ,moz(k)     ,h2osoi(1,k),h2ocan(k)  ,
     &  h2osno(k) ,soot(k)    ,eah(k)     ,fsh(k)     ,frg(k)     ,
     &  fire(k)   ,albd(1,k)  ,albi(1,k)  ,albgrd(1,k),albgri(1,k),
     &  fabd(1,k) ,fabi(1,k)  ,ftdd(1,k)  ,ftid(1,k)  ,ftii(1,k)  ,
     &  fsun(k)   ,htop(k)    ,elai(k)    ,esai(k)    ,tlai(k)    ,
     &  fcev(k)   ,fctr(k)    ,fgev(k)    ,qsoil(k)   ,qvege(k)   ,
     &  taux(k)   ,tauy(k)    ,trad(k)    ,tsa(k)     ,stemb(k)   ,
     &  fpsn(k)   ,frm(k)     ,fco2(k)    ,tsai(k)    ,rootb(k)   ,
     &  fmicr(k)  ,root(1,k)  ,watsat(k)  ,hksat(k)   ,soilc(k)   ,
     &  smpsat(k) ,bch(k)     ,watdry(k)  ,watopt(k)  ,csol(k)    ,
     &  tksol(k)  ,tkdry(k)   ,qvegt(k)   ,kmo        ,kda        ,
     &  calday    ,hydro      ,pergro     ,conchk     ,qover(k)   ,
     &  qdrai(k)  ,i          ,k          ,
     &  z0m(k)    ,z0h(k)     ,zpd(k)     ,wet1(k)    ,wet(k)     ,
     &  tskin4(k) ,trad4(k)   )
      endif
      end do

* ----------------------------------------------------------------------
* update basin runoff info
* ----------------------------------------------------------------------

      call basindrv(beglatspmd, endlatspmd, qover, qdrai)

* ----------------------------------------------------------------------
* accumulate counters. write lsm history and restart files
* ----------------------------------------------------------------------

      call histhan (nstep, nestep, beglatspmd, endlatspmd)

* ----------------------------------------------------------------------
* return required surface fields to atmospheric model:
* [kpt] vector of subgrid points -> [lpt] vector of land points -> 
* [lsmlon] x [lsmlat] grid. this mapping is for land points only.
* non-land points are undefined.
* ----------------------------------------------------------------------

* need to initialize 2-d fields to zero for subgrid averaging, but 
* only for land points because other surface models may have already 
* calculated fields for non-land points.
 
      do j = beglatspmd, endlatspmd
         do i = 1, numlon(j)
            if (surf2d(i,j) .gt. 0) then 
               shxy(i,j)    = 0.
               lhxy(i,j)    = 0.
               cfxy(i,1,j)  = 0.
               tauxxy(i,j)  = 0.
               tauyxy(i,j)  = 0.
               tsxy(i,j)    = 0.
               trefxy(i,j)  = 0.
               lwupxy(i,j)  = 0.
               asdirxy(i,j) = 0.
               asdifxy(i,j) = 0.
               aldirxy(i,j) = 0.
               aldifxy(i,j) = 0.
               snowxy(i,j)  = 0.
! JDC modify begin
               z0mxy(i,j)   = 0.
               z0hxy(i,j)   = 0.
               zpdxy(i,j)   = 0.
               wet1xy(i,j)  = 0.
               wetxy(i,j)   = 0.
               tradxy(i,j) = 0.
               emsfcxy(i,j) = 0.
! JDC modify end
            end if
         end do
      end do

      do m = 2, pcnst
         do j = beglatspmd, endlatspmd
            do i = 1, numlon(j)
               if (surf2d(i,j) .gt. 0) cfxy(i,m,j) = 0.
            end do
         end do
      end do

* map fields from subgrid vector, with length [kpt], to 
* [lsmlon] x [lsmlat] surface, weighting by subgrid fraction. 
* for each land point on [lsmlon] x [lsmlat] surface, process 
* [msub] subgrid points. if the subgrid point is not valid, a 
* dummy index to the subgrid vector is used and the weight is zero.

      do m = 1, msub              !subgrid points for each land point
!! fix this
!$omp parallel do default(shared) private(i,j,k,l,wt)

         do l = beglptspmd(beglatspmd),endlptspmd(endlatspmd)
            i = ixy(l)            !longitude index 
            j = jxy(l)            !latitude index 
            k = kvec(l,m)         !lsm subgrid vector index
            wt = wsg2g(l,m)       !subgrid weights 
c
c The following if test is so SPMD code will not use uninitialized stack
c memory values for arrays like taux.
c
            if (wt.ne.0.) then

            shxy(i,j)    = shxy(i,j)   + fsh(k)*wt
            lhxy(i,j)    = lhxy(i,j)   + (fgev(k)+fcev(k)+fctr(k))*wt
            cfxy(i,1,j)  = cfxy(i,1,j) + (qsoil(k)+qvege(k)+qvegt(k))*wt
            tauxxy(i,j)  = tauxxy(i,j) + taux(k)*wt
            tauyxy(i,j)  = tauyxy(i,j) + tauy(k)*wt
            tradxy(i,j)  = tradxy(i,j) + trad(k)*wt
            trefxy(i,j)  = trefxy(i,j) + tsa(k)*wt
            lwupxy(i,j)  = lwupxy(i,j) + fire(k)*wt
            asdirxy(i,j) = asdirxy(i,j)+ albd(1,k)*wt
            asdifxy(i,j) = asdifxy(i,j)+ albi(1,k)*wt
            aldirxy(i,j) = aldirxy(i,j)+ albd(2,k)*wt
            aldifxy(i,j) = aldifxy(i,j)+ albi(2,k)*wt
            snowxy(i,j)  = snowxy(i,j) + h2osno(k)*wt
! JDC modify begin
            z0mxy(i,j)   = z0mxy(i,j)  + z0m(k)*wt
            z0hxy(i,j)   = z0hxy(i,j)  + z0h(k)*wt
            zpdxy(i,j)   = zpdxy(i,j)  + zpd(k)*wt
            wet1xy(i,j)  = wet1xy(i,j) + wet1(k)*wt
            wetxy(i,j)   = wetxy(i,j)  + wet(k)*wt
            tsxy(i,j)    = tsxy(i,j)   + tskin4(k)*wt
            emsfcxy(i,j) = emsfcxy(i,j) + trad4(k)*wt
! JDC modify end

* set cfxy = co2 fluxes for pcnst > 1, multiplying by 4.4e-08
* to convert from umol co2/m**2/s to kg co2/m**2/s

c           cfxy(i,2,j) = cfxy(i,2,j) + fco2(k)*wt * 4.4e-08
c           cfxy(i,3,j) = cfxy(i,3,j) + fpsn(k)*wt * 4.4e-08
c           cfxy(i,4,j) = cfxy(i,4,j) + (frm(k)+frg(k))*wt * 4.4e-08
c           cfxy(i,5,j) = cfxy(i,5,j) + fmicr(k)*wt * 4.4e-08
            end if
         end do
      end do
!JDC modify begin
!$omp parallel do default(shared) private(i,j,l)

       do l = beglptspmd(beglatspmd),endlptspmd(endlatspmd)
          i = ixy(l)            !longitude index
          j = jxy(l)            !latitude index
          emsfcxy(i,j) = emsfcxy(i,j) / tsxy(i,j)
          tsxy(i,j) = tsxy(i,j) ** 0.25
       end do
!JDC modify end

!
! The following added by SJL 
!
      deallocate (pgcm)
      deallocate (tgcm)
      deallocate (qgcm)
      deallocate (firgcm)
      deallocate (solad)
      deallocate (solai) 
      deallocate (qprecc)
      deallocate (qprecl)
      deallocate (ugcm)
      deallocate (vgcm)
      deallocate (hgcm)
      deallocate (thgcm)
      deallocate (egcm)
      deallocate (rhogcm)
      deallocate (co2gcm)
      deallocate (o2gcm)
      deallocate (taux)
      deallocate (tauy)
      deallocate (fire)
      deallocate (fsh)
      deallocate (fcev)
      deallocate (fgev)
      deallocate (fctr)
      deallocate (trad)
      deallocate (tsa)
      deallocate (qsoil)
      deallocate (qvege)
      deallocate (qvegt)
      deallocate (fpsn)
      deallocate (frm)
      deallocate (frg)
      deallocate (fmicr)
      deallocate (fco2)
      deallocate (qover)
      deallocate (qdrai)
! JDC modify begin
      deallocate (z0m)
      deallocate (z0h)
      deallocate (zpd)
      deallocate (wet1)
      deallocate (wet)
      deallocate (tskin4)
      deallocate (trad4)
! JDC modify end

      return
      end
