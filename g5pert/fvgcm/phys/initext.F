#include <params.h>
      subroutine initext()

! Initialize external models and/or boundary dataset information

#if defined( CLM2 )
! JDR Added clm2 modules
      use precision
      use shr_orb_mod
      use atm_lndMod
#endif
 
      implicit none

!------------------------------Parameters-------------------------------
#include <pmgrid.h>
#include <pagrid.h>
!------------------------------Commons----------------------------------
#include <comcon.h>
#include <comlun.h>
#include <comctl.h>
#include <commap.h>
#include <commss.h>
#include <comsrf.h>
#include <comsgh.h>
#include <comtim.h>
#include <comhst.h>

#if defined( CLM2 )
! JDR Added clm2 include file for orbital parameters
#include <comsol.h>
#endif

      common/icemask/icemask(plon,plat)

!--------------------------Local Variables------------------------------
      real icemask
      integer i,lat
      real coszrs(plond)          ! Cosine solar zenith angle
      integer numlon(plat)        ! Number of longitudes per latitude
      real tfrz
      data tfrz   /273.16   / !freezing point (kelvin)

! JDR Added clm2 variables
      logical log_print       ! Flag to print out log information or not

!-----------------------------------------------------------------------
 
      if (masterproc) call attachnc(ncid_sst, bndtvs  )
      if (masterproc) call attachnc(ncid_oz,  bndtvo  )
 
! Read in ozone boundary dataset
      call oznini

! Set flag values for fields. These things will be initialized by 
! the land model. Do not initialize TS over sea ice to flag value
! since its value is set by INIDAT which has already been called.
 
      if (nstep == 0) then
!       do lat=beglat,endlat
        do lat=begj,endj
          do i=1,plon
            asdir(i,lat) = 1.e36
            asdif(i,lat) = 1.e36
            aldir(i,lat) = 1.e36
            aldif(i,lat) = 1.e36
          end do
        end do
      end if
 
! Initialize land surface model variables: land albedos, land 
! surface temperature, land oro flag and land snowh.
! NOTE: On restart, oro, lwup, ts, albedos and snowh, come from the 
! atm *restart data*.  They are computed in LSMINI *only* at nstep=0.
! NOTE: The fractional land data needed by the LSM model
! orography flags obtained from the ORO field on the initial 
! dataset. The lsm may then change this array, for example if the 
! Ross ice shelf is exteneded. 
 
      if (nstep == 0) then
         do lat=1,plat
            numlon(lat) = plon
         end do
      end if

#if ( !defined AQUA )
        irt = 0
#if defined ( fvPSAS_ON )
        nhtfrq(1)  = -6           ! write frequency for history time samples
        ninavg(1)  = ' '
#else
        nhtfrq(1)  = -12           ! write frequency for history time samples
        ninavg(1)  = 'q'
#endif

#if defined( CLM2 )
! JDR Added calls for clm initialization 
! Setup the characteristics of the orbit
! (Based on the namelist parameters)
!
        if (masterproc) then
           log_print = .true.
        else
           log_print = .false.
        end if
        call shr_orb_params (iyear_AD, eccen , obliq , mvelp, obliqr, 
     &                   lambm0, mvelpp, log_print)
!
! Initialize land model. This involves initializing land 
! albedos, surface temperature, lwup and snowh.  NOTE: On restart, 
! lwup, ts, albedos and snowh, come from the atm restart data.  
!
!   if (.not. adiabatic .and. .not. ideal_phys .and. .not. aqua_planet) then
!        if (.not. adiabatic)then
           call atmlnd_ini()
!        endif

! JDR
#else
        call lsmini(plond    ,plat     ,begj     ,endj   ,
     $              numlon   ,oro      ,mdbase   ,msbase   ,
     $              mbdate   ,mbsec    ,nstep    ,dtime    ,
     $              nsrest   ,nestep   ,nhtfrq(1),mfilt(1) ,
     $              ninavg(1),irt      ,iradsw   ,lwup     , ts ,
     $              asdir    ,asdif    ,aldir    ,aldif    ,
     $              snowh ,   mcdate,   mcsec   )
#endif
#endif

! Ocean/Sea ice surface datasets: Determine sea surface temperature 
! from input dataset.
  
      call sstini
      if(nstep == 0) call sstint

      if (anncyc .and. nstep==0) then
 
!       do lat=beglat,endlat
        do lat=begj,endj
 
          if (interp) then  ! using interpolated initial data
            do i=1,plon
              if (nint(oro(i,lat)).ne.1) then !non-land
                if (tssub(i,1,lat).le.tfrz) then
                  snowh(i,lat) = 0.005
                  sicthk(i,lat) = 2.0
                  tssub(i,2,lat)=tssub(i,1,lat)-
     $                 .25*(tssub(i,1,lat)-tfrz)
                  tssub(i,3,lat)=tssub(i,1,lat)-
     $                 .50*(tssub(i,1,lat)-tfrz)
                  tssub(i,4,lat)=tssub(i,1,lat)-
     $                 .75*(tssub(i,1,lat)-tfrz)
                  oro(i,lat)=2.
                end if
              end if
            end do
          else   !using non interpolated initial data
 
! Determine differences between LSM surface and initial dataset.
 
            do i=1,plon
              if (nint(oro(i,lat)).ne.1) then ! LSM sfc data says non-land
                if (nint(icemask(i,lat)).eq.2) then ! initial data says sea ice
                  oro(i,lat) = 2.0
                end if
              else              ! LSM sfc data says land
                if (nint(icemask(i,lat)).eq.2) then ! initial data says sea ice
                  if(masterproc) write(6,*)'LSM sfc data says land, initial data ',
     $                  'says sea ice at i,lat=',i,lat,
     $                  ' LSM data will be used'
                else if(nint(icemask(i,lat)).eq.0) then ! init data says ocean
                  if(masterproc) write(6,*)'LSM sfc data says land, initial data ',
     $                  'says ocean i,lat=',i,lat,
     $                  ' LSM data will be used'
                end if
              end if
            end do
          end if
 
! set new sea ice concentrations
 
          call sstan(lat, oro(1,lat), ts(1,lat), tssub(1,1,lat),
     $               snowh(1,lat),    sicthk(1,lat))
 
! Compute longwave up over non-land
 
            do i=1,plon
              if (nint(oro(i,lat)).ne.1) then
                 lwup(i,lat) = stebol*(ts(i,lat)**4)
              end if
            end do
        end do
      end if
 
! Initialize non-land albedos at NSTEP = 0.  At NSTEP = 1 and 
! beyond, albedos will be computed for the *next* timestep to 
! accomodate coupling with a single interface.
 
      if (nstep == 0) then
!       do lat=beglat,endlat
        do lat=begj,endj
          call zenith  (calday  ,dodiavg ,clat(lat), coszrs  )
          call albocean(oro(1,lat), snowh(1,lat), coszrs  ,
     $                  asdir(1,lat), aldir(1,lat),
     $                  asdif(1,lat), aldif(1,lat))
        end do
      end if	

      return
      end
