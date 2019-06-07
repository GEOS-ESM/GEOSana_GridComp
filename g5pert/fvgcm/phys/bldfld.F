#include <params.h>
      subroutine bldfld
C-----------------------------------------------------------------------
C
C Build Master Field List of all possible fields in a history file.
C Below are listed the field name and description of all the fields
C in the Master Field List:
C
C 'PHIS    ': Surface geopotential
C 'PS      ': Surface pressure
C 'T       ': Temperature
C 'U       ': U component (east-west) of horizontal velocity
C 'V       ': V component (north-south) of horizontal velocity
C 'ETADOT  ': Vertical velocity (vertical "eta" velocity)
C tracnam   : Constituent fields (1 = moisture)
C hadvnam   : Horizontal advection of constituent fields
C vadvnam   : Vertical advection of constituent fields
C fixcnam   : Conservation fixer for constituent fields
C tottnam   : SLT constituent tendency (advection + fixer)
C vdiffnam  : Vertical diffusion of constituents
C dcconnam  : Convective adjustment of constituents
C tendnam   : Total constituent tendency
C srcnam    : Constituent sources/sinks
C sflxnam   : Constituent surface fluxes
C 'DUH     ': U horizontal diffusive heating
C 'DVH     ': V horizontal diffusive heating
C 'DTH     ': T horizontal diffusive heating
C 'ORO     ': Land, ocean, sea ice mask
C 'SNOWH   ': Water equivalent snow depth
C 'SICTHK  ': Sea-ice thickness
C 'PRECL   ': Large-scale (stable) precipitation rate
C 'PRECC   ': Convective precipitation rate
C 'PRECSL  ': Large-scale (stable) snow rate (water equivalent)
C 'PRECSC  ': Convective snow rate (water equivalent)
C 'SHFLX   ': Surface sensible heat flux
C 'LHFLX   ': Surface latent heat flux
C 'QFLX    ': Surface water flux
C 'PBLH    ': PBL height
C 'USTAR   ': Surface friction velocity
C 'TREFHT  ': Reference height temperature
C 'CGH     ': Counter-gradient term for heat in PBL
C 'CGQ     ': Counter-gradient term for moisture in PBL
C 'CGS     ': Counter-gradient coefficient on surface kinematic fluxes
C 'TPERT   ': Perturbation temperature (eddies in PBL)
C 'QPERT   ': Perturbation specific humidity (eddies in PBL)
C 'KVH     ': Vertical diffusion diffusivities (heat/moisture)
C 'KVM     ': Vertical diffusion diffusivities (momentum)
C 'DUV     ': U vertical diffusion
C 'DVV     ': V vertical diffusion
C 'DTV     ': T vertical diffusion
C 'FSNS    ': Net solar flux at surface
C 'FLNS    ': Net longwave flux at surface
C 'FLNT    ': Net longwave flux at top
C 'FSNT    ': Net solar flux at top
C 'CLOUD   ': Cloud fraction
C 'SETLWP  ': Specified liquid water path lengths
C 'CLDLWP  ': Cloud liq water path length (accounting for cloud fraction)
C 'EFFCLD  ': "Effective" cloud fraction
C 'FLNTC   ': Clearsky net longwave flux at top
C 'FSNTC   ': Clearsky net solar flux at top
C 'FLNSC   ': Clearsky net longwave flux at surface
C 'FSNSC   ': Clearsky net solar flux at surface
C 'OMEGA   ': Vertical pressure velocity
C 'DQP     ': Specific humidity tendency due to precipitation
C 'TAUX    ': X-component (east-west) of surface stress
C 'TAUY    ': Y-component (north-south) of surface stress
C 'SRFRAD  ': Net radiative flux at surface
C 'QRS     ': Solar heating rate
C 'QRL     ': Longwave heating rate
C 'CLDTOT  ': Vertically-integrated  random overlap, total cloud amount
C 'CLDLOW  ': Vertically-integrated, random overlap, low cloud amount
C 'CLDMED  ': Vertically-integ., random overlap, mid-level cloud amount
C 'CLDHGH  ': Vertically-integrated, random overlap, high cloud amount
C 'TOTLWP  ': Vertically-integrated actual liquid water path length
C 'LWSH    ': Liquid water scale height
C tsnam     : Subsoil temperatures - top level is surface temperature
C 'TS      ': Surface temp.  In uncoupled mode identical to tsnam(1).
C             In coupled mode, part land, part ocean points will differ.
C 'SOLIN   ': Solar insolation
C 'UTEND   ': U tendency
C 'VTEND   ': V tendency
C 'TTEND   ': T tendency
C 'LPSTEN  ': Surface pressure tendency
C 'UTGW    ': U tendency - gravity wave drag
C 'VTGW    ': V tendency - gravity wave drag
C 'TAUGWX  ': East-west gravity wave drag surface stress
C 'TAUGWY  ': North-south gravity wave drag surface stress
C 'DTCOND  ': T tendency - convective adjustment
C 'CMFDT   ': T tendency - hack moist convection
C 'CMFDQ   ': Q tendency - hack moist convection
C 'ZMDT    ': T tendency - zhang moist convection
C 'ZMDQ    ': Q tendency - zhang moist convection
C 'CMFDQR  ': Q tendency - moist convective rainout
C 'CMFMC   ': Moist convection mass flux
C 'CMFSL   ': Moist convection liquid water static energy flux
C 'CMFLQ   ': Moist convection total water flux
C 'CNVCLD  ': Random overlap total convective cloud amount
C 'FSDS    ': Flus Shortwave Downwelling Surface
C 'VT      ': Meridional heat transport
C 'VZ      ': Meridional transport
C 'VQ      ': Meridional water transport
C 'VVPUU   ': Kinetic Energy
C 'RELHUM  ': Relative humidity
C 'Z3      ': Geopoten. Ht. (above sea level) -- new hydrostatic matrix
C 'MQ      ': water mass
C 'PSL     ': Sea level pressure
C 'OMEGAUP ': Average of Omega when omega is < 0 (up)
C 'NUMOMGUP': Percentage of time omega is < 0 (up)
C 'CLOUDUP ': Average of Cloud when omega is < 0 (up)
C 'DPSLON  ': Longitude Pressure Gradient
C 'DPSLAT  ': Latitude Pressure Gradient
C 'O3VMR   ': Ozone volume mixing ratio
C 'T10M    ': 10 meter temperature
C 'Q10M    ': 10 meter specified humidity
C 'U10M    ': 10 meter U wind
C 'V10M    ': 10 meter V wind
C 'T2M     ': 2 meter temperature
C 'Q2M     ': 2 meter specified humidity
C 'U2M     ': 2 meter U wind
C 'V2M     ': 2 meter V wind
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath, April 1992
C Standardized:      L. Bath, June 1992
C                    T. Acker, March 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C------------------------------Commons----------------------------------
#include <comhst.h>
C-----------------------------------------------------------------------
#include <comtrcnm.h>
C---------------------------Local variables-----------------------------
C
      integer nf         ! Field counter
      integer m,j        ! Index
      integer nsingl     ! Flag for single-level field
      integer nmulth     ! Flag for multilevel field on interfaces
      integer nmultf     ! Flag for multilevel field on Model levels
      integer active     ! Flag to write field to history file by default
      integer inactive   ! Flag not to write field to history file
      integer flag       ! Active/inactive flag for mon-ave hist file 
C
      logical lprint     ! Flag for print of full table at end
C
      character*8 fieldns(2,pflds) ! Sorted list of field names
C
      integer isloccc(0:pichsum)   ! Current location in list
      integer islocdum(0:pichsum)  ! Local work array for islocc
      integer nfldscdu(0:pichsum)  ! Local work array for nfldsc
      integer ifldss(3,pflds)      ! Sorted list of "iflds"
      integer ipos                 ! Ichar sum
      integer i                    ! Index
C
C-----------------------------------------------------------------------
C
C Set lprint true to get a complete printout of the Master Field List
C
      lprint = .false.
C
C Give flags integer values
C
      nsingl = 0
      nmulth = 1
      nmultf = 2
      active = 1
      inactive = 0
      nf = 0

C
C Call subroutine addfld to add a field to the Master Field List.
C
      call addfld('PHIS    ','M2/S2   ',nf,nsingl,active)
      call addfld('PS      ','PA      ',nf,nsingl,active)
      call addfld('T       ','K       ',nf,nmultf,active)
      call addfld('U       ','M/S     ',nf,nmultf,active)
      call addfld('V       ','M/S     ',nf,nmultf,active)
      call addfld('ETADOT  ','1/S     ',nf,nmulth,active)
      call addfld('SGH     ','        ',nf,nsingl,inactive)
C
C Constituent tracers
C
      do m=1,pcnst
        call addfld(tracnam(m),'KG/KG   ',nf,nmultf,active)
        call addfld(hadvnam(m),'KG/KGS  ',nf,nmultf,inactive)
        call addfld(vadvnam(m),'KG/KGS  ',nf,nmultf,inactive)
        call addfld(fixcnam(m),'KG/KGS  ',nf,nmultf,inactive)
        call addfld(tendnam(m),'KG/KGS  ',nf,nmultf,inactive)
        if (m.eq.1) then
           call addfld(tottnam(m),'KG/KGS  ',nf,nmultf,active)
           call addfld(vdiffnam(m),'KG/KGS  ',nf,nmultf,active)
           call addfld(dcconnam(m),'KG/KGS  ',nf,nmultf,active)
        else
           call addfld(tottnam(m),'KG/KGS  ',nf,nmultf,inactive)
           call addfld(vdiffnam(m),'KG/KGS  ',nf,nmultf,inactive)
           call addfld(dcconnam(m),'KG/KGS  ',nf,nmultf,inactive)
        endif
      end do
      do m=2,pcnst
        call addfld(srcnam(m),'KG/KGS  ',nf,nmultf,inactive)
        call addfld(sflxnam(m),'KG/M2/S ',nf,nsingl,inactive)
      end do
c
c Non-advected species
c
      do m=pcnst+1,pcnst+pnats
        call addfld(tracnam(m),'UNKNOWN ',nf,nmultf,active)
      end do
C
      call addfld('DUH     ','K/S     ',nf,nmultf,inactive)
      call addfld('DVH     ','K/S     ',nf,nmultf,inactive)
      call addfld('DTH     ','K/S     ',nf,nmultf,active)
      call addfld('ORO     ','FLAG    ',nf,nsingl,active)
      call addfld('SNOWH   ','M       ',nf,nsingl,active)
#if ( defined COUP_SOM )
      call addfld('SICTHK  ','M       ',nf,nsingl,active)
#endif
      call addfld('PRECL   ','M/S     ',nf,nsingl,active)
      call addfld('PRECC   ','M/S     ',nf,nsingl,active)
      call addfld('PRECSL  ','M/S     ',nf,nsingl,inactive)
      call addfld('PRECSC  ','M/S     ',nf,nsingl,inactive)
      call addfld('SHFLX   ','W/M2    ',nf,nsingl,active)
      call addfld('LHFLX   ','W/M2    ',nf,nsingl,active)
      call addfld('QFLX    ','KG/M2/S ',nf,nsingl,active)
      call addfld('PBLH    ','M       ',nf,nsingl,active)
      call addfld('USTAR   ','M/S     ',nf,nsingl,active)
c      call addfld('TREFHT  ','K       ',nf,nsingl,active)
      call addfld('CGH     ','K/M     ',nf,nmulth,inactive)
      call addfld('CGQ     ','1/M     ',nf,nmulth,inactive)
      call addfld('CGS     ','S/M2    ',nf,nmulth,active)
      call addfld('TPERT   ','K       ',nf,nsingl,active)
      call addfld('QPERT   ','KG/KG   ',nf,nsingl,active)
      call addfld('KVH     ','M2/S    ',nf,nmulth,active)
      call addfld('KVM     ','M2/S    ',nf,nmulth,inactive)
      call addfld('DUV     ','M/S2    ',nf,nmultf,inactive)
      call addfld('DVV     ','M/S2    ',nf,nmultf,inactive)
      call addfld('DTV     ','K/S     ',nf,nmultf,active)
      call addfld('FSNS    ','W/M2    ',nf,nsingl,active)
      call addfld('FLNS    ','W/M2    ',nf,nsingl,active)
      call addfld('FLNT    ','W/M2    ',nf,nsingl,active)
      call addfld('FSNT    ','W/M2    ',nf,nsingl,active)
      call addfld('CLOUD   ','FRACTION',nf,nmultf,active)
      call addfld('SETLWP  ','G/M2    ',nf,nmultf,inactive)
      call addfld('CLDLWP  ','G/M2    ',nf,nmultf,inactive)
      call addfld('EFFCLD  ','FRACTION',nf,nmultf,active)
      call addfld('FLNTC   ','W/M2    ',nf,nsingl,active)
      call addfld('FSNTC   ','W/M2    ',nf,nsingl,active)
      call addfld('FLNSC   ','W/M2    ',nf,nsingl,active)
      call addfld('FSNSC   ','W/M2    ',nf,nsingl,active)
      call addfld('OMEGA   ','PA/S    ',nf,nmultf,active)
      call addfld('DQP     ','KG/KGS  ',nf,nmultf,inactive)
      call addfld('TAUX    ','N/M2    ',nf,nsingl,active)
      call addfld('TAUY    ','N/M2    ',nf,nsingl,active)
      call addfld('SRFRAD  ','W/M2    ',nf,nsingl,active)
      call addfld('QRS     ','K/S     ',nf,nmultf,active)
      call addfld('QRL     ','K/S     ',nf,nmultf,active)
      call addfld('CLDTOT  ','FRACTION',nf,nsingl,active)
      call addfld('CLDLOW  ','FRACTION',nf,nsingl,active)
      call addfld('CLDMED  ','FRACTION',nf,nsingl,active)
      call addfld('CLDHGH  ','FRACTION',nf,nsingl,active)
      call addfld('TOTLWP  ','FRACTION',nf,nsingl,inactive)
      call addfld('LWSH    ','M       ',nf,nsingl,inactive)
c
c Put TS1 on history file since it will be different than TS in a coupled
c run
c
      do m=1,plevmx
        call addfld(tsnam(m),'K       ',nf,nsingl,active)
      end do
      call addfld('TS      ','K       ',nf,nsingl,active)
      call addfld('SOLIN   ','W/M2    ',nf,nsingl,active)
      call addfld('UTEND   ','M/S2    ',nf,nmultf,inactive)
      call addfld('VTEND   ','M/S2    ',nf,nmultf,inactive)
      call addfld('TTEND   ','K/S     ',nf,nmultf,inactive)
      call addfld('LPSTEN  ','PA/S    ',nf,nsingl,inactive)
      call addfld('UTGW    ','M/S2    ',nf,nmultf,active)
      call addfld('VTGW    ','M/S2    ',nf,nmultf,active)
      call addfld('TAUGWX  ','N/M2    ',nf,nsingl,active)
      call addfld('TAUGWY  ','N/M2    ',nf,nsingl,active)
      call addfld('DTCOND  ','K/S     ',nf,nmultf,active)
      call addfld('CMFDT   ','K/S     ',nf,nmultf,active)
      call addfld('CMFDQ   ','KG/KGS  ',nf,nmultf,active)
      call addfld('ZMDT    ','K/S     ',nf,nmultf,active)
      call addfld('ZMDQ    ','KG/KGS  ',nf,nmultf,active)
      call addfld('CMFDQR  ','KG/KGS  ',nf,nmultf,inactive)
      call addfld('CMFMC   ','KG/M2S  ',nf,nmulth,active)
      call addfld('CMFSL   ','W/M2    ',nf,nmulth,inactive)
      call addfld('CMFLQ   ','W/M2    ',nf,nmulth,inactive)
      call addfld('CNVCLD  ','FRACTION',nf,nsingl,active)
      call addfld('FSDS    ','W/M2    ',nf,nsingl,active)
      call addfld('O3VMR   ','FRACTION',nf,nmultf,inactive)
      call addfld('T10M    ','K       ',nf,nsingl,active)
      call addfld('Q10M    ','KG/KG   ',nf,nsingl,active)
      call addfld('U10M    ','M/S     ',nf,nsingl,active)
      call addfld('V10M    ','M/S     ',nf,nsingl,active)
      call addfld('T2M     ','K       ',nf,nsingl,active)
      call addfld('Q2M     ','KG/KG   ',nf,nsingl,active)
      call addfld('U2M     ','M/S     ',nf,nsingl,active)
      call addfld('V2M     ','M/S     ',nf,nsingl,active)
C
C Fields which go out active monthly average history files
C
      if (ninavg(1).eq.'Q') then
        flag = active
      else
        flag = inactive
      end if
      call addfld('VT      ','KM/S    ',nf,nmultf,flag)
      call addfld('VZ      ','M3/S3   ',nf,nmultf,flag)
      call addfld('VQ      ','M/S     ',nf,nmultf,flag)
      call addfld('VVPUU   ','M2/S2   ',nf,nmultf,flag)
      call addfld('RELHUM  ','FRACTION',nf,nmultf,flag)
      call addfld('Z3      ','M       ',nf,nmultf,flag)
      call addfld('MQ      ','KG/M2   ',nf,nmultf,flag)
      call addfld('PSL     ','PASCALS ',nf,nsingl,flag)
      call addfld('OMEGAUP ','PA/S    ',nf,nmultf,flag)
      call addfld('NUMOMGUP','FRACTION',nf,nmultf,flag)
      call addfld('CLOUDUP ','FRACTION',nf,nmultf,flag)
      call addfld('DPSLON  ','M/S2    ',nf,nmultf,flag)
      call addfld('DPSLAT  ','M/S2    ',nf,nmultf,flag)
C
C Sort the master field list by summation of ichar's of 1st 4 characters
C pichsum = 127*4.  Minimum ichar output is 0. Maximum is 127.
C
      islocdum(0) = 0
      isloccc(0) = 0
      do i=0,pichsum
         nfldscdu(i) = 0
      end do
C
C Count up number of fields with a given sum of 4 ichar's
C
      do i=1,nf
         ipos = ichar(fieldn(1,i)(1:1)) + ichar(fieldn(1,i)(2:2)) +
     $          ichar(fieldn(1,i)(3:3)) + ichar(fieldn(1,i)(4:4))
         nfldscdu(ipos) = nfldscdu(ipos) + 1
      end do
C
C Define starting location for each sorted sum
C
      do i=1,pichsum
         islocdum(i) = islocdum(i-1) + nfldscdu(i-1)
         isloccc(i) = islocdum(i)
      end do
C
C Sort the hashed field names by increasing sum
C
      do i=1,nf
         ipos = ichar(fieldn(1,i)(1:1)) + ichar(fieldn(1,i)(2:2)) +
     $          ichar(fieldn(1,i)(3:3)) + ichar(fieldn(1,i)(4:4))
         isloccc(ipos) = isloccc(ipos) + 1
         fieldns(1,isloccc(ipos)) = fieldn(1,i)
         fieldns(2,isloccc(ipos)) = fieldn(2,i)
         ifldss(1,isloccc(ipos)) = iflds(1,i)
         ifldss(3,isloccc(ipos)) = iflds(3,i)
      end do
C
C Copy the local sorted list into the permanent list
C
      do i=1,nf
         fieldn(1,i) = fieldns(1,i)
         fieldn(2,i) = fieldns(2,i)
         iflds(1,i) = ifldss(1,i)
         iflds(3,i) = ifldss(3,i)
      end do
C
C Check number of fields against max for master list.
C
      if (nf.gt.pflds) then
         write(6,*)'BLDFLD: Too many fields -- pflds,nf=',pflds,nf
         call endrun
      end if
      if (lprint) then
         write(6,*)' ******* MASTER FIELD LIST *******'
         do j=1,nf
            write(6,9000)j,fieldn(1,j),iflds(1,j),iflds(3,j),
     $                   fieldn(2,j)
         end do
      end if
C
      return
C
 9000 format (i5,8x,a8,i3,i3,5x,a8)
      end
 
