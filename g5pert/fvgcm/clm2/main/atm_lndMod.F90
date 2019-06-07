#include <misc.h>
#include <preproc.h>
#include <params.h>


module atm_lndMod

#if (defined COUP_CAM)

!-----------------------------------------------------------------------
!
! Purpose:
! Atm - Land interface module
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

!JDR Replace CCM4 modules with fvCCM3 include files (comsrf & pmgrid)
! Reformat and rename for clm2 some fvCCM3 includes written in fixed form

  use precision
!  use pmgrid, only: plon, plond, plat
!  use tracers, only: pcnst, pnats
!  use rgrid, only: nlon
!  use ppgrid, only: pcols, begchunk, endchunk
!  use phys_grid
!  use comsrf, only : oro, ts, asdir, aldir, asdif, aldif, snowh, lwup, &
!                     wsx, wsy, lhf, shf, cflx, tref, &
!                     zbot, ubot, vbot, tbot, thbot, qbot, pbot, flwds, &
!                     precsc, precsl, precc, precl, soll, sols, solld, solsd 
!  use history, only : caseid, ctitle, inithist, nhtfrq, mfilt
  use shr_const_mod, only: SHR_CONST_PI
!JDR added for clm2 calendar
!  use clm_varctl, only: ymdsec
  implicit none

#include <pmgrid.h>
#include <comsrf.h>


  private              ! By default make data private
  integer :: landmask(plon,plat) !2d land mask
  integer :: nlon(plat)   ! from rgrid
  
  character(len=17) caseid      ! Case identifier (max 16 chars)
  character(len=80) ctitle      ! Case title
  character(len=8)  inithist    ! If set to 'MONTHLY' or 'YEARLY' then write IC file 
  integer, parameter :: ptapes = 6             ! max number of tapes
  integer :: i
  integer :: nhtfrq(ptapes)     ! history write frequency (0 monthly)
  integer :: mfilt(ptapes)      ! number of time samples per tape

  public atmlnd_ini, atmlnd_drv, is_lsm   ! Public interfaces


!===============================================================================
CONTAINS
!===============================================================================

  subroutine atmlnd_ini()

    use initializeMod, only : initialize           !initialization of clm 
    use error_messages, only: alloc_err
    use time_manager, only: get_nstep
#if ( defined SPMD )
    use mpishorthand
#endif

#include <commap.h>
#include <comsol.h>
#include <comctl.h>
#include <commss.h>

!-----------------------------------------------------------------------
! Initialize land surface model and obtain relevant atmospheric model 
! arrays back from (i.e. albedos, surface temperature and snow cover over land)
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer  :: i,lat,n !indices
    integer  :: istat               !error return
    real(r8) :: landfrac(plon,plat) !2d fractional land
    real(r8) :: latixy(plon,plat)   !2d latitude  grid (degrees)
    real(r8) :: longxy(plon,plat)   !2d longitude grid (degrees)
    real(r8) :: pi,dp
    integer  :: nstep
    
!-----------------------------------------------------------------------

! Time management variables.

    nstep = get_nstep()

#if (defined SPMD)
    call mpibcast (oro_glob, size(oro_glob), mpir8, 0, mpicom)
#endif

! JDR Assign ccm4 timing variables to equivalent fvccm3 variables

    inithist='NONE'

    do lat = 1,plat
       nlon(lat)=plon
    enddo

    pi = SHR_CONST_PI
    longxy(:,:) = 1.e36
    dp = 180./(plat-1)
    do lat = 1,plat
       do i = 1,nlon(lat)
          longxy(i,lat) = (i-1)*360.0/nlon(lat)
           latixy(i,lat) = -90.0 + (lat-1)*dp
          if (nint(oro(i,lat)) == 1) then
             landmask(i,lat) = 1
             landfrac(i,lat) = 1.0
          else
             landmask(i,lat) = 0
             landfrac(i,lat) = 0.
          endif
       end do
    end do

    call  initialize(eccen    , obliqr   , lambm0  , mvelpp  , caseid  , &
                     ctitle   , nsrest   , nstep   , iradsw  , inithist, &
                     nhtfrq(1), mfilt(1) , longxy  , latixy  , nlon    , &
                     landmask , landfrac , irt)

! For initial run only - get 2d data back from land model (Note that 
! in SPMD case, only masterproc contains valid recv2d data) and 
! split 2d data into appropriate arrays contained in module comsrf. 

    if (nstep == 0) then

!       call lnd_to_ccm_mapping_ini(recv2d)
!       call scatter_field_to_chunk(1,nrecv_ccm,1,plon,recv2d,recv2d_chunk)

!       do lchnk=begchunk,endchunk
!          ncols = get_ncols_p(lchnk)
!          do i=1,ncols
!             if (landmask_chunk(i,lchnk) == 1) then
!                ts(i,lchnk)    = recv2d_chunk(i, 1,lchnk) 
!                asdir(i,lchnk) = recv2d_chunk(i, 2,lchnk) 
!                aldir(i,lchnk) = recv2d_chunk(i, 3,lchnk) 
!                asdif(i,lchnk) = recv2d_chunk(i, 4,lchnk) 
!                aldif(i,lchnk) = recv2d_chunk(i, 5,lchnk) 
!                snowh(i,lchnk) = recv2d_chunk(i, 6,lchnk) 
!                lwup(i,lchnk)  = recv2d_chunk(i,11,lchnk) 
!             endif
!          end do
!       end do

    endif
    
    return
  end subroutine atmlnd_ini

!===============================================================================

  subroutine atmlnd_drv (iradsw, eccen, obliqr, lambm0, mvelpp)

!-----------------------------------------------------------------------
! Pack data to be sent to land model into a single array. 
! Send data to land model and call land model driver. 
! Receive data back from land model in a single array.
! Unpack this data into component arrays. 
! NOTE: component arrays are contained in module comsrf.
!-----------------------------------------------------------------------

#if ( defined SPMD )
    use mpishorthand
#endif
    use clm_varder
    use clm_varcon
    use clm_varsur
    use clm_varmap , only : numpatch,patchvec
    use time_manager, only: get_nstep, advance_timestep

!---------------------------Arguments----------------------------------- 
    integer , intent(in) :: iradsw   !Iteration frequency for shortwave radiation
    real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
    real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
    real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox longitude of perihelion + pi (radians)
    integer              :: nstep
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer :: i,lat,m !indices
    logical doalb          !true if surface albedo calculation time step
!JDR Added variables for mapping
    integer  :: j,k
    real(r8) :: wt                          !remapping weight
    real(r8) :: wtsum(plon,plat)

!-----------------------------------------------------------------------

! -----------------------------------------------------------------
! Determine doalb
! [doalb] is a logical variable that is true when the next time
! step is a radiation time step. This allows for the fact that
! an atmospheric model may not do the radiative calculations 
! every time step. For example:
!      nstep dorad doalb
!        1     F     F
!        2     F     T
!        3     T     F
!        4     F     F
!        5     F     T
!        6     T     F
! The following expression for doalb is for example only (it is 
! specific to the NCAR CAM). This variable must be calculated
! appropriately for the host atmospheric model
! -----------------------------------------------------------------
    nstep = get_nstep()
    
    doalb = iradsw==1 .or. (mod(nstep,iradsw)==0 .and. nstep+1/=1)

! Condense the 2d atmospheric data needed by the land surface model into 
! one array. Note that precc and precl precipitation rates are in units 
! of m/sec. They are turned into fluxes by multiplying by 1000 kg/m^3.
! JDR Added atm to lnd mapping of forcing

    do k = 1, numpatch
       
       i=patchvec%ixy(k)
       j=patchvec%jxy(k)

       clm(k)%forc_hgt      = zbot(i,j)          !zgcmxy  Atm state m
       clm(k)%forc_u        = ubot(i,j)          !forc_uxy  Atm state m/s
       clm(k)%forc_v        = vbot(i,j)          !forc_vxy  Atm state m/s
       clm(k)%forc_th       = thbot(i,j)         !forc_thxy Atm state K
       clm(k)%forc_q        = qbot(i,j)          !forc_qxy  Atm state kg/kg
       clm(k)%forc_pbot     = pbot(i,j)          !ptcmxy  Atm state Pa
       clm(k)%forc_t        = tbot(i,j)          !forc_txy  Atm state K
       clm(k)%forc_lwrad    = flwds(i,j)         !flwdsxy Atm flux  W/m^2
       clm(k)%forc_solad(2) = soll(i,j)          !forc_sollxy  Atm flux  W/m^2
       clm(k)%forc_solad(1) = sols(i,j)          !forc_solsxy  Atm flux  W/m^2 
       clm(k)%forc_solai(2) = solld(i,j)         !forc_solldxy Atm flux  W/m^2
       clm(k)%forc_solai(1) = solsd(i,j)         !forc_solsdxy Atm flux  W/m^2

       ! determine derived quantities

       clm(k)%forc_hgt_u = clm(k)%forc_hgt   !observational height of wind [m]            !derived
       clm(k)%forc_hgt_t = clm(k)%forc_hgt   !observational height of temperature [m]     !derived
       clm(k)%forc_hgt_q = clm(k)%forc_hgt   !observational height of humidity [m]        !derived
       clm(k)%forc_vp    = clm(k)%forc_q*clm(k)%forc_pbot / (0.622+0.378*clm(k)%forc_q)   !derived
       clm(k)%forc_rho   = (clm(k)%forc_pbot-0.378*clm(k)%forc_vp) / (rair*clm(k)%forc_t) !derived 
       clm(k)%forc_co2   = pco2*clm(k)%forc_pbot                                          !derived
       clm(k)%forc_o2    = po2*clm(k)%forc_pbot                                           !derived

       ! Determine precipitation needed by clm
! Set upper limit of air temperature for snowfall at 275.65K.
! This cut-off was selected based on Fig. 1, Plate 3-1, of Snow
! Hydrology (1956).

       if (tbot(i,j) > (tfrz + tcrit)) then
          clm(k)%forc_rain = (precc(i,j) + precl(i,j))*1000.    !mm/s
          clm(k)%forc_snow = 0.0_r8
       else
          clm(k)%forc_snow = (precc(i,j) + precl(i,j))*1000.    !mm/s 
          clm(k)%forc_rain = 0.0_r8
       endif

       if ( clm(k)%forc_snow > 0.0_r8  .and. clm(k)%forc_rain > 0.0_r8 ) then
          write(6,*) 'kpatch= ',k,' snow= ',clm(k)%forc_snow,' rain= ',clm(k)%forc_rain, &
               ' CLM cannot currently handle both non-zero rain and snow'
          call endrun
       elseif (clm(k)%forc_rain > 0.0_r8) then
          clm(k)%itypprc = 1
       elseif (clm(k)%forc_snow > 0.0_r8) then
          clm(k)%itypprc = 2
       else
          clm(k)%itypprc = 0
       endif

    end do

! Call land model driver
    call driver (doalb, eccen, obliqr, lambm0, mvelpp)

! increment time step

    call advance_timestep()
  
! Convert one dimensional land model output data to two dimensional atm data 
! JDR Added lnd to atm mapping
    do k= 1, numpatch          
       if (patchvec%wtxy(k) /= 0.) then
          i  = patchvec%ixy(k)    
          j  = patchvec%jxy(k)    
          wt = patchvec%wtxy(k) 
         
          ts(i,j)      = 0.              !tsxy 
          asdir(i,j)   = 0.              !asdir
          aldir(i,j)   = 0.              !aldir
          asdif(i,j)   = 0.              !asdif
          aldif(i,j)   = 0.              !aldif
          snowh(i,j)   = 0.              !snow (convert mm->m)
          wsx(i,j)     = 0.              !taux 
          wsy(i,j)     = 0.              !tauy
          lhf(i,j)     = 0.              !lhflx 
          shf(i,j)     = 0.              !shflx 
          lwup(i,j)    = 0.              !lwup
          cflx(i,1,j)  = 0.              !qflx 
          tref(i,j)    = 0.              !tref
       endif
    end do
    
    do k= 1, numpatch          
       if (patchvec%wtxy(k) /= 0.) then
          i  = patchvec%ixy(k)    
          j  = patchvec%jxy(k)    
          wt = patchvec%wtxy(k) 
         
          ts(i,j)      = ts(i,j)     + clm(k)%t_rad*wt              !tsxy 
          asdir(i,j)   = asdir(i,j)  + clm(k)%albd(1)*wt            !asdir
          aldir(i,j)   = aldir(i,j)  + clm(k)%albd(2)*wt            !aldir
          asdif(i,j)   = asdif(i,j)  + clm(k)%albi(1)*wt            !asdif
          aldif(i,j)   = aldif(i,j)  + clm(k)%albi(2)*wt            !aldif
          snowh(i,j)   = snowh(i,j)  + (clm(k)%h2osno/1000.)*wt     !snow (convert mm->m)
          wsx(i,j)     = wsx(i,j)    + clm(k)%taux*wt               !taux 
          wsy(i,j)     = wsy(i,j)    + clm(k)%tauy*wt               !tauy
          lhf(i,j)     = lhf(i,j)    + clm(k)%eflx_lh_tot*wt        !lhflx 
          shf(i,j)     = shf(i,j)    + clm(k)%eflx_sh_tot*wt        !shflx 
          lwup(i,j)    = lwup(i,j)   + clm(k)%eflx_lwrad_out*wt     !lwup
          cflx(i,1,j)  = cflx(i,1,j) + clm(k)%qflx_evap_tot*wt      !qflx 
          tref(i,j)    = tref(i,j)   + clm(k)%t_ref2m*wt            !tref
       endif
    end do
    
    return
  end subroutine atmlnd_drv

   logical function is_lsm ( )
     is_lsm = .false.
     return
   end function is_lsm


!===============================================================================

#endif        

end module atm_lndMod

