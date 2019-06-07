#include <misc.h>
#include <preproc.h>

subroutine iniTimeVar (readini, eccen, obliqr, lambm0 , mvelpp)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize the following time varying variables:
!
!   water      : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
!   snow       : snowdp, snowage, snl, dz, z, zi 
!   temperature: t_soisno, t_veg, t_grnd
!
! Note - h2osoi_vol is needed by clm_soilalb -this is not needed on 
! restart since it is computed before the soil albedo computation is 
! called
! 
! Note -  remaining variables are initialized by calls to ecosystem 
! dynamics and albedo subroutines. 
!
! Method: 
! Initial data is saved to instantaneous initial data files for each of
! the [maxpatch] subgrid patches for each of the [numland] land points. 
! If a subgrid patch is not active (e.g., 3 patches rather 
! than [maxpatch]), the inactive subgrid patches have data values for the 
! first subgrid patch. 
! This way, as long as the land mask DOES NOT change among runs 
! (i.e., [numland] is the same), an initial data file can be used in 
! numerous experiments even if the surface types (and hence [numpatch]) 
! differ
! 
!
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use clm_varder
  use clm_varpar  , only : nlevsoi, nlevsno, nlevlak
  use clm_varmap  , only : begpatch, endpatch
  use clm_varcon  , only : bdsno, istice, istwet, istsoil, denice, denh2o, tfrz, spval 
  use inicFileMod , only : type_inidat, inicrd, histrd 
  use shr_sys_mod , only : shr_sys_abort
  use spmdMod     , only : masterproc
  use time_manager, only : get_nstep, get_curr_calday    
#if (defined SPMD)
  use mpishorthand, only : mpicom, mpichar
#endif
  implicit none

! ------------------------ arguments ---------------------------------
  logical , intent(in) :: readini  !true if read in initial data set
  real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox long. of perihelion + pi (radians)
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  integer  :: i,j,k,l,m,n               !loop indices
  integer  :: ier                       !MPI return code
  real(r8) :: pi                        !3.14159...
  logical  :: doalb                     !true => albedo time step
  character(len=16) :: initype          !type of initial dataset
  real(r8) :: calday                    !calendar day
! --------------------------------------------------------------------

! ----------------------------------------------------------------------
! Initialize water and temperature based on:
! o readini = true : read initial data set -- requires netCDF codes
! o readini = false: arbitrary initialization
! ----------------------------------------------------------------------

  if (readini) then

     if ( masterproc ) write (6,*) 'Reading initial data '

! determine format of initial data 

     call type_inidat(initype)
#if (defined SPMD)
    call mpi_bcast (initype, len(initype), mpichar, 0, mpicom, ier)
#endif

! open and read data from netCDF initial dataset or history file

     if (trim(initype) == 'INICFILE') then
        call inicrd () 
     else if (trim(initype) == 'HISTFILE') then
        call histrd ()
     else
        call shr_sys_abort('initial data type is limited to INIC or HIST file only')
     endif

! Determine volumetric soil water

     do k = begpatch,endpatch
        do j = 1,nlevsoi
           clm(k)%h2osoi_vol(j) = clm(k)%h2osoi_liq(j)/(clm(k)%dz(j)*denh2o) &
                                + clm(k)%h2osoi_ice(j)/(clm(k)%dz(j)*denice)
        end do
     end do

  else

     if ( masterproc ) write (6,*) 'Setting initial data to non-spun up values'

! ========================================================================
! Set snow water 
! ========================================================================

! NOTE: h2ocan, h2osno, snowdp and snowage has valid values everywhere

     do k = begpatch,endpatch
        clm(k)%h2ocan = 0.
        if (clm(k)%itypwat == istice) then
           clm(k)%h2osno = 1000.
        else
           clm(k)%h2osno = 0.
        endif
        clm(k)%snowdp  = clm(k)%h2osno/bdsno
        clm(k)%snowage = 0.                
     end do
        
! ========================================================================
! Set snow layer number, depth and thickiness 
! ========================================================================

     call snowdp2lev ()

! ========================================================================
! Set snow/soil temperature
! ========================================================================
        
! NOTE: 
! t_soisno only has valid values over non-lake
! t_lake   only has valid values over lake
! t_grnd has valid values over all land
! t_veg  has valid values over all land

     do k = begpatch,endpatch

        clm(k)%t_soisno(-nlevsno+1:nlevsoi) = spval
        clm(k)%t_lake(1:nlevlak) = spval
        clm(k)%t_veg = 283.

        if (.not. clm(k)%lakpoi) then  !not lake
           clm(k)%t_soisno(-nlevsno+1:0) = spval
           if (clm(k)%snl < 0) then    !snow layer temperatures
              do i = clm(k)%snl+1, 0
                 clm(k)%t_soisno(i) = 250.
              enddo
           endif
           do i = 1, nlevsoi
              if (clm(k)%itypwat == istice) then
                 clm(k)%t_soisno(i) = 250.
              else if (clm(k)%itypwat == istwet) then
                 clm(k)%t_soisno(i) = 277.
              else
                 clm(k)%t_soisno(i) = 283.
              endif
           end do
           clm(k)%t_grnd = clm(k)%t_soisno(clm(k)%snl+1)
        else                           !lake
           clm(k)%t_lake(1:nlevlak) = 277.
           clm(k)%t_grnd = clm(k)%t_lake(1)
        endif

     end do

! ========================================================================
! Set snow/soil ice and liquid mass
! ========================================================================
        
! volumetric water is set first and liquid content and ice lens are
! then obtained
! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values 
! over soil

     do k = begpatch, endpatch
        clm(k)%h2osoi_vol(         1:nlevsoi) = spval
        clm(k)%h2osoi_liq(-nlevsno+1:nlevsoi) = spval
        clm(k)%h2osoi_ice(-nlevsno+1:nlevsoi) = spval

        if (.not. clm(k)%lakpoi) then  !not lake
           ! volumetric water
           do i = 1,nlevsoi
              if (clm(k)%itypwat == istsoil) then 
                 clm(k)%h2osoi_vol(i) = 0.3_r8
              else
                 clm(k)%h2osoi_vol(i) = 1.0_r8
              endif
              clm(k)%h2osoi_vol(i) = min(clm(k)%h2osoi_vol(i),clm(k)%watsat(i))
           end do
  
           ! liquid water and ice (note that dz is in meters below)
           if (clm(k)%snl < 0) then    !snow 
              do i = clm(k)%snl+1, 0
                 clm(k)%h2osoi_ice(i) = clm(k)%dz(i)*250.
                 clm(k)%h2osoi_liq(i) = 0.
              enddo
           endif
           do i = 1, nlevsoi           !soil layers
              if (clm(k)%t_soisno(i) <= tfrz) then
                 clm(k)%h2osoi_ice(i) = clm(k)%dz(i)*denice*clm(k)%h2osoi_vol(i) 
                 clm(k)%h2osoi_liq(i) = 0.
              else
                 clm(k)%h2osoi_ice(i) = 0.
                 clm(k)%h2osoi_liq(i) = clm(k)%dz(i)*denh2o*clm(k)%h2osoi_vol(i)
              endif
           enddo
        endif
     end do

  end if  ! end of arbitrary initialization if-block

! ========================================================================
! Remaining variables are initialized by calls to ecosystem dynamics and
! albedo subroutines. 
! Note: elai, esai, frac_veg_nosno are computed in Ecosysdyn and needed
! by Fwet and SurfaceAlbedo
! Note: fwet is needed in routine clm_twostream (called by clm_surfalb)
! ========================================================================

  calday = get_curr_calday()
  doalb = .true.

#if (defined DGVM)
  call iniTimeConstDGVM()
#endif

!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch, endpatch

     clm(k)%nstep = get_nstep()

     call EcosystemDyn (clm(k), doalb, .false.)

     clm(k)%frac_sno = clm(k)%snowdp/(10.*clm(k)%zlnd + clm(k)%snowdp)  

     clm(k)%frac_veg_nosno = clm(k)%frac_veg_nosno_alb

     call Fwet(clm(k))

     call SurfaceAlbedo (clm(k), calday, eccen, obliqr, lambm0, mvelpp)

  end do
!$OMP END PARALLEL DO

  return
end subroutine iniTimeVar
