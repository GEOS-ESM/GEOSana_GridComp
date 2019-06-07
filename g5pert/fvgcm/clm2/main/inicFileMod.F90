#include <misc.h>
#include <preproc.h>

module inicFileMod

!----------------------------------------------------------------------- 
! Purpose: 
! read and writes initial data netCDF history files
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use clm_varder
  use clm_varmap  , only : begpatch, endpatch, numland, numpatch
  use clm_varpar  , only : nlevsno, nlevsoi, nlevlak, maxpatch, rtmlon, rtmlat
  use clm_varmap  , only : landvec, patchvec, begland, endland
  use clm_varcon  , only : spval
  use fileutils   , only : getfil
#if (defined SPMD)
  use spmdMod     , only : masterproc, npes, compute_mpigs_patch, compute_mpigs_land
  use mpishorthand, only : mpir8, mpiint, mpilog, mpicom
#else
  use spmdMod     , only : masterproc
#endif
#if (defined RTM)
  use RtmMod      , only : volr
#endif
  implicit none

! loop variables

  integer, private  :: i,j,k,l,m,n         !loop indices

! netcdf data

  integer, private  :: ncid                !netCDF dataset id
  integer, private  :: dimid               !netCDF dimension id 
  integer, private  :: varid               !netCDF variable id

! input dataset dimensions

  integer, private  :: numland_dim         !value for [numland] from dataset
  integer, private  :: maxpatch_dim        !value for [maxpatch] from dataset
  integer, private  :: nlevsoi_dim         !value for [nlevsoi] from dataset
  integer, private  :: nlevsno_dim         !value for [nlevsno] from dataset
  integer, private  :: nlevtot_dim         !number of total (snow+soil) levels from dataset  
  integer, private  :: rtmlon_dim          !number of RTM longitudes
  integer, private  :: rtmlat_dim          !number of RTM latitudes

! methods

  public  :: do_inicwrite 
  private :: patch_to_land 
  private :: land_to_patch
  private :: hist_to_clm
  private :: set_init_filename

  INTERFACE patch_to_land
     MODULE procedure patch_to_land_1d_int
     MODULE procedure patch_to_land_1d_real
     MODULE procedure patch_to_land_2d_real
  END INTERFACE

  INTERFACE land_to_patch
     MODULE procedure land_to_patch_1d_int
     MODULE procedure land_to_patch_1d_real
     MODULE procedure land_to_patch_2d_real
  END INTERFACE

  INTERFACE hist_to_clm
     MODULE procedure hist_to_clm_1d
     MODULE procedure hist_to_clm_2d
  END INTERFACE

  SAVE

!=======================================================================
CONTAINS
!=======================================================================

  subroutine type_inidat (initype)
    
!----------------------------------------------------------------------- 
! 
! Purpose: 
! determine type of initial dataset (instantaneous inic file or
! history file)
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
  
    use precision
    use clm_varctl , only : finidat
    implicit none
    include 'netcdf.inc'

! ------------------------ arguments ---------------------------------
    character(len=16), intent(out) :: initype  !type of initial dataset
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
   integer ret                       !NetCDF return code
   character(len=256) :: locfn       !local file name
! --------------------------------------------------------------------

    if (masterproc) then

! open netCDF initial data file and read data 

       call getfil (finidat, locfn, 0)
       call wrap_open (locfn, NF_NOWRITE, ncid)

! check if file has dimensional history file output

       ret = nf_inq_dimid (ncid, 'patch', dimid)
       if (ret==NF_NOERR) then
          initype = 'HISTFILE'
          RETURN
       endif

! check file has numland x numpatch instantaneous output

       ret = nf_inq_dimid (ncid, 'numland', dimid)
       if (ret==NF_NOERR) then
          initype = 'INICFILE'
          RETURN
       endif

! file has non-supported form - stop the run 

       write(6,*)'TYPE_INIDAT: initial dataset has non-supported form'
       call endrun

! close the file

       call wrap_close(ncid)

    endif

  end subroutine type_inidat

!=======================================================================

  subroutine histrd()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! read initial data from netCDF history file for variables
! SNOWDP, SNOWAGE, TV, TV, TSOI, TSNOW, SOILLIQ, SOILICE, 
! SNOWLIQ, SNOWICE, H2OCAN, H2OSNO
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use precision
    use clm_varmap, only : begpatch, endpatch, numland, numpatch
    use clm_varcon, only : istsoil, tfrz, denice, denh2o
    use clm_varctl, only : finidat
    implicit none
    include 'netcdf.inc'

! ------------------------ local variables -----------------------------
    character(len=256) :: locfn            !local file name
    integer ndim                           !temporary  
    real(r8) mlevini(numpatch,nlevsoi)     !temporary
    real(r8) slevini(numpatch)             !temporary 
    integer begsing(2),lensing(2)          !single level array section spec
    integer begmult(3),lenmult(3)          !multi  level array section spec
    integer ntim                           !number of time samples on tape 
    real(r8) dzsum                         !vertical sum of snow layers
! --------------------------------------------------------------------

! Open netCDF data file and read data and check input dimensions

    if (masterproc) then
       call getfil (finidat, locfn, 0)
       call wrap_open (locfn, nf_nowrite, ncid)

       call wrap_inq_dimid (ncid, 'patch', dimid)
       call wrap_inq_dimlen (ncid, dimid, ndim)
       if (ndim /= numpatch) then
          write (6,*) 'HISTRD error: numpatch values disagree'
          write (6,*) 'finidat numpatch = ',ndim,' model numpatch = ',numpatch
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'levsoi', dimid)
       call wrap_inq_dimlen (ncid, dimid, ndim)
       if (ndim /= nlevsoi) then
          write (6,*) 'HISTRD error: nlevsoi values disagree'
          write (6,*) 'finidat nlevsoi = ',ndim,' model nlevsoi = ',nlevsoi
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'time', dimid)
       call wrap_inq_dimlen (ncid, dimid, ntim)
    endif

! obtain data

    begsing(1) = 1    ; lensing(1) = numpatch
    begsing(2) = ntim ; lensing(2) = 1

    begmult(1) = 1    ; lenmult(1) = numpatch
    begmult(2) = 1    ; lenmult(2) = nlevsoi
    begmult(3) = ntim ; lenmult(3) = 1

    if (masterproc) then
       call wrap_inq_varid (ncid, 'SOILLIQ', varid)
       call wrap_get_vara_realx (ncid, varid, begmult, lenmult, mlevini)
    endif
    call hist_to_clm(mlevini, nlevsoi)
    do j = 1,nlevsoi
       clm(begpatch:endpatch)%h2osoi_liq(j) = mlevini(begpatch:endpatch,j)
    end do

    if (masterproc) then
       call wrap_inq_varid (ncid, 'SOILICE', varid)
       call wrap_get_vara_realx (ncid, varid, begmult, lenmult, mlevini)
    endif
    call hist_to_clm(mlevini, nlevsoi)
    do j = 1,nlevsoi
       clm(begpatch:endpatch)%h2osoi_ice(j) = mlevini(begpatch:endpatch,j)
    end do

    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNOWLIQ', varid)
       call wrap_get_vara_realx (ncid, varid, begsing, lensing, slevini)
    endif
    call hist_to_clm(slevini)
    clm(begpatch:endpatch)%snowliq = slevini(begpatch:endpatch)

    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNOWICE', varid)
       call wrap_get_vara_realx (ncid, varid, begsing, lensing, slevini)
    endif
    call hist_to_clm(slevini)
    clm(begpatch:endpatch)%snowice = slevini(begpatch:endpatch)

    if (masterproc) then
       call wrap_inq_varid (ncid, 'TSOI', varid)
       call wrap_get_vara_realx (ncid, varid, begmult, lenmult, mlevini)
    endif
    call hist_to_clm(mlevini, nlevsoi)
    do j = 1,nlevsoi
       clm(begpatch:endpatch)%t_soisno(j) = mlevini(begpatch:endpatch,j)
    end do

    if (masterproc) then
       call wrap_inq_varid (ncid, 'TLAKE', varid)
       call wrap_get_vara_realx (ncid, varid, begmult, lenmult, mlevini)
    endif
    call hist_to_clm(mlevini, nlevlak)
    do j = 1,nlevlak
       clm(begpatch:endpatch)%t_lake(j) = mlevini(begpatch:endpatch,j)
    end do

    if (masterproc) then
       call wrap_inq_varid (ncid, 'TSNOW', varid)
       call wrap_get_vara_realx (ncid, varid, begsing, lensing, slevini)
    endif
    call hist_to_clm(slevini)
    clm(begpatch:endpatch)%t_snow = slevini(begpatch:endpatch)
    
    if (masterproc) then
       call wrap_inq_varid (ncid, 'TV', varid)
       call wrap_get_vara_realx (ncid, varid, begsing, lensing, slevini)
    endif
    call hist_to_clm(slevini)
    clm(begpatch:endpatch)%t_veg = slevini(begpatch:endpatch)

    if (masterproc) then
       call wrap_inq_varid (ncid, 'TG', varid)
       call wrap_get_var_realx (ncid, varid, slevini)
    endif
    call hist_to_clm(slevini)
    clm(begpatch:endpatch)%t_grnd = slevini(begpatch:endpatch)

    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OCAN', varid)
       call wrap_get_vara_realx (ncid, varid, begsing, lensing, slevini)
    endif
    call hist_to_clm(slevini)
    clm(begpatch:endpatch)%h2ocan = slevini(begpatch:endpatch)
    
    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OSNO', varid)
       call wrap_get_vara_realx (ncid, varid, begsing, lensing, slevini)
    endif
    call hist_to_clm(slevini)
    clm(begpatch:endpatch)%h2osno = slevini(begpatch:endpatch)

    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNOWDP', varid)
       call wrap_get_vara_realx (ncid, varid, begsing, lensing, slevini)
    endif
    call hist_to_clm(slevini)
    clm(begpatch:endpatch)%snowdp = slevini(begpatch:endpatch)

    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNOWAGE', varid)
       call wrap_get_vara_realx (ncid, varid, begsing, lensing, slevini)
    endif
    call hist_to_clm(slevini)
    clm(begpatch:endpatch)%snowage = slevini(begpatch:endpatch)

    if (masterproc) then
       call wrap_close (ncid)
    endif

! create the snow layers using the input snowdp

    call snowdp2lev

! put the average t_snow, snow_ice and snow_liq into the soil snow layers
! WHERE THERE IS snow, REPLACE h2osno with sum of snowice and snowliq to
! have consistency between h2osno and possibly averaged input values

    do k = begpatch,endpatch
       clm(k)%t_soisno  (-nlevsno+1:0) = spval
       clm(k)%h2osoi_liq(-nlevsno+1:0) = spval
       clm(k)%h2osoi_ice(-nlevsno+1:0) = spval
       if ((.not. clm(k)%lakpoi) .and. (clm(k)%snl < 0)) then  
          dzsum = 0.
          do j = clm(k)%snl+1, 0
             dzsum = dzsum + clm(k)%dz(j)
          end do
          do j = clm(k)%snl+1, 0
             clm(k)%t_soisno(j)   = clm(k)%t_snow
             clm(k)%h2osoi_ice(j) = clm(k)%snowice*clm(k)%dz(j)/dzsum
             clm(k)%h2osoi_liq(j) = clm(k)%snowliq*clm(k)%dz(j)/dzsum
          enddo
          clm(k)%h2osno = clm(k)%snowice + clm(k)%snowliq
       endif
    end do

    return
  end subroutine histrd

!=======================================================================

  subroutine inicrd ()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! read initial data from netCDF instantaneous initial data history file 
! for variables:
!   snlsno, dzsno, zsno, zisno, h2ocan, h2osno, snowdp, snowage, 
!   h2osoi_liq, h2osoi_ice, t_veg, t_grnd, t_soisno, t_lake
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use precision
    use clm_varctl, only : finidat
    implicit none
    include 'netcdf.inc'

! ------------------------ local variables -----------------------------
    character(len=256) :: locfn   !local file name
    integer :: ndim               !input dimension       
    integer :: ret                !netcdf return code
#if (defined SPMD)
    integer :: numrecvv(0:npes-1) !vector of items to be received  
    integer :: displsv(0:npes-1)  !displacement vector
    integer :: numsend            !number of items to be sent
    integer :: ier                !mpi return code
#endif
    integer , allocatable :: ibuf1dl(:,:)
    integer , allocatable :: ibuf1dp(:)
    real(r8), allocatable :: rbuf1dl(:,:)
    real(r8), allocatable :: rbuf1dp(:)
    real(r8), allocatable :: rbuf2dl(:,:,:)
    real(r8), allocatable :: rbuf2dp(:,:)
! --------------------------------------------------------------------

! Open netCDF data file and read data

    if (masterproc) then

       call getfil (finidat, locfn, 0)
       call wrap_open (locfn, nf_nowrite, ncid)

! check input dimensions

       call wrap_inq_dimid (ncid, 'numland', dimid)
       call wrap_inq_dimlen (ncid, dimid, numland_dim)
       if (numland_dim /= numland) then
          write (6,*) 'INICRD error: numland values disagree'
          write (6,*) 'finidat numland = ',numland_dim,' model numland = ',numland
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'maxpatch', dimid)
       call wrap_inq_dimlen (ncid, dimid, maxpatch_dim)
       if (maxpatch_dim /= maxpatch) then
          write (6,*) 'INICRD error: maxpatch values disagree'
          write (6,*) 'finidat maxpatch = ',maxpatch_dim,' model maxpatch = ',maxpatch
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'nlevsno', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlevsno_dim)
       if (nlevsno_dim /= nlevsno) then
          write (6,*) 'INICRD error: nlevsno values disagree'
          write (6,*) 'finidat levsno = ',nlevsno_dim,' model nlevsno = ',nlevsno
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'nlevsoi', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlevsoi_dim)
       if (nlevsoi_dim /= nlevsoi) then
          write (6,*) 'INICRD error: nlevsoi values disagree'
          write (6,*) 'finidat nlevsoi = ',nlevsoi_dim,' model nlevsoi = ',nlevsoi
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'nlevtot', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlevtot_dim)
       if (nlevtot_dim /= nlevsoi+nlevsno) then
          write (6,*) 'INICRD error: nlevtot values disagree'
          write (6,*) 'finidat nlevtot = ',nlevtot_dim,' model nlevtot = ',nlevsno+nlevsoi
          call endrun 
       end if

#if (defined RTM)
       ret = nf_inq_dimid (ncid, 'rtmlon', dimid)
       if (ret == NF_NOERR) then
          call wrap_inq_dimlen (ncid, dimid, rtmlon_dim)
          if (rtmlon_dim /= rtmlon) then
             write (6,*) 'INICRD error: rtmlon values disagree'
             write (6,*) 'finidat rtmlon = ',rtmlon_dim,' model rtmlon = ',rtmlon
             call endrun
          end if
       endif

       ret = nf_inq_dimid (ncid, 'rtmlat', dimid)
       if (ret == NF_NOERR) then
          call wrap_inq_dimlen (ncid, dimid, rtmlat_dim)
          if (rtmlat_dim /= rtmlat) then
             write (6,*) 'INICRD error: rtmlat values disagree'
             write (6,*) 'finidat rtmlat = ',rtmlat_dim,' model rtmlat = ',rtmlat
             call endrun
          end if
       endif
#endif

    endif ! if-masterproc block

! Obtain data - for the snow interfaces, are only examing the snow 
! interfaces above zi=0 which is why zisno and zsno have the same 
! level dimension below

    allocate (rbuf1dl(numland,maxpatch))
    allocate (ibuf1dl(numland,maxpatch))
    allocate (rbuf1dp(begpatch:endpatch))
    allocate (ibuf1dp(begpatch:endpatch))
    
    ! Read in zisno
    ! NOTE: zi(0) is set to 0 in routine iniTimeConst
    allocate (rbuf2dp(-nlevsno+0:-1,begpatch:endpatch))
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+0:-1)) 
    if (masterproc) then
       call wrap_inq_varid (ncid, 'ZISNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno) 
    do k = begpatch,endpatch
       clm(k)%zi(-nlevsno+0:-1) = rbuf2dp(-nlevsno+0:-1,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in zsno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:0))
    allocate (rbuf2dp(-nlevsno+1: 0,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'ZSNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno)
    do k = begpatch,endpatch
       clm(k)%z (-nlevsno+1:0) = rbuf2dp(-nlevsno+1:0,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in dzsno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:0))
    allocate (rbuf2dp(-nlevsno+1: 0,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'DZSNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno)
    do k = begpatch,endpatch
       clm(k)%dz(-nlevsno+1:0) = rbuf2dp(-nlevsno+1:0,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Read in h2osoi_liq
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OSOI_LIQ_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno+nlevsoi)
    do k = begpatch,endpatch
       clm(k)%h2osoi_liq(-nlevsno+1:nlevsoi) = rbuf2dp(-nlevsno+1:nlevsoi,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in h2osoi_ice
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OSOI_ICE_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno+nlevsoi)
    do k = begpatch,endpatch
       clm(k)%h2osoi_ice(-nlevsno+1:nlevsoi) = rbuf2dp(-nlevsno+1:nlevsoi,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in t_soisno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'T_SOISNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno+nlevsoi)
    do k = begpatch,endpatch
       clm(k)%t_soisno(-nlevsno+1:nlevsoi) = rbuf2dp(-nlevsno+1:nlevsoi,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in t_lake
    allocate(rbuf2dl(numland,maxpatch,1:nlevlak))
    allocate(rbuf2dp(1:nlevlak,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'T_LAKE_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevlak)
    do k = begpatch,endpatch
       clm(k)%t_lake(1:nlevlak) = rbuf2dp(1:nlevlak,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in t_veg
    if (masterproc) then
       call wrap_inq_varid (ncid, 'T_VEG_INI', varid)
       call wrap_get_var_realx (ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%t_veg  = rbuf1dp(k)
    end do
    
    ! Read in t_grnd
    if (masterproc) then
       call wrap_inq_varid (ncid, 'T_GRND_INI', varid)
       call wrap_get_var_realx (ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%t_grnd = rbuf1dp(k)
    end do
    
    ! Read in h2ocan
    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OCAN_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%h2ocan = rbuf1dp(k)
    end do

    ! Read in h2osno
    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OSNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%h2osno = rbuf1dp(k)
    end do
    
    ! Read in snowdp
    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNOWDP_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%snowdp = rbuf1dp(k)
    end do
    
    ! Read in snowage
    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNOWAGE_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%snowage= rbuf1dp(k)
    end do
    
    ! Read in snlsno
    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNLSNO_INI', varid)
       call wrap_get_var_int(ncid, varid, ibuf1dl)
    endif
    call land_to_patch (ibuf1dl, ibuf1dp)
    do k = begpatch,endpatch
       clm(k)%snl = ibuf1dp(k)
    end do

#if (defined RTM)

    if (masterproc) then
       ret = nf_inq_varid (ncid, 'RTMVOLR', varid)
       if (ret == NF_NOERR) then
          write(6,*)'INICFILE: reading in river volr'
          call wrap_get_var_realx(ncid, varid, volr)
       endif
    endif

#endif

#if (defined DGVM)

    ! Read in itypveg
    if (masterproc) then
       call wrap_inq_varid (ncid, 'ITYPVEG_INI', varid)
       call wrap_get_var_int(ncid, varid, ibuf1dl)
    endif
    call land_to_patch (ibuf1dl, ibuf1dp)
    do k = begpatch,endpatch
       clm(k)%itypveg = ibuf1dp(k)
    end do

    ! Read in fpcgrid
    if (masterproc) then
       call wrap_inq_varid (ncid, 'FPCGRID_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%fpcgrid = rbuf1dp(k)
    end do

    ! Read in lai_ind
    if (masterproc) then
       call wrap_inq_varid (ncid, 'LAI_IND_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%lai_ind = rbuf1dp(k)
    end do

    ! Read in crownarea
    if (masterproc) then
       call wrap_inq_varid (ncid, 'CROWNAREA_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%crownarea = rbuf1dp(k)
    end do

    ! Read in litterag
    if (masterproc) then
       call wrap_inq_varid (ncid, 'LITTERAG_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%litterag = rbuf1dp(k)
    end do

    ! Read in litterbg
    if (masterproc) then
       call wrap_inq_varid (ncid, 'LITTERBG_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%litterbg = rbuf1dp(k)
    end do

    ! Read in cpool_fast
    if (masterproc) then
       call wrap_inq_varid (ncid, 'CPOOL_FAST_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%cpool_fast = rbuf1dp(k)
    end do

    ! Read in cpool_slow
    if (masterproc) then
       call wrap_inq_varid (ncid, 'CPOOL_SLOW_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%cpool_slow = rbuf1dp(k)
    end do

    ! Read in present
    if (masterproc) then
       call wrap_inq_varid (ncid, 'PRESENT_INI', varid)
       call wrap_get_var_int(ncid, varid, ibuf1dl)
    endif
    call land_to_patch (ibuf1dl, ibuf1dp)
    do k = begpatch,endpatch
       clm(k)%present = .false.
       if (ibuf1dp(k) == 1) clm(k)%present = .true.
    end do

    ! Read in nind
    if (masterproc) then
       call wrap_inq_varid (ncid, 'NIND_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%nind = rbuf1dp(k)
    end do

    ! Read in lm_ind
    if (masterproc) then
       call wrap_inq_varid (ncid, 'LM_IND_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%lm_ind = rbuf1dp(k)
    end do

    ! Read in sm_ind
    if (masterproc) then
       call wrap_inq_varid (ncid, 'SM_IND_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%sm_ind = rbuf1dp(k)
    end do

    ! Read in hm_ind
    if (masterproc) then
       call wrap_inq_varid (ncid, 'HM_IND_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%hm_ind = rbuf1dp(k)
    end do

    ! Read in rm_ind
    if (masterproc) then
       call wrap_inq_varid (ncid, 'RM_IND_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%rm_ind = rbuf1dp(k)
    end do

    ! Read in htop
    if (masterproc) then
       call wrap_inq_varid (ncid, 'HTOP_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%htop = rbuf1dp(k)
    end do

    ! Read in wtxy
    if (masterproc) then
       call wrap_inq_varid (ncid, 'WTXY_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    if (masterproc) then
       do k = 1,numpatch
          l = patchvec%land(k)
          m = patchvec%mxy(k)
          patchvec%wtxy(k) = rbuf1dl(l,m)
          landvec%wtxy(l,m) = rbuf1dl(l,m)
       end do
    else
       do k = begpatch,endpatch
          l = patchvec%land(k)
          m = patchvec%mxy(k)
          patchvec%wtxy(k) = rbuf1dp(k)
          landvec%wtxy(l,m) = rbuf1dp(k)
       end do
    endif
       
#endif

    deallocate (ibuf1dl)
    deallocate (rbuf1dl)
    deallocate (ibuf1dp)
    deallocate (rbuf1dp)

    return
  end subroutine inicrd

!=======================================================================

  subroutine inicwrt ()

!----------------------------------------------------------------------- 
! Purpose: 
! write initial data to netCDF history file
!
! Method: 
! 
! Author: Mariana Vertenstein
!-----------------------------------------------------------------------

    use precision
    use clm_varctl, only : caseid, ctitle, version, fsurdat
    use time_manager, only : get_nstep
    use fileutils, only : set_filename, putfil
    use clm_varctl, only : archive_dir, mss_wpass, mss_irt
    implicit none
    include 'netcdf.inc'

! ------------------------ local variables ---------------------------
    integer :: dim1_id(1)                      !netCDF dimension id for 1-d variables   
    integer :: dim2_id(2)                      !netCDF dimension id for 2-d variables
    integer :: dim3_id(3)                      !netCDF dimension id for 3-d variables
    integer :: omode                           !netCDF dummy variable
    integer :: status                          !netCDF error status
    character(len=256) :: loc_fn               !local 
    character(len=256) :: rem_dir              !remote (archive) directory
    character(len=256) :: rem_fn               !remote (archive) filename
    character(len=256) :: str                  !global attribute string 
    character(len=256) :: name                 !name of attribute
    character(len=256) :: unit                 !units of attribute
    character(len=256) :: mode                 !field mode (aver, inst, max, min, etc)
    character(len=  8) :: curdate              !current date
    character(len=  8) :: curtime              !current time 
    integer :: snlsno_id                       !netCDF variable id
    integer :: dzsno_id                        !netCDF variable id 
    integer :: zsno_id                         !netCDF variable id
    integer :: zisno_id                        !netCDF variable id
    integer :: h2osoi_liq_id                   !netCDF variable id 
    integer :: h2osoi_ice_id                   !netCDF variable id
    integer :: t_soisno_id                     !netCDF variable id  
    integer :: t_lake_id                       !netCDF variable id  
    integer :: t_veg_id                        !netCDF variable id
    integer :: t_grnd_id                       !netCDF variable id
    integer :: h2ocan_id                       !netCDF variable id
    integer :: h2osno_id                       !netCDF variable id  
    integer :: snowdp_id                       !netCDF variable id
    integer :: snowage_id                      !netCDF variable id    
#if (defined RTM)
    integer :: volr_id                         !netCDF variable id
#endif
#if (defined DGVM)
    integer :: itypveg_id                      !netCDF variable id 
    integer :: fpcgrid_id                      !netCDF variable id
    integer :: lai_ind_id                      !netCDF variable id
    integer :: crownarea_id                    !netCDF variable id
    integer :: litterag_id                     !netCDF variable id
    integer :: litterbg_id                     !netCDF variable id
    integer :: cpool_fast_id                   !netCDF variable id
    integer :: cpool_slow_id                   !netCDF variable id
    integer :: present_id                      !netCDF variable id 
    integer :: nind_id                         !netCDF variable id
    integer :: lm_ind_id                       !netCDF variable id
    integer :: sm_ind_id                       !netCDF variable id
    integer :: hm_ind_id                       !netCDF variable id
    integer :: rm_ind_id                       !netCDF variable id 
    integer :: htop_id                         !netCDF variable id
    integer :: wtxy_id                         !netCDF variable id  
#endif
    integer , allocatable :: ibuf1dl(:,:)
    integer , allocatable :: ibuf1dp(:)
    real(r8), allocatable :: rbuf1dl(:,:)
    real(r8), allocatable :: rbuf1dp(:)
    real(r8), allocatable :: rbuf2dl(:,:,:)
    real(r8), allocatable :: rbuf2dp(:,:)
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! create initial conditions file for writing
! --------------------------------------------------------------------

    if (masterproc) then

       loc_fn = set_init_filename()
       write(6,*)
       write(6,*)'(INICFILEMOD): Writing clm2 initial conditions dataset at ',&
            trim(loc_fn), 'at nstep = ',get_nstep()
       write(6,*)

! create new netCDF file (in defined mode)

       call wrap_create (trim(loc_fn), nf_clobber, ncid)

! set fill mode to "no fill" to optimize performance

       status = nf_set_fill (ncid, nf_nofill, omode)
       if (status /= nf_noerr) then
          write (6,*) ' netCDF error = ',nf_strerror(status)
          call endrun
       end if

! define dimensions 

       call wrap_def_dim (ncid, 'numland'  , numland        ,numland_dim)
       call wrap_def_dim (ncid, 'maxpatch' , maxpatch       ,maxpatch_dim)
       call wrap_def_dim (ncid, 'nlevsno'  , nlevsno        ,nlevsno_dim)
       call wrap_def_dim (ncid, 'nlevsoi'  , nlevsoi        ,nlevsoi_dim)
       call wrap_def_dim (ncid, 'nlevtot'  , nlevsno+nlevsoi,nlevtot_dim)
#if (defined RTM)
       call wrap_def_dim (ncid, 'rtmlon'   , rtmlon         ,rtmlon_dim)
       call wrap_def_dim (ncid, 'rtmlat'   , rtmlat         ,rtmlat_dim)
#endif

! define global attributes

       str = 'CF1.0'
       call wrap_put_att_text (ncid, NF_GLOBAL, 'conventions', trim(str))
       
       call datetime(curdate, curtime)
       str = 'created on ' // curdate // ' ' // curtime
       call wrap_put_att_text (ncid, NF_GLOBAL,'history', trim(str))
       
       call getenv ('LOGNAME', str)
       call wrap_put_att_text (ncid, NF_GLOBAL, 'logname',trim(str))
       
       call getenv ('HOST', str)
       call wrap_put_att_text (ncid, NF_GLOBAL, 'host', trim(str))
       
       str = 'Community Land Model: CLM2'
       call wrap_put_att_text (ncid, NF_GLOBAL, 'source', trim(str))
       
       str = '$Name$' 
       call wrap_put_att_text (ncid, NF_GLOBAL, 'version', trim(str))
       
       str = '$Id$'
       call wrap_put_att_text (ncid, NF_GLOBAL, 'revision_id', trim(str))
       
       str = ctitle 
       call wrap_put_att_text (ncid,nf_global,'case_title',trim(str))

       str = caseid
       call wrap_put_att_text (ncid,nf_global,'case_id',trim(str))

! define current date

       mode = 'instantaneous initial conditions'

       name = 'current date as 8 digit integer (YYYYMMDD)'
       call wrap_def_var (ncid, 'mcdate', nf_int, 0, 0, varid)
       call wrap_put_att_text (ncid, varid, 'long_name',name)
       call wrap_put_att_text (ncid, varid, 'mode'     ,mode)

       name = 'current seconds of current date'
       unit = 's'
       call wrap_def_var (ncid, 'mcsec',  nf_int, 0, 0, varid)
       call wrap_put_att_text (ncid, varid, 'long_name',name)
       call wrap_put_att_text (ncid, varid, 'units'    ,unit)
       call wrap_put_att_text (ncid, varid, 'mode'     ,mode)

! define single-level fields (numland x maxpatch)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'vegetation temperature (T_VEG)'
       unit = 'K'
       call wrap_def_var (ncid, 'T_VEG_INI', nf_double, 2, dim2_id, t_veg_id)
       call wrap_put_att_text (ncid, t_veg_id, 'long_name',name)
       call wrap_put_att_text (ncid, t_veg_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, t_veg_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'ground temperature (T_GRND)'
       unit = 'K'
       call wrap_def_var (ncid, 'T_GRND_INI', nf_double, 2, dim2_id, t_grnd_id)
       call wrap_put_att_text (ncid, t_grnd_id, 'long_name',name)
       call wrap_put_att_text (ncid, t_grnd_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, t_grnd_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'canopy water (H2OCAN)'
       unit = 'kg/m2'
       call wrap_def_var (ncid, 'H2OCAN_INI', nf_double, 2, dim2_id, h2ocan_id)
       call wrap_put_att_text (ncid, h2ocan_id, 'long_name',name)
       call wrap_put_att_text (ncid, h2ocan_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, h2ocan_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'snow water (H2OSNO)'
       unit = 'kg/m2'
       call wrap_def_var (ncid, 'H2OSNO_INI', nf_double, 2, dim2_id, h2osno_id)
       call wrap_put_att_text (ncid, h2osno_id, 'long_name',name)
       call wrap_put_att_text (ncid, h2osno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, h2osno_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'snow depth (SNOWDP)'
       unit = 'm'
       call wrap_def_var (ncid, 'SNOWDP_INI', nf_double, 2, dim2_id, snowdp_id)
       call wrap_put_att_text (ncid, snowdp_id, 'long_name',name)
       call wrap_put_att_text (ncid, snowdp_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, snowdp_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'snow age (SNOWAGE)'
       call wrap_def_var (ncid, 'SNOWAGE_INI', nf_double, 2, dim2_id, snowage_id)
       call wrap_put_att_text (ncid, snowage_id, 'long_name',name)
       call wrap_put_att_text (ncid, snowage_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'number of snow layers (SNLSNO)'
       call wrap_def_var (ncid, 'SNLSNO_INI', nf_int, 2, dim2_id, snlsno_id)
       call wrap_put_att_text (ncid, snlsno_id, 'long_name',name)
       call wrap_put_att_text (ncid, snlsno_id, 'mode'     ,mode)

! define multi-level fields (numland x maxpatch x numlev)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevtot_dim
       name = 'soil-snow temperature'
       unit = 'K'
       call wrap_def_var (ncid, 'T_SOISNO_INI', nf_double, 3, dim3_id, t_soisno_id)
       call wrap_put_att_text (ncid, t_soisno_id, 'long_name',name)
       call wrap_put_att_text (ncid, t_soisno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, t_soisno_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevsoi_dim
       name = 'lake temperature'
       unit = 'K'
       call wrap_def_var (ncid, 'T_LAKE_INI', nf_double, 3, dim3_id, t_lake_id)
       call wrap_put_att_text (ncid, t_lake_id, 'long_name',name)
       call wrap_put_att_text (ncid, t_lake_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, t_lake_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevtot_dim
       name = 'liquid water'
       unit = 'kg/m2'
       call wrap_def_var (ncid, 'H2OSOI_LIQ_INI', nf_double, 3, dim3_id, h2osoi_liq_id)
       call wrap_put_att_text (ncid, h2osoi_liq_id, 'long_name',name)
       call wrap_put_att_text (ncid, h2osoi_liq_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, h2osoi_liq_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevtot_dim
       name = 'ice lens'
       unit = 'kg/m2'                                                       
       call wrap_def_var (ncid, 'H2OSOI_ICE_INI', nf_double, 3, dim3_id, h2osoi_ice_id)
       call wrap_put_att_text (ncid, h2osoi_ice_id, 'long_name',name)
       call wrap_put_att_text (ncid, h2osoi_ice_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, h2osoi_ice_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevsno_dim
       name = 'snow layer depth'
       unit = 'm'
       call wrap_def_var (ncid, 'ZSNO_INI', nf_double, 3, dim3_id, zsno_id)
       call wrap_put_att_text (ncid, zsno_id, 'long_name',name)
       call wrap_put_att_text (ncid, zsno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, zsno_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevsno_dim
       name = 'snow layer thickness'
       unit = 'm'
       call wrap_def_var (ncid, 'DZSNO_INI', nf_double, 3, dim3_id, dzsno_id)
       call wrap_put_att_text (ncid, dzsno_id, 'long_name',name)
       call wrap_put_att_text (ncid, dzsno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, dzsno_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevsno_dim
       name = 'snow interface depth'
       unit = 'm'
       call wrap_def_var (ncid, 'ZISNO_INI', nf_double, 3, dim3_id, zisno_id)
       call wrap_put_att_text (ncid, zisno_id, 'long_name',name)
       call wrap_put_att_text (ncid, zisno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, zisno_id, 'mode'     ,mode)

#if (defined RTM)
       dim2_id(1) = rtmlon_dim ; dim2_id(2) = rtmlat_dim

       name = 'water volumn in cell (volr)'
       unit = 'm3'
       call wrap_def_var (ncid, 'RTMVOLR', nf_double, 2, dim2_id, volr_id)
       call wrap_put_att_text (ncid, volr_id, 'long_name',name)
       call wrap_put_att_text (ncid, volr_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, volr_id, 'mode'     ,mode)
#endif

#if (defined DGVM)
       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim

       name = 'plant functional type' 
       call wrap_def_var (ncid, 'ITYPVEG_INI', nf_int, 2, dim2_id, itypveg_id)
       call wrap_put_att_text (ncid, itypveg_id, 'long_name',name)
       call wrap_put_att_text (ncid, itypveg_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, itypveg_id, 'mode'     ,mode)

       name = 'plant functional type cover - fraction of vegetated area'
       call wrap_def_var (ncid, 'FPCGRID_INI', nf_double, 2, dim2_id, fpcgrid_id)
       call wrap_put_att_text (ncid, fpcgrid_id, 'long_name',name)
       call wrap_put_att_text (ncid, fpcgrid_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, fpcgrid_id, 'mode'     ,mode)

       name = 'LAI per individual'
       unit = 'm2/m2'
       call wrap_def_var (ncid, 'LAI_IND_INI', nf_double, 2, dim2_id, lai_ind_id)
       call wrap_put_att_text (ncid, lai_ind_id, 'long_name',name)
       call wrap_put_att_text (ncid, lai_ind_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, lai_ind_id, 'mode'     ,mode)

       name = 'crown area'
       unit = 'm2'
       call wrap_def_var (ncid, 'CROWNAREA_INI', nf_double, 2, dim2_id, crownarea_id)
       call wrap_put_att_text (ncid, crownarea_id, 'long_name',name)
       call wrap_put_att_text (ncid, crownarea_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, crownarea_id, 'mode'     ,mode)

       name = 'above ground litter'
       unit = 'grams C/m2 gridcell area'
       call wrap_def_var (ncid, 'LITTERAG_INI', nf_double, 2, dim2_id, litterag_id)
       call wrap_put_att_text (ncid, litterag_id, 'long_name',name)
       call wrap_put_att_text (ncid, litterag_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, litterag_id, 'mode'     ,mode)
       
       name = 'below ground litter'
       unit = 'grams C/m2 gridcell area'
       call wrap_def_var (ncid, 'LITTERBG_INI', nf_double, 2, dim2_id, litterbg_id)
       call wrap_put_att_text (ncid, litterbg_id, 'long_name',name)
       call wrap_put_att_text (ncid, litterbg_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, litterbg_id, 'mode'     ,mode)

       name = 'fast carbon pool'
       unit = 'grams C/m2 gridcell area'
       call wrap_def_var (ncid, 'CPOOL_FAST_INI', nf_double, 2, dim2_id, cpool_fast_id)
       call wrap_put_att_text (ncid, cpool_fast_id, 'long_name',name)
       call wrap_put_att_text (ncid, cpool_fast_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, cpool_fast_id, 'mode'     ,mode)

       name = 'slow carbon pool'
       unit = 'grams C/m2 gridcell area'
       call wrap_def_var (ncid, 'CPOOL_SLOW_INI', nf_double, 2, dim2_id, cpool_slow_id)
       call wrap_put_att_text (ncid, cpool_slow_id, 'long_name',name)
       call wrap_put_att_text (ncid, cpool_slow_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, cpool_slow_id, 'mode'     ,mode)

       name = 'if pft is present is patch, 1=>yes,=>no'
       call wrap_def_var (ncid, 'PRESENT_INI', nf_int, 2, dim2_id, present_id)
       call wrap_put_att_text (ncid, present_id, 'long_name',name)
       call wrap_put_att_text (ncid, present_id, 'mode'     ,mode)

       name = 'number of individuals'
       unit = 'number of individuals/m2 gridcell area'
       call wrap_def_var (ncid, 'NIND_INI', nf_double, 2, dim2_id, nind_id)
       call wrap_put_att_text (ncid, nind_id, 'long_name',name)
       call wrap_put_att_text (ncid, nind_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, nind_id, 'mode'     ,mode)

       name = 'individual leaf mass'
       unit = 'grams C'
       call wrap_def_var (ncid, 'LM_IND_INI', nf_double, 2, dim2_id, lm_ind_id)
       call wrap_put_att_text (ncid, lm_ind_id, 'long_name',name)
       call wrap_put_att_text (ncid, lm_ind_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, lm_ind_id, 'mode'     ,mode)

       name = 'individual stem mass'
       unit = 'grams C'
       call wrap_def_var (ncid, 'SM_IND_INI', nf_double, 2, dim2_id, sm_ind_id)
       call wrap_put_att_text (ncid, sm_ind_id, 'long_name',name)
       call wrap_put_att_text (ncid, sm_ind_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, sm_ind_id, 'mode'     ,mode)

       name = 'individual heartwood mass'
       unit = 'grams C'
       call wrap_def_var (ncid, 'HM_IND_INI', nf_double, 2, dim2_id, hm_ind_id)
       call wrap_put_att_text (ncid, hm_ind_id, 'long_name',name)
       call wrap_put_att_text (ncid, hm_ind_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, hm_ind_id, 'mode'     ,mode)

       name = 'individual root mass'
       unit = 'grams C'
       call wrap_def_var (ncid, 'RM_IND_INI', nf_double, 2, dim2_id, rm_ind_id)
       call wrap_put_att_text (ncid, rm_ind_id, 'long_name',name)
       call wrap_put_att_text (ncid, rm_ind_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, rm_ind_id, 'mode'     ,mode)

       name = 'canopy top'
       unit = 'm'
       call wrap_def_var (ncid, 'HTOP_INI', nf_double, 2, dim2_id, htop_id)
       call wrap_put_att_text (ncid, htop_id, 'long_name',name)
       call wrap_put_att_text (ncid, htop_id, 'unit'     ,unit)
       call wrap_put_att_text (ncid, htop_id, 'mode'     ,mode)

       name = 'subgrid weight'
       call wrap_def_var (ncid, 'WTXY_INI', nf_double, 2, dim2_id, wtxy_id)
       call wrap_put_att_text (ncid, wtxy_id, 'long_name',name)
       call wrap_put_att_text (ncid, wtxy_id, 'mode'     ,mode)

#endif
! finish creating netCDF file (end define mode)

       status = nf_enddef(ncid)

    endif  !end of if-masterproc block

! --------------------------------------------------------------------
! Write single-level variables: [numland x maxpatch] and 
! multi-level variables: [numland x maxpatch x lev]
! NOTE: for non-spmd begpatch=1 and endpatch=numpatch
! --------------------------------------------------------------------

! Convert clm derived type components to patch vectors
! Convert initial data from subgrid patch form to land point form

    allocate (ibuf1dl(numland,maxpatch))
    allocate (rbuf1dl(numland,maxpatch))
    allocate (ibuf1dp(begpatch:endpatch))
    allocate (rbuf1dp(begpatch:endpatch))

    ! Write out zisno
    allocate (rbuf2dp(-nlevsno+0:-1,begpatch:endpatch))
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+0:-1)) 
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+0:-1,k) = clm(k)%zi(-nlevsno+0:-1)  
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno)
    if (masterproc) call wrap_put_var_realx (ncid, zisno_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out zsno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:0))
    allocate (rbuf2dp(-nlevsno+1: 0,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1:0,k) = clm(k)%z(-nlevsno+1:0)  
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno)
    if (masterproc) call wrap_put_var_realx (ncid, zsno_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out dzsno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:0))
    allocate (rbuf2dp(-nlevsno+1: 0,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1: 0,k) = clm(k)%dz(-nlevsno+1: 0)  
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno)
    if (masterproc) call wrap_put_var_realx (ncid, dzsno_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out h2osoi_liq
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1:nlevsoi,k) = clm(k)%h2osoi_liq(-nlevsno+1:nlevsoi) 
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno+nlevsoi)
    if (masterproc) call wrap_put_var_realx (ncid, h2osoi_liq_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out h2osoi_ice
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1:nlevsoi,k) = clm(k)%h2osoi_ice(-nlevsno+1:nlevsoi) 
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno+nlevsoi)
    if (masterproc) call wrap_put_var_realx (ncid, h2osoi_ice_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out t_soisno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1:nlevsoi,k) = clm(k)%t_soisno(-nlevsno+1:nlevsoi) 
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno+nlevsoi)
    if (masterproc) call wrap_put_var_realx (ncid, t_soisno_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out t_lake
    allocate(rbuf2dl(numland,maxpatch,1:nlevlak))
    allocate(rbuf2dp(1:nlevlak,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(1:nlevlak,k) = clm(k)%t_lake(1:nlevlak)
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevlak)
    if (masterproc) call wrap_put_var_realx (ncid, t_lake_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out t_veg
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%t_veg   
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, t_veg_id, rbuf1dl)

    ! Write out t_grnd
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%t_grnd  
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, t_grnd_id, rbuf1dl)

    ! Write out h2ocan
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%h2ocan  
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, h2ocan_id, rbuf1dl)

    ! Write out h2osno
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%h2osno  
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, h2osno_id, rbuf1dl)

    ! Write out snowdp
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%snowdp
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, snowdp_id, rbuf1dl)

    ! Write out snowage
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%snowage
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, snowage_id, rbuf1dl)

    ! Write out snlsno
    do k = begpatch,endpatch
       ibuf1dp(k) = clm(k)%snl     
    end do
    call patch_to_land (ibuf1dp, ibuf1dl)
    if (masterproc) call wrap_put_var_int (ncid, snlsno_id, ibuf1dl)

#if (defined RTM)
    ! Write out volr
    if (masterproc) call wrap_put_var_realx (ncid, volr_id, volr)
#endif

#if (defined DGVM)
    ! Write out itypveg    
    do k = begpatch,endpatch
       ibuf1dp(k) = clm(k)%itypveg
    end do
    call patch_to_land (ibuf1dp, ibuf1dl)
    if (masterproc) call wrap_put_var_int (ncid, itypveg_id, ibuf1dl)

    ! write out fpcgrid
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%fpcgrid
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, fpcgrid_id, rbuf1dl)

    ! write out lai_ind
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%lai_ind
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, lai_ind_id, rbuf1dl)

    ! write out crownarea
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%crownarea
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, crownarea_id, rbuf1dl)

    ! write out litterag
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%litterag
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, litterag_id, rbuf1dl)

    ! write out litterbg
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%litterbg
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, litterbg_id, rbuf1dl)

    ! write out cpool_fast
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%cpool_fast
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, cpool_fast_id, rbuf1dl)

    ! write out cpool_slow
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%cpool_slow
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, cpool_slow_id, rbuf1dl)

    ! write out present
    do k = begpatch,endpatch
       ibuf1dp(k) = 0
       if (clm(k)%present) ibuf1dp(k) = 1
    end do
    call patch_to_land (ibuf1dp, ibuf1dl)
    if (masterproc) call wrap_put_var_int (ncid, present_id, ibuf1dl)

    ! write out nind
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%nind
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, nind_id, rbuf1dl)

    ! write out lm_ind
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%lm_ind
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, lm_ind_id, rbuf1dl)

    ! write out sm_ind
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%sm_ind
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, sm_ind_id, rbuf1dl)

    ! write out hm_ind
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%hm_ind
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, hm_ind_id, rbuf1dl)

    ! write out rm_ind
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%rm_ind
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, rm_ind_id, rbuf1dl)

    ! write out htop
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%htop
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, htop_id, rbuf1dl)

    ! write out wtxy (master processor already has all the weights - even on restart)
    ! in this case only want to write out valid values for weights
    ! do not want to use in "filled" in values for those patches with zero weights
    ! whose value is set to the first patch with non-zero weight (see routine clm_map)
    if (masterproc) then
       rbuf1dl(:,:) = 0._r8
       do k = 1,numpatch
          l = patchvec%land(k)
          m = patchvec%mxy(k)
          rbuf1dl(l,m) = patchvec%wtxy(k)
       end do
       call wrap_put_var_realx (ncid, wtxy_id, rbuf1dl)
    endif

#endif

    deallocate (ibuf1dl)
    deallocate (rbuf1dl)
    deallocate (ibuf1dp)
    deallocate (rbuf1dp)

! archive initial conditions dataset (Mass Store currently)

    if (masterproc) then
       call wrap_close (ncid)
       if (mss_irt > 0) then 
          rem_dir = trim(archive_dir) // '/init/'
          rem_fn = set_filename(rem_dir, loc_fn)
          call putfil (loc_fn, rem_fn, mss_wpass, mss_irt, .true.)
       endif
    endif

    return
  end subroutine inicwrt

!=======================================================================
! BEGIN GENERIC PROCEDURE DEFINITIONS
!=======================================================================

  logical function do_inicwrite()

    use time_manager, only : get_curr_date, get_prev_date
    use clm_varctl, only : hist_crtinic

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Determine if initial dataset is to be written at this time step
!
! Method: 
!
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

! ------------------------ local variables ------------------------
    integer :: yr         !nstep year (0 -> ...)
    integer :: yrm1       !nstep-1 year (0 -> ...)
    integer :: daym1      !nstep-1 day (1 -> 31)
    integer :: day        !nstep day (1 -> 31)
    integer :: mon        !nstep month (1 -> 12)
    integer :: monm1      !nstep-1 month (1 -> 12)
    integer :: mcsec      !nstep time of day [seconds] 
    integer :: mcsecm1    !nstep-1 time of day [seconds]
! -----------------------------------------------------------------

    ! Set calendar for current time step and previous time step

    call get_curr_date (yr, mon, day, mcsec) 
    call get_prev_date (yrm1, monm1, daym1, mcsecm1)

    ! Determine if time to write out initial dataset

    do_inicwrite = .false.
    if (hist_crtinic /= 'NONE') then
       if (hist_crtinic == 'MONTHLY') then
          if (mon /= monm1 .and. monm1 /= -1) do_inicwrite = .true.
       else if (hist_crtinic == 'YEARLY') then
          if (monm1 == 12 .and. mon == 1)  do_inicwrite = .true.
       endif
    endif

  end function do_inicwrite

!=======================================================================

  character(len=256) function set_init_filename ()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Determine initial dataset filenames
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use clm_varctl  , only : caseid
    use time_manager, only : get_curr_date

! ------------------------ local variables ------------------------
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
! -----------------------------------------------------------------

    call get_curr_date (yr, mon, day, sec) 
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_init_filename = "./"//trim(caseid)//".clm2.i."//trim(cdate)//".nc"

  end function set_init_filename

!=======================================================================

!----------------------------------------------------------------------- 
! 
! Purpose: 
! [numland] x [maxpatch] array from 1d subgrid patches
!
! Method: 
! Map a subgrid input vector [fldin] of length [numpatch] to a 2-d
! [numland] x [maxpatch] output array [fldout]. Not all land points have
! [maxpatch] subgrid patches. Many have less. [numpatch] is some number <=
! [numland]*[maxpatch], i.e., is valid subgrid patches only. This routine
! converts a field from its [numpatch] representation to a [numland] x 
! [maxpatch] representation, setting values for non-valid subgrid patches 
! to that of the first valid subgrid patch using mapping from clm_map
! o numland  = number of land grid cells
! o maxpatch = maximum number of subgrid patches per grid cell
! o numpatch = actual number of subgrid patches (<= numland*maxpatch)
! 
!-----------------------------------------------------------------------

  subroutine patch_to_land_1d_int (fldin, fldout)
! ------------------------ arguments ---------------------------------
    integer, intent(in)  :: fldin(begpatch:endpatch)          
    integer, intent(out) :: fldout(numland,maxpatch) 
! --------------------------------------------------------------------

! ------------------------ local variables ----------------------
    integer l,m,k                   !indices
#if (defined SPMD)
    integer :: ier                   !MPI error status
    integer :: ibuf1d(numpatch)      !MPI temporary buffer 
    integer :: numrecvv(0:npes-1)    !vector of items to be received  
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numsend               !number of items to be sent
! ---------------------------------------------------------------
    ibuf1d(begpatch:endpatch) = fldin(begpatch:endpatch)
    call compute_mpigs_patch(1, numsend, numrecvv, displsv)
    call mpi_gatherv (ibuf1d(begpatch), numsend , mpiint,  &
                      ibuf1d          , numrecvv, displsv, mpiint, 0, mpicom, ier)
    if (masterproc) then
       do m = 1, maxpatch           !subgrid patches for each land point
          do l = 1, numland         !land point index for [lsmlon] x [lsmlat] grid
             k = landvec%patch(l,m) !subgrid patch index: [1] to [numpatch]
             fldout(l,m) = ibuf1d(k)
          end do
       end do
    endif
#else 
    do m = 1, maxpatch              !subgrid patches for each land point
       do l = 1, numland            !land point index for [lsmlon] x [lsmlat] grid
          k = landvec%patch(l,m)    !subgrid patch index: [1] to [numpatch]
          fldout(l,m) = fldin(k)
       end do
    end do
#endif
    return
  end subroutine patch_to_land_1d_int

!=======================================================================

  subroutine patch_to_land_1d_real (fldin, fldout)
! ------------------------ arguments ---------------------------------
    real(r8), intent(in)  :: fldin(begpatch:endpatch)          
    real(r8), intent(out) :: fldout(numland,maxpatch) 
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k                   !indices
#if (defined SPMD)
    integer  :: ier                 !MPI error status
    real(r8) :: buf1d(numpatch)     !MPI temporary buffer 
    integer  :: numrecvv(0:npes-1)  !vector of items to be received  
    integer  :: displsv(0:npes-1)   !displacement vector
    integer  :: numsend             !number of items to be sent
! ---------------------------------------------------------------
    buf1d(begpatch:endpatch) = fldin(begpatch:endpatch)
    call compute_mpigs_patch(1, numsend, numrecvv, displsv)
    call mpi_gatherv (buf1d(begpatch), numsend , mpir8, &
                      buf1d          , numrecvv, displsv, mpir8  , 0, mpicom, ier)
    if (masterproc) then
       do m = 1, maxpatch           !subgrid patches for each land point
          do l = 1, numland         !land point index for [lsmlon] x [lsmlat] grid
             k = landvec%patch(l,m) !subgrid patch index: [1] to [numpatch]
             fldout(l,m) = buf1d(k)
          end do
       end do
    endif
#else 
    do m = 1, maxpatch               !subgrid patches for each land point
       do l = 1, numland             !land point index for [lsmlon] x [lsmlat] grid
          k = landvec%patch(l,m)     !subgrid patch index: [1] to [numpatch]
          fldout(l,m) = fldin(k)
       end do
    end do
#endif
    return
  end subroutine patch_to_land_1d_real

!=======================================================================

  subroutine patch_to_land_2d_real (fldin, fldout, nlev)
! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: nlev
    real(r8), intent(in)  :: fldin(nlev,begpatch:endpatch)
    real(r8), intent(out) :: fldout(numland,maxpatch,nlev) 
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k,n                    !indices
#if (defined SPMD)                    
    integer  :: ier                   !MPI error status
    real(r8) :: buf2d(nlev,numpatch)  !MPI temporary buffer
    integer  :: numrecvv(0:npes-1)    !vector of items to be received  
    integer  :: displsv(0:npes-1)     !displacement vector
    integer  :: numsend               !number of items to be sent
! ---------------------------------------------------------------
    do k=begpatch,endpatch
       do n=1,nlev
          buf2d(n,k) = fldin(n,k)
       end do
    end do
    call compute_mpigs_patch(nlev, numsend, numrecvv, displsv)
    call mpi_gatherv (buf2d(1,begpatch), numsend , mpir8, &
                      buf2d            , numrecvv, displsv, mpir8, 0, mpicom, ier)
    if (masterproc) then
       do m = 1, maxpatch           !subgrid patches for each land point
          do l = 1, numland         !land point index for [lon] x [lat] grid
             k = landvec%patch(l,m) !subgrid patch index: [1] to [numpatch]
             do n = 1,nlev          !level index 
                fldout(l,m,n) = buf2d(n,k)
             end do
          end do
       end do
    endif
#else
    do m = 1,maxpatch               !subgrid patches for each land point
       do l = 1,numland             !land point index for [lon] x [lat] grid
          k = landvec%patch(l,m)    !subgrid patch index: [1] to [numpatch]
          do n = 1,nlev             !level index 
             fldout(l,m,n) = fldin(n,k)
          end do
       end do
    end do
#endif
    return
  end subroutine patch_to_land_2d_real

!=======================================================================

  subroutine land_to_patch_1d_int (fldin, fldout)
! ------------------------ arguments ---------------------------------
    integer, intent(in)  :: fldin(numland,maxpatch)
    integer, intent(out) :: fldout(begpatch:endpatch)
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k                   !indices
#if (defined SPMD)
    integer :: ier                  !MPI error status
    integer :: ibuf1d(numpatch)     !MPI temporary buffer
    integer :: numsendv(0:npes-1)   !vector of items to be sent
    integer :: displsv(0:npes-1)    !displacement vector
    integer :: numrecv              !number of items to be received
! ---------------------------------------------------------------
    if (masterproc) then
       do m = 1, maxpatch              !subgrid patches for each land point
          do l = 1, numland            !land point index 
             if (landvec%wtxy(l,m) > 0.) then
                k = landvec%patch(l,m) !subgrid patch index
                ibuf1d(k) = fldin(l,m)
             end if
          end do
       end do
    endif
    call compute_mpigs_patch(1, numrecv, numsendv, displsv)
    call mpi_scatterv (ibuf1d          , numsendv, displsv, mpiint, &
                       ibuf1d(begpatch), numrecv , mpiint, 0, mpicom, ier)
    fldout(begpatch:endpatch) = ibuf1d(begpatch:endpatch)
#else
    do m = 1, maxpatch                !subgrid patches for each land point
       do l = 1, numland              !land point index 
          if (landvec%wtxy(l,m) > 0.) then
             k = landvec%patch(l,m)   !subgrid patch index
             fldout(k) = fldin(l,m)
          endif
       end do
    end do
#endif
    return
  end subroutine land_to_patch_1d_int

!=======================================================================

  subroutine land_to_patch_1d_real (fldin, fldout)
! ------------------------ arguments ---------------------------------
    real(r8), intent(in)  :: fldin(numland,maxpatch)
    real(r8), intent(out) :: fldout(begpatch:endpatch)
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k  !indices
#if (defined SPMD)
    integer  :: ier                  !MPI error status
    real(r8) :: buf1d(numpatch)      !MPI temporary buffer
    integer  :: numsendv(0:npes-1)   !vector of items to be sent
    integer  :: displsv(0:npes-1)    !displacement vector
    integer  :: numrecv              !number of items to be received
! ---------------------------------------------------------------
    if (masterproc) then
       do m = 1, maxpatch             !subgrid patches for each land point
          do l = 1, numland           !land point index 
             if (landvec%wtxy(l,m) > 0.) then
                k = landvec%patch(l,m) !subgrid patch index
                buf1d(k) = fldin(l,m)
             end if
          end do
       end do
    endif
    call compute_mpigs_patch(1, numrecv, numsendv, displsv)
    call mpi_scatterv (buf1d          , numsendv, displsv, mpir8, &
                       buf1d(begpatch), numrecv , mpir8  , 0, mpicom, ier)
    fldout(begpatch:endpatch) = buf1d(begpatch:endpatch)
#else
    do m = 1, maxpatch              !subgrid patches for each land point
       do l = 1, numland            !land point index 
          if (landvec%wtxy(l,m) > 0.) then
             k = landvec%patch(l,m) !subgrid patch index
             fldout(k) = fldin(l,m)
          endif
       end do
    end do
#endif
    return
  end subroutine land_to_patch_1d_real

!=======================================================================

  subroutine land_to_patch_2d_real (fldin, fldout, nlev)
! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: nlev
    real(r8), intent(in)  :: fldin (numland,maxpatch,nlev) 
    real(r8), intent(out) :: fldout(nlev,begpatch:endpatch)
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k,n                   !indices
#if (defined SPMD)                   
    integer  :: ier                   !MPI error status
    real(r8) :: buf2d(nlev,numpatch)  !MPI temporary buffer
    integer  :: numsendv(0:npes-1)    !vector of items to be sent
    integer  :: displsv(0:npes-1)     !displacement vector
    integer  :: numrecv               !number of items to be received
! ---------------------------------------------------------------
    if (masterproc) then
       do m = 1, maxpatch             !subgrid patches for each land point
          do l = 1, numland           !land point index 
             if (landvec%wtxy(l,m) > 0.) then
                k = landvec%patch(l,m) !subgrid patch index
                do n = 1,nlev
                   buf2d(n,k) = fldin(l,m,n)
                end do
             end if
          end do
       end do
    endif
    call compute_mpigs_patch(nlev, numrecv, numsendv, displsv)
    call mpi_scatterv (buf2d            , numsendv, displsv, mpir8, &
                       buf2d(1,begpatch), numrecv , mpir8  , 0, mpicom, ier)
    fldout(1:nlev,begpatch:endpatch) = buf2d(1:nlev,begpatch:endpatch)
#else
    do m = 1, maxpatch               !subgrid patches for each land point
       do l = 1, numland             !land point index 
          k = landvec%patch(l,m)  !subgrid patch index
          if (landvec%wtxy(l,m) > 0.) then
             do n = 1,nlev
                fldout(n,k) = fldin(l,m,n)
             end do
          endif
       end do
    end do
#endif
    return
  end subroutine land_to_patch_2d_real

!=======================================================================

  subroutine hist_to_clm_1d(fldin)
! ------------------------ arguments ---------------------------------
    real(r8), intent(inout)  :: fldin(numpatch)
! --------------------------------------------------------------------
#if (defined SPMD)
! ------------------------ local variables ----------------------
    integer  :: ier                  !MPI error status
    real(r8) :: buf1d(numpatch)      !MPI temporary buffer
    integer  :: numsendv(0:npes-1)    !vector of items to be sent
    integer  :: displsv(0:npes-1)     !displacement vector
    integer  :: numrecv               !number of items to be received
! ---------------------------------------------------------------
    if (masterproc) buf1d(1:numpatch) = fldin(1:numpatch)
    call compute_mpigs_patch(1, numrecv, numsendv, displsv)
    call mpi_scatterv (buf1d          , numsendv, displsv, mpir8, &
                       buf1d(begpatch), numrecv , mpir8, 0, mpicom, ier)
    fldin(begpatch:endpatch) = buf1d(begpatch:endpatch)
#endif
    return
  end subroutine hist_to_clm_1d

!=======================================================================

  subroutine hist_to_clm_2d(fldin, nlev)
! ------------------------ arguments ---------------------------------
    integer , intent(in)    :: nlev
    real(r8), intent(inout) :: fldin(numpatch,nlev)
! --------------------------------------------------------------------
#if (defined SPMD)
! ------------------------ local variables ----------------------
    integer  :: j                    !index  
    integer  :: ier                  !MPI error status
    real(r8) :: buf2d(nlev,numpatch) !MPI temporary buffer
    integer  :: numsendv(0:npes-1)   !vector of items to be sent
    integer  :: displsv(0:npes-1)    !displacement vector
    integer  :: numrecv              !number of items to be received
! ---------------------------------------------------------------
    if (masterproc) then
       do j = 1,nlev
          buf2d(j,1:numpatch) = fldin(1:numpatch,j)
       end do
    end if
    call compute_mpigs_patch(nlev, numrecv, numsendv, displsv)
    call mpi_scatterv (buf2d            , numsendv, displsv, mpir8, &
                       buf2d(1,begpatch), numrecv , mpir8  , 0, mpicom, ier)
    do j = 1,nlev
       fldin(begpatch:endpatch,j) = buf2d(j,begpatch:endpatch)
    end do
#endif
    return
  end subroutine hist_to_clm_2d

!=======================================================================
! END GENERIC PROCEDURE DEFINITIONS
!=======================================================================

end module inicFileMod
