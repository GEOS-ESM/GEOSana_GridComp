#include <misc.h>
#include <preproc.h>

subroutine lpj ()
 
  use precision
  use clmtype
  use clm_varder
  use clm_varpar, only : maxpatch_pft, lsmlon , lsmlat
  use clm_varmap, only : numpatch, patchvec, begpatch, endpatch, &
                         numland, landvec, begland, endland
  use pft_varcon, only : pftpar, tree, sla
  use time_manager, only : get_curr_date, get_ref_date
  use histFileDGVM, only : ncid, beg3d, len3d, beg4d, len4d, &
                           afirefrac_id, acfluxfire_id, bmfm_id, afmicr_id, &
	                   histcrt_dgvm, histwrt_dgvm, spval
#if (defined SPMD)
  use spmdMod, only : masterproc, npes, compute_mpigs_patch, compute_mpigs_land, iam
  use mpishorthand, only : mpir8, mpilog, mpicom
#else
  use spmdMod, only : masterproc
#endif

! ----------------------------------------------------------------------
! purpose           : to drive the annual portion of lpj; called 1 per year
! date first created: October 2000 - lsm version 2 
! by whom           : Sam Levis
! date last revised : 
! by whom           : 
! ----------------------------------------------------------------------
! $Id$
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
  integer  :: i,j,k,l,m,ii,ji    !indices
  integer  :: yr                 !year (0 -> ...)
  integer  :: mon                !month (1 -> 12)
  integer  :: day                !day (1 -> 31)
  integer  :: sec                !seconds 
  integer  :: kyr                !used in routine climate20' below
  integer  :: ncdate             !current date
  integer  :: nbdate             !base date (reference date)
  real(r8) :: wscal
  real(r8) :: mort_max
  real(r8) :: turnover_ind(begpatch:endpatch)
  real(r8) :: acfluxfire_patch
  real(r8) :: afirefrac_patch
  real(r8) :: afirefrac_land(numland)
  real(r8) :: acfluxfire_land(numland)
  real(r8) :: bminc_patch(numpatch)
  real(r8) :: bmfm_land(maxpatch_pft,numland)
  real(r8) :: afmicr_land(maxpatch_pft,numland)
  real(r8) :: rbuf2d_xy(lsmlon,lsmlat)              ! for history write
  real(r8) :: rbuf3d_xy(lsmlon,lsmlat,maxpatch_pft) ! for history write
#if (defined SPMD)
  integer :: numrecvv(0:npes-1)   !vector of items to be received  
  integer :: displsv(0:npes-1)    !displacement vector
  integer :: numsend              !number of items to be sent
  integer :: ier                  !MPI error status
#endif
! ----------------------------------------------------------------------

! initialize land output variables

  bminc_patch(:) = 0.0 

  afirefrac_land(:) = 0.0 
  acfluxfire_land(:) = 0.0 
  bmfm_land(:,:) = 0.0
  afmicr_land(:,:) = 0.0

! kyr used in routine climate20' below
! NB: at end of first year, kyr = 2

  call get_curr_date (yr, mon, day, sec)  
  ncdate = yr*10000 + mon*100 + day
  call get_ref_date (yr, mon, day, sec)
  nbdate = yr*10000 + mon*100 + day
  kyr = ncdate/10000 - nbdate/10000 + 1   
  write(6,*)'LPJ: ncdate= ',ncdate,' nbdate= ',nbdate,' kyr= ',kyr

! loop over patches

  do k = begpatch,endpatch

! ******************************************************************************
! My version of LPJ's routine climate20 - 'Returns' tmomin20 and agdd20 for
! use in routine bioclim, which I have placed in routine Establishment
! Instead of 20-yr running mean of coldest monthly temperature,
! use 20-yr running mean of minimum 10-day running mean

     if (kyr == 2) clm(k)%tmomin20 = clm(k)%t_mo_min
     if (kyr == 2) clm(k)%agdd20 = clm(k)%agdd
     clm(k)%tmomin20 = (19.0 * clm(k)%tmomin20 + clm(k)%t_mo_min) / 20.0
     clm(k)%agdd20   = (19.0 * clm(k)%agdd20   + clm(k)%agdd    ) / 20.0

! ******************************************************************************

     bminc_patch(k) = clm(k)%bm_inc                 ![gC/m2 patch] for output
     clm(k)%bm_inc = clm(k)%bm_inc * clm(k)%fpcgrid ![gC/m2 cell vegetated area]

! determine grid values of npp and microbial respiration

     l = patchvec%land(k)         !land point index
     m = patchvec%mxy(k)          !patch index
     if (m <= maxpatch_pft) then
        bmfm_land(m,l) = bminc_patch(k)                   ![gC/m2 patch] for output 
        afmicr_land(m,l) = clm(k)%afmicr * clm(k)%fpcgrid ![gC/m2 cell veg'd area]
     endif

     if (clm(k)%itypveg > 0) then

! returns updated bm_inc, litterag

        call Reproduction (clm(k)%bm_inc, clm(k)%litterag, clm(k)%present)

! returns turnover_ind and updated litterag,bg, l,s,h,rm_ind

        call Turnover (pftpar(clm(k)%itypveg,:), clm(k)%litterag, clm(k)%litterbg, &
                       clm(k)%lm_ind           , clm(k)%sm_ind  , clm(k)%hm_ind  , &
                       clm(k)%rm_ind           , clm(k)%nind    , clm(k)%present , &
                       turnover_ind(k))

! returns updated litterag,bg, and present

        call Kill (clm(k)%bm_inc, clm(k)%litterag, clm(k)%litterbg, &
                   clm(k)%lm_ind, clm(k)%sm_ind  , clm(k)%hm_ind  , &
                   clm(k)%rm_ind, clm(k)%nind    , clm(k)%present , &
                   tree(clm(k)%itypveg))

! returns lai_ind, lai_inc, updates crownarea, htop, l,s,h,rm_ind, litterag,bg

        if (clm(k)%annpsnpot > 0.0) then
           wscal = clm(k)%annpsn/clm(k)%annpsnpot
        else
           wscal = 1.0
        end if

        call Allocation (pftpar(clm(k)%itypveg,:), clm(k)%nind, &
                         clm(k)%bm_inc  , tree(clm(k)%itypveg), &
                         clm(k)%htop    , sla(clm(k)%itypveg) , &
                         clm(k)%lm_ind  , clm(k)%sm_ind       , &
                         clm(k)%hm_ind  , clm(k)%rm_ind       , &
                         clm(k)%present , clm(k)%crownarea    , &
                         clm(k)%litterag, clm(k)%litterbg     , &
                         clm(k)%lai_ind , clm(k)%fpcgrid      , &
                         clm(k)%fpcinc  , wscal               )
     end if

  end do

! returns lm,rm_ind, fpcgrid, nind, litterag,bg via modules
! reason for different set up (ie, no external patch loop):
! in this routine sub-grid patches (k) communicate at the grid cell level (i,j)

  call Light ()

! loop over patches

  do k = begpatch,endpatch

! returns updated present, nind, litterag,bg

     if (clm(k)%itypveg > 0) then

        if (clm(k)%itypveg == 3 .or. clm(k)%itypveg == 8) then
           mort_max = 0.03 !testing diff values (add to Ecosystemini?)
        else
           mort_max = 0.01 !orig. value for all pfts
        end if
        call Mortality (clm(k)%bm_inc  , clm(k)%nind    , turnover_ind(k) , &
                        clm(k)%lm_ind  , clm(k)%sm_ind  , clm(k)%hm_ind   , &
                        clm(k)%rm_ind  , sla(clm(k)%itypveg), clm(k)%litterag , &
                        clm(k)%litterbg, clm(k)%present , tree(clm(k)%itypveg), &
                        clm(k)%agddtw  , mort_max       )

! returns updated litterag, nind

        call Fire (clm(k)%firelength   , clm(k)%litterag         , clm(k)%present, &
                   tree(clm(k)%itypveg), pftpar(clm(k)%itypveg,8), clm(k)%nind   , &
                   clm(k)%lm_ind       , clm(k)%sm_ind           , clm(k)%hm_ind , &
                   clm(k)%rm_ind       , afirefrac_patch         , clm(k)%fpcgrid, &
                   acfluxfire_patch    )

     else

        afirefrac_patch = 0.0
        acfluxfire_patch = 0.0

     end if

! determine land point value

     l = patchvec%land(k)         
     afirefrac_land(l) = afirefrac_land(l) + afirefrac_patch * clm(k)%fpcgrid
     acfluxfire_land(l) = acfluxfire_land(l) + acfluxfire_patch

  end do

! returns updated present, nind, *m_ind, crownarea, fpcgrid, htop, litter*g
! reason for different set up (ie, no external patch loop):
! in this routine sub-grid patches (k) communicate at the grid cell level (i,j)

  call Establishment ()

! -----------------------------------------------------------------------
! Write grid info, time info and first set of variables to  to DGVM
! history file
! -----------------------------------------------------------------------

! Create DGVM history file

  if (masterproc) call histcrt_dgvm()

! Write grid and time information to DGVM history file

  call histwrt_dgvm()	

! First write of first set of output fields to DGVM history file
! Rest of fields are written out in lpjreset2.F90
! First map land or patch points to xy arrays and then write to history file

! Write out afirefrac

#if (defined SPMD)
  call compute_mpigs_land(1, numsend, numrecvv, displsv)
  call mpi_gatherv (afirefrac_land(begland) , numsend , mpir8, &
                    afirefrac_land          , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
  if (masterproc) then
     rbuf2d_xy(:,:)= spval
     do l = 1, numland
        i = landvec%ixy(l)    
        j = landvec%jxy(l)    
        rbuf2d_xy(i,j) = afirefrac_land(l)
     end do
     call wrap_put_vara_realx (ncid, afirefrac_id, beg3d, len3d, rbuf2d_xy)
  endif

! Write out acfluxfire
  
#if (defined SPMD)
  call compute_mpigs_land(1, numsend, numrecvv, displsv)
  call mpi_gatherv (acfluxfire_land(begland), numsend , mpir8, &
                    acfluxfire_land         , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
  if (masterproc) then
     rbuf2d_xy(:,:)= spval
     do l = 1, numland
        i = landvec%ixy(l)    
        j = landvec%jxy(l)    
        rbuf2d_xy(i,j) = acfluxfire_land(l)
     end do
     call wrap_put_vara_realx (ncid, acfluxfire_id, beg3d, len3d, rbuf2d_xy)
  endif

! Write out bmfm

#if (defined SPMD)
  call compute_mpigs_land(maxpatch_pft    , numsend, numrecvv, displsv)
  call mpi_gatherv (bmfm_land(1,begland)  , numsend , mpir8, &
                    bmfm_land             , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
  if (masterproc) then
     rbuf3d_xy(:,:,:) = spval
     do l = 1, numland
        i = landvec%ixy(l)    
        j = landvec%jxy(l)    
        do m = 1, maxpatch_pft
           rbuf3d_xy(i,j,m) = bmfm_land(m,l)
        end do
     end do
     call wrap_put_vara_realx (ncid, bmfm_id, beg4d, len4d, rbuf3d_xy)
  endif  

! Write out afmicr

#if (defined SPMD)
  call compute_mpigs_land(maxpatch_pft    , numsend, numrecvv, displsv)
  call mpi_gatherv (afmicr_land(1,begland), numsend , mpir8, &
                    afmicr_land           , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif           
  if (masterproc) then
     rbuf3d_xy(:,:,:) = spval
     do l = 1, numland
        i = landvec%ixy(l)    
        j = landvec%jxy(l)    
        do m = 1, maxpatch_pft
           rbuf3d_xy(i,j,m) = afmicr_land(m,l)
        end do
     end do
     call wrap_put_vara_realx (ncid, afmicr_id, beg4d, len4d, rbuf3d_xy)
  endif  

  return
end subroutine lpj
