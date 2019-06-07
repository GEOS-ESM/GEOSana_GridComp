module geos_nstmod

! Santha: where are the prologues of these routines????
!
! 23Mar2015 Todling rid of nst_mask_full

use kinds,       only : r_kind

implicit none
private
public geos_nst_init
public geos_nst_set
public geos_nst_final

interface geos_nst_init
  module procedure init_
end interface
!----------------------------------------

interface geos_nst_set
  module procedure set_
end interface
!----------------------------------------

interface geos_nst_final
  module procedure final_
end interface
!----------------------------------------
 
contains

subroutine init_

use kinds,            only: i_kind, r_kind
use mpimod,           only: mype
use gsi_nstcouplermod,only: tref_full,dt_cool_full,dt_warm_full,z_c_full,z_w_full
use gsi_nstcouplermod,only: nstinfo
use guess_grids,      only: nfldsfc,nfldnst,ntguesnst  
use gridmod,          only: nlat,nlon
use mpeu_util,        only: die, perr

implicit none

! allocation of _full variables 
  if(.not.allocated(z_c_full))     allocate(z_c_full     (nlat,nlon,nfldsfc))
  if(.not.allocated(z_w_full))     allocate(z_w_full     (nlat,nlon,nfldsfc))
  if(.not.allocated(dt_cool_full)) allocate(dt_cool_full (nlat,nlon,nfldsfc))
  if(.not.allocated(dt_warm_full)) allocate(dt_warm_full (nlat,nlon,nfldsfc))
  if(.not.allocated(tref_full))    allocate(tref_full    (nlat,nlon,nfldsfc))

  nstinfo=4  ! number of fields appended to diag files
  if( ntguesnst < 1 .or. ntguesnst > nfldnst ) then
      call perr('geos_nstinit','ntguesnst = ',ntguesnst)
      call  die('geos_nstinit')
  endif
       
  if(mype==0) then
     write(6,*)'geos_nst_init_: default nstinfo:', nstinfo
  endif
end subroutine init_
!*******************************************************************************************

subroutine set_

use kinds,            only: i_kind, r_kind, r_single
use constants,        only: zero
use MAPL_Mod,         only: MAPL_UNDEF
use mpimod,           only: mpi_comm_world,ierror,mpi_rtype,mype
use gsi_nstcouplermod,only: tref_full,dt_cool_full,dt_warm_full,z_c_full,z_w_full
use satthin,          only: isli_full
use guess_grids,      only: nfldsfc,nfldnst,ntguesnst   !, sfct
use gridmod,          only: nlat,nlon,lat2,lon2,lat1,lon1,iglobal,itotsub,ijn,displs_g,strip
use general_commvars_mod, only: ltosi,ltosj
use gsi_metguess_mod, only: gsi_metguess_bundle
use gsi_bundlemod,    only: gsi_bundlegetpointer
!use m_rerank, only: rerank
implicit none
! local variables
integer(i_kind) it,ier
real(r_kind),pointer,dimension(:,:)  :: ptr2d=>NULL()
real                                 :: spval    ! spval should be set to 0.0 for all variables except tref. use tskin for tref. do this later. 

real(r_single),allocatable, dimension(:,:,:) :: tdel_full
real(r_kind),parameter                     :: TOL_dt_warm = 1.e-10_r_kind

allocate(tdel_full(nlat,nlon,nfldsfc))

!integer(i_kind) mold2(2,2)

  do it=1, nfldsfc

     call gsi_bundlegetpointer(gsi_metguess_bundle(it),'z_c',ptr2d,ier) ! get pointer to distrubuted 2d field

     if(ier==0) then
        call     dist2full_(z_c_full(:,:,it))
        spval = zero; call apply_ls_mask_(z_c_full(:,:,it), spval)
     else
        if(mype == 0)write(6,*)'trouble getting z_c in geos_nst_set_'
        CALL STOP2(999)
     endif

     call gsi_bundlegetpointer(gsi_metguess_bundle(it),'z_w',ptr2d,ier)
     if(ier==0) then
        call     dist2full_(z_w_full(:,:,it))
        spval = 2.0_r_kind; call apply_ls_mask_(z_w_full(:,:,it), spval)
     else
        if(mype == 0)write(6,*)'trouble getting z_w in geos_nst_set_'
        CALL STOP2(999)
     endif

     call gsi_bundlegetpointer(gsi_metguess_bundle(it),'dt_cool',ptr2d,ier)      
     if(ier==0) then
        call     dist2full_(dt_cool_full(:,:,it))
        spval = zero; call apply_ls_mask_(dt_cool_full(:,:,it), spval)
     else
        if(mype == 0)write(6,*)'trouble getting dt_cool in geos_nst_set_'
        CALL STOP2(999)
     endif

     call gsi_bundlegetpointer(gsi_metguess_bundle(it),'tdel',ptr2d,ier)      
     if(ier==0) then
        call     dist2full_(tdel_full(:,:,it))
        spval = 290.0_r_kind; call apply_ls_mask_(tdel_full(:,:,it), spval)
     else
        if(mype == 0)write(6,*)'trouble getting tdel in geos_nst_set_'
        CALL STOP2(999)
     endif

     call gsi_bundlegetpointer(gsi_metguess_bundle(it),'tref',ptr2d,ier)
     if(ier==0) then
        call     dist2full_(tref_full(:,:,it))
        spval = 290.0_r_kind; call apply_ls_mask_(tref_full(:,:,it), spval)
     else
        if(mype == 0)write(6,*)'trouble getting tref in geos_nst_set_'
        CALL STOP2(999)
     endif

     where (isli_full == 0)
         dt_warm_full(:,:,it) = tdel_full(:,:,it) - tref_full(:,:,it)
     else where
         dt_warm_full(:,:,it) = zero
     end where

!    dt_warm MUST be >= 0. There should BE NO negative values. Make sure this is the case.
     where( dt_warm_full(:,:,it) < TOL_dt_warm) 
        dt_warm_full(:,:,it) = zero
     end where
     
!    If you need Tskin, you have to do following for sfct
!    ptr2d => rerank(sfct(:,:,it),mold2,(/size(sfct,1),size(sfct,2)/)) 
!    call dist2full_(tskin_full(:,:,it))

  enddo

  deallocate(tdel_full)

  contains

  subroutine dist2full_(var_full)
  ! Gather variable to all processors (based on handling of sfc fields in satthin[getsfc])

  real(r_single),intent(inout),dimension(:,:) :: var_full

  real(r_kind),allocatable,  dimension(:)   :: aux
  real(r_kind),allocatable,  dimension(:)   :: work1
  real(r_kind),allocatable,  dimension(:,:) :: work2

  integer(i_kind) i,j,k,mm1

     if(.not.allocated(aux))   allocate(aux(lon1*lat1))
     if(.not.allocated(work1)) allocate(work1(itotsub))
     if(.not.allocated(work2)) allocate(work2(lat2,lon2))

     mm1=mype+1

     do j=1,lon1*lat1
        aux(j)=zero
     end do

     do j=1,lon2
        do i=1,lat2
           work2(i,j)=ptr2d(i,j)
        end do
     end do

     call strip(work2,aux)

     call mpi_allgatherv(aux,ijn(mm1),mpi_rtype,      &
                         work1,ijn,displs_g,mpi_rtype,&
                         mpi_comm_world,ierror)

     if(ier/=0) then
          PRINT *, 'dist2full_ mpi_allgatherv ierror:', ierror
          CALL STOP2(999)
     endif     

     do k=1,iglobal
        i=ltosi(k) ; j=ltosj(k)
        var_full(i,j)= work1(k)  
     end do

     if(allocated(aux))       deallocate(aux)
     if(allocated(work1))   deallocate(work1)
     if(allocated(work2))   deallocate(work2)

  end subroutine dist2full_

  subroutine apply_ls_mask_(var_full, spval)
  ! GEOS nst fields are computed over ocean ONLY. They are MAPL_UNDEF over lakes and inland water bodies.
  ! gsi has a single mask for all water (ocean, lakes, etc.) -> islip.  Here nst fields are made consistent with gsi's isli 

  real(r_single),intent(inout),dimension(:,:) :: var_full
  real,          intent(in)                   :: spval

  real, parameter                           :: tol         = 1.e-2_r_kind
! set Large_value to 1e4, because no nst_gsi bkg sfc z_w ~1e3, pick a value > 1e3.
  real, parameter                           :: Large_value = 1.e+4_r_kind
  integer                                   :: nUNDEF


  where( var_full >= Large_value)
        var_full = spval
  end where

  nUNDEF = count(abs(var_full/MAPL_UNDEF) >= tol)
  if( nUNDEF /=0) then
          PRINT *, 'MAPL_UNDEF not filled in geos_nstmod: apply_ls_mask_, num. of problem pts: ', nUNDEF
          CALL STOP2(999)
  endif


  end subroutine apply_ls_mask_

end subroutine set_
!*******************************************************************************************

subroutine final_
  use gsi_nstcouplermod, only: tref_full,dt_cool_full,dt_warm_full,z_c_full,z_w_full
  implicit none
!    clean up
     if(allocated(z_c_full))        deallocate(z_c_full)
     if(allocated(z_w_full))        deallocate(z_w_full)
     if(allocated(dt_cool_full))    deallocate(dt_cool_full)
     if(allocated(dt_warm_full))    deallocate(dt_warm_full)
     if(allocated(tref_full))       deallocate(tref_full)
end subroutine final_
!*******************************************************************************************

end module geos_nstmod

