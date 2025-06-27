module m_revBens
use kinds, only: r_kind,i_kind
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: self_add
implicit none
private

public :: revBens_ensmean_overwrite
public :: revBens_init
public :: revBens_final
public :: revBens_ensloc_refactor
public :: spread2d,spread3d
public :: update_spread

logical :: update_spread=.false.
real(r_kind),allocatable,dimension(:,:,:) :: spread2d
real(r_kind),allocatable,dimension(:,:,:,:) :: spread3d

! All internals (not to be made public)
character(len=*), parameter ::  myname='m_revBens'
logical,parameter:: verbose =.true.
integer(i_kind)  :: miter=-1
type(gsi_bundle),save,allocatable :: xincaccum(:)
real(r_kind),allocatable::sprdfactors2d(:,:)
real(r_kind),allocatable::sprdfactors3d(:,:,:)
contains

subroutine revBens_init(miter_in)
  use mpimod,only: mype
  use gridmod,only: nsig
  use gsi_4dvar, only: ens_fmnlevs
  use control_vectors, only: nc2d,nc3d
  use hybrid_ensemble_parameters, only: upd_ens_spread,write_ens_sprd
  implicit none
  integer(i_kind),intent(in) :: miter_in
  integer(i_kind) ii,nt
  character(len=*), parameter :: myname_ = myname//":init"
  miter=miter_in
  nt=size(ens_fmnlevs)
  if (.not.write_ens_sprd) then
    if(mype==0) then
       write(6,*) myname_, '; WARNING: spread calc not active, will not calculate localization (re)factors'
    endif
    return
  endif
  if(.not.allocated(spread2d)) then
     allocate(spread2d(nc2d,nt,miter))
  endif
  if(.not.allocated(spread3d)) then
     allocate(spread3d(nsig,nc3d,nt,miter))
  endif
  if(.not.allocated(sprdfactors2d)) allocate(sprdfactors2d(nc2d,nt))
  if(.not.allocated(sprdfactors3d)) allocate(sprdfactors3d(nsig,nc3d,nt))
end subroutine revBens_init

subroutine revBens_ensmean_overwrite (en_bar,ibin)
  use constants, only: zero,half,one,rd_over_cp,fv,qcmin,qmin,tgmin
  use mpimod, only: npe,mype
  use mpeu_util, only: die
  use gridmod, only: idsl5
  use hybrid_ensemble_parameters, only: grd_ens
  use hybrid_ensemble_parameters, only: q_hyb_ens
  use hybrid_ensemble_parameters, only: oz_univ_static
  !use hybrid_ensemble_parameters, only: sst_staticB,pblh_staticB
  use hybrid_ensemble_parameters, only: sst_staticB,pblri_staticB,pblrf_staticB,pblkh_staticB
  use gsi_bundlemod, only: gsi_bundlecreate
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_bundlemod, only: gsi_bundlegetvar
  use gsi_bundlemod, only: gsi_bundledestroy
  use gsi_bundlemod, only: assignment(=)
  use gsi_4dcouplermod, only: gsi_4dcoupler_getpert_set
  use gsi_4dcouplermod, only: gsi_4dcoupler_getpert
  use gsi_4dcouplermod, only: gsi_4dcoupler_getpert_unset
  use gsi_enscouplermod, only: gsi_enscoupler_get_user_ens
  use gsi_enscouplermod, only: gsi_enscoupler_create_sub2grid_info
  use gsi_enscouplermod, only: gsi_enscoupler_destroy_sub2grid_info
  use general_sub2grid_mod, only: sub2grid_info
  use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
  use jfunc, only: jiter
  use gsi_4dvar, only: ens_fmnlevs
  use state_vectors,only: dot_product

  implicit none

  type(gsi_bundle) :: en_bar
  integer(i_kind),intent(in) :: ibin

! local variables
  character(len=*), parameter :: myname_ = myname//":ensmean_overwrite"
  character(len=8)  :: xincfile
  character(len=40) :: evar
  logical :: hydrometeor
  integer(i_kind) :: im,jm,km,ii,ierges,ierbar,ierinc,istatus,ier
  integer(i_kind) :: nt,myjiter
  real(r_kind) :: ddot
  real(r_kind),pointer,dimension(:,:)   :: iptr2d =>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: iptr3d =>NULL()
  real(r_kind),pointer,dimension(:,:)   :: gptr2d =>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: gptr3d =>NULL()
  real(r_kind),pointer,dimension(:,:)   :: mptr2d =>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: mptr3d =>NULL()
  type(gsi_bundle)    :: guess
  type(gsi_bundle)    :: xinc
  type(sub2grid_info) :: grd_tmp

  ierinc=0
  im=en_bar%grid%im
  jm=en_bar%grid%jm
  km=en_bar%grid%km
  myjiter=jiter
  nt=size(ens_fmnlevs)

  hydrometeor = trim(cvars3d(ii))=='cw' .or. trim(cvars3d(ii))=='ql' .or. &
                trim(cvars3d(ii))=='qi' .or. trim(cvars3d(ii))=='qr' .or. &
                trim(cvars3d(ii))=='qs' .or. trim(cvars3d(ii))=='qg' .or. &
                trim(cvars3d(ii))=='qh'

  ! Define grid for local entities
  call gsi_enscoupler_create_sub2grid_info(grd_tmp,km,npe,grd_ens)

  if( update_spread ) then
     ! Create bundle to hold increment
     call gsi_bundlecreate(xinc,en_bar%grid,'inc',istatus,names2d=cvars2d,names3d=cvars3d)
     if (istatus/=0)then
        call die(myname_,': cannot create inc bundle, aborting ...',90)
     endif
     ! Read in increment
     if (mype==0) then
         write(6,*) myname_,': reading analysis increment at time slot ',ibin
         write(6,*) myname_,': resolution = ', grd_tmp%nlon, grd_tmp%nlat, grd_tmp%nsig
     endif
     !call GSI_EnsCoupler_get_user_ens(grd_tmp,-(myjiter-1),ibin,xinc,istatus)
     xincfile='xinc.ZZZ'
     write(xincfile(6:8),'(i3.3)') myjiter-1
     call GSI_4dCoupler_getpert_set(grd_tmp%nlat,grd_tmp%nlon,grd_tmp%lat2,grd_tmp%lon2,ibin)
     call GSI_4dcoupler_getpert(xinc,ibin,'inc',trim(xincfile))
     call GSI_4dCoupler_getpert_unset()
!    Accummulate increment when needed
     if (miter>2) then
        ! create bundle to accummulate increment
        if(.not.allocated(xincaccum)) then
           allocate(xincaccum(nt))
           do ii=1,nt
             !call gsi_bundlecreate(xincaccum(ii),xinc,'accum-inc',istatus)
              call gsi_bundlecreate(xincaccum(ii),en_bar%grid,'accum-inc',istatus,names2d=cvars2d,names3d=cvars3d)
              xincaccum(ii)=zero
           enddo
        endif
        ddot = dot_product(xinc,xinc)
        if(mype==0) then
          write(6,'(2a,1p,e25.18)') myname_, ": before accum, dot(xinc,xinc) = ", ddot
        endif
        call self_add(xincaccum(ibin),xinc)
        xinc=xincaccum(ibin)
     endif
     ddot = dot_product(xinc,xinc)
     if(mype==0) then
       write(6,'(2a,1p,e25.18)') myname_, ": after  accum, dot(xinc,xinc) = ", ddot
     endif
  else
     ! Read in state to recenter around
     if (mype==0) then
         write(6,*) myname_,': reading background at time slot ',ibin
         write(6,*) myname_,': resolution = ', grd_tmp%nlon, grd_tmp%nlat, grd_tmp%nsig
     endif
!    Create bundle to hold guess or analysis increment
     call gsi_bundlecreate(guess,en_bar%grid,'bkg',istatus,names2d=cvars2d,names3d=cvars3d)
     if (istatus/=0)then
        call die(myname_,': cannot create bundle, aborting ...',91)
     endif
     call GSI_EnsCoupler_get_user_ens(grd_tmp,0,-1,ibin,guess,istatus)
     if (istatus/=0)then
        call die(myname_,': cannot create guess bundle, aborting ...',92)
     endif
  endif

!  Rank-2 fields
   ierbar=0;ierges=0;ierinc=0
   do ii=1,nc2d
      evar=trim(cvars2d(ii)) ! in general output name same as input (but not always!)
      call gsi_bundlegetpointer (en_bar,evar,mptr2d,ierbar)
      if (update_spread) then
         call gsi_bundlegetpointer (xinc ,evar,iptr2d,ierinc)
      else
         call gsi_bundlegetpointer (guess,evar,gptr2d,ierges)
      endif
      if(ierges==0 .and. ierbar==0 .and. ierinc==0) then
         if (update_spread) then
            mptr2d=mptr2d+iptr2d
         else
            mptr2d=gptr2d
         endif
         if (evar=='sst' .and. sst_staticB) then
            mptr2d=zero
         endif
         !if ((evar=='pblri' .or. evar=='pblrf' .or. evar=='pblkh') .and. pblh_staticB) then
         if ( (evar=='pblri' .and. pblri_staticB) .or. &
              (evar=='pblrf' .and. pblrf_staticB) .or. &
              (evar=='pblkh' .and. pblkh_staticB) ) then
            mptr2d=zero
         endif
         if (mype==0) then
            write(6,*) myname_,': mean ens field ',trim(evar),' overwritten'
         endif
      else
         if (ierges/=0.and.mype==0) then
            write(6,*) myname_,': field ',evar,' not in replacement mean CV, skipping'
         endif
         if (ierbar/=0.and.mype==0) then
            write(6,*) myname_,': field ',evar,' not in ens mean CV, skipping'
         endif
      endif
   enddo

!  Rank-3 fields
   ierbar=0;ierges=0;ierinc=0
   do ii=1,nc3d
      evar=trim(cvars3d(ii)) ! in general output name same as input (but not always!)
      call gsi_bundlegetpointer (en_bar,evar,mptr3d,ierbar)
      if (update_spread) then
         call gsi_bundlegetpointer (xinc ,evar,iptr3d,ierinc)
      else
         call gsi_bundlegetpointer (guess,evar,gptr3d,ierges)
      endif
      if(ierges==0 .and. ierbar==0 .and. ierinc==0) then
         if (update_spread) then
            if (evar=='q') then
               ! TODO:
               ! the update should be something like this:
               ! mptr3d = iptr3d/qs - (dqs/qs)*mptr3d
               ! where mptr3d=rh and iptr3d=dq
               ! but I don''t know how to get qs (since en_bar carries rh),
               ! and I don''t know how to get dqs - would have to develop 
               ! tlm for genqsat - don''t feel like it.
               ! So, this simply leaves the ensemble mean of q untouched.
            else
               mptr3d=mptr3d+iptr3d
            endif
         else
            mptr3d=gptr3d
            if(evar=='q') call positive_bounds_(qmin)
         endif
         if(evar=='oz')  call positive_bounds_(tgmin)
         if(hydrometeor) call positive_bounds_(qcmin)
         if ( evar == 'oz' .and. oz_univ_static ) then
            mptr3d = zero
            cycle
         end if
         if (mype==0) then
            write(6,*) myname_,': mean ens field ',trim(evar),' overwritten'
         endif
      else
         if (ierges/=0.and.mype==0) then
            write(6,*) myname_,': field ',evar,' not in replacement mean CV, skipping'
         endif
         if (ierbar/=0.and.mype==0) then
            write(6,*) myname_,': field ',evar,' not in ens mean CV, skipping'
         endif
         if (ierinc/=0.and.mype==0) then
            write(6,*) myname_,': field ',evar,' not in increment CV, skipping'
         endif
      endif

   enddo

!  Now adjust for variable redefinition/overload
   if(.not.update_spread) then
      if(.not.q_hyb_ens) call q2rh_
   endif

!  Clean up
   if(update_spread) then
     call gsi_bundledestroy(xinc ,istatus)
   else
     call gsi_bundledestroy(guess,istatus)
   endif
   call gsi_enscoupler_destroy_sub2grid_info(grd_tmp)

contains
 
  subroutine positive_bounds_(threshold)

  real(r_kind),intent(in) :: threshold
  integer(i_kind) i,j,k

!$omp parallel do schedule(dynamic,1) private(i,j,k)
  do k=1,km
     do j=1,jm
        do i=1,im
           mptr3d(i,j,k) = max(mptr3d(i,j,k),threshold)
        end do
     end do
  end do

  end subroutine positive_bounds_

  subroutine q2rh_

! Really bad to repeat code ...
  character(len=*), parameter :: myname__ = myname_//":q2rh_"

! Note that, in most cases, the ensemble mean is not physical, so doing 
! calculations such as the one here might not be meaningful, still ...

  real(r_kind),pointer,dimension(:,:)   :: ps
  real(r_kind),pointer,dimension(:,:,:) :: tv
  real(r_kind),pointer,dimension(:,:,:) :: q
  real(r_kind),allocatable,dimension(:,:,:) :: tsen,prsl,pri,qs

  real(r_kind) kap1,kapr,rh
  integer(i_kind) i,j,k,iderivative
  logical ice

  kap1=rd_over_cp+one
  kapr=one/rd_over_cp
 
  istatus=0
  call gsi_bundlegetpointer(en_bar,'ps',ps,ier);istatus=ier
  call gsi_bundlegetpointer(en_bar,'t' ,tv,ier);istatus=istatus+ier
  call gsi_bundlegetpointer(en_bar,'q' ,q ,ier);istatus=istatus+ier
  if (istatus/=0) then
     call die(myname_,': trouble retriving pntr from en_bar, istatus= ',istatus)
  endif

! Compute RH

! Get 3d pressure field now on interfaces
  allocate(pri(im,jm,km+1))
  call general_getprs_glb(ps,tv,pri)
  allocate(prsl(im,jm,km),tsen(im,jm,km),qs(im,jm,km))

! Get sensible temperature and 3d layer pressure
  if (idsl5 /= 2) then
!$omp parallel do schedule(dynamic,1) private(k,j,i)
     do k=1,km
        do j=1,jm
           do i=1,im
              prsl(i,j,k)=((pri(i,j,k)**kap1-pri(i,j,k+1)**kap1)/&
                      (kap1*(pri(i,j,k)-pri(i,j,k+1))))**kapr
              tsen(i,j,k)=tv(i,j,k)/(one+fv*max(zero,q(i,j,k)))
           end do
        end do
     end do
  else
!$omp parallel do schedule(dynamic,1) private(k,j,i)
     do k=1,km
        do j=1,jm
           do i=1,im
               prsl(i,j,k)=(pri(i,j,k)+pri(i,j,k+1))*half
               tsen(i,j,k)=tv(i,j,k)/(one+fv*max(zero,q(i,j,k)))
           end do
        end do
     end do
  end if
  deallocate(pri)

  ice=.true.
  iderivative=0
  call genqsat(qs,tsen,prsl,im,jm,km,ice,iderivative)
  deallocate(tsen,prsl)

  do k=1,km
     do j=1,jm
        do i=1,im
           q(i,j,k)=q(i,j,k)/qs(i,j,k)
        end do
     end do
  end do
  deallocate(qs)

  if (mype==0) then
     write(6,'(2a)') myname__,': converted q to rh after overwriting ensemble mean'
  endif
  end subroutine q2rh_

end subroutine revBens_ensmean_overwrite

subroutine revBens_ensloc_refactor(hlocs,vlocs)
  use mpimod, only: mype
  use constants, only: zero,one
  use gsi_4dvar, only: ens_fmnlevs
  use gridmod,only: nsig
  use hybrid_ensemble_parameters, only: upd_ens_spread,upd_ens_localization,write_ens_sprd
  use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
  use jfunc, only: jiter
  implicit none
  real(r_kind),intent(inout) :: hlocs(nsig)
  real(r_kind),intent(inout) :: vlocs(nsig)

  integer(i_kind) n,k,myjiter,nt,ibin
  real(r_kind) scaleit
  character(len=*), parameter :: myname_ = "revBens_ensloc_refactor"

  if(.not.upd_ens_spread) return  ! nothing to do
  if(.not.write_ens_sprd) then
    if(mype==0) then
       write(6,*) myname_, '; WARNING: average spreads not available, cannot calculate spread factors'
    endif
    return
  endif
  myjiter=jiter
  nt=size(ens_fmnlevs)

! Write out summary of spreads
  if (verbose) then
     if (mype==0) then
        do ibin=1,nt
           write(6,*)'GLO_AVE_SPREAD_DUALRES_2D_FIELDS, TIME BIN: ', ibin
           do n=1,nc2d
              write(6,'(i5,2x,a10,2x,1p,e13.6)') myjiter, trim(cvars2d(n)), spread2d(n,ibin,myjiter)
           enddo
           write(6,*)'GLO_AVE_SPREAD_DUALRES_3D_FIELDS, TIME BIN:', ibin
           write(6,'(14x,9(a11))') (trim(cvars3d(n)),n=1,nc3d)
           do k=1,nsig
              write(6,'(2(2x,i5),1p,9(e11.4))') myjiter, k, (spread3d(k,n,ibin,myjiter),n=1,nc3d)
           enddo
        enddo
     endif
  endif

  if(myjiter<=1) return

! Calculate spread change factors

  do ibin=1,nt

!    2D
     do n=1,nc2d
        if(spread2d(n,ibin,1)>zero) then
           sprdfactors2d(n,ibin) = spread2d(n,ibin,myjiter)/spread2d(n,ibin,1)
        else
           sprdfactors2d(n,ibin) = one
        endif
     enddo
     if (mype==0.and.verbose) then
        write(6,*)'GLO_AVE_SPREAD_FACTORS_2D, TIME BIN: ', ibin
        do n=1,nc2d
           write(6,'(1x,a10,2x,1p,e13.6)') trim(cvars2d(n)), sprdfactors2d(n,ibin)
        enddo
     endif
!    3D
     if (mype==0.and.verbose) then
        write(6,*)'GLO_AVE_SPREAD_FACTORS_3D, TIME BIN: ', ibin
        write(6,'(11x,i5,9(a11))') myjiter, (trim(cvars3d(n)),n=1,nc3d)
     endif
     do n=1,nc3d
        do k=1,nsig
           if(spread3d(k,n,ibin,1)>zero) then
              sprdfactors3d(k,n,ibin) = spread3d(k,n,ibin,myjiter)/spread3d(k,n,ibin,1)
           else
              sprdfactors3d(k,n,ibin) = one
           endif
        enddo
     enddo
     if (mype==0.and.verbose) then
        do k=1,nsig
           write(6,'(2(1x,i5),8x,1p,9(e11.4))') myjiter, k, (sprdfactors3d(k,n,ibin),n=1,nc3d)
        enddo
     endif

  enddo ! nt

  scaleit = one
  do n=1,nc2d
     if(trim(cvars2d(n)) == 'ps' ) then
        scaleit = sum(sprdfactors2d(n,1:nt))/nt
     endif
  enddo
  if (mype==0) then
      write(6,'(2a,1p,e10.3)') myname_, ': localization scaling factor: ', scaleit
  endif
  if (upd_ens_localization) then
      hlocs = scaleit*hlocs
  endif

end subroutine revBens_ensloc_refactor

subroutine revBens_final
  use kinds, only: i_kind
  use gsi_bundlemod, only: gsi_bundledestroy
  implicit none
  integer(i_kind) ii,istatus
  if (allocated(xincaccum)) then
     do ii=size(xincaccum),1,-1
        call gsi_bundledestroy(xincaccum(ii),istatus)
     enddo
     deallocate(xincaccum)
  endif
  if(allocated(spread2d)) deallocate(spread2d)
  if(allocated(spread3d)) deallocate(spread3d)
  if(allocated(sprdfactors2d)) deallocate(sprdfactors2d)
  if(allocated(sprdfactors3d)) deallocate(sprdfactors3d)
  
end subroutine revBens_final
end module m_revBens
