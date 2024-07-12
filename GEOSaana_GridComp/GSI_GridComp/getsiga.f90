subroutine getsiga ()
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getsiga
!   prgmmr: todling
!
! abstract:  Calculate analysis errors from Lanczos-CG results
!
! program history log:
!   2010-03-16  todling  - initial code
!   2010-05-14  todling  - update to use gsi_bundle
!   2010-05-27  todling  - gsi_4dcoupler; remove all user-specific TL-related references
!   2010-08-19  lueken   - add only to module use;no machine code, so use .f90
!   2022-08-10  zhu      - add treatment for log(pbl*)
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use kinds, only: i_kind, r_kind
use mpimod, only: mype
use constants, only: zero,one
use gsi_4dvar, only: ibdate,lsqrtb,jsiga
use jfunc, only: jiter
use lanczos, only : congrad_siga
use state_vectors, only: allocate_state,deallocate_state
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_set
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_final
use gsi_bundlemod, only: gsi_bundle
implicit none
! declare local variables
character(len=*),parameter:: myname_ = "getsiga"
type(gsi_bundle)     :: siga                      ! vector to analysis errors
integer(i_kind)      :: nymd                      ! date as in YYYYMMDD
integer(i_kind)      :: nhms                      ! time as in HHMMSS
integer(i_kind)      :: nvecs
integer(i_kind)      :: ier

! consistency checks
if (jiter/=jsiga) return
if (.not.lsqrtb) then
   write(6,*)trim(myname_),': must set lsqrt=.t. to get analysis errors'
   call stop2(331)
end if

nymd = 10000*ibdate(1)+ibdate(2)*100+ibdate(3)
nhms = 10000*ibdate(4)
if(mype==0) write(6,'(2a,i8.8,2x,i6.6)')trim(myname_),': starting to calculate analysis errors at ',&
             nymd, nhms

! allocate memory for working arrays
call allocate_state(siga)

! calculate estimate of analysis errors
call congrad_siga(siga,nvecs,ier)

! write out analysis errors
if(ier==0) then
   call gsi_4dcoupler_putpert_set (nymd,nhms,ier)
   call gsi_4dcoupler_putpert (siga,nymd,nhms,'tlm','siga')
   if(mype==0) write(6,'(2a,i5,a)')trim(myname_),': complete calculating analysis errors using ',&
                                   nvecs, ' eigenvectors'
   call gsi_4dcoupler_putpert_final(ier)
else
   if(mype==0) write(6,'(2a,i6)')trim(myname_),': failed to calculate analysis errors, ier= ', ier
endif

! clean up
call deallocate_state(siga)

return
end subroutine getsiga

subroutine view_cv (xhat,mydate,filename,writecv)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    view_cv
!   prgmmr: todling
!
! abstract:  this allow writing CV to file for visualization
!
! program history log:
!   2011-02-23  todling  - initial code
!                          (not sure we'll keep this here)
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use kinds, only: i_kind, r_kind
use mpimod, only: mype
use constants, only: zero,one
use gsi_4dvar, only: nsubwin,lsqrtb
use state_vectors, only: allocate_state,deallocate_state
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_set
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_final
use gsi_bundlemod, only: gsi_bundle
use control_vectors, only: control_vector,write_cv
use state_vectors, only: allocate_state,deallocate_state,prt_state_norms
use bias_predictors, only: predictors,allocate_preds,deallocate_preds
use bias_predictors, only: write_preds
implicit none
type(control_vector)        :: xhat
integer(i_kind), intent(in) :: mydate(5) ! as in iadate or ibdate, or similar
character(len=*),intent(in) :: filename
logical,         intent(in) :: writecv   ! when .t., simply write out CV directly
! declare local variables
character(len=*),parameter:: myname_ = "view_cv"
integer(i_kind)      :: nymd                      ! date as in YYYYMMDD
integer(i_kind)      :: nhms                      ! time as in HHMMSS
integer(i_kind)      :: ii,ier
type(gsi_bundle) :: mval(nsubwin)
type(predictors) :: sbias

! in case CV not required to be transformed ...
if (writecv) then
   call write_cv(xhat,filename)
   return
else
!  if(mype==0) write(6,*) trim(myname_),': not writing CV to disk for now'
!  return
endif

! otherwise, transform CV to state-space and write out ...
nymd = 10000*mydate(1)+mydate(2)*100+mydate(3)
nhms = 10000*mydate(4)
if(mype==0) write(6,'(2a,i8.8,2x,i6.6)')trim(myname_),': start writing state on ',&
             nymd, nhms

! Allocate local variables
do ii=1,nsubwin
   call allocate_state(mval(ii))
end do
call allocate_preds(sbias)

if (lsqrtb) then
   call control2model(xhat,mval,sbias)
else
   call control2state(xhat,mval,sbias)
endif

! write out analysis errors
call gsi_4dcoupler_putpert_set (nymd,nhms,ier)
do ii=1,nsubwin
   call gsi_4dcoupler_putpert (mval(ii),nymd,nhms,'tlm',filename) ! will need to be smart for nsubwin>1
   call prt_state_norms(mval(ii),'output-state')
enddo
call gsi_4dcoupler_putpert_final (ier)
call write_preds(sbias,'preds_'//trim(filename),mype)

! Allocate local variables
call deallocate_preds(sbias)
do ii=1,nsubwin
   call deallocate_state(mval(ii))
end do

if(mype==0) write(6,'(3a)')trim(myname_),': complete writing state ', trim(filename)

return
end subroutine view_cv

subroutine view_cv_ad (xhat,mydate,filename,readcv)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    view_cv_ad
!   prgmmr: todling
!
! abstract:  this allows reading state/control vector and transforming to CV
!
! program history log:
!   2012-05-22  todling  - based on view_ad, performs opposite of view_cv
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use kinds, only: i_kind, r_kind
use mpimod, only: mype
use constants, only: zero,one
use gsi_4dvar, only: nsubwin,lsqrtb
use state_vectors, only: allocate_state,deallocate_state
use gsi_4dcouplermod, only: gsi_4dcoupler_getpert
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: gsi_bundlegetpointer
use gsi_bundlemod, only: gsi_bundledup
use gsi_bundlemod, only: assignment(=)
use control_vectors, only: control_vector,read_cv,assignment(=)
use state_vectors, only: allocate_state,deallocate_state,prt_state_norms
use bias_predictors, only: predictors,allocate_preds,deallocate_preds,assignment(=)
use bias_predictors, only: read_preds
implicit none
type(control_vector)        :: xhat
integer(i_kind), intent(in) :: mydate(5) ! as in iadate or ibdate, or similar
character(len=*),intent(in) :: filename
logical,         intent(in) :: readcv    ! when .t. simply read in CV
! declare local variables
character(len=*),parameter:: myname_ = "view_cv_ad"
integer(i_kind)      :: nymd                      ! date as in YYYYMMDD
integer(i_kind)      :: nhms                      ! time as in HHMMSS
integer(i_kind)      :: ii
type(gsi_bundle) :: mval(nsubwin)
type(predictors) :: sbias

! in case CV not required to be transformed ...
if (readcv) then
   call read_cv(xhat,filename)
   return
else ! for now only
   xhat=zero
   if(mype==0) write(6,*) trim(myname_),': input vector set to zero for now'
   return
endif

! otherwise read state-vector and transform to control-space ...
nymd = 10000*mydate(1)+mydate(2)*100+mydate(3)
nhms = 10000*mydate(4)
if(mype==0) write(6,'(2a,i8.8,2x,i6.6)')trim(myname_),': start reading state on ',&
             nymd, nhms

! Allocate local variables
do ii=1,nsubwin
   call allocate_state(mval(ii))
   mval(ii)=zero
end do
call allocate_preds(sbias)
sbias=zero
xhat=zero

if (lsqrtb) then
   call control2model(xhat,mval,sbias) ! dirty trick
endif
! read in (model/state) vector
do ii=1,nsubwin
   mval(ii)=zero
   call gsi_4dcoupler_getpert (mval(ii),'tlm',filename) ! will need better for nsubwin>1
   call prt_state_norms(mval(ii),'input-state')
enddo
call read_preds(sbias,'preds_'//trim(filename))

! convert to control vector
if (lsqrtb) then
   call control2model_ad(mval,sbias,xhat)
else
   call control2state_ad(mval,sbias,xhat)
endif

! Allocate local variables
call deallocate_preds(sbias)
do ii=1,nsubwin
   call deallocate_state(mval(ii))
end do

if(mype==0) write(6,'(3a)')trim(myname_),': complete reading state ', trim(filename)
return
end subroutine view_cv_ad

subroutine view_st (ssval,filename)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    view_st
!   prgmmr: todling
!
! abstract:  this allow writing CV to file for visualization
!
! program history log:
!   2011-02-23  todling  - initial code
!                          (not sure we'll keep this here)
!   2020-02-26  todling   obsbin time now in minute
!   2020-03-03  todling   split set/final from put
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use kinds, only: i_kind, r_kind
use mpimod, only: mype
use constants, only: zero,one
use constants, only: max_varname_length
use gsi_4dvar, only: ibdate,nobs_bins,nmn_obsbin,mn_obsbin
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_set
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_final
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: gsi_bundlegetpointer
use gsi_bundlemod, only: gsi_bundledup
use gsi_metguess_mod, only: gsi_metguess_bundle
use gsi_metguess_mod, only: gsi_metguess_get
use mpeu_util, only: getindex
use guess_grids, only: nfldsig,hrdifsig
use state_vectors, only: svars2d
use gridmod, only: lat2,lon2
implicit none
type(gsi_bundle)            :: ssval(nobs_bins)
character(len=*),intent(in) :: filename
! declare local variables
character(len=*),parameter:: myname_ = "view_st"
character(max_varname_length),allocatable,dimension(:) :: guess
type(gsi_bundle)     :: sval(nobs_bins)
integer(i_kind)      :: nymd                      ! date as in YYYYMMDD
integer(i_kind)      :: nhms                      ! time as in HHMMSS
integer(i_kind)      :: ii,status
integer(i_kind)      :: mydate(5)
integer(i_kind)      :: i,j,it,nguess,istatus,id,ic

real(r_kind),pointer,dimension(:,:  ) :: ptr2dinc =>NULL()
real(r_kind),pointer,dimension(:,:  ) :: ptr2dges =>NULL()
integer(i_kind),dimension(8) :: ida,jda
real(r_kind),dimension(5)    :: fha
real(r_kind) zt,dtmp,wktmp


! Inquire about guess fields
  call gsi_metguess_get('dim',nguess,istatus)
  if (nguess>0) then
     allocate(guess(nguess))
     call gsi_metguess_get('gsinames',guess,istatus)
  endif

  do it=1,nfldsig
     if (nobs_bins>1) then
        zt = hrdifsig(it)
        ii = NINT(zt*60/mn_obsbin)+1
     else
        ii = 1
     endif

     call gsi_bundledup (ssval(ii), sval(ii), ' copy of bundle ', istatus)
     do ic=1,nguess     
        id=getindex(svars2d,guess(ic))
        if (id>0) then
           if (trim(guess(ic))=='pblri' .or. trim(guess(ic))=='pblrf' .or. trim(guess(ic))=='pblkh') then
              call gsi_bundlegetpointer (sval(ii), guess(ic),ptr2dinc,istatus)
              if (istatus.ne.0) then
                 write(6,*) "yeg_view_st: no guess(ic) ",guess(ic)," in sval(ii) (increment)"
              end if
              call gsi_bundlegetpointer (gsi_metguess_bundle(it),guess(ic),ptr2dges,istatus)
              if (istatus.ne.0) then
                 write(6,*) "yeg_view_st: no guess(ic) ",guess(ic)," in gsi_metguess_bundle(it)"
              end if
              do j=1,lon2
                 do i=1,lat2
                    ! Convert Increment from log space to pysical space
                    !if (ptr2dinc(i,j) > 0.005_r_kind) then
                    !if (abs(ptr2dinc(i,j)) > 0.0005_r_kind) then ! log
                    if (abs(ptr2dinc(i,j)) > 0.5_r_kind) then ! physical space
                       write(6,*) i,j,'view_st_1:fg and inc=',ptr2dges(i,j),ptr2dinc(i,j)
                    else
                       ptr2dinc(i,j)=zero
                    end if
                 end do
              end do
              cycle
           end if          
        end if
     end do
  end do

! initial date/time and writer
mydate = ibdate
nymd = 10000*mydate(1)+mydate(2)*100+mydate(3)
nhms = 10000*mydate(4)+mydate(5)*100

call GSI_4dCoupler_putpert_set(nymd,nhms,status)

! write out analysis errors/increments
do ii=1,nobs_bins
   nymd = 10000*mydate(1)+mydate(2)*100+mydate(3)
   nhms = 10000*mydate(4)+mydate(5)*100
   ! iwrtinc ...

   if(mype==0) write(6,'(2a,i8.8,2x,i6.6)')trim(myname_),': start writing state on ', nymd, nhms
   call gsi_4dcoupler_putpert (sval(ii),nymd,nhms,'tlm',filename)
!  call prt_state_norms(sval(ii),'output-state')

   ! increment mydate ...
   fha(:)=0.0; ida=0; jda=0
   fha(3)=nmn_obsbin! relative time interval in minutes
   ida(1)=mydate(1) ! year
   ida(2)=mydate(2) ! month
   ida(3)=mydate(3) ! day
   ida(4)=0         ! time zone
   ida(5)=mydate(4) ! hour
   ida(6)=mydate(5) ! hour
   ! Move date-time forward by nhr_assimilation hours
   call w3movdat(fha,ida,jda)
   mydate(1)=jda(1)
   mydate(2)=jda(2)
   mydate(3)=jda(3)
   mydate(4)=jda(5)
   mydate(5)=jda(6)
enddo
! finalize writer
call GSI_4dCoupler_putpert_final(status)

if(mype==0) write(6,'(3a)')trim(myname_),': complete writing state ', trim(filename)

return
end subroutine view_st
