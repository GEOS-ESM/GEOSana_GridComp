module obs_sensitivity
!$$$ module documentation block
!           .      .    .                                       .
! module:   obs_sensitivity
!   prgmmr: tremolet
!
! abstract: Contains variables and routines for computation of
!           forecast sensitivity to observations.
!
! program history log:
!   2007-06-26 tremolet
!   2007-07-19 tremolet - increment sensitivity to observations
!   2009-08-07 lueken   - updated documentation
!   2010-04-30 tangborn - add pointer to carbon monoxide
!   2010-05-27 todling  - remove all user-specific TL-related references
!   2010-07-16 todling  - add reference to aero and aerol
!   2010-08-19 lueken   - add only to module use;no machine code, so use .f90
!   2011-03-29 todling  - add reference to pm2_5
!   2012-04-15 todling  - add reference to gust, vis, pblh
!   2015-07-10 pondeca  - add reference to wspd10m, td2m ,mxtm ,mitm ,pmsl,
!                         howv ,tcamt, lcbas, cldch
!   2016-02-20 pagowski - add pm10
!   2016-05-05 pondeca  - add reference to uwnd10m, vwnd10m
!   2017-05-06 todling  - reload ensemble when FSOI calc doing EFSOI-like
!   2017-05-21 todling  - implement ability to do 2nd EFSOI-like calc
!
! Subroutines Included:
!   init_fc_sens  - Initialize computations
!
! Variable Definitions:
!   fcsens - forecast sensitivity gradient
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
! ------------------------------------------------------------------------------
use kinds, only: r_kind,i_kind,r_quad
use constants, only: zero, zero_quad, two
use gsi_4dvar, only: nobs_bins, l4dvar, lsqrtb, nsubwin
use gsi_4dvar, only: tau_fcst,efsoi_order,efsoi_afcst,efsoi_ana
use jfunc, only: jiter, miter, niter, iter
use obsmod, only: cobstype, nobs_type, obscounts, &
                  i_ps_ob_type, i_t_ob_type, i_w_ob_type, i_q_ob_type, &
                  i_spd_ob_type, i_rw_ob_type, i_dw_ob_type, &
                  i_sst_ob_type, i_pw_ob_type, i_pcp_ob_type, i_oz_ob_type, &
                  i_o3l_ob_type, i_gps_ob_type, i_rad_ob_type, i_tcp_ob_type, &
                  i_lag_ob_type, i_colvk_ob_type, i_aero_ob_type, i_aerol_ob_type, &
                  i_pm2_5_ob_type, i_gust_ob_type, i_vis_ob_type, i_pblh_ob_type, &
                  i_wspd10m_ob_type, i_td2m_ob_type, i_mxtm_ob_type, i_mitm_ob_type, &
                  i_pmsl_ob_type, i_howv_ob_type, i_tcamt_ob_type, i_lcbas_ob_type, &
                  i_cldch_ob_type, i_uwnd10m_ob_type, i_vwnd10m_ob_type, i_pm10_ob_type, &
                  i_swcp_ob_type, i_lwcp_ob_type

use mpimod, only: mype
use control_vectors, only: control_vector,allocate_cv,read_cv,deallocate_cv, &
    dot_product,assignment(=)
use state_vectors, only: allocate_state,deallocate_state
use gsi_bundlemod, only: self_add
use gsi_bundlemod, only: assignment(=)
use gsi_bundlemod, only: gsi_bundle
use bias_predictors, only: predictors,allocate_preds,deallocate_preds, &
    assignment(=)
use mpl_allreducemod, only: mpl_allreduce
use gsi_4dcouplermod, only: gsi_4dcoupler_getpert
use hybrid_ensemble_parameters,only : l_hyb_ens,ntlevs_ens,destroy_hybens_localization_parameters
use hybrid_ensemble_isotropic, only: create_ensemble,load_ensemble,destroy_ensemble
use hybrid_ensemble_isotropic, only: hybens_localization_setup
! ------------------------------------------------------------------------------
implicit none
save
private
public lobsensfc,lobsensjb,lobsensincr,lobsensadj,&
       lobsensmin,iobsconv,llancdone,lsensrecompute, &
       fcsens, sensincr, &
       init_obsens, init_fc_sens, save_fc_sens, dot_prod_obs
public efsoi_o2_update

logical lobsensfc,lobsensjb,lobsensincr, &
        lobsensadj,lobsensmin,llancdone,lsensrecompute
integer(i_kind) :: iobsconv

! ------------------------------------------------------------------------------
type(control_vector) :: fcsens
real(r_kind), allocatable :: sensincr(:,:,:)
character(len=5) :: cobtype(nobs_type)
integer(i_kind):: my_nobs_type
integer(i_kind):: orig_tau_fcst=-1
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
subroutine init_obsens
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_obsens
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-07  lueken - added subprogram doc block
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
implicit none

lobsensfc=.false.
lobsensjb=.false.
lobsensincr=.false.
lobsensadj=.false.
lobsensmin=.false.
lsensrecompute=.false.
llancdone=.false.
iobsconv=0

end subroutine init_obsens
! ------------------------------------------------------------------------------
subroutine init_fc_sens
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_fc_sens
!   prgmmr:      tremolet
!
! abstract: Read forecast sensitivity gradient
!
! program history log:
!   2007-06-26  tremolet - initial code
!   2009-08-07  lueken - added subprogram doc block
!   2010-05-27  todling - gsi_4dcoupler; remove dependence on GMAO specifics
!   2012-05-22  todling - update interface to getpert
!   2015-12-01  todling - add several obs-types that Pondeca forget to add here
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

implicit none

character(len=12) :: clfile
type(gsi_bundle) :: fcgrad(nsubwin)
type(gsi_bundle) :: eval(ntlevs_ens)
type(predictors) :: zbias
type(control_vector) :: xwork
real(r_kind) :: zjx
integer(i_kind) :: ii
character(len=80),allocatable,dimension(:)::fname

! Define short name for obs types
my_nobs_type=0
cobtype( i_ps_ob_type)   ="spr  "; my_nobs_type=my_nobs_type+1
cobtype(  i_t_ob_type)   ="tem  "; my_nobs_type=my_nobs_type+1
cobtype(  i_w_ob_type)   ="uv   "; my_nobs_type=my_nobs_type+1
cobtype(  i_q_ob_type)   ="hum  "; my_nobs_type=my_nobs_type+1
cobtype(i_spd_ob_type)   ="spd  "; my_nobs_type=my_nobs_type+1
cobtype( i_rw_ob_type)   ="rw   "; my_nobs_type=my_nobs_type+1
cobtype( i_dw_ob_type)   ="dw   "; my_nobs_type=my_nobs_type+1
cobtype(i_sst_ob_type)   ="sst  "; my_nobs_type=my_nobs_type+1
cobtype( i_pw_ob_type)   ="pw   "; my_nobs_type=my_nobs_type+1
cobtype(i_pcp_ob_type)   ="pcp  "; my_nobs_type=my_nobs_type+1
cobtype( i_oz_ob_type)   ="oz   "; my_nobs_type=my_nobs_type+1
cobtype(i_o3l_ob_type)   ="o3l  "; my_nobs_type=my_nobs_type+1
cobtype(i_gps_ob_type)   ="gps  "; my_nobs_type=my_nobs_type+1
cobtype(i_rad_ob_type)   ="rad  "; my_nobs_type=my_nobs_type+1
cobtype(i_tcp_ob_type)   ="tcp  "; my_nobs_type=my_nobs_type+1
cobtype(i_lag_ob_type)   ="lag  "; my_nobs_type=my_nobs_type+1
cobtype(i_colvk_ob_type) ="colvk"; my_nobs_type=my_nobs_type+1
cobtype(i_aero_ob_type)  ="aero "; my_nobs_type=my_nobs_type+1
cobtype(i_aerol_ob_type) ="aerol"; my_nobs_type=my_nobs_type+1
cobtype(i_pm2_5_ob_type) ="pm2_5"; my_nobs_type=my_nobs_type+1
cobtype(i_pm10_ob_type)  ="pm10 "; my_nobs_type=my_nobs_type+1
cobtype(i_gust_ob_type)  ="gust "; my_nobs_type=my_nobs_type+1
cobtype(i_vis_ob_type)   ="vis  "; my_nobs_type=my_nobs_type+1
cobtype(i_pblh_ob_type)  ="pblh "; my_nobs_type=my_nobs_type+1
cobtype(i_wspd10m_ob_type)  ="ws10m"; my_nobs_type=my_nobs_type+1
cobtype(i_td2m_ob_type)  ="td2m "; my_nobs_type=my_nobs_type+1
cobtype(i_mxtm_ob_type)  ="mxtm "; my_nobs_type=my_nobs_type+1
cobtype(i_mitm_ob_type)  ="mitm "; my_nobs_type=my_nobs_type+1
cobtype(i_pmsl_ob_type)  ="pmsl "; my_nobs_type=my_nobs_type+1
cobtype(i_howv_ob_type)  ="howv "; my_nobs_type=my_nobs_type+1
cobtype(i_tcamt_ob_type)  ="tcamt"; my_nobs_type=my_nobs_type+1
cobtype(i_lcbas_ob_type)  ="lcbas"; my_nobs_type=my_nobs_type+1
cobtype(i_cldch_ob_type)  ="cldch"; my_nobs_type=my_nobs_type+1
cobtype(i_uwnd10m_ob_type) ="u10m "; my_nobs_type=my_nobs_type+1
cobtype(i_vwnd10m_ob_type) ="v10m "; my_nobs_type=my_nobs_type+1
cobtype(i_swcp_ob_type)  ="swcp "; my_nobs_type=my_nobs_type+1
cobtype(i_lwcp_ob_type)  ="lwcp "; my_nobs_type=my_nobs_type+1

if(my_nobs_type/=nobs_type) then
   if (mype==0) then
      write(6,*)'init_fc_sens: inconsistent nobs_types, code needs update'
   endif
   call stop2(999)
endif
if (mype==0) then
   write(6,*)'init_fc_sens: lobsensincr,lobsensfc,lobsensjb=', &
                            lobsensincr,lobsensfc,lobsensjb
   write(6,*)'init_fc_sens: lobsensadj,lobsensmin,iobsconv=', &
                            lobsensadj,lobsensmin,iobsconv
   write(6,*)'init_fc_sens: lsensrecompute=',lsensrecompute
endif

call allocate_cv(fcsens)
fcsens=zero

if (lobsensadj.and.lobsensmin) then
   if (mype==0) then
      write(6,*)'init_fc_sens: unknown method',lobsensadj,lobsensmin
   endif
   call stop2(155)
end if

if (iobsconv>=2) then
   allocate(sensincr(nobs_bins,nobs_type,niter(jiter)))
else
   allocate(sensincr(nobs_bins,nobs_type,1))
endif
sensincr=zero

! Initialize fcsens
if (lobsensfc) then
   if (lobsensincr) then
      clfile='xhatsave.ZZZ'
      write(clfile(10:12),'(I3.3)') jiter
      call read_cv(fcsens,clfile)
      if (jiter>1) then
         clfile='xhatsave.ZZZ'
         write(clfile(10:12),'(I3.3)') jiter-1
         call allocate_cv(xwork)
         call read_cv(xwork,clfile)
         do ii=1,fcsens%lencv
            fcsens%values(ii) = fcsens%values(ii) - xwork%values(ii)
         end do
         call deallocate_cv(xwork)
      endif
   else
      if (jiter==miter) then
         if (lobsensjb) then
            clfile='xhatsave.ZZZ'
            write(clfile(10:12),'(I3.3)') miter
            call read_cv(fcsens,clfile)
         else
!           read and convert output of GCM adjoint
            allocate(fname(nsubwin))
            fname='NULL'
            if (tau_fcst>0) then
               fname='B'
            endif
            do ii=1,nsubwin
               call allocate_state(fcgrad(ii))
            end do
            call allocate_preds(zbias)
            zbias=zero
            call gsi_4dcoupler_getpert(fcgrad,nsubwin,'adm',fname)
            if (lsqrtb) then
               call control2model_ad(fcgrad,zbias,fcsens)
            else
               if (l_hyb_ens) then
                  do ii=1,ntlevs_ens
                     call allocate_state(eval(ii))
                  end do
                  eval(1)=fcgrad(1)
                  fcgrad(1)=zero
                  call ensctl2state_ad(eval,fcgrad(1),fcsens)
                  call control2state_ad(fcgrad,zbias,fcsens)
                  do ii=1,ntlevs_ens
                     call deallocate_state(eval(ii))
                  end do
                  if (tau_fcst>0) then
                     call destroy_hybens_localization_parameters
                     call destroy_ensemble
                     call create_ensemble
                     ! now load actual ens background
                     efsoi_afcst=.false.
                     orig_tau_fcst=tau_fcst
                     tau_fcst=-1
                     call load_ensemble
                     call hybens_localization_setup
                  endif
               else
                  call control2state_ad(fcgrad,zbias,fcsens)
               end if
            endif
            do ii=1,nsubwin
               call deallocate_state(fcgrad(ii))
            end do
            call deallocate_preds(zbias)
            deallocate(fname)
         endif
      else
!        read gradient from outer loop jiter+1
         clfile='fgsens.ZZZ'
         WRITE(clfile(8:10),'(I3.3)') jiter+1
         call read_cv(fcsens,clfile)
      endif
   endif
   zjx=dot_product(fcsens,fcsens)
   if (mype==0) write(6,888)'init_fc_sens: Norm fcsens=',sqrt(zjx)
endif
888 format(A,3(1X,ES25.18))


return
end subroutine init_fc_sens
! ------------------------------------------------------------------------------
subroutine efsoi_o2_update(sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    order2_fc_sens
!   prgmmr:      todling
!
! abstract: Read forecast sensitivity gradient
!
! program history log:
!   2007-06-26  tremolet - initial code
!   2009-08-07  lueken - added subprogram doc block
!   2010-05-27  todling - gsi_4dcoupler; remove dependence on GMAO specifics
!   2012-05-22  todling - update interface to getpert
!   2015-12-01  todling - add several obs-types that Pondeca forget to add here
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

implicit none
type(gsi_bundle),intent(inout)  :: sval(nobs_bins)

character(len=12) :: clfile
character(len=80),allocatable,dimension(:)::fname
type(gsi_bundle) :: fcgrad(nsubwin)
type(gsi_bundle) :: eval(ntlevs_ens)
type(gsi_bundle) :: mval(nsubwin)
type(predictors) :: zbias
type(gsi_bundle) :: rval(nobs_bins)
real(r_kind) :: zjx
integer(i_kind) :: ii

tau_fcst=orig_tau_fcst
if (tau_fcst<=0) return
if (efsoi_order/=2) return
if (.not.l_hyb_ens) return

if (mype==0) then
   print *, 'in order2_fc_sens ...'
endif
if(my_nobs_type/=nobs_type) then
   write(6,*)'order2_fc_sens: inconsistent nobs_types, code needs update'
   call stop2(999)
endif

if (lobsensadj.and.lobsensmin) then
   write(6,*)'order2_fc_sens: unknown method',lobsensadj,lobsensmin
   call stop2(155)
end if

! Zero out whatever is in fcsens
 fcsens=zero

! read and convert output of GCM adjoint
 allocate(fname(nsubwin))
 fname='A'
 do ii=1,nsubwin
    call allocate_state(fcgrad(ii))
 end do
 call allocate_preds(zbias)
 zbias=zero
 do ii=1,ntlevs_ens
    call allocate_state(eval(ii))
 end do
 do ii=1,nsubwin
    fcgrad(ii)=zero
 end do
 call gsi_4dcoupler_getpert(fcgrad,nsubwin,'adm',fname)
! Wipe out loaded ensemble members to allow reloading ensemble of analysis
 call destroy_ensemble
! Read in ens forecasts issues from "analysis" 
 efsoi_afcst=.true.
 call create_ensemble
 call load_ensemble
 do ii=1,ntlevs_ens
    call allocate_state(eval(ii))
 end do
 eval(1)=fcgrad(1)
 fcgrad(1)=zero
 call ensctl2state_ad(eval,fcgrad(1),fcsens)
 call control2state_ad(fcgrad,zbias,fcsens)
 do ii=1,ntlevs_ens
    call deallocate_state(eval(ii))
 end do
 do ii=1,nsubwin
    call deallocate_state(fcgrad(ii))
 end do
 call deallocate_preds(zbias)
 deallocate(fname)
! Wipe  outensemble from memory
 call destroy_ensemble

! Report magnitude of input vector
 zjx=dot_product(fcsens,fcsens)
 if (mype==0) write(6,888)'order2_fc_sens: Norm fcsens=',sqrt(zjx)
888 format(A,3(1X,ES25.18))

! Read in ensemble of analysis (these are EnKF, not GSI, analyses obviously)
 efsoi_ana=.true.
 tau_fcst=-1
 call create_ensemble
 call load_ensemble

! Allocate local variables
 do ii=1,nsubwin
    call allocate_state(mval(ii))
 end do
 do ii=1,ntlevs_ens
    call allocate_state(eval(ii))
 end do
 call allocate_preds(zbias)

! Convert from control variable to state space
 call control2state(fcsens,mval,zbias)

! Convert ensemble control variable to state space and update from previous value
 call ensctl2state(fcsens,mval(1),eval)
 do ii=1,ntlevs_ens
    call self_add(sval(ii),eval(ii))
 end do

! Release memory
call deallocate_preds(zbias)
do ii=1,ntlevs_ens
   call deallocate_state(eval(ii))
end do
do ii=1,nsubwin
   call deallocate_state(mval(ii))
end do

return
end subroutine efsoi_o2_update
subroutine save_fc_sens
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    save_fc_sens
!   prgmmr:      tremolet
!
! abstract: Compute and save forecast sensitivity to observations
!
! program history log:
!   2007-06-26  tremolet - initial code
!   2009-08-07  lueken - added subprogram doc block
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
implicit none

real(r_kind) :: zz
integer(i_kind) :: ii,jj,kk

! Save statistics
if (mype==0) then

!  Full stats
   do jj=1,nobs_type
      write(6,'(A,2X,I3,2X,A)')'Obs types:',jj,cobtype(jj)
   enddo
   write(6,'(A,2X,I4)')'Obs bins:',nobs_bins
   write(6,*)'Obs Count Begin'
   if (.not.allocated(obscounts)) then
      write(6,*)'save_fc_sens: obscounts not allocated'
      call stop2(156)
   end if
   do jj=1,nobs_type
      write(6,'((1X,I12))')(obscounts(jj,ii),ii=1,nobs_bins)
   enddo
   write(6,*)'Obs Count End'

   write(6,*)'Obs Impact Begin'
   do kk=1,SIZE(sensincr,3)
      if (SIZE(sensincr,3)==1) then
         write(6,'(A,I4)')'Obs Impact iteration= ',niter(jiter)
      else
         write(6,'(A,I4)')'Obs Impact iteration= ',kk
      endif
      do jj=1,nobs_type
         write(6,'((1X,ES12.5))')(sensincr(ii,jj,kk),ii=1,nobs_bins)
      enddo
   enddo
   write(6,*)'Obs Impact End'

   kk=SIZE(sensincr,3)
!  Summary by obs type
   do jj=1,nobs_type
      zz=zero
      do ii=1,nobs_bins
         zz=zz+sensincr(ii,jj,kk)
      enddo
      if (zz/=zero) write(6,'(A,2X,A3,2X,ES12.5)')'Obs Impact type',cobtype(jj),zz
   enddo

!  Summary by obs bins
   do ii=1,nobs_bins
      zz=zero
      do jj=1,nobs_type
         zz=zz+sensincr(ii,jj,kk)
      enddo
      if (zz/=zero) write(6,'(A,2X,I3,2X,ES12.5)')'Obs Impact bin',ii,zz
   enddo

endif

deallocate(sensincr)
call deallocate_cv(fcsens)

return
end subroutine save_fc_sens
! ------------------------------------------------------------------------------
real(r_kind) function dot_prod_obs()
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_fc_sens
!   prgmmr:      tremolet
!
! abstract: Computes scalar product in observation space
!           (based on evaljo)
!
! program history log:
!   2007-06-27  tremolet
!   2009-01-18  todling - carry summations in quad precision
!   2009-08-07  lueken - added subprogram doc block
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
use obsmod, only: obs_diag
use obsmod, only: obsdiags
implicit none

integer(i_kind) :: ii,jj,ij,it
real(r_quad)    :: zzz
real(r_quad)    :: zprods(nobs_type*nobs_bins)
type(obs_diag),pointer:: obsptr
! ----------------------------------------------------------

zprods(:)=zero_quad

ij=0
do ii=1,nobs_bins
   do jj=1,nobs_type
      ij=ij+1

      obsptr => obsdiags(jj,ii)%head
      do while (associated(obsptr))
         if (obsptr%luse.and.obsptr%muse(jiter)) then
            zprods(ij) = zprods(ij) + obsptr%nldepart(jiter) * obsptr%obssen(jiter)
         endif
         obsptr => obsptr%next
      enddo

   enddo
enddo

! Gather contributions
call mpl_allreduce(nobs_type*nobs_bins,qpvals=zprods)

! Save intermediate values
it=-1
if (iobsconv>=2) then
   if (iter>=1.and.iter<=niter(jiter)) it=iter
else
   it=1
endif

if (it>0) then
   ij=0
   do ii=1,nobs_bins
      do jj=1,nobs_type
         ij=ij+1
         sensincr(ii,jj,it)=zprods(ij)
      enddo
   enddo
endif

! Sum
zzz=zero_quad

ij=0
do ii=1,nobs_bins
   do jj=1,nobs_type
      ij=ij+1
      zzz=zzz+zprods(ij)
   enddo
enddo

dot_prod_obs=zzz

return
end function dot_prod_obs
! ------------------------------------------------------------------------------
end module obs_sensitivity
