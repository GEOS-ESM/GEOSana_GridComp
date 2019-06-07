!
! !REVISION HISTORY:
!
!  14May2007 Todling Introduced dyn_prog; global change.
!  17May2007 Todling Largely revampped from original code:
!                    - turned into model
!                    - interfaces 1 and 2
!
!  31May2007 Todling Add hooks to handle g5 perturbations
!
!
module m_model_tad

!==============================================
! referencing used modules
!==============================================
use precision
use prognostics
use prognostics_q
use stepon, only : nymd,nhms,nouter,ninner,nstep,stepon_do,stepon_set
use stepon, only : stepon_tape_rec
use fvcore, only : dynpkg_n2,dynpkg_nsplit
use mod_comm, only : gid, mp_bcst_n_real, mpi_bcst_n_real_ad
use stepon_tad, only : stepon_do_tad, stepon_domd
use m_zeit, only : zeit_ci
use m_zeit, only : zeit_co

use c_sw_tad_store, only : c_sw_tape1_a6_5h,c_sw_tape1_c_sw,c_sw_tape1_crx_9h,c_sw_tape1_tm2_7h,c_sw_tape1_v2_6h,c_sw_tape1_xfx_1h,&
&c_sw_tape1_yfx_2h,c_sw_tape2_a6_1h,c_sw_tape2_c_sw,c_sw_tape2_slope_2h
use cd_core_do_tad_store, only : cd_core_tape_cd_core_do,cd_core_tape_delpf_3h,cd_core_tape_pt_16h,cd_core_tape_pt_5h,&
&cd_core_tape_u_12h,cd_core_tape_u_1h,cd_core_tape_uc_14h,cd_core_tape_v_13h,cd_core_tape_v_2h,cd_core_tape_vc_15h
use d_sw_tad_store, only : d_sw_tape1_d_sw,d_sw_tape1_delpf_1h,d_sw_tape1_uc_7h,d_sw_tape1_xfx_2h,d_sw_tape1_yfx_3h,&
&d_sw_tapej_a6_1h,d_sw_tapej_d_sw,d_sw_tapej_slope_2h
use fvcore_do_tad_store, only : dynpkg_1_fvcore_do,dynpkg_1_pe_4h,dynpkg_1_ps_5h,dynpkg_1_q_3h,dynpkg_n2_cx_3h,dynpkg_n2_cy_4h,&
&dynpkg_n2_dp0_1h,dynpkg_n2_fvcore_do,dynpkg_n2_mfx_5h,dynpkg_n2_mfy_6h,dynpkg_n2_q_2h
use stepon_do_tad_store, only : inner_omega_8h,inner_pe_14h,inner_stepon_do
use te_map_tad_store, only : te_map_tape_dz_8h,te_map_tape_te_9h,te_map_tape_te_map,te_map_tape_u_6h,te_map_tape_v_7h
use tpcc_tad_store, only : tpcc1_tape_a6_1h,tpcc1_tape_dm_2h,tpcc1_tape_tpcc,tpcc2_tape_a6_1h,tpcc2_tape_dm_2h,tpcc2_tape_tpcc

implicit none

PRIVATE

PUBLIC qmodel_tad
PUBLIC initial_tad
PUBLIC final_tad

interface qmodel_tad; module procedure &
          model_ad1_,&
          model_ad2_
end interface qmodel_tad

logical, save :: first_  = .true.

contains

subroutine model_ad1_( pert_ad, nymdi, nhmsi, ntsteps )

implicit none

!==============================================
! declare arguments
!==============================================
type(dynq) :: pert_ad
integer, optional, intent(in) :: nymdi
integer, optional, intent(in) :: nhmsi
integer, optional, intent(in) :: ntsteps

type(dynq) :: xpert
type(dynq) :: ypert
integer i,nymds,nhmss,nouters,nstepsv
logical, save :: setup = .true.
logical reset

!----------------------------------------------
! RESET TIME IN ADM
!----------------------------------------------
reset = present(nymdi) .and. present(nhmsi) .and. present(ntsteps)
if ( reset ) then
     nstepsv= ninner+(nouter-1)*ninner
     nymds  = nymd ; nhmss=nhms ; nouters= nouter
     nymd   = nymdi; nhms =nhmsi; nouter = ntsteps
     setup  = .false.
endif

!----------------------------------------------
! DEFFINE AUXILIAR VECTORS
!----------------------------------------------
call prognostics_q_initial ( xpert )  
call prognostics_q_initial ( ypert )  

!----------------------------------------------
! ROUTINE BODY
!----------------------------------------------

call prognostics_q_dup ( pert_ad, ypert )  

call model_ad2_ ( xpert, ypert, setup=setup )

call prognostics_q_dup ( xpert, pert_ad )  

!----------------------------------------------
! CLEAN UP
!----------------------------------------------
call prognostics_q_final ( ypert )  
call prognostics_q_final ( xpert )  

!----------------------------------------------
! RESET TIME BACK TO ORIGINAL IN TLM
!----------------------------------------------
if ( reset ) then
     nstep = nstepsv
     nymd  = nymds; nhms =nhmss; nouter =nouters
     setup = .true.
endif

end subroutine model_ad1_

subroutine model_ad2_( xpert, ypert, setup )

implicit none

!==============================================
! declare arguments
!==============================================
type(dynq) :: ypert  ! input  perturbation
type(dynq) :: xpert  ! output perturbation
logical,optional,intent(in) :: setup

type(dyn_prog) :: myprog
real(r8),allocatable::y_ad(:)
logical, save :: setup_  = .true.
integer i,ndim

if (present(setup)) then
    setup_ = setup
endif

!----------------------------------------------
! RESET ADJOINT AND TRAJECTORY VARIABLES
!----------------------------------------------
call prognostics_initial ( myprog )
if(setup_) call initial_tad
call prognostics_q_zero ( xpert )

!----------------------------------------------
! TAPE RECORDING
!----------------------------------------------
stepon_tape_rec = 0
!_RTif(setup_) call stepon_set ( myprog )
call stepon_set ( myprog )
call stepon_domd ( first_, myprog )

!----------------------------------------------
! ADJOINT COMPUTATIONS
!----------------------------------------------
   call zeit_ci('stpqadm')
call stepon_do_tad ( myprog, xpert, ypert )
   call zeit_co('stpqadm')

!----------------------------------------------
! FINALIZE ADJOINT AND TRAJECTORY VARIABLES
!----------------------------------------------
if ( setup_ ) then
     call final_tad
endif
first_=.false.
call prognostics_final ( myprog )

end subroutine model_ad2_

subroutine initial_tad
!******************************************************************
!******************************************************************
!** This routine was generated by Automatic differentiation.     **
!** FastOpt: Transformation of Algorithm in Fortran, TAF 1.6.1   **
!******************************************************************
!******************************************************************
!==============================================
! referencing used modules
!==============================================

!==============================================
! all entries are defined explicitly
!==============================================
implicit none

!----------------------------------------------
! OPEN TAPE c_sw_tape1
!----------------------------------------------
c_sw_tape1_c_sw = ninner*dynpkg_nsplit*dynpkg_n2*nl*1

!----------------------------------------------
! OPEN TAPE c_sw_tape2
!----------------------------------------------
c_sw_tape2_c_sw = ninner*dynpkg_nsplit*dynpkg_n2*nl*(jnp-1)

!----------------------------------------------
! OPEN TAPE cd_core_tape
!----------------------------------------------
cd_core_tape_cd_core_do = ninner*dynpkg_nsplit*dynpkg_n2

!----------------------------------------------
! OPEN TAPE d_sw_tape1
!----------------------------------------------
d_sw_tape1_d_sw = ninner*dynpkg_nsplit*dynpkg_n2*nl*1

!----------------------------------------------
! OPEN TAPE d_sw_tapej
!----------------------------------------------
d_sw_tapej_d_sw = ninner*dynpkg_nsplit*dynpkg_n2*nl*(jnp-1)

!----------------------------------------------
! OPEN TAPE dynpkg_1
!----------------------------------------------
dynpkg_1_fvcore_do = ninner

!----------------------------------------------
! OPEN TAPE dynpkg_n2
!----------------------------------------------
dynpkg_n2_fvcore_do = ninner*dynpkg_n2

!----------------------------------------------
! OPEN TAPE inner
!----------------------------------------------
inner_stepon_do = ninner

!----------------------------------------------
! OPEN TAPE te_map_tape
!----------------------------------------------
te_map_tape_te_map = ninner

!----------------------------------------------
! OPEN TAPE tpcc1_tape
!----------------------------------------------
tpcc1_tape_tpcc = ninner*dynpkg_nsplit*dynpkg_n2*nl*(jnp-1)

!----------------------------------------------
! OPEN TAPE tpcc2_tape
!----------------------------------------------
tpcc2_tape_tpcc = ninner*dynpkg_nsplit*dynpkg_n2*nl*(jnp-1)


end subroutine initial_tad


subroutine final_tad
!******************************************************************
!******************************************************************
!** This routine was generated by Automatic differentiation.     **
!** FastOpt: Transformation of Algorithm in Fortran, TAF 1.6.1   **
!******************************************************************
!******************************************************************
!==============================================
! referencing used modules
!==============================================

!==============================================
! all entries are defined explicitly
!==============================================
implicit none

!----------------------------------------------
! CLOSE TAPE c_sw_tape1
!----------------------------------------------
if (allocated(c_sw_tape1_xfx_1h)) then
  deallocate( c_sw_tape1_xfx_1h )
endif
if (allocated(c_sw_tape1_yfx_2h)) then
  deallocate( c_sw_tape1_yfx_2h )
endif
if (allocated(c_sw_tape1_a6_5h)) then
  deallocate( c_sw_tape1_a6_5h )
endif
if (allocated(c_sw_tape1_v2_6h)) then
  deallocate( c_sw_tape1_v2_6h )
endif
if (allocated(c_sw_tape1_tm2_7h)) then
  deallocate( c_sw_tape1_tm2_7h )
endif
if (allocated(c_sw_tape1_crx_9h)) then
  deallocate( c_sw_tape1_crx_9h )
endif

!----------------------------------------------
! CLOSE TAPE c_sw_tape2
!----------------------------------------------
if (allocated(c_sw_tape2_a6_1h)) then
  deallocate( c_sw_tape2_a6_1h )
endif
if (allocated(c_sw_tape2_slope_2h)) then
  deallocate( c_sw_tape2_slope_2h )
endif

!----------------------------------------------
! CLOSE TAPE cd_core_tape
!----------------------------------------------
if (allocated(cd_core_tape_u_1h)) then
  deallocate( cd_core_tape_u_1h )
endif
if (allocated(cd_core_tape_v_2h)) then
  deallocate( cd_core_tape_v_2h )
endif
if (allocated(cd_core_tape_delpf_3h)) then
  deallocate( cd_core_tape_delpf_3h )
endif
if (allocated(cd_core_tape_pt_5h)) then
  deallocate( cd_core_tape_pt_5h )
endif
if (allocated(cd_core_tape_u_12h)) then
  deallocate( cd_core_tape_u_12h )
endif
if (allocated(cd_core_tape_v_13h)) then
  deallocate( cd_core_tape_v_13h )
endif
if (allocated(cd_core_tape_uc_14h)) then
  deallocate( cd_core_tape_uc_14h )
endif
if (allocated(cd_core_tape_vc_15h)) then
  deallocate( cd_core_tape_vc_15h )
endif
if (allocated(cd_core_tape_pt_16h)) then
  deallocate( cd_core_tape_pt_16h )
endif

!----------------------------------------------
! CLOSE TAPE d_sw_tape1
!----------------------------------------------
if (allocated(d_sw_tape1_delpf_1h)) then
  deallocate( d_sw_tape1_delpf_1h )
endif
if (allocated(d_sw_tape1_xfx_2h)) then
  deallocate( d_sw_tape1_xfx_2h )
endif
if (allocated(d_sw_tape1_yfx_3h)) then
  deallocate( d_sw_tape1_yfx_3h )
endif
if (allocated(d_sw_tape1_uc_7h)) then
  deallocate( d_sw_tape1_uc_7h )
endif

!----------------------------------------------
! CLOSE TAPE d_sw_tapej
!----------------------------------------------
if (allocated(d_sw_tapej_a6_1h)) then
  deallocate( d_sw_tapej_a6_1h )
endif
if (allocated(d_sw_tapej_slope_2h)) then
  deallocate( d_sw_tapej_slope_2h )
endif

!----------------------------------------------
! CLOSE TAPE dynpkg_1
!----------------------------------------------
if (allocated(dynpkg_1_q_3h)) then
  deallocate( dynpkg_1_q_3h )
endif
if (allocated(dynpkg_1_pe_4h)) then
  deallocate( dynpkg_1_pe_4h )
endif
if (allocated(dynpkg_1_ps_5h)) then
  deallocate( dynpkg_1_ps_5h )
endif

!----------------------------------------------
! CLOSE TAPE dynpkg_n2
!----------------------------------------------
if (allocated(dynpkg_n2_dp0_1h)) then
  deallocate( dynpkg_n2_dp0_1h )
endif
if (allocated(dynpkg_n2_q_2h)) then
  deallocate( dynpkg_n2_q_2h )
endif
if (allocated(dynpkg_n2_cx_3h)) then
  deallocate( dynpkg_n2_cx_3h )
endif
if (allocated(dynpkg_n2_cy_4h)) then
  deallocate( dynpkg_n2_cy_4h )
endif
if (allocated(dynpkg_n2_mfx_5h)) then
  deallocate( dynpkg_n2_mfx_5h )
endif
if (allocated(dynpkg_n2_mfy_6h)) then
  deallocate( dynpkg_n2_mfy_6h )
endif

!----------------------------------------------
! CLOSE TAPE inner
!----------------------------------------------
if (allocated(inner_omega_8h)) then
  deallocate( inner_omega_8h )
endif
if (allocated(inner_pe_14h)) then
  deallocate( inner_pe_14h )
endif

!----------------------------------------------
! CLOSE TAPE te_map_tape
!----------------------------------------------
if (allocated(te_map_tape_u_6h)) then
  deallocate( te_map_tape_u_6h )
endif
if (allocated(te_map_tape_v_7h)) then
  deallocate( te_map_tape_v_7h )
endif
if (allocated(te_map_tape_dz_8h)) then
  deallocate( te_map_tape_dz_8h )
endif
if (allocated(te_map_tape_te_9h)) then
  deallocate( te_map_tape_te_9h )
endif

!----------------------------------------------
! CLOSE TAPE tpcc1_tape
!----------------------------------------------
if (allocated(tpcc1_tape_a6_1h)) then
  deallocate( tpcc1_tape_a6_1h )
endif
if (allocated(tpcc1_tape_dm_2h)) then
  deallocate( tpcc1_tape_dm_2h )
endif

!----------------------------------------------
! CLOSE TAPE tpcc2_tape
!----------------------------------------------
if (allocated(tpcc2_tape_a6_1h)) then
  deallocate( tpcc2_tape_a6_1h )
endif
if (allocated(tpcc2_tape_dm_2h)) then
  deallocate( tpcc2_tape_dm_2h )
endif

end subroutine final_tad

end module m_model_tad
