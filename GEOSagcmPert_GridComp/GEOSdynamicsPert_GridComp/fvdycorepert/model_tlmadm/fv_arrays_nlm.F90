!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
module fv_arrays_nlm_mod

#include <fms_platform.h>

implicit none
public

logical :: fv_timing_onoff

!Compiler definitions that are regular if statements for the tlm/adm, should still match how the code is compiled
type fpp_type
  logical :: FPP_MAPL_MODE   = .true.    !GMAO Mode
  logical :: FPP_OVERLOAD_R4 = .false.
end type fpp_type

type(fpp_type) :: fpp

type fv_flags_pert_type

! Traj = pert, true is faster
   logical :: split_hord   = .true. !Force traj to use pert options if false

! TLM/ADM options for horizontal transport
   integer :: hord_mt_pert = 2      ! Linear options:
   integer :: hord_vt_pert = 2      !   1 (all) first order scheme
   integer :: hord_tm_pert = 2      !   6 (all) quasi fifth order linear scheme
   integer :: hord_dp_pert = 2      !  -5 (vt,tm,dp,tr) linear scheme based on PPM 4th order interpolation
   integer :: hord_tr_pert = 2      ! 333 (vt,tm,dp,tr) linear third order cheme

! Number of sponge layers for the perturbations
   integer :: n_sponge_pert = 9

!Perturbation damping in the sponge
   logical :: hord_ks_pert    = .true. !Use the below for the trajectory in the sponge layer
   integer :: hord_mt_ks_pert = 1
   integer :: hord_vt_ks_pert = 1
   integer :: hord_tm_ks_pert = 1
   integer :: hord_dp_ks_pert = 1
   integer :: hord_tr_ks_pert = 1

!Trajectory damping in the sponge
   logical :: hord_ks_traj    = .true. !Use the below for the trajectory in the sponge layer
   integer :: hord_mt_ks_traj = 1
   integer :: hord_vt_ks_traj = 1
   integer :: hord_tm_ks_traj = 1
   integer :: hord_dp_ks_traj = 1
   integer :: hord_tr_ks_traj = 1

! TLM/ADM options for vertical remapping
   logical :: split_kord   = .true. !Force traj to use pert options if false
   integer :: kord_mt_pert = 17     ! |kord|>16 is linear remapping, not shape preserving or
   integer :: kord_wz_pert = 17     ! positive definite, anything else not recommended
   integer :: kord_tm_pert = 17     ! for perturbation quantities
   integer :: kord_tr_pert = 17 

! TLM/ADM options for damping
   logical :: split_damp = .true.   !Force traj to use pert options if false
   integer ::  nord_pert = 1 
   real    :: dddmp_pert = 0.2
   real    :: d2_bg_pert = 0.015
   real    :: d4_bg_pert = 0.150
   logical :: do_vort_damp_pert = .true.
   real    :: vtdm4_pert = 0.0005
   real    :: d2_bg_k1_pert = 4.         ! factor for d2_bg (k=1)
   real    :: d2_bg_k2_pert = 2.         ! factor for d2_bg (k=2)
   real    :: d2_bg_ks_pert = 2.         ! factor for d2_bg (k=1)

! TLM/ADM options for tracer damping
   logical :: split_damp_tr = .true. !Force traj to use pert options if false
   integer :: nord_tr_pert=0         ! Damping order
   real    :: trdm2_pert = 0.0       ! coefficient for damping

end type fv_flags_pert_type

type fv_atmos_pert_type

   type(fv_flags_pert_type) :: flagstruct

   real, allocatable ::     up(:,:,:)   ! D grid zonal wind (m/s)
   real, allocatable ::     vp(:,:,:)   ! D grid meridional wind (m/s)
   real, allocatable ::    ptp(:,:,:)   ! temperature (K)
   real, allocatable ::  delpp(:,:,:)   ! pressure thickness (pascal)
   real, allocatable ::     qp(:,:,:,:) ! specific humidity and constituents

   real, allocatable ::     wp(:,:,:)   ! cell center vertical wind (m/s)
   real, allocatable ::  delzp(:,:,:)   ! layer thickness (meters)
   real, allocatable ::   ze0p(:,:,:)   ! height at layer edges for remapping
   real, allocatable :: q_conp(:,:,:)   ! total condensates

   real, allocatable ::    psp(:,:)     ! Surface pressure (pascal)
   real, allocatable ::    pep(:,:,:)   ! edge pressure (pascal)
   real, allocatable ::    pkp(:,:,:)   ! pe**cappa
   real, allocatable ::  pelnp(:,:,:)   ! ln(pe)
   real, allocatable ::   pkzp(:,:,:)   ! finite-volume mean pk

   real, allocatable ::  omgap(:,:,:)   ! Vertical pressure velocity (pa/s)

   real, allocatable ::    uap(:,:,:)   ! (ua, va) are mostly used as the A grid winds
   real, allocatable ::    vap(:,:,:)     
   real, allocatable ::    ucp(:,:,:)   ! (uc, vc) are mostly used as the C grid winds
   real, allocatable ::    vcp(:,:,:)     

   real, allocatable ::   mfxp(:,:,:)   ! Mass fluxes
   real, allocatable ::   mfyp(:,:,:)  

   real, allocatable ::    cxp(:,:,:)  
   real, allocatable ::    cyp(:,:,:)  

end type fv_atmos_pert_type

!---- version number -----
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

contains

subroutine allocate_fv_atmos_pert_type(AtmP, isd, ied, jsd, jed, is, ie, js, je, npz, ncnst)

  implicit none
  type(fv_atmos_pert_type), intent(INOUT), target :: AtmP
  integer, intent(IN) :: isd, ied, jsd, jed, is, ie, js, je, npz, ncnst

  ! Allocate the perturbation variables
  allocate (    AtmP%up(isd:ied  ,jsd:jed+1,npz) )
  allocate (    AtmP%vp(isd:ied+1,jsd:jed  ,npz) )
  allocate (   AtmP%ptp(isd:ied  ,jsd:jed  ,npz) )
  allocate ( AtmP%delpp(isd:ied  ,jsd:jed  ,npz) )
  allocate (    AtmP%qp(isd:ied  ,jsd:jed  ,npz, ncnst) )
  allocate (    AtmP%wp(isd:ied, jsd:jed  ,npz  ) )
  allocate ( AtmP%delzp(isd:ied, jsd:jed  ,npz) )
  allocate (  AtmP%ze0p(is:ie  , js:je    ,npz+1) )
  allocate (AtmP%q_conp(isd:ied,jsd:jed,1:npz) ) 
  allocate (   AtmP%psp(isd:ied  ,jsd:jed) )
  allocate (   AtmP%pep(is-1:ie+1, npz+1,js-1:je+1) )
  allocate (   AtmP%pkp(is:ie    ,js:je  , npz+1) )
  allocate ( AtmP%pelnp(is:ie,npz+1,js:je) )
  allocate (  AtmP%pkzp(is:ie,js:je,npz) )
  allocate ( AtmP%omgap(isd:ied  ,jsd:jed  ,npz) )
  allocate (   AtmP%uap(isd:ied  ,jsd:jed  ,npz) )
  allocate (   AtmP%vap(isd:ied  ,jsd:jed  ,npz) )
  allocate (   AtmP%ucp(isd:ied+1,jsd:jed  ,npz) )
  allocate (   AtmP%vcp(isd:ied  ,jsd:jed+1,npz) )
  allocate (  AtmP%mfxp(is:ie+1, js:je,  npz) )
  allocate (  AtmP%mfyp(is:ie  , js:je+1,npz) )
  allocate (   AtmP%cxp(is:ie+1, jsd:jed, npz) )
  allocate (   AtmP%cyp(isd:ied ,js:je+1, npz) )

     AtmP%up = 0.0
     AtmP%vp = 0.0
    AtmP%ptp = 0.0
  AtmP%delpp = 0.0
     AtmP%qp = 0.0
     AtmP%wp = 0.0
  AtmP%delzp = 0.0
   AtmP%ze0p = 0.0
 AtmP%q_conp = 0.0
    AtmP%psp = 0.0
    AtmP%pep = 0.0
    AtmP%pkp = 0.0
  AtmP%pelnp = 0.0
   AtmP%pkzp = 0.0
  AtmP%omgap = 0.0
    AtmP%uap = 0.0
    AtmP%vap = 0.0
    AtmP%ucp = 0.0
    AtmP%vcp = 0.0
   AtmP%mfxp = 0.0
   AtmP%mfyp = 0.0
    AtmP%cxp = 0.0
    AtmP%cyp = 0.0

end subroutine allocate_fv_atmos_pert_type

subroutine deallocate_fv_atmos_pert_type(AtmP)

  implicit none
  type(fv_atmos_pert_type), intent(INOUT) :: AtmP

  deallocate ( AtmP%up    )
  deallocate ( AtmP%vp    )
  deallocate ( AtmP%ptp   )
  deallocate ( AtmP%delpp )
  deallocate ( AtmP%qp    )
  deallocate ( AtmP%wp    )
  deallocate ( AtmP%delzp )
  deallocate ( AtmP%ze0p  )
  deallocate ( AtmP%psp   )
  deallocate ( AtmP%pep   )
  deallocate ( AtmP%pkp   )
  deallocate ( AtmP%pelnp )
  deallocate ( AtmP%pkzp  )
  deallocate ( AtmP%omgap )
  deallocate ( AtmP%uap   )
  deallocate ( AtmP%vap   )
  deallocate ( AtmP%ucp   )
  deallocate ( AtmP%vcp   )
  deallocate ( AtmP%mfxp  )
  deallocate ( AtmP%mfyp  )
  deallocate ( AtmP%cxp   )
  deallocate ( AtmP%cyp   )

end subroutine deallocate_fv_atmos_pert_type

end module fv_arrays_nlm_mod
