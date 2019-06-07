!  $Id$

#include "MAPL_Generic.h"
#define MAPL_FieldBundleGetPointer ESMFL_BundleGetPointerToData

!-----------------------------------------------------------------------
!              ESMA - Earth System Modeling Applications
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: DynCorePert_GridCompMod --- Linearized version of FV_StateMod
! !INTERFACE:

   module GEOS_DynCorePertGridCompMod

! ESMF/MAPL uses

   use ESMF                ! ESMF base class
   use MAPL_Mod            ! GEOS base class

! GEOS_PertShared uses
   use GEOS_PertSharedMod

! fms uses
   use fms_mod,              only: fms_init, set_domain, nullify_domain
   use mpp_domains_mod,      only: mpp_update_domains, DGRID_NE, mpp_get_boundary, mpp_get_boundary_ad

! fv_core uses
   use fv_mp_mod,            only: is_master
   use fv_control_mod,       only: fv_init1, fv_init2
   use fv_control_nlm_mod,   only: fv_init_pert, fv_end_pert
   use fv_arrays_mod ,       only: FVPRC
   use fv_arrays_mod ,       only: fv_atmos_type
   use fv_arrays_nlm_mod ,   only: fv_atmos_pert_type, fv_timing_onoff
   use fv_dynamics_mod,      only: fv_dynamics
   use fv_dynamics_tlm_mod,  only: fv_dynamics_nlm => fv_dynamics
   use fv_dynamics_tlm_mod,  only: fv_dynamics_tlm
   use fv_dynamics_adm_mod,  only: fv_dynamics_fwd, fv_dynamics_bwd
   use fv_diagnostics_mod,   only: prt_minmax
   use fv_timing_mod,        only: timing_on, timing_off
   use fv_sg_tlm_mod,        only: fv_subgrid_z_tlm
   use fv_sg_adm_mod,        only: fv_subgrid_z_fwd, fv_subgrid_z_bwd
   use fv_pressures_tlm_mod, only: compute_pressures_tlm
   use fv_pressures_adm_mod, only: compute_pressures_fwd, compute_pressures_bwd

! Tapenade uses
   use tapenade_iter,        only: cp_iter_controls, cp_iter, cp_mod_ini, cp_mod_mid, cp_mod_end, &
                                   pushrealarray, poprealarray

! !PUBLIC MEMBER FUNCTIONS:

  implicit none
  private

  public SetServices
  public FV_Atm

! !DESCRIPTION: This module implements the TLM and adjoint model operators.
!               The components natural grid is the grid of the analysis, from
!               which we assume the operators are being called. Currently, this
!               grid is a Lat-Lon grid and the TLM and adjoint operate on a
!               Cubed-Sphere grid.    

!EOP
!----------------------------------------------------------------------
!BOC

  type(fv_atmos_type)     , allocatable, save :: FV_Atm(:)
  type(fv_atmos_pert_type), allocatable, save :: FV_AtmP(:)

  logical, allocatable, save             :: grids_on_this_pe(:)

! Wrapper for extracting internal state
! -------------------------------------

  integer :: proc_id
  integer :: cp_dyn_ind

  !Constants
  real(FVPRC) :: pi
  real(FVPRC) :: omega    ! angular velocity of earth's rotation
  real(FVPRC) :: cp       ! heat capacity of air at constant pressure
  real(FVPRC) :: radius   ! radius of the earth (m)
  real(FVPRC) :: rgas     ! Gas constant of the air
  real(FVPRC) :: rvap     ! Gas constant of vapor
  real(FVPRC) :: kappa    ! kappa
  real(FVPRC) :: grav     ! Gravity
  real(FVPRC) :: hlv      ! latent heat of evaporation
  real(FVPRC) :: zvir     ! RWV/RAIR-1
  real(FVPRC) :: p00

  integer             :: CustomNSplit
  real(FVPRC)         :: DT

  !Variables for Test Cases
  logical, save :: firststep = .true.
  real(FVPRC) :: lambda_c, theta_c, rr, ptop_r
  real(FVPRC), parameter   :: dry_mass = 98290., deg2rad = 0.01745329251
  real(FVPRC), allocatable :: lambda(:,:), theta(:,:)
  real, pointer :: LATS (:,:), LONS (:,:)
  real(FVPRC), save, allocatable, dimension(:,:,:) :: u_r, v_r, w_r, delz_r, pt_r, delp_r, pe_r, pk_r, peln_r, pkz_r
  real(FVPRC), save, allocatable, dimension(:,:,:) :: omga_r, ua_r, va_r, uc_r, vc_r, mfx_r, mfy_r, cx_r, cy_r
  real(FVPRC), save, allocatable, dimension(:,:,:,:) :: q_r
  real(FVPRC), save, allocatable, dimension(:,:) :: ps_r, phis_r
  real(FVPRC), save, allocatable, dimension(:) :: ak_r, bk_r
  real(FVPRC), save, allocatable, dimension(:,:,:) :: u_p, v_p, w_p, delz_p, pt_p, delp_p, pe_p, pk_p, peln_p, pkz_p
  real(FVPRC), save, allocatable, dimension(:,:,:) :: omga_p, ua_p, va_p, uc_p, vc_p, mfx_p, mfy_p, cx_p, cy_p
  real(FVPRC), save, allocatable, dimension(:,:,:,:) :: q_p
  real(FVPRC), save, allocatable, dimension(:,:) :: ps_p

contains

!----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices

! !DESCRIPTION:  SetServices registers Initialize, Run, and Finalize
!   methods for FV. Two stages of the FV run method are registered. The
!   first one does TLM calculations, and the second does the adjoint.
!   SetServices also creates a private internal state in which we
!   keep the surface geopotential and the ak, bk.

! !INTERFACE:

   subroutine SetServices ( gc, rc )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! gridded component
    integer, intent(out), optional     :: rc     ! return code
   

! !DESCRIPTION: Set services (register) for the FVCAM Dynamical Core
!               Grid Component.
!         
!EOP         
!----------------------------------------------------------------------
  
    integer                          :: status
    character(len=ESMF_MAXSTR)       :: IAm
    character(len=ESMF_MAXSTR)       :: COMP_NAME

    type (MAPL_MetaComp), pointer    :: MAPL
    type (ESMF_VM)                   :: VM
    integer                          :: p_split=1
    integer                          :: ks
    integer                          :: ndt
    integer                          :: comm
    integer                          :: tmp
    integer                          :: nx, ny, npxa, npya, npza

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "SetServices"

!BOP

! !IROUTINE: State Descriptions

! !DESCRIPTION: The component uses all three states (Import, Export
!  and Internal), in addition to a Private (non-ESMF) Internal state.
!  All three are managed by MAPL. 
!
!  The Import State contains perturbation increments to be added
!  by either of the run phases. If these are not to be used, the
!  Import state is not registered to save space.

!
!EOP   

    !FV_Setup
    !--------
    if ( ObserverMode == 0) then

       call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
       VERIFY_(STATUS)

       call ESMF_VMGetCurrent(VM, rc=STATUS)
       VERIFY_(STATUS)

       call MAPL_MemUtilsWrite(VM, trim(IAm), RC=STATUS )
       VERIFY_(STATUS)

       !Start up FMS
       call MAPL_TimerOn(MAPL,"--FMS_INIT")
       call ESMF_VMGet(VM,mpiCommunicator=comm,localPet=proc_id,rc=status)
       VERIFY_(STATUS)
       call fms_init(comm)
       call MAPL_TimerOff(MAPL,"--FMS_INIT")
       call MAPL_MemUtilsWrite(VM, 'FV_StateMod: FMS_INIT', RC=STATUS )
       VERIFY_(STATUS)

       !Model time step
       call MAPL_GetResource( MAPL, ndt, 'RUN_DT:', default=0, RC=STATUS )
       VERIFY_(STATUS)
       DT = real(ndt,FVPRC)

       !Start up FV                   
       call MAPL_TimerOn(MAPL,"--FV_INIT")
       call fv_init1(FV_Atm, DT, grids_on_this_pe, p_split)
       call MAPL_TimerOff(MAPL,"--FV_INIT")
       call MAPL_MemUtilsWrite(VM, 'FV_StateMod: FV_INIT', RC=STATUS )
       VERIFY_(STATUS)

       !FV grid dimensions setup from MAPL
       call MAPL_GetResource( MAPL, npxa, 'AGCM_IM:', default=180 , RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource( MAPL, npya, 'AGCM_JM:', default=1080, RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource( MAPL, npza, 'AGCM_LM:', default=72  , RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource( MAPL, nx  , 'NX:'     , default=0   , RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource( MAPL, ny  , 'NY:'     , default=0   , RC=STATUS )
       VERIFY_(STATUS)

       FV_Atm(1)%layout(1) = nx
       FV_Atm(1)%layout(2) = nx
       FV_Atm(1)%flagstruct%npx = npxa + 1
       FV_Atm(1)%flagstruct%npy = npxa + 1 
       FV_Atm(1)%flagstruct%npz = npza
       FV_Atm(1)%flagstruct%ntiles = 6

       !Some assertions incase of user mistake
       ASSERT_( 6*FV_Atm(1)%layout(2) == ny )
       ASSERT_( 6*(FV_Atm(1)%flagstruct%npy-1) == npya )
       ASSERT_(    FV_Atm(1)%flagstruct%ntiles == 6 )
       ASSERT_( DT > 0.0_FVPRC)
       ASSERT_( FV_Atm(1)%ks <= FV_Atm(1)%flagstruct%npz+1 )

       !Constants
       pi     = real(MAPL_PI_R8,FVPRC)
       omega  = real(MAPL_OMEGA,FVPRC)
       cp     = real(MAPL_CP,FVPRC)
       radius = real(MAPL_RADIUS,FVPRC)
       rgas   = real(MAPL_RGAS,FVPRC)
       rvap   = real(MAPL_RVAP,FVPRC)
       kappa  = real(MAPL_KAPPA,FVPRC)
       grav   = real(MAPL_GRAV,FVPRC)
       hlv    = real(MAPL_ALHL,FVPRC)
       zvir   = real(MAPL_RVAP/MAPL_RGAS - 1.,FVPRC)
       p00    = real(MAPL_P00,FVPRC)

       call MAPL_TimerOn(MAPL,"--FV_INIT")
       call fv_init2(FV_Atm, DT, grids_on_this_pe, p_split)
       call MAPL_TimerOff(MAPL,"--FV_INIT")
       call MAPL_MemUtilsWrite(VM, 'FV_StateMod: FV_INIT', RC=STATUS )
       VERIFY_(STATUS)

       !Write some infomation about the setup
       call WRITE_PARALLEL("Dynamics PE Layout ")
       call WRITE_PARALLEL(FV_Atm(1)%layout(1)    ,format='("NPES_X  : ",(   I3))')
       call WRITE_PARALLEL(FV_Atm(1)%layout(2)    ,format='("NPES_Y  : ",(   I3))')
       call WRITE_PARALLEL((/FV_Atm(1)%flagstruct%npx,FV_Atm(1)%flagstruct%npy,FV_Atm(1)%flagstruct%npz/), &
                           format='("Resolution of dynamics restart     =",3I5)'  )
       call WRITE_PARALLEL(FV_Atm(1)%ks, format='("Number of true pressure levels =", I5)'   )

       FV_Atm(1)%flagstruct%fv_debug = .false.
       prt_minmax = FV_Atm(1)%flagstruct%fv_debug

       call MAPL_MemUtilsWrite(VM, trim(Iam), RC=STATUS )
       VERIFY_(STATUS)

       !Set up the tracers
       !------------------
       call MAPL_GetResource( MAPL, FV_Atm(1)%flagstruct%ncnst  , 'ncnst:' , default=4,   rc = STATUS )
       VERIFY_(STATUS)

       !If doing moist physics add 1 tracer for the anvil fraction.
       if (DO_MOIST_PHYS .ne. 0) then
          FV_Atm(1)%flagstruct%ncnst = FV_Atm(1)%flagstruct%ncnst + 1
       endif

       !Reallocate based on GEOS tracers, rather than field_table
       if (allocated(FV_Atm(1)%q)) deallocate( FV_Atm(1)%q  )
       allocate  ( FV_Atm(1)%q (FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,FV_Atm(1)%flagstruct%npz, FV_Atm(1)%flagstruct%ncnst) )

       !Construct perturbation variables
       !--------------------------------
       call fv_init_pert(FV_Atm, FV_AtmP)


       !Reallocate non-hydrostatic variables
       deallocate( FV_Atm(1)%w )
       deallocate( FV_Atm(1)%delz )
       deallocate( FV_AtmP(1)%wp )
       deallocate( FV_AtmP(1)%delzp )

       allocate  ( FV_Atm(1)%w (FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,FV_Atm(1)%flagstruct%npz) )
       allocate  ( FV_Atm(1)%delz (FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,FV_Atm(1)%flagstruct%npz) )
       allocate  ( FV_AtmP(1)%wp (FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,FV_Atm(1)%flagstruct%npz) )
       allocate  ( FV_AtmP(1)%delzp (FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,FV_Atm(1)%flagstruct%npz) )

       !Set up the checkpointing options for this module
       !------------------------------------------------

       if (cp_iter_controls%cp_i .ne. 0) then

          !Name for the module, 3 character max
          call MAPL_GetResource( MAPL, cp_dyn_ind, 'CP_dyn_ind:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)

          !Name for this module in checkpointing
          cp_iter(cp_dyn_ind)%my_name(1:3) = 'dyn'

          !Run in test mode
          cp_iter(cp_dyn_ind)%cp_test = .false.
          call MAPL_GetResource( MAPL, tmp       , 'CP_dyn_test:', default = 0, RC = STATUS )
          VERIFY_(STATUS)
          if (tmp==1) cp_iter(cp_dyn_ind)%cp_test = .true.
   
          !Write reports on checkpointing
          cp_iter(cp_dyn_ind)%cp_rep = .false.
          call MAPL_GetResource( MAPL, tmp       , 'CP_dyn_rep:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          if (tmp==1) cp_iter(cp_dyn_ind)%cp_rep = .true.
  
          !Check status of checkpoints?
          cp_iter(cp_dyn_ind)%check_st_control = .false.
          cp_iter(cp_dyn_ind)%check_st_integer = .false.
          cp_iter(cp_dyn_ind)%check_st_real_r4 = .false.
          cp_iter(cp_dyn_ind)%check_st_real_r8 = .false.
          call MAPL_GetResource( MAPL, tmp, 'CP_dyn_chk_st_control:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          if (tmp==1) cp_iter(cp_dyn_ind)%check_st_control = .true.
          call MAPL_GetResource( MAPL, tmp, 'CP_dyn_chk_st_integer:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          if (tmp==1) cp_iter(cp_dyn_ind)%check_st_integer = .true.
          call MAPL_GetResource( MAPL, tmp, 'CP_dyn_chk_st_real_r4:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          if (tmp==1) cp_iter(cp_dyn_ind)%check_st_real_r4 = .true.
          call MAPL_GetResource( MAPL, tmp, 'CP_dyn_chk_st_real_r8:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          if (tmp==1) cp_iter(cp_dyn_ind)%check_st_real_r8 = .true.
    
          !Number of checkpoints in test mode
          call MAPL_GetResource( MAPL, cp_iter(cp_dyn_ind)%test_dim_st_control, 'CP_dyn_dim_st_control:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          call MAPL_GetResource( MAPL, cp_iter(cp_dyn_ind)%test_dim_st_integer, 'CP_dyn_dim_st_integer:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          call MAPL_GetResource( MAPL, cp_iter(cp_dyn_ind)%test_dim_st_real_r4, 'CP_dyn_dim_st_real_r4:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          call MAPL_GetResource( MAPL, cp_iter(cp_dyn_ind)%test_dim_st_real_r8, 'CP_dyn_dim_st_real_r8:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
   
          !Size of arrays to hold checkpoints in test mode
          call MAPL_GetResource( MAPL, cp_iter(cp_dyn_ind)%test_dim_cp_control, 'CP_dyn_dim_cp_control:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          call MAPL_GetResource( MAPL, cp_iter(cp_dyn_ind)%test_dim_cp_integer, 'CP_dyn_dim_cp_integer:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          call MAPL_GetResource( MAPL, cp_iter(cp_dyn_ind)%test_dim_cp_real_r4, 'CP_dyn_dim_cp_real_r4:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)
          call MAPL_GetResource( MAPL, cp_iter(cp_dyn_ind)%test_dim_cp_real_r8, 'CP_dyn_dim_cp_real_r8:' , default = 0, RC = STATUS )
          VERIFY_(STATUS)

       endif

       call MAPL_GetResource( MAPL, CustomNSplit, 'CustomNSplit:', default = 0, RC = STATUS ); VERIFY_(STATUS)
    
       call MAPL_GetResource( MAPL, tmp, 'fv_timers:', default = 0, RC = STATUS ); VERIFY_(STATUS)
       fv_timing_onoff = .false.
       if (tmp==1) fv_timing_onoff = .true.

       FV_Atm(1)%idiag%id_ws = 0
       FV_Atm(1)%idiag%id_te = 0
       FV_Atm(1)%idiag%id_amdt = 0
       FV_Atm(1)%idiag%id_divg = 0
       FV_Atm(1)%idiag%id_aam = 0
       FV_Atm(1)%idiag%zxg = 0
       FV_Atm(1)%idiag%id_mdt = 0

    endif

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="INITIALIZE"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--FMS_INIT"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--FV_INIT"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUNTLM"      ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--TLM_CORE"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUNADJ"      ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--ADJ_CORE"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="FINALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE,   Initialize, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,          Run,        rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,          Run,        rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_FINALIZE, Finalize,   rc=status)
    VERIFY_(STATUS)
 
! Generic SetServices
!--------------------

! Register prototype of cubed sphere grid and associated regridders
!------------------------------------------------------------------
    call register_grid_and_regridders()


    call MAPL_GenericSetServices( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine Initialize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: gc       ! composite gridded component 
  type(ESMF_State),    intent(inout) :: import   ! import state
  type(ESMF_State),    intent(inout) :: export   ! export state
  type(ESMF_Clock),    intent(inout) :: clock    ! the clock
  integer, optional,   intent(  out) :: rc       ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error
! Locals
!-------
  integer                            :: STATUS
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME
  type (MAPL_MetaComp),      pointer :: MAPL 
  type (ESMF_State)                  :: INTERNAL

! Locals
  integer                            :: i,j
  real(FVPRC)                        :: f_coriolis_angle
  

! Begin
! -----

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the MAPL object
! ---------------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Call Generic Initialize
! -----------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    !FV_InitState
    !------------

    f_coriolis_angle = 0.0

    if (FV_Atm(1)%flagstruct%grid_type == 4) then
       FV_Atm(1)%gridstruct%fC(:,:) = 2.*omega*sin(FV_Atm(1)%flagstruct%deglat/180.*pi)
       FV_Atm(1)%gridstruct%f0(:,:) = 2.*omega*sin(FV_Atm(1)%flagstruct%deglat/180.*pi)
    else
       if (f_coriolis_angle == -999) then
          FV_Atm(1)%gridstruct%fC(:,:) = 0.0
          FV_Atm(1)%gridstruct%f0(:,:) = 0.0
       else
          do j=FV_Atm(1)%bd%jsd,FV_Atm(1)%bd%jed+1
             do i=FV_Atm(1)%bd%isd,FV_Atm(1)%bd%ied+1
                FV_Atm(1)%gridstruct%fC(i,j) = 2.*omega*( -COS(FV_Atm(1)%gridstruct%grid(i,j,1))*&
                                               COS(FV_Atm(1)%gridstruct%grid(i,j,2))*SIN(f_coriolis_angle) + &
                                               SIN(FV_Atm(1)%gridstruct%grid(i,j,2))*COS(f_coriolis_angle) )
             enddo
          enddo
          do j=FV_Atm(1)%bd%jsd,FV_Atm(1)%bd%jed
             do i=FV_Atm(1)%bd%isd,FV_Atm(1)%bd%ied
                FV_Atm(1)%gridstruct%f0(i,j) = 2.*omega*( -COS(FV_Atm(1)%gridstruct%agrid(i,j,1))*&
                                               COS(FV_Atm(1)%gridstruct%agrid(i,j,2))*SIN(f_coriolis_angle) + &
                                               SIN(FV_Atm(1)%gridstruct%agrid(i,j,2))*COS(f_coriolis_angle) )
             enddo
          enddo
       endif
    endif

! All done
! --------

    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_TimerOff(MAPL,"INITIALIZE")

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize


!---------------------------------------------------------------------

!BOP

! !IROUTINE: Run

! !DESCRIPTION: This is the first Run stage of FV.
!
! !INTERFACE:

subroutine Run(gc, import, export, clock, rc)

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: gc
  type (ESMF_State),   intent(inout) :: import
  type (ESMF_State),   intent(inout) :: export
  type (ESMF_Clock),   intent(inout) :: clock
  integer,  optional,  intent(  out) :: rc 

!EOP

! !Local Variables:
  
    type (MAPL_MetaComp),      pointer :: MAPL 
  
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: IAm
    character(len=ESMF_MAXSTR)         :: COMP_NAME
    
    type (ESMF_FieldBundle)            :: Traj3, Traj2, Pert     
    type (T_3DVar), allocatable        :: Trajwrk3(:)
    type (T_3DVar), allocatable        :: Trajwrk2(:)
    type (T_3DVar), allocatable        :: Pertwrk(:)

    !Local convenience
    integer  :: isc,iec,jsc,jec
    integer  :: isd,ied,jsd,jed
    integer  :: npz, ncnst

    integer  :: i,j,k,l

    !Pointers to phis and ps
    real(4), pointer :: TRJ_PHIS(:,:), TRJ_PS(:,:)

    !For computing delp from Ps
    real(FVPRC), allocatable :: PrT_PS(:,:,:)

    !For compting the cloud ratio
    real(FVPRC), allocatable, dimension(:,:,:) :: Ph, PIh, TEMP, fQi
    real(FVPRC), allocatable, dimension(:,:)   :: ebuffery, nbufferx, wbuffery, sbufferx

    character(len=ESMF_MAXSTR)         :: PHASE_NAME
    integer                            :: PHASE
    type (ESMF_State)                  :: INTERNAL
    integer  :: idmodel
    logical, save:: first=.true.

    type (ESMF_Time) :: CURRENTTIME
    integer :: YY, MM, DD, HH, MN, nn

    real(FVPRC) :: tempp
    real(FVPRC), allocatable, dimension(:,:,:) :: u_dt, v_dt, t_dt

! Begin
!------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, idmodel, 'AGCM_PERT_IDMODEL:', default=0, RC = STATUS )
    VERIFY_(STATUS)

! Option to skip the dynamics
! ---------------------------

    IF (DO_DYNAMICS == 0) then ! nothing to do
        if(first) call write_parallel("Skip dynamics")
        first=.false.
        RETURN_(ESMF_SUCCESS)
        return
    endif

! Start the timer
! ---------------
    call MAPL_TimerOn(MAPL,"TOTAL")

! Model time & path for internal FV trajectories
! ----------------------------------------------

    call ESMF_ClockGet(CLOCK, currTIME=CURRENTTIME,       RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet (CurrentTime, YY=YY, MM=MM, DD=DD, H=HH, M=MN, RC=STATUS )
    VERIFY_(STATUS)

! Get the TLM/ADM phase
! ---------------------

    call ESMF_GridCompGet(GC, currentPhase=Phase, rc=STATUS)
    VERIFY_(STATUS)

    select case(phase)
    case(TLMphase)
       PHASE_NAME = "RUNTLM"
    case(ADJphase)
       PHASE_NAME = "RUNADJ"
    case default
       ASSERT_(.false.)
    end select

! Start phase timer
! -----------------

    if (phase==TLMPhase) then
       call timing_on('DYNCORE_TLM')
    elseif (phase==ADJPhase) then
       call timing_on('DYNCORE_ADM')
    endif

    call MAPL_TimerOn(MAPL,PHASE_NAME)

! Get working trajectory and pertubation from IMPORT  
!---------------------------------------------------

    call ESMF_StateGet(Import, 'TRAJWRK3', Traj3, rc=STATUS); VERIFY_(STATUS)
    call ESMF_StateGet(Import, 'TRAJWRK2', Traj2, rc=STATUS); VERIFY_(STATUS)
    call ESMF_StateGet(Import, 'PERTWRK',  Pert,  rc=STATUS); VERIFY_(STATUS)

! Incedes of local and ghosted domains
!-------------------------------------

    isc = FV_Atm(1)%bd%isc
    iec = FV_Atm(1)%bd%iec
    jsc = FV_Atm(1)%bd%jsc
    jec = FV_Atm(1)%bd%jec
    isd = FV_Atm(1)%bd%isd
    ied = FV_Atm(1)%bd%ied
    jsd = FV_Atm(1)%bd%jsd
    jed = FV_Atm(1)%bd%jed
    npz = FV_Atm(1)%npz
    ncnst = FV_Atm(1)%ncnst

    if (CustomNSplit > 0) then
       FV_Atm(1)%flagstruct%n_split = CustomNSplit
       if ( is_master() ) write(*,"(A,I3.1)") 'WARNING, running with custom n-splits: ',CustomNSplit
    endif
 
! Allocate FV vars for trajectory
!--------------------------------
    allocate(Trajwrk3(NUM_GCM3DVars))
    allocate(Trajwrk2(NUM_GCM2Dvars))
    allocate(Pertwrk(NUM_GCM3DVars))

    !Variables for pressure calculations
    allocate(PrT_PS(isd  :ied  ,jsd:jed  ,npz+1 ), stat=status); VERIFY_(STATUS)
    allocate(Ph    (isd  :ied  ,jsd:jed  ,npz   ), stat=status); VERIFY_(STATUS)
    allocate(PIh   (isd  :ied  ,jsd:jed  ,npz   ), stat=status); VERIFY_(STATUS)
    allocate(TEMP  (isd  :ied  ,jsd:jed  ,npz   ), stat=status); VERIFY_(STATUS)
    allocate(fQi   (isd  :ied  ,jsd:jed  ,npz   ), stat=status); VERIFY_(STATUS)

    allocate(ebuffery(jsd:jed,npz))
    allocate(wbuffery(jsd:jed,npz))
    allocate(nbufferx(isd:ied,npz))
    allocate(sbufferx(isd:ied,npz))
    ebuffery = 0.0_FVPRC
    wbuffery = 0.0_FVPRC
    nbufferx = 0.0_FVPRC
    sbufferx = 0.0_FVPRC

! Copy ESMF bundles into simple bundles of the working copies of Traj and Pert
!-----------------------------------------------------------------------------

    do i=1, NUM_GCM3Dvars
       TrajWrk3(i)%name = GCM3Dvars(i)
       PertWrk(i)%name  = GCM3Dvars(i)
    enddo
    TrajWrk2(1)%name = 'PHIS'
    TrajWrk2(2)%name = 'PS'

    call PERT_RefVarsFromBundle(Traj3, TrajWrk3, rc=STATUS); VERIFY_(STATUS)
    call PERT_RefVarsFromBundle(Pert,  PertWrk,  rc=STATUS); VERIFY_(STATUS)

    call MAPL_FieldBundleGetPointer(Traj2, "PHIS", TRJ_PHIS, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_FieldBundleGetPointer(Traj2, "PS", TRJ_PS, rc=STATUS)
    VERIFY_(STATUS)

    if (DynTestCase == 0) then
       FV_Atm(1)%phis = 0.0_FVPRC
       FV_Atm(1)%ps   = 0.0_FVPRC
       FV_Atm(1)%phis(isc:iec,jsc:jec) = TRJ_PHIS
       FV_Atm(1)%ps  (isc:iec,jsc:jec) = TRJ_Ps
       call mpp_update_domains(FV_Atm(1)%phis, FV_Atm(1)%domain)
       call mpp_update_domains(FV_Atm(1)%ps  , FV_Atm(1)%domain)
    endif

! Caclulate DPT from STATE_PS
!----------------------------

    FV_Atm(1)%ak(1:npz+1) = pert_ak(0:npz)
    FV_Atm(1)%bk(1:npz+1) = pert_bk(0:npz)
    FV_Atm(1)%ptop = FV_Atm(1)%ak(1)

    !Pressure at the half levels from FV_Atm(1)%ps
    PrT_PS = 0.0_FVPRC
    DO L = 1,npz+1
       PrT_PS(isc:iec,jsc:jec,L) = FV_ATM(1)%ak(L) + FV_ATM(1)%bk(L)*FV_Atm(1)%ps(isc:iec,jsc:jec)
    ENDDO 

! Set up the reference trajectory variables
!------------------------------------------

    if (DynTestCase == 0) then

       !Make sure halo is zero
       FV_Atm(1)%u  = 0.0_FVPRC
       FV_Atm(1)%v  = 0.0_FVPRC
       FV_Atm(1)%pt = 0.0_FVPRC
       FV_Atm(1)%delp = 0.0_FVPRC
       FV_Atm(1)%q  = 0.0_FVPRC

       !Trajectory
       FV_Atm(1)%u   (isc:iec,jsc:jec,:)   = TrajWrk3(PERT_GetIndex(TrajWrk3,"U" ))%X
       FV_Atm(1)%v   (isc:iec,jsc:jec,:)   = TrajWrk3(PERT_GetIndex(TrajWrk3,"V" ))%X
       FV_Atm(1)%delp(isc:iec,jsc:jec,:)   = PrT_PS(isc:iec,jsc:jec,2:npz+1) - PrT_PS(isc:iec,jsc:jec,1:npz) !From State_PS
       FV_Atm(1)%pt  (isc:iec,jsc:jec,:)   = TrajWrk3(PERT_GetIndex(TrajWrk3,"PT"))%X
       FV_Atm(1)%q   (isc:iec,jsc:jec,:,1) = TrajWrk3(PERT_GetIndex(TrajWrk3,"QV"))%X
       FV_Atm(1)%q   (isc:iec,jsc:jec,:,2) = TrajWrk3(PERT_GetIndex(TrajWrk3,"O3"))%X
       if (DO_MOIST_PHYS == 0) then
          FV_Atm(1)%q   (isc:iec,jsc:jec,:,3) = TrajWrk3(PERT_GetIndex(TrajWrk3,"QL"))%X
          FV_Atm(1)%q   (isc:iec,jsc:jec,:,4) = TrajWrk3(PERT_GetIndex(TrajWrk3,"QI"))%X
       else
          Ph   = 0.0_FVPRC
          PIh  = 0.0_FVPRC
          TEMP = 0.0_FVPRC
          fQi  = 0.0_FVPRC
          Ph(isc:iec,jsc:jec,:) = 0.01_FVPRC * 0.5_FVPRC * (PrT_PS(isc:iec,jsc:jec,1:npz) +  PrT_PS(isc:iec,jsc:jec,2:npz+1  ) ) 
          PIh(isc:iec,jsc:jec,:) = (Ph(isc:iec,jsc:jec,:)/1000.0_FVPRC)**(rgas/cp)
          TEMP(isc:iec,jsc:jec,:) = FV_Atm(1)%pt(isc:iec,jsc:jec,:)*(p00**kappa)*PIh(isc:iec,jsc:jec,:)
          DO i = isc,iec
             DO j = jsc,jec
                DO l = 1,npz
                   call IceFraction( TEMP(i,j,l), fQi(i,j,l) )
                enddo
             enddo
          enddo

          !Read in QLLS & QLCN/QILS & QICN and combine for QL/QI
          FV_Atm(1)%q (isc:iec,jsc:jec,:,3) = (TrajWrk3(PERT_GetIndex(TrajWrk3,"QLS"))%X + &
                                      TrajWrk3(PERT_GetIndex(TrajWrk3,"QCN"))%X) * (1 - fQi(isc:iec,jsc:jec,:))
          FV_Atm(1)%q (isc:iec,jsc:jec,:,4) = (TrajWrk3(PERT_GetIndex(TrajWrk3,"QLS"))%X + &
                                      TrajWrk3(PERT_GetIndex(TrajWrk3,"QCN"))%X) *      fQi(isc:iec,jsc:jec,:)

          !Cloud fractions
          FV_Atm(1)%q (isc:iec,jsc:jec,:,5) = TrajWrk3(PERT_GetIndex(TrajWrk3,"CFCN"))%X
       endif

       !These are computed internally in the core, are outputs or are superfluous variables
       FV_Atm(1)%pe   = 0.0_FVPRC
       FV_Atm(1)%peln = 0.0_FVPRC
       FV_Atm(1)%pk   = 0.0_FVPRC
       FV_Atm(1)%pkz  = 0.0_FVPRC
       FV_Atm(1)%ua   = 0.0_FVPRC
       FV_Atm(1)%va   = 0.0_FVPRC
       FV_Atm(1)%uc   = 0.0_FVPRC
       FV_Atm(1)%vc   = 0.0_FVPRC
       FV_Atm(1)%omga = 0.0_FVPRC
       FV_Atm(1)%mfx  = 0.0_FVPRC
       FV_Atm(1)%mfy  = 0.0_FVPRC
       FV_Atm(1)%cx   = 0.0_FVPRC
       FV_Atm(1)%cy   = 0.0_FVPRC
       FV_Atm(1)%ze0  = 0.0_FVPRC
       FV_Atm(1)%w    = 0.0_FVPRC
       FV_Atm(1)%delz = 0.0_FVPRC
       FV_Atm(1)%q_con = 0.0_FVPRC

       !Remapping things 
       FV_Atm(1)%flagstruct%reproduce_sum = .false.

       !Never fill tracers
       FV_Atm(1)%flagstruct%fill = .false.

       !Never debug
       FV_Atm(1)%flagstruct%fv_debug = .false.

       !Never adiabatic
       FV_Atm(1)%flagstruct%adiabatic = .false.

       !Never do saturation adjustment
       FV_Atm(1)%flagstruct%do_sat_adj = .false.

       !Never do this
       FV_Atm(1)%flagstruct%breed_vortex_inline = .false.

    endif

! Set up the perturbation variables
!----------------------------------

    if (DynTestCase == 0) then
   
       !Make sure halo is zero
       FV_AtmP(1)%up  = 0.0_FVPRC
       FV_AtmP(1)%vp  = 0.0_FVPRC
       FV_AtmP(1)%ptp = 0.0_FVPRC
       FV_AtmP(1)%delpp = 0.0_FVPRC
       FV_AtmP(1)%qp  = 0.0_FVPRC

       !Perturbation trajectory                                                            
       FV_AtmP(1)%up   (isc:iec,jsc:jec,:)   = PertWrk(PERT_GetIndex(PertWrk,"U" ))%X
       FV_AtmP(1)%vp   (isc:iec,jsc:jec,:)   = PertWrk(PERT_GetIndex(PertWrk,"V" ))%X
       FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)   = PertWrk(PERT_GetIndex(PertWrk,"PT"))%X
       FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)   = PertWrk(PERT_GetIndex(PertWrk,"DP"))%X
       FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,1) = PertWrk(PERT_GetIndex(PertWrk,"QV"))%X
       FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,2) = PertWrk(PERT_GetIndex(PertWrk,"O3"))%X
       FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,3) = PertWrk(PERT_GetIndex(PertWrk,"QL"))%X
       FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,4) = PertWrk(PERT_GetIndex(PertWrk,"QI"))%X
       if   (DO_MOIST_PHYS .ne. 0) then
          FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,5) = PertWrk(PERT_GetIndex(PertWrk,"CFCN"))%X
       endif
   
       !Zero out FV3 perturbation variables
       FV_AtmP(1)%pep   = 0.0_FVPRC
       FV_AtmP(1)%pelnp = 0.0_FVPRC
       FV_AtmP(1)%pkp   = 0.0_FVPRC
       FV_AtmP(1)%pkzp  = 0.0_FVPRC
       FV_AtmP(1)%uap   = 0.0_FVPRC
       FV_AtmP(1)%vap   = 0.0_FVPRC
       FV_AtmP(1)%ucp   = 0.0_FVPRC
       FV_AtmP(1)%vcp   = 0.0_FVPRC
       FV_AtmP(1)%omgap = 0.0_FVPRC
       FV_AtmP(1)%mfxp  = 0.0_FVPRC
       FV_AtmP(1)%mfyp  = 0.0_FVPRC
       FV_AtmP(1)%cxp   = 0.0_FVPRC
       FV_AtmP(1)%cyp   = 0.0_FVPRC
       FV_AtmP(1)%wp    = 0.0_FVPRC
       FV_AtmP(1)%delzP = 0.0_FVPRC
       FV_AtmP(1)%q_conP = 0.0_FVPRC

    endif

    call set_domain(FV_Atm(1)%domain)
 
    if(phase==TLMPhase) then

       if (DynTestCase > 0) then

          !Run the chosen test case
          call DynCoreTestCases
 
       else
 
          call mpp_get_boundary( FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%domain, &
                                 wbuffery=wbuffery, ebuffery=ebuffery, &
                                 sbufferx=sbufferx, nbufferx=nbufferx, &
                                 gridtype=DGRID_NE, complete=.true. )
          do k=1,npz
             do i=isc,iec
                FV_Atm(1)%u(i,jec+1,k) = nbufferx(i,k)
             enddo
          enddo
          do k=1,npz
             do j=jsc,jec
                FV_Atm(1)%v(iec+1,j,k) = ebuffery(j,k)
             enddo
          enddo 
   
          call mpp_get_boundary( FV_AtmP(1)%up, FV_AtmP(1)%vp, FV_Atm(1)%domain, &
                                 wbuffery=wbuffery, ebuffery=ebuffery, &
                                 sbufferx=sbufferx, nbufferx=nbufferx, &
                                 gridtype=DGRID_NE, complete=.true. )
          do k=1,npz
             do i=isc,iec
                FV_AtmP(1)%up(i,jec+1,k) = nbufferx(i,k)
             enddo
          enddo
          do k=1,npz
             do j=jsc,jec
                FV_AtmP(1)%vp(iec+1,j,k) = ebuffery(j,k)
             enddo
          enddo

          call timing_on(' FV_DYNAMICS_TLM')

          call compute_pressures_tlm( FV_Atm(1)%bd, FV_Atm(1)%npz, kappa, FV_Atm(1)%ptop, &
                                      FV_Atm(1)%delp, FV_AtmP(1)%delpp, FV_Atm(1)%pe, FV_AtmP(1)%pep, FV_Atm(1)%pk, FV_AtmP(1)%pkp, &
                                      FV_Atm(1)%pkz, FV_AtmP(1)%pkzp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp )

          !Convert potential temperature to dry temperature (GEOS to FV)
          do k=1,FV_Atm(1)%npz
            do j=jsc,jec
              do i=isc,iec
                FV_AtmP(1)%ptp(i,j,k) = FV_AtmP(1)%ptp(i,j,k)*FV_Atm(1)%pkz(i,j,k) + FV_Atm(1)%pt(i,j,k)*FV_AtmP(1)%pkzp(i,j,k)
                FV_Atm(1)%pt(i,j,k) = FV_Atm(1)%pt(i,j,k)*FV_Atm(1)%pkz(i,j,k)
              end do
            end do
          end do

          call fv_dynamics_tlm(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
                               DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                                    &
                               FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
                               cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
                               FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
                               FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,              &
                               FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                              &
                               FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                  &
                               FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,          &
                               FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,  &
                               FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga, FV_AtmP(1)%omgap,                               &
                               FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                      &
                               FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                      &
                               FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
                               FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                  &
                               FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp,                                     &
                               FV_Atm(1)%ze0,                                                                                   &
                               FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct,                                             &
                               FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct, FV_Atm(1)%neststruct,                               &
                               FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain                            )

          if (FV_Atm(1)%flagstruct%fv_sg_adj .gt. 0) then

             allocate(u_dt(isd:ied,jsd:jed,npz)); u_dt = 0.0
             allocate(v_dt(isd:ied,jsd:jed,npz)); v_dt = 0.0
             allocate(t_dt(isc:iec,jsc:jec,npz)); t_dt = 0.0

             call fv_subgrid_z_tlm( isd, ied, jsd, jed, isc, iec, jsc, jec, npz, ncnst, &
                                    dt, FV_Atm(1)%flagstruct%fv_sg_adj, FV_Atm(1)%flagstruct%nwat, &
                                    FV_Atm(1)%delp, FV_AtmP(1)%delpp, FV_Atm(1)%pe, FV_AtmP(1)%pep, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, &
                                    FV_Atm(1)%pkz, FV_AtmP(1)%pkzp, FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%q, FV_AtmP(1)%qp, &
                                    FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap, &
                                    FV_Atm(1)%flagstruct%hydrostatic, &
                                    FV_Atm(1)%w, FV_AtmP(1)%wp, FV_Atm(1)%delz, FV_AtmP(1)%delzp, &
                                    u_dt, v_dt, t_dt, FV_Atm(1)%flagstruct%n_zfilter )

             deallocate(u_dt, v_dt, t_dt)

          end if

          !Convert back to potential temperature
          do k=1,FV_Atm(1)%npz
            do j=jsc,jec
              do i=isc,iec
                FV_AtmP(1)%ptp(i,j,k) = ( FV_AtmP(1)%ptp(i,j,k)*FV_Atm(1)%pkz(i,j,k) - &
                                          FV_Atm(1)%pt(i,j,k)*FV_AtmP(1)%pkzp(i,j,k) )/FV_Atm(1)%pkz(i,j,k)**2
                FV_Atm(1)%pt(i,j,k) = FV_Atm(1)%pt(i,j,k)/FV_Atm(1)%pkz(i,j,k)
              end do
            end do
          end do

          call timing_off(' FV_DYNAMICS_TLM')
   
       endif

    else
 
       !Initilize the this iterative step
       if (cp_iter_controls%cp_i .ne. 0) then
          call cp_mod_ini(cp_dyn_ind)
       endif

       call mpp_get_boundary( FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%domain, &
                              wbuffery=wbuffery, ebuffery=ebuffery, sbufferx=sbufferx, nbufferx=nbufferx, &
                              gridtype=DGRID_NE, complete=.true. )
       do k=1,npz
          do i=isc,iec
             FV_Atm(1)%u(i,jec+1,k) = nbufferx(i,k)
          enddo
       enddo
       do k=1,npz
          do j=jsc,jec
             FV_Atm(1)%v(iec+1,j,k) = ebuffery(j,k)
          enddo
       enddo

       call timing_on(' FV_DYNAMICS_FWD')

       if (cp_iter_controls%cp_i <= 3) then

          call compute_pressures_fwd( FV_Atm(1)%bd, FV_Atm(1)%npz, kappa, FV_Atm(1)%ptop, &
                                      FV_Atm(1)%delp, FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%pkz, FV_Atm(1)%peln )

          !Convert potential temperature to dry temperature
          call pushrealarray(FV_Atm(1)%pt(isc:iec,jsc:jec,:),(iec-isc+1)*(jec-jsc+1)*npz)
          do k=1,FV_Atm(1)%npz
            do j=jsc,jec
              do i=isc,iec
                FV_Atm(1)%pt(i, j, k) = FV_Atm(1)%pt(i, j, k)*FV_Atm(1)%pkz(i, j, k)
              end do
            end do
          end do

          !Forward sweep of the dynamics with saving of checkpoints for use in backward sweep.
          call fv_dynamics_fwd(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
                               DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                                    &
                               FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
                               cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
                               FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
                               FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w,                                                           &
                               FV_Atm(1)%delz, FV_Atm(1)%flagstruct%hydrostatic,                                                &
                               FV_Atm(1)%pt, FV_Atm(1)%delp,                                                                    &
                               FV_Atm(1)%q, FV_Atm(1)%ps, FV_Atm(1)%pe,                                                         &
                               FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,                                                     &
                               FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga,                                                 &
                               FV_Atm(1)%ua, FV_Atm(1)%va,                                                                      &
                               FV_Atm(1)%uc, FV_Atm(1)%vc,                                                                      &
                               FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
                               FV_Atm(1)%mfx, FV_Atm(1)%mfy,                                                                    &
                               FV_Atm(1)%cx, FV_Atm(1)%cy,                                                                      &
                               FV_Atm(1)%ze0,                                                                                   &
                               FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct,                                             &
                               FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct, FV_Atm(1)%neststruct,                               &
                               FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain                            )
 
          if (FV_Atm(1)%flagstruct%fv_sg_adj .gt. 0) then
             allocate(u_dt(isd:ied,jsd:jed,npz)); u_dt = 0.0
             allocate(v_dt(isd:ied,jsd:jed,npz)); v_dt = 0.0
             allocate(t_dt(isc:iec,jsc:jec,npz)); t_dt = 0.0
             call fv_subgrid_z_fwd( isd, ied, jsd, jed, isc, iec, jsc, jec, npz, ncnst, &
                                    dt, FV_Atm(1)%flagstruct%fv_sg_adj, FV_Atm(1)%flagstruct%nwat, &
                                    FV_Atm(1)%delp, FV_Atm(1)%pe, FV_Atm(1)%peln, &
                                    FV_Atm(1)%pkz, FV_Atm(1)%pt, FV_Atm(1)%q, &
                                    FV_Atm(1)%ua, FV_Atm(1)%va, &
                                    FV_Atm(1)%flagstruct%hydrostatic, &
                                    FV_Atm(1)%w, FV_Atm(1)%delz, &
                                    u_dt, v_dt, t_dt, FV_Atm(1)%flagstruct%n_zfilter )
             deallocate(u_dt, v_dt, t_dt)
          end if
 
          call pushrealarray(FV_Atm(1)%cy, (ied-isd+1)*(jec-jsc+2)*npz)
          call pushrealarray(FV_Atm(1)%cx, (iec-isc+2)*(jed-jsd+1)*npz)
          call pushrealarray(FV_Atm(1)%pk, (iec-isc+1)*(jec-jsc+1)*(npz+1))
          call pushrealarray(FV_Atm(1)%pe, (iec-isc+3)*(npz+1)*(jec-jsc+3))
          call pushrealarray(FV_Atm(1)%pkz, (iec-isc+1)*(jec-jsc+1)*npz)
          call pushrealarray(FV_Atm(1)%va, (ied-isd+1)*(jed-jsd+1)*npz)
          call pushrealarray(FV_Atm(1)%omga, (ied-isd+1)*(jed-jsd+1)*npz)
          call pushrealarray(FV_Atm(1)%mfy, (iec-isc+1)*(jec-jsc+2)*npz)
          call pushrealarray(FV_Atm(1)%mfx, (iec-isc+2)*(jec-jsc+1)*npz)
          call pushrealarray(FV_Atm(1)%ua, (ied-isd+1)*(jed-jsd+1)*npz)
          call pushrealarray(FV_Atm(1)%peln, (iec-isc+1)*(npz+1)*(jec-jsc+1))
          call pushrealarray(FV_Atm(1)%phis, (ied-isd+1)*(jed-jsd+1))

          if (cp_iter_controls%cp_i .ne. 0) then
             !Push end of timestep trajectory to stack
             call pushrealarray(FV_Atm(1)%u   ,(ied-isd+1)*(jed-jsd+2)*npz)
             call pushrealarray(FV_Atm(1)%v   ,(ied-isd+2)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%w   ,(ied-isd+1)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%delz,(ied-isd+1)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%pt  ,(ied-isd+1)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%delp,(ied-isd+1)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%q   ,(ied-isd+1)*(jed-jsd+1)*npz*ncnst)
             call pushrealarray(FV_Atm(1)%ps  ,(ied-isd+1)*(jed-jsd+1))
             call pushrealarray(FV_Atm(1)%pe  ,(iec-isc+3)*(jec-jsc+3)*(npz+1))
             call pushrealarray(FV_Atm(1)%pk  ,(iec-isc+1)*(jec-jsc+1)*(npz+1))
             call pushrealarray(FV_Atm(1)%peln,(iec-isc+1)*(jec-jsc+1)*(npz+1))
             call pushrealarray(FV_Atm(1)%pkz ,(iec-isc+1)*(jec-jsc+1)*npz)
             call pushrealarray(FV_Atm(1)%phis,(ied-isd+1)*(jed-jsd+1))
             call pushrealarray(FV_Atm(1)%omga,(ied-isd+1)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%ua  ,(ied-isd+1)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%va  ,(ied-isd+1)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%uc  ,(ied-isd+2)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%vc  ,(ied-isd+1)*(jed-jsd+2)*npz)
             call pushrealarray(FV_Atm(1)%mfx ,(iec-isc+2)*(jec-jsc+1)*npz)
             call pushrealarray(FV_Atm(1)%mfy ,(iec-isc+1)*(jec-jsc+2)*npz)
             call pushrealarray(FV_Atm(1)%cx  ,(iec-isc+2)*(jed-jsd+1)*npz)
             call pushrealarray(FV_Atm(1)%cy  ,(ied-isd+1)*(jec-jsc+2)*npz)
             !Trick checkpoint schemes into not considering these superfluous checkpoints,
             !about to recover with the pop anyway.
             FV_Atm(1)%u    = 2.0_FVPRC*FV_Atm(1)%u
             FV_Atm(1)%v    = 2.0_FVPRC*FV_Atm(1)%v
             FV_Atm(1)%w    = 2.0_FVPRC*FV_Atm(1)%w
             FV_Atm(1)%delz = 2.0_FVPRC*FV_Atm(1)%delz
             FV_Atm(1)%pt   = 2.0_FVPRC*FV_Atm(1)%pt
             FV_Atm(1)%delp = 2.0_FVPRC*FV_Atm(1)%delp
             FV_Atm(1)%q    = 2.0_FVPRC*FV_Atm(1)%q
             FV_Atm(1)%ps   = 2.0_FVPRC*FV_Atm(1)%ps
             FV_Atm(1)%pe   = 2.0_FVPRC*FV_Atm(1)%pe
             FV_Atm(1)%pk   = 2.0_FVPRC*FV_Atm(1)%pk
             FV_Atm(1)%peln = 2.0_FVPRC*FV_Atm(1)%peln
             FV_Atm(1)%pkz  = 2.0_FVPRC*FV_Atm(1)%pkz
             FV_Atm(1)%phis = 2.0_FVPRC*FV_Atm(1)%phis
             FV_Atm(1)%omga = 2.0_FVPRC*FV_Atm(1)%omga
             FV_Atm(1)%ua   = 2.0_FVPRC*FV_Atm(1)%ua
             FV_Atm(1)%va   = 2.0_FVPRC*FV_Atm(1)%va
             FV_Atm(1)%uc   = 2.0_FVPRC*FV_Atm(1)%uc
             FV_Atm(1)%vc   = 2.0_FVPRC*FV_Atm(1)%vc
             FV_Atm(1)%mfx  = 2.0_FVPRC*FV_Atm(1)%mfx
             FV_Atm(1)%mfy  = 2.0_FVPRC*FV_Atm(1)%mfy
             FV_Atm(1)%cx   = 2.0_FVPRC*FV_Atm(1)%cx
             FV_Atm(1)%cy   = 2.0_FVPRC*FV_Atm(1)%cy
          endif

       endif

       call timing_off(' FV_DYNAMICS_FWD')

       !Checkpoint mid point, reset counters etc
       if (cp_iter_controls%cp_i .ne. 0) then 
          call cp_mod_mid
       endif

       call timing_on(' FV_DYNAMICS_BWD')

       if (cp_iter_controls%cp_i .ne. 0) then
          !Populate end of timestep trajectory from stack
          call poprealarray(FV_Atm(1)%cy  ,(ied-isd+1)*(jec-jsc+2)*npz)
          call poprealarray(FV_Atm(1)%cx  ,(iec-isc+2)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%mfy ,(iec-isc+1)*(jec-jsc+2)*npz)
          call poprealarray(FV_Atm(1)%mfx ,(iec-isc+2)*(jec-jsc+1)*npz)
          call poprealarray(FV_Atm(1)%vc  ,(ied-isd+1)*(jed-jsd+2)*npz)
          call poprealarray(FV_Atm(1)%uc  ,(ied-isd+2)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%va  ,(ied-isd+1)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%ua  ,(ied-isd+1)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%omga,(ied-isd+1)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%phis,(ied-isd+1)*(jed-jsd+1))
          call poprealarray(FV_Atm(1)%pkz ,(iec-isc+1)*(jec-jsc+1)*npz)
          call poprealarray(FV_Atm(1)%peln,(iec-isc+1)*(jec-jsc+1)*(npz+1))
          call poprealarray(FV_Atm(1)%pk  ,(iec-isc+1)*(jec-jsc+1)*(npz+1))
          call poprealarray(FV_Atm(1)%pe  ,(iec-isc+3)*(jec-jsc+3)*(npz+1))
          call poprealarray(FV_Atm(1)%ps  ,(ied-isd+1)*(jed-jsd+1))
          call poprealarray(FV_Atm(1)%q   ,(ied-isd+1)*(jed-jsd+1)*npz*ncnst)
          call poprealarray(FV_Atm(1)%delp,(ied-isd+1)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%pt  ,(ied-isd+1)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%delz,(ied-isd+1)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%w   ,(ied-isd+1)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%v   ,(ied-isd+2)*(jed-jsd+1)*npz)
          call poprealarray(FV_Atm(1)%u   ,(ied-isd+1)*(jed-jsd+2)*npz)
       endif

       call poprealarray(FV_Atm(1)%phis, (ied-isd+1)*(jed-jsd+1))
       call poprealarray(FV_Atm(1)%peln, (iec-isc+1)*(npz+1)*(jec-jsc+1))
       call poprealarray(FV_Atm(1)%ua, (ied-isd+1)*(jed-jsd+1)*npz)
       call poprealarray(FV_Atm(1)%mfx, (iec-isc+2)*(jec-jsc+1)*npz)
       call poprealarray(FV_Atm(1)%mfy, (iec-isc+1)*(jec-jsc+2)*npz)
       call poprealarray(FV_Atm(1)%omga, (ied-isd+1)*(jed-jsd+1)*npz)
       call poprealarray(FV_Atm(1)%va, (ied-isd+1)*(jed-jsd+1)*npz)
       call poprealarray(FV_Atm(1)%pkz, (iec-isc+1)*(jec-jsc+1)*npz)
       call poprealarray(FV_Atm(1)%pe, (iec-isc+3)*(npz+1)*(jec-jsc+3))
       call poprealarray(FV_Atm(1)%pk, (iec-isc+1)*(jec-jsc+1)*(npz+1))
       call poprealarray(FV_Atm(1)%cx, (iec-isc+2)*(jed-jsd+1)*npz)
       call poprealarray(FV_Atm(1)%cy, (ied-isd+1)*(jec-jsc+2)*npz)

       !Backward adjoint of temperature to potential temperature
       FV_AtmP(1)%pkzp = 0.0
       DO k=FV_Atm(1)%npz,1,-1
         DO j=jec,jsc,-1
           DO i=iec,isc,-1
             tempp = FV_AtmP(1)%ptp(i, j, k)/FV_Atm(1)%pkz(i, j, k)
             FV_AtmP(1)%pkzp(i, j, k) = FV_AtmP(1)%pkzp(i, j, k) - FV_Atm(1)%pt(i, j, k)*tempp/FV_Atm(1)%pkz(i, j, k)
             FV_AtmP(1)%ptp(i, j, k) = tempp
           END DO
         END DO
       END DO

       !Backward adjoint of subgrid mixing
       if (FV_Atm(1)%flagstruct%fv_sg_adj .gt. 0) then
          allocate(u_dt(isd:ied,jsd:jed,npz)); u_dt = 0.0
          allocate(v_dt(isd:ied,jsd:jed,npz)); v_dt = 0.0
          allocate(t_dt(isc:iec,jsc:jec,npz)); t_dt = 0.0
          call fv_subgrid_z_bwd( isd, ied, jsd, jed, isc, iec, jsc, jec, npz, ncnst, &
                                 dt, FV_Atm(1)%flagstruct%fv_sg_adj, FV_Atm(1)%flagstruct%nwat, &
                                 FV_Atm(1)%delp, FV_AtmP(1)%delpp, FV_Atm(1)%pe, FV_AtmP(1)%pep, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, &
                                 FV_Atm(1)%pkz, FV_AtmP(1)%pkzp, FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%q, FV_AtmP(1)%qp, &
                                 FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap, &
                                 FV_Atm(1)%flagstruct%hydrostatic, &
                                 FV_Atm(1)%w, FV_AtmP(1)%wp, FV_Atm(1)%delz, FV_AtmP(1)%delzp, &
                                 u_dt, v_dt, t_dt, FV_Atm(1)%flagstruct%n_zfilter )
          deallocate(u_dt, v_dt, t_dt)
       endif

       !Backward adjoint sweep of the dynamics.
       call fv_dynamics_bwd(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
                            DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                                    &
                            FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
                            cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
                            FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
                            FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,              &
                            FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                              &
                            FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                  &
                            FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,          &
                            FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,  &
                            FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga, FV_AtmP(1)%omgap,                               &
                            FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                      &
                            FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                      &
                            FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
                            FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                  &
                            FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp,                                     &
                            FV_Atm(1)%ze0,                                                                                   &
                            FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct,                                             &
                            FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct, FV_Atm(1)%neststruct,                               &
                            FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain                            )

       !Backward adjoint of potential temperature to temperature
       call poprealarray(FV_Atm(1)%pt(isc:iec,jsc:jec,:),(iec-isc+1)*(jec-jsc+1)*npz) 
       DO k=FV_Atm(1)%npz,1,-1
         DO j=jec,jsc,-1
           DO i=iec,isc,-1
             FV_AtmP(1)%pkzp(i, j, k) = FV_AtmP(1)%pkzp(i, j, k) + FV_Atm(1)%pt(i, j, k)*FV_AtmP(1)%ptp(i, j, k)
             FV_AtmP(1)%ptp(i, j, k) = FV_Atm(1)%pkz(i, j, k)*FV_AtmP(1)%ptp(i, j, k)
           END DO
         END DO
       END DO

       !Backward adjoint of compute pressures
       call compute_pressures_bwd( FV_Atm(1)%bd, FV_Atm(1)%npz, kappa, FV_Atm(1)%ptop, &
                                   FV_Atm(1)%delp, FV_AtmP(1)%delpp, FV_Atm(1)%pe, FV_AtmP(1)%pep, FV_Atm(1)%pk, FV_AtmP(1)%pkp, &
                                   FV_Atm(1)%pkz, FV_AtmP(1)%pkzp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp )

       call timing_off(' FV_DYNAMICS_BWD')

       !Finalize iterative step and get ready for next iteration
       if (cp_iter_controls%cp_i .ne. 0) then
          call cp_mod_end
       endif


       nbufferx = 0.0_FVPRC
       do k=1,npz
          do i=isc,iec
             nbufferx(i,k) = FV_AtmP(1)%up(i,jec+1,k)
          enddo
       enddo
       ebuffery = 0.0_FVPRC
       do k=1,npz
          do j=jsc,jec
             ebuffery(j,k) = FV_AtmP(1)%vp(iec+1,j,k)
          enddo
       enddo

       call mpp_get_boundary_ad( FV_AtmP(1)%up, FV_AtmP(1)%vp, FV_Atm(1)%domain, &
                                 wbuffery=wbuffery, ebuffery=ebuffery, sbufferx=sbufferx, nbufferx=nbufferx, &
                                 gridtype=DGRID_NE, complete=.true. )

    endif
 
    call nullify_domain()

! Copy out of r8/FVPRC, ghosted into the working bundles
!-------------------------------------------------------
    if (DynTestCase == 0) then
       PertWrk(PERT_GetIndex(PertWrk,"U" ))%X = FV_AtmP(1)%up   (isc:iec,jsc:jec,:)
       PertWrk(PERT_GetIndex(PertWrk,"V" ))%X = FV_AtmP(1)%vp   (isc:iec,jsc:jec,:)
       PertWrk(PERT_GetIndex(PertWrk,"PT"))%X = FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)
       PertWrk(PERT_GetIndex(PertWrk,"DP"))%X = FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)
       PertWrk(PERT_GetIndex(PertWrk,"QV"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,1)
       PertWrk(PERT_GetIndex(PertWrk,"O3"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,2)
       PertWrk(PERT_GetIndex(PertWrk,"QL"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,3)
       PertWrk(PERT_GetIndex(PertWrk,"QI"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,4)
       if (DO_MOIST_PHYS .ne. 0) then
          PertWrk(PERT_GetIndex(PertWrk,"CFCN"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,5)
          ncnst = ncnst - 1
       endif
    endif

! Clean up
!---------
    deallocate(Trajwrk3)
    deallocate(Trajwrk2)
    deallocate(Pertwrk)

    !Variables for pressure calculations
    deallocate(PrT_PS)
    deallocate(Ph)
    deallocate(PIh)
    deallocate(TEMP)
    deallocate(fQi)

    deallocate(ebuffery)
    deallocate(wbuffery)
    deallocate(nbufferx)
    deallocate(sbufferx)

! All done
! --------

    if (phase==TLMPhase) then
       call timing_off('DYNCORE_TLM')
    elseif (phase==ADJPhase) then
       call timing_off('DYNCORE_ADM')
    endif

    call MAPL_TimerOff(MAPL,PHASE_NAME)
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  contains

   subroutine DynCoreTestCases

    use test_cases_mod, only: init_case, test_case

! 1 - Run Jablonoski-Willaimson NL perturbation
! 2 - Run Jablonoski-Willaimson TLM perturbation
! 3 - Run Tracer perturbation with NLM
! 4 - Run Tracer perturbation with TLM

! Jablonowski-Williamson Baroclinic Instability Test Case
! -------------------------------------------------------
 
    if (DynTestCase == 1 .or. DynTestCase == 2) then
 
       !1: Generate perturbation from NL difference and run TLM
       !2: Generate perturbed trajectory and integrate it with the NLM
 
       if (firststep) then
 
          !Dont need GEOS-5 to FV3 conversions
          FV_Atm(1)%u    = 0.0_FVPRC
          FV_Atm(1)%v    = 0.0_FVPRC
          FV_Atm(1)%pt   = 0.0_FVPRC
          FV_Atm(1)%delp = 0.0_FVPRC
          FV_Atm(1)%q    = 0.0_FVPRC
          FV_Atm(1)%ps   = 0.0_FVPRC
          FV_Atm(1)%pe   = 0.0_FVPRC
          FV_Atm(1)%peln = 0.0_FVPRC
          FV_Atm(1)%pk   = 0.0_FVPRC
          FV_Atm(1)%pkz  = 0.0_FVPRC
          FV_Atm(1)%phis = 0.0_FVPRC
          FV_Atm(1)%w    = 0.0_FVPRC
          FV_Atm(1)%delz = 0.0_FVPRC
          FV_Atm(1)%uc   = 0.0_FVPRC
          FV_Atm(1)%vc   = 0.0_FVPRC
          FV_Atm(1)%ua   = 0.0_FVPRC
          FV_Atm(1)%va   = 0.0_FVPRC
   
          !Generate the unperturbed state
          test_case = 12
          call init_case(FV_Atm(1)%u,FV_Atm(1)%v,FV_Atm(1)%w,FV_Atm(1)%pt,FV_Atm(1)%delp,FV_Atm(1)%q, &
                         FV_Atm(1)%phis, FV_Atm(1)%ps,FV_Atm(1)%pe, FV_Atm(1)%peln,FV_Atm(1)%pk,FV_Atm(1)%pkz, &
                         FV_Atm(1)%uc,FV_Atm(1)%vc, FV_Atm(1)%ua,FV_Atm(1)%va,        & 
                         FV_Atm(1)%ak, FV_Atm(1)%bk, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,&
                         FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ng, &
                         ncnst, FV_Atm(1)%flagstruct%nwat,  &
                         FV_Atm(1)%flagstruct%ndims, FV_Atm(1)%flagstruct%ntiles, &
                         FV_Atm(1)%flagstruct%dry_mass, &
                         FV_Atm(1)%flagstruct%mountain,       &
                         FV_Atm(1)%flagstruct%moist_phys, FV_Atm(1)%flagstruct%hydrostatic, &
                         FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%delz, FV_Atm(1)%ze0, &
                         FV_Atm(1)%flagstruct%adiabatic, FV_Atm(1)%ks, FV_Atm(1)%neststruct%npx_global, &
                         FV_Atm(1)%ptop, FV_Atm(1)%domain, FV_Atm(1)%tile, FV_Atm(1)%bd)!, 0.0_FVPRC)

          !Create copyies of the reference state as they will be evolved by the TLM
          allocate(u_r   (isd  :ied  ,jsd  :jed+1,npz  ))
          allocate(v_r   (isd  :ied+1,jsd  :jed  ,npz  )) 
          allocate(pt_r  (isd  :ied  ,jsd  :jed  ,npz  )) 
          allocate(delp_r(isd  :ied  ,jsd  :jed  ,npz  ))
          allocate(q_r   (isd  :ied  ,jsd  :jed  ,npz  , ncnst))
          allocate(ps_r  (isd  :ied  ,jsd  :jed        ))
          allocate(pe_r  (isc-1:iec+1,npz+1,jsc-1:jec+1))
          allocate(peln_r(isc  :iec  ,npz+1,jsc  :jec  ))
          allocate(pk_r  (isc  :iec  ,jsc  :jec  ,npz+1))
          allocate(pkz_r (isc  :iec  ,jsc  :jec  ,npz  ))
          allocate(phis_r(isd  :ied  ,jsd  :jed        ))
          allocate(w_r   (isd  :ied  ,jsd  :jed  ,npz  ))
          allocate(delz_r(isd  :ied  ,jsd  :jed  ,npz  ))
          allocate(uc_r  (isd  :ied+1,jsd  :jed  ,npz  ))
          allocate(vc_r  (isd  :ied  ,jsd  :jed+1,npz  ))
          allocate(ua_r  (isd  :ied  ,jsd  :jed  ,npz  ))
          allocate(va_r  (isd  :ied  ,jsd  :jed  ,npz  ))

          u_r    = FV_Atm(1)%u
          v_r    = FV_Atm(1)%v   
          pt_r   = FV_Atm(1)%pt  
          delp_r = FV_Atm(1)%delp
          q_r    = FV_Atm(1)%q   
          ps_r   = FV_Atm(1)%ps  
          pe_r   = FV_Atm(1)%pe  
          peln_r = FV_Atm(1)%peln
          pk_r   = FV_Atm(1)%pk  
          pkz_r  = FV_Atm(1)%pkz 
          phis_r = FV_Atm(1)%phis
          w_r    = FV_Atm(1)%w   
          delz_r = FV_Atm(1)%delz
          uc_r   = FV_Atm(1)%uc  
          vc_r   = FV_Atm(1)%vc  
          ua_r   = FV_Atm(1)%ua  
          va_r   = FV_Atm(1)%va  

          FV_AtmP(1)%up    = 0.0_FVPRC
          FV_AtmP(1)%vp    = 0.0_FVPRC
          FV_AtmP(1)%ptp   = 0.0_FVPRC
          FV_AtmP(1)%delpp = 0.0_FVPRC
          FV_AtmP(1)%qp    = 0.0_FVPRC
          FV_AtmP(1)%psp   = 0.0_FVPRC
          FV_AtmP(1)%pep   = 0.0_FVPRC
          FV_AtmP(1)%pelnp = 0.0_FVPRC
          FV_AtmP(1)%pkp   = 0.0_FVPRC
          FV_AtmP(1)%pkzp  = 0.0_FVPRC
          FV_Atm(1)%phis  = 0.0_FVPRC
          FV_AtmP(1)%wp    = 0.0_FVPRC
          FV_AtmP(1)%delzp = 0.0_FVPRC
          FV_AtmP(1)%ucp   = 0.0_FVPRC
          FV_AtmP(1)%vcp   = 0.0_FVPRC
          FV_AtmP(1)%uap   = 0.0_FVPRC
          FV_AtmP(1)%vap   = 0.0_FVPRC
   
          !Generate the unperturbed state
          test_case = 13
          call init_case(FV_AtmP(1)%up,FV_AtmP(1)%vp,FV_AtmP(1)%wp,FV_AtmP(1)%ptp,FV_AtmP(1)%delpp,FV_AtmP(1)%qp, &
                         FV_Atm(1)%phis, FV_AtmP(1)%psp,FV_AtmP(1)%pep, FV_AtmP(1)%pelnp,FV_AtmP(1)%pkp,FV_AtmP(1)%pkzp, &
                         FV_AtmP(1)%ucp,FV_AtmP(1)%vcp, FV_AtmP(1)%uap,FV_AtmP(1)%vap,        & 
                         FV_Atm(1)%ak, FV_Atm(1)%bk, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,&
                         FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ng, &
                         ncnst, FV_Atm(1)%flagstruct%nwat,  &
                         FV_Atm(1)%flagstruct%ndims, FV_Atm(1)%flagstruct%ntiles, &
                         FV_Atm(1)%flagstruct%dry_mass, &
                         FV_Atm(1)%flagstruct%mountain,       &
                         FV_Atm(1)%flagstruct%moist_phys, FV_Atm(1)%flagstruct%hydrostatic, &
                         FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%delz, FV_Atm(1)%ze0, &
                         FV_Atm(1)%flagstruct%adiabatic, FV_Atm(1)%ks, FV_Atm(1)%neststruct%npx_global, &
                         FV_Atm(1)%ptop, FV_Atm(1)%domain, FV_Atm(1)%tile, FV_Atm(1)%bd)!, 0.0_FVPRC)
 
          !Test reference state only, comment save of difference
          !FV_AtmP(1)%up    = FV_Atm(1)%u   
          !FV_AtmP(1)%vp    = FV_Atm(1)%v   
          !FV_AtmP(1)%ptp   = FV_Atm(1)%pt  
          !FV_AtmP(1)%delpp = FV_Atm(1)%delp
          !FV_AtmP(1)%qp    = FV_Atm(1)%q   
          !FV_AtmP(1)%psp   = FV_Atm(1)%ps  
          !FV_AtmP(1)%pep   = FV_Atm(1)%pe  
          !FV_AtmP(1)%pelnp = FV_Atm(1)%peln
          !FV_AtmP(1)%pkp   = FV_Atm(1)%pk  
          !FV_AtmP(1)%pkzp  = FV_Atm(1)%pkz 
          !FV_Atm(1)%phis  = FV_Atm(1)%phis
          !FV_AtmP(1)%wp    = FV_Atm(1)%w   
          !FV_AtmP(1)%delzp = FV_Atm(1)%delz
          !FV_AtmP(1)%ucp   = FV_Atm(1)%uc  
          !FV_AtmP(1)%vcp   = FV_Atm(1)%vc  
          !FV_AtmP(1)%uap   = FV_Atm(1)%ua  
          !FV_AtmP(1)%vap   = FV_Atm(1)%va  

          if (DynTestCase == 1) then
             FV_AtmP(1)%up    = FV_AtmP(1)%up    - FV_Atm(1)%u   
             FV_AtmP(1)%vp    = FV_AtmP(1)%vp    - FV_Atm(1)%v   
             FV_AtmP(1)%ptp   = FV_AtmP(1)%ptp   - FV_Atm(1)%pt  
             FV_AtmP(1)%delpp = FV_AtmP(1)%delpp - FV_Atm(1)%delp
             FV_AtmP(1)%qp    = FV_AtmP(1)%qp    - FV_Atm(1)%q   
             FV_AtmP(1)%psp   = FV_AtmP(1)%psp   - FV_Atm(1)%ps  
             FV_AtmP(1)%pep   = FV_AtmP(1)%pep   - FV_Atm(1)%pe  
             FV_AtmP(1)%pelnp = FV_AtmP(1)%pelnp - FV_Atm(1)%peln
             FV_AtmP(1)%pkp   = FV_AtmP(1)%pkp   - FV_Atm(1)%pk  
             FV_AtmP(1)%pkzp  = FV_AtmP(1)%pkzp  - FV_Atm(1)%pkz 
             FV_AtmP(1)%wp    = FV_AtmP(1)%wp    - FV_Atm(1)%w   
             FV_AtmP(1)%delzp = FV_AtmP(1)%delzp - FV_Atm(1)%delz
             FV_AtmP(1)%ucp   = FV_AtmP(1)%ucp   - FV_Atm(1)%uc  
             FV_AtmP(1)%vcp   = FV_AtmP(1)%vcp   - FV_Atm(1)%vc  
             FV_AtmP(1)%uap   = FV_AtmP(1)%uap   - FV_Atm(1)%ua  
             FV_AtmP(1)%vap   = FV_AtmP(1)%vap   - FV_Atm(1)%va              
          endif

          if ( is_master() ) write(*,*) ' '
          if ( is_master() ) write(*,*) 'Jablonowski-Williamson test case initialized'
          if ( is_master() ) write(*,*) ' '

          firststep = .false.

       endif !Firststep
 
       if (DynTestCase == 1) then

          !Reset the reference state
          FV_Atm(1)%u    = u_r   
          FV_Atm(1)%v    = v_r   
          FV_Atm(1)%pt   = pt_r  
          FV_Atm(1)%delp = delp_r
          FV_Atm(1)%q    = q_r   
          FV_Atm(1)%ps   = ps_r  
          FV_Atm(1)%pe   = pe_r  
          FV_Atm(1)%peln = peln_r
          FV_Atm(1)%pk   = pk_r  
          FV_Atm(1)%pkz  = pkz_r 
          FV_Atm(1)%w    = w_r   
          FV_Atm(1)%delz = delz_r
          FV_Atm(1)%uc   = uc_r  
          FV_Atm(1)%vc   = vc_r  
          FV_Atm(1)%ua   = ua_r  
          FV_Atm(1)%va   = va_r  

          call timing_on(' FV_DYNAMICS_TLM_JWT')
!          call fv_dynamics_tlm(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
!                               DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                                    &
!                               FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
!                               cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
!                               FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
!                               FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,              &
!                               FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                              &
!                               FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                  &
!                               FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,          &
!                               FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,  &
!                               FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_AtmP(1)%q_conp, FV_Atm(1)%omga, FV_AtmP(1)%omgap,            &
!                               FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                      &
!                               FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                      &
!                               FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
!                               FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                  &
!                               FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp,                                     &
!                               FV_Atm(1)%ze0,                                                                                   &
!                               FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct,                                             &
!                               FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct, FV_Atm(1)%neststruct,                               &
!                               FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain                            ) 
          call timing_off(' FV_DYNAMICS_TLM_JWT')
 
          PertWrk(PERT_GetIndex(PertWrk,"U" ))%X = FV_AtmP(1)%up(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"V" ))%X = FV_AtmP(1)%vp(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"PT"))%X = FV_AtmP(1)%ptp(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"DP"))%X = FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"QV"))%X = FV_AtmP(1)%qp(isc:iec,jsc:jec,:,1)
          PertWrk(PERT_GetIndex(PertWrk,"O3"))%X = FV_AtmP(1)%qp(isc:iec,jsc:jec,:,2)
          PertWrk(PERT_GetIndex(PertWrk,"QL"))%X = FV_AtmP(1)%qp(isc:iec,jsc:jec,:,3)
          PertWrk(PERT_GetIndex(PertWrk,"QI"))%X = FV_AtmP(1)%qp(isc:iec,jsc:jec,:,4)

       elseif (DynTestCase == 2) then
 
          call timing_on(' FV_DYNAMICS_NLM_JWT')
          call fv_dynamics_nlm(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,               &
                           DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                                 &
                           FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                    &
                           cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                           &
                           FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                   &
                           FV_AtmP(1)%up, FV_AtmP(1)%vp, FV_AtmP(1)%wp, FV_AtmP(1)%delzp,                                    &
                           FV_Atm(1)%flagstruct%hydrostatic, FV_AtmP(1)%ptp, FV_AtmP(1)%delpp, FV_AtmP(1)%qp, FV_AtmP(1)%psp,&
                           FV_AtmP(1)%pep, FV_AtmP(1)%pkp, FV_AtmP(1)%pelnp, FV_AtmP(1)%pkzp,                                &
                           FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_AtmP(1)%omgap,                                             &
                           FV_AtmP(1)%uap, FV_AtmP(1)%vap, FV_AtmP(1)%ucp, FV_AtmP(1)%vcp,                                   &
                           FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                   &
                           FV_AtmP(1)%mfxp, FV_AtmP(1)%mfyp, FV_AtmP(1)%cxp, FV_AtmP(1)%cyp,                                 &
                           FV_Atm(1)%ze0,                                                                                &
                           FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct,                    &
                           FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain   )
          call timing_off(' FV_DYNAMICS_NLM_JWT')

          PertWrk(PERT_GetIndex(PertWrk,"U" ))%X = FV_AtmP(1)%up   (isc:iec,jsc:jec,:)   - FV_Atm(1)%u   (isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"V" ))%X = FV_AtmP(1)%vp   (isc:iec,jsc:jec,:)   - FV_Atm(1)%v   (isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"PT"))%X = FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)   - FV_Atm(1)%pt  (isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"DP"))%X = FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)   - FV_Atm(1)%delp(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"QV"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,1) - FV_Atm(1)%q   (isc:iec,jsc:jec,:,1)
          PertWrk(PERT_GetIndex(PertWrk,"O3"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,2) - FV_Atm(1)%q   (isc:iec,jsc:jec,:,2)
          PertWrk(PERT_GetIndex(PertWrk,"QL"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,3) - FV_Atm(1)%q   (isc:iec,jsc:jec,:,3)
          PertWrk(PERT_GetIndex(PertWrk,"QI"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,4) - FV_Atm(1)%q   (isc:iec,jsc:jec,:,4)

          !Test NL state only
          !PertWrk(PERT_GetIndex(PertWrk,"U" ))%X = FV_AtmP(1)%up   (isc:iec,jsc:jec,:)
          !PertWrk(PERT_GetIndex(PertWrk,"V" ))%X = FV_AtmP(1)%vp   (isc:iec,jsc:jec,:)
          !PertWrk(PERT_GetIndex(PertWrk,"PT"))%X = FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)
          !PertWrk(PERT_GetIndex(PertWrk,"DP"))%X = FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)
          !PertWrk(PERT_GetIndex(PertWrk,"QV"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,1)
          !PertWrk(PERT_GetIndex(PertWrk,"O3"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,2)
          !PertWrk(PERT_GetIndex(PertWrk,"QL"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,3)
          !PertWrk(PERT_GetIndex(PertWrk,"QI"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,4)
 
       endif
 
! Tracer Perturbation Test Case
! -----------------------------

    elseif (DynTestCase == 3 .or. DynTestCase == 4) then
 
       if (firststep) then
 
          allocate(q_r   (isd  :ied  ,jsd:jed  ,npz, ncnst), stat=status); VERIFY_(STATUS)
          q_r = 0.0
  
          allocate(lambda(isd:ied,jsd:jed), stat=status); VERIFY_(STATUS)
          allocate( theta(isd:ied,jsd:jed), stat=status); VERIFY_(STATUS)
  
          call MAPL_Get(MAPL, LATS=LATS, LONS=LONS, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
  
          lambda_c = deg2rad*TrT_lonc
           theta_c = deg2rad*TrT_latc
          lambda(isc:iec,jsc:jec) = lons !Fill internal points with lat-lon
           theta(isc:iec,jsc:jec) = lats
   
          do i = isc,iec
             do j = jsc,jec
                do k = 1,npz
                   rr = acos(sin(theta_c)*sin(theta(i,j)) + cos(theta_c)*cos(theta(i,j))*cos(lambda(i,j)-lambda_c))  
                   if (rr < 0.5) then
                      if (TrT_ver == 1) then
                         q_r(i,j,k,2) = (TrT_h0/2) * (1 + cos(pi*rr/TrT_R))
                      elseif (TrT_ver == 2) then
                         q_r(i,j,k,2) = TrT_h0
                      elseif (TrT_ver == 3) then
                         q_r(i,j,k,2) = TrT_h0
                      endif
                   endif
                   if (TrT_ver == 3) then
                      if ( abs(lambda(i,j)-lambda_c) .lt. pi/16.0 .and. theta(i,j) .gt. theta_c ) then
                         q_r(i,j,k,2) = 0.0
                      endif
                   endif
                enddo
             enddo
          enddo
  
          deallocate(lambda,theta)
  
          if ( is_master() ) write(*,*) ' '
          if ( is_master() ) write(*,*) 'Advection perturbation initialized'
          if ( is_master() ) write(*,*) ' '
  
          firststep = .false.
  
       endif
  
       if (DynTestCase == 3) then
  
          !Store the test tracer in the ozone slot
          FV_Atm(1)%q(:,:,:,2) = q_r(:,:,:,2)

          !Reset superfluous variables
          FV_Atm(1)%w    = 0.0
          FV_Atm(1)%delz = 0.0
          FV_Atm(1)%pe   = 0.0
          FV_Atm(1)%peln = 0.0
          FV_Atm(1)%pk   = 0.0
          FV_Atm(1)%pkz  = 0.0
          FV_Atm(1)%omga = 0.0
          FV_Atm(1)%uc   = 0.0
          FV_Atm(1)%vc   = 0.0
          FV_Atm(1)%ua   = 0.0
          FV_Atm(1)%va   = 0.0
          FV_Atm(1)%mfx  = 0.0
          FV_Atm(1)%mfy  = 0.0
          FV_Atm(1)%cx   = 0.0
          FV_Atm(1)%cy   = 0.0
          FV_Atm(1)%ze0  = 0.0
  
          call timing_on(' FV_DYNAMICS_NLM_TRTEST')
 
          call mpp_get_boundary( FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%domain, &
                                 wbuffery=wbuffery, ebuffery=ebuffery, &
                                 sbufferx=sbufferx, nbufferx=nbufferx, &
                                 gridtype=DGRID_NE, complete=.true. )
          do k=1,npz
             do i=isc,iec
                FV_Atm(1)%u(i,jec+1,k) = nbufferx(i,k)
             enddo
          enddo
          do k=1,npz
             do j=jsc,jec
                FV_Atm(1)%v(iec+1,j,k) = ebuffery(j,k)
             enddo
          enddo 
 
          call fv_dynamics(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                   &
                           DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                                 &
                           FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                    &
                           cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                           &
                           FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                   &
                           FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w, FV_Atm(1)%delz,                                        &
                           FV_Atm(1)%flagstruct%hydrostatic, FV_Atm(1)%pt, FV_Atm(1)%delp, FV_Atm(1)%q, FV_Atm(1)%ps,    &
                           FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,                                    &
                           FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga,                                              &
                           FV_Atm(1)%ua, FV_Atm(1)%va, FV_Atm(1)%uc, FV_Atm(1)%vc,                                       &
                           FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                   &
                           FV_Atm(1)%mfx, FV_Atm(1)%mfy, FV_Atm(1)%cx, FV_Atm(1)%cy,                                     &
                           FV_Atm(1)%ze0,                                                                                &
                           FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,                    &
                           FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain   )
  
          call timing_off(' FV_DYNAMICS_NLM_TRTEST')
  
          !Save for next timestep
          q_r = FV_Atm(1)%q
 
          !Save for output
          FV_AtmP(1)%qp(:,:,:,2) = FV_Atm(1)%q(:,:,:,2)
 
       elseif (DynTestCase == 4) then
  
          !Reset superfluous trajectory variables
          FV_Atm(1)%q(:,:,:,2:FV_Atm(1)%ncnst) = 0.0

          FV_Atm(1)%w    = 0.0
          FV_Atm(1)%delz = 0.0
          FV_Atm(1)%pe   = 0.0
          FV_Atm(1)%peln = 0.0
          FV_Atm(1)%pk   = 0.0
          FV_Atm(1)%pkz  = 0.0
          FV_Atm(1)%omga = 0.0
          FV_Atm(1)%uc   = 0.0
          FV_Atm(1)%vc   = 0.0
          FV_Atm(1)%ua   = 0.0
          FV_Atm(1)%va   = 0.0
          FV_Atm(1)%mfx  = 0.0
          FV_Atm(1)%mfy  = 0.0
          FV_Atm(1)%cx   = 0.0
          FV_Atm(1)%cy   = 0.0
          FV_Atm(1)%ze0  = 0.0

          !Reset superfluous perturbation variables
          FV_AtmP(1)%up    = 0.0
          FV_AtmP(1)%vp    = 0.0
          FV_AtmP(1)%ptp   = 0.0
          FV_AtmP(1)%delpp = 0.0
          FV_AtmP(1)%qp    = q_r

          FV_AtmP(1)%wp    = 0.0
          FV_AtmP(1)%delzp = 0.0
          FV_AtmP(1)%pep   = 0.0
          FV_AtmP(1)%pelnp = 0.0
          FV_AtmP(1)%pkp   = 0.0
          FV_AtmP(1)%pkzp  = 0.0
          FV_AtmP(1)%omgap = 0.0
          FV_AtmP(1)%ucp   = 0.0
          FV_AtmP(1)%vcp   = 0.0
          FV_AtmP(1)%uap   = 0.0
          FV_AtmP(1)%vap   = 0.0
          FV_AtmP(1)%mfxp  = 0.0
          FV_AtmP(1)%mfyp  = 0.0
          FV_AtmP(1)%cxp   = 0.0
          FV_AtmP(1)%cyp   = 0.0
  
          call timing_on(' FV_DYNAMICS_TLM_TRTEST')
  
          call mpp_get_boundary( FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%domain, &
                                 wbuffery=wbuffery, ebuffery=ebuffery, &
                                 sbufferx=sbufferx, nbufferx=nbufferx, &
                                 gridtype=DGRID_NE, complete=.true. )
          do k=1,npz
             do i=isc,iec
                FV_Atm(1)%u(i,jec+1,k) = nbufferx(i,k)
             enddo
          enddo
          do k=1,npz
             do j=jsc,jec
                FV_Atm(1)%v(iec+1,j,k) = ebuffery(j,k)
             enddo
          enddo 

        !  call fv_dynamics_tlm(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
        !                       DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                                    &
        !                       FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
        !                       cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
        !                       FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
        !                       FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,              &
        !                       FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                              &
        !                       FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                  &
        !                       FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,          &
        !                       FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,  &
        !                       FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_AtmP(1)%q_conp, FV_Atm(1)%omga, FV_AtmP(1)%omgap,            &
        !                       FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                      &
        !                       FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                      &
        !                       FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
        !                       FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                  &
        !                       FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp,                                     &
        !                       FV_Atm(1)%ze0,                                                                                   &
        !                       FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct,                                             &
        !                       FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct, FV_Atm(1)%neststruct,                               &
        !                       FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain                            )
  
          call timing_off(' FV_DYNAMICS_TLM_TRTEST')
  
          !Save for next timestep
          q_r = FV_AtmP(1)%qp
  
       endif

! Gravity wave perturbation test case
! -----------------------------------

    elseif (DynTestCase == 5 .or. DynTestCase == 6) then
 
       !1: Generate perturbation from NL difference and run TLM
       !2: Generate perturbed trajectory and integrate it with the NLM
 
       if (firststep) then
 
          !Dont need GEOS-5 to FV3 conversions
        !  TLMIdealTest = .true.
           
          FV_Atm(1)%u    = 0.0_FVPRC
          FV_Atm(1)%v    = 0.0_FVPRC
          FV_Atm(1)%pt   = 0.0_FVPRC
          FV_Atm(1)%delp = 0.0_FVPRC
          FV_Atm(1)%q    = 0.0_FVPRC
          FV_Atm(1)%ps   = 0.0_FVPRC
          FV_Atm(1)%pe   = 0.0_FVPRC
          FV_Atm(1)%peln = 0.0_FVPRC
          FV_Atm(1)%pk   = 0.0_FVPRC
          FV_Atm(1)%pkz  = 0.0_FVPRC
          FV_Atm(1)%phis = 0.0_FVPRC
          FV_Atm(1)%w    = 0.0_FVPRC
          FV_Atm(1)%delz = 0.0_FVPRC
          FV_Atm(1)%uc   = 0.0_FVPRC
          FV_Atm(1)%vc   = 0.0_FVPRC
          FV_Atm(1)%ua   = 0.0_FVPRC
          FV_Atm(1)%va   = 0.0_FVPRC
   
          !Generate the unperturbed state
          test_case = 16
          call init_case(FV_Atm(1)%u,FV_Atm(1)%v,FV_Atm(1)%w,FV_Atm(1)%pt,FV_Atm(1)%delp,FV_Atm(1)%q, &
                         FV_Atm(1)%phis, FV_Atm(1)%ps,FV_Atm(1)%pe, FV_Atm(1)%peln,FV_Atm(1)%pk,FV_Atm(1)%pkz, &
                         FV_Atm(1)%uc,FV_Atm(1)%vc, FV_Atm(1)%ua,FV_Atm(1)%va,        & 
                         FV_Atm(1)%ak, FV_Atm(1)%bk, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,&
                         FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ng, &
                         ncnst, FV_Atm(1)%flagstruct%nwat,  &
                         FV_Atm(1)%flagstruct%ndims, FV_Atm(1)%flagstruct%ntiles, &
                         FV_Atm(1)%flagstruct%dry_mass, &
                         FV_Atm(1)%flagstruct%mountain,       &
                         FV_Atm(1)%flagstruct%moist_phys, FV_Atm(1)%flagstruct%hydrostatic, &
                         FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%delz, FV_Atm(1)%ze0, &
                         FV_Atm(1)%flagstruct%adiabatic, FV_Atm(1)%ks, FV_Atm(1)%neststruct%npx_global, &
                         FV_Atm(1)%ptop, FV_Atm(1)%domain, FV_Atm(1)%tile, FV_Atm(1)%bd)!, 0.0_FVPRC)

          !Create copyies of the reference state as they will be evolved by the TLM
          allocate(u_r   (isd  :ied  ,jsd  :jed+1,npz  ))
          allocate(v_r   (isd  :ied+1,jsd  :jed  ,npz  )) 
          allocate(pt_r  (isd  :ied  ,jsd  :jed  ,npz  )) 
          allocate(delp_r(isd  :ied  ,jsd  :jed  ,npz  ))
          allocate(q_r   (isd  :ied  ,jsd  :jed  ,npz  , ncnst))
          allocate(ps_r  (isd  :ied  ,jsd  :jed        ))
          allocate(pe_r  (isc-1:iec+1,npz+1,jsc-1:jec+1))
          allocate(peln_r(isc  :iec  ,npz+1,jsc  :jec  ))
          allocate(pk_r  (isc  :iec  ,jsc  :jec  ,npz+1))
          allocate(pkz_r (isc  :iec  ,jsc  :jec  ,npz  ))
          allocate(phis_r(isd  :ied  ,jsd  :jed        ))
          allocate(w_r   (isd  :ied  ,jsd  :jed  ,npz  ))
          allocate(delz_r(isd  :ied  ,jsd  :jed  ,npz  ))
          allocate(uc_r  (isd  :ied+1,jsd  :jed  ,npz  ))
          allocate(vc_r  (isd  :ied  ,jsd  :jed+1,npz  ))
          allocate(ua_r  (isd  :ied  ,jsd  :jed  ,npz  ))
          allocate(va_r  (isd  :ied  ,jsd  :jed  ,npz  ))
          allocate(ak_r  ( npz+1 ))
          allocate(bk_r  ( npz+1 ))

          u_r    = FV_Atm(1)%u
          v_r    = FV_Atm(1)%v   
          pt_r   = FV_Atm(1)%pt  
          delp_r = FV_Atm(1)%delp
          q_r    = FV_Atm(1)%q   
          ps_r   = FV_Atm(1)%ps  
          pe_r   = FV_Atm(1)%pe  
          peln_r = FV_Atm(1)%peln
          pk_r   = FV_Atm(1)%pk  
          pkz_r  = FV_Atm(1)%pkz 
          phis_r = FV_Atm(1)%phis
          w_r    = FV_Atm(1)%w   
          delz_r = FV_Atm(1)%delz
          uc_r   = FV_Atm(1)%uc  
          vc_r   = FV_Atm(1)%vc  
          ua_r   = FV_Atm(1)%ua  
          va_r   = FV_Atm(1)%va  

          !Save grid incase gravity wave test
          ak_r   = FV_Atm(1)%ak
          bk_r   = FV_Atm(1)%bk
          ptop_r = FV_Atm(1)%ptop

          FV_AtmP(1)%up    = 0.0_FVPRC
          FV_AtmP(1)%vp    = 0.0_FVPRC
          FV_AtmP(1)%ptp   = 0.0_FVPRC
          FV_AtmP(1)%delpp = 0.0_FVPRC
          FV_AtmP(1)%qp    = 0.0_FVPRC
          FV_AtmP(1)%psp   = 0.0_FVPRC
          FV_AtmP(1)%pep   = 0.0_FVPRC
          FV_AtmP(1)%pelnp = 0.0_FVPRC
          FV_AtmP(1)%pkp   = 0.0_FVPRC
          FV_AtmP(1)%pkzp  = 0.0_FVPRC
          FV_Atm(1)%phis  = 0.0_FVPRC
          FV_AtmP(1)%wp    = 0.0_FVPRC
          FV_AtmP(1)%delzp = 0.0_FVPRC
          FV_AtmP(1)%ucp   = 0.0_FVPRC
          FV_AtmP(1)%vcp   = 0.0_FVPRC
          FV_AtmP(1)%uap   = 0.0_FVPRC
          FV_AtmP(1)%vap   = 0.0_FVPRC
   
          !Generate the unperturbed state
          test_case = 16
          call init_case(FV_AtmP(1)%up,FV_AtmP(1)%vp,FV_AtmP(1)%wp,FV_AtmP(1)%ptp,FV_AtmP(1)%delpp,FV_AtmP(1)%qp, &
                         FV_Atm(1)%phis, FV_AtmP(1)%psp,FV_AtmP(1)%pep, FV_AtmP(1)%pelnp,FV_AtmP(1)%pkp,FV_AtmP(1)%pkzp, &
                         FV_AtmP(1)%ucp,FV_AtmP(1)%vcp, FV_AtmP(1)%uap,FV_AtmP(1)%vap,        & 
                         FV_Atm(1)%ak, FV_Atm(1)%bk, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,&
                         FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ng, &
                         ncnst, FV_Atm(1)%flagstruct%nwat,  &
                         FV_Atm(1)%flagstruct%ndims, FV_Atm(1)%flagstruct%ntiles, &
                         FV_Atm(1)%flagstruct%dry_mass, &
                         FV_Atm(1)%flagstruct%mountain,       &
                         FV_Atm(1)%flagstruct%moist_phys, FV_Atm(1)%flagstruct%hydrostatic, &
                         FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%delz, FV_Atm(1)%ze0, &
                         FV_Atm(1)%flagstruct%adiabatic, FV_Atm(1)%ks, FV_Atm(1)%neststruct%npx_global, &
                         FV_Atm(1)%ptop, FV_Atm(1)%domain, FV_Atm(1)%tile, FV_Atm(1)%bd)!, 1.0_FVPRC)
 
          !Test reference state only, comment save of difference
          !FV_AtmP(1)%up    = FV_Atm(1)%u   
          !FV_AtmP(1)%vp    = FV_Atm(1)%v   
          !FV_AtmP(1)%ptp   = FV_Atm(1)%pt  
          !FV_AtmP(1)%delpp = FV_Atm(1)%delp
          !FV_AtmP(1)%qp    = FV_Atm(1)%q   
          !FV_AtmP(1)%psp   = FV_Atm(1)%ps  
          !FV_AtmP(1)%pep   = FV_Atm(1)%pe  
          !FV_AtmP(1)%pelnp = FV_Atm(1)%peln
          !FV_AtmP(1)%pkp   = FV_Atm(1)%pk  
          !FV_AtmP(1)%pkzp  = FV_Atm(1)%pkz 
          !FV_Atm(1)%phis  = FV_Atm(1)%phis
          !FV_AtmP(1)%wp    = FV_Atm(1)%w   
          !FV_AtmP(1)%delzp = FV_Atm(1)%delz
          !FV_AtmP(1)%ucp   = FV_Atm(1)%uc  
          !FV_AtmP(1)%vcp   = FV_Atm(1)%vc  
          !FV_AtmP(1)%uap   = FV_Atm(1)%ua  
          !FV_AtmP(1)%vap   = FV_Atm(1)%va  

          if (DynTestCase == 5) then
             FV_AtmP(1)%up    = FV_AtmP(1)%up    - FV_Atm(1)%u   
             FV_AtmP(1)%vp    = FV_AtmP(1)%vp    - FV_Atm(1)%v   
             FV_AtmP(1)%ptp   = FV_AtmP(1)%ptp   - FV_Atm(1)%pt  
             FV_AtmP(1)%delpp = FV_AtmP(1)%delpp - FV_Atm(1)%delp
             FV_AtmP(1)%qp    = FV_AtmP(1)%qp    - FV_Atm(1)%q   
             FV_AtmP(1)%psp   = FV_AtmP(1)%psp   - FV_Atm(1)%ps  
             FV_AtmP(1)%pep   = FV_AtmP(1)%pep   - FV_Atm(1)%pe  
             FV_AtmP(1)%pelnp = FV_AtmP(1)%pelnp - FV_Atm(1)%peln
             FV_AtmP(1)%pkp   = FV_AtmP(1)%pkp   - FV_Atm(1)%pk  
             FV_AtmP(1)%pkzp  = FV_AtmP(1)%pkzp  - FV_Atm(1)%pkz 
             FV_AtmP(1)%wp    = FV_AtmP(1)%wp    - FV_Atm(1)%w   
             FV_AtmP(1)%delzp = FV_AtmP(1)%delzp - FV_Atm(1)%delz
             FV_AtmP(1)%ucp   = FV_AtmP(1)%ucp   - FV_Atm(1)%uc  
             FV_AtmP(1)%vcp   = FV_AtmP(1)%vcp   - FV_Atm(1)%vc  
             FV_AtmP(1)%uap   = FV_AtmP(1)%uap   - FV_Atm(1)%ua  
             FV_AtmP(1)%vap   = FV_AtmP(1)%vap   - FV_Atm(1)%va              
          endif

          if ( is_master() ) write(*,*) ' '
          if ( is_master() ) write(*,*) 'Test case initialized'
          if ( is_master() ) write(*,*) ' '

          firststep = .false.

       endif !Firststep
 
       if (DynTestCase == 5) then

          !Reset the reference state
          FV_Atm(1)%u    = u_r   
          FV_Atm(1)%v    = v_r   
          FV_Atm(1)%pt   = pt_r  
          FV_Atm(1)%delp = delp_r
          FV_Atm(1)%q    = q_r   
          FV_Atm(1)%ps   = ps_r  
          FV_Atm(1)%pe   = pe_r  
          FV_Atm(1)%peln = peln_r
          FV_Atm(1)%pk   = pk_r  
          FV_Atm(1)%pkz  = pkz_r 
          FV_Atm(1)%w    = w_r   
          FV_Atm(1)%delz = delz_r
          FV_Atm(1)%uc   = uc_r  
          FV_Atm(1)%vc   = vc_r  
          FV_Atm(1)%ua   = ua_r  
          FV_Atm(1)%va   = va_r  

          FV_Atm(1)%ak   = ak_r
          FV_Atm(1)%bk   = bk_r
          FV_Atm(1)%ptop = ptop_r

          call timing_on(' FV_DYNAMICS_TLM_JWT')
     !     call fv_dynamics_tlm(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
     !                          DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                                    &
     !                          FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
     !                          cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
     !                          FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
     !                          FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,              &
     !                          FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                              &
     !                          FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                  &
     !                          FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,          &
     !                          FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,  &
     !                          FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_AtmP(1)%q_conp, FV_Atm(1)%omga, FV_AtmP(1)%omgap,            &
     !                          FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                      &
     !                          FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                      &
     !                          FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
     !                          FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                  &
     !                          FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp,                                     &
     !                          FV_Atm(1)%ze0,                                                                                   &
     !                          FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct,                                             &
     !                          FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct, FV_Atm(1)%neststruct,                               &
     !                          FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain                            ) 
          call timing_off(' FV_DYNAMICS_TLM_JWT')
 
          PertWrk(PERT_GetIndex(PertWrk,"U" ))%X = FV_AtmP(1)%up(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"V" ))%X = FV_AtmP(1)%vp(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"PT"))%X = FV_AtmP(1)%ptp(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"DP"))%X = FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"QV"))%X = FV_AtmP(1)%qp(isc:iec,jsc:jec,:,1)
          PertWrk(PERT_GetIndex(PertWrk,"O3"))%X = FV_AtmP(1)%qp(isc:iec,jsc:jec,:,2)
          PertWrk(PERT_GetIndex(PertWrk,"QL"))%X = FV_AtmP(1)%qp(isc:iec,jsc:jec,:,3)
          PertWrk(PERT_GetIndex(PertWrk,"QI"))%X = FV_AtmP(1)%qp(isc:iec,jsc:jec,:,4)

       elseif (DynTestCase == 6) then
 
          FV_Atm(1)%ak   = ak_r
          FV_Atm(1)%bk   = bk_r
          FV_Atm(1)%ptop = ptop_r

          call timing_on(' FV_DYNAMICS_TESTCASE')
!          call tap_fv_dynamics(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                        &
!                               DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                                      &
!                               FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                         &
!                               cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                                &
!                               FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                        &
!                               FV_AtmP(1)%up, FV_AtmP(1)%vp, FV_AtmP(1)%wp, FV_AtmP(1)%delzp,                                     &
!                               FV_Atm(1)%flagstruct%hydrostatic, FV_AtmP(1)%ptp, FV_AtmP(1)%delpp, FV_AtmP(1)%qp, FV_AtmP(1)%psp, &
!                               FV_AtmP(1)%pep, FV_AtmP(1)%pkp, FV_AtmP(1)%pelnp, FV_AtmP(1)%pkzp,                                 &
!                               FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_AtmP(1)%omgap,                                                 &
!                               FV_AtmP(1)%uap, FV_AtmP(1)%vap, FV_AtmP(1)%ucp, FV_AtmP(1)%vcp,                                    &
!                               FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                        &
!                               FV_AtmP(1)%mfxp, FV_AtmP(1)%mfyp, FV_AtmP(1)%cxp, FV_AtmP(1)%cyp,                                  &
!                               FV_Atm(1)%ze0,                                                                                     &
!                               FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct,  &
!                               FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain        )
          call timing_off(' FV_DYNAMICS_TESTCASE')

          PertWrk(PERT_GetIndex(PertWrk,"U" ))%X = FV_AtmP(1)%up   (isc:iec,jsc:jec,:)   - FV_Atm(1)%u   (isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"V" ))%X = FV_AtmP(1)%vp   (isc:iec,jsc:jec,:)   - FV_Atm(1)%v   (isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"PT"))%X = FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)   - FV_Atm(1)%pt  (isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"DP"))%X = FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)   - FV_Atm(1)%delp(isc:iec,jsc:jec,:)
          PertWrk(PERT_GetIndex(PertWrk,"QV"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,1) - FV_Atm(1)%q   (isc:iec,jsc:jec,:,1)
          PertWrk(PERT_GetIndex(PertWrk,"O3"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,2) - FV_Atm(1)%q   (isc:iec,jsc:jec,:,2)
          PertWrk(PERT_GetIndex(PertWrk,"QL"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,3) - FV_Atm(1)%q   (isc:iec,jsc:jec,:,3)
          PertWrk(PERT_GetIndex(PertWrk,"QI"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,4) - FV_Atm(1)%q   (isc:iec,jsc:jec,:,4)

          !Test NL state only
          !PertWrk(PERT_GetIndex(PertWrk,"U" ))%X = FV_AtmP(1)%up   (isc:iec,jsc:jec,:)
          !PertWrk(PERT_GetIndex(PertWrk,"V" ))%X = FV_AtmP(1)%vp   (isc:iec,jsc:jec,:)
          !PertWrk(PERT_GetIndex(PertWrk,"PT"))%X = FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)
          !PertWrk(PERT_GetIndex(PertWrk,"DP"))%X = FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)
          !PertWrk(PERT_GetIndex(PertWrk,"QV"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,1)
          !PertWrk(PERT_GetIndex(PertWrk,"O3"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,2)
          !PertWrk(PERT_GetIndex(PertWrk,"QL"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,3)
          !PertWrk(PERT_GetIndex(PertWrk,"QI"))%X = FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,4)
 
       endif

    endif !Test cases

   end subroutine DynCoreTestCases

  end subroutine Run

  subroutine Finalize ( gc, import, export, clock, rc )

    use fv_control_mod, only : fv_end

! !ARGUMENTS:
    
    type(ESMF_GridComp), intent(inout) :: gc       ! composite gridded component 
    type(ESMF_State),    intent(inout) :: import   ! import state
    type(ESMF_State),    intent(inout) :: export   ! export state
    type(ESMF_Clock),    intent(inout) :: clock    ! the clock
    integer, optional,   intent(  out) :: rc       ! Error code:

! Locals
!-------
    type (MAPL_MetaComp),      pointer :: MAPL 

    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: IAm
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Begin
!------

    Iam = "Finalize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"FINALIZE")

! Deallocate saved DynCorePert variables
! --------------------------------------

    !Dyn Core Test Case Variables
    if (allocated(   u_r)) deallocate(   u_r)
    if (allocated(   v_r)) deallocate(   v_r)
    if (allocated(  pt_r)) deallocate(  pt_r)
    if (allocated(delp_r)) deallocate(delp_r)
    if (allocated(   q_r)) deallocate(   q_r)
    if (allocated(  ps_r)) deallocate(  ps_r)
    if (allocated(  pe_r)) deallocate(  pe_r)
    if (allocated(peln_r)) deallocate(peln_r)
    if (allocated(  pk_r)) deallocate(  pk_r)
    if (allocated( pkz_r)) deallocate( pkz_r)
    if (allocated(phis_r)) deallocate(phis_r)
    if (allocated(   w_r)) deallocate(   w_r)
    if (allocated(delz_r)) deallocate(delz_r)
    if (allocated(  uc_r)) deallocate(  uc_r)
    if (allocated(  vc_r)) deallocate(  vc_r)
    if (allocated(  ua_r)) deallocate(  ua_r)
    if (allocated(  va_r)) deallocate(  va_r)
    if (allocated(  ak_r)) deallocate(  ak_r)
    if (allocated(  bk_r)) deallocate(  bk_r)

! Retrieve the pointer to the internal state
! ------------------------------------------
    call fv_end_pert(FV_AtmP)
    deallocate(FV_AtmP)

    call fv_end(FV_Atm, grids_on_this_pe, .false.)
    deallocate(FV_Atm)

! Deallocate memory used by internal dyn core checkpoints
! -------------------------------------------------------

    if (cp_iter_controls%cp_i .ne. 0 .and. .not.cp_iter(cp_dyn_ind)%cp_test) then
       if (is_master()) write(*,*), 'Deallocating FV checkpoints'
       call cp_mod_end
    endif

! Call Generic Finalize
! ---------------------

    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL")

    call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call WRITE_PARALLEL("DynPert: Finalize Stubbed")

    RETURN_(ESMF_SUCCESS)

  end subroutine Finalize

  subroutine register_grid_and_regridders()
    use MAPL_GridManagerMod, only: grid_manager
    use CubedSphereGridFactoryMod, only: CubedSphereGridFactory
    use MAPL_RegridderManagerMod, only: regridder_manager
    use MAPL_RegridderSpecMod, only: REGRID_METHOD_BILINEAR
    use LatLonToCubeRegridderMod
    use CubeToLatLonRegridderMod
    use CubeToCubeRegridderMod

    type (CubedSphereGridFactory) :: factory

    type (CubeToLatLonRegridder) :: cube_to_latlon_prototype
    type (LatLonToCubeRegridder) :: latlon_to_cube_prototype
    type (CubeToCubeRegridder) :: cube_to_cube_prototype

    call grid_manager%add_prototype('Cubed-Sphere',factory)
    associate (method => REGRID_METHOD_BILINEAR, mgr => regridder_manager)
      call mgr%add_prototype('Cubed-Sphere', 'LatLon', method, cube_to_latlon_prototype)
      call mgr%add_prototype('LatLon', 'Cubed-Sphere', method, latlon_to_cube_prototype)
      call mgr%add_prototype('Cubed-Sphere', 'Cubed-Sphere', method, cube_to_cube_prototype)
    end associate

  end subroutine register_grid_and_regridders

end module GEOS_DynCorePertGridCompMod
