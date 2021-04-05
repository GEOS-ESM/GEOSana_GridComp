#include "MAPL_ErrLog.h" 
#ifdef _REAL8_
#define _GMAO_FVGSI_
#endif
!#define PRINT_STATES
!#define SFCverbose
!#define UPAverbose
!#define VERBOSE
!#define DEBUG_TRACE
#include "mytrace.H"
!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOP

! !MODULE: GSI_GridCompMod -- Implements ESMF wrapper to invoke GSI
!
! !DESCRIPTION: 
! {\tt GSI\_GridComp} is the ESMF Gridded Component module for NCEP's GSI
! It defines the ESMF Initialize/Run/Finalize methods as well as ancillary routines.
! {\tt GSI\_GridComp} acts as a wrapper around NCEP's GSI. The GSI internal fields
! are defined from the component's import state which are, on import,
! already decomposed as required by the GSI.
! An important function of the {\tt GSI\_GridComp} is to take the import fields
! and use them to populate the GSI internal fields. Since these latter are defined
! on a different grid and with different units some additional manipulations
! are required. These require the need to access the GSI's global data segments
! guess_grids, grid_mod, and mpimod.
!

!
! !INTERFACE:

   module GSI_GridCompMod

! !USES:

   use ESMF                  ! ESMF
   use MAPL_Mod              ! MAPL Generic
   use gsimod                ! GSI original interface
   use mpeu_util, only: StrTemplate         ! grads style templates
   use m_chars,   only: lowercase
   use m_chars,   only: uppercase
   use constants, only : grav,fv,t0c,r0_05
   use m_tick, only: tick
   use obsmod, only: lobserver

   use gsi_bundlemod, only : GSI_BundleGetPointer
   use gsi_bundlemod, only : GSI_Bundle
   use gsi_bundlemod, only : GSI_BundlePrint
   use MAPL_LatLonGridFactoryMod
   use MAPL_GridManagerMod

! Access to GSI's global data segments

   ! guess fields and related variables

   	!>>> sfc_grids <<< all of them
   use guess_grids, only: isli,         & ! snow/land/ice mask
                          fact10,       & ! 10 meter wind factor
                          sfct,         & ! guess skin temp
                          dsfct,        & ! delta guess skin temp
                          sno,          & ! snow-ice mask
                          veg_type,     & ! vegetation type  
                          veg_frac,     & ! vegetation frac
                          soil_type,    & ! soil type
                          soil_temp,    & ! soil temperature
                          soil_moi,     & ! soil moisture
                          sfc_rough       ! surface roughness

   use guess_grids, only: nfldsig,      & ! number of guess sigma times
                          nfldsfc,      & ! number of guess surface times
                          nfldnst,      & ! 
                          hrdifsig,     & ! times for guess sigma fields
                          hrdifnst,     & ! times for guess nst fields 
                          hrdifsfc        ! times for guess surface fields 
                                          ! atmosphere in kPa
   use guess_grids, only: ntguessig,    & ! These variables are used  
                          ntguessfc,    & ! to handle file names and are
                          ntguesnst,    & ! nst analysis time handle
                          ifilenst,     & ! included for compatibility
                          ifilesfc,     & ! included for compatibility
                          ifilesig        ! with the legacy code.

   use guess_grids, only: nfldsig_all,	& ! all background times count
   			  nfldsfc_all,	&
   			  nfldnst_all,	& ! 
			  nfldsig_now,	& ! filled guess_grids count
			  nfldsfc_now,	&
			  nfldnst_now,	& ! 
   			  hrdifsig_all,	& ! hour list of all backgrounds
			  hrdifsfc_all,	&
			  hrdifnst_all,	&
			  extrap_intime	  ! do time-extrapolation or not

   use guess_grids, only: ntguessig_ref	! a fixed reference ntguessig value
   use guess_grids, only: ntguessfc_ref	! a fixed reference ntguessfc value
   use guess_grids, only: ntguesnst_ref	! a fixed reference ntguesnst value

   use guess_grids, only: guess_grids_stats

   use guess_grids, only: create_chemges_grids, &
                          destroy_chemges_grids, &
                          create_metguess_grids, &
                          destroy_metguess_grids

   ! routines from gridmod
   use gridmod,   only : create_grid_vars,   &
                         destroy_grid_vars,  &
   ! variables for create_grid_vars
                         rlats,    & ! grid latitudes (radians)
                         rlons,    & ! grid longitudes (radians)
                         coslon,   & ! cos(grid longitudes (radians))
                         sinlon,   & ! sin(grid longitudes (radians))
                         wgtlats,  & ! gaussian integration weights
                         rbs2,     & ! 1./sin(grid latitudes))**2
   ! variables for create_mapping and deter_subdomain
                         nlat_sfc, & ! no. of latitudes surface files
                         nlon_sfc, & ! no. of longitudes surface files
                         nlat,     & ! no. of analysis grid latitudes
                         nlon,     & ! no. of analysis grid longitudes
                         nsig,     & ! no. of levels
                         istart,   & ! start lat of the whole array on each pe
                         jstart,   & ! start lon of the whole array on each pe
                         ilat1,    & ! no. of lats for each subdomain (no buffer)
                         jlon1,    & ! no. of lons for each subdomain (no buffer)
                         lat1,     & ! no. of lats on subdomain (no buffer)
                         lon1,     & ! no. of lons on subdomain (no buffer)
                         lat2,     & ! lat2: no. of lats on subdomain 
                                     ! (buffer points on ends)
                         lon2,     & ! lon2: no. of lons on subdomain 
                                     ! (buffer points on ends)
                         ak5,bk5,ck5, & ! coefficients for hybrid vertical coordinate
                         idvc5,    & ! vertical coordinate identifier
   ! spectral transform grid info
                         sp_a

   ! from mpimod
   use mpimod,     only: npe,nxpe,nype,  & ! num of MPI tasks (total, along x, along y)
                                           ! (nxpe, nype new for ESMF)
                         mpi_comm_world    ! MPI communicator

   use gsi_4dvar, only: time_4dvar
   use gsi_4dvar, only: ibdate, iedate, iadatebgn, iadateend, iwinbgn, &
                         nhr_assimilation,& ! size of assimilation window (hrs)
                         min_offset,&       ! offset minutes from analysis time
                         iwrtinc,&          ! when .t., writes out increment
			 idmodel,&          ! use identity perturbation model
                         l4dvar, &          ! when .t., will run 4d-var
                         tau_fcst           ! interface of forecast (error; wrt to ana date/time)

   ! others...
   use constants, only: pi, rearth, zero, one, half
   use kinds,     only: r_kind,r_single,i_kind,r_double
   use obsmod,    only: iadate, ianldate, ndat, ndat_times, time_offset

   ! meteorological guess
!  use gsi_metguess_mod, only: gsi_metguess_create_grids
!  use gsi_metguess_mod, only: gsi_metguess_destroy_grids
   use gsi_metguess_mod, only: gsi_metguess_get
   use gsi_metguess_mod, only: gsi_metguess_bundle
   use gsi_metguess_mod, only: gsi_metguess_init

   ! chem trace gases
!  use gsi_chemguess_mod, only: gsi_chemguess_create_grids
!  use gsi_chemguess_mod, only: gsi_chemguess_destroy_grids
   use gsi_chemguess_mod, only: gsi_chemguess_get
   use gsi_chemguess_mod, only: gsi_chemguess_bundle
   use gsi_chemguess_mod, only: gsi_chemguess_init

   implicit none

   private

! !PUBLIC ROUTINES:

   public GSI_GridCompSetServices
   public GSI_GridCompSetupSpecs
   public GSI_GridCompSetAlarms
   public GSI_clockShow

   public GSI_ExpId
   public GSI_AgcmPertGrid
   public GSI_aClock
   public PERTMOD_RUN_DT
   public GEOS_MU_SKIN
   public GSI_HALOWIDTH
   public GSI_Time0
   public GSI_RefTime
   public GSI_FcsTime
   public GSI_bkg_fname_tmpl
   public GSI_ensbkg_fname_tmpl
   public GSI_ensana_fname_tmpl
   public GSI_ensprgA_fname_tmpl
   public GSI_ensprgB_fname_tmpl
   public GSI_fsens_fname_tmpl
   public GSI_ferrA_fname_tmpl
   public GSI_ferrB_fname_tmpl
   public GSI_ensread_blocksize

   public PPMV2GpG
!
! !REVISION HISTORY:
!
!   19Mar2007 Todling  Removed sigi,sigl; added module w/ 4dvar timers
!   05Apr2007 CCruz/RT Handle for many background fields
!   14Apr2007 Todling  Calculating idate4 internally; removed from rc file
!   17Apr2007 Todling  Added swap of vertical levels as internal feature;
!                      changed reference times in hrdifsig/sfc - per YT
!   25Apr2007 Todling  Added ts and lwi as export for FileSpec compliance
!   08May2007 Todling  Removed reference to global surface arrays
!   05Jul2007 Todling  Analysis time set in rc file (defining ntguess)
!   10Jul2007 Todling  Adjustment to cope w/ write out of increment
!   08Jul2008 Todling  Merge fdda-b1p3 w/ das-215 (MAPL update)
!   18Nov2008 Todling  Merge with NCEP-May-2008; add sfc_rough
!   09Oct2009 Wu       replace nhr_offset with min_offset since it's 1.5 hr for regional
!   10Jul2009 Todling  Remove vertical halo
!   19Aug2009 Guo      Added hrdifsig_all etc. and other changes for
!		       multi-pass observer.
!   31Mar2010 Treadon  replace specmod routines with general_specmod
!   21Apr2010 Todling  - Add chem tracers capability
!                      - Rename Initalize/Run/Finalize
!   20May2010 Todling  Change initialization of chem
!   11Feb2011 Todling  Add GOCART aerosols to import state
!                      REMARK: for some reason this produces roundoff diffs
!   07Apr2011 Guo      Initialized all tracer pointers to null();
!                      Added a loop in Scale_import() for undefined tracer ptrs;
!                      Added additional parameters for GEOSgsi_Coupler;
!                      Reordered the reference to the gsi_4dcoupler.
!   20May2011 Guo      - Moved alarm setting ahead of gsi_4dcouple_init_traj().
!		       - Added a public interface GSI_clockShow() to desplay
!			 clock information.
!   05Jul2011 Todling  Adjustments to run GSI(r_single); back to using GSI-div/vor
!   15Apr2012 Todling  - GMAO-GSI now must have use_sp_eqspace set to true in namelist
!                      - Remove reference to vars not used here
!   18Sep2012 Akella   - Add SwapJI2i_ to GSI_GridCompSwapJI_ for use in writing out TSKIN when doing NST analysis
!   19Oct2013 Todling - metguess now holds background
!
!EOP
!-------------------------------------------------------------------------

! GLOBAL scope vars

   character(len=*),parameter:: myname='GSI_GridComp'

   ! parameters
   integer, parameter  :: HW = 1    ! horizontal halowidth of distributed GSI fields
   real, parameter     :: UNDEF_SSI_      = -9.9899991E33 
   real, parameter     :: UNDEF_SOIL_TEMP = 1.E+15
   real, parameter     :: UNDEF_SNOW_DEP  = 1.E+12
   real, parameter     :: PPMV2GpG        = 1.6571E-6         ! ((47.9982 g/mol)/((28.9644 g/mol))*1e-6(mol/mol)-> g/g
   real, parameter     :: KGpKG2PPBV      = (28./28.97)*1.E+9 ! mol/mol to ppbv
   real, parameter     :: KGpKG2PPBVaero  = 1.E+9             ! kg/kg to ppbv (mol/mol to ppbv)
   real, parameter     :: KGpKG2ppmvCO2   = 1.E+6             ! kg/kg to ppmv (mol/mol to ppmv)
   real, parameter     :: kPa_per_Pa      = 0.001
   real, parameter     :: Pa_per_kPa      = 1000.
   logical, parameter  :: verbose         = .false.

   integer, save       :: nbkgfreq        = 1
   integer, save       :: nthTimeIndex    = 0
   integer, save       :: BKGfreq_hr,BKGfreq_mn,BKGfreq_sc
   integer, save       :: ANAfreq_hr,ANAfreq_mn,ANAfreq_sc
   integer, save       :: ANAwndw_hr
   integer(i_kind), save :: MYIDATE(4)
   integer(i_kind), save :: MYHOURG = 0.
   real(r_kind), save  :: BKGfrac_hr

   logical, save       :: doVflip = .true.
   character(len=ESMF_MAXSTR),save :: GSI_ExpId
   character(len=ESMF_MAXSTR),save :: GSI_bkg_fname_tmpl
   character(len=ESMF_MAXSTR),save :: GSI_ensbkg_fname_tmpl
   character(len=ESMF_MAXSTR),save :: GSI_ensana_fname_tmpl
   character(len=ESMF_MAXSTR),save :: GSI_ensprgA_fname_tmpl
   character(len=ESMF_MAXSTR),save :: GSI_ensprgB_fname_tmpl
   character(len=ESMF_MAXSTR),save :: GSI_fsens_fname_tmpl
   character(len=ESMF_MAXSTR),save :: GSI_ferrA_fname_tmpl
   character(len=ESMF_MAXSTR),save :: GSI_ferrB_fname_tmpl
   integer(i_kind)     :: GSI_ensread_blocksize

   logical             :: idco

   integer(i_kind) :: nmguess                                   ! number of extra meteol fields
   character(len=ESMF_MAXSTR),save,allocatable:: mguess_gsi(:)  ! names of extra met-fields per GSI
   character(len=ESMF_MAXSTR),save,allocatable:: mguess_usr(:)  ! names of extra met-fields per input

   integer(i_kind) :: ntgases                                   ! number of tracer gases (namelist)
   character(len=ESMF_MAXSTR),save,allocatable:: tgases_gsi(:)  ! names of tracer gases (internal names)
   character(len=ESMF_MAXSTR),save,allocatable:: tgases_usr(:)  ! names of tracer gases input names

   character(len=ESMF_MAXSTR),save,allocatable :: obstab(:,:)
   character(len=ESMF_MAXSTR),save,allocatable :: obscls(:,:)

   integer,parameter :: NFLDSLOT=2 ! in lobserver mode, internal guess grid size is controlled

   integer(i_kind),save :: GSI_HALOWIDTH = HW
   integer(i_kind),save :: PERTMOD_RUN_DT = HUGE(PERTMOD_RUN_DT)
   real(r_kind)   ,save :: GEOS_MU_SKIN
   type(ESMF_Grid),save :: GSI_AgcmPertGrid
   type(ESMF_Clock),save :: GSI_aClock
   type(ESMF_Time),save :: GSI_Time0
   type(ESMF_Time),save :: GSI_RefTime
   type(ESMF_Time),save :: GSI_FcsTime

! local utilities

   ! transpose IJ-JI
   interface GSI_GridCompSwapIJ_
     module procedure SwapIJ2r8_
     module procedure SwapIJ2i_
     module procedure SwapIJ3r8_
   end interface
   ! transpose JI-IJ
   interface GSI_GridCompSwapJI_
     module procedure SwapJI2r_
     module procedure SwapJI2i_
     module procedure SwapJI3r8_
   end interface
   ! swap vertical levels
   interface GSI_GridCompSwapV_
     module procedure SwapVr8_
     module procedure SwapVr4_
   end interface
   ! swap vertical levels
   interface GSI_GridCompFlipLons_
     module procedure GSI_GridCompFlipLons2_
   end interface

   interface GSI_clockShow
     module procedure myClock_show_
   end interface

   logical, save :: rc_obstable_initialized_=.false.

!---------------------------------------------------------------------------

   CONTAINS

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: GSI_GridCompSetServices -- Sets ESMF services for the GSI

! !INTERFACE:

   subroutine GSI_GridCompSetServices ( gc, rc )

!
! !USES:
!
  use gsi_4dcouplermod, only: gsi_4dCoupler_setServices
  use mpeu_util, only: tell
  implicit NONE

! !ARGUMENTS:

   type(ESMF_GridComp)               :: gc  ! gridded component
   integer            , intent(out ) :: rc

! !DESCRIPTION: This version uses the MAPL_GenericSetServices,
!        which sets the Run, Initialize, and Finalize services, 
!       as well as allocating our instance of a generic state and putting it in the 
!       gridded component (GC). 
!
!EOPI
!-------------------------------------------------------------------------

   integer                               :: STATUS
   logical                               :: IamRoot    
   character(len=*), parameter :: IAm='GSI_GridComp.SetServices'

! start
_ENTRY_(trim(Iam))

   IamRoot = MAPL_AM_I_ROOT()

   if(IamRoot) print *, Iam,": Start ",trim(Iam)

! Set the Initialize, Run, Finalize entry points

   call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE, Initialize, &
                                                       STATUS)
   VERIFY_(STATUS)

   call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN, Run,         &
                                                       STATUS)
   VERIFY_(STATUS)

   call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_FINALIZE, Finalize,  &
                                                       STATUS)
   VERIFY_(STATUS)

! Set Import/Export Coupling SPECS through the Generic Internal State

   call GSI_GridCompSetupSpecs (GC, rc=STATUS) 
   VERIFY_(STATUS)

! Set analysis timer

   call MAPL_TimerAdd(gc, name='Analysis time', rc=STATUS)
   VERIFY_(STATUS)

! Generic SetServices

   call MAPL_GenericSetServices ( gc, RC=STATUS)
   VERIFY_(STATUS)

! A non-standard serServices() call for non-standard pertmod component
   call gsi_4dCoupler_setServices (RC=STATUS)
   VERIFY_(STATUS)

   if(IamRoot.and.verbose) print *,Iam,": End ",trim(Iam)

   _EXIT_(trim(Iam))
   RETURN_(ESMF_SUCCESS)

   end subroutine GSI_GridCompSetServices

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: Initialize -- Initialize the GSI Gridded Comp

! !INTERFACE:

   subroutine Initialize ( gc, import, export, clock, rc )

!
! !USES:
!

   use gsi_4dcouplermod, only: gsi_4dCoupler_init_traj
   use mpeu_util, only: tell
   implicit NONE

! !ARGUMENTS:

   type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component
   type(ESMF_State),    intent(INOUT) :: import ! Import state
   type(ESMF_State),    intent(INOUT) :: export ! Export state
   type(ESMF_Clock),    intent(INOUT) :: clock  ! The clock
   integer, optional,   intent(  OUT) :: rc     ! Error code:

! !DESCRIPTION: This function initializes the GSI gridded component.
!       It  creates its ESMF grid, defines the import/export specs, and 
!       initializes other variables associated with the GSI legacy code.
!
! !REVISION HISTORY:
!
!   05Apr2007 CCruz/RT  Removed clock advance; add to AlarmSet interface
!   19Apr2007 Todling   Placed call to CompGetVertParms a bit sooner
!   12Nov2011 Todling   Add fix4hybrid
!
!EOPI
!-------------------------------------------------------------------------

   integer                          :: STATUS    ! error code STATUS
   type(MAPL_MetaComp), pointer     :: GENSTATE  ! GEOS Generic state
   type(ESMF_Config)                :: CF        ! configuration data
   type(ESMF_VM)                    :: VM        ! virtual machine
   type(ESMF_Grid)                  :: GSIGrid   ! this component's grid
   logical                          :: IamRoot
   integer                          :: GSIGridType
   type(ESMF_Time)                  :: AnaTime
   type(ESMF_Time)                  :: CurrTime
   type(ESMF_TimeInterval)          :: alarmInterval
   type (MAPL_VarSpec), pointer     :: import_spec(:)
   type (MAPL_VarSpec), pointer     :: export_spec(:)
   integer                          :: i, mype, atime
   integer                          :: yy, mm, dd, h, m
   character(len=ESMF_MAXSTR)       :: aname
   character(len=*), parameter :: IAm='GSI_GridComp.Initialize'

! start
_ENTRY_(trim(Iam))

   IamRoot =  MAPL_AM_I_ROOT()
   if(IamRoot.and.verbose) print *, Iam,": Start ",trim(Iam)

   if(IamRoot) call myClock_show_(clock,Iam,'entry-clock')

! Get vm, config objects

   call ESMF_GridCompGet( gc, config=CF, RC=STATUS )
   VERIFY_(STATUS)
   call ESMF_VMGetCurrent(vm=vm, rc=STATUS)
   VERIFY_(STATUS)
   call ESMF_VMGet(vm, petCount=npe, localPet=mype, &
        mpiCommunicator=mpi_comm_world, &
        rc=STATUS)
   VERIFY_(STATUS)

!                       ---------------------
!                       Initialize Legacy GSI
!                       ---------------------

! initialize internal GSI variables and read namelists
   
   if(IamRoot.and.verbose) print *, trim(Iam),': gsimain_Initialize'

   call GSI_GridCompDeterPEpartition_()

   call gsimain_Initialize

   call fix4hybrid_()

! REMARKS:
! 1) In the legacy code, domain decomposition, grid setup, etc, are
!    performed in routine gsisub(), which is now our run() method.
! 2) These are now done here during initialization as we need this
!    information to implement the couplers.
   
! domain decomposition

    call GSI_GridCompDeterSubdomain_()
     
! get GSI parameters that are usually obtained from a file header
! (vertical grid parameters ak, bk, sigma levels, etc.)

   call GSI_GridCompGetVertParms_(rc=STATUS)

! create GSI ESMF_GRID - using GSI domain decomposition
! Also set grid lons/lats.

   call GSI_GridCompGridCreate_(gsigrid,GSI_AgcmPertGrid)

! Make sure GSI G.C. gets this new grid

   call ESMF_GridCompSet(gc, grid=gsigrid, rc=STATUS)
   VERIFY_(STATUS)

! Determine how many global atmospheric and surface guess times are
! used (no. of first guess time levels, interval).

   call GSI_GridCompGetBKGTimes_()

! Allocate memory for GSI internal variables (first guess time levels)

   call GSI_GridCompAlloc_()

!                       -----------------------
!                       Initialize MAPL Generic
!                       -----------------------

! GEOS generic initialize
   
   call MAPL_GenericInitialize(gc, import, export, clock, rc=STATUS)
   VERIFY_(STATUS)

!                       --------------
!                       Set GSI alarms
!                       --------------

   call GSI_GridCompSetAlarms (GENSTATE, GC, cf, clock)

   GSI_aClock = ESMF_ClockCreate(clock, rc=STATUS)
   VERIFY_(STATUS)

   if(l4dvar) then
     call gsi_4dcoupler_init_traj(idmodel, rc=STATUS)
     VERIFY_(STATUS)
   endif

   if(IamRoot.and.verbose) print *, Iam,": End ",trim(Iam)

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( import, rc=STATUS )
    call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( export, rc=STATUS )
#endif

   if(IamRoot) call myClock_show_(clock,Iam,'exit-clock')

   _EXIT_(trim(Iam))
   RETURN_(STATUS)

   CONTAINS

   subroutine fix4hybrid_()
   use hybrid_ensemble_parameters, only: l_hyb_ens
   use hybrid_ensemble_parameters, only: nlat_ens
!  the following is a big hack, but the only way
!  dual resolution works ...
   if(l_hyb_ens) then  ! when doing hybrid GSI ...
      if(nlat_ens/=nlat) then ! and dual resolution ...
         if(mod(nlat_ens,2)/=0) nlat_ens=nlat_ens-1 ! make sure nlat is even.
         if(IamRoot) then
            print*, 'GSI_GridComp*fix4hybrid: WARNING WARNING WARNING WARNING WARNING   '
            print*, 'GSI_GridComp*fix4hybrid: WARNING, redefining nlat_ens to nlat_ens-1'
            print*, '==================================================================='
         endif
      endif
   endif
   end subroutine fix4hybrid_
!-------------------------------------------------------------------------
   subroutine GSI_GridCompDeterPEpartition_()
!-------------------------------------------------------------------------
   character(len=*), parameter :: IAm='GSI_GridCompDeterSubdomain_'
   integer :: default_nxpe, default_nype

   if(IamRoot.and.verbose) print *,trim(Iam),': determine GSI domain decomposition'

   call ESMF_ConfigGetAttribute( CF, nxpe, label ='NX:', rc = STATUS )
   VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute( CF, nype, label ='NY:', rc = STATUS )
   VERIFY_(STATUS)

!ALT: the next few lines are needed for concurrent (2 VMs) execution
   default_nxpe = nxpe
   default_nype = nype
   call ESMF_ConfigGetAttribute( CF, nxpe, label ='ANA_NX:', default=default_nxpe, rc = STATUS )
   VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute( CF, nype, label ='ANA_NY:', default=default_nype, rc = STATUS )
   VERIFY_(STATUS)

   end subroutine GSI_GridCompDeterPEpartition_

!-------------------------------------------------------------------------
   subroutine GSI_GridCompDeterSubdomain_()
!-------------------------------------------------------------------------
   character(len=*), parameter :: IAm='GSI_GridCompDeterSubdomain_'
   integer :: default_nxpe, default_nype

   if(IamRoot.and.verbose) print *,trim(Iam),': determine GSI domain decomposition'

   ! Note
   ! nlat: no. of latitudes
   ! nlon: no. of longitudes
   ! nsig: no. of levels
   ! npe:  no. of PEs
   ! These are defined in gsimain_Initialize
   
   ASSERT_(mod(nlon,8)==0)	! as required by GSI (e.g. 192,288,512,544)

   nlat_sfc = nlat
   nlon_sfc = nlon

   ! the following routines are called in gsisub:

   call create_grid_vars()
   
   end subroutine GSI_GridCompDeterSubdomain_

!-------------------------------------------------------------------------
   subroutine GSI_GridCompGridCreate_ ( grid, agcmPertGrid )

!-------------------------------------------------------------------------
! !REVISION HISTORY:
!
!  19Apr2007 Todling  Added ak and bk as attribute to file header
!  06Jan2009 RT/Guo   Redefine spectral coefficients on lat-lon grid
!  12Sep2011 Todling  Revisit calculation of grid for consistency w/ rest
!  15Apr2012 Todling  - Remove re-initialization of spectral routines
!                     - Extra pert-grid must be created all the time
!  05May2013 Todling  Bug fix: proper initialization of wgtlats
!
!-------------------------------------------------------------------------

   type (ESMF_Grid),  intent(inout)   :: grid    ! component grid
   type (ESMF_Grid),  intent(inout)   :: agcmPertGrid	! for 4dvar pertmod

! Local vars

   logical  :: flip_poles
   logical  :: flip_lons
   type(ESMF_DELayout) :: LAYOUT    ! DE layout
   integer, allocatable:: imxy(:), jmxy(:), lmxy(:)
   integer(kind(nlon)) :: i,j,k,i1,nzpe
   integer(i_kind)     :: ROOT
   real                :: lon0, lat0
   real                :: lon0d, lat0d
   real(kind(half))    :: pi, d2r
   real(kind(half))    :: dlon,dlat,pih
   character(len=*), parameter :: IAm='GSI_GridCompGridCreate'

   character(len=30) GSIGRIDNAME
   real, allocatable,  dimension(:)   :: ak5r4(:),bk5r4(:)
   type(LatLonGridFactory) :: ll_factory

! start

    rc  = 0
    pi  = 4.0_r_kind * atan ( 1.0_r_kind ) 
   d2r  = pi / 180._r_kind

! Horizontal grid : uniform (e.g. Equally spaced) or non-uniform (e.g. Gaussian)
! 0 = uniform, 1 = non-uniform

   call ESMF_ConfigGetAttribute( CF, GSIGridType, label ='GSIGridType:', rc = STATUS )
   VERIFY_(STATUS)

   if(IamRoot.and.verbose) print *,trim(Iam),': Will create GSI grid of TYPE  ',GsiGridType

! Query start longitude and latitude

   call ESMF_ConfigGetAttribute( CF, LON0, label ='ORIGIN_CENTER_LON:', rc = STATUS )
   VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute( CF, LAT0, label ='ORIGIN_CENTER_LAT:', rc = STATUS )
   VERIFY_(STATUS)

   lon0d= lon0; lat0d=lat0
   lon0 = lon0 * d2r
   lat0 = lat0 * d2r
   dlon=(pi+pi)/nlon	! in radians
   dlat=pi/(nlat-1)

   if(GsiGridType==0) then  ! equally spaced dgrid


! Set grid longitude array used by GSI.
      do i=1,nlon			! from 0 to 2pi
         rlons (i)=(i-one)*dlon
         coslon(i)=cos(rlons(i))
         sinlon(i)=sin(rlons(i))
      end do

! Set grid latitude array used by GSI.
      pih =half*pi
      do j=1,nlat			! from -pi/2 to +pi/2
         rlats(j)=(j-one)*dlat - pih
      end do

! wgtlats is used by spectral code. The values are used as divisor in the
! compact_diffs::inisph() routine.  Therefore, set to TINY instead of ZERO.
!     wgtlats(:)=TINY(wgtlats)
      wgtlats=zero
      do i=sp_a%jb,sp_a%je
         i1=i+1
         wgtlats(i1)=sp_a%wlat(i) !sp_a%clat(i)
         i1=nlat-i
         wgtlats(i1)=sp_a%wlat(i) !sp_a%clat(i)
      end do

! rbs2=1/cos^2(rlats)) is used in pcp.  polar points are set to zeroes.
      rbs2(1       )=zero
      rbs2(2:nlat-1)=cos(rlats(2:nlat-1))
      rbs2(2:nlat-1)=one/(rbs2(2:nlat-1)*rbs2(2:nlat-1))
      rbs2(  nlat  )=zero

   else                      ! Gaussian grid

! Set Gaussian grid lon/lat arrays used by GSI.
      call gengrid_vars  

  end if
 
! Re-Define South-West Corner of First Grid-Box
! ---------------------------------------------
!
!   NW---------NE
!   |           |
!   |     C     |
!   |           |
!   SW---------SE
!

! Create Grid. The deltas define whether is is uniform or non-uniform

   call MAPL_DefGridName (nlon,nlat,GSIGRIDNAME,MAPL_am_I_root())
   ll_factory = LatLonGridFactory(grid_name=GSIGRIDNAME, nx=nxpe, ny=nype, &
                im_world=nlon, jm_world=nlat, lm=nsig, pole='PC', dateline='DC', &
                rc=status)
   VERIFY_(status)
   grid = grid_manager%make_grid(ll_factory,rc=status)
   VERIFY_(STATUS)

   call ESMF_ConfigGetAttribute( CF, lon0d, label='AgcmPert_ORIGIN_CENTER_LON:', rc = STATUS )
  	if(STATUS/=ESMF_SUCCESS) lon0d=-180.
   call ESMF_ConfigGetAttribute( CF, lat0d, label='AgcmPert_ORIGIN_CENTER_LAT:', rc = STATUS )
     	if(STATUS/=ESMF_SUCCESS) lat0d=-90.

   ll_factory = LatLonGridFactory(grid_name=GSIGRIDNAME, nx=nxpe, ny=nype, &
                im_world=nlon, jm_world=nlat, lm=nsig, pole='PC', dateline='DC', &
                rc=status)
   VERIFY_(status)
   agcmPertGrid = grid_manager%make_grid(ll_factory,rc=status)
   VERIFY_(STATUS)

   if(IamRoot) then
     call tell(trim(Iam),'Creating agcmPertGrid with BegLon =',lon0d)
     call tell(trim(Iam),'                           BegLat =',lat0d)
   endif


! distribute grid
! ---------------

   call ESMF_VMBroadcast(vm, rlons, size(rlons), MAPL_Root, RC=STATUS); VERIFY_(STATUS)
   call ESMF_VMBroadcast(vm, rlats, size(rlats), MAPL_Root, RC=STATUS); VERIFY_(STATUS)

   ! Fields (bkg files) on a uniform grid have longitude range on [-pi,pi]
   ! and the latitudes are "correctly" oriented, so:
   if(GsiGridType==0) then 
      flip_lons  = .true.    ! flip to [0,2*pi]
      flip_poles = .false.
   end if   
   ! Fields on a Gaussian grid have longitude range on [0,2*pi]
   ! and the latitudes are oriented in reverse, so:
   if(GsiGridType==1) then 
      flip_lons  = .false.
      flip_poles = .true.
   end if   
   call ESMF_AttributeSet(grid, "FLIP_LONS", flip_lons, rc=STATUS)
   VERIFY_(STATUS)
   call ESMF_AttributeSet(grid, "FLIP_POLES", flip_poles, rc=STATUS)
   VERIFY_(STATUS)

!!! TODO: Set some grid attributes to perform these actions during coupling:
!!! GSI distributed fields are transposed from IJK to JIK
!   index_order = 2-1-3  ! e.g.
!   call ESMF_AttributeSet(grid, "GRID_INDEX_ORDER", index_order , rc= STATUS)
!   VERIFY_(STATUS)

   allocate(ak5r4(nsig),bk5r4(nsig))
   ak5r4=ak5
   call ESMF_AttributeSet(grid, name='ak', valuelist=ak5r4, rc=STATUS)
   VERIFY_(STATUS)
   bk5r4=bk5
   call ESMF_AttributeSet(grid, name='bk', valuelist=bk5r4, rc=STATUS)
   VERIFY_(STATUS)
!  ck5r4=ck5
!  call ESMF_AttributeSet(grid, name='ck', value=nsig+1, valuelist=ck5r4, rc=STATUS)
!  VERIFY_(STATUS)
   deallocate(ak5r4,bk5r4)

   end subroutine GSI_GridCompGridCreate_

!-------------------------------------------------------------------------
   subroutine GSI_GridCompGetVertParms_(rc)
   use m_set_eta, only: set_eta
   use mpimod, only: mpi_real8
!-------------------------------------------------------------------------
!
! !REVISION HISTORY:
!
!  21Mar2007 Todling  Added ck5
!
!-------------------------------------------------------------------------
   integer, optional,   intent(  OUT) :: rc     ! Error code:
   integer                            :: k
   real, allocatable,  dimension(:)   :: ak5r4(:),bk5r4(:)
   real(r_double), allocatable,  dimension(:):: ak(:),bk(:)
   real(r_double) :: ptop,pint
   integer ks
   integer, allocatable, dimension(:) :: mylevs(:)
   character(len=16)                  :: vgridlabl
   character(len=3)                   :: cnsig
   character(len=*), parameter :: IAm='GSI_GridCompGetVertParms'

! start

   if(IamRoot.and.verbose) print *,trim(Iam),': Get GSI g.c. parameters '

! Create the label to be searched for in the RC file based on nsig

   allocate ( ak (nsig+1), bk (nsig+1), stat=STATUS )
   if(IamRoot) then
      call set_eta ( nsig,ks,ptop,pint,ak,bk ) ! in GEOS orientation
   endif
   call mpi_bcast(ak, nsig+1,MPI_REAL8,MAPL_root,mpi_comm_world,rc)
   call mpi_bcast(bk, nsig+1,MPI_REAL8,MAPL_root,mpi_comm_world,rc)
   do i=1,nsig+1
      k=nsig-i+2              ! in gsi orientation
      ak5(k)=ak(i)*kPa_per_Pa ! convert to cb
      bk5(k)=bk(i)
      ck5(k)=zero
   enddo
   deallocate(ak,bk)

   if(IamRoot.and.verbose) then
      print *,trim(IAm),' - lev, ak, bk - '
      do i=1,nsig+1
         write(*,'(1x,i3,2f16.6)') i,ak5(i),bk5(i)
      end do
   end if

! Since sigma = P/Ps => P = sigma*Ps
! but P = ak + bk*Ps
! hence for Gaussian grid type set ak=0 => sigma=bk

   idvc5 = 2    ! set to sigma-pressure
   if(GsiGridType==1) then
      idvc5 = 1 ! reset to sigma-only
      if(nsig/=64) then
         print *,trim(IAm),' GSIsa_Gaussian only supports 64 levels; nsig = ',nsig
         RETURN_(ESMF_FAILURE)
      end if
   end if

   end subroutine GSI_GridCompGetVertParms_

!-------------------------------------------------------------------------
   subroutine GSI_GridCompGetBKGTimes_()
!-------------------------------------------------------------------------
   use mpeu_util, only: tell
   implicit none
   
! local variables

   real(r_kind)     :: bkgbits
   integer(i_kind)  :: i
   integer(i_kind)  :: JOB_SGMT(2)
   integer(i_kind)  :: BKGfreq
   integer(i_kind)  :: ANAfreq
   integer(i_kind)  :: RUN_DT
   character(len=*),parameter :: IAm='GSI_GridCompGetBKGTimes'

! Begin...

   if(IamRoot.and.verbose) print *,trim(Iam),': Get GSI background times'

   ! proceed as in GSI's read_files.f90 

   call ESMF_ConfigGetAttribute( CF, GSI_ExpId, label ='expid:', rc = STATUS )
   VERIFY_(STATUS)
   
   call ESMF_ConfigGetAttribute(CF, JOB_SGMT, count=2, label='JOB_SGMT:', rc=STATUS)
        VERIFY_(STATUS)

   ANAwndw_hr  = JOB_SGMT(2)/10000 + JOB_SGMT(1)*24

   CALL ESMF_ConfigGetAttribute(CF, BKGfreq, label = 'BKG_FREQUENCY:', rc=STATUS)
        VERIFY_(STATUS)

   BKGfreq_hr = BKGfreq/10000
   BKGfreq_mn = mod(BKGfreq,10000)/100
   BKGfreq_sc = mod(BKGfreq,100)
   BKGfreq_sc = BKGfreq_hr*3600 + BKGFreq_mn*60 + BKGfreq_sc
   BKGfrac_hr = BKGfreq_hr + real(BKGfreq_mn)/60
 
   call ESMF_ConfigGetAttribute(CF, RUN_DT, label='RUN_DT:', rc=STATUS)
        VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute(CF, PERTMOD_RUN_DT, label='PERTMOD_RUN_DT:', default=RUN_DT, rc=STATUS)
        VERIFY_(STATUS)  
   ! if GSI_DT specified have it overwrite run_dt, but leave pertmod_run_dt untouched
   call ESMF_ConfigGetAttribute(CF, RUN_DT, label='GSI_DT:', default=RUN_DT, rc=STATUS)
        VERIFY_(STATUS)

   nbkgfreq  = BKGfreq_sc / RUN_DT    ! bkg per dt

   CALL ESMF_ConfigGetAttribute(CF, ANAfreq, label = 'ANA_FREQUENCY:', rc=STATUS)
        VERIFY_(STATUS)

   ANAfreq_hr = ANAfreq / 10000
   ANAfreq_mn = mod(ANAfreq,10000)/100
   ANAfreq_sc = mod(ANAfreq,100)
   ANAfreq_sc = ANAfreq_hr*3600 + ANAFreq_mn*60 + ANAfreq_sc

   if ( mod(BKGfreq_sc,6*3600)==0) then
        nfldsig_all = nhr_assimilation*3600 / BKGfreq_sc
   else
        nfldsig_all = nhr_assimilation*3600 / BKGfreq_sc + 1
   endif
   nfldsfc_all = nfldsig_all   ! In GEOS-5 upper-air bkg and surface bkg come equal numbers
   nfldnst_all = nfldsig_all   ! make it simple for now

	! to use gsimod in a clock-driving mode, only upto two time slots are
	! needed for guess_grid

   nfldsig = nfldsig_all
   nfldsfc = nfldsfc_all
   nfldnst = nfldnst_all
   if(lobserver) then
     nfldsig = min(NFLDSLOT,nfldsig_all)
     nfldsfc = min(NFLDSLOT,nfldsfc_all)
     nfldnst = min(NFLDSLOT,nfldnst_all)
   endif

   extrap_intime = .not.lobserver .or. nfldsig_all==1

   nfldsig_now = 0	! no data has filled into guess_grids
   nfldsfc_now = 0

   ! these are internal GSI global variables
   ! and are deallocated in Finalize
   allocate(hrdifsig    (1:nfldsig    ), &
   	    hrdifsfc    (1:nfldsfc    ), &
   	    hrdifnst    (1:nfldnst    ), &
	    hrdifsig_all(1:nfldsig_all), &
	    hrdifnst_all(1:nfldnst_all), &
   	    hrdifsfc_all(1:nfldsfc_all), stat=STATUS)
   VERIFY_(STATUS)
   allocate(ifilesig(1:nfldsig_all), ifilesfc(1:nfldsfc_all), ifilenst(1:nfldnst_all), stat=STATUS)
   VERIFY_(STATUS)

   ifilesig(1:nfldsig_all)=6    ! it doesn't matter what goes in ifilesig
   ifilesfc(1:nfldsig_all)=6    ! it doesn't matter what goes in ifilesfc
   ifilenst(1:nfldsig_all)=6    ! it doesn't matter what goes in ifilesfc

#ifdef VERBOSE
   call tell(Iam,'nfldsig_all=',nfldsig_all)
   call tell(Iam,'nfldsfc_all=',nfldsfc_all)
   call tell(Iam,'BKGfrac_hr =',BKGfrac_hr )
#endif
   bkgbits = 0.0
   do i = 1, nfldsig_all
      hrdifsig_all(i)  = bkgbits
      bkgbits      = bkgbits + bkgfrac_hr
!_RT#ifdef VERBOSE
      if (IamRoot) then
         call tell(Iam,'             i =',i)
         call tell(Iam,'hrdifsig_all(i)=',hrdifsig_all(i))
      endif
!_RT#endif
   enddo

   bkgbits = 0.0
   do i = 1, nfldsfc_all
      hrdifsfc_all(i)  = bkgbits
      bkgbits      = bkgbits + bkgfrac_hr
#ifdef VERBOSE
      call tell(Iam,'             i =',i)
      call tell(Iam,'hrdifsfc_all(i)=',hrdifsfc_all(i))
#endif
   enddo
   hrdifnst_all =  hrdifsfc_all ! make it simple

   call ESMF_ConfigGetAttribute(CF, GSI_ensread_blocksize, label='ensemble_read_blocksize:', &
                                default=4,rc=STATUS)

!  Read in filename template for background files
!  ----------------------------------------------
   call ESMF_ConfigGetAttribute(CF, GSI_bkg_fname_tmpl, label='upper-air_bkg_filename:', &
                                default='%s.bkg.eta.%y4%m2%d2_%h2z.nc4',rc=STATUS)
        VERIFY_(STATUS)

!  Read in filename template for ensemble of backgrounds
!  -----------------------------------------------------
   call ESMF_ConfigGetAttribute(CF, GSI_ensbkg_fname_tmpl, label='ensemble_upabkg_filename:', &
                                default='%s.bkg.eta.%y4%m2%d2_%h2z.nc4',rc=STATUS)
        VERIFY_(STATUS)

!  Read in filename template for ensemble of analysis
!  --------------------------------------------------
   call ESMF_ConfigGetAttribute(CF, GSI_ensana_fname_tmpl, label='ensemble_upaana_filename:', &
                                default='%s.ana.eta.%y4%m2%d2_%h2z.nc4',rc=STATUS)
        VERIFY_(STATUS)
!  Read in filename template for ensemble of forecasts
!  ---------------------------------------------------
   call ESMF_ConfigGetAttribute(CF, GSI_ensprgA_fname_tmpl, label='ensemble_upaprga_filename:', &
                                default='%s.proga.eta.%y4%m2%d2_%h2z.nc4',rc=STATUS)
        VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute(CF, GSI_ensprgB_fname_tmpl, label='ensemble_upaprgb_filename:', &
                                default='%s.progb.eta.%y4%m2%d2_%h2z.nc4',rc=STATUS)
        VERIFY_(STATUS)
!  Define template filename of forecast sensitivity vector
!  -------------------------------------------------------
   call ESMF_ConfigGetAttribute(CF, GSI_fsens_fname_tmpl, label='forecast_sensitivity_filename:', &
                                default='fsens.eta.nc4',rc=STATUS)
        VERIFY_(STATUS)
!  Define template filename of forecast error vectors
!  -------------------------------------------------------
   call ESMF_ConfigGetAttribute(CF, GSI_ferrA_fname_tmpl, label='forecast_errorA_filename:', &
                                default='ferrA.eta.nc4',rc=STATUS)
        VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute(CF, GSI_ferrB_fname_tmpl, label='forecast_errorB_filename:', &
                                default='ferrB.eta.nc4',rc=STATUS)
        VERIFY_(STATUS)

!  Read in the temperature profile exponent: mu_skin used for skin SST analysis
!  SA: PERHAPS THIS resource parameter SHOULD NOT BE in GSI_GridCompGetBKGTimes_. 
!  Ricardo, plz suggest suitable place.
!  ----------------------------------------------------------------------------
   CALL ESMF_ConfigGetAttribute(CF, GEOS_MU_SKIN, label = 'mu_skin:', default=0.2_r_kind, rc=STATUS)
        VERIFY_(STATUS)
   if(IamRoot) print *,trim(Iam),': Set MU_SKIN= ', GEOS_MU_SKIN
  
   end subroutine GSI_GridCompGetBKGTimes_

!-------------------------------------------------------------------------
   subroutine GSI_GridCompAlloc_()
!-------------------------------------------------------------------------

! allocation of GSI internal state variables

   character(len=*), parameter :: IAm='GSI_GridCompAlloc'
   integer(i_kind) ii,jj,kk,nn

   allocate( isli     (lat2,lon2,nfldsfc),&
             fact10   (lat2,lon2,nfldsfc),&
             sfct     (lat2,lon2,nfldsfc),&
             dsfct    (lat2,lon2,nfldsfc),&
             sno      (lat2,lon2,nfldsfc),&
             veg_type (lat2,lon2,nfldsfc),&
             veg_frac (lat2,lon2,nfldsfc),&
             soil_type(lat2,lon2,nfldsfc),&
             soil_temp(lat2,lon2,nfldsfc),&
             soil_moi (lat2,lon2,nfldsfc),&
             sfc_rough(lat2,lon2,nfldsfc),&
             stat=STATUS)
   VERIFY_(STATUS)
   do nn=1,nfldsfc
      do jj=1,lon2
         do ii=1,lat2
            isli(ii,jj,nn)=0
            fact10(ii,jj,nn)=zero
            sfct(ii,jj,nn)=zero
            dsfct(ii,jj,nn)=zero
            sno(ii,jj,nn)=zero
            veg_type(ii,jj,nn)=zero
            veg_frac(ii,jj,nn)=zero
            soil_type(ii,jj,nn)=zero
            soil_temp(ii,jj,nn)=zero
            soil_moi(ii,jj,nn)=zero
         end do
      end do
   end do

!  When proper connection to ESMF is complete,
!  the following will not be needed here
!  ------------------------------------------
   if (nmguess>0) then
       call create_metguess_grids(mype,status)
       VERIFY_(STATUS)
   endif

   if (ntgases>0) then
       call create_chemges_grids(mype,status)
       VERIFY_(STATUS)
   endif

   RETURN_(ESMF_SUCCESS)

   end subroutine GSI_GridCompAlloc_

   end subroutine Initialize

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: Run -- Run the GSI Gridded Component

! !INTERFACE:

   subroutine Run ( gc, import, export, clock, rc )

!
! !USES:
!
   use mpeu_util, only: tell
   use ncepgfs_io, only: read_gfs_chem
   implicit NONE

! !ARGUMENTS:

   type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component
   type(ESMF_State),    intent(INOUT) :: import ! Import state
   type(ESMF_State),    intent(INOUT) :: export ! Export state
   type(ESMF_Clock),    intent(INOUT) :: clock  ! The clock
   integer, optional,   intent(  OUT) :: rc     ! Error code:

! !DESCRIPTION: set up and run the GSI analysis.
!
!  05Apr2007 Cruz    Handle for many background fields  
!  20Apr2007 Todling Added phis to export state
!  23Apr2007 Todling Added export of delp
!  07Apr2011 Guo     Nullify chem pointers
!  07May2011 Todling/Merkova - add cloud-fraction for radiation to import
!  11Jun2013 Todling Allow overwrite of trace gases in background
!
!EOPI
!-------------------------------------------------------------------------

! local variables

   integer                           :: i, j, k, L
   integer                           :: im,jm,km
   integer                           :: mype
   integer                           :: GSIGridType
   integer                           :: STATUS    ! error code STATUS
   type(MAPL_MetaComp), pointer      :: GENSTATE  ! GEOS Generic state
   type(ESMF_Config)                 :: CF        ! coneiguration data
   type(ESMF_VM)                     :: VM        ! virtual machine
   type(ESMF_Grid)                   :: GSIGrid   ! this component's grid
   type(ESMF_Alarm)                  :: GSIALARM
   logical                           :: do_analysis
   logical                           :: do_observer
   logical                           :: get_background
   logical                           :: IamRoot
   ! import state upper air pointers
   real(r_single),dimension(:,:  ), pointer :: hsp  ! terrain
   real(r_single),dimension(:,:  ), pointer :: psp  ! surf. pressure
   real(r_single),dimension(:,:,:), pointer :: up   ! u wind
   real(r_single),dimension(:,:,:), pointer :: vp   ! v wind
   real(r_single),dimension(:,:,:), pointer :: tp   ! virtual temp.
   real(r_single),dimension(:,:,:), pointer :: qp   ! spec. hum.
   real(r_single),dimension(:,:,:), pointer :: ozp  ! ozone
   real(r_single),dimension(:,:,:), pointer :: qimr ! cloud ice    mixing ratio
   real(r_single),dimension(:,:,:), pointer :: qlmr ! cloud liquid mixing ratio
   real(r_single),dimension(:,:,:), pointer :: qrmr ! rain mixing ratio
   real(r_single),dimension(:,:,:), pointer :: qsmr ! snow mixing ratio
   real(r_single),dimension(:,:,:), pointer :: clfr ! cloud fraction for radiation
   ! import chem tracers ... preliminary (dummy implementation)
   real(r_single),dimension(:,:,:), pointer :: cop   =>NULL() ! carbone monoxide
   real(r_single),dimension(:,:,:), pointer :: co2p  =>NULL() ! carbone dioxide
   ! import aerosol optical depth
   real(r_single),dimension(:,:),   pointer :: aodp =>NULL()  ! aerosol optical depth
   ! import aerosols: dust
   real(r_single),dimension(:,:,:), pointer :: du001p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: du002p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: du003p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: du004p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: du005p =>NULL()  ! 
   ! import aerosols: sea-salt
   real(r_single),dimension(:,:,:), pointer :: ss001p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: ss002p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: ss003p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: ss004p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: ss005p =>NULL()  ! 
   ! import aerosols: sulfates
   real(r_single),dimension(:,:,:), pointer :: dmsp =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: so2p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: so4p =>NULL()  ! 
   real(r_single),dimension(:,:,:), pointer :: msap =>NULL()  ! 
   ! import carbonaceous
   real(r_single),dimension(:,:,:), pointer :: bcphobicp =>NULL()  ! dry black carbon
   real(r_single),dimension(:,:,:), pointer :: bcphilicp =>NULL()  ! wet black carbon
   real(r_single),dimension(:,:,:), pointer :: ocphobicp =>NULL()  ! dry organic carbon
   real(r_single),dimension(:,:,:), pointer :: ocphilicp =>NULL()  ! wet organic carbon
   ! import state surface pointers
   real(r_single),dimension(:,:  ), pointer :: f10p ! 10m winf factors
   real(r_single),dimension(:,:  ), pointer :: tskp ! skin Temp.
   real(r_single),dimension(:,:  ), pointer :: snop ! snow depth
   real(r_single),dimension(:,:  ), pointer :: sotp ! soil Temp.
   real(r_single),dimension(:,:  ), pointer :: soqp ! soil moist. 
   real(r_single),dimension(:,:  ), pointer :: frland    ! land fraction
   real(r_single),dimension(:,:  ), pointer :: frlandice ! land-ice fraction
   real(r_single),dimension(:,:  ), pointer :: frlake    ! lake fraction
   real(r_single),dimension(:,:  ), pointer :: frocean   ! ocean fraction
   real(r_single),dimension(:,:  ), pointer :: frseaice  ! sea-ice fraction
   real(r_single),dimension(:,:  ), pointer :: vtyp ! veg. type
   real(r_single),dimension(:,:  ), pointer :: styp ! soil type
   real(r_single),dimension(:,:  ), pointer :: vfrp ! veg. frac.
   real(r_single),dimension(:,:  ), pointer :: sz0p ! surf roughness
   ! import state skin-layer ocean variables
   real(r_single),dimension(:,:  ), pointer :: z_cp    =>NULL()  ! depth of cool layer
   real(r_single),dimension(:,:  ), pointer :: z_wp    =>NULL()  ! depth of warm layer
   real(r_single),dimension(:,:  ), pointer :: dt_coolp=>NULL()  ! temperature drop from skin to base of cool layer
   real(r_single),dimension(:,:  ), pointer :: tdelp   =>NULL()  ! temperature at top of warm layer
   real(r_single),dimension(:,:  ), pointer :: trefp   =>NULL()  ! foundation sea surface temperature- with no diurnal variation 
   ! u10m and v10m used to calculate 10m wind factors
   real(r_single),dimension(:,:  ), pointer :: u10p, v10p
   ! export state pointers - tendencies
   real(r_single),dimension(:,:  ), pointer :: df10  ! factor 10m
   real(r_single),dimension(:,:  ), pointer :: dsli  ! land/water/ice mask
   real(r_single),dimension(:,:  ), pointer :: dts   ! skin/surface temperature
   real(r_single),dimension(:,:  ), pointer :: dhs   ! terrain
   real(r_single),dimension(:,:  ), pointer :: dps   ! surf pressure
   real(r_single),dimension(:,:,:), pointer :: ddp   ! del pressure
   real(r_single),dimension(:,:,:), pointer :: du    ! u wind
   real(r_single),dimension(:,:,:), pointer :: dv    ! v wind
   real(r_single),dimension(:,:,:), pointer :: dt    ! virtual Temp.
   real(r_single),dimension(:,:,:), pointer :: dq    ! spec. hum.
   real(r_single),dimension(:,:,:), pointer :: doz   ! ozone
   real(r_single),dimension(:,:,:), pointer :: dqimr ! cloud ice    mixing ratio
   real(r_single),dimension(:,:,:), pointer :: dqlmr ! cloud liquid mixing ratio
   real(r_single),dimension(:,:,:), pointer :: dqrmr ! rain mixing ratio
   real(r_single),dimension(:,:,:), pointer :: dqsmr ! snow mixing ratio
   real(r_single),dimension(:,:  ), pointer :: dfrland    ! land fraction
   real(r_single),dimension(:,:  ), pointer :: dfrlandice ! land-ice fraction
   real(r_single),dimension(:,:  ), pointer :: dfrlake    ! lake fraction
   real(r_single),dimension(:,:  ), pointer :: dfrocean   ! ocean fraction
   real(r_single),dimension(:,:  ), pointer :: dfrseaice  ! sea-ice fraction
   ! export chem tracers ... preliminary (dummy implementation)
   real(r_single),dimension(:,:,:), pointer :: dcop  ! carbone monoxide
   real(r_single),dimension(:,:,:), pointer :: dco2p ! carbone dioxide
   ! export AOD
   real(r_single),dimension(:,:),   pointer :: daodp ! AOD
   !
   character(len=ESMF_MAXSTR)        :: aname
   character(len=ESMF_MAXSTR)        :: opt
   character(len=*), parameter :: IAm='Gsi_GridComp.Run'
   integer :: ier
   integer :: atime

! start

   IamRoot = MAPL_AM_I_ROOT() 
   if(IamRoot.and.verbose) print *, Iam,": Start ",trim(Iam)

   call ESMF_VMGetCurrent(vm=vm, rc=STATUS)
   VERIFY_(STATUS)
   call ESMF_VMGet(vm, petCount=npe, localPet=mype, rc=STATUS)
   VERIFY_(STATUS)
   call ESMF_GridCompGet( gc, config=CF, RC=STATUS )
   VERIFY_(STATUS)
   call ESMF_GridCompGet( gc, grid=GSIgrid, RC=STATUS )
   VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute( CF, GSIGridType, label ='GSIGridType:', rc = STATUS )
   VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute( cf, opt, label ='ANALYSIS_INCREMENT:', rc = STATUS )
   VERIFY_(STATUS)
   opt = adjustl(opt)
   if ( scan(opt(1:1),"Yy") /= 0 ) then
        if(iwrtinc<0) then
            if(IamRoot) print *,trim(Iam)//': overwrite GSI namelist setting, iwrtinc = ', iwrtinc
        endif
        iwrtinc = 1 ! other than being >0, this index does not control output in case of GMAO interface
   endif



! Get ALARMs
!-----------

   call MAPL_GetObjectFromGC ( GC, GENSTATE, RC=STATUS); VERIFY_(STATUS)

   nthTimeIndex = nthTimeIndex + 1
   if(IamRoot) print *,trim(Iam)//': time index = ',nthTimeIndex, ' bkgfreq = ', nbkgfreq

!  Determine if it's time do gather background
!  -------------------------------------------
!  if ( mod(nthTimeIndex-1,nbkgfreq) .ne. 0 ) then
!      if(IamRoot) print *,trim(Iam),': no time for bkg'
!      RETURN_(ESMF_SUCCESS)
!  endif
!  L = (nthTimeIndex-1)/nbkgfreq + 1
   L =  nthTimeIndex

!  Determine if it's time to actually do the analysis
!  --------------------------------------------------
   do_analysis = .false.
   call MAPL_StateAlarmGet(GENSTATE, GSIALARM, NAME='last-bkg', RC=STATUS); VERIFY_(STATUS)
   do_analysis = ESMF_AlarmIsRinging(GSIALARM, RC=STATUS); VERIFY_(STATUS)
   do_analysis = do_analysis .and. (.not.lobserver)

   do_observer = L>=nfldsig  .and. lobserver

#ifdef VERBOSE
   if(IamRoot) then
     call tell(Iam,'lobserver =',lobserver)
     call tell(Iam,'do_analysis =',do_analysis)
     call tell(Iam,'do_observer =',do_observer)
     call tell(Iam,'nthTimeIndex =',nthTimeIndex)
     call tell(Iam,"L =",L)
     call tell(Iam,"nfldsig_all",nfldsig_all)
     call tell(Iam,"nfldsig_now",nfldsig_now)
     call tell(Iam,"nfldsig",nfldsig)
     call tell(Iam,'init_pass(L==nfldsig) =',(L==nfldsig))
     call tell(Iam,'last_pass(L==nfldsig_all) =',(L==nfldsig_all))
   endif
#endif

!  Set alarm
!  ---------

   call GSI_GridCompSetAnaTime_()
#ifdef VERBOSE
   call tell(Iam,"returned from GSI_GridCompSetAnaTime_()")
#endif

   call GSI_GridCompGetPointers_()
   call GSI_GridCompCopyImportDyn2Internal_(L)
   call GSI_GridCompComputeVorDiv_(L)
   call GSI_GridCompCopyImportSfc2Internal_(L)
   call GSI_GridCompGetNCEPsfcFromFile_(L)

!  Set observations input
!  ----------------------
   if(L==1)then
     call mpi_barrier(mpi_comm_world,ier)
     if(IamRoot) call GSI_GridCompSetObsNames_(L)
     call mpi_barrier(mpi_comm_world,ier)
   endif
#ifdef VERBOSE
   call tell(Iam,"returned from GSI_GridCompSetObsNames_(L) at L =",L)
#endif

!  Run observer or analysis
!  ------------
   if(lobserver) then
     if(.not. do_observer) then
#ifdef VERBOSE
       if(IamRoot) call tell(Iam,'skip back to AANA with do_observer =',do_observer)
#endif
       RETURN_(ESMF_SUCCESS)
     endif

     call gsimain_run(  init_pass=(L==nfldsig), &
     			last_pass=(L==nfldsig_all)) 	! if in "observer" mode
#ifdef VERBOSE
     call tell(Iam,"returned from gsimain_run()")
#endif

   else

     if(.not. do_analysis) then
#ifdef VERBOSE
       if(IamRoot) call tell(Iam,'skip back to AANA with do_analysis =',do_analysis)
#endif
       RETURN_(ESMF_SUCCESS)
     endif

     call gsimain_run( init_pass=.true.,last_pass=.true.)
#ifdef VERBOSE
     call tell(Iam,"returned from gsimain_run()")
#endif

!  Copy analysis to export state
!  -----------------------------
     call GSI_GridCompCopyInternal2Export_(ntguessig)
#ifdef VERBOSE
     call tell(Iam,"returned from GSI_GridCompCopyInternal2Export_(), ntguessig=",ntguessig)
#endif

   endif

!  Done
!  ----
#ifdef VERBOSE
   if(IamRoot) call tell(Iam,'Exiting ...')
#endif

   RETURN_(ESMF_SUCCESS)

   CONTAINS

!-------------------------------------------------------------------------
   subroutine GSI_GridCompGetAlarms_(GENSTATE,nfldsig_all)
!-------------------------------------------------------------------------
   type(MAPL_MetaComp),intent(inout)   :: GENSTATE  ! GEOS Generic state
   integer,intent(in)               :: nfldsig_all

   type(ESMF_Alarm), pointer        :: ALARM(:)
   type(ESMF_Time)                  :: CurrTime
   type(ESMF_Time)                  :: AlaTime
   integer                          :: atime
   integer                          :: yy, mm, dd, hh, mn, sec, n
   character*2                      :: idxtime
   character(len=*), parameter :: IAm='GSI_GridCompGetAlarms'
                                                                                                                                
! start
                                                                                                                                
   call ESMF_ClockGet (clock, currTime=currTime, rc=STATUS)
   VERIFY_(STATUS)
   call ESMF_TimeGet(currTime, yy=YY, mm=MM, dd=DD, h=HH, m=MN, s=SEC, rc=status)
   VERIFY_(STATUS)
   if(IamRoot) then
      write(6,'(a,1x,i4.4,5(a,i2.2))') trim(Iam)//" The GSIgc RUN current TIME is ", &
                                      YY, "/", MM, "/", DD, " ", HH, ":", MN, ":", SEC
   end if
                                                                                                                                
!!$ Is this part of code doing anything useful?  It seems ALARM is first allocated
!!$ and (shallow) assigned to 'some-bkg' and 'lat-bkg' alarms of GENSTATE.  Then
!!$ throughed away after geting one of alarm for its ESMF_Time value, known only
!!$ locally, as AlaTime.
!!$ 

!  get alarm at ana-t

   allocate(ALARM(nfldsig_all), stat=STATUS); VERIFY_(STATUS)
                                                                                                                                
   do atime=1,nfldsig_all-1
     write(idxtime,'(i2)',iostat=STATUS) atime
     aname = "some-bkg " // idxtime
     VERIFY_(STATUS)
     call MAPL_StateAlarmGet(GENSTATE, ALARM(atime), NAME=aname, RC=STATUS)
     VERIFY_(STATUS)
   end do

   call MAPL_StateAlarmGet(GENSTATE, ALARM(nfldsig_all), NAME='last-bkg', RC=STATUS)
   VERIFY_(STATUS)
   call ESMF_AlarmGet(ALARM(nfldsig_all), ringTime=AlaTime, rc=status)
   VERIFY_(STATUS)
                                                                                                                                
   deallocate(ALARM, stat=STATUS); VERIFY_(STATUS)
   RETURN_(ESMF_SUCCESS)


   end subroutine GSI_GridCompGetAlarms_

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: GSI_GridCompGetPointers_  -- get pointers

! !INTERFACE:

   subroutine GSI_GridCompGetPointers_()

   implicit none

! !DESCRIPTION: get pointers to fields in import and export states
!
! !REVISION HISTORY:
!
!  19Mar2007 Todling  GSI no longer handles log(ps)
!  20Apr2007 Todling  Added phis(dhs) to export state
!  22Apr2007 Todling  Properly naming u/v/tv fields in file
!  07May2011 Todling/Merkova - add cloud-fraction for radiation to import
!
!EOPI
!-------------------------------------------------------------------------
!   type(ESMF_Grid)                   :: GEOSGrid 
!   type(ESMF_Field)                  :: GEOSFld

   character(len=*), parameter :: IAm='GSI_GridCompGetPointers_'

   integer(i_kind) :: nt,i4crtm
   character(len=ESMF_MAXSTR) :: cvar
    

! imports
! -------

   call ESMFL_StateGetPointerToData(import, tskp,     'ts',   rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(import, snop,  'SNOWDP',  rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(import, soqp,  'GWETTOP', rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(import, sotp,  'TSOIL1',  rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(import, frland,   'frland',     rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(import, frlandice,'frlandice',  rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(import, frlake,   'frlake',     rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(import, frocean,  'frocean',    rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(import, frseaice,'frseaice',    rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(import, sz0p,     'Z0M',  rc=STATUS)
   VERIFY_(STATUS)
   if(GsiGridType==0) then
      call ESMFL_StateGetPointerToData(import, u10p,  'U10M', rc=STATUS)
      VERIFY_(STATUS)
      call ESMFL_StateGetPointerToData(import, v10p,  'V10M', rc=STATUS)
      VERIFY_(STATUS)
   else ! f10m is provided in a file and stored in V10M
      call ESMFL_StateGetPointerToData(import, f10p,  'V10M', rc=STATUS)
      VERIFY_(STATUS)
   end if

! Will bootstrap these variables ...
! ----------------------------------
!  call ESMFL_StateGetPointerToData(import, vfrp, 'NCEP_VEGFRAC', rc=STATUS)
!  if(status==0) ncbootstrap(1)=.false.
!  call ESMFL_StateGetPointerToData(import, vtyp, 'NCEP_VEGTYPE', rc=STATUS)
!  if(status==0) ncbootstrap(2)=.false.
!  call ESMFL_StateGetPointerToData(import, styp, 'NCEP_SOITYPE', rc=STATUS)
!  if(status==0) ncbootstrap(3)=.false.
!  status=0

! Extra Meteorological fields
! ---------------------------
   do nt=1,nmguess
      cvar = trim(mguess_usr(nt))
      select case (uppercase(cvar))
!        case ( 'qctot' )
!           call ESMFL_StateGetPointerToData(import, qcmr, trim(cvar), rc=STATUS)
!           VERIFY_(STATUS)
!           if(mype==0) write(6,*) trim(Iam), ': this variable is local: ', trim(cvar)
         case ( 'PS' )
            call ESMFL_StateGetPointerToData(import, psp,trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ( 'PHIS' )
            call ESMFL_StateGetPointerToData(import, hsp,trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ( 'U' )
            call ESMFL_StateGetPointerToData(import, up, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ( 'V' )
            call ESMFL_StateGetPointerToData(import, vp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ( 'TV' )
            call ESMFL_StateGetPointerToData(import, tp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ( 'SPHU' )
            call ESMFL_StateGetPointerToData(import, qp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ( 'OZONE' )
            call ESMFL_StateGetPointerToData(import, ozp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ( 'QITOT' )
            call ESMFL_StateGetPointerToData(import, qimr, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('QLTOT')
            call ESMFL_StateGetPointerToData(import, qlmr, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('QRTOT')
            call ESMFL_StateGetPointerToData(import, qrmr, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('QSTOT')
            call ESMFL_StateGetPointerToData(import, qsmr, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('CLOUD')
            call ESMFL_StateGetPointerToData(import, clfr, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('DCOOL')
            call ESMFL_StateGetPointerToData(import, z_cp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('DWARM')
            call ESMFL_StateGetPointerToData(import, z_wp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('TDROP')
            call ESMFL_StateGetPointerToData(import, dt_coolp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('TDEL')
            call ESMFL_StateGetPointerToData(import, tdelp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('TS_FOUND')
            call ESMFL_StateGetPointerToData(import, trefp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case default
            if(mype==0) write(6,*) trim(Iam), ': ', trim(cvar), ' not in import state ...'
!           status = 1
!           VERIFY_(STATUS)
      end select
   enddo

! Chemistry tracer imports
! When proper connection w/ Tracer Bundel is made
! the following won't be necessary
! -----------------------------------------------
   do nt=1,ntgases
      cvar = trim(tgases_usr(nt))
      select case (uppercase(cvar))
         case ( 'AOD' )
            call ESMFL_StateGetPointerToData(import, aodp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ( 'CO' )
            call ESMFL_StateGetPointerToData(import, cop, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('CO2')
            call ESMFL_StateGetPointerToData(import, co2p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('DU001')
            call ESMFL_StateGetPointerToData(import, du001p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('DU002')
            call ESMFL_StateGetPointerToData(import, du002p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('DU003')
            call ESMFL_StateGetPointerToData(import, du003p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('DU004')
            call ESMFL_StateGetPointerToData(import, du004p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('DU005')
            call ESMFL_StateGetPointerToData(import, du005p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('SS001')
            call ESMFL_StateGetPointerToData(import, ss001p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('SS002')
            call ESMFL_StateGetPointerToData(import, ss002p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('SS003')
            call ESMFL_StateGetPointerToData(import, ss003p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('SS004')
            call ESMFL_StateGetPointerToData(import, ss004p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('SS005')
            call ESMFL_StateGetPointerToData(import, ss005p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('DMS')
            call ESMFL_StateGetPointerToData(import, dmsp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('SO2')
            call ESMFL_StateGetPointerToData(import, so2p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('SO4')
            call ESMFL_StateGetPointerToData(import, so4p, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('MSA')
            call ESMFL_StateGetPointerToData(import, msap, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('BCPHOBIC')
            call ESMFL_StateGetPointerToData(import, bcphobicp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('BCPHILIC')
            call ESMFL_StateGetPointerToData(import, bcphilicp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('OCPHOBIC')
            call ESMFL_StateGetPointerToData(import, ocphobicp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case ('OCPHILIC')
            call ESMFL_StateGetPointerToData(import, ocphilicp, trim(cvar), rc=STATUS)
            VERIFY_(STATUS)
         case default
            if(mype==0) write(6,*) trim(Iam), ': ', trim(cvar), ' no such chem (imp) variable, aborting ...'
            status = 1
            VERIFY_(STATUS)
      end select
   enddo

! exports
! -------

   call ESMFL_StateGetPointerToData(export, dts,               'ts',  alloc=.true., rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(export, dfrland,        'frland', alloc=.true., rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(export, dfrlandice,  'frlandice', alloc=.true., rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(export, dfrlake,        'frlake', alloc=.true., rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(export, dfrocean,      'frocean', alloc=.true., rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(export, dfrseaice,    'frseaice', alloc=.true., rc=STATUS)
   VERIFY_(STATUS)
   call ESMFL_StateGetPointerToData(export, ddp,             'delp',  alloc=.true., rc=STATUS)
   VERIFY_(STATUS)

! Extra meteorological fields
! ---------------------------
   do nt=1,nmguess
      cvar = trim(mguess_usr(nt))
      select case (cvar)
         case ('phis')
            call ESMFL_StateGetPointerToData(export,  dhs,  trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('ps')
            call ESMFL_StateGetPointerToData(export,  dps,  trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('u')
            call ESMFL_StateGetPointerToData(export,    du, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('v')
            call ESMFL_StateGetPointerToData(export,    dv, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('tv')
            call ESMFL_StateGetPointerToData(export,    dt, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('sphu')
            call ESMFL_StateGetPointerToData(export,    dq, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('qitot')
            call ESMFL_StateGetPointerToData(export, dqimr, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('qltot')
            call ESMFL_StateGetPointerToData(export, dqlmr, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('qrtot')
            call ESMFL_StateGetPointerToData(export, dqrmr, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('qstot')
            call ESMFL_StateGetPointerToData(export, dqsmr, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('ozone')
           call ESMFL_StateGetPointerToData(export, doz,    trim(cvar), alloc=.true., rc=STATUS)
           VERIFY_(STATUS)
         case default
            if(mype==0) write(6,*) trim(Iam), ': ', trim(cvar), ' no such met-guess (ex) variable, skipping ...'
!_RT        if(mype==0) write(6,*) trim(Iam), ': no such met-guess (ex) variable, aborting ...'
!_RT        status = 1
!_RT        VERIFY_(STATUS)
      end select
   enddo

! Chemistry tracer exports
! When proper connection w/ Tracer Bundel is made
! the following won't be necessary
! -----------------------------------------------
   do nt=1,ntgases
      cvar = trim(tgases_usr(nt))
!_RT      call gsi_chemguess_get ( 'i4crtm::'//trim(cvar), i4crtm, status )
!_RT      if(i4crtm>=10) cycle ! will ignore aerosols for now
      select case (cvar)
         case ('AOD')
            call ESMFL_StateGetPointerToData(export, daodp, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('CO')
            call ESMFL_StateGetPointerToData(export, dcop, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case ('CO2')
            call ESMFL_StateGetPointerToData(export, dco2p, trim(cvar), alloc=.true., rc=STATUS)
            VERIFY_(STATUS)
         case default
            if(mype==0) write(6,*) trim(Iam), ': ', trim(cvar), ' not available in chem (ex), skipping ...'
      end select
   enddo


   end subroutine GSI_GridCompGetPointers_

   subroutine Scale_Import_()
   character(len=*), parameter :: IAm='GSI_GridComp.Scale_Import_'
   integer(i_kind) :: nt
   character(len=ESMF_MAXSTR) :: cvar

   where ( psp /= MAPL_UNDEF )
           psp = psp * kPa_per_Pa       ! convert ps to cb
   endwhere
   if(GsiGridType==0) then
      if(associated(hsp)) then
         where ( hsp /= MAPL_UNDEF )
                 hsp = hsp / grav       ! convert geop h to h
         endwhere
      endif
      if(associated(ozp)) then
         where ( ozp /= MAPL_UNDEF )
                 ozp = ozp * PPMV2GpG  ! convert from ppmv to g/g
         endwhere
      endif
      if(associated(cop)) then
         where ( cop /= MAPL_UNDEF )
                 cop = cop * KGpKG2PPBV  ! convert carbon monoxide unit (need gen. way of handling chemistry)
         endwhere
      endif
      if(associated(co2p)) then
         where ( co2p /= MAPL_UNDEF )
                 co2p = co2p * KGpKG2ppmvCO2  ! convert carbon dioxide to ppmv
         endwhere
      endif
      if(associated(du001p)) then
         where ( du001p /= MAPL_UNDEF )
                 du001p = du001p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(du002p)) then
         where ( du002p /= MAPL_UNDEF )
                 du002p = du002p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
         print *, 'DEBUG doing the right thing ...'
      endif
      if(associated(du003p)) then
         where ( du003p /= MAPL_UNDEF )
                 du003p = du003p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(du004p)) then
         where ( du004p /= MAPL_UNDEF )
                 du004p = du004p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(du005p)) then
         where ( du005p /= MAPL_UNDEF )
                 du005p = du005p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(ss001p)) then
         where ( ss001p /= MAPL_UNDEF )
                 ss001p = ss001p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(ss002p)) then
         where ( ss002p /= MAPL_UNDEF )
                 ss002p = ss002p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(ss003p)) then
         where ( ss003p /= MAPL_UNDEF )
                 ss003p = ss003p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(ss004p)) then
         where ( ss004p /= MAPL_UNDEF )
                 ss004p = ss004p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(ss005p)) then
         where ( ss005p /= MAPL_UNDEF )
                 ss005p = ss005p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(dmsp)) then
         where ( dmsp /= MAPL_UNDEF )
                 dmsp = dmsp * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(so2p)) then
         where ( so2p /= MAPL_UNDEF )
                 so2p = so2p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(so4p)) then
         where ( so4p /= MAPL_UNDEF )
                 so4p = so4p * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(msap)) then
         where ( msap /= MAPL_UNDEF )
                 msap = msap * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(bcphobicp)) then
         where ( bcphobicp /= MAPL_UNDEF )
                 bcphobicp = bcphobicp * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(bcphilicp)) then
         where ( bcphilicp /= MAPL_UNDEF )
                 bcphilicp = bcphilicp * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(ocphobicp)) then
         where ( ocphobicp /= MAPL_UNDEF )
                 ocphobicp = ocphobicp * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
      if(associated(ocphilicp)) then
         where ( ocphilicp /= MAPL_UNDEF )
                 ocphilicp = ocphilicp * KGpKG2PPBVaero ! convert from kg/kg to ppbv
         endwhere
      endif
   endif

   end subroutine Scale_Import_ 
!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: GSI_GridCompCopyImportDyn2Internal  -- get vars from impGSI

! !INTERFACE:

   subroutine GSI_GridCompCopyImportDyn2Internal_(lit)

   use mpeu_util, only: tell
   implicit none

   integer, intent(in) :: lit ! logical index of first guess time level to copy

! !DESCRIPTION: Extract all upper air fields from GSI import state to
!                define the corresponding GSI internal variables
!
! !REVISION HISTORY:
!
!  17Apr2007 Todling  Swap vertical on way in
!  10Mar2009 Todling  Handle qi/qlmr instead of qctot
!  07May2011 Todling/Merkova - add cloud-fraction for radiation to import
!  24Oct2011 Todling/Akella  - add NST variables
!
!EOPI
!-------------------------------------------------------------------------

! local variables

   character(len=*), parameter    :: &
            IAm='GSI_GridCompCopyImportDyn2Internal_'
   integer(i_kind) :: it,nn,nt,itguessig
   integer(i_kind) :: irank, ipnt, ier,istatus
   character(len=ESMF_MAXSTR) :: cvar
   real(r_kind),pointer,dimension(:,:,:):: ptr_cwmr,ptr_qimr,ptr_qlmr


! start   

   if(IamRoot) print *,trim(Iam),': Copy contents of import to internal state, it= ', lit

! The upper air fields from the GEOS import state are already decomposed 
! for the GSI. A final operation before updating the corresponding GSI 
! internal fields is to transpose the horizontal grid  from IJK (GEOS) 
! to JIK (GSI). TODO: this could be done during coupling.

   call Scale_Import_()

   nfldsig_now = lit
   if(lit > nfldsig) then	! rotate the storage

     nfldsig_now = nfldsig	! take over the last storage slot

     do nn = 1,nfldsig-1	! by first moving data slots up

       hrdifsig(nn) = hrdifsig(nn+1)

       if(nmguess>0) GSI_MetGuess_bundle(nn)= GSI_MetGuess_bundle(nn+1)
       if(ntgases>0) GSI_chemguess_bundle(nn)= GSI_chemguess_bundle(nn+1)
     enddo
   endif

   it = nfldsig_now
   ntguessig = ntguessig_ref - lit + it
   itguessig = max(1,min(ntguessig,nfldsig))
   if(IamRoot) call tell(Iam, &
   	'reset ntguessig, (from, to) =',(/ntguessig,itguessig/))
   ntguessig = itguessig

   ntguesnst = ntguessig ! make life simple for now
   hrdifnst  = hrdifsig  ! same here

#ifdef VERBOSE
   call tell(Iam,'nfldsig.it  =',it)
   call tell(Iam,'nfldsig.lit =',lit)
   call tell(Iam,'nfldsig     =',nfldsig)
#endif
   hrdifsig(it)=hrdifsig_all(lit)
#ifdef VERBOSE
   call tell(Iam,'hrfldsig(it) =',hrdifsig(it))
   call tell(Iam,"ntguessig =",ntguessig)
#endif

!  Handle extra meteological guess fields
!  -------------------------------------
   do nt = 1,nmguess
      cvar=trim(mguess_gsi(nt))
      call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), cvar, ipnt, status, irank=irank )
      VERIFY_(STATUS)
      if(irank==2) then
         select case (cvar)
            case ('z')
              call GSI_GridCompSwapIJ_(hsp,GSI_MetGuess_bundle(it)%r2(ipnt)%q)
            case ('ps')
              call GSI_GridCompSwapIJ_(psp,GSI_MetGuess_bundle(it)%r2(ipnt)%q)
            case ('z_c')
              call GSI_GridCompSwapIJ_(z_cp,GSI_MetGuess_bundle(it)%r2(ipnt)%q)
            case ('z_w')
              call GSI_GridCompSwapIJ_(z_wp,GSI_MetGuess_bundle(it)%r2(ipnt)%q)
            case ('dt_cool')
              call GSI_GridCompSwapIJ_(dt_coolp,GSI_MetGuess_bundle(it)%r2(ipnt)%q)
            case ('tdel')
              call GSI_GridCompSwapIJ_(tdelp,GSI_MetGuess_bundle(it)%r2(ipnt)%q)
            case ('tref')
              call GSI_GridCompSwapIJ_(trefp ,GSI_MetGuess_bundle(it)%r2(ipnt)%q)
         end select
#ifdef SFCverbose
         call guess_grids_stats(cvar, GSI_MetGuess_Bundle(it)%r2(ipnt)%q, mype)
#endif
      endif
      if(irank==3) then
         select case (cvar)
            case ('u')
              call GSI_GridCompSwapIJ_(up,  GSI_MetGuess_bundle(it)%r3(ipnt)%q)
            case ('v')
              call GSI_GridCompSwapIJ_(vp,  GSI_MetGuess_bundle(it)%r3(ipnt)%q)
            case ('tv')
              call GSI_GridCompSwapIJ_(tp,  GSI_MetGuess_bundle(it)%r3(ipnt)%q)
            case ('q')
              call GSI_GridCompSwapIJ_(qp,  GSI_MetGuess_bundle(it)%r3(ipnt)%q)
            case ('oz')
              call GSI_GridCompSwapIJ_(ozp, GSI_MetGuess_bundle(it)%r3(ipnt)%q)
            case ('qi')
              call GSI_GridCompSwapIJ_(qimr,GSI_MetGuess_bundle(it)%r3(ipnt)%q)
            case ('ql')
              call GSI_GridCompSwapIJ_(qlmr,GSI_MetGuess_bundle(it)%r3(ipnt)%q)
            case ('qr')
              call GSI_GridCompSwapIJ_(qrmr,GSI_MetGuess_bundle(it)%r3(ipnt)%q)
            case ('qs')
              call GSI_GridCompSwapIJ_(qsmr,GSI_MetGuess_bundle(it)%r3(ipnt)%q)
            case ('cf')
              call GSI_GridCompSwapIJ_(clfr,GSI_MetGuess_bundle(it)%r3(ipnt)%q)
         end select
#ifdef UPAverbose
         call guess_grids_stats(cvar, GSI_MetGuess_Bundle(it)%r3(ipnt)%q, mype)
#endif
      endif
   enddo

! Cloud liquid water content and cloud ice water content are added into
! cloud condensate mixing ratio - internal GSI field
   ier=0
   call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'cw', ptr_cwmr, istatus );ier=ier+istatus
   call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qi', ptr_qimr, istatus );ier=ier+istatus
   call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'ql', ptr_qlmr, istatus );ier=ier+istatus
   if (ier==0) then
      where(ptr_qimr.ne.MAPL_UNDEF .and. ptr_qlmr.ne.MAPL_UNDEF)
            ptr_cwmr=ptr_qimr+ptr_qlmr
      elsewhere
            ptr_cwmr=MAPL_UNDEF
      endwhere 
   endif

!  Handle trace gases
!  ------------------
   do nt = 1,ntgases
      cvar=trim(tgases_gsi(nt))
      call GSI_BundleGetPointer ( GSI_chemguess_bundle(it), cvar, ipnt, status, irank=irank )
      VERIFY_(STATUS)
      if (irank==2) then
         select case (cvar)
           case ('aod')
             call GSI_GridCompSwapIJ_(aodp,GSI_chemguess_bundle(it)%r2(ipnt)%q)
         end select
      endif
      if (irank==3) then
         select case (cvar)
           case ('co')
             call GSI_GridCompSwapIJ_(cop,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('co2')
             call GSI_GridCompSwapIJ_(co2p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('du001')
             call GSI_GridCompSwapIJ_(du001p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('du002')
             call GSI_GridCompSwapIJ_(du002p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('du003')
             call GSI_GridCompSwapIJ_(du003p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('du004')
             call GSI_GridCompSwapIJ_(du004p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('du005')
             call GSI_GridCompSwapIJ_(du005p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('ss001')
             call GSI_GridCompSwapIJ_(ss001p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('ss002')
             call GSI_GridCompSwapIJ_(ss002p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('ss003')
             call GSI_GridCompSwapIJ_(ss003p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('ss004')
             call GSI_GridCompSwapIJ_(ss004p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('ss005')
             call GSI_GridCompSwapIJ_(ss005p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('dms')
             call GSI_GridCompSwapIJ_(dmsp,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('so2')
             call GSI_GridCompSwapIJ_(so2p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('so4')
             call GSI_GridCompSwapIJ_(so4p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('msa')
             call GSI_GridCompSwapIJ_(msap,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('bcphobic')
             call GSI_GridCompSwapIJ_(bcphobicp,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('bcphilic')
             call GSI_GridCompSwapIJ_(bcphilicp,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('ocphobic')
             call GSI_GridCompSwapIJ_(ocphobicp,GSI_chemguess_bundle(it)%r3(ipnt)%q)
           case ('ocphilic')
             call GSI_GridCompSwapIJ_(ocphilicp,GSI_chemguess_bundle(it)%r3(ipnt)%q)
         end select
      end if ! rank 3
#ifdef UPAverbose
       call guess_grids_stats(cvar, GSI_chemguess_bundle(it)%r3(ipnt)%q, mype)
#endif
   enddo

!  If so, this allows overwriting CO2 (and other trace gases) in background 
!  ------------------------------------------------------------------------
   call read_gfs_chem (iadate(1),iadate(2),iadate(3),it)

! simple statistics

#ifdef VERBOSE
   if(IamRoot) print *,trim(Iam),': Complete copy contents of import to internal state, it= ', lit
#endif

   end subroutine GSI_GridCompCopyImportDyn2Internal_

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: GSI_GridCompComputeVorDiv_

! !INTERFACE:

   subroutine GSI_GridCompComputeVorDiv_(lit)

#ifdef _GMAO_FVGSI_
   use m_ggGradient,only : ggGradient
   use m_ggGradient,only : ggGradient_init,clean
   use m_ggGradient,only : ggDivo
#else /* _GMAO_FVGSI_ */
   use compact_diffs, only: cdiff_created
   use compact_diffs, only: cdiff_initialized
   use compact_diffs, only: create_cdiff_coefs
   use compact_diffs, only: inisph
   use compact_diffs, only: uv2vordiv
   use xhat_vordivmod, only: xhat_vordiv_calc2
#endif /* _GMAO_FVGSI_ */

   implicit none

   integer, intent(in) :: lit ! logical index of first guess time level to copy

! !DESCRIPTION: vorticity-divergence calculation 
!
! !REVSION HISTORY:
!
!  18Feb2009 Todling - let GSI calc of vor/div when GMAO util not available
!
!EOPI
!-------------------------------------------------------------------------

! local variables

   character(len=*), parameter :: IAm='GSI_GridCompComputeVorDiv_'
   integer(i_kind) :: it,nn,ier,istatus
#ifdef _GMAO_FVGSI_
   type(ggGradient) :: gr
   integer :: iGdim,iGloc,iGlen
   integer :: jGdim,jGloc,jGlen
   real(r_double),dimension(:,:,:), pointer :: up8
   real(r_double),dimension(:,:,:), pointer :: vp8
#endif /* _GMAO_FVGSI_ */
   real(r_kind),dimension(:,:,:),pointer :: ges_u_it,ges_v_it
   real(r_kind),dimension(:,:,:),pointer :: ges_div_nnn,ges_div_np1
   real(r_kind),dimension(:,:,:),pointer :: ges_vor_nnn,ges_vor_np1

! start

   if(IamRoot.and.verbose) print *,trim(Iam),': Compute vorticity and divergence'

   ier=0
   it = nfldsig_now
   if(lit > nfldsig) then	! rotate the storage
     do nn=1,nfldsig-1
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(nn)  , 'div', ges_div_nnn, istatus );ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(nn)  , 'vor', ges_vor_nnn, istatus );ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(nn+1), 'div', ges_div_np1, istatus );ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(nn+1), 'vor', ges_vor_np1, istatus );ier=ier+istatus
       if(ier/=0) exit ! get out of here ...
       ges_div_nnn = ges_div_np1
       ges_vor_nnn = ges_vor_np1
     enddo
   endif
   if(ier/=0) return ! nothing to do

#ifdef _GMAO_FVGSI_
   allocate(up8 (lat2,lon2,nsig), &
            vp8 (lat2,lon2,nsig), &
            stat=STATUS)
   VERIFY_(STATUS)

! define -real8 and swapped- uwnd,vwnd needed by ggDivo
  
   if(associated(up).and.associated(vp)) then
      call GSI_GridCompSwapIJ_(up,up8)
      call GSI_GridCompSwapIJ_(vp,vp8)
   else
      deallocate(up8,vp8)
      return
   endif

   call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'div', ges_div_nnn, istatus );ier=ier+istatus
   call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'vor', ges_vor_nnn, istatus );ier=ier+istatus
   if(ier/=0) then
      deallocate(up8,vp8)
      return  ! nothing to do
   endif
   ges_vor_nnn = zero
   ges_div_nnn = zero

   iGloc=jstart(myPE+1)
   iGlen=jlon1 (myPE+1)
   jGloc=istart(myPE+1)
   jGlen=ilat1 (myPE+1)

   call ggGradient_init(gr,nlon,real(rlats(1:nlat),r_double),	&
  	iGloc,iGlen,jGloc,jGlen,nsig, mpi_comm_world)

   if(IamRoot.and.verbose) print *,Iam,': call ggDivo '
   call ggDivo(gr,		&
	up8  (:,:,:),	        &	! in: u
	vp8  (:,:,:),	        &	! in: v
	ges_div_nnn,		&	! div(u,v)
	ges_vor_nnn,		&	! vor(u,v)
	mpi_comm_world)		        ! communicator

   call clean(gr)
   deallocate(up8, vp8, stat=STATUS)
   VERIFY_(STATUS)

#else /* _GMAO_FVGSI_ */
   if(.not.cdiff_created()) call create_cdiff_coefs()
   if(.not.cdiff_initialized()) call inisph(rearth,rlats(2),wgtlats(2),nlon,nlat-2)
   ier=0
   call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'u', ges_u_it, istatus );ier=ier+istatus
   call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'v', ges_v_it, istatus );ier=ier+istatus
   if(ier==0) then
      call xhat_vordiv_calc2 (ges_u_it,ges_v_it,ges_vor_nnn,ges_div_nnn)
   endif
!!   call destroy_cdiff_coefs
#endif /* _GMAO_FVGSI_ */

! simple statistics

#ifdef UPAverbose 
   call guess_grids_stats('GCges_vor', ges_vor_nnn, mype)
   call guess_grids_stats('GCges_div', ges_div_nnn, mype)
#endif

   end subroutine GSI_GridCompComputeVorDiv_

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: GSI_GridCompCopyImportSfc2Internal_

! !INTERFACE:

   subroutine GSI_GridCompCopyImportSfc2Internal_(lit)

   implicit none

   integer, intent(in) :: lit ! logical index of first guess time level to copy

! !DESCRIPTION: Extract all surface fields from GSI import state to
!                define the corresponding GSI internal variables
! !REVISION HISTORY:
! 
!  24Apr2007 Todling  Renamed tskin to ts (getting from upa bkg file)
!  08May2007 Todling  Corrected calculation of fact10; do not replace
!                     soil temp w/ tskin - leave undefs in place; revised
!                     calculation of isli
!
!EOPI
!-------------------------------------------------------------------------

! local variables

   character(len=*), parameter :: &
            IAm='GSI_GridCompCopyImportSfc2Internal_'
   real(r_kind)                      :: wspd
   integer,allocatable,dimension(:,:):: islip
   integer(i_kind) :: it,nn,itguessfc,ier,istatus
   real(r_kind),dimension(:,:,:),pointer :: ges_u_it,ges_v_it

! start

   if(IamRoot) print *,trim(Iam),': Copy contents of import to internal state, it= ',lit

   nfldsfc_now = lit
   if(lit > nfldsfc) then
     nfldsfc_now = nfldsfc	! take over the last storage slot

     do nn=1,nfldsfc-1
       hrdifsfc(nn) = hrdifsfc(nn+1)

       dsfct(:,:,nn) = dsfct(:,:,nn+1)
        sfct(:,:,nn) =  sfct(:,:,nn+1)
	sno (:,:,nn) =  sno (:,:,nn+1)

       soil_moi (:,:,nn) = soil_moi (:,:,nn+1)
       soil_temp(:,:,nn) = soil_temp(:,:,nn+1)
       sfc_rough(:,:,nn) = sfc_rough(:,:,nn+1)

        isli(:,:,nn) =   isli(:,:,nn+1)
      fact10(:,:,nn) = fact10(:,:,nn+1)
     enddo
   endif

   it = nfldsfc_now
   ntguessfc = ntguessfc_ref - lit + it
   itguessfc = max(1,min(ntguessfc,nfldsfc))
   if(IamRoot) call tell(Iam, &
   	'reset ntguessfc, (from, to) =',(/ntguessfc,itguessfc/))
   ntguessfc = itguessfc

#ifdef VERBOSE
   call tell(Iam,'nfldsfc.it  =',it)
   call tell(Iam,'nfldsfc.lit =',lit)
#endif
   hrdifsfc(it)=hrdifsfc_all(lit)
#ifdef VERBOSE
   call tell(Iam,'hrfldsfc(it) =',hrdifsfc(it))
   call tell(Iam,"ntguessfc =",ntguessfc)
#endif

! The surface fields from the GEOS import state are already decomposed 
! for the GSI. A final operation before updating the corresponding GSI 
! internal fields is to transpose the horizontal grid  from IJK (GEOS) 
! to JIK (GSI).

! For consistency with the GMAO_INTFC results the surface fields have been
! "processed" in a similar fashion using undef_2ssi

! Skin temperature

   dsfct(:,:,it) = zero

   call GSI_GridCompSwapIJ_(tskp,sfct(:,:,it))

! Snow depth
   call undef_2ssi(snop,MAPL_UNDEF,trim(Iam),	&
        verb=.true.,vname='SNOWDP')
   call undef_2ssi(snop,UNDEF_SNOW_DEP,trim(Iam),	&
        verb=.false.,vname='SNOWDP')
   where(snop==UNDEF_SSI_)
      snop=0.
   elsewhere(snop <zero)
      snop=0.
   endwhere
   call GSI_GridCompSwapIJ_(snop,sno(:,:,it))

! Soil moisture
   call undef_2ssi(soqp,MAPL_UNDEF,trim(Iam),	&
        verb=.false.,vname='GWETTOP')
   where(soqp == UNDEF_SSI_)	! in case of any
      soqp=1.
   elsewhere(soqp <=zero)		! why not <.05?
      ! This block would undo any :=UNDEF_SSI_ == -9.99e33
      soqp=.05
   elsewhere(soqp > one)
      soqp=1.
   endwhere
   call GSI_GridCompSwapIJ_(soqp,soil_moi(:,:,it))

! Soil temperature
   call undef_2ssi(sotp,MAPL_UNDEF,trim(Iam),	&
        verb=.false.,vname='TSOIL1')
   call undef_2ssi(sotp,UNDEF_SOIL_TEMP,trim(Iam),	&
        verb=.true.,vname='TSOIL1')
! A hack solution to work around possible undesired undef values.
   where(sotp == UNDEF_SSI_)
      sotp=tskp
   endwhere

   call GSI_GridCompSwapIJ_(sotp,soil_temp(:,:,it))

! Surface roughness
   call undef_2ssi(sz0p,MAPL_UNDEF,trim(Iam),	&
        verb=.false.,vname='Z0M')
! A hack solution to work around possible undesired undef values.
   where(sz0p == UNDEF_SSI_ .or. sz0p<zero)
      sz0p=0.
   endwhere

   call GSI_GridCompSwapIJ_(sz0p,sfc_rough(:,:,it))

! SLI mask  (adaptation from L. Takacs change elsewhere)
   allocate(islip(lon2,lat2), stat=STATUS)
   VERIFY_(STATUS)
                                            islip = 1  ! Land
   where (  frocean+frlake >= 0.6         ) islip = 0  ! Water
   where (  islip==0 .and. frseaice > 0.5 ) islip = 2  ! Ice
   where (  islip==0 .and.   tskp < 271.4 ) islip = 2  ! Ice

! Set output (export) land/water/etc information
! RT: Is there a more decent way to set outputs from inputs in ESMF?
!     Why can't I use fr-arrays directly?
!   if ( it==ntguessfc ) then
!
! Check stored hrdifsfc(:) for expected time (J.G.)
   if ( hrdifsfc(it)==hrdifsfc_all(ntguessfc_ref) ) then
       dts        = tskp
       dfrland    = frland
       dfrlandice = frlandice
       dfrlake    = frlake
       dfrocean   = frocean
       dfrseaice  = frseaice
   endif

   call GSI_GridCompSwapIJ_(islip,isli(:,:,it))

   deallocate(islip, stat = STATUS)
   VERIFY_(STATUS)

   if(GSIGridType==0) then ! we need to calculate 10m wind factors
                           ! else they are provided in the import state

! F10M is derived from U10M,V10m,U,V and defined through GSI's
! internal variable fact10 (declared in guess_grids.f90)

      call undef_2ssi(u10p,MAPL_UNDEF,trim(Iam),	&
           verb=.false.,vname='U10M')
      call undef_2ssi(v10p,MAPL_UNDEF,trim(Iam),	&
           verb=.false.,vname='V10M')
#ifdef SFCverbose
      print *,Iam,':  u10m (PE,min,max,sum): ',mype,minval(u10p),maxval(u10p),sum(u10p)
      print *,Iam,':  v10m (PE,min,max,sum): ',mype,minval(v10p),maxval(v10p),sum(v10p)
#endif
      allocate(f10p(lon2,lat2), stat=STATUS)
      VERIFY_(STATUS)
   
      where(u10p == UNDEF_SSI_ .or. v10p == UNDEF_SSI_)	        ! in case there is any
        f10p(:,:)=one ! RT: this would happen over the halo only anyway 
      elsewhere
        f10p(:,:)=sqrt(u10p(:,:)*u10p(:,:)+v10p(:,:)*v10p(:,:))
      endwhere

      call GSI_GridCompSwapIJ_(f10p,fact10(:,:,it))

!     zero somewhere.  To compute factor at 10m, wind speed at the
!     bottom of the atmosphere is needed, which is level 1 of the
!     GSI grid.  It is debateable what fact10 should be where the
!     denominator is zero.  I don't have any good solution, but
!     simply mimic what has been done before.
                                                                                                                                      
      ier=0
      call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'u', ges_u_it, istatus );ier=ier+istatus
      call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'v', ges_v_it, istatus );ier=ier+istatus
      if (ier==0) then
         do j=lbound(fact10,2),ubound(fact10,2)
            do i=lbound(fact10,1),ubound(fact10,1)
               wspd=sqrt(ges_u_it(i,j,1)*ges_u_it(i,j,1) + &
                         ges_v_it(i,j,1)*ges_v_it(i,j,1))
                                                                                                                                    
!              Reset fact10 if fact10 > zero
               if(zero < fact10(i,j,it)) then
                  if(fact10(i,j,it) <= wspd) then
                     fact10(i,j,it)=fact10(i,j,it)/wspd
                  else
                     fact10(i,j,it)=one
                  endif
               endif
            end do
         end do
      endif
      deallocate(f10p, stat = STATUS)
      VERIFY_(STATUS)

   end if ! GridType

! simple statistics

#ifdef SFCverbose
   call guess_grids_stats('GCsfct',      sfct     (:,:,it), mype)
   call guess_grids_stats('GCsno',       sno      (:,:,it), mype)
   call guess_grids_stats('GCsoil_moi',  soil_moi (:,:,it), mype)
   call guess_grids_stats('GCsfc_rough', sfc_rough(:,:,it), mype)
   call guess_grids_stats('GCsoil_temp', soil_temp(:,:,it), mype)
   call guess_grids_stats('GCfact10',    fact10   (:,:,it), mype)
#endif


   end subroutine GSI_GridCompCopyImportSfc2Internal_

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: GSI_GridCompGetNCEPsfcFromFile_

! !INTERFACE:

   subroutine GSI_GridCompGetNCEPsfcFromFile_(lit)

   use sfcio_module
   implicit none

   integer, intent(in) :: lit ! logical index of first guess time level to copy

! !DESCRIPTION: read NCEP global surface fields, vegetation fraction,
!               vegetation and soil type,  from a file.
!   17Mar2016  Todling  Revisit handling of extra surface fields: use sfcio;
!                       add option to set veg-type to constant for sensitivity
!                       evalulation.
!
!EOPI
!-------------------------------------------------------------------------

! local variables

   integer                              :: glon,glat,glat2
   real(r_single), allocatable, dimension(:,:) :: vfrbuf
   real(r_single), allocatable, dimension(:,:) :: vtybuf
   real(r_single), allocatable, dimension(:,:) :: stybuf
   real(r_single), allocatable, dimension(:,:) :: sfcbuf
   real(r_single), allocatable, dimension(:,:) :: sfcbufpp
   real(r_single) yhour
   integer version
   integer,dimension(4):: igdate

   type(sfcio_head):: head
   type(sfcio_data):: data

   character(len=*), parameter    :: &
            IAm='GSI_GridCompGetNCEPsfcFromFile_'

   integer(i_kind) :: it,nn,iret
   integer(i_kind) :: iset_veg_type
   logical,parameter :: wrt_ncep_sfc_grd=.false.
   logical,parameter :: use_sfcio=.true.

! start

   if(IamRoot.and.verbose) print *,trim(Iam),': read NCEP global surface fields from a file'

   call ESMF_ConfigGetAttribute( CF, iset_veg_type, label ='THIS_VEG_TYPE:', default=-1, rc = STATUS )
   VERIFY_(STATUS)

   if(lit > nfldsfc) then
     do nn=1,nfldsfc-1
        veg_type(:,:,nn) = veg_type(:,:,nn+1)
       soil_type(:,:,nn) =soil_type(:,:,nn+1)   
        veg_frac(:,:,nn) = veg_frac(:,:,nn+1)
     enddo
   endif
   it = nfldsfc_now

! Now we use hinterp to regrid ncepsfc fields as well as swap lats, add poles
! and flip longitudes before using hinterp (also flip lons after hinterp).
!
   if(IamRoot) then

      if (wrt_ncep_sfc_grd) then
         call baopenwt(36,'ncepsfc.grd',iret)
      endif
      ! NCEP surface fields are defined on a Gaussian grid...
      if (use_sfcio) then
         print *,trim(Iam),': Using sfcio to read NCEP surface file (as BC)'
         call sfcio_srohdc(34,'ncepsfc',head,data,iret)
         yhour=head%fhour
         igdate=head%idate
         glat2=head%latb
         glon =head%lonb
         version = head%ivs
         glat=glat2
         if (iset_veg_type>0) then
            data%vfrac=1.0 ! test sens of CRTM to veg fractions
            data%vtype=13.0! test sens of CRTM to veg type
            print*, trim(Iam), ', Veg-fractions reset: vfrac=1.0'
            print*, trim(Iam), ', Veg-type reset by request to: ',iset_veg_type
         endif
      else
         open(34,file='ncepsfc', form='unformatted')
         read (34) 
         read(34) yhour,igdate,glon,glat2,version
         close(34)
         glat=glat2+2 
      endif
      print *,trim(Iam),': ncepsfc hour/date  : ',yhour,' / ',igdate
      print *,trim(Iam),': ncepsfc fields dims, version: ',glon,' x ',glat2, ' v', version

      allocate(sfcbuf(glon,glat2),  & ! ncepsfc fields buffer
               sfcbufpp(glon,glat), & ! same but with poles added
               stat=STATUS)  
      VERIFY_(STATUS)
      allocate(vtybuf(nlon,nlat), &
               stybuf(nlon,nlat), &
               vfrbuf(nlon,nlat), & 
               stat=STATUS)
      VERIFY_(STATUS)

      if (.not.use_sfcio ) then
         call ncep_rwsurf_ ( .false., 'ncepsfc', glat2, glon, &
                             STATUS, jrec=12, fld=sfcbuf)
         VERIFY_(STATUS)
      else
         sfcbuf=data%vfrac
      endif
      call GSI_GridCompSP2NP_(sfcbufpp,sfcbuf)
      call GSI_GridCompFlipLons_(sfcbufpp)
      call hinterp2 ( sfcbufpp ,glon, glat   , &
                      vfrbuf, nlon, nlat, 1, &
                      UNDEF_SSI_)
      call GSI_GridCompFlipLons_(vfrbuf)
      where(vfrbuf<0.0)
        vfrbuf=0.0
      end where
      where(vfrbuf>1.0)
        vfrbuf=1.0
      end where
      if (wrt_ncep_sfc_grd) then
         call wryte(36,4*nlat*nlon,vfrbuf)
      endif

      if (.not.use_sfcio ) then
         call ncep_rwsurf_ ( .false., 'ncepsfc', glat2, glon, &
                             STATUS, jrec=15, fld=sfcbuf)
         VERIFY_(STATUS)
      else
         sfcbuf=data%vtype
      endif
      call GSI_GridCompSP2NP_(sfcbufpp,sfcbuf)
      call GSI_GridCompFlipLons_(sfcbufpp)
      call hinterp2 ( sfcbufpp ,glon, glat   , &
                      vtybuf, nlon, nlat, 1, &
                      UNDEF_SSI_)
      call GSI_GridCompFlipLons_(vtybuf)
      vtybuf=real(nint(vtybuf),r_single) 
      where(vtybuf <  0.0_r_single) vtybuf =  0.0_r_single
      where(vtybuf > 13.0_r_single) vtybuf = 13.0_r_single
      if (wrt_ncep_sfc_grd) then
         call wryte(36,4*nlat*nlon,vtybuf)
      endif

      if (.not.use_sfcio ) then
         call ncep_rwsurf_ ( .false., 'ncepsfc', glat2, glon, &
                             STATUS, jrec=16, fld=sfcbuf)
         VERIFY_(STATUS)
      else
         sfcbuf=data%stype
      endif
      call GSI_GridCompSP2NP_(sfcbufpp,sfcbuf)
      call GSI_GridCompFlipLons_(sfcbufpp)
      call hinterp2 ( sfcbufpp ,glon, glat   , &
                      stybuf, nlon, nlat, 1, &
                      UNDEF_SSI_)
      call GSI_GridCompFlipLons_(stybuf)
      stybuf=real(nint(stybuf),r_single)
      where(stybuf <  0.0_r_single) stybuf = 0.0_r_single
      where(stybuf >  9.0_r_single) stybuf = 9.0_r_single
      if (wrt_ncep_sfc_grd) then
         call wryte(36,4*nlat*nlon,stybuf)
         call baclose(36,iret) 
      endif

      deallocate(sfcbuf, sfcbufpp, stat=STATUS)
      VERIFY_(STATUS)

   end if
   
   allocate(vfrp(lon2,lat2), stat=STATUS)
   VERIFY_(STATUS)
   call ArrayScatter(vfrp, vfrbuf, GSIGrid, hw=hw, rc=STATUS)
   VERIFY_(STATUS)
   call GSI_GridCompSwapIJ_(vfrp, veg_frac(:,:,it))
   deallocate(vfrp, stat = STATUS)
   VERIFY_(STATUS)

   allocate(vtyp(lon2,lat2), stat=STATUS)
   VERIFY_(STATUS)
   call ArrayScatter(vtyp, vtybuf, GSIGrid, hw=hw, rc=STATUS)
   VERIFY_(STATUS)
!  Dirty fix to get rid of complete insanity (but still not correct)
!  -----------------------------------------
   where(vtyp <= zero) vtyp = zero
   call GSI_GridCompSwapIJ_(vtyp, veg_type(:,:,it))
   deallocate(vtyp, stat = STATUS)
   VERIFY_(STATUS)

   allocate(styp(lon2,lat2), stat=STATUS)
   VERIFY_(STATUS)
   call ArrayScatter(styp, stybuf, GSIGrid, hw=hw, rc=STATUS)
   VERIFY_(STATUS)
!  Dirty fix to get rid of complete insanity (but still not correct)
!  -----------------------------------------
   where(styp <= zero) styp = zero
   call GSI_GridCompSwapIJ_(styp,soil_type(:,:,it))   
   deallocate(styp, stat = STATUS)
   VERIFY_(STATUS)

   if(IamRoot) then
      deallocate(vtybuf, stybuf, vfrbuf, stat=STATUS)
      VERIFY_(STATUS)
   end if

#ifdef SFCverbose
   call guess_grids_stats('GCsoil_type', soil_type(:,:,it), mype) 
   call guess_grids_stats('GCveg_type',  veg_type (:,:,it), mype)
   call guess_grids_stats('GCveg_frac',  veg_frac (:,:,it), mype)
#endif

   end subroutine GSI_GridCompGetNCEPsfcFromFile_
   
!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: GSI_GridCompCopyInternal2Export_  -- copy vars to expGSI

! !INTERFACE:

   subroutine GSI_GridCompCopyInternal2Export_(lit)

   use constants,only: one
   use gsi_nstcouplermod,  only: nst_gsi
   implicit none

   integer, intent(in) :: lit ! logical index of first guess time level to copy

! !DESCRIPTION: Extract GSI internal variables and put into GSI
!                export state.
!
! !TO DO: for GSIsa export state contains analysis fields. However it can 
!         also contain increments. This is controlled in GSI routine 
!         update_guess.f90
!
! !REVISION HISTORY:
!
!   19Mar2007 Todling  GSI no longer uses log(ps)
!   17Apr2007 Todling  Swap vertical on way out
!   20Apr2007 Todling  Added phis(ges_z) to export state
!   05Dec2008 Todling  Update SST before writing out
!   10Mar2009 Todling  Place cwmr into qimr; zero out qltot
!   12Mar2009 Todling  Revamp skin temperature imp/exp
!   20Jan2010 Todling  Fix for dealing w/ delp when writing out ainc
!   21Apr2010 Todling  Add trace gases - eventually own GridComp
!   01Sep2012 Akella   Update TSKIN over ocean when doing NSST (nst_gsi > 2)
!   06Dec2013 M.-J.Kim Split cw into qi/ql
!
!EOPI
!-------------------------------------------------------------------------

! local variables

   character(len=*), parameter    :: &
            IAm='GSI_GridCompCopyInternal2Export_'
   integer kk
   real(r_kind),allocatable, dimension(:,:) :: wrk
   integer,     allocatable, dimension(:,:) :: ls_mask
   integer :: nt,it,i4crtm,istatus
   integer(i_kind) :: irank, ipnt, icw
   real(r_kind),pointer:: tv(:,:,:) =>NULL()
   real(r_kind),pointer::  q(:,:,:) =>NULL()
   real(r_kind),pointer:: ql(:,:,:) =>NULL()
   real(r_kind),pointer:: qi(:,:,:) =>NULL()
   real(r_kind),pointer:: qr(:,:,:) =>NULL()
   real(r_kind),pointer:: qs(:,:,:) =>NULL()
   real(r_kind),pointer:: cw(:,:,:) =>NULL()
   character(len=ESMF_MAXSTR) :: cvar

   if(IamRoot) print *,trim(Iam),': Copy contents of internal state to export, it= ',lit

   it = 1
   do it=1,nfldsig
     if (hrdifsig(it) == hrdifsig_all(lit)) exit
   enddo

! simple statistics

!  If so, save full skin temperature field to work array
!  -----------------------------------------------------
   if ( iwrtinc<=0 ) then
       allocate(wrk(size(dts,1),size(dts,2)))
       wrk = dts                                  ! wrk is BKG skin temp ( also called tskp)

       if( nst_gsi > 2) then
          allocate(ls_mask(size(dts,1),size(dts,2)))
          ls_mask  =  zero
          call GSI_GridCompSwapJI_(ls_mask,  isli   (:,:,  it))

          ! gather bkg skin temp into wrk
          ! there is no reason why gsi-isli mask should agree with GCM fractional mask.
          ! Here tdelp has not passed through geos_nstmod- so we need to check for MAPL_UNDEF below.
          where( (ls_mask == zero) .AND. ((tdelp-MAPL_UNDEF) > 0.1*MAPL_UNDEF))  
             wrk = tdelp - dt_coolp     ! over open ocean, and sea-ice (not 100% ice)
          else where  
             wrk = dts                  ! over land, and 100% sea ice
          end where
          deallocate(ls_mask)
       end if
   endif

   call GSI_GridCompSwapJI_(dts,  dsfct   (:,:,  it))

!  Split cloud water fields when needed
!  -------------------------------------
   call gsi_metguess_get ( 'i4crtm::cw', i4crtm , status )
   if (i4crtm==10) then
      call GSI_BundleGetPointer ( GSI_MetGuess_bundle(it), 'cw', cw, status)
      call GSI_BundleGetPointer ( GSI_MetGuess_bundle(it), 'ql', ql, ier); status=status+ier
      call GSI_BundleGetPointer ( GSI_MetGuess_bundle(it), 'qi', qi, ier); status=status+ier
      call GSI_BundleGetPointer ( GSI_MetGuess_bundle(it), 'tv', tv, ier); status=status+ier
      call GSI_BundleGetPointer ( GSI_MetGuess_bundle(it), 'q' ,  q, ier); status=status+ier
      if (status==0) then
         call cloudwater_split_ (tv,q,cw,ql,qi)
      else
         if(IamRoot) print *,trim(Iam),': failed to split cw into qi/ql, aborting ...'
         VERIFY_(STATUS)
      endif
   endif

!  Now handle extra met-guess
!  --------------------------
   do nt=1,nmguess
      cvar = trim(mguess_gsi(nt))
      call GSI_BundleGetPointer ( GSI_MetGuess_bundle(it), cvar, ipnt, status, irank=irank )
      VERIFY_(STATUS)
      if (irank==2) then
          select case (cvar)
            case('ps')
               call GSI_GridCompSwapJI_(dps,GSI_MetGuess_Bundle(it)%r2(ipnt)%q)
            case('z')
               call GSI_GridCompSwapJI_(dhs,GSI_MetGuess_Bundle(it)%r2(ipnt)%q)
          end select
      endif
      if (irank==3) then
          select case (cvar)
            case('u')
               call GSI_GridCompSwapJI_(du,GSI_MetGuess_Bundle(it)%r3(ipnt)%q)
            case('v')
               call GSI_GridCompSwapJI_(dv,GSI_MetGuess_Bundle(it)%r3(ipnt)%q)
            case('tv')
               call GSI_GridCompSwapJI_(dt,GSI_MetGuess_Bundle(it)%r3(ipnt)%q)
            case('q')
               call GSI_GridCompSwapJI_(dq,GSI_MetGuess_Bundle(it)%r3(ipnt)%q)
            case('oz')
               call GSI_GridCompSwapJI_(doz,GSI_MetGuess_Bundle(it)%r3(ipnt)%q)
!_RT        case('cw')   ! GMAO does not save total cloud condensate anymore
            case('qi')
               call GSI_GridCompSwapJI_(dqimr,GSI_MetGuess_Bundle(it)%r3(ipnt)%q)
            case ('ql')
               CALL GSI_GridCompSwapJI_(dqlmr,GSI_MetGuess_Bundle(it)%r3(ipnt)%q)
            case ('qr')
               CALL GSI_GridCompSwapJI_(dqrmr,GSI_MetGuess_Bundle(it)%r3(ipnt)%q)
            case ('qs')
               CALL GSI_GridCompSwapJI_(dqsmr,GSI_MetGuess_Bundle(it)%r3(ipnt)%q)
          end select
      endif
   enddo

!  Now handle trace gases
!  ----------------------
   idco=.false.
   do nt=1,ntgases
      cvar = trim(tgases_gsi(nt))
      call GSI_BundleGetPointer ( GSI_chemguess_bundle(it), cvar, ipnt, status, irank=irank )
      VERIFY_(STATUS)
      if (irank==2) then
         select case (cvar)
           case('aod')
              call GSI_GridCompSwapJI_(daodp,GSI_chemguess_bundle(it)%r2(ipnt)%q)
         end select
      endif
      if (irank==3) then
         select case (cvar)
           case('co')
              call GSI_GridCompSwapJI_(dcop,GSI_chemguess_bundle(it)%r3(ipnt)%q)
              idco=.true.
           case ('co2')
              call GSI_GridCompSwapJI_(dco2p,GSI_chemguess_bundle(it)%r3(ipnt)%q)
         end select
      endif
   enddo

!  SST needs to be updated before being written out 
!  ------------------------------------------------
   if ( iwrtinc<=0 ) then
      if( nst_gsi > 2) then
          dts =  wrk + dts
      else
          where (dts/=MAPL_UNDEF) dts =  wrk + dts             ! Original way of handling sst [no nst]
      end if
      deallocate(wrk)
   endif
#ifdef SFCverbose
   print *,Iam,':  dts (PE,min,max,sum): ',mype,minval(dts),maxval(dts),sum(dts)
#endif

   call UnScale_Export_()

   if ( iwrtinc>0 ) then
        do k=1,nsig
           kk = nsig-k+1
           ddp(:,:,kk)=(bk5(k+1)-bk5(k))*dps(:,:)
        end do
   else
        do k=1,nsig
           kk = nsig-k+1
           ddp(:,:,kk)=Pa_per_kPa*(ak5(k+1)-ak5(k)) +  &
                                  (bk5(k+1)-bk5(k))*dps(:,:)
        end do
   endif

!  If ak(1)+bk(1)*ps is the surface, then ak(2)+bk(2)*ps would
!  be smaller, then dp would be negtive and should be reset to -dp.

   if(sfcFirst_(Pa_per_kPa*ak5,bk5,maxval(dps),iwrtinc>0)) ddp(:,:,:)=-ddp(:,:,:)

   if(IamRoot) print *,trim(Iam),': Copy contents of internal state to export, it= ',lit

   end subroutine GSI_GridCompCopyInternal2Export_

   subroutine cloudwater_split_(tv,q,cw,ql,qi)
   real(r_kind),dimension(:,:,:),intent(in):: tv,q,cw
   real(r_kind),dimension(:,:,:),intent(inout)::ql,qi
   real(r_kind) :: ges_tsen, work 
   integer(i_kind) :: ii,jj,kk
   real(r_kind) :: dcw
 
    do kk=1,nsig
      do jj=1,lon2
        do ii=1,lat2
           dcw = cw(ii,jj,kk) - (ql(ii,jj,kk)+qi(ii,jj,kk)) 
           ges_tsen = tv(ii,jj,kk)/(one+fv*q(ii,jj,kk))
           work = -r0_05*(ges_tsen-t0c) 
           work = max(zero,work)   
           work = min(one,work)   
           ql(ii,jj,kk) = ql(ii,jj,kk) + dcw*(one-work)
           qi(ii,jj,kk) = qi(ii,jj,kk) + dcw*work 
        enddo
      enddo
    enddo
    if(IamRoot) print *,trim(Iam),': Complete splitting cw into qi/ql'

   end subroutine cloudwater_split_

   subroutine UnScale_Export_()
   character(len=*), parameter :: IAm='GSI_GridComp.UnScale_Export_'
   if(associated(dps)) dps = dps  * Pa_per_kPa
   if(GsiGridType==0) then
      if(associated(dhs)) dhs = grav * dhs
      if(associated(doz)) doz = doz  / PPMV2GpG
      if(idco.and.associated(dcop)) then
         dcop = dcop  / KGpKG2PPBV
      endif
   endif
   end subroutine UnScale_Export_

!-------------------------------------------------------------------------
   subroutine GSI_GridCompSetAnaTime_()
!-------------------------------------------------------------------------
   character(len=*), parameter :: IAm='GSI_GridCompSetAnaTime'
   integer(i_kind)              :: nmin_an
   integer(i_kind),dimension(8) :: ida,jda
   integer(i_kind)              :: iyr,ihourg
   integer(i_kind),dimension(4) :: idate4
   character(len=8)             :: cymd
   character(len=6)             :: chms
   real(r_single)               :: hourg4
   real(r_kind),dimension(5)    :: fha

! start

   if(IamRoot) print *,trim(Iam),': Set GSI analysis time '

   ihourg = MYHOURG
   idate4 = MYIDATE

   iyr=idate4(4)
   if ( iyr>=0 .and. iyr<=99 ) then
      if(iyr>51) then
         iyr=iyr+1900
      else
         iyr=iyr+2000
      end if
   end if
   fha(:)=0.0; ida=0; jda=0
   fha(2)=ihourg    ! relative time interval in hours
   ida(1)=iyr       ! year
   ida(2)=idate4(2) ! month
   ida(3)=idate4(3) ! day
   ida(4)=0         ! time zone
   ida(5)=idate4(1) ! hour
   ! Move date-time forward by nhr_assimilation hours
   call w3movdat(fha,ida,jda)

!  For the given time tag, nhr_assimilation, read input guess
!  data header to generate required grid specifications for the
!  GSI analysis grid data.
                                                                                                                                
   ! iadate(:) is for obsmod.
   iadate(1)=jda(1) ! year
   iadate(2)=jda(2) ! mon
   iadate(3)=jda(3) ! day
   iadate(4)=jda(5) ! hour
   iadate(5)=0      ! minute
   ianldate =jda(1)*1000000+jda(2)*10000+jda(3)*100+jda(5)
                                                                                                                                
!  Determine date and time at start of assimilation window
   ida(:)=0
   jda(:)=0
   fha(:)=0.0
   fha(2)=-real(min_offset/60)
   fha(3)=-real(mod(min_offset,60))
   ida(1:3)=iadate(1:3)
   ida(5:6)=iadate(4:5)
   call w3movdat(fha,ida,jda)

   ibdate(1:5)=(/jda(1),jda(2),jda(3),jda(5),jda(6)/)
   iadatebgn=jda(1)*1000000+jda(2)*10000+jda(3)*100+jda(5)

! Set the analysis time - this is output info...
! w3fs21(NCEP-w3) converts analysis time to minutes relative to a fixed date.
   call w3fs21(ibdate,nmin_an)
   iwinbgn = nmin_an

!  Determine date and time at end of assimilation window
   ida(:)=jda(:)
   jda(:)=0
   fha(:)=0.0
   if ( min_offset == 0 ) then
        fha(2)=0
   else
        fha(2)=nhr_assimilation
   endif
   call w3movdat(fha,ida,jda)
 
   iedate(1:5)=(/jda(1),jda(2),jda(3),jda(5),jda(6)/)
   iadateend=jda(1)*1000000+jda(2)*10000+jda(3)*100+jda(5)

!  Get time offset
   call time_4dvar(ianldate,time_offset)

   if(IamRoot) then ! print out some information.

     write(6,*)trim(IAm),':  Guess date is ',idate4,ihourg
     write(6,*)trim(IAm),':  Analysis date is ',iadate
     write(6,*)trim(IAm),':  Start of assimilation window ',ibdate,iadatebgn
     write(6,*)trim(IAm),':  End   of assimilation window ',iedate,iadateend

     write(cymd,'(i4.4,i2.2,i2.2)',iostat=STATUS) iadate(1),iadate(2),iadate(3)
     VERIFY_(STATUS)
     write(chms,'(i2.2,i2.2,i2.2)',iostat=STATUS) iadate(4),0,0
     VERIFY_(STATUS)
     write(6,*) trim(IAm),':        analysis yyyymmdd:hhmmss = '//cymd//':'//chms
     write(6,*) trim(IAm),':        analysis time in minutes =',nmin_an
     write(6,*) trim(IAm),':  time since start of var window =',iwinbgn
     write(6,*) trim(IAm),':                     time offset =',time_offset
      
   endif

   end subroutine GSI_GridCompSetAnaTime_

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: GSI_GridCompSetObsNames  -- 

! !INTERFACE:

   subroutine GSI_GridCompSetObsNames_(lit)

   use mpeu_util, only: tell
   implicit none
   integer, intent(in) :: lit ! logical index of first guess time level

! !DESCRIPTION: temporary routine to set names of GSI obs data files
!
!EOPI
!-------------------------------------------------------------------------

! local variables
   integer             :: nymd, nhms, i, indx

! obs file name templates defined in GSI_GridComp.rc
! ...are there any more than 9?

   character(len=*), parameter :: IAm='GSI_GridCompSetObsNames_'
   character(len=ESMF_MAXSTR) :: DateStamp
   character(len=ESMF_MAXSTR) :: ObsFileDir
   character(len=ESMF_MAXSTR) :: gsiobfname
   character(len=8)           :: cymd
   character(len=6)           :: chms
   integer(i_kind)            :: nobfiles
   integer(i_kind)            :: OBSfreq, OBSfreq_sc, OBSfreq_mn, OBSfreq_hr

! start
#ifdef VERBOSE
   call tell(Iam,'entered with lit =',lit)
#endif

   call ESMF_GridCompGet( gc, config=CF, RC=STATUS )
   VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute( CF, OBSfreq, label ='OBS_FREQUENCY:', rc = STATUS )
   VERIFY_(STATUS)

#ifdef VERBOSE
   call tell(Iam,'returned from ESMF_ConfiGetAttribute()')
#endif

   OBSfreq_hr = OBSfreq/10000
   OBSfreq_mn = mod(OBSfreq,10000)/100
   OBSfreq_sc = mod(OBSfreq,100)
   OBSfreq_sc = OBSfreq_hr*3600 + OBSFreq_mn*60 + OBSfreq_sc
#ifdef VERBOSE
   call tell(Iam,'OBSfreq_hr =',OBSfreq_hr)
   call tell(Iam,'OBSfreq_mn =',OBSfreq_mn)
   call tell(Iam,'OBSfreq_sc =',OBSfreq_sc)
   call tell(Iam,'OBSfreq_sc =',OBSfreq_sc)
#endif

!_RT#ifdef USING_GEOS_FILE_MODEL
! temporarily we get it from iadate...
   write(cymd,'(i4.4,i2.2,i2.2)',iostat=STATUS) iadate(1),iadate(2),iadate(3)
   VERIFY_(STATUS)
   write(chms,'(i2.2,i2.2,i2.2)',iostat=STATUS) iadate(4),0,0
   VERIFY_(STATUS)
#ifdef VERBOSE
   call tell(Iam,'cymd =',trim(cymd))
   call tell(Iam,'chms =',trim(chms))
#endif
!_RT#else
! This is the way to do it if the GSI clock's current time is the analysis 
! time (should be available during coupling to model):
!_RT   call get_DateStamp ( clock, DateStamp, expid, rc=STATUS )
!_RT   VERIFY_(STATUS)
!_RT#endif


   read(cymd,'(i8.8)',iostat=STATUS) nymd
   VERIFY_(STATUS)
   read(chms,'(i6.6)',iostat=STATUS) nhms
   VERIFY_(STATUS)
#ifdef VERBOSE
   call tell(Iam,'nymd =',nymd)
   call tell(Iam,'nhms =',nhms)
#endif

!  Read in observation files table
!  -------------------------------

      call ESMF_ConfigGetAttribute( CF, ObsFileDir, label ='OBS_FILE_DIRECTORY:', default="Obsloc/", rc=STATUS)
      VERIFY_(STATUS)

      allocate( obscls(ndat,2) )
      call get_obsclass_

!  Allocate obs table array
!  ------------------------
   call tell(Iam,'ndat =',ndat)

   call RC_obstable_ ( .false., nobfiles )
   allocate( obstab(nobfiles,3) )
   obstab(:,:)=""
   call RC_obstable_ ( .true., nobfiles )
#ifdef VERBOSE
   call tell(Iam,'returned from RC_obstable_()')

   call tell(Iam,'ndat_times =',ndat_times)
#endif
   if ( ndat_times > 1 ) then
#ifdef VERBOSE
      call tell(Iam,'ndat_times > 1')
#endif
      call tell(Iam,'observation time yyyymmdd:hhmmss = '//cymd//':'//chms)
      do indx = 1, ndat_times
#ifdef VERBOSE
         call tell(Iam,'indx =',indx)
         call tell(Iam,'nobfiles =',nobfiles)
#endif
         ! print *, trim(IAm),':  observation time yyyymmdd:hhmmss = '//cymd//':'//chms
         do i=1,nobfiles
#ifdef VERBOSE
            call tell(Iam,'iobfiles =',i)
#endif
            call match_obclass_ ( trim(obstab(i,1)), gsiobfname )
            if (trim(gsiobfname)/='null') then
                call setObsFile_( trim(obstab(i,3)), trim(gsiobfname), ObsFileDir, nymd, nhms, indx )
            endif
#ifdef VERBOSE
            call tell(Iam,'returned from setObsFile_(), iobfiles =',i)
#endif
         enddo
         call tick ( nymd, nhms, obsfreq_sc ) ! _RT: now this is cheating ... needs generalization
#ifdef VERBOSE
         call tell(Iam,'returned from tick()',i)
#endif
      enddo
   else
#ifdef VERBOSE
      call tell(Iam,'ndat_times !> 1')
#endif
      call tell(Iam,'observation time yyyymmdd:hhmmss = '//cymd//':'//chms)
      ! print *, trim(IAm),':  observation time yyyymmdd:hhmmss = '//cymd//':'//chms
#ifdef VERBOSE
      call tell(Iam,'nobfiles =',nobfiles)
#endif
      do i=1,nobfiles
#ifdef VERBOSE
         call tell(Iam,'iobfiles =',i)
#endif
         call match_obclass_ ( trim(obstab(i,1)), gsiobfname )
         if (trim(gsiobfname)/='null') then
             call setObsFile_( trim(obstab(i,3)), trim(gsiobfname), ObsFileDir, nymd, nhms )
         endif
#ifdef VERBOSE
         call tell(Iam,'returned from setObsFile_(), iobfiles =',i)
#endif
      enddo
   endif ! < ndat_times >
#ifdef VERBOSE
   call tell(Iam,'endif(ndat_times>1)')
   call tell(Iam,'exiting..')
#endif

   deallocate(obstab)
   deallocate(obscls)
   end subroutine GSI_GridCompSetObsNames_

   subroutine RC_obstable_ ( fillin_table, nobfiles )

!  01Oct2007 Todling initial code
!  17Feb2009 Jing Guo	- removed mpeu dependency using m_inpak90
!			- removed VERIFY_() reference, since this
!			  code does not meet the requirements for
!			  VERIFY_() macro and VERIFY(rc/=0) will
!			  not produce an abort().
!  18Feb2009 Jing Guo	- replaced mpeu dependency using m_strtemplate
!			  with local module mpeu_util.

      use mpeu_util,only: warn,tell
      implicit none

      logical,         intent(in)  :: fillin_table
      integer(i_kind), intent(out) :: nobfiles
 
      character(len=*), parameter :: IAm='GSI_GridComp.RC_obstable_'
      character(len=*), parameter :: myname_   = IAm
      character(len=*), parameter :: tablename = 'observation_files::'

      integer iret, ientry, iamroot
      character(len=ESMF_MAXSTR) :: token

!     Initialize
!     ----------
      iamroot = MAPL_AM_I_ROOT()
      iret    = 0

!     Read table with paired datatypes:
!     --------------------------------
      call ESMF_ConfigFindLabel(CF,trim(tablename),rc=iret)
         if (iret/=0) then
            write(6,'(2a,i5,3a)') myname_, ': ESMF_ConfigFindLabel error, iret=', iret, &
                                           ': trying to read "', trim(tablename),'"'
            call MAPL_Abort()
         end if

      ientry   = 0
      nobfiles = 0

      call getnexttoken1_(CF,token,myname_,iret)
      do while (iret==0)     ! read table entries

         ientry=ientry+1
	    if(ientry>ndat) then
               write(6,'(2a,i5)') myname_, ': number of entries in tabel exceed max allowed ', ndat
               write(6,'(2a,2a)') myname_, ': field #1 = "',trim(token),'"'
               iret = 99
	       call MAPL_Abort()
	    endif

         if(fillin_table) then
	 	! Expect three tokens (ignore second) from this line
   	    obstab(ientry,1)=token                           ! ob-class
	    call getnexttoken2_(CF,obstab(ientry,2),myname_) ! ignore this token
	    call getnexttoken2_(CF,obstab(ientry,3),myname_) ! filename template

	     	! echo the input
             !if (iamroot) write(6,'(5a)') myname_, &
             !                            ': observation data entries: ', &
             !                            trim(obstab(ientry,1)), ' ',   &
             !     			  trim(obstab(ientry,3))
         endif

	 	! try the next token
	 call getnexttoken1_(CF,token,myname_,iret)
      enddo

      nobfiles = ientry

!     Release resource file:
!     ---------------------

   end subroutine RC_obstable_

   subroutine match_obclass_ ( current_obclass, gsiobfname )
   character(len=*),intent(in)  :: current_obclass
   character(len=*),intent(out) :: gsiobfname
   integer(i_kind)  ii
   gsiobfname='null'
   do ii=1,ndat
      if (trim(obscls(ii,2))==trim(current_obclass)) then
         gsiobfname=trim(obscls(ii,1))
         exit
      endif
   enddo
   end subroutine match_obclass_

   subroutine getnexttoken1_(cf,token,myname,rc)
     implicit none
     type(ESMF_Config),intent(inout) :: cf
     character(len=*),intent(out) :: token
     character(len=*),intent(in ) :: myname
     integer         ,intent(out) :: rc
     character(len=*), parameter :: IAm='GSI_GridComp.getnexttoken1_'

     token=""
     call ESMF_ConfigNextLine(cf, rc=rc)
       if(rc/=0) return           ! end-of-file is expected.  No error message is produced.

     call ESMF_ConfigGetAttribute(cf, token, rc=rc)
       if(rc/=0) then
         write(6,'(2a,i5)') myname, ': ESMF_ConfigGetAttribute(field #1) error, rc=', rc
         call MAPL_Abort()
       endif

     if(token=='::') rc=-1	! end-of-table is expected.  No error message is produced
   end subroutine getnexttoken1_

   subroutine getnexttoken2_(cf,token,myname)
     implicit none
     type(ESMF_Config),intent(inout) :: cf
     character(len=*),intent(out) :: token
     character(len=*),intent(in ) :: myname
     character(len=*), parameter :: IAm='GSI_GridComp.getnexttoken2_'

     integer :: rc

     token=""
     call ESMF_ConfigGetAttribute(cf, token, rc=rc)
       if(rc/=0) then
         write(6,'(2a,i5)') myname, ': ESMF_ConfigGetAttribute(field #2) error, rc=', rc
         call MAPL_Abort()
       endif
   end subroutine getnexttoken2_


!-------------------------------------------------------------------------
   subroutine setObsFile_(tmpl,gsiname,ObsFileDir,nymd,nhms,indx)
!-------------------------------------------------------------------------
   use mpeu_util, only: tell,warn
! args
   character(len=*), intent(in)           :: tmpl, gsiname
   character(len=*), intent(in)           :: ObsFileDir
   integer(i_kind),  intent(in)           :: nymd, nhms
   integer(i_kind),  intent(in), optional :: indx
! locals
   integer                               :: rc
   logical                               :: fexist
   character(len=ESMF_MAXSTR)            :: syscmd1,syscmd2
   character(len=ESMF_MAXSTR)            :: obsfile, obs_file
   character(len=*), parameter :: Iam='GSI_GridComp.setObsFile_'

   call StrTemplate ( obsfile, tmpl, nymd=nymd, nhms=nhms, stat=STATUS )
   VERIFY_(STATUS)
   inquire(file=trim(ObsFileDir)//trim(obsfile), exist=fexist)
   syscmd1 = '/bin/rm -f '//trim(gsiname)
   syscmd2 = '/bin/ln -s '//trim(ObsFileDir)//trim(obsfile)//' '//trim(gsiname)
   if(fexist) then
      if(present(indx)) then
         write(syscmd2,'(2a,i2.2)') trim(syscmd2), '.', indx
      endif
      ! RT: commented out due to MPT mpi issue in handling system calls 
      ! RT: this is now controlled at script level
      !call system(trim(syscmd1))
      !call system(trim(syscmd2))
      !call EXECUTE_COMMAND_LINE(trim(syscmd1),WAIT=.TRUE.)
      !call EXECUTE_COMMAND_LINE(trim(syscmd2),WAIT=.TRUE.)
      call tell(Iam,trim(syscmd2))
   else
      call warn(Iam,"observation file does not exist, obsfile =",trim(obsfile))
   end if

   end subroutine setObsFile_

!-------------------------------------------------------------------------
   subroutine get_obsclass_
!-------------------------------------------------------------------------
!  28Sep2013  Todling  GMAO-specific routine mimic of obsmod routine to get
!                      GMAO extra column in gsiparm.nml file w/ obsclass
!  11Mar2013  Sienkiewicz  update 'xdfile' size to match recent 'dfile' size
!                          change for aircraft BC profile file name
!-------------------------------------------------------------------------
   use file_utility, only : get_lun
   use mpeu_util, only: gettablesize
   use mpeu_util, only: gettable
   use mpeu_util, only: getindex
   implicit none

   character(len=*),parameter::myname_=myname//'*get_obsclass_'
   character(len=*),parameter::rcname='gsiparm.anl'
   character(len=*),parameter:: tbname='OBS_INPUT::'
   integer(i_kind) luin,i,ii,nrows,ntot
   character(len=ESMF_MAXSTR),allocatable,dimension(:):: utable
   character(10) xdtype,xditype,xdplat
   character(15) xdfile
   character(20) xdsis
   character(32) xcolumn
   real(r_kind)  xdval
   integer(i_kind) xdsfcalc,xdthin

! load file
  luin=get_lun()
  open(luin,file=trim(rcname),form='formatted')

  call gettablesize(tbname,luin,ntot,nrows)
  if(nrows>ndat) then
     call tell(myname_,"inconsistent obs table")
     call MAPL_Abort()
  endif

! Get contents of table
  allocate(utable(nrows))
  call gettable(tbname,luin,ntot,nrows,utable)

! release file unit
  close(luin)

! Retrieve each token of interest from table and define
! variables participating in state vector
   do ii=1,nrows
      read(utable(ii),*) xdfile,   & ! local file name from which to read observatinal data
                         xdtype,   & ! character string identifying type of observatio
                         xdplat,   & ! currently contains satellite id (no meaning for non-sat data)
                         xdsis,    & ! sensor/instrument/satellite identifier for info files
                         xdval,    & ! 
                         xdthin,   & ! thinning flag (1=thinning on; otherwise off)
                         xdsfcalc, & ! use orig bilinear FOV surface calculation (routine deter_sfc)
                         xcolumn     ! use orig bilinear FOV surface calculation (routine deter_sfc)
      obscls(ii,1)=trim(xdfile)
      obscls(ii,2)=trim(xcolumn)
   enddo

   deallocate(utable)

   end subroutine get_obsclass_

   function sfcFirst_(ak,bk,ps,amIinc)
     implicit none
     real(r_kind),dimension(:),intent(in) :: ak,bk
     real,intent(in) :: ps
     logical :: sfcFirst_
     logical :: amIinc
     integer :: m
     m=size(ak)
     if (amIinc) then
             !   p_sfc=bk(1) > p_top=bk(m)
       sfcFirst_= bk(1)      > bk(m)    ! simplified test when dealing w/ increment
     else
             !   p_sfc=pres(1)   >  p_top=pres(m)
       sfcFirst_= ak(1)+bk(1)*ps > ak(m)+bk(m)*ps
     endif
   end function sfcFirst_


   end subroutine Run

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: Finalize -- Finalizes the GSI gridded component

! !INTERFACE:

  subroutine Finalize ( gc, import, export, clock, rc )

!
! !USES:
!
  use gsi_4dcouplermod, only: gsi_4dcoupler_final_traj
  implicit NONE

! !ARGUMENTS:

   type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
   type(ESMF_State),    intent(INOUT) :: import ! Import state
   type(ESMF_State),    intent(INOUT) :: export ! Export state
   type(ESMF_Clock),    intent(INOUT) :: clock  ! The clock
   integer, optional,   intent(  OUT) :: rc     ! Error code:

! !DESCRIPTION:
!
! !TO DO:
!  1. revisit gsi_4dcoupler_final_traj
!
! !REVISION HISTORY:
!
!   14Apr2007 Todling Handling ifile arrays as allocatables 
!
!EOPI
!-------------------------------------------------------------------------

   integer                           :: STATUS    ! error code STATUS
   logical                           :: IamRoot    
   character(len=*), parameter :: IAm='Gsi_GridComp.Finalize'

! start

   IamRoot = MAPL_AM_I_ROOT()
   if(IamRoot.and.verbose) print *, Iam,": Start ",trim(Iam)

!_RT  if(l4dvar) then
!_RT    call gsi_4dcoupler_final_traj(rc=STATUS)
!_RT    VERIFY_(STATUS)
!_RT  endif

   call ESMF_ClockDestroy(GSI_aClock,rc=STATUS)
   VERIFY_(STATUS)

! GSI_AgcmPertGrid, is a grid defined for the AgcmPert grid in 4DVAR
   call ESMF_GridDestroy(GSI_AgcmPertGrid,  RC=STATUS)
   VERIFY_(STATUS)

! Call GEOS Generic Finalize

   call MAPL_GenericFinalize ( gc, import, export, clock,  RC=STATUS)
   VERIFY_(STATUS)

! Deallocate memory for GSI variables

   call GSI_GridCompDealloc_()

   deallocate(ifilesig, ifilesfc, ifilenst, stat=STATUS)
   VERIFY_(STATUS)
   deallocate(hrdifsig, hrdifsfc, hrdifnst, stat=STATUS)
   VERIFY_(STATUS)
   deallocate(hrdifsig_all, hrdifsfc_all, hrdifnst_all, stat=STATUS)
   VERIFY_(STATUS)
   
   if(IamRoot) print *, Iam,": End ",trim(Iam)

   RETURN_(ESMF_SUCCESS)

   CONTAINS

!-------------------------------------------------------------------------
   subroutine GSI_GridCompDealloc_()
!-------------------------------------------------------------------------
   character(len=*), parameter :: IAm='GSI_GridCompDealloc'

   integer(i_kind) mype
   type(ESMF_VM)  :: VM

   call ESMF_VMGetCurrent(vm=vm, rc=STATUS)
   VERIFY_(STATUS)
   call ESMF_VMGet(vm, localPet=mype, rc=STATUS)
   VERIFY_(STATUS)

!  Deallocate memory for guess files for trace gases
!  ------------------------------------------------
   call destroy_chemges_grids(status)
   VERIFY_(STATUS)

   call destroy_metguess_grids(mype,status)
   VERIFY_(STATUS)

   deallocate( isli     ,&
               fact10   ,&
               sfct     ,&
               dsfct    ,&
               sno      ,&
               veg_type ,&
               veg_frac ,&
               soil_type,&
               soil_temp,&
               soil_moi ,&
               sfc_rough,&
               stat=STATUS)
   VERIFY_(STATUS)

   call destroy_grid_vars

   end subroutine GSI_GridCompDealloc_

   end subroutine Finalize	

!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: GSI_GridCompSetupSpecs -- Sets import/export specs

! !INTERFACE:

   subroutine GSI_GridCompSetupSpecs (GC, opthw, rc)

!
! !USES:
!
  implicit NONE

! !ARGUMENTS:

   type(ESMF_GridComp), intent(INOUT) :: gc    ! gridded component
   integer, optional, intent(IN)      :: opthw ! optional halowidth
   integer, optional, intent(OUT)     :: rc    ! return code

! !DESCRIPTION: Sets the import and export specs. The GSI expects (imports) 
!       the following fields:
!
!       UPA: Hs, Ps, Tsk, u, v, Tv, vor, div, q, oz, qi, ql
!       SFC: F10m, Sno, qg, Tg, SLI, Vfrac, Vtyp, Styp, Sfc_rough
! 
! GEOS does not export vor, div - calculated from u,v by GSI
!                      F10m     - calculated from u,v,u10,v10 (exported by GEOS)
!                      Vfrac, Vtyp, Styp - read from a file 
!
!  Whenever possible we have used the same SHORT_NAME between corresponding
!  GEOS-5 and GSI fields. In some cases, as when coupling, it clearly changes the
!  nature of the field, we have used GSI specific names, to emphasize 
!  the change.
!
! !REVISION HISTORY:
!
!  22Apr2007 Todling Properly naming u/v/tv fields
!  24Apr2007 Todling Added delp to export state; tskin to ts from bkg.eta
!  10Mar2009 Todling Add fraction (land/ice/etc) and replace qc w/ qi/qltot
!
!EOPI
!-------------------------------------------------------------------------

   integer                               :: STATUS, local_hw
   character(len=*), parameter :: IAm='GSI_GridCompSetupSpecs'

   integer(i_kind) ii, iii
   logical iamroot
   character(len=ESMF_MAXSTR) :: cvar

!  Declare import 2d-fields
!  ------------------------
   integer, parameter :: nin2d=14
   character(len=16), parameter :: insname2d(nin2d) = (/  &
                                   'phis            ',    &
                                   'ps              ',    &
                                   'ts              ',    &
!                                  'NCEP_VEGFRAC    ',    &
!                                  'NCEP_VEGTYPE    ',    &
!                                  'NCEP_SOITYPE    ',    &
                                   'U10M            ',    &
                                   'V10M            ',    &
                                   'SNOWDP          ',    &
                                   'GWETTOP         ',    &
                                   'TSOIL1          ',    &
                                   'Z0M             ',    &
                                   'frland          ',    &
                                   'frlandice       ',    &
                                   'frlake          ',    &
                                   'frocean         ',    &
                                   'frseaice        '    /)
   character(len=32), parameter :: inlname2d(nin2d) = (/  &
                      'geopotential height             ', &
                      'surface pressure                ', &
                      'skin temperature                ', &
!                     'NCEP(CRTM-like) veg fraction    ', &
!                     'NCEP(CRTM-like) veg type        ', &
!                     'NCEP(CRTM-like) soil type       ', &
                      'u 10 m-wind                     ', &
                      'v 10 m-wind                     ', &
                      'snow depth                      ', &
                      'surface soil wetness            ', &
                      'soil temperature                ', &
                      'roughness length                ', &
                      'fraction of land                ', &
                      'fraction of land ice            ', &
                      'fraction of lake                ', &
                      'fraction of ocean               ', &
                      'fraction of sea ice             ' /)
   character(len=16), parameter :: inunits2d(nin2d) = (/  &
                                   'm**2/s**2       ',    &
                                   'Pa              ',    &
                                   'K               ',    &
!                                  '1               ',    &
!                                  '1               ',    &
!                                  '1               ',    &
                                   'm/s             ',    &
                                   'm/s             ',    &
                                   'm               ',    &
                                   '1               ',    &
                                   'K               ',    &
                                   'm               ',    &
                                   '1               ',    &
                                   '1               ',    &
                                   '1               ',    &
                                   '1               ',    &
                                   '1               '    /)

!  Declare import 3d-fields
!  ------------------------
   integer, parameter :: nin3d=5
   character(len=16), parameter :: insname3d(nin3d) = (/  &
                                   'u               ',    &
                                   'v               ',    &
                                   'tv              ',    &
                                   'sphu            ',    &
                                   'ozone           '    /)
   character(len=40), parameter :: inlname3d(nin3d) = (/  &
                      'eastward wind                           ', &
                      'northward wind                          ', &
                      'air virtual temperature                 ', &
                      'specific humidity                       ', &
                      'ozone mass mixing ratio                 ' /)
   character(len=16), parameter :: inunits3d(nin3d) = (/  &
                                   'm/s             ',    &
                                   'm/s             ',    &
                                   'K               ',    &
                                   'g/g             ',    &
                                   'ppmv            '    /)

!  Declare import 2d-fields (extra met-guess)
!  ------------------------
   integer, parameter :: nin2dx=5
   character(len=16), parameter :: insname2dx(nin2dx) = (/  &
                                   'DCOOL              ',      &
                                   'DWARM              ',      &
                                   'TDROP              ',      &
                                   'TDEL               ',      &
                                   'TS_FOUND           ' /)
   character(len=40), parameter :: inlname2dx(nin2dx) = (/  &
                      'depth of cool layer                     ', &
                      'depth of warm layer                     ', &
                      'temperature drop cool layer             ', &
                      'temperature at top of warm layer        ', &
                      'foundation SST                          ' /)
   character(len=16), parameter :: inunits2dx(nin2dx) = (/&
                                   'm               ',    &
                                   'm               ',    &
                                   'K               ',    &
                                   'K               ',    &
                                   'K               '    /)

!  Declare import 3d-fields (extra met-guess)
!  ------------------------
   integer, parameter :: nin3dx=5
   character(len=16), parameter :: insname3dx(nin3dx) = (/  &
                                   'qitot           ',      &
                                   'qltot           ',      &
                                   'qrtot           ',      &
                                   'qstot           ',      &
                                   'cloud           '      /)
   character(len=40), parameter :: inlname3dx(nin3dx) = (/  &
                      'mass fraction of cloud ice water        ', &
                      'mass fraction of cloud liquid water     ', &
                      'mass fraction of rain             r     ', &
                      'mass fraction of snow             r     ', &
                      'cloud fraction for radiation            ' /)
   character(len=16), parameter :: inunits3dx(nin3dx) = (/&
                                   '1               ',    &
                                   '1               ',    &
                                   '1               ',    &
                                   '1               ',    &
                                   '1               '    /)

!  Declare import 3d-fields for trace gases - this needs 
!  to come from gsi_chemguess_mod
!  -----------------------------------------------------
   integer, parameter :: nin2dg=1
   character(len=16), parameter :: insname2dg(nin2dg) = (/ &
                                   'AOD             '/)
   character(len=32), parameter :: inlname2dg(nin2dg) = (/ &
                      'aerosol optical depth           '/)
   character(len=16), parameter :: inunits2dg(nin2dg) = (/ &
                                   '1               '/)

!  Declare import 3d-fields for trace gases - this needs 
!  to come from gsi_chemguess_mod
!  -----------------------------------------------------
   integer, parameter :: nin3dg=20
   character(len=16), parameter :: insname3dg(nin3dg) = (/ &
                                   'CO              ',     &
                                   'CO2             ',     &
                                   'DU001           ',     &
                                   'DU002           ',     &
                                   'DU003           ',     &
                                   'DU004           ',     &
                                   'DU005           ',     &
                                   'SS001           ',     &
                                   'SS002           ',     &
                                   'SS003           ',     &
                                   'SS004           ',     &
                                   'SS005           ',     &
                                   'DMS             ',     &
                                   'SO2             ',     &
                                   'SO4             ',     &
                                   'MSA             ',     &
                                   'BCPHOBIC        ',     &
                                   'BCPHILIC        ',     &
                                   'OCPHOBIC        ',     &
                                   'OCPHILIC        '/)
   character(len=32), parameter :: inlname3dg(nin3dg) = (/ &
                      'carbon monoxide                 ',  &
                      'carbon dioxide                  ',  &
                      'dust                            ',  &
                      'dust                            ',  &
                      'dust                            ',  &
                      'dust                            ',  &
                      'dust                            ',  &
                      'sea salt                        ',  &
                      'sea salt                        ',  &
                      'sea salt                        ',  &
                      'sea salt                        ',  &
                      'sea salt                        ',  &
                      'sulfate                         ',  &
                      'sulfate                         ',  &
                      'sulfate                         ',  &
                      'sulfate                         ',  &
                      'dry black carbon                ',  &
                      'wet black carbon                ',  &
                      'dry organic carbon              ',  &
                      'wet organic carbon              '/)
   character(len=16), parameter :: inunits3dg(nin3dg) = (/ &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             ',     &
                                   'g/g             '     /)


!  Declare export 2d-fields
!  ------------------------
   integer, parameter :: nex2d=8
   character(len=16), parameter :: exsname2d(nex2d) = (/  &
                                   'phis            ',    &
                                   'ps              ',    &
                                   'ts              ',    &
                                   'frland          ',    &
                                   'frlandice       ',    &
                                   'frlake          ',    &
                                   'frocean         ',    &
                                   'frseaice        '    /)
   character(len=32), parameter :: exlname2d(nex2d) = (/  &
                      'geopotential height             ', &
                      'surface pressure inc            ', &
                      'skin temperature inc            ', &
                      'fraction of land                ', &
                      'fraction of land ice            ', &
                      'fraction of lake                ', &
                      'fraction of ocean               ', &
                      'fraction of sea ice             ' /)
   character(len=16), parameter :: exunits2d(nex2d) = (/  &
                                   'm**2/s**2       ',    &
                                   'Pa              ',    &
                                   'K               ',    &
                                   '1               ',    &
                                   '1               ',    &
                                   '1               ',    &
                                   '1               ',    &
                                   '1               '    /)

!  Declare export 3d-fields
!  ------------------------
   integer, parameter :: nex3d=6
   character(len=16), parameter :: exsname3d(nex3d) = (/  &
                                   'u               ',    &
                                   'v               ',    &
                                   'tv              ',    &
                                   'delp            ',    &
                                   'sphu            ',    &
                                   'ozone           '    /)
   character(len=40), parameter :: exlname3d(nex3d) = (/  &
                      'u-wind inc                              ', &
                      'v-wind inc                              ', &
                      'virtual temperature inc                 ', &
                      'delta pressure inc                      ', &
                      'specific humidity inc                   ', &
                      'ozone mass mixing ratio inc             ' /)
   character(len=16), parameter :: exunits3d(nex3d) = (/  &
                                   'm/s             ',    &
                                   'm/s             ',    &
                                   'K               ',    &
                                   'Pa              ',    &
                                   'g/g             ',    &
                                   'g/g             '    /)

!  Declare export 3d-fields (extra met-guess)
!  ------------------------
   integer, parameter :: nex3dx=4
   character(len=16), parameter :: exsname3dx(nex3dx) = (/  &
                                   'qitot           ',      &
                                   'qltot           ',      &
                                   'qrtot           ',      &
                                   'qstot           '      /)
   character(len=40), parameter :: exlname3dx(nex3dx) = (/  &
                      'mass fraction of cloud ice water inc    ', &
                      'mass fraction of cloud liquid water inc ', &
                      'mass fraction of rain inc               ', &
                      'mass fraction of snow inc               ' /)
   character(len=16), parameter :: exunits3dx(nex3dx) = (/  &
                                   '1               ',      &
                                   '1               ',      &
                                   '1               ',      &
                                   '1               '      /)

!  Declare export 2d-fields for trace gases - this needs 
!  to come from gsi_chemguess_mod
!  -----------------------------------------------------
   integer, parameter :: nex2dg=1
   character(len=16), parameter :: exsname2dg(nex2dg) = (/ &
                                   'AOD             '/)
   character(len=32), parameter :: exlname2dg(nex2dg) = (/ &
                      'aerosol optical depth           '/)
   character(len=16), parameter :: exunits2dg(nex2dg) = (/ &
                                   '1               '     /)

!  Declare export 3d-fields for trace gases - this needs 
!  to come from gsi_chemguess_mod
!  -----------------------------------------------------
   integer, parameter :: nex3dg=2
   character(len=16), parameter :: exsname3dg(nex3dg) = (/ &
                                   'CO              ',     &
                                   'CO2             '/)
   character(len=32), parameter :: exlname3dg(nex3dg) = (/ &
                      'carbon monoxide inc             ',  &
                      'carbon dioxide inc              '/)
   character(len=16), parameter :: exunits3dg(nex3dg) = (/ &
                                   'g/g             ',     &
                                   'g/g             '     /)


! Begin

   if(present(opthw)) then
      local_hw=opthw
   else
      local_hw=HW
   end if

!  Inquire about extra-guess and chemistry fields

   iamroot = MAPL_AM_I_ROOT()

   call gsi_metguess_init(iamroot=iamroot)  ! should not be done here
   call gsi_metguess_get('dim',nmguess,status)
   VERIFY_(STATUS)
   if (nmguess>0.and.(.not.allocated(mguess_gsi))) then
       allocate(mguess_gsi(nmguess))
       call gsi_metguess_get('gsinames',mguess_gsi,status)
   endif
   if (nmguess>0.and.(.not.allocated(mguess_usr))) then
       allocate(mguess_usr(nmguess))
       call gsi_metguess_get('usrnames',mguess_usr,status)
   endif

   call gsi_chemguess_init(iamroot=iamroot)  ! should not be done here
   call gsi_chemguess_get('dim',ntgases,status)
   if (ntgases>0.and.(.not.allocated(tgases_gsi))) then
       allocate(tgases_gsi(ntgases))
       call gsi_chemguess_get('gsinames',tgases_gsi,status)
   endif
   if (ntgases>0.and.(.not.allocated(tgases_usr))) then
       allocate(tgases_usr(ntgases))
       call gsi_chemguess_get('usrnames',tgases_usr,status)
   endif

! Get variables from chem
! call gsi_chemguess_get('gsinames')
! call gsi_chemguess_get('longnames')
! call gsi_chemguess_get('units')

! Imports


! Begin

   if(present(opthw)) then
      local_hw=opthw
   else
      local_hw=HW
   end if

! Imports

    do ii = 1, nin2d 
       call MAPL_AddImportSpec(GC,        &
         SHORT_NAME= trim(insname2d(ii)), &
         LONG_NAME = trim(inlname2d(ii)), &
         UNITS     = trim(inunits2d(ii)), &
         DIMS      = MAPL_DimsHorzOnly,   &
         VLOCATION = MAPL_VLocationNone,  &
         HALOWIDTH = local_hw,            &
         RC=STATUS  ); VERIFY_(STATUS)
    enddo
    do ii = 1, nin3d 
     call MAPL_AddImportSpec(GC,           &
         SHORT_NAME= trim(insname3d(ii)),  &
         LONG_NAME = trim(inlname3d(ii)),  &
         UNITS     = trim(inunits3d(ii)),  &
         DIMS      = MAPL_DimsHorzVert,    &
         VLOCATION = MAPL_VLocationCenter, &
         HALOWIDTH = local_hw,             &
         RC=STATUS  ); VERIFY_(STATUS)
    enddo
    do iii = 1, nmguess 
       cvar = trim(mguess_usr(iii))
       do ii = 1, nin2dx 
          if ( trim(lowercase(cvar))==trim(lowercase(insname2dx(ii))) ) then
               call MAPL_AddImportSpec(GC,           &
                   SHORT_NAME= trim(insname2dx(ii)), &
                   LONG_NAME = trim(inlname2dx(ii)), &
                   UNITS     = trim(inunits2dx(ii)), &
                   DIMS      = MAPL_DimsHorzOnly,    &
                   VLOCATION = MAPL_VLocationNone,   &
                   HALOWIDTH = local_hw,             &
                   RC=STATUS  ); VERIFY_(STATUS)
          endif
       enddo
       do ii = 1, nin3dx 
          if ( trim(lowercase(cvar))==trim(lowercase(insname3dx(ii))) ) then
               call MAPL_AddImportSpec(GC,           &
                   SHORT_NAME= trim(insname3dx(ii)), &
                   LONG_NAME = trim(inlname3dx(ii)), &
                   UNITS     = trim(inunits3dx(ii)), &
                   DIMS      = MAPL_DimsHorzVert,    &
                   VLOCATION = MAPL_VLocationCenter, &
                   HALOWIDTH = local_hw,             &
                   RC=STATUS  ); VERIFY_(STATUS)
          endif
       enddo
    enddo
    do iii = 1, ntgases
       cvar = trim(tgases_usr(iii))
       do ii = 1, nin2dg 
          if ( trim(lowercase(cvar))==trim(lowercase(insname2dg(ii))) ) then
               call MAPL_AddImportSpec(GC,           &
                  SHORT_NAME= trim(insname2dg(ii)), &
                  LONG_NAME = trim(inlname2dg(ii)), &
                  UNITS     = trim(inunits2dg(ii)), &
                  DIMS      = MAPL_DimsHorzOnly,    &
                  VLOCATION = MAPL_VLocationNone,   &
                  HALOWIDTH = local_hw,             &
                  RC=STATUS  ); VERIFY_(STATUS)
          endif
       enddo
       do ii = 1, nin3dg 
          if ( trim(lowercase(cvar))==trim(lowercase(insname3dg(ii))) ) then
               call MAPL_AddImportSpec(GC,           &
                  SHORT_NAME= trim(insname3dg(ii)), &
                  LONG_NAME = trim(inlname3dg(ii)), &
                  UNITS     = trim(inunits3dg(ii)), &
                  DIMS      = MAPL_DimsHorzVert,    &
                  VLOCATION = MAPL_VLocationCenter, &
                  HALOWIDTH = local_hw,             &
                  RC=STATUS  ); VERIFY_(STATUS)
          endif
       enddo
    enddo

! Exports

    do ii = 1, nex2d
       call MAPL_AddExportSpec(GC,         &
         SHORT_NAME= trim(exsname2d(ii)),  &
         LONG_NAME = trim(exlname2d(ii)),  &
         UNITS     = trim(exunits2d(ii)),  &
         DIMS      = MAPL_DimsHorzOnly,    &
         VLOCATION = MAPL_VLocationNone,   &
         HALOWIDTH = local_hw,             &
         RC=STATUS  ); VERIFY_(STATUS)
    enddo
    do ii = 1, nex3d
       call MAPL_AddExportSpec(GC,         &
         SHORT_NAME= trim(exsname3d(ii)),  &
         LONG_NAME = trim(exlname3d(ii)),  &
         UNITS     = trim(exunits3d(ii)),  &
         DIMS      = MAPL_DimsHorzVert,    &
         VLOCATION = MAPL_VLocationCenter, &
         HALOWIDTH = local_hw,             &
         RC=STATUS  );  VERIFY_(STATUS)
    enddo
    do iii = 1, nmguess
       cvar = trim(mguess_usr(iii))
       do ii = 1, nex3dx
          if ( trim(lowercase(cvar))==trim(lowercase(exsname3dx(ii))) ) then
               call MAPL_AddExportSpec(GC,         &
                 SHORT_NAME= trim(exsname3dx(ii)),  &
                 LONG_NAME = trim(exlname3dx(ii)),  &
                 UNITS     = trim(exunits3dx(ii)),  &
                 DIMS      = MAPL_DimsHorzVert,    &
                 VLOCATION = MAPL_VLocationCenter, &
                 HALOWIDTH = local_hw,             &
                 RC=STATUS  );  VERIFY_(STATUS)
          endif
       enddo
    enddo
    do iii = 1, ntgases
       cvar = trim(tgases_usr(iii))
       do ii = 1, nex2dg
          if ( trim(lowercase(cvar))==trim(lowercase(exsname2dg(ii))) ) then
               call MAPL_AddExportSpec(GC,         &
                 SHORT_NAME= trim(exsname2dg(ii)), &
                 LONG_NAME = trim(exlname2dg(ii)), &
                 UNITS     = trim(exunits2dg(ii)), &
                 DIMS      = MAPL_DimsHorzOnly,    &
                 VLOCATION = MAPL_VLocationNone,   &
                 HALOWIDTH = local_hw,             &
                 RC=STATUS  );  VERIFY_(STATUS)
          endif
       enddo
       do ii = 1, nex3dg
          if ( trim(lowercase(cvar))==trim(lowercase(exsname3dg(ii))) ) then
               call MAPL_AddExportSpec(GC,         &
                 SHORT_NAME= trim(exsname3dg(ii)), &
                 LONG_NAME = trim(exlname3dg(ii)), &
                 UNITS     = trim(exunits3dg(ii)), &
                 DIMS      = MAPL_DimsHorzVert,    &
                 VLOCATION = MAPL_VLocationCenter, &
                 HALOWIDTH = local_hw,             &
                 RC=STATUS  );  VERIFY_(STATUS)
          endif
       enddo
    enddo

   end subroutine GSI_GridCompSetupSpecs

!-------------------------------------------------------------------------
   subroutine GSI_GridCompSetAlarms (GENSTATE, GC, cf, clock)
!-------------------------------------------------------------------------
   use mpeu_util,only : tell,perr,die
   type(MAPL_MetaComp), pointer     :: GENSTATE  ! GEOS Generic state
   type(ESMF_GridComp)              :: GC
   type(ESMF_Config)                :: cf
   type(ESMF_Clock)                 :: clock

!  locals
!  ------
   character(len=*), parameter :: IAm='GSI_GridCompSetAlarms'

   type(ESMF_Alarm), pointer        :: ALARM(:)  ! this component's alarms
   type(ESMF_Time)                  :: CurrTime
   type(ESMF_Time)                  :: AlarmTime
   type(ESMF_Time)                  :: AlarmTime0
   type(ESMF_Time)                  :: AlaTime
   type(ESMF_Time)                  :: AnaTime
   type(ESMF_Time)                  :: RefTime
   type(ESMF_Time)                  :: FcsTime
   type(ESMF_TimeInterval)          :: FrqBKG
   type(ESMF_TimeInterval)          :: FrqANA
   type(ESMF_TimeInterval)          :: AnaOST
   type(ESMF_TimeInterval)          :: SmallTimeStep
   integer                          :: atime, i
   integer                          :: REF_TIME(6), RREF_DATE, RREF_TIME
   integer                          :: FCS_TIME(6), nymd_fcst,nhms_fcst
   integer                          :: status, rc
   integer                          :: yy, mm, dd, hh, mn, sec
   logical                          :: IamRoot
   character(len=ESMF_MAXSTR)       :: aname
   character(len=2)                 :: idxtime
   type(ESMF_Alarm), allocatable    :: AlarmList(:)
   type(ESMF_Time),  allocatable    :: AlarmRingTime(:)
   logical,          allocatable    :: ringingState(:)
   integer                          :: nalarms

!  start
!  -----
   IamRoot = MAPL_AM_I_ROOT()

   if(IamRoot) call myclock_show_(clock,IAm,'initial-clock')

   call MAPL_GetObjectFromGC ( GC, GENSTATE, RC=STATUS ); VERIFY_(STATUS)

!ALT: since we are manipulating the clock we need to do some bookeeping
!     first we will save the state of all alarms already in the clock
!     and after the clock manipulation we will restore them

   ! save alarms' states
   call ESMF_ClockGet(clock, alarmCount = nalarms, rc=status)
   VERIFY_(STATUS)
   allocate (alarmList(nalarms), alarmRingTime(nalarms), ringingState(nalarms), stat = status)
   VERIFY_(STATUS)
   call ESMF_ClockGetAlarmList(clock, alarmListFlag=ESMF_ALARMLIST_ALL, &
        alarmList=alarmList, alarmCount = nalarms, rc=status)
   VERIFY_(STATUS)
   DO I = 1, nalarms
      call ESMF_AlarmGet(alarmList(I), ringTime=alarmRingTime(I), ringing=ringingState(I), rc=status)
      VERIFY_(STATUS)
   END DO

   call ESMF_ClockAdvance(clock, rc=status)
   VERIFY_(STATUS)

   if(IamRoot) call myclock_show_(clock,IAm,'advanced-clock')

!  Define analysis time from reference time in rc file
!  ---------------------------------------------------
   CALL ESMF_ConfigGetAttribute( cf, RREF_DATE, label = 'RECORD_REF_DATE:', rc=status )
     VERIFY_(status)
   CALL ESMF_ConfigGetAttribute( cf, RREF_TIME, label = 'RECORD_REF_TIME:', rc=status )
     VERIFY_(status)

   REF_TIME(1) =     RREF_DATE/10000
   REF_TIME(2) = mod(RREF_DATE,10000)/100
   REF_TIME(3) = mod(RREF_DATE,100)
   REF_TIME(4) =     RREF_TIME/10000
   REF_TIME(5) = mod(RREF_TIME,10000)/100
   REF_TIME(6) = mod(RREF_TIME,100)

!  in case necessary, initial time tau-days from now
!  ------------------------------------------------- 
   if (tau_fcst>0) then
      nymd_fcst = RREF_DATE
      nhms_fcst = RREF_TIME
      call tick (nymd_fcst,nhms_fcst,tau_fcst*3600)
      FCS_TIME(1) =     nymd_fcst/10000
      FCS_TIME(2) = mod(nymd_fcst,10000)/100
      FCS_TIME(3) = mod(nymd_fcst,100)
      FCS_TIME(4) =     nhms_fcst/10000
      FCS_TIME(5) = mod(nhms_fcst,10000)/100
      FCS_TIME(6) = mod(nhms_fcst,100)
      call ESMF_TimeSet(  GSI_FcsTime, YY =  FCS_TIME(1), &
                                       MM =  FCS_TIME(2), &
                                       DD =  FCS_TIME(3), &
                                       H  =  FCS_TIME(4), &
                                       M  =  FCS_TIME(5), &
                                       S  =  FCS_TIME(6), rc=status ); VERIFY_(STATUS)
   endif

!  initialize current time
!  -----------------------
   call ESMF_TimeSet(  RefTime, YY =  REF_TIME(1), &
                                MM =  REF_TIME(2), &
                                DD =  REF_TIME(3), &
                                H  =  REF_TIME(4), &
                                M  =  REF_TIME(5), &
                                S  =  REF_TIME(6), rc=status ); VERIFY_(STATUS)

   GSI_RefTime = RefTime

!  Create num. FG time levels + 1 alarms centered on the analysis time
!  -------------------------------------------------------------------
   call ESMF_ClockGet (clock, currTime=currTime, rc=STATUS)
     VERIFY_(STATUS)
   call ESMF_TimeGet(currTime, yy=YY, mm=MM, dd=DD, h=HH, m=MN, s=SEC, rc=rc)
     if(IamRoot) write(6,'(a,1x,i4.4,5(a,i2.2))') trim(Iam)//" The current TIME is ", &
                                       YY, "/", MM, "/", DD, " ", HH, ":", MN, ":", SEC

   call ESMF_TimeSet(AlarmTime, yy=yy, mm=mm, dd=dd, h=hh, rc=STATUS)
     VERIFY_(STATUS)

   call ESMF_TimeIntervalSet(FrqBKG, s=BKGfreq_sc, rc=STATUS)
     VERIFY_(STATUS)

   call ESMF_TimeIntervalSet(FrqANA, h=ANAfreq_sc, rc=STATUS)
     VERIFY_(STATUS)

   call ESMF_TimeIntervalSet(AnaOST, h=min_offset/60,m=mod(min_offset,60), rc=STATUS)
     VERIFY_(STATUS)


!
!  n+1 alarms (n=number of BKG time levels) are set around the analysis time:
!
! n=1                      *
! n=3                *     *     *
! n=5          *     *     *     *     *
! etc                      |
!                        ana-t
!
! n+1 th alarm rings when all BKG time levels have been collected
!

   allocate(ALARM(nfldsig_all), stat=STATUS); VERIFY_(STATUS)

   AnaTime    = AlarmTime + AnaOST  ! alarm to indicate analysis time
   AlarmTime0 = AlarmTime
   GSI_Time0  = AlarmTime0
   do atime = 1,nfldsig_all-1
      write(idxtime,'(i2)',iostat=STATUS) atime;                                     VERIFY_(STATUS)
      aname = "some-bkg " // idxtime
      ALARM(atime) = ESMF_AlarmCreate(name=aname, clock=clock, ringTime=alarmTime0, rc=rc);     VERIFY_(STATUS)
      call MAPL_StateAlarmAdd(GENSTATE,ALARM(atime),RC=status);                      VERIFY_(STATUS)
      call ESMF_AlarmGet(ALARM(atime), ringTime=AlaTime, rc=status);                 VERIFY_(STATUS)
      call ESMF_TimeGet(AlaTime, yy=YY, mm=MM, dd=DD, h=HH, m=MN, s=SEC, rc=status); VERIFY_(STATUS)
      if(IamRoot) print *, Iam,": ",trim(aname)," background time set to ", YY, "/", MM, "/", DD, " ", HH, ":", MN, ":", SEC
      if( AlaTime == RefTime ) then
          ntguessig_ref = atime
          ntguessfc_ref = atime
          ntguesnst_ref = atime
          if (IamRoot) print *, Iam,": Upper-air updated pointer (ntguessig_ref): ", ntguessig_ref
          if (IamRoot) print *, Iam,": Surface   updated pointer (ntguessfc_ref): ", ntguessfc_ref
          if (IamRoot) print *, Iam,": Surface   updated pointer (ntguesnst_ref): ", ntguesnst_ref
      endif
#ifdef VERBOSE
      call tell(Iam,"if(Alatime==AnaTime)")
#endif
      if( AlaTime == AnaTime ) then
          MYIDATE = (/ HH, MM, DD, YY /)
      endif
#ifdef VERBOSE
      call tell(Iam,"endif(Alatime==AnaTime)")
#endif
      alarmTime0 = alarmTime0 + FrqBKG
   end do
#ifdef VERBOSE
   call tell(Iam,"enddo")
#endif

   	! locate where the slot is available
#ifdef VERBOSE
   call tell(Iam,"nfldsig =",nfldsig)
   call tell(Iam,"ntguessig_ref =",ntguessig_ref)
   call tell(Iam,"allocated(hrdifsig) =",allocated(hrdifsig))
   call tell(Iam,"allocated(hrdifsig_all) =",allocated(hrdifsig_all))
   call tell(Iam,"size(hrdifsig,1) =",size(hrdifsig,1))
   call tell(Iam,"size(hrdifsig_all,1) =",size(hrdifsig_all,1))
   call tell(Iam,"hrdifsig_all(ntguessig_ref) =",hrdifsig_all(ntguessig_ref))
#endif

   ! the last alarm notifies aana that analysis is done

   atime = nfldsig_all
   aname = "last-bkg"
   ALARM(atime) = ESMF_AlarmCreate(name=aname, clock=clock, ringTime=alarmTime0, rc=rc);     VERIFY_(STATUS)
   call MAPL_StateAlarmAdd(GENSTATE,ALARM(atime),RC=status);                      VERIFY_(STATUS)
   call ESMF_AlarmGet(ALARM(atime), ringTime=AlaTime, rc=status);                 VERIFY_(STATUS)
   call ESMF_TimeGet(AlaTime, yy=YY, mm=MM, dd=DD, h=HH, m=MN, s=SEC, rc=status); VERIFY_(STATUS)
   if(IamRoot) print *,Iam,": ",trim(aname)," at ana time set to ", YY, "/", MM, "/", DD, " ", HH, ":", MN, ":", SEC
   if( AlaTime == RefTime ) then
       ntguessig_ref = atime
       ntguessfc_ref = atime
       ntguesnst_ref = atime
       if (IamRoot) print *, Iam,": Upper-air updated pointer (ntguessig_ref): ", ntguessig_ref
       if (IamRoot) print *, Iam,": Surface   updated pointer (ntguessig_ref): ", ntguessfc_ref
       if (IamRoot) print *, Iam,": Surface   updated pointer (ntguesnst_ref): ", ntguesnst_ref

   	! locate where the slot is available
#ifdef VERBOSE
   call tell(Iam,"nfldsig =",nfldsig)
   call tell(Iam,"ntguessig_ref =",ntguessig_ref)
   call tell(Iam,"allocated(hrdifsig) =",allocated(hrdifsig))
   call tell(Iam,"allocated(hrdifsig_all) =",allocated(hrdifsig_all))
   call tell(Iam,"size(hrdifsig,1) =",size(hrdifsig,1))
   call tell(Iam,"size(hrdifsig_all,1) =",size(hrdifsig_all,1))
   call tell(Iam,"hrdifsig_all(ntguessig_ref) =",hrdifsig_all(ntguessig_ref))
#endif
   endif

   deallocate(ALARM, stat=STATUS); VERIFY_(STATUS)

   call ESMF_ClockSet(clock, direction=ESMF_DIRECTION_REVERSE, rc=status); VERIFY_(STATUS)
   call ESMF_ClockAdvance(clock, rc=STATUS)                         ; VERIFY_(STATUS)
   call ESMF_ClockSet(clock, direction=ESMF_DIRECTION_FORWARD, rc=status); VERIFY_(STATUS) 

   if(IamRoot) call myclock_show_(clock,IAm,'reversed-clock')


   ! restore the state of the alarms
   DO I = 1, nalarms
      call ESMF_AlarmSet(alarmList(I), ringTime=alarmRingTime(I), ringing=ringingState(I), rc=status)
      VERIFY_(STATUS)
   END DO
          
   deallocate(alarmList, alarmRingTime, ringingState)


   RETURN_(ESMF_SUCCESS)

   end subroutine GSI_GridCompSetAlarms

! Utilities: 

!-------------------------------------------------------------------------
  subroutine GSI_GridCompSP2NP_(rbufr,sfcin)
!-------------------------------------------------------------------------
! NCEP native Gaussian gridded fields have latitudes reversed and the poles
! excluded. This routine reverses latitudes and adds poles.
  implicit none
  real(r_single),dimension(:,:),intent(inout) :: rbufr
  real(r_single),dimension(:,:),intent(in ) :: sfcin
  integer :: im,jm,j
  character(len=*), parameter :: IAm='GSI_GridCompSP2NP_'

  im=size(rbufr,1)
  jm=size(rbufr,2)
  if(jm>size(sfcin,2))then
     rbufr(:, 1)=sum(sfcin(:,jm-2))/im	! add South Pole points
  else
     rbufr(:, 1)=sfcin(:,jm)		! add South Pole points
  endif
  do j=2,jm-1
				! for j :=    2,   3, ..., jm-2,jm-1
				!  jm-j == jm-2,jm-3, ...,    2,   1
    rbufr(:,j)=sfcin(:,jm-j) 
  end do
  if(jm>size(sfcin,2))then
     rbufr(:,jm)=sum(sfcin(:,   1))/im	! add North Pole points
  else
     rbufr(:,jm)=sfcin(:,   1)		! add North Pole points
  endif
  end subroutine GSI_GridCompSP2NP_

!-------------------------------------------------------------------------
   subroutine SwapVr8_(fld)
!-------------------------------------------------------------------------
   implicit none
   real(r_double),intent(inout) ::  fld(:,:,:)
   real(r_double),allocatable   :: work(:,:,:)
   character(len=*), parameter :: IAm='GSI_GridComp.SwapVr8_'
   integer im, jm, km
   if(.not.doVflip)return
   im   = size(fld,1)
   jm   = size(fld,2)
   km   = size(fld,3)
   allocate (work(im,jm,km))
   work = fld
   fld(:,:,km:1:-1) = work(:,:,1:km:+1)
   deallocate (work)
   end subroutine SwapVr8_
   subroutine SwapVr4_(fld)
!-------------------------------------------------------------------------
   implicit none
   real(r_single),intent(inout) ::  fld(:,:,:)
   real(r_single),allocatable   :: work(:,:,:)
   character(len=*), parameter :: IAm='GSI_GridComp.SwapVr4_'
   integer im, jm, km
   if(.not.doVflip)return
   im   = size(fld,1)
   jm   = size(fld,2)
   km   = size(fld,3)
   allocate (work(im,jm,km))
   work = fld
   fld(:,:,km:1:-1) = work(:,:,1:km:+1)
   deallocate (work)
   end subroutine SwapVr4_

!-------------------------------------------------------------------------
   subroutine SwapIJ3r8_(aij,aji)
!-------------------------------------------------------------------------
! transpose IJK-ordered array to JIK-ordered array
   implicit none
   real(r_single),dimension(:,:,:),intent(in) :: aij
   real(r_kind),  dimension(:,:,:),intent(inout) :: aji
   character(len=*), parameter :: IAm='GSI_GridComp.SwapIJ3r8_'
   integer :: i,j,k,isz,jsz,ksz
!
   isz=size(aij,1)
   jsz=size(aij,2)
   ksz=size(aij,3)
   do k=1,ksz
      do i=1,isz
         aji(1:jsz,i,k)=aij(i,1:jsz,k)
      end do
   end do
   call Halo3d_Undef_OutJIr_(aji(:,:,1:ksz))
   call GSI_GridCompSwapV_(aji(:,:,1:ksz))
   end subroutine SwapIJ3r8_

!-------------------------------------------------------------------------
   subroutine SwapJI3r8_(aij,aji)
!-------------------------------------------------------------------------
! transpose JI-ordered array to IJ-ordered array
   implicit none
   real(r_single),dimension(:,:,:),intent(inout) :: aij
   real(r_kind),dimension(:,:,:),intent(in) :: aji
   character(len=*), parameter :: IAm='GSI_GridComp.SwapJI3r8_'
   integer :: i,j,k,isz,jsz,ksz
!
   isz=size(aji,2)
   jsz=size(aji,1)
   ksz=size(aji,3)
   do k=1,ksz
      do i=1,isz
         aij(i,1:jsz,k)=aji(1:jsz,i,k)
      end do
   end do
   call GSI_GridCompSwapV_(aij(:,:,1:ksz))

   end subroutine SwapJI3r8_

!-------------------------------------------------------------------------
   subroutine SwapIJ2r8_(aij,aji)
!-------------------------------------------------------------------------
! transpose IJ-ordered array to JI-ordered array
   implicit none
   real(r_single),dimension(:,:),intent(in) :: aij
   real(r_kind),dimension(:,:),intent(inout) :: aji
   character(len=*), parameter :: IAm='GSI_GridComp.SwapIJ2r8_'
!
   integer :: i,j,isz,jsz
!
   isz=size(aij,1)
   jsz=size(aij,2)
   do i=1,isz
      aji(1:jsz,i)=aij(i,1:jsz)
   end do
   call Halo2d_Undef_OutJIr_(aji)

   end subroutine SwapIJ2r8_

!-------------------------------------------------------------------------
   subroutine SwapIJ2i_(aij,aji)
!-------------------------------------------------------------------------
! transpose IJ-ordered int-array to JI-ordered int-array
   implicit none
   integer, dimension(:,:),intent(in   ) :: aij
   integer, dimension(:,:),intent(inout) :: aji
!
   character(len=*), parameter :: IAm='GSI_GridComp.SwapIJ2i_'
   integer :: i,j,isz,jsz
!
     isz=size(aij,1)
     jsz=size(aij,2)
     do i=1,isz
       aji(1:jsz,i)=aij(i,1:jsz)
     end do

   end subroutine SwapIJ2i_

!-------------------------------------------------------------------------
   subroutine SwapJI2r_(aij,aji)
!-------------------------------------------------------------------------
! transpose IJ-ordered array to JI-ordered array
   implicit none
   real(r_single),dimension(:,:),intent(inout) :: aij
   real(r_kind), dimension(:,:),intent(in) :: aji
!
   character(len=*), parameter :: IAm='GSI_GridComp.SwapJI2r_'
   integer :: i,j,isz,jsz
!
   jsz=size(aji,1)
   isz=size(aji,2)
   do i=1,isz
     aij(i,1:jsz)=aji(1:jsz,i)
   end do

   end subroutine SwapJI2r_

!-------------------------------------------------------------------------
   subroutine SwapJI2i_(aij,aji)
!-------------------------------------------------------------------------
! transpose IJ-ordered array to JI-ordered array
   implicit none
   integer, dimension(:,:),intent(inout) :: aij
   integer, dimension(:,:),intent(in)    :: aji
!
   character(len=*), parameter :: IAm='GSI_GridComp.SwapJI2i_'
   integer :: i,j,isz,jsz
!
   jsz=size(aji,1)
   isz=size(aji,2)
   do i=1,isz
     aij(i,1:jsz)=aji(1:jsz,i)
   end do

   end subroutine SwapJI2i_

!-------------------------------------------------------------------------
   subroutine Halo2d_Undef_OutJIr_(aji)
!-------------------------------------------------------------------------
!  !REVISION HISTORY:
!    21Aug2007 Todling Implemented to remove MAPL_UNDEF from halo
!                      This is a temporary fix that needs better handling
!                      to differentiate between scalar and vector fields
!  THIS ASSUMES horizontal HW=1
!-------------------------------------------------------------------------
   implicit none
   real(r_kind),  dimension(:,:),intent(inout) :: aji
   character(len=*), parameter :: IAm='GSI_GridComp.Halo2d_Undef_OutJIr_'
   integer :: j,i,k,isz,jsz
   jsz=size(aji,1)
   isz=size(aji,2)
   do i=1,isz
      if(abs(aji(1  ,i)-MAPL_UNDEF).lt.1.e-2) aji(1  ,i) = aji(2    ,i)
      if(abs(aji(jsz,i)-MAPL_UNDEF).lt.1.e-2) aji(jsz,i) = aji(jsz-1,i)
   end do
   end subroutine Halo2d_Undef_OutJIr_

!-------------------------------------------------------------------------
   subroutine Halo3d_Undef_OutJIr_(aji)
!-------------------------------------------------------------------------
!  !REVISION HISTORY:
!    21Aug2007 Todling Implemented to remove MAPL_UNDEF from halo
!                      This is a temporary fix that needs better handling
!                      to differentiate between scalar and vector fields
!  THIS ASSUMES horizontal HW=1
!-------------------------------------------------------------------------
   implicit none
   real(r_kind),  dimension(:,:,:),intent(inout) :: aji
   character(len=*), parameter :: IAm='GSI_GridComp.Halo3d_Undef_OutJIr_'
   integer :: j,i,k,isz,jsz,ksz
   jsz=size(aji,1)
   isz=size(aji,2)
   ksz=size(aji,3)
   do k = 1,ksz
      do i=1,isz
         if(abs(aji(1  ,i,k)-MAPL_UNDEF).lt.1.e-2) aji(1  ,i,k) = aji(2    ,i,k)
         if(abs(aji(jsz,i,k)-MAPL_UNDEF).lt.1.e-2) aji(jsz,i,k) = aji(jsz-1,i,k) 
      end do
   end do
   end subroutine Halo3d_Undef_OutJIr_
                                                                                                                      
!-------------------------------------------------------------------------
   subroutine GSI_GridCompFlipLons2_(q)
!-------------------------------------------------------------------------
     implicit none
     real,dimension(:,:),intent(inout) :: q
     character(len=*), parameter :: IAm='GSI_GridCompFlipLons2_'
     integer :: j
     do j=1,size(q,2)
        call GSI_GridCompFlipLonsi_(q(:,j))
     end do
   end subroutine GSI_GridCompFlipLons2_

!-------------------------------------------------------------------------
   subroutine GSI_GridCompFlipLonsi_(q)
!-------------------------------------------------------------------------
     implicit none
     real,dimension(:),intent(inout) :: q
     character(len=*), parameter :: IAm='GSI_GridCompFlipLonsi_'
     integer :: im
     real,dimension(size(q,1)/2) :: d
     im=size(q,1)
     d(     1:im/2) = q(     1:im/2)
     q(     1:im/2) = q(im/2+1:im  )
     q(im/2+1:im  ) = d(     1:im/2)
   end subroutine GSI_GridCompFlipLonsi_

! The following routines exist somewhere else , but are not readily
! available

!-------------------------------------------------------------------
   subroutine get_DateStamp (clock, DateStamp, expid, offset, rc)
!-------------------------------------------------------------------

! This can be found in GEOShistory.

    type (ESMF_Clock)                 :: clock
    character(len=ESMF_MAXSTR)        :: DateStamp
    character(len=ESMF_MAXSTR)        :: expid
    type(ESMF_TimeInterval), optional :: offset
    integer, optional                 :: rc

    character(len=*), parameter :: IAm='GSI_GridComp.get_DateStamp'

    type(ESMF_Time)                   :: currentTime
    type(ESMF_TimeInterval)           :: TimeStep
    character(len=ESMF_MAXSTR)        :: TimeString
    integer                           :: secs
    character(len=ESMF_MAXSTR)        :: TimeStamp
    character                         :: String(ESMF_MAXSTR)

    character*4 year
    character*2 month
    character*2 day
    character*2 hour
    character*2 minute
    character*2 second

    equivalence ( string(01),TimeString )
    equivalence ( string(01),year       )
    equivalence ( string(06),month      )
    equivalence ( string(09),day        )
    equivalence ( string(12),hour       )
    equivalence ( string(15),minute     )
    equivalence ( string(18),second     )

    call ESMF_ClockGet (clock, currTime=currentTime, rc=rc)
    if (present(offset)) then
       currentTime = currentTime - offset
    end if
    call ESMF_TimeGet  (currentTime, timeString=TimeString, rc=rc)
    call ESMF_ClockGet (clock, timeStep=TimeStep, rc=rc)
    call ESMF_TimeIntervalGet (TimeStep, S=secs, rc=rc)

    DateStamp = year // month // day // '_' // hour // minute // second // 'z'
    TimeStamp = '   Date: ' // year // '/' // month // '/' // day
    TimeStamp = trim(TimeStamp) // '  Time: ' // timestring(12:19)

    if (.not. present(OFFSET)) then
       call WRITE_PARALLEL ( 'Expid: ' // trim(expid) // trim(TimeStamp) )
    endif

   end subroutine get_DateStamp

!-------------------------------------------------------------------------
!   NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, DAS  !
!-------------------------------------------------------------------------
!BOPI
!
! !ROUTINE: ncep_rwsurf_ - read and/or create SSI surface background file
!
! !INTERFACE:

   subroutine ncep_rwsurf_ ( verbose_, infile, glat2, glon, rc, & 
                             wout, jrec, fld )

! !USES:

   implicit none

! INPUT PARAMETERS:

   logical :: verbose_
   integer :: glat2, glon
   character(len=*), intent(in)           :: infile
   logical,          intent(in), optional :: wout
   integer,                      optional :: jrec     ! set to record 
                                                      ! of fld to be returned
! OUTPUT PARAMETERS:

   real(r_single),intent(inout),optional :: fld(glon,glat2) ! requested output field from file
   integer, intent(out) :: rc                        ! return error code

! !DESCRIPTION: This module handles surface fields form the NCEP surface files.
!
!EOPI
!-------------------------------------------------------------------------

   character(len=*), parameter :: IAm='GSI_GridComp.get_ncep_rwsurf_'

   real(r_single) :: sdum4
   real(r_single), allocatable :: fdum4(:,:)
   real(r_single), allocatable :: fdum42(:,:,:)
   real(r_single), allocatable :: albedo4(:,:,:)
   real(r_single), allocatable :: zenith2(:,:,:)

   logical :: writeout
   integer :: icount, irec, status
   character(8),dimension(4):: labfix
   real(r_single) yhour
   integer latd,lonl,version
   integer,dimension(4):: igdate,iadate
   integer,allocatable,dimension(:):: lonsperlat

   rc       = 0
   irec     = 0
   writeout = .false.
   if (present(jrec)) then
       if(.not.present(fld))then
          print *, trim(Iam), ': Error, missing fld array'
          rc = 1
          return
       endif
       irec = jrec
   endif

   ! unit 34 is a temporary file that contains surface fields
   ! not available from GMAO and so the corresponding NCEP fields are used

   open(34,file=infile,form='unformatted')
   if(present(wout)) then
      if(wout)then
         writeout = .true.
      endif
   endif
   allocate ( fdum4  (glon,glat2),   &
              fdum42 (glon,glat2,2), & 
              albedo4(glon,glat2,4), &
              zenith2(glon,glat2,2), &
              stat=STATUS)
   VERIFY_(STATUS)
   read(34) 
   read(34) yhour,IgDATE,LONL,LATD,version
   allocate ( lonsperlat(latd/2), stat=STATUS)
   VERIFY_(STATUS)

   rewind(34)
   read(34) labfix
   read(34) yhour,IgDATE,LONL,LATD,version,lonsperlat
   icount = 0
   read(34) fdum4                      ! record 1  tsf
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum42                     ! record 2  soilm
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum42
          fld = fdum42(:,:,1)
       end if
   read(34) fdum4                      ! record 3  snow
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum42                     ! record 4  soilt
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum42
          fld = fdum42(:,:,1)
       end if
   read(34) fdum4                      ! record 5  tg3
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 6  zor
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 7  cv
       icount = icount + 1
       if(verbose_) print*,Iam,': ',icount, ' ', maxval(fdum4)
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 8  cvb
       icount = icount + 1
       if(verbose_) print*,Iam,': ',icount,' ', maxval(fdum4)
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 9  cvt
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) albedo4                    ! record 10 albedo
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) albedo4
          fld = albedo4(:,:,1)
       end if
   read(34) fdum4                      ! record 11 slimsk
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 12 vegetation cover
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 13 canopy water
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 14 f10m
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 15 vegetation type
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 16 soil type
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) zenith2                    ! record 17 zenith
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) zenith2
          fld = zenith2(:,:,1)
       end if
   read(34) fdum4                      ! record 18 ustar
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 19 ffmm
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if
   read(34) fdum4                      ! record 20 ffhh
       icount = icount + 1
       if(irec==icount) then
          if(verbose_) print *,Iam,': reading record ',irec
          if(writeout) write(35) fdum4
          fld = fdum4
       end if

   close(34)

   deallocate ( lonsperlat, stat=STATUS)
   VERIFY_(STATUS) 
   deallocate ( fdum4,fdum42,albedo4,zenith2 , stat=STATUS)
   VERIFY_(STATUS)

   if(verbose_) print * , trim(Iam), ': done'
   end subroutine ncep_rwsurf_

! undef_2ssi exists in transf?
! Unfortunately present code needs it r4; transf has it r8.

!-------------------------------------------------------------------------
   subroutine undef_2ssi (var, undef_in, cwhere, reltol, verb, vname )
!-------------------------------------------------------------------------

     implicit none

     real,dimension(:,:),intent(inout) ::     var
     real             ,intent(in) :: undef_in
     character(len=*) ,intent(in) :: cwhere
     real,optional    ,intent(in) :: reltol
     logical,optional ,intent(in) :: verb
     character(len=*),optional,intent(in) :: vname

     character(len=*), parameter :: IAm='GSI_GridComp.undef_2ssi_'

     real :: reltol_

     reltol_=.01
     if(present(reltol)) reltol_=reltol

     if(verb==.true.) then
        call lookfor_(fundef_ssi() ,reltol_,var,6,cwhere,vname)	! -.999e33
        call lookfor_(undef_in     ,reltol_,var,6,cwhere,vname)	! given
        call lookfor_(1.e12        ,reltol_,var,6,cwhere,vname)	! 1.e12
     endif
     
     call udf2udf2d_(var, undef_in,fundef_ssi(),cwhere,	&
          reltol=reltol, verb=verb, vname=vname )

   contains

     function fundef_ssi()
       implicit none
       real :: fundef_ssi
       fundef_ssi=UNDEF_SSI_
     end function fundef_ssi

     subroutine lookfor_(spval,reltol,v,lu,cwhere,vname)
       implicit none
       real,intent(in) :: spval
       real,intent(in) :: reltol
       real,dimension(:,:),intent(in) :: v
       integer,intent(in) :: lu
       character(len=*),intent(in) :: cwhere,vname
       character(len=*), parameter :: IAm='GSI_GridComp.lookfor_'

       integer :: n,m
       n=count(abs(v-spval)<=abs(reltol*spval))
       m=size(v)
!       if(n/=0) write(lu,'(2a,e12.6,3a,i6,a,i6)') trim(where),	&
!            ': find ',spval,' in variable "',trim(vname),'", ',	&
!            n,' out of ',m
     end subroutine lookfor_

   end subroutine undef_2ssi


!-------------------------------------------------------------------------
   subroutine udf2udf2d_ (var,undef_in,undef_out,cwhere,	&
      	reltol,verb,vname)
!-------------------------------------------------------------------------
     implicit none

     real,intent(in) :: undef_in
     real,intent(in) :: undef_out
     character(len=*),intent(in) :: cwhere

     real,dimension(:,:),intent(inout) ::     var

     real,optional,intent(in) :: reltol
     logical,optional,intent(in) :: verb
     character(len=*),optional,intent(in) :: vname

     character(len=*), parameter :: IAm='GSI_GridComp.udf2udf2d_'

     integer  i,j
     integer  nonmiss
     real     tol
     real     val_max_in, val_max_out
     logical :: verbose_

     verbose_=.false.
     if(present(verb)) verbose_=verb
     
      tol   = 0.01
      if(present(reltol)) tol=reltol
      tol   = tol*undef_in

!     Make sure there is no precision problem with undef's
!     ----------------------------------------------------
      nonmiss = 0
      val_max_in  = maxval(var)
      do j = 1, size(var,2)
         do i = 1, size(var,1)
            if ( abs(var(i,j)-undef_in) .lt. tol ) then
               nonmiss = nonmiss + 1
               var(i,j) = undef_out
            end if
         end do
      end do
      val_max_out = maxval(abs(var))
      
      if (nonmiss.ne.0) then
         if(verbose_) then
	   if(present(vname)) then
              if(MAPL_AM_I_ROOT().and.verbose) write(6,*) &
                   'No. of UNDEF_in fixed to UNDEF_out for "'//	&
                   trim(vname)//'"',nonmiss
	   else
              if(MAPL_AM_I_ROOT().and.verbose) write(6,*) &
                   'No. of UNDEF_in fixed to UNDEF_out',nonmiss
	   endif
        endif
        
        if ( val_max_out .gt. abs(undef_out) ) then
           if (MAPL_AM_I_ROOT()) then
               write(6,*) &
	       	     ' Largest  abs(value) on  input: ',  val_max_in
               write(6,*) &
	   	     ' Largest  abs(value) on output: ',  val_max_out
               write(6,*) &
		     ' Undef_in     value spec.: ',  undef_in
               write(6,*) &
		     ' Undef_out    value spec.: ',  undef_out
           endif
        end if
     end if

     return
   end subroutine udf2udf2d_

subroutine myClock_show_(aClock,who,what,rc)
  use ESMF, only: ESMF_SUCCESS
  use ESMF, only: ESMF_Clock
  use ESMF, only: ESMF_ClockGet
  use ESMF, only: ESMF_Time
  use ESMF, only: ESMF_TimeGet
  use ESMF, only: ESMF_TimeInterval
  use ESMF, only: ESMF_TimeIntervalGet
  use ESMF, only: ESMF_Direction_Flag
  use ESMF, only: ESMF_DIRECTION_FORWARD
  use ESMF, only: ESMF_DIRECTION_REVERSE
  use mpeu_util, only: tell,warn,perr,die
  implicit none
  type(ESMF_Clock),intent(in):: aClock
  character(len=*),intent(in):: who
  character(len=*),intent(in):: what
  integer(i_kind),optional,intent(out):: rc

  character(len=*),parameter:: myname_=myname//"::myClock_show_"
  character(len=20):: timestr
  type(ESMF_Time):: bTime,cTime,eTime
  type(ESMF_TimeInterval):: dTime
  integer(i_kind):: d,h,m,s,ier
! type(ESMF_Direction):: direction
#define MYNAME_ myname_//'/'//trim(who)//'/'//trim(what)
_ENTRY_(MYNAME_)

  if(present(rc)) rc=ESMF_SUCCESS
	return

  call ESMF_ClockGet(aClock, startTime=btime,stopTime=etime, &
  			      currTime=ctime,timeStep=dTime, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMF_ClockGet(), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(MYNAME_)
	  return
	endif

  call ESMF_TimeGet(btime,timeStringISOFrac=timestr)
  call tell(who,'startTime('//trim(what)//') = '//trim(timestr))
  call ESMF_TimeGet(ctime,timeStringISOFrac=timestr)
  call tell(who,' currTime('//trim(what)//') = '//trim(timestr))
  call ESMF_TimeGet(etime,timeStringISOFrac=timestr)
  call tell(who,' stopTime('//trim(what)//') = '//trim(timestr))

  call ESMF_TimeIntervalGet(dtime,d=d,h=h,m=m,s=s)
  call tell(who,' timeStep('//trim(what)//') =',(/((d*100+h)*100+m)*100+s/),format='(i10.8)')
_EXIT_(MYNAME_)
end subroutine myClock_show_

subroutine echo_time_ ( CLOCK, msg )
    use m_StrTemplate
    implicit none
    type(ESMF_Clock) :: Clock
    character(len=*) :: msg

    type(ESMF_Time)  :: Time
    integer :: YY, MM, DD, HH, MN, SC, nymd, nhms
    integer :: rc,status

    call ESMF_ClockGet(Clock, CurrTime=Time)
    call ESMF_TimeGet(Time, yy=YY, mm=MM, dd=DD, h=HH, m=MN, s=SC, rc=rc)
    if(MAPL_AM_I_ROOT()) write(6,'(a,1x,i4.4,4(a,i2.2))') trim(msg)//" The current TIME is ", &
                                        YY, "/", MM, "/", DD, " ", HH, ":", MN
end subroutine echo_time_

   end module GSI_GridCompMod
