!#define DEBUG_TRACE
!#define DEBUG_METACOMP
!#define DBLE_STATE
module stub_GEOSagcmPert_GridCompMod
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module stub_GEOSagcmPert_GridCompMod
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2011-02-11
!
! abstract: this is a stub of GEOSagcmPert_GridCompMod
!
! program history log:
!   2011-02-11  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! module interface:

  use kinds, only: i_kind
#ifdef DBLE_STATE
  use ESMF, only: r_ESkind => ESMF_KIND_R8
# define _PRECISION_ PRECISION=ESMF_KIND_R8,
#else
  use ESMF, only: r_ESkind => ESMF_KIND_R4
# define _PRECISION_
#endif
  use mpeu_util, only: tell, warn, perr, die
  implicit none
  private	! except
  public :: setServices		! the minimum interface

  public :: TLMphase		! flag of operation is in TL phase (TLMphase_addImport_getExport)
  public :: ADMphase		! flag of operation is in AD phase (ADMphase_addImport_getExport)

  public :: phase		! currPhase = phase(adjoint=.false.,addIMPORT=.true.,getEXPORT=.true.)
  public :: phase_opAdjoint	! opAdjoint = phase_opAdjoint(phase)
  public :: phase_addImport	! addImport = phase_addimport(phase)
  public :: phase_getExport	! getExport = phase_getexport(phase)
  
  public :: gcdiag

  interface phase
  	module procedure phase_  ; end interface
  interface phase_opAdjoint
  	module procedure adjoint_; end interface
  interface phase_addImport
  	module procedure import_ ; end interface
  interface phase_getExport
  	module procedure export_ ; end interface

  interface gcdiag
  	module procedure show_   ; end interface

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='stub_GEOSagcmPert_GridCompMod'

  integer, parameter :: TLMPhase_addIMPORT_getEXPORT = 1	! 000 +1
  integer, parameter :: ADMPhase_addIMPORT_getEXPORT = 2	! 001 +1
  integer, parameter :: TLMPhase_notIMPORT_getEXPORT = 3	! 010 +1
  integer, parameter :: ADMPhase_notIMPORT_getEXPORT = 4	! 011 +1
  integer, parameter :: TLMPhase_addIMPORT_notEXPORT = 5	! 100 +1
  integer, parameter :: ADMPhase_addIMPORT_notEXPORT = 6	! 101 +1
  integer, parameter :: TLMPhase_notIMPORT_notEXPORT = 7	! 110 +1
  integer, parameter :: ADMPhase_notIMPORT_notEXPORT = 8	! 111 +1

  integer(i_kind),parameter:: TLMphase = TLMphase_addIMPORT_getEXPORT
  integer(i_kind),parameter:: ADMphase = ADMphase_addIMPORT_getEXPORT

!! all inState_ components
  real(r_ESkind),allocatable,dimension(:,:,:),save::  u_
  real(r_ESkind),allocatable,dimension(:,:,:),save::  v_
  real(r_ESkind),allocatable,dimension(:,:,:),save:: tv_
  real(r_ESkind),allocatable,dimension(:,:,:),save:: dp_
  real(r_ESkind),allocatable,dimension(:,:,:),save:: qv_
  real(r_ESkind),allocatable,dimension(:,:,:),save:: ql_
  real(r_ESkind),allocatable,dimension(:,:,:),save:: qi_
  real(r_ESkind),allocatable,dimension(:,:,:),save:: o3_
  integer(i_kind),save:: last_phase_   = -1
  integer(i_kind),save:: call_counter_ = -1

  logical,save:: inState_created_=.false.

  	! These are expected from a PGCM state
  character(len=*),parameter::    u_GCMpert='U'		! = u
  character(len=*),parameter::    v_GCMpert='V'		! = v
  character(len=*),parameter::   tv_GCMpert='TV'	! = tv
  character(len=*),parameter::   dp_GCMpert='DP'	! = d(p3d,ps)/dk
  character(len=*),parameter::   qv_GCMpert='QV'	! = q
  character(len=*),parameter::   ql_GCMpert='QL'	! = cw - qi
  character(len=*),parameter::   qi_GCMpert='QI'	! = qi (=0)
  character(len=*),parameter::   o3_GCMpert='O3'	! = oz

  integer(i_kind),parameter:: GENERIC_ERROR=-9999

#include "MAPL_ErrLog.h"
#include "mytrace.H"
contains
subroutine setServices(gc,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine setServices
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2011-02-11
!
! abstract: the minimum interface entry
!
! program history log:
!   2011-02-11  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_GridComp
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: MAPL_GenericSetServices
  use MAPL_Mod, only: MAPL_GridCompSetEntryPoint
  use MAPL_Mod, only: MAPL_AM_I_ROOT
#ifdef DEBUG_METACOMP
  use MAPL_Mod, only: MAPL_MetaCompPrint
#endif

  use ESMF, only: ESMF_METHOD_INITIALIZE
  use ESMF, only: ESMF_METHOD_RUN
  use ESMF, only: ESMF_METHOD_FINALIZE

  use mpeu_util, only: perr,die
  implicit none
  type(ESMF_GridComp), intent(inout):: gc
  integer            , intent(  out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer:: phase0, phase1
  character(len=*),parameter :: myname_=myname//'::setServices'
_ENTRY_(myname_)

  rc=ESMF_SUCCESS

  call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE , Initialize, rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_INITIALIZE), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  	! register for run(), phase 1
  call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN  , Run       , rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_RUN,1), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  	! register for run(), phase 2
  call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN  , Run       , rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_RUN,2), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  phase0=phase(import=.true. )
  phase1=phase(import=.false.)

  if(phase0/=phase1) then
  	! register for run(), phase 3
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN  , Run       , rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_RUN,3), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  	! register for run(), phase 4
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN  , Run       , rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_RUN,4), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

    phase0=phase(export=.true. )
    phase1=phase(export=.false.)

    if(phase0/=phase1) then

  	! register for run(), phase 5
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN  , Run       , rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_RUN,5), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  	! register for run(), phase 6
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN  , Run       , rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_RUN,6), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  	! register for run(), phase 7
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN  , Run       , rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_RUN,7), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  	! register for run(), phase 8
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN  , Run       , rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_RUN,8), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

      if(MAPL_AM_I_ROOT()) write(*,'(1x,2a)') trim(myname_),' -- 8 RUN() phases set'
    else
      if(MAPL_AM_I_ROOT()) write(*,'(1x,2a)') trim(myname_),' -- 4 RUN() phases set'
    endif
  else
      if(MAPL_AM_I_ROOT()) write(*,'(1x,2a)') trim(myname_),' -- 2 RUN() phases set'
  endif

  call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE, Finalize  , rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GridCompSetEntryPoint(ESMF_METHOD_FINALIZE), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  call addSpecs_( GC, rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'addSpecs_(), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  call MAPL_GenericSetServices ( GC, RC=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GenericSetServices(), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

#ifdef DEBUG_METACOMP
call MAPL_MetaCompPrint(GC,header='before exiting '//myname_)
#endif
_EXIT_(myname_)
end subroutine setServices
subroutine addSpecs_(gc, rc)
  use ESMF, only: ESMF_SUCCESS
  use ESMF, only: ESMF_GridComp
  use ESMF, only: ESMF_KIND_R4
  use ESMF, only: ESMF_KIND_R8
  use MAPL_Mod, only: MAPL_AddImportSpec
  use MAPL_Mod, only: MAPL_AddExportSpec
  use MAPL_Mod, only: MAPL_DimsHorzVert
  use MAPL_Mod, only: MAPL_VLocationCenter
  
  use mpeu_util, only: perr,die
  implicit none

  type(ESMF_GridComp),intent(inout):: gc
  integer(i_kind),intent(out):: rc

  character(len=*),parameter:: myname_=myname//'::addSpecs_'

_ENTRY_(myname_)
    rc=ESMF_SUCCESS

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME =  u_GCMpert,                                  &
         LONG_NAME  = 'eastward_wind_analysis_increment',          &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddImportSpec("'// u_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME =  v_GCMpert,                                  &
         LONG_NAME  = 'northward_wind_analysis_increment',         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddImportSpec("'// v_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = tv_GCMpert,                                  &
         LONG_NAME  = 'virtual_temperature_analysis_increment',    &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddImportSpec("'//tv_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = dp_GCMpert,                                  &
         LONG_NAME  = 'pressure_thickness_analysis_increment',     &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddImportSpec("'//dp_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = qv_GCMpert,                                  &
         LONG_NAME  = 'specific_humidity_analysis_increment',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddImportSpec("'//qv_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = ql_GCMpert,                                  &
         LONG_NAME  = 'cloud_liquid_water_analysis_increment',     &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddImportSpec("'//ql_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = qi_GCMpert,                                  &
         LONG_NAME  = 'cloud_ice_analysis_increment',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddImportSpec("'//qi_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = o3_GCMpert,                                  &
         LONG_NAME  = 'ozone_analysis_increment',                  &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddImportSpec("'//o3_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

! !EXPORT STATE:

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME =  u_GCMpert,                                  &
         LONG_NAME  = 'perturbation_eastward_wind_on_A-grid',      &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddExportSpec("'// u_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME =  v_GCMpert,                                  &
         LONG_NAME  = 'perturbation_northward_wind_on_A-grid',     &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddExportSpec("'// v_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = tv_GCMpert,                                  &
         LONG_NAME  = 'perturbation_air_virtual_temperature',      &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddExportSpec("'//tv_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = dp_GCMpert,                                  &
         LONG_NAME  = 'perturbation_air_pressure_thickness',       &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddExportSpec("'//dp_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = qv_GCMpert,                                  &
         LONG_NAME  = 'perturbation_specific_humidity',            &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddExportSpec("'//qv_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = ql_GCMpert,                                  &
         LONG_NAME  = 'perturbation_cloud_liquid_water_mixing_ration',&
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddExportSpec("'//ql_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = qi_GCMpert,                                  &
         LONG_NAME  = 'perturbation_cloud_ice_mixing_ration',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddExportSpec("'//qi_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = o3_GCMpert,                                  &
         LONG_NAME  = 'perturbation_ozone_mole_mixing_ration',     &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert, _PRECISION_               &
         VLOCATION  = MAPL_VLocationCenter,             RC=rc      )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_AddExportSpec("'//o3_GCMpert//'"), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif
_EXIT_(myname_)
end subroutine addSpecs_

subroutine show_(gc,rc,where)
    use mpeu_util, only: tell,perr
    use ESMF , only: ESMF_GridComp
    use ESMF , only: ESMF_GridCompPrint
    use ESMF , only: ESMF_GridCompGet
    use ESMF , only: ESMF_State
    use ESMF , only: ESMF_MAXSTR
    use ESMF , only: ESMF_StateGet
    use ESMF , only: ESMF_SUCCESS
    use MAPL_Mod , only: MAPL_MetaComp
    use MAPL_Mod , only: MAPL_Get
    use MAPL_Mod , only: MAPL_GetObjectFromGC

    implicit none
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    integer, optional,   intent(  out) :: RC     ! Error code
    character(len=*),optional,intent(in):: where

    character(len=ESMF_MAXSTR)          :: IAm 
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME,GCname,IMname

    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_GridComp),      pointer  :: GCS(:)
    character(len=ESMF_MAXSTR),pointer  :: GCNames(:)

    character(len=*),parameter:: myname_=myname//"::show_"
    integer:: i

_ENTRY_(myname_)
    if(present(where)) call tell(myname_,where)

    call ESMF_GridCompPrint(GC, rc=rc)
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_GridCompPrint(), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=rc)
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_GridCompGet(name), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    Iam = trim(COMP_NAME) // '-' //myname_
    call tell(myname_,'%comp_name =',trim(comp_name))

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=rc)
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_GetObjectFromGC(mapl), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GCNames=GCNames, RC=rc )
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'MAPL_Get(GCS, GIM, ...), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif
    nullify(MAPL)

    call tell(myname_,trim(comp_name)//'.MAPL.#GCNames =',size(GCNames))
    call tell(myname_,trim(comp_name)//'.MAPL.#GCS     =',size(GCS    ))
    call tell(myname_,trim(comp_name)//'.MAPL.#GIM     =',size(GIM    ))
    do i=1,size(GCNames)
	GCname=GCNames(i)
        call tell(myname_,trim(comp_name)//'.MAPL.       i =',i)
	call tell(myname_,trim(comp_name)//'.MAPL.GCname_i =',trim(GCname))
	IMname=""
        call ESMF_StateGet(GIM(i),name=IMname, RC=rc)
	 	if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_StateGet(GIM,name), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif
	call tell(myname_,trim(comp_name)//'.MAPL.IMname_i = ',trim(IMname))
    enddo
    call tell(myname_,trim(comp_name))

    nullify(GCS)
    nullify(GIM)
    nullify(GCNames)

_EXIT_(myname_)
end subroutine show_

subroutine Initialize(gc,import,export,clock,rc)
  use ESMF, only: ESMF_GridComp
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_Clock
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: MAPL_GenericInitialize
  use kinds   , only: i_kind
  use mpeu_util, only: perr,die

  implicit none
  type(ESMF_GridComp ),intent(inout):: gc
  type(ESMF_State    ),intent(inout):: import,export
  type(ESMF_Clock    ),intent(inout):: clock
  integer(kind=i_kind),intent(  out):: rc

  character(len=*),parameter:: myname_=myname//"::Initialize"

_ENTRY_(myname_)
  rc=ESMF_SUCCESS

  call MAPL_GenericInitialize(gc, import, export, clock, rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GenericInitialize(), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  call inState_create_(import)		! make inState_
  call inState_init_  (0._r_ESkind)	! inState_ = 0.

  last_phase_ = 0
_EXIT_(myname_)
end subroutine Initialize

subroutine Run(gc,import,export,clock,rc)
  use ESMF, only: ESMF_GridComp
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_Clock
  use ESMF, only: ESMF_GridCompGet
  use ESMF, only: ESMF_SUCCESS
  use kinds   , only: i_kind

  implicit none
  type(ESMF_GridComp ),intent(inout):: gc
  type(ESMF_State    ),intent(inout):: import,export
  type(ESMF_Clock    ),intent(inout):: clock
  integer(kind=i_kind),intent(  out):: rc

  character(len=*),parameter:: myname_=myname//"::Run"

  integer(i_kind):: phase, phase_
  logical:: addincs, makeexp

_ENTRY_(myname_)
  rc=ESMF_SUCCESS

  call ESMF_GridCompGet( GC, currentPhase=phase, RC=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMF_GridCompGet(currentPhase), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  phase_=TLMphase
  if(phase_opAdjoint(phase)) phase_=ADMphase
  addincs = phase_addImport(phase)
  makeexp = phase_getExport(phase)

  select case(phase_)
	! Verify the call phase
  case(TLMphase)	! valid only if last_phase_ is 0, TLMphase, or 2
    if(last_phase_/=0        .and. &
       last_phase_/=TLMphase .and. &
       last_phase_/=ADMphase) then
      call perr(myname_,'out of phase, phase =',phase )
      call perr(myname_,'    operation phase_=',phase_)
      call perr(myname_,'         last_phase_=',last_phase_)
      _EXIT_(myname_)
      return
    endif

  case(ADMphase)	! valid only if last_phase_ is 1 or 2
    if(last_phase_/=TLMphase .and. &
       last_phase_/=ADMphase) then
      call perr(myname_,'out of phase, phase =',phase )
      call perr(myname_,'    operation phase_=',phase_)
      call perr(myname_,'         last_phase_=',last_phase_)
      _EXIT_(myname_)
      return
    endif

  case default
    call perr(myname_,'unknown value, phase =',phase )
    call perr(myname_,'     operation phase_=',phase_)
    _EXIT_(myname_)
    return
  endselect

  if(phase_/= last_phase_) then
    	! reset call_counter_ and save this phase as last_phase_
    last_phase_ = phase_
    call_counter_ = 0
  endif

  if(call_counter_ == 0) call inState_init_(0._r_ESkind)	! inState_ = 0.

  call_counter_ = call_counter_ + 1

  if(addincs) call inState_addst_(import)		! inState_ = inState_+import
  if(makeexp) call inState_setst_(export,alloc=.true.)	! export   = inState_

!  call pertClock_tick(clock)

_EXIT_(myname_)
end subroutine Run

function myClock_isAtInitialTime_(clock,phase) result(yes)
  use ESMF, only: ESMF_Clock
  use ESMF, only: ESMF_Time
  use ESMF, only: ESMF_ClockGet
  use ESMF, only: ESMF_ClockPrint
  use ESMF, only: ESMF_SUCCESS
  use ESMF, only: operator(==)
  implicit none
  type(ESMF_Clock),intent(in):: clock
  integer(i_kind) ,intent(in):: phase
  logical:: yes

  character(len=*),parameter:: myname_=myname//"::myClock_isAtInitialTime_"
  type(ESMF_Time):: startTime,stoptime,currTime
  integer(i_kind):: ier
  call ESMF_ClockGet(clock, startTime=startTime, &
  			     stopTime= stopTime, &
			     currTime= currTime, &
			rc=ier)
	if(ier/=ESMF_SUCCESS) &
	  call die(myname_,'ESMF_ClockGet(), rc =',ier)

  yes=.false.
  select case(phase)
  case(TLMphase)
    yes = currTime==startTime
  case(ADMphase)
    yes = currTime== stopTime

  case default
#ifdef DEBUG_VERBOSE
    call ESMF_ClockPrint(clock,rc=ier)
#endif
    call die(myname_,'an unknown value, phase =',phase)
  endselect

end function myClock_isAtInitialTime_

subroutine copy3dvar_(var,s,r,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine get3darIFx_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-08-10
!
! abstract: Let r%v = s%v, if associated(x).
!
! program history log:
!   2010-08-10  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use kinds, only: i_kind
  use mpeu_util, only: tell,perr,die
  implicit none
  character(len=*)        ,intent(in   ):: var	! variable name in s and r
  type(ESMF_State),target ,intent(inout):: s	!
  type(ESMF_State),target ,intent(inout):: r	! r%v = s%v
  integer(i_kind),optional,intent(  out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//"::copy3dvar_"
  integer(i_kind):: ier,km
  real(r_ESkind),pointer,dimension(:,:,:):: sptr,rptr

#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'('//trim(var)//')'
_ENTRY_(_MYNAME_)

  if(present(rc)) rc=ESMF_SUCCESS

  call ESMFL_StateGetPointerToData(s,sptr,var,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'ESMFL_StateGetPointerToData(s%"'//trim(var)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(_MYNAME_)
	  return
	endif

  call ESMFL_StateGetPointerToData(r,rptr,var,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'ESMFL_StateGetPointerToData(r%"'//trim(var)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(_MYNAME_)
	  return
	endif

  rptr(:,:,:)=sptr(:,:,:)

  nullify(rptr)
  nullify(sptr)
_EXIT_(_MYNAME_)
end subroutine copy3dvar_

subroutine Finalize(gc,import,export,clock,rc)
  use ESMF, only: ESMF_GridComp
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_Clock
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: MAPL_GenericFinalize
  use kinds   , only: i_kind
  implicit none
  type(ESMF_GridComp ),intent(inout):: gc
  type(ESMF_State    ),intent(inout):: import,export
  type(ESMF_Clock    ),intent(inout):: clock
  integer(kind=i_kind),intent(  out):: rc

  character(len=*),parameter :: myname_=myname//"::Finalize"

_ENTRY_(myname_)
  rc=ESMF_SUCCESS

  call MAPL_GenericFinalize(gc,import,export,clock,rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'MAPL_GenericFinalize(), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  call inState_destroy_(rc=rc)		! destroy inState_
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'inState_destroy_(), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

_EXIT_(myname_)
end subroutine Finalize


subroutine inState_create_(template,rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use kinds, only: i_kind

  implicit none
  type(ESMF_State),intent(inout):: template
  integer(i_kind),optional,intent(out):: rc

  character(len=*),parameter:: myname_=myname//"::inState_create_"
  character(len=*),parameter:: varn="DP"
  real(r_ESkind),pointer,dimension(:,:,:):: ptr3
  integer(i_kind):: ier,l,m,n

_ENTRY_(myname_)
  if(present(rc)) rc=ESMF_SUCCESS

  if(inState_created_) then
    call perr(myname_,'inState_, already created, inState_created_ =', &
      inState_created_)
    if(.not.present(rc)) call die(myname_)
    rc=GENERIC_ERROR
    _EXIT_(myname_)
    return
  endif

  call ESMFL_StateGetPointerToData(template,ptr3,varn,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMFL_StateGetPointerToData("'//varn//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

	if(.not.associated(ptr3)) then
	  call perr(myname_,'.not.associated(template%"'//varn//'")')
	  if(.not.present(rc)) call die(myname_)
	  rc=GENERIC_ERROR
	  _EXIT_(myname_)
	  return
	endif

  l=size(ptr3,1)
  m=size(ptr3,2)
  n=size(ptr3,3)

  allocate( u_(l,m,n), v_(l,m,n),tv_(l,m,n),dp_(l,m,n), &
  	   qv_(l,m,n),ql_(l,n,n),qi_(l,m,n),o3_(l,m,n), stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allocate(inState_), stat =',ier)
	  call perr(myname_,'                 size(1) =',l)
	  call perr(myname_,'                 size(2) =',m)
	  call perr(myname_,'                 size(3) =',n)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  inState_created_=.true.
_EXIT_(myname_)
end subroutine inState_create_

subroutine inState_destroy_(rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use kinds, only: i_kind

  implicit none
  integer(i_kind),optional,intent(out):: rc

  character(len=*),parameter:: myname_=myname//"::inState_destroy_"
  character(len=*),parameter:: varn="DP"
  real(r_ESkind),pointer,dimension(:,:,:):: ptr3
  integer(i_kind):: ier,l,m,n

_ENTRY_(myname_)
  if(present(rc)) rc=ESMF_SUCCESS

  if(.not.inState_created_) then
    call perr(myname_,'inState_, not created, inState_created_ =', &
      inState_created_)
    if(.not.present(rc)) call die(myname_)
    rc=GENERIC_ERROR
    _EXIT_(myname_)
    return
  endif

  deallocate( u_, v_,tv_,dp_,qv_,ql_,qi_,o3_, stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate(inState_), stat =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  inState_created_=.false.
_EXIT_(myname_)
end subroutine inState_destroy_

subroutine inState_init_(v,rc)
  use ESMF, only: ESMF_SUCCESS
  use kinds, only: i_kind
  implicit none
  real(r_ESkind),intent(in):: v
  integer(i_kind),optional,intent(out):: rc

  character(len=*),parameter:: myname_=myname//"::inState_init_"

_ENTRY_(myname_)
  if(present(rc)) rc=ESMF_SUCCESS

  if(.not.inState_created_) then
    call perr(myname_,'inState_, not created, inState_created_ =', &
      inState_created_)
    if(.not.present(rc)) call die(myname_)
    rc=GENERIC_ERROR
    _EXIT_(myname_)
    return
  endif

   u_(:,:,:)=v
   v_(:,:,:)=v
  tv_(:,:,:)=v
  dp_(:,:,:)=v
  qv_(:,:,:)=v
  ql_(:,:,:)=v
  qi_(:,:,:)=v
  o3_(:,:,:)=v
_EXIT_(myname_)
end subroutine inState_init_

subroutine inState_getst_(state,rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use kinds, only:  i_kind
  implicit none
  type(ESMF_State),intent(inout):: state	! (in)
  integer(i_kind),optional,intent(out):: rc

  character(len=*),parameter:: myname_=myname//"::inState_getst_"
  integer:: ier

_ENTRY_(myname_)
  if(present(rc)) rc=ESMF_SUCCESS

  if(.not.inState_created_) then
    call perr(myname_,'inState_, not created, inState_created_ =', &
      inState_created_)
    if(.not.present(rc)) call die(myname_)
    rc=GENERIC_ERROR
    _EXIT_(myname_)
    return
  endif

  call get3d_( u_,state, u_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'get3d_("'// u_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call get3d_( v_,state, v_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'get3d_("'// v_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call get3d_(tv_,state,tv_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'get3d_("'//tv_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call get3d_(dp_,state,dp_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'get3d_("'//dp_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call get3d_(qv_,state,qv_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'get3d_("'//qv_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call get3d_(ql_,state,ql_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'get3d_("'//ql_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call get3d_(qi_,state,qi_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'get3d_("'//qi_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call get3d_(o3_,state,o3_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'get3d_("'//o3_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
_EXIT_(myname_)
end subroutine inState_getst_

subroutine inState_addst_(state,rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use kinds, only:  i_kind
  implicit none
  type(ESMF_State),intent(inout):: state	! (in)
  integer(i_kind),optional,intent(out):: rc

  character(len=*),parameter:: myname_=myname//"::inState_addst_"
  integer:: ier

_ENTRY_(myname_)
  if(present(rc)) rc=ESMF_SUCCESS

  if(.not.inState_created_) then
    call perr(myname_,'inState_, not created, inState_created_ =', &
      inState_created_)
    if(.not.present(rc)) call die(myname_)
    rc=GENERIC_ERROR
    _EXIT_(myname_)
    return
  endif

  call add3d_( u_,state, u_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'add3d_("'// u_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call add3d_( v_,state, v_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'add3d_("'// v_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call add3d_(tv_,state,tv_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'add3d_("'//tv_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call add3d_(dp_,state,dp_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'add3d_("'//dp_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call add3d_(qv_,state,qv_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'add3d_("'//qv_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call add3d_(ql_,state,ql_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'add3d_("'//ql_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call add3d_(qi_,state,qi_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'add3d_("'//qi_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call add3d_(o3_,state,o3_GCMpert,rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'add3d_("'//o3_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
_EXIT_(myname_)
end subroutine inState_addst_
subroutine inState_setst_(state,alloc,rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use kinds, only:  i_kind
  implicit none
  type(ESMF_State),intent(inout):: state
  logical,optional,intent(in):: alloc
  integer(i_kind),optional,intent(out):: rc

  character(len=*),parameter:: myname_=myname//"::inState_setst_"
  integer:: ier

_ENTRY_(myname_)
  if(present(rc)) rc=ESMF_SUCCESS

  if(.not.inState_created_) then
    call perr(myname_,'inState_, not created, inState_created_ =', &
      inState_created_)
    if(.not.present(rc)) call die(myname_)
    rc=GENERIC_ERROR
    _EXIT_(myname_)
    return
  endif

  call set3d_( u_,state, u_GCMpert,rc=ier,alloc=alloc)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'set3d_("'// u_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call set3d_( v_,state, v_GCMpert,rc=ier,alloc=alloc)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'set3d_("'// v_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call set3d_(tv_,state,tv_GCMpert,rc=ier,alloc=alloc)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'set3d_("'//tv_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call set3d_(dp_,state,dp_GCMpert,rc=ier,alloc=alloc)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'set3d_("'//dp_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call set3d_(qv_,state,qv_GCMpert,rc=ier,alloc=alloc)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'set3d_("'//qv_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call set3d_(ql_,state,ql_GCMpert,rc=ier,alloc=alloc)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'set3d_("'//ql_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call set3d_(qi_,state,qi_GCMpert,rc=ier,alloc=alloc)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'set3d_("'//qi_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
  call set3d_(o3_,state,o3_GCMpert,rc=ier,alloc=alloc)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'set3d_("'//o3_GCMpert//'")',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif
_EXIT_(myname_)
end subroutine inState_setst_

subroutine get3d_(x,state,vnm,rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  implicit none
  real(r_ESkind)  ,dimension(:,:,:),intent(out):: x
  type(ESMF_State),intent(inout):: state	! (in)
  character(len=*),intent(in   ):: vnm
  integer(i_kind) ,intent(out  ):: rc

  character(len=*),parameter:: myname_=myname//"::get3d_"
  real(r_ESkind),pointer,dimension(:,:,:):: ptr3_
_ENTRY_(myname_)

  ptr3_ => null()
  call ESMFL_StateGetPointerToData(state,ptr3_,vnm,rc=rc)
	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMFL_StateGetPointerToData("'//vnm//'"), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

	if(.not.associated(ptr3_)) then
	  call perr(myname_,'.not.associated(ptr => "'//vnm//'")')
	  rc=GENERIC_ERROR
	  _EXIT_(myname_)
	  return
	endif

  x(:,:,:) = ptr3_(:,:,:)
  nullify(ptr3_)
_EXIT_(myname_)
end subroutine get3d_
subroutine add3d_(x,state,vnm,rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  implicit none
  real(r_ESkind)  ,dimension(:,:,:),intent(inout):: x
  type(ESMF_State),intent(inout):: state	! (in)
  character(len=*),intent(in   ):: vnm
  integer(i_kind) ,intent(out  ):: rc

  character(len=*),parameter:: myname_=myname//"::add3d_"
  real(r_ESkind),pointer,dimension(:,:,:):: ptr3_
_ENTRY_(myname_)

  ptr3_ => null()
  call ESMFL_StateGetPointerToData(state,ptr3_,vnm,rc=rc)
	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMFL_StateGetPointerToData("'//vnm//'"), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

	if(.not.associated(ptr3_)) then
	  call perr(myname_,'.not.associated(ptr => "'//vnm//'")')
	  rc=GENERIC_ERROR
	  _EXIT_(myname_)
	  return
	endif

  x(:,:,:) = x(:,:,:)+ ptr3_(:,:,:)
  nullify(ptr3_)
_EXIT_(myname_)
end subroutine add3d_
subroutine set3d_(x,state,vnm,rc,alloc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  implicit none
  real(r_ESkind)  ,dimension(:,:,:),intent(in):: x
  type(ESMF_State),intent(inout):: state	! copy x state%vnm
  character(len=*),intent(in   ):: vnm
  integer(i_kind) ,intent(out  ):: rc
  logical,optional,intent(in   ):: alloc	! do alloc(), if not already so

  character(len=*),parameter:: myname_=myname//"::set3d_"
  real(r_ESkind),pointer,dimension(:,:,:):: ptr3_
  logical:: alloc_
_ENTRY_(myname_)
  alloc_=.false.
  if(present(alloc)) alloc_=alloc

  rc=ESMF_SUCCESS

  	! This get_pointer(vnm) will ensure the presence of pointer VNM
  ptr3_ => null()
  call ESMFL_StateGetPointerToData(state,ptr3_,vnm,rc=rc)
	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMFL_StateGetPointerToData("'//vnm//'"), rc =',rc)
	  _EXIT_(myname_)
	  return
	endif

  	! this is an implicit rule of GEOS-5, a variable in the export
	! state will not be returned, unless its storage in the export
	! state has been allocated by the user.
  if(associated(ptr3_)) then
    ptr3_(:,:,:) = x(:,:,:)
    nullify(ptr3_)
  endif
_EXIT_(myname_)
end subroutine set3d_

function phase_(adjoint,import,export)
  implicit none
  logical,optional,intent(in):: adjoint
  logical,optional,intent(in):: import
  logical,optional,intent(in):: export
  integer:: phase_

  logical:: adjoint_,import_,export_

  adjoint_=.false.; if(present(adjoint)) adjoint_=adjoint
  import_ =.true. ; if(present(import )) import_ =import
  export_ =.true. ; if(present(export )) export_ =export

  phase_=0
  if(    adjoint_) phase_=ior(phase_,1)
  if(.not.import_) phase_=ior(phase_,2)
  if(.not.export_) phase_=ior(phase_,4)
  phase_=phase_+1
end function phase_

function adjoint_(phase)
  implicit none
  integer,intent(in):: phase
  logical:: adjoint_
  adjoint_ = iand(phase-1,1)==1
end function adjoint_
function import_(phase)
  implicit none
  integer,intent(in):: phase
  logical:: import_
  import_  = iand(phase-1,2)==0
end function import_
function export_(phase)
  implicit none
  integer,intent(in):: phase
  logical:: export_
  export_  = iand(phase-1,4)==0
end function export_
end module stub_GEOSagcmPert_GridCompMod
