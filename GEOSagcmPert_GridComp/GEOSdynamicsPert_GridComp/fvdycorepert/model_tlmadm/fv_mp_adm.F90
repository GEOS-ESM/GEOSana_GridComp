module fv_mp_adm_mod

use fv_arrays_mod,   only : R_GRID
use fv_arrays_nlm_mod,only : fv_timing_onoff
use mpp_domains_mod, only : domain2D
use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
use mpp_domains_mod, only : mpp_group_update_initialized
use mpp_domains_mod, only : mpp_create_group_update,mpp_reset_group_update_field
use mpp_domains_mod, only : group_halo_update_type => mpp_group_update_type

#ifdef OLDMPP
use mpp_domains_mod, only : mpp_update_domains
#endif
use mpp_domains_mod, only : mpp_global_sum

use mpp_domains_mod, only : mpp_update_domains_ad, mpp_get_boundary_ad

use fv_mp_mod, only : XDir, YDir, ng
use fv_mp_mod, only : is, ie, js, je
use fv_mp_mod, only : isd, ied, jsd, jed

use fv_timing_mod, only: timing_on, timing_off

use tapenade_iter, only: pushcontrol, popcontrol, pushinteger, popinteger, &
                         pushrealarray, poprealarray, pushrealarray_adm, poprealarray_adm

implicit none
private

#include "mpif.h"

integer :: commglobal, ierror, npes

!fv_mp_mod routines
public complete_group_halo_update
public start_group_halo_update, start_group_halo_update_adm
public fill_corners_adm

!mpp_domains interface
public mpp_update_domains_adm, mpp_get_boundary_adm, mpp_global_sum_adm


! Regular fv_mp_mod routines
! --------------------------
interface start_group_halo_update
  module procedure start_var_group_update_2d
  module procedure start_var_group_update_3d
  module procedure start_var_group_update_4d
  module procedure start_vector_group_update_2d
  module procedure start_vector_group_update_3d
end interface start_group_halo_update

interface start_group_halo_update_adm
  module procedure start_var_group_update_2d_adm
  module procedure start_var_group_update_3d_adm
  module procedure start_var_group_update_4d_adm
  module procedure start_vector_group_update_2d_adm
  module procedure start_vector_group_update_3d_adm
end interface 

interface fill_corners_adm
    module procedure fill_corners_2d_r4_adm
    module procedure fill_corners_2d_r8_adm
    module procedure fill_corners_xy_2d_r4_adm
    module procedure fill_corners_xy_2d_r8_adm
end interface

interface fill_corners_agrid_adm
    module procedure fill_corners_agrid_r4_adm
    module procedure fill_corners_agrid_r8_adm
end interface

interface fill_corners_cgrid_adm
    module procedure fill_corners_cgrid_r4_adm
    module procedure fill_corners_cgrid_r8_adm
end interface

interface fill_corners_dgrid_adm
    module procedure fill_corners_dgrid_r4_adm
    module procedure fill_corners_dgrid_r8_adm
end interface


! These are invented interfaces to mpp_domains
! --------------------------------------------

interface mpp_global_sum_adm
    module procedure mpp_global_sum_2d_adm
end interface

interface mpp_update_domains_adm
    module procedure mpp_update_domain2d_2d_adm
    module procedure mpp_update_domain2d_3d_adm
    module procedure mpp_update_domain2d_4d_adm
    module procedure mpp_update_domain2d_5d_adm
    module procedure mpp_update_domain2d_2dv_adm
    module procedure mpp_update_domain2d_3dv_adm
    module procedure mpp_update_domain2d_4dv_adm
    module procedure mpp_update_domain2d_5dv_adm
end interface

interface mpp_get_boundary_adm
    module procedure mpp_get_boundary_2d_adm
    module procedure mpp_get_boundary_3d_adm
    module procedure mpp_get_boundary_2dv_adm
    module procedure mpp_get_boundary_3dv_adm
end interface

contains
 

! start_group_halo_update
! -----------------------

subroutine start_var_group_update_2d(group, &
#ifndef OLDMPP
                                     groupp, &
#endif
                                     array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)

  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real, dimension(:,:),         intent(inout) :: array
  type(domain2D),               intent(inout) :: domain
  integer,      optional,       intent(in)    :: flags
  integer,      optional,       intent(in)    :: position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete
  real                                        :: d_type
  logical                                     :: is_complete
! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      domain - contains domain information.
!  (in)      flags  - An optional integer indicating which directions the
!                       data should be sent.  
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (fv_timing_onoff) call timing_on('  FWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(array, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#else

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  else
    call mpp_create_group_update(group, array, domain, flags=flags, position=position, &
             whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete) then 
     call mpp_start_group_update(group, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  FWD_COMM_TOTAL')

end subroutine start_var_group_update_2d


subroutine start_var_group_update_3d(group, &
#ifndef OLDMPP
                                     groupp, &
#endif
                                     array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)

  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real, dimension(:,:,:),       intent(inout) :: array
  type(domain2D),               intent(inout) :: domain
  integer,           optional,  intent(in)    :: flags
  integer,           optional,  intent(in)    :: position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete
  real                                        :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      domain - contains domain information.
!  (in)      flags  - An optional integer indicating which directions the
!                       data should be sent.  
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (fv_timing_onoff) call timing_on('  FWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(array, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#else

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  else
    call mpp_create_group_update(group, array, domain, flags=flags, position=position, &
          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete) then
     call mpp_start_group_update(group, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  FWD_COMM_TOTAL')

end subroutine start_var_group_update_3d

subroutine start_var_group_update_4d(group, &
#ifndef OLDMPP
                                     groupp, &
#endif
                                     array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real, dimension(:,:,:,:),     intent(inout) :: array
  type(domain2D),               intent(inout) :: domain
  integer,           optional,  intent(in)    :: flags
  integer,           optional,  intent(in)    :: position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete
  real                                        :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      domain - contains domain information.
!  (in)      flags  - An optional integer indicating which directions the
!                       data should be sent.  
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  integer :: dirflag

  if (fv_timing_onoff) call timing_on('  FWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(array, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#else

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  else
    call mpp_create_group_update(group, array, domain, flags=flags, position=position, &
              whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete) then
     call mpp_start_group_update(group, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  FWD_COMM_TOTAL')

end subroutine start_var_group_update_4d



subroutine start_vector_group_update_2d(group, &
#ifndef OLDMPP
                                        groupp, &
#endif
                                        u_cmpt, v_cmpt, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real,       dimension(:,:),   intent(inout) :: u_cmpt, v_cmpt
  type(domain2d),               intent(inout) :: domain
  integer,            optional, intent(in)    :: flags
  integer,            optional, intent(in)    :: gridtype
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete
  real                                        :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged. 
!  (in)      domain - Contains domain decomposition information.
!  (in)      flags - An optional integer indicating which directions the
!                        data should be sent. 
!  (in)      gridtype - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      CGRID_NE or DGRID_NE, indicating where the two components of the
!                      vector are discretized. 
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (fv_timing_onoff) call timing_on('  FWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(u_cmpt,    v_cmpt,    domain, flags=flags, gridtype=gridtype, &
                                                   whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#else

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, domain, &
            flags=flags, gridtype=gridtype, &
            whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete) then
     call mpp_start_group_update(group, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  FWD_COMM_TOTAL')

end subroutine start_vector_group_update_2d

subroutine start_vector_group_update_3d(group, &
#ifndef OLDMPP
                                        groupp, &
#endif
                                        u_cmpt, v_cmpt, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real,       dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt
  type(domain2d),               intent(inout) :: domain
  integer,            optional, intent(in)    :: flags
  integer,            optional, intent(in)    :: gridtype
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete
  real                                        :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged. 
!  (in)      domain - Contains domain decomposition information.
!  (in)      flags - An optional integer indicating which directions the
!                        data should be sent. 
!  (in)      gridtype - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      CGRID_NE or DGRID_NE, indicating where the two components of the
!                      vector are discretized. 
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (fv_timing_onoff) call timing_on('  FWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(u_cmpt,    v_cmpt,    domain, flags=flags, gridtype=gridtype, &
                                                   whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#else

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, domain, &
            flags=flags, gridtype=gridtype, &
            whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete) then
     call mpp_start_group_update(group, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  FWD_COMM_TOTAL')

end subroutine start_vector_group_update_3d

! complete_group_halo_update
! --------------------------

subroutine complete_group_halo_update(group,&
#ifndef OLDMPP
                                      groupp, &
#endif
                                      domain)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  type(domain2d),               intent(inout) :: domain
  real                                        :: d_type

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!  (in)      domain - Contains domain decomposition information.

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

#ifndef OLDMPP

     call mpp_complete_group_update(group, domain, d_type)
     call mpp_complete_group_update(groupp, domain, d_type)

#endif

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine complete_group_halo_update


! start_group_halo_update_adm
! ---------------------------

 subroutine start_var_group_update_2d_adm(group, &
#ifndef OLDMPP
                                          groupp, &
#endif
                                          array, arrayp, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real, dimension(:,:),         intent(inout) :: array, arrayp
  real(8) :: array8(10,10)
  type(domain2D),               intent(inout) :: domain
  integer,      optional,       intent(in)    :: flags
  integer,      optional,       intent(in)    :: position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete
  real                                 :: d_type
  logical                                     :: is_complete
! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      domain - contains domain information.
!  (in)      flags  - An optional integer indicating which directions the
!                       data should be sent.  
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains_ad(arrayp, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#endif

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine start_var_group_update_2d_adm

subroutine start_var_group_update_3d_adm(group, &
#ifndef OLDMPP
                                         groupp, &
#endif
                                         array, arrayp, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real, dimension(:,:,:),       intent(inout) :: array, arrayp
  type(domain2D),               intent(inout) :: domain
  integer,           optional,  intent(in)    :: flags
  integer,           optional,  intent(in)    :: position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete 
  real                                        :: d_type
  logical                                     :: is_complete 

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      domain - contains domain information.
!  (in)      flags  - An optional integer indicating which directions the
!                       data should be sent.  
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains_ad(arrayp, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#endif

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine start_var_group_update_3d_adm

subroutine start_var_group_update_4d_adm(group, &
#ifndef OLDMPP
                                         groupp, &
#endif
                                         array, arrayp, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real, dimension(:,:,:,:),     intent(inout) :: array, arrayp
  type(domain2D),               intent(inout) :: domain
  integer,           optional,  intent(in)    :: flags
  integer,           optional,  intent(in)    :: position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete
  real                                        :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      domain - contains domain information.
!  (in)      flags  - An optional integer indicating which directions the
!                       data should be sent.  
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  integer :: dirflag

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains_ad(arrayp, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#endif

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine start_var_group_update_4d_adm



subroutine start_vector_group_update_2d_adm(group, &
#ifndef OLDMPP
                                            groupp, &
#endif
                                            u_cmpt, u_cmptp, v_cmpt, v_cmptp, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real,       dimension(:,:),   intent(inout) :: u_cmpt, v_cmpt, u_cmptp, v_cmptp
  type(domain2d),               intent(inout) :: domain
  integer,            optional, intent(in)    :: flags
  integer,            optional, intent(in)    :: gridtype
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete
  real                                 :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged. 
!  (in)      domain - Contains domain decomposition information.
!  (in)      flags - An optional integer indicating which directions the
!                        data should be sent. 
!  (in)      gridtype - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      CGRID_NE or DGRID_NE, indicating where the two components of the
!                      vector are discretized. 
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains_ad(u_cmptp, v_cmptp, domain, flags=flags, gridtype=gridtype, &
                                                  whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#endif

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine start_vector_group_update_2d_adm

subroutine start_vector_group_update_3d_adm(group, &
#ifndef OLDMPP
                                            groupp, &
#endif
                                            u_cmpt, u_cmptp, v_cmpt, v_cmptp, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: groupp
#endif
  real,       dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt, u_cmptp, v_cmptp
  type(domain2d),               intent(inout) :: domain
  integer,            optional, intent(in)    :: flags
  integer,            optional, intent(in)    :: gridtype
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete
  real                                 :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged. 
!  (in)      domain - Contains domain decomposition information.
!  (in)      flags - An optional integer indicating which directions the
!                        data should be sent. 
!  (in)      gridtype - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      CGRID_NE or DGRID_NE, indicating where the two components of the
!                      vector are discretized. 
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains_ad(u_cmptp, v_cmptp, domain, flags=flags, gridtype=gridtype, &
                                                  whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

#endif

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine start_vector_group_update_3d_adm



! fill_corners_adm
! ----------------

  SUBROUTINE FILL_CORNERS_2D_R4_ADM(q, q_ad, npx, npy, fill, agrid, &
&   bgrid)
    IMPLICIT NONE
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: q
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: q_ad
    INTEGER, INTENT(IN) :: npx, npy
! X-Dir or Y-Dir 
    INTEGER, INTENT(IN) :: fill
    LOGICAL, OPTIONAL, INTENT(IN) :: agrid, bgrid
    INTEGER :: i, j
    INTRINSIC PRESENT
    REAL(kind=4) :: tmp
    REAL(kind=4) :: tmp_ad
    REAL(kind=4) :: tmp0
    REAL(kind=4) :: tmp_ad0
    REAL(kind=4) :: tmp1
    REAL(kind=8) :: tmp_ad1
    REAL(kind=8) :: tmp2
    REAL(kind=8) :: tmp_ad2
    REAL(kind=4) :: tmp3
    REAL(kind=4) :: tmp_ad3
    REAL(kind=4) :: tmp4
    REAL(kind=4) :: tmp_ad4
    REAL(kind=4) :: tmp5
    REAL(kind=4) :: tmp_ad5
    REAL(kind=4) :: tmp6
    REAL(kind=4) :: tmp_ad6
    REAL(kind=4) :: tmp7
    REAL(kind=4) :: tmp_ad7
    REAL(kind=4) :: tmp8
    REAL(kind=4) :: tmp_ad8
    REAL(kind=4) :: tmp9
    REAL(kind=4) :: tmp_ad9
    REAL(kind=4) :: tmp10
    REAL(kind=4) :: tmp_ad10
    REAL(kind=4) :: tmp11
    REAL(kind=4) :: tmp_ad11
    REAL(kind=4) :: tmp12
    REAL(kind=4) :: tmp_ad12
    REAL(kind=4) :: tmp13
    REAL(kind=4) :: tmp_ad13
    REAL(kind=4) :: tmp14
    REAL(kind=4) :: tmp_ad14
    REAL(kind=4) :: tmp15
    REAL(kind=4) :: tmp_ad15
    REAL(kind=4) :: tmp16
    REAL(kind=4) :: tmp_ad16
    REAL(kind=4) :: tmp17
    REAL(kind=4) :: tmp_ad17
    REAL(kind=4) :: tmp18
    REAL(kind=4) :: tmp_ad18
    REAL(kind=4) :: tmp19
    REAL(kind=4) :: tmp_ad19
    REAL(kind=4) :: tmp20
    REAL(kind=4) :: tmp_ad20
    REAL(kind=4) :: tmp21
    REAL(kind=4) :: tmp_ad21
    REAL(kind=4) :: tmp22
    REAL(kind=4) :: tmp_ad22
    INTEGER :: branch
    IF (PRESENT(bgrid)) THEN
      IF (bgrid) THEN
        SELECT CASE  (fill) 
        CASE (xdir) 
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad2 = q_ad(npx+i, npy+j)
                q_ad(npx+i, npy+j) = 0.0
                q_ad(npx+j, npy-i) = q_ad(npx+j, npy-i) + tmp_ad2
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad1 = q_ad(npx+i, 1-j)
                q_ad(npx+i, 1-j) = 0.0
                q_ad(npx+j, i+1) = q_ad(npx+j, i+1) + tmp_ad1
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad0 = q_ad(1-i, npy+j)
                q_ad(1-i, npy+j) = 0.0
                q_ad(1-j, npy-i) = q_ad(1-j, npy-i) + tmp_ad0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad = q_ad(1-i, 1-j)
                q_ad(1-i, 1-j) = 0.0
                q_ad(1-j, i+1) = q_ad(1-j, i+1) + tmp_ad
              END IF
            END DO
          END DO
        CASE (ydir) 
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad6 = q_ad(npx+j, npy+i)
                q_ad(npx+j, npy+i) = 0.0
                q_ad(npx-i, npy+j) = q_ad(npx-i, npy+j) + tmp_ad6
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad5 = q_ad(npx+j, 1-i)
                q_ad(npx+j, 1-i) = 0.0
                q_ad(npx-i, 1-j) = q_ad(npx-i, 1-j) + tmp_ad5
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad4 = q_ad(1-j, npy+i)
                q_ad(1-j, npy+i) = 0.0
                q_ad(i+1, npy+j) = q_ad(i+1, npy+j) + tmp_ad4
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad3 = q_ad(1-j, 1-i)
                q_ad(1-j, 1-i) = 0.0
                q_ad(i+1, 1-j) = q_ad(i+1, 1-j) + tmp_ad3
              END IF
            END DO
          END DO
        CASE DEFAULT
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad10 = q_ad(npx+i, npy+j)
                q_ad(npx+i, npy+j) = 0.0
                q_ad(npx+j, npy-i) = q_ad(npx+j, npy-i) + tmp_ad10
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad9 = q_ad(npx+i, 1-j)
                q_ad(npx+i, 1-j) = 0.0
                q_ad(npx+j, i+1) = q_ad(npx+j, i+1) + tmp_ad9
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad8 = q_ad(1-i, npy+j)
                q_ad(1-i, npy+j) = 0.0
                q_ad(1-j, npy-i) = q_ad(1-j, npy-i) + tmp_ad8
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad7 = q_ad(1-i, 1-j)
                q_ad(1-i, 1-j) = 0.0
                q_ad(1-j, i+1) = q_ad(1-j, i+1) + tmp_ad7
              END IF
            END DO
          END DO
        END SELECT
      END IF
    ELSE IF (PRESENT(agrid)) THEN
      IF (agrid) THEN
        SELECT CASE  (fill) 
        CASE (xdir) 
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad14 = q_ad(npx-1+i, npy-1+j)
                q_ad(npx-1+i, npy-1+j) = 0.0
                q_ad(npx-1+j, npy-1-i+1) = q_ad(npx-1+j, npy-1-i+1) + &
&                 tmp_ad14
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad13 = q_ad(npx-1+i, 1-j)
                q_ad(npx-1+i, 1-j) = 0.0
                q_ad(npx-1+j, i) = q_ad(npx-1+j, i) + tmp_ad13
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad12 = q_ad(1-i, npy-1+j)
                q_ad(1-i, npy-1+j) = 0.0
                q_ad(1-j, npy-1-i+1) = q_ad(1-j, npy-1-i+1) + tmp_ad12
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad11 = q_ad(1-i, 1-j)
                q_ad(1-i, 1-j) = 0.0
                q_ad(1-j, i) = q_ad(1-j, i) + tmp_ad11
              END IF
            END DO
          END DO
        CASE (ydir) 
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad18 = q_ad(npx-1+j, npy-1+i)
                q_ad(npx-1+j, npy-1+i) = 0.0
                q_ad(npx-1-i+1, npy-1+j) = q_ad(npx-1-i+1, npy-1+j) + &
&                 tmp_ad18
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad17 = q_ad(npx-1+j, 1-i)
                q_ad(npx-1+j, 1-i) = 0.0
                q_ad(npx-1-i+1, 1-j) = q_ad(npx-1-i+1, 1-j) + tmp_ad17
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad16 = q_ad(1-j, npy-1+i)
                q_ad(1-j, npy-1+i) = 0.0
                q_ad(i, npy-1+j) = q_ad(i, npy-1+j) + tmp_ad16
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad15 = q_ad(1-j, 1-i)
                q_ad(1-j, 1-i) = 0.0
                q_ad(i, 1-j) = q_ad(i, 1-j) + tmp_ad15
              END IF
            END DO
          END DO
        CASE DEFAULT
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad22 = q_ad(npx-1+j, npy-1+i)
                q_ad(npx-1+j, npy-1+i) = 0.0
                q_ad(npx-1-i+1, npy-1+j) = q_ad(npx-1-i+1, npy-1+j) + &
&                 tmp_ad22
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad21 = q_ad(npx-1+j, 1-i)
                q_ad(npx-1+j, 1-i) = 0.0
                q_ad(npx-1-i+1, 1-j) = q_ad(npx-1-i+1, 1-j) + tmp_ad21
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad20 = q_ad(1-j, npy-1+i)
                q_ad(1-j, npy-1+i) = 0.0
                q_ad(i, npy-1+j) = q_ad(i, npy-1+j) + tmp_ad20
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad19 = q_ad(1-j, 1-i)
                q_ad(1-j, 1-i) = 0.0
                q_ad(i, 1-j) = q_ad(i, 1-j) + tmp_ad19
              END IF
            END DO
          END DO
        END SELECT
      END IF
    END IF
  END SUBROUTINE FILL_CORNERS_2D_R4_ADM

  SUBROUTINE FILL_CORNERS_2D_R8_ADM(q, q_ad, npx, npy, fill, agrid, &
&   bgrid)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: q
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: q_ad
    INTEGER, INTENT(IN) :: npx, npy
! X-Dir or Y-Dir 
    INTEGER, INTENT(IN) :: fill
    LOGICAL, OPTIONAL, INTENT(IN) :: agrid, bgrid
    INTEGER :: i, j
    INTRINSIC PRESENT
    REAL(kind=8) :: tmp
    REAL(kind=8) :: tmp_ad
    REAL(kind=8) :: tmp0
    REAL(kind=8) :: tmp_ad0
    REAL(kind=8) :: tmp1
    REAL(kind=8) :: tmp_ad1
    REAL(kind=8) :: tmp2
    REAL(kind=8) :: tmp_ad2
    REAL(kind=8) :: tmp3
    REAL(kind=8) :: tmp_ad3
    REAL(kind=8) :: tmp4
    REAL(kind=8) :: tmp_ad4
    REAL(kind=8) :: tmp5
    REAL(kind=8) :: tmp_ad5
    REAL(kind=8) :: tmp6
    REAL(kind=8) :: tmp_ad6
    REAL(kind=8) :: tmp7
    REAL(kind=8) :: tmp_ad7
    REAL(kind=8) :: tmp8
    REAL(kind=8) :: tmp_ad8
    REAL(kind=8) :: tmp9
    REAL(kind=8) :: tmp_ad9
    REAL(kind=8) :: tmp10
    REAL(kind=8) :: tmp_ad10
    REAL(kind=8) :: tmp11
    REAL(kind=8) :: tmp_ad11
    REAL(kind=8) :: tmp12
    REAL(kind=8) :: tmp_ad12
    REAL(kind=8) :: tmp13
    REAL(kind=8) :: tmp_ad13
    REAL(kind=8) :: tmp14
    REAL(kind=8) :: tmp_ad14
    REAL(kind=8) :: tmp15
    REAL(kind=8) :: tmp_ad15
    REAL(kind=8) :: tmp16
    REAL(kind=8) :: tmp_ad16
    REAL(kind=8) :: tmp17
    REAL(kind=8) :: tmp_ad17
    REAL(kind=8) :: tmp18
    REAL(kind=8) :: tmp_ad18
    REAL(kind=8) :: tmp19
    REAL(kind=8) :: tmp_ad19
    REAL(kind=8) :: tmp20
    REAL(kind=8) :: tmp_ad20
    REAL(kind=8) :: tmp21
    REAL(kind=8) :: tmp_ad21
    REAL(kind=8) :: tmp22
    REAL(kind=8) :: tmp_ad22
    INTEGER :: branch
    IF (PRESENT(bgrid)) THEN
      IF (bgrid) THEN
        SELECT CASE  (fill) 
        CASE (xdir) 
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad2 = q_ad(npx+i, npy+j)
                q_ad(npx+i, npy+j) = 0.0
                q_ad(npx+j, npy-i) = q_ad(npx+j, npy-i) + tmp_ad2
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad1 = q_ad(npx+i, 1-j)
                q_ad(npx+i, 1-j) = 0.0
                q_ad(npx+j, i+1) = q_ad(npx+j, i+1) + tmp_ad1
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad0 = q_ad(1-i, npy+j)
                q_ad(1-i, npy+j) = 0.0
                q_ad(1-j, npy-i) = q_ad(1-j, npy-i) + tmp_ad0
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad = q_ad(1-i, 1-j)
                q_ad(1-i, 1-j) = 0.0
                q_ad(1-j, i+1) = q_ad(1-j, i+1) + tmp_ad
              END IF
            END DO
          END DO
        CASE (ydir) 
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad6 = q_ad(npx+j, npy+i)
                q_ad(npx+j, npy+i) = 0.0
                q_ad(npx-i, npy+j) = q_ad(npx-i, npy+j) + tmp_ad6
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad5 = q_ad(npx+j, 1-i)
                q_ad(npx+j, 1-i) = 0.0
                q_ad(npx-i, 1-j) = q_ad(npx-i, 1-j) + tmp_ad5
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad4 = q_ad(1-j, npy+i)
                q_ad(1-j, npy+i) = 0.0
                q_ad(i+1, npy+j) = q_ad(i+1, npy+j) + tmp_ad4
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad3 = q_ad(1-j, 1-i)
                q_ad(1-j, 1-i) = 0.0
                q_ad(i+1, 1-j) = q_ad(i+1, 1-j) + tmp_ad3
              END IF
            END DO
          END DO
        CASE DEFAULT
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad10 = q_ad(npx+i, npy+j)
                q_ad(npx+i, npy+j) = 0.0
                q_ad(npx+j, npy-i) = q_ad(npx+j, npy-i) + tmp_ad10
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad9 = q_ad(npx+i, 1-j)
                q_ad(npx+i, 1-j) = 0.0
                q_ad(npx+j, i+1) = q_ad(npx+j, i+1) + tmp_ad9
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad8 = q_ad(1-i, npy+j)
                q_ad(1-i, npy+j) = 0.0
                q_ad(1-j, npy-i) = q_ad(1-j, npy-i) + tmp_ad8
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad7 = q_ad(1-i, 1-j)
                q_ad(1-i, 1-j) = 0.0
                q_ad(1-j, i+1) = q_ad(1-j, i+1) + tmp_ad7
              END IF
            END DO
          END DO
        END SELECT
      END IF
    ELSE IF (PRESENT(agrid)) THEN
      IF (agrid) THEN
        SELECT CASE  (fill) 
        CASE (xdir) 
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad14 = q_ad(npx-1+i, npy-1+j)
                q_ad(npx-1+i, npy-1+j) = 0.0
                q_ad(npx-1+j, npy-1-i+1) = q_ad(npx-1+j, npy-1-i+1) + &
&                 tmp_ad14
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad13 = q_ad(npx-1+i, 1-j)
                q_ad(npx-1+i, 1-j) = 0.0
                q_ad(npx-1+j, i) = q_ad(npx-1+j, i) + tmp_ad13
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad12 = q_ad(1-i, npy-1+j)
                q_ad(1-i, npy-1+j) = 0.0
                q_ad(1-j, npy-1-i+1) = q_ad(1-j, npy-1-i+1) + tmp_ad12
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad11 = q_ad(1-i, 1-j)
                q_ad(1-i, 1-j) = 0.0
                q_ad(1-j, i) = q_ad(1-j, i) + tmp_ad11
              END IF
            END DO
          END DO
        CASE (ydir) 
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad18 = q_ad(npx-1+j, npy-1+i)
                q_ad(npx-1+j, npy-1+i) = 0.0
                q_ad(npx-1-i+1, npy-1+j) = q_ad(npx-1-i+1, npy-1+j) + &
&                 tmp_ad18
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad17 = q_ad(npx-1+j, 1-i)
                q_ad(npx-1+j, 1-i) = 0.0
                q_ad(npx-1-i+1, 1-j) = q_ad(npx-1-i+1, 1-j) + tmp_ad17
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad16 = q_ad(1-j, npy-1+i)
                q_ad(1-j, npy-1+i) = 0.0
                q_ad(i, npy-1+j) = q_ad(i, npy-1+j) + tmp_ad16
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad15 = q_ad(1-j, 1-i)
                q_ad(1-j, 1-i) = 0.0
                q_ad(i, 1-j) = q_ad(i, 1-j) + tmp_ad15
              END IF
            END DO
          END DO
        CASE DEFAULT
          DO j=1,ng
            DO i=1,ng
!SW Corner 
              IF (is .EQ. 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NW Corner
              IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!SE Corner
              IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
!NE Corner
              IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          DO j=ng,1,-1
            DO i=ng,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                tmp_ad22 = q_ad(npx-1+j, npy-1+i)
                q_ad(npx-1+j, npy-1+i) = 0.0
                q_ad(npx-1-i+1, npy-1+j) = q_ad(npx-1-i+1, npy-1+j) + &
&                 tmp_ad22
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad21 = q_ad(npx-1+j, 1-i)
                q_ad(npx-1+j, 1-i) = 0.0
                q_ad(npx-1-i+1, 1-j) = q_ad(npx-1-i+1, 1-j) + tmp_ad21
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad20 = q_ad(1-j, npy-1+i)
                q_ad(1-j, npy-1+i) = 0.0
                q_ad(i, npy-1+j) = q_ad(i, npy-1+j) + tmp_ad20
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                tmp_ad19 = q_ad(1-j, 1-i)
                q_ad(1-j, 1-i) = 0.0
                q_ad(i, 1-j) = q_ad(i, 1-j) + tmp_ad19
              END IF
            END DO
          END DO
        END SELECT
      END IF
    END IF
  END SUBROUTINE FILL_CORNERS_2D_R8_ADM

  SUBROUTINE FILL_CORNERS_XY_2D_R4_ADM(x, x_ad, y, y_ad, npx, npy, dgrid&
&   , agrid, cgrid, vector)
    IMPLICIT NONE
!(isd:ied  ,jsd:jed+1)
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x_ad
!(isd:ied+1,jsd:jed  )
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y_ad
    INTEGER, INTENT(IN) :: npx, npy
    LOGICAL, OPTIONAL, INTENT(IN) :: dgrid, agrid, cgrid, vector
    INTEGER :: i, j
    REAL(kind=4) :: mysign
    INTRINSIC PRESENT
    INTEGER :: branch
    mysign = 1.0
    IF (PRESENT(vector)) THEN
      IF (vector) THEN
        CALL PUSHCONTROL1B(0)
        mysign = -1.0
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (PRESENT(dgrid)) THEN
      CALL FILL_CORNERS_DGRID_ADM(x, x_ad, y, y_ad, npx, npy, mysign)
    ELSE IF (PRESENT(cgrid)) THEN
      CALL FILL_CORNERS_CGRID_ADM(x, x_ad, y, y_ad, npx, npy, mysign)
    ELSE IF (PRESENT(agrid)) THEN
      CALL FILL_CORNERS_AGRID_ADM(x, x_ad, y, y_ad, npx, npy, mysign)
    ELSE
      CALL FILL_CORNERS_AGRID_ADM(x, x_ad, y, y_ad, npx, npy, mysign)
    END IF
    CALL POPCONTROL1B(branch)
  END SUBROUTINE FILL_CORNERS_XY_2D_R4_ADM

  SUBROUTINE FILL_CORNERS_XY_2D_R8_ADM(x, x_ad, y, y_ad, npx, npy, dgrid&
&   , agrid, cgrid, vector)
    IMPLICIT NONE
!(isd:ied  ,jsd:jed+1)
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x_ad
!(isd:ied+1,jsd:jed  )
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y_ad
    INTEGER, INTENT(IN) :: npx, npy
    LOGICAL, OPTIONAL, INTENT(IN) :: dgrid, agrid, cgrid, vector
    INTEGER :: i, j
    REAL(kind=8) :: mysign
    INTRINSIC PRESENT
    INTEGER :: branch
    mysign = 1.0
    IF (PRESENT(vector)) THEN
      IF (vector) THEN
        CALL PUSHCONTROL1B(0)
        mysign = -1.0
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (PRESENT(dgrid)) THEN
      CALL FILL_CORNERS_DGRID_ADM(x, x_ad, y, y_ad, npx, npy, mysign)
    ELSE IF (PRESENT(cgrid)) THEN
      CALL FILL_CORNERS_CGRID_ADM(x, x_ad, y, y_ad, npx, npy, mysign)
    ELSE IF (PRESENT(agrid)) THEN
      CALL FILL_CORNERS_AGRID_ADM(x, x_ad, y, y_ad, npx, npy, mysign)
    ELSE
      CALL FILL_CORNERS_AGRID_ADM(x, x_ad, y, y_ad, npx, npy, mysign)
    END IF
    CALL POPCONTROL1B(branch)
  END SUBROUTINE FILL_CORNERS_XY_2D_R8_ADM

! fill_corners_agrid_adm
! ----------------------

  SUBROUTINE FILL_CORNERS_AGRID_R4_ADM(x, x_ad, y, y_ad, npx, npy, &
&   mysign)
    IMPLICIT NONE
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x_ad
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y_ad
    INTEGER, INTENT(IN) :: npx, npy
    REAL(kind=4), INTENT(IN) :: mysign
    INTEGER :: i, j
    INTEGER :: branch
    DO j=1,ng
      DO i=1,ng
!SW Corner
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=1,ng
      DO i=1,ng
!SW Corner
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          x_ad(npx-1-i+1, npy-1+j) = x_ad(npx-1-i+1, npy-1+j) + mysign*&
&           y_ad(npx-1+j, npy-1+i)
          y_ad(npx-1+j, npy-1+i) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(npx-1-i+1, 1-j) = x_ad(npx-1-i+1, 1-j) + y_ad(npx-1+j, 1-&
&           i)
          y_ad(npx-1+j, 1-i) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(i, npy-1+j) = x_ad(i, npy-1+j) + y_ad(1-j, npy-1+i)
          y_ad(1-j, npy-1+i) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(i, 1-j) = x_ad(i, 1-j) + mysign*y_ad(1-j, 1-i)
          y_ad(1-j, 1-i) = 0.0
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          y_ad(npx-1+j, npy-1-i+1) = y_ad(npx-1+j, npy-1-i+1) + mysign*&
&           x_ad(npx-1+i, npy-1+j)
          x_ad(npx-1+i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(npx-1+j, i) = y_ad(npx-1+j, i) + x_ad(npx-1+i, 1-j)
          x_ad(npx-1+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(1-j, npy-1-i+1) = y_ad(1-j, npy-1-i+1) + x_ad(1-i, npy-1+&
&           j)
          x_ad(1-i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(1-j, i) = y_ad(1-j, i) + mysign*x_ad(1-i, 1-j)
          x_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
  END SUBROUTINE FILL_CORNERS_AGRID_R4_ADM

  SUBROUTINE FILL_CORNERS_AGRID_R8_ADM(x, x_ad, y, y_ad, npx, npy, &
&   mysign)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x_ad
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y_ad
    INTEGER, INTENT(IN) :: npx, npy
    REAL(kind=8), INTENT(IN) :: mysign
    INTEGER :: i, j
    INTEGER :: branch
    DO j=1,ng
      DO i=1,ng
!SW Corner
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=1,ng
      DO i=1,ng
!SW Corner
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je .EQ. npy - 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie .EQ. npx - 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie .EQ. npx - 1 .AND. je .EQ. npy - 1) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          x_ad(npx-1-i+1, npy-1+j) = x_ad(npx-1-i+1, npy-1+j) + mysign*&
&           y_ad(npx-1+j, npy-1+i)
          y_ad(npx-1+j, npy-1+i) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(npx-1-i+1, 1-j) = x_ad(npx-1-i+1, 1-j) + y_ad(npx-1+j, 1-&
&           i)
          y_ad(npx-1+j, 1-i) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(i, npy-1+j) = x_ad(i, npy-1+j) + y_ad(1-j, npy-1+i)
          y_ad(1-j, npy-1+i) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(i, 1-j) = x_ad(i, 1-j) + mysign*y_ad(1-j, 1-i)
          y_ad(1-j, 1-i) = 0.0
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          y_ad(npx-1+j, npy-1-i+1) = y_ad(npx-1+j, npy-1-i+1) + mysign*&
&           x_ad(npx-1+i, npy-1+j)
          x_ad(npx-1+i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(npx-1+j, i) = y_ad(npx-1+j, i) + x_ad(npx-1+i, 1-j)
          x_ad(npx-1+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(1-j, npy-1-i+1) = y_ad(1-j, npy-1-i+1) + x_ad(1-i, npy-1+&
&           j)
          x_ad(1-i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(1-j, i) = y_ad(1-j, i) + mysign*x_ad(1-i, 1-j)
          x_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
  END SUBROUTINE FILL_CORNERS_AGRID_R8_ADM

! fill_corners_cgrid_adm
! ----------------------

  SUBROUTINE FILL_CORNERS_CGRID_R4_ADM(x, x_ad, y, y_ad, npx, npy, &
&   mysign)
    IMPLICIT NONE
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x_ad
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y_ad
    INTEGER, INTENT(IN) :: npx, npy
    REAL(kind=4), INTENT(IN) :: mysign
    INTEGER :: i, j
    INTEGER :: branch
    DO j=1,ng
      DO i=1,ng
!SW Corner 
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie + 1 .EQ. npx .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie + 1 .EQ. npx .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=1,ng
      DO i=1,ng
!SW Corner 
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie + 1 .EQ. npx .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie + 1 .EQ. npx .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          x_ad(npx+j, npy-i) = x_ad(npx+j, npy-i) + y_ad(npx-1+i, npy+j)
          y_ad(npx-1+i, npy+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(npx+j, i) = x_ad(npx+j, i) + mysign*y_ad(npx-1+i, 1-j)
          y_ad(npx-1+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(1-j, npy-i) = x_ad(1-j, npy-i) + mysign*y_ad(1-i, npy+j)
          y_ad(1-i, npy+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(1-j, i) = x_ad(1-j, i) + y_ad(1-i, 1-j)
          y_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          y_ad(npx-j, npy+i) = y_ad(npx-j, npy+i) + x_ad(npx+i, npy-1+j)
          x_ad(npx+i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(npx-j, 1-i) = y_ad(npx-j, 1-i) + mysign*x_ad(npx+i, 1-j)
          x_ad(npx+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(j, npy+i) = y_ad(j, npy+i) + mysign*x_ad(1-i, npy-1+j)
          x_ad(1-i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(j, 1-i) = y_ad(j, 1-i) + x_ad(1-i, 1-j)
          x_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
  END SUBROUTINE FILL_CORNERS_CGRID_R4_ADM

  SUBROUTINE FILL_CORNERS_CGRID_R8_ADM(x, x_ad, y, y_ad, npx, npy, &
&   mysign)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x_ad
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y_ad
    INTEGER, INTENT(IN) :: npx, npy
    REAL(kind=8), INTENT(IN) :: mysign
    INTEGER :: i, j
    INTEGER :: branch
    DO j=1,ng
      DO i=1,ng
!SW Corner 
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie + 1 .EQ. npx .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie + 1 .EQ. npx .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=1,ng
      DO i=1,ng
!SW Corner 
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie + 1 .EQ. npx .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie + 1 .EQ. npx .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          x_ad(npx+j, npy-i) = x_ad(npx+j, npy-i) + y_ad(npx-1+i, npy+j)
          y_ad(npx-1+i, npy+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(npx+j, i) = x_ad(npx+j, i) + mysign*y_ad(npx-1+i, 1-j)
          y_ad(npx-1+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(1-j, npy-i) = x_ad(1-j, npy-i) + mysign*y_ad(1-i, npy+j)
          y_ad(1-i, npy+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(1-j, i) = x_ad(1-j, i) + y_ad(1-i, 1-j)
          y_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          y_ad(npx-j, npy+i) = y_ad(npx-j, npy+i) + x_ad(npx+i, npy-1+j)
          x_ad(npx+i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(npx-j, 1-i) = y_ad(npx-j, 1-i) + mysign*x_ad(npx+i, 1-j)
          x_ad(npx+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(j, npy+i) = y_ad(j, npy+i) + mysign*x_ad(1-i, npy-1+j)
          x_ad(1-i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(j, 1-i) = y_ad(j, 1-i) + x_ad(1-i, 1-j)
          x_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
  END SUBROUTINE FILL_CORNERS_CGRID_R8_ADM

! fill_corners_dgrid_adm
! ----------------------

  SUBROUTINE FILL_CORNERS_DGRID_R4_ADM(x, x_ad, y, y_ad, npx, npy, &
&   mysign)
    IMPLICIT NONE
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x_ad
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y
    REAL(kind=4), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y_ad
    INTEGER, INTENT(IN) :: npx, npy
    REAL(kind=4), INTENT(IN) :: mysign
    INTEGER :: i, j
    INTEGER :: branch
    DO j=1,ng
      DO i=1,ng
!   if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) =        y(j+1  ,1-i    )  !SW Corner 
!   if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) = mySign*y(j+1  ,npy-1+i)  !NW Corner
!   if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) = mySign*y(npx-j,1-i    )  !SE Corner
!   if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) =        y(npx-j,npy-1+i)  !NE Corner
!SW Corner 
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie + 1 .EQ. npx .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie + 1 .EQ. npx .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=1,ng
      DO i=1,ng
!  if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) =        x(1-j    ,i+1  )  !SW Corner 
!  if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) = mySign*x(1-j    ,npy-i)  !NW Corner
!  if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) = mySign*x(npx-1+j,i+1  )  !SE Corner
!  if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) =        x(npx-1+j,npy-i)  !NE Corner
!SW Corner 
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie + 1 .EQ. npx .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie + 1 .EQ. npx .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          x_ad(npx-j, npy+i) = x_ad(npx-j, npy+i) + mysign*y_ad(npx+i, &
&           npy-1+j)
          y_ad(npx+i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(npx-j, 1-i) = x_ad(npx-j, 1-i) + y_ad(npx+i, 1-j)
          y_ad(npx+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(j, npy+i) = x_ad(j, npy+i) + y_ad(1-i, npy-1+j)
          y_ad(1-i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(j, 1-i) = x_ad(j, 1-i) + mysign*y_ad(1-i, 1-j)
          y_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          y_ad(npx+j, npy-i) = y_ad(npx+j, npy-i) + mysign*x_ad(npx-1+i&
&           , npy+j)
          x_ad(npx-1+i, npy+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(npx+j, i) = y_ad(npx+j, i) + x_ad(npx-1+i, 1-j)
          x_ad(npx-1+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(1-j, npy-i) = y_ad(1-j, npy-i) + x_ad(1-i, npy+j)
          x_ad(1-i, npy+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(1-j, i) = y_ad(1-j, i) + mysign*x_ad(1-i, 1-j)
          x_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
  END SUBROUTINE FILL_CORNERS_DGRID_R4_ADM

  SUBROUTINE FILL_CORNERS_DGRID_R8_ADM(x, x_ad, y, y_ad, npx, npy, &
&   mysign)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: x_ad
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y
    REAL(kind=8), DIMENSION(isd:, jsd:), INTENT(INOUT) :: y_ad
    INTEGER, INTENT(IN) :: npx, npy
    REAL(kind=8), INTENT(IN) :: mysign
    INTEGER :: i, j
    INTEGER :: branch
    DO j=1,ng
      DO i=1,ng
!   if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) =        y(j+1  ,1-i    )  !SW Corner 
!   if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) = mySign*y(j+1  ,npy-1+i)  !NW Corner
!   if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) = mySign*y(npx-j,1-i    )  !SE Corner
!   if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) =        y(npx-j,npy-1+i)  !NE Corner
!SW Corner 
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie + 1 .EQ. npx .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie + 1 .EQ. npx .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=1,ng
      DO i=1,ng
!  if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) =        x(1-j    ,i+1  )  !SW Corner 
!  if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) = mySign*x(1-j    ,npy-i)  !NW Corner
!  if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) = mySign*x(npx-1+j,i+1  )  !SE Corner
!  if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) =        x(npx-1+j,npy-i)  !NE Corner
!SW Corner 
        IF (is .EQ. 1 .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NW Corner
        IF (is .EQ. 1 .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!SE Corner
        IF (ie + 1 .EQ. npx .AND. js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!NE Corner
        IF (ie + 1 .EQ. npx .AND. je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          x_ad(npx-j, npy+i) = x_ad(npx-j, npy+i) + mysign*y_ad(npx+i, &
&           npy-1+j)
          y_ad(npx+i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(npx-j, 1-i) = x_ad(npx-j, 1-i) + y_ad(npx+i, 1-j)
          y_ad(npx+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(j, npy+i) = x_ad(j, npy+i) + y_ad(1-i, npy-1+j)
          y_ad(1-i, npy-1+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          x_ad(j, 1-i) = x_ad(j, 1-i) + mysign*y_ad(1-i, 1-j)
          y_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
    DO j=ng,1,-1
      DO i=ng,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          y_ad(npx+j, npy-i) = y_ad(npx+j, npy-i) + mysign*x_ad(npx-1+i&
&           , npy+j)
          x_ad(npx-1+i, npy+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(npx+j, i) = y_ad(npx+j, i) + x_ad(npx-1+i, 1-j)
          x_ad(npx-1+i, 1-j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(1-j, npy-i) = y_ad(1-j, npy-i) + x_ad(1-i, npy+j)
          x_ad(1-i, npy+j) = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          y_ad(1-j, i) = y_ad(1-j, i) + mysign*x_ad(1-i, 1-j)
          x_ad(1-i, 1-j) = 0.0
        END IF
      END DO
    END DO
  END SUBROUTINE FILL_CORNERS_DGRID_R8_ADM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MPP_DOMAINS INTERFACES THAT ARE NOT NORMALLY IN FV_MP_MOD !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! mpp_global_sum_adm
! ------------------

  real(kind=r_grid) function mpp_global_sum_2d_adm(domain, field, field_ad, flags, position, tile_count)

    implicit none
    type(domain2d), intent(in) :: domain
    real(r_grid), intent(in) :: field(:, :)
    real(r_grid), intent(in) :: field_ad(:, :)
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: position
    integer, intent(in), optional :: tile_count

    mpp_global_sum_2d_adm = 0.0

  end function mpp_global_sum_2d_adm


! mpp_update_domains_adm
! ----------------------

subroutine mpp_update_domain2d_2d_adm(array, arrayp, domain, flags, complete, position, &
                                             whalo, ehalo, shalo, nhalo, name, tile_count)

  real, dimension(:,:), intent(inout) :: array, arrayp
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

   call mpp_update_domains_ad(arrayp, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_update_domain2d_2d_adm

subroutine mpp_update_domain2d_3d_adm(array, arrayp, domain, flags, complete, position, &
                                             whalo, ehalo, shalo, nhalo, name, tile_count)

  real, dimension(:,:,:), intent(inout) :: array, arrayp 
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

   call mpp_update_domains_ad(arrayp, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_update_domain2d_3d_adm

subroutine mpp_update_domain2d_4d_adm(array, arrayp, domain, flags, complete, position, &
                                             whalo, ehalo, shalo, nhalo, name, tile_count)

  real, dimension(:,:,:,:), intent(inout) :: array, arrayp
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')
   call mpp_update_domains_ad(arrayp, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_update_domain2d_4d_adm

subroutine mpp_update_domain2d_5d_adm(array, arrayp, domain, flags, complete, position, &
                                             whalo, ehalo, shalo, nhalo, name, tile_count)

  real, dimension(:,:,:,:,:), intent(inout) :: array, arrayp
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

   call mpp_update_domains_ad(arrayp, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_update_domain2d_5d_adm

subroutine mpp_update_domain2d_2dv_adm( u_cmpt, u_cmptp, v_cmpt, v_cmptp, domain, flags, gridtype, complete, &
                                                        whalo, ehalo, shalo, nhalo, name, tile_count )

  real, dimension(:,:), intent(inout) :: u_cmpt, v_cmpt
  real, dimension(:,:), intent(inout) :: u_cmptp, v_cmptp
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

   call mpp_update_domains_ad(u_cmptp, v_cmptp, domain, flags=flags, gridtype=gridtype, &
                                              whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_update_domain2d_2dv_adm

subroutine mpp_update_domain2d_3dv_adm( u_cmpt, u_cmptp, v_cmpt, v_cmptp, domain, flags, gridtype, complete, &
                                                        whalo, ehalo, shalo, nhalo, name, tile_count )

  real, dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt
  real, dimension(:,:,:), intent(inout) :: u_cmptp, v_cmptp
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

   call mpp_update_domains_ad(u_cmptp, v_cmptp, domain, flags=flags, gridtype=gridtype, &
                                              whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_update_domain2d_3dv_adm

subroutine mpp_update_domain2d_4dv_adm( u_cmpt, u_cmptp, v_cmpt, v_cmptp, domain, flags, gridtype, complete, &
                                                        whalo, ehalo, shalo, nhalo, name, tile_count )

  real, dimension(:,:,:,:), intent(inout) :: u_cmpt, v_cmpt
  real, dimension(:,:,:,:), intent(inout) :: u_cmptp, v_cmptp
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

   call mpp_update_domains_ad(u_cmptp, v_cmptp, domain, flags=flags, gridtype=gridtype, &
                                              whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_update_domain2d_4dv_adm

subroutine mpp_update_domain2d_5dv_adm( u_cmpt, u_cmptp, v_cmpt, v_cmptp, domain, flags, gridtype, complete, &
                                                        whalo, ehalo, shalo, nhalo, name, tile_count )

  real, dimension(:,:,:,:,:), intent(inout) :: u_cmpt, v_cmpt
  real, dimension(:,:,:,:,:), intent(inout) :: u_cmptp, v_cmptp
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

   call mpp_update_domains_ad(u_cmptp, v_cmptp, domain, flags=flags, gridtype=gridtype, &
                                              whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_update_domain2d_5dv_adm


! mpp_get_boundary_adm
! --------------------

 subroutine mpp_get_boundary_2d_adm( array, arrayp, domain, &
                                            ebuffer, sbuffer, wbuffer, nbuffer, &
                                            ebuffer_ad, sbuffer_ad, wbuffer_ad, nbuffer_ad, &
                                            flags, position, complete, tile_count )

  real, dimension(:,:), intent(in) :: array, arrayp
  type(domain2d), intent(in) :: domain
  real, intent(inout), optional :: ebuffer(:), sbuffer(:), wbuffer(:), nbuffer(:)
  real, intent(inout), optional :: ebuffer_ad(:), sbuffer_ad(:), wbuffer_ad(:), nbuffer_ad(:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

   call mpp_get_boundary_ad( arrayp, domain,ebuffer=ebuffer_ad,sbuffer=sbuffer_ad,wbuffer=wbuffer_ad,nbuffer=nbuffer_ad,&
                                     flags = flags, position = position, complete = complete, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

 end subroutine mpp_get_boundary_2d_adm

 subroutine mpp_get_boundary_3d_adm( array, arrayp, domain, &
                                            ebuffer, sbuffer, wbuffer, nbuffer, &
                                            ebuffer_ad, sbuffer_ad, wbuffer_ad, nbuffer_ad, &
                                            flags, position, complete, tile_count )

  real, dimension(:,:,:), intent(in) :: array, arrayp
  type(domain2d), intent(in) :: domain
  real, intent(inout), optional :: ebuffer(:,:), sbuffer(:,:), wbuffer(:,:), nbuffer(:,:)
  real, intent(inout), optional :: ebuffer_ad(:,:), sbuffer_ad(:,:), wbuffer_ad(:,:), nbuffer_ad(:,:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

   call mpp_get_boundary_ad( arrayp, domain,ebuffer=ebuffer_ad,sbuffer=sbuffer_ad,wbuffer=wbuffer_ad,nbuffer=nbuffer_ad,&
                                     flags = flags, position = position, complete = complete, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_get_boundary_3d_adm

 subroutine mpp_get_boundary_2dv_adm( u_cmpt, u_cmptp, v_cmpt, v_cmptp, domain, &
                                                      ebufferx, sbufferx, wbufferx, nbufferx, &
                                                      ebuffery, sbuffery, wbuffery, nbuffery, &
                                                      ebufferx_ad, sbufferx_ad, wbufferx_ad, nbufferx_ad, &
                                                      ebuffery_ad, sbuffery_ad, wbuffery_ad, nbuffery_ad, &
                                                      flags, gridtype, complete, tile_count )

  real, dimension(:,:), intent(in) :: u_cmpt, v_cmpt
  real, dimension(:,:), intent(in) :: u_cmptp, v_cmptp
  type(domain2d), intent(in) :: domain
  real, intent(inout), optional :: ebufferx(:), sbufferx(:), wbufferx(:), nbufferx(:)
  real, intent(inout), optional :: ebuffery(:), sbuffery(:), wbuffery(:), nbuffery(:)
  real, intent(inout), optional :: ebufferx_ad(:), sbufferx_ad(:), wbufferx_ad(:), nbufferx_ad(:)
  real, intent(inout), optional :: ebuffery_ad(:), sbuffery_ad(:), wbuffery_ad(:), nbuffery_ad(:)
  integer,     intent(in),    optional :: flags, gridtype, tile_count
  logical,     intent(in),    optional :: complete

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

  call mpp_get_boundary_ad( u_cmptp, v_cmptp, domain, &
                            ebufferx = ebufferx_ad, sbufferx = sbufferx_ad, wbufferx = wbufferx_ad, nbufferx = nbufferx_ad, &
                            ebuffery = ebuffery_ad, sbuffery = sbuffery_ad, wbuffery = wbuffery_ad, nbuffery = nbuffery_ad, &
                            flags = flags, gridtype = gridtype, &
                            complete = complete, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_get_boundary_2dv_adm

 subroutine mpp_get_boundary_3dv_adm( u_cmpt, u_cmptp, v_cmpt, v_cmptp, domain, &
                                                      ebufferx, sbufferx, wbufferx, nbufferx, &
                                                      ebuffery, sbuffery, wbuffery, nbuffery, &
                                                      ebufferx_ad, sbufferx_ad, wbufferx_ad, nbufferx_ad, &
                                                      ebuffery_ad, sbuffery_ad, wbuffery_ad, nbuffery_ad, &
                                                      flags, gridtype, complete, tile_count )

  real, dimension(:,:,:), intent(in) :: u_cmpt, v_cmpt 
  real, dimension(:,:,:), intent(in) :: u_cmptp, v_cmptp
  type(domain2d), intent(in) :: domain
  real, intent(inout), optional :: ebufferx(:,:), sbufferx(:,:), wbufferx(:,:), nbufferx(:,:)
  real, intent(inout), optional :: ebuffery(:,:), sbuffery(:,:), wbuffery(:,:), nbuffery(:,:)
  real, intent(inout), optional :: ebufferx_ad(:,:), sbufferx_ad(:,:), wbufferx_ad(:,:), nbufferx_ad(:,:)
  real, intent(inout), optional :: ebuffery_ad(:,:), sbuffery_ad(:,:), wbuffery_ad(:,:), nbuffery_ad(:,:)
  integer,     intent(in),    optional :: flags, gridtype, tile_count
  logical,     intent(in),    optional :: complete

  if (fv_timing_onoff) call timing_on('  BWD_COMM_TOTAL')

  call mpp_get_boundary_ad( u_cmptp, v_cmptp, domain, &
                            ebufferx = ebufferx_ad, sbufferx = sbufferx_ad, wbufferx = wbufferx_ad, nbufferx = nbufferx_ad, &
                            ebuffery = ebuffery_ad, sbuffery = sbuffery_ad, wbuffery = wbuffery_ad, nbuffery = nbuffery_ad, &
                            flags = flags, gridtype = gridtype, &
                            complete = complete, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  BWD_COMM_TOTAL')

end subroutine mpp_get_boundary_3dv_adm


end module fv_mp_adm_mod
