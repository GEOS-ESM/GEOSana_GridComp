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
module fv_mp_tlm_mod

use fv_arrays_mod,   only: R_GRID
use fv_arrays_nlm_mod,only: fv_timing_onoff
use mpp_domains_mod, only : domain2D
use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
use mpp_domains_mod, only : mpp_group_update_initialized
use mpp_domains_mod, only : mpp_create_group_update,mpp_reset_group_update_field
use mpp_domains_mod, only : group_halo_update_type => mpp_group_update_type

use fv_mp_mod, only : XDir, YDir, ng
use fv_mp_mod, only : is, ie, js, je
use fv_mp_mod, only : isd, ied, jsd, jed
use mpp_domains_mod, only : mpp_update_domains, mpp_get_boundary
use mpp_domains_mod, only : mpp_global_sum

use fv_timing_mod, only: timing_on, timing_off

implicit none
private

#include "mpif.h"

integer :: commglobal, ierror, npes

!fv_mp_mod routines
public complete_group_halo_update
public start_group_halo_update, start_group_halo_update_tlm
public fill_corners_tlm, mp_reduce_sum_tlm

!mpp_domains interface
public mpp_update_domains_tlm, mpp_get_boundary_tlm, mpp_global_sum_tlm


! Regular fv_mp_mod routines
! --------------------------
interface start_group_halo_update
  module procedure start_var_group_update_2d
  module procedure start_var_group_update_3d
  module procedure start_var_group_update_4d
  module procedure start_vector_group_update_2d
  module procedure start_vector_group_update_3d
end interface start_group_halo_update

interface start_group_halo_update_tlm
  module procedure start_var_group_update_2d_tlm
  module procedure start_var_group_update_3d_tlm
  module procedure start_var_group_update_4d_tlm
  module procedure start_vector_group_update_2d_tlm
  module procedure start_vector_group_update_3d_tlm
end interface 

interface fill_corners_tlm
    module procedure fill_corners_2d_r4_tlm
    module procedure fill_corners_2d_r8_tlm
    module procedure fill_corners_xy_2d_r4_tlm
    module procedure fill_corners_xy_2d_r8_tlm
end interface

interface fill_corners_agrid_tlm
    module procedure fill_corners_agrid_r4_tlm
    module procedure fill_corners_agrid_r8_tlm
end interface

interface fill_corners_cgrid_tlm
    module procedure fill_corners_cgrid_r4_tlm
    module procedure fill_corners_cgrid_r8_tlm
end interface

interface fill_corners_dgrid_tlm
    module procedure fill_corners_dgrid_r4_tlm
    module procedure fill_corners_dgrid_r8_tlm
end interface

interface mp_reduce_sum_tlm
    module procedure mp_reduce_sum_r4_tlm
    module procedure mp_reduce_sum_r4_1d_tlm
    module procedure mp_reduce_sum_r8_tlm
    module procedure mp_reduce_sum_r8_1d_tlm
end interface


! These are invented interfaces to mpp_domains
! --------------------------------------------

interface mpp_global_sum_tlm
    module procedure mpp_global_sum_2d_tlm
end interface

interface mpp_update_domains_tlm
    module procedure mpp_update_domain2d_2d_tlm
    module procedure mpp_update_domain2d_3d_tlm
    module procedure mpp_update_domain2d_4d_tlm
    module procedure mpp_update_domain2d_5d_tlm
    module procedure mpp_update_domain2d_2dv_tlm
    module procedure mpp_update_domain2d_3dv_tlm
    module procedure mpp_update_domain2d_4dv_tlm
    module procedure mpp_update_domain2d_5dv_tlm
end interface

interface mpp_get_boundary_tlm
    module procedure mpp_get_boundary_2d_tlm
    module procedure mpp_get_boundary_3d_tlm
    module procedure mpp_get_boundary_2dv_tlm
    module procedure mpp_get_boundary_3dv_tlm
end interface

contains
 

subroutine start_var_group_update_2d(group, array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

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

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_var_group_update_2d


subroutine start_var_group_update_3d(group, array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

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

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_var_group_update_3d

subroutine start_var_group_update_4d(group, array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

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

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_var_group_update_4d



subroutine start_vector_group_update_2d(group, u_cmpt, v_cmpt, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

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

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_vector_group_update_2d

subroutine start_vector_group_update_3d(group, u_cmpt, v_cmpt, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group
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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

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

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_vector_group_update_3d


! start_group_halo_update_tlm
! ---------------------------

 subroutine start_var_group_update_2d_tlm(group, &
#ifndef OLDMPP
                                          group_tl, &
#endif
                                          array, array_tl, domain, flags, position, whalo, ehalo, shalo, nhalo, complete, complete_tl)

  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: group_tl
#endif
  real, dimension(:,:),         intent(inout) :: array, array_tl
  type(domain2D),               intent(inout) :: domain
  integer,      optional,       intent(in)    :: flags
  integer,      optional,       intent(in)    :: position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete, complete_tl
  real                                        :: d_type
  logical                                     :: is_complete, is_complete_tl
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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(array   , domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)
     call mpp_update_domains(array_tl, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

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


  if (mpp_group_update_initialized(group_tl)) then
    call mpp_reset_group_update_field(group_tl,array_tl)
  else
    call mpp_create_group_update(group_tl, array_tl, domain, flags=flags, position=position, &
             whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete_tl = .TRUE.
  if(present(complete_tl)) is_complete_tl = complete_tl
  if(is_complete_tl) then 
     call mpp_start_group_update(group_tl, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_var_group_update_2d_tlm

subroutine start_var_group_update_3d_tlm(group, &
#ifndef OLDMPP
                                         group_tl, &
#endif
                                         array, array_tl, domain, flags, position, whalo, ehalo, shalo, nhalo, complete, complete_tl)

  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: group_tl
#endif
  real, dimension(:,:,:),       intent(inout) :: array, array_tl
  type(domain2D),               intent(inout) :: domain
  integer,           optional,  intent(in)    :: flags
  integer,           optional,  intent(in)    :: position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete, complete_tl
  real                                        :: d_type
  logical                                     :: is_complete, is_complete_tl

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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(array   , domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)
     call mpp_update_domains(array_tl, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

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

  if (mpp_group_update_initialized(group_tl)) then
    call mpp_reset_group_update_field(group_tl,array_tl)
  else
    call mpp_create_group_update(group_tl, array_tl, domain, flags=flags, position=position, &
          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete_tl = .TRUE.
  if(present(complete_tl)) is_complete_tl = complete_tl
  if(is_complete_tl) then
     call mpp_start_group_update(group_tl, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_var_group_update_3d_tlm

subroutine start_var_group_update_4d_tlm(group, &
#ifndef OLDMPP
                                         group_tl, &
#endif
                                         array, array_tl, domain, flags, position, whalo, ehalo, shalo, nhalo, complete, complete_tl)

  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: group_tl
#endif
  real, dimension(:,:,:,:),     intent(inout) :: array, array_tl
  type(domain2D),               intent(inout) :: domain
  integer,           optional,  intent(in)    :: flags
  integer,           optional,  intent(in)    :: position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete, complete_tl
  real                                        :: d_type
  logical                                     :: is_complete, is_complete_tl

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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(array   , domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)
     call mpp_update_domains(array_tl, domain, flags=flags, position=position, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)

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

  if (mpp_group_update_initialized(group_tl)) then
    call mpp_reset_group_update_field(group_tl,array_tl)
  else
    call mpp_create_group_update(group_tl, array_tl, domain, flags=flags, position=position, &
              whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete_tl = .TRUE.
  if(present(complete_tl)) is_complete_tl = complete_tl
  if(is_complete_tl) then
     call mpp_start_group_update(group_tl, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_var_group_update_4d_tlm



subroutine start_vector_group_update_2d_tlm(group,&
#ifndef OLDMPP
                                            group_tl, &
#endif
                                            u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete, complete_tl)

  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: group_tl
#endif
 real,       dimension(:,:),   intent(inout) :: u_cmpt, v_cmpt
  real,       dimension(:,:),   intent(inout) :: u_cmpt_tl, v_cmpt_tl
  type(domain2d),               intent(inout) :: domain
  integer,            optional, intent(in)    :: flags
  integer,            optional, intent(in)    :: gridtype
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete, complete_tl
  real                                        :: d_type
  logical                                     :: is_complete, is_complete_tl

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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(u_cmpt,    v_cmpt,    domain, flags=flags, gridtype=gridtype, &
                                                   whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)
     call mpp_update_domains(u_cmpt_tl, v_cmpt_tl, domain, flags=flags, gridtype=gridtype, &
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


  if (mpp_group_update_initialized(group_tl)) then
    call mpp_reset_group_update_field(group_tl,u_cmpt_tl, v_cmpt_tl)
  else
    call mpp_create_group_update(group_tl, u_cmpt_tl, v_cmpt_tl, domain, &
            flags=flags, gridtype=gridtype, &
            whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete_tl = .TRUE.
  if(present(complete_tl)) is_complete_tl = complete_tl
  if(is_complete_tl) then
     call mpp_start_group_update(group_tl, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_vector_group_update_2d_tlm

subroutine start_vector_group_update_3d_tlm(group, &
#ifndef OLDMPP
                                            group_tl, &
#endif
                                            u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete, complete_tl)

  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: group_tl
#endif
  real,       dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt
  real,       dimension(:,:,:), intent(inout) :: u_cmpt_tl, v_cmpt_tl
  type(domain2d),               intent(inout) :: domain
  integer,            optional, intent(in)    :: flags
  integer,            optional, intent(in)    :: gridtype
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete, complete_tl
  real                                        :: d_type
  logical                                     :: is_complete, is_complete_tl

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

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

#ifdef OLDMPP

     call mpp_update_domains(u_cmpt,    v_cmpt,    domain, flags=flags, gridtype=gridtype, &
                                                   whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, complete=complete)
     call mpp_update_domains(u_cmpt_tl, v_cmpt_tl, domain, flags=flags, gridtype=gridtype, &
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

  if (mpp_group_update_initialized(group_tl)) then
    call mpp_reset_group_update_field(group_tl,u_cmpt_tl, v_cmpt_tl)
  else
    call mpp_create_group_update(group_tl, u_cmpt_tl, v_cmpt_tl, domain, &
            flags=flags, gridtype=gridtype, &
            whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete_tl = .TRUE.
  if(present(complete_tl)) is_complete_tl = complete_tl
  if(is_complete_tl) then
     call mpp_start_group_update(group_tl, domain, d_type)
  endif

#endif

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine start_vector_group_update_3d_tlm



! complete_group_halo_update
! --------------------------

subroutine complete_group_halo_update(group,&
#ifndef OLDMPP
                                      group_tl, &
#endif
                                      domain)
  type(group_halo_update_type), intent(inout) :: group
#ifndef OLDMPP
  type(group_halo_update_type), intent(inout) :: group_tl
#endif
  type(domain2d),               intent(inout) :: domain
  real                                        :: d_type

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!  (in)      domain - Contains domain decomposition information.

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

#ifndef OLDMPP

     call mpp_complete_group_update(group, domain, d_type)
     call mpp_complete_group_update(group_tl, domain, d_type)

#endif

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine complete_group_halo_update


! fill_corners_tlm
! ----------------

  subroutine fill_corners_2d_r4_tlm(q, q_tl, npx, npy, fill, agrid,  bgrid)
    implicit none
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: q
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: q_tl
    integer, intent(in) :: npx, npy
! x-dir or y-dir 
    integer, intent(in) :: fill
    logical, optional, intent(in) :: agrid, bgrid
    integer :: i, j
    intrinsic present
    if (present(bgrid)) then
      if (bgrid) then
        select case  (fill) 
        case (xdir) 
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-i, 1-j) = q_tl(1-j, i+1)
                q(1-i, 1-j) = q(1-j, i+1)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-i, npy+j) = q_tl(1-j, npy-i)
                q(1-i, npy+j) = q(1-j, npy-i)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx+i, 1-j) = q_tl(npx+j, i+1)
                q(npx+i, 1-j) = q(npx+j, i+1)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx+i, npy+j) = q_tl(npx+j, npy-i)
                q(npx+i, npy+j) = q(npx+j, npy-i)
              end if
            end do
          end do
        case (ydir) 
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-j, 1-i) = q_tl(i+1, 1-j)
                q(1-j, 1-i) = q(i+1, 1-j)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-j, npy+i) = q_tl(i+1, npy+j)
                q(1-j, npy+i) = q(i+1, npy+j)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx+j, 1-i) = q_tl(npx-i, 1-j)
                q(npx+j, 1-i) = q(npx-i, 1-j)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx+j, npy+i) = q_tl(npx-i, npy+j)
                q(npx+j, npy+i) = q(npx-i, npy+j)
              end if
            end do
          end do
        case default
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-i, 1-j) = q_tl(1-j, i+1)
                q(1-i, 1-j) = q(1-j, i+1)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-i, npy+j) = q_tl(1-j, npy-i)
                q(1-i, npy+j) = q(1-j, npy-i)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx+i, 1-j) = q_tl(npx+j, i+1)
                q(npx+i, 1-j) = q(npx+j, i+1)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx+i, npy+j) = q_tl(npx+j, npy-i)
                q(npx+i, npy+j) = q(npx+j, npy-i)
              end if
            end do
          end do
        end select
      end if
    else if (present(agrid)) then
      if (agrid) then
        select case  (fill) 
        case (xdir) 
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-i, 1-j) = q_tl(1-j, i)
                q(1-i, 1-j) = q(1-j, i)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-i, npy-1+j) = q_tl(1-j, npy-1-i+1)
                q(1-i, npy-1+j) = q(1-j, npy-1-i+1)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx-1+i, 1-j) = q_tl(npx-1+j, i)
                q(npx-1+i, 1-j) = q(npx-1+j, i)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx-1+i, npy-1+j) = q_tl(npx-1+j, npy-1-i+1)
                q(npx-1+i, npy-1+j) = q(npx-1+j, npy-1-i+1)
              end if
            end do
          end do
        case (ydir) 
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-j, 1-i) = q_tl(i, 1-j)
                q(1-j, 1-i) = q(i, 1-j)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-j, npy-1+i) = q_tl(i, npy-1+j)
                q(1-j, npy-1+i) = q(i, npy-1+j)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx-1+j, 1-i) = q_tl(npx-1-i+1, 1-j)
                q(npx-1+j, 1-i) = q(npx-1-i+1, 1-j)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx-1+j, npy-1+i) = q_tl(npx-1-i+1, npy-1+j)
                q(npx-1+j, npy-1+i) = q(npx-1-i+1, npy-1+j)
              end if
            end do
          end do
        case default
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-j, 1-i) = q_tl(i, 1-j)
                q(1-j, 1-i) = q(i, 1-j)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-j, npy-1+i) = q_tl(i, npy-1+j)
                q(1-j, npy-1+i) = q(i, npy-1+j)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx-1+j, 1-i) = q_tl(npx-1-i+1, 1-j)
                q(npx-1+j, 1-i) = q(npx-1-i+1, 1-j)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx-1+j, npy-1+i) = q_tl(npx-1-i+1, npy-1+j)
                q(npx-1+j, npy-1+i) = q(npx-1-i+1, npy-1+j)
              end if
            end do
          end do
        end select
      end if
    end if
  end subroutine fill_corners_2d_r4_tlm

  subroutine fill_corners_2d_r8_tlm(q, q_tl, npx, npy, fill, agrid,  bgrid)
    implicit none
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: q
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: q_tl
    integer, intent(in) :: npx, npy
! x-dir or y-dir 
    integer, intent(in) :: fill
    logical, optional, intent(in) :: agrid, bgrid
    integer :: i, j
    intrinsic present
    if (present(bgrid)) then
      if (bgrid) then
        select case  (fill) 
        case (xdir) 
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-i, 1-j) = q_tl(1-j, i+1)
                q(1-i, 1-j) = q(1-j, i+1)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-i, npy+j) = q_tl(1-j, npy-i)
                q(1-i, npy+j) = q(1-j, npy-i)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx+i, 1-j) = q_tl(npx+j, i+1)
                q(npx+i, 1-j) = q(npx+j, i+1)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx+i, npy+j) = q_tl(npx+j, npy-i)
                q(npx+i, npy+j) = q(npx+j, npy-i)
              end if
            end do
          end do
        case (ydir) 
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-j, 1-i) = q_tl(i+1, 1-j)
                q(1-j, 1-i) = q(i+1, 1-j)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-j, npy+i) = q_tl(i+1, npy+j)
                q(1-j, npy+i) = q(i+1, npy+j)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx+j, 1-i) = q_tl(npx-i, 1-j)
                q(npx+j, 1-i) = q(npx-i, 1-j)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx+j, npy+i) = q_tl(npx-i, npy+j)
                q(npx+j, npy+i) = q(npx-i, npy+j)
              end if
            end do
          end do
        case default
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-i, 1-j) = q_tl(1-j, i+1)
                q(1-i, 1-j) = q(1-j, i+1)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-i, npy+j) = q_tl(1-j, npy-i)
                q(1-i, npy+j) = q(1-j, npy-i)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx+i, 1-j) = q_tl(npx+j, i+1)
                q(npx+i, 1-j) = q(npx+j, i+1)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx+i, npy+j) = q_tl(npx+j, npy-i)
                q(npx+i, npy+j) = q(npx+j, npy-i)
              end if
            end do
          end do
        end select
      end if
    else if (present(agrid)) then
      if (agrid) then
        select case  (fill) 
        case (xdir) 
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-i, 1-j) = q_tl(1-j, i)
                q(1-i, 1-j) = q(1-j, i)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-i, npy-1+j) = q_tl(1-j, npy-1-i+1)
                q(1-i, npy-1+j) = q(1-j, npy-1-i+1)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx-1+i, 1-j) = q_tl(npx-1+j, i)
                q(npx-1+i, 1-j) = q(npx-1+j, i)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx-1+i, npy-1+j) = q_tl(npx-1+j, npy-1-i+1)
                q(npx-1+i, npy-1+j) = q(npx-1+j, npy-1-i+1)
              end if
            end do
          end do
        case (ydir) 
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-j, 1-i) = q_tl(i, 1-j)
                q(1-j, 1-i) = q(i, 1-j)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-j, npy-1+i) = q_tl(i, npy-1+j)
                q(1-j, npy-1+i) = q(i, npy-1+j)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx-1+j, 1-i) = q_tl(npx-1-i+1, 1-j)
                q(npx-1+j, 1-i) = q(npx-1-i+1, 1-j)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx-1+j, npy-1+i) = q_tl(npx-1-i+1, npy-1+j)
                q(npx-1+j, npy-1+i) = q(npx-1-i+1, npy-1+j)
              end if
            end do
          end do
        case default
          do j=1,ng
            do i=1,ng
!sw corner 
              if (is .eq. 1 .and. js .eq. 1) then
                q_tl(1-j, 1-i) = q_tl(i, 1-j)
                q(1-j, 1-i) = q(i, 1-j)
              end if
!nw corner
              if (is .eq. 1 .and. je .eq. npy - 1) then
                q_tl(1-j, npy-1+i) = q_tl(i, npy-1+j)
                q(1-j, npy-1+i) = q(i, npy-1+j)
              end if
!se corner
              if (ie .eq. npx - 1 .and. js .eq. 1) then
                q_tl(npx-1+j, 1-i) = q_tl(npx-1-i+1, 1-j)
                q(npx-1+j, 1-i) = q(npx-1-i+1, 1-j)
              end if
!ne corner
              if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
                q_tl(npx-1+j, npy-1+i) = q_tl(npx-1-i+1, npy-1+j)
                q(npx-1+j, npy-1+i) = q(npx-1-i+1, npy-1+j)
              end if
            end do
          end do
        end select
      end if
    end if
  end subroutine fill_corners_2d_r8_tlm

  subroutine fill_corners_xy_2d_r4_tlm(x, x_tl, y, y_tl, npx, npy, dgrid&
&   , agrid, cgrid, vector)
    implicit none
!(isd:ied  ,jsd:jed+1)
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: x
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: x_tl
!(isd:ied+1,jsd:jed  )
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: y
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: y_tl
    integer, intent(in) :: npx, npy
    logical, optional, intent(in) :: dgrid, agrid, cgrid, vector
    integer :: i, j
    real(kind=4) :: mysign
    intrinsic present
    mysign = 1.0
    if (present(vector)) then
      if (vector) mysign = -1.0
    end if
    if (present(dgrid)) then
      call fill_corners_dgrid_tlm(x, x_tl, y, y_tl, npx, npy, mysign)
    else if (present(cgrid)) then
      call fill_corners_cgrid_tlm(x, x_tl, y, y_tl, npx, npy, mysign)
    else if (present(agrid)) then
      call fill_corners_agrid_tlm(x, x_tl, y, y_tl, npx, npy, mysign)
    else
      call fill_corners_agrid_tlm(x, x_tl, y, y_tl, npx, npy, mysign)
    end if
  end subroutine fill_corners_xy_2d_r4_tlm

  subroutine fill_corners_xy_2d_r8_tlm(x, x_tl, y, y_tl, npx, npy, dgrid&
&   , agrid, cgrid, vector)
    implicit none
!(isd:ied  ,jsd:jed+1)
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: x
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: x_tl
!(isd:ied+1,jsd:jed  )
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: y
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: y_tl
    integer, intent(in) :: npx, npy
    logical, optional, intent(in) :: dgrid, agrid, cgrid, vector
    integer :: i, j
    real(kind=8) :: mysign
    intrinsic present
    mysign = 1.0
    if (present(vector)) then
      if (vector) mysign = -1.0
    end if
    if (present(dgrid)) then
      call fill_corners_dgrid_tlm(x, x_tl, y, y_tl, npx, npy, mysign)
    else if (present(cgrid)) then
      call fill_corners_cgrid_tlm(x, x_tl, y, y_tl, npx, npy, mysign)
    else if (present(agrid)) then
      call fill_corners_agrid_tlm(x, x_tl, y, y_tl, npx, npy, mysign)
    else
      call fill_corners_agrid_tlm(x, x_tl, y, y_tl, npx, npy, mysign)
    end if
  end subroutine fill_corners_xy_2d_r8_tlm

  subroutine fill_corners_agrid_r4_tlm(x, x_tl, y, y_tl, npx, npy, &
&   mysign)
    implicit none
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: x
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: x_tl
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: y
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: y_tl
    integer, intent(in) :: npx, npy
    real(kind=4), intent(in) :: mysign
    integer :: i, j
    do j=1,ng
      do i=1,ng
!sw corner
        if (is .eq. 1 .and. js .eq. 1) then
          x_tl(1-i, 1-j) = mysign*y_tl(1-j, i)
          x(1-i, 1-j) = mysign*y(1-j, i)
        end if
!nw corner
        if (is .eq. 1 .and. je .eq. npy - 1) then
          x_tl(1-i, npy-1+j) = y_tl(1-j, npy-1-i+1)
          x(1-i, npy-1+j) = y(1-j, npy-1-i+1)
        end if
!se corner
        if (ie .eq. npx - 1 .and. js .eq. 1) then
          x_tl(npx-1+i, 1-j) = y_tl(npx-1+j, i)
          x(npx-1+i, 1-j) = y(npx-1+j, i)
        end if
!ne corner
        if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
          x_tl(npx-1+i, npy-1+j) = mysign*y_tl(npx-1+j, npy-1-i+1)
          x(npx-1+i, npy-1+j) = mysign*y(npx-1+j, npy-1-i+1)
        end if
      end do
    end do
    do j=1,ng
      do i=1,ng
!sw corner
        if (is .eq. 1 .and. js .eq. 1) then
          y_tl(1-j, 1-i) = mysign*x_tl(i, 1-j)
          y(1-j, 1-i) = mysign*x(i, 1-j)
        end if
!nw corner
        if (is .eq. 1 .and. je .eq. npy - 1) then
          y_tl(1-j, npy-1+i) = x_tl(i, npy-1+j)
          y(1-j, npy-1+i) = x(i, npy-1+j)
        end if
!se corner
        if (ie .eq. npx - 1 .and. js .eq. 1) then
          y_tl(npx-1+j, 1-i) = x_tl(npx-1-i+1, 1-j)
          y(npx-1+j, 1-i) = x(npx-1-i+1, 1-j)
        end if
!ne corner
        if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
          y_tl(npx-1+j, npy-1+i) = mysign*x_tl(npx-1-i+1, npy-1+j)
          y(npx-1+j, npy-1+i) = mysign*x(npx-1-i+1, npy-1+j)
        end if
      end do
    end do
  end subroutine fill_corners_agrid_r4_tlm

  subroutine fill_corners_agrid_r8_tlm(x, x_tl, y, y_tl, npx, npy, &
&   mysign)
    implicit none
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: x
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: x_tl
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: y
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: y_tl
    integer, intent(in) :: npx, npy
    real(kind=8), intent(in) :: mysign
    integer :: i, j
    do j=1,ng
      do i=1,ng
!sw corner
        if (is .eq. 1 .and. js .eq. 1) then
          x_tl(1-i, 1-j) = mysign*y_tl(1-j, i)
          x(1-i, 1-j) = mysign*y(1-j, i)
        end if
!nw corner
        if (is .eq. 1 .and. je .eq. npy - 1) then
          x_tl(1-i, npy-1+j) = y_tl(1-j, npy-1-i+1)
          x(1-i, npy-1+j) = y(1-j, npy-1-i+1)
        end if
!se corner
        if (ie .eq. npx - 1 .and. js .eq. 1) then
          x_tl(npx-1+i, 1-j) = y_tl(npx-1+j, i)
          x(npx-1+i, 1-j) = y(npx-1+j, i)
        end if
!ne corner
        if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
          x_tl(npx-1+i, npy-1+j) = mysign*y_tl(npx-1+j, npy-1-i+1)
          x(npx-1+i, npy-1+j) = mysign*y(npx-1+j, npy-1-i+1)
        end if
      end do
    end do
    do j=1,ng
      do i=1,ng
!sw corner
        if (is .eq. 1 .and. js .eq. 1) then
          y_tl(1-j, 1-i) = mysign*x_tl(i, 1-j)
          y(1-j, 1-i) = mysign*x(i, 1-j)
        end if
!nw corner
        if (is .eq. 1 .and. je .eq. npy - 1) then
          y_tl(1-j, npy-1+i) = x_tl(i, npy-1+j)
          y(1-j, npy-1+i) = x(i, npy-1+j)
        end if
!se corner
        if (ie .eq. npx - 1 .and. js .eq. 1) then
          y_tl(npx-1+j, 1-i) = x_tl(npx-1-i+1, 1-j)
          y(npx-1+j, 1-i) = x(npx-1-i+1, 1-j)
        end if
!ne corner
        if (ie .eq. npx - 1 .and. je .eq. npy - 1) then
          y_tl(npx-1+j, npy-1+i) = mysign*x_tl(npx-1-i+1, npy-1+j)
          y(npx-1+j, npy-1+i) = mysign*x(npx-1-i+1, npy-1+j)
        end if
      end do
    end do
  end subroutine fill_corners_agrid_r8_tlm

  subroutine fill_corners_cgrid_r4_tlm(x, x_tl, y, y_tl, npx, npy, &
&   mysign)
    implicit none
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: x
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: x_tl
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: y
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: y_tl
    integer, intent(in) :: npx, npy
    real(kind=4), intent(in) :: mysign
    integer :: i, j
    do j=1,ng
      do i=1,ng
!sw corner 
        if (is .eq. 1 .and. js .eq. 1) then
          x_tl(1-i, 1-j) = y_tl(j, 1-i)
          x(1-i, 1-j) = y(j, 1-i)
        end if
!nw corner
        if (is .eq. 1 .and. je + 1 .eq. npy) then
          x_tl(1-i, npy-1+j) = mysign*y_tl(j, npy+i)
          x(1-i, npy-1+j) = mysign*y(j, npy+i)
        end if
!se corner
        if (ie + 1 .eq. npx .and. js .eq. 1) then
          x_tl(npx+i, 1-j) = mysign*y_tl(npx-j, 1-i)
          x(npx+i, 1-j) = mysign*y(npx-j, 1-i)
        end if
!ne corner
        if (ie + 1 .eq. npx .and. je + 1 .eq. npy) then
          x_tl(npx+i, npy-1+j) = y_tl(npx-j, npy+i)
          x(npx+i, npy-1+j) = y(npx-j, npy+i)
        end if
      end do
    end do
    do j=1,ng
      do i=1,ng
!sw corner 
        if (is .eq. 1 .and. js .eq. 1) then
          y_tl(1-i, 1-j) = x_tl(1-j, i)
          y(1-i, 1-j) = x(1-j, i)
        end if
!nw corner
        if (is .eq. 1 .and. je + 1 .eq. npy) then
          y_tl(1-i, npy+j) = mysign*x_tl(1-j, npy-i)
          y(1-i, npy+j) = mysign*x(1-j, npy-i)
        end if
!se corner
        if (ie + 1 .eq. npx .and. js .eq. 1) then
          y_tl(npx-1+i, 1-j) = mysign*x_tl(npx+j, i)
          y(npx-1+i, 1-j) = mysign*x(npx+j, i)
        end if
!ne corner
        if (ie + 1 .eq. npx .and. je + 1 .eq. npy) then
          y_tl(npx-1+i, npy+j) = x_tl(npx+j, npy-i)
          y(npx-1+i, npy+j) = x(npx+j, npy-i)
        end if
      end do
    end do
  end subroutine fill_corners_cgrid_r4_tlm

  subroutine fill_corners_cgrid_r8_tlm(x, x_tl, y, y_tl, npx, npy, &
&   mysign)
    implicit none
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: x
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: x_tl
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: y
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: y_tl
    integer, intent(in) :: npx, npy
    real(kind=8), intent(in) :: mysign
    integer :: i, j
    do j=1,ng
      do i=1,ng
!sw corner 
        if (is .eq. 1 .and. js .eq. 1) then
          x_tl(1-i, 1-j) = y_tl(j, 1-i)
          x(1-i, 1-j) = y(j, 1-i)
        end if
!nw corner
        if (is .eq. 1 .and. je + 1 .eq. npy) then
          x_tl(1-i, npy-1+j) = mysign*y_tl(j, npy+i)
          x(1-i, npy-1+j) = mysign*y(j, npy+i)
        end if
!se corner
        if (ie + 1 .eq. npx .and. js .eq. 1) then
          x_tl(npx+i, 1-j) = mysign*y_tl(npx-j, 1-i)
          x(npx+i, 1-j) = mysign*y(npx-j, 1-i)
        end if
!ne corner
        if (ie + 1 .eq. npx .and. je + 1 .eq. npy) then
          x_tl(npx+i, npy-1+j) = y_tl(npx-j, npy+i)
          x(npx+i, npy-1+j) = y(npx-j, npy+i)
        end if
      end do
    end do
    do j=1,ng
      do i=1,ng
!sw corner 
        if (is .eq. 1 .and. js .eq. 1) then
          y_tl(1-i, 1-j) = x_tl(1-j, i)
          y(1-i, 1-j) = x(1-j, i)
        end if
!nw corner
        if (is .eq. 1 .and. je + 1 .eq. npy) then
          y_tl(1-i, npy+j) = mysign*x_tl(1-j, npy-i)
          y(1-i, npy+j) = mysign*x(1-j, npy-i)
        end if
!se corner
        if (ie + 1 .eq. npx .and. js .eq. 1) then
          y_tl(npx-1+i, 1-j) = mysign*x_tl(npx+j, i)
          y(npx-1+i, 1-j) = mysign*x(npx+j, i)
        end if
!ne corner
        if (ie + 1 .eq. npx .and. je + 1 .eq. npy) then
          y_tl(npx-1+i, npy+j) = x_tl(npx+j, npy-i)
          y(npx-1+i, npy+j) = x(npx+j, npy-i)
        end if
      end do
    end do
  end subroutine fill_corners_cgrid_r8_tlm

  subroutine fill_corners_dgrid_r4_tlm(x, x_tl, y, y_tl, npx, npy, &
&   mysign)
    implicit none
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: x
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: x_tl
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: y
    real(kind=4), dimension(isd:, jsd:), intent(inout) :: y_tl
    integer, intent(in) :: npx, npy
    real(kind=4), intent(in) :: mysign
    integer :: i, j
    do j=1,ng
      do i=1,ng
!   if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) =        y(j+1  ,1-i    )  !sw corner 
!   if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) = mysign*y(j+1  ,npy-1+i)  !nw corner
!   if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) = mysign*y(npx-j,1-i    )  !se corner
!   if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) =        y(npx-j,npy-1+i)  !ne corner
!sw corner 
        if (is .eq. 1 .and. js .eq. 1) then
          x_tl(1-i, 1-j) = mysign*y_tl(1-j, i)
          x(1-i, 1-j) = mysign*y(1-j, i)
        end if
!nw corner
        if (is .eq. 1 .and. je + 1 .eq. npy) then
          x_tl(1-i, npy+j) = y_tl(1-j, npy-i)
          x(1-i, npy+j) = y(1-j, npy-i)
        end if
!se corner
        if (ie + 1 .eq. npx .and. js .eq. 1) then
          x_tl(npx-1+i, 1-j) = y_tl(npx+j, i)
          x(npx-1+i, 1-j) = y(npx+j, i)
        end if
!ne corner
        if (ie + 1 .eq. npx .and. je + 1 .eq. npy) then
          x_tl(npx-1+i, npy+j) = mysign*y_tl(npx+j, npy-i)
          x(npx-1+i, npy+j) = mysign*y(npx+j, npy-i)
        end if
      end do
    end do
    do j=1,ng
      do i=1,ng
!  if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) =        x(1-j    ,i+1  )  !sw corner 
!  if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) = mysign*x(1-j    ,npy-i)  !nw corner
!  if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) = mysign*x(npx-1+j,i+1  )  !se corner
!  if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) =        x(npx-1+j,npy-i)  !ne corner
!sw corner 
        if (is .eq. 1 .and. js .eq. 1) then
          y_tl(1-i, 1-j) = mysign*x_tl(j, 1-i)
          y(1-i, 1-j) = mysign*x(j, 1-i)
        end if
!nw corner
        if (is .eq. 1 .and. je + 1 .eq. npy) then
          y_tl(1-i, npy-1+j) = x_tl(j, npy+i)
          y(1-i, npy-1+j) = x(j, npy+i)
        end if
!se corner
        if (ie + 1 .eq. npx .and. js .eq. 1) then
          y_tl(npx+i, 1-j) = x_tl(npx-j, 1-i)
          y(npx+i, 1-j) = x(npx-j, 1-i)
        end if
!ne corner
        if (ie + 1 .eq. npx .and. je + 1 .eq. npy) then
          y_tl(npx+i, npy-1+j) = mysign*x_tl(npx-j, npy+i)
          y(npx+i, npy-1+j) = mysign*x(npx-j, npy+i)
        end if
      end do
    end do
  end subroutine fill_corners_dgrid_r4_tlm

  subroutine fill_corners_dgrid_r8_tlm(x, x_tl, y, y_tl, npx, npy, &
&   mysign)
    implicit none
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: x
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: x_tl
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: y
    real(kind=8), dimension(isd:, jsd:), intent(inout) :: y_tl
    integer, intent(in) :: npx, npy
    real(kind=8), intent(in) :: mysign
    integer :: i, j
    do j=1,ng
      do i=1,ng
!   if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) =        y(j+1  ,1-i    )  !sw corner 
!   if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) = mysign*y(j+1  ,npy-1+i)  !nw corner
!   if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) = mysign*y(npx-j,1-i    )  !se corner
!   if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) =        y(npx-j,npy-1+i)  !ne corner
!sw corner 
        if (is .eq. 1 .and. js .eq. 1) then
          x_tl(1-i, 1-j) = mysign*y_tl(1-j, i)
          x(1-i, 1-j) = mysign*y(1-j, i)
        end if
!nw corner
        if (is .eq. 1 .and. je + 1 .eq. npy) then
          x_tl(1-i, npy+j) = y_tl(1-j, npy-i)
          x(1-i, npy+j) = y(1-j, npy-i)
        end if
!se corner
        if (ie + 1 .eq. npx .and. js .eq. 1) then
          x_tl(npx-1+i, 1-j) = y_tl(npx+j, i)
          x(npx-1+i, 1-j) = y(npx+j, i)
        end if
!ne corner
        if (ie + 1 .eq. npx .and. je + 1 .eq. npy) then
          x_tl(npx-1+i, npy+j) = mysign*y_tl(npx+j, npy-i)
          x(npx-1+i, npy+j) = mysign*y(npx+j, npy-i)
        end if
      end do
    end do
    do j=1,ng
      do i=1,ng
!  if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) =        x(1-j    ,i+1  )  !sw corner 
!  if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) = mysign*x(1-j    ,npy-i)  !nw corner
!  if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) = mysign*x(npx-1+j,i+1  )  !se corner
!  if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) =        x(npx-1+j,npy-i)  !ne corner
!sw corner 
        if (is .eq. 1 .and. js .eq. 1) then
          y_tl(1-i, 1-j) = mysign*x_tl(j, 1-i)
          y(1-i, 1-j) = mysign*x(j, 1-i)
        end if
!nw corner
        if (is .eq. 1 .and. je + 1 .eq. npy) then
          y_tl(1-i, npy-1+j) = x_tl(j, npy+i)
          y(1-i, npy-1+j) = x(j, npy+i)
        end if
!se corner
        if (ie + 1 .eq. npx .and. js .eq. 1) then
          y_tl(npx+i, 1-j) = x_tl(npx-j, 1-i)
          y(npx+i, 1-j) = x(npx-j, 1-i)
        end if
!ne corner
        if (ie + 1 .eq. npx .and. je + 1 .eq. npy) then
          y_tl(npx+i, npy-1+j) = mysign*x_tl(npx-j, npy+i)
          y(npx+i, npy-1+j) = mysign*x(npx-j, npy+i)
        end if
      end do
    end do
  end subroutine fill_corners_dgrid_r8_tlm


! mp_reduce_sum_tlm
! -----------------

  subroutine mp_reduce_sum_r4_tlm(mysum,mysum_tl)
     real(kind=4), intent(INOUT)  :: mysum, mysum_tl

     real(kind=4) :: gsum, gsum_tl

     call MPI_ALLREDUCE( mysum, gsum, 1, MPI_REAL, MPI_SUM, &
                         commglobal, ierror )
     call MPI_ALLREDUCE( mysum_tl, gsum_tl, 1, MPI_REAL, MPI_SUM, &
                         commglobal, ierror )

     mysum = gsum
     mysum_tl = 0.0!gsum_tl

  end subroutine mp_reduce_sum_r4_tlm

  subroutine mp_reduce_sum_r8_tlm(mysum,mysum_tl)
     real(kind=8), intent(INOUT)  :: mysum, mysum_tl

     real(kind=8) :: gsum, gsum_tl

     call MPI_ALLREDUCE( mysum, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         commglobal, ierror )
     call MPI_ALLREDUCE( mysum_tl, gsum_tl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         commglobal, ierror )

     mysum = gsum
     mysum_tl = 0.0!gsum_tl

  end subroutine mp_reduce_sum_r8_tlm

  subroutine mp_reduce_sum_r4_1d_tlm(mysum, mysum_tl, sum1d, sum1d_tl, npts)
     integer, intent(in)  :: npts
     real(kind=4), intent(in)     :: sum1d(npts), sum1d_tl(npts)
     real(kind=4), intent(INOUT)  :: mysum, mysum_tl

     real(kind=4) :: gsum, gsum_tl
     integer :: i

     mysum = 0.0
     mysum_tl = 0.0
     do i=1,npts
        mysum = mysum + sum1d(i)
        mysum_tl = mysum_tl + sum1d_tl(i)
     enddo 

     call MPI_ALLREDUCE( mysum, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         commglobal, ierror )
     call MPI_ALLREDUCE( mysum_tl, gsum_tl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         commglobal, ierror )

     mysum = gsum
     mysum_tl = 0.0!gsum_tl

  end subroutine mp_reduce_sum_r4_1d_tlm

  subroutine mp_reduce_sum_r8_1d_tlm(mysum, mysum_tl, sum1d, sum1d_tl, npts)
     integer, intent(in)  :: npts
     real(kind=8), intent(in)     :: sum1d(npts), sum1d_tl(npts)
     real(kind=8), intent(INOUT)  :: mysum, mysum_tl

     real(kind=8) :: gsum, gsum_tl
     integer :: i

     mysum = 0.0
     mysum_tl = 0.0
     do i=1,npts
        mysum = mysum + sum1d(i)
        mysum_tl = mysum_tl + sum1d_tl(i)
     enddo 

     call MPI_ALLREDUCE( mysum, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         commglobal, ierror )
     call MPI_ALLREDUCE( mysum_tl, gsum_tl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         commglobal, ierror )


     mysum = gsum
     mysum_tl = 0.0!gsum_tl

  end subroutine mp_reduce_sum_r8_1d_tlm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MPP_DOMAINS INTERFACES THAT ARE NOT NORMALLY IN FV_MP_MOD !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! mpp_global_sum_tlm
! ------------------

  real(kind=r_grid) function mpp_global_sum_2d_tlm(domain, field, field_tl, flags, position, tile_count, mpp_global_sum_2d)

    implicit none
    type(domain2d), intent(in) :: domain
    real(r_grid), intent(in) :: field(:, :)
    real(r_grid), intent(in) :: field_tl(:, :)
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: position
    integer, intent(in), optional :: tile_count
    real(kind=r_grid) :: mpp_global_sum_2d

    mpp_global_sum_2d     = mpp_global_sum(domain,field,flags=flags,position=position,tile_count=tile_count)
    mpp_global_sum_2d_tlm = 0.0!mpp_global_sum(domain,field_tl,flags=flags,position=position,tile_count=tile_count)

  end function mpp_global_sum_2d_tlm


! mpp_update_domains_tlm
! ----------------------

subroutine mpp_update_domain2d_2d_tlm(array, array_tl, domain, flags, complete, position, &
                                                       whalo, ehalo, shalo, nhalo, name, tile_count)

  real, dimension(:,:), intent(inout) :: array, array_tl
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  call mpp_update_domains( array, domain, &
                           flags = flags, complete = complete, position = position, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )
  call mpp_update_domains( array_tl, domain, &
                           flags = flags, complete = complete, position = position, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_update_domain2d_2d_tlm

subroutine mpp_update_domain2d_3d_tlm(array, array_tl, domain, flags, complete, position, &
                                                       whalo, ehalo, shalo, nhalo, name, tile_count)

  real, dimension(:,:,:), intent(inout) :: array, array_tl
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  call mpp_update_domains( array, domain, &
                           flags = flags, complete = complete, position = position, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )
  call mpp_update_domains( array_tl, domain, &
                           flags = flags, complete = complete, position = position, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_update_domain2d_3d_tlm

subroutine mpp_update_domain2d_4d_tlm(array, array_tl, domain, flags, complete, position, &
                                                       whalo, ehalo, shalo, nhalo, name, tile_count)

  real, dimension(:,:,:,:), intent(inout) :: array, array_tl
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  call mpp_update_domains( array, domain, &
                           flags = flags, complete = complete, position = position, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )
  call mpp_update_domains( array_tl, domain, &
                           flags = flags, complete = complete, position = position, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_update_domain2d_4d_tlm

subroutine mpp_update_domain2d_5d_tlm(array, array_tl, domain, flags, complete, position, &
                                                       whalo, ehalo, shalo, nhalo, name, tile_count)

  real, dimension(:,:,:,:,:), intent(inout) :: array, array_tl
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  call mpp_update_domains( array, domain, &
                           flags = flags, complete = complete, position = position, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )
  call mpp_update_domains( array_tl, domain, &
                           flags = flags, complete = complete, position = position, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_update_domain2d_5d_tlm

subroutine mpp_update_domain2d_2dv_tlm( u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl, domain, flags, gridtype, complete, &
                                        whalo, ehalo, shalo, nhalo, name, tile_count )

  real, dimension(:,:), intent(inout) :: u_cmpt, v_cmpt, u_cmpt_tl, v_cmpt_tl
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  call mpp_update_domains( u_cmpt,v_cmpt,domain,flags = flags, gridtype = gridtype, complete = complete, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )
  call mpp_update_domains( u_cmpt_tl,v_cmpt_tl,domain,flags = flags, gridtype = gridtype, complete = complete, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_update_domain2d_2dv_tlm

subroutine mpp_update_domain2d_3dv_tlm( u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl, domain, flags, gridtype, complete, &
                                        whalo, ehalo, shalo, nhalo, name, tile_count )

  real, dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt, u_cmpt_tl, v_cmpt_tl
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  call mpp_update_domains( u_cmpt,v_cmpt,domain,flags = flags, gridtype = gridtype, complete = complete, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )
  call mpp_update_domains( u_cmpt_tl,v_cmpt_tl,domain,flags = flags, gridtype = gridtype, complete = complete, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_update_domain2d_3dv_tlm

subroutine mpp_update_domain2d_4dv_tlm( u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl, domain, flags, gridtype, complete, &
                                        whalo, ehalo, shalo, nhalo, name, tile_count )

  real, dimension(:,:,:,:), intent(inout) :: u_cmpt, v_cmpt, u_cmpt_tl, v_cmpt_tl
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  call mpp_update_domains( u_cmpt,v_cmpt,domain,flags = flags, gridtype = gridtype, complete = complete, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )
  call mpp_update_domains( u_cmpt_tl,v_cmpt_tl,domain,flags = flags, gridtype = gridtype, complete = complete, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_update_domain2d_4dv_tlm

subroutine mpp_update_domain2d_5dv_tlm( u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl, domain, flags, gridtype, complete, &
                                        whalo, ehalo, shalo, nhalo, name, tile_count )

  real, dimension(:,:,:,:,:), intent(inout) :: u_cmpt, v_cmpt, u_cmpt_tl, v_cmpt_tl
  type(domain2d), intent(inout) :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  call mpp_update_domains( u_cmpt,v_cmpt,domain,flags = flags, gridtype = gridtype, complete = complete, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )
  call mpp_update_domains( u_cmpt_tl,v_cmpt_tl,domain,flags = flags, gridtype = gridtype, complete = complete, &
                           whalo = whalo, ehalo = ehalo, shalo = shalo, nhalo = nhalo, &
                           name = name, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_update_domain2d_5dv_tlm


! mpp_get_boundary_tlm
! --------------------

 subroutine mpp_get_boundary_2d_tlm( array, array_tl, domain, &
                                     ebuffer, ebuffer_tl, &
                                     sbuffer, sbuffer_tl, &
                                     wbuffer, wbuffer_tl, &
                                     nbuffer, nbuffer_tl, &
                                     flags, position, complete, tile_count )

  real, dimension(:,:), intent(in) :: array, array_tl
  type(domain2d), intent(in) :: domain
  real, intent(inout), optional :: ebuffer(:), sbuffer(:), wbuffer(:), nbuffer(:)
  real, intent(inout), optional :: ebuffer_tl(:), sbuffer_tl(:), wbuffer_tl(:), nbuffer_tl(:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  if (present(wbuffer)) wbuffer = 0.0
  if (present(sbuffer)) sbuffer = 0.0
  if (present(ebuffer)) ebuffer = 0.0
  if (present(nbuffer)) nbuffer = 0.0
  if (present(wbuffer_tl)) wbuffer_tl = 0.0
  if (present(sbuffer_tl)) sbuffer_tl = 0.0
  if (present(ebuffer_tl)) ebuffer_tl = 0.0
  if (present(nbuffer_tl)) nbuffer_tl = 0.0

  call mpp_get_boundary( array,domain,ebuffer=ebuffer,sbuffer=sbuffer,wbuffer=wbuffer,nbuffer=nbuffer,&
                         flags = flags, position = position, complete = complete, tile_count = tile_count )
  call mpp_get_boundary( array_tl,domain,ebuffer=ebuffer_tl,sbuffer=sbuffer_tl,wbuffer=wbuffer_tl,nbuffer=nbuffer_tl,&
                         flags = flags, position = position, complete = complete, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_get_boundary_2d_tlm

 subroutine mpp_get_boundary_3d_tlm( array, array_tl, domain, &
                                     ebuffer, ebuffer_tl, &
                                     sbuffer, sbuffer_tl, &
                                     wbuffer, wbuffer_tl, &
                                     nbuffer, nbuffer_tl, &
                                     flags, position, complete, tile_count )

  real, dimension(:,:,:), intent(in) :: array, array_tl
  type(domain2d), intent(in) :: domain
  real, intent(inout), optional :: ebuffer(:,:), sbuffer(:,:), wbuffer(:,:), nbuffer(:,:)
  real, intent(inout), optional :: ebuffer_tl(:,:), sbuffer_tl(:,:), wbuffer_tl(:,:), nbuffer_tl(:,:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete

  if (present(wbuffer)) wbuffer = 0.0
  if (present(sbuffer)) sbuffer = 0.0
  if (present(ebuffer)) ebuffer = 0.0
  if (present(nbuffer)) nbuffer = 0.0
  if (present(wbuffer_tl)) wbuffer_tl = 0.0
  if (present(sbuffer_tl)) sbuffer_tl = 0.0
  if (present(ebuffer_tl)) ebuffer_tl = 0.0
  if (present(nbuffer_tl)) nbuffer_tl = 0.0

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  call mpp_get_boundary( array,domain,ebuffer=ebuffer,sbuffer=sbuffer,wbuffer=wbuffer,nbuffer=nbuffer,&
                         flags = flags, position = position, complete = complete, tile_count = tile_count )
  call mpp_get_boundary( array_tl,domain,ebuffer=ebuffer_tl,sbuffer=sbuffer_tl,wbuffer=wbuffer_tl,nbuffer=nbuffer_tl,&
                         flags = flags, position = position, complete = complete, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_get_boundary_3d_tlm

 subroutine mpp_get_boundary_2dv_tlm( u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl, domain, &
                                      ebufferx, ebufferx_tl, &
                                      sbufferx, sbufferx_tl, &
                                      wbufferx, wbufferx_tl, &
                                      nbufferx, nbufferx_tl, &
                                      ebuffery, ebuffery_tl, &
                                      sbuffery, sbuffery_tl, &
                                      wbuffery, wbuffery_tl, &
                                      nbuffery, nbuffery_tl, &
                                      flags, gridtype, complete, tile_count )

  real, dimension(:,:), intent(in) :: u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl
  type(domain2d), intent(in) :: domain
  real, intent(inout), optional :: ebufferx(:), sbufferx(:), wbufferx(:), nbufferx(:)
  real, intent(inout), optional :: ebufferx_tl(:), sbufferx_tl(:), wbufferx_tl(:), nbufferx_tl(:)
  real, intent(inout), optional :: ebuffery(:), sbuffery(:), wbuffery(:), nbuffery(:)
  real, intent(inout), optional :: ebuffery_tl(:), sbuffery_tl(:), wbuffery_tl(:), nbuffery_tl(:)
  integer,     intent(in),    optional :: flags, gridtype, tile_count
  logical,     intent(in),    optional :: complete

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  if (present(wbufferx)) wbufferx = 0.0
  if (present(sbufferx)) sbufferx = 0.0
  if (present(ebufferx)) ebufferx = 0.0
  if (present(nbufferx)) nbufferx = 0.0
  if (present(wbufferx_tl)) wbufferx_tl = 0.0
  if (present(sbufferx_tl)) sbufferx_tl = 0.0
  if (present(ebufferx_tl)) ebufferx_tl = 0.0
  if (present(nbufferx_tl)) nbufferx_tl = 0.0
  if (present(wbuffery)) wbuffery = 0.0
  if (present(sbuffery)) sbuffery = 0.0
  if (present(ebuffery)) ebuffery = 0.0
  if (present(nbuffery)) nbuffery = 0.0
  if (present(wbuffery_tl)) wbuffery_tl = 0.0
  if (present(sbuffery_tl)) sbuffery_tl = 0.0
  if (present(ebuffery_tl)) ebuffery_tl = 0.0
  if (present(nbuffery_tl)) nbuffery_tl = 0.0

  call mpp_get_boundary( u_cmpt, v_cmpt, domain, &
                         ebufferx = ebufferx, &
                         sbufferx = sbufferx, &
                         wbufferx = wbufferx, &
                         nbufferx = nbufferx, &
                         ebuffery = ebuffery, &
                         sbuffery = sbuffery, &
                         wbuffery = wbuffery, &
                         nbuffery = nbuffery, &
                         flags = flags, gridtype = gridtype, &
                         complete = complete, tile_count = tile_count )

  call mpp_get_boundary( u_cmpt_tl, v_cmpt_tl, domain, &
                         ebufferx = ebufferx_tl, &
                         sbufferx = sbufferx_tl, &
                         wbufferx = wbufferx_tl, &
                         nbufferx = nbufferx_tl, &
                         ebuffery = ebuffery_tl, &
                         sbuffery = sbuffery_tl, &
                         wbuffery = wbuffery_tl, &
                         nbuffery = nbuffery_tl, &
                         flags = flags, gridtype = gridtype, &
                         complete = complete, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_get_boundary_2dv_tlm

 subroutine mpp_get_boundary_3dv_tlm( u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl, domain, &
                                      ebufferx, ebufferx_tl, &
                                      sbufferx, sbufferx_tl, &
                                      wbufferx, wbufferx_tl, &
                                      nbufferx, nbufferx_tl, &
                                      ebuffery, ebuffery_tl, &
                                      sbuffery, sbuffery_tl, &
                                      wbuffery, wbuffery_tl, &
                                      nbuffery, nbuffery_tl, &
                                      flags, gridtype, complete, tile_count )

  real, dimension(:,:,:), intent(in) :: u_cmpt, u_cmpt_tl, v_cmpt, v_cmpt_tl
  type(domain2d), intent(in) :: domain
  real, intent(inout), optional :: ebufferx(:,:), sbufferx(:,:), wbufferx(:,:), nbufferx(:,:)
  real, intent(inout), optional :: ebufferx_tl(:,:), sbufferx_tl(:,:), wbufferx_tl(:,:), nbufferx_tl(:,:)
  real, intent(inout), optional :: ebuffery(:,:), sbuffery(:,:), wbuffery(:,:), nbuffery(:,:)
  real, intent(inout), optional :: ebuffery_tl(:,:), sbuffery_tl(:,:), wbuffery_tl(:,:), nbuffery_tl(:,:)
  integer,     intent(in),    optional :: flags, gridtype, tile_count
  logical,     intent(in),    optional :: complete

  if (fv_timing_onoff) call timing_on('  TLM_COMM_TOTAL')

  if (present(wbufferx)) wbufferx = 0.0
  if (present(sbufferx)) sbufferx = 0.0
  if (present(ebufferx)) ebufferx = 0.0
  if (present(nbufferx)) nbufferx = 0.0
  if (present(wbufferx_tl)) wbufferx_tl = 0.0
  if (present(sbufferx_tl)) sbufferx_tl = 0.0
  if (present(ebufferx_tl)) ebufferx_tl = 0.0
  if (present(nbufferx_tl)) nbufferx_tl = 0.0
  if (present(wbuffery)) wbuffery = 0.0
  if (present(sbuffery)) sbuffery = 0.0
  if (present(ebuffery)) ebuffery = 0.0
  if (present(nbuffery)) nbuffery = 0.0
  if (present(wbuffery_tl)) wbuffery_tl = 0.0
  if (present(sbuffery_tl)) sbuffery_tl = 0.0
  if (present(ebuffery_tl)) ebuffery_tl = 0.0
  if (present(nbuffery_tl)) nbuffery_tl = 0.0

  call mpp_get_boundary( u_cmpt, v_cmpt, domain, &
                         ebufferx = ebufferx, &
                         sbufferx = sbufferx, &
                         wbufferx = wbufferx, &
                         nbufferx = nbufferx, &
                         ebuffery = ebuffery, &
                         sbuffery = sbuffery, &
                         wbuffery = wbuffery, &
                         nbuffery = nbuffery, &
                         flags = flags, gridtype = gridtype, &
                         complete = complete, tile_count = tile_count )

  call mpp_get_boundary( u_cmpt_tl, v_cmpt_tl, domain, &
                         ebufferx = ebufferx_tl, &
                         sbufferx = sbufferx_tl, &
                         wbufferx = wbufferx_tl, &
                         nbufferx = nbufferx_tl, &
                         ebuffery = ebuffery_tl, &
                         sbuffery = sbuffery_tl, &
                         wbuffery = wbuffery_tl, &
                         nbuffery = nbuffery_tl, &
                         flags = flags, gridtype = gridtype, &
                         complete = complete, tile_count = tile_count )

  if (fv_timing_onoff) call timing_off('  TLM_COMM_TOTAL')

end subroutine mpp_get_boundary_3dv_tlm

end module fv_mp_tlm_mod
