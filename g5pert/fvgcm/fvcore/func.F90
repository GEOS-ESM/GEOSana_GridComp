!  14May2007  Todling   Introduced dyn_prog; global change.
!  28Sep2007  Todling   Upgraded to latest interfaces.
subroutine func( prog )

  use precision
  use prognostics, only : dyn_prog
  use stepon, only: stepon_set, stepon_do, nouter, ninner
  use stepon, only: stepon_tape_rec

  use timingModule

  implicit none

  type(dyn_prog) :: prog

!-----------------------------------------------------------------------
! FastOpt: start taping here
!-----------------------------------------------------------------------

#ifdef CHECKPOINTING

!$taf init outer   = static, nouter
!$taf init onetape = static, 1

#ifdef    EXTERNAL_TRAJECTORY
!$taf init dummytape = user, 'xx'
#endif /* EXTERNAL_TRAJECTORY */

#else  /* CHECKPOINTING */

!$taf init inner   = static, ninner
!$taf init onetape = static, 1

#ifndef  TAPING_IN_DYNPKG
!$taf init dynpkg_1     = static, ninner
!$taf init dynpkg_n2    = static, ninner*dynpkg_n2
!$taf init cd_core_tape = static, ninner*dynpkg_nsplit*dynpkg_n2
!$taf init te_map_tape  = static, ninner
#endif

#ifndef  D_SW_TAPING
!$taf init d_sw_tape1   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*1
!$taf init d_sw_tapej   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*(jnp-1)
#endif

#ifndef  SW_CORE_STORING
!$taf init c_sw_tape1   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*1
!$taf init c_sw_tape2   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*(jnp-1)
#endif

#ifndef  TP_CORE_STORING
!$taf init tpcc1_tape   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*(jnp-1)
!$taf init tpcc2_tape   = static, ninner*dynpkg_nsplit*dynpkg_n2 * nl*(jnp-1)
#endif
#endif /* CHECKPOINTING */

  stepon_tape_rec = 0

! Invoke driving routine for time integration

  call timing_on('stepon')
  call stepon_do ( prog )
  call timing_off('stepon')

end subroutine func

