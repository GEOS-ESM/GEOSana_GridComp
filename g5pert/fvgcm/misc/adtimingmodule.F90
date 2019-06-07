subroutine adtiming_on( fn )
#ifdef TIMING
  use timingmodule, only: timing_off
#endif
  implicit none
  character*(*) :: fn
#ifdef TIMING
  call timing_off( 'ad'//fn )
#endif
end subroutine adtiming_on


subroutine adtiming_off( fn )
#ifdef TIMING
  use timingmodule, only: timing_on
#endif
  implicit none
  character*(*) :: fn
#ifdef TIMING
  call timing_on( 'ad'//fn )
#endif
end subroutine adtiming_off
