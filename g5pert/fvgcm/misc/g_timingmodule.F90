subroutine g_timing_on( fn )
#ifdef TIMING
  use timingmodule, only: timing_on
#endif
  implicit none
  character*(*) :: fn
#ifdef TIMING
  call timing_on( 'g_'//fn )
#endif
end subroutine g_timing_on


subroutine g_timing_off( fn )
#ifdef TIMING
  use timingmodule, only: timing_off
#endif
  implicit none
  character*(*) :: fn
#ifdef TIMING
  call timing_off( 'g_'//fn )
#endif
end subroutine g_timing_off
