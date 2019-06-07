      subroutine d_split ( nsplit,  im,   jm,  dt )

#if defined ( SPMD )
      use mod_comm, only : gid
#define CPP_PRT_PREFIX  if(gid == 0)
#else
#define CPP_PRT_PREFIX
#endif

! Purpose:
!    If nsplit=0 then determine a good value for nsplit based on
!    resolution and the large-time-step (pdt). The user may have to
!    set this manually if instability occurs.

! Algorithm developer: S.-J. Lin

      use precision
      implicit none

      integer,  intent(in)    ::  im       ! E-W dimension
      integer,  intent(in)    ::  jm       ! N-S dimension
      real(r8), intent(in)    ::  dt       ! Time-step in seconds
                                            ! Negative dt (backward in time
                                            ! integration) is allowed

      integer, intent(inout) ::  nsplit     ! Number of time splits for
                                            ! the fast Lagrangian dynamics

! Local
      real(r8)  dimx
      real(r8)  dim0                      ! base dimension
      real(r8)  dt0                       ! base time step
      real(r8)  ns0                       ! base nsplit for base dimension

      parameter ( dim0 = 180.  )
      parameter ( dt0  = 1800. )
      parameter ( ns0  = 4.    )

      if ( nsplit == 0 ) then
   
          dimx   = max ( im, 2*(jm-1) )
          nsplit = int ( ns0*abs(dt)*dimx/(dt0*dim0) + 0.75 )
          nsplit = max ( 1, nsplit )     ! for cases in which dt or dimx is too small
          CPP_PRT_PREFIX write(6,*) 'nsplit is set to ', nsplit
      endif

      return
      end
