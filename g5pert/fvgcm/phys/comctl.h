!
! $Id$
! $Author$
!
!
! Model control variables
!
      common/comctc/settrace
      character*4 settrace

      common/comctl/itsst   ,nsrest  ,iradsw  ,iradlw  ,iradae  
      common/comctl/anncyc  ,nlend   ,nlres   ,nlhst   ,lbrnch  
      common/comctl/aeres   ,ozncyc  ,sstcyc  ,dodiavg ,aeregen 
      common/comctl/cpuchek ,incorhst,incorbuf,incorrad,adiabatic
      common/comctl/flxave  ,interp

      integer itsst    ! Sea surf. temp. update freq. (iters)
      integer nsrest   ! Restart flag
      integer iradsw   ! Iteration frequency for shortwave radiation
      integer iradlw   ! Iteration frequency for longwave radiation
      integer iradae   ! Iteration freq. for absorptivity/emissivity

      logical anncyc   ! Do annual cycle (otherwise perpetual)
      logical nlend    ! Flag for end of run
      logical nlres    ! If true, continuation run
      logical nlhst    ! If true, regeneration run
      logical lbrnch   ! If true, branch run
      logical aeres    ! If true, a/e data will be stored on restart file
      logical ozncyc   ! If true, cycle ozone dataset
      logical sstcyc   ! If true, cycle sst dataset
      logical dodiavg  ! true => diurnal averaging
      logical aeregen  ! true => absor/emis part of regeneration data
      logical cpuchek  ! If true, check remaining cpu time each writeup
      logical incorhst ! true => keep history buffer in-core
      logical incorbuf ! true => keep model buffers in-core
      logical incorrad ! true => keep abs/ems buffer in-core
      logical adiabatic! true => no physics
      logical flxave   ! true => send to coupler only on radiation time steps
      logical interp   ! if true interpolate initial conditions


 
