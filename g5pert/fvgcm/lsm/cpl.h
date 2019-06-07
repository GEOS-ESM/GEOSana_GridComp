* MPI info

      integer     info          ! msg lib error return code
      integer     tid_cpl       ! cpl task identifier
      parameter ( tid_cpl = 0)

* Msg id

      integer msg_id            ! The message type. Special 
			        ! message numbers are used on initalization.  

* Ibuff information received from coupler

      integer nibuff                  
      parameter (nibuff=100)
      integer ibuff(nibuff)

* Coupler timing information 

      integer irtc               ! real time clock
      integer irtc_w             ! rc ticks when waiting for msg
      integer irtc_r             ! rtc tics when msg recved
      integer irtc_s             ! rtc tics when msg sent

      logical csm_timing
      common /pvm_tim/ csm_timing

* Arrays sent/recvd to/from coupler

      integer
     $     nsnd,                 ! number of send variables
     $     nrcv,                 ! number of recv variables
     $     nsizs,                ! size of send array
     $     nsizr                 ! size of recv array
      parameter(
     $     nsnd=12,
     $     nrcv=14,
     $     nsizs=lsmlon*lsmlat*nsnd,
     $     nsizr=lsmlon*lsmlat*nrcv )
      real
     $     arput(lsmlon,lsmlat,nsnd),! send array
     $     arget(lsmlon,lsmlat,nrcv) ! recv array


 
