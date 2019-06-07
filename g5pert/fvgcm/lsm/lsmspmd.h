c ------------------------ code history ---------------------------
c source file:       lsmspmd.h
c purpose:           lsmspmd common block for the 
c                    distributed memory version of the code.
c                    initialized in lsmini.F 
c date last revised: July 1995 - lsm version 1
c author:            John Truesdale
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
c -----------------------------------------------------------------
      common/lsmspmd/beglptspmd(lsmlat),endlptspmd(lsmlat),
     $     begkptspmd(lsmlat),endkptspmd(lsmlat),numks,npes,
     $     tids(0:lsmlat-1),masterproc
C
      integer beglptspmd  ! beg lpt index for this processor
      integer endlptspmd  ! beg lpt index for this processor
      integer begkptspmd  ! beginning kpt array for all lats
      integer endkptspmd  ! ending kpt array for all lats
      integer numks       ! number of k pts on processor
      integer npes        ! number of processors
      integer tids        ! proc id array
      logical masterproc  ! proc 0 logical for printing msgs
 
