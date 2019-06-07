c
c $Id$
c $Author$
c
C
C Character variables for file maintenance
C
      common /comloc/njusn, rest_pfile, rgpath, home_dir
C
      character njusn*8        ! Logon name
      character rest_pfile*80  ! File name for restart dataset
      character rgpath*80      ! Pathname for regeneration dataset nsrest=2
      character home_dir*80    ! Pathname for regeneration dataset nsrest=2
C
 
