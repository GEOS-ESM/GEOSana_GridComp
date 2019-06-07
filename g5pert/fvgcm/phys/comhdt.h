c
c $Id$
c $Author$
c
C
C Save certain history tape header values
C
      common /comhdt/ nstepht(ptapes)
C
      integer nstepht       ! Save time step to calc nitslf
C
      common /comhdtc/
     $             ldhstvss ,lthstvss ,lshstvss ,
     $             ldhstvos ,lthstvos ,lshstvos 
C
      character*8 ldhstvss  ! Date last header record written on sst data
      character*8 lthstvss  ! Time last header record written on sst data
      character*8 lshstvss  ! Seq. no. of last header record write, sst
      character*8 ldhstvos  ! Date last header record written, ozone data
      character*8 lthstvos  ! Time last header record written, ozone data
      character*8 lshstvos  ! Seq. no. of last header record write, ozone
C
