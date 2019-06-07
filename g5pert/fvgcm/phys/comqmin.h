c
c $Id$
c $Author$
c
C
C Minimum mass mixing ratio for constituents (water vapor is first)
C
      common /comqmin/ qmin(pcnst),qmincg(pcnst)
C
      real qmin      ! Global minimum constituent concentration
      real qmincg    ! Min. constituent concentration counter-gradient term
C
 
