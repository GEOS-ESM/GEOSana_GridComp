c
c $Id$
c $Author$
c
C
C Global integrals for moisture and mass conservation
C
      common/comqfl/tmass(plat), tqf(plev), dqf(plev), tmass0, tmassf ,
     $              qmassf ,fixmas ,qmass1, qmass2, pdela(plond,plev)
C
      real tmass    ! Mass integral for each latitude pair
      real tqf      ! Global moisture integral for each level
      real dqf      ! Global integral of neg moisture for each level
      real tmass0   ! Specified dry mass of atmosphere
      real tmassf   ! Global mass integral
      real qmassf   ! Global moisture integral
      real fixmas   ! Proportionality factor for ps in dry mass fixer
      real qmass1   ! Contribution to global moisture integral (mass
C                   !  weighting is based upon the "A" part of the 
C                   !  hybrid grid)
      real qmass2   ! Contribution to global moisture integral (mass
C                   !  weighting is based upon the "B" part of the 
C                   !  hybrid grid)
      real pdela    ! pressure difference between interfaces (pressure
C                   !  defined using the "A" part of hybrid grid only)
C
 
