c
c $Id$
c $Author$
c
C
C Gravity wave drag constants
C
      common/comgw/cappa   ,gravit  ,vmin    ,v0      ,hmin    ,
     $             fc2     ,fudge   ,tnlim   ,bv2min  ,g2or    ,
     $             g2ocp   ,fgor    ,cpair   ,cpvir
C
      real cappa     ! R/cp
      real gravit    ! Acceleration of gravity
      real vmin      ! Threshold value lowest lev wind for gwd to apply
      real v0        ! Approximately zero wind
      real hmin      ! Threshold value of orographic standard deviation
      real fc2       ! Critical froude number squared
      real fudge     ! Tunable parameter e*mu/2 = e*pi/l
      real tnlim     ! Max |v| tendency, 250 m/s/day
      real bv2min    ! Min value of Brundt-Vaisala frequency squared
      real g2or      ! -(gravit**2)/rair
      real g2ocp     !  (gravit**2)/cpair
      real fgor      ! Fudge*gravit/rair
      real cpair     ! Specific heat of dry air
      real cpvir     ! Derived constant for cp of moist air
C
 
