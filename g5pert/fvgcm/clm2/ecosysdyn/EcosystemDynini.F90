subroutine EcosystemDynini ()

! ------------------------ code history ------------------------------
! purpose           : LPJ and other DGVM related initializations
! date first created: October 2000 - lsm version 2
! by whom           : Sam Levis (adapted from LPJ initialization subroutines)
! date last revised :
! by whom           :
! --------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use clm_varder
  use infnan
  use clm_varpar, ONLY: numpft, npftpar, maxpatch_pft
  use clm_varcon, ONLY: tfrz
  use clm_varmap, ONLY: begpatch,endpatch
  use pft_varcon, ONLY: pftpar , tree   , summergreen, raingreen  , sla     , &
                        lm_sapl, sm_sapl, hm_sapl    , rm_sapl    , latosa  , &
                        allom1 , allom2 , allom3     , reinickerp , wooddens, &
                        noveg
  use shr_const_mod, ONLY: SHR_CONST_PI
  implicit none

! ------------------------ local variables ---------------------------
  integer  :: n
  integer  :: j
  integer  :: k
  integer  :: pft
  real(r8) :: pi
  real(r8) :: x
  real(r8) :: stemdiam
  real(r8) :: height_sapl
  real(r8) :: table(0:numpft,1:npftpar)
! --------------------------------------------------------------------

! PFT PARAMETERS (follows LPJ subroutine pftparameters)

!  1  fraction of roots in upper soil layer
!  2  plants with C4 (1) or C3 (0) photosynthetic pathway
!  3  water scalar value at which leaves shed by drought deciduous PFT
!  4  canopy conductance component (gmin, mm/s) not associated with
!     photosynthesis (Haxeltine & Prentice 1996, Table 4)
!  5  maintenance respiration coefficient
!  6  flammability threshold
!  7  maximum foliar N content (mg/g)
!     (Haxeltine & Prentice 1996a, Fig 4)
!  8  fire resistance index
!  9  leaf turnover period (years)
! 10  leaf longevity (years)
! 11  sapwood turnover period (sapwood converted to heartwood) (years)
! 12  root turnover period (years)
! 13  leaf C:N mass ratio
! 14  sapwood C:N mass ratio
! 15  root C:N mass ratio
! 16  leaf type: broadleaved (1), needleleaved (2) or grass (3)
! 17  phenology type: evergreen (1), summergreen (2), raingreen (3),
!     any type (4)
! 18  leaf to root ratio under non-water stressed conditions
! 19  summergreen phenology ramp, GDD5 requirement to grow full leaf canopy
! 20  tree maximum crown area (m2)
! 21  sapling (or grass on initialisation) LAI
! 22  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
! 23  boreal pft (1), non-boreal pft (0)
! 24  low temperature limit for CO2 uptake
! 25  lower range of temperature optimum for photosynthesis
! 26  upper range of temperature optimum for photosynthesis
! 27  high temperature limit for CO2 unptake

!BIOCLIMATIC LIMITS

! 28 minimum coldest monthly mean temperature
! 29 maximum coldest monthly mean temperature
! 30 minimum growing degree days (at or above 5 deg C)
! 31 upper limit of temperature of the warmest month (twmax)
! 32 lower limit of growth efficiency (g/m2)


! ---------------------------------------------------------------------
!      1      2      3      4      5      6      7      8          PFT
! ---------------------------------------------------------------------

  data ((table(pft,n),n=1,8),pft=0,numpft) /               &
     inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      !  0
    0.70,   0.0,  0.00,   0.3,  1.20,  0.15, 100.0,  0.12, &      !  1
    0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      !  2 was 1.20
    0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      !  3 was 1.20
    0.85,   0.0,  0.00,   0.5,  0.50,  0.15, 100.0,  0.12, &      !  4 was 0.20
    0.70,   0.0,  0.00,   0.5,  1.20,  0.15, 100.0,  0.50, &      !  5
    0.70,   0.0,  0.35,   0.5,  0.50,  0.15, 100.0,  0.50, &      !  6 was 0.20
    0.80,   0.0,  0.00,   0.5,  1.20,  0.15, 120.0,  0.12, &      !  7
    0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      !  8 was 1.20
    0.85,   0.0,  0.00,   0.5,  1.20,  0.15, 100.0,  0.12, &      !  9 (=4or5)
    0.80,   0.0,  0.00,   0.5,  1.20,  0.15, 120.0,  0.12, &      ! 10
    0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      ! 11 was 1.20
    0.90,   0.0,  0.35,   0.5,  0.60,  0.15, 100.0,  1.00, &      ! 12 was 1.20
    0.90,   0.0,  0.35,   0.5,  0.60,  0.15, 100.0,  1.00, &      ! 13 was 1.20
    0.90,   1.0,  0.35,   0.5,  1.20,  0.15, 100.0,  1.00, &      ! 14
    0.90,   1.0,  0.35,   0.5,  1.20,  0.15, 100.0,  1.00, &      ! 15 (.not.
    0.90,   1.0,  0.35,   0.5,  1.20,  0.15, 100.0,  1.00/        ! 16 present)


! ---------------------------------------------------------------------
!      9     10     11     12     13     14     15     16     17   PFT
! ---------------------------------------------------------------------

  data ((table(pft,n),n=9,17),pft=0,numpft) /                     &
     inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 0
     2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, & ! 1
     2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, & ! 2
     1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! 3
     2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! 4
     1.0,  1.00,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! 5
     1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, & ! 6
     1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, & ! 7
     1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! 8
     2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! 9
     1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, & !10
     1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & !11
     1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & !12
     1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & !13
     1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & !14
     1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & !15
     1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0/   !16


! ------------------------------------------------------
!       18      19     20      21    22     23     PFT
! ------------------------------------------------------

  data ((table(pft,n),n=18,23),pft=0,numpft) /  &
       inf,    inf,   inf,    inf,  inf,   inf, &  ! 0
       1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 1
       1.0, 1000.0,  15.0,  1.500,  1.2,   1.0, &  ! 2
       1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  ! 3
       1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 4
       1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 5
       1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 6
       1.0,  200.0,  15.0,  1.500,  1.2,   0.0, &  ! 7
       1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  ! 8
       1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 9
       1.0,  200.0,  15.0,  1.500,  1.2,   0.0, &  !10
       1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  !11
      0.75,  100.0,   0.0,  0.001,  1.2,   1.0, &  !12
      0.75,  100.0,   0.0,  0.001,  1.2,   1.0, &  !13
      0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  !14
      0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  !15
      0.75,  100.0,   0.0,  0.001,  1.2,   0.0/    !16

! -------------------------------------
!      24     25     26      27    PFT
! -------------------------------------
  data ((table(pft,n),n=24,27),pft=0,numpft) / &
      inf,   inf,   inf,    inf, & ! 0
     -4.0,  20.0,  30.0,   42.0, & ! 1
     -4.0,  15.0,  25.0,   38.0, & ! 2
     -4.0,  15.0,  25.0,   38.0, & ! 3
      2.0,  25.0,  30.0,   55.0, & ! 4
     -4.0,  20.0,  30.0,   42.0, & ! 5
      2.0,  25.0,  30.0,   55.0, & ! 6
     -4.0,  20.0,  25.0,   38.0, & ! 7
     -4.0,  15.0,  25.0,   38.0, & ! 8
      2.0,  25.0,  30.0,   55.0, & ! 9
     -4.0,  20.0,  25.0,   38.0, & !10
     -4.0,  15.0,  25.0,   38.0, & !11
     -4.0,  10.0,  30.0,   45.0, & !12
     -4.0,  10.0,  30.0,   45.0, & !13
      6.0,  20.0,  45.0,   55.0, & !14
      6.0,  20.0,  45.0,   55.0, & !15
      6.0,  20.0,  45.0,   55.0/   !16


! --------------------------------------------------------
!      28       29      30       31      32     PFT
! --------------------------------------------------------
  data ((table(pft,n),n=28,npftpar),pft=0,numpft) / &
      inf,     inf,    inf,  1000.0,    inf, & !  0
     -2.0,    22.0,  900.0,  1000.0,    0.0, & !  1
    -32.5,    -2.0,  600.0,    23.0,    0.0, & !  2
   9999.9,    -2.0,  350.0,    23.0,    0.0, & !  3 (was -1000.0)
     15.5,  1000.0,    0.0,  1000.0,    0.0, & !  4
      3.0,    18.8, 1200.0,  1000.0,    0.0, & !  5
     15.5,  1000.0,    0.0,  1000.0,    0.0, & !  6
    -17.0,    15.5, 1200.0,  1000.0,    0.0, & !  7
  -1000.0,    -2.0,  350.0,    23.0,    0.0, & !  8
   9999.9,  1000.0,    0.0,  1000.0,    0.0, & !  9 (was    15.5)
   9999.9,    15.5, 1200.0,  1000.0,    0.0, & ! 10 (was   -17.0)
   9999.9,    -2.0,  350.0,    23.0,    0.0, & ! 11 (was -1000.0)
  -1000.0,   -17.0,    0.0,  1000.0,    0.0, & ! 12 (an LSM type)
    -17.0,    15.5,    0.0,  1000.0,    0.0, & ! 13 (was -1000.0)
     15.5,  1000.0,    0.0,  1000.0,    0.0, & ! 14
   9999.9,  1000.0,    0.0,  1000.0,    0.0, & ! 15 (was    15.5)
   9999.9,  1000.0,    0.0,  1000.0,    0.0/   ! 16 (was    15.5)
!----------------------------------------------------------------------------

  pi = SHR_CONST_PI

  pftpar(noveg,:) = table(noveg,:)
  sla(noveg)      = inf
  lm_sapl(noveg)  = inf
  rm_sapl(noveg)  = inf
  sm_sapl(noveg)  = inf
  hm_sapl(noveg)  = inf
  tree(noveg)        = .false.
  summergreen(noveg) = .false.
  raingreen(noveg)   = .false.

  do pft = 1, numpft

! Transfer parameter values to array pftpar

     do n = 1, npftpar
        pftpar(pft,n) = table(pft,n)
     enddo

! Assign leaf and phenology logicals

     if (pftpar(pft,16) <= 2.0) then     !woody vegetation: trees, shrubs
        tree(pft) = .true.
     else                                !non woody vegetation: grasses
        tree(pft) = .false.
     endif

     if     (pftpar(pft,17) == 1.0) then !evergreen
        summergreen(pft) = .false.
        raingreen(pft)   = .false.
     elseif (pftpar(pft,17) == 2.0) then !summergreen
        summergreen(pft) = .true.
        raingreen(pft)   = .false.
     elseif (pftpar(pft,17) == 3.0) then !raingreen
        summergreen(pft) = .false.
        raingreen(pft)   = .true.
     else                                !any of the above
        summergreen(pft) = .true.
        raingreen(pft)   = .true.
     endif

! Calculate specific leaf area (SLA) for each PFT from leaf longevity
! Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
! Equation based on Reich et al 1997, Fig 1f:

! SLA = 2e-4 * exp(6.15 - 0.46 ln (leaf_longevity * 12))

! SLA in m2/gC, leaf_longevity in years

     sla(pft) = 2.0e-4 * exp(6.15 - 0.46*log(pftpar(pft,10)*12.0))

! Define initial mass structure

     if (tree(pft)) then !woody PFTs

! Calculate leafmass for a sapling individual
!  (1) lai = leafmass * sla / (crown area)
!  (2) (leaf area) = latosa * (sapwood xs area)
!         (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
!  (3) (crown area) = allom1 * (stem diameter) ** reinickerp
!         (Reinickes theory)
! From (1),
!  (4) leafmass = lai * (crown area) / sla
! From (1) & (3),
!  (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
! From (2),
!  (6) leafmass = latosa * (sapwood xs area) / sla
!  (7) (sapwood xs area) = pi * (sapwood diameter)**2 / 4
! From (6) and (7),
!  (8) leafmass = latosa * pi * (sapwood diameter)**2 / 4 / sla
! From (8),
!  (9) (sapwood diameter) = [ 4 * leafmass * sla / pi / latosa ]**0.5
! (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
! Define x,
! (11) x = [ (sapwood diameter)+(heartwood diameter) ] /
!          (sapwood diameter)
! From (10) & (11),
! (12) (stem diameter) = x * (sapwood diameter)
! From (5), (9) & (12),
! (13) leafmass = lai * allom1 * x**reinickerp *
!               (4*leafmass*sla/pi/latosa)**(reinickerp*0.5) / sla
! From (13),
! (14) leafmass = [ lai * allom1 * x**reinickerp *
!      (4*sla/pi/latosa)**(reinickerp*0.5) / sla ]**(2/(2-reinickerp))

        x = pftpar(pft,22)

        lm_sapl(pft) = (pftpar(pft,21) * allom1 * x**reinickerp *          &
          (4.0 * sla(pft) / pi / latosa)**(reinickerp * 0.5) / sla(pft))** &
          (2.0/(2.0-reinickerp)) !eqn 14

! Calculate sapling stem diameter
! From (9) & (12),
! (15) (stem diameter) = x * [ 4 * leafmass * sla / pi / latosa ]**0.5

        stemdiam = x * (4.0*lm_sapl(pft)*sla(pft)/pi/latosa)**0.5 !Eqn 15

! Calculate sapling height
! (16) height = allom2 * (stem diameter)**allom3 (source?)

        height_sapl = allom2 * stemdiam**allom3 !Eqn 16

! Calculate sapling sapwood mass
! (17) (sapwood volume) = height * (sapwood xs area)
! (18) (sapwood xs area) = leafmass * sla / latosa
! From (17) & (18),
! (19) (sapwood volume) = height * leafmass * sla / latosa
! (20) (sapwood mass) = (wood density) * (sapwood volume)
! From (19) & (20),
! (21) (sapwood mass) = (wood density) * height * leafmass * sla / latosa

        sm_sapl(pft)=wooddens*height_sapl*lm_sapl(pft)*sla(pft)/latosa !Eqn 21

! Calculate sapling heartwood mass
! From (11),
! (22) (heartwood mass) = (x-1) * (sapwood mass)

        hm_sapl(pft) = (x-1.0) * sm_sapl(pft) !Eqn 22

     else !grass PFTs

        lm_sapl(pft) = pftpar(pft,21) / sla(pft)

     endif

! Calculate sapling or initial grass rootmass
! (23) lmtorm = (leafmass) / (rootmass)
! where lmtorm=pftpar(pft,18)

     rm_sapl(pft) = lm_sapl(pft) / pftpar(pft,18) !From Eqn 23

  enddo ! pft loop

! ---------------------------------------------------------------
! Some of the following comes from LPJ subroutine initgrid
! ---------------------------------------------------------------

  do k = begpatch,endpatch

! used in Phenology

     clm(k)%t10min = 1.0e+36
     clm(k)%lai_ind = 0.0

! updated in Phenology

     clm(k)%dphen = 0.0
     clm(k)%leafon = 0.0
     clm(k)%leafof = 0.0

! accumulated in FireSeason

     clm(k)%firelength = 0.0       !must reset at end of every year

! used in FireSeason; updated in LitterSOM; updated in annual portion of LPJ

     clm(k)%litterag = 0.0
     clm(k)%litterbg = 0.0

! updated in LitterSOM

     clm(k)%cpool_fast = 0.0
     clm(k)%cpool_slow = 0.0
     clm(k)%k_fast_ave = 0.0
     clm(k)%k_slow_ave = 0.0
     clm(k)%litter_decom_ave = 0.0
     clm(k)%fmicr = 0.0 !initialize b/c use in Biogeochemistry before LitterSOM

! used and updated in annual portion of LPJ

     clm(k)%present  = .false.
     clm(k)%nind     = 0.
     clm(k)%lm_ind   = 0.
     clm(k)%sm_ind   = 0.
     clm(k)%hm_ind   = 0.
     clm(k)%rm_ind   = 0.
     clm(k)%tmomin20 = 0.
     clm(k)%agdd20   = 0.
     clm(k)%t_mo_min = 1.0e+36

! already a variable in LSM but now updated in annual portion of LPJ

     clm(k)%htop    = 0.
     clm(k)%tsai    = 0.
     clm(k)%fpcgrid = 1.0 / maxpatch_pft !not relevant where ist/=istsoil
                                         !consistent with subr. surfrd

! accumulated in Biogeochemistry and used/reset in annual portion of LPJ

     clm(k)%bm_inc = 0.0
     clm(k)%afmicr = 0.0

! accumulated in Stomata and used/reset in annual portion of LPJ

     clm(k)%annpsn = 0.0
     clm(k)%annpsnpot = 0.0

  end do

end subroutine EcosystemDynini
