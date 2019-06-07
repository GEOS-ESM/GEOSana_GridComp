#include <misc.h>
#include <preproc.h>

subroutine Establishment ()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: Called once per year
! 
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. establishment)
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use clmtype
  use clm_varder
  use clm_varmap, only : patchvec, begpatch, endpatch, begland, endland
  use clm_varpar, only : numpft
  use clm_varcon, only : tfrz, istsoil
  use pft_varcon, only : pftpar, latosa, allom1, allom2, allom3,   &
                         wooddens, reinickerp, tree, sla, lm_sapl, &
                         sm_sapl, hm_sapl, rm_sapl, noveg
  use shr_const_mod, only : SHR_CONST_CDAY,SHR_CONST_PI
  implicit none

! ------------------------ local variables -----------------------------
  integer  :: i,j,k,l,m
  integer  :: ngrass(begland:endland)
  integer  :: npft_estab(begland:endland)
  logical  :: present(begland:endland,0:numpft)
  logical  :: survive(begland:endland,0:numpft)
  logical  :: estab(begland:endland,0:numpft)
  real(r8) :: sumwt(begland:endland)
  real(r8) :: fpc_tree_total(begland:endland)
  real(r8) :: fpc_grass_total(begland:endland)
  real(r8) :: fpc_total(begland:endland)
  real(r8) :: fpc_total_new(begland:endland)
  real(r8) :: tmomin20(begland:endland)
  real(r8) :: agdd20(begland:endland)
  real(r8) :: agddtw(begland:endland)
  real(r8) :: fpc_ind
  real(r8) :: fpcgridtemp
  real(r8) :: estab_rate
  real(r8) :: estab_grid
  real(r8) :: crownarea_max
  real(r8) :: bare_max
  real(r8) :: bare
  real(r8) :: nind_old
  real(r8) :: sm_ind_temp
  real(r8) :: stemdiam
  real(r8) :: pi
  real(r8) :: tcmin
  real(r8) :: tcmax
  real(r8) :: gddmin

  !minimum individual density for persistence of PFT (indiv/m2)
  real(r8), parameter :: nind_min = 1.0e-10 
                                            
  !minimum precip. for establishment (mm/s)
  real(r8), parameter :: prec_min_estab = 100./(365.*SHR_CONST_CDAY) 

  !maximum sapling establishment rate (indiv/m2)
  real(r8), parameter :: estab_max = 0.24   
! ----------------------------------------------------------------------

! **********************************************************************
! My version of LPJ's subr. bioclim

! Limits based on 20-year running averages of coldest-month mean
! temperature and growing degree days (5 degree base).
! For SURVIVAL, coldest month temperature and GDD should be
! at least as high as PFT-specific limits.
! For REGENERATION, PFT must be able to survive AND coldest month
! temperature should be no higher than a PFT-specific limit.

  tmomin20(:) = 0.
  agdd20(:) = 0.
  agddtw(:) = 0.
  sumwt(:) = 0.

  do k = begpatch,endpatch
     if (patchvec%wtxy(k) > 0.) then
        l = patchvec%land(k)
        tmomin20(l) = tmomin20(l) + patchvec%wtxy(k) *clm(k)%tmomin20
        agdd20(l) = agdd20(l) + patchvec%wtxy(k) * clm(k)%agdd20
        agddtw(l) = agddtw(l) + patchvec%wtxy(k) * clm(k)%agddtw
        sumwt(l) = sumwt(l) + patchvec%wtxy(k)
     endif
  end do
  do l = begland,endland
     if (sumwt(l) /= 0.) then
        tmomin20(l) = tmomin20(l)/sumwt(l)
        agdd20(l) = agdd20(l)/sumwt(l)
        agddtw(l) = agddtw(l)/sumwt(l)
     end if
  end do

! Must go thru all 16 pfts and decide which can/cannot establish or survive
! Determine present, survive, estab
! note - if tmomin20 > tcmax then  crops, shrubs and 2nd boreal
! summergreen tree cannot exist yet (see EcosystemDynini) because 
! they often coexist using up all pft patches and allowing for no bare
! ground. Thus they make fpc_grid_total < 1. Solve by allowing one 
! one more pft patch per grid cell than the number of pfts.

  do l = begland,endland
     do m = 0,numpft
        present(l,m) = .false.
        survive(l,m) = .false.
        estab(l,m)   = .false.
     end do
  end do

  do m = 1, numpft
     tcmin  = pftpar(m,28) + tfrz !PFT-specific minimum coldest-month temp.
     tcmax  = pftpar(m,29) + tfrz !PFT-specific maximum coldest-month temp.
     gddmin = pftpar(m,30)        !PFT-specific minimum GDD
     do l = begland,endland
        if (tmomin20(l) >= tcmin) then
           survive(l,m) = .true.
           if (  tmomin20(l) <= tcmax &
                .and. agdd20(l) >= gddmin &
                .and. nint(agddtw(l)) == 0) then
              estab(l,m) = .true.
           endif
        endif
     end do
  end do

! **********************************************************************

  pi = SHR_CONST_PI

  do k = begpatch,endpatch
     l = patchvec%land(k)
     present(l,clm(k)%itypveg) = clm(k)%present
     if (.not.clm(k)%present) clm(k)%itypveg = 0
  end do

! Kill PFTs not adapted to current climate, introduce newly "adapted" PFTs

! slevis: case 1 -- pft ceases to exist

  do k = begpatch,endpatch
     l = patchvec%land(k)
     if (clm(k)%present .and. &
         (.not.survive(l,clm(k)%itypveg) .or. clm(k)%nind<nind_min)) then

        clm(k)%present = .false.
        present(l,clm(k)%itypveg) = .false.

        ! Add killed biomass to litter

        if (tree(clm(k)%itypveg)) then
           clm(k)%litterag = clm(k)%litterag + clm(k)%nind * &
                             (clm(k)%lm_ind + clm(k)%sm_ind + clm(k)%hm_ind)
        else !if grass
           clm(k)%litterag = clm(k)%litterag + clm(k)%nind * clm(k)%lm_ind
        endif
        clm(k)%litterbg = clm(k)%litterbg + clm(k)%nind * clm(k)%rm_ind

        clm(k)%itypveg = 0             !need to set anything other than this?
        clm(k)%fpcgrid = 0.0
     end if
  end do 

! slevis: case 2 -- pft begins to exist

  do k = begpatch,endpatch
     l = patchvec%land(k)
     if (clm(k)%itypwat == istsoil) then
        if (.not.clm(k)%present .and. clm(k)%prec365 >= prec_min_estab) then
           do m = 1, numpft
              if (.not.clm(k)%present .and. clm(k)%itypveg /= m) then
                 if (.not.present(l,m) .and. estab(l,m)) then

                    clm(k)%present = .true.
                    present(l,m) = .true.
                    clm(k)%itypveg = m
                    
                    if (tree(m)) then
                       clm(k)%nind = 0.0
                    else
                       clm(k)%nind = 1.0 !each grass PFT = 1 "individual"
                    endif

                    clm(k)%lm_ind = 0.0
                    clm(k)%sm_ind = 0.0
                    clm(k)%rm_ind = 0.0
                    clm(k)%hm_ind = 0.0
                    clm(k)%fpcgrid = 0.0
                    
                    if (.not.tree(m)) clm(k)%crownarea = 1.0
                    
                 end if !conditions suitable for establishment
              end if    !no pft present and pft 'm' was not here until now
           end do       !numpft
        end if          !more conditions for establishment
     end if             !if soil
  end do  

! slevis: case 3 -- some pfts continue to exist (no change) and
!                   some pfts continue to not exist (no change)

! SAPLING AND GRASS ESTABLISHMENT

! Initialize; then...
! Calculate total woody FPC and number of woody PFTs present and
! able to establish

  do l = begland,endland
     ngrass(l) = 0
     npft_estab(l) = 0
     fpc_tree_total(l) = 0.0
     fpc_grass_total(l) = 0.0
     fpc_total(l)=0.0
     fpc_total_new(l) = 0.0
  end do

  do k = begpatch,endpatch
     l = patchvec%land(k)
     if (clm(k)%present) then
        if (tree(clm(k)%itypveg)) then
           fpc_tree_total(l) = fpc_tree_total(l) + clm(k)%fpcgrid
           if (estab(l,clm(k)%itypveg)) npft_estab(l) = npft_estab(l) + 1
        else if (.not.tree(clm(k)%itypveg) .and. clm(k)%itypveg > 0) then !grass
           ngrass(l) = ngrass(l) + 1
           fpc_grass_total(l) = fpc_grass_total(l) + clm(k)%fpcgrid
        endif
        fpc_total(l) = fpc_total(l) + clm(k)%fpcgrid
     endif
  enddo

! Prohibit establishment under extreme temperature or water stress.

  do k = begpatch,endpatch
     l = patchvec%land(k)

     if (clm(k)%prec365 >= prec_min_estab .and. npft_estab(l) > 0) then

! Calculate establishment rate over available space, per tree PFT
! Maximum establishment rate reduced by shading as tree FPC approaches 1
! Total establishment rate partitioned equally among regenerating woody PFTs

        estab_rate = estab_max * (1.0-exp(5.0*(fpc_tree_total(l)-1.0))) / &
                     real(npft_estab(l))

! Calculate grid-level establishment rate per woody PFT
! Space available for woody PFT establishment is proportion of grid cell
! not currently occupied by woody PFTs

        estab_grid = estab_rate * (1.0-fpc_tree_total(l))

     else !if unsuitable climate for establishment

        estab_grid = 0.0

     endif

     if (clm(k)%present .and. tree(clm(k)%itypveg) .and. estab(l,clm(k)%itypveg)) then

        crownarea_max = pftpar(clm(k)%itypveg,20)

! Add new saplings to current population

        nind_old = clm(k)%nind
        clm(k)%nind = nind_old + estab_grid

        clm(k)%lm_ind = (clm(k)%lm_ind * nind_old + &
                         lm_sapl(clm(k)%itypveg) * estab_grid) / clm(k)%nind
        sm_ind_temp   = (clm(k)%sm_ind * nind_old + &
                         sm_sapl(clm(k)%itypveg) * estab_grid) / clm(k)%nind
        clm(k)%hm_ind = (clm(k)%hm_ind * nind_old + &
                         hm_sapl(clm(k)%itypveg) * estab_grid) / clm(k)%nind
        clm(k)%rm_ind = (clm(k)%rm_ind * nind_old + &
                         rm_sapl(clm(k)%itypveg) * estab_grid) / clm(k)%nind

! Calculate height, diameter and crown area for new average
! individual such that the basic allometric relationships (A-C below)
! are satisfied.

! (A) (leaf area) = latosa * (sapwood xs area)
!        (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
! (B) (leaf mass) = lmtorm * (root mass)
! (C) height = allom2 * (stem diameter)**allom3  (source?)
! (D) (crown area) = min (allom1 * (stem diameter)**reinickerp,
!                                 crownarea_max)

! From (A),
!  (1) sap_xsa = lm_ind * sla / latosa
!  (2) wooddens = (sm_ind + hm_ind) / stemvolume
!  (3) stemvolume = stem_xsa * height
! From (1), (2) & (3),
!  (4) stem_xsa = (sm_ind + hm_ind) / wooddens / height
!  (5) stem_xsa = pi * (stemdiam**2) / 4
! From (5),
!  (6) stemdiam = ( 4 * stem_xsa / pi )**0.5
! From (4) & (6),
!  (7) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens / height /
!        pi )**0.5
! From (C) & (7),
!  (8) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens /
!        ( allom2 * stemdiam**allom3 ) / pi )**0.5
! From (8),
!  (9) stemdiam = ( 4 * (sm_ind + hm_ind ) / wooddens / pi /
!        allom2 )**( 1 / (2 + allom3) )

        stemdiam = (4.0*(sm_ind_temp+clm(k)%hm_ind)/wooddens/pi/allom2)** &
                   (1.0/(2.0+allom3))              !Eqn 9
        clm(k)%htop = allom2 * stemdiam**allom3  !Eqn C
        clm(k)%crownarea = min(crownarea_max, allom1*stemdiam**reinickerp)!Eqn D

! Recalculate sapwood mass, transferring excess sapwood to heartwood
! compartment, if necessary to satisfy Eqn A

        clm(k)%sm_ind = clm(k)%lm_ind * clm(k)%htop * wooddens * &
                        sla(clm(k)%itypveg) / latosa
        clm(k)%hm_ind = clm(k)%hm_ind + (sm_ind_temp-clm(k)%sm_ind)

! Update LAI and FPC

        if (clm(k)%crownarea > 0.0) then
           clm(k)%lai_ind = clm(k)%lm_ind * sla(clm(k)%itypveg) / clm(k)%crownarea
        else
           clm(k)%lai_ind = 0.0
        endif

        fpc_ind = 1.0 - exp(-0.5*clm(k)%lai_ind)
        clm(k)%fpcgrid = clm(k)%crownarea *clm(k)%nind * fpc_ind

        fpc_total_new(l) = fpc_total_new(l) + clm(k)%fpcgrid

     endif
  enddo

  do k = begpatch,endpatch
     l = patchvec%land(k)
     if (fpc_total_new(l) > 0.95) then
        if (tree(clm(k)%itypveg) .and. clm(k)%present) then
           nind_old = clm(k)%nind
           clm(k)%nind = clm(k)%nind / (fpc_total_new(l)/0.95)
           clm(k)%fpcgrid = clm(k)%fpcgrid / (fpc_total_new(l)/0.95)
           clm(k)%litterag = clm(k)%litterag + (nind_old-clm(k)%nind) * &
                             (clm(k)%lm_ind+clm(k)%sm_ind+clm(k)%hm_ind)
           clm(k)%litterbg = clm(k)%litterbg + (nind_old-clm(k)%nind) * &
                             clm(k)%rm_ind
        endif
        fpc_total(l) = 0.95
     endif
  enddo

! SECTION FOR GRASSES
 
  do k = begpatch,endpatch
     l = patchvec%land(k)
     if (clm(k)%present .and. .not.tree(clm(k)%itypveg)) then

        if (estab(l,clm(k)%itypveg)) then
           ! Grasses can establish in non-vegetated areas
           if (ngrass(l) > 0) then
              bare = (1.0-fpc_total(l)) / real(ngrass(l))
           else
              bare = 0.0
           endif

           bare_max = (-2.0 * clm(k)%crownarea *                     &
                       log(max(1.0-bare-clm(k)%fpcgrid, 0.000001)) / &
                       sla(clm(k)%itypveg) -                             &
                       clm(k)%lm_ind) /                              &
                       lm_sapl(clm(k)%itypveg)

           bare = max(0.0, min(bare, bare_max))

           clm(k)%lm_ind = clm(k)%lm_ind + bare * lm_sapl(clm(k)%itypveg)
           clm(k)%rm_ind = clm(k)%rm_ind + bare * rm_sapl(clm(k)%itypveg)

        endif

        if (clm(k)%lm_ind <= 0.0) then
           clm(k)%present = .false.
           clm(k)%litterbg = clm(k)%litterbg + clm(k)%rm_ind * clm(k)%nind
        endif
     endif
  enddo

! recalculate fpc's and do error check

  do l = begland,endland
     fpc_total(l)=0.0
     fpc_total_new(l) = 0.0
  end do

  do k = begpatch,endpatch
     l = patchvec%land(k)
     if (clm(k)%present) then
        if (clm(k)%crownarea > 0.0) then
           clm(k)%lai_ind = clm(k)%lm_ind * sla(clm(k)%itypveg) / clm(k)%crownarea
        else
           clm(k)%lai_ind = 0.0
        endif

        fpc_ind = 1.0 - exp(-0.5*clm(k)%lai_ind)
        clm(k)%fpcgrid = clm(k)%crownarea * clm(k)%nind * fpc_ind
        fpc_total(l) = fpc_total(l) + clm(k)%fpcgrid
     else
        clm(k)%fpcgrid = 0.0
     endif
  enddo

  do k = begpatch,endpatch    
     l = patchvec%land(k)
     if (fpc_total(l) < 0.0) then
        write(6,*) 'Error in Establishment: fpc_total is',fpc_total(l),l
        call endrun
     end if
  end do

! Adjustment b/c fpc_total > 1. This can happen, because grasses are allowed
! to grow on the bare soil as determined before updating the trees
! This may imply that LPJ allows a brief coexistence of some
! trees and grasses. LSM cannot allow such coexistence. (slevis)

  do k = begpatch,endpatch
     l = patchvec%land(k)
     if (fpc_total(l) > 1.0) then
        if (fpc_total(l) > 1.1) then
           write(6,*) 'Error in Establishment: fpc_total is',fpc_total(l),l
           call endrun
        end if
        if (clm(k)%itypveg >= 12 .and. clm(k)%fpcgrid > 0.0) then
           fpcgridtemp = clm(k)%fpcgrid
           clm(k)%fpcgrid = max(0.0, clm(k)%fpcgrid-(fpc_total(l)-1.0))
           if (clm(k)%fpcgrid == 0.0) then
              clm(k)%present = .false.
              clm(k)%litterag = clm(k)%litterag + clm(k)%nind * clm(k)%lm_ind
              clm(k)%litterbg = clm(k)%litterbg + clm(k)%nind * clm(k)%rm_ind
           end if
           fpc_total(l) = fpc_total(l) - fpcgridtemp + clm(k)%fpcgrid
        end if
     end if
  end do

! Next lines avoid fpcgrid=0 with present=.T. and
! itypveg=0 with lai>0. Should C pools be reset, to balance carbon?

  do k = begpatch,endpatch
     l = patchvec%land(k)
     if (clm(k)%fpcgrid == 0.) clm(k)%present = .false.
     if (.not. clm(k)%present) clm(k)%lai_ind = 0.0
     if (fpc_total(l) < 1.0) then
        if (clm(k)%itypveg == noveg) then
           clm(k)%fpcgrid = 1.0 - fpc_total(l)
           fpc_total(l) = 1.0
        end if
     end if
     fpc_total_new(l) = fpc_total_new(l) + clm(k)%fpcgrid
  end do

  do l = begland,endland
     if (abs(fpc_total_new(l) - 1.0) > 1.0e-6) then
        write(6,*) 'Error in Establishment: fpc_total_new =',fpc_total_new(l),l
        call endrun
     end if
  end do

  return
end subroutine Establishment
