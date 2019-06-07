#include <misc.h>
#include <preproc.h>

subroutine Light ()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: Called once per year
! 
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine light)
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use clmtype
  use clm_varder
  use clm_varmap, ONLY : patchvec, begpatch, endpatch, begland, endland
  use pft_varcon, ONLY : tree, sla
  implicit none

! ------------------------ local variables -----------------------------
  integer  :: l,k
  real(r8) :: fpc_tree_total(begland:endland)
  real(r8) :: fpc_inc_tree(begland:endland)
  real(r8) :: fpc_grass_total(begland:endland)
  integer  :: ntree(begland:endland)
  real(r8) :: excess
  real(r8) :: nind_kill
  real(r8) :: lm_old
  real(r8) :: lm_kill
  real(r8) :: rm_kill
  real(r8) :: lai_ind
  real(r8) :: fpc_ind
  real(r8), parameter :: fpc_tree_max = 0.95  !maximum total tree FPC
! ----------------------------------------------------------------------

! Initialize; then ...
! Calculate total woody FPC, FPC increment and grass cover (= crown area)

  do l = begland,endland
     fpc_tree_total(l) = 0.0
     fpc_inc_tree(l) = 0.0
     fpc_grass_total(l) = 0.0
     ntree(l) = 0
  end do

  do k = begpatch,endpatch
     l = patchvec%land(k)
     if (clm(k)%present) then
        if (tree(clm(k)%itypveg)) then
           ntree(l) = ntree(l) + 1
           fpc_tree_total(l) = fpc_tree_total(l) + clm(k)%fpcgrid
           fpc_inc_tree(l) = fpc_inc_tree(l) + clm(k)%fpcinc
        else    !if grass
           fpc_grass_total(l) = fpc_grass_total(l) + clm(k)%fpcgrid
        endif
     endif
  enddo

! LIGHT COMPETITION

  do k = begpatch,endpatch
     l = patchvec%land(k)

     if (clm(k)%present) then
        if (tree(clm(k)%itypveg)) then

           if (fpc_tree_total(l) > fpc_tree_max) then    ! case (1)

              if (fpc_inc_tree(l) > 0.0) then
                 excess = (fpc_tree_total(l) - fpc_tree_max) * &
                      clm(k)%fpcinc / fpc_inc_tree(l)
              else
                 excess = (fpc_tree_total(l) - fpc_tree_max) / &
                          real(ntree(l))
              endif

              ! Reduce individual density (and thereby gridcell-level biomass)
              ! so that total tree FPC reduced to 'fpc_tree_max'

              nind_kill = clm(k)%nind * excess / clm(k)%fpcgrid
              clm(k)%nind = clm(k)%nind - nind_kill

              ! Transfer lost biomass to litter

              clm(k)%litterag = clm(k)%litterag + nind_kill * &
                   (clm(k)%lm_ind + clm(k)%sm_ind + clm(k)%hm_ind)
              clm(k)%litterbg = clm(k)%litterbg + nind_kill * clm(k)%rm_ind

           endif

        else !if grass

           if (fpc_grass_total(l) > &
                (1.0-min(fpc_tree_total(l), fpc_tree_max))) then

              ! grass competes with itself if total fpc exceeds 1
              ! NEEDS COMMENTS !!!!!!????!!!!

              excess = (min(fpc_tree_total(l), fpc_tree_max) + &
                   fpc_grass_total(l) - 1.0) *            &
                   clm(k)%fpcgrid / fpc_grass_total(l)
              lm_old = clm(k)%lm_ind
              clm(k)%lm_ind = -2.0 * alog(1.0-(clm(k)%fpcgrid - excess)) / &
                   sla(clm(k)%itypveg)
              lm_kill = lm_old - clm(k)%lm_ind
              rm_kill = clm(k)%rm_ind * lm_kill/lm_old
              clm(k)%rm_ind = clm(k)%rm_ind - rm_kill
              
              ! Transfer lost biomass to litter

              clm(k)%litterag = clm(k)%litterag + lm_kill
              clm(k)%litterbg = clm(k)%litterbg + rm_kill
           endif

        endif

        ! update fpc (for establishment routine)
        ! slevis: lai_ind is local here

        if (clm(k)%crownarea > 0.0) then
           lai_ind = clm(k)%lm_ind * sla(clm(k)%itypveg) / clm(k)%crownarea
        else
           lai_ind = 0.0
        endif
        
        fpc_ind = 1.0 - exp(-0.5*lai_ind)
        clm(k)%fpcgrid = clm(k)%crownarea * clm(k)%nind * fpc_ind
        
     endif ! present
  enddo   !patch loop

  return
end subroutine Light








