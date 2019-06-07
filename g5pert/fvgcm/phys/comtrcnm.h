c
c $Id$
c $Author$
c
C
C Lists of tracer names and diagnostics
C
      common /comtrcnm/ tracnam(pcnst+pnats), hadvnam(pcnst),
     $                  vadvnam(pcnst), vdiffnam(pcnst),
     $                  dcconnam(pcnst), fixcnam(pcnst),
     $                  tendnam(pcnst), sflxnam(pcnst), srcnam(pcnst),
     $                  tottnam(pcnst), tsnam(plevmx)
C
      character*8 tracnam   ! tracer names (including any non-advected)
      character*8 hadvnam   ! names of horizontal advection tendencies
      character*8 vadvnam   ! names of vertical advection tendencies
      character*8 vdiffnam  ! names of v-diff tendencies
      character*8 dcconnam  ! names of convection tendencies
      character*8 fixcnam   ! names of species slt fixer tendencies
      character*8 tendnam   ! names of total tendencies of species
      character*8 sflxnam   ! names of surface fluxes of species
      character*8 srcnam    ! names of source/sink tendencies of species
      character*8 tottnam   ! names for horz + vert + fixer tendencies
      character*8 tsnam     ! names of sub-surface temperature fields
C
 
