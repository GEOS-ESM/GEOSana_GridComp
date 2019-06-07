#include <params.h>
      subroutine vdinti(cpwvx   ,cpairx  ,gravx   ,rairx, ak, bk )
C-----------------------------------------------------------------------
C
C Initialization of time independent fields for vertical diffusion.
C Call initialization routine for boundary layer scheme.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
C Reviewed:          P. Rasch, B. Boville, August 1992
C Reviewed:          P. Rasch, April 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comvd.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real cpwvx        ! spec. heat of water vapor at const. pressure
      real cpairx       ! specific heat of dry air
      real gravx        ! acceleration due to gravity
      real rairx        ! gas constant for dry air

C   SJL: 10/25/98
      real ak(plev+1)
      real bk(plev+1)
C---------------------------Local workspace-----------------------------
C
      integer k         ! vertical loop index
      real    pe_ref    ! reference edge pressure for determine pbl ht.
C
C-----------------------------------------------------------------------
C
C Hard-wired numbers.
C zkmin = minimum k = kneutral*f(ri)
C
      zkmin = 0.01
C
C Set physical constants for vertical diffusion and pbl:
C
      cpair  = cpairx
      gravit = gravx
      rair   = rairx
C
C Derived constants
C ntopfl = top level to which v-diff is applied
C npbl = max number of levels (from bottom) in pbl
C
      npbl   = 9                ! original CCM3 L18 setting
C Modified by SJL to accomodate higher vertical resolution.
      do k=plev,1,-1
            pe_ref = ak(k) + bk(k)*100000.
            if(masterproc) write(6,*) k, ' pe_ref=',pe_ref, ak(k), bk(k)
         if(pe_ref .le. 40000.) then
            npbl = plevp - k
            goto 10
         endif
      enddo
10    continue
      if(masterproc) write(*,*) 'npbl= ', npbl

      ntopfl = 1                ! original CCM3 L18 setting
C SJL     
c     ntopfl = 2
      if (plev.eq.1) ntopfl = 0
      if(masterproc) write(*,*) 'ntopfl= ', ntopfl

      cpvir = cpwvx/cpairx - 1.
C
C Set the square of the mixing lengths.
C
      ml2(1) = 0.
      do k=2,plev
        ml2(k) = 30.0**2
      end do
      ml2(plevp) = 0.
C
C Initialize pbl variables
C
      call pbinti(gravx)
C
      return
      end
 
