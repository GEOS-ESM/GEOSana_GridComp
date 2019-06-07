#include <params.h>
      subroutine kdpfnd(pkdim   ,pmap    ,sig     ,sigdp   ,kdpmap  ,
     $                  kdp     )
C-----------------------------------------------------------------------
C
C Determine vertical departure point indices that point into a grid
C containing the full or half sigma levels.  Use an artificial evenly 
C spaced vertical grid to map into the true model levels.
C
C Note:
C Indices are computed assuming the the sigdp values have
C been constrained so that sig(1) .le. sigdp(i,j) .lt. sig(pkdim).
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C Reviewed:          D. Williamson, P. Rasch, March 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
      integer pgls
      parameter ( pgls = plon*plev ) ! length of latitude slice
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer pkdim             ! dimension of "sig"
      integer pmap              ! dimension of "kdpmap"
C
      real sig  (pkdim)         ! vertical grid coordinates
      real sigdp(pgls)          ! vertical coords. of departure points
C
      integer kdpmap(pmap)      ! array of model grid indices which
C                               ! are mapped into the artificial grid.
C
C Output arguments
C
      integer kdp(pgls)         ! vertical index for each dep. pt.
C
C---------------------------Local variables-----------------------------
C
      integer i                 ! index
      integer ii                ! index
C
      real rdel                 ! recip. of interval in artificial grid
      real sig1ln               ! ln (sig(1))
C
C-----------------------------------------------------------------------
C
      rdel   = float(pmap)/( log(sig(pkdim)) - log(sig(1)) )
      sig1ln = log( sig(1) )
C
      do i = 1,pgls
C
C First guess of the departure point's location in the model grid
C
         ii = max0(1,min0(pmap,int((log(sigdp(i))-sig1ln)*rdel + 1.)))
         kdp(i) = kdpmap(ii)
C
C Determine if location is in next interval
C
         if(sigdp(i) .ge. sig(kdp(i)+1)) then
            kdp(i) = kdp(i) + 1
         end if
      end do
C
      return
      end
 
