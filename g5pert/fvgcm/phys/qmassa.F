#include <params.h>
      subroutine qmassa(cwava   ,w       ,q3      ,pdel    ,hw1lat  )
C-----------------------------------------------------------------------
C
C Calculate contribution of current latitude to mass of constituents
C being advected by slt.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          P. Rasch, D. Williamson, August 1992
C Modified:          B. Eaton and J. Olson, 1994 and 1995
C Reviewed:          P. Rasch, D. Williamson, April 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real cwava                ! normalization factor    l/(g*plon)
      real w                    ! gaussian weight this latitude
      real q3(plond,plev,pcnst) ! constituents
      real pdel(plond,plev)     ! pressure diff between interfaces
C
C Output arguments
C
      real hw1lat(pcnst)        ! accumulator
C
C---------------------------Local variables-----------------------------
C
      integer i,k,m             ! longitude, level, constituent indices
      real const
C
C-----------------------------------------------------------------------
C
C Integration factor (the 0.5 factor arises because gaussian weights
C sum to 2)
C
      const = cwava*w*0.5
C
C Compute mass integral for water ONLY
C
      do k=1,plev
         do i=1,plon
            hw1lat(1) = hw1lat(1) + q3(i,k,1)*pdel(i,k)
         end do
      end do
C
C Compute mass integral for non-water constituents (on a DRY basis)
C
      do m=2,pcnst
         do k=1,plev
            do i=1,plon
               hw1lat(m) = hw1lat(m) +
     $                       q3(i,k,m)*(1. - q3(i,k,1))*pdel(i,k)
            end do
         end do
      end do
C
      do m = 1,pcnst
         hw1lat(m) = hw1lat(m)*const
      end do
      return
      end
 
