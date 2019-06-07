#include <params.h>
      subroutine radozn(lat,pmid,o3vmr)
C--------------------------------------------------------------------------
C
C Interpolate ozone from current time-interpolated values to model levels
C
C Bruce Briegleb  29 September 1992
C Reviewed: J. Kiehl April 1996
C Reviewed: B. Briegleb  May 1996
C
C--------------------------------------------------------------------------
c
c $Id$
c $Author$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comozp.h>
C--------------------------------------------------------------------------
C
C input parameters
C
      integer lat             ! latitude index
      real pmid(plond,plev)   ! level pressures (mks)
C
C output parameters
C
      real o3vmr(plond,plev)  ! ozone volume mixing ratio
C
C local storage
C
      integer i               ! longitude index
      integer k,kk,kkstart    ! level indices
      integer kupper(plond)   ! Level indices for interpolation
      integer kount           ! Counter

      real dpu                ! upper level pressure difference
      real dpl                ! lower level pressure difference
C
C Initialize index array
C
      do i=1,plon
        kupper(i) = 1
      end do
      do k=1,plev
C
C Top level we need to start looking is the top level for the previous k
C for all longitude points
C
        kkstart = pozlev
        do i=1,plon
          kkstart = min0(kkstart,kupper(i))
        end do
        kount = 0
C
C Store level indices for interpolation
C
        do kk=kkstart,pozlev-1
          do i=1,plon
            if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
              kupper(i) = kk
              kount = kount + 1
            end if
          end do
C
C If all indices for this level have been found, do the interpolation and
C go to the next level
C
          if (kount.eq.plon) then
            do i=1,plon
              dpu = pmid(i,k) - pin(kupper(i))
              dpl = pin(kupper(i)+1) - pmid(i,k)
              o3vmr(i,k) = (ozmix(kupper(i)  ,lat)*dpl +
     $                      ozmix(kupper(i)+1,lat)*dpu)/(dpl + dpu)
            end do
            goto 35
          end if
        end do
C
C If we've fallen through the kk=1,pozlev-1 loop, we cannot interpolate and
C must extrapolate from the bottom or top ozone data level for at least some
C of the longitude points.
C
        do i=1,plon
          if (pmid(i,k) .lt. pin(1)) then
            o3vmr(i,k) = ozmix(1,lat)*pmid(i,k)/pin(1)
          else if (pmid(i,k) .gt. pin(pozlev)) then
            o3vmr(i,k) = ozmix(pozlev,lat)
          else
            dpu = pmid(i,k) - pin(kupper(i))
            dpl = pin(kupper(i)+1) - pmid(i,k)
            o3vmr(i,k) = (ozmix(kupper(i)  ,lat)*dpl +
     $                    ozmix(kupper(i)+1,lat)*dpu)/(dpl + dpu)
          end if
        end do
        if (kount.gt.plon) then
          write(6,*)'RADOZN: Bad ozone data: non-monotonicity suspected'
          call endrun
        end if
   35   continue
      end do
      return
      end
 
