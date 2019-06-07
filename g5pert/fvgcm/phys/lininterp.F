      subroutine lininterp(arrin, yin, nlev, nlatin,
     &                     arrout, yout, nlatout)
C-----------------------------------------------------------------------
C
C Do a linear interpolation from input mesh defined by yin to output
C mesh defined by yout.  Where extrapolation is necessary, values will
C be copied from the extreme edge of the input grid.  Vectorization is over
C the vertical (nlev) dimension.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
!
! Input arguments
!
      integer nlev               ! number of vertical levels
      integer nlatin             ! number of input latitudes
      integer nlatout            ! number of output latitudes
      real arrin(nlev,nlatin)    ! input array of values to interpolate
      real yin(nlatin)           ! input mesh
      real yout(nlatout)         ! output mesh
!
! Output arguments
!
      real arrout(nlev,nlatout)  ! interpolated array
!
! Local workspace
!
      integer j, jj              ! latitude indices
      integer js, jn, jjprev     ! latitude indices
      integer k                  ! level index
      integer icount             ! number of values

      real extrap                ! percent grid non-overlap
!
! Dynamic
!
      pointer (pjjm,jjm)
      pointer (pjjp,jjp)
      pointer (pwgts,wgts)
      pointer (pwgtn,wgtn)

      integer jjm(nlatout)       ! interp. indices south
      integer jjp(nlatout)       ! interp. indices north
      real wgts(nlatout)         ! interp. weight south
      real wgtn(nlatout)         ! interp. weight norty
!
!------------------------------------------------------------------------------
!
      call getmem('LININTERP', nlatout, pjjm)
      call getmem('LININTERP', nlatout, pjjp)
      call getmem('LININTERP', nlatout, pwgts)
      call getmem('LININTERP', nlatout, pwgtn)
!
! Check validity of input coordinate arrays: must be monotonically increasing,
! and have a total of at least 2 elements
!
      if (nlatin.lt.2) then
        write(6,*)'LININTERP: Must have at least 2 input points for',
     &            ' interpolation'
        call endrun
      end if
      icount = 0
      do j=1,nlatin-1
        if (yin(j).gt.yin(j+1)) icount = icount + 1
      end do
      do j=1,nlatout-1
        if (yout(j).gt.yout(j+1)) icount = icount + 1
      end do
      if (icount.gt.0) then
        write(6,*)'Non-monotonic coordinate array(s) found'
        call endrun
      end if
!
! Initialize index arrays for later checking
!
      do j=1,nlatout
        jjm(j) = 0
        jjp(j) = 0
      end do
!
! For values which extend beyond N and S boundaries, set weights
! such that values will just be copied.
!
      do js=1,nlatout
        if (yout(js).gt.yin(1)) goto 10
        jjm(js) = 1
        jjp(js) = 1
        wgts(js) = 1.
        wgtn(js) = 0.
      end do
   10 do jn=nlatout,1,-1
        if (yout(jn).le.yin(nlatin)) goto 20
        jjm(jn) = nlatin
        jjp(jn) = nlatin
        wgts(jn) = 1.
        wgtn(jn) = 0.
      end do
!
! Loop though output indices finding input indices and weights
!
   20 jjprev = 1
      do j=js,jn
        do jj=jjprev,nlatin-1
          if (yout(j).gt.yin(jj) .and. yout(j).le.yin(jj+1)) then
            jjm(j) = jj
            jjp(j) = jj + 1
            wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
            wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
            goto 30
          end if
        end do
        write(6,*)'LININTERP: Failed to find interp values'
   30   jjprev = jj
      end do
!
! Check grid overlap
!
      extrap = 100.*((js - 1.) + float(nlatout - jn))/nlatout
      if (extrap.gt.30.) then
         write(6,*)'********LININTERP WARNING:',extrap,' % of output',
     &             ' grid will have to be extrapolated********'
      end if
!
! Check that interp/extrap points have been found for all outputs
!
      icount = 0
      do j=1,nlatout
        if (jjm(j).eq.0 .or. jjp(j).eq.0) icount = icount + 1
      end do
      if (icount.gt.0) then
        write(6,*)'LININTERP: Point found without interp indices'
        call endrun
      end if
!
! Do the interpolation
!
      do j=1,nlatout
        do k=1,nlev
          arrout(k,j) = arrin(k,jjm(j))*wgts(j) + 
     &                  arrin(k,jjp(j))*wgtn(j)
        end do
      end do

      call freemem(pjjm)
      call freemem(pjjp)
      call freemem(pwgts)
      call freemem(pwgtn)

      return
      end
 
