#include <params.h>
      subroutine mgwinti (cpx, cpwv, gx, rx)
!-----------------------------------------------------------------------
! Time independent initialization for multiple gravity wave parameterization.
!-----------------------------------------------------------------------
      implicit none

#include <pmgrid.h>
#include <mgw.h>

! Input variables
      real 
     $     cpx,                 ! specific heat of dry air (constant p)
     $     cpwv,                ! specific heat of water vapor (constant p)
     $     gx,                  ! acceleration of gravity
     $     rx                   ! gas constant for dry air
 
! Local variables
 
      integer k
      real oroe                 ! efficiency factor for orographic term

      integer ks
      real ptop, pint, pedge
      real ak(plev+1), bk(plev+1)
 
! Copy model constants
      cp     = cpx
      cpvir  = cpwv / cp - 1
      g      = gx
      r      = rx
c+
c Set MGWD constants
c-
! SJL kwv    = 6.28e-5          ! 100 km wave length
      kwv    = 6.28e-5 * min( 4., max(1.,  (plat-1)/90.))
      dback  = 0.05             ! background diffusivity
      orofc2 = OROFC2           ! critical froude number squared
      oroe   = OROE             ! efficiency for orographic generation
      effgw  = EFFGW		! efficiency of the tendencies
      tauscal= TAUSCAL          ! scale factor for background stress
      taubgnd= TAUBGND          ! background stress amplitude

!SJL's mod: 10/30/2001
!=========================
#if ( !defined GW_PRECIP )
      taubgnd= 0.7*taubgnd
#endif
!=========================

      fracldv= FRACLDV          ! fraction of tau0 diverged in low level region
      zldvcon= 10.              ! constant for determining zldv
c+
c Set phase speeds 
c-
      do k = -pgwv, pgwv
         c(k)   = 10. * k       ! 0, +/- 10, +/- 20, ... m/s
      end do

      if( masterproc ) then
          print *, ' '
          print *, 'mgwinti: pgwv = ', pgwv
          print *, 'mgwinti: c(l) = ', c
          print *, ' '
      endif
c+
c Set radiative damping times
c-
      do k = 0, plev
         alpha(k) = 1.e-6       ! about 10 days.
      end do
c+
c Min and max values to keep things reasonable
c-
      mxasym = MXASYM           ! max factor of 10 from |tau(c)| to |tau(-c)|
      mxrange= MXRANGE          ! factor of 100 from max to min |tau(c)|
      n2min  = 1.e-8
      orohmin= 10.              ! min sdv for orographic wave drag
      orovmin=  2.              ! min wind speed for orographic wave drag
      taumin = 1.e-10
      tndmin = 1.e-10
      tndmax = 500. / 86400.    ! max permitted tendency (500 m/s/day)
      ubmc2mn= 0.01

! Determine efficiency factor from gravity wave parameters

      oroeko2 = 0.5 * oroe * kwv 
      efcncy  = oroe * orofc2
c+
c Determine other derived constants
c-
      rog    = r/g
      efkw   = efcncy * kwv
c+
c Determine the bounds of the background and orographic stress regions
c-

      ktopbg  = 0
      kbotoro = plev

      call set_eta(plev, ks, ptop, pint, ak, bk)

      do k = 0, plev
         pedge = ak(k+1) + bk(k+1)*1.e5
         if ( pedge < 10000.) kbotbg  = k       ! 100 mb
         if ( pedge <   100.) ktoporo = k
      end do

      if ( ak(1) >= 30. ) ktoporo = 0

      if( masterproc ) then
	write (6,*) 'Top layer for background GWD     =', ktopbg
	write (6,*) 'Bottom layer for background GWD  =', kbotbg
	write (6,*) 'Top layer for ORO GWD            =', ktoporo
	write (6,*) 'Bottom layer for ORO GWD         =', kbotoro
	write (6,*) ' '
      endif

      return
      end
