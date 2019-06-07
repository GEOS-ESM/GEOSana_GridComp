
      subroutine tropopause (alpha,im,jm,levels,ple,pl,
     .                   Ptop,TropP,TmpU,TropT,SphU,TropQ )

C********************************************************************
C!                     Subroutine Tropfind
C********************************************************************
C
C!ROUTINE:  Tropfind
C
C!DESCRIPTION:
C
C     This routine finds the tropopause pressure and temperature from
C     model temperature profiles.  It finds these values for a
C     given 3-dimensional grid.  The algorithm is based on the
C     temperature profile only.  It is similar to visually inspecting
C     a Skew-T Log P diagram for the left-most kink in the temperature
C     profile.  In a standard Skew-T Log P diagram the isotherms
C     intersect the Log P scale at a 45 degree angle.  To more
C     effectively isolate the tropopause,  this angle can be adjusted.
C     ThatÕs how this program works.  It adjusts the angle at which
C     the isotherms intersect the Log P axis via the coeffiecient
C     Òalpha.Ó It simply looks for the furthest point to the left on
C     the profile.  The routine defines the leftward position of each
C     temperature point as xfact:
C
C         xfact = (alpha * temperature) - log (pressure)
C
C     The tropopause for a given profile is where xfact is a minimum.
C     Uncertainty can occur when the kink is not very distinct.  This
C     is where the selection of alpha becomes important.  Optimal
C     values for alpha appear to be between .02 and .04.  For lower
C     values of alpha, the tropopause selection will favor lower
C     levels (higher P).  For higher values of alpha, the selection
C     will favor higher levels (lower P).  A value of alpha = .03
C     appears to be optimal in generating tropopause values within
C     the range indicated in the Handbook of Geophysics and the Space
C     Environment, AFGL, chapt 14, 1985.
C
C!INPUT PARAMETERS:
C
C     alpha   = see discussion above ((log mb)/deg K)
C     im      = number of longitude grid points
C     jm      = number of latitude grid points
C     levels  = number of model levels
C     ple     = pressures at model edges
C     pl      = pressures at model mid-layers
C     Ptop    = pressure at the top model edge (mb)
C     TmpU    = 3-d array of gridded temperature on model level (deg K)
C
C!OUTPUT PARAMETERS:
C
C     TropP = array of gridded tropopause pressures (mb)
C     TropT = array of gridded tropopause temperatures (deg K)
C
C!REVISION HISTORY:
C
C     Created 25 Jun 97 by Jim Stobie
C     Mods:   23 Jul 97 by L.Takacs
C          1) Send in Ps-Ptop (rather than Ps)
C          2) Use  log10 (average pressure) rather than average 
C             (log10(pressure) )
C    05 Feb 2002 by da Silva
C             Eliminated getcon() dependecy; ps_top not needed here
C             Change order of pl/ple indices (i,j,k)->i,k,j)
C
C********************************************************************

      implicit none

      integer    levmax
      parameter (levmax=1000)


c Passed Variables

      integer, intent(in) :: im
      integer, intent(in) :: jm
      integer, intent(in) :: levels

      real, intent(in)  :: alpha
      real, intent(in)  :: Ptop
      real, intent(in)  :: pl(im,levels,jm)    ! Pressure at model mid-layers
      real, intent(in)  :: ple(im,levels+1,jm)  ! Pressure at model edges
      real, intent(in)  :: TmpU(im,jm,levels)
      real, intent(in)  :: SphU(im,jm,levels)
      real, intent(out) :: TropP(im,jm)
      real, intent(out) :: TropT(im,jm)
      real, intent(out) :: TropQ(im,jm)
!!!      integer, intent(out), OPTIONAL :: TropK(im,jm) ! tropopause index

c Local Variables

      integer i,j,k          !loop variables
      integer kend           !end index for model level search
      integer kstart         !start index for model level search
      integer ktrop          !index for tropopause level

      real phigh             !highest pressure (mb) for search
      real plow              !lowest pressure (mb) for search
      real undef             !value for undefined variables
      real xfacmn            !minimum x-factor, see prologue
      real xfact(levmax)     !x-factor, see prologue

      real getcon            !constant function

c----------------------------------------------------------------
c Get value for undefined variables from funciton getcon
c---------------------------------------------------------------

!!!      undef = getcon('UNDEF')
       undef = 1.E15

c----------------------------------------------------------------
c Set vertical limits on search.  Tropopause search will be
c limited to the range between plow and phigh (mb).
c According to Handbook of Geophysics and the Space Environment,
c AFGL, 1985, pg 14-6, the lowest tropopause values are near 8 km
c in the polar winter (approx 350 mb) and the highest near 18 km
c in the tropics (approx 80 mb).
c----------------------------------------------------------------

        plow  = 40.
        phigh = 550.

c----------------------------------------------------------------
c If this particular run does not go up to plow, then write a
c warning message and fill TropT and TropP with undefined data
c values.
c----------------------------------------------------------------

      if (plow.LT.Ptop) then

!!!	write(6,1000) Ptop,plow

        do j = 1, jm
	 do i = 1, im
!!!	    if ( present(tropk) ) TropK(i,j) = -1
	    TropT(i,j) = undef
	    TropP(i,j) = undef
            TropQ(i,j) = undef

        end do
        end do

      else

c--------------------
c Loop over lat/lon Grid
c--------------------

        do j = 1, jm
        do i = 1, im

c-------------------------------------------------------
c Find pressure range for search.  Search will begin
c at first model level edge above phigh, and end
c at the first model level edge below plow.
c-------------------------------------------------------

               kend = levels
            do while (ple(i,kend,j).GE.phigh)
               kend = kend-1
            enddo

               kstart = 1
            do while (ple(i,kstart,j).le.plow)
               kstart = kstart+1
            enddo

c-----------------------------------------------------
c Calculate pressure of the model layer midpoints.
c Then calculate xfact for these points.  See prologue
c for description of xfact.
c-----------------------------------------------------

            do k = kstart, kend
              xfact(k) = alpha * TmpU(i,j,k) - log10(pl(i,k,j))
            end do

c-----------------------------------------------
c Tropopause is level for which xfact is minimum
c-----------------------------------------------

            xfacmn = 100000.

            do k = kstart, kend

              if (xfact(k).LT.xfacmn) then
                xfacmn = xfact(k)
                ktrop = k
              end if

            end do
c-------------------------------------------------------
c If the minimum value of xfact is at the upper or lower
c boundary of the search, then the tropopause has not
c been sucessfully isolated.  In this case a warning
c message is printed and the grid point value filled
c with the undefined value.
c-------------------------------------------------------

            if(ktrop.EQ.kstart.or.ktrop.EQ.kend) then

c             write(6,1100) i,j
              !!! if ( present(tropk) ) TropK(i,j) = -1
              tropp(i,j) = undef
              tropt(i,j) = undef
              TropQ(i,j) = undef

            else

c------------------------------------------------------
c If the tropopause has been successfully isolated
c     store tropopause pressure    in TropP
c and store tropopause temperature in TropT.
c------------------------------------------------------

              !!! if ( present(tropk) ) tropk(i,j) =  ktrop
              tropp(i,j) =   pl(i,ktrop,j)
              tropt(i,j) = TmpU(i,j,ktrop)
              TropQ(i,j) = SphU(i,j,ktrop)

            end if

          end do

        end do

      end if

      return

 1000 format("WARNING.  Can't find tropopause.  Top of model",
     *       " is too low (",f4.0," mb).  Tropfind requires",
     *       " model top to be ",f4.0," mb")

 1100 format("WARNING.  Can't find tropopause at gridpoint:",
     *       "  i = ",i4,"  j = ",i4)

      end


