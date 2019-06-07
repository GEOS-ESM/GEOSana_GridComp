#ifdef linux
#  define REAL16 double precision
#else
#  define REAL16 real*16
#endif
      subroutine bsslzr(bes,n)
C-----------------------------------------------------------------------
C
C Return n zeros (or if n>50, approximate zeros), of the Bessel function
C j0,in the array bes. The first 50 zeros will be given exactly, and the
C remaining zeros are computed by extrapolation,and therefore not exact.
C
C Modified 1/23/97 by Jim Rosinski to use REAL16 arithmetic
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, D. Williamson, August 1992
C Reviewed:          J. Hack, D. Williamson, April 1996
C
C-----------------------------------------------------------------------
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer n              ! Number of zeros to return
C
C Output arguments
C
      REAL16 bes(n)         ! Array containing zeros of j0
C
C---------------------------Local workspace-----------------------------
C
      REAL16 one
      parameter (one = 1.)

      integer j,nn           ! loop indices
      REAL16 pi             ! 3.14.......
      REAL16 bz(50)         ! table of first 50 zeros
      save bz                ! ensure re-entrancy
C
      data bz           / 2.4048255577,   5.5200781103,
     $    8.6537279129,  11.7915344391,  14.9309177086,  18.0710639679,
     $   21.2116366299,  24.3524715308,  27.4934791320,  30.6346064684,
     $   33.7758202136,  36.9170983537,  40.0584257646,  43.1997917132,
     $   46.3411883717,  49.4826098974,  52.6240518411,  55.7655107550,
     $   58.9069839261,  62.0484691902,  65.1899648002,  68.3314693299,
     $   71.4729816036,  74.6145006437,  77.7560256304,  80.8975558711,
     $   84.0390907769,  87.1806298436,  90.3221726372,  93.4637187819,
     $   96.6052679510,  99.7468198587, 102.8883742542, 106.0299309165,
     $  109.1714896498, 112.3130502805, 115.4546126537, 118.5961766309,
     $  121.7377420880, 124.8793089132, 128.0208770059, 131.1624462752,
     $  134.3040166383, 137.4455880203, 140.5871603528, 143.7287335737,
     $  146.8703076258, 150.0118824570, 153.1534580192, 156.2950342685/
C
      pi = 4.*atan(one)
      nn = n
      if (n.gt.50) then
        bes(50) = bz(50)
        do j=51,n
          bes(j) = bes(j-1) + pi
        end do
        nn = 49
      end if
      do j=1,nn
        bes(j) = bz(j)
      end do
      return
      end
