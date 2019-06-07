#include <params.h>
      subroutine mtdlss(a       ,b       ,c       ,rhs     ,x       ,
     $                  n       ,ld      ,num     ,indx    ,ws      ,
     $                  lws     ,ier     )
C------------------------------------------------------------------------
C
C     solves multiple tridiagonal linear systems ax=b for the vector x
C
C     the input vectors are:
C        a  : subdiagonal elements   -
C        b  : diagonal elements       |- from matrix a above (see below)
C        c  : superdiagonal elements -
C        rhs: right hand side (vector b above)
C        indx: Index array for gather
C        npts: Number of points in indx
C
C     other inputs:
C        n  : order of linear system (length of diagonal element vector)
C        ld : leading dimension of subdiagonal matrices (see below)
C        num: number of tridiagonal systems to be solved
C        ws : workspace declared in calling routine
C             must be at least of length (ld*n)
C        lws: length of workspace ws as declared in calling program.
C
C     outputs:
C        x  : solution vector x
C        ier: error status flag
C             .eq. 0 => no errors
C             .lt. 0 => insufficient storage declared for array ws
C                       or num .gt. ld (both declared storage errors)
C             .gt. 0 => numerical error (singular matrix, etc.)
C                       value is equal to location of linear system
C                       that failed where 1 .le. ier .le. num
C
C     nb: the parameter ld is the leading dimension of all input and
C         output arrays, while n, the order of the linear systems, is
C         the second dimension.  this storage strategy is adopted for
C         efficient memory (i.e., stride 1) referencing when solving
C         multiple systems on a vector machine.  thus, num, the number
C         of systems to be solved, can take any value between 1 and ld.
C         Care should be taken defining a and c in the calling program.
C         Although the diagonal information a, b, and c, all are of
C         length n, the subdiagonal coefficients, a, are only defined
C         for n = 2, 3, 4, ... n, while the superdiagonal coefficients,
C         c, are only defined for n = 1, 2, 3, ... n-1 (i.e., the index
C         refers to the "row number" of the diagonal element and as such
C         the first element of a and last element of c do not exist).
C
C         since there is no pivoting, it is possible, but unlikely that
C         this routine could fail due to numerical instability, even if
C         the matrix a is nonsingular.  there is no attempt to diagnose
C         the nature of a failure if one is encountered in the procedure
C
C----------------------------Code History-------------------------------
C
C Original version:  J. Hack
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, August 1992
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
#include <implicit.h>
C------------------------------Arguments--------------------------------
C
      integer n  
      integer ld  
      integer num  
      integer lws 
      integer indx(ld)
C
      real a(ld,n)  
      real b(ld,n)  
      real c(ld,n)  
      real rhs(ld,n)  
      real x(ld,n)  
      real ws(ld,n) 
C
      integer ier
C
C---------------------------Local variables-----------------------------
C
      integer i,ii
      integer j
C
C-----------------------------------------------------------------------
C
      ier = 0
C
C Check for sufficient working storage in workspace ws and for
C logically consistant storage declaration for input arrays
C
      if ((lws .lt. (ld*n)) .or. (num.gt.ld)) then
        ier = -1
        return
      end if
C
C Decomposition and forward substitution loops
C Note: references to ws(i,1) denote special use of workspace ws
C
CDIR$ IVDEP
      do ii=1,num
        i = indx(ii)
        if (b(i,1) .eq. 0.0) then
          ier = i
          return
        end if
      end do
C
      do ii=1,num
        i = indx(ii)
        ws(i,1) = b(i,1)
        x (i,1) = rhs(i,1)/ws(i,1)
      end do
C
      do j=2,n
CDIR$ IVDEP
        do ii=1,num
          i = indx(ii)
          ws(i,j) = c(i,j-1)/ws(i,1)
          ws(i,1) = b(i,j) - a(i,j)*ws(i,j)
        end do
CDIR$ IVDEP
        do ii=1,num
          i = indx(ii)
          if (ws(i,1) .eq. 0.0) then
            ier = i
            return
          end if
        end do
CDIR$ IVDEP
        do ii=1,num
          i = indx(ii)
          x(i,j) = (rhs(i,j) - a(i,j)*x(i,j-1))/ws(i,1)
        end do
      end do
C
C Backsubstitution loop
C
      do j=n-1,1,-1
CDIR$ IVDEP
        do ii=1,num
          i = indx(ii)
          x(i,j) = x(i,j) - ws(i,j+1)*x(i,j+1)
        end do
      end do
C
      return
      end
 
