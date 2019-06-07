#include <params.h>
      integer function getunit(iunit)
C-----------------------------------------------------------------------
C
C Return index of next available Fortran unit number.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      L. Bath, Jun 1992
C                    L. Buja, Feb 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Commons----------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C-----------------------------------------------------------------------
#include <comlun.h>
C---------------------------Local variables-----------------------------
C
      integer iunit
C
C-----------------------------------------------------------------------
C
      if (.not.lutag(iunit)) then
         lutag(iunit) = .true.
         getunit = iunit
         return
      end if
      write(6,*) 'GETUNIT: Unit ',iunit,' is in use.  Pick another one'

      call endrun
      end
 
