#include <params.h>
      integer function navu()
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
#include <pagrid.h>
#include <comlun.h>
C---------------------------Local variables-----------------------------
C
      integer itst      ! Index
C
C-----------------------------------------------------------------------
C
      do itst=1,99
        if(.not.lutag(itst)) then
          navu = itst
          lutag(itst) = .true.
          return
        end if
      end do
      write(6,*) 'NAVU: Ran out of Fortran unit numbers'
      call endrun       ! Abnormal termination
      end
 
