#include <params.h>
      subroutine relunit(iunit)
C-----------------------------------------------------------------------
C
C Utility function to return a unit number no longer in use
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C                    T. Acker, March 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C--------------------------Parameters-----------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C------------------------------Commons----------------------------------
#include <comlun.h>
C-----------------------------Arguments---------------------------------
C
C Input arguments
C
      integer iunit    ! unit number
C
C-----------------------------------------------------------------------
C
C Print error messages when appropriate
C
      if (.not.lutag(iunit)) then
        write(6,*)'Warning: Unit ',iunit,' is not flagged as in use'
      end if
      if (iunit.lt.1 .or. iunit.gt.99) then
        write(6,*)'RELUNIT: Attempt to return out of range unit=',iunit
        call endrun
      end if
C
C Close file associated with unit number
C
      close(iunit)
C
C Release unit number (return to pool of unused numbers)
C
      lutag(iunit) = .false.
C
      return
C
      end
 
