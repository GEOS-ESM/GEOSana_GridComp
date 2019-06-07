      subroutine blowup
!-----------------------------------------------------------------------
!
! Model is blowing up or something went wrong reading the boundary
! datasets.  Dispose current history tape and quit.  This routine was
! cannibalized from WRAPUP.  It was not possible to simply include the
! history dispose logic in ENDRUN because of potential recursion
! problems.
!
!     write(6,*) 'Blowup!!!!'
      stop
      return
      end
