module m_zeitTimer
! This module is an implementation of abstractTimer extension

  use m_abstractTimer, only: anytimer => abstractTimer
  use m_zeit, only: PUTIME,PSTIME,MWTIME
  implicit none
  private
  public:: timer
  public:: timer_typemold

  type, extends(anytimer):: timer
    private
  contains
    procedure:: on
    procedure:: off
    procedure:: reset
    procedure:: flush
    procedure:: allflush
    procedure, nopass:: mytype => typename
  end type timer

  character(len=*),parameter:: myname="m_zeitTimer"
  logical,parameter:: verbose=.false.
  !logical,parameter:: verbose=.true.

  type(timer),target:: mold_
  type(timer),parameter:: const_ = timer()

        ! UMTIME represents "user-mask", a customized timer mask.

  !integer,parameter:: UWTIME = PUTIME+PSTIME+MWTIME
  integer,parameter:: UWTIME = MWTIME

        ! where, individual timer masks are for
        !       PUTIME: process-user-time (~ CPU time)
        !       PSTIME: process-system-time (~ I/O time)
        !       MWTIME: MPI wall-clock time

contains

function timer_typemold()
  implicit none
  type(timer),pointer:: timer_typemold
  timer_typemold => mold_
end function timer_typemold
!-------------------------------------------------------
function typename()
  implicit none
  character(len=:),allocatable:: typename
  typename='['//myname//'::timer]'
end function typename

subroutine on(tm,name)
  use m_zeit, only: zeit_ci
  implicit none
  class(timer), intent(inout):: tm
  character(len=*),intent(in):: name
  call zeit_ci(name)
end subroutine on

subroutine off(tm,name)
  use m_zeit, only: zeit_co
  implicit none
  class(timer), intent(inout):: tm
  character(len=*),intent(in):: name
  call zeit_co(name)
end subroutine off

subroutine reset(tm)
  use m_zeit, only: zeit_reset
  implicit none
  class(timer), intent(inout):: tm
  call zeit_reset()
end subroutine reset

subroutine flush(tm,lu)
  use m_zeit, only: zeit_flush
  implicit none
  class(timer), intent(in):: tm
  integer     , intent(in):: lu

  call zeit_flush(lu,umask=UWTIME,subname_at_end=.true.)
end subroutine flush

subroutine allflush(tm,lu,comm,root)
  use m_zeit, only: zeit_allflush 
  implicit none
  class(timer), intent(in):: tm
  integer     , intent(in):: lu
  integer     , intent(in):: comm
  integer     , intent(in):: root

  call zeit_allflush(comm,root,lu,umask=UWTIME,subname_at_end=.true.)
end subroutine allflush

end module m_zeitTimer
