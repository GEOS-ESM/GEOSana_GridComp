!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ozTmpl - ozinfo template file IO
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ozTmpl
      implicit none
      private	! except

      public :: ozTmpl		! data type
      public :: ozTmpl_get	! get from a file
      public :: ozTmpl_put	! put to a file

      public :: ptr_nusis	! reference to tmpl(:)%nusis
      public :: ptr_nulev	! reference to tmpl(:)%nulev
      public :: ptr_iuse	! reference to tmpl(:)%iuse

      public :: ozTmpl_alloc	! the initial allocation
      public :: ozTmpl_realloc	! an incremental allocation, if needed
      public :: ozTmpl_dealloc	! simple deallocate

      public :: NAMESIZE

    interface ozTmpl_get; module procedure get_; end interface
    interface ozTmpl_put; module procedure put_; end interface

    interface ozTmpl_alloc  ; module procedure   alloc_; end interface
    interface ozTmpl_realloc; module procedure realloc_; end interface
    interface ozTmpl_dealloc; module procedure dealloc_; end interface

    integer,parameter:: PAGESIZE=64
    integer,parameter:: NAMESIZE=20

    type ozTmpl
      private
      character(len=NAMESIZE):: nusis
      integer                 :: nulev
      integer                 :: iuse
      real   , dimension(5)   :: params
    end type ozTmpl


! !REVISION HISTORY:
! 	08Jun07	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ozTmpl'

#include "assert.H"
contains
function ptr_nusis(tmpl,lbnd,ubnd,strd) result(ptr)
  implicit none
  type(ozTmpl),target,dimension(:),intent(in) :: tmpl    ! the target object
  integer,optional,intent(in):: lbnd,ubnd,strd
  character(len=NAMESIZE),pointer,dimension(:):: ptr	! returned pointer 
  integer:: lbnd_,ubnd_
  lbnd_=lbound(tmpl,1); if(present(lbnd)) lbnd_=lbnd
  ubnd_=ubound(tmpl,1); if(present(ubnd)) ubnd_=ubnd
  if(present(strd)) then
    ptr => tmpl(lbnd_:ubnd_:strd)%nusis
  else
    ptr => tmpl(lbnd_:ubnd_     )%nusis
  endif
end function ptr_nusis
function ptr_nulev(tmpl,lbnd,ubnd,strd) result(ptr)
  implicit none
  type(ozTmpl),target,dimension(:),intent(in) :: tmpl    ! assumed null()
  integer,optional,intent(in):: lbnd,ubnd,strd
  integer,pointer,dimension(:):: ptr
  integer:: lbnd_,ubnd_
  lbnd_=lbound(tmpl,1); if(present(lbnd)) lbnd_=lbnd
  ubnd_=ubound(tmpl,1); if(present(ubnd)) ubnd_=ubnd
  if(present(strd)) then
    ptr => tmpl(lbnd_:ubnd_:strd)%nulev
  else
    ptr => tmpl(lbnd_:ubnd_     )%nulev
  endif
end function ptr_nulev
function ptr_iuse(tmpl,lbnd,ubnd,strd) result(ptr)
  implicit none
  type(ozTmpl),target,dimension(:),intent(in) :: tmpl    ! assumed null()
  integer,optional,intent(in):: lbnd,ubnd,strd
  integer,pointer,dimension(:) :: ptr
  integer:: lbnd_,ubnd_
  lbnd_=lbound(tmpl,1); if(present(lbnd)) lbnd_=lbnd
  ubnd_=ubound(tmpl,1); if(present(ubnd)) ubnd_=ubnd
  if(present(strd)) then
    ptr => tmpl(lbnd_:ubnd_:strd)%iuse
  else
    ptr => tmpl(lbnd_:ubnd_     )%iuse
  endif
end function ptr_iuse

subroutine alloc_(ptr,nptr,nominal,page)
    ! allocate an INTEGER pointer(:), which is assumed null()
  implicit none
  type(ozTmpl),pointer,dimension(:) :: ptr    ! assumed null()
  integer,intent(out) :: nptr            ! content size, set to 0
  integer,optional,intent(in) :: nominal ! nominal minimum size
  integer,optional,intent(in) :: page    ! block size of this alloc.

  integer :: psz,msz

  psz=PAGESIZE
  if(present(page)) psz=max(page,PAGESIZE)
  msz=psz
  if(present(nominal)) msz=max(1,(nominal+psz-1)/psz)*psz
  allocate(ptr(msz))
  nptr=0
end subroutine alloc_

subroutine realloc_(ptr,nptr,incr,page)
    ! reallocate an INTEGER pointer(:), which is associated().
  use satinfo_util, only: assert_
  implicit none
  type(ozTmpl),pointer,dimension(:):: ptr	! is assumed associated()
  integer,             intent(inout):: nptr	! record count, present but unchanged
  integer,             intent(in   ):: incr	! an estimated size increment
  integer, optional,   intent(in   ):: page	! block size of this increase.

  type(ozTmpl),pointer,dimension(:):: tmp
  integer :: msz,psz

  msz=size(ptr)
    ASSERT(msz>=nptr)
  if( nptr+incr <= msz) return

  psz=PAGESIZE
  if(present(page)) psz=page

  msz=msz+max(0,((incr+psz-1)/psz)*psz)
  tmp => ptr
  allocate(ptr(msz))
  ptr(1:nptr)=tmp(1:nptr)
  deallocate(tmp)
end subroutine realloc_

subroutine dealloc_(ptr,nptr)
    ! deallocate an INTEGER pointer(:).
  implicit none
  type(ozTmpl),pointer,dimension(:):: ptr	! is assumed associated()
  integer,intent(out) :: nptr         		! content size, set to -1
  deallocate(ptr)
  nptr=-1
end subroutine dealloc_

subroutine get_(fname,jpch,tmpl,append)
!--  Read in info template file here
  use satinfo_util,only : luavail,stdout
  use satinfo_util,only : die,perr,tell
  implicit none
  character(len=*),    intent(in   ):: fname	! input info.rc.tmpl
  integer         ,    intent(inout):: jpch	! a count of entries
  type(ozTmpl),pointer,dimension(:):: tmpl	! a table of entries
  logical,optional,    intent(in   ):: append
  
  character(len=*),parameter :: myname_=myname//"::get_"

  integer :: lu,ios,j,ir
  character(len=512) :: line
  character(len=1) :: key
  logical :: append_,skipline
  append_=.false.; if(present(append)) append_=append

!! Open ...
  lu=luavail()
  open(lu,file=fname,form='formatted',status='old',iostat=ios)
    if (ios /= 0) call die(myname_, &
      'open("'//trim(fname)//'") for input, iostat =',ios)

  call tell(myname_,'reading "'//trim(fname)//'" for input')

  if(.not.append_) call alloc_(tmpl ,jpch) 	!! note jpch is set to 0 by alloc()
  j=jpch

	ir=0
	read(lu,'(a)',iostat=ios) line
	do while(ios==0)
	  ir=ir+1

          call realloc_(tmpl,j,incr=1)

	  tmpl(j+1) = ozTmpl('.undef.',HUGE(j),HUGE(j),HUGE(1.))

	  key='?'
	  read(line,*,iostat=ios) key	! use "read() key" instead of "key=..".
	  select case(key(1:1))
	  case('?')
	    call die(myname_,'cann''t read, rec # =',ir)

	  case('!')
	    skipline=.true.

	  case default		! in the "new" format
	    skipline=.false.
	    j=j+1
            read(line,*,iostat=ios) tmpl(j)
		if(ios/=0) exit		! It has reached a blank record.

	  	! If the input line contains unexpected values,
		! which is often the result of a failed "list-
		! directed" read (i.e. fmt=*), echo the input
		! then die().

	    if(     tmpl(j)%nusis == '.undef.' .or. &
	    	    tmpl(j)%nulev == HUGE(j)   .or. &
	    	    tmpl(j)%iuse  == HUGE(j)   .or. &
	    	any(tmpl(j)%params== HUGE(1.)) ) then

	      call perr(myname_,'tmpl_j =',j)
	      call perr(myname_,'%nusis =',tmpl(j)%nusis )
	      call perr(myname_,'%nulev =',tmpl(j)%nulev )
	      call perr(myname_,'%iuse  =',tmpl(j)%iuse  )
	      call perr(myname_,'%params(1)=',tmpl(j)%params(1))
	      call perr(myname_,'%params(2)=',tmpl(j)%params(2))
	      call perr(myname_,'%params(3)=',tmpl(j)%params(3))
	      call perr(myname_,'%params(4)=',tmpl(j)%params(4))
	      call perr(myname_,'%params(5)=',tmpl(j)%params(5))
	      call perr(myname_,'"'//trim(line)//'"')
	      call die(myname_,'unexpected info entry at rec # =',j)
	    endif

	  	! Count in the record as a good table entry, and
		! set the use_rad flag to _off_.
	    jpch=j
	  end select

		! Read the next record
	  read(lu,'(a)',iostat=ios) line
        enddo
        close(lu)

	call tell(myname_,'number of total nusis-by-levels, jpch =',jpch)

    if(jpch==0) call die(myname_, &
      'no coefficient found in "'//trim(fname)//'", jpch =',jpch)
end subroutine get_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: put_ - write out satinfo table for GSI
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine put_(fname,jpch,tmpl,nymd,nhms)
      use satinfo_util,only : luavail
      use satinfo_util,only : die,tell
      implicit none
      character(len=*),intent(in) :: fname      ! output xxinfo.txt
      integer,intent(in) :: jpch		! a count of entries
      type(ozTmpl),pointer,dimension(:):: tmpl	! a table of entries
      integer,intent(in) :: nymd,nhms

! !REVISION HISTORY:
! 	08Jun07	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::put_'
  integer :: lu,ios,j,l
  character(len=max(14,len(tmpl%nusis))) :: nusisj

        ! Modify original satinfo file according to the satinfo_base.rc
        !===============================================================
        
  call tell(myname_,'writing "'//trim(fname)//'" for output')

	lu=luavail()
        open(lu,file=fname,form='formatted',status='unknown',iostat=ios)
	  if(ios/=0) call die(myname_, &
	    'open("'//trim(fname)//'") for output, iostat =',ios)

	! Write header lines

	write(lu,'(a)') "! For mls data, pressure and obs errors are pulled from bufr, so not listed here"
	write(lu,'(a)') "! sens/instr/sat lev  use pressure gross   obs    b_oz  pg_oz"
	write(lu,'(a)') "!                                  error  error variational qc"

        do j = 1,jpch
	  nusisj=tmpl(j)%nusis
	  l=max(len_trim(nusisj),14)
	  write(lu,'(1x,a,i4,i5,1x,f8.3,4f7.3)',iostat=ios) nusisj(1:l), &
	  	tmpl(j)%nulev, tmpl(j)%iuse, tmpl(j)%params(:)

	  if(ios/=0) call die(myname_, &
	      'write("'//trim(fname)//'") for output, iostat =',ios)
	enddo
end subroutine put_
end module m_ozTmpl
