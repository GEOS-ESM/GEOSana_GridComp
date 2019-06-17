!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_convTmpl - convinfo template file IO
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_convTmpl
      implicit none
      private	! except

      public :: convTmpl	! data type
      public :: convTmpl_get	! get from a file
      public :: convTmpl_put	! put to a file

      public :: ptr_otype	! reference to tmpl(:)%otype, obs. type
      				!   e.g. ps, pw, uv, pm2_5, etc.
      public :: ptr_type	! reference to tmpl(:)%type, platform type
      public :: ptr_sub		! reference to tmpl(:)%sub, instrument class
      public :: ptr_iuse	! reference to tmpl(:)%iuse

      public :: convTmpl_alloc	! the initial allocation
      public :: convTmpl_realloc	! an incremental allocation, if needed
      public :: convTmpl_dealloc	! simple deallocate

      public :: NAMESIZE

    interface convTmpl_get; module procedure get_; end interface
    interface convTmpl_put; module procedure put_; end interface

    interface convTmpl_alloc  ; module procedure   alloc_; end interface
    interface convTmpl_realloc; module procedure realloc_; end interface
    interface convTmpl_dealloc; module procedure dealloc_; end interface

    integer,parameter:: PAGESIZE=64
    integer,parameter:: NAMESIZE=10
    integer,parameter:: BUFRSIZE=256

    type convTmpl
      private
      character(len=NAMESIZE) :: otype
      integer                 :: type
      integer                 :: sub
      integer                 :: iuse
      character(len=BUFRSIZE) :: rest
    end type convTmpl

! !REVISION HISTORY:
! 	08Jun07	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_convTmpl'

#include "assert.H"
contains
function ptr_otype(tmpl,lbnd,ubnd,strd) result(ptr)
  implicit none
  type(convTmpl),target,dimension(:),intent(in) :: tmpl    ! the target object
  integer,optional,intent(in):: lbnd,ubnd,strd
  character(len=NAMESIZE),pointer,dimension(:):: ptr	! returned pointer 
  integer:: lbnd_,ubnd_
  lbnd_=lbound(tmpl,1); if(present(lbnd)) lbnd_=lbnd
  ubnd_=ubound(tmpl,1); if(present(ubnd)) ubnd_=ubnd
  if(present(strd)) then
    ptr => tmpl(lbnd_:ubnd_:strd)%otype
  else
    ptr => tmpl(lbnd_:ubnd_     )%otype
  endif
end function ptr_otype
function ptr_type(tmpl,lbnd,ubnd,strd) result(ptr)
  implicit none
  type(convTmpl),target,dimension(:),intent(in) :: tmpl    ! assumed null()
  integer,optional,intent(in):: lbnd,ubnd,strd
  integer,pointer,dimension(:):: ptr
  integer:: lbnd_,ubnd_
  lbnd_=lbound(tmpl,1); if(present(lbnd)) lbnd_=lbnd
  ubnd_=ubound(tmpl,1); if(present(ubnd)) ubnd_=ubnd
  if(present(strd)) then
    ptr => tmpl(lbnd_:ubnd_:strd)%type
  else
    ptr => tmpl(lbnd_:ubnd_     )%type
  endif
end function ptr_type
function ptr_sub(tmpl,lbnd,ubnd,strd) result(ptr)
  implicit none
  type(convTmpl),target,dimension(:),intent(in) :: tmpl    ! assumed null()
  integer,optional,intent(in):: lbnd,ubnd,strd
  integer,pointer,dimension(:) :: ptr
  integer:: lbnd_,ubnd_
  lbnd_=lbound(tmpl,1); if(present(lbnd)) lbnd_=lbnd
  ubnd_=ubound(tmpl,1); if(present(ubnd)) ubnd_=ubnd
  if(present(strd)) then
    ptr => tmpl(lbnd_:ubnd_:strd)%sub
  else
    ptr => tmpl(lbnd_:ubnd_     )%sub
  endif
end function ptr_sub
function ptr_iuse(tmpl,lbnd,ubnd,strd) result(ptr)
  implicit none
  type(convTmpl),target,dimension(:),intent(in) :: tmpl    ! assumed null()
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
  type(convTmpl),pointer,dimension(:) :: ptr    ! assumed null()
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
  type(convTmpl),pointer,dimension(:):: ptr	! is assumed associated()
  integer,             intent(inout):: nptr	! record count, present but unchanged
  integer,             intent(in   ):: incr	! an estimated size increment
  integer, optional,   intent(in   ):: page	! block size of this increase.

  type(convTmpl),pointer,dimension(:):: tmp
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
  type(convTmpl),pointer,dimension(:):: ptr	! is assumed associated()
  integer,intent(out) :: nptr         		! content size, set to -1
  deallocate(ptr)
  nptr=-1
end subroutine dealloc_

function iseq_(a,b)
  implicit none
  type(convTmpl),intent(in):: a,b
  logical:: iseq_
  iseq_ = a%otype==b%otype .and. &
  	  a% type==b% type .and. &
	  a%  sub==b%  sub .and. &
	  a% iuse==b% iuse .and. &
	  a% rest==b% rest
end function iseq_

subroutine get_(fname,jpch,tmpl,append)
!--  Read in info template file here
  use satinfo_util,only : luavail,stdout
  use satinfo_util,only : die,perr,tell
  implicit none
  character(len=*),    intent(in   ):: fname	! input info.rc.tmpl
  integer         ,    intent(inout):: jpch	! a count of entries
  type(convTmpl),pointer,dimension(:):: tmpl	! a table of entries
  logical,optional,    intent(in   ):: append
  
  character(len=*),parameter :: myname_=myname//"::get_"

  type(convTmpl),parameter:: udeftmpl=convTmpl('.undef.',HUGE(1),HUGE(1),HUGE(1),'')
  integer :: lu,ios,j,ir
  character(len=BUFRSIZE) :: line
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

	  tmpl(j+1) = udeftmpl

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
	    tmpl(j)%rest=line

	  	! If the input line contains unexpected values,
		! which is often the result of a failed "list-
		! directed" read (i.e. fmt=*), echo the input
		! then die().

	    if(iseq_(tmpl(j),udeftmpl)) then
	      call perr(myname_,'%otype =',tmpl(j)%otype )
	      call perr(myname_,'%type  =',tmpl(j)%type  )
	      call perr(myname_,'%sub   =',tmpl(j)%sub   )
	      call perr(myname_,'%iuse  =',tmpl(j)%iuse  )
	      call perr(myname_,'%rest  =',tmpl(j)%rest  )
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
      type(convTmpl),pointer,dimension(:):: tmpl	! a table of entries
      integer,intent(in) :: nymd,nhms

! !REVISION HISTORY:
! 	08Jun07	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::put_'
  integer :: lu,ios,j,l
  type(convTmpl):: tmpl_j

        ! Modify original satinfo file according to the satinfo_base.rc
        !===============================================================
        
  call tell(myname_,'writing "'//trim(fname)//'" for output')

	lu=luavail()
        open(lu,file=fname,form='formatted',status='unknown',iostat=ios)
	  if(ios/=0) call die(myname_, &
	    'open("'//trim(fname)//'") for output, iostat =',ios)

	! Write header lines

  write(lu,'(a,i8.8,a,i6.6,a)') "! -- created by make_convinfo, for ",nymd,":",nhms," --"
  write(lu,'(a)') "! otype   = observation type (a7, t, uv, q, etc.)"
  write(lu,'(a)') "! type    = prepbufr observation type (if available)"
  write(lu,'(a)') "! sub     = prepbufr subtype (not yet available)"
  write(lu,'(a)') "! iuse    = flag if to use/not use / monitor data"
  write(lu,'(a)') "!         = 1  use data"
  write(lu,'(a)') "!         = 0  do not use data"
  write(lu,'(a)') "!         = -1 monitor data"
  write(lu,'(a)') "! twindow = time window (+/- hours)"
  write(lu,'(a)') "! numgrp  = cross validation parameter - number of groups"
  write(lu,'(a)') "! ngroup  = cross validation parameter - group to remove from data use"
  write(lu,'(a)') "! nmiter  = cross validation parameter - external iteration to introduce removed data"
  write(lu,'(a)') "! gross   = gross error parameter - gross error"
  write(lu,'(a)') "! ermax   = gross error parameter - max error"
  write(lu,'(a)') "! ermin   = gross error parameter - min error"
  write(lu,'(a)') "! var_b   = variational quality control parameter -  b parameter"
  write(lu,'(a)') "! var_pg ithin rmesh npred  = variational quality control parameter -  pg parameter"
  write(lu,'(a)') "!otype   type  sub iuse twindow numgrp ngroup nmiter gross ermax ermin var_b    var_pg ithin rmesh  pmesh  npred"

        do j = 1,jpch
	  l=max(len_trim(tmpl(j)%otype),8)
	  write(lu,'(1x,a,i4.3,2i5,5x,a)',iostat=ios) tmpl(j)%otype(1:l), &
	  	tmpl(j)%type, tmpl(j)%sub, tmpl(j)%iuse, trim(tmpl(j)%rest(29:))

	  if(ios/=0) call die(myname_, &
	      'write("'//trim(fname)//'") for output, iostat =',ios)
	enddo
end subroutine put_
end module m_convTmpl
