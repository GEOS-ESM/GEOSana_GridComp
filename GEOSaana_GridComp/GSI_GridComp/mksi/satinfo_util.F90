!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: satinfo_util - utilities
!
! !DESCRIPTION:
!   to be compiled with -Dsys`uname -s`
!
! !INTERFACE:

    module satinfo_util
      implicit none
      private	! except

      public :: die, perr, warn, tell, assert_
      public :: luavail
      public :: stdin, stdout, stderr
      public :: alloc, realloc, dealloc
      public :: getarec
      public :: locate

    integer,parameter :: STDIN = 5
    integer,parameter :: STDOUT= 6
#ifdef sysHP_UX
    integer,parameter :: STDERR= 7
#else
    integer,parameter :: STDERR= 0
#endif
    interface luavail; module procedure luavail_; end interface

    interface die; module procedure &
      die_bul_, &
      die_chr_, &
      die_int_, &
      die_vint_, &
      die_flt_, &
      die_dbl_, &
      die2_,	&
      die_; end interface
    interface perr; module procedure &
      perr_bul_, &
      perr_chr_, &
      perr_int_, &
      perr_vint_, &
      perr_flt_, &
      perr_dbl_, &
      perr_; end interface
    interface warn; module procedure &
      warn_bul_, &
      warn_chr_, &
      warn_int_, &
      warn_vint_, &
      warn_flt_, &
      warn_dbl_, &
      warn_; end interface
    interface tell; module procedure &
      tell_bul_, &
      tell_chr_, &
      tell_int_, &
      tell_vint_, &
      tell_flt_, &
      tell_dbl_, &
      tell_; end interface

    interface mprint; module procedure &
      mprint_bul_, &
      mprint_chr_, &
      mprint_int_, &
      mprint_vint_, &
      mprint_flt_, &
      mprint_dbl_, &
      mprint_; end interface

    interface alloc; module procedure &
      flt_alloc_, &
      dbl_alloc_, &
      int_alloc_, &
      chr_alloc_; end interface
    interface realloc; module procedure &
      flt_realloc_, &
      dbl_realloc_, &
      int_realloc_, &
      chr_realloc_; end interface
    interface dealloc; module procedure &
      flt_dealloc_, &
      dbl_dealloc_, &
      int_dealloc_, &
      chr_dealloc_; end interface

    interface getarec; module procedure getarec_; end interface
    interface locate; module procedure &
      int_locate_, &
      chr_locate_; end interface

! !REVISION HISTORY:
! 	02May07	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname="satinfo_util"

  real*4,parameter :: R4=1.E+0
  real*8,parameter :: R8=1.D+0
  integer,parameter :: SP=kind(R4)
  integer,parameter :: DP=kind(R8)

  integer,parameter :: PAGESIZE_=64
  integer,parameter :: MAX_LUNIT=1024

#include "assert.H"

!!$ Usecases:
!!$
!!$ Usecase 1: store data into a reallocatable array.
!!$
!!$	use satinfo_util,only : alloc,realloc,dealloc
!!$	integer,pointer,dimensioni(:) :: iptr =>null()
!!$	integer :: ilen,itmp,ierr
!!$	call alloc(iptr,ilen)
!!$	read(*,*,iostat=ierr) itmp
!!$	do while(ierr==0)
!!$	  ![.. can iptr(:) store one more?  If not, realloc it ..]
!!$	  call realloc(iptr,ilen,incr=1)
!!$	  ilen=ilen+1
!!$	  iptr(ilen)=itmp
!!$	  read(*,*,iostat=ierr) itmp
!!$	end do
!!$	![.. array iptr(1:ilen) contains all values of itmp ..]
!!$	call dealloc(iptr,ilen)
!!$
!!$ Usecase 2: get a table record from an ascii file containing comments
!!$
!!$	use satinfo_util,only: getarec
!!$	character(len=MXRECLEN) :: arec
!!$	character(len=32) :: t
!!$	character(len=32),pointer,dimension(:) :: tags=>null()
!!$	call alloc(tags,ntag)	! initialize tags
!!$	call getarec(5,arec,ier)
!!$	do while(ier==0)
!!$	  read(arec,*,iostat=ier) t
!!$	    if(ier/=0) [...]
!!$	    call realloc(tags,ntag,incr=1)
!!$	    tags(ntag)=t
!!$	  call getarec(5,arec,ier)
!!$	enddo
!!$	[.. tags(1:ntag) is a list of tags ..]
!!$	call dealloc(tags,ntag)

contains
subroutine int_alloc_(ptr,nptr,nominal,page)
    ! allocate an INTEGER pointer(:), which is assumed null()
  implicit none
  integer,pointer,dimension(:) :: ptr    ! assumed null()
  integer,intent(out) :: nptr            ! content size, set to 0
  integer,optional,intent(in) :: nominal ! nominal minimum size
  integer,optional,intent(in) :: page    ! block size of this alloc.

  integer :: psz,msz

  psz=PAGESIZE_
  if(present(page)) psz=max(page,PAGESIZE_)
  msz=psz
  if(present(nominal)) msz=max(1,(nominal+psz-1)/psz)*psz
  allocate(ptr(msz))
  nptr=0
end subroutine int_alloc_

subroutine int_realloc_(ptr,nptr,incr,page)
    ! reallocate an INTEGER pointer(:), which is associated().
  implicit none
  integer,pointer,dimension(:) :: ptr ! assumed assoicated()
  integer,intent(in) :: nptr          ! content size, unchanged
  integer,intent(in) :: incr          ! required size increment
  integer,optional,intent(in) :: page ! block size of this increase.

  integer,pointer,dimension(:) :: tmp
  integer :: msz,psz

  msz=size(ptr)
    ASSERT(msz>=nptr)
  if( nptr+incr <= msz) return

  psz=PAGESIZE_
  if(present(page)) psz=page

  msz=msz+max(0,((incr+psz-1)/psz)*psz)
  tmp => ptr
  allocate(ptr(msz))
  ptr(1:nptr)=tmp(1:nptr)
  deallocate(tmp)
end subroutine int_realloc_

subroutine int_dealloc_(ptr,nptr)
    ! deallocate an INTEGER pointer(:).
  implicit none
  integer,pointer,dimension(:) :: ptr ! assumed associated()
  integer,intent(out) :: nptr         ! content size, set to -1
  deallocate(ptr)
  nptr=-1
end subroutine int_dealloc_

subroutine chr_alloc_(ptr,nptr,nominal,page)
    ! allocate a CHARACTER(LEN=*) pointer(:), assumed null().
  implicit none
  character(len=*),pointer,dimension(:) :: ptr ! assumed null()
  integer,intent(out) :: nptr                  ! content size, set to 0
  integer,optional,intent(in) :: nominal ! nominal minimum size
  integer,optional,intent(in) :: page    ! block size of this alloc.

  integer :: psz,msz

  psz=PAGESIZE_
  if(present(page)) psz=max(page,PAGESIZE_)
  msz=psz
  if(present(nominal)) msz=max(1,(nominal+psz-1)/psz)*psz
  allocate(ptr(msz))
  nptr=0
end subroutine chr_alloc_

subroutine chr_realloc_(ptr,nptr,incr,page)
    ! reallocate a CHARACTER(LEN=*) pointer(:), assumed associated().
  implicit none
  character(len=*),pointer,dimension(:) :: ptr ! assumed associated()
  integer,intent(in) :: nptr          ! content size, unchanged
  integer,intent(in) :: incr          ! required size increment
  integer,optional,intent(in) :: page ! block size of this increase.

  character(len=len(ptr)),pointer,dimension(:) :: tmp
  integer :: msz,psz

  msz=size(ptr)
    ASSERT(msz>=nptr)
  if( nptr+incr <= msz) return

  psz=PAGESIZE_
  if(present(page)) psz=page

  msz=msz+max(0,((incr+psz-1)/psz)*psz)
  tmp => ptr
  allocate(ptr(msz))
  ptr(1:nptr)=tmp(1:nptr)
  deallocate(tmp)
end subroutine chr_realloc_

subroutine chr_dealloc_(ptr,nptr)
    ! deallocate a CHARACTER(LEN=*) pointer(:)
  implicit none
  character(len=*),pointer,dimension(:) :: ptr ! assumed associated()
  integer,intent(out) :: nptr                  ! content size, set to -1
  deallocate(ptr)
  nptr=-1
end subroutine chr_dealloc_

subroutine flt_alloc_(ptr,nptr,nominal,page)
    ! allocate a REAL pointer(:), assumed null().
  implicit none
  real(SP),pointer,dimension(:) :: ptr       ! assumed null()
  integer,intent(out) :: nptr            ! content size, set to 0
  integer,optional,intent(in) :: nominal ! nominal minimum size
  integer,optional,intent(in) :: page    ! block size of this alloc.

  integer :: psz,msz

  psz=PAGESIZE_
  if(present(page)) psz=max(page,PAGESIZE_)
  msz=psz
  if(present(nominal)) msz=max(1,(nominal+psz-1)/psz)*psz
  allocate(ptr(msz))
  nptr=0
end subroutine flt_alloc_

subroutine flt_realloc_(ptr,nptr,incr,page)
    ! reallocate a REAL pointer(:), assumed associated().
  implicit none
  real(SP),pointer,dimension(:) :: ptr    ! assumed associated()
  integer,intent(in) :: nptr          ! content size, unchanged
  integer,intent(in) :: incr          ! required size increment
  integer,optional,intent(in) :: page ! block size of this increase.

  real(kind(ptr)),pointer,dimension(:) :: tmp
  integer :: msz,psz

  msz=size(ptr)
    ASSERT(msz>=nptr)
  if( nptr+incr <= msz) return

  psz=PAGESIZE_
  if(present(page)) psz=page

  msz=msz+max(0,((incr+psz-1)/psz)*psz)
  tmp => ptr
  allocate(ptr(msz))
  ptr(1:nptr)=tmp(1:nptr)
  deallocate(tmp)
end subroutine flt_realloc_

subroutine flt_dealloc_(ptr,nptr)
    ! deallocate a REAL pointer(:).
  implicit none
  real(SP),pointer,dimension(:) :: ptr    ! assumed associated()
  integer,intent(out) :: nptr         ! content size, set to -1
  deallocate(ptr)
  nptr=-1
end subroutine flt_dealloc_

subroutine dbl_alloc_(ptr,nptr,nominal,page)
    ! allocate a DOUBLE PRECISION pointer(:), assumed null().
  implicit none
  real(DP),pointer,dimension(:) :: ptr ! assumed null()
  integer,intent(out) :: nptr            ! content size, set to 0
  integer,optional,intent(in) :: nominal ! nominal minimum size
  integer,optional,intent(in) :: page    ! block size of this alloc.

  integer :: psz,msz

  psz=PAGESIZE_
  if(present(page)) psz=max(page,PAGESIZE_)
  msz=psz
  if(present(nominal)) msz=max(1,(nominal+psz-1)/psz)*psz
  allocate(ptr(msz))
  nptr=0
end subroutine dbl_alloc_

subroutine dbl_realloc_(ptr,nptr,incr,page)
    ! reallocate a DOUBLE PRECISION pointer(:), assumed assoicated().
  implicit none
  real(DP),pointer,dimension(:) :: ptr ! assumed associated()
  integer,intent(in) :: nptr          ! content size, unchanged
  integer,intent(in) :: incr          ! required size increment
  integer,optional,intent(in) :: page ! block size of this increase.

  real(kind(ptr)),pointer,dimension(:) :: tmp
  integer :: msz,psz

  msz=size(ptr)
    ASSERT(msz>=nptr)
  if( nptr+incr <= msz) return

  psz=PAGESIZE_
  if(present(page)) psz=page

  msz=msz+max(0,((incr+psz-1)/psz)*psz)
  tmp => ptr
  allocate(ptr(msz))
  ptr(1:nptr)=tmp(1:nptr)
  deallocate(tmp)
end subroutine dbl_realloc_

subroutine dbl_dealloc_(ptr,nptr)
    ! deallocate a DOUBLE PRECISION  pointer(:).
  implicit none
  real(DP),pointer,dimension(:) :: ptr ! assumed assoicated().
  integer,intent(out) :: nptr         ! content size, set to -1
  deallocate(ptr)
  nptr=-1
end subroutine dbl_dealloc_

subroutine getarec_(lu,line,ier,nrec,commchar)
  implicit none
  integer,intent(in) :: lu
  character(len=*),intent(out) :: line
  integer,intent(out) :: ier
  integer,optional,intent(out) :: nrec ! count of record readings
  character(len=*),optional,intent(in) :: commchar ! set of comment chars

character(len=1) :: c
character(len=*),parameter :: SPC=achar(32),TAB=achar(09)
character(len=*),parameter :: NUL=achar(00),COM='#'

if(present(nrec)) nrec=0

    ! Read records until a line of non-blank and any non-comment text.
    ! A pure commont text record is a record with non-block content
    ! lead by a "#".
  read(lu,'(a)',iostat=ier) line
  do while(ier==0)
    if(present(nrec)) nrec=nrec+1
    c=leadchar_(line)
    if(present(commchar)) then
      if(c/=SPC .and. c/=TAB .and. index(commchar,c)/=1) exit
    else
      if(c/=SPC .and. c/=TAB .and. c/=COM) exit
    endif
    read(lu,'(a)',iostat=ier) line
  enddo
contains
function leadchar_(line) result(c)
  implicit none
  character(len=*),intent(in) :: line
  character(len=1) :: c
  integer :: i,l
  i=0
  l=len(line)
  c=SPC
  do while(i<l)
    i=i+1
    if(line(i:i)==SPC .or. line(i:i)==TAB) cycle
    c=line(i:i)
    return
  enddo
end function leadchar_
end subroutine getarec_

function int_locate_(key,list) result(l)
  implicit none
  integer,intent(in) :: key
  integer,dimension(:),intent(in) :: list
  integer :: l
  integer :: i
  l=0
  do i=1,size(list)
    if(key/=list(i)) cycle
    l=i
    return
  enddo
end function int_locate_
function chr_locate_(key,list) result(l)
  implicit none
  character(len=*),intent(in) :: key
  character(len=*),dimension(:),intent(in) :: list
  integer :: l
  integer :: i
  l=0
  do i=1,size(list)
    if(key/=list(i)) cycle
    l=i
    return
  enddo
end function chr_locate_

function myid_(who)
  implicit none
  character(len=*),intent(in) :: who
  character(len=len(who)) :: myid_
  myid_=adjustl(who)
end function myid_

function luavail_() result(lu)
  implicit none
  integer :: lu

  character(len=*),parameter :: myname_=myname//'::luavail_'
  integer ios
  logical inuse
  character*8 attr

  lu=-1
  ios=0
  inuse=.true.

  do while(ios.eq.0.and.inuse)
    lu=lu+1

	! Test #1, reserved units

    inuse = lu.eq.stdout .or. lu.eq.stdin .or. lu.eq.stderr

#ifdef sysSunOS
	! Reserved units under SunOS
    inuse = lu.eq.100 .or. lu.eq.101 .or. lu.eq.102
#endif

	! Test #2, in-use

    if(.not.inuse) inquire(unit=lu,opened=inuse,iostat=ios)

    if(lu >= MAX_LUNIT) ios=-1
  end do
  if(ios.ne.0) lu=-1
end function luavail_

subroutine mprint_(lu,who,what)
  implicit none
  integer,intent(in) :: lu
  character(len=*),intent(in) :: who,what
  write(lu,'(3a)') &
    trim(myid_(who)),'(): ',trim(what)
end subroutine mprint_
subroutine mprint_chr_(lu,who,what,val)
  implicit none
  integer,intent(in) :: lu
  character(len=*),intent(in) :: who,what,val
  write(lu,'(3a)',advance='no') &
    trim(myid_(who)),'(): ',what
  write(lu,'(1x,3a)') '"',trim(val),'"'
end subroutine mprint_chr_
subroutine mprint_bul_(lu,who,what,val)
  implicit none
  integer,intent(in) :: lu
  character(len=*),intent(in) :: who,what
  logical,intent(in) :: val
  write(lu,'(3a)',advance='no') &
    trim(myid_(who)),'(): ',what
  write(lu,'(1x,l1)') val
end subroutine mprint_bul_
subroutine mprint_int_(lu,who,what,val,base)
  implicit none
  integer,intent(in) :: lu
  character(len=*),intent(in) :: who,what
  integer,intent(in) :: val
  character(len=*),optional,intent(in) :: base
        ! base is either z(hex), o(oct), b(binary), default(decimal)
  write(lu,'(3a)',advance='no') &
    trim(myid_(who)),'(): ',what
  if(present(base)) then
    select case(base)
    case('z','Z')
      write(lu,'(1x,z8.8)') val
    case('o','O')
      write(lu,'(1x,o11.11)') val
    case('b','B')
      write(lu,'(1x,b32.32)') val
    case default
      write(lu,'(1x,i0)') val
    endselect
  else
    write(lu,'(1x,i0)') val
  endif
end subroutine mprint_int_
subroutine mprint_vint_(lu,who,what,vals,base,format,sum)
  implicit none
  integer,intent(in) :: lu
  character(len=*),intent(in) :: who,what
  integer,dimension(:),intent(in) :: vals
  character(len=*),optional,intent(in) :: base
  character(len=*),optional,intent(in) :: format
  logical,optional,intent(in) :: sum
  integer:: i
  write(lu,'(3a)',advance='no') &
    trim(myid_(who)),'(): ',what
  if(present(sum)) then
    if(sum) then
      if(present(format)) then
        write(lu,fmt=format,advance='no') sumof_(vals(:))
      else
        write(lu,'(1x,i0)',advance='no') sumof_(vals(:))
      endif
    endif
  endif
  do i=1,size(vals)
  if(present(base)) then
    select case(base)
    case('z','Z')
      write(lu,'(1x,z8.8)',advance='no') vals(i)
    case('o','O')
      write(lu,'(1x,o11.11)',advance='no') vals(i)
    case('b','B')
      write(lu,'(1x,b32.32)',advance='no') vals(i)
    case default
      if(present(format)) then
        write(lu,fmt=format,advance='no') vals(i)
      else
        write(lu,'(1x,i0)',advance='no') vals(i)
      endif
    endselect
  else
    if(present(format)) then
      write(lu,fmt=format,advance='no') vals(i)
    else
      write(lu,'(1x,i0)',advance='no') vals(i)
    endif
  endif
  enddo
  write(lu,'()',advance='yes')
end subroutine mprint_vint_
function sumof_(vals)
  implicit none
  integer,dimension(:),intent(in) :: vals
  integer:: sumof_
  sumof_=sum(vals(:))
end function sumof_
subroutine mprint_flt_(lu,who,what,val)
  implicit none
  integer,intent(in) :: lu
  character(len=*),intent(in) :: who,what
  real(SP),intent(in) :: val
  write(lu,'(3a)',advance='no') &
    trim(myid_(who)),'(): ',what
  write(lu,'(1x,g15.7)') val
end subroutine mprint_flt_
subroutine mprint_dbl_(lu,who,what,val)
  implicit none
  integer,intent(in) :: lu
  character(len=*),intent(in) :: who,what
  real(DP),intent(in) :: val
  write(lu,'(3a)',advance='no') &
    trim(myid_(who)),'(): ',what
  write(lu,'(1x,g23.15)') val
end subroutine mprint_dbl_
subroutine mprint_err_(lu,who,what,errm)
  implicit none
  integer,intent(in) :: lu
  character(len=*),intent(in) :: who,what,errm
  write(lu,'(3a)',advance='no') &
    trim(myid_(who)),'(): >>> ERROR <<< ',what
  write(lu,'(1x,3a)') '"',trim(errm),'"'
end subroutine mprint_err_

subroutine perr_(who,what)
  implicit none
  character(len=*),intent(in) :: who,what
  call mprint(stderr,who,'>>> ERROR <<< '//what)
  if(stderr==stdout) return
  call mprint(stdout,who,'>>> ERROR <<< '//what)
end subroutine perr_

subroutine perr_chr_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what,val
  call mprint(stderr,who,'>>> ERROR <<< '//what,val)
  if(stderr==stdout) return
  call mprint(stdout,who,'>>> ERROR <<< '//what,val)
end subroutine perr_chr_
subroutine perr_bul_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  logical,intent(in) :: val
  call mprint(stderr,who,'>>> ERROR <<< '//what,val)
  if(stderr==stdout) return
  call mprint(stdout,who,'>>> ERROR <<< '//what,val)
end subroutine perr_bul_
subroutine perr_int_(who,what,val,base)
  implicit none
  character(len=*),intent(in) :: who,what
  integer,intent(in) :: val
  character(len=*),optional,intent(in) :: base
  call mprint(stderr,who,'>>> ERROR <<< '//what,val,base=base)
  if(stderr==stdout) return
  call mprint(stdout,who,'>>> ERROR <<< '//what,val,base=base)
end subroutine perr_int_
subroutine perr_vint_(who,what,vals,base,format,sum)
  implicit none
  character(len=*),intent(in) :: who,what
  integer,dimension(:),intent(in) :: vals
  character(len=*),optional,intent(in) :: base
  character(len=*),optional,intent(in) :: format
  logical,optional,intent(in) :: sum
  call mprint(stderr,who,'>>> ERROR <<< '//what,vals,base=base,format=format,sum=sum)
  if(stderr==stdout) return
  call mprint(stdout,who,'>>> ERROR <<< '//what,vals,base=base,format=format,sum=sum)
end subroutine perr_vint_
subroutine perr_flt_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  real(SP),intent(in) :: val
  call mprint(stderr,who,'>>> ERROR <<< '//what,val)
  if(stderr==stdout) return
  call mprint(stdout,who,'>>> ERROR <<< '//what,val)
end subroutine perr_flt_
subroutine perr_dbl_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  real(DP),intent(in) :: val
  call mprint(stderr,who,'>>> ERROR <<< '//what,val)
  if(stderr==stdout) return
  call mprint(stdout,who,'>>> ERROR <<< '//what,val)
end subroutine perr_dbl_

subroutine warn_(who,what)
  implicit none
  character(len=*),intent(in) :: who,what
  call mprint(stdout,who,'>>> WARNING <<< '//what)
end subroutine warn_
subroutine warn_chr_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what,val
  call mprint(stdout,who,'>>> WARNING <<< '//what,val)
end subroutine warn_chr_
subroutine warn_bul_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  logical,intent(in) :: val
  call mprint(stdout,who,'>>> WARNING <<< '//what,val)
end subroutine warn_bul_
subroutine warn_int_(who,what,val,base)
  implicit none
  character(len=*),intent(in) :: who,what
  integer,intent(in) :: val
  character(len=*),optional,intent(in) :: base
  call mprint(stdout,who,'>>> WARNING <<< '//what,val,base=base)
end subroutine warn_int_
subroutine warn_vint_(who,what,vals,base,format,sum)
  implicit none
  character(len=*),intent(in) :: who,what
  integer,dimension(:),intent(in) :: vals
  character(len=*),optional,intent(in) :: base
  character(len=*),optional,intent(in) :: format
  logical,optional,intent(in) :: sum
  call mprint(stdout,who,'>>> WARNING <<< '//what,vals,base=base,format=format,sum=sum)
end subroutine warn_vint_
subroutine warn_flt_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  real(SP),intent(in) :: val
  call mprint(stdout,who,'>>> WARNING <<< '//what,val)
end subroutine warn_flt_
subroutine warn_dbl_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  real(DP),intent(in) :: val
  call mprint(stdout,who,'>>> WARNING <<< '//what,val)
end subroutine warn_dbl_

subroutine tell_(who,what)
  implicit none
  character(len=*),intent(in) :: who,what
  call mprint(stdout,who,what)
end subroutine tell_
subroutine tell_chr_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what,val
  call mprint(stdout,who,what,val)
end subroutine tell_chr_
subroutine tell_bul_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  logical,intent(in) :: val
  call mprint(stdout,who,what,val)
end subroutine tell_bul_
subroutine tell_int_(who,what,val,base)
  implicit none
  character(len=*),intent(in) :: who,what
  integer,intent(in) :: val
  character(len=*),optional,intent(in) :: base
  call mprint(stdout,who,what,val,base)
end subroutine tell_int_
subroutine tell_vint_(who,what,vals,base,format,sum)
  implicit none
  character(len=*),intent(in) :: who,what
  integer,dimension(:),intent(in) :: vals
  character(len=*),optional,intent(in) :: base
  character(len=*),optional,intent(in) :: format
  logical,optional,intent(in) :: sum
  call mprint(stdout,who,what,vals,base=base,format=format,sum=sum)
end subroutine tell_vint_
subroutine tell_flt_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  real(SP),intent(in) :: val
  call mprint(stdout,who,what,val)
end subroutine tell_flt_
subroutine tell_dbl_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  real(DP),intent(in) :: val
  call mprint(stdout,who,what,val)
end subroutine tell_dbl_

subroutine dropdead_()
#ifdef USE_MPI_ABORT
  use mpeu_mpif, only: mpi_comm_world
#endif
  implicit none
  integer:: ier
  integer,parameter:: myer=2

  character(len=08):: cdate
  character(len=10):: ctime
  character(len=05):: czone

  call date_and_time(date=cdate,time=ctime,zone=czone)
  call mprint(stdout,'dropdead',ctime//'('//czone//') '//cdate)
  call mprint(stderr,'dropdead',ctime//'('//czone//') '//cdate)

#ifdef USE_MPI_ABORT
  call mpi_abort(mpi_comm_world,myer,ier)
#else
  ! this is a GSI_GridComp::abort()
# ifdef USE_GSI_ABOR1
    call abor1(myname//"::dropdead_()")
# else
    call exit(myer)
# endif
#endif
end subroutine dropdead_

subroutine die_(who)
  implicit none
  character(len=*),intent(in) :: who
  call perr_(who,'terminated')
  call dropdead_()
end subroutine die_
subroutine die2_(who,what)
  implicit none
  character(len=*),intent(in) :: who,what
  call perr_(who,what)
  call dropdead_()
end subroutine die2_
subroutine die_chr_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what,val
  call perr_chr_(who,what,val)
  call dropdead_()
end subroutine die_chr_
subroutine die_bul_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  logical,intent(in) :: val
  call perr_bul_(who,what,val)
  call dropdead_()
end subroutine die_bul_
subroutine die_int_(who,what,val,base)
  implicit none
  character(len=*),intent(in) :: who,what
  integer,intent(in) :: val
  character(len=*),optional,intent(in) :: base
  call perr_int_(who,what,val,base=base)
  call dropdead_()
end subroutine die_int_
subroutine die_vint_(who,what,vals,base,format,sum)
  implicit none
  character(len=*),intent(in) :: who,what
  integer,dimension(:),intent(in) :: vals
  character(len=*),optional,intent(in) :: base
  character(len=*),optional,intent(in) :: format
  logical,optional,intent(in) :: sum
  call perr_vint_(who,what,vals,base=base,format=format,sum=sum)
  call dropdead_()
end subroutine die_vint_
subroutine die_flt_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  real(SP),intent(in) :: val
  call perr_flt_(who,what,val)
  call dropdead_()
end subroutine die_flt_
subroutine die_dbl_(who,what,val)
  implicit none
  character(len=*),intent(in) :: who,what
  real(DP),intent(in) :: val
  call perr_dbl_(who,what,val)
  call dropdead_()
end subroutine die_dbl_
subroutine assert_(str,from,line)
  implicit none
  character(len=*),intent(in) :: str	! a message of assert_()
  character(len=*),intent(in) :: from	! where assert_() is invoked.
  integer,         intent(in) :: line	! where assert_() is invoked.
  character(len=*),parameter :: myname_='ASSERT_'
  call perr(myname_,'failed',str)
  call perr(myname_,'file =',from)
  call perr(myname_,'line #',line)
  call die(myname_)
end subroutine assert_
subroutine assert_GE_(m,n,who,str)
  implicit none
  integer,intent(in) :: m,n
  character(len=*),intent(in) :: who	! where assert_GE_() is invoked.
  character(len=*),intent(in) :: str	! a message of assert_GE_()
  character(len=*),parameter :: myname_='ASSERT_GE_'
  if(.not.(m>=n)) then
    call perr(myname_,'test failed',str)
    call perr(myname_,'operand 1 = ',m)
    call perr(myname_,'operand 2 = ',n)
    call die(myname_,who)
  endif
end subroutine assert_GE_
end module satinfo_util
