!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_actvchan - Input module of active-channels table
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_actvchan
      implicit none
      private	! except

      public :: actvchan_select
      public :: actvchan_clean
      public :: actvchan_show

      interface actvchan_select; module procedure &
        select_; end interface
      interface actvchan_clean; module procedure &
        dealloc_; end interface
      interface actvchan_show; module procedure &
        show_; end interface

! !REVISION HISTORY:
! 	06Jul09	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- found a suspected compiler error, and changed
!		  the application of "append=" optional argument
!		  to have a work arround.
! 	04Jun07	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_actvchan'

!!$ Usecase:
!!$
!!$   character(len=LEN_WORD),pointer,dimension(:) :: p_nm
!!$   integer,pointer,dimension(:) :: p_ch
!!$   integer :: n
!!$   call actvchan_select(p_nm,p_sn,p_ch,n, &
!!$                        from="active_channels.tbl",	 &
!!$                        where_dt=(/19990330,060000/), &
!!$                        where_nm=(/"f10","f11"/), &
!!$			   ...)

!!$   [.. p_nm(1:n),p_ch(1:n) ..]
!!$   call actvchan_clean(p_nm,p_ch,n)
!!$

contains

subroutine select_(p_nm,p_sn,p_ch,n, from,where_dt,	&
		   where_nm,where_sn, dbase,bufrsize,	&
		   append,verbose)
  use satinfo_util,only : alloc
  implicit none
  character(len=*),pointer,dimension(:) :: p_nm
  character(len=*),pointer,dimension(:) :: p_sn
  integer         ,pointer,dimension(:) :: p_ch
  integer,intent(inout) :: n

  character(len=*),intent(in) :: from

  integer         ,dimension(:),         intent(in) :: where_dt
  character(len=*),dimension(:),optional,intent(in) :: where_nm
  character(len=*),dimension(:),optional,intent(in) :: where_sn

  character(len=*),optional,intent(in) :: dbase  ! default ".".
  integer         ,optional,intent(in) :: bufrsize ! default 256 bytes
  logical         ,optional,intent(in) :: append !
  logical         ,optional,intent(in) :: verbose ! verbal messages

  integer :: lbufr
  logical :: append_

  append_=.false.; if(present(append)) append_=append

  if(.not.append_) then
    call alloc(p_nm,n)
    call alloc(p_sn,n)
    call alloc(p_ch,n)
  endif

  lbufr=256
  if(present(bufrsize)) lbufr=bufrsize
  if(present(dbase)) then
    call selectb_(p_nm,p_sn,p_ch,n,from,dbase,lbufr,		&
      where_dt=where_dt,where_nm=where_nm,where_sn=where_sn,	&
      verbose=verbose)
  else
    call selectb_(p_nm,p_sn,p_ch,n,from,'.'  ,lbufr,		&
    where_dt=where_dt,where_nm=where_nm,where_sn=where_sn,	&
      verbose=verbose)
  endif
end subroutine select_

subroutine selectb_(p_nm,p_sn,p_ch,n,from,dbase,bufrsize, &
		    where_dt,where_nm,where_sn, verbose)

  use satinfo_util,only : perr,die,tell
  use satinfo_util,only : luavail,getarec
  use satinfo_util,only : realloc
  use satinfo_util,only : locate
  implicit none
  character(len=*),pointer,dimension(:) :: p_nm
  character(len=*),pointer,dimension(:) :: p_sn
  integer         ,pointer,dimension(:) :: p_ch
  integer,intent(inout) :: n

  character(len=*),intent(in) :: from
  character(len=*),intent(in) :: dbase
  integer         ,intent(in) :: bufrsize

  character(len=*),dimension(:),optional,intent(in) :: where_nm
  integer,         dimension(:),optional,intent(in) :: where_dt
  character(len=*),dimension(:),optional,intent(in) :: where_sn

  logical,optional,intent(in) :: verbose ! verbal messages

  character(len=*),parameter :: myname_=myname//"::selectb_"

  integer :: lu,ier,j,k,m
  integer :: iymd,ihms
  integer :: nymd,nhms
  integer :: lymd,lhms
  character(len=len(p_nm)) :: nm
  character(len=len(p_sn)) :: sn
  character(len=bufrsize) :: bufr

  logical :: nm_matched, nm_any
  logical :: dt_matched, dt_any
  logical :: sn_matched, sn_any

  logical :: verb_
  verb_  =.false.; if(present(verbose))  verb_=verbose

  nm_any=.true.; if(present(where_nm)) nm_any=.false.
  dt_any=.true.; if(present(where_dt)) dt_any=.false.
  sn_any=.true.; if(present(where_sn)) sn_any=.false.

  lu=luavail()
  open(lu,file=trim(dbase)//'/'//trim(from),status='old',iostat=ier)
    if(ier/=0) call die(myname_, &
      'open("'//trim(dbase)//'/'//trim(from)//'"), iostat =',ier)
    if(verb_) call tell(myname_, &
      'open("'//trim(dbase)//'/'//trim(from)//'") for input')

  nymd=-1
  nhms=-1
  if(.not.dt_any) then
    nymd=where_dt(1)
    nhms=where_dt(2)
  endif

  call getarec(lu,bufr,ier)
  do while(ier==0)
    read(bufr,*,iostat=ier) nm,iymd,ihms,lymd,lhms,sn,m
      if(ier/=0) then
        call perr(myname_,'bufr="'//trim(bufr)//'"')
        call die(myname_, &
	  'stem-read("'//trim(dbase)//'/'//trim(from)//'"), iostat =',ier)
      endif

    nm_matched = nm_any
    dt_matched = dt_any
    sn_matched = sn_any

    	! xx_matched can be modified only if present(where_xx).

    if(.not.nm_any) then
      nm_matched = locate(nm,where_nm(:))>0
    endif
    if(.not.dt_any) dt_matched = &
       (iymd<nymd .or. (iymd==nymd.and.ihms<=nhms)) .and. &
       (nymd<lymd .or. (nymd==lymd.and.nhms<=lhms))
    if(.not.sn_any) then
      sn_matched = locate(sn,where_sn(:))>0
    endif

    if( nm_matched .and. dt_matched .and. sn_matched ) then

      call realloc(p_nm,n,incr=m)
      call realloc(p_sn,n,incr=m)
      call realloc(p_ch,n,incr=m)

      p_nm(n+1:n+m)=nm
      p_sn(n+1:n+m)=sn
      read(bufr,*,iostat=ier) nm,iymd,ihms,lymd,lhms,sn,k, &
                              p_ch(n+1:n+m),k
        if(ier==0) then
          call perr(myname_,'bufr="'//trim(bufr)//'"')
          call die(myname_, &
	    'leaf-read("'//trim(dbase)//'/'//trim(from)//'"), too many channels',m+1)
	endif

      read(bufr,*,iostat=ier) nm,iymd,ihms,lymd,lhms,sn,k, &
                              p_ch(n+1:n+m)
        if(ier/=0) then
          call perr(myname_,'bufr="'//trim(bufr)//'"')
          call die(myname_, &
	    'leaf-read("'//trim(dbase)//'/'//trim(from)//'"), iostat =',ier)
        endif

      n=n+m
    endif

        ! potential read() errors are ignored
    call getarec(lu,bufr,ier)
  enddo

  close(lu,iostat=ier)
    if(ier/=0) call die(myname_, &
      'close("'//trim(dbase)//'/'//trim(from)//'"), iostat =',ier)

!  if(verb_) call show_(p_nm,p_sn,p_ch,n,name=trim(dbase)//'/'//trim(from))
end subroutine selectb_

subroutine dealloc_(p_nm,p_sn,p_ch,n)
  use satinfo_util,only : dealloc
  implicit none
  character(len=*),pointer,dimension(:) :: p_nm
  character(len=*),pointer,dimension(:) :: p_sn
  integer         ,pointer,dimension(:) :: p_ch
  integer,intent(out) :: n

  call dealloc(p_nm,n)
  call dealloc(p_sn,n)
  call dealloc(p_ch,n)
end subroutine dealloc_

subroutine show_(p_nm,p_sn,p_ch,n,name)
  use satinfo_util,only : tell
  implicit none
  character(len=*),dimension(:),intent(in) :: p_nm
  character(len=*),dimension(:),intent(in) :: p_sn
  integer         ,dimension(:),intent(in) :: p_ch
  integer,intent(in) :: n
  character(len=*),optional,intent(in) :: name

  integer :: i

  if(present(name)) call tell(myname,'name = "'//trim(name)//'"')
  call tell(myname,'n =',n)
  if(n>0) call tell(myname,'i, sn, nm, ch')

  do i=1,n
    print'(a20,i8.6)', trim(p_sn(i))//' '//trim(p_nm(i)),p_ch(i)
  end do
end subroutine show_
end module m_actvchan
