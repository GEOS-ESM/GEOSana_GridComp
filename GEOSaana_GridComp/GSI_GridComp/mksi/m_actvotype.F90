!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_actvotype - Input module of active-obstypes table
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_actvotype
      implicit none
      private	! except

      public :: actvotype_select
      public :: actvotype_clean
      public :: actvotype_show

      interface actvotype_select; module procedure &
        select_; end interface
      interface actvotype_clean; module procedure &
        dealloc_; end interface
      interface actvotype_show; module procedure &
        show_; end interface

! !REVISION HISTORY:
! 	06Jul09	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- found a suspected compiler error, and changed
!		  the application of "append=" optional argument
!		  to have a work arround.
! 	04Jun07	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_actvotype'

!!$ Usecase:
!!$
!!$   integer,pointer,dimension(:) :: p_nm		! "type" in convinfo.txt
!!$   integer,pointer,dimension(:) :: p_sb		! "sub" in convinfo.txt
!!$   character(len=_LEN_),pointer,dimension(:) :: p_ot	! "otype" in convinfo.txt
!!$   integer :: n
!!$   call actvotype_select(p_nm,p_sb,p_ot,n, &
!!$                        from="active_otypes.tbl",	 &
!!$                        where_dt=(/19990330,060000/), &
!!$                        where_nm=(/120,243,253,254/), &
!!$			   ...)
!!$
!!$   [.. p_nm(1:n),p_ot(1:n) ..]
!!$   call actvotype_clean(p_nm,p_sb,p_ot,n)
!!$

contains

subroutine select_(p_nm,p_sb,p_ot,n, from,where_dt,	&
		   where_nm,where_sb, dbase,bufrsize,	&
		   append,verbose)
  use satinfo_util,only : alloc
  implicit none
  integer         ,pointer,dimension(:) :: p_nm
  integer         ,pointer,dimension(:) :: p_sb
  character(len=*),pointer,dimension(:) :: p_ot
  integer,intent(inout) :: n

  character(len=*),intent(in) :: from

  integer         ,dimension(:),         intent(in) :: where_dt
  integer         ,dimension(:),optional,intent(in) :: where_nm
  integer         ,dimension(:),optional,intent(in) :: where_sb

  character(len=*),optional,intent(in) :: dbase  ! default ".".
  integer         ,optional,intent(in) :: bufrsize ! default 256 bytes
  logical         ,optional,intent(in) :: append !
  logical         ,optional,intent(in) :: verbose ! verbal messages

  integer :: lbufr
  logical :: append_

  append_=.false.; if(present(append)) append_=append

  if(.not.append_) then
    call alloc(p_nm,n)
    call alloc(p_sb,n)
    call alloc(p_ot,n)
  endif

  lbufr=256
  if(present(bufrsize)) lbufr=bufrsize
  if(present(dbase)) then
    call selectb_(p_nm,p_sb,p_ot,n,from,dbase,lbufr,		&
      where_dt=where_dt,where_nm=where_nm,where_sb=where_sb,	&
      verbose=verbose)
  else
    call selectb_(p_nm,p_sb,p_ot,n,from,'.'  ,lbufr,		&
    where_dt=where_dt,where_nm=where_nm,where_sb=where_sb,	&
      verbose=verbose)
  endif
end subroutine select_

subroutine selectb_(p_nm,p_sb,p_ot,n,from,dbase,bufrsize, &
		    where_dt,where_nm,where_sb, verbose)

  use satinfo_util,only : perr,die,tell
  use satinfo_util,only : luavail,getarec
  use satinfo_util,only : realloc
  use satinfo_util,only : locate
  implicit none
  integer         ,pointer,dimension(:) :: p_nm
  integer         ,pointer,dimension(:) :: p_sb
  character(len=*),pointer,dimension(:) :: p_ot
  integer,intent(inout) :: n

  character(len=*),intent(in) :: from
  character(len=*),intent(in) :: dbase
  integer         ,intent(in) :: bufrsize

  integer,         dimension(:),optional,intent(in) :: where_dt
  integer,         dimension(:),optional,intent(in) :: where_nm
  integer,         dimension(:),optional,intent(in) :: where_sb

  logical,optional,intent(in) :: verbose ! verbal messages

  character(len=*),parameter :: myname_=myname//"::selectb_"

  integer :: lu,ier,j,k,m,nrec,irec
  integer :: iymd,ihms
  integer :: nymd,nhms
  integer :: lymd,lhms
  integer :: nm,sb
  character(len=len(p_ot)) :: ot
  character(len=bufrsize) :: bufr

  logical :: nm_matched, nm_any
  logical :: dt_matched, dt_any
  logical :: sb_matched, sb_any

  logical :: verb_
  verb_  =.false.; if(present(verbose))  verb_=verbose

  nm_any=.true.; if(present(where_nm)) nm_any=.false.
  dt_any=.true.; if(present(where_dt)) dt_any=.false.
  sb_any=.true.; if(present(where_sb)) sb_any=.false.

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

  irec=0
  call getarec(lu,bufr,ier,nrec=nrec)
  do while(ier==0)
    irec=irec+nrec
    read(bufr,*,iostat=ier) sb,iymd,ihms,lymd,lhms,nm,m
      if(ier/=0) then
        call perr(myname_,'stem-read(bufr) error, iostat =',ier)
        call perr(myname_,'                     at # rec =',irec)
        call perr(myname_,'   bufr =',trim(bufr))
	call perr(myname_,'   file =',trim(dbase)//'/'//trim(from))
        call perr(myname_,'expected fields are, "sb,iymd,ihms,lymd,lhms,nm,m"')
        call  die(myname_)
      endif

    nm_matched = nm_any
    dt_matched = dt_any
    sb_matched = sb_any

    	! xx_matched can be modified only if present(where_xx).

    if(.not.nm_any) then
      nm_matched = locate(nm,where_nm(:))>0
    endif
    if(.not.dt_any) dt_matched = &
       (iymd<nymd .or. (iymd==nymd.and.ihms<=nhms)) .and. &
       (nymd<lymd .or. (nymd==lymd.and.nhms<=lhms))
    if(.not.sb_any) then
      sb_matched = locate(sb,where_sb(:))>0
    endif

    if( nm_matched .and. dt_matched .and. sb_matched ) then

      call realloc(p_nm,n,incr=m)
      call realloc(p_sb,n,incr=m)
      call realloc(p_ot,n,incr=m)

      p_nm(n+1:n+m)=nm
      p_sb(n+1:n+m)=sb
      read(bufr,*,iostat=ier) sb,iymd,ihms,lymd,lhms,nm,k, &
                              p_ot(n+1:n+m),ot
        if(ier==0.and.(ot(1:1)/='#'.and.ot(1:1)/='!')) then
          call perr(myname_,'leaf-read(bufr) error, too many otypes available, at rec# =',irec)
          call perr(myname_,'                               allowed total otypes count =',m) 
          call perr(myname_,'   bufr =',trim(bufr))
	  call perr(myname_,'   file =',trim(dbase)//'/'//trim(from))
          call  die(myname_)
	endif

      read(bufr,*,iostat=ier) sb,iymd,ihms,lymd,lhms,nm,k, &
                              p_ot(n+1:n+m)
        if(ier/=0) then
          call perr(myname_,'leaf-read(bufr) error, not enough otypes available, at rec# =',irec)
          call perr(myname_,'                              expected minimum otypes count =',m) 
          call perr(myname_,'   bufr =',trim(bufr))
	  call perr(myname_,'   file =',trim(dbase)//'/'//trim(from))
          call  die(myname_)
        endif

      n=n+m
    endif

        ! potential read() errors are ignored
    call getarec(lu,bufr,ier,nrec=nrec)
  enddo

  close(lu,iostat=ier)
    if(ier/=0) call die(myname_, &
      'close("'//trim(dbase)//'/'//trim(from)//'"), iostat =',ier)

!  if(verb_) call show_(p_nm,p_sb,p_ot,n,name=trim(dbase)//'/'//trim(from))
end subroutine selectb_

subroutine dealloc_(p_nm,p_sb,p_ot,n)
  use satinfo_util,only : dealloc
  implicit none
  integer         ,pointer,dimension(:) :: p_nm
  integer         ,pointer,dimension(:) :: p_sb
  character(len=*),pointer,dimension(:) :: p_ot
  integer,intent(out) :: n

  call dealloc(p_nm,n)
  call dealloc(p_sb,n)
  call dealloc(p_ot,n)
end subroutine dealloc_

subroutine show_(p_nm,p_sb,p_ot,n,name)
  use satinfo_util,only : tell
  implicit none
  integer         ,dimension(:),intent(in) :: p_nm
  integer         ,dimension(:),intent(in) :: p_sb
  character(len=*),dimension(:),intent(in) :: p_ot
  integer,intent(in) :: n
  character(len=*),optional,intent(in) :: name

  integer :: i

  if(present(name)) call tell(myname,'name = "'//trim(name)//'"')
  call tell(myname,'n =',n)

  if(n>0) call tell(myname,'i, sub, itype, otype')
  do i=1,n
    print'(i4,i5,a,i3.3,1x,a20)', i, p_sb(i),'_',p_nm(i),trim(p_ot(i))
  end do
end subroutine show_
end module m_actvotype
