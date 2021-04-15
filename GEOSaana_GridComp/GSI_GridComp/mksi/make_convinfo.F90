!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: make_convinfo - create a convinfo file for a given date-time
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"
#include "assert.H"

    program make_convinfo
      use m_actvotype  , only : actv_select => actvotype_select
      use m_actvotype  , only : actv_clean  => actvotype_clean
      use m_convtmpl   , only : NAMESIZE
      use m_convtmpl   , only : convtmpl
      use m_convtmpl   , only : convtmpl_get,convtmpl_put,convtmpl_dealloc
      use m_convtmpl   , only : ptr_iuse
      use satinfo_util, only : tell,warn,perr,die,stdin,stdout
      implicit none

! !REVISION HISTORY:
!	09Nov11	- Jing Guo <jing.guo@nasa.gov>
!		- Adapted from make_satinfo, for convinfo management.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='make_convinfo'
  integer,parameter :: PATHSIZE=512
  integer,parameter :: BUFRSIZE=132

  integer,parameter :: default_NOTUSED   = -9
  integer,parameter :: default_AVAILABLE = -1
  integer,parameter :: default_MONITOR   = +0
  integer,parameter :: default_ACTIVE    = +1

  character(len=*),parameter :: default_AVAILABLE_TBL  = 'available.tbl'
  character(len=*),parameter :: default_ACTIVE_TBL     =    'active.tbl'
  character(len=*),parameter :: default_MONITOR_TBL    =   'monitor.tbl'

! Elements of convinfo database tables, an in-stream
  integer                ,pointer,dimension(:) :: xxxx_nm => null()
  integer                ,pointer,dimension(:) :: xxxx_sb => null()
  character(len=NAMESIZE),pointer,dimension(:) :: xxxx_ot => null()
  integer :: nxxxx=-1

! A given convinfo tempalte table, an inout-stream
  type(convTmpl),pointer,dimension(:):: tmpl => null()
  integer :: ntmpl=-1

! Configuration (setup) parameters.

    ! [nymd, nhms] must be defined explicitly.
  integer :: nymd
  integer :: nhms

    ! optional parameters.  Default values can be used.
  character(len=PATHSIZE) :: info_tmpl = "convinfo.tmpl"
  character(len=PATHSIZE) :: info_outf = "convinfo.txt"
  character(len=PATHSIZE) :: dbname        = "./convinfo.db"

  character(len=PATHSIZE) :: available_tbl = default_AVAILABLE_TBL
  character(len=PATHSIZE) ::    active_tbl = default_ACTIVE_TBL
  character(len=PATHSIZE) ::   monitor_tbl = default_MONITOR_TBL

  logical :: verbose=.false.
  logical :: nowarn =.false.

  integer ::   NOTUSED_FLAG = default_NOTUSED
  integer :: AVAILABLE_FLAG = default_AVAILABLE
  integer ::    ACTIVE_FLAG = default_ACTIVE
  integer ::   MONITOR_FLAG = default_MONITOR

! The port of /setup/ NAMELIST
  namelist/setup/ nymd,nhms,dbname,	&
    info_tmpl,info_outf,	&
    available_tbl,	&
    active_tbl,	&
    monitor_tbl,	&
    verbose,		&
    NOTUSED_FLAG,	&
    AVAILABLE_FLAG,	&
    ACTIVE_FLAG,	&
    MONITOR_FLAG,	&
    nowarn

! Internal working variables
  integer :: itmpl
  integer :: iver
  integer :: ios
  logical :: warning

  integer :: iuse_tbl
  integer :: itbl
  logical :: report
  logical :: isOptional
  character(len=PATHSIZE) :: mtag_tbl
  character(len=PATHSIZE) :: this_tbl
  integer,pointer,dimension(:) :: tmpl_iu

!!$ Get required and optional control arguments:
!!$			
  nymd = -HUGE(nymd)
  nhms = -HUGE(nhms)

  read(stdin,setup,iostat=ios)
  	if(ios/=0) call die(myname,'read(/setup/), iostat',ios)

    !!$ Verify /setup/ input for required arguments:
    !!$		nymd: <yyyymmdd>, e.g. 20050101 for 01/01/05
    !!$		nhms:   <hhmmss>, e.g.   083000 for 08:30:00

      if( nymd == -HUGE(nymd) .or. nhms == -HUGE(nhms) ) then
	if( nymd==-HUGE(nymd) ) call perr(myname,'undefined nymd')
	if( nhms==-HUGE(nhms) ) call perr(myname,'undefined nhms')
  	call die(myname,'read(/setup/)')
      endif

  call tell(myname,'processing for [nymd, nhms] = ',(/nymd,nhms/),format='(i9.6)')

!!$ Mark tmpl-entry list according to the active-entry list.

  call convTmpl_get(info_tmpl, ntmpl,tmpl)

!!$ <<< Initialization pass >>> All levels in the template are turned off
  tmpl_iu => ptr_iuse (tmpl,ubnd=ntmpl)
  do itmpl=1,ntmpl
    tmpl_iu(itmpl)=min(tmpl_iu(itmpl),NOTUSED_FLAG)
  enddo

!!$ <<< Pass 0:n >>> All levels are reset if matched by a table entry

  ! iuse_flag values and priorities (latter override earliers)
  ! < -1, keep as is (not used)                          everything in convinfo.tmpl
  ! = -1, O-F diagnostics produced                       available entries ("passive")
  ! =  1, use data with complete quality control         active entries
  ! =  0, monitoring and use in QC only                  monitoring entries

  do itbl=0,2
    report=.false.
    warning=.not.nowarn
    isOptional=.false.
    select case(itbl)
    case(0)
	!!$ <<< Pass 1 >>> Set all levels PASSIVE
      iuse_tbl=AVAILABLE_FLAG
      warning =.true.
      mtag_tbl="passive"
      this_tbl=available_tbl

    case(1)
      iuse_tbl=ACTIVE_FLAG
      mtag_tbl="active"
      this_tbl=active_tbl

    case(2)
      iuse_tbl=MONITOR_FLAG
      report  =.true.
      isOptional=.true.
      mtag_tbl='monitor'
      this_tbl=monitor_tbl

    case default
      call die(myname,'unknown table file, #',itbl)
    end select
    call tell(myname,'    itbl =',itbl)
    call tell(myname,'mtag_tbl =',mtag_tbl)

    if( isOptional .and. &
    	.not. file_exist_(this_tbl,dbase=dbname)) cycle
    
	!!$ Read the channel list
    call actv_select(xxxx_nm,xxxx_sb,xxxx_ot,nxxxx,	&
      where_dt=(/nymd,nhms/), from=this_tbl,		&
      dbase=dbname,verbose=.true.)

	!!$ Set all levels to iuse_tbl
    call tmpl_xcheck(iuse_tbl,tmpl,ntmpl, &
      xxxx_nm,xxxx_sb,xxxx_ot,nxxxx, tag=mtag_tbl, &
      warning=warning,report=report)

	!!$ reset working arrays
    call actv_clean(xxxx_nm,xxxx_sb,xxxx_ot,nxxxx)
  enddo

  nullify(tmpl_iu )

!!$ write out the convinfo table.
  call convtmpl_put(info_outf,ntmpl,tmpl, nymd=nymd,nhms=nhms)
  call convtmpl_dealloc(tmpl,ntmpl)

contains
function file_exist_(name,dbase)
  implicit none
  character(len=*),intent(in) :: name
  character(len=*),optional,intent(in) :: dbase
  logical :: file_exist_
  file_exist_=.false.
  if(trim(name)/="") then
    if(present(dbase)) then
      !!$ print*,0,trim(dbase)//'/'//trim(name)
      inquire(file=trim(dbase)//'/'//trim(name),exist=file_exist_)
    else
      !!$ print*,1,trim(name)
      inquire(file=name,exist=file_exist_)
    endif
  endif
end function file_exist_
subroutine tmpl_xcheck(inuse,tmpl,ntmpl, &
    xxxx_nm,xxxx_sb,xxxx_ot,nxxxx, tag, report, warning)

  use satinfo_util,only : assert_,warn,tell
  use m_convTmpl, only: ptr_type, ptr_sub, ptr_otype, ptr_iuse
  implicit none
  integer,intent(in) :: inuse	! one of channel-in-use flags:
  				! 	  NotUsed_FLAG
				!	AVAILABLE_FLAG
				!	   ACTIVE_FLAG
				!	  MONITOR_FLAG

  type(convtmpl),target,dimension(:),intent(inout) :: tmpl
  integer,intent(in) :: ntmpl

  integer         ,dimension(:),intent(in) :: xxxx_nm
  integer         ,dimension(:),intent(in) :: xxxx_sb
  character(len=*),dimension(:),intent(in) :: xxxx_ot
  integer,intent(in) :: nxxxx
  
  integer                ,pointer,dimension(:) :: tmpl_nm
  integer                ,pointer,dimension(:) :: tmpl_sb
  character(len=NAMESIZE),pointer,dimension(:) :: tmpl_ot
  integer                ,pointer,dimension(:) :: tmpl_iu

  character(len=*),intent(in) :: tag
  logical,optional,intent(in) :: report
  logical,optional,intent(in) :: warning
  
  character(len=BUFRSIZE) :: xxxx_mesg
  logical :: xxxx_matched,report_,warning_
  integer :: ixxxx,itmpl

  report_   =.false.; if(present(report)) report_=report
  warning_  =.false.; if(present(warning)) warning_=warning

  call tell(myname,'nxxxx =',nxxxx)
  call tell(myname,'ntmpl =',ntmpl)

  tmpl_nm =>ptr_type (tmpl,ubnd=ntmpl)
  tmpl_sb =>ptr_sub  (tmpl,ubnd=ntmpl)
  tmpl_ot =>ptr_otype(tmpl,ubnd=ntmpl)
  tmpl_iu =>ptr_iuse (tmpl,ubnd=ntmpl)

  do ixxxx=1,nxxxx

    if(report_.or.warning_) then
      write(xxxx_mesg,'(a,i2,3a,i2,a,i3.3,3a)') 'iuse ',inuse,	&
      	' ("',trim(tag),'") -> [',xxxx_sb(ixxxx),'_',		&
	xxxx_nm(ixxxx),', ',trim(xxxx_ot(ixxxx)),']'
    endif

        ! With a given (sb,nm,ot) entry in the active list, search the
        ! template list for a match.  If a match is found, set the
	! given inuse value to the matching entry in the template list.

    xxxx_matched=.false.
    do itmpl=1,ntmpl
      xxxx_matched = xxxx_sb(ixxxx)==tmpl_sb(itmpl) .and. &
		     xxxx_nm(ixxxx)==tmpl_nm(itmpl) .and. &
		     xxxx_ot(ixxxx)==tmpl_ot(itmpl) 

      if(xxxx_matched) then
        tmpl_iu(itmpl)=inuse ! mark this tmpl[%sb,%nm,%ot] to tmpl[inuse]

	if(report_) call tell(myname,'set '//trim(xxxx_mesg))

        exit	! no need to continue this loop of (itmpl=1,ntmpl)
      endif
    enddo

    if((.not.xxxx_matched) .and. warning_) then
      call warn(myname,'no entry found, '//trim(xxxx_mesg))
      cycle	! continue this loop of (ixxxx=1,nxxxx)
    endif
  enddo

  nullify(tmpl_nm)
  nullify(tmpl_sb)
  nullify(tmpl_ot)
  nullify(tmpl_iu)
  call tell(myname,'exiting, tag =',tag)
end subroutine tmpl_xcheck
end program make_convinfo
