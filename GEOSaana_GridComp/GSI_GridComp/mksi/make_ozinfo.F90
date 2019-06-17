!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: make_ozinfo - create a ozinfo file for a given date-time
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"
#include "assert.H"

    program make_ozinfo
      use m_actvchan  , only : actvlevs_select => actvchan_select
      use m_actvchan  , only : actvlevs_clean  => actvchan_clean
      use m_oztmpl    , only : NUSIS_SIZE => NAMESIZE
      use m_oztmpl    , only : oztmpl
      use m_oztmpl    , only : oztmpl_get,oztmpl_put,oztmpl_dealloc
      use m_oztmpl    , only : ptr_nusis,ptr_nulev
      use m_oztmpl    , only : ptr_inuse => ptr_iuse
      use satinfo_util, only : tell,warn,perr,die,stdin,stdout
      implicit none

! !REVISION HISTORY:
!	09Nov11	- Jing Guo <jing.guo@nasa.gov>
!		- Adapted from make_satinfo, for ozinfo management.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='make_ozinfo'
  integer,parameter :: PATHSIZE=512
  integer,parameter :: DESCSIZE=80

  character(len=NUSIS_SIZE),pointer,dimension(:) :: xxxx_nm => null()
  character(len=NUSIS_SIZE),pointer,dimension(:) :: xxxx_sn => null()
  integer,pointer,dimension(:) :: xxxx_ch => null()
  integer :: nchan=-1

  type(ozTmpl),pointer,dimension(:):: tmpl_ozi   => null()
  integer     ,pointer,dimension(:):: tmpl_inuse => null()
  integer :: ntmpl=-1

    ! key parameters.  Must be defined explicitly.
  integer :: nymd
  integer :: nhms

  integer,parameter :: default_NOTUSED   = -9
  integer,parameter :: default_AVAILABLE = -1
  integer,parameter :: default_MONITOR   = +0
  integer,parameter :: default_ACTIVE    = +1

  character(len=*),parameter :: default_AVAILABLE_TBL = 'available.tbl'
  character(len=*),parameter :: default_ACTIVE_TBL    =    'active.tbl'
  character(len=*),parameter :: default_MONITOR_TBL   =   'monitor.tbl'

    ! optional parameters.  Default values can be used.
  character(len=PATHSIZE) :: info_tmpl =  "ozinfo.tmpl"
  character(len=PATHSIZE) :: info_outf =  "ozinfo.txt"
  character(len=PATHSIZE) :: dbname    ="./ozinfo.db"
  character(len=PATHSIZE) :: available_tbl = default_AVAILABLE_TBL
  character(len=PATHSIZE) ::    active_tbl = default_ACTIVE_TBL
  character(len=PATHSIZE) ::   monitor_tbl = default_MONITOR_TBL

  logical :: verbose=.false.
  logical :: nowarn =.false.

  integer ::   NOTUSED_FLAG = default_NOTUSED
  integer :: AVAILABLE_FLAG = default_AVAILABLE
  integer ::    ACTIVE_FLAG = default_ACTIVE
  integer ::   MONITOR_FLAG = default_MONITOR

  integer :: itmpl
  integer :: iver
  integer :: ios
  logical :: warning

  integer :: iuse_tbl
  integer :: itbl
  logical :: report
  logical :: isOptional
  character(len=PATHSIZE) :: mtag_tbl
  character(len=PATHSIZE) :: levs_tbl

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
      warning=.not.nowarn

  call tell(myname,'processing for [nymd, nhms] = ',(/nymd,nhms/),format='(i9.6)')

!!$ Mark tmpl-levels list according to the active-levels list.
!!$   for the same sisname-sischan, set rad_flag = yes

  call ozTmpl_get(info_tmpl,ntmpl, tmpl_ozi)

  tmpl_inuse  => ptr_inuse (tmpl_ozi,ubnd=ntmpl)

!!$ <<< Initialization pass >>> All levels in the template are turned off
!  where(tmpl_inuse(1:ntmpl)>NOTUSED_FLAG)
  do itmpl=1,ntmpl
    tmpl_inuse(itmpl)=min(tmpl_inuse(itmpl),NOTUSED_FLAG)
  enddo
!  endwhere

!!$ <<< Pass 0:n >>> All levels are reset if matched by a table entry

  ! iuse_flag values and priorities (latter override earliers)
  ! = -2 do not use                                             everything in ozinfo.tmpl
  ! = -1 monitor if diagnostics produced                        available levels
  ! =  1 use data with complete quality control                 active levels
  ! =  0 monitor and use in QC only                             monitoring levels

  do itbl=0,2
    report=.false.
    isOptional=.false.
    select case(itbl)
    case(0)
	!!$ <<< Pass 1 >>> Set all levels PASSIVE
      iuse_tbl=AVAILABLE_FLAG
      mtag_tbl="passive"
      levs_tbl=available_tbl

    case(1)
      iuse_tbl=ACTIVE_FLAG
      mtag_tbl="active"
      levs_tbl=active_tbl

    case(2)
      iuse_tbl=MONITOR_FLAG
      report=.true.
      isOptional=.true.
      mtag_tbl='monitor'
      levs_tbl=monitor_tbl

    case default
      call die(myname,'unknown table file, #',itbl)
    end select
    call tell(myname,'    itbl =',itbl)
    call tell(myname,'mtag_tbl =',mtag_tbl)

    if( isOptional .and. &
    	.not. file_exist_(levs_tbl,dbase=dbname)) cycle
    
	!!$ Read the channel list
    call actvlevs_select(xxxx_nm,xxxx_sn,xxxx_ch,nchan, &
      where_dt=(/nymd,nhms/),	&
      from=levs_tbl,dbase=dbname,verbose=.true.)

	!!$ Set all levels to iuse_tbl
    call tmpl_xcheck(iuse_tbl,ptr_nusis(tmpl_ozi,ubnd=ntmpl),ptr_nulev(tmpl_ozi,ubnd=ntmpl),tmpl_inuse,ntmpl, &
      xxxx_nm,xxxx_sn,xxxx_ch,nchan, tag=mtag_tbl, &
      warning=warning,report=report)

	!!$ reset working arrays
    call actvlevs_clean(xxxx_nm,xxxx_sn,xxxx_ch,nchan)
  enddo

  nullify(tmpl_inuse )

!!$ write out the ozinfo table.
  call oztmpl_put(info_outf,ntmpl, tmpl_ozi, nymd=nymd,nhms=nhms)
  call oztmpl_dealloc(tmpl_ozi,ntmpl)

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
subroutine tmpl_xcheck(inuse,tmpl_nusis,tmpl_nuchn,tmpl_inuse,ntmpl, &
    xxxx_nm,xxxx_sn,xxxx_ch,nchan, tag, report, warning)

  use satinfo_util,only : assert_,warn,tell
  implicit none
  integer,intent(in) :: inuse	! one of channel-in-use flags:
  				! 	  NotUsed_FLAG
				!	AVAILABLE_FLAG
				!	   ACTIVE_FLAG
				!	  MONITOR_FLAG

  character(len=*),dimension(:),intent(in   ) :: tmpl_nusis
  integer         ,dimension(:),intent(in   ) :: tmpl_nuchn
  integer         ,dimension(:),intent(inout) :: tmpl_inuse
  integer,intent(in) :: ntmpl

  character(len=*),dimension(:),intent(in) :: xxxx_nm
  character(len=*),dimension(:),intent(in) :: xxxx_sn
  integer         ,dimension(:),intent(in) :: xxxx_ch
  integer,intent(in) :: nchan
  
  character(len=*),intent(in) :: tag
  logical,optional,intent(in) :: report
  logical,optional,intent(in) :: warning
  
  character(len=len(tmpl_nusis)) :: xxxx_sis
  integer :: xxxx_chn

  logical :: xxxx_matched,report_,warning_
  integer :: ichan,itmpl

  report_   =.false.; if(present(report)) report_=report
  warning_  =.false.; if(present(warning)) warning_=warning

  call tell(myname,'nchan =',nchan)
  call tell(myname,'ntmpl =',ntmpl)

  do ichan=1,nchan

        ! With a given (nm,ch) entry in the active-levels list,
        ! search the full-levels list for a match, to find the
        ! corresponding (sis,chn) value.

    xxxx_sis=trim(xxxx_sn(ichan))//'_'//trim(xxxx_nm(ichan))
    xxxx_chn=xxxx_ch(ichan)

        ! With the (si,ch) of the matching (nm,ch), locate the same
        ! entry in the ozinfo-template list, and set its _inuse flag.

    xxxx_matched=.false.
    do itmpl=1,ntmpl
      xxxx_matched = xxxx_sis==tmpl_nusis(itmpl) .and. &
                     xxxx_chn==tmpl_nuchn(itmpl)

      if( xxxx_matched ) then

        tmpl_inuse(itmpl)=inuse ! mark this channel to [inuse]
        xxxx_matched=.true.

	if(report_) call tell(myname,'flag "'//trim(tag)//'" is set to '// &
		'"'//trim(xxxx_sis)//'" chan#',xxxx_chn)
        exit
      endif
    enddo

    if((.not.xxxx_matched) .and. warning_) then
      call warn(myname,trim(tag)//'_sis =',xxxx_sis)
      call warn(myname,trim(tag)//'_chn =',xxxx_chn)
      call warn(myname,'['//trim(tag)//'_sis:chn] entry not seen in the ozinfo-template')
      cycle
    endif
  enddo

  call tell(myname,'exiting, tag =',tag)
end subroutine tmpl_xcheck
end program make_ozinfo
