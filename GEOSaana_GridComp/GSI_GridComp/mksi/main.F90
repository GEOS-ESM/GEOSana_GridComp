!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: make_satinfo - create a satinfo file for a given date-time
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"
#include "assert.H"

    program make_satinfo
      use m_actvchan  , only : actvchan_select,actvchan_clean
      use m_sitmpl    , only : sitmpl_get,sitmpl_put,sitmpl_clean
      use satinfo_util, only : tell,warn,perr,die,stdin,stdout
      use satinfo_util, only : alloc,realloc,dealloc
      implicit none

! !REVISION HISTORY:
!	01Sep10	- Jing Guo <jing.guo@nasa.gov>
!		- Changed to process additional channel tables, for 
!		  observation-bias-correction classes, 2, 3, and 4.
!	06Jul09	- Jing Guo <jing.guo@nasa.gov>
!		- changed code comments
!	20Apr09	- Jing Guo <jing.guo@nasa.gov>
!		- removed "satops.tbl" and "satinst.tbl"
!		- renamed "usable_channels.tbl" to "available_channels.tbl"
!		- renamed usableXXXX related variables to availableXXXX
!		- added "monitor.tbl" and related variables
!		- implemented a new algorithm flowing through
!		  1) avaialble-channels table,
!		  2) active-channels table, then, optionally,
!		  3) monitor-channels table.
! 	21Sep07	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- added "usable_channels" table for "passive" data.
!		- added user-definable values of INUSE flag. (setup.nml)
! 	24May07	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- redesigned with a more explicit database concept
!		- added additional tables to the existing database
!		- reengineered with modules
!		- removed dependency to GMAO libraries
!       04Jan07 - Jing Guo <guo@gmao.gsfc.nasa.gov>
!               - restructured the original &SETUP namelist to getjpch_
!		  to avoid crashing on every variable added to the GSI
!		  &SETUP namelist.
!		- modified exit messages to better inform users.
!	14Feb07 - R. Todling
!		- add usage routine
! 	before  - E. Liu
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='make_satinfo'
  integer,parameter :: PATHSIZE=512
  integer,parameter :: NAMESIZE=20
  integer,parameter :: BUFRSIZE=256

  character(len=NAMESIZE),pointer,dimension(:) :: xxxx_nm => null()
  character(len=NAMESIZE),pointer,dimension(:) :: xxxx_sn => null()
  integer,pointer,dimension(:) :: xxxx_ch => null()
  integer :: nchan=-1

  character(len=NAMESIZE),pointer,dimension(:) :: tmpl_nusis => null()
  integer,pointer,dimension(:) :: tmpl_nusat => null()
  integer,pointer,dimension(:) :: tmpl_nuchn => null()
  integer,pointer,dimension(:) :: tmpl_inuse => null()
  character(len=BUFRSIZE),pointer,dimension(:) :: tmpl_rest => null()
  integer :: ntmpl=-1

    ! key parameters.  Must be defined explicitly.
  integer :: nymd
  integer :: nhms

  integer,parameter :: default_NOTUSED   = -9
  integer,parameter :: default_AVAILABLE = -1
  integer,parameter :: default_MONITOR   = +0
  integer,parameter :: default_ACTIVE    = +1

  character(len=*),parameter :: default_AVAILABLECHAN  = 'available_channels.tbl'
  character(len=*),parameter :: default_ACTIVECHAN     =    'active_channels.tbl'
  character(len=*),parameter :: default_MONITORCHAN    =   'monitor_channels.tbl'

    ! optional parameters.  Default values can be used.
  character(len=PATHSIZE) ::    info_tmpl = "satinfo.tmpl"
  character(len=PATHSIZE) ::    info_outf = "satinfo.txt"
  character(len=PATHSIZE) :: satinfo_tmpl = "satinfo.tmpl"	! backward compatable
  character(len=PATHSIZE) :: satinfo_outf = "satinfo.txt"	! backward compatable
  character(len=PATHSIZE) :: dbname        ='./sidb'
  character(len=PATHSIZE) :: availablechan_tbl = default_AVAILABLECHAN
  character(len=PATHSIZE) ::    activechan_tbl = default_ACTIVECHAN
  character(len=PATHSIZE) ::   monitorChan_tbl = default_MONITORCHAN

  logical :: verbose=.false.
  logical :: nowarn =.false.
  logical :: samever=.false.

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
  character(len=PATHSIZE) :: chan_tbl

  namelist/setup/ nymd,nhms,dbname,	&
       info_tmpl,   info_outf,	&
    satinfo_tmpl,satinfo_outf,	&
    availablechan_tbl,	&
    activechan_tbl,	&
    monitorchan_tbl,	&
    verbose,samever,	&
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

  	! if the values between variables with the old name and the new name do
	! not agree, make both to the same new value.
  if(satinfo_tmpl/=info_tmpl .and. satinfo_tmpl/="satinfo.tmpl") info_tmpl=satinfo_tmpl
  if(satinfo_outf/=info_outf .and. satinfo_outf/="satinfo.txt" ) info_outf=satinfo_outf

  call tell(myname,'processing for [nymd, nhms] = ',(/nymd,nhms/),format='(i9.6)')

!!$ Mark tmpl-channels list according to the active-channels list.
!!$   for the same sisname-sischan, set rad_flag = yes

  call sitmpl_get(satinfo_tmpl,ntmpl, &
    tmpl_nusis,tmpl_nusat,tmpl_nuchn,tmpl_inuse,tmpl_rest, vern=iver )

!!$ <<< Initialization pass >>> All channels in the template are turned off
  where(tmpl_inuse(1:ntmpl)>NOTUSED_FLAG)
    tmpl_inuse(1:ntmpl)=NOTUSED_FLAG
  endwhere

!!$ <<< Pass 0:n >>> All channels are reset if matched by a table entry

  ! iuse_flag values and priorities (latter override earliers)
  ! = -2 do not use                                             everything in satinfo.tmpl
  ! = -1 monitor if diagnostics produced                        available channels
  ! =  1 use data with complete quality control                 active channels
  ! =  2 use data with no airmass bias correction               iuseClass-2 channels
  ! =  3 use data with no angle dependent bias correction       iuseClass-3 channels
  ! =  4 use data with no bias correction                       iuseClass-4 channels
  ! =  0 monitor and use in QC only                             monitoring channels

  do itbl=0,5
    report=.false.
    isOptional=.false.
    select case(itbl)
    case(0)
	!!$ <<< Pass 1 >>> Set all channels PASSIVE
      iuse_tbl=AVAILABLE_FLAG
      mtag_tbl="passive"
      chan_tbl=availablechan_tbl

    case(1)
      iuse_tbl=ACTIVE_FLAG
      mtag_tbl="active"
      chan_tbl=activechan_tbl

    case(2)
      iuse_tbl=2
      report=.true.
      isOptional=.true.
      mtag_tbl="obcClass-2"
      chan_tbl="iuseClass-2_channels.tbl"
    case(3)
      iuse_tbl=3
      report=.true.
      isOptional=.true.
      mtag_tbl="obcClass-3"
      chan_tbl="iuseClass-3_channels.tbl"
    case(4)
      iuse_tbl=4
      report=.true.
      isOptional=.true.
      mtag_tbl="obcClass-4"
      chan_tbl="iuseClass-4_channels.tbl"

    case(5)
      iuse_tbl=MONITOR_FLAG
      report=.true.
      isOptional=.true.
      mtag_tbl='monitor'
      chan_tbl=monitorchan_tbl

    case default
      call die(myname,'unknown table file, #',itbl)
    end select

    if( isOptional .and. &
    	.not. file_exist_(chan_tbl,dbase=dbname)) cycle
    
	!!$ Read the channel list
    call actvchan_select(xxxx_nm,xxxx_sn,xxxx_ch,nchan, &
      where_dt=(/nymd,nhms/),	&
      from=chan_tbl,dbase=dbname,verbose=.true.)

	!!$ Set all channels to iuse_tbl
    call tmpl_xcheck(iuse_tbl,tmpl_nusis,tmpl_nuchn,tmpl_inuse,ntmpl, &
      xxxx_nm,xxxx_sn,xxxx_ch,nchan, tag=mtag_tbl, &
      warning=warning,report=report)

	!!$ reset working arrays
    call actvchan_clean(xxxx_nm,xxxx_sn,xxxx_ch,nchan)
  enddo

!!$ write out the satinfo table.
  if(.not.samever) iver=max(3,iver) ! use the later version by the default

  call sitmpl_put(satinfo_outf,ntmpl, &
    tmpl_nusis,tmpl_nusat,tmpl_nuchn,tmpl_inuse,tmpl_rest, &
    vern=iver,nymd=nymd,nhms=nhms)

  call sitmpl_clean(ntmpl,tmpl_nusis,tmpl_nusat,tmpl_nuchn,tmpl_inuse,tmpl_rest)

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

  do ichan=1,nchan

        ! With a given (nm,ch) entry in the active-channels list,
        ! search the full-channels list for a match, to find the
        ! corresponding (sis,chn) value.

    xxxx_sis=trim(xxxx_sn(ichan))//'_'//trim(xxxx_nm(ichan))
    xxxx_chn=xxxx_ch(ichan)

        ! With the (si,ch) of the matching (nm,ch), locate the same
        ! entry in the satinfo-template list, and set its _inuse flag.

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
      call warn(myname,'['//trim(tag)//'_sis:chn] entry not seen in the satinfo-template')
      cycle
    endif
  enddo

end subroutine tmpl_xcheck
end program make_satinfo
