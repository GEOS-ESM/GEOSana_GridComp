!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_sitmpl - satinfo si template file IO
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_sitmpl
      implicit none
      private   ! except

      public :: sitmpl_get
      public :: sitmpl_put
      public :: sitmpl_clean

    interface sitmpl_get; module procedure get_; end interface
    interface sitmpl_put; module procedure put_; end interface
    interface sitmpl_clean; module procedure clean_; end interface

! !REVISION HISTORY:
!       08Jun07 - Jing Guo <guo@gmao.gsfc.nasa.gov>
!               - initial prototype/prolog/code
!       10Jul12 - Jing Guo <jguo@nasa.gov>
!               . fixed incorrect name conversion to n05 (NOAA-5) back
!                 to tirosn (TIROS-N).
!       22Sep16 - Jing Guo <jing.guo@nasa.gov>
!               . Used token_shift() to keep unknown tails in entries.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_sitmpl'

contains
subroutine get_(fname,jpch,nusis,nusat,nuchan,iuse_rad,rest,vern,append)
!--  Read in satinfo template file here
  use satinfo_util,only : alloc,realloc
  use satinfo_util,only : luavail,stdout
  use satinfo_util,only : die,perr,tell
  use satinfo_util,only : token_shift
  implicit none
  character(len=*),intent(in) :: fname  ! input satinfo.rc.tmpl
        ! satinfo table
  integer,intent(inout) :: jpch             ! entry count
  character(len=*),pointer,dimension(:) :: nusis ! Sensor-Instr.-Satil.
  integer,pointer,dimension(:) :: nusat     ! old sensor_sat key
  integer,pointer,dimension(:) :: nuchan    ! satellite channel
  integer,pointer,dimension(:) :: iuse_rad  ! channel in-use flag
  character(len=*),pointer,dimension(:):: rest  ! the full input lines
  integer,intent(out) :: vern   ! version of the input template file format
  logical,optional,intent(in) :: append
  
  character(len=*),parameter :: myname_=myname//"::get_"
  character(len=*),parameter :: proc_=myname_//'(): '

  integer :: lu,ios,j,ir
  character(len=512) :: line
  character(len=1) :: key
  logical :: append_,skipline
  append_=.false.; if(present(append)) append_=append

  vern=-1
!! Open ...
  lu=luavail()
  open(lu,file=fname,form='formatted',status='old',iostat=ios)
    if (ios /= 0) call die(myname_, &
      'open("'//trim(fname)//'") for input, iostat =',ios)

  call tell(myname_,'reading "'//trim(fname)//'" for input')

  if(.not.append_) then !! note jpch is set to 0 by alloc()
    call alloc(nusis    ,jpch)
    call alloc(nusat    ,jpch)
    call alloc(nuchan   ,jpch)
    call alloc(iuse_rad ,jpch)
    call alloc(rest     ,jpch)
  endif
  j=jpch

           nusis(j+1:)='.undef.'
           nusat(j+1:)=HUGE(j)
          nuchan(j+1:)=HUGE(j)
        iuse_rad(j+1:)=HUGE(j)
            rest(j+1:)=""

        ir=0
        read(lu,'(a)',iostat=ios) line
        do while(ios==0)
          ir=ir+1

          call realloc(nusis    ,j,incr=1)
          call realloc(nusat    ,j,incr=1)
          call realloc(nuchan   ,j,incr=1)
          call realloc(iuse_rad ,j,incr=1)
          call realloc(rest     ,j,incr=1)

          key='?'
          skipline=.false.
          read(line,*,iostat=ios) key
          select case(key(1:1))
          case('?')
            call die(myname_,'cann''t read, rec # =',ir)

          case('!')
            skipline=.true.
            if(vern==-1) vern=3

          case default          ! in the "new" format
            if(vern==-1) vern=2
            j=j+1
            read(line,*,iostat=ios) nusis(j),nuchan(j),iuse_rad(j)

                ! If this record is an expected end-of-table mark,
                ! end the read loop.
            if( ios/=0 .or. nusis(j)=='sensor'.or.nusis(j)=='sat' ) exit
            nusat(j)=-1
            rest(j)=token_shift(line,3) ! 3 tokens, for nusis, nuchan, and iuse_rad
          end select

                ! If the input line contains unexpected values,
                ! which is often the result of a failed 
                ! "list-directed" read (i.e. fmt=*), echo the
                ! input then die().
          if(.not.skipline) then
                ! Count in the record as a good table entry, and
                ! set the use_rad flag to _off_.

            jpch=j
            iuse_rad(j)=-1
          endif

                ! Read the next record
          read(lu,'(a)',iostat=ios) line
        enddo
        close(lu)

        call tell(myname_,'number of channels, jpch =',jpch)
        call tell(myname_,'input format version, vern =',vern)

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

    subroutine put_(fname,jpch,nusis,nusat,nuchan,iuse_rad,rest,vern,nymd,nhms)
      use satinfo_util,only : luavail
      use satinfo_util,only : die,tell
      implicit none
      character(len=*),intent(in) :: fname      ! output satinfo.txt
                                                ! satinfo table
      integer,intent(in) :: jpch                ! entry count
      character(len=*),dimension(:),intent(in) :: nusis ! Sensor-Instr.-Satil.
      integer,dimension(:),intent(in) :: nusat     ! old sensor_sat key
      integer,dimension(:),intent(in) :: nuchan    ! satellite channel
      integer,dimension(:),intent(in) :: iuse_rad  ! channel in-use flag
      character(len=*),dimension(:),intent(in) :: rest
      integer,intent(in) :: vern ! vertion of the output template file format
      integer,intent(in) :: nymd,nhms

! !REVISION HISTORY:
!       08Jun07 - Jing Guo <guo@gmao.gsfc.nasa.gov>
!               - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::put_'
  integer :: lu,ios,j,l
  integer :: nusatj
  character(len=max(20,len(nusis))) :: nusisj

        ! Modify original satinfo file according to the satinfo_base.rc
        !===============================================================
        
  call tell(myname_,'writing "'//trim(fname)//'" for output')

        lu=luavail()
        open(lu,file=fname,form='formatted',status='unknown',iostat=ios)
          if(ios/=0) call die(myname_, &
            'open("'//trim(fname)//'") for output, iostat =',ios)

  call tell(myname_,'output format version, vern =',vern)

        select case(vern)
        case(3)
          write(lu,'(a,i5,a,i8.8,1x,i6.6,a)') &
            '!sensor/instr/sat   chan iuse   error error_cld ermax   var_b  var_pg  cld_det  ### ',&
                    jpch, ' total (',nymd,nhms,')'
        endselect
        do j = 1,jpch
          nusisj=nusis(j)
          nusatj=nusat(j)

          select case(vern)
          case(1)
            call die(myname_,'no longer supported, vern = 1')

          case(2:)
            l=max(18,len_trim(nusisj))
            write(lu,110,iostat=ios) nusisj(1:l),nuchan(j),iuse_rad(j),trim(rest(j))
            if(ios/=0) call die(myname_, &
                'write-110("'//trim(fname)//'") for output, iostat =',ios)
          endselect
          
        enddo
        select case(vern)
        case(1)
          call die(myname_,'no longer supported, vern = 1')
        case(2)
          write(lu,'(a,i5,a,i8.8,1x,i6.6,a)') &
            '!sensor/instr/sat   chan iuse  error/error_cld ermax  var_b  var_pg # ', &
                    jpch, ' total (',nymd,nhms,')'
        endselect
        close(lu)

110     format(1x,a,i5,i5,a)
115     format(i6,i5,i5,f7.3,f7.3,f7.3,f7.2,f8.3) 
111     format(' Channel Index :',/,30(10i6/)) 
112     format(' Channel Number:',/,30(10i6/)) 
113     format(' Sat/Inst ID   :',/,30(10i6/)) 
end subroutine put_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean pointers
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(jpch,nusis,nusat,nuchan,iuse_rad,rest)
      use satinfo_util,only : dealloc
      implicit none
      integer,intent(out) :: jpch                   ! entry count
      character(len=*),pointer,dimension(:) :: nusis ! Sensor-Instr.-Satil.
      integer,pointer,dimension(:) :: nusat     ! old sensor_sat key
      integer,pointer,dimension(:) :: nuchan    ! satellite channel
      integer,pointer,dimension(:) :: iuse_rad  ! channel in-use flag
      character(len=*),pointer,dimension(:) :: rest

! !REVISION HISTORY:
!       08Jun07 - Jing Guo <guo@gmao.gsfc.nasa.gov>
!               - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: j

  j=jpch
  call dealloc(nusis    ,j)
  call dealloc(nusat    ,j)
  call dealloc(nuchan   ,j)
  call dealloc(iuse_rad ,j)
  call dealloc(rest     ,j)
  jpch=-1

end subroutine clean_
subroutine sat2sis_(nusat,sisname)
!-- map old numerical nusat values to new sis values.
  implicit none
  integer,intent(in) :: nusat
  character(len=*),intent(out) :: sisname

  integer :: isat

select case(nusat)
case(5)
  isat=nusat
  sisname='hirs2_tirosn'
case(6:12,14)
  isat=nusat
  write(sisname,'(a,i2.2)') 'hirs2_n',isat
  !nusat=-1
case(15:17)
  isat=nusat
  write(sisname,'(a,i2.2)') 'hirs3_n',isat
  !nusat=-1
case(18)
  isat=nusat
  write(sisname,'(a,i2.2)') 'hirs4_n',isat
  !nusat=-1
case(49)
  isat=nusat
  sisname='airs281SUBSET_aqua'
  !nusat=-1
case(58,60,62)
  isat=nusat-50
  write(sisname,'(a,i2.2)') 'sndr_g',isat
  !nusat=-1
case(205)
  isat=nusat-200
  sisname='msu_tirosn'
case(206:212,214)
  isat=nusat-200
  write(sisname,'(a,i2.2)') 'msu_n',isat
  !nusat=-1
case(258,260,262)
  isat=nusat-250
  write(sisname,'(a,i2.2)') 'imgr_g',isat
  !nusat=-1
case(305)
  isat=nusat-300
  sisname='ssu_tirosn'
case(306:312,314)
  isat=nusat-300
  write(sisname,'(a,i2.2)') 'ssu_n',isat
  !nusat=-1
case(315:317)
  isat=nusat-300
  write(sisname,'(a,i2.2)') 'amsua_n',isat
  !nusat=-1
case(318)
  isat=nusat-300
  write(sisname,'(a,i2.2)') 'amsua_n',isat
  !nusat=-1
case(349)
  isat=nusat-300
  sisname='amsua_aqua'
  !nusat=-1
case(415:417)
  isat=nusat-400
  write(sisname,'(a,i2.2)') 'amsub_n',isat
  !nusat=-1
case(418)
  isat=nusat-400
  write(sisname,'(a,i2.2)') 'mhs_n',isat
  !nusat=-1
case(449)
  isat=nusat-400
  sisname='hsb_aqua'
  !nusat=-1
case(516)
  isat=nusat-500
  write(sisname,'(a,i2.2)') 'ssmis_f16'
  !nusat=-1
case(549)
  isat=nusat-500
  sisname='amsre_aqua'
  !nusat=-1
case(616:618)
  isat=nusat-600
  write(sisname,'(a,i2.2)') 'avhrr3_n',isat
  !nusat=-1
case(708,710,711,713:715)
  isat=nusat-700
  write(sisname,'(a,i2.2)') 'ssmi_f',isat
  !nusat=-1
case default
  isat=nusat
  write(sisname,'(a,i3.3)') 'unknown_',isat
  !nusat=-1
end select
end subroutine sat2sis_
end module m_sitmpl
