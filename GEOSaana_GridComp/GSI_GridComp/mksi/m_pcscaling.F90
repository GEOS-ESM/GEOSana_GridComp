module m_pcscaling
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module m_pcscaling
! abstract: apply scale factors to varx of predictor coefficients
!
!$$$  end subprogram documentation block

! module interface:

  use satinfo_util, only: DP
  use satinfo_util, only: perr, die, tell

  implicit none
  private       ! except
  public :: pcscaling_yes       ! data structure
  public :: pcscaling_nmlread   ! read a namelist to initialize
    interface pcscaling_nmlread; module procedure nmlread_; end interface
  public :: pcscaling_apply     ! apply a table of scale-factors to varx(:)
    interface pcscaling_apply; module procedure apply_; end interface

  integer,parameter:: MSIZE_=100
! integer,parameter:: MSIZE_=4
  !logical,parameter:: PCSCALING_TRUE_=.false.    ! a convinient way to turn off this feature.
  logical,parameter:: PCSCALING_TRUE_=.true.     ! a convinient way to turn off this feature.

  logical,save:: pcscaling_yes = .false.
  logical,save:: pcscaling_verbose = .false.
  integer,save:: pcscaling_count=-1

  character(len=20),dimension(MSIZE_),save:: pcscaling_isis  = ".undefined."
  integer          ,dimension(MSIZE_),save:: pcscaling_ichan = -1
  integer          ,dimension(MSIZE_),save:: pcscaling_nymdh = -1
  real(kind=DP)    ,dimension(MSIZE_),save:: pcscaling_scale = 1._DP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='m_pcscaling'

contains
subroutine nmlread_(pcscaling_nmlfile)
  implicit none
  character(len=*),intent(in):: pcscaling_nmlfile

! &PCSCALING_NML
!   pcscaling=.true.    ! turn pcscaling on/off (default is off)
!   verbose  =.true.    ! in case there is a need to show additional information
!   pcscales =          ! -- nymdh, isis, ichan, scale --
!       2020020200  hirs3_n15  14  1.6
!       2020122200  hirs3_n15  14 .625
!       ...
! /

  character(len=*),parameter:: myname_=myname//"::nmlread_"
  integer:: ip,ios,n
  integer:: inunit

  type:: pcscale
    integer          :: nymdh
    character(len=20):: isis
    integer          :: ichan
    real(kind=DP)    :: scale
  end type pcscale
  type(pcscale):: pcscale_default = pcscale(-1,".undefined.",-1,1._DP)

  logical:: pcscaling=.false.
  logical:: verbose  =.false.
  type(pcscale),dimension(MSIZE_):: pcscales

  namelist/pcscaling_NML/pcscaling,verbose,pcscales

  pcscaling=.false.
  verbose  =.false.
  pcscales(:)=pcscale_default

  if(.not.PCSCALING_TRUE_) return

  open(newunit=inunit,file=pcscaling_nmlfile,form='formatted',status="old",iostat=ios)
  if(ios==0) then
    read(inunit,pcscaling_NML,iostat=ios)
        if(ios/=0) then
          pcscaling=.false.
          call perr(myname_,'NAMELIST/PCSCALING/ reading error, ios =',ios)
          call  die(myname_)
        endif
  endif
  close(inunit)

  pcscaling_yes=pcscaling
  pcscaling_count=0

  n=0
  do ip=1,size(pcscales)
    if(pcscales(ip)%isis==".undefined." ) cycle

    n=n+1
    ! Copy all defined entries
    pcscaling_isis (n)=pcscales(ip)%isis
    pcscaling_ichan(n)=pcscales(ip)%ichan
    pcscaling_nymdh(n)=pcscales(ip)%nymdh
    pcscaling_scale(n)=pcscales(ip)%scale
  enddo

  pcscaling_count=n
  pcscaling_yes=pcscaling_yes.and.pcscaling_count>0

  call tell(myname_,'  pcscaling_yes =',pcscaling_yes)
  call tell(myname_,'pcscaling_count =',pcscaling_count)
  !do ip=1,pcscaling_count
  !  write(*,'(i4,i12,2x,a12,i6,f12.4)') ip,&
  !    pcscaling_nymdh(ip),pcscaling_isis(ip),pcscaling_ichan(ip),pcscaling_scale(ip)
  !enddo
end subroutine nmlread_

subroutine apply_(nymdh,isis,ichan,varx,changed)
! Where (nymdh, isis, and ichan) keys match an entry in the table, apply
! a scale factor to the given varx(:) argument.

  implicit none

  integer,intent(in):: nymdh    ! ((yr*100+mo)*100+dy)+hr
  character(len=*),intent(in):: isis
  integer,intent(in):: ichan
  real(kind=DP),dimension(:),intent(inout):: varx(:)
  logical,optional,intent(inout):: changed

  character(len=*),parameter:: myname_=myname//"::apply_"
  integer:: ip,ionce_
  logical:: matched_

  if(.not.PCSCALING_TRUE_) return
  if(.not.pcscaling_yes) return

  ionce_=0
  do ip=1,pcscaling_count
    matched_= pcscaling_nymdh(ip)==nymdh
    if(matched_) matched_= pcscaling_isis(ip)==isis
    if(matched_) matched_= pcscaling_ichan(ip)==0 .or. pcscaling_ichan(ip)==ichan
    !if(isis=='hirs3_n15') then
    !  write(*,'(l2,i2,2i12,2x,2a10,2i4)') matched_,ip,nymdh,pcscaling_nymdh(ip),isis,pcscaling_isis(ip),ichan,pcscaling_ichan(ip)
    !endif

    if( matched_ ) then
      select case(ionce_)
      case(0)
        varx(:)=pcscaling_scale(ip)*varx(:)
        if(present(changed)) changed=.true.
        ionce_=ip

        if(pcscaling_verbose) then
          call tell(myname_,'pcscaling applied to isis =',trim(isis))
          call tell(myname_,'                    ichan =',ichan)
          call tell(myname_,'                    nymdh =',nymdh)
          call tell(myname_,'                    scale =',pcscaling_scale(ip))
        endif

      case default
        call perr(myname_,'redendent pcscaling entries, pcscaling_count =',pcscaling_count)
        call perr(myname_,'                              previous entry #',ionce_)
        call perr(myname_,'                        pcscaling_isis(prev) =',trim(pcscaling_isis (ionce_)))
        call perr(myname_,'                       pcscaling_ichan(prev) =',     pcscaling_ichan(ionce_) )
        call perr(myname_,'                       pcscaling_nymdh(prev) =',     pcscaling_nymdh(ionce_) )
        call perr(myname_,'                       pcscaling_scale(prev) =',     pcscaling_scale(ionce_) )
        call perr(myname_,'                               current entry #',ip)
        call perr(myname_,'                        pcscaling_isis(curr) =',trim(pcscaling_isis (ip)))
        call perr(myname_,'                       pcscaling_ichan(curr) =',     pcscaling_ichan(ip) )
        call perr(myname_,'                       pcscaling_nymdh(curr) =',     pcscaling_nymdh(ip) )
        call perr(myname_,'                       pcscaling_scale(curr) =',     pcscaling_scale(ip) )
        call  die(myname_)
      end select

    endif
  enddo
end subroutine apply_
end module m_pcscaling
