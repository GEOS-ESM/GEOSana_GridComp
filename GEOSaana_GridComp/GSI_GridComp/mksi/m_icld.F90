module m_icld
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module m_icld
! abstract: update icld number
!
!$$$  end subprogram documentation block

! module interface:

  use satinfo_util, only: perr, die, tell

  implicit none
  private       ! except
  public:: icld_TRUE_   ! hard off/on switch
  public:: icld_yes     ! soft off/on switch

  public:: icld_init    ! initialize icldtable(:) through a namelist input
        interface icld_init   ; module procedure nmlinit_; end interface

  public:: icld_nmlread   ! read a namelist to initialize
        interface icld_nmlread; module procedure nmlread_; end interface

  public:: icld_update    ! reset %icld_det it is a matching sitmpl entry
        interface icld_update ; module procedure update_ ; end interface

  public:: icld_clean     ! clean icldtable to its initial state
        interface icld_clean  ; module procedure clean_  ; end interface

!USECASE: 
!       use m_icld,only: icld_yes,icld_nmlread,icld_update
!       call icld_nmlread([<nml_filename>])
!       if(icld_yes) then
!         do i=1,sitmpl_count
!           call icld_update(sitmpl(i)%icld_det, &
!                            sitmpl(i)%isis, &
!                            sitmpl(i)%ichan )
!         enddo
!       endif
!       call icld_clean()

  integer,parameter:: MSIZE_=500
  !logical,parameter:: ICLD_TRUE_=.false.    ! a convinient way to turn off this feature.
  logical,parameter:: ICLD_TRUE_=.true.     ! a convinient way to turn off this feature.
  !logical,parameter:: ICLD_DEBUG=.true.     ! a convinient way to turn off this feature.
  logical,parameter:: ICLD_DEBUG=.false.    ! a convinient way to turn off this feature.

  logical,save:: icld_yes     = .true.
  logical,save:: icld_verbose = .false.
  integer,save:: icld_count=-1

  type:: icldtuple
    integer          :: lymdh=-1   ! lower yyyymmddhh
    integer          :: uymdh=-1   ! upper yyyymmddhh
    character(len=20):: isis =".undef."
    integer          :: ichan=-1   ! channel number
    integer          :: icld =1    ! corresponding icld_det values
  end type icldtuple

  type(icldtuple),dimension(MSIZE_),save:: icldTable=icldtuple()
  character(len=*),parameter:: ICLD_NMLFILE_default="sidb/icldtable.nml"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='m_icld'

contains
subroutine nmlinit_(icld_nmlfile)
  implicit none
  character(len=*),optional,intent(in):: icld_nmlfile
  if(present(icld_nmlfile)) then
    call nmlread_(icld_nmlfile)
  else
    call nmlread_(icld_nmlfile_default)
  endif
end subroutine nmlinit_

subroutine clean_()
  implicit none
  icld_count=0
  icldtable(:)=icldtuple()
  icld_yes=.false.
end subroutine clean_

subroutine nmlread_(icld_nmlfile)
  implicit none
  character(len=*),intent(in):: icld_nmlfile

! &icldtable_NML
!   verbose      =.true.    ! verbose output if additional runtime information is needed
!   icldtable(:) =          ! -- lymdh, uymdh, isis, ichan, icld --
!       2010010100      2100010100      iasi_metop-a       1671   31
!       2010010100      2100010100      airs_aqua          1012   31
!       ...
! /

  character(len=*),parameter:: myname_=myname//"::nmlread_"
  integer:: ip,ios,n
  integer:: inunit

  logical:: verbose  =.false.

  namelist/icldtable_NML/verbose,icldtable

  verbose  =.false.
  icldtable(:)=icldtuple()

  if(.not.ICLD_TRUE_) return

  open(newunit=inunit,file=icld_nmlfile,form='formatted',status="old",iostat=ios)
        if(ios/=0) then
          call perr(myname_,'open("'//trim(icld_nmlfile)//'"), iostat =',ios)
          call  die(myname_)
        endif

  read(inunit,icldtable_NML,iostat=ios)
        if(ios/=0) then
          call perr(myname_,'read(nml=icldtable_NML), iostat =',ios)
          call perr(myname_,'                   icld_nmlfile =',trim(icld_nmlfile))
          call  die(myname_)
        endif

  close(inunit)

        ! trim icldtable(:) off all ".undef." entries
  icld_count=0
  n=0
  do ip=1,size(icldtable)
    if(icldtable(ip)%isis==".undef." ) cycle

    n=n+1
    ! Copy all defined entries
    if(n/=ip) icldtable(n)=icldtable(ip)
  enddo

  icld_count  = n
  icld_yes    = icld_count>0
  icld_verbose= verbose

  call tell(myname_,'/icldtable_nml/ loaded, nmlfile =',icld_nmlfile)
  if(ICLD_DEBUG.or.icld_verbose) then
    call tell(myname_,'icld_yes     =',icld_yes)
    call tell(myname_,'icld_count   =',icld_count)
    call tell(myname_,'icld_verbose =',icld_verbose)
  endif
end subroutine nmlread_

subroutine update_(nymdh,isis,ichan,icld)
! Where (nymdh, isis, and ichan) keys match an entry in the table, update
! a scale factor to the given varx(:) argument.

  implicit none

  integer,intent(in):: nymdh    ! ((yr*100+mo)*100+dy)+hr
  character(len=*),intent(in   ):: isis
  integer         ,intent(in   ):: ichan
  integer         ,intent(inout):: icld

  character(len=*),parameter:: myname_=myname//"::update_"
  integer:: ip,im
  logical:: matched_

  if(.not.icld_TRUE_) return
  if(.not.icld_yes  ) return

  if(ICLD_DEBUG.and.isis=="iasi_metop-a".and.ichan==1427) then
    call tell(myname_,'given nymdh =',nymdh)
    call tell(myname_,'       isis =',isis)
    call tell(myname_,'      ichan =',ichan)
    call tell(myname_,'       icld =',icld)
  endif

  im=0
  do ip=1,icld_count
    matched_= icldtable(ip)%lymdh<=nymdh .and. nymdh<=icldtable(ip)%uymdh
    if(matched_) matched_= icldtable(ip)%isis ==isis
    if(matched_) matched_= icldtable(ip)%ichan==0 .or. icldtable(ip)%ichan==ichan
    if(ICLD_DEBUG.and.isis=="iasi_metop-a".and.ichan==1427) then
      call tell(myname_,'entry matched =',merge("YES","NOP",matched_))
      call tell(myname_,'        entry #',ip)
      call tell(myname_,'       [ymdh] =',[icldtable(ip)%lymdh,icldtable(ip)%uymdh])
      call tell(myname_,'       %isis  =',trim(icldtable(ip)%isis ))
      call tell(myname_,'       %ichan =',     icldtable(ip)%ichan )
      call tell(myname_,'       %icld  =',     icldtable(ip)%icld  )
    endif

    if(.not.matched_) cycle

    select case(im)
      case(0)
        im=ip

        if(icld_verbose) then
          call tell(myname_,'new icld applied to "'//trim(isis)// &
            '", [nymhh,ichan,icld->icld] =',[nymdh,ichan,icld,icldtable(ip)%icld])
        endif
        icld=icldtable(ip)%icld ! reset icld value according to the database

      case default
        call perr(myname_,'redendent icldtable entries, isis =',trim(isis))
        call perr(myname_,'                            ichan =',ichan)
        call perr(myname_,'                            nymdh =',nymdh)
        call perr(myname_,'                       icld_count =',icld_count)
        call perr(myname_,'                   previous entry #',im)
        call perr(myname_,'            icldtable(prev)%isis  =',trim(icldtable(im)%isis ))
        call perr(myname_,'            icldtable(prev)%ichan =',     icldtable(im)%ichan )
        call perr(myname_,'            icldtable(prev)%lymdh =',     icldtable(im)%lymdh )
        call perr(myname_,'            icldtable(prev)%uymdh =',     icldtable(im)%uymdh )
        call perr(myname_,'            icldtable(prev)%icld  =',     icldtable(im)%icld  )
        call perr(myname_,'                    current entry #',ip)
        call perr(myname_,'            icldtable(curr)%isis  =',trim(icldtable(ip)%isis ))
        call perr(myname_,'            icldtable(curr)%ichan =',     icldtable(ip)%ichan )
        call perr(myname_,'            icldtable(curr)%lymdh =',     icldtable(ip)%lymdh )
        call perr(myname_,'            icldtable(curr)%uymdh =',     icldtable(ip)%uymdh )
        call perr(myname_,'            icldtable(curr)%icld  =',     icldtable(ip)%icld  )
        call  die(myname_)
    endselect
  enddo
end subroutine update_
end module m_icld
