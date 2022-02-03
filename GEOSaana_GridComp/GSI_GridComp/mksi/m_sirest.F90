module m_sirest
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module m_sirest
!   prgmmr:      j guo <jing.guo-1@nasa.gov>
!      org:      NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:      2021-01-26
!
! abstract: for remaining fields of sitmpl
!
! program history log:
!   2021-01-26  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! module interface:

  use satinfo_util, only: perr,die
  use satinfo_util, only: token_shift
  implicit none
  private       ! except
  public:: sirest_update
        interface sirest_update; module procedure update_; end interface

!> In a satinfo file, ...      |-- sirest starts here and may be writen as "5(f9.3),i8"
!>
!>!sensor/instr/sat   chan iuse   error error_cld ermax  var_b   var_pg  icld_det
!> amsua_n15             1    1   3.000   9.100   4.500  10.000   0.000      1
!> amsua_n15             2    1   2.000  13.500   4.500  10.000   0.000      1
!> ...
  
  integer,parameter:: BUFRSIZE=256
  character(len=*),parameter:: FMT_sirest="(5(f8.3),i7,a)"
  type:: sirest
    private
    real   :: varch     =0.
    real   :: varch_cld =0.
    real   :: ermax_rad =0.
    real   :: b_rad     =0.
    real   :: pg_rad    =0.
    integer:: icld      =-1
  end type sirest

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='m_sirest'
  logical,parameter:: IUSE_ALL=.false.

contains
subroutine update_(jpch,nusis,nuchan,iuse,rest,vern,nymd,nhms)
! abstract: update additional fields in rest

  use m_icld, only: icld_yes,icld_update
  implicit none
  integer         ,intent(in):: jpch                   ! entry count
  character(len=*),dimension(:),intent(in):: nusis     ! Sensor-Instr.-Satil.
  integer         ,dimension(:),intent(in):: nuchan    ! satellite channel
  integer         ,dimension(:),intent(in):: iuse      ! channel in-use flag
  character(len=*),dimension(:),intent(inout):: rest   ! additonal fields
  integer         ,intent(in):: vern            ! expected version of the output format
  integer         ,intent(in):: nymd,nhms       ! date-time tag: (yyyymmdd,hhmms)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//'::update_'

  type(sirest):: sir
  character(len=BUFRSIZE):: tail_
  integer:: ios,j,nymdh

  if(.not.icld_yes) return

  select case(vern)
  case(3:)
    continue
  case default
    call die(myname_,'not supported, vern =',vern)
  end select

  ios=0
  do j=1,jpch
    read(rest(j),*,iostat=ios) sir
        if(ios/=0) then
          call perr(myname_,'read(rest(j)), iostat =',ios)
          call perr(myname_,'              rest(j) =',rest(j))
          call perr(myname_,'                    j =',j)
          call  die(myname_)
        endif
    tail_=token_shift(rest(j),6)

    nymdh=nymd*100+nhms/10000
    if(IUSE_ALL.or.iuse(j)>=0) call icld_update(nymdh,nusis(j),nuchan(j),sir%icld)
    write(rest(j),fmt=FMT_sirest) sir,trim(tail_)
  enddo

end subroutine update_
end module m_sirest
!.
