module m_fsi_weight
use kinds, only: i_kind,r_kind
use constants, only: one
use mpimod, only: mype
use mpeu_util, only: gettablesize
use mpeu_util, only: gettable
use mpeu_util, only: die,tell
use mpeu_util, only: luavail
implicit none
private

public :: fsi_init
public :: fsi_weight
public :: fsi_apply_weight

interface fsi_init
  module procedure init_
end interface
interface fsi_apply_weight
  module procedure apply_lev_weight_
  module procedure apply_chn_weight_
end interface

integer, parameter :: MAXSTR=256
character(len=*),parameter :: myname = 'm_fsi_weight'
logical :: fsi_weight
logical :: iamroot_

type fsi_type
  integer :: ncri ! number of criteria
  character(len= 4),allocatable :: criterium(:)
  character(len=20),allocatable :: variable(:)
  integer,     allocatable :: instrument(:)
  real(r_kind),allocatable :: lata(:), latb(:)
  real(r_kind),allocatable :: levl(:), levh(:)
  real(r_kind),allocatable :: sfactor(:)
end type

logical :: initialized_=.false.
type(fsi_type) :: fsi

contains

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  init_ --- Initialize info about FSI weighting factors
!
! !INTERFACE:
!
subroutine init_ (iamroot)
! !USES:
use mpeu_util, only: die
implicit none
! !INPUT PARAMETERS:
   logical,optional,intent(in) :: iamroot 
! !DESCRIPTION: Define scale factors based on offline evaluation
!               of FSI.
!
! !REVISION HISTORY:
!   2025-02-24  todling  initial code
!
! !REMARKS:
!   language: f90
!   machine:  discover
!
! !AUTHOR:
!   Ricardo Todling  org: gmao      date: 2025-02-24
!
!EOP
!-------------------------------------------------------------------------
!BOC
character(len=*),parameter:: rcname='anavinfo'  ! filename should have extension
character(len=*),parameter:: tbname='fsi_weights::'
integer :: luin,ii,ntot,nrows,instrument
character(len=MAXSTR),allocatable,dimension(:):: utable
character(len=4) criterium
character(len=20) var
real(r_kind) lata,latb
real(r_kind) levl,levh
real(r_kind) sfactor
character(len=*),parameter::myname_=myname//'*init_'

fsi_weight = .false.
if(initialized_) return

iamroot_=mype==0
if(present(iamroot)) iamroot_=iamroot 

! load file
luin=luavail()
open(luin,file=rcname,form='formatted')

! Scan file for desired table first
! and get size of table
call gettablesize(tbname,luin,ntot,nrows)
if(nrows==0) then
   close(luin)
   return
endif
fsi%ncri=nrows

! Get contents of table
allocate(utable(fsi%ncri))
call gettable(tbname,luin,ntot,fsi%ncri,utable)

! release file unit
close(luin)

! Retrieve each token of interest from table and define
! variables participating in state vector

call fsi_type_init_(fsi)

! Count variables first
if(iamroot_) write(6,*) myname_,': FSI scaling factors'
do ii=1,fsi%ncri
   read(utable(ii),*) criterium, var, instrument, lata, latb, levl, levh, sfactor
   fsi%criterium(ii) = trim(criterium)
   fsi%variable(ii) = trim(var)
   fsi%instrument(ii) = instrument
   fsi%lata(ii) = lata
   fsi%latb(ii) = latb
   fsi%levl(ii) = levl
   fsi%levh(ii) = levh
   fsi%sfactor(ii) = sfactor
   if(iamroot_) then
      write(6,'(1x,2(a10,1x),i4,1x,4f15.3,1x,f7.2)') trim(criterium), trim(var), instrument, lata, latb, levl, levh, sfactor
   endif
enddo

! release table
deallocate(utable)

fsi_weight = .true.
initialized_=.true.

end subroutine init_
!EOC

subroutine fsi_type_init_(this) 
 type(fsi_type) this
 integer n
 n=this%ncri
 allocate(this%criterium(n))
 allocate(this%variable(n))
 allocate(this%instrument(n))
 allocate(this%lata(n))
 allocate(this%latb(n))
 allocate(this%levl(n))
 allocate(this%levh(n))
 allocate(this%sfactor(n))
end subroutine fsi_type_init_

subroutine fsi_type_clean_(this) 
 type(fsi_type) this
 deallocate(this%sfactor)
 deallocate(this%levh)
 deallocate(this%levl)
 deallocate(this%latb)
 deallocate(this%lata)
 deallocate(this%instrument)
 deallocate(this%variable)
 deallocate(this%criterium)
end subroutine fsi_type_clean_

subroutine apply_lev_weight_(sfactor,var,kx,lat,lon,lev)
implicit none
real(r_kind), intent(inout) :: sfactor
character(len=*), intent(in) :: var
real(r_kind), intent(in) :: lat,lon
real(r_kind), intent(in) :: lev ! mb
integer, intent(in) :: kx
!
logical :: done
integer :: ii,iiall
character(len=4) :: criterium

! if not intialized, nothing to do ...
if(.not. initialized_ ) return

sfactor = one
iiall=-1
done=.false.
do ii=1,fsi%ncri
   if(trim(fsi%variable(ii))=='all') then
      iiall=ii
      cycle ! see special handle below
   endif
   if(trim(fsi%variable(ii)) == trim(var)) then
     if(fsi%instrument(ii) == kx) then
       if( lat>fsi%lata(ii) .and. lat<fsi%latb(ii) .and. &
           lev<fsi%levl(ii) .and. lev>fsi%levh(ii) ) then
           criterium = fsi%criterium(ii)
           sfactor = fsi%sfactor(ii)
           done = .true.
!          print *, 'DEBUG_RT: ', trim(var), kx, sfactor
           exit
       endif
     endif
   endif
enddo
! special handle case of all others
if (.not. done) then
  if(iiall>0) then
    ii=iiall
    if( .not. done .and. &
        lat>fsi%lata(ii) .and. lat<fsi%latb(ii) .and. &
        lev<fsi%levl(ii) .and. lev>fsi%levh(ii) ) then
        criterium = fsi%criterium(ii)
        sfactor = fsi%sfactor(ii)
        done = .true.
    endif
  endif
endif

end subroutine apply_lev_weight_

subroutine apply_chn_weight_(sfactor,dplat,obstype,ich,lat,lon)
implicit none
real(r_kind),     intent(inout) :: sfactor
character(len=*), intent(in) :: dplat ! metop-c - not yet implemented
character(len=*), intent(in) :: obstype ! iasi - not yet implemented
real(r_kind),     intent(in) :: lat,lon
integer(i_kind),  intent(in) :: ich
!
logical :: done
integer :: ii,iirad
character(len=4) :: criterium

! if not intialized, nothing to do ...
if(.not. initialized_ ) return

sfactor = one
iirad=-1
done=.false.
do ii=1,fsi%ncri
   if(trim(fsi%variable(ii))=='rad') then
      iirad=ii
      cycle ! see special handle below
   endif
!   if(trim(fsi%variable(ii)) == trim(var)) then
!     if(fsi%instrument(ii) == kx) then
!       if( lat>fsi%lata(ii) .and. lat<fsi%latb(ii) .and. &
!           lev<fsi%levl(ii) .and. lev>fsi%levh(ii) ) then
!           criterium = fsi%criterium(ii)
!           sfactor = fsi%sfactor(ii)
!           done = .true.
!           exit
!       endif
!     endif
!   endif
enddo
! special handle case of all others
if (.not. done) then
  if(iirad>0) then
    ii=iirad
    if( lat>fsi%lata(ii) .and. lat<fsi%latb(ii) ) then
        criterium = fsi%criterium(ii)
        sfactor = fsi%sfactor(ii)
        done = .true.
    endif
  endif
endif

end subroutine apply_chn_weight_

end module m_fsi_weight
