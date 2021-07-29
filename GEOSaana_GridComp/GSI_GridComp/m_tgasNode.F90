module m_tgasNode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_tgasNode
!   prgmmr:	 
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type tgasNode (trace gas)
!
! program history log:
!   2021-04-25  j guo   - added this document block for the initial polymorphic
!                         implementation.
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
  use m_obsdiagNode, only: obs_diag,aofp_obs_diag => fptr_obsdiagNode
  use m_obsdiagNode, only: obs_diags
  use kinds , only: i_kind,r_kind
  use mpeu_util, only: assert_,die,perr,warn,tell
  use m_obsNode, only: obsNode
  implicit none
  private

  public:: tgasNode

  type,extends(obsNode):: tgasNode
     type(aofp_obs_diag), dimension(:),   pointer :: diags   => NULL()

     real(r_kind),        dimension(:),   pointer :: res     => NULL()          ! residual
     real(r_kind),        dimension(:),   pointer :: err2    => NULL()          ! error squared
     real(r_kind),        dimension(:),   pointer :: raterr2 => NULL()          ! square of ratio of final obs error to original obs error
     real(r_kind),        dimension(:,:), pointer :: avgker  => NULL()          ! averaging kernel
     real(r_kind),        dimension(:,:), pointer :: avgwgt  => NULL()          ! averaging weights
     integer(i_kind),     dimension(:),   pointer :: ipos    => NULL()

     integer(i_kind)    :: nchanl         ! number of channels for this profile
     integer(i_kind)    :: npro           ! number of elements for this profile

     integer(i_kind)    :: ij(4)          ! horizontal locations
     real(r_kind)       :: wij(4)         ! horizontal interpolation weights
     character(len=256) :: obstype        ! observation type of int/stp procedures

  contains
    procedure,nopass::  mytype
    procedure::  setHop => obsNode_setHop_
    procedure::   xread => obsNode_xread_
    procedure::  xwrite => obsNode_xwrite_
    procedure:: isvalid => obsNode_isvalid_
    procedure::  getTLDdp => getTLDdp_

    ! procedure, nopass:: headerRead  => obsHeader_read_
    ! procedure, nopass:: headerWrite => obsHeader_write_
    ! procedure:: init  => obsNode_init_
    procedure:: clean => obsNode_clean_

  end type tgasNode

  public:: tgasNode_typecast
  public:: tgasNode_nextcast
        interface tgasNode_typecast; module procedure typecast_ ; end interface
        interface tgasNode_nextcast; module procedure nextcast_ ; end interface

  public:: tgasNode_appendto
        interface tgasNode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: MYNAME="m_tgasNode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(aNode) result(ptr_)
!-- cast a class(obsNode) to a type(tgasNode)
  use m_obsNode, only: obsNode
  implicit none
  type(tgasNode),pointer:: ptr_
  class(obsNode ),pointer,intent(in):: aNode
  ptr_ => null()
  if(.not.associated(aNode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(aNode)
  type is(tgasNode)
    ptr_ => aNode
  end select
return
end function typecast_

function nextcast_(aNode) result(ptr_)
!-- cast an obsNode_next(obsNode) to a type(tgasNode)
  use m_obsNode, only: obsNode,obsNode_next
  implicit none
  type(tgasNode),pointer:: ptr_
  class(obsNode ),target ,intent(in):: aNode

  class(obsNode),pointer:: inode_
  inode_ => obsNode_next(aNode)
  ptr_ => typecast_(inode_)
return
end function nextcast_

subroutine appendto_(aNode,oll)
  use m_obsNode , only: obsNode
  use m_obsLList, only: obsLList,obsLList_appendNode
  implicit none
  type(tgasNode),pointer,intent(in):: aNode
  type(obsLList),intent(inout):: oLL

  class(obsNode),pointer:: inode_
  inode_ => aNode
  call obsLList_appendNode(oLL,inode_)
  inode_ => null()
end subroutine appendto_

! obsNode implementations

function mytype()
  implicit none
  character(len=:),allocatable:: mytype
  mytype="[tgasNode]"
end function mytype

subroutine obsNode_clean_(aNode)
  implicit none
  class(tgasNode),intent(inout):: aNode

  character(len=*),parameter:: myname_=MYNAME//'::obsNode_clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',aNode%mytype())
    if(associated(aNode%diags  )) deallocate(aNode%diags  )
    if(associated(aNode%res    )) deallocate(aNode%res    )
    if(associated(aNode%err2   )) deallocate(aNode%err2   )
    if(associated(aNode%raterr2)) deallocate(aNode%raterr2)
    if(associated(aNode%avgker )) deallocate(aNode%avgker )
    if(associated(aNode%avgwgt )) deallocate(aNode%avgwgt )
    if(associated(aNode%ipos   )) deallocate(aNode%ipos   )
_EXIT_(myname_)
return
end subroutine obsNode_clean_

subroutine obsNode_xread_(aNode,iunit,istat,diagLookup,skip)
  use gridmod, only: nsig
  use m_obsdiagNode, only: obsdiagLookup_locate
  implicit none
  class(tgasNode) , intent(inout):: aNode
  integer(i_kind) , intent(in   ):: iunit
  integer(i_kind) , intent(  out):: istat
  type(obs_diags) , intent(in   ):: diagLookup
  logical,optional, intent(in   ):: skip

  character(len=*),parameter:: myname_=MYNAME//'::obsNode_xread_'
  integer(i_kind),allocatable,dimension(:):: ich
  integer(i_kind):: k,nchanl,npro,nsig_read
  logical:: skip_
_ENTRY_(myname_)
  skip_=.false.
  if(present(skip)) skip_=skip

  istat=0
  if(skip_) then
    read(iunit,iostat=istat)
                if (istat/=0) then
                  call perr(myname_,'skipping read(%nchanl, ...), istat =',istat)
                  _EXIT_(myname_)
                  return
                end if
    read(iunit,iostat=istat)
                if(istat/=0) then
                  call perr(myname_,'skipping read(ich, %(res, ...)), istat =',istat)
                  _EXIT_(myname_)
                  return
                endif

  else
    read(iunit,iostat=istat) aNode%nchanl,aNode%npro,nsig_read,aNode%obstype	! nchanl,npro,nsig,obstype
                if (istat/=0) then
                  call perr(myname_,'read(%nchanl, ...), istat =',istat)
                  _EXIT_(myname_)
                  return
                end if

	if(nsig/=nsig_read) then
          call perr(myname_,'read(%nchanl, ,,,), expecting nsig =',nsig)
          call perr(myname_,'                     but nsig_read =',nsig_read)
          _EXIT_(myname_)
          return
	endif

                ! re-allocation is needed, because a node may have been read in
                ! but excluded later.
        if(associated(aNode%diags  )) deallocate(aNode%diags  )
        if(associated(aNode%res    )) deallocate(aNode%res    )
        if(associated(aNode%err2   )) deallocate(aNode%err2   )
        if(associated(aNode%raterr2)) deallocate(aNode%raterr2)
        if(associated(aNode%ipos   )) deallocate(aNode%ipos   )
        if(associated(aNode%avgker )) deallocate(aNode%avgker )
        if(associated(aNode%avgwgt )) deallocate(aNode%avgwgt )

    nchanl = aNode%nchanl
    npro   = aNode%npro

        allocate( aNode%res    (nchanl) &
                 ,aNode%diags  (nchanl) &
                 ,aNode%err2   (nchanl) &
                 ,aNode%raterr2(nchanl) &
                 ,aNode%ipos   (nchanl) & 
                 ,aNode%avgker (nchanl,npro) & 
                 ,aNode%avgwgt (npro  ,nsig) &
		)

        allocate(       ich    (nchanl) )

    read(iunit,iostat=istat)          ich    , & !(nchanl)
                                aNode%res    , & !(nchanl)
                                aNode%err2   , & !(nchanl)
                                aNode%raterr2, & !(nchanl)
                                aNode%ipos   , & !(nchanl)
                                aNode%avgker , & !(nchanl,npro)
                                aNode%avgwgt     !(npro  ,nsig)
                if (istat/=0) then
                  call perr(myname_,'read(ich,%(...)), istat =',istat)
                  _EXIT_(myname_)
                  return
                end if

    do k=1,nchanl
      ! Not sure if ich=k or ich=ich(k) for now.  More verification is needed.
      aNode%diags(k)%ptr => obsdiagLookup_locate(diagLookup,aNode%idv,aNode%iob,ich(k))
                if(.not.associated(aNode%diags(k)%ptr)) then
                  call perr(myname_,'obsdiagLookup_locate(k), k =',k)
                  call perr(myname_,'                      %idv =',aNode%idv)
                  call perr(myname_,'                      %iob =',aNode%iob)
                  call perr(myname_,'                       ich =',ich(k))
                  call  die(myname_)
                  istat=-1
                  _EXIT_(myname_)
                  return
                endif
    enddo
    deallocate(ich)
  endif
_EXIT_(myname_)
return
end subroutine obsNode_xread_

subroutine obsNode_xwrite_(aNode,junit,jstat)
  implicit none
  class(tgasNode),intent(in):: aNode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=MYNAME//'::obsNode_xwrite_'
  integer(i_kind):: k
_ENTRY_(myname_)

  jstat=0
  write(junit,iostat=jstat) aNode%nchanl,aNode%npro,size(aNode%avgwgt,2),aNode%obstype	! nchanl,npro,nsig
                if (jstat/=0) then
                  call perr(myname_,'write(%nchanl, ...), jstat =',jstat)
                  _EXIT_(myname_)
                  return
                end if

  write(junit,iostat=jstat) (/ (aNode%diags(k)%ptr%ich, k=1,aNode%nchanl) /), &
                                aNode%res    , & !(%nchanl)
                                aNode%err2   , & !(%nchanl)
                                aNode%raterr2, & !(%nchanl)
                                aNode%ipos   , & !(%nchanl)
				aNode%avgker , & !(%nchanl,%npro)
				aNode%avgwgt     !(%npro,nsig)
                if (jstat/=0) then
                  call perr(myname_,'write(ich,%(res,...)), istat =',jstat)
                  _EXIT_(myname_)
                  return
                end if
_EXIT_(myname_)
return
end subroutine obsNode_xwrite_

subroutine obsNode_setHop_(aNode)
  !use m_obsNode, only: obstype_getHop => obsNode_getHop
  use m_cvgridLookup, only: cvgridLookup_getiw
  implicit none
  class(tgasNode),intent(inout):: aNode

  character(len=*),parameter:: myname_=MYNAME//'::obsNode_setHop_'
_ENTRY_(myname_)
  call cvgridLookup_getiw(aNode%elat,aNode%elon,aNode%ij,aNode%wij)
_EXIT_(myname_)
return
end subroutine obsNode_setHop_

function obsNode_isvalid_(aNode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(tgasNode),intent(in):: aNode

  character(len=*),parameter:: myname_=MYNAME//'::obsNode_isvalid_'
  integer(i_kind):: k
_ENTRY_(myname_)
  isvalid_ = all( (/ (associated(aNode%diags(k)%ptr), k=1,aNode%nchanl) /) )
_EXIT_(myname_)
return
end function obsNode_isvalid_

pure subroutine getTLDdp_(aNode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(tgasNode), intent(in):: aNode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  integer(kind=i_kind):: k
  do k=1,aNode%nchanl
    tlddp = tlddp + aNode%diags(k)%ptr%tldepart(jiter)*aNode%diags(k)%ptr%tldepart(jiter)
  enddo
  if(present(nob)) nob=nob+aNode%nchanl
return
end subroutine getTLDdp_

end module m_tgasNode
