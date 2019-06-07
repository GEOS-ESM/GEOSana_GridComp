
!NOTE: To do original code defined _MAX_ORIGINAL_ and _TIME_LOOP_
!=============================================================================
!BOP

! !MODULE: 

! !INTERFACE:

module GEOS_AgcmPertADTester
#include "MAPL_Generic.h"

! !USES:

!  use ESMF
!  use MAPL_Mod
  use GEOS_PertSharedMod, only: T_3Dvar

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public:: ADTester_test

  interface ADTester_test; module procedure test_; end interface

! !DESCRIPTION:
!
!EOP

  real   *4 , parameter:: RS=0.
  real   *8 , parameter:: RD=0.
  real   *16, parameter:: RQ=0.
  integer*4 , parameter:: I4=0
  integer   , parameter:: SP=kind(RS)
  integer   , parameter:: DP=kind(RD)
  integer   , parameter:: QP=kind(RQ)
  integer   , parameter:: IK=kind(I4)
  integer   , parameter:: MAXNAMELEN = 64

  character(len=*), parameter:: myname="GEOS_AgcmPertADTester"
  
  character(len=MAXNAMELEN), save:: name_adtester="TLM"
  real(SP), target, allocatable, dimension(:,:,:,:,:,:),save :: save_adtester
  real(QP), save:: dpp_adtester
  real(QP), save:: dqq_adtester
  real(QP), save:: dpq_adtester
  integer , parameter:: lbnd_adtester = 1
  integer , parameter:: ubnd_adtester = 4

#include "mpif.h"
contains

  subroutine test_(pertWrk, icall,interval,phase,name)
    use MAPL_Mod, only: MAPL_ASRT, MAPL_AM_I_ROOT
    use ESMF,     only: ESMF_SUCCESS
    use ESMF,     only: ESMF_VM
    use ESMF,     only: ESMF_VMgetcurrent
    use ESMF,     only: ESMF_VMget
    use GEOS_PertSharedMod, only: T_3DVar,TLMphase,ADJphase
    implicit none
    type(T_3DVar), target, dimension(:), intent(in) :: pertWrk
    integer,intent(in):: icall
    integer,intent(in):: interval
    integer,intent(in):: phase
    character(len=*),optional,intent(in):: name

    integer:: RC
    character(len=*),parameter:: myname_=myname//'::test_'
    character(len=*),parameter:: IAm=myname_
    type(ESMF_VM), save :: VM


    integer(IK):: mcomm
    integer(IK):: mx,my,mz,mv,iv,ier
    real(kind=QP):: dpp,dqq,dpq,d(3)
    real(kind=QP):: cpq,rpq
    integer(IK):: ipq,mpq

    real(kind(save_adtester)),pointer,dimension(:,:,:):: p,q
    
    if( interval<lbnd_adtester .or. interval>ubnd_adtester) then
      if(MAPL_AM_I_ROOT()) then
	write(*,'(1x,a,3(a,i2.2))') trim(Iam), &
      	  '() --                   p.i.i = ',phase,'.',interval,'.',icall
      endif
      !! tlm
      return
    endif

    mv=size(pertWrk)
    if(.not.allocated(save_adtester)) then
      mx=0; my=0; mz=0
      if(mv>0) then
        mx=size(pertwrk(1)%x,1)
        my=size(pertwrk(1)%x,2)
        mz=size(pertwrk(1)%x,3)
      endif

      allocate(save_adtester(mx,my,mz,mv, 0:1, lbnd_adtester:ubnd_adtester))
      if(MAPL_AM_I_ROOT()) then
	write(*,'(1x,2a,7(i4))') myname_, &
      	  '() -- shape(save_adtester) = ' ,shape(save_adtester)
	write(*,'(1x,2a,3a   )') myname_, &
      	  '() --       name_adtester  = ',' "',trim(name_adtester),'"'
	write(*,'(1x,2a,i4   )') myname_, &
      	  '() --       lbnd_adtester  = ',lbnd_adtester
	write(*,'(1x,2a,i4   )') myname_, &
      	  '() --       ubnd_adtester  = ',ubnd_adtester
      endif
    endif

    select case(phase)
    case(TLMphase)
!      if(MAPL_AM_I_ROOT()) then
!	write(*,'(1x,a,3(a,i2.2))') trim(Iam), &
!      	  '() --          saved at p.i.i = ',phase,'.',interval,'.',icall
!      endif

      ASSERT_(allocated(save_adtester))
      do iv=1,mv
        save_adtester(:,:,:,iv,icall,interval) = pertwrk(iv)%x(:,:,:)
      enddo

    case(ADJphase)
      dpp=0.
      dqq=0.
      dpq=0.

      ASSERT_(allocated(save_adtester))
      do iv=1,mv
        p => save_adtester(:,:,:,iv,icall,interval)	! this is the original TL pert-vector
        q => pertwrk(iv)%x(:,:,:)			! this is current AD pert-vector

	dpp = dpp+sum(p(:,:,:)*p(:,:,:))
	dqq = dqq+sum(q(:,:,:)*q(:,:,:))
	dpq = dpq+sum(p(:,:,:)*q(:,:,:))
      enddo

      nullify(p,q)

      call ESMF_VMgetcurrent(VM,RC=RC)
      			ASSERT_(RC==ESMF_SUCCESS)
      call ESMF_VMget(VM,mpicommunicator=mcomm,RC=RC)
      			ASSERT_(RC==ESMF_SUCCESS)

      d(1)=dpp; d(2)=dqq; d(3)=dpq
      call mpi_allreduce((d),d,size(d),MPI_REAL16,MPI_SUM,mcomm,RC)
      			ASSERT_(RC==MPI_SUCCESS)
      dpp=d(1); dqq=d(2); dpq=d(3)

      select case(icall)
      case(1)
	if(MAPL_AM_I_ROOT()) then
	  write(*,'(1x,4a,i1,a,i2.2,a,i1,2x,4(1x,e16.8),2i6)') trim(Iam), &
		'M=',trim(name_adtester),			    &
		"() -- (u    ,v=Mq), uu, vv = ", phase,'.',interval,'.',icall, dpq, dqq, dpp
	endif

	dpp_adtester=dpp
	dqq_adtester=dqq
	dpq_adtester=dpq

      case(0)
	cpq=1._QP
	if(abs(dpq_adtester)>0._QP) cpq=dpq/dpq_adtester
	rpq=cpq-1._QP
	! ipq=int(-log(abs(rpq)+tiny(rpq))/log(10._QP))
	ipq=int(-log(max(abs(rpq),tiny(1.)))/log(10.))

	rpq=tiny(1.)
	rpq=max(dpp_adtester*dqq_adtester,dpp*dqq,rpq*rpq)
	rpq=rpq**0.5_QP
	if (rpq.ne.0.) rpq=(dpq-dpq_adtester)/rpq

	mpq=int(-log(max(abs(rpq),tiny(1.)))/log(10.))

	if(MAPL_AM_I_ROOT()) then
	  write(*,'(1x,4a,i1,a,i2.2,a,i1,2x,4(1x,e16.8),2i6)') trim(Iam), &
		'M=',trim(name_adtester),			    &
		"() -- (p=M'u,q   ), qq, pp = ", phase,'.',interval,'.',icall, dpq, dqq, dpp, rpq, ipq, mpq
	endif

       if(interval==lbnd_adtester) then
         ASSERT_(allocated(save_adtester))
         deallocate(save_adtester)
       endif

      case default
	ASSERT_(.false.)
      endselect

    case default
      ASSERT_(.false.)
    endselect
  end subroutine test_

  function dot_pert_(pertWrk) result(dot_)
    use GEOS_PertSharedMod, only: T_3DVar
    use MAPL_Mod, only: MAPL_ASRT
    implicit none
    type(T_3DVar), target, dimension(:), intent(in) :: pertWrk
    real(QP):: dot_
    integer:: RC,iv
    character(len=*),parameter:: myname_=myname//'::dot_pert_'
    character(len=*),parameter:: IAm=myname_

    dot_=0._QP
    do iv=1,size(pertWrk)
      dot_ = dot_+sum(pertWrk(iv)%X(:,:,:)*pertWrk(iv)%X(:,:,:))
    enddo

    call mpi_allreduce((dot_),dot_,1,MPI_REAL16,MPI_SUM,MPI_comm_World,RC)
							ASSERT_(RC==0)
  end function dot_pert_
end module GEOS_AgcmPertADTester
