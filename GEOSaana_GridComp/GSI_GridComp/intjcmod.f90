module intjcmod
!$$$ module documentation block
!           .      .    .                                       .
! module:   intjcmod    module for weak constraint int routines
!  pgrmmr:  kleist
!
! abstract: module for Jc int routines
!
! program history log:
!   2012-01-21  kleist - consolidation of Jc int routines into single module
!   2013-10-25  todling - nullify work pointers
!
! subroutines included:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
use kinds, only: r_kind,i_kind,r_quad
use constants, only: zero,two,one,half,zero_quad,one_quad,two_quad
use gsi_bundlemod, only: gsi_bundle,gsi_bundlegetpointer

implicit none

PRIVATE
PUBLIC intlimq,intlimg,intlimp,intlimv,intjcdfi,intjcpdry

contains

subroutine intlimq(rval,sval,itbin)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intlimq
!   prgmmr: derber           org: np23                date: 1996-11-19
!
! abstract: limit negative q as a weak constraint
!
! program history log:
!   1996-11-19  derber
!   1998-07-10  weiyu yang
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   2004-03-15  kleist, d., derber, j., treadon, r., use negative q only
!   2004-06-02  kleist, add penalty for excess moisture
!   2004-08-02  treadon - add only to module use, add intent in/out
!   2007-02-13  derber - modify to use rh rather than q
!   2008-06-02  safford - rm unused vars
!   2010-05-13  todling - update to use gsi_bundle
!   2011-12-27  kleist - add multiple time level capability (for 4densvar option)
!
!   input argument list:
!     sq       - increment in grid space
!     itbin    - observation bin (time level)
!
!   output argument list:
!     rq       - results from limiting operator                 
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use gridmod, only: lat2,lon2,nsig,lat1,lon1
  use jfunc, only: factqmin,factqmax
  use derivsmod, only: qgues,qsatg
  use guess_grids, only: ges_qsat
  use gsi_metguess_mod, only: gsi_metguess_bundle
  implicit none

! Declare passed variables
  type(gsi_bundle),intent(in   ) :: sval
  type(gsi_bundle),intent(inout) :: rval
  integer, intent(in)            :: itbin

! Declare local variables
  integer(i_kind) i,j,k,ier,istatus
  real(r_kind) q
  real(r_kind),pointer,dimension(:,:,:) :: sq=>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: rq=>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: ges_q_it=>NULL()

  if (factqmin==zero .and. factqmax==zero) return

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'q',sq,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'q',rq,istatus);ier=istatus+ier
  if(ier/=0)return

  call gsi_bundlegetpointer(gsi_metguess_bundle(itbin),'q',ges_q_it,ier)
  if(ier/=0)return
 
!$omp parallel do  schedule(dynamic,1) private(k,j,i,q)
  do k = 1,nsig
     do j = 2,lon1+1
        do i = 2,lat1+1
           q = ges_q_it(i,j,k) + sq(i,j,k)
           
!          Lower constraint limit
           if (q < zero) then
              rq(i,j,k) = rq(i,j,k) + factqmin*q/(ges_qsat(i,j,k,itbin)*ges_qsat(i,j,k,itbin))

!          Upper constraint limit
           else if (q > ges_qsat(i,j,k,itbin)) then
              rq(i,j,k) = rq(i,j,k) + factqmax*(q-ges_qsat(i,j,k,itbin))/ &
                          (ges_qsat(i,j,k,itbin)*ges_qsat(i,j,k,itbin))
           
           end if
        end do
     end do
  end do
  
  return
end subroutine intlimq

subroutine intlimg(rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intlimg
!   prgmmr: zhu           org: np23                date: 2011-02-20
!
! abstract: limit negative gust as a weak constraint
!
! program history log:
!   2011-02-20  zhu
!
!   input argument list:
!     sg       - increment in grid space
!
!   output argument list:
!     rg       - results from limiting operator                 
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: zero
  use gridmod, only: lat2,lon2,nsig,lat1,lon1
  use jfunc, only: factg
  use derivsmod, only: ggues
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  implicit none

! Declare passed variables
  type(gsi_bundle),intent(in   ) :: sval
  type(gsi_bundle),intent(inout) :: rval

! Declare local variables
  integer(i_kind) i,j,ier,istatus
  real(r_kind) gust
  real(r_kind),pointer,dimension(:,:) :: sg=>NULL()
  real(r_kind),pointer,dimension(:,:) :: rg=>NULL()

  if (factg==zero) return

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'gust',sg,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'gust',rg,istatus);ier=istatus+ier
  if(ier/=0)return
 
  do j = 2,lon1+1
     do i = 2,lat1+1
        gust = ggues(i,j) + sg(i,j)
           
!       Lower constraint limit
        if (gust < zero) then
           rg(i,j) = rg(i,j) + factg*gust/(ggues(i,j)*ggues(i,j))
        end if
     end do
  end do
  
  return
end subroutine intlimg

subroutine intlimp(rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intlimp
!   prgmmr: zhu           org: np23                date: 2011-02-20
!
! abstract: limit negative pblh as a weak constraint
!
! program history log:
!   2011-02-20  zhu
!
!   input argument list:
!     sp       - increment in grid space
!
!   output argument list:
!     rp       - results from limiting operator                 
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use gridmod, only: lat2,lon2,nsig,lat1,lon1
  use jfunc, only: factp
  use derivsmod, only: pgues
  implicit none

! Declare passed variables
  type(gsi_bundle),intent(in   ) :: sval
  type(gsi_bundle),intent(inout) :: rval

! Declare local variables
  integer(i_kind) i,j,ier,istatus
  real(r_kind) pblh
  real(r_kind),pointer,dimension(:,:) :: sp=>NULL()
  real(r_kind),pointer,dimension(:,:) :: rp=>NULL()

  if (factp==zero) return

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblh',sp,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblh',rp,istatus);ier=istatus+ier
  if(ier/=0)return
 
  do j = 2,lon1+1
     do i = 2,lat1+1
        pblh = pgues(i,j) + sp(i,j)
           
!       Lower constraint limit
        if (pblh < zero) then
           rp(i,j) = rp(i,j) + factp*pblh/(pgues(i,j)*pgues(i,j))
        end if
     end do
  end do
  
  return
end subroutine intlimp

subroutine intlimv(rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intlimv
!   prgmmr: zhu           org: np23                date: 2011-02-20
!
! abstract: limit negative vis as a weak constraint
!
! program history log:
!   2011-02-20  zhu
!
!   input argument list:
!     sv       - increment in grid space
!
!   output argument list:
!     rv       - results from limiting operator                 
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use gridmod, only: lat2,lon2,nsig,lat1,lon1
  use jfunc, only: factv
  use derivsmod, only: vgues
  implicit none

! Declare passed variables
  type(gsi_bundle),intent(in   ) :: sval
  type(gsi_bundle),intent(inout) :: rval

! Declare local variables
  integer(i_kind) i,j,ier,istatus
  real(r_kind) vis
  real(r_kind),pointer,dimension(:,:) :: sv=>NULL()
  real(r_kind),pointer,dimension(:,:) :: rv=>NULL()

  if (factv==zero) return

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'vis',sv,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'vis',rv,istatus);ier=istatus+ier
  if(ier/=0)return
 
  do j = 2,lon1+1
     do i = 2,lat1+1
        vis = vgues(i,j) + sv(i,j)
           
!       Lower constraint limit
        if (vis < zero) then
           rv(i,j) = rv(i,j) + factv*vis/(vgues(i,j)*vgues(i,j))
        end if
     end do
  end do
  
  return
end subroutine intlimv

subroutine intjcpdry(rval,sval,pjc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intjcpdry   adjoint for mean dry ps conservation
!   prgmmr: kleist           org: np23                date: 2009-07-07
!
! abstract: calculate stepsize contribution and penalty for limiting changes to dry mass
!
! program history log:
!   2009-07-07  kleist
!   2010-05-13  todling - update to use gsi_bundle
!   2010-05-25  derber  - modify to minimize number of communications
!   2010-08-18  hu      - added qpvals= to mpl_allreduce call
!   2010-11-03  treadon - correct i,j loop limits for rq,rc update
!   2011-11-01  eliu    - add handling for ql & qi increments and search directions
!   2013-05-05  todling - separate dry mass from the rest (zero-diff change)
!                         collapse two verions of this routine into one (add opt arg)
!
!   input argument list:
!     rq       - q search direction
!     rc       - cloud water search direction
!     rp       - surface pressure search direction
!     sq       - q increment
!     sc       - cloud water increment
!     sp       - increment in grid space
!     mype     - integer PE
!
!   output argument list:
!     rq       - q search direction
!     rc       - cloud water search direction
!     rp       - surface pressure search direction
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use mpimod, only: mype
  use gridmod, only: lat2,lon2,nsig,wgtlats,nlon,istart
  use guess_grids, only: ges_prsi,ntguessig
  use mpl_allreducemod, only: mpl_allreduce
  use jcmod, only: bamp_jcpdry
  use gsi_metguess_mod,  only: gsi_metguess_get
  implicit none

! Declare passed variables
  type(gsi_bundle),intent(in   ) :: sval
  type(gsi_bundle),intent(inout) :: rval
  real(r_quad)    ,intent(  out),optional :: pjc

! Declare local variables
  real(r_quad),dimension(2) :: mass ! 1=dry;2=wv
  real(r_quad),dimension(nsig) :: mass2
  real(r_quad) rcon,con,dmass
  integer(i_kind) i,j,k,it,ii,mm1,ier,icw,iql,iqi,istatus
  real(r_kind),pointer,dimension(:,:,:) :: sq =>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: sc =>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: sql=>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: sqi=>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: rq =>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: rc =>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: rql=>NULL()
  real(r_kind),pointer,dimension(:,:,:) :: rqi=>NULL()
  real(r_kind),pointer,dimension(:,:)   :: sp =>NULL()
  real(r_kind),pointer,dimension(:,:)   :: rp =>NULL()
  
  it=ntguessig

! Retrieve pointers
! Simply return if any pointer not found
  ier=0; icw=0; iql=0; iqi=0
  call gsi_bundlegetpointer(sval,'q' ,sq, istatus);ier=istatus+ier
  call gsi_bundlegetpointer(sval,'cw',sc, istatus);icw=istatus+icw
  call gsi_bundlegetpointer(sval,'ql',sql,istatus);iql=istatus+iql
  call gsi_bundlegetpointer(sval,'qi',sqi,istatus);iqi=istatus+iqi
  call gsi_bundlegetpointer(sval,'ps',sp, istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'q' ,rq, istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'cw',rc, istatus);icw=istatus+icw
  call gsi_bundlegetpointer(rval,'ql',rql,istatus);iql=istatus+iql
  call gsi_bundlegetpointer(rval,'qi',rqi,istatus);iqi=istatus+iqi
  call gsi_bundlegetpointer(rval,'ps',rp, istatus);ier=istatus+ier
  if(ier+icw*(iql+iqi)/=0)then
    if (mype==0) write(6,*)'intjcpdry: checking ier+icw*(iql+iqi)=', ier+icw*(iql+iqi)
    return
  end if

  mass(:)=zero_quad
  rcon=one_quad/(two_quad*float(nlon))
  mm1=mype+1

! Calculate mean surface pressure contribution in subdomain
  do j=2,lon2-1
    do i=2,lat2-1
      ii=istart(mm1)+i-2
      mass(1)=mass(1)+sp(i,j)*wgtlats(ii)
    end do
  end do

  mass2(:)=zero_quad
! Calculate water-vapor contribution to total mass
!$omp parallel do  schedule(dynamic,1) private(k,j,i,ii,con)
  do k=1,nsig
     do j=2,lon2-1
        do i=2,lat2-1
           ii=istart(mm1)+i-2
           con = (ges_prsi(i,j,k,it)-ges_prsi(i,j,k+1,it))*wgtlats(ii)
           mass2(k)=mass2(k)+sq(i,j,k)*con
           if (icw==0) then
              mass2(k)=mass2(k)+sc(i,j,k)*con
           else
              mass2(k)=mass2(k)+(sql(i,j,k)+sqi(i,j,k))*con
           endif
        end do
     end do
  end do
  do k=1,nsig
     mass(2)=mass(2)+mass2(k)
  end do

! First, use MPI to get global mean increment
  call mpl_allreduce(2,qpvals=mass)

! Remove water-vapor contribution to get incremental dry ps
! if (mype==0) write(6,*)'intjcpdry: total mass =', mass(1)
! if (mype==0) write(6,*)'intjcpdry: wv    mass =', mass(2)
  dmass=mass(1)-mass(2)
  dmass=bamp_jcpdry*dmass*rcon*rcon
  if(present(pjc)) then
     pjc = dmass*dmass
  endif

! Calculate mean surface pressure contribution in subdomain
  do j=2,lon2-1
    do i=2,lat2-1
      ii=istart(mm1)+i-2
      rp(i,j)=rp(i,j)+dmass*wgtlats(ii)
    end do
  end do
! Remove water to get incremental dry ps
!$omp parallel do  schedule(dynamic,1) private(k,j,i,ii,con)
  do k=1,nsig
     do j=2,lon2-1
        do i=2,lat2-1
           ii=istart(mm1)+i-2
           con = dmass*(ges_prsi(i,j,k,it)-ges_prsi(i,j,k+1,it))*wgtlats(ii)
           rq(i,j,k)=rq(i,j,k)-con
           if (icw==0)then
              rc(i,j,k)=rc(i,j,k)-con
           else
              rql(i,j,k)=rql(i,j,k)-con
              rqi(i,j,k)=rqi(i,j,k)-con
           endif
        end do
     end do
  end do

  return
end subroutine intjcpdry

subroutine intjcdfi(rval,sval,pjc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intjcdfi    calculate Jc DFI terms and contribution to gradient
!   prgmmr: tremolet
!
! program history log:
!   2007-10-18  tremolet - initial version
!   2009-01-18  todling  - carry summation in quad precision
!   2009-08-14  lueken   - update documentation
!   2010-05-14  todling  - update to use gsi_bundle
!   2011-08-01  lueken   - replace F90 with f90 (no machine logic) and replaced izero/ione with 0/1
!   2012-01-19  kleist   - adaptation of evaljcdfi
!   2013-05-18  todling  - code consolidation; remove evaljcdfi and renamed to intjcdfi
!
!   input argument list:
!    sval
!    rval
!
!   output argument list:
!    rval
!    pjc
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use jcmod, only: wgtdfi,alphajc
use gsi_4dvar, only: nobs_bins
use mpimod, only: mype
use state_vectors, only : allocate_state,deallocate_state
use gsi_bundlemod, only : self_add,self_mul,assignment(=)
implicit none

! Declare passed variables
type(gsi_bundle),dimension(nobs_bins),intent(in   ) :: sval
type(gsi_bundle),dimension(nobs_bins),intent(inout) :: rval
real(r_quad),               optional, intent(  out) :: pjc

! Declare local variables
integer(i_kind) :: jj,idfi
real(r_quad),parameter :: half_quad=0.5_r_quad
type(gsi_bundle) :: sfilter,afilter
real(r_quad) :: cost

!************************************************************************************  

  idfi = (nobs_bins-1)/2+1
  call allocate_state(sfilter)
  call allocate_state(afilter)

! Compute filtered state
  sfilter=zero
  do jj=1,nobs_bins
     call self_add(sfilter,wgtdfi(jj),sval(jj))
  enddo

! Compute difference from filtered state
  call self_add(sfilter,-one,sval(idfi))

! Apply Jc multiplicative factor
  call self_mul(sfilter,alphajc)

! Compute Jc (norm of difference)
! Jc = 1/2 * wgt * sfilter *sfilter
! afilter = wgt * sfilter
  call enorm_state(sfilter,cost,afilter)
  if(present(pjc))then
     pjc=half_quad*cost
     if (mype==0) write(6,*)'Jc DFI=',pjc
  endif

! Adjoint Jc multiplicative factor
  call self_mul(afilter,alphajc)

! Adjoint of difference from filtered state
  call self_add(rval(idfi),-one,afilter)

! Compute filtered state
  do jj=1,nobs_bins
     call self_add(rval(jj),wgtdfi(jj),afilter)
  enddo

  call deallocate_state(sfilter)
  call deallocate_state(afilter)

  return
end subroutine intjcdfi

end module intjcmod
