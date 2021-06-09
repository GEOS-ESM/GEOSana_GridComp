!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  geos_StateIO --- IO of GSI states
!
! !INTERFACE:
!
  module geos_StateIO


  use ESMF, only: ESMF_MAXGRIDDIM
  use MAPL

  use GSI_GridCompMod, only: GSI_bkg_fname_tmpl
  use GSI_GridCompMod, only: GSI_ensbkg_fname_tmpl
  use GSI_GridCompMod, only: GSI_ensana_fname_tmpl
  use GSI_GridCompMod, only: GSI_ensprga_fname_tmpl
  use GSI_GridCompMod, only: GSI_ensprgb_fname_tmpl
  use GSI_GridCompMod, only: GSI_ensread_blocksize
  use GSI_GridCompMod, only: PPMV2GpG

  use kinds,       only : r_kind,r_single,i_kind
  use mpimod,      only : mype,mpi_rtype,mpi_comm_world
  use gridmod,     only : nsig           ! no. levels
  use gridmod,     only : bk5         
  use gsi_4dvar,   only : evfsoi_afcst    ! might want to have this as opt arg of interface        
  use constants,   only : zero,one,tiny_r_kind,grav
  use state_vectors,only: dot_product

  use m_tick, only: tick
  use mpeu_util, only: tell,warn,perr,die
  use timermod, only: timer_ini,timer_fnl

  implicit none
  private
  public:: State_get
  public:: State_put

  interface State_get
    module procedure get_1State_
    module procedure get_NState_
  end interface 
  interface State_put; module procedure&
    put_1State_; end interface 

  character(len=*),parameter :: myname="geos_StateIO"
#include "MAPL_ErrLog.h"


   character(len=*),parameter:: fnpert = 'fsens.eta.nc4'
   character(len=*),parameter:: fnxgsi = 'xxgsi.eta'

   integer,save :: mycount = 0

   integer(i_kind),parameter:: ROOT =0
   real(r_kind), parameter :: kPa_per_Pa = 0.001_r_kind
   real(r_kind), parameter :: Pa_per_kPa = 1000._r_kind
   real(r_kind), parameter :: R3600      = 3600.0_r_kind

contains
!------------------------------------------------------------------------------------
subroutine get_1State_(xx,sgrid,nymd,nhms,iwhat,tau)
!
use constants, only:zero
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: assignment(=)
use general_sub2grid_mod, only: sub2grid_info
use GSI_GridCompMod, only: GSI_ExpId
implicit none
type(gsi_bundle),   intent(inout) :: xx
type(sub2grid_info),intent(in   ) :: sgrid ! internal subdomain grid
integer(i_kind),    intent(in   ) :: nymd,nhms
integer(i_kind),    intent(in   ) :: iwhat
integer(i_kind),optional,intent(in ) :: tau   ! time interval in hours

character(len=*),parameter:: myname_=myname//'.get_1State_'
character(len=255) :: fname
integer(i_kind):: nymdb,nhmsb
integer(i_kind):: ierr
integer(i_kind):: tau_

  tau_  = -one
  nymdb = nymd
  nhmsb = nhms
  if (present(tau)) then
     tau_ = tau
  endif
  if (tau_>zero) then
     tau_ = 3600*tau
     call tick (nymdb,nhmsb,tau_)
  endif
  if (tau_>0) then ! read forecast fields
     if (evfsoi_afcst) then
        write(fname,'(a,i3.3,2a)') 'mem',iwhat,'/',trim(GSI_ensprga_fname_tmpl)
     else
        write(fname,'(a,i3.3,2a)') 'mem',iwhat,'/',trim(GSI_ensprgb_fname_tmpl)
     endif
  else             ! read background fields
#ifdef _ALSO_READ_XINC_
     ! the stuff inside this ifdef is for testing only; this code should never
     ! deal w/ perturbations; see geos_pertStateIO for that
     if (iwhat<0) then
       ! in this case, abs(iwhat) refers to outer iteration index (this should
       write(fname,'(3a,i3.3,a)')  trim(GSI_ExpId), '.', 'xinc.', abs(iwhat), '.eta.%y4%m2%d2_%h2%n2z.nc4'
     else
#endif /* _ALSO_READ_XINC_ */
       if (evfsoi_afcst) then
          write(fname,'(a,i3.3,2a)') 'mem',iwhat,'/',trim(GSI_ensana_fname_tmpl)
       else
          write(fname,'(a,i3.3,2a)') 'mem',iwhat,'/',trim(GSI_ensbkg_fname_tmpl)
       endif
#ifdef _ALSO_READ_XINC_
     endif
#endif /* _ALSO_READ_XINC_ */
  endif

  xx=zero
  call gcm2gsi0_(xx,sgrid,ierr,filename=fname,nymd=nymdb,nhms=nhmsb)
        if(ierr/=0) then
          call perr(myname_,'trouble reading guess-like vector')
          call die(myname_)
        endif
end subroutine get_1State_

subroutine get_NState_(xx,sgrid,nymd,nhms,tau)
!
use constants, only:zero
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: assignment(=)
use general_sub2grid_mod, only: sub2grid_info
implicit none
type(gsi_bundle),   intent(inout) :: xx(:)
type(sub2grid_info),intent(in   ) :: sgrid ! internal subdomain grid
integer(i_kind),    intent(in   ) :: nymd,nhms
integer(i_kind),optional,intent(in ) :: tau   ! time interval in hours

character(len=*),parameter:: myname_=myname//'.get_NState_'
character(len=255),allocatable :: fname(:)
integer(i_kind):: nymdb,nhmsb
integer(i_kind):: ierr
integer(i_kind):: tau_
integer nf

  tau_  = -1
  nymdb = nymd
  nhmsb = nhms
  if (present(tau)) then
     tau_ = tau
  endif
  if (tau_>0) then
     tau_ = 3600*tau
     call tick (nymdb,nhmsb,tau_)
  endif
  allocate(fname(size(xx)))
  do nf=1,size(xx)
     if (tau_>0) then ! read forecast fields
        if (evfsoi_afcst) then
           write(fname(nf),'(a,i3.3,2a)') 'mem',nf,'/',trim(GSI_ensprga_fname_tmpl)
        else
           write(fname(nf),'(a,i3.3,2a)') 'mem',nf,'/',trim(GSI_ensprgb_fname_tmpl)
        endif
     else             ! read background fields
        if (evfsoi_afcst) then
           write(fname(nf),'(a,i3.3,2a)') 'mem',nf,'/',trim(GSI_ensana_fname_tmpl)
        else
           write(fname(nf),'(a,i3.3,2a)') 'mem',nf,'/',trim(GSI_ensbkg_fname_tmpl)
        endif
     endif
     xx(nf)=zero
!    call gcm2gsi0_(xx(nf),sgrid,ierr,filename=fname(nf),nymd=nymdb,nhms=nhmsb)
!       if(ierr/=0) then
!         call perr(myname_,'trouble reading guess-like vector')
!         call die(myname_)
!       endif
  enddo

  call gcm2gsi0N_(xx,sgrid,ierr,filename=fname,nymd=nymdb,nhms=nhmsb)
        if(ierr/=0) then
          call perr(myname_,'trouble reading guess-like vector')
          call die(myname_)
        endif

  deallocate(fname)
end subroutine get_NState_

subroutine put_1State_(xx,sgrid,nymd,nhms,iwhat)
!
use constants, only:zero
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: assignment(=)
use general_sub2grid_mod, only: sub2grid_info
implicit none
type(gsi_bundle),   intent(inout) :: xx
type(sub2grid_info),intent(in   ) :: sgrid ! internal subdomain grid
integer(i_kind),    intent(in   ) :: nymd,nhms
integer(i_kind),    intent(in   ) :: iwhat

character(len=*),parameter:: myname_=myname//'.put_1State_'
character(len=255) :: fname
integer(i_kind):: ierr

  write(fname,'(a,i3.3,a)') '%s.gsi',iwhat,'.eta.%y4%m2%d2_%h2%n2z.nc4'

  call gsi2pgcm0_ (nymd,nhms,xx,sgrid,'tlm',ierr,filename=fname)
        if(ierr/=0) then
          call perr(myname_,'trouble writing guess-like vector')
          call die(myname_)
        endif
end subroutine put_1State_
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gsi2pgcm0_:  Convert GSI increments to GCM perturbations
!
! !INTERFACE:

      subroutine gsi2pgcm0_ ( nymd, nhms, xx, sgrid, which, stat, &
                              xp, filename )  ! optionals

! !USES:

      use ESMF

      use GSI_GridCompMod, only: GSI_ExpId
      use GSI_GridCompMod, only: GSI_RefTime

      use mpimod,    only: nxpe,nype
      use general_sub2grid_mod, only: sub2grid_info

      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: assignment(=)

      use m_StrTemplate, only: StrTemplate

      implicit none

! !INPUT PARAMETERS:
                                                                                                                           
      integer(i_kind),    intent(in)    :: nymd   ! date as in YYYYMMDD
      integer(i_kind),    intent(in)    :: nhms   ! time as in HHMMSS
      type(gsi_bundle),   intent(inout) :: xx     ! GSI increment
      type(sub2grid_info),intent(in)    :: sgrid  ! subdomain GSI grid
      character(len=*),   intent(in)    :: which  ! adm or tlm

      character(len=*),optional, intent(in) :: filename ! output filename

! !OUTPUT PARAMETERS:

      type(MAPL_SimpleBundle),optional,intent(out) :: xp
      integer(i_kind),                 intent(out) :: stat

! !DESCRIPTION: Convert GSI increment vector to GEOS-5 perturbation vector
!               (as gsi2pgcm1_, but output GCM perturbation to file)
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  30Sep2007  Todling   Updated interface to putpert.
!  17Jan2010  Todling   Update interface to putpert (add forceflip).
!  15May2010  Todling   Update to use GSI_Bundle
!
!EOP
!-----------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'*gsi2pgcm0_'
     character(len=*), parameter :: Iam = myname_
!    character(len=*), parameter :: only_vars='ps,delp,u,v,tv,sphu,qitot,qltot,ozone'
     type(ESMF_FieldBundle)  :: WBundle
     type(MAPL_SimpleBundle) :: xpert
     type(MAPL_CFIO)         :: cfio
     type(ESMF_Clock)        :: Clock
     type(ESMF_Time)         :: OutTime
     type(ESMF_TimeInterval) :: TimeStep
     type(ESMF_VM)           :: VM
     type(ESMF_Grid)         :: GSI_oGrid ! grid used for write out

     character(len=255) fname,bkgfname,tmpl
     character(len=30) GSIGRIDNAME
     integer(i_kind) thistime(6)
     integer(i_kind) rc,status,ierr
     integer(i_kind) idim,jdim,kdim,i,j,k,i_dp,i_ps
     type(LatLonGridFactory) :: ll_factory

     stat = 0

!    Define working bundle on a GEOS-5 (oriented) grid
!    -------------------------------------------------
     call ESMF_VMGetCurrent(vm=VM, rc=STATUS)
     VERIFY_(STATUS)
     call MAPL_DefGridName (sgrid%nlon,sgrid%nlat,GSIGRIDNAME,MAPL_am_I_root())
     
     ll_factory = LatLonGridFactory(grid_name=GSIGRIDNAME, nx=nxpe, ny=nype, &
                 im_world=sgrid%nlon, jm_world=sgrid%nlat, lm = nsig, &
                 pole='PC', dateline='DC', rc=status)
     VERIFY_(status)
     GSI_oGrid = grid_manager%make_grid(ll_factory,rc=status)
     VERIFY_(STATUS)

     WBundle = ESMF_FieldBundleCreate ( name='Work Bundle', rc=status )
               VERIFY_(status)
     call ESMF_FieldBundleSet ( WBundle, grid=GSI_oGrid, rc=status )
     VERIFY_(status)

!    Set file to be read
!    -------------------
!    fname = fnpert
!    if(present(filename)) fname = trim(filename)

!    Get memory for output bundle by reading inital bkg file
!    NOTE: 1. this needs revision since not general enough
!          2. filanme template wired-in for now
!    -------------------------------------------------------
     thistime(1) =     nymd/10000
     thistime(2) = mod(nymd,10000)/100
     thistime(3) = mod(nymd,100)
     thistime(4) =     nhms/10000
     thistime(5) = mod(nhms,10000)/100
     thistime(6) = mod(nhms,100)

!   initialize current time
!   -----------------------
     call ESMF_TimeSet( OutTime, YY =  thistime(1), &
                                 MM =  thistime(2), &
                                 DD =  thistime(3), &
                                 H  =  thistime(4), &
                                 M  =  thistime(5), &
                                 S  =  thistime(6), rc=status ); VERIFY_(STATUS)

!    Just to get room for ESMF bundle ...
!    ------------------------------------
     call StrTemplate ( bkgfname, GSI_bkg_fname_tmpl, 'GRADS', & 
                        xid=trim(GSI_ExpId), nymd=nymd, nhms=nhms, &
                        stat=status )
          VERIFY_(STATUS)

     call MAPL_CFIORead  ( trim(bkgfname), OutTime, WBundle, &
                           TIME_IS_CYCLIC=.false., verbose=.true., &
                                                noread=.true., rc=status )
!                          only_vars=only_vars, noread=.true., rc=status )
          VERIFY_(status)

!    Link Bundle to Simple-Bundle
!    ----------------------------
     xpert = MAPL_SimpleBundleCreate ( WBundle, rc=status )
             VERIFY_(status)

!    Convert to GSI perturbation vector
!    ----------------------------------
     call gsi2pgcm1_ ( xx, xpert, sgrid, GSI_oGrid, which, status )
          VERIFY_(status)

!    Create clock set at reference time
!    ----------------------------------
     call ESMF_TimeIntervalSet( TimeStep, h=0, m=0, s=0, rc=status )
          VERIFY_(STATUS)
     Clock = ESMF_ClockCreate ( name="statewrite", timeStep=TimeStep, &
                                startTime=OutTime, rc=status )
             VERIFY_(STATUS)

!    Write out perturbation
!    ----------------------
     mycount = mycount + 1
     if (present(filename)) then
       tmpl = trim(filename)
     else
       write(tmpl,'(4a,i3.3,1a)') trim(GSI_ExpId), '.', trim(fnxgsi), '_', mycount, '.%y4%m2%d2_%h2%n2z.nc4'
     endif
  
     call StrTemplate ( fname, tmpl, 'GRADS', xid=GSI_ExpId, nymd=nymd, nhms=nhms, &
                        stat=status )
          VERIFY_(STATUS)

!    Before writing to the file, we create a CFIO object
!    Here we define a file with reduced resolution
!    ---------------------------------------------------
     call MAPL_CFIOCreate ( cfio, fname, clock, WBundle, rc=status )
          VERIFY_(STATUS)

     call MAPL_CFIOWrite ( cfio, Clock, WBundle, verbose=.true., rc=status )
       if(status/=0)then
           stat = 90
           if(mype==ROOT) print*, trim(myname_), ': Error retrieving perturbation'
           return
       endif

     call MAPL_cfioDestroy ( cfio )

     if (present(xp)) then
!        Should copy simple bundle xpert to xp
         call die(myname_, 'feature not ready ...')
     endif

!    Release GCM perturbation vector
!    -------------------------------
!    call MAPL_SimpleBundleDestroy ( xpert, rc=STATUS )  ! _RT is this a nullify?? it better be
!         VERIFY_(STATUS)
     call ESMF_FieldBundleDestroy ( WBundle, rc=STATUS )
          VERIFY_(STATUS)


     end subroutine gsi2pgcm0_

!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1    !
!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gsi2pgcm1_:  Convert GSI increments to GCM perturbations
!
! !INTERFACE:

      subroutine gsi2pgcm1_ ( xx, xpert, sg, GSI_xGrid, which, stat )

! !USES:

      use ESMF
      use general_sub2grid_mod, only: sub2grid_info
      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: gsi_bundlegetpointer
      use gsi_bundlemod, only: assignment(=)

      implicit none

! !INPUT PARAMETERS:

      type(gsi_bundle),   intent(inout) :: xx    ! GSI increment vector
      character(len=*),   intent(in)    :: which ! adm or tlm
      type(ESMF_Grid),    intent(in)    :: GSI_xGrid

! !OUTPUT PARAMETERS:

      type(MAPL_SimpleBundle),intent(inout) :: xpert ! GCM perturbation vector
      type(sub2grid_info),    intent(in)    :: sg    ! subdomain GSI grid
      integer(i_kind),        intent(out)   :: stat  ! return error code

! !DESCRIPTION: Converts GSI increments in to ADM/TLM perturbations.
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  21Sep2007  Todling   Handles for O3 and CW.
!  19Nov2008  Todling   Update to use gsi-3d pressure instead of ps.
!  15May2010  Todling   Update to use GSI_Bundle
!  21Feb2011  Todling   Adapt to work with MAPL-SimpleBundle
!  28Oct2013  Todling   Rename p3d to prse
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*gsi2pgcm_'

      real(r_kind), allocatable, dimension(:,:,:) :: sub_tv,sub_u,sub_v,sub_q,sub_delp
      real(r_kind), allocatable, dimension(:,:,:) :: sub_oz,sub_cw
      real(r_kind), allocatable, dimension(:,:,:) :: sub_qi,sub_ql,sub_qr,sub_qs
      real(r_kind), allocatable, dimension(:,:)   :: sub_ps

      integer(i_kind)  i,j,k,kk,ijk,ij
      integer(i_kind)  i_u,i_v,i_t,i_q,i_oz,i_cw,i_prse,i_p,i_dp
      integer(i_kind)  i_qi,i_ql,i_qr,i_qs
      integer(i_kind)  ierr,istatus,rc,status
      character(len=255) :: whatin

      stat = 0

!     Initializes internal arrays
!     ---------------------------
      allocate (sub_tv  (sg%lat2,sg%lon2,nsig), sub_u (sg%lat2,sg%lon2,nsig),&
                sub_v   (sg%lat2,sg%lon2,nsig), sub_q (sg%lat2,sg%lon2,nsig),&
                sub_oz  (sg%lat2,sg%lon2,nsig), & 
                sub_delp(sg%lat2,sg%lon2,nsig), sub_ps(sg%lat2,sg%lon2), stat=ierr )
        if ( ierr/=0 ) then
            stat = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_)'
            return
        end if

!     Get poiners to GSI state-vector
!     -------------------------------
      ierr=0
      call gsi_bundlegetpointer(xx,'sf',  i_u,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'vp',  i_v,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'t' ,  i_t,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'q' ,  i_q,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'oz',  i_oz , istatus);ierr=ierr+istatus
!     call gsi_bundlegetpointer(xx,'prse',i_prse,istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'ps',  i_p,   istatus);ierr=ierr+istatus
      if(ierr/=0) then
         write(6,*) myname_, ': trouble getting pointers'
         call stop2(999)
      endif

!     Fill in subdomain arrays
!     ------------------------
      do k=1,nsig
         do j=1,sg%lon2
            do i=1,sg%lat2
               sub_u (i,j,k) = xx%r3(i_u)%q(i,j,k)
               sub_v (i,j,k) = xx%r3(i_v)%q(i,j,k)
               sub_tv(i,j,k) = xx%r3(i_t)%q(i,j,k)
               sub_q (i,j,k) = xx%r3(i_q)%q(i,j,k)
               sub_oz(i,j,k) = xx%r3(i_oz)%q(i,j,k) / PPMV2GpG
            enddo
         enddo
      enddo

      call gsi_bundlegetpointer(xx,'cw',  i_cw , istatus)
      if (i_cw>0) then
         allocate(sub_cw(sg%lat2,sg%lon2,nsig), stat=ierr )
         if ( ierr/=0 ) then
             stat = 91
             if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_cw)'
             return
         end if
         do k=1,nsig
            do j=1,sg%lon2
               do i=1,sg%lat2
                  sub_cw(i,j,k) = xx%r3(i_cw)%q(i,j,k)
               enddo
            enddo
         enddo
      endif

      call gsi_bundlegetpointer(xx,'qi',  i_qi , istatus)
      if (i_qi>0) then
         allocate(sub_qi(sg%lat2,sg%lon2,nsig), stat=ierr )
         if ( ierr/=0 ) then
             stat = 91
             if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_qi)'
             return
         end if
         do k=1,nsig
            do j=1,sg%lon2
               do i=1,sg%lat2
                  sub_qi(i,j,k) = xx%r3(i_qi)%q(i,j,k)
               enddo
            enddo
         enddo
      endif

      call gsi_bundlegetpointer(xx,'ql',  i_ql , istatus)
      if (i_ql>0) then
         allocate(sub_ql(sg%lat2,sg%lon2,nsig), stat=ierr )
         if ( ierr/=0 ) then
             stat = 91
             if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_ql)'
             return
         end if
         do k=1,nsig
            do j=1,sg%lon2
               do i=1,sg%lat2
                  sub_ql(i,j,k) = xx%r3(i_ql)%q(i,j,k)
               enddo
            enddo
         enddo
      endif

      call gsi_bundlegetpointer(xx,'qr',  i_qr , istatus)
      if (i_qr>0) then
         allocate(sub_qr(sg%lat2,sg%lon2,nsig), stat=ierr )
         if ( ierr/=0 ) then
             stat = 91
             if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_qr)'
             return
         end if
         do k=1,nsig
            do j=1,sg%lon2
               do i=1,sg%lat2
                  sub_qr(i,j,k) = xx%r3(i_qr)%q(i,j,k)
               enddo
            enddo
         enddo
      endif

      call gsi_bundlegetpointer(xx,'qs',  i_qs , istatus)
      if (i_qs>0) then
         allocate(sub_qs(sg%lat2,sg%lon2,nsig), stat=ierr )
         if ( ierr/=0 ) then
             stat = 91
             if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_qs)'
             return
         end if
         do k=1,nsig
            do j=1,sg%lon2
               do i=1,sg%lat2
                  sub_qs(i,j,k) = xx%r3(i_qs)%q(i,j,k)
               enddo
            enddo
         enddo
      endif

      do j=1,sg%lon2
         do i=1,sg%lat2
            sub_ps(i,j) = xx%r2(i_p)%q(i,j) * Pa_per_kPa
         enddo
      enddo
      sub_delp=zero
      do k=1,nsig
         sub_delp(:,:,k) = sub_ps ! just a trick to use gsi2pert_ w/o interface change
      enddo

!     Gather from GSI subdomains/Scatter to GCM
!     -----------------------------------------
      i_p = MAPL_SimpleBundleGetIndex ( xpert, 'ps'  , 2, rc=status )
      i_dp= MAPL_SimpleBundleGetIndex ( xpert, 'delp', 3, rc=status )
      i_u = MAPL_SimpleBundleGetIndex ( xpert, 'u'   , 3, rc=status )
      i_v = MAPL_SimpleBundleGetIndex ( xpert, 'v'   , 3, rc=status )
      i_t = MAPL_SimpleBundleGetIndex ( xpert, 'tv'  , 3, rc=status )
      i_q = MAPL_SimpleBundleGetIndex ( xpert, 'sphu', 3, rc=status )
      i_oz= MAPL_SimpleBundleGetIndex ( xpert, 'ozone',3, rc=status )
      call gsi2pert_ ( sub_u,    xpert%r3(i_u)%qr4,  ierr )
      call gsi2pert_ ( sub_v,    xpert%r3(i_v)%qr4,  ierr )
      call gsi2pert_ ( sub_tv,   xpert%r3(i_t)%qr4,  ierr )
      call gsi2pert_ ( sub_q ,   xpert%r3(i_q)%qr4,  ierr )
      call gsi2pert_ ( sub_oz,   xpert%r3(i_oz)%qr4, ierr )
      call gsi2pert_ ( sub_delp, xpert%r3(i_dp)%qr4, ierr )

      if (i_cw>0) then
         i_cw= MAPL_SimpleBundleGetIndex ( xpert, 'qltot',3, rc=status )
         call gsi2pert_ ( sub_cw,   xpert%r3(i_cw)%qr4, ierr )
      endif
      if (i_ql>0) then ! in case, both cw and ql are present ...
                       ! too bad: ql will overwrite cw
         i_ql= MAPL_SimpleBundleGetIndex ( xpert, 'qltot',3, rc=status )
         call gsi2pert_ ( sub_ql,   xpert%r3(i_ql)%qr4, ierr )
      endif
      if (i_qi>0) then
         i_qi= MAPL_SimpleBundleGetIndex ( xpert, 'qitot',3, rc=status )
         call gsi2pert_ ( sub_qi,   xpert%r3(i_qi)%qr4, ierr )
      endif
      if (i_qr>0) then
         i_qr= MAPL_SimpleBundleGetIndex ( xpert, 'qrtot',3, rc=status )
         call gsi2pert_ ( sub_qr,   xpert%r3(i_qr)%qr4, ierr )
      endif
      if (i_qs>0) then
         i_qs= MAPL_SimpleBundleGetIndex ( xpert, 'qstot',3, rc=status )
         call gsi2pert_ ( sub_qs,   xpert%r3(i_qs)%qr4, ierr )
      endif

!     Build surface pressure perturbation - output purposes
!     -----------------------------------------------------
      xpert%r2(i_p )%qr4 = xpert%r3(i_dp)%qr4(:,:,1) ! it has been flipped
      do k=1,nsig ! the following is not correct for full fields
         kk = nsig-k+1
         xpert%r3(i_dp)%qr4(:,:,kk)=(bk5(k+1)-bk5(k))*xpert%r2(i_p)%qr4
      end do

      deallocate (sub_tv, sub_u, sub_v, sub_q, sub_delp, sub_ps, &
                  sub_oz, stat=ierr )
        if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_)'
            return
        end if

      if (allocated(sub_cw)) deallocate(sub_cw)
      if (allocated(sub_qi)) deallocate(sub_qi)
      if (allocated(sub_ql)) deallocate(sub_ql)
      if (allocated(sub_qr)) deallocate(sub_qr)
      if (allocated(sub_qs)) deallocate(sub_qs)

      CONTAINS

      subroutine gsi2pert_ ( sub, fld, stat_ )

      real(r_kind),   intent(in)  :: sub(:,:,:)
      real(4),      intent(inout) :: fld(:,:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*gsi2pert_'

      real(r_kind), allocatable :: fldsm(:,:)
      real(r_kind), allocatable :: work3d(:,:,:)     ! auxliar 3d array
      real(r_kind), allocatable :: work(:)

      integer(i_kind) imr,jnp,nl
      integer(i_kind) i,j,k,mm1,ierr

      mm1 = mype+1 
      stat_= 0
      imr  =sg%nlon
      jnp  =sg%nlat
      nl   =nsig

      allocate ( work3d(sg%nlon,sg%nlat,nsig), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work3d)'
            return
        end if
      allocate ( work(max(sg%iglobal,sg%itotsub)), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work)'
            return
        end if
      allocate ( fldsm(sg%lat1*sg%lon1,nsig), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(fldsm)'
            return
        end if

!     Strip off boundary points from subdomains
!     -----------------------------------------
      call strip_(sg,sub,fldsm,nsig)

!     Gather GSI perturbations to root processor
!     ------------------------------------------
      do k=1,nsig
         call mpi_gatherv(fldsm(1,k),sg%ijn(mm1),mpi_rtype,&
              work,sg%ijn,sg%displs_g,mpi_rtype,&
              ROOT,mpi_comm_world,ierr)
         if (mype==ROOT) then
            call reorder12_(sg,work,work3d(:,:,k))
         endif
      end do

      deallocate ( fldsm, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(fldsm)'
            return
        end if
      deallocate ( work, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work)'
            return
        end if

!     Flip horizontal and vertical
!     ----------------------------
      if ( mype==ROOT ) then
           call hflip3_ ( work3d, imr,jnp,nl )
           call SwapV_ ( work3d )
      endif

!     Scatter perturbations to GCM decomposition
!     ------------------------------------------
      call myscatter_ ( which, GSI_xGrid, imr, jnp, nl, work3d, fld )

!     Swap work memory
!     ----------------
      deallocate ( work3d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work3d)'
            return
        end if

      end subroutine gsi2pert_

      end subroutine gsi2pgcm1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gcm2gsi0_:  Convert gcm vector to gsi vector
!
! !INTERFACE:

      subroutine gcm2gsi0_ ( xx, sgrid, stat, &
                             filename, skiptraj, nymd, nhms )

! !USES:

      use kinds, only: i_kind,r_kind
      use constants, only: zero
      use m_StrTemplate, only: StrTemplate

      use ESMF

      use GSI_GridCompMod, only: GSI_RefTime
      use GSI_GridCompMod, only: GSI_ExpId

      use mpimod,    only: nxpe,nype
      use general_sub2grid_mod, only: sub2grid_info
      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: assignment(=)

      implicit none

! !INPUT PARAMETERS:

      character(len=*), optional, intent(in)  :: filename ! name of file w/ perturbation
      logical,          optional, intent(in)  :: skiptraj ! allows skip of trajectory/ics 

! !OUTPUT PARAMETERS:
			
      type(gsi_bundle),         intent(inout) :: xx       ! GSI increment
      type(sub2grid_info),      intent(in)    :: sgrid    ! subdomain GSI grid

      integer(i_kind),            intent(out) :: stat

      integer(i_kind),optional,   intent(in ) :: nymd
      integer(i_kind),optional,   intent(in ) :: nhms

! !DESCRIPTION: Convert GEOS-5 perturbation vector to GSI increment vector
!               (as gcm2gsi1_, but reads GCM perturbation from file)
!
! !REMARKS:
!
!  31Oct2007  Todling   De-activated interpolation capability; not fully coded yet
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  17Jul2007  Todling   Add ability to read in perturbation at diff resolution
!  17Jul2008  Todling   Add filename as optional argument
!  15May2010  Todling   Update to use GSI_Bundle
!  20Feb2011  Todling   Adapt to work with MAPL-SimpleBundle
!  01Dec2014  Todling   Add do-parallel per Ben Auer
!
!EOP
!-----------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'*gcm2gsi0_'
     character(len=*), parameter :: Iam = myname_
     type(ESMF_FieldBundle)  :: WBundle
     type(MAPL_SimpleBundle) :: xpert
     type(ESMF_Time)         :: OutTime
     type(ESMF_VM)           :: VM
     type(ESMF_Grid)         :: GSI_oGrid ! grid used for write out
     type(LatLonGridFactory) :: ll_factory

     character(len=255) :: fname
     character(len=30)  :: GSIGRIDNAME
     integer(i_kind) myimr,myjnp,mynl,mync
     integer(i_kind) ierr,yy,mm,dd,hh,mn,sec
     integer(i_kind) thistime(6)
     integer         rc,status
     real(r_kind)    dmodel,dgsi

     stat = 0
     xx   = zero

     call timer_ini ('gcm2gsi0_')

!    Set file to be read
!    -------------------
     if(present(nymd) .and. present(nhms) ) then
        call StrTemplate ( fname, filename, 'GRADS', & 
                           xid=trim(GSI_ExpId), nymd=nymd, nhms=nhms, &
                           stat=status )
                           VERIFY_(STATUS)
    else
        print*, 'gcm2gsi: parameters are really not optional'
        call stop2(999)
    endif

!   Set ESMF clock to proper time
!   -----------------------------
    thistime(1) =     nymd/10000
    thistime(2) = mod(nymd,10000)/100
    thistime(3) = mod(nymd,100)
    thistime(4) =     nhms/10000
    thistime(5) = mod(nhms,10000)/100
    thistime(6) = mod(nhms,100)

!   initialize current time
!   -----------------------
     call ESMF_TimeSet( OutTime, YY =  thistime(1), &
                                 MM =  thistime(2), &
                                 DD =  thistime(3), &
                                 H  =  thistime(4), &
                                 M  =  thistime(5), &
                                 S  =  thistime(6), rc=status ); VERIFY_(STATUS)

!    Define working bundle on a GEOS-5 (oriented) grid
!    -------------------------------------------------
     call ESMF_VMGetCurrent(vm=VM, rc=STATUS)
     VERIFY_(STATUS)
     call MAPL_DefGridName (sgrid%nlon,sgrid%nlat,GSIGRIDNAME,MAPL_am_I_root())

     ll_factory = LatLonGridFactory(grid_name=GSIGRIDNAME, nx=nxpe, ny=nype, &
                 im_world=sgrid%nlon, jm_world=sgrid%nlat, lm = nsig, &
                 pole='PC', dateline='DC', rc=status)
     VERIFY_(status)
     GSI_oGrid = grid_manager%make_grid(ll_factory,rc=status)
     VERIFY_(STATUS)

     WBundle = ESMF_FieldBundleCreate ( name='Work Bundle', rc=status )
               VERIFY_(status)
     call ESMF_FieldBundleSet ( WBundle, grid=GSI_oGrid, rc=status )
     VERIFY_(status)

!    Read perturbation from file into ESMF Bundle
!    --------------------------------------------
     call timer_ini('g2g_read')
     call MAPL_CFIORead  ( fname, OutTime, WBundle, & !only_vars=only_vars, &
                           TIME_IS_CYCLIC=.false., verbose=.false., &
                           doParallel=.true., gsiMode=.true., rc=status )
          VERIFY_(status)
     call timer_fnl('g2g_read')

!    Link Bundle to Simple-Bundle
!    ----------------------------
     xpert = MAPL_SimpleBundleCreate ( WBundle, rc=status )
             VERIFY_(status)

!    Need to convert ESMF time to conventional time
!    ----------------------------------------------
!    if (present(nymd_in)) then
!        call ESMF_TimeGet(GSI_RefTime, yy=YY, mm=MM, dd=DD, rc=status)
!        nymd_in = 10000*yy + 100*mm + dd
!    endif
!    if (present(nhms_in)) then
!        call ESMF_TimeGet(GSI_RefTime, h=HH, m=MN, s=SEC, rc=status)
!        nhms_in = 10000*hh + 100*mn + sec
!    endif

!    Simply convert in GSI fields
!    ----------------------------
     call gcm2gsi1_ ( xpert, xx, sgrid, GSI_oGrid, stat )

     call MAPL_SimpleBundleDestroy ( xpert, rc=STATUS )  ! _RT is this a nullify?? it better be
          VERIFY_(STATUS)
     call ESMF_FieldBundleDestroy ( WBundle, rc=STATUS )
          VERIFY_(STATUS)

     call timer_fnl ('gcm2gsi0_')
     end subroutine gcm2gsi0_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gcm2gsi0N_:  Convert N gcm vectors to N gsi vectors
!
! !INTERFACE:

      subroutine gcm2gsi0N_ ( xx, sgrid, stat, &
                             filename, skiptraj, nymd, nhms )

! !USES:

      use kinds, only: i_kind,r_kind
      use constants, only: zero
      use m_StrTemplate, only: StrTemplate

      use ESMF

      use GSI_GridCompMod, only: GSI_RefTime
      use GSI_GridCompMod, only: GSI_ExpId

      use mpimod,    only: nxpe,nype
      use general_sub2grid_mod, only: sub2grid_info
      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: assignment(=)

      implicit none

! !INPUT PARAMETERS:

      character(len=*), optional, intent(in)  :: filename(:) ! name of file w/ perturbation
      logical,          optional, intent(in)  :: skiptraj    ! allows skip of trajectory/ics 

! !OUTPUT PARAMETERS:
			
      type(gsi_bundle),         intent(inout) :: xx(:)    ! GSI increment
      type(sub2grid_info),      intent(in)    :: sgrid    ! subdomain GSI grid

      integer(i_kind),            intent(out) :: stat

      integer(i_kind),optional,   intent(in ) :: nymd
      integer(i_kind),optional,   intent(in ) :: nhms

! !DESCRIPTION: Convert GEOS-5 perturbation vector to GSI increment vector
!               (as gcm2gsi1_, but reads GCM perturbation from file)
!
! !REMARKS:
!
!  31Oct2007  Todling   De-activated interpolation capability; not fully coded yet
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  17Jul2007  Todling   Add ability to read in perturbation at diff resolution
!  17Jul2008  Todling   Add filename as optional argument
!  15May2010  Todling   Update to use GSI_Bundle
!  20Feb2011  Todling   Adapt to work with MAPL-SimpleBundle
!  01Dec2014  Todling   Add do-parallel per Ben Auer
!  19Jun2017  Todling   Add to hook-up with MAPL multi-bundle read of B. Auer
!
!EOP
!-----------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'*gcm2gsi0N_'
     character(len=*), parameter :: Iam = myname_
     type(ESMF_FieldBundle),pointer :: WBundle(:)
     type(MAPL_SimpleBundle) :: xpert
     type(ESMF_Time)         :: OutTime
     type(ESMF_VM)           :: VM
     type(ESMF_Grid)         :: GSI_oGrid ! grid used for write out

     character(len=255),allocatable :: fname(:)
     character(len=30)  :: GSIGRIDNAME
     integer(i_kind) myimr,myjnp,mynl,mync
     integer(i_kind) ierr,yy,mm,dd,hh,mn,sec
     integer(i_kind) thistime(6)
     integer         nf,rc,status
     real(r_kind)    dmodel,dgsi
     type(LatLonGridFactory) :: ll_factory

     stat = 0

     call timer_ini ('gcm2gsi0N_')

!    Set file to be read
!    -------------------
     allocate(fname(size(xx)))
     if(present(nymd) .and. present(nhms) ) then
        do nf=1,size(xx)
           call StrTemplate ( fname(nf), filename(nf), 'GRADS', & 
                              xid=trim(GSI_ExpId), nymd=nymd, nhms=nhms, &
                              stat=status )
                              VERIFY_(STATUS)
        enddo
    else
        print*, 'gcm2gsi: parameters are really not optional'
        call stop2(999)
    endif

!   Set ESMF clock to proper time
!   -----------------------------
    thistime(1) =     nymd/10000
    thistime(2) = mod(nymd,10000)/100
    thistime(3) = mod(nymd,100)
    thistime(4) =     nhms/10000
    thistime(5) = mod(nhms,10000)/100
    thistime(6) = mod(nhms,100)

!   initialize current time
!   -----------------------
     call ESMF_TimeSet( OutTime, YY =  thistime(1), &
                                 MM =  thistime(2), &
                                 DD =  thistime(3), &
                                 H  =  thistime(4), &
                                 M  =  thistime(5), &
                                 S  =  thistime(6), rc=status ); VERIFY_(STATUS)

!    Define working bundle on a GEOS-5 (oriented) grid
!    -------------------------------------------------
     call ESMF_VMGetCurrent(vm=VM, rc=STATUS)
     VERIFY_(STATUS)
     call MAPL_DefGridName (sgrid%nlon,sgrid%nlat,GSIGRIDNAME,MAPL_am_I_root())
     ll_factory = LatLonGridFactory(grid_name=GSIGRIDNAME, nx=nxpe, ny=nype, &
                 im_world=sgrid%nlon, jm_world=sgrid%nlat, lm = nsig, &
                 pole='PC', dateline='DC', rc=status)
     VERIFY_(status)
     GSI_oGrid = grid_manager%make_grid(ll_factory,rc=status)
     VERIFY_(STATUS)

     allocate(WBundle(size(xx)))
     do nf=1,size(xx)
        WBundle(nf) = ESMF_FieldBundleCreate ( name='Work Bundle', rc=status )
                      VERIFY_(status)
        call ESMF_FieldBundleSet ( WBundle(nf), grid=GSI_oGrid, rc=status )
        VERIFY_(status)
     enddo

!    Read perturbation from file into ESMF Bundle
!    --------------------------------------------
     call timer_ini('g2g_read')
     call MAPL_CFIOReadParallel(WBundle,fname,OutTime,blocksize=GSI_ensread_blocksize,gsiMode=.true.,rc=status)
     if (status/=0) then
        call die(myname_,': failed reading ensemble, aborting ...',99)
     endif
     call timer_fnl('g2g_read')

!    Link Bundle to Simple-Bundle
!    ----------------------------
     do nf=1,size(xx)
        xpert = MAPL_SimpleBundleCreate ( WBundle(nf), rc=status )
                VERIFY_(status)

!       Simply convert in GSI fields
!       ----------------------------
        xx(nf) = zero
        call gcm2gsi1_ ( xpert, xx(nf), sgrid, GSI_oGrid, stat )
     enddo

     call MAPL_SimpleBundleDestroy ( xpert, rc=STATUS )  ! _RT is this a nullify?? it better be
          VERIFY_(STATUS)
     do nf=size(xx),1,-1
        call ESMF_FieldBundleDestroy ( WBundle(nf), rc=STATUS )
             VERIFY_(STATUS)
     enddo
     deallocate(WBundle)
     deallocate(fname)

     call timer_fnl ('gcm2gsi0N_')
     end subroutine gcm2gsi0N_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gcm2gsi1_:  Convert gcm-adm/tlm vector to gsi vector
!
! !INTERFACE:

      subroutine gcm2gsi1_ ( xpert, xx, sg, GSI_xGrid, stat )

! !USES:

      use ESMF
      use constants, only: one,zero
      use gsi_bundlemod,only: gsi_bundle     ! GSI state vector
      use gsi_bundlemod,only: gsi_bundlegetpointer
      use gsi_bundlemod,only: assignment(=)
      use general_sub2grid_mod, only: sub2grid_info
 
      implicit none

! !INPUT PARAMETERS:

      type(MAPL_SimpleBundle), intent(inout)  :: xpert  ! GCM perturbation vector 
      type(sub2grid_info),     intent(in)     :: sg     ! subdomain grid
      type(ESMF_Grid),         intent(in)     :: GSI_xGrid

! !OUTPUT PARAMETERS:

      type(gsi_bundle),         intent(inout) :: xx     ! GSI increment

      integer(i_kind),            intent(out) :: stat

! !DESCRIPTION: Convert GEOS-5 perturbation vector to GSI increment vector
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  21Sep2007  Todling   Handles for O3 and CW.
!  22Feb2008  Todling   Handle for forecast (error) gradient vector.
!  19Nov2008  Todling   Update to use gsi-3d pressure instead of ps.
!  15May2010  Todling   Update to use GSI_Bundle
!  20Feb2011  Todling   Adapt to work with MAPL-SimpleBundle
!  30Nov2014  Todling   Add some variable exceptions (qi/ql/ph/ts)
!  06May2020  Todling   Shuffled pointer check and alloc for increased flexibility
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*gcm2gsi1_'
      character(len=*), parameter :: which='tlm'

      real(r_kind),  allocatable, dimension(:,:,:) :: sub_u,sub_v,sub_q,sub_tv
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_oz,sub_ci,sub_cl
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_cr,sub_cs
      real(r_kind),  allocatable, dimension(:,:)   :: sub_ps,sub_ph,sub_ts
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_3d

      integer(i_kind) i,j,k,ij,ijk,id,ng_d,ng_s,ni,nj
      integer(i_kind) i_u,i_v,i_t,i_q,i_oz,i_cw,i_ps,i_ts,i_tsen,i_ph
      integer(i_kind) i_ci,i_cl,i_cr,i_cs
      integer(i_kind) j_cr,j_cs
      integer(i_kind) ierr,istatus,status

      stat = 0

      if ( sg%nsig/=nsig ) then
          stat = 90
          if(mype==ROOT) print*, trim(myname_), ': unnacceptable (sg%nsig,nsig) = ', sg%nsig, nsig
          return
      end if

!     Gather from GCM/Scatter to GSI subdomains
!     -----------------------------------------
      i_u  = MAPL_SimpleBundleGetIndex ( xpert, 'u'    , 3, rc=status )
      i_v  = MAPL_SimpleBundleGetIndex ( xpert, 'v'    , 3, rc=status )
      i_t  = MAPL_SimpleBundleGetIndex ( xpert, 'tv'   , 3, rc=status )
      i_q  = MAPL_SimpleBundleGetIndex ( xpert, 'sphu' , 3, rc=status )
      i_oz = MAPL_SimpleBundleGetIndex ( xpert, 'ozone', 3, rc=status )
      i_cl = MAPL_SimpleBundleGetIndex ( xpert, 'qltot', 3, rc=status )
      i_ci = MAPL_SimpleBundleGetIndex ( xpert, 'qitot', 3, rc=status )
      i_ps = MAPL_SimpleBundleGetIndex ( xpert, 'ps'   , 2, rc=status )
      i_ts = MAPL_SimpleBundleGetIndex ( xpert, 'ts'   , 2, rc=status )
      i_ph = MAPL_SimpleBundleGetIndex ( xpert, 'phis' , 2, rc=status )
      if (i_u>0.and.i_v>0) then
         allocate(sub_u(sg%lat2,sg%lon2,nsig), stat=ierr )
         call pert2gsi_ ( xpert%r3(i_u)%qr4 , sub_u  , ierr )
         allocate(sub_v(sg%lat2,sg%lon2,nsig), stat=ierr )
         call pert2gsi_ ( xpert%r3(i_v)%qr4 , sub_v  , ierr )
      endif
      if (i_t>0) then
         allocate(sub_tv(sg%lat2,sg%lon2,nsig), stat=ierr )
         call pert2gsi_ ( xpert%r3(i_t)%qr4 , sub_tv , ierr )
      endif
      if (i_q>0) then
         allocate(sub_q(sg%lat2,sg%lon2,nsig), stat=ierr )
         call pert2gsi_ ( xpert%r3(i_q)%qr4 , sub_q  , ierr )
      endif
      if (i_oz>0) then
         allocate(sub_oz(sg%lat2,sg%lon2,nsig), stat=ierr )
         call pert2gsi_ ( xpert%r3(i_oz)%qr4, sub_oz , ierr )
      endif
      if (i_cl>0) then
         allocate(sub_cl(sg%lat2,sg%lon2,nsig), stat=ierr )
         call pert2gsi_ ( xpert%r3(i_cl)%qr4, sub_cl , ierr )
      endif
      if (i_ci>0) then
         allocate(sub_ci(sg%lat2,sg%lon2,nsig), stat=ierr )
         call pert2gsi_ ( xpert%r3(i_ci)%qr4, sub_ci , ierr )
      endif
      if (i_ps>0) then
         allocate(sub_ps(sg%lat2,sg%lon2))
         call pert2gsi2d_ ( xpert%r2(i_ps)%qr4, sub_ps, ierr )
      endif
      if (i_ph>0) then
         allocate(sub_ph(sg%lat2,sg%lon2))
         call pert2gsi2d_ ( xpert%r2(i_ph)%qr4, sub_ph, ierr )
      endif
      if (i_ts>0) then
         allocate(sub_ts(sg%lat2,sg%lon2))
         call pert2gsi2d_ ( xpert%r2(i_ts)%qr4, sub_ts, ierr )
      endif
      if ( ierr/=0 ) then
          stat = 99
          if(mype==ROOT) print*, trim(myname_), ': unfinished convertion ...'
          return
      end if

!     Get pointers to GSI state-vector
!     --------------------------------
      ierr=0
      if (allocated(sub_u) .and. allocated(sub_v)) then
         call gsi_bundlegetpointer(xx,'sf' ,  i_u,  istatus);ierr=ierr+istatus
         call gsi_bundlegetpointer(xx,'vp' ,  i_v,  istatus);ierr=ierr+istatus
         if (i_u>0 .and. i_v>0) then
            do k=1,nsig
               do j=1,sg%lon2
                  do i=1,sg%lat2
                     xx%r3(i_u)%q(i,j,k) = sub_u (i,j,k)
                     xx%r3(i_v)%q(i,j,k) = sub_v (i,j,k)
                  enddo
               enddo
            enddo
         endif
         deallocate(sub_v)
         deallocate(sub_u)
      endif
      if (allocated(sub_tv)) then
         call gsi_bundlegetpointer(xx,'t'  ,  i_t,  istatus);ierr=ierr+istatus
         if (i_t>0) then
            do k=1,nsig
               do j=1,sg%lon2
                  do i=1,sg%lat2
                     xx%r3(i_t)%q(i,j,k) = sub_tv(i,j,k)
                  enddo
               enddo
            enddo
         endif
         deallocate(sub_tv)
      endif
      if (allocated(sub_q)) then
         call gsi_bundlegetpointer(xx,'q',  i_q,  istatus);ierr=ierr+istatus
         if (i_q>0) then
            do k=1,nsig
               do j=1,sg%lon2
                  do i=1,sg%lat2
                     xx%r3(i_q)%q(i,j,k) = sub_q (i,j,k)
                  enddo
               enddo
            enddo
         endif
         deallocate(sub_q)
      endif
      if (allocated(sub_oz)) then
         call gsi_bundlegetpointer(xx,'oz' ,  i_oz, istatus);ierr=ierr+istatus
         if (i_oz>0) then
            do k=1,nsig
               do j=1,sg%lon2
                  do i=1,sg%lat2
                     xx%r3(i_oz)%q(i,j,k)= sub_oz(i,j,k) * PPMV2GpG
                  enddo
               enddo
            enddo
         endif
         deallocate(sub_oz)
      endif

      if (allocated(sub_ci).and.allocated(sub_cl)) then
         call gsi_bundlegetpointer(xx,'cw', i_cw, istatus)
         if (i_cw>0 .and. i_ci>0 .and. i_cl>0) then
            do k=1,nsig
               do j=1,sg%lon2
                  do i=1,sg%lat2
                     xx%r3(i_cw)%q(i,j,k)= sub_ci(i,j,k)+sub_cl(i,j,k)
                  enddo
               enddo
            enddo
         endif
      endif

      if (allocated(sub_ci)) then
         call gsi_bundlegetpointer(xx,'qi', i_ci, istatus)
         if (i_ci>0) then ! present in xx
            do k=1,nsig
               do j=1,sg%lon2
                  do i=1,sg%lat2
                     xx%r3(i_ci)%q(i,j,k)= sub_ci(i,j,k)
                  enddo
               enddo
            enddo
         endif
         deallocate(sub_ci)
      endif

      if (allocated(sub_cl)) then
         call gsi_bundlegetpointer(xx,'ql', i_cl, istatus)
         if (i_cl>0) then ! present in xx
            do k=1,nsig
               do j=1,sg%lon2
                  do i=1,sg%lat2
                     xx%r3(i_cl)%q(i,j,k)= sub_cl(i,j,k)
                  enddo
               enddo
            enddo
         endif
         deallocate(sub_cl)
      endif

      call gsi_bundlegetpointer(xx,'qr', i_cr, istatus)
      if (i_cr>0) then ! present in xx
         j_cr = MAPL_SimpleBundleGetIndex ( xpert, 'qrtot', 3, rc=status )
         if(j_cr>0) then
           allocate( sub_cr(sg%lat2,sg%lon2,nsig),stat=ierr )
           if ( ierr/=0 ) then
               stat = 91
               if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_cr)'
               return
           end if
           call pert2gsi_ ( xpert%r3(j_cr)%qr4, sub_cr , ierr )
         endif
         if (i_cr>0) then ! present in xpert
            do k=1,nsig
               do j=1,sg%lon2
                  do i=1,sg%lat2
                     xx%r3(i_cr)%q(i,j,k)= sub_cr(i,j,k)
                  enddo
               enddo
            enddo
            deallocate(sub_cr)
         endif
      endif

      call gsi_bundlegetpointer(xx,'qs', i_cs, istatus)
      if (i_cs>0) then ! present in xx
         j_cs = MAPL_SimpleBundleGetIndex ( xpert, 'qstot', 3, rc=status )
         if(j_cs>0) then
           allocate( sub_cs(sg%lat2,sg%lon2,nsig),stat=ierr)
           if ( ierr/=0 ) then
               stat = 91
               if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_cs)'
               return
           end if
           call pert2gsi_ ( xpert%r3(j_cs)%qr4, sub_cs , ierr )
         endif
         if (i_cs>0) then ! present in xpert
            do k=1,nsig
               do j=1,sg%lon2
                  do i=1,sg%lat2
                     xx%r3(i_cs)%q(i,j,k)= sub_cs(i,j,k)
                  enddo
               enddo
            enddo
            deallocate(sub_cs)
         endif
      endif

      if (allocated(sub_ps)) then
         call gsi_bundlegetpointer(xx,'ps'  , i_ps, istatus);ierr=ierr+istatus
         if (i_ps>0) then
            do j=1,sg%lon2
               do i=1,sg%lat2
                  xx%r2(i_ps)%q(i,j)= sub_ps(i,j) * kPa_per_Pa
               enddo
            enddo
         endif
         deallocate(sub_ps)
      endif

      if (allocated(sub_ph)) then
         call gsi_bundlegetpointer(xx,'z'  ,  i_ph, ierr)
         if (i_ph>0) then
            do j=1,sg%lon2
               do i=1,sg%lat2
                  xx%r2(i_ph)%q(i,j)= sub_ph(i,j) / grav
               enddo
            enddo
         endif
         deallocate(sub_ph)
      endif

      if (allocated(sub_ts)) then
         call gsi_bundlegetpointer(xx,'sst',  i_ts, ierr)
         if (i_ts>0) then
            do j=1,sg%lon2
               do i=1,sg%lat2
                  xx%r2(i_ts)%q(i,j)= sub_ts(i,j)
               enddo
            enddo
         endif
         deallocate(sub_ts)
      endif

!     check on essential vars
      if ( ierr/=0 ) then
          stat = 99
          if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_)'
          return
      end if
      if(allocated(sub_cl)) deallocate(sub_cl)
      if(allocated(sub_ci)) deallocate(sub_ci)

      CONTAINS

      subroutine pert2gsi2d_ ( fld, sub, stat_ )

      real(4),        intent(in)  :: fld(:,:)
      real(r_kind),   intent(out) :: sub(:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*pert2gsi2d_'

      integer(i_kind) i,j

      stat_= 0

      call timer_ini('pert2gsi_')

      do j=1,sg%lon2
         do i=1,sg%lat2
            sub(i,j) = fld(j,i)
         end do
      end do

      call timer_fnl('pert2gsi_')
      end subroutine pert2gsi2d_

      subroutine pert2gsi_ ( fld, sub, stat_ )

      real(4),        intent(in)  :: fld(:,:,:)
      real(r_kind),   intent(out) :: sub(:,:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*pert2gsi_'

      integer(i_kind) i,j,k,mm1,ierr
      integer(i_kind) imr,jnp,nl

      mm1 = mype+1 
      stat_= 0
      imr  =sg%nlon
      jnp  =sg%nlat
      nl   =nsig
      if(size(fld,3)==1) nl=1

      call timer_ini('pert2gsi_')

      do k=1,nl
         do j=1,sg%lon2
            do i=1,sg%lat2
               sub(i,j,k) = fld(j,i,k)
            end do
         end do
      enddo

      call timer_fnl('pert2gsi_')
      end subroutine pert2gsi_

   end subroutine gcm2gsi1_

!------ BELOW THIS POINT: Routines of general (internal-only) use --------

!-------------------------------------------------------------------------
   subroutine SwapV_(fld)
!-------------------------------------------------------------------------
   implicit none
   real(r_kind),intent(inout) ::  fld(:,:,:)
   real(r_kind),allocatable   :: work(:,:,:)
   integer(i_kind) im, jm, km
   im   = size(fld,1)
   jm   = size(fld,2)
   km   = size(fld,3)
   allocate (work(im,jm,km))
   work = fld
   fld(:,:,km:1:-1) = work(:,:,1:km:+1)
   deallocate (work)
   end subroutine SwapV_
!-------------------------------------------------------------------------
   subroutine SwapIJK_(aij,aji)
!-------------------------------------------------------------------------
! transpose IJK-ordered array to JIK-ordered array
   implicit none
   real(r_kind),dimension(:,:,:),intent(in)    :: aij
   real(r_kind),dimension(:,:,:),intent(inout) :: aji
   integer(i_kind) :: i,k,isz,jsz,ksz,kk
!
   isz=size(aij,1)
   jsz=size(aij,2)
   ksz=size(aij,3)
   kk=1
   do k=1,ksz
      do i=1,isz
         aji(1:jsz,i,kk)=aij(i,1:jsz,k)
      end do
      kk=kk+1
   end do
   call SwapV_(aji)
   end subroutine SwapIJK_

   subroutine hflip3_ ( q,im,jm,km )
      implicit none
      integer(i_kind)  im,jm,km,i,j,k
      real(r_kind), intent(inout) :: q(im,jm,km)
      real(r_kind), allocatable   :: dum(:)
      allocate ( dum(im) )
      do k=1,km
      do j=1,jm
      do i=1,im/2
         dum(i) = q(i+im/2,j,k)
         dum(i+im/2) = q(i,j,k)
      enddo
         q(:,j,k) = dum(:)
      enddo
      enddo
      deallocate ( dum )
   end subroutine hflip3_
  
   subroutine hflip2_ ( q,im,jm )
      implicit none
      integer(i_kind)  im,jm,i,j
      real(r_kind), intent(inout) :: q(im,jm)
      real(r_kind), allocatable   :: dum(:)
      allocate ( dum(im) )
      do j=1,jm
      do i=1,im/2
         dum(i) = q(i+im/2,j)
         dum(i+im/2) = q(i,j)
      enddo
         q(:,j) = dum(:)
      enddo
      deallocate ( dum )
   end subroutine hflip2_

   subroutine myscatter_ ( what, GSI_xGrid, im_world, jm_world, km_world, vari, varo )
!  NOTE: there is not need for Halo here ... this is going onto a "model" field
   use ESMF
   implicit none
   character(len=*),intent(in)    :: what
   type(ESMF_Grid), intent(in)    :: GSI_xGrid
   integer,         intent(in)    :: im_world,jm_world,km_world
   real(r_kind),    intent(in)    :: vari(im_world,jm_world,km_world)
   real(4),         intent(inout) :: varo(:,:,:)

   character(len=*),parameter :: Iam=myname//'myscatter_'
   real(4),allocatable :: work2d(:,:)
   integer  k, rc,status

    allocate(work2d(im_world,jm_world))
    if (what=='adm') then
        do k=1,km_world  ! this would be a good place to swapV
           work2d=vari(:,:,k)
           call ArrayScatter(varo(:,:,k), work2d, GSI_xGrid, rc=status)   ! _RT: ideally this should be the AD of gather
                VERIFY_(STATUS)
        enddo
    else
        do k=1,km_world  ! this would be a good place to swapV
           work2d=vari(:,:,k)
           call ArrayScatter(varo(:,:,k), work2d, GSI_xGrid, rc=status)
                VERIFY_(STATUS)
        enddo
    endif
    deallocate(work2d)
    
   end subroutine myscatter_

   subroutine mygather_ ( what, GSI_xgrid, im_world, jm_world, km_world, vari, varo )
!  NOTE: there is not need for Halo here ... this is coming from a "model" field
   use ESMF
   implicit none
   character(len=*), intent(in)    :: what
   type(ESMF_Grid),  intent(in)    :: GSI_xGrid
   integer,          intent(in)    :: im_world,jm_world,km_world
   real(4),          intent(in)    :: vari(:,:,:)
   real(r_kind),     intent(inout) :: varo(im_world,jm_world,km_world)

   character(len=*),parameter :: Iam=myname//'mygather_'
   real(4),allocatable :: work2d(:,:)
   integer  k, rc,status

    allocate(work2d(im_world,jm_world))
    if (what=='adm') then
        do k=1,km_world  ! this would be a good place to swapV
           work2d=zero
           call ArrayGather(vari(:,:,k), work2d, GSI_xGrid, rc=status)   ! _RT: ideally this should be the AD of gather
                VERIFY_(STATUS)
           varo(:,:,k) = work2d
        enddo
    else
        do k=1,km_world  ! this would be a good place to swapV
           work2d=zero
           call ArrayGather(vari(:,:,k), work2d, GSI_xGrid, rc=status)
                VERIFY_(STATUS)
           varo(:,:,k) = work2d
        enddo
    endif
    deallocate(work2d)
    
   end subroutine mygather_

   subroutine reorder21_(sg,grid_in,grid_out)
   use kinds, only: r_kind,r_double
   use general_sub2grid_mod, only: sub2grid_info
   implicit none

   type(sub2grid_info),                    intent(in  ) :: sg        ! subdomain info
   real(r_kind),dimension(sg%nlon,sg%nlat),intent(in  ) :: grid_in   ! input grid
   real(r_kind),dimension(sg%itotsub)     ,intent( out) :: grid_out  ! output gridend subroutine reorder21_

!  Declare local variables
   integer(i_kind) i,j,k

!  Transfer input 2d array to output 1d array
   do k=1,sg%itotsub
      i=sg%ltosi_s(k)
      j=sg%ltosj_s(k)
      grid_out(k)=grid_in(j,i)
   end do

   end subroutine reorder21_

   subroutine reorder12_(sg,grid_in,grid_out)
   use kinds, only: r_kind,r_double
   use general_sub2grid_mod, only: sub2grid_info
   implicit none

   type(sub2grid_info),                               intent(in   ) :: sg        ! subdomain info
   real(r_kind),dimension(max(sg%iglobal,sg%itotsub)),intent(in   ) :: grid_in   !  input grid
   real(r_kind),dimension(sg%nlon,sg%nlat)           ,intent(  out) :: grid_out  !  input grid

!  Declare local variables
   integer(i_kind) i,j,k

!  Transfer input 1d array to output 2d array
   do k=1,sg%iglobal
      i=sg%ltosi(k)
      j=sg%ltosj(k)
      grid_out(j,i) = grid_in(k)
   end do

   end subroutine reorder12_

   subroutine strip_(sg,field_in,field_out,nz)

! !USES:

    use kinds, only: r_kind
   use general_sub2grid_mod, only: sub2grid_info
    implicit none

! !INPUT PARAMETERS:

    type(sub2grid_info)                       , intent(in)    :: sg         !  subdomain info
    integer(i_kind)                           , intent(in   ) :: nz         !  number of levs in subdomain array
    real(r_kind),dimension(sg%lat2,sg%lon2,nz), intent(in   ) :: field_in   ! full subdomain 
                                                                            !    array containing 
                                                                            !    buffer points
! !OUTPUT PARAMETERS:

    real(r_kind),dimension(sg%lat1,sg%lon1,nz), intent(  out) :: field_out ! subdomain array
                                                                           !   with buffer points
                                                                           !   stripped off

! !DESCRIPTION: strip off buffer points froms subdomains for mpi comm
!               purposes
!
! !REVISION HISTORY:
!
!   2004-01-25  kleist
!   2004-05-14  kleist, documentation
!   2004-07-15  todling, protex-compliant prologue
!   2011-11-11  todling, adapt to work generally with subgrid info
!
! !REMARKS:
!
!   language: f90
!   machine:  ibm rs/6000 sp; sgi origin 2000; compaq/hp
!
! !AUTHOR: 
!    kleist           org: np20                date: 2004-01-25
!
!EOP
!-------------------------------------------------------------------------

    integer(i_kind) i,j,k,jp1

    do k=1,nz
       do j=1,sg%lon1
          jp1 = j+1
          do i=1,sg%lat1
             field_out(i,j,k)=field_in(i+1,jp1,k)
          end do
       end do
    end do

    return
   end subroutine strip_
end module geos_StateIO
