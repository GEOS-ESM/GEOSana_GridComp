!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_IOState --- GFIO-type Input/Output of GCM state vector
!
! !INTERFACE:
!
      module m_IOState

! !USES:
#if defined( SPMD )
      use mod_comm, only : mp_gather4d,mp_scatter2d,mp_scatter4d,mp_bcst_int
#endif

      use precision

      use control, only : ibeg_u   , iend_u
      use control, only : ibeg_v   , iend_v
      use control, only : ibeg_delp, iend_delp
      use control, only : ibeg_pt  , iend_pt
      use control, only : ibeg_q   , iend_q
      use control, only : ibeg_ps  , iend_ps
      use control, only : nvct
      use control, only : mod2vect

      use stepon, only : ng_d, ng_s
      use stepon, only : nsecf
      use stepon, only : job, trajfrq, rmaptrj, trajtmpl
      use stepon, only : nouter, ninner, pdt, ptop, ks

      use m_FileResolv    		! file resolver
      use m_strtemplate

      use m_mpif90, only : MP_comm_world
      use m_mpif90, only : MP_comm_rank

      use mod_comm, only : imr      ! no. of grid points in the zonal direction
      use mod_comm, only : jnp      ! no. of grid points in the meridional direction
      use mod_comm, only : nl       ! no. of levels
      use mod_comm, only : nc       ! no. of tracers

      use prognostics, only : jfirst   ! initial lat index for this pe
      use prognostics, only : jlast    ! last    lat index for this pe
      use prognostics, only : kfirst   ! 
      use prognostics, only : klast    ! 

      use prognostics, only : dyn_prog

      use prognostics, only : sgh      ! Topography height stdv (meter)
      use prognostics, only : phis     ! topography geopotential (g*z_s)

      use prognostics, only : ps       ! surface pressure
      use prognostics, only : pe
      use prognostics, only : pk
      use prognostics, only : pkz
      use prognostics, only : peln

      use timingModule

      use m_dyn

      use m_inpak90

      use m_die, only : MP_die
      use m_die, only : die
      use m_die, only : perr

      implicit none

! !PUBLIC MEMBER FUNCTIONS:

      PUBLIC PutState
      PUBLIC GetState

      PUBLIC GetState_init
      PUBLIC GetState_clean
 
      interface PutState; module procedure
     .          PutState1_,
     .          PutState2_
      end interface

      interface GetState; module procedure 
     .        getstate_,
     .        retrieve_,
     .        all_
      end interface

      interface GetState_init
        module procedure init_
      end interface

      interface GetState_clean
        module procedure clean_
      end interface

!
! !DESCRIPTION: Controls I/O of GFIO-type of GCM state.
!
! !REVISION HISTORY:
!
!  19Nov2002  Todling   Initial design/interfaces.
!  25Feb2004  Todling   Added interface to allow param specs via rc file
!  22Nov2006  Oloso     Changes for SPMD compatibility
!  20Apr2007  Kim       Changes for GetSatet_all_ to read in all trajectories
!  03May2007  Todling   Option to read trajectory via i/o
!  07May2007  Todling   Cleared some redundancies
!  12May2007  Todling   All use to the top; introduced dyn_prog
!  14May2007  Todling   Introduced dyn_prog; global change.
!  28Oct2009  Todling   All dyn-get calls deal w/ geos-4 type vectors here
!
!EOP
!-------------------------------------------------------------------------

      PRIVATE
      save

      real(r8), allocatable :: u_all(:,:,:,:)   ! zonal wind on D-grid
      real(r8), allocatable :: v_all(:,:,:,:)   ! meridional wind
      real(r8), allocatable :: pt_all(:,:,:,:)  ! virtual potential temperature
      real(r8), allocatable :: ps_all(:,:,:)    ! surface pressure (pascal)
      real(r8), allocatable :: delp_all(:,:,:,:)! pressure thickness (pascal)
      real(r8), allocatable :: q_all(:,:,:,:,:) ! specific humidity & tracer mixing ratios

      real(r8), pointer ::    u(:,:,:)   ! zonal wind on D-grid
      real(r8), pointer ::    v(:,:,:)   ! meridional wind
      real(r8), pointer ::   pt(:,:,:)   ! virtual potential temperature
      real(r8), pointer :: delp(:,:,:)   ! pressure thickness (pascal)
      real(r8), pointer ::    q(:,:,:,:) ! specific humidity & tracer mixing ratios

      integer, allocatable :: outer_loop_nymd(:)
      integer, allocatable :: outer_loop_nhms(:)

!     real(r8), pointer ::   ps(:,:)     ! surface pressure (pascal)

      integer, parameter :: ROOT = 0     ! should come from above

      character(len=*), parameter :: myname = 'm_IOState'
      character(len=*), parameter :: trajfn_def = 'trajectory'

      logical, save :: MEMTRJ_  = .true.
      logical, save :: verbose_ = .false.

      integer :: w_f_grid_lm
      integer :: nouter_trj
      integer :: nymd_t0, nhms_t0


      CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PutState1_ --- 
!
! !INTERFACE:

      subroutine PutState1_ ( outfile, 
     &                        job, nymd, nhms, prog, freq, nstep,
     &                        ak, bk, Ts, oro, prec, stat )

! !USES:


      implicit none

! !INPUT PARAMETERS:

      character(len=*)     outfile

      character(len=*), intent(in) :: job
      integer,          intent(in) :: nymd   ! Date: year-month-day
      integer,          intent(in) :: nhms   ! Time: hour-min-sec
      type(dyn_prog),   intent(in) :: prog   ! prognostic variables
      integer,          intent(in) :: freq   ! Analysis frequency (HHMMSS)
      integer,          intent(in) :: nstep  ! model time step

      real,    intent(in)  ::  ak(nl+1)      ! vertical grid a-coefficient
      real,    intent(in)  ::  bk(nl+1)      ! vertical grid b-coefficient
      real,    intent(in)  ::  Ts(imr,jnp)   ! Sea surface temperature (K)
      real,    intent(in)  ::  oro(imr,jnp)  ! Land-water-ice mask:
                                             !   oro = 0  over ocean
                                             !   oro = 1  over land
                                             !   oro = 2  over sea ice

      integer, intent(in), optional ::  prec ! 0 = 32 bits
                                             ! 1 = 64 bits

! !OUTPUT PARAMETERS:

      integer, intent(out), optional :: stat

! !DESCRIPTION: Writes out state vector as a GFIO file.
!
! !TO DO: 
!
!  (1) Generalize to output at any specified frequency other
!      than model time step (pdt)
!  (2) MPI version
!
! !REVISION HISTORY:
!
!  19Nov2002  Todling    Initial code.
!  11Mar2003  Todling    Removed circular dependence on stepon.
!  26Jan2007  Todling    SPMD implementation
!  06Mar2009  Todling    Update dyn_init interface
!  19Jan2010  Todling    Gather phis and sgh locally
!
!EOP
!-------------------------------------------------------------------------


      character(len=*), parameter :: myname_ = myname//'*PutState1_'

      type (dyn_vect) :: w_f

      integer  prec_out
      integer  myid, ierr
      integer, parameter :: prec_def = 1  ! no. of tracers (where in GCM?)
      real(r8) ::  myphis(imr,jnp),mysgh(imr,jnp)
      real(r8), allocatable ::  xvect(:)  ! u,v,pt,q,ps
      real(r8), allocatable ::  aux2d(:,:)! temp 2d array for land/etc fractions

      if(present(stat)) stat = 0

      call MP_comm_rank(MP_comm_world,myID,ierr)
         if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)

      if ( freq == 0 ) then
           print *, ' Invalid frequency to do I/O ... skip it '
           return
      endif

!     Allocate local arrays
!     ---------------------
      allocate( xvect(nvct), stat=ierr )

      call mod2vect( nvct, xvect, prog )

      call mp_gather4d(phis, myphis, imr, jnp, 1, 1, jfirst, jlast, 1, 1, 0, 0, 0)
      call mp_gather4d(sgh , mysgh , imr, jnp, 1, 1, jfirst, jlast, 1, 1, 0, 0, 0)

      if ( myid==ROOT ) then

!     Allocate aux array
!     ------------------
      allocate(aux2d(imr,jnp))
      aux2d=0.d0  ! can't do much better for land/etc fractions for now

!     Make sure w_f does not point to anything
!     ----------------------------------------
      call dyn_null ( w_f )

!     Assign pointers in w_f input arrays
!     -----------------------------------
      call Dyn_Init ( imr, jnp, nl, nc, ptop, ks, ak, bk,
     &                myphis, mysgh, Ts, oro,
     &                aux2d, aux2d, aux2d, aux2d, aux2d,
     &                reshape(xvect(ibeg_ps  :iend_ps)  ,(/imr,jnp/)),
     &                reshape(xvect(ibeg_delp:iend_delp),(/imr,jnp,nl/)),
     &                reshape(xvect(ibeg_u   :iend_u   ),(/imr,jnp,nl/)), 
     &                reshape(xvect(ibeg_v   :iend_v   ),(/imr,jnp,nl/)),
     &                reshape(xvect(ibeg_pt  :iend_pt  ),(/imr,jnp,nl/)), 
     &                reshape(xvect(ibeg_q   :iend_q   ),(/imr,jnp,nl,nc/)),
     &                w_f, ierr  )
      if ( ierr/=0 ) then
           if (present(stat)) then
               stat = 2
               return
           else
               call MP_die(myname_,'error from Dyn_Init()', ierr) 
           endif
      end if

      prec_out = prec_def
      if(present(prec)) prec_out = prec

!     Write out
!     ---------
      call dyn_put ( outfile, nymd, nhms, prec_out, w_f, ierr, freq=freq )
	if(ierr==0) then
           if (myid==ROOT) then
               print *, myname_//': wrote dynamics state file '//trim(outfile)
               print *, myname_//': nymd, nhms = ', nymd, nhms
               print *
           endif
        else
             call MP_die (myname_,'Error from dyn_put',ierr)
        endif

!     Release pointers
!     ----------------
      call dyn_null ( w_f )
      deallocate(aux2d)

      end if ! < myid >

      deallocate( xvect )


      end subroutine PutState1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetState_:  Read in model state
!
! !INTERFACE:

      subroutine GetState_ ( infile, job, nymd, nhms, prog, stat )

! !USES:

      implicit none

      include "comcon.h"               ! only : cappa

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: infile
      character(len=*), intent(in) :: job

      integer,          intent(in) :: nymd   ! Date: year-month-day
      integer,          intent(in) :: nhms   ! Time: hour-min-sec

      type(dyn_prog),   intent(inout) :: prog ! prognostic variables

! !OUTPUT PARAMETERS:

      integer, intent(out), optional :: stat

! !DESCRIPTION: Reads in model state vector.
!
! !TO DO: 
!
!  (1) Generalize to allow the option of reading the nearest-time 
!      available in file to the input (nymd,nhms) time.
!
! !REVISION HISTORY:
!
!  19Nov2002  Todling   Initial code (as per ana.f).
!  25Sep2003  RG/RT     Bug fix: added dyn_null.
!  29Jan2004  Todling   Fixed exit for "present(stat)" and field not found
!  30Mar2005  Gelaro/RT Fixed handle for tracers>1
!  28Mar2005  Todling   Taking orography from input traj; zero out d-grid u-wind
!                       large (fill) values just in case.
!  22Nov2006  Oloso     Changes for SPMD compatibility
!  30Jun2007  Todling   Read on root and scatter
!  19Jan2010  Todling   Bug fix: q_all%q was being used instead of prog%q
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*GetState_'


      real(r8) ptopk
      integer  nymd_in, nhms_in
      integer  myid, i, j, k, n, w_f_grid_lm
      integer  ier
      character(len=255) :: fname

      real(r8), allocatable ::   pstmp(:,:),    utmp(:,:,:), vtmp(:,:,:)
      real(r8), allocatable :: delptmp(:,:,:), pttmp(:,:,:), qtmp(:,:,:,:)

      type (dyn_vect) :: w_f

      if(present(stat)) stat=0

      call MP_comm_rank(MP_comm_world,myID,ier)
         if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

      fname = trim(infile)

!     Get memory for temporary arrays
!     -------------------------------
      if ( myid==ROOT ) then
          allocate (   pstmp(imr,jnp),    utmp(imr,jnp,nl),  vtmp(imr,jnp,nl),
     .               delptmp(imr,jnp,nl), pttmp(imr,jnp,nl), qtmp(imr,jnp,nl,nc) )
      else
          allocate (   pstmp(0,0),    utmp(0,0,0), vtmp(0,0,0),
     .               delptmp(0,0,0), pttmp(0,0,0), qtmp(0,0,0,0) )
      endif

!     Set date/time of trajectory
!     ---------------------------
      nymd_in = nymd
      nhms_in = nhms

      if ( myid== ROOT ) then

!          Make sure w_f does not point to anything
!          ----------------------------------------
           call dyn_null ( w_f )

!          Load dynamics vector
!          -------------------
    	   call timing_on ( 'dyn_get' )

           call dyn_get ( trim(infile), nymd_in, nhms_in, w_f, ier, timidx=0, vectype=4 )
	   if ( ier .eq. 0 ) then
	        print *, myname_//': read dynamics state file '//trim(fname)
	        print *, myname_//': nymd, nhms = ', nymd, nhms
     	        print *
           else
              if (present(stat)) then
                 stat = 1  ! cannot read state
	         print *, myname_//': read dynamics failed to read state'
                 call dyn_clean ( w_f )
                 return
              else
                call die(myname_,'could not read dynamics state file '//trim(fname))
              end if
           end if
  	   call timing_off ( 'dyn_get' )

             pstmp = w_f%ps
              utmp = w_f%u
              vtmp = w_f%v
             pttmp = w_f%pt
           delptmp = w_f%delp
              qtmp = w_f%q

           w_f_grid_lm = w_f%grid%lm

!          Release root memory
!          -------------------
           call dyn_clean ( w_f )

      endif ! < myid=ROOT >

      call mp_bcst_int (w_f_grid_lm)

      call mp_scatter2d(  pstmp,        ps, imr,  jnp, jfirst, jlast,                            0)
      call mp_scatter4d(   utmp,    prog%u, imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_s, 0)
      call mp_scatter4d(   vtmp,    prog%v, imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_s, ng_d, 0)
      call mp_scatter4d(  pttmp,   prog%pt, imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
      call mp_scatter4d(delptmp, prog%delp, imr,  jnp, nl,  1, jfirst, jlast, 1, nl,    0,    0, 0)
      call mp_scatter4d(   qtmp,    prog%q, imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)

      if ( nc .gt. 1 ) then
         if ( nc .le. w_f_grid_lm ) then
            do n = 2, nc
               call mp_scatter4d (qtmp(:,:,:,n),  prog%q(:,:,:,n), imr,  jnp, nl, 1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
             enddo
         else
            call die(myname_,'not enough tracers ')
         end if
      end if

!     Release temporary memory
!     ------------------------
      deallocate ( pstmp, utmp, vtmp, delptmp, pttmp, qtmp )

!     Recalculate phz, pe, pk, and pkz
!     NOTE1: line below extracted from fvgcm.F
!     NOTE2: openMP intentionally deactivated
!     ---------------------------------------
      ptopk = ptop ** cappa

!___$omp parallel do private(i,j,k)
      do j=jfirst,jlast
        do i=1,imr
           pe(i,1,j) = ptop
        enddo

        do k=2,nl+1
        do i=1,imr
           pe(i,k,j) = pe(i,k-1,j) + prog%delp(i,j,k-1)
        enddo
        enddo
      enddo

!___$omp parallel do default(shared) private(i,j,k)
        do j=jfirst,jlast
           do i=1,imr
              pk(i,j,1) = ptopk
           enddo

           do k=2,nl+1
              do i=1,imr
                 pk(i,j,k) = pe(i,k,j)**cappa
              enddo
           enddo
        enddo

        call pkez(1, imr, nl, jfirst, jlast, 
     &            pe, pk, cappa, ks, peln, pkz, .false.)

      end subroutine GetState_


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetState77_:  Read in model state
!
! !INTERFACE:

      subroutine GetState77_ ( infile,
     &                         nymd, nhms, prog,
     &                         imr, jnp, km, nq, jfirst, jlast, ng_s, ng_d,
     &                         ps, delp, u, v, pt, q,
     &                         stat )

! !USES:

      implicit none

      include "comcon.h"               ! only : cappa

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: infile

      integer :: nymd     ! date
      integer :: nhms     ! time
      type(dyn_prog), intent(inout) :: prog
      integer :: imr      ! no. of grid points in the zonal direction
      integer :: jnp      ! no. of grid points in the meridional direction
      integer :: km       ! no. of levels
      integer :: nq       ! no. of tracers
      integer :: jfirst   ! 
      integer :: jlast    ! 
      integer :: ng_s     ! ghosting parameter 
      integer :: ng_d     ! ghosting parameter

! !OUTPUT PARAMETERS:

      real(r8) :: ps(imr,jfirst:jlast) 
      real(r8) :: u(imr,jfirst-ng_d:jlast+ng_s,km) 
      real(r8) :: v(imr,jfirst-ng_s:jlast+ng_d,km) 
      real(r8) :: pt(imr,jfirst-ng_d:jlast+ng_d,km) 
      real(r8) :: delp(imr,jfirst:jlast,km)
      real(r8) :: q(imr,jfirst-ng_d:jlast+ng_d,km,nq)

      integer, intent(out), optional :: stat

! !DESCRIPTION: Reads in model state vector.
!
! !TO DO: 
!
!  (1) Generalize to allow the option of reading the nearest-time 
!      available in file to the input (nymd,nhms) time.
!  (2) MPI conterpart
!  (3) Do remap
!
! !REVISION HISTORY:
!
!  11Mar2003  Todling   Initial code NOT COMPLETED: needs work
!  22Nov2006  Oloso     Changes for SPMD compatibility
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*GetState77_'

      real(r8) ptopk
      integer  nymd_in, nhms_in
      integer  myid, i, j, k
      integer  ier
      character(len=255) :: fname

      type (dyn_vect) :: w_f

      if(present(stat)) stat=0

      call MP_comm_rank(MP_comm_world,myID,ier)
         if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

      fname = trim(infile)

!     Set date/time of trajectory
!     ---------------------------
      nymd_in = nymd
      nhms_in = nhms

!     Load dynamics vector
!     -------------------
	call timing_on ( 'dyn_get' )
      call dyn_get ( trim(infile), nymd_in, nhms_in, w_f, ier, timidx=0, vectype=4 )
	if ( ier .eq. 0 ) then
           if ( myid==ROOT ) then
	       print *, myname_//': read dynamics state file '//trim(fname)
	       print *, myname_//': nymd, nhms = ', nymd, nhms
  	       print *
           endif
        else
          if (present(stat)) then
             stat = 1  ! cannot read state
	     print *, myname_//': read dynamics failed to read state'
          else
            call die(myname_,'could not read dynamics state file '//trim(fname))
          end if
        end if
	call timing_off ( 'dyn_get' )

!     Copy state to model arrays and release w_f allocation
!     -----------------------------------------------------
      prog%delp(:,jfirst:jlast,:)   = w_f%delp(:,jfirst:jlast,:)
      prog%u   (:,jfirst:jlast,:)   = w_f%u   (:,jfirst:jlast,:)
      prog%v   (:,jfirst:jlast,:)   = w_f%v   (:,jfirst:jlast,:)
      prog%pt  (:,jfirst:jlast,:)   = w_f%pt  (:,jfirst:jlast,:)
      prog%q   (:,jfirst:jlast,:,1) = w_f%q   (:,jfirst:jlast,:,1)

!     prog%ps  (:,jfirst:jlast)     = w_f%ps  (:,jfirst:jlast)
           ps  (:,jfirst:jlast)     = w_f%ps  (:,jfirst:jlast)


      call dyn_clean ( w_f )

!     Recalculate phz, pe, pk, and pkz
!     NOTE1: line below extracted from fvgcm.F
!     NOTE2: openMP intentionally deactivated
!     ---------------------------------------
      ptopk = ptop ** cappa

!___$omp parallel do private(i,j,k)
      do j=jfirst,jlast
        do i=1,imr
           pe(i,1,j) = ptop
        enddo

        do k=2,km+1
        do i=1,imr
           pe(i,k,j) = pe(i,k-1,j) + prog%delp(i,j,k-1)
        enddo
        enddo
      enddo

!___$omp parallel do default(shared) private(i,j,k)
        do j=jfirst,jlast
           do i=1,imr
              pk(i,j,1) = ptopk
           enddo

           do k=2,km+1
              do i=1,imr
                 pk(i,j,k) = pe(i,k,j)**cappa
              enddo
           enddo
        enddo

        call pkez(1, imr, km, jfirst, jlast, 
     &            pe, pk, cappa, ks, peln, pkz, .false.)

      end subroutine GetState77_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PutState2_ --- 
!
! !INTERFACE:

      subroutine PutState2_ ( job, nymd, nhms, prog, nstep,
     &                        ak, bk, Ts, oro,
     &                        rcfname, prec, stat )


! !USES:


      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: job
      integer,          intent(in) :: nymd   ! Date: year-month-day
      integer,          intent(in) :: nhms   ! Time: hour-min-sec
      type(dyn_prog),   intent(inout) :: prog

      integer,          intent(in) :: nstep  ! model time step

      real,    intent(in)  ::  ak(nl+1)      ! vertical grid a-coefficient
      real,    intent(in)  ::  bk(nl+1)      ! vertical grid b-coefficient
      real,    intent(in)  ::  Ts(imr,jnp)   ! Sea surface temperature (K)
      real,    intent(in)  ::  oro(imr,jnp)  ! Land-water-ice mask:
                                             !   oro = 0  over ocean
                                             !   oro = 1  over land
                                             !   oro = 2  over sea ice
 
      character(len=*), intent(in) :: rcfname

      integer, intent(in), optional ::  prec ! 0 = 32 bits
                                             ! 1 = 64 bits

! !OUTPUT PARAMETERS:

      integer, intent(out), optional :: stat

! !DESCRIPTION: Writes out state vector as a GFIO file.
!
! !REVISION HISTORY:
!
!  25Feb2004  Todling    Added this interface
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*PutState2_'

      character(len=255) :: fname          ! trajectory filename
      integer            :: freq           ! trajectory output freq

      integer strajfrq, snhms, ierror

      if(present(stat)) stat = 0

!     Get info on filename and frequency from resouce file
!     ----------------------------------------------------
      call putinit_(ierror)
         if(ierror/=0) return

      if ( freq .gt. 0 ) then
        strajfrq = nsecf(freq)
        snhms    = nsecf(nhms)
        if ( mod(snhms,strajfrq) == 0 ) then
          call putstate1_ ( fname,
     &                      job, nymd, nhms, prog, freq, nstep,
     &                      ak, bk, Ts, oro,
     &                      stat=ierror )
        endif
      endif

      CONTAINS

      subroutine putinit_ (stat)

      integer stat

      character(len=*), parameter :: DefRCfn = 'iostate.rc' ! default RC file name
      integer, parameter :: freq_def = 003000  ! default i/o freq
      character(len=255) :: fname_def          ! default filename

      character(len=255) template
      character(len=8)   cymd

      integer envin, ier
      logical exist

      stat = 0
      inquire(file=trim(rcfname),exist=exist)
      if (.not.exist) then
         stat = 99
         return
      endif

!     Default analysis-related filenames
!     ----------------------------------
      write(cymd,'(i8.8)') nymd
      fname_def = trim(job) // '.nlms.eta.' // cymd // '.hdf' 

      fname  = trim(fname_def)
      freq   = freq_def

!     Load resources
!     -------------
      call i90_loadf ( trim(rcfname), ier )
        if ( ier/=0 ) return   ! if rc not found, return silently

        call i90_label ( 'nlm_statevector_filename:', ier )
        if ( ier .ne. 0 ) then
           print *, 'taking nlm state-vec  default filename '
        else
          call i90_gtoken ( template, ier )
             if(ier.ne.0) call die (myname, 'premature end of rc file')
             call strTemplate ( fname, template, 'GRADS',
     &                          xid=trim(job), nymd=nymd, nhms=nhms )
        endif

        call i90_label ( 'nlm_statevector_frequency:', ier )
        if ( ier .ne. 0 ) then
           print *, 'taking nlm state-vec i/o default default '
        else
          freq = i90_GInt( ier )
            if(ier.ne.0) call die (myname, 'premature end of rc file')
        endif

        call i90_release()

      end subroutine putinit_

      end subroutine putstate2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: all_:  Read in all trajectory model states
!
! !INTERFACE:

      subroutine all_ ( nymd_in, nhms_in )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: nymd_in   ! Date: year-month-day
      integer, intent(in) :: nhms_in   ! Time: hour-min-sec

! !DESCRIPTION: Reads in model state vector.
!
! !REVISION HISTORY:
!
!  20Apr2007  Kim     Initial code (as per GetState_).
!  05Jun2007  Todling Bug fix; was not handling case trajfrq<0 right
!  30Jun2007  Todling Read on root and scatter
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*all_'
      integer  myid, i, j, k, n, iouter,iinner, i_trj, i_test
      integer  ier
      integer  nymd_trj, nhms_trj
      integer  nymd_all, nhms_all
      integer  strajfrq, snhms
      integer  rmaptfrq
      character(len=255) :: fname
      character*255 :: fnametrj

      real(r8), allocatable ::   pstmp(:,:),    utmp(:,:,:), vtmp(:,:,:)
      real(r8), allocatable :: delptmp(:,:,:), pttmp(:,:,:), qtmp(:,:,:,:)

      type (dyn_vect) :: w_f

!     If trajectory not in memory ...
!     -------------------------------
      if ( .not. MEMTRJ_ ) return

      call MP_comm_rank(MP_comm_world,myID,ier)
         if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

      if (trajfrq<0) return

      nymd_all = nymd_in
      nhms_all = nhms_in

!     Read state trajectories at initial time step (iouter=0)
!     ------------------------------------------------------
      nymd_trj = nymd_all
      nhms_trj = nhms_all
      call strTemplate( fnametrj, trajtmpl,
     &                  xid=trim(job), nymd=nymd_trj, nhms=nhms_trj )
      fname = trim(fnametrj)

!     Get memory for temporary arrays
!     -------------------------------
      if ( myid==ROOT ) then
          allocate (   pstmp(imr,jnp),    utmp(imr,jnp,nl),  vtmp(imr,jnp,nl),
     .               delptmp(imr,jnp,nl), pttmp(imr,jnp,nl), qtmp(imr,jnp,nl,nc) )
      else
          allocate (   pstmp(0,0),    utmp(0,0,0), vtmp(0,0,0),
     .               delptmp(0,0,0), pttmp(0,0,0), qtmp(0,0,0,0) )
      endif

!     Read state on root
!     ------------------
      if ( myid==ROOT ) then

          call dyn_null ( w_f )
          call dyn_get ( trim(fnametrj), nymd_trj, nhms_trj, w_f, ier, timidx=0, vectype=4 )
          if ( ier .eq. 0 ) then
               print *, myname_//': read dynamics state file '//trim(fname)
               print *, myname_//': nymd, nhms = ', nymd_all, nhms_all
               print *
          endif

          delptmp = w_f%delp
            pstmp = w_f%ps
             utmp = w_f%u
             vtmp = w_f%v
            pttmp = w_f%pt
             qtmp(:,:,:,1) = w_f%q(:,:,:,1)

          w_f_grid_lm = w_f%grid%lm

          if ( nc .gt. 1 ) then
            if ( nc .le. w_f_grid_lm ) then
                 do n = 2, nc
                    qtmp(:,:,:,n) = w_f%q(:,:,:,n)
                 enddo
            else 
                 ier = 1
            endif
          endif

!         Release root memory
!         -------------------
          call dyn_clean ( w_f )

      endif ! < myid==ROOT >

      call mp_bcst_int (ier)
         if(ier/=0) call MP_die(myname_,'not enough tracers',ier)
      call mp_bcst_int (w_f_grid_lm)

      call mp_scatter2d(  pstmp,   ps_all(:,:,    0), imr,  jnp, jfirst, jlast,                            0)
      call mp_scatter4d(   utmp,    u_all(:,:,:,  0), imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_s, 0)
      call mp_scatter4d(   vtmp,    v_all(:,:,:,  0), imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_s, ng_d, 0)
      call mp_scatter4d(  pttmp,   pt_all(:,:,:,  0), imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
      call mp_scatter4d(delptmp, delp_all(:,:,:,  0), imr,  jnp, nl,  1, jfirst, jlast, 1, nl,    0,    0, 0)
      call mp_scatter4d(   qtmp(:,:,:,1), q_all(:,:,:,:,0), imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)

      if ( nc .gt. 1 ) then
         if ( nc .le. w_f_grid_lm ) then
            do n = 2, nc
               call mp_scatter4d (qtmp(:,:,:,n),  q_all(:,:,:,n,0), imr,  jnp, nl, 1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
             enddo
         else
            call die(myname_,'not enough tracers ')
         end if
      end if

      outer_loop_nymd(0) = nymd_trj
      outer_loop_nhms(0) = nhms_trj

!     Read the rest of the trajectory files
!     -------------------------------------
      do iouter = 1, nouter
        nymd_trj = nymd_all
        nhms_trj = nhms_all
        call tick( nymd_trj,nhms_trj,pdt )

        if (trajfrq>0) then
          strajfrq = nsecf(trajfrq)
          snhms    = nsecf(nhms_trj)
          if (mod(snhms,strajfrq)==0) then
            call strTemplate( fnametrj, trajtmpl,
     &                        xid=trim(job), nymd=nymd_trj, nhms=nhms_trj )
            fname = trim(fnametrj)

            if ( myid==ROOT ) then

                call dyn_null ( w_f )

!               Load single instance of dynamics vector
!               ----------------------------------------
                call dyn_get ( trim(fnametrj), nymd_trj, nhms_trj, w_f, ier, timidx=0, vectype=4 )
    
                if ( ier .eq. 0 ) then
                     print *, myname_//': read dynamics state file '//trim(fname)
                     print *, myname_//': nymd, nhms = ',nymd_trj,nhms_trj
                     print *
                endif

                delptmp = w_f%delp
                  pstmp = w_f%ps
                   utmp = w_f%u
                   vtmp = w_f%v
                  pttmp = w_f%pt
                   qtmp(:,:,:,1) = w_f%q(:,:,:,1)
 
               w_f_grid_lm = w_f%grid%lm

                if ( nc .gt. 1 ) then
                  if ( nc .le. w_f_grid_lm ) then
                       do n = 2, nc
                          qtmp(:,:,:,n) = w_f%q(:,:,:,n)
                       enddo
                  else 
                       ier = 1
                  endif
                endif

!              Release root memory
!              -------------------
               call dyn_clean ( w_f )

            endif  ! < myid==ROOT >

            call mp_bcst_int (w_f_grid_lm)

!           Copy state to model arrays and release w_f allocation
!           -----------------------------------------------------
            call Trj_idx_(nymd_trj,nhms_trj,.true.,i)

            call mp_scatter2d(  pstmp,   ps_all(:,:,    i), imr,  jnp, jfirst, jlast,                            0)
            call mp_scatter4d(   utmp,    u_all(:,:,:,  i), imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_s, 0)
            call mp_scatter4d(   vtmp,    v_all(:,:,:,  i), imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_s, ng_d, 0)
            call mp_scatter4d(  pttmp,   pt_all(:,:,:,  i), imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
            call mp_scatter4d(delptmp, delp_all(:,:,:,  i), imr,  jnp, nl,  1, jfirst, jlast, 1, nl,    0,    0, 0)
            call mp_scatter4d(   qtmp(:,:,:,1), q_all(:,:,:,:,i), imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)

            if ( nc .gt. 1 ) then
               if ( nc .le. w_f_grid_lm ) then
                    do n = 2, nc 
                       call mp_scatter4d (qtmp(:,:,:,n), q_all(:,:,:,n,i), imr,  jnp, nl, 1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
                    enddo
               else
                  call die(myname_,'not enough tracers ')
               end if
            end if

            outer_loop_nymd(i) = nymd_trj
            outer_loop_nhms(i) = nhms_trj

          endif ! < strajfrq >

        endif ! < trajfrq>0 >

!       Mimic internal loop of dynamics
!       -------------------------------
        do iinner = 1, ninner
          call tick( nymd_all,nhms_all,pdt )
        end do

      enddo

!     Release temporary memory
!     ------------------------
      deallocate ( pstmp, utmp, vtmp, delptmp, pttmp, qtmp )

      end subroutine all_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Init_:  Allocate arrays for all trajectory model states
!
! !INTERFACE:

      subroutine Init_ ( nymd_in, nhms_in,
     &                   memtrj, verbose )         ! optionals

! !USES:
      
      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: nymd_in   ! Date: year-month-day
      integer, intent(in) :: nhms_in   ! Time: hour-min-sec

      logical, optional, intent(in) :: memtrj   ! Determines trajectory status
      logical, optional, intent(in) :: verbose  ! Set verbose on/off

! !DESCRIPTION: Allocate arrays for all trajectory model state vectors.
!
! !REVISION HISTORY:
!
!  20Apr2007  Kim  Initial code (as per GetState_).
!  03May2007  Todling Option to read trajectory via i/o
!  05Jun2007  Todling Handle for trajfrq<0
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//"*init_"
      integer iouter,iinner,myid,ier
      integer nymd_trj, nhms_trj
      integer nymd_all, nhms_all
      integer nymd_fix, nhms_fix
      integer strajfrq, snhms

      if (trajfrq<0) return
 
      call MP_comm_rank(MP_comm_world,myID,ier)
         if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

      if (present(verbose) ) then
          verbose_ = verbose
      endif

!     If so, reset default trajectory status
!     --------------------------------------
      MEMTRJ_ = .true.
      if (present(memtrj) ) then
          MEMTRJ_ = memtrj
      endif

      if(.not.MEMTRJ_) return

!     Dummy parameter setups for test
!     -------------------------------
      nymd_t0 = nymd_in
      nhms_t0 = nhms_in

      nouter_trj=0
      nymd_all = nymd_in
      nhms_all = nhms_in
      do iouter = 1, nouter
        nymd_trj = nymd_all
        nhms_trj = nhms_all

        call tick(nymd_trj,nhms_trj,pdt)

        if (trajfrq>0) then
          strajfrq = nsecf(trajfrq)
          snhms    = nsecf(nhms_trj)
          if (mod(snhms,strajfrq)==0) then
              nouter_trj=nouter_trj+1
          endif
        endif

        do iinner = 1, ninner
          call tick(nymd_all,nhms_all,pdt)
        end do
      enddo

      if(myid==ROOT) print *,'The number of outer loop trajectories to be read in =', nouter_trj+1

!     Allocate five state variables for GetState_all
!     ----------------------------------------------
      allocate(u_all(imr,jfirst-ng_d:jlast+ng_s,kfirst:klast,0:nouter_trj))
      allocate(v_all(imr,jfirst-ng_s:jlast+ng_d,kfirst:klast,0:nouter_trj)) 
      allocate(pt_all(imr,jfirst-ng_d:jlast+ng_d,kfirst:klast,0:nouter_trj))  
      allocate(ps_all(imr,jfirst:jlast,0:nouter_trj) )
      allocate(delp_all(imr,jfirst:jlast,kfirst:klast,0:nouter_trj))          
      allocate(q_all(imr,jfirst-ng_d:jlast+ng_d,kfirst:klast,nc,0:nouter_trj))

      allocate(outer_loop_nymd(0:nouter_trj),outer_loop_nhms(0:nouter_trj))

      end subroutine Init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Retrieve_:  Read in model state with Getstate_all_
!
! !INTERFACE:

      subroutine Retrieve_ ( nymd_in, nhms_in, prog, doremap )

! !USES:

      implicit none

      include "comcon.h"               ! only : cappa

! !INPUT PARAMETERS:

      integer, intent(in) :: nymd_in, nhms_in

      type(dyn_prog), intent(inout) :: prog

! !OUTPUT PARAMETERS:

      logical, intent(out) :: doremap

! !DESCRIPTION: Reads in model state vector.
!
!  20Apr2007  Kim     Initial code (as per GetState_).
!  03May2007  Todling Option to read trajectory via i/o
!  18May2007  Todling Add verbose
!  05Jun2007  Todling Handle for trajfrq<0 case
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'::retrieve_'

      real(r8) ptopk
      integer  i, j, k, n, myid
      integer  ier, i_test, rmaptfrq
      integer  strajfrq, snhms
      character*255 :: fnametrj
      character(len=255) :: fname

      if (trajfrq<0) return

      call MP_comm_rank(MP_comm_world,myID,ier)
         if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

!     If trajectory not in memory ...
!     -------------------------------
      if ( .not. MEMTRJ_ ) then
         strajfrq = nsecf(trajfrq)
         snhms    = nsecf(nhms_in)
         if ( trajfrq>0 ) then
             if (mod(snhms,strajfrq)==0) then
                 call strTemplate( fnametrj, trajtmpl,
     &                             xid=trim(job), nymd=nymd_in, nhms=nhms_in )
                 fname = trim(fnametrj)
                 call GetState_ ( fname, job, nymd_in, nhms_in, prog, stat=ier )
                    if (ier/=0) call die (myname, 'cannot find trajectory file')
             endif
         endif

      else ! < MEMTRJ_ >

!        Identify trajectory locations      
!        -----------------------------
         call Trj_idx_(nymd_in,nhms_in,.false.,i)

         if ( verbose_ ) then
             if(myid==ROOT)  print *, trim(myname_), ': ', i, outer_loop_nymd(i), outer_loop_nhms(i)
         endif

!        Copy trajectory state variables
!        -------------------------------
         prog%delp(:,jfirst:jlast,:)             = delp_all(:,jfirst:jlast,:,i)
         prog%u   (:,jfirst-ng_d:jlast+ng_s,:)   = u_all   (:,jfirst-ng_d:jlast+ng_s,:,i)
         prog%v   (:,jfirst-ng_s:jlast+ng_d,:)   = v_all   (:,jfirst-ng_s:jlast+ng_d,:,i)
         prog%pt  (:,jfirst-ng_d:jlast+ng_d,:)   = pt_all  (:,jfirst-ng_d:jlast+ng_d,:,i)
         prog%q   (:,jfirst-ng_d:jlast+ng_d,:,1) = q_all   (:,jfirst-ng_d:jlast+ng_d,:,1,i)
              ps  (:,jfirst:jlast)               = ps_all  (:,jfirst:jlast,i)

         if ( nc .gt. 1 ) then
            if ( nc .le. w_f_grid_lm ) then
                do n=2,nc
                   prog%q(:,jfirst-ng_d:jlast+ng_d,:,n)  = q_all(:,jfirst-ng_d:jlast+ng_d,:,n,i)
                enddo
            else
                call die(myname_,'not enough tracers ')
            endif
         endif
      
      endif ! < MEMTRJ_ >

!     Decide on whether or not to remap
!     ---------------------------------
      rmaptfrq = nsecf(rmaptrj)
      snhms    = nsecf(nhms_in)
      doremap = .false.
      if (rmaptrj>0) then
          if (mod(snhms,rmaptfrq)==0) doremap = .true.
      endif

!     Recalculate phz, pe, pk, and pkz
!     NOTE2: openMP intentionally deactivated
!     ---------------------------------------
      ptopk = ptop ** cappa

!___$omp parallel do private(i,j,k)
      do j=jfirst,jlast
        do i=1,imr
           pe(i,1,j) = ptop
        enddo

        do k=2,nl+1
        do i=1,imr
           pe(i,k,j) = pe(i,k-1,j) + prog%delp(i,j,k-1)
        enddo
        enddo
      enddo

!___$omp parallel do default(shared) private(i,j,k)
        do j=jfirst,jlast
           do i=1,imr
              pk(i,j,1) = ptopk
           enddo

           do k=2,nl+1
              do i=1,imr
                 pk(i,j,k) = pe(i,k,j)**cappa
              enddo
           enddo
        enddo

        call pkez(1, imr, nl, jfirst, jlast, pe, pk, cappa, ks, peln, pkz, .false.)

      end subroutine Retrieve_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Clean_: Deallocate arrays w.r.t. model states of Getstate_all_
!
! !INTERFACE:

      subroutine clean_

      implicit none
! !DESCRIPTION: Deaalocate arrays w.r.t. model states of Getstate_all_
!
!  20Apr2007  Kim  Initial code (as per GetState_).
!  05Jun2007  Todling Handle for trajfrq<0
!
!EOP
!-----------------------------------------------------------------------

      if (trajfrq<0) return

      if(.not.MEMTRJ_) return

      deallocate ( u_all )
      deallocate ( v_all )
      deallocate ( pt_all )
      deallocate ( ps_all )
      deallocate ( delp_all )
      deallocate ( q_all )
      deallocate ( outer_loop_nymd, outer_loop_nhms )

      end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Trj_idx_: Counting trajectory location in memory
!
! !INTERFACE:

      subroutine Trj_idx_ ( nymd_in, nhms_in, load, frq_idx )

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: nymd_in, nhms_in
      logical, intent(in) :: load

! !OUTPUT PARAMETERS:

      integer, intent(out) :: frq_idx

! !DESCRIPTION: Reads in model state vector.
!
!  10May2007  Kim     Initial code
!  29May2007  Todling Revised to really tap on date/time to extract traj
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*trj_idx'
      integer    i,iouter,iinner
      integer :: nymd_all, nhms_all,nymd_trj, nhms_trj, strajfrq, snhms
      real(r8) :: calday_in, calday_trj

!     Initial time
!     ------------
      nymd_all = nymd_t0
      nhms_all = nhms_t0
      frq_idx=0

      if (nymd_in == nymd_t0 .and. nhms_in == nhms_t0) return

      if ( load ) then

           do iouter = 1, nouter
             nymd_trj = nymd_all
             nhms_trj = nhms_all
             call tick( nymd_trj,nhms_trj,pdt )

             if (trajfrq>0) then
               strajfrq = nsecf(trajfrq)
               snhms    = nsecf(nhms_trj)
               if (mod(snhms,strajfrq)==0) then
                   frq_idx=frq_idx+1
                   if (nymd_in == nymd_trj .and. nhms_in == nhms_trj) return
               endif
             endif

!          Mimic internal loop of dynamics
!          -------------------------------
           do iinner = 1, ninner
              call tick( nymd_all,nhms_all,pdt )
           end do
           enddo

      else

!          Identifying the location of current trajectory by Julian day
!          ------------------------------------------------------------
           frq_idx = -1
           do i = 0, nouter_trj
              call mcalday(nymd_in, nhms_in, calday_in)
              call mcalday(outer_loop_nymd(i), outer_loop_nhms(i), calday_trj)
              if ( calday_in .ge. calday_trj ) frq_idx=i
           enddo
           if(frq_idx<0)then
              print*, 'Date/Time/nouter searched for: ', nymd_in, nhms_in, nouter_trj
              call die(myname_,' cannot determine dyn trajectory index',99)
           endif

      endif



      end subroutine Trj_idx_

      end module m_IOState
