!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_trjphys --- Read-In Physics Trajectories
!
! !INTERFACE:
!
module m_trjphys

! !USES:

  use precision
  use prognostics

  use stepon, only : nsecf
  use stepon, only : tdt
  use stepon, only : nouter, ninner, pdt, ptrjtmpl, ptrjfrq, job

  use m_const, only : zvir

  use m_gfio_getfld
  use m_StrTemplate, only : strTemplate

  use m_mpif90, only : MP_comm_world
  use m_mpif90, only : MP_comm_rank

  use mod_comm, only :imr, jnp, nl
  use mod_comm, only : mp_scatter4d

  use m_die, only   : mp_die
  use m_die, only   : die
 
  implicit none

  PRIVATE

! !PUBLIC MEMBER FUNCTIONS:

  PUBLIC  physdrv1_get_all
  PUBLIC  physdrv1_get_init
  PUBLIC  physdrv1_get_clean
  PUBLIC  physdrv1_get_state

  interface physdrv1_get_all   ; module procedure all_      ; end interface
  interface physdrv1_get_init  ; module procedure init_     ; end interface
  interface physdrv1_get_clean ; module procedure clean_    ; end interface
  interface physdrv1_get_state ; module procedure retrieve_ ; end interface

!
! !DESCRIPTION: Initial code for reading physics trajectories all in once
!
! !REVISION HISTORY:
!  25Apr2007  Kim       read in all physics trajectories
!  12May2007  Todling   Fixed routine name convension 
!  01Jul2007  Todling   Some clean up
!
!EOP
!-------------------------------------------------------------------------

   integer, PUBLIC,  parameter :: nvars = 4   ! number of diffusion coefs
   character(len=3), parameter :: vname(nvars)=(/'cah','cam','cch','ccm'/)

   character(len=*), parameter :: myname = 'm_trjphys'

   real(kind=r8), PUBLIC, allocatable :: vcoef3d_all(:,:,:,:,:)  ! upper and lower diagonal elements of matrix

   integer, allocatable :: outer_nymd_phys(:)
   integer, allocatable :: outer_nhms_phys(:)
   integer :: nouter_trjphys

   integer :: State_ptrjfrq, nymd_t0, nhms_t0

   integer, parameter :: ROOT = 0        ! should come from above
   logical, save      :: MEMPHYS_ = .true.
   logical, save      :: verbose_ = .false.

CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: all_: Read all Physics Trajectories 
!
! !INTERFACE:
 
 subroutine all_ (nymd,nhms)
 
! !USES:

      implicit none
 
! !INPUT PARAMETERS:

      integer, intent(in) :: nymd ! Year-month-day, e.g.,  19971012
      integer, intent(in) :: nhms ! hour-minute-sec, e.g.,   120000
      

! !DESCRIPTION: 
                                                        
! !REVISION HISTORY:
!
!  23Apr2007  Kim      Initial code 
!  05Jun2007  Todling  Handle for ptrjfrq<0
!  01Jul2007  Todling  Read on root and scatter
!
!EOP
!-------------------------------------------------------------------------

! !LOCAL VARIALBLES

      character(len=*), parameter :: myname_=myname//'::all_'

      character*255 :: vFile ! name of file containing coefs
      integer i,j,k, j1, iouter, iinner, nymd_trj, nhms_trj, myid
      integer snhms, sptrjfrq
      integer ierr, iouter_iinner, i_test
      integer nymd_all, nhms_all      
      real(r8), allocatable :: aux(:,:,:,:)

      if(ptrjfrq<0) return

!     If trajectory not in memory ...
!     -------------------------------
      if ( .not. MEMPHYS_ ) return

      call MP_comm_rank(MP_comm_world,myID,ierr)
         if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)

!     Allocate temporary space
!     ------------------------
      if ( myid==ROOT ) then
           allocate( aux(imr,jnp,nl,nvars) )
      else
           allocate( aux(0,0,0,nvars) )
      endif

!     read state trajectories at initial time step (iouter=0)
!     -------------------------------------------------------
      nymd_all=nymd
      nhms_all=nhms
      call strTemplate( vFile, ptrjtmpl, &
                        xid=trim(job), nymd=nymd_all, nhms=nhms_all )

!     Read trajectory on root processor
!     ---------------------------------
      if ( myid==ROOT ) then
           call GFIO_GetFld( vFile, nymd_all, nhms_all, imr, jnp, nl, &
     		             nvars, vname, aux, stat=ierr )
           aux = tdt * aux
      endif ! < myid==root >

!     Scatter trajectory
!     ------------------    
      do k = 1, nvars
         call mp_scatter4d(aux(:,:,:,k:k), vcoef3d_all(:,:,:,k:k,0), imr, jnp, nl, 1, jfirst, jlast, 1, nl, 0, 0, 0)
      enddo
      outer_nymd_phys(0)=nymd_all
      outer_nhms_phys(0)=nhms_all

!     read the rest of trajectories
!     -----------------------------
      nymd_all=nymd
      nhms_all=nhms
      iouter_iinner=0
      do iouter = 1, nouter
        do iinner = 1, ninner
	   iouter_iinner=iouter_iinner+1
           call tick( nymd_all,nhms_all,pdt )

	   if ( ptrjfrq == 0 ) then
              print *, trim(myname_),': this option not yet implemented '
	      call exit(7)
	   else
	      sptrjfrq = nsecf(ptrjfrq)
	      snhms    = nsecf(nhms_all)
	      if (mod(snhms,sptrjfrq)==0) then
                  call strTemplate( vFile, ptrjtmpl, &
                              xid=trim(job), nymd=nymd_all, nhms=nhms_all )

!                 Read trajectory on root processor
!                 ---------------------------------
                  if ( myid==ROOT ) then

                       call GFIO_GetFld ( vFile, nymd_all, nhms_all, imr, jnp, nl, &
             		 		  nvars, vname, aux, stat=ierr )
               	            if (ierr/=0) call mp_die ( myname_, 'cannot read phys trajectory',99 )

	               aux = tdt * aux

                  endif ! < myid==root >

                  call Phys_idx_(nymd_all,nhms_all,.true.,i)

!                 Scatter trajectory
!                 ------------------
                  do k = 1, nvars
                     call mp_scatter4d(aux(:,:,:,k:k), vcoef3d_all(:,:,:,k:k,i), imr, jnp, nl, 1, jfirst, jlast, 1, nl, 0, 0, 0)
                  enddo

		  outer_nymd_phys(i)=nymd_all
		  outer_nhms_phys(i)=nhms_all
              endif
           endif
 
        end do !do iinner
      enddo !do iouter

!     Release temporary space
!     -----------------------
      deallocate ( aux )

end subroutine all_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Init_:  Initialization
!
! !INTERFACE:

   subroutine init_( nymd_in, nhms_in, memphys, verbose )

! !USES:

! !INPUT PARAMETERS:

      integer,           intent(in) :: nymd_in  ! Year-month-day, e.g.,  19971012
      integer,           intent(in) :: nhms_in  ! hour-minute-sec, e.g.,   120000
      logical, optional, intent(in) :: memphys  ! Determines trajectory status
      logical, optional, intent(in) :: verbose  ! Sets verbose on/off

! !DESCRIPTION: Reads in model state vector.
!
! !REVISION HISTORY:
!
!  20Apr2007  Kim     Initial code
!  18May2007  Todling Add verbose handle
!  05Jun2007  Todling Handle for ptrjfrq<0
!  01Jul2007  Todling Simplified: 2nd dim of vcoef need only jf to jl
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_=myname//'::init_'

      integer :: nymd_all, nhms_all      
      integer :: nymd_fix, nhms_fix
      integer snhms, sptrjfrq
      integer ierr, iouter, iinner, nouter_trj

      if(ptrjfrq<0) return

!     Get initial time for counting trajectories in memory
!     ----------------------------------------------------
      nymd_t0 = nymd_in
      nhms_t0 = nhms_in

      if (present(verbose) ) then
          verbose_ = verbose
      endif

!     Set/Reset trajectory status
!     ---------------------------
      if (present(memphys) ) then
          MEMPHYS_ = memphys
      endif
      if(.not.MEMPHYS_) return

      nouter_trjphys=0

      nymd_all=nymd_in
      nhms_all=nhms_in

      do iouter = 1, nouter
        do iinner = 1, ninner
           call tick( nymd_all,nhms_all,pdt )

	   if ( ptrjfrq == 0 ) then
              print *, trim(myname_), ': this option not yet implemented '
	      call exit(7)
	   else
	      sptrjfrq = nsecf(ptrjfrq)
	      snhms    = nsecf(nhms_all)
	      if (mod(snhms,sptrjfrq)==0) then
                  nouter_trjphys=nouter_trjphys+1
              endif
           endif
 
        end do !do iinner
      enddo !do iouter

!     Allocate the trajectory in memory
!     ---------------------------------
      allocate ( vcoef3d_all(imr, jfirst:jlast, nl, nvars, 0:nouter_trjphys), stat=ierr )
      allocate(outer_nymd_phys(0:nouter_trjphys))
      allocate(outer_nhms_phys(0:nouter_trjphys))
      if(ierr/=0) call die(myname,' Alloc(vcoef3d_all) error',ierr)

   end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Retrieve_: Retrieving trajectrories from memory
!
! !INTERFACE:

   subroutine retrieve_(vcoef3d,nymd_in,nhms_in)

! !USES:

! !INPUT PARAMETERS:
      integer, intent(in) :: nymd_in,nhms_in

! !OUTPUT PARAMETERS:
      real(kind=r8), intent(out) :: vcoef3d(imr, jfirst:jlast, nl, nvars)

! !DESCRIPTION: Reads in physics coefficients as part of trajectory.
!
! !REVISION HISTORY:
!
!  20Apr2007  Kim      Initial code
!  05Jun2007  Todling  Handle for ptrjfrq<0
!  01Jul2007  Todling  Simplified: 2nd dim of vcoef need only jf to jl
!  11Dec2007  Todling  Fix for case of when traj not in memory
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'::retrieve_'

      character*255 :: vFile        ! name of file containing coefs
      integer :: iouter_phys
      integer :: ierr, i_test, myid, k
      logical :: verb
      real(r8), allocatable :: aux(:,:,:,:)

      if(ptrjfrq<0) return

      call MP_comm_rank(MP_comm_world,myID,ierr)
         if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)

!     If trajectory not in memory ...
!     -------------------------------
      verb=.false.; if(myid==ROOT) verb=.true.
      if ( .not. MEMPHYS_ ) then

!        Allocate temporary space
!        ------------------------
         if ( myid==ROOT ) then
              allocate( aux(imr,jnp,nl,nvars) )
         else
              allocate( aux(0,0,0,nvars) )
         endif

         call strTemplate( vFile, ptrjtmpl, &
                           xid=trim(job), nymd=nymd_in, nhms=nhms_in )

!        Read trajectory on root processor
!        ---------------------------------
         if ( myid==ROOT ) then
              call GFIO_GetFld( vFile, nymd_in, nhms_in, imr, jnp, nl, nvars, vname, aux, stat=ierr )
              aux = tdt * aux
         endif ! < myid==root >
                                                                                                                      
!        Scatter trajectory
!        ------------------
         do k = 1, nvars
            call mp_scatter4d(aux(:,:,:,k:k), vcoef3d(:,:,:,k:k), imr, jnp, nl, 1, jfirst, jlast, 1, nl, 0, 0, 0)
         enddo

         deallocate ( aux )
         return
      endif

      call Phys_idx_(nymd_in,nhms_in,.false.,iouter_phys)

      if ( verbose_ ) then
           if(myid==ROOT) print*,trim(myname_), ': ',iouter_phys,&
                                 outer_nymd_phys(iouter_phys),outer_nhms_phys(iouter_phys)
      endif

      vcoef3d(:,:,:,:)=vcoef3d_all(:,:,:,:,iouter_phys)

   end subroutine retrieve_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Clean_: Deallocate arrays
!
! !INTERFACE:

   subroutine clean_
! !DESCRIPTION: Deaalocate arrays
!
!  20Apr2007  Kim  Initial code
!  05Jun2007  Todling Handle for ptrjfrq<0
!
!EOP
!-----------------------------------------------------------------------

      if(ptrjfrq<0) return

!     If trajectory not in memory ...
!     -------------------------------
      if ( .not. MEMPHYS_ ) return

      deallocate (vcoef3d_all)
      deallocate (outer_nymd_phys)
      deallocate (outer_nhms_phys)
   end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Phys_idx: Counting trajectry location in memory
!
! !INTERFACE:

   subroutine Phys_idx_(nymd_in,nhms_in,load,frq_idx)

!USES:

   implicit none

! !INPUT PARAMETERS:

   integer, intent(in)  :: nymd_in, nhms_in
   logical, intent(in)  :: load

! !OUTPUT PARAMETERS:

   integer, intent(out) :: frq_idx

! !DESCRIPTION: Reads in model state vector.
!
!  10May2007  Kim     Initial code (as per GetState_).
!  29May2007  Todling Revised to really tap on date/time to extract traj
!
!EOP
!-----------------------------------------------------------------------

   character(len=*), parameter :: myname_ = myname//'::Phys_idx'

   integer nymd_all, nhms_all,nymd_trj, nhms_trj, sptrjfrq, snhms
   integer i,iouter, iinner
   real(r8) calday_in,calday_trj

!  Initial time
!  ------------
   nymd_all = nymd_t0
   nhms_all = nhms_t0
   frq_idx=0

   if (nymd_in == nymd_t0 .and. nhms_in == nhms_t0) return

   if ( load ) then
        do iouter = 1, nouter
           do iinner = 1, ninner
              call tick( nymd_all,nhms_all,pdt )
              if ( ptrjfrq == 0 ) then
                   print *, trim(myname_), ': this option not yet implemented '
     	           call exit(7)
              else
	           sptrjfrq = nsecf(ptrjfrq)
	           snhms    = nsecf(nhms_all)
	           if (mod(snhms,sptrjfrq)==0) then
                       frq_idx=frq_idx+1
                       if (nymd_in == nymd_all .and. nhms_in == nhms_all) return     
                   endif
              endif
           end do !do iinner
        enddo !do iouter

   else

        frq_idx = -1
        do i = 1, nouter_trjphys
           call mcalday(nymd_in, nhms_in, calday_in)
           call mcalday(outer_nymd_phys(i), outer_nhms_phys(i), calday_trj)
           if ( calday_in .ge. calday_trj ) frq_idx=i
        enddo
        if(frq_idx<0) call die(myname_,' cannot determine phys trajectory index',99)

   endif

   end subroutine Phys_idx_

   end module m_trjphys
