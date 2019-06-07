!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_TrajMng --- Trajectory Management Module
!
! !INTERFACE:
!
      module m_trajmng

! !USES:

      use mod_comm, only : gid,mp_gather4d,mp_scatter2d,mp_scatter4d,mp_bcst_int
      use mod_comm, only : imr         ! no. of grid points in the zonal direction
      use mod_comm, only : jnp         ! no. of grid points in the meridional direction
      use mod_comm, only : nl          ! no. of levels
      use mod_comm, only : nc          ! no. of tracers

      use stepon, only : ng_d, ng_s
 
      use precision
      use prognostics, only : jfirst   ! mpi begin lat for this pe
      use prognostics, only : jlast    ! mpi end   lat for this pe
      use prognostics, only : sgh      ! topography height stdv (meter)
      use prognostics, only : phis     ! topography geopotential (g*z_s)

      use prognostics, only : dyn_prog
      use prognostics_q, only : dynq

      use m_dyn
      use m_inpak90
      use m_die, only : MP_die

      use m_GFIO_Putfld, only : GFIO_PutFld
      use m_GFIO_Getfld, only : GFIO_GetFld

      use timingModule
                                                                                                                       
      use stepon,  only : ptop
      use stepon,  only : ks
      use stepon,  only : ak
      use stepon,  only : bk
      use stepon,  only : ts
      use stepon,  only : oro
                                                                                                                       
      use m_mpif90, only : MP_comm_world
      use m_mpif90, only : MP_comm_rank

      use m_dyn
      use m_FileResolv                 ! file resolver
      use m_die,  only : die

      implicit none


! !PUBLIC MEMBER FUNCTIONS:

      PRIVATE

      PUBLIC PutPert
      PUBLIC GetPert
      PUBLIC InqPert_dims
 
      interface PutPert; module procedure
     .          PutPert1_,
     .          PutPert2_,
     .          PutPert3_
      end interface

      interface GetPert; module procedure
     .          GetPert1_,
     .          GetPert2_,
     .          GetPert3_
      end interface
      interface inqpert_dims; module procedure
     .          dims_
      end interface
!
! !DESCRIPTION: Controls I/O of GFIO-type trajectory used in TLM
!               and ADM codes. 
!
! !REVISION HISTORY:
!
!  19Nov2002  Todling   Initial design/interfaces.
!  09Jun2006  Elena N.  Bug fixed in GetPert2: Added Dyn_Init call to define vector size
!                                  Replaced Dyn_Clean by Dyn_Null to release pointers
!  10Jan2007  Todling   Merged in Oloso's changes
!  14May2007  Todling   Introduced dyn_prog; global change.
!  17Jul2007  Todling   Add routine to inquire dimension of perturbation fields from file
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'm_TrajMng'
      character(len=*), parameter :: trajfn_def = 'trajectory'

#if !defined( SPMD )
      integer, parameter :: gid = 0
#endif
      integer, parameter :: ROOT = 0     ! should come from above

      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: PutPert1_:  Write out perturbations to GFIO file
!
! !INTERFACE:

      subroutine PutPert1_ ( job, nymd, nhms, prog, freq, nstep,
     &                       ak, bk, Ts, oro, ps,
     &                       fntmpl, fntype, 
     &                       nymd_b, nhms_b, nymd_e, nhms_e,
     &                       stat, vectype, forceflip )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      integer, optional :: nymd_b
      integer, optional :: nhms_b
      integer, optional :: nymd_e
      integer, optional :: nhms_e

      character(len=*), intent(in) :: job
      integer,          intent(in) :: nymd   ! Date: year-month-day
      integer,          intent(in) :: nhms   ! Time: hour-min-sec
 
      type(dyn_prog),   intent(inout) :: prog ! prognostic variables

      integer,          intent(in) :: freq   ! Analysis frequency (HHMMSS)
      integer,          intent(in) :: nstep  ! model time step

      real,    intent(in)  ::  ak(nl+1)      ! vertical grid a-coefficient
      real,    intent(in)  ::  bk(nl+1)      ! vertical grid b-coefficient
      real,    intent(in)  ::  Ts(imr,jfirst:jlast)   ! Sea surface temperature (K)
      real,    intent(in)  ::  oro(imr,jfirst:jlast)  ! Land-water-ice mask:
                                             !   oro = 0  over ocean
                                             !   oro = 1  over land
                                             !   oro = 2  over sea ice
      real,    intent(in)  ::  ps(imr,jfirst:jlast) ! surface pressure

      integer, optional :: stat
      character(len=*), intent(in), optional :: fntmpl
      character(len=*), intent(in), optional :: fntype
      integer, intent(in), optional :: vectype
      logical, intent(in), optional :: forceflip

! !DESCRIPTION: Reads in model state vector.
!
! !REVISION HISTORY:
!
!  22Nov2002  Todling   Initial code (as per ana.f).
!  22Nov2006  Oloso     Ifdef'd timinig routines for SPMD
!  28Aug2007  Todling   - Add vectype in support of geos-5 perturbations
!                       - Updated filename convention to GEOS-5
!  30Oct2007  Todling   Changed template definition
!  29Jan2009  Todling   Replace mcalday with julian
!  06Mar2009  Todling   Update dyn_init interface
!  16Jan2010  Todling   Add forceflip
!  19Jan2010  Todling   Simplified put for SPMD
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*PutPert1_'

      character(len=*), parameter :: fname_def = 'perts.eta'

      character(len=80)  :: fname
      character(len=2)   :: mark
      integer, parameter :: prec = 0  ! no need to output pert in 64 bits
      integer b_date,e_date
      integer, external :: julday

      real(r8), allocatable ::   ts_x(:,:)      ! Sea surface temperature (K)
      real(r8), allocatable ::  oro_x(:,:)      ! Land-water-ice mask:
      real(r8), allocatable ::  sgh_x(:,:)      ! topography height stdv (meter)
      real(r8), allocatable :: phis_x(:,:)      ! topography geopotential (g*z_s)
      real(r8), allocatable ::   ps_x(:,:)      ! surface pressure (pascal)
      real(r8), allocatable :: delp_x(:,:,:)    ! pressure thickness
      real(r8), allocatable ::    u_x(:,:,:)    ! zonal wind
      real(r8), allocatable ::    v_x(:,:,:)    ! meridional wind
      real(r8), allocatable ::   pt_x(:,:,:)    ! scaled virtual potential temperature
      real(r8), allocatable ::    q_x(:,:,:,:)  ! specific humidity & tracers

      type (dyn_vect) :: w_f
      real(r8), allocatable :: aux2d(:,:)  ! aux array for land/etc fractions

      integer  j,lm, ierr
      integer  nymd_v, nhms_v
      integer  hour, hour_b, hour_e
      integer  yrb,mmb,ddb
      integer  yre,mme,dde
      character(len=8) str

      if(present(stat)) stat = 0

      lm = max(size(prog%q,4),nc)

!     Allocate local arrays
!     ---------------------
      allocate(   ts_x(imr,jnp)   , oro_x(imr,jnp) ,  sgh_x(imr,jnp),
     &          phis_x(imr,jnp)   ,  ps_x(imr,jnp) , delp_x(imr,jnp,nl),
     &             u_x(imr,jnp,nl), v_x(imr,jnp,nl),   pt_x(imr,jnp,nl),
     &             q_x(imr,jnp,nl,lm), stat=ierr )
        if ( ierr/=0 ) then
           if(present(stat))then
              stat = 1
              return
           else
             call MP_die (myname_,'Error from Alloc',ierr)
           endif
        end if

!     Assign pointers in w_f input arrays
!     -----------------------------------
#if defined( SPMD )
           call mp_gather4d(phis, phis_x, imr, jnp,  1,  1, jfirst, jlast, 1,  1, 0, 0, 0)
           call mp_gather4d(sgh ,  sgh_x, imr, jnp,  1,  1, jfirst, jlast, 1,  1, 0, 0, 0)
           call mp_gather4d(ts  ,   ts_x, imr, jnp,  1,  1, jfirst, jlast, 1,  1, 0, 0, 0)
           call mp_gather4d(oro ,  oro_x, imr, jnp,  1,  1, jfirst, jlast, 1,  1, 0, 0, 0)
           call mp_gather4d(ps  ,   ps_x, imr, jnp,  1,  1, jfirst, jlast, 1,  1, 0, 0, 0)
           call mp_gather4d(prog%u   ,    u_x, imr, jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_s, 0)
           call mp_gather4d(prog%v   ,    v_x, imr, jnp, nl,  1, jfirst, jlast, 1, nl, ng_s, ng_d, 0)
           call mp_gather4d(prog%pt  ,   pt_x, imr, jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
           call mp_gather4d(prog%delp, delp_x, imr, jnp, nl,  1, jfirst, jlast, 1, nl,    0,    0, 0)
           call mp_gather4d(prog%q   ,    q_x, imr, jnp, nl, lm, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
#else
           phis_x = phis
           sgh_x  = sgh
           ts_x   = ts
           oro_x  = oro
           ps_x   = ps
           delp_x = prog%delp
           u_x    = prog%u
           v_x    = prog%v
           pt_x   = prog%pt
           q_x    = prog%q
#endif

      if ( gid == 0 ) then

!     Allocate aux array
!     ------------------
      allocate(aux2d(imr,jnp))
      aux2d=0.d0  ! can't do much better for land/etc fractions for now

!     Make sure w_f does not point to anything
!     ----------------------------------------
      call dyn_null ( w_f )

!     Assign pointers in w_f input arrays
!     -----------------------------------
      call Dyn_Init ( imr, jnp, nl, lm, ptop, ks, ak, bk,
     &                phis_x, sgh_x, Ts_x, oro_x, ps_x,
     &                aux2d, aux2d, aux2d, aux2d, aux2d,
     &                delp_x, u_x,v_x, pt_x,q_x,
     &                w_f, ierr, vectype=vectype  )
        if ( ierr .ne. 0 ) then
           if(present(stat))then
              stat = 1
              return
           else
             call MP_die (myname_,'Error from Dyn_Init',ierr)
           endif
        end if

      nymd_v = nymd
      nhms_v = nhms

!     Name output file
!     ----------------
      if (present(fntmpl)) then
          if (present(nymd_b).and.present(nhms_b).and.
     &        present(nymd_e).and.present(nhms_e)    ) then
              write(str,'(i8.8)') nymd_b
              read(str,'(i4,2i2)') yrb,mmb,ddb
              write(str,'(i8.8)') nymd_e
              read(str,'(i4,2i2)') yre,mme,dde
              b_date = julday ( mmb, ddb, yrb )
              e_date = julday ( mme, dde, yre )
              mark = 'z+'
              if(b_date>e_date) mark = 'z-'
              hour_b = nhms_b/10000
              hour_e = nhms_e/10000
              if ( present(fntype) ) then
                 if( gid == 0 ) 
     &           write(fname,'(a,2(a,i8.8,a,i2.2),3a)') trim(fntmpl   ),
     &                                                   '.',       nymd_b, '_', hour_b,
     &                                                  trim(mark), nymd_e, '_', hour_e, 'z.',
     &                                                  trim(fntype), '.nc4' 
              else
                 if( gid == 0 ) 
     &           write(fname,'(a,2(a,i8.8,a,i2.2),a)')  trim(fntmpl   ),
     &                                                   '.',       nymd_b, '_', hour_b,
     &                                                  trim(mark), nymd_e, '_', hour_e, 'z.nc4' 
              endif
          else
              hour = nhms_v/10000
              if ( present(fntype) ) then
                   if ( gid==0 ) 
     &             write(fname,'(2a,i8.8,a,i2.2,3a)') trim(fntmpl), '.',
     &                                                nymd_v, '_', hour, 'z.',
     &                                                trim(fntype), '.nc4' 
              else
                   if ( gid==0 ) 
     &             write(fname,'(2a,i8.8,a,i2.2,a)')  trim(fntmpl   ), '.',
     &                                                nymd_v, '_', hour, 'z.nc4' 
              endif
          end if
      else     
!         hour = nhms_v/10000
!         write(fname,'(2a,i8.8,a)') trim(fname_def), 
!    &                                      '.', nymd_v, '.nc4' 
          fname = trim(job)
      end if

!     Write out
!     ---------
      call dyn_put ( fname, nymd, nhms, prec, w_f, ierr, freq=freq, nstep=nstep, vectype=vectype, forceflip=forceflip )
	if(ierr/=0) call MP_die (myname_,'Error from dyn_put',ierr)
        print *, myname_//': wrote perturbation file '//trim(fname)
        print *

!     Release pointers
!     ----------------
      call dyn_null ( w_f )
      deallocate(aux2d)

      endif  ! < gid >

      deallocate( ts_x, oro_x, sgh_x, phis_x, ps_x,
     &            delp_x, u_x, v_x, pt_x, q_x, stat=ierr ) 

      end subroutine PutPert1_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PutPert2_ --- 
!
! !INTERFACE:

      subroutine PutPert2_ ( job, nymd, nhms, prog, nstep,
     &                       ak, bk, Ts, oro, ps,
     &                       rcfname, stat, vectype, forceflip )


! !USES:

      use precision
      use prognostics, only : imr      ! no. of grid points in the zonal direction
      use prognostics, only : jnp      ! no. of grid points in the meridional direction
      use prognostics, only : nl       ! no. of levels

      use m_StrTemplate, only : strTemplate  ! grads style template
      use m_inpak90
      use m_die, only : die

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: job
      integer,          intent(in) :: nymd   ! Date: year-month-day
      integer,          intent(in) :: nhms   ! Time: hour-min-sec

      type(dyn_prog),   intent(inout) :: prog

      integer,          intent(in) :: nstep  ! model time step

      real,    intent(in)  ::  ak(nl+1)      ! vertical grid a-coefficient
      real,    intent(in)  ::  bk(nl+1)      ! vertical grid b-coefficient
      real,    intent(in)  ::  Ts(imr,jfirst:jlast)   ! Sea surface temperature (K)
      real,    intent(in)  ::  oro(imr,jfirst:jlast)  ! Land-water-ice mask:
                                             !   oro = 0  over ocean
                                             !   oro = 1  over land
                                             !   oro = 2  over sea ice
      real,    intent(in)  ::  ps(imr,jfirst:jlast)  ! surface pressure

      character(len=*), intent(in) :: rcfname
      integer, intent(in),  optional :: vectype
      logical, intent(in),  optional :: forceflip


! !OUTPUT PARAMETERS:

      integer, intent(out), optional :: stat

! !DESCRIPTION: Writes out state vector as a GFIO file.
!
! !REVISION HISTORY:
!
!  25Feb2004  Todling    Added this interface
!  28Aug2007  Todling    Added vectype in support of geos-5 perturbations
!  16Jan2010  Todling    Add forceflip
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*PutPert2_'

      character(len=255) :: fname          ! trajectory filename
      integer            :: freq           ! trajectory output freq

      integer strajfrq, snhms, ierror

      integer nsecf, n
      nsecf(n) = n/10000*3600 + mod(n,10000)/100* 60 + mod(n,100)

      if(present(stat)) stat = 0

!     Get info on filename and frequency from resouce file
!     ----------------------------------------------------
      call putinit_ (ierror)
           if(ierror/=0)return

      if ( freq .gt. 0 ) then
        strajfrq = nsecf(freq)
        snhms    = nsecf(nhms)
        if ( mod(snhms,strajfrq) == 0 ) then
          call PutPert1_ ( fname, nymd, nhms, prog, freq, nstep,
     &                     ak, bk, Ts, oro, ps, vectype=vectype,
     &                     forceflip=forceflip )
        endif
      endif

      CONTAINS

      subroutine putinit_ ( stat )

      integer  stat

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
      fname_def = trim(job) // '.pts.eta.' // cymd // '.nc4'

      fname  = trim(fname_def)
      freq   = freq_def

!     Load resources
!     -------------
      call i90_loadf ( trim(rcfname), ier )
        if ( ier/=0 ) return   ! if rc not found, return silently

        call i90_label ( 'pt_statevector_filename:', ier )
        if ( ier .ne. 0 ) then
           if(gid==0) print *, 'taking nlm state-vec  default filename '
        else
          call i90_gtoken ( template, ier )
             if(ier.ne.0) call die (myname, 'premature end of rc file')
             call strTemplate ( fname, template, 'GRADS',
     &                          xid=trim(job), nymd=nymd, nhms=nhms )
        endif

        call i90_label ( 'pt_statevector_frequency:', ier )
        if ( ier .ne. 0 ) then
           if(gid==0) print *, 'taking nlm state-vec i/o default default '
        else
          freq = i90_GInt( ier )
            if(ier.ne.0) call die (myname, 'premature end of rc file')
        endif

        call i90_release()

      end subroutine putinit_

      end subroutine PutPert2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: PutPert3_:  Write out tracer perturbations to GFIO file
!
! !INTERFACE:

      subroutine PutPert3_ ( job, nymd, nhms, qpert, freq, nstep,
     &                       fntmpl, fntype, 
     &                       nymd_b, nhms_b, nymd_e, nhms_e,
     &                       stat, forceflip, RCfile )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      integer, optional :: nymd_b
      integer, optional :: nhms_b
      integer, optional :: nymd_e
      integer, optional :: nhms_e

      character(len=*), intent(in) :: job
      integer,          intent(in) :: nymd   ! Date: year-month-day
      integer,          intent(in) :: nhms   ! Time: hour-min-sec
 
      type(dynq),   intent(inout) :: qpert   ! prognostic variables

      integer,          intent(in) :: freq   ! Analysis frequency (HHMMSS)
      integer,          intent(in) :: nstep  ! model time step

      integer, intent(out), optional :: stat
      logical, intent(in),  optional :: forceflip
      character(len=*), intent(in), optional :: fntmpl
      character(len=*), intent(in), optional :: fntype
      character(len=*), intent(in), optional :: RCfile

! !DESCRIPTION: Writes out tracer perturbations
!
! !REVISION HISTORY:
!
!  25Sep2007  Todling   Overload Put to handle qpert.
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*PutPert3_'

      character(len=*), parameter :: fname_def = 'qpert.eta'

      character(len=80)  :: fname
      integer, parameter :: prec = 0  ! not need to output pert in 64 bits

      real(r8), allocatable :: q_x(:,:,:,:)  ! specific humidity & tracers

      integer  j,n,lm, ierr
      integer  nymd_v, nhms_v
      integer  hour, hour_b, hour_e

      character(len=255) vname (nc)
      character(len=255) vunits(nc)
      character(len=255) vtitle(nc)

      character(len=*), parameter :: title  = "Tracers File"
      character(len=*), parameter :: source = "Global Modeling and Assimilation Office"

      if(present(stat)) stat = 0

      lm = max(size(qpert%q,4),nc)

!     Allocate local arrays
!     ---------------------
      allocate( q_x(imr,jnp,nl,nc), stat=ierr )
        if ( ierr/=0 ) then
           if(present(stat))then
              stat = 1
              return
           else
             call MP_die (myname_,'Error from Alloc',ierr)
           endif
        end if

#if defined( SPMD )
      call mp_gather4d(qpert%q, q_x, imr, jnp, nl, lm, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
#else
      q_x = qpert%q
#endif

!     Fill auxiliar informational arrays
!     ----------------------------------
      call readRC_ ( nc, vname, vunits, vtitle, ierr, RCfile )
        if ( ierr/=0 ) then
           if(present(stat))then
              stat = ierr
              return
           else
             call MP_die (myname_,'Error from readRC_',ierr)
           endif
        end if

      if ( gid == 0 ) then

      nymd_v = nymd
      nhms_v = nhms

!     Name output file
!     ----------------
      if (present(fntmpl)) then
          if (present(nymd_b).and.present(nhms_b).and.
     &        present(nymd_e).and.present(nhms_e)    ) then
              hour_b = nhms_b/10000
              hour_e = nhms_e/10000
              if ( present(fntype) ) then
                 if( gid == 0 ) 
     &           write(fname,'(a,2(a,i8.8,a,i2.2),3a)') trim(fntmpl   ),
     &                                                  '.', nymd_b, '_', hour_b,
     &                                                  'z+', nymd_e, '_', hour_e, 'z.',
     &                                                  trim(fntype), '.nc4' 
              else
                 if( gid == 0 ) 
     &           write(fname,'(a,2(a,i8.8,a,i2.2),a)')  trim(fntmpl   ),
     &                                                  '.', nymd_b, '_', hour_b,
     &                                                  'z+', nymd_e, '_', hour_e, 'z.nc4' 
              endif
          else
              hour = nhms_v/10000
              if ( present(fntype) ) then
                   if ( gid==0 ) 
     &             write(fname,'(2a,i8.8,a,i2.2,3a)') trim(fntmpl), '.',
     &                                                nymd_v, '_', hour, 'z.',
     &                                                trim(fntype), '.nc4' 
              else
                   if ( gid==0 ) 
     &             write(fname,'(2a,i8.8,a,i2.2,a)')  trim(fntmpl   ), '.',
     &                                                nymd_v, '_', hour, 'z.nc4' 
              endif
          end if
      else     
!         hour = nhms_v/10000
!         write(fname,'(2a,i8.8,a)') trim(fname_def), 
!    &                                      '.', nymd_v, '.nc4' 
          fname = trim(job)
      end if

!     Write out
!     ---------
      call GFIO_PutFld ( title,  source,  fname,
     &                   nymd,   nhms,     freq,
     &                   imr,     jnp,       nl,
     &                   ptop,     ks,  ak,  bk,
     &                   nc, vname, vtitle, vunits, q_x, 
     &                   ierr, untag=.true., forceflip=forceflip )
	 if(ierr/=0) call MP_die (myname_,'Error from GFIO_PutFld',ierr)
         print *, myname_//': wrote perturbation file '//trim(fname)
         print *

      endif  ! < gid >

      deallocate( q_x, stat=ierr ) 

      end subroutine PutPert3_

      subroutine readRC_ ( nvars, vname, vunits, vtitle, ierr, RCfile )

      implicit none
      integer,          intent(in)  :: nvars
      character(len=*), intent(out) :: vname(nvars)
      character(len=*), intent(out) :: vunits(nvars)
      character(len=*), intent(out) :: vtitle(nvars)
      character(len=*), optional :: RCfile
      integer          ierr
      integer          iret,n
      character*3      ntracer
      character(len=255) mylabel
      ierr = 0

!     Load resources
!     --------------
      if ( present(RCfile) ) then
         call i90_loadf (RCfile, iret)
           if (iret/=0) then
               write(6,'(2a,i5,2a)') myname, ': I90_loadf error, iret =', iret,
     .                                       ': trying to load ', trim(RCfile)
               ierr=iret
           end if

!        Read tracer information in a dummy way - transform to table later
!        -----------------------------------------------------------------
         do n = 1, nvars
            write(ntracer, '(i3.3)') n
   
!           Read tracer name
!           ----------------
            mylabel = 'short_name_tracer'//ntracer//':'
            call I90_label(trim(mylabel), iret)
            if (iret/=0) then
               write(6,'(2a,i5,2a)') myname, ': I90_label error, iret=', iret,
     .                                       ': trying to read ', 'tracer name:'
               ierr=iret
            end if
            call I90_GToken(vname(n),iret)
            if (iret/=0) then
               write(6,'(2a,i5)') myname, ': I90_GToken error, iret=', iret
               ierr=iret
            end if

!           Read tracer description
!           -----------------------
            mylabel = 'long_name_tracer'//ntracer//':'
            call I90_label(trim(mylabel), iret)
            if (iret/=0) then
               write(6,'(2a,i5,2a)') myname, ': I90_label error, iret=', iret,
     .                                       ': trying to read ', 'tracer name:'
               ierr=iret
            end if
            call I90_GToken(vtitle(n),iret)
            if (iret/=0) then
               write(6,'(2a,i5)') myname, ': I90_GToken error, iret=', iret
               ierr=iret
            end if

!           Read tracer description
!           -----------------------
            mylabel = 'units_tracer'//ntracer//':'
            call I90_label(trim(mylabel), iret)
            if (iret/=0) then
               write(6,'(2a,i5,2a)') myname, ': I90_label error, iret=', iret,
     .                                        ': trying to read ', 'tracer name:'
               ierr=iret
            end if
            call I90_GToken(vunits(n),iret)
            if (iret/=0) then
               write(6,'(2a,i5)') myname, ': I90_GToken error, iret=', iret
               ierr=iret
            end if

         enddo ! < nc >

!        Release resource file
!        ---------------------
         call i90_release()

      else ! Case without RCfile

         do n = 1, nvars
            write(ntracer,'(i3.3)') n
            vname(n)  = 'q'//ntracer
            vunits(n) = "kg/kg"
            vtitle(n) = 'Tracer '//ntracer
         enddo
 
      endif

      end subroutine readRC_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetPert1_:  Read in perturbation state vector
!
! !INTERFACE:

      subroutine GetPert1_ ( infile, 
     &                      job, nymd, nhms,
     &                      im, jm, km, lm, jfirst, jlast,
     &                      delp_pt, pt_pt, u_pt, v_pt, q_pt,
     &                      stat )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      integer,          intent(in) :: im
      integer,          intent(in) :: jm
      integer,          intent(in) :: km
      integer,          intent(in) :: lm
      integer,          intent(in) :: jfirst
      integer,          intent(in) :: jlast

      real,             intent(out) :: delp_pt(im,jm,km)
      real,             intent(out) :: pt_pt  (im,jm,km)
      real,             intent(out) :: u_pt   (im,jm,km)
      real,             intent(out) :: v_pt   (im,jm,km)
      real,             intent(out) :: q_pt   (im,jm,km,lm)

      character(len=*), intent(in) :: infile
      character(len=*), intent(in) :: job

      integer,          intent(in) :: nymd   ! Date: year-month-day
      integer,          intent(in) :: nhms   ! Time: hour-min-sec

! !OUTPUT PARAMETERS:

      integer, intent(out), optional :: stat

! !DESCRIPTION: Reads in model state vector.
!
! !REVISION HISTORY:
!
!  19Nov2002  Todling   Initial code (as per ana.f).
!  25Sep2003  RG/RT     Bug fix: added dyn_null.
!  29Jan2004  Todling   Fixed exit for "present(stat)" and field not found
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*GetPert1_'

      integer  nymd_in, nhms_in
      integer  i, j, k
      integer  ier
      character(len=255) :: fname

      type (dyn_vect) :: w_f

      if(present(stat)) stat=0

      fname = trim(infile)

!     Make sure w_f does not point to anything
!     ----------------------------------------
      call dyn_null ( w_f )

!     Set date/time of trajectory
!     ---------------------------
      nymd_in = nymd
      nhms_in = nhms

!     Load dynamics vector
!     -------------------
	call timing_on ( 'dyn_get' )
      call dyn_get ( trim(infile), nymd_in, nhms_in, w_f, ier, timidx=0 )
	if ( ier .eq. 0 ) then
            if ( gid==0 ) then
	         print *, myname_//': read perturbation file '//trim(fname)
	         print *, myname_//': nymd, nhms = ', nymd, nhms
                 print *
            endif
        else
          if (present(stat)) then
             stat = 1  ! cannot read state
	     if(gid==0) print *, myname_//': read perturbation failed'
             call dyn_clean ( w_f )
             return
          else
            call die(myname_,'could not read perturbation file '//trim(fname))
          end if
        end if
	call timing_off ( 'dyn_get' )

!     Check dimensions
!     ----------------
      if(im/=w_f%grid%im .or. jm/=w_f%grid%jm .or.
     &   km/=w_f%grid%km .or. lm/=w_f%grid%lm    ) then
         call dyn_clean ( w_f )
         if(present(stat)) then
            stat = 1 
            return
          else
            call die (myname_,'Imcompatible dimension in state vector')
          endif
      endif

      call mp_scatter4d(w_f%u,       u_pt, imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_s, 0)
      call mp_scatter4d(w_f%v,       v_pt, imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_s, ng_d, 0)
      call mp_scatter4d(w_f%pt,     pt_pt, imr,  jnp, nl,  1, jfirst, jlast, 1, nl, ng_d, ng_d, 0)
      call mp_scatter4d(w_f%delp, delp_pt, imr,  jnp, nl,  1, jfirst, jlast, 1, nl,    0,    0, 0)
      call mp_scatter4d(w_f%q,       q_pt, imr,  jnp, nl, nc, jfirst, jlast, 1, nl, ng_d, ng_d, 0)

      call dyn_clean ( w_f )

      end subroutine GetPert1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetPert_:  Read in perturbation state vector
!
! !INTERFACE:

      subroutine GetPert2_ ( infile, nymd, nhms, prog,
     &                       pick, stat, vectype, forceflip )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: infile

      integer,       intent(inout) :: nymd   ! Date: year-month-day
      integer,       intent(inout) :: nhms   ! Time: hour-min-sec

      type(dyn_prog), intent(inout) :: prog

      logical,          intent(in), optional :: pick      ! allows getting last time in file
      integer,          intent(in), optional :: vectype   ! determines whether pert is g4 or g5
      logical,          intent(in), optional :: forceflip ! allows dyn_get to forcefully flip lon's

! !OUTPUT PARAMETERS:

      integer, intent(out), optional :: stat

! !DESCRIPTION: Reads in model state vector.
!
! !REVISION HISTORY:
!
!  19Nov2002  Todling   Initial code (as per ana.f).
!  25Sep2003  RG/RT     Bug fix: added dyn_null.
!  29Jan2004  Todling   Fixed exit for "present(stat)" and field not found
!  26Aug2004  Todling   Added pick option to allow taking last time in file
!  09Jun2006  Elena N.  Bug fixed: Added Dyn_Init call to define vector size
!                                  Replaced Dyn_Clean by Dyn_Null to release pointers
!  19Jun2006  Elena N.  Added tlm Init component
!  22Nov2006  Oloso     Ifdef'd timinig routines for SPMD
!  28Aug2007  Todling   Added vectype and forceflip options
!  20Jan2010  Todling   Only need to read on root and then scatter
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*GetPert2_'


      integer  nymd_in, nhms_in
      integer  i, j, k, n, myID, w_f_grid_lm
      integer  ier
      logical  pick_time
      character(len=255) :: fname

      real(r8), allocatable ::    utmp(:,:,:), vtmp(:,:,:)
      real(r8), allocatable :: delptmp(:,:,:), pttmp(:,:,:), qtmp(:,:,:,:)

      type (dyn_vect) :: w_f

      pick_time=.true.
      if(present(stat)) stat = 0
      if(present(pick)) pick_time = pick

!     Make sure w_f does not point to anything
!     ----------------------------------------
      fname = trim(infile)

      call MP_comm_rank(MP_comm_world,myID,ier)
         if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

!     Get memory for temporary arrays
!     -------------------------------
      if ( myid==ROOT ) then
          allocate (    utmp(imr,jnp,nl),  vtmp(imr,jnp,nl),
     .               delptmp(imr,jnp,nl), pttmp(imr,jnp,nl), qtmp(imr,jnp,nl,nc) )
      else           
          allocate (    utmp(0,0,0), vtmp(0,0,0),
     .               delptmp(0,0,0), pttmp(0,0,0), qtmp(0,0,0,0) )
      endif

      if ( myid== ROOT ) then

      call dyn_null ( w_f )

!     Load dynamics vector
!     -------------------
      call timing_on ( 'dyn_get' )
      if ( pick_time ) then
           nymd_in = nymd
           nhms_in = nhms
           call dyn_get ( trim(infile), nymd_in, nhms_in, w_f, ier, timidx=0, vectype=vectype, forceflip=forceflip )
      else
           call dyn_get ( trim(infile), nymd_in, nhms_in, w_f, ier, vectype=vectype, forceflip=forceflip )
           nymd = nymd_in
           nhms = nhms_in
      endif
	if ( ier .eq. 0 ) then
           if ( gid==0 ) then
	      print *, myname_//': read perturbation file '//trim(fname)
	      print *, myname_//': nymd, nhms = ', nymd, nhms
  	      print *
           endif
        else
          if (present(stat)) then
             stat = 1  ! cannot read state
	     if(gid==0) print *, myname_//': read perturbation failed'
             call dyn_clean ( w_f )
             return
          else
            call die(myname_,'could not read perturbation state file '//trim(fname))
          end if
        end if
	call timing_off ( 'dyn_get' )

              utmp = w_f%u
              vtmp = w_f%v
             pttmp = w_f%pt
           delptmp = w_f%delp
              qtmp = w_f%q

           w_f_grid_lm = w_f%grid%lm

      call dyn_null ( w_f )

      endif ! < myid=ROOT >

      call mp_bcst_int (w_f_grid_lm)

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

      deallocate ( utmp, vtmp, pttmp, delptmp, qtmp )

      end subroutine GetPert2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetPert3_:  Read in tracer perturbations
!
! !INTERFACE:

      subroutine GetPert3_ ( infile, nymd, nhms, qpert, stat, RCfile )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: infile

      integer,          intent(in) :: nymd   ! Date: year-month-day
      integer,          intent(in) :: nhms   ! Time: hour-min-sec
      character(len=*), intent(in), optional :: RCfile

! !OUTPUT PARAMETERS:

      type(dynq)                     :: qpert

      integer, intent(out), optional :: stat

! !DESCRIPTION: Reads in tracer perturbations.
!
! !REVISION HISTORY:
!
!  25Sep2007  Todling  Overload get to handle qpert
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*GetPert3_'

      integer  nymd_in, nhms_in
      integer  i, j, k
      integer  ier
      character(len=255) :: fname
      character(len=255) :: vname (nc)
      character(len=255) :: vunits(nc)
      character(len=255) :: vtitle(nc)
      real(kind=r8), allocatable :: q_x(:,:,:,:)

      if(present(stat)) stat=0

      fname = trim(infile)

!     Set date/time of trajectory
!     ---------------------------
      nymd_in = nymd
      nhms_in = nhms

!     Allocate auxiliar arrays
!     ------------------------
!     if (gid==0) then
          allocate( q_x(imr,jnp,nl,nc) )
!     else
!         allocate( q_x(0,0,0,0) )
!     endif

!     Fill auxiliar informational arrays
!     ----------------------------------
      call readRC_ ( nc, vname, vunits, vtitle, ier, RCfile )
        if ( ier/=0 ) then
           if(present(stat))then
              stat = ier
              return
           else
             call MP_die (myname_,'Error from readRC_',ier)
           endif
        end if

!     Load dynamics vector
!     -------------------
      if ( gid==0 ) then
           call GFIO_GetFld ( fname, nymd, nhms, imr, jnp, nl, nc,
     &                        vname, q_x, verbose=.true., stat=ier )
           if ( ier .eq. 0 ) then
               if ( gid==0 ) then
	            print *, myname_//': read perturbation file '//trim(fname)
	            print *, myname_//': nymd, nhms = ', nymd, nhms
                    print *
               endif
           else
             if (present(stat)) then
                stat = 1  ! cannot read state
	        if(gid==0) print *, myname_//': read perturbation failed'
                return
             else
               call die(myname_,'could not read perturbation file '//trim(fname))
             end if
           end if
      endif ! < gid >

!     Check dimensions
!     ----------------
      call mp_scatter4d(q_x, qpert%q, imr,  jnp, nl, nc, jfirst, jlast, 1, nl, ng_d, ng_d, 0)

!     Release local memory
!     --------------------
      deallocate (q_x)

      end subroutine GetPert3_

      subroutine dims_ ( fname, myimr, myjnp, mynl, mync, stat )

      implicit none

      character(len=*), intent(in) :: fname 
      integer,          intent(out) :: myimr, myjnp, mynl, mync
      integer,          intent(out) :: stat

      integer, parameter :: READ_ONLY = 1
      integer  nvars, ngatts, mylm
      integer  fid,ierr

!     Open the file
!     -------------
      call GFIO_Open ( trim(fname), READ_ONLY, fid, ierr )
      if ( ierr .ne. 0 ) then
           stat = 91
           return
      endif
  
!     Get dimensions
!     --------------
      call GFIO_DimInquire ( fid, myimr, myjnp, mynl, mylm, nvars, ngatts, ierr )
      if ( ierr .ne. 0 ) then
         stat = 92
         call GFIO_close ( fid, ierr )
         return
      end if

!     Now, this is ugly, should be tied to m_dyn
!     ------------------------------------------
      mync = nvars - 9
      if ( mync<1 ) then
         stat = 93
          call GFIO_close ( fid, ierr )
          return
      endif
 
!    Close GFIO file
!     ---------------
      call GFIO_close ( fid, ierr )
      if ( ierr .ne. 0 ) then
           stat = 94
           return
      endif

      end subroutine dims_

      end module m_trajmng
