!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_fvPSAS --- Implements analysis interface to FVGCM
! 
! !INTERFACE:
!
      MODULE  m_fvpsas
            
! !USES:

      use  m_dyn                ! dynamics state vector type & methods
      use  m_stdio,  only: stdout
      use  m_ioutil, only: luflush

      Implicit NONE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  fvPSAS            ! Performs statiscal analysis and updates
                                !  model state.
      PUBLIC  fvPSAS_out        ! Output GCM state vector
!
! !DESCRIPTION: This module implements the main interface between 
!               the PSAS based statistical analsysis and FVGCM. 
!
! !REVISION HISTORY: 
!
! 23dec1999  da Silva  Initial specs and prologues.
! 03feb2000  da Silva  Added nstep, fvPSAS_out().
! 08feb2000  da Silva  Built actual ana.x interface.
! 10apr2000  da Silva  Replace ana.x interface with perl script fvana 
!                      interface.
! 31Oct2002  Todling   Made fvPSAS_out public
!
!EOP
!-------------------------------------------------------------------------


      CONTAINS


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  fvPSAS --- Statistical Analysis/FVGCM Interface
! 
! !INTERFACE:
!

     subroutine fvPSAS ( expid,   nymd,    nhms,   freq,   nstep,    &
                         im,      jm,      km,                       &
                         jfirst,  jlast,   ng_s,   ng_d,             &
                         ptop,    pint,    ks,     ak,     bk,       &
                         phis,    hs_stdv, Ts,     lwi,    ps,       &
                         delp,    u,       v,      pt,     q,        &
                         delpmin, nremap,  prec,   rc)

!
! !USES:
!

      use  m_stdio,  only: stdout
      use  m_ioutil, only: luflush

#if defined (SPMD)
      use mod_comm
#endif

!
! !INPUT PARAMETERS: 
!

      character(len=*), intent(in)  :: expid   ! Experiment identifier

      integer,          intent(in)   :: nymd   ! Date: year-month-day
      integer,          intent(in)   :: nhms   ! Time: hour-min-sec

      integer,          intent(in)   :: freq   ! Analysis frequency (HHMMSS)
      integer,          intent(in)   :: nstep  ! model time step

                                               ! First guess dimensions
      integer, intent(in)  :: im               !  zonal
      integer, intent(in)  :: jm               !  meridional
      integer, intent(in)  :: km               !  vertical

                                               ! SPMD domain information
      integer, intent(in)  :: jfirst           !  first latitude
      integer, intent(in)  :: jlast            !  last latitude
      integer, intent(in)  :: ng_s             !  ghost zones
      integer, intent(in)  :: ng_d             !  ghost zones

      real,    intent(in)  ::  ptop            ! Pressure at top (Pa)
      real,    intent(in)  ::  pint            ! Pressure at interface (Pa)
      integer, intent(in)  ::  ks              ! no. of pure pressure layers
      real,    intent(in)  ::  ak(km+1)        ! vertical grid a-coefficient
      real,    intent(in)  ::  bk(km+1)        ! vertical grid b-coefficient


      real,    intent(in)  ::  phis(im,jm)     ! Topography geopotential 
                                               ! (meter2/sec2)

      real,    intent(in)  ::  hs_stdv(im,jm)  ! Topography height
                                               !  stdv (meter)

      real,    intent(in)  ::  Ts(im,jm)       ! Sea surface temperature (K)

      real,    intent(in)  ::  lwi(im,jm)      ! Land-water-ice mask:
                                               !   lwi = 0  over ocean
                                               !   lwi = 1  over land
                                               !   lwi = 2  over sea ice
                                               ! NOTE: same as ORO in FVGCM.

      integer, intent(in), OPTIONAL ::  prec   ! precision for GFIO:
                                                   ! 0 = 32 bits
                                                   ! 1 = 64 bits
!
! !INPUT/OUTPUT PARAMETERS:
!
                                               ! First Guess/Analysis fields:
      real,    intent(inout)  ::  ps(im,jfirst:jlast)      
                                               ! Surface pressure (Pa)
      real,    intent(inout)  ::  delp(im,jfirst:jlast,km) 
                                               ! Delta pressure (Pa)
      real,    intent(inout)  ::  u(im,jfirst-ng_d:jlast+ng_s,km) 
                                               ! u-wind (m/s)
      real,    intent(inout)  ::  v(im,jfirst-ng_s:jlast+ng_d,km) 
                                               ! v-wind (m/s)
      real,    intent(inout)  ::  pt(im,jfirst-ng_d:jlast+ng_d,km) 
                                               ! scaled virtual potential
                                               !  temperature (T p**kappa)
      real,    intent(inout)  ::  q(im,jfirst-ng_d:jlast+ng_d,km,1)
                                               ! Specific humdity (kg/kg)
!
! !OUTPUT PARAMETERS:
!
      real,    intent(out)  :: delpmin         ! minimum depth of the lowest layer
      integer, intent(out)  :: nremap          ! number of points need to remap
      integer, intent(out)  :: rc              ! Error return code:
                                                 !  0   all is well
                                                 !  1   ...

! !DESCRIPTION: Given the first guess fields ({\tt delp\_f, u\_f, v\_f, pt\_f,}
!               and {\tt q\_f}) from FVGCM this routine invokes the PSAS
!  based statistical analysis to compute the analysis fields {\tt w\_a});
!  the analysis fields are returned in the same arraysbas the first guess.
!  This version implements a {\em Nested Executable} interface. 
!
! !REVISION HISTORY: 
!
! 30Nov2001  Chern     First crack
! 11Dec2001  da Silva  Moved from psasdrv.F here, renamed
!
!EOP
!-------------------------------------------------------------------------

#if defined (SPMD)

#define  CPP_PS   pstmp
#define  CPP_U    utmp
#define  CPP_V    vtmp
#define  CPP_DELP delptmp
#define  CPP_PT   pttmp
#define  CPP_Q    qtmp
#else
#define CPP_PS   ps
#define CPP_U    u
#define CPP_V    v
#define CPP_DELP delp
#define CPP_PT   pt
#define CPP_Q    q
   integer :: gid = 0
#endif

! Local Variables
      integer :: i, j                 ! index
#if defined (SPMD)
! Gathered first guess fields
      real    :: pstmp( im,jm ),    utmp( im,jm,km ) 
      real    :: vtmp( im,jm,km ),  delptmp( im,jm,km )
      real    :: pttmp( im,jm,km ), qtmp(im,jm,km,1)
#endif

#if defined (SPMD)           
      call mp_gather4d(ps,   pstmp,   im,   jm,   1,  1, jfirst, jlast, &
                       1,    1,       0,    0,    0)
      call mp_gather4d(u,    utmp,    im,   jm,   km, 1, jfirst, jlast, &
                       1,    km,      ng_d, ng_s, 0)
      call mp_gather4d(v,    vtmp,    im,   jm,   km, 1, jfirst, jlast, &
                       1,    km,      ng_s, ng_d, 0)
      call mp_gather4d(pt,   pttmp,   im,   jm,   km, 1, jfirst, jlast, &
                       1,    km,      ng_d, ng_d, 0)
      call mp_gather4d(delp, delptmp, im,   jm,   km, 1, jfirst, jlast, &
                       1,    km,      0,    0,    0)
      call mp_gather4d(q,    qtmp,    im,   jm,   km, 1, jfirst, jlast, &
                       1,    km,      ng_d, ng_d, 0)
#endif

      if (gid == 0 ) then

        call  fvPSAS_ ( expid,  nymd,     nhms,  freq,  nstep,        &
                        im,     jm,       km,                         &
                        ptop,   pint,     ks,    ak,    bk,           &
                        phis,   hs_stdv,  Ts,    lwi,                 &
                        CPP_PS, CPP_DELP, CPP_U, CPP_V, CPP_PT,       &
                        CPP_Q,  rc,       prec )

        if ( rc .ne. 0 ) then
          print *, 'fvpsas_: error from fvPSAS_(), rc ', rc
          call exit(1)
        else
          write(6,*) 'fvPSAS_ successfully called at ',nymd,'/',nhms
        endif
        call luFlush ( stdout )

! Check whether remapping is necessary

        delpmin = 100000.0
        nremap  = 0
        do j = 1, jm
           do i = 1, im
              ! check if the lowest layer thickness < 0.1 mb
              if( CPP_DELP (i,j,km) < 10.0 ) then
                 nremap  = nremap + 1
                 delpmin = min( delpmin, CPP_DELP (i,j,km) )
              endif
           enddo
        enddo

      endif

#if defined (SPMD)           
      call mp_barrier()

! Scatter (broadcast)
      call mp_bcst_real(delpmin)
      call mp_bcst_int(nremap)
      call mp_scatter2d(pstmp,   ps,   im,   jm,   jfirst, jlast, 0) 
      call mp_scatter4d(utmp,    u,    im,   jm,   km, 1, jfirst, jlast, &
                        1,       km,   ng_d, ng_s, 0)
      call mp_scatter4d(vtmp,    v,    im,   jm,   km, 1, jfirst, jlast, &
                        1,       km,   ng_s, ng_d, 0)
      call mp_scatter4d(pttmp,   pt,   im,   jm,   km, 1, jfirst, jlast, &
                        1,       km,   ng_d, ng_d, 0)
      call mp_scatter4d(delptmp, delp, im,   jm,   km, 1, jfirst, jlast, &
                        1,       km,   0,    0,    0)
      call mp_scatter4d(qtmp,    q,    im,   jm,   km, 1, jfirst, jlast, &
                        1,       km,   ng_d, ng_d, 0)
#endif          

     return
     end subroutine fvpsas

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
!
! !IROUTINE:  fvPSAS_ --- Statistical Analysis/FVGCM INTERNAL Interface
! 
! !INTERFACE:
!
   subroutine fvPSAS_ ( expid, nymd,    nhms,  freq,  nstep,         &
                        im,    jm,      km,                          &
                        ptop,  pint,    ks,    ak,    bk,            &
                        phis,  hs_stdv, Ts,    lwi,                  &
                        ps,    delp,    u,     v,     pt,            &
                        q,     rc,      precision ) 

!
! !USES:
!
   Implicit NONE

!
! !INPUT PARAMETERS: 
!

      character(len=*), intent(in)  :: expid   ! Experiment identifier

      integer,          intent(in)   :: nymd   ! Date: year-month-day
      integer,          intent(in)   :: nhms   ! Time: hour-min-sec

      integer,          intent(in)   :: freq   ! Analysis frequency (HHMMSS)
      integer,          intent(in)   :: nstep  ! model time step

                                               ! First guess dimensions
      integer, intent(in)  :: im               !  zonal
      integer, intent(in)  :: jm               !  meridional
      integer, intent(in)  :: km               !  vertical


      real,    intent(in)  ::  ptop            ! Pressure at top (Pa)
      real,    intent(in)  ::  pint            ! Pressure at interface (Pa)
      integer, intent(in)  ::  ks              ! no. of pure pressure layers
      real,    intent(in)  ::  ak(km+1)        ! vertical grid a-coefficient
      real,    intent(in)  ::  bk(km+1)        ! vertical grid b-coefficient


      real,    intent(in)  ::  phis(im,jm)     ! Topography geopotential 
                                               ! (meter2/sec2)

      real,    intent(in)  ::  hs_stdv(im,jm)  ! Topography height
                                               !  stdv (meter)

      real,    intent(in)  ::  Ts(im,jm)       ! Sea surface temperature (K)

      real,    intent(in)  ::  lwi(im,jm)      ! Land-water-ice mask:
                                               !   lwi = 0  over ocean
                                               !   lwi = 1  over land
                                               !   lwi = 2  over sea ice
                                               ! NOTE: same as ORO in FVGCM.

      integer, intent(in), OPTIONAL ::  precision  ! precision:
                                                   ! 0 = 32 bits
                                                   ! 1 = 64 bits

!
! !INPUT/OUTPUT PARAMETERS:
!
                                                  ! First Guess/Analysis fields:
                                                  ! [0, 360]
      real,    intent(inout)  ::    ps(im,jm)     ! Surface pressure (Pa)
      real,    intent(inout)  ::  delp(im,jm,km)  ! Delta pressure (Pa)
      real,    intent(inout)  ::     u(im,jm,km)  ! u-wind (m/s)
      real,    intent(inout)  ::     v(im,jm,km)  ! v-wind (m/s)
      real,    intent(inout)  ::    pt(im,jm,km)  ! scaled virtual potential 
                                                  !  temperature (T p**kappa)
      real,    intent(inout)  ::     q(im,jm,km,1)! Specific humdity (kg/kg)


      integer, intent(inout)  :: rc               ! Error return code:
                                                  !  0   all is well
                                                  !  1   ...
                                                  !  

! !DESCRIPTION: Given the first guess fields ({\tt delp\_f, u\_f, v\_f, pt\_f,}
!               and {\tt q\_f}) from FVGCM this routine invokes the PSAS
!  based statistical analysis to compute the analysis fields {\tt w\_a});
!  the analysis fields are returned in the same arraysbas the first guess.
!  This version implements a {\em Nested Executable} interface. 
!
! !REVISION HISTORY: 
!
! 23dec1999  da Silva  Initial code.
! 10apr2000  da Silva  Interface to fvana instead of ana.x.  
! 15mar2001  da Silva  Made fvCCM arrays in/out
! ??Apr2001  Chern     Converted (ps_a,ps_f) -> ps, etc., changing 
!                      argument list; set precision to 32 bits
! 04May2001  da Silva  Set default precision back to 64 bits.
! 11May2001  da Silva  Added 'precision' to argument list (optional)
! 11Jan2002  da Silva  Added replay mode.
! 31Jan2002  da Silva  Added myenv handling.
! 06Mar2009  Todling   Update to dyn_init API; put anything in fractions
!
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'fvPSAS'

      integer, parameter :: lm = 1    ! tracer dimension: sphu only, no tracers

      character(len=255) :: fname_f   ! first guess file name 
      character(len=255) :: fname_a   ! analysis file name
      character(len=255) :: fname_o   ! post-analysis ODS file name
      character(len=255) :: cmd       ! command to spawn analysis process
   
  
      integer            :: prec = 1  !  Precision of output files:
                                      !      0 = 32 bits
                                      !   else = 64 bits

      character(len=2) chour      
      character(len=8) cymd
      character(len=6) chms
      character(len=1) cprec

      type(dyn_vect) w_f   ! first guess
      type(dyn_vect) w_a   ! analysis

      integer ier, nymd_a, nhms_a, system
      logical replay, myenv


      write(*,'(a,2i10)') myname // ': starting fvPSAS on ', nymd, nhms


!     Set precision
!     -------------
      if ( present(precision) ) prec = precision
      if ( prec .eq. 0 ) then
           cprec = '0'
      else
           prec  = 1      
           cprec = '1'
      end if

!     Make sure w_f/w_a do not point to anything
!     ------------------------------------------
      call dyn_null ( w_f )
      call dyn_null ( w_a )

!     Assign pointers in w_f input arrays
!     -----------------------------------
      call Dyn_Init ( im, jm, km, lm, ptop, ks, ak, bk,      &
                      phis, hs_stdv, Ts, lwi, ps,            &
                      lwi, lwi, lwi, lwi, lwi,               &
                      delp, u, v, pt, q,                     &
                      w_f, ier  )
      if ( ier .ne. 0 ) then
         rc = 1
         return
      end if


!     Assign pointers in w_a to output arrays
!     ---------------------------------------
      call Dyn_Init ( im, jm, km, lm, ptop, ks, ak, bk,      &
                      phis, hs_stdv, Ts, lwi, ps,            &
                      lwi, lwi, lwi, lwi, lwi,               &
                      delp, u, v, pt, q,                     &
                      w_a, ier  )
      if ( ier .ne. 0 ) then
         rc = 2
         return
      end if


!     Form file names
!     ---------------
      if ( nhms .gt. 240000 .or. nymd .gt. 29999999 ) then
         rc = 3
         return
      end if
      write(chour,'(I2.2)') nhms / 10000
      write(cymd,'(I8.8)')  nymd
      write(chms,'(i6.6)')  nhms
      fname_f = trim(expid) // '.bkg.eta.' // cymd // '.hdf' 
      fname_a = trim(expid) // '.ana.eta.' // cymd // '.hdf' 
      fname_o = trim(expid) // '.ana.obs.' // cymd // '.ods' 


!     Either replay mode...
!     ---------------------
      inquire ( file=trim('replay.acq'), exist=replay )
      if ( replay ) then

        write(*,'(a)') myname//': replaying analysis from GFIO file: ' &
                      //trim(fname_a)
        call luFlush ( stdout )
        nymd_a = nymd
        nhms_a = nhms
        call  dyn_get ( fname_a, nymd_a, nhms_a, w_a, ier, timidx=0 )
        if ( ier .ne. 0 ) then
           rc = 6
           return
        end if

!     ... or regular analysis cycle
!     -----------------------------
      else

!       Save first guess to disk
!       ------------------------
        write(*,'(a)') myname//': saving GFIO file: '//trim(fname_f)
        call dyn_put ( fname_f, nymd, nhms, prec, w_f, ier, &
                       nstep=nstep, freq=freq )
        if ( ier .ne. 0 ) then
           print *, 'ier = ', ier
           rc = 4
           return
        end if

!       ---------------------------------------------------------------
!       NOTE: If file .myenv exixts, then we the take environment from 
!       this file. This hack is needed because some early mpirun 
!       implementation did not export the environment to its children.
!       In this case, we save the environment in the parent script
!       with "printenv > .myenv"
!       ---------------------------------------------------------------

!       Spawn analysis process
!       ----------------------
        cmd = 'fvana -o ' // trim(fname_o)  &
            //     ' -a ' // trim(fname_a)  &
            //     ' -prec ' // cprec       &
            // ' ' // cymd // ' ' // chms // ' ' // trim(fname_f) 
        write(*,'(a)') myname//': starting analysis process: '//trim(cmd)
        call luFlush ( stdout )
        inquire ( file=trim('.myenv'), exist=myenv ) ! see NOTE above
        if ( myenv ) then
	     ier =  system ( 'env "`cat .myenv`" ' // trim(cmd) ) 
        else
	     ier =  system ( trim(cmd) )
        end if
        call luFlush ( stdout )
        if ( ier .ne. 0 ) then
           rc = 5
           return
        end if


!       Retrieve analysis fields
!       ------------------------
        write(*,'(a)') myname//': retrieving GFIO file: '//trim(fname_a)
        call luFlush ( stdout )
        call  dyn_get ( fname_a, nymd_a, nhms_a, w_a, ier )
        if ( ier .ne. 0 ) then
           rc = 6
           return
        end if
        if ( nymd_a .ne. nymd .or. nhms_a .ne. nhms ) then
           rc = 7
           return
        end if

     end if


!     Nullify pointers
!     ----------------
      call dyn_null ( w_f )
      call dyn_null ( w_a )


!     All done
!     --------
      rc = 0
      return

      end subroutine fvPSAS_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  fvPSAS_out - Writes out dynamics state vector
! 
! !INTERFACE:
!
   subroutine fvPSAS_out ( expid, nymd, nhms, nstep, freq,          &
                           im, jm, km, ptop, pint, ks, ak, bk,      &
                           phis, hs_stdv, Ts, lwi,                  &
                           ps_f, delp_f, u_f, v_f, pt_f, q_f, rc,   &
                           usrfname )

!
! !USES:
!
   Implicit NONE

!
! !INPUT PARAMETERS: 
!

      character(len=*), intent(in)   :: expid  ! Experiment identifier

      integer,          intent(in)   :: nymd   ! Date: year-month-day
      integer,          intent(in)   :: nhms   ! Time: hour-min-sec

      integer,          intent(in)   :: nstep  ! model time step
      integer,          intent(in)   :: freq   ! frequency (HHMMSS) of
                                               !  output fields

                                               ! First guess dimensions
      integer, intent(in)  :: im               !  zonal
      integer, intent(in)  :: jm               !  meridional
      integer, intent(in)  :: km               !  vertical


      real,    intent(in)  ::  ptop            ! Pressure at top (Pa)
      real,    intent(in)  ::  pint            ! Pressure at interface (Pa)
      integer, intent(in)  ::  ks              ! no. of pure pressure layers
      real,    intent(in)  ::  ak(km+1)        ! vertical grid a-coefficient
      real,    intent(in)  ::  bk(km+1)        ! vertical grid b-coefficient


      real,    intent(in)  ::  phis(im,jm)     ! Topography geopotential 
                                               ! (meter2/sec2)

      real,    intent(in)  ::  hs_stdv(im,jm)  ! Topography height
                                               !  stdv (meter)

      real,    intent(in)  ::  Ts(im,jm)       ! Sea surface temperature (K)

      real,    intent(in)  ::  lwi(im,jm)      ! Land-water-ice mask:
                                               !   lwi = 0  over ocean
                                               !   lwi = 1  over land
                                               !   lwi = 2  over sea ice
                                               ! NOTE: same as ORO in FVGCM.

                                               ! Dynamics fields:
                                               ! [0, 360]
      real,    intent(in)  ::    ps_f(im,jm)     ! Surface pressure (Pa)
      real,    intent(in)  ::  delp_f(im,jm,km)  ! Delta pressure (Pa)
      real,    intent(in)  ::     u_f(im,jm,km)  ! u-wind (m/s)
      real,    intent(in)  ::     v_f(im,jm,km)  ! v-wind (m/s)
      real,    intent(in)  ::    pt_f(im,jm,km)  ! scaled virtual potential 
                                                 !  temperature (T p**kappa)
      real,    intent(in)  ::     q_f(im,jm,km,1)! Specific humdity (kg/kg)

      character(len=*), intent(in), optional :: usrfname ! Template for file name

!
! !OUTPUT PARAMETERS:
!

      integer, intent(out)  :: rc               ! Error return code:
                                                !  0   all is well
                                                !  1   ...
                                                !  

! !DESCRIPTION: This routine writes out the dynamics vector fields 
!               ({\tt delp\_f, u\_f, v\_f, pt\_f,} and {\tt q\_f}) 
!  using GFIO. This routine is nothing more than a FORTRAN 77 style
!  interface to dyn_put(). The output file is in single precision
!  (32 bits), HDF/COARDS compliant.
!
! !REVISION HISTORY: 
!
!  03feb2000 da Silva  Initial code.
!  31oct2002 Todling   Added usrfname.
!  06Mar2009 Todling   Update to dyn_init API; put anything in fractions
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'fvPSAS_out'

      integer, parameter :: lm = 1    ! tracer dimension: sphu only, no tracers

      character(len=255) :: fname_f   ! first guess file 
  
      integer, parameter :: prec = 0  !  Precision of output files:
                                      !      0 = 32 bits
                                      !      1 = 64 bits

      character(len=2) chour      
      character(len=8) cymd

      type(dyn_vect) w_f   ! state vector

      integer ier, nymd_a, nhms_a, ios
      real, allocatable :: epv(:,:,:)


!     Make sure w_f does not point to anything
!     ----------------------------------------
      call dyn_null ( w_f )

!     Assign pointers in w_f input arrays
!     -----------------------------------
      call Dyn_Init   ( im, jm, km, lm, ptop, ks, ak, bk,      &
                        phis, hs_stdv, Ts, lwi, ps_f,        &
                        lwi, lwi, lwi, lwi, lwi,               &
                        delp_f, u_f, v_f, pt_f, q_f, w_f, ier  )
      if ( ier .ne. 0 ) then
         rc = 1
         return
      end if


!     Form file name
!     --------------
      if ( nhms .gt. 240000 .or. nymd .gt. 29999999 ) then
         rc = 3
         return
      end if
      write(chour,'(I2.2)') nhms / 10000
      write(cymd,'(I8.8)')  nymd
      if (present(usrfname)) then
         fname_f = trim(usrfname) // '.prog.eta.' // cymd // '.hdf' 
      else
         fname_f = trim(expid) // '.prog.eta.' // cymd // '.hdf' 
      end if


!     Compute Ertel's PV
!     ------------------
      allocate ( epv(im,jm,km), stat=ios )
      if ( ios .eq. 0 ) then
         call epvd (im,jm,km,1,jm,u_f,v_f,pt_f,delp_f,epv)
      else
         rc = 4
         return
      end if


!     Write the GFIO file
!     -------------------
      call dyn_put ( fname_f, nymd, nhms, prec, w_f, ier, &
                     nstep=nstep, verbose=.true., freq=freq, epv=epv )
      deallocate ( epv )
      if ( ier .ne. 0 ) then
         print *, myname//': ier = ', ier
         rc = 5
         return
      end if


!     Nullify pointers
!     ----------------
      call dyn_null ( w_f )


!     All done
!     --------
      rc = 0
      return

      end subroutine fvPSAS_out

      end module m_fvPSAS
