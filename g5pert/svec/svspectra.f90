
 program svspectra 

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !PROGRAM: svspectra - energy spectra calculations for singular vector
!
! !USES:
                                                                                                                           
  use m_mpif90,only : MP_init
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world

  use m_svspectra, only : svspectra_Process
  use m_svspectra, only : svspectra_Clean

  use m_die, only : MP_die

  implicit none
                                                                                                                           
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  25Aug2005   Elena N./Ron Gelaro  Coding based on sv similarity code
!  02Feb2007  Todling    Placed mpi hooks (only runs in 1 pe)
!
!EOP
!-------------------------------------------------------------------------
 
  character(len=*), parameter :: myname= 'svspectra_'


  integer :: nvecs
  integer :: nymd(2)   ! Date of svecs: 1 - ini, 2 - final
  integer :: nhms(2)   ! Time of svecs: 1 - ini, 2 - final
  integer :: etype     ! energy type to be plotted
  integer ::  ifwd, ibak    ! Number of hours of forward/backward integration
  real :: svec_fac      ! Factor to multiply ISVEC's for plotting
  character*255 :: expid     ! experiment id
  character*255 :: svecnormI ! norm used to calculate SVecs ( 'ke', 'te', or else)
  character*255 :: res       ! resolution of the experiment to be shown at plots labels
  logical :: TE_scale
  integer :: myID, ier

! MPI Initialization
! ------------------
  call MP_init(ier)
        if(ier/=0) call MP_die(myname,'MP_init()',ier)
  call MP_comm_rank(MP_comm_world,myID,ier)
        if(ier/=0) call MP_die(myname,'MP_comm_rank()',ier)

!     Defaults
!     --------
      nvecs = 1
      svecnormI =  'te'
      expid = 'EXP_ID'
      res = 'b32'
      etype = 1
      TE_scale = .false.

      nymd(1) = 20050210
      nhms(1) = 00
      nymd(2) = 20050212
      nhms(2) = 00
      svec_fac = 40.0
      ifwd = 48
      ibak = 48


! Parse command line
! ------------------
  call Init_ ( expid, nymd, nhms, svecnormI, nvecs, etype, svec_fac, ifwd, ibak, res, TE_scale )

! Generate initial perturbation
! -----------------------------
  call svspectra_Process ( expid, nymd, nhms, svecnormI, nvecs, etype, svec_fac, ifwd, ibak, res, TE_scale)

! Clean up
! --------
  call svspectra_Clean( )

  call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)

  contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Init_ --- Initialize 
!
! !DESCRIPTION: parses command line.
!
! !INTERFACE:
!
      subroutine Init_ ( expid, nymd, nhms, svecnormI, nvecs, etype, svec_fac, ifwd, ibak, res, TE_scale )

      implicit NONE

      integer,       intent(inout) :: nvecs     ! Total number of vector pairs (ini - fin) to process
      integer,       intent(inout) :: nymd(2)   ! Date of svecs: 1 - ini, 2 - final
      integer,       intent(inout) :: nhms(2)   ! Time of svecs: 1 - ini, 2 - final
      integer,       intent(inout) :: etype     ! energy type to be plotted
      integer,       intent(inout) :: ifwd, ibak    ! Number of hours of forward/backward integration
      real,      intent(inout) :: svec_fac      ! Factor to multiply ISVEC's for plotting 
      character(len=255), intent(inout) :: expid     ! experiment id 
      character(len=255), intent(inout) :: svecnormI ! norm used to calculate SVecs ( 'ke', 'te', or else)
      character(len=255), intent(inout) :: res           !  Resolution of the experiment to be shown at plot labels
      logical, intent(inout) :: TE_scale
!
! !REVISION HISTORY:
!       15Nov2004  Winslow/RT   Initial code, based on similar codes
!
!EOP
!BOC

      character*4, parameter :: myname = 'init_svspectra_'

      integer iret, i, iarg, argc, iargc
      character(len=255) :: argv

      integer mfiles, nfiles

!     Defaults
!     --------
      mfiles = 7               ! number of required input values
                               ! -tbeg 2 -tend 2 -expid 1 -res 1 -nvecs 1  

!     Parse command line
!     ------------------
      argc =  iargc()
      if ( argc .lt. 1 ) call usage()
                                                                                                                                                                                                                                                           
      iarg = 0
      nfiles = 0
                                                                                                                                                                                                                                                           
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) exit
         call GetArg ( iarg, argv )
         if (index(argv,'-expid' ) .gt. 0 ) then
            iarg = iarg + 1
            call GetArg ( iArg, expid )
            nfiles = nfiles + 1
         else if (index(argv,'-E_type' ) .gt. 0 ) then
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*)  etype
         else if (index(argv,'-TE_scale' ) .gt. 0 ) then
            TE_scale = .true.
         else if (index(argv,'-tbeg' ) .gt. 0 ) then
            if ( iarg+2 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) nymd(1)
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) nhms(1)
            nfiles = nfiles + 2
         else if (index(argv,'-tend' ) .gt. 0 ) then
            if ( iarg+2 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) nymd(2)
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) nhms(2)
            nfiles = nfiles + 2
         else if (index(argv,'-svecnorm' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, svecnormI )
         else if (index(argv,'-res' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, res )
            nfiles = nfiles + 1
         else if (index(argv,'-nvecs' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*)  nvecs 
            nfiles = nfiles + 1
         else if (index(argv,'-fac' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) svec_fac 
         else if (index(argv,'-ifwd' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) ifwd 
         else if (index(argv,'-ibak' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) ibak 
         end if
      end do

!      if ( TE_scale ) then
!          svec_fac = 1.0
!          print *, 'Evolved SVECs will be scaled by the total energy norm'
!      endif
                                                                                                                                                                                                                                                           
             if ( nfiles .lt. mfiles) then
                call usage()
                print *, 'More input values are required'
             endif

      return
      end subroutine Init_

!.................................................................

      subroutine usage()
      print *
      print *, '  -------------------------------------------------------------------'
      print *, '  svspectra - Energy spectra calculation for singular vectors                  '
      print *, '  -------------------------------------------------------------------'
! ********************************************************************
! ************************ RUNTIME VARIABLES **************************
!
!  define the variable to be plotted
!
!     ivar= 1 total energy (default)
!           2 kinetic energy
!           3 potential energy (temp only)
!           4 rotational kinetic energy (vor only)
!           5 divergent kinetic energy (div only)
!
!
! number of singular vectors to plot
      print *
      print *
      print *,'Usage: '
      print *,'  svspectra.x -expid EXP_ID -tbeg YYYYMMDD HHMMSS -tend YYYYMMDD HHMMSS'
      print *,'             -nvecs NVECS -res RES [-TE_scale]'
      print *,'             [-E_type etype] [-svecnorm SVECNORM]'
      print *,'             [-fac FACTOR] [-ifwd HOURS] [-ibak HOURS]'
      print *
      print *, 'where'
      print *
      print *, '-expid    EXPERIMENT ID  '
      print *
      print *, '-tbeg    DATE/TIME for initial singular vectors'
      print *
      print *, '-tend    DATE/TIME for evolved singular vectors'
      print *
      print *, '-nvecs     Number of singular vector pairs to process'
      print *
      print *, '-res       Resolution to be used for plotting labels'
      print *
      print *, '-TE_scale  OPTIONAL: scale evolved SVECs by total energy norm'
      print *, '                     before plotting'
      print *
      print *, '-E_type    OPTIONAL: Energy type to be plotted '
      print *, '              (default: 1 = total energy)'
      print *, '              (         2 kinetic energy)'
      print *, '              (         3 potential energy (temp only))'
      print *
      print *, '-fac       OPTIONAL: Factor to multiply ISVECs for plotting'
      print *, '              (default: 40.0)'
      print *
      print *, '-ifwd      OPTIONAL: Number of hours in forward integration'
      print *, '              (default: 48)'
      print *
      print *, '-ibak      OPTIONAL: Number of hours in backward integration'
      print *, '              (default: 48)'
      print *
      print *, '-svecnorm  OPTIONAL: NORM used to calculate'
      print *, '                     singular vectors, defines SVEC filenames'
      print *, '              (default: te ) '
      print *
      print *, ' Last updated: 02 Feb 2007; Todling '
      print *
      call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)
      call exit(1)
      return
      end subroutine usage
      
!---------------------------------------------------------------------------

 end program svspectra 
