
  program postStats

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !PROGRAM: postStats - Statistical post processing
!
! !USES:

  use m_mpif90,only : MP_init
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world
                                                                                                                           
  use m_postStats, only : postStats_set
  use m_postStats, only : postStats_process
  use m_postStats, only : postStats_clean

  use m_pertutil, only : pertutil_setparam

  use m_die, only : MP_die, die

  implicit none
                                                                                                                           
! !DESCRIPTION: This routine post processes output from various combinations 
!               of TLM, NLM, and/or ADM runs
!
! To do:
!  \begin{itemize}
!      \item MPI part
!      \item double check this is really correct
!  \end{itemize}
!
! !REVISION HISTORY:
!
!  18Oct2004  Winslow     initial code based on similar code
!  23Sep2005  RG/Elena N. Modified usage/comments for clarity 
!  02Feb2007  Todling    Placed mpi hooks (only runs in 1 pe)
!  01Nov2007  Todling    Removed inconsistency in command line
!
!EOP
!-------------------------------------------------------------------------
 
  character(len=*), parameter :: myname= 'postStats'

  integer, parameter :: mfiles = 20       ! max no. of input dyn-vector files
  integer, parameter :: COMM = 0          ! Fake mpi
  integer, parameter :: ROOT = 0          ! MPI ROOT PE

  character(len=255) :: dynfiles(mfiles)   ! input dyn filename(s)
  character(len=255) :: RCfile
  logical :: pert_with_LPO

  integer :: nfiles               ! actual no. of dyn files read in
  integer :: ntimes               !  number of dates and times examined
  integer, allocatable :: nymd(:) ! date (YYYYMMDD) of stats
  integer, allocatable :: nhms(:) ! time   (HHMMSS) of stats
  integer :: corr_kind(2)         ! slots to be correlated
  logical :: rescale_adj         ! true if adm field comparison
                                  ! should extend through the column
  logical :: e_calc   ! true if energy norm to be calc.
  integer :: myID, ier

! MPI Initialization
! ------------------
  call MP_init(ier)
        if(ier/=0) call MP_die(myname,'MP_init()',ier)
  call MP_comm_rank(MP_comm_world,myID,ier)
        if(ier/=0) call MP_die(myname,'MP_comm_rank()',ier)

! Parse command line
! ------------------
  call init_ ( mfiles, nfiles, dynfiles, RCfile, ntimes, pert_with_LPO )

! Set up program's options
! ------------------------
  allocate ( nymd(ntimes) ); nymd = 0
  allocate ( nhms(ntimes) ); nhms = 0
  call postStats_set ( MP_comm_world, root, &
                       ntimes, nymd, nhms, corr_kind, rescale_adj, e_calc, &
                       RCfile, ier )
     if(ier/=0) call die(myname,'could not setup run')

! Generate initial perturbation
! -----------------------------
  call postStats_Process ( comm, root, nfiles, corr_kind, rescale_adj, e_calc,      &
                           dynfiles, &
                           ntimes, nymd, nhms, RCfile, pert_with_LPO )

! Clean up
! --------
  call postStats_Clean( )
  deallocate ( nhms )
  deallocate ( nymd )

  call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Init_ --- Initialize post processor
!
! !DESCRIPTION: parses command line.
!
! !INTERFACE:
!
      subroutine Init_ ( mfiles, nfiles, dynfiles, RCfile, ntimes, pert_with_LPO )

      implicit NONE

      integer,       intent(in)  :: mfiles     ! max. number of eta files
                                               ! dynamics file names (eta)
      character(len=*), intent(out) :: dynfiles(mfiles) 
      integer,          intent(out) :: nfiles     ! actual no. of input eta files
      character(len=*), intent(out) :: RCfile     ! resource file
      integer,          intent(out) :: ntimes     ! actual no. of output eta files
      logical,          intent(out) :: pert_with_LPO  !defines whether to output perturbation with imposed LPO
      
!
! !REVISION HISTORY:
!       15Oct2004  Winslow   Initial code, based on similar codes
!
!EOP
!BOC

      character(len=*), parameter :: myname = 'init'

      integer iret, i, iarg, argc, iargc
      integer vectype
      character(len=255) :: argv, res
      character*10 str

!     Defaults
!     --------
      RCfile   ='postStats.rc'   ! default not to use an RCfile
      pert_with_LPO = .false.

!     Set defaults
!     ------------
      nfiles = 0
      ntimes = 1
      vectype = 4

!     Parse command line
!     ------------------
      argc =  iargc()
      if ( argc .lt. 1 ) call usage()
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) exit
         call GetArg ( iarg, argv )
         select case (argv)
           case ("-g5")
             vectype = 5
           case ("-ntimes")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iarg, str )
             read(str,*) ntimes
           case ("-rc")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iarg, RCfile )
           case ("-pert_with_LPO")
             pert_with_LPO = .true.
           case default
             nfiles = nfiles + 1
             if ( nfiles .gt. mfiles ) call die(myname,'too many eta files')
             dynfiles(nfiles) = trim(argv)
         end select
      end do
      if ( nfiles .lt. 2 ) call usage()

                       call pertutil_setparam ( 'vectype', vectype )
      if( vectype==5 ) call pertutil_setparam ( 'tname', 'vT_' )

      end subroutine Init_

!.................................................................

      subroutine usage()
      print *
      print *, '  -------------------------------------------------------------------'
      print *, '  stats5 - post processing of output                                 '
      print *, '  -------------------------------------------------------------------'
      print *
      print *
      print *,'Usage: '
      print *,'  postStats.x [-nfiles NFILES] [-ntimes NTIMES] [-rc RCfile] dynfiles '
      print *
      print *, 'where'
      print *
      print *, '-ntimes   NTIMES   number of times to be examined '
      print *
      print *, '-nfiles   NFILES   number of input dynamic files '
      print *, '                   (see postStats.rc to follow an order of files) '
      print *
      print *, '-rc       RCfile   resource filename '
      print *, '                   (default: postStats.rc)'
      print *
      print *, '-pert_with_LPO     output of perturbation with imposed LPO'
      print *, '                   as defined in resource file RCfile'
      print *, '                   (default: no output of perturbation with LPO)'
      print *
      print *, ' Last updated: 01 Nov 2007; Todling '
      print *
      call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)
      call exit(1)
      end subroutine usage
      
!.................................................................

 end program postStats
