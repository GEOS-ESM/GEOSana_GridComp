
  program simSV

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !PROGRAM: simSV - similarity calculations for simgular vector
!
! !USES:

  use m_mpif90,only : MP_init
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world
                                                                                                                           
  use m_simsv, only : simsv_set
  use m_simsv, only : simsv_process
  use m_simsv, only : simsv_clean

  use m_die, only : MP_die, die

  implicit none
                                                                                                                           
! !DESCRIPTION: This routine makes the actual perturbation vector
!
! To do:
!  \begin{itemize}
!      \item MPI part
!      \item double check this is really correct
!  \end{itemize}
!
! !REVISION HISTORY:
!
!  18Nov2004  Winslow    initial code based on similar code
!  02Feb2007  Todling    Placed mpi hooks (only runs in 1 pe)
!
!EOP
!-------------------------------------------------------------------------
 
  character(len=*), parameter :: myname= 'simSV'

  integer, parameter :: ROOT = 0 ! MPI ROOT PE

  character(len=255) :: RCfile
  character(len=255) :: jobR

  integer :: nvecsA
  integer :: nvecsB
  integer :: ntimes               !  number of dates and times examined
  integer :: nymd ! date (YYYYMMDD) of stats
  integer, allocatable :: nhms(:) ! time   (HHMMSS) of stats
  character(len=255), allocatable :: jobA(:)
  character(len=255), allocatable :: jobB(:)
  logical :: details
  integer :: myID, ier

! MPI Initialization
! ------------------
  call MP_init(ier)
        if(ier/=0) call MP_die(myname,'MP_init()',ier)
  call MP_comm_rank(MP_comm_world,myID,ier)
        if(ier/=0) call MP_die(myname,'MP_comm_rank()',ier)

! Parse command line
! ------------------
  call Init_ ( nvecsA, nvecsB, ntimes, RCfile )


! Set up program's options
! ------------------------
  allocate ( nhms(ntimes) )
  allocate ( jobA(nvecsA) )
  allocate ( jobB(nvecsB) )
  call simsv_Set ( MP_comm_world, root,  ntimes, nhms, nymd, details,  &
                     nvecsA, nvecsB, jobR, jobA, jobB,          &
                     RCfile, ier )
     if(ier/=0) call die(myname,'could not setup run')

! Generate initial perturbation
! -----------------------------
  call simsv_Process ( MP_comm_world, ROOT, ntimes, nhms, nymd, details, nvecsA, nvecsB, &
                         jobR, jobA, jobB, RCfile)

! Clean up
! --------
  deallocate ( nhms )
  deallocate ( jobA )
  deallocate ( jobB )
  call simsv_Clean( )

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
      subroutine Init_ ( nvecsA, nvecsB, ntimes, RCfile )

      implicit NONE

      integer,       intent(out) :: nvecsA
      integer,       intent(out) :: nvecsB
      integer,       intent(out) :: ntimes     ! actual no. of output eta files
      character*255, intent(out) :: RCfile     ! resource file
      
!
! !REVISION HISTORY:
!       15Nov2004  Winslow/RT   Initial code, based on similar codes
!
!EOP
!BOC

      character*4, parameter :: myname = 'init'

      integer iret, i, iarg, argc, iargc
      logical verb
      character(len=255) :: argv, res
      character*10 str

!     Defaults
!     --------
      RCfile   ='fvsimsv.rc'   ! default not to use an RCfile

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
           case ("-A")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iarg, str )
             read(str,*) nvecsA
           case ("-B")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iarg, str )
             read(str,*) nvecsB
           case ("-t")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iarg, str )
             read(str,*) ntimes
           case ("-rc")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iarg, RCfile )
           case default
             call die(myname,'not a valid option')
         end select
      end do

      end subroutine Init_

!.................................................................

      subroutine usage()
      print *
      print *, '  -------------------------------------------------------------------'
      print *, '  simSV - similarity transform for singular vectors                  '
      print *, '  -------------------------------------------------------------------'
      print *
      print *
      print *,'Usage: '
      print *,'  simSV.x [-A nvecsA] [-B nvecsB] [-t ntimes] [-rc RCfile] '
      print *
      print *, 'where'
      print *
      print *, '-A   ntimes   no. vectors in set A '
      print *
      print *, '-B   ntimes   no. vectors in set B '
      print *
      print *, '-t   ntimes   number of times to be examined '
      print *
      print *, '-rc   RCfile   resource filename '
      print *, '              (default: fvsimsv.rc)'
      print *
      print *, ' Last updated: 02 Feb 2007; Todling '
      print *
      call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)
      call exit(1)
      end subroutine usage
      
!.................................................................

 end program simSV
