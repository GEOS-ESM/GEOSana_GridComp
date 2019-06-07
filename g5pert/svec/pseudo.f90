
  program pseudo 

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !PROGRAM: pseudo - Pseudo-inverse fields based on forecast error projection
!
! !USES:
                                                                                                                           
  use m_mpif90,only : MP_init
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world

  use m_pseudo, only : pseudo_set
  use m_pseudo, only : pseudo_process
  use m_pseudo, only : pseudo_clean

  use m_die, only : MP_die, die

  implicit none
                                                                                                                           
! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  05Dec2005  RG/EN   Adopted similarity calculation code to create this program 
!  02Feb2007  Todling    Placed mpi hooks (only runs in 1 pe)
!
!EOP
!-------------------------------------------------------------------------
 
  character(len=*), parameter :: myname= 'pseudoInverse'

  integer, parameter :: ROOT = 0 ! MPI ROOT PE

  character(len=255) :: RCfile
  character(len=255) :: jobR

  integer :: nvecsA
  integer :: nvecsB
  integer :: ntimes               !  number of dates and times examined
  integer :: nymd                 ! date (YYYYMMDD) of stats
  integer :: nhms                 ! time   (HHMMSS) of stats
  character(len=255) :: FSVEC
  character(len=255) :: ISVEC 
  character(len=255), allocatable :: jobB(:)
  character(len=255) :: svalu_file
  character(len=255) :: dynfiles(4)             ! intent to output 4 dynamic vectors:
                                                ! 1. forecast error
                                                ! 2. forecast error projected on nvecsA EVOLVED SVECs
                                                ! 3. pseudo-inverse fields
                                                ! 4. sensitivity projected on nvecsA INITIAL SVECs
  integer :: myID, ier

! MPI Initialization
! ------------------
  call MP_init(ier)
        if(ier/=0) call MP_die(myname,'MP_init()',ier)
  call MP_comm_rank(MP_comm_world,myID,ier)
        if(ier/=0) call MP_die(myname,'MP_comm_rank()',ier)

  nvecsB = 2

! Parse command line
! ------------------
  call Init_ ( nvecsA, dynfiles, RCfile )

! Set up program's options
! ------------------------
  allocate ( jobB(nvecsB) )
  call pseudo_Set ( MP_comm_world, root,    &
                     nvecsA, nvecsB, jobR, FSVEC, ISVEC, svalu_file, jobB,          &
                     RCfile, ier )
     if(ier/=0) call die(myname,'could not setup run')

! Generate initial perturbation
! -----------------------------
  call pseudo_Process ( nvecsA, nvecsB, &
                         jobR, FSVEC, ISVEC, svalu_file, jobB, dynfiles, RCfile)

! Clean up
! --------
  deallocate ( jobB )
  call pseudo_Clean( )

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
      subroutine Init_ ( nvecsA, dynfiles, RCfile )

      implicit NONE

      integer,       intent(out) :: nvecsA
      character*255, intent(out) :: dynfiles(4)           
      character*255, intent(out) :: RCfile     ! resource file
      
!
! !REVISION HISTORY:
!      10Dec2005    EN   Adopted simsv.f90 to create this application 
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
      RCfile   ='pseudo.rc'   ! default not to use an RCfile
      dynfiles(1) = 'FCSTERR.nc4'
      dynfiles(2) = 'PROJFERR.nc4'
      dynfiles(3) = 'PSEUDO.nc4'
      dynfiles(4) = 'PROJSENS.nc4'
      nvecsA = 1

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
           case ("-nsvecs")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iarg, str )
             read(str,*) nvecsA
           case ("-rc")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iarg, RCfile )
           case ("-fcsterr")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, dynfiles(1) )
           case ("-projferr")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, dynfiles(2) )
           case ("-pseudo")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, dynfiles(3) )
           case ("-projsens")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, dynfiles(4) )
           case default
             call die(myname,'not a valid option')
         end select
      end do

      if (nvecsA .lt. 1) then
         print *
         print *, ' Error: No SVECs in a list '
         print *, '     Exiting on error    '
         print *
         stop
      endif

      end subroutine Init_

!.................................................................

      subroutine usage()
      print *
      print *, '  -------------------------------------------------------------------'
      print *, '  pseudo.x - computes full and projected forecast error, pseudo-inverse'
      print *, '                     and projected sensitivity' 
      print *, '  -------------------------------------------------------------------'
      print *
      print *
      print *,'Usage: '
      print *,'  pseudo.x -nsvecs NSVECS -rc RCfile '
      print *,'          [-fcsterr dynfile1] [-projferr dynfile2] [-pseudo dynfile3]'
      print *,'          [-projsens dynfile4] '
      print *
      print *, 'where'
      print *
      print *, '-nsvecs   NSVECS   number of SVECS to be used for projection'
      print *
      print *, '-rc       RCfile   resource filename '
      print *, '              (default: pseudo.rc)'
      print *
      print *, '-fcsterr  dynfile1   Filename to write out forecast error'
      print *, '              (default:  FCSTERR.nc4)'
      print *
      print *, '-projferr dynfile2   Filename to write out projected forecast error'
      print *, '              (default:  PROJFERR.nc4)'
      print *
      print *, '-pseudo   dynfile3   Filename to write out pseudo-inverse fileds'
      print *, '              (default:  PSEUDO.nc4)'
      print *
      print *, '-projsens dynfile4   Filename to write out projected sensitivity'
      print *, '              (default:  PROJSENS.nc4)'
      print *
      print *, ' Last updated: 05 Mar 2009 Todling '
      print *
      call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)
      call exit(1)
      end subroutine usage
      
!.................................................................

 end program pseudo 
