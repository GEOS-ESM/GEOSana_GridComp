!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOI
!
! !TITLE: The GEOS-4/5 Adjoint-based Library
!
! !AUTHORS: 
!
! !AFFILIATION: Global Modeling and Assimilation Office, NASA/GSFC, Greenbelt, MD 20771
!
! !DATE: 03 August 2006
!
! !INTRODUCTION: Adjoint Initialization Program
!
!  This program prepares initial conditions for model adjoint integrations.
!
!EOI
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: initadm - Initializes perturbations for ADM/SENS runs
!
! !INTERFACE:

  program initadm

! !USES:
                                                                                                                           
  use m_mpif90,only : MP_init
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world

  use m_initadj, only : initadj_set
  use m_initadj, only : initadj_iadm
  use m_initadj, only : initadj_clean

  use m_ioutil, only : luavail
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
!  20Nov2002  Todling    Initial code: projection matrix assumed diagonal.
!  20Dec2002  Gelaro     Properly defined projection operator.
!  02Feb2007  Todling    Placed mpi hooks (only runs in 1 pe)
!  07May2007  Todling    Add option -pick
!  27Oct2007  Todling    Support to GEOS-5-type input/output vectors
!  13Dec2007  Todling    Calculate various gradients at "once"
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname= 'InitADM'

  integer, parameter :: mfiles = 2        ! max no. of dyn-vector files
  integer, parameter :: ROOT = 0          ! MPI ROOT PE

  logical :: diff_file                    ! true if difference of two files
  character(len=255) :: dynfiles(mfiles)  ! input dyn filename(s)
  character(len=255) :: jtype(10), init_field
  character(len=255) :: ofname, RCfile
  character(len=255), allocatable :: knames(:)

  integer :: nfiles  ! actual no. of dyn files read in
  integer :: knmax   ! total no. of dyn-like arrays
  integer :: nymd, nhms
  integer :: vectype
  integer :: j,lu,ier
  integer :: idims(2)

! MPI Initialization
! ------------------
  call MP_init(ier)
        if(ier/=0) call MP_die(myname,'MP_init()',ier)

! Parse command line
! ------------------
  call init_ ( mfiles, nfiles, dynfiles, RCfile, nymd, nhms, idims, vectype )

! Set up program's options
! ------------------------
  knmax = nfiles + 1
  allocate ( knames(knmax) )
  call initadj_set ( MP_comm_world, root, diff_file, jtype, init_field, knmax, knames, RCfile, vectype, ier )
     if(ier/=0) call die(myname,'could not setup run')

  lu = luavail()
  open (lu, file='Jnormf.txt', form='formatted')

! Generate initial perturbation
! -----------------------------
  ofname ='Jgradf_XXX.eta.nc4'  ! wired in for now
  j  = 1
  do while ( trim(jtype(j)) .ne. 'NONE' )
     write(ofname(8:10),'(a3)') trim(jtype(j))
     call initadj_iadm ( MP_comm_world, root, lu, trim(ofname), jtype(j), diff_file, init_field, nfiles, dynfiles, & 
                         knames, nymd, nhms, idims, vectype, RCfile )
     j=j+1
  enddo

! Clean up
! --------
  close (lu)
  deallocate ( knames )
  call initadj_clean( )

  call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_ --- Initialize initadj
!
! !INTERFACE:

      subroutine Init_ ( mfiles, nfiles, dynfiles, RCfile, nymd, nhms, &
                         idims, vectype )

      implicit NONE

      integer,       intent(in)  :: mfiles  ! max. number of eta files
                                            ! dynamics file names (eta)
      character*255, intent(out) :: dynfiles(mfiles) 
      integer,       intent(out) :: nfiles  ! actual no. of eta files
      character*255, intent(out) :: RCfile  ! resource file

      integer,       intent(out) :: nymd    ! date of intial perturbation
      integer,       intent(out) :: nhms    ! time of intial perturbation
      integer,       intent(out) :: vectype ! GEOS-4/5 input/output vector-type
      integer,       intent(out) :: idims(2)! im,jm of output, when applicable
      
! !DESCRIPTION: parses command line.
!
! !REVISION HISTORY:
!       15oct2004  Winslow/RT   Initial code, based on similar codes
!       07May2007  Todling      Date/time of pert defined here
!       13Dec2007  Todling      Removed output filename as an optional parameter
!       23Feb2008  Todling      Add -res option
!       21Apr2009  Todling      Updated default resolutions of GEOS-5
!
!EOP

      character*4, parameter :: myname = 'init'

      integer iret, i, iarg, argc, iargc, ires
      logical verb
      character(len=255) :: argv
      character*10 str

      character(len=1) :: res
      integer, dimension(5), parameter :: IMS4 = (/ 72, 144, 288, 576, 1000 /)
      integer, dimension(5), parameter :: IMS5 = (/ 72, 144, 288, 576, 1152 /)
      integer, dimension(5), parameter :: JMS  = (/ 46,  91, 181, 361,  721 /)
      logical :: interp
 

!     Defaults
!     --------
      RCfile   ='initadj.rc'      ! default not to use an RCfile

!     Parse command line
!     ------------------
      argc =  iargc()
      if ( argc .lt. 1 ) call usage()

      iarg   = 0
      nfiles = 0
      nymd   = 0
      nhms   = 0
      vectype = 4    ! default: GEOS-4 input/output vector
      idims   =-1    ! default: don't interpolate output
      interp = .false.

      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) exit
         call GetArg ( iarg, argv )
         select case (argv)
           case ("-g5")
             vectype = 5
           case ('-pick')
               if ( iarg+2 .gt. argc ) call usage()
               iarg = iarg + 1
               call GetArg ( iarg, argv )
               read(argv,*) nymd
               iarg = iarg + 1
               call GetArg ( iarg, argv )
               read(argv,*) nhms
           case ("-rc")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iarg, RCfile )
           case ("-res" )
              if ( iarg+1 .gt. argc ) call usage()
              iarg = iarg + 1
              call GetArg ( iArg, res )
              select case (res)
                case ("a")
                      ires=1
                case ("b")
                      ires=2
                case ("c")
                      ires=3
                case ("d")
                      ires=4
                case ("e")
                      ires=5
                case default
                      print *, 'Sorry this resolution not supported'
                      call exit(1)
              end select
              interp = .true.
           case default
             nfiles = nfiles + 1
             if ( nfiles .gt. mfiles ) call die(myname,'too many eta files')
             dynfiles(nfiles) = argv
         end select
      end do
 
      if ( nfiles .lt. 1 ) call usage()

      if ( interp ) then
           if ( vectype==4 ) then
             idims(1) = ims4(ires)
           else
             idims(1) = ims5(ires)
           endif
           idims(2) = jms (ires)
      endif

      end subroutine Init_
!EOC

!.................................................................

      subroutine usage()
      print *
      print *, '  -------------------------------------------------------------------'
      print *, '  initadj - generates initial conditions for adj/sensitivity program '
      print *, '  -------------------------------------------------------------------'
      print *
      print *
      print *,'Usage: '
      print *,'  initadj.x [-rc RCfile] dynfiles '
      print *
      print *, 'where'
      print *
      print *, '-rc   RCfile   resource filename '
      print *, '              (default: fvpert.rc)'
      print *, '-g5   specify when input/output vectors are GEOS-5 compliant'
      print *
      print *, '        Writes out J-norm into a file Jnormf.txt ' 
      print *, ' (values for total, u, v, T, and ps components of J-norm '
      print *, '              in format 5f20.12)'
      print *
      print *, ' Last updated: 02 Feb 2007; Todling '
      print *
      call MP_finalize(ier)
            if(ier/=0) call MP_die(myname,'MP_finalized()',ier)
      call exit(1)
      end subroutine usage
      
!.................................................................

 end program initadm
