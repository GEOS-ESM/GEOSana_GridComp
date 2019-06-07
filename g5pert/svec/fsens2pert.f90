
  program fsens2pert 

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !PROGRAM: fsens2pert - Performs transformation of ADM/SENS run output 
!                        into a perturbation dynamic vector
!
! !USES:
                                                                                                                           
  use m_mpif90,only : MP_init
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world

  use m_pertutil, only : pertutil_setparam
  use m_fsens2pert, only : fsens2pert_set
  use m_fsens2pert, only : fsens2pert_run

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
!  01Oct2004  Winslow    Hacked to work with adjoint initialization
!  08Aug2005  RG/EN      Adopted initadj.f90 initial code to create
!                        a stand-alone post-processing transformation 
!                        of ADM/SENS outputs into perturbation vectors
!  02Feb2007  Todling    Placed mpi hooks (only runs in 1 pe)
!  08Dec2008  Todling    Unification of codes
!
!EOP
!-------------------------------------------------------------------------
 
  character(len=*), parameter :: myname= 'FSENS2PERT'

  integer, parameter :: mfiles = 4        ! max no. of dyn-vector and txt files
  integer, parameter :: ROOT = 0 ! MPI ROOT PE

  character(len=255) :: dynfiles(mfiles)  ! input dyn filename(s)
  character(len=255) :: pertfile          ! output perturbation filename
  character(len=255), allocatable :: knames(:)
  character(len=255) :: RCfile

  integer :: nfiles  ! actual no. of dyn files read in
  integer :: knmax   ! total no. of dyn-like arrays
  integer :: nymd, nhms
  integer :: myID, ier
  integer :: vectype

! MPI Initialization
! ------------------
  call MP_init(ier)
        if(ier/=0) call MP_die(myname,'MP_init()',ier)
  call MP_comm_rank(MP_comm_world,myID,ier)
        if(ier/=0) call MP_die(myname,'MP_comm_rank()',ier)

! Parse command line
! ------------------
  call init_ ( mfiles, nfiles, dynfiles, pertfile, RCfile, vectype, nymd, nhms )


! Set up program's options
! ------------------------
  call fsens2pert_set ( MP_comm_world, root, RCfile, vectype, ier )
     if(ier/=0) call die(myname,'could not setup run')

  knmax = nfiles + 1
  allocate ( knames(knmax) )

! Generate initial perturbation
! -----------------------------
  call fsens2pert_run ( nfiles, dynfiles, pertfile, knames, nymd, nhms )

! Clean up
! --------
  deallocate ( knames )

  call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)
 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Init_ --- Initialize fsens2pert 
!
! !DESCRIPTION: parses command line.
!
! !INTERFACE:
!
      subroutine Init_ ( mfiles, nfiles, dynfiles, pertfile, RCfile, &
                         vectype, nymd, nhms )

      implicit NONE

      integer,       intent(in)  :: mfiles  ! max. number of eta files
                                            ! dynamics file names (eta)
      character*255, intent(out) :: dynfiles(mfiles) 
      character*255, intent(out) :: pertfile! perturbation filename
      character*255, intent(out) :: RCfile
      integer,       intent(out) :: nfiles  ! actual no. of eta files
      integer,       intent(out) :: nymd
      integer,       intent(out) :: nhms
      integer,       intent(out) :: vectype
      
!
! !REVISION HISTORY:
!       15oct2004  Winslow/RT   Initial code, based on similar codes
!       08Aug2005  RG/EN        Adopted initadj.f90 for this stand-alone application
!
!EOP
!BOC

      character*4, parameter :: myname = 'init'

      integer iret, i, iarg, argc, iargc
      character(len=255) :: argv, res
      character(len=255) :: sensfile, ref_state_file, Jnorm_file
      character*10 str
      logical      pick

!     Defaults
!     --------
      pertfile ='fpert.nc4'             ! default filename for output pert
      sensfile = 'fsens.nc4'            ! default input perturbation
      rcfile = 'NONE'
      ref_state_file = 'ref.nc4'        ! default reference state
      Jnorm_file = 'Jnormf.txt'
      pick = .false.
      nymd = 0
      nhms = 0
      vectype = 4

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

         if (index(argv,'-g5' ) .gt. 0 ) then
            vectype = 5
         else if (index(argv,'-o' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, pertfile )
            nfiles = nfiles + 1
         else if (index(argv,'-sens' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, sensfile )
            nfiles = nfiles + 1
         else if (index(argv,'-ref' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, ref_state_file )
            nfiles = nfiles + 1
         else if (index(argv,'-rc' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, RCfile )
         else if (index(argv,'-Jnorm' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, Jnorm_file )
            nfiles = nfiles + 1
         else if (index(argv,'-pick' ) .gt. 0 ) then
            pick = .true.
            if ( iarg+2 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) nymd
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) nhms
         end if

             if ( nfiles .gt. mfiles ) call die(myname,'too many eta files')
      end do

             if ( nfiles .lt. 2 ) then
                call usage()
                call die(myname,'BOTH SENS and REF files are required')
             endif

            dynfiles(1) = trim(sensfile)
            dynfiles(2) = trim(ref_state_file)
            dynfiles(3) = trim(pertfile)
            dynfiles(4) = trim(Jnorm_file)


      end subroutine Init_

!.................................................................

      subroutine usage()
      print *
      print *, '  -------------------------------------------------------------------'
      print *, '  fsens2pert - generates a perturbation vector from adj/sensitivity outputs '
      print *, '  -------------------------------------------------------------------'
      print *
      print *
      print *,'Usage: '
      print *,'  fsens2pert.x [-o pertfile ] -sens Sens_file -ref Ref_State [-pick nymd nhms] '
      print *,'               [-Jnorm Jnorm_file] [-Temp]'
      print * 
      print *, 'where'
      print *
      print *, '-o      pertfile     output perturbation vector filename'
      print *, '              (default: fsens2pert.nc4)'
      print *
      print *, '-ref    Ref_State    reference state dynamic vector filename'
      print *, '              (default: ref_state.nc4)'
      print *
      print *, '-pick nymd nhms   date and time at the reference state vector'
      print *, '                    (default: use latest on file)'
      print *
      print *, '-sens   Sens_file   input ADM/SENS vector filename'
      print *, '              (default: fsens.nc4)'
      print *
      print *, '-Jnorm  Jnorm_file  filename where a value of J-norm for forecast error is stored'
      print *, '                    in format (e12.4). This value will be used to calculate scale factor.'
      print *, '              (default: the scale factor is set to 1.0)'
      print *
      print *, '-Temp             logical flag if conversion Theta to T is required'
      print *, '                    in output pertfile'
      print *, '              (default: output file will be in THETA)'
      print * 
      print *, ' Last updated: 06 Mar 2009; Todling '
      print *
      call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)
      call exit(1)
      end subroutine usage
      
!.................................................................

 end program fsens2pert 
