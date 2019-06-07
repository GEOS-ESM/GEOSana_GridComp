
  program pertenergy 

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !PROGRAM: pertenergy - POstprocess perturbation vector(s) to output energy field 
!                        as a dynamic vector (in *.nc4 file).  If multiple perturbation
!                        vectors are processed, then result is the average energy field.
!
! !USES:

  use precision

  use m_mpif90,only : MP_init
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world

  use m_pertenergy, only : pertenergy_Set
  use m_pertenergy, only : pertenergy_Process
  use m_pertenergy, only : pertenergy_clean

  use m_pertutil, only : pertutil_setparam

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
!  16Aug2005  RG/EN      Adopted postStats.f90 initial code to create
!                           a stand-alone post-processing utility
!  03May2007  Gelaro     Added rcfile capability to process multiple perturbations
!                        a stand-alone post-processing utility
!  02Feb2007  Todling    Placed mpi hooks (only runs in 1 pe)
!  10Jan2007  Todling    - Remove energy-type opt(generalized)
!                        - Add end date/time to command line arg list
!  09Aug2017  Todling    Opt not to time-tag output file
!
!EOP
!-------------------------------------------------------------------------
 
  character(len=*), parameter :: myname= 'pertenergy'

  integer, parameter :: mfiles = 3        ! max no. of dyn-vector and txt files
  integer, parameter :: ROOT = 0 ! MPI ROOT PE

  character(len=255) :: dynfiles(mfiles)  ! input dyn filename(s)
  character(len=255) :: expid             ! name of experiment
  character(len=255) :: rcfile            ! input resource file

  integer :: knmax        ! total no. of dyn-like arrays
  logical :: pick
  logical :: notag
  integer :: myID, ier
  integer :: nymd         ! beginning (or only) date (YYYYMMDD) to process
  integer :: nhms         ! beginning (or only) time (HHMMSS) to process
  integer :: nymd_end     ! ending date (YYYYMMDD) to process
  integer :: nhms_end     ! ending time (HHMMSS) to process


! MPI Initialization
! ------------------
  call MP_init(ier)
        if(ier/=0) call MP_die(myname,'MP_init()',ier)
  call MP_comm_rank(MP_comm_world,myID,ier)
        if(ier/=0) call MP_die(myname,'MP_comm_rank()',ier)

! Parse command line
! ------------------
  call init_ ( mfiles, dynfiles, pick, &
               nymd, nhms, nymd_end, nhms_end, expid, notag, rcfile )


! Set up program options from rcfile
! ----------------------------------
  call pertenergy_Set ( MP_comm_world, root, rcfile, ier )


! Process perturbation file(s)
! ----------------------------
  call pertenergy_Process ( mfiles, dynfiles, expid, notag, nymd, nhms, nymd_end, nhms_end ) 

! Clean up
! --------
  call pertenergy_clean( )

  print *, 'Program done'
  call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)
 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Init_ --- Initialize pertenergy 
!
! !DESCRIPTION: parses command line.
!
! !INTERFACE:
!
      subroutine Init_ ( mfiles, dynfiles, pick, &
                         nymd, nhms, nymd_end, nhms_end, expid, notag, rcfile )

      implicit NONE

      integer,       intent(in)  :: mfiles  ! max. number of eta files
                                            ! dynamics file names (eta)
      character*255, intent(out) :: dynfiles(mfiles) 
      character*255, intent(out) :: expid
      logical,       intent(out) :: pick
      integer,       intent(out) :: nymd
      integer,       intent(out) :: nhms
      integer,       intent(out) :: nymd_end
      integer,       intent(out) :: nhms_end
      logical,       intent(out) :: notag
      
!
! !REVISION HISTORY:
!       15oct2004  Winslow/RT  Initial code, based on similar codes
!       08Aug2005  RG/EN       Adopted initadj.f90 for this stand-alone application
!       10Aug2008  Todling     Add various handles and safe-guards
!
!EOP
!BOC

      character*4, parameter :: myname = 'init'

      integer iret, i, iarg, argc, iargc, vectype, nfiles
      character(len=255) :: argv, res
      character(len=255) :: pertfile, ref_state_file, energy_file, rcfile
      character*10 str
      real eps_eer

!     Defaults
!     --------
      ref_state_file = 'ref_state.nc4'
      energy_file ='pertenergy.eta'    ! default filename for output pert
      pertfile    = 'pert.nc4'
      pick        = .false.
      nymd        = 0
      nhms        = 0
      nymd_end    = -1
      nhms_end    = -1
      rcfile      = 'null'
      expid       = 'null'
      notag       = .false.
      eps_eer     = -999.0

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
            call GetArg ( iArg, energy_file )
            nfiles = nfiles + 1
         else if (index(argv,'-pert' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, pertfile )
            nfiles = nfiles + 1
         else if (index(argv,'-ref' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, ref_state_file )
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
         else if (index(argv,'-eps_eer' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) eps_eer
         else if (index(argv,'-epick' ) .gt. 0 ) then
            pick = .true.
            if ( iarg+2 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) nymd_end
            iarg = iarg + 1
            call GetArg ( iArg, argv )
            read(argv,*) nhms_end
            rcfile = 'pertenergy.rc'   ! reset default name for rc file
         else if (index(argv,'-expid' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, expid )
         else if (index(argv,'-notag' ) .gt. 0 ) then
            notag = .true.
         else if (index(argv,'-rc' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, rcfile )
         end if

             if ( nfiles .gt. mfiles ) call die(myname,'too many eta files')
      end do

      if ( nfiles .lt. 2 .and. rcfile=='null' ) then
           print *
           print *, myname,': BOTH SENS and REF files are required'
           print *
           call usage()
      endif

      if ( expid=='null' .and. rcfile/='null' ) call usage()

      call pertutil_setparam ( 'vectype', vectype )
      if( vectype==5 ) then
          call pertutil_setparam ( 'tname', 'vT_' )
          call pertutil_setparam ( 'wgrid', 'a'   )
      endif

      if ( nymd_end==nhms_end .and. nymd_end==-1 ) then
           nymd_end= nymd
           nhms_end= nhms
      endif
      if ( rcfile=='null' ) then
          dynfiles(1) = pertfile
          dynfiles(2) = ref_state_file
          if (eps_eer < 0.0) then
             print* 
             print*, myname, ': when no RC-file command line eps_eer must be used'
             print* 
             call usage()
          else
             call pertutil_setparam ( 'eps_eer', real(eps_eer,r8) )
          endif
      endif
      dynfiles(3) = energy_file

      end subroutine Init_

!.................................................................

      subroutine usage()
      print *
      print *, '  -------------------------------------------------------------------'
      print *, '  pertenergy - generates a perturbation vector from adj/sensitivity outputs '
      print *, '  -------------------------------------------------------------------'
      print *
      print *
      print *,'Usage: '
      print *,'  To process single file'
      print *,'    pertenergy.x [-o energy ] -pert pert -ref refstate -pick nymd nhms '
      print *
      print *,'  To process multiple files'
      print *,'    pertenergy.x [-o energy ] -expid expid -rc RCfile '
      print *, '                             -pick  nymdb nhmsb '
      print *, '                            -epick  nymde nhmse '
      print * 
      print *, 'where'
      print *
      print *, '-o      energy       prefix for output nc4 filename'
      print *
      print *, '-ref    refstate    reference state dynamic vector filename'
      print *, '                    (default: ref_state.nc4)'
      print *
      print *, '-pick   nymd nhms   date and time at the reference state vector'
      print *, '                    (default: read single-time files in arg list)'
      print *
      print *, '-epick  nymde nhmse end date and time of files to process '
      print *, '                    (default: ignore it)'
      print *
      print *, '-pert   pert        filename of perturbation file'
      print *, '                    (default: pert.nc4)'
      print *, '-eps_eer EPS_EER    must be specified when RC file not given'
      print *
      print *, '-rc     rcfile      resource file to control processing, overrides -ref and -pert'
      print *, '                    REQUIRED to process multiple perturbation files'
      print *, '                    OPTIONAL to process a single perturbation file'
      print *, '                    (default name: pertenegy.rc)'
      print *
      print *, ' Last updated: 14 Sep 2013; Todling '
      print *
      call MP_finalize(ier)
        if(ier/=0) call MP_die(myname,'MP_finalized()',ier)
      call exit(1)
      end subroutine usage
      
!.................................................................

 end program pertenergy 
