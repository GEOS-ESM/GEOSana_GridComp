!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS        !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fvpert: Main program for fvTLM integration

! !INTERFACE:
 
program fvpert

! !USES:

#ifdef SPMD
  use mod_comm,only : mp_init
  use mod_comm,only : mp_exit
#else
  use m_mpif90,only : MP_init
#endif
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world

  use m_gAtlm, only : gAtlm_init
  use m_gAtlm, only : gAtlm_set
  use m_gAtlm, only : gAtlm_run
  use m_gAtlm, only : gAtlm_dotp
  use m_gAtlm, only : gAtlm_wrap
  use m_gAtlm, only : gAtlm_clean

  use m_ioutil, only : luavail

  use m_zeit, only: zeit_ci,zeit_co,zeit_flush
  use m_zeit, only: zeit_allflush

  use m_die,  only : MP_die
  use m_die,  only : die

  implicit none

! !DESCRIPTION: Calculates singular vectors of FVGCM.
!
! !REVISION HISTORY:
!
!  27Jul2004  Todling   Initial code, based on FastOpt "prgadmtlm".
!  07Jun2007  Todling   Add option for checking ADM/TLM
!  28Aug2007  Todling   Various opts to control g5-pert evolution
!  11Dec2007  Todling   Add option to load trajectory to memory
!  18Jul2008  Todling   Add checkpoint opts; add differencing integration
!
!EOP
!-----------------------------------------------------------------------

! Define parameters
! -----------------
  character(len=*), parameter :: myname = 'FVPert'
  integer, parameter :: ROOT = 0

  character(len=255) :: oprefix

  integer n    ! number of independent (initial state) variables
  integer m    ! number of   dependent (final   state) variables
  integer ierr
  integer myID,seed
  integer lu
  integer icase,ndotp
  integer vectype, freq
  logical donlm,oneshot,memtraj,differencing
  integer checkpoint(2)

! MPI Initialization
! ------------------
#ifdef SPMD
  call mp_init
#else
  call MP_init(ierr)
	if(ierr/=0) call MP_die(myname,'MP_init()',ierr)
#endif
  call MP_comm_rank(MP_comm_world,myID,ierr)
	if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

  call init_ ( donlm, ndotp, seed, vectype, oneshot, memtraj, freq, oprefix, differencing, checkpoint )

    call zeit_ci ('fvpert')
  call gAtlm_init  ( ierr, verbose=.true., oneshot=oneshot, freq=freq, memtraj=memtraj, differencing=differencing )
  call gAtlm_set   ( ierr, vectype=vectype )
	if(ierr/=0) call MP_die(myname,'set()',ierr)
  if ( ndotp>0 ) then
       do icase = 0, ndotp-1
          call gAtlm_dotp  ( seed+icase, ierr, checkpoint_tl=checkpoint(1), checkpoint_ad=checkpoint(2) )
       enddo
  else
       call gAtlm_run   ( ierr, nonlinear=donlm )
  endif
  if(trim(oprefix)=='x') then
     call gAtlm_wrap  ( ierr )
  else
     call gAtlm_wrap  ( ierr, fname=oprefix )
  endif
  call gAtlm_clean ( ierr )
    call zeit_co ('fvpert')

    if(myid==ROOT) call zeit_flush(6,subname_at_end=.true.)
                   call zeit_allflush(MP_comm_world,0,6,subname_at_end=.true.)


! Echo run status and close MPI
! -----------------------------
  if(myid==ROOT) print *, myname, ': Normal Execution.'
  lu = luavail()
  open(lu,file='status.log')
  write(lu,*) 'Normal Execution.'
  close(lu)
#if defined ( SPMD )
  call mp_exit
#else
  call MP_finalize(ierr)
	if(ierr/=0) call MP_die(myname,'MP_finalized()',ierr)
#endif


contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: init: initialize odsshuffle
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine init_ ( nlm, ndotp, seed, vectype, oneshot, memtraj, &
                         freq, oprefix, differencing, checkpoint )

! !INPUT/OUTPUT PARAMETERS:
!
      implicit NONE
      logical, intent(out)  :: nlm     ! 0=TLM; 1=two nonlinear integrations
      integer, intent(out)  :: ndotp   ! when !=0, performs test of ADM/TLM
      integer, intent(out)  :: seed    ! random number seed
      integer, intent(out)  :: vectype ! type of input perturbation (g4/g5)
      logical, intent(out)  :: oneshot ! defines whether to run TLM in single shot or not
      logical, intent(out)  :: memtraj ! when .t., load traj to memory
      integer, intent(out)  :: freq    ! freq of output perturbation
      character(len=*), intent(out) :: oprefix ! name of output file
      logical, intent(out)  :: differencing ! when .t., replaces TLM by finite-diff NLM
      integer, intent(out)  :: checkpoint(2)   ! checkpoint:
                                               !  checkpoint(1): controls TLM
                                               !  0=    no checkpoint
                                               !  1=       checkpoint of initial fwrd run
                                               !  checkpoint(2): controls ADM
                                               !  0=    no checkpoint of outer traj
                                               !  1= apprx checkpoint of outer traj
                                               !  2= exact checkpoint of outer traj
!
!
! !REVISION HISTORY:
!
!	16Sep2006 Todling  Initial code.
!       07Jun2007 Todling  Add flag to allow dot product check of ADM/TLM
!       10Jul2007 Todling  Add seed opt
!       28Aug2007 Todling  Add support for g5 perturbation; add oneshot
!
!EOP

      character*4, parameter :: myname = 'init_'

      integer iret, i, lv, iarg, argc, iargc
      character*255 argv,str
      character*10  SS

!     Set defaults
!     ------------
      nlm  = .false.
      ndotp = 0
      seed = -1
      vectype = 4
      oneshot = .true.
      oprefix =  'x'
      freq  = -1
      memtraj  = .false.
      checkpoint = 0
      differencing = .false.

!     Parse command line
!     ------------------
      argc =  iargc()
      if ( argc .lt. 1 ) return
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
         if (index(argv,'-h' ) .gt. 0 ) then
             call usage_()
         endif
         if (index(argv,'-nlm' ) .gt. 0 ) then
            nlm = .true.
         endif
         if (index(argv,'-dotp' ) .gt. 0 ) then
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) ndotp
         endif
         if (index(argv,'-memtraj' ) .gt. 0 ) then
            memtraj = .true.
         endif
         if (index(argv,'-fdiff' ) .gt. 0 ) then
            differencing = .true.
         endif
         if (index(argv,'-g5' ) .gt. 0 ) then
            vectype = 5
         endif
         if (index(argv,'-inc' ) .gt. 0 ) then
            oneshot = .false.
         endif
         if (index(argv,'-checkpoint') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage_()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) checkpoint(1)
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) checkpoint(2)
         endif
         if (index(argv,'-freq') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage_()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) freq
         endif
         if (index(argv,'-o' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage_()
            iarg = iarg + 1
            call GetArg ( iArg, oprefix )
         endif
         if (index(argv,'-seed' ) .gt. 0 ) then
            iarg = iarg + 1
            call GetArg ( iarg, str )
            read(str,*) seed
         endif
      end do
 111  continue

      return

      end subroutine init_

      subroutine usage_()
      integer ierr
      if(myid==ROOT) then
      print *
      print *, 'Usage:  fvpert [-h] [-nlm] [-dotp #] [-g5] [-inc] [-checkpoint TLMLEV ADMLEV]'
      print *
      print *, 'where'
      print *
      print *,   '-h         echoes this Usage statement'
      print *,   '-nlm       evolves perturbation with two nonlinear model integrations'
      print *,   '-dotp N    performs dop-product N-times check for ADM/TLM validity'
      print *,   '-seed #    seed to random # generator; seed=-1, use input perturbation'
      print *,   '-g5        to be specified when input perturbation is GEOS-5-like'
      print *,   '-inc       runs TLM incrementally (step-by-step)'
      print *,   '-memtraj   load trajectory to memory'
      print *,   '-fdiff     replaces TLM with finite-differencing NLM (not tangent)'
      print *,   '-checkpoint TLMLEV ADMLEV  controls checkpoint levels for TLM and ADM'
      print *,   '                           TLMLEV=0, non-zero; ADMLEV=0,1,2'
      print *,   '-freq HHMMSS frequency for output perturbation'
      print *,   '-o PREFIXFN  prefix of filename for propagated perturbation'
      print *
      endif
      call MP_finalize(ierr)
	if(ierr/=0) call MP_die(myname,'MP_finalized()',ierr)
      stop
      end subroutine usage_

end program fvpert
