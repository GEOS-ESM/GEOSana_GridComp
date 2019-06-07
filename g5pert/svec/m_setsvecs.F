!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_SetsVecs --- Set up parameters for SV and Norms
!
! !INTERFACE:

      module m_SetsVecs

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type

      use m_stdio
      use m_inpak90
      use m_die, only: MP_die, die

      use m_sVecDef

      implicit NONE
 
! !PUBLIC MEMBER FUNCTIONS:
 
      PRIVATE

      PUBLIC Set_sVecs

      interface Set_sVecs ; module procedure
     .   SetsVecs0_,
     .   SetsVecs1_,
     .   SetsVecs2_,
     .   SetsVecs3_
      end interface

!
! !DESCRIPTION: Set up for singular vector calculation.
!
! !REVISION HISTORY:
!
!  24Oct2002  Todling   Modularized; split arpack-related initialization.
!  24Jun2002  Gelaro    Added nag-related initialization.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'm_SetsVecs'

      CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: SetsVecs0_ --- General paramater setup for singular vector calc
! 
! !INTERFACE:
!
      subroutine SetsVecs0_ ( comm, root )
 
! !USES: 
    
      Implicit None

! !INCLUDES:

! !INPUT PARAMETERS: 
!
      integer, intent(in) ::  comm	! MPI communicator
      integer, intent(in) ::  root	! MPI ROOT PE

! !OUTPUT PARAMETERS:
!

! !FILES USED:  fvsvec.rc
!
! !DESCRIPTION:  Initializes basic parameters used during the 
!                calculation of analysis errors.
!  Tests:
!     test\_norm=0  : perform no tests of norm
!     test\_norm=1  : compare ynorm $=y^T~S_2~y$ with znorm $=z^T~O~z$
!                    These should be identical except for roundoff. Values of
!                    diff\_norm=ynorm-znorm that are 10*(-10) smaller than
!                    ynorm are acceptable
!     test\_norm=2  : like test\_norm=1, except ADM and TLM replaced by identity
!
! !REVISION HISTORY: 
!
!   24Mar1998  Todling    Initial code.
!   28Oct2002  Todling    Adapted from PSAS-sigA code.
!   14Mar2003  Todling    Moved test\_norm def from svecdrv here.
!   22Apr2004  Todling    Created after original setup.
!   24Jun2004  Gelaro     Made setup generic, with NAG,ARPACK as choice
!   26Jul2004  Todling    Moved proj-coeffs setting to m_fvproj.
!   19Dec2005  Gelaro     Modified for Lanczos restart capability
!   13Dec2007  Todling    Add Ehrendorfer et al. (1999) parameter eps_eer
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'SetsVecs0_'
      character(len=*), parameter :: def_svecrc = 'fvsvec.rc'

      character(len=255) token
      character(len=255) svecrc
      character(len=3  ) evolve_svec

      integer       i, j, k, iret, ierr
      integer       ival
      integer       myID

      call MP_comm_rank(comm,myID,ierr)
        if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)

!     Load resources from fvsvec.rc
!     ---------------------------
      svecrc = ' '
      call getenv('SIGA_RC',svecrc)		! Unix binding
      if(svecrc.eq.' ') svecrc=def_svecrc	! default name
      call i90_loadf (trim(svecrc), iret)
      if( iret .ne. 0) then
      write(stderr,'(2a,i5)') myname_,': I90_loadf error, iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) then
         write(stdout,'( a  )') '---------------------------------'
         write(stdout,'(2a  )') myname_, ': Reading resource file'
         write(stdout,'( a,/)') '---------------------------------'
      end if


!     Select eigensolver: 1=ARPACK, 2=NAG
!     -------------------------------------
      call I90_label('eigensolver_package:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,': I90_label error, iret =',iret
        call die(myname_)
      end if
      isolve = I90_GInt(iret)
      if( iret .ne. 0) then
         write(stderr,'(2a,i5)') myname_,': I90_GInt error, iret =',iret
         call die(myname_)
      end if
      if ( isolve .eq. 1 ) then
        if(myID==ROOT) write(stdout,'(a,i2,a)') 'Eigensolver: ',isolve, ' ARPACK'
      else if ( isolve .eq. 2 ) then
        if(myID==ROOT) write(stdout,'(a,i2,a)') 'Eigensolver: ',isolve, ' NAG'
      else if ( isolve .eq. 3 ) then
        if(myID==ROOT) write(stdout,'(a,i2,a)') 'Will solve CNOP problem'
      else
        if(myID==ROOT) write(stdout,'(2a,i2,a)') myname_, ': Invalid eigensolver: ',isolve, ' ...must choose 1=ARPACK or 2=NAG'
        call die(myname_)
      end if

!     Define type of eigen-problem
!     ----------------------------
      svecnormF = 'KE'   ! default: KE
      call I90_label('svec_norm:', iret)
      if (iret .eq. 0) then
        call I90_Gtoken ( token, iret )
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,
     .                 ': I90_Gtoken error, iret =',iret
        else
          svecnormF = trim(token)
        end if
      end if
      if(myID==ROOT) write(stdout,'(2a)') 'Final svec norm: ', trim(svecnormF)

!     Define initial normalization
!     ----------------------------
      svecnormI = svecnormF   ! default: same as final norm
      call I90_label('initial_svec_norm:', iret)
      if (iret .eq. 0) then
        call I90_Gtoken ( token, iret )
        if (iret .ne. 0) then
          write(stdout,'(2a,i5)') myname_,
     .                 ': not defined, taking default ...'
        else
          svecnormI = trim(token)
        end if
      end if
      if(myID==ROOT) write(stdout,'(2a)') 'Initial svec norm: ', trim(svecnormI)

!     Determine whether to do norm test
!     ---------------------------------
      test_norm = 0  ! default: no testing of norms
      call I90_label('do_norm_test:', iret)
      if (iret .eq. 0) then
        ival = I90_GInt(iret)
        if (iret .ne. 0) then
          write(stdout,'(2a,i5)') myname_,
     .                 ': not defined, taking default ...'
        else
          test_norm = ival
        end if
      end if
      if(myID==ROOT) write(stdout,'(a,i5)') 'Test norm value: ', test_norm


!     Read option to calculate eigenvectors
!     -------------------------------------
      call I90_label('calculate_eigenvectors:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,
     .               ': I90_label error, iret =',iret
        call die(myname_)
      end if
      call I90_Gtoken ( token, iret )
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,
     .               ': I90_Gtoken error, iret =',iret
        call die(myname_)
      end if
      rvec = .false.
      if( trim(token) .eq. 'yes' ) rvec = .true.
      if(myID==ROOT) write(stdout,'(2a)')
     .             'Calculate eigenvectors: ', trim(token)


!     Decide whether to evolve svecs or not
!     -------------------------------------
      propsvec = .false.   ! default: no propagation
      call I90_label('evolve_svec:', iret)
      if (iret .eq. 0) then
        call I90_Gtoken ( token, iret )
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,
     .                 ': I90_Gtoken error, iret =',iret
        else
          evolve_svec = trim(token)
        end if
      end if
      if ( trim(evolve_svec)=='yes' .or. trim(evolve_svec)=='y' ) propsvec = .true.
      if(myID==ROOT) write(stdout,'(2a)') 'Evolve svecs: ', trim(evolve_svec)

!     Read maximum number of Lanczos iterations per job
!     -------------------------------------------------
      call I90_label('maximum_iterations_per_job:', iret)
      if( iret .ne. 0) then
         maxitr = maxd
         if (myID==0) then
           write(stderr,'(2a)') myname_,': I90_label, label not found'
           write(stdout,'(a,i8)') 'Using default maxitr: ', maxitr
         endif ! < ROOT >
      else
        maxitr = I90_GInt(iret)
        if ( iret .ne. 0) then
           maxitr = maxd
           if (myID==ROOT) then
             write(stderr,'(2a)') myname_,': I90_GInt error, value not found'
             write(stdout,'(a,i8)') 'Using default maxitr: ', maxitr
           endif ! < ROOT >
        else
           if ( maxitr .lt. 0 ) maxitr = maxd
           if(myID==ROOT) write(stdout,'(a,i7)') 'Max iterations for this job: ',maxitr
        end if
      end if

!     Read Ehrendorfer, Errico, and Raeder's epsilon factor
!     -----------------------------------------------------
      eps_eer = 1.0d0
      call I90_label('ehrendorfer_errico_raedder_eps:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_, ': I90_label error, iret =',iret
      else
        eps_eer = I90_GFloat(iret)
        if( iret .ne. 0) then
           write(stderr,'(3a,i5)') myname_,': I90_GFloat error, ', ' iret =',iret
           call die(myname_)
        end if
      end if
      if(myID==ROOT) write(stdout,'(a,e13.6)') 'Ehrendorfer, Errico, and Raeder eps: ',eps_eer


!     release resource file:
!     ---------------------
      call I90_release()

      return
      end subroutine SetsVecs0_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: SetsVecs1_ --- Initialize ARPACK-specific paramters
! 
! !INTERFACE:
!
      subroutine SetsVecs1_ ( comm, root, nev, ncv, which )
 
! !USES: 
    
      Implicit None

! !INCLUDES:

      include 'debug.h'                  ! ARPACK include for debugging/output purposes

! !INPUT PARAMETERS: 
!
      integer, intent(in) ::  comm	 ! MPI communicator
      integer, intent(in) ::  root	 ! MPI ROOT PE

! !OUTPUT PARAMETERS:
!

      integer,          intent(out) :: nev     ! number of eigenvalues/eigenvectors to calculate
      integer,          intent(out) :: ncv     ! number of Lanczos basis vectors
      character(len=2), intent(out) :: which   ! leading/trailing/both eigenvalues

! !FILES USED:  fvsvec.rc
!
! !DESCRIPTION:  Initializes ARPACK-specific parameters.
!
! !REVISION HISTORY: 
!
!   24Mar1998  Todling    Initial code.
!   28Oct2002  Todling    Adapted from PSAS-sigA code.
!   14Mar2003  Todling    Moved test\_norm def from svecdrv here.
!   22Apr2004  Todling    Created after split of original setup.
!   19Dec2005  Gelaro     Modified for Lanczos restart capability
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'SetsVecs1_'
      character(len=*), parameter :: def_svecrc = 'fvsvec.rc'

      character(len=255) token
      character(len=255) svecrc

      integer       i, j, k, iret, ierr
      integer       maxupd
      integer       ival
      integer       myID

      call MP_comm_rank(comm,myID,ierr)
        if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)

!     Setup debug options of ARPACK
!     -----------------------------
      ndigit = -3
      logfil = 6
      msgets = 0
      msaitr = 0 
      msapps = 0
      msaupd = 1
      msaup2 = 0
      mseigt = 0
      mseupd = 0
 
      pflabel = 'svec'  ! perturbation file tag name

!     Load resources from fvsvec.rc
!     -----------------------------
      svecrc = ' '
      call getenv('SIGA_RC',svecrc)		! Unix binding
      if(svecrc.eq.' ') svecrc=def_svecrc	! default name
      call i90_loadf (trim(svecrc), iret)
      if( iret .ne. 0) then
      write(stderr,'(2a,i5)') myname_,': I90_loadf error, iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) then
         write(stdout,'( a  )') '---------------------------------'
         write(stdout,'(2a  )') myname_, ': Reading resource file'
         write(stdout,'( a,/)') '---------------------------------'
      end if


!     Read number of eigenvalues to calculate
!     ---------------------------------------
      call I90_label('number_eigenvectors:', iret)
      if( iret .ne. 0) then
         nev = nevd
         if (myID==0) then
           write(stderr,'(2a)') myname_,': I90_label, label not found'
           write(stdout,'(a,i8)') 'Using default (nev): ', nev
         endif ! < ROOT >
      else
        nev = I90_GInt(iret)
        if ( iret .ne. 0) then
           nev = nevd
           if (myID==ROOT) then
             write(stderr,'(2a)') myname_,': I90_GInt error, value not found'
             write(stdout,'(a,i8)') 'Using default (nev): ', nev
           endif ! < ROOT >
        else
           if ( nev .lt. 0 ) nev = nevd
           if(myID==ROOT) write(stdout,'(a,i7)') 'Number of singular values(nev): ',nev
        end if
      end if

!     Read number of vectors for Lanczos basis
!     ----------------------------------------
      call I90_label('number_lanczos_basis_vectors:', iret)
      if (iret .ne. 0) then
        ncv = ncvd
        write(stderr,'(2a)') myname_,
     .                 ': I90_label error, label not found '
        write(stdout,'(a,i8)') 'Using default (ncv): ', ncv
      else
        ival = I90_GInt(iret)
        if( iret .ne. 0 ) then
           ncv = ncvd
           write(stderr,'(2a,i5)') myname_,': I90_GInt error, iret =',iret
           write(stdout,'(a,i8)') 'Using default (ncv): ', ncv
        else
           ncv = max(ival,nev+1)
           if(myID==ROOT)
     .     write(stdout,'(a,i8)') 'Number of Lanczos basis vectors (ncv): ', ncv
        end if
      endif

!     Which eigenvalues: leading(LM), trailing(SM) or both(BE) 
!     --------------------------------------------------------
      which = 'LM'   ! default is leading eigenvalues
      call I90_label('which_eigenvalues:', iret)

      if (iret .eq. 0) then
        call I90_Gtoken ( token, iret )

        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,
     .                 ': I90_Gtoken error, iret =',iret
          write(stdout,'(a)') 'Using default (which=LM): Leading eigenvalues'
        else
          if (( trim(token) == 'LM' .or. trim(token) == 'lm') .or.
     .        ( trim(token) == 'SM' .or. trim(token) == 'sm') .or.
     .        ( trim(token) == 'BE' .or. trim(token) == 'be') ) then
             which = trim(token)
             if(myID==ROOT) write(stdout,'(2a)') 'Which eigenvalues: ', which
          else
             write(stderr,'(2a)') myname_,': Invalid choice (which): Must be LM,SM or BE'
             write(stdout,'(a)') 'Using default (which=LM): Leading eigenvalues'
          end if
        end if

      else
         write(stderr,'(2a)') myname_,': Which eigenvalues (LM,SM or BE) not specified'
         write(stdout,'(a)') 'Using default (which=LM): Leading eigenvalues'
      end if

!     Read relative accuracy of eigenvalue calculation
!     ------------------------------------------------
      call I90_label('eigenvalue_relative_accuracy:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,
     .                 ': I90_label error, iret =',iret
        call die(myname_)
      end if
      tol = I90_GFloat(iret)
      if( iret .ne. 0) then
         write(stderr,'(3a,i5)') myname_,': I90_GFloat error, ',
     .                                  ' iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) write(stdout,'(a,e13.6)') 'Tolerance: ',tol


!     Read maximum allowed number of Arnoldi update iterations
!     --------------------------------------------------------
      call I90_label('maximum_arnoldi_iterations:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,
     .                 ': I90_label error, iret =',iret
        call die(myname_)
      end if
      maxupd = I90_GInt(iret)
      if( iret .ne. 0) then
         write(stderr,'(2a,i5)') myname_,': I90_GInt error, iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) write(stdout,'(a,i7)') 'Maximum allowed number of Arnoldi update iterations: ',maxupd


!     Read option to check on accuracy of eigendecomposition
!     ------------------------------------------------------
      call I90_label('eigen_decomposition_accuracy:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,
     .               ': I90_label error, iret =',iret
        call die(myname_)
      end if
      call I90_Gtoken ( token, iret )
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,
     .               ': I90_Gtoken error, iret =',iret
        call die(myname_)
      end if
      eigchk = .false.
      if( trim(token) .eq. 'yes' ) eigchk = .true.
      if(myID==ROOT) write(stdout,'(2a)') 
     .             'Check on accuracy of eigen-decomposition: ',
     .              trim(token)


!     release resource file:
!     ---------------------
      call I90_release()

!     now do some setting up:
!     ----------------------
      do j = 1, 11
         iparam(j) = 0
      end do
      iparam(1) = 1		!  shifting option; set to exact shifts
      iparam(3) = maxupd	!  maximum number of Arnoldi update iterations
      if ( trim(svecnormI) == 'Pa' .or. trim(svecnormI) == 'pa' ) then
           iparam(7) =  2	!  type of eigen-problem: generalized eigen-problem
           bmat      = 'G'
           if ( myID==ROOT ) 
     .          write(stdout,'(a)') '*** Will use Pa-based norm '
      else if (( trim(svecnormI) == 'L2' .or. trim(svecnormI) == 'l2') .or.
     .         ( trim(svecnormI) == 'S2' .or. trim(svecnormI) == 's2') .or.
     .         ( trim(svecnormI) == 'Da' .or. trim(svecnormI) == 'da') .or.
     .         ( trim(svecnormI) == 'KE' .or. trim(svecnormI) == 'ke') .or.
     .         ( trim(svecnormI) == 'TE' .or. trim(svecnormI) == 'te')    ) then
           iparam(7) =  1	!  type of eigen-problem: standard eigen-problem
           bmat      = 'I'
      else
           write(stderr,'(2a)') myname_,
     .                   ': No such scheme implemented. Aborting ...'
           call die(myname_)
      end if
 
      if(myID==ROOT)then
      write(stdout,'(/,2a)') myname_, ': calculation setup completed '
      write(stdout,'(a,/)')  '-------------------------------------- '
      endif
 
      return
      end subroutine SetsVecs1_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: SetsVecs2_ --- Initialize NAG-specific parameters
! 
! !INTERFACE:
!
      subroutine SetsVecs2_ ( comm, root, lanmax )
 
! !USES: 
    
      Implicit None

! !INCLUDES:


! !INPUT PARAMETERS: 
!
      integer, intent(in) ::  comm	 ! MPI communicator
      integer, intent(in) ::  root	 ! MPI ROOT PE

! !OUTPUT PARAMETERS:
!

      integer, intent(out) :: lanmax	 ! number of Lanczos steps to be taken

! !FILES USED:  svec.rc
!
! !DESCRIPTION:  Initializes NAG-specific parameters.
!
! !REVISION HISTORY: 
!
!   04May2004  Gelaro     Created to initialize NAG solver, based on SetsVecs1_
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'SetsVecs2_'
      character(len=*), parameter :: def_svecrc = 'fvsvec.rc'

      character(len=255) token
      character(len=255) svecrc

      integer       i, j, k, iret, ierr
      integer       ival
      integer       myID

      call MP_comm_rank(comm,myID,ierr)
        if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)


      pflabel = 'svec'  ! perturbation file tag name

!     Load resources from svec.rc
!     ---------------------------
      svecrc = ' '
      call getenv('SIGA_RC',svecrc)		! Unix binding
      if(svecrc.eq.' ') svecrc=def_svecrc	! default name
      call i90_loadf (trim(svecrc), iret)
      if( iret .ne. 0) then
      write(stderr,'(2a,i5)') myname_,': I90_loadf error, iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) then
         write(stdout,'( a  )') '---------------------------------'
         write(stdout,'(2a  )') myname_, ': Reading resource file'
         write(stdout,'( a,/)') '---------------------------------'
      end if


!     Read number of Lanczos iterations (applications of Oseledec operator)
!     ---------------------------------------------------------------------
      call I90_label('number_lanczos_iterations:', iret)
      if (iret .ne. 0) then
        lanmax = lanmaxd
        write(stderr,'(2a)') myname_,
     .                 ': I90_label error, label not found '
        write(stdout,'(a,i8)') 'Using default (lanmax): ', lanmax
      else
        lanmax = I90_GInt(iret)
        if( iret .ne. 0 ) then
           lanmax = lanmaxd
           write(stderr,'(2a,i5)') myname_,': I90_GInt error, iret =',iret
           write(stdout,'(a,i8)') 'Using default (lanmax): ', lanmax
        else
           if ( lanmax .lt. 1 ) lanmax = lanmaxd
           if(myID==ROOT)
     .     write(stdout,'(a,i8)') 'Number of Lanczos steps (lanmax): ',lanmax
        end if
      endif

!     Read accuracy value for eigenvector calculation
!     ----------------------------------------------
      call I90_label('eigenvector_accuracy:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,
     .                 ': I90_label error, iret =',iret
        call die(myname_)
      end if
      kappa = I90_GFloat(iret)
      if( iret .ne. 0) then
         write(stderr,'(3a,i5)') myname_,': I90_GFloat error, ',
     .                                  ' iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) write(stdout,'(a,e13.6)') 'Eigenvector accuracy: ',kappa


      if(myID==ROOT)then
      write(stdout,'(/,2a)') myname_, ': calculation setup completed '
      write(stdout,'(a,/)')  '-------------------------------------- '
      endif
 
!     release resource file:
!     ---------------------
      call I90_release()

      return
      end subroutine SetsVecs2_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: SetsVecs3_ --- Initialize CNOP-specific parameters
! 
! !INTERFACE:
!
      subroutine SetsVecs3_ ( comm, root, what, ncv )
 
! !USES: 
    
      Implicit None

! !INCLUDES:


! !INPUT PARAMETERS: 
!
      integer, intent(in) ::  comm	 ! MPI communicator
      integer, intent(in) ::  root	 ! MPI ROOT PE

! !OUTPUT PARAMETERS:
!

      integer,          intent(out) :: ncv  ! number of CNOPs to get
      character(len=4), intent(out) :: what ! dummy to trick interface

! !FILES USED:  fvsvec.rc
!
! !DESCRIPTION:  Initializes CNOP-specific parameters.
!
! !REVISION HISTORY: 
!
!   29Aug2006  Todling     Initial code.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'SetsVecs3_'
      character(len=*), parameter :: def_svecrc = 'fvsvec.rc'

      character(len=255) token
      character(len=255) svecrc

      integer       i, j, k, iret, ierr
      integer       ival
      integer       myID

      call MP_comm_rank(comm,myID,ierr)
        if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)

      what    = 'dummy'
      pflabel = 'cnop'  ! perturbation file tag name

!     Load resources from svec.rc
!     ---------------------------
      svecrc = ' '
      call getenv('SIGA_RC',svecrc)		! Unix binding
      if(svecrc.eq.' ') svecrc=def_svecrc	! default name
      call i90_loadf (trim(svecrc), iret)
      if( iret .ne. 0) then
      write(stderr,'(2a,i5)') myname_,': I90_loadf error, iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) then
         write(stdout,'( a  )') '---------------------------------'
         write(stdout,'(2a  )') myname_, ': Reading resource file'
         write(stdout,'( a,/)') '---------------------------------'
      end if


!     Read number of CNOPs to calculate
!     ---------------------------------
      call I90_label('number_of_cnops:', iret)
      if (iret .ne. 0) then
        ncv = 1
        write(stderr,'(2a)') myname_,
     .                 ': I90_label error, label not found '
        write(stdout,'(a,i8)') 'Using default number of CNOPs to get: ', ncv
      else
        ncv = I90_GInt(iret)
        if( iret .ne. 0 ) then
           ncv = 1
           write(stderr,'(2a,i5)') myname_,': I90_GInt error, iret =',iret
           write(stdout,'(a,i8)') 'Using default (ncv): ', ncv
        else
           if(myID==ROOT) write(stdout,'(a,i8)') 'Number of CNOPs to get: ',ncv
        end if
      endif


!     Read number of iterations in SPG method
!     ---------------------------------------
      call I90_label('number_spg_iterations:', iret)
      if (iret .ne. 0) then
        spg_miter = spg_miterd
        write(stderr,'(2a)') myname_,
     .                 ': I90_label error, label not found '
        write(stdout,'(a,i8)') 'Using default (spg_miter): ', spg_miter
      else
        spg_miter = I90_GInt(iret)
        if( iret .ne. 0 ) then
           spg_miter = spg_miterd
           write(stderr,'(2a,i5)') myname_,': I90_GInt error, iret =',iret
           write(stdout,'(a,i8)') 'Using default (spg_miter): ', spg_miter
        else
           if ( spg_miter .lt. 1 ) spg_miter = spg_miterd
           if(myID==ROOT)
     .     write(stdout,'(a,i8)') 'Number of SPG max iterations (spg_miter): ',spg_miter
        end if
      endif

!     Read in tolerance for norm comparison in SPG method
!     ---------------------------------------------------
      call I90_label('cnop_tolerance:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,
     .                 ': I90_label error, iret =',iret
        call die(myname_)
      end if
      cnop_tol = I90_GFloat(iret)
      if( iret .ne. 0) then
         write(stderr,'(3a,i5)') myname_,': I90_GFloat error, ',
     .                                  ' iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) write(stdout,'(a,e13.6)') 'CNOP norm tolerance: ',cnop_tol

!     Read in magnitude of norm to be acheived by CNOP
!     ------------------------------------------------
      call I90_label('cnop_norm_value:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,
     .                 ': I90_label error, iret =',iret
        call die(myname_)
      end if
      cnop_sigma = I90_GFloat(iret)
      if( iret .ne. 0) then
         write(stderr,'(3a,i5)') myname_,': I90_GFloat error, ',
     .                                  ' iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) write(stdout,'(a,e13.6)') 'CNOP target norm value: ',cnop_sigma
                                                                                                               


      if(myID==ROOT)then
      write(stdout,'(/,2a)') myname_, ': calculation setup completed '
      write(stdout,'(a,/)')  '-------------------------------------- '
      endif
 
!     release resource file:
!     ---------------------
      call I90_release()

      return
      end subroutine SetsVecs3_

      end module m_SetsVecs
