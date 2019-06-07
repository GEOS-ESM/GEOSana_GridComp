        module m_sVecDef

! History:
! -------
!
!  28Oct2002 Todling  - Created.
!  23Jun2004 Gelaro   - Modified for use with NAG or ARPACK solver
!  20Dec2005 Gelaro   - Modified for Lanczos restart capability
!  14Aug2006 Todling  - Added parameter related to CNOP calculation
!  13Dec2007 Todling  - Add eps_eer
!
! -----------------------------------------------------------------


        implicit none

!       Define eigen-problem case
!       -------------------------

!       0. Generic parameters
!       ---------------------
        integer, save ::   isolve        ! selects eignen solver: 1=ARPACK, 2=NAG

	logical, save ::   rvec          ! calculate eigenvectors
	logical, save ::   eigchk        ! check on accuracy of eigen-decomposition

        logical, save ::   propsvec      ! determines whether or not to evolve svecs

        character*2, save :: svecnormI   ! initial state norm
        character*2, save :: svecnormF   ! final state norm for eigenvectors
                                         !   L2  - for L2-norm
                                         !   KE  - for Kinetic Energy norm
                                         !   TE  - for Total   Energy norm
                                         !   Pa  - for analysis error cov-based norm
                                         !   Da  - for analysis error variances norm

                                         ! local projection operator box limits
        logical, save :: lclproj         !   .t. when local proj is applied
        real   , save :: projlon(2)      !   longitudes
        real   , save :: projlat(2)      !   latitudes
        integer, save :: projlev(2)      !   vertical levels

        integer, save :: test_norm       !   levs of norm tests: 0, 1, or 2

        integer, save ::   perc_var      ! percentage variance captured in innov cov

        integer, save ::   maxitr        ! maximum no. of Lanczos iterations per job
        integer, save ::   maxd = 1      ! default maximum no. Lanczos iterations per job

        integer, save ::   ncalls        ! counter for lanczos iterations
        integer, save ::   lanunit       ! unit number for lanczos restart data

        character(len=4), save :: pflabel! perturbation file tag name

        real, save    :: eps_eer  ! eps from Ehrendorfer, Errico and Raeder (1999)
                                  ! this controls the extent to which q
                                  ! influences the "wet-energy" norm

!       1. ARPACK-specific parameters
!       -----------------------------
        integer, save ::   nevd = 1      ! default number of eigenvalues/vectors
        integer, save ::   ncvd = 2      ! default number of Lanczos basis vectors

	character, save ::   bmat	 ! define eigen-problem type in ARPACK sense
 
	integer, save ::   iparam(11)    ! array of specific options for ARPACK
	real   , save ::   tol           ! relative accuracy of eigenvalues


!       2. NAG-specific parameters
!       --------------------------
        integer, save :: lanmaxd = 1     ! default total number of Lanczos iterations

        real   , save ::   kappa         ! relative accuracy of eigenvalues

!       3. CNOP-specific parameters
!       ---------------------------
        integer, save :: spg_mfls  = 1000! number of function evaluations used in SPG line search
	integer, save :: spg_maxfc = 2000! max number of function evaluations
	integer, save :: spg_miterd= 100 ! default max number of iterations
	real   , save :: spg_einf  = 0.0 ! SPG inf norm criterium
        logical, save :: spg_verb  = .true.
	
	integer, save :: spg_miter       ! max number of iterations
        real   , save :: cnop_sigma      ! magnitude of norm of growing perturbation;
                                         ! typical LSV total energy
        real   , save :: cnop_tol        ! tolerance for convergence of CNOPs
	
        end module m_sVecDef
