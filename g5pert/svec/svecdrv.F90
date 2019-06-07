!--------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!--------------------------------------------------------------------------
!BOI
!
! !TITLE: The GEOS-4/5 Adjoint-based Library
!
! !AUTHORS: R. Todling, R. Gelaro, and R. Errico 
!
! !AFFILIATION: Global Modeling and Assimilation Office, NASA/GSFC, Greenbelt, MD 20771
!
! !DATE: 03 August 2006
!
! !INTRODUCTION: Singular Vector Program
!
!  This document describes the program developed by GMAO to calculate
!  singular values and singular vectors for GEOS-4 and GEOS-5 model
!  applications.
!
!EOI
!-------------------------------------------------------------------------

!--------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!--------------------------------------------------------------------------
!BOP
!
! !ROUTINE: svdrv --- main driver for either NAG- or ARPACK-based singular
!                     vector calculation
!
! !INTERFACE:

subroutine svdrv( ROOT, comm )

! !USES:

  use precision
  use stepon,  only : stepon_set
  use stepon,  only : job
  use stepon,  only : nymd
  use stepon,  only : nhms
  use stepon,  only : fvpsasdt
  use stepon,  only : nstep
  use stepon,  only : ptop
  use stepon,  only : imr,jfirst,jlast,nl ! temporary
  use stepon,  only : ks
  use stepon,  only : ak
  use stepon,  only : bk
  use stepon,  only : ts
  use stepon,  only : oro
  use stepon,  only : ptrjtmpl
  use control, only : mod2cont

  use prognostics, only : prognostics_initial
  use prognostics, only : prognostics_final
  use prognostics, only : dyn_prog

  use control,     only : cont2mod

  use m_model_tl,  only : amodel_tl

  use m_setsvecs

  use m_admtlm
! use m_admtlm,  only : admtlm_init
! use m_admtlm,  only : admtlm_set
! use m_admtlm,  only : admtlm_clean

  use m_sVecDef, only : eigchk
  use m_sVecDef, only : rvec
  use m_sVecDef, only : bmat
  use m_sVecDef, only : maxitr
  use m_sVecDef, only : ncalls
  use m_sVecDef, only : iparam
  use m_sVecDef, only : tol        ! Eigenvalue accuracy ARPACK
  use m_sVecDef, only : svecnormI  ! Initial state norm
  use m_sVecDef, only : svecnormF  ! Final state norm
  use m_sVecDef, only : propsvec   ! Evolve svecs or not?
  use m_sVecDef, only : isolve     ! Eigen solver package
  use m_sVecDef, only : kappa      ! Eigenvalue accuracy NAG
  use m_sVecDef, only : pflabel

  use m_svnorms, only : norm_zvecsize
  use m_svnorms, only : norm_svecA
  use m_svnorms, only : norm_svecBC

  use m_svprojs, only : proj_init
  use m_svprojs, only : proj_clean

  use m_delp2ps, only : delp2ps

  use m_iostate, only : PutState
  use m_iostate, only : getstate_init
  use m_iostate, only : getstate
  use m_iostate, only : getstate_clean

  use m_trajmng, only : PutPert

  use m_trjphys, only : physdrv1_get_init
  use m_trjphys, only : physdrv1_get_all
  use m_trjphys, only : physdrv1_get_clean

  use m_zeit, only : zeit_ci
  use m_zeit, only : zeit_co

  use m_ioutil, only : luavail

  use m_mpif90,only : MP_comm_rank
  use m_die, only: die, MP_die
  use m_stdio, only : stdout

  implicit none

! !INPUT PARAMETERS:

#include "lapack.H"
#include "arpack.H"
  integer, intent (in) :: ROOT 
  integer, intent (in) :: comm

! !OUTPUT PARAMETERS:

! !DESCRIPTION: 
!
! To do:
! \begin{itemize}
!    \item Output both state and pert as GFIO files (use dyn\_put)
!    \item Set controls for trajectory update
!    \item Implement norms 
!    \item Use m\_dyn way of defining single state vector x
!    \item Extend m\_dyn to define perturbations
! \end{itemize}
!
! !REVISION HISTORY:
!
!  25Oct2002 Todling - Initial code, half based on FastOpt demo, 
!  20Nov2000 Todling - Added opt to evolve svecs.
!  21Jan2003 Errico  - Corrected zvector size and use, 
!                      add TE and KE norms as only correct options  
!  11Mar2003 Todling - Redefined svec files naming convention
!  14Mar2003 Todling - Moved definition of test\_norm to rcfile level
!                      Updated interfaces to PutPert/PutState
!                      Redefined name of pert files for pesto purposes
!  26Apr2004 Todling - Moved admtlm to m_admtlm
!  25Jun2004 Gelaro  - Added option to use NAG eigensolver
!  25Apr2005 Gelaro  - Changed indexing to output SV info in descending order
!  23Jul2005 Gelaro  - Added output of svalues and convergence info
!  23Dec2005 Gelaro  - Added Lanczos restart capability
!  30Aug2006 Todling - Added CNOP interface and capability
!  07May2007 Todling - Added Kim's trajectory-in-mem knobs
!  13May2007 Todling - Introduced dyn_prog; global change
!  08Jun2007 Kim     - Various modifications for MPI implemetation
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname  = 'svdrv'
  integer,          parameter :: prec    = 0      ! 32-bits
  logical,          parameter :: verbose = .true.
  logical    :: memtraj = .true.

  integer    n           ! number of independent (initial state) variables
  integer    m           ! DUMMY 
  integer    nzvecsize   ! size of z vector 
  integer    nev, ncv, nconv
  integer    nymd_b, nhms_b, nymd_e, nhms_e
  integer    ido
  integer    lworkl
  integer    info, j
  integer    ierr
  integer    ipntr(11)
  integer    lanmax, lanstrt
  integer    lunit1
  integer    myID

  character(len=2)   :: which
  character(len=4)   :: what

  character(len=255) :: fname
  character(len=255) :: ftype
  character(len=255) :: svalu

  real       asigma       ! shifts for ARPACK (when applied)
  real(r8)   ynorm        ! value of norm at beginning or end 


! Allocatable arrays
! ------------------
  real(r8), allocatable :: ps(:,:)      ! surface pressure perturbation
  real(r8), allocatable :: Eval(:)      ! singular values
  real(r8), allocatable :: Lvec(:,:)    ! singular vectors
  real(r8), allocatable :: errbnd(:)    ! error bounds for singular values
  real(r8), allocatable :: resid(:)     ! all residuals
  real(r8), allocatable :: workd(:)     ! workspace for reverse communication
  real(r8), allocatable :: workl(:)     ! internal ARPACK workspace
  logical,  allocatable :: select(:)    ! specify how many/which Ritz 
                                        !   vectors to calculate
  type(dyn_prog) :: prog 
  type(dyn_prog) :: prog_tl
   
  call MP_comm_rank(comm,myID,ierr)
      if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

! Get number of the independent and dependent variables
! -----------------------------------------------------
  call setfunc( n, m, prog )
        if(n/=m) call die(myname,'invalid dimensions(n,m)')

! Sanity:
! ------
  if(n/=m) return  ! not allowed
  ierr = 0

! Initialize options, model, and set the first guess of the control variables
! ---------------------------------------------------------------------------
  what   = ' LSV'                               ! problem identifier 
  ncalls = 0                                    ! initialize lanczos counter
  call stepon_set ( prog )			! read in (initial nonlinear) state
  call set_svecs  ( comm, ROOT )         	! generic options for svec calc.
  if ( isolve .eq. 1 ) &                        ! ARPACK-specific options for SV calc.
       call set_svecs ( comm, ROOT, nev, ncv, which )
  if ( isolve .eq. 2 ) &                        ! NAG-specific options for SV calc.
       call set_svecs ( comm, ROOT, lanmax )
  if ( isolve .eq. 3 ) &                        ! SPG-based CNOP 
       call set_svecs ( comm, ROOT, what, ncv )

  write(fname,'(2a,i8.8,a)') trim(job), '.prog.eta.', nymd, '.nc4'
  if(fvpsasdt/=0) &
  call PutState ( fname,  &
                  job, nymd, nhms, prog, fvpsasdt, nstep, &
                  ak, bk, Ts, oro, &
                  prec=prec, stat=ierr ) ! output (initial) state GFIO file
	if(ierr/=0) call MP_die(myname,'Error from PutState',ierr)

  if ( isolve .ne. 3 ) then
       call norm_zvecsize ( comm, ROOT, svecnormI, nzvecsize, stat=ierr ) ! determine z vector size 
	  if(ierr/=0) call MP_die(myname,'Error from norm_zvecsize',ierr)
  else
       nzvecsize = n  ! needs revision when SPMD
  endif

! Read starting Lanczos iteration number from svalu file
! ------------------------------------------------------
  write( svalu, '(2a,i8.8,a)' ) trim(job), '.svalu.', nymd, '.txt'
  lunit1= luavail()
  open( unit=lunit1, file=svalu, form='formatted' )
  read( lunit1, '(a49,i3)', iostat=ierr ) fname, lanstrt 
      if(ierr/=0) lanstrt = 0
  close( lunit1 )
  if( isolve == 2 ) ncv= min( lanmax, lanstrt+maxitr )
  if( isolve <= 2 ) then
     if(myID==ROOT) write(stdout,'(1x,2a,i5)') myname, ': Starting iteration for Lanczos= ', lanstrt
  endif

! Initialize norm/projector
! -------------------------
  call proj_init ( comm, root, ierr )
     if(ierr/=0) call MP_die ( myname, 'Error from Proj_Init:', ierr )

! Allocate workspace for eigenpairs
! ---------------------------------
  allocate ( Eval(ncv), Lvec(nzvecsize,ncv), errbnd(ncv), stat=ierr )
        if(ierr/=0) call MP_die(myname,'Alloc(Eval,Lvec,errbnd)',ierr)

! Initialize dynamics trajectory handle
! -------------------------------------
  call getstate_init ( nymd, nhms, memtrj=memtraj, verbose=verbose )
  call getstate      ( nymd, nhms )

! Initialize physics trajectory handle
! -------------------------------------
  call physdrv1_get_init ( nymd, nhms, memphys=memtraj, verbose=verbose )
  call physdrv1_get_all  ( nymd, nhms )

! Choose NAG or ARPACK routines for Lanczos algorithm
! ---------------------------------------------------
  if ( isolve == 1 ) call lanc_ARPACK
  if ( isolve == 2 ) call lanc_PNAG
  if ( isolve == 3 ) call cnop ( comm, ROOT )

  
! Postprocessing
! --------------

! Write eigenvalues & convergence info to svalu file
! --------------------------------------------------
  if (isolve==2) ncalls = ncalls+1  ! Account for 'seed' iteration in NAG
  lanstrt = min( ncalls, lanstrt+maxitr )
  lunit1= luavail()
  if ( myID==ROOT ) then

    open( unit=lunit1, file=svalu, form='formatted', status='replace' )
    write( lunit1, '(a,i8.8,a,i6.6,a,i3)' )  &
       ' num      eigval     errbnd      ', nymd, '_', nhms, ' ', lanstrt

    write(stdout,'(1x,2a,i5)') myname, 'Ending iteration for Lanczos= ', lanstrt

    do j = 1, nconv  ! loop over converged eigenvalues  ! flip index for descending order
       write(lunit1, '(i4,2e14.6)') j, Eval(nconv-j+1), errbnd(nconv-j+1) 
    end do

  endif ! < ROOT >

  
! Eigenvectors 
! ------------
  if ( isolve == 1 ) rvec = rvec .and. ido==99         ! Postprocess SVs only
  if ( isolve == 2 ) rvec = rvec .and. ncalls>=lanmax  ! after final segment

  if ( rvec ) then

!   Allocate perturbation arrays 
!   ----------------------------
    call prognostics_initial ( prog_tl )

    do j = nconv, 1, -1  ! loop over converged SVs; reverse index for descending order
     
!      Unscale and remap singular vector z to gridded field array xpert
!      ----------------------------------------------------------------
       if ( isolve .ne. 3 ) then
          call norm_svecA  ( comm, ROOT, svecnormI, prog, prog_tl, nzvecsize, Lvec(:,j), ierr )
            if (ierr/=0) call MP_die(myname,'Error from norm_svecA',ierr)
!_RT   else
!_RT       xpert(:) = Lvec(:,j)
       endif

!      Write out singular vector
!      -------------------------
       allocate (ps(imr,jfirst:jlast))
       call delp2ps ( prog_tl%delp, ps )
       write(fname,'(4a)')     trim(job), '.i', trim(pflabel), '.eta'
       write(ftype,'(a,i3.3)') trim(svecnormI), nconv-j+1

       if(fvpsasdt/=0) &
       call PutPert ( job, nymd, nhms, prog_tl, fvpsasdt, nstep, &
                      ak, bk, Ts, oro, ps, fname, ftype )
       deallocate (ps)

!      Check normalization of SVs (as checked here, beginning time and end
!      time norms must be the same).
!      -------------------------------------------------------------------
       call norm_svecBC ( comm, ROOT, svecnormI, prog, prog_tl, ynorm, 'NA', ierr )
       if(ierr/=0) call MP_die(myname,'Error from norm_svecBC C2',ierr)
       if(myID==ROOT) write(stdout,'(1x,2a,1p,e13.6,3a,i5)') myname,' initial ynorm= ',ynorm,' for ', what,'-',nconv-j+1


    end do
    call prognostics_final ( prog_tl )

!   If so, propagate singular vectors
!   NOTE: need to do separate appl 
!         for this part
!   ---------------------------------
    if ( propsvec ) then

    call zeit_ci('propv')
    call prognostics_initial ( prog_tl )

       do j = nconv, 1, -1  ! loop over evolved SVs; reverse index for descending order

!         Unscale and remap singular vector z to gridded field array xpert
!         ----------------------------------------------------------------
          if ( isolve .ne. 3 ) then
             call norm_svecA  ( comm, ROOT, svecnormI, prog, prog_tl, nzvecsize, Lvec(:,j), ierr )
               if (ierr/=0) call MP_die(myname,'Error from norm_svecA (2)',ierr)
!_RT      else
!_RT         xpert(:) = Lvec(:,j)
          endif

!         Apply TLM
!         ---------
          call amodel_tl ( prog_tl )

!         Write out evolved vector(s)
!         ---------------------------
          allocate (ps(imr,jfirst:jlast))
          call delp2ps ( prog_tl%delp, ps )
          write(fname,'(4a)')     trim(job), '.f', trim(pflabel), '.eta'
          write(ftype,'(a,i3.3)') trim(svecnormI), nconv-j+1

	  if(fvpsasdt/=0) &
          call PutPert ( job, nymd, nhms, prog_tl, fvpsasdt, nstep, &
                         ak, bk, Ts, oro, ps, fname, ftype )
          deallocate (ps)

!         Evaluate norm of y and write it out
!         -----------------------------------
          call norm_svecBC ( comm, ROOT, svecnormF, prog, prog_tl, ynorm, 'NA', ierr )
	       if(ierr/=0) call MP_die(myname,'Error from norm_svecBC',ierr)
          write(stdout,'(1x,2a,1p,e13.6,3a,i5)') myname,' final ynorm=',ynorm,' for ', what,'-',nconv-j+1

       end do

       call prognostics_final ( prog_tl )
       call zeit_ci('propv')

    endif  ! < propsvec >

  endif  ! < rvec >

  call prognostics_final ( prog )

! Clean norm/projector
! --------------------
  call proj_clean ( ierr )
        if(ierr/=0) call MP_die(myname,'Error from Proj_Clean',ierr)

  deallocate ( Eval, Lvec, errbnd, stat=ierr )
	if(ierr/=0) call MP_die(myname,'Dealloc(Eval,Lvec,errbnd)',ierr)

  call getstate_clean
  call physdrv1_get_clean

!_RT  deallocate ( x, y, stat=ierr )

  return

CONTAINS
subroutine lanc_ARPACK

!---------------------------------------------------------------------------
! Routine for ARPACK-based Lanczos algorithm
!---------------------------------------------------------------------------
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  25Oct2002 Todling - Initial code.
!  28Sep2007 Todling - Add timers
!
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::lanc_ARPACK'
!
! Allocate workspace for eigen-decomposition
! ------------------------------------------
  lworkl = ncv*(ncv+8)
  allocate ( resid(nzvecsize), workd(3*nzvecsize), workl(lworkl), stat=ierr )
        if(ierr/=0) call MP_die(myname_,'Alloc(resid,etc)',ierr)
                                                                                           
  call zeit_ci ('lanc')
  call admtlm_init ( n, m, nzvecsize )
  call admtlm_set  ( prog )

! M A I N   L O O P
! -----------------
                                                                                           
  ido  = 0
  info = 0
  if(myid==root) write(stdout,'(2a)') myname_, ': Entering SV/ARPACK Loop ...'
                                                                                           
  do while ( (ido==0 .or. ido==-1 .or. ido==1) .and. (ncalls<maxitr+lanstrt) )
                                                                                           
!    %----------------------------------------------%
!    | Repeatedly call the routine _SAUPD and take  |
!    | actions indicated by parameter IDO until     |
!    | either convergence is achieved or maximum    |
!    | number of iterations has been exceeded.      |
!    %----------------------------------------------%
                                                                                           
#if defined ( SPMD )
        call zeit_ci ( '_PSAUPD' )
     call _PSAUPD ( comm, ido, bmat, nzvecsize, which, nev, tol, resid, &
                     ncv, Lvec, nzvecsize, iparam, ipntr,          &
                     workd, workl, lworkl, info )
        call zeit_co( '_PSAUPD' )
#else
        call zeit_ci ( '_SAUPD' )
      call _SAUPD ( ido, bmat, nzvecsize, which, nev, tol, resid, &
                    ncv, Lvec, nzvecsize, iparam, ipntr,          &
                    workd, workl, lworkl, info )
        call zeit_co ( '_SAUPD' )
#endif
                     
!    %------------------------------------------------%
!    | Perform matrix vector multiplication           |
!    |         z_output = S^T * S * z_input           |
!    |         S= N_2^(1/2) * TLM * N_1^(-1/2)        |
!    | where z_input=workd(ipntr(1))                  |
!    | where z_output=workd(ipntr(2))                 |
!    %------------------------------------------------%
                     
     if ( ido /= 99 ) then  ! convergence flag from _SAUPD
                     
        workd(ipntr(2):ipntr(2)+nzvecsize-1) = workd(ipntr(1):ipntr(1)+nzvecsize-1)
        call admtlm ( workd(ipntr(2):ipntr(2)+nzvecsize-1), nzvecsize )
        
     end if
 
  end do
 
! %--------------------------------------------------%
! | Either we have convergence or there is an error. |
! %--------------------------------------------------%
!
  if ( info<0 ) then
 
!    %---------------------------------------------------%
!    | Error message. Check the documentation in SSAUPD. |
!    %---------------------------------------------------%
 
     call MP_die ( myname_, 'Error with _saupd, info = ', info )
 
  else
 
!    %-------------------------------------------%
!    | No fatal errors occurred.                 |
!    | Post-Process using SSEUPD.                |
!    |                                           |
!    | Eigenvalues can be retrieved.             |
!    |                                           |
!    | Eigenvectors may also be calculated if    |
!    | rvec = .true.                             |
!    |                                           |
!    %-------------------------------------------%

     nconv =  iparam(5)         !  This indicates how many evals are
                                !  accurate to the requested tolerance

! Make eval convergence info (error bounds) available for output
! --------------------------------------------------------------
     errbnd(1:nconv) = workl(3*ncv+1:3*ncv+nconv)
 
     allocate ( select(ncv), stat=ierr )
        if(ierr/=0) call MP_die(myname_,'Alloc(select)',ierr)
        
#if defined ( SPMD )
        call zeit_ci ( '_PSEUPD' )
     call _PSEUPD ( comm, rvec, 'A', select, Eval, Lvec, nzvecsize, asigma, &
                   bmat, nzvecsize, which, nev, tol, resid, ncv, Lvec, &
                   nzvecsize, iparam, ipntr, workd, workl, lworkl, ierr )
        call zeit_co ( '_PSEUPD' )
#else
        call zeit_ci ( '_SEUPD' )
     call _SEUPD ( rvec, 'A', select, Eval, Lvec, nzvecsize, asigma, &
                   bmat, nzvecsize, which, nev, tol, resid, ncv, Lvec, &
                   nzvecsize, iparam, ipntr, workd, workl, lworkl, ierr )
        call zeit_co ( '_SEUPD' )
#endif
        
!    %----------------------------------------------%
!    | Eigenvalues are returned in Eval and their   |
!    | corresponding eigenvectors are returned in   |
!    | the first NCONV (=IPARAM(5)) columns of the  |
!    | two dimensional array Lvec if requested.     |
!    | Otherwise, an orthogonal basis for the       |
!    | invariant subspace corresponding to the      |
!    | eigenvalues in Eval is returned in Lvec.     |
!    %----------------------------------------------%
 
     if ( ierr/=0 ) then
      
!       %-----------------------------------------------------%
!       | Error condition: Check the documentation of SSEUPD. |
!       %-----------------------------------------------------%
      
        call MP_die(myname_,'Error in _seupd, info',ierr)
      
     else if ( eigchk ) then
      
        call chk_svec ( root, comm, nzvecsize, ncv, nconv, Eval, Lvec )
      
     end if ! < info,eigchk >
      
!    %-------------------------------------------%
!    | Print additional convergence information. |
!    %-------------------------------------------%
      
     call sumry_ARPACK ( comm, ROOT, info, nzvecsize, nev, ncv, nconv, which, Eval )
      
!    Release space
!    -------------
     deallocate ( select )
      
  end if ! < info/=0 >

  call admtlm_clean()
  call zeit_co ('lanc')
      
  deallocate ( resid, workd, workl, stat=ierr )
        if(ierr/=0) call MP_die(myname_,'Dealloc(resid,etc)',ierr)

end subroutine lanc_ARPACK

subroutine lanc_PNAG

!---------------------------------------------------------------------------
! Routine for Parallel NAG-based Lanczos algorithm
!---------------------------------------------------------------------------

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  25Jun2004 Gelaro - Initial code.
!  23Jul2005 Gelaro - Made svalues and convergence info available for output
!  23Dec2005 Gelaro - Modified for Lanczos restart capability
!  09Jul2007 Kim    - Modified for Parallel NAG-based Lanczos algorithm
!  29Aug2007 Kim    - Got rid of evecs module; updated interface to clanc
!  28Sep2007 Todling- Revised timers
!
!-------------------------------------------------------------------------

  use m_ioutil, only : luavail

  character(len=*), parameter :: myname_ = myname//'::lanc_NAG'

  integer    nw
  integer    lsteps
  integer    lunit2
  integer    i, J, N_lanso, MAXPRS, CONDM, ENDL, ENDR, lunit, MSGLVL

  integer    indev(ncv)
  integer    indvec(ncv)

  real(r8)   bndev(ncv)

!
! Allocate workspace for eigen-decomposition
! ------------------------------------------
  nw = 6*nzvecsize+1+4*ncv+ncv*ncv
  allocate ( workd(nw), stat=ierr )
        if(ierr/=0) call MP_die(myname_,'Alloc(workd)',ierr)


! Initialize internal state vector and get basic state
! ----------------------------------------------------
  call zeit_ci ('lanc')
  call admtlm_init ( n, m, nzvecsize )
  call admtlm_set  ( prog )

! Set flag to generate eigenvectors
! ---------------------------------
  lunit2= -1
  if ( rvec ) lunit2= luavail()
  if(myid==root) write(stdout,'(1x,2a,i5)') myname_,': NAG eigenvector indicator unit= ',lunit2

  write(stdout,'(1x,2a,i5)') myname_,': PNAG eigenvector indicator unit= ',lunit2

! Primary interface for Parallel NAG Lanczos package
! --------------------------------------------------
  call zeit_ci ( 'pLANSO' )

  CALL clanc(nzvecsize,ncv,lunit2,KAPPA, &
             lsteps,nev,Lvec,Eval,bndev,workd,nw,nconv,indev,indvec,ierr,comm)
             ! lsteps number of Lanczos steps actually taken,
             ! Nev    number of Ritz values stabilized(see KAPPA),
             ! Eval   array that holds the Ritz values,
             ! BNDEV  array that holds the error bounds

  call zeit_ci ( 'pLANSO' )

! Print convergence information
! -----------------------------
  call sumry_NAG ( comm, ROOT, ierr, nzvecsize, ncv, lunit2, lsteps, nev, &
                   nconv, Eval, bndev, indev, indvec )

! Make converged eigenvectors, etc. accessible for postprocessing
! NOTE: Lvec's already come out in the proper order
! ---------------------------------------------------------------
  do i = 1, ncv
     if (workd(i+ncv) .le. kappa*abs(workd(i))) then
         Eval(i)  = workd(i)
         errbnd(i)= workd(i+ncv)
     endif
  enddo

! Release workspace
! -----------------
  deallocate ( workd, stat=ierr )
        if(ierr/=0) call MP_die(myname_,'Dealloc(workd)',ierr)
  call zeit_co ('lanc')
                                                                                 
end subroutine lanc_PNAG

subroutine cnop ( comm, ROOT )

!---------------------------------------------------------------------------
! Routine to calculate conditional nonlinear optimal perturbations
!---------------------------------------------------------------------------
!
! !USES:
!
  use m_sVecDef, only : spg_miter
  use m_sVecDef, only : spg_mfls
  use m_sVecDef, only : spg_einf
  use m_sVecDef, only : cnop_tol
  use m_sVecDef, only : spg_maxfc
  use m_sVecDef, only : spg_verb
  use m_sVecDef, only : cnop_sigma

  use m_cnop, only : cnop_init
  use m_cnop, only : cnop_set
  use m_cnop, only : cnop_clean
  use m_cnop, only : cnop_fcnt
  
! !INPUT PARAMETERS: 

  integer, intent(in) :: comm
  integer, intent(in) :: ROOT
   
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  14Aug2006 Todling - Initial code.
!
!-------------------------------------------------------------------------
 
  character(len=*), parameter :: myname_ = myname//'::cnop'

  integer  nfevals  ! number of function evaluations
  integer  ngevals  ! number of gradient evaluations
  integer  niter    ! number of iterations
  real(r8) fval     ! value of function at the minimum
  real(r8) cginfn   !
  real(r8) cgtwon   !
  
! Initialization and setup
! ------------------------
  call cnop_init ( n, ierr )
     if(ierr/=0) call MP_die(myname_,'Failed in cnop_init ',ierr)
  if(n/=nzvecsize) call MP_die(myname_,'Dimension inconsistency in cnop ',ierr)
!_RT  call cnop_set  ( comm, ROOT, n, x, ierr, xpert=Lvec(:,1) )
!_RT     if(ierr/=0) call MP_die(myname_,'Failed in cnop_set ',ierr)
  
! Calculate optimal perturbation
! ------------------------------
  call dspg ( n, Lvec(:,1), spg_mfls, spg_einf, cnop_tol, & 
              spg_miter, spg_maxfc, spg_verb, &
              fval, cginfn, cgtwon, niter, nfevals, ngevals, ierr )

     if ( ierr <= 1 ) then
          nconv = 1
          Eval(1) = cnop_sigma  ! for now to fill up the whole
          if(myid==root) write(stdout,'(2a)') trim(myname_), ' successfully calculated CNOP'
     else if ( ierr == 2 ) then
          nconv = 1
          Eval(1) = cnop_sigma  ! for now to fill up the whole
          if(myid==root) write(stdout,'(2a)') '--------------------------------------------'
          if(myid==root) write(stdout,'(2a)') '         CAUTION  CAUTION  CAUTION          '
          if(myid==root) write(stdout,'(2a)') trim(myname_), ' CNOP has not fully converged'
          if(myid==root) write(stdout,'(2a)') '--------------------------------------------'
     else
          nconv = -1
          if(myid==root) write(stdout,'(2a,i3)') trim(myname_), ' was unsuccessful in finding CNOP, ierr = ', ierr
     endif
     if(myid==root) write(stdout,'(2a,i7)') trim(myname_), ' Number of function calls (no grad): ', cnop_fcnt
  
! Wrap up
! -------
  call cnop_clean ()
   
end subroutine cnop

end subroutine svdrv
!EOC

!--------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!--------------------------------------------------------------------------
!BOP
!
! !ROUTINE: chk_svec --- Check eigenvectors/singular vectors
!
! !INTERFACE:

subroutine chk_svec ( root, comm, nzvecsize, ncv, nconv, Eval, Lvec ) 

! !USES:

  use prognostics, only : dyn_prog
  use prognostics, only : prognostics_initial
  use prognostics, only : prognostics_final
  use m_admtlm
  use m_die,   only: MP_die
  use m_stdio, only: stdout

  implicit none

! !INPUT PARAMETERS:

  integer, intent(in) :: nzvecsize
  integer, intent(in) :: root
  integer, intent(in) :: comm
  integer, intent(in) :: ncv
  integer, intent(in) :: nconv
  real,    intent(in) :: Eval(ncv)
  real,    intent(in) :: Lvec(nzvecsize,ncv)

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  25Oct2002 Todling - Initial code.
!  18May2007 Todling - Reimplemented call to operator.
!
!EOP
!-------------------------------------------------------------------------
  character(len=*), parameter :: myname = 'chk_svec'

  integer  ierr, j
#if defined ( SPMD )
  real     _PSNRM2 
#else
  real     _SNRM2
#endif
  real,    allocatable :: ax(:)           ! array use for residual check
  real,    allocatable :: Xval(:)         ! computed eval residuals
  integer  n, m
  type(dyn_prog)       :: prog

  call admtlm_set  ( prog )

!       Allocate workspace
!       ------------------
        allocate ( ax(nzvecsize), Xval(nconv), stat=ierr )
            if(ierr/=0) call MP_die(myname,'Alloc(ax,xval)',ierr)

        do j = 1, nconv

!          %--------------------------------------%
!          | Compute the residual norm            |
!          |                                      |
!          |   ||  Op*x - lambda*x ||             |
!          |                                      |
!          | for the NCONV accurately computed    |
!          | eigenvalues and eigenvectors.        |
!          %--------------------------------------%

                         ax = Lvec(1:nzvecsize,j)
           call  admtlm ( ax, nzvecsize )
           call _SAXPY (nzvecsize, -Eval(j), Lvec(1,j), 1, ax, 1)
#if defined ( SPMD )
                Xval(j) = _PSNRM2(comm, nzvecsize, ax, 1)
#else
                Xval(j) = _SNRM2(nzvecsize, ax, 1)
#endif
                Xval(j) = Xval(j) / abs(Eval(j))

        end do  ! < j >

!       %-----------------------------%
!       | Display computed residuals. |
!       %-----------------------------%
!
        call _MOUT(stdout, nconv, 1, Eval, ncv, -6, 'Ritz values')
        call _MOUT(stdout, nconv, 1, Xval, ncv, -6, 'Computed residuals')

!       Release memory
!       --------------
        deallocate ( ax, Xval )

end subroutine chk_svec


!--------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!--------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sumry_ARPACK --- Summarizes ARPACK calculations
!
! !INTERFACE:

subroutine sumry_ARPACK ( comm, ROOT, info, nzvecsize, nev, ncv, nconv, which, Eval )

! !USES:

  use m_svecdef
  use m_stdio, only: stdout
  use m_mpif90,only : MP_comm_rank
  use m_die, only : MP_die

  implicit none

! !INPUT PARAMETERS:

  integer, intent(in) :: comm
  integer, intent(in) :: ROOT
  integer, intent(in) :: info
  integer, intent(in) :: nzvecsize
  integer, intent(in) :: nev
  integer, intent(in) :: ncv
  integer, intent(in) :: nconv
  real   , intent(in) :: Eval(ncv)

  character(len=2), intent(in) :: which

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  25Oct2002 FastOpt - Initial code.
!  15Jun2007 Gelaro - Root write protection.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'sumry_ARPACK'
  integer  j,myID,ier

  call MP_comm_rank(comm,myID,ier)
      if(ier/=0) call MP_die(myname,'MP_comm_rank()',ier)

  if ( info==1 ) then
       if(myID==ROOT) write(stdout,'(2a)') myname, &
                    ': Maximum number of iterations reached.'
  else if ( info .eq. 3 ) then
       if(myID==ROOT) write(stdout,'(3a)') myname, &
                    ' No shifts could be applied during implicit', &
                    ' Arnoldi update, try increasing NCV.'
  end if  ! < info==1 >

  if ( myID==ROOT ) then
       write(stdout,'(2a,/)') myname,  &
                    ': Summary of Eigen-Decomposition'
       write(stdout,'(a,i4,a,i7)') ' On PE: ', myID,' Size of the matrix is ', nzvecsize
       write(stdout,'(a,i7)') ' The number of Ritz values requested is ', nev
       write(stdout,'(2a,i7)') &
                   ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
       write(stdout,'(2a)') ' What portion of the spectrum: ', which
       write(stdout,'(a,i7)') ' The number of converged Ritz values is ', nconv
       write(stdout,'(2a,i7)') ' The number of Implicit Arnoldi update', &
                   ' iterations taken is ', iparam(3)
       write(stdout,'(a,i7)') ' The number of OP*x is ', iparam(9)
       write(stdout,'(a,e9.2)') ' The convergence criterion is ', tol
    
       write(stdout,'(2a)') myname, ': converged eigenvalues follow'
       write(stdout,'(1p,7e11.4)') (Eval(j),j=1,nconv)
       call flush(stdout)
  endif

end subroutine sumry_ARPACK



!--------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!--------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sumry_NAG --- Summarizes NAG calculations
!
! !INTERFACE:

subroutine sumry_NAG ( comm, ROOT, ierr, nzvecsize, lanmax, lunit2, lsteps, nev, &
                       nconv, Eval, bndev, indev, indvec )

! !USES:

  use m_svecdef,  only : kappa
  use m_stdio,    only : stdout
  use m_mpif90,only : MP_comm_rank
  use m_die, only : MP_die

  implicit none

! !INPUT PARAMETERS:

  integer, intent(in) :: comm
  integer, intent(in) :: ROOT
  integer, intent(in) :: ierr
  integer, intent(in) :: nzvecsize
  integer, intent(in) :: lanmax
  integer, intent(in) :: lunit2
  integer, intent(in) :: lsteps
  integer, intent(in) :: nev
  integer, intent(in) :: nconv
  integer, intent(in) :: indev(lanmax)
  integer, intent(in) :: indvec(lanmax)

  real   , intent(in) :: Eval(nconv)
  real   , intent(in) :: bndev(nconv)

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  25Jun2004 Gelaro - Initial code.
!  15Jun2007 Gelaro - Root write protection.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'sumry_NAG'
  integer  j,myID,ier

  call MP_comm_rank(comm,myID,ier)
      if(ier/=0) call MP_die(myname,'MP_comm_rank()',ier)

  if ( myID==ROOT ) then

     write(stdout,'(2a,/)') myname,  &
                  ': Summary of Eigen-Decomposition'

     write (stdout,'(a,i10)')    'Dimension of problem  (nzvecsize) = ', nzvecsize
     write (stdout,'(a,i10)')    'Max number iterations (lanmax)    = ', lanmax
     write (stdout,'(a,i10)')    'Lanczos output unit   (lunit2)    = ', lunit2
     write (stdout,'(a,e22.15)') 'Accuracy of e-vectors (kappa)     = ', kappa
     write (stdout,'(a,i10)')    'Lanczos steps taken   (lsteps)    = ', lsteps
     write (stdout,'(a,i10)')    'E-values accepted     (nev)       = ', nev
     write (stdout,'(a,i10)')    'E-vectors accepted    (nconv)     = ', nconv
     write (stdout,'(a,i10)')    'Error condition       (ierr)      = ', ierr
!
     write (stdout,201) ( Eval   (j) , j = 1 , lanmax )
     write (stdout,202) ( bndev  (j) , j = 1 , lanmax )
  201     format (/2x,'e-valus',5e15.7/40(9x,5e15.7/))
  202     format (/2x,'e-bnds ',5e15.7/40(9x,5e15.7/))
!
     write (stdout,204) ( indev   (j) , j = 1 , lanmax )
     write (stdout,205) ( indvec  (j) , j = 1 , lanmax )
  204     format (/2x,'indev  ',5i15/40(9x,5i15/))
  205     format (/2x,'indvec ',5i15/40(9x,5i15/))

     call flush(stdout)
  endif ! < ROOT >

end subroutine sumry_NAG




