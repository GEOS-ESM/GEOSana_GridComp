!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_physdrv1_ad --- Adjoint of simplified physics
!
! !INTERFACE:
!
module m_physdrv1_ad

! !USES:

      use precision 
      use m_const, only : zvir
      use M_GFIO_GETFLD
      use m_StrTemplate, only : strTemplate
      use stepon, only : nsecf
      use stepon, only : tdt
      use m_die, only : die
      use m_trjphys, only: physdrv1_get_state
      use m_trjphys, only: nvars

#ifdef SPMD
      use m_mpif, only:mpi_double_precision,mpi_comm_world
      use mod_comm, only : mp_send_s, mp_recv_n,             &
                           mp_send3d_ns, mp_recv3d_ns,        &
                           mp_barrier,mp_send_s_ad,mp_recv_n_ad, &
                           mp_send3d_ns_ad,mp_recv3d_ns_ad,gid, &
                           mp_gather4d,mp_scatter4d,mp_send_s_ad,mp_recv_n_ad
#endif

  implicit none
! !PUBLIC MEMBER FUNCTIONS:
 
  PRIVATE

  PUBLIC  physdrv1_ad

  interface physdrv1_ad ; module procedure physdrv1_ad_ ; end interface

!
! !DESCRIPTION: This module implements the first version of the
!               adjoint simplified physics (vertical diffusion).
!
! !REVISION HISTORY:
!
!  13Feb2004  Todling   Modularized package.
!  28Jul2004  Errico/RT Scaling coefficients w/ dt.
!  08Jul2005  Errico/EN Changed calls to horiz diff routines
!  15Jul2005  Elena N.  Modified USES for this entire module
!  19Apr2006  Errico/Elena N. Adjusted sponge parameters for GEOS-5
!  22Nov2006  Oloso     Various SPMD-related fixes throughout
!  10Jan2007  Todling   Cleaned lines left commented by Oloso
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'm_physdrv1_ad'
   integer,          parameter :: ROOT = 0

!
!  WARNING **** The horiz diff  code will only work correctly for jfirst=1
!  jlast=jnp; The problem is that extra latitudes are required to
!  compute the gradients at the 2 edges of each sub domain.  An MPI version
!  will therefore require modification.
!
!  Also Note that the parameters in the hdiff_tl routines must be set the same
!  as in the hdiff_tl routines until the code is restructured.

       real(8), parameter :: earthr=6.37d6
       real(8), parameter :: fcrit=1.0d0/8.0d0
       real(8), parameter :: tdamp_days=2.0d0 
       real(8), parameter :: tdamp_dt=4.d0 ! limiting damping time in *dt
       real(8), parameter :: pi =  3.1415926535897931
       real(8), parameter :: twopi = 2.d0*pi


CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: physdrv1_ad_ --- adjoint of version 1 of driver for simplified physics
!
! !INTERFACE:

subroutine physdrv1_ad_ ( ptrjtmpl, ptrjfrq, job, imr, jnp, nl, nc, jfirst, jlast, &
                          coslon, sinlon,  q, nymd, nhms,  &
                          u_ad,   v_ad,   pt_ad, ng_d, ng_s)

! !USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in)  :: ptrjtmpl   ! filename template for vdiff fields
      integer, intent(in)  :: ptrjfrq             ! frequecy (HHMMSS) of trajectory update
                                                  !   < 0 returns without exercising routine
                                                  !   = 0 doesn't update trajectory
                                                  !   > 0 updates trajectory
      character(len=*), intent(in)  :: job        ! name of experiment. DO THIS SOME OTHER WAY!!!
      integer, intent(in)  :: imr                 ! number of grid points in long. direction
      integer, intent(in)  :: jnp                 ! number of grid points in lat.  direction
      integer, intent(in)  :: nl                  ! number of grid points in vert. direction
      integer, intent(in)  :: nc                  ! number of constituent fields  
      integer, intent(in)  :: jfirst              !
      integer, intent(in)  :: jlast               !
      integer, intent(in)  :: ng_d              !
      integer, intent(in)  :: ng_s              !
      real(kind=r8), intent(in) :: coslon(imr)         ! cosine of longitude
      real(kind=r8), intent(in) :: sinlon(imr)         ! sine   of longitude
      real(kind=r8), intent(in) :: q(imr,jfirst-ng_d:jlast+ng_d,nl,nc)  

      integer   nymd                              !  Year-month-day, e.g.,  19971012
      integer   nhms                              !  hour-minute-sec, e.g.,   120000

! !INPUT/OUTPUT PARAMETERS:

      real(kind=r8), intent(inout) :: u_ad(imr,jfirst-ng_d:jlast+ng_s,nl)
      real(kind=r8), intent(inout) :: v_ad(imr,jfirst-ng_s:jlast+ng_d,nl)
      real(kind=r8), intent(inout) :: pt_ad(imr,jfirst-ng_d:jlast+ng_d,nl)   ! virt pt


! !DESCRIPTION: This is the adjoint of subroutine physdrv1.  Note that 
!               the coding does not follow the TAF standard.  The location
!               of some adjoint calculations have been shifted to increase 
!               efficiency. 
!
! !REVISION HISTORY:
!
!  06Jan2004  Errico  Initial code and implementation
!  07Jan2004  Winslow modified variable names, and placed into NLM/ADM context
!  10Feb2004  Todling Added template ref for i/o traj filename;
!                     added feature to allow skipping physics (simply don't
!                     provide a phys. trajectory file)
!  13Feb2004  Todling Set explicitly  ng_d, ng_s to zero when calling 
!                     d2a3d_ad since this code only works in the non-SPMD
!                     case - need revision otherwise (see cal
!  10Oct2005  Errico  corrections added to Hdiff
!  14Jan2007  Todling Bug fix, after Oloso's change (utotal wasn't def in non SPMD case)
!  01Jul2007  Todling Simplified: 2nd dim of vcoef need only jf to jl
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname_ = myname//':physdrv1_ad_'

! !LOCAL VARIALBLES
!
! The following are lower (a) and upper (c) diagonals (for a matrix ordered top! to bottom)
! of tri-diagonal matrix of the vertical diffusion, pbl, and
! surface drag equations for momentum (m) and heat (h):

      character*255 :: vFile                      ! name of file containing coefs
      real(kind=r8), allocatable :: vcoef3d(:,:,:,:)  ! upper and lower diagonal elements of matrix
                                                  ! used in vdiff (1=cah,2=cam,3=cch,4=ccm)
      integer, parameter :: nfields2=2 ! max number of simultaneous field slices 
      integer, parameter :: ntopfl=1   ! topmost level for vertical diffusion
                                       ! copied from vdinti.F
      integer            :: i          ! longitudinal index
      integer            :: j          ! latitudinal index
      integer            :: k          ! vertical index
      integer            :: nfields    ! number of simultaneous field slices 
      integer            :: ierr, iouter_phys
!
      real(kind=r8) :: ua_ad(imr,jfirst-ng_d:jlast,nl) ! u or its physics change at pt points
      real(kind=r8) :: va_ad(imr,jfirst:jlast,nl) ! v or its physics change at pt points
      real(kind=r8) :: ta_ad(imr,jfirst:jlast,nl) ! physics change of pt 
      real(kind=r8) :: fm_ad(imr,nl,nfields2) ! 2-d slice of fields before physics
      real(kind=r8) :: fp_ad(imr,nl,nfields2) ! 2-d slice of fields after physics
      real(kind=r8) :: ca(imr,nl)      ! upper diagonal of matrix (for 1 lat) 
      real(kind=r8) :: cc(imr,nl)      ! lower diagonal of matrix (for 1 lat)
      real(r8) :: utotal(imr,jnp,nl)

!     Return when no simplified physics not wanted
!     --------------------------------------------
      if (ptrjfrq<0) return

      allocate ( vcoef3d(imr, jfirst:jlast, nl, nvars), stat=ierr )
          if(ierr/=0) call die(myname,' Alloc(vcoef3d) error',ierr)

!     Get the diffusion coefficents
!     -----------------------------
      if ( ptrjfrq == 0 ) then
         if(gid==ROOT) print *, myname_, ': this option not yet implemented '
         call exit(7)
      else
         call physdrv1_get_state(vcoef3d,nymd,nhms) 
      endif

!
! adjoint of horizontal diffusion

    call hdiff_ad (imr,jnp,nl,jfirst,jlast,tdt,pt_ad,u_ad,v_ad,ng_d,ng_s)
!
!  Apply sponge in top portion of model
     call sponge_tl_ ( imr, jnp, jfirst, jlast, nl, tdt, pt_ad(:,jfirst:jlast,:), 'pt' )
     call sponge_tl_ ( imr, jnp, jfirst, jlast, nl, tdt, u_ad(:,jfirst:jlast,:),  'u'  )
     call sponge_tl_ ( imr, jnp, jfirst, jlast, nl, tdt, v_ad(:,jfirst:jlast,:),  'v'  )

! Adjoint of remapping of u changes back to D grid
! ------------------------------------------------
#ifdef SPMD
    call mp_gather4d(u_ad,utotal,imr,jnp,nl,1,jfirst,jlast, &
                    1,nl,ng_d,ng_d,0)
    call mpi_bcast(utotal,imr*jnp*nl,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
#else
    utotal = u_ad
#endif
      do j=max(2,jfirst),min(jnp-1,jlast)
        ua_ad(1:imr,j,ntopfl:nl)=0.5d0*(utotal(1:imr,j,ntopfl:nl)+utotal(1:imr,j+1,ntopfl:nl))
      enddo
      if (jfirst == 1) then
        ua_ad(1:imr,1,ntopfl:nl)=0.5d0*u_ad(1:imr,2,ntopfl:nl)
      endif
      if (jlast == jnp) then 
        ua_ad(1:imr,jnp,ntopfl:nl)=0.5d0*u_ad(1:imr,jnp,ntopfl:nl)
      endif
 
! Adjoint of remapping of v changes back to D grid
! ------------------------------------------------
      va_ad(:,jfirst,:) = 0.0d0
      va_ad(:,jlast ,:) = 0.0d0
      do j=max(2,jfirst),min(jnp-1,jlast)
        va_ad(imr,j,ntopfl:nl)=0.5d0*(v_ad(imr,j,ntopfl:nl)+v_ad(1,j,ntopfl:nl))       
        do i=1,imr-1
          va_ad(i,j,ntopfl:nl)=0.5d0*(v_ad(i,j,ntopfl:nl)+v_ad(i+1,j,ntopfl:nl))  
        enddo
      enddo

! Compute adjoint of tridiagonal solver
! -------------------------------------
      do j=jfirst, jlast

! Peel off 2-d slices of ua and va
! --------------------------------
        nfields=1
        if ( (j>1) .and. (j<jnp) ) then
          nfields=2
        endif 
        ca(1:imr,ntopfl:nl) = vcoef3d(1:imr,j,ntopfl:nl,2) ! cam   
        cc(1:imr,ntopfl:nl) = vcoef3d(1:imr,j,ntopfl:nl,4) ! ccm

        fp_ad(1:imr,ntopfl:nl,1)=ua_ad(1:imr,j,ntopfl:nl)
        if (nfields == 2) then 
          fp_ad(1:imr,ntopfl:nl,2)=va_ad(1:imr,j,ntopfl:nl) 
        endif

! apply tri-diagonal solver to wind fields
! ----------------------------------------
        call tridiag_ad_ ( imr, nl, ntopfl, nfields2, nfields, ca, cc, &
                           fp_ad, fm_ad )
        ua_ad(1:imr,j,ntopfl:nl)=fm_ad(1:imr,ntopfl:nl,1)-ua_ad(1:imr,j,ntopfl:nl)  
        if ( nfields == 2) then 
          va_ad(1:imr,j,ntopfl:nl)=fm_ad(1:imr,ntopfl:nl,2)-va_ad(1:imr,j,ntopfl:nl)
        endif 
      enddo ! loop over latitudes

! Compute adjoints of u and v on D grid from A-grid values
! --------------------------------------------------------
      call d2a3d_ad( u_ad, v_ad, ua_ad, &
      va_ad, imr, jnp, nl, jfirst, jlast, ng_d, ng_s, coslon, sinlon )

! adjoint for pt field
! The polavg routine is self-adjoint
! ----------------------------------
!     call polavg(  pt_ad(1,jfirst,k), imr, jnp, jfirst, jlast)
      ta_ad(1:imr,jfirst:jlast,ntopfl:nl)=pt_ad(1:imr,jfirst:jlast,ntopfl:nl)
      do j=jfirst, jlast

! Peel off 2-d slices of pt (and convert virt potential temp. to pot. T)
! ----------------------------------------------------------------------
        nfields=1
        ca(1:imr,ntopfl:nl) = vcoef3d(1:imr,j,ntopfl:nl,1)  ! cah 
        cc(1:imr,ntopfl:nl) = vcoef3d(1:imr,j,ntopfl:nl,3)  ! cch
        fp_ad(1:imr,ntopfl:nl,1)=ta_ad(1:imr,j,ntopfl:nl)*(1.d0+zvir*q(1:imr,j,ntopfl:nl,1))
        pt_ad(1:imr,j,ntopfl:nl)=pt_ad(1:imr,j,ntopfl:nl)-ta_ad(1:imr,j,ntopfl:nl)

! apply tri-diagonal solver to pt fields
! --------------------------------------
        call tridiag_ad_ ( imr, nl, ntopfl, nfields2, nfields, ca, cc, &
                           fp_ad, fm_ad )
        pt_ad(1:imr,j,ntopfl:nl)=pt_ad(1:imr,j,ntopfl:nl)+ &
                          fm_ad(1:imr,ntopfl:nl,1)/(1.d0+zvir*q(1:imr,j,ntopfl:nl,1))
      enddo ! loop over latitudes

      deallocate(vcoef3d)

end subroutine physdrv1_ad_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tridiag_ad_ --- solver for adjoint of implicit tri-diag. diffusion
!
! !INTERFACE:

subroutine tridiag_ad_ ( imr, nl, ntopfl, nfields2, nfields, &
                         ca, cc, fm_ad, fp_ad )

! !USES:
  implicit none

! !INPUT PARAMETERS:

  integer,  intent(in) :: imr      ! number of grid points in long. direction
  integer,  intent(in) :: nl       ! number of grid points in vert. direction
  integer,  intent(in) :: ntopfl   !
  integer,  intent(in) :: nfields2 ! max. number of fields to apply
  integer,  intent(in) :: nfields  ! number of fields to apply this call
  real(8), intent(in) :: ca(imr,nl) ! lower diagonal of matrix at imr points
  real(8), intent(in) :: cc(imr,nl) ! upper diagonal of matrix at imr points

! !INPUT/OUTPUT PARAMETERS:

  real(8), intent(inout) :: fm_ad(imr,nl,nfields2) ! fields before diffusion 

! !OUTPUT PARAMETERS:

  real(8), intent(out) :: fp_ad(imr,nl,nfields2) ! fields after diffusion 

! !DESCRIPTION: This is the adjoint of the tri-diagonal solver (subroutine 
!               tridiag).  Since the solver is identical except for 
!               transposing the matrix, the matrix and its solution are 
!               not described here (see instead subroutine tridiag). The 
!               transposition of the matrix is performed by using ca and cc
!               in place of each other, with adjusted indices, except 
!               for the use of ca and cc in defining the main diagonal, which
!               is unchanged. 
!
! !REMARKS:
!  Note that the routine is desined to input an array fm (same as the original
!               routine) and output an array fp.  So, for adjoint applications,
!               the arrays fm and fp will typically be reversed in the 
!               subroutine call.  Note also that the input array is altered 
!               (essentially used as a work space).

! !REVISION HISTORY:
!
!  6Jan2004  Errico  Initial code and implimentation
!
!EOP

!qqq  character(len=*), parameter :: myname_ = myname//'::tridiag_tl' 
  integer   :: i             ! longitudinal index  
  integer   :: k             ! vertical index  
  integer   :: km1           ! k-1    
  integer   :: kp1           ! k+1              
  integer   :: nf            ! field indicator   
  real(8)   :: work(imr,nl)  ! b_k - c_kE_k-1 where E_k=a_k/work_k  
  real(8)   :: b             ! diagonal element of matrix  
  real(8)   :: alpha         ! 

  do i=1,imr
    work(i,ntopfl)=1.d0+ca(i,ntopfl)
  enddo

  do k=ntopfl+1,nl
    km1=k-1
    do i=1,imr
      b=1.d0+ca(i,k)+cc(i,k)
      alpha=ca(i,km1)/work(i,km1)
      work(i,k)=b-alpha*cc(i,k)
      do nf=1,nfields
        fm_ad(i,k,nf)=fm_ad(i,k,nf)+alpha*fm_ad(i,km1,nf)
      enddo
    enddo 
  enddo

  do i=1,imr
    do nf=1,nfields
      fp_ad(i,nl,nf)=fm_ad(i,nl,nf)/work(i,nl)
    enddo
  enddo  

  do k=nl-1,ntopfl,-1
    kp1=k+1
    do i=1,imr
      do nf=1,nfields
        fp_ad(i,k,nf)=(fm_ad(i,k,nf)+cc(i,kp1)*fp_ad(i,kp1,nf))/work(i,k)
      enddo
    enddo
  enddo  

end subroutine tridiag_ad_

! 
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sponge_tl_ --- adjoint of solver for implicit sponge (same as TLM)
!
! !INTERFACE:
                                                                                              
subroutine sponge_tl_ ( imr, jnp, jfirst, jlast, nl, tdt, field, field_name )
                                                                                              
! !USES:
  implicit none
                                                                                              
! !INPUT PARAMETERS:
                                                                                              
  integer,  intent(in) :: imr      ! number of grid points in long. direction
  integer,  intent(in) :: jnp      ! number of grid points in lat. direction
  integer,  intent(in) :: jfirst, jlast   ! range of latitude indeces
  integer,  intent(in) :: nl       ! number of grid points in vert. direction
  real(8),  intent(in) :: tdt      ! time step in seconds
  character(len=*), intent(in) :: field_name ! to denote grid index
                                                                                              
! !INPUT/OUTPUT PARAMETERS:
                                                                                              
  real(8), intent(inout) :: field(imr,jfirst:jlast,nl) ! field to be damped
                                                                                              
! !OUTPUT PARAMETERS:
                                                                                              
                                                                                              
! !DESCRIPTION:
!    This routine applies a sponge to prognostic fields near the model's
!    top (levels k=1,...,kmax=nl/3).  The damping is linear df/dt=-cf.  It is
!    is solved as f(t+delta t)=f(t)*exp(-c*delta t).  The e-folding
!    damping rate c decreases geometrically as lower levels are considered,
!    from 1/coef1 at k=1 to 1/coef2 at k=kamx.
                                                                         
! !REMARKS:
                                                                                              
! !REVISION HISTORY:
!
!  20Feb2005  Errico  Initial code and implimentation
!  19Apr2006  Errico/Elena N. Adjusted sponge parameters for GEOS-5
!  17Sep2007  Errico/Elena N. Re-adjusted coeffs based on strat. SVs
!
!EOP

  character(len=*), parameter :: myname_ = myname//'::sponge_tl_'
  integer   :: k             ! vertical index
  integer   :: kmax          ! max k at which sponge is applied
  integer   :: j1,j2         ! max and min of lat index to consider
  real(8)   :: coefs(nl)     ! damping factors for each level
  real(8)   :: coef1         ! e-folding damping time (days) at k=1
  real(8)   :: coef2         ! e-folding damping time (days) at k=kmax
  real(8)   :: coef          ! e-folding damping rate at level k
  real(8)   :: coefx         ! factor by which the damping rate decreases

                             ! for each successive level k  
  coef1=0.7d0
  coef2=30.d0

  kmax=nl/2

  coefx=(coef1/coef2)**(1./real(kmax-1))
  coef=1.d0/(coef1*coefx*86400.d0)
 
  do k=1,kmax
    coef=coef*coefx
    coefs(k)=exp(-tdt*coef)
  enddo
!
! account for extent of field on D-grid
  if (field_name == 'u') then
    j1=max(jfirst,2)
    j2=jlast
  elseif (field_name == 'v') then
    j1=max(jfirst,2)
    j2=min(jlast,jnp-1)
  else
    j1=jfirst
    j2=jlast
  endif
 
  do k=1,kmax
    field(1:imr,j1:j2,k)=coefs(k)*field(1:imr,j1:j2,k)
  enddo
 
end subroutine sponge_tl_
! 
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine hdiff_ad (imr,jnp,nl,jfirst,jlast,dt,pt_ad,u_ad,v_ad,ng_d,ng_s)

      implicit none

      integer,  intent(in)  :: imr
      integer,  intent(in)  :: jnp
      integer,  intent(in)  :: nl
      integer,  intent(in)  :: jfirst
      integer,  intent(in)  :: jlast
      integer,  intent(in)  :: ng_d
      integer,  intent(in)  :: ng_s
      real(r8), intent(in)  :: dt 

      real(r8), intent(inout)  :: pt_ad(imr,jfirst-ng_d:jlast+ng_d,nl) 
      real(r8), intent(inout)  :: u_ad(imr,jfirst-ng_d:jlast+ng_s,nl) 
      real(r8), intent(inout)  :: v_ad(imr,jfirst-ng_s:jlast+ng_d,nl) 
      
      integer  :: k, i, imh, ntime
      integer  :: ifax(13)         ! factorization of imr 
      real(r8) :: trigs(3*imr/2+1) ! trig factors for FFT
      real(r8) :: zam5, zamda, dp, dl
      real(r8) :: cosp(jnp), cose(jnp), sinp(jnp), sine(jnp)
      real(r8) :: cosl5(imr), sinl5(imr), coslon(imr), sinlon(imr)

        call setrig (imr, jnp, dp, dl, cosp, cose, sinp, sine)

        imh=imr/2
        do i=1,imh
           zam5          = (dble(i-1)-0.5d0) * dl
           cosl5(i)      = -dcos(zam5)
           cosl5(i+imh)  = -cosl5(i)
           sinl5(i)      = -dsin(zam5)
           sinl5(i+imh)  = -sinl5(i)
           zamda         = dble(i-1)*dl
           coslon(i)     = -dcos(zamda)
           coslon(i+imh) = -coslon(i)
           sinlon(i)     = -dsin(zamda)
           sinlon(i+imh) = -sinlon(i)
        enddo

!
! compute factors required for the FFT
        call fftfax (imr,ifax,trigs)
!
! assume that there is no sensitivity to pole values of v here
        if (jfirst==1) then 
          v_ad(:,  1,:)=0.d0
        endif
        if (jlast==jnp) then          
          v_ad(:,jnp,:)=0.d0
        endif
!
        call  hdiff_tend_ad (imr,jnp,jfirst,jlast,nl,cosp,cose, &
                            coslon,sinlon,dt,ifax,trigs,pt_ad,'S',ng_d,ng_d)
! 
        call  hdiff_tend_ad (imr,jnp,jfirst,jlast,nl,cose,cosp, &
                            coslon,sinlon,dt,ifax,trigs,u_ad,'u',ng_d,ng_s)
!
        call  hdiff_tend_ad (imr,jnp,jfirst,jlast,nl,cosp,cose, &
                            coslon,sinlon,dt,ifax,trigs,v_ad,'v',ng_s,ng_d)
      end subroutine hdiff_ad
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine hdiff_fac_ad (imr,dt,xnu,cosp,damp)
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
!  ! DESCRIPTION:
!     Compute coefficients for del**4 implicit damping
!
      implicit none

      integer,  intent(in) :: imr
      real(r8), intent(in) :: dt        ! physics time step in seconds
      real(r8), intent(in) :: xnu       ! Laplacian damping factor
      real(r8), intent(in) :: cosp      ! cosine at field points

      real(r8), intent(out) :: damp(imr)
!
! Compute coefficients for Fourier filter of initial fields in the
! tangent linear model or final fields in the adjoint model.
!
! The filtering is done based on a function of the zonal wavenumber k,
! accounting for the cosine factor in the length of the latitude circle.
! The function here is:
!
!      df/dt = xnu*(d^4 c \ d x^4)
!      c(t+dt)=c(t)-dt*xnu*m**4*c(t+dt)
!      k=m*2*pi/L ; L=2*pi*R*cos(lat);  k=m/(R*cos(lat))      
!      c(t+dt)=(1/(1+dt*xnu*m**4))c(t)
!
! The indexes for the array damp(i) are i for
! dimensionless zonal wavenumber m=int((i+1)/2). Odd and even i correspond
! respectively to real and imaginary components of the Fourier
! coefficients (those damping coefs should be identical). m is the 
! wavenumber measured in multiples of the fundamental wavenumber (that 
! for the largest zonal wavelength at each particular latitude).
!    
! damp varies between 0 and 1. Smaller values indicate greater filtering.

      integer :: i,j,m
      real(r8) :: fac,faccos,xx,rxx
!
      fac=sqrt(xnu*dt)/earthr**2
      if (cosp /= 0.d0 ) then
        faccos=fac/cosp**2
        do i=2,imr,2
          m=i/2
          xx = faccos*m*m
          rxx = 1.0d0/(1.0d0+xx*xx)
          damp(i)=rxx
          damp(i-1)=damp(i)
        enddo
      else
          damp=1.d0
      endif
!
      end subroutine hdiff_fac_ad 
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
     subroutine hdiff_zonal_ad (imr,jfirst,jlast,j1,j2,nl, &
                            dt,xnu,cosp,ifax,trigs,f)

!
! !DESCRIPTION: This routine applies a filter at latitudes close to poles
!               for calculation of singular vectors to
!               a single 2-d (lat-lon) field.
! Selection of near-pole-latitudes is based on 'fcrit' value setup in module's
! global constants.
! All filtering performed by applying zonal FFT to fields, multiplying 
! zonal spectral coefficients by a supplied factor 0 <= damp <= 1, and 
! projecting back onto the grid, replacing the original field.
!
 
     implicit none

! Input
     integer,  intent(in) :: imr  
     integer,  intent(in) :: jfirst,jlast,j1,j2
     integer,  intent(in) :: nl
     integer,  intent(in) :: ifax(13)        ! factorization of im
     real(r8), intent(in) :: trigs(3*imr/2+1) ! trig factors for FFT
     real(r8), intent(in) :: dt,xnu
     real(r8), intent(in) :: cosp(j1:j2)

! Input/Output
     real(r8), intent(inout) :: f(imr,jfirst:jlast,nl)

! Local
     integer   :: i,j,k
     integer   :: imrm1,imrm2,imrp1,imrp2
     integer   :: mmax
     real(r8)  :: facx
     real(r8)  :: damp(imr+2)    ! implicit spectral damping factors
     real(r8)  :: fs(imr+2,nl)   ! dt*tendency, or array for copying fields
     real(r8)  :: work(imr+2,nl) ! work array for FFT
     real(r8)  :: dlon_eq, afac, afacx
!
     imrm1=imr-1
     imrm2=imr-2
     imrp1=imr+1 
     imrp2=imr+2
!
      dlon_eq=earthr*twopi/real(imr)    ! long grid spacing at equator (km)
      afac=xnu*dt
      afacx=afac/dlon_eq**4
      mmax = imr/2

     do j=j1,j2
!
! facx i
!
       facx=(xnu*dt*imr**4.)/(pi*earthr*cosp(j))**4.

       if (facx.ge.fcrit ) then       ! compute using discrete implicit scheme

! Determine spectral damping factors
!         
         call hdiff_fac_ad (imr,dt,xnu,cosp(j),damp)
!
! Copy 1 lat of the field to new array structure 
         do k=1,nl
           fs(1:imr,k)=f(1:imr,j,k)
           fs(imrp1,k)=0.d0    
           fs(imrp2,k)=0.d0    
         enddo
!
! Replace field values in q1 by Fourier coefficients
! These coefs are ordered real, imaginary, real, imaginary, ...
! for dimensionless wavenumber 0,0, 1,1, 2,2, ..., im/2, im/2)
!
         call rfftmlt(fs, work, trigs, ifax, 1, imrp2, imr, nl, -1)
!
! Damp spectral coefficients.
! Note that the array damp is presumed to be ordered the same as 
! the spectral coefficients except that the values start corresponding 
! to wavenumber 1, rather than at wavenumber 0 as for q1.
!

        do k=1,nl
          do i=3,imrp2
            fs(i,k) = fs(i,k)*damp(i-2)
          enddo
        enddo

!
! Compute filtered fields from filtered coefficients

        call rfftmlt(fs, work, trigs, ifax, 1, imrp2, imr, nl, 1)

!  Adjoint to Laplacian:

        do k=1,nl       ! Loop for k-levels
             f(1:imr,j,k) = fs(1:imr,k)
        enddo
!
      endif  ! test on whether implicit spectral
    enddo    ! loop over j
!
 end subroutine hdiff_zonal_ad
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
      subroutine hdiff_del2_ad (imr,jnp,jfirst,jlast,nl,cosp,cose,   &
                     coslon, sinlon, xnu, dt, f, del2f,name,ng_d1,ng_d2)
!__________________________________
!
!   ! DESCRIPTION: 
!       Compute explicit Laplacian of a field f
!
      implicit none

      integer,  intent(in) :: imr
      integer,  intent(in) :: jnp
      integer,  intent(in) :: jfirst
      integer,  intent(in) :: jlast
      integer,  intent(in) :: nl
      real(r8), intent(in) :: cosp(jnp) ! cosine at field points
      real(r8), intent(in) :: cose(jnp) ! cosine at edges between points
      real(r8), intent(in) :: coslon(imr)  !
      real(r8), intent(in) :: sinlon(imr)  !
      real(r8), intent(in) :: xnu, dt
      integer :: ng_d1
      integer :: ng_d2
      real(r8), intent(inout) :: f(imr,jfirst-ng_d1:jlast+ng_d2,nl)
! Note: this nomenclature for cosp,e differs from FVGCM: for
! u field points, it is the reverse of the FVGCM specs, which refers to 
! T field points. 
      character(len=*), intent(in) :: name  ! 'u', 'v', or 'S'
     
      real(r8), intent(in) :: del2f(imr,jfirst-ng_d1:jlast+ng_d2,nl)

! local variables
      integer  :: j1, j2, ja, jb, jp, jm, jnume
      integer  :: i, j, k
      integer  :: nmax
      real(r8) :: afac, afacx, afacy
      real(r8) :: dlon_eq, dle, dlp, dlat, area, sine
      real(r8) :: fac, fac0, facm, facp, facx, fac1
      real(r8) :: facm_glob(jnp), facp_glob(jnp), fac0_glob(jnp)
      real(r8) :: sum_diff
      real(r8) :: cvr, cvi 
      real(r8) :: fcheck
      integer  :: jstart, jend
      character(len=128) :: myrankname,myrankname1
      logical, save :: openghost = .true.
!
! Set indecies for special treatment of first and last rows of fields
      if (name=='u') then
        j1=max(jfirst,3)
      else
        j1=max(jfirst,2)
      endif
      j2=min(jlast,jnp-1)
!
      dlon_eq=earthr*twopi/real(imr)  ! long grid spacing at equator (km)
      dlat=earthr*pi/(jnp-1)          ! lat grid spacing in km 
      afac=sqrt(xnu*dt)               ! by appling this routine 2x, this gets squared
      afacx=afac/dlon_eq**2
      afacy=afac/dlat**2


! Compute del**2 for points where pole is not an issue.
!
#ifdef SPMD
      call mp_barrier
      call mp_send3d_ns(imr,jnp,jfirst,jlast,1,nl,ng_d1,ng_d2,del2f,1)
      call mp_barrier
      call mp_recv3d_ns(imr,jnp,jfirst,jlast,1,nl,ng_d1,ng_d2,del2f,1)
#endif
      do j=2,jnp-1
        jp=j+1
        jm=j-1
        if (name == 'u') then
          facp_glob(j)=afacy*cose( j)/cosp(j)
          facm_glob(j)=afacy*cose(jm)/cosp(j)
        else
          facp_glob(j)=afacy*cose(jp)/cosp(j)
          facm_glob(j)=afacy*cose( j)/cosp(j)
        endif
        fac0_glob(j)=-(facp_glob(j)+facm_glob(j))
      enddo

      if(j1 == 2 .or. j1==3) then
        f(:,j1-1,:)=f(:,j1-1,:)+facm_glob(j1)*del2f(:,j1,:)
        f(:,j1,:)  =f(:,j1,:)  +fac0_glob(j1)*del2f(:,j1,:)+facm_glob(j1+1)*del2f(:,j1+1,:)
        if (name=='u') then
          jstart = 4
        else
          jstart = 3
        endif
      else
        jstart = j1 
      endif
      if(j2 == jnp-1) then
        f(:,j2+1,:)=f(:,j2+1,:)+facp_glob(j2)*del2f(:,j2,:)
        f(:,j2,:)  =f(:,j2,:)+facp_glob(j2-1)*del2f(:,j2-1,:)+fac0_glob(j2)*del2f(:,j2,:)
        jend = jnp-2
      else
        jend = j2
      endif

      do j=jstart,jend
        jp=j+1
        jm=j-1
         do k=1,nl
          do i=1,imr
            f(i,j,k)  = f(i,j,k)  + facp_glob(jm)*del2f(i,jm,k)   &
                         + fac0_glob(j)*del2f(i,j ,k)   &
                         + facm_glob(jp)*del2f(i,jp,k)
          enddo
        enddo
      enddo
!
! Consider contribution by longitudinal gradients
       
      if (name=='u') then
        j1=max(jfirst,2)
        j2=min(jlast,jnp)
      endif
       
      do j=j1,j2
        fcheck = (xnu*dt*imr**4.)/(pi*earthr*cosp(j))**4.
!
! Do not apply near poles where delta x is too small and will lead to a numerical
! instability: imr/(pi*earthr*cosp)**4 is the factor you get if you apply 
! d**4/dx**4 to an f that is a 2 delta x wave using centered finite differences
        if ( fcheck < fcrit )  then
!                        
          facx=afacx/cosp(j)**2  ! this is sqrt(xnu*dt)(1/delta x)**2 at this latitude
          fac0=-2.0d0*facx
          do k=1,nl
            f(  1,j,k)  =  f(  1,j,k)    + fac0*del2f(  1,j,k)
            f(  2,j,k)  =  f(2,j,k)      + facx*del2f(  1,j,k)
            f(  imr,j,k)=  f(  imr,j,k)  + facx*del2f(  1,j,k) 
            f(  imr,j,k)=  f(  imr,j,k)  + fac0*del2f(imr,j,k)
            f(imr-1,j,k)=  f(imr-1,j,k)  + facx*del2f(imr,j,k)
            f(  1,j,k)  =  f(  1,j,k)    + facx*del2f(imr,j,k)   
            do i=2,imr-1
              f(  i,j,k) = f(  i,j,k) + fac0*del2f(i,j,k)
              f(i+1,j,k) = f(i+1,j,k) + facx*del2f(i,j,k)
              f(i-1,j,k) = f(i-1,j,k) + facx*del2f(i,j,k)
            enddo
          enddo

        endif   ! for fcheck
      enddo
!
      if ((name=='S') .and. ((j1>jfirst) .or. (j2<jlast)) ) then
        sine=sqrt(1.0d0-cose(2)**2)
        area=earthr*dlon_eq*(1.0d0-sine)    ! area of cap / imr 
        dle=cose(2)*dlon_eq                 ! length of cap / imr
        fac=afac*dle/(area*dlat*imr)    
!
        if (j1 > jfirst) then ! input range of lats includes S.P.
          ja=1
          jb=2
          do k=1,nl
            sum_diff=0.d0   
            do i=1,imr
              sum_diff=sum_diff+del2f(i,ja,k)  
            enddo
            sum_diff=fac*sum_diff
            f(:,jb,k)=f(:,jb,k) + sum_diff    
            f(:,ja,k)=f(:,ja,k) - sum_diff
          enddo 
        endif  ! test if S.P. included 
!
        if (j2 < jlast) then ! input range of lats includes N.P.
          ja=jnp
          jb=jnp-1
          do k=1,nl
            sum_diff=0.d0   
            do i=1,imr
              sum_diff=sum_diff+del2f(i,ja,k)     
            enddo
            sum_diff=fac*sum_diff
            f(:,jb,k)=f(:,jb,k) + sum_diff
            f(:,ja,k)=f(:,ja,k) - sum_diff
          enddo 
        endif  ! test if N.P. included 
!
      endif    ! test if Scalar field at poles included
!
      if ((name=='v') .and. ((j1>jfirst) .or. (j2<jlast)) ) then
!
        if (j1 > jfirst) then ! input range of lats includes S.P.
          ja=1
          jb=2
          do k=1,nl
            cvr=0.d0
            cvi=0.d0
            do i=1,imr
              cvr=cvr+coslon(i)*f(i,ja,k) 
              cvi=cvi+sinlon(i)*f(i,ja,k)
              f(i,ja,k)=0.d0
            enddo
            cvr=2.0d0*cvr/real(imr, r8)
            cvi=2.0d0*cvi/real(imr, r8)  
            do i=1,imr
              f(i,jb,k)=f(i,jb,k) + cvr*coslon(i)+cvi*sinlon(i)
            enddo
          enddo
!
        endif  ! test if next to S.P. included
!
        if (j2 < jlast) then ! input range of lats includes N.P.
!                      
          ja=jnp
          jb=jnp-1
          do k=1,nl
            cvr=0.d0
            cvi=0.d0
            do i=1,imr
              cvr=cvr+coslon(i)*f(i,ja,k)
              cvi=cvi+sinlon(i)*f(i,ja,k)
              f(i,ja,k)=0.d0
            enddo
            cvr=2.0d0*cvr/real(imr, r8)
            cvi=2.0d0*cvi/real(imr, r8)
            do i=1,imr
              f(i,jb,k)=f(i,jb,k) + cvr*coslon(i)+cvi*sinlon(i)
            enddo
          enddo
                 
        endif  ! test if next to N.P. included
!
      endif    ! test if v field next to poles included
!
!
! reset j1 and j2 properly for u field
      if (name=='u') then
        j1=max(jfirst,3)
        j2=min(jlast,jnp-1)
      endif
!
      if ((name=='u') .and. ((j1>jfirst) .or. (j2<jlast)) ) then
!
        sine=sqrt(1.0d0-cose(2)**2)
        area=earthr*dlon_eq*(1.0d0-sine)    ! area of sector
        dle=cose(2)*dlon_eq                 ! length of sector's base
        fac0=afac*dle/(dlat*area)
        fac1=-fac0
!
        if (j1 > jfirst) then ! input range of lats includes S.P.
          ja=2
          jb=3
          do k=1,nl
            f(:,ja,k)=f(:,ja,k)+fac1*del2f(:,ja,k)
            f(:,jb,k)=f(:,jb,k)+fac0*del2f(:,ja,k)
          enddo
        endif  ! test if next to S.P. included
!
        if (j2 < jlast) then ! input range of lats includes N.P.
          ja=jnp
          jb=jnp-1
          do k=1,nl
            f(:,ja,k)=f(:,ja,k)+fac1*del2f(:,ja,k)
            f(:,jb,k)=f(:,jb,k)+fac0*del2f(:,ja,k)
          enddo
        endif  ! test if next to N.P. included
!
      endif    ! test if u field next to poles included
!
      end subroutine hdiff_del2_ad



!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine hdiff_tend_ad (imr,jnp,jfirst,jlast,nl,cosp,cose, &
                             coslon,sinlon,dt,ifax,trigs,f,name,ng_d1,ng_d2)
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
                                                                                                                             
      implicit none

      integer,  intent(in) :: imr
      integer,  intent(in) :: jnp
      integer,  intent(in) :: jfirst
      integer,  intent(in) :: jlast
      integer,  intent(in) :: nl
      integer,  intent(in) :: ifax(13)         ! factorization of imr
      real(r8), intent(in):: trigs(3*imr/2+1) ! trig factors for FFT
      real(r8), intent(in) :: cosp(jnp) ! cosine at field points
      real(r8), intent(in) :: cose(jnp) ! cosine at edges between points
      real(r8), intent(in) :: coslon(imr),sinlon(imr) ! only used for v 
      real(r8), intent(in) :: dt        ! step for time marching
! Note: this nomenclature for cosp,e differs from FVGCM: for
! u field points, it is the reverse of the FVGCM specs, which refers to
! T field points.
      character(len=*), intent(in) :: name  ! 'u', 'v', or 'S'
      integer :: ng_d1
      integer :: ng_d2
!
      real(r8), intent(inout) :: f(imr,jfirst-ng_d1:jlast+ng_d2,nl)
!
      real(r8) :: del4f(imr,jfirst-ng_d1:jlast+ng_d2,nl)
      real(r8) :: work(imr,jfirst-ng_d1:jlast+ng_d2,nl)
!
! local variables
!
      integer  :: j1, j2
      integer  :: nmax
      real(r8) :: tdamp_secs, xnu
      integer  :: j
      real(r8) :: facx

!
!    Compute the tendency of del**4 by applying hdiff_del2 twice   
!    Then adjust the tendency replacing zonal part of explicit del**4 
!    by implicit zonal damping at latitudes "close" to the poles 
!    "Close" criteria is specified within routine 'new_hdiff_zonal_tl' itself 
! 
!
      nmax=imr/2
      tdamp_secs=max(dt*tdamp_dt,tdamp_days*86400.d0)
      xnu=earthr**4/(tdamp_secs*(nmax**2+nmax)**2)
!
      xnu = xnu*(0.5*pi)**4        ! The factor (pi/2)^4 accounts for the effective 
                                   ! scalar w that arises from wf=(d**4 /dx**4)f where
                                   ! f=A*sin((i-1)*delta x/nmax) and the derivative
                                   ! is computed as discrete centered differences. 
                                   ! Without this factor of pi/2, w is 1/(n*(n-1))**2,
                                   ! which is the factor for del**4 applied to a spherical   
                                   ! harmonic of scale n=nmax.
!
      j1=max(jfirst,2)
      if (name=='u') then
        j2=min(jlast,jnp)
      else
        j2=min(jlast,jnp-1)
      endif
!
      del4f(:,:,:) = -f(:,:,:)

      call hdiff_zonal_ad(imr,jfirst,jlast,j1,j2,nl, &
                            dt,xnu,cosp(j1:j2),ifax,trigs,f(:,jfirst:jlast,:))
      work(:,:,:)  = 0.d0 
      call hdiff_del2_ad(imr,jnp,jfirst,jlast,nl,cosp,cose, &
                      coslon,sinlon,xnu,dt,work,del4f, name,ng_d1,ng_d2)
!
      call hdiff_del2_ad(imr,jnp,jfirst,jlast,nl,cosp,cose, &
                      coslon,sinlon,xnu,dt,f,work, name,ng_d1,ng_d2)
!
!
      end subroutine hdiff_tend_ad
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
end module m_physdrv1_ad
