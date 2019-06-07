!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_physdrv1_tl --- Tangent linear of simplified physics
!
! !INTERFACE:
!
module m_physdrv1_tl

! !USES:

      use precision
      use m_const, only : zvir
      use m_StrTemplate, only : strTemplate
      use stepon, only : nsecf
      use stepon, only : tdt
      use m_die, only : die
      use m_trjphys, only: physdrv1_get_state

#ifdef SPMD
      use mod_comm, only : mp_send_s,mp_recv_n,mp_send3d_ns,mp_recv3d_ns, &
             mp_barrier,mp_send2_n,mp_recv2_s
      use mod_comm, only : gid
#endif
!
  implicit none

! !PUBLIC MEMBER FUNCTIONS:

  PRIVATE

  PUBLIC  physdrv1_tl

  interface physdrv1_tl ; module procedure physdrv1_tl_ ; end interface

!
! !DESCRIPTION: This module implements the first version of the
!               tangent linear simplified physics (vertical diffusion).
!
! !REMARKS:
!   The following routines have been hand-modified from their original
!   automatic generation: d2a2_tl1_, d2a3d_tl1_,
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

   character(len=*), parameter :: myname = 'm_physdrv1_tl'
   integer,          parameter :: ROOT = 0

!
!  WARNING **** The horiz diff  code will only work correctly for jfirst=1
!  jlast=jnp; The problem is that extra latitudes are required to
!  compute the gradients at the 2 edges of each sub domain.  An MPI version
!  will therefore require modification.
!
!  Also Note that the parameters in the hdiff_tl routines must be set the same
!  as in the hdiff_ad routines until the code is restructured. 

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
! !ROUTINE: physdrv1_tl_ --- version 1 of driver for simplified physics
!
! !INTERFACE:
 
subroutine physdrv1_tl_ ( ptrjtmpl, ptrjfrq, job, imr,  jnp,   nl,   nc, jfirst, jlast,  &
                           coslon,   sinlon,    q, nymd, nhms, &
                            u_tl,      v_tl, pt_tl, ng_d, ng_s)
 
! !USES:

      use m_trjphys, only: nvars

      implicit none
 
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: ptrjtmpl    ! filename template for vdiff fields
      character(len=*), intent(in) :: job         ! name of experiment. DO THIS SOME OTHER WAY!!!
      integer, intent(in)  :: ptrjfrq             ! frequecy (HHMMSS) of trajectory update
                                                  !   < 0 returns without exercising routine
                                                  !   = 0 doesn't update trajectory
                                                  !   > 0 updates trajectory
      integer, intent(in)  :: imr                 ! number of grid points in long. direction
      integer, intent(in)  :: jnp                 ! number of grid points in lat.  direction
      integer, intent(in)  :: nl                  ! number of grid points in vert. direction
      integer, intent(in)  :: nc                  ! number of constituent fields
      integer, intent(in)  :: jfirst              !
      integer, intent(in)  :: jlast               !
      integer, intent(in)  :: ng_d               !
      integer, intent(in)  :: ng_s               !
      real(kind=r8), intent(in) :: coslon(imr)    ! cosine of longitude
      real(kind=r8), intent(in) :: sinlon(imr)    ! sine   of longitude
      real(kind=r8), intent(inout) :: q(imr,jfirst-ng_d:jlast+ng_d,nl,nc)
      integer   nymd                              !  Year-month-day, e.g.,  19971012
      integer   nhms                              !  hour-minute-sec, e.g.,   120000

! !INPUT/OUTPUT PARAMETERS:
      real(kind=r8), intent(inout) :: u_tl(imr,jfirst-ng_d:jlast+ng_s,nl)
      real(kind=r8), intent(inout) :: v_tl(imr,jfirst-ng_s:jlast+ng_d,nl)
      real(kind=r8), intent(inout) :: pt_tl(imr,jfirst-ng_d:jlast+ng_d,nl)   ! virt pt

! !OUTPUT PARAMETERS:
                     

! !DESCRIPTION: This routine computes u, v, and pt fields modified by
!               vertical diffusion applied as a temporally implicit scheme.
!               Input are fields after adiabatic step.  The modified fields
!               overwrite those values.  This simplified physics uses
!               prescribed vertical diffusion factors described in terms of
!               their corresponding matrix elements (the factors are products
!               of eddy diffusion coefs, air density, 1/delta z, and time
!               step).  These coefficients are 3D varying.
                                                        
! !REMARKS:  The polar filter used in high resolution versions of the model
!            has not been inserted here.  What is called the upper diagonal
!            here is called the lower diagonal in CCM3 and vice versa, because
!            here it is considered that the matrix rows are ordered with
!            k=1 (the topmost atmospheric level) at the top of the matrix.
!            It is assumed that pt is virtual potential temperature on input,
!            as it is on output. q is only used to create pot. temp. from pt.
                                                        
                                                        
! !REVISION HISTORY:
!
!  10Mar2003  Errico  Initial code and implementation
!  07Jan2004  Winslow modified variable names, rearranged some loops, 
!                     and placed into TLM context
!  10Feb2004  Todling Added template ref for i/o traj filename;
!                     added feature to allow skipping physics (simply don't
!                     provide a phys. trajectory file); added param to 
!                     control trajectory update.
!  10Oct2005  Errico  Corrections added to Hdiff
!  20Apr2007  Kim     Updated interface to physics trajectory package
!  01Jul2007  Todling Simplified: 2nd dim of vcoef need only jf to jl
!
!EOP
!-------------------------------------------------------------------------

! !LOCAL VARIALBLES

      character(len=*), parameter :: myname_ = myname//':physdrv1_tl_'

      character*255 :: vFile                      ! name of file containing coefs
      real(kind=r8), allocatable :: vcoef3d(:,:,:,:)  ! upper and lower diagonal elements of matrix
                                                  ! used in vdiff (1=cah,2=cam,3=cch,4=ccm)
      real(kind=r8) ua_tl(imr,jfirst-1:jlast, nl)   ! u field on a grid (tangent linear)
      real(kind=r8) va_tl(imr,jfirst:jlast, nl)   ! v field on a grid (tangent linear)
      real(kind=r8) ta_tl(imr,jfirst:jlast, nl)   ! change in potential temperature
      real(kind=r8) ca(imr,nl)                    ! -upper diag for momentum or heat and constituts
      real(kind=r8) cc(imr,nl)                    ! -lower diag for momentum or heat and constits
      integer, parameter :: nfields2=2            ! max number of simultaneous field slices
      real(kind=r8) fm_tl(imr,nl,nfields2)        ! wind or potential temperature fields before diffusion
      real(kind=r8) fp_tl(imr,nl,nfields2)        ! wind fields after diffusion
      integer       nfields                       !
      integer :: ntopfl
      integer i,j,k, j1
      integer ierr, iouter_phys

!     Return when no simplified physics not wanted 
!     --------------------------------------------
      if (ptrjfrq<0) return

!     Allocate vertical specify diffusion coefficients array
!     ------------------------------------------------------
      allocate ( vcoef3d(imr, jfirst:jlast, nl, nvars), stat=ierr )
          if(ierr/=0) call die(myname,' Alloc(vcoef3d) error',ierr)

!     Get the diffusion coefficents for this time
!     -----------------------------
      if ( ptrjfrq == 0 ) then
         if(gid==ROOT) print *, myname_, ': this option not yet implemented '
         call exit(7)
      else
         call physdrv1_get_state(vcoef3d,nymd,nhms) 
      endif
 
!     convert winds to the a grid
      call d2a3d_tl1( u_tl(:,jfirst:jlast+1,:), v_tl, ua_tl(:,jfirst:jlast,:), va_tl, imr, jnp, nl, &
                     jfirst, jlast, ng_d, ng_s, coslon, sinlon )
 
!     loop over latitude
      ntopfl = 1
      do j = jfirst, jlast
 
!     Peel off 2-d slices of u and v (do not need pole values of v)
        fm_tl(1:imr,ntopfl:nl,1) = ua_tl(1:imr,j,ntopfl:nl)
        nfields=1
        if ( (j>1) .and. (j<jnp) ) then
          fm_tl(1:imr,ntopfl:nl,2) = va_tl(1:imr,j,ntopfl:nl)
          nfields=2
        endif
        ca(1:imr,ntopfl:nl) = vcoef3d(1:imr,j,ntopfl:nl,2)
        cc(1:imr,ntopfl:nl) = vcoef3d(1:imr,j,ntopfl:nl,4)
            
!     apply tri-diagonal solver to wind fields
        call tridiag_tl_ ( imr, nl, ntopfl, nfields2, nfields, ca, cc, &
                           fm_tl, fp_tl )
            
!     compute change in winds due to physics
        ua_tl(1:imr,j,ntopfl:nl) = fp_tl(1:imr,ntopfl:nl,1) - ua_tl(1:imr,j,ntopfl:nl)
        if (nfields == 2) then
          va_tl(1:imr,j,ntopfl:nl) = fp_tl(1:imr,ntopfl:nl,2) - va_tl(1:imr,j,ntopfl:nl)
        endif

!     Peel off 2-d slices of pt (and convert virt potential temp. to pot. T)
        nfields = 1
        fm_tl(1:imr,ntopfl:nl,1) = pt_tl(1:imr,j,ntopfl:nl) / &
                         (1.d0+zvir*q(1:imr,j,ntopfl:nl,1))
        ca(1:imr,ntopfl:nl) = vcoef3d(1:imr,j,ntopfl:nl,1)
        cc(1:imr,ntopfl:nl) = vcoef3d(1:imr,j,ntopfl:nl,3)
            
!     apply tri-diagonal solver to pt fields
        call tridiag_tl_ ( imr, nl, ntopfl, nfields2, nfields, ca, cc, &
                           fm_tl, fp_tl )
            
!     compute change in pt due to physics (virtual pot. temperature here)
        ta_tl(1:imr,j,ntopfl:nl) =  fp_tl(1:imr,ntopfl:nl,1) * &
                            (1.d0+zvir*q(1:imr,j,ntopfl:nl,1)) &
                           - pt_tl(1:imr,j,ntopfl:nl)
      enddo ! loop over latitudes

!     The tendencies then get filtered and mapped back to the wind form on the D-grid
!     Filter changes to the fields on D grid using routine pft2D
!     NO FILTERING DONE HERE (in physdrv, filter applied only if imr>144)
!     This filter should otherwise be applied to ua, va, ta
            
!     Add (possibly filtered) changes of pt to pt
      pt_tl(1:imr,jfirst:jlast,ntopfl:nl) = pt_tl(1:imr,jfirst:jlast,ntopfl:nl) + &
                        ta_tl(1:imr,jfirst:jlast,ntopfl:nl)
            
!     Add change to u by remapping changes from D grid
!     Compute u at jfirst later; need ua(j=jfirst-1)
#ifdef SPMD
      call mp_barrier
      call mp_send2_n(imr,jnp,jfirst,jlast,1,nl,1,0,ua_tl,ua_tl)
      call mp_barrier
      call mp_recv2_s(imr,jnp,jfirst,jlast,1,nl,1,0,ua_tl,ua_tl)
#endif
      j1=max(2,jfirst)
!      do j=jfirst+1,jlast
      do j=j1,jlast
          u_tl(1:imr,j,ntopfl:nl) = u_tl(1:imr,j,ntopfl:nl) + 0.5d0 * &
                        (ua_tl(1:imr,j,ntopfl:nl)+ua_tl(1:imr,j-1,ntopfl:nl))
      enddo

!     Add change to v by remapping changes from D grid
      do j=max(2,jfirst),min(jnp-1,jlast)
        v_tl(1,j,ntopfl:nl) = v_tl(1,j,ntopfl:nl) + 0.5d0 * &
                        (va_tl(1,j,ntopfl:nl)+va_tl(imr,j,ntopfl:nl))
        do i=2,imr
          v_tl(i,j,ntopfl:nl) = v_tl(i,j,ntopfl:nl) + 0.5d0 * &
                        (va_tl(i,j,ntopfl:nl)+va_tl(i-1,j,ntopfl:nl))
        enddo
      enddo

      deallocate ( vcoef3d )

!
!  Apply sponge in top portion of model
     call sponge_tl_ ( imr, jnp, jfirst, jlast, nl, tdt, pt_tl(:,jfirst:jlast,:), 'pt' )
     call sponge_tl_ ( imr, jnp, jfirst, jlast, nl, tdt, u_tl(:,jfirst:jlast,:),  'u'  )
     call sponge_tl_ ( imr, jnp, jfirst, jlast, nl, tdt, v_tl(:,jfirst:jlast,:),  'v'  )

!
! Horizontal diffusion
    call hdiff_tl (imr,jnp,nl,jfirst,jlast,tdt,pt_tl,u_tl,v_tl,ng_d,ng_s)

end subroutine physdrv1_tl_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tridiag_tl_ --- solver for implicit tri-diagonal diffusion
!
! !INTERFACE:

subroutine tridiag_tl_ ( imr, nl, ntopfl, nfields2, nfields, &
                         ca, cc, fm_tl, fp_tl )

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

  real(8), intent(inout) :: fm_tl(imr,nl,nfields2) ! fields before diffusion 

! !OUTPUT PARAMETERS:

  real(8), intent(out) :: fp_tl(imr,nl,nfields2) ! fields after diffusion 

! !DESCRIPTION: This routine computes x=(I+K)**(-1) y, given vector y and 
!               tridiagonal matrix K.  I is the identity matrix. The elements
!               of K are the flux coefficients K * density * d/dz determined 
!               for the combination of surface drag, vertical diffusion, 
!               and pbl scheme. In the notation of the NCAR CCM3 description 
!               (Kiehl, et al. 1996, page 95):
!         ca(i,k)   = a_k at longitudinal point i
!         cc(i,k)   = c_k at longitudinal point i
!         b         = 1.+a_k +c_k  (the diagonal of I+K)
!         work(i,k) = b_k + c_k*E_(k-1)
!         fm_tl     = field on input, but then modified to be
!                     Fu_k * (b_k-c_k*E(k-1))
!   The rescursive excpressions for work and f_m are:
!         fm(k)= fieldin(k) + c_k*fm(k-1)/work(k-1)
!         work(k)=b_k - c(k)*a(k-1)/work(k-1)
!   This follows the reformulation of the CCM3 formulation of 
!   Richtmeyer and Morton (1967, page 198-201) as used in version 2 of
!   the Mesoscale Adjoint Modeling System (MAMS2) by R. M. Errico. In
!   addition to the formulation being different (the arrays work and 
!   fm being calculated rather than Fu and E), the surface momentum drag
!   ia also treated as linear and temporally implicit here, rather than 
!   as a forcing and temporally explicit as in CCM3.

! !REMARKS:

! !REVISION HISTORY:
!
!  10Mar2003  Errico  Initial code and implimentation
!
!EOP

  character(len=*), parameter :: myname_ = myname//'::tridiag_tl_' 
  integer   :: i             ! longitudinal index  
  integer   :: k             ! vertical index  
  integer   :: km1           ! k-1    
  integer   :: kp1           ! k+1              
  integer   :: nf            ! field indicator   
  real(8)  :: work(imr,nl)   ! b_k - c_kE_k-1 where E_k=a_k/work_k  
  real(8)  :: b              ! diagonal element of matrix  
  real(8)  :: alpha          
  

  do i=1,imr
    work(i,ntopfl)=1.d0+ca(i,ntopfl)
  enddo

  do k=ntopfl+1,nl
    km1=k-1
    do i=1,imr
      b=1.d0+ca(i,k)+cc(i,k)
      alpha=cc(i,k)/work(i,km1)
      work(i,k)=b-alpha*ca(i,km1)
      do nf=1,nfields
        fm_tl(i,k,nf)=fm_tl(i,k,nf)+alpha*fm_tl(i,km1,nf)
      enddo
    enddo 
  enddo

  do i=1,imr
    do nf=1,nfields
      fp_tl(i,nl,nf)=fm_tl(i,nl,nf)/work(i,nl)
    enddo
  enddo  

  do k=nl-1,ntopfl,-1
    kp1=k+1
    do i=1,imr
      do nf=1,nfields
        fp_tl(i,k,nf)=(fm_tl(i,k,nf)+ca(i,k)*fp_tl(i,kp1,nf))/work(i,k)
      enddo
    enddo
  enddo  


end subroutine tridiag_tl_


subroutine d2a2_tl1( u_tm, v_tm, ua_tl, va_tl, im, jm, jfirst, jlast, ng_d, ng_s, coslon, sinlon )
!******************************************************************
!******************************************************************
!** This routine was generated by Automatic differentiation.     **
!** FastOpt: Transformation of Algorithm in Fortran, TAF 1.5.2   **
!******************************************************************
!******************************************************************
!==============================================
! referencing used modules
!==============================================
use precision

!==============================================
! all entries are defined explicitly
!==============================================
implicit none

!==============================================
! declare arguments
!==============================================
integer :: im
real(kind=r8) :: coslon(im)
integer :: jfirst
integer :: jlast
integer :: jm
integer :: ng_d
integer :: ng_s
real(kind=r8) :: sinlon(im)
!real(kind=r8) :: u_tm(im,jfirst-ng_d:jlast+ng_s)
real(kind=r8) :: u_tm(im,jfirst:jlast+1)
real(kind=r8) :: ua_tl(im,jfirst:jlast)
real(kind=r8) :: v_tm(im,jfirst:jlast)
real(kind=r8) :: va_tl(im,jfirst:jlast)

!==============================================
! declare local variables
!==============================================
integer :: i
integer :: imh
integer :: j
integer :: jn2g0
integer :: js2g0
real(kind=r8) :: un
real(kind=r8) :: un_tl
real(kind=r8) :: us
real(kind=r8) :: us_tl
real(kind=r8) :: vn
real(kind=r8) :: vn_tl
real(kind=r8) :: vs
real(kind=r8) :: vs_tl

!----------------------------------------------
! TANGENT LINEAR AND FUNCTION STATEMENTS
!----------------------------------------------
imh = im/2
jn2g0 = min(jlast,jm-1)
js2g0 = max(jfirst,2)
do j = js2g0, jn2g0
  do i = 1, im
    ua_tl(i,j) = 0.5*u_tm(i,j+1)+0.5*u_tm(i,j)
  end do
end do
do j = js2g0, jn2g0
  do i = 1, im-1
    va_tl(i,j) = 0.5*v_tm(i+1,j)+0.5*v_tm(i,j)
  end do
  va_tl(im,j) = 0.5*v_tm(im,j)+0.5*v_tm(1,j)
end do
if (jfirst .eq. 1) then
  us_tl = 0.0d0
  vs_tl = 0.0d0
  do i = 1, imh
    us_tl = ua_tl(i+imh,2)*sinlon(i)-ua_tl(i,2)*sinlon(i)+us_tl-va_tl(i+imh,2)*coslon(i)+va_tl(i,2)*coslon(i)
    vs_tl = ua_tl(i+imh,2)*coslon(i)-ua_tl(i,2)*coslon(i)+va_tl(i+imh,2)*sinlon(i)-va_tl(i,2)*sinlon(i)+vs_tl
  end do
  us_tl = us_tl/dble(im)
  vs_tl = vs_tl/dble(im)
  do i = 1, imh
    ua_tl(i,1) = (-(us_tl*sinlon(i)))-vs_tl*coslon(i)
    va_tl(i,1) = us_tl*coslon(i)-vs_tl*sinlon(i)
    ua_tl(i+imh,1) = -ua_tl(i,1)
    va_tl(i+imh,1) = -va_tl(i,1)
  end do
endif
if (jlast .eq. jm) then
  un_tl = 0.0d0
  vn_tl = 0.0d0
  do i = 1, imh
    un_tl = ua_tl(i+imh,jm-1)*sinlon(i)-ua_tl(i,jm-1)*sinlon(i)+un_tl+va_tl(i+imh,jm-1)*coslon(i)-va_tl(i,jm-1)*coslon(i)
    vn_tl = (-(ua_tl(i+imh,jm-1)*coslon(i)))+ua_tl(i,jm-1)*coslon(i)+va_tl(i+imh,jm-1)*sinlon(i)-va_tl(i,jm-1)*sinlon(i)+vn_tl
  end do
  un_tl = un_tl/dble(im)
  vn_tl = vn_tl/dble(im)
  do i = 1, imh
    ua_tl(i,jm) = (-(un_tl*sinlon(i)))+vn_tl*coslon(i)
    va_tl(i,jm) = (-(un_tl*coslon(i)))-vn_tl*sinlon(i)
    ua_tl(i+imh,jm) = -ua_tl(i,jm)
    va_tl(i+imh,jm) = -va_tl(i,jm)
  end do
endif

end subroutine d2a2_tl1

subroutine d2a3d_tl1( u_tm, v_tm,  ua_tl,  va_tl, im, jm, km, jfirst, jlast, ng_d, ng_s, coslon, sinlon )
!******************************************************************
!******************************************************************
!** This routine was generated by Automatic differentiation.     **
!** FastOpt: Transformation of Algorithm in Fortran, TAF 1.5.2   **
!******************************************************************
!******************************************************************
!==============================================
! referencing used modules
!==============================================
use precision

!==============================================
! all entries are defined explicitly
!==============================================
implicit none

!==============================================
! declare arguments
!==============================================
integer :: im
real(kind=r8) :: coslon(im)
integer :: jfirst
integer :: jlast
integer :: jm
integer :: km
integer :: ng_d
integer :: ng_s
real(kind=r8) :: sinlon(im)
!real(kind=r8) :: u_tm(im,jfirst-ng_d:jlast+ng_s,km)
real(kind=r8) :: u_tm(im,jfirst:jlast+1,km)
real(kind=r8) :: ua_tl(im,jfirst:jlast,km)
real(kind=r8) :: v_tm(im,jfirst-ng_s:jlast+ng_d,km)
real(kind=r8) :: va_tl(im,jfirst:jlast,km)

!==============================================
! declare local variables
!==============================================
integer :: k

!----------------------------------------------
! TANGENT LINEAR AND FUNCTION STATEMENTS
!----------------------------------------------

#ifdef SPMD
call mp_barrier
call mp_send_s(im,jm,jfirst,jlast,1,km,0,1,u_tm)
call mp_barrier
call mp_recv_n(im,jm,jfirst,jlast,1,km,0,1,u_tm)
#endif  

#ifdef USE_OPENMP
!$omp parallel do private(k)
#endif /* not USE_OPENMP */
do k = 1, km
  call d2a2_tl1( u_tm(1:im,jfirst:jlast+1,k),v_tm(1,jfirst,k), &
    ua_tl(1,jfirst,k),va_tl(1,jfirst,k),im,jm,jfirst,jlast,ng_d,ng_s, &
    coslon,sinlon)
end do

end subroutine d2a3d_tl1

! 
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sponge_tl_ --- solver for implicit sponge
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
      subroutine hdiff_tl (imr,jnp,nl,jfirst,jlast,dt,pt_tl,u_tl,v_tl,ng_d,ng_s)

! !DESCRIPTION: Horizontal diffusion on eta surfaces using a 5 point
!               discret Laplacian, except at the poles where the 
!               Laplacian is computed as a polar cap (for T) or polar 
!               sector (for u) area-average of the flux divergence out of
!               the defined area, computed using the line integral of the 
!               noraml gradiant (Gauss's theorem). A time splittng 
!               algorithm is used: first the effect of meridional gradients
!               are determined, then that of the zonal gradients using 
!               the previous meridionally diffused values. For numerical 
!               stability of the latter calculations near the poles, where 
!               delta x and thus the dimensional zonal wavenumbers can 
!               become too large for the model time step, the algorithm 
!               uses an implicit time scheme.  This is done by computing
!               the diffusion in zonal spectral space, using FFTs for the
!               transforms.

      implicit none

      integer,  intent(in)  :: imr
      integer,  intent(in)  :: jnp
      integer,  intent(in)  :: nl
      integer,  intent(in)  :: jfirst
      integer,  intent(in)  :: jlast
      integer,  intent(in)  :: ng_d
      integer,  intent(in)  :: ng_s
      real(r8), intent(in)  :: dt 

      real(kind=r8), intent(inout) :: pt_tl(imr,jfirst-ng_d:jlast+ng_d,nl)
      real(kind=r8), intent(inout) :: u_tl(imr,jfirst-ng_d:jlast+ng_s,nl)
      real(kind=r8), intent(inout) :: v_tl(imr,jfirst-ng_s:jlast+ng_d,nl)
      
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
        call  hdiff_tend_tl (imr,jnp,jfirst,jlast,nl,cosp,cose, &
                          coslon,sinlon,dt,ifax,trigs,pt_tl,'S',ng_d,ng_d)
!
        call  hdiff_tend_tl (imr,jnp,jfirst,jlast,nl,cose,cosp, &
                          coslon,sinlon,dt,ifax,trigs,u_tl,'u',ng_d,ng_s)
!         
        call  hdiff_tend_tl (imr,jnp,jfirst,jlast,nl,cosp,cose, &
                          coslon,sinlon,dt,ifax,trigs,v_tl,'v',ng_s,ng_d)
!
! The pole values of v_tl have been are overwritten in hdiff_tend_tl. It is 
! assumed that they in fact are not referenced after this point, unless
! re-assigned in further time steps. Therefore they are formally set to
! zero, as the adjoint will assume.
!  
        if (jfirst==1) then 
          v_tl(:,  1,:)=0.d0
        endif
        if (jlast==jnp) then          
          v_tl(:,jnp,:)=0.d0
        endif
!
      end subroutine hdiff_tl 
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine hdiff_fac_tl (imr,dt,xnu,cosp,damp)
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
!
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
      end subroutine hdiff_fac_tl 
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
     subroutine hdiff_zonal_tl (imr,jfirst,jlast,j1,j2,nl, &
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
!
!     integer, parameter :: r8=8   ! use precision instead

! Input
     integer,  intent(in) :: imr  
     integer,  intent(in) :: jfirst,jlast,j1,j2
     integer,  intent(in) :: nl
     real(r8), intent(in) :: dt,xnu
     integer,  intent(in) :: ifax(13)        ! factorization of im
     real(r8), intent(in) :: trigs(3*imr/2+1) ! trig factors for FFT
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
!
     do j=j1,j2
!
       facx=(xnu*dt*imr**4.)/(pi*earthr*cosp(j))**4.     ! 'facx' includes (PI/2) adjustment
!
       if ( facx.ge.fcrit ) then       ! compute using discrete implicit scheme
!
! Determine spectral damping factors
!         
         call hdiff_fac_tl (imr,dt,xnu,cosp(j),damp)
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

        do k=1,nl
          do i=3,imrp2
            fs(i,k) = fs(i,k)*damp(i-2)
          enddo
        enddo

!
! Compute filtered fields from filtered coefficients

        call rfftmlt(fs, work, trigs, ifax, 1, imrp2, imr, nl, 1)
      
! Copy field into an array with two extra elements for each j
! as required for the FFT routine
        do k=1,nl
          f(1:imr,j,k)=fs(1:imr,k) 
        enddo
!
      endif  ! test on whether implicit spectral
    enddo    ! loop over j

 end subroutine hdiff_zonal_tl
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine hdiff_del2_tl (imr,jnp,jfirst,jlast,nl,cosp,cose,   &
                     coslon, sinlon, xnu, dt, f, del2f,name,ng_d1,ng_d2)
!____________________________________________________________________
!
!   ! DESCRIPTION: 
!       Compute explicit Laplacian of a field f
!
      implicit none

      integer,  intent(in) :: imr
      integer,  intent(in) :: jnp
      integer,  intent(in) :: jfirst
      integer,  intent(in) :: jlast
      integer,  intent(in) :: ng_d1
      integer,  intent(in) :: ng_d2
      integer,  intent(in) :: nl
      real(r8), intent(in) :: cosp(jnp) ! cosine at field points
      real(r8), intent(in) :: cose(jnp) ! cosine at edges between points
      real(r8), intent(in) :: coslon(imr)  !
      real(r8), intent(in) :: sinlon(imr)  !
      real(r8), intent(in) :: xnu, dt
! Note: this nomenclature for cosp,e differs from FVGCM: for
! u field points, it is the reverse of the FVGCM specs, which refers to 
! T field points. 
      character(len=*), intent(in) :: name  ! 'u', 'v', or 'S'
     
      real(r8), intent(inout) :: f(imr,jfirst-ng_d1:jlast+ng_d2,nl)
      real(r8), intent(out)   :: del2f(imr,jfirst-ng_d1:jlast+ng_d2,nl)

! local variables
      integer  :: j1, j2, ja, jb, jp, jm, jnume
      integer  :: i, j, k
      integer  :: nmax
      real(r8) :: afac, afacx, afacy
      real(r8) :: dlon_eq, dle, dlp, dlat, area, sine
      real(r8) :: fac, fac0, facm, facp, facx, fac1
      real(r8) :: sum_diff
      real(r8) :: cvr, cvi 
      real(r8) :: fcheck
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
!
!
      if ((name=='v') .and. ((j1>jfirst) .or. (j2<jlast)) ) then
!
! Special treatment for v at the poles. In order to take a gradiant between 
! v at the pole and at the adjacent latitude, values of v at the pole are
! specified as identical to the wavenumber 1 component of the v field at the 
! adjacent latitude.  This over-writes any pole values that may have been 
! specified for earlier calculations.
!
        if (j1 > jfirst) then ! input range of lats includes S.P.
          ja=1
          jb=2
          do k=1,nl
            cvr=0.d0
            cvi=0.d0
            do i=1,imr
              cvr=cvr+coslon(i)*f(i,jb,k)
              cvi=cvi+sinlon(i)*f(i,jb,k)
            enddo
            cvr=2.0d0*cvr/real(imr, r8)  ! this is  2* real part of wavenumber 1
            cvi=2.0d0*cvi/real(imr, r8)  ! this is -2* imag part of wavenumber 1
            do i=1,imr
              f(i,ja,k)=cvr*coslon(i)+cvi*sinlon(i)
            enddo
          enddo
!
        endif  ! test if next to S.P. included
!
        if (j2 < jlast) then ! input range of lats includes N.P.
                      
          ja=jnp
          jb=jnp-1
          do k=1,nl
            cvr=0.d0
            cvi=0.d0
            do i=1,imr
              cvr=cvr+coslon(i)*f(i,jb,k)
              cvi=cvi+sinlon(i)*f(i,jb,k)
            enddo
            cvr=2.0d0*cvr/real(imr, r8)  ! this is  2* real part of wavenumber 1
            cvi=2.0d0*cvi/real(imr, r8)  ! this is -2* imag part of wavenumber 1
            do i=1,imr
              f(i,ja,k)=cvr*coslon(i)+cvi*sinlon(i)
            enddo
          enddo
                 
        endif  ! test if next to N.P. included
!
      endif    ! test if v field next to poles included
!
      if ((name=='S') .and. ((j1>jfirst) .or. (j2<jlast)) ) then
!
! Special treatment for Scalar field at poles.  The Lplacian is computed 
! as the area mean of the divergence of the gradiant vector out of the 
! polar cap centered on the pole.  The assumption here is that the grid 
! spacing is symmetric with respect to the equator, so e.g. only cose(2) 
! need be ref.
!
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
              sum_diff=sum_diff+f(i,jb,k)-f(i,ja,k)    
            enddo
            del2f(1:imr,ja,k)=fac*sum_diff
          enddo 
        endif  ! test if S.P. included 
!
        if (j2 < jlast) then ! input range of lats includes N.P.
          ja=jnp
          jb=jnp-1
          do k=1,nl
            sum_diff=0.d0   
            do i=1,imr
              sum_diff=sum_diff+f(i,jb,k)-f(i,ja,k)    
            enddo
            del2f(1:imr,ja,k)=fac*sum_diff
          enddo 
        endif  ! test if N.P. included 
!
      endif    ! test if Scalar field at poles included
!
      if ((name=='u') .and. ((j1>jfirst) .or. (j2<jlast)) ) then
!
! Special treatment for u field at poles. Compute as 4-point Laplacian 
! on sector using Gauss's law:  The Lplacian is computed as the area 
! mean of the divergence of the gradiant vector out the three sides of 
! the sector whose vertex is at the pole. The assumption here is that 
! the grid spacing is symmetric with respect to the equator, so e.g. 
! only cose(2) need be ref.
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
            del2f(:,ja,k)=fac0*f(:,jb,k)+fac1*f(:,ja,k)
          enddo
        endif  ! test if next to S.P. included
!
        if (j2 < jlast) then ! input range of lats includes N.P.
          ja=jnp
          jb=jnp-1
          do k=1,nl
            del2f(:,ja,k)=fac0*f(:,jb,k)+fac1*f(:,ja,k)
          enddo
        endif  ! test if next to N.P. included
!
      endif    ! test if u field next to poles included
!
!
! Compute del**2 for points where pole is not an issue.
!
! Consider contribution by latitudinal gradients first
#ifdef SPMD
      call mp_barrier
      call mp_send3d_ns(imr,jnp,jfirst,jlast,1,nl,ng_d1,ng_d2,f,2)
      call mp_barrier
      call mp_recv3d_ns(imr,jnp,jfirst,jlast,1,nl,ng_d1,ng_d2,f,2)
#endif
      do j=j1,j2
        jp=j+1
        jm=j-1
        if (name == 'u') then
          facp=afacy*cose( j)/cosp(j)
          facm=afacy*cose(jm)/cosp(j)
        else
          facp=afacy*cose(jp)/cosp(j)
          facm=afacy*cose( j)/cosp(j)
        endif
        fac0=-(facp+facm)
         do k=1,nl 
          do i=1,imr
            del2f(i,j,k)=facp*f(i,jp,k)+fac0*f(i,j,k)+facm*f(i,jm,k) 
          enddo
        enddo
      enddo
!
! Consider contribution by longitudinal gradients 
! The change of j1 and j2 here for the u field now allows u at latitudes 
! adjacent to the poles to be treated.  But note that this treatment here
! has not been modified to account for the difference between the area of
! the polar sector for this u latitude and the area of u grid boxes not 
! adjacent to the pole.  
!
      if (name=='u') then
        j1=max(jfirst,2)
        j2=min(jlast,jnp)
      endif
!
      do j=j1,j2
        fcheck = (xnu*dt*imr**4.)/(pi*earthr*cosp(j))**4.

        if ( fcheck < fcrit )  then
! Do not apply near poles where delta x is too small and will lead to a numerical
! instability: imr/(pi*earthr*cosp)**4 is the factor you get if you apply 
! d**4/dx**4 to an f that is a 2 delta x wave using centered finite differences
!
          facx=afacx/cosp(j)**2  ! this is sqrt(xnu*dt)(1/delta x)**2 at this latitude
          fac0=-2.0d0*facx
          do k=1,nl 
            del2f(  1,j,k)=del2f(  1,j,k) &
                       +fac0*f(  1,j,k)+facx*(f(2,j,k)+f(  imr,j,k))
            del2f(imr,j,k)=del2f(imr,j,k) &
                       +fac0*f(imr,j,k)+facx*(f(1,j,k)+f(imr-1,j,k))
            do i=2,imr-1
              del2f(i,j,k)=del2f(i,j,k) &
                       +fac0*f(  i,j,k)+facx*(f(i+1,j,k)+f(i-1,j,k))
            enddo
          enddo
        endif   ! for fcheck
      enddo 
!
      end subroutine hdiff_del2_tl
!
!  x x x x x x x x x x   x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine hdiff_tend_tl (imr,jnp,jfirst,jlast,nl,cosp,cose, &
                             coslon,sinlon,dt,ifax,trigs,f,name,ng_d1,ng_d2)
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
                                                                                                                             
      implicit none
          
      integer,  intent(in) :: imr
      integer,  intent(in) :: jnp
      integer,  intent(in) :: jfirst
      integer,  intent(in) :: jlast
      integer,  intent(in) :: ng_d1
      integer,  intent(in) :: ng_d2
      integer,  intent(in) :: nl
      integer,  intent(in) :: ifax(13)        ! factorization of im
      real(r8), intent(in) :: trigs(3*imr/2+1) ! trig factors for FFT
      real(r8), intent(in) :: cosp(jnp) ! cosine at field points
      real(r8), intent(in) :: cose(jnp) ! cosine at edges between points
      real(r8), intent(in) :: coslon(imr),sinlon(imr) ! only used for v 
      real(r8), intent(in) :: dt        ! step for time marching

! Note: this nomenclature for cosp,e differs from FVGCM: for
! u field points, it is the reverse of the FVGCM specs, which refers to
! T field points.
      character(len=*), intent(in) :: name  ! 'u', 'v', or 'S'
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
!    "Close" criteria is specified within routine 'hdiff_zonal_tl' itself 
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
      work(:,:,:)=0.d0
      call hdiff_del2_tl(imr,jnp,jfirst,jlast,nl,cosp,cose, &
                      coslon,sinlon,xnu,dt,f,work,name,ng_d1,ng_d2)
!
      del4f(:,:,:)=0.d0
      call hdiff_del2_tl(imr,jnp,jfirst,jlast,nl,cosp,cose, &
                      coslon,sinlon,xnu,dt,work,del4f, name,ng_d1,ng_d2)
!
! Note: these j1 and j2 are not the same as in hdiff_del2_tl
      j1=max(jfirst,2)
      if (name=='u') then
        j2=min(jlast,jnp)
      else
        j2=min(jlast,jnp-1)
      endif

      call hdiff_zonal_tl(imr,jfirst,jlast,j1,j2,nl, &
                            dt,xnu,cosp(j1:j2),ifax,trigs,f(:,jfirst:jlast,:))

      f(:,:,:) = f(:,:,:) - del4f(:,:,:)
!
      end subroutine hdiff_tend_tl
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
   end module m_physdrv1_tl
