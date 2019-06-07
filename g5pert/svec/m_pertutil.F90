!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_pertutil --- Utility programs to handle perturbations
!
! !INTERFACE:

      module m_pertutil

      use precision
      
      use m_const, only : Cp     => cpm
      use m_const, only : R      => rgas
      use m_const, only : kappa
      use m_const, only : zvir
      use m_const, only : pstd
      use m_const, only : tstd
      use m_const, only : alhl       ! latent heat condensation

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type

      use m_mapz, only : set_eta

      use m_dyn
      use m_maph_pert

      use m_stdio
      use m_inpak90
      use m_die, only: MP_die, die

      use m_daInterp, only : daInterp_vatod
      use m_daInterp, only : daInterp_vdtoa
      use m_daInterp, only : dainterp_vdtoa_ad
      use m_daInterp, only : daInterp_vatod_ad

      use m_ioutil, only : luavail
      use m_chars, only : lowercase

      implicit NONE
 
! !PUBLIC MEMBER FUNCTIONS:
 
      PRIVATE

      PUBLIC pertutil_find_name
      PUBLIC pertutil_calculate_ps
      PUBLIC pertutil_calculate_ps_ad
      PUBLIC pertutil_transform_q2r
      PUBLIC pertutil_transform_q2r_ad
      PUBLIC pertutil_transform_T
      PUBLIC pertutil_transform_T_tl
      PUBLIC pertutil_transform_T_ad
      PUBLIC pertutil_transform_T_adinv
      PUBLIC pertutil_enorm
      PUBLIC pertutil_inverse_grad_enorm
      PUBLIC pertutil_horiz_grid
      PUBLIC pertutil_vert_grid
      PUBLIC pertutil_pt2vpt
      PUBLIC pertutil_ptv2pt
      PUBLIC pertutil_pt2t
      PUBLIC pertutil_pt2t_ad
      PUBLIC pertutil_pt2t_tl
      PUBLIC pertutil_t2pt
      PUBLIC pertutil_t2pt_ad
      PUBLIC pertutil_t2pt_tl
      PUBLIC pertutil_pkappa
      PUBLIC pertutil_read_data
      PUBLIC pertutil_write_data
      PUBLIC pertutil_setparam
      PUBLIC pertutil_getparam

!     PUBLIC g5tog4
!     PUBLIC g4tog5_pert

      PUBLIC e_cross_product
      PUBLIC comp_kdomain

      interface pertutil_transform_T; module procedure &
                         transform_T0_, &
                         transform_T2_
      end interface

      interface pertutil_calculate_ps;       module procedure calculate_ps_;       end interface
      interface pertutil_transform_q2r;      module procedure transform_q2r_;      end interface
      interface pertutil_transform_q2r_ad;   module procedure transform_q2r_ad_;   end interface
      interface pertutil_enorm;              module procedure enorm_;              end interface
      interface pertutil_inverse_grad_enorm; module procedure inverse_grad_enorm_; end interface
      interface pertutil_find_name;          module procedure find_name_;          end interface
      interface pertutil_horiz_grid;         module procedure horiz_grid_;         end interface
      interface pertutil_vert_grid;          module procedure vert_grid_;          end interface
      interface pertutil_pt2vpt;             module procedure pt2vpt_;             end interface
      interface pertutil_ptv2pt;             module procedure ptv2pt_;             end interface
      interface pertutil_pt2t;               module procedure pt2t_;               end interface
      interface pertutil_t2pt;               module procedure t2pt_;               end interface
      interface pertutil_pkappa;             module procedure pkappa_;             end interface
      interface pertutil_read_data;          module procedure read_data_;          end interface
      interface pertutil_write_data;         module procedure write_data_;         end interface

      interface pertutil_transform_T_tl;     module procedure &
                         transform_T1_tl_
      end interface
      interface pertutil_transform_T_ad;     module procedure &
                         transform_T0_ad_, &
                         transform_T1_ad_
      end interface
      interface pertutil_transform_T_adinv; module procedure &
                         transform_T0_adinv_
      end interface
      interface pertutil_setparam;           module procedure &
                         setparamI_, &
                         setparamL_, &
                         setparamR_, &
                         setparamC_
      end interface
      interface pertutil_getparam;           module procedure &
                         getparamI_, &
                         getparamL_, &
                         getparamR_, &
                         getparamC_
      end interface

      interface pertutil_calculate_ps_ad;    module procedure calculate_ps_ad_;    end interface
      interface pertutil_pt2t_ad;            module procedure pt2t_ad_;            end interface
      interface pertutil_t2pt_ad;            module procedure t2pt_ad_;            end interface

      interface pertutil_pt2t_tl;            module procedure pt2t_tl_;            end interface
      interface pertutil_t2pt_tl;            module procedure t2pt_tl_;            end interface

      interface g5tog4;                      module procedure g5tog4_;             end interface
      interface g4tog5_pert;                 module procedure g4tog5_pert_;        end interface

      interface e_cross_product;             module procedure e_cross_product_;    end interface
      interface comp_kdomain;                module procedure comp_kdomain_;       end interface
!
! !DESCRIPTION: Energy unscaling of sensitivity fields 
!
! !REVISION HISTORY:
!
!  05Jan2008 Todling  Trying to organize this mess
!  16Mar2012 Todling  Add wgrid_
!  15Oct2015 Holdaway Support for moist available enthalpy calculation
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'm_pertutil'
      real(r8), parameter :: pref = 100.d0*pstd
      real(r8), parameter :: tref = tstd

      integer,          save :: vectype_ = 4     ! Default i/o vectors are GEOS-4 compliant
      character(len=1), save :: wgrid_   = 'd'   ! Default grid for winds: d-grid
      character(len=3), save :: tname_   = 'vpT' ! Default temperature in file
      real(r8),         save :: eps_eer_ = 1.0d0 ! Ehrendorfer, Errico, and Raedder q-energy coefficient
      integer,          save :: vnorm_   = 0     ! 0=mass weight energy norm
                                                 ! 1=height weighted
                                                 ! 2=hybrid mean mass+height weighted

      CONTAINS

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: calculate_ps --- calculate surface pressure as vertical sum of delp
!
!
! !INTERFACE:
!
      subroutine calculate_ps_ (im1, jn1, nl1, ptop, delp, ps)

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: im1,jn1,nl1
      real(r8), intent(in)  :: ptop
      real(r8), intent(in)  :: delp(im1,jn1,nl1) ! p-thickness of layers

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: ps(im1,jn1)      ! surface pressure
     
! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!
! calculate surface pressure as vertical sum of delp 
!    

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  27Jul2005  Elena N.   Added ptop to surface pressure value
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer  :: i,j,k
   
      do j=1,jn1
        do i=1,im1
          ps(i,j)=ptop
          do k=1,nl1 
            ps(i,j)=ps(i,j)+delp(i,j,k)
          enddo
        enddo
      enddo

      end subroutine calculate_ps_ 

!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: calculate_ps_ad --- Compute adjoint ps field from adjoint delp field.
!
!
! !INTERFACE:
!
      subroutine calculate_ps_ad_ (imr, jnp, nl, delp, delbk, ps)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: imr,jnp,nl
      real(r8), intent(in)  :: delbk(nl)
      real(r8), intent(in)  :: delp(imr,jnp,nl) ! dJ/d(p-thickness layers)

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: ps(imr,jnp)      ! dJ/d(surface pressure)

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!
! Compute adjoint ps field from adjoint delp field.
! This computation is based on an assumed relation for calculating
!  initial conditions of delp from ps.  That relationship is
!        delp(k) = deltab(k) * ps
!  where deltab is the thickness in b, the coefficient that multiplies ps
!  to define the p values at interfaces in the FVGCM eta coordinate.
!  Note that near the model top, where a pure p-coordinate is used, b=0
!  so delbk=0.  That is, changes to ps do not change initial values of
!  delp in that region, so the adjoint of delp there does not affect ps.
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  08Aug2005  Elena N.   Adopted this routine from m_postStats.F90 
!                        to work in the fsens2pert.x stand-alone program
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer  :: i,j,k
   
      do j=1,jnp
        do i=1,imr
          ps(i,j)=0.d0
          do k=1,nl 
            ps(i,j)=ps(i,j)+delbk(k)*delp(i,j,k)
          enddo
        enddo
      enddo

      end subroutine calculate_ps_ad_ 

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: transform_q2r_ --- transform specific humidity to water vapour mixing ratio
!
!
! !INTERFACE:
!
      subroutine transform_q2r_ (im1,jn1,nl1,nfields,fnames,fields)

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: im1,jn1,nl1,nfields

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

      character(len=*), intent(inout) :: fnames(nfields) 
      real(r8), intent(inout) :: fields(im1,jn1,nl1,nfields) 

! !DESCRIPTION:
!
! transform specific humidity to water vapour mixing ratio 
!    

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  15Oct2015  D. Holdaway  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer  :: i,j,k
      integer :: indexQ
   
      call find_name_ (nfields,fnames,'q__',indexQ)

      fields(:,:,:,indexQ) = fields(:,:,:,indexQ) / (1 - fields(:,:,:,indexQ))

      fnames(indexQ) = 'r__'    ! now converted to r


      end subroutine transform_q2r_ 

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: transform_q2r_ad_ --- adjoint of transform specific humidity to water vapour mixing ratio
!
!
! !INTERFACE:
!
      subroutine transform_q2r_ad_ (im1,jn1,nl1,nfields,fnames,fields_ad,fields_rf)

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: im1,jn1,nl1,nfields

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

      character(len=*), intent(inout) :: fnames(nfields) 
      real(r8), intent(inout) :: fields_ad(im1,jn1,nl1,nfields,1) 
      real(r8), intent(inout) :: fields_rf(im1,jn1,nl1,nfields,1) 

! !DESCRIPTION:
!
! adjoint of transform specific humidity to water vapour mixing ratio 
!    

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  15Oct2015  D. Holdaway  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer  :: i,j,k
      integer :: indexR
   
      call find_name_ (nfields,fnames,'r__',indexR)

      !Transform the reference back to q
      fields_rf(:,:,:,indexR,1) = fields_rf(:,:,:,indexR,1) / (1 + fields_rf(:,:,:,indexR,1))

      !Adjoint of transform q to r (actually takes r_ad to q_ad ready for pert model)
      fields_ad(:,:,:,indexR,1) = ( fields_rf(:,:,:,indexR,1)/(1-fields_rf(:,:,:,indexR,1))**2 &
                                  + 1/(1-fields_rf(:,:,:,indexR,1)) ) * fields_ad(:,:,:,indexR,1)

      fnames(indexR) = 'q__'    ! now converted to q

      print *, 'Input ad-moisture field reset to  ', fnames(indexR)

      end subroutine transform_q2r_ad_ 


!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: enorm --- calculates energy norm
!
!
! !INTERFACE:
!
      subroutine enorm_ (im1,jn1,nl1,nfields,                 &
                   jweights,ptop,                             &
                   fields,delp_ref,psfield,fnames,            &
                   esum) 

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1            ! grid dimensions
      integer, intent(in) :: nfields                ! number of fields 
      character(len=*), intent(in) :: fnames(nfields) 
      real(r8), intent(in)  :: jweights(jn1,2)      ! dlat(radians)*cos(lat)/2
      real(r8), intent(in)  :: ptop
      real(r8), intent(in)  :: fields(im1,jn1,nl1,nfields)  ! ref. or diff fields
      real(r8), intent(in)  :: delp_ref(im1,jn1,nl1)        ! ref. delta p 
      real(r8), intent(in)  :: psfield(im1,jn1)             ! ref. or diff ps

! !OUTPUT PARAMETERS:
 
      real(r8), intent(out) :: esum

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
  
! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  11Apr2005  Elena N.   Corrected a bug with indexing
!                        replaced  'do i=i,im1' to be  'do i=1,im1' in KEV looping 
!                        Corrected 0.5 factors for sum(jweights)=1
!  09Jan2008  Todling    Add moisture to energy norm
!  11Sep2013  Todling    Update to hande A-grid winds
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      character(len=*), parameter :: myname_ = myname //'*enorm_'

      integer  :: i,j,k,kk,nlp1,ierr
      integer  :: idu,idv,idt,idq
      integer  :: lunit1   
      real(r8) :: wfactor
      real(r8) :: dbleimr
      real(r8) :: tfac, ufac, pfac, qfac
      real(r8) :: eps_eer
      real(r8) :: dsigma(im1,jn1,nl1)  ! sigma-thickness of layer
      real(r8) :: ps_ref(im1,jn1)     ! reference surface pressure
      real(r8) :: jfac(jn1)           ! fractional area of globe in T grid box 
      real(r8) :: fa(im1,jn1) 
      real(r8) :: fb(im1,jn1)         ! (interpolated) squared values of field at T points
      real(r8) :: fb_ad(im1,jn1)      ! d(E)/d(fb) without ufac,pfac,tfac 
      real(r8) :: eu,ev,et,ep,eq,e(5)
      real(r8), allocatable :: pe(:,:,:)
      character(len=1) wgrid
      integer vnorm
      integer kz(1)

! Compute fractional areas
      dbleimr=1.d0/dble(im1)
      do j=1,jn1
        jfac(j)=jweights(j,2)*dbleimr  ! fractional area of T grid points
      enddo

! Assumption for division by area where sum(jweights)=1
      call pertutil_getparam ( 'vnorm', vnorm )
      call pertutil_getparam ( 'eps_eer', eps_eer )
      tfac=0.5d0*Cp/tref
      ufac=0.5d0
      pfac=0.5d0*R*tref/pref**2
      qfac=0.5d0*eps_eer*alhl*alhl/(Cp*tref)
      e=0.d0

      call pertutil_getparam ( 'wgrid', wgrid )
      print *, 'Energy operator treating winds as on ',wgrid,'-grid'

      call find_name_ (nfields,fnames,'u__',idu)      
      call find_name_ (nfields,fnames,'v__',idv)      
      call find_name_ (nfields,fnames,'T__',idt)  ! it is assumed that vpt->T has already happened 
      call find_name_ (nfields,fnames,'q__',idq)

      if (idt<=0) call MP_die (myname_,'incorrect temperature index',99)
!
! compute 3-D varying dsigma

      call calculate_ps_ (im1,jn1,nl1,ptop,delp_ref,ps_ref)

      if ( vnorm>0 ) then

          print *, 'Using V-norm'
          nlp1 = nl1+1
          allocate( pe(im1,jn1,nlp1), stat=ierr )
              if(ierr/=0) call MP_die(myname_,'Alloc(pe)',ierr)

          pe(:,:,1) = ptop
          do k=2,nlp1
             pe(:,:,k) = pe(:,:,k-1) + delp_ref(:,:,k-1)
          enddo
          do k=1,nl1
             dsigma(:,:,k) = log(pe(:,:,k+1)/pe(:,:,k))/log(pe(:,:,nlp1)/pe(:,:,1))
          enddo
          if(vnorm>1) then
             pe=0.0 ! pe used as auxiliar array now
             do k=1,nl1
                pe(:,:,k)=delp_ref(:,:,k)/ps_ref(:,:)
             enddo
             wfactor=0.5
             do k=1,nl1
                dsigma(:,:,k)=wfactor*(pe(:,:,k)+dsigma(:,:,k))
             enddo
          endif
          deallocate( pe, stat=ierr )
              if(ierr/=0) call MP_die(myname_,'DeAlloc(pe)',ierr)

      else

         do k=1,nl1
           dsigma(:,:,k)=delp_ref(:,:,k)/ps_ref(:,:)
         enddo    

      endif

!
!  Calculate total energy by looping over levels, then fields
!
      do k=1,nl1

        do j=1,jn1
          do i=1,im1
            fb_ad(i,j)=dsigma(i,j,k)*jfac(j)
          enddo
        enddo
!     
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              fa(i,j)=fields(i,j,k,idu)
            enddo
          enddo

          if (wgrid=='d') then
             do i=1,im1
               fb(i,  1)   =fa(i,  2)*fa(i,  2)
               fb(i,jn1)   =fa(i,jn1)*fa(i,jn1)
               fb(i,2)     =0.5*(fa(i,2)*fa(i,2) + fa(i,3)*fa(i,3))
               do j=3,jn1-1     
                 fb(i,j)   =0.5*(fa(i,j)*fa(i,j) + fa(i,j+1)*fa(i,j+1))
               enddo
             enddo
             eu=0.d0
             do j=1,jn1
               do i=1,im1
                 eu=eu+fb_ad(i,j)*fb(i,j)
               enddo
             enddo
          else
             eu=0.d0 
             do j=1,jn1
               do i=1,im1
                 eu=eu+fb_ad(i,j)*fa(i,j)*fa(i,j)
               enddo
             enddo
          endif
          
          e(1)=e(1)+ufac*eu
 

!  Energy calculations for v field e(2)
!  Note that fb grid is assumed to be offset to east of v grid point
          fa = 0.0
          do i=1,im1
            do j=2,jn1-1
              fa(i,j)=fields(i,j,k,idv)
            enddo
          enddo

          if (wgrid=='d') then
             do j=2,jn1-1
               fb(im1,j)=0.5*(fa(1,j)*fa(1,j) + fa(im1,j)*fa(im1,j))
               do i=1,im1-1  
                 fb(i,j)=0.5*(fa(i,j)*fa(i,j) + fa(i+1,j)*fa(i+1,j))
               enddo  
             enddo
             do i=1,im1
               fb(i,  1)=fa(i,    2)*fa(i,    2)
               fb(i,jn1)=fa(i,jn1-1)*fa(i,jn1-1)
             enddo
 
             ev=0.d0
             do j=1,jn1
               do i=1,im1
                 ev=ev+fb_ad(i,j)*fb(i,j)
               enddo
             enddo
          else
             ev=0.d0 
             do j=1,jn1
               do i=1,im1
                 ev=ev+fb_ad(i,j)*fa(i,j)*fa(i,j)
               enddo
             enddo
          endif
          e(2)=e(2)+ufac*ev

!     
!  Energy calculations for T field e(3)
!  T form is assumed
!     
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              fa(i,j)=fields(i,j,k,idt)
            enddo
          enddo

          et=0.d0 
          do j=1,jn1
            do i=1,im1
              et=et+fb_ad(i,j)*fa(i,j)*fa(i,j)
            enddo
          enddo
          e(3)=e(3)+tfac*et  

!     
!  Energy calculations for Q field e(5)
!     
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              fa(i,j)=fields(i,j,k,idq)
            enddo
          enddo

          eq=0.d0 
          do j=1,jn1
            do i=1,im1
              eq=eq+fb_ad(i,j)*fa(i,j)*fa(i,j)
            enddo
          enddo
          e(5)=e(5)+qfac*eq  

      enddo ! loop over k
 
!
! Energy calculations for ps field e(4)
!
        fa = 0.0
        kk = 1
        do i=1,im1
          do j=1,jn1
            fa(i,j)=psfield(i,j)
          enddo
        enddo

        ep=0.d0
        do j=1,jn1
          do i=1,im1
            ep=ep+jfac(j)*fa(i,j)*fa(i,j)
          enddo
        enddo
        e(4)=pfac*ep

      esum=sum(e)
      print *,' ENERGY NORM '
      print *,' TOTAL ENERGY   ',esum      
      print *,' ENERGY BY COMPONENT (u,v,T,ps,q)  ',e

      lunit1= luavail()
      open(lunit1, file='Jnormi.txt', form='formatted')
      write (lunit1, '(6f20.12)') esum, e(1:5)
      close(lunit1)
      print *
      print *,' Wrote sensitivity norm into a file `Jnormi.txt`'
      print *

      end subroutine enorm_
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: inverse_grad_enorm --- computes perturbation fields from ADM/SENS run result 
!
!
! !INTERFACE:
!
      subroutine inverse_grad_enorm_ (im1,jn1,nl1,nfields,bk, &
                   jweights,ptop,                             &
                   fields,delp_ref,psfield,fnames,            &
                   fields_ad,psfield_ad) 

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1            ! grid dimensions
      real(r8), intent(in)  :: bk(nl1+1) 
      integer, intent(in) :: nfields                ! number of fields 
      character(len=*), intent(in) :: fnames(nfields) 
      real(r8), intent(in)  :: jweights(jn1,2)      ! dlat(radians)*cos(lat)/2
      real(r8), intent(in)  :: ptop
      real(r8), intent(in)  :: fields(im1,jn1,nl1,nfields)  ! ref. or diff fields
      real(r8), intent(in)  :: delp_ref(im1,jn1,nl1)        ! ref. delta p 
      real(r8), intent(in)  :: psfield(im1,jn1)            ! ref. or diff ps

! !OUTPUT PARAMETERS:
 
      real(r8), intent(out) :: fields_ad(im1,jn1,nl1,nfields)  ! adjoint u,v,T  
      real(r8), intent(out) :: psfield_ad(im1,jn1)            ! adjoint ps

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
  
! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  11Apr2005  Elena N.   Corrected a bug with indexing
!                        replaced  'do i=i,im1' to be  'do i=1,im1' in KEV looping 
!                        Corrected 0.5 factors for sum(jweights)=1
!  11Sep2013  Todling    Update to hande A-grid winds
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      character(len=*), parameter :: myname_ = myname //'*inverse_grad_enorm_'

      integer  :: i,j,k,kk,nlp1,ierr
      integer  :: idu,idv,idt,idp,idq
      integer  :: number_init         ! number of fields initialized
      real(r8) :: wfactor
      real(r8) :: dbleimr
      real(r8) :: tfac, ufac, pfac, qfac
      real(r8) :: eps_eer
      real(r8) :: dsigma(im1,jn1,nl1)  ! sigma-thickness of layer
      real(r8) :: ps_ref(im1,jn1)     ! reference surface pressure
      real(r8) :: jfac(jn1)           ! fractional area of globe in T grid box 
      real(r8) :: fa(im1,jn1) 
      real(r8) :: fb(im1,jn1)         ! (interpolated) squared values of field at T points
      real(r8) :: fa_ad(im1,jn1)      ! d(E)/d(fa) without ufac,pfac,tfac 
      real(r8) :: fb_ad(im1,jn1)      ! d(E)/d(fb) without ufac,pfac,tfac 
      real(r8), allocatable :: pe(:,:,:)
      character(len=1) wgrid
      integer vnorm
      integer kz(1)


! Compute fractional areas
      dbleimr=1.d0/dble(im1)
      do j=1,jn1
        jfac(j)=jweights(j,2)*dbleimr  ! fractional area of T grid points
      enddo

! Assumption for division by area where sum(jweights)=1
      call pertutil_getparam ( 'vnorm', vnorm )
      call pertutil_getparam ( 'eps_eer', eps_eer )
      ufac=0.5d0
      tfac=0.5d0*Cp/tref
      pfac=0.5d0*R*tref/pref**2
      qfac=0.5d0*eps_eer*alhl*alhl/(Cp*tref)

      ufac = 1.d0/ufac
      tfac = 1.d0/tfac
      pfac = 1.d0/pfac
      qfac = 1.d0/qfac

      call pertutil_getparam ( 'wgrid', wgrid )
      print *, 'Inverse energy operator treating winds as on ',wgrid,'-grid'

      call find_name_ (nfields,fnames,'u__',idu)      
      call find_name_ (nfields,fnames,'v__',idv)      
      call find_name_ (nfields,fnames,'T__',idt)     ! it is assumed that vpt->T has already happened 
      call find_name_ (nfields,fnames,'q__',idq)
!
! compute 3-D varying dsigma

      call calculate_ps_ (im1,jn1,nl1,ptop,delp_ref,ps_ref)

      if ( vnorm>0 ) then

          print *, 'Using V-norm'
          nlp1 = nl1+1
          allocate( pe(im1,jn1,nlp1), stat=ierr )
              if(ierr/=0) call MP_die(myname_,'Alloc(pe)',ierr)

          pe(:,:,1) = ptop
          do k=2,nlp1
             pe(:,:,k) = pe(:,:,k-1) + delp_ref(:,:,k-1)
          enddo
          do k=1,nl1
             dsigma(:,:,k) = log(pe(:,:,k+1)/pe(:,:,k))/log(pe(:,:,nlp1)/pe(:,:,1))
          enddo
          if(vnorm>1) then
             pe=0.0 ! pe used as auxiliar array now
             do k=1,nl1
                pe(:,:,k)=delp_ref(:,:,k)/ps_ref(:,:)
             enddo
             wfactor=0.5
             do k=1,nl1
                dsigma(:,:,k)=wfactor*(pe(:,:,k)+dsigma(:,:,k))
             enddo
          endif
          deallocate( pe, stat=ierr )
              if(ierr/=0) call MP_die(myname_,'DeAlloc(pe)',ierr)

      else

         do k=1,nl1
           dsigma(:,:,k)=delp_ref(:,:,k)/ps_ref(:,:)
         enddo    

      endif
!
!  Calculate total energy by looping over levels, then fields
!
      do k=1, nl1

        do j=1,jn1
          do i=1,im1
            fb_ad(i,j)=dsigma(i,j,k)*jfac(j)
          enddo
        enddo
!     
!  Energy calculations for u and v fields e(1) and e(2) respectively
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              fa(i,j)=fields(i,j,k,idu)
            enddo
          enddo

          if (wgrid=='d') then
             do i=1,im1
               fa_ad(i,  2)=fa(i,  2)/( 2.0d0*fb_ad(i,  1)+fb_ad(i,    2) )
               fa_ad(i,jn1)=fa(i,jn1)/( 2.0d0*fb_ad(i,jn1)+fb_ad(i,jn1-1) )
               do j=3,jn1-1     
                 fa_ad(i,j)=fa(i,j)/( fb_ad(i,j)+fb_ad(i,j-1) )
               enddo
             enddo
          else
             do j=1,jn1
               do i=1,im1
                 fa_ad(i,j)=0.5d0*fa(i,j)/fb_ad(i,j)
               enddo
             enddo
          endif

          fa_ad(:,1)=0.d0
          fields_ad(:,:,k,idu)=ufac*fa_ad(:,:)

!  Energy calculations for v field e(2)
!  Note that fb grid is assumed to be offset to east of v grid point
          fa = 0.0
          number_init=number_init+1       

          if (wgrid=='d') then
             do i=1,im1
               do j=2,jn1-1
                 fa(i,j)=fields(i,j,k,idv)
               enddo
             enddo
             do j=3,jn1-2
               fa_ad(1,j)=fa(1,j)/( fb_ad(im1,j)+fb_ad(1,j) )
               do i=1,im1-1
                 fa_ad(i+1,j)=fa(i+1,j)/( fb_ad(i,j)+ fb_ad(i+1,j) )
               enddo
             enddo
             fa_ad(1,    2)=fa(1,    2)/   &
                        ( fb_ad(im1,    2)+fb_ad(1,    2)+2.0*fb_ad(1,  1))
             fa_ad(1,jn1-1)=fa(1,jn1-1)/   &
                        ( fb_ad(im1,jn1-1)+fb_ad(1,jn1-1)+2.0*fb_ad(1,jn1))
             do i=2,im1
               fa_ad(i,    2)=fa(i,    2)/ &
                        ( fb_ad(i-1,    2)+fb_ad(i,    2)+2.0*fb_ad(i,  1))
               fa_ad(i,jn1-1)=fa(i,jn1-1)/ &
                        ( fb_ad(i-1,jn1-1)+fb_ad(i,jn1-1)+2.0*fb_ad(i,jn1))
             enddo
             fa_ad(:,  1)=0.d0
             fa_ad(:,jn1)=0.d0
          else
             do i=1,im1
               do j=1,jn1
                 fa(i,j)=fields(i,j,k,idv)
               enddo
             enddo
             do j=1,jn1
               do i=1,im1
                 fa_ad(i,j)=0.5d0*fa(i,j)/fb_ad(i,j)
               enddo
             enddo
          endif
          fields_ad(:,:,k,idv)=ufac*fa_ad(:,:)

!  Energy calculations for T field e(3)
!  T form is assumed
!     
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              fa(i,j)=fields(i,j,k,idt)
            enddo
          enddo

          do j=1,jn1
            do i=1,im1
              fa_ad(i,j)=0.5d0*fa(i,j)/fb_ad(i,j)
            enddo
          enddo

          fields_ad(:,:,k,idT)=tfac*fa_ad(:,:)

!  Energy calculations for Q field e(5)
!     
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              fa(i,j)=fields(i,j,k,idq)
            enddo
          enddo

          do j=1,jn1
            do i=1,im1
              fa_ad(i,j)=0.5d0*fa(i,j)/fb_ad(i,j)
            enddo
          enddo

          fields_ad(:,:,k,idq)=qfac*fa_ad(:,:)

      enddo ! loop over k

!
! Energy calculations for ps field e(4)
!
! Compute adjoint ps field from adjoint delp field
! Note: for global computation, where LPO is not applied,
!       difference (bk(nl+1) - bk(1)) is = 1.d0. Therefore,
!       this denominator is omitted here.

      fa = 0.0

        do i=1,im1
          do j=1,jn1
            fa(i,j)=psfield(i,j) 
          enddo
        enddo

        do j=1,jn1
          do i=1,im1
            fa_ad(i,j)=0.5d0*fa(i,j)/jfac(j)
          enddo
        enddo

        psfield_ad(:,:)=pfac*fa_ad(:,:)


      call find_name_ (nfields,fnames,'dp_',idp)
!
!  Update delp in output vector by distributing its updated ps through levels
!  Note: since (bk(nl1+1) - bk(1))= 1.d0, this denominator is omitted
!
      do k=1, nl1
         fields_ad(:,:,k,idp) = psfield_ad(:,:)*(bk(k+1) - bk(k))
         fields_ad(:,  1,k,idp) = sum(fields_ad(:,  1,k,idp))/dble(im1)
         fields_ad(:,jn1,k,idp) = sum(fields_ad(:,jn1,k,idp))/dble(im1)
      enddo

      end subroutine inverse_grad_enorm_
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: find_name --- 
!
!
! !INTERFACE:
!
      subroutine find_name_ (nums,all_names,name,index1)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)  :: nums     ! number of names to examine               
      character(len=*), intent(in)  :: all_names(nums)
      character(len=*), intent(in)  :: name

! !OUTPUT PARAMETERS:

      integer, intent(out) :: index1   ! if > 0, index of unique name found
                                       ! if = 0, no name found
                                       ! if < 0, - of number of names matching
     
! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
  
! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: num, nchar, nchars, indexa
      character(len=1) :: cchar
      logical :: test(nums)
       
      test(:)=.true.      
      nchars=len(name)
      
      do nchar=1,nchars
        cchar=name(nchar:nchar)
        if (cchar /= '*') then 
          do num=1,nums
            if (all_names(num)(nchar:nchar) /= cchar ) then 
              test(num)=.false.
            endif
          enddo
        endif
      enddo

      index1=0
      do num=1,nums
        if (test(num)) then
          if (index1 == 0) then 
            index1=num
          else if (index1 > 0) then
            index1=-2
          else
            index1=index1-1
          endif
        endif 
      enddo

      if (index1 < 1) then
        print *,' ***  WARNING **** '
        if (index1 == 0) then
          print *,' name = ',name,' not found in list = ',all_names
        else
          indexa=abs(index1)
          print *,indexa,' occurances of name = ',name,  &
                  ' found in list = ',all_names
        endif
      endif

      end subroutine find_name_
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: horiz_grid_ --- Determine some D-grid information required
!
!
! !INTERFACE:
!
      subroutine horiz_grid_ (im1, jn1, jweights, glats)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: im1, jn1    

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: jweights(jn1,2) ! area weights (1) u, (2) T   
      real(r8), intent(out) :: glats(jn1,2)    ! degrees lats (1) u, (2) T   

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!
!  Determine some D-grid information required
!
!  i=1 on v field corresponds to long=0
!  i=1 on T,u fields corresponds to long=0+360/(2*im1)
!  v defined on same lats as T, excluding poles
!  u defined between T lats, but on T lons
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer   :: i,j
      real(r8)  :: slats(jn1,2)  ! sines of latitude for (1)u and (2)T
      real(r8)  :: pi, pi180
      real(r8)  :: rlat, rlat_half
      real(r8)  :: tlat, ulat
      
      pi=4.d0*datan(1.d0)
      pi180=180.d0/pi
      rlat=pi/dble(jn1-1)
      rlat_half=0.5d0*rlat

      tlat=-0.5d0*pi   ! latitude of south pole
      glats(1,1)=pi180*rlat  ! a special value since otherwise not used
      glats(1,2)=pi180*tlat  
      slats(1,1)=0.d0  ! a value not used
      slats(1,2)=-1.d0
      do j=2,jn1
        ulat=tlat+rlat_half
        tlat=tlat+rlat
        glats(j,1)=pi180*ulat
        glats(j,2)=pi180*tlat
        slats(j,1)=dsin(ulat)
        slats(j,2)=dsin(tlat)
      enddo
!
       jweights(1,1)=0.d0  ! not used
       jweights(1,2)=0.5d0*(1.d0+slats(2,1))
       do j=2,jn1-1
         jweights(j,1)=0.5d0*(slats(j,2)-slats(j-1,2))
         jweights(j,2)=0.5d0*(slats(j+1,1)-slats(j,1))
       enddo
       jweights(jn1,1)=0.5d0*(slats(jn1,2)-slats(jn1-1,2))
       jweights(jn1,2)=0.5d0*(1.d0-slats(jn1,1))
    
       end subroutine horiz_grid_
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pkappa --- Compute (p/po)**kappa from delp 
!
!
! !INTERFACE:
!
      subroutine pkappa_ (im1,jn1,nl1,ptop,delp_r,pk_r,delp_x,pk_x,compute)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)      :: im1,jn1,nl1
      real(r8), intent(in)     :: ptop
      real(r8), intent(in)     :: delp_r(im1,jn1,nl1)
      character(len=*), optional, intent(in) :: compute

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout)  :: pk_r(im1,jn1,nl1)
      real(r8), optional, intent(inout)  :: delp_x(im1,jn1,nl1)
      real(r8), optional, intent(inout)  :: pk_x(im1,jn1,nl1)

! !DESCRIPTION:
!
!  Compute (p/po)**kappa from delp for reference field.
!  If requested also compute the corresponding TLM field from pert delp.
!  Or, if requested update the adjoint delp field by an input adjoint
!  field of (p/po)**kappa. Note that pk\_x and delp\_x are not referenced if
!  only the nonlinear calculation is performed.  The \_x refers to either the
!  \_tl or \_ad fields.
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  27Oct2007  Todling    Turned perturbation fields to optionals
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k,ier               
      real(r8)            :: f1
      real(r8)            :: f2
      real(r8)            :: pe_r(nl1+1)
      real(r8)            :: pke_r(nl1+1)
      real(r8)            :: peln_r(nl1+1)
      real(r8)            :: pe_x(nl1+1)
      real(r8)            :: pke_x(nl1+1)
      real(r8)            :: peln_x(nl1+1)      

! sanity check
      if ( present(compute) ) then
           ier=0
           if(.not.present(delp_x)) ier=1
           if(.not.present(pk_x)  ) ier=1
           if (ier/=0) then
               print*, 'inconsistent call to pkappa_'
               call exit(ier)
           endif
      endif

      do j=1,jn1
        do i=1,im1
!
! compute reference values
          pe_r(1)=ptop
          pke_r(1)=pe_r(1)**kappa
          peln_r(1)=log(pe_r(1))
          do k=2,nl1+1
            pe_r(k)=pe_r(k-1)+delp_r(i,j,k-1)
            pke_r(k)=pe_r(k)**kappa
            peln_r(k)=log(pe_r(k))
          enddo
          do k=1,nl1
            pk_r(i,j,k) = (pke_r(k+1)-pke_r(k)) / &
                          (kappa*(peln_r(k+1)-peln_r(k)))
          enddo

!
! compute tlm values
          if ( present(compute) ) then

          if (compute == 'TLM') then
            pe_x(1)=0.d0
            pke_x(1)=0.d0
            peln_x(1)=0.d0
            do k=2,nl1+1
              pe_x(k)=pe_x(k-1)+delp_x(i,j,k-1)
              pke_x(k)=pe_x(k)*kappa*pke_r(k)/pe_r(k)
              peln_x(k)=pe_x(k)/pe_r(k)
            enddo
            do k=1,nl1
              pk_x(i,j,k)=(pke_x(k+1)-pke_x(k)) / &
                          (kappa*(peln_r(k+1)-peln_r(k))) &
                         -(kappa*(peln_x(k+1)-peln_x(k)))*(pk_r(i,j,k)**2) / &
                                (pke_r(k+1)-pke_r(k))
            enddo
          endif
!
! compute adm values
          if (compute == 'ADM') then
            pke_x(:)=0.d0
            peln_x(:)=0.d0
            pe_x(:)=0.d0
            do k=1,nl1
              f1=pk_x(i,j,k)/(kappa*(peln_r(k+1)-peln_r(k)))
              f2=kappa*pk_x(i,j,k)*(pk_r(i,j,k)**2) / & 
                 (pke_r(k+1)-pke_r(k)) 
              pke_x(k+1)=pke_x(k+1)+f1
              pke_x(k  )=pke_x(k  )-f1
              peln_x(k+1)=peln_x(k+1)-f2
              peln_x(k  )=peln_x(k  )+f2
            enddo
            do k=2,nl1+1
              pe_x(k)=pke_x(k)*kappa*pke_r(k)/pe_r(k) +peln_x(k)/pe_r(k)
            enddo
            do k=nl1+1,2,-1 
              pe_x(k-1)=pe_x(k-1)+pe_x(k)
              delp_x(i,j,k-1)=delp_x(i,j,k-1)+pe_x(k)
            enddo
          endif

          endif ! <compute>
!
        enddo
      enddo

      end subroutine pkappa_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pt2vpt ---Replace potential temperature by virtual potential temperature.
!
!
! !INTERFACE:
!
      subroutine  pt2vpt_ (im1,jn1,nl1,pt,q)

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)  :: im1,jn1,nl1
      real(r8), intent(in)    :: q(im1,jn1,nl1)  ! specific humidity

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: pt(im1,jn1,nl1) ! see above

! !DESCRIPTION:
!
! Replace potential temperature by virtual potential temperature.
! This will also work for replacing temperature by virtual temperature.
! This will also work if pt is a _ad field of either virtual pot. temp. or
! virt. temp.as long as the input q is the reference q field and no q_ad
! needs to be computed.
!
!  pt explanation:
!      if input pt = pot. temp.   then output pt = virt. pot. temp.
!      if input pt =      temp.   then output pt = virt.      temp.
!      if input pt = dJ/d(virt. pot. temp.) then output pt = dJ/d(pot. temp.)
!      if input pt = dJ/d(virt.      temp.) then output pt = dJ/d(     temp.)
!  These latter are based on the need for adjoint calculations to be based on
!      a calculation of virt T being computed from T, the latter being considered
!      as the initial field, the former as a derived field. They assume that there
!      is no prior dependence on virt. pot. temp. (or virt. temp.) and that no
!      dJ/dq needs to be computed.

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k               

      do k=1,nl1
        do j=1,jn1
          do i=1,im1 
            pt(i,j,k)=pt(i,j,k)*(1.d0+zvir*q(i,j,k))
          enddo
        enddo
      enddo

      end subroutine pt2vpt_

!
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ptv2pt ---Replace virtual potential temperature by potential temperature
!
!
! !INTERFACE:
!
      subroutine ptv2pt_ (imr,jnp,nl,pt,q)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)  :: imr,jnp,nl
      real(r8), intent(in)    :: q(imr,jnp,nl)   ! specific humidity

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: pt(imr,jnp,nl)  ! ptv on input; pt on output

! !DESCRIPTION:
!
! Replace virtual potential temperature by potential temperature.
! This will also work for replacing virtual temperature by temperature.
! This will also work if pt is a _tl field as long as q_tl=0 and the input q
! is the reference q field.
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k               

      do k=1,nl
        do j=1,jnp
          do i=1,imr 
            pt(i,j,k)=pt(i,j,k)/(1.d0+zvir*q(i,j,k))
          enddo
        enddo
      enddo

      end subroutine ptv2pt_
!
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pt2t --- Replace potential temperature by temperature.
!
!
! !INTERFACE:
!
      subroutine pt2t_ (im1,jn1,nl1,pt,pk)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)  :: im1,jn1,nl1
      real(r8), intent(in)    :: pk(im1,jn1,nl1) ! (p/po)**kappa 

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: pt(im1,jn1,nl1) ! pt on input, t on output

! !DESCRIPTION:
!
! Replace potential temperature by temperature.
! This will also work if pt is a _tl field if pk_tl=0.

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k               

       do k=1,nl1
        do j=1,jn1
          do i=1,im1 
            pt(i,j,k)=pt(i,j,k)*pk(i,j,k)
          enddo
        enddo
      enddo

      end subroutine pt2t_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pt2t_tl --- Replace perturbation potential temperature by 
!                       perturbation temperature      
!
!
! !INTERFACE:
!
subroutine  pt2t_tl_ (im1,jn1,nl1,pt,pt_tl,pk,pk_tl)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)  :: im1,jn1,nl1
      real(r8), intent(in)    :: pt(im1,jn1,nl1)     ! reference theta
      real(r8), intent(in)    :: pk(im1,jn1,nl1)     ! reference (p/po)**kappa
      real(r8), intent(in)    :: pk_tl(im1,jn1,nl1)  ! perturbed pk

! !OUTPUT PARAMETERS:

      real(r8), intent(out)   :: pt_tl(im1,jn1,nl1)  ! theta on input, t on output

! !DESCRIPTION:
!
! Replace perturbation potential temperature by perturbation temperature,
! accounting for possible perturbations of pk
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k               

      do k=1,nl1
        do j=1,jn1
          do i=1,im1 
            pt_tl(i,j,k)=pt(i,j,k)*pk_tl(i,j,k)+pt_tl(i,j,k)*pk(i,j,k)
          enddo
        enddo
      enddo

      end subroutine pt2t_tl_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pt2t_ad --- Adjoint of routine pt2t_ad
!
!
! !INTERFACE:
!
      subroutine pt2t_ad_ (im1,jn1,nl1,pt,pt_ad,pk,pk_ad)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)    :: im1,jn1,nl1
      real(r8), intent(in)    :: pt(im1,jn1,nl1)     ! reference pot. temp.
                                                     ! adjoint theta on output
      real(r8), intent(in)    :: pk(im1,jn1,nl1)     ! reference (p/po)**kappa

! !OUTPUT PARAMETERS:

      real(r8), intent(out)   :: pk_ad(im1,jn1,nl1)  ! adjoint pk

! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: pt_ad(im1,jn1,nl1)  ! adjoint T on input

! !DESCRIPTION:
!
! Adjoint of routine pt2t\_tl.  This replaces dJ/dT by dJ/d(pt) and
! creates a resulting sensitivity dJ/d(pk).  Note that it is assumed
! that there is no pre-existing value of dJ/d(pk) that needs to be added to.
! This routine is necessary when J is defined by T, but the model is begun
! by a field that depends on pt.
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k               

      do k=1,nl1
        do j=1,jn1
          do i=1,im1  
            pk_ad(i,j,k)=pt_ad(i,j,k)*pt(i,j,k)
          enddo
          do i=1,im1 
            pt_ad(i,j,k)=pt_ad(i,j,k)*pk(i,j,k)
          enddo
        enddo
      enddo

      end subroutine pt2t_ad_
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: read_data --- Read in GEOS4 dynamics vector an input into 
!                          format for init\_adj code
!
!
! !INTERFACE:
!
      subroutine read_data_ (im,jm,km,nfields,fnames,nymd,nhms,job, &
                             fields3d,what,rc)

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: im, jm, km, nfields
      character(len=*), intent(in) :: fnames(nfields)
      character(len=*), intent(in) :: job
      character(len=*), intent(in) :: what

! !INPUT/OUTPUT PARAMETERS:

      integer, intent(inout) :: nymd, nhms

! !OUTPUT PARAMETERS:
 
      real(r8), intent(out) :: fields3d(:,:,:,:,:)
      integer,  intent(out) :: rc

! !DESCRIPTION:
!  
! !REMARKS:
!     The support to GEOS-5 like vectors provided here is only 
!  a quick fix. The way currently implemented the code will convert
!  the GEOS-5 input vector to GEOS-4-like and all will proceed as
!  the original code. There is quite a lot of redundancy doing things
!  this way and ultimately this code should be upgraded to handle 
!  GEOS-5 vectors properly. (RT)
!
! !SEE ALSO: write_data
!    
! !REVISION HISTORY:
!
!  04Oct2004  Winslow    Initial algorithm
!  27Oct2007  Todling    Support to GEOS-5 input vector
!  08Jan2008  Todling    Add opt to pick date
!  14Sep2013  Todling    No longer calls g5tog4 since norms handle A-grid winds
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: ier, indexF
      logical :: pick

      type(dyn_vect) w

      rc = 0
      pick = .false.
      if ( nymd/=0 .and. nhms/=0 ) pick=.true.

      call dyn_null ( w )
      if ( vectype_==5 ) then
         ! Note: in case the input vector is GEOS-5 compliant we will read it as such, but will
         !       need to flip it to make it GEOS-4 orientation-compliant ... before converting the
         !       whole vector to GEOS-4 - ultimately, we should make proper changes in this code
         !       to fully support GEOS-5 vectors and not have to convert the input to GEOS-4 type.
         if ( pick ) then
              call dyn_get ( trim(job), nymd, nhms, w, ier, vectype=vectype_, forceflip=.true., timidx=0 )
         else
              call dyn_get ( trim(job), nymd, nhms, w, ier, vectype=vectype_, forceflip=.true. )
         endif

         ! Now convert vector to GEOS-4
         !_RT no longer need this since norms can now handle A-grid winds
         ! call g5tog4_ ( w, what )
      else
         if ( pick ) then
              call dyn_get ( trim(job), nymd, nhms, w, ier, timidx=0 )
         else
              call dyn_get ( trim(job), nymd, nhms, w, ier )
         endif
      endif
      rc = ier
      if ( rc/=0 ) then
        call dyn_clean ( w )
        return
      endif
      if (w%grid%im.ne.im .or. w%grid%jm.ne.jm .or. w%grid%km.ne.km ) then
        print*,'***ERROR READING DYNAMIC VECTOR***'
        print*,'read_data: dimensions do not agree'
        print*,w%grid%im,w%grid%jm,w%grid%km
        print*,im,jm,km
        print*,'***EXIT ON READING ERROR***'
        rc = 2
        call dyn_clean ( w )
        return
      endif
      if (ier.ne.0) then
        print*,'failed to read data from ',job
      else
        print*,'read data from ',job
      endif
                                                                       
      call find_name_ (nfields,fnames,'u__',indexF)
      if (indexF .gt. 0) fields3d(:,:,:,indexF,1)=reshape(w%u,(/w%grid%im,w%grid%jm,w%grid%km/))
 
      call find_name_ (nfields,fnames,'v__',indexF)
      if (indexF .gt. 0) fields3d(:,:,:,indexF,1)=reshape(w%v,(/w%grid%im,w%grid%jm,w%grid%km/))
 
      call find_name_ (nfields,fnames,tname_,indexF)
      if (indexF .gt. 0) fields3d(:,:,:,indexF,1)=reshape(w%pt,(/w%grid%im,w%grid%jm,w%grid%km/))
 
      call find_name_ (nfields,fnames,'q__',indexF)
      if (indexF .gt. 0) fields3d(:,:,:,indexF,1)=reshape(w%q(:,:,:,1),(/w%grid%im,w%grid%jm,w%grid%km/))
 
      call find_name_ (nfields,fnames,'dp_',indexF)
      if (indexF .gt. 0) fields3d(:,:,:,indexF,1)=reshape(w%delp,(/w%grid%im,w%grid%jm,w%grid%km/))

      call find_name_ (nfields,fnames,'ps_',indexF)
      if (indexF .gt. 0) fields3d(:,:,1,indexF,1)=reshape(w%ps,(/w%grid%im,w%grid%jm/))

      call dyn_clean ( w )
 
      end subroutine read_data_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: transform_T0_ --- 
!
!
! !INTERFACE:
!
      subroutine transform_T0_ (im1,jn1,nl1,ptop,nfields,fnames,  &
                                fields_tl,fields_rf)

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1,nfields
      real(r8), intent(in) :: ptop

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
      character(len=*), intent(inout) :: fnames(nfields) 
      real(r8), intent(inout) :: fields_tl(im1,jn1,nl1,nfields)    
      real(r8), intent(inout) :: fields_rf(im1,jn1,nl1,nfields)    

! !DESCRIPTION:
!
! The last index of fields refers to (1) _tl fields or (2) ref fields.
! The input reference 'T__' field becomes 'vpT' fields on output 
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  07Jan2008  Todling    Handle either vT or vpT
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      character(len=*), parameter :: myname_ = myname//'*transform_T0_'
      real(r8)  :: pk(im1,jn1,nl1), dum1(im1,jn1,nl1)
      real(r8)  :: pk_tl(im1,jn1,nl1)
      integer :: indexT, indexPT, indexQ, indexdp
      integer :: i,j,k,nf

      call find_name_ (nfields,fnames,'T__',indexT)
      call find_name_ (nfields,fnames,'q__',indexQ)
      call find_name_ (nfields,fnames,'dp_',indexdP)
!
      if (indexT.gt. 0) then                         
  
           if ( tname_ == 'vpT' ) then

!  This first call pkappa simply computes pk and pk_tl using option 'TLM'
!
              call pkappa_ (im1,jn1,nl1,ptop,fields_rf(:,:,:,indexdP),pk,   &
                            delp_x=fields_tl(:,:,:,indexdP),pk_x=pk_tl,compute='TLM')

!  Replace reference T by reference pT first
!
              call t2pt_ (im1,jn1,nl1,fields_rf(:,:,:,indexT),pk)
!
!  Replace T_tl by pT_tl then
!
              call t2pt_tl_ (im1,jn1,nl1,fields_rf(:,:,:,indexT),  &
                             fields_tl(:,:,:,indexT),pk,pk_tl)
          endif
!  Also for the pert and reference state, replace pT by vpT for output. 
!  Identical fields_rf(:,:,:,indexQ) fields are used here because the assumption is 
!  that there is no q perturbation.
!
          call pt2vpt_ (im1,jn1,nl1,                                       &
                       fields_tl(:,:,:,indexT),fields_rf(:,:,:,indexQ))
          call pt2vpt_ (im1,jn1,nl1,                                       &
                       fields_rf(:,:,:,indexT),fields_rf(:,:,:,indexQ))
!
!  Finish conversion stating vpT in temperature fields name
!
          fnames(indexT) = tname_
          print *, 'Temperature field reset to ', fnames(indexT)

      else
	  call MP_die (myname_,'failed to convert temperature',99)
      endif
!
      end subroutine transform_T0_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: transform_T0_adinv_ --- The adjoint of routine transform_T_tl.
!
!
! !INTERFACE:
!
      subroutine transform_T0_adinv_ (im1,jn1,nl1,ptop,nfields,fnames,  &
                                      fields_ad,fields_rf)

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1,nfields
      real(r8), intent(in) :: ptop

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
      character(len=*), intent(inout) :: fnames(nfields) 
      real(r8), intent(inout) :: fields_ad(im1,jn1,nl1,nfields,1)    
      real(r8), intent(inout) :: fields_rf(im1,jn1,nl1,nfields,1)    

! !DESCRIPTION:
! The adjoint of routine= transform\_T\_tl.
! This sequence the operations required to initialize the adjoint virtual
! potential temperature field from the adjoint temperature field.  An adjoint delp
! field is also computed. NOTE: it is assumed that the input reference field
! is temperature.
!
! The last index of fields refers to (1) \_ad fields or (2) ref fields.
! The input reference vpT field is restored on output (within precision of
! inverse calculation)
!
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  27Oct2007  Todling    Updated interface to pkappa
!  07Jan2008  Todling    Handle either vT or vpT
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      character(len=*), parameter :: myname_ = myname//'*transform_T0_adinv_'
      real(r8)  :: pk(im1,jn1,nl1)
      real(r8)  :: pk_ad(im1,jn1,nl1)
      integer :: indexT, indexPT, indexVT, indexQ, indexdp
      integer :: i,j,k,nf
      real(r8) :: summ

      call find_name_ (nfields,fnames,'T__',indexT)
      call find_name_ (nfields,fnames,'q__',indexQ)
      call find_name_ (nfields,fnames,'dp_',indexdP)
      call find_name_ (nfields,fnames,'vT_',indexVT)

      if ( indexT.gt.0 .and. indexVT.gt.0 .and. indexT==indexVT ) then
	   call MP_die (myname_,'incorrect index setting',99)
      endif

      if (indexT>0 .and. indexVT==0 .and. tname_=='vpT') then                         
  
!        This first call pkappa simply computes pk
         call pkappa_ (im1,jn1,nl1,ptop,fields_rf(:,:,:,indexdP,1),pk)

!        Replace reference T by reference pT
         call t2pt_ (im1,jn1,nl1,fields_rf(:,:,:,indexT,1),pk)

!        transform adjoint T to adjoint pT field
         call pt2t_ad_ (im1,jn1,nl1,fields_rf(:,:,:,indexT,1),     &
                          fields_ad(:,:,:,indexT,1),pk,pk_ad)
         call pkappa_ (im1,jn1,nl1,ptop,fields_rf(:,:,:,indexdp,1),pk,   &
                       delp_x=fields_ad(:,:,:,indexdp,1),pk_x=pk_ad,compute='ADM')

!
!  Compute adjoint pTv field from adjoint pT field
!  Also for the reference state, replace pT by pTv as originally input. 
!  Identical routines are called here because the assumption is 
!  that there is no q perturbation.
          call pt2vpt_ (im1,jn1,nl1,                                       &
                        fields_ad(:,:,:,indexT,1),fields_rf(:,:,:,indexQ,1))

          call pt2vpt_ (im1,jn1,nl1,                                       &
                        fields_rf(:,:,:,indexT,1),fields_rf(:,:,:,indexQ,1))

          fnames(indexT) = tname_
          print *, 'Input ad-temperature field reset to ', fnames(indexT)
      endif

!  Compute adjoint Tv field from adjoint T field
      if (indexT>0 .and. indexVT==0 .and. tname_=='vT_') then                         
          call pt2vpt_ (im1,jn1,nl1,                                       &
                        fields_ad(:,:,:,indexT,1),fields_rf(:,:,:,indexQ,1))

          call pt2vpt_ (im1,jn1,nl1,                                       &
                        fields_rf(:,:,:,indexT,1),fields_rf(:,:,:,indexQ,1))
          fnames(indexT) = tname_
          print *, 'Input ad-temperature field reset to ', fnames(indexT)
      endif

      print *, 'Complete ', myname_
      end subroutine transform_T0_adinv_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: transform_T2 --- virtual (potential) temperature to dry temperature
!
!
! !INTERFACE:
!
      subroutine transform_T2_ (im1,jn1,nl1,ptop,nfields,fnames,fields)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: im1,jn1,nl1,nfields
      real(r8), intent(in) :: ptop

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
      character(len=*), intent(inout) :: fnames(nfields) 
      real(r8), intent(inout) :: fields(im1,jn1,nl1,nfields)    

! !DESCRIPTION: Convert virtual (potential) temperature into dry temperature
!    
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  27Oct2007  Todling    Updated interface to pkappa
!  07Jan2008  Todling    Handle either vT or vpT
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      real(r8) :: pk(im1,jn1,nl1)
      integer :: indexvpT, indexQ, indexP
 
      call find_name_ (nfields,fnames,tname_,indexvpT)  ! vpT or vT 
      call find_name_ (nfields,fnames,'q__',indexQ)
!
! compute T from vPT 
!
      if (indexvpT .gt. 0 ) then
          call ptv2pt_ (im1,jn1,nl1,fields(:,:,:,indexvpT),fields(:,:,:,indexQ))
          if ( tname_ == 'vpT' ) then
              call find_name_ (nfields,fnames,'dp_',indexP)
              call pkappa_ (im1,jn1,nl1,ptop,fields(:,:,:,indexP),pk)
              call pt2t_ (im1,jn1,nl1,fields(:,:,:,indexvpT),pk)
          endif
          fnames(indexvpT) = 'T__'                    ! now converted to T
          print *, 'Temperature field reset to ', 'T__'
      endif

      end subroutine transform_T2_

 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: transform_T1_tl --- Transform _tl field of virtual potential 
!                              temperature to either potential temperature
!                              or temperature
! !INTERFACE:
!
      subroutine transform_T1_tl_ (imr,jnp,nl,ptop,nfields,fnames,fields)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: imr,jnp,nl,nfields
      real(r8), intent(in) :: ptop
      character(len=*), intent(inout) :: fnames(nfields) 

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: fields(imr,jnp,nl,nfields,2)    

! !DESCRIPTION:
!
!  Transform _tl field of virtual potential temperature to either
!  potential temperature or temperature, depending on specification of fnames
!
!  The last index of fields refers to (1) _tl fields or (2) ref fields.
!  Note that both _tl and reference fields are transformed by this routine.
!
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  07Jan2008  Todling    Generalized to handle either vT and vpT
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      real(r8), allocatable  :: pk(:,:,:)
      real(r8), allocatable  :: pk_tl(:,:,:)
      integer :: indexT, indexPT, indexVT, indexVPT, indexQ, indexP

      call find_name_ (nfields,fnames,'vT_',indexVT)
      call find_name_ (nfields,fnames,'vpT',indexVPT)
      call find_name_ (nfields,fnames,'PT_',indexPT)
      call find_name_ (nfields,fnames,'q__',indexQ)

!
! compute PT from VPT if requested or T from VT
!
      if ( indexPT.gt.0 ) then    
        indexT = indexPT
        call ptv2pt_ (imr,jnp,nl,                                         &
                      fields(:,:,:,indexPT,2),fields(:,:,:,indexQ,2))
        call ptv2pt_ (imr,jnp,nl,                                         &
                      fields(:,:,:,indexPT,1),fields(:,:,:,indexQ,2))
        fnames(indexT) = 'PT_'
        print *, 'Virtual potential temperature field reset to ', 'PT_'
      endif
!
! Or calculate T from VPT
! The call to pkappa here computes both ref. and _tl values of (p/po)**kappa
!
      if ( indexVPT.gt.0 ) then                         
        indexT = indexVPT
        call ptv2pt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexT,2),fields(:,:,:,indexQ,2))
        call ptv2pt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexT,1),fields(:,:,:,indexQ,2))

        call find_name_ (nfields,fnames,'dp_',indexP)
        allocate (pk(imr,jnp,nl))
        allocate (pk_tl(imr,jnp,nl))
        call pkappa_ (imr,jnp,nl,ptop,fields(:,:,:,indexP,2),pk,   &
                      delp_x=fields(:,:,:,indexP,1),pk_x=pk_tl,compute='TLM')

        pk_tl = 0.d0   ! (EN) SVEC calculation in SVEC run does not count for pk_tl
                       ! (EN) To match the norm calculation in SVEC pk_tl is set to ZERO here 

        call pt2t_tl_ (imr,jnp,nl,fields(:,:,:,indexT,2),     &
                       fields(:,:,:,indexT,1),pk,pk_tl)
    
        call pt2t_ (imr,jnp,nl,fields(:,:,:,indexT,2),pk)
        deallocate (pk,pk_tl)
        fnames(indexT) = 'T__'
        print *, 'Virtual potential temperature field reset to ', 'T__'
      endif

! Or replace vT with T

      if ( indexVT.gt.0 ) then                         
        indexT = indexVT
        call ptv2pt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexVT,2),fields(:,:,:,indexQ,2))
        call ptv2pt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexVT,1),fields(:,:,:,indexQ,2))
        fnames(indexT) = 'T__'
        print *, 'Virtual temperature field reset to ', 'T__'
      endif
!
      end subroutine transform_T1_tl_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: transform_T0_ad_ --- Transform _ad field of virtual potential 
!                              temperature to either potential temperature
!                              or temperature.
!
!
! !INTERFACE:
!
      subroutine transform_T0_ad_ (imr,jnp,nl,ptop,nfields,fnames,   &
                                   fields_ad, fields_ref)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: imr,jnp,nl,nfields
      real(r8), intent(in) :: ptop
      character(len=*), intent(inout) :: fnames(nfields) 

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: fields_ad (imr,jnp,nl,nfields,1)
      real(r8), intent(inout) :: fields_ref(imr,jnp,nl,nfields,1)

! !DESCRIPTION:
!
!  Transform _ad field of virtual potential temperature to either
!  potential temperature or temperature, depending on specification of fnames
!
!  The last index of fields refers to (1) _ad fields or (2) ref fields.
!  Note that both _ad and reference fields are transformed
!  DELP is also transformed if pt_ad is replaced by t_ad: in the first case
!  delp is dJ/d(delp) with theta held fixed; in the latter case, with T held
!  fixed (the d/d here are actually partial derivitives).
!
! !SEE ALSO:
!
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  09Jan2008  Todling    - Generalized to handle either vT and vpT
!                        - Revisited calculation of vpT to T
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      character(len=*), parameter :: myname_ = myname//'*transform_T0_ad_'
      real(r8), allocatable  :: pk(:,:,:)
      real(r8), allocatable  :: pk_ad(:,:,:)

      integer :: indexT, indexQ, indexP
                                                                                                                   
      call find_name_ (nfields,fnames,'q__',indexQ)
      call find_name_ (nfields,fnames,tname_,indexT)
      call find_name_ (nfields,fnames,'dp_',indexP)

! If so, compute T from vpT

      if (tname_=='vpT') then                         
        call ptv2pt_ (imr,jnp,nl,                                         &
                     fields_ref(:,:,:,indexT,1),fields_ref(:,:,:,indexQ,1))
        call pt2vpt_ (imr,jnp,nl,                                         &
                      fields_ad(:,:,:,indexT,1),fields_ref(:,:,:,indexQ,1))

        allocate (pk(imr,jnp,nl),pk_ad(imr,jnp,nl))
        pk_ad = 0.d0
        call pkappa_  (imr,jnp,nl,ptop,fields_ref(:,:,:,indexP,1),pk)
        call t2pt_ad_ (imr,jnp,nl,fields_ref(:,:,:,indexT,1), &
                                  fields_ad (:,:,:,indexT,1),pk,pk_ad )
        call pt2t_    (imr,jnp,nl,fields_ref(:,:,:,indexT,1),pk)
        deallocate (pk,pk_ad)

        fnames(indexT) = 'T__'
        print *, 'Input ad-virtual potential temperature field reset to ', fnames(indexT)
      endif

! If so, compute T from vT

      if (tname_=='vT_') then                         
        call ptv2pt_ (imr,jnp,nl,                                         &
                     fields_ref(:,:,:,indexT,1),fields_ref(:,:,:,indexQ,1))
        call pt2vpt_ (imr,jnp,nl,                                         &
                     fields_ad (:,:,:,indexT,1),fields_ref(:,:,:,indexQ,1))
        fnames(indexT) = 'T__'
        print *, 'Input ad-virtual temperature field reset to ', fnames(indexT)
       endif
 
      end subroutine transform_T0_ad_
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: transform_T1_ad --- Transform _ad field of virtual potential 
!                              temperature to either potential temperature
!                              or temperature.
!
!
! !INTERFACE:
!
      subroutine transform_T1_ad_ (imr,jnp,nl,ptop,nfields,fnames,fields)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: imr,jnp,nl,nfields
      real(r8), intent(in) :: ptop
      character(len=*), intent(in) :: fnames(nfields) 

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: fields(imr,jnp,nl,nfields,2)    

! !DESCRIPTION:
!
!  Transform _ad field of virtual potential temperature to either
!  potential temperature or temperature, depending on specification of fnames
!
!  The last index of fields refers to (1) _ad fields or (2) ref fields.
!  Note that both _ad and reference fields are transformed
!  DELP is also transformed if pt_ad is replaced by t_ad: in the first case
!  delp is dJ/d(delp) with theta held fixed; in the latter case, with T held
!  fixed (the d/d here are actually partial derivitives).
!
! !SEE ALSO:
!
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  09Jan2008  Todling    Knob to calculate T from vT
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      character(len=*), parameter :: myname_ = myname//'*transform_T1_ad_'
      real(r8), allocatable  :: pk(:,:,:)
      real(r8), allocatable  :: pk_ad(:,:,:)

      integer :: indexT, indexPT, indexvpT, indexvT, indexQ, indexP
                                                                                                                   
      call find_name_ (nfields,fnames,'T__',indexT)
      call find_name_ (nfields,fnames,'vpT',indexvpT)
      call find_name_ (nfields,fnames,'vT_',indexvT)
      call find_name_ (nfields,fnames,'PT_',indexPT)
      call find_name_ (nfields,fnames,'q__',indexQ)

! Compute PT_ad from VPT_ad if requested (using adjoint of VPT=PT*F(Q))
! Also replace reference VPT by PT 
!
      if (indexPT.gt.0) then    
        call ptv2pt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexPT,2),fields(:,:,:,indexQ,2))
        call pt2vpt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexPT,1),fields(:,:,:,indexQ,2))
      endif
!
! also compute T from PT if requested (using adjoint of PT=T/PK)
!
      if (indexT.gt.0) then                         
        call ptv2pt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexT,2),fields(:,:,:,indexQ,2))
        call pt2vpt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexT,1),fields(:,:,:,indexQ,2))

        call find_name_ (nfields,fnames,'dp_',indexP)
        allocate (pk(imr,jnp,nl))
        allocate (pk_ad(imr,jnp,nl))
! first need reference (p/po)**kappa (f(...,1) and pk_ad dummy here)
        call pkappa_ (imr,jnp,nl,ptop,fields(:,:,:,indexP,2),pk,   &
                      delp_x=fields(:,:,:,indexP,1),pk_x=pk_ad,compute='NLM')
        call t2pt_ad_ (imr,jnp,nl,fields(:,:,:,indexT,2),      &
                       fields(:,:,:,indexT,1),pk,pk_ad )
        call pkappa_ (imr,jnp,nl,ptop,fields(:,:,:,indexP,2),pk,   &
                      delp_x=fields(:,:,:,indexP,1),pk_x=pk_ad,compute='ADM')
        call pt2t_ (imr,jnp,nl,fields(:,:,:,indexT,2),pk)
        deallocate (pk,pk_ad)
      endif
!
! If so, replace vT with T
!
      if (indexVT.gt.0) then                         
        call ptv2pt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexVT,2),fields(:,:,:,indexQ,2))
        call pt2vpt_ (imr,jnp,nl,                                         &
                     fields(:,:,:,indexVT,1),fields(:,:,:,indexQ,2))
      endif
!
      print *, 'Complete ', myname_
      end subroutine transform_T1_ad_
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: t2pt --- Routine that replaces temperature by potential temperature
!
!
! !INTERFACE:
!
      subroutine  t2pt_ (im1,jn1,nl1,pt,pk)

!USES:

      implicit none

! !INPUT PARAMETERS:
 
      integer,  intent(in)    :: im1,jn1,nl1
      real(r8), intent(in)    :: pk(im1,jn1,nl1)     ! (p/po)**kappa

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: pt(im1,jn1,nl1)     ! potential temp.

! !DESCRIPTION:
  
! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k               

      do k=1,nl1
        do j=1,jn1
          do i=1,im1 
            pt(i,j,k)= pt(i,j,k)/pk(i,j,k)
          enddo
        enddo
      enddo

      end subroutine t2pt_
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: t2pt_ad --- Adjoint of a routine that replaces temperature by 
!                       potential temperature
!
!
! !INTERFACE:
!
      subroutine  t2pt_ad_ (imr,jnp,nl,pt,pt_ad,pk,pk_ad)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)  :: imr,jnp,nl
      real(r8), intent(in)    :: pt(imr,jnp,nl)     ! reference potential temp.
      real(r8), intent(in)    :: pk(imr,jnp,nl)     ! reference (p/po)**kappa

! !OUTPUT PARAMETERS:

      real(r8), intent(out)   :: pk_ad(imr,jnp,nl)  ! adjoint of pk field

! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: pt_ad(imr,jnp,nl)  ! adj. pt in, adj. T out

! !DESCRIPTION:
!
!  Adjoint of a routine that replaces temperature by potential temperature 
!  This assumes that there is no prior dependence on T or pk. 
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k               

      do k=1,nl
        do j=1,jnp
          do i=1,imr 
            pk_ad(i,j,k)=-pt_ad(i,j,k)*pt(i,j,k)/pk(i,j,k)
            pt_ad(i,j,k)= pt_ad(i,j,k)/pk(i,j,k)
          enddo
        enddo
      enddo

      end subroutine t2pt_ad_
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: t2pt_tl_ --- replaces pert temperature by potential temperature
!
!
! !INTERFACE:
!
      subroutine  t2pt_tl_ (imr,jnp,nl,pt,pt_tl,pk,pk_tl)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)  :: imr,jnp,nl
      real(r8), intent(in)    :: pt(imr,jnp,nl)     ! reference potential temp.
      real(r8), intent(in)    :: pk(imr,jnp,nl)     ! reference (p/po)**kappa
      real(r8), intent(in)   :: pk_tl(imr,jnp,nl)  ! pert of pk field

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: pt_tl(imr,jnp,nl)  ! pert: T in, theta out

! !DESCRIPTION:
!
!  Routine that replaces temperature by potential temperature 
!  Used only for testing corresponding adjoint
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k               

      do k=1,nl
        do j=1,jnp
          do i=1,imr 
            pt_tl(i,j,k)= pt_tl(i,j,k)/pk(i,j,k)    &
               -pk_tl(i,j,k)*pt(i,j,k)/pk(i,j,k)
          enddo
        enddo
      enddo

      end subroutine t2pt_tl_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: vert_grid --- compute vertical weights and set eta
!
!
! !INTERFACE:
!
      subroutine vert_grid_ (nl1, ks, ak, bk, kweights, ptop, &
                            pmid, delbk)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: nl1

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: ptop
      real(r8), intent(out) :: kweights(nl1)
      real(r8), intent(out) :: pmid(nl1)
      real(r8), intent(out) :: delbk(nl1)
      integer,intent(out)    :: ks
      real(r8),intent(out)   :: ak(nl1+1)
      real(r8),intent(out)   :: bk(nl1+1)

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
  
! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  20Jul2005  Elena N.   Linked constants in vert_grid (pref) to common constants
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer    :: k
      real(r8)   :: p(nl1+1)
      real(r8)   :: pint
      real(r8)   :: ps
     
      call set_eta(nl1, ks, ptop, pint, ak, bk)

      ps=pref
      do k=1,nl1+1
        p(k)=ak(k)+bk(k)*ps
      enddo

      do k=1,nl1
        kweights(k)=(p(k+1)-p(k))/p(nl1+1)
        pmid(k)=0.5*(p(k)+p(k+1))
        delbk(k)=bk(k+1) - bk(k) 
      enddo

      end subroutine vert_grid_  

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: write_data --- Write out perturbation/adjoint initialization
!
!
! !INTERFACE:
!
      subroutine write_data_ (im1,jn1,nl1,lm1,nfields2d,nfields3d,ptop,ks, &
                              ak,bk,ps,fields3d,fnames,ofname,nymd,nhms,what,&
                              idims)  

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: im1, jn1, nl1, lm1, nfields2d, nfields3d
      real(r8), intent(in)  :: ptop
      integer,  intent(in)  :: ks
      real(r8), intent(in)  :: ak(nl1+1)
      real(r8), intent(in)  :: bk(nl1+1)
      character(len=3), intent(in) :: fnames(nfields2d+nfields3d)
      character(len=*), intent(in) :: ofname  ! output file name
      integer , intent(in)  :: nhms, nymd
      character(len=*), intent(in) :: what  ! specify type of pert/field

      integer,          optional, intent(in) :: idims(2) ! alternative horiz. dim to interpolate to
                                                                 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: ps(im1,jn1)
      real(r8), intent(inout) :: fields3d(im1,jn1,nl1,nfields3d)
 
! !DESCRIPTION:
! 
! !SEE ALSO: read_data
!
! !REVISION HISTORY:
!
!  04Oct2004  Winslow    Initial algorithm
!  07Apr2005  Elena N.   Changed output precision to 32-bit
!  21Jul2005  Elena N.   Removed extraneous parameters
!  27Oct2007  Todling    Support to GEOS-5-type ouput vector
!  25Feb2008  Todling    Add option to interpolate perturbation on way out
!  05Mar2009  Todling    Update interface to dyn_init
!  27Oct2009  Todling    Pass vectype to dyn_init
!  19Jan2010  Todling    Flip fields on way out
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables
 
      character(len=*), parameter :: myname_ = 'write_data'

      integer :: ier
      integer, parameter :: frq = 3000
      logical, parameter :: fnew = .true.
      integer :: im2_, jn2_
      integer :: vectype
      integer :: fnhms, fnymd
      integer :: f_id                            ! field id
      integer :: ierr                            ! return code
      integer :: prec_das                        ! precision for DAS output
                                                 ! 0: 32 bits
                                                 ! 1: 64 bits (default)
      logical  forceflip_
      real(r8) :: aux2d(im1,jn1)
      real(r8) :: delp_ad(im1,jn1,nl1)
      real(r8) :: u_ad(im1,jn1,nl1)
      real(r8) :: v_ad(im1,jn1,nl1)
      real(r8) :: pt_ad(im1,jn1,nl1)
      real(r8) :: q_ad(im1,jn1,nl1,lm1)
      type (dyn_vect) :: w_f
      type (dyn_vect) :: w_i
!
! --------------- Initialize local arrays
!
      aux2d = 0.d0
!
      call getparamI_ ( 'vectype', vectype )

!     If so, convert resulting perturbation vector to GOES-5-type
!     -----------------------------------------------------------
      if ( vectype_==5 ) then
          if( what/='adm' .and. what/='tlm') then
              call MP_die (myname_,'this routine to be call for perturbations only',99)
          endif
          call g4tog5_pert_ ( im1,jn1,nl1,lm1,nfields2d,nfields3d,ptop,bk, &
                              ps,fields3d,fnames,what )
      endif

      call find_name_ (nfields3d,fnames,'dp_',f_id )
      delp_ad(:,:,:)=reshape(fields3d(:,:,:,f_id),(/im1,jn1,nl1/))
      call find_name_ (nfields3d,fnames,'u__',f_id )
      u_ad(:,:,:)=reshape(fields3d(:,:,:,f_id),(/im1,jn1,nl1/))
      call find_name_ (nfields3d,fnames,'v__',f_id )
      v_ad(:,:,:)=reshape(fields3d(:,:,:,f_id),(/im1,jn1,nl1/))
      call find_name_ (nfields3d,fnames,tname_,f_id )
      pt_ad(:,:,:)=reshape(fields3d(:,:,:,f_id),(/im1,jn1,nl1/))
      call find_name_ (nfields3d,fnames,'q__',f_id )
      q_ad(:,:,:,1)=fields3d(:,:,:,f_id)
      if (lm1>1) q_ad(:,:,:,2:lm1) = 0.d0  ! Assign zeros to other tracers for NOW 

      call dyn_null ( w_f )
      call Dyn_Init ( im1, jn1, nl1, lm1, ptop, ks, ak, bk, &
                      aux2d, aux2d, aux2d, aux2d, ps,       &
                      aux2d, aux2d, aux2d, aux2d, aux2d,    &
                      delp_ad, u_ad, v_ad, pt_ad, q_ad,     &
                      w_f, ierr, vectype=vectype )
      if ( ierr .ne. 0 ) then
         call MP_die (myname_,'Error from Dyn_Init',ierr)
      end if

      prec_das = 0
      fnhms = nhms
      fnymd = nymd

      im2_=im1; jn2_=jn1
      if ( present(idims) ) then
           if(idims(1)>0) im2_ = idims(1)
           if(idims(2)>0) jn2_ = idims(2)
      endif

!     Check whether or not to interpolate, and write out
!     --------------------------------------------------
      if ( im2_==im1 .and. jn2_==jn1  ) then

!        Write out perturbation in original grid
!        ---------------------------------------
         if ( vectype_==5 ) then
              forceflip_=.true.
         else
              forceflip_=.false.
         endif
         call dyn_put ( trim(ofname), fnymd, fnhms, prec_das, w_f, ierr, & 
                        new=fnew, freq=frq, vectype=vectype, forceflip=forceflip_ )
           if(ierr/=0) then
             call MP_die (myname_,'Error from dyn_put',ierr)
           else
             print*,'Wrote perturbation to ',trim(ofname)
           endif

      else

!          Initialize dimension of vertically interpolated vector
!          -------------------------------------------------------
           call dyn_null ( w_i )
           call dyn_init ( im2_, jn2_, w_f%grid%km, w_f%grid%lm, w_i, ierr, &
                           w_f%grid%ptop, w_f%grid%ks, w_f%grid%ak, w_f%grid%bk, vectype=vectype )
                   if (ierr/=0) then
                       call die ( myname, 'error initializing perturbation on way out ' )
                   endif

!          Interpolate perturbation now
!          ----------------------------
           call h_map_pert ( w_f, w_i, what, ierr )
                  if (ierr/=0) then
                      call die ( myname, 'error interpolating perturbation on way out ' )
                  else
                      print *, trim(myname), ': interpolated perturbation on way out'
                  endif

!          Write out interpolated perturbation
!          -----------------------------------
           call dyn_put ( trim(ofname), fnymd, fnhms, prec_das, w_i, ierr, & 
                          new=fnew, freq=frq, vectype=vectype )
             if(ierr/=0) then
               call MP_die (myname_,'Error from dyn_put',ierr)
             else
               print*,'Wrote interpolated perturbation to ',trim(ofname)
             endif

             call dyn_null ( w_i )

      endif

      call dyn_null ( w_f )
      end subroutine write_data_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: e_cross_product --- Cross product of two state vectors computed
!
!
! !INTERFACE:
!
      subroutine e_cross_product_ (imr,jnp,nl,nfields,ptop,        &
                   udomain, vdomain, ptdomain, kdomain, jweights,  &
                   fieldsA,fieldsB,psfieldA,psfieldB,              &
                   delp_ref,fnames, e) 

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: imr,jnp,nl,nfields
      real(r8), intent(in) :: ptop
      real(r8), intent(in) :: udomain(imr,jnp,nl)
      real(r8), intent(in) :: vdomain(imr,jnp,nl)
      real(r8), intent(in) :: ptdomain(imr,jnp,nl)
      integer, intent(in) :: kdomain(2)
      character(len=*), intent(in) :: fnames(nfields) 
      real(r8), intent(in)  :: jweights(jnp)
      real(r8), intent(in)  :: fieldsA(imr,jnp,nl,nfields)
      real(r8), intent(in)  :: fieldsB(imr,jnp,nl,nfields)
      real(r8), intent(in)  :: psfieldA(imr,jnp)   
      real(r8), intent(in)  :: psfieldB(imr,jnp)   
      real(r8), intent(in)  :: delp_ref(imr,jnp,nl)

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: e(nl+2,3)   

! !INPUT/OUTPUT PARAMETERS:
                                                                                                                     
! !DESCRIPTION:

!
! Cross product of two state vectors computed using the
! same total energy norm as defined for the FVGCM Singular
! Vector calculations. Separate contributions
! by kinetic and potential energies for each pressure layer are determined.
!

! !SEE ALSO:
!
                                                                                                                     
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  15Apr2005  Elena N.   Corrected loop indexing in V-energy calculations 
!                        Corrected sum(jweights)=1
!                        Corrected dsigma indexing  in V-energy calculations
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer  :: i,j,k
      integer  :: idu,idv,idt
      real(r8) :: dbleimr
      real(r8) :: tfac, ufac, pfac
      real(r8) :: dsigma(imr,jnp,nl)
      real(r8) :: ps_ref(imr,jnp)
      real(r8) :: jfac(jnp)      
      real(r8) :: faA(imr,jnp) , faB(imr,jnp) 
      real(r8) :: fb(imr,jnp) 
      real(r8) :: eu,ev,et,ep,eut,evt,ett

      dbleimr=dble(imr)

      do j=1,jnp
        jfac(j)=jweights(j)/dbleimr
      enddo

! One factor of 0.5 is for definition of energy as 1/2 * ...
! Assumption for division by area where sum(jweights)=1
      tfac=0.5d0*Cp/tref
      ufac=0.5d0
      pfac=0.5d0*R*tref/pref**2

      call find_name_ (nfields,fnames,'u__',idu)      
      call find_name_ (nfields,fnames,'v__',idv)      
      call find_name_ (nfields,fnames,'T__',idt)      
!
! compute 3-D varying dsigma
      call calculate_ps_ (imr,jnp,nl,ptop,delp_ref,ps_ref)
      do k=1,nl
        dsigma(:,:,k)=delp_ref(:,:,k)/ps_ref(:,:)
      enddo    
!
! Initialize e to zero
      e(:,:)=0.d0
!
!  Loop over levels
!
      eut = 0.0
      evt = 0.0
      ett = 0.0
      do k=kdomain(1),kdomain(2)
!     
!  Calculations for u field
!     
        faA = 0.0
        faB = 0.0
        do i=1,imr
          do j=1,jnp
            if (udomain(i,j,k) .ne. 0.0) then
              faA(i,j)=fieldsA(i,j,k,idu)
              faB(i,j)=fieldsB(i,j,k,idu)
            endif
          enddo
        enddo

        do i=1,imr
          fb(i,  1)=faA(i,  2)*faB(i,  2)
          fb(i,jnp)=faA(i,jnp)*faB(i,jnp)
          do j=2,jnp-1     
            fb(i,j)=0.5*(faA(i,j)*faB(i,j) + faA(i,j+1)*faB(i,j+1))
          enddo
        enddo

        eu=0.d0
        do j=1,jnp
          do i=1,imr
            eu=eu+dsigma(i,j,k)*jfac(j)*fb(i,j)
          enddo
        enddo
        eut = eut + eu*ufac
!     
!  Calculations for v field
!     
        faA = 0.0
        faB = 0.0
        do i=1,imr
          do j=1,jnp
            if (vdomain(i,j,k) .ne. 0.0) then
              faA(i,j)=fieldsA(i,j,k,idv)
              faB(i,j)=fieldsB(i,j,k,idv)
            endif
          enddo
        enddo

        do i=1,imr
          fb(i,  1)=faA(i,    2)*faB(i,    2)
          fb(i,jnp)=faA(i,jnp-1)*faB(i,jnp-1)
        enddo
        do j=2,jnp-1
          fb(imr,j)=0.5*(faA(1,j)*faB(1,j) + faA(imr,j)*faB(imr,j))
          do i=1,imr-1     
            fb(i,j)=0.5*(faA(i,j)*faB(i,j) + faA(i+1,j)*faB(i+1,j))
          enddo
        enddo

        ev=0.d0
        do j=1,jnp
          do i=1,imr
            ev=ev+dsigma(i,j,k)*jfac(j)*fb(i,j)
          enddo
        enddo
        evt = evt + ev*ufac

        e(k,1)=ufac*(eu+ev)
        e(nl+1,1)=e(nl+1,1)+e(k,1)
!     
!  Calculations for T field
!     
        faA = 0.0
        faB = 0.0
        do i=1,imr
          do j=1,jnp
            if (ptdomain(i,j,k) .ne. 0.0) then
              faA(i,j)=fieldsA(i,j,k,idt)
              faB(i,j)=fieldsB(i,j,k,idt)
            endif
          enddo
        enddo

        et=0.d0 
        do j=1,jnp
          do i=1,imr
            et=et+dsigma(i,j,k)*jfac(j)*faA(i,j)*faB(i,j)
          enddo
        enddo

        e(k,2)=tfac*et
        e(k,3)=e(k,1)+e(k,2)
        e(nl+1,2)=e(nl+1,2)+e(k,2)
!
!
      enddo ! loop over k
      ett = e(nl+1,2)
!
!
! Calculations for ps field
!
      faA = 0.0
      faB = 0.0
      k = kdomain(1)
      do i=1,imr
        do j=1,jnp
          if (ptdomain(i,j,k) .ne. 0.0) then
            faA(i,j)=psfieldA(i,j)
            faB(i,j)=psfieldB(i,j)
          endif
        enddo
      enddo
      ep=0.d0
      do j=1,jnp
        do i=1,imr
          ep=ep+jfac(j)*faA(i,j)*faB(i,j)
        enddo
      enddo

      e(nl+1,3)=e(nl+1,1)+e(nl+1,2)
      e(nl+2,1)=0.d0                            ! not used
      e(nl+2,2)=pfac*ep
      e(nl+2,3)=e(nl+2,2)+e(nl+1,3)             ! total E
      print*,'KEU, KEV, APET, APEP=',eut,evt,ett,e(nl+2,2)

      end subroutine e_cross_product_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: comp_kdomain --- compute the vertical indices which delineate 
!                            the vertical domain
!
!
! !INTERFACE:
!
      subroutine comp_kdomain_(im1, jn1, nl1, ptdomain, kdomain)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer , intent(in)  :: im1, jn1, nl1
      real(r8), intent(in) :: ptdomain(im1,jn1,nl1)

! !OUTPUT PARAMETERS:

      integer, intent(out) :: kdomain(2)

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
  
! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Oct2004  Winslow    Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: i,j,k
      logical  :: kflag

      kflag = .true.
      do j=1,jn1
        do i=1,im1
          do k=1,nl1
            if (ptdomain(i,j,k) .ne. 0.0) then
              kdomain(2) = k
              if ((kflag)) then
                kflag = .false.
                kdomain(1) = k
              endif
            endif
          enddo
        enddo
      enddo

      end subroutine comp_kdomain_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: g5tog4_ --- Converts G-5 vectors to G-4 (NEVER make this public)
!
!
! !INTERFACE:

      subroutine g5tog4_ ( w, what ) ! CAUTION: never ever export this routine

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: what

! !INPUT/OUTPUT PARAMETERS:

      type(dyn_vect) :: w
  
! !DESCRIPTION: Convert GEOS-5 full vector to GEOS-4 like
!
! !REVISION HISTORY:
!
!  18Nov2007  Todling    Initial (temporary) code
!  18Dec2008  Errico/RT  Bug fix: zero out u/v before adjoint calculation
!
!EOP
!-------------------------------------------------------------------------


      character(len=255) cwhat
      integer :: im, jm, km, jfirst, jlast, ng_d, ng_s
      integer :: i
      real(r8) :: pi,dl
      real(r8),allocatable :: pk(:,:,:)
      real(r8),allocatable :: ua(:,:,:)
      real(r8),allocatable :: va(:,:,:)
      real(r8),allocatable :: coslon(:)
      real(r8),allocatable :: sinlon(:)

      if ( vectype_/=5 ) return

!     Assume this is not a legitimate MPI program (meaning, run as serial)
!     -------------------------------------------
      cwhat  = lowercase(what)
      im     = w%grid%im
      jfirst = 1
      jm     = w%grid%jm
      jlast  = jm
      km     = w%grid%km
      ng_d   = 0
      ng_s   = 0

      allocate( pk(im,jm,km) )
      allocate( ua(im,jm,km) )
      allocate( va(im,jm,km) )
      allocate( coslon(im), sinlon(im) )

!     Define logitude at the center of the finite-volume i=1, Zamda = -pi
!     -------------------------------------------------------------------
      pi = 4.d0*datan(1.0d0)
      dl = (2*pi)/dble(im)
      do i=1,im/2
         coslon(i)      =  cos((i-1)*dl)
         coslon(i+im/2) = -coslon(i)
         sinlon(i)      =  sin((i-1)*dl)
         sinlon(i+im/2) = -sinlon(i)
      enddo

!     Convert winds from A- to D-grid
!     -------------------------------
      ua = w%u
      va = w%v
      if ( cwhat == 'adm' ) then
           w%u = 0.d0
           w%v = 0.d0
           call daInterp_vdtoa_ad ( w%u, w%v, ua, va )
      else
!          call a2d3d ( ua, va, w%u, w%v, &
!                       im, jm, km, jfirst, jlast, ng_d, ng_s, coslon, sinlon )
         call daInterp_vatod ( ua,va,w%u,w%v )
      endif
    
      deallocate( coslon, sinlon )
      deallocate( va )
      deallocate( ua )
      deallocate( pk )

      print *, 'Complete conversion of state from GEOS-5 to GEOS-4 type'

      end subroutine g5tog4_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: g4tog5_pert_ --- Converts G-4 pert to G-5 (NEVER make this public)
!
!
! !INTERFACE:
!
      subroutine g4tog5_pert_ ( im, jm, km, lm, n2d, n3d, ptop, bk, &
                                fld2d, fld3d, fnames, what )

! !USES:

      implicit none

! !INPUT PARAMETERS:
 
      integer,  intent(in) :: im,jm,km,lm,n2d,n3d
      real(r8), intent(in) :: ptop
      real(r8), intent(in) :: bk(km+1)
      character(len=3), intent(in) :: fnames(n2d+n3d)
      character(len=*), intent(in) :: what

! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: fld2d(im,jm)
      real(r8), intent(inout) :: fld3d(im,jm,km,n3d)

! !DESCRIPTION: Convert GEOS-4 perturbation vector to GEOS-5 like
!
! !REVISION HISTORY:
!
!  18Nov2007  Todling    Initial (temporary) code
!  18Dec2008  Errico/RT  Bug fix: zero out u/v before adjoint calculation
!  14Sep2013  Todling    A-grid winds can now be properly handled in norms
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*g4tog5_pert_'

      integer :: jfirst,jlast,ng_d, ng_s
      integer :: i,j,k,idd,idu,idv,idt
      real(r8),allocatable :: delp_ad(:,:,:) ! 
      real(r8),allocatable :: u_ad(:,:,:)    ! U-Wind on A-Grid
      real(r8),allocatable :: v_ad(:,:,:)    ! V-Wind on A-Grid
      real(r8),allocatable :: coslon(:)      ! cosine of longitudes
      real(r8),allocatable :: sinlon(:)      ! sine   of longitudes
      real(r8) pi,dl,bkweight
      character(len=255) :: cwhat

      if ( vectype_/=5 ) return

!     Assume this is not a legitimate MPI program (meaning, run as serial)
!     -------------------------------------------
      cwhat  = lowercase(what)
      jfirst = 1
      jlast  = jm
      ng_d   = 0
      ng_s   = 0

      allocate( delp_ad(im,jm,km),u_ad(im,jm,km),v_ad(im,jm,km))
      allocate( coslon(im), sinlon(im) )

!     Define logitude at the center of the finite-volume i=1, Zamda = -pi
!     -------------------------------------------------------------------
      pi = 4.d0*datan(1.0d0)
      dl = (2*pi)/dble(im)
      do i=1,im/2
         coslon(i)      =  cos((i-1)*dl)
         coslon(i+im/2) = -coslon(i)
         sinlon(i)      =  sin((i-1)*dl)
         sinlon(i+im/2) = -sinlon(i)
      enddo

      call find_name_ (n3d,fnames,'dp_',idd )
      call find_name_ (n3d,fnames,'u__',idu )
      call find_name_ (n3d,fnames,'v__',idv )
      delp_ad(:,:,:)=reshape(fld3d(:,:,:,idd),(/im,jm,km/))

!     Convert wind fields from D-grid to A-grid
!     ------------------------------------------
!_RT no longer need this, since norms now properly account for A-grid winds
!_RT  if ( cwhat == 'adm' ) then
!_RT       u_ad(:,:,:)=0.d0
!_RT       v_ad(:,:,:)=0.d0
!_RT       call daInterp_vatod_ad (u_ad,v_ad,fld3d(:,:,:,idu),fld3d(:,:,:,idv))
!_RT  else
!_RT       call daInterp_vdtoa (fld3d(:,:,:,idu),fld3d(:,:,:,idv),u_ad,v_ad)
!_RT  endif
!_RT
!_RT  fld3d(:,:,:,idu) = u_ad(:,:,:)
!_RT  fld3d(:,:,:,idv) = v_ad(:,:,:)

!     Also create gradient in ps based on the delp gradient 
!     -----------------------------------------------------
      fld2d(:,:) = 0._r8
      do k=1,km
         bkweight = ( bk(k+1) - bk(k) ) / ( bk(km) - bk(1) )
         do j=jfirst,jlast
            do i=1,im
               fld2d(i,j) = fld2d(i,j) + bkweight * delp_ad(i,j,k) ! this is an adj expression
            enddo
         enddo
      enddo

      deallocate( coslon, sinlon )
      deallocate( delp_ad, u_ad, v_ad )

      print *, myname_, 'Complete conversion of perturbation from GEOS-4 to GEOS-5 type'
      end subroutine g4tog5_pert_

      subroutine setparamC_ ( name, val )
      implicit none
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: val 
      if ( trim(name) == 'tname' ) tname_ = val(1:3)
      if ( trim(name) == 'wgrid' ) wgrid_ = val(1:1)
      end subroutine setparamC_
 
      subroutine setparamI_ ( name, val )
      implicit none
      character(len=*), intent(in) :: name
      integer,          intent(in) :: val 
      if ( trim(name) == 'vectype' ) vectype_ = val
      if ( trim(name) == 'vnorm' ) vnorm_ = val
      end subroutine setparamI_

      subroutine setparamL_ ( name, val )
      implicit none
      character(len=*), intent(in) :: name
      logical,          intent(in) :: val 
      end subroutine setparamL_

      subroutine setparamR_ ( name, val )
      implicit none
      character(len=*), intent(in) :: name
      real(r8),         intent(in) :: val 
      if ( trim(name) == 'eps_eer' ) eps_eer_ = val
      end subroutine setparamR_

      subroutine getparamC_ ( name, val )
      implicit none
      character(len=*), intent(in) :: name
      character(len=*), intent(out) :: val 
      if ( trim(name) == 'tname' ) val = tname_
      if ( trim(name) == 'wgrid' ) val = wgrid_
      end subroutine getparamC_
 
      subroutine getparamI_ ( name, val )
      implicit none
      character(len=*), intent(in) :: name
      integer,          intent(out) :: val 
      if ( trim(name) == 'vectype' ) val = vectype_
      if ( trim(name) == 'vnorm' ) val = vnorm_
      end subroutine getparamI_

      subroutine getparamL_ ( name, val )
      implicit none
      character(len=*), intent(in) :: name
      logical,         intent(out) :: val 
      val=.false.
      end subroutine getparamL_

      subroutine getparamR_ ( name, val )
      implicit none
      character(len=*), intent(in) :: name
      real(r8),         intent(out) :: val 
      if ( trim(name) == 'eps_eer' ) val = eps_eer_
      end subroutine getparamR_

      end module m_pertutil
