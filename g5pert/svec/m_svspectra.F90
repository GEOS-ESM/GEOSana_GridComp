!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_svspectra ---
!
! !INTERFACE:
                                                                                                       
      module m_svspectra

      use precision
                                                                                                       
      use m_const, only : Cp     => cpm   
      use m_const, only : R      => rgas  
      use m_const, only : kappa
      use m_const, only : pstd
      use m_const, only : tstd

      use m_pertutil, only : ptv2pt => pertutil_ptv2pt
      use m_pertutil, only : pt2t   => pertutil_pt2t
      use m_pertutil, only : pt2t_tl => pertutil_pt2t_tl
      use m_pertutil, only : calculate_ps => pertutil_calculate_ps
      use m_pertutil, only : pkappa => pertutil_pkappa
      use m_pertutil, only : horiz_grid => pertutil_horiz_grid
      use m_pertutil, only : vert_grid => pertutil_vert_grid
      use m_pertutil, only : find_name => pertutil_find_name
      use m_pertutil, only : transform_t_tl => pertutil_transform_t_tl

      use m_pertutil, only : e_cross_product 

      use m_mapz, only : set_eta
                                                                                                       
      use m_shtrans_DG
                                                                                                       
      use m_dyn

      use m_stdio
      use m_inpak90

      use m_ioutil, only : luavail
                                                                                                       
      implicit NONE
                                                                                                       
! !PUBLIC MEMBER FUNCTIONS:
                                                                                                       
      PRIVATE
                                                                                                       
      PUBLIC svspectra_Process
      PUBLIC svspectra_Clean
                                                                                                       
      interface svspectra_Process; module procedure svspectraProcess_; end interface
      interface svspectra_Clean; module procedure svspectraClean_; end interface
 
!
! !DESCRIPTION: Singular Vector Energy Calculation and Plotting routines
!
! !REVISION HISTORY:
!
!  25Aug2005   Elena N./Ron Gelaro  Coding based on SV similarity code
!                                   and on Ron's plotting program   
!  24Apr2006   Elena N.    Modified write into a binary file to work at PALM
!                          Added 72-levels pressure values for pl_vert_multi plotting routine
!  05Jan2008   Todling     Trying to clean up:
!                            - calc_ps using ptop now
!                            - using find_name w/ error check on
!  06Mar2009   Todling     .hdf to .nc4
!
!EOP
!-------------------------------------------------------------------------
 
      character(len=*), parameter :: myname = 'm_svEnergy_'
      real(r8), parameter :: pref = 100.d0*pstd
      real(r8), parameter :: tref = tstd
 
      CONTAINS
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: svspectra - compute vertical energy distribution and 
!                       energy spectra of two sets of singular vectors 
!
!
! !INTERFACE:
!
! NOTE to Ron E.: In SV code, 
!   check calc of E for v at NP:  should be cosp(j) in SV code
!   check  pe(:,:,k) = pe(:,:,k-1) + delp(:,:,k)
!
!
      subroutine svspectraProcess_  ( expid, nymd, nhms, svecnormI, nvecs, etype,    &
                                         isvec_fac, ifwd, ibak, res, TE_scale )
                                                            
      implicit NONE
                   
! !INPUT PARAMETERS:

      integer, intent(in) :: nvecs     ! Total number of vector pairs (ini - fin) to process 
      integer, intent(in) :: nymd(2)   ! Date of svecs: 1 - ini, 2 - final 
      integer, intent(in) :: nhms(2)   ! Time of svecs: 1 - ini, 2 - final 
      integer, intent(in) :: etype     ! define the energy type to be plotted
      integer, intent(in) :: ifwd, ibak    ! Number of hours of forward/backward integration  
      real, intent(in) :: isvec_fac ! Factor for ISVEC's to be plotted
      character(len=*), intent(in) :: expid     ! experiment id
      character(len=*), intent(in) :: svecnormI ! norm used to calculate SVecs ( 'ke', 'te', or else) 
      character(len=*), intent(in) :: res       ! resolution of the experiment
      logical, intent(in) :: TE_scale

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:
       
! !DESCRIPTION:
! compute vertical energy distribution
!    and energy spectra of two sets of singular vectors (pairs of initial and evolved vectors)

! !SEE ALSO:
!
       
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!  04Nov2004  Winslow    Placed into modular form
!  24Apr2006  Elena N.   Modified write into a binary file to work at PALM
!  17May2006  Elena N.   Added an option to scale evolved SVECs by teh global total energy norm
!                        for plotting purposes
!  05Jan2008  Todling    Updated e_cross_prod interface for unification of codes; 
!                        for that, add (k)domain and set manually to global
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: im1        ! grid points in zonal direction
      integer :: jn1        ! grid points in lat. direction
      integer :: nl1        ! grid points in vert. direction
      integer :: lm1        ! number of tracers from reference file(s)
!
!
! Specify the specific form of temperature to be examined.  This is the
! third argument in the character array fnames. The FVGCM field is 'VPT' 
! for virtual potential temperature. Instead, the derived potential 
! temperature will be examined if ' PT' is specified in its place.  
! '  T' means temperature itself will be examined in place of 'PVT'.
! ' ps' is the surface pressure derived from the model ' dp' (delp) field.
! The remaining values here generally should not be changed, but allow 
! easy generalization to include additional model or derived fields.
!
      integer, parameter :: nplot = 3     ! number of plots per page for GrADS script
      integer, parameter :: nfields3d=5   ! number of 3D fields to be examined
      integer, parameter :: nfields2d=1   ! number of 2D fields to be examined
      integer, parameter :: nfields=nfields3d+nfields2d

      character(len=3), dimension(nfields) :: fnames
      character(len=3), dimension(nfields), parameter :: fnames_def = &
                  (/ 'u__','v__','T__','q__','dp_','ps_' /)
!
      integer :: i,j,k,nf         ! looping indices

! Counters for do-loops and for other array indecies
      integer :: ndT, ndp, nps, kdp  ! indicates particular fields
      integer :: nvecA
      integer :: nvecB

      real(r8), allocatable  :: fields3dA(:,:,:,:,:)
      real(r8), allocatable  :: fields3dB(:,:,:,:,:)
      real(r8), allocatable  :: psfieldA(:,:)
      real(r8), allocatable  :: psfieldB(:,:)

      integer  :: ks                           ! level of strat
      real(r8), allocatable  :: jweights(:,:)  ! area weights (1) u, (2) T   
      real(r8), allocatable  :: kweights(:)    ! vertical layer weights
      real(r8), allocatable  :: glats(:,:)     ! degrees lats (1) u, (2) T   
      real(r8) :: ptop                         ! pressure at top of grid
      real(r8), allocatable  :: pmid(:)        ! p at data levels for standard grid
      real(r8), allocatable  :: delbk(:)       ! "thickness" of ps dependent portion of  
                                               ! standard (CCM) hybrid-coord. vert. grid
      real(r8), allocatable  :: ak(:)          ! hybrid coordinates
      real(r8), allocatable  :: bk(:)          ! hybrid coordinates

      integer :: kdomain(2)
      real(r8), allocatable :: domain(:,:,:)
      real(r8), allocatable :: energy_vert(:,:,:,:)
      real(r8), allocatable :: energy_spec(:,:,:,:)
      real(r8), allocatable :: V_norm(:)
      real(r8), allocatable  :: energy(:,:) 

      integer  :: kmax, mmax, nmax, nspects
      integer, allocatable   :: ntrunc(:)
      integer, allocatable   :: mtrunc(:)
      integer, allocatable   :: index(:,:)

      integer  :: ier
      integer  :: lunit1, lunit2, lunit3
      integer  :: e_index 
                                                                                                                                               
      real*4, allocatable :: aux_vert(:)
      real*4, allocatable :: aux_spec(:)
      real*4, allocatable :: growth(:)

      character(len=255) :: filename, svalu_file
      character(len=255) :: grads_string, header
! 
!                      BEGIN EXECUTION
!

      fnames = fnames_def
      call find_name (nfields3d,fnames,'dp_',ndp)

!    Form filename for initial reference state vector
!    -------------------------
     write(filename,'(2a,i8.8,a)')               &
              trim(expid), '.traj.lcv.',nymd(1),'.nc4'

!    Set some grid information from initial reference state
!    -------------------------

      call dyn_getdim   (filename,im1,jn1,nl1,lm1,ier)

      allocate (fields3dA(im1,jn1,nl1,nfields,2))
      allocate (fields3dB(im1,jn1,nl1,nfields,2))
      allocate (psfieldA(im1,jn1))
      allocate (psfieldB(im1,jn1))
      allocate (jweights(jn1,2)) ! area weights (1) u, (2) T
      allocate (kweights(nl1))   ! vertical layer weights
      allocate (glats(jn1,2))    ! degrees lats (1) u, (2) T
      allocate (pmid(nl1))       ! p at data levels for standard grid
      allocate (delbk(nl1))      ! "thickness" of ps dependent portion of
      allocate (ak(nl1+1))    
      allocate (bk(nl1+1))    
      allocate (energy(nl1+2,3))

!    Snag some grid information
!    -------------------------
      call horiz_grid (im1, jn1, jweights, glats)
      call vert_grid  (nl1, ks, ak, bk, kweights, ptop, pmid, delbk)

!    Determine number of modes for spectral analysis (corresponding to D - grid)
!    -------------------------
                                                                                                                                               
      mmax = im1/2 - 1
      nmax = im1/2 - 1
      kmax = im1/2 - 1
                                                                                                                                               
      allocate (ntrunc(0:mmax), mtrunc(0:nmax), index(0:mmax,0:nmax))
                                                                                                                                               
      call sh_initial_D (nmax,mmax,kmax,im1,jn1,nspects,   &
                           ntrunc,mtrunc,index,'both',.true.,  &
                           .true.,ier)
!
      if (ier /= 0) print *,' sh_initial_D returns ierror=',ier

      allocate(energy_vert(nl1+2,3,nvecs,2))
      allocate(energy_spec(0:nmax,3,nvecs,2))
!
      allocate(aux_vert(nl1))
      allocate(aux_spec(nmax+1))
      allocate(growth(nvecs+1))
      allocate(V_norm(nvecs))
      
!    Read reference state into arrays fields3dA(:,:,:,:,2) and fields3dB(:,:,:,:,2) 
!    -------------------------

      call read_data2_ (im1,jn1,nl1,lm1,nfields,fnames, &
                        trim(filename),nhms(1),nymd(1), &
                        fields3dA(:,:,:,:,2))

!    Form filename for final reference state vector
!    -------------------------

     write(filename,'(2a,i8.8,a,i2.2,2a,i3.3,a)')               &
              trim(expid), '.traj.lcv.',nymd(2),'.nc4'

      call read_data2_ (im1,jn1,nl1,lm1,nfields,fnames, &
                        trim(filename),nhms(2),nymd(2), &
                        fields3dB(:,:,:,:,2))

       lunit1 = luavail()
       open(lunit1,file='evert.dat',form='unformatted',          &
             access='direct', recl=nl1)

       lunit2 = luavail()
       open(lunit2,file='espec.dat',form='unformatted',          &
             access='direct', recl=(nmax+1))
!      
!   Loop over set of vectors A
!    -------------------------
      aux_spec(:) = 0.d0
      aux_vert(:) = 0.d0
      

!
!   Define energy type to be plotted
!    -------------------------
     if (etype.eq.1) e_index = 3
     if (etype.eq.2) e_index = 1
     if (etype.eq.3) e_index = 2
     if (etype.gt.3) e_index = 3

     allocate(domain(im1,jn1,nl1) )
     domain = 1.0d0
     kdomain(1)=1; kdomain(2)=nl1
 
      do nvecA=1,nvecs

!    Form filename for initial singular vector name
!    -------------------------
     write(filename,'(2a,i8.8,a,i2.2,2a,i3.3,a)')               &
              trim(expid), '.isvec.eta.',nymd(1),'_',nhms(1),'z.',trim(svecnormI),nvecA,'.nc4'


!
!    Read vector number nvecA from set A into fields3dA(:,:,:,:,1)  
!    -------------------------

        call read_data2_ (im1,jn1,nl1,lm1,nfields,fnames, &
                          trim(filename),nhms(1),nymd(1),      &
                          fields3dA(:,:,:,:,1))
!
!    Transform T fields 
!    -------------------------
        call transform_T_tl (im1,jn1,nl1,ptop,nfields,fnames, &
                                  fields3dA(:,:,:,:,1:2) )
!
!    Compute ps fields
!    -------------------------
        call calculate_ps (im1,jn1,nl1,ptop,fields3dA(:,:,:,ndp,1),psfieldA(:,:))

!    Calculate spectral energy and store in 'energy_vert' and 'energy_spec', respectively
!    -------------------------

                call compute_spectra_ (im1,jn1,nl1,nfields3d,       &
                   kweights,mmax, nmax, kmax, nspects,    &
                   ntrunc, mtrunc, index,                        &
                   fields3dA(:,:,:,:,1),fields3dA(:,:,:,ndp,2),     &
                   psfieldA(:,:),fnames(1:nfields3d),               &
                   energy_vert(:,:,nvecA,1), energy_spec(:,:,nvecA,1))

	       do j = 1, nl1 
        	   aux_vert(j) = energy_vert(j,e_index,nvecA,1) 
	       enddo
!
  	       write (lunit1, rec=nvecA) (aux_vert(j), j = 1, nl1)
!
!    Output spectral energy data of an initial vector in a binary file 
!    -------------------------
	       do j = 0, nmax   ! Start from first harmonic
        	   aux_spec(j+1) = energy_spec(j,e_index,nvecA,1)
	       enddo
!
	        write(lunit2, rec=nvecA) (aux_spec(j), j=1, nmax+1)

       enddo    ! End loop for initial SVECs
!    -------------------------
!
!    Writing means energy values for initial vectors
!    -------------------------
       do j = 1, nl1
           aux_vert(j) = sum(energy_vert(j,e_index,1:nvecs,1))/dble(nvecs)
       enddo
!
       write (lunit1, rec=nvecs+1) (aux_vert(j), j = 1, nl1)
!
       do j = 0, nmax   ! Start from first harmonic
            aux_spec(j+1) = sum(energy_spec(j,e_index,1:nvecs,1))/dble(nvecs)
       enddo

       write(lunit2, rec=nvecs+1) (aux_spec(j), j=1, nmax+1)
!
!    Loop over set of evolved vectors B
!    -------------------------
       do nvecB=1,nvecs

!    Form filename for final singular vector name
!    -------------------------
     write(filename,'(2a,i8.8,a,i2.2,2a,i3.3,a)')               &
              trim(expid), '.fsvec.eta.',nymd(2),'_',nhms(2),'z.',trim(svecnormI),nvecB,'.nc4'
!
!    Read vector number nvecB from set B into fields3dB(:,:,:,:,1)  
!    -------------------------
         call read_data2_ (im1,jn1,nl1,lm1,nfields,fnames, &
                          trim(filename),nhms(2),nymd(2),      &
                          fields3dB(:,:,:,:,1))
!
!    Transform T fields 
!    -------------------------
         call transform_T_tl (im1,jn1,nl1,ptop,nfields,fnames, &
                                  fields3dB(:,:,:,:,1:2) )
!
!    Compute ps fields
!    -------------------------
         call calculate_ps (im1,jn1,nl1,ptop,fields3dB(:,:,:,ndp,1),  &
                             psfieldB(:,:))
!
        call e_cross_product (im1,jn1,nl1,nfields,ptop,       &
           domain, domain, domain, kdomain, reshape(jweights(:,2),(/jn1/)),      &
           reshape(fields3dB(:,:,:,:,1),(/im1,jn1,nl1,nfields/)),   &
           reshape(fields3dB(:,:,:,:,1),(/im1,jn1,nl1,nfields/)),   &
           psfieldB, psfieldB, &
           reshape(fields3dB(:,:,:,ndp,2),(/im1,jn1,nl1/)),    &
           fnames(1:nfields3d), energy )

        V_norm(nvecB) = 1.d0/sqrt(energy(nl1+2,3))

! Normalize FSVECS if requested

        if ( TE_scale ) fields3dB(:,:,:,:,1) = fields3dB(:,:,:,:,1)*V_norm(nvecB)

              call compute_spectra_ (im1,jn1,nl1,nfields3d,     &
                   kweights, mmax, nmax, kmax, nspects,    &
                   ntrunc, mtrunc, index,                    &    
                   fields3dB(:,:,:,:,1),fields3dB(:,:,:,ndp,2), &
                   psfieldB(:,:),fnames(1:nfields3d),           &
                   energy_vert(:,:,nvecB,2), energy_spec(:,:,nvecB,2))

!
	       do j = 1, nl1
        	   aux_vert(j) = energy_vert(j,e_index,nvecB,2)
	       enddo
!
	       write (lunit1, rec=nvecs+1+nvecB) (aux_vert(j), j = 1, nl1)
!
!    Plot vertical profiles
!    -------------------------
	     do j = 0, nmax   ! Start from first harmonic
          	 aux_spec(j+1) = energy_spec(j,e_index,nvecB,2)
	     enddo
!
	     write (lunit2, rec=nvecs+1+nvecB) (aux_spec(j), j=1, nmax+1)             
!
      enddo   ! End loop over final SVEC's 
!    -------------------------
!
!     Writing mean energy values for final vectors
!    -------------------------
!
       do j = 1, nl1
           aux_vert(j) = sum(energy_vert(j,e_index,1:nvecs,2))/dble(nvecs)
       enddo
!
       write (lunit1, rec=2*nvecs+2) (aux_vert(j), j = 1, nl1)
!
       do j = 0, nmax   ! Start from first harmonic
            aux_spec(j+1) = sum(energy_spec(j,e_index,1:nvecs,2))/dble(nvecs)
       enddo
!
       write(lunit2, rec=2*nvecs+2) (aux_spec(j), j=1, nmax+1)
!

      close(lunit1)
      close(lunit2)
!    -------------------------
! ---- RON's GRAPHICS ROUTINES
! ----
! ---- When ready, growth(1:nvecs) will be read from file, and growth(nvecs+1) will be their mean
! ----
!    -------------------------                             
      growth(:) = 1.0
      growth(nvecs+1)= 0.0

      write(svalu_file,'(2a,i8.8,a)') trim(expid), '.svalu.', nymd(1), '.txt'
      lunit3 = luavail()
      open (lunit3, file = svalu_file, form='formatted')
      read (lunit3, *) header
      do j=1,nvecs
        read (lunit3, '(i4,e14.6)') i, growth(j) 
        growth(j) = sqrt(growth(j))
        growth(nvecs+1) = growth(nvecs+1)+growth(j)
      enddo
      close (lunit3)

      growth(nvecs+1)= growth(nvecs+1)/nvecs

      write(grads_string,'(i8.8,i2.2)')               &
              nymd(1), nhms(1)/10000
                                                                                                                                               
!    Plot vertical profiles
!    -------------------------
 
      call pl_vert_multi(nvecs+1,nplot,nl1,etype,isvec_fac,grads_string,    &
                          trim(expid),.true., growth, ifwd, ibak, res)

!    Plot spectral energy distribution
!    -------------------------
                                                                                                                                               
      call pl_spec_multi(nvecs+1,nplot,nmax+1,etype,isvec_fac,grads_string,  &
                          trim(expid),1,.true., growth, ifwd, ibak, res)

!    Deallocate memory 
!    -------------------------
!
      call sh_clean
      
      deallocate ( domain )
      deallocate (growth, V_norm)
      deallocate (ntrunc, mtrunc, index, stat = ier)
      deallocate (fields3dA,fields3dB)
      deallocate (psfieldA,psfieldB)
      deallocate(jweights,kweights,glats,pmid,delbk,ak,bk)
      deallocate(energy_vert, energy_spec)
      deallocate (aux_vert, aux_spec)
      deallocate ( energy)
!
      print *,' '
      print *,' '
      print *,' '
      print *,' PROGRAM COMPLETED'
!
      end subroutine svspectraProcess_  
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: compute_spectra --- Compute contributions to the total energy norm E.
!
!
! !INTERFACE:
!
      subroutine compute_spectra_ (im1,jn1,nl1,nfields,    &
                   kweights, mmax, nmax, kmax, nspects,    &    
                   ntrunc, mtrunc, index,                  &
                   fields,delp_ref,psfield,fnames,e_vert,e_spec) 

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1,nfields
      character(len=*), intent(in) :: fnames(nfields) 
      real(r8), intent(in)  :: kweights(jn1)
      real(r8), intent(in)  :: fields(im1,jn1,nl1,nfields)
      real(r8), intent(in)  :: delp_ref(im1,jn1,nl1)
      real(r8), intent(in)  :: psfield(im1,jn1)   
      integer, intent(in)   :: nspects, mmax, nmax, kmax
      integer, intent(in)   :: ntrunc(0:mmax)
      integer, intent(in)   :: mtrunc(0:nmax)
      integer, intent(in)   :: index(0:mmax,0:nmax)

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: e_vert(nl1+2,3)  
      real(r8), intent(out) :: e_spec(0:nmax,3) 

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
! 
! Compute contributions to the total energy norm E.  This is the same E as
! defined for the FVGCM Singular Vector calculations. Separate contributions
! by kinetic and potential energies for each pressure layer are determined. 
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm to compute energy norm
!  15Apr2005  Elena N.   Added calls to spectral analysis routines
!                        and put them in this separate subroutine   
!  05Jul2005  R. Errico  Modified spectral power routines: 
!                        their calls have more input values
!                        to allow for different truncation 

!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer  :: i,j,k
      integer  :: idu,idv,idt
      real(r8) :: tfac, ufac, pfac
      real(r8) :: ps_ref(im1,jn1)
!                                                                 
      real(r8), allocatable     :: u(:,:,:)
      real(r8), allocatable     :: v(:,:,:)
      real(r8), allocatable     :: T(:,:,:)
      real(r8), allocatable     :: Ps(:,:,:)
 
      integer  :: m,n    ! m=zonal wavenumber, n=n-m
      integer  :: ier
      integer  :: level
!      integer  :: nspects
!      integer  :: mmax
!      integer  :: nmax 
!      integer  :: kmax 
!      integer, allocatable   :: ntrunc(:)
!      integer, allocatable   :: mtrunc(:)
!      integer, allocatable   :: index(:,:)
 
      real(r8), allocatable :: powerT(:,:,:)
      real(r8), allocatable :: powerV(:,:,:)
      real(r8), allocatable :: powerPs(:,:,:)
 
      real(r8), allocatable  :: powerTOTAL(:,:,:)
 
 
      complex(r8), allocatable :: scoefsT(:,:,:)
      complex(r8), allocatable :: scoefsV(:,:,:)
      complex(r8), allocatable :: scoefsPs(:,:,:)
 
      real(r8) :: total_sum

! One factor of 0.5 is for defenition of energy as 1/2 * ...
! Assumption here is that sum(jweights)=1
      tfac=0.5d0*Cp/tref
      ufac=0.5d0
      pfac=0.5d0*R*tref/pref**2

      call find_name (nfields,fnames,'u__',idu)      
      call find_name (nfields,fnames,'v__',idv)      
      call find_name (nfields,fnames,'T__',idt)      
!
! ------- Compute PS of the referenced state 

      ps_ref(:,:) = sum(delp_ref(:,:,:),3)       
!
      allocate ( u(im1,jn1,nl1), v(im1,jn1,nl1), &
          T(im1,jn1,nl1), Ps(im1,jn1, 1), stat = ier)
!
      if (ier /= 0) print *, 'Cannot allocate u,v,T,dsigma etc' 

!      mmax = im1/2
!      nmax = im1/2
!      kmax = im1/2

      allocate ( powerT(0:nmax,nl1,1), stat = ier)
      allocate ( powerV(0:nmax,nl1,2), stat = ier)
      allocate ( powerPs(0:nmax,1,1), stat = ier)
!
      if (ier /= 0) print *, 'Cannot allocate powerSpec'
 
      allocate ( powerTOTAL(0:nmax,nl1,1), stat = ier)
 
      u(:,:,:) = fields(:,:,:,idu)
      v(:,:,:) = fields(:,:,:,idv)
      T(:,:,:) = fields(:,:,:,idt)
      Ps(:,:,1) = psfield(:,:)

! ------- Find nspects for grid_D with known im1, jn1, nl1 values:
!
!      allocate (ntrunc(0:mmax), mtrunc(0:nmax), index(0:mmax,0:nmax))
!                                         
!      call sh_initial_D (nmax,mmax,kmax,im1,jn1,nspects,   &
!                           ntrunc,mtrunc,index,'both',.true.,  &
!                           .true.,ier)
!                                    
!      if (ier /= 0) print *,' sh_initial_D returns ierror=',ier
!
! ------- Now that nspects is known for this truncation, allocate spectral coefs
! 
      allocate (scoefsT(nspects,nl1,1), stat = ier)
      allocate (scoefsV(nspects,nl1,2), stat = ier)
      allocate (scoefsPs(nspects,1, 1), stat = ier)
!
      if (ier /= 0) print *, 'Cannot allocate scoefs'
!
!
! ------- Initialize coefs with zeros
!
      scoefsT(:,:,:)=cmplx(0.d0, 0.d0, r8)
      scoefsV(:,:,:)=cmplx(0.d0, 0.d0, r8)
      scoefsPs(:,:,:)=cmplx(0.d0, 0.d0, r8)
!
      call sh_transforms (im1,jn1,mmax,nmax,kmax,nl1, &
                                nspects,ntrunc,index,       &
                                1,T,scoefsT,'S','f2s',ier)
!
      if (ier /= 0) print *,' ier=',ier,'Error returned from T sh_transform'
!

      call sh_transforms (im1,jn1,mmax,nmax,kmax,nl1, &
                                nspects,ntrunc,index,       &
                                2,u,scoefsV,'U','f2s',ier)
!
      if (ier /= 0) print *,' ier=',ier,'Error returned from U sh_transform'
!        
      call sh_transforms (im1,jn1,mmax,nmax,kmax,nl1, &
                                nspects,ntrunc,index,       &
                                2,v,scoefsV,'V','f2s',ier)
!
      if (ier /= 0) print *,' ier=',ier,'Error returned from V sh_transform'
!
      call sh_transforms (im1,jn1,mmax,nmax,kmax,1, &
                                nspects,ntrunc,index,       &
                                1,Ps,scoefsPs,'S','f2s',ier)
!
      if (ier /= 0) print *,' ier=',ier,'Error returned from Ps sh_transform'
!
! ------- Compute power spectra from complex coefficients
!
      powerT(:,:,1) = 0.d0 
      powerV(:,:,1:2) = 0.d0 
      powerPs(:,1,1) = 0.d0
! 
! ------- NEW VERSION OF SH_CALC_POWER IS USED (07/18/2005)
!
      call sh_calc_power (mmax,nmax,kmax,im1,nspects,nl1,1,    &
                          mtrunc,index,powerT(:,:,1:1),scoefsT,'S')
      call sh_calc_power (mmax,nmax,kmax,im1,nspects,nl1,2,    &
                          mtrunc,index,powerV(:,:,1:2),scoefsV,'VV')
      call sh_calc_power (mmax,nmax,kmax,im1,nspects,1,1,          &
                          mtrunc,index,powerPs(:,1:1,1:1),scoefsPs,'S')
 
!----------------------------------------------------------------------
!
      total_sum = 0.0
      e_vert(:,:) = 0.d0
                                                                 
      do level=1,nl1
                                                                 
          powerTOTAL(0:nmax,level,1) = kweights(level)*              &
                                  (ufac*(powerV(0:nmax,level,1) +    &
                                   powerV(0:nmax,level,2)) +         &
                                   tfac*powerT(0:nmax,level,1)  +    &  
                                   pfac*powerPs(0:nmax,1,1)) 

          total_sum = total_sum + sum(powerTOTAL(:,level,1))

          e_vert(level,1) = kweights(level)*ufac*(sum(powerV(0:nmax,level,1)) + sum(powerV(0:nmax,level,2)))

          e_vert(level,2) = kweights(level)*tfac*sum(powerT(0:nmax,level,1)) 

          e_vert(level,3) = e_vert(level,1) + e_vert(level,2) 

          e_vert(nl1+1,:) = e_vert(nl1+1,:) + e_vert(level,:)
                                                                 
     enddo

     e_vert(nl1+2,3) = e_vert(nl1+1,3) + pfac*sum(powerPs(0:nmax,1,1))

     e_spec(:,:) = 0.d0

     do j = 0, nmax   ! Start from first harmonic

          e_spec(j,1) = ufac*(sum(kweights(:)*powerV(j,:,1))+sum(kweights(:)*powerV(j,:,2)))
          e_spec(j,2) = tfac*sum(kweights(:)*powerT(j,:,1))
          e_spec(j,3) = sum(powerTOTAL(j,:,1))  

     enddo

! ------- Call clean to deallocate saved module internal arrays
                                                                 
!      call sh_clean
                                                                 
      deallocate (scoefsT,scoefsV,scoefsPs, stat = ier)
      deallocate (u, v, T, Ps, stat = ier)
      deallocate (powerT, powerV, powerPs, powerTOTAL, stat = ier)
!      deallocate (ntrunc, mtrunc, index, stat = ier)

! 
      if (ier /= 0) print *, 'Cannot deallocate u,v,T,dsigma etc'
!
      end subroutine compute_spectra_ 
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: read_data2 --- Read in GEOS 4 vector for placement into stats
!                          configuration
!
!
! !INTERFACE:
!
      subroutine read_data2_ (im4,jn4,nl4,lm4,nfields,fnames, &
                            job,nhms,nymd,fields3d)
                                                                                                         
!USES:
                                                                                                         
      implicit none
                                                                                                         
                                                                                                         
! !INPUT PARAMETERS:
                                                                                                         
      integer,  intent(in)  :: im4, jn4, nl4, lm4, nfields
      character(len=*), intent(in) :: fnames(nfields)
      character(len=*), intent(in) :: job
      integer , intent(in)  :: nhms, nymd
                                                                                                         
! !OUTPUT PARAMETERS:
                                                                                                         
      real(r8), intent(out) :: fields3d(im4,jn4,nl4,nfields)
                                                                                                         
! !INPUT/OUTPUT PARAMETERS:
                                                                                                         
! !DESCRIPTION:
             
! !SEE ALSO:
!
            
! !REVISION HISTORY:
!
!  04Oct2004  Winslow    Original code
!  15Apr2005  Elena N.   Change if/else options for T/PT/VPT read to ensure
!                        that a corresponding fname label is found 
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables
             
      integer :: ier, nhms_in, nymd_in, indexF
      type(dyn_vect) w
             
      nhms_in = nhms
      nymd_in = nymd
      call dyn_null ( w )
      call dyn_get ( trim(job), nymd_in, nhms_in, w, ier, timidx=0 )
      if (w%grid%im.ne.im4 .or. w%grid%jm.ne.jn4 .or. w%grid%km.ne.nl4 ) then
        print*,'***ERROR READING DYNAMIC VECTOR***'
        print*,'read_data2: dimensions do not agree'
        print*,w%grid%im,w%grid%jm,w%grid%km
        print*,im4,jn4,nl4
        print*,'***EXIT ON READING ERROR***'
        stop
      endif

      if (ier.ne.0) then
        print*,'failed to read data from ',trim(job)
      endif
             
             
      call find_name (nfields,fnames,'u__',indexF)
      if (indexF .eq. 0) then
        print*,'Error read_data2: did not fill u'
        stop
      else
        fields3d(:,:,:,indexF)=reshape(w%u,(/w%grid%im,w%grid%jm,w%grid%km/))
      endif
             
      call find_name (nfields,fnames,'v__',indexF)
      if (indexF .eq. 0) then
        print*,'Error read_data2: did not fill v'
        stop
      else
        fields3d(:,:,:,indexF)=reshape(w%v,(/w%grid%im,w%grid%jm,w%grid%km/))
      endif
             
      call find_name (nfields,fnames,'T__',indexF)
      if (indexF .eq. 0) then
        call find_name (nfields,fnames,'vpT',indexF)
        if (indexF .eq. 0) then
          call find_name (nfields,fnames,'PT_',indexF)
          if (indexF .eq. 0) then
            print*,'Error read_data: did not read T'
            stop
          else
            fields3d(:,:,:,indexF)=reshape(w%pt,(/w%grid%im,w%grid%jm,w%grid%km/))
          endif
        else
          fields3d(:,:,:,indexF)=reshape(w%pt,(/w%grid%im,w%grid%jm,w%grid%km/))
        endif
      else
      fields3d(:,:,:,indexF)=reshape(w%pt,(/w%grid%im,w%grid%jm,w%grid%km/))
      endif

      call find_name (nfields,fnames,'q__',indexF)
      if (indexF .eq. 0) then
        print*,'Error read_data2: did not fill q'
        stop
      else
        fields3d(:,:,:,indexF)=reshape(w%q(:,:,:,1),(/w%grid%im,w%grid%jm,w%grid%km/))
      endif
             
      call find_name (nfields,fnames,'dp_',indexF)
      if (indexF .eq. 0) then
        print*,'Error read_data2: did not fill dp'
        stop
      else
        fields3d(:,:,:,indexF)=reshape(w%delp,(/w%grid%im,w%grid%jm,w%grid%km/))
      endif

      call dyn_clean ( w )
      call dyn_null ( w )
             
      end subroutine read_data2_
             
      subroutine svspectraClean_ ()
      end subroutine svspectraClean_

   end module m_svspectra 
