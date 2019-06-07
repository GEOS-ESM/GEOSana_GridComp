! -------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_pertenergy --- Computes energy field of a perturbation vector 
!
! !INTERFACE:

      module m_pertenergy

      use precision

      use m_const, only : Cp     => cpm
      use m_const, only : R      => rgas
      use m_const, only : pstd    
      use m_const, only : alhl    
      use m_const, only : tref   => tstd
      use m_const, only : kappa  

      use m_pertutil, only : ptv2pt => pertutil_ptv2pt
      use m_pertutil, only : pt2t => pertutil_pt2t
      use m_pertutil, only : pt2t_tl => pertutil_pt2t_tl
      use m_pertutil, only : calculate_ps => pertutil_calculate_ps
      use m_pertutil, only : pkappa => pertutil_pkappa
      use m_pertutil, only : find_name => pertutil_find_name
      use m_pertutil, only : horiz_grid => pertutil_horiz_grid
      use m_pertutil, only : vert_grid => pertutil_vert_grid
      use m_pertutil, only : transform_t_tl => pertutil_transform_t_tl
      use m_pertutil, only : read_data => pertutil_read_data
      use m_pertutil, only : getparam  => pertutil_getparam
      use m_pertutil, only : setparam  => pertutil_setparam

      use m_GFIO_PutFld,only : GFIO_PutFld

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type

      use m_stdio, only : stdout
      use m_stdio, only : stderr

      use m_strTemplate, only : strTemplate
      use m_inpak90

      use m_die, only: MP_die, die

      use m_mapz, only : set_eta

      use m_dyn


      implicit NONE
                                                                                                      
! !PUBLIC MEMBER FUNCTIONS:
                                                                                                   
      PRIVATE
                                                                                                           
      PUBLIC pertenergy_Set
      PUBLIC pertenergy_Process
      PUBLIC pertenergy_Clean
                                                                                                          
      interface pertenergy_Set; module procedure pertengy_set_; end interface
      interface pertenergy_Process; module procedure energy_field_Process_; end interface
      interface pertenergy_Clean; module procedure StatsClean_; end interface
                                                                                                      
!
! !DESCRIPTION: Computes energy field of a perurbation vector 
!
! !REVISION HISTORY:
!
!  18Oct2004  Winslow   Initial code based on similar code
!  10May2005  Elena N.  Linked constants of this module to common constants 
!  24Apr2006  Elena N.  Added 72-lev pressure info for Grads output
!  01May2006  Elena N.  Adjusted Grads output options to work at PALM/Altix
!  03May2007  R.Gelaro  Added rcfile capability to process multiple perturbations
!  05Jan2008  Todling   Clean up:
!                        - calc_ps using ptop now 
!                        - find_name doing error check
!                        - revisited energy calculation
!                        - replaced output from binary grads to nc4 file
!  22Jan2012  Todling   Allow for reset of vertical norm (weights)
!
!EOP
!-------------------------------------------------------------------------
 
      character(len=*), parameter :: myname = 'm_pertenergy'

      real(r8), parameter :: pref = 100.d0*pstd 

      logical, save :: readRC_ = .false.

!     Define setup parameters for pertenergy calculation
!     --------------------------------------------------
       character(len=255), save :: pertmpl  ! template for pert filenames
       character(len=255), save :: reftmpl  ! template for ref state filenames

       integer, save :: nvecs        ! number of perts (e.g. SVs) per case/date
       integer, save :: tinc         ! time icrement between cases/dates (sec)
       integer, save :: ifwd         ! length of forward model integration (sec)
       integer, save :: ibak         ! length of adjoint model integration (sec)

       integer, parameter :: nkxe = 1   ! index for kinetic energy contribution
       integer, parameter :: nape = 2   ! index for avail. pot. energy contribution
       integer, parameter :: ngpe = 3   ! index for geopotential. energy contribution
       integer, parameter :: nqxe = 4   ! index for moist energy contribution
       integer, parameter :: ntxe = 5   ! index for total dry energy
       integer, parameter :: ntwe = 6   ! index for total wet energy
       integer, parameter :: nall = 6   ! max index from those above

       character(len=*), dimension(nall), parameter :: vnames = &
                          (/'kxe','ape','pse','qxe','txe','twe'/)
       character(len=*), dimension(nall), parameter :: vunits = &
                          (/'J/Kg','J/Kg','J/Kg','J/Kg','J/Kg','J/Kg'/)
       character(len=*), dimension(nall), parameter :: vtitles = &
                           (/'Kinetic Energy            ', &
                             'Available Potential Energy', &
                             'Potential Energy',           &
                             'Latent Heat Energy        ', &
                             'Total Dry Energy          ', &
                             'Total Wet Energy          '/)

      CONTAINS
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: energy_field_Process_ --- Compute basic stats and do T transform
!
!
! !INTERFACE:
!
      subroutine energy_field_Process_ ( nfiles, dynfiles, expid, notag, &
                                         nymd_beg, nhms_beg, nymd_end, nhms_end )
 

!USES:
                                                                                                                                               
      implicit none
 
! !INPUT PARAMETERS:
 
      integer, intent(in) :: nfiles
      character(len=*), intent(inout) :: dynfiles(nfiles)
      character(len=*), intent(in)    :: expid
      logical, intent(in)    :: notag
      integer, intent(inout) :: nymd_beg
      integer, intent(inout) :: nhms_beg
      integer, intent(inout) :: nymd_end
      integer, intent(inout) :: nhms_end

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!  16Sep2005  Elena N./RG Adopted m_pertenergy.F90 to process perturbation vector
!                         and to produce energy fields
!  10Jan2008  Todling    - Replaced binary grads output with nc4 output
!                        - Unification of post-proc routines
!  16Mar2012  Todling    - handle d/a grids
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer, parameter :: nfields3d=5   ! number of 3D fields to be examined
      integer, parameter :: nfields2d=1   ! number of 2D fields to be examined
      integer, parameter :: nfields=nfields3d+nfields2d
      character(len=3), dimension(nfields) :: fnames
      character(len=3), dimension(nfields), parameter :: fnames_def = &
                  (/ 'u__','v__','T__','q__','dp_','ps_' /)
      integer :: im1   ! grid points in zonal direction from reference file(s)
      integer :: jn1   ! grid points in lat. direction from reference file(s)
      integer :: nl1   ! grid points in vert. direction from reference file(s)
      integer :: lm1   ! number of tracers from reference file(s)
      integer :: im2,jn2,nl2,lm2

      logical :: first ! flag for first pass through time loop

! Counters for do-loops and for other array indicies
      integer :: ndT, ndp, nps, kdp  ! indicates particular fields
      character(len=255) :: energy_filename, filename
      character(len=3)   :: tname

      real(r8), allocatable :: fields2d(:,:,:,:,:)
      real(r8), allocatable :: fields3d(:,:,:,:,:)

      real(r8), allocatable :: energy          (:,:)
      real(r8), allocatable :: energy_field    (:,:,:,:)
      real(r8), allocatable :: energy_field_ave(:,:,:,:)
       
      real(r8), allocatable :: jweights(:,:) ! area weights (1) u, (2) T   
      real(r8), allocatable :: kweights(:)    ! vertical layer weights
      real(r8), allocatable :: glats(:,:)    ! degrees lats (1) u, (2) T   
      real(r8) :: ptop            ! pressure at top of grid
      real(r8), allocatable :: pmid(:)        ! p at data levels for standard grid
      real(r8), allocatable :: delbk(:)       ! "thickness" of ps dependent portion of  
                                  ! standard (CCM) hybrid-coord. vert. grid
      integer  :: ks              ! index for strat
      real(r8), allocatable :: ak(:)        ! values for vertical mapping
      real(r8), allocatable :: bk(:)        ! values for vertical mapping

      integer :: nymd,nhms
      integer :: ncount          ! counts total perts processed (nvecs * dates/times)
      integer :: ier
      logical :: pick

! 
!                      BEGIN EXECUTION
!
      fnames = fnames_def
      call find_name (nfields3d,fnames,'T__',ndt)
      call getparam ( 'tname', tname )

      call find_name (nfields3d,fnames,'dp_',ndp)

      nps = nfields2d 

      first = .true.
      ncount = 0

      nymd = nymd_beg
      nhms = nhms_beg
      if(nymd==nhms .and. nymd==0) pick = .false.
      if(.not.readRC_) tinc = 999999
!      
! Loop over time
!
      do while ( nymd.le.nymd_end  .and.  nhms.le.nhms_end )

         if (readRC_) call fill_templ ( dynfiles, expid, nfiles, nymd, nhms )

         if ( first ) then 

!          Read in dimensions of the reference field
           call dyn_getdim (trim(dynfiles(1)),im1,jn1,nl1,lm1,ier)
           call dyn_getdim (trim(dynfiles(2)),im2,jn2,nl2,lm2,ier)
          
             if (nl1.ne.nl2.or.im1.ne.im2.or.jn1.ne.jn2) then
                 call die(myname,'Resolution of pert inconsistent w/ that of ref state')
             endif
 
!          Compute weighting terms

           allocate (jweights(jn1,2))
           allocate (kweights(nl1))
           allocate (glats(jn1,2))
           allocate (pmid(nl1))
           allocate (delbk(nl1))
           allocate (ak(nl1+1))
           allocate (bk(nl1+1))

           call horiz_grid (im1, jn1, jweights, glats)
           call vert_grid  (nl1, ks, ak, bk, kweights, ptop, pmid, delbk)

           allocate(energy                   (nl1+2,3))
           allocate(energy_field    (im1,jn1,nl1,nall))
           allocate(energy_field_ave(im1,jn1,nl1,nall))

           allocate(fields3d(im1,jn1,nl1,nfields3d,2))
           allocate(fields2d(im1,jn1,  1,nfields2d,2))
      
           energy_field    = 0.d0
           energy_field_ave= 0.d0
  
           first = .false. 

         endif

! Read data from reference state
! NOTE: For some unknown reason, this program gives a memory corruption error
!       when reading fields directly in their proper slot - so the first slot
!       is used to read both perturbation and ref state and the result is 
!       placed in proper slot (Todling)
 
         fnames(ndt) = tname
         print *, 'Temperature field on input set to: ', tname

         call read_data (im1,jn1,nl1,nfields,fnames,nymd,nhms, &
                         trim(dynfiles(2)),fields3d(:,:,:,:,1:1),'nlm',ier)
         fields3d(:,:,:,:,2:2)=fields3d(:,:,:,:,1:1)
         if (.not.pick .and. .not. readRC_) then
              nymd_end = nymd
              nhms_end = nhms
              pick = .true.
         endif

! Read data from a perturbation vector 
         call read_data (im1,jn1,nl1,nfields,fnames,nymd,nhms, &
                         trim(dynfiles(1)),fields3d(:,:,:,:,1:1),'tlm',ier)
!
         call transform_T_tl (im1,jn1,nl1,ptop,nfields3d,fnames(1:nfields3d), &
                              fields3d(:,:,:,:,1:2) )

! Compute ps fields
         call calculate_ps (im1, jn1, nl1, ptop, fields3d(:,:,:,ndp,1),  & 
                                                 fields2d(:,:,1,nps,1))

! reference fields always follow perturbed fields
         call calculate_ps (im1, jn1, nl1, ptop, fields3d(:,:,:,ndp,2),  & 
                                                 fields2d(:,:,1,nps,2))
 
         energy_filename = trim(dynfiles(3)) 
!
         call compute_e_field_ (im1,jn1,nl1,nfields3d,nall,   &
                   jweights(:,2),                             &
                   fields3d(:,:,:,:,1),fields3d(:,:,:,ndp,2), &
                   fields2d(:,:,1,nps,1),fnames(1:nfields3d), &
                   energy, energy_field)

         ncount = ncount + 1
         energy_field_ave = energy_field_ave + energy_field
         call tick ( nymd, nhms, tinc )

      end do
      print *,'Total number of perturbations processed = ', ncount

      energy_field_ave = energy_field_ave / ncount

      if ( notag ) then
          filename = trim(energy_filename)
      else
          write(filename, '(2a,i8,a,i2.2,a)') trim(energy_filename),'.', &
                                              nymd_beg,'_',nhms_beg/10000, 'z.nc4'
      endif
      call GFIO_PutFld ( 'gmao_pert_energy',  'GMAO',  filename,  &
                          nymd_beg,   nhms_beg, 240000, &
                          im1,jn1,nl1,ptop,ks,ak,bk, &
                          nall, vnames, vtitles, vunits, & 
                          energy_field_ave, ier, untag=.true. )

      deallocate(fields3d)
      deallocate(fields2d)
      deallocate(energy_field_ave)
      deallocate(energy_field)
      deallocate(energy)
      deallocate(bk)
      deallocate(ak)
      deallocate(delbk)
      deallocate(pmid)
      deallocate(glats)
      deallocate(kweights)
      deallocate(jweights)
 
      end subroutine energy_field_Process_
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: compute_e_field_ --- Compute contributions to the total energy norm E.
!
!
! !INTERFACE:
!
      subroutine compute_e_field_ (im1,jn1,nl1,nfields,nall,&
                   jweights,                                &
                   fields,delp_ref,psfield,fnames,          &
                   e, energy_field) 

!USES:

      use m_pertutil, only: pertutil_getparam
      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1,nfields,nall
      character(len=*), intent(in) :: fnames(nfields) 
      real(r8), intent(in)  :: jweights(jn1)
      real(r8), intent(in)  :: fields(im1,jn1,nl1,nfields)
      real(r8), intent(in)  :: delp_ref(im1,jn1,nl1)
      real(r8), intent(in)  :: psfield(im1,jn1)   

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: e(nl1+2,3)   
      real(r8), intent(out) :: energy_field(im1,jn1,nl1,nall)

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
!  04Jun2003  R. Errico  Initial algorithm
!  08Apr2005  R. Errico  Corrected assumption that sum(jweights)=1
!  10Jan2008  Todling    - Add wet(q) energy component  
!                        - Generalized to output all energy partitions
!  22Jan2012  Todling    Allow for vertical norm change
!  16Mar2012  Todling    Handle d/a grids
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      character(len=*),parameter:: myname_='compute_e_field_'

      integer  :: i,j,k
      integer  :: idu,idv,idt,idq
      real(r8) :: dbleimr
      real(r8) :: tfac, ufac, pfac, qfac
      real(r8) :: dsigma(im1,jn1,nl1)
      real(r8) :: ps_ref(im1,jn1)
      real(r8) :: jfac(jn1)      
      real(r8) :: fa(im1,jn1) 
      real(r8) :: fb(im1,jn1) 
      real(r8) :: value
      real(r8) :: eu,ev,et,ep,eq
      real(r8) :: euu,evv,ett,eqq

      character(len=1) :: wgrid
      integer    :: nlp1,ierr
      integer    :: ks
      real(r8)   :: ak(nl1+1)
      real(r8)   :: bk(nl1+1)
      real(r8)   :: pint
      real(r8)   :: ptop
      real(r8)   :: wfactor
      real(r8)   :: eps_eer
      real(r8), allocatable :: pe(:,:,:)
      integer vnorm
      integer kz(1)
        
      call set_eta(nl1, ks, ptop, pint, ak, bk)

      dbleimr=1.d0/dble(im1)
      do j=1,jn1
        jfac(j)=jweights(j)*dbleimr
      enddo

! One factor of 0.5 is for defenition of energy as 1/2 * ...
! Assumption here is that sum(jweights)=1
      call getparam ( 'vnorm', vnorm )
      call getparam ( 'eps_eer', eps_eer )
      print *, 'enorm using eps_eer: ', eps_eer
      tfac=0.5d0*Cp/tref
      ufac=0.5d0
      pfac=0.5d0*R*tref/pref**2
      qfac=0.5d0*eps_eer*alhl*alhl/(Cp*tref)

      call find_name (nfields,fnames,'u__',idu)      
      call find_name (nfields,fnames,'v__',idv)      
      call find_name (nfields,fnames,'T__',idt)      
      call find_name (nfields,fnames,'q__',idq)      

      call pertutil_getparam ( 'wgrid', wgrid )

! compute 3-D varying dsigma

      call calculate_ps (im1,jn1,nl1,ptop,delp_ref,ps_ref) 

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

         print *, 'Using E-norm'
         do k=1,nl1
           dsigma(:,:,k)=delp_ref(:,:,k)/ps_ref(:,:)
         enddo    

      endif
!
! Initialize e to zero
!
      e(:,:)=0.d0
      energy_field=0.d0
      euu=0.d0 
      evv=0.d0
      ett=0.d0 
      eqq=0.d0 
!
!  Loop over levels
!
      do k=1, nl1
!     
!  Calculations for u field
!     
        fa(:,:)=0.d0
        fb(:,:)=0.d0
        do j=1,jn1
          do i=1,im1
            fa(i,j)=fields(i,j,k,idu)   
          enddo
        enddo


        eu=0.d0
        if ( wgrid=='d' ) then
           do i=1,im1
             fb(i,  1)=fa(i,  2)*fa(i,  2)
             fb(i,jn1)=fa(i,jn1)*fa(i,jn1)
             do j=2,jn1-1     
               fb(i,j)=0.5*(fa(i,j)*fa(i,j) + fa(i,j+1)*fa(i,j+1))
             enddo
           enddo

           do j=1,jn1
             do i=1,im1
               value=dsigma(i,j,k)*jfac(j)*fb(i,j)
               eu=eu+value
               energy_field(i,j,k,nkxe)=energy_field(i,j,k,nkxe) + ufac*value
             enddo
           enddo
        else if ( wgrid=='a' ) then
           do j=1,jn1
             do i=1,im1
               value=dsigma(i,j,k)*jfac(j)*fa(i,j)*fa(i,j)
               eu=eu+value
               energy_field(i,j,k,nkxe)=energy_field(i,j,k,nkxe) + ufac*value
             enddo
           enddo
        else
           call die(myname_,'grid for winds not properly defined, aborting')
        endif
        euu=euu+eu
!     
!  Calculations for v field
!    
        fa(:,:)=0.d0
        fb(:,:)=0.d0
        do j=1,jn1
          do i=1,im1
            fa(i,j)=fields(i,j,k,idv)     
          enddo
        enddo
        ev=0.d0
        if ( wgrid=='d' ) then
           do i=1,im1
             fb(i,  1)=fa(i,    2)*fa(i,    2)
             fb(i,jn1)=fa(i,jn1-1)*fa(i,jn1-1)
           enddo
           do j=2,jn1-1
             fb(im1,j)=0.5*(fa(1,j)*fa(1,j) + fa(im1,j)*fa(im1,j))
             do i=1,im1-1     
               fb(i,j)=0.5*(fa(i,j)*fa(i,j) + fa(i+1,j)*fa(i+1,j)) 
             enddo
           enddo

           do j=1,jn1
             do i=1,im1
               value=dsigma(i,j,k)*jfac(j)*fb(i,j)   
               ev=ev+value
               energy_field(i,j,k,nkxe)=energy_field(i,j,k,nkxe) + ufac*value
             enddo
           enddo
        else if ( wgrid=='a') then
           do j=1,jn1
             do i=1,im1
               value=dsigma(i,j,k)*jfac(j)*fa(i,j)*fa(i,j)
               ev=ev+value
               energy_field(i,j,k,nkxe)=energy_field(i,j,k,nkxe) + ufac*value
             enddo
           enddo
        else ! will never make it here
           call die(myname_,'grid for winds not properly defined, aborting')
        endif
        evv=evv+ev

        e(k,1)=ufac*(eu+ev)   

        e(nl1+1,1)=e(nl1+1,1)+e(k,1)

!     
!  Calculations for T field
!     
        fa(:,:)=0.d0
        do j=1,jn1
          do i=1,im1
            fa(i,j)=fields(i,j,k,idt) 
          enddo
        enddo

        et=0.d0
        do j=1,jn1
          do i=1,im1
            value=dsigma(i,j,k)*jfac(j)*fa(i,j)*fa(i,j)
            et=et+value
            energy_field(i,j,k,nape)=energy_field(i,j,k,nape) + tfac*value
          enddo
        enddo
        ett=ett+et

        e(k,2)=tfac*et
        e(k,3)=e(k,1)+e(k,2)
        e(nl1+1,2)=e(nl1+1,2)+e(k,2)
!     
!  Calculations for Q field
!     
        fa(:,:)=0.d0
        do j=1,jn1
          do i=1,im1
            fa(i,j)=fields(i,j,k,idq) 
          enddo
        enddo

        eq=0.d0
        do j=1,jn1
          do i=1,im1
            value=dsigma(i,j,k)*jfac(j)*fa(i,j)*fa(i,j)
            eq=eq+value
            energy_field(i,j,k,nqxe)=energy_field(i,j,k,nqxe) + qfac*value
          enddo
        enddo
        eqq=eqq+eq

!
!
      enddo ! loop over k
!
! Calculations for ps field
!
      fa(:,:)=0.d0
      do j=1,jn1
        do i=1,im1
          fa(i,j)=psfield(i,j)
        enddo
      enddo

      ep=0.d0
      do j=1,jn1
        do i=1,im1
          value=jfac(j)*fa(i,j)*fa(i,j)
          ep=ep+value
          do k=1, nl1
             energy_field(i,j,k,ngpe)=energy_field(i,j,k,ngpe)+pfac*dsigma(i,j,k)*value
          enddo
        enddo
      enddo

      print*, 'Energy scaled fields (u,v,t,p,q): ', ufac*euu,ufac*evv,tfac*ett,pfac*ep,qfac*eqq
      energy_field(:,:,:,nape) = energy_field(:,:,:,nape) + energy_field(:,:,:,ngpe)
      energy_field(:,:,:,ntxe) = energy_field(:,:,:,nape) + energy_field(:,:,:,nkxe)
      energy_field(:,:,:,ntwe) = energy_field(:,:,:,nqxe) + energy_field(:,:,:,ntxe)

      e(nl1+1,3)=e(nl1+1,1)+e(nl1+1,2)
      e(nl1+2,1)=0.d0                            ! not used
      e(nl1+2,2)=pfac*ep
      e(nl1+2,3)=e(nl1+2,2)+e(nl1+1,3)           ! total E

      print *, 'Energy partition: '
      print *, 'Kinetic          Energy: ', sum(energy_field(:,:,:,nkxe))
      print *, 'Potential        Energy: ', sum(energy_field(:,:,:,ngpe))
      print *, 'Avail. Potential Energy: ', sum(energy_field(:,:,:,nape))
      print *, 'Total dry        Energy: ', sum(energy_field(:,:,:,ntxe))
      print *, 'Total wet        Energy: ', sum(energy_field(:,:,:,ntwe))

      end subroutine compute_e_field_  

!-------------------------------------------------------------------------
!         NASA/GSFC, GMAO  Code 610.1 GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pertengy_set_ --- General paramater setup for pertenergy calc
!
! !INTERFACE:
!
      subroutine pertengy_set_ ( comm, root, rcfile, stat )
                                                                                                 
      Implicit None
                                                                                                 
! !INPUT PARAMETERS:
!
      integer, intent(in) ::  comm      ! MPI communicator
      integer, intent(in) ::  root      ! MPI ROOT PE

      character(len=*), intent(in) :: rcfile  ! resource filename
                                                                                                 
! !OUTPUT PARAMETERS:

      integer, intent(out) :: stat       ! Return error code:
                                         !   stat=1 ...
                                         !   stat=2 ...
                                         !   stat=3 ...

! !FILES USED:  pert_energy.rc
!
! !DESCRIPTION:  Reads resource file, does basic setup for calculation of
!                 3-D perturbation energy field (for sensitivty vectors or SVs)
!
! !REVISION HISTORY:
!
!   26Apr2007  Gelaro      Initial code
!   05Jan2008  Todling     Revisited
!   22Jan2012  Todling     Add vnorm option
!
!EOP
!-------------------------------------------------------------------------
                                                                                                 
      character(len=*), parameter :: myname_ = myname//'*pertenergy_Set'

      character(len=255) token
      character(len=255) pertrc

      integer     iret, ierr
      integer     myID
      real(r8)    eps_eer
      integer     vnorm
   
      stat    = 0

!     Initialize parameters independently of whether RCfile used or not
!     -----------------------------------------------------------------
      eps_eer = 1.0d0
      call setparam ( 'eps_eer', eps_eer )

      if (rcfile=='null') return 

      call MP_comm_rank(comm,myID,ierr)
        if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)

!     Load resources from pert_energy.rc
!     ----------------------------------
      pertrc = ' '
      call getenv('FVPERT_RC',pertrc)           ! Unix binding
      if(pertrc.eq.' ') pertrc=rcfile           ! default name
      call i90_loadf (trim(pertrc), iret)
      if( iret .ne. 0) then
          if(myID==ROOT) then
            write(stdout,'( a  )') '----------------------------------'
            write(stdout,'(2a  )') myname_, ': No resource file found'
            write(stdout,'( a,/)') '----------------------------------'
            return
          end if
      end if
      if(myID==ROOT) then
         write(stdout,'( a  )') '---------------------------------'
         write(stdout,'(2a  )') myname_, ': Reading resource file'
         write(stdout,'( a,/)') '---------------------------------'
      end if
                                                                                                 
!     Read template for pert filename
!     -------------------------------
      call I90_label('pert_template:', iret)
      if (iret .eq. 0) then
        call I90_Gtoken ( token, iret )
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_, ': I90_Gtoken error, iret =',iret
          call die(myname_)
        else
          pertmpl= trim(token)
        end if
      end if
      if(myID==ROOT) write(stdout,'(2a)') 'Perturbation template: ', trim(pertmpl)
                                                                                                 
!     Read template for ref state filename
!     -----------------------------------
      call I90_label('ref_template:', iret)
      if (iret .eq. 0) then
        call I90_Gtoken ( token, iret )
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_, ': I90_Gtoken error, iret =',iret
          call die(myname_)
        else
          reftmpl= trim(token)
        end if
      end if
      if(myID==ROOT) write(stdout,'(2a)') 'Reference state template: ', trim(reftmpl)

!     Read time increment between cases
!     ---------------------------------
      call I90_label('time_increment:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,': I90_label error(time_increment), iret =',iret
         call die(myname_)
      end if
      tinc = I90_GInt(iret)
      if( iret .ne. 0) then
         write(stderr,'(2a,i5)') myname_,': I90_GInt error(tinc), iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) write(stdout,'(a,i4)') 'Time inc between perturbations (hrs): ', tinc
      tinc = tinc * 3600

!     Read forward forecast length
!     ----------------------------
      call I90_label('tau_forward:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,': I90_label error(tau_forward), iret =',iret
         call die(myname_)
      end if
      ifwd = I90_GInt(iret)
      if( iret .ne. 0) then
         write(stderr,'(2a,i5)') myname_,': I90_GInt error(ifwd), iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) write(stdout,'(a,i4)') 'Forecast integrations time (hrs): ', ifwd
      ifwd = ifwd * 3600

!     Read backward (adjoint) forecast length
!     ---------------------------------------
      call I90_label('tau_backward:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,': I90_label error(tau_backward), iret =',iret
         call die(myname_)
      end if
      ibak = I90_GInt(iret)
      if( iret .ne. 0) then
         write(stderr,'(2a,i5)') myname_,': I90_GInt error(ibak), iret =',iret
         call die(myname_)
      end if
      if(myID==ROOT) write(stdout,'(a,i4)') 'Adjoint integrations time (hrs): ', ibak
      ibak = ibak * 3600

!     Read Ehrendorfer, Errico, and Raeder's epsilon factor (apply to Q-component)
!     ----------------------------------------------------------------------------
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
      if(myID==ROOT) write(stdout,'(a,1p,e13.6)') 'Ehrendorfer, Errico, and Raeder eps: ',eps_eer
      call setparam ( 'eps_eer', eps_eer )

!     Decide between E- and V-norms
!     -----------------------------
      call I90_label('vnorm:', iret)
      if (iret .ne. 0) then
        write(stdout,'(2a)') myname_,    &
                     ': I90_label not found will use default '
      else
        vnorm = I90_GInt(iret)
        call I90_Gtoken ( token, iret )
        if (iret .ne. 0) then
            write(stderr,'(2a,i5)') myname_,    &
                         ': I90_Gtoken error, iret =',iret
            call die(myname_)
        end if
        if( vnorm > 0 ) then
            call setparam ( 'vnorm', vnorm )
            if(myID==ROOT) write(stdout,'(a,i5)')   &
                   'Using V-weights option: ', vnorm
         endif
      end if

!     release resource file:
!     ---------------------
      call I90_release()

      readRC_ = .true.
                                                                                                 
      return
      end subroutine pertengy_set_

!-------------------------------------------------------------------------
!         NASA/GSFC, GMAO  Code 610.1 GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fill_templ --- Advance date/time and fill in templates when 
!                            processing multiple files
!
! !INTERFACE:
!
      subroutine  fill_templ ( dynfiles, expid, nfiles, nymd, nhms )

      Implicit None

! !INPUT PARAMETERS:

      integer, intent(in) ::  nfiles
      character(len=*), intent(in) :: expid

! !OUTPUT/OUTPUT PARAMETERS:

      character(len=*), intent(inout) :: dynfiles(nfiles)

      integer, intent(inout) :: nymd
      integer, intent(inout) :: nhms

! !DESCRIPTION:
!
! Advance date/time and fill in templates when processing multiple files
!
! !REVISION HISTORY:
!
!  02May2007  R. Gelaro    Initial algorithm
!  10Jan2008  Todling      Revisited time ticking
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer nymdb, nhmsb, nymde, nhmse

      character(len=255)  str, tmp
      character(len=255)  pertfile, refsfile

! Compute beginning and ending times in filenames
! -----------------------------------------------

      nymde = nymd
      nhmse = nhms

      call tick ( nymde, nhmse, ibak )

      nymdb = nymde
      nhmsb = nhmse

      call tick ( nymdb, nhmsb, -ifwd )

! Fill template for perturbation file
! -----------------------------------

! initial date of forward integration
      call strTemplate ( str, pertmpl, nymd=nymdb, nhms=nhmsb, xid=expid )
! ending date of adjoint integration
      call strTemplate ( tmp, str, nymd=nymde, nhms=nhmse )
! beginning date of adjoint integration
      call strTemplate ( pertfile, tmp, nymd=nymd, nhms=nhms )


! Fill template for ref state file
! --------------------------------

! initial date of forward integration
      call strTemplate ( str, reftmpl, nymd=nymdb, nhms=nhmsb, xid=expid )
! ending date of adj integration
      call strTemplate ( refsfile, str, nymd=nymd, nhms=nhms )
                                                                                              
      print *, 'Perturbation file = ', trim(pertfile)
      print *, 'Reference    file = ', trim(refsfile)

      if(nfiles>=1) dynfiles(1) = pertfile
      if(nfiles>=2) dynfiles(2) = refsfile

      return
      end subroutine fill_templ

      subroutine StatsClean_ ()
      end subroutine StatsClean_
             
      end module m_pertenergy

