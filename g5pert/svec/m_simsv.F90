!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_simsv ---
!
! !INTERFACE:
                                                                                                       
      module m_simsv

      use precision
                                                                                                       
      use prognostics, only : prognostics_initial
      use prognostics, only : prognostics_final
      use prognostics, only : dyn_prog

      use m_const, only : Cp     => cpm   
      use m_const, only : R      => rgas  
      use m_const, only : kappa
      use m_const, only : pstd
      use m_const, only : tstd
   
      use m_pertutil, only : ptv2pt => pertutil_ptv2pt
      use m_pertutil, only : pt2t => pertutil_pt2t
      use m_pertutil, only : pt2t_tl => pertutil_pt2t_tl
      use m_pertutil, only : find_name => pertutil_find_name
      use m_pertutil, only : calculate_ps => pertutil_calculate_ps
      use m_pertutil, only : pkappa => pertutil_pkappa
      use m_pertutil, only : horiz_grid => pertutil_horiz_grid
      use m_pertutil, only : vert_grid => pertutil_vert_grid
      use m_pertutil, only : transform_t_tl => pertutil_transform_t_tl

      use m_pertutil, only : e_cross_product
      use m_pertutil, only : comp_kdomain
                                                                                                       
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type
                                                                                                       
      use m_mapz, only : set_eta
                                                                                                       
      use m_svprojs, only : proj_init
      use m_svprojs, only : proj_clean
      use m_svprojs, only : proj_svec
                                                                                                       
      use m_dyn

      use m_stdio
      use m_inpak90
      use m_die, only: MP_die, die
                                                                                                       
      implicit NONE
                                                                                                       
! !PUBLIC MEMBER FUNCTIONS:
                                                                                                       
      PRIVATE
                                                                                                       
      PUBLIC simsv_Set
      PUBLIC simsv_Process
      PUBLIC simsv_Clean
                                                                                                       
      interface simsv_Set ; module procedure SimSet_   ; end interface
      interface simsv_Process; module procedure simSV; end interface
      interface simsv_Clean; module procedure SimClean_; end interface
 
!
! !DESCRIPTION: Set up for adjoint initialization
!
! !REVISION HISTORY:
!
!  18Oct2004  Winslow   Initial code based on similar code
!  11Apr2005  Elena N.  Linked constants of this module to constants in /shared/hermes/
!  15Jul2005  Elena N.  Corrected bug in routine calls, where fields/fields3d were messed up
!                       This bug resulted to incorrectly accessed fileds in a dynamic vector   
!  18Jul2005  Elena N.  Corrected BUG in e_cross_product calculation: psfields were incorrect 
!  16May2007  Todling   Introduced dyn_prog
!
!EOP
!-------------------------------------------------------------------------
 
      character(len=*), parameter :: myname = 'm_postStats'
      real(r8), parameter :: pref = 100.d0*pstd
      real(r8), parameter :: tref = tstd
 
      CONTAINS
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: simsv - compute similarity of two sets of singular vectors 
!
!
! !INTERFACE:
!
! NOTE to Ron E.: In SV code, 
!   check calc of E for v at NP:  should be cosp(j) in SV code
!   check  pe(:,:,k) = pe(:,:,k-1) + delp(:,:,k)
!
!
      subroutine simSV (comm,ROOT, &
                        ntimes,times,nymd, details, nvecsA, nvecsB,jobR,jobA,jobB,RCfile)

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: comm
      integer, intent(in) :: ROOT
      integer, intent(in) :: ntimes         ! number of model times to be exam.
      integer, intent(in) :: times(ntimes)  ! times to be examined
      integer, intent(in) :: nymd           ! date to be examined
      logical, intent(in) :: details
      integer, intent(in) :: nvecsA      
      integer, intent(in) :: nvecsB      
      character(len=*), intent(in) :: jobR
      character(len=*), intent(in) :: jobA(nvecsA)
      character(len=*), intent(in) :: jobB(nvecsB)
      character(len=*), intent(in) :: RCfile

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:
       
! !DESCRIPTION:
! compute similarity of two sets of singular vectors

! !SEE ALSO:
!
       
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!  04Nov2004  Winslow    Placed into modular form
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
      integer,  dimension(2), parameter :: kdomainps=(/1,1/) ! Do not change

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
      integer, parameter :: nfields3d=5   ! number of 3D fields to be examined
      integer, parameter :: nfields2d=1   ! number of 2D fields to be examined
      integer, parameter :: nfields=nfields3d+nfields2d

      character(len=3), dimension(nfields) :: fnames
      character(len=3), dimension(nfields), parameter :: fnames_def = &
                  (/ 'u__','v__','T__','q__','dp_','ps_' /)
!
! Stuff for LPO
      type(dyn_prog)        :: ypert           ! perturbation array corres. to x
      real(r8), allocatable :: domain(:,:,:)   ! mask from local projection operator
      real(r8), allocatable :: udomain(:,:,:)  ! mask from local projection operator
      real(r8), allocatable :: vdomain(:,:,:)  ! mask from local projection operator
      real(r8), allocatable :: psdomain(:,:,:)  ! mask from local projection operator
      integer  :: m                          ! size of long vector
      integer ::  ndim            ! size of state vector
      integer ::  dims(5)         ! array of grid parameters
      integer :: nev              ! dummy argument (number of eigenvectors)
      integer :: ncv              ! dummy argument (number of Lanczos vectors)
      integer :: projrc           ! return flag for projection operator
      integer :: ibeg_u, iend_u
      integer :: ibeg_v, iend_v
      integer :: ibeg_delp, iend_delp
      integer :: ibeg_pt, iend_pt
      integer :: ibeg_q, iend_q
      integer :: i,j,k,nf         ! looping indices

! Counters for do-loops and for other array indecies
      integer :: kind
      integer :: ndT, ndp, nps, kdp  ! indicates particular fields
      integer :: nvecA
      integer :: nvecB
      integer :: nvecsC

      real(r8), allocatable  :: fields3dA(:,:,:,:,:)
      real(r8), allocatable  :: fields3dB(:,:,:,:,:)
      real(r8), allocatable  :: fields3dC(:,:,:,:)
      real(r8), allocatable  :: psfieldA(:,:)
      real(r8), allocatable  :: psfieldB(:,:)
      real(r8), allocatable  :: psfieldC(:,:)
      real(r8), allocatable  :: energyA(:,:)
      real(r8), allocatable  :: energyB(:,:)
      real(r8), allocatable  :: energyAB(:,:)
      real(r8), allocatable  :: simAB(:,:)
      real(r8), allocatable  :: simC(:)

      integer  :: kdomain(2)                   ! indexes corresp to pdomain
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

! 
!                      BEGIN EXECUTION
!

      fnames = fnames_def
      call find_name (nfields3d,fnames,'dp_',ndp)


! Set some grid information
      call read_dim_   (im1,jn1,nl1,lm1,jobR)
      ibeg_u  = 1
      iend_u  = ibeg_u  + im1 * jn1 * nl1- 1
      ibeg_v  = iend_u  + 1
      iend_v  = ibeg_v  + im1 * jn1 * nl1- 1
      ibeg_pt = iend_v  + 1
      iend_pt = ibeg_pt + im1 * jn1 * nl1- 1
      ibeg_q  = iend_pt + 1
      iend_q  = ibeg_q  + im1 * jn1 * nl1- 1
      ibeg_delp = iend_q    + 1
      iend_delp = ibeg_delp + im1 * jn1 * nl1- 1

!
!  construct grid mask of local projection operator (LPO)
!
      projrc = 1
                                                       
      dims(1) = im1
      dims(2) = jn1
      dims(3) = nl1
      dims(4) = 1
      dims(5) = jn1
                                                                                     
      m = 5*im1*jn1*nl1
      call prognostics_initial (ypert,im1,jn1,nl1,1)
      ypert%u = 1.0
      ypert%v = 1.0
      ypert%pt = 1.0

      call proj_init (  comm, root, stat=projrc, dims=dims, RCfile=RCfile )
      if(projrc/=0) print*,'projrc = ',projrc,' proj_init'
                                                                                      
      call proj_svec ( comm, ROOT, ypert, projrc )   ! Apply local projection operator
      if(projrc/=0) print*,'projrc = ',projrc,' proj_svec'
                 
      call proj_clean( projrc )
      if(projrc/=0) print*,'projrc = ',projrc,' proj_clean'
                 
      allocate (udomain   (im1,jn1,nl1))
      udomain(:,:,:)    =  ypert%u          ! <-- looks fishy _RT
      allocate (vdomain   (im1,jn1,nl1))
      vdomain(:,:,:)    =  ypert%v          ! <-- looks fishy _RT
      allocate (domain    (im1,jn1,nl1))
      domain(:,:,:)   =  ypert%pt           ! <-- looks fishy _RT
      allocate (psdomain    (im1,jn1,nl1))
      psdomain(:,:,:)   =  ypert%pt         ! <-- looks fishy _RT

      call prognostics_final (ypert)

      allocate (fields3dA(im1,jn1,nl1,nfields,2))
      allocate (fields3dB(im1,jn1,nl1,nfields,2))
      allocate (fields3dC(im1,jn1,nl1,nfields))
      allocate (psfieldA(im1,jn1))
      allocate (psfieldB(im1,jn1))
      allocate (psfieldC(im1,jn1))
      allocate (energyA(nl1+2,3))
      allocate (energyB(nl1+2,3))
      allocate (energyAB(nl1+2,3))
      allocate (simAB(1+nvecsA,1+nvecsB))
      allocate (jweights(jn1,2)) ! area weights (1) u, (2) T
      allocate (kweights(nl1))   ! vertical layer weights
      allocate (glats(jn1,2))    ! degrees lats (1) u, (2) T
      allocate (pmid(nl1))       ! p at data levels for standard grid
      allocate (delbk(nl1))      ! "thickness" of ps dependent portion of
      allocate (ak(nl1+1))    
      allocate (bk(nl1+1))    

!   Compute kdomain to maintain as much code history as possible
      call comp_kdomain (im1, jn1, nl1, domain, kdomain)
      k = kdomain(1)
      psdomain(:,:,1)   =  domain(:,:,k)

!    Snag some grid information
      call horiz_grid (im1, jn1, jweights, glats)
      call vert_grid  (nl1, ks, ak, bk, kweights, ptop, pmid, delbk)

!
! Read reference state into array fields3dA(:,:,:,:,2) and copy into B
      call read_data2_ (im1,jn1,nl1,lm1,nfields,fnames, &
                        jobR,times(1),nymd, &
                        fields3dA(:,:,:,:,2))

      fields3dB(:,:,:,:,2)=fields3dA(:,:,:,:,2)

!      
! Loop over set of vectors A
!
      do nvecA=1,nvecsA

!
! Read vector number nvecA from set A into fields3dA(:,:,:,:,1)  
        call read_data2_ (im1,jn1,nl1,lm1,nfields,fnames, &
                          jobA(nvecA),times(1),nymd,      &
                          fields3dA(:,:,:,:,1))

!
! Transform T fields 
        call transform_T_tl (im1,jn1,nl1,ptop,nfields,fnames, &
                                  fields3dA(:,:,:,:,1:2) )
!
! Compute ps fields
        call calculate_ps (im1,jn1,nl1,ptop,fields3dA(:,:,:,ndp,1),psfieldA(:,:))
!
! Compute E-norm for A as normalization factor 
!
        fields3dC(:,:,:,:)=fields3dA(:,:,:,:,1)
!       psfieldC(:,:,:)  =psfieldA(:,:,:)
        psfieldC(:,:)  =psfieldA(:,:)
        print*
        print*,'nvecA = ',nvecA
        print*
        call e_cross_product (im1,jn1,nl1,nfields,ptop,      &
           udomain, vdomain, domain, kdomain, reshape(jweights(:,2),(/jn1/)),   &
           reshape(fields3dA(:,:,:,:,1),(/im1,jn1,nl1,nfields/)),fields3dC(:,:,:,:),  &
           psfieldA(:,:),psfieldC(:,:), &
           reshape(fields3dA(:,:,:,ndp,2),(/im1,jn1,nl1/)), &
           fnames(1:nfields), energyA(:,:) )

!
! DEBUG
!       print ('(a,i3,3x,e20)'),'Energy for A-vector ', &
!              nvecA,energyA(nl1+2,3)
!
! Print table of contributions to E-norm if desired
!        if (details) then
! FIX QQQQQQ   call print_e_ (nl1,1,1,knames,             &
!                       lldomain,pdomain,times,pmid,energy)
!        endif
!      
! Loop over set of vectors B
!
        do nvecB=1,nvecsB

!
! Read vector number nvecB from set B into fields3dB(:,:,:,:,1)  
        call read_data2_ (im1,jn1,nl1,lm1,nfields,fnames, &
                          jobB(nvecB),times(1),nymd,      &
                          fields3dB(:,:,:,:,1))

!
! Transform T fields 
          call transform_T_tl (im1,jn1,nl1,ptop,nfields,fnames, &
                                  fields3dB(:,:,:,:,1:2) )
!
! Compute ps fields
         call calculate_ps (im1,jn1,nl1,ptop,fields3dB(:,:,:,ndp,1),  &
                             psfieldB(:,:))
!
! Compute E-norm for B as normalization factor 
!
          fields3dC(:,:,:,:)=fields3dB(:,:,:,:,1)
          psfieldC(:,:)  =psfieldB(:,:)
          call e_cross_product (im1,jn1,nl1,nfields,ptop,       &
             udomain, vdomain, domain,  kdomain, reshape(jweights(:,2),(/jn1/)),   &
             reshape(fields3dB(:,:,:,:,1),(/im1,jn1,nl1,nfields/)),fields3dC(:,:,:,:), &
             psfieldB(:,:),psfieldC(:,:), &
             reshape(fields3dB(:,:,:,ndp,2),(/im1,jn1,nl1/)),  &
             fnames(1:nfields), energyB(:,:) )

!
! Compute cross product A*B as E-norm 
!
          call e_cross_product (im1,jn1,nl1,nfields,ptop,      &
             udomain, vdomain, domain, kdomain, reshape(jweights(:,2),(/jn1/)),   &
             reshape(fields3dA(:,:,:,:,1),(/im1,jn1,nl1,nfields/)),fields3dB(:,:,:,:,1),          &
             psfieldA(:,:),psfieldB(:,:),&
             reshape(fields3dB(:,:,:,ndp,2),(/im1,jn1,nl1/)), &
             fnames(1:nfields), energyAB(:,:) )
             simAB(nvecA,nvecB)=energyAB(nl1+2,3) /              &
                     sqrt(energyA(nl1+2,3)*energyB(nl1+2,3))

!
!
!          if (nvecA==1) then
!           print ('(a,i3,3x,e20)'),'Energy for B-vector ', &
!              nvecB,energyB(nl1+2,3)
!          endif 
        enddo   ! loop over B vectors
      enddo     ! loop over A vectors

!
! Replace elements of simAB by their squared values, sum rows and 
! columns, and compute similarity indices for various subsets
      nvecsC=min(nvecsA,nvecsB)
      allocate (simC(nvecsC))
      call sim_index_ (nvecsA,nvecsB,nvecsC,simAB,simC)
!
! Print matrix of squared projections 
      call print_simAB_ (nvecsA,nvecsB,simAB)
!
! Print similarity matrix 
      call print_simC_ (nvecsC,simC)

      deallocate (simAB,simC)
      deallocate (fields3dA,fields3dB,fields3dC)
      deallocate (psfieldA,psfieldB,psfieldC)
      deallocate (energyA,energyB,energyAB)
      deallocate(jweights,kweights,glats,pmid,delbk,ak,bk)
      deallocate(domain,udomain,vdomain,psdomain)

      print *,' '
      print *,' '
      print *,' '
      print *,' PROGRAM COMPLETED'
!

      end subroutine simsv

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
             
      call dyn_null ( w )
             
      end subroutine read_data2_
             
!
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: read_dim --- get dimensions for input vector to be analyzed
!
!
! !INTERFACE:
!
      subroutine read_dim_ (im1,jn1,nl1,lm1,job1)
                                                                                                         
!USES:
                                                                                                         
      implicit none
                                                                                                         
! !INPUT PARAMETERS:
                                                                                                         
      character(len=*), intent(in) :: job1
                                                                                                         
! !OUTPUT PARAMETERS:
                                                                                                         
      integer,  intent(out)  :: im1, jn1, nl1, lm1
                                                                                                         
! !INPUT/OUTPUT PARAMETERS:
                                                                                                         
! !DESCRIPTION:
                                                                                                         
! !SEE ALSO:
!
                                                                                                         
! !REVISION HISTORY:
!
!  04Oct2004  Winslow    Original code
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables
       
       
      type(dyn_vect) w
      integer :: ier, nymd_out, nhms_out
       
       
      call dyn_null ( w )
      call dyn_get ( trim(job1), nymd_out, nhms_out, w, ier )
      if(ier.ne.0) then
        print*,'trouble getting ',trim(job1)
      endif
      im1 = w%grid%im
      jn1 = w%grid%jm
      nl1 = w%grid%km
      lm1 = w%grid%lm
      call dyn_null ( w )
      end subroutine read_dim_
       
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: print_e --- print table of contributions by total energy norm
!
!
! !INTERFACE:
!
      subroutine print_e_ (nl,nkinds,ntimes,knames,            &
                          lldomain,pdomain,times,pmid,energy)


!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: nl,nkinds,ntimes
      character(len=*), intent(in) :: knames(nkinds)
      real(r8), intent(in)  :: lldomain(4)
      real(r8), intent(in)  :: pdomain(2)
      real(r8), intent(in)  :: times(ntimes)
      real(r8), intent(in)  :: pmid(nl)
      real(r8), intent(in)  :: energy(nl+2,3,nkinds,ntimes)
                                                                                                                     
! !OUTPUT PARAMETERS:
                                                                                         
! !INPUT/OUTPUT PARAMETERS:
                                                                                                                     
! !DESCRIPTION:

!
!  print table of contributions by total energy norm as a function of
!  data kind, data field, data time, and model pressure layer 
!

! !SEE ALSO:
!
                                                                                                                     
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: k, kind, nt

      print *,' '
      print *,' '
      print *,' '
      print ('(30x,a)'),'E - N O R M    S T A T I S T I C S'
      print ('(4(a,f6.1),2(a,f9.1))')                           &
               ,'Longitudes: ',lldomain(1),' to ',lldomain(2),  &
               ', Latitudes:  ',lldomain(3),' to ',lldomain(4), &
               ', pressures:  ',pdomain(1),' to ',pdomain(2) 
        
      do nt=1,ntimes
        do kind=1,nkinds

! Energy has only been computed for particular kinds
          if ( (knames(kind) (1:3) == 'TLM') .or.        &
               (knames(kind) (1:3) == 'DNN') .or.        &
               (knames(kind) (1:3) == 'DTT') .or.        &
               (knames(kind) (1:3) == 'DTN')      ) then
            print *,' '
            print *,' KIND=',knames(kind),'  TIME=',times(nt) 
            print ('(a,5x,a,8x,a,8x,a)'),'  k p-level','  KE',' APE','TOTE'  
            do k=1,nl
              print ('(i3,f8.1,1p3e12.2)'),k,pmid(k),energy(k,1:3,kind,nt)
            enddo
            print ('(a,1p3e12.2)'),'      sums:',energy(nl+1,1:3,kind,nt)
            print ('(a,1p3e12.2)'),'0, psE, TOT',energy(nl+2,1:3,kind,nt)
          endif

        enddo
      enddo

      end subroutine print_e_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sim_index --- Compute similarity index 
!
!
! !INTERFACE:
!
      subroutine sim_index_ (nvecsA,nvecsB,nvecsC,simAB,simC)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)    :: nvecsA, nvecsB, nvecsC

! !OUTPUT PARAMETERS:

      real,    intent(out)   :: simC(nvecsC)

! !INPUT/OUTPUT PARAMETERS:

      real,    intent(inout) :: simAB(1+nvecsA,1+nvecsB)

! !DESCRIPTION:
! 
!  Compute similarity index by summing rows and columns of the
!  simiarity matrix.  
!
!

! !SEE ALSO:
!
                                                                                                                     
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!  20Apr2005  R. Errico/Elena N.  Corrected bug in similarity index computation
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables
      integer :: nvecA,nvecB,nvecC
      real    :: sum
!
! Replace matrix of projection values by matrix of squares of
! projection values, and sum results over each row and column.
!
      simAB(1+nvecsA,:)=0.
      simAB(:,1+nvecsB)=0.
      do nvecA=1,nvecsA
        do nvecB=1,nvecsB
          simAB(nvecA,nvecB)=simAB(nvecA,nvecB)*simAB(nvecA,nvecB)
          simAB(1+nvecsA,nvecB)=simAB(1+nvecsA,nvecB)+simAB(nvecA,nvecB)
          simAB(nvecA,1+nvecsB)=simAB(nvecA,1+nvecsB)+simAB(nvecA,nvecB)
	enddo
      enddo
!
! Compute similarity index 
!
     do nvecC=1,nvecsC
       sum=0.
       do nvecA=1,nvecC
         do nvecB=1,nvecC
           sum=sum+simAB(nvecA,nvecB)
         enddo
       enddo
       simC(nvecC)=sum/real(nvecC)
     enddo
!
    end subroutine sim_index_ 

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: print_simAB --- Print matrix of squared projection values
!
!
! !INTERFACE:
!
      subroutine print_simAB_ (nvecsA,nvecsB,simAB)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: nvecsA, nvecsB
      real,    intent(in) :: simAB(1+nvecsA,1+nvecsB)

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:
                                                                                                                     
! !DESCRIPTION:
                                                                                                                     
! !SEE ALSO:
!
                                                                                                                     
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables
  
!
      integer :: ifac=100
      integer :: n1,n2
      integer :: isim(nvecsB+1)
!
      print *,' '
      print *,' '
      print *,' Matrix of squared projection values and their partial sums:'
      print *,' (values multiplied by', ifac,' for printing purposes)'
      print *,' Row index refers to vectors from set A'
      print *,' Column index refers to vectors from set B '
      print *,' Last column is the sum over all B-vectors for each A-vector'
      print *,' Last row is the sum over all A-vectors for each B-vector'
      print *,' '
!
      do n2=1,nvecsB
        isim(n2)=n2
      enddo
      print ('(a,25i4)'),'    ',isim(1:nvecsB)
      print *,' '
!
      do n1=1,nvecsA
        do n2=1,nvecsB+1
          isim(n2)=int(ifac*simAB(n1,n2))
        enddo
        print ('(i3,1x,(25i4))'),n1,isim(1:nvecsB+1)
      enddo 
!
      do n2=1,nvecsB
        isim(n2)=int(ifac*simAB(nvecsA+1,n2))
      enddo
      print *,' '
      print ('(4x,(25i4))'),isim(1:nvecsB)
!
     end subroutine print_simAB_ 

!
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: print_simC --- Print similarity index
!
!
! !INTERFACE:
!
      subroutine print_simC_ (nvecsC,simC)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: nvecsC
      real,    intent(in) :: simC(nvecsC)
  
! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:
                                                                                                                     
! !DESCRIPTION:
                                                                                                                     
! !SEE ALSO:
!
                                                                                                                     
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

!
      integer :: n
!
      print *,' '
      print *,' '
      print *,' Similarity index for spaces of leading n-vectors'
      print *,' '

      do n=1,nvecsC
        print ('(a,i4,a,f10.5)'),'n=',n,'  sim=',simC(n) 
      enddo 
!
     end subroutine print_simC_ 

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: SimSet_ --- General paramater setup for singular vector calc
!
! !INTERFACE:
!
      subroutine SimSet_ ( comm, root, ntimes, nhms, nymd, details, &
                               nvecsA, nvecsB, jobR, jobA, jobB,        &
                               RCfile, stat )
                                                                                                       
! !USES:
                                                                                                       
      Implicit None
                                                                                                       
! !INPUT PARAMETERS:
!
      integer, intent(in) ::  comm      ! MPI communicator
      integer, intent(in) ::  root      ! MPI ROOT PE
      integer, intent(in) ::  ntimes
      integer, intent(in) ::  nvecsA
      integer, intent(in) ::  nvecsB
      character(len=*), intent(in) :: RCfile  ! resource filename
                                                                                                       
! !OUTPUT PARAMETERS:
                                                                                                       
      logical, intent(out) :: details
      character(len=*), intent(out) :: jobR
      character(len=*), intent(out) :: jobA(nvecsA)
      character(len=*), intent(out) :: jobB(nvecsB)
      integer, intent(out) :: nymd            ! date (YYYYMMDD) of output perturbation
      integer, intent(out) :: nhms(ntimes)    ! time   (HHMMSS) of output perturbation
      integer, intent(out) :: stat            ! Return error code:
                                              !  stat=1 cannot determined proc
                                              !  stat=2 RC file not found
                                              !  stat=3 missing RC name indicator
                                              !  stat=4 missing RC entry
                                                                                                       
                                                                                                       
! !FILES USED:  fvsimsv.rc
!
! !DESCRIPTION:  Initializes parameters for similarity calculations
!
! !REVISION HISTORY:
!
!   30Nov2004  Winslow    Initial code, based on similar codes
!
!EOP
!-------------------------------------------------------------------------
      character(len=*), parameter :: myname_ = myname//'SimSet_'
                                                                                                       
      character(len=255) token
      character(len=255) simrc
 
      integer       i, j, k, iret, ierr
      integer       ifiles
      integer       myID
 
      stat = 0
 
      call MP_comm_rank(comm,myID,ierr)
        if (ierr/=0) then
            stat = 1  ! cannot grab proc id
            return
        endif
 
!     Load resources from fvsimsv.rc
!     ------------------------------
      simrc = ' '
      call getenv('FVPERT_RC',simrc)           ! Unix binding
      if(simrc.eq.' ') simrc=RCfile           ! default name
      call i90_loadf (trim(simrc), iret)
      if( iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_loadf error, iret =',iret
          stat = 2 ! missing RC file
          return
      end if
      if(myID==ROOT) then
         write(stdout,'( a  )') '---------------------------------'
         write(stdout,'(2a  )') myname_, ': Reading resource file'
         write(stdout,'( a,/)') '---------------------------------'
      end if
 
 
!     Define Reference file
!     ---------------------
         call I90_label('file_nameR:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobR = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for reference file: ', trim(jobR)
             
!     Define files for group A
!     ------------------------
      if ( nvecsA .ge. 1 ) then
         call I90_label('file_nameA1:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA(1) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file A group: ', trim(jobA(1))
      endif

      if ( nvecsA .ge. 2 ) then
         call I90_label('file_nameA2:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA(2) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file A group: ', trim(jobA(2))
      endif
                                                                                                                
      if ( nvecsA .ge. 3 ) then
         call I90_label('file_nameA3:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA(3) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file A group: ', trim(jobA(3))
      endif
                                                                                                                
      if ( nvecsA .ge. 4 ) then
         call I90_label('file_nameA4:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA(4) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file A group: ', trim(jobA(4))
      endif
                                                                                                                
      if ( nvecsA .ge. 5 ) then
         call I90_label('file_nameA5:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA(5) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file A group: ', trim(jobA(5))
      endif
                                                                                                                
      if ( nvecsA .ge. 6 ) then
         call I90_label('file_nameA6:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA(6) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file A group: ', trim(jobA(6))
      endif
                                                                                                                
      if ( nvecsA .ge. 7 ) then
         call I90_label('file_nameA7:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA(7) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file A group: ', trim(jobA(7))
      endif
                                                                                                                
      if ( nvecsA .ge. 8 ) then
         call I90_label('file_nameA8:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA(8) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file A group: ', trim(jobA(8))
      endif
                                                                                                                
      if ( nvecsA .ge. 9 ) then
         call I90_label('file_nameA9:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA(9) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file A group: ', trim(jobA(9))
      endif
                                                                                                                
!     Define files for group B
!     ------------------------
      if ( nvecsB .ge. 1 ) then
         call I90_label('file_nameB1:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobB(1) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(1))
      endif
                                                                                                                
      if ( nvecsB .ge. 2 ) then
         call I90_label('file_nameB2:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobB(2) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(2))
      endif
                                                                                                                
      if ( nvecsB .ge. 3 ) then
         call I90_label('file_nameB3:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobB(3) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(3))
      endif
                                                                                                                
      if ( nvecsB .ge. 4 ) then
         call I90_label('file_nameB4:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobB(4) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(4))
      endif
                                                                                                                
      if ( nvecsB .ge. 5 ) then
         call I90_label('file_nameB5:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobB(5) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(5))
      endif
                                                                                                                
      if ( nvecsB .ge. 6 ) then
         call I90_label('file_nameB6:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobB(6) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(6))
      endif
                                                                                                                
      if ( nvecsB .ge. 7 ) then
         call I90_label('file_nameB7:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobB(7) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(7))
      endif
                                                                                                                
      if ( nvecsB .ge. 8 ) then
         call I90_label('file_nameB8:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobB(8) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(8))
      endif
                                                                                                                
      if ( nvecsB .ge. 9 ) then
         call I90_label('file_nameB9:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobB(9) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(9))
      endif
                                                                                                                
!     Define output times
!     -------------------
      if ( ntimes .ge. 1 ) then
        call I90_label('time1:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_label error(time1), iret =',iret
          stat = 3 ! missing RC name indicator
          return
        end if
        nhms(1) = I90_GInt(iret)
        if( iret .ne. 0) then
           write(stderr,'(2a,i5)') myname_,': I90_GInt error(time1), iret =',iret
           stat = 4 ! missing RC entry
           return
        end if
      endif
             
      if ( ntimes .ge. 2 ) then
        call I90_label('time2:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_label error(time2), iret =',iret
          stat = 3 ! missing RC name indicator
          return
        end if
        nhms(2) = I90_GInt(iret)
        if( iret .ne. 0) then
           write(stderr,'(2a,i5)') myname_,': I90_GInt error(time2), iret =',iret
           stat = 4 ! missing RC entry
           return
        end if
      endif
             
      if ( ntimes .ge. 3 ) then
        call I90_label('time3:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_label error(time3), iret =',iret
          stat = 3 ! missing RC name indicator
          return
        end if
        nhms(3) = I90_GInt(iret)
        if( iret .ne. 0) then
           write(stderr,'(2a,i5)') myname_,': I90_GInt error(time3), iret =',iret
           stat = 4 ! missing RC entry
           return
        end if
      endif
             
      if ( ntimes .ge. 4 ) then
        call I90_label('time4:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_label error(time4), iret =',iret
          stat = 3 ! missing RC name indicator
          return
        end if
        nhms(4) = I90_GInt(iret)
        if( iret .ne. 0) then
           write(stderr,'(2a,i5)') myname_,': I90_GInt error(time4), iret =',iret
           stat = 4 ! missing RC entry
           return
        end if
      endif

!     Define output date
!     ------------------
        call I90_label('date:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_label error(date), iret =',iret
          stat = 3 ! missing RC name indicator
          return
        end if
        nymd = I90_GInt(iret)
        if( iret .ne. 0) then
           write(stderr,'(2a,i5)') myname_,': I90_GInt error(nymd), iret =',iret
           stat = 4 ! missing RC entry
           return
        end if

!     Decide if you want to print out the details of the similarity calculation
!     -------------------------------------------------------------------------
        call I90_label('details:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,    &
                       ': I90_label error, iret =',iret
          call die(myname_)
        end if
        call I90_Gtoken ( token, iret )
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,    &
                       ': I90_Gtoken error, iret =',iret
          call die(myname_)
        end if
        details = .false.
        if( trim(token) .eq. 'yes' ) details = .true.
        if(myID==ROOT) write(stdout,'(2a)')   &
                     'Print out the details: ', trim(token)
                                                                                                                     
!     release resource file:
!     ---------------------
      call I90_release()
 
      return
      end subroutine SimSet_
 
      subroutine SimClean_ ()
      end subroutine SimClean_
 
      end module m_simsv

