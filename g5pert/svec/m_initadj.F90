!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_initadj --- generates initial perturbation for ADM runs
!
! !INTERFACE:

      module m_initadj

      use precision

      use prognostics, only : prognostics_initial
      use prognostics, only : prognostics_final
      use prognostics, only : dyn_prog
      
      use m_const, only : Cp     => cpm
      use m_const, only : R      => rgas
      use m_const, only : kappa
      use m_const, only : pstd
      use m_const, only : tstd
      use m_const, only : alhl       ! latent heat condensation
      use m_const, only : radius_earth

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type

      use m_mapz, only : set_eta

      use m_svprojs, only : proj_init
      use m_svprojs, only : proj_clean
      use m_svprojs, only : proj_svec

      use m_pertutil, only : pt2vpt  => pertutil_pt2vpt
      use m_pertutil, only : t2pt    => pertutil_t2pt
      use m_pertutil, only : t2pt_ad => pertutil_t2pt_ad
      use m_pertutil, only : pt2t    => pertutil_pt2t
      use m_pertutil, only : pt2t_ad => pertutil_pt2t_ad
      use m_pertutil, only : pt2t_tl => pertutil_pt2t_tl
      use m_pertutil, only : ptv2pt  => pertutil_ptv2pt
      use m_pertutil, only : horiz_grid => pertutil_horiz_grid
      use m_pertutil, only : vert_grid => pertutil_vert_grid
      use m_pertutil, only : find_name => pertutil_find_name
      use m_pertutil, only : pkappa => pertutil_pkappa
      use m_pertutil, only : calculate_ps => pertutil_calculate_ps
      use m_pertutil, only : transform_q2r => pertutil_transform_q2r
      use m_pertutil, only : transform_q2r_ad => pertutil_transform_q2r_ad
      use m_pertutil, only : transform_T => pertutil_transform_T
      use m_pertutil, only : transform_t_tl => pertutil_transform_t_tl
      use m_pertutil, only : transform_T_adinv => pertutil_transform_T_adinv
      use m_pertutil, only : read_data => pertutil_read_data
      use m_pertutil, only : write_data => pertutil_write_data
      use m_pertutil, only : setparam => pertutil_setparam
      use m_pertutil, only : getparam => pertutil_getparam


      use m_pertutil, only : comp_kdomain

      use m_dyn

      use m_stdio
      use m_inpak90
      use m_chars,   only: lowercase
      use m_die, only: MP_die, die

      implicit NONE
 
! !PUBLIC MEMBER FUNCTIONS:
 
      PRIVATE

      PUBLIC InitAdj_Set
      PUBLIC InitAdj_iADM
      PUBLIC InitAdj_Clean

      interface InitAdj_Set ; module procedure PertSet_   ; end interface
      interface InitAdj_iADM; module procedure init_adjoint_; end interface
      interface InitAdj_Clean; module procedure clean_; end interface

!
! !DESCRIPTION: Set up for adjoint initialization
!
! !REVISION HISTORY:
!
!  24Oct2002  Todling   Modularized; split arpack-related initialization.
!  24Jun2002  Gelaro    Added nag-related initialization.
!  30Sep2002  Winslow   Hacked to work for adjinit
!  11Apr2005  Elena N.  Linked constants of this module 
!                       to already defined constants in shared/hermes 
!  10Jul2005  Elena N.  Changed a reference state to use a forecast vector,
!                       rather than an analysis vector 
!  15Aug2005  Elena N.  Output of Jnorm into a file 'Jnormf.txt' 
!  22Mar2006  Elena N.  Modified to work with GEOS-5
!  14Apr2006  Elena N.  Fixed bug with imposed conditions on v-wind in forecast error 
!  25Apr2006  Elena N.  Fixed LPO applied to DELP_AD
!  27Oct2007  Todling   Support to GEOS-5-type vectors
!  14Nov2007  Todling   Various spots to allow multiple tracers to pass through
!  13Dec2007  Todling   Add moisture to total energy norm
!  08Jan2008  Todling   Unification of post-proc routines
!  25Feb2008  Todling   - Correction for input perturbation handle
!                       - Add interpolation capability
!  27Oct2009  Todling   Pass vectype for GEOS-540 dyn-vect compatibility
!  15Oct2015  Holdaway  Added support for moist available enthalpy norm
!
!
! !REMARKS:
!
!  Relevant References:
!      
!   1) Errico, R. N., R. Gelaro, E. Novakovskaia, and R. Todling, 2007:
!         General characteristics of stratospheric singular vectors.
!         Meteorologische Zeitschrift, Vol. 16, 621-634.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'm_initadj'
      real(r8), parameter :: pref = 100.d0*pstd
      real(r8), parameter :: tref = tstd

      CONTAINS

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_adjoint_ --- Generate a perturbation/adjoint initialization
!
!
! !INTERFACE:
!
      subroutine init_adjoint_(comm,ROOT,&
                               lu,ofname,jtype,diff_file,init_field,nfiles,dynfiles,&
                               knames,nymd,nhms,idims,vectype,RCfile)
                                                                                                                           
!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: comm
      integer, intent(in) :: ROOT
      integer, intent(in) :: lu    ! logical unit for text output file
      character(len=*), intent(in)  :: ofname
      character(len=*), intent(in)  :: jtype
      logical, intent(in) :: diff_file
      character(len=*), intent(in)  :: init_field 
      integer, intent(in) :: nfiles
      character(len=*), intent(in) :: dynfiles(nfiles)
      character(len=*), intent(in) :: knames  (nfiles+1)
      integer, intent(inout) :: nymd 
      integer, intent(inout) :: nhms 
      integer, intent(in)    :: idims(2) 
      integer, intent(in)    :: vectype
      character(len=*), intent(in) :: RCfile

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION:
 
! !SEE ALSO:
!
 
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  22Mar2006  Elena N.   Modified to work with GEOS-5 
!  18Nov2007  Todling    Support for GEOS-5 vector-type
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: im2,im3   ! grid points in zonal direction from reference file(s)
      integer :: jn2,jn3   ! grid points in lat. direction from reference file(s)
      integer :: nl2,nl3   ! grid points in vert. direction from reference file(s)
      integer :: lm2,lm3   ! number of tracers from reference file(s)
      integer :: lm1

!  Any of the following options can be used for slots 2 and 3: RTL1, RTL2, RAD1, RAD2, NLM1, NLM2 

      integer, parameter :: nfields3d=5   ! number of 3D fields 
      integer, parameter :: nfields2d=1   ! number of 2D fields 
      integer, parameter :: nfields=nfields3d+nfields2d
      character(len=3), dimension(nfields) :: fnames
      character(len=3), dimension(nfields) :: fnames_def = &
                  (/ 'u__','v__','vpT','q__','dp_','ps_' /)
       
      real(r8), allocatable :: domain(:,:,:) ! mask from local projection operator
      real(r8), allocatable :: udomain(:,:,:) ! mask from local projection operator
      real(r8), allocatable :: vdomain(:,:,:) ! mask from local projection operator
      real(r8), allocatable :: ptdomain(:,:,:) ! mask from local projection operator
      real(r8), allocatable :: qdomain(:,:,:) ! mask from local projection operator
      real(r8), allocatable :: delpdomain(:,:,:) ! mask from local projection operator
      integer  :: kdomain(2)         ! indexes corresp to pdomain
      integer,  dimension(2), parameter :: kdomainps=(/1,1/)   ! Do not change
      real(r8), dimension(1), parameter :: kweightps=(/1.0d0/) ! Do not change

      integer  :: k_ref ! index for reference field used to define J
                        ! Generally, this is either reference field or diff 
      integer  :: ndT, ndp, nps, kdp  ! indicates particular fields
      integer  :: id_field        ! index for single field used to define J

      real(r8), allocatable :: fields2d1(:,:,:,:,:)
      real(r8), allocatable :: fields3d1(:,:,:,:,:)
      real(r8), allocatable :: aux4d (:,:,:,:) ! auxillary array
      real(r8), allocatable :: aux4d2(:,:,:,:) ! auxillary array
      real(r8), allocatable :: aux2d (:,:)     ! auxillary array

       
      real(r8), allocatable :: jweights(:,:) ! area weights (1) u, (2) T
      real(r8), allocatable :: kweights(:)   ! vertical layer weights
      real(r8), allocatable :: glats(:,:)    ! degrees lats (1) u, (2) T
      real(r8), allocatable :: pmid(:)       ! p at data levels for standard grid
      real(r8), allocatable :: delbk(:)      ! "thickness" of ps dependent portion of
                                             ! standard (CCM) hybrid-coord. vert. grid
      real(r8), allocatable :: ak(:)         ! hybrid coordinates
      real(r8), allocatable :: bk(:)         ! hybrid coordinates
      type(dyn_prog)        :: ypert         ! perturbation array corres. to x
      integer  :: m                          ! size of long vector
      real(r8) :: ptop                       ! pressure at top of grid
      real(r8) :: eps_eer
      integer  :: ks                         ! level of stratosphere

      logical  :: trans_T         ! flag that transforms between T and theta_v required
      logical  :: trans_ps        ! flag that transforms between ps and delp required
      logical  :: trans_q         ! flag that transforms between q and r = q/(1-q) (moist available enthalpy)
      integer :: ibeg_u, iend_u
      integer :: ibeg_v, iend_v
      integer :: ibeg_delp, iend_delp
      integer :: ibeg_pt, iend_pt
      integer :: ibeg_q, iend_q
      integer :: nu, nv           ! pointers for field data
      integer :: i,j,k,nf         ! looping indices
      integer :: ier

      integer ::  ndim            ! size of state vector
      integer ::  dims(5)         ! array of grid parameters
      integer :: nev              ! dummy argument (number of eigenvectors)
      integer :: ncv              ! dummy argument (number of Lanczos vectors)
      integer :: projrc           ! return flag for projection operator
      integer :: kinds            ! overall number of state-vect-like arrays
      real(r8) :: summ
      character(len=3) tname

      real(r8), allocatable :: Tza(:,:), rza(:,:)
      integer :: indexT, indexR

      fnames = fnames_def                          ! set default field names
      call getparam ( 'tname', tname )
      call find_name (nfields3d,fnames,'vpT',ndt)
      fnames(ndt) = tname                          ! reset temperature to what's in file
      print *, 'Temperature field on input set to: ', fnames(ndt)

      call find_name (nfields3d,fnames,'dp_',ndp)
      call find_name (nfields3d,fnames,'u__',nu )
      call find_name (nfields3d,fnames,'v__',nv )
      call find_name (nfields,  fnames,'ps_',nps)
      nps=nps-nfields3d      
      kinds = nfiles + 1

      call getparam ( 'eps_eer', eps_eer )

!  Read in reference field(s)
      call dyn_SetVecTyp ( ier, vectype=vectype )
      call dyn_getdim   (trim(dynfiles(1)),im2,jn2,nl2,lm2,ier)
      allocate(fields2d1(im2,jn2,  1,nfields,kinds))
      allocate(fields3d1(im2,jn2,nl2,nfields,kinds))
      allocate(aux4d(im2,jn2,nl2,nfields))

      fields3d1(:,:,:,:,:)=0.d0
      fields2d1(:,:,:,:,:)=0.d0

      call read_data (im2,jn2,nl2,nfields,fnames,nymd,nhms, &
                      trim(dynfiles(1)),fields3d1(:,:,:,:,2:2),'nlm',ier)
      if(ier/=0) call MP_die (myname,'Failed to read ref state',99)
!     fields2d1(:,:,1,nps,2) = fields3d1(:,:,1,nps,2)
      print*

! Set some grid information
      ibeg_u  = 1
      iend_u  = ibeg_u  + im2 * jn2 * nl2- 1
      ibeg_v  = iend_u  + 1
      iend_v  = ibeg_v  + im2 * jn2 * nl2- 1
      ibeg_pt = iend_v  + 1
      iend_pt = ibeg_pt + im2 * jn2 * nl2- 1
      ibeg_q  = iend_pt + 1
      iend_q  = ibeg_q  + im2 * jn2 * nl2- 1
      ibeg_delp = iend_q    + 1
      iend_delp = ibeg_delp + im2 * jn2 * nl2- 1

!  Read in a 2nd field if required, such as for computing forecast differences.
      if (nfiles>1) then
        call dyn_getdim   (trim(dynfiles(2)),im3,jn3,nl3,lm3,ier)
        if (nl3.ne.nl2.or.im3.ne.im2.or.jn3.ne.jn2) then
          print*,'interpolating for?',trim(dynfiles(2))
          print*,'interpolator runs externally'
          stop
        elseif (lm3.ne.lm2) then
          print*,' NUMBER OF TRACERS IN PROVIDED STATE VECTORS DIFFERS'
          print *, ' read_data lm2 = ', lm2
          print *, ' read_data lm3 = ', lm3
        endif
        call read_data (im3,jn3,nl3,nfields,fnames,nymd,nhms, &
                        trim(dynfiles(2)),fields3d1(:,:,:,:,3:3),'nlm',ier)
        if(ier/=0) call MP_die (myname,'Failed to read 2nd ref state',99)
        print*
      endif

      lm1 = lm2

!   Compute weighting terms
      allocate (jweights(jn2,2))
      allocate (kweights(nl2))
      allocate (glats(jn2,2))
      allocate (pmid(nl2))
      allocate (delbk(nl2))
      allocate (ak(nl2+1))
      allocate (bk(nl2+1))
      call horiz_grid (im2, jn2, jweights, glats)
      call vert_grid  (nl2, ks, ak, bk, kweights, ptop, pmid, delbk)

!
!  Determine flags to indicate that transforms between  T or ps and
!  theta_v or delp are required
      trans_T=.false.
      trans_ps=.false.
      trans_q=.false.
      if ( ( trim(jtype) == 'txe' )  .or. &
           ( trim(jtype) == 'twe' )  .or. &
           ( trim(jtype) == 'ave' )  .or. &
           ( trim(jtype) == 'mae' )  .or. &
           ( trim(jtype) == 'ate' )   ) then
        trans_T=.true.
      endif
      if ( ( trim(jtype) == 'txe' ) .or. &
           ( trim(jtype) == 'twe' ) .or. &
           ( trim(jtype) == 'mae' ) .or. &
           ( trim(jtype) == 'ape' ))  then
        trans_ps=.true.
      endif
      if ( (trim(jtype) == 'bmx' ) .or. &
           (trim(jtype) == 'bms' ) .or. &
           (trim(jtype) == 'rms' ) ) then
        if (trim(init_field) == 'T__' ) trans_T=.true.
        if (trim(init_field) == 'ps_' ) trans_ps=.true.
      endif
      if ( ( trim(jtype) == 'mae' ) ) then
        trans_q=.true.
      endif


!  Reference field has vpT in GEOS-4
!  Reference field has  vT in GEOS-5
!  Transform reference field(s) to dry T field if required
      if (trans_T) then   
        if ( knames(2)(1:3) == 'NLM' .and. knames(3)(1:3) == 'NLM' ) then
             call transform_T (im2,jn2,nl2,ptop,nfields,fnames,fields3d1(:,:,:,:,2))
             if (nfiles>1) then     
                fnames(ndT) = tname  ! reset to original
                call transform_T (im2,jn2,nl2,ptop,nfields,fnames,fields3d1(:,:,:,:,3))
             endif         
        else if ( knames(2)(1:3) == 'TLM' .and. knames(3)(1:3) == 'NLM' ) then
                     call transform_T_tl (im2,jn2,nl2,ptop,nfields,fnames,fields3d1(:,:,:,:,2:3) )
        else
             call MP_die (myname,'invalid choice of knames',99)
        endif         
      endif
      
!
!  Compute reference ps field(s) if required
      if (trans_ps) then   
        call calculate_ps (im2,jn2,nl2,ptop,fields3d1(:,:,:,ndp,2),fields2d1(:,:,1,nps,2))
        if (nfiles>1) then
          call calculate_ps (im2,jn2,nl2,ptop,fields3d1(:,:,:,ndp,3),fields2d1(:,:,1,nps,3))
        endif
      endif

!
!  Compute the water vapour mixing ratio from specific humidity
      if (trans_q) then
        if ( knames(2)(1:3) == 'NLM' .and. knames(3)(1:3) == 'NLM' ) then
             call transform_q2r (im2,jn2,nl2,nfields,fnames,fields3d1(:,:,:,:,2))
             if (nfiles>1) then
                call transform_q2r (im2,jn2,nl2,nfields,fnames,fields3d1(:,:,:,:,3))
             endif
        else
             call MP_die (myname,'moist enthalpy norm not compatible',99)
        endif     
      endif
      

!
!  Replace 2nd reference state by its difference with the 1st reference
!  state if required. Note this will be for transformed variables.
!  NOTE: aux4d(:,:,:,:) is necessary to keep FORECAST as the reference field
!  in case when 'diff_file = .true.' (DO NOT use analysis field as reference)
!
      if (diff_file) then
        aux4d(:,:,1,:) = fields2d1(:,:,1,:,2)     ! 4th index in fields2d is always 1
                                                  ! (see array allocate statement above)
        fields2d1(:,:,:,:,2)=fields2d1(:,:,:,:,2)-fields2d1(:,:,:,:,3)
        fields2d1(:,:,1,:,3) = aux4d(:,:,1,:)
                                                                                                                                               
        aux4d(:,:,:,:) = fields3d1(:,:,:,:,2)
        fields3d1(:,:,:,:,2)=fields3d1(:,:,:,:,2)-fields3d1(:,:,:,:,3)
        fields3d1(:,:,:,:,3) = aux4d(:,:,:,:)
      endif

!
!  construct grid mask of local projection operator (LPO)
!
      print*,'Constructing LPO'
      projrc = 1

      dims(1) = im2
      dims(2) = jn2
      dims(3) = nl2
      dims(4) = 1   ! jfirst      
      dims(5) = jn2 ! jlast

      m = 5*im2*jn2*nl2
      call prognostics_initial (ypert,im2,jn2,nl2,lm2)
      ypert%u  = 1.0
      ypert%v  = 1.0
      ypert%pt = 1.0
      ypert%q  = 1.0
      ypert%delp  = 1.0

      print*,'Building projector'
      call proj_init ( comm, root, stat=projrc, dims=dims,  RCfile=RCfile )
      if(projrc/=0) print*,'projrc = ',projrc,' proj_init'

      print*,'Projecting'
      call proj_svec ( comm, ROOT, ypert, projrc )   ! Apply local projection operator
      if(projrc/=0) print*,'projrc = ',projrc,' proj_svec'

      call proj_clean( projrc )
      if(projrc/=0) print*,'projrc = ',projrc,' proj_clean'

      allocate (domain    (im2,jn2,nl2))
      allocate (udomain   (im2,jn2,nl2))
      udomain(:,:,:)    =  ypert%u
      allocate (vdomain   (im2,jn2,nl2))
      vdomain(:,:,:)    =  ypert%v
      allocate (ptdomain  (im2,jn2,nl2))
      ptdomain(:,:,:)   =  ypert%pt
      allocate (qdomain   (im2,jn2,nl2))
      qdomain(:,:,:)    =  ypert%q(:,:,:,1)
      allocate (delpdomain(im2,jn2,nl2))
      delpdomain(:,:,:) =  ypert%delp

      call prognostics_final (ypert)

!   Compute kdomain to maintain as much code history as possible
      call comp_kdomain (im2, jn2, nl2, ptdomain, kdomain) 

!  Compute gradient of selected J
      print*
      if      (init_field(1:3) == 'u__' ) then
        domain(:,:,:) = udomain(:,:,:)
      else if (init_field(1:3) == 'v__' ) then
        domain(:,:,:) = vdomain(:,:,:)
      else if (init_field(1:3) == 'vpT' ) then
        domain(:,:,:) = ptdomain(:,:,:)
      else if (init_field(1:3) == 'vT_' ) then
        domain(:,:,:) = ptdomain(:,:,:)
      else if (init_field(1:3) == 'q__' ) then
        domain(:,:,:) = qdomain(:,:,:)
      else if (init_field(1:3) == 'dp_' ) then
        domain(:,:,:) = delpdomain(:,:,:)
      else if (init_field(1:3) == 'T__' ) then
        domain(:,:,:) = ptdomain(:,:,:)
      else if (init_field(1:3) == 'ps_' ) then
        domain(:,:,1:nl2-1) = 0.0
        domain(:,:,nl2) = delpdomain(:,:,nl2)
        kweights(nl2) = 1.0
        kweights(1:nl2-1) = 0.0
      endif

!  Gradients are analytic and related by a sum of the weights for these 3 cases
      if ( (trim(jtype) == 'bmx')     .or.           &  ! e.g. sum(Tw)    -> sum(w)    
           (trim(jtype) == 'bms' )    .or.           &  ! e.g. sum(Tw**2) -> 2*sum(T*w)  + sum(w)
           (trim(jtype) == 'rms' ) ) then               ! see subroutine 
        call find_name (nfields,fnames,trim(init_field),id_field)      ! determine field to compute gradient 
        if (id_field <= nfields3d ) then  ! field is 3-dimensional
          call grad_box_ (im2,jn2,nl2,fields3d1(:,:,:,id_field,2),      &
                          domain,kdomain,jweights,kweights,trim(jtype),       &
                         trim(init_field), fields3d1(:,:,:,id_field,1) )
        else                              ! field is 2-dimensional
          id_field=id_field-nfields3d
          call grad_box_ (im2,jn2,  1,fields2d1(:,:,:,id_field,2),      &
                          domain(:,:,nl2),kdomainps,jweights,           &
                         kweightps,trim(jtype),                          &
                         trim(init_field), fields2d1(:,:,:,id_field,1) )
        endif               

!  Gradients are computed using the energy norm for these five cases
!  Moist available enthalpy also computed in grad_enorm as u,v,p the same
      elseif ( ( trim(jtype) == 'txe' ) .or.        &
               ( trim(jtype) == 'twe' ) .or.        &
               ( trim(jtype) == 'mae' ) .or.        &
               ( trim(jtype) == 'kxe' ) .or.        &
               ( trim(jtype) == 'ape' ) .or.        &
               ( trim(jtype) == 'ate' )   ) then
        if (init_field(1:3) == '___' ) then

        !Grab zonal average of reference temperature and water 
        !vapor mixing ratio for available moist enthalpy norm
        if ( trim(jtype) == 'mae' ) then
          allocate(Tza(jn2,nl2) )      
          allocate(rza(jn2,nl2) )      

          !Get index for T and r
          call find_name (nfields,fnames,'T__',indexT)
          call find_name (nfields,fnames,'r__',indexR)

          if (nfiles>1) then
            Tza = sum(fields3d1(:,:,:,indexT,3),1)/im2
            rza = sum(fields3d1(:,:,:,indexR,3),1)/im2
          else
            Tza = sum(fields3d1(:,:,:,indexT,2),1)/im2
            rza = sum(fields3d1(:,:,:,indexR,2),1)/im2
          endif

        endif

        if (nfiles>1) then
          call grad_enorm_ (nymd,nhms,im2,jn2,nl2,nfields,ptop,R,Cp,eps_eer, &
                            udomain,vdomain,ptdomain,kdomain,jweights,  &
                            fields3d1(:,:,:,:  ,2),                     &
                            fields3d1(:,:,:,ndp,3),                     &
                            fields2d1(:,:,1,nps,2),fnames,              &
                            fields3d1(:,:,:,:  ,1),                     &
                            fields2d1(:,:,1,nps,1),trim(jtype),lu,      &
                            Tza,rza                                     )
        else
          call grad_enorm_ (nymd,nhms,im2,jn2,nl2,nfields,ptop,R,Cp,eps_eer, &
                            udomain,vdomain,ptdomain,kdomain,jweights,  &
                            fields3d1(:,:,:,:  ,2),                     &
                            fields3d1(:,:,:,ndp,2),                     &
                            fields2d1(:,:,1,nps,2),fnames,              &
                            fields3d1(:,:,:,:  ,1),                     &
                            fields2d1(:,:,1,nps,1),trim(jtype),lu,      &
                            Tza,rza                                     )
        endif

        if ( trim(jtype) == 'mae' ) then
          deallocate(Tza,rza)
        endif

        else
          print*,'init_field set incorrectly'
          stop
        endif
       
!  Gradient is computed based on mean circulation in the box
      elseif (trim(jtype) == 'mcx' ) then
        if (init_field(1:3) == '___' ) then
          call grad_circulation_(im2,jn2,nl2,                           &
                                udomain,vdomain,kdomain,jweights,      &
                                kweights,                              &
                                fields3d1(:,:,:,nu,1),                 &
                                fields3d1(:,:,:,nv,1))
        else
          print*,'init_field set incorrectly'
          stop
        endif

!  User made a boo-boo
      else
        print *,' **************************************************'
        print *,' NO VALID OPTION FOR JTYPE SELECTED '
        print *,' selected jtype = ',jtype 
      endif

!  Transform adjoint T field to adjoint virtual (potential) temperature field
      if (trans_T) then 
        if (nfiles>1) then
          call transform_T_adinv (im2,jn2,nl2,ptop,nfields, fnames,  &
                                  fields3d1(:,:,:,:,1:1),fields3d1(:,:,:,:,3:3) )
        else
          call transform_T_adinv (im2,jn2,nl2,ptop,nfields, fnames,  &
                                  fields3d1(:,:,:,:,1:1),fields3d1(:,:,:,:,2:2) )
        endif
      endif
 
!  Add contribution of adjoint ps field to adjoint delp field
      if (trans_ps) then 
        call adjoint_ps2delp_ (im2, jn2, nl2, fields3d1(:,:,:,ndp,1), &
                                              fields2d1(:,:,1,nps,1) )
      endif

!  Transform adjoint r water vapour mixing ratio field to adjoint q specific humidity field
      if (trans_q) then 
        if (nfiles>1) then
          call transform_q2r_ad (im2,jn2,nl2,nfields, fnames,fields3d1(:,:,:,:,1:1),fields3d1(:,:,:,:,3:3) )
        else
          call transform_q2r_ad (im2,jn2,nl2,nfields, fnames,fields3d1(:,:,:,:,1:1),fields3d1(:,:,:,:,2:2) )
        endif
      endif

! Apply LPO to DELP_AD
      if (kdomain(1).gt.  1) fields3d1(:,:,  1:kdomain(1)-1,ndp,1) = 0.d0
      if (kdomain(2).lt.nl2) fields3d1(:,:,kdomain(2)+1:nl2,ndp,1) = 0.d0

! Write output file with initialization fields
      summ = sum(fields3d1(:,:,:,1:nfields3d,1))

      allocate(aux2d(im2,jn2))
      aux4d (:,:,:,:) = fields3d1(:,:,:,:,1)
      aux2d (:,:)     = fields2d1(:,:,1,nps,1)
      call write_data (im2,jn2,nl2,lm1,nfields2d,nfields3d,ptop,ks,ak,bk, &
                       aux2d,aux4d,fnames,trim(ofname),nymd,nhms,'adm',idims=idims )
      deallocate(aux2d)

      summ = sum(fields3d1(:,:,:,1:nfields3d,1))
      print*,'wi=',summ

      deallocate(domain,udomain,vdomain,ptdomain,qdomain,delpdomain)
      deallocate(jweights,kweights,glats,pmid,delbk,ak,bk)
      deallocate(fields2d1,fields3d1)
      deallocate(aux4d)

      end subroutine init_adjoint_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: grad_box --- The adjoint fields are prescribed as the gradient 
!                        of the mean, mean square, or rms value of a particular
!                        field (specified as init\_field) where the mean is 
!                        computed as a horizontal area and delta\_sigma within
!                        some grid-rectangular volume.
!                        Gradients are analytic. ONLY Three options are available:
!                        'box mean'            dy/dx = factor1 + factor2
!                        'box mean square'     dy/dx = factor1 + factor2
!                        'root mean square'    dy/dx = factor1 + factor2 + factor3
!                        increasing down via the Chain Rule in complexity
!
!
!
! !INTERFACE:
!
      Subroutine grad_box_ (im1,jn1,nl1,field,          &
                           domain, kdomain,          &
                           jweights,kweights,jtype,   &
                           fname,field_ad             )

!USES:
 
      implicit none
                                                                                                              
! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1        ! grid dimensions
      integer, intent(in) :: kdomain(2)         ! indices defining range of k
      character(len=*), intent(in) :: fname     ! name of field that defines J
      character(len=*), intent(in) :: jtype     ! type of forecast measure J
      real(r8), intent(in)  :: domain(im1,jn1,nl1)   ! =1 if in box, =0 if out
      real(r8), intent(in)  :: jweights(jn1,2)       ! dlat*cos(lat)/2.
      real(r8), intent(in)  :: kweights(nl1)         ! layer weights (depths)
      real(r8), intent(in)  :: field(im1,jn1,nl1)    ! ref. or diff. field

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: field_ad(im1,jn1,nl1) ! adjoint field to create
                                                                                                                
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

      integer  :: i,j,k,jfirst,jlast
      integer  :: id1
      real(r8) :: kweight_sum   ! delta sigma of vertical span of box
      real(r8) :: afac          ! (2*pi/im1)/(4*pi) so that area=fraction of globe
      real(r8) :: area          ! fraction of global area covered by selected domain
      real(r8) :: fractn        ! fraction of box volume in current grid box
      real(r8) :: volume        ! total area*dsigma volume of box / area of globe
      real(r8) :: f2int         ! volume mean square of field in box
      real(r8) :: rms           ! rms of field in box
      real(r8) :: r_rms2        ! d(sqrt(x))/dx for x=f2int
      real(r8) :: fc
      logical  :: aflag(nl1)   ! flag that tells if lat lon area has been computed for a given level

!
      id1=2           ! pointer to jweights
      if (fname(1:1) == 'u' ) id1=1
                  
!  Factor 1: Compute total area and volume in box using the Chain Rule.
!  All 3 options have this first summation
      afac=1.d0/dble(im1)
      area=0.d0
      do j=1,jn1
        do i=1,im1
           if (domain(i,j,kdomain(1)) .ne. 0.0 ) then
             area=area+jweights(j,id1)*afac
           endif
        enddo
      enddo

      kweight_sum=0.d0
      do k=kdomain(1),kdomain(2)
        kweight_sum=kweight_sum+kweights(k)
      enddo
      volume=kweight_sum*area
                  
      f2int=0.d0
      do k=kdomain(1),kdomain(2)
        do j=1,jn1
          fractn=afac*kweights(k)*jweights(j,id1)/volume
          do i=1,im1
            if (domain(i,j,k) .ne. 0.0) then
                  
! Factor 2: This form of the second factor is computed only for 'box mean square' and 'root mean square'
              if  ( (trim(jtype) == 'bms') .or. &
                    (trim(jtype) == 'rms') )then
                fc=fractn*field(i,j,k)
                field_ad(i,j,k)=2.0d0*fc
                f2int=f2int+fc*field(i,j,k)
                  
! Factor 2: The second factor is different when computed only for 'box mean'
              elseif (trim(jtype) == 'bmx') then
                field_ad(i,j,k)=fractn
              endif
            endif
          enddo
        enddo
      enddo
             
! Factor 3: Only 'root mean square' has this third factor
      if (trim(jtype) == 'rms') then
        if (f2int > 0.d0 ) then
          rms=sqrt(f2int)
          r_rms2=0.5d0/rms
        else                                            !!! jtype = 'box mean square ' or 'box mean'
          rms=0.d0
          r_rms2=0.d0
        endif
        field_ad(:,:,:)=r_rms2*field_ad(:,:,:)
      endif
                  
      print *,' J set to gradient of ',jtype
      print *,'    field = ',fname
      print *,'    fractional area  = ',area
      print *,'    fractional depth = ',kweight_sum
      if (trim(jtype) == 'rms') then
        print *,'    J = ',rms
      endif
      if (trim(jtype) == 'bms') then
        print *,'    J = ',f2int
      endif
                  
      end subroutine grad_box_

!
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: grad_circulation --- compute the gradient of the mean vorticity
!
!
! !INTERFACE:
!
      subroutine grad_circulation_ (im1,jn1,nl1,                        &
                                   udomain,vdomain,kdomain,            &
                                   jweights,kweights,                  &
                                   u_ad,v_ad)

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: im1,jn1,nl1          ! grid dimensions
      integer, intent(in) :: kdomain(2)            ! indices defining range of k
      real(r8), intent(in) :: udomain(im1,jn1,nl1) ! 1 or 0 defines box
      real(r8), intent(in) :: vdomain(im1,jn1,nl1) ! 1 or 0 defines box
      real(r8), intent(in)   :: jweights(jn1,2)    ! dlat*cos(lat)/2
      real(r8), intent(in)   :: kweights(nl1)      ! layer weights (depths)

! !OUTPUT PARAMETERS:
 
      real(r8), intent(out)  :: u_ad(im1,jn1,nl1)  ! d(mean vort)/du
      real(r8), intent(out)  :: v_ad(im1,jn1,nl1)  ! d(mean vort)/dv
 

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!  The adjoint fields are prescribed as the gradient of the mean vorticity
!  within a volume.  Vorticity is computed on the T grid.
!
!  vorticity = 1/(r*cosp) ( (v\_e-v\_w)/dlon - (u*cose\_n - u*cos\_s)/dlat)
!  area of T grid point = r**2 cosp*dlat*dlon
!  where cosp=cos value of T grid point, cose is value on edge of grid box.
!  e,w,n,s aabove refer to edges of box
!  vorticity*grid-point-area = r*((v\_e-v\_w)*dlat - (u*cose\_n - u*cos\_s)*dlon )
!
!  Note that jweight(*,1)=0.5*cose*dlat, jweight(*,2)=0.5*cosp*dlat

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  24Jan2014  Todling    Use hermes to define constants (earth radius)
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer  :: i,j,k,kk
      real(r8) :: kweight_sum
      real(r8) :: area        ! horiz. area of box as fraction of earth's surface
      real(r8) :: pi
      real(r8) :: dlat        ! latitudinal spacing in radians
      real(r8) :: dlon        ! longitudinal spacing in radians
      real(r8) :: latfac
      real(r8) :: lonfac
      real(r8) :: afac
      real(r8) :: fac
      real(r8) :: earthr
      real(r8) :: volume      ! volume of box / volume of atmosphere
      logical  :: aflag(nl1)   ! flag that tells if lat lon area has been computed for a given level
                  
!
      earthr = radius_earth
      pi=4.d0*datan(1.d0)
      dlat=pi/dble(jn1-1)
      dlon=2.0d0*pi/dble(im1)
      latfac=earthr*dlat
      lonfac=earthr*dlon/dlat ! remove dlat factor already in jweight.
      afac=latfac*lonfac
                  
      kweight_sum=0.d0
      do k=kdomain(1),kdomain(2)
        kweight_sum=kweight_sum+kweights(k)
      enddo
             
! South Pole
      area = 0.0
      k=kdomain(1)
      do i=1,im1
        if (udomain(i,1,k) .ne. 0.0) area=area+jweights(1,2)
        do kk=kdomain(1),kdomain(2)
          if (udomain(i,1,k) .ne. 0.0 ) then
            u_ad(i,2,k)=-lonfac*jweights(2,1)*kweights(k) 
          endif
        enddo
      enddo
             
! North pole
      k=kdomain(1)
      do i=1,im1
        if (udomain(i,jn1,k) .ne. 0.0) area=area+jweights(jn1,2)
        do kk=kdomain(1),kdomain(2)
          if (udomain(i,jn1,k) .ne. 0.0 ) then
            u_ad(i,jn1,k)=lonfac*jweights(jn1,1)*kweights(k)
          endif
        enddo
      enddo
             
! Non-Pole points
      kk = kdomain(1)
      do j=2,jn1-1

!  Special treatment of i=im1 since next eastward point is i=1
        if (udomain(im1,j,kk) .ne. 0.d0) then
          area=area+jweights(j,2)*afac
          do k=kdomain(1),kdomain(2)
            u_ad(im1,j,k)  = u_ad(im1,j,k)   + lonfac*kweights(k)*jweights(j,  1)
            u_ad(im1,j+1,k)= u_ad(im1,j+1,k) - lonfac*kweights(k)*jweights(j+1,1)
            v_ad(im1,j,k)  = v_ad(im1,j,k)   - latfac*kweights(k)
            v_ad(1,j,k)    = v_ad(1,j,k)     + latfac*kweights(k)
          enddo
        endif

!  rest of the i's
        do i=1,im1-1
          if (udomain(i,j,kk) .ne. 0.d0) then
            area=area+jweights(j,2)*afac
            do k=kdomain(1),kdomain(2)
              u_ad(i,j,k)  = u_ad(i,j,k)   + lonfac*kweights(k)*jweights(j,  1)
              u_ad(i,j+1,k)= u_ad(i,j+1,k) - lonfac*kweights(k)*jweights(j+1,1)
              v_ad(i,j,k)  = v_ad(i,j,k)   - latfac*kweights(k)
              v_ad(i+1,j,k)= v_ad(i+1,j,k) + latfac*kweights(k)
            enddo
          endif
        enddo

      enddo  ! loop over j
!
! divide by total volume of box
      volume=area*kweight_sum
      if (volume > 0.d0 ) then
        u_ad(:,:,:)=u_ad(:,:,:)/volume
        v_ad(:,:,:)=v_ad(:,:,:)/volume
      endif
                  
      print *,' J set to derivative of mean circulation in a box'
      print *,'    vert. range = ',kdomain(1),kdomain(2)
      print *,'    horizontal area = ',area
                  
      end subroutine grad_circulation_
                  
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: grad_enorm --- adjoint fields are prescribed as the gradient of 
!                          some kind of energy norm within some volume.
!
!
! !INTERFACE:
!
      subroutine grad_enorm_ (nymd,nhms,im1,jn1,nl1,nfields,ptop,R,Cp,eps_eer, &
                              udomain,vdomain,ptdomain,kdomain,jweights, &
                              fields,delp_ref,psfield,fnames,            &
                              fields_ad,psfield_ad,jtype,lunit,Tza,rza) 

!USES:
 
      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: nymd,nhms               ! date and time
      integer, intent(in) :: im1,jn1,nl1             ! grid dimensions
      real(r8), intent(in)  :: ptop
      real(r8), intent(in)  :: R
      real(r8), intent(in)  :: Cp
      real(r8), intent(in)  :: eps_eer
      integer, intent(in) :: nfields                ! number of fields 
      integer, intent(in) :: lunit 
      integer, intent(in) :: kdomain(2)             ! range of k to consider
      character(len=*), intent(in) :: fnames(nfields) 
      character(len=*), intent(in) :: jtype         ! type of forecast measure J
      real(r8), intent(in)  :: jweights(jn1,2)      ! dlat(radians)*cos(lat)/2
      real(r8), intent(in)  :: udomain(im1,jn1,nl1)   ! 1s or 0s that define box
      real(r8), intent(in)  :: vdomain(im1,jn1,nl1)   ! 1s or 0s that define box
      real(r8), intent(in)  :: ptdomain(im1,jn1,nl1)   ! 1s or 0s that define box
      real(r8), intent(in)  :: fields(im1,jn1,nl1,nfields)  ! ref. or diff fields
      real(r8), intent(in)  :: delp_ref(im1,jn1,nl1)        ! ref. delta p 
      real(r8), intent(in)  :: psfield(im1,jn1)            ! ref. or diff ps
      real(r8), optional, intent(in) :: Tza(jn1,nl1)      ! ref zonal average of T
      real(r8), optional, intent(in) :: rza(jn1,nl1)      ! ref zonal average of r_v

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
!                        replaced  $do i=i,im1$ to be  $do i=1,im1$ in KEV looping 
!                        Corrected 0.5 factors for sum(jweights)=1
!  15Aug2005  Elena N.   Output of Jnorm in a file is included now
!  12Dec2007  Todling    Add moisture to energy norm
!  10Oct2008  Todling    Add V-norm opt (as from Errico et al., 2007)
!  11Sep2013  Todling    Capability to handle winds as on A-grid
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      character(len=*), parameter :: myname_ = myname//'*grad_enorm_'

      integer  :: i,j,k,kk,nlp1,ierr
      integer  :: idu,idv,idt,idq
      real(r8) :: wfactor
      real(r8) :: dbleimr
      real(r8) :: tfac, ufac, pfac, qfac
      real(r8) :: dsigma(im1,jn1,nl1)  ! sigma-thickness of layer
      real(r8) :: ps_ref(im1,jn1)     ! reference surface pressure
      real(r8) :: jfac(jn1)           ! fractional area of globe in T grid box 
      real(r8) :: fa(im1,jn1) 
      real(r8) :: fb(im1,jn1)         ! (interpolated) squared values of field at T points
      real(r8) :: fa_ad(im1,jn1)      ! d(E)/d(fa) without ufac,pfac,tfac 
      real(r8) :: fb_ad(im1,jn1)      ! d(E)/d(fb) without ufac,pfac,tfac 
      real(r8) :: eu,ev,et,ep,eq,e(5),esum
      real(r8), allocatable :: pe(:,:,:)
      real(r8), allocatable :: tfac_mae(:,:,:), qfac_mae(:,:,:)
      character(len=1) wgrid
      integer vnorm
      integer kz(1)

! Compute fractional areas
      dbleimr=1.d0/dble(im1)
      do j=1,jn1
        jfac(j)=jweights(j,2)*dbleimr  ! fractional area of T grid points
      enddo

! Assumption for division by area where sum(jweights)=1
      tfac=0.5d0*Cp/tref
      ufac=0.5d0
      qfac=0.5d0*eps_eer*alhl*alhl/(Cp*tref)
      pfac=0.5d0*R*tref/pref**2
      e=0.d0 

      if ( trim(jtype) == 'mae' ) then

        allocate( tfac_mae(im1,jn1,nl1), stat=ierr )
        allocate( qfac_mae(im1,jn1,nl1), stat=ierr )

        do i=1,im1
          tfac_mae(i,:,:) = 0.5d0*(Cp/Tref)*(Tref/Tza)**2
          qfac_mae(i,:,:) = 0.5d0*461.52*Tref/rza
        enddo
      endif

      call getparam ( 'wgrid', wgrid )
      call getparam ( 'vnorm', vnorm )
      print *, trim(myname),' operator treating winds as on ',wgrid,'-grid'

      call find_name (nfields,fnames,'u__',idu)      
      call find_name (nfields,fnames,'v__',idv)      
      if ( ( trim(jtype) == 'txe' ) .or. &
           ( trim(jtype) == 'twe' ) .or. &
           ( trim(jtype) == 'mae' ) .or. &
           ( trim(jtype) == 'ape' ) .or. &
           ( trim(jtype) == 'ate' ) ) then
        call find_name (nfields,fnames,'T__',idt)     ! it is assumed that vpt->T has already happened 
      endif
      if ( trim(jtype) == 'twe'  ) then
        call find_name (nfields,fnames,'q__',idq)
      endif
      if ( trim(jtype) == 'mae' ) then
        call find_name (nfields,fnames,'r__',idq)
      endif!
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
!  Calculate total energy by looping over levels, then fields
!
      do k=kdomain(1),kdomain(2)

        do j=1,jn1
          do i=1,im1
            fb_ad(i,j)=dsigma(i,j,k)*jfac(j)
          enddo
        enddo
!     
!  Energy calculations for u and v fields e(1) and e(2) respectively
!  Contributions from winds are only computed for 'total energy' and
!  'kinetic energy'. All other options: 'available potential energy'  
!  and 'available thermal energy' do not consider wind fields
        if ( ( trim(jtype) == 'txe' ) .or. & 
             ( trim(jtype) == 'twe' ) .or. &
             ( trim(jtype) == 'mae' ) .or. &
             ( trim(jtype) == 'kxe' ))  then
      
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              if (udomain(i,j,k) .ne. 0.0) fa(i,j)=fields(i,j,k,idu)
            enddo
          enddo

          if(wgrid=='d') then
             do i=1,im1
               fa_ad(i,  2)=fa(i,  2)*( 2.0d0*fb_ad(i,  1)+fb_ad(i,    2) )
               fa_ad(i,jn1)=fa(i,jn1)*( 2.0d0*fb_ad(i,jn1)+fb_ad(i,jn1-1) )
               fb(i,  1)   =fa(i,  2)*fa(i,  2)
               fb(i,jn1)   =fa(i,jn1)*fa(i,jn1)
               fb(i,2)     =0.5*(fa(i,2)*fa(i,2) + fa(i,3)*fa(i,3))
               do j=3,jn1-1     
                 fa_ad(i,j)=fa(i,j)*( fb_ad(i,j)+fb_ad(i,j-1) )
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
                 fa_ad(i,j)=2.0d0*fb_ad(i,j)*fa(i,j)
                 eu=eu+fb_ad(i,j)*fa(i,j)*fa(i,j)
               enddo
             enddo
          endif
          e(1)=e(1)+ufac*eu
 
          fa_ad(:,1)=0.d0
          fields_ad(:,:,k,idu)=ufac*fa_ad(:,:)

!  Energy calculations for v field e(2)
!  Note that fb grid is assumed to be offset to east of v grid point
          fa = 0.0
          do i=1,im1
            do j=2,jn1-1
              if (vdomain(i,j,k) .ne. 0.0) fa(i,j)=fields(i,j,k,idv)
            enddo
          enddo

          if(wgrid=='d') then
             do j=2,jn1-1
               fa_ad(1,j)=fa(1,j)*( fb_ad(im1,j)+fb_ad(1,j) )
               fb(im1,j)=0.5*(fa(1,j)*fa(1,j) + fa(im1,j)*fa(im1,j))
               do i=1,im1-1  
                 fa_ad(i+1,j)=fa(i+1,j)*( fb_ad(i,j)+ fb_ad(i+1,j) ) 
                 fb(i,j)=0.5*(fa(i,j)*fa(i,j) + fa(i+1,j)*fa(i+1,j))
               enddo  
             enddo
             do i=1,im1
               fa_ad(i,    2)=2.0d0*fa(i,    2)*fb_ad(i,  1)+fa_ad(i,    2)
               fa_ad(i,jn1-1)=2.0d0*fa(i,jn1-1)*fb_ad(i,jn1)+fa_ad(i,jn1-1)
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
                 fa_ad(i,j)=2.0d0*fb_ad(i,j)*fa(i,j)
                 ev=ev+fb_ad(i,j)*fa(i,j)*fa(i,j)
               enddo
             enddo
          endif
          e(2)=e(2)+ufac*ev

          fa_ad(:,  1)=0.d0
          fa_ad(:,jn1)=0.d0
          fields_ad(:,:,k,idv)=ufac*fa_ad(:,:)
        endif
!     
!  Energy calculations for T field e(3)
!  T form is assumed
!     
        if ( ( trim(jtype) == 'txe' ) .or.   & 
             ( trim(jtype) == 'twe' ) .or.   &
             ( trim(jtype) == 'ape' ) .or.   &
             ( trim(jtype) == 'ate' ) ) then
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              if (ptdomain(i,j,k) .ne. 0.0) fa(i,j)=fields(i,j,k,idt)
            enddo
          enddo

          et=0.d0 
          do j=1,jn1
            do i=1,im1
              fa_ad(i,j)=2.0d0*fb_ad(i,j)*fa(i,j)
              et=et+fb_ad(i,j)*fa(i,j)*fa(i,j)
            enddo
          enddo
          e(3)=e(3)+tfac*et  

          fields_ad(:,:,k,idT)=tfac*fa_ad(:,:)
        endif

!     
!  Energy calculations for Q field e(5)
!     
        if ( trim(jtype) == 'twe'  ) then
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              if (ptdomain(i,j,k) .ne. 0.0) fa(i,j)=fields(i,j,k,idq)
            enddo
          enddo

          eq=0.d0 
          do j=1,jn1
            do i=1,im1
               fa_ad(i,j)=2.0d0*fb_ad(i,j)*fa(i,j)
               eq=eq+fb_ad(i,j)*fa(i,j)*fa(i,j)
            enddo
          enddo
          e(5)=e(5)+qfac*eq  

          fields_ad(:,:,k,idq)=qfac*fa_ad(:,:)
        endif ! test on whether 'Q' is included in this version of J

!
!  Moist available enthalpy calculations
! 
        if ( trim(jtype) == 'mae'  ) then

          !Temperature part
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              if (ptdomain(i,j,k) .ne. 0.0) fa(i,j)=fields(i,j,k,idt)
            enddo
          enddo

          !et=0.d0 
          do j=1,jn1
            do i=1,im1
              fa_ad(i,j)=2.0d0*fb_ad(i,j)*fa(i,j)
              !et=et+fb_ad(i,j)*fa(i,j)*fa(i,j)
              e(3) = e(3) + tfac_mae(i,j,k)*fb_ad(i,j)*fa(i,j)*fa(i,j)
            enddo
          enddo

          fields_ad(:,:,k,idT)=tfac_mae(:,:,k)*fa_ad(:,:)

          !Water vapour mixing ratio part
          fa = 0.0
          do i=1,im1
            do j=1,jn1
              if (ptdomain(i,j,k) .ne. 0.0) fa(i,j)=fields(i,j,k,idq)
            enddo
          enddo

          !eq=0.d0 
          do j=1,jn1
            do i=1,im1
               fa_ad(i,j)=2.0d0*fb_ad(i,j)*fa(i,j)
               !eq=eq+fb_ad(i,j)*fa(i,j)*fa(i,j)
               e(5)=e(5)+qfac_mae(i,j,k)*fb_ad(i,j)*fa(i,j)*fa(i,j)
            enddo
          enddo

          fields_ad(:,:,k,idq)=qfac_mae(:,:,k)*fa_ad(:,:)

        endif

      enddo ! loop over k


      if ( trim(jtype) == 'mae' ) then

        deallocate( tfac_mae, stat=ierr )
        deallocate( qfac_mae, stat=ierr )

      endif


!
! Energy calculations for ps field e(4)
!
      if ( ( trim(jtype) == 'txe' ) .or. & 
           ( trim(jtype) == 'twe' ) .or. &
           ( trim(jtype) == 'mae' ) .or. &
           ( trim(jtype) == 'ape' ))  then

        fa = 0.0
        kk = kdomain(1)
        do i=1,im1
          do j=1,jn1
            if (ptdomain(i,j,kk) .ne. 0.0) fa(i,j)=psfield(i,j)
          enddo
        enddo

        ep=0.d0
        do j=1,jn1
          do i=1,im1
            fa_ad(i,j)=2.0d0*jfac(j)*fa(i,j)
            ep=ep+jfac(j)*fa(i,j)*fa(i,j)
          enddo
        enddo
        e(4)=pfac*ep

        psfield_ad(:,:)=pfac*fa_ad(:,:)
      endif

      esum=sum(e)

      if ( trim(jtype) /= 'mae'  ) then

      print *,' J set to derivative of ',jtype
      print *,'    vert. range = ',kdomain(1),kdomain(2)
      print *,'    J of ref or dif field = ',esum      
      print *,'    J of ref or dif components (u,v,T,ps,q) = ',e(1:5)
!
      write (lunit, '(a,i8.8,a,i6.6,2a)') 'Date: ', nymd, ' Time: ', nhms, ' Norm type: ', trim(jtype)
      write (lunit, '(a)') 'sum, u, v, t, ps, q'
      write (lunit, '(6f20.12)') esum, e(1:5)
   
      else

      print *,' J set to derivative of ',jtype
      print *,'    vert. range = ',kdomain(1),kdomain(2)
      print *,'    J of ref or dif field = ',esum      
      print *,'    J of ref or dif components (u,v,T,ps,r) = ',e(1:5)
!
      write (lunit, '(a,i8.8,a,i6.6,2a)') 'Date: ', nymd, ' Time: ', nhms, ' Norm type: ', trim(jtype)
      write (lunit, '(a)') 'sum, u, v, t, ps, r'
      write (lunit, '(6f20.12)') esum, e(1:5)

      endif

!
      end subroutine grad_enorm_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: adjoint_ps2delp --- adjoint of routine calculate ps
!
!
! !INTERFACE:
!
      subroutine adjoint_ps2delp_ (im1, jn1, nl1, delp_ad, ps_ad)
!
! add contribution to adjoint delp field by adjoint ps field
! This is the adjoint of routine calculate ps
!    

!USES:
 
      implicit none

! !INPUT PARAMETERS:
 
      integer,  intent(in)     :: im1,jn1,nl1
      real(r8), intent(in)     :: ps_ad(im1,jn1)      ! surface pressure
     
! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
      real(r8), intent(inout)  :: delp_ad(im1,jn1,nl1) ! p-thickness of layers

! !DESCRIPTION:
!
! add contribution to adjoint delp field by adjoint ps field
! This is the adjoint of routine calculate ps
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
 
      integer  :: i,j,k

      do j=1,jn1
        do i=1,im1
          do k=1,nl1 
            delp_ad(i,j,k)=ps_ad(i,j)+delp_ad(i,j,k)
          enddo
        enddo
      enddo

      end subroutine adjoint_ps2delp_ 

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PertSet_ --- General paramater setup for singular vector calc
! 
! !INTERFACE:
!
      subroutine PertSet_ ( comm, root, diff_file, jtype, init_field, &
                                 knmax, knames, RCfile, vectype, stat )
 
! !USES: 
    
      Implicit None

! !INPUT PARAMETERS: 
!
      integer, intent(in) ::  comm	! MPI communicator
      integer, intent(in) ::  root	! MPI ROOT PE
      integer, intent(in) ::  knmax	! dim(knames)
      integer, intent(in) ::  vectype   ! GEOS-4/5 vector specifier

      character(len=*), intent(in) :: RCfile  ! resource filename

! !OUTPUT PARAMETERS:
 
      logical, intent(out) :: diff_file     ! true if difference of two files
      character(len=*), intent(inout) :: knames(knmax)! field: U, V, PT, etc
      character(len=*), intent(out) :: jtype(:)     ! norm type
      character(len=*), intent(out) :: init_field   ! field: U, V, PT, etc
      integer, intent(out) :: stat            ! Return error code:
                                              !  stat=1 cannot determined proc
                                              !  stat=2 RC file not found
                                              !  stat=3 missing RC name indicator
                                              !  stat=4 missing RC entry


! !FILES USED:  fvpert.rc
!
! !DESCRIPTION:  Initializes perturbation used for sensitivity runs.
!
! !REVISION HISTORY: 
!
!   30Jul2004  Winslow    Initial code, based on similar codes
!   07May2007  Todling    No longer reads date and time.
!   27Oct2007  Todling    Support to GEOS-5-type input/output vectors
!   10Oct2008  Todling    Add vnorm opt
!   11Sep2013  Todling    Set a-grid winds when applicable
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'PertSet_'

      character(len=255) token
      character(len=255) pertrc
      character(len=3  ) evolve_svec

      real(r8)      eps_eer
      integer       i, j, k, jcnt, maxjtype, iret, ierr
      integer       vnorm
      integer       myID

      stat = 0

      call MP_comm_rank(comm,myID,ierr)
        if (ierr/=0) then
            stat = 1  ! cannot grab proc id
            return
        endif

!     Set input/output vector type
!     ----------------------------
      call setparam ( 'vectype', vectype )
      if( vectype==5 ) then
          call setparam ( 'tname', 'vT_' )
          call setparam ( 'wgrid', 'a'   )
      endif


!     Load resources from fvpert.rc
!     -----------------------------
      pertrc = ' '
      call getenv('FVPERT_RC',pertrc)		! Unix binding
      if(pertrc.eq.' ') pertrc=RCfile	        ! default name
      call i90_loadf (trim(pertrc), iret)
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

      if ( knmax >= 1 ) then
           knames(1) = 'ADM1'
      else
        if(myID==ROOT) write(stdout,'(2a,i2)') myname_, ': unacceptable knmax ', knmax
        stat = 6
        return
      endif

!     Decide if you want to difference two files
!     ------------------------------------------
        call I90_label('diff_file:', iret)
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
        diff_file = .false.
        if( trim(token) .eq. 'yes' ) diff_file = .true.
        if(myID==ROOT) write(stdout,'(2a)')   &
                     'Difference of two files: ', trim(token)

!     Decide between E- and V-norms
!     -----------------------------
      call I90_label('vnorm:', iret)
      if (iret .ne. 0) then
        write(stdout,'(2a)') myname_,    &
                     ': I90_label not found will use default '
      else
        vnorm = I90_GInt(iret)
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

!     Define type of File1
!     --------------------
      if ( knmax >= 2 ) then
         knames(2) = 'NLM1'   ! default: NLM1
         call I90_label('file_type1:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             knames(2) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File type for first file: ', trim(knames(2))
      endif

!     Define type of File2
!     --------------------
      if ( knmax == 3 ) then
         knames(3) = 'NLM2'   ! default: NLM1
         call I90_label('file_type2:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             knames(3) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File type for second file: ', trim(knames(3))
      endif
                                                                                                                            
!     Define type of norm
!     -------------------
      maxjtype = size(jtype)
      jcnt      = 0
      jtype(1)  = 'txe'   ! default: TE
                          ! other possibilities are:
                          !  twe - total wet energy (q contribution)
                          !  mae - moist available enthalpy (T and q different from txe/twe)
                          !  ape - available potential energy
                          !  ate - available thermal energy
                          !  mcx - mean circulation
                          !  bmx - box mean
                          !  bms - box mean square
                          !  rms - root mean square
      jtype(2:maxjtype) = 'NONE'
      call I90_label('pert_norm:', iret)
      if (iret .eq. 0) then
          do while ( iret==0 .and. jcnt<maxjtype )
             call I90_Gtoken ( token, iret )
             if ( iret==0 ) then
                  jcnt = jcnt + 1
                  jtype(jcnt) = lowercase(trim(token))
                  if(myID==ROOT) write(stdout,'(2a)') 'Norm for perturbation: ', trim(jtype(jcnt))
             endif
          enddo
      end if

!     Read Ehrendorfer, Errico, and Raeder's epsilon factor (apply to Q-component)
!     ----------------------------------------------------------------------------
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
      if(myID==ROOT) write(stdout,'(a,1p,e13.6)') 'Ehrendorfer, Errico, and Raeder eps: ',eps_eer
      call setparam ( 'eps_eer', eps_eer )

!     Define field used by norm
!     -------------------------
      init_field = '___' ! default
      call I90_label('init_field:', iret)
      if (iret .eq. 0) then
        call I90_Gtoken ( token, iret )
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_, &
                       ': I90_Gtoken error, iret =',iret
        else
          init_field = trim(token)
        end if
      end if
      if(myID==ROOT) write(stdout,'(2a)') 'Field for perturbation: ', trim(init_field)


!     release resource file:
!     ---------------------
      call I90_release()

      return
      end subroutine PertSet_

      subroutine clean_ ()
      end subroutine clean_

      end module m_initadj
