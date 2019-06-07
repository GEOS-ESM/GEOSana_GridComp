! -------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_postStats --- 
!
! !INTERFACE:

      module m_postStats

      use precision

      use prognostics, only : prognostics_initial
      use prognostics, only : prognostics_final
      use prognostics, only : prognostics_cnst
      use prognostics, only : dyn_prog

      use m_const, only : Cp     => cpm
      use m_const, only : R      => rgas
      use m_const, only : pstd    
      use m_const, only : tref   => tstd
      use m_const, only : kappa  
      use m_const, only : zvir
                                                                                                         
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type
                                                                                                        
      use m_mapz, only : set_eta

      use m_svprojs, only : proj_init
      use m_svprojs, only : proj_clean
      use m_svprojs, only : proj_svec

      use m_pertutil, only : pt2vpt => pertutil_pt2vpt
      use m_pertutil, only : ptv2pt => pertutil_ptv2pt
      use m_pertutil, only : pt2t   => pertutil_pt2t
      use m_pertutil, only : pt2t_tl => pertutil_pt2t_tl
      use m_pertutil, only : t2pt_ad => pertutil_t2pt_ad
      use m_pertutil, only : calculate_ps => pertutil_calculate_ps
      use m_pertutil, only : calculate_ps_ad => pertutil_calculate_ps_ad
      use m_pertutil, only : find_name => pertutil_find_name
      use m_pertutil, only : horiz_grid => pertutil_horiz_grid
      use m_pertutil, only : vert_grid => pertutil_vert_grid
      use m_pertutil, only : pkappa => pertutil_pkappa
      use m_pertutil, only : comp_kdomain
      use m_pertutil, only : transform_t => pertutil_transform_t
      use m_pertutil, only : transform_t_tl => pertutil_transform_t_tl
      use m_pertutil, only : transform_t_ad => pertutil_transform_t_ad
      use m_pertutil, only : read_data => pertutil_read_data
      use m_pertutil, only : getparam  => pertutil_getparam
      use m_pertutil, only : setparam  => pertutil_setparam

      use m_dyn

      use m_gfio_putfld
                                                                                                         
      use m_stdio
      use m_inpak90
      use m_die, only: MP_die, die
                                                                                                         
      implicit NONE
                                                                                                      
! !PUBLIC MEMBER FUNCTIONS:
                                                                                                   
      PRIVATE
                                                                                                           
      PUBLIC postStats_Set
      PUBLIC postStats_Process
      PUBLIC postStats_Clean

      interface postStats_Set ; module procedure StatsSet_   ; end interface
      interface postStats_Process; module procedure field_stats_; end interface
      interface postStats_Clean; module procedure StatsClean_; end interface
                                                                                                      
!
! !DESCRIPTION: Set up for adjoint initialization
!
! !REVISION HISTORY:
!
!  18Oct2004  Winslow   Initial code based on similar code
!  10May2005  Elena N.  Linked constants of this module to common constants 
!  05Jan2008  Todling   Trying to clean up; make use of m_pertutil:
!                        - calc_ps using version w/ ptop
!                        - using find_name w/ error check on
!
!EOP
!-------------------------------------------------------------------------
 
      character(len=*), parameter :: myname = 'm_postStats'
      real(r8), parameter :: pref = 100.d0*pstd 

      integer, save ::  nOutfiles                          ! dim(knames)
      character(len=255), allocatable, save :: knames(:)   ! field types (NLM, TLM)
      character(len=255), allocatable, save :: outNames(:) ! file names

      CONTAINS
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: field_stats_ --- Compute basic stats and do T transform
!
!
! !INTERFACE:
!
      subroutine field_stats_ (comm,ROOT, &
                               kindsIn,corr_kind,rescale_adj,e_calc, jobIn, &
                               ntimes,datesOut,timesOut,RCfile,pert_with_LPO)  

!USES:


      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: comm
      integer, intent(in) :: ROOT
      integer, intent(in) :: kindsIn          ! number of kinds of data read in
      integer, intent(in) :: corr_kind(2)     ! slots to be correlated
      logical, intent(in) :: rescale_adj      ! true if adm field comparison 
                                              ! should extend through the column
      logical, intent(in) :: e_calc ! true if energy norm to be calc.
      character(len=*), dimension(kindsIn), intent(in) :: jobIn 
      integer, intent(in) :: ntimes           ! number of output model times to be exam.
      integer, dimension(ntimes), intent(inout) :: datesOut 
      integer, dimension(ntimes), intent(inout) :: timesOut 
      character(len=*), intent(in) :: RCfile
      logical, intent(in) :: pert_with_LPO


! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
! kinds is the number of kinds of fields to be read or computed
! Note that a TLM# kind is always followed by a RTL# kind, that is the
! corresponding reference state about which the linearization has been
! performed.  These are fields defined at the same time. Similarly, an
! ADM# kind is followed by a RAD# kind. Lastly, the difference fields always go on the end.
! See the options below
!
!
! Options:
! (nOutfiles=2) (kindsIn=2) (/ 'TLM1', 'RTL1' /)  Examine only TLM and its ref. fields
! (nOutfiles=2) (kindsIn=2) (/ 'ADM1', 'RAD1' /)  Examine only ADM and its ref. fields
! (nOutfiles=3) (kindsIn=2) (/ 'NLM1', 'NLM2', 'DNN ' /) Compare 2 NLM results
! (nOutfiles=5) (kindsIn=4) (/ 'TLM1', 'RTL1', 'TLM2', 'RTL2', 'DTT ' /)  Compare 2 TLM results
! (nOutfiles=5) (kindsIn=4) (/ 'ADM1', 'RAD1', 'ADM2', 'RAD2', 'DAA '  /) Compare 2 ADM results
! (nOutfiles=6) (kindsIn=4) (/ 'TLM1', 'RTL1', 'NLM2', 'NLM3', 'DNN ', 'DTN ' /)
!   This last compares the difference of 2 NLM results with a TLM result
!
! Set the 2-valued integer parameter corr_kind
! If not 0,0 this will cause the corrleations of kinds value 1, value 2 to
! be computed (See the character list knames for what these 2 values mean)
!
! Set e_calc=.true. if the energy norm is to be computed for all TLM or
! difference fields.  In this case, the field-name specification for the
! temperature field should be 'T__' (see (7) below).
!
!
! Specify experiments and model times to examine
!
!       'u004_a18.tlm1.eta.20170123.nc4', &     ! file name for output of kname 1 
!       'u004_a18.rtl1.eta.20170123.nc4'  &     ! file name for output of kname 2 
!       'u004_a18.dtt3.eta.20170123.nc4', &     ! file name for output of kname 3 
!       'u004_a18.dtt4.eta.20170123.nc4', &     ! file name for output of kname 4 
!       'u004_a18.dtt5.eta.20170123.nc4', &     ! file name for output of kname 5 
!       'u004_a18.dtt6.eta.20170123.nc4', &     ! file name for output of kname 6 
!       /)
!       'u000_a18.traj.eta.20170123.nc4', &     ! file name for kname reference field 1
!       'u000_a18.traj.eta.20170123.nc4'  &     ! file name for kname reference field 2
!       'u000_a18.traj.eta.20170123.nc4', &     ! file name for kname reference field 3
!       'u000_a18.traj.eta.20170124.nc4', &     ! file name for kname reference field 4
!       /)
!
! Specify the specific form of temperature to be examined.  
! The FVGCM field is 'vpT' 
! for virtual potential temperature. If the derived potential 
! temperature will be examined a 'PT_' is specified.  
! 'T__' means temperature itself will be examined.
! 'ps_' is the surface pressure derived from the model 'dp_' (delp) field.
! The remaining values here generally should not be changed, but allow 
! easy generalization to include additional model or derived fields.
!
!  The code now computes 5 measures only (Min, max, mean, rms, stdv)
!  in addition to the optional energy norm and field correlations. So,
!  measures =5.  This parameter allows easy generalization to additional
!  metrics.
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!  16May2007  Todling    Introduced dyn_prog.
!  05Nov2007  Todling    Quick fix to handle knames/outNames
!  06Mar2009  Todling    .hdf to .nc4
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer, parameter :: nfields3d=5   ! number of 3D fields to be examined
      integer, parameter :: nfields2d=1   ! number of 2D fields to be examined
      integer, parameter :: nfields=nfields3d+nfields2d
      character(len=3) :: tname_orig, tname 
      character(len=3), dimension(nfields) :: fnames 
      character(len=3), dimension(nfields), parameter :: fnames_def = &
                  (/ 'u__','v__','vpT','q__','dp_','ps_' /)
      integer :: im1   ! grid points in zonal direction from reference file(s)
      integer :: jn1   ! grid points in lat. direction from reference file(s)
      integer :: nl1   ! grid points in vert. direction from reference file(s)
      integer :: lm1   ! number of tracers from reference file(s)
      integer,  dimension(2), parameter :: kdomainps=(/1,1/) ! Do not change

      integer, parameter :: measures=5  ! number of kinds of stats examined
      character(len=4), dimension(measures), parameter :: mnames = &
               (/' MAX',' MIN','MEAN',' RMS','STDV'/)

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

      integer :: i,j,k,nf         ! looping indices
      integer :: ier

! Counters for do-loops and for other array indicies
      integer :: ntimex
      integer :: kindx
      integer :: countKT
      integer :: nkind0
      integer :: nkind1
      integer :: nlp1
      integer :: ndT, ndp, nps, kdp  ! indicates particular fields

      real(r8), allocatable :: fields2d(:,:,:,:,:)
      real(r8), allocatable :: fields3d(:,:,:,:,:)
      real(r8), allocatable :: fields3d_LPO(:,:,:,:)

      real(r8), allocatable :: stats3d(:,:,:,:,:)
      real(r8), allocatable :: stats2d(:,:,:,:,:)
      real(r8), allocatable :: energy(:,:,:,:)
      real(r8), allocatable :: correls3d(:,:,:)
      real(r8), allocatable :: correls2d(:,:,:)
       
      integer  :: kdomain(2)      ! indexes corresp to pdomain
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

! 
!                      BEGIN EXECUTION
!
!  Read in reference field(s)
      call dyn_getdim   (jobIn(1),im1,jn1,nl1,lm1,ier)

      nlp1 = nl1 + 1

!
!  construct grid mask of local projection operator (LPO)
!
      print*
      print*,'Constructing LPO:'
      projrc = 1
                                                                                                                            
      dims(1) = im1
      dims(2) = jn1
      dims(3) = nl1
      dims(4) = 1
      dims(5) = jn1
                                                                                                                            
      m = 5*im1*jn1*nl1
      call prognostics_initial (ypert,im1,jn1,nl1,1)
      call prognostics_cnst ( 1.d0,ypert )

      call proj_init (  comm, root, stat=projrc, dims=dims, RCfile=RCfile )
      if(projrc/=0) print*,'projrc = ',projrc,' proj_init'

      call proj_svec ( comm, ROOT, ypert, projrc )   ! Apply local projection operator
      if(projrc/=0) print*,'projrc = ',projrc,' proj_svec'

      call proj_clean( projrc )
      if(projrc/=0) print*,'projrc = ',projrc,' proj_clean'

      allocate (udomain   (im1,jn1,nl1))
      udomain(:,:,:) =  ypert%u
      allocate (vdomain   (im1,jn1,nl1))
      vdomain(:,:,:) =  ypert%v
      allocate (domain    (im1,jn1,nl1))
      domain(:,:,:)  = ypert%pt
      allocate (psdomain    (im1,jn1,nl1))
      psdomain(:,:,:)=  ypert%pt

      call prognostics_final (ypert)

!   Compute kdomain to maintain as much code history as possible
      call comp_kdomain (im1, jn1, nl1, domain, kdomain)
      k = kdomain(1)
      psdomain(:,:,1) =  domain(:,:,k)
 
!   Compute weighting terms
      allocate (jweights(jn1,2))
      allocate (kweights(nl1))
      allocate (glats(jn1,2))
      allocate (pmid(nl1))
      allocate (delbk(nl1))
      allocate (ak(nl1+1))
      allocate (bk(nl1+1))
      call horiz_grid (im1, jn1, jweights, glats)
      call vert_grid  (nl1, ks, ak, bk, kweights, ptop, pmid, delbk)

      fnames = fnames_def                          ! set default field names
      call find_name (nfields3d,fnames,'vpT',ndt)
      call getparam ( 'tname', tname_orig )
      fnames(ndt) = tname_orig
      call find_name (nfields3d,fnames,'dp_',ndp)
      call find_name (nfields,  fnames,'ps_',nps)
      nps=nps-nfields3d      
!      
! Loop over time
!
      allocate(stats3d(nlp1,measures,nfields3d,nOutfiles,ntimes))
      allocate(stats2d(    1,measures,nfields2d,nOutfiles,ntimes))
      allocate(energy( nl1+2,3,nOutfiles,ntimes))
      allocate(correls3d(nlp1,nfields3d ,ntimes))
      allocate(correls2d(  1,nfields2d ,ntimes))

      countKT = 0
      do ntimex=1,ntimes
!      
! Loop over kinds of data
!

        allocate(fields3d(im1,jn1,nl1,nfields3d,nOutfiles))
        allocate(fields3d_LPO(im1,jn1,nl1,nfields3d))
        allocate(fields2d(im1,jn1,  1,nfields2d,nOutfiles))

        call getparam ( 'tname', tname_orig )
        do kindx=1,nOutfiles
          countKT = countKT + 1 

          if ((knames(kindx)(1:3) == 'TLM') .or.      &
              (knames(kindx)(1:3) == 'ADM')     ) then
            nkind0=2
          else
            nkind0=1
          endif
          nkind1=kindx+nkind0-1

! Read data from TLM, ADM, or NLM data set for one time     
! Note that if nkind0=2, 2 kinds of data are read
          if (knames(kindx)(1:1) .ne. 'D' ) then
          if (knames(kindx)(1:1) .ne. 'R' ) then
            fnames(ndt) = tname_orig   ! reset temperature to what's in file
            print *, 'Temperature field on input set to: ', fnames(ndt)

            call read_data (im1,jn1,nl1,nfields,fnames, &
                            datesOut(ntimex),timesOut(ntimex),trim(jobIn(kindx)), &
                            fields3d(:,:,:,:,kindx:kindx),'nlm',ier)
            if (knames(kindx)(1:3) == 'TLM' .or. knames(kindx)(1:3) == 'ADM' ) then
            call read_data (im1,jn1,nl1,nfields,fnames, &
                            datesOut(ntimex),timesOut(ntimex),trim(jobIn(kindx+1)), &
                            fields3d(:,:,:,:,kindx+1:kindx+1),'nlm',ier)
            endif
          endif
          endif
          if (pert_with_LPO) fields3d_LPO(:,:,:,:) = fields3d(:,:,:,:,kindx)
          print *, '_RT fields3d 1 after lpo', sum(fields3d(:,:,:,:,kindx:kindx))
          print *, '_RT fields3d 2 after lpo', sum(fields3d(:,:,:,:,nkind1:nkind1))
!
! Transform T fields if requested 
! (In ADM calculations, mofify the delp field also)
! ***POTENTIAL BUG***
! Descriptions of transform_T_tl and transform_T_ad say that pert and ref 
! fields are transformed, yet a simple else would result in the transforms of ref fields a second time.
          if (knames(kindx)(1:3) == 'TLM' ) then
            print*
            print*,'Transforming TLM and REF temperature '
            call transform_T_tl (im1,jn1,nl1,ptop,nfields3d,fnames(1:nfields3d), &
                                 fields3d(:,:,:,:,kindx:kindx+1) )
          else if (knames(kindx)(1:3) == 'ADM' ) then
            print*
            print*,'Transforming ADM and REF vpT temperature'
            call transform_T_ad (im1,jn1,nl1,ptop,nfields3d,fnames(1:nfields3d), &
                                 fields3d(:,:,:,:,kindx:kindx+1))
          else if (knames(kindx)(1:3) == 'NLM' ) then    
            print*
            print*,'Transforming NLM temperature '
            call transform_T (im1,jn1,nl1,ptop,nfields3d,fnames(1:nfields3d),    &
                              fields3d(:,:,:,:,kindx))
          endif

          print *, '_RT fields3d 1 before ps', sum(fields3d(:,:,:,:,kindx:kindx))
          print *, '_RT fields3d 2 before ps', sum(fields3d(:,:,:,:,nkind1:nkind1))
!
! Compute ps fields
! TLM and NLM compute ps in same way
          if ((knames(kindx)(1:3) == 'TLM' ) .or.                          & 
              (knames(kindx)(1:3) == 'NLM' )      ) then
            call calculate_ps (im1, jn1, nl1, ptop, fields3d(:,:,:,ndp,kindx),    & 
                               fields2d(:,:,1,nps,kindx))
          endif
! reference fields always follow perturbed fields
          if ((knames(kindx)(1:3) == 'TLM' ) .or.                          & 
              (knames(kindx)(1:3) == 'ADM' )      ) then
            call calculate_ps (im1, jn1, nl1, ptop, fields3d(:,:,:,ndp,kindx+1),  & 
                               fields2d(:,:,1,nps,kindx+1))
          endif
! Adjoint ps calculation is different
          if (knames(kindx)(1:3) == 'ADM' ) then
            call calculate_ps_ad (im1, jn1, nl1, fields3d(:,:,:,ndp,kindx), & 
                                  delbk, fields2d(:,:,1,nps,kindx) )
          endif

          print *, '_RT fields3d 1 before diff', sum(fields3d(:,:,:,:,kindx:kindx))
          print *, '_RT fields3d 2 before diff', sum(fields3d(:,:,:,:,nkind1:nkind1))
!
!  Compute appropriate differences          
          if ((knames(kindx)(1:4) == 'DTT_') .or. (knames(kindx)(1:4) == 'DAA_') ) then
            fields3d(:,:,:,:,kindx)=fields3d(:,:,:,:,1)-fields3d(:,:,:,:,3)
            fields2d(:,:,:,:,kindx)=fields2d(:,:,:,:,1)-fields2d(:,:,:,:,3)
          else if ((knames(kindx)(1:4) == 'DNN_') ) then
            fields3d(:,:,:,:,kindx)=fields3d(:,:,:,:,kindx-2)-  &
                                    fields3d(:,:,:,:,kindx-1)
            fields2d(:,:,:,:,kindx)=fields2d(:,:,:,:,kindx-2)-  &
                                    fields2d(:,:,:,:,kindx-1)
          else if ((knames(kindx)(1:4) == 'DTN_') ) then
            fields3d(:,:,:,:,kindx)=fields3d(:,:,:,:,1)-  &
                                    fields3d(:,:,:,:,kindx-1)
            fields2d(:,:,:,:,kindx)=fields2d(:,:,:,:,1)-  &
                                    fields2d(:,:,:,:,kindx-1)
          endif
          print *, '_RT fields3d 1 before ad', sum(fields3d(:,:,:,:,kindx:kindx))
          print *, '_RT fields3d 2 before ad', sum(fields3d(:,:,:,:,nkind1:nkind1))

          if (knames(kindx)(1:3) == 'ADM' .and. (rescale_adj)                &
              .and. nOutfiles .eq. 2 ) then

                call scale_adj_ (im1,jn1,nl1,nfields2d,nfields3d,fnames,     &
                                jweights,kweights,fields2d(:,:,:,:,kindx),   &
                                fields3d(:,:,:,:,kindx))
          endif

!     
! compute horizontal statistics
          if (knames(kindx)(1:1) /= 'R' ) then
          print *, '_RT fields3d 1', sum(fields3d(:,:,:,:,kindx:kindx))
          print *, '_RT fields3d 2', sum(fields3d(:,:,:,:,nkind1:nkind1))
            call horiz_stats_ (im1, jn1, nl1, nfields3d, measures,     &
                              nkind0, nlp1,                            &
                              domain, udomain, vdomain, kdomain,       &
                              jweights, kweights,                      &
                              fields3d(:,:,:,:,kindx:nkind1),          &
                              fnames(1:nfields3d), mnames(:),          &
                              stats3d(:,:,:,kindx:nkind1,ntimex) )
            print *, '_RT sum of stats 3d', sum(stats3d(:,:,:,kindx:nkind1,ntimex))
            call horiz_stats_ (im1, jn1, 1, nfields2d, measures,       &
                              nkind0, 1,                               &
                              psdomain, udomain, vdomain, kdomainps,   &
                              jweights, kweights,                      &
                              fields2d(:,:,:,:,kindx:nkind1),          &
                              fnames(nfields3d+1:nfields), mnames(:),  & 
                              stats2d(:,:,:,kindx:nkind1,ntimex) ) 
            print *, '_RT sum of stats 2d', sum(stats2d(:,:,:,kindx:nkind1,ntimex))
          endif  

! GFIO version
          call write_data_ (im1,jn1,nl1,lm1,nps,nfields2d,nfields3d,         &
                           ptop,ks,ak,bk,                                    &
                           kdomain,psdomain,udomain,vdomain,domain,          &
                           fields2d(:,:,1,:,kindx),fields3d(:,:,:,:,kindx),  &
                           trim(knames(kindx)),outNames(kindx),fnames,         &
                           timesOut(ntimex),datesOut(ntimex))

          if ( (pert_with_LPO).and.(knames(kindx)(1:3) == 'TLM') ) &
          call write_pert_data_ (im1,jn1,nl1,lm1,nps,nfields2d,nfields3d,         &
                           ptop,ks,ak,bk,                                    &
                           kdomain,psdomain,udomain,vdomain,domain,          &
                           fields2d(:,:,1,:,kindx+1),fields3d_LPO(:,:,:,:),  &
                           trim(knames(kindx)),'pert_with_LPO.nc4',fnames,         &
                           timesOut(ntimex),datesOut(ntimex))

        enddo ! loop over kinds
       
!
! Compute E-norm if desired
!
        if (e_calc) then 
!
!  Only can compute E if T field exists 
!
          if (fnames(ndT) == 'T__') then
            do kindx=1,nOutfiles      
              kdp=0
              if (knames(kindx) (1:3) == 'TLM') kdp=kindx+1
              if (knames(kindx) (1:3) == 'DNN') kdp=kindx-1
              if (knames(kindx) (1:3) == 'DTT') kdp=kindx-1
              if (knames(kindx) (1:3) == 'DTN') kdp=kindx-2
              if (kdp /= 0) then
                call compute_e_ (im1,jn1,nl1,nfields3d,        &
                   udomain, vdomain, domain,                        &
                   kdomain, jweights(:,2),                          &
                   fields3d(:,:,:,:,kindx),fields3d(:,:,:,ndp,kdp), &
                   fields2d(:,:,1,nps,kindx),fnames(1:nfields3d),   &
                   energy(:,:,kindx,ntimex))
              endif  
            enddo
          else
            print *,' '
            print *,' ******  ERROR DETECTED *******'
            print *,' E-norm could not be calculated'
            print *,'fnames(ndT)= not present'
            print *,'*******************************'
          endif
        endif
!
! Compute correlations of specified pair of kinds if desired
!
        if ((corr_kind(1) /= 0 ) .and. (corr_kind(2) /= 0 ) ) then
          call compute_corr_ (im1, jn1, nl1, nfields3d, nlp1,       &
                             udomain, vdomain, domain, kdomain,     &
                             jweights, kweights,                    &
                             fields3d(:,:,:,:, corr_kind(1)),       &
                             fields3d(:,:,:,:, corr_kind(2)),       &
                             fnames, correls3d(:,:,ntimex) )
     
          call compute_corr_ (im1, jn1, 1, nfields2d, 1,            &
                             psdomain, psdomain, psdomain,          &
                             kdomainps,                             &
                             jweights, kweights,                    &
                             fields2d(:,:,:,:, corr_kind(1)),       &
                             fields2d(:,:,:,:, corr_kind(2)),       &
                             fnames, correls2d(:,:,ntimex) )
        endif     

        deallocate(fields3d)
        deallocate(fields2d)

      enddo   ! loop over times

!
!
! Print out all tables previously computed
! Print stats table
      call print_stats_ (im1, jn1, nl1, nfields2d, nfields3d, measures, nlp1, nOutfiles, &
                         ntimes, fnames, mnames,                                          &
                         domain, kdomain, timesOut, pmid, stats2d, stats3d )
!
! Print E-norm table
      if (e_calc) then
        call print_e_ (im1,jn1,nl1,nOutfiles,ntimes,     &
                      domain,kdomain,timesOut,pmid,energy) 
      endif
!
! Print correlation table
      if ((corr_kind(1) /= 0 ) .and. (corr_kind(2) /= 0 ) ) then
        call print_correls_ (im1,jn1,nl1,nfields2d,nfields3d,nlp1,ntimes,fnames(:), &
                            trim(knames(corr_kind(1))),trim(knames(corr_kind(2))),  &
                            domain,kdomain,timesOut,pmid,                           &
                            correls2d,correls3d )
      endif

      print *,' '
      print *,' '
      print *,' '
      print *,' PROGRAM COMPLETED'

      deallocate(domain,udomain,vdomain,psdomain)
      deallocate(energy)
      deallocate(stats3d)
      deallocate(stats2d)
      deallocate(correls3d)
      deallocate(correls2d)
      deallocate(jweights)
      deallocate(kweights)
      deallocate(glats)
      deallocate(pmid)
      deallocate(delbk)
      deallocate(ak)
      deallocate(bk)

      end subroutine field_stats_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: write_data --- write out arrays in stats format to GFIO format
!
!
! !INTERFACE:
!
      Subroutine write_data_ (im1,jn1,nl1,lm1,nps,nfields2d,nfields3d,ptop,ks, &
                             ak,bk, kdomain,psdomain,udomain,vdomain,domain,   &
                             fields2d,fields3d,kname,job1,fnames,              &
                             nhms,nymd)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: im1, jn1, nl1, lm1, nps, nfields2d, nfields3d
      real(r8), intent(in) :: ptop
      integer,  intent(in) :: ks
      real(r8), intent(in) :: ak(nl1+1)
      real(r8), intent(in) :: bk(nl1+1)
      integer, intent(in)  :: kdomain(2)
      real(r8), intent(in) :: psdomain(im1,jn1,nl1)
      real(r8), intent(in) :: udomain(im1,jn1,nl1)
      real(r8), intent(in) :: vdomain(im1,jn1,nl1)
      real(r8), intent(in) :: domain(im1,jn1,nl1)
      real(r8), intent(in) :: fields2d(im1,jn1, 1,nfields2d)
      real(r8), intent(in) :: fields3d(im1,jn1,nl1,nfields3d)
      character(len=4), intent(in) :: kname
      character(len=*), intent(in) :: job1
      character(len=3), intent(in) :: fnames(nfields2d+nfields3d)
      integer,  intent(in) :: nhms, nymd

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
  
! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Oct2004  Winslow    Original code
!  07Apr2005  Elena N.   Changed output precision to 32-bit  
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables
                                                                                                                    
      integer :: ier
      character(len=18), parameter :: fileTitle  =                     &
                                             'post proc. fields '
      character(len=10), parameter :: fileSource = 'stats5.f90'
      logical, parameter :: funtag = .true.
      logical, parameter :: fverbose = .true.
      logical, parameter :: fnew = .true.
      integer, parameter :: frq = 060000
      integer :: dpid,u_id,v_id,T_id,q_id        ! field id
      integer :: ierr                            ! return code
      integer :: prec_das                        ! precision for DAS output
                                                 ! 0: 32 bits
                                                 ! 1: 64 bits (default)
      integer :: nAllFields
      integer :: i,j,k
      real(r8),         allocatable :: allFields(:,:,:,:)
      character(len=256), allocatable :: allFieldNames(:)
      character(len=256), allocatable :: allFieldTitles(:)
      character(len=256), allocatable :: allFieldUnits(:)
                  
                  
      nAllFields = (nfields2d+nfields3d)
      allocate (allFields(im1,jn1,nl1,nAllFields))
      allocate (allFieldNames(nAllFields))
      allocate (allFieldUnits(nAllFields))
      allocate (allFieldTitles(nAllFields))

!     Apply the LPO and transfer to a GFIO array
      allFields = 0.0 
      call find_name (nfields3d,fnames,'dp_',dpid )
      call find_name (nfields3d,fnames,'u__',u_id )
      call find_name (nfields3d,fnames,'v__',v_id )
      call find_name (nfields3d,fnames,'T__',T_id )
      call find_name (nfields3d,fnames,'q__',q_id )
      do k=kdomain(1),kdomain(2)
      do j=1,jn1
      do i=1,im1
          if ( domain(i,j,k) .ne. 0.0) allFields(i,j,nl1,1)=fields2d(i,j,1,nps )
          if ( domain(i,j,k) .ne. 0.0) allFields(i,j,k  ,2)=fields3d(i,j,k,dpid)
          if (udomain(i,j,k) .ne. 0.0) allFields(i,j,k  ,3)=fields3d(i,j,k,u_id)
          if (vdomain(i,j,k) .ne. 0.0) allFields(i,j,k  ,4)=fields3d(i,j,k,v_id)
          if ( domain(i,j,k) .ne. 0.0) allFields(i,j,k  ,5)=fields3d(i,j,k,T_id)
          if ( domain(i,j,k) .ne. 0.0) allFields(i,j,k  ,6)=fields3d(i,j,k,q_id)
      enddo
      enddo
      enddo

!     Set the variable names
      if (kname(1:3) == 'TLM' .or.  kname(1:2) == 'DT') then
        allFieldNames (1) = 'ps_tl'
        allFieldUnits (1) = 'Pa'
        allFieldTitles(1) = 'pert_ps'
        allFieldNames (2) = 'delp_tl'
        allFieldUnits (2) = 'Pa'
        allFieldTitles(2) = 'pert_delp'
        allFieldNames (3) = 'u_tl'
        allFieldUnits (3) = 'm/s'
        allFieldTitles(3) = 'pert_u'
        allFieldNames (4) = 'v_tl'
        allFieldUnits (4) = 'm/s'
        allFieldTitles(4) = 'pert_v'
        allFieldNames (5) = 'T_tl'
        allFieldUnits (5) = 'unknown'
        allFieldTitles(5) = 'pert_T'
        allFieldNames (6) = 'q_tl'
        allFieldUnits (6) = 'kg/kg'
        allFieldTitles(6) = 'pert_q'
      endif
      if (kname(1:3) == 'ADM' .or.  kname(1:2) == 'DA') then
        allFieldNames (1) = 'ps_ad'
        allFieldUnits (1) = 'Pa'
        allFieldTitles(1) = 'dJ/dps'
        allFieldNames (2) = 'delp_ad'
        allFieldUnits (2) = 'J/Pa'
        allFieldTitles(2) = 'dJ/ddelp'
        allFieldNames (3) = 'u_ad'
        allFieldUnits (3) = 'J*s/m'
        allFieldTitles(3) = 'dJ/du'
        allFieldNames (4) = 'v_ad'
        allFieldUnits (4) = 'J*s/m'
        allFieldTitles(4) = 'dJ/dv'
        allFieldNames (5) = 'T_ad'
        allFieldUnits (5) = 'unknown'
        allFieldTitles(5) = 'dJ/dT'
        allFieldNames (6) = 'q_ad'
        allFieldUnits (6) = 'J*kg/kg'
        allFieldTitles(6) = 'dJ/dq'
      endif
      if (kname(1:1) == 'R' .or. kname(1:1) == 'N' .or.  kname(1:2) == 'DN') then
        allFieldNames (1) = 'ps'
        allFieldUnits (1) = 'Pa'
        allFieldTitles(1) = 'Surface Pressure'
        allFieldNames (2) = 'delp'
        allFieldUnits (2) = 'Pa'
        allFieldTitles(2) = 'Pressure Thickness'
        allFieldNames (3) = 'u'
        allFieldUnits (3) = 'm/s'
        allFieldTitles(3) = 'Zonal Wind'
        allFieldNames (4) = 'v'
        allFieldUnits (4) = 'm/s'
        allFieldTitles(4) = 'Meridional Wind'
        allFieldNames (5) = 'T'
        allFieldUnits (5) = 'unknown'
        allFieldTitles(5) = 'Temperature'
        allFieldNames (6) = 'q'
        allFieldUnits (6) = 'kg/kg'
        allFieldTitles(6) = 'Specific Humidity'
      endif
                  
!     Write in GFIO format           
      prec_das = 0 
      call GFIO_PutFld( fileTitle,fileSource,trim(job1),                        &
                        nymd,   nhms,  frq,                                     &
                        im1,   jn1,    nl1,                                     &
                        ptop,     ks,        ak,  bk,                           &
                        nAllFields, allFieldNames,allFieldTitles,allFieldUnits, &
                        allFields, ierr,                                        &
                        iprec=prec_das,  untag=funtag, verbose=fverbose )
                  
      deallocate(allFields,allFieldNames,allFieldTitles,allFieldUnits)
                  
      end subroutine write_data_

!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: write_pert_data --- write out arrays in stats format to GFIO format
!
!
! !INTERFACE:
!
      Subroutine write_pert_data_ (im1,jn1,nl1,lm1,nps,nfields2d,nfields3d,ptop,ks, &
                             ak,bk, kdomain,psdomain,udomain,vdomain,domain,   &
                             fields2d,fields3d,kname,job1,fnames,              &
                             nhms,nymd)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: im1, jn1, nl1, lm1, nps, nfields2d, nfields3d
      real(r8), intent(in) :: ptop
      integer,  intent(in) :: ks
      real(r8), intent(in) :: ak(nl1+1)
      real(r8), intent(in) :: bk(nl1+1)
      integer, intent(in)  :: kdomain(2)
      real(r8), intent(in) :: psdomain(im1,jn1,nl1)
      real(r8), intent(in) :: udomain(im1,jn1,nl1)
      real(r8), intent(in) :: vdomain(im1,jn1,nl1)
      real(r8), intent(in) :: domain(im1,jn1,nl1)
      real(r8), intent(in) :: fields2d(im1,jn1, 1,nfields2d)
      real(r8), intent(in) :: fields3d(im1,jn1,nl1,nfields3d)
      character(len=4), intent(in) :: kname
      character(len=*), intent(in) :: job1
      character(len=3), intent(in) :: fnames(nfields2d+nfields3d)
      integer,  intent(in) :: nhms, nymd

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
  
! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Oct2004  Winslow    Original code
!  07Apr2005  Elena N.   Changed output precision to 32-bit  
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables
                                                                                                                    
      integer :: ier
      character(len=18), parameter :: fileTitle  =                     &
                                             'post proc. fields '
      character(len=10), parameter :: fileSource = 'stats5.f90'
      logical, parameter :: funtag = .true.
      logical, parameter :: fverbose = .true.
      logical, parameter :: fnew = .true.
      integer, parameter :: frq = 060000
      integer :: dpid,u_id,v_id,T_id,q_id        ! field id
      integer :: ierr                            ! return code
      integer :: prec_das                        ! precision for DAS output
                                                 ! 0: 32 bits
                                                 ! 1: 64 bits (default)
      integer :: nAllFields
      integer :: i,j,k
      real(r8),         allocatable ::  phis(:,:), sgh(:,:), Ts(:,:), oro(:,:), ps(:,:)
      real(r8),         allocatable :: delp_tl(:,:,:), u_tl(:,:,:), v_tl(:,:,:), pt_tl(:,:,:)
      real(r8),         allocatable :: q_tl(:,:,:,:)
      type (dyn_vect) :: w_f
      
      allocate( phis(im1,jn1), sgh(im1,jn1),  Ts(im1,jn1), oro(im1,jn1), ps(im1,jn1),stat=ierr)
        if(ierr/=0) call MP_die (myname,'Failed in alloc(phis,etc)',ierr)
      allocate( delp_tl(im1,jn1,nl1), u_tl(im1,jn1,nl1), v_tl(im1,jn1,nl1), pt_tl(im1,jn1,nl1), q_tl(im1,jn1,nl1,lm1),stat=ierr )
        if(ierr/=0) call MP_die (myname,'Failed in alloc(delp_ad,etc)',ierr)
!
! Initialize extraneous 2D fields of output vector to zeros
!
      phis = 0.d0
      sgh  = 0.d0
      Ts   = 0.d0
      oro  = 0.d0
      ps   = 0.d0
      delp_tl = 0.d0
         u_tl = 0.d0
         v_tl = 0.d0
        pt_tl = 0.d0
         q_tl = 0.d0
                  
!     Apply the LPO and transfer to a GFIO array
      call find_name (nfields3d,fnames,'dp_',dpid )
      call find_name (nfields3d,fnames,'u__',u_id )
      call find_name (nfields3d,fnames,'v__',v_id )
      call find_name (nfields3d,fnames,'T__',T_id )
      if (T_id.le.0) then
         call find_name (nfields3d,fnames,'vpT',T_id )
      endif
      call find_name (nfields3d,fnames,'q__',q_id )
      do j=1,jn1
      do i=1,im1
          if ( domain(i,j,1) .ne. 0.0) ps(i,j)=fields2d(i,j,1,nps )
      do k=kdomain(1),kdomain(2)
          if ( domain(i,j,k) .ne. 0.0) delp_tl(i,j,k) =fields3d(i,j,k,dpid)
          if (udomain(i,j,k) .ne. 0.0) u_tl(i,j,k)    =fields3d(i,j,k,u_id)
          if (vdomain(i,j,k) .ne. 0.0) v_tl(i,j,k)    =fields3d(i,j,k,v_id)
          if ( domain(i,j,k) .ne. 0.0) pt_tl(i,j,k)   =fields3d(i,j,k,T_id)
          if ( domain(i,j,k) .ne. 0.0) q_tl(i,j,k,1)  =fields3d(i,j,k,q_id)
      enddo
      enddo
      enddo
!
! --- 01Nov2005 --- Adjoint runs with a single tracer FOR NOW
!
      if ( lm1 > 1 ) q_tl(:,:,:,2:lm1) = 0.d0
!
      call dyn_null ( w_f )
      call Dyn_Init ( im1, jn1, nl1, lm1, w_f, ierr, ptop, ks, ak, bk )
      w_f%phis      = phis
      w_f%hs_stdv   = sgh
      w_f%Ts        = Ts
      w_f%lwi       = oro
      w_f%ps        = ps
      w_f%delp      = delp_tl
      w_f%u         = u_tl
      w_f%v         = v_tl
      w_f%pt        = pt_tl
      w_f%q         = q_tl

      if ( ierr .ne. 0 ) then
         call MP_die (myname,'Error from Dyn_Init',ierr)
      end if
      prec_das = 0
      call dyn_put ( trim(job1), nymd, nhms, prec_das, w_f, ierr, &
                     new=fnew, freq=frq )
        if(ierr/=0) then
          call MP_die (myname,'Error from dyn_put',ierr)
        else
          print*,'Wrote out perturbation to ',trim(job1)
        endif
      call dyn_null ( w_f )
      deallocate( delp_tl, u_tl, v_tl, pt_tl, q_tl )
      deallocate( phis, sgh,  Ts, oro, ps)

      end subroutine write_pert_data_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: compute_e --- Compute contributions to the total energy norm E.
!
!
! !INTERFACE:
!
      subroutine compute_e_ (im1,jn1,nl1,nfields,    &
                   udomain, vdomain, domain,             &
                   kdomain, jweights,                    &
                   fields,delp_ref,psfield,fnames,e) 

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1,nfields
      real(r8), intent(in) :: udomain(im1,jn1,nl1)
      real(r8), intent(in) :: vdomain(im1,jn1,nl1)
      real(r8), intent(in) :: domain(im1,jn1,nl1)
      integer, intent(in) :: kdomain(2)
      character(len=*), intent(in) :: fnames(nfields) 
      real(r8), intent(in)  :: jweights(jn1)
      real(r8), intent(in)  :: fields(im1,jn1,nl1,nfields)
      real(r8), intent(in)  :: delp_ref(im1,jn1,nl1)
      real(r8), intent(in)  :: psfield(im1,jn1)   

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: e(nl1+2,3)   

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
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!  08Apr2005  Elena N.   Corrected a bug with indexing in KEV calculation
!  08Apr2005  R. Errico  Corrected assumption that sum(jweights)=1

!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer  :: i,j,k
      integer  :: idu,idv,idt
      real(r8) :: dbleimr
      real(r8) :: tfac, ufac, pfac
      real(r8) :: dsigma(im1,jn1,nl1)
      real(r8) :: ps_ref(im1,jn1)
      real(r8) :: jfac(jn1)      
      real(r8) :: fa(im1,jn1) 
      real(r8) :: fb(im1,jn1) 
      real(r8) :: eu,ev,et,ep

      integer  :: ks
      real(r8)   :: ak(nl1+1)
      real(r8)   :: bk(nl1+1)
      real(r8)   :: pint
      real(r8)   :: ptop
        
      real(r8)   :: KEU, KEV                                                                                                             
      call set_eta(nl1, ks, ptop, pint, ak, bk)

      dbleimr=dble(im1)

      do j=1,jn1
        jfac(j)=jweights(j)/dbleimr
      enddo



! One factor of 0.5 is for defenition of energy as 1/2 * ...
! Assumption here is that sum(jweights)=1
      tfac=0.5d0*Cp/tref
      ufac=0.5d0
      pfac=0.5d0*R*tref/pref**2

      call find_name (nfields,fnames,'u__',idu)      
      call find_name (nfields,fnames,'v__',idv)      
      call find_name (nfields,fnames,'T__',idt)      

!
! compute 3-D varying dsigma

      ps_ref(:,:) = sum(delp_ref(:,:,:),3)       

      do k=1,nl1
        dsigma(:,:,k)=delp_ref(:,:,k)/ps_ref(:,:)
      enddo    
!
! Initialize e to zero
!
      e(:,:)=0.d0
!
!  Loop over levels
!
      do k=kdomain(1),kdomain(2)
!     
!  Calculations for u field
!     
        fa(:,:)=0.d0
        fb(:,:)=0.d0
        do j=1,jn1
          do i=1,im1
            if (udomain(i,j,k) .ne. 0.0) fa(i,j)=fields(i,j,k,idu)   
          enddo
        enddo


        do i=1,im1
          fb(i,  1)=fa(i,  2)*fa(i,  2)
          fb(i,jn1)=fa(i,jn1)*fa(i,jn1)
          do j=2,jn1-1     
            fb(i,j)=0.5*(fa(i,j)*fa(i,j) + fa(i,j+1)*fa(i,j+1))
          enddo
        enddo

        eu=0.d0
        do j=1,jn1
          do i=1,im1
            eu=eu+dsigma(i,j,k)*jfac(j)*fb(i,j)
          enddo
        enddo
!     
!  Calculations for v field
!    
        fa(:,:)=0.d0
        fb(:,:)=0.d0
        do j=1,jn1
          do i=1,im1
            if (vdomain(i,j,k) .ne. 0.0) fa(i,j)=fields(i,j,k,idv)     
          enddo
        enddo
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

        ev=0.d0
        do j=1,jn1
          do i=1,im1
            ev=ev+dsigma(i,j,k)*jfac(j)*fb(i,j)   
          enddo
        enddo

        e(k,1)=ufac*(eu+ev)   

        e(nl1+1,1)=e(nl1+1,1)+e(k,1)
!     
!  Calculations for T field
!     
        fa(:,:)=0.d0
        do j=1,jn1
          do i=1,im1
            if (domain(i,j,k) .ne. 0.0) fa(i,j)=fields(i,j,k,idt) 
          enddo
        enddo

        et=0.d0 
        do j=1,jn1
          do i=1,im1
            et=et+dsigma(i,j,k)*jfac(j)*fa(i,j)*fa(i,j)
          enddo
        enddo

        e(k,2)=tfac*et
        e(k,3)=e(k,1)+e(k,2)
        e(nl1+1,2)=e(nl1+1,2)+e(k,2)
!
!
      enddo ! loop over k
!
! Calculations for ps field
!
      fa(:,:)=0.d0
      do j=1,jn1
        do i=1,im1
          if (domain(i,j,kdomain(1)) .ne. 0.0) fa(i,j)=psfield(i,j)
        enddo
      enddo

      ep=0.d0
      do j=1,jn1
        do i=1,im1
          ep=ep+jfac(j)*fa(i,j)*fa(i,j)
        enddo
      enddo

      e(nl1+1,3)=e(nl1+1,1)+e(nl1+1,2)
      e(nl1+2,1)=0.d0                            ! not used
      e(nl1+2,2)=pfac*ep
      e(nl1+2,3)=e(nl1+2,2)+e(nl1+1,3)             ! total E

      end subroutine compute_e_ 
!
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
      subroutine print_e_ (im1,jn1,nl,nkinds,ntimes,      &
                          domain,kdomain,times,pmid,energy) 

!USES:
      use m_ioutil, only : luavail

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl,nkinds,ntimes
      real(r8), intent(in)  :: domain(im1,jn1,nl)
      integer, intent(in)  :: kdomain(2)
      integer, intent(in)  :: times(ntimes)
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

      integer :: i, j, k, kind, nt, imin, imax, jmin, jmax
      integer :: lunit1
      real(r4) :: aux_stats3d(nl), pe(nl+1) 
      integer  :: ks
      real(r8)   :: ak(nl+1)
      real(r8)   :: bk(nl+1)
      real(r8)   :: pint
      real(r8)   :: ptop

      call set_eta(nl, ks, ptop, pint, ak, bk)
     
      do i = 1, nl + 1
         pe(i) = ak(i) + bk(i)*pref
      enddo

      imin = im1
      imax = 1
      jmin = jn1
      jmax = 1
      k = kdomain(1)
      do i=1,im1
        do j=1,jn1
          if( domain(i,j,k) .ne. 0.0) then
            if (i .ge. imax) imax = i
            if (i .le. imin) imin = i
            if (j .ge. jmax) jmax = j
            if (j .le. jmin) jmin = j
          endif
        enddo
      enddo

      open(14,file='stats_ecalc.txt')

      print *,' '
      write(14,*)
      print *,' '
      write(14,*)
      print *,' '
      write(14,*)
      print ('(30x,a)'),'E - N O R M    S T A T I S T I C S'
      write(14,*) '                             E - N O R M    S T A T I S T I C S'
      print ('(4(a,I6),2(a,I6))')                           &
               ,'Longitudes: i = ',imin,' to ',imax,  &
               ', Latitudes: j =  ',jmin,' to ',jmax, &
               ', pressures: k =  ',kdomain(1),' to ',kdomain(2) 
      write(14,100) imin, imax, jmin, jmax, kdomain(1), kdomain(2) 
        
      do nt=1,ntimes
        do kind=1,nkinds

! Energy has only been computed for particular kinds
          if ( (knames(kind) (1:3) == 'TLM') .or.        &
               (knames(kind) (1:3) == 'DNN') .or.        &
               (knames(kind) (1:3) == 'DTT') .or.        &
               (knames(kind) (1:3) == 'DTN')      ) then
            print *,' '
            write(14,*)
            print *,' KIND=',trim(knames(kind)),'  TIME=',times(nt) 
            write(14,*) ' KIND=',trim(knames(kind)),'  TIME=',times(nt) 
            print ('(a,5x,a,8x,a,8x,a)'),'  k p-level','  KE',' APE','TOTE'  
            write(14,*) '  k p-level       KE         APE        TOTE'
            do k=kdomain(1),kdomain(2)
              print ('(i3,f8.1,1p3e12.2)'),k,pmid(k),energy(k,1:3,kind,nt)
              write(14,101) k,pmid(k),energy(k,1:3,kind,nt)
            enddo
            print ('(a,1p3e12.2)'),'      sums:',energy(nl+1,1:3,kind,nt)
            write(14,102) '      sums:',energy(nl+1,1:3,kind,nt)
            print ('(a,1p3e12.2)'),'0, psE, TOT',energy(nl+2,1:3,kind,nt)
            write(14,103) '0, psE, TOT',energy(nl+2,1:3,kind,nt)
          endif

        enddo
      enddo
      close(14)

!
      lunit1 = luavail()
!
! --- Write out binary data
!
      open ( lunit1, file = 'energy_stats.dat', form = 'unformatted',  &
             access = 'direct', recl = nl)

      k = 1
      do nt=1,ntimes
         do kind=1,nkinds
            aux_stats3d(1:nl) = 0.01*pmid(1:nl)
            write(lunit1, rec=k) aux_stats3d
            k = k + 1
            aux_stats3d(1:nl) = log10(0.01*pmid(1:nl))
            write(lunit1, rec=k) aux_stats3d
            k = k + 1            
            do i = 1, nl
                aux_stats3d(i) = log10(pe(i+1)/pe(i))
            enddo
            write(lunit1, rec=k) aux_stats3d
            k = k + 1
            aux_stats3d(1:nl) = pe(2:nl+1)
            write(lunit1, rec=k) aux_stats3d
            k = k + 1
            do j = 1, 3
               aux_stats3d(1) = energy(1, j, kind, nt)
               do i = 2, nl
                   aux_stats3d(i) = aux_stats3d(i-1) + energy(i, j, kind, nt)
               enddo
               write(lunit1, rec=k) aux_stats3d
               k = k + 1
            enddo
            do i = 1, 3
                  aux_stats3d(1:nl) = energy(1:nl, i, kind, nt)
                  write(lunit1, rec=k) aux_stats3d
                  k = k + 1
            enddo
         enddo
      enddo

      close(lunit1)
!
! --- Finish writing data

      100 format('Longitudes: i = ',i6,' to ',i6,', Latitudes: j =  ',i6, &
                 ' to ',i6,', pressures: k =  ',i6,' to ',i6)
      101 format(i3,f8.1,1p3e12.2)
      102 format(a,1p3e12.2)
      103 format(a,1p3e12.2)

      end subroutine print_e_

!
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: print_stats --- Print table of field statistics
!
!
! !INTERFACE:
!
      subroutine print_stats_ (im1,jn1,nl,nfields2d,nfields3d,measures,nlp1,nkinds, &
                              ntimes,fnames,mnames,                                 &
                              domain,kdomain,times,pmid,stats2d,stats3d )

!USES:

      use m_ioutil, only : luavail

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: im1,jn1,nl,nfields3d,nfields2d,measures
      integer,  intent(in)  :: nlp1,nkinds,ntimes
      character(len=*), intent(in) :: fnames(nfields3d+nfields2d) 
      character(len=*), intent(in) :: mnames(measures)
      real(r8), intent(in)  :: domain(im1,jn1,nl)
      integer, intent(in)  :: kdomain(2)
      integer,  intent(in)  :: times(ntimes)
      real(r8), intent(in)  :: pmid(nl)
      real(r8), intent(in)  :: stats3d(nlp1,measures,nfields3d,nkinds,ntimes)   
      real(r8), intent(in)  :: stats2d(   1,measures,nfields2d,nkinds,ntimes)   

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!
!  Print table of field statistics as a function of
!  data kind, data field, data time, and model pressure layer
!  Print tables for both 3-D and 2_D fields.
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

      integer :: i,j,k, nf, kind, m, nt, imin, imax, jmin, jmax
      integer :: lunit1
      real(r4) :: aux_stats3d(nl)

      imin = im1
      imax = 1
      jmin = jn1
      jmax = 1
      k = kdomain(1)
      do i=1,im1
        do j=1,jn1
          if( domain(i,j,k ) .ne. 0.0) then
            if (i .ge. imax) imax = i
            if (i .le. imin) imin = i
            if (j .ge. jmax) jmax = j
            if (j .le. jmin) jmin = j
          endif
        enddo
      enddo
!
      open(14,file='stats_stats.txt')
      print *,' '
      write(14,*)
      print *,' '
      write(14,*)
      print ('(30x,a)'),'F I E L D   S T A T I S T I C S'
      write(14,*) '                             F I E L D   S T A T I S T I C S'
!     print ('(4(a,f6.1),2(a,f9.1))')                           &
      print ('(4(a,I6),2(a,I6))')                           &
               ,'Longitudes: i = ',imin,' to ',imax,  &
               ', Latitudes: j =  ',jmin,' to ',jmax, &
               ', pressures: k =  ',kdomain(1),' to ',kdomain(2) 
      write(14,300) imin, imax, jmin, jmax, kdomain(1), kdomain(2)
        
      do nt=1,ntimes
        print *,' '
        write(14,*)
        print *,' TIME=',times(nt) 
        write(14,*) ' TIME=',times(nt) 
        do kind=1,nkinds
          print *,' '
          write(14,*)
          print *,' KIND=',trim(knames(kind))
          write(14,*) ' KIND=',trim(knames(kind))
          do nf=1,nfields3d
            print *,' '
            write(14,*)
            print ('(5a,i6)'),' Field=',fnames(nf),                      &
                                 '  kind=',trim(knames(kind)),'  time=',times(nt)
            write(14,301) ' Field=',fnames(nf),'  kind=',trim(knames(kind)),'  time=',times(nt)
            print ('(a,6a12)'),'  k p-level',mnames(1:measures)
            write(14,302) '  k p-level',mnames(1:measures)
            do k=kdomain(1),kdomain(2)
              print ('(i3,f8.1,1p6e12.2)'),                               &
                    k,pmid(k),stats3d(k,1:measures,nf,kind,nt)
              write(14,303) k,pmid(k),stats3d(k,1:measures,nf,kind,nt)
            enddo
            print ('(a,1p6e12.2)'),                                       &
                    ' v-integral',stats3d(k,1:measures,nf,kind,nt)
            write(14,304) ' v-integral',stats3d(k,1:measures,nf,kind,nt)
          enddo
          print *,' '
          write(14,*)
          print ('(a11,6a12)'),'  2-D field',mnames(1:measures)
          write(14,305) '  2-D field',mnames(1:measures)
          do nf=1,nfields2d
            print ('(8x,a3,1p6e12.2)'),                            &
                    fnames(nfields3d+nf),stats2d(1,1:measures,nf,kind,nt)
            write(14,306) fnames(nfields3d+nf),stats2d(1,1:measures,nf,kind,nt)
          enddo
!    
        enddo  ! loop over kinds
      enddo    ! loop over times      
      close(14)

!
      lunit1 = luavail()
!
! --- Write out binary data
!
      open ( lunit1, file = 'stats3d.dat', form = 'unformatted',  &
             access = 'direct', recl = nl)

      k = 1
      do nt=1,ntimes
         do kind=1,nkinds
            do nf=1,nfields3d
               aux_stats3d(1:nl) = max(abs(stats3d(1:nl, 1, nf, kind, nt)),abs(stats3d(1:nl, 2, nf, kind, nt)))
               write(lunit1, rec=k) aux_stats3d
               k = k + 1
               do m = 3, measures, 2
                  aux_stats3d(1:nl) = stats3d(1:nl, m, nf, kind, nt)
                  write(lunit1, rec=k) aux_stats3d
                  k = k + 1
               enddo
            enddo
         enddo
      enddo

      close(lunit1)
!
! --- Finish writing data
!
      300 format('Longitudes: i = ',i6,' to ',i6,', Latitudes: j =  ',i6, &
                 ' to ',i6,', pressures: k =  ',i6,' to ',i6)
      301 format(5a,i6)
      302 format(a,6a12)
      303 format(i3,f8.1,1p6e12.2)
      304 format(a,1p6e12.2)
      305 format(a11,6a12)
      306 format(8x,a3,1p6e12.2)


      end subroutine print_stats_

!
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: print_correls --- Print table of corresponding field correlations
!
!
! !INTERFACE:
!
      subroutine print_correls_ (im1,jn1,nl,nfields2d,nfields3d,nlp1,    &
                                ntimes,fnames,kname1,kname2,    &
                                domain,kdomain,times,pmid,    &
                                correls2d,correls3d )

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: im1,jn1,nl,nfields3d,nfields2d,nlp1,ntimes
      character(len=*), intent(in) :: fnames(nfields3d+nfields2d) 
      character(len=*), intent(in) :: kname1,kname2
      real(r8), intent(in)  :: domain(im1,jn1,nl)
      integer, intent(in)  :: kdomain(2)
      integer , intent(in)  :: times(ntimes)
      real(r8), intent(in)  :: pmid(nl)
      real(r8), intent(in)  :: correls3d(nlp1,nfields3d,ntimes) 
      real(r8), intent(in)  :: correls2d(   1,nfields2d,ntimes)

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!
!  Print table of corresponding field correlations afor 2 kinds of data,
!  as a function of field, time, and model pressure layer.
!  Print tables for both 3-D and 2_D fields.
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
 
      integer :: i ,j, k, nf, nt,imin,imax,jmin,jmax

      imin = im1
      imax = 1
      jmin = jn1
      jmax = 1
      k = kdomain(1)
      do i=1,im1
        do j=1,jn1
          if(domain(i,j,k) .ne. 0.0) then
            if (i .ge. imax) imax = i
            if (i .le. imin) imin = i
            if (j .ge. jmax) jmax = j
            if (j .le. jmin) jmin = j
          endif
        enddo
      enddo

      open(14,file='stats_correls.txt')
      print *,' '
      write(14,*)
      print *,' '
      write(14,*)
      print ('(30x,a)'),'C O R R E L A T I O N S'
      write(14,*) '                             C O R R E L A T I O N S'
      print *,' '
      write(14,*)
      print *,' For kinds = ',trim(kname1),' and ',trim(kname2)
      write(14,*) ' For kinds = ',trim(kname1),' and ',trim(kname2)
!     print ('(4(a,f6.1),2(a,f9.1))')                           &
      print ('(4(a,I6),2(a,I6))')                           &
               ,'Longitudes: i = ',imin,' to ',imax,  &
               ', Latitudes: j =  ',jmin,' to ',jmax, &
               ', pressures: k =  ',kdomain(1),' to ',kdomain(2) 
      write(14,200) imin, imax, jmin, jmax, kdomain(1), kdomain(2)

  
      do nt=1,ntimes
        print *,' '
        write(14,*)
        print *,' TIME = ',times(nt) 
        write(14,*) ' TIME = ',times(nt) 
        print ('(a,8(5x,a))'),'  k p-level',fnames(1:nfields3d)
        write(14,201) '  k p-level',fnames(1:nfields3d)
        do k=kdomain(1),kdomain(2)
          print ('(i3,f8.1,8f8.3)'),                         &
                    k,pmid(k),correls3d(k,1:nfields3d,nt)
          write(14,202) k,pmid(k),correls3d(k,1:nfields3d,nt)
        enddo
        print ('(a,8f8.3)'),'  3d corr  ',correls3d(nlp1,1:nfields3d,nt)
        write(14,203) '  3d corr  ',correls3d(nlp1,1:nfields3d,nt)

        if (nfields2d > 0) then
          print *,' '
          write(14,*) 
          print ('(a11,8(5x,a))'),' 2D fields ',             &
                    fnames(nfields3d+1:nfields3d+nfields2d)
          write(14,204) ' 2D fields ', fnames(nfields3d+1:nfields3d+nfields2d)
          print ('(11x,8f8.3)'),correls2d(1,1:nfields2d,nt)
          write(14,205) correls2d(1,1:nfields2d,nt)
        endif  
  
        print *,' '
        write(14,*)
        print *,' '
        write(14,*)

      enddo    ! loop over times      
      close(14)
      200 format('Longitudes: i = ',i6,' to ',i6,', Latitudes: j =  ',i6, &
                 ' to ',i6,', pressures: k =  ',i6,' to ',i6)
      201 format(a,8(5x,a))
      202 format(i3,f8.1,8f8.3)
      203 format(a,8f8.3)
      204 format(a11,8(5x,a))
      205 format(11x,8f8.3)

      end subroutine print_correls_


!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: horiz_stats --- Compute some statistics of 2-D or 3-D fields.
!
!
! !INTERFACE:
!
      subroutine horiz_stats_ (im1, jn1, nl1, nfields, measures, nx, nlp1,    &
                              domain, udomain, vdomain, kdomain,             &
                              jweights, kweights,                            &
                              fields, fnames, mnames, stats )

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1,nfields,measures,nx,nlp1
      real(r8), intent(in) :: domain(im1,jn1,nl1)
      real(r8), intent(in) :: udomain(im1,jn1,nl1)
      real(r8), intent(in) :: vdomain(im1,jn1,nl1)
      integer, intent(in) :: kdomain(2)
      real(r8), intent(in)  :: jweights(jn1,2)
      real(r8), intent(in)  :: kweights(nl1)
      real(r8), intent(in)  :: fields(im1,jn1,nl1,nfields,nx)    
      character(len=*), intent(in) :: fnames(nfields) 
      character(len=*), intent(in) :: mnames(measures)

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: stats(nlp1,measures,nfields,nx)   

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!
!  Compute some statistics of 2-D or 3-D fields.
!  For 3-D fields, stats for each model layer are determined separately, but
!  also full 3-D stats are determined.  The 2-D and 3-D single layer stats
!  use a proper latitudinal weighting, but no vertical (layer) weighting. The
!  full 3-D stats uses a weighting by an effective dsigma defined by the model
!  standard grid, with ps assumed to be 10**5 Pa. No vertical weighting is
!  applied to the max and min statistics.

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

      real(r8), parameter :: rbig=1.0d30
      integer  :: idmax,idmin,idmean,idrms,idstdv
      integer  :: i,j,k,nf,n
      integer  :: id1,id2
      integer  :: ifirst,ilast,jfirst,jlast
      integer  :: Idomain(im1,jn1,nl1)
      real(r8) :: kweight_sum
      real(r8) :: area
      real(r8) :: rmean
      real(r8) :: rsqr
      real(r8) :: rmean_sum
      real(r8) :: rsqr_sum
      real(r8) :: rmax      
      real(r8) :: rmin  
      real(r8) :: rms
      real(r8) :: rstdv

      call find_name (measures,mnames,' MAX',idmax)      
      call find_name (measures,mnames,' MIN',idmin)      
      call find_name (measures,mnames,'MEAN',idmean)      
      call find_name (measures,mnames,' RMS',idrms)      
      call find_name (measures,mnames,'STDV',idstdv)      
   
      kweight_sum=0.d0
      do k=kdomain(1),kdomain(2)
        kweight_sum=kweight_sum+kweights(k)
      enddo

      do n=1,nx
      do nf=1,nfields
        id1=2           ! pointer to jweights
        id2=3           ! pointer to ijdomain
        Idomain = domain
        if (fnames(nf)(1:1) == 'u' ) then
          Idomain = nint(udomain)
          id1=1      
          id2=1
        else if (fnames(nf)(1:1) == 'v' ) then 
          Idomain = nint(vdomain)
          id2=2
        endif     

        area=0.d0
        k = kdomain(1)
        do j=1,jn1
          do i=1,im1
             if (Idomain(i,j,k) .ne. 0.0 ) then
               area=area+jweights(j,id1)
               Idomain(i,j,k) = 1.0
             endif
          enddo
        enddo

        rmean_sum=0.d0
        rsqr_sum =0.d0
        stats(nl1+1,idmax,nf,n)=0.d0 !_RT -rbig
        stats(nl1+1,idmin,nf,n)=0.d0 !_RT rbig
 
        do k=kdomain(1),kdomain(2)

          rmax =-rbig
          rmin = rbig
          rmean= 0.d0
          rsqr = 0.d0
          do j=1,jn1
            do i=1,im1
              if (fields(i,j,k,nf,n) > rmax .and. Idomain(i,j,k) .ne. 0.0) rmax=fields(i,j,k,nf,n)
              if (fields(i,j,k,nf,n) < rmin .and. Idomain(i,j,k) .ne. 0.0) rmin=fields(i,j,k,nf,n)
              if (Idomain(i,j,k) .ne. 0.0) rmean=rmean+jweights(j,id1)*fields(i,j,k,nf,n)
              if (Idomain(i,j,k) .ne. 0.0) rsqr=rsqr+jweights(j,id1)*fields(i,j,k,nf,n)*fields(i,j,k,nf,n)
            enddo
          enddo

          if (area>0.d0) then
              rmean=rmean/area
              rsqr=rsqr/area 
          else
              print *, 'SEVERE WARNING: Suspect trouble in stats: area is zero'
          endif

          if ( rsqr > 0.d0 ) then
            rms=sqrt(rsqr)
          else
            rms=0.d0
          endif 
          rstdv=rsqr-rmean*rmean
          if ( rstdv > 0.d0 ) then
            rstdv=sqrt(rstdv)
          else
            rstdv=0.d0
          endif

          stats(k,idmax ,nf,n)=rmax   
          stats(k,idmin ,nf,n)=rmin   
          stats(k,idmean,nf,n)=rmean
          stats(k,idrms ,nf,n)=rms   
          stats(k,idstdv,nf,n)=rstdv

          if (nlp1 > 1) then
            if ( stats(nlp1,idmax,nf,n) < stats(k,idmax,nf,n) ) &
                 stats(nlp1,idmax,nf,n) = stats(k,idmax,nf,n)  
            if ( stats(nlp1,idmin,nf,n) > stats(k,idmin,nf,n) ) &
                 stats(nlp1,idmin,nf,n) = stats(k,idmin,nf,n)  
            rmean_sum=rmean_sum+kweights(k)*rmean
            rsqr_sum=rsqr_sum+kweights(k)*rsqr
          endif

        enddo  ! loop over levels

        if (nlp1 > 1) then
          rmean_sum=rmean_sum/kweight_sum
          rsqr_sum =rsqr_sum /kweight_sum

          stats(nlp1,idmean,nf,n)=rmean_sum
          if ( rsqr_sum > 0.d0 ) then
            stats(nlp1,idrms,nf,n)=sqrt(rsqr_sum)
          else
            stats(nlp1,idrms,nf,n)=0.d0
          endif 
          rstdv=rsqr_sum-rmean_sum*rmean_sum
          if ( rstdv > 0.d0 ) then
            stats(nlp1,idstdv,nf,n)=sqrt(rstdv)
          else
            stats(nlp1,idstdv,nf,n)=0.d0
          endif
        endif

     enddo     ! loop over fields
     enddo     ! loop over nx sets of fields

     end subroutine horiz_stats_      


!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: compute_corr --- Compute horizontal correlations of 
!                            corresponding 2-D or 3-D fields
!
!
! !INTERFACE:
!
      subroutine compute_corr_ (im1, jn1, nl1, nfields, nlp1,         &
                              udomain, vdomain, domain, kdomain,     &
                              jweights, kweights,                    &
                              fields1, fields2, fnames, correls )

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im1,jn1,nl1,nfields,nlp1
      real(r8), intent(in) :: udomain(im1,jn1,nl1)
      real(r8), intent(in) :: vdomain(im1,jn1,nl1)
      real(r8), intent(in) :: domain(im1,jn1,nl1)
      integer, intent(in) :: kdomain(2)
      character(len=*), intent(in) :: fnames(nfields) 
      real(r8), intent(in)  :: jweights(jn1,2)
      real(r8), intent(in)  :: kweights(nl1)
      real(r8), intent(in)  :: fields1(im1,jn1,nl1,nfields)    
      real(r8), intent(in)  :: fields2(im1,jn1,nl1,nfields)    

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: correls(nlp1,nfields)   

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!
!  Compute horizontal correlations of corresponding 2-D or 3-D fields.
!  For 3-D fields, correls for each model layer are determined separately, but
!  also full 3-D correls are determined.  The 2-D and 3-D single layer stats
!  use a proper latitudinal weighting, but no vertical (layer) weighting. The
!  full 3-D stats uses a weighting by an effective dsigma defined by the model
!  standard grid, with ps assumed to be 10**5 Pa.

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

      integer  :: i,j,k,nf
      integer  :: id1,id2
      integer  :: ifirst,ilast,jfirst,jlast
      integer  :: Idomain(im1,jn1,nl1)
      real(r8) :: kweight_sum
      real(r8) :: area
      real(r8) :: f1f1, f2f2
      real(r8) :: f1f1_3d,f2f2_3d
      real(r8) :: f1f2
      real(r8) :: f1f2_3d

      kweight_sum=0.d0
      do k=kdomain(1),kdomain(2)
        kweight_sum=kweight_sum+kweights(k)
      enddo

      do nf=1,nfields

        id1=2           ! pointer to jweights
        id2=3           ! pointer to ijdomain
        Idomain = domain
        if (fnames(nf)(1:1) == 'u' ) then
          id1=1      
          id2=1
          Idomain = udomain
        else if (fnames(nf)(1:1) == 'v' ) then 
          id2=2
          Idomain = vdomain
        endif     

        area=0.d0
        k = kdomain(1)
        do j=1,jn1
          do i=1,im1
            if (Idomain(i,j,k) .ne. 0.0) then
              area=area+jweights(j,id1)
            endif
          enddo
        enddo

        f1f2_3d=0.d0
        f1f1_3d=0.d0
        f2f2_3d=0.d0
 
        do k=kdomain(1),kdomain(2)

          f1f2=0.d0
          f1f1=0.d0
          f2f2=0.d0

          do j=1,jn1
            do i=1,im1
              if (Idomain(i,j,k) .ne. 0.0) then
                f1f2=f1f2+jweights(j,id1)*fields1(i,j,k,nf)*fields2(i,j,k,nf)
                f1f1=f1f1+jweights(j,id1)*fields1(i,j,k,nf)*fields1(i,j,k,nf)
                f2f2=f2f2+jweights(j,id1)*fields2(i,j,k,nf)*fields2(i,j,k,nf)
              endif
            enddo
          enddo
     
          f1f2=f1f2/area
          f1f1=f1f1/area
          f2f2=f2f2/area

          f1f2_3d=f1f2_3d+kweights(k)*f1f2
          f1f1_3d=f1f1_3d+kweights(k)*f1f1
          f2f2_3d=f2f2_3d+kweights(k)*f2f2

          if (f1f1 <= 0.d0 ) then
            correls(k,nf)=2.d0
            if (f2f2 <= 0.d0 ) then
              correls(k,nf)=0.d0
            endif
          else if (f2f2 <= 0.d0 ) then
             correls(k,nf)=3.d0
          else
             correls(k,nf)=f1f2/sqrt(f1f1*f2f2)
          endif 

        enddo  ! loop over levels

        if (nlp1 > 1) then

          f1f2_3d=f1f2_3d/kweight_sum
          f1f1_3d=f1f1_3d/kweight_sum
          f2f2_3d=f2f2_3d/kweight_sum

          if (f1f1_3d <= 0.d0 ) then
            correls(nlp1,nf)=2.d0
            if (f2f2_3d <= 0.d0 ) then
              correls(nlp1,nf)=0.d0
            endif
          else if (f2f2_3d <= 0.d0 ) then
             correls(nlp1,nf)=3.d0
          else
             correls(nlp1,nf)=f1f2_3d/sqrt(f1f1_3d*f2f2_3d)
          endif 
       endif

     enddo     ! loop over fields

     end subroutine compute_corr_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: scale_adj --- Scale adjoint fields so they are physical, 
!                         grid-independent, sensitivities rather than
!                         grid-sensitive sensitivities.
!
!
! !INTERFACE:
!
      subroutine scale_adj_ (imr,jnp,nl,nfields2d,nfields3d,fnames, &
                            jweights,kweights,fields2d,fields3d) 

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)      :: imr,jnp,nl,nfields2d,nfields3d
      character(len=*), intent(in) :: fnames(nfields3d+nfields2d) 
      real(r8), intent(in)     :: jweights(jnp,2), kweights(nl)

! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout)  :: fields2d(imr,jnp,1 ,nfields2d)
      real(r8), intent(inout)  :: fields3d(imr,jnp,nl,nfields3d)

! !DESCRIPTION:
!
!  Scale adjoint fields so they are physical, grid-independent,
!  sensitivities rather than grid-sensitive sensitivities.
!
!  The scaling is based on arguments in the Appendix of the paper
!  Lewis, J.M., K. D. Raeder, and R. M. Errico, Tellus 2001, 53A,
!  pages 74-93.
!
!  It attempts to remove most of the potential mis-representation of
!  the adjoint (sensitivity) fields due to grid irregularity. This is
!  important if general physically-based interpretations of the adjoint
!  fields are to be made, as opposed to specific model-grid-based
!  interpretations.  The latter requires only a Euclidean dot product
!  of perturbation and corresponding sensitivity gradient vectors
!  defined on the model grid with no explicit consideration of varying
!  areas and depths of grid boxes. The former instead considers an L2-norm
!  computed as an integral of the product of perturbation and
!  corresponding physically-based sensitivity vectors integrated over
!  "volume" defined by the fractional area of the earth and vertical sigma
!  coordinate.  The latter is defined using the standard vertical grid
!  for the eta-coordinate model assuming a surface pressure of 10**5
!  pascals.
!
!  To transform the sensitivities to this physical representation, it
!  is necessary to divide each of the model-derived sensitivities by
!  the fractional volume (Delta area *delta sigma/Total area) defined
!  for each grid box.  The units and magnitudes of the new scaled fields
!  will be markedly different than the original ones, and the inner
!  products used to compute first--order variations must include a
!  fractional grid-volume factor at each grid point. The tendency of
!  model-computed adjoint fields to be inversely proportional to
!  fractional volumes will be greatly mitigated by this rescaling.  Also,
!  the newly scaled sensitivities will be relatively grid independent,
!  reflecting only how the resolution impacts the dynamics, not the
!  definition of the fields themselves. For further discussion see the
!  referenced paper and other citations within it.
!
!  It is assumed that values of scalar fields at the poles output by the
!  adjoint model have already been adjusted by averaging them over the
!  index i and re-setting them all to that value. Pole values for all i are
!  therefore identical and each value represents the sensitivity by a single
!  polar sector that has 1/imr the area of the polar cap. For this reason, 
!  the rescaling here has the same dependence on imr at polar points and 
!  non-polar points. The true sensitivity to the pole value is actually the 
!  sum of the imr pole values if the pole is instead considered as a single 
!  polar cap. The averaging of pole values that has been assumed to 
!  have already been applied and the optional rescaling here should ensure 
!  that both the unscaled and rescaled scalar adjoint fields should be 
!  continuous at the poles and not reflect the difference between the grid's
!  single polar cap and the longitudinal array of grid boxes at all other 
!  latitudes. 
! 

! !SEE ALSO: m_jacobian for further description of the interpretation of 
! the treatment of redundant values at the poles.
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  Winslow    Modified to work with GEOS4 I/O
!  20Oct2005  R. Errico  Corected typo in hfac calculation ('*jfacsum')
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer  :: i,j,k,n   ! grid point and field indexes
      integer  :: idu,idv   ! indicators for u,v fields
      integer  :: j1,j2     ! starting and ending values of j 
      integer  :: jgrid     ! =1 for u field, =2 all other fields
      real(r8) :: dbleimr, jfacsum
      real(r8) :: hfac(jnp,2)  ! fraction of earth's area in grid box    
      real(r8) :: hkfac        ! fraction of sigma and area in grid volume
!
!  Specify fractional area of each gridpoint
      dbleimr=dble(imr)
      jfacsum=0.d0
      do j=1,jnp
        jfacsum=jfacsum+jweights(j,2)
      enddo
      jfacsum=1.d0/(dbleimr*jfacsum)
      do j=1,jnp
        hfac(j,1:2)=jweights(j,1:2)*jfacsum
      enddo
!
!  Identify u and v fields since they are not on T,delp,q etc. grid
      call find_name (nfields3d,fnames,'u__',idu)      
      call find_name (nfields3d,fnames,'v__',idv)      
!
!  Loop over 3-D fields
!
      do n=1,nfields3d
!
!  Specify start and end values for grid and which weights to use
        j1=1
        j2=jnp       
        jgrid=2
        if (n == idu) then
          j1=1
          jgrid=1
        else if (n == idv) then
          j1=2
          j2=jnp-1
        endif

        do k=1,nl
          do j=j1,j2
            hkfac=1.d0/(hfac(j,jgrid)*kweights(k))
            do i=1,imr
              fields3d(i,j,k,n)=hkfac*fields3d(i,j,k,n)
            enddo   ! loop over i
          enddo     ! loop over j
        enddo       ! loop over k
      enddo         ! loop over n
!
!  Loop over 2-D fields  (assumption is all are on T,delp,q grid)
!
      j1=1
      j2=jnp       
      jgrid=2

      do n=1,nfields2d
        do j=j1,j2
          hkfac=1.d0/hfac(j,jgrid)
          do i=1,imr
            fields2d(i,j,1,n)=hkfac*fields2d(i,j,1,n)
          enddo   ! loop over i
        enddo     ! loop over j
      enddo         ! loop over n

      end subroutine scale_adj_
!
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: StatsSet_ --- General paramater setup for singular vector calc 
!
! !INTERFACE:
!
      subroutine StatsSet_ ( comm, root, &
                             ntimes, nymd, nhms, corr_kind, rescale_adj, &
                             e_calc,                                     &
                             RCfile, stat )
                                                                                                    
! !USES:
                                                                                                    
      Implicit None
                                                                                                    
! !INPUT PARAMETERS:
!
      integer, intent(in) ::  comm      ! MPI communicator
      integer, intent(in) ::  root      ! MPI ROOT PE
      integer, intent(in) ::  ntimes    ! dim(nymd,nhms)
                                                                                                    
      character(len=*), intent(in) :: RCfile  ! resource filename
                                                                                                    
! !OUTPUT PARAMETERS:
                                                                                                    
      integer, intent(out) :: nymd(ntimes)                 ! date (YYYYMMDD) of stats
      integer, intent(out) :: nhms(ntimes)                 ! time   (HHMMSS) of stats
      integer, intent(out) :: corr_kind(2)                 ! slots to be correlated
      logical, intent(out) :: rescale_adj     ! true if adm field comparison
                                              ! should extend through the column
      logical, intent(out) :: e_calc          ! true if energy norm to be calc.
      integer, intent(out) :: stat            ! Return error code:
                                              !  stat=1 cannot determined proc
                                              !  stat=2 RC file not found
                                              !  stat=3 missing RC name indicator
                                              !  stat=4 missing RC entry
             
             
! !FILES USED:  Resource file RCfile
!
! !DESCRIPTION:  Initializes perturbation used for sensitivity runs.
!
! !REVISION HISTORY:
!
!   30Jul2004  Winslow    Initial code, based on similar codes
!   05Nov2007  Todling    Add iopt to remove redundancy in command line
!
!EOP
!-------------------------------------------------------------------------
      character(len=*), parameter :: myname_ = myname//':StatsSet_'
             
      character(len=255) token
      character(len=255) statsrc
             
      integer       i, j, k, iret, ierr
      integer       iopt
      integer       myID
             
      stat = 0
             
      call MP_comm_rank(comm,myID,ierr)
        if (ierr/=0) then
            stat = 1  ! cannot grab proc id
            return
        endif
             
!     Load resources from RC file
!     ---------------------------
      statsrc = ' '
      call getenv('STATS_RC',statsrc)           ! Unix binding
      if(statsrc.eq.' ') statsrc=RCfile           ! default name
      call i90_loadf (trim(statsrc), iret)
      if( iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_loadf error, iret =',iret
          stat = 2 ! missing RC file
          return
      end if
      if(myID==ROOT) then
         write(stdout,'( a  )') '--------------------------------------------'
         write(stdout,'(2a  )') myname_, ': Reading resource file'
         write(stdout,'( a,/)') '--------------------------------------------'
      end if

      call I90_label('iopt_stats:', iret)
      if (iret .eq. 0) then
           iopt = I90_GInt(iret)
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
             stat = 2 ! missing RC file
             return
           endif
      else
            iopt = 2 
            write(stdout,'()')  myname_, 'Using default case, iopt_stat = ', iopt
      end if
      if ( iopt/=2 .and. iopt/=3 .and. iopt/=5 .and. iopt/=6 ) then
          write(stderr,'(2a,i5)') myname_,': invalid option, iopt_stats =',iopt
          stat = 9 ! messed up option
          return
      endif
      nOutfiles = iopt
      allocate( knames(nOutfiles), outNames(nOutfiles) )

!     Define type of File1 and 2
!     --------------------------
      if ( nOutfiles .ge. 2 ) then
         call I90_label('file_type1:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             knames(1) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File type for first output file: ', trim(knames(1))

         call I90_label('file_type2:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             knames(2) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File type for second output file: ', trim(knames(2))
      endif
             
!     Define type of File3
!     --------------------
      if ( nOutfiles .ge. 3 ) then
         call I90_label('file_type3:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             knames(3) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File type for third output file: ', trim(knames(3))
      endif
                                                                                                 
!     Define type of File4
!     --------------------
      if ( nOutfiles .ge. 4 ) then
         call I90_label('file_type4:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             knames(4) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File type for fourth output file: ', trim(knames(4))
      endif
                                                                                                 
!     Define type of File5
!     --------------------
      if ( nOutfiles .ge. 5 ) then
         call I90_label('file_type5:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             knames(5) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File type for fifth output file: ', trim(knames(5))
      endif
                                                                                                 
!     Define type of File6
!     --------------------
      if ( nOutfiles .ge. 6 ) then
         call I90_label('file_type6:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             knames(6) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File type for sixth output file: ', trim(knames(6))
      endif
                                                                                                 
             
!     Define output File1 and 2
!     -------------------------
      if ( nOutfiles .ge. 2 ) then
         call I90_label('file_name1:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             outnames(1) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for first file: ', trim(outnames(1))
                                                                                                 
         call I90_label('file_name2:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             outnames(2) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for second file: ', trim(outnames(2))
      endif

      if ( nOutfiles .ge. 3 ) then
         call I90_label('file_name3:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             outnames(3) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for third file: ', trim(outnames(3))
      endif

      if ( nOutfiles .ge. 4 ) then
         call I90_label('file_name4:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             outnames(4) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for fourth file: ', trim(outnames(4))
      endif

      if ( nOutfiles .ge. 5 ) then
         call I90_label('file_name5:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             outnames(5) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for fifth file: ', trim(outnames(5))
      endif

      if ( nOutfiles .eq. 6 ) then
         call I90_label('file_name6:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             outnames(6) = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'File name for sixth file: ', trim(outnames(6))
      endif

!     Define correlation slots
!     ------------------------
      call I90_label('corr_kind1:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,': I90_label error(corr_kind1), iret =',iret
        stat = 3 ! missing RC name indicator
        return
      end if
      corr_kind(1) = I90_GInt(iret)
      if( iret .ne. 0) then
         write(stderr,'(2a,i5)') myname_,': I90_GInt error(corr_kind(1)), iret =',iret
         stat = 4 ! missing RC entry
         return
      end if

      call I90_label('corr_kind2:', iret)
      if (iret .ne. 0) then
        write(stderr,'(2a,i5)') myname_,': I90_label error(corr_kind2), iret =',iret
        stat = 3 ! missing RC name indicator
        return
      end if
      corr_kind(2) = I90_GInt(iret)
      if( iret .ne. 0) then
         write(stderr,'(2a,i5)') myname_,': I90_GInt error(corr_kind(2)), iret =',iret
         stat = 4 ! missing RC entry
         return
      end if

                                                                                                 
!     Define output dates
!     -------------------
      if ( ntimes .ge. 1 ) then
        call I90_label('date1:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_label error(date1), iret =',iret
          stat = 3 ! missing RC name indicator
          return
        end if
        nymd(1) = I90_GInt(iret)
        if( iret .ne. 0) then
           write(stderr,'(2a,i5)') myname_,': I90_GInt error(nymd(1)), iret =',iret
           stat = 4 ! missing RC entry
           return
        end if
      endif
             
      if ( ntimes .ge. 2 ) then
        call I90_label('date2:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_label error(date2), iret =',iret
          stat = 3 ! missing RC name indicator
          return
        end if
        nymd(2) = I90_GInt(iret)
        if( iret .ne. 0) then
           write(stderr,'(2a,i5)') myname_,': I90_GInt error(nymd(2)), iret =',iret
           stat = 4 ! missing RC entry
           return
        end if
      endif
                                                                                                 
      if ( ntimes .ge. 3 ) then
        call I90_label('date3:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_label error(date3), iret =',iret
          stat = 3 ! missing RC name indicator
          return
        end if
        nymd(3) = I90_GInt(iret)
        if( iret .ne. 0) then
           write(stderr,'(2a,i5)') myname_,': I90_GInt error(nymd(3)), iret =',iret
           stat = 4 ! missing RC entry
           return
        end if
      endif

      if ( ntimes .ge. 4 ) then
        call I90_label('date4:', iret)
        if (iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_label error(date4), iret =',iret
          stat = 3 ! missing RC name indicator
          return
        end if
        nymd(4) = I90_GInt(iret)
        if( iret .ne. 0) then
           write(stderr,'(2a,i5)') myname_,': I90_GInt error(date4), iret =',iret
           stat = 4 ! missing RC entry
           return
        end if
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

!     Decide if you want to do the adjoint rescaling
!     ----------------------------------------------
        call I90_label('rescale_adj:', iret)
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
        rescale_adj = .false.
        if( trim(token) .eq. 'yes' ) rescale_adj = .true.
        if(myID==ROOT) write(stdout,'(2a)')   &
                     'Rescale the adjoint: ', trim(token)

!     Decide if you want to calculatethe e-norm
!     -----------------------------------------
        call I90_label('e_calc:', iret)
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
        e_calc = .false.
        if( trim(token) .eq. 'yes' ) e_calc = .true.
        if(myID==ROOT) write(stdout,'(2a)')   &
                     'Calculate the energy norm: ', trim(token)


!     release resource file:
!     ---------------------
      call I90_release()
             
      return
      end subroutine StatsSet_
             
      subroutine StatsClean_ ()
      deallocate (knames,outNames)
      end subroutine StatsClean_
             
      end module m_postStats

!
