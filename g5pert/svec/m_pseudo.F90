!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_pseudo --- Pseudo-inverse fields based on forecast error projection
!                        into a set of final singular vectors
!
! !INTERFACE:
                                                                                                       
      module m_pseudo

      use precision

      use prognostics, only : prognostics_initial
      use prognostics, only : prognostics_final
      use prognostics, only : dyn_prog
                                                                                                       
      use m_const, only : Cp     => cpm   
      use m_const, only : R      => rgas  
      use m_const, only : kappa
      use m_const, only : pstd
      use m_const, only : tstd
                                                                                                       
      use m_pertutil, only : pt2t => pertutil_pt2t
      use m_pertutil, only : ptv2pt => pertutil_ptv2pt
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

      use m_ioutil, only : luavail
                                                                                                       
      implicit NONE
                                                                                                       
! !PUBLIC MEMBER FUNCTIONS:
                                                                                                       
      PRIVATE
                                                                                                       
      PUBLIC pseudo_Set
      PUBLIC pseudo_Process
      PUBLIC pseudo_Clean
                                                                                                       
      interface pseudo_Set ; module procedure PseudoSet_   ; end interface
      interface pseudo_Process; module procedure PseudoInverse_; end interface
      interface pseudo_Clean; module procedure PseudoClean_; end interface
 
!
! !DESCRIPTION: calculation of pseudo-inverse fields 
!
! !REVISION HISTORY:
!
!  07Dec2005  EN/RG     Adopted m_simsv.F90 to create this application, 
!                       which computes forecast error, pseudo-inverse,
!                       projected forecast error, and projected sensitivity.   
!  01May2006  Elena N.  Fixed output format for PALM/altix
!                       Fixed typo bug in reading pseudo.rc file 
!  18May2007  Todling   Fake mpi for now
!  05Jan2008  Todling   Clean up:
!                         - find_name doing error check
!
!EOP
!-------------------------------------------------------------------------
 
      character(len=*), parameter :: myname = 'm_PseudoInverse'
      real(r8), parameter :: pref = 100.d0*pstd
      real(r8), parameter :: tref = tstd

      integer, parameter :: COMM = 0  ! fake MPI
      integer, parameter :: ROOT = 0  ! fake MPI
 
      CONTAINS
 
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PseudoInverse - compute similarity of two sets of singular vectors 
!
!
! !INTERFACE:
!
      subroutine PseudoInverse_ ( nvecsA, nvecsB,     &
                         jobR, FSVEC_name, ISVEC_name, svalu_file, jobB, dynfiles, RCfile)

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: nvecsA         ! number of SVECs to be used for projection
      integer, intent(in) :: nvecsB      
      character(len=*), intent(in) :: jobR
      character(len=*), intent(in) :: FSVEC_name    ! name of the first FSVEC is used as a template 
      character(len=*), intent(in) :: ISVEC_name    ! name of the first ISVEC is used as a template
      character(len=*), intent(in) :: jobB(nvecsB)
      character(len=*), intent(in) :: svalu_file 
      character(len=*), intent(in) :: dynfiles(4)   ! intent to output 4 dynamic vectors:
                                                    ! 1. forecast error
                                                    ! 2. forecast error projected on nvecsA EVOLVED SVECs
                                                    ! 3. pseudo-inverse fields
      character(len=*), intent(in) :: RCfile

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:
       
! !DESCRIPTION:
! compute similarity of two sets of singular vectors

! !SEE ALSO:
!
       
! !REVISION HISTORY:
!
!  05Dec2005  EN/RG  Adopted similarity calculation code to create this application
!  01May2006  Elena N. Fixed formatted output for PALM/altix
!  05Jan2008  Todling  Updated interface to e_cross_prod to simplify programs
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
! --- Dimensions of the vectors on set B
!
      integer :: im2        ! grid points in zonal direction
      integer :: jn2        ! grid points in lat. direction
      integer :: nl2        ! grid points in vert. direction
      integer :: lm2        ! number of tracers from reference file(s)

      integer :: nhms, nymd
      integer :: nhms_ISVEC, nymd_ISVEC
      integer :: nhms_FSVEC, nymd_FSVEC
!
      integer,  dimension(2), parameter :: kdomainps=(/1,1/) ! Do not change
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
      real(r8), allocatable :: GL_domain(:,:,:)  ! mask for global norm calculation
      
      integer  :: m                          ! size of long vector
      integer ::  ndim            ! size of state vector
      integer ::  dims(5)         ! array of grid parameters
      integer :: nev              ! dummy argument (number of eigenvectors)
      integer :: ncv              ! dummy argument (number of Lanczos vectors)
      integer :: projrc           ! return flag for projection operator
      integer :: stat 
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

      real(r8), allocatable  :: FSVEC(:,:,:,:)
      real(r8), allocatable  :: ISVEC(:,:,:,:)

      real(r8), allocatable  :: forecast(:,:,:,:)
      real(r8), allocatable  :: analysis(:,:,:,:)
      real(r8), allocatable  :: reference(:,:,:,:)

      real(r8), allocatable  :: aux1(:,:,:,:,:)
      real(r8), allocatable  :: aux2(:,:,:,:,:)
      real(r8), allocatable  :: psfieldA(:,:)
      real(r8), allocatable  :: psfieldB(:,:)
      real(r8), allocatable  :: V_norm(:)
      real(r8), allocatable  :: C_coefs(:)
      real(r8), allocatable  :: V_eigen(:)
      real(r8), allocatable  :: projferr(:,:,:,:)
      real(r8), allocatable  :: pseudo(:,:,:,:)
      real(r8), allocatable  :: sens_trun(:,:,:,:)
      real(r8) :: fcst_norm

      integer  :: kdomain(2)                   ! indexes corresp to pdomain
      integer  :: VGL_domain(2)                ! vertical domain for global norm calculation
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
      real(r8), allocatable  :: energy(:,:)    !

      character(len=256) :: header
      logical :: exists
      integer :: ier

! 
!                      BEGIN EXECUTION
!

      fnames = fnames_def
      call find_name (nfields,fnames,'dp_',ndp)


! Set some grid information
      call dyn_getdim   (jobR,im1,jn1,nl1,lm1, ier)
         print *
         print *, ' DIMENSIONS in reference state vector: ', im1,jn1,nl1,lm1
         print *
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
      call prognostics_initial ( ypert )
      ypert%u  = 1.0
      ypert%v  = 1.0
      ypert%pt = 1.0

      call proj_init (  comm, root, stat=projrc, dims=dims, RCfile=RCfile )
      if(projrc/=0) print*,'projrc = ',projrc,' proj_init'
                                                                                      
      call proj_svec ( COMM, ROOT, ypert, projrc )   ! Apply local projection operator
      if(projrc/=0) print*,'projrc = ',projrc,' proj_svec'
                 
      call proj_clean( projrc )
      if(projrc/=0) print*,'projrc = ',projrc,' proj_clean'
                 
      allocate (udomain   (im1,jn1,nl1))
      udomain(:,:,:)    =  ypert%u          ! ---> fishy
      allocate (vdomain   (im1,jn1,nl1))
      vdomain(:,:,:)    =  ypert%v          ! ---> fishy
      allocate (domain    (im1,jn1,nl1))
      domain(:,:,:)   =  ypert%pt           ! ---> fishy
      allocate (psdomain    (im1,jn1,nl1))
      psdomain(:,:,:)   =  ypert%pt         ! ---> fishy

      call prognostics_final ( ypert )

      allocate (GL_domain (im1,jn1,nl1))
      GL_domain(:,:,:) = 1.d0
      VGL_domain(1) = 1
      VGL_domain(2) = nl1

      allocate (FSVEC(im1,jn1,nl1,nfields))
      allocate (ISVEC(im1,jn1,nl1,nfields))

      allocate (forecast(im1,jn1,nl1,nfields))
      allocate (analysis(im1,jn1,nl1,nfields))
      allocate (reference(im1,jn1,nl1,nfields))

      allocate (projferr(im1,jn1,nl1,nfields))
      allocate (pseudo(im1,jn1,nl1,nfields))
      allocate (sens_trun(im1,jn1,nl1,nfields))

      allocate (aux1(im1,jn1,nl1,nfields,2))
      allocate (aux2(im1,jn1,nl1,nfields,2))
      allocate (psfieldA(im1,jn1))
      allocate (psfieldB(im1,jn1))
      allocate (V_norm(nvecsA))
      allocate (C_coefs(nvecsA))
      allocate (V_eigen(nvecsA)) 
      allocate (jweights(jn1,2)) ! area weights (1) u, (2) T
      allocate (kweights(nl1))   ! vertical layer weights
      allocate (glats(jn1,2))    ! degrees lats (1) u, (2) T
      allocate (pmid(nl1))       ! p at data levels for standard grid
      allocate (delbk(nl1))      ! "thickness" of ps dependent portion of
      allocate (ak(nl1+1))    
      allocate (bk(nl1+1))    
      allocate (energy(nl1+2,3))

!   Compute kdomain to maintain as much code history as possible
      call comp_kdomain (im1, jn1, nl1, domain, kdomain)
      k = kdomain(1)
      psdomain(:,:,1)   =  domain(:,:,k)

!    Snag some grid information
      call horiz_grid (im1, jn1, jweights, glats)
      call vert_grid  (nl1, ks, ak, bk, kweights, ptop, pmid, delbk)
!
! Read reference state into array fields3dA(:,:,:,:,2) and copy into B
!
      call dyn_getdim  (FSVEC_name,im2,jn2,nl2,lm2, ier)
      call compare_dim_(im1,jn1,nl1,lm1,im2,jn2,nl2,lm2, ' Reference state and SVECs ')
      call dyn_getdim  (ISVEC_name,im2,jn2,nl2,lm2, ier)

      call read_data2_ (im1,jn1,nl1,lm1,nfields,fnames, &
                        jobR, .true., nhms_FSVEC, nymd_FSVEC, &
                        reference(:,:,:,:), stat )

      if (stat .ne. 0) then
         print *, ' Reference state vector does not have a time stamp'
         print *, '      matching the time stamp of the FSVEC ' 
         print *, ' Program will read the last time stamp of the reference vector'
         print *, '             and proceed further '
         call read_data2_ (im1,jn1,nl1,lm1,nfields,fnames, &
                        jobR, .false., nhms_FSVEC, nymd_FSVEC, &
                        reference(:,:,:,:), stat )
      endif
!
! Vectors B is a FORECAST ERROR VECTOR
!
        nvecB=1

!
! Read vector number nvecB from set B into fields3dB(:,:,:,:,1)  
!
         call dyn_getdim (jobB(1),im2,jn2,nl2,lm2, ier)
         call compare_dim_(im1,jn1,nl1,lm1,im2,jn2,nl2,lm2, ' Reference state and forecast ')
!
! --- (EN) This vector is a FORECAST vector
!
        call read_data2_ (im2,jn2,nl2,lm2,nfields,fnames, &
                          jobB(nvecB), .false., nhms_FSVEC, nymd_FSVEC,      &
                          forecast(:,:,:,:), stat )

        aux2(:,:,:,:,1) = forecast(:,:,:,:)

        call transform_T2_ (im2, jn2, nl2, ptop, nfields, fnames,   &
                          aux2(:,:,:,:,1))  

!
! --- (EN) This vector is ANALYSIS vector
!
        if (nvecsB > 1) then
!
       		call read_data2_ (im2,jn2,nl2,lm2,nfields,fnames, &
       	                  jobB(nvecB+1), .false., nhms_FSVEC, nymd_FSVEC,      &
               	          analysis(:,:,:,:), stat )

        aux1(:,:,:,:,1) = analysis(:,:,:,:)
!
        call transform_T2_ (im2, jn2, nl2, ptop, nfields, fnames,   &
                          aux1(:,:,:,:,1))
!
        endif
!
        aux2(:,:,:,:,1) = aux2(:,:,:,:,1) - aux1(:,:,:,:,1) 
        aux1(:,:,:,:,2) = reference(:,:,:,:)
        aux2(:,:,:,:,2) = reference(:,:,:,:)

!
! Compute ps fields
        call calculate_ps (im1,jn1,nl1,ptop,aux2(:,:,:,ndp,1),  &
                             psfieldB(:,:))
!
! Compute E-norm for (FORECAST - ANALYSIS) as normalization factor
! It is an energy weighted mean square forecast error 
!
        call e_cross_product (im1,jn1,nl1,nfields,ptop,        &
             udomain, vdomain, domain,  kdomain, reshape(jweights(:,2),(/jn1/)),&
             reshape(aux2(:,:,:,:,1),(/im1,jn1,nl1,nfields/)), &
             reshape(aux2(:,:,:,:,1),(/im1,jn1,nl1,nfields/)), &
             psfieldB(:,:),psfieldB(:,:), &
             reshape(aux2(:,:,:,ndp,2),(/im1,jn1,nl1/)),        &
             fnames(1:nfields), energy )
             fcst_norm = energy(nl1+2,3)
!
        write (*,'(1a,f18.6)') ' ENERGY WEIGHTED MEAN SQUARED FORECAST ERROR: ', fcst_norm
!
! Find dimensions of SVECS using the first one in the set 
! 
        exists = .true.
        call get_svalues_(nvecsA, V_eigen, svalu_file, exists)
        if (exists) V_eigen(:) = sqrt(V_eigen(:))

! Loop over set of vectors A

     do i=1,nvecsA

! Read initial SVEC #i and just keep it to calculate PROJSENS later 
        call make_svec_name(i, ISVEC_name, header)
        call read_data2_ (im1,jn1,nl1,lm2,nfields,fnames, &
                          trim(header), .false., nhms, nymd,      &
                          ISVEC(:,:,:,:), stat )

!
! Read final SVEC #i, keep it and make a working copy to calculate coefs  
        call make_svec_name(i, FSVEC_name, header)
        call read_data2_ (im1,jn1,nl1,lm2,nfields,fnames, &
                          trim(header), .false., nhms, nymd,      &
                          FSVEC(:,:,:,:), stat )

        aux1(:,:,:,:,1) = FSVEC(:,:,:,:) 
!
! Transform T fields 
        call transform_T_tl (im1,jn1,nl1,ptop,nfields,fnames, &
                             aux1(:,:,:,:,1:2) )
!
! Compute ps fields for FSVEC #i
        call calculate_ps (im1,jn1,nl1,ptop,aux1(:,:,:,ndp,1),psfieldA(:,:))
!
! Compute E-norm for A as normalization factor 
        if (.not.exists) then
	        call e_cross_product (im1,jn1,nl1,nfields,ptop,       &
        	   udomain, vdomain, domain, kdomain, reshape(jweights(:,2),(/jn1/)),   &
           	   reshape(aux1(:,:,:,:,1),(/im1,jn1,nl1,nfields/)), & 
                   reshape(aux1(:,:,:,:,1),(/im1,jn1,nl1,nfields/)), &
                   psfieldA(:,:),psfieldA(:,:),&
                   reshape(aux1(:,:,:,ndp,2),(/im1,jn1,nl1/)),       &
          	   fnames(1:nfields), energy )
                   V_eigen(i) = sqrt(energy(nl1+2,3))
        endif

        call e_cross_product (im1,jn1,nl1,nfields,ptop,       &
           GL_domain, GL_domain, GL_domain, VGL_domain, reshape(jweights(:,2),(/jn1/)),   &
       	   reshape(aux1(:,:,:,:,1),(/im1,jn1,nl1,nfields/)), & 
           reshape(aux1(:,:,:,:,1),(/im1,jn1,nl1,nfields/)), &
           psfieldA(:,:),psfieldA(:,:), &
           reshape(aux1(:,:,:,ndp,2),(/im1,jn1,nl1/)),        &
           fnames(1:nfields), energy )
           V_norm(i) = 1.d0/sqrt(energy(nl1+2,3))

! Normalize FSVECS both in T and in THETA 

        aux1(:,:,:,:,1) = aux1(:,:,:,:,1)*V_norm(i)
        FSVEC(:,:,:,:)  = FSVEC(:,:,:,:)*V_norm(i)
      
        call e_cross_product (im1,jn1,nl1,nfields,ptop,       &
           udomain, vdomain, domain, kdomain, reshape(jweights(:,2),(/jn1/)),   &
       	   reshape(aux1(:,:,:,:,1),(/im1,jn1,nl1,nfields/)), & 
           reshape(aux2(:,:,:,:,1),(/im1,jn1,nl1,nfields/)), &
           psfieldA(:,:), psfieldB(:,:), &
           reshape(aux1(:,:,:,ndp,2),(/im1,jn1,nl1/)),        &
           fnames(1:nfields), energy )
           C_coefs(i) = energy(nl1+2,3)

        projferr(:,:,:,:) = projferr(:,:,:,:) +    &
                     FSVEC(:,:,:,:)*C_coefs(i) 

        pseudo(:,:,:,:) = pseudo(:,:,:,:) +        &
                     ISVEC(:,:,:,:)*C_coefs(i)/V_eigen(i)

        sens_trun(:,:,:,:) = sens_trun(:,:,:,:) +  &
                     ISVEC(:,:,:,:)*C_coefs(i)*V_eigen(i)

      enddo     ! loop over A vectors
!
! Print out coefficients
!
      print *
      print *, '             V_NORMS ' 
      write (*,'(12i8)') (i, i=1, nvecsA)
      write (*,'(12f8.4)') (1.d0/(V_norm(i)*V_norm(i)), i=1, nvecsA)  
      print *
      print *, '             C_COEFS '
      write (*,'(12i8)') (i, i=1, nvecsA)
      write (*,'(12f8.4)') (C_COEFS(i), i=1, nvecsA)
      print *
      print *, '             PSEUDO '
      write (*,'(12i8)') (i, i=1, nvecsA)
      write (*,'(12f8.4)') (C_COEFS(i)/V_eigen(i), i=1, nvecsA)
      print *
      print *, '             SENS_TRUNC '
      write (*,'(12i8)') (i, i=1, nvecsA)
      write (*,'(12f8.4)') (C_COEFS(i)*V_eigen(i), i=1, nvecsA)
      print *

! Write out file

      call find_name (nfields, fnames, 'ps_', nps)
!
      aux1(:,:,:,:,1) = forecast(:,:,:,:) - analysis(:,:,:,:)
      psfieldA(:,:) = sum(aux1(:,:,:,ndp,1), 3) 
      call write_data_ (im1,jn1,nl1,lm1,nfields,ptop,ks,ak,bk,   &
                       psfieldA(:,:),aux1(:,:,:,:,1),           &
                       trim(dynfiles(1)),fnames,nymd_FSVEC,nhms_FSVEC)

      psfieldA(:,:) = sum(projferr(:,:,:,ndp), 3)
      call write_data_ (im1,jn1,nl1,lm1,nfields,ptop,ks,ak,bk,   &
                       psfieldA(:,:),projferr,           &
                       trim(dynfiles(2)),fnames,nymd_FSVEC,nhms_FSVEC)

      psfieldA(:,:) = sum(pseudo(:,:,:,ndp), 3)
      call write_data_ (im1,jn1,nl1,lm1,nfields,ptop,ks,ak,bk,   &
                      psfieldA(:,:),pseudo,           &
                      trim(dynfiles(3)),fnames,nymd_ISVEC,nhms_ISVEC)

      psfieldA(:,:) = sum(sens_trun(:,:,:,ndp), 3)
      call write_data_ (im1,jn1,nl1,lm1,nfields,ptop,ks,ak,bk,   &
                       psfieldA(:,:),sens_trun,           &
                       trim(dynfiles(4)),fnames,nymd_ISVEC,nhms_ISVEC)

!
!
! Replace elements of simAB by their squared values, sum rows and 
! columns, and compute similarity indices for various subsets

      deallocate ( energy )
      deallocate (reference, forecast, analysis)
      deallocate (FSVEC, ISVEC, projferr, pseudo, sens_trun)
      deallocate (psfieldA,psfieldB)
      deallocate (aux1, aux2)
      deallocate ( V_norm, C_coefs, V_eigen)
      deallocate(jweights,kweights,glats,pmid,delbk,ak,bk)
      deallocate(domain,udomain,vdomain,psdomain)
      deallocate(GL_domain)

      print *,' '
      print *,' '
      print *,' '
      print *,' PROGRAM COMPLETED'
!
      end subroutine PseudoInverse_ 

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
! !ROUTINE: read svalues from file
!
! !INTERFACE:
!
      subroutine get_svalues_(nvecsA, V_eigen, filename, exists)
!
!USES:
                                                                                                                             
      implicit none
                                                                                                                             
! !INPUT PARAMETERS:

      integer, intent(in) :: nvecsA
      character(len=*), intent(in) :: filename

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: V_eigen(nvecsA) 
                                                                                                                             
! !INPUT/OUTPUT PARAMETERS:

      logical, intent(inout) :: exists

! !DESCRIPTION:
                                                                                                                             
! !SEE ALSO:
!
! !REVISION HISTORY:
!
! 05Dec2005     EN    Reading SVALUES from file
!EOP
!-------------------------------------------------------------------------
!
      character(len=255) :: header
      integer :: i, isvec, lunit

      inquire( file=trim(filename), exist=exists )
      if (exists) then
         lunit= luavail()
	 open (lunit, file = trim(filename))
         read (lunit,*) header
         do i=1,nvecsA
            read (lunit,*) isvec, V_eigen(i)
         enddo
         close(lunit)
       else
         print *, ' File with eigen values is not found '
         print *, '  Will calculate SVALUES here '
         print *, ' Program will process further       '
       endif

      return
      end subroutine get_svalues_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
! !ROUTINE: make a name for a singular vector
!
! !INTERFACE:
!
   subroutine make_svec_name(svec_index, example, svec_name)
!
! !SEE ALSO:
!
! !REVISION HISTORY:
!
! 05Dec2005     EN   This routine forms names of singular vector files
!                    using template for a first file 
!
!USES:
                                                                                                                             
      implicit none
                                                                                                                             
! !INPUT PARAMETERS:

   integer, intent(in) :: svec_index
   character(len=*), intent(in) :: example     ! Template for a first SVEC 

! !OUTPUT PARAMETERS:

   character(len=*), intent(out) :: svec_name 
!
!EOP
!---------------------------------------------------------------------------

   integer :: i
   character(len=3) :: ch1

   ch1 = '000'

   svec_name = trim(example)

   if (svec_index .lt. 10) write ( ch1, '(a2,i1)') '00', svec_index
   if ((svec_index .ge. 10) .and. (svec_index .lt. 100) ) write ( ch1, '(a1,i2)') '0', svec_index
   if (svec_index .ge. 100) write ( ch1, '(i3)') svec_index

   svec_name(INDEX(svec_name, '.nc4')-3 : INDEX(svec_name, '.nc4')-1) = ch1 

   return

   end subroutine make_svec_name
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
! !ROUTINE: compare dimensions of two dynamic vectors 
!
! !INTERFACE:
!
  subroutine compare_dim_(im1,jn1,nl1,lm1,im2,jn2,nl2,lm2, types)
! !SEE ALSO:
!
! !REVISION HISTORY:
!
! 05Dec2005     EN   
!
!USES:
                                                                                                                                                                                                                                                           
      implicit none
                                                                                                                                                                                                                                                           
! !INPUT PARAMETERS:

  integer, intent(in) :: im1,jn1,nl1, im2,jn2,nl2,lm2
  character(len=*), intent(in) :: types
!
! !INPUT/OUTPUT PARAMETERS:
!
  integer, intent(inout) :: lm1


         if ((im1.ne.im2).or.(jn1.ne.jn2).or.(nl1.ne.nl2)) then
            print *
            print *, types,  ': DIMENSIONS DO NOT AGREE. '
            print *, '      Exit on error.                      '
            stop
         else
            lm1 = max(lm1, lm2)
            print *
            print *, ' Ignoring difference in number of tracers '
            print *, '   in ', types
            print *, ' DIMENSIONS of output files will correspond '
            print *, '       to the maximum number of tracers ', lm1
            print *
         endif

   return
   end subroutine compare_dim_

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
                            job,pick,nhms,nymd,fields3d,ier)
                                                                                                         
!USES:
                                                                                                         
      implicit none
                                                                                                         
                                                                                                         
! !INPUT PARAMETERS:
                                                                                                         
      integer,  intent(in)  :: im4, jn4, nl4, lm4, nfields
      character(len=*), intent(in) :: fnames(nfields)
      character(len=*), intent(in) :: job
      integer , intent(in)  :: nhms, nymd
      logical, intent(in)  :: pick
                                                                                                         
! !OUTPUT PARAMETER
                                                                                                         
      real(r8), intent(out) :: fields3d(im4,jn4,nl4,nfields)
      integer,  intent(out) :: ier 
                                                                                                         
! !INPUT/OUTPUT PARAMETERS:
                                                                                                         
! !DESCRIPTION:
             
! !SEE ALSO:
!
            
! !REVISION HISTORY:
!
!  04Oct2004  Winslow    Original code
!  15Apr2005  Elena N.   Change if/else options for T/PT/VPT reading to ensure
!                        that temperature field is filled in any case 
!  10Dec2005  Elena N.   Introduced "pick" to have an option whether the requested 
!                        or the last position at the file is read
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables
             
      integer :: nhms_in, nymd_in, indexF
      type(dyn_vect) w
             
      nhms_in = nhms
      nymd_in = nymd
      call dyn_null ( w )

      if (pick) then        ! --- (EN) Read the requested position in the file
           call dyn_get ( trim(job), nymd_in, nhms_in, w, ier, timidx=0 ) 
      else                  ! --- (EN) Read the last position in the file
           call dyn_get ( trim(job), nymd_in, nhms_in, w, ier)   
      endif

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
             

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: transform_T2 --- Transform reference field of virtual potential
!                            to temperature
!
!
! !INTERFACE:
!
      subroutine transform_T2_ (im1,jn1,nl1,ptop,nfields,fnames,fields)
                                                                                                                                                                                                         
!USES:
                                                                                                                                                                                                         
      implicit none
                                                                                                                                                                                                         
! !INPUT PARAMETERS:
                                                                                                                                                                                                         
      integer, intent(in) :: im1,jn1,nl1,nfields
      real(r8), intent(in) :: ptop
      character(len=*), intent(in) :: fnames(nfields)
                                                                                                                                                                                                         
! !OUTPUT PARAMETERS:
                                                                                                                                                                                                         
! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: fields(im1,jn1,nl1,nfields)
                                                                                                                                                                                                         
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
                                                                                                                                                                                                         
      real(r8) :: pk(im1,jn1,nl1)
      real(r8) :: dum1(im1,jn1,nl1), dum2(im1,jn1,nl1)
      integer :: indexT, indexQ, indexP
                                                                                                                                                                                                         
!
! compute T from vPT
!
      call find_name (nfields,fnames,'T__',indexT)
      if (indexT.le.0) then
            call find_name (nfields,fnames,'vpT',indexT)
            if (indexT.le.0) then
      		call find_name (nfields,fnames,'PT_',indexT)
            endif
      endif

      if (indexT .gt. 0 ) then
        call find_name (nfields,fnames,'q__',indexQ)
        call ptv2pt (im1,jn1,nl1,                                      &
                     fields(:,:,:,indexT),fields(:,:,:,indexQ) )
        call find_name (nfields,fnames,'dp_',indexP)
        call pkappa (im1,jn1,nl1,ptop,fields(:,:,:,indexP), &
                     dum1,pk,dum2,'NLM')
        call pt2t (im1,jn1,nl1,fields(:,:,:,indexT),pk)
      endif
      end subroutine transform_T2_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PseudoSet_ --- General paramater setup for singular vector calc
!
! !INTERFACE:
!
      subroutine PseudoSet_ ( comm, root,  &
                               nvecsA, nvecsB, jobR, jobA, jobI, svalu_file, jobB,        &
                               RCfile, stat )
                                                                                                       
! !USES:
                                                                                                       
      Implicit None
                                                                                                       
! !INPUT PARAMETERS:
!
      integer, intent(in) ::  comm      ! MPI communicator
      integer, intent(in) ::  root      ! MPI ROOT PE
      integer, intent(in) ::  nvecsA
      integer, intent(inout) ::  nvecsB
      character(len=*), intent(in) :: RCfile  ! resource filename
                                                                                                       
! !OUTPUT PARAMETERS:
                                                                                                       
      character(len=*), intent(out) :: jobR
      character(len=*), intent(out) :: jobA
      character(len=*), intent(out) :: jobI
      character(len=*), intent(out) :: jobB(nvecsB)
      character(len=*), intent(out) :: svalu_file
      integer, intent(out) :: stat            ! Return error code:
                                              !  stat=1 cannot determined proc
                                              !  stat=2 RC file not found
                                              !  stat=3 missing RC name indicator
                                              !  stat=4 missing RC entry
                                                                                                       
                                                                                                       
! !FILES USED:  
!
! !DESCRIPTION:  
!
! !REVISION HISTORY:
!
!  05Dec2005    EN    Modified a similar routine to have any number of singular vectors
!                     as provided in pseudo.x command line argument -nsvecs NSVECS 
!  01May2006    EN    Fixed bug in svalu_file filename
!
!   
!
!EOP
!-------------------------------------------------------------------------
      character(len=*), parameter :: myname_ = myname//'PseudoSet_'
                                                                                                       
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
             
!     Define filename template for FSVECs
!     ------------------------
      if ( nvecsA .ge. 1 ) then
         call I90_label('FSVEC_filename:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobA = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'Filename template for FSVECs: ', trim(jobA)
      endif

             
!     Define filename template for FSVECs
!     ------------------------
      if ( nvecsA .ge. 1 ) then
         call I90_label('ISVEC_filename:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             jobI = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'Filename template for ISVECs: ', trim(jobI)
      endif

!     Define filename for SVALUES files
!     ---------------------
         call I90_label('SVALUES_file:', iret)
         if (iret .eq. 0) then
           call I90_Gtoken ( token, iret )
           if (iret .ne. 0) then
             write(stderr,'(2a,i5)') myname_, &
                          ': I90_Gtoken error, iret =',iret
           else
             svalu_file = trim(token)
           end if
         end if
         if(myID==ROOT) write(stdout,'(2a)') 'SVALUES are in a file: ', trim(svalu_file)

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
           if(myID==ROOT) write(stdout,'(2a)') 'File name for first file B group: ', trim(jobB(2))
         else
           if(myID==ROOT) write(stdout,'(2a)') 'Only one file for the second group '
           nvecsB = 1
         endif
      endif

                                                                                                                     
!     release resource file:
!     ---------------------
      call I90_release()
 
      return
      end subroutine PseudoSet_
 
      subroutine PseudoClean_ ()
      end subroutine PseudoClean_
 
!
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
      subroutine write_data_ (im1,jn1,nl1,lm1,nfields3d,ptop,ks, &
                             ak,bk, fields2d,fields3d,job1,fnames,      &
                             nymd,nhms)  

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: im1, jn1, nl1, lm1, nfields3d
      real(r8), intent(in)  :: ptop
      integer,  intent(in)  :: ks
      real(r8), intent(in)  :: ak(nl1+1)
      real(r8), intent(in)  :: bk(nl1+1)
      real(r8), intent(in) :: fields2d(im1,jn1)
      real(r8), intent(in) :: fields3d(im1,jn1,nl1,nfields3d)
      character(len=*), intent(in) :: job1
      character(len=3), intent(in) :: fnames(nfields3d)
      integer , intent(in)  :: nhms, nymd
                                                                 
! !OUTPUT PARAMETERS:
 
! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
  
! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Oct2004  Winslow    Initial algorithm
!  07Apr2005  Elena N.   Changed output precision to 32-bit
!  21Jul2005  Elena N.   Removed extraneous parameters
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables
 
      character(len=*), parameter :: myname = 'write_data'
      integer :: ier
      integer, parameter :: frq = 3000
      logical, parameter :: fnew = .true.
      integer :: fnhms, fnymd
      integer :: f_id                            ! field id
      integer :: ierr                            ! return code
      integer :: prec_das                        ! precision for DAS output
                                                 ! 0: 32 bits
                                                 ! 1: 64 bits (default)
      type (dyn_vect) :: w_f
                                                                 
      call dyn_null ( w_f )

      call Dyn_Init ( im1, jn1, nl1, lm1, w_f, ierr, ptop, ks, ak, bk )

      w_f%phis      = 0.d0
      w_f%hs_stdv   = 0.d0
      w_f%Ts        = 0.d0
      w_f%lwi       = 0.d0
      w_f%ps        = fields2d 
      call find_name (nfields3d,fnames,'dp_',f_id )
      w_f%delp      = fields3d(:,:,:,f_id)
      call find_name (nfields3d,fnames,'u__',f_id )
      w_f%u         = fields3d(:,:,:,f_id)
      call find_name (nfields3d,fnames,'v__',f_id )
      w_f%v         = fields3d(:,:,:,f_id)
      call find_name (nfields3d,fnames,'vpT',f_id )
      if (f_id .le. 0) then
       call find_name (nfields3d,fnames,'T__',f_id )
        if (f_id .le. 0) print *, ' No field name vpT or T__ '
      endif
      w_f%pt        = fields3d(:,:,:,f_id)
      call find_name (nfields3d,fnames,'q__',f_id )
      w_f%q         = fields3d(:,:,:,f_id:f_id+lm1-1)


      if ( ierr .ne. 0 ) then
         call MP_die (myname,'Error from Dyn_Init',ierr)
      end if
      prec_das = 0
      fnhms = nhms
      fnymd = nymd

      call dyn_put ( trim(job1), fnymd, fnhms, prec_das, w_f, ierr, & 
                     new=fnew, freq=frq )
        if(ierr/=0) then
          call MP_die (myname,'Error from dyn_put',ierr)
        else
          print*,'Wrote out perturbation to ',trim(job1)
        endif

      call dyn_null ( w_f )

      end subroutine write_data_
                                                                  
!
!
!---------------------------------------------------------------------------
     end module m_pseudo
