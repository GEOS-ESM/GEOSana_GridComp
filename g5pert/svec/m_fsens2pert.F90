!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_fsens2pert --- unscales output from ADM runs
!
! !INTERFACE:

      module m_fsens2pert

      use precision
      
      use m_const, only : Cp     => cpm
      use m_const, only : R      => rgas
      use m_const, only : kappa
      use m_const, only : zvir
      use m_const, only : pstd
      use m_const, only : tstd

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type

      use m_mapz, only : set_eta

      use m_dyn

      use m_pertutil, only: find_name          => pertutil_find_name
      use m_pertutil, only: calculate_ps       => pertutil_calculate_ps
      use m_pertutil, only: calculate_ps_ad    => pertutil_calculate_ps_ad
      use m_pertutil, only: transform_T        => pertutil_transform_T
      use m_pertutil, only: transform_T_ad     => pertutil_transform_T_ad
      use m_pertutil, only: enorm              => pertutil_enorm
      use m_pertutil, only: inverse_grad_enorm => pertutil_inverse_grad_enorm
      use m_pertutil, only: horiz_grid         => pertutil_horiz_grid
      use m_pertutil, only: vert_grid          => pertutil_vert_grid
      use m_pertutil, only: read_data          => pertutil_read_data
      use m_pertutil, only: write_data         => pertutil_write_data
      use m_pertutil, only: getparam           => pertutil_getparam
      use m_pertutil, only: setparam           => pertutil_setparam

      use m_stdio
      use m_inpak90
      use m_die, only: MP_die, die

      use m_ioutil, only : luavail

      implicit NONE
 
! !PUBLIC MEMBER FUNCTIONS:
 
      PRIVATE

      PUBLIC fsens2pert_set
      PUBLIC fsens2pert_run

      interface fsens2pert_run; module procedure fsens2pert_; end interface
      interface fsens2pert_set; module procedure pertset_; end interface

!
! !DESCRIPTION: Energy unscaling of sensitivity fields 
!
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm for m_initadj.F90
!  08Aug2005  RG/EN      Adopted m_initadj.F90 to create this stand-alone application
!                        Purpose: transform ADM/SENS outputs into a perturbation dynamic vector
!  08Jan2008  Todling    Major unification of post-proc codes
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'm_fsens2pert'
 
      CONTAINS

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fsens2pert_ --- Convert adjoint fields into a perturbation fileds
!
!
! !INTERFACE:
!
      subroutine fsens2pert_(nfiles,dynfiles,pertfile,knames,nymd,nhms)
                                                                                                                           
!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: nfiles
      character(len=*), intent(in) :: dynfiles(nfiles)
      character(len=*), intent(in) :: pertfile
      character(len=*), intent(in) :: knames  (nfiles+1)
      integer, intent(inout) :: nymd 
      integer, intent(inout) :: nhms 

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION:
 
! !SEE ALSO:
!
 
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  08Aug2005  RG/EN      Adopted m_initadj.F90 to create this stand-alone application
!                        Purpose: transform ADM/SENS outputs into a perturbation dynamic vector 
!  01Nov2005  RG/EN      Changed #of tracers in output perturbation vector to match #of tracers
!                                           in the reference state vector
!                        Writing out the result on the grid 
!                        with the dimensions of a reference state vector
!                        Note: Match of FSENS and REF vector dimensions
!                        have been checked, when files are read
!  08Jan2008 Todling     Revisited most of the interfaces
!  18Dec2008 Errico      Bug fix call arg of read_data (adm instead of tlm) 
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer :: im2,im3   ! grid points in zonal direction from reference file(s)
      integer :: jn2,jn3   ! grid points in lat. direction from reference file(s)
      integer :: nl2,nl3   ! grid points in vert. direction from reference file(s)
      integer :: lm2,lm3   ! number of tracers from reference file(s)

!  Any of the following options can be used for slots 2 and 3: RTL1, RTL2, RAD1, RAD2, NLM1, NLM2 

      integer, parameter :: nfields3d=5   ! number of 3D fields 
      integer, parameter :: nfields2d=1   ! number of 2D fields 
      integer, parameter :: nfields=nfields3d+nfields2d
      character(len=3), dimension(nfields) :: fnames 
      character(len=3), dimension(nfields) :: fnames_def = &
                  (/ 'u__','v__','vpT','q__','dp_','ps_' /)
       
      integer  :: k_ref ! index for reference field used to define J
                        ! Generally, this is either reference field or diff 
      integer  :: ndT, ndp, nps, kdp  ! indicates particular fields
      integer  :: id_field        ! index for single field used to define J

      real(r8), allocatable :: fields2d1(:,:,:,:,:)
      real(r8), allocatable :: fields3d1(:,:,:,:,:)

      real(r8), allocatable :: aux4d (:,:,:,:)
      real(r8), allocatable :: aux2d (:,:)
 
      real(r8), allocatable :: jweights(:,:) ! area weights (1) u, (2) T
      real(r8), allocatable :: kweights(:)   ! vertical layer weights
      real(r8), allocatable :: glats(:,:)    ! degrees lats (1) u, (2) T
      real(r8), allocatable :: pmid(:)       ! p at data levels for standard grid
      real(r8), allocatable :: delbk(:)      ! "thickness" of ps dependent portion of
                                             ! standard (CCM) hybrid-coord. vert. grid
      real(r8), allocatable :: ak(:)         ! hybrid coordinates
      real(r8), allocatable :: bk(:)         ! hybrid coordinates
      integer  :: m                          ! size of long vector
      real(r8) :: ptop                       ! pressure at top of grid
      integer  :: ks                         ! level of stratosphere

      character(len=3) tname
      logical  :: trans_T         ! flag that transforms between T and theta_v required
      integer :: nu, nv           ! pointers for field data
      integer :: nq
      integer :: i,j,k,nf         ! looping indices

      integer, parameter :: comm = 0 ! MPI communicator
      integer, parameter :: ROOT = 0 ! MPI ROOT PE
      integer :: ier              ! I/O error code
      integer :: kinds            ! overall number of state-vect-like arrays
      integer :: lunit1           ! logical unit number for ASCII output
      real(r8) :: lambda
      real(r8) ::  e(2)
      logical  :: exists

      fnames = fnames_def
      call find_name (nfields3d,fnames,'vpT',ndt)
      call getparam ( 'tname', tname )
      fnames(ndt) = tname
      print *, 'Temperature field on input set to: ', tname

      call find_name  (nfields3d,fnames,'dp_',ndp)
      call find_name  (nfields3d,fnames,'u__',nu )
      call find_name  (nfields3d,fnames,'v__',nv )
      call find_name  (nfields,  fnames,'ps_',nps)

      kinds = nfiles + 1

!  Read in reference field(s)

      call dyn_getdim  (trim(dynfiles(1)),im2,jn2,nl2,lm2,ier)  ! DIM of ADM/SENS vector
      call dyn_getdim  (trim(dynfiles(2)),im3,jn3,nl3,lm3,ier)  ! DIM of ref state vector
       
        if (nl3.ne.nl2.or.im3.ne.im2.or.jn3.ne.jn2) then
          print*,'nl2,im2,jn2', nl2,im2,jn2
          print*,'nl3,im3,jn3', nl3,im3,jn3
          print*,'interpolator runs externally'
          call die ('FSENS and REF vector dimensions do not match')
        endif
        print*, 'Will force to usual number of 4'
        print*
        lm2 = 4
        lm3 = lm2

!   Having grid info available, allocate weighting terms

      allocate (jweights(jn2,2))
      allocate (kweights(nl2))
      allocate (glats(jn2,2))
      allocate (pmid(nl2))
      allocate (delbk(nl2))
      allocate (ak(nl2+1))
      allocate (bk(nl2+1))

!   Compute weighting terms

      call horiz_grid (im2, jn2, jweights, glats)
      call vert_grid  (nl2, ks, ak, bk, kweights, ptop, pmid, delbk)

!   Allocate arrays for fields

      allocate(fields2d1(im2,jn2,  1,nfields,kinds))
      allocate(fields3d1(im2,jn2,nl2,nfields,kinds))
      fields3d1(:,:,:,:,:)=0.d0
      fields2d1(:,:,:,:,:)=0.d0

! Read in ADM/SENS run output from fsens.nc4 file

       call read_data (im2,jn2,nl2,nfields,fnames,nymd,nhms,trim(dynfiles(1)), &
                       fields3d1(:,:,:,:,1:1),'adm',ier)
            if(ier/=0) call die(myname,'Error reading SENS/ADM file')

!  Read in reference state field

       call read_data (im2,jn2,nl2,nfields,fnames,nymd,nhms,trim(dynfiles(2)), &
                       fields3d1(:,:,:,:,3:3),'nlm',ier)
            if(ier/=0) call die(myname,'Error reading REF. STATE file')

! Convert VPT to T in sens run result:

! For adjoint fnames(ndT) should be set to 'T__' on output

       call transform_T_ad (im2,jn2,nl2,ptop,nfields,fnames,   &
                            fields3d1(:,:,:,:,1:1), fields3d1(:,:,:,:,3:3))

! Calculate REF STATE SURFACE PRESSURE: call to calculate_ps is replaced by sum

      call calculate_ps (im2, jn2, nl2, ptop, fields3d1(:,:,:,ndp,3), fields2d1(:,:,1,nps,3))

! Calculate adjoint fields PERTURBED SURFACE PRESSURE

      call calculate_ps_ad (im2, jn2, nl2,fields3d1(:,:,:,ndp,1),delbk,fields2d1(:,:,1,nps,1))
 
! Convert a vector of adjoint fields to perturbation fields by unscaling with E^(-1) matrix
 
      call inverse_grad_enorm (im2,jn2,nl2,nfields,bk,                &
                               jweights,ptop,                         &
                               fields3d1(:,:,:,:  ,1),                &
                               fields3d1(:,:,:,ndp,3),                &
                               fields2d1(:,:,1,nps,1),fnames,         &
                               fields3d1(:,:,:,:  ,2),                &
                               fields2d1(:,:,1,nps,2))
 
! Calculate total energy
 
      call enorm (im2,jn2,nl2,nfields,                        &
                  jweights,ptop,                              &
                  fields3d1(:,:,:,:  ,2),                     &
                  fields3d1(:,:,:,ndp,3),                     &
                  fields2d1(:,:,1,nps,2),fnames,              &
                  e(2))

!
      inquire( file=trim(dynfiles(4)), exist=exists )
      if (exists) then
     	 lunit1= luavail()
         open(lunit1, file=trim(dynfiles(4)))
         read(lunit1, *) 
         read(lunit1, *) 
         read(lunit1, '(f20.12)') e(1)
         close(lunit1)      
!
         print *, ' Energy in final   perturbation: e(1) = ', e(1)
         print *, ' Energy in initial perturbation: e(2) = ', e(2)
         if ( e(2) > 0.0 ) then
             lambda = 0.5*e(1)/e(2)     
         else
	     call MP_die (myname,'unacceptable energy',99)
         endif
!
         lunit1= luavail()
         open(lunit1, file='lambda.txt', form='formatted')
         write (lunit1, '(e16.4,6x,A)') lambda, trim(dynfiles(1)) 
         close(lunit1)
!
         print *
         print *,' LAMBDA has been calculated:  ', lambda
         print *,' Output fields have been scaled by LAMBDA' 
         print *

     else

         lambda = 1.d0
         print *
         print *, ' File with final Jnorm value is not found'
         print *, '   Value of lambda is set to 1.d0'
         print *, '    Program will process further '
         print *

      endif

!
      fields3d1(:,:,:,:,2) = lambda*fields3d1(:,:,:,:,2)
      fields2d1(:,:,:,:,2) = lambda*fields2d1(:,:,:,:,2)

!  Transform adjoint T field to adjoint virtual (potential) temperature field

      call transform_T (im2,jn2,nl2,ptop,nfields, fnames,  &
                        fields3d1(:,:,:,:,2),fields3d1(:,:,:,:,3))

! Write output file
      allocate(aux4d(im2,jn2,nl2,nfields3d))
      allocate(aux2d(im2,jn2))
      aux4d(:,:,:,:) = fields3d1(:,:,:,:,2)
      do j = 1,jn2
      do i = 1,im2
      aux2d(i,j) = fields2d1(i,j,1,nps,2)
      enddo
      enddo
      call write_data (im2,jn2,nl2,lm2,nfields2d,nfields3d,ptop,ks,ak,bk,    &
                       aux2d,aux4d,fnames,trim(dynfiles(3)),nymd,nhms,'adm')
      deallocate(aux2d)
      deallocate(aux4d )

      deallocate(jweights,kweights,glats,pmid,delbk,ak,bk)
      deallocate(fields2d1,fields3d1)

      end subroutine fsens2pert_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PertSet_ --- General paramater setup for singular vector calc
! 
! !INTERFACE:
!
      subroutine PertSet_ ( comm, root, RCfile, vectype, stat )
 
! !USES: 
    
      Implicit None

! !INPUT PARAMETERS: 
!
      integer, intent(in) ::  comm	! MPI communicator
      integer, intent(in) ::  root	! MPI ROOT PE
      integer, intent(in) ::  vectype	! Determine dyn-vect type

      character(len=*), intent(in) :: RCfile  ! resource filename

! !OUTPUT PARAMETERS:
 
      integer, intent(out) :: stat            ! Return error code:
                                              !  stat=1 cannot determined proc
                                              !  stat=2 RC file not found
                                              !  stat=3 missing RC name indicator
                                              !  stat=4 missing RC entry


! !FILES USED:  fsens2pert.rc
!
! !DESCRIPTION:  Initializes parameters for fsens2pert program.
!
! !REVISION HISTORY: 
!
!   12Sep2013  Todling   Initial code.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'PertSet_'

      character(len=255) token
      character(len=255) pertrc
      character(len=3  ) evolve_svec

      real(r8)      eps_eer
      integer       i, j, k, jcnt, maxjtype, iret, ierr
      integer       myID
      integer       vnorm

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

      if ( trim(RCfile) == "NONE" ) return

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

!     release resource file:
!     ---------------------
      call I90_release()

      return
      end subroutine PertSet_

 
      end module m_fsens2pert
