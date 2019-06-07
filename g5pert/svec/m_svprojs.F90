!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_svprojs --- Implements Singular Vector Local Projection
!
! !INTERFACE:
!
module m_svprojs

! !USES:

  use precision
  use prognostics, only : dyn_prog
  use prognostics, only : def_imr     => imr
  use prognostics, only : def_jnp     => jnp
  use prognostics, only : def_nl      => nl
  use prognostics, only : nc
  use prognostics, only : def_jfirst  => jfirst
  use prognostics, only : def_jlast   => jlast

  use control, only : control_number

  use m_sVecDef, only : lclproj
  use m_sVecDef, only : projlon
  use m_sVecDef, only : projlat
  use m_sVecDef, only : projlev

  use m_mpif90,only : MP_comm_rank

  use m_inpak90
  use m_stdio,   only: stdout
  use m_stdio,   only: stderr
  use m_die,     only: mp_die

  implicit none

! !PUBLIC MEMBER FUNCTIONS:
 
  PRIVATE

  PUBLIC  proj_init
  PUBLIC  proj_clean
  PUBLIC  proj_svec

  interface proj_init ; module procedure init_     ; end interface
  interface proj_clean; module procedure clean_    ; end interface
  interface proj_svec ; module procedure proj_svec_; end interface

!
! !DESCRIPTION: This module implements the application of the
!               local projection operator.
!
! !REVISION HISTORY:
!
!  20Nov2002  Todling   Initial design/interfaces.
!  12Jul2004  Todling   Redefined setting of imr,jnp,nl,etc.
!
!EOP
!-------------------------------------------------------------------------
  character(len=*), parameter :: myname = 'm_SVnorms'

  integer, dimension(5), parameter :: rc = (/99, &  ! Invalid option
                                             98, &  ! Already initialized
                                             97, &  ! Not yet implemented
                                             96, &  ! Undefined normalization
                                             95 /)  ! Packg. not initialized
  character(len=32), dimension(5), parameter :: rcmsg = (/ &
                     'Invalid option                  ',   &
                     'Already initialized             ',   &
                     'Not yet implemented             ',   &
                     'Undefined normalization         ',   &
                     'Package not initializedn        ' /)

  logical,  save ::  doproj
  logical,  save ::  letinit = .true.

  real(r8), allocatable, save :: coefft(:,:)
  real(r8), allocatable, save :: coeffu(:,:)
  real(r8), allocatable, save :: coeffv(:,:)
  real(r8), allocatable, save :: coeffk(:)

  integer, save :: imr
  integer, save :: jnp
  integer, save :: nl
  integer, save :: jfirst
  integer, save :: jlast

CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_ --- initialize norms
!
! !INTERFACE:

subroutine init_ ( comm, root, stat, dims, verbose, RCfile )

! !USES:

  implicit none

! !INPUT PARAMETERS:

   integer, intent(in)           :: comm
   integer, intent(in)           :: root
   logical, intent(in), optional :: verbose
   integer, intent(in), optional :: dims(5)
   character(len=*), intent(in), optional :: RCfile ! specify resource filename

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:

  integer,  intent(out), optional ::  stat

! !DESCRIPTION: Initialized weights and whatever else for norm utilization.
!
! !REVISION HISTORY:
!
!  24Nov2002  Todling    Initial code.
!  20Dec2002  Gelaro     Properly defined projection operator.
!  13Oct2004  Todling    Added opt rcfile.
!  04Jul2007  Todling    MPI-related details.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::init_'

  real(r8), allocatable ::  sine(:)
  real(r8), allocatable ::  cosp(:)
  real(r8), allocatable ::  sinp(:)
  real(r8), allocatable ::  cose(:)

  real(r8) dp, dl
  real(r8) pi, clon, clat
  real(r8) rlon1, rlon2, rlat1, rlat2
  integer  ndim
  integer  i, j, k, myID
  integer  ierr
  logical  verb

  if(present(stat)) stat=0

  if (.not.letinit) then
      if (present(stat)) then
          stat = rc(2)
          return
      else
          call mp_die(myname_,trim(rcmsg(2)),rc(2)) 
      endif
  endif

  verb = .false.
  if (present(verbose) ) then
      verb = verbose
  endif

  call MP_comm_rank(comm,myID,ierr)
    if(ierr/=0) call mp_die(myname_,'MP_comm_rank()',ierr)


! Read in projection parameters
! -----------------------------
  call set_ ( comm, root, RCfile=RCfile )

! Make sure that ibeg's and iend's are indeed defined 
! ----------------------------------------------------
!_RT  call control_number ( ndim )

! Set local dimensions 
! --------------------
  if (present(dims)) then
     imr     = dims(1)
     jnp     = dims(2)
     nl      = dims(3)
     jfirst  = dims(4)
     jlast   = dims(5)
  else
     imr     = def_imr
     jnp     = def_jnp
     nl      = def_nl
     jfirst  = def_jfirst
     jlast   = def_jlast
  endif

! Calculate cosines; this is called here just to get dl/dp
! --------------------------------------------------------
  allocate( sine(jnp), cosp(jnp), sinp(jnp), cose(jnp), stat=ierr )
	if(ierr/=0) call mp_die(myname_,'Alloc(sine,...)',ierr)

  call setrig(imr, jnp, dp, dl, cosp, cose, sinp, sine)

  deallocate( sine, cosp, sinp, cose, stat=ierr )
	if(ierr/=0) call mp_die(myname_,'Dealloc(sine,...)',ierr)

! Set coefficients 0/1's for local projection operator
! ----------------------------------------------------
  allocate( coefft(imr,jnp), coeffu(imr,jnp), coeffv(imr,jnp),  &
            coeffk(nl), stat=ierr )
    if(ierr/=0) call mp_die(myname_,'Alloc(coefft,...)',ierr)

  if(projlon(1) < 0.) projlon(1) = projlon(1) + 360.
  if(projlon(2) < 0.) projlon(2) = projlon(2) + 360.

  if(projlev(2) > nl) projlev(2) = nl

  pi    = 4. * atan(1.0)
  rlon1 = pi * projlon(1)/180.
  rlon2 = pi * projlon(2)/180.
  rlat1 = pi * projlat(1)/180.
  rlat2 = pi * projlat(2)/180.

! Coefficients for vertical levels
! --------------------------------
  coeffk = 0.0

  do k=1,nl
  if(k >= projlev(1) .and. k <= projlev(2) ) coeffk(k) = 1.0
  end do

! D-grid requires different horizontal coefficents for different variables
! ------------------------------------------------------------------------

! Coefficents for temperature, delta_p and moisture
! -------------------------------------------------
  coefft = 0.0

  if (rlon1 < rlon2) then

     do j = 1, jnp
        clat = -pi*0.5 + (j-1.)*dp
        if ( clat >= rlat1 .and. clat <= rlat2 ) then
          do i = 1, imr
             clon = dl*0.5 + (i-1.)*dl
             if ( clon >= rlon1 .and. clon <= rlon2 ) coefft(i,j) = 1.0
          end do
        end if
     end do

  else

     do j = 1, jnp
        clat = -pi*0.5 + (j-1.)*dp
        if ( clat >= rlat1 .and. clat <= rlat2 ) then
          do i = 1, imr
             clon = dl*0.5 + (i-1.)*dl
             if ( clon >= rlon1 .or. clon <= rlon2 ) coefft(i,j) = 1.0
          end do
        end if
     end do

  endif ! (rlon1 < rlon2 for t)

! Coefficents for u-component
! ---------------------------
  coeffu = 0.0

  if (rlon1 < rlon2) then

     do j = 2, jnp
        clat = -pi*0.5 + dp*0.5 + (j-2.)*dp
        if ( clat >= rlat1 .and. clat <= rlat2 ) then
          do i = 1, imr
             clon = dl*0.5 + (i-1.)*dl
             if ( clon >= rlon1 .and. clon <= rlon2 ) coeffu(i,j) = 1.0
          end do
        end if
     end do

  else

     do j = 2, jnp
        clat = -pi*0.5 + dp*0.5 + (j-2.)*dp
        if ( clat >= rlat1 .and. clat <= rlat2 ) then
          do i = 1, imr
             clon = dl*0.5 + (i-1.)*dl
             if ( clon >= rlon1 .or. clon <= rlon2 ) coeffu(i,j) = 1.0
          end do
        end if
     end do

  endif ! (rlon1 < rlon2 for u) 


! Coefficents for v-component
! ---------------------------
  coeffv = 0.0

  if (rlon1 < rlon2) then

     do j = 2, jnp - 1
        clat = -pi*0.5 + (j-1.)*dp
        if ( clat >= rlat1 .and. clat <= rlat2 ) then
          do i = 1, imr
             clon = (i-1.)*dl
             if ( clon >= rlon1 .and. clon <= rlon2 ) coeffv(i,j) = 1.0
          end do
        end if
     end do

  else

    do j = 2, jnp - 1
        clat = -pi*0.5 + (j-1.)*dp
        if ( clat >= rlat1 .and. clat <= rlat2 ) then
          do i = 1, imr
             clon = (i-1.)*dl
             if ( clon >= rlon1 .or. clon <= rlon2 ) coeffv(i,j) = 1.0
          end do
        end if
     end do

  endif ! (rlon1 < rlon2 for v)

  if (myid==root .and. verb) then

!    Print test
!    ----------
     print *
     print *,'LPO Verical coefficients'
     do k=1,nl
     print '(a,i2,a,f3.1)', 'Lev= ',k,'  Proj Coeff= ',coeffk(k)
     end do

     print *
     print *,'LPO coefficients for T'
     print '(20x,4(i2,18x))', (j, j=10,jnp,10)
     do i=1,imr
     print '(i2,1x,46f2.0)', i,(coefft(i,j), j=1,jnp)
     end do
     print '(20x,4(i2,18x))', (j, j=10,jnp,10)
  
     print *
     print *,'LPO coefficients for u'
     print '(20x,4(i2,18x))', (j, j=10,jnp,10)
     do i=1,imr
     print '(i2,1x,46f2.0)', i,(coeffu(i,j), j=1,jnp)
     end do
     print '(20x,4(i2,18x))', (j, j=10,jnp,10)

     print *
     print *,'LPO coefficients for v'
     print '(20x,4(i2,18x))', (j, j=10,jnp,10)
     do i=1,imr
     print '(i2,1x,46f2.0)', i,(coeffv(i,j), j=1,jnp)
     end do
     print '(20x,4(i2,18x))', (j, j=10,jnp,10)

  endif


! Initialization done
! -------------------
  doproj  = lclproj
  letinit = .false.

end subroutine init_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_ --- Read in projection boundary parameters
!
! !INTERFACE:
!
      subroutine set_ ( comm, root, RCfile )

! !USES:
   
      Implicit None

! !INCLUDES:

! !INPUT PARAMETERS:

      integer, intent(in)                    :: comm
      integer, intent(in)                    :: root
      character(len=*), intent(in), optional :: RCfile ! rc filename

! !OUTPUT PARAMETERS:
!

! !FILES USED:  fvsvec.rc
!
! !DESCRIPTION:  Read in projection boundary parameters.
!
! !REVISION HISTORY:
!
!   26Jul2004  Todling    Moved RC-read from m_setsvecs here; modular.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'set_'
      character(len=*), parameter :: def_projrc = 'fvsvec.rc'

      character(len=255) token
      character(len=255) svecrc
      character(len=3  ) evolve_svec

      integer       myID, ierr, iret
      real          pnext

      character(len=255) projrc
      logical            invalid,prt

      call MP_comm_rank(comm,myID,ierr)
        if(ierr/=0) call MP_die(myname_,'MP_comm_rank()',ierr)

      prt = myID==root

!     Load resources from resource file
!     ---------------------------------
      if (present(RCfile)) then
         projrc = trim(RCfile)
      else
         projrc = ' '
         call getenv('PROJ_RC',projrc)             ! Unix binding
         if(projrc.eq.' ') projrc=def_projrc       ! default name
      endif
      call i90_loadf (trim(projrc), iret)
      if(iret .ne. 0) call mp_die(myname_,': I90_loadf error, iret =',iret)
      if (prt) then
          write(stdout,'( a  )') '---------------------------------'
          write(stdout,'(2a  )') myname_, ': Reading resource file'
          write(stdout,'( a,/)') '---------------------------------'
      endif

!     Read in values defining lat/lon box for local projection operator
!     -----------------------------------------------------------------
      lclproj    = .false.
      projlon(1) =    0.
      projlon(2) =  360.
      projlat(1) =  -90.
      projlat(2) =   90.
      projlev(1) =    1
      projlev(2) =   99
      call I90_label('local_svec_projection:', iret)
      if (iret .ne. 0) then
        if(prt) write(stdout,'(2a)') myname_, &
                     ': Using default local projection operator, i.e., none'
      else
          invalid = .false.
          pnext = I90_GFloat(iret)
          if (iret/=0) invalid = .true.
          if (.not.invalid) then
             projlon(1) = pnext
          end if
          pnext = I90_GFloat(iret)
          if (iret/=0) invalid = .true.
          if (.not.invalid) then
             projlon(2) = pnext
          end if
          pnext = I90_GFloat(iret)
          if (iret/=0) invalid = .true.
          if (.not.invalid) then
             projlat(1) = pnext
          end if
          pnext = I90_GFloat(iret)
          if (iret/=0) invalid = .true.
          if (.not.invalid) then
             projlat(2) = pnext
          end if
          pnext = I90_GFloat(iret)
          if (iret/=0) invalid = .true.
          if (.not.invalid) then
             projlev(1) = pnext
          end if
          pnext = I90_GFloat(iret)
          if (iret/=0) invalid = .true.
          if (.not.invalid) then
             projlev(2) = pnext
          end if
!         Quick sanity check
!         ------------------
          if ( projlat(1) > projlat(2) ) invalid = .true.
          if ( projlev(1) > projlev(2) ) invalid = .true.

          if ( invalid ) then
               projlon(1) =    0.
               projlon(2) =  360.
               projlat(1) =  -90.
               projlat(2) =   90.
               projlev(1) =    1
               projlev(2) =   99
               write(stdout,'(2a,/,a)') myname_, &
                  ': Something went wrong while setting projection box.', &
                  '  Taking default local projection (-180,180) (-90,90) (All Levels)'
          else
               lclproj = .true.
               if(prt) write(stdout,'(a,/,2(a,f7.2,a,f7.2,/),2(a,i7),/)')    &
                  'User specified local projection: ',               &
                  '  From Lon ', projlon(1), ' to Lon ', projlon(2), &
                  '  From Lat ', projlat(1), ' to Lat ', projlat(2), &
                  '  From Lev ', projlev(1), ' to Lev ', projlev(2)
          end if

      endif

!     release resource file:
!     ---------------------
      call I90_release()

end subroutine set_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: clean_ --- clean norms
!
! !INTERFACE:

subroutine clean_ ( stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:

  integer,  intent(out), optional ::  stat

! !DESCRIPTION: Clear whatever needs to be cleared from norm usage.
!
! !REVISION HISTORY:
!
!  24Nov2002  Todling    Initial code.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::clean_'

  integer  ierr

  if(present(stat)) stat=0

  deallocate( coefft, coeffu, coeffv, coeffk, stat=ierr )
    if(ierr/=0) call mp_die(myname_,'Alloc(coefft...)',ierr)

  letinit = .true.

end subroutine clean_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Proj_Svec_ --- Apply local projection operator
!
! !INTERFACE:

subroutine proj_svec_ ( comm, ROOT, pert, stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:

  integer,  intent(in) :: comm
  integer,  intent(in) :: ROOT

! !INPUT/OUTPUT PARAMETERS:

  type(dyn_prog), intent(inout) :: pert

! !OUTPUT PARAMETERS:

  integer,  intent(out), optional ::  stat

! !DESCRIPTION: This routine applies a local projection operator (LPO)
!               to perturbation state vector
!
! To do:
!  \begin{itemize}
!      \item MPI part
!      \item double check this is really correct
!  \end{itemize}
!
! !REVISION HISTORY:
!
!  20Nov2002  Todling    Initial code: projection matrix assumed diagonal.
!  20Dec2002  Gelaro     Properly defined projection operator.
!  16May2007  Todling    Introduced dyn_prog
!  15Oct2015  Holdaway   Bug fix on q part, nc not known to this module
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::proj_svec_'

  integer  i, j, k, m
  integer  ierr
  integer  myID

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ierr)
      if(ierr/=0) call mp_die(myname,'MP_comm_rank()',ierr)

  if(letinit) then
      if (present(stat)) then
          stat = rc(5)
          return
      else
          call mp_die(myname_,trim(rcmsg(5)),rc(5))
      endif
  endif

  if(.not.doproj) return

! Apply horizontal and vertical projection coeffs to perturbations
! ----------------------------------------------------------------

  do j = jfirst, jlast
     do k = 1, nl
        do i = 1, imr
           pert%u(i,j,k) = pert%u(i,j,k) * coeffu(i,j) * coeffk(k)
        end do
     end do
  end do

  do j = jfirst, jlast
     do k = 1, nl
        do i = 1, imr
           pert%v(i,j,k) = pert%v(i,j,k) * coeffv(i,j) * coeffk(k)
        end do
     end do
  end do

  do j = jfirst, jlast
     do k = 1, nl
        do i = 1, imr
           pert%pt(i,j,k) = pert%pt(i,j,k) * coefft(i,j) * coeffk(k)
        end do
     end do
  end do

  do m = 1, size(pert%q,4)
     do j = jfirst, jlast
        do k = 1, nl
           do i = 1, imr
              pert%q(i,j,k,m) = pert%q(i,j,k,m) * coefft(i,j) * coeffk(k)
           end do
        end do
     end do
  end do

  do j = jfirst, jlast
     do k = 1, nl
        do i = 1, imr
           pert%delp(i,j,k) = pert%delp(i,j,k) * coefft(i,j) * coeffk(k)
        end do
     end do
  end do

  if(myID==ROOT) write(stdout,'(2a)') myname_, ': Applied local projection operator'

end subroutine proj_svec_

end module m_svprojs

