subroutine advect_cv(mype,mydate,tau,cvector)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    advect_cv    advect cv using guess winds
!   prgmmr: todling          org: np20                date: 2005-09-29
!
! abstract: advect cv using guess winds
!           both AB and RK4 are not good for time reversal (i.e., 
!           change of sign in t; will require adjoint).
!  >>       I think we need to do this using spectral decomposition
!           of the gradient operator.
!
! program history log:
!   2013-10-28  todling - stripped off from calctends
!
! usage:
!   input argument list:
!     mype     - task id
!
!   output argument list:
!     p_t      - time tendency of 3d prs
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds,only: r_kind,i_kind
  use gsi_4dvar, only: nsubwin, lsqrtb
  use gridmod, only: lat2,lon2,nsig,istart,rlats,nlat,&
     jstart,nlon,nthreads,jtstart,jtstop,bk5
  use guess_grids, only: ntguessig
  use constants, only: zero,one,two,three,half,max_varname_length
  use control_vectors, only: control_vector,cvars3d,cvars2d  
  use derivsmod, only: dvars2d,dvars3d
  use derivsmod, only: init_anadv

  use mp_compact_diffs_mod2, only: init_mp_compact_diffs2
  use mp_compact_diffs_mod2, only: destroy_mp_compact_diffs2

  use gsi_metguess_mod, only: gsi_metguess_get,gsi_metguess_bundle

  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlecreate
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_bundlemod, only: gsi_bundledestroy

  use mpeu_util, only: die
  use mpeu_util, only: getindex
  implicit none

! Declare passed variables
  integer(i_kind),intent(in) :: mype
  integer(i_kind),intent(in) :: mydate(5)
  integer(i_kind),intent(in) :: tau
  type(control_vector) :: cvector

! Declare local variables
  character(len=*), parameter :: myname = "advect_cv"

  logical, parameter :: debug = .false.
  logical, parameter :: showcv = .true.
  logical,allocatable :: adv2d(:),adv3d(:)
  character(max_varname_length) :: varname
  character(len=4) label
  integer(i_kind) i,j,k,it,ii,jj,kk,idx,n2d,n3d,nc,norder,ntimes,istatus 
  integer(i_kind) order
  real(r_kind) mytau, dt

  character(len=3) :: scheme
  integer(i_kind), parameter :: aborder = 1
! real(r_kind), parameter :: advrate = 0.75 ! pass via namelist
  real(r_kind), parameter :: advrate = 1.00 ! pass via namelist
  logical, parameter :: filter = .true.
  logical, parameter :: cnstwind_test = .false.

  type(gsi_bundle) :: f
  type(gsi_bundle) :: aux
  type(gsi_bundle), allocatable :: fn(:)
  type(gsi_bundle) :: xderivative
  type(gsi_bundle) :: yderivative

  real(r_kind), dimension(:), allocatable :: alpha,beta

  real(r_kind),dimension(:,:,:),pointer :: uwind  =>NULL()
  real(r_kind),dimension(:,:,:),pointer :: vwind  =>NULL()

  real(r_kind),dimension(:,:,:),pointer :: avptr3d=>NULL()
  real(r_kind),dimension(:,:,:),pointer :: cvptr3d=>NULL()
  real(r_kind),dimension(:,:,:),pointer :: dxptr3d=>NULL()
  real(r_kind),dimension(:,:,:),pointer :: dyptr3d=>NULL()

  real(r_kind),dimension(:,:  ),pointer :: avptr2d=>NULL()
  real(r_kind),dimension(:,:  ),pointer :: cvptr2d=>NULL()
  real(r_kind),dimension(:,:  ),pointer :: dxptr2d=>NULL()
  real(r_kind),dimension(:,:  ),pointer :: dyptr2d=>NULL()

  if (lsqrtb) then
      call die(myname,': should not be here, aborting ... ',99)
  endif
! if (tau<0) return

  scheme = "rk4"

! Mark sure derivative variables are defined
  call init_anadv
  call init_mp_compact_diffs2(nsig,mype,.false.)

  allocate(adv2d(size(dvars2d)),adv3d(size(dvars3d)))
  adv2d=.false.
  adv3d=.false.

! Initialize coefficients of Adams-Bashforth scheme
  read(scheme(3:3),*) norder
  allocate(alpha(norder),beta(norder))

! Debug (advect blob)
  if (debug) call blob_ (cvector%step(1))

! Write out input CV for testing
  if (showcv) call view_cv (cvector,mydate,'cv_before_adv',.false.)

! Create temporary bundle to hold two instance of fields at two time steps
  call gsi_bundlecreate(f,cvector%grid_step,'f',istatus,names2d=dvars2d,names3d=dvars3d)
  allocate(fn(norder))
  do jj=1,norder
     write(label,'(a,i1,a)') 'f(', jj, ')'
     call gsi_bundlecreate(fn(jj),cvector%grid_step,label,istatus,names2d=dvars2d,names3d=dvars3d)
  enddo
  if (scheme(1:2)=='rk') then
     call gsi_bundlecreate(aux,cvector%grid_step,label,istatus,names2d=dvars2d,names3d=dvars3d)
  endif

! Create vectors to hold x and y derivatives
  call gsi_bundlecreate(xderivative,cvector%grid_step,'dx',istatus,names2d=dvars2d,names3d=dvars3d)
  call gsi_bundlecreate(yderivative,cvector%grid_step,'dy',istatus,names2d=dvars2d,names3d=dvars3d)

! Get pointer for guess wind components
  call gsi_bundlegetpointer (gsi_metguess_bundle(ntguessig),'u',uwind,istatus);
  call gsi_bundlegetpointer (gsi_metguess_bundle(ntguessig),'v',vwind,istatus);

! Get dt based on some rough estimate of CFL
  mytau=tau
  mytau=1.0
  call check_cfl (dt,mytau,uwind,vwind)

  ntimes = nint(advrate*mytau*3600._r_kind / dt)
  if (mype==0) then
     write(6,'(2a,(2f7.3,i5))') myname, " tau_fcst, dt, ntimes  = ", mytau, dt, ntimes
  endif

  do jj=1,nsubwin

!    Fill auxiliar CV as needed/able
     call from_cv_(jj)

!    Loop over integration time
     do it=1,ntimes

        order = min(it,norder)
        if (order <= norder) call abcoeff_(scheme,order)

        select case(scheme(1:2))
           case ('ab') ! A class of Adams-Bashforth schemes

!            Calculate horizontal advection term at iteration it
             call hadvect_(dt,f,fn(1))
     
!            Adams-Bashforth N-th order update
             do ii=1,order
                call axpy_(alpha(ii),f,fn(ii))
             enddo
  
!            Store evaluated function (advection term) as term at time-step it-1
             if(it<ntimes) call rotate_(it,fn)

           case ('rk') ! A class of Runge-Kutta schemes

!            Calculate horizontal advection term at iteration it
             call hadvect_  (dt,f,fn(1))

             call    axpy_ (beta(2),f,fn(1),gx=aux)
             call hadvect_ (dt,aux,fn(2))

             if (norder>2) then
                call    axpy_ (beta(3),f,fn(2),gx=aux)
                call hadvect_ (dt,aux,fn(3))
             endif

             if (norder>3) then
                call    axpy_ (beta(4),f,fn(3),gx=aux)
                call hadvect_ (dt,aux,fn(4))
             endif

             do ii=1,norder
                call axpy_(alpha(ii),f,fn(ii))
             enddo

        end select

     enddo

!    Retrieve advected fields back to CV
     call to_cv_(jj)

  enddo

! Clean up
  if (scheme(1:2)=='rk') call gsi_bundledestroy(aux,istatus)
  call gsi_bundledestroy(f ,istatus)
  call gsi_bundledestroy(yderivative,istatus)
  call gsi_bundledestroy(xderivative,istatus)
  do jj=norder,1,-1
     call gsi_bundledestroy(fn(jj),istatus)
  enddo
  deallocate(fn)
  call destroy_mp_compact_diffs2
  deallocate(adv2d,adv3d)
  deallocate(alpha,beta)

! Write out advected CV
  if (showcv) call view_cv (cvector,mydate,'cv_after_adv',.false.)

  return

  contains

  subroutine from_cv_(jjj)

  character(len=*), parameter :: myname_ = myname//"*from_cv_"
  integer(i_kind), intent(in) :: jjj
  integer ier

  ier=0
  idx=getindex(dvars3d,'qi')
  call gsi_bundlegetpointer (cvector%step(jjj),'sf',cvptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (f                ,'qi',avptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  adv3d(idx)=.true.
  avptr3d=cvptr3d

  ier=0
  idx=getindex(dvars3d,'ql')
  call gsi_bundlegetpointer (cvector%step(jjj),'vp',cvptr3d,istatus);ier=ier+istatus 
  call gsi_bundlegetpointer (f                ,'ql',avptr3d,istatus);ier=ier+istatus 
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  adv3d(idx)=.true.
  avptr3d=cvptr3d

  ier=0
  idx=getindex(dvars3d,'tv')
  call gsi_bundlegetpointer (cvector%step(jjj),'t' ,cvptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (f                ,'tv',avptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  adv3d(idx)=.true.
  avptr3d=cvptr3d

  ier=0
  idx=getindex(dvars3d,'q')
  call gsi_bundlegetpointer (cvector%step(jjj),'q' ,cvptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (f                ,'q' ,avptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  adv3d(idx)=.true.
  avptr3d=cvptr3d

  ier=0
  idx=getindex(dvars2d,'ps')
  call gsi_bundlegetpointer (cvector%step(jjj),'ps',cvptr2d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (f                ,'ps',avptr2d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  adv2d(idx)=.true.
  avptr2d=cvptr2d
 
  end subroutine from_cv_

  subroutine to_cv_(jjj)
  
  integer(i_kind) jjj

  character(len=*), parameter :: myname_ = myname//"*to_cv_"
  integer ier

  ier=0
  call gsi_bundlegetpointer (f                ,'ql',avptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (cvector%step(jjj),'sf',cvptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  cvptr3d=avptr3d

  ier=0
  call gsi_bundlegetpointer (f                ,'ql',avptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (cvector%step(jjj),'vp',cvptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  cvptr3d=avptr3d

  ier=0
  call gsi_bundlegetpointer (f                ,'tv',avptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (cvector%step(jjj),'t' ,cvptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  cvptr3d=avptr3d

  ier=0
  call gsi_bundlegetpointer (f                ,'q' ,avptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (cvector%step(jjj),'q' ,cvptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  cvptr3d=avptr3d

  ier=0
  call gsi_bundlegetpointer (f                ,'ps',avptr2d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (cvector%step(jjj),'ps',cvptr2d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  cvptr2d=avptr2d
 
  end subroutine to_cv_

  subroutine hadvect_(deltime,gi,go)

  use tendsmod, only: what9,prsth9,prdif9,r_prdif9

  implicit none

  real(r_kind) :: deltime
  type(gsi_bundle) :: gi,go
  integer(i_kind) ix

  real(r_kind) :: tmp,localx,localy
  integer k1
  real(r_kind), pointer ::  func3d(:,:,:),dx3d(:,:,:),dy3d(:,:,:), x3d(:,:,:)
  real(r_kind), pointer ::  func2d(:,:  ),dx2d(:,:  ),dy2d(:,:  )

  real(r_kind),dimension(lat2,lon2,nsig):: div
  real(r_kind),dimension(lat2,lon2,nsig):: prdif9u,prdif9v
  real(r_kind),dimension(lat2,lon2,nsig):: prdif9u_x,prdif9v_y2

  do k=1,nsig
    do j=1,lon2
      do i=1,lat2
        prdif9u(i,j,k)=prdif9(i,j,k)*uwind(i,j,k)
        prdif9v(i,j,k)=prdif9(i,j,k)*vwind(i,j,k)
      end do
    end do
  end do

  call mp_compact_dlon2(prdif9u,prdif9u_x,.false.,nsig,mype)
  call mp_compact_dlat2(prdif9v,prdif9v_y2,.true.,nsig,mype)

  div(:,:,:)=zero
  do k=1,nsig
    do j=1,lon2
      do i=1,lat2
        div(i,j,k)=div(i,j,k)+prdif9u_x(i,j,k)+prdif9v_y2(i,j,k)
      end do
    end do
  end do

  do j=1,lon2
    do i=1,lat2
      prsth9(i,j,nsig+1)=zero
    end do
  end do
  do k=nsig,1,-1
    do j=1,lon2
      do i=1,lat2
         prsth9(i,j,k)=prsth9(i,j,k+1) - div(i,j,k)
      end do
    end do
  end do

  do k=2,nsig
     do j=1,lon2
       do i=1,lat2
          what9(i,j,k)=prsth9(i,j,k)-bk5(k)*prsth9(i,j,1)
       end do
     end do
  end do

! Calculate derivatives
  call get_derivatives (gi,xderivative,yderivative)
! call get_derivatives_(gi,xderivative,yderivative)

! Loop over 3d-fields
  n3d = size(dvars3d)
  do nc = 1,n3d
    if(.not.adv3d(nc)) cycle ! don''t waste time

    call gsi_bundlegetpointer (gi         ,dvars3d(nc),   x3d,istatus)
    call gsi_bundlegetpointer (go         ,dvars3d(nc),func3d,istatus)
    call gsi_bundlegetpointer (xderivative,dvars3d(nc),  dx3d,istatus)
    call gsi_bundlegetpointer (yderivative,dvars3d(nc),  dy3d,istatus)

    if (filter) then ! the better thing to do to wipe small scale stuff would be to do a 
                     ! back/forth spectral decomposition eliminating the small scales in between
       localx = maxval(abs(dx3d))
       localy = maxval(abs(dy3d))
       where(abs(dx3d)<1e-10*localx) dx3d=zero
       where(abs(dy3d)<1e-10*localy) dy3d=zero
    endif

! Loop over threads
!$omp parallel do schedule(dynamic,1) private(i,j,k,kk,k1)
   do kk=1,nthreads

    k1=1
!   Horizontal advection of "tracer" quantities
    if ( cnstwind_test ) then
       k=1
      do j=jtstart(kk),jtstop(kk)
        do i=1,lat2
           ix=istart(mype+1)+i-2
           if (ix==1 .OR. ix==nlat) then
              func3d(i,j,k) = zero
           else
              func3d(i,j,k) = - deltime * ( 20._r_kind*dx3d(i,j,k) ) ! constant zonal wind
           endif
        end do
      end do
       k=2
      do j=jtstart(kk),jtstop(kk)
        do i=1,lat2
           ix=istart(mype+1)+i-2
           if (ix==1 .OR. ix==nlat) then
              func3d(i,j,k) = zero
           else
              func3d(i,j,k) = - deltime * ( 20._r_kind*dy3d(i,j,k) ) ! constant meridional wind
           endif
        end do
      end do
      k1=3
    endif
    do k=k1,nsig
      do j=jtstart(kk),jtstop(kk)
        do i=1,lat2
           ix=istart(mype+1)+i-2
           if (ix==1 .OR. ix==nlat) then
              func3d(i,j,k) = zero
           else
              func3d(i,j,k) = - deltime * ( uwind(i,j,k)*dx3d(i,j,k) + vwind(i,j,k)*dy3d(i,j,k) )
           endif
        end do
      end do
    end do

   end do ! end of threading loop

!  call hdiff_(deltime,avptr3d,dxptr3d,dyptr3d,mype)

! vertical flux terms
  do k=1,nsig
    do j=1,lon2
      do i=1,lat2
        if (k.gt.1) then
          tmp = half*what9(i,j,k)*r_prdif9(i,j,k)
          func3d(i,j,k) = func3d(i,j,k) - tmp*(x3d(i,j,k-1)-x3d(i,j,k))
        end if
        if (k.lt.nsig) then
          tmp = half*what9(i,j,k+1)*r_prdif9(i,j,k)
          func3d(i,j,k) = func3d(i,j,k) - tmp*(x3d(i,j,k)-x3d(i,j,k+1))
        end if
      end do  !end do i
    end do    !end do j
  end do      !end do k

 end do ! end loop over 3d-fields

! Loop over 2d-fields
  n2d = size(dvars2d)
  do nc = 1,n2d
    if(.not.adv2d(nc)) cycle ! don''t waste time

    call gsi_bundlegetpointer (go         ,dvars2d(nc),func2d,istatus)
    call gsi_bundlegetpointer (xderivative,dvars2d(nc),  dx2d,istatus)
    call gsi_bundlegetpointer (yderivative,dvars2d(nc),  dy2d,istatus)

    if (filter) then ! the better thing to do to wipe small scale stuff would be to do a 
                     ! back/forth spectral decomposition eliminating the small scales in between
       localx = maxval(abs(dx2d))
       localy = maxval(abs(dy2d))
       where(abs(dx2d)<1e-10*localx) dx2d=zero
       where(abs(dy2d)<1e-10*localy) dy2d=zero
    endif

! Loop over threads
!$omp parallel do schedule(dynamic,1) private(i,j,k,kk)
   do kk=1,nthreads

!   Horizontal advection of "tracer" quantities
    k=1
    do j=jtstart(kk),jtstop(kk)
      do i=1,lat2
         ix=istart(mype+1)+i-2
         if (ix==1 .OR. ix==nlat) then
            func2d(i,j) = zero
         else
            func2d(i,j) = - deltime *( uwind(i,j,k)*dx2d(i,j) + vwind(i,j,k)*dy2d(i,j) )
          endif
      end do
    end do

   end do ! end of threading loop

 end do ! end loop over 2d-fields

 end subroutine hadvect_

 subroutine rotate_(istep,gn)
  
  integer istep
  type(gsi_bundle) gn(:)

  character(len=*), parameter :: myname_ = myname//"*rotate_"
  integer ii,order,idxi,idxo,iv,ier

  integer(i_kind),  parameter :: mynvars3d = 4
  character(len=2), parameter :: myvars3d(mynvars3d) = (/ 'qi', 'ql', 'tv', 'q ' /)
  
  order = min(istep,norder)

  do ii=order,1,-1
     idxo=min(ii,norder)
     idxi=max(ii-1,1)
     !if(mype==0) then
     !  print *, myname_, ' ', idxi, idxo
     !endif
     if(idxo==idxi) cycle

     do iv=1,mynvars3d
        ier=0
        call gsi_bundlegetpointer (gn(idxi),trim(myvars3d(iv)),avptr3d,istatus);ier=ier+istatus
        call gsi_bundlegetpointer (gn(idxo),trim(myvars3d(iv)),cvptr3d,istatus);ier=ier+istatus
        if(ier/=0) call die(myname_,'invalid pointer =',ier)
        cvptr3d=avptr3d
     enddo

     ier=0
     call gsi_bundlegetpointer (gn(idxi),'ps',avptr2d,istatus);ier=ier+istatus
     call gsi_bundlegetpointer (gn(idxo),'ps',cvptr2d,istatus);ier=ier+istatus
     if(ier/=0) call die(myname_,'invalid pointer =',ier)
     cvptr2d=avptr2d

  enddo
 
 end subroutine rotate_

 subroutine get_derivatives_(x,dx,dy)

! use mp_compact_diffs_mod2, only: mp_compact_dlon2,mp_compact_dlat2
  implicit none

  type(gsi_bundle)  x
  type(gsi_bundle)  dx
  type(gsi_bundle)  dy

  character(len=*), parameter :: myname_ = myname//"*get_derivatives_"
  character(len=20) varname
  logical, parameter :: vector=.false.
  integer ier

  real(r_kind), pointer ::  xptr3d(:,:,:) 
  real(r_kind), pointer :: dxptr3d(:,:,:) 
  real(r_kind), pointer :: dyptr3d(:,:,:) 
  real(r_kind), pointer ::  xptr2d(:,:  ) 
  real(r_kind), pointer :: dxptr2d(:,:  ) 
  real(r_kind), pointer :: dyptr2d(:,:  ) 
  real(r_kind), allocatable :: xaux3d(:,:,:) 
  real(r_kind), allocatable :: daux3d(:,:,:) 
  
  ier=0
  call gsi_bundlegetpointer (x ,'qi', xptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dx,'qi',dxptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dy,'qi',dyptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  call mp_compact_dlon2(xptr3d,dxptr3d,vector,nsig,mype)
  call mp_compact_dlat2(xptr3d,dyptr3d,vector,nsig,mype)

  ier=0
  call gsi_bundlegetpointer (x ,'ql', xptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dx,'ql',dxptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dy,'ql',dyptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid pointer =',ier)
  call mp_compact_dlon2(xptr3d,dxptr3d,vector,nsig,mype)
  call mp_compact_dlat2(xptr3d,dyptr3d,vector,nsig,mype)

  ier=0
  varname='tv'
  call gsi_bundlegetpointer (x ,trim(varname), xptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dx,trim(varname),dxptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dy,trim(varname),dyptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid '// varname // ' pointer =',ier)
  call mp_compact_dlon2(xptr3d,dxptr3d,vector,nsig,mype)
  call mp_compact_dlat2(xptr3d,dyptr3d,vector,nsig,mype)

  ier=0
  call gsi_bundlegetpointer (x ,'q', xptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dx,'q',dxptr3d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dy,'q',dyptr3d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid (q) pointer =',ier)
  call mp_compact_dlon2(xptr3d,dxptr3d,vector,nsig,mype)
  call mp_compact_dlat2(xptr3d,dyptr3d,vector,nsig,mype)

  ier=0
  call gsi_bundlegetpointer (x ,'ps', xptr2d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dx,'ps',dxptr2d,istatus);ier=ier+istatus
  call gsi_bundlegetpointer (dy,'ps',dyptr2d,istatus);ier=ier+istatus
  if(ier/=0) call die(myname_,'invalid (ps)  =',ier)
  allocate(xaux3d(size(xptr2d,1),size(xptr2d,2),nsig))
  allocate(daux3d(size(xptr2d,1),size(xptr2d,2),nsig))
  xaux3d=zero
  xaux3d(:,:,1) =  xptr2d
  call mp_compact_dlon2(xaux3d,daux3d,vector,nsig,mype)
  dxptr2d = daux3d(:,:,1)
  call mp_compact_dlat2(xaux3d,daux3d,vector,nsig,mype)
  dyptr2d = daux3d(:,:,1)
  deallocate(daux3d)
  deallocate(xaux3d)

 end subroutine get_derivatives_

 subroutine abcoeff_(sch,n)
 implicit none
 character(len=*), intent(in):: sch
 integer(i_kind) :: n

 character(len=*), parameter :: myname_=myname//"*abcoeff"

 select case(n)
   case(1)
      alpha(1) = one
      beta(1)  = one
   case(2)
      if (trim(sch)=='rk') then ! Heun''s Method
         alpha(1) =  half
         alpha(2) =  half
         beta(1)  =  zero
         beta(2)  =  one
      else
         alpha(1) =  three/two
         alpha(2) = -half
      endif
   case(3)
      alpha(1) =  23./12._r_kind
      alpha(2) = -16./12._r_kind
      alpha(3) =   5./12._r_kind
   case(4)
      if (trim(sch)=='rk') then
         alpha(1) =  one/6._r_kind
         alpha(2) =  two/6._r_kind
         alpha(3) =  two/6._r_kind
         alpha(4) =  one/6._r_kind
         beta(1)  =  zero
         beta(2)  =  half
         beta(3)  =  half
         beta(4)  =  one
      else
         alpha(1) =  55./24._r_kind
         alpha(2) = -59./24._r_kind
         alpha(3) =  37./24._r_kind
         alpha(4) =  -9./24._r_kind
      endif
   case(5)
      alpha(1) =  1901./720._r_kind
      alpha(2) = -2774./720._r_kind
      alpha(3) =  2616./720._r_kind
      alpha(4) = -1274./720._r_kind
      alpha(5) =   251./720._r_kind
   case(6)
      alpha(1) =  4277./1440._r_kind
      alpha(2) = -7923./1440._r_kind
      alpha(3) =  9982./1440._r_kind
      alpha(4) = -7298./1440._r_kind
      alpha(5) =  2877./1440._r_kind
      alpha(6) =  -475./1440._r_kind
   case default
      call die(myname_,'Unacceptable Adams-Bashforth order, aborting',99)   
   end select

 end subroutine abcoeff_

 subroutine axpy_(coef,g,gn,gx)

! Adams–Bashforth

 real(r_kind)    :: coef
 type(gsi_bundle) g
 type(gsi_bundle) gn
 type(gsi_bundle),optional :: gx

 real(r_kind), pointer, dimension(:,:,:) :: ptr3d,ptr3dn,ptr3dx
 real(r_kind), pointer, dimension(:,:  ) :: ptr2d,ptr2dn,ptr2dx
 integer(i_kind) ii

! Loop over 3d-fields
  n3d = size(dvars3d)
  do nc = 1,n3d
    if(.not.adv3d(nc)) cycle ! don''t waste time

    call gsi_bundlegetpointer (g ,dvars3d(nc),ptr3d ,istatus)
    call gsi_bundlegetpointer (gn,dvars3d(nc),ptr3dn,istatus)
    if (present(gx)) then
       call gsi_bundlegetpointer (gx,dvars3d(nc),ptr3dx,istatus)
       ptr3dx = ptr3d + coef*ptr3dn
    else
       ptr3d  = ptr3d + coef*ptr3dn
    endif

 enddo

! Loop over 2d-fields
  n2d = size(dvars2d)
  do nc = 1,n2d
    if(.not.adv2d(nc)) cycle ! don''t waste time

    call gsi_bundlegetpointer (g ,dvars2d(nc),ptr2d ,istatus)
    call gsi_bundlegetpointer (gn,dvars2d(nc),ptr2dn,istatus)
    if (present(gx)) then
       call gsi_bundlegetpointer (gx,dvars2d(nc),ptr2dx,istatus)
       ptr2dx = ptr2d + coef*ptr2dn
    else
       ptr2d  = ptr2d + coef*ptr2dn
    endif

 enddo

 end subroutine axpy_

 subroutine blob_(blob)

! set a blob on each subdomain

 type(gsi_bundle) blob

 real(r_kind),pointer :: ptr2d(:,:)
 integer(i_kind), parameter :: imany=4 ! choose an even number
 integer(i_kind) imid,jmid

 call gsi_bundlegetpointer (blob,'ps',ptr2d   ,istatus)

 imid = nint(lat2/2.)-imany/2
 jmid = nint(lon2/2.)-imany/2

 ptr2d = zero
 if (mod(mype,8)==0) then
    ptr2d(imid:imid+imany,jmid:jmid+imany) = one
 endif

 end subroutine blob_

 subroutine handle_poles_(xx,mype)
 use constants, only: zero
 use gridmod, only: lat2,istart,nlat
 implicit none
 real(r_kind), pointer :: xx(:,:,:)
 integer(i_kind) mp1,ix,mype

!mp1=mype+1
!if(nlat+1-istart(mp1)+2==lat2) then
!   xx(lat2,:,:)=zero
!end if
!if(2-istart(mm1)==1) then
!   xx(1,:,:)=zero
!end if
    do k=1,nsig
      do j=1,lon2
        do i=1,lat2
          ix=istart(mype+1)+i-2
          if (ix == 1 .OR. ix == nlat) then
            xx(i,j,k)=zero
          end if
        end do
      end do
    end do  !end do k
 end subroutine handle_poles_
 
 subroutine hdiff_(time_step,q_t,q_x,q_y,mype)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    hdiff           calculate horiziontal diffusion
!   prgmmr: rancic                                    
!
! abstract: compute horizontal diffusion passive tracer
!
! program history log:
!   2010-02-25  rancic
!   2020-11-20  todling - adapt to case in point 
!
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind  
  use constants, only: zero,half,one,two
  use gridmod, only: lat2,lon2,nsig
  implicit none

! Declare passed variables
  integer(i_kind)                       ,intent(in   ) :: mype
  real(r_kind)                          ,intent(in   ) :: time_step
  real(r_kind),dimension(lat2,lon2,nsig),intent(in   ) :: q_x,q_y
  real(r_kind),dimension(lat2,lon2,nsig),intent(inout) :: q_t

! Declare local variables
  real(r_kind),dimension(lat2,lon2,nsig):: qc_x,qc_y
  real(r_kind),dimension(lat2,lon2,nsig):: q_xx,q_yy
  real(r_kind) uave,vave,up,vp
  real(r_kind) khdff,kht,facdec
  integer(i_kind) i,j,k

  khdff=1.73e04_r_kind
  kht=khdff*time_step

  q_xx=zero; q_yy=zero

  do k=1,nsig
    if(k<nsig-2) then
      facdec=exp(-3.*(k-1)/(nsig-1))   
    else
      facdec=one
    end if
    facdec=facdec*kht
    do j=1,lon2
    do i=1,lat2
      qc_x(i,j,k)=q_x(i,j,k)*facdec
      qc_y(i,j,k)=q_y(i,j,k)*facdec
    end do
    end do
  end do

  call mp_compact_dlon2(qc_x,q_xx,.false.,nsig,mype)
  call mp_compact_dlat2(qc_y,q_yy,.false.,nsig,mype)     
      
  do k=1,nsig
    do j=1,lon2
    do i=1,lat2
       q_t(i,j,k)=q_t(i,j,k)+q_xx(i,j,k)+q_yy(i,j,k)
    end do
    end do
  end do

 end subroutine hdiff_

 subroutine hdiff_ad_(time_step,q_t,q_x,q_y,mype)
  use kinds,only: r_kind,i_kind  
  use constants, only: zero,half,one,two
  use gridmod, only: lat2,lon2,nsig
  implicit none

! Declare passed variables
  integer(i_kind)                       ,intent(in   ) :: mype
  real(r_kind)                          ,intent(in   ) :: time_step
  real(r_kind),dimension(lat2,lon2,nsig),intent(  out) :: q_x,q_y
  real(r_kind),dimension(lat2,lon2,nsig),intent(inout) :: q_t

! Declare local variables
  integer(i_kind) i,j,k
  real(r_kind),dimension(lat2,lon2,nsig):: qc_x,qc_y
  real(r_kind),dimension(lat2,lon2,nsig):: q_xx,q_yy
  real(r_kind) khdff,kht,facdec
!
! Preliminaries
!

  khdff=1.73e04_r_kind
  kht=khdff*time_step

!
! Start adjoint
!
  q_xx=q_t   ;   q_yy=q_t
  qc_x=zero  ;   qc_y=zero

  call mp_compact_dlon2_ad(qc_x,q_xx,.false.,nsig,mype)  
  call mp_compact_dlat2_ad(qc_y,q_yy,.false.,nsig,mype)     


  do k=nsig,1,-1
    if(k<nsig-2) then
      facdec=exp(-3.*(k-1)/(nsig-1))   
    else
      facdec=one
    end if
    facdec=facdec*kht
    do j=1,lon2
    do i=1,lat2
      q_x(i,j,k)=qc_x(i,j,k)*facdec
      q_y(i,j,k)=qc_y(i,j,k)*facdec
    end do
    end do
  end do

 end subroutine hdiff_ad_

 subroutine check_cfl(deltat,fcstlen,u,v) ! rough CFL estimate
 use constants, only: ten,rearth
 use mpimod, only: mpi_comm_world,mpi_rtype,mpi_max
 use guess_grids, only: get_ref_gesprs
 use gridmod, only: nlon,nlat
 implicit none
! RT: needless to say this is not really CFL since CFL varies from cell to cell
! (depending in dlon,dlat on spherical coordinates) and getting impossible at 
! the poles
 real(r_kind), intent(inout) :: deltat
 real(r_kind), intent(in) :: fcstlen
 real(r_kind), intent(in) :: u(:,:,:),v(:,:,:)
 real(r_kind), allocatable:: ple(:)
 real(r_kind) dt0,wfc,mindt,dlon,dlat,umax,vmax,cfl,maxu,maxv
 integer k,ier
 dlon=rearth/nlon
 dlat=rearth/(nlat-1)
 allocate(ple(size(u,3)+1))
 call get_ref_gesprs(ple)
 deltat = 180.0
 dt0 = deltat
 if(mype==0) then
   print *, 'on way in (dt): ', deltat
   print *
   write(6,'(3a)') '    prs     cfl      dt' 
   write(6,'( a)') '-----------------------'
 endif
 do k=1,size(u,3)
    maxu = maxval( abs(u(:,:,k)) )
    maxv = maxval( abs(v(:,:,k)) )
    call mpi_allreduce(maxu,umax,1,mpi_rtype,mpi_max,mpi_comm_world,ier)
    call mpi_allreduce(maxv,vmax,1,mpi_rtype,mpi_max,mpi_comm_world,ier)
    wfc = (umax/dlon+vmax/dlat)
    cfl = deltat*wfc
    mindt = one/wfc
    if(mype==0) then
       write(6,'(3(f7.2,1x))')  ten*ple(k), cfl, mindt
    endif
    if(mindt<deltat) deltat=mindt
 enddo
 if (deltat<dt0) then
    deltat = 0.8*deltat
    if(mype==0) print *, 'on way out (dt): ', deltat
    ntimes = nint(fcstlen*3600._r_kind/deltat)
    if(mype==0) print *, 'on way out (ntimes): ', ntimes
    deltat = fcstlen*3600_r_kind/ntimes
 endif
 if(mype==0) print *, 'final (dt): ', deltat
 deallocate(ple)
 end subroutine check_cfl

end subroutine advect_cv
