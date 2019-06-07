module control

! !REVISION HISTORY:
!  
!  05Feb2007  Todling  Some cleaning; add vect2mod/mod2vect expanded; to get ps
!  14May2007  Todling  Introduced dyn_prog; global change.
!  20Aug2007  Todling  Bug fix end_ps was too long
!  05Sep2007  Todling  Generalized the number of q (tracer)
!

#if defined( SPMD )
    use mod_comm,  only : gid
    use stepon, only: ng_d, ng_s
    use mod_comm, only: mp_gather4d,mp_bcst_n_real
    use mod_comm, only: mp_scatter4d
#define CPP_PRT_PREFIX  if(gid == 0)
#else
#define CPP_PRT_PREFIX
#endif

  use precision
  use prognostics, only : dyn_prog
  use prognostics, only : imr,jnp,nl,nc
  use prognostics, only : jfirst,jlast
  use prognostics, only : kfirst,klast
  
  use prognostics, only : ps

  use m_die, only : mp_die

  implicit none

  PRIVATE

! public mapvars
  public cont2mod
  public control_number 
  public mod2cont
  public mod2vect
  public vect2mod

  public ibeg_u   , iend_u
  public ibeg_v   , iend_v
  public ibeg_delp, iend_delp
  public ibeg_pt  , iend_pt
  public ibeg_q   , iend_q
  public ibeg_ps  , iend_ps
  public ibeg_tr, iend_tr
  public ncnt                 ! no. of control vars: u, v, pt, q, and delp
  public nvct                 ! no. of control vars + ps

  interface mapvars ; module procedure mapvars_    ; end interface

  integer, save :: ibeg_u   , iend_u
  integer, save :: ibeg_v   , iend_v
  integer, save :: ibeg_delp, iend_delp
  integer, save :: ibeg_pt  , iend_pt
  integer, save :: ibeg_q   , iend_q
  integer, save :: ibeg_ps  , iend_ps
  integer, save,allocatable :: ibeg_tr(:)  , iend_tr(:)
  integer, save :: ncnt
  integer, save :: nvct

  logical, save :: iamset = .false.

  real(r8), pointer ::    u(:,:,:)   ! zonal wind on D-grid
  real(r8), pointer ::    v(:,:,:)   ! meridional wind
  real(r8), pointer ::   pt(:,:,:)   ! virtual potential temperature
  real(r8), pointer :: delp(:,:,:)   ! pressure thickness (pascal)
  real(r8), pointer ::    q(:,:,:,:) ! specific humidity & tracer mixing ratios


contains
subroutine mapvars_ ( imr, jnp, nl, nc )

  integer, intent(in) :: imr
  integer, intent(in) :: jnp
  integer, intent(in) :: nl 
  integer, intent(in) :: nc 

  integer m, iend_tr0

  ibeg_u  = 1
  iend_u  = ibeg_u  + imr * jnp * nl - 1

  ibeg_v  = iend_u  + 1
  iend_v  = ibeg_v  + imr * jnp * nl - 1

  ibeg_pt = iend_v  + 1
  iend_pt = ibeg_pt + imr * jnp * nl - 1

  ibeg_q  = iend_pt + 1
  iend_q  = ibeg_q  + imr * jnp * nl - 1

  ibeg_delp = iend_q    + 1
  iend_delp = ibeg_delp + imr * jnp * nl - 1

  iend_tr0 = iend_delp
  if ( nc .gt. 1 ) then
       allocate ( ibeg_tr(nc-1), iend_tr(nc-1) )
       ibeg_tr(1) = iend_delp + 1 
       iend_tr(1) = ibeg_tr(1) + imr * jnp * nl - 1
       do m = 2, nc-1
          ibeg_tr(m) = iend_tr(m-1) + 1
          iend_tr(m) = ibeg_tr(m) + imr * jnp * nl - 1
       enddo
       iend_tr0 = iend_tr(nc-1)
  endif

  ibeg_ps = iend_tr0   + 1
  iend_ps = ibeg_ps + imr * jnp - 1

end subroutine mapvars_

subroutine cont2mod( n, x, prog )

  implicit none
  integer  :: n,m
  real(r8) :: x(n)
  type(dyn_prog), TARGET :: prog
  real(r8) :: qhelp(imr,jfirst     :jlast     ,nl)
  integer i

! Set pointers
! ------------
       u   => prog%u
       v   => prog%v
      pt   => prog%pt
      delp => prog%delp
      q    => prog%q

#ifdef SPMD

  call mp_scatter4d( x(ibeg_u)   , qhelp, imr, jnp, nl, 1, jfirst, jlast, &
                              kfirst, klast, 0   , 0   , 0 )
  u   (:,jfirst:jlast,:) = qhelp(:,jfirst:jlast,:)
   
  call mp_scatter4d( x(ibeg_v)   , qhelp, imr, jnp, nl, 1, jfirst, jlast, &
                              kfirst, klast, 0   , 0   , 0 )
  v   (:,jfirst:jlast,:) = qhelp(:,jfirst:jlast,:)
   
  call mp_scatter4d( x(ibeg_pt)  , qhelp, imr, jnp, nl, 1, jfirst, jlast, &
                              kfirst, klast, 0   , 0   , 0 )
  pt  (:,jfirst:jlast,:) = qhelp(:,jfirst:jlast,:)

  call mp_scatter4d( x(ibeg_q)   , qhelp, imr, jnp, nl, 1, jfirst, jlast, &
                              kfirst, klast, 0   , 0   , 0 )
  q   (:,jfirst:jlast,:,1) = qhelp(:,jfirst:jlast,:)

  call mp_scatter4d( x(ibeg_delp), qhelp, imr, jnp, nl, 1, jfirst, jlast, &
                              kfirst, klast, 0   , 0   , 0 )
  delp(:,jfirst:jlast,:) = qhelp(:,jfirst:jlast,:)

  if (nc .gt. 1) then
      do m = 2, nc
         call mp_scatter4d( x(ibeg_tr(m-1)), qhelp, imr, jnp, nl, 1, jfirst, jlast, &
                            kfirst, klast, 0   , 0   , 0 )
         q   (:,jfirst:jlast,:,m) = qhelp(:,jfirst:jlast,:)
      enddo
  endif

#else
  u(:,jfirst:jlast,:)   = reshape(x(ibeg_u :iend_u ),(/ imr,jnp,nl /))

  v(:,jfirst:jlast,:)   = reshape(x(ibeg_v :iend_v ),(/ imr,jnp,nl /))

  pt(:,jfirst:jlast,:)  = reshape(x(ibeg_pt:iend_pt),(/ imr,jnp,nl /))

  q(:,jfirst:jlast,:,1) = reshape(x(ibeg_q :iend_q ),(/ imr,jnp,nl /))

  delp(:,jfirst:jlast,:) = reshape(x(ibeg_delp:iend_delp),(/ imr,jnp,nl  /))
#endif

end subroutine cont2mod

subroutine vect2mod( n, x, prog )

  implicit none
  integer  :: n
  real(r8) :: x(n)
  real(r8) :: qhelp(imr,jfirst:jlast)
  type(dyn_prog), TARGET :: prog

  if ( n/=nvct ) &
       call mp_die('vect2mod',': Inconsistent dim ... aborting ',n)

  call cont2mod( n, x, prog )

#ifdef SPMD
  call mp_scatter4d( x(ibeg_ps), qhelp, imr, jnp, 1, 1, jfirst, jlast, &
                     1, 1, 0   , 0   , 0 )
  ps  (:,jfirst:jlast) = qhelp(:,jfirst:jlast)
#else
  ps(:,jfirst:jlast) = reshape(x(ibeg_ps:iend_ps),(/ imr,jnp  /))
#endif

end subroutine vect2mod

subroutine mod2cont( n, x, prog )

  implicit none
  integer  :: n, m
  real(r8) :: x(n)
  type(dyn_prog), TARGET :: prog

! Set pointers
! ------------
       u   => prog%u
       v   => prog%v
      pt   => prog%pt
      delp => prog%delp
      q    => prog%q

  !-----------------------------------------------------------------------

#ifdef SPMD
  call mp_gather4d( u   , x(ibeg_u)   , imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_d, ng_s, 0 )
  call mp_gather4d( v   , x(ibeg_v)   , imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_s, ng_d, 0 )
  call mp_gather4d( pt  , x(ibeg_pt)  , imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_d, ng_d, 0 )
  call mp_gather4d( q(:,:,:,1), x(ibeg_q), &
                    imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_d, ng_d, 0 )
  call mp_gather4d( delp, x(ibeg_delp), imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, 0   , 0   , 0 )

  if (nc .gt. 1) then
      do m = 2, nc
         call mp_gather4d( q(:,:,:,m), x(ibeg_tr(m-1)), &
                           imr, jnp, nl, 1, jfirst, jlast, &
                           kfirst, klast, ng_d, ng_d, 0 )
      enddo
  endif

! broadcast control vector
! ------------------------
! call mp_bcst_n_real( x, n )

#else
  x(ibeg_delp:iend_delp)=reshape(delp(:,jfirst:jlast,:),(/1+iend_delp-ibeg_delp/))
  x(ibeg_q :iend_q )   = reshape(q (:,jfirst:jlast,:,1),(/1+iend_q -ibeg_q /))
  x(ibeg_pt:iend_pt)   = reshape(pt(:,jfirst:jlast,:)  ,(/1+iend_pt-ibeg_pt/))
  x(ibeg_v :iend_v )   = reshape(v (:,jfirst:jlast,:)  ,(/1+iend_v -ibeg_v /))
  x(ibeg_u :iend_u )   = reshape(u (:,jfirst:jlast,:)  ,(/1+iend_u -ibeg_u /))
#endif

end subroutine mod2cont

subroutine mod2vect( n, x, prog )
  implicit none
  integer  :: n
  real(r8) :: x(n)
  type(dyn_prog), TARGET :: prog

  if ( n/=nvct ) &
       call mp_die('mod2vect','Inconsistent dim ... aborting ',n)

#ifdef SPMD
  call mp_gather4d( ps, x(ibeg_ps), imr, jnp, 1, 1, jfirst, jlast, &
                    1, 1, 0, 0, 0 )
#else
  x(ibeg_ps :iend_ps )   = reshape(ps (:,jfirst:jlast)  ,(/1+iend_ps -ibeg_ps /))
#endif
  call mod2cont( n, x, prog )
 
end subroutine mod2vect

subroutine control_number( n )
  implicit none
  integer, optional, intent(out) :: n
  integer :: myN,m

  if (iamset) then
      if(present(n)) n = ncnt
      return
  endif

  myN = 0
  myN = myN + imr * jnp * nl
  myN = myN + imr * jnp * nl
  myN = myN + imr * jnp * nl
  myN = myN + imr * jnp * nl
  myN = myN + imr * jnp * nl
  do m = 2,nc
     myN = myN + imr * jnp * nl
  enddo

  ncnt = myN
  nvct = myN + imr * jnp
  CPP_PRT_PREFIX print *,'       #control = ', ncnt
  CPP_PRT_PREFIX print *,' #ps + #control = ', nvct
  if(present(n)) n = myN

  call mapvars_ ( imr, jnp, nl, nc )

  iamset = .true.

end subroutine control_number

end module control
