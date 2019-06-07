!-----------------------------------------------------------------------
! FastOpt:
!   splitted into module and module procedures
!   changed cy to cyy to avoid name crash
!-----------------------------------------------------------------------

!BOP
module cd_core

! !USES:
   use precision

      implicit none

! Declare permanent local arrays
      integer  ifax(13)                      !ECMWF fft
      real(r8), allocatable :: trigs(:)
      real(r8), allocatable :: dc(:,:), de(:,:), sc(:), se(:)
      real(r8), allocatable :: fc(:), f0(:)
      real(r8), allocatable :: cdx(:,:), cdy(:,:)
      real(r8), allocatable :: dtdx(:), dtdxe(:), txe5(:), dtxe5(:)
      real(r8), allocatable :: dyce(:),   dx(:) ,  rdx(:),    cyy(:)
      real(r8), allocatable :: dtdx2(:), dtdx4(:),  dxdt(:), dxe(:)
      real(r8), allocatable :: cye(:),    dycp(:),  rdxe(:)

      real(r8) rdy, dtdy, dydt, dtdy5, tdy5

      real(r8) zt_c   
      real(r8) zt_d  
      real(r8) dt0, dt5

      integer :: cd_core_tape_rec         ! FastOpt

!
!EOP
!---------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: cd_core_allocate --- Dynamical core allocation
!
! !INTERFACE:
  subroutine cd_core_allocate( im,   jm,   km,               &
                               jfirst, jlast,                &
                               ng_c, ng_d, ng_s              &
                             )

! !USES:
   use precision
   use sw_core
!  use tp_core

   implicit none

! Input paraterers:
  integer im, jm, km
  integer jfirst
  integer jlast
  integer ng_c
  integer ng_d               ! Max NS dependencies
  integer ng_s

!EOP
!---------------------------------------------------------------------
!BOC
! Local
  integer js1g1, js2g0, js2g1, js2gc
  integer jn2g0, jn1g1, jn1gc
  integer ktot, ktotp

!-----------------------------------------------------------------------

! Set general loop limits
! jfirst >= 1; jlast <= jm
      js1g1 = max(1,jfirst-1)
      js2g0 = max(2,jfirst)
      js2g1 = max(2,jfirst-1)
      jn2g0 = min(jm-1,jlast)
      jn1g1 = min(jm,jlast+1)

! Construct C-grid dependent loop limits
      js2gc  = max(2,jfirst-ng_c)     ! NG latitudes on S (starting at 2)
!
! WS: 00.04.13 : An ugly hack to handle JORD>1 JCD=1
!
      if ( ng_c == 1 .AND. ng_d  > 1 ) THEN
        js2gc  = max(2,jfirst-2) 
      endif

      jn1gc  = min(jm,jlast+ng_c)     ! NG latitudes on N (ending at jm)

      allocate(dtdx(jm),dtdx2(jm), dtdx4(jm), dtdxe(jm), dxdt(jm),    &
                   dxe(jm),  cye(jm),dycp(jm),rdxe(jm),                   &
                   txe5(jm), dtxe5(jm),dyce(jm),                          &
                   dx(jm),rdx(jm),cyy(jm) )

      allocate( trigs(3*im/2+1) )
      allocate( sc(js2g0:jn2g0),    se(js2g0:jn1g1)    )
      allocate( dc(im,js2g0:jn2g0), de(im,js2g0:jn1g1) )

      allocate( cdx(js2g0:jn1g1,km) )
      allocate( cdy(js2g0:jn1g1,km) )

      allocate( f0(jfirst-ng_s:jlast+ng_d) ) ! 000304 bug fix: ng_s not ng_d
      allocate( fc(js2gc:jn1gc) )


end subroutine cd_core_allocate


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: cd_core_deallocate --- Dynamical core deallocation
!
! !INTERFACE:
  subroutine cd_core_deallocate

! !USES:
   use precision
   use sw_core
!  use tp_core

   implicit none

!EOP
!---------------------------------------------------------------------
!BOC

!-----------------------------------------------------------------------
  deallocate(dtdx,dtdx2, dtdx4, dtdxe, dxdt, &
           dxe,  cye,dycp,rdxe,              &
           txe5, dtxe5,dyce,                 &
           dx,rdx,cyy )

  deallocate( trigs )
  deallocate( sc, se )
  deallocate( dc, de )

  deallocate( cdx )
  deallocate( cdy )

  deallocate( f0 )
  deallocate( fc )

end subroutine cd_core_deallocate


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: cd_core_initialize --- Dynamical core initialization
!
! !INTERFACE:
  subroutine cd_core_initialize( im,   jm,   km,                &
                                 jfirst, jlast,                 &
                                 ng_c, ng_d, ng_s,              &
                                 dt, ae, om,  ptop, umax,       &
                                 sinp, cosp, cose, acosp, akap  &
                               )

! !USES:
   use precision
   use sw_core
!  use tp_core

#if defined( SPMD )
   use mod_comm, only: gid

#define CPP_PRT_PREFIX  if(gid==0)
#else
#define CPP_PRT_PREFIX
#endif

   implicit none

! Input paraterers:
  integer im, jm, km
  integer jfirst
  integer jlast
  integer ng_c
  integer ng_d               ! Max NS dependencies
  integer ng_s
  real(r8) dt            !small time step in seconds
  real(r8) ae                ! Radius of the Earth (m)
  real(r8) om                ! rotation rate
  real(r8) ptop
  real(r8) umax
  real(r8) cosp(jm)
  real(r8) cose(jm)
  real(r8) sinp(jm)
  real(r8) acosp(jm)
  real(r8) akap

!EOP
!---------------------------------------------------------------------
!BOC
! Local
      integer i, j, k
      integer ks
      integer js1g1, js2g0, js2g1, js2gc
      integer jn2g0, jn1g1, jn1gc
      integer iord , jord
      integer ktot, ktotp

      double precision pi
      real(r8)  rat, ycrit
      real(r8)  dt5

      real(r8) ak(km+1), bk(km+1)
      real(r8) ptmp, pint, press
      real(r8)  tau, fac

      real(r8) p1, p2
      real(r8) dl, dp

!-----------------------------------------------------------------------
! Set general loop limits
! jfirst >= 1; jlast <= jm
      js1g1 = max(1,jfirst-1)
      js2g0 = max(2,jfirst)
      js2g1 = max(2,jfirst-1)
      jn2g0 = min(jm-1,jlast)
      jn1g1 = min(jm,jlast+1)

! Construct C-grid dependent loop limits
      js2gc  = max(2,jfirst-ng_c)     ! NG latitudes on S (starting at 2)
!
! WS: 00.04.13 : An ugly hack to handle JORD>1 JCD=1
!
      if ( ng_c == 1 .AND. ng_d  > 1 ) THEN
        js2gc  = max(2,jfirst-2) 
      endif

      jn1gc  = min(jm,jlast+ng_c)     ! NG latitudes on N (ending at jm)

      call fftfax(im, ifax, trigs)

! Determine ycrit such that effective DX >= DY
      pi  = 4.d0 * datan(1.d0)
      rat = float(im)/float(2*(jm-1))
      ycrit = acos( min(0.81_r8, rat) ) * (180./pi)

      call pft_cf(im, jm, js2g0, jn2g0, jn1g1, sc, se, dc, de,       &
                  cosp, cose, ycrit)

!     do j=1,jm
      do j=max(1,jfirst-ng_s),min(jm,jlast+ng_d)        ! 000304 bug fix
         f0(j) = (om+om)*sinp(j)
      enddo

! Compute coriolis parameter at cell corners.
!     do j=2,jm
      do j=js2gc, jn1gc
         fc(j) = 0.5*(f0(j) + f0(j-1))
      enddo

        dt0 = dt
        dt5 = 0.5*dt

        pi = 4.d0 * datan(1.d0)
        dl = (pi+pi)/im
        dp = pi/(jm-1)

        rdy   = 1./(ae*dp)
        dtdy  = dt *rdy
        dtdy5 = dt5*rdy
        dydt  = (ae*dp) / dt
        tdy5  = 0.5/dtdy

        do j=2,jm-1
          dx(j)    = dl*ae*cosp(j)
          rdx(j)   = 1./dx(j)
          dtdx(j)  = dt / dx(j)
          dxdt(j)  = dx(j) / dt
          dtdx2(j) = 0.5*dtdx(j)
          dtdx4(j) = 0.5*dtdx2(j)
          dycp(j)  = ae*dp/cosp(j)
          cyy(j)   =  rdy * acosp(j)
        enddo

        do j=2,jm
          dxe(j)   = ae*dl*cose(j)
          rdxe(j)  = 1. / dxe(j)
          dtdxe(j) = dt / dxe(j)
          dtxe5(j) = 0.5*dtdxe(j)
          txe5(j)  = 0.5/dtdxe(j)
           cye(j)  =  1. / (ae*cose(j)*dp)
          dyce(j)  = ae*dp/cose(j)
        enddo

! C-grid
        zt_c = abs(umax*dt5) / (dl*ae)
!       CPP_PRT_PREFIX write(6,*) 'c-core: ', (180./pi)*acos(zt_c)
! D-grid
        zt_d = abs(umax*dt) / (dl*ae)
!       CPP_PRT_PREFIX write(6,*) 'd-coret: ', (180./pi)*acos(zt_d)

        call set_eta( km, ks, ptmp, pint, ak, bk )

        if ( ptop /= ptmp) then
             write(6,*) 'PTOP as input to cd_core != ptop from set_eat'
             stop
        endif

!---------------------------------------
! Divergence damping coeff. unit: m**2/s
!---------------------------------------
 
!_RT        CPP_PRT_PREFIX write(6,*) 'Divergence coefficient (m**2/sec):'

          do k=1,km
             press = 0.5 * ( ak(k)+ak(k+1) + (bk(k)+bk(k+1))*1.E5 )
             tau = 8. * (1.+ tanh(1.0*log(max(ptop,2._r8)/press)) )
! The last factor of the following eqn is to adjust for higher resolution --
             tau = max(1._r8, tau) / (128.* max(abs(dt), 225._r8))
            do j=js2g0,jn1g1
              fac = tau * ae / cose(j)
              cdx(j,k) = fac*dp
              cdy(j,k) = fac*dl
            enddo
          enddo

end subroutine cd_core_initialize




!-----------------------------------------------------------------------
!BOP
! !ROUTINE: cd_core_do --- Dynamical core for both C- and D-grid Lagrangian
!                          dynamics; tracer advection is optional
!                          Default: no tracer advection with small-time-step
!
! !INTERFACE:

 subroutine cd_core_do(im,  jm,  km,  nq,  nx,                     &
                       jfirst, jlast,    u,  v,  pt,               &
                       delp,    pe,    pk,  ns,                    &
                       dt, ptop , umax  , fill, filter,            &
                       acap,  ae,  rcap,  cp,  akap, iord_c,       &
                       jord_c, iord_d, jord_d, ng_c, ng_d, ng_s,   &
                       ipe, om, hs, sinp, cosp, cose, acosp,       &
                       sinlon, coslon, cosl5, sinl5,               &
                       cx3  , cy3, mfx, mfy,                       &
                       delpf, uc, vc, ptc, dpt, ptk,               &
                       wz3, pkc, wz )

! !USES:

   use precision
   use sw_core
!  use tp_core

   use timingModule

#if defined( SPMD )
   use mod_comm, only: mp_send3d_ns, mp_recv3d_ns,                     &
                       mp_send3d_ns2, mp_recv3d_ns2,                   &
                       mp_send2_n, mp_recv2_s, mp_send_s, mp_recv_n,   &
                       mp_send2_ns, mp_recv2_ns, gid, gsize, mp_barrier

#define CPP_PRT_PREFIX  if(gid==0)
#else
#define CPP_PRT_PREFIX
#endif

   implicit none

! !INPUT PARAMETERS:

! Input paraterers:
      integer im, jm, km
      integer nq
      integer ns                 ! # of time splits
      integer nx                 ! # of split pieces in longitude direction
      integer jfirst
      integer jlast
      integer ipe                ! ipe=1:  end of cd_core()
                                 ! ipe=-1: start of cd_core()
                                 ! ipe=0 :
      real(r8) ae                ! Radius of the Earth (m)
      real(r8) om                ! rotation rate

      real(r8) ptop
      real(r8) umax

      real(r8) dt            !small time step in seconds
      real(r8) acap
      real(r8) rcap
      real(r8) cp
      real(r8) akap

      logical fill
      logical filter

      integer iord_c, jord_c
      integer iord_d, jord_d
      integer ng_c
      integer ng_d            ! Max NS dependencies
      integer ng_s

! Input time independent arrays:
      real(r8) hs(im,jfirst:jlast)   !surface geopotential
      real(r8) sinp(jm)
      real(r8) cosp(jm)
      real(r8) acosp(jm)
      real(r8) cose(jm)

      real(r8) sinlon(im)
      real(r8) coslon(im)
      real(r8) sinl5(im)
      real(r8) cosl5(im)

! !INPUT/OUTPUT PARAMETERS:
      real(r8)  u(im,jfirst-ng_d:jlast+ng_s,km)      ! u-Wind (m/s)
      real(r8)  v(im,jfirst-ng_s:jlast+ng_d,km)      ! v-wind (m/s)
      real(r8) pt(im,jfirst-ng_d:jlast+ng_d,km)      ! Potential temperature
      real(r8) delp(im,jfirst:jlast,km)        ! Delta pressure

! Input/output: accumulated winds & mass fluxes on c-grid for large-
!               time-step transport
      real(r8) cx3(im,jfirst-ng_d:jlast+ng_d,km)!Accumulated Courant number in X
      real(r8) cy3(im,jfirst:jlast+1,km)        !Accumulated Courant number in Y
      real(r8) mfx(im,jfirst:jlast,km)          !Mass flux in X  (unghosted)
      real(r8) mfy(im,jfirst:jlast+1,km)        !Mass flux in Y

! !OUTPUT PARAMETERS:
      real(r8) pe(im,km+1,jfirst:jlast)         ! Edge pressure
      real(r8) pk(im,jfirst:jlast,km+1)         ! Pressure to the kappa

! Input work arrays:
      real(r8) delpf(im,jfirst-ng_d:jlast+ng_d,km)
      real(r8)   uc(im,jfirst-ng_d:jlast+ng_d,km)
      real(r8)   vc(im,jfirst-2:   jlast+2,   km)

      real(r8) ptc(im,jfirst:jlast,km)
      real(r8) ptk(im,jfirst:jlast,km)

      real(r8) dpt(im,jfirst-1:jlast+1,km)
      real(r8) wz3(im,jfirst-1:jlast  ,km+1)
      real(r8) pkc(im,jfirst-1:jlast+1,km+1) 
      real(r8) wz (im,jfirst-1:jlast+1,km+1)

#if defined (HIGH_P)
      real(r8) wzz(im,jfirst:jlast+1,km+1)
#endif

! ! !DESCRIPTION:
!    Perform a dynamical update for one small time step; the small
!    time step is limitted by the fastest wave within the Lagrangian control-
!    volume 
!
! !REVISION HISTORY:
!     SJL  99.01.01:   Original SMP version
!     WS   99.04.13:   Added jfirst:jlast concept
!     SJL  99.07.15:   Merged c_core and d_core to this routine
!     WS   99.09.07:   Restructuring, cleaning, documentation
!     WS   99.10.18:   Walkthrough corrections; frozen for 1.0.7
!     WS   99.11.23:   Pruning of some 2-D arrays
!     SJL  99.12.23:   More comments; general optimization; reduction
!                      of redundant computation & communication
!     WS   00.05.14:   Modified ghost indices per Kevin's definition
!     WS   00.07.13:   Changed PILGRIM API
!
!EOP
!---------------------------------------------------------------------
!BOC
! Local 2D arrays:
      real(r8)  wk(im,jfirst:jlast+2)
      real(r8) wk2(im,jfirst:jlast+1)

      real(r8) wk1(im,jfirst-1:jlast+1)
      real(r8) wk3(im,jfirst-1:jlast+1)

      real(r8) pkchelp(im,jfirst-1:jlast+1)      ! FastOpt

! Local 1D
      real(r8) p1d(im)
      real(r8) ak(km+1), bk(km+1)
      integer ks
      real(r8) ptmp, pint, press
      real(r8) p1, p2

!Local scalars
      real(r8) rat, ycrit

      integer i, j, k
      integer js1g1, js2g0, js2g1, js2gc
      integer jn2g0, jn1g1, jn1gc
      integer iord , jord

      double precision pi
      real(r8) tau, fac, pk4
      real(r8) tiny
      parameter (tiny = 1.e-10)

      integer  c_sw_tape_rec      ! FastOpt
      integer  d_sw_tape_rec      ! FastOpt
      integer, parameter :: store_kind = 8   ! FastOpt precision of tape

! Set general loop limits
! jfirst >= 1; jlast <= jm
      js1g1 = max(1,jfirst-1)
      js2g0 = max(2,jfirst)
      js2g1 = max(2,jfirst-1)
      jn2g0 = min(jm-1,jlast)
      jn1g1 = min(jm,jlast+1)

! Construct C-grid dependent loop limits
      js2gc  = max(2,jfirst-ng_c)     ! NG latitudes on S (starting at 2)
!
! WS: 00.04.13 : An ugly hack to handle JORD>1 JCD=1
!
      if ( ng_c == 1 .AND. ng_d  > 1 ) THEN
        js2gc  = max(2,jfirst-2) 
      endif

      jn1gc  = min(jm,jlast+ng_c)     ! NG latitudes on N (ending at jm)


!FastOpt begin
	wz = 0.
	wz3 = 0.
	uc = 0.
	vc = 0.
!Fastopt end

#if defined( SPMD )
!-------------------------------
! Send (u, v) on the way
!-------------------------------
    call timing_on('Send_uv')
      call mp_barrier
      call mp_send3d_ns( im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u,  1)
      call mp_send3d_ns2(im, jm, jfirst, jlast, 1, km, ng_s, ng_d, v,  2)
    call timing_off('Send_uv')
#endif

      if ( ipe == -1 .or. ns == 1 ) then          ! starting cd_core
!$taf loop = parallel
!$omp parallel do private(i, j, k, wk, wk2)
      do k=1,km
         do j=jfirst,jlast
            do i=1,im
               delpf(i,j,k) = delp(i,j,k)
            enddo
         enddo
       call pft2d(delpf(1,js2g0,k),  sc(js2g0), dc(1,js2g0),     &
                  im, jn2g0-js2g0+1, ifax, trigs, wk, wk2)
      enddo
      endif

#if defined( SPMD )
!-------------------------------
! Receive (u, v)
!-------------------------------
      call timing_on('Recv_uv')
        call mp_barrier
        call mp_recv3d_ns( im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u,  1) 
        call mp_recv3d_ns2(im, jm, jfirst, jlast, 1, km, ng_s, ng_d, v,  2) 
      call timing_off('Recv_uv')

      if ( ipe == -1 .or. ns == 1 ) then          ! starting cd_core
!-------------------------------
! Send/Recv pt and delpf
!-------------------------------
        call timing_on('Ghost_pt&delpf')
          call mp_barrier
          call mp_send3d_ns( im, jm, jfirst, jlast, 1, km, ng_d, ng_d, pt,    1)
          call mp_send3d_ns2(im, jm, jfirst, jlast, 1, km, ng_d, ng_d, delpf, 2)

          call mp_barrier
          call mp_recv3d_ns( im, jm, jfirst, jlast, 1, km, ng_d, ng_d, pt,    1) 
          call mp_recv3d_ns2(im, jm, jfirst, jlast, 1, km, ng_d, ng_d, delpf, 2) 
        call timing_off('Ghost_pt&delpf')
       endif   
#endif

       call timing_on('C_CORE')

!$taf store u,v            = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind
!$taf store delpf,delp,pt  = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind

!$omp parallel do default(shared) private(k, iord, jord, c_sw_tape_rec)

      do k=1,km

        if ( k < km/8 ) then
             iord = 1
             jord = 1
        else
             iord = iord_c
             jord = jord_c
        endif

!-----------------------------------------------------------------
! Call the vertical independent part of the dynamics on the C-grid
!-----------------------------------------------------------------
	c_sw_tape_rec =  k-1 + cd_core_tape_rec*km

        call c_sw( u(1,jfirst-ng_d,k),   v(1,jfirst-ng_s,k),         &
                  pt(1,jfirst-ng_d,k),                               &
                  delp(1,jfirst,k),       uc(1,jfirst-ng_d,k),       &
                  vc(1,jfirst-2,k),       ptc(1,jfirst,k),           &
                  delpf(1,jfirst-ng_d,k), ptk(1:im,jfirst:jlast,k),  &
                  cosp,   acosp,   cose,   coslon,   sinlon,         &
                  dxdt,   dxe,      dtdx2,                           &
                  dtdx4,  dtxe5,   rdxe,   dycp,     dydt, dtdy5,    &
                  cye,    fc,      ifax,   trigs,    dc(1,js2g0),    &
                  sc,     zt_c,    tiny,   rcap,     im,             &
                  jm,     jfirst,  jlast,  ng_c,     ng_d,           &
                  ng_s,   js2g0,  jn2g0,   js2gc,    jn1gc,          &
                  iord,   jord,   cosl5,   sinl5,                    &
		  c_sw_tape_rec                                      &
		  )
      enddo

      call timing_off('C_CORE')

!$taf store ptc,ptk = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind

! Need to ghost for D_CORE
! uc(jfirst-ng_d:jlast+ng_d)
! vc(jfirst     :jlast+1   )

    call timing_on('C_GEOP')
    call geopk(ptop, pe, ptk, pkc, wz, hs, ptc, im, jm, km,   &
               jfirst, jlast, 0, cp, akap, nx, 0, .false.)
    call timing_off('C_GEOP')

#if defined( SPMD )
    call timing_on('Send_pkc_wz')
    call mp_send2_n(im, jm, jfirst, jlast, 1, km+1, 1, 1, pkc, wz)
    call timing_off('Send_pkc_wz')
#endif

!$taf store pkc,wz = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind

#if defined ( cHIGH_P )
     call timing_on('HIGHP')
        call highp2(pkc,  wz,  wz3,  wzz,  dpt, im,             &
                    jm,   km,  jfirst, jlast, nx )
     call timing_off('HIGHP')

     call timing_on('C_U_LOOP')

!$omp parallel do private(i, j, k, p1d, wk, wk2)
      do k=1,km
         do j=js2g0,jn2g0
             do i=1,im
                p1d(i) = pkc(i,j,k+1) - pkc(i,j,k)
             enddo

             uc(1,j,k) = uc(1,j,k) + dtdx2(j) / (p1d(1)+p1d(im)) *    &
                (dpt(im,j,k)-dpt(1,j,k)-wz3(1,j,k)+wz3(1,j,k+1))
            do i=2,im
              uc(i,j,k) = uc(i,j,k) + dtdx2(j) / (p1d(i)+p1d(i-1)) *  &
                 (dpt(i-1,j,k)-dpt(i,j,k)-wz3(i,j,k)+wz3(i,j,k+1))
            enddo
         enddo          ! j-loop

         call pft2d(uc(1,js2g0,k), sc(js2g0), dc(1,js2g0), im,       &
                    jn2g0-js2g0+1, ifax,   trigs, wk, wk2 )
      enddo                 ! end k paralle loop  
    call timing_off('C_U_LOOP')
#else

    call timing_on('C_U_LOOP')

!$omp parallel do private(i, j, k, p1d, wk, wk2)
      do k=1,km
         do j=js2g0,jn2g0
             do i=1,im
                p1d(i) = pkc(i,j,k+1) - pkc(i,j,k)
             enddo
            uc(1,j,k) = uc(1,j,k) + dtdx2(j) * (                      &
                (wz(im,j,k+1)-wz(1,j,k))*(pkc(1,j,k+1)-pkc(im,j,k))   &
              + (wz(im,j,k)-wz(1,j,k+1))*(pkc(im,j,k+1)-pkc(1,j,k)))  &
                / (p1d(1)+p1d(im))

            do i=2,im
              uc(i,j,k) = uc(i,j,k) + dtdx2(j) * (                     &
              (wz(i-1,j,k+1)-wz(i,j,k))*(pkc(i,j,k+1)-pkc(i-1,j,k))    &
            + (wz(i-1,j,k)-wz(i,j,k+1))*(pkc(i-1,j,k+1)-pkc(i,j,k)) )  &
                                 / (p1d(i)+p1d(i-1))
            enddo
         enddo

         call pft2d(uc(1,js2g0,k), sc(js2g0), dc(1,js2g0), im,        &
                    jn2g0-js2g0+1,  ifax,   trigs, wk, wk2)
      enddo                 ! end k paralle loop  
    call timing_off('C_U_LOOP')
#endif

#if defined( SPMD )
    call timing_on('Recv_pkc_wz')
      call mp_barrier
      call mp_recv2_s(im, jm, jfirst, jlast, 1, km+1, 1, 1, pkc, wz)
    call timing_off('Recv_pkc_wz')
!  Send uc on the way
    call timing_on('Send_uc')
      call mp_send3d_ns2(im, jm, jfirst, jlast, 1, km, ng_d, ng_d, uc, 2)
    call timing_off('Send_uc')
#endif

!$taf store pkc,wz = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind

    call timing_on('C_V_LOOP')

!$omp parallel do private(i, j, k, wk, wk1)
      do 2500  k=1,km
         do j=js1g1,jlast
             do i=1,im
                wk1(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
             enddo
         enddo

         do j=js2g0,jlast
            do i=1,im
              vc(i,j,k) = vc(i,j,k) + dtdy5/(wk1(i,j)+wk1(i,j-1)) *   &
#if defined (cHIGH_P)
              ( dpt(i,j-1,k)-dpt(i,j,k)-wzz(i,j,k)+wzz(i,j,k+1))
#else
            ((wz(i,j-1,k+1)-wz(i,j,k))*(pkc(i,j,k+1)-pkc(i,j-1,k))    &
           + (wz(i,j-1,k)-wz(i,j,k+1))*(pkc(i,j-1,k+1)-pkc(i,j,k)))
#endif
          enddo
        enddo

          call pft2d(vc(1,js2g0,k), se(js2g0), de(1,js2g0), im,      &
                     jlast-js2g0+1,  ifax, trigs, wk, wk1 )
2500  continue

    call timing_off('C_V_LOOP')

#if defined( SPMD )
!------------
! Receive uc
!------------
    call timing_on('Recv_uc')
      call mp_barrier
      call mp_recv3d_ns2(im, jm, jfirst, jlast, 1, km, ng_d, ng_d, uc, 2) 
    call timing_off('Recv_uc')

!------------
! Send/recv vc
!------------
    call timing_on('Ghost_vc')
      call mp_send_s(im, jm, jfirst, jlast, 1, km, 2, 2, vc)
      call mp_barrier
      call mp_recv_n(im, jm, jfirst, jlast, 1, km, 2, 2, vc)
    call timing_off('Ghost_vc')
#endif

! Construct D-grid dependent loop limits

    call timing_on('D_CORE')

!$taf store u,v,uc,vc,pt   = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind

!$omp parallel do default(shared) private(k, iord, jord, d_sw_tape_rec)

    do k=1,km

       if( k < km/8 ) then
           if( k == 1 ) then
              iord = 1
              jord = 1
           else
              iord = min(2, iord_d)
              jord = min(2, jord_d)
           endif
       else
          iord = iord_d
          jord = jord_d
       endif

!-----------------------------------------------------------------
! Call the vertical independent part of the dynamics on the D-grid
!-----------------------------------------------------------------
     d_sw_tape_rec = k-1 + cd_core_tape_rec*km

     call d_sw( u(1,jfirst-ng_d,k),      v(1,jfirst-ng_s,k),     &
                uc(1,jfirst-ng_d,k),    vc(1,jfirst-2,k),        &
                pt(1,jfirst-ng_d,k),   delp(1,jfirst,k),         &
                delpf(1,jfirst-ng_d,k), cx3(1,jfirst-ng_d,k),    &
                cy3(1,jfirst,k),        mfx(1,jfirst,k),         &
                mfy(1,jfirst,k), cdx(js2g0,k),  cdy(js2g0,k),    &
                dtdx,   dtdxe,  dtxe5,  txe5,  dyce,  rdx,  cyy, &
                dx,  f0(jfirst-ng_d), js2g0,  jn1g1, im,  jm,    &
                jfirst, jlast,  ng_d,  ng_s,   nq,    iord,      &
                jord,   zt_d,   rcap,  tiny,   dtdy,             &
                dtdy5,  tdy5,   rdy,    cosp,  acosp, cose,      &
                coslon, sinlon, cosl5, sinl5,                    &
		d_sw_tape_rec                                    &
		)
    enddo

    call timing_off('D_CORE')

!$taf store delp,pt = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind

    call timing_on('D_GEOP')
    call geopk(ptop, pe, delp, pkc, wz, hs, pt, im, jm, km,     &
               jfirst, jlast, ng_d, cp, akap, nx, ipe, .true.)
    call timing_off('D_GEOP')

#if defined( SPMD )
!-------------------------------------
!  Send pkc and wz NS boundary regions
!-------------------------------------
    call timing_on('Send_pkc_wz')
      call mp_send2_ns(im, jm, jfirst, jlast, 1, km+1, 1, pkc, wz)
    call timing_off('Send_pkc_wz')
#endif

!$taf store pkc,wz = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind

    if ( ipe /= 1 ) then          !  not the last call
!$omp parallel do private(i, j, k, wk, wk2)
      do k=1,km
         do j=jfirst,jlast
            do i=1,im
               delpf(i,j,k) = delp(i,j,k)
             enddo
          enddo
        call pft2d(delpf(1,js2g0,k), sc(js2g0), dc(1,js2g0),   &
                   im, jn2g0-js2g0+1, ifax, trigs, wk, wk2)
      enddo

    else

! Last call
!$omp parallel do private(i, j, k)
      do k=1,km+1
        do j=jfirst,jlast
           do i=1,im
             pk(i,j,k) = pkc(i,j,k)
           enddo
        enddo
      enddo
    endif

#if defined( SPMD )

! Recv pkc  & wz
    call timing_on('Recv_pkc_wz')
      call mp_barrier
      call mp_recv2_ns(im, jm, jfirst, jlast, 1, km+1, 1, pkc, wz)
    call timing_off('Recv_pkc_wz')

!$taf store pkc,wz = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind

    if ( ipe /= 1 ) then          !  not the last call
      call timing_on('Send_pt&delpf')
        call mp_send3d_ns( im, jm, jfirst, jlast, 1, km, ng_d, ng_d, pt   , 1)
        call mp_send3d_ns2(im, jm, jfirst, jlast, 1, km, ng_d, ng_d, delpf, 2)
      call timing_off('Send_pt&delpf')
    endif
#endif

#if defined (HIGH_P)
    call timing_on('HIGHP')
      call highp2(pkc,  wz,  wz3,  wzz,  dpt, im,         &
                  jm,   km,  jfirst, jlast, nx )
    call timing_off('HIGHP')
#else

!$omp parallel do private(i, j, k)
      do k=1,km
        do j=js1g1,jn1g1                  ! dpt needed NS
          do i=1,im                       ! wz, pkc ghosted NS
            dpt(i,j,k)=(wz(i,j,k+1)+wz(i,j,k))*(pkc(i,j,k+1)-pkc(i,j,k))
          enddo
        enddo
      enddo
#endif

      call timing_on('D-4500')

!$omp parallel do private(i, j, k, wk1, wk3, pkchelp)

#if defined (HIGH_P)
      do 4500 k=1,km+1
#else
      do 4500 k=2,km+1
#endif

!!!     do j=2,jm-1
        do j=js2g1,jn2g0
#if defined (HIGH_P)
          do i=1,im
             wk3(i,j) = wz3(i,j,k) 
          enddo
#else
! i=1

            wk3(1,j) = (wz(1,j,k)+wz(im,j,k)) *        &
                     (pkc(1,j,k)-pkc(im,j,k))
          do i=2,im
            wk3(i,j) = (wz(i,j,k)+wz(i-1,j,k)) *       &
                       (pkc(i,j,k)-pkc(i-1,j,k))
          enddo
#endif
        enddo

      do j=js2g1,jn2g0                       ! wk1 needed S
         do i=1,im-1
            wk1(i,j) = wk3(i,j) + wk3(i+1,j)
         enddo
            wk1(im,j) = wk3(im,j) + wk3(1,j)
      enddo

      if ( jfirst == 1 ) then
          do i=1,im
            wk1(i,  1) = 0.
          enddo
      endif

      if ( jlast == jm ) then
          do i=1,im
            wk1(i,jm) = 0.
          enddo
      endif

!!!     do j=2,jm
        do j=js2g0,jlast                          ! wk1 ghosted S
          do i=1,im
            wz3(i,j,k) = wk1(i,j) + wk1(i,j-1)
          enddo
        enddo

! N-S walls

        do j=js2g0,jn1g1                         ! wk1 needed N
#if defined (HIGH_P)
          do i=1,im                            ! wz, pkc ghosted NS
             wk1(i,j) = wzz(i,j,k)
          enddo
#else
          do i=1,im                            ! wz, pkc ghosted NS
             wk1(i,j) = (wz(i,j,k)+wz(i,j-1,k))*(pkc(i,j,k)-pkc(i,j-1,k))
          enddo
#endif

        enddo

        do j=js2g0,jn1g1                         ! wk3 needed N
            wk3(1,j) = wk1(1,j) + wk1(im,j)      ! wk1 ghosted N
          do i=2,im
            wk3(i,j) = wk1(i,j) + wk1(i-1,j)   ! wk1 ghosted N
          enddo
        enddo

        do j=js2g0,jn2g0
          do i=1,im
            wz(i,j,k) = wk3(i,j) + wk3(i,j+1)  ! wk3 ghosted N
          enddo
        enddo

! FastOpt begin
        do j=jfirst-1,jlast+1
          do i=1,im
            pkchelp(i,j) = pkc(i,j,k)
          enddo
        enddo

!        call avgc( pkc(1,jfirst-1,k), pkc(1,jfirst,k), im, jm,       &
!                  jfirst, jlast, wk1)
        call avgc( pkchelp, pkc(1,jfirst,k), im, jm,       &
                  jfirst, jlast, wk1)
! FastOpt end

4500  continue

      call timing_off('D-4500')

#if ( !defined HIGH_P )
      do j=js2g0,jlast
        do i=1,im
          wz3(i,j,1) = 0.
           wz(i,j,1) = 0.
        enddo
      enddo

          pk4 = 4.*ptop**akap
      do j=js2g0,jn1g1
        do i=1,im
          pkc(i,j,1) = pk4
        enddo
      enddo
#endif

      call timing_on('D-6000')

!$taf store pkc,wz,wz3 = cd_core_tape, key=cd_core_tape_rec+1, kind=store_kind

!$omp parallel do private(i, j, k, wk, wk2, wk1, wk3)
      do 6000 k=1,km

        call avgc(dpt(1,jfirst-1,k), wk2(1,jfirst),im,jm,       &
                  jfirst,jlast,wk1)

!!!     do j=2,jm
        do j=js2g0,jn1g1
          do i=1,im
             wk(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
          enddo
        enddo

        do j=js2g0,jlast
          do i=1,im-1
            wk3(i,j) = uc(i,j,k) + dtdxe(j)/(wk(i,j) + wk(i+1,j))  &
                    * (wk2(i,j)-wk2(i+1,j)+wz3(i,j,k+1)-wz3(i,j,k))
          enddo
          wk3(im,j) = uc(im,j,k) + dtdxe(j)/(wk(im,j) + wk(1,j))   &
                     * (wk2(im,j)-wk2(1,j)+wz3(im,j,k+1)-wz3(im,j,k))
        enddo

        do j=js2g0,jn2g0                  ! Assumes wk2 ghosted on N
          do i=1,im
            wk1(i,j) = vc(i,j,k) + dtdy/(wk(i,j)+wk(i,j+1)) *     &
                      (wk2(i,j)-wk2(i,j+1)+wz(i,j,k+1)-wz(i,j,k))
          enddo
        enddo

#if ( !defined ALT_PFT )
          call pft2d(wk3(1,js2g0), se(js2g0), de(1,js2g0), im,     &
                     jlast-js2g0+1,  ifax, trigs, wk, wk2 )
          call pft2d(wk1(1,js2g0), sc(js2g0), dc(1,js2g0), im,     &
                     jn2g0-js2g0+1,  ifax, trigs, wk, wk2 )
#endif

         do j=js2g0,jlast
           do i=1,im
             u(i,j,k) = u(i,j,k) + wk3(i,j)
           enddo
         enddo

         do j=js2g0,jn2g0
           do i=1,im
             v(i,j,k) = v(i,j,k) + wk1(i,j)
           enddo
         enddo

#if defined ( ALT_PFT )
       call pft2d(u(1,js2g0,k), se(js2g0), de(1,js2g0), im,       &
                  jlast-js2g0+1,  ifax, trigs, wk, wk2)
       call pft2d(v(1,js2g0,k), sc(js2g0), dc(1,js2g0), im,       &
                  jn2g0-js2g0+1,  ifax, trigs, wk, wk2)
#endif

6000  continue
      call timing_off('D-6000')

#if defined( SPMD )
      if ( ipe /= 1 ) then
        call timing_on('Recv_pt&delpf')
        call mp_barrier
        call mp_recv3d_ns( im, jm, jfirst, jlast, 1, km, ng_d, ng_d, pt   , 1)
        call mp_recv3d_ns2(im, jm, jfirst, jlast, 1, km, ng_d, ng_d, delpf, 2)
        call timing_off('Recv_pt&delpf')
      endif
#endif
   return
!EOC
 end subroutine cd_core_do
!-----------------------------------------------------------------------


end module cd_core
!-----------------------------------------------------------------------
