!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  epvd --- Calculate absolute potential vorticity
!
! !INTERFACE:
      subroutine epvd(im, jm, km, jfirst, jlast, u, v, pt, delp,  &
                      epv, ng_s, ng_d)
! !USES:
      use precision
#if defined( SPMD )
      use mod_comm, only :  mp_send_s, mp_recv_n, mp_barrier
#endif
      implicit none

! !INPUT PARAMETERS:
      integer im, jm, km                   ! Horizontal dimensions
      integer jfirst, jlast                ! Latitude strip
      integer ng_s, ng_d

      real (r8) :: u(im,jfirst-ng_d:jlast+ng_s,km) 
      real (r8) :: v(im,jfirst-ng_s:jlast+ng_d,km) 
      real (r8) :: pt(im,jfirst-ng_d:jlast+ng_d,km) 
      real (r8) :: delp(im,jfirst:jlast,km)

! !OUTPUT PARAMETERS:
      real(r8) epv(im,jfirst:jlast,km)

! !DESCRIPTION:
!     Compute absolute vorticity on the D grid
!        epv = -g * (vort+f0)*dpt/dp
!
! !REVISION HISTORY:
!   WS  99.11.02   Documentation; indentation; jfirst:jlast
!   WS  00.07.08   Use precision module; Kevin's ghost indices
!
!EOP
!---------------------------------------------------------------------
!BOC

      real(r8) te(im,jm,km+1), t2(im,km),delp2(im,km)
      real(r8) vort(im,jfirst:jlast),fx(im,jm),fy(im,jm),te2(im,km+1)
! Geometric arrays
      real(r8) sine(jm),cosp(jm),sinp(jm)
      real(r8), allocatable, save :: f0(:), rdx(:), cy(:), cose(:)

      integer i, j, k,  js2g0, jn2g0
      real(r8) grav
      real(r8) ae, pi, dl, dp, rcap, rotsec, omega, rdy
      real(r8) c1, c2
      logical first
      data first /.true./
      save rdy, rcap

      js2g0 = max(2,jfirst)
      jn2g0 = min(jm-1,jlast)

! Geometric factors

      grav = 9.81
      if(first) then
        allocate( f0(jm), rdx(jm), cy(jm), cose(jm) )
        ae = 6.37122e6
        pi = 4. * atan(1.)
        call setrig(im,jm,dp,dl,cosp,cose,sinp,sine)
        rcap = dp / ( im*(1.+sine(2)) )
        rotsec = 86164.09
        omega  = 2.*pi / rotsec

	do j=1,jm
          f0(j)  = 2.*omega*sinp(j)
	enddo

        rdy   = 1./(ae*dp)

        do j=2,jm-1
          rdx(j) = 1./(dl*ae*cosp(j))
          cy(j) =  rdy / cosp(j)
        enddo

        first = .false.
      endif

#if defined( SPMD )
      call mp_barrier
      call mp_send_s(im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u)
#endif

! Compute PT at layer edges.

!$omp  parallel do                 &
!$omp  default(shared)             &
!$omp  private(i,j,k,t2,delp2,te2)

      do 1000 j=jfirst,jlast

        do k=1,km
          do i=1,im
	    t2(i,k) =   pt(i,j,k)
            delp2(i,k) = delp(i,j,k)
          enddo
        enddo

        call ppme(t2,te2,delp2,im,km)

        do k=1,km+1
          do i=1,im
	    te(i,j,k) = te2(i,k)
          enddo
        enddo

!     do k=2,KM
!     do i=1,IM
!       te(i,j,k) = 0.5*(PT(i,j,k)+PT(i,j,k-1))
!     enddo
!     enddo

!     do i=1,IM
!       te(i,j,   1) = PT(i,j,   1)
!       te(i,j,KM+1) = PT(i,j,KM+1)
!     enddo

1000  continue

#if defined( SPMD )
        call mp_barrier
        call mp_recv_n(im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u)
#endif

!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,j,k,fx,fy,vort,c1,c2)

      do 2000 k=1,km

! Compute relative vorticity
        do j=js2g0,min(jlast+1,jm)
          do i=1,im
            fx(i,j) = u(i,j,k)*cose(j)
          enddo
        enddo

        do j=js2g0,jn2g0
          do i=1,im-1
            fy(i,j) =  v(i+1,j,k) - v(i,j,k)
          enddo
        enddo

        do j=js2g0,jn2g0
          fy(im,j) = v(1,j,k) - v(im,j,k)
        enddo

        do j=js2g0,jn2g0
          do i=1,im
            vort(i,j) = (fx(i,j)-fx(i,j+1))*cy(j) + fy(i,j)*rdx(j)
          enddo
        enddo

! Vort at poles computed by circulation theorem

        if ( jfirst .eq. 1 ) then
          c1 = -SUM(fx(1:im,2))*rdy*rcap
          do i=1,im
            vort(i,  1) = c1
          enddo
        endif 
        if ( jlast .eq. jm )  then
          c2 = SUM(fx(1:im,jm))*rdy*rcap
          do i=1,im
            vort(i,jm) = c2
          enddo
        endif

        do j=jfirst,jlast
          do i=1,im
! Entropy is the thermodynamic variable in the following formulation.
            epv(i,j,k) = grav*(vort(i,j)+f0(j))*(te(i,j,k)-te(i,j,k+1))  &
                       / (pt(i,j,k)*delp(i,j,k))
          enddo
        enddo
2000  continue

      return
      end
