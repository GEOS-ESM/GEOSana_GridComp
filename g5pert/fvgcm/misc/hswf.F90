!-----------------------------------------------------------------------
!BOP
module hswf

! !USES:
      use precision

      implicit none

      integer ks
      real (r8) , allocatable :: sinp2(:),cosp2(:),cosp4(:), rf(:)
      integer    :: klow         ! FastOpt allow using klow for allocation

      public :: hswf_initialize, hswf_do, hswf_finalize

contains
  subroutine hswf_initialize(im, jm, km, pdt, rayf, sinp, cosp, sine, cose )

      implicit none

! !INPUT PARAMETERS:
      integer im, jm, km
      integer pdt
      logical rayf
      real (r8) sinp(jm), cosp(jm), sine(jm), cose(jm)

      real (r8) ak(km+1), bk(km+1)
      real (r8) ptop, pint
      real (r8) pc, c1
      integer j, k
      real (r8) tmp

        allocate( sinp2(jm), cosp2(jm), cosp4(jm), rf(km) )
        do j=2,jm-1
          sinp2(j) = sinp(j)**2
          cosp2(j) = cosp(j)**2
        enddo

        sinp2(1) = ( 0.5*(-1.+sine(2)) )**2
        sinp2(jm) = sinp2(1)
        cosp2(1) = ( 0.5*cose(2) ) **2
        cosp2(jm) = cosp2(1)

        do j=1,jm
          cosp4(j) = cosp2(j)**2
        enddo

        call set_eta(km,ks,ptop,pint,ak,bk)

        if ( rayf ) then
	  c1 = 1. / (12.*3600)
          pc = 1.
          do k=1,ks
            tmp = 0.5*(ak(k) + ak(k+1))
            rf(k) = c1*(1.+tanh(1.5*log10(pc/tmp)))
	    rf(k) = 1./(1.+pdt*rf(k))
          enddo
        endif
 
      klow = max(km/16, 1)

  end subroutine hswf_initialize


  subroutine hswf_finalize
  implicit none

  deallocate( sinp2, cosp2, cosp4, rf )
  end subroutine hswf_finalize


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: hswf --- Held-Suarez forcing
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine hswf_do(im, jm, km, jfirst, jlast,                    &
                         u,v,pt,pe,pkz,pdt,akap,grav,rg, dcaf, strat,  &
                         rayf, sinp,cosp,sine,cose, coslon, sinlon,    &
                         ng_s, ng_d)
!****6***0*********0*********0*********0*********0*********0**********72
! !USES:


#if defined( SPMD )
      use mod_comm, only :  gid, gsize, mp_barrier, mp_send_pe, mp_recv_pe, &
                            mp_send_s, mp_recv_n, mp_send_ua, mp_recv_ua
     
#endif
      implicit none

! !INPUT PARAMETERS:
      integer im, jm, km
      integer jfirst, jlast
      integer ng_s, ng_d
      integer pdt
      real(r8) akap, grav, rg
      logical strat
      logical rayf
      logical dcaf
      real(r8) cosp(jm),sinp(jm),cose(jm),sine(jm)
      real(r8) coslon(im), sinlon(im)

! !INPUT/OUTPUT PARAMETERS:

      real(r8) :: u(im,jfirst-ng_d:jlast+ng_s,km)
      real(r8) :: v(im,jfirst-ng_s:jlast+ng_d,km)
      real(r8) :: pt(im,jfirst-ng_d:jlast+ng_d,km)
      real(r8) :: pe(im,km+1,jfirst:jlast)
      real(r8) :: pkz(im,jfirst:jlast,km)


! !DESCRIPTION:
!    Author: Shian-Jiann Lin, JCET, UMBC and NASA/GSFC
!
! !REVISION HISTORY:
!   SJL 99.09.30:  Delivery
!   WS  99.10.28:  jfirst:jlast; documentation
!   WS  99.11.07:  pruned arrays
!   WS  00.07.10:  Incorporated simplfications to PILGRIM
!
!EOP
!-----------------------------------------------------------------------
!BOC
!

      real (r8) p0, t0, sday, rkv, rka, rks, rkt, sigb, rsgb
      real (r8) tmp
      real (r8) ap0k, algpk, rfc
      real (r8) tey, tez, fac, pw, sigl, tmin
      real (r8) pl(im,km+1)             ! FastOpt : fixed bug km+1
      real (r8) frac(im,jm)
      real (r8) teq(im,km)
      real (r8)  h0, dz
      real (r8) dt_tropic
      real (r8) rmr, rms
      real (r8) relx, tau
      real (r8) t_st, t_ms
      real (r8) f1
      real (r8) ptop, pint
      real (r8) ak(km+1), bk(km+1)
      real (r8) t2(im, km)
      real (r8) dp(im, km)

      real (r8) ua(im,jfirst:jlast,km)
      real (r8) va(im,jfirst:jlast,km)
      real (r8) u2(im,km)
      real (r8) v2(im,km)
      real (r8) fu(im,km)
      real (r8) fv(im,km)

#if defined( SPMD )
      real (r8) , allocatable :: pesouth(:,:), uasouth(:,:)    !will be removed later
!      real (r8) :: pesouth(im, km+1), uasouth(im, km)    !will be removed later
#endif
      real(r8) rdt
 
      integer i, j, k, js2g0, js2gm1, jn2g0, count, ierror

      js2g0  = max(2,jfirst)
      js2gm1 = max(2,jfirst+1)
      jn2g0  = min(jm-1,jlast)

      p0 = 100000.
      t0 = 200.
      h0 = 7.
      sday = 24*3600
      rkv = 0.5*pdt/sday
      rka = pdt/ (40.*sday)      ! was 40 days
      rfc = 1./(1.+rka)
      rks = pdt/ (4.*sday)       ! was 4 days

! For strat-mesosphere
      t_ms = 10.
      t_st = 40.
      tau = (t_st - t_ms) / log(100.)
      rms = pdt/(t_ms*sday)
      rmr =  1./(1.+rms)

      sigb = 0.7
      rsgb = 1./(1.-sigb)
      ap0k = 1./p0**akap
      algpk = log(ap0k)


      if (dcaf) then
!
! WS 99.10.26 : Replaced loop 500 with one call to d2a3d
!
        call d2a3d(u,v,ua,va,im,jm,km,jfirst,jlast,ng_d, ng_s, coslon,sinlon)
      endif

!$omp  parallel do                                      &
!$omp  default(shared)                                  &
!$omp  private(i,j,k,pl,tey,tez,dz,relx,dt_tropic, rdt) &
!$omp  private(teq,tmin,sigl,f1,rkt, tmp, u2, v2, fu, fv, t2, dp)

      do 1000 j=jfirst,jlast

!$taf init hswf_tape     = static, km
!$taf init dry_adj_tape  = static, im*klow

        tey = ap0k*(315.-60.*sinp2(j))
        tez = ap0k*10./akap*cosp2(j)

        do k=1,km
          do i=1,im
            pl(i,k) = 0.5*(pe(i,k,j)+pe(i,k+1,j))
          enddo
        enddo

        do i=1,im                 ! FastOpt : quick bug fix must be changed
          pl(i,km+1) = pl(i,km)   ! FastOpt : quick bug fix must be changed
        enddo                     ! FastOpt : quick bug fix must be changed

        teq = 0.                  ! FastOpt dummy initialization

        do k=km,1,-1

!$taf store teq(:,k)  = hswf_tape, key = k, kind=4
!$taf store pt(:,j,k) = hswf_tape, key = k, kind=4
                                     ! FastOpt : the last one can be avoided
          do i=1,im
            if (strat .and. pl(i,k) .lt. 10000.    &
                      .and. pl(i,k) .gt. 100.  )  then
              dz = h0 * log(pl(i,k+1)/pl(i,k))
!
! Lapse rate above tropic stratopause is 2.25 deg/km
! Relaxation time is t_st days at 100 mb (as H-S) and gradually
! decreases to t_ms Days at and above the stratopause
!
              relx =  t_ms + tau*log(0.01*pl(i,k))
              relx = pdt/(relx*sday)
              dt_tropic = 2.25*cosp(j) * dz
              teq(i,k) = (teq(i,k+1)*pkz(i,j,k+1) + &
                          dt_tropic)/pkz(i,j,k)
              pt(i,j,k) = (pt(i,j,k)+relx*teq(i,k))/(1.+relx)
!!!              pt(i,j,k) = teq(i,k)
            elseif (strat .and. pl(i,k) .le. 100.)  then
!
! Mesosphere
!
              dz = h0 * log(pl(i,k+1)/pl(i,k))
              dt_tropic = -2.25*cosp(j) * dz
              tmp = teq(i,k+1)*pkz(i,j,k+1) + dt_tropic
              teq(i,k) =  tmp / pkz(i,j,k)
!!!              teq(i,k) = max(200., tmp) / pkz(i,j,k)
              pt(i,j,k) = (pt(i,j,k)+rms*teq(i,k))*rmr
!!!              pt(i,j,k) = teq(i,k)
            else
!
! Trop:  strictly Held-Suarez
!
              sigl = pl(i,k)/pe(i,km+1,j)
              f1 = max(0._r8, (sigl-sigb) * rsgb )
              tmin = t0/pkz(i,j,k)
              teq(i,k) = tey - tez*(log(pkz(i,j,k))+algpk)
              teq(i,k) = max(tmin, teq(i,k))
              rkt = rka + (rks-rka)*f1*cosp4(j)
              pt(i,j,k) = (pt(i,j,k)+rkt*teq(i,k))/(1.+rkt)
            endif
          enddo     !i-loop
        enddo     !k-loop

! Do dry_con
        if (dcaf) then
          do k=1, km
            do i=1,im
              dp(i,k) = pe(i,k+1,j) - pe(i,k,j)
              fu(i,k) = 0.
              fv(i,k) = 0.
              u2(i,k) = ua(i,j,k)
              v2(i,k) = va(i,j,k)
              t2(i,k) = pt(i,j,k)
            enddo
          enddo

          rdt = 1.  / pdt
          call dry_adj(im, km, rdt, t2, fu, fv,   &
                       u2, v2, dp, j)

          do k=1, km
! Adjust D-grid v-winds
            v(1,j,k) = v(1,j,k) + 0.5*(fv(1,k)+fv(im,k))
            do i=2,im
              v(i,j,k) = v(i,j,k) + 0.5*(fv(i,k)+fv(i-1,k))
            enddo

            do i=1,im
              ua(i,j,k) = fu(i,k)
              pt(i,j,k) = t2(i,k)
            enddo
          enddo
        endif

1000  continue

#if defined( SPMD )
!
! Communication might include ua and/or pe on the south only

      allocate (uasouth(im, km)) 
      allocate (pesouth(im, km+1)) 

      if ( dcaf ) then
!           goto 1001
           call mp_barrier
           call mp_send_ua(im, jm, jfirst, jlast, 1, km, ua)
           call mp_barrier
           call mp_recv_ua(im, jm, jfirst, jlast, 1, km, uasouth)
           call mp_barrier
      endif
1001  continue
!      call mp_barrier 
      call mp_send_pe(im, jm, jfirst, jlast, 1, km+1, pe)
      call mp_barrier
      call mp_recv_pe(im, jm, jfirst, jlast, 1, km+1, pesouth)
      call mp_barrier

!
! Communication finished
!
#endif

!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,j,k,sigl,fac,frac)


      do 2000 k=1,km

        if (dcaf) then
          do j=js2gm1,jlast                  
            do i=1,im
              u(i,j,k) = u(i,j,k) + 0.5*(ua(i,j,k)+ua(i,j-1,k))
            enddo
          enddo
#if defined( SPMD )
          if ( jfirst .gt. 1 ) then
            do i=1,im
              u(i,jfirst,k) = u(i,jfirst,k)            &
                              + 0.5*(ua(i,jfirst,k)+uasouth(i,k))
            enddo
          endif
#endif
        endif

        if (rayf .and. k.le. ks) then
! Apply Rayleigh friction
          do j=js2g0,jlast
            do i=1,im
              u(i,j,k) = u(i,j,k)*rf(k)
            enddo
          enddo

          do j=js2g0,jn2g0
            do i=1,im
              v(i,j,k) = v(i,j,k)*rf(k)
            enddo
          enddo
        else
! Surface Rayleigh friction according to Held-Suarez
          do j=jfirst,jlast
            do i=1,im
              sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,km+1,j)
              frac(i,j) = max(0._r8, (sigl-sigb)*rsgb )
            enddo
          enddo
#if defined( SPMD )
          if ( jfirst .gt. 1 ) then
            do i=1,im
              sigl = 0.5*(pesouth(i,k)+pesouth(i,k+1)) / pesouth(i,km+1)
              frac(i,jfirst-1) = max(0._r8, (sigl-sigb)*rsgb )
            enddo
          endif
#endif

! Backward adjustment
          do j=js2g0,jlast
            do i=1,im
              fac = frac(i,j)+frac(i,j-1)
              if (fac .gt. 0.) then
                u(i,j,k) = u(i,j,k)/(1.+rkv*fac)
              endif
            enddo
          enddo

          do j=js2g0,jn2g0
            do i=2,im
              fac = frac(i,j)+frac(i-1,j)
              if (fac .gt. 0.) then
                v(i,j,k) = v(i,j,k)/(1.+rkv*fac)
              endif
            enddo
          enddo

          do j=js2g0,jn2g0
            fac = frac(1,j)+frac(im,j)
            if (fac .gt. 0.) then
              v(1,j,k) = v(1,j,k)/(1.+rkv*fac)
            endif
          enddo
        endif
2000  continue


#if defined (SPMD)
      deallocate (pesouth)
      deallocate (uasouth)
#endif

      return
!EOP
      end subroutine hswf_do
!-----------------------------------------------------------------------

end module hswf


