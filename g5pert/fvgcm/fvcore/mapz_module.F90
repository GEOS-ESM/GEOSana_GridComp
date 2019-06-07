
! FastOpt separated initialization and finalization
! FastOpt make variable diag a parameter
! FastOpt set dtmp in else part
! FastOpt renamed cosp to te_cosp to avoid name conflict
! FastOpt renamed acap to te_acap to avoid name conflict
! FastOpt renamed ptop to te_ptop to avoid name conflict

module mapz_module

   use precision

   implicit none

! Geometric arrays
      real(r8) , allocatable :: te_cosp(:), ak(:), bk(:)
      real(r8) dp, dl, te_acap

      real(r8) te_ptop
      integer ks

#if defined( SPMD )
      real(r8), allocatable :: pesouth(:,:)
#endif

 public te_map_initialize, te_map, te_map_finalize

contains

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: te_map_initialize --- Allocation and initialization for te_map
!
! !INTERFACE:

subroutine te_map_initialize( im, jm, km )
   implicit none

! !INPUT PARAMETERS:
   integer im, jm, km            ! x, y, z dimensions

!EOP
!-----------------------------------------------------------------------
!BOC
! Local arrays:
      real(r8) pint
      real(r8) sine(jm)
      real(r8) sinp(jm)
      real(r8) cose(jm)

      allocate(te_cosp(jm),ak(km+1),bk(km+1))
#if defined( SPMD )
      allocate( pesouth(im,km+1) )
#endif

      call setrig(im, jm, dp, dl, te_cosp, cose, sinp, sine)
      te_acap = im*(1.+sine(2))/dp
      call set_eta(km, ks, te_ptop, pint, ak, bk)

end subroutine te_map_initialize
!-----------------------------------------------------------------------
 

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: te_map_finalize --- Deallocation
!
! !INTERFACE:

subroutine te_map_finalize

   implicit none

!EOP
!-----------------------------------------------------------------------
!BOC

   deallocate( te_cosp, ak, bk )
#if defined( SPMD )
   deallocate( pesouth )
#endif

end subroutine te_map_finalize

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: te_map --- Map vertical Lagrangian coordinates to normal grid
!
! !INTERFACE:
 subroutine te_map(consv, convt, ps, omga, pe, delp, pkz, pk, mdt,   &
                   im, jm, km, nx, jfirst, jlast, nq,  u,  v,        &
                   pt, q, hs, cp, akap, kord, peln, te0,             &
                   ng_d, ng_s, te_map_tape_rec )
!
! !USES:

  use precision

#if defined( SPMD )
      use mod_comm, only: mp_send_pe, mp_recv_pe, mp_send_s, mp_recv_n, mp_barrier
#endif

      implicit none

! !INPUT PARAMETERS:
      logical consv                 ! flag to force TE conservation
      logical convt                 ! flag to control pt output (see below)
      integer mdt                   ! mapping time step (same as phys)
      integer im, jm, km            ! x, y, z dimensions
      integer nq                    ! number of tracers (including h2o)
      integer nx                    ! number of SMP "decomposition" in x
      integer ng_d
      integer ng_s
      integer jfirst, jlast         ! starting & ending latitude index
      real(r8) hs(im,jfirst:jlast)  ! surface geopotential
      real(r8) cp
      real(r8) te0
      integer  te_map_tape_rec           ! FastOpt

! !INPUT/OUTPUT PARAMETERS:
      real(r8) pk(im,jfirst:jlast,km+1) ! pe to the kappa
      real(r8) q(im,jfirst-ng_d:jlast+ng_d,km,nq)
      real(r8) delp(im,jfirst:jlast,km) ! pressure thickness
      real(r8) pe(im,km+1,jfirst:jlast) ! pressure at layer edges
      real(r8) ps(im,jfirst:jlast)      ! surface pressure

      real(r8) u(im,jfirst-ng_d:jlast+ng_s,km)   ! u-wind (m/s)
      real(r8) v(im,jfirst-ng_s:jlast+ng_d,km)   ! v-wind (m/s)
      real(r8) pt(im,jfirst-ng_d:jlast+ng_d,km)  ! virtual potential temperature as input
                                    ! Output: virtual temperature if convt is true
                                    ! false: output is (virtual) potential temperature 
 
! !OUTPUT PARAMETERS:
      real(r8) omga(im,km,jfirst:jlast)    ! vertical press. velocity (pascal/sec)
      real(r8) peln(im,km+1,jfirst:jlast)  ! log(pe)
      real(r8) pkz(im,jfirst:jlast,km)     ! layer-mean pk for converting t to pt

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! WS 99.05.19 : Replaced IMR, JMR, JNP and NL with IM, JM-1, JM and KM
! WS 99.05.25 : Revised conversions with IMR and JMR; removed fvcore.h
! WS 99.07.29 : Reduced computation region to jfirst:jlast
! WS 99.07.30 : Tested concept with message concatenation of te_map calls
! WS 99.10.01 : Documentation; indentation; cleaning
! SJL 99.12.31: SMP "decomposition" in-E-W direction
! WS 00.05.14 : Renamed ghost indices as per Kevin's definitions
! WS 00.07.13 : Changed PILGRIM API
!
!EOP
!-----------------------------------------------------------------------
!BOC
! Local arrays:
      real(r8) te(im,jfirst:jlast,km)  ! Work array (cache performance)
      real(r8) dz(im,jfirst:jlast,km)  ! Work array (cache performance)
      real(r8) rmin(nx*jm), rmax(nx*jm)
      real(r8) tte(jfirst:jlast)
! x-y
      real(r8)  u2(im,jfirst:jlast+1)
      real(r8)  v2(im,jfirst:jlast)
      real(r8)  t2(im,jfirst:jlast)
! y-z
      real(r8)  pe0(im,km+1)
      real(r8)  pe1(im,km+1)
      real(r8)  pe2(im,km+1)
      real(r8)  pe3(im,km+1)
      real(r8) phis(im,km+1)
! x
      real(r8)     gz(im)
      real(r8)  ratio(im)
      real(r8)    bte(im)
! z
      real(r8) pe1w(km+1)
      real(r8) pe2w(km+1)

      integer i, j, k, js2g0, jn2g0
      integer kord
      integer krd

      real(r8) akap, dak, bkh, qmax, qmin
      real(r8) te_sp,  te_np
      real(r8) xsum, ysum, tsum
      real(r8) dtmp
      real(r8) rdt5
      real(r8) rg
      real(r8) tvm
      real(r8) te1
      real(r8) dlnp

      integer ixj, jp, it, i1, i2

      logical, parameter  :: diag = .false.       ! FastOpt
      integer, parameter :: store_kind = 8        ! FastOpt precision of tape

      js2g0  = max(2,jfirst)
      jn2g0  = min(jm-1,jlast)

#if defined( SPMD )
      call mp_send_s(im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u)
#endif
      call pkez(nx, im, km, jfirst, jlast,             &
                pe, pk, akap, ks, peln, pkz, .false.)
 
!$taf store pk,pkz    = te_map_tape, key=te_map_tape_rec+1, kind=store_kind
!$taf store dz,te     = te_map_tape, key=te_map_tape_rec+1, kind=store_kind
!$taf store peln      = te_map_tape, key=te_map_tape_rec+1, kind=store_kind

#if defined( SPMD )
      call mp_barrier
      call mp_recv_n(im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u)
      call mp_send_pe(im, jm, jfirst, jlast, 1, km+1, pe)
#endif

!$taf store u,v       = te_map_tape, key=te_map_tape_rec+1, kind=store_kind

!$omp parallel do private(i, j, k, u2, v2, t2, te_sp, te_np)
      do 1000 k=1,km
! Compute cp*T + KE
        do j=js2g0,min(jlast+1,jm)
           do i=1,im
              u2(i,j) = u(i,j,k)**2
           enddo
        enddo

        do j=js2g0,jn2g0
          do i=1,im
             v2(i,j) = v(i,j,k)**2
          enddo
        enddo

        do j=jfirst,jlast
           do i=1,im
              t2(i,j) = cp*pt(i,j,k)
           enddo
        enddo

        do j=js2g0,jn2g0
          do i=1,im-1
            te(i,j,k) = 0.25 * ( u2(i,j) + u2(i,j+1) +       &
                                 v2(i,j) + v2(i+1,j)  ) +    &
                        t2(i,j)*pkz(i,j,k)
          enddo
! i=im
            te(im,j,k) = 0.25 * ( u2(im,j) + u2(im,j+1) +    &
                                v2(im,j) + v2(1,j)  ) +      &
                         t2(im,j)*pkz(im,j,k)
        enddo

        if( jfirst == 1 ) then
            te_sp = 0.
          do i=1,im
            te_sp = te_sp + u2(i,  2) + v2(i,  2)
          enddo
            te_sp = 0.5*te_sp/float(im) + t2(1,1)*pkz(1,1,k)
          do i=1,im
            te(i,1,k) = te_sp
          enddo
        endif

        if ( jlast == jm ) then
            te_np = 0.
          do i=1,im
            te_np = te_np + u2(i,jm) + v2(i,jm-1)
          enddo
            te_np = 0.5*te_np/float(im) + t2(1,jm)*pkz(1,jm,k)
          do i=1,im
            te(i,jm,k) = te_np
          enddo
        endif

! Compute dz; geo-potential increments
        do j=jfirst,jlast
           do i=1,im
              dz(i,j,k) = t2(i,j)*(pk(i,j,k+1)-pk(i,j,k))
           enddo
        enddo
1000  continue

!$taf store dz,te = te_map_tape, key=te_map_tape_rec+1, kind=store_kind

#if defined( SPMD )
      call mp_barrier
      call mp_recv_pe(im, jm, jfirst, jlast, 1, km+1, pesouth)
#endif
      it = im / nx
      jp = nx * ( jlast - jfirst + 1 )

!$omp parallel do default(shared)                &
!$omp private(i,j,k,pe0,pe1,pe2,pe3,ratio)       &
!$omp private(dak,bkh,rdt5,phis,krd,ixj,i1,i2)   &
!$omp private( pe1w, pe2w )

!     do 2000 j=jfirst,jlast
      do 2000 ixj=1, jp

         j  = jfirst + (ixj-1) / nx
         i1 = 1 + it * mod(ixj-1, nx)
         i2 = i1 + it - 1

! Copy data to local 2D arrays.
        do k=1,km+1
           do i=i1,i2
              pe1(i,k) = pe(i,k,j)
           enddo

! Ghosting for v mapping
           if( i1 == 1 ) then
               pe1w(k) = pe(im,k,j)
           else
               pe1w(k) = pe(i1-1,k,j)
           endif
        enddo

        do k=1,ks+1
           do i=i1,i2
              pe0(i,k) = ak(k)
              pe2(i,k) = ak(k)
              pe3(i,k) = ak(k)
            enddo
        enddo

        do k=ks+2,km
           do i=i1,i2
              pe0(i,k) = ak(k) + bk(k)* ps(i,j)
              pe2(i,k) = ak(k) + bk(k)*pe1(i,km+1)
           enddo
        enddo

        do i=i1,i2
           pe0(i,km+1) =  ps(i,j)
           pe2(i,km+1) = pe1(i,km+1)
        enddo

! Ghosting for v mapping
        do k=ks+2,km
           pe2w(k) = ak(k) + bk(k)*pe1w(km+1)
         enddo
           pe2w(km+1) = pe1w(km+1)

! Compute omga (dp/dt)
          rdt5 = 0.5 / float(mdt)
        do k=2,km+1
           do i=i1,i2
              pe0(i,k) = pe1(i,k) - pe0(i,k)
           enddo
        enddo

        do i=i1,i2
! update ps
            ps(i,j)   = pe1(i,km+1)
          omga(i,1,j) = rdt5 * pe0(i,2)
        enddo

        do k=2,km
          do i=i1,i2
             omga(i,k,j) = rdt5 * ( pe0(i,k) + pe0(i,k+1) )
          enddo
        enddo

        if(ks /= 0) then
           do k=1,ks
             dak = ak(k+1) - ak(k)
             do i=i1,i2
                delp(i,j,k) = dak
             enddo
           enddo
        endif

        do k=ks+1,km
          do i=i1,i2
             delp(i,j,k) = pe2(i,k+1) - pe2(i,k)
          enddo
        enddo

! Compute correction terms to Total Energy
        do i=i1,i2
           phis(i,km+1) = hs(i,j)      
        enddo

        do k=km,1,-1
          do i=i1,i2
             phis(i,k) = phis(i,k+1) + dz(i,j,k)   
          enddo
        enddo

        do k=1,km+1
          do i=i1,i2
             phis(i,k) = phis(i,k) * pe1(i,k)
          enddo
        enddo

! <<< Compute Total Energy >>>
        do k=1,km
          do i=i1,i2
            te(i,j,k) =  te(i,j,k) + (phis(i,k+1) - phis(i,k)) /   &
                         (pe1(i,k+1) - pe1(i,k) )
          enddo
        enddo

! Map Total Energy
        call map1_ppm ( km,   pe1,   te,                           &
                        km,   pe2,   te,  0,  0,                   &
                        im, i1, i2, j, jfirst, jlast, 1, kord)

! Map constituents

       if(nq /= 0) then
          if(kord == 8) then
             krd = 8
          else
             krd = 7
          endif

          call mapn_ppm ( km,   pe1,   q, nq,                     &
                          km,   pe2,   q, ng_d, ng_d,             &
                          im, i1, i2, j, jfirst, jlast, 0, krd)
       endif

! map u
        if(j /= 1) then

! WS 99.07.29 : protect j==jfirst case
          if (j > jfirst) then
            do k=2,km+1
              do i=i1,i2
                pe0(i,k) = 0.5*(pe1(i,k)+pe(i,k,j-1))
              enddo
            enddo

            do k=ks+2,km+1
              bkh = 0.5*bk(k)
              do i=i1,i2
                pe3(i,k) = ak(k) + bkh*(pe1(i,km+1)+pe(i,km+1,j-1))
              enddo
            enddo

#if defined( SPMD )
          else
            do k=2,km+1
              do i=i1,i2
                pe0(i,k) = 0.5*(pe1(i,k)+pesouth(i,k))
              enddo
            enddo

            do k=ks+2,km+1
              bkh = 0.5*bk(k)
              do i=i1,i2
                pe3(i,k) = ak(k) + bkh*(pe1(i,km+1)+pesouth(i,km+1))
              enddo
            enddo
#endif
          endif

          call map1_ppm( km,   pe0,   u,                         &
                         km,   pe3,   u, ng_d, ng_s,             &
                         im, i1, i2, j, jfirst, jlast, -1, kord)
        endif

! map v
        if(j .ne. 1 .and. j .ne. jm) then
          do k=2,km+1
! pe1(i1-1,1:km+1) must be ghosted
               pe0(i1,k) = 0.5*(pe1(i1,k)+pe1w(k))
            do i=i1+1,i2
               pe0(i ,k) = 0.5*(pe1(i,k)+pe1(i-1,k))
            enddo
          enddo

          do k=ks+2,km+1
! pe2(i1-1,ks+2:km+1) must be ghosted
               pe3(i1,k) = 0.5*(pe2(i1,k)+pe2w(k))
            do i=i1+1,i2
               pe3(i,k) = 0.5*(pe2(i,k)+pe2(i-1,k))
            enddo
          enddo

          call map1_ppm ( km,   pe0,   v,                       &
                          km,   pe3,   v, ng_s, ng_d,           &
                          im, i1, i2, j, jfirst, jlast, -1, kord)
        endif

! Save new PE to temp storage peln
        do k=2,km
          do i=i1,i2
             peln(i,k,j) = pe2(i,k)
          enddo
        enddo

! Check deformation.
       if( diag ) then
          rmax(ixj) = 0.
          rmin(ixj) = 1.
          do k=1,km
             do i=i1,i2
              ratio(i) = (pe1(i,k+1)-pe1(i,k)) / (pe2(i,k+1)-pe2(i,k))
             enddo

             do i=i1,i2
              if(ratio(i) > rmax(ixj)) then
                 rmax(ixj) = ratio(i)
              elseif(ratio(i) .lt. rmin(ixj)) then
                 rmin(ixj) = ratio(i)
              endif
            enddo
          enddo
       endif
2000  continue

#if defined( SPMD )
      call mp_send_s(im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u)
#endif

!$taf store v,te = te_map_tape, key=te_map_tape_rec+1, kind=store_kind

      if( diag ) then
            qmin = rmin(1)
        do ixj=2, jp
          if(rmin(ixj) .lt. qmin) then
            qmin = rmin(ixj)
          endif
        enddo
        write(6,*) 'rmin=', qmin

            qmax = rmax(1)
        do ixj=2, jp
          if(rmax(ixj) > qmax) then
            qmax = rmax(ixj)
          endif
        enddo
        write(6,*) 'rmax=', qmax
      endif

!$omp parallel do private(i,j,k)
      do j=jfirst,jlast
        do k=2,km
          do i=1,im
            pe(i,k,j) = peln(i,k,j)
          enddo
        enddo
      enddo

      call pkez(nx, im, km, jfirst, jlast,          &
                pe, pk, akap, ks, peln, pkz, .true.)

!$taf store peln,pkz,pe,delp = te_map_tape, key=te_map_tape_rec+1, kind=store_kind

! ((((((((((((((((( compute globally integrated TE >>>>>>>>>>>>>>>>

      if( consv ) then
!$omp parallel do private(i,j,k)
      do k=1,km
        do j=jfirst,jlast
           do i=1,im
              dz(i,j,k) = te(i,j,k) * delp(i,j,k)
           enddo
        enddo
      enddo

!$omp parallel do private(i, j, k, bte, xsum)

! FastOpt specify that dz is not broadcasts
!$taf omp broadcast()

      do 4000 j=jfirst,jlast
! Perform vertical integration
        if ( j == 1 ) then
! SP
          tte(1) = 0.

          do k=1,km
            tte(1) = tte(1) + dz(1,1,k)
          enddo
          tte(1)  = te_acap * tte(1)

        elseif ( j == jm) then
! NP
          tte(jm) = 0.

          do k=1,km
            tte(jm) = tte(jm) + dz(1,jm,k)
          enddo
          tte(jm) = te_acap * tte(jm)

        else
! Interior
          do i=1, im
            bte(i) = 0.
          enddo

          do k=1,km
            do i=1,im
              bte(i) = bte(i) + dz(i,j,k)
            enddo
          enddo

          xsum = 0.
          do i=1,im
            xsum = xsum + bte(i)
          enddo
          tte(j) = xsum*te_cosp(j)

        endif
4000  continue

      call par_vecsum(jm, jfirst, jlast, tte, te1)

!$omp parallel do private(i, j, xsum, ysum)
! FastOpt specify that peln,ps are not broadcast
!$taf omp broadcast()
      do j=jfirst,jlast

       if( j == 1 ) then
           tte(1) = te_acap*cp * (ps(1,1) - 2.*te_ptop -               &
                    akap*te_ptop*(peln(1,km+1,1) - peln(1,1,1) ) )
       elseif( j == jm ) then
           tte(jm)= te_acap*cp * (ps(1,jm) -                        &
                    akap*te_ptop*(peln(1,km+1,jm) - peln(1,1,jm) ) )
       else
          xsum = 0.
          ysum = 0.
        do i=1,im
          xsum = xsum + ps(i,j)
          ysum = ysum + peln(i,km+1,j)
        enddo
        tte(j) = cp*te_cosp(j)*(xsum - te_ptop*im -                   &
                             akap*te_ptop*(ysum - peln(1,1,j)*im) )
       endif
      enddo

      call par_vecsum(jm, jfirst, jlast, tte, tsum)
        
!$taf store tsum = te_map_tape, key=te_map_tape_rec+1, kind=store_kind

      dtmp = (te0 - te1) / tsum
      if( diag ) write(6,*) 'te=',te0,                         &
                            ' Energy deficit in T = ', dtmp
      else
         dtmp = 0.       ! FastOpt define dtmp to be store
      endif        ! end consv check

#if defined( SPMD )
      call mp_barrier
      call mp_recv_n(im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u)
#endif

!$taf store u,dtmp = te_map_tape, key=te_map_tape_rec+1, kind=store_kind

!$omp parallel do private(i, j, k, u2, v2, te_sp, te_np)
      do 8000 k=1,km
! Compute KE
        do j=js2g0,min(jlast+1,jm)
          do i=1,im
             u2(i,j) = u(i,j,k)**2
          enddo
        enddo

        do j=js2g0,jn2g0
          do i=1,im
             v2(i,j) = v(i,j,k)**2
          enddo
        enddo

        do j=js2g0,jn2g0
          do i=1,im-1
            te(i,j,k) = te(i,j,k) - 0.25 * ( u2(i,j) + u2(i,j+1)  &
                                            +v2(i,j) + v2(i+1,j) )
          enddo
           te(im,j,k) = te(im,j,k) - 0.25*( u2(im,j) + u2(im,j+1) &
                                           +v2(im,j) + v2(1,j) )
        enddo

! poles
        if ( jfirst == 1 ) then
          te_sp = 0.
          do i=1,im
            te_sp = te_sp + u2(i,2) + v2(i,2)
          enddo
            te_sp = te(1,1,k) - 0.5*te_sp/float(im)

          do i=1,im
             te(i,1,k) = te_sp
          enddo
        endif

        if ( jlast == jm ) then
            te_np = 0.
          do i=1,im
            te_np = te_np + u2(i,jm) + v2(i,jm-1)
          enddo

            te_np = te(1,jm,k) - 0.5*te_np/float(im)
          do i=1,im
            te(i,jm,k) = te_np
          enddo
        endif
8000  continue

!$taf store te = te_map_tape, key=te_map_tape_rec+1, kind=store_kind

! Recover (virtual) temperature
!$omp parallel do private(ixj, i1, i2, i, j, k, rg, gz, tvm, dlnp)
! FastOpt specify that delp,pe,peln,pkz,te are not broadcast
!$taf omp broadcast()
!     do 9000 j=jfirst,jlast
      do 9000 ixj=1,jp

         j  = jfirst + (ixj-1) / nx
         i1 = 1 + it * mod(ixj-1, nx)
         i2 = i1 + it - 1

!$taf init te_map_gz = static, km
         rg = akap * cp
         do i=i1,i2
            gz(i) = hs(i,j)      
         enddo

        do k=km,1,-1
!$taf store gz = te_map_gz, kind=store_kind
          do i=i1,i2
            dlnp = rg*(peln(i,k+1,j) - peln(i,k,j))
            tvm  = delp(i,j,k)*(te(i,j,k) - gz(i)) /      &
                  ( cp*delp(i,j,k) - pe(i,k,j)*dlnp )
! Update phis
            gz(i) = gz(i) + dlnp*tvm
            pt(i,j,k) = tvm         ! pt is now (virtual) temperature
          enddo

          if( consv ) then
              do i=i1,i2
                 pt(i,j,k) = pt(i,j,k) + dtmp
              enddo
          endif

          if( .not. convt ) then
              do i=i1,i2
                 pt(i,j,k) = pt(i,j,k) / pkz(i,j,k)
              enddo
           endif
        enddo           ! end k-loop
9000  continue

      return
!EOC
  end subroutine te_map


  subroutine map1_ppm( km,   pe1,   q1,                      &
                       kn,   pe2,   q2, ng_s, ng_n,          &
                       im, i1, i2, j, jfirst, jlast, iv, kord)

  use precision

! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

      implicit none

      real(r8)       r3, r23
      parameter (r3 = 1./3., r23 = 2./3.)

! Input:
      integer i1, i2
      integer im, jfirst, jlast, iv, kord
      integer km                             ! Original vertical dimension
      integer kn                             ! Target vertical dimension
      integer ng_s, ng_n

      real(r8) pe1(im,km+1)
      real(r8) pe2(im,kn+1)
      real(r8)  q1(im,jfirst-ng_s:jlast+ng_n,km)
! Output
      real(r8)  q2(im,jfirst-ng_s:jlast+ng_n,kn)

! Local
      real(r8)   dp1(i1:i2,km)
      real(r8)  q4(4,i1:i2,km)

      integer i, j, k, l, ll, k0
      real(r8)    pl, pr, qsum, delp, esl
      integer ikind                       ! FastOpt
      integer lf                          ! FastOpt

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

! Mapping
!$taf loop = parallel
      do 1000 i=i1,i2
         k0 = 1
!$taf loop = parallel
      do 555 k=1,kn

      ikind = 0
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
         if(pe2(i,k) .ge. pe1(i,l) .and. pe2(i,k) .le. pe1(i,l+1)) then
! entire new grid is within the original grid
            if(pe2(i,k+1) .le. pe1(i,l+1)) then
               k0 = l
               ikind = 1
               exit
            else
! Fractional area...
               ikind = 2
               exit
            endif
         endif
100   continue

      if (ikind .eq. 0) then
         stop 'error'
      else if (ikind .eq. 1) then
         pl = (pe2(i,k  )-pe1(i,l)) / dp1(i,l)
         pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
         q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l)) &
                      *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
      else
         pl = (pe2(i,k  )-pe1(i,l)) / dp1(i,l)
         qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+     &
                 q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*             &
                  (r3*(1.+pl*(1.+pl))))
         lf = 0
         do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
            if(pe2(i,k+1) .gt. pe1(i,ll+1) ) then
! Whole layer..
            else
               lf = ll
               k0 = ll
               exit
            endif
         enddo

         if (lf .eq. 0) then
            do ll=l+1,km
               qsum = qsum + dp1(i,ll)*q4(1,i,ll)
            enddo
         else
            do ll=l+1,lf-1
               qsum = qsum + dp1(i,ll)*q4(1,i,ll)
            enddo
            ll = lf
            delp = pe2(i,k+1)-pe1(i,ll)
            esl = delp / dp1(i,ll)
            qsum = qsum + delp*(q4(2,i,ll)+0.5*esl*               &
                    (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
         endif

         q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
      end if

555   continue
1000  continue

      return
 end subroutine map1_ppm

 subroutine mapn_ppm(km,   pe1,   q1, nq,                  &
                     kn,   pe2,   q2, ng_s, ng_n,          &
                     im, i1, i2, j, jfirst, jlast, iv, kord)

 use precision

! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

      implicit none

      real(r8)       r3, r23
      parameter (r3 = 1./3., r23 = 2./3.)

! Input:
      integer i1, i2
      integer im, jfirst, jlast, iv, kord
      integer km                             ! Original vertical dimension
      integer kn                             ! Target vertical dimension
      integer nq
      integer ng_s, ng_n

      real(r8) pe1(im,km+1)
      real(r8) pe2(im,kn+1)
      real(r8)  q1(im,jfirst-ng_s:jlast+ng_n,km, nq)
! Output
      real(r8)  q2(im,jfirst-ng_s:jlast+ng_n,kn, nq)

! Local
      real(r8)   dp1(i1:i2,km)
      real(r8)  q4(4,i1:i2,km)

      integer i, j, k, l, ll, k0, iq
      real(r8)    pl, pr, qsum, delp, esl
      integer ikind                       ! FastOpt
      integer lf                          ! FastOpt

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         enddo
      enddo

      do 2000 iq=1,nq


      do k=1,km
         do i=i1,i2
            q4(1,i,k) = q1(i,j,k,iq)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

! Mapping
!$taf loop = parallel
      do 1000 i=i1,i2
         k0 = 1

!$taf loop = parallel
      do 555 k=1,kn

      ikind = 0
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) .ge. pe1(i,l) .and. pe2(i,k) .le. pe1(i,l+1)) then
         if(pe2(i,k+1) .le. pe1(i,l+1)) then
! entire new grid is within the original grid
            k0 = l
            ikind = 1
            exit
         else
! Fractional area...
            ikind = 2
            exit
         endif
      endif
100   continue

      if (ikind .eq. 0) then
         stop 'error'
      else if (ikind .eq. 1) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
         q2(i,j,k,iq) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)


      else
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
! Fractional area...
         qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+     &
                 q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*             &
                  (r3*(1.+pl*(1.+pl))))

         lf = 0
         do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
            if(pe2(i,k+1) > pe1(i,ll+1) ) then
! Whole layer..
            else
               lf = ll
               k0 = ll
               exit
            endif
         enddo

         if (lf .eq. 0) then
            do ll=l+1,km
               qsum = qsum + dp1(i,ll)*q4(1,i,ll)
            enddo
         else
            do ll=l+1,lf-1
               qsum = qsum + dp1(i,ll)*q4(1,i,ll)
            enddo
            ll = lf
            delp = pe2(i,k+1)-pe1(i,ll)
            esl  = delp / dp1(i,ll)
            qsum = qsum + delp*(q4(2,i,ll)+0.5*esl*                 &
                    (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
         endif

         q2(i,j,k,iq) = qsum / ( pe2(i,k+1) - pe2(i,k) )
      endif

555   continue
1000  continue
2000  continue

      return
 end subroutine mapn_ppm


 subroutine ppm2m(a4, delp, km, i1, i2, iv, kord)

! iv =-1: winds
! iv = 0: positive definite scalars stored in the q-array
! iv = 1: others

 use precision

      implicit none
! Input
      integer km, lmt, iv
      integer i1, i2
      integer kord
      real(r8)    delp(i1:i2,km)
      real(r8)    a4(4,i1:i2,km)

! local arrays.
      real(r8)  dc(i1:i2,km)
      real(r8)  h2(i1:i2,km)
      real(r8) delq(i1:i2,km)

      real(r8) a1, a2, c1, c2, c3, d1, d2
      real(r8) qmax, qmin, cmax, cmin
      real(r8) qm, dq, tmp

! Local scalars:
      integer i, k, km1
      real(r8) qmp, pmp
      real(r8) lac
      integer it
      real(r8) df2(i1:i2,km)
      real(r8) d4(i1:i2,km)
      real(r8) fac

      integer, parameter :: store_kind = 8      ! FastOpt precision of tape

      km1 = km - 1
       it = i2 - i1 + 1

!$taf init ppm2m = static, 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

!$taf store a4(1,i1:i2,:) = ppm2m, kind=store_kind

      do k=2,km1
         do i=i1,i2
            c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
             dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
         enddo
      enddo

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

      do k=3,km1
      do i=i1,i2
      c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
      a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
      a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
      a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

!$taf store a4(1:2,i1:i2,:) = ppm2m, kind=store_kind

      call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top

!$taf store a4(2,i1:i2,1:3), a4(1,i1:i2,1:2) = ppm2m, kind=store_kind

      do i=i1,i2
      d1 = delp(i,1)
      d2 = delp(i,2)
      qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
      dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
      c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
      a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
      a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
      dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
      cmax = max(a4(1,i,1), a4(1,i,2))
      cmin = min(a4(1,i,1), a4(1,i,2))
      a4(2,i,2) = max(cmin,a4(2,i,2))
      a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo

!$taf store a4(1,i1:i2,1), a4(2,i1:i2,1:2) = ppm2m, kind=store_kind

      if( iv == 0 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      elseif ( iv == -1 ) then
! Winds:
        if( km > 32 ) then
          do i=i1,i2
! More dampping: top layer as the sponge
             a4(2,i,1) = a4(1,i,1)
             a4(3,i,1) = a4(1,i,1)
          enddo
        else
          do i=i1,i2
             if( a4(1,i,1)*a4(2,i,1) <=  0. ) then
                 a4(2,i,1) = 0.
             else
                 a4(2,i,1) = sign(min(abs(a4(1,i,1)),    &
                                      abs(a4(2,i,1))),   &
                                          a4(1,i,1)  )
            endif
          enddo
        endif
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface

!$taf store a4(1,i1:i2,km1:km), a4(2,i1:i2,km1:km), a4(3,i1:i2,km) &
!$taf&      = ppm2m, kind=store_kind

      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

! Enforce constraint at the surface

!$taf store a4(1,i1:i2,km), a4(3,i1:i2,km) = ppm2m, kind=store_kind

      if ( iv == 0 ) then
! Positive definite scalars:
           do i=i1,i2
              a4(3,i,km) = max(0._r8, a4(3,i,km))
           enddo
      elseif ( iv == -1 ) then
! Winds:
           do i=i1,i2
              if( a4(1,i,km)*a4(3,i,km) <=  0. ) then
                  a4(3,i,km) = 0.
              else
                  a4(3,i,km) = sign( min(abs(a4(1,i,km)),   &
                                         abs(a4(3,i,km))),  &
                                             a4(1,i,km)  )
              endif
           enddo
      endif

      do k=1,km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo
 
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
 
! Top 2 and bottom 2 layers always use monotonic mapping
!$taf store a4(:,i1:i2,1:2) = ppm2m, kind=store_kind

!$taf loop = parallel
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

!$taf store a4(:,i1:i2,3:km-2) = ppm2m, kind=store_kind

      if(kord .ge. 7) then
!****6***0*********0*********0*********0*********0*********0**********72
! Huynh's 2nd constraint
!****6***0*********0*********0*********0*********0*********0**********72
      do k=2, km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord .eq. 7 ) then
         fac = 1.5           ! original quasi-monotone
      else
         fac = 0.125         ! full monotone
      endif

!$taf loop = parallel
      do k=3, km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
	enddo                                          ! FastOpt

        do i=i1,i2                                     ! FastOpt
         pmp   = 2.*dc(i,k)                            ! FastOpt
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to prevent negatives when kord=7
         if (iv .eq. 0 .and. kord .eq. 7) then
             call kmppm(dc(i1,k), a4(1,i1,k), it, 2)
         endif
      enddo

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv .eq. 0) lmt = min(2, lmt)

!$taf loop = parallel
      do k=3, km-2
      if( kord .ne. 4) then
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
      endif
         call kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

!$taf store a4(:,i1:i2,km1:km) = ppm2m, kind=store_kind

!$taf loop = parallel
      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      return
 end subroutine ppm2m


 subroutine kmppm(dm, a4, im, lmt)

 use precision

 implicit none

      real(r8)       r12
      parameter (r12 = 1./12.)

      integer im, lmt
      real(r8)    a4(4,im)  ! FastOpt : avoid assumed-size dummy argument
      real(r8)    dm(im)    ! FastOpt : avoid assumed-size dummy argument
      real(r8) qmp

! PPM notations:
! AA <-- a4(1,i)
! AL <-- a4(2,i)
! AR <-- a4(3,i)
! A6 <-- a4(4,i)

      integer i
      real(r8) da1, da2, a6da
      real(r8) fmin

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt .eq. 3 ) return

      if(lmt .eq. 0) then
! Standard PPM constraint
      do i=1,im
      if(dm(i) .eq. 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da .lt. -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da .gt. da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt .eq. 1) then

! Improved full monotonicity constraint (Lin)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, im
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt .eq. 2) then

! Positive definite constraint
      do i=1,im
      if( abs(a4(3,i)-a4(2,i)) .lt. -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin .lt. 0. ) then
         if(a4(1,i).lt.a4(3,i) .and. a4(1,i).lt.a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) .gt. a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

      return
 end subroutine kmppm

 subroutine steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)

! !USES:
   use precision
   implicit none

! !INPUT PARAMETERS:
      integer km                        ! Total levels
      integer i1                        ! Starting longitude
      integer i2                        ! Finishing longitude
      real(r8)  dp(i1:i2,km)            ! grid size
      real(r8)  dq(i1:i2,km)            ! backward diff of q
      real(r8)  d4(i1:i2,km)            ! backward sum:  dp(k)+ dp(k-1) 
      real(r8) df2(i1:i2,km)            ! first guess mismatch
      real(r8)  dm(i1:i2,km)            ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real(r8)  a4(4,i1:i2,km)          ! first guess/steepened

!
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, k
      real(r8) alfa(i1:i2,km)
      real(r8)    f(i1:i2,km)
      real(r8)  rat(i1:i2,km)
      real(r8)  dg2

! Compute ratio of dq/dp
      do k=2,km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=2,km-1
         do i=i1,i2
            f(i,k) = (rat(i,k+1) - rat(i,k))                             &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=3,km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1).lt.0. .and. df2(i,k).ne.0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0._r8, min(0.5_r8, -0.1875*dg2/df2(i,k))) 
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=4,km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

      return
!EOC
 end subroutine steepz

end module mapz_module
