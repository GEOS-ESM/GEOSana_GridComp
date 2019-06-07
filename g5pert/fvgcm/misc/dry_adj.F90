      subroutine dry_adj(im,  km,  rdt, pt, fu, fv,      &   
                         u,   v,   dp, j) 

      implicit none

      integer i, j, k, im, km
      real pt(im,km)
      real dp(im,km)
      real rdt
      real u(im,km)
      real v(im,km)
      real fu(im,km)
      real fv(im,km)

! local
      real ut(im,km)
      real vt(im,km)
      real dp1, dp2
      real ptm, um, vm
      integer klow
      logical mixm

      mixm = .true.

      klow = max(km/16, 1)

      do k=1,klow+3
         do i=1,im
            ut(i,k) = u(i,k)
            vt(i,k) = v(i,k)
         enddo
      enddo

! FastOpt additional statements to break data dependence
! FastOpt write tape because intermediate results are lost

      do k=1,klow

!$taf store pt(:,k:k+1) = dry_adj_tape, rec=k

!$taf loop = parallel
      do i=1,im
         if( pt(i,k) < pt(i,k+1) ) then
!           write(*,*) "unstable points:",j, k,i

!$taf store dp(i,k:k+3) = dry_adj_tape, rec=k
!$taf store u (i,k:k+3) = dry_adj_tape, rec=k
!$taf store v (i,k:k+3) = dry_adj_tape, rec=k

            dp1 = dp(i,k) + dp(i,k+1)
            ptm = (pt(i,k)*dp(i,k)+pt(i,k+1)*dp(i,k+1)) / dp1
            pt(i,k)   = ptm
            pt(i,k+1) = ptm
! winds
            if ( mixm ) then
             um = (u(i,k)*dp(i,k)+u(i,k+1)*dp(i,k+1)) / dp1
             u(i,k  ) = um
             u(i,k+1) = um
             vm = (v(i,k)*dp(i,k)+v(i,k+1)*dp(i,k+1)) / dp1
             v(i,k  ) = vm
             v(i,k+1) = vm
            endif

            if( ptm < pt(i,k+2) ) then
!               write(*,*) "2nd stage in dry_adj"
                dp2 = dp1 + dp(i,k+2)
                ptm = (ptm*dp1 + pt(i,k+2)*dp(i,k+2)) / dp2
                pt(i,k  ) = ptm
                pt(i,k+1) = ptm
                pt(i,k+2) = ptm
! Winds:
             if ( mixm ) then
                um = (u(i,k)*dp(i,k)+u(i,k+1)*dp(i,k+1)) / dp1    ! FastOpt
                um = (um*dp1 + u(i,k+2)*dp(i,k+2)) / dp2
                u(i,k  ) = um
                u(i,k+1) = um
                u(i,k+2) = um
                vm = (v(i,k)*dp(i,k)+v(i,k+1)*dp(i,k+1)) / dp1    ! FastOpt
                vm = (vm*dp1 + v(i,k+2)*dp(i,k+2)) / dp2
                v(i,k  ) = vm
                v(i,k+1) = vm
                v(i,k+2) = vm
             endif

             if( ptm < pt(i,k+3) ) then
!               write(*,*) "3rd stage in dry_adj"
                ptm = (ptm*dp2 + pt(i,k+3)*dp(i,k+3)) / (dp2 + dp(i,k+3)) 
                pt(i,k  ) = ptm
                pt(i,k+1) = ptm
                pt(i,k+2) = ptm
                pt(i,k+3) = ptm
! winds:
              if ( mixm ) then
                um = (u(i,k)*dp(i,k)+u(i,k+1)*dp(i,k+1)) / dp1    ! FastOpt
                um = (um*dp1 + u(i,k+2)*dp(i,k+2)) / dp2          ! FastOpt
                um = (um*dp2 + u(i,k+3)*dp(i,k+3)) / (dp2 + dp(i,k+3)) 
                u(i,k  ) = um
                u(i,k+1) = um
                u(i,k+2) = um
                u(i,k+3) = um
                vm = (v(i,k)*dp(i,k)+v(i,k+1)*dp(i,k+1)) / dp1    ! FastOpt
                vm = (vm*dp1 + v(i,k+2)*dp(i,k+2)) / dp2          ! FastOpt
                vm = (vm*dp2 + v(i,k+3)*dp(i,k+3)) / (dp2 + dp(i,k+3)) 
                v(i,k  ) = vm
                v(i,k+1) = vm
                v(i,k+2) = vm
                v(i,k+3) = vm
              endif
             endif
            endif
         endif
      enddo
      enddo

      if ( mixm ) then
      do k=1,klow+3
         do i=1,im
            fu(i,k) = fu(i,k) + (u(i,k) - ut(i,k)) * rdt
            fv(i,k) = fv(i,k) + (v(i,k) - vt(i,k)) * rdt
         enddo
      enddo
      endif

      return
      end
