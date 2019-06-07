#include <params.h>
      subroutine moistdrv( 
     &           nstep,   bdt,     coslon,  sinlon, 
     &           ptop,    phis,    pblht,   ps,      ptold,
     &           tpert,   qpert,   delp,    u,       v,
     &           pt,      q    ,
     &           jfirst,  jlast,   ng_d,    ng_s)
!-----------------------------------------------------------------------
!
! Purpose:
! Moist physics driver to reduce the excessive precipitation right after fvPSAS
!
! Method:
! Call model moist physics to adjust pt and q
!
! Author: Jiundar Chern, April 2001
!
!-----------------------------------------------------------------------

      implicit none


!----------------- common blocks ---------------------------------------
#include <pmgrid.h>  ! use only: plond, plat, plev, plevp, pcnst, pnats, beglat, endlat
#include <comcon.h>  ! use only: cappa, zvir, gravit, rair, cpair, latvap, rga, tmelt
#include <comsrf.h>  ! use only: ts, wsx, wsy, shf, cflx
!---------------------------------------------------------------------
!
! Input Arguments
! 
      integer, intent(in) :: nstep                   ! Current time step
      integer, intent(in) :: jfirst, jlast
      integer, intent(in) ::  ng_d, ng_s
      real, intent(in) :: bdt                        ! 2 delta t (model time increment)
      real, intent(in) :: coslon                     ! cosine of latitude
      real, intent(in) :: sinlon                     ! sine of latitude
      real, intent(in) :: ptop                       ! pressure at model top
      real, intent(in) :: phis(plond,jfirst:jlast)   ! Surface geopotential (m2/s2)
      real, intent(in) :: ps(plond,jfirst:jlast)     ! Surface pressure (pa)
      real, intent(in) :: ptold(plond,jfirst-ng_d:jlast+ng_d,plev)   ! scaled virtual potential 
                                                                     ! temperature (K) before fvPSAS
      real, intent(in) :: delp(plond,jfirst:jlast,plev)    ! Delta pressure (Pa) after fvPSAS
      real, intent(in) :: u(plond,jfirst-ng_d:jlast+ng_s,plev)       ! u wind (m/s) after fvPSAS
      real, intent(in) :: v(plond,jfirst-ng_s:jlast+ng_d,plev)       ! v wind (m/s) after fvPSAS
      real, intent(in) :: pblht(plond,plat)        ! Planetary boundary layer height
      real, intent(in) :: tpert(plond,plat)        ! Thermal temperature excess (K)
      real, intent(in) :: qpert(plond,pcnst,plat)  ! Thermal humidity & constituent excess
!
! Input/Output Arguments
!
      real, intent(inout) :: pt(plond,jfirst-ng_d:jlast+ng_d,plev) ! scaled virtual potential 
                                                                   ! temperature (K)
      real, intent(inout) :: q(plond,jfirst-ng_d:jlast+ng_d,plev,pcnst+pnats)   ! Specific humidity (kg/kg)

! Local Variables

      real ua(plond,jfirst:jlast,plev)   ! u wind at A grid
      real va(plond,jfirst:jlast,plev)   ! v wind at A grid
      real wtadj(plond)                  ! weighting for moist physics adjustment
      real wkadj(plond)                  ! working space for moist physics adjustment

      real pm(plond,plev)                ! Pressure at middle points (Pa)
      real pe2d(plond,plevp)             ! Pressure at layer edge (Pa)      
      real pk2d(plond,plevp)             ! pe2d**cappa
      real peln2d(plond,plevp)           ! log( pe2d )
      real pmln(plond,plev)              ! 0.5*(peln2d(i,k) + peln2d(i,k+1))
      real pkz(plond,plev)               ! layer-mean pk for converting temperature
                                         !  to potential temperature
      real psomc(plond,plev)
      real tvm(plond,plev)               ! Virtual temperature (K)
      real told(plond,plev)              ! Temperature (K) before adjustment
      real qold(plond,plev,pcnst+pnats)  ! Specific humidity (kg/kg) before adjustment
      real pdel(plond,plev)              ! Delta pressure (Pa)
      real tnew(plond,plev)              ! adjusted temperature (K)
      real qnew(plond,plev,pcnst+pnats)  ! adjusted specific humidity (kg/kg)
      real pblht1(plond)                 ! Planetary boundary layer height
      real tpert1(plond)                 ! Thermal temperature excess (K)
      real qpert1(plond,pcnst)           ! Thermal humidity & constituent excess
      real um1(plond,plev)               ! u-wind before vertical diffusion
      real vm1(plond,plev)               ! v-wind before vertical diffusion
      real lnptop                        ! log( ptop )
      real ptopk                         ! ptop**cappa
      real tmp1                          ! dummy variable

      logical filz
      integer i, j, k, ic

      lnptop      = log( ptop)
      ptopk       = ptop ** cappa
      filz        = .true.  
!
! Convert u and v from D-Grid to A-Grid
!
      call d2a3d(u, v, ua, va, plond, plat, plev, jfirst, jlast,
     &           ng_d, ng_s, coslon, sinlon)

!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j,k,ic)
!$omp& private(pe2d, peln2d, pk2d, pkz, pm, pmln, pdel, psomc)
!$omp& private(tpert1, qpert1, pblht1, qold, tvm, told)
!$omp& private(um1, vm1, qnew, tnew, wtadj, wkadj, tmp1)

      do j = jfirst, jlast      ! Parallel j-loop
!
! Compute variables needed to call tphysaa
!
        do i = 1, plond
          tpert1(i)     = tpert(i,j)
          pblht1(i)     = pblht(i,j)
        end do
   
        do ic = 1, pcnst
          do i = 1, plond
            qpert1(i,ic) = qpert(i,ic,j)
          end do
        end do

        do ic = 1, pcnst+pnats
          do k = 1, plev
            do i = 1, plond
              qold(i,k,ic) = q(i,j,k,ic)
            end do
          end do
        end do

        do i = 1, plond
          pe2d(i,1)     = ptop
          peln2d(i,1)   = lnptop
          pk2d(i,1)     = ptopk
        enddo

        do k = 2, plevp
          do i = 1, plond
            pe2d(i,k)   = pe2d(i,k-1) + delp(i,j,k-1)
            peln2d(i,k) = log( pe2d(i,k) )
            pk2d(i,k)   = pe2d(i,k) ** cappa
          end do
        end do
    
        do k = 1, plev
          do i = 1, plond
            pm(i,k)     = 0.5 * ( pe2d(i,k) + pe2d(i,k+1) )
            pmln(i,k)   = 0.5 * ( peln2d(i,k) + peln2d(i,k+1) )
            pdel(i,k)   = delp(i,j,k)
            pkz(i,k)    = ( pk2d(i,k) - pk2d(i,k+1) ) / cappa
     &                  / ( peln2d(i,k) - peln2d(i,k+1) )
            psomc(i,k)  = pk2d(i,plevp) / pkz(i,k)
            tvm(i,k)    = pt(i,j,k) * pkz(i,k)
            told(i,k)   = tvm(i,k) / ( 1.0 + zvir * qold(i,k,1) )

            um1(i,k)    = ua(i,j,k)
            vm1(i,k)    = va(i,j,k)
          end do
        end do

        if(filz) call fillz(plond, plev, pcnst+pnats, qold, pdel)


! Call Moist physics after analysis

!****6***0*********0*********0*********0*********0*********0**********72
        call tphysaa(j     ,nstep    ,bdt        ,
     &           gravit    ,rair     ,zvir       ,
     &           cpair     ,latvap   ,rga        ,tmelt    ,
     &           phis(1,j) ,ts(1,j)  ,pblht1     ,tpert1   ,qpert1  ,
     &           pe2d      ,peln2d   ,pdel       ,pm       ,pmln    ,
     &           told      ,qold     ,tvm        , 
     &           um1       ,vm1      ,psomc      ,
     &           wsx(1,j)  ,wsy(1,j) ,shf(1,j)   ,cflx(1,1,j)       ,
     &           tnew      ,qnew                       )
!****6***0*********0*********0*********0*********0*********0**********72
!
! compute adjustment weighting
!
        do i = 1, plond
          wtadj(i) = 0.0
          wkadj(i) = 0.0
        end do

        do k = plev, 1, -1
          do i = 1, plond
            wkadj(i) = wkadj(i)
     &               + abs( pt(i,j,k) - ptold(i,j,k) )
     &               / ptold(i,j,k) * delp(i,j,k)
          enddo
        enddo

        do i = 1, plond
          wkadj(i) = wkadj(i) / ( ps(i,j) - ptop )
          wtadj(i) = exp( -2.5e3 * wkadj(i) )
        enddo
!
! only update pt and specific humidity q(i,j,k,1)
!
! no data, no harm
!
        do k = 1, plev
          do i = 1, plond
            tmp1        = tnew(i,k) * (1.0 + zvir * qnew(i,k,1) ) 
     &                  / pkz(i,k)
            pt(i,j,k)   = wtadj(i) * pt(i,j,k)
     &                  + ( 1.0 - wtadj(i) ) * tmp1
            q(i,j,k,1)  = wtadj(i) * q(i,j,k,1)
     &                  + ( 1.0 - wtadj(i) ) * qnew(i,k,1)
          enddo
        end do

      enddo               ! End Parallel j-loop
      return
      end
