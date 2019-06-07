!-----------------------------------------------------------------------
!BOP
! !ROUTINE: geopk --- Calculate geopotential to the kappa
!
! !INTERFACE:
      subroutine geopk(ptop, pe, delp, pk, wz, hs, pt, im, jm, km,  &
                       jfirst, jlast, nd, cp, akap, nx, id, dp_check)

      use precision
      implicit none

! !INPUT PARAMETERS:

      integer im, jm, km, jfirst, jlast, id
      integer nx                        ! # of pieces in longitude direction
      integer nd
      real(r8)    akap, cp, ptop
      logical dp_check
      real(r8) hs(im,jfirst:jlast)

! !INPUT/OUTPUT PARAMETERS:
      real(r8)  pt(im,jfirst-nd:jlast+nd,km)  ! only altered if dp_check
      real(r8) delp(im,jfirst:jlast,km)      ! only altered if dp_check

! !OUTPUT PARAMETERS
      real(r8) wz(im,jfirst-1:jlast+1,km+1)  ! space N*1 S*1
      real(r8) pk(im,jfirst-1:jlast+1,km+1)  ! space N*1 S*1
      real(r8) pe(im,km+1,jfirst:jlast)      ! only if id .eq. 1

! !DESCRIPTION:
!     Calculates geopotential and pressure to the kappa.  This is an expensive
!     operation and several out arrays are kept around for future use.
!
! !REVISION HISTORY:
!
!  WS  99.10.22: MPIed SJ's original SMP version
!  SJL 00.01.01: Merged C-core and D-core computation
!                SMP "decmposition" in E-W by combining i and j loops
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
      integer i, j, k
      real(r8)    p1d(im)
      integer ixj, jp, it, i1, i2

      real(r8) ptk
      real(r8) dpmin
      real(r8) dp
      parameter (dpmin = 0.1)              ! unit: pascal

      integer, parameter :: store_kind = 8   ! FastOpt precision of tape

      it = im / nx
      jp = nx * ( jlast - jfirst + 1 )

! FastOpt begin
      do k = 1,km+1
        do j = jfirst-1,jlast+1
          do i = 1,im
            pk(i,j,k) = 0.
          enddo
        enddo
      enddo
! FastOpt end

!$omp parallel do default(shared) private(i1, i2, ixj, i, j, k, p1d, ptk, dp)

!     do 2000 j=jfirst,jlast
      do 2000 ixj=1, jp

!$taf init geopk_tape_1 = static, 1
!$taf init geopk_tape   = static, km

         j  = jfirst + (ixj-1)/nx
         i1 = 1 + it * mod(ixj-1, nx)
         i2 = i1 + it - 1

           ptk  = ptop ** akap
        do i=i1,i2
           p1d(i) = ptop
           pk(i,j,1) = ptk
           wz(i,j,km+1) = hs(i,j)
        enddo

        if(id .eq. 1) then
           do i=i1,i2
              pe(i,1,j) = ptop
           enddo
        endif

        if( dp_check ) then

          do k=1, km-1

!$taf store delp(:,j,k),pt(:,j,k:k+1) = geopk_tape, kind=store_kind

             do i=i1,i2
              if(delp(i,j,k) < dpmin) then
! Remap from below and mix pt
                dp = dpmin - delp(i,j,k)
                pt(i,j,k) = (pt(i,j,k)*delp(i,j,k) + pt(i,j,k+1)*dp) / dpmin
                delp(i,j,k) = dpmin
                delp(i,j,k+1) = delp(i,j,k+1) - dp
              endif
            enddo
          enddo

!$taf store delp(:,j,k),pt(:,j,km-1:km) = geopk_tape_1, kind=store_kind

! Bottom (k=km):
          do i=i1,i2
            if(delp(i,j,km) < dpmin) then
! Remap from above and mix pt
              dp = dpmin - delp(i,j,km)
              pt(i,j,km) = (pt(i,j,km)*delp(i,j,km) + pt(i,j,km-1)*dp)/dpmin
              delp(i,j,km) = dpmin
              delp(i,j,km-1) = delp(i,j,km-1) - dp
            endif
          enddo
        endif

! Top down
        do k=2,km+1

!$taf store p1d = geopk_tape, kind=store_kind

          do i=i1,i2
            p1d(i)  = p1d(i) + delp(i,j,k-1)
            pk(i,j,k) = p1d(i) ** akap
          enddo
        if(id .eq. 1) then
           do i=i1,i2
              pe(i,k,j) = p1d(i)
           enddo
        endif
        enddo

! Bottom up
        do k=km,1,-1
          do i=i1,i2
            wz(i,j,k) = wz(i,j,k+1)+cp*pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
          enddo
        enddo
2000  continue

      return
      end
