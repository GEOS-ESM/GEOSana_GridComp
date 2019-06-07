      subroutine zout(iou,iavg,jm,km,jfirst,jlast,psx,ux,vx,tx, qx,   &
                      wx, pvx,usus,usvs,tsts,vsts,vspv)

      use precision
      implicit none

      integer jm, km, jfirst, jlast
      integer iou,iavg
      real(r8)   psx(jfirst:jlast)
      real(r8)    ux(jfirst:jlast,km)
      real(r8)    vx(jfirst:jlast,km)
      real(r8)    tx(jfirst:jlast,km)
      real(r8)    qx(jfirst:jlast,km)
      real(r8)    wx(jfirst:jlast,km)
      real(r8)   pvx(jfirst:jlast,km)
      real(r8)  usus(jfirst:jlast,km)
      real(r8)  usvs(jfirst:jlast,km)
      real(r8)  vsts(jfirst:jlast,km)
      real(r8)  vspv(jfirst:jlast,km)
      real(r8)  tsts(jfirst:jlast,km)
      real(r8) rt
      real(r8) rw
      integer j,k

      rt = 1./float(iavg)
      rw = 864./float(iavg)

        do j=jfirst,jlast
            psx(j) = psx(j) * rt
        enddo

!$omp  parallel do           &
!$omp  default(shared)       &
!$omp  private(j,k)

      do k=1,km
        do j=jfirst,jlast
            ux(j,k) = ux(j,k) * rt
            vx(j,k) = vx(j,k) * rt
            tx(j,k) = tx(j,k) * rt
            qx(j,k) = qx(j,k) * rt
            wx(j,k) = wx(j,k) * rw
           pvx(j,k) = pvx(j,k) * rt
          usus(j,k) = usus(j,k) * rt
          usvs(j,k) = usvs(j,k) * rt
          tsts(j,k) = tsts(j,k) * rt
          vsts(j,k) = vsts(j,k) * rt
          vspv(j,k) = vspv(j,k) * rt
        enddo
      enddo

      call wrt3dr(iou,1,jm, 1,psx, jfirst, jlast)
      call wrt3dr(iou,1,jm,km,ux,  jfirst, jlast)
      call wrt3dr(iou,1,jm,km,vx,  jfirst, jlast)
      call wrt3dr(iou,1,jm,km,wx,  jfirst, jlast)
      call wrt3dr(iou,1,jm,km,tx,  jfirst, jlast)
      call wrt3dr(iou,1,jm,km,qx,  jfirst, jlast)
      call wrt3dr(iou,1,jm,km,pvx, jfirst, jlast)

      call wrt3dr(iou,1,jm,km,usus, jfirst, jlast)
      call wrt3dr(iou,1,jm,km,usvs, jfirst, jlast)
      call wrt3dr(iou,1,jm,km,tsts, jfirst, jlast)
      call wrt3dr(iou,1,jm,km,vsts, jfirst, jlast)
      call wrt3dr(iou,1,jm,km,vspv, jfirst, jlast)

      return
      end
