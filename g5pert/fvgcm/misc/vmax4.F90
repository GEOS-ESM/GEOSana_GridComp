      real*4 function vmax4(a, qmin4, im, jm, jfirst, jlast)

#if defined (SPMD )
      use mod_comm
#endif

      use precision
      implicit none

      integer im, jm, jfirst, jlast
      real(r4) qmin4
      real(r4) a(im, jfirst:jlast)

      integer i, j
      real qmax, qmin

      vmax4 = a(1,jfirst)
      qmin4 = vmax4

      do j=jfirst, jlast
         do i=1,im
            qmin4 = min(qmin4, a(i,j))
            vmax4 = max(vmax4, a(i,j))
         enddo
      enddo

#if defined (SPMD )
      qmin = qmin4
      qmax = vmax4
      call mp_minmax(qmin, qmax)
      qmin4 = qmin
      vmax4 = qmax
#endif
      return
      end
