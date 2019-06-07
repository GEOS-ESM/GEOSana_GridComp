      real function vmax(a, qmin, im, jm, jfirst, jlast)

#if defined (SPMD )
      use mod_comm
#endif

      use precision
      implicit none

      integer im, jm, jfirst, jlast
!!!      real a(im, jfirst:jlast)
      real a(im,jm)

      integer i, j
      real qmax, qmin

      qmax = a(1,jfirst)
      qmin = qmax

      do j=jfirst, jlast
         do i=1,im
            qmin = min(qmin, a(i,j))
            qmax = max(qmax, a(i,j))
         enddo
      enddo

#if defined (SPMD )
      call mp_minmax(qmin, qmax)
#endif

      vmax = qmax
      return
      end
