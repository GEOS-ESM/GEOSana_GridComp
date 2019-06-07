      subroutine wrt3dr(iout, im, jm, km, a3, jfirst, jlast)

      use precision

#if defined( SPMD )
      use mod_comm, only : gid, La_2_Ga
#define CPP_ARRAY  global
#else
#define CPP_ARRAY  a3
#endif

      implicit none
      integer iout, im, jm, km
      integer i, j, k
      integer jfirst, jlast

      real(r8) a3(im,jfirst:jlast,km)     ! local slice
      real(r4) a2(im,jm)

#if defined( SPMD  )
      real(r8) global(im,jm,km)

      call La_2_Ga(a3, global, im, jm, km, jfirst, jlast, 1, km)
      if ( gid .eq. 0 ) then
#endif

        call zflip(CPP_ARRAY, im*jm, km, 1)
        do 50 k=1,km
          do j=1,jm
            do i=1,im
              if(abs(CPP_ARRAY(i,j,k)) .lt. 1.e-25) then
                a2(i,j) = 0.
              else
                a2(i,j) = CPP_ARRAY(i,j,k)
              endif
            enddo
          enddo
          write(iout) a2
50      continue

#if defined( SPMD )
      endif
#endif

      return
      end
