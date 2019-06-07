      integer function memget (nbytes)

      implicit none

      integer malloc
      integer nbytes
      integer nwords
      integer ier

#if (defined CRAY)
      integer ptr
      nwords = nbytes / 4
      call hpalloc(ptr,nwords,ier,0)
      memget = ptr

#elif (defined AIX)
      memget = malloc(%VAL(nbytes))

#else

      memget = malloc(nbytes)
#endif

      return
      end


      subroutine memfree (ptr)

      implicit none

      real space
      integer ier

      pointer (ptr,space)

#if (defined CRAY)
      call hpdeallc(ptr,ier,0)

#elif (defined AIX)
      call free (%VAL(ptr))

#else
      call free (ptr)
#endif

      return
      end

