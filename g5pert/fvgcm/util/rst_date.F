      program rst_date
c
c ... Print out the time stamp embedded in dynamics
c     restart file (d_rst).
c
      implicit none

#include "fvgcm.h"

      character*254 flnm
      character*8  buf

      integer argnum
      integer iargc
      integer ierr
      integer len_trim

      integer ir, n, i, n2, n3
ccc   integer nstep, nymd, nhms
      integer*4 nstep, nymd, nhms

      logical verbose, stats
      logical exist

      real, allocatable :: ps(:), delp(:), u(:), v(:), pt(:), q(:,:) 
      integer, parameter :: im = FVGCM_LON
      integer, parameter :: jm = FVGCM_LAT
      integer, parameter :: km = FVGCM_LEV
      integer, parameter :: nq = FVGCM_TRACER

c
c ... Nstep will be printed out in verbose mode
c
      stats = .false.
      verbose = .false.
      n = iargc()
      if (n .eq. 0 .or. n .gt. 2) then
        print *
        print *, 'Usage: rst_date [-v|-s] filename'
        print *
      else
        verbose = .false.
        if (n .eq. 2) then
          argnum = 1
          call Getarg (argnum, buf) 
          if (buf(:2) .eq. '-v') then
            verbose = .true.
          endif
          if (buf(:2) .eq. '-s') then ! Print stats
            stats = .true.
          endif
          argnum = 2
          call Getarg (argnum, flnm) 
        else
          argnum = 1
          call Getarg (argnum, flnm) 
        endif

        inquire (file=flnm, exist=exist)

        if (exist) then
          ir = 55

#if defined( CRAY_T3E )
          call asnunit( ir, '-R', ierr )
          call asnunit( ir, '-F f77', ierr )
#elif defined( CRAY )
          call asnunit( ir, '-R', ierr )
          call asnunit( ir, '-F f77 -N ieee_dp', ierr )
#endif

          open (ir, file=flnm, form='unformatted', status='old')
          read (ir) nstep, nymd, nhms
          if (verbose .or. stats) then
            print '(a12,i8.8)', ' YYYYMMDD = ', nymd
            print '(a12,i6.6)', ' HHMMSS   =   ', nhms
            print *, 'NSTEP    =', nstep
          else
            print '(i8.8,1x,i6.6)', nymd, nhms
          endif
          if ( stats ) then
              n2 = im*jm
              n3 = im*jm*km
              allocate ( ps(n2), delp(n3), u(n3), v(n3), pt(n3), q(n3,nq), 
     .                   stat=ierr)
              if ( ierr .ne. 0 ) then
                 print *, 'rst_date: cannot allocate memory'
                 call exit(1)
              end if
              read(ir, iostat=ierr) ps, delp, u, v, pt
              if ( ierr .ne. 0 ) then
                 print *, 'rst_date: cannot read ps, u, v, pt'
                 call exit(1)
              end if
              read(ir, iostat=ierr) q
              if ( ierr .ne. 0 ) then
                 print *, 'rst_date: cannot read q'
                 call exit(1)
              end if
              print *, 'Min/Max Values:'
              print *, '   im = ', im
              print *, '   jm = ', jm
              print *, '   km = ', km
              print *, '   nq = ', nq
              print *, '   PS ', ' max = ', maxval(ps), ' min = ', minval(ps)
              print *, ' delP ', ' max = ', maxval(delp), 
     .                   ' min = ', minval(delp)
              print *, '    U ', ' max = ', maxval(u), ' min = ', minval(u)
              print *, '    V ', ' max = ', maxval(v), ' min = ', minval(v)
              print *, '   PT ', ' max = ', maxval(pt), ' min = ', minval(pt)
              do i = 1, nq
              print *, '    Q ', ' max = ', maxval(q(1:n3,i)), 
     .                 ' min = ', minval(q(1:n3,i))
              end do
              deallocate ( ps, delp, u, v, pt, q )
          end if
          close (ir)

#if defined( CRAY ) || defined( CRAY_T3E )
          call asnunit( ir, '-R', ierr )
#endif

        else
          if (verbose) then
            print '(/a,a,a/)', '  Error --> File ''', 
     .            flnm(:len_trim(flnm)), ''' does not exist.'
          endif
        endif
      endif

      end


      integer function len_trim (string)

************************************************************************
*                                                                      * 
*     len_trim.f - return the length of string without the trailing    * 
*                  blanks and tabs                                     * 
*                                                                      * 
*     Last Modified:  Thu Mar 13 17:47:34 EST 1997                     * 
*                                                                      * 
************************************************************************

      implicit         none

      integer          tab, blank
      parameter        (tab   =  9)
      parameter        (blank = 32)

      character*(*)    string
      integer          n, code

      len_trim = 0
      n = len(string)
      if (n .eq. 0) return

      do while (n .gt. 0)
        code = ichar(string(n:n))
        if (code .eq. blank  .or. code .eq. tab) then
          len_trim = n - 1
	  n = n - 1
        else
	  len_trim = n
          n = 0
        end if
      end do

      return 
      end


#if defined( CRAY_T3E )

      subroutine Getarg (argnum, arg)

      implicit none

      integer       :: argnum
      character*(*) :: arg

      integer :: ierror
      integer :: ilen

      call Pxfgetarg (argnum, arg, ilen, ierror)

      if (ierror /= 0) then
        Write (6,*) 'Problem in Getarg/CRAY_T3E.'
        Stop
      end if

      return

      end

#endif

