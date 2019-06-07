      module timingModule

      use precision

!
! ... Use system etime() function for timing
!
#if  (defined  TIMING)

#if defined( SPMD )
         use mod_comm, only : gid, gsize, mp_bcst_int, mp_reduce_max
#endif

#if defined( CRAY ) || defined( CRAY_T3E )
#define CPP_REAL r8
#define CPP_INT  i8
#else
#define CPP_REAL r4
#define CPP_INT  i4
#endif

      implicit none
      integer, private      :: nblks
      parameter  (nblks   = 100)

      character*20, private :: blkname(nblks)



      integer (CPP_INT), private      :: tblk

      real (CPP_REAL), private       :: etime
      real (CPP_REAL), private       :: totim
      real (CPP_REAL), private       :: tarray(2)
      type tms
           private
           real (CPP_REAL) :: usr, sys
      end type tms


      type (tms), private   :: accum(nblks), last(nblks)

      real (CPP_REAL), private       :: us_tmp1(nblks,2)
      real (CPP_REAL), private       :: us_tmp2(nblks,2)

      contains
         subroutine timing_init
!
! init
!
         implicit none

# if defined( CRAY ) || defined( CRAY_T3E )
         real (CPP_REAL)   :: real8c, real8r
# endif
      
         integer (CPP_INT) :: C, R, M
         real (CPP_REAL)   :: wclk

         integer  n


         tblk=0
         do n = 1, nblks
            accum(n)%usr = 0.
            accum(n)%sys = 0.
            last(n)%usr  = 0.
            last(n)%sys  = 0.
         end do
!
! ... To reduce the overhead for the first call
!

# if defined( IRIX64 ) || ( defined FFC )
         totim = etime(tarray)
# else
! FastOpt begin replace named arguments
!        CALL SYSTEM_CLOCK(Count=C, Count_Rate=R, Count_Max=M)
         CALL SYSTEM_CLOCK(C, R, M)
! FastOpt end

#  if defined( CRAY ) || defined( CRAY_T3E )
         real8c = C
         real8r = R
         wclk   = real8c / real8r
         totim  = wclk
#  else
         wclk =  REAL(C) / REAL(R)
         totim = wclk
#  endif
# endif

         end subroutine timing_init


         subroutine timing_on(blk_name)
!
! timing_on
!

         implicit none

         character*(*)  blk_name



         character*20   ctmp
         integer i
         integer iblk
# if defined( CRAY ) || defined( CRAY_T3E )
         real (CPP_REAL)   :: real8c, real8r
# endif

         integer (CPP_INT) :: C, R, M
         real (CPP_REAL)   :: wclk

         integer ierr


!         call upper(blk_name)
!c         ctmp=blk_name(:len_trim(blk_name))
         ctmp=trim(blk_name)

!         write(*,*) 'timing_on ', ctmp
         iblk=0
         do i=1, tblk
            if ( ctmp .EQ. blkname(i) ) then
               iblk =i
            endif
         enddo
      
         if ( iblk .eq. 0 ) then

# if defined ( SPMD )
            if (gid .eq. 0 ) then
# endif
               tblk=tblk+1
# if defined ( SPMD )
            endif
            call mp_bcst_int(tblk)
!            write(*,*) 'tblk ', tblk, 'gid ', gid, ctmp
# endif
            iblk=tblk
!            call upper(blk_name)
!C            blkname(iblk)=blk_name(:len_trim(blk_name))
            blkname(iblk)=trim(blk_name)
        endif

# if defined( IRIX64 ) || ( defined FFC )
        totim = etime(tarray)
        last(iblk)%usr = tarray(1)
        last(iblk)%sys = tarray(2)
# else
! FastOpt begin replace named arguments
!        CALL SYSTEM_CLOCK(Count=C, Count_Rate=R, Count_Max=M)
         CALL SYSTEM_CLOCK(C, R, M)
! FastOpt end

#  if defined( CRAY ) || defined( CRAY_T3E )
        real8c = C
        real8r = R
        wclk   = real8c / real8r
#  else
        wclk = REAL(C) / REAL(R)
#  endif
        last(iblk)%usr = wclk
        last(iblk)%sys = 0.0
# endif

        end subroutine timing_on


        subroutine timing_off(blk_name)
!
! Timing_off
!

        implicit none
        character*(*) blk_name

        character*20  ctmp
        integer i
# if defined( CRAY ) || defined( CRAY_T3E )
        real (CPP_REAL)   :: real8c, real8r
# endif

        integer (CPP_INT) :: C, R, M
        real (CPP_REAL)   :: wclk

        integer  iblk


!        call upper(blk_name)
!v        ctmp=blk_name(:len_trim(blk_name))
        ctmp=trim(blk_name)

        iblk=0
        do i=1, tblk
           if ( ctmp .EQ. blkname(i) ) then
              iblk =i
           endif
        enddo
      
!         write(*,*) 'timing_off ', ctmp, tblk, tblk
        if ( iblk .eq. 0 ) then
!SJL           write(*,*) ' timing_on missing ', ctmp, 'gid', gid
!           write(*,*) 'stop in timing off in ', ctmp
!           stop 
        endif


# if defined( IRIX64 ) || ( defined FFC ) 
        totim = etime(tarray)
        accum(iblk)%usr = accum(iblk)%usr +           &
                        tarray(1) - last(iblk)%usr
        accum(iblk)%sys = accum(iblk)%sys +           &
                        tarray(2) - last(iblk)%sys
        last(iblk)%usr = tarray(1)
        last(iblk)%sys = tarray(2)
# else
! FastOpt begin replace named arguments
!        CALL SYSTEM_CLOCK(Count=C, Count_Rate=R, Count_Max=M)
         CALL SYSTEM_CLOCK(C, R, M)
! FastOpt end

#  if defined( CRAY ) || defined( CRAY_T3E )
        real8c = C
        real8r = R
        wclk   = real8c / real8r
#  else
        wclk = REAL(C) / REAL(R)
#  endif
        accum(iblk)%usr = accum(iblk)%usr + wclk - last(iblk)%usr
        accum(iblk)%sys = 0.0
        last(iblk)%usr  = wclk
        last(iblk)%sys  = 0.0
# endif

        end subroutine timing_off



        subroutine timing_prt
!
! Timing_prt
!
        implicit none

        integer  n

        type (tms)   :: others, tmp(nblks)
        real         :: tmpmax(1)

#if defined( SPMD )

        do n = 1, nblks                   !will clean these later
           tmpmax(1) = accum(n)%usr
           call mp_reduce_max(1, tmpmax(1))
           tmp(n)%usr = tmpmax(1)
           tmpmax(1) = accum(n)%sys
           call mp_reduce_max(1, tmpmax(1))
           tmp(n)%sys = tmpmax(1)
        enddo
        if ( gid .eq. 0 ) then
#else
        do n = 1, nblks
           tmp(n)%usr = accum(n)%usr
           tmp(n)%sys = accum(n)%sys
        enddo
#endif

        print *
        print *,                                  &
        '  -----------------------------------------------------'
        print *,                                  &
        '     Block          User time  System Time   Total Time'
        print *,                                  &
        '  -----------------------------------------------------'

        do n = 1, tblk
           print '(3x,a20,2x,3(1x,f12.4))', blkname(n),     &
               tmp(n)%usr, tmp(n)%sys, tmp(n)%usr + tmp(n)%sys
        end do


        print *
#if defined( SPMD )
        print *
        endif
#endif

        end subroutine timing_prt

#else  /* endif of (defined TIMING) */

      contains
        subroutine timing_init
        end subroutine timing_init

        subroutine timing_on ( a )
        character(len=*) :: a
        end subroutine timing_on

        subroutine timing_off ( a )
        character(len=*) :: a
        end subroutine timing_off

        subroutine timing_prt
        end subroutine timing_prt

#endif  /* endif of (defined TIMING) */

      end module TimingModule
