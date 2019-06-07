!-----------------------------------------------------------------------
!BOP
! !ROUTINE: writerst --- Write the restart file
!
! !INTERFACE:
      subroutine writerst(iuic, dflnm, im, jm, km, nq, jfirst, jlast,
     .                    nstep, nymd, nhms, ps, delp, u, v, pt, q,
     .                    ng_s, ng_d)
! !USES:

#if defined (SPMD)
      use mod_comm
#endif

      implicit none

#if defined( SPMD )
#define CPP_PS   pstmp
#define CPP_U    utmp
#define CPP_V    vtmp
#define CPP_DELP delptmp
#define CPP_PT   pttmp
#define CPP_Q    qtmp
#else
#define CPP_PS   ps
#define CPP_U    u
#define CPP_V    v
#define CPP_DELP delp
#define CPP_PT   pt
#define CPP_Q    q
#endif

! !INPUT PARAMETERS:
      integer iuic                         ! Unit number
      character*(*) dflnm                  ! dynamics restart filename
      integer im, jm, km, nq               ! Dimensions of arrays
      integer jfirst, jlast                ! Latitude strip
      integer nymd, nhms, nstep            ! Dates 
      integer ng_s, ng_d                   ! info for ghost grids

      real :: ps(im,jfirst:jlast) 
      real :: u(im,jfirst-ng_d:jlast+ng_s,km) 
      real :: v(im,jfirst-ng_s:jlast+ng_d,km) 
      real :: pt(im,jfirst-ng_d:jlast+ng_d,km) 
      real :: delp(im,jfirst:jlast,km) 
      real :: q(im,jfirst-ng_d:jlast+ng_d,km,nq)

! !REVISION HISTORY:
!   10Feb2006 Todling change to be compatible w/ Oloso changes to readrst routine;
!                     rst file assumed to be open

#if defined( SPMD )
      real pstmp(im*jm), utmp(im*jm*km), vtmp(im*jm*km)
      real delptmp(im*jm*km), pttmp(im*jm*km), qtmp(im,jm,km,nq) 
      integer i, j, k, l, ij, n

! --- Gather from all subdomains to GID=0 ---
      call mp_barrier
      call mp_gather4d(   u,    utmp, im, jm, km,  1, jfirst, jlast, 1, km, ng_d, ng_s, 0)
      call mp_gather4d(   v,    vtmp, im, jm, km,  1, jfirst, jlast, 1, km, ng_s, ng_d, 0)
      call mp_gather4d(  pt,   pttmp, im, jm, km,  1, jfirst, jlast, 1, km, ng_d, ng_d, 0)
      call mp_gather4d(delp, delptmp, im, jm, km,  1, jfirst, jlast, 1, km,    0,    0, 0)
      call mp_gather4d(  ps,   pstmp, im, jm,  1,  1, jfirst, jlast, 1,  1,    0,    0, 0)
      call mp_gather4d(   q,    qtmp, im, jm, km, nq, jfirst, jlast, 1, km, ng_d, ng_d, 0)

      if (gid == 0) then
#endif
!_RT     open(iuic,file=dflnm,form='unformatted',access='sequential')
         rewind(iuic)
         write(iuic) nstep, nymd, nhms
         write(iuic) CPP_PS, CPP_DELP, CPP_U, CPP_V, CPP_PT
         if(nq .ne. 0) write(iuic) CPP_Q
         close(iuic)
         write(*,*) "Dynamic restart file written and closed"

#if defined( SPMD )
      endif
#endif
      return
!EOC
      end subroutine writerst
!---------------------------------------------------------------------
