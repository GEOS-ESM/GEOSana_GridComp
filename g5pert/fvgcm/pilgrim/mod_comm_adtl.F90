#ifdef SPMD

!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOI
!
! !TITLE: MPI Adjoint-Related Interfaces
!
! !AUTHORS: Ralf Giering and Thomas Kaminski
!
! !AFFILIATION:  FastOpt, Germany
!
! !DATE: Feb 2002
!
! !INTRODUCTION: Description of adjoint-related MPI utilities
!
! \begin{verbatim}
!
! TAF flow directives and hand written 
! tangent linear and adjoint subroutines to relevant
! communication subroutines of the NASA-NCAR fvgcm
!
! tangent linear routines: g_<routine name>
! adjoint routines       : ad_<routine name>
! routines treated       : mp_send3d_ns
!                        : mp_recv3d_ns
!                        : mp_send2_s
!                        : mp_recv2_s
!                        : mp_send2_n
!                        : mp_recv2_n
!                        : mp_send2_ns
!                        : mp_recv2_ns
!                        : mp_send_s
!                        : mp_recv_n
!
! mp_send4d_ns  --  mp_recv4d_ns
! mp_send3d_ns  --  mp_recv3d_ns
! mp_send2_ns   --  mp_recv2_ns
! 
! mp_send2_s  --  mp_recv2_n
! mp_send2_n  --  mp_recv2_s
!
! mp_send_s   --  mp_recv_n
! mp_send_pe  --  mp_recv_pe
! mp_send_ua  --  mp_recv_ua
!
! \end{verbatim}
!
!EOI
!-----------------------------------------------------------------------
!*********************************************************************
#define DEBUG_COM_NO

!*********************************************************************
! mp_send2_ns
!*********************************************************************
!$taf module mod_comm subroutine mp_send2_ns   input = 1,2,3,4,5,6,7,8,9
!$taf module mod_comm subroutine mp_send2_ns  output =
!$taf module mod_comm subroutine mp_send2_ns  active = 8,9
!$taf module mod_comm subroutine mp_send2_ns  depend = 1,2,3,4,5,6,7
!$taf module mod_comm subroutine mp_send2_ns  adname = mod_comm::mp_send2_ns_ad
!$taf module mod_comm subroutine mp_send2_ns ftlname = mod_comm::mp_send2_ns
!$taf module mod_comm subroutine mp_send2_ns common mp_2_ns output = 1
!$taf module mod_comm subroutine mp_send2_ns common mp_2_ns active = 1

!*********************************************************************
! mp_recv2_ns
!*********************************************************************
!$taf module mod_comm subroutine mp_recv2_ns   input = 1,2,3,4,5,6,7
!$taf module mod_comm subroutine mp_recv2_ns  output = 8,9
!$taf module mod_comm subroutine mp_recv2_ns  active = 8,9
!$taf module mod_comm subroutine mp_recv2_ns  depend = 1,2,3,4,5,6,7
!$taf module mod_comm subroutine mp_recv2_ns  adname = mod_comm::mp_recv2_ns_ad
!$taf module mod_comm subroutine mp_recv2_ns ftlname = mod_comm::mp_recv2_ns
!$taf module mod_comm subroutine mp_recv2_ns common mp_2_ns  input = 1
!$taf module mod_comm subroutine mp_recv2_ns common mp_2_ns active = 1


!*********************************************************************
! mp_send3d_ns
!*********************************************************************
!$taf module mod_comm subroutine mp_send3d_ns   input = 1,2,3,4,5,6,7,8,9,10
!$taf module mod_comm subroutine mp_send3d_ns  output =
!$taf module mod_comm subroutine mp_send3d_ns  active = 9
!$taf module mod_comm subroutine mp_send3d_ns  depend = 1,2,3,4,5,6,7,8  ,10
!$taf module mod_comm subroutine mp_send3d_ns  adname = mod_comm::mp_send3d_ns_ad
!$taf module mod_comm subroutine mp_send3d_ns ftlname = mod_comm::mp_send3d_ns
!$taf module mod_comm subroutine mp_send3d_ns common mp_3d_ns output = 1
!$taf module mod_comm subroutine mp_send3d_ns common mp_3d_ns active = 1

!*********************************************************************
! mp_recv3d_ns
!*********************************************************************
!$taf module mod_comm subroutine mp_recv3d_ns   input = 1,2,3,4,5,6,7,8  ,10
!$taf module mod_comm subroutine mp_recv3d_ns  output = 9
!$taf module mod_comm subroutine mp_recv3d_ns  active = 9
!$taf module mod_comm subroutine mp_recv3d_ns  depend = 1,2,3,4,5,6,7,8  ,10
!$taf module mod_comm subroutine mp_recv3d_ns  adname = mod_comm::mp_recv3d_ns_ad
!$taf module mod_comm subroutine mp_recv3d_ns ftlname = mod_comm::mp_recv3d_ns
!$taf module mod_comm subroutine mp_recv3d_ns common mp_3d_ns  input = 1
!$taf module mod_comm subroutine mp_recv3d_ns common mp_3d_ns active = 1

!*********************************************************************
! mp_send3d_ns2
! mp_recv3d_ns2
!*********************************************************************
! we need these wrappers to include three call in cd_core_ad
! that would otherwise be missing
!*********************************************************************
!$taf module mod_comm subroutine mp_send3d_ns2   input = 1,2,3,4,5,6,7,8,9,10
!$taf module mod_comm subroutine mp_send3d_ns2  output =
!$taf module mod_comm subroutine mp_send3d_ns2  active = 9
!$taf module mod_comm subroutine mp_send3d_ns2  depend = 1,2,3,4,5,6,7,8  ,10
!$taf module mod_comm subroutine mp_send3d_ns2  adname = mp_send3d_ns2_ad
!$taf module mod_comm subroutine mp_send3d_ns2 ftlname = mp_send3d_ns2
!$taf module mod_comm subroutine mp_send3d_ns2 common mp_3d_ns2 output = 1
!$taf module mod_comm subroutine mp_send3d_ns2 common mp_3d_ns2 active = 1

!$taf module mod_comm subroutine mp_recv3d_ns2   input = 1,2,3,4,5,6,7,8  ,10
!$taf module mod_comm subroutine mp_recv3d_ns2  output = 9
!$taf module mod_comm subroutine mp_recv3d_ns2  active = 9
!$taf module mod_comm subroutine mp_recv3d_ns2  depend = 1,2,3,4,5,6,7,8  ,10
!$taf module mod_comm subroutine mp_recv3d_ns2  adname = mp_recv3d_ns2_ad
!$taf module mod_comm subroutine mp_recv3d_ns2 ftlname = mp_recv3d_ns2
!$taf module mod_comm subroutine mp_recv3d_ns2 common mp_3d_ns2  input = 1
!$taf module mod_comm subroutine mp_recv3d_ns2 common mp_3d_ns2 active = 1

!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: mp_send3d_ns2_ad --- Adjoint of mp_send3d_ns2
!
! !INTERFACE:
 
subroutine mp_send3d_ns2_ad(im,jm,jfirst,jlast,kfirst,klast,ng_s,ng_n,q,iq)

! !USES:

  use mod_comm, only: mp_send3d_ns_ad
  implicit none

! !INPUT PARAMETERS:
!   to do
! !OUTPUT PARAMETERS:
!   to do
!
! !DESCRIPTION: Adjoint of mp\_send3d\_ns2
!
! !REVISION HISTORY:
!
!  01Jun2002  FastOpt   Initial code.
!  13Nov2002  Todling   Added Protex prologue.
!
!EOP
!-----------------------------------------------------------------------
  integer im, jm
  integer jfirst, jlast
  integer kfirst, klast
  integer ng_s      ! southern zones to ghost 
  integer ng_n      ! northern zones to ghost 
  integer iq        ! Counter
  real    q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)
  call mp_send3d_ns_ad(im,jm,jfirst,jlast,kfirst,klast,ng_s,ng_n,q,iq)
end subroutine mp_send3d_ns2_ad

!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: mp_recv3d_ns2_ad --- Adjoint of mp_recv3d_ns2
!
! !INTERFACE:
 
subroutine mp_recv3d_ns2_ad(im,jm,jfirst,jlast,kfirst,klast,ng_s,ng_n,q,iq)

! !USES:

  use mod_comm, only: mp_recv3d_ns_ad
  implicit none

! !INPUT PARAMETERS:

  integer im, jm
  integer jfirst, jlast
  integer kfirst, klast
  integer ng_s      ! southern zones to ghost 
  integer ng_n      ! northern zones to ghost 

! !OUTPUT PARAMETERS:
 
! !DESCRIPTION: Adjoint of mp\_recv3d\_ns2
!
! !REVISION HISTORY:
!
!  01Jun2002  FastOpt   Initial code.
!  13Nov2002  Todling   Added Protex prologue.
!
!EOP
!-----------------------------------------------------------------------
  real    q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)
  integer iq
  call mp_recv3d_ns_ad(im,jm,jfirst,jlast,kfirst,klast,ng_s,ng_n,q,iq)
end subroutine mp_recv3d_ns2_ad


!*********************************************************************
! mp_send4d_ns
!*********************************************************************
!$taf module mod_comm subroutine mp_send4d_ns   input = 1,2,3,4,5,6,7,8,9,10
!$taf module mod_comm subroutine mp_send4d_ns  output =
!$taf module mod_comm subroutine mp_send4d_ns  active = 10
!$taf module mod_comm subroutine mp_send4d_ns  depend = 1,2,3,4,5,6,7,8,9
!$taf module mod_comm subroutine mp_send4d_ns  adname = mod_comm::mp_send4d_ns_ad
!$taf module mod_comm subroutine mp_send4d_ns ftlname = mod_comm::mp_send4d_ns
!$taf module mod_comm subroutine mp_send4d_ns common mp_4d_ns output = 1
!$taf module mod_comm subroutine mp_send4d_ns common mp_4d_ns active = 1

!*********************************************************************
! mp_recv4d_ns
!*********************************************************************
!$taf module mod_comm subroutine mp_recv4d_ns   input = 1,2,3,4,5,6,7,8,9
!$taf module mod_comm subroutine mp_recv4d_ns  output = 10
!$taf module mod_comm subroutine mp_recv4d_ns  active = 10
!$taf module mod_comm subroutine mp_recv4d_ns  depend = 1,2,3,4,5,6,7,8,9
!$taf module mod_comm subroutine mp_recv4d_ns  adname = mod_comm::mp_recv4d_ns_ad
!$taf module mod_comm subroutine mp_recv4d_ns ftlname = mod_comm::mp_recv4d_ns
!$taf module mod_comm subroutine mp_recv4d_ns common mp_4d_ns  input = 1
!$taf module mod_comm subroutine mp_recv4d_ns common mp_4d_ns active = 1


!*********************************************************************
! mp_send2_s
!*********************************************************************
!$taf module mod_comm subroutine mp_send2_s   input = 1,2,3,4,5,6,7,8,9,10
!$taf module mod_comm subroutine mp_send2_s  output =
!$taf module mod_comm subroutine mp_send2_s  active = 9,10
!$taf module mod_comm subroutine mp_send2_s  depend = 1,2,3,4,5,6,7,8
!$taf module mod_comm subroutine mp_send2_s  adname = mod_comm::mp_send2_s_ad
!$taf module mod_comm subroutine mp_send2_s ftlname = mod_comm::mp_send2_s
!$taf module mod_comm subroutine mp_send2_s common mp_2_s output = 1
!$taf module mod_comm subroutine mp_send2_s common mp_2_s active = 1

!*********************************************************************
! mp_send2_n
!*********************************************************************
!$taf module mod_comm subroutine mp_send2_n   input = 1,2,3,4,5,6,7,8,9,10
!$taf module mod_comm subroutine mp_send2_n  output =
!$taf module mod_comm subroutine mp_send2_n  active = 9,10
!$taf module mod_comm subroutine mp_send2_n  depend = 1,2,3,4,5,6,7,8
!$taf module mod_comm subroutine mp_send2_n  adname = mod_comm::mp_send2_n_ad
!$taf module mod_comm subroutine mp_send2_n ftlname = mod_comm::mp_send2_n
!$taf module mod_comm subroutine mp_send2_n common mp_2_n output = 1
!$taf module mod_comm subroutine mp_send2_n common mp_2_n active = 1

!*********************************************************************
! mp_recv2_s
!*********************************************************************
!$taf module mod_comm subroutine mp_recv2_s   input = 1,2,3,4,5,6,7,8
!$taf module mod_comm subroutine mp_recv2_s  output = 9,10
!$taf module mod_comm subroutine mp_recv2_s  active = 9,10
!$taf module mod_comm subroutine mp_recv2_s  depend = 1,2,3,4,5,6,7,8
!$taf module mod_comm subroutine mp_recv2_s  adname = mod_comm::mp_recv2_s_ad
!$taf module mod_comm subroutine mp_recv2_s ftlname = mod_comm::mp_recv2_s
!$taf module mod_comm subroutine mp_recv2_s common mp_2_n  input = 1
!$taf module mod_comm subroutine mp_recv2_s common mp_2_n active = 1

!*********************************************************************
! mp_recv2_n
!*********************************************************************
!$taf module mod_comm subroutine mp_recv2_n   input = 1,2,3,4,5,6,7,8
!$taf module mod_comm subroutine mp_recv2_n  output = 9,10
!$taf module mod_comm subroutine mp_recv2_n  active = 9,10
!$taf module mod_comm subroutine mp_recv2_n  depend = 1,2,3,4,5,6,7,8
!$taf module mod_comm subroutine mp_recv2_n  adname = mod_comm::mp_recv2_n_ad
!$taf module mod_comm subroutine mp_recv2_n ftlname = mod_comm::mp_recv2_n
!$taf module mod_comm subroutine mp_recv2_n common mp_2_s  input = 1
!$taf module mod_comm subroutine mp_recv2_n common mp_2_s active = 1



!*********************************************************************
! mp_send_s
!*********************************************************************
!$taf module mod_comm subroutine mp_send_s   input = 1,2,3,4,5,6,7,8,9
!$taf module mod_comm subroutine mp_send_s  output =
!$taf module mod_comm subroutine mp_send_s  active = 9
!$taf module mod_comm subroutine mp_send_s  depend = 1,2,3,4,5,6,7,8
!$taf module mod_comm subroutine mp_send_s  adname = mod_comm::mp_send_s_ad
!$taf module mod_comm subroutine mp_send_s ftlname = mod_comm::mp_send_s
!$taf module mod_comm subroutine mp_send_s common mp_sn output = 1
!$taf module mod_comm subroutine mp_send_s common mp_sn active = 1

!*********************************************************************
! mp_recv_n
!*********************************************************************
!$taf module mod_comm subroutine mp_recv_n   input = 1,2,3,4,5,6,7,8
!$taf module mod_comm subroutine mp_recv_n  output = 9
!$taf module mod_comm subroutine mp_recv_n  active = 9
!$taf module mod_comm subroutine mp_recv_n  depend = 1,2,3,4,5,6,7,8
!$taf module mod_comm subroutine mp_recv_n  adname = mod_comm::mp_recv_n_ad
!$taf module mod_comm subroutine mp_recv_n ftlname = mod_comm::mp_recv_n
!$taf module mod_comm subroutine mp_recv_n common mp_sn  input = 1
!$taf module mod_comm subroutine mp_recv_n common mp_sn active = 1



!*********************************************************************
! mp_send_pe
!*********************************************************************
!$taf module mod_comm subroutine mp_send_pe   input = 1,2,3,4,5,6,7
!$taf module mod_comm subroutine mp_send_pe  output =
!$taf module mod_comm subroutine mp_send_pe  active = 7
!$taf module mod_comm subroutine mp_send_pe  depend = 1,2,3,4,5,6
!$taf module mod_comm subroutine mp_send_pe  adname = mod_comm::mp_send_pe_ad
!$taf module mod_comm subroutine mp_send_pe ftlname = mod_comm::mp_send_pe
!$taf module mod_comm subroutine mp_send_pe common mp_pe output = 1
!$taf module mod_comm subroutine mp_send_pe common mp_pe active = 1

!*********************************************************************
! mp_recv_pe
!*********************************************************************
!$taf module mod_comm subroutine mp_recv_pe   input = 1,2,3,4,5,6
!$taf module mod_comm subroutine mp_recv_pe  output = 7
!$taf module mod_comm subroutine mp_recv_pe  active = 7
!$taf module mod_comm subroutine mp_recv_pe  depend = 1,2,3,4,5,6
!$taf module mod_comm subroutine mp_recv_pe  adname = mod_comm::mp_recv_pe_ad
!$taf module mod_comm subroutine mp_recv_pe ftlname = mod_comm::mp_recv_pe
!$taf module mod_comm subroutine mp_recv_pe common mp_pe  input = 1
!$taf module mod_comm subroutine mp_recv_pe common mp_pe active = 1


!*********************************************************************
! mp_send_ua
!*********************************************************************
!$taf module mod_comm subroutine mp_send_ua   input = 1,2,3,4,5,6,7
!$taf module mod_comm subroutine mp_send_ua  output =
!$taf module mod_comm subroutine mp_send_ua  active = 7
!$taf module mod_comm subroutine mp_send_ua  depend = 1,2,3,4,5,6
!$taf module mod_comm subroutine mp_send_ua  adname = mod_comm::mp_send_ua_ad
!$taf module mod_comm subroutine mp_send_ua ftlname = mod_comm::mp_send_ua
!$taf module mod_comm subroutine mp_send_ua common mp_ua output = 1
!$taf module mod_comm subroutine mp_send_ua common mp_ua active = 1

!*********************************************************************
! mp_recv_ua
!*********************************************************************
!$taf module mod_comm subroutine mp_recv_ua   input = 1,2,3,4,5,6
!$taf module mod_comm subroutine mp_recv_ua  output = 7
!$taf module mod_comm subroutine mp_recv_ua  active = 7
!$taf module mod_comm subroutine mp_recv_ua  depend = 1,2,3,4,5,6
!$taf module mod_comm subroutine mp_recv_ua  adname = mod_comm::mp_recv_ua_ad
!$taf module mod_comm subroutine mp_recv_ua ftlname = mod_comm::mp_recv_ua
!$taf module mod_comm subroutine mp_recv_ua common mp_ua  input = 1
!$taf module mod_comm subroutine mp_recv_ua common mp_ua active = 1


!*********************************************************************
! mp_reduce_sum
! is linear and self adjoint
!*********************************************************************
!$taf module mod_comm subroutine mp_reduce_sum   input = 1
!$taf module mod_comm subroutine mp_reduce_sum  output = 1
!$taf module mod_comm subroutine mp_reduce_sum  active = 1
!$taf module mod_comm subroutine mp_reduce_sum  depend =
!$taf module mod_comm subroutine mp_reduce_sum  adname = mod_comm::mp_reduce_sum


!*********************************************************************
! mp_sum1d
!*********************************************************************
!$taf module mod_comm subroutine mp_sum1d   input = 1,2,3,4
!$taf module mod_comm subroutine mp_sum1d  output =         5
!$taf module mod_comm subroutine mp_sum1d  active =       4,5
!$taf module mod_comm subroutine mp_sum1d  depend = 1,2,3
!$taf module mod_comm subroutine mp_sum1d  adname = admp_sum1d
!$taf module mod_comm subroutine mp_sum1d ftlname = mod_comm::mp_sum1d
!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: admp_sum1d --- Adjoint of mp_reduce_sum
!
! !INTERFACE:

      subroutine admp_sum1d( jm, jfirst, jlast, adqin, adsum0 )

! !USES:

      use mod_comm, only: gid
      use mod_comm, only: mp_reduce_sum
      implicit none

! !INPUT PARAMETERS:

      integer jm, jfirst, jlast

! !INPUT/OUTPUT PARAMETERS:

      real  adqin(jfirst:jlast)
      real  adsum0
 
! !DESCRIPTION: Adjoint of mp\_reduce\_sum 
!
! !REVISION HISTORY:
!
!  01Jun2002  FastOpt   Initial code.
!  13Nov2002  Todling   Added Protex prologue.
!
!EOP
!-----------------------------------------------------------------------
      real  sum0
      integer j

#ifdef DEBUG_COM
       print *, 'admp_sum1d: ', gid
       call flush(6)
#endif
      call mp_reduce_sum( adsum0 )
      do j = jfirst, jlast
	adqin(j) = adqin(j) + adsum0
      enddo
      adsum0 = 0.0

      end subroutine admp_sum1d

!*********************************************************************
! mp_reduce_max
!*********************************************************************
!$taf module mod_comm subroutine mp_reduce_max   input = 1,2
!$taf module mod_comm subroutine mp_reduce_max  output =   2
!$taf module mod_comm subroutine mp_reduce_max  active =   2
!$taf module mod_comm subroutine mp_reduce_max  depend = 1,2
!$taf module mod_comm subroutine mp_reduce_max  adname = admp_reduce_max
!$taf module mod_comm subroutine mp_reduce_max ftlname = g_mp_reduce_max
!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: g_mp_reduce_max --- Tangent of mp_reduce_max
!
! !INTERFACE:

      subroutine g_mp_reduce_max(km, cymax, g_cymax)

! !USES:

      use mod_comm, only: gid, commglobal

      implicit none

#if !defined(USE_MLP)
#include "mpif.h"
#define mp_2precision MPI_2DOUBLE_PRECISION
#endif

! !DESCRIPTION: Tangent of mp_reduce_max
!
! !REVISION HISTORY:
!
!  01Jun2002  FastOpt   Initial code.
!  13Nov2002  Todling   Added Protex prologue.
!
!EOP
!-----------------------------------------------------------------------
      integer km
      real      cymax(km)
      real    g_cymax(km)

      real maxin (2,km)
      real maxout(2,km)
      integer k
      integer ierror

#ifdef DEBUG_COM
       print *, 'g_mp_reduce_max: ', gid
       call flush(6)
#endif

#if !defined(USE_MLP)
!$omp parallel do private(k)
      do k=1,km
        maxin(1,k) = cymax(k)
        maxin(2,k) = gid
      enddo
			  
      call MPI_Allreduce( maxin, maxout, km, mp_2precision, MPI_MAXLOC, &
                          commglobal, ierror )

      do k=1,km
        cymax(k) = maxout(1,k)
        if (gid .ne. int(maxout(2,k)) ) then
          g_cymax(k) = 0.
	endif
      enddo
			  
#else
	stop 'g_mp_reduce_max: USE_MLP not implemented'
#endif
      end subroutine g_mp_reduce_max

!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: admp_reduce_max --- Adjoint of
!
! !INTERFACE:

      subroutine admp_reduce_max(km, cymax, adcymax)

! !USES:

      use mod_comm, only: gid, commglobal

      implicit none

#if !defined(USE_MLP)
#include "mpif.h"
#define mp_2precision MPI_2DOUBLE_PRECISION
#endif

! !DESCRIPTION: Adjoint of
!
! !REVISION HISTORY:
!
!  01Jun2002  FastOpt   Initial code.
!  13Nov2002  Todling   Added Protex prologue.
!
!EOP
!-----------------------------------------------------------------------
      integer km
      real      cymax(km)
      real    adcymax(km)

      real    maxin (2,km)
      real    maxout(2,km)
      real    adcysum(km)
      integer k
      integer ierror

#ifdef DEBUG_COM
       print *, 'admp_reduce_max: ', gid
       call flush(6)
#endif

#if !defined(USE_MLP)
      call MPI_Allreduce( adcymax, adcysum, km, MPI_DOUBLE_PRECISION, MPI_SUM, &
                          commglobal, ierror )
!$omp parallel do private(k)
      do k=1,km
        maxin(1,k) = cymax(k)
        maxin(2,k) = gid
      enddo
			  
      call MPI_Allreduce( maxin, maxout, km, mp_2precision, MPI_MAXLOC, &
                          commglobal, ierror )

      do k=1,km
        adcymax(k) = 0.
        if (gid .eq. int(maxout(2,k)) ) then
          adcymax(k) = adcysum(k)
	endif
      enddo
#else
	stop 'admp_reduce_max: USE_MLP not implemented'
#endif

      end subroutine admp_reduce_max

!*********************************************************************
! La_2_Ga
!*********************************************************************
!$taf module mod_comm subroutine La_2_Ga   input = 1  ,3,4,5,6,7,8,9
!$taf module mod_comm subroutine La_2_Ga  output =   2
!$taf module mod_comm subroutine La_2_Ga  active = 1,2
!$taf module mod_comm subroutine La_2_Ga  depend =     3,4,5,6,7,8,9
!$taf module mod_comm subroutine La_2_Ga  adname = adLa_2_Ga
!$taf module mod_comm subroutine La_2_Ga ftlname = mod_comm::La_2_Ga
!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: adLa_2_Ga --- Adjoint of
!
! !INTERFACE:

subroutine adLa_2_Ga(addata_l, addata_g, im, jm, km,&
                     jfirstm, jlastp, kfirst, klast )
! !USES:

  use mod_comm, only: gid
  implicit none

! !DESCRIPTION: Adjoint of
!
! !REVISION HISTORY:
!
!  01Jun2002  FastOpt   Initial code.
!  13Nov2002  Todling   Added Protex prologue.
!
!EOP
!-----------------------------------------------------------------------
  integer :: i, j, k
  integer :: im, jm, km
  integer :: jfirstm, jlastp, kfirst, klast
  real    :: addata_l(1:im, jfirstm:jlastp, kfirst:klast)
  real    :: addata_g(im, jm, km)

  stop 'adLa_2_Ga not needed yet'

end subroutine adLa_2_Ga


!*********************************************************************
! mp_barrier
! is empty for MPI but not for USE_MLP !!!
! to be implemented
!*********************************************************************
!$taf module mod_comm subroutine mp_barrier   input =
!$taf module mod_comm subroutine mp_barrier  output =


!*********************************************************************
!                mp_bcst_real
!*********************************************************************
!$taf module mod_comm subroutine mp_bcst_real    input = 1
!$taf module mod_comm subroutine mp_bcst_real   output = 1
!$taf module mod_comm subroutine mp_bcst_real   active = 1
!$taf module mod_comm subroutine mp_bcst_real   depend =
!$taf module mod_comm subroutine mp_bcst_real  ftlname = mod_comm::mp_bcst_real
!$taf module mod_comm subroutine mp_bcst_real   adname = mod_comm::mpi_bcst_real_ad


!*********************************************************************
!                mp_bcst_n_real
!*********************************************************************
!$taf module mod_comm subroutine mp_bcst_n_real    input = 1,2
!$taf module mod_comm subroutine mp_bcst_n_real   output = 1
!$taf module mod_comm subroutine mp_bcst_n_real   active = 1
!$taf module mod_comm subroutine mp_bcst_n_real   depend =   2
!$taf module mod_comm subroutine mp_bcst_n_real  ftlname = mod_comm::mp_bcst_n_real
!$taf module mod_comm subroutine mp_bcst_n_real   adname = mod_comm::mpi_bcst_n_real_ad


!*********************************************************************
!                mp_scatter4d
!*********************************************************************
!$taf module mod_comm subroutine mp_scatter4d   input = 1  ,3,4,5,6,7,8,9
!$taf module mod_comm subroutine mp_scatter4d   input +=    10,11,12,13
!$taf module mod_comm subroutine mp_scatter4d  output =   2
!$taf module mod_comm subroutine mp_scatter4d  active = 1,2
!$taf module mod_comm subroutine mp_scatter4d  depend =     3,4,5,6,7,8,9
!$taf module mod_comm subroutine mp_scatter4d  depend +=    10,11,12,13
!$taf module mod_comm subroutine mp_scatter4d ftlname = mod_comm::mp_scatter4d
!$taf module mod_comm subroutine mp_scatter4d  adname = mod_comm::mp_scatter4d_ad

!*********************************************************************
!                mp_gather4d
!*********************************************************************
!$taf module mod_comm subroutine mp_gather4d   input = 1  ,3,4,5,6,7,8,9
!$taf module mod_comm subroutine mp_gather4d   input +=    10,11,12,13
!$taf module mod_comm subroutine mp_gather4d  output =   2
!$taf module mod_comm subroutine mp_gather4d  active = 1,2
!$taf module mod_comm subroutine mp_gather4d  depend =     3,4,5,6,7,8,9
!$taf module mod_comm subroutine mp_gather4d  depend +=    10,11,12,13
!$taf module mod_comm subroutine mp_gather4d ftlname = mod_comm::mp_gather4d
!$taf module mod_comm subroutine mp_gather4d  adname = mod_comm::mp_gather4d_ad

#else  /* SPMD */

subroutine mpi_dummy_empty()
end subroutine mpi_dummy_empty

#endif /* SPMD */
