      module mod_comm
! !REVISION HISTORY:
!  21Nov2007  Todling  Resolution-free mod_comm
!  06Mar2009  Todling  rst file suffix change from hdf to nc4
!
      use m_dyn, only : dyn_getdim
      use m_inpak90
      implicit none

#if defined ( SPMD ) 
!!    PRIVATE

      logical, public, save :: masterproc

      integer, public, save :: PLON
      integer, public, save :: PLAT
      integer, public, save :: PLEV
      integer, public, save :: PCNST
      integer, public, parameter :: PNATS = 0

      integer, public, save :: imr, jnp, nl, nc

      integer, save :: maxpro       ! Max no. of MLP PE allowed
      integer, save :: max_nq       ! Be carefiul: max_nq = max(nc, 2)
                                    ! nc is the total # of advected tracers

      integer, parameter :: nbuf = 2
      integer, parameter :: nghost = 3

#  if !defined(USE_MLP)
#    include "mpif.h"
#    define mp_precision MPI_DOUBLE_PRECISION

      integer, save :: pkgs_per_pro
      integer, save :: idimsize
      integer, parameter :: max_call = 2*2    ! FastOpt tangent need twice as much
      integer, parameter :: igosouth = 0
      integer, parameter :: igonorth = 1
#    if defined(AIX) && defined(MPI2)
      integer(kind=MPI_ADDRESS_KIND) intptr
      pointer (buff_r_ptr, buff_r)
      pointer (buff_s_ptr, buff_s)
      pointer (buff4d_ptr, buff4d)
      pointer (buff4d_r4_ptr, buff4d_r4)
      pointer (buff3d_i_ptr, buff3d_i)
      real :: buff_r
      real :: buff_s
      real :: buff4d
      real*4 :: buff4d_r4
      integer :: buff3d_i
#    else /* AIX & MPI2 */
      real, allocatable, SAVE:: buff_r(:)
      real, allocatable, SAVE:: buff_s(:)
      real, allocatable, SAVE:: buff4d(:)
      real*4, allocatable, SAVE:: buff4d_r4(:)
      integer, allocatable, SAVE:: buff3d_i(:)
#    endif /* AIX & MPI2 */
      integer, SAVE:: ncall_s
      integer, SAVE:: ncall_r

#    if defined(MPI2)
#      if defined(LINUX)
#        define MPI_ADDRESS_KIND 8
        integer, parameter :: MPI_MODE_NOCHECK             = 0
        integer, parameter :: MPI_MODE_NOSTORE             = 0
        integer, parameter :: MPI_MODE_NOPUT               = 0
        integer, parameter :: MPI_MODE_NOPRECEDE           = 0
        integer, parameter :: MPI_MODE_NOSUCCEED           = 0
#      endif /* LINUX */
      integer(kind=MPI_ADDRESS_KIND) bsize, tdisp
      integer, SAVE:: buffwin      ! Communication window
      integer, SAVE:: buff4dwin    ! Communication window
      integer, SAVE:: buff4d_r4win ! Communication window
      integer, SAVE:: buff3d_iwin  ! Communication window
#    else /* MPI2 */
      integer, SAVE:: tdisp
      integer, SAVE:: nsend                   ! Number of messages out-going
      integer, SAVE:: nrecv                   ! Number of messages in-coming
      integer, SAVE:: nread                   ! Number of messages read
      integer, SAVE, allocatable:: sqest(:)
      integer, SAVE, allocatable:: rqest(:)
#    endif /* MPI2 */

      integer, SAVE,public:: commglobal   ! Global Communicator
      integer, SAVE:: Status(MPI_STATUS_SIZE)
      integer, SAVE, allocatable:: Stats(:)
      integer ierror

#  else /* USE_MLP */

#    if defined (LAHEY)

#      define PTR_INT TRUE
#      define NOT_ASSIGNED

#      include "mlp_ptr.h"

#      undef  PTR_INT
#      undef NOT_ASSIGNED

#    else /* LAHEY */
!
! Main vars:
!
      pointer (wing_4d, g_4d)
      real, allocatable :: g_4d(:, :, :, :)

! Other work arrays:
!
! Type 1: For variables defined at layer edge (wz & pk)
!
      pointer (wing_t1, g_t1)
      real, allocatable :: g_t1(:, :, :, :)

!
! Type 2: For edge pressure (pe)
!
      pointer (wing_t2, g_t2)
      real, allocatable :: g_t2(:, :, :)
!
! Type 3: 
!
      pointer (wing_t3, g_t3)
      real, allocatable :: g_t3(:, :)

!
! General purpose 2D (x-y) array
!
      pointer (wing_2d, g_2d)
      real, allocatable :: g_2d(:,:)

!
! General purpose 1D array
!
      pointer (wing_1d, g_1d)
      real, allocatable :: g_1d(:)
#    endif /* LAHEY */
#  endif /* USE_MLP */

      integer, SAVE, allocatable, public:: numcps(:)
      integer, SAVE, public:: nowpro,numpro,numcpu  
      integer, SAVE, public :: gid
      integer, SAVE, public:: gsize
      integer, allocatable, SAVE, public:: yfirst(:)  ! First latitude
      integer, allocatable, SAVE, public:: ylast(:)   ! Last latitude
      integer, allocatable, SAVE, public:: zfirst(:)  ! First level
      integer, allocatable, SAVE, public:: zlast(:)   ! Last level

      public mp_init, mp_exit, y_decomp, set_decomp
      public mp_reduce_sum
      public mp_send3d_ns_ad
      public mp_recv3d_ns_ad
      public mp_recv3d_ns
      public mp_send3d_ns
      public mp_barrier
      public mp_scatter2d
      public mp_bcst_int
      public mp_scatter4d
      public mp_gather4d
      public mp_bcst_n_real
      public mpi_bcst_n_real_ad
      public mp_recv2_s
      public mp_send2_n
      public mp_recv_n
      public mp_send_s
      public mp_recv_n_ad
      public mp_send_s_ad
      public mp_scatter4d_ad
      public mp_bcst_r2d
      public mp_bcst_real
      public mp_add1d
      public mp_gather4d_ad
      public mp_recv_s
      public mp_send_n
      public mp_send3d_ns2, mp_recv3d_ns2
      public mp_recv2_ns
      public mp_send_pe,mp_recv_pe
      public mp_wrt3d
      public mp_gath_r2d
      public La_2_Ga
      public mp_send2_ns
      public mp_send_pe_ad
      public mp_recv_pe_ad
      public mp_send_ua, mp_recv_ua
      public mp_sum1d
      public mp_reduce_max
      public mp_send4d_ns,mp_recv4d_ns, mp_recv2_n,mp_send2_s
#  if defined(USE_MLP)
      public Ga_Put
      public Ga_Get
#  endif /* USE_MLP */

#  if defined (SEMA)
      integer semid
#  endif /* SEMA */

!.................

      contains

      subroutine mod_comm_init ( )

      implicit none

      character(len=*),  parameter :: rcname = 'fort.811'
      character(len=80) :: rstname, token
      integer :: ntracers, iret
      logical :: exists

!     Get dimensions from restart file
!     --------------------------------
      ierror = 0
      if ( gid==0 ) then

           inquire( file=rcname, exist=exists )
           if ( exists ) then
              call i90_loadf (rcname, iret)
              if( iret .eq. 0) then
                 call I90_label('g4_restart_name:', iret)
                 if (iret .eq. 0) then
                     call I90_Gtoken(token,iret)
                     if( iret .eq. 0) then
                         rstname = trim(token)
                         print *, 'g4 restart filename: ', trim(rstname)
                     endif
                 endif
              endif
              call i90_release ()
           endif

           call dyn_getdim ( trim(rstname), PLON, PLAT, PLEV, PCNST, ierror )
           if(ierror/=0) print *, 'dims from file (im,jm,km,nc): ', PLON, PLAT, PLEV, PCNST

! This is a little hack to allow an rst file w/ many tracers to be
! used in conjunction with a trajectory w/ smaller number of tracers.

           inquire( file=rcname, exist=exists )
           if ( exists ) then
              call i90_loadf (rcname, iret)
              if( iret .eq. 0) then
                 call I90_label('number_of_tracers:', iret)
                 if (iret .eq. 0) then
                     ntracers = I90_GInt(iret)
                     if( iret .eq. 0) PCNST = ntracers
                     print *, 'reset number of tracers to ', ntracers
                 endif
              endif
              call i90_release ()
           endif

      endif ! < gid=0 >

      call mp_bcst_int ( PLON  )
      call mp_bcst_int ( PLAT  )
      call mp_bcst_int ( PLEV  )
      call mp_bcst_int ( PCNST )
      call mp_bcst_int ( ierror)

      if(ierror/=0)then
        print *, 'Trouble reading drst file. Aborting ...'
        call exit(1)
      endif

      imr = PLON
      jnp = PLAT
      nl  = PLEV
      nc  = PCNST + PNATS
      
      maxpro = PLAT/4      ! This is the max 1D decomp
      max_nq = PCNST + 1

#  if !defined(USE_MLP)

      idimsize = PLON*nghost*PLEV*PCNST
#    if defined(AIX) && defined(MPI2)
      integer(kind=MPI_ADDRESS_KIND) intptr
      allocate (buff_r(idimsize*nbuf*max_call))
      allocate (buff_s(idimsize*nbuf*max_call))
      allocate (buff4d(PLON*PLAT*(PLEV+1)*PCNST))
      allocate (buff4d_r4(PLON*PLAT*(PLEV+1)*PCNST))
      allocate (buff3d_i(PLON*PLAT*PLEV))
#    else /* AIX & MPI2 */
      allocate( buff_r(idimsize*nbuf*max_call) )
      allocate ( buff_s(idimsize*nbuf*max_call) )
      allocate ( buff4d(PLON*PLAT*(PLEV+1)*PCNST) )
      allocate ( buff4d_r4(PLON*PLAT*(PLEV+1)*PCNST) )
      allocate ( buff3d_i(PLON*PLAT*PLEV) )
#    endif /* AIX & MPI2 */

#    if defined(MPI2)
#    else /* MPI2 */
      allocate (sqest(nbuf*max_call))
      allocate (rqest(nbuf*max_call))
#    endif /* MPI2 */

      allocate (Stats(nbuf*max_call*MPI_STATUS_SIZE))

#  else /* USE_MLP */

#    if defined (LAHEY)

#    else /* LAHEY */
!
! Main vars:
!
      allocate( g_4d(PLON, PLAT, PLEV, max_nq) )

! Other work arrays:
!
! Type 1: For variables defined at layer edge (wz & pk)
!
      allocate( g_t1(PLON, PLAT, PLEV+1, nbuf))

!
! Type 2: For edge pressure (pe)
!
      allocate( g_t2(PLON, PLEV+1, PLAT) )
!
! Type 3: 
!
      allocate( g_t3(PLEV+PLAT, maxpro) )

!
! General purpose 2D (x-y) array
!
      allocate( g_2d(PLON,PLAT) )

!
! General purpose 1D array
!
      allocate( g_1d(PLAT) )
#    endif /* LAHEY */
#  endif /* USE_MLP */

      allocate (numcps(maxpro))

      end subroutine mod_comm_init

      subroutine mp_init

#  if !defined(USE_MLP)
        integer idimBuff, idimBuff4d
        integer n, nowpro, nowcpu
        integer npthreads
        character*80 evalue
        integer info
        logical flag

        integer mp_size
#if defined (USE_OPENMP)
#if !defined (SET_CPUS)
#if defined (IRIX64)
        integer mp_suggested_numthreads
#else
        integer omp_get_num_threads
#endif
#endif
#endif /* USE_OPENMP */


#if defined(MPI2) && !defined(AIX) && (!defined LINUX)
        call MPI_INITIALIZED( flag, ierror )
        if ( .not. flag ) then
          call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, npthreads, ierror)
        endif
        call MPI_QUERY_THREAD(npthreads, ierror)
#if !defined(MT_OFF)
        if (npthreads == MPI_THREAD_SINGLE) then
          write(*,*) gid, 'Provided MPI_THREAD_SINGLE. ', &
                'Change to MPI_THREAD_MULTIPLE with MPI_INIT_THREAD ', &
                'for multi-threading MPI2'
          call MPI_FINALIZE(ierror)
          stop
        elseif (npthreads == MPI_THREAD_FUNNELED) then
          write(*,*) gid, 'Provided MPI_THREAD_FUNNELED. ', &
                'Change to MPI_THREAD_MULTIPLE with MPI_INIT_THREAD ', &
                'for multi-threading MPI2'
          call MPI_FINALIZE(ierror)
          stop
        elseif (npthreads == MPI_THREAD_SERIALIZED) then
          write(*,*) gid, 'Provided MPI_THREAD_SERIALIZED. ', &
                'Change to MPI_THREAD_MULTIPLE with MPI_INIT_THREAD ', &
                'for multi-threading MPI2'
          call MPI_FINALIZE(ierror)
          stop
        elseif (npthreads == MPI_THREAD_MULTIPLE) then
!          write(*,*) 'Provided MPI_THREAD_MULTIPLE on', gid
        else
          write(*,*) gid,': Error in MPI_INIT_THREAD. ', &
                'Change to MPI_THREAD_MULTIPLE with MPI_INIT_THREAD ', &
                'for multi-threading MPI2'
          call MPI_FINALIZE(ierror)
          stop
        endif
#endif
#else
        call MPI_INITIALIZED( flag, ierror )
        if ( .not. flag ) then
          call MPI_INIT( ierror )
        endif
#endif
        call MPI_COMM_RANK (MPI_COMM_WORLD, gid, ierror)
        call MPI_COMM_SIZE (MPI_COMM_WORLD, numpro, ierror)
        call MPI_COMM_DUP (MPI_COMM_WORLD, commglobal, ierror)

        call mod_comm_init()

#if defined(MPI2)
        call MPI_INFO_CREATE(info, ierror)
        call MPI_INFO_SET(info, "no_locks", "true", ierror)
#if defined(AIX)
        call MPI_TYPE_EXTENT(mp_precision, mp_size, ierror)
        bsize=idimsize*nbuf*max_call*mp_size
        call MPI_ALLOC_MEM(bsize, MPI_INFO_NULL, intptr, ierror)
        buff_r_ptr = intptr
        call MPI_WIN_CREATE(buff_r, bsize, mp_size, info, commglobal, &
                            buffwin, ierror)

        call MPI_ALLOC_MEM(bsize, MPI_INFO_NULL, intptr, ierror)
        buff_s_ptr = intptr

        bsize=PLON*PLAT*(PLEV+1)*PCNST*mp_size
        call MPI_ALLOC_MEM(bsize, MPI_INFO_NULL, intptr, ierror)
        buff4d_ptr = intptr
        call MPI_WIN_CREATE(buff4d, bsize, mp_size, info, commglobal, &
                            buff4dwin, ierror)

        call MPI_TYPE_EXTENT(MPI_REAL, mp_size, ierror)
        bsize=PLON*PLAT*(PLEV+1)*PCNST*mp_size
        call MPI_ALLOC_MEM(bsize, MPI_INFO_NULL, intptr, ierror)
        buff4d_r4_ptr = intptr
        call MPI_WIN_CREATE(buff4d_r4, bsize, mp_size, info, commglobal, &
                            buff4d_r4win, ierror)

        call MPI_TYPE_EXTENT(MPI_INTEGER, mp_size, ierror)
        bsize=PLON*PLAT*PLEV*mp_size
        call MPI_ALLOC_MEM(bsize, MPI_INFO_NULL, intptr, ierror)
        buff3d_i_ptr = intptr
        call MPI_WIN_CREATE(buff3d_i, bsize, mp_size, info, commglobal, &
                            buff3d_iwin, ierror)
#else
        call MPI_TYPE_EXTENT(mp_precision, mp_size, ierror)
        bsize=idimsize*nbuf*max_call
        call MPI_WIN_CREATE(buff_r, bsize, mp_size, info, commglobal, &
                            buffwin, ierror)

        bsize=PLON*PLAT*(PLEV+1)*PCNST
        call MPI_WIN_CREATE(buff4d, bsize, mp_size, info, commglobal, &
                            buff4dwin, ierror)

        bsize=PLON*PLAT*(PLEV+1)*PCNST
        call MPI_TYPE_EXTENT(MPI_REAL, mp_size, ierror)
        call MPI_WIN_CREATE(buff4d_r4, bsize, mp_size, info, commglobal, &
                            buff4d_r4win, ierror)

        bsize=PLON*PLAT*PLEV
        call MPI_TYPE_EXTENT(MPI_INTEGER, mp_size, ierror)
        call MPI_WIN_CREATE(buff3d_i, bsize, mp_size, info, commglobal, &
                            buff3d_iwin, ierror)
#endif
        call MPI_INFO_FREE(info, ierror)
#else
        nsend = 0
        nrecv = 0
        nread = 0
#endif

        ncall_r = 0
        ncall_s = 0

#if defined(SET_CPUS)
        call getenv('NUMBER_CPUS_PER_MLP_PROCESS',evalue)
        if (gid == 0) then
          read(evalue,*) numcpu
        endif
        call mp_bcst_int(numcpu)
#if defined (IRIX64)
       call mp_set_numthreads(numcpu)  !keep it for a while, :)
#else
       call omp_set_num_threads(numcpu)
#endif
#if  defined( IRIX64 ) && defined(PIN_CPUS)
!$omp parallel do private(n,nowcpu)
        nowpro = gid
        do n=1,numcpu
          nowcpu = n + (nowpro) * numcpu-1
          call mp_assign_to_cpu(nowcpu)
        enddo
#endif
#else
#if defined (USE_OPENMP)
#if defined (IRIX64)
        numcpu = mp_suggested_numthreads(0)
#else
        numcpu = omp_get_num_threads()
#endif
#else  /* USE_OPENMP */
        numcpu = 1
#endif /* USE_OPENMP */
#endif

#if defined(MT_OFF)
        pkgs_per_pro = 1
#else
        pkgs_per_pro = numcpu
#endif
#else
        if ( max_nq < PCNST ) then
           write(*,*) "Buffer size for MLP is NOT large enough!"
           stop
        endif

        call gotmem
        call forkit
#if defined (SEMA)
        call semcreate(semid)
#endif

#endif
        allocate( yfirst( numpro ) )
        allocate( ylast( numpro ) )
        allocate( zfirst( numpro ) )
        allocate( zlast( numpro ) )
      end subroutine mp_init

      subroutine mp_exit
        deallocate( yfirst )
        deallocate( ylast )
        deallocate( zfirst )
        deallocate( zlast )
#if !defined(USE_MLP)
#if defined(MPI2)
        call MPI_WIN_FREE( buffwin, ierror )
        call MPI_WIN_FREE( buff4dwin, ierror )
        call MPI_WIN_FREE( buff4d_r4win, ierror )
#endif
        call MPI_FINALIZE (ierror)
#endif
        return
      end subroutine mp_exit


#if defined(USE_MLP)
      subroutine gotmem

#define NOT_ASSIGNED
#include "mlp_ptr.h"
#undef  NOT_ASSIGNED


      integer n_svar
      integer*8 numvar       ! Total # of shared vars     
      parameter (n_svar=100)
      integer*8 isize(n_svar),ipnt(n_svar)
      integer n

      numvar    =  6
      isize(1)  =  PLON*PLAT*PLEV*max_nq
      isize(2)  =  PLON*PLAT*(PLEV+1)*nbuf
      isize(3)  =  PLON*PLAT*(PLEV+1)
      isize(4)  = (PLEV+PLAT)*maxpro
      isize(5)  =  PLON*PLAT
      isize(6)  =  PLAT

      do n=1,numvar
         isize(n) = isize(n) * 8
      enddo

      call mlp_getmem(numvar,isize,ipnt)

      wing_4d  = ipnt(1)
      wing_t1  = ipnt(2)
      wing_t2  = ipnt(3)
      wing_t3  = ipnt(4)
      wing_2d  = ipnt(5)
      wing_1d  = ipnt(6)

#if defined (LAHEY)
      ptrg_4d  = wing_4d
      ptrg_t1  = wing_t1
      ptrg_t2  = wing_t2
      ptrg_t3  = wing_t3
      ptrg_2d  = wing_2d
      ptrg_1d  = wing_1d
#endif

      return
      end subroutine gotmem


      subroutine forkit

#if defined(IRIX64)
#include <ulocks.h>
#endif
      integer fork,getpid
      integer master, n, nowpid, ierror, nowcpu
      character*80 evalue

!-----create mp environment
      call getenv('NUMBER_MLP_PROCESSES',evalue)
      read(evalue,*) numpro
      call getenv('NUMBER_CPUS_PER_MLP_PROCESS',evalue)
      read(evalue,*) numcpu

!-----get master pid
      master = getpid()
      nowpro = 1

!-----print fork message
!!!      write(*,510) numpro

#if defined(IRIX64)
!-----destroy mp environment
      call mp_destroy
#endif

!-----spawn the processes - manual forks
      do n=2,numpro
         nowpid = getpid()
         if(nowpid == master) then 
                              ierror=fork()
                              endif
         nowpid = getpid()
         if(nowpid /= master) then
                              nowpro=n
                              go to 200
                              endif
      enddo

!-----write note
  200 if(nowpro == 1) nowpid = master
!!!     write(*,500) nowpro,nowpid

      call omp_set_num_threads(numcpu)

#if  defined( IRIX64 ) && defined(PIN_CPUS)

!$omp parallel do private(n,nowcpu)
      do n=1,numcpu
         nowcpu = n+(nowpro-1)*numcpu-1
         call mp_assign_to_cpu(nowcpu)
      enddo
#endif

!******************************
!*    I/O formats             *
!******************************
  500 format('FORKIT: Current process:',i3,'    PID:',i10)
  510 format('FORKIT: Total active processes spawned:',i3)

      gid = nowpro-1
      gsize = numpro

      return
      end subroutine forkit
#endif

      subroutine y_decomp(jm, km, jfirst, jlast, kfirst, klast, myid)

      implicit none
      integer jm     ! Dimensions
      integer km     ! Levels
      integer myid

! OUTPUT PARAMETERS:
      integer jfirst, jlast, kfirst, klast

! Local
      integer p, p1, p2, lats, pleft
      integer, allocatable:: ydist(:)

      if (myid == 0) print *, "numpro", numpro, "numcpu", numcpu
      allocate( ydist( numpro ) )

      lats = jm / numpro
      pleft = jm - lats * numpro

      if( lats < 3 ) then
         write(*,*) 'Number of Proc is too large for jm=',jm
         stop
      endif

      do p=1,numpro
         ydist(p) = lats
      enddo

      if ( pleft .ne. 0 ) then
          p1 = (numpro+1) / 2 
          p2 = p1 + 1
        do while ( pleft .ne. 0 )
           if( p1 .eq. 1 ) p1 = numpro
               ydist(p1) = ydist(p1) + 1
               pleft = pleft - 1
               if ( pleft .ne. 0 ) then
                    ydist(p2) = ydist(p2) + 1
                    pleft = pleft - 1
               endif
               p2 = p2 + 1
               p1 = p1 - 1
        enddo
      endif

! Safety check:
      lats = 0
      do p = 1, numpro
         lats = lats + ydist(p)
      enddo

      if ( lats .ne. jm ) then
         print *, "Decomp: big trouble sum(ydist) = ", lats, "!=", jm
      endif
 
      jfirst = 1
      jlast  = ydist(1)
      yfirst(1) = jfirst
      ylast(1) = jlast
      kfirst = 1
      klast = km
      zfirst(1) = kfirst
      zlast(1) = klast

      do p = 1,numpro-1
         yfirst(p+1) = ylast(p) + 1
         ylast(p+1) = ylast(p) + ydist(p+1) 
         if( p == myid ) then
            jfirst = yfirst(p+1)
            jlast  = ylast (p+1)
         endif
         zfirst(p+1) = kfirst
         zlast(p+1) = klast
      enddo

      deallocate (ydist)
      end subroutine y_decomp


!-----
      subroutine set_decomp(nprocs, jm, km, ydist, zdist)
      integer nprocs
      integer jm, km
      integer ydist(nprocs)
      integer zdist(nprocs)   ! Currently not used

!
! Set the decomposition if it is defined external to mod_comm
!
      integer lats, p

! Safety check:
      lats = 0
      do p = 1, nprocs
         lats = lats + ydist(p)
      enddo

      if ( lats .ne. jm ) then
         print *, "Decomp: big trouble sum(ydist) = ", lats, "!=", jm
      endif

      yfirst(1) = 1
      ylast(1) = ydist(1)
      zfirst(1) = 1
      zlast(1) = km

      do p = 1,nprocs-1
         yfirst(p+1) = ylast(p) + 1
         ylast(p+1) = ylast(p) + ydist(p+1)
         zfirst(p+1) = 1
         zlast(p+1) = km
      enddo

      end subroutine set_decomp

!-----
      subroutine mp_send4d_ns(im, jm, jfirst, jlast, kfirst, klast, &
                              nq, ng_s, ng_n, q)
!-----
      implicit none
      integer im, jm
      integer jfirst, jlast
      integer kfirst, klast
      integer nq
      integer ng_s      ! southern zones to ghost 
      integer ng_n      ! noruthern zones to ghost 
      real q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
! Local:
      integer iq        ! Counter
      integer i, j, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif
! Send to south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        src = gid - 1
        recv_tag = src
        qsize = im*ng_s*(klast-kfirst+1)*nq
        nrecv = nrecv + 1
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        dest = gid - 1
        qsize = im*ng_n*(klast-kfirst+1)*nq
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack4d(q, 1, im, jfirst-ng_s, jlast+ng_n, &
                             kfirst, klast, 1, nq, &
                             1, im, jfirst, jfirst+ng_n-1, &
                             kfirst, klast, 1 , nq, &
                             buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(MPI2)
        src = gid + 1
        recv_tag = src
        qsize = im*ng_n*(klast-kfirst+1)*nq
        nrecv = nrecv + 1
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        dest = gid + 1
        qsize = im*ng_s*(klast-kfirst+1)*nq
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack4d(q, 1, im, jfirst-ng_s, jlast+ng_n, &
                             kfirst, klast, 1, nq, &
                             1, im, jlast-ng_s+1, jlast, &
                             kfirst, klast, 1, nq, &
                             buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif
#else
#include "mlp_ptr.h"
      do iq=1,nq
!$omp parallel do private(i,j,k)
      do k=kfirst,klast
         if ( jfirst > 1 ) then
! Send to south
            do j=jfirst,jfirst+ng_n-1
               do i=1,im
                  g_4d(i,j,k,iq) = q(i,j,k,iq)
               enddo
            enddo
         endif
         if ( jlast < jm ) then
! Send to north
            do j=jlast-ng_s+1,jlast
               do i=1,im
                  g_4d(i,j,k,iq) = q(i,j,k,iq)
               enddo
            enddo
         endif
      enddo
      enddo
#endif
      end subroutine mp_send4d_ns
 
!-----
      subroutine mp_recv4d_ns(im, jm, jfirst, jlast, kfirst, klast, &
                              nq, ng_s, ng_n, q)
!-----
      implicit none
      integer im, jm
      integer jfirst, jlast
      integer kfirst, klast
      integer nq
      integer ng_s      ! southern zones to ghost 
      integer ng_n      ! noruthern zones to ghost 
      real q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
! Local:
      integer iq        ! Counter
      integer i, j, k

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
! Recv from south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack4d(q, 1, im, jfirst-ng_s, jlast+ng_n, &
                               kfirst, klast, 1, nq, &
                               1, im, jfirst-ng_s, jfirst-1, &
                               kfirst, klast, 1, nq, &
                               buff_r(tdisp+1))
      endif
! Recv from north
      if ( jlast < jm ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack4d(q, 1, im, jfirst-ng_s, jlast+ng_n, &
                               kfirst, klast, 1, nq, &
                               1, im, jlast+1, jlast+ng_n, &
                               kfirst, klast, 1, nq, &
                               buff_r(tdisp+1))
      endif
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#else
#include "mlp_ptr.h"
      do iq=1,nq

!$omp parallel do private(i,j,k)
      do k=kfirst,klast
         if ( jfirst > 1 ) then
! Recv from south
            do j=jfirst-ng_s,jfirst-1
               do i=1,im
                  q(i,j,k,iq) = g_4d(i,j,k,iq)
               enddo
            enddo
         endif

         if ( jlast < jm ) then
! Recv from north
            do j=jlast+1,jlast+ng_n
               do i=1,im
                  q(i,j,k,iq) = g_4d(i,j,k,iq)
               enddo
            enddo
         endif
      enddo
      enddo
#endif
      end subroutine mp_recv4d_ns

!-----
      subroutine mp_send3d_ns(im, jm, jfirst, jlast, kfirst, klast, & 
                              ng_s, ng_n, q, iq)
!-----
      implicit none
      integer im, jm
      integer jfirst, jlast
      integer kfirst, klast
      integer ng_s      ! southern zones to ghost 
      integer ng_n      ! noruthern zones to ghost 
      real q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)
      integer iq
! Local:
      integer i,j,k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif
! Send to south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        src = gid - 1
        recv_tag = src
        qsize = im*ng_s*(klast-kfirst+1)
        nrecv = nrecv + 1
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif

        dest = gid - 1
        qsize = im*ng_n*(klast-kfirst+1)
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf 
        call BufferPack3d(q, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                             1, im, jfirst, jfirst+ng_n-1, kfirst, klast, &
                             buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(MPI2)
        src = gid + 1
        recv_tag = src
        qsize = im*ng_n*(klast-kfirst+1)
        nrecv = nrecv + 1
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        dest = gid + 1
        qsize = im*ng_s*(klast-kfirst+1)
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf 
        call BufferPack3d(q, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                             1, im, jlast-ng_s+1, jlast, kfirst, klast, &
                             buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif
#else
#include "mlp_ptr.h"
!$omp parallel do private(i,j,k)
      do k=kfirst,klast
         if ( jfirst > 1 ) then
! Send to south
            do j=jfirst,jfirst+ng_n-1
               do i=1,im
                  g_4d(i,j,k,iq) = q(i,j,k)
               enddo
            enddo
         endif
         if ( jlast < jm ) then
! Send to north
            do j=jlast-ng_s+1,jlast
               do i=1,im
                  g_4d(i,j,k,iq) = q(i,j,k)
               enddo
            enddo
         endif
      enddo
#endif
      end subroutine mp_send3d_ns
 
!-----
      subroutine mp_recv3d_ns(im, jm, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q, iq)
!-----
      implicit none
      integer im, jm
      integer jfirst, jlast
      integer kfirst, klast
      integer ng_s      ! southern zones to ghost 
      integer ng_n      ! noruthern zones to ghost 
      integer iq        ! Counter
      real q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)

! Local:
      integer i,j,k
      integer src
      integer recv_tag

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
! Recv from south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3d(q, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                               1, im, jfirst-ng_s, jfirst-1,   kfirst, klast, &
                               buff_r(tdisp+1))
      endif
! Recv from north
      if ( jlast < jm ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3d(q, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                               1, im, jlast+1,     jlast+ng_n, kfirst, klast, &
                               buff_r(tdisp+1))
      endif
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#else
#include "mlp_ptr.h"
!$omp parallel do private(i,j,k)
      do k=kfirst,klast
         if ( jfirst > 1 ) then
! Recv from south
            do j=jfirst-ng_s,jfirst-1
               do i=1,im
                  q(i,j,k) = g_4d(i,j,k,iq)
               enddo
            enddo
         endif

         if ( jlast < jm ) then
! Recv from north
            do j=jlast+1,jlast+ng_n
               do i=1,im
                  q(i,j,k) = g_4d(i,j,k,iq)
               enddo
            enddo
         endif
      enddo
#endif
      end subroutine mp_recv3d_ns

!-----
      subroutine mp_send2_n(im, jm, jfirst, jlast, kfirst, klast, &
                            ng_s, ng_n, q1, q2)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast  ! careful: might be klast+1 outside
      integer, intent(in):: ng_s
      integer, intent(in):: ng_n
      real, intent(in):: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      real, intent(in):: q2(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      integer i, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
      integer displ
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#else
! Start receive from south
      if ( jfirst > 1 ) then
         src = gid - 1
         recv_tag = src
         qsize = im*2*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
#endif

! Send data to north
      if ( jlast < jm ) then
#if !defined(USE_MLP)
        dest = gid + 1
        qsize = im*(klast-kfirst+1)*2
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3d(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                              1, im, jlast,       jlast,      kfirst, klast, &
                              buff_s(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferPack3d(q2, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                              1, im, jlast,       jlast,      kfirst, klast, &
                              buff_s(displ+tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif

#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_t1(i,jlast,k,1) = q1(i,jlast,k)
               g_t1(i,jlast,k,2) = q2(i,jlast,k)
            enddo
         enddo
#endif
      endif
      end subroutine mp_send2_n

!-----
      subroutine mp_send2_s(im, jm, jfirst, jlast, kfirst, klast, &
                            ng_s, ng_n, q1, q2)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast   !careful: klast might be klast+1
      integer, intent(in):: ng_s
      integer, intent(in):: ng_n
      real, intent(in):: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      real, intent(in):: q2(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 

! Local:
      integer i, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
      integer displ
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#else
! Start recv from north
      if ( jlast < jm ) then
         src = gid + 1
         recv_tag = src
         qsize = im*2*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
#endif

! Send data to south
      if ( jfirst > 1 ) then
#if !defined(USE_MLP)
        dest = gid - 1
        qsize = im*(klast-kfirst+1)*2
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3d(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                              1, im, jfirst,      jfirst,     kfirst, klast, &
                              buff_s(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferPack3d(q2, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                              1, im, jfirst,      jfirst,     kfirst, klast, &
                              buff_s(displ+tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend  = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif

#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_t1(i,jfirst,k,1) = q1(i,jfirst,k)
               g_t1(i,jfirst,k,2) = q2(i,jfirst,k)
            enddo
         enddo
#endif
      endif
      end subroutine mp_send2_s

!-----
      subroutine mp_recv2_s(im, jm, jfirst, jlast, kfirst, klast, &
                            ng_s, ng_n, q1, q2)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast     ! careful: klast might be klast+1
      integer, intent(in):: ng_s
      integer, intent(in):: ng_n
      real, intent(inout):: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      real, intent(inout):: q2(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 

! Local:
      integer i,j, k, n
      integer displ

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv data from south
      if ( jfirst > 1 ) then
        j = jfirst - 1
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3d(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                                1, im, j,           j,          kfirst, klast, &
                                buff_r(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferUnPack3d(q2, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                                1, im, j,           j,          kfirst, klast, &
                                buff_r(displ+tdisp+1))
#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               q1(i,j,k) = g_t1(i,j,k,1) 
               q2(i,j,k) = g_t1(i,j,k,2) 
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_recv2_s

!-----
      subroutine mp_recv2_n(im, jm, jfirst, jlast, kfirst, klast, &
                            ng_s, ng_n, q1, q2)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast     !careful: klast might be klast+1
      integer, intent(in):: ng_s
      integer, intent(in):: ng_n
      real, intent(inout):: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      real, intent(inout):: q2(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 

! Local:
      integer i,j, k, n
      integer displ

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv data from north
      if ( jlast < jm ) then
        j = jlast + 1
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3d(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                                1, im, j,           j,          kfirst, klast, &
                                buff_r(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferUnPack3d(q2, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                                1, im, j,           j,          kfirst, klast, &
                                buff_r(displ+tdisp+1))
#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               q1(i,j,k) = g_t1(i,j,k,1) 
               q2(i,j,k) = g_t1(i,j,k,2) 
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_recv2_n

!-----
      subroutine mp_send_s(im, jm, jfirst, jlast, kfirst, klast, &
                           nd_s, nd_n, q)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: nd_s, nd_n
      real, intent(in):: q(im,jfirst-nd_s:jlast+nd_n,kfirst:klast) 

! Local:
      integer i, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#else
! Start receive from north
      if ( jlast < jm ) then
         src = gid + 1
         recv_tag = src
         qsize = im*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
#endif

! Send data to south
      if ( jfirst > 1 ) then
#if !defined(USE_MLP)
        dest = gid - 1
        qsize = im*(klast-kfirst+1)
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3d(q, 1, im, jfirst-nd_s, jlast+nd_n, kfirst, klast, &
                             1, im, jfirst,      jfirst,     kfirst, klast, &
                             buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
         send_tag = gid
         nsend = nsend + 1
         call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                        send_tag, commglobal, sqest(nsend), ierror)
#endif

#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_4d(i,jfirst,k,1) = q(i,jfirst,k)
            enddo
         enddo
#endif
      endif
      end subroutine mp_send_s

!-----
      subroutine mp_recv_n(im, jm, jfirst, jlast, kfirst, klast, nd_s, nd_n, q)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: nd_s, nd_n
      real, intent(inout):: q(im,jfirst-nd_s:jlast+nd_n,kfirst:klast) 

! Local:
      integer i, j, k, n

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv data from north
      if ( jlast < jm ) then
        j = jlast + 1
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3d(q, 1, im, jfirst-nd_s, jlast+nd_n, kfirst, klast, &
                               1, im, j,           j,          kfirst, klast, &
                               buff_r(tdisp+1))
#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               q(i,j,k) = g_4d(i,j,k,1) 
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_recv_n

!-----
      subroutine mp_send_n(im, jm, jfirst, jlast, kfirst, klast, &
                           ng_s, ng_n, q1)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast  ! careful: might be klast+1 outside
      integer, intent(in):: ng_s
      integer, intent(in):: ng_n
      real, intent(in):: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      integer i, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
      integer displ
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#else
! Start receive from south
      if ( jfirst > 1 ) then
         src = gid - 1
         recv_tag = src
         qsize = im*1*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
#endif

! Send data to north
      if ( jlast < jm ) then
#if !defined(USE_MLP)
        dest = gid + 1
        qsize = im*(klast-kfirst+1)*2
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3d(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                              1, im, jlast,       jlast,      kfirst, klast, &
                              buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif

#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_t1(i,jlast,k,1) = q1(i,jlast,k)
            enddo
         enddo
#endif
      endif
      end subroutine mp_send_n

!-----
      subroutine mp_recv_s(im, jm, jfirst, jlast, kfirst, klast, &
                            ng_s, ng_n, q1)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast     ! careful: klast might be klast+1
      integer, intent(in):: ng_s
      integer, intent(in):: ng_n
      real, intent(inout):: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 

! Local:
      integer i,j, k, n
      integer displ

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv data from south
      if ( jfirst > 1 ) then
        j = jfirst - 1
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3d(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                                1, im, j,           j,          kfirst, klast, &
                                buff_r(tdisp+1))
#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               q1(i,j,k) = g_t1(i,j,k,1) 
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_recv_s


!-----
      subroutine mp_send2_ns(im, jm, jfirst, jlast, kfirst, klast, nd, q1, q2)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
      real, intent(in):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real, intent(in):: q2(im,jfirst-nd:jlast+nd,kfirst:klast) 

! Local:
      integer i, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
      integer displ
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif

! Send to south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
! Start recv to north
        src = gid - 1
        recv_tag = src
        qsize = im*2*(klast-kfirst+1)
        nrecv = nrecv + 1
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        dest = gid - 1
        qsize = im*(klast-kfirst+1)*2
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3d(q1, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                              1, im, jfirst,    jfirst,   kfirst, klast, &
                              buff_s(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferPack3d(q2, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                              1, im, jfirst,    jfirst,   kfirst, klast, &
                              buff_s(displ+tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif

! Send to north
      if ( jlast < jm ) then
#if !defined(MPI2)
! Start recv to south
        src = gid + 1
        recv_tag = src
        qsize = im*2*(klast-kfirst+1)
        nrecv = nrecv + 1
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        dest = gid + 1
        qsize = im*(klast-kfirst+1)*2
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3d(q1, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                              1, im, jlast,     jlast,    kfirst, klast, &
                              buff_s(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferPack3d(q2, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                              1, im, jlast,     jlast,    kfirst, klast, &
                              buff_s(displ+tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif

#else
#include "mlp_ptr.h"
!$omp parallel do private(i, k)
      do k=kfirst,klast
! Send to south
      if ( jfirst > 1 ) then
           do i=1,im
              g_t1(i,jfirst,k,1) = q1(i,jfirst,k)
              g_t1(i,jfirst,k,2) = q2(i,jfirst,k)
           enddo
      endif

! Send to north
      if ( jlast < jm ) then
           do i=1,im
              g_t1(i,jlast,k,1) = q1(i,jlast,k)
              g_t1(i,jlast,k,2) = q2(i,jlast,k)
           enddo
      endif
      enddo
#endif
      end subroutine mp_send2_ns

!-----
      subroutine mp_recv2_ns(im, jm, jfirst, jlast, kfirst, klast, nd, q1, q2)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
      real, intent(inout):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real, intent(inout):: q2(im,jfirst-nd:jlast+nd,kfirst:klast) 

! Local:
      integer i,j, k, n
      integer displ

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
! Recv from south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        j = jfirst - 1
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3d(q1, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                                1, im, j,         j,        kfirst, klast, &
                                buff_r(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferUnPack3d(q2, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                                1, im, j,         j,        kfirst, klast, &
                                buff_r(displ+tdisp+1))
      endif
! Recv from north
      if ( jlast < jm ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        j = jlast + 1
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3d(q1, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                                1, im, j,         j,        kfirst, klast, &
                                buff_r(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferUnPack3d(q2, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                                1, im, j,         j,        kfirst, klast, &
                                buff_r(displ+tdisp+1))
      endif
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#else
#include "mlp_ptr.h"
!$omp parallel do private(i, j, k)
      do k=kfirst,klast
! Recv data from south
      if ( jfirst > 1 ) then
            j = jfirst - 1
            do i=1,im
               q1(i,j,k) = g_t1(i,j,k,1) 
               q2(i,j,k) = g_t1(i,j,k,2) 
            enddo
      endif

! Recv data from north
      if ( jlast < jm ) then
            j = jlast + 1
            do i=1,im
               q1(i,j,k) = g_t1(i,j,k,1) 
               q2(i,j,k) = g_t1(i,j,k,2) 
            enddo
      endif
      enddo
#endif
      end subroutine mp_recv2_ns

      subroutine mp_send_pe(im, jm, jfirst, jlast, kfirst, klast, p)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast, kfirst, klast
      real, intent(in):: p(im,kfirst:klast,jfirst:jlast) 
      integer i, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer n, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif
#endif

#if !defined(USE_MLP) && !defined(MPI2)
! Start recv from south
      if ( jfirst > 1 ) then
         src = gid - 1
         recv_tag = src
         qsize = im*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
! Send data to North
      if ( jlast < jm ) then
#if !defined(USE_MLP)
        dest = gid + 1
        qsize = im*(klast-kfirst+1)
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack2d(p(1,1,jlast), 1, im, kfirst, klast, &
                                        1, im, kfirst, klast, &
                                        buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(n,tmpsize,mysize,mydisp)
        do n=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(n-1)),0))
          mydisp = tdisp + (n-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif

#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_t2(i,k,jlast) = p(i,k,jlast)
            enddo
         enddo
#endif
      endif
      end subroutine mp_send_pe

!-----
      subroutine mp_recv_pe(im, jm, jfirst, jlast, kfirst, klast, pesouth)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast, kfirst, klast
      real, intent(inout):: pesouth(im,kfirst:klast) 
      integer i, j, k, n
#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv from south
      if ( jfirst > 1 ) then
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack2d(pesouth, 1, im, kfirst, klast, &
                                     1, im, kfirst, klast, &
                                     buff_r(tdisp+1))
#else
!$omp parallel do private(i, j, k)
         do k=kfirst,klast
            j = jfirst - 1
            do i=1,im
               pesouth(i,k) = g_t2(i,k,j)
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_recv_pe

!-----
      subroutine mp_send_ua(im, jm, jfirst, jlast, kfirst, klast, p)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast, kfirst, klast
      real, intent(in):: p(im,jfirst:jlast,kfirst:klast) 
      integer i, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer n, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif
#endif

#if !defined(USE_MLP) && !defined(MPI2)
! Start recv from north
      if ( jfirst > 1 ) then
         src = gid - 1
         recv_tag = src
         qsize = im*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
! Send data to North
      if ( jlast < jm ) then
#if !defined(USE_MLP)
        dest = gid + 1
        qsize = im*(klast-kfirst+1)
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3d(p, 1, im, jfirst, jlast, kfirst, klast, &
                             1, im, jlast,  jlast, kfirst, klast, &
                             buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(n,tmpsize,mysize,mydisp)
        do n=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(n-1)),0))
          mydisp = tdisp + (n-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_4d(i,jlast,k,1) = p(i,jlast,k)
            enddo
         enddo
#endif
      endif
      end subroutine mp_send_ua

!-----
      subroutine mp_recv_ua(im, jm, jfirst, jlast, kfirst, klast, uasouth)
!-----
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast, kfirst, klast
      real, intent(inout):: uasouth(im, kfirst:klast) 
      integer i, j, k, n
#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv from south
      if ( jfirst > 1 ) then
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack2d(uasouth, 1, im, kfirst, klast, &
                                     1, im, kfirst, klast, &
                                     buff_r(tdisp+1))
#else
!$omp parallel do private(i, j, k)
         do k=kfirst,klast
            j = jfirst - 1
            do i=1,im
               uasouth(i,k) = g_4d(i,j,k,1)
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_recv_ua

!-----
      subroutine mp_reduce_max(km, cymax)
!-----
      implicit none
      integer k, km, n
      real maxin(km)
      real cymax(km)

#if !defined(USE_MLP)
!$omp parallel do private(k)
      do k=1,km
        maxin(k) = cymax(k)
      enddo
      call mpi_allreduce( maxin, cymax, km, mp_precision, MPI_MAX, &
                          commglobal, ierror )
#else
#include "mlp_ptr.h"
      do k=1,km
         g_t3(k,nowpro) = cymax(k)
      enddo
      call mlp_barrier(gid, gsize)

      do n=1,numpro
         do k=1,km
            cymax(k) = max(g_t3(k,n), cymax(k))
         enddo
      enddo
      call mlp_barrier(gid, gsize) !may not be necessay, test, BW
#endif
      end subroutine mp_reduce_max

!-----
      subroutine mp_minmax(qmin, qmax)
!-----
      implicit none
      real, intent(inout):: qmin, qmax
      real minin, maxin
      integer n

#if !defined(USE_MLP)
      maxin = qmax
      call mpi_allreduce(maxin, qmax, 1, mp_precision, MPI_MAX, &
                         commglobal, ierror)
      minin = qmin
      call mpi_allreduce(minin, qmin, 1, mp_precision, MPI_MIN, &
                         commglobal, ierror)
#else
#include "mlp_ptr.h"
      g_t3(1,nowpro) = qmin
      g_t3(2,nowpro) = qmax
      call mlp_barrier(gid, gsize)

      do n=1,numpro
         qmin = min(g_t3(1,n), qmin)
         qmax = max(g_t3(2,n), qmax)
      enddo
      call mlp_barrier(gid, gsize)
#endif
      end subroutine mp_minmax

!-----
      subroutine mp_sum1d(jm, jfirst, jlast, qin, sum0)
!-----
      implicit none
      integer jm
      integer jfirst, jlast
      real  qin(jfirst:jlast)
! Output:
      real  sum0

#ifdef FASTOPT_FASTER_BUT_RUNDOFF
      sum0 = 0.0
      do j = jfirst, jlast
        sum0 = sum0 + qin(j)
      enddo

      call mp_reduce_sum( sum0 )
#else
! Local:
      integer j, n
      real qout(jm)
#if !defined(USE_MLP)
      call mp_allgather1d(jm, jfirst, jlast, qin, qout)
      sum0 = 0.
      do j=1,jm
        sum0 = sum0 + qout(j)
      enddo
#else
#include "mlp_ptr.h"
! Gather all subdomain vector from all PEs to a global array
      do j=jfirst,jlast
         g_1d(j) = qin(j)
      enddo
      call mlp_barrier(gid, gsize)

! Compute the sum if "Master"
      if ( gid == 0 ) then
            sum0 = 0.
         do j=1,jm
            sum0 = sum0 + g_1d(j)
         enddo
      endif
      call mp_bcst_real(sum0)
#endif
#endif
      end subroutine mp_sum1d

!-----
      subroutine mp_bcst_real(val)
!-----
! Send real "val" from Process=id to All other Processes
      real val

#if !defined(USE_MLP)
      call mpi_bcast(val, 1, mp_precision, 0, commglobal, ierror)
#else
#include "mlp_ptr.h"
      if ( gid == 0 ) then
          g_1d(1) = val
      endif

      call mlp_barrier(gid, gsize)

      if ( gid /= 0 ) then
          val = g_1d(1)
      endif
      call mlp_barrier(gid, gsize)   !may not be necessary, BW
#endif
      end subroutine mp_bcst_real


!-----
      subroutine mp_bcst_int(intv)
!-----
! Send integer "intv" from Process=id to All other Processes
      integer intv

#if !defined(USE_MLP)
      call mpi_bcast(intv, 1, MPI_INTEGER, 0, commglobal, ierror)
#else
#include "mlp_ptr.h"

      if ( gid == 0 ) then
          g_1d(1) = intv
      endif

      call mlp_barrier(gid, gsize)

      if ( gid /= 0 ) then
          intv = nint(g_1d(1))
      endif

      call mlp_barrier(gid, gsize)   !may not be necessary, BW
#endif
      end subroutine mp_bcst_int


!-----
      subroutine mp_bcst_i2d(im, jm, jfirst, jlast, qin, id)
!-----
! Send 2D array qin from Process=id to All other Processes
      integer im, jm
      integer id        ! source ID
      integer jfirst, jlast
      integer qin(im,jm)
      integer i, j, n
      integer j1, j2
      integer qsize_s
      integer qsize_r
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP)
      integer rqst(numpro), rq_stats(numpro*MPI_STATUS_SIZE)
#if defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#endif

#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff3d_iwin, ierror)
#endif
      if (gid == id) then
         do n=1,numpro
           qsize_s = im*jm
           dest = n-1
           tdisp = 0
#if defined(MPI2)
           call BufferPack2d_i(qin, 1, im, 1, jm, 1, im, 1, jm, &
                               buff3d_i(tdisp+1))
!$omp parallel do private(p,tmpsize,mysize,mydisp)
           do p=1,pkgs_per_pro
             tmpsize = ceiling(real(qsize_s)/real(pkgs_per_pro))
             mysize = MIN(tmpsize, MAX(qsize_s-(tmpsize*(p-1)),0))
             mydisp = tdisp + (p-1)*tmpsize
             call MPI_PUT(buff3d_i(mydisp+1), mysize, MPI_INTEGER, dest, &
                          mydisp, mysize, MPI_INTEGER, buff3d_iwin, ierror)
           enddo
#else
           send_tag = gid
           nsend = nsend + 1
           call mpi_isend(qin, qsize_s, MPI_INTEGER, dest, send_tag, &
                          commglobal, rqst(nsend), ierror)
#endif
         enddo
      endif
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, buff3d_iwin, &
                         ierror)
#else
      qsize_r = im*jm
      src = id
      recv_tag = src
      call mpi_recv(buff3d_i, qsize_r, MPI_INTEGER, src, recv_tag, &
                    commglobal, Status, ierror)
#endif
      tdisp = (jfirst-1)*im
      call BufferUnPack2d_i(qin, 1, im, 1, jm, 1, im, jfirst, jlast, &
                            buff3d_i(tdisp+1))

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, rqst, rq_stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
      if ( gid == id ) then
!$omp parallel do private(i, j)
          do j=1,jm
             do i=1,im
                g_2d(i,j) = qin(i,j)
             enddo
          enddo
      endif
      call mlp_barrier(gid, gsize)

      if ( gid /= id ) then
!$omp parallel do private(i, j)
          do j=jfirst,jlast
             do i=1,im
                qin(i,j) = nint(g_2d(i,j))
             enddo
          enddo
      endif
      call mlp_barrier(gid, gsize)
#endif
      end subroutine mp_bcst_i2d


!-----
      subroutine mp_bcst_r2d(im, jm, jfirst, jlast, qin, id)
!-----
! Send 2D array qin from Process=id to All other Processes
      integer im, jm
      integer id        ! source ID
      integer jfirst, jlast
      real qin(im,jm)
      integer i, j, n 
      integer j1, j2
      integer qsize_s
      integer qsize_r
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP)
      integer rqst(numpro), rq_stats(numpro*MPI_STATUS_SIZE)
#if defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#endif

#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4dwin, ierror)
#endif
      if (gid == id) then
         do n=1,numpro
           qsize_s = im*jm
           dest = n-1
           tdisp = 0
#if defined(MPI2)
           call BufferPack2d(qin, 1, im, 1, jm, 1, im, 1, jm, buff4d(tdisp+1))
!$omp parallel do private(p,tmpsize,mysize,mydisp)
           do p=1,pkgs_per_pro
             tmpsize = ceiling(real(qsize_s)/real(pkgs_per_pro))
             mysize = MIN(tmpsize, MAX(qsize_s-(tmpsize*(p-1)),0))
             mydisp = tdisp + (p-1)*tmpsize
             call MPI_PUT(buff4d(mydisp+1), mysize, mp_precision ,dest, &
                          mydisp, mysize, mp_precision, buff4dwin, ierror)
           enddo
#else
           send_tag = gid
           nsend = nsend + 1
           call mpi_isend(qin, qsize_s, mp_precision, dest, &
                          send_tag, commglobal, rqst(nsend), ierror)
#endif
         enddo
      endif
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4dwin, ierror)
#else
      qsize_r = im*jm
      src = id
      recv_tag = src
      call mpi_recv(buff4d, qsize_r, mp_precision, src, recv_tag, &
                    commglobal, Status, ierror)
#endif
      tdisp = (jfirst-1)*im
      call BufferUnPack2d(qin, 1, im, 1, jm, 1, im, jfirst, jlast, &
                          buff4d(tdisp+1))

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, rqst, rq_stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
      if ( gid == id ) then
!$omp parallel do private(i, j)
          do j=1,jm
             do i=1,im
                g_2d(i,j) = qin(i,j)
             enddo
          enddo
      endif
      call mlp_barrier(gid, gsize)

      if ( gid /= id ) then
!$omp parallel do private(i, j)
          do j=jfirst,jlast
             do i=1,im
                qin(i,j) = g_2d(i,j)
             enddo
          enddo
      endif
      call mlp_barrier(gid, gsize)
#endif
      end subroutine mp_bcst_r2d


!-----
      subroutine mp_gath_r2d(im, jm, jfirst, jlast, qin, id)
!-----
      integer im, jm
      integer id        ! source ID
      integer jfirst, jlast
      real qin(im,jm)

! Local:
      integer i, j, k, iq
      integer j1, j2
      integer n
      integer qsize
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
  
#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4dwin, ierror)
#endif
      dest = id
      tdisp = (jfirst-1)*im
      call BufferPack2d(qin, 1, im, 1, jm, 1, im, jfirst, jlast, &
                        buff4d(tdisp+1))
      qsize = (jlast-jfirst+1)*im
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
      do p=1,pkgs_per_pro
        tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = tdisp + (p-1)*tmpsize
        call MPI_PUT(buff4d(mydisp+1), mysize, mp_precision, dest, &
                     mydisp, mysize, mp_precision, buff4dwin, ierror)
      enddo
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4dwin, ierror)
#else
      send_tag = gid
      nsend = nsend + 1
      call mpi_isend(buff4d(tdisp+1), qsize, mp_precision, dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
      if (gid==id) then
        do n=0,numpro-1
          j1 = yfirst(n+1)
          j2 = ylast(n+1)
          qsize = im*(j2-j1+1)
          src = n
          recv_tag = src
          tdisp = (j1-1)*im
          call mpi_recv(buff4d(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, Status, ierror)
        enddo
      endif
#endif
      if (gid == id) then
        call BufferUnPack2d(qin, 1, im, 1, jm, 1, im, 1, jm, buff4d)
      endif

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
!$omp parallel do private(i, j)
          do j=jfirst,jlast
             do i=1,im
                g_2d(i,j) = qin(i,j)
             enddo
          enddo
      call mlp_barrier(gid, gsize)

      if ( gid == id ) then
!$omp parallel do private(i, j)
          do j=1,jm
             do i=1,im
                qin(i,j) = g_2d(i,j) 
             enddo
          enddo
      endif
      call mlp_barrier(gid, gsize)
#endif
      end subroutine mp_gath_r2d

!-----
      subroutine mp_gath_3d(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, id)
!-----

      integer idim, jdim, kdim
      integer id
      integer i1, i2, j1, j2, k1, k2
      real    q(idim,jdim,kdim)
      integer i, j, k, n, p
      integer nsend, rqst1, rqst2
      integer qsize
      integer dest, src
      integer recv_tag, send_tag
      integer ir1, ir2, jr1, jr2, kr1, kr2
      integer dims_snd(6)
      integer dims_rcv(6)
      real buff_snd(i1:i2,j1:j2,k1:k2)
      real buff_rcv(idim*jdim*kdim)

#if !defined(USE_MLP)
      nsend = 0
      if (gid /= id) then
        nsend = nsend + 1
        dims_snd(1) = i1
        dims_snd(2) = i2
        dims_snd(3) = j1
        dims_snd(4) = j2
        dims_snd(5) = k1
        dims_snd(6) = k2
        qsize = 6
        dest = id
        send_tag = gid+1
        call mpi_isend(dims_snd, qsize, MPI_INTEGER, dest, send_tag, &
                       commglobal, rqst1, ierror)

        nsend = nsend + 1
        call BufferPack3d(q, 1,  idim, 1,  jdim, 1,  kdim, &
                             i1, i2,   j1, j2,   k1, k2,   &
                             buff_snd)
        qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        dest = id
        send_tag = gid
        call mpi_isend(buff_snd, qsize, mp_precision, dest, send_tag, &
                       commglobal, rqst2, ierror)

      else

        do n=1,numpro-1
          qsize = 6
          src = n
          recv_tag = src + 1
          call mpi_recv(dims_rcv, qsize, MPI_INTEGER, src, recv_tag, &
                        commglobal, Status, ierror)
          ir1 = dims_rcv(1)
          ir2 = dims_rcv(2)
          jr1 = dims_rcv(3)
          jr2 = dims_rcv(4)
          kr1 = dims_rcv(5)
          kr2 = dims_rcv(6)
          qsize = (ir2-ir1+1)*(jr2-jr1+1)*(kr2-kr1+1)
          src = n
          recv_tag = src
          call mpi_recv(buff_rcv, qsize, mp_precision, src, recv_tag, &
                        commglobal, Status, ierror)
          call BufferUnPack3d(q, 1,   idim, 1,   jdim, 1,   kdim, &
                                 ir1, ir2,  jr1, jr2,  kr1, kr2,  &
                                 buff_rcv)
        enddo
      endif

      if (nsend /= 0) then
         call mpi_wait(rqst1, Status, ierror)
         call mpi_wait(rqst2, Status, ierror)
         nsend = 0
      endif

#else
#include "mlp_ptr.h"
      call mp_cp_3d(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, wing_t1, 1)
      call mp_barrier

      if ( gid == id ) then
        call mp_cp_3d(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, wing_t1, -1)
      endif
      call mp_barrier
#endif
      end subroutine mp_gath_3d

#if defined(USE_MLP)
!-----
      subroutine mp_cp_3d(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, &
                          wing, dir)
!-----

      integer idim, jdim, kdim
      integer i1, i2, j1, j2, k1, k2
      real    q(idim,jdim,kdim)
      integer dir         ! direction: 1 -> local to global
      integer i, j, k

      pointer (wing, b3d)
      real b3d(idim,jdim,kdim)

      if ( dir == 1 ) then

!$omp parallel do private(i, j, k)
      do k=k1,k2
         do j=j1,j2
            do i=i1,i2
               b3d(i,j,k) = q(i,j,k)
            enddo
         enddo
      enddo

      else

!$omp parallel do private(i, j, k)
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               q(i,j,k) = b3d(i,j,k)
            enddo
         enddo
      enddo

      endif

      end subroutine mp_cp_3d

!-----
      subroutine mp_cp_3d_int(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, &
                              wing, dir)
!-----

      integer idim, jdim, kdim
      integer i1, i2, j1, j2, k1, k2
      integer q(idim,jdim,kdim)
      integer dir         ! direction: 1 -> local to global
      integer i, j, k

      pointer (wing, b3d)
      integer b3d(idim,jdim,kdim)

      if ( dir == 1 ) then

!$omp parallel do private(i, j, k)
      do k=k1,k2
         do j=j1,j2
            do i=i1,i2
               b3d(i,j,k) = q(i,j,k)
            enddo
         enddo
      enddo

      else

!$omp parallel do private(i, j, k)
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               q(i,j,k) = b3d(i,j,k)
            enddo
         enddo
      enddo

      endif

      end subroutine mp_cp_3d_int
#endif

!-----
      subroutine mp_scat3d_int(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, id)
!-----

      integer idim, jdim, kdim
      integer id
      integer i1, i2, j1, j2, k1, k2
      integer q(idim,jdim,kdim)
      integer i, j, k, n, p
#if !defined(USE_MLP)
      integer nsend, rqst1(numpro-1), rqst2(numpro-1)
      integer rq_stats((numpro-1) * MPI_STATUS_SIZE)
      integer qsize
      integer dest, src
      integer recv_tag, send_tag
      integer ir1, ir2, jr1, jr2, kr1, kr2
      integer dims_snd(6)
      integer dims_rcv(6)
      integer buff_snd(i1:i2,j1:j2,k1:k2)
      integer buff_rcv(idim*jdim*kdim)
#endif

#if !defined(USE_MLP)
      nsend = 0
      if (gid == id) then
        nsend = nsend + 1
        dims_snd(1) = i1
        dims_snd(2) = i2
        dims_snd(3) = j1
        dims_snd(4) = j2
        dims_snd(5) = k1
        dims_snd(6) = k2
        qsize = 6
        send_tag = gid+1
        do n=1,numpro-1
          dest = n
          call mpi_isend(dims_snd, qsize, MPI_INTEGER, dest, send_tag, &
                         commglobal, rqst1(n), ierror)
        enddo

        call BufferPack3d_i(q, 1,  idim, 1,  jdim, 1,  kdim, &
                               i1, i2,   j1, j2,   k1, k2,   &
                               buff_snd)
        qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        send_tag = id
        do n=1,numpro-1
          nsend = nsend + 1
          dest = n
          call mpi_isend(buff_snd, qsize, MPI_INTEGER, dest, send_tag, &
                         commglobal, rqst2(n), ierror)
        enddo

      else

        qsize = 6
        src = id
        recv_tag = src + 1
        call mpi_recv(dims_rcv, qsize, MPI_INTEGER, src, recv_tag, &
                      commglobal, Status, ierror)
        ir1 = dims_rcv(1)
        ir2 = dims_rcv(2)
        jr1 = dims_rcv(3)
        jr2 = dims_rcv(4)
        kr1 = dims_rcv(5)
        kr2 = dims_rcv(6)
        qsize = (ir2-ir1+1)*(jr2-jr1+1)*(kr2-kr1+1)
        src = id
        recv_tag = src
        call mpi_recv(buff_rcv, qsize, MPI_INTEGER, src, recv_tag, &
                      commglobal, Status, ierror)
        call BufferUnPack3d_i(q, 1,   idim, 1,   jdim, 1,   kdim, &
                                 ir1, ir2,  jr1, jr2,  kr1, kr2,  &
                                 buff_rcv)

      endif

      if (nsend /= 0) then
         call mpi_waitall(numpro-1, rqst1, rq_stats, ierror)
         call mpi_waitall(numpro-1, rqst2, rq_stats, ierror)
         nsend = 0
      endif

#else
#include "mlp_ptr.h"
      if ( gid == id ) then
        call mp_cp_3d_int(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, &
                          wing_t1, 1)
      endif
      call mp_barrier

      if ( gid /= id ) then
        call mp_cp_3d_int(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, &
                          wing_t1, -1)
      endif
      call mp_barrier
#endif
      end subroutine mp_scat3d_int

!-----
      subroutine mp_scat3d(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, id)
!-----

      integer idim, jdim, kdim
      integer id
      integer i1, i2, j1, j2, k1, k2
      real    q(idim,jdim,kdim)
      integer i, j, k, n, p
#if !defined(USE_MLP)
      integer nsend, rqst1(numpro-1), rqst2(numpro-1)
      integer rq_stats((numpro-1) * MPI_STATUS_SIZE)
      integer qsize
      integer dest, src
      integer recv_tag, send_tag
      integer ir1, ir2, jr1, jr2, kr1, kr2
      integer dims_snd(6)
      integer dims_rcv(6)
      real buff_snd(i1:i2,j1:j2,k1:k2)
      real buff_rcv(idim*jdim*kdim)
#endif

#if !defined(USE_MLP)
      nsend = 0
      if (gid == id) then
        nsend = nsend + 1
        dims_snd(1) = i1
        dims_snd(2) = i2
        dims_snd(3) = j1
        dims_snd(4) = j2
        dims_snd(5) = k1
        dims_snd(6) = k2
        qsize = 6
        send_tag = gid+1
        do n=1,numpro-1
          dest = n
          call mpi_isend(dims_snd, qsize, MPI_INTEGER, dest, send_tag, &
                         commglobal, rqst1(n), ierror)
        enddo

        call BufferPack3d(q, 1,  idim, 1,  jdim, 1,  kdim, &
                             i1, i2,   j1, j2,   k1, k2,   &
                             buff_snd)
        qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        send_tag = id
        do n=1,numpro-1
          nsend = nsend + 1
          dest = n
          call mpi_isend(buff_snd, qsize, mp_precision, dest, send_tag, &
                         commglobal, rqst2(n), ierror)
        enddo

      else

        qsize = 6
        src = id
        recv_tag = src + 1
        call mpi_recv(dims_rcv, qsize, MPI_INTEGER, src, recv_tag, &
                      commglobal, Status, ierror)
        ir1 = dims_rcv(1)
        ir2 = dims_rcv(2)
        jr1 = dims_rcv(3)
        jr2 = dims_rcv(4)
        kr1 = dims_rcv(5)
        kr2 = dims_rcv(6)
        qsize = (ir2-ir1+1)*(jr2-jr1+1)*(kr2-kr1+1)
        src = id
        recv_tag = src
        call mpi_recv(buff_rcv, qsize, mp_precision, src, recv_tag, &
                      commglobal, Status, ierror)
        call BufferUnPack3d(q, 1,   idim, 1,   jdim, 1,   kdim, &
                               ir1, ir2,  jr1, jr2,  kr1, kr2,  &
                               buff_rcv)

      endif

      if (nsend /= 0) then
         call mpi_waitall(numpro-1, rqst1, rq_stats, ierror)
         call mpi_waitall(numpro-1, rqst2, rq_stats, ierror)
         nsend = 0
      endif

#else
#include "mlp_ptr.h"
      if ( gid == id ) then
        call mp_cp_3d(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, wing_t1, 1)
      endif
      call mp_barrier

      if ( gid /= id ) then
        call mp_cp_3d(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, wing_t1, -1)
      endif
      call mp_barrier
#endif
      end subroutine mp_scat3d

!-----
      subroutine mp_add1d(jdim, qin)
!-----
      implicit none
      integer jdim
      real  qin(jdim)
! Local:
      integer n, j
      integer dest, qsize
      integer inc
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif
#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4dwin, ierror)
      tdisp = jdim*gid
      qsize = jdim
!$omp parallel do private(n,dest)
      do n = 1,numpro
        dest = n-1
        call MPI_PUT(qin, qsize, mp_precision, dest, tdisp, &
                     qsize, mp_precision, buff4dwin, ierror)
      enddo
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4dwin, ierror)
#else
      call mpi_allgather(qin, jdim, mp_precision, buff4d, jdim, &
                         mp_precision, commglobal, ierror)
#endif
      do j=1,jdim
           qin(j) = 0.
        do n=1, numpro
           inc = (n-1)*jdim + j
           qin(j) = qin(j) + buff4d(inc)
        enddo
      enddo
#else
        do j=1,jdim
           g_t3(j,gid+1) = qin(j)       ! gid starts from 0!
        enddo

        call mlp_barrier(gid, gsize)

        do j=1,jdim
             qin(j) = 0.
          do n=1, numpro
             qin(j) = qin(j) + g_t3(j,n)
          enddo
        enddo

        call mlp_barrier(gid, gsize)
#endif
      end subroutine mp_add1d

!-----
      subroutine mp_gath4_r2d(im, jm, jfirst, jlast, qin, id)
!-----

! 4-BYTE version of mp_gath_r2d

      integer im, jm
      integer id        ! source ID
      integer jfirst, jlast
      real*4, intent(inout):: qin(im,jm)

      integer i, j, n
      integer j1, j2
      integer qsize
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4d_r4win, ierror)
#endif
      dest = id
      tdisp = (jfirst-1)*im
      call BufferPack2d_r4(qin, 1, im, 1, jm, 1, im, jfirst, jlast, &
                           buff4d_r4(tdisp+1))
      qsize = im*(jlast-jfirst+1)
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
      do p=1,pkgs_per_pro
        tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = tdisp + (p-1)*tmpsize
        call MPI_PUT(buff4d_r4(mydisp+1), mysize, MPI_REAL, dest, &
                     mydisp, mysize, MPI_REAL, buff4d_r4win, ierror)
      enddo
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4d_r4win, ierror)
#else
      send_tag = gid
      nsend = nsend + 1
      call mpi_isend(buff4d_r4(tdisp+1), qsize, MPI_REAL, dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
      if (gid==id) then
        do n=0,numpro-1
          j1 = yfirst(n+1)
          j2 = ylast(n+1)
          qsize = im*(j2-j1+1)
          src = n
          recv_tag = src
          tdisp = (j1-1)*im
          call mpi_recv(buff4d_r4(tdisp+1), qsize, MPI_REAL, src, &
                        recv_tag, commglobal, Status, ierror)
        enddo
      endif
#endif
      if (gid == id) then
        call BufferUnPack2d_r4(qin, 1, im, 1, jm, 1, im, 1, jm, buff4d_r4)
      endif

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
!$omp parallel do private(i, j)
          do j=jfirst,jlast
             do i=1,im
                g_2d(i,j) = qin(i,j)
             enddo
          enddo
      call mlp_barrier(gid, gsize)


      if ( gid == id ) then
!$omp parallel do private(i, j)
          do j=1,jm
             do i=1,im
                qin(i,j) = g_2d(i,j)
             enddo
          enddo
      endif
      call mlp_barrier(gid, gsize)
#endif
      end subroutine mp_gath4_r2d

!-----
      subroutine mp_scat4_r2d(qin, qout, im, jm, jfirst, jlast, id)
!-----
! Send 2D array qin from Process=id to All other Processes
      integer, intent(in):: im, jm
      integer, intent(in):: id         ! Source ID
      integer, intent(in):: jfirst, jlast
      real*4,  intent(in):: qin(im,jm)
      real*4, intent(out):: qout(im,jfirst:jlast)
 
! Local:
      integer i, j
      integer j1, j2
      integer n
      integer qsize
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP)
      integer rqst(numpro), rq_stats(numpro*MPI_STATUS_SIZE)
#if defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#endif
 
#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4d_r4win, ierror)
#endif
      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
           qsize = im*(j2-j1+1)
           dest = n-1
           tdisp = (j1-1)*im
           call BufferPack2d_r4(qin, 1, im, 1, jm, 1, im, j1, j2, &
                                buff4d_r4(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
           do p=1,pkgs_per_pro
             tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
             mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
             mydisp = tdisp + (p-1)*tmpsize
             call MPI_PUT(buff4d_r4(mydisp+1), mysize, MPI_REAL, dest, &
                          mydisp, mysize, MPI_REAL, buff4d_r4win, ierror)
           enddo
#else
           send_tag = gid
           nsend = nsend + 1
           call mpi_isend(buff4d_r4(tdisp+1), qsize, MPI_REAL, dest, &
                          send_tag, commglobal, rqst(nsend), ierror)
#endif
         enddo
      endif
      tdisp = (jfirst-1)*im
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4d_r4win, ierror)
#else
      qsize = im*(jlast-jfirst+1)
      src = id
      recv_tag = src
      call mpi_recv(buff4d_r4(tdisp+1), qsize, MPI_REAL, src, &
                    recv_tag, commglobal, Status, ierror)
#endif
      call BufferUnPack2d_r4(qout, 1, im, jfirst, jlast, &
                                   1, im, jfirst, jlast, &
                                   buff4d_r4(tdisp+1))

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, rqst, rq_stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
      if ( gid == id ) then
          do j=1,jm
             do i=1,im
                g_2d(i,j) = qin(i,j)
             enddo
          enddo
      endif
      call mlp_barrier(gid, gsize)
 
      do j=jfirst,jlast
         do i=1,im
            qout(i,j) = g_2d(i,j)
         enddo
      enddo
      call mlp_barrier(gid, gsize)
#endif
      end subroutine mp_scat4_r2d

!-----
      subroutine mp_allgather1d(jm, jfirst, jlast, qin, qout)
!-----
      implicit none
      integer jm
      integer jfirst, jlast
      real  qin(jfirst:jlast)
! Output:
      real  qout(jm)
! Local:
      integer j, n
      integer qsize
      integer dest
      integer qrsize(numpro), displ(numpro)
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4dwin, ierror)
      tdisp = (jfirst-1)
      qsize = jlast-jfirst+1
!$omp parallel do private(n,dest)
      do n=1,numpro
        dest = n-1
        call MPI_PUT(qin, qsize, mp_precision, dest, tdisp, qsize, &
                     mp_precision, buff4dwin, ierror)
      enddo
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4dwin, ierror)
      call BufferUnPack2d(qout, 1, 1, 1, jm, 1, 1, 1, jm, buff4d)
#else
      qsize = jlast - jfirst + 1
      do n=1,numpro
        qrsize(n) = (ylast(n) - yfirst(n) + 1)
        displ(n) = (yfirst(n) - 1)
      enddo
      call mpi_allgatherv(qin, qsize, mp_precision, qout, qrsize, &
                          displ, mp_precision, commglobal, ierror)
#endif

#else
#include "mlp_ptr.h"

      do j=jfirst,jlast
         g_1d(j) = qin(j)
      enddo
      call mlp_barrier(gid, gsize)

      do j=1,jm
         qout(j) = g_1d(j)
      enddo
#endif
      end subroutine mp_allgather1d
 
!-----
      subroutine mp_barrier
!-----
#if !defined(USE_MLP)
!BMP        call MPI_WIN_FENCE(0, buffwin, ierror)
!BMP        call mpi_barrier(commglobal, ierror)
#else
#include "mlp_ptr.h"
        call mlp_barrier(gid, gsize)
#endif
      end subroutine mp_barrier

!BW added

#if defined(USE_MLP)
#if defined ( PGI ) || defined ( LAHEY )
!-----
      subroutine mlp_barrier(id, isize)
!-----

      integer id, isize

#if defined (SEMA)
      call u_barrier(numpro, semid)
#else
      call mlp_barrier_(numpro)
#endif

      end subroutine mlp_barrier
#endif
#endif



#if defined(USE_MLP)
!.......
      subroutine Ga_put(datain, im, jm, km, jfirstm, jlastp, &
                        kfirst, klast, win)
      integer i, j, k
      integer im, jm, km, lm
      integer jfirstm, jlastp, kfirst, klast
      real :: datain(1:im, jfirstm:jlastp, kfirst:klast)
      pointer (win, data_g)
      real :: data_g(im, jm, km)
#include "mlp_ptr.h"



!$omp parallel do private(i,j,k)
      do k=kfirst, klast
         do j=jfirstm, jlastp
            do i=1, im
              data_g(i,j,k)=datain(i,j,k)
            enddo
         enddo
      enddo

      return
      end subroutine Ga_put



!.......
      subroutine Ga_get(dataout, im, jm, km, jfirstm, jlastp, &
                        kfirst, klast, win)
      integer i, j, k
      integer im, jm, km
      integer jfirstm, jlastp, kfirst, klast
      real :: dataout(1:im, jfirstm:jlastp, kfirst:klast)
      pointer (win, data_g)
      real :: data_g(im, jm, km)
#include "mlp_ptr.h"

!$omp parallel do private(i,j,k)
      do k=kfirst, klast
         do j=jfirstm, jlastp
            do i=1, im
              dataout(i,j,k)=data_g(i,j,k)
            enddo
         enddo
      enddo

      return
      end subroutine Ga_get
#endif


!......
      subroutine Ga_2_La(data_g, data_l, im, jm, km, jfirstm, jlastp, &
                         kfirst, klast)

      integer i, j, k
      integer im, jm, km
      integer jfirstm, jlastp, kfirst, klast
      real :: data_l(1:im, jfirstm:jlastp, kfirst:klast)
      real :: data_g(im, jm, km)

#if !defined(USE_MLP)
      call mp_scatter4d(data_g, data_l, im, jm, km, 1, jfirstm, jlastp,&
                        kfirst, klast, 0, 0, 0)
#else
#include "mlp_ptr.h"
      if (gid .eq. 0 ) then
         call Ga_put(data_g, im, jm, km, 1, jm, 1, km, wing_t1)
      endif
      call mp_barrier            !to a global array

      call Ga_get(data_l, im,  jm, km, jfirstm,jlastp, &
                  kfirst, klast, wing_t1)
      call mp_barrier            !to get its own copy for each PE
#endif
      return
      end subroutine Ga_2_La

!
      subroutine La_2_Ga(data_l, data_g, im, jm, km, jfirstm, jlastp,&
                         kfirst, klast)

      integer i, j, k
      integer im, jm, km
      integer jfirstm, jlastp, kfirst, klast
      real :: data_l(1:im, jfirstm:jlastp, kfirst:klast)
      real :: data_g(im, jm, km)

#if !defined(USE_MLP)
      call mp_gather4d(data_l, data_g, im, jm, km, 1, jfirstm, jlastp,&
                       kfirst, klast, 0, 0, 0)
#else
#include "mlp_ptr.h"
      call Ga_put(data_l, im, jm, km, jfirstm, jlastp, &
                  kfirst, klast, wing_t1)
      call mp_barrier            !to a global array

      if (gid .eq. 0 ) then
         call Ga_get(data_g, im,  jm, km, 1, jm, 1, km, wing_t1)
      endif
      call mp_barrier            !to get a whole array for each PE
#endif
      return
      end subroutine La_2_Ga

!-----
      subroutine mp_scatter2d(qin, qout, im, jm, jfirst, jlast, id)
!-----
! Send 2D array qin from Process=id to All other Processes
      integer, intent(in):: im, jm
      integer, intent(in):: id         ! Source ID
      integer, intent(in):: jfirst, jlast
      real, intent(in):: qin(im,jm)
      real, intent(out):: qout(im,jfirst:jlast)

! Local:
      integer i, j, k, n, iq
      integer j1, j2
      integer qsize_s
      integer qsize_r
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP)
      integer rqst(numpro), rq_stats(numpro*MPI_STATUS_SIZE)
#if defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#endif

#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4dwin, ierror)
#endif
      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
           qsize_s = im*(j2-j1+1)
           dest = n-1
           tdisp = (j1-1)*im
#if defined(MPI2)
           call BufferPack2d(qin, 1, im, 1, jm, 1, im, j1, j2, buff4d(tdisp+1))
!$omp parallel do private(p,tmpsize,mysize,mydisp)
           do p=1,pkgs_per_pro
             tmpsize = ceiling(real(qsize_s)/real(pkgs_per_pro))
             mysize = MIN(tmpsize, MAX(qsize_s-(tmpsize*(p-1)),0))
             mydisp = tdisp + (p-1)*tmpsize
             call MPI_PUT(buff4d(mydisp+1), mysize, mp_precision, dest, &
                          mydisp, mysize, mp_precision, buff4dwin, ierror)
           enddo
#else
           send_tag = gid
           nsend = nsend + 1
           call mpi_isend(qin(1,j1), qsize_s, mp_precision, dest, &
                          send_tag, commglobal, rqst(nsend), ierror)
#endif
         enddo
      endif
      tdisp = (jfirst-1)*im
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4dwin, ierror)
#else
      qsize_r = im*(jlast-jfirst+1)
      src = id
      recv_tag = src
      call mpi_recv(buff4d(tdisp+1), qsize_r, mp_precision, src, &
                    recv_tag, commglobal, Status, ierror)
#endif
      call BufferUnPack2d(qout, 1, im, jfirst, jlast, 1, im, jfirst, jlast, &
                          buff4d(tdisp+1))

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, rqst, rq_stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
      if ( gid == id ) then
          do j=1,jm
             do i=1,im
                g_2d(i,j) = qin(i,j)
             enddo
          enddo
      endif
      call mlp_barrier(gid, gsize)

      do j=jfirst,jlast
         do i=1,im
            qout(i,j) = g_2d(i,j)
         enddo
      enddo
      call mlp_barrier(gid, gsize)
#endif
      end subroutine mp_scatter2d

!-----
      subroutine mp_scatter4d(qin, qout, im, jm, km, nq, jfirst, jlast, &
                              kfirst, klast, ng_s, ng_n, id)
!-----

      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast, kfirst, klast
      integer, intent(in):: ng_s, ng_n
      integer, intent(in):: nq
      integer, intent(in):: id           ! Source (usually 0)
      real, intent(in)::   qin(im,jm,km,nq)
      real, intent(out):: qout(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)

! Local:
      integer i, j, k, n, iq
      integer j1, j2, k1, k2
      integer qsize_s
      integer qsize_r
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP)
      integer rqst(numpro), rq_stats(numpro*MPI_STATUS_SIZE)
#if defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#endif

#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4dwin, ierror)
#endif
      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
           k1 = zfirst(n)
           k2 = zlast(n)
           qsize_s = im*(j2-j1+1)*(k2-k1+1)*nq
           dest = n-1
           tdisp = ((j1-1)*km+(k1-1))*im*nq
           call BufferPack4d(qin, 1, im, 1,  jm, 1, km, 1, nq, &
                                  1, im, j1, j2, k1, k2, 1, nq, &
                                  buff4d(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
           do p=1,pkgs_per_pro
             tmpsize = ceiling(real(qsize_s)/real(pkgs_per_pro))
             mysize = MIN(tmpsize, MAX(qsize_s-(tmpsize*(p-1)),0))
             mydisp = tdisp + (p-1)*tmpsize
             call MPI_PUT(buff4d(mydisp+1), mysize, mp_precision, dest, &
                          mydisp, mysize, mp_precision, buff4dwin, ierror)
           enddo
#else
           send_tag = gid
           nsend = nsend + 1
           call mpi_isend(buff4d(tdisp+1), qsize_s, mp_precision, dest, &
                          send_tag, commglobal, rqst(nsend), ierror)
#endif
         enddo
      endif
      tdisp = ((jfirst-1)*km+(kfirst-1))*im*nq
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4dwin, ierror)
#else
      qsize_r = im*(jlast-jfirst+1)*(klast-kfirst+1)*nq
      src = id
      recv_tag = src
      call mpi_recv(buff4d(tdisp+1), qsize_r, mp_precision, src, &
                    recv_tag, commglobal, Status, ierror)
#endif
      call BufferUnPack4d(qout, 1, im, jfirst-ng_s, jlast+ng_n, kfirst,&
                          klast, 1, nq, 1, im, jfirst, jlast, kfirst,  &
                          klast, 1, nq, buff4d(tdisp+1))

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, rqst, rq_stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
      if ( gid == id ) then
      do iq=1,nq
!$omp parallel do private(i,j,k)
         do k=1,km
            do j=1,jm
               do i=1,im
                  g_4d(i,j,k,iq) = qin(i,j,k,iq)
               enddo
            enddo
         enddo
      enddo
      endif
      call mp_barrier

      do iq=1,nq
!$omp parallel do private(i,j,k)
         do k=kfirst,klast
            do j=jfirst,jlast
               do i=1,im
                  qout(i,j,k,iq) = g_4d(i,j,k,iq)
               enddo
            enddo
         enddo
      enddo
      call mp_barrier
#endif
      end subroutine mp_scatter4d

!-----
      subroutine mp_gather4d(qin, qout, im, jm, km, nq, jfirst, jlast, &
                             kfirst, klast, ng_s, ng_n, id)
!-----

      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast, kfirst, klast
      integer, intent(in):: ng_s, ng_n
      integer, intent(in):: nq
      integer, intent(in):: id         ! process ID to gather data to
      real, intent(in)  :: qin(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
      real, intent(out) :: qout(im,jm,km,nq)

! Local:
      integer i, j, k, iq
      integer j1, j2, k1, k2
      integer n
      integer qsize
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
 
#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4dwin, ierror)
#endif
      dest = id
      tdisp = ((jfirst-1)*km+(kfirst-1))*im*nq
      call BufferPack4d(qin, 1, im, jfirst-ng_s, jlast+ng_n, &
                        kfirst, klast, 1, nq, 1, im, jfirst, jlast,&
                        kfirst, klast, 1, nq, buff4d(tdisp+1))
      qsize = (jlast-jfirst+1)*im*(klast-kfirst+1)*nq
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
      do p=1,pkgs_per_pro
        tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = tdisp + (p-1)*tmpsize
        call MPI_PUT(buff4d(mydisp+1), mysize, mp_precision, dest, &
                     mydisp, mysize, mp_precision, buff4dwin, ierror)
      enddo
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4dwin, ierror)
#else
      send_tag = gid
      nsend = nsend + 1
      call mpi_isend(buff4d(tdisp+1), qsize, mp_precision, dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
#endif
      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
!
! BMP: fvGCM edit
! Rough hack to account for gathering arrays k=1:1 or k=1:km+1
!
           if (km == PLEV) then
             k1 = zfirst(n)
             k2 = zlast(n)
           else
             k1 = 1
             k2 = km
           endif
           tdisp = ((j1-1)*km+(k1-1))*im*nq

#if !defined(MPI2)
           qsize = im*(j2-j1+1)*(k2-k1+1)*nq
           src = n-1
           recv_tag = src
           call mpi_recv(buff4d(tdisp+1), qsize, mp_precision, src, &
                         recv_tag, commglobal, Status, ierror)
#endif
           call BufferUnPack4d(qout, 1, im, 1,  jm, 1, km, 1, nq, &
                                     1, im, j1, j2, k1, k2, 1, nq, &
                                     buff4d(tdisp+1))
         enddo
      endif

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
      do iq=1,nq
!$omp parallel do private(i,j,k)
         do k=kfirst,klast
            do j=jfirst,jlast
               do i=1,im
                  g_4d(i,j,k,iq) = qin(i,j,k,iq)
               enddo
            enddo
         enddo
      enddo
      call mp_barrier

      if ( gid == id ) then
      do iq=1,nq
!$omp parallel do private(i,j,k)
         do k=1,km
            do j=1,jm
               do i=1,im
                  qout(i,j,k,iq) = g_4d(i,j,k,iq)
               enddo
            enddo
         enddo
      enddo
      endif
      call mp_barrier
#endif
      end subroutine mp_gather4d
 
!-----
      subroutine mp_reduce_sum (sum)
!-----

      implicit none
      integer  n
      real sum
      real sumin
      real sumg(numpro)

#if !defined(USE_MLP)
      call mpi_gather(sum, 1, mp_precision, sumg, 1, mp_precision, 0, &
                      commglobal, ierror)
      if (gid == 0) then
        sum=0.
        do n=1,numpro
          sum = sum + sumg(n)
        enddo
      endif
      call mp_bcst_real(sum)
#else
#include "mlp_ptr.h"
      g_1d(gid+1) = sum           
      call mlp_barrier(gid, gsize)

      sum=0.0
      do n=1,numpro
         sum = sum + g_1d(n)
      enddo
      call mlp_barrier(gid, gsize) !may not be necessay, test, BW
#endif
      end subroutine mp_reduce_sum

!-----
      subroutine mp_wrt3d(iout, nrec, r_zero, qin, im, jm, km, &
                          jfirst, jlast, kfirst, klast, id )
!-----
      integer, intent(in):: iout
      integer, intent(in):: nrec
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast, kfirst, klast
      integer, intent(in):: id         ! process ID
      real, intent(in) :: r_zero
      real, intent(in) :: qin(im,jfirst:jlast,kfirst:klast)
! Local:
      integer i, j, k, iq
      real*4 a2(im,jm)
      integer n
      integer qsize
      integer src, dest
      integer send_tag, recv_tag
      real qout(im,jm,km)

#if !defined(USE_MLP)
      call mp_gather4d(qin, qout, im, jm, km, 1, jfirst, jlast, &
                       kfirst, klast, 0, 0, 0)
      if (gid == id) then
        do k=1,km
          do j=1,jm
            do i=1,im
              if( abs(qout(i,j,k)) < r_zero ) then
                 qout(i,j,k) = 0.
              endif
              a2(i,j) = qout(i,j,k)
            enddo
          enddo
        write(iout,rec=nrec+k) a2
        enddo
      endif
#else
#include "mlp_ptr.h"
!$omp parallel do private(i,j,k)
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=1,im
               if( abs(qin(i,j,k)) < r_zero ) then
                   g_t1(i,j,k,1) = 0.
               else
                   g_t1(i,j,k,1) = qin(i,j,k)
               endif
            enddo
         enddo
      enddo
      call mp_barrier

      if ( gid == id ) then
!-----------------------------------
! This is a crude paralle I/O attemp
!-----------------------------------
! It is actually slower than sequential version
!**** !$omp parallel do private(i,j,k,a2)
        do k=1,km
           do j=1,jm
              do i=1,im
! Convert to 32-bit
                 a2(i,j) = g_t1(i,j,k,1)
              enddo
           enddo
           write(iout,rec=nrec+k) a2
         enddo
      endif
      call mp_barrier
#endif
      end subroutine mp_wrt3d

#if !defined(USE_MLP)
      subroutine BufferPack4d ( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                   nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                   nq1, nq2, buff )
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real, intent(in)     :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      real, intent(out)    :: buff(i1:i2,j1:j2,k1:k2,nq1:nq2) ! Packed Buffer
! Local
      integer i, j, k, iq

      do iq = nq1, nq2
!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              buff(i,j,k,iq) = q(i,j,k,iq)
            enddo
          enddo
        enddo
      enddo
      end subroutine BufferPack4d

      subroutine BufferUnPack4d ( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                     nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                     nq1, nq2, buff )
      integer, intent(in) :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real, intent(in)    :: buff(i1:i2,j1:j2,k1:k2,nq1:nq2) ! Packed Buffer
      integer, intent(in) :: i1, i2, j1, j2, k1, k2, nq1, nq2
      real, intent(out)   :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
! Local
      integer i, j, k, iq

      do iq = nq1, nq2
!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              q(i,j,k,iq) = buff(i,j,k,iq)
            enddo
          enddo
        enddo
      enddo
      end subroutine BufferUnPack4d

      subroutine BufferPack3d ( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                   i1,    i2,  j1,    j2,  k1,    k2,  &
                                   buff )
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto 
      real, intent(in)     :: q(ifrom:ito,jfrom:jto,kfrom:kto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2
      real, intent(out)    :: buff(i1:i2,j1:j2,k1:k2) ! Packed Buffer
! Local
      integer i, j, k

!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              buff(i,j,k) = q(i,j,k)
            enddo
          enddo
        enddo
      end subroutine BufferPack3d

      subroutine BufferUnPack3d ( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                     i1,    i2,  j1,    j2,  k1,    k2,  &
                                     buff )
      integer, intent(in) :: ifrom, ito, jfrom, jto, kfrom, kto
      real, intent(in)    :: buff(i1:i2,j1:j2,k1:k2) ! Packed Buffer
      integer, intent(in) :: i1, i2, j1, j2, k1, k2
      real, intent(out)   :: q(ifrom:ito,jfrom:jto,kfrom:kto)
! Local
      integer i, j, k

!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              q(i,j,k) = buff(i,j,k)
            enddo
          enddo
        enddo
      end subroutine BufferUnPack3d

      subroutine BufferPack3d_i ( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                     i1,    i2,  j1,    j2,  k1,    k2,  &
                                     buff)
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto
      integer, intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2
      integer, intent(out) :: buff(i1:i2,j1:j2,k1:k2) ! Packed Buffer
! Local
      integer i, j, k

!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              buff(i,j,k) = q(i,j,k)
            enddo
          enddo
        enddo
      end subroutine BufferPack3d_i

      subroutine BufferUnPack3d_i ( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                       i1,    i2,  j1,    j2,  k1,    k2,  &
                                       buff)
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto
      integer, intent(in)  :: buff(i1:i2,j1:j2,k1:k2) ! Packed Buffer
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2
      integer, intent(out) :: q(ifrom:ito,jfrom:jto,kfrom:kto)
! Local
      integer i, j, k

!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              q(i,j,k) = buff(i,j,k)
            enddo
          enddo
        enddo
      end subroutine BufferUnPack3d_i

      subroutine BufferPack2d ( q, ifrom, ito, kfrom, kto, &
                                   i1,    i2,  k1,    k2,  &
                                   buff)
      integer, intent(in)  :: ifrom, ito, kfrom, kto
      real, intent(in)     :: q(ifrom:ito,kfrom:kto)
      integer, intent(in)  :: i1, i2, k1, k2
      real, intent(out)    :: buff(i1:i2,k1:k2) ! Packed Buffer
! Local
      integer i, k

!$omp parallel do private(i,k)
        do k = k1, k2
            do i = i1, i2
              buff(i,k) = q(i,k)
            enddo
        enddo
      end subroutine BufferPack2d

      subroutine BufferUnPack2d ( q, ifrom, ito, kfrom, kto, &
                                     i1,    i2,  k1,    k2,  &
                                     buff ) 
      integer, intent(in) :: ifrom, ito, kfrom, kto
      real, intent(in)    :: buff(i1:i2,k1:k2) ! Packed Buffer
      integer, intent(in) :: i1, i2, k1, k2
      real, intent(out)   :: q(ifrom:ito,kfrom:kto)
! Local
      integer i, k

!$omp parallel do private(i,k)
        do k = k1, k2
            do i = i1, i2
              q(i,k) = buff(i,k)
            enddo
        enddo
      end subroutine BufferUnPack2d

      subroutine BufferPack2d_i ( q, ifrom, ito, kfrom, kto, &
                                     i1,    i2,  k1,    k2,  &
                                     buff)
      integer, intent(in)  :: ifrom, ito, kfrom, kto
      integer, intent(in)  :: q(ifrom:ito,kfrom:kto)
      integer, intent(in)  :: i1, i2, k1, k2
      integer, intent(out) :: buff(i1:i2,k1:k2) ! Packed Buffer
! Local
      integer i, k

!$omp parallel do private(i,k)
        do k = k1, k2
            do i = i1, i2
              buff(i,k) = q(i,k)
            enddo
        enddo
      end subroutine BufferPack2d_i

      subroutine BufferUnPack2d_i ( q, ifrom, ito, kfrom, kto, &
                                       i1,    i2,  k1,    k2,  &
                                       buff)
      integer, intent(in)  :: ifrom, ito, kfrom, kto
      integer, intent(in)  :: buff(i1:i2,k1:k2) ! Packed Buffer
      integer, intent(in)  :: i1, i2, k1, k2
      integer, intent(out) :: q(ifrom:ito,kfrom:kto)
! Local
      integer i, k

!$omp parallel do private(i,k)
        do k = k1, k2
            do i = i1, i2
              q(i,k) = buff(i,k)
            enddo
        enddo
      end subroutine BufferUnPack2d_i

      subroutine BufferPack2d_r4 ( q, ifrom, ito, kfrom, kto, &
                                      i1,    i2,  k1,    k2,  &
                                      buff)
      integer, intent(in)  :: ifrom, ito, kfrom, kto
      real*4, intent(in)   :: q(ifrom:ito,kfrom:kto)
      integer, intent(in)  :: i1, i2, k1, k2
      real*4, intent(out)  :: buff(i1:i2,k1:k2) ! Packed Buffer
! Local
      integer i, k

!$omp parallel do private(i,k)
        do k = k1, k2
            do i = i1, i2
              buff(i,k) = q(i,k)
            enddo
        enddo
      end subroutine BufferPack2d_r4

      subroutine BufferUnPack2d_r4 ( q, ifrom, ito, kfrom, kto, &
                                        i1,    i2,  k1,    k2,  &
                                        buff)
      integer, intent(in)  :: ifrom, ito, kfrom, kto
      real*4, intent(in)   :: buff(i1:i2,k1:k2) ! Packed Buffer
      integer, intent(in)  :: i1, i2, k1, k2
      real*4, intent(out)  :: q(ifrom:ito,kfrom:kto)
! Local
      integer i, k

!$omp parallel do private(i,k)
        do k = k1, k2
            do i = i1, i2
              q(i,k) = buff(i,k)
            enddo
        enddo
      end subroutine BufferUnPack2d_r4
#endif

!=========================================================================
! FastOpt
! These two routines are wrappers to include three calls in cd_core_ad
! that would otherwise be missing.
!=========================================================================
subroutine mp_send3d_ns2(im, jm, jfirst, jlast, kfirst, klast, ng_s, ng_n, q, iq)
  implicit none
  integer im, jm
  integer jfirst, jlast
  integer kfirst, klast
  integer ng_s      ! southern zones to ghost 
  integer ng_n      ! noruthern zones to ghost 
  real q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)
  integer iq
  call mp_send3d_ns(im, jm, jfirst, jlast, kfirst, klast, ng_s, ng_n, q, iq)
end subroutine mp_send3d_ns2

subroutine mp_recv3d_ns2(im, jm, jfirst, jlast, kfirst, klast, ng_s, ng_n, q, iq)
  implicit none
  integer im, jm
  integer jfirst, jlast
  integer kfirst, klast
  integer ng_s      ! southern zones to ghost 
  integer ng_n      ! noruthern zones to ghost 
  integer iq        ! Counter
  real q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)
  call mp_recv3d_ns(im, jm, jfirst, jlast, kfirst, klast, ng_s, ng_n, q, iq)
end subroutine mp_recv3d_ns2

!=========================================================================
!=========================================================================
! adjoint
!=========================================================================
!=========================================================================
#if !defined(USE_MLP)

      subroutine BufferPack2dx( q, ifrom, ito, kfrom, kto, &
                                   i1,    i2,  k1,    k2,  &
                                   buff)
      integer, intent(in)    :: ifrom, ito, kfrom, kto
      real   , intent(inout) :: q(ifrom:ito,kfrom:kto)
      integer, intent(in)    :: i1, i2, k1, k2
      real   , intent(out)   :: buff(i1:i2,k1:k2) ! Packed Buffer
! Local
      integer i, k

!$omp parallel do private(i,k)
        do k = k1, k2
            do i = i1, i2
              buff(i,k) = q(i,k)
	      q(i,k)    = 0.0                     ! FastOpt: reset
            enddo
        enddo
      end subroutine BufferPack2dx


      subroutine BufferUnPack2dx( q, ifrom, ito, kfrom, kto, &
                                     i1,    i2,  k1,    k2,  &
                                     buff ) 
      integer, intent(in)  :: ifrom, ito, kfrom, kto
      real   , intent(in)  :: buff(i1:i2,k1:k2) ! Packed Buffer
      integer, intent(in)  :: i1, i2, k1, k2
      real   , intent(out) :: q(ifrom:ito,kfrom:kto)
! Local
      integer i, k

!$omp parallel do private(i,k)
        do k = k1, k2
            do i = i1, i2
              q(i,k) = q(i,k) + buff(i,k)         ! FastOpt: add
            enddo
        enddo
      end subroutine BufferUnPack2dx


      subroutine BufferPack3dx ( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                    i1,    i2,  j1,    j2,  k1,    k2,  &
                                    buff )
!=========================================================================
      integer, intent(in)    :: ifrom, ito, jfrom, jto, kfrom, kto 
      real   , intent(inout) :: q(ifrom:ito,jfrom:jto,kfrom:kto)
      integer, intent(in)    :: i1, i2, j1, j2, k1, k2
      real   , intent(out)   :: buff(i1:i2,j1:j2,k1:k2) ! Packed Buffer
! Local
      integer i, j, k

!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              buff(i,j,k) = q(i,j,k)
              q(i,j,k)    = 0.0                   ! FastOpt: reset
            enddo
          enddo
        enddo
      end subroutine BufferPack3dx

!=========================================================================
      subroutine BufferUnPack3dx ( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                      i1,    i2,  j1,    j2,  k1,    k2,  &
                                      buff )
!=========================================================================
      integer, intent(in)    :: ifrom, ito, jfrom, jto, kfrom, kto
      real   , intent(in)    :: buff(i1:i2,j1:j2,k1:k2) ! Packed Buffer
      integer, intent(in)    :: i1, i2, j1, j2, k1, k2
      real   , intent(inout) :: q(ifrom:ito,jfrom:jto,kfrom:kto)
! Local
      integer i, j, k

!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              q(i,j,k) = q(i,j,k) + buff(i,j,k)   ! FastOpt: add
            enddo
          enddo
        enddo
      end subroutine BufferUnPack3dx


      subroutine BufferPack4dx( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                   nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                   nq1, nq2, buff )
      integer, intent(in)    :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real   , intent(inout) :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)    :: i1, i2, j1, j2, k1, k2, nq1, nq2
      real   , intent(out)   :: buff(i1:i2,j1:j2,k1:k2,nq1:nq2) ! Packed Buffer
! Local
      integer i, j, k, iq

      do iq = nq1, nq2
!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              buff(i,j,k,iq) = q(i,j,k,iq)
              q(i,j,k,iq)    = 0.0                   ! FastOpt: reset
            enddo
          enddo
        enddo
      enddo
      end subroutine BufferPack4dx


      subroutine BufferUnPack4dx( q, ifrom, ito, jfrom, jto, kfrom, kto, &
                                     nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                     nq1, nq2, buff )
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real   , intent(in)  :: buff(i1:i2,j1:j2,k1:k2,nq1:nq2) ! Packed Buffer
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      real   , intent(out) :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
! Local
      integer i, j, k, iq

      do iq = nq1, nq2
!$omp parallel do private(i,j,k)
        do k = k1, k2
          do j = j1, j2
            do i = i1, i2
              q(i,j,k,iq) = q(i,j,k,iq) + buff(i,j,k,iq)   ! FastOpt: add
            enddo
          enddo
        enddo
      enddo
      end subroutine BufferUnPack4dx


#endif /* !defined(USE_MLP) */


!=========================================================================
      subroutine mp_recv2_ns_ad(im, jm, jfirst, jlast, kfirst, klast, nd, q1, q2)
!=========================================================================
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
      real   , intent(inout):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real   , intent(inout):: q2(im,jfirst-nd:jlast+nd,kfirst:klast) 

! Local:
      integer i, j, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
      integer displ
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif

! Send to south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
! Start recv to north
        src = gid - 1
        recv_tag = src
        qsize = im*2*(klast-kfirst+1)
        nrecv = nrecv + 1
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        j = jfirst - 1
        dest = gid - 1
        qsize = im*(klast-kfirst+1)*2
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3dx(q1, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                               1, im, j,         j,        kfirst, klast, &
                               buff_s(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferPack3dx(q2, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                               1, im, j,         j,        kfirst, klast, &
                               buff_s(displ+tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif

! Send to north
      if ( jlast < jm ) then
#if !defined(MPI2)
! Start recv to south
        src = gid + 1
        recv_tag = src
        qsize = im*2*(klast-kfirst+1)
        nrecv = nrecv + 1
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        j = jlast + 1
        dest = gid + 1
        qsize = im*(klast-kfirst+1)*2
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3dx(q1, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                               1, im, j,         j,        kfirst, klast, &
                               buff_s(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferPack3dx(q2, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                               1, im, j,         j,        kfirst, klast, &
                               buff_s(displ+tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif

#else
#include "mlp_ptr.h"
!$omp parallel do private(i, k)
      do k=kfirst,klast
! Send to south
      if ( jfirst > 1 ) then
        j = jfirst - 1
           do i=1,im
              g_t1(i,j,k,1) = q1(i,j,k)
	      q1(i,j,k)     = 0.               ! FastOpt: reset
              g_t1(i,j,k,2) = q2(i,j,k)
	      q2(i,j,k)     = 0.               ! FastOpt: reset
           enddo
      endif

! Send to north
      if ( jlast < jm ) then
        j = jlast + 1
           do i=1,im
              g_t1(i,j,k,1) = q1(i,j,k)
	      q1(i,j,k)     = 0.               ! FastOpt: reset
              g_t1(i,j,k,2) = q2(i,j,k)
	      q2(i,j,k)     = 0.               ! FastOpt: reset
           enddo
      endif
      enddo
#endif
      end subroutine mp_recv2_ns_ad


!=========================================================================
      subroutine mp_send2_ns_ad(im, jm, jfirst, jlast, kfirst, klast, nd, q1, q2)
!=========================================================================
      implicit none
      integer, intent(in) :: im, jm, jfirst, jlast
      integer, intent(in) :: kfirst, klast !careful: klast might be klast+1
      integer, intent(in) :: nd
      real, intent(inout) :: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real, intent(inout) :: q2(im,jfirst-nd:jlast+nd,kfirst:klast) 

! Local:
      integer i,j, k, n
      integer displ

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
! Recv from south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        j = jfirst
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3dx(q1, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                                 1, im, j,         j,        kfirst, klast, &
                                 buff_r(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferUnPack3dx(q2, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                                 1, im, j,         j,        kfirst, klast, &
                                 buff_r(displ+tdisp+1))
      endif
! Recv from north
      if ( jlast < jm ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        j = jlast
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3dx(q1, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                                 1, im, j,         j,        kfirst, klast, &
                                buff_r(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferUnPack3dx(q2, 1, im, jfirst-nd, jlast+nd, kfirst, klast, &
                                 1, im, j,         j,        kfirst, klast, &
                                 buff_r(displ+tdisp+1))
      endif
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#else
#include "mlp_ptr.h"
!$omp parallel do private(i, j, k)
      do k=kfirst,klast
! Recv data from south
      if ( jfirst > 1 ) then
            j = jfirst
            do i=1,im
               q1(i,j,k) = q1(i,j,k) + g_t1(i,j,k,1)   ! FastOpt: add 
               q2(i,j,k) = q2(i,j,k) + g_t1(i,j,k,2)   ! FastOpt: add 
            enddo
      endif

! Recv data from north
      if ( jlast < jm ) then
            j = jlast
            do i=1,im
               q1(i,j,k) = q1(i,j,k) + g_t1(i,j,k,1)    ! FastOpt: add
               q2(i,j,k) = q2(i,j,k) + g_t1(i,j,k,2)    ! FastOpt: add
            enddo
      endif
      enddo
#endif
      end subroutine mp_send2_ns_ad


!=========================================================================
      subroutine mp_recv3d_ns_ad(im, jm, jfirst, jlast, kfirst, klast, & 
                                ng_s, ng_n, q, iq)
!=========================================================================
      implicit none
      integer im, jm
      integer jfirst, jlast
      integer kfirst, klast
      integer ng_s      ! southern zones to ghost 
      integer ng_n      ! noruthern zones to ghost 
      real    q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)
      integer iq
! Local:
      integer i,j,k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif
! Send to south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        src = gid - 1
        recv_tag = src
        qsize = im*ng_s*(klast-kfirst+1)
        nrecv = nrecv + 1
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif

        dest = gid - 1
        qsize = im*ng_n*(klast-kfirst+1)
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf 
        call BufferPack3dx(q, 1, im, jfirst-ng_s, jlast +ng_n, kfirst, klast, &
                              1, im, jfirst-ng_s, jfirst-1   , kfirst, klast, &
                              buff_s(tdisp+1))

#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(MPI2)
        src = gid + 1
        recv_tag = src
        qsize = im*ng_n*(klast-kfirst+1)
        nrecv = nrecv + 1
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        dest = gid + 1
        qsize = im*ng_s*(klast-kfirst+1)
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf 
        call BufferPack3dx(q, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                              1, im, jlast+1,     jlast+ng_n, kfirst, klast, &
                              buff_s(tdisp+1))

#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif
#else
#include "mlp_ptr.h"
!$omp parallel do private(i,j,k)
      do k=kfirst,klast
         if ( jfirst > 1 ) then
! Send to south
            do j=jfirst-ng_s, jfirst-1
               do i=1,im
                  g_4d(i,j,k,iq) = q(i,j,k)
		  q(i,j,k)       = 0.        ! FastOpt: reset
               enddo
            enddo
         endif
         if ( jlast < jm ) then
! Send to north
            do j=jlast+1,jlast+ng_n
               do i=1,im
                  g_4d(i,j,k,iq) = q(i,j,k)
		  q(i,j,k)       = 0.        ! FastOpt: reset
               enddo
            enddo
         endif
      enddo
#endif
      end subroutine mp_recv3d_ns_ad
 
!=========================================================================
      subroutine mp_send3d_ns_ad(im, jm, jfirst, jlast, kfirst, klast, &
                                 ng_s, ng_n, q, iq)
!=========================================================================
      implicit none
      integer im, jm
      integer jfirst, jlast
      integer kfirst, klast
      integer ng_s      ! southern zones to ghost 
      integer ng_n      ! noruthern zones to ghost 
      integer iq        ! Counter
      real    q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)

! Local:
      integer i,j,k
      integer src
      integer recv_tag

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
! Recv from south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3dx(q, 1, im, jfirst-ng_s, jlast +ng_n  , kfirst, klast, &
                                1, im, jfirst     , jfirst+ng_n-1, kfirst, klast, &
                                buff_r(tdisp+1))
      endif
! Recv from north
      if ( jlast < jm ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3dx(q, 1, im, jfirst-ng_s  , jlast+ng_n, kfirst, klast, &
                                1, im, jlast -ng_s+1, jlast     , kfirst, klast, &
                                buff_r(tdisp+1))
      endif
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#else
#include "mlp_ptr.h"
!$omp parallel do private(i,j,k)
      do k=kfirst,klast
         if ( jfirst > 1 ) then
! Recv from south
            do j=jfirst, jfirst+ng_n-1
               do i=1,im
                  q(i,j,k) = q(i,j,k) + g_4d(i,j,k,iq)   ! FastOpt: add
               enddo
            enddo
         endif

         if ( jlast < jm ) then
! Recv from north
            do j=jlast-ng_s+1, jlast
               do i=1,im
                  q(i,j,k) = q(i,j,k) + g_4d(i,j,k,iq)   ! FastOpt: add
               enddo
            enddo
         endif
      enddo
#endif
      end subroutine mp_send3d_ns_ad

!============================================================================
!============================================================================

!============================================================================
      subroutine mp_recv4d_ns_ad(im, jm, jfirst, jlast, kfirst, klast, &
                                 nq, ng_s, ng_n, q)
!============================================================================
      implicit none
      integer im, jm
      integer jfirst, jlast
      integer kfirst, klast
      integer nq
      integer ng_s      ! southern zones to ghost 
      integer ng_n      ! noruthern zones to ghost 
      real    q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
! Local:
      integer iq        ! Counter
      integer i, j, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif
! Send to south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        src = gid - 1
        recv_tag = src
        qsize = im*ng_s*(klast-kfirst+1)*nq
        nrecv = nrecv + 1
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        dest = gid - 1
        qsize = im*ng_n*(klast-kfirst+1)*nq
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack4dx(q, 1, im, jfirst-ng_s, jlast+ng_n, &
                              kfirst, klast, 1, nq, &
                              1, im, jfirst-ng_s, jfirst-1, &
                              kfirst, klast, 1 , nq, &
                              buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(MPI2)
        src = gid + 1
        recv_tag = src
        qsize = im*ng_n*(klast-kfirst+1)*nq
        nrecv = nrecv + 1
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
#endif
        dest = gid + 1
        qsize = im*ng_s*(klast-kfirst+1)*nq
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack4dx(q, 1, im, jfirst-ng_s, jlast+ng_n, &
                              kfirst, klast, 1, nq, &
                              1, im, jlast+1, jlast+ng_n, &
                              kfirst, klast, 1, nq, &
                              buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
      endif
#else
#include "mlp_ptr.h"
      do iq=1,nq
!$omp parallel do private(i,j,k)
      do k=kfirst,klast
         if ( jfirst > 1 ) then
! Send to south
            do j=jfirst-ng_s, jfirst-1
               do i=1,im
                  g_4d(i,j,k,iq) = q(i,j,k,iq)
		  q(i,j,k,iq)    = 0.        ! FastOpt: reset
               enddo
            enddo
         endif
         if ( jlast < jm ) then
! Send to north
            do j=jlast+1, jlast+ng_n
               do i=1,im
                  g_4d(i,j,k,iq) = q(i,j,k,iq)
		  q(i,j,k,iq)    = 0.        ! FastOpt: reset
               enddo
            enddo
         endif
      enddo
      enddo
#endif
      end subroutine mp_recv4d_ns_ad
 
!=========================================================================
      subroutine mp_send4d_ns_ad(im, jm, jfirst, jlast, kfirst, klast, &
                                 nq, ng_s, ng_n, q)
!=========================================================================
      implicit none
      integer im, jm
      integer jfirst, jlast
      integer kfirst, klast
      integer nq
      integer ng_s      ! southern zones to ghost 
      integer ng_n      ! noruthern zones to ghost 
      real q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
! Local:
      integer iq        ! Counter
      integer i, j, k

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
! Recv from south
      if ( jfirst > 1 ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack4dx(q, 1, im, jfirst-ng_s, jlast+ng_n, &
                                kfirst, klast, 1, nq, &
                                1, im, jfirst, jfirst+ng_n-1, &
                                kfirst, klast, 1, nq, &
                                buff_r(tdisp+1))
      endif
! Recv from north
      if ( jlast < jm ) then
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack4dx(q, 1, im, jfirst-ng_s, jlast+ng_n, &
                                kfirst, klast, 1, nq, &
                                1, im, jlast-ng_s+1, jlast, &
                                kfirst, klast, 1, nq, &
                                buff_r(tdisp+1))
      endif
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#else
#include "mlp_ptr.h"
      do iq=1,nq

!$omp parallel do private(i,j,k)
      do k=kfirst,klast
         if ( jfirst > 1 ) then
! Recv from south
            do j=jfirst, jfirst+ng_n-1
               do i=1,im
                  q(i,j,k,iq) = q(i,j,k,iq) + g_4d(i,j,k,iq)   ! FastOpt: add
               enddo
            enddo
         endif

         if ( jlast < jm ) then
! Recv from north
            do j=jlast-ng_s+1, jlast
               do i=1,im
                  q(i,j,k,iq) = q(i,j,k,iq) + g_4d(i,j,k,iq)   ! FastOpt: add
               enddo
            enddo
         endif
      enddo
      enddo
#endif
      end subroutine mp_send4d_ns_ad

!============================================================================
!============================================================================


!============================================================================
      subroutine mp_recv2_n_ad(im, jm, jfirst, jlast, kfirst, klast, &
                               ng_s, ng_n, q1, q2)
!============================================================================
      implicit none
      integer, intent(in) :: im, jm, jfirst, jlast
      integer, intent(in) :: kfirst, klast  ! careful: might be klast+1 outside
      integer, intent(in) :: ng_s
      integer, intent(in) :: ng_n
      real, intent(inout) :: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      real, intent(inout) :: q2(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      integer i, j, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
      integer displ
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#else
! Start receive from south
      if ( jfirst > 1 ) then
         src = gid - 1
         recv_tag = src
         qsize = im*2*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
#endif

! Send data to north
      if ( jlast < jm ) then
        j = jlast + 1
#if !defined(USE_MLP)
        dest = gid + 1
        qsize = im*(klast-kfirst+1)*2
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3dx(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                               1, im, j,           j,          kfirst, klast, &
                               buff_s(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferPack3dx(q2, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                               1, im, j,           j,          kfirst, klast, &
                               buff_s(displ+tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif

#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_t1(i,j,k,1) = q1(i,j,k)
	       q1(i,j,k)     = 0.              ! FastOpt reset
               g_t1(i,j,k,2) = q2(i,j,k)
	       q2(i,j,k)     = 0.              ! FastOpt reset
            enddo
         enddo
#endif
      endif
      end subroutine mp_recv2_n_ad


!============================================================================
      subroutine mp_recv2_s_ad(im, jm, jfirst, jlast, kfirst, klast, &
                             ng_s, ng_n, q1, q2)
!============================================================================
      implicit none
      integer, intent(in) :: im, jm, jfirst, jlast
      integer, intent(in) :: kfirst, klast   !careful: klast might be klast+1
      integer, intent(in) :: ng_s
      integer, intent(in) :: ng_n
      real, intent(inout) :: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      real, intent(inout) :: q2(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 

! Local:
      integer i, j, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
      integer displ
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#else
! Start recv from north
      if ( jlast < jm ) then
         src = gid + 1
         recv_tag = src
         qsize = im*2*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
#endif

! Send data to south
      if ( jfirst > 1 ) then
        j = jfirst - 1
#if !defined(USE_MLP)
        dest = gid - 1
        qsize = im*(klast-kfirst+1)*2
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3dx(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                               1, im, j,           j,          kfirst, klast, &
                               buff_s(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferPack3dx(q2, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                               1, im, j,           j,          kfirst, klast, &
                               buff_s(displ+tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend  = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif

#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_t1(i,j,k,1) = q1(i,jt,k)
	       q1(i,j,k)     = 0.              ! FastOpt reset
               g_t1(i,j,k,2) = q2(i,j,k)
	       q2(i,j,k)     = 0.              ! FastOpt reset
            enddo
         enddo
#endif
      endif
      end subroutine mp_recv2_s_ad

!============================================================================
      subroutine mp_send2_s_ad(im, jm, jfirst, jlast, kfirst, klast, &
                             ng_s, ng_n, q1, q2)
!============================================================================
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast     ! careful: klast might be klast+1
      integer, intent(in):: ng_s
      integer, intent(in):: ng_n
      real, intent(inout):: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      real, intent(inout):: q2(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 

! Local:
      integer i,j, k, n
      integer displ

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv data from south
      if ( jfirst > 1 ) then
        j = jfirst
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3dx(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                                 1, im, j,           j,          kfirst, klast, &
                                 buff_r(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferUnPack3dx(q2, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                                 1, im, j,           j,          kfirst, klast, &
                                 buff_r(displ+tdisp+1))
#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               q1(i,j,k) = q1(i,j,k) + g_t1(i,j,k,1)   ! FastOpt add
               q2(i,j,k) = q2(i,j,k) + g_t1(i,j,k,2)   ! FastOpt add
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_send2_s_ad


!============================================================================
      subroutine mp_send2_n_ad(im, jm, jfirst, jlast, kfirst, klast, &
                               ng_s, ng_n, q1, q2)
!============================================================================
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast     !careful: klast might be klast+1
      integer, intent(in):: ng_s
      integer, intent(in):: ng_n
      real, intent(inout):: q1(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 
      real, intent(inout):: q2(im,jfirst-ng_s:jlast+ng_n,kfirst:klast) 

! Local:
      integer i,j, k, n
      integer displ

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv data from north
      if ( jlast < jm ) then
        j = jlast
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3dx(q1, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                                 1, im, j,           j,          kfirst, klast, &
                                 buff_r(tdisp+1))
        displ = im*(klast-kfirst+1)
        call BufferUnPack3dx(q2, 1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, &
                                 1, im, j,           j,          kfirst, klast, &
                                 buff_r(displ+tdisp+1))
#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               q1(i,j,k) = g_t1(i,j,k,1) 
               q2(i,j,k) = g_t1(i,j,k,2) 
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_send2_n_ad


!============================================================================
!============================================================================

!============================================================================
      subroutine mp_recv_n_ad(im, jm, jfirst, jlast, kfirst, klast, &
                              nd_s, nd_n, q)
!============================================================================
      implicit none
      integer, intent(in)    :: im, jm, jfirst, jlast
      integer, intent(in)    :: kfirst, klast
      integer, intent(in)    :: nd_s, nd_n
      real   , intent(inout) :: q(im,jfirst-nd_s:jlast+nd_n,kfirst:klast) 

! Local:
      integer i, j, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#else
! Start receive from south
      if ( jfirst > 1 ) then
         src = gid - 1
         recv_tag = src
         qsize = im*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
#endif

! Send data to north
      if ( jlast < jm ) then
        j = jlast + 1
#if !defined(USE_MLP)
        dest = gid + 1
        qsize = im*(klast-kfirst+1)
        tdisp = igosouth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack3dx(q, 1, im, jfirst-nd_s, jlast+nd_n, kfirst, klast, &
                              1, im, j,           j,          kfirst, klast, &
                              buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
        do p=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
          mydisp = tdisp + (p-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
         send_tag = gid
         nsend = nsend + 1
         call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                        send_tag, commglobal, sqest(nsend), ierror)
#endif

#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_4d(i,j,k,1) = q(i,j,k)
	       q(i,j,k)      = 0.              ! FastOpt reset
            enddo
         enddo
#endif
      endif
      end subroutine mp_recv_n_ad

!============================================================================
      subroutine mp_send_s_ad(im, jm, jfirst, jlast, kfirst, klast, nd_s, nd_n, q)
!============================================================================
      implicit none
      integer, intent(in):: im, jm, jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: nd_s, nd_n
      real, intent(inout):: q(im,jfirst-nd_s:jlast+nd_n,kfirst:klast) 

! Local:
      integer i, j, k, n

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv data from south
      if ( jfirst > 1 ) then
        j = jfirst
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igosouth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3dx(q, 1, im, jfirst-nd_s, jlast+nd_n, kfirst, klast, &
                                1, im, j,           j,          kfirst, klast, &
                                buff_r(tdisp+1))
#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               q(i,j,k) = q(i,j,k) + g_4d(i,j,k,1)     ! FastOpt add
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_send_s_ad

!============================================================================
!============================================================================

!============================================================================
      subroutine mp_recv_pe_ad(im, jm, jfirst, jlast, kfirst, klast, pesouth)
!============================================================================
      implicit none
      integer, intent(in)    :: im, jm, jfirst, jlast, kfirst, klast
      real   , intent(inout) :: pesouth(im,kfirst:klast) 
      integer i, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer n, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif
#endif

#if !defined(USE_MLP) && !defined(MPI2)
! Start recv from North
      if ( jlast < jm ) then
         src = gid + 1
         recv_tag = src
         qsize = im*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
! Send data to South
      if ( jfirst > 1 ) then
#if !defined(USE_MLP)
        dest = gid - 1
        qsize = im*(klast-kfirst+1)
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack2dx(pesouth, 1, im, kfirst, klast, &
                                    1, im, kfirst, klast, &
                                    buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(n,tmpsize,mysize,mydisp)
        do n=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(n-1)),0))
          mydisp = tdisp + (n-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif

#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_t2(i,k,jfirst-1) = pesouth(i,k)
	       pesouth(i,k)       = 0.              ! FastOpt reset
            enddo
         enddo
#endif
      endif
      end subroutine mp_recv_pe_ad


!============================================================================
      subroutine mp_send_pe_ad(im, jm, jfirst, jlast, kfirst, klast, p)
!============================================================================
      implicit none
      integer, intent(in)    :: im, jm, jfirst, jlast, kfirst, klast
      real   , intent(inout) :: p(im,kfirst:klast,jfirst:jlast) 
      integer i, j, k, n
#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv from North
      if ( jlast < jm ) then
         j = jlast
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack2dx(p(1,1,j), 1, im, kfirst, klast, &
                                       1, im, kfirst, klast, &
                                       buff_r(tdisp+1))
#else
!$omp parallel do private(i, j, k)
         do k=kfirst,klast
            do i=1,im
	       p(i,k,j) = p(i,k,j) + g_t2(i,k,j)  ! FastOpt add
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_send_pe_ad

!============================================================================
!============================================================================

!============================================================================
      subroutine mp_recv_ua_ad(im, jm, jfirst, jlast, kfirst, klast, uasouth)
!============================================================================
      implicit none
      integer, intent(in)    :: im, jm, jfirst, jlast, kfirst, klast
      real   , intent(inout) :: uasouth(im, kfirst:klast) 
      integer i, k
      integer src, dest
      integer qsize
      integer recv_tag, send_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer n, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif

#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_s = ncall_s + 1
#if defined(MPI2)
      if (ncall_s == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buffwin, ierror)
      endif
#endif
#endif

#if !defined(USE_MLP) && !defined(MPI2)
! Start recv from north
      if ( jlast < jm ) then
         src = gid + 1
         recv_tag = src
         qsize = im*(klast-kfirst+1)
         nrecv = nrecv + 1
         tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
         call mpi_irecv(buff_r(tdisp+1), qsize, mp_precision, src, &
                        recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif
! Send data to South
      if ( jfirst > 1 ) then
#if !defined(USE_MLP)
        dest = gid - 1
        qsize = im*(klast-kfirst+1)
        tdisp = igonorth*idimsize + (ncall_s-1)*idimsize*nbuf
        call BufferPack2dx(uasouth, 1, im, kfirst, klast, &
                                    1, im, kfirst, klast, &
                                    buff_s(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(n,tmpsize,mysize,mydisp)
        do n=1,pkgs_per_pro
          tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
          mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(n-1)),0))
          mydisp = tdisp + (n-1)*tmpsize
          call MPI_PUT(buff_s(mydisp+1), mysize, mp_precision, dest, &
                       mydisp, mysize, mp_precision, buffwin, ierror)
        enddo
#else
        send_tag = gid
        nsend = nsend + 1
        call mpi_isend(buff_s(tdisp+1), qsize, mp_precision, dest, &
                       send_tag, commglobal, sqest(nsend), ierror)
#endif
#else
!$omp parallel do private(i, k)
         do k=kfirst,klast
            do i=1,im
               g_4d(i,jfirst-1,k,1) = uasouth(i,k)
	       uasouth(i,k)         = 0.              ! FastOpt reset
            enddo
         enddo
#endif
      endif
      end subroutine mp_recv_ua_ad


!============================================================================
      subroutine mp_send_ua_ad(im, jm, jfirst, jlast, kfirst, klast, p)
!============================================================================
      implicit none
      integer, intent(in)    :: im, jm, jfirst, jlast, kfirst, klast
      real   , intent(inout) :: p(im,jfirst:jlast,kfirst:klast) 
      integer i, j, k, n
#if defined(USE_MLP)
#include "mlp_ptr.h"
#endif

#if !defined(USE_MLP)
      ncall_r = ncall_r + 1
#if defined(MPI2)
      if (ncall_r == 1) then
        call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                           buffwin, ierror)
      endif
#endif
#endif

! Recv from North
      if ( jlast < jm ) then
         j = jlast
#if !defined(USE_MLP)
#if !defined(MPI2)
        nread = nread + 1
        call mpi_wait(rqest(nread), Status, ierror)
#endif
        tdisp = igonorth*idimsize + (ncall_r-1)*idimsize*nbuf
        call BufferUnPack3dx(p, 1, im, jfirst, jlast, kfirst, klast, &
                                1, im, j     , j    , kfirst, klast, &
                                buff_r(tdisp+1))
#else
!$omp parallel do private(i, j, k)
         do k=kfirst,klast
            do i=1,im
               p(i,j,k) = p(i,j,k) + g_4d(i,j,k,1)
            enddo
         enddo
#endif
      endif

#if !defined(USE_MLP)
      if (ncall_r == ncall_s) then
#if !defined(MPI2)
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nrecv = 0
        nread = 0
        nsend = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif
#endif
      end subroutine mp_send_ua_ad


      subroutine mpi_bcst_real_ad( buf_ad )
      implicit none
      real     buf_ad
      real     res_ad
      integer  myid, ierr

      call MPI_Reduce( buf_ad,res_ad,1,mp_precision,MPI_SUM,0,commglobal,ierr )
      if (gid .eq. 0) then
        buf_ad = res_ad
      else
        buf_ad = 0.
      endif

      end subroutine mpi_bcst_real_ad


      subroutine mpi_bcst_n_real_ad( buf_ad, n )
      implicit none
      integer  n
      real     buf_ad(n)
      real     res_ad(n)
      integer  myid, ierr

      call MPI_Reduce( buf_ad,res_ad,n,mp_precision,MPI_SUM,0,commglobal,ierr )
      if (gid .eq. 0) then
        buf_ad(:) = res_ad(:)
      else
        buf_ad(:) = 0.
      endif

      end subroutine mpi_bcst_n_real_ad


!-----
      subroutine mp_bcst_n_real(val,n)
!-----
! Send real "val" from Process=id to All other Processes
      real val(n)
      integer n

#if !defined(USE_MLP)
      call mpi_bcast(val, n, mp_precision, 0, commglobal, ierror)
#else
#include "mlp_ptr.h"
! FastOpt: missing variable gnn to hold val(n)
      if ( gid == 0 ) then
          gnn(1:n) = val(1:n)
      endif

      call mlp_barrier(gid, gsize)

      if ( gid /= 0 ) then
          val(1:n) = gnn(1:n)
      endif
      call mlp_barrier(gid, gsize)   !may not be necessary, BW
#endif
      end subroutine mp_bcst_n_real


!-----
      subroutine mp_scatter4d_ad(qout, qin, im, jm, km, nq, jfirst, jlast, &
                                 kfirst, klast, ng_s, ng_n, id)
!-----

      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast, kfirst, klast
      integer, intent(in):: ng_s, ng_n
      integer, intent(in):: nq
      integer, intent(in):: id         ! process ID to gather data to
      real, intent(inout) :: qin(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
      real, intent(inout) :: qout(im,jm,km,nq)

! Local:
      integer i, j, k, iq
      integer j1, j2, k1, k2
      integer n
      integer qsize
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP) && defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
 
#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4dwin, ierror)
#endif
      dest = id
      tdisp = ((jfirst-1)*km+(kfirst-1))*im*nq
      call BufferPack4dx(qin, 1, im, jfirst-ng_s, jlast+ng_n, &
                         kfirst, klast, 1, nq, 1, im, jfirst, jlast,&
                         kfirst, klast, 1, nq, buff4d(tdisp+1))
      qsize = (jlast-jfirst+1)*im*(klast-kfirst+1)*nq
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
      do p=1,pkgs_per_pro
        tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = tdisp + (p-1)*tmpsize
        call MPI_PUT(buff4d(mydisp+1), mysize, mp_precision, dest, &
                     mydisp, mysize, mp_precision, buff4dwin, ierror)
      enddo
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4dwin, ierror)
#else
      send_tag = gid
      nsend = nsend + 1
      call mpi_isend(buff4d(tdisp+1), qsize, mp_precision, dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
#endif
      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
!
! BMP: fvGCM edit
! Rough hack to account for gathering arrays k=1:1 or k=1:km+1
!
           if (km == PLEV) then
             k1 = zfirst(n)
             k2 = zlast(n)
           else
             k1 = 1
             k2 = km
           endif
           tdisp = ((j1-1)*km+(k1-1))*im*nq

#if !defined(MPI2)
           qsize = im*(j2-j1+1)*(k2-k1+1)*nq
           src = n-1
           recv_tag = src
           call mpi_recv(buff4d(tdisp+1), qsize, mp_precision, src, &
                         recv_tag, commglobal, Status, ierror)
#endif
           call BufferUnPack4dx(qout, 1, im, 1,  jm, 1, km, 1, nq, &
                                      1, im, j1, j2, k1, k2, 1, nq, &
                                      buff4d(tdisp+1))
         enddo
      endif

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, sqest, Stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
      do iq=1,nq
!$omp parallel do private(i,j,k)
         do k=kfirst,klast
            do j=jfirst,jlast
               do i=1,im
                  g_4d(i,j,k,iq) = qin(i,j,k,iq)
               enddo
            enddo
         enddo
      enddo
      call mp_barrier

      if ( gid == id ) then
      do iq=1,nq
!$omp parallel do private(i,j,k)
         do k=1,km
            do j=1,jm
               do i=1,im
                  qout(i,j,k,iq) = qout(i,j,k,iq) + g_4d(i,j,k,iq) ! FastOpt
               enddo
            enddo
         enddo
      enddo
      endif
      call mp_barrier
#endif
      end subroutine mp_scatter4d_ad


!-----
      subroutine mp_gather4d_ad(qout, qin, im, jm, km, nq, jfirst, jlast, &
                                 kfirst, klast, ng_s, ng_n, id)
!-----

      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast, kfirst, klast
      integer, intent(in):: ng_s, ng_n
      integer, intent(in):: nq
      integer, intent(in):: id           ! Source (usually 0)
      real, intent(inout):: qin(im,jm,km,nq)
      real, intent(inout):: qout(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)

! Local:
      integer i, j, k, n, iq
      integer j1, j2, k1, k2
      integer qsize_s
      integer qsize_r
      integer src, dest
      integer send_tag, recv_tag
#if !defined(USE_MLP)
      integer rqst(numpro), rq_stats(numpro*MPI_STATUS_SIZE)
#if defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#endif

#if !defined(USE_MLP)
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, buff4dwin, ierror)
#endif
      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
           k1 = zfirst(n)
           k2 = zlast(n)
           qsize_s = im*(j2-j1+1)*(k2-k1+1)*nq
           dest = n-1
           tdisp = ((j1-1)*km+(k1-1))*im*nq
           call BufferPack4dx(qin, 1, im, 1,  jm, 1, km, 1, nq, &
                                   1, im, j1, j2, k1, k2, 1, nq, &
                                   buff4d(tdisp+1))
#if defined(MPI2)
!$omp parallel do private(p,tmpsize,mysize,mydisp)
           do p=1,pkgs_per_pro
             tmpsize = ceiling(real(qsize_s)/real(pkgs_per_pro))
             mysize = MIN(tmpsize, MAX(qsize_s-(tmpsize*(p-1)),0))
             mydisp = tdisp + (p-1)*tmpsize
             call MPI_PUT(buff4d(mydisp+1), mysize, mp_precision, dest, &
                          mydisp, mysize, mp_precision, buff4dwin, ierror)
           enddo
#else
           send_tag = gid
           nsend = nsend + 1
           call mpi_isend(buff4d(tdisp+1), qsize_s, mp_precision, dest, &
                          send_tag, commglobal, rqst(nsend), ierror)
#endif
         enddo
      endif
      tdisp = ((jfirst-1)*km+(kfirst-1))*im*nq
#if defined(MPI2)
      call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                         buff4dwin, ierror)
#else
      qsize_r = im*(jlast-jfirst+1)*(klast-kfirst+1)*nq
      src = id
      recv_tag = src
      call mpi_recv(buff4d(tdisp+1), qsize_r, mp_precision, src, &
                    recv_tag, commglobal, Status, ierror)
#endif
      call BufferUnPack4dx(qout, 1, im, jfirst-ng_s, jlast+ng_n, kfirst,&
                           klast, 1, nq, 1, im, jfirst, jlast, kfirst,  &
                           klast, 1, nq, buff4d(tdisp+1))

#if !defined(MPI2)
      if (nsend /= 0) then
        call mpi_waitall(nsend, rqst, rq_stats, ierror)
        nsend = 0
      endif
#endif

#else
#include "mlp_ptr.h"
      if ( gid == id ) then
      do iq=1,nq
!$omp parallel do private(i,j,k)
         do k=1,km
            do j=1,jm
               do i=1,im
                  g_4d(i,j,k,iq) = qin(i,j,k,iq)
               enddo
            enddo
         enddo
      enddo
      endif
      call mp_barrier

      do iq=1,nq
!$omp parallel do private(i,j,k)
         do k=kfirst,klast
            do j=jfirst,jlast
               do i=1,im
                  qout(i,j,k,iq) = qout(i,j,k,iq) + g_4d(i,j,k,iq) ! FastOpt add
               enddo
            enddo
         enddo
      enddo
      call mp_barrier
#endif
      end subroutine mp_gather4d_ad


#endif /*  defined ( SPMD )  */

      end module mod_comm
