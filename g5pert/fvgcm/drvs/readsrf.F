!-----------------------------------------------------------------------
!BOP
! !ROUTINE: readsrf --- Read the surface data file (and distributed data)
!
! !INTERFACE:
      subroutine readsrf(iuhs, srfdata, im, jm, jfirst, jlast, phis, gz)

#if defined( SPMD )
      use mod_comm, only : gid, mp_barrier, mp_scatter2d
#endif

      use m_ioutil, only : luavail
      use m_dyn
      use m_inpak90

      implicit none

! !INPUT PARAMETERS:
      integer iuhs
      character*80 srfdata                 ! Name of surface data file
      integer im, jm                       ! Dimensions of arrays
      integer jfirst, jlast

! !OUTPUT PARAMETERS:
      real  phis(im,jfirst:jlast)        ! Topography
      real  gz(im,jfirst:jlast)   

! !DESCRIPTION:
!     Read the surface data file and distribute the data if necessary
!
! !REVISION HISTORY:
!   WS  99.09.07:  Creation
!   WS  00.07.10:  Cleaned up MPI usage
!   SJL 08/29/01:  Bug fix for MLP
!   Elena N. 06.14.06: Fix to read surface fields from rst.hdf file
!   Oloso    11.20.06: Added barrier for SPMD option
!   Todling  01.16.07  Merged Oloso/Elena's changes
!   Todling  03.13.09  generalized name of rst file
!   Todling  10.28.09  pass vectype to dyn_get (always geos-4 type)
!   Todling  01.19.10  temp arrays can be zero-dim except on root
!
!EOP
!---------------------------------------------------------------------
!BOC
      logical exist
      real,allocatable::  g_phis(:,:)
      real,allocatable::  g_gz(:,:)   

      character(len=*),  parameter :: rcname = 'ccmrun.namelist'
      character(len=80) token

      integer ::  iuhs2
      integer ::  rc
      integer ::  nymd, nhms 
      character(len=80) :: rstfname
      type(dyn_vect) :: w_f     ! fv dynamics vector in eta

      nymd = 0
      nhms = 0

#if defined( SPMD )
      call mp_barrier
      if ( gid == 0 ) then
         allocate(g_phis(im,jm),g_gz(im,jm))
#endif

#if defined (BINRST)
        inquire (file=srfdata, exist=exist)
        if (.not. exist) then
          print *, 'GET_SRF: Cannot find surface data file ---> ',
     &             srfdata
          print *, 'Stop.'
          stop
        endif
        open (iuhs, file=srfdata, form='unformatted', status='old')

#if defined( SPMD )
        read(iuhs) g_phis
        read(iuhs) g_gz
#else
        read(iuhs) phis
        read(iuhs) gz
#endif
        close (iuhs)

#else

        inquire( file=rcname, exist=exist )
        if ( exist ) then  
           call i90_loadf (rcname, rc)
           if( rc .eq. 0) then
              call I90_label('g4_restart_name:', rc)
              if (rc .eq. 0) then
                  call I90_Gtoken(token,rc)
                  if( rc .eq. 0) then
                      rstfname = trim(token) 
                  endif
              endif 
           endif 
           call i90_release ()
        endif

        call dyn_null(w_f)
        call dyn_get ( trim(rstfname), nymd, nhms, w_f, rc, vectype=4 )

        if ((w_f%grid%im .ne. im) .or. (w_f%grid%jm .ne. jm)) then 
           print *, ' Resolution of Dynamic restart vector (2D fields): ', w_f%grid%im, w_f%grid%jm
           print *, ' Required resolution: ', im, jm
           print *, ' Exiting on READRST error'
           stop
        endif

        if (rc .ne. 0) then
               print *, 'could not read dynamics state file to get restarts'
               stop
        endif

#if defined( SPMD )
        g_phis          = w_f%phis
        g_gz            = w_f%hs_stdv
#else
        phis            = w_f%phis
        gz              = w_f%hs_stdv
#endif

        call dyn_clean ( w_f )

        print *
        print *,'USING RST_NC TO READ HS and HS_STDV'
        print *

#endif


#if defined( SPMD )
      else  ! <gid=0>
         allocate(g_phis(0,0),g_gz(0,0))
      endif ! <gid=0>
      call mp_barrier
      call mp_scatter2d(g_phis, phis, im, jm, jfirst, jlast, 0)
      call mp_scatter2d(g_gz,   gz,   im, jm, jfirst, jlast, 0)

      deallocate(g_phis,g_gz)
#endif

      return
!EOC
      end subroutine readsrf
!---------------------------------------------------------------------

