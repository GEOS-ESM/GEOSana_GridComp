      module m_pmgrid

      integer, save :: begj       ! actual computation begins
      integer, save :: endj       ! computation ends
 
      integer, save :: plevmx  = 4
      integer, save :: plevp
      integer, save :: plond
       
      integer, save :: beglat
      integer, save :: endlat

      contains
      subroutine pmgrid_set ( plon, plat, plev, jfirst, jlast ) 
      implicit none
      integer, intent(in) :: plon, plat, plev, jfirst, jlast
      beglat = 1
      endlat = PLAT
      plevp  = PLEV + 1
      plond  = PLON  ! per SJL
      begj = jfirst
      endj = jlast
      end subroutine pmgrid_set
      end module m_pmgrid
