      integer plon       ! number of longitudes
      integer plev       ! number of vertical levels
      integer plat       ! number of latitudes
      integer pcnst      ! number of constituents (including water vapor)
      integer pnats      ! number of non-advected trace species
      integer plevmx     ! number of subsurface levels
      integer plevp      ! plev + 1
      integer plond      ! slt extended domain longitude
      integer plevd      ! fold plev,pcnst indices into one

      integer beglat     ! beg. index for 2D phys arrays
      integer endlat     ! end. index for 2D phys arrays

      integer begj       ! actual computation begins
      integer endj       ! computation ends
 
      logical masterproc ! Flag for (iam eq 0)
 
      parameter (plon    = PLON) 
      parameter (plev    = PLEVR)
      parameter (plat    = PLAT)
      parameter (pcnst   = PCNST)
      parameter (pnats   = PNATS)
      parameter (plevmx  = 4)
      parameter (plevp   = PLEV + 1)
      parameter (plond   = PLON)    ! SJL
       
      parameter ( beglat = 1 )
      parameter ( endlat = PLAT )

      common/spmdlats/begj, endj
      common/spmdlats/masterproc
   
      
