      subroutine binf2c(flon   , flat  ,nflon   ,nflat   ,fine   ,
     $                  clon   , clat  ,nclon   ,nclat   ,cmean  )

c-----------------------------------------------------------------------------
c       Bin going from a fine grid to a coarse grid.
c       A schematic for the coarse and fine grid systems is shown in
c       Figure 1.  This code assumes that each data point is represent
c       it's surrounding area, called a cell.  The first grid data point
c       for both grids is assumed to be located at 0E (GM).  This
c       implies that the 1st cell for both the fine and the coarse grids
c       strattles the Greenwich Meridian (GM).  This code also assumes
c       that there is no data wraparound (last data value is located at
c       360-dx).
c
c       FIGURE 1:  Overview of the coarse (X) and fine (@) grids
c                  longitudinal structure where:
c                  X = location of each coarse grid data point
c                  @ = location of each fine   grid data point
c
c           Greenwich                                     Greenwich
c              0          Coarse cells                       360
c              :          v                                   :
c       clon(1):  clon(2) v clon(3)                clon(nclon):
c         v    :    v     v   v                        v      :
c         xxxxxxxxxxxxxxxxxxxxxxxxxxxx..xxxxxxxxxxxxxxxx      :
c         x         x         x              x         x      :
c         x         x         x              x         x      :
c         x  c(1)   x  c(2)   x              x c(nclon)x      :
c         x    X    x    X    x              x    X    x      :
c         x   ___ ___ ___ ___ ___ ___    ___ ___ ___ ___ ___  :
c         x  |   |   |   |   |   |   |  |   |   |   |   |   | :
c         x  | @ | @ | @ | @ | @ | @ |..| @ | @ | @ | @ | @ | :
c         xxx|___|___|___|___|___|___|  |___|___|___|___|___| :
c              v       v       v                      v   v   :
c          flon(1)   flon(3)   v            flon(nflon-1) flon(nflon)
c              :               v                              :
c              :               Fine cells                     :
c              0                                             360
c
c    The Longitude/Latitude search:
c    ------------------------------
c
c      Given a coarse grid cell with west and east boundaries of cWest
c      and cEast and south and north boundaries of cSouth and cNorth
c      (outlined by "x" in figure 2), find the indices of the fine grid
c      points which are contained within the coarse grid cell. imin and 
c      imax are the indices fine grid points which overlap the western
c      and eastern boundary of the coarse cell.  jmin and jmax are the
c      corresponding indices in the S-N direction.  Bin these overlapping
c      values to generate coarse(n), the coarse grid data values.
c
c        FIGURE 2: Detail of Coarse and Fine cell overlap.
c                   @ = fine   grid data point
c                   X = coarse grid data point
c
c                             cWest             cEast
c                   |       |   x   |       |     x |
c                  -@-------@---x---@-------@-----x-@-
c                   |       |  x*xxxxxxxxxxxxxxxxx*x|xx cNorth
c                   |       |   x   |       |     x |
c                   |       |   x   |       |     x |
c                   @-------@---x---@-------@-----x-@- jmax
c                   |       |   x   |  c(n) |     x |
c                   |   @   |   x   |       |     x |
c                   |       |   x   |       |     x |
c                   @-------@---x---@-------@-----x-@- jmin
c                   |       |   x   |       |     x |
c                   |   @   |  x*xxxxxxx@xxxxxxxxx*x|xx cSouth
c                   |       |   x   |       |     x |
c                  -@-------@---x---@-------@-----x-@-  
c                   |             imin    imax      |
c
c
c      When a cell coarse cell strattles the Greenwich Meridian
c      ---------------------------------------------------------
c
c      The first coarse grid cell strattles the GM, so when the western
c      boundary of the coarse cell is < 0, an additional search is carried out. 
c      It ASSUMES that the easternmost fine grid point overlaps and searches
c      westward from nflon, looking for a grid point west of clon(1)
c      This generates a second set of longitudinal indices, imin1 and imax1.   
c      See Figure 3.
c
c       Figure 3:  Detail of Coarse cell strattling GM:
c       -----------------------------------------------
c
c              Greenwich                                   Greenwich
c                 0                                           360
c          cWest  :  cEast                            cWest    :
c          clon(1):  clon(2)               clon(nclon+1)=clon(1)
c            v    :    v                                v      :
c            xxxxxxxxxxxxxxxxxxxxxxx ... xxxxxxxxxxxxxxxx      :
c            x         x        x            x          x      :
c            x         x        x            x          x      :
c            x  c(1)   x        x            x  c(nclon)x      :
c            x    X    x        x            x     X    x      :
c            x   ___ ___ ___ _                    ___ ___ ___  :
c            x  |   |   |   |                        |   |   | :
c            x  | @ | @ | @ |                      @ | @ | @ | :
c            xxx|___|___|___|_                    ___|___|___| :
c                 ^ : ^   ^                            ^   ^   :
c            flon(1): ^ flon(3)              flon(nflon-1) ^   :
c                 ^ : ^                                ^   ^   :
c                 ^ :flon(2)                           ^ flon(nflon)
c                 ^ : ^                                ^   ^   :
c              imin : imax                         imin1 imax1 :
c                   :                                          :
c
c
c         In this case, imin=1, imax=2, imin1=nflon-1 and imax1=nflon.
c         because the last two cells of the fine grid will have some
c         contribution the the 1st cell of the coarse grid.
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------Arguments---------------------------------

      integer nflon              ! Input: number of fine longitude points
      integer nflat              ! Input: number of fine latitude points
      real    flon(nflon)        ! Input: fine grid lons, centers (deg)
      real    flat(nflat)        ! Input: fine grid lats, centers (deg)
      real    fine(nflon,nflat)  ! Input: Fine grid data array

      integer nclon              ! Input: number of coarse longitude points
      integer nclat              ! Input: number of coarse latitude points
      real    clon(nclon+1)      ! Input: coarse grid cell lons, west  edge (deg)
      real    clat(nclat+1)      ! Input: coarse grid cell lat,  south edge (deg)
      real    cmean(nclon,nclat) ! Output: mean of fine grid points over coarse cell

c--------------------------Local variables------------------------------

      real    cWest              ! Coarse cell longitude, west edge (deg)
      real    cEast              ! Coarse cell longitude, east edge (deg)
      real    cSouth             ! Coarse cell latitude, south edge (deg)
      real    cNorth             ! Coarse cell latitude, notrh edge (deg)
      real    sum               ! coarse tmp value
      
      integer i,j                ! Indices
      integer imin ,imax         ! Max/Min E-W indices of intersecting fine cell.
      integer imin1,imax1        ! Additional fine E-W indices when coarse cell strattles GM
      integer jmin ,jmax         ! Max/Min N-S indices of intersecting fine cell.
      integer iclon,jclat        ! coarse grid indices
      integer num                ! increment 

c-----------------------------------------------------------------------------

      do 200 jclat= 1,nclat         ! loop over coarse latitudes
         cSouth = clat(jclat)
         cNorth = clat(jclat+1)

         do 100 iclon=1,nclon       ! loop over coarse longitudes
            cWest  = clon(iclon)
            cEAST  = clon(iclon+1)

c  1. Normal longitude search:  Find imin and imax

            imin = 0
            imax = 0
            do i=1,nflon-1                     ! loop over fine lons, W -> E
               if (flon(i) .gt. cEast) goto 10 ! fine grid point is E of coarse box
               if (flon(i) .ge. cWest .and. imin.eq.0) imin=i
               imax=i
            enddo

c  2. If cWest < 0, then coarse cell strattles GM.  Hunt westward
c     from the end to find indices of any overlapping fine grid cells:
c     imin1 and imax1.

   10       imin1 = 0           ! borders for cWest, cEast
            imax1 = -1          ! borders for cWest, cEast
            if (cWest .lt. 0) then
               cWest = cWest + 360.
               imax1 = nflon
               do i=nflon,1,-1                     ! loop over fine lons, E -> W
                  imin1=i
                  if (flon(i) .le. cWest) goto 20 ! fine grid point is W of coarse box 
               enddo
            endif

c  3. Do the latitude search S -> N for jmin and jmax

   20       jmin = 0
            jmax = 0
            do  j=1,nflat                       ! loop over fine lats, S -> N
               if (flat(j) .gt. cNorth) goto 30 ! fine grid point is N of coarse box
               if (flat(j) .ge. cSouth .and. jmin.eq.0) jmin=j
               jmax=j
            enddo
   30       continue

c  4. Sum 

          sum = 0.                        ! Initialize coarse data value
          num  = 0

          do j=jmin,jmax                   ! loop over fine lats, S -> N
             do i=imin,imax                ! loop over fine lons, W -> E
                sum = sum + fine(i,j) 
                num  = num + 1
             enddo
             do i=imin1,imax1              ! If coarse cell strattles GM
                sum = sum + fine(i,j) 
                num  = num + 1 
             enddo
          enddo

          if (num .gt. 0) then
             cmean(iclon,jclat) = sum/num
          else
             cmean(iclon,jclat) = 1.e30
          endif

  100   continue
  200 continue
      
      return
      end

