#include <params.h>
      subroutine sstan(lat     ,oro     ,ts      ,tssub   ,snowh   ,
     $                 sicthk  )
C-----------------------------------------------------------------------
C
C Update sea surface temperatures (sst's) and sea ice distribution
C
C Assume that the sst data exists in a two dimensional field
C encoded as follows:
C
C     Land               values where oro field says so ("valid sst's
C                        are provided globally, the model's land mask
C                        determines whether the sst is used or not)
C
C     Ocean without      values degrees celcius (greater than tsice)
C      sea ice
C
C     Ocean with         values less than tsice
C      sea ice
C
C New sea ice has a constant 0.5 cm value for snow cover prescribed
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Buja, Jun 1992, Feb 1996
C Reviewed:          J. Hack, B. Boville, Aug 1992, Feb 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
#include <parsst.h>
C------------------------------Commons----------------------------------
#include <comsst.h>
C-----------------------------------------------------------------------
#include <comtsc.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat               ! Latitude index
C
C input/output arguments
C
      real oro(plond)           ! Surface type flag array
      real ts(plond)            ! Surface temperature
      real tssub(plond,plevmx)  ! Sub/srf tmps over ocean set to sst.
      real snowh(plond)         ! Snow cover (liquid water equivalent)
      real sicthk(plond)        ! Sea-ice thickness (m)
C
C---------------------------Local variables-----------------------------
C
      integer i,ii              ! Longitude indices
      integer n                 ! Sub/srf level index
      integer nvalsi            ! Number of sea-ice points
      integer nvalo             ! Number of open ocean points
      integer indxsi(plond)     ! Indices of sea-ice points
      integer indxo(plond)      ! Indices of open ocean points
C
C--------------------------Statement Functions--------------------------
C
      logical land,ocean
      land(i) = nint(oro(i)).eq.1
      ocean(i) = nint(oro(i)).eq.0
C
C-----------------------------------------------------------------------
C
C Find where dataset says sea ice.
C
      call whenfle(plon,sst(1,lat),1,tsice,indxsi,nvalsi)
C
C       Change things only where dataset *and* model say not land
C
      do ii=1,nvalsi
        i=indxsi(ii)
        if (.not.land(i)) then
          snowh(i) = 0.005
          sicthk(i) = 2.0
        end if
      end do
C
C New sea ice.  Set oro field and initialize surface and subsurface
C temperatures.
C
      do ii=1,nvalsi
        i = indxsi(ii)
        if (ocean(i)) then
          oro(i) = 2.0
          ts(i) = sst(i,lat) + 273.16
          do n=1,plevmx
            tssub(i,n) = sst(i,lat) + 273.16
          end do
        end if
      end do
C
C Open ocean
C
      nvalo = 0
      do i=1,plon
        if (.not.land(i) .and. sst(i,lat).gt.tsice) then
          nvalo = nvalo + 1
          indxo(nvalo) = i
          ts(i) = sst(i,lat) + 273.16
          oro(i) = 0.
          snowh(i) = 0.
          sicthk(i) = 0.
        end if
      end do
      do n=1,plevmx
        do ii=1,nvalo
          i=indxo(ii)
          tssub(i,n) = sst(i,lat) + 273.16
        end do
      end do
C
      return
      end
 
