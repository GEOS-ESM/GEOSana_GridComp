#include <preproc.h>
      subroutine soialb (npt    ,
     &                   nsol   ,nband  ,istice ,ist    ,isc    ,
     &                   indsol ,fsno   ,h2osoi ,albsnd ,albsni ,
     &                   albsat ,albdry ,albice ,albgrd ,albgri ,
     &                   tg     ,tfrz   ,alblak ,coszen )
 
#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       soialb.F
* purpose:           ground surface albedo, accounting for snow
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------
 
* ------------------------- input/output variables ----------------
* input
      integer npt              !number of points
      integer nband            !number of solar radiation waveband classes
      integer nsol             !number of points with coszen > 0
      integer istice           !value of ist for land ice
      real tfrz                !freezing temperature (kelvin)

      integer ist(npt)         !"soil texture" class
      integer isc(npt)         !soil color class
      integer indsol(npt)      !index for points with coszen > 0

      real fsno(npt)           !fraction of surface covered with snow
      real h2osoi(msl,npt)     !volumetric soil water content
      real tg(npt)             !ground temperature (kelvin)
      real coszen(npt)         !cosine solar zenith angle

      real albsnd(mband,npt)   !direct beam snow albedo (vis, nir)
      real albsni(mband,npt)   !diffuse snow albedo (vis, nir)

      real albsat(msc,mband)   !saturated soil albedo (vis, nir)
      real albdry(msc,mband)   !dry soil albedo (vis, nir)

      real albice(mband)       !albedo of land ice (vis, nir)
      real alblak(mband)       !albedo of frozen water (vis, nir)

* output
      real albgrd(mband,npt)   !ground albedo (direct beam: vis, nir)
      real albgri(mband,npt)   !ground albedo (diffuse: vis, nir)
* -----------------------------------------------------------------

* ------------------------- local variables -----------------------
      integer i,j              !loop/array index
      integer ib               !waveband number (1=vis, 2=nir)

      real inc                 !soil water correction factor for soil albedo
      real albsod              !soil albedo (direct)
      real albsoi              !soil albedo (diffuse)
* -----------------------------------------------------------------

      do ib = 1, nband
         do j = 1, nsol
            i = indsol(j)
            inc = max(0.11-0.40*h2osoi(1,i), 0.)
            if (ist(i) .eq. 1)  then                     !soil
               albsod = min(albsat(isc(i),ib)+inc,albdry(isc(i),ib))
               albsoi = albsod
            else if (ist(i) .eq. istice)  then           !land ice
               albsod = albice(ib)
               albsoi = albsod
            else if (tg(i) .gt. tfrz) then               !unfrozen lake, wetland
               albsod = 0.06/(max(0.01,coszen(i))**1.7 + 0.15)
               albsoi = 0.06
            else                                         !frozen lake, wetland
               albsod = alblak(ib)
               albsoi = albsod
            end if

* increase desert and semi-desert albedos

            if (ist(i) .eq. 1 .and. isc(i) .eq. 9) then
               albsod = albsod + 0.10
               albsoi = albsoi + 0.10
            end if

            albgrd(ib,i) = albsod*(1.-fsno(i)) + albsnd(ib,i)*fsno(i)
            albgri(ib,i) = albsoi*(1.-fsno(i)) + albsni(ib,i)*fsno(i)
         end do
      end do

      return
      end
 
