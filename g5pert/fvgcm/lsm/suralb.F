#include <preproc.h>
      subroutine suralb (npt      ,ivt      ,ist      ,isc    ,
     &                   lati     ,long     ,calday   ,soot   ,
     &                   elai     ,esai     ,tam      ,tg     ,
     &                   tv       ,h2osno   ,h2osoi   ,fwet   ,
     &                   albgrd   ,albgri   ,albd     ,albi   ,
     &                   fabd     ,fabi     ,ftdd     ,ftid   ,
     &                   ftii     ,fsun     ,fsno     ,kpti   ,
     &                   loopi    ,nstep    )
 
#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ----------------------------------
* source file:       suralb.F
* purpose:           surface albedo and two-stream fluxes
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* ------------------------------------------------------------------------

* ------------------------ notes -----------------------------------------
* surface albedos. also fluxes (per unit incoming direct and diffuse
* radiation) reflected, transmitted, and absorbed by vegetation. 
* also sunlit fraction of the canopy. 

* the calling sequence is:
*   -> suralb:           albedos for next time step
*        -> lsmzen:      cosine solar zenith angle for next time step
*        -> snoalb:      snow albedos: direct beam
*        -> snoalb:      snow albedos: diffuse
*        -> soialb:      soil/lake albedos
*        -> twostr:      absorbed, reflected, transmitted solar fluxes (vis dir)
*        -> twostr:      absorbed, reflected, transmitted solar fluxes (vis dif)
*        -> twostr:      absorbed, reflected, transmitted solar fluxes (nir dir)
*        -> twostr:      absorbed, reflected, transmitted solar fluxes (nir dif)
* ------------------------------------------------------------------------

* ------------------------ input/output variables ------------------------
* input
      integer npt           !number of points
      integer nstep         !time step index
      integer loopi         !"little" vector index (1 -> numlv)
      integer kpti          !first point in "big" kpt vec for loopi "little" vec
      integer ivt(npt)      !vegetation type
      integer ist(npt)      !soil "texture" type
      integer isc(npt)      !soil color type

      real calday           !calendar day at greenwich (1.00, ..., 365.99)

      real h2osno(npt)      !liquid water snow (kg/m**2)
      real h2osoi(msl,npt)  !volumetric soil water
      real soot(npt)        !soot concentration
      real tam(npt)         !"aerodynamic" surface air temperature (kelvin)
      real tg(npt)          !ground temperature (kelvin)
      real tv(npt)          !vegetation temperature (kelvin)
      real lati(npt)        !latitude  (+ = NH, - = SH)
      real long(npt)        !longitude (+ = EH, - = WH
      real elai(npt)        !leaf area, one-sided, adjusted for burying by snow
      real esai(npt)        !stem area, one-sided, adjusted for burying by snow
      real fsno(npt)        !fraction of grid covered by snow
      real fwet(npt)        !fraction of canopy that is wet

* output
      real albgrd(mband,npt)  !ground albedo (direct)
      real albgri(mband,npt)  !ground albedo (diffuse)
      real albd(mband,npt)    !surface albedo (direct)
      real albi(mband,npt)    !surface albedo (diffuse)
      real fabd(mband,npt)    !flux abs by veg (per unit direct flux) 
      real fabi(mband,npt)    !flux abs by veg (per unit diffuse flux) 
      real ftdd(mband,npt)    !down direct flux below veg (per unit dir flux) 
      real ftid(mband,npt)    !down diffuse flux below veg (per unit dir flux) 
      real ftii(mband,npt)    !down diffuse flux below veg (per unit dif flux) 

      real fsun(npt)          !sunlit fraction of canopy
* ------------------------------------------------------------------------

* ------------------------ common block variables ------------------------
#include <phycon.h>
#include <radcon.h>
#include <soicon.h>
#include <vegcon.h>
* ------------------------------------------------------------------------

* ------------------------ local variables -------------------------------
      integer i,k,ik           !indices
      integer nband            !number of solar radiation wave bands
      integer ic               !direct beam: ic=0; diffuse: ic=1
      integer nsol             !number of grid points with coszen > 0
      integer indsol(npt)      !index of grid points with coszen > 0

      real wl                  !fraction of lai+sai that is lai
      real ws                  !fraction of lai+sai that is sai
      real mpe                 !prevents overflow for division by zero 
 
      real rho(mband,npt)      !leaf/stem refl weighted by fraction lai and sai
      real tau(mband,npt)      !leaf/stem tran weighted by fraction lai and sai
      real ftdi(mband,npt)     !down direct flux below veg per unit dif flux = 0
      real albsnd(mband,npt)   !snow albedo (direct)
      real albsni(mband,npt)   !snow albedo (diffuse)

      real vai(npt)         !elai+esai
      real gdir(npt)        !average projected leaf/stem area in solar direction
      real ext (npt)        !optical depth direct beam per unit leaf+stem area
      real coszen (npt)     !cosine solar zenith angle for next time step
* ------------------------------------------------------------------------

      if (npt .eq. 0) return

      nband = mband
      mpe = 1.e-06

* cosine solar zenith angle for next time step

      call lsmzen (npt, calday, long, lati, coszen)

* initialize output because solar radiation only done if coszen > 0

      do i = 1, nband
         do k = 1, npt
            albd(i,k) = 0.
            albi(i,k) = 0.
            albgrd(i,k) = 0.
            albgri(i,k) = 0.
            fabd(i,k) = 0.
            fabi(i,k) = 0.
            ftdd(i,k) = 0.
            ftid(i,k) = 0.
            ftii(i,k) = 0.
            if (i.eq.1) fsun(k) = 0.
         end do
      end do

* index points with positive coszen for subsequent calculations 
 
      nsol = 0                                            
      do k = 1, npt
         if (coszen(k) .gt. 0.) then                       
            nsol = nsol + 1                              
            indsol(nsol) = k
         endif                                        
      end do
      if (nsol .eq. 0) then
         return
      end if
* weight reflectance/transmittance by lai and sai
      
      do i = 1, nband
         do ik = 1, nsol
            k = indsol(ik)
            vai(k) = elai(k) + esai(k)
            wl = elai(k) / max(vai(k),mpe)
            ws = esai(k) / max(vai(k),mpe)
            rho(i,k) = max( rhol(ivt(k),i)*wl+rhos(ivt(k),i)*ws, mpe)
            tau(i,k) = max( taul(ivt(k),i)*wl+taus(ivt(k),i)*ws, mpe)
         end do
      end do

* snow albedos: only if coszen > 0 and h2osno > 0

         ic=0
         call snoalb (npt, nband, h2osno, tam, coszen, ic, soot, albsnd)
         ic=1
         call snoalb (npt, nband, h2osno, tam, coszen, ic, soot, albsni)

* ground surface albedos

      call soialb (npt    ,
     &             nsol   ,nband  ,istice ,ist    ,isc    ,
     &             indsol ,fsno   ,h2osoi ,albsnd ,albsni ,
     &             albsat ,albdry ,albice ,albgrd ,albgri ,
     &             tg     ,tfrz   ,alblak ,coszen )

* loop over nband wavebands to calculate surface albedos and solar 
* fluxes for unit incoming direct (ic=0) and diffuse flux (ic=1)

      do i = 1, nband
         ic = 0
         call twostr (npt    ,
     &                fabd   ,albd   ,ftdd ,ftid   ,albgrd ,
     &                albgri ,coszen ,i    ,ic     ,nsol   ,
     &                indsol ,ivt    ,vai  ,fwet   ,rho    ,
     &                tau    ,tv     ,tfrz ,omegas ,betads ,
     &                betais ,xl     ,gdir )
         ic = 1
         call twostr (npt    ,
     &                fabi   ,albi   ,ftdi ,ftii   ,albgrd ,
     &                albgri ,coszen ,i    ,ic     ,nsol   ,
     &                indsol ,ivt    ,vai  ,fwet   ,rho    ,
     &                tau    ,tv     ,tfrz ,omegas ,betads ,
     &                betais ,xl     ,gdir )
      end do

* sunlit fraction of canopy. set fsun = 0 if fsun < 0.01.
* code looks strange but vectorizes

      do ik = 1, nsol
         k = indsol(ik)
         ext(k) = gdir(k)/coszen(k) * sqrt(1.-rho(1,k)-tau(1,k))
         fsun(k) = (1.-exp(-ext(k)*vai(k))) / max(ext(k)*vai(k),mpe)
         ext(k) = fsun(k)                                        !temporary fsun
      end do
      do ik = 1, nsol
         k = indsol(ik)
         if (ext(k) .lt. 0.01) then 
            wl = 0.                                              !temporary fsun
         else
            wl = ext(k)                                          !temporary fsun
         end if
         fsun(k) = wl
      end do

      return
      end
 
