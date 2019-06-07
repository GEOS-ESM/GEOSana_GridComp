#include <preproc.h>
      subroutine twostr (npt     ,
     &                   fab    ,fre    ,ftd   ,fti    ,albgrd ,
     &                   albgri ,coszen ,ib    ,ic     ,nsol   ,
     &                   indsol ,ivt    ,vai   ,fwet   ,rho    ,
     &                   tau    ,t      ,tfrz  ,omegas ,betads ,
     &                   betais ,xl     ,gdir  )
 
#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       twostr.F
* purpose:           two-stream fluxes for canopy radiative transfer
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------
 
* ------------------------ notes ----------------------------------
* use two-stream approximation of Dickinson (1983) Adv Geophysics
* 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
* to calculate fluxes absorbed by vegetation, reflected by vegetation,
* and transmitted through vegetation for unit incoming direct or diffuse 
* flux given an underlying surface with known albedo.
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer npt             !number of points
      integer ib              !waveband number 
      integer ic              !0=unit incoming direct; 1=unit incoming diffuse
      integer nsol            !number of points with coszen > 0

      real tfrz               !freezing temperature (kelvin)
      real betads             !betad for snow
      real betais             !betai for snow

      integer indsol(npt)     !index of points with coszen > 0
      integer ivt(npt)        !vegetation type

      real coszen(npt)        !cosine of direct zenith angle
      real vai(npt)           !one-sided leaf+stem area index
      real fwet(npt)          !fraction of lai, sai that is wetted
      real t(npt)             !surface temperature (kelvin)

      real albgrd(mband,npt)  !direct  albedo of underlying surface 
      real albgri(mband,npt)  !diffuse albedo of underlying surface
      real rho(mband,npt)     !leaf+stem reflectance
      real tau(mband,npt)     !leaf+stem transmittance

      real omegas(mband)      !omega for snow

      real xl(mvt)            !departure of leaf angles from spherical distrib

* output
      real fab(mband,npt)  !flux abs by veg layer (per unit incoming flux)
      real fre(mband,npt)  !flux refl above veg layer (per unit incoming flux)
      real ftd(mband,npt)  !down dir flux below veg layer (per unit in flux)
      real fti(mband,npt)  !down dif flux below veg layer (per unit in flux)
      real gdir(npt)       !relative projected leaf+stem area in solar direction
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      real omega(npt)   !fraction of intercepted radiation that is scattered
      real omegal(npt)  !omega for leaves
      real betai(npt)   !upscatter parameter for diffuse radiation 
      real betail(npt)  !betai for leaves
      real betad(npt)   !upscatter parameter for direct beam radiation 
      real betadl(npt)  !betad for leaves
      real ext(npt)     !optical depth of direct beam per unit leaf area 
      real avmu(npt)    !average diffuse optical depth

      integer j,i       !array index 

      real cosz         !0.001 <= coszen <= 1.000
      real asu          !single scattering albedo
      real chil         ! -0.4 <= xl <= 0.6

      real tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
      real p1,p2,p3,p4,s1,s2,u1,u2,u3
      real b,c,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
      real phi1,phi2,sigma
      real ftds,ftis,fres
* -----------------------------------------------------------------

      if (nsol .eq. 0) return


* calculate two-stream parameters omega, betad, betai, avmu, gdir, ext.
* omega, betad, betai are adjusted for snow. values for omega*betad 
* and omega*betai are calculated and then divided by the new omega
* because the product omega*betai, omega*betad is used in solution. 
* also, the transmittances and reflectances (tau, rho) are linear 
* weights of leaf and stem values.

      do j = 1, nsol
         i = indsol(j)
         cosz = max(0.001, coszen(i))
         chil = min( max(xl(ivt(i)), -0.4), 0.6)
         if (abs(chil) .le. 0.01) chil = 0.01
         phi1 = 0.5 - 0.633*chil - 0.330*chil*chil
         phi2 = 0.877 * (1.-2.*phi1)
         gdir(i) = phi1 + phi2*cosz
         ext(i) = gdir(i)/cosz
         avmu(i) = ( 1. - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
         omegal(i) = rho(ib,i) + tau(ib,i)
         tmp0 = gdir(i) + phi2*cosz
         tmp1 = phi1*cosz
         asu = 0.5*omegal(i)*gdir(i)/tmp0 * ( 1. - tmp1/tmp0 *
     &         log((tmp1+tmp0)/tmp1) )
         betadl(i) = (1.+avmu(i)*ext(i))/(omegal(i)*avmu(i)*ext(i))*asu
         betail(i) = 0.5 * ( rho(ib,i)+tau(ib,i) + (rho(ib,i)-tau(ib,i)) 
     &               * ((1.+chil)/2.)**2 ) / omegal(i)
      end do

* adjust omega, betad, and betai for intercepted snow

      do j = 1, nsol
         i = indsol(j)
         if (t(i) .gt. tfrz) then                                !no snow
            tmp0 = omegal(i)           
            tmp1 = betadl(i) 
            tmp2 = betail(i)  
         else
            tmp0 = (1.-fwet(i))*omegal(i) + fwet(i)*omegas(ib)           
            tmp1 = ( (1.-fwet(i))*omegal(i)*betadl(i) +    
     &              fwet(i)*omegas(ib)*betads ) / tmp0
            tmp2 = ( (1.-fwet(i))*omegal(i)*betail(i) +  
     &              fwet(i)*omegas(ib)*betais ) / tmp0
         end if
         omega(i) = tmp0           
         betad(i) = tmp1 
         betai(i) = tmp2  
      end do

* absorbed, reflected, transmitted fluxes per unit incoming radiation

      do j = 1, nsol
         i = indsol(j)
         b = 1. - omega(i) + omega(i)*betai(i)
         c = omega(i)*betai(i)
         tmp0 = avmu(i)*ext(i)
         d = tmp0 * omega(i)*betad(i)
         f = tmp0 * omega(i)*(1.-betad(i))
         tmp1 = b*b - c*c
         h = sqrt(tmp1) / avmu(i)
         sigma = tmp0*tmp0 - tmp1
         p1 = b + avmu(i)*h
         p2 = b - avmu(i)*h
         p3 = b + tmp0
         p4 = b - tmp0
         s1 = exp(-h*vai(i))
         s2 = exp(-ext(i)*vai(i))
         if (ic .eq. 0) then
            u1 = b - c/albgrd(ib,i)
            u2 = b - c*albgrd(ib,i)
            u3 = f + c*albgrd(ib,i)
         else
            u1 = b - c/albgri(ib,i)
            u2 = b - c*albgri(ib,i)
            u3 = f + c*albgri(ib,i)
         end if
         tmp2 = u1 - avmu(i)*h
         tmp3 = u1 + avmu(i)*h
         d1 = p1*tmp2/s1 - p2*tmp3*s1
         tmp4 = u2 + avmu(i)*h
         tmp5 = u2 - avmu(i)*h
         d2 = tmp4/s1 - tmp5*s1
         h1 = -d*p4 - c*f
         tmp6 = d - h1*p3/sigma
         tmp7 = ( d - c - h1/sigma*(u1+tmp0) ) * s2
         h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
         h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
         h4 = -f*p3 - c*d
         tmp8 = h4/sigma
         tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
         h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
         h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
         h7 = (c*tmp2) / (d1*s1)
         h8 = (-c*tmp3*s1) / d1
         h9 = tmp4 / (d2*s1)
         h10 = (-tmp5*s1) / d2

* downward direct and diffuse fluxes below vegetation
 
         if (ic .eq. 0) then
            ftds = s2
            ftis = h4*s2/sigma + h5*s1 + h6/s1
         else
            ftds = 0.
            ftis = h9*s1 + h10/s1
         end if
         ftd(ib,i) = ftds
         fti(ib,i) = ftis

* flux reflected by vegetation
 
         if (ic .eq. 0) then
            fres = h1/sigma + h2 + h3
         else
            fres = h7 + h8
         end if
         fre(ib,i) = fres

* flux absorbed by vegetation

         fab(ib,i) = 1. - fre(ib,i) - (1.-albgrd(ib,i))*ftd(ib,i) 
     &               - (1.-albgri(ib,i))*fti(ib,i)
 
      end do

      return
      end
 
