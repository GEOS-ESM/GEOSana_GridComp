module SORADMOD

IMPLICIT NONE

PRIVATE
PUBLIC :: sorad, deledd, getvistau1, getnirtau1

contains

 subroutine sorad ( m               , & !Number of soundings
                    np              , & !Number of model levels
                    nb              , & !Number of bands
                    cosz_dev        , & !Cosine of solar zenith angle
                    pl_dev          , & !Pressure (Pa)
                    ta_dev          , & !Temperature (K)
                    wa_dev          , & !Specific humidity (kgkg^-1)
                    oa_dev          , & !Ozone (kgkg^-1)
                    co2             , & !CO2 (pppv)
                    cwc_dev         , & !Cloud properties (kgkg^-1)
                    fcld_dev        , & !Cloud fractions (1)
                    ict             , & !Level index separating high and middle clouds
                    icb             , & !Level index separating middle and low clouds
                    reff_dev        , & !Moisture refflectivity properties
                    hk_uv           , & !Solar UV constant
                    hk_ir           , & !Solar IR constant
                    taua_dev        , & !Aerosol optical thickness 
                    ssaa_dev        , & !Aerosol single scattering albedo
                    asya_dev        , & !Aerosol asymmetry factor
                    rsuvbm_dev      , & !Surface reflectivity in the UV+par for beam insolation
                    rsuvdf_dev      , & !Surface reflectivity in the UV+par for diffuse insolation
                    rsirbm_dev      , & !Surface reflectivity in the near-ir region for beam insolation
                    rsirdf_dev      , & !Surface reflectivity in the near-ir region for diffuse insolation
!Outputs
                    flx_dev         , & !
!Constants            
                    CONS_GRAV,                   &
                    wk_uv, zk_uv, ry_uv,         &
                    xk_ir, ry_ir,                &
                    cah, coa,                    &
                    aig_uv, awg_uv, arg_uv,      &
                    aib_uv, awb_uv, arb_uv,      &
                    aib_nir, awb_nir, arb_nir,   &
                    aia_nir, awa_nir, ara_nir,   &
                    aig_nir, awg_nir, arg_nir,   &
                    caib, caif                   ) 

   IMPLICIT NONE

   ! Parameters
   ! ----------

   integer, parameter :: nu = 43
   integer, parameter :: nw = 37
   integer, parameter :: nx = 62
   integer, parameter :: ny = 101
   integer, parameter :: nband_uv = 5
   integer, parameter :: nk_ir = 10
   integer, parameter :: nband_ir = 3

   integer, parameter :: nband = nband_uv + nband_ir

   real(8),    parameter :: dsm = 0.602


!-----input values

!-----input parameters

      integer m,np,ict,icb,nb
      real(8) cosz_dev(m),pl_dev(m,np+1),ta_dev(m,np),wa_dev(m,np),oa_dev(m,np),co2
      real(8) cwc_dev(m,np,4),fcld_dev(m,np),reff_dev(m,np,4), hk_uv(5),hk_ir(3,10)
      real(8) rsuvbm_dev(m),rsuvdf_dev(m),rsirbm_dev(m),rsirdf_dev(m)

      real(8) taua_dev(m,np,nb)
      real(8) ssaa_dev(m,np,nb)
      real(8) asya_dev(m,np,nb)

      logical :: overcast

! Constants
real(8), intent(in) :: wk_uv(5), zk_uv(5), ry_uv(5)
real(8), intent(in) :: xk_ir(10), ry_ir(3)
real(8), intent(in) :: cah(43,37), coa(62,101)

real(8), intent(in) :: aig_uv(3), awg_uv(3), arg_uv(3)
real(8), intent(in) :: aib_uv, awb_uv(2), arb_uv(2)
real(8), intent(in) :: aib_nir, awb_nir(3,2), arb_nir(3,2)
real(8), intent(in) :: aia_nir(3,3), awa_nir(3,3), ara_nir(3,3)
real(8), intent(in) :: aig_nir(3,3), awg_nir(3,3), arg_nir(3,3)
real(8), intent(in) :: caib(11,9,11), caif(9,11)

real(8), intent(in) :: CONS_GRAV

!-----output parameters

      real(8) flx_dev(m,np+1),flc_dev(m,np+1)
      real(8) flxu_dev(m,np+1),flcu_dev(m,np+1)
      real(8) fdiruv_dev (m),fdifuv_dev (m)
      real(8) fdirpar_dev(m),fdifpar_dev(m)
      real(8) fdirir_dev (m),fdifir_dev (m)
      real(8) flx_sfc_band_dev(m,nband)

!-----temporary arrays

      integer :: i,j,k,l,in,ntop

      real(8) :: dp(np),wh(np),oh(np)
      real(8) :: scal(np)
      real(8) :: swh(np+1),so2(np+1),df(0:np+1)
      real(8) :: scal0, wvtoa, o3toa, pa
      real(8) :: snt,cnt,x,xx4,xtoa
      real(8) :: dp_pa(np)

!-----parameters for co2 transmission tables

      real(8) :: w1,dw,u1,du

      integer :: ib,rc
      real(8) :: tauclb(np),tauclf(np),asycl(np)
      real(8) :: taubeam(np,4),taudiff(np,4)
      real(8) :: fcld_col(np)
      real(8) :: cwc_col(np,4)
      real(8) :: reff_col(np,4)
      real(8) :: taurs,tauoz,tauwv
      real(8) :: tausto,ssatau,asysto
      real(8) :: tautob,ssatob,asytob
      real(8) :: tautof,ssatof,asytof
      real(8) :: rr(0:np+1,2),tt(0:np+1,2),td(0:np+1,2)
      real(8) :: rs(0:np+1,2),ts(0:np+1,2)
      real(8) :: fall(np+1),fclr(np+1),fsdir,fsdif
      real(8) :: fupa(np+1),fupc(np+1)
      real(8) :: cc1,cc2,cc3
      real(8) :: rrt,ttt,tdt,rst,tst

      integer :: iv,ik
      real(8) :: ssacl(np)

      integer :: im

      integer :: ic,iw 
      real(8) :: ulog,wlog,dc,dd,x0,x1,x2,y0,y1,y2,du2,dw2

      integer :: ih

      !if (overcast == true) then
      !real(8) :: rra(0:np+1),rxa(0:np+1)
      !real(8) :: ttaold,tdaold,rsaold
      !real(8) :: ttanew,tdanew,rsanew 
      !else
      real(8) :: rra(0:np+1,2,2),tta(0:np,2,2)
      real(8) :: tda(0:np,2,2)
      real(8) :: rsa(0:np,2,2),rxa(0:np+1,2,2)
      !endif

      real(8) :: flxdn
      real(8) :: fdndir,fdndif,fupdif
      real(8) :: denm,yy

      integer :: is
      real(8) :: ch,cm,ct

      integer :: foundtop
      real(8) :: dftop

!-----Variables for aerosols

      integer :: II, JJ, irhp1, an

      real(8) :: dum


      RUN_LOOP: do i=1,m

         ntop = 0
         fdndir = 0.0
         fdndif = 0.0

!-----Beginning of sorad code

!-----wvtoa and o3toa are the water vapor and o3 amounts of the region 
!     above the pl(1) level.
!     snt is the secant of the solar zenith angle

         snt    = 1.0/cosz_dev(i)
         xtoa   = max(pl_dev(i,1),1.e-3)
         scal0  = xtoa*(0.5*xtoa/300.)**.8
         o3toa  = 1.02*oa_dev(i,1)*xtoa*466.7 + 1.0e-8
         wvtoa  = 1.02*wa_dev(i,1)*scal0 * (1.0+0.00135*(ta_dev(i,1)-240.)) + 1.0e-9
         swh(1)  = wvtoa

         do k=1,np

!-----compute layer thickness. indices for the surface level and
!     surface layer are np+1 and np, respectively.

            dp(k) = pl_dev(i,k+1)-pl_dev(i,k)
            dp_pa(k) = dp(k) * 100. ! dp in pascals
 
!-----compute scaled water vapor amount following Eqs. (3.3) and (3.5) 
!     unit is g/cm**2
!
            pa   = 0.5*(pl_dev(i,k)+pl_dev(i,k+1))
            scal(k) = dp(k)*(pa/300.)**.8
            wh(k)   = 1.02*wa_dev(i,k)*scal(k) * (1.+0.00135*(ta_dev(i,k)-240.)) + 1.e-9
            swh(k+1)= swh(k)+wh(k)

!-----compute ozone amount, unit is (cm-atm)stp
!     the number 466.7 is the unit conversion factor
!     from g/cm**2 to (cm-atm)stp

            oh(k)   = 1.02*oa_dev(i,k)*dp(k)*466.7 + 1.e-8

!-----Fill the reff, cwc, and fcld for the column

            fcld_col(k) = fcld_dev(i,k)
            do l = 1, 4
               reff_col(k,l) = reff_dev(i,k,l)
               cwc_col(k,l) = cwc_dev(i,k,l)
            end do

         end do

!-----Initialize temporary arrays to zero to avoid UMR

         rr = 0.0
         tt = 0.0
         td = 0.0
         rs = 0.0
         ts = 0.0

         rra = 0.0
         rxa = 0.0

!if( OVERCAST == .false. ) then
         tta = 0.0
         tda = 0.0
         rsa = 0.0
!endif

!-----initialize fluxes for all-sky (flx), clear-sky (flc), and
!     flux reduction (df)
!
         do k=1,np+1
            flx_dev(i,k)=0.
            flc_dev(i,k)=0.
            flxu_dev(i,k)=0.
            flcu_dev(i,k)=0.
         end do

!-----Initialize new per-band surface fluxes

         do ib = 1, nband
            flx_sfc_band_dev(i,ib) = 0.
         end do

!-----Begin inline of SOLUV

!-----compute solar uv and par fluxes

!-----initialize fdiruv, fdifuv, surface reflectances and transmittances.
!     the reflectance and transmittance of the clear and cloudy portions
!     of a layer are denoted by 1 and 2, respectively.
!     cc is the maximum cloud cover in each of the high, middle, and low
!     cloud groups.
!     1/dsm=1/cos(53) = 1.66

         fdiruv_dev(i)=0.0
         fdifuv_dev(i)=0.0
         rr(np+1,1)=rsuvbm_dev(i)
         rr(np+1,2)=rsuvbm_dev(i)
         rs(np+1,1)=rsuvdf_dev(i)
         rs(np+1,2)=rsuvdf_dev(i)
         td(np+1,1)=0.0
         td(np+1,2)=0.0
         tt(np+1,1)=0.0
         tt(np+1,2)=0.0
         ts(np+1,1)=0.0
         ts(np+1,2)=0.0
         rr(0,1)=0.0
         rr(0,2)=0.0
         rs(0,1)=0.0
         rs(0,2)=0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
         tt(0,1)=1.0
         tt(0,2)=1.0
         ts(0,1)=1.0
         ts(0,2)=1.0
         cc1=0.0
         cc2=0.0
         cc3=0.0

!-----options for scaling cloud optical thickness

!if ( OVERCAST == .true. ) then

!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

!         call getvistau1(np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
!                        taubeam,taudiff,asycl,                                  &
!                         aig_uv, awg_uv, arg_uv,                                 &
!                         aib_uv, awb_uv, arb_uv,                                 &
!                         aib_nir, awb_nir, arb_nir,                              &
!                         aia_nir, awa_nir, ara_nir,                              &
!                         aig_nir, awg_nir, arg_nir,                              &
!                         caib, caif,                                             &
!                         CONS_GRAV                                               )

!else

!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

         call getvistau1(np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,ict,icb, &
                         taubeam,taudiff,asycl,                                  &
                         aig_uv, awg_uv, arg_uv,                                 &
                         aib_uv, awb_uv, arb_uv,                                 &
                         aib_nir, awb_nir, arb_nir,                              &
                         aia_nir, awa_nir, ara_nir,                              &
                         aig_nir, awg_nir, arg_nir,                              &
                         caib, caif,                                             &
                         CONS_GRAV                                               )

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
!     The cc1,2,3 are still needed in the flux calculations below

!MAT---DO NOT FUSE THIS LOOP
!MAT---Loop must run to completion so that cc[1,2,3] are correct.
         do k=1,np
            if (k.lt.ict) then
               cc1=max(cc1,fcld_dev(i,k))
            elseif (k.lt.icb) then
               cc2=max(cc2,fcld_dev(i,k))
            else
               cc3=max(cc3,fcld_dev(i,k))
            end if
         end do
!MAT---DO NOT FUSE THIS LOOP

!endif !overcast

         do k=1,np
            tauclb(k)=taubeam(k,1)+taubeam(k,2)+taubeam(k,3)+taubeam(k,4)
            tauclf(k)=taudiff(k,1)+taudiff(k,2)+taudiff(k,3)+taudiff(k,4)
         end do

!-----integration over spectral bands

!-----Compute optical thickness, single-scattering albedo and asymmetry
!     factor for a mixture of "na" aerosol types. [Eqs. (4.16)-(4.18)]

         do ib=1,nband_uv

!-----compute direct beam transmittances of the layer above pl(1)

            td(0,1)=exp(-(wvtoa*wk_uv(ib)+o3toa*zk_uv(ib))/cosz_dev(i))
            td(0,2)=td(0,1)

            do k=1,np

!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor (Eqs. 6.2-6.4)

               taurs=ry_uv(ib)*dp(k)
               tauoz=zk_uv(ib)*oh(k)
               tauwv=wk_uv(ib)*wh(k)

               tausto=taurs+tauoz+tauwv+taua_dev(i,k,ib)+1.0e-7
               ssatau=ssaa_dev(i,k,ib)+taurs
               asysto=asya_dev(i,k,ib)

               tautob=tausto
               asytob=asysto/ssatau
               ssatob=ssatau/tautob+1.0e-8
               ssatob=min(ssatob,0.999999)

!-----for direct incident radiation

               call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

               call deledd(tautob,ssatob,asytob,dsm,rst,tst,dum)

               rr(k,1)=rrt
               tt(k,1)=ttt
               td(k,1)=tdt
               rs(k,1)=rst
               ts(k,1)=tst

!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer

!-----for direct incident radiation
!     The effective layer optical properties. Eqs. (6.2)-(6.4)

               tautob=tausto+tauclb(k)
               ssatob=(ssatau+tauclb(k))/tautob+1.0e-8
               ssatob=min(ssatob,0.999999)
               asytob=(asysto+asycl(k)*tauclb(k))/(ssatob*tautob)

!-----for diffuse incident radiation

               tautof=tausto+tauclf(k)
               ssatof=(ssatau+tauclf(k))/tautof+1.0e-8
               ssatof=min(ssatof,0.999999)
               asytof=(asysto+asycl(k)*tauclf(k))/(ssatof*tautof)

!-----for direct incident radiation
!     note that the cloud optical thickness is scaled differently 
!     for direct and diffuse insolation, Eqs. (7.3) and (7.4).

               call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation 
!     with an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

               call deledd(tautof,ssatof,asytof,dsm,rst,tst,dum)

               rr(k,2)=rrt
               tt(k,2)=ttt
               td(k,2)=tdt
               rs(k,2)=rst
               ts(k,2)=tst
            end do

!-----flux calculations
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)

            do k=1,np+1
               fclr(k)=0.0
               fall(k)=0.0
               fupa(k)=0.0
               fupc(k)=0.0
            end do

            fsdir=0.0
            fsdif=0.0

!if ( OVERCAST == .true. ) then

!-----Inline CLDFLXY

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)

!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

!-----ih=1 for clear sky; ih=2 for cloudy sky.

!-----First set is ih = 1
!            rra(np+1)=rr(np+1,1)
!            rxa(np+1)=rs(np+1,1)
!
!            do k=np,0,-1
!               denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
!               rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
!               rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
!            end do
!
!            do k=1,np+1
!               if (k <= np) then
!                  if (k == 1) then
!                     tdaold = td(0,1)
!                     ttaold = tt(0,1)
!                     rsaold = rs(0,1)
!
!                     tdanew = 0.0
!                     ttanew = 0.0
!                     rsanew = 0.0
!                  end if
!                  denm=ts(k,1)/(1.-rsaold*rs(k,1))
!                  tdanew=tdaold*td(k,1)
!                  ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
!                  rsanew=rs(k,1)+ts(k,1)*rsaold*denm
!               end if
!
!               denm=1./(1.-rsaold*rxa(k))
!               fdndir=tdaold
!               xx4=tdaold*rra(k)
!               yy=ttaold-tdaold
!               fdndif=(xx4*rsaold+yy)*denm
!               fupdif=(xx4+yy*rxa(k))*denm
!               flxdn=fdndir+fdndif-fupdif
!               fupc(k)=fupdif 
!               fclr(k)=flxdn
!
!               tdaold = tdanew
!               ttaold = ttanew
!               rsaold = rsanew
!
!               tdanew = 0.0
!               ttanew = 0.0
!               rsanew = 0.0
!            end do
!
!!-----Second set is ih = 2
!
!            rra(np+1)=rr(np+1,2)
!            rxa(np+1)=rs(np+1,2)
!
!            do k=np,0,-1
!               denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
!               rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
!               rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
!            end do
!
!            do k=1,np+1
!               if (k <= np) then
!                  if (k == 1) then
!                     tdaold = td(0,2)
!                     ttaold = tt(0,2)
!                     rsaold = rs(0,2)

!                     tdanew = 0.0
!                     ttanew = 0.0
!                     rsanew = 0.0
!                  end if
!                  denm=ts(k,2)/(1.-rsaold*rs(k,2))
!                  tdanew=tdaold*td(k,2)
!                  ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
!                  rsanew=rs(k,2)+ts(k,2)*rsaold*denm
!               end if
!
!               denm=1./(1.-rsaold*rxa(k))
!               fdndir=tdaold
!               xx4=tdaold*rra(k)
!               yy=ttaold-tdaold
!               fdndif=(xx4*rsaold+yy)*denm
!               fupdif=(xx4+yy*rxa(k))*denm
!               flxdn=fdndir+fdndif-fupdif
!
!               fupa(k)=fupdif
!               fall(k)=flxdn
!
!               tdaold = tdanew
!               ttaold = ttanew
!               rsaold = rsanew
!
!               tdanew = 0.0
!               ttanew = 0.0
!               rsanew = 0.0
!            end do
!
!            fsdir=fdndir
!            fsdif=fdndif
!
!!-----End CLDFLXY inline
!
!else

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)

!-----Inline CLDFLX

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.

!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition

            do ih=1,2
               tda(0,ih,1)=td(0,ih)
               tta(0,ih,1)=tt(0,ih)
               rsa(0,ih,1)=rs(0,ih)
               tda(0,ih,2)=td(0,ih)
               tta(0,ih,2)=tt(0,ih)
               rsa(0,ih,2)=rs(0,ih)

               do k=1,ict-1
                  denm=ts(k,ih)/(1.-rsa(k-1,ih,1)*rs(k,ih))
                  tda(k,ih,1)=tda(k-1,ih,1)*td(k,ih)
                  tta(k,ih,1)=tda(k-1,ih,1)*tt(k,ih)+(tda(k-1,ih,1)*rsa(k-1,ih,1)&
                        *rr(k,ih)+tta(k-1,ih,1)-tda(k-1,ih,1))*denm
                  rsa(k,ih,1)=rs(k,ih)+ts(k,ih)*rsa(k-1,ih,1)*denm
                  tda(k,ih,2)=tda(k,ih,1)
                  tta(k,ih,2)=tta(k,ih,1)
                  rsa(k,ih,2)=rsa(k,ih,1)
               end do ! k loop

!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition

               do k=ict,icb-1
                  do im=1,2
                     denm=ts(k,im)/(1.-rsa(k-1,ih,im)*rs(k,im))
                     tda(k,ih,im)=tda(k-1,ih,im)*td(k,im)
                     tta(k,ih,im)=tda(k-1,ih,im)*tt(k,im)+(tda(k-1,ih,im)&
                           *rsa(k-1,ih,im)*rr(k,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                     rsa(k,ih,im)=rs(k,im)+ts(k,im)*rsa(k-1,ih,im)*denm
                  end do ! im loop
               end do ! k loop
            end do ! ih loop

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.

!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition

            do is=1,2
               rra(np+1,1,is)=rr(np+1,is)
               rxa(np+1,1,is)=rs(np+1,is)
               rra(np+1,2,is)=rr(np+1,is)
               rxa(np+1,2,is)=rs(np+1,is)

               do k=np,icb,-1
                  denm=ts(k,is)/(1.-rs(k,is)*rxa(k+1,1,is))
                  rra(k,1,is)=rr(k,is)+(td(k,is)*rra(k+1,1,is)+(tt(k,is)-td(k,is))&
                        *rxa(k+1,1,is))*denm
                  rxa(k,1,is)=rs(k,is)+ts(k,is)*rxa(k+1,1,is)*denm
                  rra(k,2,is)=rra(k,1,is)
                  rxa(k,2,is)=rxa(k,1,is)
               end do ! k loop

!-----for middle clouds

               do k=icb-1,ict,-1
                  do im=1,2
                     denm=ts(k,im)/(1.-rs(k,im)*rxa(k+1,im,is))
                     rra(k,im,is)=rr(k,im)+(td(k,im)*rra(k+1,im,is)+(tt(k,im)-td(k,im))&
                           *rxa(k+1,im,is))*denm
                     rxa(k,im,is)=rs(k,im)+ts(k,im)*rxa(k+1,im,is)*denm
                  end do ! im loop
               end do ! k loop
            end do ! is loop

!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.

            do ih=1,2

!-----clear portion 
               if(ih.eq.1) then
                  ch=1.0-cc1
!-----cloudy portion
               else
                  ch=cc1
               end if

               do im=1,2
!-----clear portion 
                  if(im.eq.1) then
                     cm=ch*(1.0-cc2)
!-----cloudy portion
                  else
                     cm=ch*cc2 
                  end if

                  do is=1,2
!-----clear portion 
                     if(is.eq.1) then
                        ct=cm*(1.0-cc3) 
!-----cloudy portion
                     else
                        ct=cm*cc3
                     end if

!-----add one layer at a time, going down.

                     do k=icb,np
                        denm=ts(k,is)/(1.-rsa(k-1,ih,im)*rs(k,is))
                        tda(k,ih,im)=tda(k-1,ih,im)*td(k,is)
                        tta(k,ih,im)=tda(k-1,ih,im)*tt(k,is)+(tda(k-1,ih,im)*rr(k,is)&
                              *rsa(k-1,ih,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                        rsa(k,ih,im)=rs(k,is)+ts(k,is)*rsa(k-1,ih,im)*denm
                     end do ! k loop

!-----add one layer at a time, going up.

                     do k=ict-1,0,-1
                        denm=ts(k,ih)/(1.-rs(k,ih)*rxa(k+1,im,is))
                        rra(k,im,is)=rr(k,ih)+(td(k,ih)*rra(k+1,im,is)+(tt(k,ih)-td(k,ih))&
                              *rxa(k+1,im,is))*denm
                        rxa(k,im,is)=rs(k,ih)+ts(k,ih)*rxa(k+1,im,is)*denm
                     end do ! k loop

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)

!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

                     do k=1,np+1
                        denm=1./(1.-rsa(k-1,ih,im)*rxa(k,im,is))
                        fdndir=tda(k-1,ih,im)
                        xx4=tda(k-1,ih,im)*rra(k,im,is)
                        yy=tta(k-1,ih,im)-tda(k-1,ih,im)
                        fdndif=(xx4*rsa(k-1,ih,im)+yy)*denm
                        fupdif=(xx4+yy*rxa(k,im,is))*denm
                        flxdn=fdndir+fdndif-fupdif

!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)

                        if(ih.eq.1 .and. im.eq.1 .and. is.eq.1) then
                           fupc(k)=fupdif
                           fclr(k)=flxdn
                        end if
                        fupa(k)=fupa(k)+fupdif*ct
                        fall(k)=fall(k)+flxdn*ct
                     end do ! k loop
                     fsdir=fsdir+fdndir*ct
                     fsdif=fsdif+fdndif*ct
                  end do ! is loop
               end do ! im loop
            end do ! ih loop

!-----End CLDFLX inline

!endif !overcast

!-----flux integration, Eq. (6.1)

            do k=1,np+1
               flx_dev(i,k)=flx_dev(i,k)+fall(k)*hk_uv(ib)
               flc_dev(i,k)=flc_dev(i,k)+fclr(k)*hk_uv(ib)
               flxu_dev(i,k)=flxu_dev(i,k)+fupa(k)*hk_uv(ib)
               flcu_dev(i,k)=flcu_dev(i,k)+fupc(k)*hk_uv(ib)
            end do

!-----get surface flux for each band
            flx_sfc_band_dev(i,ib)=flx_sfc_band_dev(i,ib)+fall(np+1)*hk_uv(ib)

!-----compute direct and diffuse downward surface fluxes in the UV
!     and par regions

            if(ib.lt.5) then
               fdiruv_dev(i)=fdiruv_dev(i)+fsdir*hk_uv(ib)
               fdifuv_dev(i)=fdifuv_dev(i)+fsdif*hk_uv(ib)
            else
               fdirpar_dev(i)=fsdir*hk_uv(ib)
               fdifpar_dev(i)=fsdif*hk_uv(ib)
            end if
         end do

!-----Inline SOLIR

!-----compute and update solar ir fluxes

         fdirir_dev(i)=0.0
         fdifir_dev(i)=0.0
         rr(np+1,1)=rsirbm_dev(i)
         rr(np+1,2)=rsirbm_dev(i)
         rs(np+1,1)=rsirdf_dev(i)
         rs(np+1,2)=rsirdf_dev(i)
         td(np+1,1)=0.0
         td(np+1,2)=0.0
         tt(np+1,1)=0.0
         tt(np+1,2)=0.0
         ts(np+1,1)=0.0
         ts(np+1,2)=0.0
         rr(0,1)=0.0
         rr(0,2)=0.0
         rs(0,1)=0.0
         rs(0,2)=0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
         tt(0,1)=1.0
         tt(0,2)=1.0
         ts(0,1)=1.0
         ts(0,2)=1.0
         cc1=0.0
         cc2=0.0
         cc3=0.0

!-----integration over spectral bands

!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10)
!     The indices 1, 2, 3 are for ice, water, rain particles,
!     respectively.

         do ib=1,nband_ir
            iv=ib+5

!-----options for scaling cloud optical thickness

!if ( OVERCAST == .true. ) then

!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

!            call getnirtau1(ib,np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
!                           taubeam,taudiff,asycl,ssacl,                               &
!                            aig_uv, awg_uv, arg_uv,                                    &
!                            aib_uv, awb_uv, arb_uv,                                    &
!                            aib_nir, awb_nir, arb_nir,                                 &
!                            aia_nir, awa_nir, ara_nir,                                 &
!                            aig_nir, awg_nir, arg_nir,                                 &
!                            caib, caif,                                                &
!                            CONS_GRAV                                                  )

!else

!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

            call getnirtau1(ib,np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,ict,icb, &
                            taubeam,taudiff,asycl,ssacl,                               &
                            aig_uv, awg_uv, arg_uv,                                    &
                            aib_uv, awb_uv, arb_uv,                                    &
                            aib_nir, awb_nir, arb_nir,                                 &
                            aia_nir, awa_nir, ara_nir,                                 &
                            aig_nir, awg_nir, arg_nir,                                 &
                            caib, caif,                                                &
                            CONS_GRAV                                                  )

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

!MAT--DO NOT FUSE THIS LOOP
!MAT  Loop must run to completion so that cc[1,2,3] are correct.
            do k=1,np
               if (k.lt.ict) then
                  cc1=max(cc1,fcld_dev(i,k))
               elseif (k.lt.icb) then
                  cc2=max(cc2,fcld_dev(i,k))
               else
                  cc3=max(cc3,fcld_dev(i,k))
               end if
            end do
!MAT--DO NOT FUSE THIS LOOP

!endif !overcast

            do k=1,np
               tauclb(k)=taubeam(k,1)+taubeam(k,2)+taubeam(k,3)+taubeam(k,4)
               tauclf(k)=taudiff(k,1)+taudiff(k,2)+taudiff(k,3)+taudiff(k,4)
            end do


!-----integration over the k-distribution function

            do ik=1,nk_ir

!-----compute direct beam transmittances of the layer above pl(1)

               td(0,1)=exp(-wvtoa*xk_ir(ik)/cosz_dev(i))
               td(0,2)=td(0,1)

               do k=1,np
                  taurs=ry_ir(ib)*dp(k)
                  tauwv=xk_ir(ik)*wh(k)

!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor. Eqs.(6.2)-(6.4)

                  tausto=taurs+tauwv+taua_dev(i,k,iv)+1.0e-7
                  ssatau=ssaa_dev(i,k,iv)+taurs+1.0e-8
                  asysto=asya_dev(i,k,iv)
                  tautob=tausto
                  asytob=asysto/ssatau
                  ssatob=ssatau/tautob+1.0e-8
                  ssatob=min(ssatob,0.999999)

!-----Compute reflectance and transmittance of the clear portion 
!     of a layer

!-----for direct incident radiation

                  call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

                  call deledd(tautob,ssatob,asytob,dsm,rst,tst,dum)

                  rr(k,1)=rrt
                  tt(k,1)=ttt
                  td(k,1)=tdt
                  rs(k,1)=rst
                  ts(k,1)=tst

!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer

!-----for direct incident radiation. Eqs.(6.2)-(6.4)

                  tautob=tausto+tauclb(k)
                  ssatob=(ssatau+ssacl(k)*tauclb(k))/tautob+1.0e-8
                  ssatob=min(ssatob,0.999999)
                  asytob=(asysto+asycl(k)*ssacl(k)*tauclb(k))/(ssatob*tautob)

!-----for diffuse incident radiation

                  tautof=tausto+tauclf(k)
                  ssatof=(ssatau+ssacl(k)*tauclf(k))/tautof+1.0e-8
                  ssatof=min(ssatof,0.999999)
                  asytof=(asysto+asycl(k)*ssacl(k)*tauclf(k))/(ssatof*tautof)

!-----for direct incident radiation

                  call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs.(6.5) and (6.6)

                  call deledd(tautof,ssatof,asytof,dsm,rst,tst,dum)

                  rr(k,2)=rrt
                  tt(k,2)=ttt
                  td(k,2)=tdt
                  rs(k,2)=rst
                  ts(k,2)=tst
               end do

!-----FLUX CALCULATIONS

!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)

               do k=1,np+1
                  fclr(k)=0.0
                  fall(k)=0.0
                  fupc(k)=0.0
                  fupa(k)=0.0
               end do

               fsdir=0.0
               fsdif=0.0

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.

!if ( OVERCAST == .true. ) then

!-----Inline CLDFLXY

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)

!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

!-----ih=1 for clear sky; ih=2 for cloudy sky.

!-----First set is ih = 1
!               rra(np+1)=rr(np+1,1)
!               rxa(np+1)=rs(np+1,1)
!
!               do k=np,0,-1
!                  denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
!                  rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
!                  rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
!               end do
!
!               do k=1,np+1
!                  if (k <= np) then
!                     if (k == 1) then
!                        tdaold = td(0,1)
!                        ttaold = tt(0,1)
!                        rsaold = rs(0,1)
!
!                        tdanew = 0.0
!                        ttanew = 0.0
!                        rsanew = 0.0
!                     end if
!                     denm=ts(k,1)/(1.-rsaold*rs(k,1))
!                     tdanew=tdaold*td(k,1)
!                     ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
!                     rsanew=rs(k,1)+ts(k,1)*rsaold*denm
!                  end if
!
!                  denm=1./(1.-rsaold*rxa(k))
!                  fdndir=tdaold
!                  xx4=tdaold*rra(k)
!                  yy=ttaold-tdaold
!                  fdndif=(xx4*rsaold+yy)*denm
!                  fupdif=(xx4+yy*rxa(k))*denm
!                  flxdn=fdndir+fdndif-fupdif
!
!                  fupc(k)=fupdif
!                  fclr(k)=flxdn
!
!                  tdaold = tdanew
!                  ttaold = ttanew
!                  rsaold = rsanew
!
!                  tdanew = 0.0
!                  ttanew = 0.0
!                  rsanew = 0.0
!               end do
!
!!-----Second set is ih = 2
!
!               rra(np+1)=rr(np+1,2)
!               rxa(np+1)=rs(np+1,2)
!
!               do k=np,0,-1
!                  denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
!                  rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
!                  rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
!               end do
!
!               do k=1,np+1
!                  if (k <= np) then
!                     if (k == 1) then
!                        tdaold = td(0,2)
!                        ttaold = tt(0,2)
!                        rsaold = rs(0,2)
!
!                        tdanew = 0.0
!                        ttanew = 0.0
!                        rsanew = 0.0
!                     end if
!                     denm=ts(k,2)/(1.-rsaold*rs(k,2))
!                     tdanew=tdaold*td(k,2)
!                     ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
!                     rsanew=rs(k,2)+ts(k,2)*rsaold*denm
!                  end if
!
!                  denm=1./(1.-rsaold*rxa(k))
!                  fdndir=tdaold
!                  xx4=tdaold*rra(k)
!                  yy=ttaold-tdaold
!                  fdndif=(xx4*rsaold+yy)*denm
!                  fupdif=(xx4+yy*rxa(k))*denm
!                  flxdn=fdndir+fdndif-fupdif
!
!                  fupa(k)=fupdif
!                  fall(k)=flxdn
!
!                  tdaold = tdanew
!                  ttaold = ttanew
!                  rsaold = rsanew
!
!                  tdanew = 0.0
!                  ttanew = 0.0
!                  rsanew = 0.0
!               end do
!
!               fsdir=fdndir
!               fsdif=fdndif
!
!!-----End CLDFLXY inline
!
!else

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)

!-----Inline CLDFLX

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.

!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition

               do ih=1,2
                  tda(0,ih,1)=td(0,ih)
                  tta(0,ih,1)=tt(0,ih)
                  rsa(0,ih,1)=rs(0,ih)
                  tda(0,ih,2)=td(0,ih)
                  tta(0,ih,2)=tt(0,ih)
                  rsa(0,ih,2)=rs(0,ih)

                  do k=1,ict-1
                     denm=ts(k,ih)/(1.-rsa(k-1,ih,1)*rs(k,ih))
                     tda(k,ih,1)=tda(k-1,ih,1)*td(k,ih)
                     tta(k,ih,1)=tda(k-1,ih,1)*tt(k,ih)+(tda(k-1,ih,1)*rsa(k-1,ih,1)&
                           *rr(k,ih)+tta(k-1,ih,1)-tda(k-1,ih,1))*denm
                     rsa(k,ih,1)=rs(k,ih)+ts(k,ih)*rsa(k-1,ih,1)*denm
                     tda(k,ih,2)=tda(k,ih,1)
                     tta(k,ih,2)=tta(k,ih,1)
                     rsa(k,ih,2)=rsa(k,ih,1)
                  end do ! k loop

!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition

                  do k=ict,icb-1
                     do im=1,2
                        denm=ts(k,im)/(1.-rsa(k-1,ih,im)*rs(k,im))
                        tda(k,ih,im)=tda(k-1,ih,im)*td(k,im)
                        tta(k,ih,im)=tda(k-1,ih,im)*tt(k,im)+(tda(k-1,ih,im)&
                              *rsa(k-1,ih,im)*rr(k,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                        rsa(k,ih,im)=rs(k,im)+ts(k,im)*rsa(k-1,ih,im)*denm
                     end do ! im loop
                  end do ! k loop
               end do ! ih loop

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.

!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition

               do is=1,2
                  rra(np+1,1,is)=rr(np+1,is)
                  rxa(np+1,1,is)=rs(np+1,is)
                  rra(np+1,2,is)=rr(np+1,is)
                  rxa(np+1,2,is)=rs(np+1,is)

                  do k=np,icb,-1
                     denm=ts(k,is)/(1.-rs(k,is)*rxa(k+1,1,is))
                     rra(k,1,is)=rr(k,is)+(td(k,is)*rra(k+1,1,is)+(tt(k,is)-td(k,is))&
                           *rxa(k+1,1,is))*denm
                     rxa(k,1,is)=rs(k,is)+ts(k,is)*rxa(k+1,1,is)*denm
                     rra(k,2,is)=rra(k,1,is)
                     rxa(k,2,is)=rxa(k,1,is)
                  end do ! k loop

!-----for middle clouds

                  do k=icb-1,ict,-1
                     do im=1,2
                        denm=ts(k,im)/(1.-rs(k,im)*rxa(k+1,im,is))
                        rra(k,im,is)=rr(k,im)+(td(k,im)*rra(k+1,im,is)+(tt(k,im)-td(k,im))&
                              *rxa(k+1,im,is))*denm
                        rxa(k,im,is)=rs(k,im)+ts(k,im)*rxa(k+1,im,is)*denm
                     end do ! im loop
                  end do ! k loop
               end do ! is loop

!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.

               do ih=1,2
!-----clear portion 
                  if(ih.eq.1) then
                     ch=1.0-cc1
!-----cloudy portion
                  else
                     ch=cc1
                  end if

                  do im=1,2
!-----clear portion 
                     if(im.eq.1) then
                        cm=ch*(1.0-cc2)
!-----cloudy portion
                     else
                        cm=ch*cc2 
                     end if

                     do is=1,2
!-----clear portion 
                        if(is.eq.1) then
                           ct=cm*(1.0-cc3) 
!-----cloudy portion
                        else
                           ct=cm*cc3
                        end if

!-----add one layer at a time, going down.

                        do k=icb,np
                           denm=ts(k,is)/(1.-rsa(k-1,ih,im)*rs(k,is))
                           tda(k,ih,im)=tda(k-1,ih,im)*td(k,is)
                           tta(k,ih,im)=tda(k-1,ih,im)*tt(k,is)+(tda(k-1,ih,im)*rr(k,is)&
                                 *rsa(k-1,ih,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                           rsa(k,ih,im)=rs(k,is)+ts(k,is)*rsa(k-1,ih,im)*denm
                        end do ! k loop

!-----add one layer at a time, going up.

                        do k=ict-1,0,-1
                           denm=ts(k,ih)/(1.-rs(k,ih)*rxa(k+1,im,is))
                           rra(k,im,is)=rr(k,ih)+(td(k,ih)*rra(k+1,im,is)+(tt(k,ih)-td(k,ih))&
                                 *rxa(k+1,im,is))*denm
                           rxa(k,im,is)=rs(k,ih)+ts(k,ih)*rxa(k+1,im,is)*denm
                        end do ! k loop

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

                        do k=1,np+1
                           denm=1./(1.-rsa(k-1,ih,im)*rxa(k,im,is))
                           fdndir=tda(k-1,ih,im)
                           xx4=tda(k-1,ih,im)*rra(k,im,is)
                           yy=tta(k-1,ih,im)-tda(k-1,ih,im)
                           fdndif=(xx4*rsa(k-1,ih,im)+yy)*denm
                           fupdif=(xx4+yy*rxa(k,im,is))*denm
                           flxdn=fdndir+fdndif-fupdif

!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)

                           if(ih.eq.1 .and. im.eq.1 .and. is.eq.1) then
                              fupc(k)=fupdif
                              fclr(k)=flxdn
                           end if
                           fupa(k)=fupa(k)+fupdif*ct
                           fall(k)=fall(k)+flxdn*ct
                        end do ! k loop
                        fsdir=fsdir+fdndir*ct
                        fsdif=fsdif+fdndif*ct
                     end do ! is loop
                  end do ! im loop
               end do ! ih loop

!-----End CLDFLX inline
!endif !overcast

!-----flux integration following Eq. (6.1)

               do k=1,np+1
                  flx_dev(i,k)=flx_dev(i,k)+fall(k)*hk_ir(ib,ik)
                  flc_dev(i,k)=flc_dev(i,k)+fclr(k)*hk_ir(ib,ik)
                  flxu_dev(i,k)=flxu_dev(i,k)+fupa(k)*hk_ir(ib,ik)
                  flcu_dev(i,k)=flcu_dev(i,k)+fupc(k)*hk_ir(ib,ik)
               end do

!-----compute downward surface fluxes in the ir region

               fdirir_dev(i)=fdirir_dev(i)+fsdir*hk_ir(ib,ik)
               fdifir_dev(i)=fdifir_dev(i)+fsdif*hk_ir(ib,ik)

!-----tabulate surface flux at ir bands
               flx_sfc_band_dev(i,iv)=flx_sfc_band_dev(i,iv)+fall(np+1)*hk_ir(ib,ik)

            end do ! ik loop
         end do

!-----compute pressure-scaled o2 amount following Eq. (3.5) with f=1.
!     unit is (cm-atm)stp. 165.22 = (1000/980)*23.14%*(22400/32)
!     compute flux reduction due to oxygen following Eq. (3.18). 0.0633 is the
!     fraction of insolation contained in the oxygen bands

         df(0) = 0.0
         cnt = 165.22*snt
         so2(1) = scal0*cnt
! LLT increased parameter 145 to 155 to enhance effect
         df(1) = 0.0633*(1.-exp(-0.000155*sqrt(so2(1))))

         do k=1,np
            so2(k+1) = so2(k) + scal(k)*cnt
! LLT increased parameter 145 to 155 to enhance effect
            df(k+1) = 0.0633*(1.0 - exp(-0.000155*sqrt(so2(k+1)))) 
         end do

!-----for solar heating due to co2 scaling follows Eq(3.5) with f=1.
!     unit is (cm-atm)stp. 789 = (1000/980)*(44/28.97)*(22400/44)

         so2(1) = (789.*co2)*scal0

         do k=1,np
            so2(k+1) = so2(k) + (789.*co2)*scal(k)
         end do

!-----The updated flux reduction for co2 absorption in Band 7 where absorption due to
!     water vapor and co2 are both moderate. df is given by the second term on the
!     right-hand-side of Eq. (3.24) divided by So. so2 and swh are the co2 and
!     water vapor amounts integrated from the top of the atmosphere

         u1 = -3.0
         du = 0.15
         w1 = -4.0
         dw = 0.15

!-----Inline RFLX
         du2=du*du
         dw2=dw*dw
         x0=u1+real(nu)*du
         y0=w1+real(nw)*dw

         x1=u1-0.5*du
         y1=w1-0.5*dw

         do k= 1, np+1
            ulog=min(log10(so2(k)*snt),x0)
            wlog=min(log10(swh(k)*snt),y0)
            ic=int((ulog-x1)/du+1.)
            iw=int((wlog-y1)/dw+1.)
            if(ic.lt.2)  ic=2
            if(iw.lt.2)  iw=2
            if(ic.gt.nu) ic=nu
            if(iw.gt.nw) iw=nw
            dc=ulog-real(ic-2)*du-u1
            dd=wlog-real(iw-2)*dw-w1   
            x2=cah(ic-1,iw-1)+(cah(ic-1,iw)-cah(ic-1,iw-1))/dw*dd
            y2=x2+(cah(ic,iw-1)-cah(ic-1,iw-1))/du*dc
            y2=max(y2,0.0)
            df(k)=df(k)+1.5*y2 ! LLT increase CO2 effect to help reduce cold tropopause bias
         end do

!-----df is the updated flux reduction for co2 absorption
!     in Band 8 where the co2 absorption has a large impact
!     on the heating of middle atmosphere. From the table
!     given by Eq. (3.19)

         u1 = 0.000250
         du = 0.000050
         w1 = -2.0
         dw = 0.05

!-----Inline RFLX
         du2=du*du
         dw2=dw*dw
         x0=u1+real(nx)*du
         y0=w1+real(ny)*dw

         x1=u1-0.5*du
         y1=w1-0.5*dw

         do k= 1, np+1
            ulog=min(co2*snt,x0)
            wlog=min(log10(pl_dev(i,k)),y0)
            ic=int((ulog-x1)/du+1.)
            iw=int((wlog-y1)/dw+1.)
            if(ic.lt.2)  ic=2
            if(iw.lt.2)  iw=2
            if(ic.gt.nx) ic=nx
            if(iw.gt.ny) iw=ny
            dc=ulog-real(ic-2)*du-u1
            dd=wlog-real(iw-2)*dw-w1   
            x2=coa(ic-1,iw-1)+(coa(ic-1,iw)-coa(ic-1,iw-1))/dw*dd
            y2=x2+(coa(ic,iw-1)-coa(ic-1,iw-1))/du*dc
            y2=max(y2,0.0)
            df(k)=df(k)+1.5*y2 ! LLT increase CO2 effect to help reduce cold tropopause bias
         end do

!-----adjust the o2-co2 reduction below cloud top following Eq. (6.18)

         foundtop = 0

         do k=1,np
            if (fcld_dev(i,k) > 0.02.and.foundtop.eq.0) then
               foundtop = 1
               ntop = k
            end if
         end do

         if (foundtop.eq.0) ntop=np+1

         dftop = df(ntop)

         do k=1,np+1
            if (k .gt. ntop) then
               xx4   = (flx_dev(i,k)/flx_dev(i,ntop))
               df(k) = dftop + xx4 * (df(k)-dftop)
            end if
         end do

!-----update the net fluxes

         do k=1,np+1
            df(k) = min(df(k),flx_dev(i,k)-1.0e-8)
!           df(k) = 0.0
            flx_dev(i,k) = flx_dev(i,k) - df(k)
            flc_dev(i,k) = flc_dev(i,k) - df(k)
         end do

!-----update the downward surface fluxes 

!        xx4 = fdirir (i) + fdifir (i) +&
!              fdiruv (i) + fdifuv (i) +&
!              fdirpar(i) + fdifpar(i)

         xx4 = flx_dev(i,np+1) + df(np+1)

         if ( abs(xx4) > epsilon(1.0) ) then
            xx4 = max(min(1.0 - df(np+1)/xx4,1.),0.)
         else
            xx4 = 0.0
         end if

         fdirir_dev(i)  = xx4*fdirir_dev(i) 
         fdifir_dev(i)  = xx4*fdifir_dev(i) 
         fdiruv_dev(i)  = xx4*fdiruv_dev(i) 
         fdifuv_dev(i)  = xx4*fdifuv_dev(i) 
         fdirpar_dev(i) = xx4*fdirpar_dev(i)
         fdifpar_dev(i) = xx4*fdifpar_dev(i)

         do ib = 1, nband
            flx_sfc_band_dev(i,ib) = xx4*flx_sfc_band_dev(i,ib)
         end do

      end do RUN_LOOP

   end subroutine sorad


!*********************************************************************

   subroutine deledd(tau1,ssc1,g01,cza1,rr1,tt1,td1)

!*********************************************************************
!
!-----uses the delta-eddington approximation to compute the
!     bulk scattering properties of a single layer
!     coded following King and Harshvardhan (JAS, 1986)
!
!  inputs:
!     tau1:  optical thickness
!     ssc1:  single scattering albedo
!     g01:   asymmetry factor
!     cza1:  cosine o the zenith angle
!
!  outputs:
!
!     rr1:  reflection of the direct beam
!     tt1:  total (direct+diffuse) transmission of the direct beam
!     td1:  direct transmission of the direct beam
!
!*********************************************************************

      implicit none

      integer,parameter :: REAL_DE = 8 ! 8 byte real
      !integer,parameter :: REAL_SP = 4 ! 4 byte real

!-----input parameters

      real(8), intent(in) :: tau1,ssc1,g01,cza1

!-----output parameters

      real(8), intent(out) :: rr1, tt1, td1

!-----temporary parameters

      real(8), parameter :: zero = 0.0_REAL_DE
      real(8), parameter :: one = 1.0_REAL_DE
      real(8), parameter :: two = 2.0_REAL_DE
      real(8), parameter :: three = 3.0_REAL_DE
      real(8), parameter :: four = 4.0_REAL_DE
      real(8), parameter :: fourth = 0.25_REAL_DE
      real(8), parameter :: seven = 7.0_REAL_DE
      real(8), parameter :: thresh = 1.e-8_REAL_DE

      real(8) ::  tau,ssc,g0,rr,tt,td 
      real(8) ::  zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2
      real(8) ::  all,bll,st7,st8,cll,dll,fll,ell,st1,st2,st3,st4

      !zth = real(cza1,kind=REAL_DE)
      !g0  = real(g01 ,kind=REAL_DE)
      !tau = real(tau1,kind=REAL_DE)
      !ssc = real(ssc1,kind=REAL_DE)

      zth = dble(cza1)
      g0  = dble(g01)
      tau = dble(tau1)
      ssc = dble(ssc1)

      ff  = g0*g0
      xx  = one-ff*ssc
      taup= tau*xx
      sscp= ssc*(one-ff)/xx
      gp  = g0/(one+g0)

      xx  = three*gp 
      gm1 = (seven-sscp*(four+xx))*fourth
      gm2 =-(one  -sscp*(four-xx))*fourth

      akk = sqrt((gm1+gm2)*(gm1-gm2))

      xx  = akk*zth
      st7 = one-xx
      st8 = one+xx
      st3 = st7*st8

      if (abs(st3) .lt. thresh) then
         zth = zth+0.0010
         if(zth > 1.0) zth = zth-0.0020 
         xx  = akk*zth
         st7 = one-xx
         st8 = one+xx
         st3 = st7*st8
      end if

      td=exp(-taup/zth)

      gm3 = (two-zth*three*gp)*fourth
      xx  = gm1-gm2
      alf1= gm1-gm3*xx
      alf2= gm2+gm3*xx

      xx  = akk*two
      all = (gm3-alf2*zth    )*xx*td
      bll = (one-gm3+alf1*zth)*xx

      xx  = akk*gm3
      cll = (alf2+xx)*st7
      dll = (alf2-xx)*st8

      xx  = akk*(one-gm3)
      fll = (alf1+xx)*st8
      ell = (alf1-xx)*st7

      st2 = exp(-akk*taup)
      st4 = st2*st2

      st1 = sscp/((akk+gm1+(akk-gm1)*st4)*st3)

      rr = ( cll-dll*st4     - all*st2)*st1
      tt =-((fll-ell*st4)*td - bll*st2)*st1

      rr = max(rr,zero)
      tt = max(tt,zero)

      tt = tt+td

      !td1 = real(td,kind=REAL_SP)
      !rr1 = real(rr,kind=REAL_SP)
      !tt1 = real(tt,kind=REAL_SP)
      td1 = real(td)
      rr1 = real(rr)
      tt1 = real(tt)

   end subroutine deledd


   subroutine getvistau1(nlevs,cosz,dp,fcld,reff,hydromets,ict,icb,  &
                         taubeam,taudiff,asycl,                      &
                         aig_uv, awg_uv, arg_uv,                     &
                         aib_uv, awb_uv, arb_uv,                     &
                         aib_nir, awb_nir, arb_nir,                  &
                         aia_nir, awa_nir, ara_nir,                  &
                         aig_nir, awg_nir, arg_nir,                  &
                         caib, caif,                                 &
                         CONS_GRAV                                   )

! !USES:

      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: nlevs              !  Number of levels
      real(8),    intent(IN ) :: cosz               !  Cosine of solar zenith angle
      real(8),    intent(IN ) :: dp(nlevs)          !  Delta pressure (Pa)
      real(8),    intent(IN ) :: fcld(nlevs)        !  Cloud fraction (used sometimes)
      real(8),    intent(IN ) :: reff(nlevs,4)      !  Effective radius (microns)
      real(8),    intent(IN ) :: hydromets(nlevs,4) !  Hydrometeors (kg/kg)
      integer, intent(IN ) :: ict, icb           !  Flags for various uses 
!                 ict  = 0   Indicates that in-cloud values have been given
!                            and are expected
!                 ict != 0   Indicates that overlap computation is needed, and:
!                               ict is the level of the mid-high boundary
!                               icb is the level of the low-mid  boundary
!                
! !OUTPUT PARAMETERS:
      real(8),    intent(OUT) :: taubeam(nlevs,4)    !  Optical Depth for Beam Radiation
      real(8),    intent(OUT) :: taudiff(nlevs,4)    !  Optical Depth for Diffuse Radiation
      real(8),    intent(OUT) ::   asycl(nlevs  )    !  Cloud Asymmetry Factor
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for visible wavelengths
!  In general will compute in-cloud - will do grid mean when called
!  for diagnostic use only. ict flag will indicate which to do.
!  Slots for reff, hydrometeors, taubeam, taudiff, and asycl are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.10.27   Molod moved to Radiation_Shared and revised arg list, units
!    2011.11.16   MAT: Generalized to a call that is per-column
!
!EOP
!------------------------------------------------------------------------------
!BOC

      integer            :: k,in,im,it,ia,kk
      real(8)               :: fm,ft,fa,xai,tauc,asyclt
      real(8)               :: cc(3)
      real(8)               :: taucld1,taucld2,taucld3,taucld4
      real(8)               :: g1,g2,g3,g4

      real(8)               :: reff_snow

      integer, parameter :: nm=11,nt=9,na=11
      real(8),    parameter :: dm=0.1,dt=0.30103,da=0.1,t1=-0.9031

      real(8), intent(in) :: aig_uv(3), awg_uv(3), arg_uv(3)
      real(8), intent(in) :: aib_uv, awb_uv(2), arb_uv(2)
      real(8), intent(in) :: aib_nir, awb_nir(3,2), arb_nir(3,2)
      real(8), intent(in) :: aia_nir(3,3), awa_nir(3,3), ara_nir(3,3)
      real(8), intent(in) :: aig_nir(3,3), awg_nir(3,3), arg_nir(3,3)
      real(8), intent(in) :: caib(11,9,11), caif(9,11)
      real(8), intent(in) :: CONS_GRAV

      taubeam = 0.0
      taudiff = 0.0

      if (ict.ne.0) then

!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

         cc = 0.0

         do k = 1, ict-1
             cc(1)=max(cc(1),fcld(k))
         end do
         do k = ict, icb-1
             cc(2)=max(cc(2),fcld(k))
         end do
         do k = icb, nlevs
             cc(3)=max(cc(3),fcld(k))
         end do

      end if

!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow

      do k = 1, nlevs

         if (reff(k,1) <= 0.) then
            taucld1=0.
         else
            taucld1=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,1))*aib_uv/reff(k,1)
         end if

         if (reff(k,2) <= 0.) then
            taucld2=0.
         else
            taucld2=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,2))*(awb_uv(1)+awb_uv(2)/reff(k,2))
         end if

            taucld3=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,3))*arb_uv(1)

!-----In the IR optical thickness calculation (getirtau.code), it was
!     found that using the table of coefficients tabulated for suspended
!     cloud ice particles (aib_ir) for falling snow lead to unphysical 
!     (negative) values of cloud optical thickness for effective radii 
!     greater than 113 microns. By restricting the effective radius of 
!     snow to 112 microns, we prevent unphysical optical thicknesses. 
!     For consistency's sake, we limit snow effective radius similarly here.

         reff_snow = min(reff(k,4),112.0)

         if (reff_snow <= 0.) then
            taucld4=0.
         else
            taucld4=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,4))*aib_uv/reff_snow
         endif

         if ( ict .ne. 0 ) then

!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

            if (k.lt.ict) then
               kk=1
            else if (k.ge.ict .and. k.lt.icb) then
               kk=2
            else
               kk=3
            end if
 
            tauc=taucld1+taucld2+taucld3+taucld4

            if (tauc.gt.0.02 .and. fcld(k).gt.0.01) then

!-----normalize cloud cover following Eq. (7.8)

               fa=fcld(k)/cc(kk)

!-----table look-up
 
               tauc=min(tauc,32.)

               fm=cosz/dm 
               ft=(log10(tauc)-t1)/dt
               fa=fa/da

               im=int(fm+1.5) 
               it=int(ft+1.5)
               ia=int(fa+1.5)
  
               im=max(im,2)
               it=max(it,2)
               ia=max(ia,2)
     
               im=min(im,nm-1)
               it=min(it,nt-1)
               ia=min(ia,na-1)
 
               fm=fm-real(im-1)
               ft=ft-real(it-1)
               fa=fa-real(ia-1)

!-----scale cloud optical thickness for beam radiation following 
!     Eq. (7.3).
!     the scaling factor, xai, is a function of the solar zenith
!     angle, optical thickness, and cloud cover.

               xai=    (-caib(im-1,it,ia)*(1.-fm)+&
                         caib(im+1,it,ia)*(1.+fm))*fm*.5+caib(im,it,ia)*(1.-fm*fm)

               xai=xai+(-caib(im,it-1,ia)*(1.-ft)+&
                         caib(im,it+1,ia)*(1.+ft))*ft*.5+caib(im,it,ia)*(1.-ft*ft)

               xai=xai+(-caib(im,it,ia-1)*(1.-fa)+&
                         caib(im,it,ia+1)*(1.+fa))*fa*.5+caib(im,it,ia)*(1.-fa*fa)

               xai=xai-2.*caib(im,it,ia)

               xai=max(xai,0.0)
               xai=min(xai,1.0)

               taubeam(k,1)=taucld1*xai
               taubeam(k,2)=taucld2*xai
               taubeam(k,3)=taucld3*xai
               taubeam(k,4)=taucld4*xai

!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
 
               xai=    (-caif(it-1,ia)*(1.-ft)+&
                         caif(it+1,ia)*(1.+ft))*ft*.5+caif(it,ia)*(1.-ft*ft)
 
               xai=xai+(-caif(it,ia-1)*(1.-fa)+&
                         caif(it,ia+1)*(1.+fa))*fa*.5+caif(it,ia)*(1.-fa*fa)
 
               xai=xai-caif(it,ia)
 
               xai=max(xai,0.0)
               xai=min(xai,1.0)
 
               taudiff(k,1)=taucld1*xai
               taudiff(k,2)=taucld2*xai
               taudiff(k,3)=taucld3*xai
               taudiff(k,4)=taucld4*xai
            end if
         else
         ! Overlap calculation scaling not needed
            taubeam(k,1)=taucld1
            taubeam(k,2)=taucld2
            taubeam(k,3)=taucld3
            taubeam(k,4)=taucld4

            taudiff(k,1)=taucld1
            taudiff(k,2)=taucld2
            taudiff(k,3)=taucld3
            taudiff(k,4)=taucld4
         end if

!-----cloud asymmetry factor for a mixture of liquid and ice particles.
!     unit of reff is micrometers. Eqs. (4.8) and (6.4)

         asyclt=1.0
         tauc=taucld1+taucld2+taucld3+taucld4

         if (tauc.gt.0.02 .and. fcld(k).gt.0.01) then
            g1=(aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff(k,1))*reff(k,1))*taucld1
            g2=(awg_uv(1)+(awg_uv(2)+awg_uv(3)*reff(k,2))*reff(k,2))*taucld2
            g3= arg_uv(1)                                           *taucld3
            g4=(aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff_snow)*reff_snow)*taucld4
            asyclt=(g1+g2+g3+g4)/tauc
         end if

         asycl(k)=asyclt

      end do

      return

!EOC
   end subroutine getvistau1



   subroutine getnirtau1(ib,nlevs,cosz,dp,fcld,reff,hydromets,ict,icb, &
                         taubeam,taudiff,asycl,ssacl,                  &
                         aig_uv, awg_uv, arg_uv,                       &
                         aib_uv, awb_uv, arb_uv,                       &
                         aib_nir, awb_nir, arb_nir,                    &
                         aia_nir, awa_nir, ara_nir,                    &
                         aig_nir, awg_nir, arg_nir,                    &
                         caib, caif,                                   &
                         CONS_GRAV                                     )


      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: ib                 !  Band number
      integer, intent(IN ) :: nlevs              !  Number of levels
      real(8),    intent(IN ) :: cosz               !  Cosine of solar zenith angle
      real(8),    intent(IN ) :: dp(nlevs)          !  Delta pressure in Pa
      real(8),    intent(IN ) :: fcld(nlevs)        !  Cloud fraction (used sometimes)
      real(8),    intent(IN ) :: reff(nlevs,4)      !  Effective radius (microns)
      real(8),    intent(IN ) :: hydromets(nlevs,4) !  Hydrometeors (kg/kg)
      integer, intent(IN ) :: ict, icb           !  Flags for various uses 

      real(8), intent(in) :: aig_uv(3), awg_uv(3), arg_uv(3)
      real(8), intent(in) :: aib_uv, awb_uv(2), arb_uv(2)
      real(8), intent(in) :: aib_nir, awb_nir(3,2), arb_nir(3,2)
      real(8), intent(in) :: aia_nir(3,3), awa_nir(3,3), ara_nir(3,3)
      real(8), intent(in) :: aig_nir(3,3), awg_nir(3,3), arg_nir(3,3)
      real(8), intent(in) :: caib(11,9,11), caif(9,11)
      real(8), intent(in) :: CONS_GRAV

! !OUTPUT PARAMETERS:
      real(8),    intent(OUT) :: taubeam(nlevs,4)   !  Optical depth for beam radiation
      real(8),    intent(OUT) :: taudiff(nlevs,4)   !  Optical depth for diffuse radiation
      real(8),    intent(OUT) ::   ssacl(nlevs  )   !  Cloud single scattering albedo
      real(8),    intent(OUT) ::   asycl(nlevs  )   !  Cloud asymmetry factor

      integer            :: k,in,im,it,ia,kk
      real(8)               :: fm,ft,fa,xai,tauc,asyclt,ssaclt
      real(8)               :: cc(3)
      real(8)               :: taucld1,taucld2,taucld3,taucld4
      real(8)               :: g1,g2,g3,g4
      real(8)               :: w1,w2,w3,w4

      real(8)               :: reff_snow

      integer, parameter :: nm=11,nt=9,na=11
      real(8),    parameter :: dm=0.1,dt=0.30103,da=0.1,t1=-0.9031

      taubeam = 0.0
      taudiff = 0.0

      if (ict.ne.0) then

!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

         cc = 0.0

         do k = 1, ict-1
             cc(1)=max(cc(1),fcld(k))
         end do
         do k = ict, icb-1
             cc(2)=max(cc(2),fcld(k))
         end do
         do k = icb, nlevs
             cc(3)=max(cc(3),fcld(k))
         end do

      end if

!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow

      do k = 1, nlevs

         if (reff(k,1) <= 0.) then
            taucld1=0.
         else
            taucld1=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,1))*aib_nir/reff(k,1)
         end if

         if (reff(k,2) <= 0.) then
            taucld2=0.
         else
            taucld2=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,2))*(awb_nir(ib,1)+awb_nir(ib,2)/reff(k,2))
         end if

            taucld3=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,3))*arb_nir(ib,1)

!-----In the IR optical thickness calculation (getirtau.code), it was 
!     found that using the table of coefficients tabulated for suspended
!     cloud ice particles (aib_ir) for falling snow lead to unphysical 
!     (negative) values of cloud optical thickness for effective radii 
!     greater than 113 microns. By restricting the effective radius of  
!     snow to 112 microns, we prevent unphysical optical thicknesses. 
!     For consistency's sake, we limit snow effective radius similarly here.

         reff_snow = min(reff(k,4),112.0)

         if (reff_snow <= 0.) then
            taucld4=0.
         else
            taucld4=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,4))*aib_nir/reff_snow
         endif

         if ( ict .ne. 0 ) then

!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

            if (k.lt.ict) then
               kk=1
            else if (k.ge.ict .and. k.lt.icb) then
               kk=2
            else
               kk=3
            end if
 
            tauc=taucld1+taucld2+taucld3+taucld4

            if (tauc.gt.0.02 .and. fcld(k).gt.0.01) then

!-----normalize cloud cover following Eq. (7.8)
             if (cc(kk).ne.0.0) then
                 fa=fcld(k)/cc(kk)
	     else
	         fa=0.0
	     end if

!-----table look-up
 
               tauc=min(tauc,32.)

               fm=cosz/dm 
               ft=(log10(tauc)-t1)/dt
               fa=fa/da

               im=int(fm+1.5) 
               it=int(ft+1.5)
               ia=int(fa+1.5)
  
               im=max(im,2)
               it=max(it,2)
               ia=max(ia,2)
     
               im=min(im,nm-1)
               it=min(it,nt-1)
               ia=min(ia,na-1)
 
               fm=fm-real(im-1)
               ft=ft-real(it-1)
               fa=fa-real(ia-1)

!-----scale cloud optical thickness for beam radiation following 
!     Eq. (7.3).
!     the scaling factor, xai, is a function of the solar zenith
!     angle, optical thickness, and cloud cover.

               xai=    (-caib(im-1,it,ia)*(1.-fm)+&
                         caib(im+1,it,ia)*(1.+fm))*fm*.5+caib(im,it,ia)*(1.-fm*fm)

               xai=xai+(-caib(im,it-1,ia)*(1.-ft)+&
                         caib(im,it+1,ia)*(1.+ft))*ft*.5+caib(im,it,ia)*(1.-ft*ft)

               xai=xai+(-caib(im,it,ia-1)*(1.-fa)+&
                         caib(im,it,ia+1)*(1.+fa))*fa*.5+caib(im,it,ia)*(1.-fa*fa)

               xai=xai-2.*caib(im,it,ia)

               xai=max(xai,0.0)
               xai=min(xai,1.0)

               taubeam(k,1)=taucld1*xai
               taubeam(k,2)=taucld2*xai
               taubeam(k,3)=taucld3*xai
               taubeam(k,4)=taucld4*xai

!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
 
               xai=    (-caif(it-1,ia)*(1.-ft)+&
                         caif(it+1,ia)*(1.+ft))*ft*.5+caif(it,ia)*(1.-ft*ft)
 
               xai=xai+(-caif(it,ia-1)*(1.-fa)+&
                         caif(it,ia+1)*(1.+fa))*fa*.5+caif(it,ia)*(1.-fa*fa)
 
               xai=xai-caif(it,ia)
 
               xai=max(xai,0.0)
               xai=min(xai,1.0)
 
               taudiff(k,1)=taucld1*xai
               taudiff(k,2)=taucld2*xai
               taudiff(k,3)=taucld3*xai
               taudiff(k,4)=taucld4*xai
            end if
         else
         ! Overlap calculation scaling not needed
            taubeam(k,1)=taucld1
            taubeam(k,2)=taucld2
            taubeam(k,3)=taucld3
            taubeam(k,4)=taucld4

            taudiff(k,1)=taucld1
            taudiff(k,2)=taucld2
            taudiff(k,3)=taucld3
            taudiff(k,4)=taucld4
         end if

!-----compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)

         ssaclt=0.99999
         asyclt=1.0
         tauc=taucld1+taucld2+taucld3+taucld4

         if (tauc.gt.0.02 .and. fcld(k).gt.0.01) then

            w1=(1.-(aia_nir(ib,1)+(aia_nir(ib,2)+aia_nir(ib,3)*reff(k,1))*reff(k,1)))*taucld1
            w2=(1.-(awa_nir(ib,1)+(awa_nir(ib,2)+awa_nir(ib,3)*reff(k,2))*reff(k,2)))*taucld2
            w3=(1.- ara_nir(ib,1))                                                   *taucld3
            w4=(1.-(aia_nir(ib,1)+(aia_nir(ib,2)+aia_nir(ib,3)*reff_snow)*reff_snow))*taucld4
            ssaclt=(w1+w2+w3+w4)/tauc

            g1=(aig_nir(ib,1)+(aig_nir(ib,2)+aig_nir(ib,3)*reff(k,1))*reff(k,1))*w1
            g2=(awg_nir(ib,1)+(awg_nir(ib,2)+awg_nir(ib,3)*reff(k,2))*reff(k,2))*w2
            g3= arg_nir(ib,1)                                                   *w3

            g4=(aig_nir(ib,1)+(aig_nir(ib,2)+aig_nir(ib,3)*reff(k,4))*reff(k,4))*w4
	    
	    if ((w1+w2+w3+w4).ne.0.0) then
             asyclt=(g1+g2+g3+g4)/(w1+w2+w3+w4)
	    end if
	     

         end if

         ssacl(k)=ssaclt
         asycl(k)=asyclt

      end do

      return

   end subroutine getnirtau1


end module SORADMOD
