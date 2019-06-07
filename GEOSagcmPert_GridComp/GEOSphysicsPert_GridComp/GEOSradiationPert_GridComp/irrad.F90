module IRRADMOD

IMPLICIT NONE

PRIVATE
PUBLIC :: irrad, planck, plancd, h2oexps, conexps, n2oexps, ch4exps, comexps, cfcexps, b10exps, &
          tablup, h2okdis, n2okdis, ch4kdis, comkdis, cfckdis, b10kdis, cldovlp, sfcflux, mkicx, SORTIT, getirtau1

contains

 subroutine irrad ( np          , & !Number of layers (LM)
                    ple_dev     , & !Pressure at level edges (Pa)
                    ta_dev      , & !d - Temperature (K)
                    wa_dev      , & !d - Specific humidity (g/g)
                    oa_dev      , & !d - Ozone (g/g)
                    tb_dev      , & !Surface air temperature (K)
                    co2         , & !* Carbon dioxide (pppv)
                    trace       , & !Option
                    n2o_dev     , & !* Nitrous oxide (pppv)
                    ch4_dev     , & !* Methane (pppv)
                    cfc11_dev   , & !* Trichlorofluoromethane (pppv)
                    cfc12_dev   , & !* Dichlorodifluoromethane (pppv)
                    cfc22_dev   , & !* Chlorodifluoromethane (pppv)
                    cwc_dev     , & !Cloud water mixing ratio (kg/kg) 
                    fcld_dev    , & !Cloud amount (fraction)
                    ict         , & !Level index separating high and middle clouds
                    icb         , & !Level index separating middle and low clouds
                    reff_dev    , & !Effective size of cloud particles (micron)
                    ns          , & !Number of sub-grid surface types
                    fs_dev      , & !Fractional cover of sub-grid regions
                    tg_dev      , & !Land or ocean surface temperature
                    eg_dev      , & !Land or ocean surface emissivity
                    tv_dev      , & !Vegetation temperature
                    ev_dev      , & !Vegetation emissivity
                    rv_dev      , & !Vegetation reflectivity 
                    na          , & !Number of bands
                    nb          , & !Number of bands in IRRAD calcs for Chou
                    taua_dev    , & !Aerosol optical thickness
                    ssaa_dev    , & !Aerosol single scattering albedo
                    asya_dev    , & !Aerosol asymmetry factor
                    flxu_dev    , & !Upwelling flux, all-sky
                    flxd_dev    , & !Downwelling flux, all-sky
                    dfdts_dev   , & !Sensitivity of net downward flux to surface temperature
                    aib_ir, awb_ir, aiw_ir,   &
                    aww_ir, aig_ir, awg_ir,   &
                    xkw, xke, mw, aw, bw, pm, &
                    fkw, gkw, cb, dcb,        &
                    w11, w12, w13, p11, p12,  &
                    p13, dwe, dpe,            &
                    c1,  c2,  c3,             &
                    oo1, oo2, oo3,            &
                    h11, h12, h13,            &
                    h21, h22, h23,            &
                    h81, h82, h83             )

 IMPLICIT NONE

 !Radiation constants, these need to be inputs for the autodiff tool
 real(8), intent(IN)    :: aib_ir(3,10), awb_ir(4,10), aiw_ir(4,10)
 real(8), intent(IN)    :: aww_ir(4,10), aig_ir(4,10), awg_ir(4,10)

 integer, intent(IN)    :: mw(9)
 real(8), intent(IN)    :: xkw(9), xke(9), aw(9), bw(9), pm(9)
 real(8), intent(IN)    :: fkw(6,9), gkw(6,3), cb(6,10), dcb(5,10)
 real(8), intent(IN)    :: w11, w12, w13, p11, p12
 real(8), intent(IN)    :: p13, dwe, dpe
 real(8), intent(IN)    :: c1(26,30),  c2(26,30),  c3(26,30)
 real(8), intent(IN)    :: oo1(26,21), oo2(26,21), oo3(26,21)
 real(8), intent(IN)    :: h11(26,31), h12(26,31), h13(26,31)
 real(8), intent(IN)    :: h21(26,31), h22(26,31), h23(26,31)
 real(8), intent(IN)    :: h81(26,31), h82(26,31), h83(26,31)

 !----- INPUTS -----
 INTEGER, INTENT(IN)                     :: np, ict, icb, ns, na, nb
 LOGICAL, INTENT(IN)                     :: trace 

 real(8), INTENT(IN)                        :: CO2

 !Rank 2 inputs
 real(8), INTENT(IN)          :: tb_dev

 !Rank 3 (Prognostic variables and tracers)
 real(8), DIMENSION(np),    INTENT(IN)    :: ta_dev, wa_dev, oa_dev, fcld_dev
 real(8), DIMENSION(np),    INTENT(IN)    :: n2o_dev, ch4_dev, cfc11_dev, cfc12_dev, cfc22_dev
 real(8), DIMENSION(np+1),  INTENT(IN)    :: ple_dev

 !Rank 3 (surface types)
 real(8), DIMENSION(ns),    INTENT(IN)    :: fs_dev, tg_dev, tv_dev
 real(8), DIMENSION(ns,10), INTENT(IN)    :: eg_dev, ev_dev, rv_dev

 !Rank 3 (diagnostic cloud parts)
 real(8), DIMENSION(np,4),  INTENT(IN)    :: cwc_dev, reff_dev

 !Rank 3 (aerosols)
 real(8), DIMENSION(np,nb), INTENT(INOUT) :: taua_dev, ssaa_dev, asya_dev


 !----- OUPUTS -----
 real(8), DIMENSION(np+1), INTENT(OUT)    :: flxu_dev
 real(8), DIMENSION(np+1), INTENT(OUT)    :: flxd_dev
 real(8), DIMENSION(np+1), INTENT(OUT)    :: dfdts_dev

 !----- LOCALS -----
 real(8), parameter :: CONS_GRAV = 9.80665

 integer, parameter :: nx1 = 26
 integer, parameter :: no1 = 21
 integer, parameter :: nc1 = 30
 integer, parameter :: nh1 = 31


 !Temporary arrays
 real(8) :: pa(0:np),dt(0:np)
 real(8) :: x1,x2,x3
 real(8) :: dh2o(0:np),dcont(0:np),dco2(0:np),do3(0:np)
 real(8) :: dn2o(0:np),dch4(0:np)
 real(8) :: df11(0:np),df12(0:np),df22(0:np)
 real(8) :: th2o(6),tcon(3),tco2(6)
 real(8) :: tn2o(4),tch4(4),tcom(6)
 real(8) :: tf11,tf12,tf22
 real(8) :: blayer(0:np+1),blevel(0:np+1)
 real(8) :: bd(0:np+1),bu(0:np+1)
 real(8) :: bs,dbs,rflxs
 real(8) :: dp(0:np)
 real(8) :: trant,tranal
 real(8) :: transfc(0:np+1)
 real(8) :: flxu(0:np+1),flxd(0:np+1)
 real(8) :: taerlyr(0:np)

 !OVERCAST
 integer :: ncld(3)
 integer :: icx(0:np)
 !OVERCAST

 integer :: idx, rc
 integer :: k,l,ip,iw,ibn,ik,iq,isb,k1,k2,ne

 real(8) :: enn(0:np)
 real(8) :: cldhi,cldmd,cldlw,tcldlyr(0:np),fclr
 real(8) :: x,xx,yy,p1,a1,b1,fk1,a2,b2,fk2
 real(8) :: w1,ff

 logical :: oznbnd,co2bnd,h2otable,conbnd,n2obnd
 logical :: ch4bnd,combnd,f11bnd,f12bnd,f22bnd,b10bnd
 logical :: do_aerosol 

 !Temp arrays and variables for consolidation of tables
 integer, parameter :: max_num_tables = 17
 real(8) :: exptbl(0:np,max_num_tables)
 type :: band_table
    integer :: start
    integer :: end
 end type band_table
 type(band_table) :: h2oexp
 type(band_table) :: conexp
 type(band_table) :: co2exp
 type(band_table) :: n2oexp
 type(band_table) :: ch4exp
 type(band_table) :: comexp
 type(band_table) :: f11exp
 type(band_table) :: f12exp
 type(band_table) :: f22exp

 !Variables for new getirtau routine
 real(8) :: dp_pa(np)
 real(8) :: fcld_col(np)
 real(8) :: reff_col(np,4)
 real(8) :: cwc_col(np,4)

 real(8) :: h2oexp_tmp(0:np,5),conexp_tmp(0:np),co2exp_tmp(0:np,6),n2oexp_tmp(0:np,2)

 !BEGIN CALCULATIONS ...

      !compute layer pressure (pa) and layer temperature minus 250K (dt) 
      do k=1,np
         pa(k) = 0.5*(ple_dev(k+1)+ple_dev(k))*0.01
         dp(k) =     (ple_dev(k+1)-ple_dev(k))*0.01
         dp_pa(k) =     (ple_dev(k+1)-ple_dev(k)) ! dp in Pascals for getirtau
         dt(k) = ta_dev(k)-250.0

         !compute layer absorber amount

         !dh2o : water vapor amount (g/cm**2)
         !dcont: scaled water vapor amount for continuum absorption
         !       (g/cm**2)
         !dco2 : co2 amount (cm-atm)stp
         !do3  : o3 amount (cm-atm)stp
         !dn2o : n2o amount (cm-atm)stp
         !dch4 : ch4 amount (cm-atm)stp
         !df11 : cfc11 amount (cm-atm)stp
         !df12 : cfc12 amount (cm-atm)stp
         !df22 : cfc22 amount (cm-atm)stp
         !the factor 1.02 is equal to 1000/980
         !factors 789 and 476 are for unit conversion
         !the factor 0.001618 is equal to 1.02/(.622*1013.25) 
         !the factor 6.081 is equal to 1800/296

         dh2o(k) = 1.02*wa_dev   (k)*dp(k)
         do3 (k) = 476.*oa_dev   (k)*dp(k)
         dco2(k) = 789.*co2           *dp(k)
         dch4(k) = 789.*ch4_dev  (k)*dp(k)
         dn2o(k) = 789.*n2o_dev  (k)*dp(k)
         df11(k) = 789.*cfc11_dev(k)*dp(k)
         df12(k) = 789.*cfc12_dev(k)*dp(k)
         df22(k) = 789.*cfc22_dev(k)*dp(k)

         dh2o(k) = max(dh2o(k),1.e-10)
         do3 (k) = max(do3 (k),1.e-6)
         dco2(k) = max(dco2(k),1.e-4)

         !Compute scaled water vapor amount for h2o continuum absorption
         !following eq. (4.21).

         xx=pa(k)*0.001618*wa_dev(k)*wa_dev(k)*dp(k)
         dcont(k) = xx*exp(1800./ta_dev(k)-6.081)

         !Fill the reff, cwc, and fcld for the column   
         fcld_col(k) = fcld_dev(k)
         do l = 1, 4
            reff_col(k,l) = reff_dev(k,l)
            cwc_col(k,l) = cwc_dev(k,l)
         end do

      end do

      !A layer is added above the top of the model atmosphere.
      !Index "0" is the layer above the top of the atmosphere.

      dp  (0) = max(ple_dev(1)*0.01,0.005)
      pa  (0) = 0.5*dp(0)
      dt  (0) = ta_dev(1)-250.0

      dh2o(0) = 1.02*wa_dev   (1)*dp(0)
      do3 (0) = 476.*oa_dev   (1)*dp(0)
      dco2(0) = 789.*co2           *dp(0)
      dch4(0) = 789.*ch4_dev  (1)*dp(0)
      dn2o(0) = 789.*n2o_dev  (1)*dp(0)
      df11(0) = 789.*cfc11_dev(1)*dp(0)
      df12(0) = 789.*cfc12_dev(1)*dp(0)
      df22(0) = 789.*cfc22_dev(1)*dp(0)

      dh2o(0) = max(dh2o(0),1.e-10)
      do3 (0) = max(do3(0),1.e-6)
      dco2(0) = max(dco2(0),1.e-4)

      xx=pa(0)*0.001618*wa_dev(1)*wa_dev(1)*dp(0)
      dcont(0) = xx*exp(1800./ta_dev(1)-6.081)

      !The surface (np+1) is treated as a layer filled with black clouds.
      !transfc is the transmittance between the surface and a pressure
      !level.
      transfc(np+1)=1.0

      !Initialize fluxes
      do k=1,np+1
         flxu_dev(k)  = 0.0
         flxd_dev(k)  = 0.0
         dfdts_dev(k)= 0.0
      end do

      !Integration over spectral bands
      do ibn=1,10

         if (ibn == 10 .and. .not. trace) return

         !if h2otable, compute h2o (line) transmittance using table look-up.
         !if conbnd,   compute h2o (continuum) transmittance in bands 2-7.
         !if co2bnd,   compute co2 transmittance in band 3.
         !if oznbnd,   compute  o3 transmittance in band 5.
         !if n2obnd,   compute n2o transmittance in bands 6 and 7.
         !if ch4bnd,   compute ch4 transmittance in bands 6 and 7.
         !if combnd,   compute co2-minor transmittance in bands 4 and 5.
         !if f11bnd,   compute cfc11 transmittance in bands 4 and 5.
         !if f12bnd,   compute cfc12 transmittance in bands 4 and 6.
         !if f22bnd,   compute cfc22 transmittance in bands 4 and 6.
         !if b10bnd,   compute flux reduction due to n2o in band 10.

         h2otable=ibn == 1 .or. ibn == 2 .or. ibn == 8
         conbnd  =ibn >= 2 .and. ibn <= 7
         co2bnd  =ibn == 3
         oznbnd  =ibn == 5
         n2obnd  =ibn == 6 .or. ibn == 7
         ch4bnd  =ibn == 6 .or. ibn == 7
         combnd  =ibn == 4 .or. ibn == 5
         f11bnd  =ibn == 4 .or. ibn == 5
         f12bnd  =ibn == 4 .or. ibn == 6
         f22bnd  =ibn == 4 .or. ibn == 6
         b10bnd  =ibn == 10

         do_aerosol = na > 0

         exptbl = 0.0

         !Control packing of the new exponential tables by band

         select case (ibn)
         case (2)
            conexp%start = 1
            conexp%end   = 1
         case (3)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 9
         case (4)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 7
            comexp%start = 8
            comexp%end   = 13
            f11exp%start = 14
            f11exp%end   = 14
            f12exp%start = 15
            f12exp%end   = 15
            f22exp%start = 16
            f22exp%end   = 16
         case (5)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 7
            comexp%start = 8
            comexp%end   = 13
            f11exp%start = 14
            f11exp%end   = 14
         case (6)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 7
            n2oexp%start = 8
            n2oexp%end   = 11
            ch4exp%start = 12
            ch4exp%end   = 15
            f12exp%start = 16
            f12exp%end   = 16
            f22exp%start = 17
            f22exp%end   = 17
         case (7)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 7
            n2oexp%start = 8
            n2oexp%end   = 11
            ch4exp%start = 12
            ch4exp%end   = 15
         case (9)
            h2oexp%start = 1
            h2oexp%end   = 6
         case (10)
            h2oexp%start = 1
            h2oexp%end   = 5
            conexp%start = 6
            conexp%end   = 6
            co2exp%start = 7
            co2exp%end   = 12
            n2oexp%start = 13
            n2oexp%end   = 14
         end select

         !blayer is the spectrally integrated planck flux of the mean layer
         !temperature derived from eq. (3.11)
         !The fitting for the planck flux is valid for the range 160-345 K.

         do k=1,np
            call planck(ibn,cb,ta_dev(k),blayer(k))
         end do

         !Index "0" is the layer above the top of the atmosphere.
         blayer(0)=blayer(1)
         blevel(0)=blayer(1)

         !Surface emission and reflectivity. See Section 9.
         !bs and dbs include the effect of surface emissivity.
         call sfcflux (ibn,cb,dcb,ns,fs_dev,tg_dev,eg_dev,tv_dev,ev_dev,rv_dev,&
               bs,dbs,rflxs) 

         blayer(np+1)=bs

         !interpolate Planck function at model levels (linear in p)

         do k=2,np
            blevel(k)=(blayer(k-1)*dp(k)+blayer(k)*dp(k-1))/(dp(k-1)+dp(k))
         end do

         !Extrapolate blevel(1) from blayer(2) and blayer(1)

         blevel(1)=blayer(1)+(blayer(1)-blayer(2))*dp(1)/(dp(1)+dp(2))
         blevel(0)=blevel(1)

         !If the surface air temperature tb is known, compute blevel(np+1)
         call planck(ibn,cb,tb_dev,blevel(np+1))

         !Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
         !     NOTE: dp_pa is only dims(1:np) as the 0'th level isn't needed in getirtau.
         !           Plus, the pressures in getirtau *MUST* be in Pascals.
         !     Slots for reff, hydrometeors and tauall are as follows:
         !                 1         Cloud Ice
         !                 2         Cloud Liquid
         !                 3         Falling Liquid (Rain)
         !                 4         Falling Ice (Snow)

         call getirtau1(ibn,np,dp_pa,fcld_col,reff_col,cwc_col,&
                       tcldlyr,enn,aib_ir,awb_ir, &
                       aiw_ir, aww_ir, aig_ir, awg_ir, CONS_GRAV)

         do k=0,np
            icx(k) = k
         end do

         call mkicx(np,ict,icb,enn,icx,ncld)

         !Compute optical thickness, single-scattering albedo and asymmetry
         !factor for a mixture of "na" aerosol types. Eqs. (7.1)-(7.3)
         if (do_aerosol) then
            taerlyr(0)=1.0

            do k=1,np

              !-----taerlyr is the aerosol diffuse transmittance
               taerlyr(k)=1.0
               if (taua_dev(k,ibn) > 0.001) then 
                  if (ssaa_dev(k,ibn) > 0.001) then
                     asya_dev(k,ibn)=asya_dev(k,ibn)/ssaa_dev(k,ibn)
                     ssaa_dev(k,ibn)=ssaa_dev(k,ibn)/taua_dev(k,ibn)

                     !Parameterization of aerosol scattering following
                     ff=.5+(.3739+(0.0076+0.1185*asya_dev(k,ibn))*asya_dev(k,ibn))*asya_dev(k,ibn)
                     taua_dev(k,ibn)=taua_dev(k,ibn)*(1.-ssaa_dev(k,ibn)*ff)
                  end if
                  taerlyr(k)=exp(-1.66*taua_dev(k,ibn))
               end if
            end do
         end if

         !Compute the exponential terms (Eq. 8.21) at each layer due to
         !water vapor line absorption when k-distribution is used
         if (.not. h2otable .and. .not. b10bnd) then
            call h2oexps(ibn,np,dh2o,pa,dt,xkw,aw,bw,pm,mw,exptbl(:,h2oexp%start:h2oexp%end))
         end if

         !compute the exponential terms (Eq. 4.24) at each layer due to
         !water vapor continuum absorption.

         ne=0
         if (conbnd) then
            ne=1
            if (ibn == 3) ne=3
            call conexps(ibn,np,dcont,xke,exptbl(:,conexp%start:conexp%end))
         end if


         if (trace) then

            !compute the exponential terms at each layer due to n2o absorption
            if (n2obnd) then
               call n2oexps(ibn,np,dn2o,pa,dt,exptbl(:,n2oexp%start:n2oexp%end))
            end if

            !compute the exponential terms at each layer due to ch4 absorption
            if (ch4bnd) then
               call ch4exps(ibn,np,dch4,pa,dt,exptbl(:,ch4exp%start:ch4exp%end))
            end if

            !Compute the exponential terms due to co2 minor absorption
            if (combnd) then
               call comexps(ibn,np,dco2,dt,exptbl(:,comexp%start:comexp%end))
            end if

            !Compute the exponential terms due to cfc11 absorption.
            !The values of the parameters are given in Table 7.
            if (f11bnd) then
               a1  = 1.26610e-3
               b1  = 3.55940e-6
               fk1 = 1.89736e+1
               a2  = 8.19370e-4
               b2  = 4.67810e-6
               fk2 = 1.01487e+1
               call cfcexps(ibn,np,a1,b1,fk1,a2,b2,fk2,df11,dt,&
                     exptbl(:,f11exp%start:f11exp%end))
            end if

            !Compute the exponential terms due to cfc12 absorption.
            if (f12bnd) then
               a1  = 8.77370e-4
               b1  =-5.88440e-6
               fk1 = 1.58104e+1
               a2  = 8.62000e-4
               b2  =-4.22500e-6
               fk2 = 3.70107e+1
               call cfcexps(ibn,np,a1,b1,fk1,a2,b2,fk2,df12,dt,&
                     exptbl(:,f12exp%start:f12exp%end))
            end if

            !Compute the exponential terms due to cfc22 absorption.
            if (f22bnd) then
               a1  = 9.65130e-4
               b1  = 1.31280e-5
               fk1 = 6.18536e+0
               a2  =-3.00010e-5 
               b2  = 5.25010e-7
               fk2 = 3.27912e+1
               call cfcexps(ibn,np,a1,b1,fk1,a2,b2,fk2,df22,dt,&
                     exptbl(:,f22exp%start:f22exp%end))
            end if

            !Compute the exponential terms at each layer in band 10 due to
            !h2o line and continuum, co2, and n2o absorption

            if (b10bnd) then
               call b10exps(np,dh2o,dcont,dco2,dn2o,pa,dt, &
                     h2oexp_tmp,exptbl(:,conexp%start:conexp%end),co2exp_tmp,n2oexp_tmp)

               exptbl(:,h2oexp%start:h2oexp%end) = h2oexp_tmp
               exptbl(:,co2exp%start:co2exp%end) = co2exp_tmp
               exptbl(:,n2oexp%start:n2oexp%end) = n2oexp_tmp

            end if
         end if

         !blayer(np+1) includes the effect of surface emissivity.

         bu(0)=0.0 ! ALT: this was undefined, check with Max if 0.0 is good value
         bd(0)=blayer(1)
         bu(np+1)=blayer(np+1)

         !do-loop 1500 is for computing upward (bu) and downward (bd)
         !Here, trant is the transmittance of the layer k2-1.
         do k2=1,np+1

            !for h2o line transmission
            if (.not. h2otable) then
               th2o=1.0
            end if

            !for h2o continuum transmission
            tcon=1.0

            x1=0.0
            x2=0.0
            x3=0.0
            trant=1.0

            if (h2otable) then

               !Compute water vapor transmittance using table look-up.
               if (ibn == 1) then
                  call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w11,p11,dwe,dpe,h11,h12,h13,trant)
               end if
               if (ibn == 2) then
                  call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w11,p11,dwe,dpe,h21,h22,h23,trant)
               end if
               if (ibn == 8) then
                  call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w11,p11,dwe,dpe,h81,h82,h83,trant)
               end if

               !for water vapor continuum absorption
               if (conbnd) then
                  tcon(1)=tcon(1)*exptbl(k2-1,conexp%start) ! Only the first exp
                  trant=trant*tcon(1)
               end if

            else

               !compute water vapor transmittance using k-distribution
               if (.not. b10bnd) then
                  call h2okdis(ibn,np,k2-1,fkw,gkw,ne,&
                        exptbl(:,h2oexp%start:h2oexp%end), &
                        exptbl(:,conexp%start:conexp%end), &
                        th2o,tcon,trant)
               end if

            end if

            if (co2bnd) then
               !Compute co2 transmittance using table look-up method
               call tablup(nx1,nc1,dco2(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                     w12,p12,dwe,dpe,c1,c2,c3,trant)
            end if

            !Always use table look-up to compute o3 transmittance.
            if (oznbnd) then
               call tablup(nx1,no1,do3(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                     w13,p13,dwe,dpe,oo1,oo2,oo3,trant)
            end if

            !include aerosol effect
            if (do_aerosol) then
               trant=trant*taerlyr(k2-1)
            end if

            !Compute upward (bu) and downward (bd) emission of the layer k2-1
            xx=(1.-enn(k2-1))*trant
            yy=min(0.9999,xx)
            yy=max(0.00001,yy)
            xx=(blevel(k2-1)-blevel(k2))/log(yy)
            bd(k2-1)=(blevel(k2)-blevel(k2-1)*yy)/(1.0-yy)-xx
            bu(k2-1)=(blevel(k2-1)+blevel(k2))-bd(k2-1)

         end do

         !Initialize fluxes
         flxu = 0.0
         flxd = 0.0

         !Compute upward and downward fluxes for each spectral band, ibn.
         do k1=0,np

            !initialization
            cldlw = 0.0
            cldmd = 0.0
            cldhi = 0.0
            tranal= 1.0

            !for h2o line transmission
            if (.not. h2otable) then
               th2o=1.0
            end if

            !for h2o continuum transmission
            tcon=1.0

            if (trace) then !Add trace gases contribution

               if (n2obnd) then !n2o
                  tn2o=1.0
               end if

               if (ch4bnd) then !ch4
                  tch4=1.0
               end if

               if (combnd) then !co2-minor
                  tcom=1.0
               end if

               if (f11bnd) then !cfc-11
                  tf11=1.0
               end if

               if (f12bnd) then !cfc-12
                  tf12=1.0
               end if

               if (f22bnd) then !cfc-22
                  tf22=1.0
               end if

               if (b10bnd) then !
                  th2o=1.0

                  tco2=1.0

                  tcon(1)=1.0

                  tn2o=1.0
               end if

            end if

            x1=0.0
            x2=0.0
            x3=0.0

            do k2=k1+1,np+1

               trant=1.0
               fclr =1.0

               if (h2otable) then

                  !Compute water vapor transmittance using table look-up.
                  if (ibn == 1) then
                     call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                           w11,p11,dwe,dpe,h11,h12,h13,trant)
                  end if
                  if (ibn == 2) then
                     call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                           w11,p11,dwe,dpe,h21,h22,h23,trant)
                  end if
                  if (ibn == 8) then
                     call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                           w11,p11,dwe,dpe,h81,h82,h83,trant)
                  end if

                  if (conbnd) then
                     tcon(1)=tcon(1)*exptbl(k2-1,conexp%start) ! Only the first exp
                     trant=trant*tcon(1)
                  end if

               else

                  !Compute water vapor transmittance using k-distribution
                  if (.not. b10bnd) then
                     call h2okdis(ibn,np,k2-1,fkw,gkw,ne,&
                           exptbl(:,h2oexp%start:h2oexp%end), &
                           exptbl(:,conexp%start:conexp%end), &
                           th2o,tcon,trant)
                  end if

               end if

               if (co2bnd) then

                  !Compute co2 transmittance using table look-up method.
                  call tablup(nx1,nc1,dco2(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w12,p12,dwe,dpe,c1,c2,c3,trant)
               end if

               if (oznbnd) then
                  !Always use table look-up to compute o3 transmittanc
                  call tablup(nx1,no1,do3(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w13,p13,dwe,dpe,oo1,oo2,oo3,trant)
               end if

               if (trace) then !Add trace gas effects

                  if (n2obnd) then !n2o
                     call n2okdis(ibn,np,k2-1,exptbl(:,n2oexp%start:n2oexp%end),tn2o,trant)
                  end if

                  if (ch4bnd) then !ch4
                     call ch4kdis(ibn,np,k2-1,exptbl(:,ch4exp%start:ch4exp%end),tch4,trant)
                  end if

                  if (combnd) then !co2-minor
                     call comkdis(ibn,np,k2-1,exptbl(:,comexp%start:comexp%end),tcom,trant)
                  end if

                  if (f11bnd) then !cfc11
                     call cfckdis(np,k2-1,exptbl(:,f11exp%start:f11exp%end),tf11,trant)
                  end if

                  if (f12bnd) then !cfc12
                     call cfckdis(np,k2-1,exptbl(:,f12exp%start:f12exp%end),tf12,trant)
                  end if

                  if (f22bnd) then !CFC22
                     call cfckdis(np,k2-1,exptbl(:,f22exp%start:f22exp%end),tf22,trant)
                  end if

                  if (b10bnd) then
                     call b10kdis(np,k2-1,&
                           exptbl(:,h2oexp%start:h2oexp%end),&
                           exptbl(:,conexp%start:conexp%end),&
                           exptbl(:,co2exp%start:co2exp%end),&
                           exptbl(:,n2oexp%start:n2oexp%end),&
                           th2o,tcon,tco2,tn2o,trant)
                  end if

               end if


               if (do_aerosol) then
                  tranal=tranal*taerlyr(k2-1)
                  trant=trant *tranal
               end if

               if (enn(k2-1) >= 0.001) then
                  call cldovlp (np,k1,k2,ict,icb,icx,ncld,enn,tcldlyr,cldhi,cldmd,cldlw)
               end if

               fclr=(1.0-cldhi)*(1.0-cldmd)*(1.0-cldlw)

               if (k2 == k1+1 .and. ibn /= 10) then
                  flxu(k1)=flxu(k1)-bu(k1)
                  flxd(k2)=flxd(k2)+bd(k1)
               end if

               xx=trant*(bu(k2-1)-bu(k2))
               flxu(k1)=flxu(k1)+xx*fclr

               if (k1 == 0) then !mjs  bd(-1) is not defined
                  xx=-trant*bd(k1)
               else
                  xx= trant*(bd(k1-1)-bd(k1))
               end if
               flxd(k2)=flxd(k2)+xx*fclr

            end do

            transfc(k1) =trant*fclr

            if (k1 > 0)then
               dfdts_dev(k1) =dfdts_dev(k1)-dbs*transfc(k1)
            end if

         end do

         if (.not. b10bnd) then !surface emission

            flxu(np+1)       =             -blayer(np+1)
            dfdts_dev(np+1)= dfdts_dev(np+1)-dbs

            do k=1,np+1
               flxu(k) = flxu(k)-flxd(np+1)*transfc(k)*rflxs
            end do

         end if

         do k=1,np+1 !Summation of fluxes over spectral bands

            flxu_dev(k) = flxu_dev(k) + flxu(k)
            flxd_dev(k) = flxd_dev(k) + flxd(k)

         end do

      end do

   end subroutine irrad

!***********************************************************************
   subroutine planck(ibn,cb,t,xlayer)
!***********************************************************************
!
!-----Compute spectrally integrated Planck flux
!
   implicit none

   integer, intent(in) :: ibn                   ! spectral band index
   real(8), intent(in) :: cb(6,10)                 ! Planck table coefficients
   real(8), intent(in) :: t                        ! temperature (K)
   real(8), intent(out) :: xlayer                   ! planck flux (w/m2)

   xlayer=t*(t*(t*(t*(t*cb(6,ibn)+cb(5,ibn))+cb(4,ibn))+cb(3,ibn))+cb(2,ibn))+cb(1,ibn)

   end subroutine planck


!***********************************************************************
   subroutine plancd(ibn,dcb,t,dbdt) 
!***********************************************************************
!
!-----Compute the derivative of Planck flux wrt temperature
!
   implicit none

   integer, intent(in) :: ibn               ! spectral band index
   real(8), intent(in) :: dcb(5,10)            ! Planck function derivative coefficients
   real(8), intent(in) :: t                    ! temperature (K)
   real(8), intent(out) :: dbdt                 ! derivative of Planck flux wrt temperature

   dbdt=t*(t*(t*(t*dcb(5,ibn)+dcb(4,ibn))+dcb(3,ibn))+dcb(2,ibn))+dcb(1,ibn)

   end subroutine plancd


!**********************************************************************
   subroutine h2oexps(ib,np,dh2o,pa,dt,xkw,aw,bw,pm,mw,h2oexp)
!**********************************************************************
!   Compute exponentials for water vapor line absorption
!   in individual layers using Eqs. (8.21) and (8.22).
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer water vapor amount for line absorption (dh2o) 
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!  absorption coefficients for the first k-distribution
!     function due to h2o line absorption (xkw)
!  coefficients for the temperature and pressure scaling (aw,bw,pm)
!  ratios between neighboring absorption coefficients for
!     h2o line absorption (mw)
!
!---- output parameters
!  6 exponentials for each layer  (h2oexp)
!**********************************************************************
   implicit none

   integer :: ib,np,ik,k

!---- input parameters ------

   real(8) :: dh2o(0:np),pa(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: h2oexp(0:np,6)

!---- static data -----

   integer :: mw(9)
   real(8) :: xkw(9),aw(9),bw(9),pm(9)

!---- temporary arrays -----

   real(8) :: xh

!**********************************************************************
!    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
!    and bw,  therefore, h2oexp for these sub-bands are identical.
!**********************************************************************

   do k=0,np

!-----xh is the scaled water vapor amount for line absorption
!     computed from Eq. (4.4).

      xh = dh2o(k)*(pa(k)/500.)**pm(ib) &
            * ( 1.+(aw(ib)+bw(ib)* dt(k))*dt(k) )

!-----h2oexp is the water vapor transmittance of the layer k
!     due to line absorption

      h2oexp(k,1) = exp(-xh*xkw(ib))

!-----compute transmittances from Eq. (8.22)

      do ik=2,6
         if (mw(ib) == 6) then
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)
            h2oexp(k,ik) = xh*xh*xh
         elseif (mw(ib) == 8) then
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)
            xh = xh*xh
            h2oexp(k,ik) = xh*xh
         elseif (mw(ib) == 9) then
            xh=h2oexp(k,ik-1)*h2oexp(k,ik-1)*h2oexp(k,ik-1)
            h2oexp(k,ik) = xh*xh*xh
         else
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)
            xh = xh*xh
            xh = xh*xh
            h2oexp(k,ik) = xh*xh
         end if
      end do
   end do

   end subroutine h2oexps


!**********************************************************************
   subroutine conexps(ib,np,dcont,xke,conexp)
!**********************************************************************
!   compute exponentials for continuum absorption in individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer scaled water vapor amount for continuum absorption (dcont) 
!  absorption coefficients for the first k-distribution function
!     due to water vapor continuum absorption (xke)
!
!---- output parameters
!  1 or 3 exponentials for each layer (conexp)
!**********************************************************************
   implicit none

   integer :: ib,np,k

!---- input parameters ------

   real(8) :: dcont(0:np)

!---- updated parameters -----

   real(8) :: conexp(0:np,3)

!---- static data -----

   real(8) :: xke(9)

!****************************************************************

   do k=0,np

      conexp(k,1) = exp(-dcont(k)*xke(ib))

!-----The absorption coefficients for sub-bands 3b and 3a are, respectively,
!     two and four times the absorption coefficient for sub-band 3c (Table 9).
!     Note that conexp(3) is for sub-band 3a. 

      if (ib  ==  3) then
         conexp(k,2) = conexp(k,1) *conexp(k,1)
         conexp(k,3) = conexp(k,2) *conexp(k,2)
      end if

   end do

   end subroutine conexps


!**********************************************************************
   subroutine n2oexps(ib,np,dn2o,pa,dt,n2oexp)
!**********************************************************************
!   Compute n2o exponentials for individual layers 
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer n2o amount (dn2o)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  2 or 4 exponentials for each layer (n2oexp)
!**********************************************************************
   implicit none

   integer :: ib,k,np

!---- input parameters -----

   real(8) :: dn2o(0:np),pa(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: n2oexp(0:np,4)

!---- temporary arrays -----

   real(8) :: xc,xc1,xc2

!-----Scaling and absorption data are given in Table 5.
!     Transmittances are computed using Eqs. (8.21) and (8.22).

   do k=0,np

!-----four exponential by powers of 21 for band 6.

      if (ib == 6) then

         xc=dn2o(k)*(1.+(1.9297e-3+4.3750e-6*dt(k))*dt(k))
         n2oexp(k,1)=exp(-xc*6.31582e-2)

         xc=n2oexp(k,1)*n2oexp(k,1)*n2oexp(k,1)
         xc1=xc*xc
         xc2=xc1*xc1
         n2oexp(k,2)=xc*xc1*xc2

!-----four exponential by powers of 8 for band 7

      else

         xc=dn2o(k)*(pa(k)/500.0)**0.48 &
               *(1.+(1.3804e-3+7.4838e-6*dt(k))*dt(k))
         n2oexp(k,1)=exp(-xc*5.35779e-2)

         xc=n2oexp(k,1)*n2oexp(k,1)
         xc=xc*xc
         n2oexp(k,2)=xc*xc
         xc=n2oexp(k,2)*n2oexp(k,2)
         xc=xc*xc
         n2oexp(k,3)=xc*xc
         xc=n2oexp(k,3)*n2oexp(k,3)
         xc=xc*xc
         n2oexp(k,4)=xc*xc
      end if

   end do

   end subroutine n2oexps


!**********************************************************************
   subroutine ch4exps(ib,np,dch4,pa,dt,ch4exp)
!**********************************************************************
!   Compute ch4 exponentials for individual layers
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer ch4 amount (dch4)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  1 or 4 exponentials for each layer (ch4exp)
!**********************************************************************
   implicit none

   integer :: ib,np,k

!---- input parameters -----

   real(8) :: dch4(0:np),pa(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: ch4exp(0:np,4)

!---- temporary arrays -----

   real(8) :: xc

!*****  Scaling and absorption data are given in Table 5  *****

   do k=0,np

!-----four exponentials for band 6

      if (ib == 6) then

         xc=dch4(k)*(1.+(1.7007e-2+1.5826e-4*dt(k))*dt(k))
         ch4exp(k,1)=exp(-xc*5.80708e-3)

!-----four exponentials by powers of 12 for band 7

      else

         xc=dch4(k)*(pa(k)/500.0)**0.65 &
               *(1.+(5.9590e-4-2.2931e-6*dt(k))*dt(k))
         ch4exp(k,1)=exp(-xc*6.29247e-2)

         xc=ch4exp(k,1)*ch4exp(k,1)*ch4exp(k,1)
         xc=xc*xc
         ch4exp(k,2)=xc*xc

         xc=ch4exp(k,2)*ch4exp(k,2)*ch4exp(k,2)
         xc=xc*xc
         ch4exp(k,3)=xc*xc

         xc=ch4exp(k,3)*ch4exp(k,3)*ch4exp(k,3)
         xc=xc*xc
         ch4exp(k,4)=xc*xc

      end if

   end do

   end subroutine ch4exps


!**********************************************************************
   subroutine comexps(ib,np,dcom,dt,comexp)
!**********************************************************************
!   Compute co2-minor exponentials for individual layers using 
!   Eqs. (8.21) and (8.22).
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer co2 amount (dcom)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  6 exponentials for each layer (comexp)
!**********************************************************************
   implicit none

   integer :: ib,ik,np,k

!---- input parameters -----

   real(8) :: dcom(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: comexp(0:np,6)

!---- temporary arrays -----

   real(8) :: xc

!*****  Scaling and absorpton data are given in Table 6  *****

   do k=0,np

      if (ib == 4) then
         xc=dcom(k)*(1.+(3.5775e-2+4.0447e-4*dt(k))*dt(k))
      end if

      if (ib == 5) then
         xc=dcom(k)*(1.+(3.4268e-2+3.7401e-4*dt(k))*dt(k))
      end if

      comexp(k,1)=exp(-xc*1.922e-7)

      do ik=2,6
         xc=comexp(k,ik-1)*comexp(k,ik-1)
         xc=xc*xc
         comexp(k,ik)=xc*comexp(k,ik-1)
      end do

   end do

   end subroutine comexps


!**********************************************************************
   subroutine cfcexps(ib,np,a1,b1,fk1,a2,b2,fk2,dcfc,dt,cfcexp)
!**********************************************************************
!   compute cfc(-11, -12, -22) exponentials for individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  parameters for computing the scaled cfc amounts
!             for temperature scaling (a1,b1,a2,b2)
!  the absorption coefficients for the
!     first k-distribution function due to cfcs (fk1,fk2)
!  layer cfc amounts (dcfc)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  1 exponential for each layer (cfcexp)
!**********************************************************************
   implicit none

   integer :: ib,np,k

!---- input parameters -----

   real(8) :: dcfc(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: cfcexp(0:np)

!---- static data -----

   real(8) :: a1,b1,fk1,a2,b2,fk2

!---- temporary arrays -----

   real(8) :: xf

!**********************************************************************

   do k=0,np

!-----compute the scaled cfc amount (xf) and exponential (cfcexp)

      if (ib == 4) then
         xf=dcfc(k)*(1.+(a1+b1*dt(k))*dt(k))
         cfcexp(k)=exp(-xf*fk1)
      else
         xf=dcfc(k)*(1.+(a2+b2*dt(k))*dt(k))
         cfcexp(k)=exp(-xf*fk2)
      end if

   end do

   end subroutine cfcexps


!**********************************************************************
   subroutine b10exps(np,dh2o,dcont,dco2,dn2o,pa,dt, &
      h2oexp,conexp,co2exp,n2oexp)
!**********************************************************************
!   Compute band3a exponentials for individual layers
!
!---- input parameters
!  number of layers (np)
!  layer h2o amount for line absorption (dh2o)
!  layer h2o amount for continuum absorption (dcont)
!  layer co2 amount (dco2)
!  layer n2o amount (dn2o)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!
!  exponentials for each layer (h2oexp,conexp,co2exp,n2oexp)
!**********************************************************************
   implicit none

   integer :: np,k

!---- input parameters -----

   real(8) :: dh2o(0:np),dcont(0:np),dn2o(0:np)
   real(8) :: dco2(0:np),pa(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: h2oexp(0:np,5),conexp(0:np),co2exp(0:np,6),n2oexp(0:np,2)

!---- temporary arrays -----

   real(8) :: xx,xx1,xx2,xx3

!**********************************************************************

   do k=0,np

!-----Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3).

      xx=dh2o(k)*(pa(k)/500.0) &
            *(1.+(0.0149+6.20e-5*dt(k))*dt(k))

!-----six exponentials by powers of 8

      h2oexp(k,1)=exp(-xx*0.10624)

      xx=h2oexp(k,1)*h2oexp(k,1)
      xx=xx*xx
      h2oexp(k,2)=xx*xx

      xx=h2oexp(k,2)*h2oexp(k,2)
      xx=xx*xx
      h2oexp(k,3)=xx*xx

      xx=h2oexp(k,3)*h2oexp(k,3)
      xx=xx*xx
      h2oexp(k,4)=xx*xx

      xx=h2oexp(k,4)*h2oexp(k,4)
      xx=xx*xx
      h2oexp(k,5)=xx*xx

!-----one exponential of h2o continuum for sub-band 3a (Table 9).

      conexp(k)=exp(-dcont(k)*109.0)

!-----Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).

      xx=dco2(k)*(pa(k)/300.0)**0.5 &
            *(1.+(0.0179+1.02e-4*dt(k))*dt(k))

!-----six exponentials by powers of 8

      co2exp(k,1)=exp(-xx*2.656e-5)

      xx=co2exp(k,1)*co2exp(k,1)
      xx=xx*xx
      co2exp(k,2)=xx*xx

      xx=co2exp(k,2)*co2exp(k,2)
      xx=xx*xx
      co2exp(k,3)=xx*xx

      xx=co2exp(k,3)*co2exp(k,3)
      xx=xx*xx
      co2exp(k,4)=xx*xx

      xx=co2exp(k,4)*co2exp(k,4)
      xx=xx*xx
      co2exp(k,5)=xx*xx

      xx=co2exp(k,5)*co2exp(k,5)
      xx=xx*xx
      co2exp(k,6)=xx*xx

!-----Compute the scaled n2o amount for Band 10 (Table 5).

      xx=dn2o(k)*(1.+(1.4476e-3+3.6656e-6*dt(k))*dt(k))

!-----Two exponentials by powers of 58

      n2oexp(k,1)=exp(-xx*0.25238)

      xx=n2oexp(k,1)*n2oexp(k,1)
      xx1=xx*xx
      xx1=xx1*xx1
      xx2=xx1*xx1
      xx3=xx2*xx2
      n2oexp(k,2)=xx*xx1*xx2*xx3

   end do

   end subroutine b10exps


!**********************************************************************
   subroutine tablup(nx1,nh1,dw,p,dt,s1,s2,s3,w1,p1, &
      dwe,dpe,coef1,coef2,coef3,tran)
!**********************************************************************
!   Compute water vapor, co2 and o3 transmittances between level
!   k1 and and level k2 , using table look-up.
!
!   Calculations follow Eq. (4.16).
!
!---- input ---------------------
!
!  number of pressure intervals in the table (nx1)
!  number of absorber amount intervals in the table (nh1)
!  layer absorber amount (dw)
!  layer pressure in mb (p)
!  deviation of layer temperature from 250K (dt)
!  first value of absorber amount (log10) in the table (w1) 
!  first value of pressure (log10) in the table (p1) 
!  size of the interval of absorber amount (log10) in the table (dwe)
!  size of the interval of pressure (log10) in the table (dpe)
!  pre-computed coefficients (coef1, coef2, and coef3)
!
!---- updated ---------------------
!
!  column integrated absorber amount (s1)
!  absorber-weighted column pressure (s2)
!  absorber-weighted column temperature (s3)
!  transmittance (tran)
!
!  Note: Units of s1 are g/cm**2 for water vapor and
!       (cm-atm)stp for co2 and o3.
!   
!**********************************************************************
   implicit none

   integer :: nx1,nh1

!---- input parameters -----

   real(8) :: w1,p1,dwe,dpe
   real(8) :: dw,p,dt
   real(8) :: coef1(nx1,nh1),coef2(nx1,nh1),coef3(nx1,nh1)

!---- update parameter -----

   real(8) :: s1,s2,s3,tran

!---- temporary variables -----

   real(8) :: we,pe,fw,fp,pa,pb,pc,ax,ba,bb,t1,ca,cb,t2
   real(8) :: x1,x2,x3,xx, x1c
   integer :: iw,ip

!-----Compute effective pressure (x2) and temperature (x3) following 
!     Eqs. (8.28) and (8.29)

   s1=s1+dw
   s2=s2+p*dw
   s3=s3+dt*dw

   x1=s1
   x1c=1.0/s1
   x2=s2*x1c
   x3=s3*x1c

!-----normalize we and pe

!       we=(log10(x1)-w1)/dwe
!       pe=(log10(x2)-p1)/dpe
   we=(log10(x1)-w1)*dwe
   pe=(log10(x2)-p1)*dpe

!-----restrict the magnitudes of the normalized we and pe.

   we=min(we,real(nh1-1))
   pe=min(pe,real(nx1-1))

!-----assign iw and ip and compute the distance of we and pe 
!     from iw and ip.

   iw=int(we+1.0)
   iw=min(iw,nh1-1)
   iw=max(iw, 2)
   fw=we-real(iw-1)

   ip=int(pe+1.0)
   ip=min(ip,nx1-1)
   ip=max(ip, 1)
   fp=pe-real(ip-1)

!-----linear interpolation in pressure

   pa = coef1(ip,iw-1)+(coef1(ip+1,iw-1)-coef1(ip,iw-1))*fp
   pb = coef1(ip,  iw)+(coef1(ip+1,  iw)-coef1(ip,  iw))*fp
   pc = coef1(ip,iw+1)+(coef1(ip+1,iw+1)-coef1(ip,iw+1))*fp

!-----quadratic interpolation in absorber amount for coef1

   ax = ( (pc+pa)*fw + (pc-pa) ) *fw*0.5 + pb*(1.-fw*fw)

!-----linear interpolation in absorber amount for coef2 and coef3

   ba = coef2(ip,  iw)+(coef2(ip+1,  iw)-coef2(ip,  iw))*fp
   bb = coef2(ip,iw+1)+(coef2(ip+1,iw+1)-coef2(ip,iw+1))*fp
   t1 = ba+(bb-ba)*fw

   ca = coef3(ip,  iw)+(coef3(ip+1,  iw)-coef3(ip,  iw))*fp
   cb = coef3(ip,iw+1)+(coef3(ip+1,iw+1)-coef3(ip,iw+1))*fp
   t2 = ca + (cb-ca)*fw

!-----update the total transmittance between levels k1 and k2

   xx=(ax + (t1+t2*x3) * x3)
   xx= min(xx,0.9999999)
   xx= max(xx,0.0000001)
   tran= tran*xx

   end subroutine tablup


!**********************************************************************
   subroutine h2okdis(ib,np,k,fkw,gkw,ne,h2oexp,conexp, &
      th2o,tcon,tran)
!**********************************************************************
!   compute water vapor transmittance between levels k1 and k2, using the k-distribution method.
!
!---- input parameters
!  spectral band (ib)
!  number of levels (np)
!  current level (k)
!  planck-weighted k-distribution function due to
!    h2o line absorption (fkw)
!  planck-weighted k-distribution function due to
!    h2o continuum absorption (gkw)
!  number of terms used in each band to compute water vapor
!     continuum transmittance (ne)
!  exponentials for line absorption (h2oexp) 
!  exponentials for continuum absorption (conexp) 
!
!---- updated parameters
!  transmittance between levels k1 and k2 due to
!    water vapor line absorption (th2o)
!  transmittance between levels k1 and k2 due to
!    water vapor continuum absorption (tcon)
!  total transmittance (tran)
!
!**********************************************************************
   implicit none

!---- input parameters ------

   integer :: ib,ne,np,k
   real(8) :: h2oexp(0:np,6),conexp(0:np,3)
   real(8) ::  fkw(6,9),gkw(6,3)

!---- updated parameters -----

   real(8) :: th2o(6),tcon(3),tran

!---- temporary arrays -----

   real(8) :: trnth2o

!-----tco2 are the six exp factors between levels k1 and k2 
!     tran is the updated total transmittance between levels k1 and k2

!-----th2o is the 6 exp factors between levels k1 and k2 due to
!     h2o line absorption. 

!-----tcon is the 3 exp factors between levels k1 and k2 due to
!     h2o continuum absorption.

!-----trnth2o is the total transmittance between levels k1 and k2 due
!     to both line and continuum absorption.

!-----Compute th2o following Eq. (8.23).

   th2o(1) = th2o(1)*h2oexp(k,1)
   th2o(2) = th2o(2)*h2oexp(k,2)
   th2o(3) = th2o(3)*h2oexp(k,3)
   th2o(4) = th2o(4)*h2oexp(k,4)
   th2o(5) = th2o(5)*h2oexp(k,5)
   th2o(6) = th2o(6)*h2oexp(k,6)

   if (ne == 0) then

!-----Compute trnh2o following Eq. (8.25). fkw is given in Table 4.

      trnth2o =(fkw(1,ib)*th2o(1) &
              + fkw(2,ib)*th2o(2) &
              + fkw(3,ib)*th2o(3) &
              + fkw(4,ib)*th2o(4) &
              + fkw(5,ib)*th2o(5) &
              + fkw(6,ib)*th2o(6))


      tran    = tran*trnth2o

   elseif (ne == 1) then

!-----Compute trnh2o following Eqs. (8.25) and (4.27).

      tcon(1) = tcon(1)*conexp(k,1)

      trnth2o =(fkw(1,ib)*th2o(1) &
              + fkw(2,ib)*th2o(2) &
              + fkw(3,ib)*th2o(3) &
              + fkw(4,ib)*th2o(4) &
              + fkw(5,ib)*th2o(5) &
              + fkw(6,ib)*th2o(6))*tcon(1)

      tran    = tran*trnth2o

   else

!-----For band 3. This band is divided into 3 subbands.

      tcon(1)= tcon(1)*conexp(k,1)
      tcon(2)= tcon(2)*conexp(k,2)
      tcon(3)= tcon(3)*conexp(k,3)

!-----Compute trnh2o following Eqs. (4.29) and (8.25).

      trnth2o = (  gkw(1,1)*th2o(1) &
              + gkw(2,1)*th2o(2) &
              + gkw(3,1)*th2o(3) &
              + gkw(4,1)*th2o(4) &
              + gkw(5,1)*th2o(5) &
              + gkw(6,1)*th2o(6) ) * tcon(1) &
              + (  gkw(1,2)*th2o(1) &
              + gkw(2,2)*th2o(2) &
              + gkw(3,2)*th2o(3) &
              + gkw(4,2)*th2o(4) &
              + gkw(5,2)*th2o(5) &
              + gkw(6,2)*th2o(6) ) * tcon(2) &
              + (  gkw(1,3)*th2o(1) &
              + gkw(2,3)*th2o(2) &
              + gkw(3,3)*th2o(3) &
              + gkw(4,3)*th2o(4) &
              + gkw(5,3)*th2o(5) &
              + gkw(6,3)*th2o(6) ) * tcon(3)

      tran    = tran*trnth2o

   end if

   end subroutine h2okdis


!**********************************************************************
   subroutine n2okdis(ib,np,k,n2oexp,tn2o,tran)
!**********************************************************************
!   compute n2o transmittances between levels k1 and k2, using the k-distribution method with linear
!    pressure scaling.
!
!---- input parameters
!   spectral band (ib)
!   number of levels (np)
!   current level (k)
!   exponentials for n2o absorption (n2oexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to n2o absorption
!     for the various values of the absorption coefficient (tn2o)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none
   integer :: ib,np,k

!---- input parameters -----

   real(8) :: n2oexp(0:np,4)

!---- updated parameters -----

   real(8) :: tn2o(4),tran

!---- temporary arrays -----

   real(8) :: xc

!-----tn2o is computed from Eq. (8.23). 
!     xc is the total n2o transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.

!-----band 6

   if (ib == 6) then

      tn2o(1)=tn2o(1)*n2oexp(k,1)
      xc=   0.940414*tn2o(1)

      tn2o(2)=tn2o(2)*n2oexp(k,2)
      xc=xc+0.059586*tn2o(2)

!-----band 7

   else

      tn2o(1)=tn2o(1)*n2oexp(k,1)
      xc=   0.561961*tn2o(1)

      tn2o(2)=tn2o(2)*n2oexp(k,2)
      xc=xc+0.138707*tn2o(2)

      tn2o(3)=tn2o(3)*n2oexp(k,3)
      xc=xc+0.240670*tn2o(3)

      tn2o(4)=tn2o(4)*n2oexp(k,4)
      xc=xc+0.058662*tn2o(4)

   end if

   tran=tran*xc

   end subroutine n2okdis


!**********************************************************************
   subroutine ch4kdis(ib,np,k,ch4exp,tch4,tran)
!**********************************************************************
!   compute ch4 transmittances between levels k1 and k2, using the k-distribution method with
!    linear pressure scaling.
!
!---- input parameters
!   spectral band (ib)
!   number of levels (np)
!   current level (k)
!   exponentials for ch4 absorption (ch4exp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to ch4 absorption
!     for the various values of the absorption coefficient (tch4)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none
   integer :: ib,np,k

!---- input parameters -----

   real(8) :: ch4exp(0:np,4)

!---- updated parameters -----

   real(8) :: tch4(4),tran

!---- temporary arrays -----

   real(8) :: xc

!-----tch4 is computed from Eq. (8.23). 
!     xc is the total ch4 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.

!-----band 6

   if (ib == 6) then

      tch4(1)=tch4(1)*ch4exp(k,1)
      xc= tch4(1)

!-----band 7

   else

      tch4(1)=tch4(1)*ch4exp(k,1)
      xc=   0.610650*tch4(1)

      tch4(2)=tch4(2)*ch4exp(k,2)
      xc=xc+0.280212*tch4(2)

      tch4(3)=tch4(3)*ch4exp(k,3)
      xc=xc+0.107349*tch4(3)

      tch4(4)=tch4(4)*ch4exp(k,4)
      xc=xc+0.001789*tch4(4)

   end if

   tran=tran*xc

   end subroutine ch4kdis


!**********************************************************************
   subroutine comkdis(ib,np,k,comexp,tcom,tran)
!**********************************************************************
!  compute co2-minor transmittances between levels k1 and k2, using the k-distribution method
!   with linear pressure scaling.
!
!---- input parameters
!   spectral band (ib)
!   number of levels (np)
!   current level (k)
!   exponentials for co2-minor absorption (comexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to co2-minor absorption
!     for the various values of the absorption coefficient (tcom)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none
   integer :: ib,np,k

!---- input parameters -----

   real(8) :: comexp(0:np,6)

!---- updated parameters -----

   real(8) :: tcom(6),tran

!---- temporary arrays -----

   real(8) :: xc

!-----tcom is computed from Eq. (8.23). 
!     xc is the total co2 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 6.

!-----band 4

   if (ib == 4) then

      tcom(1)=tcom(1)*comexp(k,1)
      xc=   0.12159*tcom(1)
      tcom(2)=tcom(2)*comexp(k,2)
      xc=xc+0.24359*tcom(2)
      tcom(3)=tcom(3)*comexp(k,3)
      xc=xc+0.24981*tcom(3)
      tcom(4)=tcom(4)*comexp(k,4)
      xc=xc+0.26427*tcom(4)
      tcom(5)=tcom(5)*comexp(k,5)
      xc=xc+0.07807*tcom(5)
      tcom(6)=tcom(6)*comexp(k,6)
      xc=xc+0.04267*tcom(6)

!-----band 5

   else

      tcom(1)=tcom(1)*comexp(k,1)
      xc=   0.06869*tcom(1)
      tcom(2)=tcom(2)*comexp(k,2)
      xc=xc+0.14795*tcom(2)
      tcom(3)=tcom(3)*comexp(k,3)
      xc=xc+0.19512*tcom(3)
      tcom(4)=tcom(4)*comexp(k,4)
      xc=xc+0.33446*tcom(4)
      tcom(5)=tcom(5)*comexp(k,5)
      xc=xc+0.17199*tcom(5)
      tcom(6)=tcom(6)*comexp(k,6)
      xc=xc+0.08179*tcom(6)

   end if

   tran=tran*xc

   end subroutine comkdis


!**********************************************************************
   subroutine cfckdis(np,k,cfcexp,tcfc,tran)
!**********************************************************************
!  compute cfc-(11,12,22) transmittances between levels k1 and k2, using the k-distribution method with
!   linear pressure scaling.
!
!---- input parameters
!   number of levels (np)
!   current level (k)
!   exponentials for cfc absorption (cfcexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to cfc absorption
!     for the various values of the absorption coefficient (tcfc)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none

!---- input parameters -----

   integer :: k, np
   real(8) :: cfcexp(0:np)

!---- updated parameters -----

   real(8) :: tcfc,tran

!-----tcfc is the exp factors between levels k1 and k2. 

   tcfc=tcfc*cfcexp(k)
   tran=tran*tcfc

   end subroutine cfckdis


!**********************************************************************
   subroutine b10kdis(np,k,h2oexp,conexp,co2exp,n2oexp, &
      th2o,tcon,tco2,tn2o,tran)
!**********************************************************************
!
!   compute h2o (line and continuum),co2,n2o transmittances between
!   levels k1 and k2, using the k-distribution
!   method with linear pressure scaling.
!
!---- input parameters
!   number of levels (np)
!   current level (k)
!   exponentials for h2o line absorption (h2oexp)
!   exponentials for h2o continuum absorption (conexp)
!   exponentials for co2 absorption (co2exp)
!   exponentials for n2o absorption (n2oexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to h2o line absorption
!     for the various values of the absorption coefficient (th2o)
!   transmittance between levels k1 and k2 due to h2o continuum
!     absorption for the various values of the absorption
!     coefficient (tcon)
!   transmittance between levels k1 and k2 due to co2 absorption
!     for the various values of the absorption coefficient (tco2)
!   transmittance between levels k1 and k2 due to n2o absorption
!     for the various values of the absorption coefficient (tn2o)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none

   integer :: np,k

!---- input parameters -----

   real(8) :: h2oexp(0:np,5),conexp(0:np),co2exp(0:np,6), &
         n2oexp(0:np,2)

!---- updated parameters -----

   real(8) :: th2o(6),tcon(3),tco2(6),tn2o(4), &
         tran

!---- temporary arrays -----

   real(8) :: xx

!-----For h2o line. The k-distribution functions are given in Table 4.

   th2o(1)=th2o(1)*h2oexp(k,1)
   xx=   0.3153*th2o(1)

   th2o(2)=th2o(2)*h2oexp(k,2)
   xx=xx+0.4604*th2o(2)

   th2o(3)=th2o(3)*h2oexp(k,3)
   xx=xx+0.1326*th2o(3)

   th2o(4)=th2o(4)*h2oexp(k,4)
   xx=xx+0.0798*th2o(4)

   th2o(5)=th2o(5)*h2oexp(k,5)
   xx=xx+0.0119*th2o(5)

   tran=xx

!-----For h2o continuum. Note that conexp(k,3) is for subband 3a.

   tcon(1)=tcon(1)*conexp(k)
   tran=tran*tcon(1)

!-----For co2 (Table 6)

   tco2(1)=tco2(1)*co2exp(k,1)
   xx=    0.2673*tco2(1)

   tco2(2)=tco2(2)*co2exp(k,2)
   xx=xx+ 0.2201*tco2(2)

   tco2(3)=tco2(3)*co2exp(k,3)
   xx=xx+ 0.2106*tco2(3)

   tco2(4)=tco2(4)*co2exp(k,4)
   xx=xx+ 0.2409*tco2(4)

   tco2(5)=tco2(5)*co2exp(k,5)
   xx=xx+ 0.0196*tco2(5)

   tco2(6)=tco2(6)*co2exp(k,6)
   xx=xx+ 0.0415*tco2(6)

   tran=tran*xx

!-----For n2o (Table 5)

   tn2o(1)=tn2o(1)*n2oexp(k,1)
   xx=   0.970831*tn2o(1)

   tn2o(2)=tn2o(2)*n2oexp(k,2)
   xx=xx+0.029169*tn2o(2)

   tran=tran*(xx-1.0)

   end subroutine b10kdis


!mjs

!***********************************************************************
   subroutine cldovlp (np,k1,k2,ict,icb,icx,ncld,enn,ett, &
      cldhi,cldmd,cldlw)
!***********************************************************************
!
!     update the effective superlayer cloud fractions between  levels k1
!     and k2 following Eqs.(6.18)-(6.21). This assumes that the input
!     are the fractions between k1 and k2-1.
!
! input parameters
!
!  np:      number of layers
!  k1:      index for top level
!  k2:      index for bottom level
!  ict:     the level separating high and middle clouds
!  icb:     the level separating middle and low clouds
!  icx:     indeces sorted by enn in the three superlayers
!  enn:     effective fractional cloud cover of a layer
!  ett:     transmittance of a cloud layer
!
! update paremeters
!
!  cldhi:   Effective high-level cloudiness
!  cldmd:   Effective middle-level cloudiness
!  cldlw:   Effective low-level cloudiness
!
!
!***********************************************************************

   implicit none

   integer, intent(IN   ) :: np,k1,k2,ict,icb,icx(0:np),ncld(3)
   real(8),    intent(IN   ) :: enn(0:np), ett(0:np)

   real(8),    intent(INOUT) :: cldhi,cldmd,cldlw

   integer :: j,k,km,kx

   km=k2-1

   if    (km <  ict               ) then ! do high clouds

      kx=ncld(1)

      if(kx==1 .or. cldhi==0.) then
         cldhi = enn(km)
      else
      !if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
         cldhi = 0.0
         if(kx/=0) then
            do k=ict-kx,ict-1
               j=icx(k)
               if(j>=k1 .and. j<=km) cldhi=enn(j)+ett(j)*cldhi
            end do
         end if
      end if

   elseif(km >= ict .and. km < icb) then ! do middle clouds

      kx=ncld(2)

      if(kx==1 .or. cldmd==0.) then
         cldmd = enn(km)
      else
      !if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
         cldmd = 0.0
         if(kx/=0) then
            do k=icb-kx,icb-1
               j=icx(k)
               if(j>=k1 .and. j<=km) cldmd=enn(j)+ett(j)*cldmd
            end do
         end if
      end if

   else                                  ! do low clouds

      kx=ncld(3)

      if(kx==1 .or. cldlw==0.) then
         cldlw = enn(km)
      else
      !if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
         cldlw = 0.0
         if(kx/=0) then
            do k=np+1-kx,np
               j=icx(k)
               if(j>=k1 .and. j<=km) cldlw=enn(j)+ett(j)*cldlw
            end do
         end if
      end if

   end if

   end subroutine cldovlp

!***********************************************************************
   subroutine sfcflux (ibn,cb,dcb,ns,fs,tg,eg,tv,ev,rv,bs,dbs,rflxs)
!***********************************************************************
! Compute emission and reflection by an homogeneous/inhomogeneous 
!  surface with vegetation cover.
!
!-----Input parameters
!  index for the spectral band (ibn)
!  number of sub-grid box (ns)
!  fractional cover of sub-grid box (fs)
!  sub-grid ground temperature (tg)
!  sub-grid ground emissivity (eg)
!  sub-grid vegetation temperature (tv)
!  sub-grid vegetation emissivity (ev)
!  sub-grid vegetation reflectivity (rv)
!                      
!-----Output parameters                                Unit
!  Emission by the surface (ground+vegetation) (bs)    W/m^2
!  Derivative of bs rwt to temperature (dbs)           W/m^2
!  Reflection by the surface (rflxs)                   W/m^2
!
!**********************************************************************
   implicit none

!---- input parameters -----
   integer :: ibn,ns
   real(8) :: cb(6,10)
   real(8) :: dcb(5,10)
   real(8) :: fs(ns),tg(ns),eg(ns,10)
   real(8) :: tv(ns),ev(ns,10),rv(ns,10)

!---- output parameters -----
   real(8) :: bs,dbs,rflxs

!---- temporary arrays -----
   integer :: j
   real(8) :: bg(ns),bv(ns),dbg(ns),dbv(ns)
   real(8) :: xx,yy,zz

!***************************************************************

!-----compute planck flux and the derivative wrt to temperature
!     c.f. eqs. (3.11) and (3.12).

   do j=1,ns
      call planck(ibn,cb,tg(j),bg(j))
      call planck(ibn,cb,tv(j),bv(j))
      call plancd(ibn,dcb,tg(j),dbg(j))
      call plancd(ibn,dcb,tv(j),dbv(j))
   end do

!-----

   if (fs(1) > 0.9999) then
      if (ev(1,ibn) < 0.0001 .and. rv(1,ibn) < 0.0001) then

!-----for homogeneous surface without vegetation
!     following Eqs. (9.4), (9.5), and (3.13)

         bs =eg(1,ibn)*bg(1)
         dbs=eg(1,ibn)*dbg(1)
         rflxs=1.0-eg(1,ibn)
      else

!-----With vegetation
!     following Eqs. (9.1), (9.3), and (9.13)

         xx=ev(1,ibn)*bv(1)
         yy=1.0-ev(1,ibn)-rv(1,ibn)
         zz=1.0-eg(1,ibn)
         bs=yy*(eg(1,ibn)*bg(1)+zz*xx)+xx

         xx=ev(1,ibn)*dbv(1)
         dbs=yy*(eg(1,ibn)*dbg(1)+zz*xx)+xx

         rflxs=rv(1,ibn)+zz*yy*yy/(1.0-rv(1,ibn)*zz)
      end if
   else

!-----for nonhomogeneous surface

      bs=0.0
      dbs=0.0
      rflxs=0.0

      if (ev(1,ibn) < 0.0001 .and. rv(1,ibn) < 0.0001) then

!-----No vegetation, following Eqs. (9.9), (9.10), and (9.13)

         do j=1,ns
            bs=bs+fs(j)*eg(j,ibn)*bg(j)
            dbs=dbs+fs(j)*eg(j,ibn)*dbg(j)
            rflxs=rflxs+fs(j)*(1.0-eg(j,ibn))
         end do
      else

!-----With vegetation, following Eqs. (9.6), (9.7), and (9.13)

         do j=1,ns
            xx=ev(j,ibn)*bv(j)
            yy=1.0-ev(j,ibn)-rv(j,ibn)
            zz=1.0-eg(j,ibn)
            bs=bs+fs(j)*(yy*(eg(j,ibn)*bg(j)+zz*xx)+xx)

            xx=ev(j,ibn)*dbv(j)
            dbs=dbs+fs(j)*(yy*(eg(j,ibn)*dbg(j)+zz*xx)+xx)

            rflxs=rflxs+fs(j)*(rv(j,ibn)+zz*yy*yy &
                  /(1.0-rv(j,ibn)*zz))
         end do
      end if
   end if

   end subroutine sfcflux

!mjs


!***********************************************************************
   subroutine mkicx(np,ict,icb,enn,icx,ncld)
!***********************************************************************

   implicit none

   integer, intent(IN   ) :: np, ict, icb
   real(8),    intent(IN   ) :: enn(0:np)
   integer, intent(INOUT) :: icx(0:np)
   integer, intent(  OUT) :: ncld(3)

   call sortit(enn(  0:ict-1),icx(  0:ict-1),ncld(1),  0, ict-1)
   call sortit(enn(ict:icb-1),icx(ict:icb-1),ncld(2),ict, icb-1)
   call sortit(enn(icb:np   ),icx(icb:   np),ncld(3),icb,    np)

   end subroutine mkicx

!***********************************************************************
   subroutine SORTIT(ENN,LST,NCL,ibg,iend)
!***********************************************************************

   implicit none
   
   integer, intent(  OUT) :: NCL
   integer, intent(IN   ) :: IBG, IEND
   real(8),    intent(IN   ) :: ENN(IBG:IEND)
   integer, intent(INOUT) :: LST(IBG:IEND)

   integer :: L, LL, I
   real(8) :: ENO

   if(ENN(IBG)>0) then
      NCL=1
   else
      NCL=0
   end if

   do L=IBG+1,IEND
      ENO = ENN(LST(L))

      if(ENO>0) NCL=NCL+1

      LL=LST(L)
      I=L-1
      do while(I>IBG-1)
         if (ENN(LST(I))<=ENO) exit
         LST(I+1)=LST(I)
         I=I-1
      end do
      LST(I+1)=LL
   end do

   end subroutine SORTIT


   subroutine getirtau1(ib,nlevs,dp,fcld,reff,hydromets,&
                       tcldlyr,enn, &
                       aib_ir1, awb_ir1, &
                       aiw_ir1, aww_ir1, &
                       aig_ir1, awg_ir1, CONS_GRAV)

! !USES:

      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: ib                 !  Band number
      integer, intent(IN ) :: nlevs              !  Number of levels
      real(8),    intent(IN ) :: dp(nlevs)          !  Delta pressure in Pa
      real(8),    intent(IN ) :: fcld(nlevs)        !  Cloud fraction (used sometimes)
      real(8),    intent(IN ) :: reff(nlevs,4)      !  Effective radius (microns)
      real(8),    intent(IN ) :: hydromets(nlevs,4) !  Hydrometeors (kg/kg)
      real(8),    intent(IN ) :: aib_ir1(3,10), awb_ir1(4,10), aiw_ir1(4,10)
      real(8),    intent(IN ) :: aww_ir1(4,10), aig_ir1(4,10), awg_ir1(4,10)
      real(8),    intent(IN ) :: CONS_GRAV

! !OUTPUT PARAMETERS:
      real(8),    intent(OUT) :: tcldlyr(0:nlevs  ) !  Flux transmissivity?
      real(8),    intent(OUT) ::     enn(0:nlevs  ) !  Flux transmissivity of a cloud layer?

! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for infrared wavelengths
!  Slots for reff, hydrometeors and tauall are as follows:
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
!    2011.11.18   MAT moved to Radiation_Shared and revised arg list, units
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer            :: k
      real(8)               :: taucld1,taucld2,taucld3,taucld4
      real(8)               :: g1,g2,g3,g4,gg
      real(8)               :: w1,w2,w3,w4,ww
      real(8)               :: ff,tauc

      real(8)               :: reff_snow


!-----Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!     Rain optical thickness is set to 0.00307 /(gm/m**2).
!     It is for a specific drop size distribution provided by Q. Fu.

      tcldlyr(0) = 1.0
      enn    (0) = 0.0

      do k = 1, nlevs
         if(reff(k,1)<=0.0) then
            taucld1=0.0
         else
            taucld1=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,1))*(aib_ir1(1,ib)+aib_ir1(2,ib)/&
                  reff(k,1)**aib_ir1(3,ib))
         end if

            taucld2=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,2))*(awb_ir1(1,ib)+(awb_ir1(2,ib)+&
                  (awb_ir1(3,ib)+awb_ir1(4,ib)*reff(k,2))*reff(k,2))*reff(k,2))

            taucld3=0.00307*(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,3))

!-----Below, we use the table of coefficients tabulated for suspended
!     cloud ice particles (aib_ir1) for falling snow. These coefficients
!     lead to unphysical (negative) values of cloud optical thickness 
!     for effective radii greater than 113 microns. By restricting the 
!     effective radius of snow to 112 microns, we prevent unphysical 
!     optical thicknesses.

         reff_snow = min(reff(k,4),112.0)

         if(reff_snow<=0.0) then
            taucld4=0.0
         else
            taucld4=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,4))*(aib_ir1(1,ib)+aib_ir1(2,ib)/&
                  reff_snow**aib_ir1(3,ib))
         end if

!-----Compute cloud single-scattering albedo and asymmetry factor for
!     a mixture of ice particles and liquid drops following 
!     Eqs. (6.5), (6.6), (6.15) and (6.16).
!     Single-scattering albedo and asymmetry factor of rain are set
!     to 0.54 and 0.95, respectively, based on the information provided
!     by Prof. Qiang Fu.

         tauc=taucld1+taucld2+taucld3+taucld4

         if (tauc > 0.02 .and. fcld(k) > 0.01) then

            w1=taucld1*(aiw_ir1(1,ib)+(aiw_ir1(2,ib)+(aiw_ir1(3,ib) &
                  +aiw_ir1(4,ib)*reff(k,1))*reff(k,1))*reff(k,1))
            w2=taucld2*(aww_ir1(1,ib)+(aww_ir1(2,ib)+(aww_ir1(3,ib)&
                  +aww_ir1(4,ib)*reff(k,2))*reff(k,2))*reff(k,2))
            w3=taucld3*0.54
            w4=taucld4*(aiw_ir1(1,ib)+(aiw_ir1(2,ib)+(aiw_ir1(3,ib) &
                  +aiw_ir1(4,ib)*reff_snow)*reff_snow)*reff_snow)
            ww=(w1+w2+w3+w4)/tauc


            g1=w1*(aig_ir1(1,ib)+(aig_ir1(2,ib)+(aig_ir1(3,ib)&
                  +aig_ir1(4,ib)*reff(k,1))*reff(k,1))*reff(k,1))
            g2=w2*(awg_ir1(1,ib)+(awg_ir1(2,ib)+(awg_ir1(3,ib) &
                  +awg_ir1(4,ib)*reff(k,2))*reff(k,2))*reff(k,2))
            g3=w3*0.95
            g4=w4*(aig_ir1(1,ib)+(aig_ir1(2,ib)+(aig_ir1(3,ib)&
                  +aig_ir1(4,ib)*reff_snow)*reff_snow)*reff_snow)
       
           if (abs(w1+w2+w3+w4).gt.0.0) then 
            gg=(g1+g2+g3+g4)/(w1+w2+w3+w4)
	    else 
	    gg=0.5
	   end if 

!-----Parameterization of LW scattering following Eqs. (6.11)
!     and (6.12). 

            ff=0.5+(0.3739+(0.0076+0.1185*gg)*gg)*gg

!ALT: temporary protection against negative cloud optical thickness

            tauc=max((1.-ww*ff),0.0)*tauc

!-----compute cloud diffuse transmittance. It is approximated by using 
!     a diffusivity factor of 1.66.

            tcldlyr(k) = exp(-1.66*tauc)
            enn    (k) = fcld(k)*(1.0-tcldlyr(k)) ! N in the documentation (6.13)

         else

            tcldlyr(k) = 1.0
            enn    (k) = 0.0

         end if
      end do


   end subroutine getirtau1

end module IRRADMOD
