subroutine nst_init_
      use geos_nstmod, only: geos_nst_init
      implicit none
      call geos_nst_init()
end subroutine nst_init_
!*******************************************************************************************

subroutine nst_set_ (mype,mype_io)
      use kinds, only: i_kind
      use geos_nstmod, only: geos_nst_set
      implicit none
      integer(i_kind), intent(in   ) :: mype,mype_io
      call geos_nst_set()
end subroutine nst_set_
!*******************************************************************************************

subroutine nst_final_
      use geos_nstmod, only: geos_nst_final
      implicit none
      call geos_nst_final()
end subroutine nst_final_
!*******************************************************************************************

subroutine deter_nst_(dlat_earth,dlon_earth,obstime,zob,tref,dtw,dtc,tz_tr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    deter_nst                     determine NSST variable at observation location over water
!   prgmmr: Xu Li           org: np2                date: 2011-04-08
! abstract:  determines NSST variables over water surface type based on surrounding surface types
!
! program history log:
!   2011-04-08 Li
!   2011-10-20 RT/ Akella                GEOS customizations
!   2013-01-23 parrish - change from grdcrd to grdcrd1 (to allow successful
!   2014-05-30 Akella/RT - pass mu-skin via RC
!
!   input argument list:
!     obstime                             - observation time relative to analysis time
!     dlat_earth                          - earth latitude in radians
!     dlon_earth                          - earth longitude in radians
!     zob                                 - obs. depth in the water
!
!   output argument list:
!     tref                                - oceanic foundation temperature
!      dtw                                - diurnal warming at depth zob
!      dtc                                - sublayer cooling at depth zob
!    tz_tr                                - d(Tz)/d(Tbar); is DIFFERENT to Xu Li. Xu Li: d(Tz)/d(Tr); Tr: foundation or ref. Temp.
!                                         - in GEOS context it is = d(Tz)/d(Ts)
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
     use kinds,       only: r_kind,i_kind
     use constants,   only: zero, one
     use mpimod,      only: mype
     use gridmod,     only: nlat,nlon,regional,tll2xy,nlat_sfc,nlon_sfc,rlats_sfc,rlons_sfc
     use guess_grids, only: nfldsfc,hrdifsfc
     use mpimod,      only: mype
     use gsi_nstcouplermod, only: tref_full,dt_cool_full,dt_warm_full,z_c_full,z_w_full
     use satthin,     only: isli_full
   
     use GSI_GridCompMod, only: MU_SKIN=>GEOS_MU_SKIN

!    use ESMF, only: ESMF_ConfigGetAttribute
!    use m_die, only: die
     implicit none

     real(r_kind), intent(in ) :: dlat_earth,dlon_earth,obstime,zob
     real(r_kind), intent(out) :: tref,dtw,dtc,tz_tr

!    local variables
     character(len=*), parameter:: myname_='deter_nst_'
     real(r_kind):: dt_cool,z_c,z_w,dt_warm
     integer(i_kind):: istyp00,istyp01,istyp10,istyp11
     integer(i_kind):: itnst,itnstp
     integer(i_kind):: ix,iy,ixp,iyp,j,i,k
     real(r_kind):: dx,dy,dx1,dy1,w00,w10,w01,w11,dtnst,dtnstp,wgtmin
     real(r_kind):: tref_00,tref_01,tref_10,tref_11,tr_tmp
     real(r_kind):: dt_cool_00,dt_cool_01,dt_cool_10,dt_cool_11
     real(r_kind):: z_c_00,z_c_01,z_c_10,z_c_11,z_c_tmp
     real(r_kind):: dt_warm_00,dt_warm_01,dt_warm_10,dt_warm_11
     real(r_kind):: z_w_00,z_w_01,z_w_10,z_w_11,z_w_tmp

     real(r_kind):: wgtavg,dlat,dlon
     logical outside


!  Read in the temperature profile exponent: mu_skin used for skin SST analysis
!  ----------------------------------------------------------------------------
!  CALL ESMF_ConfigGetAttribute(CF, MU_SKIN, label = 'mu_skin:', default=0.2_r_kind, rc=STATUS)
!  if (status/=0) then
!      call die ( myname_,': failed to get mu_skin' )
!  endif
!  if(IamRoot) print *,trim(Iam),': Set MU_SKIN= ', MU_SKIN

!    Initialize output.
     tref  = zero
     dtw   = zero
     dtc   = zero
     tz_tr = one

     if(regional)then
        call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
     else
        dlat=dlat_earth
        dlon=dlon_earth
        call grdcrd1(dlat,rlats_sfc,nlat_sfc,1)
        call grdcrd1(dlon,rlons_sfc,nlon_sfc,1)
     end if

     iy=int(dlon); ix=int(dlat)
     dy  =dlon-iy; dx  =dlat-ix
     dx1 =one-dx;    dy1 =one-dy
     w00=dx1*dy1; w10=dx*dy1; w01=dx1*dy; w11=one-w00-w10-w01

     ix=min(max(1,ix),nlat_sfc); iy=min(max(0,iy),nlon_sfc)
     ixp=min(nlat_sfc,ix+1); iyp=iy+1
     if(iy==0) iy=nlon_sfc
     if(iyp==nlon_sfc+1) iyp=1

!    Get time interpolation factors for nst files
     if(obstime > hrdifsfc(1) .and. obstime <= hrdifsfc(nfldsfc))then
        do j=1,nfldsfc-1
           if(obstime > hrdifsfc(j) .and. obstime <= hrdifsfc(j+1))then
              itnst=j
              itnstp=j+1
              dtnst=(hrdifsfc(j+1)-obstime)/(hrdifsfc(j+1)-hrdifsfc(j))
           end if
        end do
     else if(obstime <=hrdifsfc(1))then
        itnst=1
        itnstp=1
        dtnst=one
     else
        itnst=nfldsfc
        itnstp=nfldsfc
        dtnst=one
     end if
     dtnstp=one-dtnst

!    Set surface type flag.

     istyp00 = isli_full(ix ,iy )
     istyp10 = isli_full(ixp,iy )
     istyp01 = isli_full(ix ,iyp)
     istyp11 = isli_full(ixp,iyp)
!
!    Use the time interpolation factors for nst files
!
     tref_00    = tref_full   (ix ,iy ,itnst)*dtnst + tref_full   (ix ,iy ,itnstp)*dtnstp
     tref_01    = tref_full   (ix ,iyp,itnst)*dtnst + tref_full   (ix ,iyp,itnstp)*dtnstp
     tref_10    = tref_full   (ixp,iy ,itnst)*dtnst + tref_full   (ixp,iy ,itnstp)*dtnstp
     tref_11    = tref_full   (ixp,iyp,itnst)*dtnst + tref_full   (ixp,iyp,itnstp)*dtnstp

     dt_cool_00 = dt_cool_full(ix ,iy ,itnst)*dtnst + dt_cool_full(ix ,iy ,itnstp)*dtnstp
     dt_cool_01 = dt_cool_full(ix ,iyp,itnst)*dtnst + dt_cool_full(ix ,iyp,itnstp)*dtnstp
     dt_cool_10 = dt_cool_full(ixp,iy ,itnst)*dtnst + dt_cool_full(ixp,iy ,itnstp)*dtnstp
     dt_cool_11 = dt_cool_full(ixp,iyp,itnst)*dtnst + dt_cool_full(ixp,iyp,itnstp)*dtnstp

     z_c_00     = z_c_full    (ix ,iy ,itnst)*dtnst + z_c_full    (ix ,iy ,itnstp)*dtnstp
     z_c_01     = z_c_full    (ix ,iyp,itnst)*dtnst + z_c_full    (ix ,iyp,itnstp)*dtnstp
     z_c_10     = z_c_full    (ixp,iy ,itnst)*dtnst + z_c_full    (ixp,iy ,itnstp)*dtnstp
     z_c_11     = z_c_full    (ixp,iyp,itnst)*dtnst + z_c_full    (ixp,iyp,itnstp)*dtnstp

     dt_warm_00 = dt_warm_full(ix ,iy ,itnst)*dtnst + dt_warm_full(ix ,iy ,itnstp)*dtnstp
     dt_warm_01 = dt_warm_full(ix ,iyp,itnst)*dtnst + dt_warm_full(ix ,iyp,itnstp)*dtnstp
     dt_warm_10 = dt_warm_full(ixp,iy ,itnst)*dtnst + dt_warm_full(ixp,iy ,itnstp)*dtnstp
     dt_warm_11 = dt_warm_full(ixp,iyp,itnst)*dtnst + dt_warm_full(ixp,iyp,itnstp)*dtnstp

     z_w_00     = z_w_full    (ix ,iy ,itnst)*dtnst + z_w_full    (ix ,iy ,itnstp)*dtnstp
     z_w_01     = z_w_full    (ix ,iyp,itnst)*dtnst + z_w_full    (ix ,iyp,itnstp)*dtnstp
     z_w_10     = z_w_full    (ixp,iy ,itnst)*dtnst + z_w_full    (ixp,iy ,itnstp)*dtnstp
     z_w_11     = z_w_full    (ixp,iyp,itnst)*dtnst + z_w_full    (ixp,iyp,itnstp)*dtnstp

!    Interpolate nst variables to obs location (water surface only)

     wgtavg  = zero
     tr_tmp  = zero
     dt_cool = zero
     z_c_tmp = zero
     dt_warm = zero
     z_w_tmp = zero

     if (istyp00 == 0)then
        wgtavg  = wgtavg  + w00
        tr_tmp  = tr_tmp  + w00*tref_00
        dt_cool = dt_cool + w00*dt_cool_00
        z_c_tmp = z_c_tmp + w00*z_c_00
        dt_warm = dt_warm + w00*dt_warm_00
        z_w_tmp = z_w_tmp + w00*z_w_00
     end if
     if(istyp01 == 0)then
        wgtavg  = wgtavg  + w01
        tr_tmp  = tr_tmp  + w01*tref_01
        dt_cool = dt_cool + w01*dt_cool_01
        z_c_tmp = z_c_tmp + w01*z_c_01
        dt_warm = dt_warm + w01*dt_warm_01
        z_w_tmp = z_w_tmp + w01*z_w_01
     end if
     if(istyp10 == 0)then
        wgtavg  = wgtavg  + w10
        tr_tmp  = tr_tmp  + w10*tref_10
        dt_cool = dt_cool + w10*dt_cool_10
        z_c_tmp = z_c_tmp + w10*z_c_10
        dt_warm = dt_warm + w10*dt_warm_10
        z_w_tmp = z_w_tmp + w10*z_w_10
     end if
     if(istyp11 == 0)then
        wgtavg  = wgtavg  + w11
        tr_tmp  = tr_tmp  + w11*tref_11
        dt_cool = dt_cool + w11*dt_cool_11
        z_c_tmp = z_c_tmp + w11*z_c_11
        dt_warm = dt_warm + w11*dt_warm_11
        z_w_tmp = z_w_tmp + w11*z_w_11
     end if

     if(wgtavg > zero)then
        tr_tmp  = tr_tmp/wgtavg
        tref    = tr_tmp

        z_w_tmp = z_w_tmp/wgtavg
        z_w     = z_w_tmp
        dt_warm = dt_warm/wgtavg

        dt_cool = dt_cool/wgtavg
        z_c_tmp = z_c_tmp/wgtavg
        z_c     = z_c_tmp

! Jacobian calculation: d(T(z))/d(Ts)

        tz_tr = one

        if ( (zob > z_w) .AND. (zob > z_c) ) then
              ! should not have obs that is deeper than both z_w & z_c
              ! once you fix sfcpt(0)- make it agree with frac water & ice in BKG, 
              ! this branch of the if loop should never be exercised.
              ! For now, set diurnal fields to zero.  SA. 02/06/2012
              dtc   = zero
              dtw   = zero
              tz_tr = zero          ! gradient should NOT impact temperature analysis.
!             if(mype==0) WRITE(885,771)  tref, ix, iy, ixp, iyp
!             if(mype==0) WRITE(885,771)  zob, z_c, z_w, dt_cool, dt_warm
        else

           if (zob .le. z_c) then
              dtc   = (one-zob/z_c) * dt_cool                                  ! linear T(z) profile in cool-layer
              dtw   = dt_warm                                                  ! account for complete warm-layer temp. rise
              tz_tr = one
!             if(mype==0) WRITE(886,771)  zob, z_c, z_w, dt_cool, dt_warm

           elseif ( (zob > z_c) .and. (zob <= 0.05)) then                      ! z_c < zob < 5cm. That is, zob certainly corresponds to satellite measurement (IR, MW)
              dtc   = zero                                                     ! sensor does not feel cool-layer effects.
              dtw   = dt_warm * (one- ((zob-z_c)/(z_w-z_c))**MU_SKIN)          ! ZENG & BELJAARS warm layer T(z) profile 
              tz_tr = one                                                      ! THIS IS AN Approx. for Sattelite (MW) data. Correct way is below. --this MIGHT NEED REVISIT! 10/03/2014.

           else                                                                ! ((zob > z_c) .AND. (zob .le. z_w)) 
              dtc   = zero                                                     ! sensor does not feel cool-layer effects.
              dtw   = dt_warm * (one- ((zob-z_c)/(z_w-z_c))**MU_SKIN)          ! ZENG & BELJAARS warm layer T(z) profile 
              tz_tr = one- ((zob-z_c)/(z_w-z_c))**MU_SKIN                      ! larger zob => smaller tz_tr for MU_SKIN ~0.2- 0.3. But if MU_SKIN ~1, tz_tr ~1.
                                                                               ! Implies deeper (in water) obs do not change Ts as much. Which makes sense from an ATMOSPHERIC point of view.
!             if(mype==0) WRITE(887,771)  zob, z_c, z_w, dt_cool, dt_warm
           end if
!          call cal_tztr_(dt_warm, dt_cool, z_c, z_w, zob, tz_tr)
        end if

     end if

! z_c >=0 by definition. in GEOS. 

! keep Xu Li's code for future ref.
!       dtw = fac_dtl*dt_warm*(one-min(zob,z_w)/z_w)
!       if ( z_c > zero ) then
!         dtc = fac_tsl*dt_cool*(one-min(zob,z_c)/z_c)
!       else
!         dtc = zero
!       endif
!       call cal_tztr_(dt_warm, dt_cool, z_c, z_w, zob, tz_tr)
!--

!771 FORMAT(F10.4, 2x,4I4)
!771 FORMAT(E12.4, 2x, 4E12.4)
end subroutine deter_nst_
!*******************************************************************************************

!subroutine cal_tztr_(dt_warm, dt_cool, z_c, z_w, z, tztr)
!
! abstract: calculate d(Tz)/d(Ts) with T-Profile info from NSST Model
!
!   prgmmr: li, xu           org: np23                date: 2011-04-08
! input variables
!
! dt_warm :       diurnal warming amount at the surface
! xz      :       DTL depth                           (M)
! c_0     :       coefficint 1 to calculate d(Tc)/d(Ts)
! c_d     :       coefficint 2 to calculate d(Tc)/d(Ts)
! w_0     :       coefficint 1 to calculate d(Tw)/d(Ts)
! w_d     :       coefficint 2 to calculate d(Tw)/d(Ts)
!
! output variables
!
! tztr     :      d(Tz)/d(Tr)

! use kinds, only: r_kind
! use constants, only: one, zero
! real(kind=r_kind), intent(in)  :: dt_warm, dt_cool, z_c, z_w, z 
! real(kind=r_kind), intent(out) :: tztr

! set to some arbitrary value. figure it out later.
! tztr = one
! tztr = zero

! keep Xu Li's code for future ref.
! c1 = one-fac_dtl*w_0+fac_tsl*c_0
! c2 = one+fac_tsl*c_0
!
! tztr = one
!
! if ( dt_warm > zero .and.  c1 /= zero ) then
!   if ( z <= zc  ) then
!       tztr = (one+z*(fac_dtl*w_d-fac_tsl*c_d))/c1
!   elseif ( z > zc .and. z < zw ) then
!       tztr = (one+fac_tsl*c_0+z*fac_dtl*w_d)/c1
!   endif
!elseif ( dt_warm == zero .and. c2 /= zero ) then
!   if ( z <= zc ) then
!       tztr = (one-z*fac_tsl*c_d)/c2
!   endif
!endif
!
!if ( tztr <= one .and. tztr > half ) then
!        tztr = tztr
!else
!!   write(*,'(a,2I2,2F12.6,F9.3,5F12.6,F8.3,F9.6,F8.3)') ' cal_tztr : ',fac_dtl,fac_tsl,c1,c2,dt_warm,c_0,c_d,w_0,w_d,zc,zw,z,tztr
!endif
!--

!end subroutine cal_tztr_
!*******************************************************************************************
subroutine skindepth_(obstype, sd_rad)
!
! abstract: Get skin depth (instrument dependent). Ideally, a skin-depth model calculates the channel dependent sd
!
! input:                obstype
! output:               sd_rad-		penetration depth in meters
!
! to-do:                add channel, freq to input. Based on that table look-up to 
!                       find absroption coeff and penetration depth instead of assigning as done below.
!
! program history log:
!   2011-04-08  li
!   2012-02/16  akella. depth based on instrument type. see
!                       http://disc.sci.gsfc.nasa.gov/oceans/additional/science-focus/modis/MODIS_and_AIRS_SST_comp.shtml
!                       Table 1- Donlon et al, 2007: http://dx.doi.org/10.1175/BAMS-88-8-1197
!

  use kinds, only: r_kind
  implicit none
  character(10),     intent(in)  :: obstype
  real(kind=r_kind), intent(out) :: sd_rad

! a default value- 15 microns
  sd_rad = 0.000015_r_kind

! MW radiance
  if ( obstype=='amsua'     .OR. obstype=='amsub'     .OR. obstype=='mhs'         .OR. obstype=='msu'  .OR.&
       obstype=='ssmi'      .OR. obstype=='ssmis'     .OR. obstype=='hsb'         .OR. obstype=='atms' .OR.&
       obstype=='amsre_low' .OR. obstype=='amsre_mid' .OR. obstype=='amsre_high'                       .OR.&
       obstype=='tmi'       .OR. obstype=='gmi_low'   .OR. obstype=='gmi_hig' ) then

!   1.25 mm. range is ~1- 1.5 mm.
    sd_rad = 0.00125_r_kind 

! IR radiance
  elseif ( obstype=='airs'  .OR. obstype=='hirs3'      .OR. obstype=='hirs4'      .OR.                     &
           obstype=='sndr'  .OR. obstype=='goes_img'   .OR.                                                &
           obstype=='sndr1' .OR. obstype=='sndr2'      .OR. obstype=='sndr3'      .OR. obstype=='sndr4'.OR.&
           obstype=='avhrr' .OR. obstype=='avhrr_navy' .OR.                                                &
           obstype=='iasi'  .OR. obstype=='seviri' ) then
           
!   15 microns
    sd_rad = 0.000015_r_kind
  endif
  
end subroutine skindepth_
!*******************************************************************************************

