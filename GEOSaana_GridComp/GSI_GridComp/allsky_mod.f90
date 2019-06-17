 module allsky_mod
!$$$ module documentation block
!           .      .    .                                       .
! module:   allsky_mod
!   prgmmr: mkim      org: np2                 date: 2016-03-02
!
! abstract: This module contains routines related allsky microwave radiance da 
!
! subroutines included:
!   sub obserr_allsky_mw  - calculates observation error for cloud affected radiance
!   sub deter_cld_rbc_idx - determine cld_rbc_idx for bias correction for allsky mw radiance
!
! Variable Definitions:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block
  use kinds, only: r_kind,i_kind
  implicit none

! set default to private
  private
! set routines used externally to public
  public :: obserr_allsky_mw, deter_cld_rbc_idx 

  real(r_kind),parameter,dimension(15) :: cclr_amsua=(/0.05_r_kind, 0.03_r_kind, 0.03_r_kind, 0.02_r_kind, 0.00_r_kind, &
                                                           0.10_r_kind, 0.00_r_kind, 0.00_r_kind, 0.00_r_kind, 0.00_r_kind, &
                                                           0.00_r_kind, 0.00_r_kind, 0.00_r_kind, 0.00_r_kind, 0.03_r_kind/)
  real(r_kind),parameter,dimension(15) :: ccld_amsua=(/0.60_r_kind, 0.45_r_kind, 0.40_r_kind, 0.45_r_kind, 1.00_r_kind, &
                                                           1.50_r_kind, 0.00_r_kind, 0.00_r_kind, 0.00_r_kind, 0.00_r_kind, &
                                                           0.00_r_kind, 0.00_r_kind, 0.00_r_kind, 0.00_r_kind, 0.20_r_kind/)

  real(r_kind),parameter,dimension(13) :: Cclr_gmi=(/0.05_r_kind, 0.05_r_kind, 0.05_r_kind, 0.05_r_kind, &
                                                          0.05_r_kind, 0.05_r_kind, 0.05_r_kind, 0.05_r_kind, &
                                                          0.05_r_kind, 0.05_r_kind, 0.05_r_kind, 0.05_r_kind, &
                                                          0.05_r_kind/)
  real(r_kind),parameter,dimension(13) :: Ccld1_gmi=(/0.2_r_kind, 0.2_r_kind, 0.2_r_kind, 0.2_r_kind, &
                                                          0.2_r_kind, 0.2_r_kind, 0.2_r_kind, 0.2_r_kind, &
                                                          0.2_r_kind, 0.3_r_kind, 0.2_r_kind, 0.3_r_kind, &
                                                          0.3_r_kind/)
  real(r_kind),parameter,dimension(13) :: Ccld2_gmi=(/0.5_r_kind, 0.5_r_kind, 0.5_r_kind, 0.5_r_kind, &
                                                          0.5_r_kind, 0.5_r_kind, 0.5_r_kind, 0.50_r_kind, &
                                                          0.5_r_kind, 0.5_r_kind, 0.5_r_kind, 0.50_r_kind, &
                                                          0.5_r_kind/)

! Large error
  real(r_kind),parameter,dimension(13) :: tnoise1_gmi=(/17.0_r_kind,  23.0_r_kind, 13.0_r_kind, 25.0_r_kind, &
                                                         11.0_r_kind, 13.0_r_kind, 23.0_r_kind, 10.0_r_kind, &
                                                         20.0_r_kind, 15.0_r_kind, 20.0_r_kind,  8.0_r_kind, &
                                                         13.0_r_kind/)
! Small error
!  real(r_kind),parameter,dimension(13) :: tnoise1_gmi=(   / 3.5_r_kind,  5.0_r_kind,  7.0_r_kind, 15.0_r_kind, &
!                                                            8.0_r_kind, 11.0_r_kind,  4.0_r_kind, 10.0_r_kind, &
!                                                           12.0_r_kind,  7.0_r_kind,  7.0_r_kind,  4.0_r_kind, &
!                                                            5.0_r_kind/)
! Medium error
!  real(r_kind),parameter,dimension(13) :: tnoise1_gmi=(   / 3.5_r_kind,  5.0_r_kind, 10.0_r_kind, 17.0_r_kind, &
!                                                           10.0_r_kind, 12.0_r_kind,  4.0_r_kind, 10.0_r_kind, &
!                                                           12.0_r_kind, 12.0_r_kind, 12.0_r_kind,  6.0_r_kind, &
!                                                            8.0_r_kind/)

contains

subroutine deter_cld_rbc_idx(clwp_obs,clwp_guess,obstype,nchanl,cld_rbc_idx)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:   obserr_allsky_mw 
!   prgmmr: mkim          org: np23                date: 2016-03-02
!
! abstract: define observation error for all-sky microwave radiance data 
!
! program history log:
!   2014-01-31  mkim
!   2014-08-31  mkim : Adding obs error model for GMI and TMI
!
!   input argument list:
!     error0    - input observation error before applying changes for cloudy sky condition
!
!   output argument list:
!     error0      - observation error 
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: zero,one
  implicit none

  integer(i_kind) :: i
  integer(i_kind),                 intent(in   ) :: nchanl 
  real(r_kind), dimension(nchanl), intent(inout) :: cld_rbc_idx
  real(r_kind),                    intent(in   ) :: clwp_obs,clwp_guess
  character(10),                   intent(in   ) :: obstype 

  if(obstype == 'amsua')  then
     do i=1,nchanl
       if ((clwp_obs-cclr_amsua(i))*(clwp_guess-cclr_amsua(i))<zero .or.  &
        abs(clwp_obs-clwp_guess)>=0.005_r_kind) cld_rbc_idx(i)=zero
     end do
  else if(obstype == 'gmi') then
     do i=1,nchanl
       if (clwp_obs .le. cclr_gmi(i) .and. clwp_guess .le. cclr_gmi(i) .and. abs(clwp_obs-clwp_guess) < 0.001_r_kind) then
            cld_rbc_idx(i)=one   !clear/clear
       else
           cld_rbc_idx(i) = zero
       endif
     enddo 
  endif
end subroutine deter_cld_rbc_idx


subroutine obserr_allsky_mw(ichan,error0,tnoise,tnoise_cld,clwp_obs,clwp_guess,obstype)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:   obserr_allsky_mw 
!   prgmmr: mkim          org: np23                date: 2016-03-02
!
! abstract: define observation error for all-sky microwave radiance data 
!
! program history log:
!   2015-03-02  mkim : added observation error model for all-sky GMI radiance data 
!   2015-03-02  mkim : moved all-sky AMSU-A observation error model from setuprad.f90 
!
!   input argument list:
!     error0    - input observation error before applying changes for cloudy sky condition
!
!   output argument list:
!     error0      - observation error 
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: half,zero
  implicit none

  real(r_kind),    intent(inout) :: error0 
  real(r_kind),    intent(in   ) :: tnoise, tnoise_cld, clwp_obs,clwp_guess
  integer(i_kind), intent(in   ) :: ichan
  character(10),   intent(in   ) :: obstype 

! Declare local variables 
  real(r_kind) :: clwp_avg,cldx

if (obstype == 'amsua') then

   clwp_avg=half*(clwp_obs+clwp_guess)
   if(clwp_avg <= cclr_amsua(ichan)) then
      error0 = tnoise
   else if(clwp_avg > cclr_amsua(ichan) .and. clwp_avg < ccld_amsua(ichan)) then
      error0 = tnoise + &
              (clwp_avg-cclr_amsua(ichan))*(tnoise_cld-tnoise)/(ccld_amsua(ichan)-cclr_amsua(ichan))
   else
      error0 = tnoise_cld
   endif

else if (obstype == 'gmi') then
   clwp_avg=half*(clwp_obs+clwp_guess)

   if(clwp_avg .lt. cclr_gmi(ichan) ) then
      error0 =    tnoise
   else if (clwp_avg .ge.cclr_gmi(ichan) .and. clwp_avg .lt. Ccld1_gmi(ichan)) then
      error0 = tnoise1_gmi(ichan) + (tnoise- tnoise1_gmi(ichan))*(Ccld1_gmi(ichan)-clwp_avg)/(Ccld1_gmi(ichan)-cclr_gmi(ichan))
   else if (clwp_avg .ge. Ccld1_gmi(ichan) .and. clwp_avg .lt. Ccld2_gmi(ichan)) then
      error0 = tnoise_cld + (tnoise1_gmi(ichan) - tnoise_cld)*(Ccld2_gmi(ichan)-clwp_avg)/(Ccld2_gmi(ichan)-Ccld1_gmi(ichan))
   else
      error0 = tnoise_cld
   endif

else   ! temporary for other sensors' obs error models  
  
  clwp_avg=half*(clwp_obs+clwp_guess)
  if(clwp_avg < 0.05_r_kind) then
      error0 = tnoise
  else
      error0 = tnoise_cld
  endif

endif

end subroutine obserr_allsky_mw 

end module allsky_mod
