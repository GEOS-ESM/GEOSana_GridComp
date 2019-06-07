      subroutine diaginit(qname, fldtbl, qunit, qdesc,    qvdim,  cdiag,  &
                          inx,   fldloc, pick,  tavg,     nslice, nrec2d, &
                          nrec,  n2d,    n3d,   diagattr, jfirst, jlast)

      use precision 
      use mod_comm, only : imr
      use mod_comm, only : nlayr => nl
      implicit none

#include <diag.h>


!
! ----------------------------------------------------------------------
!
!  The diagnostic fields are mostly produced by the physics routines,
!  therefore we adopted the concept that's used in CCM3, i.e., by 
!  declaring a huge buffer to hold all the diagnostic fields. Since
!  CCM3 physics is parallelized in latitude (J) direction, this buffer
!  is filled in parallel regions, latitude by latitude  by different
!  processors. To achieve computational efficiency, J dimension is
!  declared the outermost dimension of the buffer (diagbuf). The 2nd
!  dimension of diagbuf, maxslice, is the total number of slices for
!  all fields, where a slice is a 2-D array with dimension (imr,jfirst:jlast). 
!
!  The two alphabetically pre-sorted tables of names for 2-D and 3-D 
!  diagnostic fields (stored in a single array, qname, 2-D first) is
!  sorted together by a merge sort. The original name table is left
!  intact, the new sorted table is stored in fldtbl with the indices
!  (w.r.t. the original table) stored in inx. Note that the reason
!  why we created a globally (2-D and 3-D) sorted table was because
!  we tried to improve the simple hash-like indexing of diagnostic
!  field by an efficient search. Unfortunately that approach was not
!  successful when used in parallel region. This part of the code
!  should be revised.
!
!  The subroutine outfld is responsible for filling diagnostic fields
!  into the buffer. The original NCAR design was to use the first
!  argument, which is a 8-character long string, as the hash key to
!  find the location and length of the field in the buffer. We pack
!  similar info into an integer attribe array which is passed into
!  outfld as the first argument. These attributes are (in that order) 
!
!       pick : 1 if this field is selected in diag.tbl, 0 if not
!       tavg : 1 if this field is time averaged, 0 if not
!      qvdim : 1 for 2-D field, nlayr for 3-D field
!     fldloc : location of the field in the buffer in unit of slice
!      cdiag : number of instantaneous samples of the field, used
!              for time averaging
!
! ----------------------------------------------------------------------
!
!!      integer  maxslice
!!      parameter   (maxslice    = pd2d + pd3d * nlayr)
      integer, intent(in):: jfirst, jlast
!!      real(r4) diagbuf(imr,maxslice,jfirst:jlast) ! diagnostic field buffer
      integer  cdiag(pdiag)                ! counter for time averaging
      integer  inx(pdiag)                  ! sorted field table indices
      integer  qvdim(pdiag)                ! vertical dimension of field
      integer  fldloc(pdiag)               ! location of field in diagbuf
      integer  tavg(pdiag)                 ! flag for time-averaged fields (1=true)
      integer  pick(pdiag)                 ! flag for selected fields (1=true)
      integer  diagattr(5,pdiag)           ! field attributes
      integer  hist(0:26)                  ! field name histogram
      integer  nrec2d                      ! diagnostics output record number (surface)
      integer  nrec                        ! diagnostics output record number
      integer  nslice                      ! total records of selected diag. fields
      character*8  qname(pdiag)                ! prescribed field table
      character*8  fldtbl(pdiag)               ! sorted field table
      character*16 qunit(pdiag)                ! unit for diagnostic fields
      character*80 qdesc(pdiag)                ! description for diagnostic fields

      character     comment
      character*8   fld
      character*128 buf

      logical     done
      logical     picked
      logical     inorder

      integer findid
      integer iutbl
      integer fldid, id
      integer head, tail
      integer n2d
      integer n3d
      integer ndiag
      integer i, j, n
      integer ios
!
! ... Define comment character used in diag.tbl
!
!     If the leftmost non-blank character in diag.tbl matches this 
!     character, that line is treated as comments.
!
      comment = '!'
!
! ***********  Note on kf77 bug on Digital machine ***********
!
!  For some reason, kf77 cannot parse those qdesc line if //
!  appears at the end of the first line. When it is moved
!  to the beginning of the second line, the problem is gone.
!
! ************************************************************
!
! ... 2-D diagnostic fields (in alphabetic order)
!

      qname(iALDIF   ) = 'ALDIF   '
      qdesc(iALDIF   ) = 'Albedo: longwave, diffuse '
      qunit(iALDIF   ) = 'fraction'

      qname(iALDIR   ) = 'ALDIR   '
      qdesc(iALDIR   ) = 'Albedo: longwave, direct '
      qunit(iALDIR   ) = 'fraction'

      qname(iASDIF   ) = 'ASDIF   '
      qdesc(iASDIF   ) = 'Albedo: shortwave, diffuse '
      qunit(iASDIF   ) = 'fraction'

      qname(iASDIR   ) = 'ASDIR   '
      qdesc(iASDIR   ) = 'Albedo: shortwave, direct '
      qunit(iASDIR   ) = 'fraction'

      qname(iBMA     ) = 'BMA     '
      qdesc(iBMA     ) = 'Bulk moisture avaliability '
      qunit(iBMA     ) = 'fraction'

      qname(iCAPEMX  ) = 'CAPEMX  '
      qdesc(iCAPEMX  ) = 'Maximum CAPE'
      qunit(iCAPEMX  ) = 'J/Kg'

      qname(iCLDHGH  ) = 'CLDHGH  '
      qdesc(iCLDHGH  ) = 'Vertically-integrated, random overlap, '     &
                         // 'high cloud amount'
      qunit(iCLDHGH  ) = 'fraction'

      qname(iCLDLOW  ) = 'CLDLOW  '
      qdesc(iCLDLOW  ) = 'Vertically-integrated, random overlap, '     &
                         // 'low cloud amount'
      qunit(iCLDLOW  ) = 'fraction'

      qname(iCLDMED  ) = 'CLDMED  '
      qdesc(iCLDMED  ) = 'Vertically-integrated, random overlap, '     &
                         // 'mid cloud amount'
      qunit(iCLDMED  ) = 'fraction'

      qname(iCLDTOT  ) = 'CLDTOT  '
      qdesc(iCLDTOT  ) = 'Vertically-integrated, random overlap, '     &
                         // 'total cloud cover'
      qunit(iCLDTOT  ) = 'fraction'

      qname(iCNVCLD  ) = 'CNVCLD  '
      qdesc(iCNVCLD  ) = 'Random overlap total convective cloud amount'
      qunit(iCNVCLD  ) = 'fraction'

      qname(iEMSFC   ) = 'EMSFC   '
      qdesc(iEMSFC   ) = 'Bulk surface emissivity'
      qunit(iEMSFC   ) = 'fraction'

      qname(iFLNS    ) = 'FLNS    '
      qdesc(iFLNS    ) = 'Net longwave flux at surface'
      qunit(iFLNS    ) = 'W/m2'

      qname(iFLNSC   ) = 'FLNSC   '
      qdesc(iFLNSC   ) = 'Clearsky net longwave flux at surface'
      qunit(iFLNSC   ) = 'W/m2'

      qname(iFLNT    ) = 'FLNT    '
      qdesc(iFLNT    ) = 'Net longwave flux at top'
      qunit(iFLNT    ) = 'W/m2'

      qname(iFLNTC   ) = 'FLNTC   '
      qdesc(iFLNTC   ) = 'Clearsky net longwave flux at top'
      qunit(iFLNTC   ) = 'W/m2'

      qname(iFSDS    ) = 'FSDS    '
      qdesc(iFSDS    ) = 'Flux shortwave downwelling surface'
      qunit(iFSDS    ) = 'W/m2'

      qname(iFSNS    ) = 'FSNS    '
      qdesc(iFSNS    ) = 'Net solar flux at surface'
      qunit(iFSNS    ) = 'W/m2'

      qname(iFSNSC   ) = 'FSNSC   '
      qdesc(iFSNSC   ) = 'Clearsky net solar flux at surface '
      qunit(iFSNSC   ) = 'W/m2'

      qname(iFSNT    ) = 'FSNT    '
      qdesc(iFSNT    ) = 'Net solar flux at top'
      qunit(iFSNT    ) = 'W/m2'

      qname(iFSNTC   ) = 'FSNTC   '
      qdesc(iFSNTC   ) = 'Clearsky net solar flux at top'
      qunit(iFSNTC   ) = 'W/m2'

      qname(iGWET    ) = 'GWET    '
      qdesc(iGWET    ) = 'Root zone soil wetness'
      qunit(iGWET    ) = 'fraction'

      qname(iGWET1   ) = 'GWET1   '
      qdesc(iGWET1   ) = 'Top soil layer wetness'
      qunit(iGWET1   ) = 'fraction'

      qname(iH300    ) = 'H300    '
      qdesc(iH300    ) = '300 hPa Geopotential height'
      qunit(iH300    ) = 'm'

      qname(iH500    ) = 'H500    '
      qdesc(iH500    ) = '500 hPa Geopotential height'
      qunit(iH500    ) = 'm'

      qname(iHTLCL   ) = 'HTLCL   '
      qdesc(iHTLCL   ) = ' Height above surface at LCL level'
      qunit(iHTLCL   ) = 'm'

      qname(iHTMMSE  ) = 'HTMMSE  '
      qdesc(iHTMMSE  ) = 'Height above surface at maximum moist static energy level'
      qunit(iHTMMSE  ) = 'm'

      qname(iLHFX    ) = 'LHFX    '
      qdesc(iLHFX    ) = 'Surface laten heat flux'
      qunit(iLHFX    ) = 'W/m2'

      qname(iLWSH    ) = 'LWSH    '
      qdesc(iLWSH    ) = 'Liquid water scale height'
      qunit(iLWSH    ) = 'm'

      qname(iORO     ) = 'ORO     '
      qdesc(iORO     ) = 'Surface type flag'
      qunit(iORO     ) = 'flag'

      qname(iPARDIF  ) = 'PARDIF  '
      qdesc(iPARDIF  ) = 'Diffuse photosynthetically active radiation'          &
                       // ' (0.35-0.70 um)'
      qunit(iPARDIF  ) = 'W/m2'

      qname(iPARDIR  ) = 'PARDIR  '
      qdesc(iPARDIR  ) = 'Direct photosynthetically active radiation'           &
                       //' (0.35-0.70 um)'
      qunit(iPARDIR  ) = 'W/m2'

      qname(iPBLH    ) = 'PBLH    '
      qdesc(iPBLH    ) = 'Planetary boundary layer height'
      qunit(iPBLH    ) = 'm'

      qname(iPRECC   ) = 'PRECC   '
      qdesc(iPRECC   ) = 'Convective precipitation rate'
      qunit(iPRECC   ) = 'mm/day'

      qname(iPRECL   ) = 'PRECL   '
      qdesc(iPRECL   ) = 'Large-scale precipitation rate'
      qunit(iPRECL   ) = 'mm/day'

      qname(iQ10M    ) = 'Q10M    '
      qdesc(iQ10M    ) = '10 Specific humidity'
      qunit(iQ10M    ) = 'Kg/Kg'

      qname(iQ2M     ) = 'Q2M     '
      qdesc(iQ2M     ) = '2 Specific humidity'
      qunit(iQ2M     ) = 'Kg/Kg'

      qname(iQFLX    ) = 'QFLX    '
      qdesc(iQFLX    ) = 'Surface water flux'
      qunit(iQFLX    ) = 'Kg/m2/s'

      qname(iQPERT   ) = 'QPERT   '
      qdesc(iQPERT   ) = 'Perturbation specific humidity '      &
                         // '(eddies in PBL)'
      qunit(iQPERT   ) = 'Kg/Kg'

      qname(iSHFX    ) = 'SHFX    '
      qdesc(iSHFX    ) = 'Surface sensible heat flux'
      qunit(iSHFX    ) = 'W/m2'

      qname(iSLP     ) = 'SLP     '
      qdesc(iSLP     ) = 'Sea level pressure'
      qunit(iSLP     ) = 'Pa'

      qname(iSNOWH   ) = 'SNOWH   '
      qdesc(iSNOWH   ) = 'Water equivalent snow depth'
      qunit(iSNOWH   ) = 'm'

      qname(iSOLIN   ) = 'SOLIN   '
      qdesc(iSOLIN   ) = 'Solar insolation'
      qunit(iSOLIN   ) = 'W/m2'

      qname(iSRFRAD  ) = 'SRFRAD  '
      qdesc(iSRFRAD  ) = 'Net radiative flux at surface'
      qunit(iSRFRAD  ) = 'W/m2'

      qname(iSURFP   ) = 'SURFP   '
      qdesc(iSURFP   ) = 'Surface pressure'
      qunit(iSURFP   ) = 'Pa'

      qname(iT10M    ) = 'T10M    '
      qdesc(iT10M    ) = '10 meter temperature'
      qunit(iT10M    ) = 'K'

      qname(iT200    ) = 'T200    '
      qdesc(iT200    ) = '200 hPa temperature'
      qunit(iT200    ) = 'K'

      qname(iT2M     ) = 'T2M     '
      qdesc(iT2M     ) = '2 meter temperature'
      qunit(iT2M     ) = 'K'

      qname(iT850    ) = 'T850    '
      qdesc(iT850    ) = '850 hPa temperature'
      qunit(iT850    ) = 'K'

      qname(iTAUGWX  ) = 'TAUGWX  '
      qdesc(iTAUGWX  ) = 'East-west gravity wave drag surface stress'
      qunit(iTAUGWX  ) = 'N/m2'

      qname(iTAUGWY  ) = 'TAUGWY  '
      qdesc(iTAUGWY  ) = 'North-south gravity wave drag surface stress'
      qunit(iTAUGWY  ) = 'N/m2'

      qname(iTAUX    ) = 'TAUX    '
      qdesc(iTAUX    ) = 'X-component (east-west) of surface stress'
      qunit(iTAUX    ) = 'N/m2'

      qname(iTAUY    ) = 'TAUY    '
      qdesc(iTAUY    ) = 'Y-component (north-south) of surface stress'
      qunit(iTAUY    ) = 'N/m2'

      qname(iTPERT   ) = 'TPERT   '
      qdesc(iTPERT   ) = 'Perturbation temperature (eddies in PBL)'
      qunit(iTPERT   ) = 'K'

      qname(iTQ      ) = 'TQ      '
      qdesc(iTQ      ) = 'Total precipitable water'
      qunit(iTQ      ) = 'Kg/m2'

      qname(iTRAD    ) = 'TRAD    '
      qdesc(iTRAD    ) = 'Surface brightness temperature'
      qunit(iTRAD    ) = 'K'

      qname(iTSKIN   ) = 'TSKIN   '
      qdesc(iTSKIN   ) = 'Surface skin temperature'
      qunit(iTSKIN   ) = 'K'

      qname(iU10M    ) = 'U10M    '
      qdesc(iU10M    ) = '10 meter U wind'
      qunit(iU10M    ) = 'm/s'

      qname(iU200    ) = 'U200    '
      qdesc(iU200    ) = '200 hPa U wind'
      qunit(iU200    ) = 'm/s'

      qname(iU2M     ) = 'U2M     '
      qdesc(iU2M     ) = '2 meter U wind'
      qunit(iU2M     ) = 'm/s'

      qname(iU850    ) = 'U850    '
      qdesc(iU850    ) = '850 hPa U wind'
      qunit(iU850    ) = 'm/s'

      qname(iUSTAR   ) = 'USTAR   '
      qdesc(iUSTAR   ) = 'Surface friction velocity'
      qunit(iUSTAR   ) = 'm/s'

      qname(iV10M    ) = 'V10M    '
      qdesc(iV10M    ) = '10 meter V wind'
      qunit(iV10M    ) = 'm/s'

      qname(iV200    ) = 'V200    '
      qdesc(iV200    ) = '200 hPa V wind'
      qunit(iV200    ) = 'm/s'

      qname(iV2M     ) = 'V2M     '
      qdesc(iV2M     ) = '2 meter V wind'
      qunit(iV2M     ) = 'm/s'
!
      qname(iV850    ) = 'V850    '
      qdesc(iV850    ) = '850 hPa V wind'
      qunit(iV850    ) = 'm/s'

      qname(iZ0H     ) = 'Z0H     '
      qdesc(iZ0H     ) = 'Roughness length , sensible heat'
      qunit(iZ0H     ) = 'm'

      qname(iZ0M     ) = 'Z0M     '
      qdesc(iZ0M     ) = 'Roughness length , momentum'
      qunit(iZ0M     ) = 'm'

      qname(iZMMB    ) = 'ZMMB    '
      qdesc(iZMMB    ) = 'Cloud base mass flux from Z&M scheme'
      qunit(iZMMB    ) = 'kg/m2/s'

      qname(iZMPR    ) = 'ZMPR    '
      qdesc(iZMPR    ) = 'Precipitation from Z&M scheme'
      qunit(iZMPR    ) = 'mm/day'

      qname(iZPD     ) = 'ZPD     '
      qdesc(iZPD     ) = 'Displacement height'
      qunit(iZPD     ) = 'm'

! ... 3-D diagnostic fields (in alphabetic order)
!
      qname(iCAPE    ) = 'CAPE    '
      qdesc(iCAPE    ) = 'CAPE'
      qunit(iCAPE    ) = 'J/Kg'

      qname(iCGS     ) = 'CGS     '
      qdesc(iCGS     ) = 'Counter-gradient coefficient on surface '      &
                         // 'kinematic fluxes'
      qunit(iCGS     ) = 's/m2'
 
      qname(iCLDLWP  ) = 'CLDLWP  '
      qdesc(iCLDLWP  ) = 'Actual cloud liquid water path length'
      qunit(iCLDLWP  ) = 'gram/m2'

      qname(iCLOUD   ) = 'CLOUD   '
      qdesc(iCLOUD   ) = 'Cloud fraction'
      qunit(iCLOUD   ) = 'fraction'

      qname(iCLOUDUP ) = 'CLOUDUP '
      qdesc(iCLOUDUP ) = 'Cloud fraction during omega < 0.'
      qunit(iCLOUDUP ) = 'fraction'
 
      qname(iCMFDQ   ) = 'CMFDQ   '
      qdesc(iCMFDQ   ) = 'q tendency - Hack convection'
      qunit(iCMFDQ   ) = 'Kg/Kg/s'

      qname(iCMFDQR2 ) = 'CMFDQR2 '
      qdesc(iCMFDQR2 ) = 'Rain production rate - Hack convection'
      qunit(iCMFDQR2 ) = 'Kg/Kg/s'

      qname(iCMFDT   ) = 'CMFDT   '
      qdesc(iCMFDT   ) = 'T tendency - Hack convection'
      qunit(iCMFDT   ) = 'K/s'

      qname(iCMFDTR  ) = 'CMFDTR  '
      qdesc(iCMFDTR  ) = 'Detrainment mass flux - Hack convection'
      qunit(iCMFDTR  ) = 'Pa/s'

      qname(iCMFETR  ) = 'CMFETR  '
      qdesc(iCMFETR  ) = 'Entrainment mass flux - Hack convection'
      qunit(iCMFETR  ) = 'Pa/s'

      qname(iCMFMC   ) = 'CMFMC   '
      qdesc(iCMFMC   ) = 'Total Moist convection mass flux'
      qunit(iCMFMC   ) = 'Pa/s'

      qname(iCMFMC2  ) = 'CMFMC2  '
      qdesc(iCMFMC2  ) = 'Hack Moist convection mass flux'
      qunit(iCMFMC2  ) = 'Pa/s'

      qname(iCONVCLD  ) = 'CONVCLD  '
      qdesc(iCONVCLD  ) = 'Convective cloud amount'
      qunit(iCONVCLD  ) = 'fraction'

      qname(iDCAFDT  ) = 'DCAFDT  '
      qdesc(iDCAFDT  ) = 'T Tendency - Dry convective adjustment'
      qunit(iDCAFDT  ) = 'K/s'

      qname(iDIABDT  ) = 'DIABDT  '
      qdesc(iDIABDT  ) = 'T Tendency - Total adiabatic (physics)'
      qunit(iDIABDT  ) = 'K/s'                   

      qname(iDQCOND  ) = 'DQCOND  '
      qdesc(iDQCOND  ) = 'Q tendency - moist physics'
      qunit(iDQCOND  ) = 'kg/kg/s'

      qname(iDQPBLCG ) = 'DQPBLCG '
      qdesc(iDQPBLCG ) = 'Q tendency - PBL counter gradient'
      qunit(iDQPBLCG ) = 'kg/kg/s'

      qname(iDQRL    ) = 'DQRL    '
      qdesc(iDQRL    ) = 'Rain production rate - large-scale'
      qunit(iDQRL    ) = 'kg/kg/s'

      qname(iDTCOND  ) = 'DTCOND  '
      qdesc(iDTCOND  ) = 'T tendency - moist physics'
      qunit(iDTCOND  ) = 'K/s'

      qname(iDTPBLCG ) = 'DTPBLCG '
      qdesc(iDTPBLCG ) = 'T tendency - PBL counter gradient'
      qunit(iDTPBLCG ) = 'K/s'
 
      qname(iDTV     ) = 'DTV     '
      qdesc(iDTV     ) = 'T vertical diffusion'
      qunit(iDTV     ) = 'K/s'

      qname(iDUV     ) = 'DUV     '
      qdesc(iDUV     ) = 'U tendency from vertical diffusion'
      qunit(iDUV     ) = 'm/s2'

      qname(iDVV     ) = 'DVV     '
      qdesc(iDVV     ) = 'V tendency from vertical diffusion'
      qunit(iDVV     ) = 'm/s2'

      qname(iEFFCLD  ) = 'EFFCLD  '
      qdesc(iEFFCLD  ) = 'Effective cloud fraction'
      qunit(iEFFCLD  ) = 'fraction'

      qname(iEVAPL   ) = 'EVAPL   '
      qdesc(iEVAPL   ) = 'Large-scale evaporation'
      qunit(iEVAPL   ) = 'kg/kg/s'

      qname(iH       ) = 'H       '
      qdesc(iH       ) = 'Geopotential height'
      qunit(iH       ) = 'm'

      qname(iKVH     ) = 'KVH     '
      qdesc(iKVH     ) = 'Vertical diffusion diffusivities '      &
                         // '(heat/moisture)'
      qunit(iKVH     ) = 'm2/s'

      qname(iKVM     ) = 'KVM     '
      qdesc(iKVM     ) = 'Eddy diffusivity for momentum'
      qunit(iKVM     ) = 'm2/s'

      qname(iO3VMR   ) = 'O3VMR   '
      qdesc(iO3VMR   ) = 'Ozone volume mixing ratio'
      qunit(iO3VMR   ) = 'Kg/Kg'

      qname(iOMEGA   ) = 'OMEGA   '
      qdesc(iOMEGA   ) = 'Vertical pressure velocity'
      qunit(iOMEGA   ) = 'Pa/s'

      qname(iQ       ) = 'Q       '
      qdesc(iQ       ) = 'Specified humidity'
      qunit(iQ       ) = 'Kg/Kg'

      qname(iQRL     ) = 'QRL     '
      qdesc(iQRL     ) = 'Longwave heating rate'
      qunit(iQRL     ) = 'K/s'

      qname(iQRS     ) = 'QRS     '
      qdesc(iQRS     ) = 'Shortwave heating rate'
      qunit(iQRS     ) = 'K/s'

      qname(iRAYFDT  ) = 'RAYFDT  '
      qdesc(iRAYFDT  ) = 'T Tendency - Rayleigh friction'
      qunit(iRAYFDT  ) = 'K/s'

      qname(iRELHUM  ) = 'RELHUM  '
      qdesc(iRELHUM  ) = 'Relative Humidity after cloud physics'
      qunit(iRELHUM  ) = '%'

      qname(iRHCLR   ) = 'RHCLR   '
      qdesc(iRHCLR   ) = 'Relative Humidity in clear region'
      qunit(iRHCLR   ) = '%'

      qname(iRNEVPDQ ) = 'RNEVPDQ '
      qdesc(iRNEVPDQ ) = 'Q Tendency - Rain evaporation'
      qunit(iRNEVPDQ ) = 'Kg/Kg/s'

      qname(iRNEVPDT ) = 'RNEVPDT '
      qdesc(iRNEVPDT ) = 'T Tendency - Rain evaporation'
      qunit(iRNEVPDT ) = 'K/s'         

      qname(iSETLWP  ) = 'SETLWP  '
      qdesc(iSETLWP  ) = 'Specified liquid water path lengths'
      qunit(iSETLWP  ) = 'gram/m2'

      qname(iSTRATCLD) = 'STRATCLD'
      qdesc(iSTRATCLD) = 'Stratiform cloud amount'
      qunit(iSTRATCLD) = 'fraction'

      qname(iT       ) = 'T       '
      qdesc(iT       ) = 'Temperature'
      qunit(iT       ) = 'K'

      qname(iTTMGW   ) = 'TTMGW   '
      qdesc(iTTMGW   ) = 'T tendency - gravity wave drag'
      qunit(iTTMGW   ) = 'K/s'

      qname(iU       ) = 'U       '
      qdesc(iU       ) = 'U wind'
      qunit(iU       ) = 'm/s'

      qname(iUTGW    ) = 'UTGW    '
      qdesc(iUTGW    ) = 'U tendency - gravity wave drag'
      qunit(iUTGW    ) = 'm/s2'

      qname(iV       ) = 'V       '
      qdesc(iV       ) = 'V wind'
      qunit(iV       ) = 'm/s'

      qname(iVD01    ) = 'VD01    '
      qdesc(iVD01    ) = 'Vertical diffusion tendency of water vapor'
      qunit(iVD01    ) = 'Kg/Kg/s'

      qname(iVTGW    ) = 'VTGW    '
      qdesc(iVTGW    ) = 'V tendency - gravity wave drag'
      qunit(iVTGW    ) = 'm/s2'

      qname(iZMCME   ) = 'ZMCME   '
      qdesc(iZMCME   ) = 'Condensation - evaporation from Z&M scheme'
      qunit(iZMCME   ) = 'kg/kg/s'

      qname(iZMDLF   ) = 'ZMDLF   '
      qdesc(iZMDLF   ) = 'Detrainment of cloud water from Z&M scheme'
      qunit(iZMDLF   ) = 'kg/kg/s'

      qname(iZMDQ    ) = 'ZMDQ    '
      qdesc(iZMDQ    ) = 'q tendency - Zhang-McFarlane convection'
      qunit(iZMDQ    ) = 'kg/kg/s'

      qname(iZMDQR   ) = 'ZMDQR   '
      qdesc(iZMDQR   ) = 'Rain production rate - Z&M convection'
      qunit(iZMDQR   ) = 'kg/kg/s'

      qname(iZMDT    ) = 'ZMDT    '
      qdesc(iZMDT    ) = 'T tendency - Zhang-McFarlane convection'
      qunit(iZMDT    ) = 'K/s'

      qname(iZMDU    ) = 'ZMDU    '
      qdesc(iZMDU    ) = 'Updraft detrainment mass flux - Z&M convection'
      qunit(iZMDU    ) = 'Pa/s'

      qname(iZMED    ) = 'ZMED    '
      qdesc(iZMED    ) = 'Downdraft entrainment mass flux - Z&M convection'
      qunit(iZMED    ) = 'Pa/s'

      qname(iZMEPS   ) = 'ZMEPS   '
      qdesc(iZMEPS   ) = 'Fractional entrainment - Z&M convection'
      qunit(iZMEPS   ) = '1/s'

      qname(iZMEU    ) = 'ZMEU    '
      qdesc(iZMEU    ) = 'Updraft entrainment mass flux - Z&M convection'
      qunit(iZMEU    ) = 'Pa/s'

      qname(iZMEVP   ) = 'ZMEVP   '
      qdesc(iZMEVP   ) = 'downdraft evaporation from Z&M convection'
      qunit(iZMEVP   ) = 'kg/kg/s'

      qname(iZMMD    ) = 'ZMMD    '
      qdesc(iZMMD    ) = 'Downdraft mass flux - Z&M convection'
      qunit(iZMMD    ) = 'Pa/s'

      qname(iZMMU    ) = 'ZMMU    '
      qdesc(iZMMU    ) = 'Updraft mass flux - Z&M convection'
      qunit(iZMMU    ) = 'pa/s'

      qname(iZMPFLX  ) = 'ZMPFLX  '
      qdesc(iZMPFLX  ) = 'Precipitation flux - Z&M convection'
      qunit(iZMPFLX  ) = 'kg/m2/s'

      qname(iZMQL    ) = 'ZMQL    '
      qdesc(iZMQL    ) = 'Cloud water in updraft - Z&M convection'
      qunit(iZMQL    ) = 'kg/kg'
!
! ... Convert all field names to upper case (paranoia check)
!
      do n = 1, pdiag
        call upper(qname(n)) 
      enddo
!
! ... Are 2D field names in alphabetic order?
!
      inorder = .true.
      n = 1
      do while (inorder .and. n .lt. pd2d) 
        if (qname(n) .gt. qname(n+1)) then
          write (6,*) '[diaginit]: 2D qname is not in alphabetic order.'
          write (6,'(13x,a11,i3,a4,a8)')                 &
                '---> qname(', n, ') = ', qname(n)
          write (6,'(13x,a11,i3,a4,a8)')                 &
                '---> qname(', n+1, ') = ', qname(n+1)
          inorder = .false.
        endif
        n = n +1
      enddo
      if (.not. inorder) stop
!
! ... Are 3D field names in alphabetic order?
!
      inorder = .true.
      n = pd2d + 1
      do while (inorder .and. n .lt. pdiag) 
        if (qname(n) .gt. qname(n+1)) then
          write (6,*) '[diaginit]: 3D qname is not in alphabetic order.'
          write (6,'(13x,a11,i3,a4,a8)')               &
                '---> qname(', n, ') = ', qname(n)
          write (6,'(13x,a11,i3,a4,a8)')               &
                '---> qname(', n+1, ') = ', qname(n+1)
          inorder = .false.
        endif
        n = n +1
      enddo
      if (.not. inorder) stop
!
! ... Vertical dimension for each diagnostic field
!
      do n = 1, pd2d
        qvdim(n) = 1
      enddo

      do n = pd2d+1, pdiag
        qvdim(n) = nlayr
      enddo
!
! ... Set initial values
!
      nrec = 1
      nrec2d = 1

      do n = 1, pdiag
        pick(n)  = 0
        tavg(n)  = 1
! SJL 11/07/98
        cdiag(n) = 0
      enddo

!!!$omp parallel do private(i,j,n)
!!      do j = jfirst, jlast
!!        do n = 1, maxslice
!!          do i = 1, imr
!!            diagbuf(i,n,j)  = 0.
!!           enddo
!!        enddo
!!      enddo


!
! ... Sort qname without physically rearrange qname, the 
!     sorted array is stored in fldtbl and the indices
!     of the sorted table (w.r.t. qname) is stored in inx.
!
      call merge (pd2d, pd3d, qname, fldtbl, inx)
!
! ... Get histogram of field names (for fast searching)
!
      do n = 0, 26
        hist(n) = 0
      enddo

      do n = 1, pdiag
        i = ichar(fldtbl(n)(1:1)) - ichar('A') + 1
        hist(i) = hist(i) + 1
      enddo

      do n = 1, 26
        hist(n) = hist(n) + hist(n-1)
      enddo
!     print '(1x,26(1x,a2))', (char(ichar('A')+n-1), n=1,26)
!     print '(26(i3)//)', (hist(n), n=1,26)

      do n = 0, 26
        if (hist(n) .eq. 0) hist(n) = 1
      enddo
!
! ... Read user input diagnostic field table, entries in the user
!     field table does not have to be in order. Therefore if there
!     are duplicated items, the setting for the entry encountered
!     later will precede the previous ones.
!
      iutbl = 70
      open (iutbl, file='diag.tbl', form='formatted', status='unknown') 

      n2d = 0
      n3d = 0
      done = .false.

      do while (.not. done)
        read (iutbl, '(a)', iostat=ios) buf
        if (ios .gt. 0) then
          print *,                        &
             '[diaginit]: Error in reading diagnostic field table.'
          stop
        else if (ios .lt. 0) then    ! EOF
          done = .true.
        else 
          call trimleft (buf)
          if (buf(1:1) .ne. comment) then
            read (buf, '(a8,2x,l4)') fld, picked
            call upper (fld)
            n = ichar(fld(1:1)) - ichar('A') + 1
            head = hist(n-1)
            tail = hist(n)
            if (head .eq. tail) then
              fldid = 0
              id = 0
            else
              fldid = findid (pdiag, fldtbl, fld, head, tail)
              id = inx(fldid)
            endif
            if (id .gt. 0 .and. id .le. pdiag) then
              if (picked) then
                pick(id) = 1
                if (qvdim(id) .eq. 1) then
                  n2d = n2d + 1
                else 
                  n3d = n3d + 1
                endif
              endif
            endif
          endif
        endif
      enddo
      close (iutbl)

      ndiag = n2d + n3d

      if (n2d .gt. pd2d .or. n3d .gt. pd3d) then
        write(6,*) '(diaginit): Number of diag fields exceeds capacity.'
        write(6,'(a,i4)') '            pd2d has to be >= ', n2d
        write(6,'(a,i4)') '            pd3d has to be >= ', n3d
        stop
      endif

!
! ... Build table for selected field location in diagbuf
!
      nslice = 0
      do n = 1, pdiag
        if (pick(n) .eq. 1) then
          fldloc(n) = nslice + 1
          nslice = nslice + qvdim(n)
        else
          fldloc(n) = 0
        endif
      enddo

!     print *, 'nslice = ', nslice
!
! ... Encode field attributes
!
      do n = 1, pdiag
        diagattr(1,n) = pick(n)
        diagattr(2,n) = tavg(n)
        diagattr(3,n) = qvdim(n)
        diagattr(4,n) = fldloc(n)
        diagattr(5,n) = cdiag(n)
      enddo

!      do n = 1, pdiag
!        write(6,'(a8,1x,i3,1x,i1,1x,i10)')         &
!                  qname(n), qvdim(n), pick(n), fldloc(n)
!      enddo
!      write(*,*) 'nslice *****', nslice
 
      return
      end

!***************************
      subroutine init_diagbuf(diagbuf, imr, nslice, jfirst, jlast)
      implicit none

      real*4 diagbuf(imr, nslice, jfirst:jlast)
      integer imr, nslice, jfirst, jlast
      integer i, j, n

!$omp parallel do private(i,j,n)
      do j = jfirst, jlast
        do n = 1, nslice
          do i = 1, imr
            diagbuf(i,n,j)  = 0.
           enddo
        enddo
      enddo

      return
      end 
