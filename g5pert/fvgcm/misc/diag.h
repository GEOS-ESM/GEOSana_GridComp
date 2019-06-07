!
!    -----------------------------------------------------------
!     Refer to ../README file for how to add or delete a field.
!    -----------------------------------------------------------
!
! ... 2-D diagnostic fields
!
      integer  iALDIF
      integer  iALDIR
      integer  iASDIF
      integer  iASDIR           
      integer  iBMA           
      integer  iCAPEMX
      integer  iCLDHGH
      integer  iCLDLOW
      integer  iCLDMED
      integer  iCLDTOT
      integer  iCNVCLD
      integer  iEMSFC
      integer  iFLNS
      integer  iFLNSC
      integer  iFLNT
      integer  iFLNTC
      integer  iFSDS
      integer  iFSNS
      integer  iFSNSC
      integer  iFSNT
      integer  iFSNTC
      integer  iGWET
      integer  iGWET1
      integer  iH300 
      integer  iH500 
      integer  iHTLCL
      integer  iHTMMSE
      integer  iLHFX
      integer  iLWSH
      integer  iORO
      integer  iPARDIF
      integer  iPARDIR
      integer  iPBLH
      integer  iPRECC
      integer  iPRECL
      integer  iQ10M
      integer  iQ2M
      integer  iQFLX
      integer  iQPERT
      integer  iSHFX
      integer  iSLP
      integer  iSNOWH
      integer  iSOLIN
      integer  iSRFRAD
      integer  iSURFP
      integer  iT10M
      integer  iT200 
      integer  iT2M
      integer  iT850 
      integer  iTAUGWX
      integer  iTAUGWY
      integer  iTAUX
      integer  iTAUY
      integer  iTPERT
      integer  iTQ
      integer  iTRAD
      integer  iTSKIN
      integer  iU10M
      integer  iU200 
      integer  iU2M
      integer  iU850 
      integer  iUSTAR
      integer  iV10M
      integer  iV200 
      integer  iV2M
      integer  iV850 
      integer  iZ0H 
      integer  iZ0M 
      integer  iZMMB 
      integer  iZMPR 
      integer  iZPD 
!
! ... 3-D diagnostic fields
!

      integer  iCAPE
      integer  iCGS
      integer  iCLDLWP
      integer  iCLOUD
      integer  iCLOUDUP
      integer  iCMFDQ
      integer  iCMFDQR2
      integer  iCMFDT
      integer  iCMFDTR
      integer  iCMFETR
      integer  iCMFMC
      integer  iCMFMC2
      integer  iCONVCLD
      integer  iDCAFDT
      integer  iDIABDT
      integer  iDQCOND
      integer  iDQPBLCG
      integer  iDQRL
      integer  iDTCOND
      integer  iDTPBLCG
      integer  iDTV
      integer  iDUV
      integer  iDVV
      integer  iEFFCLD
      integer  iEVAPL
      integer  iH
      integer  iKVH
      integer  iKVM
      integer  iO3VMR
      integer  iOMEGA
      integer  iQ
      integer  iQRL
      integer  iQRS
      integer  iRAYFDT
      integer  iRELHUM
      integer  iRHCLR
      integer  iRNEVPDQ
      integer  iRNEVPDT
      integer  iSETLWP
      integer  iSTRATCLD
      integer  iT
      integer  iTTMGW
      integer  iU
      integer  iUTGW
      integer  iV
      integer  iVD01
      integer  iVTGW
      integer  iZMCME
      integer  iZMDLF
      integer  iZMDQ
      integer  iZMDQR
      integer  iZMDT
      integer  iZMDU
      integer  iZMED
      integer  iZMEPS
      integer  iZMEU
      integer  iZMEVP
      integer  iZMMD
      integer  iZMMU
      integer  iZMPFLX
      integer  iZMQL

      integer  pd2d     ! Total number of 2-D diagnostic fields
      integer  pd3d     ! Total number of 3-D diagnostic fields
      integer  pdiag    ! Total number of diagnostic fields
!
! ... 2-D diagnosis fields
!
      parameter (iALDIF    =            1)
      parameter (iALDIR    = iALDIF   + 1)
      parameter (iASDIF    = iALDIR   + 1)
      parameter (iASDIR    = iASDIF   + 1)
      parameter (iBMA      = iASDIR   + 1)
      parameter (iCAPEMX   = iBMA     + 1)
      parameter (iCLDHGH   = iCAPEMX  + 1)
      parameter (iCLDLOW   = iCLDHGH  + 1)
      parameter (iCLDMED   = iCLDLOW  + 1)
      parameter (iCLDTOT   = iCLDMED  + 1)
      parameter (iCNVCLD   = iCLDTOT  + 1)
      parameter (iEMSFC    = iCNVCLD  + 1)
      parameter (iFLNS     = iEMSFC   + 1)
      parameter (iFLNSC    = iFLNS    + 1)
      parameter (iFLNT     = iFLNSC   + 1)
      parameter (iFLNTC    = iFLNT    + 1)
      parameter (iFSDS     = iFLNTC   + 1)
      parameter (iFSNS     = iFSDS    + 1)
      parameter (iFSNSC    = iFSNS    + 1)
      parameter (iFSNT     = iFSNSC   + 1)
      parameter (iFSNTC    = iFSNT    + 1)
      parameter (iGWET     = iFSNTC   + 1)
      parameter (iGWET1    = iGWET    + 1)
      parameter (iH300     = iGWET1   + 1)
      parameter (iH500     = iH300    + 1)
      parameter (iHTLCL    = iH500    + 1)
      parameter (iHTMMSE   = iHTLCL   + 1)
      parameter (iLHFX     = iHTMMSE  + 1)
      parameter (iLWSH     = iLHFX    + 1)
      parameter (iORO      = iLWSH    + 1)
      parameter (iPARDIF   = iORO     + 1)
      parameter (iPARDIR   = iPARDIF  + 1)
      parameter (iPBLH     = iPARDIR  + 1)
      parameter (iPRECC    = iPBLH    + 1)
      parameter (iPRECL    = iPRECC   + 1)
      parameter (iQ10M     = iPRECL   + 1)
      parameter (iQ2M      = iQ10M    + 1)
      parameter (iQFLX     = iQ2M     + 1)
      parameter (iQPERT    = iQFLX    + 1)
      parameter (iSHFX     = iQPERT   + 1)
      parameter (iSLP      = iSHFX    + 1)
      parameter (iSNOWH    = iSLP     + 1)
      parameter (iSOLIN    = iSNOWH   + 1)
      parameter (iSRFRAD   = iSOLIN   + 1)
      parameter (iSURFP    = iSRFRAD  + 1)
      parameter (iT10M     = iSURFP   + 1)
      parameter (iT200     = iT10M    + 1)
      parameter (iT2M      = iT200    + 1)
      parameter (iT850     = iT2M     + 1)
      parameter (iTAUGWX   = iT850    + 1)
      parameter (iTAUGWY   = iTAUGWX  + 1)
      parameter (iTAUX     = iTAUGWY  + 1)
      parameter (iTAUY     = iTAUX    + 1)
      parameter (iTPERT    = iTAUY    + 1)
      parameter (iTQ       = iTPERT   + 1)
      parameter (iTRAD     = iTQ      + 1)
      parameter (iTSKIN    = iTRAD    + 1)
      parameter (iU10M     = iTSKIN   + 1)
      parameter (iU200     = iU10M    + 1)
      parameter (iU2M      = iU200    + 1)
      parameter (iU850     = iU2M     + 1)
      parameter (iUSTAR    = iU850    + 1)
      parameter (iV10M     = iUSTAR   + 1)
      parameter (iV200     = iV10M    + 1)
      parameter (iV2M      = iV200    + 1)
      parameter (iV850     = iV2M     + 1)
      parameter (iZ0H      = iV850    + 1)
      parameter (iZ0M      = iZ0H     + 1)
      parameter (iZMMB     = iZ0M     + 1)
      parameter (iZMPR     = iZMMB    + 1)
      parameter (iZPD      = iZMPR    + 1)

      parameter (pd2d      = iZPD )
!
! ... 3-D diagnosis fields
!
      parameter (iCAPE     = pd2d     + 1)
      parameter (iCGS      = iCAPE     + 1)
      parameter (iCLDLWP   = iCGS     + 1)
      parameter (iCLOUD    = iCLDLWP  + 1)
      parameter (iCLOUDUP  = iCLOUD   + 1)
      parameter (iCMFDQ    = iCLOUDUP + 1)
      parameter (iCMFDQR2  = iCMFDQ   + 1)
      parameter (iCMFDT    = iCMFDQR2 + 1)
      parameter (iCMFDTR   = iCMFDT   + 1)
      parameter (iCMFETR   = iCMFDTR  + 1)
      parameter (iCMFMC    = iCMFETR  + 1)
      parameter (iCMFMC2   = iCMFMC   + 1)
      parameter (iCONVCLD  = iCMFMC2  + 1)
      parameter (iDCAFDT   = iCONVCLD + 1)
      parameter (iDIABDT   = iDCAFDT  + 1)
      parameter (iDQCOND   = iDIABDT  + 1)
      parameter (iDQPBLCG  = iDQCOND  + 1)
      parameter (iDQRL     = iDQPBLCG + 1)
      parameter (iDTCOND   = iDQRL    + 1)
      parameter (iDTPBLCG  = iDTCOND  + 1)
      parameter (iDTV      = iDTPBLCG + 1)
      parameter (iDUV      = iDTV     + 1)
      parameter (iDVV      = iDUV     + 1)
      parameter (iEFFCLD   = iDVV     + 1)
      parameter (iEVAPL    = iEFFCLD  + 1)
      parameter (iH        = iEVAPL   + 1)
      parameter (iKVH      = iH       + 1)
      parameter (iKVM      = iKVH     + 1)
      parameter (iO3VMR    = iKVM     + 1)
      parameter (iOMEGA    = iO3VMR   + 1)
      parameter (iQ        = iOMEGA   + 1)
      parameter (iQRL      = iQ       + 1)
      parameter (iQRS      = iQRL     + 1)
      parameter (iRAYFDT   = iQRS     + 1)
      parameter (iRELHUM   = iRAYFDT  + 1)
      parameter (iRHCLR    = iRELHUM  + 1)
      parameter (iRNEVPDQ  = iRHCLR   + 1)
      parameter (iRNEVPDT  = iRNEVPDQ + 1)
      parameter (iSETLWP   = iRNEVPDT + 1)
      parameter (iSTRATCLD = iSETLWP  + 1)
      parameter (iT        = iSTRATCLD+ 1)
      parameter (iTTMGW    = iT       + 1)
      parameter (iU        = iTTMGW   + 1)
      parameter (iUTGW     = iU       + 1)
      parameter (iV        = iUTGW    + 1)
      parameter (iVD01     = iV       + 1)
      parameter (iVTGW     = iVD01    + 1)
      parameter (iZMCME    = iVTGW    + 1)
      parameter (iZMDLF    = iZMCME   + 1)
      parameter (iZMDQ     = iZMDLF   + 1)
      parameter (iZMDQR    = iZMDQ    + 1)
      parameter (iZMDT     = iZMDQR   + 1)
      parameter (iZMDU     = iZMDT    + 1)
      parameter (iZMED     = iZMDU    + 1)
      parameter (iZMEPS    = iZMED    + 1)
      parameter (iZMEU     = iZMEPS   + 1)
      parameter (iZMEVP    = iZMEU    + 1)
      parameter (iZMMD     = iZMEVP   + 1)
      parameter (iZMMU     = iZMMD    + 1)
      parameter (iZMPFLX   = iZMMU    + 1)
      parameter (iZMQL     = iZMPFLX  + 1)

      parameter (pdiag     = iZMQL)
      parameter (pd3d      = pdiag  - pd2d)

