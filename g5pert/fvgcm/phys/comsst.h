C
C Sea-surface temperature values
C
c
c $Id$
c $Author$
c
      integer  pssttim,totsstsz
      parameter (pssttim=12,totsstsz=2000)

      common/sst/sstbdy(plon,plat,2), sst(plon,plat), 
     &             cdaysstm, cdaysstp
      common/sst/nm, np, sstid,lonsiz, levsiz, latsiz, timesiz, np1
      common/sst/date_sst,sec_sst


      real sstbdy     ! SST values on boundary dataset
      real sst        ! Interpolated model sst values
      real cdaysstm   ! Calendar day for prv. month SST values read in
      real cdaysstp   ! Calendar day for nxt. month SST values read in
      
      integer nm,np   ! Array indices for prv., nxt month sst data
      integer sstid   ! netcdf id for sst variable
      integer lonsiz  ! size of longitude dimension on sst dataset
      integer levsiz  ! size of level dimension on sst dataset
      integer latsiz  ! size of latitude dimension on sst dataset
      integer timesiz ! size of time dimension on sst dataset
      integer np1     ! current forward time index of sst dataset
      integer date_sst(totsstsz)! Date on sst dataset (YYYYMMDD)
      integer sec_sst(totsstsz) ! seconds of date on sst dataset (0-86399)

