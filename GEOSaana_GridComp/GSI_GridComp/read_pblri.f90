!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: gsimod  ---

!
! !INTERFACE:
!
  module read_pblri

! !USES:

      use kinds, only: r_kind,r_double,i_kind
      use constants, only: zero,one_tenth,one,deg2rad,rad2deg,three,r_missing
      use gridmod, only: diagnostic_reg,regional,nlon,nlat,&
           tll2xy,txy2ll,rotate_wind_ll2xy,rotate_wind_xy2ll,&
           rlats,rlons
      use convinfo, only: nconvtype,ctwind, &
           icuse,ictype,ioctype
      use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen
      use obsmod, only: iadate,offtime_data,bmiss
      use deter_sfc_mod, only: deter_sfc2
      use mpimod, only: npe

      implicit none

      private

! !PUBLIC ROUTINES:
!     public read_pblri

!  interface read_pblri
!    module procedure read_pblh_prepfits
!    module procedure read_pblh_wp_text
!    module procedure read_pblh_calipso
!    module procedure read_pblh_gnssro
!  end interface

      public read_pblri_prepfits
      public read_pblri_wp_text
      public read_pblri_calipso_cats_mplnet_raob_wp
      public read_pblri_gnssro
!---------------------------------------------------------------------------

   CONTAINS

      subroutine read_pblri_prepfits(nread,ndata,nodata,infile,obstype,lunout,twindin,&
         sis,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_pblri     read obs from msgs in PREPFITS files (rpf == read aircraft)
!
! program history log:
!   2009-06     whiting - coding rpf
!   2009-10-20    zhu   - modify rpf for reading in pblh data in GSI
!   2009-10-21  whiting - modify cnem & pblriob for reading Caterina's files
!   2013-01-26  parrish - change from grdcrd to grdcrd1 (to allow successful debug compile on WCOSS)
!   2015-02-23  Rancic/Thomas - add l4densvar to time window logical
!   2015-10-01  guo     - consolidate use of ob location (in deg)
!
!   input argument list:
!     infile   - unit from which to read BUFR data
!     obstype  - observation type to process
!     lunout   - unit to which to write data for further processing
!     prsl_full- 3d pressure on full domain grid
!
!   output argument list:
!     nread    - number of type "obstype" observations read
!     nodata   - number of individual "obstype" observations read
!     ndata    - number of type "obstype" observations retained for further processing
!     twindin  - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
      implicit none 

!     Declare passed variables
      character(len=*),intent(in):: infile,obstype
      character(20),intent(in):: sis
      integer(i_kind),intent(in):: lunout
      integer(i_kind),intent(inout):: nread,ndata,nodata
      integer(i_kind),dimension(npe),intent(inout):: nobs
      real(r_kind),intent(in):: twindin

!     Declare local parameters
      integer(i_kind),parameter:: MXNM=25                 ! max Nems, max Replications
      integer(i_kind),parameter:: MXRP=255                ! max Nems, max Replications
      real(r_kind),parameter:: r360 = 360.0_r_kind

      integer(i_kind) lunin,msgt,ICOMP,idate,nlevp,nlevo,nlevc
      integer(i_kind) iout,nc,nr,nmsg,nrtyp,nrtmax,nchanl
      integer(i_kind) msub,nmsub,ireadsb,ireadmg
      character(80) cnem
      character(8)  ctyp, ctyp0

      character(8) cval(5)
      real(r_kind) rval(5)
      equivalence (cval(1),rval(1))

      real(r_kind) hdr(MXNM)
!     real(r_kind) plv(MXNM,MXRP), olv(MXNM,MXRP,MXRP), &
!          slv(mxnm,mxrp,mxrp)
      real(r_kind) clv(MXNM,MXRP,MXRP)
      real(r_kind),allocatable,dimension(:,:):: cdata_all

      character(10) date
      logical first,outside,inflate_error
      integer(i_kind) ihh,idd,iret,im,iy
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon
      integer(i_kind) pblriqm,i,maxobs,j,idomsfc
      integer(i_kind) ntest
!     integer(i_kind),dimension(8):: obs_time,anal_time,lobs_time
!     real(r_kind) ltime,cenlon_tmp
      real(r_kind) usage
      real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
      real(r_kind) pblriob,pblrioe,pblrielev,pblbak
      real(r_kind) dlat,dlon,dlat_earth,dlon_earth,stnelev
      real(r_kind) dlat_earth_deg,dlon_earth_deg
      real(r_kind) :: tsavg,ff10,sfcr,zz
!     real(r_kind),dimension(5):: tmp_time

      integer(i_kind) idate5(5),minobs,minan
      real(r_kind) time_correction,timeobs,time,toff,t4dv,zeps

      data lunin /50/

!     Initialize variables
      nreal=14
      ntest=0
      nrtmax=0                       ! # rpts to print per msg type (0=all)

      call closbf(lunin)
      open(lunin,file=trim(infile),form='unformatted')
      call mesgbc(lunin,msgt,icomp)
      call openbf(lunin,'IN',lunin)
      call datelen(10)                          ! set dates to 8-digit

      maxobs=0
      ctyp0=' '
      first=.true.
      do while( ireadmg(lunin,ctyp,idate).EQ.0 )
         if ( ctyp .ne. ctyp0 ) then
           ctyp0=ctyp
           nrtyp=0                        ! counter - rpts per type
         else ! ctyp = ctyp0
           if ( nrtyp .ge. nrtmax .and. nrtmax.ne.0 ) then
             cycle
           endif ! nrtyp >= nrtmax
         endif ! ctyp != ctyp0

         if ( ctyp(1:6).ne.'ADPUPA' .and.  &
              ctyp(1:6).ne.'AIRCAR' .and.  &
              ctyp(1:6).ne.'PROFLR' .and.  &
              ctyp(1:6).ne.'AIRCFT' ) cycle

!        Time offset
         if (first) then
            call time_4dvar(idate,toff)
            first=.false.
         end if

         do while (ireadsb(lunin) .eq. 0)
            nrtyp=nrtyp+1
            cnem='SID XOB YOB DHR ELV TYP T29 ITP'
            call ufbint(lunin,hdr,MXNM,1,iret,cnem)
            if ( iret.ne.1 ) write(*,*)'ERROR - ufbseq(HEADR) iret=',iret
            kx=hdr(6)
            if(kx == 431 .or. kx == 531) nkx= 131
            if(kx == 433 .or. kx == 533) nkx= 133
            if(kx == 435 .or. kx == 535) nkx= 135
            if(kx == 120) nkx= 120
            if(kx == 227) nkx= 181

            loop_convinfo_test: do nc=1,nconvtype
              if (trim(ioctype(nc)) /= trim(obstype))cycle loop_convinfo_test
              if (nkx == ictype(nc)) then
                 ikx=nc
                 maxobs=maxobs+1
              end if
            end do loop_convinfo_test
         end do ! while ireadsb
      end do ! while ireadmg

      allocate(cdata_all(nreal,maxobs))
      nread=0
      nchanl=0
      ilon=2
      ilat=3
      call closbf(lunin)
      open(lunin,file=trim(infile),form='unformatted')
      call mesgbc(lunin,msgt,icomp)
      call openbf(lunin,'IN',lunin)
      call datelen(10)                          ! set dates to 8-digit

      nr=0                           ! rpt counter
      nmsg=0                         ! msg counter

      ctyp0=' '
      do while( ireadmg(lunin,ctyp,idate).eq.0 )
      nmsg=nmsg+1
      msub = nmsub(lunin)

      if ( ctyp .ne. ctyp0 ) then 
        if ( ctyp0 .ne. " " ) write(*,*) ! linefeed
        ctyp0=ctyp
        write(*,'(/,a,1x,i10,$)')'  new ctyp="'//ctyp//'" idate=',idate
        write(*,'(3x,a,1x,i3,$)') 'msub=',msub
        nrtyp=0                        ! counter - rpts per type
      else ! ctyp = ctyp0
        if ( nrtyp .ge. nrtmax .and. nrtmax.ne.0 ) then 
          cycle
        endif ! nrtyp >= nrtmax
      endif ! ctyp != ctyp0
       
      if ( ctyp(1:6).ne.'ADPUPA' .and.  &
          ctyp(1:6).ne.'AIRCAR' .and.  &
          ctyp(1:6).ne.'PROFLR' .and.  &
          ctyp(1:6).ne.'AIRCFT' ) then
         nr=nr+msub
         cycle
      endif 

      do while (ireadsb(lunin) .eq. 0) 
      nr=nr+1
      nrtyp=nrtyp+1
      nread=nread+1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c===
! prepfits processing
!  based on /meso/save/wx20ps/ucl4/prepfits.tab of 15 Oct 10:45

! dtyps == 
!   ADPUPA AIRCAR AIRCFT SATWND PROFLR ADPSFC SFCSHP VADWND 
!   SATBOG SATEMP SFCBOG SPSSMI SYNDAT ERS1DA GOESND MSONET GPSIPW 

! dtyp => HEADR {PLEVL}
! HEADR == SID XOB YOB DHR ELV TYP T29 ITP
! PLEVL(255) == CAT PRC PQM QQM TQM ZQM WQM  [OBLVL]
! OBLVL(255) == SRC FHR  <PEVN> <QEVN> <TEVN> <ZEVN> <WEVN> <CEVN> <SEVN>
!       PEVN == POB  PMO                                     ! pressure seq
!       QEVN == QOB                                          ! sp hum seq
!       TEVN == TOB                                          ! temperature seq
!       ZEVN == ZOB                                          ! hgt seq
!       WEVN == UOB VOB                                      ! wind seq
!       CEVN == CAPE CINH LI PBL TROP PWO                    ! convective ndx seq
!       SEVN == HOVI MXTM MITM TDO TOCC MXGS THI TCH CDBZ    ! sensible WX seq


! HEADR == SID XOB YOB DHR ELV TYP T29 ITP
!    SID: Station identification               (ascii,     character*8)     a8
!    XOB: Longitude (coarse accuracy)          (deg E,-180.00 -> 475.35)  f7.2
!    YOB: Latitude                             (deg N, -90.00 -> 237.67)  f6.2
!    DHR: Observation time minus cycle time    (hours, -24.00 -> 57.91)   f6.2
!    ELV: Station elevation                        (m,  -1000 -> 130072)    i6
!    TYP: Report type                       (code tbl,      0 -> 511)       i3
!    T29: Data dump report type             (code tbl,      0 -> 1023)      i4
!    ITP: Instrument type                   (code tbl,      0 -> 255)       i3

! PLEVL(255) == CAT PRC PQM QQM TQM ZQM WQM  [OBLVL]
!    CAT: Data level category               (code tbl,      0 -> 63)        i2
!    PRC: Pressure coordinate                     (mb,     .0 -> 1638.3)  f6.1
!    PQM: Pressure (quality) marker         (code tbl,      0 -> 31)        i2
!    QQM: Sp humidity (quality) marker      (code tbl,      0 -> 31)        i2
!    TQM: Temp (TOB) (quality) marker       (code tbl,      0 -> 31)        i2
!    ZQM: Height (quality) marker           (code tbl,      0 -> 31)        i2
!    WQM: u-, v- wind (quality) marker      (code tbl,      0 -> 31)        i2

! OBLVL(255) == SRC FHR  <PEVN> <QEVN> <TEVN> <ZEVN> <WEVN> <CEVN> <SEVN>
!    SRC: File name of data source             (ascii,    character*24)    a24
!    FHR: Forecast length                       (hour,     .0 -> 409.5)   f5.1

! PEVN == POB PMO
!    POB: Pressure observation                    (mb,     .0 -> 1638.3)  f6.1
!    PMO: Mean sea-level pressure observation     (mb,     .0 -> 1638.3)  f6.1

! QEVN == QOB
!    QOB: Specific humidity observation        (mg/kg,      0 -> 65535)     i5

! TEVN == TOB
!    TOB: Temp observation (drybulb or virt)   (deg c, -273.2 -> 1365.1)  f6.1

! ZEVN == ZOB
!    ZOB: Height observation                       (m,  -1000 -> 130071)    i6

! WEVN == UOB VOB
!    UOB: Wind u- component observation          (m/s,  -409.6 -> 409.5)  f6.1
!    VOB: Wind v- component observation          (m/s,  -409.6 -> 409.5)  f6.1

! CEVN == CAPE CINH LI PBL TROP PWO
!   CAPE: convective available pot energy       (j/kg,      0 -> 16383)     i5
!   CINH: convective inhibition                 (j/kg,  -3890 -> 205)       i5
!     LI: Lifted Index                             (K,   -5.0 -> 97.3)    f4.1
!c PLFTI: Lifted Index                             (K,    -20 -> 43)        i2  (v2)
!c  LFTI: PREPBUFR: Lifted Index                   (K,   -5.0 -> 97.3)    f4.1  (v2)
!    PBL: Height of planetary boundary layer       (m,   -4.0 -> 6549.5)  f6.1
!c  HPBL: Height of planetary boundary layer       (m,      0 -> 8191)      i4  (v2)
!   TROP: tropopause level                        (mb,     .0 -> 1638.3)  f5.1
!c  PRLC: Pressure                                (Pa,      0 -> 163830)    i6  (v2)
!    PWO: Tot column precipitable water (kg/m^2 or mm,     .0 -> 204.7)   f5.1

!  SEVN == HOVI MXTM MITM TDO TOCC MXGS THI TCH CDBZ
!   HOVI: Horizontal Visibility                    (m,     .0 -> 819.1)   f5.1
!   MXTM: Maximum temperature                      (K,    .00 -> 655.35)  f6.2
!   MITM: Minimum temperature                      (K,    .00 -> 655.35)  f6.2
!    TDO: Dew-point temperature observation    (deg c, -273.1 -> 1365.2)  f6.1
!   TOCC: cloud cover (total)                      (%,      0 -> 127)       i3
!   MXGS: Maximum wind gust speed                (m/s,     .0 -> 409.5)   f5.1
!    THI: Heat index                           (deg F,     .0 -> 6553.5)  f5.1
!    TCH: wind chill                           (deg F,     .0 -> 6553.5)  f5.1
!   CDBZ: ?                                        (m,     .0 -> 1638.3)  f5.1
!c  HOCB: Height of base of cloud                  (m,   -400 -> 20070)     i5  (v2)
!c===

! HEADR == SID XOB YOB DHR ELV TYP T29 ITP
      cnem='SID XOB YOB DHR ELV TYP T29 ITP'
      call ufbint(lunin,hdr,MXNM,1,iret,cnem)
      if ( iret.ne.1 ) write(*,*)'ERROR - ufbseq(HEADR) iret=',iret
      kx=hdr(6)
      if(kx == 431 .or. kx == 531) nkx= 131
      if(kx == 433 .or. kx == 533) nkx= 133
      if(kx == 435 .or. kx == 535) nkx= 135
      if(kx == 120) nkx= 120
      if(kx == 227) nkx= 181

!     Extract station id, type, date, and location information
      rval(1)=hdr(1)      ! SID

      if(hdr(2)>= r360)hdr(2)=hdr(2)-r360
      if(hdr(2) < zero)hdr(2)=hdr(2)+r360
      dlon_earth_deg=hdr(2)
      dlat_earth_deg=hdr(3)
      dlon_earth=hdr(2)*deg2rad
      dlat_earth=hdr(3)*deg2rad
      if(regional)then
         call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)    ! convert to rotated coordinate
         if(diagnostic_reg) then
            call txy2ll(dlon,dlat,rlon00,rlat00)
            ntest=ntest+1
            cdist=sin(dlat_earth)*sin(rlat00)+cos(dlat_earth)*cos(rlat00)* &
                 (sin(dlon_earth)*sin(rlon00)+cos(dlon_earth)*cos(rlon00))
            cdist=max(-one,min(cdist,one))
            disterr=acos(cdist)*rad2deg
            disterrmax=max(disterrmax,disterr)
         end if
         if(outside) cycle   ! check to see if outside regional domain
      else
         dlat = dlat_earth
         dlon = dlon_earth
         call grdcrd1(dlat,rlats,nlat,1)
         call grdcrd1(dlon,rlons,nlon,1)
      endif

      if(offtime_data) then

!       in time correction for observations to account for analysis
!                time being different from obs file time.
        write(date,'( i10)') idate
        read (date,'(i4,3i2)') iy,im,idd,ihh
        idate5(1)=iy
        idate5(2)=im
        idate5(3)=idd
        idate5(4)=ihh
        idate5(5)=0
        call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
        idate5(1)=iadate(1)
        idate5(2)=iadate(2)
        idate5(3)=iadate(3)
        idate5(4)=iadate(4)
        idate5(5)=0
        call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date

!       Add obs reference time, then subtract analysis time to get obs time relative to analysis

        time_correction=float(minobs-minan)/60._r_kind

      else
        time_correction=zero
      end if

      timeobs=real(real(hdr(4),4),8)
      time=timeobs + time_correction
      t4dv=timeobs + toff
      zeps=1.0e-8_r_kind
      if (t4dv<zero  .and.t4dv>      -zeps) t4dv=zero
      if (t4dv>winlen.and.t4dv<winlen+zeps) t4dv=winlen
      t4dv=t4dv + time_correction
      nc=ikx
      if (l4dvar.or.l4densvar) then
           if (t4dv<zero.OR.t4dv>winlen) cycle
      else
           if((real(abs(time)) > real(ctwind(nc)) .or. real(abs(time)) > real(twindin))) cycle 
      end if

      stnelev=hdr(5)
      pblrielev=stnelev


! PLEVL == CAT PRC PQM QQM TQM ZQM WQM [OBLVL]
!     cnem='CAT PRC PQM QQM TQM ZQM WQM'
!     call ufbint(lunin,plv,MXNM,MXRP, nlev,cnem)

! OBLVL == SRC FHR <PEVN> <QEVN> <TEVN> <ZEVN> <WEVN> <CEVN> <SEVN>
!       == SRC FHR POB PMO QOB    TOB    ZOB   UOB VOB ...  ! {P,Q,T,Z,W}EVN
!     cnem='SRC FHR POB PMO QOB TOB ZOB UOB VOB'
!     call ufbin3(lunin,olv,MXNM,MXRP,MXRP, nlevp,nlevo,cnem)

! CEVN == CAPE CINH LI PBL TROP PWO
!     cnem='CAPE CINH LI PBL TROP PWO'
      cnem='PBL'   ! for Caterina's files (ruc_raobs)
      clv=bmiss
      call ufbin3(lunin,clv,MXNM,MXRP,MXRP, nlevp,nlevc,cnem)
!     pblriob=clv(4,1,2)
      pblriob=clv(1,1,2)
      pblbak=clv(1,1,1)   ! model PBL; from Caterina's files
      if (abs(pblbak-bmiss).lt.10.e5 .or. abs(pblriob-bmiss).lt.10.e5) cycle ! <skip processing of this report>
      
      pblriqm=0
      if (pblriob .lt. 0.0) pblriqm=15
!     if (nkx==131 .or. nkx==133 .or. nkx==135) then
!       anal_time=0
!       obs_time=0
!       tmp_time=zero
!       tmp_time(2)=timeobs
!       anal_time(1)=iadate(1)
!       anal_time(2)=iadate(2)
!       anal_time(3)=iadate(3)
!       anal_time(5)=iadate(4)
!       call w3movdat(tmp_time,anal_time,obs_time) ! observation time

!       lobs_time=0
!       tmp_time=zero
!       cenlon_tmp=hdr(2)
!       if (hdr(2) > 180.0) cenlon_tmp=hdr(2)-360.0_r_kind
!       tmp_time(2)=cenlon_tmp/15.0_r_kind
!       call w3movdat(tmp_time,obs_time,lobs_time) ! local observation time
!       ltime = lobs_time(5)+lobs_time(6)/60.0_r_kind+lobs_time(7)/3600.0_r_kind
!       if ((ltime.gt.21.0) .and. (ltime.lt.5.0)) pblriqm=3
!     end if

      if (kx==531 .or. kx==533 .or. kx==535) pblriqm=3 ! descending profile

!     Set usage variable
      usage = 0.
      if(icuse(nc) <= 0) usage=150.
      if(pblriqm == 15 .or. pblriqm == 9) usage=150.

!     Set inflate_error logical 
      inflate_error=.false.
      if (pblriqm == 3) inflate_error=.true.

! SEVN == HOVI MXTM MITM TDO TOCC MXGS THI TCH CDBZ
!     cnem='HOVI MXTM MITM TDO TOCC MXGS THI TCH CDBZ'
!     call ufbin3(lunin,slv,MXNM,MXRP,MXRP, nlevp,nlevs,cnem)

!--Outputs

! ---HEADR
      if (nrtyp .eq. 1)  write(*,'(/,11x,a,$)') &
       'SID       XOB    YOB     DHR    ELV  TYP  T29 ITP'

      write(*,'(/,i3,i4,$)') nmsg, nr                            ! msg#, rpt#

      rval(1)=hdr(1) ; write(*,'(1x,a,$)') "'"//cval(1)//"'"     ! SID
      write(*,'(1x,f7.2,2(1x,f6.2),$)') (hdr(i),i=2,4)           ! XOB,YOB,DHR
      write(*,'(1x,i6,1x,i3,1x,i4,1x,i3,$)') (int(hdr(i)),i=5,8) ! ELV,TYP,T29,ITP

      write(*,'(1x,2(2x,a,1x,i3),a,$)') '(olv=',nlevo,'nlevp=',nlevp,')'
      write(*,'(1x,2(2x,a,1x,i3),a,$)') '(clv=',nlevc,'nlevp=',nlevp,')'

      ndata=ndata+1
      nodata=nodata+1
      iout=ndata

      if(ndata > maxobs) then
           write(6,*)'READ_PREPFITS:  ***WARNING*** ndata > maxobs for ',obstype
           ndata = maxobs
      end if

!     Get information from surface file necessary for conventional data here
      call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)

!     setup for sort of averaged obs 
      pblrioe=400.0_r_kind  ! temporarily
      if (nkx==120) pblrioe=200.
      if (nkx==181) pblrioe=300.
      if (inflate_error) pblrioe=pblrioe*1.5_r_kind

      cdata_all(1,iout)=pblrioe                  ! pblri error (cb)
      cdata_all(2,iout)=dlon                    ! grid relative longitude
      cdata_all(3,iout)=dlat                    ! grid relative latitude
      cdata_all(4,iout)=pblrielev                ! pblri obs elevation
      cdata_all(5,iout)=pblriob                  ! pblri obs
      cdata_all(6,iout)=rval(1)                 ! station id
      cdata_all(7,iout)=t4dv                    ! time
      cdata_all(8,iout)=nc                      ! type
      cdata_all(9,iout)=pblrioe*three            ! max error
      cdata_all(10,iout)=pblriqm                 ! quality mark
      cdata_all(11,iout)=usage                  ! usage parameter
      cdata_all(12,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
      cdata_all(13,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
      cdata_all(14,iout)=stnelev                ! station elevation (m)

      end do ! while ireadsb

      end do ! while ireadmg
      write(*,*) ! closing linefeed, debug?

      call closbf(lunin)
!   Normal exit

!   Write observation to scratch file
     call count_obs(ndata,nreal,ilat,ilon,cdata_all,nobs)
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((cdata_all(j,i),j=1,nreal),i=1,ndata)
     deallocate(cdata_all)
 
     if (ndata == 0) then
        call closbf(lunin)
        write(6,*)'READ_PREPFITS:  closbf(',lunin,')'
     endif

     close(lunin)
     end subroutine read_pblri_prepfits

     subroutine read_pblri_calipso_cats_mplnet_raob_wp(nread,ndata,nodata,infile,obstype,lunout,twindin,&
         sis,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_pblri_calipso_cats_raob_wp     read obs from msgs in NETCDF files
!
! program history log:
!   2022-10     Y. Zhu - adapted from read_pblh
!   2023-04     E.-G. Yang - added thinning based on distance for calipso and cats
!   2023-07     E.-G. Yang - added MPLNET and calculated hourly mean MPLNET pblh data
!   2023-08     E.-G. Yang - subtracted terrain height from CALIPSO PBLH (above mean sea level)
!                            in order to match a pblri model variable (above sfc height).
!                            1) CALIPSO: terrain height from obs file
!   2023-08     E.-G. Yang - added qc for CATS to exclude night-time data over land
!   2023-08     E.-G. Yang - inflate MPLNET obs error depending on the number of obs in DA window 
!   2023-10     E.-G. Yang - moving average for calipso and cats (within 12 km radius, 25 km thinning)
!                            Before moving average for CALIPSO, sfc height is subtracted.
!   2023-10     E.-G. Yang - added Radar Wind Profiler
!   2023-10     E.-G. Yang - changing the way of converting UTC to local time for obs
!                            Except for hourly obs, local time is obtained by adding float(lon/15) to UTC (YYYYMMDD HH:MI)
!                            For hourly obs (MPLNET, Wind Profiler), local time: int(lon/15) + UTC (YYYYMMDD HH:MI)
!   2023-10     E.-G. Yang - removed qc regarding LT for CALIPSO
!   2023-10     E.-G. Yang - added qc for WP
!   2023-11     E.-G. Yang - added mixed sfc type qc for all obs 
!   2023-11     E.-G. Yang - added high station elevation qc for raob
!   2023-11     E.-G. Yang - commented out mixed sfc type qc for raob
!   2023-11     E.-G. Yang - read in q-gradient based pblh
!   2023-12     E.-G. Yang - added qc for too low pblri when pblri_obs < 50 m and pblri_obs-pblh_qgrad > 1km
!   2023-12     E.-G. Yang - modified distance w.r.t. analysis grid (50 km) for thinning and moving average of calipso and cats
!   2023-12     E.-G. Yang - added night-time qc for MPLNET and WP (6 PM - 10 AM)
!   2024-01     E.-G. Yang - removed qc for too low pblri when pblri_obs < 50 m and pblri_obs-pblh_qgrad > 1km
!   2024-02     Y. Zhu - modelling observation error

      use netcdf
      use buddycheck_mod, only: gc_dist
      implicit none

!     Declare passed variables
      character(len=*),intent(in):: infile,obstype
      character(20),intent(in):: sis
      integer(i_kind),intent(in):: lunout
      integer(i_kind),intent(inout):: nread,ndata,nodata
      integer(i_kind),dimension(npe),intent(inout):: nobs
      real(r_kind),intent(in):: twindin
!     real(r_kind),dimension(nlat,nlon,nsig),intent(in):: hgtl_full

!     Declare local parameters
      real(r_kind),parameter:: r360 = 360.0_r_kind

      real(r_kind),allocatable,dimension(:,:):: cdata_all

      character(10) date
      logical first,outside,inflate_error,lexist,more_data
      logical thin_calipso, thin_cats, thin_mplnet
      integer(i_kind) iret,im,iy,ihh,idd,i,j,k
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      integer(i_kind) kk,klon1,klat1,klonp1,klatp1
      integer(i_kind) ntest,nchanl
      integer(i_kind) pblriqm,maxobs,idomsfc
      integer(i_kind) oltime
      integer(i_kind) hr,ltime_mm,ltime_dd,ltime_hh,ltime_min
      real(r_kind) hr_min
!     integer(i_kind),dimension(8):: obs_time,anal_time,lobs_time
!     real(r_kind) ltime,cenlon_tmp
      real(r_kind) usage
      real(r_kind) pblriob,pblrioe,pblrielev,pblriosfc
      real(r_kind) dlat,dlon,dlat_earth,dlon_earth
      real(r_kind) dlat_earth_deg,dlon_earth_deg
      real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
      real(r_kind) :: tsavg,ff10,sfcr,zz
!     real(r_kind),dimension(5):: tmp_time
!     real(r_kind) sin2,termg,termr,termrg
!     real(r_kind),dimension(nsig):: zges,hges
      real(r_kind) dx,dy,dx1,dy1,w00,w10,w01,w11

      integer(i_kind) idate5(5),minobs,minan
      real(r_kind) time_correction,timeobs,time,toff,t4dv,zeps

!     real(r_single) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblri_calipso
      real(r_kind) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblri_obs_kind

      integer(i_kind) :: ncid,ierr,dimid1,dimid2,dimid3,norbits,nheights
      integer(i_kind) :: noobs, strlen, site_maxstrlen, instrument_maxstrlen
      integer(i_kind) :: varid0,varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11,varid12,varid13
      integer(i_kind) :: varid14, varid15
      integer(i_kind) :: ioyr, iomo, iody, iohr, iomi, iodate
      integer(i_kind), dimension(1) :: ana_time
      real(r_kind), allocatable, dimension(:) :: lat, lon, pblri
      real(r_kind), allocatable, dimension(:) :: ryear, rmonth, rdate, rhour, rminute
      ! calipso
      real(r_kind), allocatable, dimension(:) :: sfc_elev, sfc_mask_calipso, liadr_data_alt
      real(r_kind), allocatable, dimension(:,:) :: ATB
      ! cats
      real(r_kind), allocatable, dimension(:) :: loc_hour, qc_flag, day_night_flag, sfc_mask_cats
      ! mplnet
      real(r_double), allocatable, dimension(:) :: rsecond
      real(r_kind), allocatable, dimension(:) :: surface_alti, flag_cloud_screen
      real(r_kind), allocatable, dimension(:) :: pblri_no_hr_in_da_window
      character(19),allocatable,dimension(:) :: site
      character(19),allocatable,dimension(:,:) :: site_org
      character(8),allocatable,dimension(:) :: instrument
      character(8),allocatable,dimension(:,:) :: instrument_org
      character(8)   :: c_instrument
      real(r_double) :: r_instrument
      integer(i_kind) :: use_cloud ! 0: consider only cloud free 1: consider all data
      ! mplnet tempoal thinning: hourly mean
      real(r_kind), allocatable, dimension(:) :: lat_hr, lon_hr, pblri_hr
      real(r_kind), allocatable, dimension(:) :: ryear_hr, rmonth_hr, rdate_hr, rhour_hr, rminute_hr
      real(r_kind), allocatable, dimension(:) :: surface_alti_hr, flag_use_cloud_hr, flag_use_cloud
      character(19),allocatable,dimension(:) :: site_hr
      character(8),allocatable,dimension(:) :: instrument_hr
      ! radiosonde
      real(r_kind), allocatable, dimension(:) :: nplevels, nplevels500, stn_elev,pblh_q_grad
      integer(i_kind), allocatable, dimension(:) :: loc_time, qc_mark_pblh
      character(8),allocatable,dimension(:) :: sid
      character(8),allocatable,dimension(:,:) :: sid_org
      character(8) :: c_station_id
      real(r_double) :: r_station_id
      integer(i_kind) :: ios
      ! wind profiler
      real(r_kind),allocatable,dimension(:) :: cloud_fraction, error, negative_snr_indicator
      character(5),allocatable,dimension(:) :: sid_wp
      character(5),allocatable,dimension(:,:) :: sid_wp_org
      character(5) :: c_station_id_wp
      real(r_double) :: r_station_id_wp
      ! Thinning
      real(r_kind) :: o1_lat, o1_lon, o2_lat, o2_lon, dist, thindist
      integer(i_kind) :: nthinobs
      ! Thinning for CALIPSO
      real(r_kind), allocatable, dimension(:) :: lat_thin_final, lon_thin_final, pblri_mvavg_thin_final
      ! Thinning for cats
      real(r_kind) :: mvavg_rad
      real(r_kind), allocatable, dimension(:) :: lat_qc_thin_final, lon_qc_thin_final, pblri_qc_mvavg_thin_final
      ! Observation error model
      real(r_kind) :: cpblh1, cpblh2, cerr1, cerr2

      logical :: selected
!     equivalence to handle character names
      equivalence(r_station_id,c_station_id)
      equivalence(r_instrument,c_instrument)
      equivalence(r_station_id_wp,c_station_id_wp)

!     Initialize obs err parameters
      cpblh1 = 2000.0_r_kind
      cpblh2 = 4000.0_r_kind
      cerr1 = 200.0_r_kind
      cerr2 = 500.0_r_kind

!     Check if pblri file exists
      inquire(file=trim(infile),exist=lexist)
      if (.not.lexist) return

!     Read data
      ierr =  NF90_OPEN(trim(infile),0,ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"open")

!     Read dimensions depending on observation types
      if (index(trim(infile),'calipso') > 0) then

         kx = 890 ! calipso
         nkx = 890 ! calipso
         ierr = NF90_INQ_DIMID(ncid,'norbits',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"norbits")
         ierr = NF90_INQ_DIMID(ncid,'nheights',dimid2)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nheights")

         ierr = nf90_inquire_dimension(ncid, dimid1, len = norbits)
         ierr = nf90_inquire_dimension(ncid, dimid2, len = nheights)
         print*, 'read_pblri,',trim(infile),': norbits=', norbits, ' nheights=', nheights

      else if (index(trim(infile),'cats') > 0) then

         kx = 891 ! cats
         nkx = 891 ! cats
         ierr = NF90_INQ_DIMID(ncid,'nobs',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nobs")

         ierr = nf90_inquire_dimension(ncid, dimid1, len = noobs)
         print*, 'read_pblri,',infile,': noobs=', noobs

      else if (index(trim(infile),'mplnet') > 0) then

         kx = 892 ! mplnet
         nkx = 892 ! mplnet
         ierr = NF90_INQ_DIMID(ncid,'nobs',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nobs")
         ierr = nf90_inquire_dimension(ncid, dimid1, len = noobs)

         ierr = NF90_INQ_DIMID(ncid,'Site_maxstrlen',dimid2)
         if (ierr /= nf90_noerr) call handle_err(ierr,"Site_maxstrlen")
         ierr = nf90_inquire_dimension(ncid, dimid2, len = site_maxstrlen)

         ierr = NF90_INQ_DIMID(ncid,'Instrument_maxstrlen',dimid3)
         if (ierr /= nf90_noerr) call handle_err(ierr,"Instrument_maxstrlen")
         ierr = nf90_inquire_dimension(ncid, dimid3, len = instrument_maxstrlen)

         print*, 'read_pblri,',infile,': noobs=', noobs,' site_maxstrlen=',site_maxstrlen, 'instrument_maxstrlen=',instrument_maxstrlen


      else if (index(trim(infile),'raob') > 0) then

         kx = 120 ! radiosonde
         nkx = 120 ! radiosonde

         ierr = NF90_INQ_DIMID(ncid,'nsid',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nsid")
         ierr = nf90_inquire_dimension(ncid, dimid1, len = noobs)

         ierr = NF90_INQ_DIMID(ncid,'Station_ID_maxstrlen',dimid2)
         if (ierr /= nf90_noerr) call handle_err(ierr,"Station_ID_maxstrlen")
         ierr = nf90_inquire_dimension(ncid, dimid2, len = strlen)

         print*, 'read_pblri,',trim(infile),': noobs=', noobs

      else if (index(trim(infile),'wp') > 0) then ! Wind profiler

         kx = 893 ! wind profiler 
         nkx = 893 ! wind profiler

         ierr = NF90_INQ_DIMID(ncid,'nobs',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nobs")
         ierr = nf90_inquire_dimension(ncid, dimid1, len = noobs)

         ierr = NF90_INQ_DIMID(ncid,'Sid_maxstrlen',dimid2)
         if (ierr /= nf90_noerr) call handle_err(ierr,"Sid_maxstrlen")
         ierr = nf90_inquire_dimension(ncid, dimid2, len = strlen)

         print*, 'read_pblri,',trim(infile),': noobs=', noobs


      end if

      select case(kx) ! Temporary setting
         case (890) ! calipso
           maxobs=norbits 
         case (891) ! cats
            maxobs=noobs
         case (892) ! mplnet
            maxobs=noobs
         case (120) ! radiosonde
            maxobs=noobs
         case (893) ! wp
            maxobs=noobs
      end select

!     Allocate
      allocate(lat(maxobs), lon(maxobs), pblri(maxobs))
      allocate(ryear(maxobs), rmonth(maxobs), rdate(maxobs), rhour(maxobs), rminute(maxobs))
      select case(kx)
         case (890) ! calipso
            allocate(sfc_elev(norbits), sfc_mask_calipso(norbits), liadr_data_alt(nheights), ATB(norbits,nheights))
         case (891) ! cats
            allocate(loc_hour(noobs), qc_flag(noobs), day_night_flag(noobs), sfc_mask_cats(noobs))
         case (892) ! mplnet
            allocate(rsecond(noobs), surface_alti(noobs), flag_cloud_screen(noobs), site_org(noobs,site_maxstrlen), site(noobs))
            allocate(instrument_org(noobs,instrument_maxstrlen), instrument(noobs))
         case (120) ! radiosonde
            !allocate(sid_org(noobs, strlen), sid(noobs), loc_time(noobs), qc_mark_pblh(noobs), nplevels(noobs), nplevels500(noobs), stn_elev(noobs))
            allocate(sid_org(noobs,strlen), sid(noobs), loc_time(noobs), qc_flag(noobs), nplevels(noobs), nplevels500(noobs), stn_elev(noobs))
            allocate(pblh_q_grad(noobs))
         case (893) ! wind profiler
            allocate(rsecond(noobs), surface_alti(noobs), cloud_fraction(noobs), error(noobs), negative_snr_indicator(noobs))
            allocate(sid_wp_org(noobs,strlen), sid_wp(noobs))
      end select

!     Read variables

!     Analysis Time
      ierr = NF90_INQ_VARID(ncid,'Ana_Time',varid0)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid0,ana_time)
      print*, "YEGREAD kx=",kx,",ana_time=", ana_time
!     Latitude: degrees
      ierr = NF90_INQ_VARID(ncid,'lat',varid1)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid1,lat)
!     Longitude: degrees
      ierr = NF90_INQ_VARID(ncid,'lon',varid2)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid2,lon)
      ierr = NF90_INQ_VARID(ncid,'Year',varid3)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid3,ryear)
      ierr = NF90_INQ_VARID(ncid,'Month',varid4)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid4,rmonth)
      ierr = NF90_INQ_VARID(ncid,'Date',varid5)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid5,rdate)
      ierr = NF90_INQ_VARID(ncid,'Hour',varid6)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid6,rhour)
      ierr = NF90_INQ_VARID(ncid,'Minute',varid7)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid7,rminute)
!     PBL_Height: meters
      select case(kx)
         case (890) ! calipso
            ierr = NF90_INQ_VARID(ncid,'PBL_Height',varid8)
         case (891) ! cats
            ierr = NF90_INQ_VARID(ncid,'PBL_Height',varid8)
         case (892) ! mplnet
            ierr = NF90_INQ_VARID(ncid,'PBL_Height',varid8)
         case (120) ! radiosonde
            ierr = NF90_INQ_VARID(ncid,'PBL_Height_ri',varid8)
         case (893) ! wind profiler
            ierr = NF90_INQ_VARID(ncid,'PBL_Height',varid8)
      end select

      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,pblri)

!     Other variables
      select case(kx)

         case (890) ! calipso

            !SurfaceElevation: meters
            ierr = NF90_INQ_VARID(ncid,'SurfaceElevation',varid9)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,sfc_elev)
            !0=shallow ocean, 1=land, 2=coastlines, 3=shallow inland water, 4=intermittent water, 
            !5=deep inland water, 6=continental ocean, 7=deep ocean
            ierr = NF90_INQ_VARID(ncid,'Land_Water_Mask',varid10)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,sfc_mask_calipso)
            !Lidar_Data_Altitude: km
            ierr = NF90_INQ_VARID(ncid,'Lidar_Data_Altitude',varid11)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,liadr_data_alt)
            !Total_Attenuated_Backscatter_532: km^-1 sr^-1
            ierr = NF90_INQ_VARID(ncid,'Total_Attenuated_Backscatter_532',varid12)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,ATB)

         case (891) ! cats

            !Local_Hour: hour
            ierr = NF90_INQ_VARID(ncid,'Local_Hour',varid9)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,loc_hour)
            ierr = NF90_INQ_VARID(ncid,'QC_Flag',varid10)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,qc_flag)
            ierr = NF90_INQ_VARID(ncid,'Day_Night_Flag',varid11)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,day_night_flag)
            ierr = NF90_INQ_VARID(ncid,'Land_Water_Mask',varid12)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,sfc_mask_cats)

         case (892) ! mplnet

            ierr = NF90_INQ_VARID(ncid,'Second',varid9)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,rsecond)
            ierr = NF90_INQ_VARID(ncid,'Surface_Altitude',varid10)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,surface_alti)
            ierr = NF90_INQ_VARID(ncid,'Flag_Cloud_Screen',varid11)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,flag_cloud_screen)

            ierr = NF90_INQ_VARID(ncid,'Site',varid12)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,site_org)
            ierr = NF90_INQ_VARID(ncid,'Instrument',varid13)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid13,instrument_org)

            do i = 1, noobs
               site(i)=site_org(i,1)
               instrument(i)=instrument_org(i,1)
               print*, "YEG_mplnet_CHECK: i=,",i,",site(i)=",site(i),",instruemnt(i)=",instrument(i)
            end do


         case (120) ! radiosonde

            ierr = NF90_INQ_VARID(ncid,'Sid',varid9)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,sid_org)

            do i = 1, noobs
               sid(i)=sid_org(i,1)
               print*, "YEG_RAOB_CHECK: i=,",i,",sid(i)=",sid(i)
            end do

         
            ierr = NF90_INQ_VARID(ncid,'Loc_Time_ri',varid10) ! YYYYMMDDHH
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,loc_time)
            ierr = NF90_INQ_VARID(ncid,'QC_Mark_PBLH',varid11)
            !if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,qc_mark_pblh)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,qc_flag)
            ierr = NF90_INQ_VARID(ncid,'Number_of_PLevels',varid12)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,nplevels)
            ierr = NF90_INQ_VARID(ncid,'Number_of_PLevels_500',varid13)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid13,nplevels500)
            ierr = NF90_INQ_VARID(ncid,'Stn_Elev',varid14)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid14,stn_elev)
            ierr = NF90_INQ_VARID(ncid,'PBL_Height_q',varid15)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid15,pblh_q_grad)

         case (893) ! wind profiler

            ierr = NF90_INQ_VARID(ncid,'Second',varid9)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,rsecond)
            ierr = NF90_INQ_VARID(ncid,'Altitude',varid10)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,surface_alti)
            ierr = NF90_INQ_VARID(ncid,'Cloud_Fraction',varid11)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,cloud_fraction)
            ierr = NF90_INQ_VARID(ncid,'Vert_Resolution_2x',varid12)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,error)
            ierr = NF90_INQ_VARID(ncid,'Negative_SNR_Indicator',varid13)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid13,negative_snr_indicator)

            ierr = NF90_INQ_VARID(ncid,'Station_ID',varid14)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid14,sid_wp_org)

            do i = 1, noobs
               sid_wp(i)=sid_wp_org(i,1)
               print*, "YEG_wp_CHECK: i=,",i,",sid_wp_org(i,1)=",sid_wp_org(i,1)
               print*, "YEG_wp_CHECK: i=,",i,",sid_wp(i)=",sid_wp(i)
            end do

      end select

      ierr = NF90_CLOSE(ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"close")

      nreal=17 ! 16 ! 15 ! Temporary setting

      !-------------------------------
      ! Spatial thinning for CALIPSO and CATS
      ! Temporal thinning for MPLNET
      !-------------------------------
      thin_calipso=.true.
      thin_cats=.true.
      thin_mplnet=.true.
      !-------------------------------
      ! COUNT for thinning: CALIPSO (subtracting sfc_elev and moving average)
      !-------------------------------
      if (thin_calipso .and. kx==890) then

         !count the number of obs after thinning 
         thindist=50.0_r_kind ! [km] GSI Grid, analysis grid (lower resolution than model grid)
         mvavg_rad=25.0_r_kind !12! [km] moving average radius 25 km
         call count_elev_mvavg_thin_obs(maxobs,lat,lon,sfc_elev,pblri,thindist,mvavg_rad,nthinobs,lat_thin_final,lon_thin_final,pblri_mvavg_thin_final)

         allocate(cdata_all(nreal,nthinobs))
         print*, "yeg_thinning_calipso: maxobs=",maxobs,",nthinobs=",nthinobs

      !-------------------------------
      ! COUNT for thinning: CATS (considering qc and moving average)
      !-------------------------------
      elseif (thin_cats .and. kx==891) then

         !count the number of obs after thinning 
         thindist=50.0_r_kind ! [km] GSI Grid
         mvavg_rad=25.0_r_kind ! [km] moving average radius 25 km
         !Thinning based on distance as well as qc
         call count_qc_mvavg_thin_obs(maxobs,lat,lon,qc_flag,pblri,thindist,mvavg_rad,nthinobs,lat_qc_thin_final,lon_qc_thin_final,pblri_qc_mvavg_thin_final)

         allocate(cdata_all(nreal,nthinobs))
         print*, "yeg_thinning_cats: maxobs=",maxobs,",nthinobs=",nthinobs

      !-------------------------------
      ! thinning for MPLNET
      !-------------------------------
      ! do thin here although other obs do not here.
      elseif (thin_mplnet .and. kx==892) then

         print*, "yeg_thinning_mplnet: before thinning: maxobs=",maxobs

         !calculate hourly mean
         use_cloud = 1 ! 0: consider only cloud free 1: consider all data (including cloud affected data)

         call hourly_mean_obs(maxobs,lat,lon,ryear,rmonth,rdate,rhour,rminute,pblri,& 
                              surface_alti,flag_cloud_screen,site,instrument,use_cloud,&
                              nthinobs,lat_hr,lon_hr,ryear_hr,rmonth_hr,rdate_hr,rhour_hr,rminute_hr,pblri_hr,&
                              surface_alti_hr,flag_use_cloud_hr,site_hr,instrument_hr, pblri_no_hr_in_da_window)

         allocate(cdata_all(nreal,nthinobs))
         maxobs=nthinobs
         print*, "yeg_thinning_mplnet:  after thinning: maxobs=",maxobs,", nthinobs=",nthinobs ! 2160 -> 36

!        Deallocate for thinned data
         deallocate(lat, lon, pblri)
         deallocate(ryear, rmonth, rdate, rhour, rminute)
         deallocate(rsecond, surface_alti,flag_cloud_screen,site,instrument)

!        Allocate for thinned data
         allocate(lat(maxobs), lon(maxobs), pblri(maxobs))
         allocate(ryear(maxobs), rmonth(maxobs), rdate(maxobs), rhour(maxobs), rminute(maxobs))
         allocate(surface_alti(maxobs), flag_use_cloud(maxobs), site(maxobs))
         allocate(instrument(maxobs))
         lat = lat_hr
         lon = lon_hr
         pblri = pblri_hr
         ryear = ryear_hr
         rmonth = rmonth_hr
         rdate = rdate_hr
         rhour = rhour_hr
         rminute = rminute_hr
         surface_alti = surface_alti_hr
         flag_use_cloud = flag_use_cloud_hr
         site = site_hr
         instrument = instrument_hr

      else

         allocate(cdata_all(nreal,maxobs))

      end if
      nchanl=0
      nread=0

      first = .true.
      !maxobs=0
      !do i = 1, norbits
      do i = 1, maxobs

         print*, " -------------------------"
!        Time offset         
!        Obs Time
         ioyr = ryear(i)
         iomo = rmonth(i)
         iody = rdate(i)
         iohr = rhour(i)
         iomi = rminute(i)

         if (ioyr==r_missing) then
            write (date,'(i4,3i2.2)') 0,0,0,0
         else
            write (date,'(i4,3i2.2)') ioyr,iomo,iody,iohr
         end if
         read (date,'( i10)') iodate

         print*, 'i=',i,', kx=',kx,',iodate=', iodate

         if (first) then
            !call time_4dvar(idate,toff) 
            call time_4dvar(ana_time(1),toff) 
            ! toff: time since start of 4D-Var window (hours)
            print*, 'YEG12: iodate,toff=',iodate,toff
            first=.false.
         end if

         nread=nread+1
         !kx=890
         if(kx == 120) nkx= 120
         if(kx == 890) nkx= 890
         if(kx == 891) nkx= 891
         if(kx == 892) nkx= 892
         if(kx == 893) nkx= 893

         !----------------------------------------------
         ! thinning calipso depending on distance
         if (thin_calipso .and. kx==890) then

            selected=.False.
            findloopp: do j = 1, nthinobs

               if (lat(i)==lat_thin_final(j) .and. lon(i)==lon_thin_final(j)) then
                  selected=.True.
                  print*, "Selected calipso: i=",i,",j=",j,",lat(i),lon(i)=",lat(i),lon(i)
                  print*, "Selected calipso: i=",i,",j=",j,",pblri(i)=",pblri(i)
                  pblri(i)=pblri_mvavg_thin_final(j)
                  print*, "Selected calipso: i=",i,",j=",j,",pblri_mvavg_thin_final(i)=",pblri(i)
                  exit findloopp
               end if

            end do findloopp

            if (.not.selected) cycle

         end if

         !----------------------------------------------
         ! thinning cats depending on distance as well as qc
         if (thin_cats .and. kx==891) then

            selected=.False.
            findloop: do j = 1, nthinobs

               if (lat(i)==lat_qc_thin_final(j) .and. lon(i)==lon_qc_thin_final(j)) then
                  selected=.True.
                  print*, "Selected cats: i=",i,",j=",j,",lat(i),lon(i)=",lat(i),lon(i)
                  print*, "Selected cats: i=",i,",j=",j,",pblri(i)=",pblri(i)
                  pblri(i)=pblri_qc_mvavg_thin_final(j)
                  print*, "Selected cats: i=",i,",j=",j,",pblri_qc_mvavg_thin_final(i)=",pblri(i)
                  exit findloop
               end if

            end do findloop

            if (.not.selected) cycle

         end if


!        Is pblri in convinfo file
         ikx=0
         loop_convinfo_pblri: do j=1,nconvtype
           if (trim(ioctype(j)) /= trim(obstype))cycle loop_convinfo_pblri
           if (nkx == ictype(j)) then
              ikx=j
              exit
           end if
         end do loop_convinfo_pblri

         print*, "YEG_read_pblri L819:nkx=",nkx,", ikx=",ikx
         if(ikx == 0) cycle

         lon_deg = lon(i)
         lat_deg = lat(i)
         if(lon_deg>= r360)lon_deg=lon_deg-r360
         if(lon_deg < zero)lon_deg=lon_deg+r360
         dlon_earth_deg=lon_deg
         dlat_earth_deg=lat_deg
         dlon_earth=lon_deg*deg2rad
         dlat_earth=lat_deg*deg2rad
         if(regional)then
            call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)    ! convert to rotated coordinate
            if(diagnostic_reg) then
               call txy2ll(dlon,dlat,rlon00,rlat00)
               ntest=ntest+1
               cdist=sin(dlat_earth)*sin(rlat00)+cos(dlat_earth)*cos(rlat00)* &
                    (sin(dlon_earth)*sin(rlon00)+cos(dlon_earth)*cos(rlon00))
               cdist=max(-one,min(cdist,one))
               disterr=acos(cdist)*rad2deg
               disterrmax=max(disterrmax,disterr)
            end if
            if(outside) cycle   ! check to see if outside regional domain
         else
            dlat = dlat_earth
            dlon = dlon_earth
            call grdcrd1(dlat,rlats,nlat,1)
            call grdcrd1(dlon,rlons,nlon,1)
         endif

!        Interpolate guess pressure profile to observation location
!        klon1= int(dlon);  klat1= int(dlat)
!        dx   = dlon-klon1; dy   = dlat-klat1
!        dx1  = one-dx;     dy1  = one-dy
!        w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy

!        klat1=min(max(1,klat1),nlat); klon1=min(max(0,klon1),nlon)
!        if (klon1==0) klon1=nlon
!        klatp1=min(nlat,klat1+1); klonp1=klon1+1
!        if (klonp1==nlon+1) klonp1=1
!        do kk=1,nsig
!           hges(kk)=w00*hgtl_full(klat1 ,klon1 ,kk) +  &
!                    w10*hgtl_full(klatp1,klon1 ,kk) + &
!                    w01*hgtl_full(klat1 ,klonp1,kk) + &
!                    w11*hgtl_full(klatp1,klonp1,kk)
!        end do
!        sin2  = sin(dlat_earth)*sin(dlat_earth)
!        termg = grav_equator * &
!           ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
!        termr = semi_major_axis /(one + flattening + grav_ratio -  &
!           two*flattening*sin2)
!        termrg = (termg/grav)*termr
!        do k=1,nsig
!           zges(k) = (termr*hges(k)) / (termrg-hges(k))
!        end do

         if(offtime_data) then

!           in time correction for observations to account for analysis
!                time being different from obs file time.
            write(date,'( i10)') iodate
            read (date,'(i4,3i2)') iy,im,idd,ihh
            !idate5(1)=iyear
            !idate5(2)=imonth
            !idate5(3)=idate
            !idate5(4)=ihour
            idate5(1)=iy
            idate5(2)=im
            idate5(3)=idd
            idate5(4)=ihh
            idate5(5)=0
            call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
            idate5(1)=iadate(1)
            idate5(2)=iadate(2)
            idate5(3)=iadate(3)
            idate5(4)=iadate(4)
            idate5(5)=0
            call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date

!           Add obs reference time, then subtract analysis time to get obs time relative to analysis

            time_correction=float(minobs-minan)/60._r_kind

         else
            time_correction=zero
         end if

!-----------------------------------------------------------------------
!        Time check
!-----------------------------------------------------------------------
!        timeobs (obs - analysis time [hrs]) should be calculated here. 

!        in time correction for observations to account for analysis
!              time being different from obs file time.
         ! obs time
         idate5(1)=ioyr
         idate5(2)=iomo
         idate5(3)=iody
         idate5(4)=iohr
         idate5(5)=iomi
         call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
         idate5(1)=iadate(1)
         idate5(2)=iadate(2)
         idate5(3)=iadate(3)
         idate5(4)=iadate(4)
         idate5(5)=0
         call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date

!        Add obs reference time, then subtract analysis time to get obs time relative to analysis
         timeobs = float(minobs-minan)/60._r_kind
!-----------------------------------------------------------------------
!        timeobs=real(rhour(i), r_double)
         time=timeobs + time_correction
         t4dv=timeobs + toff
         zeps=1.0e-8_r_kind
         if (t4dv<zero  .and.t4dv>      -zeps) t4dv=zero
         if (t4dv>winlen.and.t4dv<winlen+zeps) t4dv=winlen
         t4dv=t4dv + time_correction
         if (l4dvar.or.l4densvar) then
           if (t4dv<zero.OR.t4dv>winlen) cycle
         else
           if((real(abs(time)) > real(ctwind(ikx)) .or. real(abs(time)) > real(twindin))) cycle 
         end if
         print*, 'read_pblri: iodate, iadate, timeobs, toff, t4dv=', iodate, iadate, timeobs, toff, t4dv
         print*, 'YEG GRID: t4dv =', t4dv

         !------------------------------------------
         ! Local Time (only available for 2015/Aug)
         !------------------------------------------
         ! Converting UTC to Local Time (UTC + hr)
         ! longitude (degree East) -> East: lon>0, West: lon<0
         ! 1) hourly obs (MPLNET and Radar wind profiler)
         ! 2) non-hourly obs (CALIPSO, CATS, raob)

         select case(kx)

            case(890) ! CALIPSO
               ! Converting UTC to Local Time (UTC + hr_min)
               ! longitude (degree East) -> East: lon>0, West: lon<0
               if (lon(i)>180) then
                  hr_min = (lon(i)-360.0)/15.
               else
                  hr_min = lon(i)/15.
               end if

               print*, "yeg_read_pblri 890 L964 iodate=",iodate,", lon(i)=",lon(i),", hr_min=",hr_min
               call convert_localtime_min(ltime_mm, ltime_dd, ltime_hh, ltime_min, hr_min, iomo, iody, iohr,iomi)
               oltime = ioyr*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh
               print*, "yeg_read_pblri 890 L966 Local time, oltime=",oltime

            case(120) ! raob
!               oltime = loc_time(i) ! YYYYMMDDHH
               ! Converting UTC to Local Time (UTC + hr_min)
               ! longitude (degree East) -> East: lon>0, West: lon<0
               if (lon(i)>180) then
                  hr_min = (lon(i)-360.0)/15.
               else
                  hr_min = lon(i)/15.
               end if

               print*, "yeg_read_pblri 120 iodate=",iodate,", lon(i)=",lon(i),", hr_min=",hr_min
               call convert_localtime_min(ltime_mm, ltime_dd, ltime_hh, ltime_min, hr_min, iomo, iody, iohr,iomi)
               oltime = ioyr*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh
               print*, "yeg_read_pblri 120 Local time, oltime=",oltime

            case(891) ! CATS
               !oltime = r_missing ! using localtime from cats or newly calculate ??????????? 
               ! Converting UTC to Local Time (UTC + hr_min)
               ! longitude (degree East) -> East: lon>0, West: lon<0
               if (lon(i)>180) then
                  hr_min = (lon(i)-360.0)/15.
               else
                  hr_min = lon(i)/15.
               end if

               print*, "yeg_read_pblri 891 iodate=",iodate,", lon(i)=",lon(i),", hr_min=",hr_min
               call convert_localtime_min(ltime_mm, ltime_dd, ltime_hh, ltime_min, hr_min, iomo, iody, iohr,iomi)
               oltime = ioyr*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh
               print*, "yeg_read_pblri 891 Local time, oltime=",oltime

           case(892) ! MPLNET <- int(lon/15)
               ! Converting UTC to Local Time (UTC + hr)
               ! longitude (degree East) -> East: lon>0, West: lon<0
               if (lon(i)>180) then
                  hr = nint((lon(i)-360)/15.)
               else
                  hr = nint(lon(i)/15.)
               end if

               print*, "yeg_read_pblri 892 MPLNET iodate=",iodate,", lon(i)=",lon(i),", hr=",hr
               call convert_localtime_hr(ltime_mm, ltime_dd, ltime_hh, hr, iomo, iody, iohr)
               oltime = ioyr*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh
               print*, "yeg_read_pblri 892 MPLNET Local time, oltime=",oltime

           case(893) ! RWP <- int(lon/15)
               ! Converting UTC to Local Time (UTC + hr)
               ! longitude (degree East) -> East: lon>0, West: lon<0
               if (lon(i)>180) then
                  hr = nint((lon(i)-360)/15.)
               else
                  hr = nint(lon(i)/15.)
               end if

               print*, "yeg_read_pblri 893 WP iodate=",iodate,", lon(i)=",lon(i),", hr=",hr
               call convert_localtime_hr(ltime_mm, ltime_dd, ltime_hh, hr, iomo, iody, iohr)
               oltime = ioyr*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh
               print*, "yeg_read_pblri 893 WP Local time, oltime=",oltime

         end select

!        Get information from surface file necessary for conventional data here
         call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)
         !-----------------------
         ! Sfc Elevation from obs (CALIPSO, MPLNET, RAOB) or model at obs location (CATS)
         !-----------------------
         !pblrielev=sfc_elev(i)
         select case(kx)
            case(890)
               pblrielev=sfc_elev(i)
            case(891) ! no info from obs file (CATS)
               !pblrielev=r_missing
               pblrielev=zz ! surface height (zz) from model at obs location
            case(892)
               pblrielev=surface_alti(i) ! transceiver height above sea level 
            case(120)
               pblrielev=stn_elev(i)
            case(893)
               pblrielev=surface_alti(i)
         end select

         !--------------------------------------------------
         ! Subtract terrain height from CALIPSO
         ! in order to match model pblri (above sfc height)
         ! CALIPSO = height above mean sea level
         !--------------------------------------------------
         ! Subtract sfc height from CALIPSO (will be done in moving average and thinning)
         !if (kx==890) then 

         !   pblri_obs_kind = pblri(i) - pblrielev
         !   print*, "YEG: kx=",kx,",pblri(i)=",pblri(i)
         !   print*, "YEG: pblrielev=",pblrielev
         !   print*, "YEG: pblri-elev=",pblri_obs_kind

         !else

         !   pblri_obs_kind=pblri(i)

         !end if
         pblri_obs_kind=pblri(i)
         !------------------------------
         ! Obs sfc type
         ! Make the consistent land_water_mask for all obs types
         ! pblriosfc: 1=Land, 2=water, 3=ice and snow, 4=mixed (coastlines)
         !------------------------------

         select case(kx)

            case(890) ! CALIPSO
            !0=shallow ocean, 1=land, 2=coastlines, 3=shallow inland water, 4=intermittent water, 
            !5=deep inland water, 6=continental ocean, 7=deep ocean"

               if (sfc_mask_calipso(i) == 1) then ! land
                  pblriosfc=1
               else if (sfc_mask_calipso(i) == 2) then ! coastlines
                  pblriosfc=4
               else ! water 
                  pblriosfc=2
               end if

            case(891) ! CATS
            !; 17=water; 1=evergreen needleleaf forest, 2=evergreen broadleaf forest,
            !; 3=deciduous needleleaf forest, 4=deciduous broadleaf forest, 5=mixed forest,
            !; 6=closed shrublands, 7=open shrublands, 8=woody savannas, 9=savannas,
            !; 10=grasslands, 11=permanent wetlands, 12=croplands, 13=urban,
            !; 14=Cropland/Natural Vegetation Mosaic, 15=Permanent Snow and Ice,
            !; 16=Barren or Sparsely Vegetated
               if (sfc_mask_cats(i) == 17) then ! water
                  pblriosfc=2
               else if (sfc_mask_cats(i) == 15) then ! permanent snow and ice
                  pblriosfc=3
               else ! land (for CATS, there is no coastline category for sfc type)
                  pblriosfc=1
               end if

            case(892) ! mplnet
               pblriosfc=r_missing ! Check if this is okay !!!!!!!

            case(120) ! raob
               pblriosfc=r_missing ! Check if this is okay !!!!!!!

            case(893) ! wp
               pblriosfc=r_missing ! Check if this is okay !!!!!!!

         end select

         !----------------------
         ! QC
         !----------------------
         ! Initialize pblriqm
         pblriqm=0

         !----------
         ! Raob
         ! 1) assimilating --> qc=0,5,50
         !----------
         if (kx == 120) then
            pblriqm=qc_flag(i) ! qc_flag=0: good, qc_flag=5: # lev < 10 below 500 hPa, qc_flag=50: 2nd lev = pblh, but we will not specifically exclude data with qc_flags 5 and 50.
            ! not use obs for station elevation > 3 km
            if (pblrielev > 3000.0_r_kind) then
               pblriqm=9
            ! not use obs when pblri_obs < 50 m and (pblh_qgrad minus pblri_obs) > 1km
            !else if (0<pblri_obs_kind<50.0_r_kind .and. (pblh_q_grad(i)-pblri_obs_kind)>1000.0_r_kind ) then
            !   pblriqm=10
            end if
         end if

         !----------
         ! MPLNET
         !----------
         if (kx == 892) then
            pblriqm=flag_use_cloud(i) ! mplnet (all data = 1 for flag_use_cloud, only cloud-free=0 for flag_use_cloud)
            ! QC for local time (no DA 6PM-10AM, time >= 6 pm or < 10 am) 
            if (mod(oltime,100) < 10 .or. mod(oltime,100) > 17) then
               pblriqm=9
            end if
         end if

         !------------
         ! CATS
         !------------
         !For CATS, recommended to use Qflag > 1 (do not use flag < 2).
         if (kx == 891) then
            if (qc_flag(i) < 2) then
               pblriqm=9
            ! QC for night-time over land
            else if (mod(oltime,100) < 10 .or. mod(oltime,100) > 17) then ! night-time
               ! land
               if (pblriosfc /= 2) then ! pblriosfc = 2 (sfc_mask_cats(i)=17) : water
                                        ! pblriosfc /= 2 : land
                  pblriqm=10
                  print*, "yeg_cats_qc: sfc=",pblriosfc
                  print*, "yeg_cats_qc: mod(oltime,100)=",mod(oltime,100)
               end if
            end if

         end if

         !------------
         ! RWP
         !------------
         ! Radar wind profiler: SNR qc flag = 0 good obs
         !if (kx == 893 .and. negative_snr_indicator(i).ne.0) then ! WP QC 0: good flag
         if (kx == 893) then
            if (negative_snr_indicator(i).ne.0) then
               pblriqm=9
            ! QC for local time (no DA 6PM-10AM, time >= 6 pm or < 10 am)
            else if (mod(oltime,100) < 10 .or. mod(oltime,100) > 17) then
               pblriqm=10
            end if
         end if

         !----------
         ! CALIPSO
         !----------
         if (kx == 890) then ! calipso
            if (pblriosfc == 4) then ! coastline from obs data
               pblriqm=9
            !QC for CALIPSO: Latitude 
            !only consider -60 < lat < 60 for calipso
            else if (lat(i)>60.0_r_kind .or. lat(i) < (-1.0)*60.0_r_kind) then
               pblriqm=10
               print*, "yeg_calipso_qc: lat(i)=",lat(i)
            !We are not considering localtime for CALIPSO any more.
            end if
         end if

         !----------
         ! ALL OBS
         !----------
         !?????? Discuss 6000 is proper or not.
         if (pblri_obs_kind .le. 0.0 .or. pblri_obs_kind .gt. 6000_r_kind) then
            pblriqm=15
            print*, 'pblriqm=15, pblh<0 or pblh>6000, pblri_obs_kind=',pblri_obs_kind
         else if (idomsfc>=3) then ! dominate model sfc type >=3 <-- mixed if any surrounding grid has different sfc type (water(0),land(1),ice(2))
            if (kx.ne.120) then ! please comment out mixed sfc type QC for raob pblh for inflating ens spread experiment.
               pblriqm=20
               print*, 'pblriqm=20, idomsfc>=3 (mixed sfc type), idomsfc=',idomsfc
            end if
         end if

         if (pblri_obs_kind <= zero) then
            pblriob=r_missing
         else
            pblriob=pblri_obs_kind
         end if
         !------------------------
         ! Take natural logarithm
         !------------------------
         !if (pblri_obs_kind <= zero) then
         !   pblriob=r_missing
         !else if (pblri_obs_kind > zero .and. pblri_obs_kind <= one) then
         !   pblri_obs_kind=one
         !   pblriob=log(pblri_obs_kind) ! pblriob = 0
         !else 
         !   pblriob=log(pblri_obs_kind)
         !end if

!        if (nkx==131 .or. nkx==133 .or. nkx==135) then
!          anal_time=0
!          obs_time=0
!          tmp_time=zero
!          tmp_time(2)=timeobs
!          anal_time(1)=iadate(1)
!          anal_time(2)=iadate(2)
!          anal_time(3)=iadate(3)
!          anal_time(5)=iadate(4)
!          call w3movdat(tmp_time,anal_time,obs_time) ! observation time

!          lobs_time=0
!          tmp_time=zero
!          cenlon_tmp=hdr(2)
!          if (hdr(2) > 180.0) cenlon_tmp=hdr(2)-360.0_r_kind
!          tmp_time(2)=cenlon_tmp/15.0_r_kind
!          call w3movdat(tmp_time,obs_time,lobs_time) ! local observation time
!          ltime = lobs_time(5)+lobs_time(6)/60.0_r_kind+lobs_time(7)/3600.0_r_kind
!          if ((ltime.gt.21.0) .and. (ltime.lt.5.0)) pblriqm=3
!        end if

!        Set usage variable
         usage = 0.
         print*, 'YEG18: icuse(ikx)=',icuse(ikx), ',ikx=',ikx
         if(icuse(ikx) <= 0) usage=150.
         ! qm=15 for all obs (when pblh < 0 or pblh > 6000)
         ! qm=20 for all obs (when model sfc type = mixed)
         ! qm=9 (stn_elev>3km),10 (pblri_obs<50 m and pblh_qgrad minus pblri_obs > 1km) for RAOB
         ! qm=9 (coastlines),10 (high lat) for CALIPSO
         ! qm=9 (qc flag),10 (night-time land) for CATS
         ! qm=9 (night-time) for MPLNET
         ! qm=9 (qc flag),10 (night-tiime) for wind profiler
         if(pblriqm == 9 .or. pblriqm == 10 .or. pblriqm == 11 .or. pblriqm == 15 .or. pblriqm == 20) usage=150.
         print*, 'YEG19: usage=',usage

!        Set inflate_error logical 
         inflate_error=.false.
         if (pblriqm == 3) inflate_error=.true.

!--Outputs

         ndata=ndata+1
         nodata=nodata+1
         iout=ndata

         write(6,*) 'YEG1: maxobs, ndata, nodata, iout =',maxobs,ndata,nodata,iout

         if(ndata > maxobs) then
            write(6,*)'READ_PBLRI:  ***WARNING*** ndata > maxobs for ',obstype
            ndata = maxobs
         end if

         write(6,*) 'YEG2: maxobs, ndata, nodata, iout =',maxobs, ndata,nodata,iout

!        setup for sort of averaged obs 
         !---------------
         ! Obs Error
         !---------------
         !pblrioe=0.1_r_kind  ! temporarily 100/1400
         !pblrioe=1.31_r_kind  ! temporarily 100/1400
         !pblrioe=1.5_r_kind ! 2.0_r_kind  ! Mar/10/2023
         !pblrioe=0.075_r_kind ! 2.0_r_kind  ! Mar/10/2023

         !pblrioe=0.1_r_kind*pblriob ! 1.5_r_kind was used
         !pblrioe=0.05_r_kind*pblriob ! 1.5_r_kind was used
         if (pblriob == r_missing) then
            pblrioe = r_missing
         else
            if(pblriob <= cpblh1) then
               pblrioe = cerr1
            else if(pblriob > cpblh1 .and. pblriob < cpblh2) then
               pblrioe = cerr1 + (pblriob-cpblh1)* &
                           (cerr2-cerr1)/(cpblh2-cpblh1)
            else 
               pblrioe = cerr2
            endif         
         endif         
    
         if (nkx==120) then
            ! Inflate obs error of raob for station elevation > 1 km
            if (pblrielev > 1000.0_r_kind) then
               !pblrioe=1.44_r_kind ! inflate obs error to be 1.2 times larger (1.2*1.2)
               !pblrioe=1.2_r_kind ! inflate obs error to be 1.2 times larger (1.0*1.2)
               pblrioe=1.2_r_kind*pblrioe ! inflate obs error to be 1.2 times larger (1.0*1.2)
!            else
!               !pblrioe=1.0_r_kind ! 1.2 ! 1.02_r_kind ! 0.5_r_kind ! 1.02_r_kind ! 0.96_r_kind
            end if
         end if

         if (nkx==892) then ! MPLNET 

            !pblrioe=1.0_r_kind ! smaller than BEC PBLRI (1.2)
!            pblrioe=0.6_r_kind ! smaller than BEC PBLRI (1.2): latest
            ! inflate obs error of MPLNET according to no. of obs in DA window at each station
            if (pblri_no_hr_in_da_window(i) > 1.0) then
               !sqrt 2 to 6 -> 1.4, 1.7, 2.0, 2.2, 2.4
               pblrioe=pblrioe*sqrt(pblri_no_hr_in_da_window(i))
               print*, 'mplnet_oerror=',pblrioe
               print*, 'mplnet_no_obs=',pblri_no_hr_in_da_window(i)
            end if
         end if

         if (inflate_error) pblrioe=pblrioe*1.05_r_kind

         ! Station id for raob
         if (nkx==120) then
            c_station_id=sid(i)
         end if
         ! Site for mplnet
         if (nkx==892) then
            c_instrument=instrument(i)
         end if
         ! Station id for wp
         if (nkx==893) then
            c_station_id_wp=sid_wp(i)
         end if

         cdata_all(1,iout)=pblrioe                 ! pblri error (cb)
         cdata_all(2,iout)=dlon                    ! grid relative longitude
         cdata_all(3,iout)=dlat                    ! grid relative latitude
         cdata_all(4,iout)=pblrielev               ! pblri obs elevation
         cdata_all(5,iout)=pblriob                 ! pblri obs
         if (nkx==120) then
            cdata_all(6,iout)=r_station_id            ! index of station id
         else if (nkx==892) then
            cdata_all(6,iout)=r_instrument            ! index of instrument (e.g. MPL44104)
         else if (nkx==893) then
            cdata_all(6,iout)=r_station_id_wp         ! index of station id
         else
            cdata_all(6,iout)=r_missing               ! index of station id
         end if
         cdata_all(7,iout)=t4dv                    ! time
         cdata_all(8,iout)=ikx                     ! type
         print*, 'YEG_read_pblri L1087: ikx=',ikx
         cdata_all(9,iout)=oltime                  ! local time (YYYYMMDDHH) for all obs types!!!!!!!!!!!!!!!!
         cdata_all(10,iout)=pblriqm                ! quality mark
         print*, 'YEG_read_pblri L1090: pblriqm=',pblriqm,',usage=',usage
         cdata_all(11,iout)=usage                  ! usage parameter
         cdata_all(12,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
         cdata_all(13,iout)=dlat_earth_deg         ! earth relative latitude (degrees)

         !sfc type
         cdata_all(14,iout)=pblriosfc
         cdata_all(15,iout)=idomsfc
         if (nkx==120) then
            cdata_all(16,iout)=pblh_q_grad(i)         ! q-gradient based PBLH
         else
            cdata_all(16,iout)=r_missing
         end if
         cdata_all(17,iout)=kx                     ! 120 or 890-893

      end do
!     Normal exit

      ilat=3
      ilon=2
      write(6,*) 'YEG3: ndata, nodata, iout =',ndata,nodata,iout

!     Write observation to scratch file
      call count_obs(ndata,nreal,ilat,ilon,cdata_all,nobs)
      write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
      write(lunout) ((cdata_all(j,i),j=1,nreal),i=1,ndata)
      deallocate(cdata_all)
 
      deallocate(lat, lon, pblri)
      deallocate(ryear, rmonth, rdate, rhour, rminute)
      select case(kx)
         case(890) ! calipso
            deallocate(sfc_elev, sfc_mask_calipso, liadr_data_alt, ATB)
         case(891) ! cats
            deallocate(loc_hour, qc_flag, day_night_flag, sfc_mask_cats)
         case(892) ! mplnet
            deallocate(surface_alti, site, site_org,instrument,instrument_org)
            deallocate(lat_hr,lon_hr,pblri_hr,ryear_hr,rmonth_hr,rdate_hr,rhour_hr,rminute_hr)
            deallocate(surface_alti_hr,flag_use_cloud_hr,flag_use_cloud,site_hr,instrument_hr)
         case(120) ! raob
            !deallocate(nplevels, nplevels500, stn_elev, loc_time, qc_mark_pblh, sid, sid_org)
            deallocate(nplevels, nplevels500, stn_elev, loc_time, qc_flag, sid, sid_org)
            deallocate(pblh_q_grad)
         case(893) ! wp
            deallocate(surface_alti, rsecond)
            deallocate(cloud_fraction, error, negative_snr_indicator, sid_wp, sid_wp_org)
      end select

      end subroutine read_pblri_calipso_cats_mplnet_raob_wp

      subroutine read_pblri_gnssro(nread,ndata,nodata,infile,obstype,lunout,twindin,&
         sis,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_pblri_gnssro     read obs from msgs in gnssro data files
!
!   input argument list:
!     infile   - unit from which to read BUFR data
!     obstype  - observation type to process
!     lunout   - unit to which to write data for further processing
!     hgtl_full- 3d geopotential height on full domain grid
!
!   output argument list:
!     nread    - number of type "obstype" observations read
!     nodata   - number of individual "obstype" observations read
!     ndata    - number of type "obstype" observations retained for further processing
!     twindin  - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
      use kinds, only: r_kind,r_double,i_kind,r_single
      use constants, only: zero,one_tenth,one,deg2rad,rad2deg,three,rearth, &
           grav_equator,eccentricity,somigliana,grav_ratio,grav, &
           semi_major_axis,flattening,two
      use gridmod, only: diagnostic_reg,regional,nlon,nlat,nsig, &
           tll2xy,txy2ll,rotate_wind_ll2xy,rotate_wind_xy2ll,&
           rlats,rlons
      use convinfo, only: nconvtype,ctwind, &
           icuse,ictype,ioctype
      use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen
      use obsmod, only: iadate,offtime_data,bmiss
      use deter_sfc_mod, only: deter_sfc2
      use mpimod, only: npe
      implicit none

!     Declare passed variables
      character(len=*),intent(in):: infile,obstype
      character(20),intent(in):: sis
      integer(i_kind),intent(in):: lunout
      integer(i_kind),intent(inout):: nread,ndata,nodata
      integer(i_kind),dimension(npe),intent(inout):: nobs
      real(r_kind),intent(in):: twindin

     end subroutine read_pblri_gnssro

     subroutine read_pblri_wp_text(nread,ndata,nodata,infile,obstype,lunout,twindin,&
         sis,hgtl_full,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_pblri     read obs from msgs in wind profiler text data files
!
! program history log:
!   2021 zhu
!
!   input argument list:
!     infile   - unit from which to read BUFR data
!     obstype  - observation type to process
!     lunout   - unit to which to write data for further processing
!     hgtl_full- 3d geopotential height on full domain grid
!
!   output argument list:
!     nread    - number of type "obstype" observations read
!     nodata   - number of individual "obstype" observations read
!     ndata    - number of type "obstype" observations retained for further processing
!     twindin  - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
      use kinds, only: r_kind,r_double,i_kind,r_single
      use constants, only: zero,one_tenth,one,deg2rad,rad2deg,three,rearth, &
           grav_equator,eccentricity,somigliana,grav_ratio,grav, &
           semi_major_axis,flattening,two
      use gridmod, only: diagnostic_reg,regional,nlon,nlat,nsig, &
           tll2xy,txy2ll,rotate_wind_ll2xy,rotate_wind_xy2ll,&
           rlats,rlons
      use convinfo, only: nconvtype,ctwind, &
           icuse,ictype,ioctype
      use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen
      use obsmod, only: iadate,offtime_data,bmiss
      use deter_sfc_mod, only: deter_sfc2
      use mpimod, only: npe
      implicit none

!     Declare passed variables
      character(len=*),intent(in):: infile,obstype
      character(20),intent(in):: sis
      integer(i_kind),intent(in):: lunout
      integer(i_kind),intent(inout):: nread,ndata,nodata
      integer(i_kind),dimension(npe),intent(inout):: nobs
      real(r_kind),intent(in):: twindin
      real(r_kind),dimension(nlat,nlon,nsig),intent(in):: hgtl_full

!     Declare local parameters
      real(r_kind),parameter:: r360 = 360.0_r_kind

      integer(i_kind) lunin,idate

      real(r_kind),allocatable,dimension(:,:):: cdata_all

      character(10) date
      logical first,outside,inflate_error,lexist,more_data
      integer(i_kind) ihh,idd,iret,im,iy,istat,i,j,k
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      integer(i_kind) kk,klon1,klat1,klonp1,klatp1
      integer(i_kind) ntest,nchanl
      integer(i_kind) pblriqm,maxobs,idomsfc
!     integer(i_kind),dimension(8):: obs_time,anal_time,lobs_time
!     real(r_kind) ltime,cenlon_tmp
      real(r_kind) usage
      real(r_kind) pblriob,pblrioe,pblrielev
      real(r_kind) dlat,dlon,dlat_earth,dlon_earth,stnelev
      real(r_kind) dlat_earth_deg,dlon_earth_deg
      real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
      real(r_kind) :: tsavg,ff10,sfcr,zz
!     real(r_kind),dimension(5):: tmp_time
      real(r_kind) sin2,termg,termr,termrg
      real(r_kind),dimension(nsig):: zges,hges
      real(r_kind) dx,dy,dx1,dy1,w00,w10,w01,w11

      integer(i_kind) idate5(5),minobs,minan
      real(r_kind) time_correction,timeobs,time,toff,t4dv,zeps

!     real(r_single) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblri_wp
      real(r_kind) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblri_wp
      data lunin /50/

!     Initialize variables
      nreal=14
      ntest=0
      kx= 888

!     Check if pblri file exists
      inquire(file=trim(infile),exist=lexist)
      if (.not.lexist) return

      open(lunin,file=trim(infile),form='formatted')
      read(lunin, *) idate

      maxobs=0
      first=.true.
      more_data=.true.
      readfile1: do while(more_data)
         read(lunin,*,iostat=istat) stationid, lat_deg, lon_deg, altitude, localtime, utctime, localday, utcday, pblri_wp
         if (istat/=0) then
            more_data=.false.
         else
            maxobs=maxobs+1
         end if

!        Time offset
         if (first) then
            call time_4dvar(idate,toff)
            first=.false.
         end if
      end do readfile1

      if (maxobs == 0) return

      allocate(cdata_all(nreal,maxobs))
      nread=0
      nchanl=0
      ilon=2
      ilat=3
      close(lunin)
      open(lunin,file=trim(infile),form='formatted')
      read(lunin, *) idate
      print*, 'read_pblri idate=',idate,' maxobs=',maxobs
      do i=1,maxobs

         nread=nread+1
         if(kx == 120) nkx= 120
         if(kx == 227) nkx= 181
         if(kx == 888) nkx= 888

!        Is pblri in convinfo file
         ikx=0
         do j=1,nconvtype
            if(kx == ictype(j)) then
               ikx=j
               exit
            end if
         end do
         if(ikx == 0) cycle

         read(lunin,*) stationid, lat_deg, lon_deg, altitude, localtime, utctime, localday, utcday, pblri_wp
         print*, 'read_pblri:', stationid, lat_deg, lon_deg, altitude, localtime, utctime, localday, utcday, pblri_wp

         if(lon_deg>= r360)lon_deg=lon_deg-r360
         if(lon_deg < zero)lon_deg=lon_deg+r360
         dlon_earth_deg=lon_deg
         dlat_earth_deg=lat_deg
         dlon_earth=lon_deg*deg2rad
         dlat_earth=lat_deg*deg2rad
         if(regional)then
            call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)    ! convert to rotated coordinate
            if(diagnostic_reg) then
               call txy2ll(dlon,dlat,rlon00,rlat00)
               ntest=ntest+1
               cdist=sin(dlat_earth)*sin(rlat00)+cos(dlat_earth)*cos(rlat00)* &
                    (sin(dlon_earth)*sin(rlon00)+cos(dlon_earth)*cos(rlon00))
               cdist=max(-one,min(cdist,one))
               disterr=acos(cdist)*rad2deg
               disterrmax=max(disterrmax,disterr)
            end if
            if(outside) cycle   ! check to see if outside regional domain
         else
            dlat = dlat_earth
            dlon = dlon_earth
            call grdcrd1(dlat,rlats,nlat,1)
            call grdcrd1(dlon,rlons,nlon,1)
         endif

!        Interpolate guess pressure profile to observation location
         klon1= int(dlon);  klat1= int(dlat)
         dx   = dlon-klon1; dy   = dlat-klat1
         dx1  = one-dx;     dy1  = one-dy
         w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy

         klat1=min(max(1,klat1),nlat); klon1=min(max(0,klon1),nlon)
         if (klon1==0) klon1=nlon
         klatp1=min(nlat,klat1+1); klonp1=klon1+1
         if (klonp1==nlon+1) klonp1=1
         do kk=1,nsig
            hges(kk)=w00*hgtl_full(klat1 ,klon1 ,kk) +  &
                     w10*hgtl_full(klatp1,klon1 ,kk) + &
                     w01*hgtl_full(klat1 ,klonp1,kk) + &
                     w11*hgtl_full(klatp1,klonp1,kk)
         end do
         sin2  = sin(dlat_earth)*sin(dlat_earth)
         termg = grav_equator * &
            ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
         termr = semi_major_axis /(one + flattening + grav_ratio -  &
            two*flattening*sin2)
         termrg = (termg/grav)*termr
         do k=1,nsig
            zges(k) = (termr*hges(k)) / (termrg-hges(k))
         end do


         if(offtime_data) then

!          in time correction for observations to account for analysis
!                time being different from obs file time.
           write(date,'( i10)') idate
           read (date,'(i4,3i2)') iy,im,idd,ihh
           idate5(1)=iy
           idate5(2)=im
           idate5(3)=idd
           idate5(4)=ihh
           idate5(5)=0
           call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
           idate5(1)=iadate(1)
           idate5(2)=iadate(2)
           idate5(3)=iadate(3)
           idate5(4)=iadate(4)
           idate5(5)=0
           call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date

!          Add obs reference time, then subtract analysis time to get obs time relative to analysis

           time_correction=float(minobs-minan)/60._r_kind

         else
           time_correction=zero
         end if


!        Time check
         timeobs=real(utctime, r_double)
         time=timeobs + time_correction
         t4dv=timeobs + toff
         zeps=1.0e-8_r_kind
         if (t4dv<zero  .and.t4dv>      -zeps) t4dv=zero
         if (t4dv>winlen.and.t4dv<winlen+zeps) t4dv=winlen
         t4dv=t4dv + time_correction
         if (l4dvar.or.l4densvar) then
           if (t4dv<zero.OR.t4dv>winlen) cycle
         else
           if((real(abs(time)) > real(ctwind(ikx)) .or. real(abs(time)) > real(twindin))) cycle 
         end if
         print*, 'read_pblri: idate, iadate, timeobs, toff, t4dv=', idate, iadate, timeobs, toff, t4dv

         stnelev=altitude
         pblrielev=stnelev
         pblriob=pblri_wp
         pblriqm=0
         if (pblriob .lt. 0.0) pblriqm=15

!        if (nkx==131 .or. nkx==133 .or. nkx==135) then
!          anal_time=0
!          obs_time=0
!          tmp_time=zero
!          tmp_time(2)=timeobs
!          anal_time(1)=iadate(1)
!          anal_time(2)=iadate(2)
!          anal_time(3)=iadate(3)
!          anal_time(5)=iadate(4)
!          call w3movdat(tmp_time,anal_time,obs_time) ! observation time

!          lobs_time=0
!          tmp_time=zero
!          cenlon_tmp=hdr(2)
!          if (hdr(2) > 180.0) cenlon_tmp=hdr(2)-360.0_r_kind
!          tmp_time(2)=cenlon_tmp/15.0_r_kind
!          call w3movdat(tmp_time,obs_time,lobs_time) ! local observation time
!          ltime = lobs_time(5)+lobs_time(6)/60.0_r_kind+lobs_time(7)/3600.0_r_kind
!          if ((ltime.gt.21.0) .and. (ltime.lt.5.0)) pblriqm=3
!        end if

!        Set usage variable
         usage = 0.
         if(icuse(ikx) <= 0) usage=150.
         if(pblriqm == 15 .or. pblriqm == 9) usage=150.

!        Set inflate_error logical 
         inflate_error=.false.
         if (pblriqm == 3) inflate_error=.true.

!--Outputs

         ndata=ndata+1
         nodata=nodata+1
         iout=ndata

         if(ndata > maxobs) then
            write(6,*)'READ_PBLH:  ***WARNING*** ndata > maxobs for ',obstype
            ndata = maxobs
         end if

!        Get information from surface file necessary for conventional data here
         call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)

!        setup for sort of averaged obs 
         pblrioe=100.0_r_kind  ! temporarily
         if (nkx==120) pblrioe=50.
         if (nkx==181) pblrioe=100.
         if (inflate_error) pblrioe=pblrioe*1.5_r_kind

         cdata_all(1,iout)=pblrioe                  ! pblri error (cb)
         cdata_all(2,iout)=dlon                    ! grid relative longitude
         cdata_all(3,iout)=dlat                    ! grid relative latitude
         cdata_all(4,iout)=pblrielev                ! pblri obs elevation
         cdata_all(5,iout)=pblriob                  ! pblri obs
         cdata_all(6,iout)=stationid               ! station id
         cdata_all(7,iout)=t4dv                    ! time
         cdata_all(8,iout)=ikx                     ! type
         cdata_all(9,iout)=localtime               ! obs local time
         cdata_all(10,iout)=pblriqm                 ! quality mark
         cdata_all(11,iout)=usage                  ! usage parameter
         cdata_all(12,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
         cdata_all(13,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
         cdata_all(14,iout)=altitude               ! station elevation (m)

      end do
      close(lunin)
!   Normal exit

!   Write observation to scratch file
     call count_obs(ndata,nreal,ilat,ilon,cdata_all,nobs)
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((cdata_all(j,i),j=1,nreal),i=1,ndata)
     deallocate(cdata_all)

     if (ndata == 0) then
        close(lunin)
        write(6,*)'READ_pblh:  close(',lunin,')'
     endif

     close(lunin)
     end subroutine read_pblri_wp_text

     subroutine convert_localtime_hr(ltime_mm, ltime_dd, ltime_hh, interval, ana_mm, ana_dd, ana_hh)
     use kinds, only: i_kind
     integer(i_kind), intent(out) :: ltime_mm, ltime_dd, ltime_hh
     integer(i_kind), intent(in ) :: interval, ana_mm, ana_dd, ana_hh
     ltime_mm = ana_mm
     ltime_dd = ana_dd
     ltime_hh = ana_hh + interval

     ! localtime = UTC + hr
     if (ltime_hh >= 24) then

        ltime_hh = ltime_hh - 24
        ltime_dd = ana_dd + 1
        if (ana_mm == 8 .and. ltime_dd > 31) then
           ltime_mm = 9
           ltime_dd = 1
        end if
        if (ana_mm == 9 .and. ltime_dd > 30) then
           ltime_mm = 10
           ltime_dd = 1
        end if

     elseif (ltime_hh < 0) then

        ltime_hh = 24 + ltime_hh
        ltime_dd = ana_dd - 1
        if (ana_mm == 9 .and. ana_dd == 1) then
           ltime_mm = 8
           ltime_dd = 31
        end if

     end if
     end subroutine convert_localtime_hr

     subroutine convert_localtime_min(ltime_mm, ltime_dd, ltime_hh, ltime_min, interval_hr_min, ana_mm, ana_dd, ana_hh,ana_min)
     use kinds, only: i_kind,r_kind
     integer(i_kind), intent(out) :: ltime_mm, ltime_dd, ltime_hh, ltime_min
     integer(i_kind), intent(in ) :: ana_mm, ana_dd, ana_hh, ana_min
     real(r_kind), intent(in ) :: interval_hr_min
     real(r_kind) :: ana_hh_real, ltime_hh_real

     ! localtime = UTC + interval_hr_min
     ltime_mm = ana_mm
     ltime_dd = ana_dd
     ! 1) Convert HH:MM to HH.xx
     !ana_hh_real: ana_hh + float(ana_min) to account for minutes
     ana_hh_real = real(ana_hh) + real(ana_min)/60.0
     ! 2) UTC + interval => Local time [hr]
     !ltime_hh = ana_hh + interval_hr_min
     ltime_hh_real = ana_hh_real + interval_hr_min
     ! 3) local time [hr] -> integer HH:MM
     ! 0.9 -> 0, 1.9 -> 1
     ! for negative, -0.9 -> -1, -1.9 -> -2
     ltime_hh = floor(ltime_hh_real) ! 
     ! 4) local time [minutes] -> integer HH:MM
     ltime_min = nint((ltime_hh_real - ltime_hh)*60.0) ! Round (1.6 -> 2)

     ! only for Aug and Sep
     if (ltime_hh >= 24) then

        ltime_hh = ltime_hh - 24
        ltime_dd = ana_dd + 1
        if (ana_mm == 8 .and. ltime_dd > 31) then
           ltime_mm = 9
           ltime_dd = 1
        end if
        if (ana_mm == 9 .and. ltime_dd > 30) then
           ltime_mm = 10
           ltime_dd = 1
        end if

     elseif (ltime_hh < 0) then

        ltime_hh = 24 + ltime_hh
        ltime_dd = ana_dd - 1
        if (ana_mm == 9 .and. ana_dd == 1) then
           ltime_mm = 8
           ltime_dd = 31
        end if

     end if
     end subroutine convert_localtime_min


     subroutine count_elev_mvavg_thin_obs(nobs,lat,lon,sfc_elev,pblri,mesh,mvavg_rad,nthinobs,lat_thin,lon_thin,pblri_mvavg_thin)
     use kinds, only: i_kind,r_kind
     use buddycheck_mod, only: gc_dist
     implicit none
     integer(i_kind), intent(in)  :: nobs ! the number of obs before thinning
     real(r_kind), dimension(nobs), intent(in)  :: lat, lon
     real(r_kind), dimension(nobs), intent(in)  :: sfc_elev,pblri
     real(r_kind), intent(in)  :: mesh ! thinning size (km)
     real(r_kind), intent(in)  :: mvavg_rad! moving average radius
     integer(i_kind), intent(out) :: nthinobs ! the number of obs after thinning
     real(r_kind), dimension(:), allocatable, intent(out) :: lat_thin, lon_thin
     real(r_kind), dimension(:), allocatable, intent(out)  :: pblri_mvavg_thin

     integer(i_kind) :: i,j,k
     integer(i_kind) :: begin,ending
     integer(i_kind) :: no_obs_mvavg
     real(r_kind) :: o1_lat, o1_lon, o2_lat, o2_lon, dist
     real(r_kind) :: tobs_mvavg, sum_pblri

     ! FOR CALIPSO 
     ! Count the number of obs after
     ! 1) Thinning depending on distance (mesh)
     ! 2) Substracting sfc_elev from pblri CALIPSO
     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs
 
     ! 1) Thinning depending on distance (mesh)
     ! Count the number of obs after thinning depending on distance (mesh)
     do i = 1, nobs
        if (i == 1) then
           o1_lat = lat(i)
           o1_lon = lon(i)
           nthinobs=1
        else
           o2_lat = lat(i)
           o2_lon = lon(i)
           !function gc_dist(inlat1,inlon1,inlat2,inlon2) [meters]
           !lat1,lon1,lat2,lon2 in degrees
           ! *0.001_r_kind
           dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind
           !print*,"yeg_thinning, i=",i,",dist[km]=",dist

           if (dist >= mesh) then
              nthinobs = nthinobs+1
              !print*,"yeg_thinning, dist[km]=",dist,",nthinobs=",nthinobs
              o1_lat = o2_lat
              o1_lon = o2_lon
           end if
        end if
     end do

     allocate(lat_thin(nthinobs),lon_thin(nthinobs),pblri_mvavg_thin(nthinobs))

     ! 2) Save thinned obs and substract sfc_elev from pblri
     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs
     ! Save thinned QCed obs (lat,lon,pblri)
     no_obs_mvavg=40 !80 !40 ! moving average --> +- 40 obs along the track (12 km, 330m between obs)
     tobs_mvavg=0 ! total number of obs for moving avg
     sum_pblri=0

     do i = 1, nobs
        if (i == 1) then
           nthinobs=1
           o1_lat = lat(i)
           o1_lon = lon(i)

           lat_thin(nthinobs)=lat(i)
           lon_thin(nthinobs)=lon(i)
           sum_pblri = sum_pblri + pblri(i) - sfc_elev(i) ! subtracting sfc height
           tobs_mvavg = tobs_mvavg + 1
        else
           !--------------------------
           ! Decide thinning or not
           o2_lat = lat(i)
           o2_lon = lon(i)
           !function gc_dist(inlat1,inlon1,inlat2,inlon2) [meters]
           !lat1,lon1,lat2,lon2 in degrees
           ! *0.001_r_kind
           dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind
           !print*,"yeg_thinning, i=",i,",dist[km]=",dist

           if (dist >= mesh) then

              tobs_mvavg=0 ! total number of obs for moving avg
              sum_pblri=0
              nthinobs = nthinobs+1

              lat_thin(nthinobs)=lat(i)
              lon_thin(nthinobs)=lon(i)
              sum_pblri = sum_pblri + pblri(i) - sfc_elev(i)
              tobs_mvavg = tobs_mvavg + 1

              o1_lat = o2_lat
              o1_lon = o2_lon
           else
              cycle
           end if

        end if

        ! Consider near obs for moving average
        if (i <= no_obs_mvavg) then ! j=1,2,3,4,5
           begin=i*(-1)+1
           ending=no_obs_mvavg
        else if (i > nobs-no_obs_mvavg) then  !j=16,17,18,19,20
           begin=no_obs_mvavg*(-1)
           ending=nobs-i
        else
           begin=no_obs_mvavg*(-1)
           ending=no_obs_mvavg
        end if

        ! for +- obs near thinned obs
        do k = begin,ending
           if (k.ne.0) then
              o2_lat = lat(i+k)
              o2_lon = lon(i+k)
              dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind
              if (dist < mvavg_rad) then ! select obs for moving average within radius 12 km
                 sum_pblri = sum_pblri + pblri(i+k) - sfc_elev(i)
                 tobs_mvavg = tobs_mvavg + 1
              end if
           end if
        end do

        ! Save moving average pblh
        pblri_mvavg_thin(nthinobs) = sum_pblri / tobs_mvavg
        print*, 'YEG_CALIPSO_MVAVG: nthinobs=',nthinobs
        print*, 'YEG_CALIPSO_MVAVG: sum_pblri,tobs_mvavg,pblh_avg=',sum_pblri,tobs_mvavg,pblri_mvavg_thin(nthinobs)

     end do


     end subroutine count_elev_mvavg_thin_obs

     subroutine count_qc_mvavg_thin_obs(nobs,lat,lon,qc,pblri,mesh,mvavg_rad,n_qc_thin,lat_qc_thin,lon_qc_thin,pblri_qc_mvavg_thin)
     use kinds, only: i_kind,r_kind
     use buddycheck_mod, only: gc_dist
     implicit none
     integer(i_kind), intent(in)  :: nobs ! the number of obs before thinning
     real(r_kind), dimension(nobs), intent(in)  :: lat, lon, qc, pblri
     real(r_kind), dimension(:), allocatable :: lat_qc, lon_qc, pblri_qc
     real(r_kind), dimension(:), allocatable, intent(out) :: lat_qc_thin, lon_qc_thin
     real(r_kind), dimension(:), allocatable, intent(out)  :: pblri_qc_mvavg_thin
     real(r_kind), intent(in)  :: mesh ! thinning size (km)
     real(r_kind), intent(in)  :: mvavg_rad! moving average radius
     integer(i_kind), intent(out) :: n_qc_thin ! the number of obs after thinning

     integer(i_kind) :: i,j,k
     integer(i_kind) :: begin,ending
     integer(i_kind) :: n_qc, no_obs_mvavg
     real(r_kind) :: o1_lat, o1_lon, o2_lat, o2_lon, dist, dist_prev
     real(r_kind) :: tobs_mvavg, sum_pblri
 
     ! FOR CATS 
     ! Count the number of obs after
     ! 1) Considering qc_flag>1 (good obs)
     ! 2) Thinning depending on distance (mesh)
     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs

     n_qc=0

     ! 1) Save QC>1 (good obs)

     ! Count the number for qc>1 (n_qc), qc<=1 (n_bqc)
     do i = 1, nobs
        if (qc(i) > 1) then
           n_qc=n_qc+1
        end if 
     end do

     allocate(lat_qc(n_qc),lon_qc(n_qc),pblri_qc(n_qc))

     n_qc=0
     ! Create lat, lon, qc for good qc
     do i = 1, nobs
        if (qc(i) > 1) then
           n_qc=n_qc+1
           lat_qc(n_qc) = lat(i)
           lon_qc(n_qc) = lon(i)
           pblri_qc(n_qc) = pblri(i)
        end if 
     end do

     print*, "YEG_CATS_QC_THIN, 1) total_obs,n_qc(qc>1)=",nobs,n_qc

     ! 2) Thinning for qc > 1 obs (good obs)

     ! Count the number of thinned obs for qc>1
     do i = 1, n_qc
        if (i == 1) then
           o1_lat = lat_qc(i)
           o1_lon = lon_qc(i)
           n_qc_thin=1
        else
           o2_lat = lat_qc(i)
           o2_lon = lon_qc(i)
           !function gc_dist(inlat1,inlon1,inlat2,inlon2) [meters]
           !lat1,lon1,lat2,lon2 in degrees
           ! *0.001_r_kind
           dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind

           if (dist >= mesh) then
              n_qc_thin = n_qc_thin+1
              o1_lat = o2_lat
              o1_lon = o2_lon
           end if

        end if
     end do

     allocate(lat_qc_thin(n_qc_thin),lon_qc_thin(n_qc_thin),pblri_qc_mvavg_thin(n_qc_thin))

     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs
     ! Save thinned QCed obs (lat,lon,pblri)
     no_obs_mvavg=5 !10 !5 ! moving average --> +- 5 obs (12 km) along the track
     tobs_mvavg=0 ! total number of obs for moving avg
     sum_pblri=0

     do i = 1, n_qc

        if (i == 1) then
           n_qc_thin=1
           j=i
           o1_lat = lat_qc(i)
           o1_lon = lon_qc(i)

           lat_qc_thin(n_qc_thin) = lat_qc(i)
           lon_qc_thin(n_qc_thin) = lon_qc(i)
           sum_pblri = sum_pblri + pblri_qc(i)
           tobs_mvavg = tobs_mvavg + 1
        else
           !--------------------------
           ! Decide thinning or not
           o2_lat = lat_qc(i)
           o2_lon = lon_qc(i)
           !function gc_dist(inlat1,inlon1,inlat2,inlon2) [meters]
           !lat1,lon1,lat2,lon2 in degrees
           ! *0.001_r_kind
           dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind

           if (dist >= mesh) then

              tobs_mvavg=0 ! total number of obs for moving avg
              sum_pblri=0
              n_qc_thin = n_qc_thin+1
              j=i

              lat_qc_thin(n_qc_thin) = lat_qc(i)
              lon_qc_thin(n_qc_thin) = lon_qc(i)
              sum_pblri = sum_pblri + pblri_qc(i)
              tobs_mvavg = tobs_mvavg + 1

              o1_lat = o2_lat
              o1_lon = o2_lon

           else
              cycle
           end if
        end if

        ! Consider near obs for moving average
        if (j <= no_obs_mvavg) then ! j=1,2,3,4,5
           begin=j*(-1)+1
           ending=no_obs_mvavg
        else if (j > n_qc-no_obs_mvavg) then  !j=16,17,18,19,20
           begin=no_obs_mvavg*(-1)
           ending=n_qc-j
        else
           begin=no_obs_mvavg*(-1)
           ending=no_obs_mvavg
        end if

        ! for +- obs near thinned obs
        do k = begin,ending
           if (k.ne.0) then
              o2_lat = lat_qc(j+k)
              o2_lon = lon_qc(j+k)
              dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind
              if (dist < mvavg_rad) then ! select obs for moving average within radius 12 km
                 sum_pblri = sum_pblri + pblri_qc(j+k)
                 tobs_mvavg = tobs_mvavg + 1
              end if
           end if
        end do

        ! Save moving average pblh
        pblri_qc_mvavg_thin(n_qc_thin) = sum_pblri / tobs_mvavg
        print*, 'YEG_CATS_MVAVG: n_qc_thin=',n_qc_thin
        print*, 'YEG_CATS_MVAVG: sum_pblri,tobs_mvavg,pblh_avg=',sum_pblri,tobs_mvavg,pblri_qc_mvavg_thin(n_qc_thin)

     end do

     print*, "YEG_CATS_QC_THIN, 2) n_qc_thin=",n_qc_thin

     end subroutine count_qc_mvavg_thin_obs

     subroutine hourly_mean_obs(nobs,lat,lon,ryear,rmonth,rdate,rhour,rminute,pblri,& 
                              surface_alti,flag_cloud_screen,site,instrument,use_cloud,&
                              nthinobs,lat_hr,lon_hr,ryear_hr,rmonth_hr,rdate_hr,rhour_hr,rminute_hr,pblri_hr,&
                              surface_alti_hr,flag_use_cloud_hr,site_hr,instrument_hr,pblri_no_hr_in_da_window)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  hourly_mean_obs    calculate hourly mean obs for MPLNET PBLH
!
! program history log:
!   2023-07-04  eyang   -  
!
!   input argument list:
!     nobs          - the number of obs
!     lat
!     lon
!     ryear
!     rmonth
!     rdate
!     rhour
!     rminute
!     pblri           - pbl height
!     surface_alti    - surface altitude [meters]
!     flag_cloud_screen - 1: cloud free 2: cloud_fraction>0% 4: cloud_detection_fail 8: no_could_product
!     site            - site information (e.g., GSFC)
!     instrument      - instrument information (e.g., MPL44104)
!     use_cloud       - use all data (including cloud-affected data) or use only cloud free
!                       0: consider only cloud free 1: consider all data (including cloud affected data)
!   output argument list:
!     nthinobs        - the number of thinned obs (temporally)
!     lat_hr
!     lon_hr
!     ryear_hr
!     rmonth_hr
!     rdate_hr
!     rhour_hr
!     rminute_hr
!     pblri_hr        - pbl height
!     surface_alti_hr - surface altitude [meters]
!     flag_use_cloud_hr - 0: use only cloud free data 1: use all data available (including cloud affected data) 
!     site_hr         - site information (e.g., GSFC)
!     instrument_hr   - instrument information (e.g., MPL44104)
!     pblri_no_hr_in_da_window   - no. of obs in DA window at each station (0 - 6)
!
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

     use kinds, only: i_kind,r_kind,r_double
     use, INTRINSIC :: IEEE_ARITHMETIC
     implicit none
     integer(i_kind), intent(in)  :: nobs ! the number of obs before thinning
     real(r_kind), dimension(nobs), intent(in)  :: lat, lon, ryear,rmonth,rdate,rhour,rminute, pblri
     real(r_kind), dimension(nobs), intent(in)  :: surface_alti,flag_cloud_screen
     integer(i_kind), intent(in)  :: use_cloud
     character(19),dimension(nobs), intent(in)  :: site
     character(8), dimension(nobs), intent(in)  :: instrument
     real(r_kind), dimension(:), allocatable, intent(out)  :: lat_hr, lon_hr, ryear_hr,rmonth_hr,rdate_hr,rhour_hr
     real(r_kind), dimension(:), allocatable, intent(out)  :: rminute_hr, pblri_hr
     real(r_kind), dimension(:), allocatable, intent(out)  :: surface_alti_hr,flag_use_cloud_hr
     character(19),dimension(:), allocatable, intent(out)  :: site_hr
     character(8), dimension(:), allocatable, intent(out)  :: instrument_hr
     real(r_kind), dimension(:), allocatable, intent(out)  :: pblri_no_hr_in_da_window !no. of obs in DA window at each station (0 - 6)
     real(r_kind), dimension(:), allocatable :: lat_qc, lon_qc
     integer(i_kind), intent(out) :: nthinobs ! the number of obs after thinning

     character(19),dimension(nobs) :: site_unique_temp
     character(19),dimension(:), allocatable :: site_unique
     character(8),dimension(:), allocatable :: instrument_unique
     real(r_kind), dimension(6) :: hour_unique
     ! temp array
     real(r_kind), dimension(:,:), allocatable :: pblri_temp, pblri_no
     real(r_kind), dimension(:), allocatable :: pblri_no_hr_in_da_window_temp
     real(r_kind), dimension(:,:), allocatable :: lat_temp, lon_temp, ryear_temp, rmonth_temp, rdate_temp
     real(r_kind), dimension(:,:), allocatable :: surface_alti_temp
     integer(i_kind) :: i,j,k,m
     integer(i_kind) :: site_count, hour_count
     real(r_kind) :: o1_lat, o1_lon

     ! Calculate hourly mean pblh for each station
     ! 1) Find unique sites and hours
     ! 2) collect data for hourly mean depending on sites
     ! 3) Calculate average
     ! 4) Assign hourly mean value to 1-d array

     ! 1) Find unique sites and hours

     site_unique_temp(1)=site(1)
     site_count=1

     do i = 2, nobs

        do j = 1, site_count

           if (site(i)==site_unique_temp(j)) then
              exit
           end if

        end do

        if (j == site_count+1) then
           site_count = site_count + 1
           site_unique_temp(site_count) = site(i)
        end if

     end do

     ! Store unique sites and hours
     ! a) sites and instruments
     allocate(site_unique(site_count))
     allocate(instrument_unique(site_count))
     site_unique(1)=site(1)
     instrument_unique(1)=instrument(1)
     site_count=1

     do i = 2, nobs

        do j = 1, site_count

           if (site(i)==site_unique(j)) then
              exit
           end if

        end do

        if (j == site_count+1) then
           site_count = site_count + 1
           site_unique(site_count) = site(i)
           instrument_unique(site_count) = instrument(i)
        end if

     end do

     ! b) hours
     hour_unique(1)=rhour(1)
     hour_count=1

     do i = 2, nobs

        do j = 1, hour_count

           if (rhour(i)==hour_unique(j)) then
              exit
           end if

        end do

        if (j == hour_count+1) then
           hour_count = hour_count + 1
           hour_unique(hour_count) = rhour(i)
        end if

     end do


     ! 2) collect data for hourly mean depending on sites
     allocate(pblri_temp(site_count, 6))
     allocate(pblri_no(site_count, 6))
     allocate(pblri_no_hr_in_da_window_temp(site_count))
     allocate(lat_temp(site_count,6), lon_temp(site_count,6), ryear_temp(site_count,6),rmonth_temp(site_count,6))
     allocate(rdate_temp(site_count,6),surface_alti_temp(site_count,6))
     pblri_temp = 0.
     pblri_no = 0.

     ! read data 
     do i = 1, nobs ! all obs
        do j = 1, site_count ! sites
           do k = 1, 6 ! hour bin

              if (site(i) == site_unique(j) .and. rhour(i) == hour_unique(k)) then

                 lat_temp(j,k)=lat(i)
                 lon_temp(j,k)=lon(i)
                 ryear_temp(j,k)=ryear(i)
                 rmonth_temp(j,k)=rmonth(i)
                 rdate_temp(j,k)=rdate(i)
                 surface_alti_temp(j,k)=surface_alti(i)

                 if (IEEE_IS_NAN(pblri(i))) then
                    cycle
                 else
                    if (use_cloud==1) then ! use all data (including cloud-affected data)
                       pblri_temp(j, k) = pblri_temp(j,k) + pblri(i)
                       pblri_no(j, k) = pblri_no(j,k) + 1
                    else ! use_cloud==0 (use only cloud-free data)
                       if (flag_cloud_screen(i) == 1) then ! cloud-free
                          pblri_temp(j, k) = pblri_temp(j,k) + pblri(i)
                          pblri_no(j, k) = pblri_no(j,k) + 1
                       end if
                    end if
                 end if

              end if

           end do
        end do
     end do

     ! 3) Calculate average and 
     !    the number of obs in DA window at each station
     !    pblri_no_hr_in_da_window   - no. of obs in DA window at each station (0 - 6)
     pblri_no_hr_in_da_window_temp=0.
     do j = 1, site_count
        do k = 1, 6 ! hour bin

           if (pblri_no(j,k)==0) then
              pblri_temp(j,k) = -9999
           else
              pblri_temp(j,k) = pblri_temp(j,k)/pblri_no(j,k)
              pblri_no_hr_in_da_window_temp(j) = pblri_no_hr_in_da_window_temp(j) + 1
              print*, 'j=',j,',k=',k,',pblrinoinda(j)=',pblri_no_hr_in_da_window_temp(j)
           end if
           print*, 'j=',j,'k=',k,',pblri_no=',pblri_no(j,k)
           print*, 'j=',j,'k=',k,',pblri_temp=',pblri_temp(j,k)

        end do
        print*, 'j=',j,',pblrinoinda(j)=',pblri_no_hr_in_da_window_temp(j)
     end do
              
     ! 4) Assign hourly mean value to 1-d array 

     nthinobs=site_count*6
     allocate(lat_hr(nthinobs),lon_hr(nthinobs),ryear_hr(nthinobs),rmonth_hr(nthinobs),rdate_hr(nthinobs))
     allocate(rhour_hr(nthinobs),rminute_hr(nthinobs),pblri_hr(nthinobs),surface_alti_hr(nthinobs),flag_use_cloud_hr(nthinobs))
     allocate(site_hr(nthinobs),instrument_hr(nthinobs), pblri_no_hr_in_da_window(nthinobs))
     m=1
     do j = 1, site_count
        do k = 1, 6 ! hour bin

           lat_hr(m)=lat_temp(j,k)
           lon_hr(m)=lon_temp(j,k)
           ryear_hr(m)=ryear_temp(j,k)
           rmonth_hr(m)=rmonth_temp(j,k)
           rdate_hr(m)=rdate_temp(j,k)
           rhour_hr(m)=hour_unique(k)
           rminute_hr(m)=30.0
           pblri_hr(m)=pblri_temp(j,k)
           surface_alti_hr(m)=surface_alti_temp(j,k)
           flag_use_cloud_hr(m)=use_cloud ! 1: use all data 0: use only cloud-free data
           site_hr(m)=site_unique(j)
           instrument_hr(m)=instrument_unique(j)
           pblri_no_hr_in_da_window(m)=pblri_no_hr_in_da_window_temp(j)
           m=m+1

        end do
     end do

     end subroutine hourly_mean_obs

  end module read_pblri
