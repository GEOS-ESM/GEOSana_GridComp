c
c $Id$
c $Author$
c
C
C History tape header:
C   Record number 1 (integer values) provides description for history
C   tape output.  See the CCM3 User's Guide for description of each
C   of the individual header variables.
C
      common /comhdi/ 
     $     lenhdi  ,mftyp   ,mfilh    ,mfilth  ,nrbd    ,
     $     maxsiz  ,ndavu   ,mxxx     ,mlon    ,nlonw   ,
     $     morec   ,mlev    ,mtrm     ,mtrn    ,mtrk    ,
     $     nfldh   ,nsteph  ,nstprh   ,nitslf  ,ndbase  ,
     $     nsbase  ,ndcur   ,nscur    ,nbdate  ,nbsec   ,
     $     ncdate  ,ncsec   ,mdt      ,mhisf   ,mfstrt  ,
     $     lenhdc  ,lenhdr  ,mpsig    ,mplat   ,mpwts   ,
     $     mpflds  ,mpcfld  ,mflds(3,pflds)
C
      integer lenhdi       ! Length of header record 1
      integer  mftyp       ! Format code
      integer  mfilh       ! Logical file number
      integer  mfilth      ! Max number of files on history tape
      integer  nrbd        ! Number of records before data records
C
      integer  maxsiz      ! Length of data record for this volume
      integer  ndavu       ! Length of the data record after unpacking
      integer  mxxx        ! Horizontal domain flag
      integer  mlon        ! Number of longitude points per latitude line
      integer  nlonw       ! Number of longitude data values written
C
      integer  morec       ! Number of latitude lines or data records
      integer  mlev        ! Number of vertical levels
      integer  mtrm        ! M spectral truncation parameter
      integer  mtrn        ! N spectral truncation parameter
      integer  mtrk        ! K spectral truncation parameter
C
      integer  nfldh       ! Number of fields on the header
      integer  nsteph      ! Time step number
      integer  nstprh      ! Time step number for the start of this run
      integer  nitslf      ! Time steps since last file was written
      integer  ndbase      ! Base day number for this case
C
      integer  nsbase      ! Base number of seconds for this case
      integer  ndcur       ! Current day corresponding to  NSTEPH
      integer  nscur       ! Current seconds corresponding to NSTEPH
      integer  nbdate      ! Base date (yr mo day) as 6-digit integer
      integer  nbsec       ! Seconds to complete  NODT date
C
      integer  ncdate      ! Current date (yymmdd)
      integer  ncsec       ! Current seconds for date
      integer  mdt         ! Model timestep in seconds
      integer  mhisf       ! Frequency that history files are written
      integer  mfstrt      ! Flag to indicate type of start
C
      integer  lenhdc      ! Length of header record 2
      integer  lenhdr      ! Length of header record 3
      integer  mpsig       ! Pointer to first word of sigma value list
      integer  mplat       ! Pointer to list of latitude lines
      integer  mpwts       ! Pointer to list of Gaussian weights
C
      integer  mpflds      ! Pointer to header field info integer list
      integer  mpcfld      ! Pointer to field info chararacter list
      integer  mflds       ! Array of integer field list information
      integer comhdi(37+3*pflds)
C
#if ( defined SUN )
      static comhdi
#endif
      equivalence (comhdi,lenhdi)
C     
C History tape header, record number 2 (character values)
C
      common /comhdc/
     $     mcase      ,mcstit     ,    
     $     lnhstc     ,ldhstc     ,lthstc    ,lshstc    ,
     $     lnhstf     ,ldhstf     ,lthstf    ,lshstf    ,
     $     lnhsti     ,ldhsti     ,lthsti    ,lshsti    ,
     $     lnhstt     ,ldhstt     ,lthstt    ,lshstt    ,
     $     lnhstvs    ,ldhstvs    ,lthstvs   ,lshstvs   ,
     $     lnhstvo    ,ldhstvo    ,lthstvo   ,lshstvo   ,
     $     mcflds(2,pflds)
C
      character*8 mcase                      ! Case identifier
      character*8 ldhstc  ,lthstc  ,lshstc   ! Current hist file info
      character*8 ldhstf  ,lthstf  ,lshstf   ! First hist file info
      character*8 ldhsti  ,lthsti  ,lshsti   ! Initial hist file info
      character*8 ldhstt  ,lthstt  ,lshstt   ! Boundary hist file info
      character*8 ldhstvs ,lthstvs ,lshstvs  ! SST hist file info
      character*8 ldhstvo ,lthstvo ,lshstvo  ! Ozone hist file info
      character*8 mcflds                     ! Array of character
C                                            !    field list info
      character*80 mcstit               ! Case title MSS names:
      character*80 lnhstc               ! Current history file
      character*80 lnhstf               ! First history file
      character*80 lnhsti               ! Initial history file 
      character*80 lnhstt               ! Boundary history file
      character*80 lnhstvs              ! SST history file 
      character*80 lnhstvo              ! Ozone history file 
      character*8 comhdc(7*10 + 19)     ! /comhdc/ scalars
C
#if ( defined SUN )
      static comhdc
#endif
      equivalence (comhdc,mcase)
C
C History tape header, record number 3 (real values)
C
C Real header record contains sigma,latitude and gaussian weights
C
      common /comhdr/ 
     $     sigapb(2*plev+1) ,siga(2*plev+1) ,sigb(2*plev+1) ,
     $        hdlat(plat)      ,hdwt(plat)
C
      real sigapb      ! Hybrid A + B coefficients
      real siga        ! Hybrid A (pressure) coefficients
      real sigb        ! Hybrid B (sigma) coefficients
      real hdlat       ! Latitude list (South to North)
      real hdwt        ! Gaussian weight list.
C
 
