!
! $Id$
! $Author$
!
!
! Character variables associated with Mass Store pathnames
!
      common/commsc/nhpath  ,nrpath  ,ncdata  ,bndtvs  
      common/commsc/bndtvo  ,nrfil   ,nsmvn   ,nrmvn   ,msscom  
      common/commsc/nswrps  ,lcroot  ,hycofile,topofile
      common/commsc/topofilenc,topopath
!
      character*72  nhpath   ! MSS pathname for history tapes
      character*72  nrpath   ! MSS pathname for restart files
      character*80  ncdata   ! MSS pathname for initial dataset
!
      character*80  bndtvs   ! MSS path for time-variant sst dataset
      character*80  bndtvo   ! MSS path for time-variant ozone dataset
      character*8   nrfil    ! Current file name for regen dataset
      character*8   nsmvn    ! Virtual volume name for history tapes
!
      character*8   nrmvn    ! Virtual volume name for restart data
      character*80  msscom   ! MSS comment field
      character*8   nswrps   ! MSS write password
      character*80  lcroot   ! Prepend to MSS paths for local disk name
      character*80  hycofile ! full pathname for hybrid coef. file
      character*80  topofile ! full pathname for high res topo file
      character*80  topofilenc! full pathname for high res topo file
      character*80  topopath ! full pathname for high res topo file
!
! Non-character variables associated with Mass Store pathnames
!
      common/commss/nhpthl  ,nrpthl  ,irt     ,rirt
!
      integer  nhpthl        ! Length of nhpath,history tape pathname
      integer  nrpthl        ! Length of nrpath,restart data pathname
      integer  irt           ! Mass Store retention period, history tapes
      integer  rirt          ! Mass Store retention time, restart data
!
