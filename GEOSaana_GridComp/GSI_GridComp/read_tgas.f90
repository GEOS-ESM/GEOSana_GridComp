subroutine read_tgas(nread, npuse, nouse, jsatid, infile, gstime, lunout,      &
                     obstype, twind, sis, ithin, rmesh, nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_tgas --- Read trace gas data
!    prgrmmr:    weir             org: gmao                date: 2014-04-14
!
! abstract:      This routine reads trace gas observations. The initial code is
!                taken from read_ozone and read_co.
!
! program history log: 
!   2010-03-30   weir  - initial code
!
!   input argument list:
!     jsatid  - satellite id to read
!     infile  - file from which to read trace gas data
!     gstime  - analysis time in minutes from reference date
!     lunout  - unit to which to write data for further processing
!     obstype - observation type to process
!     twind   - input group time window (hours)
!     sis     - satellite/instrument/sensor indicator
!     ithin   - flag to thin data
!     rmesh   - thinning mesh size (km)
!
!   output argument list:
!     nread   - number of observations read
!     npuse   - number of profiles     retained for further processing
!     nouse   - number of observations retained for further processing
!     nobs    - array of observations on each subdomain for each processor
!
!$$$ end documentation block

! Declare module variables
  use kinds,     only: r_kind, i_kind
  use gridmod,   only: nlat, nlon, regional, tll2xy, rlats, rlons
  use constants, only: deg2rad, zero, r60inv, tiny_r_kind
  use gsi_4dvar, only: l4dvar, l4densvar, iwinbgn, winlen
  use netcdf,    only: nf90_open, nf90_nowrite, nf90_noerr, nf90_inq_dimid,    &
                       nf90_inquire_dimension, nf90_inq_varid, nf90_get_var,   &
                       nf90_close
  use mpimod,    only: npe
  use tgasinfo,  only: tgas_szamax,tgas_albmax,tgas_cldmax
  ! thinning
  use satthin,   only: makegrids,map2tgrid,destroygrids,finalcheck,itxmax  

  implicit none

! Declare passed variables
  character(len=*), intent(in   ) :: obstype, infile, jsatid
  character(len=*), intent(in   ) :: sis
  integer(i_kind),  intent(in   ) :: lunout, ithin
  real(r_kind),     intent(in   ) :: gstime, twind, rmesh
  integer(i_kind),  intent(inout) :: nread
  integer(i_kind),  intent(inout) :: npuse, nouse
  integer(i_kind), dimension(npe), intent(inout) :: nobs

! Declare local parameters
  logical,         parameter :: ldebug = .false.

  integer(i_kind), parameter :: ilat   = 4                   ! index of lat in tgasout
  integer(i_kind), parameter :: ilon   = 3                   ! index of lon in tgasout
  integer(i_kind), parameter :: maxobs = 16777216            ! greatest IEEE int w/ exact real

  real(r_kind),    parameter :: r90    =  90.0_r_kind
  real(r_kind),    parameter :: r360   = 360.0_r_kind

! Declare local variables
  logical         :: loutside, lskip, useak

  integer(i_kind) :: idate5(5)
  integer(i_kind) :: id_fin, id_nsound, id_nchanl, id_navg, id_nedge
  integer(i_kind) :: id_date, id_time, id_year, id_mon, id_day
  integer(i_kind) :: id_hour, id_mint, id_lat, id_lon
  integer(i_kind) :: id_prs, id_bad, id_unc, id_obs
  integer(i_kind) :: id_avgker, id_priorobs, id_priorpro
  integer(i_kind) :: id_mode, id_stype
  integer(i_kind) :: imode, istype
  integer(i_kind) :: nsound, nreal, ntgasdat, nchanl, navg, nedge
  integer(i_kind) :: nmind, j, k, n, iout, istat, ichan

  real(r_kind)    :: deglat, deglon, radlat, radlon, grdlat, grdlon
  real(r_kind)    :: grdtime

  character(len=10), allocatable :: qcflags10(:)

  integer(i_kind),   allocatable :: iyears(:), imons(:), idays(:)
  integer(i_kind),   allocatable :: ihours(:), imints(:), isecs(:)
  integer(i_kind),   allocatable :: idates(:), itimes(:), itimes6(:,:)
  integer(i_kind),   allocatable :: imodes(:), istypes(:)
  integer(i_kind),   allocatable :: isbads(:,:), isbads1(:)

  real(r_kind),      allocatable :: satlats(:), satlons(:)
  real(r_kind),      allocatable :: pchnom(:), pchanls(:,:), uncerts(:,:)
  real(r_kind),      allocatable :: psurfs(:), peavgs(:,:), avgkers(:,:,:)
  real(r_kind),      allocatable :: priorobses(:,:), priorpros(:,:)
  real(r_kind),      allocatable :: obses(:,:), tgasout(:,:)

  real(r_kind),      allocatable :: pchanls1(:), uncerts1(:)
  real(r_kind),      allocatable :: avgkers1(:,:), priorobses1(:), obses1(:)

  ! NO2/SO2 stuff
  logical                        :: isdoas
  integer(i_kind)                :: idx, ier, iflag
  integer(i_kind)                :: id_sza,  id_scwp,  id_albd, id_cldfrc, id_qav
  integer(i_kind)                :: id_vcds, id_tropp, id_amfs, id_rows
  real(r_kind)                   :: amf
  real(r_kind),      allocatable :: scwtpress(:)
  real(r_kind),      allocatable :: szas1(:), cldfrcs1(:), qflags1(:), albds1(:), qav1(:)
  real(r_kind),      allocatable :: rows(:), tropp(:), vcds(:), amfs(:)

  ! thinning
  real(r_kind)                   :: tdiff,sstime,crit1,dist1,timedif
  integer(i_kind)                :: i,kk,iidx,itx,itt,nodata
  logical                        :: iuse
  real(r_kind),parameter         :: r6   = 6.0_r_kind

! tgez = Point obs on a vertical grid of altitudes
! tgev = Point obs on a vertical grid of pressures
! tgaz = Avg kernel obs on a vertical grid of altitudes
! tgav = Avg kernel obs on a vertical grid of pressures
! tgop = ObsPack (specialized, may remove)

!$$$ -- 1. ALLOCATE AND INITIALIZE VARIABLES
!===============================================================================
! DOAS style NO2/SO2 obs?
  isdoas = (obstype == 'gomno2' .or. obstype == 'tomno2' .or.                   &
            obstype == 'omno2'  .or. obstype == 'omso2'  .or.                   &
            obstype == 'momno2' )

! Check if this is an averaging-kernel obs type
  useak = (obstype == 'tgav'   .or. obstype == 'tgaz'   .or.                   &
           obstype == 'acos'   .or. isdoas ) 

! Get dimension lengths
  if (obstype == 'tgez' .or. obstype == 'tgev') then
     call check(nf90_open(infile, nf90_nowrite,  id_fin))
     call check(nf90_inq_dimid(id_fin, 'nsound', id_nsound))
     call check(nf90_inquire_dimension(id_fin,   id_nsound, len=nsound))

!    Try to get number of channels, if not, set to 1
     istat = nf90_inq_dimid(id_fin, 'nchanl', id_nchanl)
     if (istat == nf90_noerr) then
        call check(nf90_inquire_dimension(id_fin, id_nchanl, len=nchanl))
     else
        nchanl = 1
     end if

     navg  = 1
     nreal = 6 + 3*nchanl

  else if (obstype == 'tgav' .or. obstype == 'tgaz') then
     call check(nf90_open(infile, nf90_nowrite,  id_fin))
     call check(nf90_inq_dimid(id_fin, 'nsound', id_nsound))
     call check(nf90_inquire_dimension(id_fin,   id_nsound, len=nsound))
     call check(nf90_inq_dimid(id_fin, 'navg',   id_navg))
     call check(nf90_inquire_dimension(id_fin,   id_navg,   len=navg))

!    Try to get number of channels, if not, set to 1
     istat = nf90_inq_dimid(id_fin, 'nchanl', id_nchanl)
     if (istat == nf90_noerr) then
        call check(nf90_inquire_dimension(id_fin, id_nchanl, len=nchanl))
     else
        nchanl = 1
     end if

!    Try to get number of edges, if not, set to navg + 1
     istat = nf90_inq_dimid(id_fin, 'nedge', id_nedge)
     if (istat == nf90_noerr) then
        call check(nf90_inquire_dimension(id_fin, id_nedge, len=nedge))
     else
        nedge = navg + 1
     end if

     nreal = 6 + 3*nchanl + nedge + (navg + navg*nchanl + nchanl)

  else if (obstype == 'tgop') then
     call check(nf90_open(infile, nf90_nowrite,  id_fin))
     call check(nf90_inq_dimid(id_fin, 'obs',    id_nsound))
     call check(nf90_inquire_dimension(id_fin,   id_nsound, len=nsound))

     nchanl = 1
     navg   = 1
     nreal  = 6 + 3*nchanl

  else if (obstype == 'acos') then
     call check(nf90_open(infile, nf90_nowrite, id_fin))
     call check(nf90_inq_dimid(id_fin, 'sounding_id', id_nsound))
     call check(nf90_inquire_dimension(id_fin,  id_nsound, len=nsound))

     nchanl = 1
     navg   = 20
     nedge  = navg
     nreal  = 6 + 3*nchanl + nedge + (navg + navg*nchanl + nchanl)

  ! OMI NO2 & SO2
  else if ( isdoas ) then

     ! make thinning grids
     call makegrids(rmesh,ithin)

     call check(nf90_open(infile, nf90_nowrite, id_fin))
     call check(nf90_inq_dimid(id_fin, 'nrec', id_nsound))
     call check(nf90_inquire_dimension(id_fin, id_nsound, len=nsound))
     call check(nf90_inq_dimid(id_fin, 'nlev', id_navg))
     call check(nf90_inquire_dimension(id_fin, id_navg, len=navg))

     nchanl = 1
     nedge  = navg + 1
     nreal  = 6 + 3*nchanl + nedge + (navg + navg*nchanl + nchanl)

  else
     return

  end if

! Make sure dimensions, etc. are ok, if eg the allocates below fail
  if (ldebug) then
     print *, 'obstype  = ', obstype, 'jsatid   = ', jsatid, 'sis      = ', sis
     print *, '---'
     print *, 'nsound   = ', nsound
     print *, 'nreal    = ', nreal
     print *, 'nchanl   = ', nchanl
     print *, 'ntgasdat = ', ntgasdat
     if (useak) print *, 'navg     = ', navg, 'nedge    = ', nedge
  end if

! Allocate local arrays
  allocate(idates(nsound), itimes(nsound), iyears(nsound), imons(nsound),      &
           idays(nsound),  ihours(nsound), imints(nsound), isecs(nsound))
  allocate(satlats(nsound), satlons(nsound))
  allocate(pchnom(nchanl), pchanls(nchanl,nsound), isbads(nchanl,nsound),      &
           uncerts(nchanl,nsound),  obses(nchanl,nsound))

! Fill with negative so we don't accidentally misuse
  pchanls = -999.

  if (useak) allocate(psurfs(nsound), peavgs(nedge,nsound),                    &
                      avgkers(nchanl,navg,nsound), priorobses(nchanl,nsound),  &
                      priorpros(navg,nsound))

  if (obstype == 'acos') allocate(imodes(nsound), istypes(nsound))

  if (nchanl  == 1) allocate(pchanls1(nsound),    isbads1(nsound),             &
                             uncerts1(nsound),    avgkers1(navg,nsound),       &
                             priorobses1(nsound), obses1(nsound))

  if ( isdoas ) then 
     allocate(szas1(nsound), cldfrcs1(nsound), qflags1(nsound), albds1(nsound), qav1(nsound) )
     allocate(scwtpress(navg), rows(nsound), tropp(nsound), vcds(nsound), amfs(nsound))
  endif

! Allocate output array
  ntgasdat = nreal + nchanl
  if ( isdoas ) then
     allocate(tgasout(ntgasdat,itxmax))
  else
     allocate(tgasout(ntgasdat,maxobs))
  endif

!$$$ -- 2. READ DATA
!$$$       POINT (TGEZ,TGEV) & AVG KER (TGAV,TGAZ) TEMPLATES
!          ON ALTITUDE (TG*Z) & PRESSURE (TG*V) GRIDS
!===============================================================================
  if (obstype == 'tgez' .or. obstype == 'tgev' .or. obstype == 'tgav' .or.      &
      obstype == 'tgaz') then
!    Required variables standard across all types
     call check(nf90_inq_varid(id_fin, 'date', id_date))
     call check(nf90_inq_varid(id_fin, 'time', id_time))
     call check(nf90_inq_varid(id_fin, 'lat',  id_lat))
     call check(nf90_inq_varid(id_fin, 'lon',  id_lon))

     call check(nf90_get_var(id_fin, id_date, idates))
     call check(nf90_get_var(id_fin, id_time, itimes))
     call check(nf90_get_var(id_fin, id_lat,  satlats))
     call check(nf90_get_var(id_fin, id_lon,  satlons))

!    Translate dates and times to years, months, etc.
     iyears = idates/10000
     imons  = idates/100 - iyears*100
     idays  = idates - imons*100  - iyears*10000
     ihours = itimes/10000
     imints = itimes/100 - ihours*100
     isecs  = itimes - imints*100 - ihours*10000

!    A variable defining the vertical grid must exist.  It depends on obstype,
!    can be fixed, can be altitudes (tgez) or pressures (tgev), and can
!    support averaging kernels (tgaz/tgav).
     if (obstype == 'tgez') then
        call check(nf90_inq_varid(id_fin, 'zchanl', id_prs))
!       Assume it is the nominal grid, if not, try to read as per sounding
        istat = nf90_get_var(id_fin, id_prs, pchnom)
        if (istat == nf90_noerr) then
           do n = 1,nsound
              pchanls(:,n) = pchnom
           end do
        else
           call check(nf90_get_var(id_fin, id_prs, pchanls))
        end if

     else if (obstype == 'tgev') then
        call check(nf90_inq_varid(id_fin, 'pchanl', id_prs))
!       Assume it is the nominal grid, if not, try to read as per sounding
        istat = nf90_get_var(id_fin, id_prs, pchnom)
        if (istat == nf90_noerr) then
           do n = 1,nsound
              pchanls(:,n) = pchnom
           end do
        else
           call check(nf90_get_var(id_fin, id_prs, pchanls))
        end if

     else if (obstype == 'tgav') then
        call check(nf90_inq_varid(id_fin, 'peavg', id_prs))
        call check(nf90_get_var(id_fin, id_prs, peavgs))

     else if (obstype == 'tgaz') then
        call check(nf90_inq_varid(id_fin, 'zeavg', id_prs))
        call check(nf90_get_var(id_fin, id_prs, peavgs))
     end if

!    Optional prior profile, always has same form
     priorpros = 0.
     istat = nf90_inq_varid(id_fin, 'priorpro', id_priorpro)
     if (istat == nf90_noerr) call check(nf90_get_var(id_fin, id_priorpro, priorpros))

!    Support optional nchanl dimension, not required for column retrievals
     ichan = nf90_inq_dimid(id_fin, 'nchanl', id_nchanl)
     if (ichan == nf90_noerr) then
!       Fill optional variables with nominal values
        isbads  = 0
        uncerts = 1.
        avgkers    = 1.
        priorobses = 0.

        istat = nf90_inq_varid(id_fin, 'isbad',  id_bad)
        if (istat == nf90_noerr) call check(nf90_get_var(id_fin, id_bad, isbads))

        istat = nf90_inq_varid(id_fin, 'uncert', id_unc)
        if (istat == nf90_noerr) call check(nf90_get_var(id_fin, id_unc, uncerts))

        istat = nf90_inq_varid(id_fin, 'avgker', id_avgker)
        if (istat == nf90_noerr) call check(nf90_get_var(id_fin, id_avgker, avgkers))

        istat = nf90_inq_varid(id_fin, 'priorobs', id_priorobs)
        if (istat == nf90_noerr) call check(nf90_get_var(id_fin, id_priorobs, priorobses))

        call check(nf90_inq_varid(id_fin, 'obs', id_obs))
        call check(nf90_get_var(id_fin, id_obs, obses))
     else
!       Fill optional variables with nominal values
        isbads1  = 0
        uncerts1 = 1.
        avgkers1    = 1.
        priorobses1 = 0.

        istat = nf90_inq_varid(id_fin, 'isbad',  id_bad)
        if (istat == nf90_noerr) call check(nf90_get_var(id_fin, id_bad, isbads1))

        istat = nf90_inq_varid(id_fin, 'uncert', id_unc)
        if (istat == nf90_noerr) call check(nf90_get_var(id_fin, id_unc, uncerts1))

        istat = nf90_inq_varid(id_fin, 'avgker', id_avgker)
        if (istat == nf90_noerr) call check(nf90_get_var(id_fin, id_avgker, avgkers1))

        istat = nf90_inq_varid(id_fin, 'priorobs', id_priorobs)
        if (istat == nf90_noerr) call check(nf90_get_var(id_fin, id_priorobs, priorobses1))

        call check(nf90_inq_varid(id_fin, 'obs', id_obs))
        call check(nf90_get_var(id_fin, id_obs, obses1))

        isbads(1,:)  = isbads1
        uncerts(1,:) = uncerts1
        avgkers(1,:,:)  = avgkers1
        priorobses(1,:) = priorobses1
        obses(1,:)      = obses1
     end if

!$$$       TGOP: DATA FROM NOAA/ESRL OBSPACK (MAY BE TRANSITIONED TO TGEZ)
!===============================================================================
  else if (obstype == 'tgop') then
     allocate(itimes6(6,nsound), qcflags10(nsound))

     call check(nf90_inq_varid(id_fin, 'time_components', id_time))
     call check(nf90_inq_varid(id_fin, 'latitude',        id_lat))
     call check(nf90_inq_varid(id_fin, 'longitude',       id_lon))
     call check(nf90_inq_varid(id_fin, 'value',           id_obs))

     call check(nf90_get_var(id_fin, id_time, itimes6))
     call check(nf90_get_var(id_fin, id_lat,  satlats))
     call check(nf90_get_var(id_fin, id_lon,  satlons))
     call check(nf90_get_var(id_fin, id_obs,  obses1))

!    Some datasets incorrectly report altitude as gps_altitude
     istat = nf90_inq_varid(id_fin, 'altitude', id_prs)
     if (istat /= nf90_noerr) then
        call check(nf90_inq_varid(id_fin, 'gps_altitude', id_prs))
     end if
!    Hack to hold altitude in the pressure variable
     call check(nf90_get_var(id_fin, id_prs, pchanls1))

!    All  datasets report an integer   obs_flag, which is sometimes useful
!    Some datasets report a  character qcflag,   which is always    useful
!    Some datasets report an integer   qcflag,   which is maybe?    useful
!    For now, we only use the first two fields of the character qcflag since
!    the other flags aren't consistently reported
     isbads = 0
     istat  = nf90_inq_varid(id_fin, 'qcflag', id_bad)
     if (istat == nf90_noerr) then
!       Try to read it as a character array
        istat = nf90_get_var(id_fin, id_bad, qcflags10)
        if (istat == nf90_noerr) then
           do n = 1,nsound
              if (qcflags10(n)(1:2) /= '..') isbads(1,n) = 1
           end do
        end if
     end if
!    If the dataset has a CT_MDM field, use that over everything else
     istat = nf90_inq_varid(id_fin, 'CT_ASSIM', id_unc)
     if (istat == nf90_noerr) then
        call check(nf90_get_var(id_fin, id_unc, isbads))
        isbads = max(1 - isbads, 0)
     end if

!    Some datasets report an uncertainty, but very many don't; for the ones
!    that don't, we make the very bad assumption that it is 10^-4 times the
!    obs value (FIXME)
     uncerts1 = 1.e-4_r_kind * obses1
     istat    = nf90_inq_varid(id_fin, 'value_unc', id_unc)
     if (istat == nf90_noerr) then
        call check(nf90_get_var(id_fin, id_unc, uncerts1))
     end if
!    If the dataset has a CT_MDM field, use that over everything else
     istat = nf90_inq_varid(id_fin, 'CT_MDM', id_unc)
     if (istat == nf90_noerr) then
        call check(nf90_get_var(id_fin, id_unc, uncerts1))
     end if

     pchanls(1,:) = pchanls1
     uncerts(1,:) = uncerts1
     obses(1,:)   = obses1

!    Translate dates and times to years, months, etc.
     iyears = itimes6(1,:)
     imons  = itimes6(2,:)
     idays  = itimes6(3,:)
     ihours = itimes6(4,:)
     imints = itimes6(5,:)
     isecs  = itimes6(6,:)

     deallocate(itimes6, qcflags10)

!$$$       ACOS GOSAT/OCO/GeoCarb (WILL BE TRANSITIONED TO TGAV)
!===============================================================================
  else if (obstype == 'acos') then
     call check(nf90_inq_varid(id_fin, 'sounding_date',       id_date))
     call check(nf90_inq_varid(id_fin, 'sounding_time',       id_time))
     call check(nf90_inq_varid(id_fin, 'latitude',            id_lat))
     call check(nf90_inq_varid(id_fin, 'longitude',           id_lon))
     call check(nf90_inq_varid(id_fin, 'psurf',               id_prs))
     call check(nf90_inq_varid(id_fin, 'qcflag',              id_bad))
     call check(nf90_inq_varid(id_fin, 'xco2_uncert',         id_unc))
     call check(nf90_inq_varid(id_fin, 'xco2_final',          id_obs))
     call check(nf90_inq_varid(id_fin, 'co2_profile_apriori', id_priorpro))
     call check(nf90_inq_varid(id_fin, 'xco2_avgker',         id_avgker))
     call check(nf90_inq_varid(id_fin, 'xco2_apriori',        id_priorobs))

     call check(nf90_inq_varid(id_fin, 'operation_mode',      id_mode))
     call check(nf90_inq_varid(id_fin, 'surface_type',        id_stype))

     call check(nf90_get_var(id_fin, id_date,     idates))
     call check(nf90_get_var(id_fin, id_time,     itimes))
     call check(nf90_get_var(id_fin, id_lat,      satlats))
     call check(nf90_get_var(id_fin, id_lon,      satlons))
     call check(nf90_get_var(id_fin, id_prs,      psurfs))
     call check(nf90_get_var(id_fin, id_bad,      isbads1))
     call check(nf90_get_var(id_fin, id_unc,      uncerts1))
     call check(nf90_get_var(id_fin, id_obs,      obses1))
     call check(nf90_get_var(id_fin, id_avgker,   avgkers1))
     call check(nf90_get_var(id_fin, id_priorpro, priorpros))
     call check(nf90_get_var(id_fin, id_priorobs, priorobses1))

     call check(nf90_get_var(id_fin, id_mode,     imodes))
     call check(nf90_get_var(id_fin, id_stype,    istypes))

!    Translate dates and times to years, months, etc.
     iyears = idates/10000
     imons  = idates/100 - iyears*100
     idays  = idates - imons*100  - iyears*10000
     ihours = itimes/10000
     imints = itimes/100 - ihours*100
     isecs  = itimes - imints*100 - ihours*10000

!    Construct obs edge pressures
     peavgs(1,:) = psurfs*1.e-4_r_kind
     do k = 2,nedge
        peavgs(k,:) = psurfs - real(nedge-k, r_kind)/19._r_kind*psurfs
     end do

!    Hack to deal with singleton leading dimension
     pchanls(1,:) = psurfs
     isbads(1,:)  = isbads1
     uncerts(1,:) = uncerts1

     avgkers(1,:,:)  = avgkers1
     priorobses(1,:) = priorobses1
     obses(1,:)      = obses1

!$$$       NO2 (MINDS OMI/GOME/TROPOMI and OMI SO2 
!===============================================================================
  else if ( isdoas ) then
     ! omi specific arrays

     call check(nf90_inq_varid(id_fin, 'Year',                          id_year))
     call check(nf90_inq_varid(id_fin, 'Month',                         id_mon))
     call check(nf90_inq_varid(id_fin, 'Day',                           id_day))
     call check(nf90_inq_varid(id_fin, 'Hour',                          id_hour))
     call check(nf90_inq_varid(id_fin, 'Minute',                        id_mint))
     call check(nf90_inq_varid(id_fin, 'Row',                           id_rows))
     call check(nf90_inq_varid(id_fin, 'Latitude',                      id_lat))
     call check(nf90_inq_varid(id_fin, 'Longitude',                     id_lon))
     call check(nf90_inq_varid(id_fin, 'SolarZenithAngle',              id_sza))
     call check(nf90_inq_varid(id_fin, 'TerrainPressure',               id_prs))
     call check(nf90_inq_varid(id_fin, 'ScatteringWeight',              id_avgker))
     call check(nf90_inq_varid(id_fin, 'CloudRadianceFraction',         id_cldfrc))
     ! should change to:
     !call check(nf90_inq_varid(id_fin, 'CloudFraction',                 id_cldfrc))
     if ( obstype == 'omso2' ) then
        call check(nf90_inq_varid(id_fin, 'SlantColumnAmountSO2',          id_obs))
        call check(nf90_inq_varid(id_fin, 'AlgorithmFlag_SnowIce',         id_unc))
        call check(nf90_inq_varid(id_fin, 'LayerBottomPressure',           id_scwp))
        call check(nf90_inq_varid(id_fin, 'SurfaceReflectivity',           id_albd))
        call check(nf90_inq_varid(id_fin, 'Flag_RowAnomaly',               id_bad))
        call check(nf90_inq_varid(id_fin, 'Flag_SAA',                      id_tropp))
     else
        if ( obstype == 'omno2' ) then
           call check(nf90_inq_varid(id_fin, 'SlantColumnAmountNO2Destriped', id_obs))
           call check(nf90_inq_varid(id_fin, 'ScatteringWtPressure',          id_scwp))
        else
           call check(nf90_inq_varid(id_fin, 'SlantColumnAmountNO2',          id_obs))
           call check(nf90_inq_varid(id_fin, 'ScatteringWeightPressure',      id_scwp))
        endif
        call check(nf90_inq_varid(id_fin, 'SlantColumnAmountNO2Std',       id_unc))
        call check(nf90_inq_varid(id_fin, 'TerrainReflectivity',           id_albd))
        call check(nf90_inq_varid(id_fin, 'VcdQualityFlags',               id_bad))
        call check(nf90_inq_varid(id_fin, 'TropopausePressure',            id_tropp))
!        call check(nf90_inq_varid(id_fin, 'ColumnAmountNO2Strat',          id_vcds))
!        call check(nf90_inq_varid(id_fin, 'AmfStrat',                      id_amfs))
        if ( obstype == 'tomno2' ) then
           call check(nf90_inq_varid(id_fin, 'qa_value', id_qav))
        endif
     endif

     call check(nf90_get_var(id_fin, id_year,   iyears))
     call check(nf90_get_var(id_fin, id_mon,    imons))
     call check(nf90_get_var(id_fin, id_day,    idays))
     call check(nf90_get_var(id_fin, id_hour,   ihours))
     call check(nf90_get_var(id_fin, id_mint,   imints))
     call check(nf90_get_var(id_fin, id_rows,   rows  ))
     call check(nf90_get_var(id_fin, id_lat,    satlats))
     call check(nf90_get_var(id_fin, id_lon,    satlons))
     call check(nf90_get_var(id_fin, id_obs,    obses1))
     call check(nf90_get_var(id_fin, id_unc,    uncerts1))
     call check(nf90_get_var(id_fin, id_prs,    psurfs))
     call check(nf90_get_var(id_fin, id_avgker, avgkers1))
     call check(nf90_get_var(id_fin, id_scwp,   scwtpress))
     call check(nf90_get_var(id_fin, id_sza,    szas1))
     call check(nf90_get_var(id_fin, id_albd,   albds1))
     call check(nf90_get_var(id_fin, id_cldfrc, cldfrcs1))
     call check(nf90_get_var(id_fin, id_bad,    qflags1))
     call check(nf90_get_var(id_fin, id_tropp,  tropp))
     if ( obstype == 'tomno2' ) then
        call check(nf90_get_var(id_fin, id_qav,  qav1))
     endif
     ! For OMNO2, albedo and cloud fraction are integers (0-1000). Convert to fractions here
     if ( obstype == 'omno2' ) then
        albds1(:)   = albds1(:)   / 1000.0
        cldfrcs1(:) = cldfrcs1(:) / 1000.0
     endif

!     call check(nf90_get_var(id_fin, id_vcds,   vcds))
!     call check(nf90_get_var(id_fin, id_amfs,   amfs))

     ! Hack for compatability
     isecs = 0._r_kind

     ! Set observation pressure to 1000.0
     pchanls1(:) = 1000.0
     pchanls(1,:) = pchanls1

     ! obs edge pressures. this is the same for all obs, get from data
     ! and broadcast to all observations.
     ! omi pressures are in hPa. 
     ! set surface pressure to really high value to make sure that we will
     ! cover the entire GEOS profile (will be cropped in setuptgas).
     peavgs(1,:) = 1200_r_kind
     ! SO2: bottom layer pressure. OMI SO2 vertical axis is GEOS-style, i.e., index 1 is top of 
     ! atmosphere
     if ( obstype == 'omso2' ) then
        do k = 2,navg
           peavgs(k,:) = scwtpress(navg-k+1)
        end do
     ! NO2: middle of layer
     else
        do k = 2,nedge-1
           peavgs(k,:) = 0.5_r_kind*(scwtpress(k-1)+scwtpress(k))  ! hPa 
           peavgs(k,:) = min(peavgs(k,:), peavgs(1,:))
        end do
     endif
     ! set top level to max GEOS level 
     peavgs(nedge,:) = 0.01_r_kind                              ! hPa

     ! set data quality to 'good' for every data point. Bad data 
     ! will be skipped in the writing step below (when calculating AMF, etc.)
     isbads(:,:) = 0
  end if

  call check(nf90_close(id_fin))


!$$$ -- 3. FILL & WRITE OUTPUT ARRAY
!===============================================================================
  nsound = min(nsound, maxobs)
  npuse  = 0
  do n = 1,nsound
!    Convert observation time to relative time
     idate5(1) = iyears(n)
     idate5(2) = imons(n)
     idate5(3) = idays(n)
     idate5(4) = ihours(n)
     idate5(5) = imints(n)

!    Get obs time (nmind) in minutes relative to historic date
     call w3fs21(idate5, nmind)

!    Compute grid relative time in hours (grdtime): subtract off starting time
!    in reference minutes (iwinbgn), add seconds back in, and rescale to hours
!    grdtime = (real(nmind - iwinbgn, r_kind) + isecs(n)*r60inv)*r60inv
     grdtime = real(nmind - iwinbgn, r_kind)*r60inv

!    Toss obs outside time window
     if (l4dvar .or. l4densvar) then
        if (grdtime < zero .or.          winlen <= grdtime) cycle
     else
        if (grdtime < zero .or. 2._r_kind*twind <= grdtime) cycle
     end if

!    Convert observation lat/lon to relative position
     deglat = satlats(n)
     deglon = satlons(n)
     if (deglon < zero) deglon = deglon + r360

!    Toss out any obs with suspicious locations
     if (r90 < abs(deglat) .or. deglon < zero .or. r360 < deglon) cycle

     radlat = deglat * deg2rad
     radlon = deglon * deg2rad
    
     if (regional) then
        call tll2xy(radlon, radlat, grdlon, grdlat, loutside)
        if (loutside) cycle 
     else
        grdlat = radlat
        grdlon = radlon
        call grdcrd1(grdlat, rlats, nlat, 1)
        call grdcrd1(grdlon, rlons, nlon, 1)
     end if

!    Optionally match observing mode and surface type for acos obstype
!    ---
!    Only does anything if sis specifies a valid retrieval mode.  In that
!    case, only uses retrievals with the appropriate mode.  Otherwise, all
!    modes are processed together.
     if (obstype == 'acos') then
       lskip = .false.

       if (index(jsatid,'oco') /= 0) then
          imode  = imodes(n)
          istype = istypes(n)

          if      (index(sis,'lndnd')  /= 0) then
            if (imode /= 0 .or. istype /= 1) lskip = .true.
          else if (index(sis,'lndgl')  /= 0) then
            if (imode /= 1 .or. istype /= 1) lskip = .true.
          else if (index(sis,'ocngl')  /= 0) then
            if (imode /= 1 .or. istype /= 0) lskip = .true.
          else if (index(sis,'targ')   /= 0) then
            if (imode /= 2)                  lskip = .true.
          else if (index(sis,'trans')  /= 0) then
            if (imode /= 3)                  lskip = .true.
          else if (index(sis,'sam')    /= 0) then
            if (imode /= 4)                  lskip = .true.
          end if

       else if (index(jsatid,'gosat') /= 0) then
          imode  = imodes(n)
          istype = istypes(n)

          if      (index(sis,'lndmg') /= 0) then
            if (imode /= 0 .or. istype /= 1) lskip = .true.
          else if (index(sis,'lndhg') /= 0) then
            if (imode /= 1 .or. istype /= 1) lskip = .true.
          else if (index(sis,'ocnhg') /= 0) then
            if (imode /= 1 .or. istype /= 0) lskip = .true.
          end if
       end if

       if (lskip) cycle
     end if

     ! DOAS style observations
     if ( isdoas ) then

        ! check data quality: 
        if ( obstype == 'omso2' ) then
           lskip = .false.
           if ( rows(n)<5. .or. rows(n)>55. ) lskip = .true.  ! rows 5-55 (1-based)
           if ( uncerts1(n) /= 0.0          ) lskip = .true.  ! non-zero SnowIce flag
           if ( qflags1(n)  /= 0.0          ) lskip = .true.  ! RowAnommaly flag
           if ( tropp(n)    /= 0.0          ) lskip = .true.  ! SAA flag
           if ( szas1(n)    > tgas_szamax(2)) lskip = .true.  ! SZA > value set in GSI_GridComp.rc 
           if ( albds1(n)   > tgas_albmax(2)) lskip = .true.  ! Surface reflectivity > 0.3
           if ( cldfrcs1(n) > tgas_cldmax(2)) lskip = .true.  ! Cloud radiance fraction > 0.3  
           !if ( albds1(n)   > 0.3           ) lskip = .true.  ! Surface reflectivity > 0.3
           !if ( cldfrcs1(n) > 0.3           ) lskip = .true.  ! Cloud radiance fraction > 0.3  
        else
           ! - VcdQualityFlag is an even integer 
           ! - sza < 80
           ! - reflectivity (albedo) < 0.3
           ! - cloud radiance fraction < 0.5
           lskip = .true.
           iflag = int(qflags1(n))
           if ( IAND ( iflag, 1 ) .eq. 0    ) lskip = .false.
           if ( rows(n)<5. .or. rows(n)>55. ) lskip = .true.  ! rows 5-55 (1-based)
           if ( szas1(n)    > tgas_szamax(1)) lskip = .true.
           if ( albds1(n)   > tgas_albmax(1)) lskip = .true.
           if ( cldfrcs1(n) > tgas_cldmax(1)) lskip = .true.
           !if ( albds1(n)   > 300.          ) lskip = .true.
           !if ( cldfrcs1(n) > 500.          ) lskip = .true.
           ! for TROPOMI, skip values with a qa_value <= 0.75
           if ( obstype == 'tomno2' ) then
              if ( qav1(n) <= 0.75 ) lskip = .true.
           endif

        endif
        if (lskip) cycle

        ! thin data
        if (l4dvar .or. l4densvar) then
           timedif = zero
        else
           sstime=real(nmind,r_kind)
           tdiff=(sstime-gstime)*r60inv
           timedif = r6*abs(tdiff)
        endif
        crit1 = 0.01_r_kind+timedif
        call map2tgrid(radlat,radlon,dist1,crit1,itx,ithin,itt,iuse,sis)
        if ( .not. iuse ) cycle
        call finalcheck(dist1,crit1,itx,iuse)
        if ( .not. iuse ) cycle

        ! Convert obs to 1e15 molec cm-2
        obses(1,n)      = obses1(n) / 1.0e15
        if ( obstype == 'omso2' ) then
           ! flip averaging kernel to make sure that surface is index 1
           avgkers(1,:,n)  = avgkers1(navg:1:-1,n)
        else
           avgkers(1,:,n)  = avgkers1(:,n)
        endif

        ! dummy a-priori profile 
        priorpros(:,n) = -99.0

        ! Uncertainties
        if ( obstype == 'omso2' ) then
           ! The OMI SO2 readme document reports SCD uncertainties of ~0.2 DU
           ! if SZA < 50deg, and 0.3-0.4 DU for SZAs 50-70degrees. Translate to 1.0e15 molec/cm2 here.
           if ( szas1(n) < 50.0 ) then
              uncerts(1,n)  = 5.374   ! 0.2 DU = 5.374e15 molec/cm2 
           else
              uncerts(1,n)  = 10.748  ! 0.4 DU = 10.748e15 molec/cm2
           endif
        else
           ! Pass SCD uncertainties. Will derive VCD uncertainties from this 
           uncerts(1,n)   = uncerts1(n) / 1.0e15
        endif

     endif

!    Write record to output file
     npuse = npuse + 1
     nodata = npuse

     ! index in tgasout to write into
     iidx = npuse
     if ( isdoas ) iidx = itx

     tgasout(1,iidx) = n
     tgasout(2,iidx) = grdtime
     tgasout(3,iidx) = grdlon
     tgasout(4,iidx) = grdlat
     tgasout(5,iidx) = deglon
     tgasout(6,iidx) = deglat

     iout = 6
     do k = 1,nchanl
        iout = iout + 1
        tgasout(iout,iidx) = pchanls(k,n)             ! channel pressure
     end do
     do k = 1,nchanl
        iout = iout + 1
        tgasout(iout,iidx) = isbads(k,n)              ! qc flag
     end do
     do k = 1,nchanl
        iout = iout + 1
        tgasout(iout,iidx) = uncerts(k,n)             ! channel uncertainty
     end do 

     if (useak) then
        do k = 1,nedge
           iout = iout + 1
           tgasout(iout,iidx) = peavgs(k,n)           ! averaging kernel edge pressures
        end do
        do k = 1,navg
           do j = 1,nchanl
              iout = iout + 1
              tgasout(iout,iidx) = avgkers(j,k,n)     ! averaging kernel
           end do
        end do
        do k = 1,navg
           iout = iout + 1
           tgasout(iout,iidx) = priorpros(k,n)        ! a priori profile
        end do 
        do j = 1,nchanl
           iout = iout + 1
           tgasout(iout,iidx) = priorobses(j,n)       ! a priori obs data
        end do
     end if

     do k = 1,nchanl
        iout = iout + 1
        tgasout(iout,iidx) = obses(k,n)               ! observation
     end do

!    Print diagnostic output if debugging
     if (ldebug) then
        print *,    '---', n, '---'
        print *,    'date     = ', idate5(1), idate5(2), idate5(3)
        print *,    'time     = ', idate5(4), idate5(5)
        print *,    'grdtime  = ', grdtime
        print *,    'lat      = ', satlats(n), 'lon      = ', satlons(n)
        print *,    'grdlat   = ', grdlat,     'grdlon   = ', grdlon

        print *,    'pchanl   = '
        print *,    pchanls(:,n)
        print *,    'isbad    = '
        print *,    isbads(:,n)
        print *,    'uncert   = '
        print *,    uncerts(:,n)

        if (useak) then
           print *, 'peavg    = '
           print *, peavgs(:,n)
           print *, 'avgker   = '
           print *, avgkers(:,:,n)
           print *, 'priorpro = '
           print *, priorpros(:,n)
           print *, 'priorobs = '
           print *, priorobses(:,n)
        end if

        print *,    'obs      = '
        print *,    obses(:,n)
     end if
  end do

  nread = nsound * nchanl
  nouse = npuse  * nchanl

! If thinning, compress array to thinned data
  if (isdoas) then
     kk=0
     do k=1,itxmax
        if (tgasout(1,k)>zero) then
           kk=kk+1
           do i=1,ntgasdat
              tgasout(i,kk)=tgasout(i,k)
           end do
        endif
     end do
     npuse=kk
     nouse=npuse
  endif

  if (ldebug) then
     print *, '---'
     print *, 'nread    = ', nread
     print *, 'npuse    = ', npuse
     print *, 'nouse    = ', nouse
  end if

! Write header record and data to output file for further processing
! ***FIXME*** count_obs is a subroutine from obs_para that is neither external
! nor contained in a module
  call count_obs(npuse, ntgasdat, ilat, ilon, tgasout, nobs)
  write(lunout) obstype, sis, nreal, nchanl, ilat, ilon
  write(lunout) ((tgasout(j,k),j=1,ntgasdat),k=1,npuse)

! Deallocate output array
  deallocate(tgasout)

! Deallocate local arrays
  deallocate(idates, itimes, iyears, imons, idays, ihours, imints, isecs)
  deallocate(satlats, satlons, pchnom, pchanls, isbads, uncerts, obses)

  if (useak) deallocate(psurfs, peavgs, avgkers, priorobses, priorpros)

  if (obstype == 'acos') deallocate(imodes, istypes)
  if (nchanl  == 1)      deallocate(pchanls1, isbads1, uncerts1, avgkers1,     &
                                    priorobses1, obses1)
  if ( isdoas ) then 
     deallocate(szas1,cldfrcs1,qflags1,albds1,scwtpress,rows,tropp,vcds,amfs,qav1)
     call destroygrids
  endif

  return
  
contains
  subroutine check(istat)
    use mpeu_util, only: perr, die
    use netcdf,    only: nf90_strerror

    integer, intent(in) :: istat

    character(len=*), parameter :: myname_ = 'READ_TGAS.CHECK'

    if (istat /= nf90_noerr) then 
       call perr(myname_, 'nf90_strerror = ', trim(nf90_strerror(istat)))
       call die(myname_)
    end if
  end subroutine check 
end subroutine read_tgas
