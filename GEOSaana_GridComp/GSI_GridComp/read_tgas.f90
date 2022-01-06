subroutine read_tgas(nread, ndata, nodata, jsatid, infile, gstime, lunout,     &
                     obstype, twind, sis, ithin, rmesh)
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
!   2022-01-06   weir  - whew. been awhile. touching to change timestamp
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
!     ndata   - number of profiles retained for further processing
!     nodata  - number of observations retained for further processing
!
!$$$ end documentation block

! Declare module variables
  use kinds,     only: r_kind, i_kind
  use gridmod,   only: nlat, nlon, regional, tll2xy, rlats, rlons
  use constants, only: deg2rad, zero, r60inv
  use obsmod,    only: nlmopitt, nlacos, nlflask
  use gsi_4dvar, only: iwinbgn
  use netcdf,    only: nf90_open, nf90_nowrite, nf90_noerr, nf90_strerror,     &
                       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, &
                       nf90_get_var, nf90_close

  implicit none

! Declare passed variables
  character(len=*), intent(in   ) :: obstype, infile, jsatid
  character(len=*), intent(in   ) :: sis
  integer(i_kind),  intent(in   ) :: lunout, ithin
  integer(i_kind),  intent(inout) :: nread
  integer(i_kind),  intent(inout) :: ndata, nodata
  real(r_kind),     intent(in   ) :: gstime, twind, rmesh

! Declare local parameters
  real(r_kind),    parameter :: r90    =  90.0_r_kind
  real(r_kind),    parameter :: r360   = 360.0_r_kind

  integer(i_kind), parameter :: ilat   = 4
  integer(i_kind), parameter :: ilon   = 3
  integer(i_kind), parameter :: maxobs = 1e6

  logical,         parameter :: ldebug = .false.

! Declare local variables
  real(r_kind)    :: deglat, deglon, radlat, radlon, grdlat, grdlon, grdtime
  integer(i_kind) :: idate5(5)
  integer(i_kind) :: id_fin, id_nrec, id_date, id_time, id_year, id_mon
  integer(i_kind) :: id_day, id_hour, id_minu, id_lat, id_lon, id_qcflag
  integer(i_kind) :: id_psurf, id_avgker, id_priorobs, id_priorpro, id_uncert
  integer(i_kind) :: id_obs
  integer(i_kind) :: npro, nchanl, nreal, ntgasdat, nrec, nele, nmind, j, k
  logical         :: loutside

  integer(i_kind), allocatable :: idates(:), itimes(:), iyears(:), imons(:)
  integer(i_kind), allocatable :: idays(:), ihours(:), iminus(:)
  real(r_kind),    allocatable :: satlats(:), satlons(:), qcflags(:), psurfs(:)

  real(r_kind),    allocatable :: avgkers1(:,:), avgkers(:,:,:)
  real(r_kind),    allocatable :: priorobses1(:), priorobses(:,:)
  real(r_kind),    allocatable :: priorpros(:,:)
  real(r_kind),    allocatable :: uncerts1(:), uncerts(:,:)
  real(r_kind),    allocatable :: obses1(:), obses(:,:), tgasout(:,:)

! Get dimension lengths and allocate arrays
  if (obstype == 'acos') then
     npro   = nlacos
     nchanl = 1
     nreal  = 8 + 2*npro + 2

     call check(nf90_open(infile, nf90_nowrite, id_fin))

     call check(nf90_inq_dimid(id_fin, 'soundings', id_nrec))
     call check(nf90_inquire_dimension(id_fin, id_nrec, len=nrec))

     allocate(idates(nrec), itimes(nrec), satlats(nrec), satlons(nrec),        &
              qcflags(nrec), psurfs(nrec), priorpros(npro,nrec),               &
              avgkers1(npro,nrec), priorobses1(nrec), uncerts1(nrec),          &
              obses1(nrec))

  else if (obstype == 'mopitt') then
     npro   = nlmopitt
     nchanl = npro
     nreal  = 8 + 2*npro + npro**2

     call check(nf90_open(infile, nf90_nowrite, id_fin))

     call check(nf90_inq_dimid(id_fin, 'TimeCount', id_nrec))
     call check(nf90_inquire_dimension(id_fin, id_nrec, len=nrec))

     allocate(iyears(nrec), imons(nrec), idays(nrec), ihours(nrec),            &
              iminus(nrec), satlats(nrec), satlons(nrec), psurfs(nrec),        &
              priorpros(npro,nrec), avgkers(nchanl,npro,nrec),                 &
              uncerts(nchanl,nrec), obses(nchanl,nrec))

  else if (obstype == 'flask') then
     npro   = nlflask
     nchanl = 2
     nreal  = 8 + 2*npro

     call check(nf90_open(infile, nf90_nowrite, id_fin))

     call check(nf90_inq_dimid(id_fin, 'measurements', id_nrec))
     call check(nf90_inquire_dimension(id_fin, id_nrec, len=nrec))

     allocate(idates(nrec), itimes(nrec), satlats(nrec), satlons(nrec),        &
              uncerts(nchanl,nrec), obses(nchanl,nrec))

  end if

! Allocate output array
  ntgasdat = nreal + nchanl
  nele     = min(nrec, maxobs)
  allocate(tgasout(ntgasdat,nele))

  if (ldebug) then
     print *, 'nrec     = ', nrec
     print *, 'nreal    = ', nreal
     print *, 'nchanl   = ', nchanl
     print *, 'ntgasdat = ', ntgasdat
     print *, 'nele     = ', nele
  end if

! Read data
  if (obstype == 'acos') then
     call check(nf90_inq_varid(id_fin, 'sounding_date',       id_date))
     call check(nf90_inq_varid(id_fin, 'sounding_time',       id_time))
     call check(nf90_inq_varid(id_fin, 'latitude',            id_lat))
     call check(nf90_inq_varid(id_fin, 'longitude',           id_lon))
     call check(nf90_inq_varid(id_fin, 'qcflag',              id_qcflag))
     call check(nf90_inq_varid(id_fin, 'psurf',               id_psurf))
     call check(nf90_inq_varid(id_fin, 'co2_profile_apriori', id_priorpro))
     call check(nf90_inq_varid(id_fin, 'xco2_avgker',         id_avgker))
     call check(nf90_inq_varid(id_fin, 'xco2_apriori',        id_priorobs))
     call check(nf90_inq_varid(id_fin, 'xco2_uncert',         id_uncert))
     call check(nf90_inq_varid(id_fin, 'xco2_final',          id_obs))

     call check(nf90_get_var(id_fin, id_date,     idates))
     call check(nf90_get_var(id_fin, id_time,     itimes))
     call check(nf90_get_var(id_fin, id_lat,      satlats))
     call check(nf90_get_var(id_fin, id_lon,      satlons))
     call check(nf90_get_var(id_fin, id_qcflag,   qcflags))
     call check(nf90_get_var(id_fin, id_psurf,    psurfs))
     call check(nf90_get_var(id_fin, id_priorpro, priorpros))
     call check(nf90_get_var(id_fin, id_avgker,   avgkers1))
     call check(nf90_get_var(id_fin, id_priorobs, priorobses1))
     call check(nf90_get_var(id_fin, id_uncert,   uncerts1))
     call check(nf90_get_var(id_fin, id_obs,      obses1))

  else if (obstype == 'mopitt') then
     call check(nf90_inq_varid(id_fin, 'Year',                   id_year))
     call check(nf90_inq_varid(id_fin, 'Month',                  id_mon))
     call check(nf90_inq_varid(id_fin, 'Day',                    id_day))
     call check(nf90_inq_varid(id_fin, 'Hour',                   id_hour))
     call check(nf90_inq_varid(id_fin, 'Minute',                 id_minu))
     call check(nf90_inq_varid(id_fin, 'Latitude',               id_lat))
     call check(nf90_inq_varid(id_fin, 'Longitude',              id_lon))
     call check(nf90_inq_varid(id_fin, 'SurfacePressure',        id_psurf))
     call check(nf90_inq_varid(id_fin, 'APrioriCOFullProfile',   id_priorpro))
     call check(nf90_inq_varid(id_fin, 'RetrievalAveragingKernelMatrix',       &
                               id_avgker))
     call check(nf90_inq_varid(id_fin, 'RetrievalUncertainty',   id_uncert))
     call check(nf90_inq_varid(id_fin, 'RetrievedCOFullProfile', id_obs))

     call check(nf90_get_var(id_fin, id_year,     iyears))
     call check(nf90_get_var(id_fin, id_mon,      imons))
     call check(nf90_get_var(id_fin, id_day,      idays))
     call check(nf90_get_var(id_fin, id_hour,     ihours))
     call check(nf90_get_var(id_fin, id_minu,     iminus))
     call check(nf90_get_var(id_fin, id_lat,      satlats))
     call check(nf90_get_var(id_fin, id_lon,      satlons))
     call check(nf90_get_var(id_fin, id_psurf,    psurfs))
     call check(nf90_get_var(id_fin, id_priorpro, priorpros))
     call check(nf90_get_var(id_fin, id_avgker,   avgkers))
     call check(nf90_get_var(id_fin, id_uncert,   uncerts))
     call check(nf90_get_var(id_fin, id_obs,      obses))

  else if (obstype == 'flask') then
     call check(nf90_inq_varid(id_fin, 'date', id_date))
     call check(nf90_inq_varid(id_fin, 'time', id_time))
     call check(nf90_inq_varid(id_fin, 'lat',  id_lat))
     call check(nf90_inq_varid(id_fin, 'lon',  id_lon))
     call check(nf90_inq_varid(id_fin, 'std',  id_uncert))
     call check(nf90_inq_varid(id_fin, 'obs',  id_obs))

     call check(nf90_get_var(id_fin, id_date,   idates))
     call check(nf90_get_var(id_fin, id_time,   itimes))
     call check(nf90_get_var(id_fin, id_lat,    satlats))
     call check(nf90_get_var(id_fin, id_lon,    satlons))
     call check(nf90_get_var(id_fin, id_uncert, uncerts))
     call check(nf90_get_var(id_fin, id_obs,    obses))

  end if

  call check(nf90_close(id_fin))

! Loop over the records in the data
110 continue

  if (nele < nread + 1) goto 150
  nread = nread + 1

  if (obstype == 'acos' .or. obstype == 'flask') then
     idate5(1) = idates(nread)/10000                              ! year
     idate5(2) = idates(nread)/100 - idate5(1)*100                ! month
     idate5(3) = idates(nread) - idate5(2)*100 - idate5(1)*10000  ! day
     idate5(4) = itimes(nread)/10000                              ! hour
     idate5(5) = itimes(nread)/100 - idate5(4)*100                ! minute

  else if (obstype == 'mopitt') then
     idate5(1) = iyears(nread)
     idate5(2) = imons(nread)
     idate5(3) = idays(nread)
     idate5(4) = ihours(nread)
     idate5(5) = iminus(nread)
  end if

! Convert observation time to relative time
  call w3fs21(idate5, nmind)
  grdtime = real(nmind - iwinbgn, r_kind)*r60inv

! Convert observation lat/lon to relative position
  deglat = satlats(nread)
  deglon = satlons(nread)
  if (deglon < zero) deglon = deglon + r360

! Toss out any obs with suspicious locations
  if (r90 < abs(deglat) .or. deglon < zero .or. r360 < deglon) goto 110

  radlat = deglat * deg2rad
  radlon = deglon * deg2rad
 
  if (regional) then
     call tll2xy(radlon, radlat, grdlon, grdlat, loutside)
     if (loutside) goto 110
  else
     grdlat = radlat
     grdlon = radlon
     call grdcrd1(grdlat, rlats, nlat, 1)
     call grdcrd1(grdlon, rlons, nlon, 1)
  end if

  if (ldebug) then
     print *, '---', nread, '---'
     print *, 'date     = ', idate5(1), idate5(2), idate5(3)
     print *, 'time     = ', idate5(4), idate5(5)
     print *, 'grdtime  = ', grdtime
     print *, 'lat      = ', satlats(nread), 'lon      = ', satlons(nread)
     print *, 'grdlat   = ', grdlat,         'grdlon   = ', grdlon
  end if

! Write record to output file
  ndata = ndata + 1

  if (obstype == 'acos') then
     tgasout(1,ndata) = 0.0
     tgasout(2,ndata) = grdtime
     tgasout(3,ndata) = grdlon
     tgasout(4,ndata) = grdlat
     tgasout(5,ndata) = deglon
     tgasout(6,ndata) = deglat
     tgasout(7,ndata) = qcflags(nread)
     tgasout(8,ndata) = psurfs(nread)
     do k = 1,npro
        tgasout(k+8,ndata) = priorpros(k,nread)
     end do
     do k = 1,npro
        tgasout(k+8+npro,ndata) = avgkers1(k,nread)
     end do
     tgasout(1+8+2*npro,ndata) = priorobses1(nread)
     tgasout(2+8+2*npro,ndata) = uncerts1(nread)
     tgasout(3+8+2*npro,ndata) = obses1(nread)

     if (ldebug) then
        print *, 'psurf    = ', psurfs(nread)
        print *, 'priorpro = '
        print *, priorpros(:,nread)
        print *, 'avgker   = '
        print *, avgkers1(:,nread)
        print *, 'priorobs = ', priorobses1(nread)
        print *, 'uncert   = ', uncerts1(nread)
        print *, 'obs      = ', obses1(nread)
     end if

  else if (obstype == 'mopitt') then
     tgasout(1,ndata) = 0.0
     tgasout(2,ndata) = grdtime
     tgasout(3,ndata) = grdlon
     tgasout(4,ndata) = grdlat
     tgasout(5,ndata) = deglon
     tgasout(6,ndata) = deglat
     tgasout(7,ndata) = 0.0
     tgasout(8,ndata) = psurfs(nread)
     do k = 1,npro
        tgasout(k+8,ndata) = priorpros(k,nread)
     end do 
     do k = 1,npro
        do j = 1,npro
           tgasout(k+(j-1)*npro+8+npro,ndata) = avgkers(j,k,nread)
        end do 
     end do 
     do k = 1,npro
        tgasout(k+8+npro+npro**2,ndata) = uncerts(k,nread)
     end do 
     do k = 1,npro
        tgasout(k+8+2*npro+npro**2,ndata) = obses(k,nread)
     end do

     if (ldebug) then
        print *, 'psurf    = ', psurfs(nread)
        print *, 'priorpro = '
        print *, priorpros(:,nread)
        print *, 'avgker   = '
        print *, avgkers(:,:,nread)
        print *, 'uncert   = '
        print *, uncerts(:,nread)
        print *, 'obs      = '
        print *, obses(:,nread)
     end if

  else if (obstype == 'flask') then
     tgasout(1,ndata)  = 0.0
     tgasout(2,ndata)  = grdtime
     tgasout(3,ndata)  = grdlon
     tgasout(4,ndata)  = grdlat
     tgasout(5,ndata)  = deglon
     tgasout(6,ndata)  = deglat
     tgasout(7,ndata)  = 0.0
     tgasout(8,ndata)  = 0.0

     do k = 1,nchanl
        tgasout(k+8,ndata) = uncerts(k,nread)
     end do
     do k = 1,nchanl
        tgasout(k+8+nchanl,ndata) = obses(k,nread)
     end do

     if (ldebug) then
        print *, 'uncert   = ', uncerts(:,nread)
        print *, 'obs      = ', obses(:,nread)
     end if

  end if

! Loop back to read next profile
  goto 110

! Jump here when eof detected
150 continue

  nodata = ndata*nchanl

  if (ldebug) then
     print *, '---'
     print *, 'nread    = ', nread
     print *, 'ndata    = ', ndata
     print *, 'nodata   = ', nodata
  end if

! Write header record and data to output file for further processing
  write(lunout) obstype, sis, nreal, nchanl, ilat, ilon
  write(lunout) ((tgasout(j,k),j=1,ntgasdat),k=1,ndata)

160 continue

! Deallocate local arrays
  deallocate(tgasout)

  if (obstype == 'acos') then
     deallocate(idates, itimes, satlats, satlons, qcflags, psurfs, priorpros,  &
                avgkers1, priorobses1, uncerts1, obses1)

  else if (obstype == 'mopitt') then
     deallocate(iyears, imons, idays, ihours, iminus, satlats, satlons,        &
                psurfs, priorpros, avgkers, uncerts, obses)

  else if (obstype == 'flask') then
     deallocate(idates, itimes, satlats, satlons, uncerts, obses)

  end if

  return
  
contains
  subroutine check(status)
    use mpeu_util, only: perr, die
    use netcdf,    only: nf90_noerr, nf90_strerror

    integer, intent(in) :: status

    character(len=*), parameter :: myname_ = 'READ_TGAS.CHECK'

    if (status /= nf90_noerr) then 
       call perr(myname_, 'nf90_strerror = ', trim(nf90_strerror(status)))
       call die(myname_)
    end if
  end subroutine check 
end subroutine read_tgas
