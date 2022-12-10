subroutine read_mls_molec(nread,ndata,nodata,molec_type, dfile,dtype,lunout,gstime,twind,&
                          sis)

  use netcdf, only: nf90_open
  use netcdf, only: nf90_nowrite
  use netcdf, only: nf90_noerr
  use netcdf, only: nf90_inq_dimid
  use netcdf, only: nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid
  use netcdf, only: nf90_get_var
  use netcdf, only: nf90_close
  use gridmod, only: nlat,nlon,regional,tll2xy,rlats,rlons
  use gsi_4dvar, only: l4dvar,iwinbgn,winlen,l4densvar
  use constants, only: deg2rad,zero,one,rad2deg,one_tenth,r60inv
  use tgasinfo, only: jpch_tgas,nusis_tgas,iuse_tgas
  use kinds, only: i_kind,r_kind,r_double
  use m_extOzone, only: check
  use mpeu_util, only: perr,die
 
  implicit none

  character(len=*), intent(in):: molec_type ! "h2o", "hcl", "hno3",...
  character(len=*), intent(in):: dfile   ! obs_input filename
  character(len=*), intent(in):: dtype   ! obs_input dtype
  character(len=*), intent(in):: sis     ! satellite/instrument/sensor indicator
  integer(i_kind),  intent(in   ) :: lunout

  real   (kind=r_kind), intent(in):: gstime ! analysis time (minute) from reference date
  real   (kind=r_kind), intent(in):: twind  ! input group time window (hour)
  integer(kind=i_kind), intent(inout):: nread, ndata, nodata
  character(len=*) , parameter:: myname='read_mls_constituent'

  integer(kind=i_kind):: ier, iprof, maxobs
  integer(kind=i_kind):: nsound, nreal, nchanl, nmolecdat
  integer(kind=i_kind):: i, ilev, iout, ikx, ncid, k0, nmrecs
  integer(kind=i_kind),allocatable,dimension(:):: ipos

  real   (kind=r_kind), allocatable :: molecout(:,:)
  integer(kind=i_kind):: nrecDimId,lonVarID,latVarID,yyVarID,mmVarID
  integer(kind=i_kind):: ddVarID,hhVarID,minVarID,ssVarID
  integer(kind=i_kind):: pressVarID
  integer(kind=i_kind):: qualVarID, errVarID, molecVarID
  integer(kind=i_kind):: levsDimID,levs

  integer(kind=i_kind), allocatable :: iya(:),ima(:),idda(:),ihha(:),imina(:),iseca(:)
  real   (kind=r_kind), allocatable :: slatsa(:),slonsa(:)
  real   (kind=r_kind), allocatable :: press(:), molecule(:,:), usage(:,:)
  real   (kind=r_kind), allocatable :: err(:,:)

  integer(kind=i_kind):: nmind,ilat,ilon,j,k

  real   (kind=r_kind):: slons0,slats0
  real   (kind=r_kind):: pob

  real   (kind=r_kind):: dlon,dlon_earth,dlon_earth_deg
  real   (kind=r_kind):: dlat,dlat_earth,dlat_earth_deg
  real   (kind=r_kind):: tdiff,sstime,t4dv,rsat
  integer(kind=i_kind):: idate5(5)
  real   (kind=r_kind), parameter:: rmiss = -9999.9_r_kind
!  real   (kind=r_kind), parameter:: badoz = 10000.0_r_kind
!  real   (kind=r_kind), parameter:: r6    =     6.0_r_kind
  real   (kind=r_kind), parameter:: r360  =   360.0_r_kind

  logical:: outside
  logical:: first

  rsat=999._r_kind
  maxobs = 100000

  ! Open file and read dimensions
  call check(nf90_open(trim(dfile),nf90_nowrite,ncid),stat=ier)

  ! ignore if the file is not actually present.
  if(ier/=nf90_noerr) return

  ! Get dimensions from the input file
  call check(nf90_inq_dimid(ncid, "nprofiles", nrecDimId),stat=ier)

  ! ignore if the file header is empty
  if(ier/=nf90_noerr) then
     call check(nf90_close(ncid),stat=ier)
     return
  endif

  ! Get dimensions from the input file: # of profiles and # of levels
  nsound=0
  call check(nf90_inquire_dimension(ncid, nrecDimId, len = nsound),stat=ier)

  if(nsound==0) then
     call check(nf90_close(ncid),stat=ier)
     return
  endif

  ! Continue the input
  call check(nf90_inq_dimid(ncid, "nlevs", levsDimId))
  call check(nf90_inquire_dimension(ncid, levsDimId, len = levs))

  nchanl = levs
  nreal  = 6 + 3*nchanl
  nmolecdat = nreal + nchanl 

  allocate(molecout(nmolecdat, maxobs))

  ilon=3    
  ilat=4

  ndata  = 0

  !  NOTE: Make sure that the 'tgasinfo' file has the same number of levels
  allocate(ipos(levs))
  ipos=999

  ! Process limb data in NetDCF format
  ikx = 0 
  first=.false.
  do i=1,jpch_tgas
     if( (.not. first) .and. index(nusis_tgas(i), trim(dtype))/=0) then
        k0=i
        first=.true.
     end if
     if(first.and.index(nusis_tgas(i),trim(dtype))/=0) then 
        ikx=ikx+1
        ipos(ikx)=k0+ikx-1
     end if
  end do

  if (ikx/=levs) call die(myname//': inconsistent levs for '//dtype)

  nmrecs=0
  ! Allocate space and read data
  allocate(iya(nsound),ima(nsound),idda(nsound),ihha(nsound),imina(nsound), &
       iseca(nsound),slatsa(nsound),slonsa(nsound), molecule(levs,nsound),  &
       usage(levs,nsound), err(levs,nsound), press(levs))

  ! Read variables and store them in these arrays
  call check(nf90_inq_varid(ncid, "lon", lonVarId))
  call check(nf90_get_var(ncid, lonVarId, slonsa))

  call check(nf90_inq_varid(ncid, "lat", latVarId))
  call check(nf90_get_var(ncid, latVarId, slatsa))

  call check(nf90_inq_varid(ncid, "yy", yyVarId))
  call check(nf90_get_var(ncid, yyVarId, iya))

  call check(nf90_inq_varid(ncid, "mm", mmVarId))
  call check(nf90_get_var(ncid, mmVarId, ima))

  call check(nf90_inq_varid(ncid, "dd", ddVarId))
  call check(nf90_get_var(ncid, ddVarId, idda))

  call check(nf90_inq_varid(ncid, "hh", hhVarId))
  call check(nf90_get_var(ncid, hhVarId, ihha))

  call check(nf90_inq_varid(ncid, "min", minVarId))
  call check(nf90_get_var(ncid, minVarId, imina))

  call check(nf90_inq_varid(ncid, "ss", ssVarId))
  call check(nf90_get_var(ncid, ssVarId, iseca))

  call check(nf90_inq_varid(ncid, "press", pressVarId))
  call check(nf90_get_var(ncid, pressVarId, press))

  call check(nf90_inq_varid(ncid, "oberr", errVarId))
  call check(nf90_get_var(ncid, errVarId, err))

  call check(nf90_inq_varid(ncid, trim(molec_type), molecVarId))
  call check(nf90_get_var(ncid, molecVarId, molecule))

  ! close the data file
  call check(nf90_close(ncid))

  ! 'Unpack' the data
  nmrecs = 0
  do iprof = 1,nsound
     do ilev = 1,levs
        if (press(ilev)          .lt. -900.0 .or. &
            err(ilev,iprof)      .lt. -900.0 .or. &
            molecule(ilev,iprof) .lt. -900.0 .or. &
            iuse_tgas(ipos(ilev)) < 0) then
           usage(ilev,iprof) = one
        else
           usage(ilev,iprof) = zero
        endif
        nmrecs=nmrecs+1

        !       censor very negative data
        if (molecule(ilev,iprof) < -err(ilev,iprof)) usage(ilev,iprof) = one

        !       convert observation location to radians
        slons0=slonsa(iprof)
        slats0=slatsa(iprof)
        if(abs(slats0)>90._r_kind .or. abs(slons0)>r360) cycle
        if(slons0< zero) slons0=slons0+r360
        if(slons0==r360) slons0=zero
        dlat_earth_deg = slats0
        dlon_earth_deg = slons0
        dlat_earth = slats0 * deg2rad
        dlon_earth = slons0 * deg2rad

        if(regional)then
           call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
           if(outside) cycle    
        else
           dlat = dlat_earth
           dlon = dlon_earth
           call grdcrd1(dlat,rlats,nlat,1)
           call grdcrd1(dlon,rlons,nlon,1)
        endif

        idate5(1) = iya(iprof) !year
        idate5(2) = ima(iprof) !month
        idate5(3) = idda(iprof) !day
        idate5(4) = ihha(iprof) !hour
        idate5(5) = imina(iprof) !minute
        call w3fs21(idate5,nmind)
        t4dv=real((nmind-iwinbgn),r_kind)*r60inv
        if (l4dvar.or.l4densvar) then
           if (t4dv<zero .OR. t4dv>winlen) then
              write(6,*)'read_mls_molec: ', dtype,' obs time idate5=',idate5,', t4dv=',&
                   t4dv,' is outside time window, sstime=',sstime*r60inv
              cycle
           end if
        else
           sstime=real(nmind,r_kind)
           tdiff=(sstime-gstime)*r60inv
           if(abs(tdiff) > twind)then
              write(6,*)'read_mls_molec: ',dtype,' obs time idate5=',idate5,', tdiff=',&
                   tdiff,' is outside time window=',twind
              cycle
           end if
        end if
     end do

     ndata  = ndata+1
     if(ndata<=maxobs) then
        molecout(1,ndata)=rsat
        molecout(2,ndata)=t4dv
        molecout(3,ndata)=dlon                 ! grid relative longitude
        molecout(4,ndata)=dlat                 ! grid relative latitude
        molecout(5,ndata)=dlon_earth_deg       ! earth relative longitude (degrees)
        molecout(6,ndata)=dlat_earth_deg       ! earth relative latitude (degrees)

        iout = 6
        do ilev = 1,levs
           iout = iout + 1
           molecout(iout,ndata) = press(ilev)               ! pressure 
        end do
        do ilev = 1,levs
           iout = iout + 1
           molecout(iout,ndata) = usage(ilev,iprof)         ! quality flag
        end do
        do ilev = 1,levs
           iout = iout + 1
           molecout(iout,ndata) = err(ilev,iprof)           ! mixing ratio precision in ppmv
        end do

        do ilev = 1,levs
           iout = iout + 1
           molecout(iout,ndata) = molecule(ilev,iprof)      ! mixing ratio in ppmv
        end do
     endif
  end do

  nread  = nsound * levs
  nodata = ndata  * levs

! Write header record and data to output file for further processing
  write(lunout) dtype, sis, nreal, nchanl, ilat, ilon
  write(lunout) ((molecout(j,k),j=1,nmolecdat),k=1,ndata)

  deallocate(iya,ima,idda,ihha,imina,iseca,slatsa,slonsa, molecule, &
             err, usage, press)
  deallocate(ipos,molecout)


  if (ndata > maxobs) then 
     call perr('read_mls_molec','Number of MLS obs reached maxobs = ', maxobs)
     call perr(myname,' Number of MLS molec. obs reached maxobs = ', maxobs)
     call perr(myname,'                           ndata = ', ndata)
     call perr(myname,'                          nodata = ', nodata)
     call die(myname)
  endif

end subroutine read_mls_molec

