module m_nc_ncepsfc
use netcdf
implicit none
private

public :: nc_ncepsfc_vars_init
public :: nc_ncepsfc_vars_final
public :: nc_ncepsfc_vars
public :: nc_ncepsfc_dims
public :: nc_ncepsfc_read
public :: nc_ncepsfc_write
public :: nc_ncepsfc_getpointer

type nc_ncepsfc_vars
   logical :: initialized=.false.
   integer :: nlon,nlat
   integer :: nveg
   real(4),pointer,dimension(:,:)  :: vtype,stype,vfrac
   real(4),pointer,dimension(:,:)  :: v2d
end type nc_ncepsfc_vars

character(len=*), parameter :: myname = 'm_nc_ncepsfc'

integer, parameter :: nv2dx = 3
character(len=5),parameter :: cvars2dx(nv2dx) = (/ 'stype', 'vtype', 'vfrac' /)

interface nc_ncepsfc_dims; module procedure    &
  read_dims_ ; end interface
interface nc_ncepsfc_read; module procedure    &
  read_ncepsfc_ ; end interface
interface nc_ncepsfc_write; module procedure    &
  write_ncepsfc_ ; end interface
interface nc_ncepsfc_vars_init; module procedure    &
  init_ncepsfc_vars_ ; end interface
interface nc_ncepsfc_vars_final; module procedure    &
  final_ncepsfc_vars_ ; end interface
interface nc_ncepsfc_getpointer 
  module procedure get_pointer_2d_ 
end interface

contains

subroutine read_dims_ (fname,nlat,nlon,rc, myid,root)
  implicit none
  character(len=*), intent(in)    :: fname ! input filename
  integer, intent(out) :: rc
  integer, intent(out) :: nlat,nlon
  integer, intent(in), optional :: myid, root

! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid, ier
  integer :: mype_,root_

! Local variables
  character(len=*), parameter :: myname_ = myname//"::dims_"
   
! Return code (status)
  rc=0; mype_=0; root_=0
  if(present(myid) .and. present(root) ) then
     mype_ = myid
     root_ = root
  endif
 
! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
! the file.

  call check_( nf90_open(fname, NF90_NOWRITE, ncid), rc, mype_, root_ )
  if(rc/=0) return

! Read global attributes
  call check_( nf90_inq_dimid(ncid, "lon", varid), rc, mype_, root_)
  call check_( nf90_inquire_dimension(ncid, varid, len=nlon), rc, mype_, root_ )
  call check_( nf90_inq_dimid(ncid, "lat", varid), rc, mype_, root_ )
  call check_( nf90_inquire_dimension(ncid, varid, len=nlat), rc, mype_, root_ )

! Close the file, freeing all resources.
  call check_( nf90_close(ncid), rc, mype_, root_ )

  return

end subroutine read_dims_

subroutine read_ncepsfc_ (fname,sfcvar,rc, myid,root)
  implicit none
  character(len=*), intent(in)    :: fname ! input filename
  type(nc_ncepsfc_vars),intent(inout) :: sfcvar ! background error variables
  integer, intent(out) :: rc
  integer, intent(in), optional :: myid,root ! accommodate MPI calling programs

! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid

! Local variables
  character(len=*), parameter :: myname_ = myname//"::read_"
  character(len=4) :: cindx
  integer :: nv,nlat,nlon
  integer :: ndims_, nvars_, ngatts_, unlimdimid_
  integer :: nlat_,nlon_
  integer :: mype_,root_
  real(4), allocatable :: data_in(:,:,:)
  logical :: verbose
  logical :: init_
  
   
! Return code (status)
  rc=0; mype_=0; root_=0
  verbose=.true.
  init_=.false.
  if(present(myid).and.present(root) )then
    if(myid/=root) verbose=.false.
    mype_ = myid
    root_ = root
  endif
 
! Get dimensions
  call read_dims_ (fname,nlat_,nlon_,rc, mype_,root_)

  init_ = sfcvar%initialized
  if ( init_ ) then
!   Set dims
    nlat=sfcvar%nlat
    nlon=sfcvar%nlon

!   Consistency check
    if (nlon_ /= nlon .or. nlat_ /=nlat) then
       rc=1
       if(myid==root) then
         print *, 'nlat(file) = ', nlat_, 'nlat(required) = ', nlat
         print *, 'nlon(file) = ', nlon_, 'nlon(required) = ', nlon
         print *, myname_,  'Inconsistent dimensions, aborting ... '
       endif
       return
    endif
  else
!   Set dims
    nlat=nlat_
    nlon=nlon_
    call init_ncepsfc_vars_(sfcvar,nlon,nlat)
  endif

! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
! the file.

  call check_( nf90_open(fname, NF90_NOWRITE, ncid), rc, mype_, root_ )
  if(rc/=0) return

! Read global attributes
! call check_( nf90_inquire(ncid, ndims_, nvars_, ngatts_, unlimdimid_), rc, mype_, root_ )
! call check_( nf90_inq_dimid(ncid, "lon", varid), rc, mype_, root_ )
! call check_( nf90_inquire_dimension(ncid, varid, len=nlon_), rc, mype_, root_ )
! call check_( nf90_inq_dimid(ncid, "lat", varid), rc, mype_, root_ )
! call check_( nf90_inquire_dimension(ncid, varid, len=nlat_), rc, mype_, root_ )
! call check_( nf90_inq_dimid(ncid, "lev", varid), rc, mype_, root_ )
! call check_( nf90_inquire_dimension(ncid, varid, len=nlev_), rc, mype_, root_ )

! Read in lat/lon fields
  allocate(data_in(nlon,nlat,1))
  do nv = 1, nv2dx
     call check_( nf90_inq_varid(ncid, trim(cvars2dx(nv)), varid), rc, mype_, root_ )
     call check_( nf90_get_var(ncid, varid, data_in(:,:,1)), rc, mype_, root_ )
     if(trim(cvars2dx(nv))=="vfrac"     ) then
        sfcvar%vfrac = transpose(data_in(:,:,1))
     endif
     if(trim(cvars2dx(nv))=="vtype" ) then 
        sfcvar%vtype = transpose(data_in(:,:,1))
     endif
     if(trim(cvars2dx(nv))=="stype" ) then 
        sfcvar%stype = transpose(data_in(:,:,1))
     endif
  enddo
  deallocate(data_in)

  sfcvar%nveg = nint(maxval(sfcvar%vtype))+1

! Close the file, freeing all resources.
  call check_( nf90_close(ncid), rc, mype_, root_ )

  if(verbose) print *,"*** Finish reading file: ", trim(fname)

  return

end subroutine read_ncepsfc_

subroutine write_ncepsfc_ (fname,sfcvar,plevs,lats,lons,rc, myid,root)
  implicit none
  character(len=*), intent(in)    :: fname ! input filename
  type(nc_ncepsfc_vars),intent(in)    :: sfcvar ! background error variables
  real(4), intent(in) :: lats(:)           ! latitudes  per GSI: increase index from South to North Pole
  real(4), intent(in) :: lons(:)           ! longitudea per GSI: increase index from East to West
  real(4), intent(in) :: plevs(:)
  integer, intent(out) :: rc
  integer, intent(in), optional :: myid,root        ! accommodate MPI calling programs

  character(len=*), parameter :: myname_ = myname//"::write_"
  integer, parameter :: NDIMS = 3

! When we create netCDF files, variables and dimensions, we get back
! an ID for each one.
  character(len=4) :: cindx
  integer :: ncid, dimids(NDIMS)
  integer :: x_dimid, y_dimid, z_dimid
  integer :: lon_varid, lat_varid, lev_varid
  integer :: ii,jj,nl,nv,nn,nlat,nlon
  integer :: mype_,root_
  integer, allocatable :: varid1d(:), varid2d(:), varid2dx(:), varidMLL(:)
  logical :: verbose
  
! This is the data array we will write. It will just be filled with
! a progression of integers for this example.
  real(4), allocatable :: data_out(:,:,:)

! Return code (status)
  rc=0; mype_=0; root_=0
  verbose=.true.
  if(present(myid).and.present(root) )then
    if(myid/=root) verbose=.false.
    mype_ = myid
    root_ = root
  endif

! Set dims
  nlat=sfcvar%nlat
  nlon=sfcvar%nlon

! Always check the return code of every netCDF function call. In
! this example program, wrapping netCDF calls with "call check()"
! makes sure that any return which is not equal to nf90_noerr (0)
! will print a netCDF error message and exit.

! Create the netCDF file. The nf90_clobber parameter tells netCDF to
! overwrite this file, if it already exists.
  call check_( nf90_create(fname, NF90_CLOBBER, ncid), rc, mype_, root_ )
  if(rc/=0) return

! Define the dimensions. NetCDF will hand back an ID for each. 
  call check_( nf90_def_dim(ncid, "lon", nlon, x_dimid), rc, mype_, root_ )
  call check_( nf90_def_dim(ncid, "lat", nlat, y_dimid), rc, mype_, root_ )

  call check_( nf90_def_var(ncid, "lon", NF90_REAL, x_dimid, lon_varid), rc, mype_, root_ )
  call check_( nf90_def_var(ncid, "lat", NF90_REAL, y_dimid, lat_varid), rc, mype_, root_ )
  call check_( nf90_def_var(ncid, "lev", NF90_REAL, z_dimid, lev_varid), rc, mype_, root_ )

  call check_( nf90_put_att(ncid, lon_varid, "units", "degress"), rc, mype_, root_ )
  call check_( nf90_put_att(ncid, lat_varid, "units", "degress"), rc, mype_, root_ )
  call check_( nf90_put_att(ncid, lev_varid, "units", "hPa"), rc, mype_, root_ )

! The dimids array is used to pass the IDs of the dimensions of
! the variables. Note that in fortran arrays are stored in
! column-major format.
  dimids =  (/ x_dimid, y_dimid, z_dimid /)

! Define variables.
  allocate(varid2dx(nv2dx))
  do nv = 1, nv2dx
     call check_( nf90_def_var(ncid, trim(cvars2dx(nv)), NF90_REAL, (/ x_dimid, y_dimid /), varid2dx(nv)), rc, mype_, root_ )
  enddo

! End define mode. This tells netCDF we are done defining metadata.
  call check_( nf90_enddef(ncid), rc, mype_, root_ )

! Write coordinate variables data
  call check_( nf90_put_var(ncid, lon_varid, lons ), rc, mype_, root_ )
  call check_( nf90_put_var(ncid, lat_varid, lats ), rc, mype_, root_ )

! Write out lat/lon fields
  allocate(data_out(nlon,nlat,1))
  do nv = 1, nv2dx
     if(trim(cvars2dx(nv))=="vfrac" ) then
        data_out(:,:,1) = transpose(sfcvar%vfrac)
     endif
     if(trim(cvars2dx(nv))=="vtype" ) then 
        data_out(:,:,1) = transpose(sfcvar%vtype)
     endif
     if(trim(cvars2dx(nv))=="stype" ) then 
        data_out(:,:,1) = transpose(sfcvar%stype)
     endif
     call check_( nf90_put_var(ncid, varid2dx(nv), data_out(:,:,1)), rc, mype_, root_ )
  enddo
  deallocate(data_out)

! Close file
  call check_( nf90_close(ncid), rc, mype_, root_ )

  deallocate(varidMLL)
  deallocate(varid2d)
  deallocate(varid1d)

  if(verbose) print *,"*** Finish writing file: ", trim(fname)

  return

end subroutine write_ncepsfc_

subroutine init_ncepsfc_vars_(vr,nlon,nlat)

  integer,intent(in) :: nlon,nlat
  type(nc_ncepsfc_vars) vr

  if(vr%initialized) return

  vr%nlon=nlon 
  vr%nlat=nlat

! allocate arrays
  allocate(vr%vfrac(nlat,nlon),vr%vtype(nlat,nlon),vr%stype(nlat,nlon))
  vr%initialized=.true.
  end subroutine init_ncepsfc_vars_

  subroutine final_ncepsfc_vars_(vr)
  type(nc_ncepsfc_vars) vr
! deallocate arrays
  if(.not. vr%initialized) return
  deallocate(vr%vtype,vr%vfrac,vr%stype)
  vr%initialized=.false.
end subroutine final_ncepsfc_vars_

subroutine get_pointer_2d_ (vname, sfcvar, ptr, rc )
implicit none
character(len=*), intent(in) :: vname
type(nc_ncepsfc_vars) sfcvar
real(4),pointer,intent(inout) :: ptr(:,:)
integer,intent(out) :: rc
rc=-1
if(trim(vname)=='vtype') then
  ptr => sfcvar%vtype
  rc=0
endif
if(trim(vname)=='vfrac') then
  ptr => sfcvar%vfrac
  rc=0
endif
if(trim(vname)=='stype') then
  ptr => sfcvar%stype
  rc=0
endif
end subroutine get_pointer_2d_

subroutine check_(status,rc, myid, root)
    integer, intent ( in) :: status
    integer, intent (out) :: rc
    integer, intent ( in) :: myid, root
    rc=0
    if(status /= nf90_noerr) then 
      if(myid==root) print *, trim(nf90_strerror(status))
      rc=999
    end if
end subroutine check_  

end module m_nc_ncepsfc
