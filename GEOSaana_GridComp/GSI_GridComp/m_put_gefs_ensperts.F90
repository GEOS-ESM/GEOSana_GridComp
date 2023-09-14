module m_put_gefs_ensperts
private
public :: ens_spread_dualres
public :: write_spread_dualres
interface ens_spread_dualres
  module procedure ens_spread_dualres_
end interface
contains
subroutine ens_spread_dualres_(ibin,prefix,en_bar)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    ens_spread_dualres  output ensemble spread for diagnostics
!   prgmmr: kleist           org: np22                date: 2010-01-05
!
! abstract: compute ensemble spread on ensemble grid, interpolate to analysis grid
!             and write out for diagnostic purposes.
!
!
! program history log:
!   2010-01-05  kleist, initial documentation
!   2010-02-28  parrish - make changes to allow dual resolution capability
!   2011-03-19  parrish - add pseudo-bundle capability
!   2011-11-01  kleist  - 4d capability for ensemble/hybrid
!   2019-07-10  todling - truly handling 4d output; and upd to out all ens c-variables
!   2020-05-12  todling - fix normalization from 1/M to 1/(M-1)
!
!   input argument list:
!     ibin  - time bin index
!   optional input argument list:
!     en_bar - ensemble mean
!     en_bar - prefix for naming output files
!
!   output argument list:
!
! NOTE: The update made by Dave to handle dual resolution leads to non-positive
!       spreads - the interpolation routines are not quarantee to preserve 
!       positiveness (Todling). A better version of this routine would
!       interpolate the errors (x(m)-xbar) and then proceed to calculate the
!       spreads; but it would be more costly.
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block
  use mpimod, only: mype
  use kinds, only: r_single,r_kind,i_kind
  use hybrid_ensemble_parameters, only: n_ens,grd_ens,grd_anl,p_e2a,uv_hyb_ens
  use hybrid_ensemble_parameters, only: en_perts,nelen
  use general_sub2grid_mod, only: sub2grid_info,general_sub2grid_create_info,general_sube2suba
  use constants, only:  zero,two,half,one
  use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
  use control_vectors, only: be2d,be3d
  use mpeu_util, only: getindex   
  use gsi_bundlemod, only: gsi_bundlecreate
  use gsi_bundlemod, only: gsi_grid
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_bundlemod, only: gsi_bundledestroy
  use gsi_bundlemod, only: gsi_gridcreate
  implicit none

  integer(i_kind),intent(in):: ibin
  character(len=*),optional,intent(in):: prefix
  type(gsi_bundle),optional,intent(in):: en_bar

  type(gsi_bundle):: sube,suba
  type(gsi_grid):: grid_ens,grid_anl
  real(r_kind) sp_norm
  type(sub2grid_info)::se,sa

  integer(i_kind) i,n,ic3,k,ipic
  logical regional
  integer(i_kind) ipc3d(nc3d),ipc2d(nc2d)
  integer(i_kind) num_fields,inner_vars,istatus
  logical,allocatable::vector(:)
  character(len=80) prefix_

  prefix_="NULL"
  if(present(prefix)) then
    prefix_ = prefix
  endif
!      create simple regular grid
  call gsi_gridcreate(grid_anl,grd_anl%lat2,grd_anl%lon2,grd_anl%nsig)
  call gsi_gridcreate(grid_ens,grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)

! create two internal bundles, one on analysis grid and one on ensemble grid

  call gsi_bundlecreate (suba,grid_anl,'ensemble work',istatus, &
                                 names2d=cvars2d,names3d=cvars3d)
  if(istatus/=0) then
     write(6,*)' ens_spread_dualres: trouble creating bundle_anl bundle'
     call stop2(999)
  endif
  call gsi_bundlecreate (sube,grid_ens,'ensemble work ens',istatus, &
                            names2d=cvars2d,names3d=cvars3d)
  if(istatus/=0) then
     write(6,*)' ens_spread_dualres: trouble creating bundle_ens bundle'
     call stop2(999)
  endif

  sp_norm=(one/max(one,n_ens-one))

  sube%values=zero
  if (present(en_bar)) then
     do n=1,n_ens
        do i=1,nelen
           sube%values(i)=sube%values(i) &
              +(en_perts(n,ibin)%valuesr4(i)-en_bar%values(i))*(en_perts(n,ibin)%valuesr4(i)-en_bar%values(i))
        end do
     end do
  else
     do n=1,n_ens
        do i=1,nelen
           sube%values(i)=sube%values(i) &
              +(en_perts(n,ibin)%valuesr4(i))*(en_perts(n,ibin)%valuesr4(i))
        end do
     end do
  endif

  do i=1,nelen
    sube%values(i) = sqrt(sp_norm*sube%values(i))
  end do

! apply anavinfo factors
  call gsi_bundlegetpointer (sube,cvars3d,ipc3d,istatus)
  if(istatus/=0) then
    write(6,*) ' ens_spread_dualres',': cannot find 3d pointers'
    call stop2(999)
  endif
  call gsi_bundlegetpointer (sube,cvars2d,ipc2d,istatus)
  if(istatus/=0) then
    write(6,*) ' ens_spread_dualres',': cannot find 2d pointers'
    call stop2(999)
  endif
  do n=1,nc2d
     ipic=ipc2d(n)
     if (be2d(ipic)<zero) cycle
     sube%r2(ipic)%q = be2d(ipic)*sube%r2(ipic)%q
  end do
  do n=1,nc3d
     ipic=ipc3d(n)
     if (be3d(ipic)<zero) cycle
     sube%r3(ipic)%q = be3d(ipic)*sube%r3(ipic)%q
  end do

  if(grd_ens%latlon1n == grd_anl%latlon1n) then
     do i=1,nelen
        suba%values(i)=sube%values(i)
     end do
  else
     regional=.false.
     inner_vars=1
     num_fields=max(0,nc3d)*grd_ens%nsig+max(0,nc2d)
     allocate(vector(num_fields))
     vector=.false.
     do ic3=1,nc3d
        if(trim(cvars3d(ic3))=='sf'.or.trim(cvars3d(ic3))=='vp') then !  RTodling: this is not meaningful
                                                                      !            since are dealing w/ spread
                                                                      !            not the physical u/v
           do k=1,grd_ens%nsig
              vector((ic3-1)*grd_ens%nsig+k)=uv_hyb_ens
           end do
        end if
     end do
     call general_sub2grid_create_info(se,inner_vars,grd_ens%nlat,grd_ens%nlon,grd_ens%nsig,num_fields, &
                                       regional,vector)
     call general_sub2grid_create_info(sa,inner_vars,grd_anl%nlat,grd_anl%nlon,grd_anl%nsig,num_fields, &
                                       regional,vector)
     deallocate(vector)
     call general_sube2suba(se,sa,p_e2a,sube%values,suba%values,regional)
  end if

  call write_spread_dualres(ibin,suba,prefix_)

  return
end subroutine ens_spread_dualres_


subroutine write_spread_dualres(ibin,bundle,prefix)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    write_spread_dualres   write ensemble spread for diagnostics
!   prgmmr: kleist           org: np22                date: 2010-01-05
!
! abstract: write ensemble spread (previously interpolated to analysis grid)
!             for diagnostic purposes.
!
!
! program history log:
!   2010-01-05  kleist, initial documentation
!   2010-02-28  parrish - make changes to allow dual resolution capability
!   2018-04-01  eliu - add hydrometeors 
!   2019-07-10  todling - generalize to write out all variables in the ensemble
!                       - also allows for print out of different time bins
!   2021-10-08  todling - name wind vars correctly in file when ens uses wind vectors
!
!   input argument list:
!     bundle -  spread bundle
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block
  use mpimod, only: mype
  use mpeu_util, only: die
  use mpimod, only: mpi_rtype,mpi_itype,mpi_comm_world
  use kinds, only: r_kind,i_kind,r_single
  use guess_grids, only: get_ref_gesprs 
  use gridmod, only: rlats
  use hybrid_ensemble_parameters, only: grd_anl,uv_hyb_ens
  use gsi_bundlemod, only: gsi_grid
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
  use constants, only: zero
  use jfunc, only: miter,jiter ! should really pass as argument
  use m_revBens, only: spread2d,spread3d
  implicit none

  integer(i_kind), intent(in) :: ibin
  type(gsi_bundle):: bundle
  character(len=*),intent(in) :: prefix

! local variables
  character(255):: grdfile,grdctl,var

  real(r_kind),allocatable,dimension(:,:,:):: work8_3d
  real(r_kind),allocatable,dimension(:,:):: work8_2d

  real(r_single),allocatable,dimension(:,:,:):: work4_3d
  real(r_single),allocatable,dimension(:,:):: work4_2d

  real(r_single),allocatable,dimension(:,:):: cosrlats2d

  real(r_kind),pointer,dimension(:,:,:):: ptr3d
  real(r_kind),pointer,dimension(:,:):: ptr2d

  integer(i_kind) iret,i,j,k,n,mem2d,mem3d,num3d,lu,istat,ifailed
  real(r_kind),dimension(grd_anl%nsig+1) :: prs

  character(len=*),parameter :: myname='write_spread_dualres'

! Initial memory used by 2d and 3d grids
  mem2d = 4*grd_anl%nlat*grd_anl%nlon
  mem3d = 4*grd_anl%nlat*grd_anl%nlon*grd_anl%nsig
  num3d=11

  allocate(work8_3d(grd_anl%nlat,grd_anl%nlon,grd_anl%nsig))
  allocate(work8_2d(grd_anl%nlat,grd_anl%nlon))
  allocate(work4_3d(grd_anl%nlon,grd_anl%nlat,grd_anl%nsig))
  allocate(work4_2d(grd_anl%nlon,grd_anl%nlat))

  if (mype==0) then
    if(trim(prefix)=="NULL") then
       write(grdfile,'(a,2(i3.3,a))') 'ens_spread_',ibin, '.iter' ,jiter, '.grd'
    else
       write(grdfile,'(2a,2(i3.3,a))') trim(prefix), 'ens_spread_',ibin, '.iter' ,jiter, '.grd'
    endif
    call baopenwt(22,trim(grdfile),iret)
    write(6,*)'WRITE_SPREAD_DUALRES:  open 22 to ',trim(grdfile),' with iret=',iret
  endif

  if(mype==0) allocate(cosrlats2d(grd_anl%nlon,grd_anl%nlat))

! Process 3d arrays
  ifailed=0
  do n=1,nc3d
    call gsi_bundlegetpointer(bundle,cvars3d(n),ptr3d,istat)
    work8_3d=zero
    do k=1,grd_anl%nsig
      call gather_stuff2(ptr3d(1,1,k),work8_3d(1,1,k),mype,0)
    end do
    if (mype==0) then
      do i=1,grd_anl%nlat
         cosrlats2d(:,i) = cos(rlats(i))
      enddo
      if(any(cosrlats2d<zero)) ifailed=1
      spread3d(:,n,ibin,jiter)=zero
      do k=1,grd_anl%nsig
        do j=1,grd_anl%nlon
          do i=1,grd_anl%nlat
            work4_3d(j,i,k) =work8_3d(i,j,k)
          end do
        end do
        spread3d(k,n,ibin,jiter)=sqrt(sum(cosrlats2d*work4_3d(:,:,k)*work4_3d(:,:,k))/(grd_anl%nlon*grd_anl%nlat)) ! not quite the proper grid weight
      end do
      call wryte(22,mem3d,work4_3d)
      write(6,*)'WRITE_SPREAD_DUALRES FOR VARIABLE ',trim(cvars3d(n))
    endif
  end do
  call mpi_bcast(ifailed,1,mpi_itype,0,mpi_comm_world,istat)
  call mpi_bcast(spread3d(1,1,ibin,jiter),grd_anl%nsig*nc3d,mpi_rtype,0,mpi_comm_world,istat)
  if(ifailed/=0) then
     call die(myname,': fishy cos(lat), aborting',ifailed)
  endif

! Process 2d array
  do n=1,nc2d
    call gsi_bundlegetpointer(bundle,cvars2d(n),ptr2d,istat)
    work8_2d=zero
    call gather_stuff2(ptr2d,work8_2d,mype,0)
    if (mype==0) then
       spread2d(n,ibin,jiter)=zero
       do j=1,grd_anl%nlon
          do i=1,grd_anl%nlat
             work4_2d(j,i)=work8_2d(i,j)
          end do
       end do
       spread2d(n,ibin,jiter)=sqrt(sum(cosrlats2d*work4_2d*work4_2d)/(grd_anl%nlon*grd_anl%nlat)) ! not quite the proper grid weight
       call wryte(22,mem2d,work4_2d)
       write(6,*)'WRITE_SPREAD_DUALRES FOR VARIABLE ',trim(cvars2d(n))
    endif
  end do
  call mpi_bcast(spread2d(1,ibin,jiter),nc2d,mpi_rtype,0,mpi_comm_world,istat)

  if(mype==0) deallocate(cosrlats2d)

! Close byte-addressable binary file for grads
  if (mype==0) then
     call baclose(22,iret)
     write(6,*)'WRITE_SPREAD_DUALRES:  close 22 with iret=',iret
  end if

! Get reference pressure levels for grads purposes
  call get_ref_gesprs(prs)

! Write out a corresponding grads control file
  if (mype==0) then
     if(trim(prefix)=="NULL") then
       write(grdctl,'(a,2(i3.3,a))') 'ens_spread_',ibin,  '.iter' ,jiter, '.ctl'
     else
       write(grdctl,'(2a,2(i3.3,a))') trim(prefix), 'ens_spread_',ibin,  '.iter' ,jiter, '.ctl'
     endif
     open(newunit=lu,file=trim(grdctl),form='formatted')
     write(lu,'(2a)') 'DSET  ^', trim(grdfile)
     write(lu,'(2a)') 'TITLE ', 'gsi ensemble spread'
     write(lu,'(a,2x,e13.6)') 'UNDEF', 1.E+15 ! any other preference for this?
     write(lu,'(a,2x,i4,2x,a,2x,f5.1,2x,f9.6)') 'XDEF',grd_anl%nlon, 'LINEAR',   0.0, 360./grd_anl%nlon
     write(lu,'(a,2x,i4,2x,a,2x,f5.1,2x,f9.6)') 'YDEF',grd_anl%nlat, 'LINEAR', -90.0, 180./(grd_anl%nlat-1.)
     write(lu,'(a,2x,i4,2x,a,100(1x,f10.5))')      'ZDEF',grd_anl%nsig, 'LEVELS', prs(1:grd_anl%nsig)
     write(lu,'(a,2x,i4,2x,a)')   'TDEF', 1, 'LINEAR 12:00Z04JUL1776 6hr' ! any date suffices
     write(lu,'(a,2x,i4)')        'VARS',nc3d+nc2d
     do n=1,nc3d
        var = trim(cvars3d(n))
        if (uv_hyb_ens .and. trim(var)=='sf') var='u'
        if (uv_hyb_ens .and. trim(var)=='vp') var='v'
        write(lu,'(a,1x,2(i4,1x),a)') trim(var),grd_anl%nsig,0,trim(var)
     enddo
     do n=1,nc2d
        write(lu,'(a,1x,2(i4,1x),a)') trim(cvars2d(n)),    1,0,trim(cvars2d(n))
     enddo
     write(lu,'(a)') 'ENDVARS'
     close(lu)
  endif

! clean up
  deallocate(work4_2d)
  deallocate(work4_3d)
  deallocate(work8_2d)
  deallocate(work8_3d)

  return
end subroutine write_spread_dualres

end module m_put_gefs_ensperts
