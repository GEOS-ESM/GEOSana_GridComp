subroutine write_bkgvars_grid(a,b,c,d,fname,mype)
!$$$  subroutine documentation block
!
! subprogram:    write_bkgvars_grid
!
!   prgrmmr:
!
! abstract:  modified routine to write out files to compare spectral computation
!            of horizontal derivatives with the derivatives that are being
!            carried around for the dynamical balance constraint
!
! program history log:
!   2008-03-27  safford -- add subprogram doc block, rm unused vars and uses
!
!   input argument list:
!     mype     - mpi task id
!     a        -
!     b        -
!     c        -
!     d        -
!
!   output argument list:
!
! attributes:
!   language:  f90
!   machine:
!
!$$$
  use kinds, only: r_kind,i_kind,r_single
  use gridmod, only: nlat,nlon,nsig,lat2,lon2
  use guess_grids, only: get_ref_gesprs
  use mpeu_util, only: get_lun => luavail
  implicit none

  integer(i_kind)                       ,intent(in   ) :: mype
  character(len=*)                      ,intent(in   ) :: fname

  real(r_kind),dimension(lat2,lon2,nsig),intent(in   ) :: a,b,c
  real(r_kind),dimension(lat2,lon2)     ,intent(in   ) :: d


  real(r_kind),dimension(nlat,nlon,nsig):: ag,bg,cg
  real(r_kind),dimension(nlat,nlon):: dg

  real(r_single),dimension(nlon,nlat,nsig):: a4,b4,c4
  real(r_single),dimension(nlon,nlat):: d4

  character(len=80) :: grdfile
  integer(i_kind) iret,i,j,k,lu
  real(r_kind),dimension(nsig+1)::prs

! gather stuff to processor 0
  do k=1,nsig
     call gather_stuff2(a(1,1,k),ag(1,1,k),mype,0)
     call gather_stuff2(b(1,1,k),bg(1,1,k),mype,0)
     call gather_stuff2(c(1,1,k),cg(1,1,k),mype,0)
  end do
  call gather_stuff2(d,dg,mype,0)

  if (mype==0) then
     write(6,*) 'WRITE OUT NEW VARIANCES'
! load single precision arrays
     do k=1,nsig
        do j=1,nlon
           do i=1,nlat
              a4(j,i,k)=ag(i,j,k)
              b4(j,i,k)=bg(i,j,k)
              c4(j,i,k)=cg(i,j,k)
           end do
        end do
     end do
     do j=1,nlon
        do i=1,nlat
           d4(j,i)=dg(i,j)
        end do
     end do

! Create byte-addressable binary file for grads
     grdfile = trim(fname)//'.grd' 
     call baopenwt(22,trim(grdfile),iret)
     call wryte(22,4*nlat*nlon*nsig,a4)
     call wryte(22,4*nlat*nlon*nsig,b4)
     call wryte(22,4*nlat*nlon*nsig,c4)
     call wryte(22,4*nlat*nlon,d4)
     call baclose(22,iret)

     call get_ref_gesprs(prs)

!    Now create corresponding grads table file
     lu=get_lun()
     open(lu,file=trim(fname)//'.ctl',form='formatted')
     write(lu,'(2a)') 'DSET  ^', trim(grdfile)
     write(lu,'(2a)') 'TITLE ', 'gsi berror variances'
     write(lu,'(a,2x,e13.6)') 'UNDEF', 1.E+15 ! any other preference for this?
     write(lu,'(a,2x,i4,2x,a,2x,f5.1,2x,f9.6)') 'XDEF',nlon, 'LINEAR',   0.0, 360./nlon
     write(lu,'(a,2x,i4,2x,a,2x,f5.1,2x,f9.6)') 'YDEF',nlat, 'LINEAR', -90.0, 180./(nlat-1.)
     write(lu,'(a,2x,i4,2x,a,1x,f8.3)')    'ZDEF',nsig, 'LEVELS', prs(1) ! prs is in mb
     do k=2,nsig  ! grads tends not to like a long line of pressures, thus split
        write(lu,'(1x,f8.3)')   prs(k) ! prs is in mb
     enddo
     write(lu,'(a,2x,i4,2x,a)')   'TDEF', 1, 'LINEAR 12:00Z04JUL1776 6hr' ! any date suffices
     write(lu,'(a,2x,i4)')        'VARS', 4
     write(lu,'(a,1x,2(i4,1x),a)') 'a',nsig,0,'a'
     write(lu,'(a,1x,2(i4,1x),a)') 'b',nsig,0,'b'
     write(lu,'(a,1x,2(i4,1x),a)') 'c',nsig,0,'c'
     write(lu,'(a,1x,2(i4,1x),a)') 'd',   1,0,'d'
     write(lu,'(a)') 'ENDVARS'
     close(lu)

  end if ! mype=0
   
  return
end subroutine write_bkgvars_grid

subroutine write_bkgvars2_grid
!$$$  subroutine documentation block
!
! subprogram:    write_bkgvars2_grid
!
!   prgrmmr:
!
! abstract:  modified routine to write out files to compare spectral computation
!            of horizontal derivatives with the derivatives that are being
!            carried around for the dynamical balance constraint
!
! program history log:
!   2008-03-27  safford -- add subprogram doc block, rm unused vars and uses
!   2010-06-18  todling -- generalized to show all variances; create ctl
!   2010-10-20  pagowski - add cmaq
!   2017-03-23  Hu      - add code to use hybrid vertical coodinate in WRF MASS
!                           core
!   2018-02-15  wu      - add code for fv3_regional
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language:  f90
!   machine:
!
!$$$
  use kinds, only: r_kind,i_kind,r_single
  use mpimod, only: mype
  use constants, only: ten
  use guess_grids, only: get_ref_gesprs
  use gridmod, only: nlat,nlon,nsig
  use control_vectors, only: nc3d,nc2d,mvars
  use control_vectors, only: cvars3d,cvars2d,cvarsmd
  use berror, only: dssv,dssvs
  use file_utility, only : get_lun
  implicit none

  character(255):: grdfile

  real(r_kind),dimension(nlat,nlon,nsig,nc3d):: ag
  real(r_kind),dimension(nlat,nlon,nc2d+mvars):: dg

  real(r_single),dimension(nlon,nlat,nsig,nc3d):: a4
  real(r_single),dimension(nlon,nlat,nc2d+mvars):: d4

  real(r_kind)   ,dimension(nsig+1)::prs
  integer(i_kind) ncfggg,iret,lu,i,j,k,n

! gather stuff to processor 0
  do n=1,nc3d
     do k=1,nsig
        call gather_stuff2(dssv(1,1,k,n),ag(1,1,k,n),mype,0)
     end do
  end do
  do n=1,nc2d
     call gather_stuff2(dssvs(1,1,n),dg(1,1,n),mype,0)
  end do
  do n=1,mvars
     call gather_stuff2(dssvs(1,1,nc2d+n),dg(1,1,nc2d+n),mype,0)
  end do

  call get_ref_gesprs(prs)

  if (mype==0) then
     write(6,*) 'WRITE OUT NEW VARIANCES'
!    Load single precision arrays
     do n=1,nc3d
        do k=1,nsig
           do j=1,nlon
              do i=1,nlat
                 a4(j,i,k,n)=ag(i,j,k,n)
              end do
           end do
        end do
     end do
     do n=1,nc2d+mvars
        do j=1,nlon
           do i=1,nlat
              d4(j,i,n)=dg(i,j,n)
           end do
        end do
     end do

!    Create byte-addressable binary file for grads
     grdfile='bkgvar_smooth.grd'
     ncfggg=len_trim(grdfile)
     lu=get_lun()
     call baopenwt(lu,grdfile(1:ncfggg),iret)
!    Loop over 3d-variances
     do n=1,nc3d
        call wryte(lu,4*nlat*nlon*nsig,a4(1,1,1,n))
     enddo
!    Loop over 2d-variances
     do n=1,nc2d+mvars
        call wryte(lu,4*nlat*nlon,d4(1,1,n))
     enddo
     call baclose(lu,iret)

!    Now create corresponding grads table file
     lu=get_lun()
     open(lu,file='bkgvar_smooth.ctl',form='formatted')
     write(lu,'(2a)') 'DSET  ^', trim(grdfile)
     write(lu,'(2a)') 'TITLE ', 'gsi berror variances'
     write(lu,'(a,2x,e13.6)') 'UNDEF', 1.E+15 ! any other preference for this?
     write(lu,'(a,2x,i4,2x,a,2x,f5.1,2x,f9.6)') 'XDEF',nlon, 'LINEAR',   0.0, 360./nlon
     write(lu,'(a,2x,i4,2x,a,2x,f5.1,2x,f9.6)') 'YDEF',nlat, 'LINEAR', -90.0, 180./(nlat-1.)
     write(lu,'(a,2x,i4,2x,a,100(1x,f10.5))')      'ZDEF',nsig, 'LEVELS', prs(1:nsig) ! prs is in cbar (convert to mb)
     write(lu,'(a,2x,i4,2x,a)')   'TDEF', 1, 'LINEAR 12:00Z04JUL1776 6hr' ! any date suffices
     write(lu,'(a,2x,i4)')        'VARS',nc3d+nc2d+mvars
     do n=1,nc3d
        write(lu,'(a,1x,2(i4,1x),a)') trim(cvars3d(n)),nsig,0,trim(cvars3d(n))
     enddo
     do n=1,nc2d
        write(lu,'(a,1x,2(i4,1x),a)') trim(cvars2d(n)),   1,0,trim(cvars2d(n))
     enddo
     do n=1,mvars
        write(lu,'(a,1x,2(i4,1x),a)') trim(cvarsmd(n)),   1,0,trim(cvarsmd(n))
     enddo
     write(lu,'(a)') 'ENDVARS'
     close(lu)

  end if ! mype=0
   
  return
end subroutine write_bkgvars2_grid

