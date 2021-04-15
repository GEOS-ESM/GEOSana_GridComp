module ncepgfs_io
!$$$ module documentation block
!           .      .    .                                       .
! module:   ncepgfs_io
!   prgmmr: treadon     org: np23                date: 2006-01-10
!
! abstract: This module contains routines which handle input/output
!           operations for NCEP GFS atmospheric and surface files.
!
! program history log:
!   2006-01-10 treadon
!   2009-08-26 li      - add write_gfs_sfc_nst,read_gfsnst, write_gfsnst
!   2010-02-20 parrish - make sigio_cnvtdv8 public so can be accessed by general_read_gfsatm, when
!                          reading in gefs sigma files at resolution different from analysis.
!   2010-03-31 treadon - add read_gfs, use sp_a and sp_b
!   2010-05-19 todling - add read_gfs_chem
!   2011-04-08 li      - (1) add integer nst_gsi to control the mode of NSST 
!                      - (2) add subroutine write_gfs_sfc_nst to save sfc and nst files
!
! Subroutines Included:
!   sub read_gfs          - driver to read ncep gfs atmospheric ("sigma") files
!   sub read_gfssfc       - read ncep gfs surface file, scatter on grid to 
!                           analysis subdomains
!   sub sfc_interpolate   - interpolate from gfs atm grid to gfs sfc grid
!   sub write_gfs         - driver to write ncep gfs atmospheric and surface
!                           analysis files
!   sub write_gfssfc      - gather/write on grid ncep surface analysis file
!   sub read_gfsnst       - read ncep nst file, scatter on grid to
!                           analysis subdomains
!   sub write_gfs_sfc_nst - gather/write on grid ncep surface & nst analysis file
!
! Variable Definitions:
!   none
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  implicit none

  private
  public read_sfc
  public read_gfs
  public read_gfs_chem
  public read_gfssfc
  public read_gfsnst
  public write_gfs
  public write_gfs_sfc_nst
  public sfc_interpolate
  public sigio_cnvtdv8

contains

  subroutine read_gfs(mype)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    read_gfs
!
!   prgrmmr:
!
! abstract:
!
! program history log:
!   2010-03-31  treadon - create routine
!   2011-05-01  todling - cwmr no longer in guess-grids; use metguess bundle now
!   2011-10-01  mkim    - add calculation of hydrometeor mixing ratio from total condensate (cw)  
!   2011-11-01  eliu    - add call to set_cloud_lower_bound (qcmin) 
!   2011-11-01  eliu    - move then calculation of hydrometeor mixing ratio from total condensate to cloud_efr;
!                         rearrange Min-Jeong's code  
!   2013-10-19  todling - update cloud_efr module name
!   2013-10-29  todling - revisit write to allow skipping vars not in MetGuess
!
!   input argument list:
!     mype               - mpi task id
!
!   output argument list:
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

    use kinds, only: i_kind,r_kind
    use gridmod, only: hires_b,sp_a,grd_a,jcap_b,nlon,nlat,lat2,lon2,nsig
    use guess_grids, only: ifilesig,nfldsig 
    use gsi_metguess_mod, only: gsi_metguess_bundle
    use gsi_bundlemod, only: gsi_bundlegetpointer
    use mpeu_util, only: die
    use cloud_efr_mod, only: cloud_calc_gfs,set_cloud_lower_bound    
    use gsi_io, only: mype_io
    use general_specmod, only: general_init_spec_vars,general_destroy_spec_vars,spec_vars
    implicit none

    integer(i_kind),intent(in   ) :: mype

    character(24) filename
    logical:: l_cld_derived
    integer(i_kind):: it,nlon_b
    integer(i_kind):: iret,iret_ql,iret_qi,istatus 

    real(r_kind),dimension(lat2,lon2  ):: aux_ps
    real(r_kind),dimension(lat2,lon2  ):: aux_z
    real(r_kind),dimension(lat2,lon2,nsig):: aux_u
    real(r_kind),dimension(lat2,lon2,nsig):: aux_v
    real(r_kind),dimension(lat2,lon2,nsig):: aux_vor
    real(r_kind),dimension(lat2,lon2,nsig):: aux_div
    real(r_kind),dimension(lat2,lon2,nsig):: aux_tv
    real(r_kind),dimension(lat2,lon2,nsig):: aux_q
    real(r_kind),dimension(lat2,lon2,nsig):: aux_oz
    real(r_kind),dimension(lat2,lon2,nsig):: aux_cwmr

    real(r_kind),pointer,dimension(:,:  ):: ges_ps_it   => NULL()
    real(r_kind),pointer,dimension(:,:  ):: ges_z_it    => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_u_it    => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_v_it    => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_vor_it  => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_div_it  => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_tv_it   => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_q_it    => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_oz_it   => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_cwmr_it => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_ql_it   => NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_qi_it   => NULL()

    type(spec_vars):: sp_b


!   If needed, initialize for hires_b transforms
    nlon_b=((2*jcap_b+1)/nlon+1)*nlon
    if (nlon_b /= sp_a%imax) then
       hires_b=.true.
       call general_init_spec_vars(sp_b,jcap_b,jcap_b,nlat,nlon_b)
       if (mype==0) &
            write(6,*)'READ_GFS:  allocate and load sp_b with jcap,imax,jmax=',&
            sp_b%jcap,sp_b%imax,sp_b%jmax
    endif

    do it=1,nfldsig

       write(filename,100) ifilesig(it)
100    format('sigf',i2.2)
       if (hires_b) then

!         If hires_b, spectral to grid transform for background
!         uses double FFT.   Need to pass in sp_a and sp_b

          call general_read_gfsatm(grd_a,sp_a,sp_b,filename,mype,.true., &
               aux_z,aux_ps,&
               aux_vor,aux_div,&
               aux_u,aux_v,&
               aux_tv,aux_q,&
               aux_cwmr,aux_oz,iret)

       else

!         Otherwise, use standard transform.  Use sp_a in place of sp_b.

          call general_read_gfsatm(grd_a,sp_a,sp_a,filename,mype,.true., &
               aux_z,aux_ps,&
               aux_vor,aux_div,&
               aux_u,aux_v,&
               aux_tv,aux_q,&
               aux_cwmr,aux_oz,iret)
       endif

!      Set values to actual MetGuess fields
       call set_guess_

       l_cld_derived = associated(ges_cwmr_it).and.&
                       associated(ges_q_it)   .and.&
                       associated(ges_ql_it)  .and.&
                       associated(ges_qi_it)  .and.&
                       associated(ges_tv_it)

!      call set_cloud_lower_bound(ges_cwmr_it)
       if (mype==0) write(6,*)'READ_GFS: l_cld_derived = ', l_cld_derived

       if (l_cld_derived) &            
       call cloud_calc_gfs(ges_ql_it,ges_qi_it,ges_cwmr_it,ges_q_it,ges_tv_it) 

    end do

    if (hires_b) call general_destroy_spec_vars(sp_b)


  contains


  subroutine set_guess_

  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'ps',ges_ps_it  ,istatus) 
  if(istatus==0) ges_ps_it = aux_ps
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'z' ,ges_z_it   ,istatus) 
  if(istatus==0) ges_z_it = aux_z
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'u' ,ges_u_it   ,istatus) 
  if(istatus==0) ges_u_it = aux_u
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'v' ,ges_v_it   ,istatus) 
  if(istatus==0) ges_v_it = aux_v
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'vor',ges_vor_it,istatus) 
  if(istatus==0) ges_vor_it = aux_vor
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'div',ges_div_it,istatus) 
  if(istatus==0) ges_div_it = aux_div
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'tv',ges_tv_it  ,istatus) 
  if(istatus==0) ges_tv_it = aux_tv
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'q' ,ges_q_it   ,istatus) 
  if(istatus==0) ges_q_it = aux_q
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'oz',ges_oz_it  ,istatus) 
  if(istatus==0) ges_oz_it = aux_oz
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'cw',ges_cwmr_it,istatus) 
  if(istatus==0) ges_cwmr_it = aux_cwmr
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'ql',ges_ql_it,  iret_ql) 
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'qi',ges_qi_it,  iret_qi)           
  if (iret_ql/=0) then 
     if (mype==0) write(6,*)'READ_GFS: cannot get pointer to ql,iret_ql= ',iret_ql 
  endif
  if (iret_qi/=0) then 
     if (mype==0) write(6,*)'READ_GFS: cannot get pointer to qi,iret_qi= ',iret_qi 
  endif

  end subroutine set_guess_

  end subroutine read_gfs

  subroutine read_gfs_chem (iyear, month,idd )
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    read_gfs_chem
!
!   prgrmmr: todling
!
! abstract: fills chemguess_bundle with GFS chemistry. 
!
! remarks: 
!    1. Right now, only CO2 is done and even this is treated 
!        as constant througout the assimialation window.
!    2. iyear and month could come from obsmod, but logically
!       this program should never depend on obsmod
! 
!
! program history log:
!   2010-04-15  hou - Initial code
!   2010-05-19  todling - Port Hou's code from compute_derived(!)
!                         into this module and linked with the chemguess_bundle
!   2011-02-01  r. yang - proper initialization of prsi
!   2011-05-24  yang    - add idd for time interpolation of co2 field
!   2011-06-29  todling - no explict reference to internal bundle arrays
!   2013-11-08  todling - revisit check for present of GHG in chem-bundle
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

    use kinds, only: i_kind,r_kind
    use mpimod, only: mype
    use gridmod, only: lat2,lon2,nsig,nlat,rlats,istart
    use ncepgfs_ghg, only: read_gfsco2,read_ch4n2oco
    use guess_grids, only: nfldsig,ntguessig
    use gsi_bundlemod, only: gsi_bundlegetpointer
    use gsi_chemguess_mod, only: gsi_chemguess_bundle
    use gsi_chemguess_mod, only: gsi_chemguess_get



    implicit none

!   Declared argument list
    integer(i_kind), intent(in):: iyear
    integer(i_kind), intent(in):: month
    integer(i_kind), intent(in):: idd

!   Declare local variables
    character(len=*),parameter :: myname='read_gfs_chem'
    integer(i_kind)            :: i,j,n,ier
    integer(i_kind)            :: ico24crtm,ich44crtm,in2o4crtm,ico4crtm
    character(len=3) :: char_ghg
    real(r_kind),dimension(lat2):: xlats
    real(r_kind),pointer,dimension(:,:,:)::p_co2=>NULL()
    real(r_kind),pointer,dimension(:,:,:)::p_ch4=>NULL()
    real(r_kind),pointer,dimension(:,:,:)::p_n2o=>NULL()
    real(r_kind),pointer,dimension(:,:,:)::p_co=>NULL()
    real(r_kind),pointer,dimension(:,:,:)::ptr3d_co2=>NULL()
    real(r_kind),pointer,dimension(:,:,:)::ptr3d_ch4=>NULL()
    real(r_kind),pointer,dimension(:,:,:)::ptr3d_n2o=>NULL()
    real(r_kind),pointer,dimension(:,:,:)::ptr3d_co=>NULL()

    if(.not.associated(gsi_chemguess_bundle)) return

!   Get subdomain latitude array
    j = mype + 1
    do i = 1, lat2
       n = min(max(1, istart(j)+i-2), nlat)
       xlats(i) = rlats(n)
    enddo
!!NOTE: NEED TO CHANGE THIS BLOCK, THE CHECK AND READ OF TRACE GASES ARE HARDWIRED !!!!!!
!!      WILL CHANGE THE CODE FOLLOWING WHAT I DID IN crtm_interface.f90            !!!!!!

! check whether CO2 exist
    call gsi_bundlegetpointer(gsi_chemguess_bundle(1),'co2',p_co2,ier)
    if (associated(p_co2)) then
       call gsi_chemguess_get ( 'i4crtm::co2', ico24crtm, ier )
       if (ico24crtm >= 0 ) then
          call read_gfsco2 (iyear,month,idd,ico24crtm,xlats,&
                          lat2,lon2,nsig,mype,  &
                          p_co2 )
! Approximation: assign three time slots (nfldsig) of ghg with same values
          do n=2,nfldsig
             call gsi_bundlegetpointer(gsi_chemguess_bundle(n),'co2',ptr3d_co2,ier)
             ptr3d_co2 = p_co2
          enddo
          char_ghg='co2'
! take comment out for printing out the interpolated tracer gas fields.
!        call write_ghg_grid (ptr3d_co2,char_ghg,mype)
       endif
    endif ! <co2>

! check whether CH4 data exist
    call gsi_bundlegetpointer(gsi_chemguess_bundle(1),'ch4',p_ch4,ier)
    if (associated(p_ch4)) then
       call gsi_chemguess_get ( 'i4crtm::ch4', ich44crtm, ier )
       if (ich44crtm > 0 ) then
          char_ghg='ch4'
          call read_ch4n2oco (iyear,month,idd,char_ghg,xlats,&
                          lat2,lon2,nsig,mype,  &
                          p_ch4 )
          do n=2,nfldsig
             call gsi_bundlegetpointer(gsi_chemguess_bundle(n),'ch4',ptr3d_ch4,ier)
             ptr3d_ch4 = p_ch4
          enddo
! take comment out for printing out the interpolated tracer gas fields.
!         call write_ghg_grid (ptr3d_ch4,char_ghg,mype)
       endif
    endif ! <ch4>

! check whether N2O data exist
    call gsi_bundlegetpointer(gsi_chemguess_bundle(1),'n2o',p_n2o,ier)
    if (associated(p_n2o)) then
       call gsi_chemguess_get ( 'i4crtm::n2o', in2o4crtm, ier )
       if (in2o4crtm > 0 ) then
          char_ghg='n2o'
          call read_ch4n2oco (iyear,month,idd,char_ghg,xlats,&
                          lat2,lon2,nsig,mype,  &
                          p_n2o )
          do n=2,nfldsig
             call gsi_bundlegetpointer(gsi_chemguess_bundle(n),'n2o',ptr3d_n2o,ier)
             ptr3d_n2o = p_n2o
          enddo
! take comment out for printing out the interpolated tracer gas fields.
!        call write_ghg_grid (ptr3d_n2o,char_ghg,mype)
       endif
    endif ! <n2o>

! check whether CO data exist
    call gsi_bundlegetpointer(gsi_chemguess_bundle(1),'co',p_co,ier)
    if (associated(p_co)) then
       call gsi_chemguess_get ( 'i4crtm::co', ico4crtm, ier )
       if (ico4crtm > 0 ) then
          char_ghg='co'
          call read_ch4n2oco ( iyear,month,idd,char_ghg,xlats,&
                          lat2,lon2,nsig,mype,  &
                          p_co )
          do n=2,nfldsig
             call gsi_bundlegetpointer(gsi_chemguess_bundle(n),'co',ptr3d_co,ier)
             ptr3d_co = p_co
          enddo
! take comment out for printing out the interpolated tracer gas fields.
!        call write_ghg_grid (ptr3d_co,char_ghg,mype)
       endif
    endif ! <co>
  end subroutine read_gfs_chem
subroutine write_ghg_grid(a,char_ghg,mype)
!$$$  subroutine documentation block
!
! subprogram:    write_ghg_grid
!
!   prgrmmr:  yang: follow write_bkgvars_grid
!
!
!   input argument list:
!     mype     - mpi task id
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
  use file_utility, only : get_lun
  implicit none

  integer(i_kind)                       ,intent(in   ) :: mype

  real(r_kind),dimension(lat2,lon2,nsig),intent(in   ) :: a
  character(len=3),intent(in) :: char_ghg

  character(255):: grdfile

  real(r_kind),dimension(nlat,nlon,nsig):: ag

  real(r_single),dimension(nlon,nlat,nsig):: a4
  integer(i_kind) ncfggg,iret,i,j,k,lu

! gather stuff to processor 0
  do k=1,nsig
     call gather_stuff2(a(1,1,k),ag(1,1,k),mype,0)
  end do
  if (mype==0) then
     write(6,*) 'WRITE OUT INTERPOLATED',char_ghg
! load single precision arrays
     do k=1,nsig
        do j=1,nlon
           do i=1,nlat
              a4(j,i,k)=ag(i,j,k)
           end do
        end do
     end do

! Create byte-addressable binary file for grads
     grdfile=trim(char_ghg)//'clim_grd'
     ncfggg=len_trim(grdfile)
     lu=get_lun()
     call baopenwt(lu,grdfile(1:ncfggg),iret)
     call wryte(lu,4*nlat*nlon*nsig,a4)
     call baclose(lu,iret)
  end if

  return
end subroutine write_ghg_grid

  subroutine read_sfc(filename,nsfc,fact10,sfct,sno,veg_type,veg_frac,soil_type,soil_temp, &
           soil_moi,sfc_rough,terrain,isli)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    read_sfc
!
!   prgrmmr: whitaker
!
! abstract: read a ncep GFS surface file on a specified task,
!           broadcast data to other tasks.
!
! program history log:
!   2012-01-24  whitaker - create routine
!
!   input argument list:
!     lunges             - unit number to use for IO
!     filename           - gfs surface file to read
!
!   output argument list:
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block
    ! read data from sfc file on a single task, bcast data to other tasks.
    use sfcio_module, only: sfcio_srohdc,sfcio_head,sfcio_data,sfcio_intkind
    use sfcio_module, only: sfcio_axdata
    use kinds, only: i_kind,r_single,r_kind
    use gridmod, only: nlat_sfc,nlon_sfc
    character(*),intent(in) :: filename
    integer(i_kind), intent(in) :: nsfc
    integer(i_kind),dimension(nlat_sfc,nlon_sfc),intent(  out) :: isli
    real(r_kind)   ,dimension(nlat_sfc,nlon_sfc),intent(  out) :: fact10,sfct,sno,&
         veg_type,veg_frac,soil_type,soil_temp,soil_moi,sfc_rough,terrain
    integer(i_kind) :: latb,lonb
    integer(i_kind) :: iret,n,i,j
    type(sfcio_head) :: sfc_head
    type(sfcio_data) :: sfc_data
    real(r_kind),allocatable,dimension(:,:):: outtmp
    integer(sfcio_intkind):: lunges = 11
! read a surface file on the task
    call sfcio_srohdc(lunges,filename,sfc_head,sfc_data,iret)
!   Check for possible problems
    if (iret /= 0) then
       write(6,*)'READ_SFC:  ***ERROR*** problem reading ',filename,&
            ', iret=',iret
       call sfcio_axdata(sfc_data,iret)
       call stop2(80)
    endif
    lonb = sfc_head%lonb
    latb = sfc_head%latb
    if ( (latb /= nlat_sfc-2) .or. (lonb /= nlon_sfc) ) then
         write(6,*)'READ_GFSSFC:  ***ERROR*** inconsistent grid dimensions.  ',&
              ', nlon,nlat-2=',nlon_sfc,nlat_sfc-2,' -vs- sfc file lonb,latb=',&
                 lonb,latb
         call sfcio_axdata(sfc_data,iret)
         call stop2(80)
    endif
!$omp parallel do private(n,i,j,outtmp)
    do n=1,nsfc
     if(n == 1)then                                  !skin temperature

       call tran_gfssfc(sfc_data%tsea,sfct,lonb,latb)                                 
     else if(n == 2) then                            ! soil moisture

       call tran_gfssfc(sfc_data%smc(1:lonb,1:latb,1),soil_moi,lonb,latb)  

     else if(n == 3) then                            ! snow depth

       call tran_gfssfc(sfc_data%sheleg,sno,lonb,latb)        

     else if(n == 4) then                            ! soil temperature

       call tran_gfssfc(sfc_data%stc(1:lonb,1:latb,1),soil_temp,lonb,latb)  

     else if(n == 5) then                            ! sea/land/ice mask

       allocate(outtmp(latb+2,lonb))
       call tran_gfssfc(sfc_data%slmsk,outtmp,lonb,latb)                       
       do j=1,lonb
         do i=1,latb+2
            isli(i,j) = nint(outtmp(i,j))
         end do
       end do
       deallocate(outtmp)

     else if(n == 6) then                             ! vegetation cover

       call tran_gfssfc(sfc_data%vfrac,veg_frac,lonb,latb)                       
     else if(n == 7) then                             ! 10m wind factor

       call tran_gfssfc(sfc_data%f10m,fact10,lonb,latb)                           
     else if(n == 8) then                             ! vegetation type

       call tran_gfssfc(sfc_data%vtype,veg_type,lonb,latb)            

     else if(n == 9) then                             ! soil type

       call tran_gfssfc(sfc_data%stype,soil_type,lonb,latb)                     

     else if(n == 10) then                            ! surface roughness length (cm)

       call tran_gfssfc(sfc_data%zorl,sfc_rough,lonb,latb)            

     else if(n == 11) then                            ! terrain

       call tran_gfssfc(sfc_data%orog,terrain,lonb,latb)            

     end if


!      End of loop over data records
    end do
    deallocate(&
       sfc_data%tsea,&
       sfc_data%smc,&
       sfc_data%sheleg,&
       sfc_data%stc,&
       sfc_data%slmsk,&
       sfc_data%zorl,&
       sfc_data%vfrac,&
       sfc_data%f10m,&
       sfc_data%vtype,&
       sfc_data%stype,&
       sfc_data%orog)
    call sfcio_axdata(sfc_data,iret)

!   Print date/time stamp
    write(6,700) latb,lonb,sfc_head%fhour,sfc_head%idate
700 format('READ_GFSSFC:  ges read/scatter, nlat,nlon=',&
         2i6,', hour=',f10.1,', idate=',4i5)
  end subroutine read_sfc


  subroutine read_gfssfc(filename,iope,mype,fact10,sfct,sno,veg_type,&
       veg_frac,soil_type,soil_temp,soil_moi,isli,sfc_rough,terrain)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_gfssfc     read gfs surface file
!   prgmmr: treadon          org: np23                date: 2003-04-10
!
! abstract: read gfs surface file
!
! program history log:
!   2003-04-10  treadon
!   2004-05-18  kleist, add global isli & documentation
!   2004-09-07  treadon fix mpi bug when npe > nsfc
!   2005-01-27  treadon - rewrite to make use of sfcio module
!   2005-03-07  todling - die gracefully when return error from sfcio
!   2006-09-28  treadon - pull out surface roughness
!   2008-05-28  safford - rm unused vars
!   2009-01-12  gayno   - add read of terrain height
!
!   input argument list:
!     filename - name of surface guess file
!     iope     - mpi task handling i/o
!     mype     - mpi task id
!
!   output argument list:
!     fact10    - 10 meter wind factor
!     sfct      - surface temperature (skin temp)
!     sno       - snow depth
!     veg_type  - vegetation type
!     veg_frac  - vegetation fraction
!     soil_type - soil type
!     soil_temp - soil temperature of first layer
!     soil_moi  - soil moisture of first layer
!     isli      - sea/land/ice mask (subdomain)
!     isli_g    - global sea/land/ice mask
!     sfc_rough - surface roughness
!     terrain   - terrain height
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
    use kinds, only: r_kind,i_kind
    use gridmod, only: nlat_sfc,nlon_sfc
    use mpimod, only: mpi_itype,mpi_rtype,mpi_comm_world
    use constants, only: zero
    implicit none

!   Declare passed variables
    character(*)                                ,intent(in   ) :: filename
    integer(i_kind)                             ,intent(in   ) :: iope
    integer(i_kind)                             ,intent(in   ) :: mype
    integer(i_kind),dimension(nlat_sfc,nlon_sfc),intent(  out) :: isli
    real(r_kind)   ,dimension(nlat_sfc,nlon_sfc),intent(  out) :: fact10,sfct,sno,&
         veg_type,veg_frac,soil_type,soil_temp,soil_moi,sfc_rough,terrain

!   Declare local parameters
    integer(i_kind),parameter:: nsfc=11

!   Declare local variables
    integer(i_kind):: iret,npts


!-----------------------------------------------------------------------------
!   Read surface file on processor iope
    if(mype == iope)then
      call read_sfc(filename,nsfc,fact10,sfct,sno,veg_type,veg_frac,soil_type,soil_temp, &
           soil_moi,sfc_rough,terrain,isli)
    end if

!     Load onto all processors

    npts=nlat_sfc*nlon_sfc

    call mpi_bcast(sfct,npts,mpi_rtype,iope,mpi_comm_world,iret)
    call mpi_bcast(fact10,npts,mpi_rtype,iope,mpi_comm_world,iret)
    call mpi_bcast(sno,npts,mpi_rtype,iope,mpi_comm_world,iret)
    call mpi_bcast(sfc_rough,npts,mpi_rtype,iope,mpi_comm_world,iret)
    call mpi_bcast(terrain,npts,mpi_rtype,iope,mpi_comm_world,iret)
    call mpi_bcast(isli,npts,mpi_itype,iope,mpi_comm_world,iret)
    call mpi_bcast(veg_type,npts,mpi_rtype,iope,mpi_comm_world,iret)
    call mpi_bcast(veg_frac,npts,mpi_rtype,iope,mpi_comm_world,iret)
    call mpi_bcast(soil_type,npts,mpi_rtype,iope,mpi_comm_world,iret)
    call mpi_bcast(soil_temp,npts,mpi_rtype,iope,mpi_comm_world,iret)
    call mpi_bcast(soil_moi,npts,mpi_rtype,iope,mpi_comm_world,iret)

    return
  end subroutine read_gfssfc

subroutine tran_gfssfc(ain,aout,lonb,latb)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    tran_gfssfc     transform gfs surface file to analysis grid
!   prgmmr: derber          org: np2                date: 2003-04-10
!
! abstract: transform gfs surface file to analysis grid
!
! program history log:
!   2012-31-38  derber  - initial routine
!
!   input argument list:
!     ain      - input surface record on processor iope
!     lonb     - input number of longitudes
!     latb     - input number of latitudes
!
!   output argument list:
!     aout     - output transposed surface record
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
    use kinds, only: r_kind,i_kind
    use constants, only: zero
    use sfcio_module, only: sfcio_realkind
    implicit none

!   Declare passed variables
    integer(i_kind)                  ,intent(in ) :: lonb,latb
    real(sfcio_realkind),dimension(lonb,latb),intent(in ) :: ain
    real(r_kind),dimension(latb+2,lonb),intent(out) :: aout

!   Declare local variables
    integer(i_kind) i,j
    real(r_kind) sumn,sums
!   of surface guess array
    sumn = zero
    sums = zero
    do i=1,lonb
       sumn = ain(i,1)    + sumn
       sums = ain(i,latb) + sums
    end do
    sumn = sumn/float(lonb)
    sums = sums/float(lonb)

!    Transfer from local work array to surface guess array
    do j = 1,lonb
       aout(1,j)=sums
       do i=2,latb+1
          aout(i,j) = ain(j,latb+2-i)
       end do
       aout(latb+2,j)=sumn
    end do

    return
    end subroutine tran_gfssfc


  subroutine read_gfsnst(filename,mype,tref,dt_cool,z_c,dt_warm,z_w,c_0,c_d,w_0,w_d)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_gfsnst     read gfs nst guess file (quadratic Gaussin grids) without scattering to tasks
!   prgmmr: li          org: np23                date: 2009-08-26
!
! abstract: read nst file
!
! program history log:
!
!   input argument list:
!     filename - name of nst guess file
!     mype     - mpi task id
!
!  output argument list:
!  tref     (:,:)                        ! oceanic foundation temperature
!  dt_cool  (:,:)                        ! sub-layer cooling amount at sub-skin layer
!  z_c      (:,:)                        ! depth of sub-layer cooling layer
!  dt_warm  (:,:)                        ! diurnal warming amount at sea surface (skin layer)
!  z_w      (:,:)                        ! depth of diurnal warming layer
!  c_0      (:,:)                        ! coefficient to calculate d(Tz)/d(tr) in dimensionless
!  c_d      (:,:)                        ! coefficient to calculate d(Tz)/d(tr) in m^-1
!  w_0      (:,:)                        ! coefficient to calculate d(Tz)/d(tr) in dimensionless
!  w_d      (:,:)                        ! coefficient to calculate d(Tz)/d(tr) in m^-1

! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
    use kinds, only: r_kind,i_kind
    use gridmod, only: itotsub,nlat_sfc,nlon_sfc
    use nstio_module, only: nstio_srohdc,nstio_head,nstio_data

    use constants, only: zero,two
    implicit none

!   Declare passed variables
    character(*),intent(in):: filename
    integer(i_kind),intent(in):: mype
    real(r_kind),dimension(nlat_sfc,nlon_sfc),intent(out):: tref,dt_cool,z_c,dt_warm,z_w,c_0,c_d,w_0,w_d
!   Declare local parameters
    integer(i_kind):: lun_nst = 13
    integer(i_kind),parameter:: n_nst=9
!   Declare local variables
    integer(i_kind) i,j,k,latb,lonb
    integer(i_kind):: irets
    real(r_kind) sumn,sums
    real(r_kind),allocatable,dimension(:,:,:):: work,nstges

    type(nstio_head):: nst_head
    type(nstio_data):: nst_data
!   Read nst file
    call nstio_srohdc(lun_nst,trim(filename),nst_head,nst_data,irets)

!   Check for possible problems
    if (irets /= 0) then
       write(6,*)'READ_GFSNST:  ***ERROR*** problem reading ',filename,&
            ', irets=',irets
       call stop2(80)
    endif
    latb=nst_head%latb                   ! e.g. 576     for T382
    lonb=nst_head%lonb                   ! e.g. 1152    for T382
    if ( (latb /= nlat_sfc-2) .or. &
         (lonb /= nlon_sfc) ) then
       write(6,*)'READ_GFSNST:  ***ERROR*** inconsistent grid dimensions.  ',&
            ', nlon_sfc,nlat_sfc-2=',nlon_sfc,nlat_sfc-2,' -vs- nst file lonb,latb=',&
            lonb,latb
       call stop2(80)
    endif

!   Load nst guess fields required in Tr analysis into local work array
    allocate(work(lonb,latb,n_nst),nstges(latb+2,lonb,n_nst))

    work = zero
    work(:,:,1)  = nst_data%tref(:,:)
    work(:,:,2)  = nst_data%dt_cool(:,:)
    work(:,:,3)  = nst_data%z_c(:,:)
    work(:,:,4)  = two*nst_data%xt(:,:)/nst_data%xz(:,:)
    work(:,:,5)  = nst_data%xz(:,:)
    work(:,:,6)  = nst_data%c_0(:,:)
    work(:,:,7)  = nst_data%c_d(:,:)
    work(:,:,8)  = nst_data%w_0(:,:)
    work(:,:,9)  = nst_data%w_d(:,:)

!     Fill nst guess array
      do k=1,n_nst

!        Compute mean for southern- and northern-most rows
!        of surface guess array
         sumn = zero
         sums = zero
         do i=1,lonb
            sumn = work(i,1,k)    + sumn
            sums = work(i,latb,k) + sums
         end do
         sumn = sumn/float(lonb)
         sums = sums/float(lonb)

!        Transfer from local work array to surface guess array
         do j = 1,lonb
            nstges(1,j,k)=sums
            nstges(latb+2,j,k)=sumn
            do i=2,latb+1
              nstges(i,j,k) = work(j,latb+2-i,k)
            end do
          end do

!     End of loop over data records
      end do


!     Deallocate local work arrays
      deallocate(work)
!   Load data into output arrays
    do j=1,lonb
     do i=1,latb+2
       tref(i,j)    = nstges(i,j,1)
       dt_cool(i,j) = nstges(i,j,2)
       z_c(i,j)     = nstges(i,j,3)
       dt_warm(i,j) = nstges(i,j,4)
       z_w(i,j)     = nstges(i,j,5)
       c_0(i,j)     = nstges(i,j,6)
       c_d(i,j)     = nstges(i,j,7)
       w_0(i,j)     = nstges(i,j,8)
       w_d(i,j)     = nstges(i,j,9)
     end do
    end do
    deallocate(nstges)

!   Print date/time stamp
    if(mype==0) then
       write(6,700) latb,lonb,nst_head%fhour,nst_head%idate
700    format('READ_GFSNST:  ges read, nlat,nlon=',&
            2i6,', hour=',f10.1,', idate=',4i5)
    end if

    return
  end subroutine read_gfsnst

  subroutine write_gfs(increment,mype,mype_atm,mype_sfc)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    write_gfs
!
!   prgrmmr:
!
! abstract:
!
! program history log:
!   2006-07-31  kleist - pass ges_ps instead of ges_lnps
!   2006-10-11  treadon - update 10m wind factor in sfc file
!   2008-05-28  safford - rm unused vars, add doc block
!   2008-12-05  todling - adjustment for dsfct time dimension addition
!   2009-08-28  li      - add nst i/o
!   2009-11-28  todling - add increment option (hook-only for now)
!   2010-03-31  treadon - add hires_b, sp_a, and sp_b
!   2011-05-01  todling - cwmr no longer in guess-grids; use metguess bundle now
!   2013-10-19  todling - update cloud_efr module name
!   2013-10-29  todling - revisit write to allow skipping vars not in MetGuess
!
!   input argument list:
!     increment          - when >0 will write increment from increment-index slot
!     mype               - mpi task id
!     mype_atm,mype_sfc  -
!
!   output argument list:
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

    use kinds, only: i_kind,r_kind
    use guess_grids, only: dsfct,isli2
    use guess_grids, only: ntguessig,ntguessfc,ifilesig,nfldsig
    use gridmod, only: hires_b,sp_a,grd_a,jcap_b,nlon,nlat,lat2,lon2,nsig
    use gsi_metguess_mod, only: gsi_metguess_bundle
    use gsi_bundlemod, only: gsi_bundlegetpointer
    use mpeu_util, only: die
    use radinfo, only: nst_gsi
    use constants, only:zero
    use general_specmod, only: general_init_spec_vars,general_destroy_spec_vars,spec_vars
    use gsi_4dvar, only: lwrite4danl

    implicit none

    integer(i_kind),intent(in   ) :: increment
    integer(i_kind),intent(in   ) :: mype,mype_atm,mype_sfc
    character(24):: filename
    integer(i_kind) itoutsig,istatus,iret_write,nlon_b,ntlevs,it

    character(24):: file_sfc,file_nst

    real(r_kind),dimension(lat2,lon2  ):: aux_ps
    real(r_kind),dimension(lat2,lon2  ):: aux_z
    real(r_kind),dimension(lat2,lon2,nsig):: aux_u
    real(r_kind),dimension(lat2,lon2,nsig):: aux_v
    real(r_kind),dimension(lat2,lon2,nsig):: aux_vor
    real(r_kind),dimension(lat2,lon2,nsig):: aux_div
    real(r_kind),dimension(lat2,lon2,nsig):: aux_tv
    real(r_kind),dimension(lat2,lon2,nsig):: aux_q
    real(r_kind),dimension(lat2,lon2,nsig):: aux_oz
    real(r_kind),dimension(lat2,lon2,nsig):: aux_cwmr

    real(r_kind),pointer,dimension(:,:  ):: ges_ps_it  =>NULL()
    real(r_kind),pointer,dimension(:,:  ):: ges_z_it   =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_u_it   =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_v_it   =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_div_it =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_vor_it =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_tv_it  =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_q_it   =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_oz_it  =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_cwmr_it=>NULL()

    type(spec_vars):: sp_b

!   Write atmospheric analysis file
    if (.not.lwrite4danl) then
       ntlevs=1
    else
       ntlevs=nfldsig
    end if

    aux_ps=zero
    aux_z=zero
    aux_u=zero
    aux_v=zero
    aux_vor=zero
    aux_div=zero
    aux_tv=zero
    aux_q=zero
    aux_oz=zero
    aux_cwmr=zero

    do it=1,ntlevs
       if (increment>0) then
          filename='siginc'
          itoutsig=increment
          if(mype==0) write(6,*) 'WRITE_GFS: writing time slot ', itoutsig
       else if (.not.lwrite4danl) then
          filename='siganl'
          itoutsig=ntguessig
          if(mype==0) write(6,*) 'WRITE_GFS: writing single analysis state for F ', itoutsig
       else
          write(filename,100) ifilesig(it)
100       format('siga',i2.2)
          itoutsig=it
          if(mype==0) write(6,*) 'WRITE_GFS: writing full analysis state for F ', itoutsig
       endif

       call set_analysis_(itoutsig)

!   If hires_b, spectral to grid transform for background
!   uses double FFT.   Need to pass in sp_a and sp_b
       nlon_b=((2*jcap_b+1)/nlon+1)*nlon
       if (nlon_b /= sp_a%imax) then
          hires_b=.true.
          call general_init_spec_vars(sp_b,jcap_b,jcap_b,nlat,nlon_b)
          if (mype==0) &
               write(6,*)'WRITE_GFS:  allocate and load sp_b with jcap,imax,jmax=',&
               sp_b%jcap,sp_b%imax,sp_b%jmax

          call general_write_gfsatm(grd_a,sp_a,sp_b,filename,mype,mype_atm, &
               aux_z,aux_ps,&
               aux_vor,aux_div,&
               aux_tv,aux_q,&
               aux_oz,aux_cwmr,it,&
               iret_write)

          call general_destroy_spec_vars(sp_b)

!   Otherwise, use standard transform.  Use sp_a in place of sp_b.
       else
          call general_write_gfsatm(grd_a,sp_a,sp_a,filename,mype,mype_atm, &
               aux_z,aux_ps,&
               aux_vor,aux_div,&
               aux_tv,aux_q,&
               aux_oz,aux_cwmr,it,&
               iret_write)
       endif
    end do ! end do over ntlevs

!   Write surface analysis file
    if (increment>0) then
       filename='sfcinc.gsi'
       call write_gfssfc(filename,mype,mype_sfc,dsfct(1,1,ntguessfc))
    else
      if ( nst_gsi > 0 ) then
        file_sfc = 'sfcanl'
        file_nst = 'nstanl'
        call write_gfs_sfc_nst(file_sfc,file_nst,mype,mype_sfc,dsfct(1,1,ntguessfc))
      else
        filename='sfcanl.gsi'
        call write_gfssfc(filename,mype,mype_sfc,dsfct(1,1,ntguessfc))
      endif
    endif

  contains


  subroutine set_analysis_(it)
  implicit none
  integer(i_kind),intent(in) :: it

  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'ps',ges_ps_it  ,istatus) 
  if(istatus==0) aux_ps = ges_ps_it
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'z' ,ges_z_it   ,istatus) 
  if(istatus==0) aux_z = ges_z_it
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'u' ,ges_u_it   ,istatus) 
  if(istatus==0) aux_u = ges_u_it
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'v' ,ges_v_it   ,istatus) 
  if(istatus==0) aux_v = ges_v_it
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'vor',ges_vor_it,istatus) 
  if(istatus==0) aux_vor = ges_vor_it
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'div',ges_div_it,istatus) 
  if(istatus==0) aux_div = ges_div_it
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'tv',ges_tv_it  ,istatus) 
  if(istatus==0) aux_tv = ges_tv_it
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'q' ,ges_q_it   ,istatus) 
  if(istatus==0) aux_q = ges_q_it
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'oz',ges_oz_it  ,istatus) 
  if(istatus==0) aux_oz = ges_oz_it
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'cw',ges_cwmr_it,istatus) 
  if(istatus==0) aux_cwmr = ges_cwmr_it

  end subroutine set_analysis_

  end subroutine write_gfs

  subroutine write_gfssfc(filename,mype,mype_sfc,dsfct)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    write_gfssfc --- Write surface analysis to file
!
!   prgrmmr:     treadon -  initial version; org: np22
!
! abstract:     This routine writes the updated surface analysis.  At
!               this point (20040615) the only surface field update by 
!               the gsi is the skin temperature.  The current (20040615)
!               GDAS setup does use the updated surface file.  Rather,
!               the output from surface cycle is used as the surface
!               analysis for subsequent GFS runs.
!
!               The routine gathers surface fields from subdomains, 
!               reformats the data records, and then writes each record
!               to the output file.  
!
!               Since the gsi only update the skin temperature, all
!               other surface fields are simply read from the guess
!               surface file and written to the analysis file.
!
!   Structure of GFS surface file  
!       data record  1    label
!       data record  2    date, dimension, version, lons/lat record
!       data record  3    tsf
!       data record  4    soilm(two layers)
!       data record  5    snow
!       data record  6    soilt(two layers)
!       data record  7    tg3
!       data record  8    zor
!       data record  9    cv
!       data record 10    cvb
!       data record 11    cvt
!       data record 12    albedo (four types)
!       data record 13    slimsk
!       data record 14    vegetation cover
!       data record 15    plantr
!       data record 16    f10m
!       data record 17    canopy water content (cnpanl)
!       data record 18    vegetation type
!       data record 19    soil type
!       data record 20    zenith angle dependent vegetation fraction (two types)
!
! program history log:
!   2004-06-15  treadon -  updated documentation
!   2004-07-15  todling -  protex-compliant prologue; added intent/only's
!   2004-12-03  treadon -  replace mpe_igatherv (IBM extension) with
!                          standard mpi_gatherv
!   2005-01-27  treadon - rewrite to make use of sfcio module
!   2005-02-09  kleist  - clean up unit number and filename for updated surface file
!   2005-03-07  todling -  die gracefully when return error from sfcio
!   2005-03-10  treadon - remove iadate from calling list, access via obsmod
!   2006-10-11  treadon - update 10m wind factor in sfc file
!   2008-05-28  safford - rm unused vars
!   2013-10-25  todling - move ltosj/s to comm_vars
!
!   input argument list:
!     filename  - file to open and write to
!     dsfct     - delta skin temperature
!     mype      - mpi task number
!     mype_sfc  - mpi task to write output file
!
!   output argument list:
!
! attributes:
!   language: f90
!   machines: ibm RS/6000 SP; SGI Origin 2000; Compaq HP
!
!$$$ end documentation block

! !USES:
    use kinds, only: r_kind,r_single,i_kind
  
    use mpimod, only: mpi_rtype
    use mpimod, only: mpi_comm_world
    use mpimod, only: ierror
    
    use gridmod, only: nlat,nlon
    use gridmod, only: lat1,lon1
    use gridmod, only: lat2,lon2
    use gridmod, only: iglobal
    use gridmod, only: ijn
    use gridmod, only: displs_g
    use gridmod, only: itotsub
    
    use general_commvars_mod, only: ltosi,ltosj

    use obsmod, only: iadate
    
    use constants, only: zero_single
    
    use sfcio_module, only: sfcio_intkind,sfcio_head,sfcio_data,&
         sfcio_srohdc,sfcio_swohdc,sfcio_axdata
    
    implicit none

! !INPUT PARAMETERS:
    character(*)                     ,intent(in   ) :: filename  ! file to open and write to

    real(r_kind),dimension(lat2,lon2),intent(in   ) :: dsfct   ! delta skin temperature

    integer(i_kind)                  ,intent(in   ) :: mype     ! mpi task number
    integer(i_kind)                  ,intent(in   ) :: mype_sfc ! mpi task to write output file

! !OUTPUT PARAMETERS:

!-------------------------------------------------------------------------

!   Declare local parameters
    character( 6),parameter:: fname_ges='sfcf06'
    integer(sfcio_intkind),parameter:: ioges = 12
    integer(sfcio_intkind),parameter:: ioanl = 52

    real(r_kind),parameter :: houra = zero_single

!   Declare local variables
    integer(sfcio_intkind):: iret
    integer(i_kind) latb,lonb,nlatm2
    integer(i_kind) i,j,ip1,jp1,ilat,ilon,jj,mm1

    real(r_kind),dimension(nlon,nlat):: buffer
    real(r_kind),dimension(lat1,lon1):: sfcsub
    real(r_kind),dimension(nlon,nlat):: grid
    real(r_kind),dimension(max(iglobal,itotsub)):: sfcall
    real(r_kind),allocatable,dimension(:,:):: buffer2

    type(sfcio_head):: head
    type(sfcio_data):: data

  
!*****************************************************************************

!   Initialize local variables
    mm1=mype+1
    nlatm2=nlat-2

!   Gather skin temperature information from all tasks.  
    do j=1,lon1
       jp1 = j+1
       do i=1,lat1
          ip1 = i+1
          sfcsub(i,j)=dsfct(ip1,jp1)
       end do
    end do
    call mpi_gatherv(sfcsub,ijn(mm1),mpi_rtype,&
         sfcall,ijn,displs_g,mpi_rtype,mype_sfc,&
         mpi_comm_world,ierror)

!   Only MPI task mype_sfc writes the surface file.
    if (mype==mype_sfc) then

!      Reorder updated skin temperature to output format
       do i=1,iglobal
          ilon=ltosj(i)
          ilat=ltosi(i)
          grid(ilon,ilat)=sfcall(i)
       end do
       do j=1,nlat
          jj=nlat-j+1
          do i=1,nlon
             buffer(i,j)=grid(i,jj)
          end do
       end do

!      For now, rather than carry around all the surface fields in memory from
!      the read in ingesfc, just read fields from surface file.  Also, for
!      now, only update the 6-hour forecast surface guess file.

!      Read surface guess file
       call sfcio_srohdc(ioges,fname_ges,head,data,iret)
       if (iret /= 0) then
          write(6,*)'WRITE_GFSSFC:  ***ERROR*** problem reading ',fname_ges,&
               ', iret=',iret
          call sfcio_axdata(data,iret)
          call stop2(80)
       endif
       latb=head%latb
       lonb=head%lonb
       allocate(buffer2(lonb,latb))
       if ( (latb /= nlatm2) .or. &
            (lonb /= nlon) ) then
          write(6,*)'WRITE_GFSSFC:  different grid dimensions analysis vs sfc. interpolating sfc temperature  ',&
               ', nlon,nlat-2=',nlon,nlatm2,' -vs- sfc file lonb,latb=',&
               lonb,latb
          call sfc_interpolate(buffer,nlon,nlat,buffer2,lonb,latb)
       else
          do j=1,latb
             do i=1,lonb
                buffer2(i,j)=buffer(i,j+1)
             end do
          end do
       endif

!      Update guess date/time to analysis date/time
       head%fhour = houra       ! forecast hour
       head%idate(1)=iadate(4)  ! hour
       head%idate(2)=iadate(2)  ! month
       head%idate(3)=iadate(3)  ! day
       head%idate(4)=iadate(1)  ! year


       do j=1,latb
          do i=1,lonb
             data%tsea(i,j) = data%tsea(i,j)+buffer2(i,j)
          end do
       end do
       deallocate(buffer2)

!      Write updated information to surface analysis file
       call sfcio_swohdc(ioanl,filename,head,data,iret)


!      Deallocate local work arrays
       call sfcio_axdata(data,iret)

       write(6,100) lonb,latb,houra,iadate(1:4),iret
100    format(' WRITE_GFSSFC:  sfc analysis written  for ',&
            2i6,1x,f4.1,4(i4,1x),' with iret=',i2)

    endif
    
!   End of routine
    return
  end subroutine write_gfssfc

  subroutine write_gfs_sfc_nst(fname_sfc,fname_nst,mype,mype_so,dsfct)

!
! abstract: write both sfc and nst analysis files (nst_gsi dependent)
!   2013-10-25  todling - move ltosj/s to comm_vars
!

!
! !USES:
!
    use kinds, only: r_kind,r_single,i_kind

    use mpimod, only: mpi_rtype,mpi_itype
    use mpimod, only: mpi_comm_world
    use mpimod, only: ierror

    use gridmod, only: nlat_sfc,nlon_sfc
    use gridmod, only: nlat,nlon
    use gridmod, only: lat1,lon1
    use gridmod, only: lat2,lon2
    use gridmod, only: iglobal
    use gridmod, only: ijn
    use gridmod, only: displs_g
    use gridmod, only: itotsub

    use general_commvars_mod, only: ltosi,ltosj

    use obsmod, only: iadate

    use constants, only: zero_single,zero,two

    use guess_grids, only: isli2
    use radinfo, only: nst_gsi
    use sfcio_module, only: sfcio_intkind,sfcio_head,sfcio_data,&
         sfcio_srohdc,sfcio_swohdc,sfcio_axdata

    use nstio_module, only: nstio_intkind,nstio_head,nstio_data,&
         nstio_srohdc,nstio_swohdc,nstio_axdata

    implicit none
!
! !INPUT PARAMETERS:
!
    character(24),intent(in):: fname_sfc,fname_nst
    real(r_kind),dimension(lat2,lon2), intent(in) :: dsfct       ! delta tr temperature
    integer(i_kind),                   intent(in) :: mype        ! mpi task number
    integer(i_kind),                   intent(in) :: mype_so     ! mpi task to write output file
!
! !OUTPUT PARAMETERS:
!

! !DESCRIPTION: This routine writes the sfc & nst analysis files and is nst_gsi dependent.
!               Tr (foundation temperature), instead of skin temperature, is the analysis variable.
!               nst_gsi >  2: Tr analysis is on
!               nst_gsi <= 2: Tr analysis is off
!
!               The routine gathers Tr field from subdomains,
!               reformats the data records, and then writes each record
!               to the output files.
!
!               Since the gsi only update the Tr temperature, all
!               other fields in surface are simply read from the guess
!               files and written to the analysis file.
!
! !REVISION HISTORY:
!
! !REMARKS:
!
!   language: f90
!   machines: ibm RS/6000 SP; SGI Origin 2000; Compaq HP
!
! !AUTHOR:
!
!   2009-08-28  xu li -  initial version; org: np22
!
!EOP
!-------------------------------------------------------------------------

!   Declare local parameters
    character(10),parameter:: fname_ges_sfc ='sfcf06'
    character(6),parameter:: fname_ges_nst ='nstf06'

    integer(sfcio_intkind),parameter:: ioges_sfc = 12
    integer(sfcio_intkind),parameter:: ioanl_sfc = 52

    integer(sfcio_intkind),parameter:: ioges_nst = 13
    integer(sfcio_intkind),parameter:: ioanl_nst = 53

    real(r_kind),parameter :: houra = zero_single



!   Declare local variables
    integer(i_kind):: iret
    integer(i_kind) latb,lonb,nlatm2
    integer(i_kind) i,j,ip1,jp1,ilat,ilon,jj,mm1

    real(r_kind),dimension(nlon,nlat):: buffer
    real(r_kind),dimension(lat1,lon1):: sosub
    real(r_kind),dimension(nlon,nlat):: grid
    real(r_kind),dimension(max(iglobal,itotsub)):: soall
    real(r_kind),allocatable,dimension(:,:):: buffer2

    integer(i_kind),dimension(nlon,nlat):: isli
    integer(i_kind),dimension(lat1,lon1):: isosub
    integer(i_kind),dimension(nlon,nlat):: igrid
    integer(i_kind),dimension(max(iglobal,itotsub)):: isoall

    type(sfcio_head):: head_sfc
    type(sfcio_data):: data_sfc

    type(nstio_head):: head_nst
    type(nstio_data):: data_nst
!*****************************************************************************

!   Initialize local variables
    mm1=mype+1
    nlatm2=nlat-2

!   Gather analysis variable (reference/foundation temperature) information from all tasks.
    do j=1,lon1
       jp1 = j+1
       do i=1,lat1
          ip1 = i+1
          sosub(i,j)=dsfct(ip1,jp1)
       end do
    end do
    call mpi_gatherv(sosub,ijn(mm1),mpi_rtype,&
         soall,ijn,displs_g,mpi_rtype,mype_so ,&
         mpi_comm_world,ierror)

!   Gather land/ice/sea mask information from all tasks.
    do j=1,lon1
       jp1 = j+1
       do i=1,lat1
          ip1 = i+1
          isosub(i,j)=isli2(ip1,jp1)
       end do
    end do
    call mpi_gatherv(isosub,ijn(mm1),mpi_itype,&
         isoall,ijn,displs_g,mpi_itype,mype_so ,&
         mpi_comm_world,ierror)

!   Only MPI task mype_so  writes the surface & nst file.
    if (mype==mype_so) then

!      Record updated skin temperature to output format
       do i=1,iglobal
          ilon=ltosj(i)
          ilat=ltosi(i)
          grid(ilon,ilat)=soall(i)
       end do
       do j=1,nlat
          jj=nlat-j+1
          do i=1,nlon
             buffer(i,j)=grid(i,jj)
          end do
       end do

!      Record updated isli to output format
       do i=1,iglobal
          ilon=ltosj(i)
          ilat=ltosi(i)
          igrid(ilon,ilat)=isoall(i)
       end do
       do j=1,nlat
          jj=nlat-j+1
          do i=1,nlon
             isli(i,j)=igrid(i,jj)
          end do
       end do

!
!      set dsfct to be zero over non-water grids
!
       do j=1,nlat
         do i=1,nlon
           if ( isli(i,j) > 0 ) then
             buffer(i,j)=zero
           endif
         end do
       end do

!      For now, rather than carry around all the sfc and nst fields in memory from
!      the read in read_gfssfc and read_gfsnst, just read fields from sfc & nst file.  Also, for
!      now, only update the 6-hour forecast sfc & nst guess file.

!      Read  nst guess file
       call nstio_srohdc(ioges_nst,fname_ges_nst,head_nst,data_nst,iret)
       if (iret /= 0) then
          write(6,*)'WRITE_GFSNST:  ***ERROR*** problem reading ',fname_ges_nst,&
               ', iret=',iret
          call nstio_axdata(data_nst,iret)
          call stop2(80)
       endif

!      Read surface guess file
       call sfcio_srohdc(ioges_sfc,fname_ges_sfc,head_sfc,data_sfc,iret)
       if (iret /= 0) then
          write(6,*)'WRITE_GFSSFC:  ***ERROR*** problem reading ',fname_ges_sfc,&
               ', iret=',iret
          call sfcio_axdata(data_sfc,iret)
          call stop2(80)
       endif

       if ( head_nst%latb /= head_sfc%latb .or. head_nst%lonb /= head_sfc%lonb ) then
          write(6,*) 'Inconsistent dimension for sfc & nst files. head_nst%latb,head_nst%lonb : ',head_nst%latb,head_nst%lonb, &
                     'head_sfc%latb,head_sfc%lonb : ',head_sfc%latb, head_sfc%lonb
          call stop2(80)
       endif

       latb=head_sfc%latb
       lonb=head_sfc%lonb
       allocate(buffer2(lonb,latb))
       if ( (latb /= nlatm2) .or. &
            (lonb /= nlon) ) then
          write(6,*)'WRITE_GFSSFC:  different grid dimensions analysis vs sfc. interpolating sfc temperature  ',&
               ', nlon,nlat-2=',nlon,nlatm2,' -vs- sfc file lonb,latb=',&
               lonb,latb
          call sfc_interpolate(buffer,nlon,nlat,buffer2,lonb,latb)
       else
          do j=1,latb
            do i=1,lonb
              buffer2(i,j)=buffer(i,j+1)
            end do
          end do
       endif


       do j=1,latb
         do i=1,lonb
           if ( nint(data_sfc%slmsk(i,j)) > 0) then
             buffer2(i,j)=zero
           endif
         end do
       end do

!      Update guess date/time to analysis date/time for surface file
       head_sfc%fhour = houra       ! forecast hour
       head_sfc%idate(1)=iadate(4)  ! hour
       head_sfc%idate(2)=iadate(2)  ! month
       head_sfc%idate(3)=iadate(3)  ! day
       head_sfc%idate(4)=iadate(1)  ! year

!
!      update tsea (in the surface file) When Tr analysis is on
!
       if ( nst_gsi > 2 ) then
         do j=1,latb
            do i=1,lonb
              if ( nint(data_sfc%slmsk(i,j)) == 0) then
                data_sfc%tsea(i,j) = max(data_nst%tref(i,j)+buffer2(i,j)+& 
                     two*data_nst%xt(i,j)/data_nst%xz(i,j)-data_nst%dt_cool(i,j),271.0_r_kind)
              endif
            end do
         end do
       else
         do j=1,latb
            do i=1,lonb
              if ( nint(data_sfc%slmsk(i,j)) == 0) then
                data_nst%tref(i,j) = data_sfc%tsea(i,j)           ! keep tref as tsea before analysis
              endif
              data_sfc%tsea(i,j) = max(data_sfc%tsea(i,j)+buffer2(i,j),271.0_r_kind)
            end do
         end do
       endif

!      Write updated information to surface analysis file
       call sfcio_swohdc(ioanl_sfc,fname_sfc,head_sfc,data_sfc,iret)

       write(6,100) lonb,latb,houra,iadate(1:4),iret
100    format(' WRITE_GFSSFC:  sfc analysis written  for ',&
            2i6,1x,f4.1,4(i4,1x),' with iret=',i2)



!      Update guess date/time to analysis date/time for nst file
       head_nst%fhour = houra       ! forecast hour
       head_nst%idate(1)=iadate(4)  ! hour
       head_nst%idate(2)=iadate(2)  ! month
       head_nst%idate(3)=iadate(3)  ! day
       head_nst%idate(4)=iadate(1)  ! year

!
!      update tref (in the nst file) When Tr analysis is on
!
       if ( nst_gsi > 2 ) then
         do j=1,latb
            do i=1,lonb
               if ( nint(data_nst%slmsk(i,j)) == 0) then
                 data_nst%tref(i,j) = max(data_nst%tref(i,j)+buffer2(i,j),271.0_r_kind)
               else
                 data_nst%tref(i,j) = data_sfc%tsea(i,j)
               endif
            end do
         end do
       endif

!      Write updated information to nst analysis file
       call nstio_swohdc(ioanl_nst,fname_nst,head_nst,data_nst,iret)

       write(6,101) lonb,latb,houra,iadate(1:4),iret
101    format(' WRITE_GFSNST:  nst analysis written  for ',&
            2i6,1x,f4.1,4(i4,1x),' with iret=',i2)

!      Deallocate local work arrays
       call sfcio_axdata(data_sfc,iret)
       call nstio_axdata(data_nst,iret)
       deallocate(buffer2)

    endif                               ! if (mype == mype_so ) then

!   End of routine
  end subroutine write_gfs_sfc_nst

  subroutine sfc_interpolate(a,na_lon,na_lat,b,ns_lon,ns_lat)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    sfc_interpolate --- interpolates from analysis grid to 
!                                    surface grid
!   prgrmmr:     derber -  initial version; org: np2
!
! abstract:      This routine interpolates a on analysis grid to b on 
!                surface grid
!
! program history log:
!   2008-02-26  derber  - original routine
!   2008-05-28  safford - add subprogram doc block, rm unused uses
!   2011-04-01  li - change kind of output field (b: single to r_kind)
!   2013-01-26  parrish - change from grdcrd to grdcrd1 (to allow successful debug compile on WCOSS)
!
!   input argument list:
!     na_lon  - number of longitude grid analysis 
!     na_lat  - number of latitude grid analysis
!     ns_lon  - number of longitude grid sfc 
!     ns_lat  - number of latitude grid sfc
!     a       - analysis values
!
!   output argument list:
!     b       - surface values
!
! attributes:
!   language: f90
!   machines: ibm RS/6000 SP; SGI Origin 2000; Compaq HP
!
!$$$ end documentation block

! !USES:
    use kinds, only: r_kind,i_kind,r_single
    use constants, only: zero,one
    use gridmod, only: rlats,rlons,rlats_sfc,rlons_sfc
    
    implicit none

! !INPUT PARAMETERS:
    integer(i_kind)                        ,intent(in   ) :: na_lon  ! number of longitude grid analysis 
    integer(i_kind)                        ,intent(in   ) :: na_lat  ! number of latitude grid analysis
    integer(i_kind)                        ,intent(in   ) :: ns_lon  ! number of longitude grid sfc 
    integer(i_kind)                        ,intent(in   ) :: ns_lat  ! number of latitude grid sfc

    real(r_kind), dimension(na_lon,na_lat),intent(in   ) :: a   ! analysis values

! !OUTPUT PARAMETERS:
    real(r_kind), dimension(ns_lon,ns_lat),intent(  out) :: b   ! surface values


!   Declare local variables
    integer(i_kind) i,j,ix,iy,ixp,iyp
    real(r_kind) dx1,dy1,dx,dy,w00,w01,w10,w11,bout,dlat,dlon

!*****************************************************************************

    b=zero
!   Loop over all points to get interpolated value
    do j=1,ns_lat
       dlat=rlats_sfc(j)
       call grdcrd1(dlat,rlats,na_lat,1)
       iy=int(dlat)
       iy=min(max(1,iy),na_lat)
       dy  =dlat-iy
       dy1 =one-dy
       iyp=min(na_lat,iy+1)


       do i=1,ns_lon
          dlon=rlons_sfc(i)
          call grdcrd1(dlon,rlons,na_lon,1)
          ix=int(dlon)
          dx  =dlon-ix
          dx=max(zero,min(dx,one))
          dx1 =one-dx
          w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy

          ix=min(max(0,ix),na_lon)
          ixp=ix+1
          if(ix==0) ix=na_lon
          if(ixp==na_lon+1) ixp=1
          bout=w00*a(ix,iy)+w01*a(ix,iyp)+w10*a(ixp,iy)+w11*a(ixp,iyp)
          b(i,j)=bout

       end do
    end do

    
!   End of routine
    return
  end subroutine sfc_interpolate


!-------------------------------------------------------------------------------
  subroutine sigio_cnvtdv8(im,ix,km,idvc,idvm,ntrac,iret,t,q,cpi,cnflg)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    sigio_cnvtdv8
!
!   prgrmmr:
!
! abstract:
!
! program history log:
!   2008-05-28  safford -- add subprogram doc block
!
!   input argument list:
!     im,ix,km,idvc,idvm,ntrac,cnflg
!     q, cpi
!     t
!
!   output argument list:
!     iret
!     t
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

    use kinds, only: i_kind,r_kind
    use constants, only: zero,one,fv
    implicit none
    integer(i_kind),intent(in   ) :: im,ix,km,idvc,idvm,ntrac,cnflg
    integer(i_kind),intent(  out) :: iret
    real(r_kind)   ,intent(in   ) :: q(ix,km,ntrac), cpi(0:ntrac)
    real(r_kind)   ,intent(inout) :: t(ix,km)
    integer(i_kind) :: thermodyn_id, n
    real(r_kind) :: xcp(ix,km), sumq(ix,km)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    iret=0
    thermodyn_id = mod(IDVM/10,10)
!
    if (thermodyn_id == 3 .and. idvc == 3) then
       xcp(1:im,:)  = zero
       sumq(1:im,:) = zero
       do n=1,NTRAC
          if( cpi(n) /= zero) then
             xcp(1:im,:)  = xcp(1:im,:)  + q(1:im,:,n) * cpi(n)
             sumq(1:im,:) = sumq(1:im,:) + q(1:im,:,n)
          endif
       enddo
       xcp(1:im,:)  = (one-sumq(1:im,:))*cpi(0) + xcp(1:im,:)   ! Mean Cp
!
    else
       xcp(1:im,:) = one + fv*Q(1:im,:,1)        ! Virt factor
    endif
    if (cnflg > 0) then
       t(1:im,:) = t(1:im,:) / xcp(1:im,:)
    else
       t(1:im,:) = t(1:im,:) * xcp(1:im,:)
    endif
!
    return
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine sigio_cnvtdv8

end module ncepgfs_io
