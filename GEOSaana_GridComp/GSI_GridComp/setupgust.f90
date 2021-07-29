module gust_setup
  implicit none
  private
  public:: setup
        interface setup; module procedure setupgust; end interface

contains
subroutine setupgust(obsLL,odiagLL,lunin,mype,bwork,awork,nele,nobs,is,conv_diagsave)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    setupgust    compute rhs for conventional surface gust
!   prgmmr: derber           org: np23                date: 2004-07-20
!
! abstract: For sea surface temperature observations
!              a) reads obs assigned to given mpi task (geographic region),
!              b) simulates obs from guess,
!              c) apply some quality control to obs,
!              d) load weight and innovation arrays used in minimization
!              e) collects statistics for runtime diagnostic output
!              f) writes additional diagnostic information to output file
!
! program history log:
!   2009-03-10  zhu
!   2011-02-18  zhu - update
!   2013-01-26  parrish - change from grdcrd to grdcrd1, 
!                          tintrp2a to tintrp2a1, tintrp2a11 (to allow successful debug compile on WCOSS)
!   2013-10-19  todling - metguess now holds background
!   2014-01-28  todling - write sensitivity slot indicator (ioff) to header of diagfile
!   2014-07-21  carley - ensure no division by 0 when calculating presw
!   2014-12-30  derber - Modify for possibility of not using obsdiag
!                          before retuning to setuprhsall.f90
!   2015-10-01  guo   - full res obvsr: index to allow redistribution of obsdiags
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!   2016-06-24  guo     - fixed the default value of obsdiags(:,:)%tail%luse to luse(i)
!                       . removed (%dlat,%dlon) debris.
!   2016-10-07  pondeca - if(.not.proceed) advance through input file first
!   2017-02-06  todling - add netcdf_diag capability; hidden as contained code
!   2017-02-09  guo     - Remove m_alloc, n_alloc.
!                       . Remove my_node with corrected typecast().
!   2018-01-08  pondeca - addd option l_closeobs to use closest obs to analysis
!                                     time in analysis
!   2020-02-26  todling - reset obsbin from hr to min
!
!   input argument list:
!     lunin    - unit from which to read observations
!     mype     - mpi task id
!     nele     - number of data elements per observation
!     nobs     - number of observations
!
!   output argument list:
!     bwork    - array containing information about obs-ges statistics
!     awork    - array containing information for data counts and gross checks
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use mpeu_util, only: die,perr
  use kinds, only: r_kind,r_single,r_double,i_kind

  use guess_grids, only: hrdifsig,nfldsig,ges_lnprsl, &
               geop_hgtl,sfcmod_gfs,sfcmod_mm5,comp_fact10     
  use m_obsdiagNode, only: obs_diag
  use m_obsdiagNode, only: obs_diags
  use m_obsdiagNode, only: obsdiagLList_nextNode
  use m_obsdiagNode, only: obsdiagNode_set
  use m_obsdiagNode, only: obsdiagNode_get
  use m_obsdiagNode, only: obsdiagNode_assert

  use obsmod, only: rmiss_single,&
                    lobsdiagsave,nobskeep,lobsdiag_allocated,time_offset,ianldate
  use m_obsNode, only: obsNode
  use m_gustNode, only: gustNode
  use m_gustNode, only: gustNode_appendto
  use m_obsLList, only: obsLList
  use obsmod, only: bmiss,luse_obsdiag
  use obsmod, only: netcdf_diag, binary_diag, dirname
  use nc_diag_write_mod, only: nc_diag_init, nc_diag_header, nc_diag_metadata, &
       nc_diag_write, nc_diag_data2d
  use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_get_dim, nc_diag_read_close
  use gsi_4dvar, only: nobs_bins,mn_obsbin,min_offset
  use oneobmod, only: magoberr,maginnov,oneobtest
  use gridmod, only: nsig
  use gridmod, only: get_ij,twodvar_regional
  use constants, only: zero,tiny_r_kind,one,one_tenth,half,wgtlim,rd,grav,&
            two,cg_term,three,four,huge_single,r1000,r3600,r100,&
            grav_ratio,flattening,grav,deg2rad,grav_equator,somigliana, &
            semi_major_axis,eccentricity
  use jfunc, only: jiter,last,miter
  use qcmod, only: dfact,dfact1,npres_print
  use convinfo, only: nconvtype,cermin,cermax,cgross,cvar_b,cvar_pg,ictype
  use convinfo, only: icsubtype
  use m_dtime, only: dtime_setup, dtime_check
  use gsi_bundlemod, only : gsi_bundlegetpointer
  use gsi_metguess_mod, only : gsi_metguess_get,gsi_metguess_bundle
  use rapidrefresh_cldsurf_mod, only: l_closeobs

  implicit none

  type(obsLList ),target,dimension(:),intent(inout):: obsLL
  type(obs_diags),target,dimension(:),intent(inout):: odiagLL

! Declare passed variables
  logical                                          ,intent(in   ) :: conv_diagsave
  integer(i_kind)                                  ,intent(in   ) :: lunin,mype,nele,nobs
  real(r_kind),dimension(100+7*nsig)               ,intent(inout) :: awork
  real(r_kind),dimension(npres_print,nconvtype,5,3),intent(inout) :: bwork
  integer(i_kind)                                  ,intent(in   ) :: is ! ndat index

! Declare external calls for code analysis
  external:: tintrp2a1,tintrp2a11
  external:: stop2

! Declare local parameters
  real(r_kind),parameter:: r0_1_bmiss=one_tenth*bmiss
  character(len=*),parameter:: myname='setupgust'

! Declare local variables
  
  real(r_double) rstation_id

  real(r_kind) gustges,dlat,dlon,ddiff,dtime,error,r0_001,thirty
  real(r_kind) scale,val2,rsig,rsigp,ratio,ressw2,ress,residual
  real(r_kind) obserrlm,obserror,val,valqc,rlow,rhgh,drpx
  real(r_kind) term,rwgt
  real(r_kind) cg_gust,wgross,wnotgross,wgt,arg,exp_arg,rat_err2
  real(r_kind) presw,factw,dpres,sfcchk
  real(r_kind) ratio_errors,tfact,fact,wflate,ten,psges,goverrd,zsges
  real(r_kind) slat,sin2,termg,termr,termrg,pobl
  real(r_kind) dz,zob,z1,z2,p1,p2,dz21,dlnp21,dstn
  real(r_kind) errinv_input,errinv_adjst,errinv_final
  real(r_kind) err_input,err_adjst,err_final,skint,sfcr
  real(r_kind),dimension(nobs):: dup
  real(r_kind),dimension(nsig)::prsltmp,zges
  real(r_kind),dimension(nele,nobs):: data
  real(r_single),allocatable,dimension(:,:)::rdiagbuf


  integer(i_kind) ier,ilon,ilat,ihgt,igust,ipres,id,itime,ikx,imaxerr,iqc
  integer(i_kind) iuse,ilate,ilone,istnelv,iprvd,isprvd
  integer(i_kind) i,nchar,nreal,k,k1,k2,ii,ikxx,nn,isli,ibin,ioff,ioff0,jj
  integer(i_kind) l,mm1
  integer(i_kind) idomsfc,iskint,iff10,isfcr
  
  logical,dimension(nobs):: luse,muse
  integer(i_kind),dimension(nobs):: ioid ! initial (pre-distribution) obs ID
  logical proceed

  character(8) station_id
  character(8),allocatable,dimension(:):: cdiagbuf
  character(8),allocatable,dimension(:):: cprvstg,csprvstg
  character(8) c_prvstg,c_sprvstg
  real(r_double) r_prvstg,r_sprvstg

  logical:: in_curbin, in_anybin
  type(gustNode),pointer:: my_head
  type(obs_diag),pointer:: my_diag
  type(obs_diags),pointer:: my_diagLL
  real(r_kind) :: hr_offset

  equivalence(rstation_id,station_id)
  equivalence(r_prvstg,c_prvstg)
  equivalence(r_sprvstg,c_sprvstg)
  
  real(r_kind),allocatable,dimension(:,:,:) :: ges_ps
  real(r_kind),allocatable,dimension(:,:,:) :: ges_z
  real(r_kind),allocatable,dimension(:,:,:) :: ges_gust
  type(obsLList),pointer,dimension(:):: gusthead
  gusthead => obsLL(:)

! Check to see if required guess fields are available
  call check_vars_(proceed)
  if(.not.proceed) then
     read(lunin)data,luse   !advance through input file
     return  ! not all vars available, simply return
  endif

! If require guess vars available, extract from bundle ...
  call init_vars_

!*********************************************************************************
! Read and reformat observations in work arrays.
  read(lunin)data,luse,ioid

!  index information for data array (see reading routine)
  ier=1       ! index of obs error
  ilon=2      ! index of grid relative obs location (x)
  ilat=3      ! index of grid relative obs location (y)
  ipres=4     ! index of pressure
  ihgt=5      ! index of observation elevation
  igust=6     ! index of gust observation
  id=7        ! index of station id
  itime=8     ! index of observation time in data array
  ikxx=9      ! index of ob type
  imaxerr=10  ! index of gust max error
  iqc=11      ! index of qulaity mark
  iuse=12     ! index of use parameter
  idomsfc=13  ! index of dominant surface type
  iskint=14   ! index of surface skin temperature
  iff10=15    ! index of 10 meter wind factor
  isfcr=16    ! index of surface roughness
  ilone=17    ! index of longitude (degrees)
  ilate=18    ! index of latitude (degrees)
  istnelv=19  ! index of station elevation (m)
  iprvd=20    ! index of provider
  isprvd=21   ! index of subprovider

  mm1=mype+1
  scale=one
  rsig=nsig
  thirty = 30.0_r_kind
  ten = 10.0_r_kind
  r0_001=0.001_r_kind
  rsigp=rsig+one
  goverrd=grav/rd

  do i=1,nobs
     muse(i)=nint(data(iuse,i)) <= jiter
  end do

! Check for missing data
  if (.not. oneobtest) then
  do i=1,nobs
    if (data(igust,i) > r0_1_bmiss)  then
       muse(i)=.false.
       data(igust,i)=rmiss_single   ! for diag output
    end if
  end do
  end if

! Check for duplicate observations at same location
  hr_offset=min_offset/60.0_r_kind
  dup=one
  do k=1,nobs
     do l=k+1,nobs
        if(data(ilat,k) == data(ilat,l) .and.  &
           data(ilon,k) == data(ilon,l) .and.  &
           data(ier,k) < r1000 .and. data(ier,l) < r1000 .and. &
           muse(k) .and. muse(l))then

           if(l_closeobs) then
              if(abs(data(itime,k)-hr_offset)<abs(data(itime,l)-hr_offset)) then
                  muse(l)=.false.
              else
                  muse(k)=.false.
              endif
           else
              tfact=min(one,abs(data(itime,k)-data(itime,l))/dfact1)
              dup(k)=dup(k)+one-tfact*tfact*(one-dfact)
              dup(l)=dup(l)+one-tfact*tfact*(one-dfact)
           endif
        end if
     end do
  end do



! If requested, save select data for output to diagnostic file
  if(conv_diagsave)then
     ii=0
     nchar=1
     ioff0=22
     nreal=ioff0
     if (lobsdiagsave) nreal=nreal+4*miter+1
     allocate(cdiagbuf(nobs),rdiagbuf(nreal,nobs))
     allocate(cprvstg(nobs),csprvstg(nobs))
     if(netcdf_diag) call init_netcdf_diag_
  end if

  call dtime_setup()
  do i=1,nobs
     dtime=data(itime,i)
     call dtime_check(dtime, in_curbin, in_anybin)
     if(.not.in_anybin) cycle

     if(in_curbin) then
        dlat=data(ilat,i)
        dlon=data(ilon,i)

        ikx  = nint(data(ikxx,i))
        error=data(ier,i)
        isli=data(idomsfc,i)
     endif

!    Link observation to appropriate observation bin
     if (nobs_bins>1) then
        ibin = NINT( dtime*60/mn_obsbin ) + 1
     else
        ibin = 1
     endif
     IF (ibin<1.OR.ibin>nobs_bins) write(6,*)mype,'Error nobs_bins,ibin= ',nobs_bins,ibin

     if (luse_obsdiag) my_diagLL => odiagLL(ibin)

!    Link obs to diagnostics structure
     if (luse_obsdiag) then
        my_diag => obsdiagLList_nextNode(my_diagLL      ,&
                create = .not.lobsdiag_allocated        ,&
                   idv = is             ,&
                   iob = ioid(i)        ,&
                   ich = 1              ,&
                  elat = data(ilate,i)  ,&
                  elon = data(ilone,i)  ,&
                  luse = luse(i)        ,&
                 miter = miter          )

        if(.not.associated(my_diag)) call die(myname, &
                'obsdiagLList_nextNode(), create =', .not.lobsdiag_allocated)
     endif

     if(.not.in_curbin) cycle

! Interpolate to get gust at obs location/time
     call tintrp2a11(ges_gust,gustges,dlat,dlon,dtime,hrdifsig,&
          mype,nfldsig)

!   Process observations with reported height
    drpx = zero
    dpres = data(ihgt,i)
    dstn = data(istnelv,i)

!   Get guess surface elevation and geopotential height profile
!   at observation location.
    call tintrp2a11(ges_z,zsges,dlat,dlon,dtime,hrdifsig,&
            mype,nfldsig)
!   Subtract off combination of surface station elevation and
!   model elevation depending on how close to surface
    fact = zero
    if(dpres-dstn > 10._r_kind)then
      if(dpres-dstn > 1000._r_kind)then
         fact = one
      else
         fact=(dpres-dstn)/990._r_kind
      end if
    end if
    dpres=dpres-(dstn+fact*(zsges-dstn))
    drpx=0.003*abs(dstn-zsges)*(one-fact)

    if (.not. twodvar_regional) then
       call tintrp2a1(geop_hgtl,zges,dlat,dlon,dtime,hrdifsig,&
               nsig,mype,nfldsig)
!      For observation reported with geometric height above sea level,
!      convert geopotential to geometric height.
!      Convert geopotential height at layer midpoints to geometric
!      height using equations (17, 20, 23) in MJ Mahoney's note
!      "A discussion of various measures of altitude" (2001).
!      Available on the web at
!      http://mtp.jpl.nasa.gov/notes/altitude/altitude.html
!
!      termg  = equation 17
!      termr  = equation 21
!      termrg = first term in the denominator of equation 23
!      zges  = equation 23

       slat = data(ilate,i)*deg2rad
       sin2  = sin(slat)*sin(slat)
       termg = grav_equator * &
            ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
       termr = semi_major_axis /(one + flattening + grav_ratio -  &
            two*flattening*sin2)
       termrg = (termg/grav)*termr
       do k=1,nsig
          zges(k) = (termr*zges(k)) / (termrg-zges(k))  ! eq (23)
       end do
    else
       zges(1) = ten
    end if

!   Given observation height, (1) adjust 10 meter wind factor if
!   necessary, (2) convert height to grid relative units, (3) compute
!   compute observation pressure (for diagnostic purposes only), and
!   (4) compute location of midpoint of first model layer above surface
!   in grid relative units

!   Convert observation height (in dpres) from meters to grid relative
!   units.  Save the observation height in zob for later use.
    zob = dpres
    call grdcrd1(dpres,zges,nsig,1)

    if (zob >= zges(1)) then
       factw=one
    else
       factw = data(iff10,i)
       if(sfcmod_gfs .or. sfcmod_mm5) then
          sfcr = data(isfcr,i)
          skint = data(iskint,i)
          call comp_fact10(dlat,dlon,dtime,skint,sfcr,isli,mype,factw)
       end if
       if (.not. twodvar_regional) then
         if (zob <= ten) then
            if(zob < ten)then
              term = max(zob,zero)/ten
              factw = term*factw
            end if
         else
            term = (zges(1)-zob)/(zges(1)-ten)
            factw = one-term+factw*term
         end if
       else
          if(zob < ten)then
             term = max(zob,zero)/ten
             factw = term*factw
          end if
       end if
       gustges=factw*gustges
    endif

!   Compute observation pressure (only used for diagnostics & for type 2**)
!   Get guess surface pressure and mid layer pressure
!   at observation location.
    if (ictype(ikx)>=280 .and. ictype(ikx)<290) then
       call tintrp2a11(ges_ps,psges,dlat,dlon,dtime,hrdifsig,&
            mype,nfldsig)
       call tintrp2a1(ges_lnprsl,prsltmp,dlat,dlon,dtime,hrdifsig,&
            nsig,mype,nfldsig)
       if ((dpres-one) < tiny_r_kind) then
          z1=zero;    p1=log(psges)
          z2=zges(1); p2=prsltmp(1)
       elseif (dpres>nsig) then
          z1=zges(nsig-1); p1=prsltmp(nsig-1)
          z2=zges(nsig);   p2=prsltmp(nsig)
          drpx = 1.e6_r_kind
       else
          k=dpres
          k1=min(max(1,k),nsig)
          k2=max(1,min(k+1,nsig))
          z1=zges(k1); p1=prsltmp(k1)
          z2=zges(k2); p2=prsltmp(k2)
       endif

       dz21     = z2-z1
       if(dz21==zero)cycle
       dlnp21   = p2-p1
       dz       = zob-z1
       pobl     = p1 + (dlnp21/dz21)*dz
       presw    = ten*exp(pobl)
    else
       presw = ten*exp(data(ipres,i))
    end if


!   Determine location in terms of grid units for midpoint of
!   first layer above surface
    sfcchk=zero
    call grdcrd1(sfcchk,zges,nsig,1)

!    Checks based on observation location relative to model surface and top
     rlow=max(sfcchk-dpres,zero)
     rhgh=max(dpres-r0_001-rsigp,zero)
     if(luse(i))then
        awork(1) = awork(1) + one
        if(rlow/=zero) awork(2) = awork(2) + one
        if(rhgh/=zero) awork(3) = awork(3) + one
     end if
 
!    Adjust observation error
!    ratio_errors=error/((data(ier,i)+adjustment)*sqrt(dup(i)))  ! qc dependent adjustment
     wflate=zero
     if (ictype(ikx)==188 .or. ictype(ikx)==288 .or. ictype(ikx)==195 .or. ictype(ikx)==295) then  !inflate Mesonet obs error for gusts<7.2m/s
       if (data(igust,i)<7.2) then
          wflate=4.0_r_kind*data(ier,i)
       else
          wflate=0.8_r_kind*data(ier,i)
       end if
     end if
     ratio_errors=error/((data(ier,i)+drpx+wflate+1.0e6*rhgh+four*rlow)*sqrt(dup(i)))
     error=one/error

!    Compute innovations
     ddiff=data(igust,i)-gustges

!    If requested, setup for single obs test.
     if (oneobtest) then
        ddiff=maginnov
        error=one/magoberr
        ratio_errors=one
     endif

!    Gross check using innovation normalized by error
     obserror = one/max(ratio_errors*error,tiny_r_kind)
     obserrlm = max(cermin(ikx),min(cermax(ikx),obserror))
     residual = abs(ddiff)
     ratio    = residual/obserrlm
     if (ratio> cgross(ikx) .or. ratio_errors < tiny_r_kind) then
        if (luse(i)) awork(6) = awork(6)+one
        error = zero
        ratio_errors=zero
     end if
     if (ratio_errors*error <=tiny_r_kind) muse(i)=.false.

     if (nobskeep>0 .and. luse_obsdiag) call obsdiagNode_get(my_diag, jiter=nobskeep, muse=muse(i))

!    Compute penalty terms (linear & nonlinear qc).
     val      = error*ddiff
     if(luse(i))then
        val2     = val*val
        exp_arg  = -half*val2
        rat_err2 = ratio_errors**2
        if (cvar_pg(ikx) > tiny_r_kind .and. error > tiny_r_kind) then
           arg  = exp(exp_arg)
           wnotgross= one-cvar_pg(ikx)
           cg_gust=cvar_b(ikx)
           wgross = cg_term*cvar_pg(ikx)/(cg_gust*wnotgross)
           term = log((arg+wgross)/(one+wgross))
           wgt  = one-wgross/(arg+wgross)
           rwgt = wgt/wgtlim
        else
           term = exp_arg
           wgt  = wgtlim
           rwgt = wgt/wgtlim
        endif
        valqc = -two*rat_err2*term

!       Accumulate statistics for obs belonging to this task
        if (muse(i)) then
           if(rwgt < one) awork(21) = awork(21)+one
           awork(4)=awork(4)+val2*rat_err2
           awork(5)=awork(5)+one
           awork(22)=awork(22)+valqc
        end if
        ress   = ddiff*scale
        ressw2 = ress*ress
        val2   = val*val
        rat_err2 = ratio_errors**2
        nn=1
        if (.not. muse(i)) then
           nn=2
           if(ratio_errors*error >=tiny_r_kind)nn=3
        end if
        if (abs(data(igust,i)-rmiss_single) >=tiny_r_kind) then
           bwork(1,ikx,1,nn)  = bwork(1,ikx,1,nn)+one           ! count
           bwork(1,ikx,2,nn)  = bwork(1,ikx,2,nn)+ress          ! (o-g)
           bwork(1,ikx,3,nn)  = bwork(1,ikx,3,nn)+ressw2        ! (o-g)**2
           bwork(1,ikx,4,nn)  = bwork(1,ikx,4,nn)+val2*rat_err2 ! penalty
           bwork(1,ikx,5,nn)  = bwork(1,ikx,5,nn)+valqc         ! nonlin qc penalty
        end if

     endif

     if (luse_obsdiag) then
        call obsdiagNode_set(my_diag, wgtjo=(error*ratio_errors)**2, &
          jiter=jiter, muse=muse(i), nldepart=ddiff)
     endif

!    If obs is "acceptable", load array with obs info for use
!    in inner loop minimization (int* and stp* routines)
     if (.not. last .and. muse(i)) then

        allocate(my_head)
        call gustNode_appendto(my_head,gusthead(ibin))

        my_head%idv = is
        my_head%iob = ioid(i)
        my_head%elat= data(ilate,i)
        my_head%elon= data(ilone,i)

!       Set (i,j) indices of guess gridpoint that bound obs location
        call get_ij(mm1,dlat,dlon,my_head%ij,my_head%wij)

        my_head%res     = ddiff
        my_head%err2    = error**2
        my_head%raterr2 = ratio_errors**2    
        my_head%time    = dtime
        my_head%b       = cvar_b(ikx)
        my_head%pg      = cvar_pg(ikx)
        my_head%luse    = luse(i)
        if(luse_obsdiag) then
           call obsdiagNode_assert(my_diag, my_head%idv,my_head%iob,1,myname,'my_diag:my_head')
           my_head%diags => my_diag
        endif
        my_head => null()
     endif


!    Save stuff for diagnostic output
     if(conv_diagsave .and. luse(i))then
        ii=ii+1
        rstation_id = data(id,i)
        err_input   = data(ier,i)
        err_adjst   = data(ier,i)
        if (ratio_errors*error>tiny_r_kind) then
           err_final = one/(ratio_errors*error)
        else
           err_final = huge_single
        endif
 
        errinv_input = huge_single
        errinv_adjst = huge_single
        errinv_final = huge_single
        if (err_input>tiny_r_kind) errinv_input = one/err_input
        if (err_adjst>tiny_r_kind) errinv_adjst = one/err_adjst
        if (err_final>tiny_r_kind) errinv_final = one/err_final

        if(binary_diag) call contents_binary_diag_(my_diag)
        if(netcdf_diag) call contents_netcdf_diag_(my_diag)
 
     end if


  end do

! Release memory of local guess arrays
  call final_vars_

! Write information to diagnostic file
  if(conv_diagsave) then
     if(netcdf_diag) call nc_diag_write
     if(binary_diag .and. ii>0)then
        write(7)'gst',nchar,nreal,ii,mype,ioff0
        write(7)cdiagbuf(1:ii),rdiagbuf(:,1:ii)
        deallocate(cdiagbuf,rdiagbuf)

        write(7)cprvstg(1:ii),csprvstg(1:ii)
        deallocate(cprvstg,csprvstg)
     end if
  end if

! End of routine

  return
  contains

  subroutine check_vars_ (proceed)
  logical,intent(inout) :: proceed
  integer(i_kind) ivar, istatus
! Check to see if required guess fields are available
  call gsi_metguess_get ('var::ps', ivar, istatus )
  proceed=ivar>0
  call gsi_metguess_get ('var::z' , ivar, istatus )
  proceed=proceed.and.ivar>0
  call gsi_metguess_get ('var::gust' , ivar, istatus )
  proceed=proceed.and.ivar>0
  end subroutine check_vars_ 

  subroutine init_vars_

  real(r_kind),dimension(:,:  ),pointer:: rank2=>NULL()
  character(len=5) :: varname
  integer(i_kind) ifld, istatus

! If require guess vars available, extract from bundle ...
  if(size(gsi_metguess_bundle)==nfldsig) then
!    get gust ...
     varname='gust'
     call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank2,istatus)
     if (istatus==0) then
         if(allocated(ges_gust))then
            write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
            call stop2(999)
         endif
         allocate(ges_gust(size(rank2,1),size(rank2,2),nfldsig))
         ges_gust(:,:,1)=rank2
         do ifld=2,nfldsig
            call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank2,istatus)
            ges_gust(:,:,ifld)=rank2
         enddo
     else
         write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
         call stop2(999)
     endif
!    get ps ...
     varname='ps'
     call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank2,istatus)
     if (istatus==0) then
         if(allocated(ges_ps))then
            write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
            call stop2(999)
         endif
         allocate(ges_ps(size(rank2,1),size(rank2,2),nfldsig))
         ges_ps(:,:,1)=rank2
         do ifld=2,nfldsig
            call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank2,istatus)
            ges_ps(:,:,ifld)=rank2
         enddo
     else
         write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
         call stop2(999)
     endif
!    get z ...
     varname='z'
     call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank2,istatus)
     if (istatus==0) then
         if(allocated(ges_z))then
            write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
            call stop2(999)
         endif
         allocate(ges_z(size(rank2,1),size(rank2,2),nfldsig))
         ges_z(:,:,1)=rank2
         do ifld=2,nfldsig
            call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank2,istatus)
            ges_z(:,:,ifld)=rank2
         enddo
     else
         write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
         call stop2(999)
     endif
  else
     write(6,*) trim(myname), ': inconsistent vector sizes (nfldsig,size(metguess_bundle) ',&
                 nfldsig,size(gsi_metguess_bundle)
     call stop2(999)
  endif
  end subroutine init_vars_

  subroutine init_netcdf_diag_
  character(len=80) string
  character(len=128) diag_conv_file
  integer(i_kind) ncd_fileid,ncd_nobs
  logical append_diag
  logical,parameter::verbose=.false. 
     write(string,900) jiter
900  format('conv_gust_',i2.2,'.nc4')
     diag_conv_file=trim(dirname) // trim(string)

     inquire(file=diag_conv_file, exist=append_diag)

     if (append_diag) then
        call nc_diag_read_init(diag_conv_file,ncd_fileid)
        ncd_nobs = nc_diag_read_get_dim(ncd_fileid,'nobs')
        call nc_diag_read_close(diag_conv_file)

        if (ncd_nobs > 0) then
           if(verbose) print *,'file ' // trim(diag_conv_file) // ' exists.  Appending.  nobs,mype=',ncd_nobs,mype
        else
           if(verbose) print *,'file ' // trim(diag_conv_file) // ' exists but contains no obs.  Not appending. nobs,mype=',ncd_nobs,mype
           append_diag = .false. ! if there are no obs in existing file, then do not try to append
        endif
     end if

     call nc_diag_init(diag_conv_file, append=append_diag)

     if (.not. append_diag) then ! don't write headers on append - the module will break?
        call nc_diag_header("date_time",ianldate )
     endif
  end subroutine init_netcdf_diag_
  subroutine contents_binary_diag_(odiag)
  type(obs_diag),pointer,intent(in):: odiag
        cdiagbuf(ii)    = station_id         ! station id
 
        rdiagbuf(1,ii)  = ictype(ikx)        ! observation type
        rdiagbuf(2,ii)  = icsubtype(ikx)     ! observation subtype
 
        rdiagbuf(3,ii)  = data(ilate,i)      ! observation latitude (degrees)
        rdiagbuf(4,ii)  = data(ilone,i)      ! observation longitude (degrees)
        rdiagbuf(5,ii)  = data(istnelv,i)    ! station elevation (meters)
        rdiagbuf(6,ii)  = presw              ! observation pressure (hPa)
        rdiagbuf(7,ii)  = data(ihgt,i)       ! observation height (meters)
        rdiagbuf(8,ii)  = dtime-time_offset  ! obs time (hours relative to analysis time)

        rdiagbuf(9,ii)  = data(iqc,i)        ! input prepbufr qc or event mark
        rdiagbuf(10,ii) = rmiss_single       ! setup qc or event mark
        rdiagbuf(11,ii) = data(iuse,i)       ! read_prepbufr data usage flag
        if(muse(i)) then
           rdiagbuf(12,ii) = one             ! analysis usage flag (1=use, -1=not used)
        else
           rdiagbuf(12,ii) = -one
        endif

        rdiagbuf(13,ii) = rwgt               ! nonlinear qc relative weight
        rdiagbuf(14,ii) = errinv_input       ! prepbufr inverse obs error (K**-1)
        rdiagbuf(15,ii) = errinv_adjst       ! read_prepbufr inverse obs error (K**-1)
        rdiagbuf(16,ii) = errinv_final       ! final inverse observation error (K**-1)
 
        rdiagbuf(17,ii) = data(igust,i)      ! GUST observation (K)
        rdiagbuf(18,ii) = ddiff              ! obs-ges used in analysis (K)
        rdiagbuf(19,ii) = data(igust,i)-gustges! obs-ges w/o bias correction (K) (future slot)
 
        rdiagbuf(20,ii) = factw              ! 10m wind reduction factor

        rdiagbuf(21,ii) = data(idomsfc,i)    ! dominate surface type
        rdiagbuf(22,ii) = zsges              ! model terrain at ob location
        r_prvstg        = data(iprvd,i)
        cprvstg(ii)     = c_prvstg           ! provider name
        r_sprvstg       = data(isprvd,i)
        csprvstg(ii)    = c_sprvstg          ! subprovider name

        if (lobsdiagsave) then
           ioff=ioff0
           do jj=1,miter 
              ioff=ioff+1 
              if (odiag%muse(jj)) then
                 rdiagbuf(ioff,ii) = one
              else
                 rdiagbuf(ioff,ii) = -one
              endif
           enddo
           do jj=1,miter+1
              ioff=ioff+1
              rdiagbuf(ioff,ii) = odiag%nldepart(jj)
           enddo
           do jj=1,miter
              ioff=ioff+1
              rdiagbuf(ioff,ii) = odiag%tldepart(jj)
           enddo
           do jj=1,miter
              ioff=ioff+1
              rdiagbuf(ioff,ii) = odiag%obssen(jj)
           enddo
        endif
  end subroutine contents_binary_diag_
  subroutine contents_netcdf_diag_(odiag)
  type(obs_diag),pointer,intent(in):: odiag
! Observation class
  character(7),parameter     :: obsclass = '   gust'
  real(r_kind),dimension(miter) :: obsdiag_iuse
           call nc_diag_metadata("Station_ID",              station_id             )
           call nc_diag_metadata("Observation_Class",       obsclass               )
           call nc_diag_metadata("Observation_Type",        ictype(ikx)            )
           call nc_diag_metadata("Observation_Subtype",     icsubtype(ikx)         )
           call nc_diag_metadata("Latitude",                data(ilate,i)          )
           call nc_diag_metadata("Longitude",               data(ilone,i)          )
           call nc_diag_metadata("Station_Elevation",       data(istnelv,i)        )
           call nc_diag_metadata("Pressure",                presw*r100             )
           call nc_diag_metadata("Height",                  data(ihgt,i)           )
           call nc_diag_metadata("Time",                    dtime-time_offset      )
           call nc_diag_metadata("Prep_QC_Mark",            data(iqc,i)            )
           call nc_diag_metadata("Prep_Use_Flag",           data(iuse,i)           )
!          call nc_diag_metadata("Nonlinear_QC_Var_Jb",     var_jb                 )
           call nc_diag_metadata("Nonlinear_QC_Rel_Wgt",    rwgt                   )                 
           if(muse(i)) then
              call nc_diag_metadata("Analysis_Use_Flag",    one                    )
           else
              call nc_diag_metadata("Analysis_Use_Flag",   -one                    )              
           endif

           call nc_diag_metadata("Errinv_Input",            errinv_input           )
           call nc_diag_metadata("Errinv_Adjust",           errinv_adjst           )
           call nc_diag_metadata("Errinv_Final",            errinv_final           )

           call nc_diag_metadata("Observation",                   data(igust,i)    )
           call nc_diag_metadata("Obs_Minus_Forecast_adjusted",   ddiff            )
           call nc_diag_metadata("Obs_Minus_Forecast_unadjusted", data(igust,i)-gustges )
 
           call nc_diag_metadata("Dominant_Sfc_Type", data(idomsfc,i)              )
           call nc_diag_metadata("Model_Terrain",     zsges                        )
           r_prvstg            = data(iprvd,i)
           call nc_diag_metadata("Provider_Name",     c_prvstg                     )    
           r_sprvstg           = data(isprvd,i)
           call nc_diag_metadata("Subprovider_Name",  c_sprvstg                    )

           if (lobsdiagsave) then
              do jj=1,miter
                 if (odiag%muse(jj)) then
                       obsdiag_iuse(jj) =  one
                 else
                       obsdiag_iuse(jj) = -one
                 endif
              enddo
   
              call nc_diag_data2d("ObsDiagSave_iuse",     obsdiag_iuse                             )
              call nc_diag_data2d("ObsDiagSave_nldepart", odiag%nldepart )
              call nc_diag_data2d("ObsDiagSave_tldepart", odiag%tldepart )
              call nc_diag_data2d("ObsDiagSave_obssen",   odiag%obssen   )             
           endif
   
  end subroutine contents_netcdf_diag_

  subroutine final_vars_
    if(allocated(ges_z   )) deallocate(ges_z   )
    if(allocated(ges_ps  )) deallocate(ges_ps  )
    if(allocated(ges_gust)) deallocate(ges_gust)
  end subroutine final_vars_

end subroutine setupgust
end module gust_setup
