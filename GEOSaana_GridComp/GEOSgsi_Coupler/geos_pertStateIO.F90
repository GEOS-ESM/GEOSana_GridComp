!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  geos_pertStateIO --- IO of GSI perturbation states
!
! !INTERFACE:
!
  module geos_pertStateIO


  use ESMF, only: ESMF_MAXGRIDDIM
  use ESMF, only: ESMF_MAXSTR
  use ESMF, only: ESMF_Field
  use ESMF, only: ESMF_FieldBundle
  use ESMF, only: ESMF_FieldBundleGet
  use ESMF, only: ESMF_FieldGet
  use MAPL

  use GSI_GridCompMod, only: GSI_bkg_fname_tmpl
  use GSI_GridCompMod, only: GSI_fsens_fname_tmpl
  use GSI_GridCompMod, only: GSI_ferrA_fname_tmpl
  use GSI_GridCompMod, only: GSI_ferrB_fname_tmpl

  use kinds,       only : r_kind,r_single,i_kind
  use mpimod,      only : mype,mpi_rtype,mpi_comm_world
  use mpimod,      only : nxpe,nype
  use gridmod,     only : strip
  use gridmod,     only : displs_s,ijn_s
  use gridmod,     only : nlat, nlon     ! no. lat/lon
  use gridmod,     only : nsig           ! no. levels
  use gridmod,     only : iglobal        ! no. of horizontal points on global grid
  use gridmod,     only : ijn            ! no. of horiz. pnts for each subdomain (no buffer)
  use gridmod,     only : displs_g       ! comm. array, displacement for receive on global grid
  use gridmod,     only : itotsub        ! no. of horizontal points of all subdomains combined
  use gridmod,     only : bk5         
  use gsi_bias,    only : reorder21,reorder12
  use gsi_4dvar,   only : ibdate,l4densvar
  use obsmod,      only : iadate
  use constants,   only : zero,one,tiny_r_kind
  use state_vectors,only: dot_product

  use m_interpack,    only : interpack_terpv
  use m_interpack_ad, only : interpack_terpv_ad

  use m_StrTemplate, only: StrTemplate

  use mpeu_util, only: tell,warn,perr,die
  use GSI_GridCompMod, only: PPMV2GpG

  implicit none
  private
  public:: pertState_get_set
  public:: pertState_get
  public:: pertState_get_unset
  public:: pertState_put
  public:: pertState_put_set
  public:: pertState_put_final

  interface pertState_get; module procedure &
    get_1pert_, &
    get_1pert_date_, &
    get_Npert_; end interface
  interface pertState_put; module procedure &
    put_1pert_, &
    put_Npert_; end interface
  interface pertState_put_set; module procedure &
    set_gsi2pgcm0_; end interface
  interface pertState_put_final; module procedure &
    final_gsi2pgcm0_; end interface
  interface pertState_get_set; module procedure &
    get_pert_set_; end interface
  interface pertState_get_unset; module procedure &
    get_pert_unset_; end interface

  character(len=*),parameter :: myname="geos_pertStateIO"
#include "MAPL_ErrLog.h"


   character(len=*),parameter:: fnxgsi = 'xxgsi.eta'

   integer,save :: mycount = 0

   logical :: gsiMode = .false.

   integer(i_kind),parameter:: ROOT =0
   real(r_kind), parameter :: kPa_per_Pa = 0.001_r_kind
   real(r_kind), parameter :: Pa_per_kPa = 1000._r_kind
   real(r_kind), parameter :: R3600      = 3600.0_r_kind

   type(ESMF_FieldBundle)  :: WBundle

   integer(i_kind) ::  nlon_set = -1
   integer(i_kind) ::  nlat_set = -1
   integer(i_kind) ::  lon2_set = -1
   integer(i_kind) ::  lat2_set = -1
   integer(i_kind) ::  nymd_set = -1
   integer(i_kind) ::  nhms_set = -1
contains
!------------------------------------------------------------------------------------
subroutine get_1pert_(xx,what,filename)
! get perturbation from user''s model and convert it to relevant gsi bundle
use constants, only: zero
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: assignment(=)
implicit none
type(gsi_bundle),intent(inout) :: xx
character(len=*),intent(in) :: what     ! indicates whether tl or ad type perturbation
character(len=*),intent(in) :: filename ! filename w/ perturbation

character(len=*),parameter:: myname_=myname//'.get_1pert_'
integer(i_kind):: ierr

  xx=zero
  call pgcm2gsi0_(xx,what,ierr,filename=filename)
 	if(ierr/=0) then
 	  call perr(myname_,'pgcm2gsi("'//trim(what)//'"), ierr =',ierr)
 	  call perr(myname_,'trouble reading analysis errors')
 	  call die(myname_)
 	endif
end subroutine get_1pert_
!------------------------------------------------------------------------------------
subroutine get_1pert_date_(xx,nt,what,filename)
! get perturbation from user''s model and convert it to relevant gsi bundle
use kinds, only: i_kind,r_kind
use gsi_4dvar, only: ens_fmnlevs
use constants, only: zero
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: assignment(=)
implicit none
type(gsi_bundle),intent(inout) :: xx
integer(i_kind),intent(in)  :: nt
character(len=*),intent(in) :: what     ! indicates whether tl or ad type perturbation
character(len=*),intent(in) :: filename ! filename w/ perturbation

character(len=*),parameter:: myname_=myname//'.get_1pert_'
integer(i_kind):: nymd,nhms
integer(i_kind):: ierr
integer(i_kind) ida(8),jda(8)
real(r_kind) fha(5)

   ida(1:3)=ibdate(1:3)
   ida(5:6)=ibdate(4:5)
   jda(:)=0
   fha(:)=0.0
   if (l4densvar) then
      fha(3)=ens_fmnlevs(nt)-180.0_r_kind ! NCEP counts time from previous syn analysis (180min=3hr)
   else
      nymd = iadate(1)*10000 + iadate(2)*100 + iadate(3)
      nhms = iadate(4)*10000 + iadate(5)*100
   endif
   call w3movdat(fha,ida,jda)
   nymd=jda(1)*10000+jda(2)*100+jda(3)
   nhms=jda(5)*10000+jda(6)*100
   if(mype==0) print *, myname_, ': will read perturbations at: ', nymd, ' ', nhms

  xx=zero
  call pgcm2gsi0_(xx,what,ierr,filename=filename)
 	if(ierr/=0) then
 	  call perr(myname_,'pgcm2gsi("'//trim(what)//'"), ierr =',ierr)
 	  call perr(myname_,'trouble reading analysis errors')
 	  call die(myname_)
 	endif
end subroutine get_1pert_date_
!------------------------------------------------------------------------------------
subroutine put_1pert_(xx,nymd,nhms,what,label)
! convert xx to the user''s model perturbation and write it out
use kinds, only: i_kind
use gsi_bundlemod, only: gsi_bundle
implicit none
type(gsi_bundle),intent(inout) :: xx     ! gsi perturbation (bundle) vector
character(len=*),intent(in)    :: what   ! indicates whether tl or ad type perturbation
character(len=*),intent(in)    :: label  ! label used to identify output filename
integer(i_kind), intent(in)    :: nymd   ! date to write out field, as in, YYYYMMDD
integer(i_kind), intent(in)    :: nhms   ! time to write out field, as in, HHMMSS

character(len=*),parameter:: myname_=myname//'.put_1pert_'
integer(i_kind):: ierr

! Write out perturbation consistent with underlying atmospheric model

  call gsi2pgcm0_(nymd,nhms,xx,what,ierr, &
  	filename=trim(label)//'.eta')
	if(ierr/=0) then
	  call perr(myname_,'gsi2pgcm0_("'//trim(what)//'"), ierr =',ierr)
	  call perr(myname_,'filename =',trim(label)//'.eta')
	  call perr(myname_,'trouble writing analysis errors')
	  call die(myname_)
	endif
end subroutine put_1pert_

!------------------------------------------------------------------------------------
subroutine get_Npert_(xx,n,what,filename)
! get perturbation from user''s model and convert it to relevant gsi bundle
use kinds, only: i_kind
use constants, only: zero
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: assignment(=)
implicit none
integer(i_kind),intent(in) :: n
type(gsi_bundle),intent(inout) :: xx(n)
character(len=*),intent(in) :: what   ! indicates whether tl or ad type perturbation
character(len=*),intent(in) :: filename(n)

integer(i_kind):: ii,ierr
character(len=*),parameter:: myname_=myname//'.get_Npert_'

  if(n>1) call die(myname_,'not ready for nsubwin>1, n =',n)

  do ii=1,n
    xx(ii)=zero
    call pgcm2gsi0_(xx(ii),what,ierr,filename=filename(ii))
     	if(ierr/=0) then
	  call perr(myname_,'pgcm2gsi0_("'//trim(what)//'"), ierr =',ierr)
	  call perr(myname_,'trouble reading analysis errors, ii =',ii)
	  call die(myname_)
	endif
  enddo

end subroutine get_Npert_
!------------------------------------------------------------------------------------
subroutine put_Npert_(xx,n,what)
! convert xx to the user''s model perturbation and write it out
use kinds, only: i_kind
use gsi_bundlemod, only: gsi_bundle
implicit none
integer(i_kind),intent(in) :: n
type(gsi_bundle),intent(in) :: xx(n)     ! gsi perturbation (bundle) vector
character(len=*),intent(in) :: what      ! indicates whether tl or ad type perturbation

integer(i_kind):: ii,ierr
character(len=*),parameter:: myname_=myname//'.put_Npert_'

  call perr(myname_,'not implemented')
  do ii=1,n
    ierr=1
    !call gsi2pgcm0_(xx(ii),what,ierr)
	if(ierr/=0) then
	  call perr(myname_,'gsi2pgcm0_("'//trim(what)//'"), ierr =',ierr)
	  call perr(myname_,'trouble writing analysis errors, ii =',ii)
	  call die(myname_)
	endif
  enddo
end subroutine put_Npert_
!------------------------------------------------------------------------------------
subroutine get_pert_set_ (nlat_in,nlon_in,lat2_in,lon2_in,nt)
   use gsi_4dvar, only: ens_fmnlevs,l4densvar
   implicit none
   integer(i_kind), intent(in) :: lat2_in,lon2_in
   integer(i_kind), intent(in) :: nlat_in,nlon_in
   integer(i_kind), intent(in) :: nt

   character(len=*),parameter:: myname_=myname//'.get_set_'
   integer(i_kind) ida(8),jda(8)
   real(r_kind) fha(5)

   nlat_set=nlat_in
   nlon_set=nlon_in
   lat2_set=lat2_in
   lon2_set=lon2_in

   if (nt>0) then
      ida(1:3)=ibdate(1:3)
      ida(5:6)=ibdate(4:5)
      jda(:)=0
      fha(:)=0.0
      if(l4densvar) then
         fha(3)=ens_fmnlevs(nt)-180.0_r_kind ! NCEP counts time from previous syn analysis (180min=3hr)
      else
         fha(3)=ens_fmnlevs(nt)
      endif
      call w3movdat(fha,ida,jda)
      nymd_set=jda(1)*10000+jda(2)*100+jda(3)
      nhms_set=jda(5)*10000+jda(6)*100
      if(mype==0) print *, myname_, ': will read perturbations at: ', nymd_set, ' ', nhms_set
   endif

end subroutine get_pert_set_
!------------------------------------------------------------------------------------
subroutine get_pert_unset_
   implicit none
   nymd_set=-1
   nhms_set=-1
   nlat_set=-1
   nlon_set=-1
   lat2_set=-1
   lon2_set=-1
end subroutine get_pert_unset_
!------------------------------------------------------------------------------------

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_gsi2pgcm0_:  Convert GSI increments to GCM perturbations
!
! !INTERFACE:

      subroutine set_gsi2pgcm0_ ( nymd, nhms, stat )

! !USES:

      use ESMF

      use GSI_GridCompMod, only: GSI_AgcmPertGrid
      use GSI_GridCompMod, only: GSI_ExpId
      use GSI_GridCompMod, only: GSI_RefTime

      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: assignment(=)

      use m_StrTemplate, only: StrTemplate

      implicit none

! !INPUT PARAMETERS:
                                                                                                                           
      integer(i_kind),    intent(in)    :: nymd   ! date as in YYYYMMDD
      integer(i_kind),    intent(in)    :: nhms   ! time as in HHMMSS

! !OUTPUT PARAMETERS:

      integer(i_kind),                 intent(out) :: stat

! !DESCRIPTION: Set work bundle for writing output.
!
! !REVISION HISTORY:
!
!  03Mar2020  Todling   Broken up from original gsi2pgcm0_
!
!EOP
!-----------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'*set_gsi2pgcm0_'
     character(len=*), parameter :: Iam = myname_
     type(MAPL_CFIO)         :: cfio
     type(ESMF_Clock)        :: Clock
     type(ESMF_Time)         :: OutTime
     type(ESMF_Time)         :: RefTime
     type(ESMF_TimeInterval) :: TimeStep

     character(len=255) fname,bkgfname,tmpl
     character(len=ESMF_MAXSTR) :: etmpl
     character(len=*), parameter :: only_vars='ps,ts,pblri,pblrf,pblkh,delp,u,v,tv,sphu,ozone,qitot,qltot,qrtot,qstot'
     integer(i_kind) thistime(6)
     integer(i_kind) rc,status,ierr
     integer(i_kind) idim,jdim,kdim,i,j,k,i_dp,i_ps
     integer nymd_ref,nhms_ref
     real(r_kind)    dmodel,dgsi

     stat = 0

!    Define working bundle
!    ---------------------
     WBundle = ESMF_FieldBundleCreate ( name='Work Bundle', rc=status )
               VERIFY_(status)
     call ESMF_FieldBundleSet ( WBundle, grid=GSI_AgcmPertGrid, rc=status )
     VERIFY_(status)

!    Reference date/time
!    -------------------
     nymd_ref = 10000*ibdate(1)+ibdate(2)*100+ibdate(3)
     nhms_ref = 10000*ibdate(4)+ibdate(5)*100
     thistime(1) =     nymd_ref/10000
     thistime(2) = mod(nymd_ref,10000)/100
     thistime(3) = mod(nymd_ref,100)
     thistime(4) =     nhms_ref/10000
     thistime(5) = mod(nhms_ref,10000)/100
     thistime(6) = mod(nhms_ref,100)

!    set reference time
!    ------------------
     call ESMF_TimeSet( RefTime, YY =  thistime(1), &
                                 MM =  thistime(2), &
                                 DD =  thistime(3), &
                                 H  =  thistime(4), &
                                 M  =  thistime(5), &
                                 S  =  thistime(6), rc=status ); VERIFY_(STATUS)

!    Get memory for output bundle by reading inital bkg file
!    NOTE: 1. this needs revision since not general enough
!          2. filanme template wired-in for now
!    -------------------------------------------------------

     call StrTemplate ( bkgfname, GSI_bkg_fname_tmpl, 'GRADS', & 
                        xid=trim(GSI_ExpId), nymd=nymd_ref, nhms=nhms_ref, &
                        stat=status )
          VERIFY_(STATUS)
     call MAPL_CFIORead  ( trim(bkgfname), RefTime, WBundle, &
                           TIME_IS_CYCLIC=.false., verbose=.true., &
                           only_vars=only_vars, rc=status )
!                          noread=.true., only_vars=only_vars, rc=status )
          VERIFY_(status)


     end subroutine set_gsi2pgcm0_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gsi2pgcm0_:  Convert GSI increments to GCM perturbations
!
! !INTERFACE:

      subroutine gsi2pgcm0_ ( nymd, nhms, xx, which, stat, &
                              xp, filename )  ! optionals

! !USES:

      use ESMF

      use GSI_GridCompMod, only: GSI_AgcmPertGrid
      use GSI_GridCompMod, only: GSI_ExpId
      use GSI_GridCompMod, only: GSI_RefTime

      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: assignment(=)

      use m_StrTemplate, only: StrTemplate

      implicit none

! !INPUT PARAMETERS:
                                                                                                                           
      integer(i_kind),    intent(in)    :: nymd   ! date as in YYYYMMDD
      integer(i_kind),    intent(in)    :: nhms   ! time as in HHMMSS
      type(gsi_bundle),   intent(inout) :: xx     ! GSI increment
      character(len=*),   intent(in)    :: which  ! adm or tlm

      character(len=*),optional, intent(in) :: filename ! output filename

! !OUTPUT PARAMETERS:

      type(MAPL_SimpleBundle),optional,intent(out) :: xp
      integer(i_kind),                 intent(out) :: stat

! !DESCRIPTION: Convert GSI increment vector to GEOS-5 perturbation vector
!               (as gsi2pgcm1_, but output GCM perturbation to file)
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  30Sep2007  Todling   Updated interface to putpert.
!  17Jan2010  Todling   Update interface to putpert (add forceflip).
!  15May2010  Todling   Update to use GSI_Bundle
!  12Jul2012  Todling   Now handling ts as well
!  13Mar2014  Todling   Allow output of qi/ql
!  03Mar2017  RT/JJin   Zero out field read in from template file
!  03Mar2020  Todling   Move set/final out into own procedures
!
!EOP
!-----------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'*gsi2pgcm0_'
     character(len=*), parameter :: Iam = myname_
     type(MAPL_SimpleBundle) :: xpert
     type(MAPL_CFIO)         :: cfio
     type(ESMF_Clock)        :: Clock
     type(ESMF_Time)         :: OutTime
     type(ESMF_Time)         :: RefTime
     type(ESMF_TimeInterval) :: TimeStep

     character(len=255) fname,bkgfname,tmpl
     character(len=ESMF_MAXSTR) :: etmpl
     character(len=*), parameter :: only_vars='ps,ts,pblri,pblrf,pblkh,delp,u,v,tv,sphu,ozone,qitot,qltot,qrtot,qstot'
     integer(i_kind) thistime(6)
     integer(i_kind) rc,status,ierr
     integer(i_kind) idim,jdim,kdim,i,j,k,i_dp,i_ps
     integer nymd_ref,nhms_ref
     real(r_kind)    dmodel,dgsi

     stat = 0

!    Make sure read in vars are zeroed out just in case
!    --------------------------------------------------
     call refresh_bundle_ ( WBundle )

!    Link Bundle to Simple-Bundle
!    ----------------------------
     xpert = MAPL_SimpleBundleCreate ( WBundle, rc=status )
             VERIFY_(status)

!    Convert to GSI perturbation vector
!    ----------------------------------
     call gsi2pgcm1_ ( xx, xpert, which, status )
          VERIFY_(status)

!    Build surface pressure perturbation - output purposes
!    -----------------------------------------------------
     i_dp = MAPL_SimpleBundleGetIndex ( xpert, 'delp', 3, rc=status )
     i_ps = MAPL_SimpleBundleGetIndex ( xpert, 'ps'  , 2, rc=status )
     idim = size(xpert%r3(i_dp)%qr4,1)
     jdim = size(xpert%r3(i_dp)%qr4,2)
     kdim = size(xpert%r3(i_dp)%qr4,3)
     xpert%r2(i_ps)%qr4 = zero
     do k=1,kdim
        do j=1,jdim
           do i=1,idim
              xpert%r2(i_ps)%qr4(i,j) = xpert%r2(i_ps)%qr4(i,j) + xpert%r3(i_dp)%qr4(i,j,k)
           end do
       end do
     end do

     dgsi = dot_product(xx,xx)
#ifdef OLD_GEOS_PERT 
!  needs to be done ...
     dmodel = DOT_PRODUCT_TO_BE_DONE(xpert,xpert)
     if(mype==ROOT) write(6,'(2a,1p,e25.18)') trim(myname_), ': magnitude of output vector in model    space ', dmodel
#endif /* OLD_GEOS_PERT */ 
     if(mype==ROOT) write(6,'(2a,1p,e25.18)') trim(myname_), ': magnitude of output vector in analysis space ', dgsi

!    Set ESMF-output date/time
!    -------------------------
     thistime(1) =     nymd/10000
     thistime(2) = mod(nymd,10000)/100
     thistime(3) = mod(nymd,100)
     thistime(4) =     nhms/10000
     thistime(5) = mod(nhms,10000)/100
     thistime(6) = mod(nhms,100)

!   initialize current time
!   -----------------------
     call ESMF_TimeSet( OutTime, YY =  thistime(1), &
                                 MM =  thistime(2), &
                                 DD =  thistime(3), &
                                 H  =  thistime(4), &
                                 M  =  thistime(5), &
                                 S  =  thistime(6), rc=status ); VERIFY_(STATUS)
!    Create clock set at reference time
!    ----------------------------------
     call ESMF_TimeIntervalSet( TimeStep, h=0, m=0, s=0, rc=status )
          VERIFY_(STATUS)
     Clock = ESMF_ClockCreate ( name="pertwrite", timeStep=TimeStep, &
                                startTime=OutTime, rc=status )
             VERIFY_(STATUS)

!    Write out perturbation
!    ----------------------
     mycount = mycount + 1
     if (present(filename)) then
       write(tmpl,'(4a)')      trim(GSI_ExpId), '.', trim(filename), '.%y4%m2%d2_%h2%n2z.nc4'
     else
       write(tmpl,'(4a,i3.3,1a)') trim(GSI_ExpId), '.', trim(fnxgsi), '_', mycount, '.%y4%m2%d2_%h2%n2z.nc4'
     endif
  
     call StrTemplate ( fname, tmpl, 'GRADS', xid=GSI_ExpId, nymd=nymd, nhms=nhms, &
                        stat=status )
          VERIFY_(STATUS)

!    Before writing to the file, we create a CFIO object
!    Here we define a file with reduced resolution
!    ---------------------------------------------------
     call MAPL_CFIOCreate ( cfio, fname, clock, WBundle, rc=status )
          VERIFY_(STATUS)

     call MAPL_CFIOWrite ( cfio, Clock, WBundle, verbose=.true., rc=status )
       if(status/=0)then
           stat = 90
           if(mype==ROOT) print*, trim(myname_), ': Error retrieving perturbation'
           return
       endif

     call MAPL_cfioDestroy ( cfio )

     if (present(xp)) then
!        Should copy simple bundle xpert to xp
         call die(myname_, 'feature not ready ...')
     endif

!    call MAPL_SimpleBundleDestroy ( xpert, rc=STATUS )  ! _RT is this a nullify?? it better be
!        VERIFY_(STATUS)

     end subroutine gsi2pgcm0_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: final_gsi2pgcm0_:  Finalizes put routine
!
! !INTERFACE:

      subroutine final_gsi2pgcm0_ ( stat )

! !USES:

      use ESMF

      use GSI_GridCompMod, only: GSI_AgcmPertGrid
      use GSI_GridCompMod, only: GSI_ExpId
      use GSI_GridCompMod, only: GSI_RefTime

      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: assignment(=)

      use m_StrTemplate, only: StrTemplate

      implicit none

! !INPUT PARAMETERS:
                                                                                                                           
! !OUTPUT PARAMETERS:

      integer(i_kind),                 intent(out) :: stat

! !DESCRIPTION: Finalizes work-bundle used for output.
!
! !REVISION HISTORY:
!
!  03Mar2020  Todling   Broken up from original gsi2pgcm0_
!
!EOP
!-----------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'*final_gsi2pgcm0_'
     character(len=*), parameter :: Iam = myname_
     type(MAPL_SimpleBundle) :: xpert
     type(MAPL_CFIO)         :: cfio
     type(ESMF_Clock)        :: Clock
     type(ESMF_Time)         :: OutTime
     type(ESMF_Time)         :: RefTime
     type(ESMF_TimeInterval) :: TimeStep

     character(len=255) fname,bkgfname,tmpl
     character(len=ESMF_MAXSTR) :: etmpl
     character(len=*), parameter :: only_vars='ps,ts,pblri,pblrf,pblkh,delp,u,v,tv,sphu,ozone,qitot,qltot,qrtot,qstot'
     integer(i_kind) thistime(6)
     integer(i_kind) rc,status,ierr
     integer(i_kind) idim,jdim,kdim,i,j,k,i_dp,i_ps
     integer nymd_ref,nhms_ref
     real(r_kind)    dmodel,dgsi

     stat = 0

!     Release GCM perturbation vector
!     -------------------------------
      call ESMF_FieldBundleDestroy ( WBundle, rc=STATUS )
          VERIFY_(STATUS)


     end subroutine final_gsi2pgcm0_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: refresh_bundle_:  Zero out contents of input bundle

     subroutine refresh_bundle_ (InBundle)
!
! !INTERFACE:

     implicit none
     type(ESMF_FieldBundle)  :: InBundle  

! !DESCRIPTION: Clears out Bundle contents
!
! !REVISION HISTORY:
!
!  03Mar2017  Todling   Add to address bug found by Jianjun Jin related to
!                       how fields not written out by GSI enherit the contents
!                       in originating read-in file from which bundle is created.
!
!EOP

     character(len=ESMF_MAXSTR) :: name
     type(ESMF_Field)        :: Field
     real, pointer, dimension(:,:)   :: ptr2d=>NULL()
     real, pointer, dimension(:,:,:) :: ptr3d=>NULL()      

     integer(i_kind) ifld,nfld,rank,status

     call ESMF_FieldBundleGet(INbundle,FieldCount=nfld, rc=status )
     do ifld=1,nfld
        call ESMF_FieldBundleGet(INbundle, ifld, Field, rc=status )
        call ESMF_FieldGet(Field, NAME=NAME, dimCount = rank, rc=status )
        if (rank==2) then
           call ESMF_FieldGet(Field, farrayPtr=ptr2d, rc=status )
           ptr2d = 0.0
        endif ! <rank=2>
        if (rank==3) then
           call ESMF_FieldGet(Field, farrayPtr=ptr3d, rc=status )
           ptr3d = 0.0
        endif ! <rank=3>
     enddo

     end subroutine refresh_bundle_

!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1    !
!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gsi2pgcm1_:  Convert GSI increments to GCM perturbations
!
! !INTERFACE:

      subroutine gsi2pgcm1_ ( xx, xpert, which, stat )

! !USES:

      use gridmod,       only: lat1, lon1     ! no. lat/lon on subdomain (no buffer)
      use gridmod,       only: lat2, lon2     ! no. lat/lon on subdomain (buffer pnts on ends)
      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: gsi_bundlegetpointer
      use gsi_bundlemod, only: assignment(=)

      implicit none

! !INPUT PARAMETERS:

      type(gsi_bundle),   intent(inout) :: xx    ! GSI increment vector
      character(len=*),   intent(in)    :: which ! adm or tlm

! !OUTPUT PARAMETERS:

      type(MAPL_SimpleBundle),intent(inout) :: xpert ! GCM perturbation vector
      integer(i_kind),        intent(out)   :: stat  ! return error code

! !DESCRIPTION: Converts GSI increments in to ADM/TLM perturbations.
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  21Sep2007  Todling   Handles for O3 and CW.
!  19Nov2008  Todling   Update to use gsi-3d pressure instead of ps.
!  15May2010  Todling   Update to use GSI_Bundle
!  21Feb2011  Todling   Adapt to work with MAPL-SimpleBundle
!  12Jul2012  Todling   Now handling ts as well
!  28Oct2013  Todling   Rename p3d to prse
!  13Mar2014  Todling   Handle for qi/ql
!  07Oct2014  Todling   More flexibility to handle variables not always available
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*gsi2pgcm1_'

      real(r_kind),   allocatable, dimension(:,:,:) :: sub_tv,sub_u,sub_v,sub_q,sub_delp
      real(r_kind),   allocatable, dimension(:,:,:) :: sub_oz,sub_cw
      real(r_kind),   allocatable, dimension(:,:,:) :: sub_qi,sub_ql
      real(r_kind),   allocatable, dimension(:,:,:) :: sub_qr,sub_qs
      real(r_kind),   allocatable, dimension(:,:)   :: sub_ps,sub_ts
      real(r_kind),   allocatable, dimension(:,:)   :: sub_pblri,sub_pblrf,sub_pblkh
      real(r_single), allocatable, dimension(:,:,:) :: aux3d

      integer(i_kind)  i,j,k,ijk,ij
      integer(i_kind)  i_u,i_v,i_t,i_q,i_oz,i_cw,i_prse,i_p,i_tsen,i_ts
      integer(i_kind)  i_pblri,i_pblrf,i_pblkh
      integer(i_kind)  i_qi,i_ql,i_qr,i_qs
      integer(i_kind)  ierr,ierr_adm,istatus,rc,status
      character(len=255) :: whatin

      stat = 0

!     Get poiners to GSI state-vector
!     -------------------------------
      ierr=0; ierr_adm=0
      call gsi_bundlegetpointer(xx,'u' ,  i_u,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'v' ,  i_v,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'tv',  i_t,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'q' ,  i_q,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'oz',  i_oz , istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'prse',i_prse,istatus);ierr_adm=ierr_adm+istatus
      call gsi_bundlegetpointer(xx,'ps',  i_p,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'tsen',i_tsen,istatus);ierr_adm=ierr_adm+istatus
      call gsi_bundlegetpointer(xx,'sst', i_ts  ,istatus);ierr=ierr+istatus
      if(ierr/=0) then
         write(6,*) myname_, ': trouble getting pointers'
         call stop2(999)
      endif
      if (which == 'adm' .and. ierr_adm/=0) then
         write(6,*) myname_, ': trouble getting pointers need for ADM transforms'
         call stop2(999)
      endif

!     Initializes internal arrays
!     ---------------------------
      allocate (sub_tv  (lat2,lon2,nsig), sub_u (lat2,lon2,nsig),&
                sub_v   (lat2,lon2,nsig), sub_q (lat2,lon2,nsig),&
                sub_delp(lat2,lon2,nsig), sub_ps(lat2,lon2), &
                sub_ts  (lat2,lon2), stat=ierr )
        if ( ierr/=0 ) then
            stat = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_)'
            return
        end if
      if ( i_oz>0 ) then
          allocate (sub_oz  (lat2,lon2,nsig), stat=ierr )
            if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_oz)'
                return
            end if
      endif

      if (which == 'adm') then
          call tv_to_tsen_ad(xx%r3(i_t)%q,xx%r3(i_q)%q,xx%r3(i_tsen)%q)
          call getprs_ad    (xx%r2(i_p)%q,xx%r3(i_t)%q,xx%r3(i_prse)%q)
      endif

!     Fill in subdomain arrays
!     ------------------------
      do k=1,nsig
         do j=1,lon2
            do i=1,lat2
               sub_u (i,j,k) = xx%r3(i_u)%q(i,j,k)
               sub_v (i,j,k) = xx%r3(i_v)%q(i,j,k)
               sub_tv(i,j,k) = xx%r3(i_t)%q(i,j,k)
               sub_q (i,j,k) = xx%r3(i_q)%q(i,j,k)
            enddo
         enddo
      enddo
      if ( i_oz>0 ) then
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    sub_oz(i,j,k) = xx%r3(i_oz)%q(i,j,k)
                 enddo
              enddo
           enddo
           sub_oz = sub_oz / PPMV2GpG
      endif
      call gsi_bundlegetpointer(xx,'cw',  i_cw , istatus)
      if ( i_cw>0 ) then
           allocate (sub_cw  (lat2,lon2,nsig), stat=ierr )
           if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_cw)'
                return
           end if
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    sub_cw(i,j,k) = xx%r3(i_cw)%q(i,j,k)
                 enddo
              enddo
           enddo
      endif
      call gsi_bundlegetpointer(xx,'qi',  i_qi,  istatus)
      if ( i_qi>0 ) then
           allocate (sub_qi  (lat2,lon2,nsig), stat=ierr )
           if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_qi)'
                return
           end if
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    sub_qi(i,j,k) = xx%r3(i_qi)%q(i,j,k)
                 enddo
              enddo
           enddo
      endif
      call gsi_bundlegetpointer(xx,'ql',  i_ql,  istatus)
      if ( i_ql>0 ) then
           allocate (sub_ql  (lat2,lon2,nsig), stat=ierr )
           if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_ql)'
                return
           end if
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    sub_ql(i,j,k) = xx%r3(i_ql)%q(i,j,k)
                 enddo
              enddo
           enddo
      endif
      call gsi_bundlegetpointer(xx,'qr',  i_qr,  istatus)
      if ( i_qr>0 ) then
           allocate (sub_qr  (lat2,lon2,nsig), stat=ierr )
           if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_qr)'
                return
           end if
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    sub_qr(i,j,k) = xx%r3(i_qr)%q(i,j,k)
                 enddo
              enddo
           enddo
      endif
      call gsi_bundlegetpointer(xx,'qs',  i_qs,  istatus)
      if ( i_qs>0 ) then
           allocate (sub_qs  (lat2,lon2,nsig), stat=ierr )
           if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_qs)'
                return
           end if
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    sub_qs(i,j,k) = xx%r3(i_qs)%q(i,j,k)
                 enddo
              enddo
           enddo
      endif
      if (i_prse>0) then
         sub_ps = zero  ! will derive ps from delp(prse)
      else
         do j=1,lon2
            do i=1,lat2
               sub_ps(i,j) = xx%r2(i_p)%q(i,j)
            enddo
         enddo
      endif
      if (i_ts>0) then
         do j=1,lon2
            do i=1,lat2
               sub_ts(i,j) = xx%r2(i_ts)%q(i,j)
            enddo
         enddo
      endif
      call gsi_bundlegetpointer(xx,'pblri', i_pblri, istatus)
      if (i_pblri>0) then
         allocate (sub_pblri(lat2,lon2), stat=ierr )
         if ( ierr/=0 ) then
             stat = 91
             if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_pblri)'
             return
         end if
         do j=1,lon2
            do i=1,lat2
               sub_pblri(i,j) = xx%r2(i_pblri)%q(i,j)
            enddo
         enddo
      endif
      call gsi_bundlegetpointer(xx,'pblrf', i_pblrf, istatus)
      if (i_pblrf>0) then
         allocate (sub_pblrf(lat2,lon2), stat=ierr )
         if ( ierr/=0 ) then
             stat = 91
             if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_pblrf)'
             return
         end if
         do j=1,lon2
            do i=1,lat2
               sub_pblrf(i,j) = xx%r2(i_pblrf)%q(i,j)
            enddo
         enddo
      endif
      call gsi_bundlegetpointer(xx,'pblkh', i_pblkh, istatus)
      if (i_pblkh>0) then
         allocate (sub_pblkh(lat2,lon2), stat=ierr )
         if ( ierr/=0 ) then
             stat = 91
             if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_pblkh)'
             return
         end if
         do j=1,lon2
            do i=1,lat2
               sub_pblkh(i,j) = xx%r2(i_pblkh)%q(i,j)
            enddo
         enddo
      endif

!     Calculate perturbation delp
!     ---------------------------
      if (which == 'adm') then
          call delp2ps_ad_ ( kPa_per_Pa ) 
      else if (which == 'tlm') then
          call ps2delp_    ( Pa_per_kPa )
      else
          call die ( myname_, ': invalid option' )
      endif

!     Gather from GSI subdomains/Scatter to GCM
!     -----------------------------------------
          i_u = MAPL_SimpleBundleGetIndex ( xpert, 'u'   , 3, rc=ierr )
          i_v = MAPL_SimpleBundleGetIndex ( xpert, 'v'   , 3, rc=ierr )
          i_t = MAPL_SimpleBundleGetIndex ( xpert, 'tv'  , 3, rc=ierr )
          i_p = MAPL_SimpleBundleGetIndex ( xpert, 'delp', 3, rc=ierr )
          i_q = MAPL_SimpleBundleGetIndex ( xpert, 'sphu', 3, rc=ierr )
          i_ts= MAPL_SimpleBundleGetIndex ( xpert, 'ts'  , 2, rc=ierr )
          call gsi2pert_ ( sub_u,    xpert%r3(i_u)%qr4,  ierr )
          call gsi2pert_ ( sub_v,    xpert%r3(i_v)%qr4,  ierr )
          call gsi2pert_ ( sub_tv,   xpert%r3(i_t)%qr4,  ierr )
          call gsi2pert_ ( sub_delp, xpert%r3(i_p)%qr4,  ierr )
          call gsi2pert_ ( sub_q ,   xpert%r3(i_q)%qr4,  ierr )
      if (i_oz>0) then
          i_oz = MAPL_SimpleBundleGetIndex ( xpert, 'ozone', 3, rc=ierr )
          call gsi2pert_ ( sub_oz,   xpert%r3(i_oz)%qr4, ierr )
      endif
      if (i_cw>0) then ! _RT for now put increment of cw onto ql increment
          i_cw = MAPL_SimpleBundleGetIndex ( xpert, 'qltot', 3, rc=ierr )
          call gsi2pert_ ( sub_cw,   xpert%r3(i_cw)%qr4, ierr )
      endif
      if (i_qi>0) then
          i_qi = MAPL_SimpleBundleGetIndex ( xpert, 'qitot', 3, rc=ierr )
          call gsi2pert_ ( sub_qi,   xpert%r3(i_qi)%qr4, ierr )
      endif
      if (i_ql>0) then
          i_ql = MAPL_SimpleBundleGetIndex ( xpert, 'qltot', 3, rc=ierr )
          call gsi2pert_ ( sub_ql,   xpert%r3(i_ql)%qr4, ierr )
      endif
      if (i_qr>0) then
          i_qr = MAPL_SimpleBundleGetIndex ( xpert, 'qrtot', 3, rc=ierr )
          call gsi2pert_ ( sub_qr,   xpert%r3(i_qr)%qr4, ierr )
      endif
      if (i_qs>0) then
          i_qs = MAPL_SimpleBundleGetIndex ( xpert, 'qstot', 3, rc=ierr )
          call gsi2pert_ ( sub_qs,   xpert%r3(i_qs)%qr4, ierr )
      endif
        if ( ierr/=0 ) then
            stat = 98
            if(mype==ROOT) print*, trim(myname_), ': unfinished convertion ...'
            return
        end if

!     Handle 2d fields
!     ----------------
          sub_q=zero
          sub_q(:,:,1)=sub_ts          ! copy to level 1
          allocate(aux3d(size(xpert%r3(i_q)%qr4,1),size(xpert%r3(i_q)%qr4,2),size(xpert%r3(i_q)%qr4,3)))
          call gsi2pert_ ( sub_q, aux3d,  ierr )
          xpert%r2(i_ts)%qr4=aux3d(:,:,nsig) ! level nsig will contain result give vertical swap
          deallocate(aux3d)

          if (i_pblri>0) then
             i_pblri = MAPL_SimpleBundleGetIndex ( xpert, 'pblri'  , 2, rc=ierr )
             sub_q=zero
             sub_q(:,:,1)=sub_pblri       ! copy to level 1
             allocate(aux3d(size(xpert%r3(i_q)%qr4,1),size(xpert%r3(i_q)%qr4,2),size(xpert%r3(i_q)%qr4,3)))
             call gsi2pert_ ( sub_q, aux3d,  ierr )
             xpert%r2(i_pblri)%qr4=aux3d(:,:,nsig) ! level nsig will contain result give vertical swap
             deallocate(aux3d)
          end if
          if (i_pblrf>0) then
             i_pblrf = MAPL_SimpleBundleGetIndex ( xpert, 'pblrf'  , 2, rc=ierr )
             sub_q=zero
             sub_q(:,:,1)=sub_pblrf       ! copy to level 1
             allocate(aux3d(size(xpert%r3(i_q)%qr4,1),size(xpert%r3(i_q)%qr4,2),size(xpert%r3(i_q)%qr4,3)))
             call gsi2pert_ ( sub_q, aux3d,  ierr )
             xpert%r2(i_pblrf)%qr4=aux3d(:,:,nsig) ! level nsig will contain result give vertical swap
             deallocate(aux3d)
          end if
          if (i_pblkh>0) then
             i_pblkh = MAPL_SimpleBundleGetIndex ( xpert, 'pblkh'  , 2, rc=ierr )
             sub_q=zero
             sub_q(:,:,1)=sub_pblkh       ! copy to level 1
             allocate(aux3d(size(xpert%r3(i_q)%qr4,1),size(xpert%r3(i_q)%qr4,2),size(xpert%r3(i_q)%qr4,3)))
             call gsi2pert_ ( sub_q, aux3d,  ierr )
             xpert%r2(i_pblkh)%qr4=aux3d(:,:,nsig) ! level nsig will contain result give vertical swap
             deallocate(aux3d)
          end if

      deallocate (sub_tv, sub_u, sub_v, sub_q, sub_delp, sub_ps, sub_ts, stat=ierr )
        if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_)'
            return
        end if
      if ( i_oz>0 ) then
          deallocate (sub_oz, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_oz)'
                return
            end if
      endif
      if ( allocated(sub_cw) ) then
            deallocate (sub_cw, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_cw)'
                return
            end if
      endif
      if ( allocated(sub_qi) ) then
          deallocate (sub_qi, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_qi)'
                return
            end if
      endif
      if ( allocated(sub_ql) ) then
          deallocate (sub_ql, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_ql)'
                return
            end if
      endif
      if ( allocated(sub_qr) ) then
          deallocate (sub_qr, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_qr)'
                return
            end if
      endif
      if ( allocated(sub_qs) ) then
          deallocate (sub_qs, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_qs)'
                return
            end if
       endif
       if ( allocated(sub_pblri) ) then
          deallocate (sub_pblri, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_pblri)'
                return
            end if
       endif
       if ( allocated(sub_pblrf) ) then
          deallocate (sub_pblrf, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_pblrf)'
                return
            end if
       endif
       if ( allocated(sub_pblkh) ) then
          deallocate (sub_pblkh, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_pblkh)'
                return
            end if
       endif



      CONTAINS

      subroutine ps2delp_ ( alpha )  ! p2delp
! tlm-only
      implicit none
      real(r_kind), intent(in) :: alpha
      real(r_kind) bkweight
      real(r_kind), allocatable, dimension(:,:,:) :: sub_aux
      integer(i_kind) i,j,k,ij,ijk
        allocate ( sub_aux(lat2,lon2,nsig+1) )
        if ( i_prse>0 ) then
           do k=1,nsig+1
              do j=1,lon2
                 do i=1,lat2
                    sub_aux(i,j,k) = alpha * xx%r3(i_prse)%q(i,j,k)
                 enddo
              enddo
           enddo
        else
           do k=1,nsig+1
              do j=1,lon2
                 do i=1,lat2
                    sub_aux(i,j,k) = alpha*bk5(k)*xx%r2(i_p)%q(i,j)
                 enddo
              enddo
           enddo
        endif
        do k=1,nsig
           do j=1,lon2
              do i=1,lat2
                 sub_delp(i,j,k) = sub_aux(i,j,k) - sub_aux(i,j,k+1)
              enddo
           enddo
        enddo
        deallocate ( sub_aux )
      end subroutine ps2delp_

      subroutine delp2ps_ad_ ( alpha )
! inverse-adm
      implicit none
      real(r_kind), intent(in) :: alpha
      real(r_kind), allocatable :: sub_aux(:,:,:)
      integer(i_kind) i,j,k
        allocate ( sub_aux(lat2,lon2,nsig+1) )
        sub_aux=zero
        do k=1,nsig+1
           do j=1,lon2
              do i=1,lat2
                 sub_aux(i,j,k) = sub_aux(i,j,k) + xx%r3(i_prse)%q(i,j,k)
                 xx%r3(i_prse)%q(i,j,k) = zero
              end do
           end do
        end do
        do k=1,nsig
           do j=1,lon2
              do i=1,lat2
                 sub_aux(i,j,k+1) = sub_aux(i,j,k+1) + sub_aux(i,j,k)
                 sub_delp(i,j,k) = sub_delp(i,j,k) + alpha * sub_aux(i,j,k)
                 sub_aux(i,j,k) = zero
              end do
           end do
        end do
        deallocate ( sub_aux )
        do k=nsig,1,-1
           do j=1,lon2
              do i=1,lat2
                 sub_delp(i,j,k) = sub_delp(i,j,k) + alpha * xx%r2(i_p)%q(i,j)
              end do
           end do
        end do
        xx%r2(i_p)%q=zero
      end subroutine delp2ps_ad_

      subroutine gsi2pert_ ( sub, fld, stat_ )

      real(r_kind),   intent(in)  :: sub(:,:,:)
      real(r_single), intent(inout) :: fld(:,:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*gsi2pert_'

      real(r_kind), allocatable :: fldsm(:,:)
      real(r_kind), allocatable :: work4d(:,:,:,:)   ! auxliar 4d array
      real(r_kind), allocatable :: work3d(:,:,:)     ! auxliar 3d array
      real(r_kind), allocatable :: work(:)

      integer(i_kind) imr,jnp,nl
      integer(i_kind) i,j,k,mm1,ierr

      mm1 = mype+1 
      stat_ = 0

      allocate ( work3d(nlon,nlat,nsig), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work3d)'
            return
        end if
      allocate ( work(max(iglobal,itotsub)), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work)'
            return
        end if
      allocate ( fldsm(lat1*lon1,nsig), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(fldsm)'
            return
        end if

!     Strip off boundary points from subdomains
!     -----------------------------------------
      call strip(sub,fldsm,nsig)

!     Gather GSI perturbations to root processor
!     ------------------------------------------
      do k=1,nsig
         call mpi_gatherv(fldsm(1,k),ijn(mm1),mpi_rtype,&
              work,ijn,displs_g,mpi_rtype,&
              ROOT,mpi_comm_world,ierr)
         if (mype==ROOT) then
            call reorder12(work,work3d(:,:,k))
         endif
      end do

      deallocate ( fldsm, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(fldsm)'
            return
        end if
      deallocate ( work, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work)'
            return
        end if

      call mygetdims_(imr,jnp,nl)
      allocate ( work4d(imr,jnp,nl,1), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work4d)'
            return
        end if

!     Flip horizontal and vertical
!     ----------------------------
      if ( mype==ROOT ) then
              call hflip3_ ( work3d, imr,jnp,nl )
              if (imr/=nlon .or. jnp/=nlat ) then
                  if (which=='adm') then
                      work4d = zero
                      call interpack_terpv_ad ( imr,jnp,nl,work4d(:,:,:,1),work4d(:,:,:,1), nlon,nlat,nsig,work3d, ierr )
                  else if (which=='tlm') then
                      work4d = zero
                      call interpack_terpv    ( nlon,nlat,nsig,work3d, imr,jnp,nl,work4d(:,:,:,1), ierr )
                  else
                      call die ( myname_,': invalid option' )
                  endif
              else
                  work4d(:,:,:,1) = work3d(:,:,:)
              endif
           call SwapV_ ( work4d(:,:,:,1) )
      endif

!     Scatter perturbations to GCM decomposition
!     ------------------------------------------
      call myscatter_ ( which, imr, jnp, nl, work4d(:,:,:,1), fld )

!     Swap work memory
!     ----------------
      deallocate ( work4d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(work4d)'
            return
        end if

      deallocate ( work3d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work3d)'
            return
        end if

      end subroutine gsi2pert_

      end subroutine gsi2pgcm1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pgcm2gsi0_:  Convert gcm-adm/tlm vector to gsi vector
!
! !INTERFACE:

      subroutine pgcm2gsi0_ ( xx, which, stat, &
                              filename, nymd_in, nhms_in )

! !USES:

      use kinds, only: i_kind,r_kind
      use constants, only: zero

      use ESMF

      use GSI_GridCompMod, only: GSI_AgcmPertGrid
      use GSI_GridCompMod, only: GSI_RefTime
      use GSI_GridCompMod, only: GSI_FcsTime
      use GSI_GridCompMod, only: GSI_ExpId

      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: assignment(=)

      use gsi_4dvar, only: tau_fcst

      implicit none

! !INPUT PARAMETERS:

      character(len=*),           intent(in)  :: which    ! adm or tlm
      character(len=*), optional, intent(in)  :: filename ! name of file w/ perturbation

! !OUTPUT PARAMETERS:

      type(gsi_bundle),         intent(inout) :: xx       ! GSI increment

      integer(i_kind),            intent(out) :: stat

      integer(i_kind),optional,   intent(out) :: nymd_in  ! Date read in
      integer(i_kind),optional,   intent(out) :: nhms_in  ! Time read in

! !DESCRIPTION: Convert GEOS-5 perturbation vector to GSI increment vector
!               (as pgcm2gsi1_, but reads GCM perturbation from file)
!
! !REMARKS:
!
!  31Oct2007  Todling   De-activated interpolation capability; not fully coded yet
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  17Jul2007  Todling   Add ability to read in perturbation at diff resolution
!  17Jul2008  Todling   Add filename as optional argument
!  15May2010  Todling   Update to use GSI_Bundle
!  20Feb2011  Todling   Adapt to work with MAPL-SimpleBundle
!  22Jun2012  Todling   Change handling of filename 
!  12Jul2012  Todling   Now handling ts as well
!  10Oct2016  M.Kim     Handle ql, qi, qr, and qs 
!  21May2017  Todling   Knobs to differentiate among different input vectors
!  06May2020  Todling   Considerable rewrite to allow reading pert into vector
!                       of arbitrary dim (not associated w/ analysis increment) 
!
!EOP
!-----------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'*pgcm2gsi0_'
     character(len=*), parameter :: Iam = myname_
     character(len=*), parameter :: fname_def ='fsens.eta.nc4' ! full forecast sensitivity from ADM integration
     character(len=*), parameter :: sens_vars='ts,delp,u,v,tv,sphu,ozone'
     character(len=*), parameter :: xinc_vars='pblri,pblrf,pblkh,ts,ps,u,v,tv,sphu,ozone,qitot,qltot,qrtot,qstot' ! not sure if this 
                                                                  ! subroutine is only used for fcst sensitivity calculation,
                                                                  ! so add PBLH here for now
     character(len=256) :: fname, ffname(1)
     character(len=256) :: my_vars
     character(len=30)  :: GSIGRIDNAME
     type(ESMF_VM)           :: VM
     type(ESMF_FieldBundle),pointer  :: WBundle(:)
     type(MAPL_SimpleBundle) :: xpert
     type(ESMF_Grid)         :: GSI_iGrid
     type(ESMF_Time)         :: This_RefTime
     type(LatLonGridFactory) :: ll_factory

     integer(i_kind) myimr,myjnp,mynl,mync
     integer(i_kind) ierr,yy,mm,dd,hh,mn,sec
     integer         rc,status
     integer         THIS_TIME(6)
     integer(i_kind) my_nlat,my_nlon
     real(r_kind)    dmodel,dgsi

     stat = 0
     if (nlat_set>0 .and. nlon_set>0) then
        my_nlon=nlon_set
        my_nlat=nlat_set
     else ! this default is somewhat dangerous and could trip the unadvised,
          ! but it avoids having to update older get calls w/ set/unset
        my_nlon=nlon
        my_nlat=nlat
     endif

!    Define working bundle
!    ---------------------
     call ESMF_VMGetCurrent(vm=VM, rc=STATUS)
     VERIFY_(STATUS)
     call MAPL_DefGridName (my_nlon,my_nlat,GSIGRIDNAME,MAPL_am_I_root())

     ll_factory = LatLonGridFactory(grid_name=GSIGRIDNAME, nx=nxpe, ny=nype, &
                 im_world=my_nlon, jm_world=my_nlat, lm = nsig, &
                 pole='PC', dateline='DC', rc=status)
     VERIFY_(status)
     GSI_iGrid = grid_manager%make_grid(ll_factory,rc=status)
     VERIFY_(STATUS)

     allocate(WBundle(1))
     WBundle(1) = ESMF_FieldBundleCreate ( name='Work Bundle', rc=status )
                  VERIFY_(status)
     call ESMF_FieldBundleSet ( WBundle(1), grid=GSI_iGrid, rc=status )
     VERIFY_(status)

!    Set file to be read
!    -------------------
     my_vars=trim(sens_vars)
     if(present(filename)) then
        fname = trim(filename)
        if(trim(fname)=='NULL') fname = GSI_fsens_fname_tmpl
        if(trim(fname)=='A')    fname = GSI_ferrA_fname_tmpl
        if(trim(fname)=='B')    fname = GSI_ferrB_fname_tmpl
        if(trim(fname(1:4))=='xinc') then
           write(fname,'(4a)')  trim(GSI_ExpId), '.', trim(filename), '.eta.%y4%m2%d2_%h2%n2z.nc4'
           my_vars= trim(xinc_vars)
        endif
     else 
        fname=fname_def
     endif

     if (nymd_set>0 .and. nhms_set>=0 ) then

!       Fill in filename template
!       -------------------------
        call StrTemplate ( ffname(1), fname, 'GRADS', &
                           xid=trim(GSI_ExpId), nymd=nymd_set, nhms=nhms_set, &
                           stat=status )
                           VERIFY_(STATUS)

         THIS_TIME(1) =     nymd_set/10000
         THIS_TIME(2) = mod(nymd_set,10000)/100
         THIS_TIME(3) = mod(nymd_set,100)
         THIS_TIME(4) =     nhms_set/10000
         THIS_TIME(5) = mod(nhms_set,10000)/100
         THIS_TIME(6) = mod(nhms_set,100)
         call ESMF_TimeSet(  This_RefTime, YY =  THIS_TIME(1), &
                                           MM =  THIS_TIME(2), &
                                           DD =  THIS_TIME(3), &
                                           H  =  THIS_TIME(4), &
                                           M  =  THIS_TIME(5), &
                                           S  =  THIS_TIME(6), rc=status ); VERIFY_(STATUS)

!        Read perturbation from file into ESMF Bundle
!        --------------------------------------------
         call MAPL_CFIORead  ( ffname(1), This_RefTime, WBundle(1), only_vars=trim(my_vars), &
                               TIME_IS_CYCLIC=.false., verbose=.false., &
                               doParallel=.true., gsiMode=gsiMode, rc=status )
         VERIFY_(status)


     else

!       When tau_fcst is set, redefine reference time of field to
!       be read in as that determined as current time ticked 
!       forward by tau_fcst
!       ---------------------------------------------------------
        if ( tau_fcst>0 ) then

!          Read perturbation from file into ESMF Bundle
!          --------------------------------------------
           call MAPL_CFIORead  ( fname, GSI_FcsTime, WBundle(1), only_vars=trim(my_vars), &
                                 TIME_IS_CYCLIC=.false., verbose=.false., &
                                 doParallel=.true., gsiMode=gsiMode, rc=status )
           VERIFY_(status)

        else

!          Read perturbation from file into ESMF Bundle
!          --------------------------------------------
           call MAPL_CFIORead  ( fname, GSI_RefTime, WBundle(1), only_vars=trim(my_vars), &
                                 TIME_IS_CYCLIC=.false., verbose=.false., &
                                 doParallel=.true., gsiMode=gsiMode, rc=status )
           VERIFY_(status)

        endif
     endif

!    Link Bundle to Simple-Bundle
!    ----------------------------
     xpert = MAPL_SimpleBundleCreate ( WBundle(1), rc=status )
             VERIFY_(status)

!    Simply convert in GSI fields
!    ----------------------------
     call pgcm2gsi1_ ( xpert, xx, which, stat, jgradf=.true. )

    dgsi = dot_product(xx,xx)
#ifdef OLD_GEOS_PERT 
! needs to be done ...
    dmodel = DOT_PRODUCT_TO_BE_DONE(xpert,xpert)
    if(mype==ROOT) write(6,'(2a,1p,e25.18)') trim(myname_), ': magnitude of input vector in model    space ', dmodel
#endif /* OLD_GEOS_PERT */ 
    if(mype==ROOT) write(6,'(2a,1p,e25.18)') trim(myname_), ': magnitude of input vector in analysis space ', dgsi

!    Need to convert ESMF time to conventional time
!    ----------------------------------------------
     if (present(nymd_in).and.present(nhms_in)) then
         call ESMF_TimeGet(GSI_RefTime, yy=YY, mm=MM, dd=DD, h=HH, m=MN, s=SEC, rc=status)
         nymd_in = 10000*yy + 100*mm + dd
         nhms_in = 10000*hh + 100*mn + sec
     endif

     call MAPL_SimpleBundleDestroy ( xpert, rc=STATUS )  ! _RT is this a nullify?? it better be
          VERIFY_(STATUS)
     call ESMF_FieldBundleDestroy ( WBundle(1), rc=STATUS )
          VERIFY_(STATUS)
     deallocate(WBundle)

     end subroutine pgcm2gsi0_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pgcm2gsi1_:  Convert gcm-adm/tlm vector to gsi vector
!
! !INTERFACE:

      subroutine pgcm2gsi1_ ( xpert, xx, which, stat, jgradf )

! !USES:

      use constants, only: one,zero
      use gsi_bundlemod,only: gsi_bundle     ! GSI state vector
      use gsi_bundlemod,only: gsi_bundlegetpointer
      use gsi_bundlemod,only: assignment(=)

      use gridmod,       only: lat2, lon2    ! no. lat/lon on subdomain (buffer pnts on ends)
      implicit none

! !INPUT PARAMETERS:

      type(MAPL_SimpleBundle), intent(inout)  :: xpert  ! GCM perturbation vector 
      character(len=*),           intent(in)  :: which  ! adm or tlm
      logical, optional,          intent(in)  :: jgradf ! specify when input is forecast (error) gradient

! !OUTPUT PARAMETERS:

      type(gsi_bundle),         intent(inout) :: xx     ! GSI increment

      integer(i_kind),            intent(out) :: stat

! !DESCRIPTION: Convert GEOS-5 perturbation vector to GSI increment vector
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  21Sep2007  Todling   Handles for O3 and CW.
!  22Feb2008  Todling   Handle for forecast (error) gradient vector.
!  19Nov2008  Todling   Update to use gsi-3d pressure instead of ps.
!  15May2010  Todling   Update to use GSI_Bundle
!  20Feb2011  Todling   Adapt to work with MAPL-SimpleBundle
!  12Jul2012  Todling   Now handling ts as well
!  28Oct2013  Todling   Rename p3d to prse
!  12Jul2016  Todling   Slightly revisit handle for tracers
!  08Oct2016  M.Kim     Allow for input of qr/qs
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*pgcm2gsi1_'

      real(r_kind),  allocatable, dimension(:,:,:) :: sub_u,sub_v,sub_delp,sub_q,sub_tv
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_oz,sub_cw
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_ql,sub_qi,sub_qr,sub_qs
      real(r_kind),  allocatable, dimension(:,:)   :: sub_ps,sub_ts
      real(r_kind),  allocatable, dimension(:,:)   :: sub_pblri,sub_pblrf,sub_pblkh
      real(r_kind),  allocatable, dimension(:,:,:) :: aux3d

      character(len=256) :: failvars
      integer(i_kind) i,j,k,ij,ijk,id,ng_d,ng_s
      integer(i_kind) i_u,i_v,i_t,i_q,i_oz,i_cw,i_prse,i_p,i_tsen,i_ts,i_ql,i_qi,i_qr,i_qs,& 
                      i_pblri,i_pblrf,i_pblkh
      integer(i_kind) slat2,slon2
      integer(i_kind) ierr,istatus,status
      logical         scaleit

      scaleit = .true.  ! default: scale input vector as original var from G-5 GCM
      stat = 0
      if ( present(jgradf) ) then
           if(jgradf) scaleit = .false. ! input vector is a gradient, don''t scale vars
      endif

      if (lon2_set>0.and.lat2_set>0) then
         slon2 = lon2_set
         slat2 = lat2_set
      else ! this default is somewhat dangerous and could trip the unadvised,
           ! but it avoids having to update older get calls w/ set/unset
         slon2 = lon2
         slat2 = lat2
      endif

      ierr=0

!     Gather from GCM/Scatter to GSI subdomains
!     -----------------------------------------
      i_u = MAPL_SimpleBundleGetIndex ( xpert, 'u'   , 3, rc=status )
      i_v = MAPL_SimpleBundleGetIndex ( xpert, 'v'   , 3, rc=status )
      i_t = MAPL_SimpleBundleGetIndex ( xpert, 'tv'  , 3, rc=status )
      i_q = MAPL_SimpleBundleGetIndex ( xpert, 'sphu', 3, rc=status )
      i_ts= MAPL_SimpleBundleGetIndex ( xpert, 'ts'  , 2, rc=status )
      if (i_u>0) then
         allocate (sub_u(slat2,slon2,nsig), stat=istatus ); ierr=ierr+istatus
         call pert2gsi_ ( xpert%r3(i_u)%qr4, sub_u   , ierr )
      endif
      if (i_v>0) then
         allocate (sub_v(slat2,slon2,nsig), stat=istatus ); ierr=ierr+istatus
         call pert2gsi_ ( xpert%r3(i_v)%qr4, sub_v   , ierr )
      endif
      if (i_t>0) then
         allocate (sub_tv  (slat2,slon2,nsig), stat=istatus ); ierr=ierr+istatus
         call pert2gsi_ ( xpert%r3(i_t)%qr4, sub_tv  , ierr )
      endif
      if (i_q>0) then
         allocate (sub_q   (slat2,slon2,nsig), stat=istatus ); ierr=ierr+istatus
         call pert2gsi_ ( xpert%r3(i_q)%qr4, sub_q   , ierr )
      endif
      if(which/='inc') then
        i_p = MAPL_SimpleBundleGetIndex ( xpert, 'delp', 3, rc=status )
        if(i_p>0) then
           allocate(sub_delp(slat2,slon2,nsig))
           call pert2gsi_ ( xpert%r3(i_p)%qr4, sub_delp, ierr )
        endif
      endif
      i_oz = MAPL_SimpleBundleGetIndex ( xpert, 'ozone', 3, quiet=.true., rc=status )
      if (i_oz>0) then
         allocate (sub_oz  (slat2,slon2,nsig), stat=ierr )
         if ( ierr/=0 ) then
            stat = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_oz)'
            return
         end if
         call pert2gsi_ ( xpert%r3(i_oz)%qr4, sub_oz , ierr )
      endif
      i_cw = MAPL_SimpleBundleGetIndex ( xpert, 'qctot', 3, quiet=.true., rc=status )
      if (i_cw>0) then
          allocate (sub_cw  (slat2,slon2,nsig), stat=ierr )
          if ( ierr/=0 ) then
              stat = 91
              if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_cw)'
              return
          end if
          call pert2gsi_ ( xpert%r3(i_cw)%qr4, sub_cw , ierr )
      endif
      i_ql = MAPL_SimpleBundleGetIndex ( xpert, 'qltot', 3, quiet=.true., rc=status )
      if (i_ql>0) then
          allocate (sub_ql  (slat2,slon2,nsig), stat=ierr )
          if ( ierr/=0 ) then
              stat = 91
              if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_ql)'
              return
          end if
          call pert2gsi_ ( xpert%r3(i_ql)%qr4, sub_ql , ierr )
      endif
      i_qi = MAPL_SimpleBundleGetIndex ( xpert, 'qitot', 3, quiet=.true., rc=status )
      if (i_qi>0) then
          allocate (sub_qi  (slat2,slon2,nsig), stat=ierr )
          if ( ierr/=0 ) then
              stat = 91
              if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_qi)'
              return
          end if
          call pert2gsi_ ( xpert%r3(i_qi)%qr4, sub_qi , ierr )
      endif
      i_qr = MAPL_SimpleBundleGetIndex ( xpert, 'qrtot', 3, quiet=.true., rc=status )
      if (i_qr>0) then
          allocate (sub_qr  (slat2,slon2,nsig), stat=ierr )
          if ( ierr/=0 ) then
              stat = 91
              if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_qr)'
              return
          end if
          call pert2gsi_ ( xpert%r3(i_qr)%qr4, sub_qr , ierr )
      endif
      i_qs = MAPL_SimpleBundleGetIndex ( xpert, 'qstot', 3, quiet=.true., rc=status )
      if (i_qs>0) then
          allocate (sub_qs  (slat2,slon2,nsig), stat=ierr )
          if ( ierr/=0 ) then
              stat = 91
              if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_qs)'
              return
          end if
          call pert2gsi_ ( xpert%r3(i_qs)%qr4, sub_qs , ierr )
      endif
      if (i_ts>0) then
         allocate (sub_ts(slat2,slon2), stat=istatus ); ierr=ierr+istatus
         call pert2gsi2d_ ( xpert%r2(i_ts)%qr4, sub_ts, ierr )
      endif
      i_pblri= MAPL_SimpleBundleGetIndex ( xpert, 'pblri'  , 2, rc=status )
      if (i_pblri>0) then
         allocate (sub_pblri(slat2,slon2), stat=istatus ); ierr=ierr+istatus
         call pert2gsi2d_ ( xpert%r2(i_ts)%qr4, sub_pblri, ierr )
      endif
      i_pblrf= MAPL_SimpleBundleGetIndex ( xpert, 'pblrf'  , 2, rc=status )
      if (i_pblrf>0) then
         allocate (sub_pblrf(slat2,slon2), stat=istatus ); ierr=ierr+istatus
         call pert2gsi2d_ ( xpert%r2(i_ts)%qr4, sub_pblrf, ierr )
      endif
      i_pblkh= MAPL_SimpleBundleGetIndex ( xpert, 'pblkh'  , 2, rc=status )
      if (i_pblkh>0) then
         allocate (sub_pblkh(slat2,slon2), stat=istatus ); ierr=ierr+istatus
         call pert2gsi2d_ ( xpert%r2(i_ts)%qr4, sub_pblkh, ierr )
      endif

      if (which=='inc') then
         i_p= MAPL_SimpleBundleGetIndex ( xpert, 'ps', 2, rc=status )
         if (i_p>0) then
            allocate (sub_ps(slat2,slon2), stat=istatus ); ierr=ierr+istatus
            call pert2gsi2d_ ( xpert%r2(i_p)%qr4, sub_ps, ierr )
         endif
      endif

      if ( ierr/=0 ) then
          stat = 99
          if(mype==ROOT) print*, trim(myname_), ': unfinished convertion ...'
          return
      end if

!     Get poiners to GSI state-vector
!     -------------------------------
      ierr=0; failvars=""
      call gsi_bundlegetpointer(xx,'u' ,  i_u,   istatus)
      if (istatus/=0) then
          istatus=0
          call gsi_bundlegetpointer(xx,'sf', i_u,istatus);ierr=ierr+istatus
          if(ierr/=0) failvars = trim(failvars)//"sf:"
      endif
      call gsi_bundlegetpointer(xx,'v' ,  i_v,   istatus)
      if (istatus/=0) then
          istatus=0
          call gsi_bundlegetpointer(xx,'vp', i_v,istatus);ierr=ierr+istatus
          if(ierr/=0) failvars = trim(failvars)//"vp:"
      endif
      call gsi_bundlegetpointer(xx,'tv',  i_t,   istatus)
      if (istatus/=0) then
          istatus=0
          call gsi_bundlegetpointer(xx,'t', i_t,istatus);ierr=ierr+istatus
          if(ierr/=0) failvars = trim(failvars)//"t:"
      endif
      call gsi_bundlegetpointer(xx,'q' ,  i_q,   istatus);ierr=ierr+istatus
          if(ierr/=0) failvars = trim(failvars)//"q:"
      call gsi_bundlegetpointer(xx,'oz',  i_oz , istatus);ierr=ierr+istatus
          if(ierr/=0) failvars = trim(failvars)//"oz:"
      call gsi_bundlegetpointer(xx,'ps',  i_p,   istatus);ierr=ierr+istatus
          if(ierr/=0) failvars = trim(failvars)//"ps:"
      if (which /= 'inc') then
         call gsi_bundlegetpointer(xx,'cw',  i_cw , istatus);ierr=ierr+istatus
          if(ierr/=0) failvars = trim(failvars)//"cw:"
         call gsi_bundlegetpointer(xx,'prse',i_prse,istatus);ierr=ierr+istatus
          if(ierr/=0) failvars = trim(failvars)//"prse:"
         call gsi_bundlegetpointer(xx,'tsen',i_tsen,istatus);ierr=ierr+istatus
          if(ierr/=0) failvars = trim(failvars)//"tsen:"
      endif
      call gsi_bundlegetpointer(xx,'sst', i_ts  ,istatus);ierr=ierr+istatus
      if(ierr/=0) then
         write(6,*) myname_, 'Could not find var(s): ',trim(failvars)
         call die (myname_, ': trouble getting pointers', ierr)
      endif

!     Calculate perturbation ps for GSI
!     ---------------------------------
      select case(which)
        case('adm')
          if ( scaleit ) then
               call ps2delp_ad_ ( Pa_per_kPa ) 
          else
               call delp2ps_    ( one ) 
          endif
        case('tlm')
          call delp2ps_    ( kPa_per_Pa ) 
        case('inc')
          do j=1,slon2
             do i=1,slat2
                xx%r2(i_p)%q(i,j) = sub_ps(i,j) * kPa_per_Pa
             enddo
          enddo
        case default
          call die ( myname_,': invalid option' )
      end select

      if (allocated(sub_ps)) then
         deallocate (sub_ps, stat=ierr )
         if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_ps)'
            return
         end if
      endif
      if ( allocated(sub_delp) ) then
          deallocate (sub_delp, stat=ierr )
          if ( ierr/=0 ) then
              stat = 99
              if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_delp)'
              return
           end if
      endif

!     Calculate all other perturbation for GSI
!     ----------------------------------------
      if (allocated(sub_tv)) then
         do k=1,nsig
            do j=1,slon2
               do i=1,slat2
                  xx%r3(i_t)%q(i,j,k) = sub_tv(i,j,k)
               enddo
            enddo
         enddo
         deallocate (sub_tv, stat=ierr )
         if ( ierr/=0 ) then
             stat = 99
             if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_tv)'
             return
         endif
      endif
      if (allocated(sub_q)) then
         do k=1,nsig
            do j=1,slon2
               do i=1,slat2
                  xx%r3(i_q)%q(i,j,k) = sub_q (i,j,k)
               enddo
            enddo
         enddo
         deallocate (sub_q, stat=ierr )
         if ( ierr/=0 ) then
             stat = 99
             if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_q)'
             return
         end if
      endif
      if (allocated(sub_u)) then
         do k=1,nsig
            do j=1,slon2
               do i=1,slat2
                  xx%r3(i_u)%q(i,j,k) = sub_u (i,j,k)
               enddo
            enddo
         enddo
         deallocate (sub_u, stat=ierr )
         if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_u)'
            return
         end if
      endif
      if (allocated(sub_v)) then
         do k=1,nsig
            do j=1,slon2
               do i=1,slat2
                  xx%r3(i_v)%q(i,j,k) = sub_v (i,j,k)
               enddo
            enddo
         enddo
         deallocate (sub_v, stat=ierr )
         if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_v)'
            return
         end if
      endif
      if ( i_oz>0 .and. allocated(sub_oz) ) then
         do k=1,nsig
            do j=1,slon2
               do i=1,slat2
                  xx%r3(i_oz)%q(i,j,k) = sub_oz (i,j,k)
               enddo
            enddo
         enddo
         if (which == 'adm') then
            if(scaleit) xx%r3(i_oz)%q  = xx%r3(i_oz)%q * PPMV2GpG
         else
            xx%r3(i_oz)%q  = xx%r3(i_oz)%q * PPMV2GpG
         endif
         deallocate (sub_oz, stat=ierr )
         if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_oz)'
            return
         end if
      endif
      if ( i_cw>0 .and. allocated(sub_cw) ) then
         do k=1,nsig
            do j=1,slon2
               do i=1,slat2
                  xx%r3(i_cw)%q(i,j,k) = sub_cw (i,j,k)
               enddo
            enddo
         enddo
         deallocate (sub_cw, stat=ierr )
         if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_cw)'
            return
         end if
      endif
      if ( allocated(sub_ql) ) then
         call gsi_bundlegetpointer(xx,'ql',  i_ql , istatus)
         if ( i_ql>0 ) then
             do k=1,nsig
                do j=1,slon2
                   do i=1,slat2
                      xx%r3(i_ql)%q(i,j,k) = sub_ql (i,j,k)
                   enddo
                enddo
             enddo
             deallocate (sub_ql, stat=ierr )
             if ( ierr/=0 ) then
                 stat = 99
                 if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_ql)'
                 return
             end if
         endif
      endif
      if ( allocated(sub_qi) ) then
         call gsi_bundlegetpointer(xx,'qi',  i_qi , istatus)
         if ( i_qi>0 ) then
            do k=1,nsig
               do j=1,slon2  
                  do i=1,slat2
                      xx%r3(i_qi)%q(i,j,k) = sub_qi (i,j,k)
                   enddo
                enddo
             enddo
             deallocate (sub_qi, stat=ierr )
             if ( ierr/=0 ) then
                 stat = 99 
                 if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_qi)'
                 return
             end if
         endif
      endif
      if ( allocated(sub_qr) ) then
         call gsi_bundlegetpointer(xx,'qr',  i_qr , istatus)
         if ( i_qr>0 ) then
            do k=1,nsig
               do j=1,slon2  
                  do i=1,slat2
                     xx%r3(i_qr)%q(i,j,k) = sub_qr (i,j,k)
                  enddo
               enddo
            enddo
            deallocate (sub_qr, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_qr)'
                return
            end if
         endif
     endif
     if ( allocated(sub_qs) ) then
        call gsi_bundlegetpointer(xx,'qs',  i_qs , istatus)
        if ( i_qs>0 ) then
           do k=1,nsig
              do j=1,slon2
                 do i=1,slat2
                    xx%r3(i_qs)%q(i,j,k) = sub_qs (i,j,k)
                 enddo
              enddo
           enddo
           deallocate (sub_qs, stat=ierr )
             if ( ierr/=0 ) then
                 stat = 99
                 if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_qs)'
                 return
             end if
           endif
     endif

      if ( which == 'tlm' ) then
          call getprs_tl   (xx%r2(i_p)%q,xx%r3(i_t)%q,xx%r3(i_prse)%q)
          call tv_to_tsen  (xx%r3(i_t)%q,xx%r3(i_q)%q,xx%r3(i_tsen)%q)
      endif

      if (allocated(sub_ts)) then
         do j=1,slon2
            do i=1,slat2
               xx%r2(i_ts)%q(i,j) = sub_ts(i,j)
            enddo
         enddo
         deallocate (sub_ts, stat=ierr )
         if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_ts)'
            return
         end if
      endif
     
      call gsi_bundlegetpointer(xx,'pblri',i_pblri,istatus)
      if ( i_pblri>0 .and. allocated(sub_pblri) ) then
         do j=1,slon2
            do i=1,slat2
               xx%r2(i_pblri)%q(i,j) = sub_pblri(i,j)
            enddo
         enddo
         deallocate (sub_pblri, stat=ierr )
         if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_pblri)'
            return
         end if
      end if

      call gsi_bundlegetpointer(xx,'pblrf',i_pblrf,istatus)
      if ( i_pblrf>0 .and. allocated(sub_pblrf) ) then
         do j=1,slon2
            do i=1,slat2
               xx%r2(i_pblrf)%q(i,j) = sub_pblrf(i,j)
            enddo
         enddo
         deallocate (sub_pblrf, stat=ierr )
         if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_pblrf)'
            return
         end if
      end if

      call gsi_bundlegetpointer(xx,'pblkh',i_pblkh,istatus)
      if ( i_pblkh>0 .and. allocated(sub_pblkh) ) then
         do j=1,slon2
            do i=1,slat2
               xx%r2(i_pblkh)%q(i,j) = sub_pblkh(i,j)
            enddo
         enddo
         deallocate (sub_pblkh, stat=ierr )
         if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_pblkh)'
            return
         end if
      end if

      CONTAINS

      subroutine delp2ps_ ( alpha )  ! delp2p
! inverse
      implicit none
      real(r_kind), intent(in) :: alpha
      real(r_kind) :: bkweight
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_aux
      integer(i_kind) i,j,k
        xx%r2(i_p)%q = zero
        do k=1,nsig
           do j=1,slon2
              do i=1,slat2
                 xx%r2(i_p)%q(i,j) = xx%r2(i_p)%q(i,j) + alpha * sub_delp(i,j,k)
              end do
           end do
        end do
        allocate ( sub_aux(slat2,slon2,nsig+1) )
           k=nsig+1
           do j=1,slon2
              do i=1,slat2
                 sub_aux(i,j,k) = zero
              end do
           end do
        do k=nsig,1,-1
           do j=1,slon2
              do i=1,slat2
                 sub_aux(i,j,k) = sub_aux(i,j,k+1) + alpha * sub_delp(i,j,k)
              end do
           end do
        end do
        do k=1,nsig+1
           do j=1,slon2
              do i=1,slat2
                 xx%r3(i_prse)%q(i,j,k) = sub_aux(i,j,k)
              end do
           end do
        end do
        deallocate ( sub_aux )
      end subroutine delp2ps_

      subroutine ps2delp_ad_ ( alpha )
! adm-only
      implicit none
      real(r_kind), intent(in) :: alpha
      real(r_kind), allocatable :: sub_aux(:,:,:) 
      integer(i_kind) i,j,k
        allocate ( sub_aux(slat2,slon2,nsig+1) )
        sub_aux = zero
        do k=1,nsig
           do j=1,slon2
              do i=1,slat2
                 sub_aux(i,j,k+1) = sub_aux(i,j,k+1) - sub_delp(i,j,k)
                 sub_aux(i,j,k)   = sub_aux(i,j,k)   + sub_delp(i,j,k)
                 sub_delp(i,j,k)  = zero
              enddo
           enddo
        enddo
        do k=1,nsig+1
           do j=1,slon2
              do i=1,slat2
                 xx%r3(i_prse)%q(i,j,k) = xx%r3(i_prse)%q(i,j,k) + alpha * sub_aux(i,j,k)
              enddo
           enddo
        enddo
!??     xx%p = zero
        deallocate ( sub_aux )
      end subroutine ps2delp_ad_

      subroutine pert2gsi2d_ ( fld, sub, stat_ )

      real(r_single), intent(in)  :: fld(:,:)
      real(r_kind),   intent(out) :: sub(:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*pert2gsi2d_'

      integer(i_kind) i,j,k,ierr
      real(r_single), allocatable :: fld3d(:,:,:)
      real(r_kind),   allocatable :: sub3d(:,:,:)

      if (gsiMode) then
         do j=1,slon2
            do i=1,slat2
               sub(i,j) = fld(j,i)
            end do
         end do
         stat_ = 0
      else
          allocate(fld3d(size(fld,1),size(fld,2),nsig))
          allocate(sub3d(size(sub,1),size(sub,2),nsig))
          fld3d=zero; sub3d=zero
          fld3d(:,:,1) = fld
          call pert2gsi_old_ ( fld3d, sub3d, ierr )
          sub(:,:) = sub3d(:,:,nsig)
          deallocate(sub3d)
          deallocate(fld3d)
          stat_ = ierr
      endif

      end subroutine pert2gsi2d_

      subroutine pert2gsi_ ( fld, sub, stat_ )

      real(r_single), intent(in)  :: fld(:,:,:)
      real(r_kind),   intent(out) :: sub(:,:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*pert2gsi_'

      integer(i_kind) i,j,k,mm1,ierr
      integer(i_kind) imr,jnp,nl

      if (gsiMode) then
          nl=nsig
          if(size(fld,3)==1) nl=1
          do k=1,nl
             do j=1,slon2
                do i=1,slat2
                   sub(i,j,k) = fld(j,i,k)
                end do
             end do
          enddo
          stat_ = 0
      else
          call pert2gsi_old_ ( fld, sub, ierr )
          stat_ = ierr
      endif

      end subroutine pert2gsi_

      subroutine pert2gsi_old_ ( fld, sub, stat_ )

      real(r_single), intent(in)  :: fld(:,:,:)
      real(r_kind),   intent(out) :: sub(:,:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*pert2gsi_old_'

      real(r_kind), allocatable :: work4d(:,:,:,:)   ! auxliar 4d array
      real(r_kind), allocatable :: work3d(:,:,:)     ! auxliar 3d array
      real(r_kind), allocatable :: work2d(:,:)       ! auxliar 2d array
      real(r_kind), allocatable :: work(:)
      integer(i_kind) i,j,k,mm1,ierr
      integer(i_kind) imr,jnp,nl

      mm1 = mype+1 
      stat_ = 0
      call mygetdims_(imr,jnp,nl)

      allocate ( work3d(nlon,nlat,nsig), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work3d)'
            return
        end if
      allocate ( work4d(imr,jnp,nl,1), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work4d)'
            return
        end if
                                                                                                                           
!     Gather GCM perturbations to root processor
!     ------------------------------------------
      call mygather_( which, imr, jnp, nl, fld, work4d(:,:,:,1) )

!     Flip horizontal and vertical
!     ----------------------------
      if ( mype==ROOT ) then
              call hflip3_ ( work4d, imr,jnp,nl )
              if (imr/=nlon .or. jnp/=nlat ) then
print*, ' WRONG PLACE TO BE ...';call flush(6)
                  if (which=='adm') then
                      work3d = zero
                      call interpack_terpv_ad ( nlon,nlat,nsig,work3d,work3d, imr,jnp,nl,work4d(:,:,:,1), ierr )
                  else if (which=='tlm') then
                      work3d = zero
                      call interpack_terpv    ( imr,jnp,nl,work4d(:,:,:,1),  nlon,nlat,nsig,work3d, ierr )
                  else
                      call die ( myname_,': invalid option' )
                  endif
              else
                  work3d(:,:,:) = work4d(:,:,:,1)
              endif
           call SwapV_ ( work3d )
      endif

!     Swap work memory
!     ----------------
      deallocate ( work4d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(work4d)'
            return
        end if
      allocate ( work(itotsub), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work)'
            return
        end if
      allocate ( work2d(lat2,lon2), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work4d)'
            return
        end if

!     Scatter to GSI subdomains
!     -------------------------
      do k=1,nsig
         if (mype==ROOT) then
             call reorder21(work3d(:,:,k),work)
         endif
         call mpi_scatterv(work,ijn_s,displs_s,mpi_rtype,&
              work2d,ijn_s(mm1),mpi_rtype,root,mpi_comm_world,ierr)
         do j=1,lon2
            do i=1,lat2
               sub(i,j,k) = work2d(i,j)
            end do
         end do
      end do

!     Release work memory
!     -------------------
      deallocate ( work2d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(work4d)'
            return
        end if
      deallocate ( work, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work)'
            return
        end if
      deallocate ( work3d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work3d)'
            return
        end if

      end subroutine pert2gsi_old_

   end subroutine pgcm2gsi1_

!------ BELOW THIS POINT: Routines of general (internal-only) use --------

!-------------------------------------------------------------------------
   subroutine SwapV_(fld)
!-------------------------------------------------------------------------
   implicit none
   real(r_kind),intent(inout) ::  fld(:,:,:)
   real(r_kind),allocatable   :: work(:,:,:)
   integer(i_kind) im, jm, km
   im   = size(fld,1)
   jm   = size(fld,2)
   km   = size(fld,3)
   allocate (work(im,jm,km))
   work = fld
   fld(:,:,km:1:-1) = work(:,:,1:km:+1)
   deallocate (work)
   end subroutine SwapV_
!-------------------------------------------------------------------------
   subroutine SwapIJK_(aij,aji)
!-------------------------------------------------------------------------
! transpose IJK-ordered array to JIK-ordered array
   implicit none
   real(r_kind),dimension(:,:,:),intent(in)    :: aij
   real(r_kind),dimension(:,:,:),intent(inout) :: aji
   integer(i_kind) :: i,k,isz,jsz,ksz,kk
!
   isz=size(aij,1)
   jsz=size(aij,2)
   ksz=size(aij,3)
   kk=1
   do k=1,ksz
      do i=1,isz
         aji(1:jsz,i,kk)=aij(i,1:jsz,k)
      end do
      kk=kk+1
   end do
   call SwapV_(aji)
   end subroutine SwapIJK_

   subroutine hflip3_ ( q,im,jm,km )
      implicit none
      integer(i_kind)  im,jm,km,i,j,k
      real(r_kind), intent(inout) :: q(im,jm,km)
      real(r_kind), allocatable   :: dum(:)
      allocate ( dum(im) )
      do k=1,km
      do j=1,jm
      do i=1,im/2
         dum(i) = q(i+im/2,j,k)
         dum(i+im/2) = q(i,j,k)
      enddo
         q(:,j,k) = dum(:)
      enddo
      enddo
      deallocate ( dum )
   end subroutine hflip3_
  
   subroutine hflip2_ ( q,im,jm )
      implicit none
      integer(i_kind)  im,jm,i,j
      real(r_kind), intent(inout) :: q(im,jm)
      real(r_kind), allocatable   :: dum(:)
      allocate ( dum(im) )
      do j=1,jm
      do i=1,im/2
         dum(i) = q(i+im/2,j)
         dum(i+im/2) = q(i,j)
      enddo
         q(:,j) = dum(:)
      enddo
      deallocate ( dum )
   end subroutine hflip2_

   subroutine mygetdims_ ( im_world, jm_world, km_world )
   use GSI_GridCompMod, only: GSI_AgcmPertGrid
   implicit none
   integer(i_kind), intent(out) :: im_world
   integer(i_kind), intent(out) :: jm_world
   integer(i_kind), intent(out) :: km_world
  
   character(len=*),parameter :: Iam=myname//'mygetdims_'
   integer  DIMS(ESMF_MAXGRIDDIM)
   integer  rc,status

   call MAPL_GridGet(GSI_AgcmPertGrid, globalCellCountPerDim=DIMS, RC=STATUS)
        VERIFY_(STATUS)

   im_world = DIMS(1)
   jm_world = DIMS(2)
   km_world = DIMS(3)
   end subroutine mygetdims_

   subroutine myscatter_ ( what, im_world, jm_world, km_world, vari, varo )
!  NOTE: there is not need for Halo here ... this is going onto a "model" field
   use GSI_GridCompMod, only: GSI_AgcmPertGrid
   implicit none
   character(len=*),intent(in)    :: what
   integer,         intent(in)    :: im_world,jm_world,km_world
   real(r_kind),    intent(in)    :: vari(im_world,jm_world,km_world)
   real(4),         intent(inout) :: varo(:,:,:)

   character(len=*),parameter :: Iam=myname//'myscatter_'
   real(4),allocatable :: work2d(:,:)
   integer  k, rc,status

    allocate(work2d(im_world,jm_world))
    if (what=='adm') then
        do k=1,km_world  ! this would be a good place to swapV
           work2d=vari(:,:,k)
           call ArrayScatter(varo(:,:,k), work2d, GSI_AgcmPertGrid, rc=status)   ! _RT: ideally this should be the AD of gather
                VERIFY_(STATUS)
        enddo
    else
        do k=1,km_world  ! this would be a good place to swapV
           work2d=vari(:,:,k)
           call ArrayScatter(varo(:,:,k), work2d, GSI_AgcmPertGrid, rc=status)
                VERIFY_(STATUS)
        enddo
    endif
    deallocate(work2d)
    
   end subroutine myscatter_

   subroutine mygather_ ( what, im_world, jm_world, km_world, vari, varo )
!  NOTE: there is not need for Halo here ... this is coming from a "model" field
   use GSI_GridCompMod, only: GSI_AgcmPertGrid
   implicit none
   character(len=*),        intent(in)    :: what
   integer,                 intent(in)    :: im_world,jm_world,km_world
   real(4),                 intent(in)    :: vari(:,:,:)
   real(r_kind),            intent(inout) :: varo(im_world,jm_world,km_world)

   character(len=*),parameter :: Iam=myname//'mygather_'
   real(4),allocatable :: work2d(:,:)
   integer  k, rc,status

    allocate(work2d(im_world,jm_world))
    if (what=='adm') then
        do k=1,km_world  ! this would be a good place to swapV
           work2d=zero
           call ArrayGather(vari(:,:,k), work2d, GSI_AgcmPertGrid, rc=status)   ! _RT: ideally this should be the AD of gather
                VERIFY_(STATUS)
           varo(:,:,k) = work2d
        enddo
    else
        do k=1,km_world  ! this would be a good place to swapV
           work2d=zero
           call ArrayGather(vari(:,:,k), work2d, GSI_AgcmPertGrid, rc=status)
                VERIFY_(STATUS)
           varo(:,:,k) = work2d
        enddo
    endif
    deallocate(work2d)
    
   end subroutine mygather_

end module geos_pertStateIO
