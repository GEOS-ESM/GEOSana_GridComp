!-----------------------------------------------------------------------
module fvcore
   use precision
#if defined( SPMD )
      use mod_comm,  only : gid
#define CPP_PRT_PREFIX  if(gid == 0)
#else
#define CPP_PRT_PREFIX
#endif
   implicit none

! Geometric arrays
   real(r8), allocatable :: sine(:),cosp(:),sinp(:),cose(:)
   real(r8), allocatable :: acosp(:)
   real(r8), allocatable :: sinlon(:),coslon(:),cosl5(:),sinl5(:)

   double precision pi, zamda, zam5
   real(r8) acap, rcap
   integer ns

   integer :: fvcore_tape_rec                   ! FastOpt
   integer :: dynpkg_n2                         ! FastOpt
   integer :: dynpkg_nsplit                     ! FastOpt
   integer :: nsplit, n2                        ! FastOpt
   real (r8) dt
   integer :: ng_c
   integer nx          ! # of split pieces in x-direction

contains

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  fvcore_initialize --- Initialization
! !Author: Ralf Giering, FastOpt
!
! !INTERFACE:
subroutine fvcore_initialize( im, jm, km, jfirst, jlast, ng_d, ng_s,  &
			      ns0, jord, pdt )

! !USES:
   use precision
   use cd_core, only: cd_core_allocate
   use benergy, only: benergy_initialize
   use mapz_module,  only: te_map_initialize

#if defined ( SPMD )
      use mod_comm, only: numcpu
#endif

   integer :: im, jm, km
   integer :: ns0
   integer :: jfirst, jlast
   integer :: jord, ng_d, ng_s
   integer :: pdt                            ! FastOpt

   integer :: jcd
   integer :: imh
   real (r8) dp, dl
   integer :: i, j

   real (r8) bdt

      if( jord <= 2 ) then
         jcd =  1
      else
         jcd =  -2
      endif

#if defined( SPMD )
        ng_c = min(abs(jcd ),2)
        if ( (jlast-jfirst+1)/numcpu >= 4 ) then
            nx = 1
        else
            nx = 4
        endif
#else
        ng_c = 0
        nx = 1
#endif

        allocate( sine(jm),cosp(jm),sinp(jm),cose(jm), acosp(jm) )
        allocate( sinlon(im), coslon(im), cosl5(im), sinl5(im) )

        pi = 4.d0 * datan(1.d0)

        call setrig(im,jm,dp,dl,cosp,cose,sinp,sine)

        acap = im*(1.+sine(2)) / dp
        rcap = 1.d0 / acap

        imh = im/2
        if(im .ne. 2*imh) then
           write(6,*) 'im must be an even integer'
           stop
        endif

! Define logitude at the center of the volume
! i=1, Zamda = -pi

        do i=1,imh
           zam5          = (dble(i-1)-0.5d0) * dl
           cosl5(i)      = -dcos(zam5)
           cosl5(i+imh)  = -cosl5(i)
           sinl5(i)      = -dsin(zam5)
           sinl5(i+imh)  = -sinl5(i)
           zamda         = dble(i-1)*dl
           coslon(i)     = -dcos(zamda)
           coslon(i+imh) = -coslon(i)
           sinlon(i)     = -dsin(zamda)
           sinlon(i+imh) = -sinlon(i)
        enddo

        do j=2,jm-1
          acosp(j) = 1.d0 / cosp(j)
        enddo
          acosp( 1) = rcap * im
          acosp(jm) = rcap * im

        ns = ns0

   call cd_core_allocate( im, jm, km, jfirst, jlast, ng_c, ng_d, ng_s )

   call benergy_initialize( im, jm )

   call te_map_initialize( im, jm, km )

! Determine splitting
   bdt = pdt
   call d_split( ns,   im,   jm,   bdt )

! Second level splitting
   n2 = max ( 1, ns/4 )
   nsplit = (ns+n2-1) / n2

   dt = bdt / float(nsplit*n2)

   dynpkg_n2     = n2
   dynpkg_nsplit = nsplit

   CPP_PRT_PREFIX print *,' dynpkg_n2     = ', dynpkg_n2
   CPP_PRT_PREFIX print *,' dynpkg_nsplit = ', dynpkg_nsplit

!EOC
end subroutine fvcore_initialize
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  fvcore_finalize --- Finalization
! !Author: Ralf Giering, FastOpt
!
! !INTERFACE:
subroutine fvcore_finalize

! !USES:
   use precision
   use cd_core, only: cd_core_deallocate
   use benergy, only: benergy_finalize
   use mapz_module, only: te_map_finalize

   implicit none

   deallocate( sine, cosp, sinp,cose, acosp )
   deallocate( sinlon, coslon, cosl5, sinl5 )

   call cd_core_deallocate

   call benergy_finalize

   call te_map_finalize

!EOC
end subroutine fvcore_finalize
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: fvcore_do --- driver for the finite-volume dynamical core
!           MPI in latitudinal direction SMP in (mostly) vertical direction
!
! FasOpt added argument nc to avoid assumed size array q
! !INTERFACE:

! Need to re-dimension u,v,pt

      subroutine fvcore_do( im, jm, km, nc, jfirst, jlast,                   &
                         ng_d, ng_s, nq, ps, pe,                             &
                         delp, u, v, pt, q, pk, pkz, phis, ns0, ndt,         &
                         ptop, om, cp, rg, ae, iord, jord, kord,             &
                         umax, omga, peln, consv, convt )
 
! !USES:
      use mapz_module, only: te_map                        ! FastOpt
      use cd_core, only: cd_core_initialize, cd_core_do    ! FastOpt
      use benergy, only: benergy_do                        ! FastOpt

      use cd_core    , only: cd_core_tape_rec              ! FastOpt


      use precision
      use timingModule
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in):: im        ! dimension in east-west
      integer, intent(in):: jm        ! dimension in North-South
      integer, intent(in):: km        ! number of Lagrangian layers
      integer, intent(in):: nc        ! number of tracer FastOpt
      integer, intent(in):: jfirst    ! starting latitude index for MPI
      integer, intent(in):: jlast     ! ending latitude index for MPI
      integer, intent(in):: ng_d
      integer, intent(in):: ng_s
      integer, intent(in):: nq        ! total # of tracers to be advected
      integer, intent(in):: ndt       ! the large time step in seconds
                                      ! Also the mapping time step in this setup
      integer, intent(in):: ns0       ! total number of splits for Lagrangian
                                      ! dynamics; a suitable value will be automatically
                                      ! determined if ns0 = 0 (i.e., if you don't know what value
                                      ! to be used try ns0=0
      real (r8), intent(in):: ptop    ! constant pressure at model top (pascal)
      real (r8), intent(in):: om      ! angular velocity of earth's rotation  
      real (r8), intent(in):: cp      ! heat capacity of air at constant pressure
      real (r8), intent(in):: rg      ! gas constant of the air
      real (r8), intent(in):: ae      ! radius of the earth (m)

      integer, intent(in):: iord      ! parameter controlling monotonicity in E-W
                                      ! recommendation: iord=4
      integer, intent(in):: jord      ! parameter controlling monotonicity in N-S 
                                      ! recommendation: jord=4
      integer, intent(in):: kord      ! parameter controlling monotonicity in mapping
                                      ! recommendation: kord=4
      real (r8), intent(in):: umax    ! upper bound of the maximum u-wind (m/s)
                                      ! (eg, 300 m/s is a good value)
      logical, intent(in):: consv     ! conserve the globally integrated Total Energy
      logical, intent(in):: convt     ! flag to control pt output (see below)
                                      ! Set convt =.T. to output (virtual) temperature
      real (r8), intent(in):: phis(im,jfirst:jlast)   ! surface geopotential (grav*zs)

! !INPUT/OUTPUT
      real (r8), intent(inout):: ps(im,jfirst:jlast)      ! surface pressure (pascal)
      real (r8), intent(inout):: delp(im,jfirst:jlast,km) ! pressure thickness (pascal)

! Ghosted prog arrays:
      real (r8), intent(inout):: u(im,jfirst-ng_d:jlast+ng_s,km)  ! u-wind (m/s)
      real (r8), intent(inout):: v(im,jfirst-ng_s:jlast+ng_d,km)  ! v-wind (m/s)
      real (r8), intent(inout):: pt(im,jfirst-ng_d:jlast+ng_d,km)
                                      ! Input: scaled (virtual) potential
                                      ! temperature  defined as
                                      ! Virtual_temperature (K)/pkz
                                      ! Output: (virtual) temperature (deg K) if convt is true
                                      ! Output is (virtual) potential temperature
                                      ! if convt is false
      real (r8), intent(inout):: q(im,jfirst-ng_d:jlast+ng_d,km,nc)
                                      ! tracers (e.g., specific humidity)
                                      ! tracer mass / moist_air_mass

! The following three arrays must be pre-computed as input to benergy(). They are NOT
! needed if consv=.F.; updated on output (to be used by physdrv)
! Please refer to routine pkez on the algorithm for computing pkz
! from pe and pk

      real (r8), intent(inout):: pe(im,km+1,jfirst:jlast)  ! pressure (pascal) at layer edges
      real (r8), intent(inout):: pk(im,jfirst:jlast,km+1)  ! pe**cappa (cappa = rg/cp)
      real (r8), intent(inout):: pkz(im,jfirst:jlast,km)   ! finite-volume mean of pk

! !OUTPUT (input values are not used):
      real (r8), intent(out):: peln(im,km+1,jfirst:jlast)  ! log pressure (pe) at layer edges
                                                           ! outputed values can be used by physdrv
      real (r8), intent(out):: omga(im,km,jfirst:jlast)    ! vertical pressure velocity (pa/sec)
                                                           ! This is the rate of change of the Lagrangian
                                                           ! interface in pascal/sec unit.

! !DESCRIPTION:
!
! Developer: Shian-Jiann Lin, NASA/GSFC; email: lin@dao.gsfc.nasa.gov
!
! D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!               u(i,j+1)
!                 |
!      v(i,j)---delp(i,j)---v(i+1,j)
!                 |
!               u(i,j)

! External routine required: the user needs to supply a subroutine to set up
!                            "Eulerian vertical coordinate" for remapping purpose.
!                             Currently this routine is named as set_eta()
!                             In principle any terrian following vertical
!                             coordinate can be used. The input to fvcore
!                             need not be on the same vertical coordinate
!                             as the output.
!
! Remarks: values at poles for both u and v need not be defined; but values for
!          all other scalars needed to be defined at both poles (as polar cap mean
!          quantities). Tracer advection is done "off-line" using the
!          large time step. Consistency is maintained by using the time accumulated
!          Courant numbers and horizontal mass fluxes for the FFSL algorithm.
!          The input "pt" can be either dry potential temperature
!          defined as T/pkz (adiabatic case) or virtual potential temperature
!          defined as T*/pkz (full phys case). IF convt is true, pt becomes
!          (virtual) temperature for input to physics driver. Therefore, 
!          pt must be converted back to (virtual) potential temperature after
!          the physics routines. pt will remain as potential temperature
!          if convt is false.
!          The user may set the value of nx to optimize the SMP performance
!          The optimal valuse of nx depends on the total number of available
!          shared memory CPUs per node.
!
! !REVISION HISTORY:
!   SJL 99.04.13:  Initial SMP version delivered to Will Sawyer
!   WS  99.10.03:  1D MPI completed and tested; 
!   WS  99.10.11:  Additional documentation
!   WS  99.10.28:  benergy and te_map added; arrays pruned
!   SJL 00.01.01:  SMP and MPI enhancements; documentation
!   WS  00.07.13:  Changed PILGRIM API
!
!EOP
!-----------------------------------------------------------------------
!BOC

! Local variables:

      integer ipe, it, k
      integer n

      real (r8)   te0
      real (r8)   cappa
      integer     te_map_tape_rec           ! FastOpt

! Local arrays
      real (r8), allocatable :: worka(:,:,:),dp0(:,:,:),cx(:,:,:),cy(:,:,:)
      real (r8), allocatable :: mfx(:,:,:), mfy(:,:,:)
      real (r8), allocatable :: delpf(:,:,:), uc(:,:,:), vc(:,:,:)
      real (r8), allocatable :: dwz(:,:,:), pkc(:,:,:), wz(:,:,:)
      real (r8), allocatable :: dpt(:,:,:)

      double precision pi, zamda, zam5
      logical filter, fill

      integer i, j, icd, jcd

      data filter /.true./              ! polar filter
      data fill   /.true./              ! perform a simple filling algorithm
                                        ! in case negatives were found

      integer cd_tape_rec_n             ! FastOpt
      integer irec                      ! FastOpt

      cappa = rg / cp

      if( iord <= 2 ) then
         icd =  1
      else
         icd = -2
      endif

      if( jord <= 2 ) then
         jcd =  1
      else
         jcd =  -2
      endif

! Allocate temporary work arrays
! Change later to use pointers for SMP performance???
! (prime candidates: uc, vc, delpf)

      allocate( worka(im,jfirst:     jlast,     km) )
      allocate(   dp0(im,jfirst:     jlast,     km) )
      allocate(   mfx(im,jfirst:     jlast,     km) )
      allocate(   mfy(im,jfirst:     jlast+1,   km) )
      allocate(    cx(im,jfirst-ng_d:jlast+ng_d,km) )
      allocate(    cy(im,jfirst:     jlast+1,   km) )
      allocate( delpf(im,jfirst-ng_d:jlast+ng_d,km) )
      allocate(    uc(im,jfirst-ng_d:jlast+ng_d,km) )
      allocate(    vc(im,jfirst-2:   jlast+2,   km) )
      allocate(   dpt(im,jfirst-1:   jlast+1,   km) )
      allocate(   dwz(im,jfirst-1:    jlast,    km+1) )
      allocate(   pkc(im,jfirst-1:   jlast+1,   km+1) ) 
      allocate(    wz(im,jfirst-1:   jlast+1,   km+1) )

      delpf = 0.                   ! FastOpt: init to be stored
      te0   = 0.                   ! FastOpt: init to be stored

#ifdef  TAPING_IN_DYNPKG
! FastOpt: start taping here
! variables n2,nsplit,ng_d need to be defined priviously

!$taf init dynpkg_1     = static, 1
!$taf init dynpkg_n2    = static, n2
!$taf init dynpkg_tape  = static, nsplit*n2
!$taf init cd_core_tape = static, nsplit*n2
!$taf init te_map_tape  = static, 1
  fvcore_tape_rec = 0
#endif

#ifdef  D_SW_TAPING_FVCORE
!$taf init d_sw_tape1   = static, nsplit*n2 * km*1
!$taf init d_sw_tapej   = static, nsplit*n2 * km*(jm-1)
#endif

#ifdef  SW_CORE_STORING_FVCORE
!$taf init c_sw_tape1   = static, nsplit*n2 * km*1
!$taf init c_sw_tape2   = static, nsplit*n2 * km*(jm-1)
!$taf init tpcc1_tape   = static, nsplit*n2 * km*(jm-1)
!$taf init tpcc2_tape   = static, nsplit*n2 * km*(jm-1)
#endif

! First touch pkc and wz??? (bufferpack is multitask in vertical but geop
! computations are parallel in j-loop)

      if ( km .gt. 1 ) then         ! not shallow water equations

      if( consv ) then
! Compute globally integrated Total Energy (te0)
      call timing_on('BENERGY')
      call benergy_do(im, jm, km, u, v, pt, delp, pe, pk, pkz, phis,  &
                   ng_d, ng_s, cp, te0, mfx, dp0, jfirst, jlast)
      call timing_off('BENERGY') 
      endif

      endif

!$taf store te0            = dynpkg_1, rec=fvcore_tape_rec+1

      do 2000 n=1, n2

      cd_tape_rec_n = n-1 + fvcore_tape_rec*n2

      if( nq .gt. 0 ) then

!$omp parallel do private(i,j,k)
      do 1000 k=1,km
         do j=jfirst,jlast
            do i=1,im
! Save initial delp field before the small-time-step
! Initialize the CFL number accumulators: (cx, cy)
! Initialize total mass fluxes: (mfx, mfy)
!
               dp0(i,j,k) = delp(i,j,k)
                cx(i,j,k) = 0.
                cy(i,j,k) = 0.
               mfx(i,j,k) = 0.
               mfy(i,j,k) = 0.
            enddo
         enddo
1000  continue

      endif

      call cd_core_initialize( im, jm, km, jfirst, jlast,     &
                               ng_c, ng_d, ng_s,              &
                               dt, ae, om,  ptop, umax,       &
                               sinp, cosp, cose, acosp, cappa &
                             )
      do 1500 it=1, nsplit

         irec = it + cd_tape_rec_n*nsplit

!$taf store delp,delpf,pt            = dynpkg_tape, rec=irec

      if(it == nsplit .and. n == n2) then
         ipe = 1                     ! end of fvcore; output pe for te_map
      elseif(it == 1 .and. n == 1) then
         ipe = -1                    ! start of cd_core
      else
         ipe = 0
      endif

! Call the Lagrangian dynamical core using small tme step

      cd_core_tape_rec = it-1 + cd_tape_rec_n*nsplit

      call timing_on('CD_CORE')
      call cd_core_do(im,  jm,  km,  nq, nx,                   &
                      jfirst, jlast,   u,  v,  pt,             &
                      delp,     pe,    pk, ns,                 &
                      dt, ptop , umax  , fill, filter, acap,   &
                      ae, rcap, cp, cappa,                     &
                      icd, jcd, iord, jord, ng_c, ng_d, ng_s,  &
                      ipe, om, phis, sinp, cosp, cose, acosp,  &
                      sinlon, coslon, cosl5, sinl5,            &
                      cx  , cy , mfx, mfy,                     &
                      delpf, uc, vc, pkz, dpt,                 &
                      worka, dwz, pkc, wz )
      call timing_off('CD_CORE')

1500  continue

      if(nq .ne. 0) then

! Perform large-tme-step scalar transport using the accumulated CFL and
! mass fluxes

!$taf store dp0,q,cx,cy,mfx,mfy = dynpkg_n2, rec=cd_tape_rec_n+1

      call timing_on('TRAC2D')
      call trac2d(dp0, q, nq, cx, cy, mfx, mfy, iord, jord,    &
                  ng_d, sine, cosp, acosp, acap, rcap, fill,   &
                  im, jm, km, jfirst, jlast, pkz, worka,       &
		  cd_tape_rec_n                                &
		  )
      call timing_off('TRAC2D')

      endif

2000  continue

      if ( km .gt. 1 ) then           ! not shallow water equations

! Peroform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.

!$taf store pk,q,pe,ps,pt = dynpkg_1, rec=fvcore_tape_rec+1

      te_map_tape_rec = fvcore_tape_rec

      call timing_on('TE_MAP')
      call te_map(consv, convt, ps, omga,  pe, delp, pkz, pk, ndt,   &
                  im, jm, km, nx, jfirst, jlast, nq,  u,  v,  pt,    &
                  q, phis,  cp, cappa, kord, peln, te0,              &
                  ng_d, ng_s, te_map_tape_rec )
      call timing_off('TE_MAP')
      endif

      deallocate( mfy )
      deallocate( mfx )
      deallocate(  cy )
      deallocate(  cx )
      deallocate( dp0 )
      deallocate( delpf )
      deallocate( uc    )
      deallocate( vc    )
      deallocate( dpt   )
      deallocate( pkc   )
      deallocate( dwz   )
      deallocate(  wz   )
      deallocate( worka )

      return
!EOC
      end subroutine fvcore_do
!-----------------------------------------------------------------------


end module fvcore
!-----------------------------------------------------------------------
