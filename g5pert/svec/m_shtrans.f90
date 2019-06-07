!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: m_shtrans_DG --- Module for spectral tansforms on G or D-grid
!                     
!
! !INTERFACE:
!    
     module m_shtrans_DG

!USES:

      use precision
      use m_eigen
      use m_const, only : earthr=>radius_earth

      implicit none

! !PUBLIC MEMBER FUNCTIONS:
 
      Private
      Public  :: sh_initial_D   ! routine to set some needed D-grid variables
      Public  :: sh_initial_G   ! routine to set some needed G-grid variables
      Public  :: sh_transforms  ! routine for forward or inverse transforms
      Public  :: sh_calc_power  ! routine to compute power spectra from scoefs
      Public  :: sh_test_coefs  ! routine for testing module
      Public  :: sh_clean       ! deallocate module-resident arrays
      Public  :: sh_prt_allc    ! for testing: print all nonzero coefs
      Public  :: sh_intsqr_grid ! for testing: compute global integral of f**2

! !DESCRIPTION: 
!
! This module is designed for performing exact spectral transforms on 
! either a Gaussian grid or the  
! GMAO FVGCM D-grid.  This means, if one starts with spherical harmonic 
! spectral coefficients for T, stream function and velocity potential, 
! determines global fields of T, u, and v from them, and then attempts to 
! re-compute the original spectral coefficients from those fields, the 
! initial and final set of spectral coefficients are identical. 
!
! See the routines sh_initial_all for information on how to specify the 
!   spectral truncation and to set up required initial arrays.
! See the routines sh_initial_D or sh_initial_G for information on how 
!   to call routines to initialize some variables and arrays for 
!   performing spectral transforms on either the D or G grids, respectively.
! See the routine sh_transforms for information on how to perform forward
!   or inverse transforms on a scalar or wind-component field.
! See the routine sh_calc_power for information on how to compute power 
!   spectra from the spherical harmonic coefficients for wind or scalar fields.
! The routine sh_clean should be called when the use of this module is 
!   concluded, so that module-resident arrays are deallocated.
! See the routine sh_test_coefs for an example of use of this module.
!
! Also, see routine sh_m1 for description of one peculiar aspect of 
! this computation. 
!
! Most of the code has been written rather generically so that it could
! be more easily modified for some lat/lon grid other than the D-grid.
! In this case, the routines sh_fm2s, sh_initial_D, and sh_transforms
! would all require modification.
!
! The code requires an FFT for forward and inverse zonal transforms.
! See routine sh_fft.  
!
! The computation here first computes intermediate results in the 
! form of spherical harmonic coefficients of u*cos and v*cos before
! computing coefficients for stream function and velocity potential.
! See Temperton, C, May 1991 Mon. Wea. Rev., vol. 119, pages 1303-1307.
!
! For those aspects of the transformations that concern the grid being 
! not a Gaussian one, see Swarztrauber, P.N. and W.F. Spotz, 2000,
! J. Comp. Physics, vol. 159, 213-230.
!
! For description of the general problem of spherical harmonic 
! transforms see the chapter by B. Machenhauer in the WMO GARP publication 
! No 17, Vol. 2, 1979: Numerical Methods used in Atmospheric Models.
!

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!
!EOP
!-------------------------------------------------------------------------
      character(len=*), parameter :: myname = 'm_shtrans'

!
!
!  The following arrays are saved because they are are created by 
!  calling sh_initial_D (or G) and the routines they call, but later 
!  required by sh_transforms and the routines it calls.  
!  Others of these variables are created in sh_transforms to be used
!  in lower routines. 
!   
      integer,  allocatable,save :: jfirstX(:)  ! first relevant latitude
      integer,  allocatable,save :: jlastX(:)   ! last relevant latitude
      integer,  allocatable,save :: nalpX(:)    ! ordering of alpX
      integer,  allocatable,save :: n4alpX(:)   ! picks sub-grids
      real(r8), allocatable,save :: cospX(:,:)  ! cosines of lats
      real(r8), allocatable,save :: sinpX(:,:)  ! sines of lats
      real(r8), allocatable,save :: abfaclatsX(:,:,:) ! factors for alp
      real(r8), allocatable,save :: epsiX(:,:)        ! factors for alp
      real(r8), allocatable,save :: wmatrixX(:,:,:,:) ! integrating weights
      real(r8), allocatable,save :: shiftX(:)         ! shift Prime Meridian
      real(r8), allocatable,save :: tableX(:,:,:)     ! trig funcs of longit.
      real(r8), allocatable,save :: table_fftX(:)     ! table for FFT
      real(r8), allocatable,save :: gweightsX(:)      ! Gaussian weights
      real(r8), allocatable,save :: alpX(:,:,:)       ! assoc. Legendr polyn.
      integer,  save             :: igridsX           ! number of sub-grids
      integer,  save             :: n3alpX            ! number of alp grids
      integer,  save             :: jalpX             ! number of alp lats
      integer,  save             :: igalpX            ! current alp grid
      real(r8), save             :: dpX               ! spacing of lats
      logical,  save             :: symmetricX        ! equatorial sym lats?
      logical,  save             :: alpsaveX          ! save alp values?
      character(len=1), save     :: grid_nameX        ! D or G grid  
!
      real(r8), parameter  :: zero=0.d0
      real(r8), parameter  :: one=1.d0
      real(r8), parameter  :: two=2.d0
      real(r8), parameter  :: three=3.d0
!
      CONTAINS
!
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_abfaclats --- Compute some arrays for Legendre polynomials
!                     
!
! !INTERFACE:
!    
      subroutine sh_abfaclats (abfaclats,sinlat,mmax,jmax,jfirst,jlast)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)   :: mmax
      integer,  intent(in)   :: jmax
      integer,  intent(in)   :: jfirst
      integer,  intent(in)   :: jlast
      real(r8), intent(in)   :: sinlat(jmax)

! !OUTPUT PARAMETERS:

      real(r8), intent(out)  :: abfaclats(1:jmax,0:mmax)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 

! Abfaclats is used to compute the Legendre polynomials.  Here the value 
! is such that the integral of the square of the polynomials from pole
! to pole is 1, given Gaussian weights that sum to 1. That is, the
! normalization is consistent with a weighting by fractional area on the 
! sphere.

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer                :: j
      integer                :: js
      integer                :: m
      real(r8)               :: cos2
      real(r8)               :: a
      real(r8)               :: b
      real(r8)               :: prod
      real(r8)               :: pi4
      real(r8)               :: approx1
!
! Value of approx1 is the sine of a latitude closer to the pole than any
! grid latitude, here taken as pi/(4*js) radians from the pole.  Use the 
! relationship sin(pi/2 -t)=cos(t)~1-0.5t**2 for small values of t.
! This is used to recognize the pole point if its sine is not precisely one
! due to the precision of calculations.  js=jmax if grid is symmetric.
      js=(jfirst-1)*2+jlast
      pi4=datan(one)
      approx1=dsqrt(two)*(one-((pi4/js)**2)/two)   
!
      abfaclats(jfirst:jlast,0)=dsqrt(one)
      do j=jfirst,jlast
        if ( abs(sinlat(j)) > approx1 ) then
          abfaclats(j,1:mmax)=zero      ! value for pole point
        else                            ! not pole point
          cos2=one-sinlat(j)*sinlat(j)
          prod=one
          a=one
          b=zero
          do m=1,mmax      
            a=a+two
            b=b+two
            prod=prod*cos2*a/b
            abfaclats(j,m)=dsqrt(prod)  
          enddo  ! loop over m
        endif    ! test on pole point
      enddo      ! loop over j   
!
      end subroutine sh_abfaclats
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_alp_all --- Pre-Compute all assoc. Legendre polynomials 
!                     
!
! !INTERFACE:
!    
      subroutine sh_alp_all (mmax,jmax,ntrunc)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: mmax,jmax
      integer,  intent(in)  :: ntrunc(0:mmax)
!                              
! n3alpX     is a saved global integer indicating the number of 
!            independent grids on which alpX need be defined.
! symmetricX is a logical variable that indicates if the grids are
!            equatorially symmetric

! !OUTPUT PARAMETERS:
!
!  This routine is saving globally-defined values of 
!     jalpX = number of lats for wich alp must be computed, allowing
!             for reduction from jmax if grid is equatorially symmetric
!     nalpX = description of ordering of alp as a function of n,m, 
!     alpX  = save associated Legendre polynomials.

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!
! This routine computes all the associated Legendre polynomials
! required for all latititudes for 1 grid. While this pre-computation
! require more storage, it greatly reduces some otherwise 
! redundant calculations that scale as N**3. 
!

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer  :: ig   ! grid index for looping over required sub-grids
      integer  :: j    ! latitude index
      integer  :: m    ! zonal wavenumber
      integer  :: nm   ! n-m of associated legendre polynomial
      integer  :: np   ! counter for values of allowed m,n index pairs
      integer  :: ntpx ! maximum n-m for this m
      integer  :: ig1  ! corresponding index of desired sub-grid
      real(r8), allocatable :: alp(:)
!
      if (symmetricX) then
        jalpX=1+jmax/2
      else
        jalpX=jmax
      endif
!
      allocate (nalpX(0:mmax+1))
!
      np=1
      do m=0,mmax
        nalpX(m)=np
        np=np+ntrunc(m)+2
      enddo
      nalpX(mmax+1)=np
!
      allocate (alpX(np,jalpX,n3alpX))      
      allocate (alp(0:ntrunc(0)+1))  ! ntrunc(0) is maximum n-m 
!
      do ig=1,n3alpX
        ig1=n4alpX(ig)
        do m=0,mmax
          ntpx=ntrunc(m)+1        
          do j=jfirstX(ig1),min(jlastX(ig1),jalpX)
            call sh_alp_1 (alp(:),ntpx,m,j,sinpX(j,ig1),epsiX(:,m), &
                           abfaclatsX(j,m,ig1),.false.)
            do nm=0,ntpx
              np=nm+nalpX(m)
              alpX(np,j,ig)=alp(nm)
            enddo

          enddo
        enddo
      enddo 
      deallocate (alp)
!
      end subroutine sh_alp_all
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_alp_1 --- Compute associated Legendre polynomials for 1 m,j
!                     
!
! !INTERFACE:
!    
      subroutine sh_alp_1 (alp,nlast,m,j,sinlat,epsi,abfac,alpsave)

!USES:

      implicit none

! !INPUT PARAMETERS:

      logical,  intent(in)   :: alpsave
      integer,  intent(in)   :: nlast      
      integer,  intent(in)   :: m
      integer,  intent(in)   :: j
      real(r8), intent(in)   :: sinlat
      real(r8), intent(in)   :: epsi(0:nlast)
      real(r8), intent(in)   :: abfac

! !OUTPUT PARAMETERS:

      real(r8), intent(out)  :: alp(0:nlast)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  Compute associated Legendre polynomials for given m and latitude.
!  The index of alp is the value n-m. The polynomials are normalized such
!  that are normalized such that the integral of their squared value from 
!  pole to pole equals 1. This renormalization is consistent with a 
!  weighting that represents a fraction of the surface area on the sphere. 
!  This normalization is specified by the value of abfac.

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer                :: nm  ! = n-m
      integer                :: np  ! indicator for order of alp storage 
      real(r8)               :: fm  ! = real(m)  
!
      if (alpsave) then  ! simply copy required alp from saved values
! 
        do nm=0,nlast
          np=nalpX(m)+nm
          alp(nm)=alpX(np,j,igalpX)
        enddo
!      
      else  ! otherwise actually compute alp
!
        fm=real(m, r8)
!
!  compute values for n-m=0 and n-m=1
        alp(0)=abfac
        if (nlast>0) then
          alp(1)=dsqrt(two*fm+three)*sinlat*alp(0)
        endif
!
! compute values for n-m>1 using recursion
        if (nlast>1) then
          do nm=2,nlast
            alp(nm)=(sinlat*alp(nm-1)-epsi(nm-1)*alp(nm-2))/epsi(nm)
          enddo
        endif
!
      endif
!
      end subroutine sh_alp_1
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_bssl --- Has zeros of Bessle functions 
!                     
!
! !INTERFACE:
!    
      subroutine sh_bssl (bes,n)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)   :: n

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: bes(n)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: These zeros of Bessel functions are used to compute
!               the Gaussian latitudes 

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      real(r8)              :: bz(50)
      integer               :: nn,j
      real(r8)              :: pi
      data pi/3.14159265358979/
      data bz / 2.4048255577, 5.5200781103, 8.6537279129, 11.7915344391, &
   14.9309177086, 18.0710639679, 21.2116366299, 24.3524715308, &
   27.4934791320, 30.6346064684, 33.7758202136, 36.9170983537, 40.0584257646, &
   43.1997917132, 46.3411883717, 49.4826098974, 52.6240518411, &
   55.7655107550, 58.9069839261, 62.0484691902, 65.1899648002, &
   68.3314693299, 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711, &
   84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819, & 
   96.6052679510, 99.7468198587, 102.8883742542, 106.0299309165, & 
   109.1714896498, 112.3130502805, 115.4546126537, 118.5961766309, &
   121.7377420880, 124.8793089132, 128.0208770059, 131.1624462752, &
   134.3040166383, 137.4455880203, 140.5871603528, 143.7287335737, &
   146.8703076258, 150.0118824570, 153.1534580192, 156.2950342685/

      nn=min0(n,50)
      do 1 j=1,nn
      bes(j)=bz(j)
    1 continue
      if (n.le.50) return
      do 2 j=51,n
      bes(j)=bes(j-1)+pi
    2 continue
!
      end subroutine sh_bssl 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_calc_power --- Routine for calc of power spectra from scoefs
!                     
!
! !INTERFACE:
!    
      subroutine sh_calc_power (mmax,nmax,kmax,imax,nspects,nlevels, &
                                nfields,mtrunc,index,power,scoefs,name)

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=nfields), intent(in)  :: name ! values of S or V 
      integer,  intent(in)  :: mmax    ! maximum zonal wave number m
      integer,  intent(in)  :: nmax    ! maximum n for m=0
      integer,  intent(in)  :: kmax    ! maximum n for any m
      integer,  intent(in)  :: imax    ! number of longitude points
      integer,  intent(in)  :: nspects ! numbers of spects per level
      integer,  intent(in)  :: nlevels ! number of horiz surfaces
      integer,  intent(in)  :: nfields ! number of fields

      integer, intent(in)     :: mtrunc(0:nmax) ! maximum m for each n-m
      integer, intent(in)     :: index(0:mmax,0:nmax) ! ordering for scoefs
      complex(r8), intent(in) :: scoefs(nspects,nlevels,nfields) 

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: power(0:kmax,nlevels,nfields)

! !DESCRIPTION: 
!  Compute power spectra from spectral coefficients, by summing over all
!  zonal wavenumbers m for each value of n. Formally, the power spectra 
!  are determined by summing over -m,...,0,...,m  but since thhe spectral
!  coefficients for m<0 are complex conjugates of those m >= 0, the 
!  summation need only be done for non-negative values, but for 
!  contributions by all m /= 0, a factor of 2 is introduced. If the 
!  spectral coefficients are identified as determined from a vector 
!  field (by specifying input name='V'), then the power is determined as 
!  if for kinetic energy from spectral coefficients of vorticity and 
!  divergence (except for the factor 1/2 in the definition of energy, 
!  as in KE = 1/2 * (u**2 +v**2)).  If spectral coefficients of vorticity 
!  are input, but the field is indicated as a scalar one (name='S' input)  
!  instead, then enstropy spectra rather than kinetic energy spectra 
!  will be calculated.
!
!  Results can be computed for multiple fields, in which case nfields
!  should be specified accordingly, and the input name should have
!  length nfields, with a character string of 'S's or 'V's indicating 
!  whether the power should be computed as for a scalar field or vector 
!  wind field, with computation for the latter as for kinetic energy.
!
!  Note that if 2*mmax=imax, then coefs for zonal wavenumber -max are
!  identical to those for mmax, and should not be counted separately
!  (i.e., the contribution by this zonal wavenumber to the power 
!  should not be multiplied by 2).  This value of mmax, however, is not
!  permitted if the vector wind is a transformed field. 
!
!  Note that the power spectra is only computed for n=0,...,nmax.
!  Thus, if thespectral truncation is not triangular, power that resides
!  in other spectral components will not be counted.
!
! !SEE ALSO: sh_dgrid_init for description of arrays mtrunc and index
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      logical lvector

      integer  :: nfield,nlev,nm,m1,n,m,j
      real(r8) :: aa,afac,xmfac,cr,ci
!      
      aa=earthr*earthr
!
      do nfield=1,nfields
        if (name(nfield:nfield).eq.'S') then
          lvector=.false.
        else
          lvector=.true.
        endif
!
        do nlev=1,nlevels
          do nm=0,nmax      ! n-m
!
!  Compute power for vector (scoefs of vort and divg) fields
            if (lvector) then
              if (nm.eq.0) then
                m1=1              ! no power in mode n=m=0 for vectors 
              else
                m1=0
              endif      
              do m=m1,mtrunc(nm) 
                j=index(m,nm)
                n=nm+m         
                afac=aa/(n*(n+1)) ! change e.g., enstropy to energy
                if (m == 0) then  ! account for absence of m=-0
                  xmfac=1.
                elseif (2*m == imax) then 
                  xmfac=1.        ! account for redundancy of m=-mmax
                else
                  xmfac=2.        ! account for m=-m conjugate waves 
                endif
                cr=real( scoefs(j,nlev,nfield))
                ci=aimag(scoefs(j,nlev,nfield))
                power(n,nlev,nfield)=power(n,nlev,nfield)+afac* &
                  xmfac*(cr*cr+ci*ci)
              enddo   
!
!  Compute power for scalar fields            
            else  
              do m=0,mtrunc(nm)
                j=index(m,nm)                   
                n=nm+m
                if (m == 0) then
                  xmfac=1.         ! account for absence of m=-0
                elseif (2*m == imax) then 
                  xmfac=1.         ! account for redundancy of m=-mmax
                else
                  xmfac=2.         ! account for m=-m conjugate waves 
                endif
                cr=real( scoefs(j,nlev,nfield))
                ci=aimag(scoefs(j,nlev,nfield))
                power(n,nlev,nfield)=power(n,nlev,nfield)+ &
                  xmfac*(cr*cr+ci*ci)
              enddo               
!
            endif  ! test on lvector
          enddo    ! end loop over nm
        enddo      ! end loop over nlev
      enddo        ! end loop over nfield
!
      end subroutine sh_calc_power
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_clean --- Routine for deallocation of some arrays
!                     
!
! !INTERFACE:
!    
      subroutine sh_clean

!USES:

      implicit none

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!

      deallocate (jfirstX,jlastX)
      deallocate (cospX,sinpX)
      deallocate (epsiX,abfaclatsX)
      deallocate (tableX,table_fftX)
      deallocate (shiftX)
      deallocate (n4alpX)
!
      if (allocated(gweightsX)) then
        deallocate (gweightsX)
      endif
!
      if (allocated(wmatrixX)) then
        deallocate (wmatrixX)
      endif
!
      if (allocated(alpX)) then
        deallocate (alpX,nalpX)
      endif
!
      end subroutine sh_clean
!
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_coslat --- Compute cosine latitude from sine latitude
!                     
!
! !INTERFACE:
!    
      subroutine sh_coslat (jmax, sinlat, coslat)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: jmax
      real(r8), intent(in)  :: sinlat(jmax)

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: coslat(jmax)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  Replace coslat that has been computed by the finite-difference 
!  delta(sin)/delta(lat) by the exact derivative d(sin)/d(lat). 
!  It is necessary to use the exact form when the wind is weighted
!  by coslat since the relationship between psi, chi and wind involves
!  derivatives of Legendre polynomials.

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
      integer   :: j
!
      do j=1,jmax
        coslat(j)=one-sinlat(j)*sinlat(j)
        if (coslat(j) <= zero) then
          coslat(j)=zero
        else
          coslat(j)=sqrt(coslat(j))
        endif 
      enddo
!
      end subroutine sh_coslat
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_c2fm --- project spectral coefficients of alp onto fields
!                     
!
! !INTERFACE:
!    
      subroutine sh_c2fm (cwork,nlast,jfirst,jlast,jalp,nlevels,alp,fm)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  ::  nlast
      integer,  intent(in)  ::  jfirst
      integer,  intent(in)  ::  jlast
      integer,  intent(in)  ::  jalp
      integer,  intent(in)  ::  nlevels
      real(r8), intent(in)  ::  alp(0:nlast,jfirst:jalp)
      complex(r8), intent(in)   :: cwork(0:nlast,nlevels)

! !OUTPUT PARAMETERS:

      complex(r8), intent(out)  :: fm(nlevels,jfirst:jlast)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!
!  Note that if the grid is symmetric about the equator, only 1/2 the 
!  alp are needed (plus the equator if it is on the grid).

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer    :: j    ! latitudinal index
      integer    :: n
      integer    :: nf
      integer    :: jsum ! sum of last and first j indexes
      real(r8)   :: afac ! alp with sign accounting for sym/asym
!
      jsum=jlast+jfirst
!
      fm(:,:)=cmplx(zero,zero,r8)   
      do j=jfirst,jlast
        do n=0,nlast

          if (j>jalp) then
            if (mod(n,2).eq.0) then
              afac= alp(n,jsum-j)
            else 
              afac=-alp(n,jsum-j)
            endif
          else
            afac=alp(n,j)
          endif

          do nf=1,nlevels
            fm(nf,j)=fm(nf,j)+afac*cwork(n,nf)
          enddo
        enddo
      enddo
!
      end subroutine sh_c2fm 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_c2s_s --- Increment spectra coefficents for scalar field
!                     
!
! !INTERFACE:
!    
      subroutine sh_c2s_s (scoefs,nspects,nmax,mmax,nlevels,m, &
                           nlast,index,cwork)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  ::  nmax
      integer,  intent(in)  ::  mmax
      integer,  intent(in)  ::  nlevels
      integer,  intent(in)  ::  nspects
      integer,  intent(in)  ::  m
      integer,  intent(in)  ::  nlast
      integer,  intent(in)  ::  index(0:mmax,0:nmax) 
      complex(r8), intent(in)  :: cwork(0:nlast,nlevels)
   
! !OUTPUT PARAMETERS:

      complex(r8), intent(inout) :: scoefs(nspects,nlevels,1)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!
! Add cwork to differently organized array of spectral coefficients 
! for scalar fields for 1 m.

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!  04Aug2008  R. Todling Bug fix: scoefs is inout
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
      integer :: n
      integer :: nf
      integer :: i

      do n=0,nlast
        i=index(m,n)
        do nf=1,nlevels           
          scoefs(i,nf,1)=scoefs(i,nf,1)+cwork(n,nf)
        enddo
      enddo
!
      end subroutine sh_c2s_s
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_c2s_v --- Increment spectra coefficents for vector field
!                     
!
! !INTERFACE:
!    

      subroutine sh_c2s_v (scoefs,nspects,nmax,mmax,nlevels,m, &
                           nlast,index,epsi,cwork,name)

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in)  :: name
      integer,  intent(in)  ::  nmax
      integer,  intent(in)  ::  mmax
      integer,  intent(in)  ::  nlevels
      integer,  intent(in)  ::  nspects
      integer,  intent(in)  ::  m
      integer,  intent(in)  ::  nlast
      integer,  intent(in)  ::  index(0:mmax,0:nmax) 
      real(r8), intent(in)  ::  epsi(0:nlast)
      complex(r8), intent(in)  :: cwork(0:nlast,nlevels)

! !OUTPUT PARAMETERS:

      complex(r8), intent(inout) :: scoefs(nspects,nlevels,2)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! Add spectral coeffients of scalar u or v field, divided by cos(theta) 
! to coefficients of vorticity and divergence for one m

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!  04Aug2008  R. Todling Bug fix: scoefs is inout
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
!
      integer :: n
      integer :: nf
      integer :: i
      integer :: truen  ! the true index n
      integer :: id1, id2, is
      integer :: nm1, np1
      real(r8)    :: rfacm, rfac
      complex(r8) :: cfacm

      rfacm=real(m,r8)/earthr
      cfacm=cmplx(zero,rfacm,r8)
      if (name == 'U') then
        id1=2
        id2=1
        is=1
      else
        id1=1
        id2=2
        is=-1
      endif

      do n=0,nlast-1

        truen=n+m
        np1=n+1
        nm1=n-1  
        i=index(m,n)

        do nf=1,nlevels           
          scoefs(i,nf,id1)=scoefs(i,nf,id1)+cfacm*cwork(n,nf)
        enddo

        if (n>0) then
          rfac=is*(truen+1)*epsi(n)/earthr
          do nf=1,nlevels           
            scoefs(i,nf,id2)=scoefs(i,nf,id2)+rfac*cwork(nm1,nf)
          enddo
        endif

        rfac=-is*truen*epsi(np1)/earthr
        do nf=1,nlevels           
          scoefs(i,nf,id2)=scoefs(i,nf,id2)+rfac*cwork(np1,nf)
        enddo

      enddo  ! loop over n
!
      end subroutine sh_c2s_v
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_emns --- Compute factorials required to compute alp
!                     
!
! !INTERFACE:
!    
      subroutine sh_emns (epsi,nmax,mmax)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: nmax
      integer,  intent(in)  :: mmax

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: epsi(0:nmax+1,0:mmax)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  compute constants used for recursive calculation of legendre
!  polynomials.  from roger daley 9/16/83
!  calculates epsilon(n-m,m)=sqrt((n**2-m**2)/(4*n**2-1))
!  for nm=n-m between 0 and nmax and m between 0 and mmax

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer               :: m
      integer               :: n
      integer               :: n1
      integer               :: nm
      real(r8)              :: fnum
      real(r8)              :: fden
!
!
      do m=0,mmax
        if (m == 0) then
          epsi(0,0)=zero    
          n1=1
        else
          n1=0 
        endif
        do nm=n1,nmax+1
          n=m+nm
          fnum=real(n**2-m**2, r8)
          fden=real(4*n**2-1, r8)
          epsi(nm,m)=dsqrt(fnum/fden)
        enddo
      enddo
!
      end subroutine sh_emns
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_fill --- Make data for testing
!                     
!
! !INTERFACE:
!    
      subroutine sh_fill (scoefs,index,ntrunc,mmax,nmax,  &
                          nspects,nlevels,ik,ikp)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)     :: mmax, nmax, nspects, nlevels, ik, ikp
      integer,  intent(in)     :: index(0:mmax,0:nmax)
      integer,  intent(in)     :: ntrunc(0:mmax)

! !OUTPUT PARAMETERS:

      complex(r8), intent(out) :: scoefs(nspects,nlevels,ik)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: Set selected spectral coefficients for testing transforms
!               Designed to set one coefficient (real or imaginary) to a 
!               non-zero value for each field and level 

! !SEE ALSO: sh_test_coefs
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
       
      integer, parameter :: mni=20
      integer :: i, j, k, m, n, ix
      integer :: mfill(mni), nfill(mni), ifill(mni)
      real(r8)    :: cmag 
!
!  cmag is magnitude of coef, mfill,nfill is m,n(=n-m) value
!  for coef to fill, and ifill=1 or to for real or imaginary
      data cmag/100./
      data mfill/20*1/
      data nfill/0,1,2,3,4,5,6,15,16,17,18,19,20,7*0/

!      data mfill/6,6,6,6,7,7,7,7,4,4,16,16,14,7*0/
!      data nfill/1,2,12,13,1,2,10,11,10,11,0,1,0,7*0/
!      data mfill/0,0,0,0,1,1,1,1,2,2,2,3,3,7*0/
!      data nfill/1,2,18,19,1,2,17,18,1,17,18,15,16,7*0/
!      data mfill/1,5,0,1,1,0,0,1,5,1,0,1,2,7*0/
!      data nfill/4,0,5,4,3,2,5,2,0,0,2,1,1,7*0/
!      data nfill/2,2,1,0,1,0,4,3,2,0,2,1,1,7*0/
!      data mfill/0,0,4,4,2,2,0,0,4,1,0,1,2,7*0/
!      data nfill/5,1,1,0,0,3,0,5,1,0,2,1,1,7*0/
      data ifill/1,1,1,1,1,1,1,1,2,1,1,1,1,7*1/
!
!  fill coef array for testing purposes.
!  values of mfill are m=0, ... ,mmax
!  values of nfill are n=0, ... ,ntrunc(m)
!  values of ifill are i=1,2 for real, imaginary part
!
!
      scoefs(:,:,:)=cmplx(zero, zero, r8)
!
!
      do k=1,nlevels
        ix=mod(k-1,mni)+1
        m=mfill(ix)
        n=nfill(ix)
        i=ifill(ix)
        if (m.gt.mmax) m=mmax
        if (n.ge.ntrunc(m)) n=ntrunc(m)
        j=index(m,n)
        if (i == 1) then
          scoefs(j,k,ikp)=cmplx( cmag, zero, r8)
        else  
          scoefs(j,k,ikp)=cmplx( zero, cmag, r8)
        endif 
      enddo
!
      end subroutine sh_fill
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_fm2c --- integrate the product of weight*alp with fm 
!                     
!
! !INTERFACE:
!    
      subroutine sh_fm2c (cwork,nlast,jfirst,jlast,jalp,nlevels,walp,fm)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  ::  nlast
      integer,  intent(in)  ::  jfirst
      integer,  intent(in)  ::  jlast
      integer,  intent(in)  ::  jalp     
      integer,  intent(in)  ::  nlevels
      real(r8), intent(in)  ::  walp(0:nlast,jfirst:jalp)  ! weight*alp
      complex(r8), intent(in)  :: fm(nlevels,jfirst:jlast)

! !OUTPUT PARAMETERS:

      complex(r8), intent(out) :: cwork(0:nlast,nlevels)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  Note that if the grid is symmetric about the equator, only 1/2 the 
!  alp are needed (plus the equator if it is on the grid).

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
!
      integer    :: j    ! latitudinal index
      integer    :: n    ! actually n-m
      integer    :: nf
      integer    :: jsum ! sum of last and first j indexes
      real(r8)   :: afac
!
      jsum=jlast+jfirst
!
      cwork=cmplx(zero, zero, r8)    
      do j=jfirst,jlast         
        do n=0,nlast
          if (j>jalp) then
            if (mod(n,2).eq.0) then
              afac= walp(n,jsum-j)
            else 
              afac=-walp(n,jsum-j)
            endif
          else
            afac=walp(n,j)
          endif

          do nf=1,nlevels          
            cwork(n,nf)=cwork(n,nf)+afac*fm(nf,j)
          enddo
        enddo
      enddo
!
      end subroutine sh_fm2c 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_fm2s --- Compute spher. harmonic coefs from zonal wave coefs
!                     
!
! !INTERFACE:
!    

      subroutine sh_fm2s (scoefs,fm,wmatrix,name,coslat,sinlat,     &
                          abfaclats,nmax,mmax,jmax,nlevels,nspects, &
                          jfirst,jlast,ntrunc,index,ik )

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: name
      integer,  intent(in)  :: nmax, mmax, jmax, nlevels, nspects
      integer,  intent(in)  :: jfirst, jlast, ik 
      integer,  intent(in)  :: ntrunc(0:mmax), index(0:mmax,0:nmax) 
      real(r8), intent(in)  :: sinlat(jmax)
      real(r8), intent(in)  :: coslat(jmax)
      real(r8), intent(in)  :: abfaclats(jmax,0:mmax)
      real(r8), intent(in)  :: wmatrix(jmax,jmax,0:1)
      complex(r8), intent(in)  :: fm(nlevels,jfirst:jlast,0:mmax)

! epsiX is also used as a global variable

! !OUTPUT PARAMETERS:

      complex(r8), intent(out) :: scoefs(nspects,nlevels,ik)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  Compute coefficients of spherical harmonic functions given a set 
!  of zonal wave coefficients.  If the input field is a component of a vector
!  wind field, then create coefficients for vorticty and divergence. 
!  Uses a latitudinal integration of the products of associated Legendre 
!  polynomials and coefficients of zonal waves, to project onto the
!  polynomials, that preserves the orthogonality property of the functions. 
!
!  Note that if the grid is symmetric about the equator, only 1/2 the 
!  alp are needed (plus the equator if it is on the grid). In this case,
!  jalp < jlast.

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer    :: m
      integer    :: j
      integer    :: nt   ! max of n-m given n
      integer    :: ntpx ! either nt (if scalar field) or nt+1 (if vector)
      integer    :: n,j1,j2,m01
      integer    :: jalp   ! last of the latitudes needed for alp values 
      integer    :: jsum   ! =jfirst+jlast
      real(r8)   :: afac
      real(r8),    allocatable :: alp(:,:)
      real(r8),    allocatable :: alp2(:,:)
      complex(r8), allocatable :: cwork(:,:)
!    
! Determine number of needed latitude values of alp
      jsum=jfirst+jlast
      if (symmetricX) then
        jalp=jfirst+(jlast-jfirst)/2
      else
        jalp=jlast
      endif
!
! Loop over zonal wave numbers m
      do m=0,mmax
!
        nt=ntrunc(m) 
        if ( name == 'S') then ! truncation for scalar field
          ntpx=nt
        else                   ! truncation for vector field
          ntpx=nt+1
        endif
!
! Compute alp for one m but all n and latitudes
        allocate (alp(0:ntpx,jfirst:jalp))
        do j=jfirst,jalp     
          call sh_alp_1 (alp(:,j),ntpx,m,j,sinlat(j),epsiX(:,m), &
                         abfaclats(j,m),alpsaveX)
        enddo
!
! Compute product of alp matrix with weight matrix
        allocate (alp2(0:ntpx,jfirst:jalp))
        if ( grid_nameX /= 'G' ) then     
!
! Use full weight matrix
          alp2=zero
          m01=mod(m,2)
          do n=0,ntpx 
            do j1=jfirst,jlast
!
              if (j1>jalp) then
                if (mod(n,2).eq.0) then
                  afac= alp(n,jsum-j1)
                else 
                  afac=-alp(n,jsum-j1)
                endif
              else
                afac=alp(n,j1)
              endif
!  
              do j2=jfirst,jalp
                alp2(n,j2)=alp2(n,j2)+afac*wmatrix(j1,j2,m01)
              enddo
            enddo
          enddo
!
        else
!    
! Use Gaussian weights as weighting vector. (For a Gaussian grid, 
! the weight matrix is diagonal; these weights being the diagonal 
! elemments.)
          do n=0,ntpx
            do j1=jfirst,jalp
              alp2(n,j1)=alp(n,j1)*gweightsX(j1)
            enddo
          enddo
!
        endif   ! test of grid_name
        deallocate (alp)          
!
! Compute integral of the product of alp with field over set of 
! latitudes Treat m=1 for wind fields in a special way
        allocate (cwork(0:ntpx,nlevels))                 
        if ( ( grid_nameX /= 'G' ) .and. (name /= 'S') &
             .and. (m==1) ) then
          call sh_m1 (ntpx,nmax,jfirst,jlast,jmax,jalp,nlevels, &
                      sinlat,epsiX(:,1),alp2,fm(:,:,m),cwork)
        else
          call sh_fm2c (cwork,ntpx,jfirst,jlast,jalp,nlevels, &
                        alp2,fm(:,:,m))
        endif
        deallocate (alp2)          

! If scalar fields,  copy cwork to properly organized array of 
!   spectral coefficients, otherwise add contribution to appropriate
!   spectral cofficients of vorticity and divergence 
        if (name == 'S') then
          call sh_c2s_s (scoefs,nspects,nmax,mmax,nlevels,m, &
                             ntpx,index,cwork)
        else 
          call sh_c2s_v (scoefs,nspects,nmax,mmax,nlevels,m, &
                             ntpx,index,epsiX(:,m),cwork,name)
        endif

        deallocate (cwork)
      enddo ! loop over m
!
      end subroutine sh_fm2s
!    
!!---------------------------------------------------------------------------
!! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!!---------------------------------------------------------------------------
!!BOP
!!
!! !ROUTINE: sh_fft --- Calls fast Fourier transform routines
!!                     
!!
!! !INTERFACE:
!!    
!      subroutine sh_fft (isign, imax, lot, x, ldx, y, ldy)
!
!!USES:
!
!      implicit none
!
!! !INPUT PARAMETERS:
!
!      integer, intent(in) :: isign
!      integer, intent(in) :: imax
!      integer, intent(in) :: lot
!      integer, intent(in) :: ldx
!      integer, intent(in) :: ldy
!
!! !OUTPUT PARAMETERS:
!
!! !INPUT/OUTPUT PARAMETERS:
!
!      real(r8)  :: x(ldx,lot)
!      complex(r8)  :: y(ldy,lot)
!
!! !DESCRIPTION: 
!!  This routine calls an FFT from a librrary (-lscs at Goddard). If an 
!!  appropriate FFT routine is not available, then call sh_fft_slow 
!!  in its place, which is a standard slow transform.
!
!! !SEE ALSO: 
!!
!
!! !REVISION HISTORY:
!!
!!  04Jun2003  R. Errico  Initial algorithm
!!  04Oct2004  R. Errico  Module development
!!
!!EOP
!!-------------------------------------------------------------------------
!!
!!    local variables
!
!      external dzfftm, zdfftm
!      integer             :: isys(0:1)
!      real(r8)            :: work(imax+2)
!      real(r8)            :: scale
!      integer             :: ntable
!
!      isys(0)=1
!
!      if (isign == 0) then      ! set a table for later use
!        ntable=imax*2+256
!        allocate (tableX(1,1,1))
!        allocate (table_fftX(ntable))
!        call dzfftm (isign, imax, lot, scale, x, ldx, y, ldy, &
!                     table_fftX, work, isys )
!      elseif (isign == 1) then  ! transform zonal coefs to fields
!        scale=one 
!        call zdfftm (isign, imax, lot, scale, y, ldy, x, ldx, &
!                     table_fftX, work, isys )
!      elseif (isign == -1) then ! transform fields to zonal coefs 
!        scale=one/real(imax, r8) 
!        call dzfftm (isign, imax, lot, scale, x, ldx, y, ldy, &
!                     table_fftX, work, isys )
!      endif
! 
!      end subroutine sh_fft
!!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_fft_slow --- Calls slow Fourier transform routines
!                     
!
! !INTERFACE:
!    
      subroutine sh_fft_slow (isign, imax, lot, x, ldx, y, ldy)
!
!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: isign
      integer, intent(in) :: imax
      integer, intent(in) :: lot
      integer, intent(in) :: ldx
      integer, intent(in) :: ldy

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

      real(r8)  :: x(ldx,lot)
      complex(r8)  :: y(ldy,lot)

! !DESCRIPTION: 
!  This routine is currently calling a regular (slow) Fourier 
!  transform.  imax can be any number, unlike in ffts

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer, save       :: imax2
      real(r8)            :: scale 

      if (isign == 0) then      ! set sine and cosine functions used later
        imax2=imax/2
        if (imax2*2 /= imax) then
          imax2=imax+1
        endif 
        allocate (table_fftX(1))
        allocate (tableX(2,0:imax2,imax))
        call sh_ft_slow_init (tableX,imax,imax2)
      elseif (isign == 1) then  ! transform zonal coefs to fields
        scale=one 
        call sh_ft_slow_m2f (imax,lot,scale,y,ldy,x,ldx,tableX,imax2)
      elseif (isign == -1) then ! transform fields to zonal coefs 
        scale=one/real(imax, r8) 
        call sh_ft_slow_f2m (imax,lot,scale,x,ldx,y,ldy,tableX,imax2)
      endif
 
      end subroutine sh_fft_slow
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_ft_slow_f2m --- A slow Fourier trans.: fields to zonal coefs 
!                     
!
! !INTERFACE:
!    
      subroutine sh_ft_slow_f2m (imax,lot,scale,x,ldx,y,ldy,table,imax2)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: imax
      integer,  intent(in) :: imax2
      integer,  intent(in) :: lot             ! number of transforms
      integer,  intent(in) :: ldx
      integer,  intent(in) :: ldy
      real(r8), intent(in) :: scale
      real(r8), intent(in) :: table(2,0:imax2,imax)  
      real(r8), intent(in) :: x(0:ldx-1,lot)   ! real fields  

! !OUTPUT PARAMETERS:

      complex(r8),  intent(out) :: y(0:ldy-1,lot) ! zonal wave coefficients

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! A slow algorithm to compute zonal transforms: from real fields to 
! complex zonal wave coefficients.
! An FFT that allows imax with factors of 2*(2**m)*(3**n)*(5**k) should
! be called in place of this slow algorithm. This routine is also called 
! from routine sh_fft if isign=0, but actually a set up routine to create 
! the array table should be called in its place.
! 
!
! !SEE ALSO: sh_ft_slow_m2f, sh_ft_slow_init, sh_fft  
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer             :: m,j,n
      complex(r8)         :: wjm
!      
      y=cmplx(zero,zero,r8)
!    
      do m=0,ldy-1
        do j=0,ldx-1
          wjm=cmplx(table(1,m,j+1),-table(2,m,j+1), r8)
          do n=1,lot
            y(m,n)=y(m,n)+x(j,n)*wjm
          enddo
        enddo
      enddo
      y=scale*y
!
      end subroutine sh_ft_slow_f2m 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_ft_slow_init --- Set trig factors for slow Fourier trans.
!                     
!
! !INTERFACE:
!    
      subroutine sh_ft_slow_init (table,imax,imax2)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: imax
      integer,  intent(in) :: imax2

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: table(2,0:imax2,imax)  

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! Pre-compute sines and cosines of k*x;  k=m*2*pi ; x=j/imax
! j=1,...,imax; m=0,..., imax2=imax/2 (+1 if imax is odd)
! for use in slow Fourier transform routine.
!
! !SEE ALSO:  sh_ft_slow_f2m, sh_ft_slow_m2f, sh_fft  
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer             :: m      ! integral zonal wavenumber
      integer             :: j      ! zonal grid point index
      real(r8)            :: pifac  ! 2*pi/imax
      real(r8)            :: pms    ! k*x
!      
      pifac=(two**3)*datan(one)/imax
      do m=0,imax2
        pms=pifac*m
        do j=1,imax
          table(1,m,j)=cos(pms*j)
          table(2,m,j)=sin(pms*j)
        enddo
      enddo
!
      end subroutine sh_ft_slow_init 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_ft_slow_m2f --- A slow Fourier trans.: zonal coefs to fields
!                     
!
! !INTERFACE:
!    
      subroutine sh_ft_slow_m2f (imax,lot,scale,y,ldy,x,ldx,table,imax2)

!USES:


! !INPUT PARAMETERS:

      integer,  intent(in) :: imax
      integer,  intent(in) :: imax2
      integer,  intent(in) :: lot
      integer,  intent(in) :: ldx
      integer,  intent(in) :: ldy
      real(r8), intent(in) :: scale
      real(r8), intent(in) :: table(2,0:imax2,imax)  
      complex(r8),  intent(in) :: y(0:ldy-1,lot)

! !OUTPUT PARAMETERS:

      real(r8),intent(out) :: x(0:ldx-1,lot)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! A slow algorithm to compute zonal transforms: from complex zonal wave 
! coefficients to real fields. 
! An FFT that allows imax with factors of 2*(2**m)*(3**n)*(5**k) should
! be called in place of this slow algorithm. 

! !SEE ALSO: sh_ft_slow_f2m, sh_ft_slow_init, sh_fft  
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer             :: m,j,n
      real(r8)            :: sfac
      complex(r8)         :: wjm
!
      x=zero
      sfac=scale
!
      do m=0,ldy-1
        do j=0,ldx-1
          wjm=cmplx(table(1,m,j+1),table(2,m,j+1), r8)
          do n=1,lot
            x(j,n)=x(j,n)+sfac*real(wjm*y(m,n))
          enddo
        enddo
        sfac=two*scale ! factor applied when m>0
      enddo
!
      end subroutine sh_ft_slow_m2f 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_gauss --- Compute Gaussian latitudes weights 
!                     
!
! !INTERFACE:
!    
      subroutine sh_gauss (a,w,k,iordlt)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)    :: k       ! number of latitudes
      integer, intent(in)    :: iordlt  ! ordering instruction

! !OUTPUT PARAMETERS:

      real(r8), intent(out)  :: a(k)    ! sines of Gaussian latitudes
      real(r8), intent(out)  :: w(k)    ! Gaussian weights

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 

! The weights are renormalized at the end of this routine such that their
! sum from pole to pole equals 1.  
! This is consistent with a change in normalization of Legendre Polynomials
! that are normalized such that the integral of their squared value from 
! pole to pole equals 1. This renormalization is consistent with a weighting
! that represents a fraction of the surface area on the sphere. 


! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  pre  1981  NCAR       Initial algorithm acquired from ECMWF or BMRC?
!  ?????1981  R. Errico  Develped for normal mode software
!  04Jun2003  R. Errico  Initial algorithm at GMAO
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

!  iordlt=0 if order is north to south;
!  iordlt=1 if order is south to north;
!  iordlt=2 if order is alternating north, south, starting at poles.

      integer     :: kk,is,iter,n,l,n1,n2
      real(r8)    :: worka(k),facsin,eps,c,fk
      real(r8)    :: xz,pkm2,pkm1,fn,pk,pkmrk,sp,avsp
!     
      facsin=45./atan(one)
      if (r8 == 4) then 
        eps=1.e-8
      else
        eps=1.d-14
      endif 
      c=(one-(two/3.14159265358979)**2)/(two*two)
      fk=k
      kk=k/2
      call sh_bssl (a,kk)

      do 3 is=1,kk
      xz=cos(a(is)/sqrt((fk+one/two)**2+c))
      iter=0
    1 pkm2=1.
      pkm1=xz
      iter=iter+1
      if (iter.gt.10) then
         print 100
         stop
      endif
      do 2 n=2,k
      fn=n
      pk=((two*fn-one)*xz*pkm1-(fn-one)*pkm2)/fn
      pkm2=pkm1
      pkm1=pk
    2 continue
      pkm1=pkm2
      pkmrk=(fk*(pkm1-xz*pk))/(one-xz**2)
      sp=pk/pkmrk
      xz=xz-sp
      avsp=abs(sp)
      if (avsp.gt.eps) goto 1
      a(is)=xz
!  original code
!           w(is)=(two*(one-xz**2))/(fk*pkm1)**2
!  code from M. Ehrendorfer:
      w(is)=two/((one-xz*xz)*pkmrk*pkmrk)

    3 continue

      if (k.ne.kk*2) then
         a(kk+1)=0.
         pk=2./fk**2
         do 4 n=2,k,2
         fn=n
         pk=pk*fn**2/(fn-one)**2
    4    continue
         w(kk+1)=pk
      endif

      do 5 n=1,kk
      l=k+1-n
! replaces sin lat by lat itself in degrees
!  a(n)=facsin*asin(a(n))
      a(l)=-a(n)
      w(l)=w(n)
      worka(n)=a(n)
      worka(l)=w(n)
    5 continue
      if (k.ne.kk*2) worka(kk+1)=w(kk+1)

      if (iordlt.eq.0) return
      if (iordlt.eq.1) then
         do 11 n=1,kk
         l=k+1-n
         a(n)=-worka(n)
         a(l)= worka(n)
         w(n)= worka(l)
         w(l)= worka(l)
   11    continue
         if (k.ne.kk*2) then
            a(kk+1)=0.
            w(kk+1)=worka(kk+1)
         endif
      elseif (iordlt.eq.2) then
         do 12 n=1,kk
         l=k+1-n
         n1=n*2
         n2=n1-1
         a(n2)= worka(n)
         a(n1)=-worka(n)
         w(n2)= worka(l)
         w(n1)= worka(l)
   12    continue
         if (k.ne.kk*2) then
            a(k)=0.
            w(k)=worka(kk+1)
         endif
      else
         print 101,iordlt
         stop
      endif
!
! Renormalize weights such that sum from pole to pole = 1
!
      w(:)=w(:)/two
!
  100 format(//,5x,'x x x x x  convergence failure in shgaus')
  101 format(//,5x,'x x x x x  bad option in shgaus      iordlt=',i4)
      end subroutine sh_gauss
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_hmns  --- Routine to compute derivatives of alp
!                     
!
! !INTERFACE:
!    
      subroutine sh_hmns (dalp,alp,nlast,m,epsi)

!USES:

      implicit none

! !INPUT PARAMETERS:
      integer,  intent(in)   :: nlast      
      integer,  intent(in)   :: m
      real(r8), intent(in)   :: epsi(0:nlast+1)
      real(r8), intent(in)   :: alp(0:nlast+1)

! !OUTPUT PARAMETERS:

      real(r8), intent(out)  :: dalp(0:nlast)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  Compute derivatives of associated Legendre polynomials divided by 
!  cos**2  (?) 
!  for given m and latitude. The index of dalp is the value n-m
!  THIS CODE HAS NOT BEEN CHECKED HERE SINCE IT IS NOT USED IN THIS
!  FORMULATION THAT TRANSFORMS U AND V DIRECTLY

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer                :: n
      real(r8)               :: alpilm, truen  
!
      do n=0,nlast
        truen=n+m
        alpilm=0.
        if(n>0) alpilm=alp(n-1)
        dalp(n)=(truen+one)*epsi(  n)*alpilm -  &
                      truen*epsi(n+1)*alp(n+1)
      enddo
!
      end subroutine sh_hmns
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_initial_all --- Set some initial arrays, called by sh_dgrid_init
!                     
!
! !INTERFACE:
!    
      subroutine sh_initial_all (mtrunc,ntrunc,index,imax,jmax,mmax, &
                                 nmax,kmax,icoft,nspects,func,ierror) 

!USES:

      implicit none

! !INPUT PARAMETERS:

!
      character(len=*), intent(in) :: func ! 's2f', 'f2s', or 'both'
!             only if func /= 's2f' is wmatrix needed  
      integer,  intent(in)  :: imax    ! number of longitudes
      integer,  intent(in)  :: jmax    ! 2nd dimension of field array
                                       ! = # lats for Scalar grid
      integer,  intent(in)  :: mmax    ! maximum zonal wave number
      integer,  intent(in)  :: nmax    ! maximum n for m=0
      integer,  intent(in)  :: kmax    ! maximum m-n for any m
      integer,  intent(in)  :: icoft   ! how coefficients will be ordered
                                       ! see sh_trunc 

! !OUTPUT PARAMETERS:

      integer,  intent(out)  :: ierror               ! error flag (0=OK)
      integer,  intent(out)  :: nspects              ! number of coefs for
!                                                      each level and field   
      integer,  intent(out)  :: mtrunc(0:nmax)       ! max m for each n-m
      integer,  intent(out)  :: ntrunc(0:mmax)       ! max n-m for each m
      integer,  intent(out)  :: index(0:mmax,0:nmax) ! order of m,n-m values
!
! Also set the global variables and arrays:
!     abfaclatsX  latitudinal dependent factors used to compute alp
!     wmatrixX    weighting matrix for projection onto alp functions
!     alpX        associated Legendre polynomials
!     nalpX       array for describing m,n ordering of alp functions
!     jalpX       number of lats needed to store required alp

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! Set some resolution-dependent constant arrays need by the 
! spectral transform software
!
! The spectral resolution is assumed to be pentagonal, defined by the three
! variables mmax, nmax, and kmax. The latter must satisfy
!                nmax <= kmax <= mmax+nmax
! For triangular truncation, mmax=nmax=kmax. For rhomboidal truncation,
! nmax=mmax=kmax/2. mmax is the maximum zonal wave number. nmax is the maximum 
! order of the Legendre polynomials for zonal wavenumber m=0. kmax is the
! maximum n-m for any m.
!
! The following is an example for mmax=nmax=5, with kax=7:
! The numbers from 1-30 in the chart below indicate the 
! ordering of the spectral coefficients from 1-30.  This ordering is 
! stored in the array index(m,n-m). In the example below, the spectral
! coefficient for m=2, n-m=2 is stored in an array with first index=15. 
!
! ntrunc=  5   5   5   4   3   2    max value of n-m in each column   
! 
! kmax  7         18  23  27  30   
!       6     12  17  22  26  29
! nmax= 5  6  11  16  21  25  28
!     n=4  5  10  15  20  24
!       3  4   9  14  19
!       2  3   8  13
!       1  2   7 
!     n=0  1
!       m= 0   1   2   3   4   5=mmax   
!
! mtrunc is the maximum m for each n-m (i.e., along diagonals).
! In the above example, mtrunc = 5, 5, 5, 4, 3, 2
!
!  If the conditions for pentagonal truncation are not met, error=1,
!  otherwise 0 means no error detected.
!
! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer      :: ig      ! grid number
      integer      :: m1n1, mnk
      real(r8)     :: x(2,1)  ! dummy variables needed by fft
      complex(r8)  :: y(2,1)  ! dummy variables needed by fft
!
      ierror=0
      if ( (kmax < nmax) .or. (kmax > mmax+nmax) ) then
        ierror=1
        print *,' '
        print *,' X X X X X X X X X X X X X X X X X X X X X X X X X'
        print *,'                 E R R O R                        '
        print *,' sh_initial_all detects incompatible values'
        print *,' kmax must be >= nmax' 
        print *,' kmax must be <= mmax+nmax' 
        print *,' mmax=',mmax,' nmax=',nmax,'   kmax=',kmax
      endif
!
! Compute number of spectral coefficients for pentagonal truncation
      mnk=mmax+nmax-kmax
      m1n1=(mmax+1)*(nmax+1)
      nspects=m1n1-mnk*(mnk+1)/2
!
! Compute arrays describing truncation and spectral ordering
      call sh_trunc (mtrunc,ntrunc,index,mmax,nmax,kmax,icoft)
!
! Compute factors required by FFT routines
      call sh_fft_slow (0, imax, 1, x, 1, y, 1)
!
! Allocate some arrays that are global to the module
! Note that the grid is Gaussian, wmatrixX is referenced as an 
! array argument, but is otherwise unused.
      allocate (epsiX(0:nmax+1,0:mmax))
      allocate (abfaclatsX(1:jmax,0:mmax,igridsX))
      allocate (wmatrixX(jmax,jmax,0:1,igridsX))
!
! Compute factorials used for computing Legendre polynomials
      call sh_emns (epsiX,nmax,mmax)
!
! Compute latitude and m factors used for computing associated 
! Legendre polynomials. 
      do ig=1,igridsX
        call sh_abfaclats (abfaclatsX(:,:,ig),sinpX(:,ig),      &
                           mmax,jmax,jfirstX(ig),jlastX(ig))
      enddo
!
! Save alp values for all subgrids if this option has been selected
! for the purpose of omitting some costly redundant calculations.
      if (alpsaveX) then
        call sh_alp_all (mmax,jmax,ntrunc)
      endif
!
! Compute weighting matrix for projections onto alp for m=0,1  
! to be used for all even and odd m.  If the grid is Gaussian, 
! this matrix reduces to a diagonal, whose values are the vector 
! of Gaussian weights, which is set elsewhere.
      if ( (grid_nameX /= 'G' ) .and. (func /= 's2f') ) then
        do ig=1,igridsX
          call sh_wmatrix (wmatrixX(:,:,0:1,ig),                &
                           sinpX(:,ig),abfaclatsX(:,:,ig),0,1,  &
                           mmax,jmax,jfirstX(ig),jlastX(ig)) 
        enddo
      endif
!
      end subroutine sh_initial_all
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_initial_D --- Set and save D-grid arrays needed later
!                     
!
! !INTERFACE:
!    
      subroutine sh_initial_D (nmax,mmax,kmax,imax,jmax,nspects, &
                               ntrunc,mtrunc,index,func,alpsave, & 
                               lprint,ierror)

!USES:

      implicit none

! !INPUT PARAMETERS:

      logical,  intent(in)  :: lprint      ! .true. if info to be printed
      logical,  intent(in)  :: alpsave     ! .true. if alps are to be saved
      character(len=*), intent(in) :: func ! 's2f', 'f2s', or 'both'
                      ! func indicates the function of the programs; i.e.,
                      ! spectra-to-fields, the reverse, or both will be done 
      integer,  intent(in)  :: imax    ! number of longitudes
      integer,  intent(in)  :: jmax    ! 2nd dimension of field array
                                       ! = # lats for Scalar grid
      integer,  intent(in)  :: mmax    ! maximum zonal wave number
      integer,  intent(in)  :: nmax    ! maximum n for m=0
      integer,  intent(in)  :: kmax    ! maximum n-m for any m

! !OUTPUT PARAMETERS:

      integer, intent(out) :: ierror          ! error flag (=0 means OK)
      integer, intent(out) :: nspects         ! number of spectral coefs
      integer, intent(out) :: mtrunc(0:nmax)  ! max m for given n-m
      integer, intent(out) :: ntrunc(0:mmax)  ! max n-m for given m
      integer, intent(out) :: index(0:mmax,0:nmax) ! ordering of coefs 

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! This routine must be called by the user before any transforms are 
! to be attempted.  If the spectral or grid resolution is changed, 
! it must be called again with the new resolution parameters.  
!
! This routine sets some arrays that describe the spectral resolution 
! and the ordering of the spectral coefficients.
! It also computes some module-resident arrays that are saved for later 
! calculation of the associated Legendre polynomials.
!
! This is specifically for the D-grid as used in the NASA/GMAO FVGCM.
! This specification applies to the setting of jfirstX, jlastX, shiftX, 
! sinpX, and symmetricX.  The other grid-dependent arrays produced here 
! use these values, but are otherwise generically produced for any 
! values input by this routine and what it itself calls. 
!
! The following conditions must be satisfied:
!     imax even and >2
!     jmax >= kmax+4
!     0 <= mmax < imax/2
!     nmax >= 0
!     kmax >= mmax    kmax >= nmax  kmax <= mmax+nmax 
! 
! ierror should be zero on return
!
! symmetricX is set to .true. here since the latitudes are spaced 
! symmetrically about the equator for all the fields. This symmetry
! requires computation of Legendre polynomials for 1 hemisphere only.
!
! alpsave = .true. causes the Legendre polynomials to be computed once
!   for all m, n, and lats, requiring more memory, but many fewer 
!   redundant calculations. Both memory and calculations scale as 
!   nmax*mmax*jmax.
 
! !SEE ALSO: sh_test_coefs for example of use
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
    
      integer   :: ier
      real(r8)  :: dl
!
      ierror=0          ! this is added to if errors detected
      igridsX=3         ! u,v,T on separate grids
      grid_nameX='D'    ! this is variable initilization for the D grid 
      alpsaveX=alpsave  ! indicates whether to pre-compute the alp
      n3alpX=2          ! number of independent sub-grids needed for alpX
      allocate (n4alpX(n3alpX))
      n4alpX(1)=1   ! specifies that alp grid 1 is from sub-grid 1 (u)
      n4alpX(2)=3   ! specifies that alp grid 2 is from sub-grid 3 (T)
                    ! because the sub-grid for v is inbedded in 3
!
! If some saved arrays have already been allocated, the assumption 
! here is that the must be re-determined here for a different 
! resolution and the former arrays can be discarded.
      if (allocated(sinpX)) then
        call sh_clean
      endif 
      allocate (jfirstX(igridsX),jlastX(igridsX))
      allocate (cospX(jmax,igridsX),sinpX(jmax,igridsX))
      allocate (shiftX(igridsX))
!
! Set sines and cosines of grid points
      symmetricX=.true.  ! this grid is equatorially symmetric 
      call sh_setrig_D (imax, jmax, dpX, dl, cospX(:,2), cospX(:,1), &
                                       sinpX(:,2), sinpX(:,1) )
! Replace finite-difference cosines by their exact counterparts
! as required for the u,v spectral transforms.
      call sh_coslat (jmax,sinpX(:,1),cospX(:,1))
      call sh_coslat (jmax,sinpX(:,2),cospX(:,2))
      jfirstX(1)=2      ! u grid
      jlastX(1)=jmax
      jfirstX(2)=2      ! v grid
      jlastX(2)=jmax-1  
      jfirstX(3)=1      ! S grid (for scalars such as T, q)
      jlastX(3)=jmax
      cospX(:,3)=cospX(:,2) 
      sinpX(:,3)=sinpX(:,2)
!
! Shift specifies how the first (i=1) longitudes of the u and T fields 
! are displaced 1/2 grid distance to the east of the Prime Meridian,
! whereaa i=1 for the v field is at 0 degrees.
      shiftX(1)=one/(imax*two)
      shiftX(2)=zero
      shiftX(3)=one/(imax*two)
!
! Set some arrays for describing truncation, order of coefficients, 
! factors used to compute associated Legender polynomials, and the
! weight matrix used to integrate the fields
      call sh_initial_all (mtrunc,ntrunc,index,imax,jmax,mmax,  &
                           nmax,kmax,0,nspects,func,ier)
      ierror=ierror+ier
!
!      if (lprint) then
!        print *,' '
!        print ('(a)'),' Spectral arrays initialized for a D grid:'
!        print ('(a,6i6)'),'   imax,jmax,mmax,nmax,kmax,nspects = ', &
!                              imax,jmax,mmax,nmax,kmax,nspects 
!        print ('(a,6i4)'),'   jfirst(1:3),jlast(1:3) = ',jfirstX,jlastX
!        print ('(a,3f10.4)'),'   shift(1:3) = ',shiftX   
!      endif
!
      if ( (imax <= 2*mmax) .or. (jmax < kmax+4) ) then
        ierror=ierror+10
        print *,' '
        print *,' X X X X X X X X X X X X X X X X X X X X X X X X X'
        print *,'                 W A R N I N G                    '
        print *,' sh_initial_D detects possibly incompatible values'
        print *,' imax must be > 2*mmax for Fourier transforms'
        print *,' jmax must be >= kmax+4 for field to scoefs of wind'
        print *,' jmax must be >= kmax+1 for field to scoefs of scalar'
        print *,' imax=',imax,'   mmax=',mmax
        print *,' jmax=',jmax,'   nmax=',nmax
      endif
!
      end subroutine sh_initial_D 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_initial_G --- Set and save G-grid arrays needed later
!                     
!
! !INTERFACE:
!    
      subroutine sh_initial_G (nmax,mmax,kmax,imax,jmax,nspects, &
                               ntrunc,mtrunc,index,func,alpsave, &
                               lprint,ierror)

!USES:

      implicit none

! !INPUT PARAMETERS:

      logical,  intent(in)  :: lprint      ! .true. if info to be printed
      logical,  intent(in)  :: alpsave     ! .true. if alps are to be saved
      character(len=*), intent(in) :: func ! 's2f', 'f2s', or 'both'
                      ! func indicates the function of the programs; i.e.,
                      ! spectra-to-fields, the reverse, or both will be done 
      integer,  intent(in)  :: imax    ! number of longitudes
      integer,  intent(in)  :: jmax    ! 2nd dimension of field array
                                       ! = # lats for Scalar grid
      integer,  intent(in)  :: mmax    ! maximum zonal wave number
      integer,  intent(in)  :: nmax    ! maximum n for m=0
      integer,  intent(in)  :: kmax    ! maximum n-m for any m

! !OUTPUT PARAMETERS:

      integer, intent(out) :: ierror          ! error flag (=0 means OK)
      integer, intent(out) :: nspects         ! number of spectral coefs
      integer, intent(out) :: mtrunc(0:nmax)  ! max m for given n-m
      integer, intent(out) :: ntrunc(0:mmax)  ! max n-m for given m
      integer, intent(out) :: index(0:mmax,0:nmax) ! ordering of coefs 

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! This routine must be called by the user before any transforms are 
! to be attempted.  If the spectral or grid resolution is changed, 
! it must be called again with the new resolution parameters.  
!
! This routine sets some arrays that describe the spectral resolution 
! and the ordering of the spectral coefficients.
! It also computes some module-resident arrays that are saved for later 
! calculation of the associated Legendre polynomials.
!
! This is specifically for the G-grid (Gaussian grid), which is assumed 
! the same for all grids.
! This specification applies to the setting of jfirstX, jlastX, shiftX, 
! sinpX, and symmetricX.  The other grid-dependent arrays produced here 
! use these values, but are otherwise generically produced for any 
! values input by this routine and what it itself calls. 
! 
! The following conditions must be satisfied:
!     imax even and >2
!     jmax >= kmax+2
!     0 <= mmax < imax/2
!     nmax >= 0
!     kmax >= mmax    kmax >= nmax  kmax <= mmax+nmax 
!
! ierror should be zero on return
!
! symmetricX is set to .true. here since the latitudes are spaced 
! symmetrically about the equator for all the fields. This symmetry
! requires computation of Legendre polynomials for 1 hemisphere only.
!
! alpsave = .true. causes the Legendre polynomials to be computed once
!   for all m, n, and lats, requiring more memory, but many fewer 
!   redundant calculations. Both memory and calculations scale as 
 
! !SEE ALSO: sh_test_coefs for example of use
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
    
      integer   :: ier
!
      ierror=0         ! this is added to if errors detected
      igridsX=1        ! all fields assumed to be on the same grid
      grid_nameX='G'   ! this is variable initilization for the G grid 
      alpsaveX=alpsave ! indicates whether to pre-compute the alp
      n3alpX=1         ! number of independent sub-grids needed for alpX
      allocate (n4alpX(n3alpX))
      n4alpX(1)=1   ! specifies that alp grid 1 is from sub-grid 1 
                    ! because there is only 1 sub-grid for the G-grid
!
! If some saved arrays have already been allocated, the assumption 
! here is that the must be re-determined here for a different 
! resolution and the former arrays can be discarded.
      if (allocated(sinpX)) then
        call sh_clean
      endif 
      allocate (jfirstX(igridsX),jlastX(igridsX))
      allocate (cospX(jmax,igridsX),sinpX(jmax,igridsX))
      allocate (shiftX(igridsX))
!
! Set sines and Gaussian weights of grid points
      symmetricX=.true.  ! this grid is equatorially symmetric 
      allocate (gweightsX(jmax))
      call sh_gauss (sinpX(:,1),gweightsX(:),jmax,1)  
!
! Replace finite-difference cosines by their exact counterparts
! as required for the u,v spectral transforms.
      call sh_coslat (jmax,sinpX(:,1),cospX(:,1))  
!
! All fields are on same grid with no shift in longitudes
      jfirstX(:)=1
      jlastX(:)=jmax 
      shiftX(:)=zero 
!
! Set some arrays for describing truncation, order of coefficients, 
! factors used to compute associated Legender polynomials, and the
! weight matrix used to integrate the fields
      call sh_initial_all (mtrunc,ntrunc,index,imax,jmax,mmax,  &
                           nmax,kmax,0,nspects,func,ier)
      ierror=ierror+ier
!
      if (lprint) then
        print *,' '
        print ('(a)'),' Spectral arrays initialized for a G grid:'
        print ('(a,6i6)'),'   imax,jmax,mmax,nmax,kmax,nspects = ', &
                              imax,jmax,mmax,nmax,kmax,nspects 
      endif
!
      if ( (imax <= 2*mmax) .or. (jmax < kmax+1) ) then
        ierror=ierror+10
        print *,' '
        print *,' X X X X X X X X X X X X X X X X X X X X X X X X X'
        print *,'                 W A R N I N G                    '
        print *,' sh_initial_G detects possibly incompatible values'
        print *,' imax must be > 2*mmax for Fourier transforms'
        print *,' jmax must be >= kmax+2 for field to scoefs of wind'
        print *,' jmax must be >= kmax+1 for field to scoefs of scalar'
        print *,' imax=',imax,'   mmax=',mmax
        print *,' jmax=',jmax,'   nmax=',nmax
      endif
!
      end subroutine sh_initial_G 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_intsqr_grid --- Integrate the square of a gridded field
!                     
!
! !INTERFACE:
!    
      subroutine sh_intsqr_grid (imax,jmax,nlevs,f,E,name)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: imax
      integer, intent(in) :: jmax
      integer, intent(in) :: nlevs
      character(len=*), intent(in) :: name  ! name of field 'U', 'V', or 'S'
      real(r8), intent(in) :: f(imax,jmax,nlevs) ! 3-D field

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: E(nlevs)   ! add integrals to this array

! !DESCRIPTION: 
! Compute the global area integral of the square of a field
! The integral is computed separately for each level k on which the 
! field is defined and added to an input array E(k).  Note that this
! E(k) is not energy since, e.g. the contribution to energy by the u field
! would be 1/2 u**2.
! 

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  01Feb2005  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      real(r8) :: S,Si,latfac,sum_w,frac
      real(r8) :: grid_weights(jmax)
      integer  :: i,j,k,ig
!
! Set up weights for grid areas
      if (grid_nameX=='D') then
        if (name=='U') then
          ig=1
          do j=2,jmax
            grid_weights(j)=sinpX(j,2)-sinpX(j-1,2)
          enddo
        else if (name=='V') then
          ig=2 
          do j=2,jmax-1
            grid_weights(j)=sinpX(j+1,1)-sinpX(j,1)
          enddo
        else
          ig=3
          grid_weights(1)=one+sinpX(2,1)
          grid_weights(jmax)=one-sinpX(jmax,1)
          do j=2,jmax-1
            grid_weights(j)=sinpX(j+1,1)-sinpX(j,1)
          enddo
        endif
!
      else if (grid_nameX=='G') then
        ig=1
        grid_weights(1:jmax)=gweightsX(1:jmax)
!
      else
        print *,' '
        print *,' ERROR in call to sh_sumsqr_grid'
        print *,' grid_nameX=',grid_nameX,' not an option'
      endif 
!
! divide by number of points in zonal direction
      sum_w=zero
      do j=jfirstX(ig),jlastX(ig) 
        sum_w=sum_w+grid_weights(j)
      enddo
      frac=one/(sum_w*real(imax))
      do j=jfirstX(ig),jlastX(ig) 
        grid_weights(j)=frac*grid_weights(j)
      enddo
!
      do k=1,nlevs
        S=zero
        do j=jfirstX(ig),jlastX(ig)
          Si=zero
          do i=1,imax
            Si=Si+f(i,j,k)*f(i,j,k)
          enddo
          S=S+Si*grid_weights(j)
        enddo
        E(k)=E(k)+S  
      enddo
!
      end subroutine sh_intsqr_grid
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_m1 --- Special handeling of wind for m=1
!                     
!
! !INTERFACE:
!    

      subroutine sh_m1 (ntpx,nmax,jfirst,jlast,jmax,jalp,nlevels, &
                        sinlat,epsi,alp2,fm,cwork)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)   ::  ntpx,nmax 
      integer, intent(in)   ::  jfirst,jlast,jmax,jalp
      integer, intent(in)   ::  nlevels
      real(r8), intent(in)  ::  sinlat(jmax)
      real(r8), intent(in)  ::  epsi(0:nmax+1)
      real(r8), intent(in)  ::  alp2(0:ntpx,jfirst:jalp)
      complex(r8), intent(in)  :: fm(nlevels,jfirst:jlast)

! !OUTPUT PARAMETERS:

      complex(r8), intent(out) :: cwork(0:ntpx,nlevels)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!
!  This routine appears necessary as a special treatment of zonal 
!  wavenumber 1 for vector fields.  Although the special nature of
!  m=1 is well known, the need for the treatment here is not described
!  anywhere. It is not necessary for a gaussian grid, nor for 
!  regular lat/lon grid that includes pole points (as stated in some
!  literature) but appears necessary for the D-grid, perhaps because 
!  the wind at the poles is not independently defined on that grid. 
!  The special treatment provided by this routine works, but without it,
!  the transforms for m=1 are otherwise not exact.
!  
!  Note that in the limit as the pole is approached, the wind vector 
!  tends to zero for all zonal wavenumbers except m=1, which implies
!  a special condition applies.
!
!  The "fix" here is to first compute the Legendre transforms for
!  u*cos (or v*cos in a separate call) rather than for u/cos as 
!  input here and transformed elsewhere for other m. This function is 
!  well behaved at the poles (tends to zero there). Using these 
!  projection coefficients, this function is interpolated exactly 
!  to a Gaussian grid with approximately the same resolution. 
!  On that new grid, the function u*cos is replaced
!  by u*cos.  Then, this new function is projected onto the Legendre
!  polynomials, relying on the fact that this works fine for the 
!  Gaussian grid. This procedure should be uneccessary and no 
!  mathematical analysis has been performed to show that it should work,
!  but extensive testing demonstrates it does.
!
!  Note that if the grid is symmetric about the equator, only 1/2 the 
!  alp are needed (plus the equator if it is on the grid).

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer     :: j, jx, n
      integer     :: jmaxB    ! number of lats for Gaussian grid
      integer     :: jmaxB2   ! =jmaxB/2 (number of alp needed for G-grid 
      real(r8),    allocatable :: abfaclatsB(:,:)
      real(r8),    allocatable :: alpB(:,:)
      real(r8),    allocatable :: sinlatB(:), gweightsB(:)
      complex(r8) :: fmcos(nlevels,jfirst:jlast)
      complex(r8), allocatable :: fmB(:,:)
!
!  Replace fm=wind/cos with fmc=wind*cos on original grid
      do j=jfirst,jlast
        fmcos(:,j)=fm(:,j)*(one-sinlat(j)*sinlat(j))
      enddo
!
!  Project wind*cos(lat) onto P_n^1
      call sh_fm2c (cwork,ntpx,jfirst,jlast,jalp,nlevels,alp2,fmcos)
!
! Create new, Gaussian grid
      jmaxB=ntpx+1
      if (mod(jmaxB,2) /= 0) then
        jmaxB=jmaxB+1
      endif
      jmaxB2=jmaxB/2
!
! Compute Gaussian weights and sines of Gaussian latitudes 
      allocate (sinlatB(jmaxB),gweightsB(jmaxB))
      call sh_gauss (sinlatB,gweightsB,jmaxB,1)
!
! Compute abfaclats for Gaussian grid
      allocate (abfaclatsB(jmaxB,0:1))
      call sh_abfaclats (abfaclatsB,sinlatB,1,jmaxB,1,jmaxB)
!
! Project wind*cos onto new, Gaussian grid
      allocate (alpB(0:ntpx,jmaxB2))
      allocate (fmB(nlevels,jmaxB))
      do j=1,jmaxB2
        call sh_alp_1 (alpB(:,j),ntpx,1,j,sinlatB(j),epsi, &
                       abfaclatsB(j,1),.false.)
      enddo
      call sh_c2fm (cwork,ntpx,1,jmaxB,jmaxB2,nlevels,alpB,fmB)      
!
!  On new grid, replace wind*cos by wind/cos 
      do j=1,jmaxB
        fmB(:,j)=fmB(:,j)/(one-sinlatB(j)*sinlatB(j))
      enddo
!
!  Compute weighted alp for new grid  
      do n=0,ntpx 
        do j=1,jmaxB2
          alpB(n,j)=alpB(n,j)*gweightsB(j)
        enddo
      enddo
!
!  Compute cwork from new Gaussian grid
      call sh_fm2c (cwork,ntpx,1,jmaxB,jmaxB2,nlevels,alpB,fmB(:,:))
!
      deallocate (sinlatB,fmb,alpB,abfaclatsB)
!
      end subroutine sh_m1
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_prt_allc --- For testing, print all significant coefs
!                     
!
! !INTERFACE:
!    
      subroutine sh_prt_allc (scoefs,nn1,nn2,nn3,sname)

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: sname
      integer,  intent(in) :: nn1,nn2,nn3
      complex(r8), intent(in) :: scoefs(nn1,nn2,nn3) 

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! Print out non-zero coefs for testing purposes
! Only coefs with cr**2+ci**2 > 1.e-8*max(cr**2+ci**2) printed

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer  :: n1,n2,n3
      real(r8) :: cmags, cmagsmax
      real(r8), parameter :: cprecision=1.e-8    
!
! detemine max magnitude
      cmagsmax=zero      
      do n3=1,nn3
        do n2=1,nn2
          do n1=1,nn1
            cmags=real(scoefs(n1,n2,n3)*conjg(scoefs(n1,n2,n3)))
            if (cmags > cmagsmax) cmagsmax=cmags
          enddo
        enddo
      enddo
!
! print values larger than cmagsmax*cprecision
      cmagsmax=cmagsmax*cprecision
      print *,' '
      print *,sname,'  dimension=',nn1,nn2,nn3
      do n3=1,nn3
        do n2=1,nn2
          do n1=1,nn1
            cmags=real(scoefs(n1,n2,n3)*conjg(scoefs(n1,n2,n3)))
            if (cmags > cmagsmax) then 
              print ('(a,3i3,2f20.13)'),'scoef ',n1,n2,n3,scoefs(n1,n2,n3)
            endif
          enddo
        enddo
      enddo
!      
      end subroutine sh_prt_allc 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_prt_c --- Print subset of coefs in table for testing
!                     
!
! !INTERFACE:
!    
      subroutine sh_prt_c (c,ni,n1,n2,n3,n4,caname,csubrt,id)

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in)   :: caname
      character(len=*), intent(in)   :: csubrt
      integer,  intent(in) :: n1,n2,n3,n4,ni,id
      complex(r8), intent(in) :: c(ni,n4,n2)

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  Print a (sub)set of coefs in a table with few significant figures
!  Used for testing, especially when truncation is very small (like T5)

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer, parameter :: nprt=21
      integer  :: iprt(nprt,2)
      integer  :: k,l,imax,nk
      integer  :: nskip,nsk3,ipr,n1dn3
      integer  :: lskip,lx,k1
      real(r8) :: a,amax,scale
!
!  this routine prints out 2-d arrays for testing purposes.
!  if n3 = 2, then every other value ignored (useful for printing only
!             real or imaginary components of complex arrays, in which
!             case a should be input as first element of real or imag.
!             component of array.
!  caname, csubrt, id for identifying field to be printed.
!
!  scale printing
!
    do nk=1,n4
      amax=0.
      do k=1,n1,n3
        do l=1,n2
          a=abs( real(c(k,nk,l), r8)) 
          if (a > amax) amax=a
          a=abs(aimag(c(k,nk,l))) 
          if (a > amax) amax=a
        enddo
      enddo
      if (amax > 0.) amax=log10(amax)
      imax=amax-3
      scale=10.**imax
!
!  select set of values to print
!
      nskip=(n1/n3-1)/nprt + 1
      nsk3=nskip*n3
      ipr=n1/nsk3
      n1dn3=n1/n3 
!
      print ('(5a,i1,a,i3,a,1p1e8.0,a,i2,a,i3,a,2i3)'),                &
        '  Test of ',caname,' called from ',csubrt,'  nk=',nk,' id=',id, &
        '  scale=',scale,'  every ',nskip,'th printed  dim1=',ni, &
        ' n1,n2=',n1,n2 
      print ('(a,21i6)'),'   k= ',(k,k=1,n1dn3,nskip)
!
      lskip=1
      if (n2.gt.30) lskip=2
      do lx=1,n2,lskip
        l=n2-lx+1
        k1=0
        do k=1,n1,nsk3
          k1=k1+1
          iprt(k1,1)=int( real(c(k,nk,l), r8)/scale)
          iprt(k1,2)=int(aimag(c(k,nk,l))/scale)
        enddo
        print ('(a,i3,a,21i6)'),' ',l,' R  ',(iprt(k,1),k=1,ipr)
        print ('(a,i3,a,21i6)'),' ',l,' I  ',(iprt(k,2),k=1,ipr)
        print *,' '
      enddo
!
    enddo  ! loop over nk
!
    end subroutine sh_prt_c

!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_prt_f --- Print portion of fields for testing software
!                     
!
! !INTERFACE:
!    
      subroutine sh_prt_f (f,ni,n1,n2,n3,cname,csubrt,id)

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in)   :: cname
      character(len=*), intent(in)   :: csubrt
      integer,  intent(in) :: n1,n2,n3,ni,id
      real(r8), intent(in) :: f(ni,n2)

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  this routine prints out 2-d arrays for testing purposes.
!  if n3 = 2, then every other value ignored (useful for printing only
!             real or imaginary components of complex arrays, in which
!             case a should be input as first element of real or imag.
!             component of array.
!  caname, csubrt, id for identifying field to be printed.
!
! Printed values are scaled, using few significant digits.
! Used when testing the code (e.g., to see if m=1 scoef creates
! a wavnumber 1 field)

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer, parameter :: nprt=21
      integer  :: iprt(nprt)
      integer  :: k,l,imax
      integer  :: nskip,nsk3,ipr,n1dn3
      integer  :: lskip,lx,k1
      real(r8) :: amax,scale
!
      amax=0.
      do k=1,n1,n3
        do l=1,n2
          if (abs(f(k,l)) .gt. amax) amax=abs(f(k,l))
        enddo
      enddo
      if (amax.gt.0.) amax=log10(amax)
      imax=amax-3
      scale=10.**imax
!
!  select set of values to print
!
      nskip=(n1/n3-1)/nprt + 1
      nsk3=nskip*n3
      ipr=n1/nsk3
      n1dn3=n1/n3 
!
      print ('(5a,i3,a,1p1e8.0,a,i2,a,i3,a,2i3)'),                &
        ' Test of ',cname,' called from ',csubrt,'  id=',id,     &
        '  scale=',scale,'  every ',nskip,'th printed  dim1=',ni, &
        ' n1,n2=',n1,n2 
      print ('(a,21i6)'),'   k= ',(k,k=1,n1dn3,nskip)
!
      lskip=1
      if (n2.gt.30) lskip=2
      do lx=1,n2,lskip
        l=n2-lx+1
        k1=0
        do k=1,n1,nsk3
          k1=k1+1
          iprt(k1)=f(k,l)/scale
        enddo
        print ('(a,i3,a,21i6)'),' ',l,'  ',(iprt(k),k=1,ipr)
      enddo
!
      end subroutine sh_prt_f
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_ptep_invert --- Invert the PT*P matrix to create wmatrix
!                     
!
! !INTERFACE:
!    
      subroutine sh_ptep_invert (ptep,nt)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)    :: nt

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

      real(r8), intent(inout) :: ptep(0:nt,0:nt) 

! !DESCRIPTION: 
! Invert a symmetric matrix by eigen-decomposition 
! The input matrix is over written.  ptep is the P^T * P matrix.

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      character(len=3), parameter :: purpose='VSI'
      real(r8)               :: weights(nt+1)
      real(r8), parameter    :: testvalue=1.d-8
      real(r8)               :: evectrs(0:nt,0:nt)
      real(r8)               :: evalues(0:nt)
      real(r8)               :: matrix(0:nt,0:nt)
      integer                :: nk

      if (nt == 0) then
        ptep(0,0)=one/ptep(0,0)
      else
        nk=nt+1
        weights=one
        matrix=ptep
        call Eigen_Main (evectrs, evalues, matrix, weights, nk, &
                         testvalue, purpose, 'ptep inverse', ptep)
      endif
!
      end subroutine sh_ptep_invert 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_reorderfm --- Routine for reordering array of half-transfom
!                     
!
! !INTERFACE:
!    
      subroutine sh_reorderfm (lots,jlats,mvals,ix,fkjm,fmjk,type)

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: type
      integer, intent(in)  :: lots,jlats,mvals,ix

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

      complex(r8), intent(inout) :: fmjk(ix,jlats,lots)
      complex(r8), intent(inout) :: fkjm(lots,jlats,0:mvals-1)

! !DESCRIPTION: 

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer     :: j,k,m

      do k=1,lots
        do j=1,jlats
 
          if (type == 'kjm2mjk' ) then

            do m=1,mvals
              fmjk(m,j,k)=fkjm(k,j,m-1)
            enddo
            if (ix > mvals) then
              do m=mvals+1,ix
                fmjk(m,j,k)=cmplx(zero,zero,r8)
              enddo
            endif

          else             
      
            do m=1,mvals
              fkjm(k,j,m-1)=fmjk(m,j,k)
            enddo

          endif 
        enddo  ! loop over j   
      enddo    ! loop over k
!
      end subroutine sh_reorderfm
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_sequence_f2s --- Sequence routines for field to spectra trans
!                     
!
! !INTERFACE:
!    
      subroutine sh_sequence_f2s (scoefs,f,name,shift,coslat,  &
                          sinlat,abfaclats,wmatrix,            &
                          nmax,mmax,jmax,nlevels,nspects,imax, &
                          jfirst,jlast,ntrunc,ik,index )
!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: name  ! name of field: 
                                            ! S (any scalar), V, or U
      integer,  intent(in)  :: nmax    ! maximum n for m=0
      integer,  intent(in)  :: mmax    ! maximum zonal wave number
      integer,  intent(in)  :: imax    ! number of longitudes
      integer,  intent(in)  :: jmax    ! 2nd dimension of field array
                                       ! = # lats for Scalar grid
      integer,  intent(in)  :: nlevels ! number of vertical levels
      integer,  intent(in)  :: jfirst  ! first latitude with data
      integer,  intent(in)  :: jlast   ! last  latitude with data
      integer,  intent(in)  :: nspects ! number of spectral coefs per level
      integer,  intent(in)  :: ik      ! =1 if scalar, 2 if vector field 
      integer,  intent(in)  :: ntrunc(0:mmax)       ! max n-m for each m
      integer,  intent(in)  :: index(0:mmax,0:nmax) ! order of m,n values
      real(r8), intent(in)  :: shift   ! offset of i=1 from Prime Meridian
      real(r8), intent(in)  :: sinlat(jmax)
      real(r8), intent(in)  :: coslat(jmax)
      real(r8), intent(in)  :: abfaclats(jmax,0:mmax) ! factors to compute alp
      real(r8), intent(in)  :: wmatrix(jmax,jmax,0:1) ! weighting matrix
      real(r8), intent(in)  :: f(imax,jmax,nlevels)   ! field to create

! !OUTPUT PARAMETERS:

      complex(r8), intent(out) :: scoefs(nspects,nlevels,ik) ! spectral coefs

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  This routine sequences all the computation to compute complex 
!  spherical harmonic spectral coefficients (scoefs) from fields on  
!  defined on a global D-grid as used by the GMAO FVGCM.
!
!  The input field can be:  a scalar field, such as T or q.
!                           the v component of a vector field 
!                           the u component of a vector field 

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer  :: jlats
      integer  :: mvals
      integer  :: nf
      integer  :: jnf
      integer  :: imax2p1
      real(r8) :: shift1     ! = -shift
      real(r8) :: vec1(jmax) ! vector whose values =1
      complex(r8) :: fmjk(imax/2+1,jfirst:jlast,nlevels)
      complex(r8) :: fkjm(nlevels,jfirst:jlast,0:mmax)
!
!
      jlats=jlast-jfirst+1
      mvals=mmax+1
      jnf=jlats*nlevels 
      imax2p1=imax/2 + 1
      vec1(:)=one
!
! Do zonal half-transform on all latitudes and fields, and reorder
      do nf=1,nlevels
        call sh_fft_slow (-1, imax, jlats, f(:,jfirst:jlast,nf), imax, &
                     fmjk(:,jfirst:jlast,nf), imax2p1)
      enddo
      call sh_reorderfm (nlevels,jlats,mvals,imax2p1,fkjm,fmjk,'mjk2kjm')
!
! If wind field, divide by coslat
     if (name /= 'S') then 
       call sh_uvcos (jfirst,jlast,nlevels,mmax,coslat(jfirst:jlast), &
                      fkjm,'R')
     endif
!
! Adjust phases of complex zonal Fourier coefficients if the first 
! longitude point (i=1) is not the Prime meridian.  For the FVDAS, 
! this concerns the grids for S and u, for which shift=1./(imax*2.).
      if (shift /= zero) then
        shift1=-shift
        call sh_shift (jnf, mmax, fkjm, shift1)
      endif
!
! Call Legendre transform with adjustment
     call sh_fm2s (scoefs,fkjm,wmatrix,name,vec1,sinlat,     &
                   abfaclats,nmax,mmax,jmax,nlevels,nspects, &
                   jfirst,jlast,ntrunc,index,ik )
!
     end subroutine sh_sequence_f2s
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_sequence_s2f --- Sequence routines for spectral to field trans
!                     
!
! !INTERFACE:
!    

      subroutine sh_sequence_s2f (scoefs,f,name,shift,       &
                          coslat,sinlat,abfaclats,           &
                          nmax,mmax,jmax,nlevels,nspects,    &
                          imax,jfirst,jlast,ntrunc,ik,index )


!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: name  ! name of field: S, V, or U
      integer,  intent(in)  :: nmax    ! maximum n for m=0
      integer,  intent(in)  :: mmax    ! maximum zonal wave number
      integer,  intent(in)  :: imax    ! number of longitudes
      integer,  intent(in)  :: jmax    ! 2nd dimension of field array
      integer,  intent(in)  :: nlevels ! number of vertical levels
      integer,  intent(in)  :: jfirst  ! first latitude with data
      integer,  intent(in)  :: jlast   ! last  latitude with data
      integer,  intent(in)  :: nspects ! number of spectral coefs per level
      integer,  intent(in)  :: ik      ! =1 if scalar, 2 if vector field 
      integer,  intent(in)  :: ntrunc(0:mmax)       ! max n-m for each m
      integer,  intent(in)  :: index(0:mmax,0:nmax) ! order of m,n values
      real(r8), intent(in)  :: shift   ! offset of i=1 from Prime Meridian
      real(r8), intent(in)  :: sinlat(jmax)
      real(r8), intent(in)  :: coslat(jmax)
      real(r8), intent(in)  :: abfaclats(jmax,0:mmax) ! factors to compute pnm
      real(r8), intent(in)  :: f(imax,jmax,nlevels)   ! field to create

! !OUTPUT PARAMETERS:

      complex(r8), intent(inout) :: scoefs(nspects,nlevels,ik) ! spectral coefs

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!
!  This routine sequences all the computation to compute real fields  
!  (f) from spherical harmonic spectral coefficients (scoefs). 
!  The fields are defined on lat/lon grid, with evenly spaced 
!  longitudes and latitudes specified by an array of sines of those 
!  latitudes.
!
!  The output field can be:  a scalar field, such as T or q.
!                            the v component of a vector field 
!                            the u component of a vector field 

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!  04Aug2008  R. Todling Bug fix: scoefs is inout
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
!
      integer  :: jlats
      integer  :: mvals
      integer  :: nf
      integer  :: jnf
      integer  :: imax2p1
      complex(r8) :: fmjk(imax/2+1,jfirst:jlast,nlevels)
      complex(r8) :: fkjm(nlevels,jfirst:jlast,0:mmax)
!
!
      jlats=jlast-jfirst+1
      mvals=mmax+1
      jnf=jlats*nlevels 
      imax2p1=imax/2 + 1
!
! Perform half transform from coefs to zonal waves on lats 
      call sh_s2fm (scoefs,fkjm,name,sinlat,abfaclats, &
                    nmax,mmax,jmax,nlevels,nspects,    &
                    jfirst,jlast,ntrunc,index,ik )
!
! If wind field, divide by coslat
     if (name /= 'S') then 
       call sh_uvcos (jfirst,jlast,nlevels,mmax,coslat(jfirst:jlast), &
                      fkjm,'R')
     endif
!
! Adjust phases of complex zonal Fourier coefficients if the first 
! longitude point (i=1) is not the Prime meridian.  For the FVDAS, 
! this concerns the grids for S and u, for which shift=1./(imax*2.).
      if (shift /= zero) then
        call sh_shift (jnf, mmax, fkjm, shift)
      endif
!
! Do zonal half-transform on all latitudes and fields, and reorder
     call sh_reorderfm (nlevels,jlats,mvals,imax2p1,fkjm,fmjk,'kjm2mjk')
     do nf=1,nlevels
       call sh_fft_slow (1, imax, jlats, f(1:imax,jfirst:jlast,nf), imax, &
                    fmjk(1:imax2p1,jfirst:jlast,nf), imax2p1)
     enddo

     end subroutine sh_sequence_s2f
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_setrig_D --- Set sines and cosines for D grid
!                     
!
! !INTERFACE:
!    
      subroutine sh_setrig_D (im,jm,dp,dl,cosp,cose,sinp,sine)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: im, jm  ! grid dimensions

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: sinp(jm) ! since of points on 'S' grid
      real(r8), intent(out) :: sine(jm) ! sines of points on 'U' grid
      real(r8), intent(out) :: cosp(jm) ! delta sine on 'S' grid
      real(r8), intent(out) :: cose(jm) ! delta sine on 'U' grid
      real(r8), intent(out) :: dp, dl 

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!
!  Copy of routine setrig from FVGCM that computes sines and 
!  cosines for D-grid.  Note that only sinp is computed as 
!  the exact sine at particular points (those where T, q, and delp 
!  are defined). The letter 'e' indicates 'edge' of 'S' grid cells
!
!  Note: the finite-difference forms of the cosines calculated 
!  here are not used for the transforms.
!  Instead, they are replaced by exact cosines of the sines.
!

! !SEE ALSO:  sh_coslat 
!

! !REVISION HISTORY:
!
!  ?????????  Copied and modified form FVGCM routine setrig
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer  :: j, jm1
      real(r8) :: pi, ph5

      jm1 = jm - 1
      pi  = 4.d0 * datan(1.d0)
      dl  = (pi+pi)/dble(im)   ! longitudinal spacing in radians
      dp  = pi/dble(jm1)       ! latitudinal spacing in radians

      do j=2,jm
        ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(jm1))
        sine(j) = dsin(ph5)
      enddo
      cosp( 1) =  0.
      cosp(jm) =  0.

      do j=2,jm1
        cosp(j) = (sine(j+1)-sine(j)) / dp
      enddo

! Define cosine at edges..

      do j=2,jm
        cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      enddo
      cose(1) = cose(2)

      sinp( 1) = -1.
      sinp(jm) =  1.

      do j=2,jm1
        sinp(j) = 0.5 * (sine(j) + sine(j+1))
      enddo

      end subroutine sh_setrig_D
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_shift --- change phase of wave to account for grid offset
!                     
!
! !INTERFACE:
!    
      subroutine sh_shift (lot, mmax, y, fraction)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: lot
      integer,  intent(in) :: mmax
      real(r8), intent(in) :: fraction

! !OUTPUT PARAMETERS:

      complex(r8), intent(inout) :: y(lot,0:mmax)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! Shift complex coefficients of zonal waves by a value to account for
! the first longitude not being the prime meridian. 

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer      :: m, nf
      complex(r8)  :: phase
      complex(r8)  :: cfac
      real(r8)     :: twopi

      twopi=(two**3)*atan(one)
      phase=cmplx(zero, fraction*twopi, r8)
      
      do m=0,mmax
        cfac=cexp(m*phase)
        do nf=1,lot
          y(nf,m)=y(nf,m)*cfac
        enddo
      enddo
! 
      end subroutine sh_shift
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_s2c_s --- Peel off scalar spectral coefs for one value of n-m
!                     
!
! !INTERFACE:
!    
      subroutine sh_s2c_s (scoefs,nspects,nmax,mmax,nlevels,m, &
                           nlast,index,cwork)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  ::  nmax
      integer,  intent(in)  ::  mmax
      integer,  intent(in)  ::  nlevels
      integer,  intent(in)  ::  nspects
      integer,  intent(in)  ::  m
      integer,  intent(in)  ::  nlast
      integer,  intent(in)  ::  index(0:mmax,0:nmax) 
      complex(r8), intent(in) :: scoefs(nspects,nlevels,1)

! !OUTPUT PARAMETERS:

      complex(r8), intent(out)  :: cwork(0:nlast,nlevels)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
!   
      integer :: n
      integer :: nf
      integer :: i

      do n=0,nlast
        i=index(m,n)
        do nf=1,nlevels           
          cwork(n,nf)=scoefs(i,nf,1)
        enddo
      enddo

      end subroutine sh_s2c_s
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_s2c_v --- Peel off vector spectral coefs for one value of n-m
!                     
!
! !INTERFACE:
!    

      subroutine sh_s2c_v (scoefs,nspects,nmax,mmax,nlevels,m, &
                           nlast,index,epsi,cwork,name)

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in)  :: name
      integer,  intent(in)  ::  nmax
      integer,  intent(in)  ::  mmax
      integer,  intent(in)  ::  nlevels
      integer,  intent(in)  ::  nspects
      integer,  intent(in)  ::  m
      integer,  intent(in)  ::  nlast
      integer,  intent(in)  ::  index(0:mmax,0:nmax) 
      real(r8), intent(in)  ::  epsi(0:nlast)
      complex(r8), intent(in)   :: scoefs(nspects,nlevels,2)

! !OUTPUT PARAMETERS:

      complex(r8), intent(out)  :: cwork(0:nlast,nlevels)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! Compute spectral coeffients of scalar u or v field, divided by 
! cos(theta) for one m from coefficientts for vorticity and divergence

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
!
      integer :: n
      integer :: nf
      integer :: i
      integer :: truen  ! the true index n
      integer :: id1, id2, is
      real(r8)    :: rfacm, rfac
      complex(r8) :: cfacm

      rfacm=-real(m,r8)*earthr
      if (name == 'U') then
        id1=2
        id2=1
        is=1
      else
        id1=1
        id2=2
        is=-1
      endif
!
      do n=0,nlast-1
        truen=n+m

        i=index(m,n)
        if (truen == 0) then
          cfacm=cmplx(zero,zero,r8)
        else
          cfacm=cmplx(zero,rfacm,r8)/((truen+1)*truen)
        endif
        do nf=1,nlevels           
          cwork(n,nf)=cfacm*scoefs(i,nf,id1)
        enddo
      
        if ( n > 0) then
          i=index(m,n-1)
          rfac=-is*epsi(n)*earthr/truen
          do nf=1,nlevels           
            cwork(n,nf)=cwork(n,nf)+rfac*scoefs(i,nf,id2)
          enddo
        endif  
      
        if ( n < nlast -1 ) then
          i=index(m,n+1)
          rfac=is*epsi(n+1)*earthr/(truen+1)
          do nf=1,nlevels           
            cwork(n,nf)=cwork(n,nf)+rfac*scoefs(i,nf,id2)
          enddo
        endif  

      enddo ! loop over n

      i=index(m,nlast-1)
      truen=nlast+m
      rfac=-is*epsi(nlast)*earthr/truen
      do nf=1,nlevels           
        cwork(nlast,nf)=rfac*scoefs(i,nf,id2)
      enddo

      end subroutine sh_s2c_v
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_s2fm --- Compute zonal wave coefs from spher. harmonic coefs
!                     
!
! !INTERFACE:
!    
      subroutine sh_s2fm (scoefs,fm,name,sinlat,abfaclats,  & 
                          nmax,mmax,jmax,nlevels,nspects,   &
                          jfirst,jlast,ntrunc,index,ik     )

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: name
      integer,  intent(in)  :: nmax, mmax, jmax, nlevels, nspects
      integer,  intent(in)  :: jfirst, jlast, ik 
      integer,  intent(in)  :: ntrunc(0:mmax), index(0:mmax,0:nmax) 
      real(r8), intent(in)  :: sinlat(jmax)
      real(r8), intent(in)  :: abfaclats(jmax,0:mmax)
      complex(r8), intent(in)   :: scoefs(nspects,nlevels,ik)

! !OUTPUT PARAMETERS:

      complex(r8), intent(out)  :: fm(nlevels,jfirst:jlast,0:mmax)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  Compute zonal wave coefficients given a set of coefficients for 
!  spherical harmonic functions. If values for a vector wind field are
!  to be created, then spectral coefficients for both vorticty and 
!  divergence are used. 
!
!  Note that if the grid is symmetric about the equator, only 1/2 the 
!  alp are needed (plus the equator if it is on the grid).  In this case,
!  jalp < jlast.

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables
!
      integer    :: m
      integer    :: j
      integer    :: nt   ! max of n-m given n
      integer    :: ntpx ! either nt (if scalar field) or nt+1 (if vector)
      integer    :: jalp ! last of the latitudes needed for alp values 
      real(r8),    allocatable :: alp(:,:)
      complex(r8), allocatable :: cwork(:,:)
!    
! Determine number of needed latitude values of alp
      if (symmetricX) then
        jalp=jfirst+(jlast-jfirst)/2
      else
        jalp=jlast
      endif
!
! Loop over zonal wave numbers m
      do m=0,mmax

        nt=ntrunc(m) 
        if ( name == 'S') then ! truncation for scalar field
          ntpx=nt
        else                   ! truncation for vector field
          ntpx=nt+1
        endif

! Compute alp for one m but all n and latitudes
        allocate (alp(0:ntpx,jfirst:jalp))
        do j=jfirst,jalp     
          call sh_alp_1 (alp(:,j),ntpx,m,j,sinlat(j),epsiX(:,m), &
                         abfaclats(j,m),alpsaveX)
        enddo

! If scalar fields, copy spectral coefficients for one m to an array
! If vector fields, compute spectral coefficients for one m from 
!    spectral coeficients of vorticity and divergence
        allocate (cwork(0:ntpx,nlevels))   
        if (name == 'S') then
          call sh_s2c_s (scoefs,nspects,nmax,mmax,nlevels,m, &
                             ntpx,index,cwork)
        else 
          call sh_s2c_v (scoefs,nspects,nmax,mmax,nlevels,m, &
                             ntpx,index,epsiX(:,m),cwork,name)
        endif
!
! Project polynomials onto latitudes for 1 m      
        call sh_c2fm (cwork,ntpx,jfirst,jlast,jalp,nlevels,alp,fm(:,:,m))

        deallocate (alp)
        deallocate (cwork)

      enddo ! loop over m

      end subroutine sh_s2fm
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_test_coefs --- Routine for testing this module
!                     
!
! !INTERFACE:
!    
      subroutine sh_test_coefs (imax,jmax,kmax,mmax,nmax,nlevels,name)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: imax
      integer,  intent(in) :: jmax
      integer,  intent(in) :: kmax
      integer,  intent(in) :: mmax
      integer,  intent(in) :: nmax
      integer,  intent(in) :: nlevels
      character(len=*),  intent(in) :: name

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! Routine for testing software by checking that transforms are exact
! All coefficients are set to some value. This also provides an example
! of use of this module.

! !SEE ALSO:  sh_fill for alternative way of setting coefs
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer  :: m,n,j  ! m=zonal wavenumber, n=n-m
      integer  :: ier
      integer  :: ifield,level
      integer  :: nspects
      integer  :: ntrunc(0:mmax)
      integer  :: mtrunc(0:nmax)
      integer  :: index(0:mmax,0:nmax)
     
      real(r8) :: ci
      real(r8) :: E_field(nlevels,4) ! area-mean of square of field on grid
      real(r8) :: E_coefs(nlevels,4) ! sum of power of coefs
      real(r8) :: E_diffs(nlevels,4) ! diff of field vs coef derived values
      real(r8) :: u(imax,jmax,nlevels)
      real(r8) :: v(imax,jmax,nlevels)
      real(r8) :: T(imax,jmax,nlevels)
      real(r8) :: power(0:kmax,nlevels,3) 
      complex(r8), allocatable :: scoefsS(:,:,:) 
      complex(r8), allocatable :: scoefsV(:,:,:) 
!
! Set arrays describing spectral truncation, ordering, and arrays required 
! to compute alps (the latter remain internal to the module).  alps are the
! associated Legendre polynomials. The first occurence of .true. in the 
! 2 argument lists below indicates that alp is to be computed once, and 
! then saved for use in the actual transforms to be called later.
!
      if (name == 'D') then
        call sh_initial_D (nmax,mmax,kmax,imax,jmax,nspects,   &
                           ntrunc,mtrunc,index,'both',.true.,  &
                           .true.,ier)
        print *,' sh_initial_D returns ierror=',ier,'  (0=OK)'
     else if (name == 'G') then
        call sh_initial_G (nmax,mmax,kmax,imax,jmax,nspects,   &
                           ntrunc,mtrunc,index,'both',.true.,  &
                           .true.,ier)
        print *,' sh_initial_G returns ierror=',ier,'  (0=OK)'
     else
        print *,' X X X X X X X X X X X X X X X X X X X X X X X X '
        print *,'  Invalid name submitted to sh_test_coefs:'
        print *,'  name=',name,'    only D or G permitted'
        stop
     endif
!
! Now that nspects is known for this truncation, allocate spectral coefs
      allocate (scoefsS(nspects,nlevels,1))
      allocate (scoefsV(nspects,nlevels,2))
!
! Set some scoefs for testing
      do m=0,mmax
        if (m == 0) then
          ci=zero          ! the imag. part for m=0 is zero
        else
          ci=44.d0
        endif 
        do n=0,ntrunc(m) 
          j=index(m,n)
          scoefsS(j,:,1)=cmplx(100.d0,ci,r8)
          scoefsV(j,:,1)=cmplx(100.d0,ci,r8)
          scoefsV(j,:,2)=cmplx( 50.d0,ci,r8)
        enddo
      enddo
      scoefsV(1,:,:)=cmplx(zero,zero,r8)  ! vort and divg have 0 global mean

!QQQQQQ
      scoefsS(:,:,:)=cmplx(zero,zero)
      scoefsS(1,:,:)=cmplx(one,zero)
      scoefsS(2,:,:)=cmplx(two,zero)
      scoefsS(4,:,:)=cmplx(3.d0,3.d0)

!
! Print out selected values of initial coefficients
      call sh_prt_c (scoefsS,nspects,nspects,1,1,nlevels,&
          'scoefs','after fill',1)
      call sh_prt_c (scoefsV,nspects,nspects,1,2,nlevels,&
          'scoefs','after fill',2)
      call sh_prt_allc (scoefsS,nspects,nlevels,1,'scoefsS start ')
      call sh_prt_allc (scoefsV,nspects,nlevels,2,'scoefsV start ')
!
! Initialize fields to zero since the projections of the spherical 
! harmonic coefficients will be added to the existing fields.
      u=zero
      v=zero
      T=zero
!
! Transform from spectral coefs to fields
      call sh_transforms (imax,jmax,mmax,nmax,kmax,nlevels, &
                                nspects,ntrunc,index,       &
                                1,T,scoefsS,'S','s2f',ier)
      if (ier /= 0) print *,' ier=',ier,' returned from transform 1 '

      call sh_transforms (imax,jmax,mmax,nmax,kmax,nlevels, &
                                nspects,ntrunc,index,       &
                                2,u,scoefsV,'U','s2f',ier)
      if (ier /= 0) print *,' ier=',ier,' returned from transform 2 '

      call sh_transforms (imax,jmax,mmax,nmax,kmax,nlevels, &
                                nspects,ntrunc,index,       &
                                2,v,scoefsV,'V','s2f',ier)
      if (ier /= 0) print *,' ier=',ier,' returned from transform 3 '
!
! Compute area mean of square of field on grid
      E_field(:,1:3)=zero 
      call sh_intsqr_grid (imax,jmax,nlevels,T,E_field(:,1),'S')
      call sh_intsqr_grid (imax,jmax,nlevels,u,E_field(:,2),'U')
      call sh_intsqr_grid (imax,jmax,nlevels,v,E_field(:,3),'V')
      E_field(:,4)=E_field(:,2)+E_field(:,3)

!QQQQQ     
      print *,' ' 
      print *,' FIELD=S'
      do n=jmax,1,-1
        print ('(i2,1p8e10.2)'), n,T(:,n,1)
      enddo
      print *,' ' 

!
! Re-set scoefs to 0 since transform will increment them otherwise
      scoefsS=cmplx(zero, zero, r8)
      scoefsV=cmplx(zero, zero, r8)
!
! Transform from fields to spectral coefficients
      call sh_transforms (imax,jmax,mmax,nmax,kmax,nlevels, &
                                nspects,ntrunc,index,       &
                                1,T,scoefsS,'S','f2s',ier)
      if (ier /= 0) print *,' ier=',ier,' returned from transform 4 '

     call sh_transforms (imax,jmax,mmax,nmax,kmax,nlevels, &
                                nspects,ntrunc,index,       &
                                2,u,scoefsV,'U','f2s',ier)
      if (ier /= 0) print *,' ier=',ier,' returned from transform 5 '

      call sh_transforms (imax,jmax,mmax,nmax,kmax,nlevels, &
                                nspects,ntrunc,index,       &
                                2,v,scoefsV,'V','f2s',ier)
      if (ier /= 0) print *,' ier=',ier,' returned from transform 6 '
!
! Print re-constituted spectral coefficients
      call sh_prt_allc (scoefsS,nspects,nlevels,1,'scoefsS end ')
      call sh_prt_allc (scoefsV,nspects,nlevels,2,'scoefsV end ')
!
! Compute power spectra from complex coefficients
      power(:,:,:)=zero
      call sh_calc_power (mmax,nmax,kmax,imax,nspects,nlevels,1,    &
                          mtrunc,index,power(:,:,1:1),scoefsS,'S')
      call sh_calc_power (mmax,nmax,kmax,imax,nspects,nlevels,2,    &
                          mtrunc,index,power(:,:,2:3),scoefsV,'VV')
!
! sum power
      E_coefs(:,1:3)=zero       
      do ifield=1,3
        do level=1,nlevels
          do n=0,kmax 
            E_coefs(level,ifield)=E_coefs(level,ifield)+ &
                                  power(n,level,ifield) 
          enddo      
        enddo
      enddo
      E_coefs(:,4)=E_coefs(:,2)+E_coefs(:,3)       
!
! print mean squared values computed from the fields vs. from the coefs.
      E_diffs(:,:)=E_coefs(:,:)-E_field(:,:)       
      print *,' '
      print *,' Mean squared values of scalar fields and total power'
      print *,' level, scalar coefs, scalar fields, diff; ' &
             ,' vector coefs, vector fields, diff '
      print *,' diff should be negligible for scalar fields'
      print *,' diff should be negligible for vector fields if G-grid'
      print *,' diff should be small for vector fields if D-grid'
      do level=1,nlevels
        print ('(i3,1p3e12.2,4x,3e12.2)'),level, &
                E_coefs(level,1),E_field(level,1),E_diffs(level,1), &
                E_coefs(level,4),E_field(level,4),E_diffs(level,4) 
      enddo
!
!  Rescale power for easier (nicer) values to print for this test.
!  This is particularly necessary for the wind-related values 
!  because the test values used to initialize the coefficients are 
!  not geophysical values. 
      power(:,:,1:1)=power(:,:,1:1)/100.d0**2
      power(:,:,2:3)=power(:,:,2:3)/(100.d0*earthr)**2
      print *,' '
      print *,' Power spectra follow '
      do ifield=1,3
        do level=1,nlevels
          print *,'  ifield = ',ifield,'  level = ',level
          print ('(2x,16f8.2)'),power(0:kmax,level,ifield)
        enddo
      enddo
!
! call clean to deallocate saved module internal arrays
      call sh_clean
      deallocate (scoefsS,scoefsV)
!      
      end subroutine sh_test_coefs 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_transforms --- Driving routine for either-way transforms
!                     
!
! !INTERFACE:
!    
      subroutine sh_transforms (imax,jmax,mmax,nmax,kmax,nlevels, &
                                nspects,ntrunc,index,             &
                                isp,f,scoefs,name,func,ierror )

!USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*), intent(in) :: func  ! instruction 'f2s' or 's2f'
      character(len=*), intent(in) :: name  ! name of field: 'U','V','S'

      integer,  intent(in)  :: imax   ! number of longitudes
      integer,  intent(in)  :: jmax   ! 2nd dimension of field array
                                      ! = # lats for Scalar grid
      integer, intent(in)  :: mmax    ! maximum zonal wave number m
      integer, intent(in)  :: nmax    ! maximum n for m=0
      integer, intent(in)  :: kmax    ! maximum n-m for any m
      integer, intent(in)  :: nlevels ! number of horz surf to transform
      integer, intent(in)  :: isp     ! 1 for scalar, 2 for vector coefs
      integer, intent(in) :: nspects         ! number of spectral coefs
      integer, intent(in) :: ntrunc(0:mmax)  ! max n-m for given m
      integer, intent(in) :: index(0:mmax,0:nmax) ! ordering of coefs 

! !OUTPUT PARAMETERS:

      integer,     intent(out)   :: ierror ! returned error (0=OK)

! !INPUT/OUTPUT PARAMETERS:

      real(r8),    intent(inout) :: f(imax,jmax,nlevels)        ! fields
      complex(r8), intent(inout) :: scoefs(nspects,nlevels,isp) ! coefs

! !DESCRIPTION: 
!  Note that fields or coefficients are added to here, based on whatever
!  their input values are.  A projection rather than incrementation therefore
!  requires they be set to zero before entering the routine.
!
!  Fields on multiple horizontal surfaces are processed simultaneously 
!  (variable nlevels).
!
!  Spectral transforms are in terms of spherical harmonic functions.
!  For the vector wind field, these coefficients are for the corresponding
!  vorticity and divergence fields.  In the latter case, this routine
!  should be called twice: once with the u-component of the wind as the field
!  and another time with the v-component. (u northward, v eastward). 
!
!  For transforms concerning scalar fields, isp=1. 
!  For transforms concerning components of vector fields, isp=2 since 
!  spectral coeeficients for both vorticity and divergence are computed.
!
!  ierror should be zero on output
!
!  It is assumed that some variables have already been tested in the 
!  routine sh_dgrid_init for their mutual compatability. Also, values for 
!  nspects, ntrunc, and index are determined in that routine.

! !SEE ALSO: 
!     sh_test_coefs for an example of use
!     sh_dgrid_init for an example of coefficient ordering index.

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer :: ig      ! indicates which subset of the D-grid is used

      ierror=0
!
! Determine which part of the should be used for this field; i.e.,
! if different fields are defined at different points.  Also, determine
! on which sub-grid any saved alp functions are defined for this field.
! Also check that the input name is one of those expected.
      if     ( grid_nameX == 'G' ) then   ! grid is Gaussian
        ig=1     ! same for all fields on a Gaussian grid   
        igalpX=1 ! alpX defined for only one grid    
      elseif ( grid_nameX == 'D' ) then   ! grid is FVGCM D
        if     (name == 'S') then
          ig=3
          igalpX=2 ! alpX defined for this field on sub-grid 2
        elseif (name == 'U') then
          ig=1
          igalpX=1 ! alpX defined for this field on sub-grid 1
        elseif (name == 'V') then
          ig=2
          igalpX=2 ! alpX defined for this field on sub-grid 2
        else    ! field name invalid
          ierror=ierror+1
        endif
      else      ! grid_nameX invalid
        ierror=ierror+2
      endif
!
      if (ierror.ne.0) then
        print *,' '
        print *,' X X X X X X X X X X X X X X X X X X X X X X X X X'
        print *,' sh_transforms called without valid variable name'
        print *,' name=',name,'  only U, V, or S permitted'
        print *,' grid_nameX=',grid_nameX,'   only G or D permitted'    
      endif
!
! If vector vields are to be transformed, check that spectral coefficient
! array has 3rd index that includes space for coefficients of both 
! vorticity and divergence.
      if ( (name /= 'S') .and. (isp == 1) ) then
        ierror=ierror+10
      endif
!
! Transform from spectra to fields or the reverse, as specified by func
      if (func == 's2f') then
        call sh_sequence_s2f (scoefs,f,name,shiftX(ig),              &
                         cospX(:,ig),sinpX(:,ig),abfaclatsX(:,:,ig), &
                         nmax,mmax,jmax,nlevels,nspects,             &
                         imax,jfirstX(ig),jlastX(ig),ntrunc,         &
                         isp,index)
!
      elseif (func == 'f2s') then
        call sh_sequence_f2s (scoefs,f,name,shiftX(ig),cospX(:,ig), &
                         sinpX(:,ig),abfaclatsX(:,:,ig),            &
                         wmatrixX(:,:,0:1,ig),                      &
                         nmax,mmax,jmax,nlevels,nspects,imax,       &
                         jfirstX(ig),jlastX(ig),ntrunc,isp,index )
      else
        ierror=ierror+100
        print *,' '
        print *,' X X X X X X X X X X X X X X X X X X X X X X X X X X'
        print *,' sh_transforms called without valid variable func'
        print *,' func=',func   
      endif
!
      end subroutine sh_transforms 
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_trunc --- Set some arrays that describe the triangular truncation
!                     
!
! !INTERFACE:
!    
      subroutine sh_trunc (mtrunc,ntrunc,index,mmax,nmax,kmax,icoft)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)  :: mmax 
      integer, intent(in)  :: nmax 
      integer, intent(in)  :: kmax 
      integer, intent(in)  :: icoft

! !OUTPUT PARAMETERS:

      integer, intent(out) :: mtrunc(0:nmax)
      integer, intent(out) :: ntrunc(0:mmax)
      integer, intent(out) :: index(0:mmax,0:nmax)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!  this routine calculates the two vector pointers mtrunc and ntrunc
!     mtrunc is the length of a diagonal in n,m space for given n-m+1
!     ntrunc is the length of a column   in n,m space for given m+1
!     icoft=0 if ordering is along columns (all n for each m)
!     icoft=1 if ordering is along diagonals (all m for each n)
!  values consistent with general pentagonal spectral truncation

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer :: j, ksubm, ksubn, m, n

      ksubm=kmax-mmax
      ksubn=kmax-nmax

      do n=0,nmax
        mtrunc(n)=mmax
        if (n > ksubm) mtrunc(n)=mtrunc(n-1)-1
      enddo

      do m=0,mmax
        ntrunc(m)=nmax
        if (m > ksubn) ntrunc(m)=ntrunc(m-1)-1
      enddo
!
!  specify ordering of modes
!     icoft=0 is ordering along vertical lines of m,n diagram
!     icoft=1 is ordering along diagonal lines of m,n diagram
!
      if (icoft == 0) then
        j=0
        do m=0,mmax
          do n=0,ntrunc(m)
            j=j+1
            index(m,n)=j
          enddo
        enddo
      elseif (icoft == 1) then  
        j=0
        do n=0,nmax
          do m=1,mtrunc(n)
            j=j+1
            index(m,n)=j
          enddo
        enddo
      endif
!
      end subroutine sh_trunc
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_uvcos --- multiply or divide field at 1 latitude by its cosine
!                     
!
! !INTERFACE:
!    
      subroutine sh_uvcos (jfirst,jlast,lot,mmax,coslat,fm,type)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in) :: jfirst,jlast,lot,mmax
      character(len=*),  intent(in) :: type
      real(r8),    intent(in) :: coslat(jfirst:jlast)

! !OUTPUT PARAMETERS:

      complex(r8), intent(inout) :: fm(lot,jfirst:jlast,0:mmax)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
! Either multiply or divide a field by the cosine of latitude

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      integer     :: j, m, k
      real(r8)    :: fac

      do j=jfirst,jlast

        if (type == 'R') then
          fac=one/coslat(j)
        else
          fac=coslat(j)
        endif 
    
        do m=0,mmax
          do k=1,lot
            fm(k,j,m)=fm(k,j,m)*fac
          enddo
        enddo

      enddo   
!
      end subroutine sh_uvcos
!    
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sh_wmatrix --- Compute weighting matrix for projections onto alp
!                     
!
! !INTERFACE:
!    
      subroutine sh_wmatrix (wmatrix,sinlat,abfaclats,mfirst,mlast, &
                             mmax,jmax,jfirst,jlast) 

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: mmax,jmax,jfirst,jlast
      integer,  intent(in)  :: mfirst,mlast
      real(r8), intent(in)  :: sinlat(jmax)
      real(r8), intent(in)  :: abfaclats(jmax,0:mmax)

! !OUTPUT PARAMETERS:

      real(r8), intent(out) :: wmatrix(jmax,jmax,mfirst:mlast)

! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION: 
!
! This routine computes the weighting matrix W defined in Eq. 2.7
! in Swartzrauber and Spotz (2000).  It is produced as the pseudo-
! inverse of the matrix P * P(transpose), where P is the matrix whose rows 
! are values of successive associated Legendre polynomials for selected 
! zonal wave numbers m at a single latitude, and rows are for 
! successive latitudes.  Generally, this matrix needs to be calculated
! for only m=0 and m=1. The matrix product is computed by using the
! simple expression for each element of the product matrix provided by
! the Christoffel-Darboux formula (see Swartzrauber and Spotz (2000)). 
!
! WARNING: This code may be appropriate only for m=0 1nd m=1

! !SEE ALSO: 
!

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!
!EOP
!-------------------------------------------------------------------------
!
!    local variables

      logical  :: pole
      integer  :: nm0, nm1
      real(r8), allocatable :: epsi(:,:)
      real(r8), allocatable :: alp(:,:)
      real(r8), allocatable :: wtemp(:,:)
      integer               :: j1, j2, jpoints, jpmax
      integer               :: i, j, m, n
      integer               :: n0, n1, ntpx
      real(r8)              :: jfac
      real(r8)              :: dalpn0, dalpn1
      real(r8)              :: cos2, cospole
!
      jpmax=jlast-jfirst+3
      allocate (epsi(0:jpmax+1,0:mlast))
      call sh_emns (epsi,jpmax,mlast)
!
!  Determine whether grid includes pole points
      cos2=one-sinlat(jfirst)*sinlat(jfirst) ! cos**2 at first grid point
      cospole=(one/(10*jmax))**2             ! value of cos**2 "near" pole
      if (cos2 < cospole ) then  ! first point "at" pole
        pole=.true.                  
      else                       ! first point not "near" pole   
        pole=.false.                  
      endif
!
!  Produce matrix only for selected zonal wave numbers m
!  Note that for m=1, the matrix would be singular if pole points
!  are present.  Therefore they are excluded and a pseudo inverse
!  is computed instead.
!
      do m=mfirst,mlast
!
        j1=jfirst
        j2=jlast 
        if (pole.and.(m==1)) then ! skip pole points
          j1=jfirst+1
          j2=jlast-1
        endif
        jpoints=j2-j1+1  
        nm0=jpoints
        nm1=nm0-1 
        n0=nm0+m
        ntpx=nm0+1        
        jfac=sqrt((n0**2-m*m)/((two*n0)**2-one))
!
! Compute required Legendre polynomials for n-m=0:jmax-1
        allocate (alp(0:ntpx,jfirst:jlast))
        do j=jfirst,jlast
          call sh_alp_1 (alp(:,j),ntpx,m,j,sinlat(j),epsi(:,m), &
                         abfaclats(j,m),.false.)
        enddo
! 
! Application of the Christoffel-Darboux formula:
! (note that for pole points, the elements are computed as they would be
! for an explicit matrix multiply rather than by use of the formula)
!
        allocate (wtemp(j1:j2,j1:j2))
        do i=j1,j2
          do j=i,j2
            if (i==j) then
              if (pole.and.(m==0).and.((i==j1).or.(i==j2))) then
                wtemp(i,i)=zero 
                do n=0,nm1
                  wtemp(i,i)=wtemp(i,i)+alp(n,i)*alp(n,i)
                enddo
              else
                dalpn0=(n0+one)*epsi(nm0,  m)*alp(nm0-1,i) -         &
                       (n0    )*epsi(nm0+1,m)*alp(nm0+1,i)
                dalpn1=(n0    )*epsi(nm1,  m)*alp(nm1-1,i) -         &
                       (n0-one)*epsi(nm1+1,m)*alp(nm1+1,i)
                wtemp(i,j)=jfac*                                     &
                       (dalpn0*alp(nm1,i)-dalpn1*alp(nm0,i))/        &
                       (one-sinlat(i)*sinlat(i))
              endif   
            else
              wtemp(i,j)=jfac*                                       &
                     (alp(nm0,i)*alp(nm1,j)-alp(nm1,i)*alp(nm0,j))/  &
                     (sinlat(i)-sinlat(j))
            endif
          enddo ! loop over j >= i
          if (i.gt.j1) then 
            do j=j1,i-1
              wtemp(i,j)=wtemp(j,i) ! since matrix is symmetric
            enddo
          endif   
        enddo  
        deallocate (alp)
!
! Invert matrix
        call sh_ptep_invert (wtemp,j2-j1)
!
! Construct entire W matrix from pseudo inverse of P*P(transpose)
        wmatrix(:,:,m)=zero
        wmatrix(j1:j2,j1:j2,m)=wtemp(j1:j2,j1:j2)
        if (j1>jfirst) then
          do j=jfirst,j1-1
            wmatrix (j,j,m)=one
          enddo 
        endif
        if (j2<jlast) then
          do j=j2+1,jlast
            wmatrix(j,j,m)=one
          enddo 
        endif
        deallocate (wtemp)

      enddo  ! loop over m
      deallocate (epsi)  
!
      end subroutine sh_wmatrix
!
!
      end module m_shtrans_DG
