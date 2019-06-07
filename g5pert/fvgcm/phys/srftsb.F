#include <params.h>
      subroutine srftsb(isrfty  ,indx    ,npts    ,fnt     ,dfntdt  ,
     $                  snowh   ,tsbsf   )
C-----------------------------------------------------------------------
C Compute surface and subsurface temperatures over sea-ice surfaces.
C
C Sea ice temperatures are specified in 'plevmx' layers of fixed
C thickness and thermal properties.  The forecast temperatures are
C determined from a backward/implicit diffusion calculation using
C linearized sensible/latent heat fluxes. The bottom ocean temperature
C is fixed at -2C, allowing heat flux exchange with underlying ocean.
C 
C Sub-surface layers are indexed 1 at the surface, increasing downwards
C to plevmx.  Layers have mid-points and interfaces between layers.
C
C Temperatures are defined at mid-points, while fluxes between layers
C and the top/bottom media are defined at layer interfaces.
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb, CCM2
C Standardized:      L. Buja, Jun 1992, Feb 1996
C Reviewed:          B. Briegleb, J. Hack, Aug 1992, March 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comtim.h>
C-----------------------------------------------------------------------
#include <comtsc.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer isrfty(plond)     ! Surface type index (1..7)
      integer indx(plond)       ! Index of points to be computed
      integer npts              ! Number of points to be computed

      real fnt(plond)           ! Top surface/atmosphere net energy flux
      real dfntdt(plond)        ! ts partial derivative of net sfc flux
      real snowh(plond)         ! Snow depth (liquid water equivalent)
C 
C Input/output arguments 
C 
      real tsbsf(plond,plevmx)  ! Initial/final surface/sub-surface tmps
C
C---------------------------Local parameters----------------------------
C
      integer    psrfty         ! Number of surface types
      parameter (psrfty=7)     
C
C---------------------------Local variables-----------------------------
C
      integer i,ii              ! Longitude indices
      integer ier               ! Error return on tridiagonal solver
      integer jndx              ! Surface type index
      integer n                 ! Sub-surface layer index
C 
      real cmass(plond,plevmx)  ! Specific heat of soil (J/kg/K)
      real rho(plond,plevmx)    ! Mass densty of sub-sfc mat (kg/m3)
      real tk(plond,plevmx)     ! Thermal conductivity (watts/m/K)
      real z(plond,0:plevmx)    ! Interface geometrical depth (m)
C 
      real diag(plond,plevmx)   ! Diagonal matrix elements
      real htsrc(plond,plevmx)  ! External heat source (W/m3)
      real rhs(plond,plevmx)    ! Rhs of tri-diagonal matrix equation
      real sbdiag(plond,plevmx) ! Sub-diagonal matrix elements
      real spdiag(plond,plevmx) ! Super-diagonal matrix elements
      real tin(plond,plevmx)    ! Initial sub-surface temperatures
      real ws(plond*plevmx)     ! Working storage for mtdlss
C
      real cmsnow               ! Snow mass heat capacity
      real cmty(plond)          ! Layer mass heat capacity
      real crt                  ! cmass*rho*rdtime
      real delz                 ! Layer thickness
      real delzmn               ! Thick from mid-point to lyr above mid-point
      real delzpl               ! Thick from mid-point to lyr below mid-point
      real fbt(plond)           ! Ocean heat flux into sea-ice
      real fmns                 ! 1/(delz*delzmn)
      real fpls                 ! 1/(delz*delzpl)
      real msnow                ! Mass path of snow
      real msoil                ! Mass path of soil
      real rdtime               ! Inverse model time step
      real rhsnow               ! Snow mass density
      real rhty(plond)          ! Layer mass density
      real rztop                ! 1/ztop
      real thck(plond)          ! Layer thickness
      real tkbot                ! Bottom layer top interf thermal conduct
      real tkmns                ! Layer bottom interface thermal conduct
      real tkpls                ! Layer top interface thermal conductivity
      real tksnow               ! Snow thermal conducitivity
      real tktop                ! Top layer bottom interface thermal conduct
      real tkty(plond)          ! Layer thermal conductivity
      real tmp                  ! crt - dfntdt(i)*rztop
      real zbot                 ! Bottom layer thickness
      real zm                   ! Present layer mid-point depth
      real zmmn                 ! Layer above mid-point depth
      real zmpl                 ! Layer below mid-point depth
      real zsnow                ! Snow geometric depth
      real ztop                 ! Top layer thickness
C
      logical scvr(plond)       ! True if surface snow covered
C
C--------------------------Data Statements------------------------------
C
C Specified (and invariant) thermal properties for surface types
C
      real cmtype(psrfty,plevmx)   ! Mass specific heat (J/kg/K)
      real rhtype(psrfty,plevmx)   ! Mass density (kg/m3)
      real thckly(psrfty,plevmx)   ! Layer thicknesses (m)
      real tktype(psrfty,plevmx)   ! Thermal conductivity (J/m/s)
      save cmtype,rhtype,thckly,tktype
C
      data cmtype /4.20e3,2.07e3,2.07e3,1.04e3,7.20e2,5.60e2,4.16e2,
     $             4.20e3,2.07e3,2.07e3,1.04e3,7.20e2,5.60e2,4.16e2,
     $             4.20e3,2.07e3,2.07e3,1.04e3,7.20e2,5.60e2,4.16e2,
     $             4.20e3,2.07e3,2.07e3,1.04e3,7.20e2,5.60e2,4.16e2/
C
      data rhtype /1.00e3,9.20e2,9.20e2,2.50e3,2.50e3,2.50e3,2.50e3,
     $             1.00e3,9.20e2,9.20e2,2.50e3,2.50e3,2.50e3,2.50e3,
     $             1.00e3,9.20e2,9.20e2,2.50e3,2.50e3,2.50e3,2.50e3,
     $             1.00e3,9.20e2,9.20e2,2.50e3,2.50e3,2.50e3,2.50e3/
C
      data thckly / 2., .500, .250, .050, .090, .080, .120,
     $              5., .500, .500, .366, .390, .435, .492,
     $             10., .500, .500,1.369,1.459,1.628,1.841,
     $             33.,.500,8.500,6.990,7.450,8.310,9.400/
C
      data tktype /15.0 ,2.200 ,2.200 ,1.408 ,1.104 ,1.071 ,1.019 ,
     $             15.0 ,2.200 ,2.200 ,1.408 ,1.104 ,1.071 ,1.019 ,
     $             15.0 ,2.200 ,2.200 ,1.408 ,1.104 ,1.071 ,1.019 ,
     $             15.0 ,2.200 ,2.200 ,1.408 ,1.104 ,1.071 ,1.019 /
C
C Properties of ice and air used to determine snow thermal properties
C
      real cmair,     ! Mass specific heat of air
     $     cmice,     ! Mass specific heat of ice
     $     frcair,    ! Fraction of air assumed in mix of ice and air
     $     rhair,     ! Mass density of surface air
     $     rhice,     ! Mass density of ice
     $     tkair,     ! Thermal conductivity of air
     $     tkice      ! Thermal conductivity of ice
      save cmair,cmice,frcair,rhair,rhice,tkair,tkice
C
      data cmair  /1.00e3/
      data cmice  /2.07e3/
      data frcair /0.90/
      data rhair  /1.25 /
      data rhice  /9.20e2/
      data tkair  /0.025/
      data tkice  /2.200/
C
C-----------------------------------------------------------------------
C
      rdtime = 1./dtime
C
C Calculate snow properties
C
      cmsnow = (1.-frcair)*cmice + frcair*cmair
      rhsnow = (1.-frcair)*rhice + frcair*rhair
      tksnow = (1.-frcair)*tkice + frcair*tkair
C
C No external heat source
C
      do n=1,plevmx
         do ii=1,npts
            i = indx(ii)
            htsrc(i,n) = 0.0
         end do
      end do
C
C Define logical for snow covered surfaces:
C
      do ii=1,npts
         i = indx(ii)
         scvr(i) = snowh(i).gt.0.
      end do
C
C Define thermal properities for each sub/surface layer, starting
C with the top layer
C
      do ii=1,npts
         i = indx(ii)
         jndx    = isrfty(i)
         thck(i) = thckly(jndx,1)
         cmty(i) = cmtype(jndx,1)
         rhty(i) = rhtype(jndx,1)
         tkty(i) = tktype(jndx,1)
      end do
C
CDIR$ IVDEP
      do ii=1,npts
         i = indx(ii)
C
C Initialize fields for no snow cover
C
         z(i,0)     = 0.0
         z(i,1)     = thck(i)
         cmass(i,1) = cmty(i)
         rho(i,1)   = rhty(i)
         tk(i,1)    = tkty(i)
C
C Modify layer 1 fields for snow cover if present
C
         if(scvr(i)) then
C
C Snow equivlnt depth times snow liquid water depth gives the physical
C depth of snow for thermal conduction computation; snow is mixed
C uniformly by mass with the top surface layer
C
           zsnow     = snowh(i)*snwedp
           msnow     = rhsnow*zsnow
           msoil     = rhty(i)*thck(i)
           rho(i,1)  = (msnow*rhsnow + msoil*rhty(i))/(msnow+msoil)
           cmass(i,1)= (msnow*cmsnow + msoil*cmty(i))/(msnow+msoil)
           tk(i,1)   = (msnow*tksnow + msoil*tkty(i))/(msnow+msoil)
           z(i,1)    = (msnow+msoil) / rho(i,1)
         end if
C
      end do
C
C Set surface thermal properties for the lower sub/surface layers:
C
      do n=2,plevmx
         do ii=1,npts
            i = indx(ii)
            jndx       = isrfty(i)
            thck(i)    = thckly(jndx,n)
            cmass(i,n) = cmtype(jndx,n)
            rho(i,n)   = rhtype(jndx,n)
            tk(i,n)    = tktype(jndx,n)
         end do
         do ii=1,npts
            i = indx(ii)
            z(i,n)  = z(i,n-1) + thck(i)
         end do
      end do
C
C Define set of linear equations for temperature
C
      do n=1,plevmx
         do ii=1,npts
            i = indx(ii)
            tin(i,n) = tsbsf(i,n)
         end do
      end do
C
C If sea ice, compute heat flux from underlying ocean, assumed to be at
C the temperature of -2C
C
CDIR$ IVDEP
      do ii=1,npts
         i = indx(ii)
         fbt(i) = 0.0
         if(isrfty(i).eq.2) then
            zbot    = 0.5*(z(i,plevmx) - z(i,plevmx-1))
            fbt(i)  = -tk(i,plevmx)*(271.16 - tin(i,plevmx))/zbot
         end if
      end do
C
C Set up linear equations
C
      do ii=1,npts
         i = indx(ii)
         sbdiag(i,1)      = 0.
         spdiag(i,plevmx) = 0.
      end do
      if (plevmx.eq.1) then
C
C Single layer
C
         do ii=1,npts
            i = indx(ii)
            rztop     = 1./(z(i,1) - z(i,0))
            crt       = (cmass(i,1)*rho(i,1)*rdtime)
            diag(i,1) = crt - dfntdt(i)*rztop
            rhs(i,1)  = diag(i,1)*tin(i,1) + fnt(i)*rztop -
     $                  fbt(i)*rztop + htsrc(i,1)
         end do
      else if (plevmx.gt.1) then
C
C More than one layer: top layer first
C
         do ii=1,npts
            i = indx(ii)
            crt         = cmass(i,1)*rho(i,1)*rdtime
            ztop        = z(i,1) - z(i,0)
            rztop       = 1./ztop
            tktop       = 0.5*(tk(i,1) + tk(i,2))
            zmpl        = 0.5*(z(i,2) + z(i,1))
            zm          = 0.5*(z(i,1) + z(i,0))
            delzpl      = zmpl - zm
            fpls        = 1./(ztop*delzpl)
            tmp         = crt - dfntdt(i)*rztop
            diag(i,1)   = tmp + tktop*fpls
            spdiag(i,1) = -tktop*fpls
            rhs(i,1)    = tmp*tin(i,1) + fnt(i)*rztop + htsrc(i,1)
         end do
C
C Intermediate layers
C
         do n=2,plevmx-1
            do ii=1,npts
               i = indx(ii)
               crt         = cmass(i,n)*rho(i,n)*rdtime
               delz        = z(i,n) - z(i,n-1)
               zmpl        = 0.5*(z(i,n+1) + z(i,n))
               zm          = 0.5*(z(i,n)   + z(i,n-1))
               zmmn        = 0.5*(z(i,n-1) + z(i,n-2))
               delzpl      = zmpl - zm
               delzmn      = zm - zmmn
               fpls        = 1./(delz*delzpl)
               fmns        = 1./(delz*delzmn)
               tkpls       = 0.5*(tk(i,n+1)+tk(i,n))
               tkmns       = 0.5*(tk(i,n)+tk(i,n-1))
               sbdiag(i,n) = -tkmns*fmns
               diag(i,n)   = crt + (tkpls*fpls + tkmns*fmns)
               spdiag(i,n) = -tkpls*fpls
               rhs(i,n)    = crt*tin(i,n) + htsrc(i,n)
            end do
         end do
C
C Bottom layer
C
         do ii=1,npts
            i = indx(ii)
            crt        = cmass(i,plevmx)*rho(i,plevmx)*rdtime
            zbot       = z(i,plevmx) - z(i,plevmx-1)
            zm         = 0.5*(z(i,plevmx)   + z(i,plevmx-1))
            zmmn       = 0.5*(z(i,plevmx-1) + z(i,plevmx-2))
            delzmn     = zm - zmmn
            tkbot      = 0.5*(tk(i,plevmx-1) + tk(i,plevmx))
            fmns       = 1./(zbot*delzmn)
            sbdiag(i,plevmx) = -tkbot*fmns
            diag(i,plevmx)   = crt + (tkbot*fmns)
            rhs(i,plevmx)    = crt*tin(i,plevmx) - fbt(i)/zbot
     $                          + htsrc(i,plevmx)
         end do
      end if
C
C For the linear equation ax = b,  a and b arrays are now set;
C solve for temperatures (x):
C
      if (plevmx.eq.1) then
         do ii=1,npts
            i = indx(ii)
            tsbsf(i,1) = rhs(i,1)/diag(i,1)
         end do
      else
         call mtdlss(sbdiag  ,diag   ,spdiag   ,rhs     ,tsbsf ,
     $               plevmx  ,plond  ,npts     ,indx    ,ws      ,
     $               plond*plevmx ,ier     )
         if(ier.ne.0) then
           write(6,*) 'SRFTSB: Error returned from mtdlss'
           call blowup
         end if
      end if
C
      return
      end
 
