#include <preproc.h>
      subroutine lsmmain(atmlon  ,atmlat  ,beglatspmd,endlatspmd,pcnst,
     $                   iradsw  ,nstep   ,mbdate  ,mdbase  ,msbase  ,
     $                   mbsec   ,dtime   ,pbot    ,thbot   ,ubot    ,
     $                   vbot    ,tbot    ,qbot    ,zbot    ,flwds   ,
     $                   precc   ,precl   ,sols    ,soll    ,solsd   ,
     $                   solld   ,oro     ,shf     ,lhf     ,cflx    ,
     $                   wsx     ,wsy     ,ts      ,tref    ,asdir   ,
     $                   asdif   ,aldir   ,aldif   ,snowh   ,lwup    ,
     $                   tssub   ,plevmx  ,doalb   ,mcdate  ,mcsec   ,
     $                   z0m     ,z0h     ,zpd     ,wet1    ,wet     ,
     $                   trad    ,emsfc   )

! doalb, mcdate, and mcsec added by SJL

!-----------------------------------------------------------------------

! Call land surface model. This is interface code. It is
! the time-stepping interface between the land surface model 
! and CCM. It changes units as needed. A similar interface is 
! needed regardless of what land surface model is used.
 
!-----------------------------------------------------------------------
! $Id$
! $Author$
* ------------------------ code history ---------------------------
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------
#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>
#include <basinpar.h>
!------------------------- Input variables -----------------------------
      integer atmlon               !atm number of longitudes
      integer atmlat               !atm number of latitudes
      integer pcnst                !atm number of constituents (water is first)
      integer iradsw               !Shortwave frequency
      integer nstep                !Current timestep
      integer mdbase               !Base day of run (e.g., 0)
      integer msbase               !Base seconds of base day (e.g., 0)
      integer mbdate               !Base date of run (yymmdd format) (e.g., 000901)
      integer mbsec                !Base seconds of base date (e.g., 0)
! SJL
      integer mcdate
      integer mcsec

      real dtime                   !length of timestep
      integer beglatspmd    !beg lat on proc, degenerates to 1 for nonspmd
      integer endlatspmd    !end lat on proc, degenerates to atmlat for nonspmd
      integer plevmx                  ! number of atm subsurface levels

      real  pbot(atmlon,beglatspmd:endlatspmd)    !bottom level pressure
      real thbot(atmlon,beglatspmd:endlatspmd)    !bottom lev potential temp 
      real  ubot(atmlon,beglatspmd:endlatspmd)    !bottom level u wind
      real  vbot(atmlon,beglatspmd:endlatspmd)    !bottom level v wind
      real  tbot(atmlon,beglatspmd:endlatspmd)    !bottom level temperature 
      real  qbot(atmlon,beglatspmd:endlatspmd)    !bottom lev specific humidity
      real  zbot(atmlon,beglatspmd:endlatspmd)    !btm lev hght above surface
      real precc(atmlon,beglatspmd:endlatspmd)    !convective precip rate (m/s)
      real precl(atmlon,beglatspmd:endlatspmd)    !large-scale precip rate(m/s)
      real flwds(atmlon,beglatspmd:endlatspmd)    !downward longwave radiation at surface
      real  sols(atmlon,beglatspmd:endlatspmd)    !direct beam solar radiation onto srf
      real  soll(atmlon,beglatspmd:endlatspmd)    !direct beam solar radiation onto srf
      real solsd(atmlon,beglatspmd:endlatspmd)    !diffuse solar radiation onto srf (sw)
      real solld(atmlon,beglatspmd:endlatspmd)    !diffuse solar radiation onto srf (lw)
      real   oro(atmlon,beglatspmd:endlatspmd)    !land/ocean/sea ice flag

!------------------------- Output variables -----------------------------

      real   lhf(atmlon,beglatspmd:endlatspmd)      ! latent heat flux
      real   shf(atmlon,beglatspmd:endlatspmd)      ! sensible heat flux
      real   wsx(atmlon,beglatspmd:endlatspmd)      ! surface u-stress (N)
      real   wsy(atmlon,beglatspmd:endlatspmd)      ! surface v-stress (N)
      real  tref(atmlon,beglatspmd:endlatspmd)     ! ref height surface air temp
      real asdir(atmlon,beglatspmd:endlatspmd)    ! albedo: shortwave, direct
      real asdif(atmlon,beglatspmd:endlatspmd)    ! albedo: shortwave, diffuse
      real aldir(atmlon,beglatspmd:endlatspmd)    ! albedo: longwave, direct
      real aldif(atmlon,beglatspmd:endlatspmd)    ! albedo: longwave, diffuse
      real  lwup(atmlon,beglatspmd:endlatspmd)    ! longwave up radiative flux
      real    ts(atmlon,beglatspmd:endlatspmd)       ! surface temperature
      real snowh(atmlon,beglatspmd:endlatspmd)    ! ccm snow depth (liquid water, m )

      real   cflx(atmlon,pcnst,beglatspmd:endlatspmd) ! constituent flux (evap)
      real tssub(atmlon,plevmx,beglatspmd:endlatspmd) ! atm surface/subsurface temperatures 
! JDC modify begin
      real   z0m(atmlon,beglatspmd:endlatspmd)     !roughness length, momentum (m)
      real   z0h(atmlon,beglatspmd:endlatspmd)     !roughness length, sensible heat (m)
      real   zpd(atmlon,beglatspmd:endlatspmd)     !displacement height (m)
      real   wet1(atmlon,beglatspmd:endlatspmd)    !top soil layer wetness (fraction)
      real   wet(atmlon,beglatspmd:endlatspmd)     !root zone soil wetness (fraction)
      real   trad(atmlon,beglatspmd:endlatspmd)    !surface brightness temperature (K)
      real   emsfc(atmlon,beglatspmd:endlatspmd)   !bulk surface emissivity (fraction)
! JDC modify end

!----------------------------Local variables----------------------------
      integer i,j,m                ! Loop index
      logical doalb                ! Is next timestep a radiation time?

      real prcxy (atmlon,beglatspmd:endlatspmd)   ! Conv precip (mm h2o/s)
      real prlxy (atmlon,beglatspmd:endlatspmd)   ! Lrge-scale precip(mm h2o/s)
      real snowxy(atmlon,beglatspmd:endlatspmd)   ! ccm snow depth (liquid water, mm)
!------------------------------Externals--------------------------------
 
      external lsmdrv 
 
! Determine if albedo calculation is done on this timestep
 
! SJL doalb is determined from outside!
!     doalb=iradsw.eq.1 .or. (mod(nstep,iradsw).eq.0 .and. nstep+1.ne.1)
 
! Change units of precipitation rate from m/s to mm/s.

!$omp parallel do default(shared) private(i,j)

      do j=beglatspmd,endlatspmd
         do i=1,atmlon
            prcxy(i,j) = precc(i,j)*1000.           
            prlxy(i,j) = precl(i,j)*1000.           
         end do
      end do

! Calculate lsm variables
 
      call lsmdrv(atmlon ,atmlat ,beglatspmd ,endlatspmd ,pcnst ,
     $            pbot   ,thbot  ,ubot       ,vbot       ,tbot  ,
     $            qbot   ,flwds  ,zbot       ,prcxy      ,prlxy ,
     $            sols   ,soll   ,solsd      ,solld      ,doalb ,
     $            nstep  ,ts     ,shf        ,lhf        ,cflx  ,
     $            wsx    ,wsy    ,tref       ,asdir      ,asdif ,
     $            aldir  ,aldif  ,lwup       ,snowxy     ,   
     $            mcdate ,mcsec  ,z0m        ,z0h        ,zpd   ,
     $            wet1   ,wet    ,trad       ,emsfc      )
 
! Set CCM3 tssub to ts for land points
      
!$omp parallel do default(shared) private(i,j,m)

      do j=beglatspmd,endlatspmd
         do m=1,plevmx
            do i=1,atmlon
               if (oro(i,j) == 1.) tssub(i,m,j) = ts(i,j)
            end do
         end do
! Set CCM3 snow depth for land points
         do i=1,atmlon
            if (oro(i,j) == 1.) snowh(i,j) = snowxy(i,j)*0.001
         end do
      end do

      return
      end
