#include <params.h>
      subroutine cond(tdt     ,pmid    ,pdel    ,t       ,q       ,
     $                qc      ,precl   ,evap)
C-----------------------------------------------------------------------
C
C Calculate large scale condensation
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Buja, Jun 1992, Feb 1996
C Reviewed:          J. Hack, G. Taylor, Aug 1992
C                    J. Hack, Feb 1996 
C
C Ice-phase latent heat associated with frezzing and melting
C added by S.-J. Lin
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
#include <comadj.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real tdt                    ! Physics time step (2 delta t)
      real pmid(plond,plev)       ! Pressure at layer midpoints
      real pdel(plond,plev)       ! Delta p at each model level
C
C Output arguments
C
      real t(plond,plev)          ! Temperature
      real q(plond,plev)          ! Specific humidity
      real qc(plond,plev)         ! Moisture tendency due to rainout
      real precl(plond)           ! Large-scale precipitation rate
      real evap(plond,plev)       ! Water evaporation rate
C
C---------------------------Local variables-----------------------------
C
      real absqs                  ! Intermediate quantity
      real denom                  ! Intermediate quantity
      real dqsdt                  ! Change of qsat with respect to temp.
      real est(plond,plev)        ! Saturation vapor pressure
      real omeps                  ! 1 - 0.622
      real qsat(plond,plev)       ! Saturation specific humidity
      real rain(plond)            ! Rain (units of kg/m^2 s)
      real rga                    ! Reciprocal gravitatnl acceleration
      real rhm1                   ! RH - saturation RH
      real zqcd(plond)            ! Intermed quantity (several actually)
      real zqdt                   ! Reciprocal of tdt
      real cndwtr(plond,plev)     ! Water condensation rate (kg/m**2/s)
      real ke                     ! `disposable parameter' in evaporation
      real relhum                 ! Relative humidity
      real dpovrg                 ! deltap/grav
      integer i                   ! Longitude index
      integer jiter               ! Iteration counter
      integer k                   ! Vertical index

C Added by SJL 
      integer lcond               ! starting layer for cond. computation
      real ice(plond)             ! Solid precipation (ice , snow, etc)
      real pice(plond)            ! Temp storage for solid precipatation
      real slcp                   ! Ice-phase latent heat / CP
      real t_ice                  ! frezzing point in Kelvin
C Temporatory 
      real dice
      real h0
      real h1

      t_ice = 273.16
      slcp = 3.34E5/cpair
C
C-----------------------------------------------------------------------
C
      rga   = 1./gravit
      zqdt  = 1./tdt
      omeps = 1. - epsilo

C SJL: Compute starting L-S condensation level
!      lcond = 1
!      do k=2, plev
!         if(pmid(1,k) .gt. 2000.) then
!            lcond = k
!            go to 111
!         endif
!      enddo
!111   continue

      lcond = 2

C
C First diagnose condensation rate due to stable processes
C Update column T and Q (evaporation process is `time-split')
C Condensation calculation is hard-wired for two iterations
C
      do k=lcond,plev
        do i=1,plon
          cndwtr(i,k) = 0.0
        end do
      end do
 
      do i=1,plon
         ice(i) = 0.
      enddo

      do jiter=1,2
 
        do i=1,plon
           pice(i) = 0.
        enddo

        call aqsat(t       ,pmid    ,est     ,qsat    ,plond   ,
     $             plon    ,plev    ,lcond   ,plev    )
        do k=lcond,plev
C
C Calculate condensation-rate and new t- and q-values
C
          do i=1,plon
C Use of critical saturation vapor pressure requires coefficient on the
C term omeps*est(i,k) in the next statement (e.g. omeps*est(i,k)*escrit)
C Corresponding changes must also be incorporated into estabv.for (e.g.,
C terms est(i,k) in qsat evaluation become escrit*est(i,k))

            denom   = (pmid(i,k) - omeps*est(i,k))*t(i,k)**2
            dqsdt   = clrh2o*qsat(i,k)*pmid(i,k)/denom
            absqs   = abs(qsat(i,k))
            rhm1    = q(i,k)/qsat(i,k) - 1.
            zqcd(i) = max(absqs*rhm1/(1. + cldcp*dqsdt),0.)
            if (q(i,k) .lt. 0.0) zqcd(i) = 0.
            q(i,k)  = q(i,k) - zqcd(i)
            t(i,k)  = t(i,k) + zqcd(i)*cldcp
            h1  = pdel(i,k)*rga*zqdt
            cndwtr(i,k) = cndwtr(i,k) + zqcd(i)*h1
            qc    (i,k) = qc(i,k)     + zqcd(i)*zqdt

          if(t(i,k) .lt. t_ice .and. cndwtr(i,k) .gt. 0.) then
C ********************************
C ***   Freezing *****
C ********************************
             h0   = slcp/h1
C The release of latent heat due to freezing should not warm
C the atmosphere above frezzing point.
             dice = min(cndwtr(i,k), (t_ice-t(i,k))/h0)
             pice(i) = pice(i) + dice
             t(i,k) = t(i,k) + dice*h0
             cndwtr(i,k) = cndwtr(i,k) - dice
        elseif(t(i,k) .gt. t_ice .and. pice(i) .gt. 0.) then
C ********************************
C ***   Melting *****
C ********************************
             h0   = slcp/h1
C The heat needed for melting should not cool the atmosphere
C below frezzing
             dice = min(pice(i), (t(i,k)-t_ice)/h0)
             cndwtr(i,k) = cndwtr(i,k) + dice
             t(i,k) = t(i,k) - dice*h0
             pice(i) = pice(i) - dice
          endif
          end do
        end do
 
      do i=1,plon
         ice(i) = ice(i) + pice(i)
      enddo

      end do
 
C Initialize rain vector (will be updated as rain falls through column)

      do i=1,plon
        rain(i) = max(cndwtr(i,lcond),0.0)
      end do

      call aqsat(t       ,pmid    ,est     ,qsat    ,plond   ,
     $           plon    ,plev    ,lcond   ,plev    )
C
C Evaporate condensate on the way down (see Sundqvist, 1988: Physically
C Based Modelling ..., pp 433-461, Schlesinger, Ed., Kluwer Academic)
C variable evap has units of 1/s; variable rain has units of kg/m**2/s
C rain is used to accumuluate unevaporated rain water on the way down
 
#if defined (ALT_CPS)
      ke = 1.2e-5
#else
      ke = 1.0e-5                     ! set in common block in final code
#endif

      do k=lcond+1,plev
        do i=1,plon
          dpovrg  = pdel(i,k)*rga
          relhum  = q(i,k)/qsat(i,k)
          evap(i,k)    = max(ke*(1.0 - relhum)*sqrt(rain(i)), 0.0)
          evap(i,k)    = min(evap(i,k), (qsat(i,k)-q(i,k))/tdt)
          evap(i,k)    = min(rain(i)/dpovrg,evap(i,k))
          qc(i,k) = qc(i,k) - evap(i,k)
          q(i,k)  = q(i,k) + evap(i,k)*tdt
          t(i,k)  = t(i,k) - evap(i,k)*tdt*cldcp
          rain(i) = max(rain(i) - evap(i,k)*dpovrg + cndwtr(i,k),0.0)
        end do
      end do
      do i=1,plon
        precl(i) = (rain(i)+ice(i)) / rhoh2o
      end do
C
      return
      end
