#include <params.h>
      subroutine ccmoce (nymd, nhms, ndt)

C-----------------------------------------------------------------------
C
C CCM3 ocean and sea ice surface fluxes.
C
C----------------------------Code History-------------------------------
C
C Original version: 
C Standardized:      L. Buja, Feb 1996
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
#include <comsrf.h>
C-----------------------------------------------------------------------
#include <commap.h>
C-----------------------------------------------------------------------
#include <comctl.h>
C-----------------------------------------------------------------------
#include <comtim.h>
C--------------------------Local Variables------------------------------
C
      real cosznext(plond)      ! cosine solar zenith angle next timestep
      real cdaynext             ! calendar day for next timestep
      integer j                 ! longitude index
      integer idum1,idum2
C
C-----------------------------------------------------------------------
C
! Added input by SJL
      integer nymd
      integer nhms
      integer ndt

! SJL
! Calendar day for next time step
      idum1 = nymd  
      idum2 = nhms  
      call tick( idum1, idum2, ndt )
      call mcalday(idum1, idum2, cdaynext)
! SJL

!$omp parallel do default(shared) private(j,cosznext)

      do j=begj,endj
 
C Update prescribed sea surface temperatures and sea ice distribution as well
C as snow cover over sea ice.
 
        if (anncyc .and. mod(nstep,itsst).eq.0) then
          call sstan(j, oro(1,j), ts(1,j), tssub(1,1,j), snowh(1,j),
     $               sicthk(1,j)  )
        end if
 
C Ocean surface fluxes and temperatures
 
         call srfoce(oro(1,j)     ,ubot(1,j)   ,vbot(1,j) ,tbot(1,j) ,
     $               qbot(1,j)    ,thbot(1,j)  ,zbot(1,j) ,pbot(1,j) ,
     $               tssub(1,1,j) ,cflx(1,1,j) ,wsx(1,j)  ,wsy(1,j)  ,
     $               ts(1,j)      ,shf(1,j)    ,lhf(1,j)  ,lwup(1,j) )
 
C Sea ice surface fluxes and temperatures
 
         call srfsice(oro(1,j) , snowh(1,j), ubot(1,j), vbot(1,j),
     $                tbot(1,j),
     $                qbot(1,j), thbot(1,j), zbot(1,j), pbot(1,j),
     $                srfrad(1,j),
     $                tssub(1,1,j), cflx(1,1,j), wsx(1,j), wsy(1,j),
     $                ts(1,j),
     $                shf(1,j), lhf(1,j), lwup(1,j) ) 
 
C JDC ADDED
C z0m,z0h and zpd diagnostics 
         call ocez0(
     &        j        ,plon      ,
     &        oro(1,j) ,tbot(1,j) ,pbot(1,j) ,wsx(1,j) ,wsy(1,j) ,
     &        z0m(1,j) ,z0h(1,j)  ,zpd(1,j)  ) 
         
C Albedos for next time step 
 
         call zenith(cdaynext,dodiavg ,clat(j) ,cosznext)
         call albocean(oro(1,j), snowh(1,j), cosznext, asdir(1,j),
     $                 aldir(1,j),
     $                 asdif(1,j), aldif(1,j) )
 
      end do

      return
      end
 
