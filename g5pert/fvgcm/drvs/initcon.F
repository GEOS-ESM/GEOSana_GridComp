#include <params.h>
      subroutine initcon(ccm3)

! This routine was taken from CCM3's initcom and inti

#if defined( SPMD )
        use mod_comm, only : gid
#endif
      implicit none


c Parameters

#include <pmgrid.h>
#include <pagrid.h>

c Commons

#include <comcon.h>
#include <comctl.h>
#include <comqmin.h>
C SJL 11/1/98
#include <commap.h>
#include <comtim.h>
#include <comtrcnm.h>

c---------------------------Local variables-----------------------
      real ak(plev+1)
      real bk(plev+1)
      real ptop, pint
      real prsl
      integer ks
      integer k

      integer j           ! Latitude index
      integer m           ! Index for tracer
      logical  ccm3
      real     clrh2o     ! Ratio of latent heat to gas constant (h2o)
      real pi
      real dp
      real edge
      character*3 trnum   ! Advected species number (Character)

#if ( !defined SPMD )
      integer gid
      gid = 0
#endif

      anncyc = .true.
      interp = .false.
      itsst = 1              ! update frequency for SST

c Climate boundary data set
c GSL 10/8/98 - move ozncyc, sstcyc to top-level namelist
c     ozncyc  = .true.       ! If true, cycle ozone dataset
c     sstcyc  = .true.       ! If true, cycle sst dataset
      incorrad  = .true.     ! keep abs/ems buffer in-core

c Set up the reference vertical coordinate system
      call set_eta(plev, ks, ptop, pint, ak, bk)
      do k=plev,1,-1
         prsl = 0.5*(ak(k) + ak(k+1) + (bk(k)+bk(k+1)) * 100000.)
     &          / 100.
c         if(gid==0) write(6,*) k, prsl
      enddo

      if(ccm3) then

         adiabatic = .false.

      do m=1,pcnst
c SJL 4/2/98
         if(gid==0) write(unit=trnum,fmt='(i3)') m+100
c        hadvnam(m)  = 'HA'//trnum(2:3)
c        vadvnam(m)  = 'VA'//trnum(2:3)
         vdiffnam(m) = 'VD'//trnum(2:3)
         if(gid==0) write(6,*) 'vdiffnam=', vdiffnam(m)
c        dcconnam(m) = 'DC'//trnum(2:3)
c        fixcnam(m)  = 'DF'//trnum(2:3)
c        tendnam(m)  = 'TE'//trnum(2:3)
c        tottnam(m)  = 'TA'//trnum(2:3)
         qmin(m)     =  0.0
      end do
         qmin(1)     =  1.e-12     ! Minimum mixing ratio for moisture

c C Set the minimum mixing ratio for the counter-gradient term.
c Normally this should be the same as qmin, but in order to match
c control case 414 use zero for water vapor.

      qmincg(1) = 0.
      if(pcnst .ge. 2) then
      do m=2,pcnst
         qmincg(m) = qmin(m)
      end do
      endif

      else
         adiabatic = .true.
      endif

c Initialize COMCON

      cpair  = 1004.64
      epsilo = 0.622
      gravit = 9.80616
      latvap = 2.5104e06
      latice = 3.336e5

      omega  = 7.292e-5
      rair   = 287.04
      rearth = 6.37122e+6
      rhoh2o = 1.e3

      stebol = 5.67e-8
      tmelt  = 273.16

      cappa  = rair/cpair
      ez     = omega/sqrt(0.375)
      ra     = 1./rearth
      rga    = 1./gravit

c In adiabatic case, set zvir and cpvir explicitly to zero instead of
c computing as (rh2o/rair - 1.) and (cpwv/cpair - 1.) respectively,
c in order to guarantee an identical zero.

      if (adiabatic) then
        rh2o   = rair
        zvir   = 0.
        cpwv   = cpair
        cpvir  = 0.
      else
        rh2o   = 4.61e2
        zvir   = rh2o/rair - 1.
        cpwv   = 1.81e3
        cpvir  = cpwv/cpair - 1.

        clrh2o = latvap/rh2o

c Latitude array (list of latitudes in radians and degrees)

      pi = 4. * atan(1.)
      dp = pi / float(plat-1)

      do j=1,plat
          clat(j) = -pi*0.5 + (j-1)*dp
        latdeg(j) = -90. + float(j-1)*180./float(plat-1)
      end do

c Poles: for physics
          edge    = 0.5 * (clat(1)   + clat(2))
          clat(1) = 0.5 * (clat(1)   + edge   )
        edge      = 0.5 * (latdeg(1) + latdeg(2))
        latdeg(1) = 0.5 * (latdeg(1) + edge     )

          edge       = 0.5 * (clat(plat)  + clat(plat-1))
          clat(plat) = 0.5 * (clat(plat)  + edge        )
        edge         = 0.5 * (latdeg(plat) + latdeg(plat-1))
        latdeg(plat) = 0.5 * (latdeg(plat) + edge          )

c Time independent initialization for multiple gravity wave
c parameterization.

      call mgwinti(cpair   ,cpwv    ,gravit  ,rair    )

c Initialization of time independent fields for vertical diffusion.
c Call initialization routine for boundary layer scheme.

      call vdinti (cpwv    ,cpair   ,gravit  ,rair , ak, bk )

c Initialize surface temperature calculation constants

      call tsinti (tmelt   ,latvap  ,rair    ,stebol  ,dtime   )

c Initialize various constants for radiation scheme; note that
c the radiation scheme uses cgs units.

      call radini (gravit  ,cpair   ,epsilo  ,stebol  )

c Initialize convective adjustment common block /comadj/
c Note that the variable nlvdry in comadj is set to its default
c value in preset.for, and updated from namelist input (if necessary)
c in subroutine DATA.
c SJL: nlvdry set in setup.

      call cainti (cappa   ,cpair   ,epsilo  ,gravit  ,latvap  ,
     $             rhoh2o  ,clrh2o  )

c Initialize es lookup tables

      call esinti (epsilo  ,latvap  ,latice  ,rh2o    ,cpair   )

c Initialize moist convective mass flux procedure common block, cmfmca
c the following made use of the hybrid coord. info.; need to do something
c about that.  !!!  SJL

      call mfinti (rair    ,cpair   ,gravit  ,latvap  ,rhoh2o, ak)

      end if

      return
      end
