#include <params.h>
      subroutine oznint
C-----------------------------------------------------------------------
C
C Time interpolate ozone mixing ratios to current time, reading
C in new monthly data if necessary
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath
C Rewritten for netcdf: J. Rosinski 4/29/97
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$

#if defined ( SPMD )
      use mod_comm, only: mp_bcst_r2d
#endif

C-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
#include <pagrid.h>
C-----------------------------------------------------------------------
#include <comctl.h>
C-----------------------------------------------------------------------
#include <comlun.h>
C-----------------------------------------------------------------------
#include <commap.h>
C-----------------------------------------------------------------------
#include <comozp.h>
C-----------------------------------------------------------------------
#include <comtim.h>
C-----------------------------------------------------------------------
#include <netcdf.inc>
!
! Local workspace
!
      integer cnt4(4)        ! array of counts for each dimension
      integer strt4(4)       ! array of starting indices
      integer yr             ! year number (only relevant if no cycling)
      integer ret            ! return code
      integer j, k           ! indices
      integer ntmp           ! temporary

      real fact1, fact2      ! time interpolation factors
      integer date_ozm 	     ! date of previous data (if ozncyc true, year of data ignored)
      integer date_ozp 	     ! date of future data
      real caldayloc         ! calendar day includes yr 
      real deltat            ! time (days) between interpolating ozone data
      logical found          ! true if ozone data bounds model date
!
      if (masterproc) then
!
!  current date has YYYYYYMMDD.sec
!
        caldayloc = float(mcdate) + float(mcsec)/86400.

        strt4(1) = 1
        strt4(2) = 1
        strt4(3) = 1
        cnt4(1)  = lonsiz
        cnt4(2)  = levsiz
        cnt4(3)  = latsiz
        cnt4(4)  = 1
!
! If model time is past current forward ozone timeslice, read until next
! timeslice is found for interpolation. 
!
        found = .false.
        do while (.not.found)

          if (cdayozm .le. caldayloc .and.
     &                     caldayloc .le. cdayozp) then
            found = .true.
 
          else if (caldayloc .gt. cdayozp) then
            if (ozncyc) then
              np1 = mod(np1,12) + 1
              cdayozm = cdayozp
              yr = mcdate/10000
              if (np1 .gt. 1) then
                cdayozp = float(yr*10000+date_oz(np1)) +          ! year after current year
     &                    float(sec_oz(np1))/86400.
     &                - float(int(date_oz(np1)/10000.))*10000.    ! strip off any year in ozone data
              else
                cdayozp = float((yr+1)*10000+date_oz(np1)) +      ! year after current year
     &                    float(sec_oz(np1))/86400.
     &                - float(int(date_oz(np1)/10000.))*10000.    ! strip off any year in ozone data
              endif

            else
              np1 = np1 + 1
              cdayozm = cdayozp
              cdayozp = float(date_oz(np1)) + float(sec_oz(np1))/86400.
            endif

            if (np1.gt.poztim) then
              if(masterproc) write(6,*)'OZNINT: Attempt to read past end of O3 dataset'
              call endrun
            endif
            if (caldayloc.le.cdayozp) then
              ntmp = nm
              nm = np
              np = ntmp
              strt4(4) = np1
              call NCVGT(ncid_oz,oznid,strt4,cnt4,oznbdy(1,1,1,np),ret)
              found = .true.
              if(masterproc) write(6,*)'OZNINT: Read ozone for date (yyyymmdd) ',
     &                  date_oz(np1),' sec ',sec_oz(np1)
            endif

          else if (caldayloc .lt. cdayozm) then
              if(masterproc) write(6,*)'OZNINT: Model date earlier than Ozone data',
     &                  ' Model =',caldayloc, ' Ozone =',cdayozm
              call endrun

          endif

        enddo
!
!
! Time interpolation.  Get linear interpolation factors.
!
        date_ozm = int(cdayozm)
        date_ozp = int(cdayozp)
        call interp_time(mcdate, mcsec, date_ozm,sec_oz(nm),
     &                   date_ozp,sec_oz(np1), fact1,fact2)
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
        if (abs(fact1+fact2-1.).gt.1.e-6 .or. 
     &      fact1.gt.1.000001 .or. fact1.lt.-1.e-6 .or.
     &      fact2.gt.1.000001 .or. fact2.lt.-1.e-6) then
          if(masterproc) write(6,*)'OZNINT: Bad fact1 and/or fact2=',fact1,fact2
          call endrun
        end if
           
        if (pozlon.eq.1) then
          call lininterp(oznbdy(1,1,1,nm), ozlat, pozlev, pozlat,
     &                   ozmixm(1,1,1,nm), latdeg, plat)
          call lininterp(oznbdy(1,1,1,np), ozlat, pozlev, pozlat,
     &                   ozmixm(1,1,1,np), latdeg, plat)
          do j=1,plat
!!!          do j=begj, endj
            do k=1,pozlev
              ozmix(k,j) = ozmixm(1,k,j,nm)*fact1 +
     &                     ozmixm(1,k,j,np)*fact2
            end do
          end do
        else
          if(masterproc) write(6,*)'OZNINT: Code to interp 3-d ozone dataset must ',
     &              'be tested!'
          call endrun
        end if

      end if

#if defined ( SPMD )
      if(masterproc) then
         if(pozlev > plon ) then
            print *, 'pozlev > plon: size of g_2d for MLP maybe too small'
            stop
         endif 
      endif

! Be careful: pozlev must be smaller than plon (see above remark)
      call mp_bcst_r2d(pozlev, plat, begj, endj, ozmix, 0)
#endif

      return
      end
