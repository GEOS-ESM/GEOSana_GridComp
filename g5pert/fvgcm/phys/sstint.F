#include <params.h>
      subroutine sstint
C-----------------------------------------------------------------------
C
C Time interpolate sea surface temps. to current time, reading in new 
C monthly data if necessary.
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
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer iunit              ! Logical unit number
C
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C-----------------------------------------------------------------------
#include <parsst.h>
C------------------------------Commons----------------------------------
#include <comctl.h>
C-----------------------------------------------------------------------
#include <comlun.h>
C-----------------------------------------------------------------------
#include <comsst.h>
C-----------------------------------------------------------------------
#include <comtim.h>
C-----------------------------------------------------------------------
#include <netcdf.inc>
C---------------------------Local variables-----------------------------
C
!
! Local workspace
!
      integer cnt3(3)        ! array of counts for each dimension
      integer strt3(3)       ! array of starting indices
      integer yr             ! year number 
      integer ret            ! return code
      integer i,j            ! indices
      integer ntmp           ! temporary
      integer date_sstm       ! date of previous data (if sstcyc true, year of data ignored)
      integer date_sstp       ! date of future data


      real fact1, fact2      ! time interpolation factors
      real caldayloc         ! calendar day YYYYMMDD.sec
      real deltat            ! time (days) between interpolating sst data
      logical found          ! true if model date has bounding sst data
C
C-----------------------------------------------------------------------
C
!
! Master does all the work.  Sends needed info to slaves
      if (masterproc) then
!
!  current date has YYYYYYMMDD.sec
!
        caldayloc = float(mcdate) + float(mcsec)/86400.

        strt3(1) = 1
        strt3(2) = 1
        strt3(3) = 1
        cnt3(1)  = lonsiz
        cnt3(2)  = latsiz
        cnt3(3)  = 1
!
! If model time is past current forward sst timeslice, read until next
! timeslice is found for interpolation. 
!
        found = .false.
        do while (.not.found)
          if (cdaysstm .le. caldayloc .and. 
     &                      caldayloc .le. cdaysstp) then
            found = .true.

          else if (caldayloc .gt. cdaysstp) then

            if (sstcyc) then
              np1 = mod(np1,timesiz) + 1			! if at end reset to 1
              cdaysstm = cdaysstp
              yr = mcdate/10000
              if (np1 .gt. 1) then
                cdaysstp = float(yr*10000+date_sst(np1)) +	! year after current year
     &                     float(sec_sst(np1))/86400.
              else
                cdaysstp = float((yr+1)*10000+date_sst(np1)) +	! year after current year
     &                     float(sec_sst(np1))/86400.
              endif

            else
              np1 = np1 + 1
              cdaysstm = cdaysstp
              cdaysstp = float(date_sst(np1)) + float(sec_sst(np1))/86400.
            endif

            if (np1.gt.timesiz) then
              if(masterproc) write(6,*)'SSTINT: Attempt to read past end of SST dataset'
              call endrun
            endif

            if (caldayloc.le.cdaysstp) then
              ntmp = nm
              nm = np
              np = ntmp
              strt3(3) = np1
              call NCVGT(ncid_sst,sstid,strt3,cnt3,sstbdy(1,1,np),ret)
              if(masterproc) write(6,*)'SSTINT: Read sst for date (yyyymmdd) ',
     &                  date_sst(np1),' sec ',sec_sst(np1)
              found = .true.
            endif

          else if (caldayloc .lt. cdaysstm) then
              if(masterproc) write(6,*)'SSTINT: Model date earlier than SST data',
     &                  ' Model =',caldayloc, ' SST =',cdaysstm
              call endrun
          endif

        enddo
!
! Time interpolation.  Get linear interpolation factors.
!
        date_sstm = int(cdaysstm)
        date_sstp = int(cdaysstp)
        call interp_time(mcdate, mcsec, date_sstm,sec_sst(nm), 
     &                   date_sstp,sec_sst(np1), fact1,fact2)
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
        if (abs(fact1+fact2-1.).gt.1.e-6 .or. 
     &      fact1.gt.1.000001 .or. fact1.lt.-1.e-6 .or.
     &      fact2.gt.1.000001 .or. fact2.lt.-1.e-6) then
          if(masterproc) write(6,*)'SSTINT: Bad fact1 and/or fact2=',fact1,fact2
          call endrun
        end if

!$omp parallel do private(i,j)
        do j=1,plat
!!!        do j=begj,endj
           do i=1,plon
              if (sstbdy(i,j,nm).le.tsice .and. fact2.lt.0.5) then
                 sst(i,j) = sstbdy(i,j,nm)
              else if (sstbdy(i,j,np).le.tsice .and. fact1.le.0.5) then
                 sst(i,j) = sstbdy(i,j,np)
              else
                 sst(i,j) = sstbdy(i,j,nm)*fact1 +
     &                sstbdy(i,j,np)*fact2
              end if
           end do
        end do
      end if

#if defined ( SPMD )
      call mp_bcst_r2d(plon, plat, begj, endj, sst, 0)      
#endif

      return
      end
