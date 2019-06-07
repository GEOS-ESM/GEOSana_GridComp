#include <params.h>
      subroutine oznini
C-----------------------------------------------------------------------
C
C Do initial read of time-variant ozone boundary dataset, containing
C ozone mixing ratios as a function of latitude and pressure.  Read two
C consecutive months between which the current date lies.  Routine
C RADOZ2 then evaluates the two path length integrals (with and without
C pressure weighting) from zero to the interfaces between the input
C levels.  It also stores the contribution to the integral from each
C layer.
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

c
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
#include <comozp.h>
C-----------------------------------------------------------------------
#include <comtim.h>
C-----------------------------------------------------------------------
#include <netcdf.inc>
!
! Local workspace
!
      character*(MAXNCNAM) lonnam   ! Name of longitude dimension 
      character*(MAXNCNAM) levnam   ! Name of level dimension 
      character*(MAXNCNAM) latnam   ! Name of latitude dimension 
      character*(MAXNCNAM) timenam  ! Name of time dimension 

      integer dateid                ! netcdf id for date variable
      integer secid                 ! netcdf id for seconds variable
      integer lonid                 ! netcdf id for longitude variable
      integer latid                 ! netcdf id for latitude variable
      integer levid                 ! netcdf id for level variable
      integer timeid                ! netcdf id for time variable
      integer cnt4(4)               ! array of counts for each dimension
      integer strt4(4)              ! array of starting indices
      integer yr                    ! year number (only relevant if no cycling)
      integer ret                   ! return code
      integer n, k                  ! indices

      real caldayloc                ! calendar day includes yr 

      nm = 1
      np = 2
!
! SPMD: Master does all the work.  Sends data to slaves
!
      if (masterproc) then
!
!  current date has YYYYYYMMDD.sec
!
        caldayloc = float(mcdate) + float(mcsec)/86400.

!
! Get and check dimension info
!
        lonid = NCDID(ncid_oz, 'lon', ret)
        levid = NCDID(ncid_oz, 'lev', ret)
        latid = NCDID(ncid_oz, 'lat', ret)
        timeid = NCDID(ncid_oz, 'time', ret)

        call NCDINQ(ncid_oz, lonid, lonnam, lonsiz, ret)
        if (lonsiz .ne. pozlon) then
          if(masterproc) write(6,*)'OZNINI: lonsiz=',lonsiz,' must = pozlon=',pozlon
          call endrun
        end if

        call NCDINQ(ncid_oz, levid, levnam, levsiz, ret)
        if (levsiz .ne. pozlev) then
          if(masterproc) write(6,*)'OZNINI: levsiz=',levsiz,' must = pozlev=',pozlev
          call endrun
        end if

        call NCDINQ(ncid_oz, latid, latnam, latsiz, ret)
        if (latsiz .ne. pozlat) then
          if(masterproc) write(6,*)'OZNINI: latsiz=',latsiz,' must = pozlat=',pozlat
          call endrun
        end if

        call NCDINQ(ncid_oz, timeid, timenam, timesiz, ret)
        if (timesiz .ne. poztim) then
          if(masterproc) write(6,*)'OZNINI: timesiz=',timesiz,' must = poztim=',poztim
          call endrun
        end if

        dateid = NCVID(ncid_oz, 'date', ret)
        secid = NCVID(ncid_oz, 'datesec', ret)
c       secid = NCVID(ncid_oz, 'time', ret)
        oznid = NCVID(ncid_oz, 'OZONE', ret)
        latid = NCVID(ncid_oz, 'lat', ret)
        levid = NCVID(ncid_oz, 'lev', ret)
!
! Retrieve latitude and level arrays for interpolation.
!
        call NCVGT(ncid_oz, latid, 1, latsiz, ozlat, ret)
        call NCVGT(ncid_oz, levid, 1, levsiz, pin, ret)

!
! Convert from millibars to pascals
!
        do k=1,pozlev
          pin(k) = pin(k)*100.
        end do
!
! Retrieve entire date and sec variables. Only 12 month clims allowed
!
        call NCVGT(ncid_oz, dateid, 1, timesiz, date_oz, ret)
        call NCVGT(ncid_oz, secid,  1, timesiz, sec_oz,  ret)
        if (ozncyc) then
          if (mod(date_oz(1),10000)/100.ne.1) then
            if(masterproc) write(6,*)'When cycling ozone, 1st month must be 1'
            call endrun
          end if
          if (mod(date_oz(poztim),10000)/100.ne.12) then
            if(masterproc) write(6,*)'When cycling ozone, last month must be 12'
            call endrun
          end if
        end if

        strt4(1) = 1
        strt4(2) = 1
        strt4(3) = 1
        cnt4(1)  = lonsiz
        cnt4(2)  = levsiz
        cnt4(3)  = latsiz
        cnt4(4)  = 1
!
! Special code for interpolation between December and January
!
        if (ozncyc) then
          n = 12
          np1 = 1
          yr = mcdate/10000
          cdayozm = float((yr-1)*10000+date_oz(n)) +
     &                float(sec_oz(n))/86400.
     &              - float(int(date_oz(n)/10000.))*10000.	   ! strip off any year in ozone data
          cdayozp = float(yr*10000+date_oz(np1)) +
     &               float(sec_oz(np1))/86400.
     &              - float(int(date_oz(np1)/10000.))*10000.	   ! strip off any year in ozone data


          if (caldayloc.le.cdayozp .and. caldayloc.gt.cdayozm) then
            strt4(4) = n
            call NCVGT(ncid_oz,oznid,strt4,cnt4,oznbdy(1,1,1,nm),ret)
            strt4(4) = np1
            call NCVGT(ncid_oz,oznid,strt4,cnt4,oznbdy(1,1,1,np),ret)
            goto 10
          end if
        end if
!
! Normal interpolation between consecutive time slices.
!
        do n=1,timesiz-1
          np1 = n + 1
          cdayozm = float(date_oz(n)) + float(sec_oz(n))/86400.
          cdayozp = float(date_oz(np1)) + float(sec_oz(np1))/86400.
          if (ozncyc) then
            yr = mcdate/10000
            cdayozm = cdayozm + float(yr)*10000. - 
     .                float(int(cdayozm/10000.))*10000.	   ! strip off any year in ozone data
            cdayozp = cdayozp + float(yr)*10000. - 
     .                float(int(cdayozp/10000.))*10000.	   ! strip off any year in ozone data
          endif
          if (caldayloc.gt.cdayozm .and. caldayloc.le.cdayozp) then
            strt4(4) = n
            call NCVGT(ncid_oz,oznid,strt4,cnt4,oznbdy(1,1,1,nm),ret)
            strt4(4) = np1
            call NCVGT(ncid_oz,oznid,strt4,cnt4,oznbdy(1,1,1,np),ret)
            goto 10
          end if
        end do

c  Last check for date after last record in climate file
        if (ozncyc) then
          n = 12
          np1 = 1
          yr = mcdate/10000
          cdayozm = float(yr*10000+ date_oz(n)) + 
     $               float(sec_oz(n))/86400.
     &              - float(int(date_oz(n)/10000.))*10000.	   ! strip off any year in ozone data
          cdayozp = float((yr+1)*10000 + date_oz(np1)) + 
     $               float(sec_oz(np1))/86400.
     &              - float(int(date_oz(np1)/10000.))*10000.	   ! strip off any year in ozone data
          if (caldayloc.gt.cdayozm .and. caldayloc.le.cdayozp) then
            strt4(4) = n
            call NCVGT(ncid_oz,oznid,strt4,cnt4,oznbdy(1,1,1,nm),ret)
            strt4(4) = np1
            call NCVGT(ncid_oz,oznid,strt4,cnt4,oznbdy(1,1,1,np),ret)
            goto 10
          end if  
        endif


        if(masterproc) write(6,*)'OZNINI: Failed to find dates bracketing ',
     $            'mcdate, mcsec=', mcdate, mcsec
        call endrun
   10   continue
        if(masterproc) write(6,*)'OZNINI: Read ozone data for dates ',date_oz(n), 
     &            sec_oz(n),' and ',date_oz(np1),sec_oz(np1)

      end if

#if defined ( SPMD )
      call mp_bcst_r2d(pozlev, 1, 1, 1, pin, 0)
#endif

      return
      end

 
