#include <params.h>
      subroutine sstini
C-----------------------------------------------------------------------
C
C Initialize the procedure for specifying sea surface temperatures
C Do initial read of time-varying sst boundary dataset, reading two
C consecutive months on either side of the current model date.
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath
C Standardized:      L. Buja, Jun 1992, Feb 1996
C Reviewed:          J. Hack, B. Boville, Aug 1992, Apr 1996
C Rewritten for netcdf: J. Truesdale 5/7/97
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <pagrid.h>
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
      integer cnt3(3)               ! array of counts for each dimension
      integer strt3(3)              ! array of starting indices
      integer year1                 ! year of first sst record 
      integer year2                 ! year of last sst record
      integer yr                    ! year number 
      integer ret                   ! return code
      integer n, k                  ! indices
      real caldayloc                ! calendar day includes yr 

      nm = 1
      np = 2

!
! Master does all the work.
      if (masterproc) then
!
!  current date has YYYYYYMMDD.sec
!
          caldayloc = float(mcdate) + float(mcsec)/86400.
!
! Get and check dimension info
!
        lonid = NCDID(ncid_sst, 'lon', ret)
        latid = NCDID(ncid_sst, 'lat', ret)
        timeid = NCDID(ncid_sst, 'time', ret)

        call NCDINQ(ncid_sst, lonid, lonnam, lonsiz, ret)
        if (lonsiz .ne. plon) then
          if(masterproc) write(6,*)'SSTINI: lonsiz=',lonsiz,' must = plon=',plon
          call endrun
        end if

        call NCDINQ(ncid_sst, latid, latnam, latsiz, ret)
        if (latsiz .ne. plat) then
          if(masterproc) write(6,*)'SSTINI: latsiz=',latsiz,' must = plat=',plat
          call endrun
        end if

        call NCDINQ(ncid_sst, timeid, timenam, timesiz, ret)
!
! Check to make sure space allocated for time variables is sufficient
!
        if (timesiz.gt.totsstsz) then
           if(masterproc) write(6,*)'SSTINI:  Allocated space for sst date ',
     $          'variables is insufficient.'
           if(masterproc) write(6,*)'Please increase parameter totsstsz to at least',
     $          timesiz,' and recompile.'
           call endrun
        end if
!
        dateid = NCVID(ncid_sst, 'date', ret)
        secid = NCVID(ncid_sst, 'datesec', ret)
        sstid = NCVID(ncid_sst, 'SST', ret)
        latid = NCVID(ncid_sst, 'lat', ret)
!
! Retrieve entire date and sec variables.
!
        call NCVGT(ncid_sst, dateid, 1, timesiz, date_sst, ret)
        call NCVGT(ncid_sst, secid,  1, timesiz, sec_sst,  ret)

C  Took out 12 month check for cyclic SST data
C  Check for cyclic SST that data is not multiyear.  
C  Warn user if years dont match and strip off year info in date_sst
        if (sstcyc) then
          year1 = int(float(date_sst(1))/10000.) 
          year2 = int(float(date_sst(timesiz))/10000.) 
          yr = mcdate/10000
          if (year1 .ne. year2) then
            if(masterproc) write(6,*) 
     .      'SSTINI: sstcyc = true but sst data has multiple years.'
            if(masterproc) write(6,*) 
     .      '        Set sstcyc = false or change sst to climatology',
     .      ' or single year.'
            call endrun
          else if (year1 .ne. 0 .and. yr .ne. year1) then
            if(masterproc) write(6,*) 'SSTINI: Warning - sst data year does not ',
     .                 'match model.'
            if(masterproc) write(6,*) '        Will continue run because sstcyc = ',
     .                 'true.'

            do n = 1, timesiz
              date_sst(n) = date_sst(n) - 
     .                      int(float(date_sst(n))/10000.)*10000
            enddo

          endif
        endif

        strt3(1) = 1
        strt3(2) = 1
        strt3(3) = 1
        cnt3(1)  = lonsiz
        cnt3(2)  = latsiz
        cnt3(3)  = 1
!
! Special code for interpolation between December and January
!
        if (sstcyc) then
C  Assume cyclic clim SST has a complete year whether monthly or weekly etc
          n = timesiz
          np1 = 1
          yr = mcdate/10000
          cdaysstm = float((yr-1)*10000+date_sst(n)) + 
     &               float(sec_sst(n))/86400.
          cdaysstp = float(yr*10000+date_sst(np1)) + 
     &               float(sec_sst(np1))/86400.
          if (caldayloc.le.cdaysstp .and. caldayloc.gt.cdaysstm) then
            strt3(3) = n
            call NCVGT(ncid_sst,sstid,strt3,cnt3,sstbdy(1,1,nm),ret)
            strt3(3) = np1
            call NCVGT(ncid_sst,sstid,strt3,cnt3,sstbdy(1,1,np),ret)
            goto 10
          end if
        end if
!
! Normal interpolation between consecutive time slices.
!
        do n=1,timesiz-1
          np1 = n + 1
          cdaysstm = float(date_sst(n)) + float(sec_sst(n))/86400.
          cdaysstp = float(date_sst(np1)) + float(sec_sst(np1))/86400.

c  add year onto sst data dates
          if (sstcyc) then
            yr = mcdate/10000
            cdaysstm = cdaysstm + yr*10000
            cdaysstp = cdaysstp + yr*10000
          endif

          if (caldayloc.gt.cdaysstm .and. caldayloc.le.cdaysstp) then
            strt3(3) = n
            call NCVGT(ncid_sst,sstid,strt3,cnt3,sstbdy(1,1,nm),ret)
            strt3(3) = np1
            call NCVGT(ncid_sst,sstid,strt3,cnt3,sstbdy(1,1,np),ret)
            goto 10
          end if
        end do

c  Last check for date after last record in climate file
        if (sstcyc) then
          n = timesiz 
          np1 = 1
          yr = mcdate/10000
          cdaysstm = float(yr*10000+ date_sst(n)) + 
     $               float(sec_sst(n))/86400.
          cdaysstp = float((yr+1)*10000 + date_sst(np1)) + 
     $               float(sec_sst(np1))/86400.
          if (caldayloc.gt.cdaysstm .and. caldayloc.le.cdaysstp) then
            strt3(3) = n
            call NCVGT(ncid_sst,sstid,strt3,cnt3,sstbdy(1,1,nm),ret)
            strt3(3) = np1
            call NCVGT(ncid_sst,sstid,strt3,cnt3,sstbdy(1,1,np),ret)
            goto 10
          end if
        endif

        if(masterproc) write(6,*)'SSTINI: Failed to find dates bracketing ',
     $            'mcdate, mcsec=', mcdate, mcsec
        call endrun
   10   continue
        if(masterproc) write(6,*)'SSTINI: Read sst data for dates ',date_sst(n), 
     &            sec_sst(n),' and ',date_sst(np1),sec_sst(np1)

      end if                 ! end masterproc

      return
      end
