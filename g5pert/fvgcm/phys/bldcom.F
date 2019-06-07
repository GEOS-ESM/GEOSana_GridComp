#include <params.h>
      subroutine bldcom(nstepb  ,itape)
C-----------------------------------------------------------------------
C
C Build Mass Store comment field for current model history file.
C The comment field holds date and time information and is stored
C in common block /commsc/ variable msscom.
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath, April 1992
C Standardized:      L. Bath, June 1992
C                    T. Acker, March 1996
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
C-----------------------------------------------------------------------
#include <pagrid.h>
C------------------------------Commons----------------------------------
#include <comhed.h>
C-----------------------------------------------------------------------
#include <comhst.h>
C-----------------------------------------------------------------------
#include <commss.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer  nstepb        ! Iteration of history file data
      integer  itape         ! History file number
C
C---------------------------Local variables-----------------------------
C
      real hour_s           ! hour corresponding to first "current" 
C                           !    seconds of date for history file
      real hour_e           ! hour corres. to current seconds, ncsec
      real frac_s           ! fraction of a day corres. to first 
C                           !    "current" seconds of date for hist. file
      real frac_e           ! fraction of a day corres. to current sec.
      integer idate         ! date of first "current" date on hist. file 
      integer mbday         ! modified base day for monthly-ave. h-file
      integer mbmnth        ! base month for monthly-ave hist. file
      integer mbyear        ! base year for monthly-ave hist. file
      integer ndm(12)       ! no. of days in each month (jan-dec)
C
C--------------------------Data Statements------------------------------
C
      save ndm
      data ndm/31,28,31,30,31,30,31,31,30,31,30,31/
C
C-----------------------------------------------------------------------
C
      hour_s = float(ncsecf(itape))/3600.
      hour_e = float(ncsec)/3600.
      if (ninavg(1).eq.'Q'.and.nstepb.gt.1.and.itape.eq.1) then
C
C  Monthly-averaged history file
C
C  We have to add one day to ncdatf to make the monthly average
C  output look like the Buja monthly average processor ie yymm02-yymm01
C
         idate =  ncdatf(itape)
         mbyear = idate/10000
         mbmnth = mod(idate,10000)/100
         mbday = mod(idate,100)
         mbday = mbday + 1
         if (mbday.gt.ndm(mbmnth)) then
           mbday = 1
           mbmnth = mbmnth + 1
         end if
         idate = mbyear*10000 + mbmnth*100 + mbday
         write(msscom,100)hour_s,idate,hour_e,ncdate
  100    format('MONTHLY-AVERAGE FILE - DATES: ',f6.3,'Z ',i8.8,' -',
     $          f6.3,'Z ',i8.8)
      else
C
C History file not monthly-averaged
C
         frac_s = float(nscurf(itape))/86400.
         frac_e = float(nscur)/86400.
         write(msscom,200)ndcurf(itape),frac_s,ndcur,frac_e,
     $                    hour_s,ncdatf(itape),hour_e,ncdate
  200    format('DAYS: ',i8.8,f4.3,'-',i8.8,f4.3,
     $          ' DATES:',f6.3,'Z ',i8.8,' -',f6.3,'Z ',i8.8)
      end if 
C
      return
C
      end

 
