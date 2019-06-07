      subroutine interp_time ( nymd ,sec , 
     .                         nymd1,sec1, nymd2,sec2, fac1,fac2 )
C***********************************************************************        
C                                                                               
C  PURPOSE:
C  ========
C    Compute interpolation factors, fac1 & fac2, to be used in the
C    calculation of the instantanious boundary conditions, ie:
C
C               q(i,j) = fac1*q1(i,j) + fac2*q2(i,j)
C    where:
C               q(i,j) => Boundary Data valid    at (nymd  , nhms )
C              q1(i,j) => Boundary Data centered at (nymd1 , nhms1)
C              q2(i,j) => Boundary Data centered at (nymd2 , nhms2)
C                                                                               
C  INPUT:                                                                       
C  ======                                                                       
C    nymd     : Date (yymmdd) of Current Timestep
C    sec      : Time (seconds) of Current Timestep
C    nymd1    : Date (yymmdd) of Boundary Data 1
C    sec1     : Time (seconds) of Boundary Data 1
C    nymd2    : Date (yymmdd) of Boundary Data 2
C    sec2     : Time (seconds) of Boundary Data 2
C                                                                               
C  OUTPUT:                                                                      
C  =======                                                                      
C    fac1     : Interpolation factor for Boundary Data 1
C    fac2     : Interpolation factor for Boundary Data 2
C                                                                               
C                                                                               
C***********************************************************************        
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *        
C***********************************************************************        
                                                                                
      INTEGER  YEAR , MONTH , DAY , SEC                                                    
      INTEGER  YEAR1, MONTH1, DAY1, SEC1
      INTEGER  YEAR2, MONTH2, DAY2, SEC2
      INTEGER  YEARMIN, YEARMAX, JDAY, I 

      real     fac1, fac2
      real     time, time1, time2
                                                                                
      REAL       MNDY(12), MNDYL(12)

      logical leap_year
      external leap_year
                                                                                
      DATA MNDYL /0,31,60,91,121,152,182,213,244,274,305,335/
      DATA MNDY  /0,31,59,90,120,151,181,212,243,273,304,334/
                                                                                
C***********************************************************************        
C*                         SET TIME BOUNDARIES                         *        
C***********************************************************************        
                                                                                
      YEAR   = NYMD / 10000                                                      
      MONTH  = MOD(NYMD,10000) / 100                                             
      DAY    = MOD(NYMD,100)                                                     
                                                                                
      YEAR1  = NYMD1 / 10000                                                      
      MONTH1 = MOD(NYMD1,10000) / 100                                             
      DAY1   = MOD(NYMD1,100)                                                     
                                                                                
      YEAR2  = NYMD2 / 10000                                                      
      MONTH2 = MOD(NYMD2,10000) / 100                                             
      DAY2   = MOD(NYMD2,100)                                                     
                                                                                
C***********************************************************************        
C*                    COMPUTE JULIAN DAYS
C***********************************************************************        
                                                                                
      yearmin = min(min(year,year1),year2)
      yearmax = max(max(year,year1),year2)

      jday = 0
      do i = yearmin, yearmax
        if (year .eq. i) then
          if (leap_year(year)) then
            time  = DAY+ MNDYL(MONTH) + float(sec)/86400. + jday
          else
            time  = DAY+ MNDY(MONTH) + float(sec)/86400. + jday
          endif
        endif
        if (year1 .eq. i) then
          if (leap_year(year1)) then
            time1 = DAY1+ MNDYL(MONTH1) + float(sec1)/86400. + jday
          else
            time1 = DAY1+ MNDY(MONTH1) + float(sec1)/86400. + jday
          endif
        endif
        if (year2 .eq. i) then
          if (leap_year(year2)) then
            time2 = DAY2+ MNDYL(MONTH2) + float(sec2)/86400. + jday
          else
            time2 = DAY2+ MNDY(MONTH2) + float(sec2)/86400. + jday
          endif
        endif

        if (leap_year(i)) then
          jday = jday + 366
        else
          jday = jday + 365
        endif

      enddo

      fac1  = (time2-time)/(time2-time1)
      fac2  = (time-time1)/(time2-time1)

      RETURN                                                                    
      END                                                                       

