      subroutine tick (nymd, nhms, ndt)

! Input:
      integer ndt                     ! TIME-STEP
! Inpuit/Output:
      integer nymd                    ! CURRENT YYYYMMDD
      integer nhms                    ! CURRENT HHMMSS
      integer incymd

! Revision:   S.-J. Lin Mar 2000

       NSECF(N)   = N/10000*3600 + MOD(N,10000)/100* 60 + MOD(N,100)
       NHMSF(N)   = N/3600*10000 + MOD(N,3600 )/ 60*100 + MOD(N, 60)

       NSEC = NSECF(NHMS) + ndt

       IF (NSEC.GT.86400)  THEN
           DO WHILE (NSEC.GT.86400)
              NSEC = NSEC - 86400
              NYMD = INCYMD (NYMD,1)
           ENDDO
       ENDIF

       IF (NSEC.EQ.86400)  THEN
           NSEC = 0
           NYMD = INCYMD (NYMD,1)
       ENDIF

       IF (NSEC .LT. 0)  THEN
           DO WHILE (NSEC .LT. 0)
               NSEC = 86400 + NSEC
               NYMD = INCYMD (NYMD,-1)
           ENDDO
        ENDIF

          NHMS = NHMSF (NSEC)
      return
      end

      integer FUNCTION INCYMD (NYMD,M)

!  PURPOSE
!     INCYMD:  NYMD CHANGED BY ONE DAY
!     MODYMD:  NYMD CONVERTED TO JULIAN DATE
!  DESCRIPTION OF PARAMETERS
!     NYMD     CURRENT DATE IN YYMMDD FORMAT
!     M        +/- 1 (DAY ADJUSTMENT)

      INTEGER NDPM(12)
      DATA    NDPM /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      logical leap_year
      DATA    NY00     / 1900 /

      NY = NYMD / 10000
      NM = MOD(NYMD,10000) / 100
      ND = MOD(NYMD,100) + M

      IF (ND.EQ.0) THEN
      NM = NM - 1
      IF (NM.EQ.0) THEN
          NM = 12
          NY = NY - 1
      ENDIF
      ND = NDPM(NM)
      IF (NM.EQ.2 .AND. leap_year(NY))  ND = 29
      ENDIF

      IF (ND.EQ.29 .AND. NM.EQ.2 .AND. leap_year(ny))  GO TO 20

      IF (ND.GT.NDPM(NM)) THEN
      ND = 1
      NM = NM + 1
      IF (NM.GT.12) THEN
          NM = 1
          NY = NY + 1
      ENDIF
      ENDIF

   20 CONTINUE
      INCYMD = NY*10000 + NM*100 + ND
      RETURN
      END
