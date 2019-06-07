MODULE qsat_util

!This module contains the subroutines used to generate the lookup table for computing saturation vapour pressure.
!Variables designated as real*8 should be kept at least as so, even if ESTBLX is required in real*4.

IMPLICIT NONE

PRIVATE
PUBLIC ESINIT

!Saturation vapor pressure table ESTBLX initialization.
!To be defined for entire module:
integer, parameter :: DEGSUBS    =  100
real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

CONTAINS

subroutine ESINIT(ESTBLX)

 IMPLICIT NONE

 !OUTPUT
 real(8), dimension(TABLESIZE) :: ESTBLX

 !LOCALS
 real(8), parameter :: ZEROC = 273.16, TMIX = -20.0

 real(8), dimension(TABLESIZE) :: ESTBLE, ESTBLW

 integer :: I
 real(8)    :: T, DELTA_T

 DELTA_T = 1.0/DEGSUBS

 do I=1,TABLESIZE

    T = (I-1)*DELTA_T + TMINTBL

    if(T>ZEROC) then
       call QSATLQU0(ESTBLE(I),T)
    else
       call QSATICE0(ESTBLE(I),T)
    end if

    call QSATLQU0(ESTBLW(I),T)

    T = T-ZEROC
    if(T>=TMIX .and. T<0.0) then
       ESTBLX(I) = ( T/TMIX )*( ESTBLE(I) - ESTBLW(I) ) + ESTBLW(I)
    else
       ESTBLX(I) = ESTBLE(I)
    end if

 end do

 end subroutine ESINIT


subroutine QSATLQU0(QS,TL)
!SUPERSATURATED AS LIQUID

 IMPLICIT NONE

 !INPUTS
 real(8) :: TL!, TMAXTBL

 !OUTPUTS
 real(8) :: QS

 !LOCALS
 real(8), parameter :: ZEROC   = 273.16
 real(8), parameter :: TMINLQU = ZEROC - 40.0

 real*8,  parameter :: B6 = 6.136820929E-11*100.0
 real*8,  parameter :: B5 = 2.034080948E-8 *100.0
 real*8,  parameter :: B4 = 3.031240396E-6 *100.0
 real*8,  parameter :: B3 = 2.650648471E-4 *100.0
 real*8,  parameter :: B2 = 1.428945805E-2 *100.0
 real*8,  parameter :: B1 = 4.436518521E-1 *100.0
 real*8,  parameter :: B0 = 6.107799961E+0 *100.0

 real(8) :: TX, EX, TI, TT

 TX = TL

 if    (TX<TMINLQU) then
    TI = TMINLQU
 elseif(TX>TMAXTBL) then
    TI = TMAXTBL
 else
    TI = TX
 end if

 TT = TI-ZEROC  !Starr polynomial fit
 EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)

 TL = TX
 QS = EX

 return

end subroutine QSATLQU0


subroutine QSATICE0(QS,TL)
!SUPERSATURATED AS ICE

 IMPLICIT NONE   
      
 !INPUTS
 real(8) :: TL

 !OUTPUTS
 real(8) :: QS

 !LOCALS
 real(8), parameter :: ZEROC = 273.16, TMINSTR = -95.0
 real(8), parameter :: TMINICE = ZEROC + TMINSTR

 real(8), parameter :: TSTARR1 = -75.0, TSTARR2 = -65.0, TSTARR3 = -50.0,  TSTARR4 = -40.0

 real*8,  parameter :: BI6= 1.838826904E-10*100.0
 real*8,  parameter :: BI5= 4.838803174E-8 *100.0
 real*8,  parameter :: BI4= 5.824720280E-6 *100.0
 real*8,  parameter :: BI3= 4.176223716E-4 *100.0
 real*8,  parameter :: BI2= 1.886013408E-2 *100.0
 real*8,  parameter :: BI1= 5.034698970E-1 *100.0
 real*8,  parameter :: BI0= 6.109177956E+0 *100.0
 real*8,  parameter :: S16= 0.516000335E-11*100.0
 real*8,  parameter :: S15= 0.276961083E-8 *100.0
 real*8,  parameter :: S14= 0.623439266E-6 *100.0
 real*8,  parameter :: S13= 0.754129933E-4 *100.0
 real*8,  parameter :: S12= 0.517609116E-2 *100.0
 real*8,  parameter :: S11= 0.191372282E+0 *100.0
 real*8,  parameter :: S10= 0.298152339E+1 *100.0
 real*8,  parameter :: S26= 0.314296723E-10*100.0
 real*8,  parameter :: S25= 0.132243858E-7 *100.0
 real*8,  parameter :: S24= 0.236279781E-5 *100.0
 real*8,  parameter :: S23= 0.230325039E-3 *100.0
 real*8,  parameter :: S22= 0.129690326E-1 *100.0
 real*8,  parameter :: S21= 0.401390832E+0 *100.0
 real*8,  parameter :: S20= 0.535098336E+1 *100.0

 real(8) :: TX, TI, TT, W, EX

 TX = TL

 if (TX<TMINICE) then
    TI = TMINICE
 elseif(TX>ZEROC  ) then
    TI = ZEROC
 else
    TI = TX
 end if

 TT = TI - ZEROC
 if (TT < TSTARR1 ) then
     EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
 elseif(TT >= TSTARR1 .and. TT < TSTARR2) then
     W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
     EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
              + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
 elseif(TT >= TSTARR2 .and. TT < TSTARR3) then
     EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
 elseif(TT >= TSTARR3 .and. TT < TSTARR4) then
     W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
     EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
              + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
 else
     EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
 endif

 QS = EX

 return
 
end subroutine QSATICE0

end MODULE qsat_util
