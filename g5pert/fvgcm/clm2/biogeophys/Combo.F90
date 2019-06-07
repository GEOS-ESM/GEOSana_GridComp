#include <misc.h>
#include <preproc.h>

subroutine Combo     (dz,    wliq,  wice,  t, dz2, &
                      wliq2, wice2, t2     ) 

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Combines two elements and returns the following combined
! variables: dz, t, wliq, wice. 
!
! Method:
! The combined temperature is based on the equation:
! the sum of the enthalpies of the two elements =  
! that of the combined element.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use clm_varcon,  only : cpice, cpliq, tfrz, hfus         
  implicit none

!----Arguments----------------------------------------------------------

  real(r8), intent(in) :: dz2   ! nodal thickness of 2 elements being combined [m]
  real(r8), intent(in) :: wliq2 ! liquid water of element 2 [kg/m2]
  real(r8), intent(in) :: wice2 ! ice of element 2 [kg/m2]
  real(r8), intent(in) :: t2    ! nodal temperature of element 2 [K]

  real(r8), intent(inout) :: dz   ! nodal thickness of 1 elements being combined [m]
  real(r8), intent(inout) :: wliq ! liquid water of element 1
  real(r8), intent(inout) :: wice ! ice of element 1 [kg/m2]
  real(r8), intent(inout) :: t    ! nodel temperature of elment 1 [K]

!----Local Variables----------------------------------------------------

  real(r8) dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
  real(r8) wliqc ! Combined liquid water [kg/m2]
  real(r8) wicec ! Combined ice [kg/m2]
  real(r8) tc    ! Combined node temperature [K]
  real(r8) h     ! enthalpy of element 1 [J/m2]
  real(r8) h2    ! enthalpy of element 2 [J/m2]
  real(r8) hc    ! temporary

!----End Variable List--------------------------------------------------

  dzc = dz+dz2
  wicec = (wice+wice2)
  wliqc = (wliq+wliq2)
  h =(cpice*wice+cpliq*wliq)* &
       (t-tfrz)+hfus*wliq
  h2=(cpice*wice2+cpliq*wliq2)* &
       (t2-tfrz)+hfus*wliq2

  hc = h + h2
  if(hc < 0.)then
     tc = tfrz + hc/(cpice* &
          wicec+cpliq*wliqc)
  else if(hc.le.hfus*wliqc)then
     tc = tfrz
  else
     tc = tfrz + (hc - hfus*wliqc)/ &
          (cpice*wicec+cpliq*wliqc)
  endif

  dz = dzc
  wice = wicec 
  wliq = wliqc
  t = tc

end subroutine Combo
