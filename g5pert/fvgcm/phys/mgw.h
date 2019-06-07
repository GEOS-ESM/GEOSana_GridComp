c
c $Id$
c $Author$
c
c These constants are documented in the gravity wave section of NCAR
c TN-420 "Description of the NCAR Community Climate Model (CCM3)",
c  Kiehl et al. (1996)
c
c
#ifndef GWPARMS_SET
#define GWPARMS_SET
#define PGWV	4
#define FRACLDV 0.0
#define TAULDV  0.0
#define TAUSCAL	.001
#define TAUNOR	0.75
#define TAUSOU	1.2
#define TAUBGND 6.4
#define TAUEQ	1.
#define TAUEQW	10.
#define OROFC2	.5
#define OROE	1.
#define EFFGW   .125
#define MXRANGE	.001
#define MXASYM	.1
#endif
c+
c common variables for multiple gravity wave parameterization
c-
      integer pgwv 
      parameter (pgwv = PGWV)  ! number of waves allowed

      integer
     $     kbotbg, kbotoro,     ! interface of gwd source
     $     ktopbg, ktoporo      ! top interface of gwd region

      real
     $     alpha(0:plev),       ! newtonian cooling coefficients
     $     c(-pgwv:pgwv),       ! list of wave phase speeds
     $     cp,                  ! specific heat of dry air (constant p)
     $     cpvir,               ! specific humidity factor for specific heat
     $     dback,               ! background diffusivity
     $     efcncy,              ! "efficiency" factor
     $     efkw,                ! efcncy * kwv
     $     effgw,               ! tendency efficiency
     $     fracldv,             ! fraction of stress deposited in low level region
     $     g,                   ! acceleration of gravity
     $     kwv                  ! effective horizontal wave number
      real
     $     mxasym,              ! max asymmetry between tau(c) and tau(-c)
     $     mxrange,             ! max range of tau for all c
     $     n2min,               ! min value of bouyancy frequency
     $     orofc2,              ! critical froude number
     $     oroeko2,             ! e*k/2, tunable parameter
     $     orohmin,             ! min sdv of height for orographic waves
     $     orovmin,             ! min wind speed for orographic waves
     $     r,                   ! gas constant for dry air
     $     rog,                 ! r / g
     $     taubgnd,             ! background source strength (/tauscal)
     $     taumin,              ! minimum (nonzero) stress
     $     tauscal,             ! scale factor for background stress source
     $     tndmin,              ! minimum wind tendency
     $     tndmax,              ! maximum wind tendency
     $     ubmc2mn,             ! min (u-c)**2
     $     zldvcon              ! constant for determining zldv from tau0

      common /mgwdcom/
     $     alpha, 
     $     c, cp, cpvir,
     $     dback,
     $     efcncy, efkw,
     $     effgw,
     $     fracldv,
     $     g,
     $     kwv,
     $     mxasym, mxrange,
     $     n2min,     
     $     oroeko2, orofc2, orohmin, orovmin,
     $     r,
     $     rog,
     $     taubgnd, taumin, tauscal,
     $     tndmin, tndmax,
     $     ubmc2mn,
     $     zldvcon,
     $     kbotbg, kbotoro, 
     $     ktopbg, ktoporo  
#ifdef LOGENTRY 
$Log$
Revision 1.1.1.1.2.1  2001/11/26 22:59:05  jchern
MLP/MPI version of fvccm

Revision 1.3  2001/10/08 18:37:58  sjlin
*** empty log message ***

#endif
