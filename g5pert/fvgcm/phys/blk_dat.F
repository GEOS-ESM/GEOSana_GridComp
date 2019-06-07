      subroutine blk_dat

C This routine initializes the "common block" for lw originally
C initialized by blkdat.F

      implicit none
 
#include <compbl.h>
#include <crdcae.h>
 
C Water vapor narrow band constants for longwave radiation computations
 
      real  doefa(3,4),doefb(4,4),doefc(3,4),doefd(4,4),
     $      doefe(3,4),doeff(6,2),doefg(2,4),doefh(2,4),
     $      doefi(6,2),doefj(3,2),doefk(3,2)

      integer i, j

C H2O EMISSIVITY AND ABSORTIVITY COEFFICIENTS
 
      data doefa/1.01400e+00, 6.41695e-03, 2.85787e-05,
     $           1.01320e+00, 6.86400e-03, 2.96961e-05,
     $           1.02920e+00, 1.01680e-02, 5.30226e-05,
     $           1.02743e+00, 9.85113e-03, 5.00233e-05/
 
      data doefb/8.85675e+00,-3.51620e-02, 2.38653e-04,-1.71439e-06,
     $           5.73841e+00,-1.91919e-02, 1.65993e-04,-1.54665e-06,
     $           6.64034e+00, 1.56651e-02,-9.73357e-05, 0.0,
     $           7.09281e+00, 1.40056e-02,-1.15774e-04, 0.0/
 
      data doefc/9.90127e-01, 1.22475e-03, 4.90135e-06,
     $           9.89753e-01, 1.97081e-03, 3.42046e-06,
     $           9.75230e-01, 1.03341e-03, 0.0,
     $           9.77366e-01, 8.60014e-04, 0.0/
 
      data doefd/7.03047e-01,-2.63501e-03,-1.57023e-06, 0.0,
     $           5.29269e-01,-3.14754e-03, 4.39595e-06, 0.0,
     $           7.88193e-02, 1.31290e-03, 4.25827e-06,-1.23982e-08,
     $           1.62744e-01, 2.22847e-03, 2.60102e-06,-4.30133e-08/
 
      data doefe/3.93137e-02,-4.34341e-05, 3.74545e-07,
     $           3.67785e-02,-3.10794e-05, 2.94436e-07,
     $           7.42500e-02, 3.97397e-05, 0.0,
     $           7.52859e-02, 4.18073e-05, 0.0/
 
      data doeff/2.2037 e-01,1.39719e-03,-7.32011e-06,
     $          -1.40262e-08,2.13638e-10,-2.35955e-13,
     $           3.07431e-01,8.27225e-04,-1.30067e-05,
     $           3.49847e-08,2.07835e-10,-1.98937e-12/
 
      data doefg/9.04489e+00,-9.56499e-03,
     $           1.80898e+01,-1.91300e-02,
     $           8.72239e+00,-9.53359e-03,
     $           1.74448e+01,-1.90672e-02/
 
      data doefh/5.46557e+01,-7.30387e-02,
     $           1.09311e+02,-1.46077e-01,
     $           5.11479e+01,-6.82615e-02,
     $           1.02296e+02,-1.36523e-01/
 
      data doefi/3.31654e-01,-2.86103e-04,-7.87860e-06,
     $           5.88187e-08,-1.25340e-10,-1.37731e-12,
     $           3.14365e-01,-1.33872e-03,-2.15585e-06,
     $           6.07798e-08,-3.45612e-10,-9.34139e-15/
 
      data doefj/2.82096e-02,2.47836e-04,1.16904e-06,
     $           9.27379e-02,8.04454e-04,6.88844e-06/
 
      data doefk/2.48852e-01,2.09667e-03,2.60377e-06,
     $           1.03594e+00,6.58620e-03,4.04456e-06/

c /PBL/
c     data vk, fak, fakn / 0.40, 8.5 , 7.2/

      vk   = 0.40       ! Von Karman's constant
      fak  = 8.5        ! Constant in surface temperature excess
      fakn = 7.2        ! Constant in turbulent prandtl number
 
c     data betam, betas, betah, sffrac / 15.0, 5.0, 15.0, 0.1 /

      betam  = 15.0     ! Constant in wind gradient expression 
      betas  = 5.0      ! Constant in surface layer gradient expression
      betah  = 15.0     ! Constant in temperature gradient expression
      sffrac = 0.1      ! Surface layer fraction of boundary layer
      ricr   = 0.3      ! Critical richardson number

C /crdcae/
C Narrow band data for H2O
C 200CM data for 800-1000 CM-1 and 1000-1200 CM-1.
C
      realk(1) = 0.18967069430426e-04
      realk(2) = 0.70172244841851e-04
 
      st(1) = 0.31930234492350e-03
      st(2) = 0.97907319939060e-03 

      a1(1) = 0.28775403075736e-01
      a1(2) = 0.23236701470511e-01

      a2(1) = -0.57966222388131e-04
      a2(2) = -0.95105504388411e-04

      b1(1) =  0.29927771523756e-01
      b1(2) =  0.21737073577293e-01

      b2(1) = -0.86322071248593e-04
      b2(2) = -0.78543550629536e-04

      do j=1,4
         do i=1,3
            coefa(i,j) = doefa(i,j)
         enddo
      enddo
 
      do j=1,4
         do i=1,4
            coefb(i,j) = doefb(i,j)
         enddo
      enddo
 
 
      do j=1,4
         do i=1,3
            coefc(i,j) = doefc(i,j)
         enddo
      enddo
 
      do j=1,4
         do i=1,4
            coefd(i,j) = doefd(i,j)
         enddo
      enddo
 
      do j=1,4
         do i=1,3
            coefe(i,j) = doefe(i,j)
         enddo
      enddo
 
      do j=1,2
         do i=1,6
            coeff(i,j) = doeff(i,j)
         enddo
      enddo
 
      do j=1,4
         do i=1,2
            coefg(i,j) = doefg(i,j)
         enddo
      enddo
 
      do j=1,4
         do i=1,2
            coefh(i,j) = doefh(i,j)
         enddo
      enddo
 
      do j=1,2
         do i=1,6
            coefi(i,j) = doefi(i,j)
         enddo
      enddo
 
      do j=1,2
         do i=1,3
            coefj(i,j) = doefj(i,j)
         enddo
      enddo
 
      do j=1,2
         do i=1,3
            coefk(i,j) = doefk(i,j)
         enddo
      enddo

      return
      end
