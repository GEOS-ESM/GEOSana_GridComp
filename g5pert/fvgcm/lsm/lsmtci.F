#include <preproc.h>
      subroutine lsmtci (veg, wt, beglatspmd)

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>

* ------------------------ code history ------------------------------
* source file:       lsmtci.F
* purpose:           set up multiprocessing and initialize time 
*                    invariant variables as subgrid vectors of 
*                    length kpt (for lsmtcdyn common block)
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* --------------------------------------------------------------------

* ------------------------ input variables ---------------------------
      integer veg(lsmlon,lsmlat,msub)  !vegetation type
      real wt(lsmlon,lsmlat,msub)      !subgrid weights
      integer beglatspmd    !beg lat on proc, degenerates to 1 for nonspmd
* --------------------------------------------------------------------

* ------------------------ common block variables --------------------
#include <dynmem.h>
#include <lsmtc.h>
#include <lsmtcdyn.h>
#include <soicon.h>
#include <vegcon.h>
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      integer i,j,k,m  !indices
      integer lvec     !number of land points to process per call to lsm
      integer mpt      !max points per call to lsm (lvec + remainder (kpt/lvec))
      real r1          !root fraction at top of soil layer
      real r2          !root fraction at bottom of soil layer
* --------------------------------------------------------------------
* set up multi-processing: the "big" vectors of [kpt] subgrid points 
* are processed as [numlv]-1 "little" vectors of [lvec] points and 
* one "little" vector of [mpt] points. [lvec] does not necessarily 
* equal [mpt] because [kpt] may not be evenly divisible by [numlv]. 
* the starting location in the "big" vector is [begkpt]. the first 
* "little" vector is longest so it gets done first. note that 
* [mpt] is the array dimension within the vector lsm code.

      lvec  = numks/numlv 
      mpt   = lvec+numks-numlv*lvec

      begkpt(1) = begkptspmd(beglatspmd)
      numkpt(1) = mpt 
      do i = 2, numlv
         begkpt(i) = begkpt(i-1) + numkpt(i-1)
         numkpt(i) = lvec
      end do

      if (masterproc) then
#if ( defined SPMD )
         write (6,*)'multi-processing:'
         write (6,*)'the big vector of ',kpt,'points will be processed'
         write (6,*)'as little vectors of land points.  The number of'
         write (6,*)'land points in each vector is determined by the'
         write (6,*)'atmospheric latitudes per forked process'
#else
         write (6,*) 'multi-processing:'
         write (6,*) 'the big vector of ',kpt,
     &        ' points will be processed as ',numlv,' little vectors'
         write (6,1002) numlv-(numlv-1),mpt
         write (6,1002) numlv-1,lvec
1002  format (5x,i4,' little vectors with length ',i7)
         write (6,*)
#endif
      endif
* ivt, isc, ist, sand, and clay. make sure k = kpt 
 
      k = 0                                
      do j = 1, lsmlat
         do i = 1, numlon(j)
            if (surf2d(i,j) .gt. 0) then                 !land point
               do m = 1, msub                            !subgrid points
                  if (wt(i,j,m) .gt. 0.) then            !valid subgrid point
                     k = k + 1                           !subgrid index
                     ivt(k) = veg(i,j,m)
                     isc(k) = soic2d(i,j)

                     if (m.eq.4) then                    !deep lake, from pctlak
                        ist(k) = istdlak
                        sand(k) = 0.
                        clay(k) = 0.
                     else if (m.eq.5) then               !wetland, from pctwet
                        ist(k) = istwet
                        sand(k) = 0.
                        clay(k) = 0.
                     else                                !soil type using surf2d
                        sand(k) = sand2d(i,j)
                        clay(k) = clay2d(i,j)
                        if (surf2d(i,j).eq.1) then       !glacier
                           ist(k) = istice
                        else if (surf2d(i,j).eq.27) then !wetland
                           ist(k) = istwet
                        else if (surf2d(i,j).eq.28) then !wetland
                           ist(k) = istwet
                        else                             !soil
                           ist(k) = 1
                        end if
                     end if
                  end if
               end do
            end if
         end do
      end do

      if (k .ne. kpt) then
         write (6,*) 'LSMTCI fatal error: discrepancy in subgrid points'
         write (6,*) 'expecting ',kpt,' points but found ',k
         call stoprun
      end if

* latitudes and longitudes

      do k = 1, kpt
         lati(k) = latixy(ixy(klnd(k)),jxy(klnd(k)))
         long(k) = longxy(ixy(klnd(k)),jxy(klnd(k)))
      end do

* soil/lake layer thickness, depths, and relative root fraction

      do k = 1, kpt
         if (ist(k) .eq. istdlak) then                !deep lake
            dzsoi(1,k) = 1.
            dzsoi(2,k) = 2.
            dzsoi(3,k) = 4.
            dzsoi(4,k) = 8.
            dzsoi(5,k) = 15.
            dzsoi(6,k) = 20.
         else if (ist(k) .eq. istslak) then           !shallow lake
            dzsoi(1,k) = 0.5
            dzsoi(2,k) = 1.0
            dzsoi(3,k) = 1.5
            dzsoi(4,k) = 2.0
            dzsoi(5,k) = 2.5
            dzsoi(6,k) = 2.5
         else                                         !soil, ice, wetland
            dzsoi(1,k) = 0.1
            dzsoi(2,k) = 0.2
            dzsoi(3,k) = 0.4
            dzsoi(4,k) = 0.8
            dzsoi(5,k) = 1.6
            dzsoi(6,k) = 3.2
         end if
      end do

      do j = 1, msl
         do k = 1, kpt
            if (j. eq. 1) then
               zsoi(j,k) = 0.5*dzsoi(j,k)
            else
               zsoi(j,k) = zsoi(j-1,k) + 0.5 * (dzsoi(j-1,k)+dzsoi(j,k))
            end if
            r1 = 1. - rdp(ivt(k)) ** (100.*(zsoi(j,k)-0.5*dzsoi(j,k)))
            r2 = 1. - rdp(ivt(k)) ** (100.*(zsoi(j,k)+0.5*dzsoi(j,k)))
            root(j,k) = r2 - r1
         end do
      end do

* soil thermal and hydraulic properties 

      do k = 1, kpt
         if (ist(k) .eq. 1) then                      !soil
            watsat(k) = 0.489 - 0.00126*sand(k) 
            hksat(k) = 0.0070556 * ( 10.**(-0.884+0.0153*sand(k)) )
            smpsat(k) = -10. * ( 10.**(1.88-0.0131*sand(k)) )
            bch(k) = 2.91 + 0.159*clay(k)
            watdry(k) = watsat(k) * (-316230./smpsat(k)) ** (-1./bch(k))
            watopt(k) = watsat(k) * (-158490./smpsat(k)) ** (-1./bch(k))
            tksol(k) =(8.80*sand(k)+2.92*clay(k)) / (sand(k)+clay(k))
            tkdry(k) =0.15
            csol(k)=(2.128*sand(k)+2.385*clay(k))/(sand(k)+clay(k))*1.e6
         else                                         !ice, lakes, wetlands
            watsat(k) = 1.
            hksat(k) = 999.
            smpsat(k) = 999.
            bch(k) = 999.
            watdry(k) = 0.
            watopt(k) = 1.
            tksol(k) = 999.
            tkdry(k) = 999.
            csol(k) = 999.
         end if
      end do

      return
      end
 
