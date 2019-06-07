!-----------------------------------------------------------------------
!BOP
! !ROUTINE: highp2 --- Compute finite-volume mean pressure forces
!                     using more accurate contour decomposition
!
! !INTERFACE:
      subroutine highp2(pk,   wz,    pzx,  pzy,         &
                       pzz,   im,    jm,   km,          &
                       jfirst,  jlast,  nx )

      implicit none

! !INPUT PARAMETERS:

      integer im, jm, km, jfirst, jlast
      integer nx                        ! # of pieces in longitude direction

      real pk(im,jfirst-1:jlast+1,km+1)
      real wz(im,jfirst-1:jlast+1,km+1) 

! !OUTPUT PARAMETERS
      real pzx(im,jfirst-1:jlast  ,km+1)  
      real pzy(im,jfirst  :jlast+1,km+1) 
      real pzz(im,jfirst-1:jlast+1,km)

! ! !DESCRIPTION:
! Note: this routine needs to be further optimized for speed.
!
! !REVISION HISTORY:
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
      integer i, j, k
      integer ixj, jp, it, i1, i2

      real pka(km+1)
      real pkb(km+1)
      real wza(km+1)
      real wzb(km+1)
      real pf(km+1)
      real pza(km)
      real pzb(km)
      real dka(km)
      real dkb(km)

      real    tmp
      integer im1

      integer js1g1
      integer js2g0
      integer jn2g0
      integer jn1g1

      js1g1  = max(1,jfirst-1)

! E-W PG
      jn2g0  = min(jm-1,jlast)

! N-S PG
      js2g0  = max(2,jfirst)
      jn1g1  = min(jm,jlast+1)

      it = im / nx
      jp = nx * ( jn1g1 - js1g1  + 1 )

!$omp parallel do private(i1, i2, ixj, i, j, k, tmp, pka, pkb, pf, pza, pzb, wza, wzb, dka, dkb, im1)

      do 1000 ixj=1, jp
! ***
! E-W
! ***
         j  = js1g1 + (ixj-1)/nx
         i1 = 1 + it * mod(ixj-1, nx)
         i2 = i1 + it - 1

        do i=i1,i2

        if( i == i1 ) then

           if ( i == 1) then
              im1 = im
           else
              im1 = i-1
           endif 

        do k=1,km+1
           pka(k) =  pk(im1,j,k)
           wza(k) =  wz(im1,j,k)
        enddo

        do k=1,km
           dka(k) =  pka(k+1) - pka(k)
           pza(k) = 0.5*(wza(k)-wza(k+1))*(pka(k+1)+pka(k))
        enddo

        else

        do k=1,km+1
           pka(k) = pkb(k)
           wza(k) = wzb(k)
        enddo

        do k=1,km
           pza(k) = pzb(k)
           dka(k) = dkb(k)
        enddo

        endif

        do k=1,km+1
           pkb(k) = pk(i,j,k)
           wzb(k) = wz(i,j,k)
        enddo

        do k=1,km
           dkb(k) = pkb(k+1) - pkb(k)
           tmp    = (wzb(k)-wzb(k+1))*(pkb(k+1)+pkb(k))
           pzb(k) = 0.5*tmp
           pzz(i,j,k) = tmp
        enddo

        if ( j <= jn2g0 ) then
        call pxt3 (km,   pf,   pka,  pkb,  wza,  wzb,    &
                   pza,  pzb,  dka,  dkb )

        do k=1,km+1
           pzx(i,j,k) = pf(k)
        enddo
        endif

! ***
! N-S
! ***

! N <-- E at (i,j)

        if ( j >= js2g0 ) then
        do k=1,km+1
           pka(k) = pk(i,j-1,k)
           wza(k) = wz(i,j-1,k)
        enddo

        do k=1,km
           dka(k) = pka(k+1) - pka(k)
           pza(k) = 0.5*(wza(k)-wza(k+1))*(pka(k+1)+pka(k)) 
        enddo

! need pzy(js2g0:jn2g1)
        call pxt3 (km,    pf,   pka,  pkb,  wza,  wzb,     &
                   pza,  pzb,   dka,  dkb )

        do k=1,km+1
           pzy(i,j,k) = pf(k)
        enddo
        endif

      enddo                  ! i-loop
1000  continue

      return
      end

      subroutine pxt3(km,  pf,  pka, pkb, wza, wzb, pza, pzb, dka, dkb)
 
! Compute pressure forcing along "horizontal" coordinate surface
! for a single column

      implicit none

! Input
      integer km
 
      real pka(km+1)     ! p**kappa at the west side of the box
      real pkb(km+1)

      real wza(km+1)
      real wzb(km+1)

      real pza(km)
      real pzb(km)

      real dka(km)
      real dkb(km)

! Output
      real  pf(km+1)

! Local
      integer k
      integer ka
      integer kb
      integer kk
      integer kt
      real wa
      real wb
      real tmp
      logical found

      pf(1) = -(pka(1) + pkb(1))*(wzb(1) - wza(1))

         ka = km/2
         kb = km/2

      do 1000 k=2,km+1

          wa = 0.
          wb = 0.

      if( wzb(k) > wza(k) ) then

!
!                     * ---->  wzb(k)
!                    /
!                   /
!                  /
!                 /
!                /
!               /
!              /
!             /
!            /
!           /
!          /
!         /
!        /
!       /
!      /
!     /
!    * ----> wza(k)
!
! Compute wa (along left edge)
!
! Top

      if(wzb(k) > wza(1)) then 
         wa = 0.5*(pkb(k)+pka(k))*(wzb(k)-wza(k))
      else

            kk = k
            found = .false.
         do while (.not. found )
            if(wzb(k) <= wza(kk-1) .and.          &
               wzb(k) >= wza(kk)       )  then
!
! find p**cappa at left edge (side A) at height wzb(k)
!
               tmp = pka(kk-1) + dka(kk-1) *       &
                    (wza(kk-1)-wzb(k))/(wza(kk-1)-wza(kk))
               wa = wa + 0.5*(wzb(k)-wza(kk))*(pka(kk)+tmp)
               found = .true.
            else
               wa  = wa + pza(kk-1)
               kk = kk - 1
            endif
         enddo
      endif

! Compute wb

         if( k .ne. (km+1)) then
            do kk=k,km
               if( wza(k) <= wzb(kk)   .and.             &
                   wza(k) >= wzb(kk+1)      ) then
                   tmp = pkb(kk) + dkb(kk)*(wzb(kk)-wza(k))   &
                                 / (wzb(kk)-wzb(kk+1))
                   wb = wb + 0.5*(tmp+pkb(kk))*(wzb(kk)-wza(k))
                   goto 222
               else
                  wb = wb + pzb(kk)
               endif
            enddo
222      continue
         endif        ! k safety check

         if( wza(k) < wzb(km+1) ) then
! Integrate down along the left edge!!

            do kk=kb,k-1
            if( wzb(km+1) <= wza(kk)   .and.                      &
                wzb(km+1) >= wza(kk+1)      ) then
                tmp = pka(kk) + dka(kk)*(wza(kk)-wzb(km+1))       &
                              / (wza(kk)-wza(kk+1))
                wb = wb + 0.5*(pka(kk+1)+tmp)*(wzb(km+1)-wza(kk+1))
                if( kk <= k-2) then
                    do kt=kk+1, k-1
                       wb = wb + pza(kt)
                    enddo
                endif
                kb = kk
                goto 444
            endif
            enddo
444      continue
         endif

      elseif( wzb(k) < wza(k) ) then

!============================
! for wzb(k) < wza(k) 
!============================

! Compute wa

         if(k .ne. km+1) then
            do kk=k,km
               if( wzb(k) <= wza(kk)   .and.                   &
                   wzb(k) >= wza(kk+1)     ) then
                   tmp = pka(kk) + dka(kk)*(wza(kk)-wzb(k))    &
                                         / (wza(kk)-wza(kk+1))
                   wa = wa - 0.5*(tmp+pka(kk))*(wza(kk)-wzb(k))
                  goto 666
               else
                  wa  = wa - pza(kk)
               endif
            enddo
666      continue
         endif        ! k safety check

         if( wzb(k) < wza(km+1) ) then
! Integrate  down along the right edge!!
            do kk=ka,k-1
              if( wza(km+1) <= wzb(kk)   .and.                   &
                  wza(km+1) >= wzb(kk+1)   ) then
                  tmp = pkb(kk) + dkb(kk)*(wzb(kk)-wza(km+1))    &
                                        / (wzb(kk)-wzb(kk+1))
                  wa = wa - 0.5*(tmp+pkb(kk+1))*(wza(km+1)-wzb(kk+1))
                  if( kk <= k-2 ) then
                     do kt=kk+1, k-1
                        wa  = wa - pzb(kt)
                     enddo
                  endif
                  ka = kk
                 goto 777
              endif
            enddo
777      continue
         endif

! Top
       if( wza(k) > wzb(1) ) then 
          wb = 0.5*(pkb(k)+pka(k))*(wzb(k)-wza(k))
       else
          do kk=k,2,-1
             if( wza(k) <= wzb(kk-1) .and.                   &
                 wza(k) >= wzb(kk)      )  then
                tmp = pkb(kk-1)+dkb(kk-1)*(wzb(kk-1)-wza(k)) &
                                       / (wzb(kk-1)-wzb(kk))
                wb = wb + 0.5*(tmp+pkb(kk))*(wzb(kk)-wza(k))
                goto 999
             else
                wb  = wb - pzb(kk-1)
             endif
         enddo
999     continue
       endif   ! end top checck
      endif

      pf(k) = -(wa + wb) 

1000  continue

      return
      end
