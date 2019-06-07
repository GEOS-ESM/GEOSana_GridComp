!-----------------------------------------------------------------------
!BOP
! !ROUTINE: pkez --- Calculate solution to hydrostatic equation
!
! !INTERFACE:
      subroutine pkez(nx, im, km, jfirst, jlast,      &
                      pe, pk, akap, ks, peln, pkz, eta)
!
! !USES:
      use precision

      implicit none

! !INPUT PARAMETERS:
      integer nx                          ! SMP decomposition in x
      integer im, km                  ! Dimensions
      integer jfirst, jlast               ! Latitude strip
      real(r8)  pe(im, km+1, jfirst:jlast)    ! Edge pressure
      integer ks
      logical eta     ! Is on ETA coordinate?
                      ! True:  input pe;     output pk, pkz, peln
                      ! False: input pe, pk; output pkz
      real(r8) akap

! !INPUT/OUTPUT PARAMETERS:
      real(r8)  pk(im,jfirst:jlast,km+1)
      real(r8) pkz(im,jfirst:jlast,km)

! !OUTPUT
      real(r8) peln(im, km+1, jfirst:jlast)   ! log pressure (pe) at layer edges

! !DESCRIPTION:
!
!
! !CALLED FROM:
!     te_map and fvgcm
!
! !REVISION HISTORY:
!
!     WS  99.05.19 : Removed fvcore.h
!     WS  99.07.27 : Limited region to jfirst:jlast
!     WS  99.10.22 : Deleted cp as argument (was not used)
!     WS  99.11.05 : Documentation; pruning of arguments
!     SJL 00.01.02: SMP decomposition in i
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local
      real(r8) pk2(im, km+1)
      real(r8) pek
      real(r8) lnp
      integer i, j, k
      integer ixj, jp, it, i1, i2

      it = im / nx
      jp = nx * ( jlast - jfirst + 1 )

!$omp  parallel do default(shared) private(ixj, i1, i2, i, j, k, pek, lnp, pk2)

!     do 1000 j=jfirst, jlast
!        i1 = 1
!        i2 = im

      do 1000 ixj=1,jp

         j  = jfirst + (ixj-1) / nx
         i1 =  1 + it * mod(ixj-1, nx)
         i2 = i1 + it - 1

        if ( eta ) then

! <<<<<<<<<<< Eta cordinate Coordinate  >>>>>>>>>>>>>>>>>>>
          pek =     pk(i1,j,1)
          lnp = log(pe(i1,1,j))

          do i=i1,i2
             pk2(i,1)   = pek
            peln(i,1,j) = lnp
          enddo

          if(ks .ne. 0) then
            do k=2, ks+1
              pek = pe(i1,k,j)**akap
              lnp = log(pe(i1,k,j))
              do i=i1,i2
                 pk2(i,k)   = pek
                peln(i,k,j) =  lnp
              enddo
            enddo

            do k=1, ks
              pek = (       pk2(i1,k+1)   - pk2(i1,k))   /  &
                    (akap*(peln(i1,k+1,j) - peln(i1,k,j)) )
              do i=i1,i2
                 pkz(i,j,k) = pek
              enddo
            enddo
          endif

          do k=ks+2,km
            do i=i1,i2
               pk2(i,k) = pe(i,k,j)**akap
            enddo
          enddo

          do i=i1,i2
             pk2(i,km+1) = pk(i,j,km+1)
          enddo

          do k=ks+2,km+1
            do i=i1,i2
               peln(i,k,j) =  log(pe(i,k,j))
            enddo
          enddo

          do k=ks+1,km
            do i=i1,i2
               pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k)) /           &
                            (akap*(peln(i,k+1,j) - peln(i,k,j)) )
            enddo
          enddo

          do k=2,km
            do i=i1,i2
               pk(i,j,k) = pk2(i,k)
            enddo
          enddo

        else

! <<<<<<<<<<< General Coordinate  >>>>>>>>>>>>>>>>>>>

          pek =     pk(i1,j,1)
          lnp = log(pe(i1,1,j))

          do i=i1,i2
             pk2(i,1) = pek
             peln(i,1,j) = lnp
          enddo

          do k=2,km+1
             do i=i1,i2
                peln(i,k,j) =  log(pe(i,k,j))
                pk2(i,k) =  pk(i,j,k)
             enddo
          enddo

          do k=1,km
             do i=i1,i2
                pkz(i,j,k) = (       pk2(i,k+1) - pk2(i,k) )  /  &
                             (akap*(peln(i,k+1,j) - peln(i,k,j)) )
             enddo
          enddo

        endif
1000  continue

      return
!EOC
      end
!-----------------------------------------------------------------------
