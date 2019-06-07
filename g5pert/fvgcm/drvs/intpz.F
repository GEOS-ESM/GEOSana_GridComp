      subroutine intpz(q, qm, jm, km, p, prk, kw_top, kw_bot, nl)
      implicit none

C Input
      integer jm
      real p(*)       ! coordinate at which the data q is defined
      real prk(*)      ! output coordinate

      real q(jm,*)    ! original data
      integer km
      integer kw_top, kw_bot
      integer nl

C Output
      real qm(nl, jm)       ! output interpolated data
 
C local variables
      integer l, j, k, k0
      real s,  q1d

C Pure p-ccordinate

      k0 = 1
      do 500 k=1, kw_bot

      if(prk(k) .lt. p(1)) then
C Top
         kw_top = max(1, k-1)

!$omp  parallel do
!$omp& default(shared)
!$omp& private(j)

      do j=1,jm
         qm(k,j) = q(j,1)
      enddo

      elseif(prk(k) .ge. p(km)) then
! Bottom

!$omp  parallel do
!$omp& default(shared)
!$omp& private(j)
      do j=1,jm
         qm(k,j) = q(j,km)
      enddo

      else
        do l=k0,km-1
        if(prk(k) .ge. p(l) .and. prk(k) .le. p(l+1)) then
              s  = (prk(k)-p(l)) / (p(l+1) -p(l))

!$omp  parallel do
!$omp& default(shared)
!$omp& private(j)

         do j=1,jm
            qm(k,j) = q(j,l) + s*(q(j,l+1) - q(j,l))
         enddo

         k0 = l
         goto 500
        endif
        enddo
      endif
500   continue

      return
      end
