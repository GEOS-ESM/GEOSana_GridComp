      subroutine initmem(im,jfirst,jlast, kfirst, klast, klastp,    &
                         u,v,pt,q,nq,pk,pe,delp,undef, ng_s, ng_d)

      use precision
      implicit none

      integer im, jfirst, jlast, kfirst, klast, klastp
      integer i, j, k, ic, nq, ng_s, ng_d

      real(r8) undef
      real(r8) :: u(im,jfirst-ng_d:jlast+ng_s,kfirst:klast)
      real(r8) :: v(im,jfirst-ng_s:jlast+ng_d,kfirst:klast)
      real(r8) :: pt(im,jfirst-ng_d:jlast+ng_d,kfirst:klast)            
      real(r8) :: q(im,jfirst-ng_d:jlast+ng_d,kfirst:klast,*)
      real(r8) :: pe(im,kfirst:klastp,jfirst:jlast)
      real(r8) :: delp(im,jfirst:jlast,kfirst:klast)
      real(r8) :: pk(im,jfirst:jlast,kfirst:klastp)

!$omp  parallel do         &
!$omp  default(shared)     &
!$omp  private(i,j,k)

      do 1000 k=kfirst,klast
        do j=jfirst-ng_d,jlast+ng_s
          do i=1,im
            u(i,j,k) = undef
          enddo
        enddo
        do j=jfirst-ng_s,jlast+ng_d
          do i=1,im
            v(i,j,k) = undef
          enddo
        enddo
        do j=jfirst-ng_d,jlast+ng_d
          do i=1,im
            pt(i,j,k) = undef
          enddo
        enddo
        do j=jfirst,jlast
          do i=1,im
            delp(i,j,k) = undef
          enddo
        enddo
1000  continue

      if(nq .ne. 0) then
        do 1500  ic=1,nq

!$omp  parallel do            &
!$omp  default(shared)        &
!$omp  private(i,j,k)

          do k=kfirst,klast
            do j=jfirst-ng_d,jlast+ng_d
              do i=1,im
                q(i,j,k,ic) = undef
              enddo
            enddo
          enddo

1500    continue

      endif

!$omp  parallel do            &
!$omp  default(shared)        &
!$omp  private(i,j,k)          

      do 2000  k=kfirst,klastp
        do j=jfirst,jlast
          do i=1,im
            pk(i,j,k) = undef
          enddo
        enddo
2000  continue

!$omp  parallel do           &
!$omp  default(shared)       &
!$omp  private(i,j,k)   

      do j=jfirst,jlast
        do k=kfirst,klastp
          do i=1,im
            pe(i,k,j) = undef
          enddo
        enddo
      enddo

      return
      end
