module m_delp2ps
      use precision
      private
      PUBLIC delp2ps
      interface delp2ps; module procedure &
             delp2psA_,&
             delp2psB_
      end interface
CONTAINS
      subroutine delp2psA_ ( delp, ps )
      use stepon, only : imr
      use stepon, only : jfirst, jlast
      use stepon, only : nl
      implicit none
      real(r8), intent(in)  :: delp(imr,jfirst:jlast,nl)
      real(r8), intent(out) :: ps(imr,jfirst:jlast)
      integer i,j,k
      ps=0.d0
      do k=1,nl
         do j=jfirst,jlast
            do i=1,imr
               ps(i,j) = ps(i,j) + delp(i,j,k)
            end do
         end do
      end do
      end subroutine delp2psA_

      subroutine delp2psB_ ( im, jfirst, jlast, nl, delp, ps )
      implicit none
      integer,  intent(in) :: im, jfirst, jlast, nl
      real(r8), intent(in) :: delp(im,jfirst:jlast,nl)
      real(r8), intent(out) :: ps(im,jfirst:jlast)
      integer i,j,k
      ps=0.d0
      do k=1,nl
         do j=jfirst,jlast
            do i=1,im
               ps(i,j) = ps(i,j) + delp(i,j,k)
            end do
         end do
      end do
      end subroutine delp2psB_
end module m_delp2ps
