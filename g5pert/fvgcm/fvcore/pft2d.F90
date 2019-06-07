!-----------------------------------------------------------------------
! FastOpt : flow directives
!$taf subroutine rfftmlt   input = 1  ,3,4,5,6,7,8,9
!$taf subroutine rfftmlt  output = 1
!$taf subroutine rfftmlt  active = 1
!$taf subroutine rfftmlt  depend =     3,4,5,6,7,8,9
!$taf subroutine rfftmlt  adkill =   2          
!$taf subroutine rfftmlt ftlname = rfftmlt
!$taf subroutine rfftmlt  adname = adrfftmlt
!-----------------------------------------------------------------------

 subroutine pft2d(p, s, damp, im,  jp, ifax, trigs, q1, q2)

 use precision
 implicit none

! Input
     integer im, jp
     integer ifax(13)
     real(r8)  s(jp)              ! 3-point algebraic filter
     real(r8)  damp(im,jp)        ! FFT damping coef
     real(r8) trigs(3*im/2+1)

! Input/Output
     real(r8) p(im,jp)            ! Array to be polar filtered

! Local
     integer i, j
     real(r8) rsc, bt 
     integer nj
     integer n

     real(r8) ptmp(0:im+1)
!    real(r8) q1( im+2, jp)
!    real(r8) q2((im+1)*jp)
     real(r8) q1( im+2, *)
     real(r8) q2(*)
     integer  jf(jp)

      nj = 0

      do 200 j=1,jp

#if defined ( ALT_PFT )
      if(s(j) > 1.01) then
#else
      if(s(j) > 1.01 .and. s(j) <= 4.) then

         rsc = 1./s(j)
         bt  = 0.5*(s(j)-1.)

         do i=1,im
            ptmp(i) = p(i,j)
         enddo
           ptmp(im+1) = p(1 ,j)
           ptmp(   0) = p(im,j)

         do i=1,im
            p(i,j) = rsc * ( ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)) )
         enddo

      elseif(s(j) > 4.) then
#endif

! Packing for FFT 
           nj  = nj + 1
         jf(nj) = j

         do i=1,im
            q1(i,nj) = p(i,j)
         enddo
            q1(im+1,nj) = 0.
            q1(im+2,nj) = 0.

      endif
200   continue

      if( nj == 0) return
      call rfftmlt(q1, q2, trigs, ifax, 1, im+2, im, nj, -1)

      do n=1,nj
         do i=5,im+2
            q1(i,n) = q1(i,n) * damp(i-2,jf(n))
         enddo
      enddo

      call rfftmlt(q1, q2, trigs, ifax, 1, im+2, im, nj, 1)

      do n=1,nj
         do i=1,im
            p(i,jf(n)) = q1(i,n)
         enddo
      enddo
 return
 end
