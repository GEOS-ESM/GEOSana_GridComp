
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine ppme(p,qe,delp,im,km)
!****6***0*********0*********0*********0*********0*********0**********72
      use precision
      implicit none
      integer(i4)  im, km, km1
      integer(i4)  i, k
      real(r8) p(im,km),qe(im,km+1),delp(im,km)
! local arrays.
      real(r8) dc(im,km),delq(im,km), a6(im,km)
      real(r8) c1, c2, c3, tmp, qmax, qmin
      real(r8) a1, a2, s1, s2, s3, s4, ss3, s32, s34, s42
      real(r8) a3, b2, sc, dm, d1, d2, f1, f2, f3, f4
      real(r8) qm, dq

      km1 = km - 1

      do 500 k=2,km
      do 500 i=1,im
500   a6(i,k) = delp(i,k-1) + delp(i,k)

      do 1000 k=1,km1
      do 1000 i=1,im
      delq(i,k) = p(i,k+1) - p(i,k)
1000  continue

      do 1220 k=2,km1
      do 1220 i=1,im
      c1 = (delp(i,k-1)+0.5*delp(i,k))/a6(i,k+1)
      c2 = (delp(i,k+1)+0.5*delp(i,k))/a6(i,k)
      tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /    &
                                    (a6(i,k)+delp(i,k+1))
      qmax = max(p(i,k-1),p(i,k),p(i,k+1)) - p(i,k)
      qmin = p(i,k) - min(p(i,k-1),p(i,k),p(i,k+1))
      dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
1220  continue

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

      do 12 k=3,km1
      do 12 i=1,im
      c1 = delq(i,k-1)*delp(i,k-1) / a6(i,k)
      a1 = a6(i,k-1) / (a6(i,k) + delp(i,k-1))
      a2 = a6(i,k+1) / (a6(i,k) + delp(i,k))
      qe(i,k) = p(i,k-1) + c1 + 2./(a6(i,k-1)+a6(i,k+1)) *        &
                ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -         &
                                delp(i,k-1)*a1*dc(i,k  ) )
12    continue

! three-cell parabolic subgrid distribution at model top

      do 10 i=1,im
! three-cell PP-distribution
! Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
! a3 = a / 3
! b2 = b / 2
      s1 = delp(i,1)
      s2 = delp(i,2) + s1
!
      s3 = delp(i,2) + delp(i,3)
      s4 = s3 + delp(i,4)
      ss3 =  s3 + s1
      s32 = s3*s3
      s42 = s4*s4
      s34 = s3*s4
! model top
      a3 = (delq(i,2) - delq(i,1)*s3/s2) / (s3*ss3)
!
      if(abs(a3) .gt. 1.e-14) then
         b2 =  delq(i,1)/s2 - a3*(s1+s2)
         sc = -b2/(3.*a3)
         if(sc .lt. 0. .or. sc .gt. s1) then
             qe(i,1) = p(i,1) - s1*(a3*s1 + b2)
         else
             qe(i,1) = p(i,1) - delq(i,1)*s1/s2
         endif
      else
! Linear
         qe(i,1) = p(i,1) - delq(i,1)*s1/s2
      endif
      dc(i,1) = p(i,1) - qe(i,1)
! compute coef. for the off-centered area preserving cubic poly.
      dm = delp(i,1) / (s34*ss3*(delp(i,2)+s3)*(s4+delp(i,1)))
      f1 = delp(i,2)*s34 / ( s2*ss3*(s4+delp(i,1)) )
      f2 = (delp(i,2)+s3) * (ss3*(delp(i,2)*s3+s34+delp(i,2)*s4)   &
            + s42*(delp(i,2)+s3+s32/s2))
      f3 = -delp(i,2)*( ss3*(s32*(s3+s4)/(s4-delp(i,2))            &
            + (delp(i,2)*s3+s34+delp(i,2)*s4))                     &
            + s42*(delp(i,2)+s3) )
      f4 = ss3*delp(i,2)*s32*(delp(i,2)+s3) / (s4-delp(i,2))
      qe(i,2) = f1*p(i,1)+(f2*p(i,2)+f3*p(i,3)+f4*p(i,4))*dm
10    continue

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do 15 i=1,im
      d1 = delp(i,km)
      d2 = delp(i,km1)
      qm = (d2*p(i,km)+d1*p(i,km1)) / (d1+d2)
      dq = 2.*(p(i,km1)-p(i,km)) / (d1+d2)
      c1 = (qe(i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
      qe(i,km  ) = qm - c1*d1*d2*(d2+3.*d1)
      qe(i,km+1) = d1*(8.*c1*d1**2-c3) + qe(i,km)
15    continue
      return
      end
