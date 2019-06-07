      subroutine fft99a(a,work,trigs,inc,jump,n,lot)

      use precision

!      dimension a(jump*lot),work((n+1)*lot),trigs(3*n/2+1)
      real(r8) a(jump*lot),work((n+1)*lot),trigs(3*n/2+1)
c
c     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1
c     (SPECTRAL TO GRIDPOINT TRANSFORM)
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     A(0) AND A(N/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
cDIR$ IVDEP
      do 10 l=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   10 continue
c
c     REMAINING WAVENUMBERS
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1
c
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
cDIR$ IVDEP
      do 20 l=1,lot
      work(ja)=(a(ia)+a(ib))-
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(jb)=(a(ia)+a(ib))+
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+
     *    (a(ia+inc)-a(ib+inc))
      work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))-
     *    (a(ia+inc)-a(ib+inc))
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   20 continue
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     WAVENUMBER N/4 (IF IT EXISTS)
      ia=iabase
      ja=jabase
cDIR$ IVDEP
      do 40 l=1,lot
      work(ja)=2.0*a(ia)
      work(ja+1)=-2.0*a(ia+inc)
      ia=ia+jump
      ja=ja+nx
   40 continue
c
   50 continue
      return
      end
