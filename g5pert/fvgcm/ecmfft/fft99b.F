      subroutine fft99b(work,a,trigs,inc,jump,n,lot)

      use precision

!      dimension work((n+1)*lot),a(jump*lot),trigs(3*n/2+1)
      real(r8) work((n+1)*lot),a(jump*lot),trigs(3*n/2+1)
c
c     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
c     (GRIDPOINT TO SPECTRAL TRANSFORM)
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     A(0) AND A(N/2)
      scale=1.0/float(n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
cDIR$ IVDEP
      do 10 l=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc)=0.0
      a(jb+inc)=0.0
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   10 continue
c
c     REMAINING WAVENUMBERS
      scale=0.5*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1
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
      a(ja)=scale*((work(ia)+work(ib))
     *   +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(jb)=scale*((work(ia)+work(ib))
     *   -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    +(work(ib+1)-work(ia+1)))
      a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    -(work(ib+1)-work(ia+1)))
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   20 continue
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     WAVENUMBER N/4 (IF IT EXISTS)
      ia=iabase
      ja=jabase
      scale=2.0*scale
cDIR$ IVDEP
      do 40 l=1,lot
      a(ja)=scale*work(ia)
      a(ja+inc)=-scale*work(ia+1)
      ia=ia+nx
      ja=ja+jump
   40 continue
c
   50 continue
      return
      end
