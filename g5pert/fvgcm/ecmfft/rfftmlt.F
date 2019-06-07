      subroutine rfftmlt(a,work,trigs,ifax,inc,jump,n,lot,isign)

      use precision

#if defined( SGI_FFT )
!
! WS 99.11.24 : Added SGI wrappers
!
      implicit none
      integer inc, jump, n, lot, isign
      real(r8)   a(jump,lot)
      real(r8)   trigs(1)
      real(r8)   work(1)                       ! Not used; here for plug reason
      integer ifax(*)
! Local
      integer*4 iisign,in,iinc,ijump,ilot
      integer i, j
      real(r8) scale

c-----convert to i4
      iisign = isign
      iinc = inc
      ijump = jump
      in = n
      ilot = lot

      if(iisign.lt.0) then
c-----forward
          call dzfftm1du (iisign,in,ilot,a,iinc,ijump,trigs)
       endif

      if(iisign.gt.0) then
c-----backward
          call zdfftm1du (iisign,in,ilot,a,iinc,ijump,trigs)

          scale = 1.0/float(n)
          do j=1,lot
             do i=1,jump
                a(i,j) = scale*a(i,j)
             enddo
          enddo
       endif
#else
!
! Default version
!
      real(r8) a(jump*lot+inc)       ! FastOpt
      real(r8) work((n+2)*lot+inc)   ! FastOpt
      real(r8) trigs(3*n/2+1)
      integer ifax(13)
c
c     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
c     FAST FOURIER TRANSFORM
c
c     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
c     THAT IN MRFFT2
c
c     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
c     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
c     (1970), 315-337)
c
c     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
c     WORK IS AN AREA OF SIZE (N+1)*LOT
c     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
c     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
c     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
c         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
c     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
c     N IS THE LENGTH OF THE DATA VECTORS
c     LOT IS THE NUMBER OF DATA VECTORS
c     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
c           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
c
c     ORDERING OF COEFFICIENTS:
c         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
c         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
c
c     ORDERING OF DATA:
c         X(0),X(1),X(2),...,X(N-1)
c
c     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
c     PARALLEL
c
c     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
c
c     DEFINITION OF TRANSFORMS:
c     -------------------------
c
c     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
c         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
c
c     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
c               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
c
c
c
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
c
c     IF NECESSARY, TRANSFER DATA TO WORK AREA
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
cDIR$ IVDEP
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
c
      igo=60
      go to 40
c
c     PREPROCESSING (ISIGN=+1)
c     ------------------------
c
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     COMPLEX TRANSFORM
c     -----------------
c
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c
      if (isign.eq.-1) go to 130
c
c     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
cDIR$ IVDEP
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     FILL IN ZEROS AT END
  110 continue
      ib=n*inc+1
cDIR$ IVDEP
      do 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 continue
      go to 140
c
c     POSTPROCESSING (ISIGN=-1):
c     --------------------------
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
#endif
      return
      end
