      subroutine fax(ifax,n,mode)

      use precision

!      dimension ifax(13)
      integer ifax(13)
      nn=n
      if (iabs(mode).eq.1) go to 10
      if (iabs(mode).eq.8) go to 10
      nn=n/2
      if ((nn+nn).eq.n) go to 10
      ifax(1)=-99
      return
   10 k=1
c     TEST FOR FACTORS OF 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
c     TEST FOR EXTRA FACTOR OF 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
c     TEST FOR FACTORS OF 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
c     NOW FIND REMAINING FACTORS
   50 l=5
      inc=2
c     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 l=l+inc
      inc=6-inc
      go to 60
   80 ifax(1)=k-1
c     IFAX(1) CONTAINS NUMBER OF FACTORS
      nfax=ifax(1)
c     SORT FACTORS INTO ASCENDING ORDER
      if (nfax.eq.1) go to 110
      do 100 ii=2,nfax
      istop=nfax+2-ii
      do 90 i=2,istop
      if (ifax(i+1).ge.ifax(i)) go to 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 continue
  100 continue
  110 continue
      return
      end
