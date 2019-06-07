!---------------------------------------------------------------------
      subroutine pl_spec_multi(nvectors,nplot,jtrun,ivar,efac,date,expid
     2                  ,itype,hc,growth,ifwd,ibak,res)

      use m_ioutil, only : luavail

      implicit none

!  INPUT VARIABLES

      integer, intent(in) :: nvectors,nplot,jtrun,ivar,itype,ifwd,ibak
      character*16, intent(in) :: expid
      character*10, intent(in) ::  date
      logical, intent(in) :: hc
      real*4, intent(in) :: growth(nvectors)
      real*8, intent(in) :: efac
      character*3, intent(in) ::  res
!
!  LOCAL VARIABLES
!
      character*48 exfilspec
      character*16 elabel,wlabel
      real xmin(3),xmax(3),ymin(3),ymax(3)
      integer ilead, iefac, imon, npix, jtrunm1
      integer i, j, k
      integer nf, npl, npg
      integer j1, j2, j3
      integer lunit1, lunit2
!
      character*4 yr
      character*3 mon(12)
      character*2 z,iday
!
      data xmin/0.00,3.50,7.00/
      data xmax/3.95,7.45,10.95/
      data ymin/2.00,2.00,2.00/
      data ymax/6.00,6.00,6.00/
!
      data mon/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP'
     *, 'OCT','NOV','DEC'/
!
      ilead= ifwd-ibak
      iefac= int(efac)
!
      read(date,'(a4,i2,a2,a2)') yr,imon,iday,z
!
c************************************************************
c
c
      npix= nvectors
      jtrunm1= jtrun-1
c
      if(ivar.eq.1) elabel='  TOTAL ENERGY  '
      if(ivar.eq.2) elabel=' KINETIC ENERGY '
      if(ivar.eq.3) elabel='POTENTIAL ENERGY'
      if(ivar.eq.4) elabel='R-KINETIC ENERGY'
      if(ivar.eq.5) elabel='D-KINETIC ENERGY'
c
      if(itype.eq.1) wlabel='TOTAL WAVENUMBER'
      if(itype.eq.2) wlabel='ZONAL WAVENUMBER'
c
c----------------------------------------------------------------------
c
!      growth(nvectors)= 0.
!      do 10 j=1,nvectors-1
!      growth(nvectors)= growth(nvectors)+growth_in(j)
!   10 continue
!      growth(nvectors)= growth(nvectors)/(nvectors-1)
c
c----------------------------------------------------------------------
c
      print*, 'description file for spectral energy = espec.ctl'
c 
      lunit1 = luavail()
      open(unit=lunit1,file='espec.ctl',form='formatted')
c
      write(lunit1,'(A)') 'DSET espec.dat'                                                           
      write(lunit1,'(A)') 'UNDEF  -9.99e33'
      write(lunit1,'(A)') 'OPTIONS big_endian'
      write(lunit1,'(A)') 'TITLE SV Energy Spectra'
      write(lunit1,'(A)') '*'
      write(lunit1,'(A,i3,A)') 'XDEF ',jtrun,' LINEAR 0.0  1.0'
c     write(lunit1,'(A)') 'XDEF 1 LINEAR 1 1'
      write(lunit1,'(A)') 'YDEF 1 LINEAR 1 1'
c     write(lunit1,'(A,i3,A)') 'YDEF ',jtrun+2,' LINEAR 0.0  1.0'
      write(lunit1,'(A)') 'ZDEF 1 LINEAR 1 1' 
      write(lunit1,444) z,iday,mon(imon),yr
  444 format('TDEF 2 LINEAR ',a2,'Z',a2,a3,a4,' 1dy')
      write(lunit1,'(A)') '*'
      write(lunit1,'(A,i3)') 'VARS ',npix                                                                          
c
      do 100 j=1,npix
      if(j.le.9)write(lunit1,'(A,i1,A,i1)') 'ev', j, '  0  99 SV',j                                                                
      if(j.gt.9)write(lunit1,'(A,i2,A,i2)') 'ev', j, '  0  99 SV',j
  100 continue
      write(lunit1,'(A)') 'ENDVARS'                                                                         
      close (lunit1)
c
c----------------------------------------------------------------------
c
      print*,'directives file for spectral energy = espec.ex'
c
      exfilspec='espec.ex'
      lunit2 = luavail()
      open(unit=lunit2,file=exfilspec,form='formatted')
c
      write(lunit2,'(A)') 'open espec.ctl'
      write(lunit2,'(A)') 'set grid off'
c
      do 200 j=1,npix
c
      nf= (j-1)/nplot
      npg= nf+1
      npl= j-nf*nplot
      print *,'nf= ',nf,'      npg= ',npg,'     npl= ',npl
c
      j1= nf*nplot+1
      j2= nf*nplot+2
      j3= nf*nplot+3
c
      if(nplot.gt.1) write(lunit2,'(A,4f6.2)') 'set vpage ',
     2  xmin(npl),xmax(npl),ymin(npl),ymax(npl)
c
      write(lunit2,'(A)') 'set grads off'
      if(hc.and.j.eq.1) then
c     if(npg.le.9) write(lunit2,'(A,i1)') 'enable print espec0',npg
c     if(npg.gt.9) write(lunit2,'(A,i2)') 'enable print espec',npg
c     write(lunit2,'(A,A,A,i2.2)') 'enable print espec',date,expid,npg
      write(lunit2,'(A)') 'enable print espec.gx'
      endif
c
      write(lunit2,'(A)') 'set ylopts 1 4 0.25'
      write(lunit2,'(A)') 'set xlopts 1 4 0.25'
      write(lunit2,'(A,i4)') 'set axlim 0 ', iefac 
      write(lunit2,'(A)') 'set gxout line'
c      write(lunit2,'(A,i4,A)') 'set xaxis 0 ',jtrunm1,' 5'
      if ((res(1:1) == 'a').or.(res(1:1) == 'A')) then
          write(lunit2,'(A)') 'set xaxis 0 36 5'
      elseif ((res(1:1) == 'b').or.(res(1:1) == 'B')) then
          write(lunit2,'(A)') 'set xaxis 0 72 10'
      elseif ((res(1:1) == 'c').or.(res(1:1) == 'C')) then
          write(lunit2,'(A)') 'set xaxis 0 144 20'
      else
          write(lunit2,'(A,i4,A)') 'set xaxis 0 ',jtrunm1,' 20'
      endif
      write(lunit2,'(A)') 'set cstyle 2'
      write(lunit2,'(A)') 'set ccolor 9'
      write(lunit2,'(A)') 'set cthick 8'
      write(lunit2,'(A)') 'set cmark 0'
      write(lunit2,'(A)') 'set t 1'
c
      if(j.le.9) write(lunit2,'(A,i1,A,f5.1)') 'd ev',j,'*',efac
      if(j.gt.9) write(lunit2,'(A,i2,A,f5.1)') 'd ev',j,'*',efac
c
      write(lunit2,'(A)') 'set cstyle 1'
      write(lunit2,'(A)') 'set ccolor 4'
      write(lunit2,'(A)') 'set cthick 8'
      write(lunit2,'(A)') 'set cmark 0'
      write(lunit2,'(A)') 'set t 2'
c
      if(j.le.9) write(lunit2,'(A,i1)') 'd ev',j
      if(j.gt.9) write(lunit2,'(A,i2)') 'd ev',j
c
c     if(j.lt.npix) write(lunit2,'(A,a8,A,i2,A,f5.1)')
c    2 'draw title ',date,'  SV=',j,'  IFAC=',efac
c     if(j.lt.npix) write(lunit2,'(A,i1,A,f4.1)')
c    2 'draw title SV#',j,'   amp=',growth(j)
c     if(j.eq.npix) write(lunit2,'(A,a8,i5,A,f5.1)')
c    2 'draw title ',date,nvectors,' MEAN  IFAC=',efac
c
c     write(lunit2,'(A,a16)') 'draw xlab ',wlabel
c     write(lunit2,'(A,a16)') 'draw ylab ',elabel
c
      if(npl.eq.nplot .or. j.eq.npix) then 
      write(lunit2,'(A)') 'set vpage off'
      write(lunit2,'(A)') 'set strsiz 0.13 0.15'
      write(lunit2,'(A)') 'set string 1 c 4 0'
      if(j1.le.npix)
     * write(lunit2,'(A)') 'draw string 2.30 2.00 TOTAL WAVENUMBER'
      if(j2.le.npix)
     * write(lunit2,'(A)') 'draw string 5.80 2.00 TOTAL WAVENUMBER'
      if(j3.le.npix)
     * write(lunit2,'(A)') 'draw string 9.30 2.00 TOTAL WAVENUMBER'
c
      write(lunit2,'(A)') 'set string 1 c 4 90'
      if(j1.le.npix)
     * write(lunit2,'(A,a16)') 'draw string 0.55 4.00 ',elabel
      if(j2.le.npix)
     * write(lunit2,'(A,a16)') 'draw string 4.05 4.00 ',elabel
      if(j3.le.npix)
     * write(lunit2,'(A,a16)') 'draw string 7.55 4.00 ',elabel
c
      write(lunit2,'(A)') 'set string 1 c 4 0'
      if(j1.lt.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 2.30 5.85 SV#',j1
     *,'   svalue=',growth(j1)
      if(j2.lt.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 5.80 5.85 SV#',j2
     *,'   svalue=',growth(j2)
      if(j3.lt.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 9.30 5.85 SV#',j3
     *,'   svalue=',growth(j3)
c
      if(j1.eq.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 2.30 5.85 Mean',j1-1
     *,'   svalue=',growth(j1)
      if(j2.eq.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 5.80 5.85 Mean',j2-1
     *,'   svalue=',growth(j2)
      if(j3.eq.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 9.30 5.85 Mean',j3-1
     *,'   svalue=',growth(j3)
c
      write(lunit2,'(A)') 'set strsiz 0.17 0.20'
      write(lunit2,'(A)') 'set string 2 c 4 0'
      write(lunit2,'(A,A)') 'draw string 5.50 7.50 '
     *,'SINGULAR VECTOR WAVENUMBER SPECTRA'
c
      write(lunit2,'(A)') 'set strsiz 0.13 0.15'
      write(lunit2,'(A)') 'set string 9 l 4 0'
      write(lunit2,'(A,A,i2,A)')'draw string 3.50 6.80 '
     *,'Energy(x',iefac,') at Initial Time  (DASH)'
      write(lunit2,'(A)') 'set string 4 l 4 0'
      write(lunit2,'(A,A)')'draw string 3.50 6.50 '
     *,'Energy at Verification Time  (SOLID)'
c
      write(lunit2,'(A)') 'set strsiz 0.15 0.17'
      write(lunit2,'(A)') 'set string 1 l 4 0'
      write(lunit2,'(A,A,A,A,i2,A,i2,A)') 'draw string 1.20 1.00 '
     *,'GEOS Singular Vectors ',res,' (+',ifwd,'h,-',ibak,'h)'
c     write(lunit2,'(A)') 'set string 2 c 4 0'
!      write(lunit2,'(A)') 'run tstamp_eng 2 4 2'
c
       if(hc) write(lunit2,'(A)') 'print'
       write(lunit2,'(A)') 'c'
c     if(hc) write(lunit2,'(A)') 'disable print'
      endif

  200 continue
c
!      write(lunit2,'(A)') 'quit'
      close (lunit2)
c
      return
      end
