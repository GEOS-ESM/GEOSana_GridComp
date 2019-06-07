!---------------------------------------------------------------------
      subroutine pl_vert_multi(nvectors,nplot,lm,ivar,efac,date,expid,hc
     *,growth,ifwd,ibak,res)
c
      use m_ioutil, only : luavail
                                                                                                                                               
      implicit none
                                                                                                                                               
!  INPUT VARIABLES
                                                                                                                                               
      integer, intent(in) :: nvectors,nplot,lm,ivar,ifwd,ibak
      character*16, intent(in) :: expid
      character*10, intent(in) ::  date
      logical, intent(in) :: hc
      real*4, intent(in) :: growth(nvectors)
      real*8, intent(in) :: efac
      character*3, intent(in) ::  res
!
!  LOCAL VARIABLES
!
      character*48 exfilvert
      character*16 elabel
      real xmin(3),xmax(3),ymin(3),ymax(3)
      integer ilead, iefac, imon
      integer npix, nf, npg, npl 
      integer j, j1, j2, j3
      integer lunit1, lunit2
c
      character*4 yr
      character*3 mon(12)
      character*2 z,iday
c
      data xmin/0.00,3.72,7.32/
      data xmax/3.80,7.32,10.98/
      data ymin/2.00,2.00,2.00/
      data ymax/6.00,6.00,6.00/
c
      data mon/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP'
     *, 'OCT','NOV','DEC'/
c
      ilead= ifwd-ibak
      iefac= int(efac)
c
      read(date,'(a4,i2,a2,a2)') yr,imon,iday,z
c
c************************************************************
c
c
      npix= nvectors
      print *,'npix= ',npix
c
      if(ivar.eq.1) elabel='  TOTAL ENERGY  '
      if(ivar.eq.2) elabel=' KINETIC ENERGY '
      if(ivar.eq.3) elabel='POTENTIAL ENERGY'
      if(ivar.eq.4) elabel='R-KINETIC ENERGY'
      if(ivar.eq.5) elabel='D-KINETIC ENERGY'
c
c----------------------------------------------------------------------
c 
!      growth(nvectors)= 0
!      do 10 j=1,nvectors-1
!      growth(nvectors)= growth(nvectors)+growth_in(j)
!   10 continue
!      growth(nvectors)= growth(nvectors)/(nvectors-1)
      print *,'Mean',nvectors-1,' growth rate= ',growth(nvectors)
c
c----------------------------------------------------------------------
c
      print*, 'description file for vert energy = evert.ctl'
c 
      lunit1 = luavail() 
      open(unit=lunit1,file='evert.ctl',form='formatted')
c
      write(lunit1,'(A)') 'DSET evert.dat'

      write(lunit1,'(A)') 'UNDEF  -9.99e33'
      write(lunit1,'(A)') 'OPTIONS big_endian'
      write(lunit1,'(A)') 'TITLE SV Energy Profile'
      write(lunit1,'(A)') '*'
      write(lunit1,'(A)') 'XDEF 1 LINEAR 1 1' 
      write(lunit1,'(A)') 'YDEF 1 LINEAR 1 1'
      if ( lm.eq.18) then
              write(lunit1,'(A,A)') 'ZDEF 18 LEVELS   8  29  57  93  138  190'
     *          ,'  251  319  398  488  583  678  769 850  916  962  991 1008'
      else if ( lm.eq.32) then
              write(lunit1,'(A,A,A,A)') 'ZDEF 32 LEVELS   0.7 1.5 2.8 5 8.1 12.6 18.4 25.3 32.9', 
     * ' 41.1 50.2 60.6 72.3 85.4 100.5 118.3 139.1 163.6 192.5',  
     * ' 226.2 265.7 312.1 366.4 430.3 505.3 593.4 687.6 776.4',  
     * ' 854.1 915.2 955.1 976.7 '
      else if ( lm.eq.55) then
              write(lunit1,'(A,A,A,A,A)') 'ZDEF 55 LEVELS 0.015 0.02 0.03 0.04 0.06 0.08 0.11 0.15 0.21',
     * ' 0.27 0.36 0.47 0.61 0.79 1.01 1.30 1.65 2.08 2.62 3.27 4.07 5.04 6.21 7.61 9.29 11.27',
     * ' 13.64 16.45 19.79 23.73 28.36 33.80 40.17 47.64 56.38 66.60 78.51 92.36 108.66 127.83',
     * ' 150.39 176.92 208.22 245.17 288.76 340.14 400.68 555.90 654.68 751.27 839.22 912.60',
     * ' 966.11 995.68 1010.73' 
      else if ( lm.eq.72) then
              write(lunit1,'(A,A,A,A,A,A,A)') 'ZDEF 72 LEVELS ',
     * ' 0.01  0.02  0.03  0.05  0.07  0.09  0.12  0.16  0.21  0.28  0.37  0.48',
     * ' 0.62  0.80  1.02  1.30  1.65  2.09  2.62  3.28  4.08  5.05  6.22  7.62',                  
     * '   9.29  11.28  13.64  16.46  19.79  23.73  28.37  33.81  40.18  47.64  56.39  66.60',
     * '  78.51  92.37 108.66 127.84 150.39 176.93 208.15 244.88 288.08 337.50 375.00 412.50',
     * ' 450.00 487.50 525.00 562.50 600.00 637.50 675.00 700.00 725.00 750.00 775.00 800.00',
     * ' 820.00 835.00 850.00 865.00 880.00 895.00 910.00 925.00 940.00 955.00 970.00 985.00'
      else 
              write(lunit1,'(A,i3,A)') 'ZDEF ',lm,' LINEAR 0.0 1.0'
      endif
      write(lunit1,444) z,iday,mon(imon),yr
  444 format('TDEF 2 LINEAR ',a2,'Z',a2,a3,a4,' 1dy')
      write(lunit1,'(A)') '*'
      write(lunit1,'(A,i3)') 'VARS ',npix                                                                          
c
      do 100 j=1,npix
      if(j.le.9) write(lunit1,'(A,i1,i5,A,i1)') 'ev',j,lm,'  99 SV',j                                                                
      if(j.gt.9) write(lunit1,'(A,i2,i5,A,i2)') 'ev',j,lm,'  99 SV',j
  100 continue
      write(lunit1,'(A)') 'ENDVARS'                                                                         
      close (lunit1)
c
c----------------------------------------------------------------------
c
      print*,'directives file for vert energy = evert.ex'
c
      exfilvert='evert.ex'
      print *,'exfilvert= ',exfilvert

      lunit2 = luavail()
      open(unit=lunit2,file=exfilvert,form='formatted')
c
      write(lunit2,'(A)') 'open evert.ctl'
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
c     if(npg.le.9) write(lunit2,'(A,i1)') 'enable print evert0',npg
c     if(npg.gt.9) write(lunit2,'(A,i2)') 'enable print evert',npg
c     write(lunit2,'(A,A,A,i2.2)') 'enable print evert',date,expid,npg
      write(lunit2,'(A)') 'enable print evert.gx'
      endif
c
      write(lunit2,'(A)') 'set ylopts 1 4 0.25'
      write(lunit2,'(A)') 'set xlopts 1 4 0.25'
      write(lunit2,'(A)') 'set yflip on'
      write(lunit2,'(A)') 'set zlog on'
      write(lunit2,'(A)') 'set lev 1 1000' 
      write(lunit2,'(A,i4)') 'set axlim 0 ', iefac
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
c     if(j.eq.npix) write(lunit2,'(A,a8,i5,A,f5.1)')
c    2 'draw title ',date,nvectors,' MEAN  IFAC=',efac
c
c     write(lunit2,'(A)') 'draw ylab LEVEL'
c     write(lunit2,'(A,a16)') 'draw xlab ',elabel
c
      if(npl.eq.nplot .or. j.eq.npix) then 
      write(lunit2,'(A)') 'set vpage off'
      write(lunit2,'(A)') 'set strsiz 0.13 0.15'
      write(lunit2,'(A)') 'set string 1 c 4 0'
      if(j1.le.npix)
     * write(lunit2,'(A,a16)') 'draw string 1.88 2.00 ',elabel
      if(j2.le.npix)
     * write(lunit2,'(A,a16)') 'draw string 5.52 2.00 ',elabel
      if(j3.le.npix)
     * write(lunit2,'(A,a16)') 'draw string 9.15 2.00 ',elabel
c
      write(lunit2,'(A)') 'set string 1 c 4 90'
      if(j1.le.npix)
     * write(lunit2,'(A)') 'draw string 0.07 4.00 PRESSURE'
      if(j2.le.npix)
     * write(lunit2,'(A)') 'draw string 3.68 4.00 PRESSURE'
      if(j3.le.npix)
     * write(lunit2,'(A)') 'draw string 7.31 4.00 PRESSURE'
c
      write(lunit2,'(A)') 'set string 1 c 4 0'
      if(j1.lt.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 1.90 5.85 SV#',j1
     *,'   svalue=',growth(j1)
      if(j2.lt.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 5.51 5.85 SV#',j2
     *,'   svalue=',growth(j2)
      if(j3.lt.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 9.14 5.85 SV#',j3
     *,'   svalue=',growth(j3)
c
      if(j1.eq.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 1.90 5.85 Mean',j1-1
     *,'   svalue=',growth(j1)
      if(j2.eq.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 5.51 5.85 Mean',j2-1
     *,'   svalue=',growth(j2)
      if(j3.eq.npix)
     * write(lunit2,'(A,i2,A,f4.1)')'draw string 9.14 5.85 Mean',j3-1
     *,'   svalue=',growth(j3)
c
      write(lunit2,'(A)') 'set strsiz 0.17 0.20'
      write(lunit2,'(A)') 'set string 2 c 4 0'
      write(lunit2,'(A,A)') 'draw string 5.50 7.50 '
     *,'SINGULAR VECTOR VERTICAL PROFILES'
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
     *,'GEOS Singular Vectors ', res, ' (+',ifwd,'h,-',ibak,'h)'
c     write(lunit2,'(A)') 'set string 2 c 4 0'
!      write(lunit2,'(A)') 'run tstamp_eng 2 4 2'
c
       if(hc) write(lunit2,'(A)') 'print'
       write(lunit2,'(A)') 'c'
c     if(hc) write(lunit2,'(A)') 'disable print'
      endif
c
  200 continue
c
!      write(lunit2,'(A)') 'quit'
      close (lunit2)
c

      return
      end
