      program mapz

! ****************
! Vertical mapping
! ****************

! this program is for mapping from an arbitrary vertical domain
! with to an arbitrary vertical domain with the same surface pressure
! and same horizontal resolution

! Developer: S.-J. Lin
! Mar 27, 2000

      integer im, jm, nl
      integer ih

! Horizontal resolution: (im,jm)
! Number of constituents including water vapor: nc
      parameter ( nc = 1)

      real ptop

! Original data with km layers

! vertical resolution of the target: nl

      write(6,*) '*************************************************'
      write(6,*) 'This program is for the vertical mapping of fvgcm'
      write(6,*) 'dynamical restart file: d_rst'
      write(6,*) '*************************************************'

      write(6,*) ' '
      write(6,*) 'Resoultion? Choose from the following:'
      write(6,*) '0: 4x5;  1: 2x2.5;  2: 1x1.25;  3: 0.5x0.625'
      read(*,*) ih

      if( ih .eq. 0) then
          im = 72
          jm = 46
      elseif( ih .eq. 1) then
          im = 144
          jm = 91
      elseif( ih .eq. 2) then
          im = 288
          jm = 181
      elseif( ih .eq. 3) then
          im = 576
          jm = 361
      else
          write(6,*) 'No suitable resolution chosen'
          stop
      endif
      
      write(6,*) 'Original vertical dimension km=?'
      read(*,*) km

      write(6,*) 'Original model top (Pa)'
      read(*,*) ptop

      write(6,*) 'Vertical dimension of the target nl=?'
      read(*,*) nl

      write(6,*) ' '
      write(6,*) 'Input file name is assumed to be d_rst'
      write(6,*) 'Output file will be d_rst_new'
      iuic = 71
      iout = 81
      open (unit=iuic,file='d_rst',form='unformatted',status='old')
      open (unit=iout,file='d_rst_new',form='unformatted',
     &      status='unknown')

      call z_mapping(iuic, iout, im, jm, km, nc, ptop, nl) 

      end

      subroutine z_mapping(iuic, iout, im, jm, km, nc, ptop_old, nl) 

      integer im, jm, nl

! Original data with km layers
! vertical resolution of the target is nl

! Original
      real     um(im,jm,km)        !zonal wind on D-grid
      real     vm(im,jm,km)        !meridional wind
      real    ptm(im,jm,km)        !scaled potential temperature
      real     qm(im,jm,km,nc)     !specific humidity & mixing ratios
      real    pem(im,km+1,jm)      !pressure at layer edges
      real  delpm(im,jm,km)        ! pressure thickness
      real    pkm(im,jm,km+1)
      real ptop_old

! New data with nl layers
      real    ps(im, jm)     
      real     u(im,jm,nl)        !zonal wind on D-grid
      real     v(im,jm,nl)        !meridional wind
      real    pt(im,jm,nl)        !scaled potential temperature
      real     q(im,jm,nl,nc)     !specific humidity & mixing ratios
      real    pe(im,nl+1,jm)      !pressure at layer edges
      real  delp(im,jm,nl)        ! pressure thickness
      real    pk(im,jm,nl+1)
      real ak(nl+1), bk(nl+1)


! Constants:
      rair   = 287.04
      rh2o   = 4.61e2
      cpair  = 1004.64
      akap  = rair/cpair
      zvir   = rh2o/rair - 1.

      write(6,*) 'zvir=',zvir

c Set up the NEW reference vertical coordinate system

      if ( nl .eq. 18 ) then
           call set_18(nl, ks, ptop, pint, ak, bk)
      elseif( nl .eq. 32 ) then
           call set_32(nl, ks, ptop, pint, ak, bk)
      elseif( nl .eq. 55 ) then
           call set_55(nl, ks, ptop, pint, ak, bk)
      else
!        write(6,*) 'No specific set_eta specified'
         call set_eta(nl, ks, ptop, pint, ak, bk)
      endif

      write(6,*) 'PTOP=', ptop, pint, ks

C get original data with km layers
	call get_org_data(ps, delpm,um,vm,ptm,qm,pem,im,jm,km,nc,
     &                    nymd, nhms, iuic, nstep, ptop_old)

      write(6,*) ' '
      write(6,*) 'NYMD=',nymd,' NHMS=',nhms
      write(6,*) 'nstep=', nstep

      if( km .eq. nl  .and. abs(ptop-ptop_old) .le. 1.e-2) then

      write(6,*) 'No mapping will be performed'
      write(6,*) 'New date/time information?'
      write(6,*) 'NYMD=?'
      read(*,*) nymd
      write(6,*) 'NHMS=?'
      read(*,*) nhms
      write(6,*) 'nstep='
      read(*,*)  nstep

! write original dyn rest file.
      call rst_dyn(im, jm, km, -1, iout, delpm, um, vm, ptm, nc,
     &             qm   , ps ,  nymd,  nhms, nstep )
      else
 
	do k=1,km+1
	do j=1,jm
	do i=1,im
           pkm(i,j,k) = pem(i,k,j) ** akap
	enddo
	enddo
	enddo

      pmax = vmax(pem,pmin,im*jm*(km+1))
      write(6,*) 'max PEm =', pmax, ' min PEm=', pmin

      pmax = vmax(delpm,pmin,im*jm*km)
      write(6,*) 'max delp =', pmax, ' min delp=', pmin

C Mapping.
      call gmap(im, jm, nc,    akap,
     &          km,  pkm,  pem,         um,  vm,  ptm,  qm,
     &          nl,  pk ,  pe , delp ,  u ,  v ,  pt ,  q)

      write(6,*) ' '
      write(6,*) 'Done mapping'
      write(6,*) 'PTOP=', pe(1,1,1)

C Fix Dry air mass
!     call drymadj( im, jm, nl, 1, jm,
!    &              .true., ak, bk, ps, delp, nc, q, 1 )

C write new dyn rest file.
      call rst_dyn(im, jm, nl, -1 , iout , delp, u, v, pt, nc,
     &             q   , ps ,  nymd,  nhms, nstep )

      pmax = vmax(ps,pmin,im*jm)
      write(6,*) 'max ps =', pmax, ' min =', pmin

      pmax = vmax(delp,pmin,im*jm*nl)
      write(6,*) 'max delp =', pmax, ' min delp=', pmin

      pmax = vmax(pe,pmin,im*jm*(nl+1))
      write(6,*) 'max PE =', pmax, ' min PE=', pmin

      pmax = vmax(u,pmin,im*jm*nl)
      write(6,*) 'max u =', pmax, ' min u=', pmin

      pmax = vmax(v,pmin,im*jm*nl)
      write(6,*) 'max v =', pmax, ' min v=', pmin

      pmax = vmax(PT,pmin,im*jm*nl)
      write(6,*) 'max PT =', pmax, ' min PT=', pmin

      pmax = vmax(q,pmin,im*jm*nl)
      write(6,*) 'max q =', pmax, ' min q=', pmin

      endif

      end

      subroutine get_org_data(ps,delp,u,v,pt,q,pe,im,jm,km,nq,
     &                        nymd, nhms, iuic, nstep, ptop_old)

      integer  im, jm
      integer  nq

      real    ps(im,jm)     
      real     u(im,jm,km)        !zonal wind on D-grid
      real     v(im,jm,km)        !meridional wind
      real    pt(im,jm,km)        !scaled potential temperature
      real     q(im,jm,km,nq)     !specific humidity & mixing ratios
      real    pe(im,km+1,jm)      !pressure at layer edges
      real  delp(im,jm,km)        ! pressure thickness
      real ptop_old

      call rst_dyn(im, jm, km, 1 , iuic , delp, u, v, pt, nq,
     &             q   , ps ,  nymd,  nhms, nstep )

      write(6,*) 'NYMD=',nymd,' NHMS=',nhms, nstep

	do j=1,jm
           do i=1,im
              pe(i,1,j) = ptop_old
           enddo
	enddo

	do 1000 k=1,km
	do j=1,jm
	do i=1,im
           pe(i,k+1,j) = pe(i,k,j) + delp(i,j,k)
	enddo
	enddo
1000    continue

	return
	end

c****6***0*********0*********0*********0*********0*********0**********72
      subroutine rst_dyn(im, jm, km, itp , iuic , delp, u, v, pt, nq,
     &                   q   , ps ,  nymd,  nhms, nstep )
c****6***0*********0*********0*********0*********0*********0**********72

      implicit none
      integer  im, jm
      integer  km

      integer nymd, nhms, nstep
      integer nq, iuic
      integer itp
      integer i, j, k

      real    u(im,jm,km),    v(im,jm,km),  pt(im,jm,km),
     &        delp(im,jm,km),
     &        q(im,jm,km,nq), ps(im,jm)

      real vmax, pmax, pmin

      if(itp.ge.0) then

      read(iuic) nstep, nymd, nhms
      read(iuic) ps, delp, u, v, pt

      if(itp .eq. 0) then
         nstep = 0
         write(6,*) 'nstep reset to zero in rst_dyn()'
      endif

      pmax = vmax(u,pmin,im*jm*km)
      write(6,*) 'max u =', pmax, ' min u=', pmin

      pmax = vmax(v,pmin,im*jm*km)
      write(6,*) 'max v =', pmax, ' min v=', pmin

      pmax = vmax(pt,pmin,im*jm*km)
      write(6,*) 'max PT =', pmax, ' min PT=', pmin

      if(nq .ne. 0) then
	 read(iuic) q
      pmax = vmax(q,pmin,im*jm*km)
      write(6,*) 'max q =', pmax, ' min q=', pmin
      endif

      write(6,*) 'done reading ic for the dycore',nymd, nhms, nstep

	if(nstep .eq. 0) then
	do k=1,km

 	do j=1,jm
 	do i=1,im
 	q(i,j,k,1) = max(1.2e-12, q(i,j,k,1))
 	enddo
 	enddo

c     	call polavg(delp(1,1,k), im, jm)
c     	call polavg(pt(1,1,k), im, jm)
c     	call polavg(q(1,1,k,1), im, jm)

	enddo
	endif

      else

c Write

      rewind iuic
      write(iuic) nstep, nymd, nhms
      write(iuic) ps, delp,u,v,pt
      if(nq .ne. 0) write(iuic) q

      endif

      return
      end


c****6***0*********0*********0*********0*********0*********0**********72
      subroutine gmap(im, jm, nq,  akap,
     &         km,  pk3d_m,  pe3d_m,          u_m,  v_m,  pt_m,  q_m,
     &         kn,  pk3d_n,  pe3d_n, delp_n,  u_n,  v_n,  pt_n,  q_n  )
c****6***0*********0*********0*********0*********0*********0**********72

      implicit none
      
      integer im, jm
      integer km, kn, nq

C Input:
C original data km-level

      real      u_m(im,jm,km)
      real      v_m(im,jm,km)
      real     pt_m(im,jm,km)
      real      q_m(im,jm,km,*)
      real   pk3d_m(im,jm,km+1)
      real   pe3d_m(im,km+1,jm)


C Output:
C New data (kn-level)
      real   delp_n(im,jm,kn)
      real      u_n(im,jm,kn)
      real      v_n(im,jm,kn)
      real     pt_n(im,jm,kn)
      real      q_n(im,jm,kn,*)
      real   pk3d_n(im,jm,kn+1)
      real   pe3d_n(im,kn+1,jm)

      real ak(kn+1)
      real bk(kn+1)

c local (private)
      integer i, j, k, ic

      real pe0(im,km+1),pe1(im,km+1),pe2(im,kn+1),pe3(im,kn+1)
      real pk1(im,km+1),pk2(im,kn+1)
      real pkz(im,kn)

      real ptop,  pint
      integer ks
      real akap
      real ple,  pek, dak, bkh
      real undef
      real big
      parameter ( undef = 1.e15 )
      parameter (   big = 1.e10 )


      do 2000 j=1,jm
! new
      if ( kn .eq. 18 ) then
           call set_18(kn, ks, ptop, pint, ak, bk)
      elseif( kn .eq. 32 ) then
           call set_32(kn, ks, ptop, pint, ak, bk)
      elseif( kn .eq. 55 ) then
           call set_55(kn, ks, ptop, pint, ak, bk)
      else
!           write(6,*) 'No specific set_eta specified'
           call set_eta(kn, ks, ptop, pint, ak, bk)
      endif

c Copy original data to local 2D arrays.

      do k=1,km+1
      do i=1,im
      pe1(i,k) = pe3d_m(i,k,j)
      pk1(i,k) = pk3d_m(i,j,k)
      enddo
      enddo

	pek = ptop**akap
      do i=1,im
         pe0(i,   1) = pe1(i,1)
         pe3(i,   1) = ptop

         pe2(i,   1) = ptop
         pe2(i,kn+1) = pe1(i,km+1)

         pk2(i,   1) = pek
         pk2(i,kn+1) = pk1(i,km+1)
      enddo

      if(ks .ne. 0) then
      do 800 k=2,ks+1
         ple = ak(k)
         pek = ple ** akap
      do 800 i=1,im
      pe2(i,k) = ple
      pe3(i,k) = ple
800   pk2(i,k) = pek
      endif

      do 900 k=ks+2,kn
      ple = ak(k)
      pek = bk(k)
         do i=1,im
            pe2(i,k) = ple + pek*pe1(i,km+1)
            pk2(i,k) = pe2(i,k) ** akap
         enddo
900   continue

c map pt

      call map1_ppm( km,   pk1,   pt_m,
     &               kn,   pk2,   pt_n,
     &               im, 1, im, j, 1, jm, 1, 3, undef)

! Fix out of bound (beyond original ptop)
      do k=1, ks+1
         do i=1,im
            pkz(i,k) = 0.5*(pk2(i,k) + pk2(i,k+1))
         enddo
      enddo

      do k=ks,1,-1
         do i=1,im
            if( pt_n(i,j,k) .gt. big ) then
! Isothermal extension
                pt_n(i,j,k) = pt_n(i,j,k+1)*pkz(i,k+1)/pkz(i,k)
            endif
         enddo
      enddo

      if(ks .ne. 0) then
      do k=1,ks
         dak = ak(k+1) - ak(k)
         do i=1,im
           delp_n(i,j,k) = dak
         enddo
      enddo
      endif

      do k=ks+1,kn
         do i=1,im
           delp_n(i,j,k) = pe2(i,k+1) - pe2(i,k)
         enddo
      enddo

c Map constituents
      if(nq .ne. 0) then

      do ic=1,nq
      call map1_ppm( km,   pe1,   q_m(1,1,1,ic),
     &               kn,   pe2,   q_n(1,1,1,ic),
     &               im, 1, im, j, 1, jm, 0, 3, undef)

! Fix out of bound (beyond original ptop)
      do k=ks,1,-1
         do i=1,im
            if( q_n(i,j,k,ic) .gt. big ) then
                q_n(i,j,k,ic) = q_n(i,j,k+1,ic)
            endif
         enddo
      enddo

      enddo
      endif

c map U
      if(j .ne. 1) then

      do k=2,km+1
         do i=1,im
            pe0(i,k) = 0.5*(pe1(i,k)+pe3d_m(i,k,j-1))
         enddo
      enddo

      do k=ks+2,kn+1
      bkh = 0.5*bk(k)
         do i=1,im
            pe3(i,k) = ak(k) + bkh*(pe1(i,km+1)+pe3d_m(i,km+1,j-1))
         enddo
      enddo

      call map1_ppm( km,   pe0,   u_m,
     &               kn,   pe3,   u_n,
     &               im, 1, im, j, 1, jm, 1, 3, undef)
! Fix out of bound (beyond original ptop)
      do k=ks,1,-1
         do i=1,im
            if( u_n(i,j,k) .gt. big ) then
                u_n(i,j,k) = u_n(i,j,k+1)
            endif
         enddo
      enddo

      endif

c map v
      if(j .ne. 1 .and. j .ne. jm) then

      do k=2,km+1
      do i=2,im
      pe0(i,k) = 0.5*(pe1(i,k)+pe1(i-1,k))
      enddo
      enddo

c      i=1
      do k=2,km+1
      pe0(1,k) = 0.5*(pe1(1,k)+pe1(im,k))
      enddo

      do k=ks+2,kn+1
         do i=2,im
            pe3(i,k) = 0.5*(pe2(i,k)+pe2(i-1,k))
         enddo
      enddo
c i=1
      do k=ks+2,kn+1
         pe3(1,k) = 0.5*(pe2(1,k)+pe2(im,k))
      enddo

      call map1_ppm( km,   pe0,   v_m,
     &               kn,   pe3,   v_n,
     &               im, 1, im, j, 1, jm, 1, 3, undef)
! Fix out of bound (beyond original ptop)
      do k=ks,1,-1
         do i=1,im
            if( v_n(i,j,k) .gt. big ) then
                v_n(i,j,k) = v_n(i,j,k+1)
            endif
         enddo
      enddo
      endif

c Update PE and PKHT
      do k=1,kn+1
         do i=1,im
            pe3d_n(i,k,j) = pe2(i,k)
            pk3d_n(i,j,k) = pk2(i,k)
         enddo
      enddo

2000  continue

      return
      end

      real function vmax(a,pmin,n)
      implicit none
      integer n, i
      real pmin, pmax
      real a(*)
      pmax =a(n)
      pmin =a(n)
      do 10 i=1,n-1
      pmax = max(pmax,a(i))
      pmin = min(pmin,a(i))
10    continue
      vmax = pmax
      return
      end

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  gmean --- Calculate the mean of a 2D field
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      real function gmean(im,jm,jfirst,jlast,q)
!****6***0*********0*********0*********0*********0*********0**********72
! !USES:
      implicit none

! !INPUT PARAMETERS:
      integer im, jm                       ! Horizontal dimensions
      integer jfirst, jlast                ! Latitude strip
      real q(im,jfirst:jlast)              ! 2D field 
! !DESCRIPTION:
!     Calculate the mean of a 2D field
!
! !REVISION HISTORY:
!   WS  99.11.02:  Documentation; cleaning; jfirst:jlast
!
!EOP
!---------------------------------------------------------------------
!BOC

      integer i, j, j2, jmm1
      real dp, acap, area, gsum, xsum
      real, allocatable, save :: sine(:),cosp(:),sinp(:),cose(:)
      real dl
      logical first
      data first /.true./

      save acap,area
 
      j2   = max( 2, jfirst )
      jmm1 = min( jm-1, jlast )
      if(first) then
        allocate( sine(jm),cosp(jm),sinp(jm),cose(jm) )
        call setrig(im,jm,dp,dl,cosp,cose,sinp,sine)

c scaled polar cap area.
        acap = im*(1.+sine(2)) / dp

        area = 2.*acap
        do j=2,jm-1
          area = area + cosp(j)*im
        enddo
 
        first = .false.
      endif

      gsum = 0.0
      if ( jfirst .eq. 1  ) gsum = gsum + q(1,1)*acap
      if ( jlast  .eq. jm ) gsum = gsum + q(1,jm)*acap

      do j=j2,jmm1
        xsum = 0.
        do i=1,im
          xsum = xsum + q(i,j)
        enddo
        gsum = gsum + xsum*cosp(j)
      enddo
      gmean = gsum/area

      return
!EOC
      end
!---------------------------------------------------------------------


c****6***0*********0*********0*********0*********0*********0**********72
      subroutine setrig(im,jm,dp,dl,cosp,cose,sinp,sine)
c****6***0*********0*********0*********0*********0*********0**********72
      implicit none

      integer im, jm
      integer j, jm1
      real sine(jm),cosp(jm),sinp(jm),cose(jm)
      real dp, dl
      double precision pi, ph5

      jm1 = jm - 1
      pi  = 4.d0 * datan(1.d0)
      dl  = (pi+pi)/dble(im)
      dp  = pi/dble(jm1)

      do 10 j=2,jm
         ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(jm1))
10    sine(j) = dsin(ph5)

      cosp( 1) =  0.
      cosp(jm) =  0.

      do 80 j=2,jm1
80    cosp(j) = (sine(j+1)-sine(j)) / dp

c Define cosine at edges..

      do 90 j=2,jm
90    cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      cose(1) = cose(2)

      sinp( 1) = -1.
      sinp(jm) =  1.

      do 100 j=2,jm1
100   sinp(j) = 0.5 * (sine(j) + sine(j+1))

      write(6,*) 'SETRIG called. ',im,jm
      return
      end


      subroutine wrt3dr(iout,im,jm,km,a3,a2)
      implicit none
      integer iout, im, jm, km
      integer i, j, k
      real a3(im,jm,*)
      real*4 a2(im,jm)

      do 50 k=km,1,-1
      do 10 j=1,jm
      do 10 i=1,im
      if(abs(a3(i,j,k)) .lt. 1.e-25) then
         a2(i,j) = 0.
      else
         a2(i,j) = a3(i,j,k)
      endif
10    continue
      write(iout) a2
50    continue
      return
      end



!-----------------------------------------------------------------------
!BOP
! !ROUTINE: drymadj --- Total dry air mass
!
! !INTERFACE:

      subroutine drymadj( im, jm, km, jfirst, jlast,
     &                    moun, ak, bk, ps, delp, nq, q, id )

! !USES:

      implicit   none

! !INPUT PARAMETERS:
      integer im, jm, km     ! Dimensions
      integer jfirst, jlast  ! Latitude strip
      integer       id       ! 0:  checking total dry mass only
                             ! 1:  checking total dry mass and adjust
      integer nq             ! Number of tracers         
      logical moun
      real    ak(km+1)
      real    bk(km+1)

! !INPUT/OUTPUT PARAMETERS:
      real delp(im,jfirst:jlast,km)     ! pressure thickness
      real    q(im,jfirst:jlast,km,nq) 
      real   ps(im,jfirst:jlast)        ! surface pressure

! !DESCRIPTION:
!  Perform adjustment of the total dry-air-mass while preserving total
!  tracer mass
!  Developer: S.-J. Lin
!
! !REVISION HISTORY:
!   WS  99.10.26:  Revision; documentation; removed fvcore.h
!   WS  99.11.02:  Limited to jfirst:jlast
!   SJL 00.03.20:
!
!EOP
!---------------------------------------------------------------------
!BOC

      real psd(im,jfirst:jlast)     ! surface pressure  due to dry air mass

      real drym                     ! dry air mass in pascals
      parameter ( drym = 98222. )

      integer   i, j, k
      real  psmo
      real  psdry
      real  gmean
      real  dpd
      integer ic

! Check global maximum/minimum
      PSMO = gmean(im,jm,jfirst,jlast,ps)
      write(6,*) 
     &        'Total (moist) surface pressure before adjustment = ',
     &         PSMO

#if ( defined OpenMP )
!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j,k)
#endif

#if ( defined CRAY )
cmic$*       private(i,j,k)
#endif

#if ( defined IRIX64 ) || ( defined PGI )
c$doacross local(i,j,k)
#endif

      do 1000 j=jfirst,jlast
        do i=1,im
          psd(i,j) = ak(1)
        enddo

        if(nq .ne. 0) then
          do k=1, km
            do i=1,im
              psd(i,j) = psd(i,j) +  delp(i,j,k) * ( 1. - q(i,j,k,1) )
            enddo
          enddo
        else
          do k=1, km
            do i=1,im
              psd(i,j) = psd(i,j) +  delp(i,j,k)
            enddo
          enddo
        endif
1000  continue

      psdry = gmean(im,jm,jfirst,jlast,psd(1,jfirst))
      write(6,*) 'mean (dry) surface pressure = ', psdry
      if( id .eq. 0) return

      if(moun) then
        dpd = drym - psdry
      else
        dpd = 1000.*100. - psdry
      endif

      write(6,*) 'dry mass to be added (pascals) =', dpd

#if ( defined OpenMP )
!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j)
#endif

#if defined( CRAY ) && !defined( CRAY_T3E )
cmic$ do all private(i,j,ic)
#endif

#if ( defined IRIX64 )
c$doacross local(i,j, ic)
#endif

      do 2000 j=jfirst,jlast
        do i=1,im
           ps(i,j) = ps(i,j) + dpd
        enddo

        do ic=1,nq
           do i=1,im
              q(i,j,km,ic) =  q(i,j,km,ic)*delp(i,j,km) /
     &                      ( delp(i,j,km) + dpd )
           enddo
        enddo

        do i=1,im
           delp(i,j,km) = delp(i,j,km) + dpd
        enddo
2000  continue

      PSMO = gmean(im,jm,jfirst,jlast,ps)
      write(6,*) 
     &            'Total (moist) surface pressure after adjustment = ',
     &            PSMO

      return
!EOC
      end
!---------------------------------------------------------------------

c For CCM3 with 18 level.
      subroutine set_18(km,ks,ptop,pint,ak,bk)
      implicit none
      integer nl, ks, k, km
      parameter (nl = 18) 
      real ak(km+1),bk(km+1)
      real a(nl+1),b(nl+1)
      real ptop, pint

      data a /291.70,  792.92,  2155.39,  4918.34,  8314.25,
     &       7993.08, 7577.38,  7057.52,  6429.63,  5698.38,
     &       4879.13, 3998.95,  3096.31,  2219.02,  1420.39,
     &        754.13,  268.38,   0.0000,   0.0000 /

      data b /0.0000,  0.0000,   0.0000,   0.0000,   0.0000,
     &        0.0380541, 0.0873088, 0.1489307, 0.2232996,
     &        0.3099406, 0.4070096, 0.5112977, 0.6182465,
     &        0.7221927, 0.8168173, 0.8957590, 0.9533137,
     &        0.9851122, 1.0  /

      if( km .ne. nl) then
          write(6,*) 'Stop in set_18()'
          stop
      endif
 
      ptop =  a(1)
      ks = 4
      pint = a(ks+1) 

      do k=1,nl+1
      ak(k) = a(k)
      bk(k) = b(k)
      enddo

      return
      end

      subroutine set_32(km, ks, ptop, pint, ak, bk)

! New 32-layer setup with top at 0.4 mb for high horizontal resolution runs. 
! pint = 176.93 mb

      implicit none
      integer nl, ks, k, km
      parameter (nl = 32) 
      real ak(km+1)
      real bk(km+1)
      real a(nl+1),b(nl+1)
      real ptop, pint

!GrITAS*Vertical_Levels: 1000 850 700 500 400 300 250 200 150 100 70 50 30 10 5 2 1 0.4

      data a/40.00000,     100.00000,     200.00000,
     &      370.00000,     630.00000,    1000.00000,
     &     1510.00000,    2160.00000,    2900.00000,
     &     3680.00000,    4535.00000,    5505.00000,
     &     6607.26750,    7851.22980,    9236.56610,
     &    10866.34270,   12783.70000,   15039.30000,
     &    17693.00000,   20119.20876,   21686.49129,
     &    22436.28749,   22388.46844,   21541.75227,
     &    19873.78342,   17340.31831,   13874.44006,
     &    10167.16551,    6609.84274,    3546.59643,
     &     1270.49390,       0.00000,       0.00000   /

      data b/0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00696,       0.02801,
     &       0.06372,       0.11503,       0.18330,
     &       0.27033,       0.37844,       0.51046,
     &       0.64271,       0.76492,       0.86783,
     &       0.94329,       0.98511,       1.00000   /

      if( km .ne. nl) then
          write(6,*) 'Stop in set_eta(km=32)'
          stop
      endif
 
      ks = 18
      ptop =  a(1)
      pint = a(ks+1)

      do k=1,nl+1
         ak(k) = a(k)
         bk(k) = b(k)
      enddo

      return
      end

      subroutine set_55(km, ks, ptop, pint, ak, bk)

! Revised 55-level eta with pint at 176.93 mb
! SJL: 2000-03-20

      implicit none
      integer nl, ks, k, km
      parameter (nl = 55)
      real ak(km+1),bk(km+1)
      real a(nl+1)
      real b(nl+1)
      real ptop, pint

      data a/ 1.00000,       2.00000,       3.27000,
     &        4.75850,       6.60000,       8.93450,
     &       11.97030,      15.94950,      21.13490,
     &       27.85260,      36.50410,      47.58060,
     &       61.67790,      79.51340,     101.94420,
     &      130.05080,     165.07920,     208.49720,
     &      262.02120,     327.64330,     407.65670,
     &      504.68050,     621.68000,     761.98390,
     &      929.29430,    1127.68880,    1364.33920,
     &     1645.70720,    1979.15540,    2373.03610,
     &     2836.78160,    3380.99550,    4017.54170,
     &     4764.39320,    5638.79380,    6660.33770,
     &     7851.22980,    9236.56610,   10866.34270,
     &    12783.70000,   15039.30000,   17693.00000,
     &    20119.20876,   21686.49129,   22436.28749,
     &    22388.46844,   21541.75227,   19873.78342,
     &    17340.31831,   13874.44006,   10167.16551,
     &     6609.84274,    3546.59643,    1270.49390,
     &        0.00000,       0.00000   /

      data b / 0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00696,       0.02801,       0.06372,
     &         0.11503,       0.18330,       0.27033,
     &         0.37844,       0.51046,       0.64271,
     &         0.76492,       0.86783,       0.94329,
     &         0.98511,       1.00000  /


      if( km .ne. nl) then
          write(6,*) 'Stop in set_eta(km=55)'
          stop
      endif

         ks = 41
         ptop =  a(1)
         pint = a(ks+1)

      do k=1,NL+1
         ak(k) = a(k)
         bk(k) = b(k)
      enddo

      return
      end

      subroutine map1_ppm( km,   pe1,   q1,
     &                     kn,   pe2,   q2, 
     &       im, i1, i2, j, jfirst, jlast, iv, kord, undef)

! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

      implicit none

      real       r3, r23
      parameter (r3 = 1./3., r23 = 2./3.)

! Input:
      integer i1, i2
      integer im, jfirst, jlast, iv, kord
      integer km                             ! Original vertical dimension
      integer kn                             ! Target vertical dimension

      real pe1(im,km+1)
      real pe2(im,kn+1)
      real  q1(im,jfirst:jlast,km)
      real undef
! Output
      real  q2(im,jfirst:jlast,kn)

! Local
      real   dp1(i1:i2,km)
      real  q4(4,i1:i2,km)

      integer i, j, k, l, ll, k0
      real    pl, pr, qsum, delp, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) .lt. pe1(i,1)) then
!        q2(i,j,k) = q4(2,i,1)
         q2(i,j,k) = undef
         goto 555
      elseif(pe2(i,k) .ge. pe1(i,l) .and. pe2(i,k) .le. pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) .le. pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))
     .                    *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+
     .              q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*
     .               (r3*(1.+pl*(1.+pl))))
              do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) .gt. pe1(i,ll+1) ) then
! Whole layer..
                 qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                 else
                 delp = pe2(i,k+1)-pe1(i,ll)
                  esl = delp / dp1(i,ll)
                 qsum = qsum + delp*(q4(2,i,ll)+0.5*esl*
     .                 (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
                 k0 = ll
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

      return
      end

      subroutine ppm2m(a4, delp, km, i1, i2, iv, kord)

! iv = 0: positive definite scalars
! iv = 1: others

      implicit none
! Input
      integer km, lmt, iv
      integer i1, i2
      integer kord
      real    delp(i1:i2,km)
      real    a4(4,i1:i2,km)

! local arrays.
      real   dc(i1:i2,km)
      real   h2(i1:i2,km)
      real delq(i1:i2,km)

      real a1, a2, c1, c2, c3, d1, d2
      real qmax, qmin, cmax, cmin
      real qm, dq, tmp

! Local scalars:
      integer i, k, km1
      real qmp
      real lac
      integer it

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
            a4(4,i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
            c1 = (delp(i,k-1)+0.5*delp(i,k))/a4(4,i,k+1)
            c2 = (delp(i,k+1)+0.5*delp(i,k))/a4(4,i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /
     &                              (a4(4,i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
            dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
         enddo
      enddo

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

      do k=3,km1
      do i=i1,i2
      c1 = delq(i,k-1)*delp(i,k-1) / a4(4,i,k)
      a1 = a4(4,i,k-1) / (a4(4,i,k) + delp(i,k-1))
      a2 = a4(4,i,k+1) / (a4(4,i,k) + delp(i,k))
      a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(a4(4,i,k-1)+a4(4,i,k+1)) *
     &          ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -
     &                          delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
      d1 = delp(i,1)
      d2 = delp(i,2)
      qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
      dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
      c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
      a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
      a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
      dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
      cmax = max(a4(1,i,1), a4(1,i,2))
      cmin = min(a4(1,i,1), a4(1,i,2))
      a4(2,i,2) = max(cmin,a4(2,i,2))
      a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo

      if(iv .eq. 0) then
         do i=i1,i2
            a4(2,i,1) = max(0.,a4(2,i,1))
            a4(2,i,2) = max(0.,a4(2,i,2))
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
!****6***0*********0*********0*********0*********0*********0**********72
! No over- and undershoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
!****6***0*********0*********0*********0*********0*********0**********72
      enddo

      if(iv .eq. 0) then
      do i=i1,i2
         a4(2,i,km) = max(0.,a4(2,i,km))
         a4(3,i,km) = max(0.,a4(3,i,km))
      enddo
      endif

      do k=1,km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo
 
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
 
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord .eq. 7) then
!****6***0*********0*********0*********0*********0*********0**********72
! Huynh's 2nd constraint
!****6***0*********0*********0*********0*********0*********0**********72
      do k=2, km1
         do i=i1,i2
            h2(i,k) = delq(i,k) - delq(i,k-1)
         enddo
      enddo

      do 1000 k=3, km-2
        do i=i1,i2
! Right edges
         qmp   = a4(1,i,k)                 + 2.0*delq(i,k-1)
         lac   = a4(1,i,k) + 1.5*h2(i,k-1) + 0.5*delq(i,k-1)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
         qmp   = a4(1,i,k)                 - 2.0*delq(i,k)
         lac   = a4(1,i,k) + 1.5*h2(i,k+1) - 0.5*delq(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to prevent negatives
         if (iv .eq. 0) then
             call kmppm(dc(i1,k), a4(1,i1,k), it, 2)
         endif
1000  continue

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv .eq. 0) lmt = min(2, lmt)

      do k=3, km-2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      return
      end

      subroutine kmppm(dm, a4, im, lmt)
      implicit none

      real       r12
      parameter (r12 = 1./12.)

      integer im, lmt
      real    a4(4,*)
      real    dm(*)

      integer i
      real da1, da2, a6da
      real fmin

      if ( lmt .eq. 3 ) return
! Full constraint

      if(lmt .eq. 0) then
      do i=1,im
      if(dm(i) .eq. 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da .lt. -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da .gt. da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt .eq. 2) then
! Positive definite constraint

      do i=1,im
      if( abs(a4(3,i)-a4(2,i)) .lt. -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
      if( fmin .lt. 0. ) then
      if(a4(1,i).lt.a4(3,i) .and. a4(1,i).lt.a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
      elseif(a4(3,i) .gt. a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
      else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
      endif
      endif
      endif
      enddo

      elseif (lmt .eq. 1) then
              write(6,*) 'semi-monotonic not implemented!'
              stop
      endif

      return
      end

      subroutine set_eta(km, ks, ptop, pint, ak, bk)

!     use precision
      implicit none

! Choices for vertical resolutions are as follows:
! NCAR: 18, 26, and 30
! NASA DAO: smoothed version of CCM4's 30-level, 32, 55, 64, and 96 
! New 32-layer setup with top at 0.4 mb for high horizontal
! resolution runs. pint = 176.93 mb
! Revised 55-level eta with pint at 176.93 mb  SJL: 2000-03-20

! NCAR specific
      real a18(19),b18(19)              ! CCM3
      real a26(27),b26(27)              ! CCM4
      real a30(31),b30(31)              ! CCM4

! NASA only
      real a30m(31),b30m(31)            ! smoothed CCM4 30-L
      real a32(33),b32(33)
      real a32old(33),b32old(33)
      real a55(56),b55(56)
      real a55old(56),b55old(56)
      real a64(65),b64(65)
      real a96(97),b96(97)

      integer ks, k, km
      real ak(km+1),bk(km+1)
      real ptop                      ! model top (Pa)
      real pint                      ! transition to p (Pa)

! *** NCAR settings ***

      data a18 /291.70,  792.92,  2155.39,  4918.34,  8314.25,
     &         7993.08, 7577.38,  7057.52,  6429.63,  5698.38,
     &         4879.13, 3998.95,  3096.31,  2219.02,  1420.39,
     &         754.13,  268.38,   0.0000,   0.0000 /

      data b18 /0.0000,    0.0000,    0.0000,   0.0000,   0.0000,
     &          0.0380541, 0.0873088, 0.1489307, 0.2232996,
     &          0.3099406, 0.4070096, 0.5112977, 0.6182465,
     &          0.7221927, 0.8168173, 0.8957590, 0.9533137,
     &          0.9851122, 1.0  /
     
      data a26 /219.4067,  489.5209,   988.2418,   1805.201,
     &          2983.724,  4462.334,   6160.587,   7851.243,
     &          7731.271,  7590.131,   7424.086,   7228.744,
     &          6998.933,  6728.574,   6410.509,   6036.322,
     &          5596.111,  5078.225,   4468.96,    3752.191,
     &          2908.949,  2084.739,   1334.443,   708.499,
     &          252.136,   0.,         0. /

      data b26 /0.,         0.,         0.,         0.,
     &          0.,         0.,         0.,         0.,
     &          0.01505309, 0.03276228, 0.05359622, 0.07810627,
     &          0.1069411,  0.14086370, 0.180772,   0.227722,
     &          0.2829562,  0.3479364,  0.4243822,  0.5143168,
     &          0.6201202,  0.7235355,  0.8176768,  0.8962153,
     &          0.9534761,  0.9851122,  1.        /

      data a30 /225.523952394724, 503.169186413288, 1015.79474285245,
     &         1855.53170740604, 3066.91229343414,  4586.74766123295,
     &         6332.34828710556, 8070.14182209969,  9494.10423636436,
     &        11169.321089983,  13140.1270627975,  15458.6806893349,
     &        18186.3352656364, 17459.799349308,   16605.0657629967,
     &        15599.5160341263, 14416.541159153,   13024.8308181763,
     &        11387.5567913055,  9461.38575673103,  7534.44507718086,
     &         5765.89405536652, 4273.46378564835,  3164.26791250706,
     &         2522.12174236774, 1919.67375576496,  1361.80268600583,
     &          853.108894079924, 397.881818935275,    0., 
     &            0.  /

      data b30 /0.,                 0.,
     &          0.,                 0.,                0., 
     &          0.,                 0.,                0.,
     &          0.,                 0.,                0.,
     &          0.,                 0.,                0.03935482725501,
     &          0.085653759539127,  0.140122056007385, 0.20420117676258,
     &          0.279586911201477,  0.368274360895157, 0.47261056303978,
     &          0.576988518238068,  0.672786951065063, 0.75362843275070,
     &          0.813710987567902,  0.848494648933411, 0.88112789392471,
     &          0.911346435546875,  0.938901245594025, 0.96355980634689,
     &          0.985112190246582,  1.   /

! *** NASA DAO settings ***

! Smoothed CCM4's 30-Level setup
      data a30m / 300.00000,     725.00000,    1500.00000,
     &       2600.00000,    3800.00000,    5050.00000,
     &       6350.00000,    7750.00000,    9300.00000,
     &      11100.00000,   13140.00000,   15458.00000,
     &      18186.33580,   20676.23761,   22275.23783,
     &      23025.65071,   22947.33569,   22038.21991,
     &      20274.24578,   17684.31619,   14540.98138,
     &      11389.69990,    8795.97971,    6962.67963,
     &       5554.86684,    4376.83633,    3305.84967,
     &       2322.63910,    1437.78398,     660.76994,
     &          0.00000 /

      data b30m / 0.00000,       0.00000,       0.00000,
     &            0.00000,       0.00000,       0.00000,
     &            0.00000,       0.00000,       0.00000,
     &            0.00000,       0.00000,       0.00000,
     &            0.00000,       0.00719,       0.02895,
     &            0.06586,       0.11889,       0.18945,
     &            0.27941,       0.38816,       0.50692,
     &            0.61910,       0.70840,       0.77037,
     &            0.81745,       0.85656,       0.89191,
     &            0.92421,       0.95316,       0.97850,
     &            1.00000 /

      data a32/40.00000,     100.00000,     200.00000,
     &      370.00000,     630.00000,    1000.00000,
     &     1510.00000,    2160.00000,    2900.00000,
     &     3680.00000,    4535.00000,    5505.00000,
     &     6607.26750,    7851.22980,    9236.56610,
     &    10866.34270,   12783.70000,   15039.30000,
     &    17693.00000,   20119.20876,   21686.49129,
     &    22436.28749,   22388.46844,   21541.75227,
     &    19873.78342,   17340.31831,   13874.44006,
     &    10167.16551,    6609.84274,    3546.59643,
     &     1270.49390,       0.00000,       0.00000   /

      data b32/0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00696,       0.02801,
     &       0.06372,       0.11503,       0.18330,
     &       0.27033,       0.37844,       0.51046,
     &       0.64271,       0.76492,       0.86783,
     &       0.94329,       0.98511,       1.00000   /

      data a32old /300.0000,  454.1491,  652.5746,  891.9637, 1159.7102,
     &       1492.8248, 1902.5026, 2400.4835, 2998.6740, 3708.6584,
     &       4541.1041, 5505.0739, 6607.2675, 7851.2298, 9236.5661,
     &      10866.3427, 12420.400, 13576.500, 14365.400, 14807.800,  
     &       14915.500, 14691.400, 14129.400, 13214.800, 11923.200,
     &       10220.700,  8062.000,  5849.500,  3777.000,  2017.200,
     &         720.600,     0.000,     0.000 /

      data b32old /0.00, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.003633 , 0.014628 , 0.033276 , 0.060071 ,
     &        0.095722 , 0.141171 , 0.197623 , 0.266571 , 0.349839 ,
     &        0.449632 , 0.568589 , 0.685887 , 0.793252 , 0.883128 ,
     &        0.948792 , 0.9851119, 1.0000000 /

      data a55/ 1.00000,       2.00000,       3.27000,
     &        4.75850,       6.60000,       8.93450,
     &       11.97030,      15.94950,      21.13490,
     &       27.85260,      36.50410,      47.58060,
     &       61.67790,      79.51340,     101.94420,
     &      130.05080,     165.07920,     208.49720,
     &      262.02120,     327.64330,     407.65670,
     &      504.68050,     621.68000,     761.98390,
     &      929.29430,    1127.68880,    1364.33920,
     &     1645.70720,    1979.15540,    2373.03610,
     &     2836.78160,    3380.99550,    4017.54170,
     &     4764.39320,    5638.79380,    6660.33770,
     &     7851.22980,    9236.56610,   10866.34270,
     &    12783.70000,   15039.30000,   17693.00000,
     &    20119.20876,   21686.49129,   22436.28749,
     &    22388.46844,   21541.75227,   19873.78342,
     &    17340.31831,   13874.44006,   10167.16551,
     &     6609.84274,    3546.59643,    1270.49390,
     &        0.00000,       0.00000   /

      data b55 / 0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00000,       0.00000,       0.00000,
     &         0.00696,       0.02801,       0.06372,
     &         0.11503,       0.18330,       0.27033,
     &         0.37844,       0.51046,       0.64271,
     &         0.76492,       0.86783,       0.94329,
     &         0.98511,       1.00000  /

      data a55old /1.0000,    2.0000,    3.2700,    4.7585,     6.6000,
     &        8.9345,   11.9703,   15.9495,   21.1349,    27.8526,
     &       36.5041,   47.5806,   61.6779,   79.5134,   101.9442,
     &      130.0508,  165.0792,  208.4972,  262.0212,   327.6433,
     &      407.6567,  504.6805,  621.6800,  761.9839,   929.2943,
     &     1127.6888, 1364.3392, 1645.7072, 1979.1554,  2373.0361,
     &     2836.7816, 3380.9955, 4017.5417, 4764.3932,  5638.7938,
     &     6660.3377, 7851.2298, 9236.5661,10866.3427, 12420.400 ,
     &    13576.500 , 14365.400, 14807.800, 14915.500 , 14691.400, 
     &    14129.400 , 13214.800, 11923.200, 10220.700 ,  8062.000,
     &     5849.500 ,  3777.000,  2017.200,   720.600,      0.000,
     &        0.000 /

      data b55old /   0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     &        0.0000000, 0.003633 , 0.014628 , 0.033276 , 0.060071 ,
     &        0.095722 , 0.141171 , 0.197623 , 0.266571 , 0.349839 ,
     &        0.449632 , 0.568589 , 0.685887 , 0.793252 , 0.883128 ,
     &        0.948792 , 0.9851119, 1.0000000 /

      data a64/1.00000,       1.80   ,       2.80086,
     &       3.93309,       5.20139,       6.77626,
     &       8.69654,      10.99483,      13.81736,
     &      17.26058,      21.43286,      26.45448,
     &      32.45730,      39.58402,      47.98678,
     &      57.82525,      69.26401,      82.46925,
     &      97.60468,     114.82686,     135.08787,
     &     158.92390,     186.96575,     219.95555,
     &     258.76633,     304.42522,     358.14053,
     &     421.33383,     495.67748,     583.13893,
     &     686.03282,     807.08215,     949.49044,
     &    1117.02644,    1314.12387,    1545.99882,
     &    1818.78771,    2139.70974,    2517.25793,
     &    2961.42386,    3483.96212,    4098.70138,
     &    4821.91034,    5672.72831,    6673.67169,
     &    7851.22983,    9236.56613,   10866.34270,
     &   12783.69059,   15039.35130,   17693.01955,
     &   20814.92310,   23609.16397,   25271.17281,
     &   25844.93368,   25345.63084,   23760.05052,
     &   21046.23129,   17132.35351,   12832.00555,
     &    8646.27815,    5012.23907,    2299.34286,
     &     781.15294,       0.00000  /

      data b64/0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00000,       0.00000,
     &       0.00000,       0.00879,       0.03537,
     &       0.08047,       0.14526,       0.23147,
     &       0.34138,       0.47789,       0.61606,
     &       0.74456,       0.85318,       0.93300,
     &       0.97730,       1.00000  /

      data a96/  1.00000,       2.32782,       3.34990,
     &         4.49484,       5.62336,       6.93048,
     &         8.41428,      10.06365,      11.97630,
     &        14.18138,      16.70870,      19.58824,
     &        22.84950,      26.52080,      30.62845,
     &        35.19588,      40.24273,      45.78375,
     &        51.82793,      58.43583,      65.62319,
     &        73.40038,      81.77154,      90.73373,
     &       100.27628,     110.82243,     122.47773,
     &       135.35883,     149.59464,     165.32764,
     &       182.71530,     201.93164,     223.16899,
     &       246.63988,     272.57922,     301.24661,
     &       332.92902,     367.94348,     406.64044,
     &       449.40720,     496.67181,     548.90723,
     &       606.63629,     670.43683,     740.94727,
     &       818.87329,     904.99493,    1000.17395,
     &      1105.36304,    1221.61499,    1350.09326,
     &      1492.08362,    1649.00745,    1822.43469,
     &      2014.10168,    2225.92627,    2460.02905,
     &      2718.75195,    3004.68530,    3320.69092,
     &      3669.93066,    4055.90015,    4482.46240,
     &      4953.88672,    5474.89111,    6050.68994,
     &      6687.04492,    7390.32715,    8167.57373,
     &      9026.56445,    9975.89648,   11025.06934,
     &     12184.58398,   13466.04785,   14882.28320,
     &     16447.46289,   18177.25781,   20088.97461,
     &     21886.89453,   23274.16602,   24264.66602,
     &     24868.31641,   25091.15430,   24935.41016,
     &     24399.52148,   23478.13281,   22162.01758,
     &     20438.00586,   18288.83984,   15693.01172,
     &     12624.54199,    9584.35352,    6736.55713,
     &      4231.34326,    2199.57910,     747.11890,
     &        0.00000 /

      data b96/0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00000,       0.00000,       0.00000,
     &        0.00315,       0.01263,       0.02853,
     &        0.05101,       0.08030,       0.11669,
     &        0.16055,       0.21231,       0.27249,
     &        0.34169,       0.42062,       0.51005,
     &        0.61088,       0.70748,       0.79593,
     &        0.87253,       0.93400,       0.97764,
     &        1.00000 /

      select case (km)

! *** Original CCM3 18-Level setup ***
        case (18)
          ks = 4
          do k=1,km+1
            ak(k) = a18(k)
            bk(k) = b18(k)
          enddo

        case (26)
! CCM4 26-Level setup ***
          ks = 7
          do k=1,km+1
            ak(k) = a26(k)
            bk(k) = b26(k)
          enddo

        case (30)
! CCM4 30-Level setup ***
          ks = 12
          do k=1,km+1
            ak(k) = a30(k)
            bk(k) = b30(k)
          enddo

! *** Revised 32-L setup with ptop at 0.4 mb ***
        case (32)
          ks = 18
          do k=1,km+1
            ak(k) = a32(k)
            bk(k) = b32(k)
          enddo

! *** Revised 55-L setup with ptop at 0.01 mb ***
        case (55)
          ks = 41
          do k=1,km+1
            ak(k) = a55(k)
            bk(k) = b55(k)
          enddo

! *** Others ***
        case (64)
          ks = 51
          do k=1,km+1
            ak(k) = a64(k)
            bk(k) = b64(k)
          enddo

        case (96)
          ks = 77
          do k=1,km+1
            ak(k) = a96(k)
            bk(k) = b96(k)
          enddo

      end select

          ptop = ak(1)
          pint = ak(ks+1) 

      return
      end
