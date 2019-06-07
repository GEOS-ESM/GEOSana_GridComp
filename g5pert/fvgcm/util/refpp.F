      program refpp

c Program to convert fvgcm refout output to NCEP pressure surfaces.

c f90 -O -n32 -mips4 -cpp -s -o refpp refpp.F
c f77 -O -n32 -mips4 -bytereclen -s -o refpp refpp.F

c Last modified:  Tue Feb 23 16:40:53 EST 1999

#include "fvgcm.h"

      integer plev
      real undef
C
#if (FVGCM_LEV == 55)
      parameter (plev = 26)        ! NCEP data levels
#else
      parameter (plev = 17)        ! NCEP data levels
#endif
c     parameter (undef = 1.E15)
      parameter (undef = 1.E25)

C Original data
      integer imr,  jnp,  nl
      integer nc,  nq, iout, i, j, k, ic

C Set parameters for original data
      parameter (imr = FVGCM_LON, jnp = FVGCM_LAT, nl=FVGCM_LEV)
      parameter (nc = FVGCM_TRACER)

      integer      iargc
      integer      m, n, it, irec, nrec
      integer      ymd(2), hms(2)
      character*8  flags
      character*40 case
      character*40 type
      character*80 iflnm
      character*80 oflnm
      real         wgt2d, wgt, weight(imr,jnp,plev)
      logical      tmean
      logical      done
      logical      exist
 
      real  wz(IMR,JNP,NL+1)          ! geopotential
      real slp(IMR,JNP)               ! sea level pressure
      real  ps(IMR,JNP)               ! surface pressure (pascal)
      real  ua(IMR,JNP,NL)            ! u-wind on A-grid
      real  va(IMR,JNP,NL)            ! v-wind on A-grid
      real  ta(IMR,JNP,NL)            ! temperature in Kelvin
      real epv(IMR,JNP,NL)            ! Ertle Potential Vorticity
      real q(IMR,JNP,NL,nc)           ! consituents
      real omga(IMR,JNP,NL)           ! omega
      real tg(IMR,JNP)                ! ground temperature
      real precp(IMR,JNP)             ! precip

C Grid info
      real pe(imr,jnp,nl+1)   !could be p, p**kappa, or log(p)
      real*8 ak(nl+1)
      real*8 bk(nl+1)
      real*8 ptop, pint

C 4-byte work array
      REAL*4 wk2d(IMR,JNP)

C Interpolated data at NCEP levels
      real wzp(imr,jnp,plev)
      real up(imr,jnp,plev)
      real vp(imr,jnp,plev)
      real tp(imr,jnp,plev)
      real epvp(imr,jnp,plev)
      real qp(imr,jnp,plev,nc)
      real omgap(imr,jnp,plev)

C Monthly mean data at NCEP levels
      real wzm(imr,jnp,plev)
      real um(imr,jnp,plev)
      real vm(imr,jnp,plev)
      real tm(imr,jnp,plev)
      real epvm(imr,jnp,plev)
      real qm(imr,jnp,plev,nc)
      real omgam(imr,jnp,plev)
      real slpm(imr,jnp)
      real psm(imr,jnp)
      real tgm(imr,jnp)
      real precpm(imr,jnp)
c     real bkem(imr,jnp)

C Local work arrays
      real dp1(imr,jnp,nl)
      real QR(imr,jnp,nl)
      real QL(imr,jnp,nl)
      real Q6(imr,jnp,nl)

C NCEP levels
      real pm(plev)              ! pressure
      real p_out(plev)           ! pressure vertical coordinate
C p_out could be anything

#if (FVGCM_LEV == 55)
      data pm /0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 3., 5.,
     .          10.   , 20.   , 30.  , 50.  , 70.   ,100.   ,150.   ,
     .         200.   ,250.   ,300.  ,400.  ,500.   ,600.   ,700.   ,
     .         850.   ,925.   ,1000.  /
#else
      data pm / 10.   , 20.   , 30.  , 50.  , 70.   ,100.   ,150.   ,
     &         200.   ,250.   ,300.  ,400.  ,500.   ,600.   ,700.   ,
     &         850.   ,925.   ,1000.  /
#endif

      n = iargc()
      if (n .lt. 1 .or. n .gt. 2) then
        print *
        print *, 'Usage: refpp [-tm] filename'
        print *
        stop
      endif
c
c     iflnm: Input filename in the form of
c
c       CASE_TYPE_YYYYMMDD.hhmmss-YYYYMMDD.hhmmss
c
c      CASE: case name
c      TYPE: rout or diag
c

      call getarg (n, iflnm)
      call name_decode (iflnm, case, type, ymd, hms)

      tmean = .false.
      if (n .gt. 1) then
        call getarg (1, flags)
        if (index(flags, '-tm') .ne. 0) tmean = .true.
      endif

      nstr = len_trim(iflnm)
      if (tmean) then
        write (oflnm, '(a,a1,a,a4,a)') 
     .         case(:len_trim(case)), '_', type(:len_trim(type)),
     .         '_tm_', iflnm(nstr-30:nstr)
      else
        write (oflnm, '(a,a1,a,a5,a)') 
     .         case(:len_trim(case)), '_', type(:len_trim(type)),
     .         '_prs_', iflnm(nstr-30:nstr)
      endif

C Initialization


      nq = nc

      do m = 1, nq
        do k = 1, plev
          do j = 1, jnp
            do i = 1, imr
            qm(i,j,k,m) = 0.
            end do
          end do
        end do
      end do

      do k = 1, plev
        do j = 1, jnp
          do i = 1, imr
            um(i,j,k) = 0.
            vm(i,j,k) = 0.
            tm(i,j,k) = 0.
            epvm(i,j,k) = 0.
            wzm(i,j,k) = 0.
            weight(i,j,k) = 0.
            omgam(i,j,k) = 0.
          end do
        end do
      end do

      do j = 1, jnp
        do i = 1, imr
          slpm(i,j) = 0.
          psm(i,j) = 0.
          tgm(i,j) = 0.
          precpm(i,j) = 0.
        end do
      end do

C Set up the Eta-coordinate


      call set_eta(nl, ks, ptop, pint, ak, bk)

      iunit = 82         ! input data
      iout =  72         ! output data

      inquire (file=iflnm, exist=exist)
      if (.not. exist) then
        print *, 'Input file doesn''t exist, stop.'
        print *, '-->', iflnm
        stop
      endif


       open (iunit, file=iflnm, form='unformatted',           
     .       status='unknown', access='direct',recl=imr*jnp*4)
       nrec = 0


      open (iout, file=oflnm, form='unformatted',
     .      status='unknown', access='direct', 
     .      recl=imr*jnp*4)

C Time Loop

      done = .false.
      it = 0
      irec = 1
      do while (.not. done)

        call get_refout (wz,  slp,  ps,  ua,  va,  ta,  epv,
     .                   q, omga, tg, precp, imr, jnp, nl, nq, 
     .                   iunit, wk2d,  nrec, ios)

        if (ios .ne. 0) done = .true.
        
        if (.not. done) then

C compute pe
          do k=1,nl+1
            do j=1,jnp
              do i=1,imr
                pe(i,j,k) = ak(k) + bk(k)*ps(i,j)
C
C Use log-P for vertical interpolation
C
                pe(i,j,k) = log10(pe(i,j,k))
              end do
            end do
          end do

C Compute coordinate increments : thickness
          do k=1,nl
            do j=1,jnp
              do i=1,imr
                dp1(i,j,k) = pe(i,j,k+1) - pe(i,j,k)
              end do
            end do
          end do

          do k=1, plev
            p_out(k) = log10(pm(k)*1.E2)
c           write(6,*) k, pm(k), p_out(k)
          end do

C Geopotential height interpolation: given values defined at edges (pe)
C pe is log(pe)   (from model top to bottom surface)

          call int_lin(imr*jnp,nl+1,plev,1,imr*jnp, pe,  p_out, 
     &                 wz,   wzp ,  undef )

C interpolation: given cell-mean values

          call int_ppm(imr*jnp,nl,plev,1,imr*jnp, pe, p_out, 
     &                 ua,   up ,  3, dp1, ql, qr, q6, undef )
          call int_ppm(imr*jnp,nl,plev,1,imr*jnp, pe, p_out, 
     &                 va,   vp ,  3,  dp1, ql, qr, q6, undef )
          call int_ppm(imr*jnp,nl,plev,1,imr*jnp, pe, p_out, 
     &                 ta,   tp ,  3, dp1, ql, qr, q6, undef )
          call int_ppm(imr*jnp,nl,plev,1,imr*jnp, pe, p_out, 
     &                 epv,  epvp ,  3, dp1, ql, qr, q6, undef )

          do ic=1,nq
            call int_ppm(imr*jnp,nl,plev,1,imr*jnp, pe, p_out, 
     &                   q(1,1,1,ic),qp(1,1,1,ic) ,  3,dp1,
     &                   ql,qr,q6, undef )
          end do

          call int_ppm(imr*jnp,nl,plev,1,imr*jnp, pe, p_out, 
     &                 omga,  omgap,  3, dp1, ql, qr, q6, undef )

          if (tmean) then
            do k = 1, plev
              do j = 1, jnp
                do i = 1, imr
                  if (up(i,j,k) .ne. undef) then
                    weight(i,j,k) = weight(i,j,k) + 1.
                    um(i,j,k) = um(i,j,k) + up(i,j,k)
                    vm(i,j,k) = vm(i,j,k) + vp(i,j,k)
                    tm(i,j,k) = tm(i,j,k) + tp(i,j,k)
                    epvm(i,j,k) = epvm(i,j,k) + epvp(i,j,k)
                    omgam(i,j,k) = omgam(i,j,k) + omgap(i,j,k)
                    wzm(i,j,k) = wzm(i,j,k) + wzp(i,j,k)
                  endif
                end do
              end do
            end do

            do m = 1, nq
              do k = 1, plev
                do j = 1, jnp
                  do i = 1, imr
                    if (up(i,j,k) .ne. undef) then
                      qm(i,j,k,m) = qm(i,j,k,m) + qp(i,j,k,m)
                    endif
                  end do
                end do
              end do
            end do

            do j = 1, jnp
              do i = 1, imr
                psm(i,j) = psm(i,j) + ps(i,j)
                slpm(i,j) = slpm(i,j) + slp(i,j)
                tgm(i,j) = tgm(i,j) + tg(i,j)
                precpm(i,j) = precpm(i,j) + precp(i,j)
              end do
            end do

          else

c ...  No time mean, output right away

            call wrt3dr(IOUT,imr,jnp,plev,wzp,wk2d,irec)
            call wrt3dr(IOUT,imr,jnp,   1,slp,wk2d,irec)
            call wrt3dr(IOUT,imr,jnp,   1, ps,wk2d,irec)
            call wrt3dr(IOUT,imr,jnp,plev,up,wk2d,irec)
            call wrt3dr(IOUT,imr,jnp,plev,vp,wk2d,irec)
            call wrt3dr(IOUT,imr,jnp,plev,tp,wk2d,irec)
            call wrt3dr(IOUT,imr,jnp,plev,epvp,wk2d,irec)

            do ic=1, nq
              call wrt3dr(IOUT,imr,jnp,plev,qp(1,1,1,ic),wk2d,irec)
            end do


            call wrt3dr(IOUT,imr,jnp,   1, tg,wk2d,irec)

            call wrt3dr(IOUT,imr,jnp,   1, precp,wk2d,irec)

            call wrt3dr(IOUT,imr,jnp,plev,omgap,wk2d,irec)

          endif
 
          it = it + 1

        endif

      end do    ! time loop

      close (iunit)

      print *, 'it = ',  it

      if (tmean) then
        wgt2d = real(it)

        do k = 1, plev
          do j = 1, jnp
            do i = 1, imr
              if (weight(i,j,k) .ne. 0.) then
                um(i,j,k) = um(i,j,k) / weight(i,j,k)
                vm(i,j,k) = vm(i,j,k) / weight(i,j,k)
                tm(i,j,k) = tm(i,j,k) / weight(i,j,k)
                epvm(i,j,k) = epvm(i,j,k) / weight(i,j,k)
                wzm(i,j,k) = wzm(i,j,k) / weight(i,j,k)
                omgam(i,j,k) = omgam(i,j,k) / weight(i,j,k)
              else
                um(i,j,k) = undef
                vm(i,j,k) = undef
                tm(i,j,k) = undef
                epvm(i,j,k) = undef
                wzm(i,j,k) = undef
              endif
            end do
          end do
        end do

        do m = 1, nq
          do k = 1, plev
            do j = 1, jnp
              do i = 1, imr
                  if (weight(i,j,k) .ne. 0.) then
                    qm(i,j,k,m) = qm(i,j,k,m) / weight(i,j,k)
                else
                  qm(i,j,k,m) = undef
                endif
              end do
            end do
          end do
        end do

        do j = 1, jnp
          do i = 1, imr
            psm(i,j) = psm(i,j) / wgt2d
            slpm(i,j) = slpm(i,j) / wgt2d
            tgm(i,j) = tgm(i,j) / wgt2d
            precpm(i,j) = precpm(i,j) / wgt2d
          end do
        end do

C Write to IOUT :   Grads format

        irec = 1
        call wrt3dr(IOUT,imr,jnp,   1, psm,wk2d,irec)
        call wrt3dr(IOUT,imr,jnp,   1,slpm,wk2d,irec)
cx        call wrt3dr(IOUT,imr,jnp,   1, tgm,wk2d,irec)
cx        call wrt3dr(IOUT,imr,jnp,   1, precpm,wk2d,irec)
        call wrt3dr(IOUT,imr,jnp,plev,wzm,wk2d,irec)
        call wrt3dr(IOUT,imr,jnp,plev,um,wk2d,irec)
        call wrt3dr(IOUT,imr,jnp,plev,vm,wk2d,irec)
        call wrt3dr(IOUT,imr,jnp,plev,tm,wk2d,irec)
  
        do ic=1, nq
          call wrt3dr(IOUT,imr,jnp,plev,qm(1,1,1,ic),wk2d,irec)
        end do
  
        call wrt3dr(IOUT,imr,jnp,plev,epvm,wk2d,irec)
        call wrt3dr(IOUT,imr,jnp,plev,omgam,wk2d,irec)
  
        print *, 'it   = ', it
        print *, 'irec = ', irec
      endif

      close (iout)
      end

      subroutine int_lin(im,km,kn,i1,i2,pm,prk,q1,q2,
     &                   undef)
C
C This routine is mainly for geopotential interpolation;
C assuming geopotential is given from model top (k=1) to
C bottom surface (km+1), use linear interpolation for
C get a vertical column at a time.
C
C Input
      real pm(im,km)       ! coordinate at which q1 is defined
      real q1(im,km)       ! original data
      real undef           ! undefined value assigned if out of bound
C Output
      real q2(im,kn)       ! output interpolated data
      real prk(kn)         ! output coordinate
 
C local variables
      real s

      do 555 i=i1,i2
         k0 = 1
      do 500 k=1,kn

C Start searching
      if(prk(k) .lt. pm(i,1) .or. prk(k) .gt. pm(i,km) ) then
         q2(i,k) = undef
      else
        do 45 l=k0,km-1
        if(prk(k) .ge. pm(i,l) .and. prk(k) .le. pm(i,l+1)) then

              s  = (prk(k)-pm(i,L)) / (pm(i,L+1) -pm(i,L))
         q2(i,k) = q1(i,l) + s*(q1(i,L+1) - q1(i,L))
         k0 = l
         GOTO 500
        endif
45      continue
      endif
500   continue
555   continue

      return
      end

      subroutine int_ppm(im,km,kn,i1,i2,pe1,pp,q1,q2,
     &                   lmt, dp1, ql, qr, q6, undef)

C This routine applies a high order monotonic interpolation algorithm
C for "inerpolating" cell-mena values to prescribed coordinate
C Input:
      real pe1(im,km+1)      !coordinate value at original grid edges
      real q1(im,km)         !original field
      real pp(kn)            !output coordinate
      real undef
C Input work arrays
      real dp1(im,km),QR(im,km), QL(im,km), Q6(im,km)
C
C Output:
      real q2(im,kn)            ! output fields
C Local variabls:
      real qmax, qmin, qtmp
      real s

      do k=1,km
            do i=i1,i2
              dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            end do
      end do

C Construct the subgrid distribution of the original data
C accordingto the chosen lmt value
      call PPM2D(q1,QL,QR,Q6,dp1,im,km,i1,i2,lmt)

      do 555 i=i1,i2
         k0 = 1
      do 500 k=1,kn

C Start searching
      if(pp(k) .lt. pe1(i,1) .or. pp(k) .gt. pe1(i,km+1) ) then
         q2(i,k) = undef
      else
        do 45 l=k0,km
        if(pp(k) .le. pe1(i,L+1) .and. pp(k) .ge. pe1(i,l)) then

         if(L .eq. 1) then
            qmax = max(q1(i,L), q1(i,L+1))
            qmin = min(q1(i,L), q1(i,L+1))
         elseif(L .eq. km) then
            qmax = max(q1(i,L), q1(i,L-1))
            qmin = min(q1(i,L), q1(i,L-1))
         else
            qmax = max(q1(i,L-1),q1(i,L), q1(i,L+1))
            qmin = min(q1(i,L-1),q1(i,L), q1(i,L+1))
         endif

           s  = (pp(k)-pe1(i,L)) / dp1(i,L)
         qtmp = ql(i,l) + s*(qr(i,l) - ql(i,l) +
     &                       q6(i,l)*(1.-s)     )

         qtmp    = max(qtmp, qmin)
         q2(i,k) = min(qtmp, qmax)
         k0 = l
         GOTO 500
        endif
45      continue
      endif
500   continue
555   continue

      return
      end

      subroutine fmap_ppm(im,km,kn,i1,i2,pe1,prk,q1,q2,
     &          dp1, ql, qr, q6, fac, undef)

C This routine is currently not used; itis for the interface between
C the analysis and the model
C Partial mapping using the piecewise parabolic finite-volume method

      parameter (R3 = 1./3., R23 = 2./3.)
C Input
      real pe1(im,km+1)      ! Cell edges of the roiginal coordinate
      real q1(im,km)         ! original cell-mena values 
      real prk(kn)           ! cell center coordinate
      real undef
      real fac               ! output cell thickness ( 0 < fac < 0.5)

C Output
      real q2(im,kn)         ! output partially mapped data
      real DPK

C Local work arrays
      real dp1(im,km),QR(im,km), QL(im,km), Q6(im,km)
C
      do k=1,kn
      do i=i1,i2
         q2(i,k) = undef
      end do
      end do

C fac should be < 0.5

      do k=1,km
            do i=i1,i2
              dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            end do
      end do

C Construct the monotonic subgrid distribution of the original data
      call PPM2D(q1,QL,QR,Q6,dp1,im,km,i1,i2,0)

      do 555 i=i1,i2
         k0 = 1
      do 555 k=1,kn

      if(prk(k) .gt. pe1(i,km+1)) go to 555  

         if(k .eq. 1) then
            dpk =  fac*min( prk(1)-pe1(i,1), prk(2)-prk(1) )
         elseif(k .eq. kn) then
            dpk = fac*min(prk(kn)-prk(kn-1), pe1(i,km) - prk(kn))
         else
            dpk = fac*min( prk(k)-prk(k-1), prk(k+1)-prk(k) )
         endif

         p_hi = prk(k) - DPK
         p_lo = prk(k) + DPK

      if(p_lo .gt. pe1(i,km+1)) then
         p_lo = pe1(i,km+1)
	   p_hi = prk(k) - max(p_lo - prk(k), 0.1)
      endif
         dp2 = p_lo - p_hi

         do 45 L=k0,km
         if( p_hi.ge.pe1(i,L) .and. p_hi.le. pe1(i,L+1)) then
             k0 = L
             PL = (p_hi-pe1(i,L)) / dp1(i,L)

              
             if(p_lo .le. pe1(i,L+1)) then
C within one cell
               PR = (p_lo-pe1(i,L)) / dp1(i,L)
               TT = R3*(PR*(PR+PL)+PL**2)
               q2(i,k) = QL(i,L) + 0.5*(Q6(i,L)+QR(i,L)
     &                   - QL(i,L))*(PR+PL) - Q6(i,L)*TT
              GOTO 555
             else
C
C Fractional area...
C
              delp = pe1(i,L+1) - p_hi
              TT   = R3*(1.+PL*(1.+PL))
              qsum = delp*(QL(i,L)+0.5*(Q6(i,L)+
     &               QR(i,L)-QL(i,L))*(1.+PL)-Q6(i,L)*TT)
              k1 = L + 1
             GOTO 111
             endif
         endif
45       continue

111      continue
         do 55 L=k1,km
         if( p_lo .gt. pe1(i,L+1) ) then
C
C Whole layer..
C
            qsum = qsum + dp1(i,L)*q1(i,L)
         else
           delp = p_lo - pe1(i,L)
           esl  = delp / dp1(i,L)
           qsum = qsum + delp * ( QL(i,L)+0.5*esl*
     &            (QR(i,L)-QL(i,L)+Q6(i,L)*(1.-R23*esl)) )
           k0 = L
           GOTO 123
         endif
55       continue
123   q2(i,k) = qsum / dp2
555   continue

      return
      end

C****6***0*********0*********0*********0*********0*********0**********72
      subroutine PPM2D(P,AL,AR,A6,delp,im,km,i1,i2,LMT)
C****6***0*********0*********0*********0*********0*********0**********72
      integer imr,  jnp,  nl
      parameter (imr = FVGCM_LON, jnp = FVGCM_LAT, nl=FVGCM_LEV)
      REAL P(im,km),AL(im,km),AR(im,km),A6(im,km),
     &     delp(im,km)
      REAL DC(imr*jnp,nl), delq(imr*jnp,nl)
 
      km1 = km - 1
 
      do 500 k=2,km
      do 500 i=i1,i2
500   A6(i,k) = delp(i,k-1) + delp(i,k)

      do 1000 k=1,km1
      do 1000 i=i1,i2
      delq(i,k) = P(i,k+1) - P(i,k)
1000  continue
 
      DO 1220 k=2,km1
      DO 1220 i=i1,i2
      c1 = (delp(i,k-1)+0.5*delp(i,k))/A6(i,k+1)
      c2 = (delp(i,k+1)+0.5*delp(i,k))/A6(i,k)
      tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /
     &                              (A6(i,k)+delp(i,k+1))
      Qmax = max(P(i,k-1),P(i,k),P(i,k+1)) - P(i,k)
      Qmin = P(i,k) - min(P(i,k-1),P(i,k),P(i,k+1))
      DC(i,k) = sign(min(abs(tmp),Qmax,Qmin), tmp)
1220  CONTINUE
 
C****6***0*********0*********0*********0*********0*********0**********72
C 4th order interpolation of the provisional cell edge value
C****6***0*********0*********0*********0*********0*********0**********72
 
      DO 12 k=3,km1
      DO 12 i=i1,i2
      c1 = delq(i,k-1)*delp(i,k-1) / A6(i,k)
      A1 = A6(i,k-1) / (A6(i,k) + delp(i,k-1))
      A2 = A6(i,k+1) / (A6(i,k) + delp(i,k))
      AL(i,k) = P(i,k-1) + c1 + 2./(A6(i,k-1)+A6(i,k+1)) *
     &          ( delp(i,k  )*(c1*(A1 - A2)+A2*DC(i,k-1)) -
     &                          delp(i,k-1)*A1*DC(i,k  ) )
12    CONTINUE
 
 
C three-cell parabolic subgrid distribution at model top
 
      DO 10 i=i1,i2
C three-cell PP-distribution
C Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
C a3 = a / 3
C b2 = b / 2
      S1 = delp(i,1)
      S2 = delp(i,2) + S1
C
      S3 = delp(i,2) + delp(i,3)
      S4 = S3 + delp(i,4)
      SS3 =  S3 + S1
      S32 = S3*S3
      S42 = S4*S4
      S34 = S3*S4
C model top
      a3 = (delq(i,2) - delq(i,1)*S3/S2) / (S3*SS3)
C
      if(abs(a3) .gt. 1.E-14) then
         b2 =  delq(i,1)/S2 - a3*(S1+S2)
         SC = -b2/(3.*a3)
         if(SC .lt. 0. .or. SC .gt. S1) then
             AL(i,1) = P(i,1) - S1*(a3*S1 + b2)
         else
             AL(i,1) = P(i,1) - delq(i,1)*S1/S2
         endif
      else
C Linear
         AL(i,1) = P(i,1) - delq(i,1)*S1/S2
      endif
      DC(i,1) = P(i,1) - AL(i,1)
C compute coef. for the off-centered area preserving cubic poly.
      DM = delp(i,1) / (S34*SS3*(delp(i,2)+S3)*(S4+delp(i,1)))
      F1 = delp(i,2)*S34 / ( S2*SS3*(S4+delp(i,1)) )
      F2 = (delp(i,2)+S3) * (SS3*(delp(i,2)*S3+S34+delp(i,2)*S4)
     &      + S42*(delp(i,2)+S3+S32/S2))
      F3 = -delp(i,2)*( SS3*(S32*(S3+S4)/(S4-delp(i,2))
     &      + (delp(i,2)*S3+S34+delp(i,2)*S4))
     &      + S42*(delp(i,2)+S3) )
      F4 = SS3*delp(i,2)*S32*(delp(i,2)+S3) / (S4-delp(i,2))
      AL(i,2) = F1*P(i,1)+(F2*P(i,2)+F3*P(i,3)+F4*P(i,4))*DM
C****6***0*********0*********0*********0*********0*********0**********72
C No over- and undershoot condition
      Cmax = max(P(i,1), P(i,2))
      Cmin = min(P(i,1), P(i,2))
      AL(i,2) = max(Cmin,AL(i,2))
      AL(i,2) = min(Cmax,AL(i,2))
10    continue

C****6***0*********0*********0*********0*********0*********0**********72

C Bottom
C Area preserving cubic with 2nd deriv. = 0 at the surface
      DO 15 i=i1,i2
      d1 = delp(i,km )
      d2 = delp(i,km1)
      qm = (d2*P(i,km)+d1*P(i,km1)) / (d1+d2)
      dq = 2.*(P(i,km1)-P(i,km)) / (d1+d2)
      c1 = (AL(i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
      AL(i,km) = qm - c1*d1*d2*(d2+3.*d1)
      AR(i,km) = d1*(8.*c1*d1**2-c3) + AL(i,km)
      DC(i,km) = AR(i,km) -  P(i,km)
C****6***0*********0*********0*********0*********0*********0**********72
C No over- and undershoot condition
      Cmax = max(P(i,km), P(i,km1))
      Cmin = min(P(i,km), P(i,km1))
      AL(i,km) = max(Cmin,AL(i,km))
      AL(i,km) = min(Cmax,AL(i,km))
C****6***0*********0*********0*********0*********0*********0**********72
15    continue

      do 20 k=1,km1
      do 20 i=i1,i2
      AR(i,k) = AL(i,k+1)
20    continue
C
C f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
C
      do 30 k=1,km
      do 30 i=i1,i2
      A6(i,k) = 3.*(P(i,k)+P(i,k) - (AL(i,k)+AR(i,k)))
30    continue
 
      if(LMT.LE.2) then
	len = i2-i1+1
         do k=1,km
         call lmppm(DC(i1,k),A6(i1,k),AR(i1,k),AL(i1,k),
     &               P(i1,k),len,LMT)
         end do
      endif
      return
      end

C****6***0*********0*********0*********0*********0*********0**********72
      subroutine lmppm(DM,A6,AR,AL,P,IM,LMT)
C****6***0*********0*********0*********0*********0*********0**********72
      parameter ( R12 = 1./12. )
      REAL A6(IM),AR(IM),AL(IM),P(IM),DM(IM)
 
C LMT = 0: full monotonicity
C LMT = 1: semi-monotonic constraint (no undershoot)
C LMT = 2: positive-definite constraint
 
      if(LMT.eq.0) then
C Full constraint
      do 100 i=1,IM
      if(DM(i).eq.0.) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      else
      da1  = AR(i) - AL(i)
      da2  = da1**2
      A6DA = A6(i)*da1
      if(A6DA .lt. -da2) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      elseif(A6DA .gt. da2) then
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
      endif
100   continue
      elseif(LMT.eq.1) then
C Semi-monotonic constraint
      do 150 i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) go to 150
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
150   continue
      elseif(LMT.eq.2) then
C Positive definite constraint
      do 250 i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) go to 250
      fmin = p(i) + 0.25*(AR(i)-AL(i))**2/A6(i) + A6(i)*R12
      if(fmin.ge.0.) go to 250
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
250   continue
      endif
      return
      end

      subroutine get_refout(wz,   slp,  ps,  ua,  va,  ta,  epv,
     &    q, omga,  tg, precp, imr, jnp, nl, nq, iout, wk2d, nrec, ios)

      implicit none

      integer imr,  jnp,  nl
      integer nq, iout, i, j, k, ic
      integer ios, nrec
      real q(IMR,JNP,NL,nq)
      real epv(IMR,JNP,NL)
      real  ua(IMR,JNP,NL),va(IMR,JNP,NL), ps(IMR,JNP)
     &     ,ta(IMR,JNP,NL),wz(IMR,JNP,NL+1)
     &     ,slp(IMR,JNP), omga(IMR,JNP,NL)
      real tg(IMR,JNP)
      real precp(IMR,JNP)

      REAL*4 wk2d(IMR,JNP)

!      nrec=0
      call r3d(IOUT,nrec, IMR,JNP,NL+1,WZ,wk2d,ios)
      call r3d(IOUT,nrec, IMR,JNP,1,slp,wk2d,ios)
      call r3d(IOUT,nrec, IMR,JNP,1,PS,wk2d,ios)
      call r3d(IOUT,nrec, IMR,JNP,NL,ua,wk2d,ios)
      call r3d(IOUT,nrec, IMR,JNP,NL,va,wk2d,ios)
      call r3d(IOUT,nrec, IMR,JNP,NL,TA,wk2d,ios)
      call r3d(IOUT,nrec, IMR,JNP,NL,epv,wk2d,ios)

      if(NQ.ne.0) then
      do IC=1,NQ
         call r3d(IOUT,nrec, IMR,JNP,NL,q(1,1,1,ic),wk2d,ios)
      end do

      endif
 
C surface temp.
       call r3d(IOUT,nrec, IMR,JNP,1,tg,wk2d,ios)

C precip
       call r3d(IOUT,nrec, IMR,JNP,1,precp,wk2d,ios)
c omega
      call r3d(IOUT,nrec, IMR,JNP,NL,omga,wk2d,ios)

      return
      end


      subroutine r3d(IOUT,nrec, im,jm,km,a3,a2,ios)
      implicit none
      integer iout, im, jm, km
      integer i, j, k
      integer ios, nrec
      real a3(im,jm,*)
      real*4 a2(im,jm)

      do 50 k=1,km
         nrec=nrec+1
         read(IOUT,rec=nrec, iostat=ios) a2
      do 10 j=1,jm
      do 10 i=1,im
         a3(i,j,k) = a2(i,j)
10    continue
50    continue
      return
      end

      subroutine wrt3dr(IOUT,im,jm,km,a3,a2,irec)
      implicit none
      integer iout, im, jm, km
      integer i, j, k
      integer irec
      real a3(im,jm,*)
      real*4 a2(im,jm)

      do 50 k=km,1,-1
      do 10 j=1,jm
      do 10 i=1,im
      if(abs(a3(i,j,k)) .lt. 1.E-25) then
         a2(i,j) = 0.
      else
         a2(i,j) = a3(i,j,k)
      endif
10    continue
      write(IOUT,rec=irec) a2
      irec = irec + 1
50    continue
      return
      end

      integer function len_trim (string)

************************************************************************
*                                                                      * 
*     len_trim.f - return the length of string without the trailing    * 
*                  blanks and tabs                                     * 
*                                                                      * 
*     Last Modified:  Thu Mar 13 17:47:34 EST 1997                     * 
*                                                                      * 
************************************************************************

      implicit         none

      integer          tab, blank
      parameter        (tab   =  9)
      parameter        (blank = 32)

      character*(*)    string
      integer          n, code

      len_trim = 0
      n = len(string)
      if (n .eq. 0) return

      do while (n .gt. 0)
        code = ichar(string(n:n))
        if (code .eq. blank  .or. code .eq. tab) then
          len_trim = n - 1
	  n = n - 1
        else
	  len_trim = n
          n = 0
        end if
      end do

      return 
      end

      subroutine name_decode (flnm, case, type, ymd, hms)

      implicit none

      character*(*)     flnm
      character*(*)     case
      character*(*)     type

      integer           index
      integer           len_trim
      integer           ymd(2)
      integer           hms(2)
      integer           nc, n
c
      nc = len_trim(flnm)
c
c     flnm: Filename in the form of
c
c           CASE_TYPE_YYYYMMDD.hhmmss-YYYYMMDD.hhmmss
c
      read(flnm(nc-5:nc),     '(i6.6)') hms(2)
      read(flnm(nc-14:nc-7),  '(i8.8)') ymd(2)
      read(flnm(nc-21:nc-16), '(i6.6)') hms(1)
      read(flnm(nc-30:nc-23), '(i8.8)') ymd(1)

      n = index(flnm(:nc-32), '_')

      case = flnm(:n-1)
      type = flnm(n+1:nc-32)
      
      return
      end

