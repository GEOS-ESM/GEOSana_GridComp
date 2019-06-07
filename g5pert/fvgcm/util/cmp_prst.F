      program cmp_prst
      implicit none
      
      integer    imr, jmr, nlayr, nc, mx
      parameter  (imr   = 144)
      parameter  (jmr   =  90)
      parameter  (nlayr =  55)
      parameter  (nc    =   1)
      parameter  (mx    =   4)

      integer  jnp, nl
      parameter(jnp= jmr+1, nl=nlayr)

      real oro1(imr,jnp)        !CCM surface flag
      real pblht1(imr,jnp)      !pbl height (m?)
      real tpert1(imr,jnp)      !temperature pert.
      real qpert1(imr,nc,jnp)

      real srfrad1(imr,jnp)
      real asdir1(imr,jnp)      ! albedo: shortwave, direct
      real asdif1(imr,jnp)      ! albedo: shortwave, diffuse
      real aldir1(imr,jnp)      ! albedo: longwave, direct
      real aldif1(imr,jnp)      ! albedo: longwave, diffuse
      real lwup1(imr,jnp)       ! longwave up radiative flux
      real ts1(imr,jnp)         ! sfc temp (merged w/ocean if coupled)
      real tssub1(imr,mx,jnp)
      real sicthk1(imr,jnp)
      real snowh1(imr,jnp)      ! ccm snow depth (liquid water)

      real fsns1(imr,jnp)
      real precst1(imr,jnp)
      real precc1(imr,jnp)
      real precl1(imr,jnp)
      real flwds1(imr,jnp)

      real sols1(imr,jnp)
      real soll1(imr,jnp)
      real solsd1(imr,jnp)
      real solld1(imr,jnp)

      real oro2(imr,jnp)        !CCM surface flag
      real pblht2(imr,jnp)      !pbl height (m?)
      real tpert2(imr,jnp)      !temperature pert.
      real qpert2(imr,nc,jnp)

      real srfrad2(imr,jnp)
      real asdir2(imr,jnp)      ! albedo: shortwave, direct
      real asdif2(imr,jnp)      ! albedo: shortwave, diffuse
      real aldir2(imr,jnp)      ! albedo: longwave, direct
      real aldif2(imr,jnp)      ! albedo: longwave, diffuse
      real lwup2(imr,jnp)       ! longwave up radiative flux
      real ts2(imr,jnp)         ! sfc temp (merged w/ocean if coupled)
      real tssub2(imr,mx,jnp)
      real sicthk2(imr,jnp)
      real snowh2(imr,jnp)      ! ccm snow depth (liquid water)

      real fsns2(imr,jnp)
      real precst2(imr,jnp)
      real precc2(imr,jnp)
      real precl2(imr,jnp)
      real flwds2(imr,jnp)

      real sols2(imr,jnp)
      real soll2(imr,jnp)
      real solsd2(imr,jnp)
      real solld2(imr,jnp)

      real*4 buf(imr,jnp)

      integer ir1, ir2
      integer iw1, iw2
      integer recno
      logical outp

      data outp /.true./

      ir1 = 41
      open (ir1, file='p_rst1', form='unformatted', status='unknown')
      call rst_phys(0, ir1, imr, jnp, oro1, snowh1, lwup1,
     &              ts1, tssub1, mx, sicthk1,
     &              asdir1, asdif1, aldir1, aldif1,
     &              pblht1, tpert1, qpert1,  srfrad1,
     &              fsns1,  precst1, flwds1, precc1, precl1,
     &              sols1,  soll1,  solsd1,  solld1,
     &              nc)
      close (ir1)
      if (outp) then
        iw1 = 51
        recno = 1
        open (iw1, file='p_rst1.bin', form='unformatted',
     .        status='unknown', access='direct', recl=imr*jnp*4)
        call outgrads (iw1, imr, jnp, nc, mx, oro1, snowh1,
     .                 lwup1, ts1, tssub1, sicthk1, asdir1,
     .                 asdif1, aldir1, aldif1, pblht1, tpert1,
     .                 qpert1, srfrad1, fsns1, precst1, flwds1,
     .                 precc1, precl1, sols1, soll1, solsd1,
     .                 solld1, buf, recno)
        close (iw1)
      endif

      ir2 = 42
      open (ir2, file='p_rst2', form='unformatted', status='unknown')
      call rst_phys(0, ir2, imr, jnp, oro2, snowh2, lwup2,
     &              ts2, tssub2, mx, sicthk2,
     &              asdir2, asdif2, aldir2, aldif2,
     &              pblht2, tpert2, qpert2,  srfrad2,
     &              fsns2,  precst2, flwds2, precc2, precl2,
     &              sols2,  soll2,  solsd2,  solld2,
     &              nc)
      close (ir2)
      if (outp) then
        iw2 = 52
        recno = 1
        open (iw2, file='p_rst2.bin', form='unformatted',
     .        status='unknown', access='direct', recl=imr*jnp*4)
        call outgrads (iw2, imr, jnp, nc, mx, oro2, snowh2,
     .                 lwup2, ts2, tssub2, sicthk2, asdir2,
     .                 asdif2, aldir2, aldif2, pblht2, tpert2,
     .                 qpert2, srfrad2, fsns2, precst2, flwds2,
     .                 precc2, precl2, sols2, soll2, solsd2,
     .                 solld2, buf, recno)
        close (iw1)
      endif

      call diff2 ('ORO     ', imr, jnp, oro1, oro2)
      call diff2 ('PBLHT   ', imr, jnp, pblht1, pblht2)
      call diff2 ('TPERT   ', imr, jnp, tpert1, tpert2)
      call diff2 ('SRFRAD  ', imr, jnp, srfrad1, srfrad2)
      call diff2 ('ASDIR   ', imr, jnp, asdir1, asdir2)
      call diff2 ('ASDIF   ', imr, jnp, asdif1, asdif2)
      call diff2 ('ALDIR   ', imr, jnp, aldir1, aldir2)
      call diff2 ('ALDIF   ', imr, jnp, aldif1, aldif2)
      call diff2 ('LWUP    ', imr, jnp, lwup1, lwup2)
      call diff2 ('TS      ', imr, jnp, ts1, ts2)
      call diff2 ('SICTHK  ', imr, jnp, sicthk1, sicthk2)
      call diff2 ('SNOWH   ', imr, jnp, snowh1, snowh2)
      call diff2 ('FSNS    ', imr, jnp, fsns1, fsns2)
      call diff2 ('PRECST  ', imr, jnp, precst1, precst2)
      call diff2 ('PRECC   ', imr, jnp, precc1, precc2)
      call diff2 ('PRECL   ', imr, jnp, precl1, precl2)
      call diff2 ('FLWDS   ', imr, jnp, flwds1, flwds2)
      call diff2 ('SOLS    ', imr, jnp, sols1, sols2)
      call diff2 ('SOLL    ', imr, jnp, soll1, soll2)
      call diff2 ('SOLSD   ', imr, jnp, solsd1, solsd2)
      call diff2 ('SOLLD   ', imr, jnp, solld1, solld2)

      call diff3 ('QPERT   ', imr, nc, jnp, qpert1, qpert2)
      call diff3 ('TSSUB   ', imr, mx, jnp, tssub1, tssub2)

      stop
      end

      subroutine minmax (m, n, l, f, fmax, fmin)
      
      implicit         none

      integer          m, n, l
      integer          i, j, k

      real             f(m,n,l)
      real             fmax
      real             fmin

      fmax = f(1,1,1)
      fmin = f(1,1,1)

      do k = 1, l
        do j = 1, n
          do i = 1, m
            fmax = max(fmax,f(i,j,k))
            fmin = min(fmax,f(i,j,k))
          end do
        end do
      end do

      return
      end


      subroutine diff2 (name, im, jm, f1, f2) 

      implicit none

      character*8 name
      integer im, jm 
      integer i, j
      integer imax, jmax
      real f1(im,jm)
      real f2(im,jm)
      real dif, dmax

      dmax = 1.e-25
      imax = 0
      jmax = 0
      do j = 1, jm
        do i = 1, im
          dif = f2(i,j) - f1(i,j)
          if (abs(dif) .gt. dmax) then
            dmax = abs(dif)
            imax = i
            jmax = j
          endif 
        end do
      end do
      if (imax .ne. 0) then
        print *, name, ' max. diff at: ', imax, jmax, dmax,
     .            f1(imax,jmax), f2(imax,jmax)
      endif
       
      return
      end
  
      subroutine diff3 (name, im, km,  jm, f1, f2) 

      implicit none

      character*8 name
      integer im, km, jm 
      integer i, k, j
      integer imax, kmax, jmax
      real f1(im,km,jm)
      real f2(im,km,jm)
      real dif, dmax

      imax = 0
      kmax = 0
      jmax = 0
      dmax = 1.e-25
      do k = 1, km
        do j = 1, jm
          do i = 1, im
            dif = f2(i,k,j) - f1(i,k,j)
            if (abs(dif) .gt. dmax) then
              dmax = abs(dif)
              imax = i
              kmax = k
              jmax = j
            endif 
          end do
        end do
      end do
      if (imax .ne. 0) then
        print *, name, ' max. diff at: ', imax, kmax, jmax, dmax,
     .            f1(imax,kmax,jmax), f2(imax,kmax,jmax)
      endif
       
      return
      end

      subroutine rst_phys(id, iu, im, jm, oro, snowh, lwup,
     &                    ts,  tssub, mx, sicthk,
     &                    asdir, asdif, aldir, aldif,
     &                    pblht, tpert, qpert,  srfrad,
     &                    fsns,  precst, flwds, precc, precl,
     &                    sols,  soll,  solsd,  solld,
     &                    nq)
c****6***0*********0*********0*********0*********0*********0**********72
c
c Read/Write restart file for physical paremeterizations
c
c ID=0: read
c ID=1: write

c NCAR flags: 
c oro = 0  --> over ocean
c oro = 1  --> over land
c oro = 2  --> over sea ice


      implicit none
      integer id              ! flag to determine read/write
      integer iu              ! unit to read/write
      integer im              ! e-w dimension
      integer jm              ! n-s diemsnion
      integer mx
      integer nq
      integer  i,  j

      real oro(im, jm)        !CCM surface flag
      real pblht(im, jm)      !pbl height (m?)
      real tpert(im, jm)      !temperature pert.
      real qpert(im, nq, jm)

      real srfrad(im,jm)
      real asdir(im, jm)      ! albedo: shortwave, direct
      real asdif(im, jm)      ! albedo: shortwave, diffuse
      real aldir(im, jm)      ! albedo: longwave, direct
      real aldif(im, jm)      ! albedo: longwave, diffuse
      real lwup(im, jm)       ! longwave up radiative flux
      real ts(im, jm)         ! sfc temp (merged w/ocean if coupled)
      real tssub(im,mx,jm)
      real sicthk(im,jm)
      real snowh(im,jm)      ! ccm snow depth (liquid water)

      real fsns(im,jm)
      real precst(im,jm)
      real precc(im,jm)
      real precl(im,jm)
      real flwds(im,jm)

      real sols(im,jm)
      real soll(im,jm)
      real solsd(im,jm)
      real solld(im,jm)

      real vmax
      real pmax, pmin

      rewind iu

      if( id .eq. 0 ) then

        if(mx .ne. 4) then
           write(6,*) 'tssub dimension incorrect!  mx=', mx
           stop
        endif

            read(iu) oro
            read(iu) pblht
            read(iu) tpert
            read(iu) qpert
            read(iu) snowh
            read(iu) lwup

C SJL added 4/7/08
C For more flexable time stepping.
            read(iu) srfrad
            read(iu) fsns
            read(iu) precst
            read(iu) flwds
            read(iu) precc
            read(iu) precl
            read(iu) sols
            read(iu) soll
            read(iu) solsd
            read(iu) solld
C SJL added 4/7/08
c
            read(iu) ts
            read(iu) tssub
            read(iu) sicthk
c
            read(iu) asdir
            read(iu) asdif
            read(iu) aldir
            read(iu) aldif
      else
            write(iu) oro
            write(iu) pblht
            write(iu) tpert
            write(iu) qpert
            write(iu) snowh
            write(iu) lwup

C SJL added 4/7/08
C For more flexable time stepping.
            write(iu) srfrad
            write(iu) fsns
            write(iu) precst
            write(iu) flwds
            write(iu) precc
            write(iu) precl

            write(iu) sols
            write(iu) soll
            write(iu) solsd
            write(iu) solld
C SJL added 4/7/08
c
            write(iu) ts
            write(iu) tssub
            write(iu) sicthk
c
            write(iu) asdir
            write(iu) asdif
            write(iu) aldir
            write(iu) aldif
      endif

c Check Max/Min of pblht
c     pmax = vmax(pblht,pmin,im*jm)
c     write(6,*) 'max pbl (m)=', pmax, ' min =', pmin

c     pmax = vmax(ts,pmin,im*jm)
c     write(6,*) 'max ts=', pmax, ' min =', pmin

      return
      end

      subroutine outgrads (iw, im, jm, nq, mx, oro, snowh,
     .                     lwup, ts, tssub, sicthk, asdir,
     .                     asdif, aldir, aldif, pblht, tpert,
     .                     qpert, srfrad, fsns, precst, flwds,
     .                     precc, precl, sols, soll, solsd,
     .                     solld, buf, recno)
      implicit none

      integer iw, im, jm, nq, mx

      real oro(im, jm)        !CCM surface flag
      real pblht(im, jm)      !pbl height (m?)
      real tpert(im, jm)      !temperature pert.
      real qpert(im, nq, jm)

      real srfrad(im,jm)
      real asdir(im, jm)      ! albedo: shortwave, direct
      real asdif(im, jm)      ! albedo: shortwave, diffuse
      real aldir(im, jm)      ! albedo: longwave, direct
      real aldif(im, jm)      ! albedo: longwave, diffuse
      real lwup(im, jm)       ! longwave up radiative flux
      real ts(im, jm)         ! sfc temp (merged w/ocean if coupled)
      real tssub(im,mx,jm)
      real sicthk(im,jm)
      real snowh(im,jm)      ! ccm snow depth (liquid water)

      real fsns(im,jm)
      real precst(im,jm)
      real precc(im,jm)
      real precl(im,jm)
      real flwds(im,jm)

      real sols(im,jm)
      real soll(im,jm)
      real solsd(im,jm)
      real solld(im,jm)

      real*4 buf(im,jm)
      integer recno

      call out2d (iw, im, jm, oro, buf, recno)
      call out2d (iw, im, jm, pblht, buf, recno)
      call out2d (iw, im, jm, tpert, buf, recno)
      call out2d (iw, im, jm, srfrad, buf, recno)
      call out2d (iw, im, jm, asdir, buf, recno)
      call out2d (iw, im, jm, asdif, buf, recno)
      call out2d (iw, im, jm, aldir, buf, recno)
      call out2d (iw, im, jm, aldif, buf, recno)
      call out2d (iw, im, jm, lwup, buf, recno)
      call out2d (iw, im, jm, ts, buf, recno)
      call out2d (iw, im, jm, sicthk, buf, recno)
      call out2d (iw, im, jm, snowh, buf, recno)
      call out2d (iw, im, jm, fsns, buf, recno)
      call out2d (iw, im, jm, precst, buf, recno)
      call out2d (iw, im, jm, precc, buf, recno)
      call out2d (iw, im, jm, precl, buf, recno)
      call out2d (iw, im, jm, flwds, buf, recno)
      call out2d (iw, im, jm, sols, buf, recno)
      call out2d (iw, im, jm, soll, buf, recno)
      call out2d (iw, im, jm, solsd, buf, recno)
      call out2d (iw, im, jm, solld, buf, recno)

      call out3d (iw, im, nq, jm, qpert, buf, recno)
      call out3d (iw, im, mx, jm, tssub, buf, recno)

      return
      end

      subroutine out2d (iw, im, jm, f, buf, recno)

      implicit none 

      integer iw, im, jm
      integer i, j
      integer recno
      real f(im,jm)
      real*4 buf(im,jm)

      do j = 1, jm
        do i = 1, im
          buf(i,j) = f(i,j)
        end do
      end do

      write (iw, rec=recno) buf
      recno = recno + 1

      return
      end

      subroutine out3d (iw, im, km, jm, f, buf, recno)

      implicit none 

      integer iw, im, km, jm
      integer i, j, k
      integer recno
      real f(im,km,jm)
      real*4 buf(im,jm)

      do k = 1, km
        do j = 1, jm
          do i = 1, im
            buf(i,j) = f(i,k,j)
          end do
        end do
        write (iw, rec=recno) buf
        recno = recno + 1
      end do

      return
      end
