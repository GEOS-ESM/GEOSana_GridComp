module sw_core
 
! This module contains vertical independent part of the Lagrangian dynamics

contains

 subroutine c_sw(u,      v,       pt,     delp,                       &
                 uc,     vc,      ptc,    delpf,     ptk,             &
                 cosp,   acosp,   cose,   coslon,   sinlon,           &
                 dxdt,   dxe,     dtdx2,                              &
                 dtdx4,  dtxe5,   rdxe,   dycp,     dydt,             & 
                 dtdy5,  cye,     fc,     ifax,     trigs,            &
                 dc,     sc,      zt_c,   tiny,     rcap,             &
                 im,     jm,      jfirst, jlast,    ng_c,             &
                 ng_d,   ng_s,    js2g0,  jn2g0,    js2gc,            &
                 jn1gc,  iord,    jord,   cosl5,    sinl5,            &
		 c_sw_tape_rec                                        &
		 )

! Routine for shallow water dynamics on the C-grid

! !USES:

  use precision
  use tp_core

  implicit none

! INPUT:
  integer, intent(in):: im
  integer, intent(in):: jm
  integer, intent(in):: jfirst
  integer, intent(in):: jlast
  integer, intent(in):: js2g0
  integer, intent(in):: jn2g0
  integer, intent(in):: js2gc
  integer, intent(in):: jn1gc
  integer, intent(in):: iord
  integer, intent(in):: jord
  integer, intent(in):: ng_c
  integer, intent(in):: ng_s
  integer, intent(in):: ng_d

  integer :: c_sw_tape_rec                  ! FastOpt

! polar filter related input arrays:
  integer, intent(in)::  ifax(13)                      !ECMWF fft
  real(r8), intent(in):: trigs(3*im/2+1)
  real(r8), intent(in):: dc(im,js2g0:jn2g0)
  real(r8), intent(in):: sc(js2g0:jn2g0)

! Prognostic variables:
  real(r8), intent(in):: u(im,jfirst-ng_d:jlast+ng_s)
  real(r8), intent(in):: v(im,jfirst-ng_s:jlast+ng_d)  ! Wind in Y
  real(r8), intent(in):: pt(im,jfirst-ng_d:jlast+ng_d)   ! Wind in Y
  real(r8), intent(in):: delp(im,jfirst:jlast)         ! Delta pressure
  real(r8), intent(in):: delpf(im,jfirst-ng_d:jlast+ng_d)

  real(r8), intent(in):: fc(js2gc:jn1gc)

  real(r8), intent(in):: cosp(jm)
  real(r8), intent(in):: acosp(jm)
  real(r8), intent(in):: cose(jm)

  real(r8), intent(in):: dxdt(jm)
  real(r8), intent(in):: dxe(jm)
  real(r8), intent(in):: rdxe(jm)
  real(r8), intent(in):: dtdx2(jm)
  real(r8), intent(in):: dtdx4(jm)
  real(r8), intent(in):: dtxe5(jm)
  real(r8), intent(in):: dycp(jm)
  real(r8), intent(in)::  cye(jm)

  real(r8), intent(in):: sinlon(im)
  real(r8), intent(in):: coslon(im)
  real(r8), intent(in):: sinl5(im)
  real(r8), intent(in):: cosl5(im)

  real(r8), intent(in):: zt_c
  real(r8), intent(in):: rcap
  real(r8), intent(in):: tiny
  real(r8), intent(in):: dydt
  real(r8), intent(in):: dtdy5

! Output:
  real(r8), intent(out):: uc(im,jfirst-ng_d:jlast+ng_d)
  real(r8), intent(out):: vc(im,jfirst-2   :jlast+2 )
  real(r8), intent(out):: ptc(im,jfirst:jlast)
  real(r8), intent(out):: ptk(im,jfirst:jlast)

!--------------------------------------------------------------
! Local 
    real(r8)   fx(im,jfirst:jlast)
    real(r8)  xfx(im,jfirst:jlast)
    real(r8)  tm2(im,jfirst:jlast)

    real(r8)   va(im,jfirst-1:jlast)

    real(r8)  wk1(im,jfirst-1:jlast+1)
    real(r8)  cry(im,jfirst-1:jlast+1)
    real(r8)   fy(im,jfirst-1:jlast+1)

    real(r8) ymass(im,jfirst: jlast+1) 
    real(r8)   yfx(im,jfirst: jlast+1)

    real(r8)  crx(im,jfirst-max(1,ng_c):jlast+max(2,ng_c))   ! FastOpt
    real(r8)   u2(im,jfirst-ng_d:jlast+ng_d)
    real(r8)   v2(im,jfirst-ng_d:jlast+ng_d)

    real(r8) fxj(im)
    real(r8) p1d(im)
    real(r8) cx1(im)

    real(r8) qtmp(-im/3:im+im/3)
    real(r8) slope(-im/3:im+im/3)
    real(r8) al(-im/3:im+im/3)
    real(r8) ar(-im/3:im+im/3)
    real(r8) a6(-im/3:im+im/3)

    real(r8) us, vs, un, vn
    real(r8) p1ke, p2ke

    logical ffsl(jm)
    logical sld

    integer i, j, im2
    integer js1g1
    integer js2g1
    integer js2gc1
    integer js2gcp1
    integer jn2gc
    integer jn1g1

    integer js2gs, jn1g0, jn1gn            ! FastOpt
    integer north, south                   ! FastOpt
    integer irec                           ! FastOpt
    integer, parameter :: store_kind = 8   ! FastOpt precision of tape

    im2 = im/2

! Set loop limits

    js1g1 = max(1,jfirst-1)
    js2g1 = max(2,jfirst-1)
    js2gcp1 = max(2,jfirst-ng_c-1)   ! NG-1 latitudes on S (starting at 2)
    jn1g1 = min(jm,jlast+1)
    jn2gc = min(jm-1,jlast+ng_c)   ! NG latitudes on N (ending at jm-1)
 
! FastOpt begin
    north = min(2,abs(jord))         ! north == 1 or 2
    south = north-1                  ! south == 0 or 1

    js2gs = max(2,jfirst-south)
    jn1gn = min(jm,jlast+north)
    jn1g0 = min(jm,jlast)

#ifdef  SW_CORE_STORING
!$taf init c_sw_tape1 = static, 1
!$taf init c_sw_tape2 = static, jn2g0-js2g1+1
    c_sw_tape_rec = 0
#endif

#ifdef  TP_CORE_STORING
!$taf init tpcc1_tape = static, jn1gn-js2gs+1
!$taf init tpcc2_tape = static, jn1g0-js2g0+1
   tpcc_tape_rec = 0
#else
   tpcc_tape_rec = c_sw_tape_rec
#endif

! FastOpt end

!
! Treat the special case of ng_c = 1
!
    if ( ng_c == 1 .AND. ng_d > 1 ) THEN
        js2gc1 = js2gc
    else
        js2gc1 = max(2,jfirst-ng_c+1)   ! NG-1 latitudes on S (starting at 2)
    endif

! Get D-grid V-wind at the poles.
    call vpol5(u(1,jfirst), v(1,jfirst), im, jm,            &
               coslon, sinlon, cosl5, sinl5, jfirst, jlast)

    do j=js2gcp1,jn2gc
       do i=1,im-1
          v2(i,j) = v(i,j) + v(i+1,j)
       enddo
          v2(im,j) = v(im,j) + v(1,j)
    enddo

    do j=js2gc,jn2gc
       do i=1,im
          u2(i,j) = u(i,j) + u(i,j+1)
       enddo
    enddo

        if ( jfirst == 1 ) then
! Projection at SP
          us = 0.
          vs = 0.

          do i=1,im2
            us = us + (u2(i+im2,2)-u2(i,2))*sinlon(i)         &
                    + (v2(i,2)-v2(i+im2,2))*coslon(i)
            vs = vs + (u2(i+im2,2)-u2(i,2))*coslon(i)         &
                    + (v2(i+im2,2)-v2(i,2))*sinlon(i)
          enddo

          us = us/im
          vs = vs/im

! SP
          do i=1,im2
            u2(i,1)  = -us*sinlon(i) - vs*coslon(i)
            v2(i,1)  =  us*coslon(i) - vs*sinlon(i)
            u2(i+im2,1)  = -u2(i,1)
            v2(i+im2,1)  = -v2(i,1)
          enddo

          p1ke = 0.125*(u2(1, 1)**2 + v2(1, 1)**2)

        endif

        if ( jlast == jm ) then

! Projection at NP
          un = 0.
          vn = 0.

          j = jm-1
          do i=1,im2
            un = un + (u2(i+im2,j)-u2(i,j))*sinlon(i)        &
                    + (v2(i+im2,j)-v2(i,j))*coslon(i)
            vn = vn + (u2(i,j)-u2(i+im2,j))*coslon(i)        &
                    + (v2(i+im2,j)-v2(i,j))*sinlon(i)
          enddo

          un = un/im
          vn = vn/im

! NP
          do i=1,im2
            u2(i,jm) = -un*sinlon(i) + vn*coslon(i)
            v2(i,jm) = -un*coslon(i) - vn*sinlon(i)
            u2(i+im2,jm) = -u2(i,jm)
            v2(i+im2,jm) = -v2(i,jm)
          enddo

          p2ke = 0.125*(u2(1,jm)**2 + v2(1,jm)**2)

        endif

! A -> C
        do j=js2gc,jn2gc                ! uc needed N*ng S*ng
! i=1
            uc(1,j) = 0.25*(u2(1,j)+u2(im,j))
          do i=2,im
            uc(i,j) = 0.25*(u2(i,j)+u2(i-1,j))
          enddo
        enddo

! FastOpt begin
        do i=1,im
           vc(i,jfirst-2) = 0.
           vc(i,jfirst-1) = 0.
           vc(i,jfirst  ) = 0.
           vc(i,jlast +1) = 0.
           vc(i,jlast +2) = 0.
        enddo
! FastOpt end

        do j=js2gc,jn1gc                ! vc needed N*ng, S*ng (for ycc)
          do i=1,im
            vc(i,j) = 0.25*(v2(i,j)+v2(i,j-1))  ! v2 needed N*ng S*(ng+1)
          enddo
        enddo

        do j=js2g1,jn1g1                   ! cry needed on NS
          do i=1,im
            cry(i,j) = dtdy5*vc(i,j)
          enddo
        enddo

        do j=js2g0,jn1g1                     ! ymass needed on NS
          do i=1,im
            ymass(i,j) = cry(i,j)*cose(j)
          enddo
        enddo

! New va definition
        do j=js2g1,jn2g0                     ! va needed on S (for YCC, iv==1)
          do i=1,im
            va(i,j) = 0.5*(cry(i,j)+cry(i,j+1))
          enddo
        enddo

! SJL: Check if FFSL integer fluxes need to be computed

! FastOpt begin
        do j=jfirst-1,jlast                ! ffsl needed on N*sg S*sg
          do i=1,im
            crx(i,j) = 0.0
          enddo
        enddo
! FastOpt end

        do 2222 j=js2gc,jn2gc                ! ffsl needed on N*sg S*sg
          do i=1,im
            crx(i,j) = uc(i,j)*dtdx2(j)
          enddo
          ffsl(j) = .false.
          if( cosp(j) < zt_c ) then
            do i=1,im
              if( abs(crx(i,j)) > 1. ) then
                ffsl(j) = .true. 
                go to 2222
              endif
            enddo
          endif
2222    continue

! 2D transport of polar filtered delp (for computing fluxes!)
! Update is done on the unfiltered delp

   call tp2c( ptk,  va(1,jfirst),  delpf(1,jfirst-ng_c),    &
              crx(1,jfirst-ng_c), cry(1,jfirst),             &
              im, jm, iord, jord, ng_c, xfx,                 &
              yfx, ffsl, rcap, acosp,                        &
              crx(1,jfirst), ymass, cosp,                    &
              0, jfirst, jlast)

!$taf store xfx,yfx      =  c_sw_tape1, rec=c_sw_tape_rec+1, kind=store_kind

   do j=js2g0,jn2g0                      ! xfx not ghosted
      if( ffsl(j) ) then
         do i=1,im
           xfx(i,j) = xfx(i,j)/sign(max(abs(crx(i,j)),tiny),crx(i,j))
         enddo
      endif
   enddo

! pt-advection using pre-computed mass fluxes
! use tm2 below as the storage for pt increment
! WS 99.09.20 : pt, crx need on N*ng S*ng, yfx on N

    call tp2c(tm2 ,va(1,jfirst), pt(1,jfirst-ng_c),       &
              crx(1,jfirst-ng_c), cry(1,jfirst),          &
              im, jm,  iord, jord, ng_c, fx,              &
              fy(1,jfirst), ffsl, rcap, acosp,            &
              xfx, yfx, cosp, 1, jfirst, jlast)

#if ( !defined ALT_PFT )
! use v2, crx as work arrays

     call pft2d(ptk(1,js2g0), sc(js2g0), dc(1,js2g0), im,   &
                jn2g0-js2g0+1,  ifax, trigs, v2, crx )
     call pft2d(tm2(1,js2g0), sc(js2g0), dc(1,js2g0),  im,  &
                jn2g0-js2g0+1,  ifax, trigs, v2, crx )
#endif

!$taf store tm2(:,:),ptk(:,:) =  c_sw_tape1, rec=c_sw_tape_rec+1, kind=store_kind

    do j=jfirst,jlast
       do i=1,im
          ptk(i,j) = delp(i,j) + ptk(i,j)
          ptc(i,j) = (pt(i,j)*delp(i,j) + tm2(i,j))/ptk(i,j)
       enddo
    enddo

!------------------
! Momentum equation
!------------------

     call ycc(im, jm, fy, vc(1,jfirst-2), va(1,jfirst-1),   &
              va(1,jfirst-1), jord, 1, jfirst, jlast)

     do j=js2g1,jn2g0

          irec = 1+j-js2g1 + c_sw_tape_rec*(jn2g0-js2g1+1)

!$taf store a6,slope = c_sw_tape2, rec=irec, kind=store_kind

          do i=1,im
            cx1(i) = dtdx4(j)*u2(i,j)
          enddo

          sld = .false.
          if( cosp(j) < zt_c ) then
            do i=1,im
              if( abs(cx1(i)) > 1. ) then
                sld = .true. 
                go to 3333
              endif
            enddo
          endif
3333      continue

          p1d(im) = uc(1,j)
          do i=1,im-1
            p1d(i) = uc(i+1,j)
          enddo

          call xtp(im,   sld, fxj, p1d, cx1, iord,   &
                   cx1,  cosp(j),  0,   slope,       &
                   qtmp, al,       ar,  a6)
 
          do i=1,im
            wk1(i,j) = dxdt(j)*fxj(i) + dydt*fy(i,j)
          enddo
     enddo

!$taf store a6 = c_sw_tape1, rec=c_sw_tape_rec+1, kind=store_kind

     if ( jfirst == 1 ) then
          do i=1,im
            wk1(i,1) = p1ke
          enddo
     endif

     if ( jlast == jm ) then
          do i=1,im
            wk1(i,jm) = p2ke
          enddo
     endif

! crx redefined
     do j=js2gc1,jn1gc
            crx(1,j) = dtxe5(j)*u(im,j)
          do i=2,im
            crx(i,j) = dtxe5(j)*u(i-1,j)
          enddo
     enddo

     do j=js1g1,jlast                       ! cry needed on S
          do i=1,im
            cry(i,j) = dtdy5*v(i,j)
          enddo
     enddo

     do j=jfirst,jlast
          do i=1,im
            ymass(i,j) = cry(i,j)*cosp(j)       ! ymass actually unghosted
          enddo
     enddo

     do j=js2g0,jlast
          do i=1,im
            tm2(i,j) = 0.5*(cry(i,j)+cry(i,j-1)) ! cry ghosted on S 
          enddo
     enddo

!    Compute absolute vorticity on the C-grid.

     if ( jfirst == 1 ) then
          do i=1,im
            u2(i,1) = 0.
          enddo
     endif

     do j=js2gc,jn2gc
          do i=1,im
            u2(i,j) = uc(i,j)*cosp(j)
          enddo
     enddo

     if ( jlast == jm ) then
          do i=1,im
            u2(i,jm) = 0.
          enddo
     endif

     do j=js2gc1,jn1gc
! The computed absolute vorticity on C-Grid is assigned to v2
          v2(1,j) = fc(j) + (u2(1,j-1)-u2(1,j))*cye(j) +     &
                    (vc(1,j) - vc(im,j))*rdxe(j)

          do i=2,im
             v2(i,j) = fc(j) + (u2(i,j-1)-u2(i,j))*cye(j) +  &
                       (vc(i,j) - vc(i-1,j))*rdxe(j)
          enddo
     enddo

!$taf store v2,tm2,ymass,crx = c_sw_tape1, rec=c_sw_tape_rec+1, kind=store_kind

     do 2233 j=js2gc1,jn1gc          ! ffsl needed on N*ng S*(ng-1)
          ffsl(j) = .false.
          if( cose(j) < zt_c ) then
            do i=1,im
              if( abs(crx(i,j)) > 1. ) then
                ffsl(j) = .true. 
                go to 2233
              endif
            enddo
          endif
2233    continue

! FastOpt begin
   call tpcc( tm2, ymass, v2(1,jfirst-ng_d), crx(1:im,jfirst-1:jlast+2),  &
              cry(1:im,jfirst:jlast), im, jm, ng_d,                       &
              iord, jord, fx, fy(1:im,jfirst:jlast), ffsl, cose,          &
              jfirst, jlast, slope, qtmp, al, ar, a6 )
! FastOpt end

   do j=js2g0,jn2g0
         uc(1,j) = uc(1,j) + dtdx2(j)*(wk1(im,j)-wk1(1,j)) + dycp(j)*fy(1,j)
      do i=2,im
         uc(i,j) = uc(i,j) + dtdx2(j)*(wk1(i-1,j)-wk1(i,j)) + dycp(j)*fy(i,j)
      enddo
   enddo

   do j=js2g0,jlast
        do i=1,im-1
           vc(i,j) = vc(i,j) + dtdy5*(wk1(i,j-1)-wk1(i,j))-dxe(j)*fx(i+1,j)
        enddo
           vc(im,j) = vc(im,j) + dtdy5*(wk1(im,j-1)-wk1(im,j))-dxe(j)*fx(1,j)
   enddo

 end subroutine c_sw

!--------------------------------------------------------------------------
 subroutine d_sw( u,     v,      uc,       vc,                       &   
                  pt,    delp,   delpf,    cx3,                      &
                  cy3,    mfx,    mfy,     cdx,    cdy,              &
                  dtdx,  dtdxe,  dtxe5,  txe5,  dyce,   rdx,         &
                  cy,    dx,    f0,     js2g0,  jn1g1,               &  
                  im,    jm,     jfirst, jlast, ng_d,   ng_s,        &
                  nq,    iord,   jord,   zt_d,  rcap,   tiny,        &
                  dtdy,  dtdy5,  tdy5,   rdy,   cosp,   acosp,       &
                  cose,  coslon, sinlon, cosl5, sinl5,               &
		  d_sw_tape_rec                                      &
		  )
!--------------------------------------------------------------------------
! Routine for shallow water dynamics on the D-grid

! !USES:

  use precision
  use tp_core

  implicit none

! INPUT:
  integer, intent(in):: im
  integer, intent(in):: jm
  integer, intent(in):: jfirst
  integer, intent(in):: jlast
  integer, intent(in):: iord
  integer, intent(in):: jord
  integer, intent(in):: js2g0
  integer, intent(in):: jn1g1
  integer, intent(in):: ng_d
  integer, intent(in):: ng_s
  integer, intent(in):: nq

  integer :: d_sw_tape_rec                  ! FastOpt

! Prognostic variables:
  real(r8), intent(in):: u(im,jfirst-ng_d:jlast+ng_s)          ! Wind in X
  real(r8), intent(in):: v(im,jfirst-ng_s:jlast+ng_d)          ! Wind in Y
  real(r8), intent(inout):: delp(im,jfirst:jlast)         ! Delta pressure
  real(r8), intent(inout):: pt(im,jfirst-ng_d:jlast+ng_d)! Potential temperature
  real(r8), intent(inout):: delpf(im,jfirst-ng_d:jlast+ng_d)

  real(r8), intent(in):: cosp(jm)
  real(r8), intent(in):: acosp(jm)
  real(r8), intent(in):: cose(jm)

  real(r8), intent(in):: sinlon(im)
  real(r8), intent(in):: coslon(im)
  real(r8), intent(in):: sinl5(im)
  real(r8), intent(in):: cosl5(im)

  real(r8), intent(in):: dtdx(jm)
  real(r8), intent(in):: dtdxe(jm)
  real(r8), intent(in):: dx(jm)
  real(r8), intent(in):: rdx(jm)
  real(r8), intent(in):: cy(jm)
  real(r8), intent(in):: dyce(jm)
  real(r8), intent(in):: dtxe5(jm)
  real(r8), intent(in):: txe5(jm)

  real(r8), intent(in):: cdx(js2g0:jn1g1)
  real(r8), intent(in):: cdy(js2g0:jn1g1)
  real(r8), intent(in):: f0(jfirst-ng_d:jlast+ng_d)

  real(r8), intent(in):: zt_d
  real(r8), intent(in):: rcap
  real(r8), intent(in):: tiny
  real(r8), intent(in):: tdy5
  real(r8), intent(in):: rdy
  real(r8), intent(in):: dtdy
  real(r8), intent(in):: dtdy5

! INPUT/OUTPUT:
  real(r8), intent(inout):: uc(im,jfirst-ng_d:jlast+ng_d)
  real(r8), intent(inout):: vc(im,jfirst-2   :jlast+2 )
  real(r8), intent(inout):: cx3(im,jfirst-ng_d:jlast+ng_d)! Accumulated Courant number in X
  real(r8), intent(inout):: cy3(im,jfirst:jlast+1)        ! Accumulated Courant number in Y
  real(r8), intent(inout):: mfx(im,jfirst:jlast)          ! Mass flux in X  (unghosted)
  real(r8), intent(inout):: mfy(im,jfirst:jlast+1)        ! Mass flux in Y

! Local 
    real(r8)   fx(im,jfirst:jlast)
    real(r8)  xfx(im,jfirst:jlast)

    real(r8)  wk1(im,jfirst-1:jlast+1)
    real(r8)  cry(im,jfirst-1:jlast+1)
    real(r8)   fy(im,jfirst-1:jlast+1)

    real(r8) ymass(im,jfirst: jlast+1) 
    real(r8)   yfx(im,jfirst: jlast+1)

    real(r8)   va(im,jfirst-1:jlast)
    real(r8)   ub(im,jfirst:  jlast+1)

    real(r8)  crx(im,jfirst-ng_d:jlast+ng_d)

    real(r8) fxj(im)
    real(r8) qtmp(-im/3:im+im/3)
    real(r8) slope(-im/3:im+im/3)
    real(r8) al(-im/3:im+im/3)
    real(r8) ar(-im/3:im+im/3)
    real(r8) a6(-im/3:im+im/3)

    real(r8)  c1, c2

    logical ffsl(jm)
    logical sld

    integer i, j
    integer js2gd, jn2g0, jn2g1, jn2gd, jn1gd

    integer, parameter :: store_kind = 8   ! FastOpt precision of tape

! Set loop limits

  jn2g0 = min(jm-1,jlast)
  jn2g1 = min(jm-1,jlast+1)
  js2gd = max(2,jfirst-ng_d)     ! NG latitudes on S (starting at 1)
  jn2gd = min(jm-1,jlast+ng_d)   ! NG latitudes on S (ending at jm-1)
  jn1gd = min(jm,jlast+ng_d)     ! NG latitudes on N (ending at jm)

#define D_SW_TAPING_NO
#ifdef  D_SW_TAPING
!$taf init d_sw_tape1   = static, 1
!$taf init d_sw_tapej   = static, jn1g1-js2g0+1
   d_sw_tape_rec = 0
#endif

! Get C-grid U-wind at poles.

   call upol5(uc(1,jfirst), vc(1,jfirst), im, jm, coslon, sinlon,  &
              cosl5, sinl5, jfirst, jlast)


  do j=js2gd,jn2gd                     ! crx needed on N*ng S*ng
     do i=1,im
        crx(i,j) = dtdx(j)*uc(i,j)
     enddo
  enddo

  do 2225 j=js2gd,jn2gd                ! ffsl needed on N*ng S*ng
     ffsl(j) = .false.
     if( cosp(j) < zt_d ) then
         do i=1,im
            if( abs(crx(i,j)) > 1. ) then
               ffsl(j) = .true. 
               go to 2225
            endif
         enddo
      endif
2225    continue

  do j=js2g0,jn1g1                       ! cry, ymass needed on N
     do i=1,im
        cry(i,j) = dtdy*vc(i,j)
        ymass(i,j) = cry(i,j)*cose(j)
     enddo
  enddo

  do j=js2g0,jn2g0                         ! No ghosting
     do i=1,im
        if( cry(i,j)*cry(i,j+1) > 0. ) then
           if( cry(i,j) > 0. ) then
              va(i,j) = cry(i,j)
           else
              va(i,j) = cry(i,j+1)         ! cry ghosted on N
           endif
        else
           va(i,j) = 0.
        endif
     enddo
  enddo

!$taf store delpf =  d_sw_tape1, key=d_sw_tape_rec+1, kind=store_kind

! transport polar filtered delp
      call tp2c(ub(1,jfirst), va(1,jfirst), delpf(1,jfirst-ng_d),   &
                crx(1,jfirst-ng_d),cry(1,jfirst),im,jm,iord,jord,   &
                ng_d, xfx, yfx, ffsl,                               &
                rcap, acosp,crx(1,jfirst), ymass,                   &
                cosp, 0, jfirst, jlast)

!$taf store xfx,yfx =  d_sw_tape1, key=d_sw_tape_rec+1, kind=store_kind

! <<< Save necessary data for large time step tracer transport >>>
      if( nq > 0 ) then
          do j=js2g0,jn2g0                       ! No ghosting needed
            do i=1,im
              cx3(i,j) = cx3(i,j) + crx(i,j)
              mfx(i,j) = mfx(i,j) + xfx(i,j)
            enddo
          enddo

          do j=js2g0,jlast                      ! No ghosting needed
            do i=1,im
              cy3(i,j) = cy3(i,j) + cry(i,j)
              mfy(i,j) = mfy(i,j) + yfx(i,j)
            enddo
          enddo
      endif

     do j=js2g0,jn2g0                         ! No ghosting needed
        if( ffsl(j) ) then
          do i=1,im
             xfx(i,j) = xfx(i,j)/sign(max(abs(crx(i,j)),tiny),crx(i,j))
          enddo
        endif
     enddo

! Update delp
        do j=jfirst,jlast
          do i=1,im
! SAVE old delp: pressure thickness ~ "air density"
            wk1(i,j) = delp(i,j)
            delp(i,j) = wk1(i,j) + ub(i,j)
          enddo
        enddo

!$taf store wk1, delp =  d_sw_tape1, key=d_sw_tape_rec+1, kind=store_kind

! pt Advection
! FastOpt arguments
  call tp2c(ub(1,jfirst),va(1,jfirst),pt(1,jfirst-ng_d),    &
            crx(1,jfirst-ng_d),cry(1,jfirst),               &
            im,jm,iord,jord,ng_d,fx,fy(1:im,jfirst:jlast+1),    &
            ffsl, rcap, acosp,                              &
            xfx, yfx(1,jfirst), cosp, 1, jfirst,jlast)

!$taf store ub        =  d_sw_tape1, key=d_sw_tape_rec+1, kind=store_kind

! Update pt.
      do j=jfirst,jlast
            do i=1,im
              pt(i,j) = (pt(i,j)*wk1(i,j)+ub(i,j)) / delp(i,j)
            enddo
      enddo

! Compute upwind biased kinetic energy at the four cell corners

! Start using ub as v (CFL) on B-grid (cell corners)
        do j=js2g0,jn1g1                          ! ub needed on N
            ub(1,j) = dtdy5*(vc(1,j) + vc(im,j))  
          do i=2,im
            ub(i,j) = dtdy5*(vc(i,j) + vc(i-1,j))
          enddo
        enddo

! FastOpt arguments
      call ytp(im, jm, fy(1:im,jfirst:jlast+1), v(1,jfirst-ng_d), ub(1,jfirst),  &
               ub(1,jfirst), ng_d, jord, 1, jfirst, jlast)
! End using ub as v (CFL) on B-grid

   do j=js2g0,jn1g1                 ! ub needed on N
       do i=1,im                
          ub(i,j) = dtxe5(j)*(uc(i,j) + uc(i,j-1))
! uc will be used as wrok array after this point
       enddo
   enddo

  do j=js2g0,jn1g1                       ! wk1 needed on N

!$taf store a6,slope = d_sw_tapej, key=1+j-js2g0 + d_sw_tape_rec*(jn1g1-js2g0+1), kind=store_kind

          sld = .false.
          if( cose(j) < zt_d ) then
            do i=1,im
              if( abs(ub(i,j)) > 1. ) then    ! ub ghosted on N
                sld = .true. 
                go to 2235
              endif
            enddo
          endif
2235      continue

     call xtp(im,  sld, fxj, u(1,j), ub(1,j),          &
              iord, ub(1,j), cose(j), 0,               &
              slope, qtmp, al, ar, a6 )

     do i=1,im
        wk1(i,j) =  txe5(j)*fxj(i) + tdy5*fy(i,j)  ! fy ghosted on N
     enddo

  enddo

! Add divergence damping to vector invariant form of the momentum eqn
! (absolute vorticity is damped by ffsl scheme, therefore divergence damping
! provides more consistent dissipation to divergent part of the flow)

!--------------------------
! Perform divergence damping 
!--------------------------

        do j=max(2,jfirst-1), jn2g1                   ! fy need on NS (below)
            do i=1,im
              fy(i,j) = v(i,j)*cosp(j)      ! v ghosted on NS at least
            enddo
        enddo

        do j=js2g0,jn1g1
! i=1
              uc(1,j) = u(im,j) - u(1,j)    ! u ghosted on N at least
            do i=2,im
              uc(i,j) = u(i-1,j) - u(i,j)
            enddo
        enddo

      if ( jfirst == 1 ) then
! j=2
           do i=1,im
              wk1(i,2) = wk1(i,2) - cdy(2)*fy(i, 2) + cdx(2)*uc(i,2)
           enddo
      endif

        do j=max(3,jfirst),jn2g1         ! wk1 needed on N (after TP2D)
            do i=1,im
              wk1(i,j) = wk1(i,j) + cdy(j)*(fy(i,j-1) - fy(i,j))  &
                                  + cdx(j)*uc(i,j)
            enddo
        enddo

      if ( jlast == jm ) then
           do i=1,im
              wk1(i,jm) = wk1(i,jm) + cdy(jm)*fy(i,jm-1) + cdx(jm)*uc(i,jm)
           enddo
      endif
!------------------------------------
! End divergence damping computation
!------------------------------------


! Compute Vorticity on the D grid
! delpf used as work array

      do j=js2gd,jn1gd
          do i=1,im
            delpf(i,j) = u(i,j)*cose(j)   ! u ghosted on N*ng S*ng
          enddo
      enddo


      if ( jfirst ==  1 ) then
          c1 = 0.
          do i=1,im
            c1 = c1 + delpf(i,2)
          end do
          c1 = -c1*rdy*rcap

          do i=1,im
            uc(i,1) = c1
          enddo
      endif

      if ( jlast == jm ) then
          c2 = 0.
          do i=1,im
            c2 = c2 + delpf(i,jm)
          end do
          c2 = c2*rdy*rcap

          do i=1,im
            uc(i,jm) = c2
          enddo
      else

! This is an attempt to avoid ghosting u on N*(ng+1)
          do i=1,im
! DEBUG
!            uc(i,jn2gd) = 0.0
! testing
             uc(i,jn2gd) = 1.e30
          enddo
      endif

      do j=js2gd, min(jm-1,jlast+ng_d-1)
          do i=1,im-1
             uc(i,j) = ( delpf(i,j) - delpf(i,j+1)) * cy(j)  +         &
                        (v(i+1,j) - v(i,j))    * rdx(j)
          enddo
           uc(im,j) = (delpf(im,j) - delpf(im,j+1)) *  cy(j) +         &
                       (v(1,j) - v(im,j)) * rdx(j)
      enddo

! uc is relative vorticity at this point

      do j=max(1,jfirst-ng_d), jn1gd
          do i=1,im
             uc(i,j) = uc(i,j) + f0(j)
! uc is absolute vorticity
          enddo
      enddo

!$taf store uc =  d_sw_tape1, key=d_sw_tape_rec+1, kind=store_kind

! FastOpt arguments
      call tp2d(va(1,jfirst), uc(1,jfirst-ng_d), crx(1,jfirst-ng_d),  &
                cry(1,jfirst), im, jm, iord, jord, ng_d, fx,          &
                fy(1:im,jfirst:jlast+1), ffsl, crx(1,jfirst),         &
                ymass, cosp, 0, jfirst, jlast)

      do j=js2g0,jlast
          do i=1,im-1
            uc(i,j) = dtdxe(j)*(wk1(i,j)-wk1(i+1,j)) + dyce(j)*fy(i,j)
          enddo
           uc(im,j) = dtdxe(j)*(wk1(im,j)-wk1(1,j)) + dyce(j)*fy(im,j)
      enddo

      do j=js2g0,jn2g0
          do i=1,im
            vc(i,j) = dtdy*(wk1(i,j)-wk1(i,j+1)) - dx(j)*fx(i,j)
          enddo
      enddo

 end subroutine d_sw

 end module sw_core
