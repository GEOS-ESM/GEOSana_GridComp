      subroutine rst_phys(id, iu, pflnm, im, jm, km, oro, snowh, lwup,
     &                    ts,  tssub, mx, sicthk,
     &                    asdir, asdif, aldir, aldif,
     &                    pblht, tpert, qpert,  srfrad,
     &                    fsns,  precst, flwds, precc, precl,
     &                    sols,  soll,  solsd,  solld, qbeg, tke,
     &                    nq, nstep, undef,  touch, jfirst, jlast )
 
! Read/Write restart file for physical paremeterizations
! ID=0: read
! ID=1: write

!-----------
! NCAR flags: 
!-----------
! oro = 0  --> over ocean
! oro = 1  --> over land
! oro = 2  --> over sea ice

#if defined (SPMD )
      use mod_comm, only: gid, mp_bcst_r2d, mp_gath_r2d, 
     .                    mp_scatter4d, mp_gather4d
#define CPP_TKE tketmp
#else
#define CPP_TKE tke
#endif

      implicit none
      integer id              ! flag to determine read/write
      integer iu              ! unit to read/write
      character*(*) pflnm     ! physics restart filename
      integer im              ! e-w dimension
      integer jm              ! n-s diemsnion
      integer km              ! Vertical dimension
      integer, intent(in):: jfirst, jlast 
      integer mx
      integer nq
      integer nstep
      integer i, j, k

      logical  touch         ! parallel initialization
      logical  qbeg          ! Cold start for PBL?
      integer  ios           !indicates if an error occured

      real undef
      real   oro(im,jm)        !CCM surface flag
      real pblht(im,jm)      !pbl height (m?)
      real tpert(im,jm)      !temperature pert.
      real qpert(im,nq,jm)
      real tke(im,jfirst:jlast,km)

#if defined (SPMD )
      real tketmp(im,jm,km)
#endif

      real srfrad(im,jm)
      real  asdir(im,jm)      ! albedo: shortwave, direct
      real  asdif(im,jm)      ! albedo: shortwave, diffuse
      real  aldir(im,jm)      ! albedo: longwave, direct
      real  aldif(im,jm)      ! albedo: longwave, diffuse
      real   lwup(im,jm)      ! longwave up radiative flux
      real     ts(im,jm)      ! sfc temp (merged w/ocean if coupled)
      real  tssub(im,mx,jm)
      real sicthk(im,jm)
      real  snowh(im,jm)      ! ccm snow depth (liquid water)

      real   fsns(im,jm)
      real precst(im,jm)
      real  precc(im,jm)
      real  precl(im,jm)
      real  flwds(im,jm)

      real  sols(im,jm)
      real  soll(im,jm)
      real solsd(im,jm)
      real solld(im,jm)

! Local
      real qt(im,jm)
      integer ic
      real vmax
      real pmax, pmin

#if (!defined SPMD )
      integer gid
      gid = 0
#endif

      if(gid==0) open(iu,file=pflnm,form='unformatted',access='sequential')

      if( nstep == 0 ) then

        if ( gid == 0 ) then
          read(iu) oro
          write(6,*) 'oro read'
          close(iu)
        endif
#if defined ( SPMD )
        call mp_bcst_r2d(im, jm, jfirst, jlast, oro, 0)
#endif

          do j=1, jm
             do i=1, im
                pblht(i,j) = 0.
                tpert(i,j) = 0.
                qpert(i,1,j) = 0.
                snowh(i,j) = 0.
                sicthk(i,j) = 0.
             enddo
          enddo
        return
      endif

      if( id == 0 ) then

!-------------------
! Read restart file 
!-------------------

        if(mx .ne. 4) then
           write(6,*) 'tssub dimension incorrect!  mx=', mx
           stop
        endif

        if( touch) then

! Parallel initialization of surface variables:

!$omp parallel do private (i,j)
 	do j=jfirst,jlast
 	do i=1,im
              oro(i,j) = undef
            pblht(i,j) = undef
	    tpert(i,j) = undef
	  qpert(i,1,j) = undef
	    snowh(i,j) = undef
	     lwup(i,j) = undef
           srfrad(i,j) = undef
             fsns(i,j) = undef
            flwds(i,j) = undef
           precst(i,j) = undef
           precc(i,j) = undef
           precl(i,j) = undef
            sols(i,j) = undef
            soll(i,j) = undef
           solsd(i,j) = undef
           solld(i,j) = undef
              ts(i,j) = undef
           tssub(i,1,j) = undef
           tssub(i,2,j) = undef
           tssub(i,3,j) = undef
           tssub(i,4,j) = undef
            sicthk(i,j) = undef
           asdir(i,j) = undef
           asdif(i,j) = undef
           aldir(i,j) = undef
           aldif(i,j) = undef
 	enddo
 	enddo

        endif

        if ( gid == 0 ) then
            read(iu) oro
            read(iu) pblht
            read(iu) tpert
            read(iu) qt
            read(iu) snowh
            read(iu) lwup
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
            read(iu) ts
            read(iu) tssub
            read(iu) sicthk
            read(iu) asdir
            read(iu) asdif
            read(iu) aldir
            read(iu) aldif
#if defined (ALT_PBL)
            read(iu,err=999,iostat=ios) CPP_TKE
 999        continue
            if (ios.ne.0) then
               write(6,*) 'Cold starting TKE (2.5 level PBL)'
               do k = 1,km
               do j =jfirst, jlast
               do i = 1,im
                  CPP_TKE(i,j,k) = 0.0
               enddo
               enddo
               enddo
               qbeg = .true.
            else
               qbeg = .false.
               write(6,*) 'Re-starting TKE (2.5 level PBL)'
            endif

#endif
            endif              ! end gid == 0 section

#if ( !defined SPMD )

!$omp  parallel do private (i,j)
            do j=1,jm
               do i=1,im
                  qpert(i,1,j) = qt(i,j)
               enddo
            enddo
#else
            call mp_bcst_r2d(im, jm, jfirst, jlast, oro, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, pblht, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, tpert, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, qt, 0)

!$omp  parallel do private (i,j)
            do j=jfirst,jlast
               do i=1,im
                  qpert(i,1,j) = qt(i,j)
               enddo
            enddo

           if ( nq > 1 ) then
           do ic=2,nq
!$omp parallel do private (i,j)
             do j=jfirst,jlast
               do i=1,im
                  qpert(i,ic,j) = undef
               enddo
             enddo
           enddo
           endif

            call mp_bcst_r2d(im, jm, jfirst, jlast, snowh, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, lwup, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, srfrad, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, fsns, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, precst, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, flwds, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, precc, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, precl, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, sols, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, soll, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, solsd, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, solld, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, ts, 0)

            do k=1,mx

            if ( gid == 0 ) then
               do j=1,jm
                  do i=1,im
                     qt(i,j) = tssub(i,k,j)
                  enddo
               enddo
            endif

            call mp_bcst_r2d(im, jm, jfirst, jlast, qt, 0)

            if ( gid /= 0 ) then
!$omp parallel do private (i,j)
               do j=jfirst,jlast
                  do i=1,im
                     tssub(i,k,j) = qt(i,j)
                  enddo
               enddo
            endif

            enddo

            call mp_bcst_r2d(im, jm, jfirst, jlast, sicthk, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, asdir, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, asdif, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, aldir, 0)
            call mp_bcst_r2d(im, jm, jfirst, jlast, aldif, 0)
#if defined (ALT_PBL)
            call mp_scatter4d(tketmp, tke, im, jm, km, 1, jfirst, jlast,
     &                        1, km, 0, 0, 0)
#endif
#endif

      else

!-------------------
! Write restart file 
!-------------------

#if defined (SPMD )

            do k=1,mx
!$omp parallel do private (i,j)
               do j=jfirst,jlast
                  do i=1,im
                     qt(i,j) = tssub(i,k,j)
                  enddo
               enddo
               call mp_gath_r2d(im, jm, jfirst, jlast, qt, 0)
               if ( gid == 0 ) then
               do j=1,jm
                  do i=1,im
                     tssub(i,k,j) = qt(i,j)
                  enddo
               enddo
               endif
            enddo

!$omp parallel do private (i,j)
            do j=jfirst,jlast
               do i=1,im
                  qt(i,j) = qpert(i,1,j)
               enddo
            enddo

            call mp_gath_r2d(im, jm, jfirst, jlast, oro, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, pblht, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, tpert, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, qt, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, snowh, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, lwup, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, srfrad, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, fsns, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, precst, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, flwds, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, precc, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, precl, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, sols, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, soll, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, solsd, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, solld, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, ts, 0)

            call mp_gath_r2d(im, jm, jfirst, jlast, sicthk, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, asdir, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, asdif, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, aldir, 0)
            call mp_gath_r2d(im, jm, jfirst, jlast, aldif, 0)
#if defined (ALT_PBL)
            call mp_gather4d(tke, tketmp, im, jm, km, 1, jfirst, jlast,
     .                       1, km, 0, 0, 0)
#endif
#else

!$omp  parallel do private (i,j)
            do j=1,jm
               do i=1,im
                  qt(i,j) = qpert(i,1,j)
               enddo
            enddo
#endif

          if ( gid == 0 ) then
            write(iu) oro
            write(iu) pblht
            write(iu) tpert
            write(iu) qt
            write(iu) snowh
            write(iu) lwup
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
            write(iu) ts
            write(iu) tssub
            write(iu) sicthk
            write(iu) asdir
            write(iu) asdif
            write(iu) aldir
            write(iu) aldif
#if defined (ALT_PBL)
            write(iu) tke
#endif
          endif
      endif

      if ( gid == 0 ) close(iu)

! Check Max/Min of pblht
      pmax = vmax(pblht, pmin, im, jm, jfirst, jlast)
      if( gid == 0 ) write(6,*) 'max pbl (m)=', pmax, ' min =', pmin

      pmax = vmax(ts, pmin, im, jm, jfirst, jlast)
      if( gid == 0 ) write(6,*) 'max ts=', pmax, ' min =', pmin

      return
      end
