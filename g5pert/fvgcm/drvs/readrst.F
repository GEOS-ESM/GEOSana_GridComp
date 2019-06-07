      subroutine readrst(iuic,dflnm,job,pick,im,jm,km,nq,jfirst,jlast,
     .                   nstep,nymd,nhms,ps,delp,u,v,pt,q, ng_s, ng_d)

#if defined (SPMD)
      use mod_comm, only : gid, mp_bcst_int, mp_barrier, mp_scatter2d, mp_scatter4d

#define  CPP_PS   pstmp
#define  CPP_U    utmp
#define  CPP_V    vtmp
#define  CPP_DELP delptmp
#define  CPP_PT   pttmp
#define  CPP_Q    qtmp
#else
#define CPP_PS   ps
#define CPP_U    u
#define CPP_V    v
#define CPP_DELP delp
#define CPP_PT   pt
#define CPP_Q    q

#endif

      use m_die, only : MP_die
      use m_ioutil, only : luavail
      use m_strtemplate
      use m_inpak90
      use m_dyn

      implicit none

! !INPUT PARAMETERS:
      integer iuic                         ! Unit number
      character*(*) dflnm                  ! dynamics restart filename
      character*(*) job                    ! job name
      logical, intent(in) :: pick          ! determines whether or not to pick particular date/time
      integer im, jm, km, nq               ! Dimensions of arrays
      integer jfirst, jlast                ! Latitude strip

      integer ng_s, ng_d                   ! info for ghost grids

! !OUTPUT PARAMETERS:
      integer nymd, nhms, nstep            ! Dates 
      real :: ps(im,jfirst:jlast) 
      real :: u(im,jfirst-ng_d:jlast+ng_s,km) 
      real :: v(im,jfirst-ng_s:jlast+ng_d,km) 
      real :: pt(im,jfirst-ng_d:jlast+ng_d,km) 
      real :: delp(im,jfirst:jlast,km)
      real :: q(im,jfirst-ng_d:jlast+ng_d,km,nq)

#if defined (SPMD)
      real pstmp( im,jm ), utmp( im,jm,km ), vtmp( im,jm,km )
      real delptmp( im,jm,km ), pttmp( im,jm,km ), qtmp(im,jm,km,nq)
#endif

      integer i, j, k, n, iret
      logical, save :: start=.true.
      logical exists
      character(len=*),  parameter :: rcname = 'ccmrun.namelist'

      integer iuic2

      character(len=255) :: rstfname, token
      integer :: rc 
      type(dyn_vect) :: w_f     ! fv dynamics vector in eta



!******************************
!*    All                     *
! 25Jan2006 Oloso Altix bug fix: can only open file once
! 13Jun2006 Elena N.  added option to read rst data directly from traj file
! 10Aug2007 Todling   generalized handling of tracers, when using dyn_get
! 13Aug2000 Todling   generalized name of rst file
!******************************

      if ( pick ) then
           call strTemplate ( rstfname, dflnm, 'GRADS',
     &                        xid=trim(job), nymd=nymd, nhms=nhms )
      else
#if defined (BINRST)
           rstfname = 'd_rst'
#else
           inquire( file=rcname, exist=exists )
           if ( exists ) then  
              call i90_loadf (rcname, iret)
              if( iret .eq. 0) then
                 call I90_label('g4_restart_name:', iret)
                 if (iret .eq. 0) then
                     call I90_Gtoken(token,iret)
                     if( iret .eq. 0) then
                         rstfname = trim(token)
                     endif
                 endif
              endif
              call i90_release ()
           endif
#endif
      endif


!******************************
!*    Master                  *
!******************************

#if defined (SPMD)
      call mp_barrier
      if( gid == 0 ) then
#endif

#if defined (BINRST)
        iuic2 = luavail()
        open(iuic2,file=dflnm,form='unformatted',access='sequential')
        read(iuic2) nstep, nymd, nhms
        read(iuic2) CPP_PS,CPP_DELP,CPP_U,CPP_V,CPP_PT
        if(nq .ne. 0) then
         print *,'readrst: nq= ', nq
            do n = 1, nq
              read(iuic2,END=9999) (((CPP_Q(i,j,k,n),i=1,im),j=1,jm),k=1,km)
            enddo
9999        continue
        endif
        close(iuic2)
        write(*,*) 'Done reading/closing dynamic restart file'
#else

        call dyn_null(w_f)
        if ( pick ) then
             call dyn_get ( rstfname, nymd, nhms, w_f, rc, timidx=0, nstep=nstep, vectype=4 )
        else
             call dyn_get ( rstfname, nymd, nhms, w_f, rc, nstep=nstep, vectype=4 )
        endif
          if (rc .ne. 0) 
     .        call MP_die ( 'readrst:', 'could not read dynamics state file to get restarts', rc )

        if ((w_f%grid%im .ne. im) .or. (w_f%grid%jm .ne. jm) .or. (w_f%grid%km .ne. km)) then
           print *, ' Resolution of Dynamic restart vector: ', w_f%grid%im, w_f%grid%jm, w_f%grid%km
           print *, ' Required resolution: ', im, jm, km
           print *, ' Exiting on READRST error'
           call MP_die ( 'readrst:', 'dims()', rc )
        endif

        CPP_PS            = w_f%ps
        CPP_DELP          = w_f%delp
        CPP_U             = w_f%u
        CPP_V             = w_f%v
        CPP_PT            = w_f%pt
        CPP_Q(:,:,:,1)    = w_f%q(:,:,:,1)

        if ( nq .le. w_f%grid%lm ) then
            CPP_Q(:,:,:,2:nq) = w_f%q(:,:,:,2:nq)
        else
            CPP_Q(:,:,:,2:w_f%grid%lm)    = w_f%q(:,:,:,2:w_f%grid%lm)
            CPP_Q(:,:,:,w_f%grid%lm+1:nq) = 0.d0    ! zero out rest of tracers (alt., could read tracer from another file)
        endif

        call dyn_clean ( w_f )

         print *
         print *,'USING RST_HDF TO READ RST FIELDS ON: ', nymd, nhms
         print *

#endif


      where ( abs(CPP_U) > 1e3 ) CPP_U = 0.0  ! care for S.Pole undef from G5GCM

#if defined (SPMD)           
      endif      
      call mp_barrier
      call mp_bcst_int(nstep)
      call mp_bcst_int(nymd)
      call mp_bcst_int(nhms)

! Scatter (broadcast)
      call mp_scatter2d(pstmp, ps, im,  jm, jfirst, jlast, 0) 
      call mp_scatter4d(utmp, u, im,  jm, km, 1, jfirst, jlast, 
     &                  1, km, ng_d, ng_s, 0)
      call mp_scatter4d(vtmp, v, im,  jm, km, 1, jfirst, jlast, 
     &                  1, km, ng_s, ng_d, 0)
      call mp_scatter4d(pttmp, pt, im,  jm, km, 1, jfirst, jlast, 
     &                  1, km, ng_d, ng_d, 0)
      call mp_scatter4d(delptmp, delp, im,  jm, km, 1, jfirst, jlast, 
     &                  1, km, 0, 0, 0)
      call mp_scatter4d(qtmp, q, im,  jm, km, nq, jfirst, jlast, 
     &                  1, km, ng_d, ng_d, 0)
#endif                     
      return
      end
