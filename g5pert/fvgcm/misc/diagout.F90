      subroutine diagout (jout,   qname,  qvdim,  cdiag, fldloc,   pick,    &
                          tavg,   pstart, pend,   nrec,  diagattr, diagbuf, &
                          jfirst, jlast,  nslice, oflnm, undef,    date,    &
                          time,   idt,    vcoord, unit,  label)

#if defined (SPMD )
      use mod_comm
#if defined ( PIO )
      use mod_pio, only : peinfo,  iptr4gout, collect_io_4gout,             &
                          allocate_shm_4gout, deallocate_shm_4gout
#endif
#endif

      use precision
#if defined (GFIO)
      use m_die, only: die
#endif
      implicit      none

#include <diag.h>

      integer  nslice
      integer  jfirst, jlast
      character*8  qname(pdiag)                ! prescribed field table
      real(r4) diagbuf(imr,nslice,jfirst:jlast) ! diagnostic field buffer
      integer  cdiag(pdiag)                ! counter for time averaging
      integer  qvdim(pdiag)                ! vertical dimension of field
      integer  fldloc(pdiag)               ! location of field in diagbuf
      integer  pick(pdiag)                 ! flag for selected fields (1=true)
      integer  tavg(pdiag)                 ! flag for time-averaged fields (1=true)
      integer  diagattr(5,pdiag)           ! field attributes
      integer  nrec                        ! diagnostics output record number
      integer  pstart                      ! diagnostic to start with (either 1 or n3ddiag)
      integer  pend                        ! diagnostic to end with (either n2ddiag or 
                                           ! n2ddiag+n3ddiag)

      character*(*)      oflnm                 ! output file name for HDF output
      character*(*)      unit(*)               ! unit for output variables
      character*(*)      label(*)              ! long name for output variables
      integer       date, time, idt, inc       ! yyyymmdd, nhms, nsecs, time increment(hhmmss)
      real          undef                      ! undef values for HDF output
      real          vcoord(nl)                 ! vertical levels
      real          valid_range(2, pdiag)      ! range values
      real          packing_range(2, pdiag)

      integer jout
      integer n
      integer i
      integer j
      integer k, kk
      integer kkfirst, kklast
      real(r8) ravg
      real(r4) precc(imr,jfirst:jlast)
      real(r4) precl(imr,jfirst:jlast)
      real     wk3(imr,jnp,nl), wk2(imr,jnp)       ! Global 2D work array
      real(r4) wk2_bin(imr,jnp)       ! Global 2D work array
      real(r4) pmx, pmin, pmean
      real(r4) vmax4, gmean4
      real ::  undef_check = 1.e15    ! undefined lower limit for one bumping with an undef

!   Yin 01.10.31:  Add an option to write HDF format output directly.
!     If the output file name is 'diag.bin' or 'diag.sfc.bin', the
!     output will be in GrADS binary format. Otherwise, it is HDF
! HDF related variables
      character*128   :: source = 'Data Assimilation Office, NASA/GSFC'
      character*128   ::   contact = 'data@dao.gsfc.nasa.gov'
      character*128   ::   oTitle = 'FVGCM Dynamics State Vector'
      character*128   ::   myname = 'diagout.F'
      integer         ::   out_fmode = 0        ! 0 for READ-WRITE
      integer         ::   out_fid              ! output file ID
      integer         ::   rc = 0               ! return error code
      integer         ::   tVars = 0            ! total variables to write out
      integer         ::   outPrec = 0          ! Output file precision:
                                                ! 0 = 32 bits,  1 = 64bits
      character(len=10), pointer  ::  oName(:)  ! output variable names
      integer, pointer            ::  oVdim(:)  ! number of vertical levels
      character(len=28), pointer  ::  oUnit(:)  ! output variable unit
      character(len=128), pointer ::  lName(:)  ! long name

      real lon(imr), lat(jnp)

      integer nrec0
#if defined ( SPMD )
#if defined ( PIO )
      pointer (ptr_diag, diag_buf)
      real*4  diag_buf(imr,nslice,jlast-jfirst+1)
#if defined (NO_PIO_BUFFER)
      call allocate_shm_4gout(jfirst, jlast)
#endif
      ptr_diag=iptr4gout
#endif
#else 
      integer gid
      gid = 0
#endif

      inc  = 10000 * (idt/60/60) + 100 * mod(idt/60, 60) + mod(idt, 60)

      precc(1,jfirst) = -999.
      precl(1,jfirst) = -999.

      do j = 1, jnp
         lat(j) = -90 + 180. / real(jnp - 1) * (j-1)
      enddo

      do i = 1, imr
         lon(i) = 360. / real(imr) * (i-1)
      enddo

      tVars = 0
      do n = pstart, pend
        if (pick(n) .eq. 1) then
           tVars = tVars + 1
        end if
      end do
      allocate ( oUnit(tVars), oVdim(tVars), oName(tVars), lName(tVars),    &
                 stat = rc )
!     if ( rc /= 0 )  call die (myname, 'can not allocate oUnit')

      tVars = 0
      do n = pstart, pend
        if (pick(n) .eq. 1) then
           tVars = tVars + 1
           oName(tVars) = trim(qname(n))
           oUnit(tVars) = unit(n)
           lName(tVars) = label(n)
           valid_range(1,tVars) = undef
           valid_range(2,tVars) = undef
           packing_range(1,tVars) = undef
           packing_range(2,tVars) = undef
           if (qvdim(n) .le. 1) then
              oVdim(tVars) = 0
           else
              oVdim(tVars) = nl
           end if
         end if
      end do

      do n = pstart, pend
        pick(n)   = diagattr(1,n)
        tavg(n)   = diagattr(2,n)
        qvdim(n)  = diagattr(3,n)
        fldloc(n) = diagattr(4,n)
        cdiag(n)  = diagattr(5,n)
      end do

!$omp parallel do                                           &
!$omp default(shared)                                       &
!$omp private(i,j,k,kk,n,ravg)

      do j = jfirst, jlast
        do n = pstart, pend
          if (pick(n) == 1 .and. tavg(n) == 1) then
            ravg = 1. / float( cdiag(n) )
            do k = 1, qvdim(n)
              kk = fldloc(n) + k - 1
              do i = 1, imr
                if (diagbuf(i,kk,j) .lt. undef_check) then
                  diagbuf(i,kk,j) = diagbuf(i,kk,j) * ravg
                else
                  diagbuf(i,kk,j) = undef
                endif
              enddo

              if(qname(n) == 'PRECC') then
                 do i=1, imr
                    precc(i,j) = diagbuf(i,kk,j)
                 enddo
              elseif(qname(n) == 'PRECL') then
                 do i=1, imr
                    precl(i,j) = diagbuf(i,kk,j)
                 enddo
              endif
            enddo
          endif
        enddo
      enddo

#if defined (GFIO)
!     If the output file names are not diag.bin or diag.sfc.bin, try to open
!     output file oflnm first. If oflnm doesn't exist, create it.
      if ( index(oflnm, 'diag.bin') .le. 0 .and.                              &
           index(oflnm, 'diag.sfc.bin') .le. 0 )  then
         call GFIO_Open(oflnm, out_fmode, out_fid, rc)
         if ( rc /= 0 )  then
            call GFIO_Create ( oflnm, oTitle, source, contact, undef,         &
                        imr, jnp, oVdim(1), lon, lat, vcoord, "hPa",          &
                          date, time, inc, tVars, oName, lName, oUnit,        &
                          oVdim, valid_range, packing_range, outPrec,         &
                          out_fid, rc )
           if (rc /= 0)  call die (myname, 'wrong in GFIO_Create or GFIO_Open')
         endif
       end if
#endif

!
! ... Write to file
!     HDF with parallel IO not supported


      if ( index(oflnm, 'diag.bin') .gt. 0                                   &
             .or. index(oflnm, 'diag.sfc.bin') .gt. 0 )  then
      nrec0 = 0                                  !for PIO
      do n=pstart,pend
      if (pick(n) == 1) then                     !in case of kk=0~qvdim(n)-1 
         kkfirst=fldloc(n)
         kklast =fldloc(n)+qvdim(n)-1
#if defined (PIO)
#if defined (CCM) && defined (SPMD)
         do j = jfirst, jlast
            k=j-jfirst+1
!$omp parallel do private(i,kk)
            do kk= kkfirst, kklast
               do i=1,imr
                  diag_buf(i, kk, k) = diagbuf(i,kk,j)
               enddo
            enddo
         enddo
         nrec0 = nrec0 + qvdim(n)
#endif
#else
         do kk= kkfirst, kklast                  !bw added
!$omp parallel do private(i,j)       !will be moved outside kk loop later
            do j = jfirst, jlast
               do i=1,imr
                  wk2_bin(i,j) = diagbuf(i,kk,j)
               enddo
            enddo
#if defined (SPMD )
            call mp_gath4_r2d(imr, jnp, jfirst, jlast, wk2_bin, 0)
#endif
            if(gid==0) write(jout,rec=nrec) wk2_bin
            nrec = nrec + 1
         enddo
#endif 
      endif
      enddo

#if defined (PIO)  && defined (SPMD)
      peinfo(3, gid+1) = nrec0          !# of records which will be adeed
      peinfo(6, gid+1) = nrec -1        !# of the current records
      call mp_barrier
      if ( gid == 0 ) then
         call collect_io_4gout()
      endif
      nrec = nrec + nrec0
#endif


#if defined (GFIO)
      else
        do n=pstart,pend
        if (pick(n) == 1) then                     !in case of kk=0~qvdim(n)-1
           kkfirst=fldloc(n)
           kklast =fldloc(n)+qvdim(n)-1
           k = 0
           do kk= kkfirst, kklast
           k = k + 1
!$omp parallel do private(i,j)
              do j = jfirst, jlast
              do i = 1, imr
                  wk2(i,j) = diagbuf(i,kk,j)
              enddo
              enddo
#if defined (SPMD )
             call mp_gath_r2d(imr, jnp, jfirst, jlast, wk2, 0)
#endif
             if(gid==0) then
               if (qvdim(n) .gt. 1) then
                 do j = 1, jnp
                 do i = 1, imr
                     wk3(i,j,nl-k+1) = wk2(i,j)
                 enddo
                 enddo
               endif
             endif
           enddo

           if ( qvdim(n) .le. 1 ) then
              call GFIO_PutVar (out_fid,trim(qname(n)),date,time,         &
                                imr, jnp, 0, 1, wk2, rc )
              if ( rc /= 0 )  call die (myname, 'wrong in GFIO_PutVar')
           else
              call GFIO_PutVar (out_fid,trim(qname(n)),date,time,         &
                                imr, jnp, 1, oVdim(1), wk3, rc )
              if ( rc /= 0 )  call die (myname, 'wrong in GFIO_PutVar')
           endif

        endif
        enddo
#endif
      endif

#if defined (GFIO)
      if ( index(oflnm, 'diag.bin') .le. 0 .and.                          &
           index(oflnm, 'diag.sfc.bin') .le. 0 ) then
         call GFIO_Close ( out_fid, rc )
      end if
#endif

      if (precc(1,jfirst) .ne. -999.) then
      pmean = gmean4(imr, jnp, jfirst, jlast, precc(1,jfirst))
      pmx = vmax4(precc, pmin, imr, jnp, jfirst, jlast)
      if(gid==0) write(6,*) ' PRECCmax=', pmx,' min=',pmin, ' Mean=', pmean
      endif

      if (precl(1,jfirst) .ne. -999.) then
      pmean = gmean4(imr, jnp, jfirst, jlast, precl(1,jfirst))
      Pmx = vmax4(precl, pmin, imr, jnp, jfirst, jlast)
      if(gid==0) write(6,*) ' PRECLmax=', pmx,' min=',pmin, ' Mean=', pmean
      endif

      do n = pstart, pend
        cdiag(n) = 0
        diagattr(1,n) = pick(n)
        diagattr(2,n) = tavg(n)
        diagattr(3,n) = qvdim(n)
        diagattr(4,n) = fldloc(n)
        diagattr(5,n) = cdiag(n)
      end do

! ... Reset history buffer and counter

!$omp  parallel do                                                  &
!$omp  default(shared)                                              &
!$omp  private(i,j,k,kk,n)

      do j = jfirst, jlast
        do n=pstart,pend
!!!        do k = 1, qvdim(n)
!!!           kk = fldloc(n) + k - 1
        if (pick(n) == 1) then
          do kk = fldloc(n), qvdim(n)+fldloc(n)-1

             do i = 1, imr
              diagbuf(i,kk,j) = 0.
             enddo

           enddo
         endif
         enddo
       enddo

      deallocate( oUnit, oVdim, oName, lName )

#if defined (NO_PIO_BUFFER)
      call deallocate_shm_4gout()
#endif

      return
      end
