      module mod_pio  
#if (defined SPMD) &&  (defined PIO )
      use mod_comm, only: mp_barrier, gid
      use mod_comm, only: imr,jnp,nl
      implicit none

      integer im, jm, km, nq
      integer numinfo
      parameter (im=FVGCM_LON, jm=FVGCM_LAT, km=FVGCM_LEV, nq=FVGCM_TRACER)
      parameter (numinfo=6)

      pointer (peinfo_ptr, peinfo)
      real*4 peinfo(numinfo, jm)       !jm is the upper limit 
    
      integer*8 iptr4gout              !defined in allocate_shm_4gout
                                       !Each PE will has its own copy
                                       !it is used to access each PE's
                                       !shared memory. The size for each PE
                                       !is different, and therefore it is not
                                       !easy to specify the shape in advance.

      integer numslices


!..................................................
      contains
!..................................................

      subroutine pio_init()

      integer*8 isize, ipnt, numvar
      character*80 evalue
      integer numpro, id
      integer nl
      integer itr


      call getenv('NUMBER_MLP_PROCESSES',evalue)
      read(evalue,*) numpro

      isize=numinfo*numpro*4
!      write(*,*)'isize for peinfo is ', isize

      id=-1          !must be negative
      itr=1
      call allocate_shm_4gout_c(ipnt, isize, id, itr)  !peinfo
      peinfo_ptr=ipnt


#if defined (CCM)
      call get_numslices(numslices, km)


      call getfd_4gout_c(numpro)           !4diag

      call spawn_iod_4gout_c()             !4diag
#endif

      return

      end subroutine pio_init

!.......................................
      subroutine pio_exit()

      if (gid == 0 ) then
!        
      endif
      return
      end subroutine pio_exit
!.......................................
      subroutine deallocate_shm_4gout()

      call deallocate_shm_4gout_c()


      end subroutine deallocate_shm_4gout

!.......................................
      subroutine allocate_shm_4gout(jfirst, jlast)
      integer   jfirst, jlast
      integer*8 iptr, isize
      integer   itr            !=1, truncate a file

      peinfo(1, gid+1)=im
      peinfo(2, gid+1)=jm
      peinfo(3, gid+1)=numslices
      peinfo(4, gid+1)=jfirst
      peinfo(5, gid+1)=jlast
      peinfo(6, gid+1)=0          !# of records, reset every time?      
      isize = im * numslices * (jlast-jfirst+1) * 4
      itr   = 1
      call allocate_shm_4gout_c(iptr, isize, gid, itr)

      iptr4gout=iptr

      end subroutine allocate_shm_4gout


!.......................................

      subroutine get_numslices(numslices, nl)

      implicit none
      integer numslices, nl
      integer n2d, n3d
      character*128 buf
      integer iutbl, ios
      character*8 fld
      character*3 dim
      character     comment
      integer n2d, n3d

      logical     done
      logical     picked
      logical     v3d

      iutbl = 70
      open (iutbl, file='diag.tbl', form='formatted', status='unknown')
      done = .false.
      comment='!'
      n2d=0
      n3d=0
      numslices=0

      v3d=.false.
      do while (.not. done)
        read (iutbl, '(a)', iostat=ios) buf
!        write(*,*)buf
        if (ios .gt. 0) then
          print *,                        &
             '[diaginit]: Error in reading diagnostic field table.'
          stop
        else if (ios .lt. 0) then    ! EOF
          done = .true.
        else
          call trimleft (buf)
          if (buf(1:1) .ne. comment) then
            read (buf, '(a8,2x,l4)') fld, picked
              if (picked) then
                if ( .not. V3D ) then
                  n2d = n2d + 1
                else
                  n3d = n3d + 1
                endif
              endif
          else
            read (buf, '(a1,1x,a3)') fld, dim
            if (dim .eq. '3-D' ) then
             V3D=.true.
!             write(*,*) '3D'
            endif
          endif
        endif
      enddo
      close (iutbl)

      numslices=n2d + n3d * nl

      end subroutine get_numslices
!***************************
      subroutine collect_io_4gout()
     
      call awake_goutd_c()
!
!      call mp_barrier
!
!      if (gid == 0 ) then
!      call awake_goutd_c()
!      endif

      end subroutine collect_io_4gout

#endif

      end module mod_pio
!*********************************************************
#if (defined SPMD) &&  (defined PIO )

      subroutine collect_io_4gout_f()

!      use mod_pio, only : peinfo    !necessary???, peinfo is shared
      use mod_pio, only : peinfo    !necessary???, peinfo is shared


      implicit none



      integer jfirst, jlast
      character*80 evalue
      integer numpro


      real*4, allocatable :: output(:, :, :)
      integer im, jm
      integer numslices
      integer   pe, petmp
      integer*8 ptr, isize
      integer jfirst, jlast
      integer i,j, n
      integer itr
      integer nrec
      integer jout 


      write(90, *) 'test'
      call getenv('NUMBER_MLP_PROCESSES',evalue)
      read(evalue,*) numpro


      write(90, *) numpro
   
      Do pe=1,  numpro
         if (pe == 1 ) then
             im =nint(peinfo(1, pe))
             jm =nint(peinfo(2, pe))
             numslices=nint(peinfo(3, pe))
             nrec=nint(peinfo(6, pe))
             allocate(output(im, jm, numslices))
         endif
         jfirst=nint(peinfo(4, pe))
         jlast =nint(peinfo(5, pe))
         isize=im*numslices*(jlast - jfirst +1) *4
         petmp=pe-1
         write(90,*)' isize gid pe ', isize, petmp, pe
         itr=0
         call allocate_shm_4gout_c(ptr, isize, petmp, itr)
         call pio_gather_4gout(output, im, jm, numslices, jfirst, jlast, ptr)
!!!         call deallocate_shm()    !will be done later
      enddo

       write(90,*)'after read'
       jout=84

       open (jout, file='diag.bin', form='unformatted',   &
             status='unknown', access='direct',           &
             recl=im*jm*4)

      do n=1, numslices
         write(jout, rec=n+nrec)((output(i,j, n), i=1, im), j=1, jm)
      enddo

      close(jout)
      close(90)

      deallocate(output)

      end subroutine collect_io_4gout_f
!................
      subroutine pio_gather_4gout(output, im, jm, km, jfirst, jlast, ptr)
      integer im, jm, km
      integer jfirst, jlast
      integer*8 ptr
      real*4 output(im, jm, km)
      integer i, j, k
      pointer (in_ptr, in)
      real*4 in(im, km, jlast-jfirst+1)
      

      in_ptr=ptr
    
      do k=1, km
         do j=jfirst, jlast
            do i=1, im
               output(i, j, k) = in(i, k, j-jfirst+1)
            enddo
          enddo
      enddo

      return
      end subroutine pio_gather_4gout

#endif
