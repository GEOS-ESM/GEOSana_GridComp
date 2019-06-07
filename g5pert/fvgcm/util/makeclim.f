C----------------------------
C Program 	makeclim
C Description 	reads a grads file and writes out a climate mean grads file
C----------------------------

      program makeclim

      integer im,jm
      parameter ( im = 144 )
      parameter ( jm = 91 )
      parameter ( ny = 6 )                      ! Total number of years
      parameter ( nplev  = 26 )        ! 55L set up
!     parameter ( nplev  = 17 )        ! model top below .1
      parameter ( nupper = 34 )        ! model top below .1
      parameter ( nsurf  = 48 )        ! model top below .1
      parameter ( nrec = nupper*nplev + nsurf)

      real q(im,jm), qc(im,jm), weight(im,jm)
      real undef
      integer nri
      character*80 fname, fnamec
 
! Output file name
      fnamec = 'diag_clim'

! Input
      fname  = 'diag.data'
! For write
      open (12, file=fnamec,access='direct',form='unformatted',
     .      recl=im*jm*4,status='new')

! Read
      open (20, file=fname,access='direct',form='unformatted',
     .          recl=im*jm*4)

      undef = 1.e25

        do n = 1, nrec
          write(6,*) ' nrec ',n

        do m = 1, 12
          do j = 1, jm
          do i = 1, im
          qc(i,j) = 0.
          weight(i,j) = 0.
          enddo
          enddo

        do nn = 1, ny


          nri = n+nrec*(12*(nn-1)+m-1)
          read(20,rec=nri) q

          do j = 1, jm
          do i = 1, im
            if (q(i,j) .ne. undef) then
              qc(i,j) = qc(i,j) + q(i,j)
              weight(i,j) = weight(i,j) + 1.
            endif
          enddo
          enddo

          if (nn.eq.ny) then

            do j = 1, jm
            do i = 1, im
              if (weight(i,j) .ne. 0.) then
                qc(i,j) = qc(i,j) / weight(i,j)
              else
                qc(i,j) = undef
              endif
            enddo
            enddo
            write(12,rec=n+(m-1)*nrec) qc

          endif



        enddo
      enddo
      enddo

        close(12)
        close(20)

      stop
      end

