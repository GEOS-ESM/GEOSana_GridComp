      subroutine outfld (diagattr ,field  ,idim  ,lat  ,diagbuf)

      use precision
      implicit none

#include <params.h>
#include <pmgrid.h>

! Input arguments
 
      integer idim         ! Longitude dimension of field array
      integer pick         ! flag for selected fields (1=true)
      integer tavg         ! flag for time-averaged fields (1=true)
      integer qvdim        ! vertical dimension of field
      integer fldloc       ! location of field in diagbuf
      integer cdiag        ! counter for time averaging
      integer lat          ! Latitude index
      integer i
      integer k
      integer kk

      integer diagattr(5)

      real(r4) diagbuf(idim,*)
      real(r8) field(idim,*)   ! Array containing field values

      pick   = diagattr(1)
      tavg   = diagattr(2)
      qvdim  = diagattr(3)
      fldloc = diagattr(4)
      cdiag  = diagattr(5)

      if (pick .eq. 1) then

        if (tavg .eq. 1) then
!
! ... Time-averaged fields
!
          do k = 1, qvdim
! ... Top-down
!           kk = fldloc + k - 1
! ... Bottom-up
            kk = fldloc + (qvdim - k)
            do i = 1, idim
              diagbuf(i,kk) = diagbuf(i,kk) + field(i,k)
            enddo
          enddo

          if ( lat == begj ) then
            cdiag = cdiag + 1
            diagattr(5) = cdiag
          endif
        else
!
! ... Instantaneous fields
!
          do k = 1, qvdim
! ... Top-down
!           kk = fldloc + k - 1
! ... Bottom-up
            kk = fldloc + (qvdim - k)
            do i = 1, idim
              diagbuf(i,kk) = field(i,k)
            enddo
          enddo

          if ( lat == begj ) then 
            cdiag = 1
            diagattr(5) = cdiag
          endif
        endif

      endif
 
      return
      end
