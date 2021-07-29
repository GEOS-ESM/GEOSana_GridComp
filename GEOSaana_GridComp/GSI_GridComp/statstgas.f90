subroutine statstgas(stats_tgas, ndata)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    statstgas   prints statistics for trace gases
!   prgmmr: weir             org: gmao                date: 2014-04-19
!
! abstract: The routine computes and prints statistics regarding the
!           use of trace gas observations.  Printed information includes
!           that about data counts, quality control decisions,
!           statistics based on the innovations, and penalties - all
!           as a function of observation and satellite type
!
! program history log:
!   2014-04-19  weir     - created file based on statsco and statsoz
!
!   input argument list:
!     stats_tgas - array holding sums from various statistical output
!     ndata(*,1) - number of observations kept for processing
!     ndata(*,2) - number of observations read
!     ndata(*,3) - number of observations kept after read
!
!   output argument list:
!     stats_tgas - array holding sums from various statistical output
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use kinds,     only: r_kind, i_kind
  use constants, only: zero, one
  use obsmod,    only: ndat, iout_tgas, dtype, dsis, dplat
  use tgasinfo,  only: error_tgas, nusis_tgas, nupro, iuse_tgas, jpch_tgas
  use jfunc,     only: jiter
  use qcmod,     only: npres_print, pboto3, ptopo3
  use convinfo,  only: nconvtype, ioctype
  use gridmod,   only: nsig

  implicit none

! Declare passed variables
  real(r_kind),    intent(inout) :: stats_tgas(9,jpch_tgas)
  integer(i_kind), intent(in  )  :: ndata(ndat,3)

! Declare local variables
  real(r_kind),    dimension(ndat) :: rpenal, qcpenal
  integer(i_kind), dimension(ndat) :: icount_asim, iqccount_asim
  logical,         dimension(ndat) :: idisplay

  integer(i_kind) :: j, iasim, i, icerr, ii
  real(r_kind)    :: svar, rsum, stdev, cpen, penalty_all, qcpenalty_all

! Compute and print statistics for trace gas data
  open(iout_tgas, position='append')
  write(iout_tgas, "('OUTER ITERATION jiter = ',i5)") jiter

! Print total penalty over all satellites
  penalty_all   = zero
  qcpenalty_all = zero
  do i = 1,jpch_tgas
     iasim = nint(stats_tgas(7,i))
     if (iasim > 0) then
        penalty_all   = penalty_all   + stats_tgas(5,i)
        qcpenalty_all = qcpenalty_all + stats_tgas(8,i)
     end if
  end do
  write(iout_tgas,*) 'trace gas total   penalty_all = ', penalty_all
  write(iout_tgas,*) 'trace gas total qcpenalty_all = ', qcpenalty_all

! Print counts, bias, rms, stndev as a function of profile number
  icount_asim   = 0
  iqccount_asim = 0
  rpenal   = zero
  qcpenal  = zero
  idisplay = .false.
  do ii = 1,ndat
     do i = 1,jpch_tgas
        iasim = nint(stats_tgas(1,i))
        if (iasim > 0 .and. nusis_tgas(i) == dsis(ii)) then
           if (iuse_tgas(i) == 1) then
              icount_asim(ii)   = icount_asim(ii)   + iasim
              iqccount_asim(ii) = iqccount_asim(ii) + nint(stats_tgas(9,i))

              rpenal(ii)  = rpenal(ii)  + stats_tgas(5,i)
              qcpenal(ii) = qcpenal(ii) + stats_tgas(8,i)
           end if
           idisplay(ii) = .true.
        end if
     end do
  end do

! Write stas to runtime output file
  write(iout_tgas,1101)
  do i = 1,jpch_tgas
     iasim = nint(stats_tgas(1,i))
     if (iasim > 0) then
        svar  = error_tgas(i)
        if (iuse_tgas(i) /= 1) svar = -svar
        rsum  = one/float(iasim)
        icerr = nint(stats_tgas(2,i))
        do j = 3,6   ! j=3: obs-mod(w_biascor)
                     ! j=4: (obs-mod(w_biascor))**2
                     ! j=5: penalty contribution
                     ! j=6: obs
           stats_tgas(j,i) = stats_tgas(j,i)*rsum
        end do
        stats_tgas(4,i) = sqrt(stats_tgas(4,i))
        if (iasim > 1) then
           stdev  = sqrt(stats_tgas(4,i)*stats_tgas(4,i) -                     &
                         stats_tgas(3,i)*stats_tgas(3,i))
        else
           stdev = zero
        end if
        write(iout_tgas,1102) i, nupro(i), nusis_tgas(i), iasim, icerr, svar,  &
                              stats_tgas(6,i), stats_tgas(6,i)-stats_tgas(3,i),&
                              stats_tgas(3,i), stats_tgas(5,i),                &
                              stats_tgas(4,i), stdev
     end if
  end do

! Write obs count to runtime output file
  write(iout_tgas,1109)
  do i = 1,ndat
     if (idisplay(i)) then
        cpen = zero
        if (icount_asim(i) > 0) cpen = rpenal(i)/float(icount_asim(i))
        write(iout_tgas,1115) jiter, dplat(i), dtype(i), ndata(i,2),           &
                              ndata(i,3), icount_asim(i), rpenal(i), cpen,     &
                              qcpenal(i), iqccount_asim(i)
     end if
  end do

1101 format(t13, 'isis', t34, '# used', t41, '# toss', t51, 'scale', t65,      &
            'obs', t77, 'ges', t89, 'omg', t99, 'omg^2', t114, 'Jo', t123,     &
            'stdev')
1102 format(1x, i4, i4, 3x, a20, 2i7, 1x, f8.3, 1x, 6(f11.7,1x))
1109 format(t5, 'it', t13, 'sat', t24, 'inst', t38, '# read', t48, '# kept',   &
            t58, '# used', t80, 'penalty', t92, 'cpen', t114, 'qcpen', t120,   &
            '# qcfail')
1115 format('o-g', 1x, i2.2, 1x, 'tgas ', a10, 1x, a10, 1x, 3(i9,1x),          &
            es22.14, 1x, f8.4, 1x, es22.14, 1x, i8)

! End of trace gas diagnostic print block

  close(iout_tgas)

! End of routine
  return
end subroutine statstgas
