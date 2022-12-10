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
  use constants, only: zero, one, varlen => max_varname_length
  use obsmod,    only: ndat, iout_tgas, dtype, dsis, dplat
  use tgasinfo,  only: jpch_tgas, nusis_tgas, ichanl_tgas, iuse_tgas,          &
                       error_tgas, vname_tgas, vunit_tgas
  use jfunc,     only: jiter
  use qcmod,     only: npres_print, pboto3, ptopo3
  use convinfo,  only: nconvtype, ioctype
  use gridmod,   only: nsig

  implicit none

! Declare passed variables
  real(r_kind),    intent(inout) :: stats_tgas(15,jpch_tgas)
  integer(i_kind), intent(in  )  :: ndata(ndat,3)

! Declare local variables
  real(r_kind),    dimension(ndat) :: rpenal, qcpenal
  integer(i_kind), dimension(ndat) :: icount_asim, iqccount_asim
  logical,         dimension(ndat) :: idisplay

  integer(i_kind) :: ii, jj, nn, iasim, icerr, numnh, numtr, numsh
  real(r_kind)    :: sclnat, svar, rsum, stdev, cpen
  real(r_kind)    :: avgnh, avgtr, avgsh, pennh, pentr, pensh
  real(r_kind)    :: penalty_all, qcpenalty_all

! Compute and print statistics for trace gas data
  open(iout_tgas,position='append')
  write(iout_tgas,*) 'statstgas: OUTER ITERATION, jiter =  ', jiter

! Hack to output stats in "natural" units
  do jj = 1,jpch_tgas
     sclnat = 1.e+00_r_kind

!    Leave log-space obs as is
     if (trim(vunit_tgas(jj)) == 'ln')  cycle
     if (trim(vunit_tgas(jj)) == 'log') cycle
     if (trim(vunit_tgas(jj)) == 'l10') cycle

!    Convert mol/m2 to mmol/m2
     if (trim(vunit_tgas(jj)) == 'molm2') then
        sclnat = 1.e+03_r_kind
     else
        if (trim(vname_tgas(jj)) == 'co2')   sclnat = 1.e+06_r_kind
        if (trim(vname_tgas(jj)) == 'co')    sclnat = 1.e+09_r_kind
        if (trim(vname_tgas(jj)) == 'ch4')   sclnat = 1.e+06_r_kind
        if (trim(vname_tgas(jj)) == 'sf6')   sclnat = 1.e+12_r_kind
        if (trim(vname_tgas(jj)) == 'n2o')   sclnat = 1.e+06_r_kind
        if (trim(vname_tgas(jj)) == 'h2o')   sclnat = 1.e+06_r_kind ! stratospheric scale
        if (trim(vname_tgas(jj)) == 'hno3')  sclnat = 1.e+09_r_kind
        if (trim(vname_tgas(jj)) == 'hcl')   sclnat = 1.e+09_r_kind
        if (trim(vname_tgas(jj)) == 'clo')   sclnat = 1.e+09_r_kind
     end if

     stats_tgas( 3,jj) = stats_tgas( 3,jj) * sclnat
     stats_tgas( 4,jj) = stats_tgas( 4,jj) * sclnat**2
     stats_tgas( 6,jj) = stats_tgas( 6,jj) * sclnat

     stats_tgas(11,jj) = stats_tgas(11,jj) * sclnat
     stats_tgas(14,jj) = stats_tgas(14,jj) * sclnat
  end do

! Print total penalty over all satellites, including diag only lines
  penalty_all   = zero
  qcpenalty_all = zero
  do jj = 1,jpch_tgas
     if (nint(stats_tgas(7,jj)) > 0) then
        penalty_all   = penalty_all   + stats_tgas(5,jj)
        qcpenalty_all = qcpenalty_all + stats_tgas(8,jj)
     end if
  end do
  write(iout_tgas,*) 'trace gas total   penalty_all = ',   penalty_all
  write(iout_tgas,*) 'trace gas total qcpenalty_all = ', qcpenalty_all

! Print counts, bias, rms, stndev as a function of profile number
  icount_asim   = 0
  iqccount_asim = 0
  rpenal   = zero
  qcpenal  = zero
  idisplay = .false.
  do ii = 1,ndat
     do jj = 1,jpch_tgas
        iasim = nint(stats_tgas(1,jj))
        if (iasim > 0 .and. nusis_tgas(jj) == dsis(ii)) then
           if (iuse_tgas(jj) == 1) then
              icount_asim(ii)   = icount_asim(ii)   + iasim
              iqccount_asim(ii) = iqccount_asim(ii) + nint(stats_tgas(9,jj))
           end if

           rpenal(ii)   = rpenal(ii)  + stats_tgas(5,jj)
           qcpenal(ii)  = qcpenal(ii) + stats_tgas(8,jj)
           idisplay(ii) = .true.
        end if
     end do
  end do

! Write stats to runtime output file
  write(iout_tgas,1101)
  do jj = 1,jpch_tgas
     iasim = nint(stats_tgas(1,jj))
     icerr = nint(stats_tgas(2,jj))
     if (iasim > 0) then
        svar = error_tgas(jj)
        if (iuse_tgas(jj) /= 1) svar = -svar

        rsum = one/float(iasim)
        do nn = 3,6  ! nn=3: obs-mod(w_biascor)
                     ! nn=4: (obs-mod(w_biascor))**2
                     ! nn=5: penalty contribution
                     ! nn=6: obs
           stats_tgas(nn,jj) = stats_tgas(nn,jj)*rsum
        end do

!       Recall multiplication by rsum above
        if (iasim > 1) then
           stdev = sqrt(stats_tgas(4,jj) - stats_tgas(3,jj)**2)
        else
           stdev = zero
        end if
        stats_tgas(4,jj) = sqrt(stats_tgas(4,jj))

!       Global totals
        write(iout_tgas,1102) 'Gl', iuse_tgas(jj), ichanl_tgas(jj),            &
                              nusis_tgas(jj), iasim, icerr, svar,              &
                              stats_tgas(6,jj),                                &
                              stats_tgas(6,jj) - stats_tgas(3,jj),             &
                              stats_tgas(3,jj), stats_tgas(5,jj),              &
                              stats_tgas(4,jj), stdev
     end if
  end do

! Display zonal bands if there's anything to show (only diagsave steps)
  if (maxval(stats_tgas(10,:) + stats_tgas(13,:)) > 0) then
     write(iout_tgas,1103)
     write(iout_tgas,1104)
     do jj = 1,jpch_tgas
        iasim = nint(stats_tgas(1,jj))
        icerr = nint(stats_tgas(2,jj))
        if (iasim > 0) then
           rsum = one/float(iasim)

!          Northern Extratropical totals
           numnh = nint(stats_tgas(10,jj))
           if (numnh > 0) then
              avgnh = stats_tgas(11,jj)/float(numnh)
              pennh = stats_tgas(12,jj)/float(numnh)
           else
              avgnh = 0.
              pennh = 0.
           end if

!          Tropical totals
           numtr = nint(stats_tgas(1,jj) - stats_tgas(10,jj)                   &
                                         - stats_tgas(13,jj))
           if (numtr > 0) then
              avgtr = (stats_tgas(3,jj)/rsum - stats_tgas(11,jj)               &
                                             - stats_tgas(14,jj))/float(numtr)
              pentr = (stats_tgas(5,jj)/rsum - stats_tgas(12,jj)               &
                                             - stats_tgas(15,jj))/float(numtr)
           else
              avgtr = 0.
              pentr = 0.
           end if

!          Southern Extratropical totals
           numsh = nint(stats_tgas(13,jj))
           if (numsh > 0) then
              avgsh = stats_tgas(14,jj)/float(numsh)
              pensh = stats_tgas(15,jj)/float(numsh)
           else
              avgsh = 0.
              pensh = 0.
           end if

           write(iout_tgas,1105) 'Zo', iuse_tgas(jj), ichanl_tgas(jj),         &
                                 nusis_tgas(jj), numnh, avgnh, pennh,          &
                                 numtr, avgtr, pentr, numsh, avgsh, pensh
        end if
     end do
  end if

! Write obs count to runtime output file
  write(iout_tgas,1109)
  do ii = 1,ndat
     if (idisplay(ii)) then
        cpen = zero
        if (icount_asim(ii) > 0) cpen = rpenal(ii)/float(icount_asim(ii))
        write(iout_tgas,1115) jiter, dsis(ii), ndata(ii,2), ndata(ii,3),       &
                              icount_asim(ii), rpenal(ii), cpen,               &
                              qcpenal(ii), iqccount_asim(ii)
     end if
  end do

1101 format(t10, 'kc', t13, 'isis', t34, '# used', t41, '# toss',              &
            t49, 'inflate', t63, '<obs>', t75, '<ges>', t87, '<omg>',          &
            t95, '<omg^2/R>', t108, 'rms(omg)', t120, 'std(omg)')
1102 format(1x, a2, 1x, i2, 1x, i4, 1x, a20, 2i7, 1x, f8.3, 1x, 6(f11.6,1x))


1103 format(t34, 'Northern Extratropics', t66, 'Tropics',                      &
            t98, 'Southern Extratropics')
1104 format(t34, '# used', t47, '<omg>', t55, '<omg^2/R>', t66, '# used',      &
            t79, '<omg>', t87, '<omg^2/R>', t98, '# used', t111, '<omg>',      &
            t119, '<omg^2/R>')
1105 format(1x, a2, 1x, i2, 1x, i4, 1x, a20, 3(i7,1x,f11.6,1x,f11.6,1x))

1109 format(t5, 'it', t13, 'isis', t38, '# read', t48, '# kept',               &
            t57, '# assim', t80, 'penalty', t92, 'cpen', t110, 'qcpenalty',    &
            t120, '# qcfail')
1115 format('o-g', 1x, i2.2, 1x, 'tgas ', a20, 2x, 3(i9,1x), es22.14, 1x,      &
            f8.4, 1x, es22.14, 1x, i8)

! End of trace gas diagnostic print block
  close(iout_tgas)

! End of routine
  return
end subroutine statstgas
