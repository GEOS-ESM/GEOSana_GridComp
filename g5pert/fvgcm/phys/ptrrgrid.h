! Define radiation vertical grid and buffer length for abs/ems
! out-of-core file
 
      integer plevr    ! Number of vertical levels
      integer plevrp   ! plevr + 1

      parameter(plevr = PLEVR, plevrp = plevr + 1)
 
