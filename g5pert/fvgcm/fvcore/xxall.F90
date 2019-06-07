!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS        !
!-------------------------------------------------------------------------
!
! !DESCRIPTION: dummy routines for user defined tape
!
! \begin{verbatim}
!   Copyright (C) 2002
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
! \end{verbatim}
!
! !REVISION HISTORY:
!
!  10Jan2003  FastOpt   Creation this file.
!
!-------------------------------------------------------------------------
subroutine xxread( str, len, k, l, z, prec, count, rec )
  implicit none
  character*(*) str
  integer :: len
  integer :: k, l
  real    :: z
  integer :: prec
  integer :: count
  integer :: rec
end subroutine xxread


subroutine xxwrite( str, len, k, l, z, prec, count, rec )
  implicit none
  character*(*) str
  integer :: len
  integer :: k, l
  real    :: z
  integer :: prec
  integer :: count
  integer :: rec
end subroutine xxwrite


subroutine xxopen( str, len, k, l, m, n )
  implicit none
  character*(*) str
  integer :: len
  integer :: k, l, m, n
end subroutine xxopen


subroutine xxclose( str, len, k, l, m, n )
  implicit none
  character*(*) str
  integer :: len
  integer :: k, l, m, n
end subroutine xxclose
