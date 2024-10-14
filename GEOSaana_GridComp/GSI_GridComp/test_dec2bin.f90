program test_dec2bin

implicit none

integer,parameter :: i_kind = 4
integer,parameter :: mdim = 5
integer,parameter :: ndim = 10
integer(i_kind) :: dec(mdim)
integer(i_kind) :: bin(ndim)

integer ii,input

dec(1) = -2
dec(2) = -1
dec(3) =  1
dec(4) =  0 
dec(5) = 31


print *, ' dec   bin '
do ii = 1,mdim
  input = dec(ii)
  call dec2bin_(dec(ii),bin,ndim)
  print*, input, dec(ii), bin
  print* 
enddo


contains
subroutine dec2bin_(dec,bin,ndim)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    dec2bin                  convert decimal number to binary
!   prgmmr: unknown             org: np23                date: 2010-04-06
!
! abstract:  This routine convert a decimal number to binary
!
! program history log:
!   2010-04-06  hliu
!   2013-02-05  guo  - STOP in dec2bin() was replaced with die() to signal an _abort_.
!
!   input argument list:
!     dec  - observation type to process
!
!   output argument list:
!     bin    - number of sbuv/omi ozone observations read
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!

!   use kinds, only: i_kind
!   use mpeu_util, only: die, perr

    implicit none

! Declare passed variables
    integer(i_kind) ,intent(inout) :: dec
    integer(i_kind) ,intent(in)    :: ndim
    integer(i_kind) ,intent(out)   :: bin(ndim)

! Declare local variables
    integer(i_kind):: bindec, i

!   Check to determine decimal # is within bounds
    i = ndim
    IF ((dec - 2**i) >= 0) THEN
       write(6,*) 'Error: Decimal Number too Large. Must be < 2^(',ndim-1,')'
       stop 99
    END IF

!   Determine the scalar for each of the decimal positions
    DO WHILE (i >= 1)
       bindec = 2**(i-1)
       IF ((dec - bindec) >= 0) THEN
          bin(i) = 1
          dec = dec - bindec
       ELSE
          bin(i) = 0
       END IF
       i = i - 1
    END DO
    RETURN
END subroutine dec2bin_

end program test_dec2bin
