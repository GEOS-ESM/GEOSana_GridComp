module  prognostics_q

use precision
use m_random, only : zufalli, normalen
use mod_comm, only : gid,mp_add1d
use stepon,   only : imr
use stepon,   only : jnp
use stepon,   only : nl
use stepon,   only : nc
use stepon,   only : jfirst
use stepon,   only : jlast
use stepon,   only : kfirst
use stepon,   only : klast
use stepon,   only : ng_d
use stepon,   only : jord
use m_die, only : die

!==============================================
! all entries are defined explicitly
!==============================================
implicit none

PRIVATE
SAVE

public prognostics_q_dotp
public prognostics_q_dup
public prognostics_q_final
public prognostics_q_initial
public prognostics_q_rand
public prognostics_q_zero
public dynq

character(len=*), parameter :: myname = 'prognostics_q'

!==============================================
! declare local variables
!==============================================
type dynq
real(kind=r8), pointer :: q(:,:,:,:)
end type dynq

interface prognostics_q_rand; module procedure &
          rand_
end interface

contains
subroutine prognostics_q_final ( prog )
!==============================================
! all entries are defined explicitly
!==============================================
implicit none

type(dynq) prog

!----------------------------------------------
! DEALLOCATE MODULE VARIABLES
!----------------------------------------------
deallocate( prog%q )

end subroutine prognostics_q_final

subroutine prognostics_q_initial ( prog )
!==============================================
! all entries are defined explicitly
!==============================================
implicit none

type(dynq) prog
character(len=*), parameter :: myname_ = myname // '*prognostics_q_initial'
integer ng_d, ng_s,klastp,ierr

#if defined( SPMD )
        ng_d = max(2,min(abs(jord),3))
        ng_s = 3
#else
        ng_d = 0
        ng_s = 0
#endif
      klastp = klast
      if (klast .eq. nl) klastp = klast + 1

!----------------------------------------------
! ALLOCATE AND RESET MODULE VARIABLES
!----------------------------------------------
allocate( prog%q(imr,jfirst-ng_d:jlast+ng_d,kfirst:klast,nc), stat=ierr )
         if(ierr/=0) call die(myname_,'Alloc(q)',ierr)
prog%q(:,:,:,:) = 0.

end subroutine prognostics_q_initial

subroutine prognostics_q_zero ( prog )
      implicit none
      type(dynq) :: prog
      prog%q    = 0.0d0
end subroutine prognostics_q_zero
subroutine prognostics_q_dup ( x, y )
implicit none
type(dynq) :: x
type(dynq) :: y
integer i,j,k,n
do n = 1, nc
do k = kfirst, klast
   do j = jfirst,jlast
      do i = 1, imr
         y%q(i,j,k,n) = x%q(i,j,k,n)
      enddo
   enddo
enddo
enddo
end subroutine prognostics_q_dup
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: rand_ --- random dyn_prog
!
! !INTERFACE:
subroutine rand_ ( prog, seed )
      implicit none
      integer, intent(in),optional :: seed
      type(dynq) :: prog
      integer seed_,n
      seed_ = 0
      if(present(seed)) seed_ = seed
      n = imr*(jlast-jfirst+1)*nl
      call zufalli ( seed_ )       ! initialize random number generator
      n = n*nc
      call normalen(n,prog%q)
end subroutine rand_

!-----------------------------------------------------------------------
real(r8) function prognostics_q_dotp ( x, y )
implicit none
type(dynq) :: x
type(dynq) :: y
integer, parameter :: ng = 1
real(r8) sum(nc), sgp(ng)
integer i,j,k,n
integer ntot,npptot
 
ntot= imr*jnp*klast*nc

sum = 0._r8
sgp = 0._r8
do n = 1, nc
do k = kfirst, klast
   do j = jfirst,jlast
      do i = 1, imr
         sum(n) = sum(n) + y%q(i,j,k,n) * x%q(i,j,k,n)
         sgp(1) = sgp(1) + 1
      enddo
   enddo
enddo
enddo
 
  call mp_add1d ( nc, sum )
  call mp_add1d ( ng, sgp )
  
  prognostics_q_dotp=0._r8
  do i=1,nc
    prognostics_q_dotp = prognostics_q_dotp + sum(i)
  enddo
  npptot = nint(sgp(1))
  if(ntot/=npptot)then
     print *, '            total # grid points*nc: ', ntot
     print *, 'counting pp total # grid points*nc: ', sgp(1)
     call die('prognostics_dopt',' inconsistent values ... aborting!')
  endif
 
end function prognostics_q_dotp

end module     prognostics_q


