!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_poles --- Handles pole values averaging for admtlm_ 
!
! !INTERFACE:
!
      module m_poles

! !USES:

      use precision
      use prognostics, only : dyn_prog

      implicit none

! !PUBLIC MEMBER FUNCTIONS:

      PRIVATE

      PUBLIC SetPoles 
 
      interface SetPoles; module procedure  &
                SetPoles_
      end interface
!
! !DESCRIPTION: This modules handles pole values averaging 
!               for TLM and ADM codes. 
!
! !REVISION HISTORY:
!
!  19Sep2005  Elena N.  Adopted from m_TrajMng module
!  14May2007  Todling   Introduced dyn_prog; global change.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'm_Poles'

      real(r8), pointer :: delp(:,:,:)   ! pressure thickness (pascal)
      real(r8), pointer ::    u(:,:,:)   ! zonal wind on D-grid
      real(r8), pointer ::    v(:,:,:)   ! meridional wind
      real(r8), pointer ::   pt(:,:,:)   ! virtual potential temperature
      real(r8), pointer ::    q(:,:,:,:) ! specific humidity & tracer mixing ratios
      real(r8), pointer ::   ps(:,:)     ! surface pressure (pascal)

      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: SetPoles_: Values at poles are set to be the same and equal to an average value 
!
! !INTERFACE:

      subroutine SetPoles_   ( prog, stat )

! !USES:

      use precision
      use prognostics, only : imr      ! no. of grid points in the zonal direction
      use prognostics, only : jnp      ! no. of grid points in the meridional direction
      use prognostics, only : nl       ! no. of levels
      use prognostics, only : nc       ! no. of tracers
      use prognostics, only : jfirst
      use prognostics, only : jlast

      use prognostics, only : ps

      use m_die, only : MP_die

      implicit none

! !INPUT PARAMETERS:

      integer, optional :: stat

      type(dyn_prog), intent(inout), TARGET  :: prog

! !DESCRIPTION: Assign the same values at both poles 
!
! !REVISION HISTORY:
!
!  19Sep2005 Elena N.  Adopted from  PutPert_ (Todling - Initial code)
!  13May2007 Todling   Introduced dyn_prog
!  30Dec2009 Kokron    Set ierr to zero at entry code
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*SetPoles_'

      real(r8) dbleimr, xsum 
      integer  j, k, l, ierr

      if(present(stat)) stat = 0
      ierr = 0
      dbleimr = dble(imr)

!     Set pointers
!     ------------
      delp => prog%delp
       u   => prog%u
       v   => prog%v
      pt   => prog%pt
      q    => prog%q

!     ps   => prog%ps


!     Assign same values at poles for all fields 
!     -----------------------------------
       if (jfirst == 1) then   
         j = 1
         do k = 1, nl
            xsum = sum(delp(:,j,k))/dbleimr
            delp(:,j,k) = xsum
            xsum = sum(pt(:,j,k))/dbleimr
            pt(:,j,k) = xsum  
            u(:,j,k) = 0.d0
            v(:,j,k) = 0.d0
            do l = 1, nc 
               xsum = sum(q(:,j,k,l))/dbleimr
               q(:,j,k,l) = xsum
            enddo
         enddo
         xsum = sum(ps(:,j))/dbleimr
         ps(:,j) = xsum
        endif
        if (jlast == jnp) then
          j = jnp
          do k = 1, nl
            xsum = sum(delp(:,j,k))/dbleimr
            delp(:,j,k) = xsum
            xsum = sum(pt(:,j,k))/dbleimr
            pt(:,j,k) = xsum
            v(:,j,k) = 0.d0
            do l = 1, nc
               xsum = sum(q(:,j,k,l))/dbleimr
               q(:,j,k,l) = xsum
            enddo
          enddo
          xsum = sum(ps(:,j))/dbleimr
          ps(:,j) = xsum
        endif
!   ---------------------------------------------

      if ( ierr .ne. 0 ) then
         if(present(stat))then
            stat = 1
            return
         else
           call MP_die (myname_,'Error from SetPoles_',ierr)
         endif
      end if

      end subroutine SetPoles_

      end module m_Poles
