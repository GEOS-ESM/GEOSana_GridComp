!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  slp_das --- Computes sea-level pressure ("UKMO" algorithm)
! 
! !INTERFACE:
!
      subroutine slp_das (im, km, ps, phis, slp, pe, tv, rgas, grav)

! !USES:
      use precision

      implicit none
      real(r8) gamma
      parameter (gamma = 6.5e-3)

! !INPUT PARAMETERS: 

      integer(i4), intent(in)   :: im, km       ! grid dimensions

      real(r8) grav     
      real(r8) rgas
      real(r8),    intent(in)   ::  tv(im,km)    ! layer-mean v. Temp
      real(r8),    intent(in)   ::  pe(im,km+1)  ! press at layer edges (Pa)
      real(r8),    intent(in)   ::    ps(im)     ! surface pressure (Pa)
      real(r8),    intent(in)   ::  phis(im)     ! surface geopotential

! !OUTPUT PARAMETERS: 

      real(r8),    intent(out)  ::   slp(im)     ! sea-level pressure (hPa)

! !DESCRIPTION: Let's assume that the virtual temperature at the {\em fictious}
!               layer under the ground is given by
!  \be
!        T(z) = T_* + \gamma z, \qquad \gamma \equiv = 6.5 K/Km
!  \ee
!  where $T_*$ is a temperature assumed to be valid at the surface, but different
!  from model's ground temperature (which is heavily affected by the diurnal cycle.)
!  Integrating the hydrostatic equation with this particular profile of $T(z)$
!  we obtain
!  \be
!       p_0 = p_s \( {T_* + \gamma z_s \over T_*} \)^{g/(R\gamma)}
!  \ee
!  where $z_s$ is the model's topographic height, and $p_0$ is the desired sea-level
!  pressure. 
!
!  In order to avoid spurious diurnal cycle in the sea-level pressure, we do not
!  use the model's ground temperature for $T_*$. Instead, we extrapolated the
!  virtual temperature from a reference level just above the PBL (say, the fifth
!  eta layer, $p ~700$ hPa) using the constant moist adiabatic lapse rate $\gamma$,
!  viz.
!  \be
!          T_* = T_{ref} \( p_s / p_{ref} \) ^{g/(R\gamma)}
!  \ee
!  where $T_{ref}$ is a reference temperature valid at level $p_{ref}$.
!  For additional information consult Woodage (1985, Meteor. Mag., 114, 1-13)
!  and Ingleby (1995, Wea. Forecasting, 10, 172-182).
!
! !REVISION HISTORY: 
!
!  20Mar2000  da Silva  Initial code.
!  23Mar2000  S.-J. Lin revised
!  08Jul2000  Sawyer    Use precision module
!
!EOP
!-------------------------------------------------------------------------
! Local
      integer(i4) i,k
      real(r8)  t_star                     ! extrapolated temperature (K) 
      real(r8) p_offset
      real(r8) p_bot
      real(r8) t_ref(im)                   ! Reference virtual temperature (K)
      real(r8) p_ref(im)                   ! Reference pressure level (Pa)
      real(r8) dp1, dp2
      integer(i4) k_bot, k1, k2
      real(r8) factor, yfactor
      real(r8) gg

      gg      = gamma / grav
      factor  = grav / ( Rgas * gamma ) 
      yfactor =  Rgas * gg

! Find reference temperature by averaging Tv in a couple of
! layers above the PBL.

      p_offset = 15000.                     ! 150 hPa above surface

      do i = 1, im

         p_bot = ps(i) - p_offset
         k_bot = -1

         do k = km, 2, -1
            if ( pe(i,k+1) .lt. p_bot ) then
                 k_bot = k
                 go to 123
            endif
          end do

123       continue

          k1 = k_bot - 1
          k2 = k_bot

          dp1 = pe(i,k_bot) - pe(i,k_bot-1)
          dp2 = pe(i,k_bot+1) - pe(i,k_bot)

          t_ref(i) = ( tv(i,k1)*dp1 + tv(i,k2)*dp2 ) / (dp1+dp2)
          p_ref(i) = 0.5 * ( pe(i,k_bot+1) + pe(i,k_bot-1) )
      end do

      do i = 1, im
         t_star  = t_ref(i) * ( ps(i)/p_ref(i) )**yfactor
         slp(i)  = ps(i) * (1.0 + gg*phis(i)/t_star) ** factor
      end do

      return
      end 
