        module cnvpert_tl

        use precision
        use prognostics
        use prognostics_tl
        use stepon, only : ng_d, ng_s, coslon, sinlon

        PRIVATE
        PUBLIC g5tog4_tl
        PUBLIC g4tog5_tl

        contains

        subroutine g5tog4_tl ( prog, pert )

        implicit none

        type(dyn_prog) :: prog
        type(dyn_prog) :: pert

        real(r8) ua_tl(imr,jfirst-1:jlast,nl) ! U-Wind on A-Grid
        real(r8) va_tl(imr,jfirst:jlast,nl)   ! V-Wind on A-Grid


!       If so, convert initial perturbation to GEOS-4 type perturbation
!       ---------------------------------------------------------------
!??         call ps2delp_tl ( im, jfirst, jlast, nl, pert%delp )
        call t2th_tl ( prog%delp, pert%delp, prog%pt, pert%pt )
        ua_tl(:,jfirst-1:jlast,:) = pert%u(:,jfirst-1:jlast,:)
        va_tl(:,jfirst  :jlast,:) = pert%v(:,jfirst  :jlast,:)
        call a2d3d  ( ua_tl, va_tl, pert%u, pert%v, &
                      imr, jnp, nl, jfirst, jlast, ng_d, ng_s, coslon, sinlon )

        end subroutine g5tog4_tl

        subroutine g4tog5_tl ( prog, pert )

        implicit none

        type(dyn_prog) :: prog
        type(dyn_prog) :: pert

        real(r8) ua_tl(imr,jfirst-1:jlast,nl) ! U-Wind on A-Grid
        real(r8) va_tl(imr,jfirst:jlast,nl)   ! V-Wind on A-Grid

!       If so, convert final perturbation to GEOS-5 type perturbation
!       ---------------------------------------------------------------
        call th2t_tl ( prog%delp, pert%delp, prog%pt, pert%pt )
        call d2a3d ( pert%u(:,jfirst:jlast+1,:), pert%v, ua_tl(:,jfirst:jlast,:), va_tl, &
                     imr, jnp, nl, jfirst, jlast, ng_d, ng_s, coslon, sinlon )
        pert%u(:,jfirst-ng_d:jfirst-ng_d+2,:) = 0._r8
        pert%v(:,jfirst-ng_s:jfirst     -2,:) = 0._r8
        pert%v(:,jlast    +1:jlast   +ng_d,:) = 0._r8
        pert%u(:,jfirst-1:jlast,:) = ua_tl(:,jfirst-1:jlast,:)
        pert%v(:,jfirst  :jlast,:) = va_tl(:,jfirst  :jlast,:)
!??       call delp2ps_tl ( im, jfirst, jlast, nl, pert%delp )

        end subroutine g4tog5_tl

        end module cnvpert_tl
