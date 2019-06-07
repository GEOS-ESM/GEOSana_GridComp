!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
module fv_restart_tlm_mod

implicit none
private

public :: d2c_setup, d2a_setup
public :: d2c_setup_tlm

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

CONTAINS
!  Differentiation of d2c_setup in forward (tangent) mode:
!   variations   of useful results: uc vc
!   with respect to varying inputs: u v uc vc
  SUBROUTINE D2C_SETUP_TLM(u, u_tl, v, v_tl, ua, va, uc, uc_tl, vc, &
&   vc_tl, dord4, isd, ied, jsd, jed, is, ie, js, je, npx, npy, &
&   grid_type, nested, se_corner, sw_corner, ne_corner, nw_corner, &
&   rsin_u, rsin_v, cosa_s, rsin2)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: dord4
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, npx, npy&
&   , grid_type
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: u_tl(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: v_tl(isd:ied+1, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: va
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(OUT) :: uc_tl
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(OUT) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(OUT) :: vc_tl
    LOGICAL, INTENT(IN) :: nested, se_corner, sw_corner, ne_corner, &
&   nw_corner
    REAL, INTENT(IN) :: rsin_u(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rsin_v(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: cosa_s(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rsin2(isd:ied, jsd:jed)
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp_tl, vtmp_tl
    REAL, PARAMETER :: t11=27./28., t12=-(13./28.), t13=3./7., t14=6./7.&
&   , t15=3./28.
    REAL, PARAMETER :: a1=0.5625
    REAL, PARAMETER :: a2=-0.0625
    REAL, PARAMETER :: c1=-(2./14.)
    REAL, PARAMETER :: c2=11./14.
    REAL, PARAMETER :: c3=5./14.
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    IF (dord4) THEN
      id = 1
    ELSE
      id = 0
    END IF
    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
    IF (nested) THEN
      utmp_tl = 0.0
      DO j=jsd+1,jed-1
        DO i=isd,ied
          utmp_tl(i, j) = a2*(u_tl(i, j-1)+u_tl(i, j+2)) + a1*(u_tl(i, j&
&           )+u_tl(i, j+1))
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      DO i=isd,ied
!j = jsd
        utmp_tl(i, jsd) = 0.5*(u_tl(i, jsd)+u_tl(i, jsd+1))
        utmp(i, jsd) = 0.5*(u(i, jsd)+u(i, jsd+1))
!j = jed
        utmp_tl(i, jed) = 0.5*(u_tl(i, jed)+u_tl(i, jed+1))
        utmp(i, jed) = 0.5*(u(i, jed)+u(i, jed+1))
      END DO
      vtmp_tl = 0.0
      DO j=jsd,jed
        DO i=isd+1,ied-1
          vtmp_tl(i, j) = a2*(v_tl(i-1, j)+v_tl(i+2, j)) + a1*(v_tl(i, j&
&           )+v_tl(i+1, j))
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
!i = isd
        vtmp_tl(isd, j) = 0.5*(v_tl(isd, j)+v_tl(isd+1, j))
        vtmp(isd, j) = 0.5*(v(isd, j)+v(isd+1, j))
!i = ied
        vtmp_tl(ied, j) = 0.5*(v_tl(ied, j)+v_tl(ied+1, j))
        vtmp(ied, j) = 0.5*(v(ied, j)+v(ied+1, j))
      END DO
      DO j=jsd,jed
        DO i=isd,ied
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        END DO
      END DO
    ELSE
!----------
! Interior:
!----------
      utmp = 0.
      vtmp = 0.
      IF (npt .LT. js - 1) THEN
        max1 = js - 1
      ELSE
        max1 = npt
      END IF
      IF (npy - npt .GT. je + 1) THEN
        min1 = je + 1
        utmp_tl = 0.0
      ELSE
        min1 = npy - npt
        utmp_tl = 0.0
      END IF
      DO j=max1,min1
        IF (npt .LT. isd) THEN
          max2 = isd
        ELSE
          max2 = npt
        END IF
        IF (npx - npt .GT. ied) THEN
          min2 = ied
        ELSE
          min2 = npx - npt
        END IF
        DO i=max2,min2
          utmp_tl(i, j) = a2*(u_tl(i, j-1)+u_tl(i, j+2)) + a1*(u_tl(i, j&
&           )+u_tl(i, j+1))
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      IF (npt .LT. jsd) THEN
        max3 = jsd
      ELSE
        max3 = npt
      END IF
      IF (npy - npt .GT. jed) THEN
        min3 = jed
        vtmp_tl = 0.0
      ELSE
        min3 = npy - npt
        vtmp_tl = 0.0
      END IF
      DO j=max3,min3
        IF (npt .LT. is - 1) THEN
          max4 = is - 1
        ELSE
          max4 = npt
        END IF
        IF (npx - npt .GT. ie + 1) THEN
          min4 = ie + 1
        ELSE
          min4 = npx - npt
        END IF
        DO i=max4,min4
          vtmp_tl(i, j) = a2*(v_tl(i-1, j)+v_tl(i+2, j)) + a1*(v_tl(i, j&
&           )+v_tl(i+1, j))
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
      END DO
!----------
! edges:
!----------
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
          DO j=jsd,npt-1
            DO i=isd,ied
              utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
          DO j=npy-npt+1,jed
            DO i=isd,ied
              utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (is .EQ. 1 .OR. isd .LT. npt) THEN
          IF (npt .LT. jsd) THEN
            max5 = jsd
          ELSE
            max5 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min5 = jed
          ELSE
            min5 = npy - npt
          END IF
          DO j=max5,min5
            DO i=isd,npt-1
              utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
          IF (npt .LT. jsd) THEN
            max6 = jsd
          ELSE
            max6 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min6 = jed
          ELSE
            min6 = npy - npt
          END IF
          DO j=max6,min6
            DO i=npx-npt+1,ied
              utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
      END IF
      DO j=js-1-id,je+1+id
        DO i=is-1-id,ie+1+id
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        END DO
      END DO
    END IF
! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
    IF (sw_corner) THEN
      DO i=-2,0
        utmp_tl(i, 0) = -vtmp_tl(0, 1-i)
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
    END IF
    IF (se_corner) THEN
      DO i=0,2
        utmp_tl(npx+i, 0) = vtmp_tl(npx, i+1)
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
    END IF
    IF (ne_corner) THEN
      DO i=0,2
        utmp_tl(npx+i, npy) = -vtmp_tl(npx, je-i)
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
    END IF
    IF (nw_corner) THEN
      DO i=-2,0
        utmp_tl(i, npy) = vtmp_tl(0, je+i)
        utmp(i, npy) = vtmp(0, je+i)
      END DO
    END IF
    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      IF (3 .LT. is - 1) THEN
        ifirst = is - 1
      ELSE
        ifirst = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        ilast = ie + 2
      ELSE
        ilast = npx - 2
      END IF
    ELSE
      ifirst = is - 1
      ilast = ie + 2
    END IF
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
    DO j=js-1,je+1
      DO i=ifirst,ilast
        uc_tl(i, j) = a1*(utmp_tl(i-1, j)+utmp_tl(i, j)) + a2*(utmp_tl(i&
&         -2, j)+utmp_tl(i+1, j))
        uc(i, j) = a1*(utmp(i-1, j)+utmp(i, j)) + a2*(utmp(i-2, j)+utmp(&
&         i+1, j))
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
! Xdir:
      IF (is .EQ. 1 .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc_tl(0, j) = c1*utmp_tl(-2, j) + c2*utmp_tl(-1, j) + c3*&
&           utmp_tl(0, j)
          uc(0, j) = c1*utmp(-2, j) + c2*utmp(-1, j) + c3*utmp(0, j)
          uc_tl(1, j) = rsin_u(1, j)*(t14*(utmp_tl(0, j)+utmp_tl(1, j))+&
&           t12*(utmp_tl(-1, j)+utmp_tl(2, j))+t15*(utmp_tl(-2, j)+&
&           utmp_tl(3, j)))
          uc(1, j) = (t14*(utmp(0, j)+utmp(1, j))+t12*(utmp(-1, j)+utmp(&
&           2, j))+t15*(utmp(-2, j)+utmp(3, j)))*rsin_u(1, j)
          uc_tl(2, j) = c1*utmp_tl(3, j) + c2*utmp_tl(2, j) + c3*utmp_tl&
&           (1, j)
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc_tl(npx-1, j) = c1*utmp_tl(npx-3, j) + c2*utmp_tl(npx-2, j) &
&           + c3*utmp_tl(npx-1, j)
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
          uc_tl(npx, j) = rsin_u(npx, j)*(t14*(utmp_tl(npx-1, j)+utmp_tl&
&           (npx, j))+t12*(utmp_tl(npx-2, j)+utmp_tl(npx+1, j))+t15*(&
&           utmp_tl(npx-3, j)+utmp_tl(npx+2, j)))
          uc(npx, j) = (t14*(utmp(npx-1, j)+utmp(npx, j))+t12*(utmp(npx-&
&           2, j)+utmp(npx+1, j))+t15*(utmp(npx-3, j)+utmp(npx+2, j)))*&
&           rsin_u(npx, j)
          uc_tl(npx+1, j) = c3*utmp_tl(npx, j) + c2*utmp_tl(npx+1, j) + &
&           c1*utmp_tl(npx+2, j)
          uc(npx+1, j) = c3*utmp(npx, j) + c2*utmp(npx+1, j) + c1*utmp(&
&           npx+2, j)
        END DO
      END IF
    END IF
!------
! Ydir:
!------
    IF (sw_corner) THEN
      DO j=-2,0
        vtmp_tl(0, j) = -utmp_tl(1-j, 0)
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
    END IF
    IF (nw_corner) THEN
      DO j=0,2
        vtmp_tl(0, npy+j) = utmp_tl(j+1, npy)
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
    END IF
    IF (se_corner) THEN
      DO j=-2,0
        vtmp_tl(npx, j) = utmp_tl(ie+j, 0)
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
    END IF
    IF (ne_corner) THEN
      DO j=0,2
        vtmp_tl(npx, npy+j) = -utmp_tl(ie-j, npy)
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1 .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vc_tl(i, 1) = rsin_v(i, 1)*(t14*(vtmp_tl(i, 0)+vtmp_tl(i, 1)&
&             )+t12*(vtmp_tl(i, -1)+vtmp_tl(i, 2))+t15*(vtmp_tl(i, -2)+&
&             vtmp_tl(i, 3)))
            vc(i, 1) = (t14*(vtmp(i, 0)+vtmp(i, 1))+t12*(vtmp(i, -1)+&
&             vtmp(i, 2))+t15*(vtmp(i, -2)+vtmp(i, 3)))*rsin_v(i, 1)
          END DO
        ELSE IF ((j .EQ. 0 .OR. j .EQ. npy - 1) .AND. (.NOT.nested)) &
&       THEN
          DO i=is-1,ie+1
            vc_tl(i, j) = c1*vtmp_tl(i, j-2) + c2*vtmp_tl(i, j-1) + c3*&
&             vtmp_tl(i, j)
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
          END DO
        ELSE IF ((j .EQ. 2 .OR. j .EQ. npy + 1) .AND. (.NOT.nested)) &
&       THEN
          DO i=is-1,ie+1
            vc_tl(i, j) = c1*vtmp_tl(i, j+1) + c2*vtmp_tl(i, j) + c3*&
&             vtmp_tl(i, j-1)
            vc(i, j) = c1*vtmp(i, j+1) + c2*vtmp(i, j) + c3*vtmp(i, j-1)
          END DO
        ELSE IF (j .EQ. npy .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vc_tl(i, npy) = rsin_v(i, npy)*(t14*(vtmp_tl(i, npy-1)+&
&             vtmp_tl(i, npy))+t12*(vtmp_tl(i, npy-2)+vtmp_tl(i, npy+1))&
&             +t15*(vtmp_tl(i, npy-3)+vtmp_tl(i, npy+2)))
            vc(i, npy) = (t14*(vtmp(i, npy-1)+vtmp(i, npy))+t12*(vtmp(i&
&             , npy-2)+vtmp(i, npy+1))+t15*(vtmp(i, npy-3)+vtmp(i, npy+2&
&             )))*rsin_v(i, npy)
          END DO
        ELSE
! 4th order interpolation for interior points:
          DO i=is-1,ie+1
            vc_tl(i, j) = a2*(vtmp_tl(i, j-2)+vtmp_tl(i, j+1)) + a1*(&
&             vtmp_tl(i, j-1)+vtmp_tl(i, j))
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
          END DO
        END IF
      END DO
    ELSE
! 4th order interpolation:
      DO j=js-1,je+2
        DO i=is-1,ie+1
          vc_tl(i, j) = a2*(vtmp_tl(i, j-2)+vtmp_tl(i, j+1)) + a1*(&
&           vtmp_tl(i, j-1)+vtmp_tl(i, j))
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
        END DO
      END DO
    END IF
  END SUBROUTINE D2C_SETUP_TLM
  SUBROUTINE D2C_SETUP(u, v, ua, va, uc, vc, dord4, isd, ied, jsd, jed, &
&   is, ie, js, je, npx, npy, grid_type, nested, se_corner, sw_corner, &
&   ne_corner, nw_corner, rsin_u, rsin_v, cosa_s, rsin2)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: dord4
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, npx, npy&
&   , grid_type
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: va
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(OUT) :: vc
    LOGICAL, INTENT(IN) :: nested, se_corner, sw_corner, ne_corner, &
&   nw_corner
    REAL, INTENT(IN) :: rsin_u(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rsin_v(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: cosa_s(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rsin2(isd:ied, jsd:jed)
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, PARAMETER :: t11=27./28., t12=-(13./28.), t13=3./7., t14=6./7.&
&   , t15=3./28.
    REAL, PARAMETER :: a1=0.5625
    REAL, PARAMETER :: a2=-0.0625
    REAL, PARAMETER :: c1=-(2./14.)
    REAL, PARAMETER :: c2=11./14.
    REAL, PARAMETER :: c3=5./14.
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    IF (dord4) THEN
      id = 1
    ELSE
      id = 0
    END IF
    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
    IF (nested) THEN
      DO j=jsd+1,jed-1
        DO i=isd,ied
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      DO i=isd,ied
!j = jsd
        utmp(i, jsd) = 0.5*(u(i, jsd)+u(i, jsd+1))
!j = jed
        utmp(i, jed) = 0.5*(u(i, jed)+u(i, jed+1))
      END DO
      DO j=jsd,jed
        DO i=isd+1,ied-1
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
!i = isd
        vtmp(isd, j) = 0.5*(v(isd, j)+v(isd+1, j))
!i = ied
        vtmp(ied, j) = 0.5*(v(ied, j)+v(ied+1, j))
      END DO
      DO j=jsd,jed
        DO i=isd,ied
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        END DO
      END DO
    ELSE
!----------
! Interior:
!----------
      utmp = 0.
      vtmp = 0.
      IF (npt .LT. js - 1) THEN
        max1 = js - 1
      ELSE
        max1 = npt
      END IF
      IF (npy - npt .GT. je + 1) THEN
        min1 = je + 1
      ELSE
        min1 = npy - npt
      END IF
      DO j=max1,min1
        IF (npt .LT. isd) THEN
          max2 = isd
        ELSE
          max2 = npt
        END IF
        IF (npx - npt .GT. ied) THEN
          min2 = ied
        ELSE
          min2 = npx - npt
        END IF
        DO i=max2,min2
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      IF (npt .LT. jsd) THEN
        max3 = jsd
      ELSE
        max3 = npt
      END IF
      IF (npy - npt .GT. jed) THEN
        min3 = jed
      ELSE
        min3 = npy - npt
      END IF
      DO j=max3,min3
        IF (npt .LT. is - 1) THEN
          max4 = is - 1
        ELSE
          max4 = npt
        END IF
        IF (npx - npt .GT. ie + 1) THEN
          min4 = ie + 1
        ELSE
          min4 = npx - npt
        END IF
        DO i=max4,min4
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
      END DO
!----------
! edges:
!----------
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
          DO j=jsd,npt-1
            DO i=isd,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
          DO j=npy-npt+1,jed
            DO i=isd,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (is .EQ. 1 .OR. isd .LT. npt) THEN
          IF (npt .LT. jsd) THEN
            max5 = jsd
          ELSE
            max5 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min5 = jed
          ELSE
            min5 = npy - npt
          END IF
          DO j=max5,min5
            DO i=isd,npt-1
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
          IF (npt .LT. jsd) THEN
            max6 = jsd
          ELSE
            max6 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min6 = jed
          ELSE
            min6 = npy - npt
          END IF
          DO j=max6,min6
            DO i=npx-npt+1,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
      END IF
      DO j=js-1-id,je+1+id
        DO i=is-1-id,ie+1+id
          ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
          va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        END DO
      END DO
    END IF
! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
    IF (sw_corner) THEN
      DO i=-2,0
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
    END IF
    IF (se_corner) THEN
      DO i=0,2
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
    END IF
    IF (ne_corner) THEN
      DO i=0,2
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
    END IF
    IF (nw_corner) THEN
      DO i=-2,0
        utmp(i, npy) = vtmp(0, je+i)
      END DO
    END IF
    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      IF (3 .LT. is - 1) THEN
        ifirst = is - 1
      ELSE
        ifirst = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        ilast = ie + 2
      ELSE
        ilast = npx - 2
      END IF
    ELSE
      ifirst = is - 1
      ilast = ie + 2
    END IF
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
    DO j=js-1,je+1
      DO i=ifirst,ilast
        uc(i, j) = a1*(utmp(i-1, j)+utmp(i, j)) + a2*(utmp(i-2, j)+utmp(&
&         i+1, j))
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
! Xdir:
      IF (is .EQ. 1 .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc(0, j) = c1*utmp(-2, j) + c2*utmp(-1, j) + c3*utmp(0, j)
          uc(1, j) = (t14*(utmp(0, j)+utmp(1, j))+t12*(utmp(-1, j)+utmp(&
&           2, j))+t15*(utmp(-2, j)+utmp(3, j)))*rsin_u(1, j)
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx .AND. (.NOT.nested)) THEN
        DO j=js-1,je+1
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
          uc(npx, j) = (t14*(utmp(npx-1, j)+utmp(npx, j))+t12*(utmp(npx-&
&           2, j)+utmp(npx+1, j))+t15*(utmp(npx-3, j)+utmp(npx+2, j)))*&
&           rsin_u(npx, j)
          uc(npx+1, j) = c3*utmp(npx, j) + c2*utmp(npx+1, j) + c1*utmp(&
&           npx+2, j)
        END DO
      END IF
    END IF
!------
! Ydir:
!------
    IF (sw_corner) THEN
      DO j=-2,0
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
    END IF
    IF (nw_corner) THEN
      DO j=0,2
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
    END IF
    IF (se_corner) THEN
      DO j=-2,0
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
    END IF
    IF (ne_corner) THEN
      DO j=0,2
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1 .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vc(i, 1) = (t14*(vtmp(i, 0)+vtmp(i, 1))+t12*(vtmp(i, -1)+&
&             vtmp(i, 2))+t15*(vtmp(i, -2)+vtmp(i, 3)))*rsin_v(i, 1)
          END DO
        ELSE IF ((j .EQ. 0 .OR. j .EQ. npy - 1) .AND. (.NOT.nested)) &
&       THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
          END DO
        ELSE IF ((j .EQ. 2 .OR. j .EQ. npy + 1) .AND. (.NOT.nested)) &
&       THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j+1) + c2*vtmp(i, j) + c3*vtmp(i, j-1)
          END DO
        ELSE IF (j .EQ. npy .AND. (.NOT.nested)) THEN
          DO i=is-1,ie+1
            vc(i, npy) = (t14*(vtmp(i, npy-1)+vtmp(i, npy))+t12*(vtmp(i&
&             , npy-2)+vtmp(i, npy+1))+t15*(vtmp(i, npy-3)+vtmp(i, npy+2&
&             )))*rsin_v(i, npy)
          END DO
        ELSE
! 4th order interpolation for interior points:
          DO i=is-1,ie+1
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
          END DO
        END IF
      END DO
    ELSE
! 4th order interpolation:
      DO j=js-1,je+2
        DO i=is-1,ie+1
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
        END DO
      END DO
    END IF
  END SUBROUTINE D2C_SETUP
  SUBROUTINE D2A_SETUP(u, v, ua, va, dord4, isd, ied, jsd, jed, is, ie, &
&   js, je, npx, npy, grid_type, nested, cosa_s, rsin2)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: dord4
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, npx, npy&
&   , grid_type
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: va
    REAL, INTENT(IN) :: cosa_s(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rsin2(isd:ied, jsd:jed)
    LOGICAL, INTENT(IN) :: nested
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, PARAMETER :: t11=27./28., t12=-(13./28.), t13=3./7., t14=6./7.&
&   , t15=3./28.
    REAL, PARAMETER :: a1=0.5625
    REAL, PARAMETER :: a2=-0.0625
    REAL, PARAMETER :: c1=-(2./14.)
    REAL, PARAMETER :: c2=11./14.
    REAL, PARAMETER :: c3=5./14.
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    IF (dord4) THEN
      id = 1
    ELSE
      id = 0
    END IF
    IF (grid_type .LT. 3 .AND. (.NOT.nested)) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
    IF (nested) THEN
      DO j=jsd+1,jed-1
        DO i=isd,ied
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      DO i=isd,ied
!j = jsd
        utmp(i, jsd) = 0.5*(u(i, jsd)+u(i, jsd+1))
!j = jed
        utmp(i, jed) = 0.5*(u(i, jed)+u(i, jed+1))
      END DO
      DO j=jsd,jed
        DO i=isd+1,ied-1
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
!i = isd
        vtmp(isd, j) = 0.5*(v(isd, j)+v(isd+1, j))
!i = ied
        vtmp(ied, j) = 0.5*(v(ied, j)+v(ied+1, j))
      END DO
    ELSE
      IF (npt .LT. js - 1) THEN
        max1 = js - 1
      ELSE
        max1 = npt
      END IF
      IF (npy - npt .GT. je + 1) THEN
        min1 = je + 1
      ELSE
        min1 = npy - npt
      END IF
!----------
! Interior:
!----------
      DO j=max1,min1
        IF (npt .LT. isd) THEN
          max2 = isd
        ELSE
          max2 = npt
        END IF
        IF (npx - npt .GT. ied) THEN
          min2 = ied
        ELSE
          min2 = npx - npt
        END IF
        DO i=max2,min2
          utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
        END DO
      END DO
      IF (npt .LT. jsd) THEN
        max3 = jsd
      ELSE
        max3 = npt
      END IF
      IF (npy - npt .GT. jed) THEN
        min3 = jed
      ELSE
        min3 = npy - npt
      END IF
      DO j=max3,min3
        IF (npt .LT. is - 1) THEN
          max4 = is - 1
        ELSE
          max4 = npt
        END IF
        IF (npx - npt .GT. ie + 1) THEN
          min4 = ie + 1
        ELSE
          min4 = npx - npt
        END IF
        DO i=max4,min4
          vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
        END DO
      END DO
!----------
! edges:
!----------
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
          DO j=jsd,npt-1
            DO i=isd,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
          DO j=npy-npt+1,jed
            DO i=isd,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (is .EQ. 1 .OR. isd .LT. npt) THEN
          IF (npt .LT. jsd) THEN
            max5 = jsd
          ELSE
            max5 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min5 = jed
          ELSE
            min5 = npy - npt
          END IF
          DO j=max5,min5
            DO i=isd,npt-1
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
        IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
          IF (npt .LT. jsd) THEN
            max6 = jsd
          ELSE
            max6 = npt
          END IF
          IF (npy - npt .GT. jed) THEN
            min6 = jed
          ELSE
            min6 = npy - npt
          END IF
          DO j=max6,min6
            DO i=npx-npt+1,ied
              utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
              vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
            END DO
          END DO
        END IF
      END IF
    END IF
    DO j=js-1-id,je+1+id
      DO i=is-1-id,ie+1+id
        ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
      END DO
    END DO
  END SUBROUTINE D2A_SETUP

end module fv_restart_tlm_mod
