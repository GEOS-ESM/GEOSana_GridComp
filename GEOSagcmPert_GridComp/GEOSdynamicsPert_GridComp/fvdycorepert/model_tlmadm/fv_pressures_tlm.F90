module fv_pressures_tlm_mod

use fv_arrays_mod, only: fv_grid_bounds_type

implicit none
private
public compute_pressures_tlm

contains

   SUBROUTINE COMPUTE_PRESSURES_TLM( bd, npz, kappa, ptop, delp, delp_tl, &
                                     pe, pe_tl, pk, pk_tl, pkz, pkz_tl, peln, peln_tl )
    IMPLICIT NONE
!Arguments
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npz
    REAL, INTENT(IN) :: kappa, ptop
    REAL, INTENT(IN) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(IN) :: delp_tl(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL, INTENT(OUT) :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(OUT) :: pe_tl(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL, INTENT(OUT) :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL, INTENT(OUT) :: pk_tl(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL, INTENT(OUT) :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL, INTENT(OUT) :: peln_tl(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL, INTENT(OUT) :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
    REAL, INTENT(OUT) :: pkz_tl(bd%is:bd%ie, bd%js:bd%je, npz)
!Locals
    INTEGER :: i, j, k, is, ie, js, je
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    pe(:, :, :) = 0.0
    pe(:, 1, :) = ptop
    pe_tl = 0.0
    DO k=2,npz+1
      DO j=js,je
        DO i=is,ie
          pe_tl(i, k, j) = pe_tl(i, k-1, j) + delp_tl(i, j, k-1)
          pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
        END DO
      END DO
    END DO
    peln_tl = 0.0
    DO k=1,npz+1
      DO j=js,je
        DO i=is,ie
          peln_tl(i, k, j) = pe_tl(i, k, j)/pe(i, k, j)
          peln(i, k, j) = LOG(pe(i, k, j))
        END DO
      END DO
    END DO
    pk_tl = 0.0
    DO k=1,npz+1
      DO j=js,je
        DO i=is,ie
          pk_tl(i, j, k) = kappa*peln_tl(i, k, j)*EXP(kappa*peln(i, k, j))
          pk(i, j, k) = EXP(kappa*peln(i, k, j))
        END DO
      END DO
    END DO
    pkz_tl = 0.0
    DO k=1,npz
      DO j=js,je
        DO i=is,ie
          pkz_tl(i, j, k) = ((pk_tl(i, j, k+1)-pk_tl(i, j, k))*kappa*(&
&           peln(i, k+1, j)-peln(i, k, j))-(pk(i, j, k+1)-pk(i, j, k))*&
&           kappa*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(kappa*(peln(i&
&           , k+1, j)-peln(i, k, j)))**2
          pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(kappa*(peln(i, k+1&
&           , j)-peln(i, k, j)))
        END DO
      END DO
    END DO
  END SUBROUTINE COMPUTE_PRESSURES_TLM

end module fv_pressures_tlm_mod
