module fv_pressures_adm_mod

use fv_arrays_mod, only: fv_grid_bounds_type

implicit none
private
public compute_pressures_fwd, compute_pressures_bwd

contains

  SUBROUTINE COMPUTE_PRESSURES_FWD(bd, npz, kappa, ptop, delp, pe, pk, &
&   pkz, peln)
    IMPLICIT NONE
!Arguments
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npz
    REAL, INTENT(IN) :: kappa, ptop
    REAL, INTENT(IN) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
!Locals
    INTEGER :: i, j, k, is, ie, js, je
    INTRINSIC LOG
    INTRINSIC EXP
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    pe(:, :, :) = 0.0
    pe(:, 1, :) = ptop
    DO k=2,npz+1
      DO j=js,je
        DO i=is,ie
          pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
        END DO
      END DO
    END DO
    DO k=1,npz+1
      DO j=js,je
        DO i=is,ie
          peln(i, k, j) = LOG(pe(i, k, j))
        END DO
      END DO
    END DO
    DO k=1,npz+1
      DO j=js,je
        DO i=is,ie
          pk(i, j, k) = EXP(kappa*peln(i, k, j))
        END DO
      END DO
    END DO
    DO k=1,npz
      DO j=js,je
        DO i=is,ie
          pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(kappa*(peln(i, k+1&
&           , j)-peln(i, k, j)))
        END DO
      END DO
    END DO
  END SUBROUTINE COMPUTE_PRESSURES_FWD

  SUBROUTINE COMPUTE_PRESSURES_BWD(bd, npz, kappa, ptop, delp, delp_ad, &
&   pe, pe_ad, pk, pk_ad, pkz, pkz_ad, peln, peln_ad)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: npz
    REAL, INTENT(IN) :: kappa, ptop
    REAL, INTENT(IN) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: delp_ad(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    REAL :: pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL :: pe_ad(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1)
    REAL :: pk(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL :: pk_ad(bd%is:bd%ie, bd%js:bd%je, npz+1)
    REAL :: peln(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL :: peln_ad(bd%is:bd%ie, npz+1, bd%js:bd%je)
    REAL :: pkz(bd%is:bd%ie, bd%js:bd%je, npz)
    REAL :: pkz_ad(bd%is:bd%ie, bd%js:bd%je, npz)
    INTEGER :: i, j, k, is, ie, js, je
    INTRINSIC LOG
    INTRINSIC EXP
    REAL :: temp
    REAL :: temp_ad
    REAL :: temp_ad0
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    DO k=npz,1,-1
      DO j=je,js,-1
        DO i=ie,is,-1
          temp = kappa*(peln(i, k+1, j)-peln(i, k, j))
          temp_ad = pkz_ad(i, j, k)/temp
          temp_ad0 = -((pk(i, j, k+1)-pk(i, j, k))*kappa*temp_ad/temp)
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp_ad
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp_ad
          peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad0
          peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad0
          pkz_ad(i, j, k) = 0.0
        END DO
      END DO
    END DO
    DO k=npz+1,1,-1
      DO j=je,js,-1
        DO i=ie,is,-1
          peln_ad(i, k, j) = peln_ad(i, k, j) + EXP(kappa*peln(i, k, j))&
&           *kappa*pk_ad(i, j, k)
          pk_ad(i, j, k) = 0.0
        END DO
      END DO
    END DO
    DO k=npz+1,1,-1
      DO j=je,js,-1
        DO i=ie,is,-1
          pe_ad(i, k, j) = pe_ad(i, k, j) + peln_ad(i, k, j)/pe(i, k, j)
          peln_ad(i, k, j) = 0.0
        END DO
      END DO
    END DO
    DO k=npz+1,2,-1
      DO j=je,js,-1
        DO i=ie,is,-1
          pe_ad(i, k-1, j) = pe_ad(i, k-1, j) + pe_ad(i, k, j)
          delp_ad(i, j, k-1) = delp_ad(i, j, k-1) + pe_ad(i, k, j)
          pe_ad(i, k, j) = 0.0
        END DO
      END DO
    END DO
  END SUBROUTINE COMPUTE_PRESSURES_BWD

end module fv_pressures_adm_mod
