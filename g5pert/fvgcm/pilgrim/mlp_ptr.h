#if defined ( LAHEY )
#if defined ( PTR_INT )
      integer ptrg_4d
      integer ptrg_t1
      integer ptrg_t2
      integer ptrg_t3
      integer ptrg_2d
      integer ptrg_1d
#else

!
! Main vars:
!
      pointer (wing_4d, g_4d)
      real :: g_4d(FVGCM_LON, FVGCM_LAT, FVGCM_LEV, max_nq)

! Other work arrays:
!
! Type 1: For variables defined at layer edge (wz & pk)
!
      pointer (wing_t1, g_t1)
      real :: g_t1(FVGCM_LON, FVGCM_LAT, FVGCM_LEV+1, nbuf)

!
! Type 2: For edge pressure (pe)
!
      pointer (wing_t2, g_t2)
      real :: g_t2(FVGCM_LON, FVGCM_LEV+1, FVGCM_LAT)
!
! Type 3: 
!
      pointer (wing_t3, g_t3)
      real :: g_t3(FVGCM_LEV, maxpro)

!
! General purpose 2D (x-y) array
!
      pointer (wing_2d, g_2d)
      real :: g_2d(FVGCM_LON,FVGCM_LAT)
!
! General purpose 1D (in y) array
!
      pointer (wing_1d, g_1d)
      real :: g_1d(FVGCM_LAT)

#endif

#if ( ! defined NOT_ASSIGNED )
      wing_4d=ptrg_4d
      wing_t1=ptrg_t1
      wing_t2=ptrg_t2
      wing_t3=ptrg_t3
      wing_2d=ptrg_2d
      wing_1d=ptrg_1d
#endif
#endif
