      module m_getModConst

      implicit none

      PRIVATE

      PUBLIC get_modconst

!     Module to allow using fixed-format *.h model files
!     containing physical constants used by the model in
!     free-format files.
!       This should/could be merged with hermes/m_const.
!     --------------------------------------------------

      CONTAINS

      subroutine get_modconst ( )
      implicit none
      include "comcon.h"
      end subroutine get_modconst

      end module m_getModConst
