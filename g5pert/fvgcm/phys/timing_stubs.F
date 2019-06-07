      subroutine t_startf(xxx)
C-----------------------------------------------------------------------
C
C Timing stubs.  Routines t_startf, t_stopf, and t_prf can be replaced by 
C rountines of the same name to initiate, terminate, and print the results of
C timing portions of the CCM3 code.
C
C-------------------------Code History----------------------------------
C
C Original version:  J. Rosinski, Oct 1995
C Standardized:      T. Acker, Feb 1996
C Reviewed:          
C 
C-----------------------------------------------------------------------
c
c $Id $
c
      character *(*) xxx
      return
      end

      subroutine t_stopf(xxx)
      character *(*) xxx
      return
      end

      subroutine t_prf(xxx)
      integer xxx
      return
      end
 
