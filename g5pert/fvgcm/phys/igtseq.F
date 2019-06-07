#include <params.h>
      character*8 function igtseq()
C-----------------------------------------------------------------------
C
C Return job sequence number as character variable
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath 
C Standardized:      L. Bath, Jun 1992
C                    L. Buja, Feb 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C---------------------------Local variables-----------------------------
C
#if ( defined CRAY )
      character*8 temps   ! Temporary character variable
#endif
C
C-----------------------------------------------------------------------
C
      igtseq = '        '
#if ( defined CRAY )
      call getenv('QSUB_REQNAME',temps)
      call ljust(temps,igtseq)         ! Left-justify character variable
#endif
      return
      end
 
