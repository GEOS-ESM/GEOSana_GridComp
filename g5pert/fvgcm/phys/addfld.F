#include <params.h>
      subroutine addfld(fname,units,nf,lvlflag,actinact)
C-----------------------------------------------------------------------
C
C Add a field to the Master Field List
C
C-------------------------Code History----------------------------------
C
C Original version:  J.Rosinski
C Standardized:      T.Acker, Feb 1996
C Reviewed:          
C 
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C------------------------------Commons----------------------------------
#include <comhst.h>
C-----------------------------------------------------------------------
c
c Input Arguments
c
      character*(*) fname    ! Input: field name--should be 8 characters
      character*(*) units    ! Input: units of fname--should be 8 chars
      integer nf             ! Input/output: field counter
      integer lvlflag        ! Input: flag for single/multih/multif
      integer actinact       ! Input: flag for active/inactive
C
C-----------------------------------------------------------------------
C
C Add field to Master Field List arrays fieldn and iflds
C
      nf = nf + 1            ! Increase number of fields on Master F.L.
      fieldn(1,nf) = fname   ! Store field name in field list array
      fieldn(2,nf) = units   ! Store units in field list array
      iflds(1,nf)  = lvlflag ! Store number of levels in iflds array
      iflds(3,nf)  = actinact! Store active/inactive flag in iflds array
C
      return
C
      end

 
