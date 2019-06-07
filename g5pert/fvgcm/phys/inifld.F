#include <params.h>
      subroutine inifld(fieldn  ,iflds   ,hacc    ,hbufpt  ,nden    ,
     $                  ktape   ,nrlen   ,nplen   ,nflds   ,mfld    ,
     $                  mcfld   ,npt     ,nppt    )
C-----------------------------------------------------------------------
C
C Update history buffer pointers, packed and unpacked data record
C lengths, and the number of fields in current history file. Also,
C store field list information for field fname, as follows:
C
C  mflds(3,nflds)  (integer array)
C  -------------------------------
C  mflds(1,if) = layer location flag and instantaneous
C                or accumulated field flag.
C  mflds(2,if) = pointer to first value for field(if) in data
C                record
C  mflds(3,if) = data packing density (1-4)
C
C
C  mcflds(2,nflds)  (character array)
C  ---------------------------------
C  mcflds(1,if) = name of field(if)
C  mcflds(2,if) = units of field(if)
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath  
C Standardized:      L. Bath, June 1992
C                    T. Acker, March 1996
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
C-----------------------------------------------------------------------
#include <comlun.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      character*8 fieldn(2)   ! Field name and units (from bldfld)
      integer iflds(3)        ! Single vs. multi-level flag and default
C                             ! h-file active/inactive flag (from bldfld)
      character*1 hacc        ! Accumulation indicator
C                             !   'I' = instantaneous, 'A' = averaged
C                             !   'M' = minimum,       'X' = maximum
      integer nden            !  Packing density, 1,2,3 or 4
      integer ktape           !  History file number
C
C Input/output arguments
C
      integer hbufpt(ptapes)  !  History file buffer pointers
      integer nrlen           !  Actual length of data record
      integer nplen           !  Length of packed data record
      integer nflds           !  Number of fields on file
C
C Output arguments
C
      integer mfld(3)         !  Field list information: no. of levels
C                             !    and field type (A,I,X,M)
      integer npt             !  Unpacked field pointer
      integer nppt            !  Packed field pointer
      character*8 mcfld(2)    !  Field name and units
C
C---------------------------Local variables-----------------------------
C
      integer ilen            ! Length of slice of unpacked field
      integer iplen           ! Length of slice of packed field.
      integer ihacc           ! Integer form of accumulation indicator
      integer itape           ! Index for history file pointers
C
C------------------------------Externals--------------------------------
C
      external endrun         !  Abnormal termination
C
C-----------------------------------------------------------------------
C
C Convert input types
C
      ihacc = 0
      if (hacc.eq.'A') ihacc = 1
      if (hacc.eq.'M') ihacc = 2
      if (hacc.eq.'X') ihacc = 3
C
C Monthly-average files are treated like time-averaged files
C (NINAVG = A)
C
      if (hacc.eq.'Q') ihacc = 1
C
C Save pointer to current field
C
      npt = nrlen + 1
      nppt = nplen + 1
C
C Determine unpacked and packed field lengths in words
C Note, iflds=0 -> single-level field; iflds>0 -> multi-level field
C
      ilen = plon   ! length of unpacked data field
      if (iflds(1).ge.1) ilen = ilen*plev
      if (nden.eq.1) then
         iplen = ilen
      else
         iplen = 2+(plon+nden-1)/nden
         if (iflds(1).ge.1) iplen = plev*iplen
      end if
C
C Update lengths of date record: unpacked and packed
C
      nrlen = nrlen + ilen   ! unpacked data record length
      nplen = nplen + iplen  ! packed data record length
C
C Increment number of fields on current history file
C
      nflds = nflds + 1   ! accumulate fields in latitude data record
C
C Set field information
C
      mcfld(1) = fieldn(1)  ! field name
      mcfld(2) = fieldn(2)  ! units of field
C
      mfld(1)  = ihacc*10 + iflds(1)  ! no. levels and averaging flag
C
C Set mfld(2) in fldlst after all fields are defined
C
      mfld(3) = nden        ! data packing density
C
C Update history buffer pointer for all files beyond this one
C
      do itape = ktape+1,ptapes
         hbufpt(itape) = hbufpt(itape) + ilen
      end do
C
C Update history buffer lengths for monthly-average restarts.
C rg3buf(ktape) is the length in words of the portion of the history
C buffer associated with history file number ktape.
C
      rg3buf(ktape) = rg3buf(ktape) + ilen
C
      return
C
      end






 
