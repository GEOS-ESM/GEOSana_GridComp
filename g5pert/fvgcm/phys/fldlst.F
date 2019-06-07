#include <params.h>
      subroutine fldlst(lhbuf)
C-----------------------------------------------------------------------
C
C Initialize each history file field, determine the size of the history
C buffer and organize its contents.  This routine combines the 
C information from the Master Field List and namelist variables 
C PRIMARY, EXCLUDE, and AUXF, to build the history buffer. The length 
C of the history buffer for a latitude slice (lhbuf) is computed, as well
C the hash field (islocc) used for sorting the fields in the history 
C files.
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath 
C Standardized:      L. Bath, June 1992
C                    T. Acker, Feb 1996
C Reviewed:          B. Boville, P. Rasch, August 1992
C Reviewed:          B. Boville, P. Rasch, April 1996
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
#include <comhed.h>
C-----------------------------------------------------------------------
#include <comhst.h>
C------------------------------Arguments--------------------------------
C
C Output argument
C
      integer lhbuf            ! Length of the hist. buffer, hbuf
C
C---------------------------Local variables-----------------------------
C
      character*1 ninst    !  Code for an instantaneous field value
      character*1 avg      !  Temporary for averaging flag
      character*8 fldtmp   !  Temporary for field name
c
      integer it           !  history file/tape index
      integer if           !  Index into master field list for each field
      integer ifex         !  Index into exclude fields
      integer ifin         !  Index into include fields (primary)
      integer i            !  Index
c
      integer kfld         !  Index into mflds,mcflds etc. for each field
      integer kkfld        !  Index like kfld  but for aux. tapes
      integer kt           !  Index over mtapes
      integer ipos         !  Index into 'ichar sum' array
c
      integer iascii       !  Index
      integer imlist       !  Index used for auxiliary history file field
      integer iaux         !  Index used for auxiliary history file field
      integer itape        !  Index to history file/tape
C
C------------------------------Externals--------------------------------
C
      external inifld      ! Load field list information with a field
      external endrun      ! Abnormal termination
      external ljust       ! Left justify a character string
C
C-----------------------------------------------------------------------
C
C Set ninst to represent an instantaneous field value.
C
      ninst = 'I'
C Initialize counters
C
      if (mtapes.eq.0) then    ! True if an initialization or branch run
C
C Initialize "number of fields per ichar sum" array to zero
C
         do it=1,ptapes
            do iascii=0,pichsum
               nfldsc(iascii,it) = 0
            end do
         end do
C
C Set current no. hist. files to 1 (primary)
C
         mtapes = 1    
C
C The following initializations are for the primary history file
C        
         nflds(mtapes) = 0     ! init. no. fields per data record
         nrlen(mtapes) = 0     ! init. length of data record
         nplen(mtapes) = 0     ! init. length of packed data record
         kfld = 1              ! index indicating fld. no. in data rec.
C
C Primary history file. Call inifld for all active fields. The exclude 
C list is used to specify which default fields on the Master Field List 
C will not be written to the primary history file.  
C
         do 10 if=1,pflds      ! pflds = no. flds. in Master Field List
C
C Search the exclude array. Fields in the exclude array are left out of
C the history buffer.
C The exclude array derives from namelist variable EXCLUDE. 
C
            do ifex = 1,nexcl 
               call ljust(exclude(ifex),fldtmp)
               if (fieldn(1,if).eq.fldtmp) go to 10
            end do
C
C Call inifld and include in the history buffer all fields in the  
C primary array. The primary array derives from namelist 
C variable PRIMARY.
C
            do ifin = 1,nincl
               call ljust(primary(ifin),fldtmp)
               if (fieldn(1,if).eq.fldtmp) then
C
C Set the field averaging flag to A, I, M, or X.  
C PHIS and ORO are always instantaneous (I).
C
                  if (fieldn(1,if).eq.'PHIS    '.or.
     $                fieldn(1,if).eq.'ORO     ') then
                     avg = ninst  
                  else
                     avg = ninavg(mtapes)  
                  end if
C
C Call inifld and load field into the history buffer. Inifld specifies
C field information list arrays mflds and mcflds, as well as pointers
C to fields and files within the history buffer.
C
                  call inifld(fieldn(1,if),iflds(1,if)  ,avg     ,
     $                          hbufpt  ,ndens(1),
     $                        mtapes  ,nrlen(mtapes) ,nplen(mtapes),
     $                          nflds(mtapes) ,mflds(1,kfld) ,
     $                        mcflds(1,kfld),nupnt(kfld,mtapes),
     $                          npnt(kfld,mtapes))
C
C Add one to the count of active field names starting with a given 
C ichar sum in the primary history file (nfldsc stores the no. of 
C fields in history file with ichar sum equal to ipos.)
C
                  ipos = ichar(mcflds(1,kfld)(1:1)) +
     $                   ichar(mcflds(1,kfld)(2:2)) +
     $                   ichar(mcflds(1,kfld)(3:3)) +
     $                   ichar(mcflds(1,kfld)(4:4))
                  nfldsc(ipos,mtapes) = nfldsc(ipos,mtapes) + 1
C
C Increment the number of fields on the primary history file.
C
                  kfld = kfld + 1
                  go to 10  ! process next field
               end if
            end do
C
C Search Master Field List for default active fields, and include
C them in the history buffer.
C
            if (iflds(3,if).eq.1) then
C
C Set the field averaging flag to A, I, M, or X.  
C
               if (fieldn(1,if).eq.'PHIS    '.or.
     $             fieldn(1,if).eq.'ORO     ') then
                  avg = ninst
               else
                  avg = ninavg(mtapes)
               end if
C
C Call inifld and load field into the history buffer.
C
               call inifld(fieldn(1,if),iflds(1,if)  ,avg     ,
     $                       hbufpt  ,ndens(1),
     $                     mtapes  ,nrlen(mtapes) ,nplen(mtapes),
     $                       nflds(mtapes) ,mflds(1,kfld),
     $                     mcflds(1,kfld),nupnt(kfld,mtapes),
     $                       npnt(kfld,mtapes))
C
C Add one to the count of active field names starting with a given 
C ichar sum in the primary history file.
C
               ipos = ichar(mcflds(1,kfld)(1:1)) +
     $                ichar(mcflds(1,kfld)(2:2)) +
     $                ichar(mcflds(1,kfld)(3:3)) +
     $                ichar(mcflds(1,kfld)(4:4))
               nfldsc(ipos,mtapes) = nfldsc(ipos,mtapes) + 1
C
C Increment the number of fields on the primary history file.
C
               kfld = kfld + 1
            end if
   10    continue
      else           ! do this branch if a restart or regen. run
C
C Continuation run, primary history file already declared
C
         kfld = 1
         do i=1,mtapes
C
C Set kfld to no. of fields on all previously declared hist. files
C
            kfld = kfld + nflds(i)
         end do
      end if
C
C Check for declaration of auxiliary history files. Note, new aux. hist.
C files can be declared on restart and regen. runs.
C
      do 20 itape=mtapes,ptapes-1 ! begin loop to define aux. hist. files
         nflds(itape+1) = 0   ! init. no. fields per aux file data record
         nrlen(itape+1) = 0   ! init. actual len. of aux file data record
         nplen(itape+1) = 0   ! init. len. of packed aux file data record
C
C If 1st "aux" is not blank, the user has specified an auxilliary file
C for the current "itape".  Increment mtapes accordingly (mtapes is the 
C no. of primary plus auxiliary history files requested).
C
         if (aux(1,itape).ne.'        ') then
            mtapes = mtapes + 1
C
C Set hist. buffer pointer for aux. hist. file to end of the fields in  
C the hist. buffer from the previous history file 
C
            hbufpt(mtapes) = hbufpt(mtapes-1) + nrlen(mtapes-1)
         endif
C
C Make sure all requested auxiliary fields are on the Master Field List 
C
         do iaux=1,pflds
            if (aux(iaux,itape).eq.'        ') go to 200
            do imlist=1,pflds
               if (aux(iaux,itape).eq.fieldn(1,imlist)) go to 100
            end do
            write(6,*) 'FLDLST:cant find auxilliary field ',
     $                 aux(iaux,itape)
  100       continue
         end do
  200    continue
C
C Set counter kkfld as no. of fields on current aux. hist. file
C
         kkfld = 1
         do 15 imlist=1,pflds
            do iaux=1,pflds
C
C Once we hit a blank "aux" entry, all auxiliary fields for this file
C have been specified
C
               if (aux(iaux,itape).eq.'        ') go to 15
               if (aux(iaux,itape).eq.fieldn(1,imlist)) then
C
C Set the field averaging flag to A, I, M, or X.  
C
                   if (aux(iaux,itape).eq.'PHIS    '.or.
     $                aux(iaux,itape).eq.'ORO     ') then
                     avg = ninst
                  else
                     avg = ninavg(mtapes)
                  end if
C
C Call inifld and load field into the history buffer.
C
                  call inifld(fieldn(1,imlist),iflds(1,imlist),avg,
     $                          hbufpt  ,ndens(itape+1),
     $                        mtapes  ,nrlen(mtapes) ,nplen(mtapes),
     $                          nflds(mtapes) ,mflds(1,kfld) ,
     $                        mcflds(1,kfld),nupnt(kkfld,mtapes),
     $                          npnt(kkfld,mtapes))
C
C Add one to the count of active field names starting with a given 
C ichar sum in the current auxilary history file.
C
                  ipos = ichar(mcflds(1,kfld)(1:1)) +
     $                   ichar(mcflds(1,kfld)(2:2)) +
     $                   ichar(mcflds(1,kfld)(3:3)) +
     $                   ichar(mcflds(1,kfld)(4:4))
                  nfldsc(ipos,mtapes) = nfldsc(ipos,mtapes) + 1
C
C Increment the number of fields in the history buffer.
C
                  kfld = kfld + 1
C
C Increment the number of fields in current aux. hist file.
C
                  kkfld = kkfld + 1
                  go to 15
               end if
            end do
   15    continue
   20 continue        ! end loop to define aux. hist. files
C
C Define islocc array.  Each element of the islocc array contains the
C no. of fields in the corresponding history file that has an "ichar
C sum" less than the value of the index.  
C
C Also, complete loading mflds array and compute the length of 
C the history buffer.
C
      lhbuf = 0
      kfld = 1
      do kt=1,mtapes
         islocc(0,kt) = 0   ! no valid field names have ichar sum of 0
         do i=1,pichsum
            islocc(i,kt) = islocc(i-1,kt) + nfldsc(i-1,kt)
         end do
         do i=1,nflds(kt)
C
C Load unpacked field pointers into mflds array. Add 2 to account for
C presence of lat. index and long. count at beginning of each latitude
C data record.
C
            mflds(2,kfld) = npnt(i,kt) + 2
C
C Compute length of history buffer
C
            if (mod(mflds(1,kfld),10).eq.0) then
               lhbuf = lhbuf + plon
            else
               lhbuf = lhbuf + plon*plev
            end if
C
C Accumulate no. of fields in all hist. files (i.e. in the hist. buffer)
C
            kfld = kfld + 1  
         end do
      end do
C
C Pad length of history tape buffer for unblocked I/O
C
      lhbuf = (lhbuf/512 + 1)*512
      return
C
      end
 
