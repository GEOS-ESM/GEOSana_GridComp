#include <preproc.h>
      subroutine histlst 

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>

* ------------------------ code history ---------------------------
* source file:       histlst.F
* purpose:           initialize active field list for history files
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* this subroutine sets:
*    o number      of active single-level and multi-level fields: n[s,m]lflds
*    o names       of active single-level and multi-level fields: [s,m]lfnam
*    o units       of active single-level and multi-level fields: [s,m]lfuni
*    o type        of active single-level and multi-level fields: [s,m]lftyp
*    o description of active single-level and multi-level fields: [s,m]lfdes

* the field types, which are set for each active field, are:
*    o average over history interval
*    o maximum in history interval
*    o minimum in history interval
*    o instantaneous when history file written

* fields can be made inactive by setting the excl variable to the
* appropriate field name via the lsmexp namelist or in subroutine lsmctli

* field type can be overridden by setting the chntyp variable to the
* appropriate field name and new field type via the lsmexp namelist
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmhis.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i,j                  !do loop index
      integer nflds                !maximum number of fields
      integer nacti                !number of active fields
      integer ind(malflds)         !index to active fields
      logical flag(malflds)        !true if field is active

      character*8 fldnam(malflds)  !field name
      character*8 flduni(malflds)  !field units
      character*8 fldlev(malflds)  !field levels: single, multi
      character*8 fldtyp(malflds)  !field type: inst, maxi, mini, aver
      character*8 nsing            !single-level field
      character*8 nmult            !multi-level field
      character*8 xnam(malflds)    !temporary field name
      character*8 xuni(malflds)    !temporary field units
      character*8 xlev(malflds)    !temporary field levels
      character*8 xtyp(malflds)    !temporary field type
      character*40 flddes(malflds) !field description
      character*40 xdes(malflds)   !temporary field description

      integer nflds1d                 !number of fields for 1d output
      integer ind1d(malflds)          !index to 1d output fields
      character*10 fldnam1d(malflds)  !1d field name
      character*8  flduni1d(malflds)  !1d field units
      character*8  fldlev1d(malflds)  !1d field levels: single, multi
      character*8  fldtyp1d(malflds)  !1d field type: inst, maxi, mini, aver
      character*40 flddes1d(malflds)  !1d field description
      character*10 name               !temporary field name
      integer lencs                   !length of character string function
* -----------------------------------------------------------------

      nflds = 0
      naver = 'average '       !average over history interval
      nmaxi = 'maximum '       !maximum in history interval
      nmini = 'minimum '       !minimum in history interval
      ninst = 'instant '       !instantaneous when history file written
      nsing = 'sing-lev'       !single-level field
      nmult = 'mult-lev'       !multi-level field

* set default fields

      nflds = nflds + 1
      fldnam(nflds) = 'TAM     '
      flduni(nflds) = 'K       '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'aerodynamic temperature at height z0m+d'

      nflds = nflds + 1
      fldnam(nflds) = 'TAH     '
      flduni(nflds) = 'K       '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'surface temperature at height z0h+d'

      nflds = nflds + 1
      fldnam(nflds) = 'TSA     '
      flduni(nflds) = 'K       '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = '2 m air temperature'

      nflds = nflds + 1
      fldnam(nflds) = 'TV      '
      flduni(nflds) = 'K       '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'vegetation temperature'

      nflds = nflds + 1
      fldnam(nflds) = 'TG      '
      flduni(nflds) = 'K       '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'ground temperature'

      nflds = nflds + 1
      fldnam(nflds) = 'TSOI    '
      flduni(nflds) = 'K       '
      fldlev(nflds) = nmult
      fldtyp(nflds) = naver
      flddes(nflds) = 'soil temperature'

      nflds = nflds + 1
      fldnam(nflds) = 'DTV     '
      flduni(nflds) = 'K       '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'change in tv last iteration'

      nflds = nflds + 1
      fldnam(nflds) = 'DTG     '
      flduni(nflds) = 'K       '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'change in tg last iteration'

      nflds = nflds + 1
      fldnam(nflds) = 'FSA     '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'absorbed solar radiation'

      nflds = nflds + 1
      fldnam(nflds) = 'FSR     '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'reflected solar radiation'

      nflds = nflds + 1
      fldnam(nflds) = 'NDVI    '
      flduni(nflds) = '        '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'surface ndvi'

      nflds = nflds + 1
      fldnam(nflds) = 'FCTR    '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'canopy transpiration'

      nflds = nflds + 1
      fldnam(nflds) = 'FCEV    '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'canopy evaporation'

      nflds = nflds + 1
      fldnam(nflds) = 'FGEV    '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'ground evaporation'

      nflds = nflds + 1
      fldnam(nflds) = 'TAUX    '
      flduni(nflds) = 'KG/MS2  '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'zonal surface stress'

      nflds = nflds + 1
      fldnam(nflds) = 'TAUY    '
      flduni(nflds) = 'KG/MS2  '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'meridional surface stress'

      nflds = nflds + 1
      fldnam(nflds) = 'FIRA    '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'net infrared (longwave) radiation'

      nflds = nflds + 1
      fldnam(nflds) = 'FIRE    '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'emitted infrared (longwave) radiation'

      nflds = nflds + 1
      fldnam(nflds) = 'FSH     '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'sensible heat'

      nflds = nflds + 1
      fldnam(nflds) = 'FGR     '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'heat flux into soil'

      nflds = nflds + 1
      fldnam(nflds) = 'FSM     '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'snow melt heat flux'

      nflds = nflds + 1
      fldnam(nflds) = 'CN      '
      flduni(nflds) = '        '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'neutral momentum drag coefficient'

      nflds = nflds + 1
      fldnam(nflds) = 'CM      '
      flduni(nflds) = '        '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'momentum drag coefficient'

      nflds = nflds + 1
      fldnam(nflds) = 'ELAI    '
      flduni(nflds) = 'M2/M2   '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'exposed one-sided leaf area index'

      nflds = nflds + 1
      fldnam(nflds) = 'ESAI    '
      flduni(nflds) = 'M2/M2   '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'exposed one-sided stem area index'

      nflds = nflds + 1
      fldnam(nflds) = 'RSSUN   '
      flduni(nflds) = 'S/M     '
      fldlev(nflds) = nsing
      fldtyp(nflds) = nmini
      flddes(nflds) = 'sunlit leaf stomatal resistance'

      nflds = nflds + 1
      fldnam(nflds) = 'RSSHA   '
      flduni(nflds) = 'S/M     '
      fldlev(nflds) = nsing
      fldtyp(nflds) = nmini
      flddes(nflds) = 'shaded leaf stomatal resistance'

      nflds = nflds + 1
      fldnam(nflds) = 'FMICR   '
      flduni(nflds) = 'UMOL/M2S'
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'microbial respiration'

      nflds = nflds + 1
      fldnam(nflds) = 'FPSN    '
      flduni(nflds) = 'UMOL/M2S'
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'photosynthesis'

      nflds = nflds + 1
      fldnam(nflds) = 'FRMS    '
      flduni(nflds) = 'UMOL/M2S'
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'stem maintenance respiration'

      nflds = nflds + 1
      fldnam(nflds) = 'FRMR    '
      flduni(nflds) = 'UMOL/M2S'
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'root maintenance respiration'

      nflds = nflds + 1
      fldnam(nflds) = 'FRMF    '
      flduni(nflds) = 'UMOL/M2S'
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'foliage maintenance respiration'

      nflds = nflds + 1
      fldnam(nflds) = 'FRG     '
      flduni(nflds) = 'UMOL/M2S'
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'growth respiration'

      nflds = nflds + 1
      fldnam(nflds) = 'FCO2    '
      flduni(nflds) = 'UMOL/M2S'
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'net co2 flux'

      nflds = nflds + 1
      fldnam(nflds) = 'DMI     '
      flduni(nflds) = 'UG/M2S  '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'net primary production'

      nflds = nflds + 1
      fldnam(nflds) = 'H2OSOI  '
      flduni(nflds) = 'MM3/MM3 '
      fldlev(nflds) = nmult
      fldtyp(nflds) = naver
      flddes(nflds) = 'volumetric soil water'

      nflds = nflds + 1
      fldnam(nflds) = 'H2OSNO  '
      flduni(nflds) = 'MM      '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'snow depth (liquid water)'

      nflds = nflds + 1
      fldnam(nflds) = 'H2OCAN  '
      flduni(nflds) = 'MM      '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'intercepted water'

      nflds = nflds + 1
      fldnam(nflds) = 'BEVAP   '
      flduni(nflds) = '        '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'ground evaporation beta factor'

      nflds = nflds + 1
      fldnam(nflds) = 'BTRAN   '
      flduni(nflds) = '        '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'transpiration beta factor'

      nflds = nflds + 1
      fldnam(nflds) = 'QINFL   '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'infiltration'

      nflds = nflds + 1
      fldnam(nflds) = 'QOVER   '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'surface runoff'

      nflds = nflds + 1
      fldnam(nflds) = 'QDRAI   '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'sub-surface drainage'

      nflds = nflds + 1
      fldnam(nflds) = 'QINTR   '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'interception'

      nflds = nflds + 1
      fldnam(nflds) = 'QDRIP   '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'throughfall'

      nflds = nflds + 1
      fldnam(nflds) = 'QMELT   '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'snow melt'

      nflds = nflds + 1
      fldnam(nflds) = 'QSOIL   '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'ground evaporation'

      nflds = nflds + 1
      fldnam(nflds) = 'QVEGE   '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'canopy evaporation'

      nflds = nflds + 1
      fldnam(nflds) = 'QVEGT   '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'canopy transpiration'

      nflds = nflds + 1
      fldnam(nflds) = 'QPRECC  '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'convective precipitation'

      nflds = nflds + 1
      fldnam(nflds) = 'QPRECL  '
      flduni(nflds) = 'MM/S    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'large-scale precipitation'

      nflds = nflds + 1
      fldnam(nflds) = 'ERRSOI  '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'soil/lake energy conservation error'

      nflds = nflds + 1
      fldnam(nflds) = 'ERRSEB  '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'surface energy conservation error'

      nflds = nflds + 1
      fldnam(nflds) = 'ERRSOL  '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'solar radiation conservation error'

      nflds = nflds + 1
      fldnam(nflds) = 'ERRH2O  '
      flduni(nflds) = 'MM      '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'total water conservation error'

      nflds = nflds + 1
      fldnam(nflds) = 'TGCM    '
      flduni(nflds) = 'K       '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'atmospheric air temperature'

      nflds = nflds + 1
      fldnam(nflds) = 'VMAGGCM '
      flduni(nflds) = 'M/S     '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'atmospheric wind velocity magnitude'

      nflds = nflds + 1
      fldnam(nflds) = 'PGCM    '
      flduni(nflds) = 'PA      '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'atmospheric pressure'

      nflds = nflds + 1
      fldnam(nflds) = 'QGCM    '
      flduni(nflds) = 'KG/KG   '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'atmospheric specific humidity'

      nflds = nflds + 1
      fldnam(nflds) = 'ZGCM    '
      flduni(nflds) = 'M       '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'atmospheric reference height'

      nflds = nflds + 1
      fldnam(nflds) = 'L_DOWN  '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'incident longwave radiation'

      nflds = nflds + 1
      fldnam(nflds) = 'S_DOWN  '
      flduni(nflds) = 'W/M2    '
      fldlev(nflds) = nsing
      fldtyp(nflds) = naver
      flddes(nflds) = 'incident solar radiation'

* override fldtyp variable based on namelist input

      do i = 1, nflds
         do j = 1, nflds
            if (fldnam(i) .eq. chntyp(1,j)) then
               if (chntyp(2,j).eq.naver .or. chntyp(2,j).eq.nmaxi .or.
     &             chntyp(2,j).eq.nmini .or. chntyp(2,j).eq.ninst) then
                  fldtyp(i) = chntyp(2,j)
               else
                  write (6,*) 'histlst error: attempting to change',
     &               ' field type for field = ',fldnam(i) 
                  write (6,*) 'to inaccurate type = ',chntyp(2,j),
     &               '. valid types are: ',naver,nmaxi,nmini,ninst
                  call stoprun
               end if
            end if
         end do
      end do

* flag active fields

      do i = 1, nflds
         flag(i) = .true.
         do j = 1, nflds
            if (fldnam(i) .eq. excl(j)) then
               flag(i) = .false.
            end if
         end do
      end do

* number of active fields (nacti <= nflds) and pointer (1 to nflds)

      nacti = 0
      do i = 1, nflds
         if (flag(i)) then
            nacti = nacti + 1
            ind(nacti) = i
         end if
      end do

* re-order fields from 1 -> nflds to 1 -> nacti

      nflds = nacti
      if (nflds .gt. malflds) then
         write (6,*) 'histlst error: number fields > maximum allowed'
         call stoprun
      end if

      do i = 1, nflds
         xnam(i) = fldnam(ind(i))
         xuni(i) = flduni(ind(i))
         xlev(i) = fldlev(ind(i))
         xtyp(i) = fldtyp(ind(i))
         xdes(i) = flddes(ind(i))
      end do

      do i = 1, nflds
         fldnam(i) = xnam(i) 
         flduni(i) = xuni(i) 
         fldlev(i) = xlev(i) 
         fldtyp(i) = xtyp(i) 
         flddes(i) = xdes(i)   
      end do

* separate single-level and multi-level fields

      nslflds = 0
      nmlflds = 0
      do i = 1, nflds
         if (fldlev(i) .eq. nsing) then
            nslflds = nslflds + 1
            slfnam(nslflds) = fldnam(i)
            slfuni(nslflds) = flduni(i)
            slftyp(nslflds) = fldtyp(i)
            slfdes(nslflds) = flddes(i)
         else if (fldlev(i) .eq. nmult) then
            nmlflds = nmlflds + 1
            mlfnam(nmlflds) = fldnam(i)
            mlfuni(nmlflds) = flduni(i)
            mlftyp(nmlflds) = fldtyp(i)
            mlfdes(nmlflds) = flddes(i)
         end if
      end do

      if (nslflds .gt. mslflds) then
         write (6,*) 'histlst error: number single-level fields',
     &               ' > maximum allowed'
         call stoprun
      end if
      if (nmlflds .gt. mmlflds) then
         write (6,*) 'histlst error: number multi-level fields',
     &               ' > maximum allowed'
         call stoprun
      end if

* echo active fields 
      if (masterproc) then
         write (6,*) 'active lsm history fields: single-level'
         write (6,1002) 
         do i = 1, nslflds
            write (6,1003) i,slfnam(i),slfuni(i),slftyp(i),slfdes(i)
         end do
         write (6,*)
         
         write (6,*) 'active lsm history fields: multi-level'
         write (6,1002) 
         do i = 1, nmlflds
            write (6,1003) i,mlfnam(i),mlfuni(i),mlftyp(i),mlfdes(i)
         end do
         write (6,*)
      end if
* flag active fields which will have 1d output

      nflds1d = 0
      do i = 1, nflds
         do j = 1, nflds
            if (flds1d(i) .eq. fldnam(j) .and. flag(j)) then
               nflds1d = nflds1d + 1
               ind1d(nflds1d) = j
            endif
         end do
      end do

      do i = 1, nflds1d
         name = '          '
         name = fldnam(ind1d(i))
         fldnam1d(i) = name(1:lencs(name)) // '1D'
         flduni1d(i) = flduni(ind1d(i))
         fldlev1d(i) = fldlev(ind1d(i))
         fldtyp1d(i) = fldtyp(ind1d(i))
         flddes1d(i) = flddes(ind1d(i))
      end do

* separate 1d single-level and multi-level fields

      nsl1dflds = 0
      nml1dflds = 0
      do i = 1, nflds1d
         if (fldlev1d(i) .eq. nsing) then
            nsl1dflds = nsl1dflds + 1
            sl1dfnam(nsl1dflds) = fldnam1d(i)
            sl1dfuni(nsl1dflds) = flduni1d(i)
            sl1dftyp(nsl1dflds) = fldtyp1d(i)
            sl1dfdes(nsl1dflds) = flddes1d(i)
         else if (fldlev1d(i) .eq. nmult) then
            nml1dflds = nml1dflds + 1
            ml1dfnam(nml1dflds) = fldnam1d(i)
            ml1dfuni(nml1dflds) = flduni1d(i)
            ml1dftyp(nml1dflds) = fldtyp1d(i)
            ml1dfdes(nml1dflds) = flddes1d(i)
         end if
      end do

* echo 1d fields 

      if (masterproc) then
         if (nsl1dflds .gt. 0) then
            write (6,*) '1d lsm history fields: single-level'
            write (6,1002) 
            do i = 1, nsl1dflds
               write (6,1003)i,sl1dfnam(i),sl1dfuni(i),sl1dftyp(i),
     $              sl1dfdes(i)
            end do
            write (6,*)
         endif

         if (nml1dflds .gt. 0) then
            write (6,*) '1d lsm history fields: multi-level'
            write (6,1002) 
            do i = 1, nml1dflds
               write (6,1003)i,ml1dfnam(i),ml1dfuni(i),ml1dftyp(i),
     $              ml1dfdes(i)
            end do
            write (6,*)
         endif
      endif

1002  format(' no',' name    ',' units   ',' type    ',' description')
1003  format((1x,i2),(1x,a8),(1x,a8),(1x,a8),(1x,a40))

      return
      end
 
