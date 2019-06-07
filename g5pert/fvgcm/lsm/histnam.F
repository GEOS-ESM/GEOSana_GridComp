#include <preproc.h>
      subroutine histnam 

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>

* ------------------------ code story ---------------------------
* source file:       histnam.F
* purpose:           generate lsm history or file names 
* date last revised: April 1997
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* generate file names by appending a character extension
* to "lsmh_". then append this to the local directory name [locpnh] 
* to create the full file name [locfnh]. the character extension is
* five characters based on the file number [nfil], e.g. 00001, or 
* seven characters based on the month [nmo] and year [nyr], e.g., 0001-01.
* in addition, history files will be given an additional 3 character 
* (.nc) extension
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <lsmhis.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      character*5 prefix  !prefix for lsm history file names
      character*16 filext !unique file extension: 
                          !e.g., lsmh_00001.nc, lsmh_0001-01.nc
      character*3 ncdext  !netcdf extension for history files
      integer lencs       !integer function gets length of character string
* -----------------------------------------------------------------

      if (ninavg.eq.'q' .or. ninavg.eq.'Q') then

* unique "nyr-nmon" file extension (7 characters). 
* nyr = 4 characters. "-" = 1 character. nmon = 2 characters.

         if (nyr .gt. 9999) then              
            if (masterproc)
     &           write (6,*) 'histnam error: nyr = ',nyr,' is > 9999'
            call stoprun
         end if

      else

* unique "nfil" file extension 

         if (nfil .gt. 99999) then              
            if (masterproc)
     &        write (6,*) 'histnam error: nfil = ',nfil,' is > 99999'
            call stoprun
         end if

      end if

* unique character extension: "prefix"+"nfil"+"ncdext" or 
* "prefix"+"nyr-nmon"+"ncdext"  

      prefix = 'lsmh_'
      ncdext = '.nc' 
      if (ninavg.eq.'q' .or. ninavg.eq.'Q') then
         write (filext,'(a5,i4.4,"-",i2.2,a3)') prefix,nyr,nmon,ncdext       
      else
         write (filext,'(a5,i5.5,a3," ")') prefix,nfil,ncdext       
      end if

* add extension to local directory to get full file name

      locfnh = locpnh(1:lencs(locpnh)) // '/'
     $     // filext(1:lencs(filext)) 

      return
      end
 
