#include <preproc.h>
      subroutine histwrt (nstep,beglatspmd)

#include <vartyp.h>
#include <lsmpar.h>
#include <basinpar.h>

* ------------------------ code history ---------------------------
* source file:       histwrt.F
* purpose:           write netcdf history file 
* date last revised: April 1997
* author:            Mariana Vertenstein
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ input variables ------------------------
      integer nstep         !current time index
      integer beglatspmd    !beg lat on proc, degenerates to 1 for nonspmd
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmhis.h>
#include <lsmhisdyn.h>
#include <lsmtc.h>
#include <lsmtcdyn.h>
#include <lsmctl.h>
#include <ncd.h>
#include <lsmspmd.h>
#include <basin.h>
* -----------------------------------------------------------------

* ------------------------ netcdf declarations---------------------
#include <netcdf.inc>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i,j,n                    !indices
      integer start1d(1)               !1-d start index
      integer count1d(1)               !1-d count value
      integer start2d(2)               !2-d start indices
      integer count2d(2)               !2-d count values
      integer start3d(3)               !3-d start indices
      integer count3d(3)               !3-d count values
      integer start4d(4)               !4-d start indices
      integer count4d(4)               !4-d count values

      real simtim                      !current simulation time
      real slxyfval(lsmlon,lsmlat)     !single level xy field values
      real mlxyfval(lsmlon,lsmlat,msl) !multi  level xy field values

      integer lencs                    !length of character string function
      real longrd(lsmlon)              !grid longitudes
      real latgrd(lsmlat)              !grid latitudes

#if ( ! defined CRAY )
      pointer ( plonsub, lonsub )     
      pointer ( platsub, latsub )   
#endif
      real lonsub(kpt)                 !subgrid longitudes
      real latsub(kpt)                 !subgrid latitudes
      real pi                          !3.14156

      real lon2d(lsmlon,lsmlat)        ! reduced grid longitude array
      real lat2d(lsmlon,lsmlat)        ! reduced grid latitude array
      integer mask2d(lsmlon,lsmlat)    ! reduced grid mask array
* -----------------------------------------------------------------

c
c Gather up the single and multifield array values from all processors
c
* -----------------------------------------------------------------
* write out time-independent grid related variables
* -----------------------------------------------------------------

      if (masterproc) then
      if (ntim .eq. 1) then
            
* write out lsm grid maximum longitude array (usually at equator)
* and latitude array
        
         pi = 4.*atan(1.)
         do i=1,lsmlon
            longrd(i) = longxy(i,1)*180./pi
            if (longrd(i) .lt. 0.) longrd(i) = longrd(i) + 360.
         end do

         start1d(1)=1
         count1d(1)=lsmlon
         call ncvpt (ncid,var_lon_id,start1d,count1d,longrd,rcode)

         do j=1,lsmlat
            latgrd(j) = latixy(1,j)*180./pi
         end do

         start1d(1)=1
         count1d(1)=lsmlat
         call ncvpt (ncid,var_lat_id,start1d,count1d,latgrd,rcode)

* write out lsm reduced grid longitude, latitude and mask arrays

         do j=1,lsmlat
            do i=1,numlon(j)
               lon2d(i,j)  = longxy(i,j)*180./pi
               if (lon2d(i,j) .lt. 0) lon2d(i,j) = lon2d(i,j) + 360.
               lat2d(i,j)  = latixy(i,j)*180./pi
               mask2d(i,j) = 1
            end do
            do i=numlon(j)+1,lsmlon
               lon2d(i,j)  = spval
               lat2d(i,j)  = spval
               mask2d(i,j) = 0
            end do
         end do

         start2d(1)=1
         count2d(1)=lsmlon
         start2d(2)=1
         count2d(2)=lsmlat
         call ncvpt (ncid,var_lon2d_id ,start2d,count2d,lon2d ,rcode)
         call ncvpt (ncid,var_lat2d_id ,start2d,count2d,lat2d ,rcode)
         call ncvpt (ncid,var_mask2d_id,start2d,count2d,mask2d,rcode)

* write out sub-grid longitudes and latitudes

         do i=1,kpt
            lonsub(i) = long(i)*180./pi
            latsub(i) = lati(i)*180./pi
         end do

         start1d(1)=1
         count1d(1)=kpt
         call ncvpt (ncid,var_sublon_id,start1d,count1d,lonsub,rcode)
         call ncvpt (ncid,var_sublat_id,start1d,count1d,latsub,rcode)

* write out land points longitude and latitude indices

         start1d(1)=1
         count1d(1)=lpt
         call ncvpt (ncid,var_ixy_id,start1d,count1d,ixy,rcode)
         call ncvpt (ncid,var_jxy_id,start1d,count1d,jxy,rcode)

* write out sub-grid land and weight indices

         start2d(1)=1
         count2d(1)=lpt
         start2d(2)=1
         count2d(2)=msub
         call ncvpt (ncid,var_subkvc_id,start2d,count2d,kvec ,rcode)
         call ncvpt (ncid,var_subwgt_id,start2d,count2d,wsg2g,rcode)

* write out land surface type

         start2d(1)=1
         count2d(1)=lsmlon
         start2d(2)=1
         count2d(2)=lsmlat
         call ncvpt (ncid,var_grdsrf_id,start2d,count2d,surf2d,rcode)

* write out time-invariant model timeing variables

         start1d(1)=1
         count1d(1)=1 
         call ncvpt (ncid,mdbase_id,start1d,count1d,mdbase,rcode)
         call ncvpt (ncid,msbase_id,start1d,count1d,msbase,rcode)
         call ncvpt (ncid,mbdate_id,start1d,count1d,mbdate,rcode)
         call ncvpt (ncid,mbsec_id ,start1d,count1d,mbsec ,rcode)
         call ncvpt (ncid,dtlsm_id ,start1d,count1d,dtlsm ,rcode)
         call ncvpt (ncid,nhtfrq_id,start1d,count1d,nhtfrq,rcode)

* write out input data file

         start1d(1)=1
         count1d(1)=lencs(finidat) 
         call ncvptc(ncid,inifil_id,start1d,count1d,finidat,
     $        lencs(finidat) ,rcode)

* note the following will be implemented when reduced grids are implemented
* in the lsm

!         start1d(1)=1
!         ncount1d(1)=lencs(flondat) 
!         call ncvptc(ncid,lonfil_id,start1d,count1d,flnddat,
!     $        lencs(flondat) ,rcode)

      endif

* -----------------------------------------------------------------
* write out time dependent variables
* -----------------------------------------------------------------

* If restart file, get variable id's

      if (nsrest.eq.1 .and. ntim.gt.1 .and. ncgetvid) then

         var_tim_id(1) = ncvid(ncid,'date'    ,rcode)
         var_tim_id(2) = ncvid(ncid,'datesec' ,rcode)
         var_tim_id(3) = ncvid(ncid,'time'    ,rcode)
         var_tim_id(4) = ncvid(ncid,'timestep',rcode)
         var_tim_id(5) = ncvid(ncid,'timecom' ,rcode)

         do i=1,nslflds
            slfld_id(i) = ncvid(ncid,slfnam(i),rcode)
         end do
         do i=1,nmlflds
            mlfld_id(i) = ncvid(ncid,mlfnam(i),rcode)
         end do
         do i=1,nsl1dflds
            sl1dfld_id(i) = ncvid(ncid,sl1dfnam(i),rcode)
         end do
         do i=1,nml1dflds
            ml1dfld_id(i) = ncvid(ncid,ml1dfnam(i),rcode)
         end do

         var_bas_id = ncvid(ncid, 'BASIN_RUNOFF',rcode)

         ncgetvid = .false.

      endif

* write out time information

      start1d(1) = ntim
      count1d(1) = 1
      simtim = float(mdcur) + float(mscur)/86400.

      call ncvpt (ncid, var_tim_id(1),start1d,count1d,mcdate,rcode)
      call ncvpt (ncid, var_tim_id(2),start1d,count1d,mcsec ,rcode)
      call ncvpt (ncid, var_tim_id(3),start1d,count1d,simtim,rcode)
      call ncvpt (ncid, var_tim_id(4),start1d,count1d,nstep ,rcode)

      start2d(1) = 1
      count2d(1) = lencs(timcom)
      start2d(2) = ntim
      count2d(2) = 1
      call ncvptc(ncid, var_tim_id(5),start2d,count2d,timcom,
     $     lencs(timcom),rcode)

* write out single-level xy fields 

      start3d(1) = 1
      count3d(1) = lsmlon
      start3d(2) = 1
      count3d(2) = lsmlat
      start3d(3) = ntim
      count3d(3) = 1

      do i=1,nslflds
         call slf2xy(slfval(1,i), slxyfval)
         call ncvpt(ncid,slfld_id(i),start3d,count3d,slxyfval,rcode)
      end do

* write out multi-level xy fields 

      start4d(1) = 1
      count4d(1) = lsmlon
      start4d(2) = 1
      count4d(2) = lsmlat
      start4d(3) = 1
      count4d(3) = msl
      start4d(4) = ntim
      count4d(4) = 1

      do i=1,nmlflds
         call mlf2xy(mlfval(1,1,i), mlxyfval)
         call ncvpt(ncid,mlfld_id(i),start4d,count4d,mlxyfval,rcode)
      end do

* write out single-level fields as 1-d vectors

      start2d(1) = 1
      count2d(1) = kpt
      start2d(2) = ntim
      count2d(2) = 1

      do i=1,nsl1dflds
         n = 0
         do j = 1, nslflds
            if (flds1d(i) .eq. slfnam(j)) then
               n = j
               go to 10
            endif
         end do
   10    continue
         if (n .eq. 0) then
            write(6,*)'(HISTWRT): could not match single level 1d flds'
            write(6,*)'(HISTWRT): 1d field name : ',flds1d(i)
            stop
         endif
         call ncvpt(ncid,sl1dfld_id(i),start2d,count2d,slfval(1,n),
     $        rcode)
      end do

* write out multi-level fields as 1-d vectors

      start3d(1) = 1
      count3d(1) = kpt
      start3d(2) = 1
      count3d(2) = msl
      start3d(3) = ntim
      count3d(3) = 1

      do i=1,nml1dflds
         n = 0
         do j = 1, nmlflds
            if (flds1d(i) .eq. mlfnam(j)) then
               n = j
               go to 20
            endif
         end do
   20    continue
         if (n .eq. 0) then
            write(6,*)'(HISTWRT): could not match multi level 1d flds'
            write(6,*)'(HISTWRT): 1d field name : ',flds1d(i)
            stop
         endif
         call ncvpt(ncid,ml1dfld_id(i),start3d,count3d,mlfval(1,1,n),
     $        rcode)
      end do

* write out basin drainage vector

      if (ncbasin) then
         start2d(1) = 1
         count2d(1) = ndrn
         start2d(2) = ntim
         count2d(2) = 1
         call ncvpt(ncid,var_bas_id,start2d,count2d,bashist,rcode)
      endif
      end if     !masterproc

      return
      end
 
