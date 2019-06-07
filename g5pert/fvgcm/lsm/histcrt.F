#include <preproc.h>
      subroutine histcrt

#include <vartyp.h>
#include <lsmpar.h>
#include <basinpar.h>

* ------------------------ code history ---------------------------
* source file:       histcrtv.F
* purpose:           creates new netcdf history tape 
* date last revised: July 1996 - lsm version 1
* author:            Mariana Vertenstein
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmhis.h>
#include <lsmtvdyn.h>
#include <ncd.h>
#include <basin.h>
* -----------------------------------------------------------------

* ------------------------ netcdf declarations---------------------
#include <netcdf.inc>
* -----------------------------------------------------------------
* ------------------------ local variables ------------------------

      integer i                 ! loop index
      integer omode

      integer dim_tim_id        ! id for time dimension
      integer dim_lon_id        ! id for grid longitude 
      integer dim_lat_id        ! id for grid latitude 
      integer dim_lev_id        ! id for soil layer dimension
      integer dim_lpt_id        ! id for lnd points, not incld subgrid
      integer dim_sub_id        ! id for max number of subgrid points
      integer dim_kpt_id        ! id for lnd points, including subgrid
      integer dim_str_id        ! id for character string variables
      
      integer dim2_id(2)        ! dim id for 2-d variables
      integer dim3_id(3)        ! dim id for 3-d variables
      integer dim4_id(4)        ! dim id for 4-d variables

      character*256 name        ! name of attribute
      character*256 unit        ! units of attribute
      character*8   curdate     ! current date
      character*8   curtime     ! current time 

      integer lencs             ! length of character string function

C SJL
C	return   ! testing !!!
* -----------------------------------------------------------------

* Create a new netcdf file with a unique history file name 
* from nfil (nmon and nyr if using monthly average option).

      ncid = nccre(locfnh(1:lencs(locfnh)), NCCLOB, rcode)
      omode = ncsfil(ncid, NCNOFILL, rcode) 
      
* Create global attributes 
      do i=1,256
         caseid(i:i) = ' '
      end do
C SJL
      call getenv ('CASE',caseid)
c gsl caseid = 'lsmrun'
      if (caseid(1:1) .eq. ' ') then
         write (6,*) 'histcrt: case not defined'
         call stoprun
      end if
      if (ctitle(1:1) .eq. ' ') then
         name = caseid(1:lencs(caseid))
      else
         name = caseid(1:lencs(caseid))//' : '//ctitle(1:lencs(ctitle))
      endif
c     print *, 'lencs(NAME)=',lencs(name)
c     print *, 'NAME=', name
c     do i=1,256
c        print *, i, '(', name(i:i), ')', ichar(name(i:i))
c     end do
      call ncaptc(ncid,ncglobal,'title'  ,NCCHAR,lencs(name),
     $     name,rcode)

      call timedate(curdate,curtime)
      write(name,1000) curdate, curtime
      call ncaptc(ncID,ncglobal,'history',NCCHAR,lencs(name),
     $     name,rcode)

      name = 'Land Surface Model: lsm1.1'
      call ncaptc(ncid, ncglobal,'Source'     ,NCCHAR,lencs(name),
     $     name,rcode)
     
* Define dimsions 
      
      dim_lon_id = ncddef(ncid,'longitude'       ,lsmlon ,rcode)

      dim_lat_id = ncddef(ncid,'latitude'        ,lsmlat ,rcode)

      dim_kpt_id = ncddef(ncid,'1d_subgrid_index',kpt    ,rcode)

      dim_sub_id = ncddef(ncid,'2d_subgrid_index',msub   ,rcode)

      dim_lpt_id = ncddef(ncid,'2d_land_index'   ,lpt    ,rcode)

      dim_lev_id = ncddef(ncid,'soil_layer_index',msl    ,rcode)

      dim_tim_id = ncddef(ncid,'time'            ,NCUNLIM,rcode)

      dim_str_id = ncddef(ncid,'strlen'          ,256    ,rcode) 

* Define coordinate variables - make them the same name as
* the corresponding netcdf dimensions
 
      var_lon_id = ncvdef(ncid,'longitude',NCPREC,1,dim_lon_id,rcode)

      var_lat_id = ncvdef(ncid,'latitude' ,NCPREC,1,dim_lat_id,rcode)

* Create attributes for coordinate variables

      name = 'longitude'
      unit = 'degrees_east'
      call ncaptc(ncid,var_lon_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_lon_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'latitude'
      unit = 'degrees_north'
      call ncaptc(ncid,var_lat_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_lat_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
* Define reduced grid 2d longitude, latitude and mask grids

      dim2_id(1) = dim_lon_id
      dim2_id(2) = dim_lat_id
      var_lon2d_id  = ncvdef(ncid,'lon2d' ,NCPREC,2,dim2_id,rcode)
      var_lat2d_id  = ncvdef(ncid,'lat2d' ,NCPREC,2,dim2_id,rcode)
      var_mask2d_id = ncvdef(ncid,'mask2d',NCLONG,2,dim2_id,rcode)

* Create attributes for reduced grid arrays

      name = 'reduced grid 2d longitude array'
      unit = 'degrees_east'
      call ncaptc(ncid,var_lon2d_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_lon2d_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'reduced grid 2d latitude array'
      unit = 'degrees_north'
      call ncaptc(ncid,var_lat2d_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_lat2d_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)

      name = 'reduced grid 2d mask array'
      unit = 'unitless'
      call ncaptc(ncid,var_mask2d_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_mask2d_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)

* Create grid variables 
      
      var_sublon_id = ncvdef(ncid,'subgrid_longitude',NCPREC,1,
     $     dim_kpt_id,rcode)
      
      var_sublat_id = ncvdef(ncid,'subgrid_latitude' ,NCPREC,1,
     $     dim_kpt_id,rcode)

      var_ixy_id    = ncvdef(ncid,'ixy'              ,NCLONG,1,
     $     dim_lpt_id,rcode)

      var_jxy_id    = ncvdef(ncid,'jxy'              ,NCLONG,1,
     $     dim_lpt_id,rcode)

      dim2_id(1) = dim_lon_id
      dim2_id(2) = dim_lat_id
      var_grdsrf_id = ncvdef(ncid,'surface_type'     ,NCLONG,2,
     $     dim2_id,rcode)

      dim2_id(1) = dim_lpt_id
      dim2_id(2) = dim_sub_id
      var_subkvc_id = ncvdef(ncid,'subgrid_vector'   ,NCLONG,2,
     $     dim2_id,rcode)
      var_subwgt_id = ncvdef(ncid,'subgrid_weight'   ,NCPREC,2,
     $     dim2_id,rcode)

* Create attributes for grid variables

      name = 'longitude for subgrid vector'
      unit = 'degrees_east'
      call ncaptc(ncid,var_sublon_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_sublon_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'latitude for subgrid vector'
      unit = 'degrees_north'
      call ncaptc(ncid,var_sublat_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_sublat_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'land longitude index'
      unit = 'unitless'
      call ncaptc(ncid,var_ixy_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_ixy_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'land latitude index'
      unit = 'unitless'
      call ncaptc(ncid,var_jxy_id,'long_name',NCCHAR,lencs(name),
     $     name ,rcode)
      call ncaptc(ncid,var_jxy_id,'units'   ,NCCHAR,lencs(unit),
     $     unit,rcode)

      name = 'surface type'
      unit = 'unitless'
      call ncaptc(ncid,var_grdsrf_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_grdsrf_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)

      name = 'subgrid vector indices for land points'
      unit = 'unitless'
      call ncaptc(ncid,var_subkvc_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_subkvc_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'subgrid weights for land points'
      unit = 'unitless'
      call ncaptc(ncid,var_subwgt_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_subwgt_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)

* Create time-invariant model timing variables
      
      mdbase_id = ncvdef(ncid,'mdbase',NCLONG,0,0,rcode)
      msbase_id = ncvdef(ncid,'msbase',NCLONG,0,0,rcode)
      mbdate_id = ncvdef(ncid,'mbdate',NCLONG,0,0,rcode)
      mbsec_id  = ncvdef(ncid,'mbsec' ,NCLONG,0,0,rcode)
      dtlsm_id  = ncvdef(ncid,'dtlsm' ,NCFLOAT,0,0,rcode)
      nhtfrq_id = ncvdef(ncid,'nhtfrq',NCLONG,0,0,rcode)

      name = 'base day for this case'
      unit = 'unitless'
      call ncaptc(ncid,mdbase_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,mdbase_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'seconds to complete base day'
      unit = 's'
      call ncaptc(ncid,msbase_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,msbase_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'base date as 6 digit integer (YYMMDD)'
      unit = 'unitless'
      call ncaptc(ncid,mbdate_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,mbdate_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'seconds to complete base date'
      unit = 's'
      call ncaptc(ncid,mbsec_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,mbsec_id,'units '   ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'model timestep'
      unit = 's'
      call ncaptc(ncid,dtlsm_id,'long_name',NCCHAR,lencs(name),
     $     name ,rcode)
      call ncaptc(ncid,dtlsm_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'frequency of model writes (timesteps)'
      unit = 'unitless'
      call ncaptc(ncid,nhtfrq_id,'long_name',NCCHAR,lencs(name),
     $     name ,rcode)
      call ncaptc(ncid,nhtfrq_id,'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
* Create time-invariant input file variables

      inifil_id = ncvdef(ncid,'inifil',NCCHAR,1,dim_str_id,rcode)
      name = 'input initial data'
      call ncaptc(ncid,inifil_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)

      flnfil_id = ncvdef(ncid,'flnfil',NCCHAR,1,dim_str_id,rcode)
      name = 'fractional land data'
      call ncaptc(ncid,flnfil_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)

      lonfil_id = ncvdef(ncid,'lonfil',NCCHAR,1,dim_str_id,rcode)
      name = 'number of longitudes per latitude band data'
      call ncaptc(ncid,lonfil_id,'long_name',NCCHAR,lencs(name),
     $     name,rcode)

* Create time-dependent model timing variables
      
      dim2_id(1) = dim_str_id
      dim2_id(2) = dim_tim_id
      var_tim_id(1) =ncvdef(ncid,'date'    ,NCLONG,1,dim_tim_id,rcode)
      var_tim_id(2) =ncvdef(ncid,'datesec' ,NCLONG,1,dim_tim_id,rcode)
      var_tim_id(3) =ncvdef(ncid,'time'    ,NCPREC,1,dim_tim_id,rcode)
      var_tim_id(4) =ncvdef(ncid,'timestep',NCLONG,1,dim_tim_id,rcode)
      var_tim_id(5) =ncvdef(ncid,'timecom' ,NCCHAR,2,dim2_id   ,rcode)

      name = 'current date as 6 digit integer (YYMMDD)'
      unit = 'unitless'
      call ncaptc(ncid,var_tim_id(1),'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_tim_id(1),'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'seconds into current date'
      unit = 's'
      call ncaptc(ncid,var_tim_id(2),'long_name',NCCHAR,lencs(name),
     $     name ,rcode)
      call ncaptc(ncid,var_tim_id(2),'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'current day'
      unit = 'unitless'
      call ncaptc(ncid,var_tim_id(3),'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_tim_id(3),'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'simulation time step'
      unit = 'unitless'
      call ncaptc(ncid,var_tim_id(4),'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      call ncaptc(ncid,var_tim_id(4),'units'    ,NCCHAR,lencs(unit),
     $     unit,rcode)
      
      name = 'time comment'
      call ncaptc(ncid,var_tim_id(5),'long_name',NCCHAR,lencs(name),
     $     name,rcode)
      
* Create time-dependent variables for single level fields written as xy vectors
      
      dim3_id(1) = dim_lon_id
      dim3_id(2) = dim_lat_id
      dim3_id(3) = dim_tim_id 
      do i=1,nslflds
         slfld_id(i) = ncvdef(ncid,slfnam(i),NCPREC,3,dim3_id,
     $        rcode)
         call ncaptc (ncid,slfld_id(i),'Title'      ,NCCHAR,40,
     $        slfdes(i),rcode)
         call ncaptc (ncid,slfld_id(i),'long_name'  ,NCCHAR,40,
     $        slfdes(i),rcode)
         call ncaptc (ncid,slfld_id(i),'units'      ,NCCHAR, 8,
     $        slfuni(i),rcode)
         call ncaptc (ncid,slfld_id(i),'field_type' ,NCCHAR, 8,
     $        slftyp(i),rcode)
         call ncaptc (ncid,slfld_id(i),'level_flag' ,NCCHAR,12,
     $        'single_level',rcode)
         call ncapt (ncid,slfld_id(i),'missing_value',NCPREC,1,
     $        spval,rcode)
      end do
      
* Create time-dependent variables for multi level fields written as xy vectors
      
      dim4_id(1) = dim_lon_id
      dim4_id(2) = dim_lat_id
      dim4_id(3) = dim_lev_id
      dim4_id(4) = dim_tim_id 
      do i=1,nmlflds
         mlfld_id(i) = ncvdef(ncid,mlfnam(i),NCPREC,4,dim4_id,
     $        rcode)
         call ncaptc (ncid,mlfld_id(i),'Title'       ,NCCHAR,40,
     $        mlfdes(i),rcode)
         call ncaptc (ncid,slfld_id(i),'long_name'   ,NCCHAR,40,
     $        slfdes(i),rcode)
         call ncaptc (ncid,mlfld_id(i),'units'       ,NCCHAR, 8,
     $        mlfuni(i),rcode)
         call ncaptc (ncid,mlfld_id(i),'field_type'  ,NCCHAR, 8,
     $        mlftyp(i),rcode)
         call ncaptc (ncid,slfld_id(i),'level_flag'  ,NCCHAR,11,
     $        'multi_level',rcode)
         call ncapt (ncid,slfld_id(i),'missing_value',NCPREC, 1,
     $        spval,rcode)
      end do
      
* Create time-dependent variables for single-level 1-d fields
      
      dim2_id(1) = dim_kpt_id
      dim2_id(2) = dim_tim_id
      do i=1,nsl1dflds
         sl1dfld_id(i) = ncvdef(ncid,sl1dfnam(i),NCPREC,2,dim2_id,
     $        rcode)
         call ncaptc (ncid,sl1dfld_id(i),'Title'     ,NCCHAR,40,
     $        sl1dfdes(i),rcode)
         call ncaptc (ncid,slfld_id(i),'long_name'   ,NCCHAR,40,
     $        slfdes(i),rcode)
         call ncaptc (ncid,sl1dfld_id(i),'units'     ,NCCHAR, 8,
     $        sl1dfuni(i),rcode)
         call ncaptc (ncid,sl1dfld_id(i),'field_type',NCCHAR, 8,
     $        sl1dftyp(i),rcode)
         call ncaptc (ncid,slfld_id(i),'level_flag'  ,NCCHAR,12,
     $        'single_level',rcode)
         call ncapt (ncid,slfld_id(i),'missing_value',NCPREC, 1,
     $        spval,rcode)
      end do
      
* Create time-dependent variables for multi-level 1-d fields
      
      dim3_id(1) = dim_kpt_id
      dim3_id(2) = dim_lev_id
      dim3_id(3) = dim_tim_id 
      do i=1,nml1dflds
         ml1dfld_id(i) = ncvdef(ncid,ml1dfnam(i),NCPREC,3,dim3_id,
     $        rcode)
         call ncaptc (ncid,ml1dfld_id(i),'Title'     ,NCCHAR,40,
     $        ml1dfdes(i),rcode)
         call ncaptc (ncid,slfld_id(i),'long_name'   ,NCCHAR,40,
     $        slfdes(i),rcode)
         call ncaptc (ncid,ml1dfld_id(i),'units'     ,NCCHAR, 8,
     $        ml1dfuni(i),rcode)
         call ncaptc (ncid,ml1dfld_id(i),'field_type',NCCHAR, 8,
     $        ml1dftyp(i),rcode)
         call ncaptc (ncid,slfld_id(i),'level_flag'  ,NCCHAR,11,
     $        'multi_level',rcode)
         call ncapt (ncid,slfld_id(i),'missing_value',NCPREC, 1,
     $        spval,rcode)
      end do
      
* create basin drainage history field

      if (ncbasin) then
         dim_bas_id = ncddef(ncid,'basins',ndrn,rcode)
         dim2_id(1) = dim_bas_id
         dim2_id(2) = dim_tim_id 
         var_bas_id = ncvdef(ncid,'BASIN_RUNOFF',NCPREC,1,
     $        dim2_id,rcode)
         name = 'total runoff to various basins'
         call ncaptc(ncid,var_bas_id,'long_name',NCCHAR,lencs(name),
     $        name,rcode)
         unit = 'kg/sec'
         call ncaptc(ncid,var_bas_id,'units'    ,NCCHAR,lencs(unit),
     $        unit,rcode)
         name = 'basin 01  Eurasian Arctic'
         call ncaptc(ncid,var_bas_id,'basin_01',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 02  N. American Arctic'
         call ncaptc(ncid,var_bas_id,'basin_02',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 03  N. American Atlantic'
         call ncaptc(ncid,var_bas_id,'basin_03',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 04  European Atlantic'
         call ncaptc(ncid,var_bas_id,'basin_04',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 05  S. American Atlantic'
         call ncaptc(ncid,var_bas_id,'basin_05',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 06  African Atlantic'
         call ncaptc(ncid,var_bas_id,'basin_06',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 07  African Indian'
         call ncaptc(ncid,var_bas_id,'basin_07',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 08  Asian Indian'
         call ncaptc(ncid,var_bas_id,'basin_08',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 09  Asian Pacific'
         call ncaptc(ncid,var_bas_id,'basin_09',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 10  N. American Pacific'
         call ncaptc(ncid,var_bas_id,'basin_10',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 11  S. American Pacific'
         call ncaptc(ncid,var_bas_id,'basin_11',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 12  Australian Southern'
         call ncaptc(ncid,var_bas_id,'basin_12',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 13  Antarctic Southern'
         call ncaptc(ncid,var_bas_id,'basin_13',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 14  Mediterranean'
         call ncaptc(ncid,var_bas_id,'basin_14',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 15  Caspian Sea'
         call ncaptc(ncid,var_bas_id,'basin_15',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 16  Black Sea'
         call ncaptc(ncid,var_bas_id,'basin_16',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 17  Red Sea / Persian Gulf'
         call ncaptc(ncid,var_bas_id,'basin_17',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 18  Baltic Sea'
         call ncaptc(ncid,var_bas_id,'basin_18',NCCHAR,lencs(name),
     $        name,rcode)
         name = 'basin 19  Hudson Bay'
         call ncaptc(ncid,var_bas_id,'basin_19',NCCHAR,lencs(name),
     $        name,rcode)
      endif


* finish creating netcdf file
      
      call ncendf(ncid,rcode)

 1000 format('Data created ',a,1x,a,' $')

      return
      end

 
