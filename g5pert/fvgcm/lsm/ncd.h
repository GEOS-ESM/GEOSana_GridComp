* ------------------------ code history ---------------------------
* source file:       ncd.h
* purpose:           common block netcdf history files 
* date last revised: August 1996 (M. Vertenstein)
* author:            M. Vertenstein
* -----------------------------------------------------------------

* netcdf stuff

      integer rcode
      integer ncid

      common /ncdinfo/ ncid 

      real       spval
      parameter (spval=1.e30)

* grid variables

      integer var_lon_id
      integer var_lat_id
      integer var_lev_id
      integer var_ixy_id
      integer var_jxy_id
      integer var_grdsrf_id
      integer var_sublat_id
      integer var_sublon_id
      integer var_subkvc_id
      integer var_subwgt_id
      integer var_lon2d_id
      integer var_lat2d_id
      integer var_mask2d_id

      common /ncdgrid/ var_lon_id   , var_lat_id   , var_lev_id   ,
     $                 var_grdsrf_id, var_ixy_id   , var_jxy_id   ,
     $                 var_sublat_id, var_sublon_id, var_subkvc_id,
     $                 var_subwgt_id, var_lon2d_id , var_lat2d_id ,
     $                 var_mask2d_id

* dataset variables

      integer lonfil_id
      integer flnfil_id   
      integer inifil_id   

      common /ncdfil/ flnfil_id, inifil_id, lonfil_id

* time variables

      integer mdbase_id
      integer msbase_id
      integer mbdate_id
      integer mbsec_id
      integer dtlsm_id
      integer nhtfrq_id
      integer var_tim_id(5)
      
      common /ncdtime/ mdbase_id, msbase_id, mbdate_id, mbsec_id,
     $                 dtlsm_id , nhtfrq_id, var_tim_id

* field variables

      integer slfld_id(mslflds)
      integer mlfld_id(mmlflds)
      integer sl1dfld_id(mslflds)
      integer ml1dfld_id(mmlflds)

      common /ncdflds/ slfld_id, mlfld_id, sl1dfld_id, ml1dfld_id
 
