* 1-d vector length

      integer kpt   ! total number of lsm land points, including subgrid points
      integer lpt   ! number of land points on lsmlon x lsmlat grid
      integer numlv ! # "little" vectors = calls to process all subgrid points

      common /lsmmem/ kpt, lpt, numlv

 
