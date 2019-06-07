DSET   ^zavg.data
TITLE  fvgcm_b55 dycore test
OPTIONS  big_endian sequential
UNDEF  0.100000E+26
XDEF    1  LINEAR      0.000000     2.500000
YDEF   91  LINEAR    -90.000000     2.000000
ZDEF   55  LEVELS
   992.5559  970.5549  929.6490  867.1610  787.7020
  696.7955  600.5240  510.4550  433.8950  368.8180
  313.5010  266.4810  226.5130  192.5395  163.6615
  139.1150  118.2502  100.5145   85.4390   72.5578
   61.4957   52.0159   43.9097   36.9927   31.0889
   26.0491   21.7610   18.1243   15.0502   12.4601
   10.2849    8.4564    6.9183    5.6318    4.5617
    3.6765    2.9483    2.3526    1.8679    1.4756
    1.1600    0.9073    0.7060    0.5463    0.4204
    0.3218    0.2449    0.1854    0.1396    0.1045
    0.0777    0.0568    0.0401    0.0263    0.0150
TDEF    12 LINEAR 00:00Z01SEP0050   1mo
VARS  12
ps     0   0 surface pressure
ux    55   0 zonal mean u-wind (m/s)
vx    55   0 zonal mean v-wind (m/s)
wx    55   0 dpdt (mb/day)
TX    55   0 t (deg K)
qx    55   0 tracer
pvx   55   0 pv
usus  55   0 u variance
usvs  55   0 meridional eddy u-momentum transport
tsts  55   0 t variance
vsts  55   0 meridional eddy heat transport
vspv  55   0 meridional eddy PV transport
ENDVARS
