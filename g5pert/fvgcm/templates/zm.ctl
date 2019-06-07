DSET   ^zavg.data
TITLE  fvgcm_b32 dycore test
OPTIONS  big_endian sequential
UNDEF  0.100000E+26
XDEF   1  LINEAR      0.000000     500000
YDEF   91  LINEAR    -90.000000    2.000000
ZDEF   26  LEVELS
  992.5561  970.5548  929.6488  867.1608  787.7020
  696.7963  600.5242  510.4552  433.8952  368.8180
  313.5013  266.4811  226.5133  192.5399  163.6621
  139.1154  118.2503  100.5147   85.4391   70.0592
   53.1146   37.2303   23.9446   13.9672    7.3888
    3.5446
TDEF  2222 LINEAR 00:00Z01SEP0050   1mo
VARS  12
psx    0   0 zonal mean surface pressure
ux    26   0 zonal mean u-wind (m/s)
vx    26   0 zonal mean v-wind (m/s)
wx    26   0 dpdt (mb/day)
TX    26   0 t (deg K)
qx    26   0 tracer
pvx   26   0 pv
usus  26   0 u variance
usvs  26   0 meridional eddy u-momentumtransport
tsts  26   0 t variance
vsts  26   0 meridional eddy heat transport
vspv  26   0 meridional eddy PV transport
ENDVARS
