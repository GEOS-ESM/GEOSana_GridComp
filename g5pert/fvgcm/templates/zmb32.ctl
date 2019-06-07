DSET   ^zavg.data
TITLE  fvgcm_b32 dycore test
OPTIONS  big_endian sequential
UNDEF  0.100000E+26
XDEF   1  LINEAR      0.000000     2.500000
YDEF   91  LINEAR    -90.000000     2.000000
ZDEF   32  LEVELS
   992.5550  970.5525  929.6455  867.1572  787.7000
  696.7930  600.5238  510.4555  433.8927  368.8161
  313.4988  266.4789  226.5135  192.5410  163.6615
  139.1150  118.2502  100.5145   85.4390   72.2925
   60.5613   50.2000   41.0750   32.9000   25.3000
   18.3500   12.5500    8.1500    5.0000    2.8500
    1.5000    0.7000
TDEF    12 LINEAR 00:00Z01SEP0050   1mo
VARS  12
ps     0   0 zonal mean surface pressure (pa)
ux    32   0 zonal mean u-wind (m/s)
vx    32   0 zonal mean v-wind (m/s)
wx    32   0 dpdt (mb/day)
TX    32   0 t (deg K)
qx    32   0 tracer
pvx   32   0 pv
usus  32   0 u variance
usvs  32   0 meridional eddy u-momentumtransport
tsts  32   0 t variance
vsts  32   0 meridional eddy heat transport
vspv  32   0 meridional eddy PV transport
ENDVARS
