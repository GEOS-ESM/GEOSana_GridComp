#!/usr/local/other/python/GEOSpyD/2019.03_py3.7/2019-04-22/bin/python
'''
Python version of ana2inc.F90. Unlike the F90 version, this code does 
not require a pre-copy of ana.eta -> inc.eta, but rather makes a copy
of ana.eta.nc4 itself. Increments are computed for all values of the 
anafile except for the ones listed in SKIPVARS below.

EXAMPLES:
python ana2inc.py -a ana.eta.nc4 -b bkg.eta.nc4 -c cbkg.eta.nc4 -o inc.eta.nc4

HISTORY: 
20240201 - christoph.a.keller@nasa.gov - initial version
'''
import sys
import argparse
import logging
import datetime as dt
import time
import numpy as np
import xarray as xr
import pandas as pd

SKIPVARS = ["phis", "ps", "ts"]

def main(args):
    '''
    Calculate increments from analysis and background fields
    '''
    log = logging.getLogger(__name__)
    # read ana, bkg, and cbkg file
    # do not decode time stamps as this screws up the units 
    # attribute.
    ana = xr.open_dataset(args.anafile,decode_times=False)
    bkg = xr.open_dataset(args.bkgfile,decode_times=False)
    cbk = xr.open_dataset(args.cbkgfile,decode_times=False)
    # make copy of anafile for the increment file
    inc = ana.copy()
    # loop over all variables in inc and calculate increments
    for v in inc:
        if v in SKIPVARS:
            continue
        if v in bkg:
            inc[v].values = inc[v].values - bkg[v].values
        else:
            assert v in cbk, "variable not found: {}".format(v)
            inc[v].values = inc[v].values - cbk[v].values

    # write out
    inc.to_netcdf(args.ofile)
    log.info('file written to {}'.format(args.ofile))
    inc.close()
    ana.close()
    bkg.close()
    cbk.close()
    return


def parse_args():
    p = argparse.ArgumentParser(description='Undef certain variables')
    p.add_argument('-a', '--anafile',type=str,help='analysis file',default="ana.eta.nc4")
    p.add_argument('-b', '--bkgfile',type=str,help='background file',default="bkg.eta.nc4")
    p.add_argument('-c', '--cbkgfile',type=str,help='chemistry background file',default="cbkg.eta.nc4")
    p.add_argument('-o', '--ofile',type=str,help='output file', default="inc.eta.nc4")
    return p.parse_args()


if __name__ == '__main__':
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    log.addHandler(handler)
    main(parse_args())

