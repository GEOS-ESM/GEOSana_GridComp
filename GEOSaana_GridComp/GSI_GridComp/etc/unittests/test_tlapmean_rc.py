#!/usr/bin/env python

import os, sys
sys.path.append(os.pardir)

import unittest
import filecmp
from tlapmean_rc import TLapMeanRc

#.......................................................................
class TlapmeanTest(unittest.TestCase):

    def test_read_write(self):
        "Test that __read() and write() are inverse functions"
        infile = "input/gmao_global_tlapmean.rc"
        outfil = "outdir/gmao_global_tlapmean.rc.copy"

        f1 = TLapMeanRc(infile)
        f1.write(outfil)
        self.assertTrue(filecmp.cmp(infile, outfil))

        f2 = TLapMeanRc(outfil)
        self.assertEquals(f1._data, f2._data)

        os.remove(outfil)

#.......................................................................
if __name__ == "__main__":
    unittest.main()
