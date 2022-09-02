#!/usr/bin/env python3

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

        if not self.debug:
            os.remove(outfil)

#.......................................................................
if __name__ == "__main__":

    # -db flag will keep outfil from being deleted
    #---------------------------------------------
    TlapmeanTest.debug = False
    for arg in sys.argv[1:len(sys.argv)]:
        if arg == "-v": continue
        if arg.lower() in ["-db", "-debug"]:
            TlapmeanTest.debug = True
        sys.argv.remove(arg)

    # check for existence of outdir
    #------------------------------
    if not os.path.isdir("outdir"):
        print("> making outdir directory")
        os.mkdir("outdir")

    # run test
    #---------
    unittest.main()
