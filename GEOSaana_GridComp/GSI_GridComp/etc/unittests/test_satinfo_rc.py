#!/usr/bin/env python3

import os, sys
sys.path.append(os.pardir)

import unittest
import filecmp
from satinfo_rc import SatInfoRc

#.......................................................................
class SatinfoTest(unittest.TestCase):

    def test_read_write(self):
        "Test that __read() and write() are inverse functions"
        infile = "input/gmao_global_satinfo.rc"
        outfil = "outdir/gmao_global_satinfo.rc.copy"

        f1 = SatInfoRc(infile)
        self.assertEqual(infile, f1.filename())

        f1.write(outfil)
        self.assertTrue(filecmp.cmp(infile, outfil))
        self.assertEqual(outfil, f1.filename())

        f2 = SatInfoRc(outfil)
        self.assertEquals(f1._data, f2._data)

        if not self.debug:
            os.remove(outfil)

    def test_change_list(self):
        "Test that self.change_list equals self.change_history"
        infile =  "input/gmao_global_satinfo.rc"
        si = SatInfoRc(infile)
        self.assertEqual(sorted(si.change_list), sorted(si.change_history))

#.......................................................................
if __name__ == "__main__":

    # -db flag will keep outfil from being deleted
    #---------------------------------------------
    SatinfoTest.debug = False
    for arg in sys.argv[1:len(sys.argv)]:
        if arg == "-v": continue
        if arg.lower() in ["-db", "-debug"]:
            SatinfoTest.debug = True
        sys.argv.remove(arg)

    # check for existence of outdir
    #------------------------------
    if not os.path.isdir("outdir"):
        print("> making outdir directory")
        os.mkdir("outdir")

    # run test
    #---------
    unittest.main()
