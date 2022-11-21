#!/usr/bin/env python3

import os, sys
sys.path.append(os.pardir)

import unittest
import filecmp
from satbang_rst import SatBangRst

#.......................................................................
class SatbangTest(unittest.TestCase):

    def test_append_duplicate(self):
        "Append_zero_record() should fail if key is already in data"
        infile = os.path.join("input", "satbang_test_rst")
        sb = SatBangRst(infile)

        # create record with duplicate key
        (key, values) = sb.key_and_values(10)
        tlapmean = values["tlapmean"].replace("0", "1")

        # attempt to append record with duplicate key
        self.assertRaises(Exception, sb.append_zero_record, key, tlapmean)

    def test_append_zero_record(self):
        "Test the append_zero_record() function"
        infile = os.path.join("input", "satbang_test_rst")
        outfil = os.path.join("outdir", "satbang_rst.append_zero_record")
        outexp = os.path.join("outexp", "satbang_rst.append_zero_record")

        sb = SatBangRst(infile)

        # create new values to load from last record
        (key, values) = sb.key_and_values(sb.len()-1)
        (instrument, channel) = key

        channel = "%s" % (int(channel) + 1)
        key = (instrument, channel)
        tlapmean = "%12.6E" % (float(values["tlapmean"]) - 0.05)

        # append record to data and write output
        sb.append_zero_record(key, tlapmean)

        self.assertEqual(sb.filename(), infile)
        sb.write(outfil)
        self.assertEqual(sb.filename(), outfil)

        # compare to expected output
        self.assertTrue(filecmp.cmp(outfil, outexp))

        if not self.debug:
            os.remove(outfil)

    def test_change(self):
        "Test the change() function for various variables"
        infile = os.path.join("input", "satbang_test_rst")

        infile1 = os.path.join("outdir", "satbang_test_rst")
        outexp1 = os.path.join("outexp", "satbang_test_rst")

        outfil = os.path.join("outdir", "satbang_test.new")
        outexp = os.path.join("outexp", "satbang_test.new")

        # rewrite input to remove comments and extraneous blanks
        # (for diff'ing purposes while verifying test)
        sb1 = SatBangRst(infile)
        sb1.write(infile1)
        self.assertTrue(filecmp.cmp(infile1, outexp1))

        # change instrument in record #1
        sb2 = SatBangRst(infile1)
        (key, values) = sb2.key_and_values(0)
        (instrument, channel) = key
        sb2.change(key, "instrument", "istoo_n09")
        
        # change channel in record #2
        (key, values) = sb2.key_and_values(1)
        (instrument, channel) = key
        sb2.change(key, "channel", "12")
        
        # change tlapmean in record #3
        (key, values) = sb2.key_and_values(2)
        sb2.change(key, "tlapmean", "0.455304E-01")
        
        # change coeffs in record #4
        (key, values) = sb2.key_and_values(3)
        coeffs = ["%7.3f"%(((-1)**n)*float(n)/8) for n in range(90)]
        sb2.change(key, "coeffs", coeffs)

        # change key in record #5
        (key, values) = sb2.key_and_values(4)
        sb2.change(key, "key", ("amsub_n61", "44"))
        
        # write output
        sb2.write(outfil)

        # compare to expected output
        self.assertTrue(filecmp.cmp(outfil, outexp))

        if not self.debug:
            os.remove(infile1)
            os.remove(outfil)

    def test_data_and_records(self):
        "Test that records() agrees with self._data[] and key_and_values()"
        infile = os.path.join("input", "satbang_full_rst")
        sb = SatBangRst(infile)

        index = 0
        for (key, values) in sb.records():
            self.assertEqual((key, values), sb._data[index])
            self.assertEqual((key, values), sb.key_and_values(index))
            index += 1

    def test_duplicate_key(self):
        "Load() should fail if satbang has duplicate keys"
        infile = os.path.join("input", "satbang_test_rst")
        infilX = os.path.join("input", "satbang_dupkey_rst")

        # attempt to add duplicate key
        sb = SatBangRst(infile)
        (key1, values) = sb.key_and_values(10)
        (key2, values) = sb.key_and_values(11)
        self.assertRaises(Exception, sb.change, key2, "key", key1)

        # attempt to load file with duplicate key
        self.assertRaises(Exception, SatBangRst, infilX, 1)

    def test_ioerror_coeffs(self):
        "Load() should fail if coeffs does not have 90 elements"
        infile = os.path.join("input", "satbang_test_rst")
        outfil = os.path.join("outdir", "satbang_rst.ioerror_coeffs")
        outexp = os.path.join("outexp", "satbang_rst.ioerror_coeffs")

        # blank last value of coeff
        sb = SatBangRst(infile)
        (key, values) = sb.key_and_values(10)
        values["coeffs"][89] = ""
        sb.change(key, "coeffs", values["coeffs"])

        # this write works because coeffs is full length but with a blank
        sb.write(outfil)
        self.assertTrue(filecmp.cmp(outfil, outexp))

        # attempt to load faulty output file
        self.assertRaises(Exception, SatBangRst, outfil)
        if not self.debug:
            os.remove(outfil)

        # attempt to replace coeffs with write file with faulty record
        values["coeffs"].pop();
        
        self.assertRaises(Exception, sb.change, key, "coeffs", values["coeffs"])

    def test_nonexist(self):
        "Load() should fail if file does not exist"        
        infile = os.path.join("input", "satbang_nonexist_rst")
        self.assertRaises(Exception, SatBangRst, infile)

    def test_read_inconsistent_recnum(self):
        "Load() should fail if input has inconsistent recnum"
        infilX = os.path.join("input", "satbang_badrecnum_rst")
        self.assertRaises(Exception, SatBangRst, infilX)

    def test_read_malformatted_rec(self):
        "Load() should fail if file has malformatted record"
        infilX = os.path.join("input", "satbang_malrec_rst")
        self.assertRaises(Exception, SatBangRst, infilX)
        
    def test_read_write(self):
        "Test that __read() and write() are inverse functions"
        infile = os.path.join("input", "satbang_full_rst")
        outfil = os.path.join("outdir", "satbang_full.copy")

        sb1 = SatBangRst(infile)
        sb1.write(outfil)

        # load output
        sb2 = SatBangRst(outfil)

        # compare private data from input and output; check for file equivalence
        self.assertEqual(sb1._data, sb2._data)
        self.assertTrue(filecmp.cmp(infile, outfil))

        if not self.debug:
            os.remove(outfil)

    def test_remove_records(self):
        "Test the remove_records() function"
        infile = os.path.join("input", "satbang_test_rst")
        outfil = os.path.join("outdir", "satbang_rst.remove_records")
        outexp = os.path.join("outexp", "satbang_rst.remove_records")

        sb = SatBangRst(infile)

        # create list of keys to remove
        keylist = []
        (key, values) = sb.key_and_values(2)
        keylist.append(key)

        (key, values) = sb.key_and_values(5)
        keylist.append(key)

        # remove records and write output
        sb.remove_records(keylist)
        sb.write(outfil)

        # compare to expected output
        self.assertTrue(filecmp.cmp(outfil, outexp))

        if not self.debug:
            os.remove(outfil)

#.......................................................................
if __name__ == "__main__":

    # -db flag will keep outfils from being deleted
    #----------------------------------------------
    SatbangTest.debug = False
    for arg in sys.argv[1:len(sys.argv)]:
        if arg == "-v": continue
        if arg.lower() in ["-db", "-debug"]:
            SatbangTest.debug = True
        sys.argv.remove(arg)

    # check for existence of outdir
    #------------------------------
    if not os.path.isdir("outdir"):
        print("> making outdir directory")
        os.mkdir("outdir")

    # run test
    #---------
    unittest.main()

## found below in documentation but unable to get it to work
#suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
#unittest.TextTestRunner(verbosity=2).run(suite)
