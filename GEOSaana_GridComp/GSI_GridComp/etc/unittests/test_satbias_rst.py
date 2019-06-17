#!/usr/bin/env python

import os, sys
sys.path.append(os.pardir)

import unittest
import filecmp
from satbias_rst import SatBiasRst

#.......................................................................
class SatbiasTest(unittest.TestCase):

    def test_append_duplicate(self):
        "Append_zero_record() should fail if key is already in data"
        infile = os.path.join("input", "satbias_full_rst")
        sb = SatBiasRst(infile)

        # create record with duplicate key
        (key, values) = sb.key_and_values(10)

        # attempt to append record with duplicate key
        self.assertRaises(Exception, sb.append_zero_record, key)

    def test_append_zero_record(self):
        "Test the append_zero_record() function"
        infile = os.path.join("input", "satbias_full_rst")
        outfil = os.path.join("outdir", "satbias_rst.append_zero_record")
        outexp = os.path.join("outexp", "satbias_rst.append_zero_record")

        sb = SatBiasRst(infile)

        # create new values to load from last record
        (key, values) = sb.key_and_values(sb.len()-1)
        (instrument, channel) = key

        channel = "%s" % (int(channel) + 1)
        key = (instrument, channel)

        # append record to data and write output
        sb.append_zero_record(key)

        self.assertEqual(sb.filename(), infile)
        sb.write(outfil)
        self.assertEqual(sb.filename(), outfil)

        # compare to expected output
        self.assertTrue(filecmp.cmp(outfil, outexp))

        os.remove(outfil)

    def test_change(self):
        "Test the change() function for various variables"
        infile = os.path.join("input", "satbias_full_rst")

        infile1 = os.path.join("outdir", "satbias_full_rst")
        outexp1 = os.path.join("outexp", "satbias_full_rst")

        outfil = os.path.join("outdir", "satbias_test.new")
        outexp = os.path.join("outexp", "satbias_test.new")

        # rewrite input to remove comments and extraneous blanks
        # (for diff'ing purposes while verifying test)
        sb1 = SatBiasRst(infile)
        sb1.write(infile1)
        self.assertTrue(filecmp.cmp(infile1, outexp1))

        # change instrument in record #1
        sb2 = SatBiasRst(infile1)
        (key, values) = sb2.key_and_values(0)
        (instrument, channel) = key
        sb2.change(key, "instrument", "istoo_n09")
        
        # change channel in record #2
        (key, values) = sb2.key_and_values(1)
        (instrument, channel) = key
        sb2.change(key, "channel", "72")
        
        # change coeffs in record #3
        (key, values) = sb2.key_and_values(2)
        coeffs = ["%7.3f"%(float(n)/2) for n in range(7)]
        sb2.change(key, "coeffs", coeffs)

        # write output
        sb2.write(outfil)

        # compare to expected output
        self.assertTrue(filecmp.cmp(outfil, outexp))

        os.remove(infile1)
        os.remove(outfil)

    def test_data_and_records(self):
        "Test that records() agrees with self._data[] and key_and_values()"
        infile = os.path.join("input", "satbias_full_rst")
        sb = SatBiasRst(infile)

        index = 0
        for (key, values) in sb.records():
            self.assertEqual((key, values), sb._data[index])
            self.assertEqual((key, values), sb.key_and_values(index))
            index += 1

    def test_duplicate_record(self):
        "Load() should fail if satbias has duplicate records"
        infile = os.path.join("input", "satbias_full_rst")
        infilX = os.path.join("input", "satbias_dupkey_rst")

        # attempt to add duplicate key
        sb = SatBiasRst(infile)
        (key1, values) = sb.key_and_values(10)
        (key2, values) = sb.key_and_values(11)
        self.assertRaises(Exception, sb.change, key2, "key", key1)

        # attempt to load file with duplicate key
        self.assertRaises(Exception, SatBiasRst, infilX, 1)

    def test_ioerror_coeffs(self):
        "Load() should fail if coeffs does not have 7 elements"
        infile = os.path.join("input", "satbias_full_rst")
        outfil = os.path.join("outdir", "satbias_rst.ioerror_coeffs")
        outexp = os.path.join("outexp", "satbias_rst.ioerror_coeffs")

        # blank last value of coeff
        sb = SatBiasRst(infile)
        (key, values) = sb.key_and_values(10)
        values["coeffs"][6] = ""
        sb.change(key, "coeffs", values["coeffs"])

        # this write works because coeffs is full length but with a blank
        sb.write(outfil)
        self.assertTrue(filecmp.cmp(outfil, outexp))

        # attempt to load faulty output file
        self.assertRaises(Exception, SatBiasRst, outfil)

        os.remove(outfil)

    def test_nonexist(self):
        "Load() should fail if file does not exist"
        infile = os.path.join("input", "satbias_nonexist_rst")
        self.assertRaises(Exception, SatBiasRst, infile)

    def test_read_inconsistent_recnum(self):
        "Load() should fail if input has inconsistent recnum"
        infilX = os.path.join("input", "satbias_badrecnum_rst")
        self.assertRaises(Exception, SatBiasRst, infilX)

    def test_read_malformatted_rec(self):
        "Load() should fail if file has malformatted record"
        infilX = os.path.join("input", "satbias_malrec_rst")
        self.assertRaises(Exception, SatBiasRst, infilX)

    def test_read_write(self):
        "Test that __read() and write() are inverse functions"
        infile = os.path.join("input", "satbias_full_rst")
        outfil = os.path.join("outdir", "satbias_full_rst.copy")

        sb1 = SatBiasRst(infile)
        sb1.write(outfil)

        # load outpout
        sb2 = SatBiasRst(outfil)

        # compare private data from input and output; check for file equivalence
        self.assertEqual(sb1._data, sb2._data)
        self.assertTrue(filecmp.cmp(infile, outfil))

        os.remove(outfil)

    def test_remove_records(self):
        "Test the remove_records() function"
        infile = os.path.join("input", "satbias_full_rst")
        outfil = os.path.join("outdir", "satbias_rst.remove_records")
        outexp = os.path.join("outexp", "satbias_rst.remove_records")

        sb = SatBiasRst(infile)

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

        os.remove(outfil)

#.......................................................................
if __name__ == "__main__":
    unittest.main()
