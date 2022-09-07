#!/usr/bin/env python3

import os, sys
sys.path.append(os.pardir)

import unittest
import filecmp, shutil

from check_satbang import check_files

#.......................................................................
class CheckSatbangTest(unittest.TestCase):

    def test_check_satbang(self):
        """
        Run check_satbang for input satbang and satbias restarts
        and compare results to expected output.
        """
        satbang = os.path.join("input", "satbang_full_rst")
        satbias = os.path.join("input", "satbias_full_rst")
        satinfo = os.path.join("input", "gmao_global_satinfo.rc")
        tlapmean = os.path.join("input", "gmao_global_tlapmean.rc")

        for sortflag in ("default", "alpha", "satinfo", "satbias"):

            fname = "satbang_{}_rst".format(sortflag)
            satbangIN = os.path.join("outdir", fname)
            satbangOUT = os.path.join("outdir", fname+".new")
            sumbangOUT = os.path.join("outdir", "check_summary."+fname)
            satbangEXP = os.path.join("outexp", fname+".new")
            sumbangEXP = os.path.join("outexp", "check_summary."+fname)

            fname = "satbias_{}_rst".format(sortflag)
            satbiasIN = os.path.join("outdir", fname)
            satbiasOUT = os.path.join("outdir", fname+".new")
            sumbiasOUT = os.path.join("outdir", "check_summary."+fname)
            satbiasEXP = os.path.join("outexp", fname+".new")
            sumbiasEXP = os.path.join("outexp", "check_summary."+fname)

            sumbangOUT_ = sumbangOUT+"_"
            sumbiasOUT_ = sumbiasOUT+"_"
            sumbangEXP_ = sumbangEXP+"_"
            sumbiasEXP_ = sumbiasEXP+"_"

            # remove leftover output and temporary input
            for file in (satbangIN, satbangOUT, sumbangOUT,
                         satbiasIN, satbiasOUT, sumbiasOUT):
                if os.path.isfile(file):
                    os.remove(file)

            # copy satbang and satbias to outdir
            # because outputs will go to same directory
            shutil.copyfile(satbang, satbangIN)
            shutil.copyfile(satbias, satbiasIN)

            # run test
            check_files(satbangIN, satbiasIN, satinfo, tlapmean, sortflag)

            # remove directory-dependent line from summary files
            for sumfil in [sumbangOUT, sumbiasOUT]:
                sumfix = sumfil+"_"
                with open(sumfil, mode="r") as input:
                    with open (sumfix, mode="w") as output:
                        for line in input:
                            if line.find("current dir") == -1:
                                output.write(line)

            # compare output to expected output
            self.assertTrue(filecmp.cmp(satbangOUT, satbangEXP))
            self.assertTrue(filecmp.cmp(satbiasOUT, satbiasEXP))

            self.assertTrue(filecmp.cmp(sumbangOUT_, sumbangEXP_))
            self.assertTrue(filecmp.cmp(sumbiasOUT_, sumbiasEXP_))

            # remove output and temporary input
            if not self.debug:
                os.remove(satbangIN)
                os.remove(satbiasIN)

                os.remove(satbangOUT)
                os.remove(satbiasOUT)

                os.remove(sumbangOUT)
                os.remove(sumbiasOUT)

                os.remove(sumbangOUT_)
                os.remove(sumbiasOUT_)

    def test_check_satbang_revert(self):
        """
        Run testcheck_satbang with an old version of satinfo.rc, where
        new instrument names and channels numbers need to be reverted to
        old names.
        """
        satbang = os.path.join("outexp", "satbang_default_rst.new")
        satbias = os.path.join("outexp", "satbias_default_rst.new")
        satinfo = os.path.join("input", "gmao_revert_satinfo.rc")
        tlapmean = os.path.join("input", "gmao_global_tlapmean.rc")

        fname = "satbang_revert_rst"
        satbangIN = os.path.join("outdir", fname)
        satbangOUT = os.path.join("outdir", fname+".new")
        sumbangOUT = os.path.join("outdir", "check_summary."+fname)
        satbangEXP = os.path.join("outexp", fname+".new")
        sumbangEXP = os.path.join("outexp", "check_summary."+fname)

        fname = "satbias_revert_rst"
        satbiasIN = os.path.join("outdir", fname)
        satbiasOUT = os.path.join("outdir", fname+".new")
        sumbiasOUT = os.path.join("outdir", "check_summary."+fname)
        satbiasEXP = os.path.join("outexp", fname+".new")
        sumbiasEXP = os.path.join("outexp", "check_summary."+fname)

        sumbangOUT_ = sumbangOUT+"_"
        sumbiasOUT_ = sumbiasOUT+"_"
        sumbangEXP_ = sumbangEXP+"_"
        sumbiasEXP_ = sumbiasEXP+"_"

        # remove leftover output and temporary input
        for file in (satbangIN, satbangOUT, sumbangOUT,
                     satbiasIN, satbiasOUT, sumbiasOUT):
            if os.path.isfile(file):
                os.remove(file)

        # copy satbang and satbias to outdir
        # because outputs will go to same directory
        shutil.copyfile(satbang, satbangIN)
        shutil.copyfile(satbias, satbiasIN)

        # run test
        check_files(satbangIN, satbiasIN, satinfo, tlapmean)

        # remove directory-dependent line from summary files
        for sumfil in [sumbangOUT, sumbiasOUT]:
            sumfix = sumfil+"_"
            with open(sumfil, mode="r") as input:
                with open (sumfix, mode="w") as output:
                    for line in input:
                        if line.find("current dir") == -1:
                            output.write(line)

        # compare output to expected output
        self.assertTrue(filecmp.cmp(satbangOUT, satbangEXP))
        self.assertTrue(filecmp.cmp(satbiasOUT, satbiasEXP))

        self.assertTrue(filecmp.cmp(sumbangOUT_, sumbangEXP_))
        self.assertTrue(filecmp.cmp(sumbiasOUT_, sumbiasEXP_))

        # remove output and temporary input
        if not self.debug:
            os.remove(satbangIN)
            os.remove(satbiasIN)

            os.remove(satbangOUT)
            os.remove(satbiasOUT)

            os.remove(sumbangOUT)
            os.remove(sumbiasOUT)

            os.remove(sumbangOUT_)
            os.remove(sumbiasOUT_)

#.......................................................................
if __name__ == "__main__":

    # -db flag will keep outfils from being deleted
    #----------------------------------------------
    CheckSatbangTest.debug = False
    for arg in sys.argv[1:len(sys.argv)]:
        if arg == "-v": continue
        if arg.lower() in ["-db", "-debug"]:
            CheckSatbangTest.debug = True
        sys.argv.remove(arg)

    # check for existence of outdir
    #------------------------------
    if not os.path.isdir("outdir"):
        print("> making outdir directory")
        os.mkdir("outdir")

    # run test
    #---------
    unittest.main()
