#!/usr/bin/env python

import glob, os, re, sys

from satbang_rst import SatBangRst, satbang_dir_and_name
from satbias_rst import SatBiasRst, satbias_dir_and_name
from satinfo_rc  import SatInfoRc
from tlapmean_rc import TLapMeanRc

def check_files(satbang_path=os.curdir,
                satbias_path=os.curdir,
                satinfo_path="default",
                tlapmean_path="default",
                sortflag="default"):
    """
    Check for and fix irregularities in the satbang and satbias restart files
    when compared to information in the satinfo.rc and tlapmean.rc files.

    satbang and satbias restarts
    ----------------------------
    * Remove all records where the instrument channel is not listed in the
      satinfo.rc file.

    * Add a new record with coeffs equal zero for each instrument channel
      found in the satinfo.rc file but not in the restart file.

    * Sort records according to sortflag, if not equal "default"

    satbang only
    ------------
    * If the tlapmean value is zero, then replace it with the value from
      tlapmean.rc, if one is found there.

    * Make a note of the instrument channels with tlapmean equal zero, and
      there is no value in tlapmean.rc to replace it.

    * Make a note of the instrument channels with a non-zero tlapmean that
      disagrees with the value in tlapmean.rc (but do not replace).

    Input parameters:
    ----------------
    => satbang_path: name of input satbang restart file or name of directory 
                     containing satbang file named "*ana_satbang_rst*.txt".

    => satbias_path: name of input satbias restart file or name of directory 
                     containing satbias file named "*ana_satbias_rst*.txt".

    => satinfo_path: name of input satinfo.rc file; if "default" ,then use
                     "../run/gmao_global_satinfo.rc" relative to the satbang
                     restart directory.

    => tlapmean_path: name of input tlapmean.rc file; if "default", then use
                      "../run/gmao_global_tlapmean.rc" relative to the satbang
                      restart directory.

    => sortflag: flag indicating what order to write output records
                 "default" = keep output records in the same order as the input
                 "satbang" = write output records in same order as satbang_rst
                 "satbias" = write output records in same order as satbias_rst
                 "satinfo" = write output records in same order as satinfo.rc
                 "alpha"   = write output records in alphabetic order
    """

    # get satbang and satbias directory and filename
    #-----------------------------------------------
    (satbang_dir, satbang_name) = satbang_dir_and_name(satbang_path)
    (satbias_dir, satbias_name) = satbias_dir_and_name(satbias_path)

    satbang_path = os.path.join(satbang_dir, satbang_name)
    satbias_path = os.path.join(satbias_dir, satbias_name)

    # get default satinfo_path and tlapmean_path, if necessary
    #---------------------------------------------------------
    if satbang_dir.split(os.sep)[-1] == "recycle":
        run_dir = re.sub(r"recycle$", "run", satbang_dir) 
    else:
        run_dir = os.path.relpath(os.path.join(satbang_dir, "..", "run"))
        if run_dir.count("..") > 1:
            run_dir = os.path.abspath(os.path.join(satbang_dir, "..", "run"))

    if satinfo_path == "default":
        satinfo_path = os.path.join(run_dir, "gmao_global_satinfo.rc")

    if tlapmean_path == "default":
        tlapmean_path = os.path.join(run_dir, "gmao_global_tlapmean.rc")

    # check for existence of input files
    #-----------------------------------
    if not os.path.isfile(satinfo_path):
        raise Exception("{0} does not exist".format(satinfo_path))

    if not os.path.isfile(tlapmean_path):
        raise Exception("{0} does not exist".format(tlapmean_path))

    # create file instances
    #----------------------
    satbang = SatBangRst(satbang_path)
    satbias = SatBiasRst(satbias_path)
    satinfo = SatInfoRc(satinfo_path)
    tlapmean = TLapMeanRc(tlapmean_path)

    # initialize data collectors and more
    #------------------------------------
    dcTL = {}  # data collector for tlapmean info (in satbang)
    dcBA = {}  # data collector for satbang records
    dcBI = {}  # data collector for satbias records

    dcTL["fill"] = []  # zero-value replaced by value from tlapmean.rc
    dcTL["diff"] = []  # satbang tlapmean differs from tlapmean.rc
    dcTL["none"] = []  # no value found in satbang or tlapmean.rc
    dcTL["what"] = []  # value found in satbang but not in tlapmean.rc

    dcBA["add"] = []   # satbang records added
    dcBA["rem"] = []   # satbang records removed
    dcBA["fix"] = []   # statistics from satbang.fixkeys()

    dcBI["add"] = []   # satbias records added
    dcBI["rem"] = []   # satbias records removed
    dcBI["fix"] = []   # statistics from satbias.fixkeys()

    fpath = { satbang: satbang_path, satbias: satbias_path }
    label = { satbang: "satbang", satbias: "satbias" }
    sb_len_0 = {}

    # check satbang and satbias for needed fixes
    #-------------------------------------------
    for (sb, dc) in [(satbang, dcBA), (satbias, dcBI)]:
        sb_len_0[sb] = sb.len()

        dc["fix"] = sb.fixkeys(satinfo.change_list)

        # loop through all records in sb
        #-------------------------------
        for (key, values) in sb.records():
            if satinfo.key_exists(key):

                # check tlap in satbang
                #----------------------
                if sb == satbang:
                    sb_tlap = values["tlapmean"]

                    if tlapmean.key_exists(key):
                        tl_tlap = tlapmean.value(key)

                        if float(sb_tlap) != float(tl_tlap):
                            if float(sb_tlap) == 0.0:
                                sb.change(key, "tlapmean", tl_tlap)
                                dcTL["fill"].append((key, sb_tlap, tl_tlap))
                            else:
                                dcTL["diff"].append((key, sb_tlap, tl_tlap))

                    # no tlap value in tlapmean
                    #--------------------------
                    else:
                        if float(sb_tlap) == 0.0:
                            dcTL["none"].append((key, sb_tlap))
                        else:
                            dcTL["what"].append((key, sb_tlap))

            # remove records not found in satinfo
            #------------------------------------
            else:
                sb.remove_records(key)
                dc["rem"].append(key)

        # add additional records from satinfo
        #------------------------------------
        for (key, values) in satinfo.records():

            if not sb.key_exists(key):
                sb.append_zero_record(key)
                dc["add"].append(key)

                if sb == satbang and tlapmean.key_exists(key):
                    sb.change(key, "tlapmean", tlapmean.value(key))

    # error checks
    #-------------
    for sb in [satbang, satbias]:

        # check for matching records
        #---------------------------
        if set(sb.keys()) != set(satinfo.keys()):
            msg = "Error. Instrument/channel lists differ: satinfo.rc and {0}"
            raise Exception(msg.format(label[sb]))

        # check for duplicate (instrument, channel) keys
        #-----------------------------------------------
        keyList = sb.keys()
        if len(keyList) != len(set(keyList)):
            for key in keyList:
                if keyList.count(key) > 1:
                    msg = "Error. Duplicate key in {0}: {1}"
                    raise Exception(msg.format(label[sb], key))

    # rewrite satbang and satbias files
    #----------------------------------
    for (sb, dc) in [(satbang, dcBA), (satbias, dcBI)]:

        if   sortflag == "satinfo": outkeys = satinfo.keys()
        elif sortflag == "satbang": outkeys = satbang.keys()
        elif sortflag == "satbias": outkeys = satbias.keys()
        elif sortflag == "alpha"  : outkeys = keysort(sb.keys())
        elif sortflag == "default": outkeys = sb.keys()
        else:
            raise Exception("Error. Unknown sortflag: {0}".format(sortflag))

        # determine whether rewrite is necessary
        #---------------------------------------
        rewrite_flag = False

        for action in ["add", "rem"]:
            if dc[action]:
                rewrite_flag = True

        for (inst, new) in dc["fix"].keys():
            if dc["fix"][(inst, new)]:
                rewrite_flag = True

        if sb == satbang and dcTL["fill"]:
            rewrite_flag = True

        if sortflag != "default" and sortflag != label[sb]:
            rewrite_flag = True

        # rewrite file
        #-------------
        if rewrite_flag:
            sb.write(fpath[sb]+".new", outkeys)

        # write summary to file?
        #-----------------------
        basename = os.path.basename(fpath[sb])
        dirname = os.path.dirname(fpath[sb])

        if not(rewrite_flag):
            if sb == satbang and (dcTL["diff"] or dcTL["none"] or dcTL["what"]):
                pass
            else:
                print "No irregularites found in {0}".format(basename)
                continue

        # write summary
        #--------------
        summary_file = os.path.join(dirname, "check_summary."+basename)
        if os.path.isfile(summary_file):
            os.rename(summary_file, summary_file+'~')

        with open(summary_file, mode='w') as sfile:
            sfile.write("\n")
            sfile.write(" =====================\n")
            sfile.write(" check_{0} summary\n".format(label[sb]))
            sfile.write(" =====================\n")
            sfile.write(" current dir: {0}\n\n".format(os.path.abspath(os.curdir)))

            sfile.write(" {0}_rst: {1}\n".format(label[sb], fpath[sb]))
            sfile.write(" satinfo.rc:  {0}\n".format(satinfo_path))

            if sb == satbang:
                sfile.write(" tlapmean.rc: {0}\n\n".format(tlapmean_path))

                msg = "%6i recs filled in missing tlapmean from tlapmean.rc\n"
                sfile.write(msg % (len(dcTL["fill"])))

                msg = "%6i recs with tlapmean different than tlapmean.rc\n"
                sfile.write(msg % (len(dcTL["diff"])))

                msg = "%6i recs with tlapmean equal zero\n"
                sfile.write(msg % (len(dcTL["none"])))

                msg = "%6i recs with tlapmean in satbang but not in tlapmean.rc\n\n"
                sfile.write(msg % (len(dcTL["what"])))
            else:
                sfile.write("\n")

            msg = "%6i recs removed because instrument/channel not in satinfo.rc\n"
            sfile.write(msg % (len(dc["rem"])))

            msg = "%6i recs added from satinfo.rc\n\n"
            sfile.write(msg % (len(dc["add"])))

            # instrument name changes
            #------------------------
            for (inst, new) in sorted(dc["fix"].keys()):
                if type(new) == str:
                    msg = "%6i recs renamed %s to %s\n"
                    sfile.write(msg % (dc["fix"][(inst, new)], inst, new))

            # instrument channel number changes
            #----------------------------------
            for (inst, new) in sorted(dc["fix"].keys()):
                if type(new) == tuple:
                    (old_nums, new_nums) = new
                    msg = "\n%6i recs where %s channels were renumbered\n"
                    msg1 = "       - old: {0}\n"
                    msg2 = "       - new: {0}\n"
                    sfile.write(msg % (dc["fix"][(inst, new)], inst))
                    sfile.write(msg1.format(old_nums))
                    sfile.write(msg2.format(new_nums))

            msg = "\n%6i {0} recs before update\n".format(label[sb])
            sfile.write(msg % (sb_len_0[sb]))

            msg = "%6i {0} recs after update\n".format(label[sb])
            sfile.write(msg % (sb.len()))

            if rewrite_flag:
                msg = "\n new {0} file has been written: {1}\n"
                sfile.write(msg.format(label[sb], sb.filename()))
                sfile.write(' sortflag = "{0}"\n'.format(sortflag))

        # display summary without details to stdout
        #------------------------------------------
        cat([summary_file])
        print(" For details, see file, {0}\n".format(summary_file))

        # write details
        #--------------
        msg1 = "* records where zero-value replaced by value from tlapmean.rc"
        msg2 = "* records where satbang tlapmean differs from tlapmean.rc"
        msg3 = "* records where no value found in satbang or tlapmean.rc"
        msg4 = "* records where value found in satbang but not in tlapmean.rc"

        msg5 = "* records added to agree with satinfo.rc (coeffs = 0.000)"
        msg6 = "* records removed; not found in satinfo.rc"

        with open(summary_file, mode='a') as sfile:

            msg = "The following details are listed below:"
            sfile.write("\n\n"+msg+"\n" + '='*len(msg)+"\n")
            if sb == satbang:
                if dcTL["fill"]: sfile.write(msg1+"\n")
                if dcTL["diff"]: sfile.write(msg2+"\n")
                if dcTL["none"]: sfile.write(msg3+"\n")
                if dcTL["what"]: sfile.write(msg4+"\n")

            if dc["add"]: sfile.write(msg5+"\n")
            if dc["rem"]: sfile.write(msg6+"\n")

            if dcTL["fill"] and sb == satbang:
                sfile.write("\n\n"+msg1+"\n" + '-'*len(msg1)+"\n")
                msg = "%-25s tlapmean:  %13s => %13s\n"
                for (key, sb_tlap, tl_tlap) in dcTL["fill"]:
                    sfile.write(msg %(str(key), sb_tlap, tl_tlap))

            if dcTL["diff"] and sb == satbang:
                sfile.write("\n\n"+msg2+"\n" + '-'*len(msg2)+"\n")
                headers = ("instrument/channel", "satbang_rst",
                          "tlapmean.rc", "difference")
                underscores = tuple('-'*len(hh) for hh in headers)
                sfile.write("%-25s %13s %13s %13s\n" % headers)
                sfile.write("%-25s %13s %13s %13s\n" % underscores)
                
                msg = "%-25s %13s %13s %13.6e\n"
                for (key, sb_tlap, tl_tlap) in dcTL["diff"]:
                    diff = float(sb_tlap) - float(tl_tlap)
                    sfile.write(msg %(str(key), sb_tlap, tl_tlap, diff))

            if dcTL["none"] and sb == satbang:
                sfile.write("\n\n"+msg3+"\n" + '-'*len(msg3)+"\n")
                msg = "%-25s tlapmean: %13s\n"
                for (key, sb_tlap) in dcTL["none"]:
                    sfile.write(msg %(str(key), tl_tlap))
    
            if dcTL["what"] and sb == satbang:
                sfile.write("\n\n"+msg4+"\n" + '-'*len(msg4)+"\n")
                msg = "%-25s tlapmean in satbang_rst: %13s\n"
                for (key, sb_tlap) in dcTL["what"]:
                    sfile.write(msg %(str(key), sb_tlap))
    
            if dc["add"]:
                sfile.write("\n\n"+msg5+"\n" + '-'*len(msg5)+"\n")
                for key in keysort(dc["add"]):
                    sfile.write(str(key)+"\n")
    
            if dc["rem"]:
                sfile.write("\n\n"+msg6+"\n" + '-'*len(msg6)+"\n")
                for key in keysort(dc["rem"]):
                    sfile.write(str(key)+"\n")
            
#.......................................................................
def keysort(keylist):
    """
    Sort a list of (instrument, channel_number) keys where the channel
    numbers are sorted numerically rather than alphabetically.
    """
    newlist = [(inst, int(chan)) for (inst, chan) in keylist]
    newlist.sort()
    return [(inst, str(chan)) for (inst, chan) in newlist]

#.......................................................................
def cat(file_arr):
    "concatenate files and print on standard output"
    for file in file_arr:
        with open(file, "r") as f:
            print(f.read())

#.......................................................................
if __name__ == "__main__":
    check_files(*sys.argv[1:])
