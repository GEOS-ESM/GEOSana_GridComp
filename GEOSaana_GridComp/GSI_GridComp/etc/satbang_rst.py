import glob, os, re, sys

class SatBangRst(object):

    #.......................................................................
    def __init__(self, filename, debug=False, num_coeffs=90):
        "Initialize instance."
        if not os.path.isfile(filename):
            raise Exception("{0} does not exist".format(filename))
        self._filename = filename
        self._debug = debug
        self._data = []
        self._num_coeffs = num_coeffs
        self._read()

    #.......................................................................
    def _read(self):
        "Read input from file into self._data array."
        recnum = None
        with open(self._filename, mode='r') as satbangfile:
            for line in satbangfile:
                if line.strip()    == ""  : continue  # skip blank lines
                if line.strip()[0] == "#" : continue  # skip commented lines
                items = line.split()

                # add previous record to list and then load new info
                #---------------------------------------------------
                if len(items)==4 and items[0].isdigit():

                    if recnum:
                        key = (instrument, channel)
                        values = { "recnum"   : recnum,
                                   "tlapmean" : tlapmean,
                                   "coeffs"   : coeffs }
                        self._append_record((key, values))

                    coeffs = []
                    (recnum, instrument, channel, tlapmean) = items[:]

                # angle correction coefficients
                #------------------------------
                elif len(items) > 0:
                    if not recnum:
                        msg = "IO Error; malformatted first line of record: {0}"
                        raise Exception(msg.format(self._filename))

                    line = line.replace("-", " -").replace("E -", "E-")
                    items = line.split()
                    coeffs.extend(items)

            # add final record to list
            #-------------------------
            if recnum:
                key = (instrument, channel)
                values = { "recnum"   : recnum,
                           "tlapmean" : tlapmean,
                           "coeffs"   : coeffs }

                self._append_record((key, values))

            # check indices
            #--------------
            n = 0
            for (key, values) in self._data:
                recnum = values["recnum"]
                n += 1

                if n != int(recnum):
                    if self._debug:
                        print("satbang[{}] = {}".format(n, self._data[n]))
                    raise Exception("inconsistent recnum (%d): %s" % (n, recnum))

    #.......................................................................
    def _append_record(self, args):
        "Append record to self._data"

        (key, values) = args

        # check for duplicate keys
        #-------------------------
        if key in dict(self._data):
            msg = "Duplicate record for {0} in {1}"
            raise Exception(msg.format(key, self._filename))

        # check for correct number of coeffs
        #-----------------------------------
        coeffs = values["coeffs"]
        if len(coeffs) != self._num_coeffs:
            msg = "{0} has {1} coeffs; should have {2}: {3}"
            params = (key, len(coeffs), self._num_coeffs, self._filename)
            raise Exception(msg.format(*params))

        self._data.append((key, values))

    #.......................................................................
    def append_zero_record(self, key, tlapmean="0.000000E+00"):
        """
        Append new record to self._data array with values contained in key
        and tlapmean and with all coeffs set to 0.000.

        input parameters:
        => key: (instrument, channel) tuple of string values
        => tlapmean: tlapmean string value associated with key
        """
        if key in dict(self._data):
            msg = "key, {0}, already in {1}"
            raise Exception(msg.format(key, self._filename))

        values = {}
        values["recnum"] = len(self._data) + 1
        values["tlapmean"] = tlapmean
        values["coeffs"] = ["0.000"] * self._num_coeffs

        self._append_record((key, values))

    #.......................................................................
    def change(self, key, variable, newvalue, noDupsFlg=True):
        """
        Change value(s) in the self._data record identified by key.

        input parameters:
        => key: (instrument, channel), string tuple key into self._data
        => variable: variable in self._data record
        => newvalue: value to assign to variable
        => noDupsFlg: do not allow duplicate records if True

        Notes:
        1. Acceptable variable are "instrument", "channel", "key", or
           any key found in the values dictionary, e.g. "tlapmean" or "coeffs",
           except "recnum"
        2. When giving a newvalue for "coeffs", a complete list of coeffs
           must be supplied, not just the coeffs which are being changed.
        """
        values = dict(self._data)[key]
        index = self._data.index((key, values))

        if variable in ["instrument", "channel", "key"]:
            (instrument, channel) = key

            if (variable == "instrument"):
                newkey = (newvalue, channel)

            elif (variable == "channel"):
                newkey = (instrument, newvalue)

            elif (variable == "key"):
                newkey = newvalue

            if noDupsFlg and newkey != key and newkey in dict(self._data):
                msg = "Error. newkey already in self._data: {0}, {1}"
                raise Exception(msg.format(key, newkey))

        elif variable in values:
            if variable == "recnum":
                msg = "Error. Cannot change the recnum value in a record."
                raise Exception(msg)

            if variable == "coeffs":
                if len(newvalue) != self._num_coeffs:
                    msg = "Error. Coeffs has {0} values instead of {1}"
                    raise Exception(msg.format(len(newvalue), self._num_coeffs))
                values["coeffs"] = newvalue

            newkey = key
            values[variable] = newvalue

        else:
            msg = "Error. Variable [{0}] was not found in file records: {1}"
            raise Exception(msg.format(variable, self.filename()))

        # replace record
        #---------------
        self._data[index] = (newkey, values)

    #.......................................................................
    def filename(self):
        "Return filename associated with instance."
        return self._filename

    #.......................................................................
    def fixkeys(self, change_list):
        """
        Make instrument name and channel number adjustments specified in an
        inputted list of changes.

        For channel number changes, first check that the set of channels for
        a specified instrument matches the set in the change list.

        input parameter
        => change_list: list of instrument name changes, (old_name, new_name),
              and instr channel number changes, (inst, ([old_nums], [new_nums])

        return value
        => fix_stats: number of fixes for each item in change_list
        """

        # check change_list for instruments with channel number changes
        #--------------------------------------------------------------
        inst_list = [inst for (inst, tp) in change_list if type(tp) is tuple]
        fix_chans = {}

        for inst in inst_list:
            old_nums = dict(change_list)[inst][0]
            new_nums = dict(change_list)[inst][1]

            found = [int(chan) for (ii, chan) in dict(self._data).keys()
                     if ii == inst]
            found.sort()

            if found == sorted(old_nums):
                fix_chans[inst] = True

            elif found == sorted(new_nums):
                fix_chans[inst] = False

            else:
                msg = "Error. {0} channel list does not equal {1} or {2}: {3}"
                raise Exception(msg.format(inst, old_nums, new_nums, found))
                
        # initialize fix statistics to zero
        #----------------------------------
        fix_stats = {}
        for (old, new) in change_list:
            fix_stats[(old, new)] = 0

        # fix keys
        #---------
        for (key, values) in self._data:
            inst = key[0]
            chan = key[1]

            if inst in dict(change_list):
                new = dict(change_list)[inst]

                if type(new) is str:
                    self.change(key, "instrument", new)
                    fix_stats[(inst, new)] += 1

                elif type(new) is tuple:
                    if not fix_chans[inst]: continue

                    (old_nums, new_nums) = new
                    index = old_nums.index(int(chan))
                    self.change(key, "channel", str(new_nums[index]), False)
                    fix_stats[(inst, new)] += 1

                else:
                    msg = "Error. change_list record for {0} contains {1}: {2}."
                    raise Exception(msg.format(inst, type(new), new))

        return fix_stats

    #.......................................................................
    def key_exists(self, key):
        """
        return True if key is found in self._data array; otherwise return False.

        input parameter:
        => key: (instrument, channel) tuple of string values
        """
        return dict(self._data).__contains__(key)

    #.......................................................................
    def key_and_values(self, index):
        """
        Return (key, values) tuple for index in the self._data array.
        key = (instrument, channel) tuple
        values = {"recnum":recnum, "tlapmean": tlapmean, "coeffs":coeffs[0:89]}
        """
        return self._data[index]

    #.......................................................................
    def keys(self):
        "Return list of (instrument, channel) values in self._data array."
        return [key for (key, values) in self._data]

    #.......................................................................
    def len(self):
        "Return the number of records in self._data array."
        return len(self._data)

    #.......................................................................
    def records(self):
        """
        Create generator to iterate through records of self._data array.

        Each record is a (key, values) tuple where
        => key is a tuple with two string values, (instrument, channel)
        => values is a dictionary with keys, {"recnum", "tlapmean", "coeffs"}
        => the value of "coeffs" is an array of 90 values
        """
        for (key, values) in self._data[:]:
            yield (key, values)

    #.......................................................................
    def remove_records(self, keylist):
        """
        Remove records from self._data array.

        input parameter:
        => keylist: key or list of keys identifying records to remove from data
                    where each key is (instrument, channel) tuple of str values
        """

        # remove records
        #---------------
        if type(keylist) == tuple:
            keylist = [keylist]

        for key in keylist:
            values = dict(self._data)[key]
            self._data.remove((key, values))

        # reset recnum values in self._data array
        #-----------------------------------------
        for index in range(len(self._data)):
            (key, values) = self._data[index]
            if values["recnum"] != str(index+1):
                values["recnum"] = str(index+1)
                self._data[index] = (key, values)

    #.......................................................................
    def write(self, newname="default", outkeys="default"):
        """
        Write data from self._data array to output file.

        input parameters:
        => newname: name of output file
        => outkeys: list of (instrument, channel) keys identifying
                 records, and order of records, to write to output
        """
        if newname == "default":
            newname = self._filename

        if outkeys == "default":
            outkeys = self.keys()

        recnum = 0
        with open(newname, mode='w') as new:
            for key in outkeys:
                values = dict(self._data)[key]

                num_coeffs = len(values["coeffs"])
                if num_coeffs != self._num_coeffs:
                    msg = "Error: {0} has {1} coeffs instead of {2}"
                    raise Exception(msg.format(key, num_coeffs, self._num_coeffs))

                recnum += 1
                str_arr = [recnum]
                str_arr.extend(key)
                str_arr.append(values["tlapmean"])
                str_arr.extend(values["coeffs"])

                new.write(("%5s %-20s%6s%15s\n")  % (tuple(str_arr[0:4])))
                new.write(("    "+"%7s"*10 +"\n") % (tuple(str_arr[4:14])))
                new.write(("    "+"%7s"*10 +"\n") % (tuple(str_arr[14:24])))
                new.write(("    "+"%7s"*10 +"\n") % (tuple(str_arr[24:34])))
                new.write(("    "+"%7s"*10 +"\n") % (tuple(str_arr[34:44])))
                new.write(("    "+"%7s"*10 +"\n") % (tuple(str_arr[44:54])))
                new.write(("    "+"%7s"*10 +"\n") % (tuple(str_arr[54:64])))
                new.write(("    "+"%7s"*10 +"\n") % (tuple(str_arr[64:74])))
                new.write(("    "+"%7s"*10 +"\n") % (tuple(str_arr[74:84])))
                new.write(("    "+"%7s"*10 +"\n") % (tuple(str_arr[84:94])))
                new.write("\n")

        self._filename = newname

#.......................................................................
def satbang_dir_and_name(pathname):
    """
    Call the dir_and_file_names() routine with the pattern
    necessary to identify the dirname and filename a satbang file.
    """
    (dirname, filename) = dir_and_file_names(pathname, "*.ana_satbang_rst.*.txt")
    return (dirname, filename)

#.......................................................................
def dir_and_file_names(pathname, pattern=False):
    """
    Return the directory location and filename from a pathname and pattern.

    input parameters
    ----------------
    => pathname: name of a file, potentially including its directory path, or
                 name of a directory which includes a file with a specified
                 pattern as part of its name
    => pattern: the pattern used to identify a file, if pathname is a directory

    returns (dirname, filename)
    """

    # if pathname is directory
    #-------------------------
    if os.path.isdir(pathname) and pattern:
        dirname = re.sub(r'/$', '', pathname.strip())

        # look for files that match the pattern
        #--------------------------------------
        filelist = glob.glob(os.path.join(dirname, pattern))
        exclude = re.compile(r"^check_summary")

        for fpath in filelist:
            bname = os.path.basename(fpath)
            if exclude.match(bname):
                filelist.remove(fpath)

        # only one file found
        #--------------------
        if len(filelist) == 1:
            filename = os.path.basename(filelist[0])

        # multiple files found; choose one
        #---------------------------------
        elif len(filelist) > 1:
            filename = None
            while filename == None:
                print("\nSelect file to check")
                nnn = 1

                for sb in filelist:
                    print("{0}: {1}".format(nnn, sb))
                    nnn += 1

                sys.stdout.write("Make selection: ")
                index = int(sys.stdin.readline())
                if index in range(1, len(filelist)+1):
                    filename = filelist[int(index) - 1]

        # no files found
        #---------------
        else:
            msg = "Error. No files found in {0} with pattern {1}."
            raise Warning(msg.format(dirname, pattern))

    # if pathname is file
    #--------------------
    else:
        if not os.path.isfile(pathname):
            msg = "file not found for pathname = {0} and pattern = {1}"
            raise Warning(msg.format(pathname, pattern))

        filename = os.path.basename(pathname)
        dirname = os.path.dirname(pathname)
        if dirname == "":
            dirname = os.curdir

    return (dirname, filename)
