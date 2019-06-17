import glob, os, re
from satbang_rst import SatBangRst, dir_and_file_names

class SatBiasRst(SatBangRst):

    #.......................................................................
    def __init__(self, filename, debug=False):
        "Initialize instance."
        num_coeffs = 7
        super(SatBiasRst, self).__init__(filename, debug, num_coeffs)

    #.......................................................................
    def _read(self):
        "Read input from file into self._data array."
        recnum = None
        with open(self._filename, mode='r') as satbiasfile:
            for line in satbiasfile:
                if line.strip()    == ""  : continue  # skip blank lines
                if line.strip()[0] == "#" : continue  # skip commented lines
                line = re.sub("(\d)-(\d)", r"\1 -\2", line)

                items = line.split()

                # check line format
                #------------------
                if len(items) != 10:
                    msg = "{0} items on line instead of 10: {1} "
                    raise Exception(msg.format(len(items), line))

                if not items[0].isdigit():
                    msg = "non-numeric record number [{0}] in {1}"
                    raise Exception(msg.format(items[0], self._filename))

                # extract values from line
                #-------------------------
                (recnum, instrument, channel) = items[:3]
                coeffs = items[3:]

                key = (instrument, channel)
                values = { "recnum" : recnum, "coeffs" : coeffs }

                # append data
                #------------
                self._append_record((key, values))

        # check indices
        #--------------
        n = 0
        for (key, values) in self._data:
            recnum = values["recnum"]
            n += 1

            if n != int(recnum):
                raise Exception("inconsistent recnum (%d): %s" % (n, recnum))

    #.......................................................................
    def append_zero_record(self, key):
        """
        Append new record to self._data array with values contained in key
        and with all coeffs set to 0.000.

        input parameters:
        => key: (instrument, channel) tuple of string values
        """
        if key in dict(self._data):
            msg = "key, {0}, already in {1}"
            raise Exception(msg.format(key, self._filename))

        values = {}
        values["recnum"] = len(self._data) + 1
        values["coeffs"] = ["0.000000"] * self._num_coeffs

        self._data.append((key, values))

    #.......................................................................
    def write(self, newname="default", keys="default"):
        """
        Write data from self._data array to output file.

        input parameter:
        => newname: name of output file to which to write data
        => keys: list of (instrument, channel) keys identifying
                 records, and order of records, to write to output
        """
        if newname == "default":
            newname = self._filename

        if keys == "default":
            keys = self.keys()

        recnum = 0
        with open(newname, mode='w') as new:
            for key in keys:
                values = dict(self._data)[key]

                num_coeffs = len(values["coeffs"])
                if num_coeffs != self._num_coeffs:
                    msg = "Error: {0} had {1} coeffs; should be {2}"
                    raise Exception(msg.format(key, num_coeffs, self._num_coeffs))

                recnum += 1
                str_arr = [recnum]
                str_arr.extend(key)
                str_arr.extend(values["coeffs"])
                new.write(("%5s %-20s%6s"+"%12s"*7+"\n") % (tuple(str_arr)))

        self._filename = newname

#.......................................................................
def satbias_dir_and_name(pathname):
    """
    Call the dir_and_file_names() routine with the pattern
    necessary to identify the dirname and filename of a satbias file.
    """
    (dirname, filename) = dir_and_file_names(pathname, "*.ana_satbias_rst.*.txt")
    return (dirname, filename)
