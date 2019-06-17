import os

class SatInfoRc(object):
    "class to create instances of gmao_global_satinfo.rc file"

    #.......................................................................
    def trange(start, stop, step=1):
        "This is a range function that returns a tuple instead of a list."
        return tuple(range(start, stop, step))


    # list of instrument name changes, (old_name, new_name), and
    # instrument channel number changes, (inst, ([old_nums], [new_nums]))
    #------------------------------------------------------------------
    change_history = [ ("airs281SUBSET",   "airs_aqua"),
                       ("airs281_aqua",    "airs_aqua"),
                       ("cris399_npp",     "cris_npp"),
                       ("iasi616_metop-a", "iasi_metop-a"),
                       ("iasi616_metop-b", "iasi_metop-b"),
                       ("seviri_m08", (trange(1,9), trange(4,12))),
                       ("seviri_m09", (trange(1,9), trange(4,12))),
                       ("seviri_m10", (trange(1,9), trange(4,12))) ]

    #.......................................................................
    def __init__(self, filename):
        "Initialize instance."
        if not os.path.isfile(filename):
            raise ValueError("file, {0}, does not exist".format(filename))

        self.change_list = []
        self._filename = filename
        self._data = []

        self._read()

    #.......................................................................
    def _read(self):
        "Read input from file into self._data array."

        channel_nums = {}

        with open(self._filename, mode='r') as satinfo_rcfile:

            # read records from file
            #-----------------------
            for line in satinfo_rcfile:
                if line.strip()    == "" : continue
                if line.strip()[0] == '!': continue

                items = line.split()
                if len(items) != 9:
                    next

                # extract key info from record
                #-----------------------------
                instrument = items[0]
                channel = items[1]
                key = (instrument, channel)

                if key in dict(self._data):
                    msg = "Duplicate record for {0} in {1}"
                    raise ValueError(msg.format(key, self._filename))

                # check for instrument in change_history
                #---------------------------------------
                for (old, new) in self.change_history:

                    # need to change (old => new) or (new => old)?
                    #---------------------------------------------
                    for (n1, n2) in [(old, new), (new, old)]:

                        if instrument == n2 and type(n1) == str:
                            if (n1, n2) not in self.change_list:
                                if n1 in dict(self.change_list):
                                    msg = "Error. Cannot rename {0} to {1} and {2}."
                                    nX = dict(self.change_list)[n1]
                                    raise IOError(msg.format(n1, n2, nX))
                                self.change_list.append((n1, n2))

                            if (n2, n1) in self.change_list:
                                msg = "Error. Both {0} and {1} found in {2}."
                                raise IOError(msg.format(n1, n2, self._filename))

                    # collect channel number information
                    #-----------------------------------
                    if instrument == old and type(new) == tuple:
                        if instrument not in channel_nums:
                            channel_nums[instrument] = set()

                        channel_nums[instrument].add(int(channel))

                # extract values from record
                #---------------------------
                values = {"iuse"     : items[2],
                          "error"    : items[3],
                          "error_cld": items[4],
                          "ermax"    : items[5],
                          "var_b"    : items[6],
                          "var_pg"   : items[7],
                          "icld_det" : items[8]}

                self._data.append((key, values))

            # check for channel number changes to add to change_list
            #-------------------------------------------------------
            for instr in channel_nums:
                chans = (channel_nums[instr])
                (old, new) = dict(self.change_history)[instr]

                if chans == set(new):
                    self.change_list.append((instr, (old, new)))

                elif chans == set(old):
                    self.change_list.append((instr, (new, old)))

                # NOTE: differences in order of numbers is not a problem
                #-------------------------------------------------------
                else:
                    msg = "{0} channel numbers should be {1} or {2}: {3}"
                    raise IOError(msg.format(instr, old, new, chans))

    #.......................................................................
    def filename(self):
        "Return filename associated with instance."
        return self._filename

    #.......................................................................
    def has_instrument(self, instr):
        """
        return True if instr is found in self._data array; otherwise return False.

        input parameter:
        => instr: instrument name to search
        """
        if instr in dict(self.keys()):
            return True
        else:
            return False

    #.......................................................................
    def key_exists(self, key):
        """
        return True if key is found in self._data array; otherwise return False.

        input parameter:
        => key: (instrument, channel) tuple of string values
        """
        return dict(self._data).has_key(key)

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
        => values is a dictionary with keys, {"iuse", "error", "error_cld",
           "ermax", "var_b", "var_pg", "icld_det"}
        """
        for (key, values) in self._data:
            yield (key, values)
            
    #.......................................................................
    def write(self, newname="default"):
        """
        Write self._data array to output file.

        input parameter:
        => newname: name of output file
        """
        if newname == "default":
            newname = self._filename

        with open(newname, mode='w') as new:
            new.write("!sensor/instr/sat   chan iuse   error error_cld "+
                      "ermax  var_b   var_pg  icld_det\n")

            for ((instrument, channel), vals) in self._data:
                str_arr = [instrument, channel]
                str_arr.append(vals["iuse"])
                str_arr.append(vals["error"])
                str_arr.append(vals["error_cld"])
                str_arr.append(vals["ermax"])
                str_arr.append(vals["var_b"])
                str_arr.append(vals["var_pg"])
                str_arr.append(vals["icld_det"])

                new.write((" %-16s%7s%5s%8s%8s%8s%8s%8s%7s\n")%(tuple(str_arr)))

        self._filename = newname
