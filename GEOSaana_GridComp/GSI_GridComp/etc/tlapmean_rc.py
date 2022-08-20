import os

class TLapMeanRc(object):
    "class to create instances of gmao_global_tlapmean.rc file"

    #.......................................................................
    def __init__(self, filename):
        "Initialize instance."
        if not os.path.isfile(filename):
            raise ValueError("file, {0}, does not exist".format(filename))
        self._filename = filename
        self._data = []
        self._read()

    #.......................................................................
    def _read(self):
        "Read input from file into self._data array."
        with open(self._filename, mode='r') as tlapmean_rcfile:
            for line in tlapmean_rcfile:
                if not line or line.isspace():  continue
                items = line.split()
                if items[0][0] == '!':  continue

                if len(items) == 3:
                    instrument = items[0]
                    channel = items[1]

                    key = (instrument, channel)
                    value = items[2]

                    self._data.append((key, value))

    #.......................................................................
    def append(self, key, value):
        "Append record to self._data"
        self._data.append((key, value))

    #.......................................................................
    def filename(self):
        return self._filename

    #.......................................................................
    def key_exists(self, key):
        return dict(self._data).__contains__(key)

    #.......................................................................
    def value(self, key):
        return dict(self._data)[key]

    #.......................................................................
    def write(self, newname):
        with open(newname, mode='w') as new:
            for ((instrument, channel), value) in self._data:
                new.write(" %-20s%6s%15s\n"%(instrument, channel, value))
