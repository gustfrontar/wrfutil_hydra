import os
import re
import time
import datetime
import operator

import numpy as np

# TODO: TITLE, SEQUENTIAL?


NUMBER = '[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'

print(NUMBER)

class CTLReader(object):
    def __init__(self, filename):
        self.dimensions = {}
        self.variables = {}

        self.filename = filename
        with open(filename, 'rb') as fp:
            self.ctl = fp.read().decode('utf-8')

            print(self.ctl)

        #self._read_data()    #Esto lee los datos, por el momento lo comento porque no es compatible con templates. 
        #Habra que desarrollar esa opcion.
        #self._read_dimensions()
        self._read_vars()

    def _read_data(self):
        self.undef = eval(re.search("undef (%s)" % NUMBER, self.ctl).group(1))
        self.yrev = bool(re.search("options.*yrev", self.ctl))
        big_endian = bool(re.search("optins.*big_endian", self.ctl))
        dset = re.search("dset (.*)", self.ctl).group(1)
        if dset.startswith('^'):
            dset = os.path.join(os.path.dirname(self.filename), dset[1:])
        data = np.fromfile(dset, 'f')
        if big_endian:
            data = data.byteswap()
        self.data = np.ma.masked_values(data, self.undef)
    
    def _read_dimensions(self):
        if 'xdef' in self.ctl:
            self.variables['longitude'] = Variable('longitude', self._parse_dimension('XDEF'))
            self.dimensions['longitude'] = len(self.variables['longitude'])
        if 'ydef' in self.ctl:
            self.variables['latitude'] = Variable('latitude', self._parse_dimension('YDEF'))
            self.dimensions['latitude'] = len(self.variables['latitude'])
        if 'zdef' in self.ctl:
            self.variables['levels'] = Variable('levels', self._parse_dimension('ZDEF'))
            self.dimensions['levels'] = len(self.variables['levels'])
        if 'tdef' in self.ctl:
            self.variables['time'] = Variable('time', self._parse_dimension('TDEF'))
            self.dimensions['time'] = len(self.variables['time'])

    def _read_vars(self):
        read = False
        for line in self.ctl.split('\n'):
            if line.startswith('ENDVARS'):
                read = False
            if read:
                p = re.compile('(\w+)\s+(\d+)\s+(\d+)\s+(.*)\((.*)\)')
                m = p.match(line)
                name = m.group(1)
                var = self.variables[name] = Variable(name)
                levels = map(int, m.group(2).split(','))
                SPACE = self.dimensions['latitude'] * self.dimensions['longitude']
                if levels[0] > 0:
                    var.dimensions = ('time', 'levels', 'latitude', 'longitude')
                    size = self.dimensions['time'] * self.dimensions['levels'] * (SPACE+2)  # account for header bytes
                else:
                    var.dimensions = ('time', 'latitude', 'longitude')
                    size = self.dimensions['time'] * (SPACE+2)  # account for header bytes

                var.shape = tuple(self.dimensions[dim] for dim in var.dimensions)
                #var.data = self.data[i:i+size].reshape(-1, SPACE+2)[:,1:-1].reshape(var.shape)  # remove header bytes
                #if self.yrev:
                #    var.data = var.data[...,::-1,:]
                #i += size

                units = int(m.group(3))
                if units != 99:
                    raise NotImplementedError('Only unit 99 implemented!')

                var.attributes = {
                    'long_name': m.group(4).strip(),
                    'units': m.group(5).strip(),
                }
            if line.startswith('VAR'):
                i = 0
                read = True

    def _parse_dimension(self, dim):
        p = re.compile("%s\s+(\d+)\s+LINEAR\s+(%s)\s+(%s)" % (dim, NUMBER, NUMBER))
        m = p.search(self.ctl)
        if m:
            length = int(m.group(1))
            start = float(m.group(2))
            increment = float(m.group(3))
            return np.arange(start, start+length*increment, increment)

        p = re.compile("%s\s+\d+\s+LEVELS((\s+%s)+)" % (dim, NUMBER))
        m = p.search(self.ctl)
        if m:
            return np.fromstring(m.group(1), sep=' ')

        p = re.compile("%s\s+(\d+)\s+LINEAR\s+([:\w]+)\s+(\d{1,2})(\w{2})" % dim)
        m = p.search(self.ctl)
        if m:
            length = int(m.group(1))
            start = parse_date(m.group(2))
            increment = parse_delta(m.group(3), m.group(4))
            return np.array([ start+i*increment for i in range(length)]) 


class Variable(object):
    def __init__(self, name, data=None):
        self.name = name
        self.data = data

    def __getitem__(self, index):
        return self.data[index] 

    def __getattr__(self, key):
        return self.attributes[key]

    def __len__(self):
        return len(self.data)

def parse_date(s):
    DATE = re.compile("""
        (?:(?P<hour>\d\d))?     # hh, default 00
        (?::(?P<minute>\d\d))?  # mm, default 00
        Z?
        (?P<day>\d\d)?          # dd, default 1
        (?P<month>\w\w\w)       # 3 char month
        (?P<year>\d\d(?:\d\d)?) # yyyy or 1950 < yy < 2049
    """, re.VERBOSE)
    d = DATE.match(s).groupdict()
    if d['hour'] is None:
        hour = 0
    else:
        hour = int(d['hour'])
    if d['minute'] is None:
        minute = 0
    else:
        minute = int(d['minute'])
    if d['day'] is None:
        day = 1
    else:
        day = int(d['day'])
    month = time.strptime(d['month'], '%b')[1]
    if len(d['year']) == 4:
        year = int(d['year'])
    else:
        year = 1950 + int(d['year'])
    return datetime.datetime(year, month, day, hour, minute)


def parse_delta(value, unit):
    value = int(value)
    if unit.lower() == 'mn':
        return datetime.timedelta(minutes=value)
    if unit.lower() == 'hr':
        return datetime.timedelta(hours=value)
    if unit.lower() == 'dy':
        return datetime.timedelta(days=value)
    if unit.lower() == 'mo':
        raise NotImplementedError('Need to implement month time step')
    if unit.lower() == 'yr':
        raise NotImplementedError('Need to implement year time step')


if __name__ == '__main__':
    import sys

    f = CTLReader(sys.argv[1])
    print(f.dimensions)
    print(f.variables['latitude'])
    #print(np.max(f.variables['TEMP'][:]))
    #print(np.min(f.variables['TEMP'][:]))
