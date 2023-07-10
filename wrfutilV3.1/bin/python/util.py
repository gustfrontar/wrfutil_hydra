# --------- CONFIGURATION FUNCTIONS ---------- #
def load_conf():
   '''
   Load configuration dictionaries containing:
     - var_info: information of assimilated variables (ID number, description, units)
     - obstype_info: information of obsevation type (ID number, description, variables)
   '''

   var_info = {'u'  : {'id':2819,  'name': {'long': 'Viento zonal', 'short': 'U'}, 'units': r'm s$^{-1}$'},
               'v'  : {'id':2820,  'name': {'long': 'Viento meridional', 'short': 'V'}, 'units': r'm s$^{-1}$'},
               't'  : {'id':3073,  'name': {'long': 'Temperatura', 'short': 'T'}, 'units': 'K'},
               'tv' : {'id':3079,  'name': {'long': '', 'short': 'TV'}, 'units': ''},
               'q'  : {'id':3330,  'name': {'long': 'Humedad especifica', 'short': 'Q'}, 'units': r'g kg$^{-1}$'},
               'rh' : {'id':3331,  'name': {'long': 'Humedad relativa', 'short': 'HR'}, 'units': '%'},
               'ref': {'id':4001,  'name': {'long': 'Reflectividad', 'short': 'Ref'}, 'units': 'dBZ'},
               'dv' : {'id':4002,  'name': {'long': 'Velocidad Doppler', 'short': 'DV'}, 'units': r'm s$^{-1}$'},
               'ps' : {'id':14593, 'name': {'long': 'Presion de superficie', 'short': 'PSFC'}, 'units': 'hPa'},
               'us' : {'id':82819, 'name': {'long': 'Viento zonal', 'short': 'U'}, 'units': r'm s$^{-1}$'},
               'vs' : {'id':82820, 'name': {'long': 'Viento meridional', 'short': 'V'}, 'units': r'm s$^{-1}$'},
               'ts' : {'id':83073, 'name': {'long': 'Temperatura', 'short': 'T'}, 'units': 'K'},
               'qs' : {'id':83330, 'name': {'long': 'Humedad especifica', 'short': 'Q'}, 'units': r'g kg$^{-1}$'},
               'rhs': {'id':83331, 'name': {'long': 'Humedad relativa', 'short': 'HR'}, 'units': '%'}
              }

   obstype_info = {'ADPUPA': {'id':1,  'name': {'long': 'Sondeos', 'short': 'Sondeos'}, 'vars':['u', 'v', 't', 'rh']},
                   'AIRCFT': {'id':3,  'name': {'long': 'Aviones', 'short': 'Aviones'}, 'vars':['u', 'v', 't']},
                   'SATWND': {'id':4,  'name': {'long': 'GOES', 'short': 'GOES'}, 'vars':['u', 'v']},
                   'ADPSFC': {'id':8,  'name': {'long': 'Est. Superficie', 'short': 'Superficie'}, 'vars':['ps', 'us', 'vs', 'ts', 'rhs']},
                   'SFCSHP': {'id':9,  'name': {'long': 'Barcos', 'short': 'Barcos'}, 'vars':['ps', 'us', 'vs', 'ts']},  #CHEQUEAR
                   'RADAR' : {'id':12, 'name': {'long': 'Radar', 'short': 'Radar'}, 'vars':['ref', 'dv']},
                   'ASCATW': {'id':20, 'name': {'long': 'ASCAT', 'short': 'ASCATW'}, 'vars':[]},
                   'AIRS'  : {'id':21, 'name': {'long': 'Sat. Polar', 'short': 'Sat. Polar'}, 'vars':['t', 'q']},
                   'ADPAUT': {'id':22, 'name': {'long': 'Est. Automatica', 'short': 'Automatica'}, 'vars':['ps', 'us', 'vs', 'ts', 'rhs']}
                   }

   return var_info, obstype_info

# --------- GET FUNCTIONS ---------- #
def get_key(d, subkey, value):
   out = [name for name, props in d.items() if props[subkey] == value]
   if len(out) == 1:
      return out[0]
   else:
      return out

# --------- DATE FUNCTIONS --------- #
def str2date(string, fmt='%Y%m%d_%H%M00'):
   from datetime import datetime as dt
   return dt.strptime(string, fmt)

def date2str(date, fmt='%Y%m%d_%H%M00'):
   ''' Date: a datetime object. '''
   from  datetime import datetime as dt
   return dt.strftime(date, fmt)

def get_dates(initial_time, nt, freq, fmt='%Y%m%d_%H%M00'):
   '''Create a list of datetime objects with 'nt' times back from 'initial_time' '''
   from datetime import timedelta

   enddate = str2date(initial_time, fmt)
   delta = timedelta(seconds=freq*3600)
   inidate = (enddate - nt//freq*delta)

   dates = []
   for date in datespan(inidate, enddate+delta, delta):
      dates.append(date)
   return dates

def datespan(startDate, endDate, delta):
   '''
   La funcion devuelve un "generator" que contiene un objecto date
   Input:
       starDate (objeto): de la clase datetime que indica la fecha inicial
       endDate (objeto): de la clase datetime que indica la fecha final
       delta (objeto): de la clase datetime que indica el intervalo temporal
   '''
   currentDate = startDate
   while currentDate < endDate:
      yield currentDate
      currentDate += delta

# --------- I/O FUNCTIONS ---------#
def read_obsdat(filename):
   '''
   Read obs.dat file and create a pandas dataframe
     wk(1) = REAL(elem(n),r_sngl)
     wk(2) = REAL(rlon(n),r_sngl)
     wk(3) = REAL(rlat(n),r_sngl)
     wk(4) = REAL(rlev(n),r_sngl)
     wk(5) = REAL(odat(n),r_sngl)
     wk(6) = REAL(oerr(n),r_sngl)
     wk(7) = REAL(otyp(n),r_sngl)
     wk(8) = REAL(omb(n),r_sngl)
     wk(9) = REAL(oma(n),r_sngl)
     wk(10) = REAL(slot(n),r_sngl)
   '''
   import pandas as pd
   import numpy as np

   # Create a dtype array with the binary data format and the desired column names
   dtype='float32'
   dt = np.dtype([('HEADER1', dtype), ('VAR', dtype), ('LON', dtype), ('LAT', dtype), ('LEV', dtype), ('OBS', dtype), ('ERROR', dtype), ('TYPE', dtype), ('OMB', dtype), ('OMA', dtype), ('SLOT', dtype), ('HEADER2', dtype)])

   '''
   try:
      # Read binary file
      data = np.fromfile(filename, dtype=dt)

      # Create pandas dataframe, set column names and remove HEADER columns
      df = pd.DataFrame(data, columns=dt.names)
      df = df.drop(['HEADER1', 'HEADER2'], axis=1)

   except IOError:
      df = pd.DataFrame(columns=dt.names)

   '''
   # Read binary file
   data = np.fromfile(filename, dtype=dt)

   # Create pandas dataframe, set column names and remove HEADER columns
   df = pd.DataFrame(data, columns=dt.names)
   df = df.drop(['HEADER1', 'HEADER2'], axis=1)

   return df

def save_df_to_pkl(dataframe, pkl_file):
   '''Save pandas dataframe to pickle file'''
   dataframe.to_pickle(pkl_file)

def load_df_from_pkl(pkl_file):
   '''Load pandas dataframe from pickle file'''
   import pandas as pd

   return pd.read_pickle(pkl_file)



