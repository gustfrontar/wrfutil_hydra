import os
import argparse
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.table import Table
from matplotlib.font_manager import FontProperties

import util

def get_data(path, ana_time):
   ''' Get data to save obs stats in dataframe '''

   # Load configuration
   var_info, obstype_info = util.load_conf() 
   obstypes = [*obstype_info]
   variables = [*var_info]

   # Read data
   obsfile = '{}/{}/obs_{}_asimiladas.dat'.format(path, ana_time, ana_time)
   df = util.read_obsdat(obsfile)

   # Keep data from OBS and remove LON, LAT, LEV, 'SLOT', ERROR, OMB, OMA columns
   df = df.drop(['LON', 'LAT', 'LEV', 'SLOT', 'ERROR','OMB', 'OMA'], axis=1)

   # Group data by observation TYPE and VARIABLE and compute observation count
   df = df.groupby(['TYPE', 'VAR'], as_index=False).agg(['count', 'mean', 'min', 'max']).reset_index()

   return df

if __name__ == '__main__':

   print('------- Hello from write_obsdat_obsstats -------')
   
   # Get environment variables and arguments
   pathin = os.environ['PATH_OBS']
   pathout = os.environ['PATH_OBS']

   parser = argparse.ArgumentParser(description='Year Month Day Hour Minute')
   parser.add_argument('Year',type=int)
   parser.add_argument('Month',type=int)
   parser.add_argument('Day',type=int)
   parser.add_argument('Hour',type=int)
   parser.add_argument('Minute',type=int)
   args = parser.parse_args()

   Y = str(args.Year).zfill(2)
   M = str(args.Month).zfill(2)
   D = str(args.Day).zfill(2)
   H = str(args.Hour).zfill(2)
   Mi = str(args.Minute).zfill(2)

   ANA_DATE = Y + M + D + '_' +  H + Mi + '00' 

   # Start measuring execution time
   start_time = time.time()

   # Get data as a dataframe 
   data = get_data(pathin, ANA_DATE)

   # Save dataframe to pkl file
   fileout = '{}/{}/obs_{}_obsstats.pkl'.format(pathout, ANA_DATE, ANA_DATE)
   util.save_df_to_pkl(data, fileout)

print('It took', round(time.time()-start_time, 5), 'seconds')
