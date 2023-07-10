import os
import argparse
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
pd.plotting.register_matplotlib_converters()

import util

def get_data(path, dates):
   ''' Get data to plot '''

   for idate, date in enumerate(dates):
      #print(date)
      current_time = util.date2str(date)

      # Read data
      obsfile = '{}/{}/obs_{}_asimiladas.dat'.format(path, current_time, current_time)
      df = util.read_obsdat(obsfile)

      # Keep data from OBS and remove LON, LAT, LEV, SLOT, ERROR, OMB, OMA columns
      df = df.drop(['LON', 'LAT', 'LEV', 'SLOT', 'ERROR','OMB', 'OMA'], axis=1)

      # Add date column
      df.insert(0, 'DATE', date)

      # Group data by observation TYPE, VARIABLE and DATE and compute statistics
      statistics = {'OBS':['count', 'mean', 'min', 'max']}
      df = df.groupby(['TYPE', 'VAR', 'DATE'], as_index=False).agg(statistics)

      # Concatenate dataframes
      if idate == 0:
         cat = df
      else:
         cat = pd.concat([cat, df]).reset_index(drop=True)

   return cat

def make_plot(data, pathout):
   ''' Make plot figure data assimilation statistics '''

   # Load configuration
   var_info, obstype_info = util.load_conf()

   # Get analysis date for title
   dates = pd.to_datetime(data.DATE.unique())
   ana_date = util.date2str(dates[-1], '%Y%m%d_%H')
   ana_date_title = util.date2str(dates[-1], '%Y-%m-%d %HZ')

   # X-axis format for dates
   days = mdates.DayLocator(interval=1)  # every day
   hours= mdates.HourLocator(byhour=[0, 6, 12, 18])  # every 6 hour
   d_fmt = mdates.DateFormatter('%b %d')
   h_fmt = mdates.DateFormatter('%H')

   # Loop over obstypes and variables (only plot SFC stations)
   obstypes_id = data.TYPE.unique()
   for obstype_id in obstypes_id:
      obstype = util.get_key(obstype_info, 'id', obstype_id)
      print(obstype)

      marker=None
      if obstype == 'ADPUPA' or obstype == 'AIRS':
         marker = 'o'

      title = '{} \n VALOR MEDIO Y RANGO OBSERVADO \n INTERVALO = 1 HR; ULTIMO DATO: {}'.format(obstype_info[obstype]['name']['long'].upper(), ana_date_title)

      # Start figure 
      variables_id = data[(data.TYPE == obstype_id)].VAR.unique() 
      fig, ax = plt.subplots(len(variables_id), 1, figsize=[10, len(variables_id)*2+1], sharex=True, squeeze=False)
      ax[0,0].set_title(title, fontsize=14)

      # Subplots
      for iplot, var_id in enumerate(variables_id):
         var = util.get_key(var_info, 'id', var_id)
         #print(var)

         # Get data to plot
         tmp = data[(data.TYPE == obstype_id) & (data.VAR == var_id)]
         x = pd.to_datetime(tmp.DATE.unique())

         obs_mean = tmp.OBS['mean']
         obs_min = tmp.OBS['min']
         obs_max = tmp.OBS['max']
         if var == 'q' or var == 'qs':
            obs_mean *= 1e3
            obs_min *= 1e3
            obs_max  *= 1e3

         ax[iplot, 0].plot(x, obs_mean, color='b', label='mean', marker=marker)
         ax[iplot, 0].fill_between(x, obs_min, obs_max, facecolor='lightblue', alpha=.4)

         # Axis labels
         ax[iplot, 0].set_ylabel('{} \n ({})'.format(var_info[var]['name']['long'], var_info[var]['units']), fontsize=12)

         # Axis format
         ax[iplot, 0].xaxis.set_major_locator(days)
         ax[iplot, 0].xaxis.set_major_formatter(d_fmt)
         ax[iplot, 0].xaxis.set_minor_locator(hours)
         ax[iplot, 0].xaxis.set_minor_formatter(h_fmt)
         ax[iplot, 0].xaxis.set_tick_params(which='major', length=9, labelsize=12, pad=20)
         ax[iplot, 0].xaxis.set_tick_params(which='minor', length=6, labelsize=12, pad=5)
         ax[iplot, 0].xaxis.remove_overlapping_locs = False
         ax[iplot, 0].xaxis.grid('True', which='both')
         ax[iplot, 0].yaxis.grid('True', which='major')

      fig.align_ylabels()       
      plt.savefig('{}/monit_obsdat_obsstats_{}.png'.format(pathout, obstype), bbox_inches='tight')
      plt.close()

if __name__ == '__main__':
   print('------- Hello from plot_obsdat_obsstats -------')  
 
   # Get environment variables and arguments
   pathin = os.environ['PATH_OBS']
   pathout = os.environ['PATH_PLOT']

   parser = argparse.ArgumentParser(description='Year Month Day Hour Minute Ntimes')
   parser.add_argument('Year',type=int)
   parser.add_argument('Month',type=int)
   parser.add_argument('Day',type=int)
   parser.add_argument('Hour',type=int)
   parser.add_argument('Minute',type=int)
   parser.add_argument('Ntimes',type=int)
   args = parser.parse_args()

   NTIMES = args.Ntimes
   Y = str(args.Year).zfill(2)
   M = str(args.Month).zfill(2)
   D = str(args.Day).zfill(2)
   H = str(args.Hour).zfill(2)
   Mi = str(args.Minute).zfill(2)

   ANA_DATE = Y + M + D + '_' +  H + Mi + '00' 

   # Start measuring execution time
   start_time = time.time()

   # Loop over dates to process
   dates = util.get_dates(ANA_DATE, NTIMES, 1)
   df = get_data(pathin, dates)

   # Make plot
   make_plot(df, pathout)

print('It took', round(time.time()-start_time, 5), 'seconds')
