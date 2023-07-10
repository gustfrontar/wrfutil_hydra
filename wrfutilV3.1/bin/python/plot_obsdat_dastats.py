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
      obsfile = ('{}/{}/obs_{}_asimiladas.dat').format(path, current_time, current_time)
      df = util.read_obsdat(obsfile)

      # Keep data from slot 7 (last) and remove LON, LAT, LEV, ERROR, SLOT columns
      df = df.loc[df['SLOT'] == 7]
      df = df.drop(['LON', 'LAT', 'LEV', 'ERROR', 'SLOT'], axis=1)

      # Add FG and ANA columns with the first guess and analysis mean 
      FG = -1 * (df.OMB - df.OBS)
      ANA = -1 * (df.OMA - df.OBS)
      df = df.assign(FG = FG, ANA = ANA)
 
      # Add date column
      df.insert(0, 'DATE', date)
 
      # Group data by observation TYPE and VARIABLE and compute statistics
      statistics = {'OBS':['count', 'mean', 'std'], 'FG':['mean'], 'ANA':['mean'], 'OMB':['mean', 'std'], 'OMA':['mean', 'std']}
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
   stats = {'plot1': {'stat': 'mean', 'vars': ['OBS', 'FG', 'ANA'], 'colors':['g', 'b', 'r'], 'labels': ['OBS', 'FG', 'ANA']},
            'plot2': {'stat': 'mean', 'vars': ['OMB', 'OMA'], 'colors': ['b', 'r'], 'labels': ['OBS-FG', 'OBS-ANA']},
            'plot3': {'stat': 'std', 'vars': ['OMB', 'OMA', 'OBS'], 'colors': ['b', 'r', 'g'], 'labels': ['std(OBS-FG)', 'std(OBS-ANA)', 'std(OBS)']},
            'plot4': {'stat': 'count', 'vars': ['OBS'], 'colors': ['g'], 'labels': ['NOBS']}}


   # Get analysis date
   dates = pd.to_datetime(data.DATE.unique())
   ana_date_title = util.date2str(dates[-1], '%Y-%m-%d %HZ')
   ana_date = util.date2str(dates[-1], '%Y%m%d_%H')

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

      variables_id = data[(data.TYPE == obstype_id)].VAR.unique() 
      for var_id in variables_id:
         var = util.get_key(var_info, 'id', var_id)

         title = '{} - {} \n ESTADISTICAS - SLOT = 7 \n INTERVALO = 1 HR; ULTIMO DATO: {}'.format(var_info[var]['name']['long'].upper(), obstype_info[obstype]['name']['long'].upper(), ana_date_title)

         # Start figure 
         fig, ax = plt.subplots(4, 1, figsize=(10,10), sharex=True)
         fig.suptitle(title, y=0.97, fontsize=14)

         # Get data to plot
         tmp = data[(data.TYPE == obstype_id) & (data.VAR == var_id)]
         x = pd.to_datetime(tmp.DATE.unique())

         # Suplots
         for iplot, plot in enumerate(stats.keys()):
            plot_opt = stats[plot]
            for idata, data_ in enumerate(plot_opt['vars']):
               stat = plot_opt['stat']
               y = tmp[data_][stat]
               if var == 'q' or var == 'qs':
                  y = y * 1e3

               if obstype == 'SFCSHP' and var == 'ts':
                  continue

               ax[iplot].plot(x, y, color=plot_opt['colors'][idata], label=plot_opt['labels'][idata], marker=marker)

            if plot == 'plot2':
               ax[iplot].axhline(y=0, linestyle=':', linewidth=2, zorder=1, color='k')

            # Legend
            ax[iplot].legend(loc='upper left', handlelength=1.5, prop={'size': 11}, edgecolor='white', ncol=3, bbox_to_anchor=(0.001, 1.22))

            # Axis labels
            ax[iplot].set_ylabel(var_info[var]['units'], fontsize=12)
            if plot == 'plot4':
               ax[iplot].set_ylabel('Numero Obs.')

            # Axis format
            ax[iplot].xaxis.set_major_locator(days)
            ax[iplot].xaxis.set_major_formatter(d_fmt)
            ax[iplot].xaxis.set_minor_locator(hours)
            ax[iplot].xaxis.set_minor_formatter(h_fmt)
            ax[iplot].xaxis.set_tick_params(which='major', length=9, labelsize=12, pad=20)
            ax[iplot].xaxis.set_tick_params(which='minor', length=6, labelsize=12, pad=5)
            ax[iplot].xaxis.remove_overlapping_locs = False
            ax[iplot].xaxis.grid('True', which='both')
            ax[iplot].yaxis.grid('True', which='major')
  
         #fig.autofmt_xdate()
         fig.align_ylabels()       
         plt.savefig('{}/monit_obsdat_dastats_{}_{}.png'.format(pathout, obstype, var), bbox_inches='tight')
         plt.close()

if __name__ == '__main__':

   print('------- Hello from plot_obsdat_dastats -------')
   
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
