import os
import argparse
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import warnings
warnings.filterwarnings("ignore")

import util 

def get_data(path, dates):
   ''' Get data to plot '''

   for idate, date in enumerate(dates):
      #print(date)
      current_time = util.date2str(date)

      # Read data
      obsfile = ('{}/{}/obs_{}_asimiladas.dat').format(path, current_time, current_time)
      df = util.read_obsdat(obsfile)

      # Keep data from TYPE, VAR, LON, LAT, OBS and remove SLOT, LEV, ERROR, OMB, OMA columns
      df = df.drop(['SLOT', 'LEV', 'ERROR','OMB', 'OMA'], axis=1)

      # Concatenate dataframes
      if idate == 0:
         cat = df
      else:
         cat = pd.concat([cat, df]).reset_index(drop=True)

   # Group data by observation TYPE, VARIABLE, LON and LAT and compute observation statistics
   return cat.groupby(['TYPE', 'VAR', 'LON', 'LAT'], as_index=False).agg(['count', 'mean', 'min', 'max']).reset_index()


def make_plot(data, pathout, dates):
   ''' Make plot figure observation mean value '''

   # Load configuration
   var_info, obstype_info = util.load_conf()

   # Basemap limits
   lat_min = -60
   lat_max = -10
   lon_min = 265
   lon_max = 325

   # Get analysis date, analysis cycle and period
   period = '{} - {}'.format(util.date2str(dates[0], '%Y-%m-%d %HZ'), util.date2str(dates[-1], '%Y-%m-%d %HZ'))
   ana_cycle = str(dates[0].hour).zfill(2)
   ana_date = util.date2str(dates[-1], '%Y%m%d_%H')

   # Loop over obstypes and variables (only plot SFC stations)
   obstypes_id = data.TYPE.unique()
   for obstype_id in obstypes_id:
      obstype = util.get_key(obstype_info, 'id', obstype_id)
      if obstype == 'ADPSFC' or obstype == 'ADPAUT':
         print(obstype)

         variables_id = data[(data.TYPE == obstype_id)].VAR.unique() 
         for var_id in variables_id:
            var = util.get_key(var_info, 'id', var_id)
            #print(var)

            # Get data to plot
            tmp = data[(data.TYPE == obstype_id) & (data.VAR == var_id)]
            lon = tmp.LON
            lat = tmp.LAT
            obs = tmp.OBS['mean']

            # Compute aereal mean, min and max 
            obs_mean = str(round(obs.mean(), 1))
            obs_min = str(round(obs.min(), 1))
            obs_max = str(round(obs.max(), 1))
 
            # Start figure with basemap
            fig = plt.figure(figsize=(15,12))
            m = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max, projection='cyl',resolution = 'l')
            m.drawcoastlines(linewidth=0.5)
            m.drawcountries(linewidth=0.5)
            m.drawstates(linewidth=0.25)
            m.drawparallels(np.arange(lat_min,lat_max,5.),labels=[1,0,0,0],fontsize=18)
            m.drawmeridians(np.arange(lon_min,lon_max+1,10.),labels=[0,0,0,1],fontsize=18)

            lons, lats = m(lon, lat)
            l1 = m.scatter(lons, lats, c=obs, s=25, cmap='viridis', marker = 'o')

            cbar = plt.colorbar(l1, shrink=0.8, extend='both', pad=0.02)
            cbar.set_label(var_info[var]['units'], fontsize=14)
            cbar.ax.tick_params(labelsize=12)

            plt.title('{} - {} \n VALOR MEDIO OBSERVADO - CICLO DE ANALISIS = {}Z \n PERIODO = {} \n MEDIA = {} {}; MIN = {} {}; MAX = {} {}'.format(var_info[var]['name']['long'].upper(), obstype_info[obstype]['name']['long'].upper(), ana_cycle, period, obs_mean, var_info[var]['units'], obs_min, var_info[var]['units'], obs_max, var_info[var]['units']), fontsize=14)

            plt.savefig('{}/monit_obsdat_obsmean_{}_{}.png'.format(pathout, obstype, var), bbox_inches='tight')
            plt.close()

if __name__ == '__main__':

   print('------- Hello from plot_obsdat_obsmean -------')   

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
   dates = util.get_dates(ANA_DATE, NTIMES, 24)
   df = get_data(pathin, dates)

   # Make plot
   make_plot(df, pathout, dates)

print('It took', round(time.time()-start_time, 5), 'seconds')
