import os
import argparse
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.table import Table
from matplotlib.font_manager import FontProperties

import util

def get_data(path, ana_time, nslots):
   ''' Get data to plot '''

   # Load configuration
   var_info, obstype_info = util.load_conf() 
   obstypes = [*obstype_info]
   variables = [*var_info]

   # Read data
   obsfile = ('{}/{}/obs_{}_asimiladas.dat').format(path, ana_time, ana_time)
   df = util.read_obsdat(obsfile)

   # Keep data from OBS and remove LON, LAT, LEV, ERROR, OMB, OMA columns
   df = df.drop(['LON', 'LAT', 'LEV', 'ERROR','OMB', 'OMA'], axis=1)

   # Group data by SLOT, observation TYPE and VARIABLE and compute observation count
   df = df.groupby(['SLOT', 'TYPE', 'VAR'], as_index=False).agg(['count']).reset_index()

   # Order data in numpy array 
   data = np.zeros((NSLOTS, len(obstype_info), [*var_info].index('ps')+1))

   # Loop over slots
   for islot in range(NSLOTS):

      # Group statistics for each observation type and variable
      for iobstype, obstype in enumerate(obstypes):
         obstype_id = obstype_info[obstype]['id']

         # Check availability of obstype to plot in the data 
         if obstype_id in df.TYPE.values:
            #print(obstype, 'is present')

            for ivar, var in enumerate(variables):
               var_id = var_info[var]['id']

               # Check availability of variable to plot in the data 
               condition = (df.SLOT == islot+1) & (df.TYPE == obstype_id) & (df.VAR == var_id)
               if any(condition):
                  #print(var, 'is present')

                  row_data = df[condition]

                  if var == 'us':
                     ivar = variables.index('u')
                  if var == 'vs':
                     ivar = variables.index('v')
                  if var == 'ts':
                     ivar = variables.index('t')
                  if var == 'qs':
                     ivar = variables.index('q')
                  if var == 'rhs':
                     ivar = variables.index('rh')

                  data[islot, iobstype, ivar] = row_data['OBS']['count'].values

   return data


def make_plot(data, pathout, ana_time):
   ''' Make plot figure data assimilation statistics '''

   # Load configuration
   var_info, obstype_info = util.load_conf()
   nslot, nobstype, nvar = data.shape
   obstypes = [*obstype_info]
   variables = [*var_info]
   width, height = 1.0/nvar, 1.0/nobstype

   # Get analysis date for title
   ana_date = util.str2date(ana_time)
   ana_date_title = util.date2str(ana_date, '%Y-%m-%d %HZ')
   ana_date = util.date2str(ana_date, '%Y%m%d_%H')

   # Start figure
   title = 'OBSERVACIONES ASIMILADAS \n FECHA DE ANALISIS: {}'.format(ana_date_title)

   fig, ax = plt.subplots(4, 2, figsize = [20, 20], squeeze=True)
   fig.suptitle(title, y=0.92, fontsize=16, fontweight='bold')

   irow = 0
   icol = 0
   for iplot, islot in enumerate(np.arange(nslot + 1)):

      # Determine axis and subplot size
      iax = ax[irow, icol]
      iax.set_axis_off()
      if icol == 0:
         idx = 0.18
      if icol == 1:
         idx = 0.1
      tb = Table(iax, bbox=[idx,0,0.9,0.9])

      if islot < nslot:
         data_slot = data[islot]
         slot = 'SLOT ' + str(islot + 1)
      else:
         data_slot = np.sum(data, axis=0)
         slot = 'TOTAL'

      # Plot values
      for (i, j), val in np.ndenumerate(data_slot):
         color = 'lightgreen' if val > 0 else 'white'

         tb.add_cell(i, j, width, height, text=str(int(val)) if val else '', loc='center', facecolor=color, fontproperties=FontProperties(weight='bold'))

      # Row name
      for i in range(data_slot.shape[0]):
         row_name = obstype_info[obstypes[i]]['name']['short'].upper()
         tb.add_cell(i, -1, width, height, text=row_name, loc='right', edgecolor='none', facecolor='none', fontproperties=FontProperties(weight='bold'))

      # Column name
      for i in range(data_slot.shape[1]):
         col_name = var_info[variables[i]]['name']['short'].upper()
         tb.add_cell(-1, i, width, height/2, text=col_name, loc='center',
              edgecolor='none', facecolor='none', fontproperties=FontProperties(weight='bold'))

      tb.auto_set_font_size(False)
      tb.set_fontsize(12)
      iax.add_table(tb)
      iax.set_title(slot, fontsize='14', fontweight='bold')

      # Uptdate row and column
      icol += 1
      if iplot in np.arange(1,nslot+1,2):
         irow += 1
         icol = 0

   #fig.align_ylabels()       
   plt.savefig('{}/monit_obsdat_obscount.png'.format(pathout), bbox_inches='tight')
   plt.close()

if __name__ == '__main__':

   print('------- Hello from plot_obsdat_obscount -------')
   
   # Get environment variables and arguments
   pathin = os.environ['PATH_OBS']
   pathout = os.environ['PATH_PLOT']

   parser = argparse.ArgumentParser(description='Year Month Day Hour Minute Nslots')
   parser.add_argument('Year',type=int)
   parser.add_argument('Month',type=int)
   parser.add_argument('Day',type=int)
   parser.add_argument('Hour',type=int)
   parser.add_argument('Minute',type=int)
   parser.add_argument('Nslots',type=int)
   args = parser.parse_args()

   NSLOTS = args.Nslots
   Y = str(args.Year).zfill(2)
   M = str(args.Month).zfill(2)
   D = str(args.Day).zfill(2)
   H = str(args.Hour).zfill(2)
   Mi = str(args.Minute).zfill(2)

   ANA_DATE = Y + M + D + '_' +  H + Mi + '00' 

   # Start measuring execution time
   start_time = time.time()

   # Get data to plot
   data = get_data(pathin, ANA_DATE, NSLOTS)

   # Make plot
   make_plot(data, pathout, ANA_DATE)

print('It took', round(time.time()-start_time, 5), 'seconds')
