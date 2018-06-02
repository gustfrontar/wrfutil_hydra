# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import binary_io as bio
import bred_vector_functions as bvf
import os

basedir='/home/jruiz/share/EXPERIMENTS/experiments_large_ensemble/data/'

expname = '/OsakaPAR_1km_control1000m_smallrandompert_new/'

plotbasedir=basedir + expname + '/plots/'

undef_in=1.0e20

nbins=50
thresholdmin=0.005

undef_out=np.nan

filetype='analgp'   

qmeanlevs=[0.001,0.0050,0.05]   #Levels for vertically averaged condensate.


#The following will be used to extract a particlar variable from the original data.
#This variables should be especified according to the data that we have in the binary files.

ctl_vars='U','V','W','T','QV','QHYD'  #Complete list of variables in ctl file.
ctl_inirecord=[0,12,24,36,48,60]        #Starting record for each variable. From 0 to N
ctl_endrecord=[11,23,35,47,59,71]       #End record for each variable. From 0 to N.

#Which variables and levels are we going to plot?

plotlevels=np.array([3,7,9])            #Which levels will be plotted (this levels are equivalent to the BV plots)
plotvars='U','V','W','T','QV','QHYD'    #Which variables will be plotted.

#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,25,00)  #Initial time.
etime = dt.datetime(2013,7,13,5,25,00)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=60)

nx=180
ny=180
nz=np.max(ctl_endrecord) + 1  #Total number of records in binary file.
nlev=12                       #Number of vertical levels for 3D variables.
ntimes=1 + np.around((etime-itime).seconds / delta.seconds)  #Total number of times.


#Define regions

my_hist=dict()

ens_mean=dict()

ctime=itime 

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]

int_liquid=np.zeros([nx,ny,nlev])

it=0

while ( ctime <= etime ):

 print( ctime )

 print ( 'Reading the histogram ')

 hist_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/histogram.grd'
 max_file =basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/ensmax.grd'
 min_file =basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/ensmin.grd'

 hist=read_histogram(hist_file,max_file,min_file,nx,ny,nbins,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='>f4',undef_in=undef_in,undef_out=undef_out)

 hist_properties=analyze_histogram_fun( my_hist , thresholdmin )

 ens_mean=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

 #Compute total integrated liquid (we will use this to identify areas associated with clouds and convection)
 tmp_int_liquid = np.nansum(ens_mean['QHYD'],2)
 for ilev in range(0,nlev)   :          #Create a fake 3D array for the vertically integrated liquid
                                        #This is because the plotting function expects a 3D array as input.
  int_liquid[:,:,ilev]=tmp_int_liquid

 for key in hist :

  #Plot moments.
  my_hist_properties=hist_properties[key]

  tmp=output[key]['min1loc']

  

  

  





 ctime = ctime + delta
 
 it = it + 1

print ( "Finish time loop" )

            



