# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""

# LECTURA Y GRAFICADO RADAR (Formato binario GVAR-SMN)

import numpy as np
import matplotlib as plt
import datetime as dt
import binary_io as bio
import bred_vector_functions as bvf
import os 

basedir='/home/jruiz/share/exp/'

expname = '/breeding_osaka_pawr_1km_bip5_local_1000m_UVT/'

plotbasedir=basedir + expname + '/plots/'

nbredvector=1    #Total number of bred vectors.
inibv=1          #Initial bred vector to plot.
endbv=1          #Final bred vector to plot.
nbipiter=1       #Total number of iterations if breeding in place is activated.
iniit=1          #Initial iter to plot.
endit=1          #End iter to plot.

plotlevels=np.array([6,13,17])   #Which levels will be plotted.
plotvars=['U','V','W','T','QV','QHYD']    #Which variables will be plotted.
#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,35,00)  #Initial time.
etime = dt.datetime(2013,7,13,5,39,30)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=30)

nx=180
ny=180
nz=20

data_pp_o=dict()

ctime=itime + delta

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat_d01z001.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon_d01z001.grd',nx,ny,1,'>f4')[:,:,0]


for ibv in range (inibv , endbv + 1):

   bvstr="%04d" % ibv

   print( ' Plotting bred vector number ' + bvstr )

   while ( ctime <= etime ):

    for iter in range ( iniit , endit + 1 ):

      iterstr="%04d" % iter
 
      ptime=ctime - delta #Data correspinding to the previous step (to compute bv growth)
 
      print ( 'The date is :', ctime )

      print ( 'Reading the positive perturbation original')

      mydir=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pp_o' + iterstr + '/'
      data_pp_o=bio.read_data_scale(mydir,expname,ctime,nx,ny,nz)

      #Plot BV
      mydir=plotbasedir + '/' + '/' + bvstr + '/' + iterstr + '/'  
      bvf.plot_state(data_pp_o,lon,lat,plotvars,plotlevels,mydir,range='maxmin',mycolorbar='coolwarm',date=ctime.strftime("%Y%m%d%H%M%S"))

      ctime = ctime + delta

print ( "Finish time loop" )





