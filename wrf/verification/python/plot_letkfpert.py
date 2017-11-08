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

basedir='/data9/jruiz/EXPERIMENTS/'

expname = '/OsakaPAR_1km_control1000m_smallrandompert_new/'

plotbasedir=basedir + expname + '/plots/'

undef_in=1.0e20
undef_out=np.nan

qmeanlevs=[0.001,0.0050,0.05]   #Levels for vertically averaged condensate.

inipert=1   #Which is the first perturbation that we will plot.
endpert=1   #Which is the last perturbation that we will plot.
npert=endpert-inipert+1 #Total number of perturbations that will be plotted.

#The following will be used to extract a particlar variable from the original data.
#This variables should be especified according to the data that we have in the binary files.

ctl_vars='U','V','W','T','QV','QHYD'  #Complete list of variables in ctl file.
ctl_inirecord=[0,12,24,36,48,60]        #Starting record for each variable. From 0 to N
ctl_endrecord=[11,23,35,47,59,71]       #End record for each variable. From 0 to N.

#Which variables and levels are we going to plot?

plotlevels=np.array([3,7,9])   #Which levels will be plotted (this levels are equivalent to the BV plots)
plotvars='U','V','W','T','QV','QHYD'    #Which variables will be plotted.


#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,10,30)  #Initial time.
etime = dt.datetime(2013,7,13,5,39,30)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=30)

nx=180
ny=180
nz=np.max(ctl_endrecord) + 1  #Total number of records in binary file.
nlev=12                       #Number of vertical levels for 3D variables.
ntimes=round( (itime-etime).seconds/delta.seconds ) + 1 #Total number of times.

data_pert_anal=dict()
data_mean_anal=dict()

data_pert_gues=dict()
data_mean_gues=dict()

ctime=itime 

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]

int_liquid=np.zeros([nx,ny,nlev])

it=0

for ipert in range (inipert , endpert + 1):

   pertstr="%04d" % ipert

   print( ' Plotting bred vector number ' + pertstr )

   while ( ctime <= etime ):
 
      ptime=ctime - delta #Data correspinding to the previous step (to compute bv growth)
 
      print ( 'The date is :', ctime )

      print ( 'Reading the perturbed gues')

      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/' + pertstr + '.grd'

      data_pert_gues=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

      print ( 'Reading the gues mean')

      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/mean.grd'

      data_mean_gues=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

      #Compute total integrated liquid (we will use this to identify areas associated with clouds and convection)
      tmp_int_liquid = np.nansum(data_mean_gues['QHYD'],2)

      for ilev in range(0,nlev)   :          #Create a fake 3D array for the vertically integrated liquid
                                              #This is because the plotting function expects a 3D array as input.
       int_liquid[:,:,ilev]=tmp_int_liquid

      #Note: In pithon when multiple variables are output from a function in a tuple, then all the variables has to be decodified. 
      #If not the reconstruction of the variables will fail.

      pert=bvf.data_diff( data_pert_gues , data_mean_gues )


      #Plot LETKF perturbation norm.
      mydir=plotbasedir + '/' + pertstr + '/'

      #Plot the perturbation structure.

      for key in plotvars :

       varname='pert_' + key
       my_range='centered'
       bvf.plot_var_levels(pert[key],lon,lat,plotlevels,mydir,varname,varcontour=int_liquid,clevels=qmeanlevs,range=my_range,date=ctime.strftime("%Y%m%d%H%M%S"))

      ctime = ctime + delta
 
      it = it + 1

print ( "Finish time loop" )





