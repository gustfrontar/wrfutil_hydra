# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""

import numpy as np
import matplotlib as plt
import datetime as dt
import binary_io as bio
import bred_vector_functions as bvf
import os 

basedir='/data9/jruiz/EXPERIMENTS/'

expname = '/OsakaPAR_1km_control1000m_smallrandompert_new/'

plotbasedir=basedir + expname + '/plots/'

#The following will be used to extract a particlar variable from the original data.
#This variables should be especified according to the data that we have in the binary files.

ctl_vars='U','V','W','T','QV','QHYD'  #Complete list of variables in ctl file.
ctl_inirecord=[0,12,24,36,48,60]        #Starting record for each variable. From 0 to N
ctl_endrecord=[11,23,35,47,59,71]       #End record for each variable. From 0 to N.

undef_in=1.0e20
undef_out=np.nan

qmeanlevs=[0.001,0.0050,0.05]   #Levels for vertically averaged condensate.

#Which variables and levels are we going to plot?

plotlevels=np.array([3,7,9])   #Which levels will be plotted (this levels are equivalent to the BV plots)
plotvars='U','V','W','T','QV','QHYD'    #Which variables will be plotted.

#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,11,00)  #Initial time.
etime = dt.datetime(2013,7,13,5,39,00)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=60)

nx=180
ny=180
nz=nz=np.max(ctl_endrecord) + 1
nlev=12

data_gues=dict()
data_anal=dict()

increment=dict()

int_liquid=np.zeros([nx,ny,nlev])

ctime=itime 

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]

while ( ctime <= etime ):
 
 print ( 'The date is :', ctime )

 print ( 'Reading the first guess')
 
 my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/mean.grd'

 data_gues=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

 print ( 'Reading the analysis')

 my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/mean.grd'

 data_anal=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

 increment=bvf.data_diff( data_anal , data_gues )

 tmp_int_liquid = np.nansum(data_anal['QHYD'],2)

 for ilev in range(0,nlev)   :          #Create a fake 3D array for the vertically integrated liquid
                                        #This is because the plotting function expects a 3D array as input.
  int_liquid[:,:,ilev]=tmp_int_liquid

 for key in increment        :

   my_increment=increment[key]

   varname='increment_' + key
   myrange='centered'

   #Plot the increments.
   bvf.plot_var_levels(my_increment,lon,lat,plotlevels,plotbasedir,varname,varcontour=int_liquid,range=myrange,clevels=qmeanlevs,date=ctime.strftime("%Y%m%d%H%M%S"))

 ctime = ctime + delta

print ( "Finish time loop" )





