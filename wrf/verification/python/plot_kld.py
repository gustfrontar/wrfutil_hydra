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

expname = '/OsakaPAR_1km_control1000m_smallrandompert_noda/'

plotbasedir=basedir + expname + '/plots/'

undef_in=1.0e20

undef_out=np.nan

filetype='analgp'   #analgp , analgz , guesgp , guesgz

qmeanlevs=[0.001,0.0050,0.05]   #Levels for vertically averaged condensate.


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
itime = dt.datetime(2013,7,13,5,39,00)  #Initial time.
etime = dt.datetime(2013,7,13,5,39,00)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=60)

nx=180
ny=180
nz=np.max(ctl_endrecord) + 1  #Total number of records in binary file.
nlev=12                       #Number of vertical levels for 3D variables.
ntimes=1 + np.around((etime-itime).seconds / delta.seconds)  #Total number of times.


#Define regions

lati=np.array([34.75,34.6])
late=np.array([35.25,34.9])
loni=np.array([135.5,135.4])
lone=np.array([136.25,135.7])

reg_name='REG_1','REG_2','TOTAL'

kld=dict()

ctime=itime 

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]

#Add the global domain as a region.
lati=np.append(lati,lat[0,0])
late=np.append(late,lat[nx-1,ny-1])
loni=np.append(loni,lon[0,0])
lone=np.append(lone,lon[nx-1,ny-1])

#Compute xi,xe,yi,ye corresponding to each region.
xi , yi = bvf.lat_lon_to_i_j(lon,lat,loni,lati)
xe , ye = bvf.lat_lon_to_i_j(lon,lat,lone,late)

nregs=xi.shape[0]

kld=dict()

ens_mean=dict()

kld_regional_mean=dict()
kld_regional_max=dict()
kld_regional_min=dict()

kld_time_mean=dict()
kld_time_std=dict()

int_liquid=np.zeros([nx,ny,nlev])

time_mean_int_liquid=np.zeros([nx,ny,nlev])

my_kld=np.zeros([nx,ny,nlev])

my_kld_sq=np.zeros([nx,ny,nlev])

for key in ctl_vars  :
 
  kld_regional_mean[key]=np.zeros([ntimes,nregs])
  kld_regional_max[key] =np.zeros([ntimes,nregs])
  kld_regional_min[key] =np.zeros([ntimes,nregs])

  kld_time_mean[key]=np.zeros([nx,ny,nlev])
  kld_time_std[key]=np.zeros([nx,ny,nlev])

it=0

while ( ctime <= etime ):

 print( ctime )

 print ( 'Reading the kld ')

 my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/kldistance.grd'

 kld=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

 #Read the ensemble mean to get the information from the storm location.
 my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/moment0001.grd'

 ens_mean=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

 #Compute total integrated liquid (we will use this to identify areas associated with clouds and convection)
 tmp_int_liquid = np.nansum(ens_mean['QHYD'],2)
 for ilev in range(0,nlev)   :          #Create a fake 3D array for the vertically integrated liquid
                                        #This is because the plotting function expects a 3D array as input.
  int_liquid[:,:,ilev]=tmp_int_liquid


 for key in kld :

  #Plot moments.
  my_kld=kld[key]
  print("plotting the kld")
  varname= 'kld_' + key
  myrange='fixed'
  scale_max=0.15
  scale_min=0

  print('Kld for Var ',key,' ',(np.nanmin(my_kld)),np.nanmax(my_kld))

  #Plot the moment.
  bvf.plot_var_levels(my_kld,lon,lat,plotlevels,plotbasedir,varname,varcontour=int_liquid,range=myrange,scale_max=scale_max,scale_min=scale_min,mycolorbar='hot_r',clevels=qmeanlevs,date=ctime.strftime("%Y%m%d%H%M%S"))


  #Generate regional averages of the moment.
  kld_regional_mean[key][it,:],kld_regional_max[key][it,:],kld_regional_min[key][it,:]=bvf.get_regional_average(my_kld,xi,xe,yi,ye)

  #Acumulate the moment statistics in time. 
      
  kld_time_mean[key]=kld_time_mean[key] + my_kld                 #Accumulate the mean.
  kld_time_std[key]=kld_time_std[key]   + np.power( my_kld , 2 ) #Accumulate the standard deviation.

 time_mean_int_liquid=time_mean_int_liquid + int_liquid  #To accumulate integrated liquid.

 ctime = ctime + delta
 
 it = it + 1

print ( "Finish time loop" )

ntimes=it

#Compute the time averaged moments and plot them.
my_plotdir= plotbasedir + '/time_independent_plots/'

#Create output directory
if not os.path.exists(my_plotdir):
    os.makedirs(my_plotdir)


for key in kld_time_mean  :
        
    my_kld=kld_time_mean[key] / ntimes

    time_mean_int_liquid = time_mean_int_liquid / ntimes


    print("plotting the kld")
    varname= 'mean_kld_' + key
    myrange='fixed'
    scale_max='1.0'
    scale_min='0.0'

    #Plot time mean of the moments.

    qmeanlevs=qmeanlevs*0.1 #I reduce the countours because this is a time mean.
    bvf.plot_var_levels(my_kld,lon,lat,plotlevels,plotbasedir,varname,range=myrange,scale_max=scale_max,scale_min=scale_min,mycolorbar='hot_r',clevels=qmeanlevs,varcontour=time_mean_int_liquid,debug=False)

    #Plot time std of the moments. 
    my_moment=np.sqrt( kld_time_std[key] / ntimes  - np.power( kld_time_mean[key], 2 ) / ntimes )

    bvf.plot_var_levels(my_kld,lon,lat,plotlevels,plotbasedir,varname,range=myrange,scale_max=scale_max,scale_min=scale_min,clevels=qmeanlevs,varcontour=time_mean_int_liquid,debug=False)
    #Plot time series of moments for different regions.

    for ireg in range(0,nregs) :
       iregstr="%04d" % ( ireg + 1 )

       fig=plt.figure(1,figsize=bvf.default_figure_size)

       plt.plot(kld_regional_mean[key][:,ireg],'-')

       plt.ylabel('Kld')
       plt.xlabel('Time')

       #plt.show()

       print( 'Generationg the following figure : ' + 'Figure_kld_reg_' + iregstr + '_var_' + key + '.png' )
       plt.savefig( my_plotdir +  'Figure_kld_reg_' + iregstr + '_var_' + key + '.png' )
            



