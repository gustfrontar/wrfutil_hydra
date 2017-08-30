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

undef_out=np.nan

filetype='analgp'   #analgp , analgz , guesgp , guesgz

nmoments=4       #Cantidad total de momentos que vamos a plotear.

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
itime = dt.datetime(2013,7,13,5,11,00)  #Initial time.
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

moments=dict()

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

moments=dict()

my_std=dict()

moment_regional_mean=dict()
moment_regional_max=dict()
moment_regional_min=dict()

moment_time_mean=dict()
moment_time_std=dict()

int_liquid=np.zeros([nx,ny,nlev])

time_mean_int_liquid=np.zeros([nx,ny,nlev])

my_moment=np.zeros([nx,ny,nlev])

my_moment_sq=np.zeros([nx,ny,nlev])

my_time_counter=dict()

for key in ctl_vars  :
 
  moment_regional_mean[key]=np.zeros([ntimes,nmoments,nregs])
  moment_regional_max[key] =np.zeros([ntimes,nmoments,nregs])
  moment_regional_min[key] =np.zeros([ntimes,nmoments,nregs])

  moment_time_mean[key]=np.zeros([nx,ny,nlev,nmoments])
  moment_time_std[key]=np.zeros([nx,ny,nlev,nmoments])

  my_std[key]=np.zeros([nx,ny,nlev])

  my_time_counter[key]=np.zeros([nx,ny,nlev,nmoments])

it=0

while ( ctime <= etime ):

  print( ctime )

  for im in range( 0 , nmoments ) :
  
    imstr="%04d" % ( im + 1 ) 

    print ( 'Reading the moment ' + imstr)

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/moment' + imstr + '.grd'

    moments=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out) 

    #Compute total integrated liquid (we will use this to identify areas associated with clouds and convection)
    if im == 0   :
      tmp_int_liquid = np.nansum(moments['QHYD'],2)
      for ilev in range(0,nlev)   :          #Create a fake 3D array for the vertically integrated liquid
                                             #This is because the plotting function expects a 3D array as input.
         int_liquid[:,:,ilev]=tmp_int_liquid



    #Plot moments.
      
    for key in moments  :
      my_moment=moments[key]
      if im == 0 :
         print("plotting the mean")
         varname= 'mean_' + key
         if key == 'W'    :
            myrange='centered'
         if key == 'QHYD' :
            myrange='positive'
         else             :
            myrange='maxmin'
         scale_max='None'
         scale_min='None'         
      elif im == 1 :
         print("plotting the standard deviation")
         my_moment=np.sqrt(my_moment)
         my_std[key]=my_moment
         my_std[key][my_std[key]==0]=np.nan
         varname= 'std_' + key
         myrange='positive'
         scale_max='None'
         scale_min='None'
      elif im == 2 : 
         print("plotting the skewness")
         my_moment=my_moment/np.power( my_std[key] , 3 ) 
         varname= 'skewness_' + key
         myrange='fixed'
         scale_max=2
         scale_min=-2
      elif im == 3 :
         print("plotting the kurtosis")
         my_moment=my_moment/np.power( my_std[key] , 4 ) - 3
         varname='kurtosis_' + key
         myrange='fixed'
         scale_max=20
         scale_min=-20

      print('Moment ',(2),' Var ',key,' ',(np.nanmin(my_std[key])),np.nanmax(my_std[key]))
      print('Moment ',(im+1),' Var ',key,' ',(np.nanmin(my_moment)),np.nanmax(my_moment))

      #Plot the moment.
      bvf.plot_var_levels(my_moment,lon,lat,plotlevels,plotbasedir,varname,varcontour=int_liquid,range=myrange,scale_max=scale_max,scale_min=scale_min,clevels=qmeanlevs,date=ctime.strftime("%Y%m%d%H%M%S"))


      #Generate regional averages of the moment.
      moment_regional_mean[key][it,im,:],moment_regional_max[key][it,im,:],moment_regional_min[key][it,im,:]=bvf.get_regional_average(my_moment,xi,xe,yi,ye)

      #Acumulate the moment statistics in time. 
      
      moment_time_mean[key][:,:,:,im]=moment_time_mean[key][:,:,:,im] + my_moment                 #Accumulate the mean.
      moment_time_std[key][:,:,:,im]=moment_time_std[key][:,:,:,im]   + np.power( my_moment , 2 ) #Accumulate the standard deviation.

      time_mean_int_liquid=time_mean_int_liquid + int_liquid 

  ctime = ctime + delta
 
  it = it + 1

print ( "Finish time loop" )

ntimes=it

#Compute the time averaged moments and plot them.
my_plotdir= plotbasedir + '/time_independent_plots/'

#Create output directory
if not os.path.exists(my_plotdir):
    os.makedirs(my_plotdir)


for key in moment_time_mean  :
        
    my_moment=moment_time_mean[key] / ntimes

    time_mean_int_liquid = time_mean_int_liquid / ntimes

    for im in range(0,nmoments) :

      if im == 0 :
         print("plotting the mean")
         varname= 'mean_mean_' + key
         myrange='maxmin'
         scale_max='None'
         scale_min='None'
      elif im == 1 :
         print("plotting the standard deviation")
         my_moment=np.sqrt(my_moment)
         varname= 'mean_std_' + key
         myrange='maxmin'
         scale_max='None'
         scale_min='None'
      elif im == 2 :
         print("plotting the skewness")
         varname= 'mean_skewness_' + key
         myrange='fixed'
         scale_max=2
         scale_min=-2
      elif im == 3 :
         print("plotting the kurtosis")
         varname='mean_kurtosis_' + key
         myrange='fixed'
         scale_max=20
         scale_min=-20

    #Plot time mean of the moments.

      qmeanlevs=[0.0015,0.0030]
      bvf.plot_var_levels(my_moment[:,:,:,im],lon,lat,plotlevels,plotbasedir,varname,range=myrange,scale_max=scale_max,scale_min=scale_min,clevels=qmeanlevs,varcontour=time_mean_int_liquid,debug=False)

    #Plot time std of the moments. 
 
      my_moment=np.sqrt( moment_time_std[key] / ntimes  - np.power( moment_time_mean[key], 2 ) / ntimes )

      bvf.plot_var_levels(my_moment[:,:,:,im],lon,lat,plotlevels,plotbasedir,varname,range=myrange,scale_max=scale_max,scale_min=scale_min,clevels=qmeanlevs,varcontour=time_mean_int_liquid,debug=False)
    #Plot time series of moments for different regions.



    for im in range(0,nmoments) :
       imstr="%04d" % ( im + 1 )

       for ireg in range(0,nregs) :
         iregstr="%04d" % ( ireg + 1 )

         fig=plt.figure(1,figsize=bvf.default_figure_size)

         plt.plot(moment_regional_mean[key][:,im,ireg],'-')

         plt.ylabel('Moment')
         plt.xlabel('Time')

         #plt.show()

         print( 'Generationg the following figure : ' + 'Figure_moment' + imstr + '_reg_' + iregstr + '_var_' + key + '.png' )
         plt.savefig( my_plotdir +  imstr + '_reg_' + iregstr + '_var_' + key + '.png' )
            



