# -*- coding: utf-8 -*-
#Este script plotea perfiles verticales de la intensidad de las covarianzas.
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
import numpy.ma as ma

basedir='/home/jruiz/share/EXPERIMENTS/experiments_large_ensemble/data/'

expname = '/OsakaPAR_1km_control1000m_smallrandompert_noda/'

plotbasedir=basedir + expname + '/plots/'

undef_in=1.0e20

buffer_zone_size=20  #Data close to the domain borders will be ignored.

tr_rain=5e-4      #QHYD above this tr means rainy grid point.
tr_norain=0.5e-4    #QHYD below this tr means no rain


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

#Covariance strength

base_vars='qhydro','u','v','tk'
cov_vars='qhydro','u','v'

#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,25,00)  #Initial time.
etime = dt.datetime(2013,7,13,5,39,00)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=60)

nx=180
ny=180
nz=np.max(ctl_endrecord) + 1  #Total number of records in binary file.
nlev=12                       #Number of vertical levels for 3D variables.
ntimes=1 + np.around((etime-itime).seconds / delta.seconds)  #Total number of times.
levels=np.array([ 1000,950,900,850,800,700,600,500,400,300,200,150])
levels_str = np.empty([np.size(levels)], dtype="U10")

times=np.zeros((ntimes))

for ii in range(np.size(levels)):
    levels_str[ii]=np.char.mod('%d',levels[ii])

#Define regions


ens_mean=dict()

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]

#for ib in 'u'  :
for ib in base_vars :
#  for ic in 'u'  :
  for ic in cov_vars :

    
    ctime=itime

    it=0

    corrindex_profile=np.zeros([nlev,ntimes,2])
    covindex_profile=np.zeros([nlev,ntimes,2])
    corrdistindex_profile=np.zeros([nlev,ntimes,2])
    ndata=np.zeros([nlev,ntimes,2])

    print ( 'Reading the data ')

    while ( ctime <= etime ):
    
      times[it]=it

      print( ctime )

      var_comb= ib + '_' + ic 

      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' +  'corrindex_' + var_comb + '.grd'

      corrindex=bio.read_data_direct(my_file,nx,ny,nlev,'>f4',undef_in=undef_in,undef_out=undef_out)

      corrindex=ma.masked_where(np.isnan(corrindex),corrindex)
      corrindex.fill_value=0.0

      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' +  'covindex_' + var_comb + '.grd'

      covindex=bio.read_data_direct(my_file,nx,ny,nlev,'>f4',undef_in=undef_in,undef_out=undef_out)
 
      covindex=ma.masked_where(np.isnan(covindex),covindex)
      covindex.fill_value=0.0

      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' +  'corrdistindex_' + var_comb + '.grd'

      corrdistindex=bio.read_data_direct(my_file,nx,ny,nlev,'>f4',undef_in=undef_in,undef_out=undef_out)

      corrdistindex=ma.masked_where(np.isnan(corrdistindex),corrdistindex)
      corrdistindex.fill_value=0.0

      #Read the ensemble mean to get the information from the storm location.
      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/moment0001.grd'

      ens_mean=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

      #Compute total integrated liquid (we will use this to identify areas associated with clouds and convection)
      int_liquid = np.nanmean(ens_mean['QHYD'][buffer_zone_size:-buffer_zone_size,buffer_zone_size:-buffer_zone_size,:],2)

      #Add the profile depending if it is a rainy point or a no rainy point
      rain_mask = int_liquid > tr_rain      #Rain mask
      norain_mask = ( int_liquid < tr_norain ) & ( int_liquid >= 0.0 )  #No rain mask 
 
      for ilev in range(0,nlev) :     
      
         my_index=corrindex[buffer_zone_size:-buffer_zone_size,buffer_zone_size:-buffer_zone_size,ilev]
    
         corrindex_profile[ilev,it,0]=np.nansum( my_index[rain_mask] )
         corrindex_profile[ilev,it,1]=np.nansum( my_index[norain_mask] )

         my_index=covindex[buffer_zone_size:-buffer_zone_size,buffer_zone_size:-buffer_zone_size,ilev]

         covindex_profile[ilev,it,0]=np.nansum( my_index[rain_mask] )
         covindex_profile[ilev,it,1]=np.nansum( my_index[norain_mask] )

         my_index=corrdistindex[buffer_zone_size:-buffer_zone_size,buffer_zone_size:-buffer_zone_size,ilev]

         corrdistindex_profile[ilev,it,0]=np.nansum( my_index[rain_mask] )
         corrdistindex_profile[ilev,it,1]=np.nansum( my_index[norain_mask] )

         ndata[ilev,it,0] = np.nansum( rain_mask )
         ndata[ilev,it,1] = np.nansum( norain_mask )

      ctime = ctime + delta
 
      it = it + 1

    data_mask=( ndata > 0 )
    corrindex_profile[data_mask]=corrindex_profile[data_mask]/ndata[data_mask]
    corrindex_profile[np.logical_not(data_mask)]=np.nan
    data_mask=( ndata > 0 )
    covindex_profile[data_mask]=covindex_profile[data_mask]/ndata[data_mask]
    covindex_profile[np.logical_not(data_mask)]=np.nan
    data_mask=( ndata > 0 )
    corrdistindex_profile[data_mask]=corrdistindex_profile[data_mask]/ndata[data_mask]
    corrdistindex_profile[np.logical_not(data_mask)]=np.nan


    corrindex_profile=ma.masked_where(np.isnan(corrindex_profile),corrindex_profile)
    corrindex_profile.fill_value=0.0
    covindex_profile=ma.masked_where(np.isnan(covindex_profile),covindex_profile)
    covindex_profile.fill_value=0.0
    corrdistindex_profile=ma.masked_where(np.isnan(corrdistindex_profile),corrdistindex_profile)  
    corrdistindex_profile.fill_value=0.0

    my_plotdir= plotbasedir + '/time_independent_plots/'

    if not os.path.exists(my_plotdir) :
       os.makedirs(my_plotdir)


    print("Plotting the CORRINDEX Vertical profile rain")

    varname= var_comb

    var_range_max=np.nanmax(abs(corrindex_profile[1:,:,0]))
    var_range_min=0

    var_range_max_2=np.nanmax(abs(corrdistindex_profile[1:,:,0]))
    var_range_min_2=0

    plt.contourf(times,-np.log(levels[1:]),corrindex_profile[1:,:,0],cmap=plt.cm.get_cmap('seismic'),vmax=var_range_max,vmin=var_range_min)
    plt.colorbar()
    plt.yticks(-np.log(levels[1:]), levels_str[1:], rotation='horizontal')
    #plt.contour(times,-np.log(levels[1:]),corrdistindex_profile[1:,:,0],vmax=var_range_max_2,vmin=var_range_min_2)
   
    plt.title( var_comb ) 
    plt.ylabel('Pressure')
    plt.xlabel('Time')

#    plt.show()

    print( 'Generationg the following figure : ' + 'Figure_corrindex_profile_rain_' + var_comb + '.png' )
    plt.savefig( my_plotdir +  'Figure_corrindex_profile_rain_' + var_comb + '.png' )

    plt.close()

    print("Plotting the CORRINDEX Vertical profile no rain")

    varname= var_comb

    var_range_max=np.nanmax(abs(corrindex_profile[1:,:,1]))
    var_range_min=0

    var_range_max_2=np.nanmax(abs(corrdistindex_profile[:,:,1]))
    var_range_min_2=0

    plt.contourf(times,-np.log(levels[1:]),corrindex_profile[1:,:,1],cmap=plt.cm.get_cmap('seismic'),vmax=var_range_max,vmin=var_range_min)
    plt.colorbar()
    plt.yticks(-np.log(levels[1:]), levels_str[1:], rotation='horizontal')
    #plt.contour(times,-np.log(levels[1:]),corrdistindex_profile[1:,:,0],vmax=var_range_max_2,vmin=var_range_min_2)

    plt.title( var_comb )
    plt.ylabel('Pressure')
    plt.xlabel('Time')

#    plt.show()

    print( 'Generationg the following figure : ' + 'Figure_corrindex_profile_norain_' + var_comb + '.png' )
    plt.savefig( my_plotdir +  'Figure_corrindex_profile_norain_' + var_comb + '.png' )

    plt.close()


print ( "Finish time loop" )


