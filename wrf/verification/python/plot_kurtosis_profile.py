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
levels=np.array([ 1000,950,900,850,800,700,600,500,400,300,200,150])
levels_str = np.empty([np.size(levels)], dtype="U10")

times=np.zeros((ntimes))

for ii in range(np.size(levels)):
    levels_str[ii]=np.char.mod('%d',levels[ii])

#Define regions

lati=np.array([34.75,34.6])
late=np.array([35.25,34.9])
loni=np.array([135.5,135.4])
lone=np.array([136.25,135.7])

reg_name='REG_1','REG_2','TOTAL'

kurt=dict()

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

kurt=dict()
sigma=dict()
vprofile=dict()
ndata=dict()

ens_mean=dict()

kurt_regional_mean=dict()
kurt_regional_max=dict()
kurt_regional_min=dict()

kurt_time_mean=dict()
kurt_time_std=dict()

int_liquid=np.zeros([nx,ny])

my_kurt=np.zeros([nx,ny,nlev])

my_kurt_sq=np.zeros([nx,ny,nlev])

for key in ctl_vars  :
 
  kurt_regional_mean[key]=np.zeros([ntimes,nregs])
  kurt_regional_max[key] =np.zeros([ntimes,nregs])
  kurt_regional_min[key] =np.zeros([ntimes,nregs])

  kurt_time_mean[key]=np.zeros([nx,ny,nlev])
  kurt_time_std[key]=np.zeros([nx,ny,nlev])

it=0

while ( ctime <= etime ):
    
    
    times[it]=it

    print( ctime )

    print ( 'Reading the kurt ')

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/moment0004.grd'

    kurt=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/moment0002.grd'

    sigma=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)


 #Read the ensemble mean to get the information from the storm location.
    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ filetype + '/' + '/moment0001.grd'

    ens_mean=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

 #Compute total integrated liquid (we will use this to identify areas associated with clouds and convection)
    int_liquid = np.nanmean(ens_mean['QHYD'][buffer_zone_size:-buffer_zone_size,buffer_zone_size:-buffer_zone_size,:],2)

    if( ctime == itime ): #Initialize vprofile
        for key in kurt :
           vprofile[key]=np.zeros([nlev,ntimes,2])
           ndata[key]=np.zeros([nlev,ntimes,2])
           wprofile=np.zeros([nlev,ntimes,2])


 #Add the profile depending if it is a rainy point or a no rainy point
    rain_mask = int_liquid > tr_rain      #Rain mask
    norain_mask = ( int_liquid < tr_norain ) & ( int_liquid >= 0.0 )  #No rain mask 
 
 #For each horizontal grid point accumulate the vertical profile discriminating between rain and no-rain
    for ilev in range(0,nlev) :
        my_w=ens_mean['W'][buffer_zone_size:-buffer_zone_size,buffer_zone_size:-buffer_zone_size,ilev]
        wprofile[ilev,it,0]=np.nansum(np.abs(my_w[rain_mask]))
        wprofile[ilev,it,1]=np.nansum(np.abs(my_w[norain_mask]))       
        
    for key in kurt :

      for ilev in range(0,nlev) :     

      
         my_kurt=kurt[key][buffer_zone_size:-buffer_zone_size,buffer_zone_size:-buffer_zone_size,ilev]
         my_sigma=sigma[key][buffer_zone_size:-buffer_zone_size,buffer_zone_size:-buffer_zone_size,ilev]   
         my_kurt=my_kurt/np.power(my_sigma,2) - 3    
    
         vprofile[key][ilev,it,0]=np.nansum(my_kurt[rain_mask])
         vprofile[key][ilev,it,1]=np.nansum(my_kurt[norain_mask])

         ndata[key][ilev,it,0] = np.nansum( rain_mask )
         ndata[key][ilev,it,1] = np.nansum( norain_mask )

    ctime = ctime + delta
 
    it = it + 1

print ( "Finish time loop" )

ntimes=it

tmp=wprofile[:,:,0]
tmp_ndata=ndata['W'][:,:,0]
data_mask = tmp_ndata > 0 
tmp[data_mask]=tmp[data_mask]/tmp_ndata[data_mask]
tmp[np.logical_not(data_mask)]=np.nan
wprofile[:,:,0]=tmp

tmp=wprofile[:,:,1]
tmp_ndata=ndata['W'][:,:,1]
data_mask = tmp_ndata > 0 
tmp[data_mask]=tmp[data_mask]/tmp_ndata[data_mask]
tmp[np.logical_not(data_mask)]=np.nan
wprofile[:,:,1]=tmp

for key in vprofile :

  #Get the average
    tmp=vprofile[key][:,:,0]
    tmp_ndata=ndata[key][:,:,0]
    data_mask= tmp_ndata > 0 
    tmp[data_mask]=tmp[data_mask]/tmp_ndata[data_mask]
    tmp[np.logical_not(data_mask)]=np.nan
    vprofile[key][:,:,0]=tmp

    tmp=vprofile[key][:,:,1]
    tmp_ndata=ndata[key][:,:,1]
    data_mask= tmp_ndata > 0        
    tmp[data_mask]=tmp[data_mask]/tmp_ndata[data_mask]
    tmp[np.logical_not(data_mask)]=np.nan
    vprofile[key][:,:,1]=tmp

#Compute the time averaged moments and plot them.
my_plotdir= plotbasedir + '/time_independent_plots/'

print( my_plotdir )

#Create output directory
if not os.path.exists(my_plotdir):
    os.makedirs(my_plotdir)


for key in vprofile  :
        
    my_profile=vprofile[key] 

    print("Plotting the kurt Vertical profile rain")
    varname= key

    var_range_max=np.nanmax(abs(my_profile[:,:,0]))
    var_range_min=-var_range_max


    plt.contourf(times,-np.log(levels),my_profile[:,:,0],cmap=plt.cm.get_cmap('YlGnBu'),vmax=var_range_max,vmin=0)
    plt.colorbar()
    plt.yticks(-np.log(levels), levels_str, rotation='horizontal')
    plt.contour(times,-np.log(levels),wprofile[:,:,0],np.arange(1,5,1),vmax=5,vmin=1,)



    plt.ylabel('Pressure')
    plt.xlabel('Time')

    plt.show()

    print( 'Generationg the following figure : ' + 'Figure_kurt_profile_rain_var_' + key + '.png' )
    plt.savefig( my_plotdir +  'Figure_kurt_profile_var_' + key + '.png' )

    my_profile=vprofile[key] 

    print("Plotting the kurt Vertical profile no rain")
    varname= key

    var_range_max=np.nanmax(abs(my_profile[:,:,0]))
    var_range_min=-var_range_max


    plt.contourf(times,-np.log(levels),my_profile[:,:,1],cmap=plt.cm.get_cmap('YlGnBu'),vmax=var_range_max,vmin=0)
    plt.colorbar()
    plt.contour(times,-np.log(levels),wprofile[:,:,1],np.arange(0.0,0.5,0.05),vmax=5,vmin=1,)
    plt.yticks(-np.log(levels), levels_str, rotation='horizontal')

    plt.ylabel('Pressure')
    plt.xlabel('Time')

    plt.show()

    print( 'Generationg the following figure : ' + 'Figure_kurt_profile_norain_var_' + key + '.png' )
    plt.savefig( my_plotdir +  'Figure_kurt_profile_var_' + key + '.png' )