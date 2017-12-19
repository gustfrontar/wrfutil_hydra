# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""

# LECTURA Y GRAFICADO RADAR (Formato binario GVAR-SMN)

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import binary_io as bio
import bred_vector_functions as bvf
import os 

basedir='/home/jruiz/share/EXPERIMENTS/'

expname='/breeding_osaka_pawr_1km_bip5_local_1000m_UVT/'

basedir_kld='/home/jruiz/share/EXPERIMENTS/experiments_large_ensemble/data/'
expname_kld='OsakaPAR_1km_control1000m_smallrandompert_new'

plotbasedir= basedir_kld + expname_kld + '/plots/'

nbredvector=1    #Total number of bred vectors.
inibv=1          #Initial bred vector to plot.
endbv=1          #Final bred vector to plot.
nbipiter=1       #Total number of iterations if breeding in place is activated.
iniit=1          #Initial iter to plot.
endit=1          #End iter to plot.

undef_in=1.0e20
undef_out=np.nan

qmeanlevs=[0.001,0.0050,0.05]   #Levels for vertically averaged condensate.

norm_type='UVT'
smooth_type='Gaussian'
smooth_sigma=1.5
smooth_cutoff=10

#The following will be used to extract a particlar variable from the original data.
#This variables should be especified according to the data that we have in the binary files.

#Variables
ctl_vars='U','V','W','T','QV','QHYD'    #Complete list of variables in ctl file.
ctl_inirecord=[0,12,24,36,48,60]        #Starting record for each variable. From 0 to N
ctl_endrecord=[11,23,35,47,59,71]       #End record for each variable. From 0 to N.

#Dimensions
nx=180
ny=180
nz=np.max(ctl_endrecord) + 1  #Total number of records in binary file.
nlev=12                       #Number of vertical levels for 3D variables.

#Region definition.
lati=np.array([34.75])
late=np.array([35.25])
loni=np.array([135.0])
lone=np.array([136.25])

reg_name='CONVECTION','TOTAL'


#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,15,00)  #Initial time.
etime = dt.datetime(2013,7,13,5,39,00)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=60)

bvfreq=dt.timedelta(seconds=30)

ntimes=round( (etime-itime).seconds/delta.seconds ) + 1

nx=180
ny=180
nz=20
nz_kld=np.max(ctl_endrecord)+1

data_pp_o=dict()
data_pn_o=dict()

data_pp_r=dict()
data_pn_r=dict()

ctime=itime 

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat_d01z001.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon_d01z001.grd',nx,ny,1,'>f4')[:,:,0]

#Add the global domain as a region.
lati=np.append(lati,lat[0,0])
late=np.append(late,lat[nx-1,ny-1])
loni=np.append(loni,lon[0,0])
lone=np.append(lone,lon[nx-1,ny-1])

#Compute xi,xe,yi,ye corresponding to each region.
xi , yi = bvf.lat_lon_to_i_j(lon,lat,loni,lati)
xe , ye = bvf.lat_lon_to_i_j(lon,lat,lone,late)

nregs=lati.size

kld_series=np.zeros((ntimes,nregs))
growthrate_series=np.zeros((ntimes,nregs))
norm_series=np.zeros((ntimes,nregs))

kld=dict()

time_mean_kld=np.zeros((nx,ny))
time_mean_kld=np.zeros((nx,ny))

it=0

for ibv in range (inibv , endbv + 1):

   bvstr="%04d" % ibv

   print( ' Plotting bred vector number ' + bvstr )

   while ( ctime <= etime ):

    for iter in range ( iniit , endit + 1 ):

      iterstr="%04d" % iter
 
      ptime=ctime - bvfreq #Data correspinding to the previous step (to compute bv growth)
 
      print ( 'The date is :', ctime )

      print ( 'Reading the positive perturbation original')

      mydir=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pp_o' + iterstr + '/'
      data_pp_o=bio.read_data_scale(mydir,expname,ctime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)

      print ( 'Reading the negative perturbation original')

      mydir=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pn_o' + iterstr + '/'
      data_pn_o=bio.read_data_scale(mydir,expname,ctime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)

      print ( 'Reading the positive perturbation rescaled')

      mydir=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pp_r' + iterstr + '/'
      data_pp_r=bio.read_data_scale(mydir,expname,ptime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)

      print ( 'Reading the negative perturbation rescaled')

      mydir=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pn_r' + iterstr + '/'
      data_pn_r=bio.read_data_scale(mydir,expname,ptime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)

      print ( 'Reading the KLD distance')        

      kld_file=basedir_kld + expname_kld + '/' +  ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/kldistance.grd'

      kld=bio.read_data_scale_2(kld_file,nx,ny,nz_kld,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

      #Note: In python when multiple variables are output from a function in a tuple, then all the variables has to be decodified. 
      #If not the reconstruction of the variables will fail.

      norm_o =bvf.norm_bv_2( data_pp_o , data_pn_o , norm_type=norm_type )
      norm_r =bvf.norm_bv_2( data_pp_r , data_pn_r , norm_type=norm_type )

      mean_norm_o=np.nanmean( norm_o , 2 )
 
      norm_series[it,:],nada,nada=bvf.get_regional_average(mean_norm_o,xi,xe,yi,ye)
 
      mean_norm_o = bvf.filter_field( mean_norm_o , smooth=smooth_type , sigma=smooth_sigma , cutoff=smooth_cutoff )

      growth_rate = norm_o[:,:,:] - norm_r[:,:,:]  #3d growing rate field.

      growthrate_series[it,:],nada,nada=bvf.get_regional_average(growth_rate,xi,xe,yi,ye)

      #Define the vertically averaged growth rate
      mean_growth_rate=np.nanmean( growth_rate , 2 )
 
      #Compute the vertically averaged kld for the same variables used to compute the BV norm. 
      mean_kld=( np.nanmean( kld['U'] , 2 ) + np.nanmean( kld['V'] , 2 ) + np.nanmean( kld['T'] ,2 ) )/3.0

      kld_series[it,:],nada,nada=bvf.get_regional_average(mean_kld,xi,xe,yi,ye)

      mean_kld = bvf.filter_field( mean_kld , smooth=smooth_type , sigma=smooth_sigma , cutoff=smooth_cutoff )

      mean_growth_rate = bvf.filter_field( mean_growth_rate , smooth=smooth_type , sigma=smooth_sigma , cutoff=smooth_cutoff )
 
      myname='growthrate_UVT_kld'

      print(np.min(mean_kld),np.max(mean_kld))

      clevels=np.arange(0.02,0.08,0.01)

      my_time=ctime.strftime("%Y%m%d%H%M%S")      

      bvf.plot_var_levels(mean_growth_rate,lon,lat,np.array([0]),plotbasedir,myname,varcontour=mean_kld,mycolorbar='YlGnBu',date=my_time,range='fixed',scale_min=0,scale_max=0.006,clevels=clevels) 

      myname='norm_UVT_kld'

      bvf.plot_var_levels(mean_norm_o,lon,lat,np.array([0]),plotbasedir,myname,varcontour=mean_kld,mycolorbar='YlGnBu',date=my_time,range='fixed',scale_min=0,scale_max=0.15,clevels=clevels)
 
      ctime = ctime + delta
 
      it = it + 1

print ( "Finish time loop" )

for ireg in range(0,nregs) :

   print(reg_name[ireg])
   fig, ax1 = plt.subplots()
   t = np.arange(0,60*(it-1),60)
   ax1.plot(t, norm_series[:t.size,ireg], 'b-',linewidth=3)
   ax1.plot(t, kld_series[:t.size,ireg], 'r-',linewidth=3)
   ax1.set_xlabel('time (s)')
   # Make the y-axis label, ticks and tick labels match the line color.
   ax1.set_ylabel('Norm and KLd', color='b')
   ax1.tick_params('y', colors='b')
   ax2 = ax1.twinx()
   ax2.plot(t, growthrate_series[:t.size,ireg], 'k--',linewidth=3)
   ax2.set_ylabel('Growth rate', color='k')
   ax2.tick_params('y', colors='k')
   fig.tight_layout()
   #plt.show()
   print( 'Generationg the following figure : ' + 'Figure_growthrate_vs_kld_series_' + reg_name[ireg] + '.png' )
   plt.savefig( plotbasedir + '/Figure_' + 'Figure_growthrate_vs_kld_series_' + reg_name[ireg] + '.png' )
   plt.close()






