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

basedir='/data9/jruiz/EXPERIMENTS/'

expname = '/OsakaPAR_1km_control1000m_smallrandompert_new/'

plotbasedir=basedir + expname + '/plots/'

inipert=1          #Initial perturbation to plot.
endpert=1          #Final perturbation to plot.
npert=endpert-inipert + 1

undef_out=np.nan #This is the undef value that we will use internally.
undef_in=1.0e20  #This is the undef value in the original data (or at least a big number that is lower than the real undef value).

#The following will be used to extract a particlar variable from the original data.
#This variables should be especified according to the data that we have in the binary files.

ctl_vars='U','V','W','T','QV','QHYD'  #Complete list of variables in ctl file.
ctl_inirecord=[0,12,24,36,48,60]        #Starting record for each variable. From 0 to N
ctl_endrecord=[11,23,35,47,59,71]       #End record for each variable. From 0 to N.

#Which variables and levels are we going to plot?

plotlevels=np.array([3,7,9])   #Which levels will be plotted (this levels are equivalent to the BV plots)
plotvars='UV','V','W','T','QV','QHYD'    #Which variables will be plotted.

#Define regions

lati=np.array([34.75,34.6])
late=np.array([35.25,34.9])
loni=np.array([135.5,135.4])
lone=np.array([136.25,135.7])

reg_name='REG_1','REG_2','TOTAL'

smooth_type='None'
smooth_sigma=np.array([1.5])

#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,10,30)  #Initial time.
etime = dt.datetime(2013,7,13,5,40,00)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=30)

nx=180
ny=180
nz=np.max(ctl_endrecord) + 1  #Total number of records in binary file.
nlev=12                       #Number of vertical levels for 3D variables.
ntimes=1 + np.around((etime-itime).seconds / delta.seconds).astype(int)  #Total number of times.

ctime = itime + delta


data_pert_gues=dict()
data_mean_gues=dict()

data_pert_anal=dict()
data_mean_anal=dict()


#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]

time_mean_growth_rate=np.zeros([nx,ny,nlev])
time_sprd_growth_rate=np.zeros([nx,ny,nlev])
time_mean_norm=np.zeros([nx,ny,nlev])
time_sprd_norm=np.zeros([nx,ny,nlev])

#Convert lat lon to the nearest grid point.

#Add the global domain as a region.
lati=np.append(lati,lat[0,0])
late=np.append(late,lat[nx-1,ny-1])
loni=np.append(loni,lon[0,0])
lone=np.append(lone,lon[nx-1,ny-1])

xi , yi = bvf.lat_lon_to_i_j(lon,lat,loni,lati)
xe , ye = bvf.lat_lon_to_i_j(lon,lat,lone,late)

nregs=xi.shape[0]

norm_mean_gues=dict()
norm_max_gues=dict()
norm_min_gues=dict()

norm_mean_anal=dict()
norm_max_anal=dict()
norm_min_anal=dict()

gr_letkfpert_mean=dict()
gr_letkfpert_min=dict()
gr_letkfpert_max=dict()


#Allocate memory for the dictionaries.
for myvar in plotvars :

    norm_mean_gues[myvar]=np.zeros([ntimes,npert,nregs])
    norm_max_gues[myvar]=np.zeros([ntimes,npert,nregs])  
    norm_min_gues[myvar]=np.zeros([ntimes,npert,nregs])

    norm_mean_anal[myvar]=np.zeros([ntimes,npert,nregs])
    norm_max_anal[myvar]=np.zeros([ntimes,npert,nregs])
    norm_min_anal[myvar]=np.zeros([ntimes,npert,nregs])

    gr_letkfpert_mean[myvar]=np.zeros([ntimes,npert,nregs])
    gr_letkfpert_min[myvar]=np.zeros([ntimes,npert,nregs])
    gr_letkfpert_max[myvar]=np.zeros([ntimes,npert,nregs])


for ipert in range (inipert , endpert + 1):

   pertstr="%04d" % ipert

   #print( ' Plotting bred vector number ' + bvstr )

   while ( ctime <= etime ):

    it = ( ctime - itime ).seconds / delta.seconds
 
    print ( 'The date is :', ctime )

    print ( 'Reading the perturbed gues' )

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/' + pertstr + '.grd'

    data_pert_gues=bio.read_data_scale_2( my_file , nx , ny , nz , ctl_vars , ctl_inirecord , ctl_endrecord ,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

    print ( 'Reading the mean gues')

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/mean.grd'

    data_mean_gues=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

    print ( 'Reading the perturbed anal' )

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/' + pertstr + '.grd'

    data_pert_anal=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

    print ( 'Reading the mean anal')

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/mean.grd'

    data_mean_anal=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

    for my_var in plotvars  :

       norm_mean_gues[my_var][it,ipert-1,:],norm_max_gues[my_var][it,ipert-1,:],norm_min_gues[my_var][it,ipert-1,:],norm_gues=bvf.norm_bv( data_pert_gues , data_mean_gues , norm_type=my_var , smooth=smooth_type , sigma=smooth_sigma , xi=xi , yi=yi , xe=xe , ye=ye )
       norm_mean_anal[my_var][it,ipert-1,:],norm_max_anal[my_var][it,ipert-1,:],norm_min_anal[my_var][it,ipert-1,:],norm_anal=bvf.norm_bv( data_pert_anal , data_mean_anal , norm_type=my_var , smooth=smooth_type , sigma=smooth_sigma , xi=xi , yi=yi , xe=xe , ye=ye )

    ctime = ctime + delta

    it = it + 1

print ( "Finish time loop" )


mybv=0

mydir=plotbasedir + '/time_independent_plots/' + '/' + pertstr + '/' 
#Plot norm time series.
#bvf.plot_norm_timeseries(norm_mean_o,norm_mean_i,norm_mean_r,plotvars,reg_name,mydir,mybv,'norm_mean',figsize='None') 
plot_reg_name = True
figsize=bvf.default_figure_size
figextension='.png'

#Create output directory
if not os.path.exists(plotbasedir):
  os.makedirs(plotbasedir)

#Create time series.
time_serie=np.nan*np.zeros([ntimes*2,npert,nregs])
times=np.nan*np.zeros(ntimes*2)

for myvar in plotvars     :
 for ii in range(0,ntimes-1)  :
  for ipert in range(inipert,endpert+1)  :
   #Creamos un grafico tipo serrucho para representar la evolucion de la norma.
   time_serie[2*ii,ipert-1,:]=norm_mean_anal[myvar][ii,ipert-1,:]
   time_serie[2*ii+1,ipert-1,:]=norm_mean_gues[myvar][ii+1,ipert-1,:]

   times[2*ii]=ii
   times[2*ii+1]=ii+1

 for ireg in range(0,nregs)      :
  iregstr="%04d" % ( ireg + 1 )

  fig=plt.figure(1,figsize=figsize)

  for ipert in range(inipert,endpert+1)             :
   plt.plot(times,time_serie[:,ipert-1,ireg],'-')

   plt.ylabel('Norm')
   plt.xlabel('Time')

   #if debug == True   :
   plt.show()

   print( 'Generationg the following figure : ' + 'Figure_' + '_pertamptimeserie_' + myvar + 'reg' + iregstr + figextension )
   plt.savefig( plotbasedir + 'Figure_' + '_pertamptimeserie_' + myvar + 'reg' + iregstr + figextension )

   plt.close(fig)















