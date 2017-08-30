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

inibv=1          #Initial bred vector to plot.
endbv=1          #Final bred vector to plot.
niter=5          #End iter to plot.


nbv=endbv-inibv+1 #Total number of bred vectors.

undef_out=np.nan #This is the undef value that we will use internally.
undef_in=1.0e20  #This is the undef value in the original data (or at least a big number that is lower than the real undef value).

#Define regions

lati=np.array([34.75,34.6])
late=np.array([35.25,34.9])
loni=np.array([135.5,135.4])
lone=np.array([136.25,135.7])

reg_name='REG_1','REG_2','TOTAL'

plotlevels=np.array([6,13,17])   #Which levels will be plotted.
plotvars='UV','W','T','QV','QHYD' #Which variables will be plotted.

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

ctime = itime + delta

ntimes=1 + np.around((etime-itime).seconds / delta.seconds)

nx=180
ny=180
nz=20

data_pp_o=dict()
data_pn_o=dict()

data_pp_r=dict()
data_pn_r=dict()

bv_o=dict()
bv_r=dict()

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat_d01z001.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon_d01z001.grd',nx,ny,1,'>f4')[:,:,0]

time_mean_growth_rate=np.zeros([nx,ny,nz])
time_sprd_growth_rate=np.zeros([nx,ny,nz])
time_mean_norm=np.zeros([nx,ny,nz])
time_sprd_norm=np.zeros([nx,ny,nz])

#Convert lat lon to the nearest grid point.

#Add the global domain as a region.
lati=np.append(lati,lat[0,0])
late=np.append(late,lat[nx-1,ny-1])
loni=np.append(loni,lon[0,0])
lone=np.append(lone,lon[nx-1,ny-1])

xi , yi = bvf.lat_lon_to_i_j(lon,lat,loni,lati)
xe , ye = bvf.lat_lon_to_i_j(lon,lat,lone,late)

nregs=xi.shape[0]


norm_mean_o=dict()
norm_max_o=dict()
norm_min_o=dict()

norm_mean_r=dict()
norm_max_r=dict()
norm_min_r=dict()

norm_mean_i=dict()
norm_max_i=dict()
norm_min_i=dict()

gr_bv_mean=dict()
gr_bv_min=dict()
gr_bv_max=dict()


#Allocate memory for the dictionaries.
for myvar in plotvars :

    norm_mean_o[myvar]=np.zeros([ntimes,nbv,nregs])
    norm_max_o[myvar]=np.zeros([ntimes,nbv,nregs])  
    norm_min_o[myvar]=np.zeros([ntimes,nbv,nregs])

    norm_mean_r[myvar]=np.zeros([ntimes,nbv,nregs])
    norm_max_r[myvar]=np.zeros([ntimes,nbv,nregs])
    norm_min_r[myvar]=np.zeros([ntimes,nbv,nregs])

    norm_mean_i[myvar]=np.zeros([ntimes,nbv,nregs])
    norm_max_i[myvar]=np.zeros([ntimes,nbv,nregs])
    norm_min_i[myvar]=np.zeros([ntimes,nbv,nregs])


    gr_bv_mean[myvar]=np.zeros([ntimes,nbv,nregs])
    gr_bv_min[myvar]=np.zeros([ntimes,nbv,nregs])
    gr_bv_max[myvar]=np.zeros([ntimes,nbv,nregs])


for ibv in range (inibv , endbv + 1):

   bvstr="%04d" % ibv

   #print( ' Plotting bred vector number ' + bvstr )

   while ( ctime <= etime ):

    it = ( ctime - itime ).seconds / delta.seconds

    iterstr="%04d" % niter
 
    ptime=ctime - delta #Data correspinding to the previous step (to compute bv growth)
 
    print ( 'The date is :', ctime )

    print ( 'Reading the positive perturbation original')

    mydir=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pp_o' + iterstr + '/'
    data_pp_o=bio.read_data_scale(mydir,expname,ctime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)
    
    #print(mydir)

    print ( 'Reading the negative perturbation original')

    mydir=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pn_o' + iterstr + '/'
    data_pn_o=bio.read_data_scale(mydir,expname,ctime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)

    #print(mydir)

    print ( 'Reading the positive perturbation rescaled')

    mydir=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pp_r' + iterstr + '/'
    data_pp_i=bio.read_data_scale(mydir,expname,ptime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)

    #print(mydir)

    print ( 'Reading the negative perturbation rescaled')

    mydir=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pn_r' + iterstr + '/'
    data_pn_i=bio.read_data_scale(mydir,expname,ptime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)

    #print(mydir)

    iterstr="%04d" % 1  #Leo el rescaling de la perturbacion 1 en el tiempo t que es el rescaling de la perturbacion de la ultima iteracion del tiempo t.
 
    print ( 'Reading the positive perturbation rescaled')

    mydir=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pp_r' + iterstr + '/'
    data_pp_r=bio.read_data_scale(mydir,expname,ctime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)

    #print(mydir)

    print ( 'Reading the negative perturbation rescaled')

    mydir=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pn_r' + iterstr + '/'
    data_pn_r=bio.read_data_scale(mydir,expname,ctime,nx,ny,nz,undef_in=undef_in,undef_out=undef_out)

    #print(mydir)
    #bv_o=bvf.data_diff( data_pp_o , data_pn_o )
    #bv_r=bvf.data_diff( data_pp_o , data_pn_o )

    for my_var in plotvars  :

       norm_mean_o[my_var][it,ibv-1,:],norm_max_o[my_var][it,ibv-1,:],norm_min_o[my_var][it,ibv-1,:],norm_o=bvf.norm_bv( data_pp_o , data_pn_o , norm_type=my_var , smooth=smooth_type , sigma=smooth_sigma , xi=xi , yi=yi , xe=xe , ye=ye )
       norm_mean_r[my_var][it,ibv-1,:],norm_max_r[my_var][it,ibv-1,:],norm_min_r[my_var][it,ibv-1,:],norm_r=bvf.norm_bv( data_pp_r , data_pn_r , norm_type=my_var , smooth=smooth_type , sigma=smooth_sigma , xi=xi , yi=yi , xe=xe , ye=ye )
       norm_mean_i[my_var][it-1,ibv-1,:],norm_max_i[my_var][it-1,ibv-1,:],norm_min_i[my_var][it-1,ibv-1,:],norm_r=bvf.norm_bv( data_pp_i , data_pn_i , norm_type=my_var , smooth=smooth_type , sigma=smooth_sigma , xi=xi , yi=yi , xe=xe , ye=ye )


    ctime = ctime + delta

    ntimes = ntimes + 1


print ( "Finish time loop" )


mybv=0

mydir=plotbasedir + '/time_independent_plots/' + '/' + bvstr + '/' 
#Plot norm time series.
bvf.plot_norm_timeseries(norm_mean_o,norm_mean_i,norm_mean_r,plotvars,reg_name,mydir,mybv,'norm_mean',figsize='None') 














