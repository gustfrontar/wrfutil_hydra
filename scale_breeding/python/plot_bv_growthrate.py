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

undef_in=1.0e20
undef_out=np.nan

qmeanlevs=[0.001,0.0050,0.05]   #Levels for vertically averaged condensate.

norm_type='UVT'
smooth_type='Gaussian'
smooth_sigma=np.array([1.5])

#Define regions

lati=np.array([34.75,34.6])
late=np.array([35.25,34.9])
loni=np.array([135.5,135.4])
lone=np.array([136.25,135.7])

reg_name='REG_1','REG_2','TOTAL'

plotlevels=np.array([6,13,17])   #Which levels will be plotted.

#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,35,00)  #Initial time.
etime = dt.datetime(2013,7,13,5,39,30)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=30)

ntimes=round( (itime-etime).seconds/delta.seconds ) + 1

nx=180
ny=180
nz=20

data_pp_o=dict()
data_pn_o=dict()

data_pp_r=dict()
data_pn_r=dict()

ctime=itime + delta

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat_d01z001.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon_d01z001.grd',nx,ny,1,'>f4')[:,:,0]

#Add the global domain as a region.
lati=np.append(lati,lat[0,0])
late=np.append(late,lat[nx-1,ny-1])
loni=np.append(loni,lon[0,0])
lone=np.append(lone,lon[nx-1,ny-1])

xi , yi = bvf.lat_lon_to_i_j(lon,lat,loni,lati)
xe , ye = bvf.lat_lon_to_i_j(lon,lat,lone,late)

nregs=xi.shape[0]

time_mean_growth_rate=np.zeros([nx,ny,nz])
time_sprd_growth_rate=np.zeros([nx,ny,nz])
time_mean_norm=np.zeros([nx,ny,nz])
time_sprd_norm=np.zeros([nx,ny,nz])

norm_mean_o=np.zeros([ntimes,nregs])
norm_mean_r=np.zeros([ntimes,nregs])
norm_max_o=np.zeros([ntimes,nregs])
norm_max_r=np.zeros([ntimes,nregs])
norm_min_o=np.zeros([ntimes,nregs])
norm_min_r=np.zeros([ntimes,nregs])

gr_bv_mean=np.zeros([ntimes,nregs])
gr_bv_max=np.zeros([ntimes,nregs])
gr_bv_min=np.zeros([ntimes,nregs])

it=0

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

      #Note: In pithon when multiple variables are output from a function in a tuple, then all the variables has to be decodified. 
      #If not the reconstruction of the variables will fail.

      norm_mean_o[it,:] , norm_max_o[it,:] , norm_min_o[it,:] , norm_o =bvf.norm_bv( data_pp_o , data_pn_o , norm_type=norm_type , smooth=smooth_type , sigma=smooth_sigma , xi=xi,yi=yi,xe=xe,ye=ye)
      norm_mean_r[it,:] , norm_max_r[it,:] , norm_min_r[it,:] , norm_r =bvf.norm_bv( data_pp_r , data_pn_r , norm_type=norm_type , smooth=smooth_type , sigma=smooth_sigma , xi=xi,yi=yi,xe=xe,ye=ye)

      gr_bv_mean[it,:] , gr_bv_max[it,:] , gr_bv_min[it,:] , gr_bv = bvf.growth_rate_bv( norm_o , norm_r , xi=xi , xe=xe , yi=yi , ye=ye ) 

      #Plot BV norm.
      mydir=plotbasedir + '/' + '/' + bvstr + '/' + iterstr + '/'

      bvf.plot_var_levels( norm_o , lon , lat , plotlevels , mydir , 'norm' + norm_type , date=ctime.strftime("%Y%m%d%H%M%S") )
      bvf.plot_var_ave( norm_o , lon , lat , mydir , 'norm' , varcontour=data_pp_o['QHYD'] , clevels=qmeanlevs , date=ctime.strftime("%Y%m%d%H%M%S") )

      #Plot BV growing rate.
      bvf.plot_var_levels( gr_bv , lon , lat , plotlevels , mydir , 'grrate' + norm_type , date=ctime.strftime("%Y%m%d%H%M%S") )
      bvf.plot_var_ave( gr_bv , lon , lat , mydir , 'grrate' , varcontour=data_pp_o['QHYD'] , clevels=qmeanlevs , date=ctime.strftime("%Y%m%d%H%M%S"))

      time_mean_growth_rate = time_mean_growth_rate + gr_bv
      time_sprd_growth_rate = time_sprd_growth_rate + np.power( gr_bv , 2 ) 
  
      time_mean_norm = time_mean_norm + norm_o
      time_sprd_norm = time_sprd_norm + np.power( norm_o , 2 )

      ctime = ctime + delta
 
      it = it + 1

print ( "Finish time loop" )

time_mean_growth_rate = time_mean_growth_rate / ntimes

time_mean_norm = time_mean_norm / ntimes

time_sprd_growth_rate = np.power( time_sprd_growth_rate / ntimes - np.power( time_mean_growth_rate , 2 ) , 0.5)

time_sprd_norm = np.power( time_sprd_norm / ntimes - np.power( time_mean_norm , 2 ) , 0.5)

#Plot mean norm
mydir=plotbasedir + '/time_independent_plots/' + '/' + bvstr + '/' + iterstr + '/'
bvf.plot_var_levels( time_mean_norm , lon , lat , plotlevels , mydir , 'tmean_norm' + norm_type )
bvf.plot_var_ave( time_mean_norm , lon , lat , mydir , 'tmean_norm' )

#Plot mean growing rate
bvf.plot_var_levels( time_mean_growth_rate , lon , lat , plotlevels , mydir , 'tmean_grrate' + norm_type )
bvf.plot_var_ave( time_mean_growth_rate , lon , lat , mydir , 'tmean_grrate' )

#Plot std norm
mydir=plotbasedir + '/time_independent_plots/' + '/' + bvstr + '/' + iterstr + '/'
bvf.plot_var_levels( time_sprd_norm , lon , lat , plotlevels , mydir , 'tstd_norm' + norm_type )
bvf.plot_var_ave( time_sprd_norm , lon , lat , mydir , 'tstd_norm' )

#Plot std growing rate
bvf.plot_var_levels( time_sprd_growth_rate , lon , lat , plotlevels , mydir , 'tstd_grrate' + norm_type  )
bvf.plot_var_ave( time_sprd_growth_rate , lon , lat , mydir , 'tstd_grrate' )


#TODO plot time evolution of norm and growing rate for selected regions. (Define regions based on different norm behaviours)














