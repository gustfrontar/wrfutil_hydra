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

expname = '/breeding_osaka_pawr_1km_bip/'

plotbasedir=basedir + expname + '/plots/'

nbredvector=1    #Total number of bred vectors.
inibv=1          #Initial bred vector to plot.
endbv=1          #Final bred vector to plot.
nbipiter=1       #Total number of iterations if breeding in place is activated.
iniit=1          #Initial iter to plot.
endit=1          #End iter to plot.

plotlevels=np.array([6,13,17])   #Which levels will be plotted.
plotvars=['U','V','W','T','QV']    #Which variables will be plotted.

norm_type='UVT'
smooth_type='Gaussian'
smooth_sigma=np.array([1.5])

#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,10,30)  #Initial time.
etime = dt.datetime(2013,7,13,5,31,00)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=30)

nx=180
ny=180
nz=20

data_pp_o=dict()
data_pn_o=dict()

data_pp_r=dict()
data_pn_r=dict()

bv_o=dict()
bv_r=dict()

ctime=itime + delta

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]


time_mean_growth_rate=np.zeros([nx,ny,nz])
time_sprd_growth_rate=np.zeros([nx,ny,nz])
time_mean_norm=np.zeros([nx,ny,nz])
time_sprd_norm=np.zeros([nx,ny,nz])

ntimes=0

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
      data_pp_o=bio.read_data_scale(mydir,expname,ctime,nx,ny,nz)

      print ( 'Reading the negative perturbation original')

      mydir=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pn_o' + iterstr + '/'
      data_pn_o=bio.read_data_scale(mydir,expname,ctime,nx,ny,nz)

      print ( 'Reading the positive perturbation rescaled')

      mydir=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pp_r' + iterstr + '/'
      data_pp_r=bio.read_data_scale(mydir,expname,ptime,nx,ny,nz)

      print ( 'Reading the negative perturbation rescaled')

      mydir=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + '/grads_pn_r' + iterstr + '/'
      data_pn_r=bio.read_data_scale(mydir,expname,ptime,nx,ny,nz)


      bv_o=bvf.data_diff( data_pp_o , data_pn_o )
      bv_r=bvf.data_diff( data_pp_o , data_pn_o )

      norm_o=bvf.norm_bv( data_pp_o , data_pn_o , norm_type=norm_type , smooth=smooth_type , sigma=smooth_sigma )

      norm_r=bvf.norm_bv( data_pp_r , data_pn_r , norm_type=norm_type , smooth=smooth_type , sigma=smooth_sigma )
 
      gr_bv_mean, gr_bv_max , gr_bv_min , gr_bv = bvf.growth_rate_bv( norm_o , norm_r , xi='None' , xe='None' , yi='None' , ye='None', zi='None', ze='None') 

      #Plot BV norm.
      mydir=plotbasedir + '/' + '/' + bvstr + '/' + iterstr + '/'  
      bvf.plot_var_levels( norm_o , lon , lat , plotlevels , mydir , 'norm' + norm_type , date=ctime.strftime("%Y%m%d%H%M%S") )
      bvf.plot_var_ave( norm_o , lon , lat , mydir , 'norm' , varcontour=np.array(data_pp_o['W']) , date=ctime.strftime("%Y%m%d%H%M%S") )

      #Plot BV growing rate.
      bvf.plot_var_levels( gr_bv , lon , lat , plotlevels , mydir , 'grrate' + norm_type , date=ctime.strftime("%Y%m%d%H%M%S") )
      bvf.plot_var_ave( gr_bv , lon , lat , mydir , 'grrate' , varcontour=np.array(data_pp_o['W']), date=ctime.strftime("%Y%m%d%H%M%S"))

      time_mean_growth_rate = time_mean_growth_rate + gr_bv
      time_sprd_growth_rate = time_sprd_growth_rate + np.power( gr_bv , 2 ) 
  
      time_mean_norm = time_mean_norm + norm_o
      time_sprd_norm = time_sprd_norm + np.power( norm_o , 2 )

      ctime = ctime + delta
 
      ntimes = ntimes + 1

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














