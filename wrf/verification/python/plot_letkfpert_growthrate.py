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

basedir='/data9/jruiz/EXPERIMENTS/'

expname = '/OsakaPAR_1km_control1000m_smallrandompert_new/'

plotbasedir=basedir + expname + '/plots/'

undef_in=1.0e20
undef_out=np.nan

qmeanlevs=[0.001,0.0050,0.05]   #Levels for vertically averaged condensate.

inipert=1   #Which is the first perturbation that we will plot.
endpert=1   #Which is the last perturbation that we will plot.
npert=endpert-inipert+1 #Total number of perturbations that will be plotted.

norm_type='UVT'
smooth_type='Gaussian'
smooth_sigma=np.array([1.5])

#The following will be used to extract a particlar variable from the original data.
#This variables should be especified according to the data that we have in the binary files.

ctl_vars='U','V','W','T','QV','QHYD'  #Complete list of variables in ctl file.
ctl_inirecord=[0,12,24,36,48,60]        #Starting record for each variable. From 0 to N
ctl_endrecord=[11,23,35,47,59,71]       #End record for each variable. From 0 to N.

#Which variables and levels are we going to plot?

plotlevels=np.array([3,7,9])   #Which levels will be plotted (this levels are equivalent to the BV plots)
plotvars='U','V','W','T','QV','QHYD'    #Which variables will be plotted.

#Define regions

lati=np.array([34.75,34.6])
late=np.array([35.25,34.9])
loni=np.array([135.5,135.4])
lone=np.array([136.25,135.7])

reg_name='REG_1','REG_2','TOTAL'

#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,10,30)  #Initial time.
etime = dt.datetime(2013,7,13,5,39,30)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=30)

ntimes=round( (itime-etime).seconds/delta.seconds ) + 1

nx=180
ny=180
nz=np.max(ctl_endrecord) + 1  #Total number of records in binary file.
nlev=12                       #Number of vertical levels for 3D variables.
ntimes=round( (itime-etime).seconds/delta.seconds ) + 1 #Total number of times.

data_pert_anal=dict()
data_mean_anal=dict()

data_pert_gues=dict()
data_mean_gues=dict()

ctime=itime 

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]

#Add the global domain as a region.
lati=np.append(lati,lat[0,0])
late=np.append(late,lat[nx-1,ny-1])
loni=np.append(loni,lon[0,0])
lone=np.append(lone,lon[nx-1,ny-1])

xi , yi = bvf.lat_lon_to_i_j(lon,lat,loni,lati)
xe , ye = bvf.lat_lon_to_i_j(lon,lat,lone,late)

nregs=xi.shape[0]

time_mean_growth_rate=np.zeros([nx,ny,nlev])
time_sprd_growth_rate=np.zeros([nx,ny,nlev])
time_mean_norm=np.zeros([nx,ny,nlev])
time_sprd_norm=np.zeros([nx,ny,nlev])

norm_mean_gues=np.zeros([ntimes,npert,nregs])
norm_mean_anal=np.zeros([ntimes,npert,nregs])
norm_max_gues=np.zeros([ntimes,npert,nregs])
norm_max_anal=np.zeros([ntimes,npert,nregs])
norm_min_gues=np.zeros([ntimes,npert,nregs])
norm_min_anal=np.zeros([ntimes,npert,nregs])

gr_pert_mean=np.zeros([ntimes,npert,nregs])
gr_pert_max=np.zeros([ntimes,npert,nregs])
gr_pert_min=np.zeros([ntimes,npert,nregs])

norm_gues=np.zeros([nx,ny,nlev])
norm_anal=np.zeros([nx,ny,nlev])

int_liquid=np.zeros([nx,ny,nlev])

it=0

for ipert in range (inipert , endpert + 1):

   pertstr="%04d" % ipert

   print( ' Plotting bred vector number ' + pertstr )

   while ( ctime <= etime ):
 
      ptime=ctime - delta #Data correspinding to the previous step (to compute bv growth)
 
      print ( 'The date is :', ctime )

      print ( 'Reading the perturbed analysis')
  
      my_file=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/' + pertstr + '.grd'

      data_pert_anal=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

      print ( 'Reading the analysis mean')

      my_file=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/mean.grd'

      data_mean_anal=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

      print ( 'Reading the perturbed gues')

      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/' + pertstr + '.grd'

      data_pert_gues=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

      print ( 'Reading the gues mean')

      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/mean.grd'

      data_mean_gues=bio.read_data_scale_2(my_file,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='f4',undef_in=undef_in,undef_out=undef_out)

      #Compute total integrated liquid (we will use this to identify areas associated with clouds and convection)
      tmp_int_liquid = np.nansum(data_mean_anal['QHYD'],2)

      for ilev in range(0,nlev)   :          #Create a fake 3D array for the vertically integrated liquid
                                              #This is because the plotting function expects a 3D array as input.
        int_liquid[:,:,ilev]=tmp_int_liquid

      #Note: In pithon when multiple variables are output from a function in a tuple, then all the variables has to be decodified. 
      #If not the reconstruction of the variables will fail.

      norm_mean_gues[it,ipert-1,:] , norm_max_gues[it,ipert-1,:] , norm_min_gues[it,ipert-1,:] , norm_gues =bvf.norm_bv( data_pert_gues , data_mean_gues , norm_type=norm_type , smooth=smooth_type , sigma=smooth_sigma , xi=xi,yi=yi,xe=xe,ye=ye)
      norm_mean_anal[it,ipert-1,:] , norm_max_anal[it,ipert-1,:] , norm_min_anal[it,ipert-1,:] , norm_anal =bvf.norm_bv( data_pert_anal , data_mean_anal , norm_type=norm_type , smooth=smooth_type , sigma=smooth_sigma , xi=xi,yi=yi,xe=xe,ye=ye)

      gr_pert_mean[it,ipert-1,:] , gr_pert_max[it,ipert-1,:] , gr_pert_min[it,ipert-1,:] , gr_pert = bvf.growth_rate_bv( norm_gues , norm_anal , xi=xi , xe=xe , yi=yi , ye=ye ) 

      #Plot LETKF perturbation norm.
      mydir=plotbasedir + '/' + pertstr + '/'

      varname='norm_' + norm_type
      my_range='centered'

      bvf.plot_var_levels( norm_gues , lon , lat , plotlevels , mydir , varname , date=ctime.strftime("%Y%m%d%H%M%S") ,varcontour=int_liquid,clevels=qmeanlevs,range=my_range)
      bvf.plot_var_ave( norm_gues , lon , lat , mydir , varname , varcontour=data_mean_anal['QHYD'] , clevels=qmeanlevs , date=ctime.strftime("%Y%m%d%H%M%S"),range=my_range)

      varname='gr_' + norm_type
      my_range='centered'

      #Plot LETKF perturbation growth rate.
      bvf.plot_var_levels( gr_pert , lon , lat , plotlevels , mydir , varname , date=ctime.strftime("%Y%m%d%H%M%S") ,varcontour=int_liquid,clevels=qmeanlevs,range=my_range) 
      bvf.plot_var_ave( gr_pert , lon , lat , mydir , varname , varcontour=data_mean_anal['QHYD'] , clevels=qmeanlevs , date=ctime.strftime("%Y%m%d%H%M%S"),range=my_range)

      time_mean_growth_rate = time_mean_growth_rate + gr_pert
      time_sprd_growth_rate = time_sprd_growth_rate + np.power( gr_pert , 2 ) 
  
      time_mean_norm = time_mean_norm + norm_gues
      time_sprd_norm = time_sprd_norm + np.power( norm_gues , 2 )

      ctime = ctime + delta
 
      it = it + 1

print ( "Finish time loop" )

time_mean_growth_rate = time_mean_growth_rate / ntimes

time_mean_norm = time_mean_norm / ntimes

time_sprd_growth_rate = np.power( time_sprd_growth_rate / ntimes - np.power( time_mean_growth_rate , 2 ) , 0.5)

time_sprd_norm = np.power( time_sprd_norm / ntimes - np.power( time_mean_norm , 2 ) , 0.5)

#Plot mean norm
mydir=plotbasedir + '/time_independent_plots/' + '/' + pertstr + '/'
bvf.plot_var_levels( time_mean_norm , lon , lat , plotlevels , mydir , 'tmean_norm' + norm_type )
bvf.plot_var_ave( time_mean_norm , lon , lat , mydir , 'tmean_norm' + norm_type )

#Plot mean growing rate
bvf.plot_var_levels( time_mean_growth_rate , lon , lat , plotlevels , mydir , 'tmean_grrate' + norm_type )
bvf.plot_var_ave( time_mean_growth_rate , lon , lat , mydir , 'tmean_grrate' + norm_type )

#Plot std norm
mydir=plotbasedir + '/time_independent_plots/' + '/' + pertstr + '/'
bvf.plot_var_levels( time_sprd_norm , lon , lat , plotlevels , mydir , 'tstd_norm' + norm_type )
bvf.plot_var_ave( time_sprd_norm , lon , lat , mydir , 'tstd_norm' + norm_type )

#Plot std growing rate
bvf.plot_var_levels( time_sprd_growth_rate , lon , lat , plotlevels , mydir , 'tstd_grrate' + norm_type  )
bvf.plot_var_ave( time_sprd_growth_rate , lon , lat , mydir , 'tstd_grrate' + norm_type )




