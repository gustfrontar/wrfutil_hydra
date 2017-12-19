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
from common_functions import common_functions as comm
from common_functions import common_data      as cd


basedir='/data9/jruiz/EXPERIMENTS/'

expname = '/OsakaPAR_1km_control1000m_smallrandompert_new/'

filetype='analgp'   #analgp , analgz , guesgp , guesgz

smooth=False #0- no smooth , 1 apply smooth
smooth_lambda=10 #Smooth length scale (in number of grid points)

undef=bio.default_undef_val

input_undef_val=1e30



enssize=1000     #Total number of bred vectors.

#Estimate number of bins in histogram. 
nbins=np.round(np.sqrt(enssize)) #Number of bins to be used in histogram computation.

nmoments=4       #Cantidad total de momentos que vamos a calcular.

get_kldistance=True #Wether we compute or not the Kullback-Leiber distance.

get_histogram=True  #Wether if histogram will be explicitelly calculated and stored.

get_moments=True    #Wether we will be computeing ensemble moments.

#plotlevels=np.array([6,13,17])   #Which levels will be plotted.
#plotvars=['U','V','W','T','QV','QHYD']    #Which variables will be plotted.

#The following will be used to extract a particlar variable from the original data.
#This variables should be especified according to the data that we have in the binary files.

ctl_vars=['U','V','W','T','QV','QHYD']  #Complete list of variables in ctl file.
ctl_inirecord=[0,12,24,36,48,60]        #Starting record for each variable. From 0 to N
ctl_endrecord=[11,23,35,47,59,71]       #End record for each variable. From 0 to N.

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,11,00)  #Initial time.
etime = dt.datetime(2013,7,13,5,39,00)  #End time.

ctime = itime

#Define the delta.
delta=dt.timedelta(seconds=60)  #Original data is every 30 seconds 

nx=180
ny=180
nz=np.max(ctl_endrecord) + 1  #Total number of records in binary file.
nlev=12                       #Number of vertical levels for 3D variables.


#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]

#my_ensemble=np.zeros([nx,ny,nz,enssize])
tmp=np.zeros([nx,ny,nz])

comm.allocate_ensemble(nx=nx,ny=ny,nz=nz,nbv=enssize) #Allocate the ensemble and the undefmask.

while ( ctime <= etime ):

  print(ctime)

  #Ensemble reading loop.

  for imem in range (0, enssize):

      memstr="%04d" % ( imem + 1 )

      print( ' Reading ensemble member ' + memstr )
 
      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/' + memstr + '.grd'
      #my_ensemble[:,:,:,imem]=bio.read_data_direct_woundef(my_file,nx,ny,nz,'f4').astype('float32')
      comm.add_to_ensemble(nx=nx,ny=ny,nz=nz,ibv=imem+1,my_data=bio.read_data_direct(my_file,nx,ny,nz,'f4',undef_in=input_undef_val,undef_out=undef).astype('float32') )

  #Get the mask.
  comm.getmask(nx=nx,ny=ny,nz=nz,nbv=enssize,undef=undef)

  if get_moments :
      print("Computing pdf moments")
      my_moments=comm.compute_moments(nx=nx,ny=ny,nz=nz,nbv=enssize,nmoments=nmoments,undef=undef) 

      for imoment in range (0,nmoments)  :
        momentstr="%04d" % ( imoment + 1 )
        my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/moment' + momentstr + '.grd'
        bio.write_data_direct_woundef(my_file,my_moments[:,:,:,imoment],'f4')

  if get_kldistance    :
       print("Computeing kl distance") 
       kldist=comm.compute_kld(nx=nx,ny=ny,nz=nz,nbv=enssize,undef=undef,nbins=nbins)
       my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/kldistance.grd'
       bio.write_data_direct_woundef(my_file,kldist,'f4')

  if get_histogram     :
       print("Computing histogram")
   #Compute histogram explicitelly
       ensmin,ensmax,histogram=comm.compute_histogram(nx=nx,ny=ny,nz=nz,nbv=enssize,nbins=nbins,undef=undef)
       #Write ensemble minimum
       my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/ensmin.grd' 
       bio.write_data_direct_woundef(my_file,ensmin,'f4')
       #Write ensemble maximum
       my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/ensmax.grd'
       bio.write_data_direct_woundef(my_file,ensmax,'f4')
       #Write histogram (check!)
       my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/histogram.grd'
       histogram=np.reshape(histogram,[nx,ny,nz*nbins])
       bio.write_data_direct_woundef(my_file,ensmax,'i2') #Using low precission integer to store histogram.

  ctime = ctime + delta

print ( "Finish time loop" )





