#!/share/anaconda3/bin/python

"""
Created on 07/11/2017
"""

import sys
sys.path.append('../../../common_python/common_modules/')
sys.path.append('../../../common_python/common_functions/')
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import binary_io as bio
import os
from common_functions import common_functions as comm
from common_functions import common_data      as cd


basedir='/home/jruiz/salidas/EXPERIMENTS/'

expname = 'ANALYSIS_PARANA_2KM_control_paranafnl_newobs_60m_radar_grib_Hydra'

undef_in=1.0e30


#Defini initial and end times using datetime module.
itime = dt.datetime(2009,11,17,18,5,00)  #Initial time.
etime = dt.datetime(2009,11,17,23,5,00)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=300)

nx=250                        #Number of grid poinst in x
ny=250                        #Number of grid points in y
nz=226                        #Total number of records in binary file.
ne=60                         #Number of ensemble members

ntimes=1 + np.around((etime-itime).seconds / delta.seconds)  #Total number of times.


for source in ['anal','gues'] :

	ctime=itime

	it=0

	while ( ctime <= etime ):       #Loop over the times.

		print( ctime )

		my_data_mean=np.zeros((nx,ny,nz))
		my_data_sprd=np.zeros((nx,ny,nz))
		my_data_coun=np.zeros((nx,ny,nz))

		for im in range( 0 , ne ) :   #Loop over the members
			imstr="%05d" % ( im + 1 ) 
			my_file=basedir + '/' + expname + '/' + source + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/plev' +  imstr + '.dat'

			print ( 'Reading file: ' + my_file)

			my_data=bio.read_data_direct(my_file,nx,ny,nz,dtypein='>f4',undef_in=undef_in,undef_out=undef,seq_acces=False)

                        #Add data to the ensemble array.
			comm.add_to_ensemble(nx=nx,ny=ny,nz=nz,ibv=im+1,my_data=my_data )

			#my_data_coun[  my_data != undef_out  ]= my_data_coun[ my_data != undef_out ] + 1

			#my_data[  my_data == undef_out  ] = 0.0

			#my_data_mean=my_data_mean + my_data 
			#my_data_sprd=my_data_sprd + np.power( my_data , 2 )

                #Generate a mask for the ensemble.
		comm.getmask(nx=nx,ny=ny,nz=nz,nbv=enssize,undef=undef)
 
		#Compute the mean and the standard deviation.
		my_moments=comm.compute_moments(nx=nx,ny=ny,nz=nz,nbv=ne,nmoments=2,undef=undef)
		#Compute ensemble mean
		#my_data_mean[ my_data_coun > 0 ]=my_data_mean[ my_data_coun > 0 ] / my_data_coun[ my_data_coun > 0 ]
		#Compute ensemble spread
		#my_data_sprd[ my_data_coun > 0 ]=np.power( my_data_sprd[ my_data_coun > 0] / my_data_coun[ my_data_coun > 0 ] - np.power(my_data_mean[ my_data_coun > 0 ],2) , 0.5 )


		#my_data_mean[ my_data_coun == 0 ]=undef_in
		#my_data_sprd[ my_data_coun == 0 ]=undef_in

		#Write the output
		my_file_mean=basedir + '/' + expname + '/' + source + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/plev_mean.dat'
		my_file_sprd=basedir + '/' + expname + '/' + source + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/plev_sprd.dat'

		bio.write_data_direct_woundef(my_file_mean,my_moments[:,:,:,0],dtypein='>f4')
		bio.write_data_direct_woundef(my_file_sprd,my_moments[:,:,:,1],dtypein='>f4')
		ctime = ctime + delta
 

	print ( 'Finish time loop for source ' + source )

print('Normal end')




