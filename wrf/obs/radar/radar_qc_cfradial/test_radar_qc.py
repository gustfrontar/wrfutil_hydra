#print __doc__

# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

# History:  Created 10-2017

import numpy as np
import matplotlib.pyplot as plt
import pyart
import netCDF4
import radar_qc_module as rqc

#======================================
# QC PARAMETERS
#======================================
options = {}

#Flags
options['ifdealias']=False
options['ifrhofilter']=False   #Rhohv filter
options['ifetfilter']=True     #Echo top filter
options['ifspfilter']=True     #Speckle filter
options['ifattfileter']=True   #Attenuation filter
options['ifblfilter']=True     #Blocking filter
options['ifmissfilter']=True   #Missing values filter

#General

options['reflectivity_var_name']='dBZ'            #Reflectivity
options['corrected_reflectivity_var_name']='CdBZ' #Corrected reflectivity (qc output)
options['dv_var_name']='V'                        #Dopper velocity
options['dvda_var_name']='VDA'                    #Dealiased doppler velocity
options['corrected_dv_var_name']='CV'             #Corrected wind (qc ouput)

options['norainrefval']=-0.1
options['undef']=-9.99e9

#Dealiasing parameters (pyart)

options['interval_split']=3
options['skip_between_rays']=10
options['skip_along_rays']=10

#Rho filter parameters

options['rhofilternx']=2
options['rhofilterny']=2
options['rhofilternz']=0
options['rhofiltertr']=0.5

#Echo top parameters

options['etfilternx']=2
options['etfilterny']=2
options['etfilternz']=0
options['etfiltertr']=0.5

#Speckle parameters

options['spfilternx']=2 
options['spfilterny']=2
options['spfilternz']=0
options['spfilterreftr']=5
options['spfiltertr']=3

#Attenuation parameters


#Blocking parameters


#Detect missing parameters


#=======================================


# read in the file, create a RadarMapDisplay object
filename = '/home/jruiz/share/DATA/OBS/OBS_REAL_PARANA_20091117_CFRADIAL/cfrad.20091117_200345.000_to_20091117_200734.001_PAR_SUR.nc3'

#Performs QC operations based on options
radar = rqc.main_qc( filename , options )


#[ref_array , ref_az , ref_level , ref_time , ref_az_exact]=rqc.order_variable( radar , options['reflectivity_var_name'] )

#ref_array[np.isnan(ref_array)]=0.0

#plt.pcolor(ref_array[:,:,5])

#plt.show()
