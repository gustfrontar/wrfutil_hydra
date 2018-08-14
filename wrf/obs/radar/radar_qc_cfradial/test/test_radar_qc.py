#print __doc__

# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

# History:  Created 10-2017

#Add path to additional modules.
import sys
sys.path.append('./src/python' )
sys.path.append('./src/fortran')

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
options['ifattfilter']=False   #Attenuation filter
options['ifetfilter']=True     #Echo top filter
options['ifedfilter']=False    #Echo depth filter
options['ifspfilter']=False    #Speckle filter
options['ifblfilter']=False    #Blocking filter
options['ifmissfilter']=False  #Missing values filter

#General

options['ref_name']='dBZ'      #Reflectivity
options['cref_name']='CdBZ'    #Corrected reflectivity (qc output)
options['v_name']='V'          #Dopper velocity
options['cv_name']='CV'        #Corrected wind (qc ouput)
options['rho_name']='RhoHV'    #Rho HV

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

#Echo top parameters           #Filter layeers with echo top below a certain threshold.

options['etfilternx']=2        #Smooth parameter in x
options['etfilterny']=2        #Smooth parameter in y
options['etfilternz']=0        #Smooth parameter in z (dummy)
options['etfiltertr']=3000     #Echo top threshold.
options['etfilter_save']=True  #Wether echo top will be included in the output

#Echo depth parameters          #Filters layers with depths below a certain threshold.

options['edfilternx']=2         #Smooth parameter in x
options['edfilterny']=2         #Smooth parameter in y
options['edfilternz']=0         #Smooth parameter in z (dummy)
options['edfiltertr']=3000      #Echo top threshold.
options['edfilter_save']=True   #Wether echo top will be included in the output

#Speckle parameters

options['spfilternx']=2           #Box size in X NX=(2*spfilternx + 1)
options['spfilterny']=2           #Box size in Y
options['spfilternz']=0           #Box size in Z
options['spfilterreftr']=5        #Reflectivity threshold
options['spfiltertr']=0.3         #Count threshold
options['spfilter_save']=True     #Save filter fields.

#Attenuation parameters

options['attfiltertr']=20.0       #Attenuation threshold in dBZ
options['attcalerror']=1.0        #Calibration error
options['attfilter_save']=True    #Save filter fields

#Blocking parameters


#Detect missing parameters


#=======================================


# read in the file, create a RadarMapDisplay object
filename = '/media/jruiz/PAWR/Dropbox/DATA/DATOS_RADAR/PARANA/PAR_20091117_120/cfrad.20091117_210346.000_to_20091117_210735.000_PAR_SUR.nc'

#Performs QC operations based on options
[radar , qc_output] = rqc.main_qc( filename , options )

print('End of QC')

qc_output['echo_top'][ qc_output['echo_top'] == options['undef'] ] = 0.0

plt.figure()
plt.pcolor(qc_output['echo_top'][:,:,1])
plt.show()


