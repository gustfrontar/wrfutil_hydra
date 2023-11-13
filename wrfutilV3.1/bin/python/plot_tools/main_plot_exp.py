import numpy as np
from netCDF4 import Dataset
from wrf import getvar
import glob
import mod_plot as mp

#Get env variables and call functions accordingly. 
conf=dict()

#conf['exp_type']      = os.getenv( 'PLOT_TYPE' ) #FCST,ANAL,DAFCST
#conf['exp_ini_date']   = os.getenv( 'FECHA_INI' )
#conf['exp_end_date']   = os.getenv( 'FECHA_FIN' )
#conf['basedir']        = os.getenv( 'BASEDIR' )   #Base directory of the experiment.
#conf['histdir']        = os.getenv( 'HISTDIR' )   #Directory containing the experiment data.
#conf['plotdir']        = os.getenv( 'PLOTDIR' )   #Directory where plots will be created.
#conf['plot_type_list'] = os.getenve( 'PLOT_TYPE_LIST' ) #TODO create a plot.conf file and add this to this file.

conf['exp_type']       = ['ANAL']
conf['exp_ini_date']   = '2019-10-10 18:00:00'
conf['exp_end_date']   = '2019-10-12 18:00:00'
conf['basedir']        = '/home/ra000007/a04037/data/data_assimilation_exps/PREVENIR_HIRESOLUTION_20191010/'
conf['histdir']        = conf['basedir'] + '/HIST/' 
conf['plotdir']        = conf['basedir'] + '/PLOT/'
conf['plot_type_list'] = ['MAXDBZ']

#Plot types
plot_types=dict()
plot_types['MAXDBZ'] = {'pvar1':'mdbz','pvar1dim':2,'pvar1vlev':None,'pvar2':'W','pvar2dim':3,'pvar2vlev':np.array([5.0]),'pvar1cmap':'gist_ncar','pvar1clev':np.arange(5., 75., 5.),'pvar2clev':np.arange(5.,75.,10.)}


  
mp.ploter_loop( conf , plot_types ) 
