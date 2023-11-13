import numpy as np
from netCDF4 import Dataset
from wrf import getvar

PLOT_TYPE = os.getenv( 'PLOT_TYPE' ) #FCST,ANAL,DAFCST

EXP_INI_DATE= os.getenv( 'FECHA_INI' )
EXP_END_DATE= os.getenv( 'FECHA_FIN' )

PLOT_BASE_PATH = os.getenv( 'PLOT_BASE_PATH' )







