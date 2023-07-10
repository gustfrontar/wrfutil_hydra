import numpy as np
import pandas as pd
import os
from datetime import datetime, timedelta
import auxiliary_functions_DF as af
import glob

FECHA_INI = os.environ['FECHA_INI']
CICLO = os.environ['CICLO']
BASEDIR = os.environ['BASEDIR']
DIR_COEF = BASEDIR + '/CALIB/horarios/'
MODELO = os.environ['FILEIN_TYPE']
PLAZOS = int(os.environ['PLAZO'])
WRFUTILDIR = os.environ['WRFUTILDIR']
DIR_TEMPLATES = WRFUTILDIR + "/templates/"

#OUTDIR = DIR_COEF + '/iniciales/' + CICLO + '/'
OUTDIR = DIR_COEF + CICLO + '/'

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

FECHA_INI_obj = datetime.strptime(FECHA_INI, '%Y/%m/%d')
ayer = FECHA_INI_obj - timedelta(days = 1)

if MODELO == 'WRF':
    FILEIN = glob.glob(FECHA_INI_obj.strftime('/data/oper/wrfutilV3.0/RUNs/deterministico/WRF/00/WRFOUT/wrfout_d01*'))[0]
    freq = int(os.environ['INTERVALO_WRF']) #minutos
elif MODELO == 'GFS':
    FILEIN = FECHA_INI_obj.strftime('/data/GFS/%Y%m%d_%H/DET/gfs.t%Hz.pgrb2.0p25.f006')
    freq = int(os.environ['INTERVALO_GFS_DET']) #minutos
elif MODELO == 'GEFS':
    FILEIN = FECHA_INI_obj.strftime('/data/GFS/%Y%m%d_%H/ENS/01/gep01.t%Hz.pgrb2.0p50.f006')
    freq = int(os.environ['INTERVALO_GFS_ENS']) #minutos
else:
    print("Modelos {} desconocido !!".format(MODELO))
    exit(10)


#Hacemos merge de las tablas para tener el CIPI, Tipo y ID de los puntos de interes que estan
#en el dominio del modelo 
OUTDIR_ij = BASEDIR + '/PostDF/'
Puntos_interes = af.merge_Puntos_Interes(DIR_TEMPLATES, OUTDIR_ij, FILEIN, MODELO, estimo_ij = True)
os.remove(OUTDIR_ij + 'CIPI2ij.csv')
#Puntos_interes = af.merge_Puntos_Interes(DIR_TEMPLATES)

NPlazo = np.array(np.arange(PLAZOS/(freq/60) + 1) * (freq/60), dtype = np.int64) #Sin el +1 no se tiene en cuenta la ultima hora de pronostico

if MODELO == 'GFS':
    ind = np.argwhere(np.logical_and(NPlazo > 120, np.mod(NPlazo, 3) != 0))
    NPlazo = np.delete(NPlazo, ind)
elif MODELO == 'GEFS':
    ind = np.argwhere(np.logical_and(NPlazo > 192, np.mod(NPlazo, 6) != 0))
    NPlazo = np.delete(NPlazo, ind)


cols = Puntos_interes[['CIPI', 'tipo', 'ID']]
cols.columns = ['CIPI', 'Tipo', 'ID']
cols.set_index('CIPI', inplace = True)

CIPI = np.sort(Puntos_interes['CIPI'].unique())

indices = ['Ciclo','Plazo','CIPI']

miindex = pd.MultiIndex.from_product([[CICLO], NPlazo, CIPI])

#Inicio b
b = pd.DataFrame(np.full([len(miindex), 2], 0.), index = miindex, columns = ['b0', 'b1'])
b.rename_axis(indices, inplace = True)
b = b.join(cols, on = 'CIPI')
b.reset_index(inplace = True) 
b = b.set_index(['Ciclo', 'Plazo', 'CIPI', 'Tipo', 'ID'])

#Inicio P
P = pd.DataFrame(np.full([len(miindex), 4], 0.), index = miindex, columns = ['a', 'b', 'c', 'd'])
P.rename_axis(indices, inplace = True)
P.loc[(slice(None)), ['a', 'd']] = [4., 4.]
P = P.join(cols, on = 'CIPI')
P.reset_index(inplace = True) 
P = P.set_index(['Ciclo', 'Plazo', 'CIPI', 'Tipo', 'ID'])

#Inicio R
R = pd.DataFrame(np.full([len(miindex), 1], 4.), index = miindex, columns = ['R'])
R.rename_axis(indices, inplace = True)
R = R.join(cols, on = 'CIPI')
R.reset_index(inplace = True) 
R = R.set_index(['Ciclo', 'Plazo', 'CIPI', 'Tipo', 'ID'])

#Inicio Q
Q = pd.DataFrame(np.full([len(miindex), 2], 0.01), index = miindex, columns = ['b0', 'b1'])
Q.rename_axis(indices, inplace = True)
Q = Q.join(cols, on = 'CIPI')
Q.reset_index(inplace = True) 
Q = Q.set_index(['Ciclo', 'Plazo', 'CIPI', 'Tipo', 'ID'])

#Inicio b_hist
indices = ['Fecha', 'Ciclo', 'Plazo', 'CIPI']
miindex = pd.MultiIndex.from_product([[ayer], [CICLO], NPlazo, CIPI])
b_hist = pd.DataFrame(np.full([len(miindex), 2], 0.), index = miindex, columns = ['b0', 'b1'])
b_hist.rename_axis(indices, inplace = True)
b_hist = b_hist.join(cols, on = 'CIPI')
b_hist.reset_index(inplace = True) 
b_hist = b_hist.set_index(['Fecha', 'Ciclo', 'Plazo', 'CIPI', 'Tipo', 'ID'])


b.to_hdf(OUTDIR + 'b.h5', key = 'b')
P.to_hdf(OUTDIR + 'P.h5', key = 'P')
R.to_hdf(OUTDIR + 'R.h5', key = 'R')
Q.to_hdf(OUTDIR + 'Q.h5', key = 'Q')
b_hist.to_hdf(OUTDIR + 'b_hist.h5', key = 'b')

os.system('mkdir -p {}/iniciales/{}'.format(DIR_COEF, CICLO))
os.system('cp -r {} {}/iniciales/ '.format(OUTDIR, DIR_COEF))

