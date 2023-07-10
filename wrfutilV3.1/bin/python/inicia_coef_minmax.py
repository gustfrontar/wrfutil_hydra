import numpy as np
import pandas as pd
import os
from datetime import datetime, timedelta
import auxiliary_functions_DF as af
import glob

FECHA_INI = os.environ['FECHA_INI']
CICLO = os.environ['CICLO']
BASEDIR = os.environ['BASEDIR']
DIR_COEF = BASEDIR + '/CALIB/minmax/'
MODELO = os.environ['FILEIN_TYPE']
PLAZOS = int(os.environ['PLAZO'])
WRFUTILDIR = os.environ['WRFUTILDIR']
DIR_TEMPLATES = WRFUTILDIR + "/templates/"

#OUTDIR = DIR_COEF + '/iniciales/' + CICLO + '/'
OUTDIR = DIR_COEF + CICLO + '/'
#print(OUTDIR)
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

FECHA_INI_obj = datetime.strptime(FECHA_INI, '%Y/%m/%d')
ayer = FECHA_INI_obj - timedelta(days = 1)

if MODELO == 'WRF':
    FILEIN = glob.glob(FECHA_INI_obj.strftime('/data/oper/wrfutilV3.0/RUNs/deterministico/WRF/00/WRFOUT/wrfout_d01*'))[0]
elif MODELO == 'GFS':
    FILEIN = FECHA_INI_obj.strftime('/data/GFS/%Y%m%d_%H/DET/gfs.t%Hz.pgrb2.0p25.f006')
elif MODELO == 'GEFS':
    FILEIN = FECHA_INI_obj.strftime('/data/GFS/%Y%m%d_%H/ENS/01/gep01.t%Hz.pgrb2.0p50.f006')
else:
    print("Modelos {} desconocido !!".format(MODELO))
    exit(10)

#Hacemos merge de las tablas para tener el CIPI, Tipo y ID de los puntos de interes que estan
#en el dominio del modelo 
OUTDIR_ij = BASEDIR + '/PostDF/'
Puntos_interes = af.merge_Puntos_Interes(DIR_TEMPLATES, OUTDIR_ij, FILEIN, MODELO, estimo_ij = True)
os.remove(OUTDIR_ij + 'CIPI2ij.csv')


resto = PLAZOS%24 + int(CICLO)

ndias = PLAZOS//24 + + resto//24

NPlazo = [x*24 for x in range(ndias + 1)]

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
b.to_hdf(OUTDIR + 'b_Tmin.h5', key = 'b')
b.to_hdf(OUTDIR + 'b_Tmax.h5', key = 'b')

#Inicio P
P = pd.DataFrame(np.full([len(miindex), 4], 0.), index = miindex, columns = ['a', 'b', 'c', 'd'])
P.rename_axis(indices, inplace = True)
P.loc[(slice(None)), ['a', 'd']] = [4., 4.]
P = P.join(cols, on = 'CIPI')
P.reset_index(inplace = True)
P = P.set_index(['Ciclo', 'Plazo', 'CIPI', 'Tipo', 'ID'])
P.to_hdf(OUTDIR + 'P_Tmax.h5', key = 'P')
P.to_hdf(OUTDIR + 'P_Tmin.h5', key = 'P')

#Inicio R
R = pd.DataFrame(np.full([len(miindex), 1], 4.), index = miindex, columns = ['R'])
R.rename_axis(indices, inplace = True)
R = R.join(cols, on = 'CIPI')
R.reset_index(inplace = True)
R = R.set_index(['Ciclo', 'Plazo', 'CIPI', 'Tipo', 'ID'])
R.to_hdf(OUTDIR + 'R_Tmin.h5', key = 'R')
R.to_hdf(OUTDIR + 'R_Tmax.h5', key = 'R')

#Inicio Q
Q = pd.DataFrame(np.full([len(miindex), 2], 0.01), index = miindex, columns = ['b0', 'b1'])
Q.rename_axis(indices, inplace = True)
Q = Q.join(cols, on = 'CIPI')
Q.reset_index(inplace = True)
Q = Q.set_index(['Ciclo', 'Plazo', 'CIPI', 'Tipo', 'ID'])
Q.to_hdf(OUTDIR + 'Q_Tmin.h5', key = 'Q')
Q.to_hdf(OUTDIR + 'Q_Tmax.h5', key = 'Q')

#Inicio b_hist
indices = ['Fecha', 'Ciclo', 'Plazo', 'CIPI']
miindex = pd.MultiIndex.from_product([[ayer], [CICLO], NPlazo, CIPI])
b_hist = pd.DataFrame(np.full([len(miindex), 2], 0.), index = miindex, columns = ['b0', 'b1'])
b_hist.rename_axis(indices, inplace = True)
b_hist = b_hist.join(cols, on = 'CIPI')
b_hist.reset_index(inplace = True)
b_hist = b_hist.set_index(['Fecha', 'Ciclo', 'Plazo', 'CIPI', 'Tipo', 'ID'])
b_hist.to_hdf(OUTDIR + 'b_hist_Tmin.h5', key = 'b')
b_hist.to_hdf(OUTDIR + 'b_hist_Tmax.h5', key = 'b')


os.system('mkdir -p {}/iniciales/{}'.format(DIR_COEF, CICLO))
os.system('cp -r {} {}/iniciales/'.format(OUTDIR, DIR_COEF, CICLO))

