import pandas as pd
from datetime import datetime
import numpy as np
import os
import sys
WRFUTILDIR = os.environ['WRFUTILDIR']
sys.path.append(WRFUTILDIR + '/bin/python/')
import auxiliary_functions_DF as af

PATHIN = '/data/calibracion/wrfutilV3.1/RUNs/HIST/POST/' #Path con los DF con el viento calibrado
PATHOUT = WRFUTILDIR+'/RUNs/HIST/POST/' #Path con los df sin viento calibrado
DIR_TEMPLATES=WRFUTILDIR+'/templates/'


FECHA_INI = datetime(2019, 11, 1, 0)
FECHA_FIN = datetime(2021, 5, 31, 0)

calibrados = {'deterministico':'WRF',
              'ensamble':'WRF',
              'CalibGFS':'GFS',
              'CalibGEFS':'GEFS'}

fechas = pd.date_range(FECHA_INI, FECHA_FIN, freq = '24H')

for NOMBRE in calibrados.keys():

    modelo = calibrados[NOMBRE]

    for dia in fechas:

        print(dia)
        #Path a los dataframes con viento calibrado y sin calibrar
        calibrated_filename = dia.strftime(PATHIN + '%Y%m%d_%H0000/{}/%Y%m%d_%H0000_horario.h5'.format(NOMBRE))
        uncalibrated_filename = dia.strftime(PATHOUT + '%Y%m%d_%H0000/{}/%Y%m%d_%H0000_horario.h5'.format(NOMBRE))

        #Si no existe el archivo con los datos sin calibrar paso a la fecha siguiente
        if not os.path.isfile(uncalibrated_filename):
            continue

        #Leo el DF sin la calibracion del viento completo
        DF_uncalibrated = pd.read_hdf(uncalibrated_filename, key = 'DF' + modelo.lower())

        if os.path.isfile(calibrated_filename):
            #Del DF calibrado leo solo el viento calibrado
            DF_calibrated = pd.read_hdf(calibrated_filename, key = 'DF' + modelo.lower())[['Umet10cal', 'Vmet10cal', 'magVientoObs_cal']]

            #Agrego las variables de viento calibrado
            DF_uncalibrated['Umet10cal'] = np.nan
            DF_uncalibrated['Vmet10cal'] = np.nan
            DF_uncalibrated['magVientoObs_cal'] = np.nan

            #Le hago al DF sin las calibraciones un update del que las tiene para completarle los valores
            DF_uncalibrated.update(DF_calibrated)
        else:
            print('No hay DF con calibracion del viento. Completo con NaNs')
            #Por continuidad completa con nan
            DF_uncalibrated['Umet10cal'] = np.nan
            DF_uncalibrated['Vmet10cal'] = np.nan

            #Leo los puntos de interes y me quedo solo con las variables CIPI, Tipo y ID
            Pinteres = af.merge_Puntos_Interes(DIR_TEMPLATES, dia.strftime(PATHIN + '%Y%m%d_%H0000/{}/'.format(NOMBRE)))
            cols = Pinteres[['CIPI', 'tipo', 'ID']]
            cols.columns = ['CIPI', 'Tipo', 'ID']
            cols.set_index('CIPI', inplace = True)

            #Hago el promedio de la observaciones en cada CIPI
            tmp = DF_uncalibrated['magVientoObs'].mean(level = ['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Miembro'], skipna = True)

            #Hago un join en CIPI para que aparezcan el tipo y el ID.
            tmp = tmp.to_frame().join(cols, on = 'CIPI')

            #Pongo iguales indices que en el DF
            tmp.reset_index(inplace = True)
            tmp.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Tipo', 'ID', 'Miembro'], inplace = True)

            #Agrego la variable obscal
            DF_uncalibrated['magVientoObs_cal'] = tmp

        #Guardo
        DF_uncalibrated.to_hdf(uncalibrated_filename, key = 'DF' + modelo.lower(), mode = 'w')



