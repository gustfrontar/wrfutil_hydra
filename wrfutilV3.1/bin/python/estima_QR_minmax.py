#!/usr/bin/env python3

import numpy as np
import pandas as pd
from datetime import datetime, date, timedelta
import os
import sys
from copy import deepcopy
import auxiliary_functions_DF as af

def main():

    #Definidas en experimento.conf
    FECHA_INI = os.environ['FECHA_INI']  
    CICLO = os.environ['CICLO'] 
    NOMBRE = os.environ['NOMBRE']
    miembro_ini = int(os.environ['MIEMBRO_INI'])
    miembro_fin = int(os.environ['MIEMBRO_FIN'])
    DIR_COEF = os.environ['BASEDIR'] + '/CALIB/minmax/' 
    OUTBASEDIR = os.environ['OUTDIR_CALIB']

    WRFUTILDIR = os.environ['WRFUTILDIR']
    DIR_TEMPLATES = WRFUTILDIR + "/templates/"


    PLAZO = os.environ['PLAZO']

    modelo = os.environ['FILEIN_TYPE'] 

    tipos_obs = eval(os.environ['TIPOS_OBS'])

    nmiembros = (miembro_fin - miembro_ini) + 1

    #Cantidad de dias previos sobre los que se estiman Q y R
    ndias = 7
    
    #Cantidad de dias para los que busco los datos
    ndays = ndias + int(PLAZO)//24

    date_today = datetime.strptime(FECHA_INI,'%Y/%m/%d') # (Dia actual)  

    date_ndias = date_today - timedelta(days = ndias)

    date_ndays = date_today - timedelta(days = ndays) 

    # Abro los DF de los dias previos de donde voy a sacar el error observado calibrar
    n = True
    fechas = pd.date_range(date_ndays, date_today)
    for ite_date in fechas: 
        filename_DF = OUTBASEDIR + ite_date.strftime("%Y%m%d_{0}0000/{1}/%Y%m%d_{0}0000_minmax.h5".format(CICLO, NOMBRE))
        if os.path.isfile(filename_DF) == False: 
            continue 
        tmp = pd.read_hdf(filename_DF, key = 'DF' + modelo.lower()) 
        if n: 
            df_t2m_prev = tmp 
            n = False
        else: 
            df_t2m_prev = pd.concat([df_t2m_prev,tmp], sort = True) 
    if n:
        sys.exit()

    df_t2m_prev = df_t2m_prev.query('Tipo == @tipos_obs')

    df_t2m_prev = df_t2m_prev.sort_index(level = 0)

    plazos = df_t2m_prev.index.get_level_values('Plazo').unique()  
    variables = ['Tmin', 'Tmax']

    #os.remove(DIR_COEF + CICLO + '/Q.h5')

    OUTDIR = OUTBASEDIR + date_today.strftime('%Y%m%d') + '_' + CICLO + '0000/' + NOMBRE + '/'
#    Puntos_interes = af.merge_Puntos_Interes(DIR_TEMPLATES, OUTDIR)
#    Puntos_interes = Puntos_interes[['CIPI', 'tipo', 'ID']].values

    #Calculo la media del ensamble, si es deterministico no hace nada
    df_t2m_prev = df_t2m_prev.mean(level = ['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI'], skipna = True)

    lista_CIPI = df_t2m_prev.index.get_level_values('CIPI').unique()


    for var in variables:
    
        DFhist_R = pd.read_hdf(DIR_COEF + CICLO + '/R_{}.h5'.format(var))
        DFhist_b = pd.read_hdf(DIR_COEF + CICLO + '/b_hist_{}.h5'.format(var))
        Q = pd.DataFrame(data = 0.01, columns = ['b0', 'b1'], index = DFhist_R.index, dtype = np.dtype('float'))


        b_last = pd.read_hdf(DIR_COEF + CICLO + '/b_{}.h5'.format(var))
        b_last = pd.concat([b_last], keys=[date_today + timedelta(days = -1)], names=['Fecha'])
        DFhist_b = DFhist_b.append(b_last)
        DFhist_b = DFhist_b.sort_index(level = 0)


        DFhist_b_ind = deepcopy(DFhist_b)
        DFhist_b_ind.index = DFhist_b.index.set_levels(np.arange(1, ndias+2), level = 'Fecha')
        #El +2 en np.arange viene de que pongo un dia mas para tener en cuenta el dato b_last que
        #concateno y otro mas porque np.arange no da el valor final

        b_last = b_last.sort_index(level = 0)

        for ip in plazos: # Loop para los plazos

            ini = (date_today - timedelta(days = ndias))
            fin = (date_today - timedelta(days = 1))

            data_obs = deepcopy(df_t2m_prev.loc[(slice(None), CICLO, ip, slice(ini, fin), slice(None)), :])

            if len(data_obs) == 0:
                continue

            valideces = data_obs.index.get_level_values('Validez').unique()
            val = [x.date() for x in valideces]

            #Tomamos solo estos puntos de interes porque puede romper cuando se agregan nuevas
            PI_CIPI = data_obs.loc[(slice(None), CICLO, ip, valideces[0], slice(None)), :].index.get_level_values('CIPI')

            fcst = data_obs.loc[(slice(None), CICLO, ip, slice(ini, fin), PI_CIPI), var + '_f']
    
            #Repito los ultimos coeficientes tantas veces como valideces haya
            b_ip = pd.concat([b_last.loc[(slice(None), CICLO, ip, PI_CIPI), :]]*len(val))

            #Le asigo a b_ip los indices de T2
            b_ip.index = fcst.index


            yp = (fcst * b_ip['b1'] + b_ip['b0'])
            data_obs.loc[(slice(None), CICLO, ip, slice(ini, fin), PI_CIPI), var + '_cal'] = (fcst - yp)

            #R
            T_obs_to2 = data_obs[var + '_o_cal']
            T_fcst_to2 = data_obs[var + '_f']
            T_cal_to2 = data_obs[var + '_cal']
            error_obs_to2 = T_obs_to2- T_fcst_to2
            error_prono_to2 = T_cal_to2 - T_fcst_to2
            delta_error = error_obs_to2 - error_prono_to2
            R = delta_error.groupby(['CIPI']).var()

            #Q
            data_b1 = DFhist_b_ind.loc[(slice(2, ndias+1), CICLO, ip, slice(None)), :].droplevel('Fecha')
            data_b_1 = DFhist_b_ind.loc[(slice(1, ndias), CICLO, ip, slice(None)), :].droplevel('Fecha')
            delta_b = data_b1 - data_b_1
            Q_tmp = delta_b.groupby(['CIPI']).var()

            for CIPI in lista_CIPI:
                if (len(R) != 0) and ((CIPI) in R.index):
                    if np.isnan(R.loc[CIPI]) or R.loc[CIPI] == 0:
                       pass
                    else:
                        DFhist_R.loc[CICLO, ip, CIPI] = R.loc[CIPI]

                if (len(Q_tmp) != 0) and ((CIPI) in Q_tmp.index):
                    if np.isnan(Q_tmp.loc[CIPI]).any():
                        pass
                    else:
                        Q.loc[(CICLO, ip, CIPI), :] = Q_tmp.loc[CIPI]


        b_update = DFhist_b.loc(axis = 0)[date_ndias:]

        DFhist_R.to_hdf(DIR_COEF + CICLO + '/R_{}.h5'.format(var), key = 'R', mode = 'w')
        Q.to_hdf(DIR_COEF + CICLO + '/Q_{}.h5'.format(var), key = 'Q', mode = 'w')
        b_update.to_hdf(DIR_COEF + CICLO + '/b_hist_{}.h5'.format(var), key = 'b', mode = 'w')


###################################################################################################  
if __name__ == "__main__":

    print ("<<<< Comienza el proceso de estimacion de Q y R ")

    main()

    print (">>>> Finaliza el proceso de estimacion de Q y R ")

