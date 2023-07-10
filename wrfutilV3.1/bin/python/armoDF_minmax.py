import pandas as pd
import numpy as np
import sys
from datetime import datetime, timedelta
import os
import auxiliary_functions_DF as af
import glob


def main():

    #Inicializacion de las variables desde el entorno
    FECHA_INI = os.environ['FECHA_INI'] #YYYY/MM/DD 
    CICLO = os.environ['CICLO'] #HH
    FECHA_DFOBS = os.environ['FECHA_DFOBS'] #YYYY/MM/DD
    CICLO_DFOBS = os.environ['CICLO_DFOBS'] #HH
    NOMBRE = os.environ['NOMBRE']
    WRFUTILDIR = os.environ['WRFUTILDIR']
    miembro_ini = int(os.environ['MIEMBRO_INI'])
    miembro_fin = int(os.environ['MIEMBRO_FIN'])
    PLAZOS = int(os.environ['PLAZO'])
    BASEDIR = os.environ['BASEDIR']
    INPDIR_OBS = os.environ['INDIR_OBS']
    OUTBASEDIR = os.environ['OUTDIR_DF']
    modelo = os.environ['FILEIN_TYPE']

    DIR_PUNTOS = BASEDIR + "/PostDF/"
    DIR_TEMPLATES = WRFUTILDIR + "/templates/"
    VARIABLES = eval(os.environ['VAR_DF'])  ## ESTO HAY QUE CAMBIARLO EN EL EXPERIMENTO

    nmiembros = (miembro_fin - miembro_ini) + 1 #Cantidad de miembros
    mems = [x for x in range(miembro_ini, miembro_fin + 1)] #Lista con los miembros

    resto = PLAZOS%24 + int(CICLO)

    PLAZOS_dias = PLAZOS//24 + resto//24#[dias]

    FECHA_obj = datetime.strptime(FECHA_INI, '%Y/%m/%d')
    FECHA_INI_obj = datetime.strptime(FECHA_INI+CICLO, '%Y/%m/%d%H')

    OUTDIR = FECHA_INI_obj.strftime(OUTBASEDIR + '%Y%m%d_%H0000/{}/'.format(NOMBRE))

    #DF con los datos horarios
    filename_DF = FECHA_INI_obj.strftime(OUTDIR + '%Y%m%d_%H0000_horario.h5'.format(NOMBRE))

    #Chequeo de la existencia del archivo
    if not os.path.isfile(filename_DF):
        print('El archivo {} no existe'.format(filename_DF))
        sys.exit(1)

    #Lectura del archivo
    #DF = pd.read_hdf(filename_DF)
    DF = pd.read_hdf(filename_DF, key = 'DF' + modelo.lower())

#    if not 'Miembro' in DF.index.names:
#        DF['Miembro'] = 0
#        DF.reset_index(inplace = True)
#        DF.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'Estacion', 'Miembro'], inplace = True)

    CIPI = np.sort(DF.index.get_level_values('CIPI').unique()) #Lista con los puntos de interes
    valideces = [FECHA_obj + timedelta(days=x) for x in range(PLAZOS_dias + 1)] #Plazos

    #Inicio el dataframe de maximas y minimas
    miindex = pd.MultiIndex.from_product([[FECHA_obj], [CICLO], valideces, CIPI, mems])
    DF_minmax = pd.DataFrame(data = None, columns = VARIABLES, index = miindex, dtype = np.dtype('float'))

    DF_minmax['Plazo'] = DF_minmax.index[0][0]  # ==> agrega una columna llena de datos iguales al 1er indice

    #Calculo los plazos con la validez
    for val in valideces:
            DF_minmax.loc[(slice(None), slice(None), val, slice(None), slice(None)), 'Plazo'] = int((val - FECHA_obj).total_seconds()/3600)

    #Se definen nombres a las columnas y se reasignan los indices para incluir el plazo
    nombres = ['Fecha', 'Ciclo', 'Validez', 'CIPI', 'Miembro']
    DF_minmax = DF_minmax.rename_axis(nombres)
    DF_minmax.reset_index(inplace = True)
    #DF_minmax.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Miembro'], inplace = True)
    DF_minmax.set_index(['Validez', 'CIPI', 'Miembro'], inplace = True)

    #Itero sobre la validez para calcular las minimas y maximas
    for val in valideces:
    
        #Funcion que define si se completan o no las maximas
        minima, maxima = af.condicion_completado(val, valideces, CICLO, resto%24)

        #Indice del DF de minima y maxima
        tupla_minmax = (FECHA_obj, CICLO, slice(None), val, slice(None), slice(None))

        #Cuando se usen los DFs deterministicos del WRF que tengan el indice Miembro borrar este IF
        if minima:
            #Indices para calcular las minimas y maximas desde el DF horario
            tupla_min = (slice(None), CICLO, slice(None),
                         slice(val, val + timedelta(hours = 12), None), slice(None),
                         slice(None), slice(None), slice(None))
            tmp_min = DF.loc[tupla_min, 'T2'].groupby(['CIPI', 'Miembro']).min()
            tmp_min = pd.concat([tmp_min], keys=[val], names=['Validez'])
            tmp_min.name = 'Tmin_f'
            #DF_minmax.loc[tupla_minmax, 'Tmin_f'] = DF.loc[tupla_min, 'T2'].groupby(['CIPI', 'Miembro']).min().values #Minima del dia
            DF_minmax.update(tmp_min)
        
        if maxima:
            tupla_max = (slice(None), CICLO, slice(None), 
                         slice(val + timedelta(hours = 12), val + timedelta(hours = 24), None),
                         slice(None), slice(None), slice(None), slice(None))

            tmp_max = DF.loc[tupla_max, 'T2'].groupby(['CIPI', 'Miembro']).max()
            tmp_max = pd.concat([tmp_max], keys=[val], names=['Validez'])
            tmp_max.name = 'Tmax_f'
            DF_minmax.update(tmp_max)

    DF_minmax.reset_index(inplace = True)
    DF_minmax.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Miembro'], inplace = True)


    Puntos_interes = af.merge_Puntos_Interes(DIR_TEMPLATES, OUTDIR)

    cols = Puntos_interes[['CIPI', 'tipo', 'ID']]
    cols.columns = ['CIPI', 'Tipo', 'ID']
    cols.set_index('CIPI', inplace = True)

    DF_minmax = DF_minmax.join(cols, on = 'CIPI')

    #Agrego las columnas de tipo y ID como indices
    DF_minmax.reset_index(inplace = True) # reseteo los indices y los pone como columnas de datos
    DF_minmax = DF_minmax.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Tipo', 'ID', 'Miembro']) # seteo los indices con este nuevo orden


    #Guardo el DF
    DF_minmax.to_hdf(FECHA_INI_obj.strftime(OUTDIR + '%Y%m%d_%H0000_minmax.h5'.format(NOMBRE)), key = 'DF' + modelo.lower(), mode = 'w')



    ####################################################################################################
    # Agrego las observaciones
    ####################################################################################################

    #Columnas del archivo
    #names_diario = ['Fecha', 'Estacion', 'Latitud', 'Longitud', 'Maxima 00', 'Maxima 12',
    #                'Maxima absoluta', 'Minima 00', 'Minima 12', 'Minima absoluta', 'PP', 'Algo']
    names_diario = ['Validez', 'ID', 'Latitud', 'Longitud', 'Tmax_o', 'Maxima 12',
                    'Maxima absoluta', 'Minima 00', 'Tmin_o', 'Minima absoluta', 'PP', 'Algo']

    FECHA_OBS = datetime.strptime(FECHA_DFOBS + CICLO_DFOBS, '%Y/%m/%d%H')
    fechas_completar = pd.date_range(FECHA_OBS, FECHA_INI_obj + timedelta(days = PLAZOS_dias), freq = '6H', closed='right')
    #closed='right' es para que no incluya FECHA_OBS en la lista

    for ite_date in fechas_completar:

        CICLO_comp = ite_date.strftime('%H')
    
        if np.mod(int(CICLO_comp), 12) != 0: #Si los ciclos no son 00 o 12 salgo porque no hay observaciones
            continue

        #Archivo con las observaciones
        filename_diario = ite_date.strftime(INPDIR_OBS + 'diario_%Y%m%d%H00.lst')

        if not os.path.isfile(filename_diario):
            print('El archivo {} no existe'.format(filename_diario))
            continue


        #Lectura del archivo
        diario = pd.read_csv(filename_diario, sep = ',', header = None, names = names_diario,
                             usecols = [0, 1, 4, 8], na_values = ['-99.0', '99.0', '99', '-99'],
                             parse_dates = ['Validez'])

        #Paso las temperaturas de celsius a kelvin
        diario[['Tmax_o', 'Tmin_o']] = diario[['Tmax_o', 'Tmin_o']] + 273.15

        diario['Tipo'] = 'ESTACION'
        diario = diario.astype({'ID': str})

        diario.set_index(['Validez', 'Tipo', 'ID'], inplace = True)

        valideces_obs = diario.index.get_level_values('Validez').unique()

        #Fechas hacia atras sobre las que se itera para completar los dataframes con observaciones
        fechas = [ite_date - timedelta(hours = x) for x in range(PLAZOS+1)] 

        if CICLO_comp == '00':
            column_minmax = 'Tmax_o'
            #column_diario = 'Maxima 00'
        else:
            column_minmax = 'Tmin_o'
            #column_diario = 'Minima 12'

        #FECHA_diario = diario.iloc[0]['Fecha']

        if np.isnan(diario[column_minmax]).all():
            print('Los datos que vamos a leer del archivo {} son todos NaN'.format(filename_diario)) #TODO aca se podria escribir un archivo
            continue

        for iter in fechas:
            INOUTFILE = iter.strftime(OUTBASEDIR + '%Y%m%d_%H0000/{}/%Y%m%d_%H0000_minmax.h5'.format(NOMBRE))

            if not os.path.isfile(INOUTFILE):
                continue
            
            DF_minmax = pd.read_hdf(INOUTFILE, key = 'DF' + modelo.lower())

            valideces_fcst = DF_minmax.index.get_level_values('Validez').unique()

            if len(set(valideces_obs).intersection(valideces_fcst)) == 0:
                print('Intentamos completar observaciones en fechas que no están en el DF') #TODO aca se podria que escriba un archivo


            DF_minmax.reset_index(inplace = True)
            DF_minmax.set_index(['Validez', 'Tipo', 'ID'], inplace = True)
            DF_minmax[column_minmax].update(diario[column_minmax])
            DF_minmax.reset_index(inplace = True)
            DF_minmax.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Tipo', 'ID', 'Miembro'], inplace = True)


            #Agregamos el promedio de los valores observados que va a ser la temperatura que se usa
            #para calibrar en caso de que haya mas de un valor observado de temperatura
            DF_mean = DF_minmax[column_minmax].mean(level = ['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI'], skipna = True)
            DF_mean.name = column_minmax + '_cal'
            DF_mean = DF_mean.to_frame()
            DF_minmax.reset_index(inplace = True)
            DF_minmax.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI'], inplace = True)
            DF_minmax.update(DF_mean)
            DF_minmax.reset_index(inplace = True)
            DF_minmax.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Tipo', 'ID', 'Miembro'], inplace = True)


            #Esta era la forma vieja de completar las obsevaciones. El update parece andar mejor
#            estaciones = DF_minmax.index.get_level_values('Estacion').unique()
#            for iest in estaciones:
#                tupla_est = (slice(None), slice(None), slice(None), FECHA_diario, iest, slice(None))
#                if iest in diario.index:
#                    DF_minmax.loc[tupla_est, column_minmax] = diario.loc[iest, column_diario] + 273.15

            DF_minmax.to_hdf(INOUTFILE, key = 'DF' + modelo.lower(), mode = 'w')


if __name__ == "__main__":

    print('Inicio generacion del DF de maximas y minimas')
    main()
    print('Fin generacion del DF de maximas y minimas')

