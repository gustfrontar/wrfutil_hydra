#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ===========================================================================
# Paso 0: LIMPIO ENTORNO Y CARGO LIBRERIAS
# ===========================================================================

# Cargo librerias
import numpy as np
import os
import pandas as pd
from datetime import datetime, timedelta, date
import sys
from multiprocessing import Pool
import auxiliary_functions_DF as af
import glob

def main():

    #Definidas en experimento.conf
    FECHA_INI = os.environ['FECHA_INI'] 
    CICLO = os.environ['CICLO'] 
    NOMBRE = os.environ['NOMBRE']
    WRFUTILDIR = os.environ['WRFUTILDIR']
    miembro_ini = int(os.environ['MIEMBRO_INI'])
    miembro_fin = int(os.environ['MIEMBRO_FIN']) 

    #Definidas en experimento.plgn
    FECHA_CALIB = os.environ['FECHA_DFOBS'] 
    CICLO_CALIB = os.environ['CICLO_DFOBS'] 
    VARIABLES = eval(os.environ['VAR_DF'])  
    INDIR_OBS = os.environ['INDIR_OBS']
    OUTBASEDIR = os.environ['OUTDIR_DF'] 

    #Definidas en config.env
    ICORE = int(os.environ['ICORE']) 
    BASEDIR = os.environ['BASEDIR']


    sys.path.append(WRFUTILDIR + "/templates/")
    import gribNames


    modelo = os.environ['FILEIN_TYPE'] 
    #DIR_PUNTOS = BASEDIR + "/PostDF/" 
    DIR_TEMPLATES = WRFUTILDIR + "/templates/"



    nmiembros = (miembro_fin - miembro_ini) + 1 #Cantidad de miembros
    mems = [x for x in range(miembro_ini, miembro_fin + 1)]


    FECHA_INI_obj = datetime.strptime(FECHA_INI + CICLO, '%Y/%m/%d%H')
    FECHA_INI = FECHA_INI_obj.strftime('%Y%m%d')
    ANIO = FECHA_INI_obj.year
    MES = FECHA_INI_obj.month
    DIA = FECHA_INI_obj.day

    PLAZOS = int(os.environ['PLAZO'])
    INPDIR = os.environ['DIRIN']

    if modelo == 'WRF':
        #FILEIN = BASEDIR + FECHA_INI_obj.strftime('/WRF/{0:02d}/WRFOUT/wrfout_d01_%Y-%m-%d_%H_%M_00'.format(mems[0])) 
        FILEIN = glob.glob(BASEDIR + '/WRF/*/WRFOUT/wrfout_d01_*')[0]
        freq = int(os.environ['INTERVALO_WRF']) #minutos
    elif modelo == 'GFS':
        #FILEIN = FECHA_INI_obj.strftime(BASEDIR + '/GFS/%Y%m%d_%H/DET/gfs.t%Hz.pgrb2.0p25.f006')
        FILEIN = glob.glob(FECHA_INI_obj.strftime(BASEDIR + '/HIST/GFS_CAL/%Y%m%d_%H/DET/gfs.t%Hz.pgrb2.0p25.f*'))[0]
        freq = int(os.environ['INTERVALO_GFS_DET']) #minutos
        var_modelo = gribNames.GFS   #diccionario con nombre de variable del WRF a nombre en el GRIB de GFS
    elif modelo == 'GEFS':
        #FILEIN = FECHA_INI_obj.strftime(BASEDIR + '/GFS/%Y%m%d_%H/ENS/{0:02d}/gep{0:02d}.t%Hz.pgrb2.0p50.f006'.format(mems[0]))
        FILEIN = glob.glob(FECHA_INI_obj.strftime(BASEDIR + '/HIST/GFS_CAL/%Y%m%d_%H/ENS/*/gep*.t%Hz.pgrb2.0p50.f*'))[0]
        freq = int(os.environ['INTERVALO_GFS_ENS']) #minutos
        var_modelo = gribNames.GEFS   #diccionario con nombre de variable del WRF a nombre en el GRIB de GEFS
    else: 
        print("Modelos {} desconocido !!".format(modelo))
        exit(10)



    OUTDIR = OUTBASEDIR + FECHA_INI_obj.strftime('%Y%m%d_%H') + '0000/'+NOMBRE+'/' 
   
    pool = Pool(processes = min(ICORE, nmiembros))  

    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)

    # ===========================================================================
    # PASO 1: APERTURA DEL ARCHIVO CON LAS CIUDADES/ESTACIONES
    # ===========================================================================

    Puntos_interes = af.merge_Puntos_Interes(DIR_TEMPLATES, OUTDIR, FILEIN, modelo, estimo_ij = True)

    cols = Puntos_interes[['CIPI', 'tipo', 'ID']]
    cols.columns = ['CIPI', 'Tipo', 'ID']
    cols.set_index('CIPI', inplace = True)

    CIPI = np.sort(Puntos_interes['CIPI'].unique())

    # Definimos indices
    DateIndex = pd.date_range(date(ANIO, MES, DIA), periods = 1,freq = 'D')
    NPlazo = np.array(np.arange(PLAZOS/(freq/60) + 1) * (freq/60), dtype = np.int64) #Sin el +1 no se tiene en cuenta la ultima hora de pronostico

    if modelo == 'GFS':
        ind = np.argwhere(np.logical_and(NPlazo > 120, np.mod(NPlazo, 3) != 0))
        NPlazo = np.delete(NPlazo, ind)
    elif modelo == 'GEFS':
        ind = np.argwhere(np.logical_and(NPlazo > 192, np.mod(NPlazo, 6) != 0))
        NPlazo = np.delete(NPlazo, ind)

    NVariables=len(VARIABLES) # cant de variables a guardar en el dataframe

    # Creamos estructura de indices 
    miindex = pd.MultiIndex.from_product([DateIndex, [CICLO], NPlazo, CIPI, mems])

    # Creamos dataFrame vacio 
    DF_fcst = pd.DataFrame(np.full((len(miindex), NVariables), np.nan), index = miindex, columns = VARIABLES)

    # Armamos en columnas los indices que no podemos generar con el multiindex
    # Validez
    DF_fcst['Validez'] = DF_fcst.index[0][0]  # ==> agrega una columna llena de datos iguales al 1er indice
    validez_ini = FECHA_INI_obj
    for ip in NPlazo.tolist():
        DF_fcst.loc[(slice(None), slice(None), ip, slice(None), slice(None)), 'Validez'] = validez_ini + timedelta(hours = ip)

    # Agregamos nombres a los indices del DF que no tenian:
    nombres = ['Fecha','Ciclo','Plazo','CIPI','Miembro']
    DF_fcst = DF_fcst.rename_axis(nombres)

    # Reorganizamos los indices
    DF_fcst.reset_index(inplace = True) # reseteo los indices y los pone como columnas de datos
    DF_fcst = DF_fcst.set_index(['Fecha','Ciclo','Plazo','Validez','CIPI','Miembro']) # seteo los indices con este nuevo orden

    # ===========================================================================
    # Paso 1: LECTURA DE LOS PRONOSTICOS
    # ===========================================================================

    #Genero una lista de Dataframes en la que cada elemento contiene 1 miembro
    DFlist = [DF_fcst.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), m), :] for m in mems]

    #Lanzo la lectura de los POST de cada miembro en paralelo 
    if modelo == 'WRF':
        pool_out = pool.starmap(af.WRF2DF, [(DF, INPDIR, FECHA_INI_obj, NPlazo, VARIABLES, Puntos_interes) for DF in DFlist])
    elif modelo == 'GFS' or  modelo == 'GEFS' :
        pool_out = pool.starmap(af.GFS2DF, [(DF, INPDIR, FECHA_INI_obj, NPlazo, var_modelo, Puntos_interes) for DF in DFlist])

    #Itero sobre la lista que devuelve la paralelizacion para completar el DFgfs con los datos de cada miembro
    for DF_mem in pool_out:
        m = DF_mem.index.get_level_values('Miembro').unique()[0]
        DF_fcst.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), m), :] = DF_mem


    #Agrego las columnas de tipo y ID
    DF_fcst = DF_fcst.join(cols, on = 'CIPI')

    #Agrego las columnas de tipo y ID como indices
    DF_fcst.reset_index(inplace = True) # reseteo los indices y los pone como columnas de datos
    DF_fcst = DF_fcst.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Tipo', 'ID', 'Miembro']) # seteo los indices con este nuevo orden

    DF_fcst.to_hdf(OUTDIR + '/' + FECHA_INI + '_' + CICLO + '0000_horario.h5', key = 'DF' + modelo.lower(), mode = 'w')

    pool.close()


    # ===========================================================================
    # PASO 3: LECTURA DE LAS OBSERVACIONES
    # ===========================================================================

    #Fecha de la ultima corrida en que se llenaron con observaciones los dataframes
    FECHA_INI_OBS = datetime.strptime(FECHA_CALIB+CICLO_CALIB,'%Y/%m/%d%H') 

    #HS_OBS= FECHA_INI_OBS -  datetime.datetime.strptime(FECHA_INI+CICLO,'%Y/%m/%d%H')
    HS_OBS= FECHA_INI_obj -  FECHA_INI_OBS
    HS_OBS=int(HS_OBS.total_seconds() //3600) + int(PLAZOS)
    DFobstot=[]
    a = True
    for ip in range(HS_OBS + 1):
        f = FECHA_INI_OBS + timedelta(hours = ip)   #datetime 
        f_str = datetime.strftime(f, '%Y%m%d')
        ANIO, MES, DIA, HORA = f.year, f.month, f.day, f.hour # enteros    
        filename_obs = INDIR_OBS + 'asm_synop_' + f.strftime('%Y%m%d%H00') + '.lst'
        try:
            tmp = af.open_obs(filename_obs, f_str, HORA) # devuelve un DF con 8 cols: estacion + 7 variables
        except:
            print("no se encontro el archivo {} .. continuando".format(filename_obs))
            continue

        tmp['Validez'] = f # agrego una columna "Validez" y calculo a partir de fecha y plazo
        print(ip, filename_obs)
        if a:
           DFobstot = tmp
           a = False
        else:
           DFobstot = pd.concat([DFobstot,tmp], axis = 0)   

    if len(DFobstot) == 0:
        print("No se encontraron observaciones. Saliendo!") #TODO aca se podria hacer para que escriba un archivo
        sys.exit()

    #Agregamos esta columna para identificar la fuente de informacion y facilitar el join/merge
    DFobstot['Tipo'] = 'ESTACION'

    DFobstot = DFobstot.astype({'ID': str})

    # definimos los indices y los ponemos en orden: 
    nombres = ['Validez', 'ID', 'Tipo']
    DFobstot = DFobstot.set_index(nombres)

    valideces_obs = DFobstot.index.get_level_values('Validez').unique()

    # ===========================================================================
    # PASO 4: APERTURA DE LOS DF PARA COMPLETARLOS CON LAS OBSERVACIONES
    # ===========================================================================
    print("abrimos los DF y vamos incluyendo las observaciones")
    fin = FECHA_INI_obj
    ini = fin - timedelta(hours = max(PLAZOS, HS_OBS))

    fechas = [ini + timedelta(hours = x) for x in range(int((fin - ini).total_seconds()//3600) + 1)]

    for iter in fechas:
        #print(iter)
        ih = datetime.strftime(iter, '%Y%m%d_%H')
        OUTOBSDIR = OUTBASEDIR + datetime.strftime(iter, '%Y%m%d_%H') + '0000/' + NOMBRE + '/'
        INPFILE = ih + '0000_horario.h5'
        try :
            DF_fcst = pd.read_hdf(OUTOBSDIR + INPFILE, key = 'DF' + modelo.lower())
        except Exception as e:
            #print(e)
            continue

        valideces_fcst = DF_fcst.index.get_level_values('Validez').unique()


        if len(set(valideces_obs).intersection(valideces_fcst)) == 0:
            print('Intentamos completar observaciones en fechas que no están en el DF') #TODO aca se podria que escriba un archivo

        #Agrega las observaciones en la columna correspondiente
        DF_fcst.reset_index(inplace = True)
        DF_fcst.set_index(nombres, inplace = True)
        DF_fcst.update(DFobstot)
        DF_fcst.reset_index(inplace = True)
        DF_fcst.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Tipo', 'ID', 'Miembro'], inplace = True)


        #Agregamos el promedio de los valores observados que va a ser la observacion que se usa
        #para calibrar en caso de que haya mas de un valor observado de temperatura o viento
        DF_mean = DF_fcst[['T2obs', 'magVientoObs']].mean(level = ['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI'], skipna = True)
        DF_mean.columns = ['T2obs_cal', 'magVientoObs_cal']
        DF_fcst.reset_index(inplace = True)
        DF_fcst.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI'], inplace = True)
        DF_fcst.update(DF_mean)
        DF_fcst.reset_index(inplace = True)
        DF_fcst.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Tipo', 'ID', 'Miembro'], inplace = True)

        # Guardamos el DF
        print(ih, INPFILE)
        DF_fcst.to_hdf(OUTOBSDIR + INPFILE, key = 'DF' + modelo.lower(), mode = 'w') # lleva el nombre del DF




if __name__ == "__main__":

    main()

