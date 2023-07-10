import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import os
import auxiliary_functions_DF as af
from copy import deepcopy
import multiprocessing as mp
import gc

#Para que numpy tire FloatingPointError cuando divide por cero
np.seterr(divide='raise')

def interpolacion_temporal(DF, dict_interpolacion, DFhist_b, freq, fcst_var, cal_var):
    """
    Funcion que interpola temporalmente los coeficientes de la calibracion
    para completar los plazos en que no se calibra para alguna estacion.
    """

    #Dia del DF
    date_today = DF.index.get_level_values('Fecha').unique()

    #Ciclo del DF
    CICLO = DF.index.get_level_values('Ciclo').unique()[0]

    #Plazos del DF
    PLAZOS = DF.index.get_level_values('Plazo').unique()

    #Genero una lista de miembros
    miembros = DF.index.get_level_values('Miembro').unique()

    nmiembros = len(miembros)

    H = np.full([nmiembros, 2], 1.)

    #Cantidad de horas entre plazos consecutivos del pronostico
    delta_fcst = PLAZOS[1:] - PLAZOS[:-1]

    DFhist_b_interp = deepcopy(DFhist_b)

    #for iest in dict_interpolacion.keys():
    for CIPI in dict_interpolacion.keys():

        #De dict_interpolacion genero las horas en que hubo observacion
        horas = np.array([x for x in PLAZOS if x not in dict_interpolacion[(CIPI)]])

        #Cantidad de horas entre plazos consecutivos en que hubo observacion
        delta_cal = np.abs(horas[1:]-horas[:-1])

        #Cantidad de horas esperadas entre plazos consecutivos de pronostico para las horas en que hubo observacion
        delta_fcst_obs = [delta_fcst[ind] for ind in range(delta_fcst.size) if PLAZOS[1:][ind] in horas[1:]]

        #Busco donde el delta_cal es distinto a la frecuencia de pronostico
        index_caso = np.argwhere(delta_cal!=freq/60)

        if index_caso.shape[0] == 1:
            index_caso = [int(index_caso[0,0])]
        else:
            index_caso = np.squeeze(index_caso)

        for ind in index_caso:
            #No interpolamos si el bache de observaciones es mayor a 3 horas o su el bache es 
            #igual al que hay entre plazos consecutivos
            if (delta_cal[ind] > 3 or delta_cal[ind] == delta_fcst_obs[ind]):
                continue
            
            #Tomo los coeficientes en los bordes del bache de observaciones
            b = DFhist_b.loc[CICLO, slice(horas[ind], horas[ind]+delta_cal[ind], delta_cal[ind]), CIPI]
            #Calculo la media
            b_mean = b.mean(axis = 0)
            #Tomo las horas de los bordes del bache de observaciones
            b_ind = b.index.get_level_values('Plazo')
            #Genero las horas del bache de observaciones
            indices = [x for x in range(b_ind[0]+1, b_ind[1])]
            #Calibro con la media de los coeficientes
            for ip in indices:
                date_plazo = date_today + timedelta(hours = int(CICLO) + ip)
                fcst_today = DF.loc[date_today, CICLO, ip, date_plazo, CIPI, :]

                #Si las variables pronosticadas son el viento calculo magnitud y direccion
                if fcst_var == ['Umet10', 'Vmet10']:
                    forecast_today = np.sqrt(fcst_today['Umet10']**2 + fcst_today['Vmet10']**2).values
                    dir_viento = np.arctan2(fcst_today['Umet10'], fcst_today['Vmet10'])
                else:
                    forecast_today = fcst_today[fcst_var].values

                H[:, 1] = forecast_today
                yp = np.dot(H, b_mean)
                fcst_today_corrected = (forecast_today - yp).reshape(1, nmiembros)

                tupla = (date_today, CICLO, ip, date_plazo, CIPI, slice(None))

                #Si la variable pronosticada es el viento descompongo la magnitud calibrada en U y V
                if fcst_var == ['Umet10', 'Vmet10']:
                    fcst_today_corrected[fcst_today_corrected<0] = 0
                    wind_dir = np.array([np.sin(dir_viento), np.cos(dir_viento)])
                    DF.loc[tupla, cal_var] = (fcst_today_corrected*wind_dir).T
                else:
                    DF.loc[tupla, cal_var] = fcst_today_corrected

                DFhist_b_interp.loc[CICLO, ip, CIPI] = b_mean

                dict_interpolacion[(CIPI)].remove(ip)

    return DF, dict_interpolacion, DFhist_b_interp



def interpolacion_espacial(DIR_TEMPLATES, OUTDIR, DF, dict_interpolacion, b, R, fcst_var, cal_var):
    """
    Funcion que interpola espacialmente los coeficientes de la calibracion
    para completar los plazos en que no se calibra para alguna estacion.
    """
    #Leo la tabla de puntos de interes y pongo CIPI como indice
    Puntos_interes = af.merge_Puntos_Interes(DIR_TEMPLATES, OUTDIR)
    Puntos_interes_CIPI = Puntos_interes.drop_duplicates(subset='CIPI')
    Puntos_interes_CIPI.set_index(['CIPI'], inplace = True)

    #Dia del DF
    date_today = DF.index.get_level_values('Fecha').unique()   

    #Ciclo del DF
    CICLO = DF.index.get_level_values('Ciclo').unique()[0]

    #Genero una lista de miembros
    miembros = DF.index.get_level_values('Miembro').unique()    

    nmiembros = len(miembros)

    #Mascara que tiene False donde el coeficiente b es 0 que indica que
    #no tiene observaciones entonces no lo uso para interpolar
    b_mask_hist = ~np.isnan(b.where(b!=0)).any(axis = 1).droplevel('Ciclo')

    #Generamos una mascara con los CIPIs donde no hubo observacion inicialmente toda con True y luego
    #agregamos los False donde corresponde.
    b_mask_tot = pd.DataFrame(np.full(len(b_mask_hist.index), True), index = b_mask_hist.index, columns = ['mask'])  

    #Aniadimos False en la mascara donde no hubo observaciones para que
    #no se usen esos CIPIs y plazos en la interpolacion de alguna estacion
    for CIPI in dict_interpolacion.keys():
        b_mask_tot.loc[(dict_interpolacion[(CIPI)], CIPI), 'mask'] = False

    #Generamos una copia
    b_mask_obs = deepcopy(b_mask_tot)

    #Mascara solo para considerar aquellos CIPIs sin observaciones hoy pero que suelen medir.
    #No se vana  considerar aquellos casos donde no suele haber observaciones, en esos casos
    #va a ser True
    b_mask_obs[b_mask_hist == False] = True

    #########################
    #Interpolacion espacial #
    #########################

    #Itero sobre los puntos de interes
    for CIPI in dict_interpolacion.keys():

        #Plazos a interpolar
        PLAZOS = dict_interpolacion[(CIPI)]
        PLAZOS_array = np.array(PLAZOS)

        #Evaluo el coeficiente b y la mascara en los plazos que quiero interpolar
        #y a esta ultima la paso a array de numpy
        b_mask_plazo = np.squeeze(b_mask_tot.loc[(PLAZOS, slice(None))].values)
        b_plazo = b.loc[(slice(None), PLAZOS, slice(None)), :]

        #Calculo la distancia del punto de interes que quiero interpolar a todos los demas
        root = (np.sin((np.deg2rad(Puntos_interes_CIPI.loc[(CIPI), 'lat']) - np.deg2rad(Puntos_interes_CIPI['lat'].values))/2)**2
                ) + np.cos(np.deg2rad(Puntos_interes_CIPI.loc[(CIPI), 'lat']))*np.cos(np.deg2rad(Puntos_interes_CIPI['lat'].values))*(
                np.sin((np.deg2rad(Puntos_interes_CIPI.loc[(CIPI), 'lon']) - np.deg2rad(Puntos_interes_CIPI['lon'].values))/2)**2)

        d = 2*6371*np.arcsin(np.sqrt(root))

#        #Calculo el peso de cada punto de interes
#        try:
#            W = (d**(-2))
#        except FloatingPointError:
#            d[d == 0] = 1e-6
#            W = (d**(-2))

        #Pongo inf si la distancia es 0 para que los puntos de interes no se interpolen consigo mismo.
        #Al poner inf el peso va a ser 0
        d[d == 0] = np.inf

        #Calculo los coeficientes de la interpolacion
        coef1, coef2 = IDW_DF(PLAZOS, d, R, b_plazo, b_mask_plazo, nmiembros)

        #Si los coeficientes son 0 es altamente probable de que no haya tenido ningun CIPI
        #con el que interpolar
        #Busco los indices donde los coeficientes son 0
        indices = np.argwhere(np.logical_and(coef1 == 0, coef2 == 0))
        if len(indices) != 0:

            PLAZOS_cero = PLAZOS_array.repeat(nmiembros)[indices[:, 0]]

            #Selecciono los plazos en que el CIPI tiene historial de observaciones y para
            #algun dia en particular no tuvo observaciones.
            #En estos casos calibramos usando sus ultimos coeficientes disponibles
            indice_plazo_obs = np.argwhere(b_mask_obs.loc[(PLAZOS_cero, CIPI), :].values.repeat(nmiembros) == False)
            PLAZOS_obs = PLAZOS_cero[indice_plazo_obs[:, 0]]

            indices_obs = [i for i, ip in enumerate(PLAZOS_array.repeat(nmiembros)) if ip in PLAZOS_obs] 

            indices_obs_cero = [i for i, ip in enumerate(PLAZOS_cero) if ip in PLAZOS_obs] 

            if len(PLAZOS_obs) != 0:

                #Caso en que la estacion no midio hoy pero suele medir
                coef1[indices_obs] = b.loc[(slice(None), PLAZOS_obs, CIPI), :]['b0'].values.repeat(nmiembros)
                coef2[indices_obs] = b.loc[(slice(None), PLAZOS_obs, CIPI), :]['b1'].values.repeat(nmiembros)

            #Selecciono los plazos en los que el CIPI no tiene historial de observaciones.
            #Con lo que los coeficientes son siempre igual a 0
            #En estos casos calibramos interpolando los ultimos coeficientes disponibles de sus
            #CIPIs cercanos
            PLAZOS_no_obs = np.delete(PLAZOS_cero, indices_obs_cero)

            if len(PLAZOS_no_obs) != 0:

                indices_no_obs = [i for i, ip in enumerate(PLAZOS_array.repeat(nmiembros)) if ip in PLAZOS_no_obs] 

                b_mask_plazo = np.squeeze(b_mask_hist.loc[(PLAZOS_no_obs, slice(None))].values)
                b_plazo = b.loc[(slice(None), PLAZOS_no_obs, slice(None)), :]

                coef1_tmp, coef2_tmp = IDW_DF(np.unique(PLAZOS_no_obs), d, R, b_plazo, b_mask_plazo, nmiembros)
                coef1[indices_no_obs] = coef1_tmp
                coef2[indices_no_obs] = coef2_tmp

        forecast_today = DF.loc[(date_today, CICLO, PLAZOS, slice(None), CIPI, slice(None)), fcst_var]       
        #Si las variables pronosticadas son el viento calculo magnitud y direccion
        if fcst_var == ['Umet10', 'Vmet10']:
            fcst_today = np.sqrt(forecast_today['Umet10']**2 + forecast_today['Vmet10']**2).values
            dir_viento = np.arctan2(forecast_today['Umet10'], forecast_today['Vmet10'])
        else:
            fcst_today = forecast_today.values

        yp = coef1 + coef2*fcst_today
        fcst_corrected = fcst_today - yp

        #Guado la calibracion en el DF
        #Si la variable pronosticada es el viento descompongo la magnitud calibrada en U y V
        if fcst_var == ['Umet10', 'Vmet10']:
            fcst_corrected[fcst_corrected<0] = 0
            wind_dir = np.array([np.sin(dir_viento), np.cos(dir_viento)])
            DF.loc[(date_today, CICLO, PLAZOS, slice(None), CIPI, slice(None)), cal_var] = (fcst_corrected*wind_dir).T
        else:
            DF.loc[(date_today, CICLO, PLAZOS, slice(None), CIPI, slice(None)), cal_var] = fcst_corrected

    return DF, dict_interpolacion


def IDW_DF(PLAZOS, d, R, b, b2, nmiembros):
    """
    Funcion que tiene la parte central de la calibracion espacial y es la que
    habria que cambiar si se cambia el metodo de interpolacion
    """
    #Longitud del array de pesos
    len_d = len(d)

    #Expando el array de distancias en funcion de la cantidad de plazos a interpolar
    d_plazos = np.tile(d, len(PLAZOS))

    #Pongo la distancia en inf para los CIPIs y plazos que no calibran para que el peso resulte en 0
    d_plazos[np.logical_or(b2 == False, d_plazos>R)] = np.inf


    #Calculo los pesos de la interpolacion
    W = (d_plazos**(-2))

    #Paso los coeficientes b0 y b1 y el peso a matrices para poder operar
    #sobre todos los plazos al mismo tiempo
    b0 = b['b0'].values.reshape(len(PLAZOS), len_d)
    b1 = b['b1'].values.reshape(len(PLAZOS), len_d)
    W = W.reshape(len(PLAZOS), len_d)
    d_plazos = d_plazos.reshape(len(PLAZOS), len_d)

    #Calculo el peso de la extrapolacion.
    peso_distancia = linear_weight(d_plazos.min(axis = 1), R)

    #Calculo la suma de los pesos para cada tiempo
    W_sum = np.nansum(W, axis = 1)

    #Si la suma es 0 pongo infinito para que al dividir en los coeficientes estos den 0
    W_sum[W_sum == 0] = np.inf

    #Calculo los coeficientes
    coef1 = ((np.nansum(b0*W, axis = 1)/W_sum) * peso_distancia).repeat(nmiembros)
    coef2 = ((np.nansum(b1*W, axis = 1)/W_sum) * peso_distancia).repeat(nmiembros)

    return coef1, coef2


def linear_weight(d, R):
    """
    Funcion que calcula el peso con decaimiento lineal
    """
    weight = np.empty_like(d)
    for i, l in enumerate(d):
        if l < R:
            weight[i] = (R-l)/R
        else:
            weight[i] = 0

    return weight


def interp_dominio(FECHA_INI_obj, DIR_TEMPLATES, OUTDIR, MODELO, PLAZOS, dict_interpolacion, b, R, fcst_var, cal_var, nmiembros, tipo):
    """
    Funcion que interpola la calibracion puntual que se hace con la metodologia
    RAFK a todos los puntos de reticula de los modelos
    """
    #Cantidad de procesos a lanzar en la paralelizacion de la interpolacion
    process = int(os.environ['ICORE'])

    #Leo la tabla de puntos de interes y pongo CIPI como indice
    Puntos_interes = af.merge_Puntos_Interes(DIR_TEMPLATES, OUTDIR)
    Puntos_interes_CIPI = Puntos_interes.drop_duplicates(subset='CIPI')
    Puntos_interes_CIPI.set_index(['CIPI'], inplace = True)

    #Si la calibracion es horaria leo los POST del WRF o GRIB del GFS/GEFS y lo guardo
    #en un archivos. En caso de que sea minmax leo este archivo que se guardo para no
    #tener que leer de nuevo todos los POST o GRIB
    if tipo == 'horario':
        lats, lons, fcst = af.read_model_vars(FECHA_INI_obj, MODELO, PLAZOS, nmiembros, fcst_var, OUTDIR)
    elif tipo == 'viento':
        lats, lons, fcst_u = af.read_model_vars(FECHA_INI_obj, MODELO, PLAZOS, nmiembros, fcst_var[0], OUTDIR)
        lats, lons, fcst_v = af.read_model_vars(FECHA_INI_obj, MODELO, PLAZOS, nmiembros, fcst_var[1], OUTDIR)

        fcst = np.sqrt(fcst_u**2 + fcst_v**2)

        #Borro de la memoria los array fcst_u y fcst_v
        fcst_u = None
        fcst_v = None
        gc.collect()
    else:
        if os.path.isfile(OUTDIR + 'fcst_tmp.npy'):
            lats = np.load(OUTDIR + 'lat_tmp.npy')
            lons = np.load(OUTDIR + 'lon_tmp.npy')
            fcst = np.load(OUTDIR + 'fcst_tmp.npy')

            #Paso de los valores horarios a minimas y maximas
            fcst = af.calcula_minmax2D(fcst, tipo, FECHA_INI_obj, OUTDIR)
        else:
            #Asumimos que siempre la calibracion horaria va a haber corrido bien y el archivo
            #con los pronosticos va a estar generado
            print('No se encontro el archivo con los pronosticos horarios para la calibracion espacial')
            return

    #Los GRIB vienen con las longitudoes de 0 a 360 entonces paso de -180 a 180
#    if MODELO != 'WRF':
#        lons = lons - 360

    nlat, nlon = lats.shape

    if tipo == 'viento':
        fcst_var = 'WindSpeed'
        cal_var = 'WindSpeed_cal'

    #Eliminamos los CIPIs y plazos en donde los coeficientes son 0 ya que significa que no midieron al menos por los
    #ultimos dias
    b_nozero = b.where(b!=0).dropna()

    #Generamos una mascara con los CIPIs donde no hubo observacion hoy, inicialmente toda con True y luego
    #agregamos los False donde corresponde.
    b_mask_obs_DF = pd.DataFrame(np.full(len(b_nozero.index), True), index = b_nozero.index, columns = ['mask'])

    #Aniadimos False en la mascara donde no hubo observaciones para que
    #no se usen esos CIPIs y plazos en la interpolacion de alguna estacion
    for CIPI in dict_interpolacion.keys():
        try:
            b_mask_obs_DF.loc[(slice(None), dict_interpolacion[(CIPI)], CIPI), 'mask'] = False
        except KeyError:
            pass

    fcst_cal = deepcopy(fcst)

    #Buscamos los plazos donde efectivamente hay plazos para calibrar. En las minimas y maximas
    #dependiendo el ciclo pasa que hay plazos en que todos los coeficientes son 0 porque no se 
    #calibran
    PLAZOS_calib = b_nozero.index.get_level_values('Plazo').unique()

    if len(PLAZOS_calib) != 0:

        indices_calib = [PLAZOS.tolist().index(plazo) for plazo in PLAZOS_calib]

        arg_list = [(ip, 
                     b_nozero.loc[(slice(None), ip, slice(None)), :], 
                     b_mask_obs_DF.loc[(slice(None), ip, slice(None)), :], 
                     lats, 
                     lons, 
                     Puntos_interes_CIPI, 
                     R, 
                     fcst[:, indices_calib[i], :, :]) for i, ip in enumerate(PLAZOS_calib)]

        with mp.Pool(processes = min(process, len(PLAZOS_calib))) as pool:

            pool_out = pool.starmap(IDW_dominio, arg_list)

            for i, cal_plazo in enumerate(pool_out):
                #Buscamos el indice en PLAZOS al que corresponden los plazos efectivamente
                #interpolados
                plazo = PLAZOS_calib[i]
                ind_fcst = PLAZOS.tolist().index(plazo)

                fcst_cal[:, ind_fcst, :, :] = cal_plazo

            pool_out = None
    else:
        print('Todos los coeficientes de la calibracion son 0, no se realizo la interpolacion del dominio y el campo calibrado es igual al pronosticado')


    #Borro el array fcst de la memoria
    fcst = None
    gc.collect()

    if fcst_var == 'WindSpeed':
        fcst_cal[fcst_cal < 0] = 0

    #Nombre del archivo de salida
    filename = FECHA_INI_obj.strftime(OUTDIR + '%Y%m%d_%H0000_{}.nc'.format(cal_var))

    #Guardo el array de la variable calibrada en un netCDF
    ini = datetime.now()
    if MODELO == 'WRF':
        af.guarda_nc_interpolacion_WRF(fcst_cal, lats, lons, filename, FECHA_INI_obj, PLAZOS, cal_var)
    else:
        lats_recorte, lons_recorte, fcst_cal_recorte = af.recorta_GFS(fcst_cal,lats,lons)
        lons_recorte += 360 # Convertimos las longitudes a la convension del WRF
        af.guarda_nc_interpolacion_GFS(fcst_cal_recorte, lats_recorte, lons_recorte, filename, FECHA_INI_obj, PLAZOS, cal_var)
    fin = datetime.now()


def IDW_dominio(ip, b_plazo, b_mask_plazo, lats, lons, Puntos_interes_CIPI, R, fcst_plazo):
    """
    Funcion que aplica el metodo de interpolacion IDW a la calibracion de todo el dominio.
    """

    fcst_cal = np.empty_like(fcst_plazo)

    #Iteramos sobre los casos en que hubo observaciones hoy (True) y en los que no (False)
    #ya que la interpolacion usa distintos CIPIs y los coeficientes resultantes los guardamos
    #en el diccionario "coeficientes"
    coeficientes = {}
    for mask in [True, False]:

        b_mask_tmp = b_mask_plazo[b_mask_plazo == mask].dropna()

        lista_CIPI = b_mask_tmp.index.get_level_values('CIPI').unique()

        W = np.zeros_like(lats)
        Wb0 = np.zeros_like(lats)
        Wb1 = np.zeros_like(lats)
        
        #Generamos un campo para tener la minima distancia de cada ij del modelo a un CIPI
        min_dist = np.full_like(lats, np.inf)

        for CIPI in lista_CIPI:

            b0 = b_plazo.loc[(slice(None), slice(None), CIPI), 'b0'].values
            b1 = b_plazo.loc[(slice(None), slice(None), CIPI), 'b1'].values

            #Calculamos la distancia del CIPI a cada ij del dominio
            root = (np.sin((np.deg2rad(lats) - np.deg2rad(Puntos_interes_CIPI.at[CIPI, 'lat']))/2)**2) + np.cos(
            np.deg2rad(lats))*np.cos(np.deg2rad(Puntos_interes_CIPI.at[CIPI, 'lat']))*(
            np.sin((np.deg2rad(lons) - np.deg2rad(Puntos_interes_CIPI.at[CIPI, 'lon']))/2)**2)

            d = 2*6371*np.arcsin(np.sqrt(root))

            #Calculamos los pesos de la interpolacion. En los casos que d>R ponemos inf para que los pesos den 0
            d[d>R] = np.inf
            W_CIPI = (1/d)**2 # NOTA, es un campo

            #Sumamos los pesos de cada CIPI y multiplicamos el peso por el coeficiente
            W += W_CIPI
            Wb0 += W_CIPI*b0
            Wb1 += W_CIPI*b1

            #Calculamos el minimo entre min_dist y d
            min_dist = np.min([min_dist, d], axis = 0)

        W_min_dist = (R-min_dist)/R
        W_min_dist[W_min_dist < 0] = 0

        #Calculamos los coeficientes interpolados
        W[W==0] = np.inf #Para evitar dividir por 0
        coeficientes[mask, 'b0'] = (Wb0/W)*W_min_dist
        coeficientes[mask, 'b1'] = (Wb1/W)*W_min_dist


        #Iteramos de nuevo sobre todos los CIPIs para poner el valor de b0 y b1 del CIPI en el i, j
        #del modelo asociado a este
        for CIPI in lista_CIPI:
            i, j = Puntos_interes_CIPI.loc[(CIPI), ['i', 'j']]
            b0 = b_plazo.loc[(slice(None), slice(None), CIPI), 'b0'].values
            b1 = b_plazo.loc[(slice(None), slice(None), CIPI), 'b1'].values

            coeficientes[mask, 'b0'][i, j] = b0
            coeficientes[mask, 'b1'][i, j] = b1

    ordenada = coeficientes[True, 'b0']
    ordenada[ordenada == 0] = coeficientes[False, 'b0'][ordenada == 0]
    pendiente = coeficientes[True, 'b1']
    pendiente[pendiente == 0] = coeficientes[False, 'b1'][pendiente == 0]

    fcst_cal[:, :, :] = fcst_plazo[:, :, :] - (ordenada + fcst_plazo[:, :, :] * pendiente)

    return fcst_cal


