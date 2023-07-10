import numpy as np
import pandas as pd
from netCDF4 import Dataset
import pygrib
import sys
import os
from datetime import datetime, timedelta
import pyproj
import glob
import gc

WRFUTILDIR = os.environ['WRFUTILDIR']
DIR_TEMPLATES = WRFUTILDIR + '/templates/'
sys.path.append(DIR_TEMPLATES)

import gribNames
from atributos_nc_interp import atributos as attr


def determino_ij(FILEIN, MODELO, Tabla_CIPI, FILEOUT):#, OUTDIR):
    """
    Funcion para determinar los i,j de cada estacion en la reticula
    del modelo
    """
    n_cercanos = 16 #Cantidad de puntos cercanos que busco

    Tabla_CIPI.set_index('CIPI', inplace = True)

    if MODELO in ['GFS', 'GEFS']:
        Tabla_CIPI['lon'] = Tabla_CIPI['lon'] + 360


    #Si el archivo de salida ya existe se borra, de otro modo los datos se escribirian al
    #final de este
    if os.path.isfile(FILEOUT):
        os.remove(FILEOUT)   

    #obtengo los array de lat, lon y mascara
    lats, lons, mascara = read_model(FILEIN, MODELO)

    #Pongo el encabezado del archivo de salida
    with open(FILEOUT, 'a', encoding = 'latin-1') as fo:
        fo.write("CIPI,i,j,Superficie\n")

    #Itero sobre el archivo con las ubicaciones
    for CIPI in Tabla_CIPI.index:

        lat, lon = Tabla_CIPI.loc[CIPI]

        if (lat < lats.min() or lat > lats.max() or lon < lons.min() or lon > lons.max()):
            continue


        dist = (lats-lat)**2+(lons-lon)**2
        for k in range(n_cercanos):
            #Busco la minima distancia
            dist_min = np.argwhere(dist == dist.min())
            if k == 0:
                primer_dist_min = dist_min

            #Chequeo si el punto cae en la tierra o en el agua
            land = mascara[dist_min[0, 0], dist_min[0, 1]]

            if land == 0:
                #Si el punto es agua le asigno inf y busco otro punto cercano
                dist[dist_min[0,0],dist_min[0,1]] = np.inf
                dist[dist_min[0,0],dist_min[0,1]] = np.inf
            else:
                superficie = 'Tierra'
                break
                
            if k == n_cercanos-1:
                dist_min = primer_dist_min
                superficie = 'Agua'


        #Guardo los datos en el archivo de salida
        with open(FILEOUT, 'a', encoding = 'latin-1') as fo:
            fo.write("{},{},{},{}\n".format(CIPI, dist_min[0,0], dist_min[0,1], superficie))


def read_model(FILEIN, MODELO):
    """
    Funcion que llama determino_ij para leer la salida de los modelos
    y obtener las matrices de latitud, longitud y mascara de tierra y agua
    """
    if MODELO == 'WRF':
        try:
            f = Dataset(FILEIN)
        except:
            print('Fallo en la lectura del NETCDF')
            sys.exit(1)

        lats = np.squeeze(f.variables['XLAT'][:])
        lons = np.squeeze(f.variables['XLONG'][:])
        mascara = np.squeeze(f.variables['LANDMASK'][:])
        f.close()
        
    else:
        variable = 'Land-sea mask'

        try: 
            f = pygrib.open(FILEIN)
        except:
            print('Fallo en la lectura del GRIB')
            sys.exit(1)
        mensaje = f.select(name = variable)[0]
        mascara = mensaje.values
        lats, lons = mensaje.latlons() 
        f.close()

    return lats, lons, mascara


def read_model_vars(FECHA_INI_obj, MODELO, PLAZOS, nmiembros, fcst_var, OUTDIR):
    """
    Funcion que llama interp_dominio  para leer la salida de los modelos
    y obtener las matrices de latitud, longitud y temperatura
    """
    if MODELO == 'GFS':
        var_names = gribNames.GFS
    elif MODELO == 'GEFS':
        var_names = gribNames.GEFS
    else:
        #Para WRF no se necesita esta variable
        pass

    #Directorio base donde buscar la salida del modelo
    INPDIR = os.environ['DIRIN']

    FECHA = FECHA_INI_obj.strftime('%Y%m%d')
    CICLO = FECHA_INI_obj.strftime('%H')

    #Genero la lista de miembros y el path a la salida del modelo
    if nmiembros == 1:
        miembros = [0]
        if MODELO == 'WRF':
            INPFILE = INPDIR + '{0:0=2d}/model.WRF_DET_4km.{1}_{2}0000.{3:0=3d}.OPERSFC.nc'
        else:
            INPFILE = INPDIR + '/gfs.t{2}z.pgrb2.0p25.f{3:0=3d}'
    else:
        miembros = [x for x in range(1, nmiembros+1)]
        if MODELO == 'WRF':
            INPFILE = INPDIR + '{0:0=2d}/model.WRF_ENS_4km.{1}_{2}0000.{3:0=3d}.M{0:0=2d}_OPERSFC.nc'
        else:
            INPFILE = INPDIR + '{0:0=2d}/gep{0:0=2d}.t{2}z.pgrb2.0p50.f{3:0=3d}'

    #Itero sobre los miembros
    for im, mem in enumerate(miembros):
        #Itero sobre los plazos
        for jp, ip in enumerate(PLAZOS):

            #Path a la salida del miembro im y plazo ip
            filename = INPFILE.format(mem, FECHA, CICLO, ip)
            if not os.path.isfile(filename):
                continue

            if MODELO == 'WRF':
                #Abro el nc
                nc = Dataset(filename)

                #En el primer miembros y plazo leo las latlons y genero el array que va a tener todos los datos
                if mem == miembros[0] and ip == PLAZOS[0]:
                    lat = nc.variables['XLAT'][:].data
                    lon = nc.variables['XLONG'][:].data
                    fcst = np.full([nmiembros, len(PLAZOS), lat.shape[0], lat.shape[1]], np.nan)
            
                #Almaceno el dato pronosticado
                fcst[im, jp, :, :] = np.squeeze(nc.variables[str(fcst_var)][:])
        
                nc.close()
            else:
                #Abro el grib
                grib = pygrib.open(filename)
                #Selecciono la variable pronosticada
                mensaje = grib.select(name = var_names[fcst_var])[0]

                #En el primer miembros y plazo leo las latlons y genero el array que va a tener todos los datos
                if mem == miembros[0] and ip == PLAZOS[0]:
                    lat, lon = mensaje.latlons()

                    #GFS y GEFS tienen como dominio el globo entonces recorto
#                    lon_max = -35 + 360 #Longitud mas alta del WRF
#                    lon_min = -95 + 360 #Longitud mas baja del WRF
#                    lat_max = -11 #Latitud mas alta del WRF
#                    lat_min = -80 #Latitud para que en el dominio esten todas las estaciones antarticas

#                    i_min = np.argwhere(lat_max == lat)[0, 0]
#                    i_max = np.argwhere(lat_min == lat)[0, 0] + 1
#                    j_min = np.argwhere(lon_min == lon)[0, 1]
#                    j_max = np.argwhere(lon_max == lon)[0, 1] + 1

#                    lat = lat[i_min:i_max, j_min:j_max]
#                    lon = lon[i_min:i_max, j_min:j_max]

                    fcst = np.full([nmiembros, len(PLAZOS), lat.shape[0], lat.shape[1]], np.nan)

                #Almaceno el dato pronosticado
                fcst[im, jp, :, :] = mensaje.values #[i_min:i_max, j_min:j_max]

                grib.close()

    #Guardo las matrices de temperatura, lat y lon para calcularlas una sola vez y luego 
    #levantar el dato guardo
    if fcst_var == 'T2':
        np.save(OUTDIR + 'fcst_tmp.npy', fcst)
        np.save(OUTDIR + 'lat_tmp.npy', lat)
        np.save(OUTDIR + 'lon_tmp.npy', lon)

    return lat, lon, fcst


def calcula_minmax2D(fcst, tipo, FECHA_INI_obj, OUTDIR):
    """
    Funcion que a partir de la matriz de temperatura calcula las minimas y las 
    maximas del dia
    """
    #Sacamos de un DF horario la validez de cada plazo
    DF_horario = pd.read_hdf(FECHA_INI_obj.strftime(OUTDIR + '%Y%m%d_%H0000_horario.h5'))
    validez = DF_horario.index.get_level_values('Validez').unique()

    #Dias para los que hay pronostico
    dias = np.unique(validez.date)

    #Creo la matriz que tiene las minimas o maximas
    fcst_tipo = np.full([fcst.shape[0], len(dias), fcst.shape[2], fcst.shape[3]], np.nan)

    for i, day in enumerate(dias):
        day_datetime = datetime(day.year, day.month, day.day)

        if tipo == 'Tmin':
            minmax_range = [day_datetime + timedelta(hours = x) for x in range(0, 13)]
        else:
            minmax_range = [day_datetime + timedelta(hours = x) for x in range(12, 25)]

        if not (minmax_range[0] in validez and minmax_range[-1] in validez):
            continue

        ini = np.argwhere(validez == minmax_range[0])[0, 0]
        fin = np.argwhere(validez == minmax_range[-1])[0, 0]

        #Calculo las minimas y maximas
        if tipo == 'Tmin':
            fcst_tipo[:, i, :, :] = fcst[:, ini:fin + 1, :, :].min(axis = 1)
        else:
            fcst_tipo[:, i, :, :] = fcst[:, ini:fin + 1, :, :].max(axis = 1)


    return fcst_tipo       

def add_fill_value(data, fill_value):
    data = np.ma.fix_invalid(data, fill_value=fill_value)
    np.ma.set_fill_value(data, fill_value)
  
    return data

def guarda_nc_interpolacion_WRF(fcst, lat, lon, filename, fecha_inicio, PLAZOS, cal_var):
    """
    Funcion que guarda el campo de temperatura calibrado en un netCDF
    """

    fill_value = 1e20

    #Path de salida
    OUTDIR = os.environ['DIRIN']

    #Leo algun POST para obtener atributos relacionados a la proyeccion del WRF
    FILEIN = glob.glob(OUTDIR + '/*/model*')[0]

    POST = Dataset(FILEIN)

    #Genero la lista de miembros
    if fcst.shape[0] == 1:
        mems = [0]
        tipo = 'Deterministic'
    else:
        mems = [x for x in range(1, fcst.shape[0] + 1)]
        tipo = 'Ensemble'


    #Vector de fechas que se va a guardar en el nc para recontruir la fecha y hora de cada plazo
    fechas = (PLAZOS*60).values.astype(float)

    #Genero los vectores x e y que representan las coordenadas en la proyeccion

    #Generacion del x e y tomado de https://fabienmaussion.info/2018/01/06/wrf-projection/
    lcc = pyproj.Proj(proj='lcc', lat_1 = POST.TRUELAT1, lat_2 = POST.TRUELAT2, # Cone intersects with the sphere
                      lat_0 = POST.CEN_LAT, lon_0 = POST.CEN_LON, a=6370000, b=6370000)

    wgs_proj = pyproj.Proj(proj='latlong', datum='WGS84')
    e, n = pyproj.transform(wgs_proj, lcc, POST.CEN_LON, POST.CEN_LAT)

    dx, dy = POST.DX, POST.DY
    nx, ny = len(POST.dimensions['x']), len(POST.dimensions['y'])

    x0 = -(nx-1) / 2. * dx + e
    y0 = -(ny-1) / 2. * dy + n

    x = (np.arange(nx) * dx + x0)
    y = (np.arange(ny) * dy + y0)


    #Traspongo para que el tiempo quede como primera dimension
    fcst = fcst.transpose((1, 0, 2, 3))

    #Creo que el NC que va a guardar el campo calibrado
    NC = Dataset(filename, 'w') 

    NC.title = attr[cal_var]["title"].format("WRF")
    NC.institution = "Servicio Meteorologico Nacional"
    NC.creation_date = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    NC.Conventions = "CF-1.8"

    #Creo las dimensiones del nc
    NC.createDimension('Member', len(mems))
    NC.createDimension('Time', None)
    NC.createDimension('x', x.size)
    NC.createDimension('y', y.size)

    #Variable con los datos de la proyeccion
    outvar = NC.createVariable('Lambert_Conformal', "c", (), fill_value = fill_value)
    outvar.grid_mapping_name = 'lambert_conformal_conic'
    outvar.standard_parallel = POST.TRUELAT1
    outvar.latitude_of_projection_origin = POST.CEN_LAT
    outvar.longitude_of_central_meridian = POST.CEN_LON

    #Variable con los miembros
    outvar = NC.createVariable('Member', np.int32, ('Member'), zlib=True)
    outvar[:] = mems
    outvar.standard_name = 'realization'
    outvar.units = '1'
    outvar.long_name = 'Ensemble members'

    #Variable con el array de fechas
    outvar = NC.createVariable('Time', np.float32, ('Time'), zlib=True, fill_value = fill_value)
    outvar[:] = add_fill_value(fechas, fill_value)
    outvar.long_name = fecha_inicio.strftime('minutes since %Y-%m-%d %H:%M:%S')
    outvar.units = fecha_inicio.strftime('minutes since %Y-%m-%d %H:%M:%S')
    outvar.standard_name = 'time'
    outvar.axis = 'T'
    outvar.calendar = 'gregorian'

    #Variable con la coordenada y de la proyeccion
    outvar = NC.createVariable('y', np.float32, ('y'), zlib=True, fill_value = fill_value)
    outvar[:] = add_fill_value(y, fill_value)
    outvar.axis = 'Y'
    outvar.units = 'm'
    outvar.standard_name = 'projection_y_coordinate'
    outvar.long_name = 'y-coordinate in projected coordinate system'

    #Variable con la coordenada x de la proyeccion
    outvar = NC.createVariable('x', np.float32, ('x'), zlib=True, fill_value = fill_value)
    outvar[:] = add_fill_value(x, fill_value)
    outvar.axis = 'X'
    outvar.units = 'm'
    outvar.standard_name = 'projection_x_coordinate'
    outvar.long_name = 'x-coordinate in projected coordinate system'

    #Variable con la matriz de latitud
    outvar = NC.createVariable('XLAT', np.float32, ('y', 'x'), zlib=True, fill_value = fill_value)
    outvar[:,:] = add_fill_value(lat[:, :], fill_value)
    outvar.long_name = "Latitude"
    outvar.units = "degrees_north"
    outvar.standard_name = 'latitude'

    #Variable con la matriz de longitud
    outvar = NC.createVariable('XLONG', np.float32, ('y', 'x'), zlib=True, fill_value = fill_value)
    outvar[:,:] = add_fill_value(lon[:, :], fill_value)
    outvar.long_name = "Longitude"
    outvar.units = "degrees_east"
    outvar.standard_name = 'longitude'

    #Variable con la matriz de temperatura calibrada
    #outvar = NC.createVariable('T2cal', np.float32, ('Member', 'Time', 'y', 'x'), zlib=True)
    outvar = NC.createVariable(cal_var, np.float32, ('Time', 'Member', 'y', 'x'), zlib=True, fill_value = fill_value)
    outvar[:, :, :, :] = add_fill_value(fcst[:, :, :, :], fill_value)
    outvar.long_name = attr[cal_var]["long_name"]
    outvar.units = attr[cal_var]["units"]
    outvar.standard_name = attr[cal_var]["standard_name"]
    outvar.coordinates = "XLAT XLONG"
    outvar.description = attr[cal_var]["description"]
    outvar.grid_mapping = "Lambert_Conformal"   

    #Cierro el nc
    NC.close()

def guarda_nc_interpolacion_GFS(fcst, lat, lon, filename, fecha_inicio, PLAZOS, cal_var):
    """
    Funcion que guarda el campo de temperatura calibrado en un netCDF si el modelo
    utilizado no es el WRF.

    La diferencia con el WRF es que en este caso las reticulas de latlon es regular,
    entonces la convencion CF pide que hay que pasarlas como vector y no generar
    las dimensiones x e y que son las coordenadas en la proyeccion.
    """

    fill_value = 1e20

    #Vectores con las latitudes y longitudes
    lat = lat[:, 0]
    lon = lon[0, :]

    #Genero la lista de miembros
    if fcst.shape[0] == 1:
        mems = [0]
        modelo = 'GFS'
    else:
        mems = [x for x in range(1, fcst.shape[0] + 1)]
        modelo = 'GEFS'

    #Vector de fechas que se va a guardar en el nc para recontruir la fecha y hora de cada plazo
    fechas = (PLAZOS*60).values.astype(float)

    #Traspongo para que el tiempo quede como primera dimension
    fcst = fcst.transpose((1, 0, 2, 3))

    #Creo que el NC que va a guardar el campo calibrado
    NC = Dataset(filename, 'w') 

    #Agrego atributos globales
    NC.title = attr[cal_var]["title"].format(modelo)
    NC.institution = "Servicio Meteorologico Nacional"
    NC.creation_date = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    NC.Conventions = "CF-1.8"

    #Creo las dimensiones del nc
    NC.createDimension('Member', len(mems))
    NC.createDimension('Time', None)
    NC.createDimension('lon', lon.size)
    NC.createDimension('lat', lat.size)

    #Variable con los miembros
    outvar = NC.createVariable('Member', np.int32, ('Member'), zlib=True)
    outvar[:] = mems
    outvar.standard_name = 'realization'
    outvar.units = '1'
    outvar.long_name = 'Ensemble members'

    #Variable con el array de fechas
    outvar = NC.createVariable('Time', np.float32, ('Time'), zlib=True, fill_value = fill_value)
    outvar[:] = add_fill_value(fechas, fill_value)
    outvar.long_name = fecha_inicio.strftime('minutes since %Y-%m-%d %H:%M:%S')
    outvar.units = fecha_inicio.strftime('minutes since %Y-%m-%d %H:%M:%S')
    outvar.standard_name = 'time'
    outvar.axis = 'T'
    outvar.calendar = 'gregorian'

    #Variable con la matriz de latitud
    outvar = NC.createVariable('lat', np.float32, ('lat'), zlib=True, fill_value = fill_value)
    outvar[:] = add_fill_value(lat[:], fill_value)
    outvar.long_name = "Latitude"
    outvar.units = "degrees_north"
    outvar.standard_name = 'latitude'

    #Variable con la matriz de longitud
    outvar = NC.createVariable('lon', np.float32, ('lon'), zlib=True, fill_value = fill_value)
    outvar[:] = add_fill_value(lon[:], fill_value)
    outvar.long_name = "Longitude"
    outvar.units = "degrees_east"
    outvar.standard_name = 'longitude'

    #Variable de temperatura calibrada
    outvar = NC.createVariable(cal_var, np.float32, ('Time', 'Member', 'lat', 'lon'), zlib=True, fill_value = fill_value)
    outvar[:, :, :, :] = add_fill_value(fcst[:, :, :, :], fill_value)
    outvar.long_name = attr[cal_var]["long_name"]
    outvar.units = attr[cal_var]["units"]
    outvar.standard_name = attr[cal_var]["standard_name"]
    outvar.description = attr[cal_var]["description"]
    
    #Cierro el nc
    NC.close()


def WRF2DF(DF_mem, INPDIR, FECHA_INI_obj, NPlazo, variables, puntos_interes):
    """
    Funcion que llama armo_DF para leer los datos del WRF e ir completando los dataframes
    """

    #Eliminamos de la lista los puntos CIPI repetidos
    puntos_interes.drop_duplicates('CIPI', inplace = True)

    m = DF_mem.index.get_level_values('Miembro').unique()[0]
    print('Procesando el miembro {}'.format(m))
    FECHA_INI = FECHA_INI_obj.strftime('%Y%m%d')
    CICLO = FECHA_INI_obj.strftime('%H')
    mem = str(m).zfill(2)
    if m == 0:
        prefix = 'model.WRF_DET_4km.'
        sufix = '.OPER'
    else:
        prefix = 'model.WRF_ENS_4km.'
        sufix = '.M' + mem + '_OPER' 

    PP_prev = np.full(puntos_interes.shape[0], 0.)
    for ip in NPlazo.tolist():
        PLAZO = str(ip).zfill(3)
        val = datetime.strftime(FECHA_INI_obj + timedelta(hours = ip), '%Y-%m-%d %H:00:00')
        # Salidas del WRF deterministico:
        INPFILE = mem + '/' + prefix + FECHA_INI + '_' + CICLO + '0000.' + PLAZO + sufix + 'SFC.nc'

        if not os.path.isfile(INPDIR + INPFILE) and m==0:
            sufix = '.HIST'
            INPFILE = mem + '/' + prefix + FECHA_INI + '_' + CICLO + '0000.' + PLAZO + sufix + 'SFC.nc'
        try:
            f = Dataset(INPDIR+INPFILE)
            lista = f.variables.keys()  # muestra la lista de variables que tiene el .nc
            for ivarname in variables:
                if (ivarname in lista): # si estan todas las variables en el archivo .nc
                    ivarval = np.squeeze(f.variables[ivarname][:])
                    for indice, row in puntos_interes.iterrows(): # loop para todas las estaciones: 
                        i = row['i']
                        j = row['j']
                        numero = row['CIPI']
                        DF_mem.loc[(FECHA_INI, CICLO, ip, val, numero, m), ivarname] = ivarval[i,j]
                    if ivarname == 'PP':
                        DF_mem.loc[(FECHA_INI, CICLO, ip, val, slice(None), m), ivarname] -= PP_prev
                        PP_prev += DF_mem.loc[(FECHA_INI, CICLO, ip, val, slice(None), m), ivarname].values


            f.close()
        except Exception as e:
            print(e)
            pass

    return DF_mem


def GFS2DF(DF_mem, INPDIR, FECHA_INI_obj, NPlazo, variables, puntos_interes):
    """
    Funcion que llama armo_DF para leer los datos del GFS e ir completando los dataframes
    """

    #Eliminamos de la lista los puntos CIPI repetidos
    puntos_interes.drop_duplicates('CIPI', inplace = True)

    m = DF_mem.index.get_level_values('Miembro').unique()[0]
    print('Procesando el miembro {}'.format(m))
    FECHA_INI = FECHA_INI_obj.strftime('%Y%m%d')
    CICLO = FECHA_INI_obj.strftime('%H')
    mem = str(m).zfill(2)

    PP_prev = np.full(puntos_interes.shape[0], 0.)
    for ip in NPlazo.tolist():
        #print(ip)
        PLAZO = str(ip).zfill(3)
        val = datetime.strftime(FECHA_INI_obj + timedelta(hours = ip), '%Y-%m-%d %H:00:00')
        if mem == '00':
            INPFILE = '/gfs.t' + CICLO + 'z.pgrb2.0p25.f' + PLAZO
        else:
            INPFILE = mem + '/gep' + mem + '.t' + CICLO + 'z.pgrb2.0p50.f' + PLAZO

        try:
            f = pygrib.open(INPDIR + INPFILE)
        except Exception as e:
            print(e)
            continue

        for var in variables.keys():
            if var == 'PP' and ip == 0:
                continue
            ivarval = f.select(name = variables[var])[0].values
            if var == 'PSFC':
                ivarval = ivarval/100
            for indice, row in puntos_interes.iterrows():
                i = row['i']
                j = row['j']
                numero = row['CIPI']
                DF_mem.loc[(FECHA_INI, CICLO, ip, val, numero, m), var] = ivarval[i,j]
            if var == 'PP':
                if m == 0: #Para GFS
                    DF_mem.loc[(FECHA_INI, CICLO, ip, val, slice(None), m), var] -= PP_prev
                    PP_prev += DF_mem.loc[(FECHA_INI, CICLO, ip, val, slice(None), m), var].values
                else: #Para GEFS
                    if np.mod(ip, 6) != 0 or ip > 192: 
                        PP_prev = DF_mem.loc[(FECHA_INI, CICLO, ip, val, slice(None), m), var].values
                    else:
                        PP = DF_mem.loc[(FECHA_INI, CICLO, ip, val, slice(None), m), var].values - PP_prev
                        PP[PP<0] = 0    #Hay casos con mas PP acumulada en 3 horas que en 6
                        DF_mem.loc[(FECHA_INI, CICLO, ip, val, slice(None), m), var] = PP

        f.close()

    return DF_mem


def open_obs(filename,date,hour):
    """ Función para leer los synop y armar un DataFrame con todas las variables del SYNOP para las estaciones de Arg
    (toma cosas programadas por Federico Cutraro)
    ####################################################################
    Keywords arguments:
    filename: Nombre del archivo
    date: Fecha actual string AAAAMMDD (!!!! VER BIEN ESTA FECHA A QUE DIA CORRESPONDE)
    hour: Hora actual
    Columnas del archivo: 
    fecha// hora // estación // latitud(-90:90) // longitud(0:360) // altura(m) // temperatura(°C) // temperatura de rocío(°C) // presión(hPa) // presión a nivel del mar(hPa) // dirección del viento(°) // intensidad del viento(kts) // humedad relativa(%)

    """

    names = ['fecha','hora','estacion','latitud','longitud','altura',
    'temperatura','temperaturaRocio','presionEstacion','presionNMM',
    'direccionViento','intensidadViento','humedadRelativa']

    data_np = np.genfromtxt(filename,dtype=float,delimiter=',')
    data = pd.DataFrame(data_np,columns = names) # Lo paso a un dataframe

    # El archivo contiene datos de las últimas 3 horas 
    # Solo quedan los datos de la hora que se quiere actualizar
    data2 = data[(data['hora'] == int(hour))]

    # Si no hay datos, sale del script
    if data2.size==0:
      raise
    # El archivo se va sobrescribiendo todos los días
    # Se chequea que los datos sean del dia en curso
    fecha_data = str(int(data2['fecha'].iloc[0]))

    # Si no corresponde la fecha, sale del script
    if(date!=fecha_data):
      raise

    # Loop para pasar por todas las estaciones
    for i in data2.index:
        k = i - data2.index[0]
        estacion = data2['estacion'][i]
        lat = data2['latitud'][i]
        lon = data2['longitud'][i]
        altura = data2['altura'][i]
        T = data2['temperatura'][i] + 273.15
        Td = data2['temperaturaRocio'][i] + 273.15
        Psup = data2['presionEstacion'][i]
        PNMM = data2['presionNMM'][i]
        dirViento = data2['direccionViento'][i]
        magViento = data2['intensidadViento'][i]*0.514444
        HR = data2['humedadRelativa'][i]
        tmp = pd.DataFrame([[int(estacion), lon, lat, altura, T, Td, Psup, PNMM, dirViento, magViento, HR]],
                           columns=['ID','Longitud','Latitud','Altura','T2obs','TD2obs','Pobs','PNMMobs','dirVientoObs','magVientoObs','HRobs'])
        if k==0:
            data3 = tmp
        else:
            data3 = pd.concat([data3,tmp],axis=0) # Agrega la nueva fila al final

    data_final = data3

    return data_final

def condicion_completado(val, valideces, CICLO, resto):
    """
    Esta funcion define si se completan observaciones de minimas y maximas
    segun el ciclo y la validez
    """
    minima = True
    maxima = True

    if (val == valideces[0] and CICLO != '00'):
        minima = False

    if (val == valideces[0] and CICLO == '18'):
        maxima = False

    if (val == valideces[-1] and resto < 12):
        minima = False
        maxima = False

    if (val == valideces[-1] and resto >= 12):
        maxima = False

    return minima, maxima


def merge_Puntos_Interes(DIR_TEMPLATES, OUTDIR = None, FILEIN = False, MODELO = False, estimo_ij = False):
    """
    Funcion que mergea las listas relacionadas a los puntos de interes y en caso de que
    se pasen los argumentos necesarios calcula los ij
    """

    #Se chequea que si estan todos los argumentos para calcular los ij
    if (estimo_ij) and ((not FILEIN) or (not MODELO) or (not OUTDIR)):
        print('Faltan datos para calular los ij')
        sys.exit(10)
        
    #Leo las tablas
    Tabla_CIPI = pd.read_csv(DIR_TEMPLATES + 'CIPI.dat', encoding = 'latin-1')
    Tabla_Tipo = pd.read_csv(DIR_TEMPLATES + 'Tipo.dat', encoding = 'latin-1')
    Tabla_PInteres = pd.read_csv(DIR_TEMPLATES + 'PInteres.dat', encoding = 'latin-1')

    #Mergeo
    Puntos_interes = Tabla_PInteres.merge(Tabla_Tipo, on = 'tipo').merge(Tabla_CIPI, on = 'CIPI')

    #Nombre del archivo que almacena los ij para cada CIPI
    filename_ij = OUTDIR + 'CIPI2ij.csv'

    #Si es True estimo los ij
    if estimo_ij:
        determino_ij(FILEIN, MODELO, Tabla_CIPI, filename_ij)#, OUTDIR)
            
    #Leo el archivo con los ij y mergeo con la tabla ya creada para quedarnos
    #solo con los Puntos de interes que estan dentro del dominio del modelo
    ij_Puntos_interes = pd.read_csv(filename_ij)
    Puntos_interes = Puntos_interes.merge(ij_Puntos_interes, on = 'CIPI')

    return Puntos_interes


def recorta_GFS(fcst,lat,lon):
    """ Recorta el dominio de GFS y GEFS 
    """
    lon_max = -35 + 360 #Longitud mas alta del WRF
    lon_min = -95 + 360 #Longitud mas baja del WRF
    lat_max = -11 #Latitud mas alta del WRF
    lat_min = -80 #Latitud para que en el dominio esten todas las estaciones antarticas

    i_min = np.argwhere(lat_max == lat)[0, 0]
    i_max = np.argwhere(lat_min == lat)[0, 0] + 1
    j_min = np.argwhere(lon_min == lon)[0, 1]
    j_max = np.argwhere(lon_max == lon)[0, 1] + 1

    lat = lat[i_min:i_max, j_min:j_max]
    lon = lon[i_min:i_max, j_min:j_max]
    
    fcst = fcst[:,:,i_min:i_max, j_min:j_max]

    return lat, lon, fcst

