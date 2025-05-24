###################################################
# It needs cfgrib >= 0.9.10.0 to work as expected #
###################################################

import xarray as xr
from metpy.units import units
from metpy import constants
import re
from datetime import datetime, timedelta
import sys
import numpy as np
#import cfgrib
import os
from netCDF4 import Dataset
import dask
import wrf
from wrf import interplevel, destagger

#Para conservar los atributos al hacer operaciones entre Datasets/DataArrays
xr.set_options(keep_attrs = True)


#Agrego dBZ al registro de unidades
#### units.define('dBZ = dBZ')


def smn_netcdf_from_list(lista, variables, catalog, levels, filetype):
   """
   Funcion que procesa varios archivos en paralelo y los concatena en un solo Dataset

   lista: lista con los nombres de los archivos a procesar
   variables: diccionario con el nombre de las variables a procesar y la coordenada vertical asociada (None si es 2D)
   catalog: catalogo de variables
   levels: diccionario con las coordenadas verticales y los niveles a extraer (None si se quieren todos los niveles)
   filetype: nombre del modelo a procesar (WRF/GFS)
   """
   self = sys.modules[__name__]

   #Seteo variables de entorno para evitar que las librerias usen multi-threading
   #https://docs.dask.org/en/stable/array-best-practices.html#avoid-oversubscribing-threads
   os.environ["OMP_NUM_THREADS"] = "1"
   os.environ["MKL_NUM_THREADS"] = "1"
   os.environ["OPENBLAS_NUM_THREADS"] = "1"

   try:
      ICORE = int(os.environ['ICORE'])
   except:
      ICORE = 1

   delayed_list = []
   for filename in lista:
#      print("===CHECK=== call dask.delayed " + filename)
      ds_delayed = dask.delayed(smn_netcdf_from_file)(filename, variables, catalog, levels, filetype)
      delayed_list.append(ds_delayed)

   dask_list = dask.compute(delayed_list, scheduler = 'multiprocessing', num_workers = ICORE)[0]

#   print("===CHECK=== ds num 0" + dask_list[0])
#   print("===CHECK=== ds.START_DATE " + dask_list[0].START_DATE)


   ds_list = []
   for ds_tmp in dask_list:
      if ds_tmp:
         ds_tmp = add_dimensions(ds_tmp)
         ds_list.append(ds_tmp)

   ds = xr.combine_by_coords(ds_list, combine_attrs = 'drop_conflicts')
   
   ds = standard_coords(ds, catalog)
   ds = standard_dataset(ds, catalog)

   if ds:
      #Hay casos en que los atributos 'FCST_HOURS', 'VALID_DATE' y 'MEMBER' no se eliminan despues del combine_by_coords,
      #los elimino en los casos necesarios
      if len(ds['XTIME']) > 1:
         rm_attr = ['FCST_HOURS', 'VALID_DATE']
         for attr in rm_attr:
            ds.attrs.pop(attr, None)

      if 'MEMBER' in ds.coords:
         if len(ds['MEMBER']) > 1:
            ds.attrs.pop('MEMBER', None)


   return ds


def add_dimensions(ds):
   """
   Funcion que agrega las dimensiones asociadas al tiempo y los miembros

   ds: dataset a agregarle las dimensiones
   """

   if 'XTIME' not in ds.dims:
      ds = ds.expand_dims(XTIME = [ds['XTIME'].values])

   if 'MEMBER' in ds.coords:
      ds = ds.expand_dims(MEMBER = [ds['MEMBER'].values])

   return ds


def smn_netcdf_from_file(filename, variables, catalog, levels, filetype):
   """
   Funcion que procesa solamente un archivo

   filename: nombre del archivo a procesar
   variables: diccionario con el nombre de las variables a procesar y la coordenada vertical asociada (None si es 2D)
   catalog: catalogo de variables
   levels: diccionario con las coordenadas verticales y los niveles a extraer (None si se quieren todos los niveles)
   filetype: nombre del modelo a procesar (WRF/GFS)
   """

   self = sys.modules[__name__]

   #Creo el dataset
   ds = getattr(self, 'smn_netcdf_from_{}'.format(filetype.lower()))(filename, variables, catalog, levels)    

   #Estandarizo las coordenadas
   ds = standard_coords(ds, catalog)

   #Estandarizo las variables
   ds = standard_vars(ds, catalog)

   #Estandarizo el dataset
   ds = standard_dataset(ds, catalog)

   return ds


def standard_array(array, catalog_var):
   """
   Funcion que cambia el nombre del array y convierte unidades

   array: DataArray a procesar
   catalog_var: catalogo de la variable que tiene el DataArray
   """

   #Renombro la variable
   array.name = catalog_var['name']

   #Conversion de unidades
   unit_in = units(array.attrs['units']) 
   unit_out = units(catalog_var['attrs']['units'])

   if unit_in != unit_out:
      array = array.metpy.convert_units(unit_out).metpy.dequantify()

   array.attrs['units'] = catalog_var['attrs']['units']

   return array


def standard_coords(ds, catalog):
   """
   Funcion que opera sobre las coordenadas del dataset

   ds: Dataset a procesar
   catalog: catalogo de variables
   """

   if 'MAP_PROJ' in ds.attrs:
      projection = re.sub(' ', '_', ds.MAP_PROJ)
   else:
      projection = None


   for var in ds.coords:
      if var ==  projection:
         continue
      catalog_var = getattr(catalog, var)

      #Elimino atributos 
      ds[var].attrs.pop('_metpy_axis', None)

      #Defino el tipo de dato de la coordenada
      if var not in ['XTIME', 'DATE', projection]:
          ds[var] = ds[var].astype(catalog_var['dtype'])

      #Actualizo los atributos de las coordenadas con el catalogo
      ds[var].attrs.update(catalog_var['attrs'])
      if var in ['XTIME']:
         ds[var].attrs['long_name'] += ' {}'.format(ds.START_DATE)

   return ds


def standard_vars(ds, catalog):
   """
   Funcion que opera sobre las variables del dataset

   ds: Dataset a procesar
   catalog: catalogo de variables
   """

   if 'MAP_PROJ' in ds.attrs:
      projection = re.sub(' ', '_', ds.MAP_PROJ)
   else:
      projection = None

   rm_keys = ['MemoryOrder', 'FieldType', 'stagger', 'projection', 'vert_units', 'missing_value', 'description', '_FillValue', 'coordinates']
   for var in ds.data_vars:

      catalog_var = getattr(catalog, var)

      #Defino el tipo de dato de la variable
      ds[var] = ds[var].astype(catalog_var['dtype'])

      #Elimino atributos 
      for key in ds[var].attrs.keys():
         if (key.split('_')[0] == 'GRIB') and (key not in rm_keys):
            rm_keys.append(key)

      for rm in rm_keys:
         ds[var].attrs.pop(rm, None)

      #Actualizo los atributos de las coordenadas con el catalogo
      ds[var].attrs.update(catalog_var['attrs'])

      if projection in ds.coords:
         ds[var].attrs['grid_mapping'] = projection
         ds[var].attrs['coordinates'] = 'XLAT XLONG'

   return ds


def standard_dataset(ds, catalog):
   """
   Funcion que opera sobre atributos del dataset

   ds: Dataset a procesar
   catalog: catalogo de variables
   """

   if 'MAP_PROJ' in ds.attrs:
      projection = re.sub(' ', '_', ds.MAP_PROJ)
   else:
      projection = None

   #Defino el FILL_VALUE
   fill_value = catalog.FILL_VALUE

   #Defino el encoding de las coordenadas
   for var in ds.coords:
      if var == projection:
         continue
      catalog_var = getattr(catalog, var)
      encoding = {'zlib':'True', 'dtype': catalog_var['dtype'] , '_FillValue': None}
      if var in ['XTIME', 'DATE']:
         encoding['units'] = catalog_var['units'] + ' {}'.format(ds.START_DATE)
      ds[var].encoding = encoding

   #Defino el encoding de las variables
   for var in ds.data_vars:
      chunksize = ()
      for dim in ds[var].dims:
         if dim not in ['x', 'y', 'XLAT', 'XLONG']:
            chunksize += (1,)
         else:
            #chunksize += (ds.dims[dim],)
            chunksize += (ds.sizes[dim],)

      catalog_var = getattr(catalog, var)
      encoding = {'zlib':'True', 'dtype': catalog_var['dtype'], 'chunksizes':chunksize, '_FillValue': fill_value}
      ds[var].encoding = encoding

   return ds


####################
# WRFOUT FUNCTIONS #
####################

def smn_netcdf_from_wrf(filename, variables, catalog, levels):
   """
   Funcion que arma un dataset a partir de un WRFOUT

   filename: nombre del archivo a procesar
   variables: diccionario con el nombre de las variables a procesar y la coordenada vertical asociada (None si es 2D)
   catalog: catalogo de variables
   levels: diccionario con las coordenadas verticales y los niveles a extraer (None si se quieren todos los niveles)
   """

   #Abro el archivo
   ncid = Dataset(filename, 'r')

   #Creo el dataset
   ds = load_wrf(ncid, variables, catalog, levels)

   #Agrego la proyeccion
   proj = ds.attrs['projection'].cf()
   ds = add_projection_wrf(ds, proj)

   #Agrego atributos globales
   ds = add_global_attrs_wrf(ds, ncid)

   #Cierro el archivo
   ncid.close()

   return ds


def load_wrf(ncid, variables, catalog, levels):
   """
   Funcion que crea un dataset con los datos del WRF

   ncid: variable que tiene un wrfout abierto con la libreria netCDF4
   variables: diccionario con el nombre de las variables a procesar y la coordenada vertical asociada (None si es 2D)
   catalog: catalogo de variables
   levels: diccionario con las coordenadas verticales y los niveles a extraer (None si se quieren todos los niveles)
   """

   v_coord = {}
   for coord in levels.keys():
      catalog_var = getattr(catalog, coord)
      #coord_name = catalog_var['wrf']['name']
      #v_coord[coord] = getvar(ncid, coord_name, squeeze = True)
      read_function = eval(catalog_var['wrf']['function'])

      #data = getvar(ncid, catalog_var['wrf']['name'], squeeze = False)
      v_coord[coord] = read_function(ncid, squeeze = False, **catalog_var['wrf']['args'])

   data_vars = []
   for var in variables.keys():
      catalog_var = getattr(catalog, var)

      read_function = eval(catalog_var['wrf']['function'])

      #data = getvar(ncid, catalog_var['wrf']['name'], squeeze = False)
      data = read_function(ncid, squeeze = False, **catalog_var['wrf']['args'])

      if 'Umet' in var:
         data = data[0]
         data = data.drop('u_v')
      if 'Vmet' in var:
         data = data[1]
         data = data.drop('u_v')
      if 'mcape_mcin_lcl_lfc' in data.coords:
         data = data.drop('mcape_mcin_lcl_lfc')

#      stagged_dim = None
#      for i, dim in enumerate(data.dims):
#         if 'stag' in dim:
#            stagged_dim = dim

      coord_name = variables[var]

      if coord_name == 'level_eta':
         data = data.rename({'bottom_top':'level_eta'})
         data = data.assign_coords({'level_eta': v_coord[coord_name].values})

      if variables[var] is not None and levels[coord_name] is not None:
         if coord_name not in ['level_eta']:
#            if stagged_dim is not None:
#               data = destagger(data, data.dims.index(stagged_dim), meta = True)
 
            data = interplevel(data, v_coord[coord_name], levels[coord_name], missing = catalog.FILL_VALUE, squeeze=False)
            if data.vert_units == 'hPa':
               data = data.rename({'level':'level_p'})
            elif data.vert_units == 'm':
               data = data.rename({'level':'level_z'})
            elif data.vert_units == 'K':
               data = data.rename({'level':'level_t'})
            else:
               print('Level not recognized')
         else:
            data = data.sel({coord_name:levels[coord_name]})

      data = standard_array(data, catalog_var)

      data_vars.append(data)

   ds = xr.merge(data_vars)

   ds = ds.drop('XTIME')
   ds = ds.rename({'Time':'XTIME'})

   return ds


def add_projection_wrf(ds, proj):
   """
   Funcion que agrega la variable asociada a la proyeccion en el dataset

   ds: dataset al que se le quieren agregar las variables de la proyeccion
   proj: diccionario con parametros asociados a la proyeccion


   ATENCION: solo probado con Lambert_Conformal
   """

   #Renombro dimensiones
   ds = ds.rename({'south_north': 'y', 'west_east': 'x'})

   # Get projection attrs
   #proj['earth_radius'] = proj.pop('semi_major_axis')

   #Defino la proyeccion
   ds = ds.metpy.assign_crs(proj)

   #Obtengo las coordenadas en la proyeccion
   ds = ds.metpy.assign_y_x(tolerance = 1000 * units.m, force = True)

   #Agrego la proyeccion como coordenada
   proj_name = 'Lambert_Conformal' #re.sub(' ', '_', ds.attrs['MAP_PROJ'])
   ds = ds.assign_coords({proj_name: str()})
   ds[proj_name].attrs = proj

   #Elimino la variable que metpy agrega asociada a la proyeccion
   ds = ds.drop('metpy_crs')

   return ds


def add_global_attrs_wrf(ds, ncid, mem = None):
   """
   Funcion que genera los atributos globales del dataset

   ds: dataset a agregarle atributos globales
   ncid: variable que tiene un wrfout abierto con la libreria netCDF4
   mem = numero de miembro de un ensamble
   """

   fcst_ini = datetime.strptime(ncid.START_DATE, '%Y-%m-%d_%H:%M:%S')
   valid_date = datetime.utcfromtimestamp(ds['XTIME'].values.tolist()[0]/1e9)
   fcst_lead = np.int32((valid_date - fcst_ini).total_seconds()/3600) 

   #Genero el diccionario de atributos
   attrs = {}
   attrs['title'] = 'Python PostProcessing for SMN WRF-ARW'
   attrs['institution'] = 'Servicio Meteorologico Nacional'
   attrs['source'] = ncid.TITLE
   attrs['creation_date'] = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
   if mem is not None:
      attrs['MEMBER'] = mem
   attrs['start_lat'] = np.float32(ds.XLAT.min().data)
   attrs['start_lon'] = np.float32(ds.XLONG.min().data)
   attrs['end_lat'] = np.float32(ds.XLAT.max().data)
   attrs['end_lon'] = np.float32(ds.XLONG.max().data)
   attrs['MAP_PROJ'] = ncid.MAP_PROJ_CHAR
   attrs['STAND_LON'] = ncid.STAND_LON
   attrs['CEN_LON'] = ncid.CEN_LON
   attrs['CEN_LAT'] = ncid.CEN_LAT
   attrs['TRUELAT1'] = ncid.TRUELAT1
   attrs['TRUELAT2'] = ncid.TRUELAT2
   attrs['DX'] = ncid.DX
   attrs['DY'] = ncid.DY
   attrs['START_DATE'] = fcst_ini.strftime("%Y-%m-%d %H:%M:%S")
   attrs['FCST_HOURS'] = fcst_lead
   attrs['VALID_DATE'] = valid_date.strftime("%Y-%m-%d %H:%M:%S")
   attrs['Conventions'] = 'CF-1.8'

   #Seteo los atributos
   ds.attrs = attrs

   return ds


######################
# GFS/GEFS FUNCTIONS #
######################

def smn_netcdf_from_gfs(filename, variables, catalog, levels):
   """
   Funcion que arma un dataset a partir de un grib de GFS/GEFS

   filename: nombre del archivo a procesar
   variables: diccionario con el nombre de las variables a procesar y la coordenada vertical asociada (None si es 2D)
   catalog: catalogo de variables
   levels: diccionario con las coordenadas verticales y los niveles a extraer (None si se quieren todos los niveles)
   """

   #Creo el dataset
   ds = load_gfs(filename, variables, catalog, levels)

   #Si no se encontro ninguna de las variables en el grib, el dataset va a estar vacio. Se omiten los siguientes pasos
   if ds:
      #Estandarizo los nombres
      name_lat = getattr(catalog, 'XLAT')['gfs']['name']
      name_lon = getattr(catalog, 'XLONG')['gfs']['name']
      name_mem = getattr(catalog, 'MEMBER')['gfs']['name']
      name_time = getattr(catalog, 'XTIME')['gfs']['name']
      name_date = getattr(catalog, 'DATE')['gfs']['name']
      name_step = getattr(catalog, 'STEP')['gfs']['name']
      
      ds = ds.rename({name_time: 'XTIME', name_lat:'XLAT', name_lon:'XLONG'})

      if name_mem in ds.coords:
         ds = ds.rename({name_mem:'MEMBER'})

      #Agrego atributos globales
      ds = add_global_attrs_gfs(ds)

      ds = ds.drop_vars([name_date, name_step])

   return ds


def load_gfs(filename, variables, catalog, levels):
   """
   Funcion que crea un dataset con los datos del WRF

   ncid: variable que tiene un wrfout abierto con la libreria netCDF4
   variables: diccionario con el nombre de las variables a procesar y la coordenada vertical asociada (None si es 2D)
   catalog: catalogo de variables
   levels: diccionario con las coordenadas verticales y los niveles a extraer (None si se quieren todos los niveles)
   """

   #cfgrib no soporta leer variables de diferentes tipos de nivel asi que se leen de a uno a la vez
   dict_typeofLevel = {}
   for var in variables.keys():
      catalog_var = getattr(catalog, var)
      typeofLevel = catalog_var['gfs']['level']
      if typeofLevel is not None:
         if typeofLevel in dict_typeofLevel.keys():
            dict_typeofLevel[typeofLevel].append(var)
         else:
            dict_typeofLevel[typeofLevel] = [var]
      else:
         print(f'{var} variable is not in the catalog for GFS/GEFS')

   data_vars = []
   for typeofLevel in dict_typeofLevel.keys():
      vars_level = dict_typeofLevel[typeofLevel]

      for var in vars_level:
         catalog_var = getattr(catalog, var)

         cfVarName = catalog_var['gfs']['name']

         #La variable de precipitacion en los GRIB no tienen definido el atributo cfVarName entonces la busco por 
         #el shortName
         if var == 'PP':
            dict_backend = {'indexpath':'', 'filter_by_keys': {'typeOfLevel':typeofLevel, 'shortName': cfVarName}}
         else:
            dict_backend = {'indexpath':'', 'filter_by_keys': {'typeOfLevel':typeofLevel, 'cfVarName': cfVarName}}

         try:
            with xr.open_dataset(filename, engine = 'cfgrib', backend_kwargs = dict_backend) as ds:

               if cfVarName in ds.data_vars:

                  ds = ds.rename({cfVarName:catalog_var['name']})

                  if var == 'PP':
                     #La unidad de precipitacion en GFS/GEFS es kg/m2 y no se puede convertir automaticamente a mm.
                     #Se la convierte dividiendo por la densidad del agua
                     ds[var] = (ds[var].metpy.quantify() / constants.rho_l).metpy.dequantify()
                  elif var == 'GEOPT':
                     ds[var] = (ds[var].metpy.quantify() * constants.g).metpy.dequantify()
                  elif var == 'LANDMASK':
                     ds[var].attrs['units'] = ''

                  data = standard_array(ds[var], catalog_var)

                  #Selecciono los niveles solicitados si la variable es 3D
                  if data[typeofLevel].size > 1:
                     if variables[var] is not None:
                        coord_name = variables[var]
                        if coord_name != 'PRESSURE':
                           print(f'Levels {levels[coord_name]} could not be selected for {var}. GFS/GEFS only accept selection of pressure levels')
                           continue
                        else:
                           data = data.sel({typeofLevel:levels[coord_name]})
                     vert_units = data[typeofLevel].attrs['units']
                     if vert_units == 'hPa':
                        level_rename = 'level_p'
                     else:
                        print('Level not recognized')
                     data = data.rename({typeofLevel:level_rename}) 
                  else:
                     data = data.drop_vars(typeofLevel)

                  data_vars.append(data)
               else:
                  print(f'Variable {var} is not in {filename}')
         except EOFError:
            print(f'{filename} is broken, it is ommited')
         except FileNotFoundError:
            print(f"{filename} doesn't exist")
            

   return xr.merge(data_vars, combine_attrs = 'drop_conflicts')


def add_global_attrs_gfs(ds):
   """
   Funcion que genera los atributos globales del dataset

   ds: dataset a agregarle atributos globales
   """

   if 'MEMBER' in ds.coords:
      modelo = 'GEFS'
      miembro = np.int32(ds['MEMBER'])
   else:
      modelo = 'GFS'
      miembro = np.int32(0)

   #Genero el diccionario de atributos
   attrs = {}
   attrs['title'] = f'Python PostProcessing for {modelo}'
   attrs['institution'] = 'Servicio Meteorologico Nacional'
   attrs['source'] = modelo
   attrs['creation_date'] = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
   attrs['MEMBER'] = miembro
   attrs['start_lat'] = float(ds['XLAT'].min().data)
   attrs['start_lon'] = float(ds['XLONG'].min().data)
   attrs['end_lat'] = float(ds['XLAT'].max().data)
   attrs['end_lon'] = float(ds['XLONG'].max().data)
   attrs['MAP_PROJ'] = ds.attrs['GRIB_gridType']
   attrs['DX'] = ds.attrs['GRIB_jDirectionIncrementInDegrees']
   attrs['DY'] = ds.attrs['GRIB_iDirectionIncrementInDegrees']
   attrs['START_DATE'] = datetime.utcfromtimestamp(ds['time'].values.tolist()/1e9).strftime("%Y-%m-%d %H:%M:%S")
   attrs['FCST_HOURS'] = np.int32(ds['step'].values.astype('timedelta64[h]')/np.timedelta64(1, 'h'))
   attrs['VALID_DATE'] = datetime.utcfromtimestamp(ds['XTIME'].values.tolist()/1e9).strftime("%Y-%m-%d %H:%M:%S")
   attrs['Conventions'] = 'CF-1.8'

   #Seteo los atributos
   ds.attrs = attrs

   return ds


####################################################
### Funciones propias de generacion de variables ###
####################################################

def get_gust10(nc, squeeze):

    if 'TKE_PBL' in nc.variables:
        # Para la parametrización MYJ que calcula TKE  (Kurbatova et al. 2018)
        TKE = wrf.getvar(nc, 'TKE_PBL', squeeze = squeeze).metpy.quantify()
        TKE0 = TKE.isel(bottom_top_stag = 0) #Tomo la TKE del nivel mas bajo
        WindSpeed10 = wrf.getvar(nc, 'wspd10', squeeze = squeeze).metpy.quantify()
        WindSpeed10 = WindSpeed10.drop('wspd_wdir')
        Gust = (WindSpeed10 + 3*np.sqrt(TKE0)).metpy.dequantify()
#    elif 'UST' in nc.variables:
    else:
        # Para la parametrizacion YSU que calcula no TKE pero si la velocidad de friccion (Born et al. 2012)
        u_fric = wrf.getvar(nc, 'UST', squeeze = squeeze)
        WindSpeed = wrf.getvar(nc, 'wspd', squeeze = squeeze).drop('wspd_wdir')
        Z = wrf.getvar(nc, 'height_agl', squeeze = squeeze)
        WindSpeed30 = interplevel(WindSpeed, Z, 30, squeeze=True)
        WindSpeed30 = WindSpeed30.drop_vars('level')
        Gust = np.abs(WindSpeed30) + 3.0*2.4*u_fric
#    else:
#        print('No se pueden calcular las rafagas')

    return Gust


#####################
# Lossy compression #
#####################

def TrimPrecision(a, keepbits):
    """
    Funcion que realiza el redondeo half-to-even para poner 0's en los ultimos bits de la mantisa.
    Codigo extraido de "A note on precision-preserving compression of scientific data", Kouznetsov 2021
    https://gmd.copernicus.org/articles/14/377/2021/

    a: array con los datos a comprimir
    keepbits: cantidad de bits de la mantisa a conservar
    """

    assert (a.dtype == np.float32)
    b = a.view(dtype=np.int32)
    maskbits = 23 - keepbits
    mask = (0xFFFFFFFF>>maskbits)<<maskbits
    half_quantum1 = ( 1<<(maskbits-1) ) - 1
    b += ((b>>maskbits) & 1) +half_quantum1
    b &= mask


def lossy_compression(ds, keepbits, variables = []):
    """
    Funcion para aplicar realizar compresion con perdida de informacion en un dataset

    ds: dataset que se quiere comprimir
    keepbits: cantidad de bits de la mantisa a conservar
    variables: lista de variables a comprimir, por default se comprimen todas
    """

    if not variables:
        variables = ds.data_vars

    for var in variables:
        if var in ds.coords:
            print(f"Variable {var} is a coordinate, should not be compressed. It's ommited")
            continue

        TrimPrecision(ds[var].values, keepbits)
        ds[var].attrs['QuantizeBitRoundNumberOfSignificantBits'] = np.int32(keepbits)
                                                                                         
    return ds
 
