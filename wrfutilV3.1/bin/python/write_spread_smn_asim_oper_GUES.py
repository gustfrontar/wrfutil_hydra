# -*- coding: utf-8 -*-

# ------------------------------------------
# Genera un netcdf con la estructura de un wrfout
# con el spread al cuadrado, es decir 1/(R-1) sum (Xme-xi)²
# de manera de poder cosntruir: 
#	campos de spread, calculando la raiz cuadrada de esta cantidad 
#	valores de spread que resuman todo el campo, o todos los tiempos
#
# (siguiendo Fortin 2104)
#
# Activar el entorno:
# source activate wrfutil
#
# Ejecutar con:
#
# 
# Tambien se necesita tener en el entorno cargadas las siguientes variables:
# PATHOUT      - Directorio donde se escribiran las salidas
# GUESFILES     - Path absoluto de los analisis de todos los miembros
#		ej: export GUESFILES=$(ls /data/yanina/wrfutilV3.1/RUNs/HIST/POST/20200102_060000/prueba/*GUE* )
# MIEMBROS     - cantidad de miembros del ensamble
# ------------------------------------------

import numpy as np
import os
from netCDF4 import Dataset
from datetime import datetime
from datetime import timedelta
import pyproj

inicio = datetime.now()
print(str(datetime.now()) + " Empezamos")

path_out = os.environ["PATHOUT"]
R = os.environ["MIEMBROS"]
R = int(R)
fcst_gues = os.environ["ANALISIS"]      # en minutos -> para completar la fecha del pronostico del first guess

###################
# Defino la descripcion del archivo para escribir las salidas
###################
description_short = 'GSP'
description_long = 'Square of First Guess Spread [ 1/(R-1) sum (Xme-Xi)2 ], i={1...R} ens members'

####################
#--- Abro los wrfout
#    generando un diccionario con clave:
#	Miembro, valor:(archivo abierto,nombre archivo, fecha,miembro)
####################

guesfiles = os.environ['GUESFILES']
guesfiles = guesfiles.split("\n")
netfiles = {}
for afile in guesfiles:
      miem = int(afile.split('.')[-2].split("OPER")[-1])
      netfiles[miem] = (Dataset(afile,'r'), afile, afile.split('.')[-4], miem)

# para obetenr la info de la fecha, leo el miembro 1 (que ya está abierto):
fecha_ini = datetime.strptime(netfiles[1][0].START_DATE, '%Y-%m-%d_%H:%M:%S')
fecha_fc = fecha_ini + timedelta(minutes=int(fcst_gues))
Yi = str(fecha_ini.year)
Mi = str(fecha_ini.month).zfill(2)
Di = str(fecha_ini.day).zfill(2)
Hi = str(fecha_ini.hour).zfill(2)
hf_out = str(int((fecha_fc - fecha_ini).total_seconds()/3600)).zfill(3)

####################
#--- Creo el primer netcdf de salida, que va a contener el spread al cuadrado de las variables en altura:
####################

print (str(datetime.now())+" Creamos el archivo de salida 1, LEV")

tipo = 'LEV'
post_dir = path_out
post_name = 'model.WRF_' + description_short + '_4km.' + Yi + Mi + Di + '_' + Hi + '0000.' + hf_out + '.OPER' + tipo +'.nc'
wrf_DPA = Dataset(post_dir + post_name, 'w')

#--- Dimensiones:

XLAT = netfiles[1][0].variables["XLAT"][:]   # Latitud
XLONG = netfiles[1][0].variables["XLONG"][:]  # Longitud
LEV  = netfiles[1][0].variables["lev"][:]    # Niveles
XTIME = netfiles[1][0].variables["XTIME"]

nlev = np.shape(LEV)[0]
nlat, nlon = XLAT.shape

wrf_DPA.createDimension("time", None)
wrf_DPA.createDimension("lev", len(LEV))
wrf_DPA.createDimension("y", nlat)
wrf_DPA.createDimension("x", nlon)

#--- Atributos globales:
wrf_ana = netfiles[1][0]

atributos_globales = {}
atributos_globales['title'] = 'Python PostProcessing for SMN WRF-ARW ' + description_long + tipo
atributos_globales["institution"] = "Servicio Meteorologico Nacional"
atributos_globales["source"] = wrf_ana.source
atributos_globales['creation_date'] = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
atributos_globales['MAP_PROJ'] = wrf_ana.MAP_PROJ
atributos_globales['MAP_PROJ_CHAR'] = wrf_ana.MAP_PROJ_CHAR
atributos_globales['STAND_LON'] = wrf_ana.STAND_LON
atributos_globales['CEN_LON'] = wrf_ana.CEN_LON
atributos_globales['CEN_LAT'] = wrf_ana.CEN_LAT
atributos_globales['TRUELAT1'] = wrf_ana.TRUELAT1
atributos_globales['TRUELAT2'] = wrf_ana.TRUELAT2
atributos_globales['DX'] = wrf_ana.DX
atributos_globales['DY'] = wrf_ana.DY
atributos_globales['Z_LEVELS'] = nlev
atributos_globales['START_DATE'] = wrf_ana.START_DATE 
atributos_globales['FCST_HOURS'] = hf_out
atributos_globales['VALID_DATE'] = fecha_fc.strftime("%Y-%m-%d_%H:%M:%S")
atributos_globales["Conventions"] = "CF-1.8"

wrf_DPA.setncatts(atributos_globales)

#--- Variable que define la proyeccion

lcc = pyproj.Proj(proj="lcc", lat_1 = wrf_ana.TRUELAT1, lat_2 = wrf_ana.TRUELAT2, # Cone intersects with the sphere
                  lat_0 = wrf_ana.CEN_LAT, lon_0 = wrf_ana.CEN_LON, a=6370000, b=6370000)


#Generacion de x e y tomado de https://fabienmaussion.info/2018/01/06/wrf-projection/

wgs_proj = pyproj.Proj(proj="latlong", datum="WGS84")
e, n = pyproj.transform(wgs_proj, lcc, wrf_ana.CEN_LON, wrf_ana.CEN_LAT)

dx, dy = wrf_ana.DX, wrf_ana.DY
nx, ny = nlon, nlat

x0 = -(nx-1) / 2. * dx + e
y0 = -(ny-1) / 2. * dy + n

x = (np.arange(nx) * dx + x0)
y = (np.arange(ny) * dy + y0)

outvar = wrf_DPA.createVariable("Lambert_Conformal", "c", ())
outvar.grid_mapping_name = "lambert_conformal_conic"
outvar.standard_parallel = [wrf_ana.TRUELAT1, wrf_ana.TRUELAT2]
outvar.latitude_of_projection_origin = wrf_ana.CEN_LAT
outvar.longitude_of_central_meridian = wrf_ana.CEN_LON

#--- Variables asociadas a las dimensiones

outvar = wrf_DPA.createVariable("XTIME", np.float32, ("time"), zlib = True)
outvar[0] = XTIME[:]/60
outvar.long_name = fecha_ini.strftime("hours since %Y-%m-%d %H:%M:%S")
outvar.units = fecha_ini.strftime("hours since %Y-%m-%d %H:%M:%S")
outvar.standard_name = "time"
outvar.axis = "T"
outvar.calendar = "gregorian"

outvar = wrf_DPA.createVariable("lev", np.float32, ("lev"), zlib = True)
outvar.long_name = "Pressure"
outvar.units = "hPa"
outvar.postive = "down"
outvar.axis = "Z"
outvar.standard_name = "air_pressure"
outvar[:] = LEV[:]

outvar = wrf_DPA.createVariable("x", np.float32, ("x"), zlib = True)
outvar[:] = x
outvar.axis = "X"
outvar.units = "m"
outvar.standard_name = "projection_x_coordinate"
outvar.long_name = "x-coordinate in projected coordinate system"

outvar = wrf_DPA.createVariable("y", np.float32, ("y"), zlib = True)
outvar[:] = y
outvar.axis = "Y"
outvar.units = "m"
outvar.standard_name = "projection_y_coordinate"
outvar.long_name = "y-coordinate in projected coordinate system"

####
# Calculamos y Escribimos el cuadrado del spread de LEV:
####

varList=[("Q",       "Q_spread2",       "Specific Humidity (spread)2", "specific_humidity"  ,"(g Kg-1)2") ,
         ("T",       "T_spread2",       "Temperature (spread)2" ,"air_temperature"  ,"K2"),
         ("Umet",    "U_spread2",       "Zonal Wind Component (spread)2"    ,"eastward_wind"    ,"(m s-1)2"),
         ("Vmet",    "V_spread2",       "Meridional Wind Component (spread)2"   ,"northward_wind"   ,"(m s-1)2"),
         ("GEOPT",   "GEOPT_spread2",   "Geopotential (spread)2"    ,"geopotential_height"  ,"(m2 s-2)2")]

#Latitud
outvar = wrf_DPA.createVariable("XLAT", np.float32, ("y", "x"), zlib = True)
outvar[:, :] = XLAT[:, :]
outvar.standard_name = "latitude"
outvar.long_name = "Latitude"
outvar.units = "degrees_north"

#Longitud
outvar = wrf_DPA.createVariable("XLONG", np.float32, ("y", "x"), zlib = True)
outvar[:, :] = XLONG[:, :]
outvar.standard_name = "longitude"
outvar.long_name = "Longitude"
outvar.units = "degrees_east"


for (var_ori, alias, lname, stdname, units) in varList:
    print(str(datetime.now()) + " Calculando y Guardando -> " + str(alias))
    x = np.full((R, nlev, nlat, nlon), np.nan, np.float32)
    dif = np.full((R, nlev, nlat, nlon), np.nan, np.float32)
    for miem in range(1, 41):
        x[miem-1, :, :, :] = netfiles[miem][0].variables[var_ori][:]
    xme = np.nanmean(x, axis=0)
#   print(np.nanmean(xme))	#debug
    for miem in range(1, 41):
        dif[miem-1, :, :, :] = xme-x[miem-1, :, :, :]

    outvar = wrf_DPA.createVariable(alias, np.float32, ("time", "lev", "y", "x"), zlib = True)
    outvar[0, :, :, :] = (np.ma.sum((dif*dif), axis=0))/(R-1)
#   print(np.nanmean(outvar))	#debug
    outvar.standard_name = stdname
    outvar.long_name = lname
    outvar.units = units
    outvar.coordinates = "XLAT XLONG"
    outvar.grid_mapping = "Lambert_Conformal"

print(str(datetime.now()) + " Cerrando archivo de salida")
wrf_DPA.close()

escritura1 = datetime.now()
dif = escritura1 - inicio
print('Para generar ' + post_name + ' tardamos ' + str(dif.seconds) + ' segundos')

####################
#--- Creo el segundo netcdf de salida:
####################

print(str(datetime.now()) + " Creamos el archivo de salida 2, SFC")

tipo = 'SFC'
post_dir = path_out
post_name = 'model.WRF_' + description_short + '_4km.' + Yi + Mi + Di + '_' + Hi + '0000.' + hf_out + '.OPER' + tipo +'.nc'
wrf_DPA = Dataset(post_dir + post_name, 'w')

#--- Dimensiones:
# Ya definimos antes las variables necesarias

wrf_DPA.createDimension("time", None)
wrf_DPA.createDimension("y", nlat)
wrf_DPA.createDimension("x", nlon)

#--- Atributos globales:
atributos_globales = {}
atributos_globales['title'] = 'Python PostProcessing for SMN WRF-ARW ' + description_long + tipo
atributos_globales["institution"] = "Servicio Meteorologico Nacional"
atributos_globales["source"] = wrf_ana.source
atributos_globales['creation_date'] = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
atributos_globales['MAP_PROJ'] = wrf_ana.MAP_PROJ
atributos_globales['MAP_PROJ_CHAR'] = wrf_ana.MAP_PROJ_CHAR
atributos_globales['STAND_LON'] = wrf_ana.STAND_LON
atributos_globales['CEN_LON'] = wrf_ana.CEN_LON
atributos_globales['TRUELAT1'] = wrf_ana.TRUELAT1
atributos_globales['TRUELAT2'] = wrf_ana.TRUELAT2
atributos_globales['DX'] = wrf_ana.DX
atributos_globales['DY'] = wrf_ana.DY
atributos_globales['Z_LEVELS'] = nlev
atributos_globales['START_DATE'] = wrf_ana.START_DATE
atributos_globales['FCST_HOURS'] = hf_out
atributos_globales['VALID_DATE'] = fecha_fc.strftime("%Y-%m-%d_%H:%M:%S")
atributos_globales["Conventions"] = "CF-1.8"

wrf_DPA.setncatts(atributos_globales)

#--- Variable que define la proyeccion

lcc = pyproj.Proj(proj="lcc", lat_1 = wrf_ana.TRUELAT1, lat_2 = wrf_ana.TRUELAT2, # Cone intersects with the sphere
                  lat_0 = wrf_ana.CEN_LAT, lon_0 = wrf_ana.CEN_LON, a=6370000, b=6370000)


#Generacion de x e y tomado de https://fabienmaussion.info/2018/01/06/wrf-projection/

wgs_proj = pyproj.Proj(proj="latlong", datum="WGS84")
e, n = pyproj.transform(wgs_proj, lcc, wrf_ana.CEN_LON, wrf_ana.CEN_LAT)

dx, dy = wrf_ana.DX, wrf_ana.DY
nx, ny = nlon, nlat

x0 = -(nx-1) / 2. * dx + e
y0 = -(ny-1) / 2. * dy + n

x = (np.arange(nx) * dx + x0)
y = (np.arange(ny) * dy + y0)

outvar = wrf_DPA.createVariable("Lambert_Conformal", "c", ())
outvar.grid_mapping_name = "lambert_conformal_conic"
outvar.standard_parallel = [wrf_ana.TRUELAT1, wrf_ana.TRUELAT2]
outvar.latitude_of_projection_origin = wrf_ana.TRUELAT1
outvar.longitude_of_central_meridian = wrf_ana.CEN_LON

#--- Variables asociadas a las dimensiones

outvar = wrf_DPA.createVariable("XTIME", np.float32, ("time"), zlib = True)
outvar[0] = XTIME[:]/60
outvar.long_name = fecha_ini.strftime("hours since %Y-%m-%d %H:%M:%S")
outvar.units = fecha_ini.strftime("hours since %Y-%m-%d %H:%M:%S")
outvar.standard_name = "time"
outvar.axis = "T"
outvar.calendar = "gregorian"

outvar = wrf_DPA.createVariable("x", np.float32, ("x"), zlib = True)
outvar[:] = x
outvar.axis = "X"
outvar.units = "m"
outvar.standard_name = "projection_x_coordinate"
outvar.long_name = "x-coordinate in projected coordinate system"

outvar = wrf_DPA.createVariable("y", np.float32, ("y"), zlib = True)
outvar[:] = y
outvar.axis = "Y"
outvar.units = "m"
outvar.standard_name = "projection_y_coordinate"
outvar.long_name = "y-coordinate in projected coordinate system"

####
# Calculamos y Escribimos el cuadrado del spread de SFC:
####

varList=[("PSFC",       "PSFC_spread2",       "Surface Pressure (spread)2", "surface_air_pressure"  ,"(hPa)2") ,
         ("T2",         "T2_spread2",         "2-m Temperature (spread)2",  "air_temperature"   ,"K2"),
         ("Q2",         "Q2_spread2",         "2-m Water Vapor mixing ratio (spread)2", 'specific_humidity' ,"(g Kg-1)2"),
         ("PP",         "PP_spread2",         "Accumulated total precipitation (spread)2",  'precipitation_amount'  ,"(mm)2"),
         ("MDBZ",       "MDBZ_spread2",       "Max reflectivity (spread)2", 'equivalent_reflectivity_factor',   "(dBZ)2"),
         ("SLP",        "SLP_spread2",        "Sea Level Pressure (spread)2",   'air_pressure_at_mean_sea_level'    ,"(hPa)2"),
         ("Umet10",     "U10_spread2",        "10-m U wind component (spread)2",    'eastward_wind' ,"(m s-1)2"),
         ("Vmet10",     "V10_spread2",        "10-m V wind component (spread)2",    'northward_wind'    ,"(m s-1)2"),
         ("REFL1KM",    "REF1KM_spread2",     "Reflectivity at 1km agl (spread)2",  'equivalent_reflectivity_factor',   "(dBZ)2"),
         ("REFL4KM",    "REF4KM_spread2",     "Reflectivity at 4km agl (spread)2",  'equivalent_reflectivity_factor',   "(dBZ)2"),]

#Latitud
outvar = wrf_DPA.createVariable("XLAT", np.float32, ("y", "x"), zlib = True)
outvar[:, :] = XLAT[:, :]
outvar.standard_name = "latitude"
outvar.long_name = "Latitude"
outvar.units = "degrees_north"

#Longitud
outvar = wrf_DPA.createVariable("XLONG", np.float32, ("y", "x"), zlib = True)
outvar[:, :] = XLONG[:, :]
outvar.standard_name = "longitude"
outvar.long_name = "Longitude"
outvar.units = "degrees_east"


for (var_ori, alias, lname, stdname, units) in varList:
    print(str(datetime.now()) + " Calculando y Guardando -> " + str(alias))
    x = np.full((R, nlat, nlon), np.nan, np.float32)
    dif = np.full((R, nlat, nlon), np.nan, np.float32)
    for miem in range(1, 41):
        x[miem-1, :, :] = netfiles[miem][0].variables[var_ori][:]
    xme = np.nanmean(x, axis=0)
    for miem in range(1, 41):
        dif[miem-1, :, :] = xme-x[miem-1, :, :]

    outvar = wrf_DPA.createVariable(alias, np.float32, ("time", "y", "x"), zlib = True)
    outvar[0, :, :] = (np.ma.sum((dif*dif), axis=0))/(R-1)
    outvar.standard_name = stdname
    outvar.long_name = lname
    outvar.units = units
    outvar.coordinates = "XLAT XLONG"
    outvar.grid_mapping = "Lambert_Conformal"

print(str(datetime.now()) + " Cerrando archivo de salida")
wrf_DPA.close()


escritura2 = datetime.now()
dif = escritura2 - escritura1
print('Para generar ' + post_name + ' tardamos ' + str(dif.seconds) + ' segundos')

