
# ------------------------------------------
# Genera un netcdf con la estructura de un wrfout
# con las siguientes variables:
#	---
# para un miembro del ensamble
# para un dado tiempo
#
# Activar el entorno:
# source activate NPP
#
# Ejecutar con:
#	python write_post_smn_ens_1mem_sfc.py 
#
# Tambien se necesita tener en el entorno cargadas las siguientes variables:
# PATHOUT - Directorio donde se escribiran las salidas
# MIER    - Miembro del ensamble que se esta procesando
# FILEIN  - Path absoluto del archivo que se esta procesando
# 
# Se le agrega PW y rafagas en superficie
# ------------------------------------------

import numpy as np
import os
from wrf import getvar, interplevel
from netCDF4 import Dataset
from datetime import datetime
from datetime import timedelta
import pyproj
from funciones_post import add_fill_value

fill_value = 1e20

inicio = datetime.now()

path_out = os.environ["PATHOUT"]
mem = os.environ["MIEM"]
filename = os.environ["FILEIN"]


####################
#--- Abro el wrfout:
####################

wrfout = Dataset(filename, 'r')


fecha_ini = datetime.strptime(wrfout.START_DATE, '%Y-%m-%d_%H:%M:%S')
fecha_fc = datetime.strptime(filename[-19:], '%Y-%m-%d_%H_%M_%S')
Yi = str(fecha_ini.year)
Mi = str(fecha_ini.month).zfill(2)
Di = str(fecha_ini.day).zfill(2)
Hi = str(fecha_ini.hour).zfill(2)
hf_out = str(int((fecha_fc - fecha_ini).total_seconds()/3600)).zfill(3)
member = str(mem).zfill(2)


####################
#--- Defino los niveles verticales de salida:
####################

LEV = [700., 500.]
AGL = [1000, 4000]

####################
#--- Obtengo variables directas del wrfout:
####################

XLAT = getvar(wrfout, "XLAT") 	# Latitud
XLONG = getvar(wrfout, "XLONG") 	# Longitud
XTIME = wrfout.variables["XTIME"]
nlat, nlon = XLAT.shape

####################
#--- Creo el netcdf de salida:
####################

tipo = 'SFC'
post_dir = path_out
post_name = 'model.WRF_ENS_4km.' + Yi + Mi + Di + '_' + Hi + '0000.' + hf_out + '.M' + member + "_OPER" + tipo +'.nc'
wrf_DPA = Dataset(post_dir + post_name, 'w')

#--- Dimensiones:

wrf_DPA.createDimension("time", None)
wrf_DPA.createDimension("lev", len(LEV))
wrf_DPA.createDimension("y", nlat)
wrf_DPA.createDimension("x", nlon)

#--- Atributos globales:

atributos_globales = {}
atributos_globales['title'] = 'Python PostProcessing for SMN WRF-ARW Ensemble ' + tipo
atributos_globales["institution"] = "Servicio Meteorologico Nacional"
atributos_globales["source"] = wrfout.TITLE
atributos_globales['creation_date'] = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
atributos_globales['start_lat'] = XLAT[0,0]
atributos_globales['start_lon'] = XLONG[0,0]
atributos_globales['end_lat'] = XLAT[nlat-1,nlon-1]
atributos_globales['end_lon'] = XLONG[nlat-1,nlon-1]
atributos_globales['MAP_PROJ'] = wrfout.MAP_PROJ
atributos_globales['MAP_PROJ_CHAR'] = wrfout.MAP_PROJ_CHAR
atributos_globales['STAND_LON'] = wrfout.STAND_LON
atributos_globales['CEN_LON'] = wrfout.CEN_LON
atributos_globales['CEN_LAT'] = wrfout.CEN_LAT
atributos_globales['TRUELAT1'] = wrfout.TRUELAT1
atributos_globales['TRUELAT2'] = wrfout.TRUELAT2
atributos_globales['DX'] = wrfout.DX
atributos_globales['DY'] = wrfout.DY
atributos_globales['START_DATE'] = wrfout.START_DATE 
atributos_globales['FCST_HOURS'] = hf_out
atributos_globales['VALID_DATE'] = fecha_fc.strftime("%Y-%m-%d_%H:%M:%S")
atributos_globales['ENS_MEMBER'] = member
atributos_globales['Conventions'] = "CF-1.8"

wrf_DPA.setncatts(atributos_globales)

#--- Variable que define la proyeccion

lcc = pyproj.Proj(proj="lcc", lat_1 = wrfout.TRUELAT1, lat_2 = wrfout.TRUELAT2, # Cone intersects with the sphere
                  lat_0 = wrfout.CEN_LAT, lon_0 = wrfout.CEN_LON, a=6370000, b=6370000)


#Generacion de x e y tomado de https://fabienmaussion.info/2018/01/06/wrf-projection/

wgs_proj = pyproj.Proj(proj="latlong", datum="WGS84")
e, n = pyproj.transform(wgs_proj, lcc, wrfout.CEN_LON, wrfout.CEN_LAT)

dx, dy = wrfout.DX, wrfout.DY
nx, ny = nlon, nlat

x0 = -(nx-1) / 2. * dx + e
y0 = -(ny-1) / 2. * dy + n

x = (np.arange(nx) * dx + x0)
y = (np.arange(ny) * dy + y0)

outvar = wrf_DPA.createVariable("Lambert_Conformal", "c", (), fill_value = fill_value)
outvar.grid_mapping_name = "lambert_conformal_conic"
outvar.standard_parallel = [wrfout.TRUELAT1, wrfout.TRUELAT2]
outvar.latitude_of_projection_origin = wrfout.CEN_LAT
outvar.longitude_of_central_meridian = wrfout.CEN_LON

#--- Atributos y Escritura de cada variable:

outvar = wrf_DPA.createVariable("XTIME", np.float32, ("time"), zlib = True, fill_value = fill_value)
outvar[0] = add_fill_value(XTIME[:]/60, fill_value)
outvar.long_name = fecha_ini.strftime("hours since %Y-%m-%d %H:%M:%S")
outvar.units = fecha_ini.strftime("hours since %Y-%m-%d %H:%M:%S")
outvar.standard_name = "time"
outvar.axis = "T"
outvar.calendar = "gregorian"

outvar = wrf_DPA.createVariable("lev", np.float32, ("lev"), zlib = True, fill_value = fill_value)
outvar[:] = add_fill_value(LEV[:], fill_value)
outvar.long_name = "Pressure"
outvar.units = "hPa"
outvar.postive = "down"
outvar.axis = "Z"
outvar.standard_name = "air_pressure"

outvar = wrf_DPA.createVariable("x", np.float32, ("x"), zlib = True, fill_value = fill_value)
outvar[:] = add_fill_value(x, fill_value)
outvar.axis = "X"
outvar.units = "m"
outvar.standard_name = "projection_x_coordinate"
outvar.long_name = "x-coordinate in projected coordinate system"

outvar = wrf_DPA.createVariable("y", np.float32, ("y"), zlib = True, fill_value = fill_value)
outvar[:] = add_fill_value(y, fill_value)
outvar.axis = "Y"
outvar.units = "m"
outvar.standard_name = "projection_y_coordinate"
outvar.long_name = "y-coordinate in projected coordinate system"



#####A partir de aca las Variables

#Latitud
outvar = wrf_DPA.createVariable("XLAT", np.float32, ("y", "x"), zlib = True, fill_value = fill_value)
outvar[:, :] = add_fill_value(XLAT[:, :], fill_value)
outvar.standard_name = "latitude"
outvar.long_name = "Latitude"
outvar.units = "degrees_north"

#Longitud
outvar = wrf_DPA.createVariable("XLONG", np.float32, ("y", "x"), zlib = True, fill_value = fill_value)
outvar[:, :] = add_fill_value(XLONG[:, :], fill_value)
outvar.standard_name = "longitude"
outvar.long_name = "Longitude"
outvar.units = "degrees_east"

# Presion de superficie [hPa]
data = getvar(wrfout, "PSFC")*1e-2
outvar = wrf_DPA.createVariable("PSFC", np.float32, ("time", "y", "x"),zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
outvar.standard_name = "surface_air_pressure"
outvar.long_name = "Surface Pressure"
outvar.units = "hPa"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

# Presion a nivel del mar [hPa]
data = getvar(wrfout, "slp")
outvar = wrf_DPA.createVariable('SLP', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
outvar.standard_name = 'air_pressure_at_mean_sea_level'
outvar.long_name = "Sea Level Pressure"
outvar.units = "hPa"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"


# T 2m [K]
data = getvar(wrfout, "T2")
outvar = wrf_DPA.createVariable("T2", np.float32, ("time", "y", "x"), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
outvar.standard_name = "air_temperature"
outvar.long_name = "2-m Temperature"
outvar.units = "K"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

# Q 2m [g kg-1]
data = getvar(wrfout, "Q2")*1e3
outvar = wrf_DPA.createVariable('Q2', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
outvar.standard_name = 'specific_humidity'
outvar.long_name = "2-m Water Vapor mixing ratio"
outvar.units = "g kg-1"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

# Precipitacion acumulada [mm]
data = getvar(wrfout, 'RAINNC')
outvar = wrf_DPA.createVariable('PP', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
outvar.standard_name = 'precipitation_amount'
outvar.long_name = "Accumulated total precipitation"
outvar.units = "mm"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

# Nieve y hielo acumulado [mm]
data = getvar(wrfout, "SNOWNC")
outvar = wrf_DPA.createVariable('SNOWNC', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
outvar.standard_name = 'sea_ice_and_surface_snow_amount'
outvar.long_name = "Accumulated total snow and ice"
outvar.units = "mm"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

#cape_2d tine CAPE, CIN, LCL y LFC
data = getvar(wrfout, "cape_2d")

# CAPE [J Kg-1]
outvar = wrf_DPA.createVariable('MCAPE', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(np.where(~np.isnan(data[0, :, :]),data[0, :, :],0), fill_value)
outvar.standard_name = 'atmosphere_convective_available_potential_energy'
outvar.long_name = "Maximum CAPE"
outvar.units = "J Kg-1"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

# CIN [J Kg-1]
outvar = wrf_DPA.createVariable('CIN', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(data[1, :, :], fill_value)
outvar.standard_name = 'atmosphere_convective_inhibition'
outvar.long_name = "Maximum Convective Inhibition"
outvar.units = "J Kg-1"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"


# Agua precipitable [kg m-2]
data = getvar(wrfout, "pw")
outvar = wrf_DPA.createVariable('PW', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
outvar.standard_name = 'atmosphere_mass_content_of_water_vapor'
outvar.long_name = "Precipitable Water"
outvar.units = "kg m-2"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

# Reflectividad 3D --> se necesita para interpolar a 1 y 4 km
DBZ = getvar(wrfout, "dbz")
hgt = getvar(wrfout, "height_agl", units = 'm')
DBZ_int = np.full((2, nlat, nlon), np.nan)
for i, agl in enumerate(AGL):
    DBZ_int[i, :, :] = interplevel(DBZ, hgt, agl)

# Reflectividad 1km [dBZ]
outvar = wrf_DPA.createVariable('REFL1KM', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(DBZ_int[0, :, :], fill_value)
outvar.standard_name = 'equivalent_reflectivity_factor'
outvar.long_name = "Reflectivity at 1km agl"
outvar.units = "dBZ"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

# Reflectividad 4km 
outvar = wrf_DPA.createVariable('REFL4KM', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(DBZ_int[1, :, :], fill_value)
outvar.standard_name = 'equivalent_reflectivity_factor'
outvar.long_name = "Reflectivity at 4km agl"
outvar.units = "dBZ"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

# Max reflectividad
data = getvar(wrfout, "mdbz")
outvar = wrf_DPA.createVariable('MDBZ', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
outvar.standard_name = 'equivalent_reflectivity_factor'
outvar.long_name = "Max Reflectivity"
outvar.units = "dBZ"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"


### Aca los vientos#

# Viento Zonal a 10m [m s-1]
viento = getvar(wrfout, "uvmet10")
outvar = wrf_DPA.createVariable('Umet10', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(viento[0][:, :], fill_value)
outvar.standard_name = 'eastward_wind'
outvar.long_name = "10-m U wind component"
outvar.units = "m s-1"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

#Viento meridional a 10m [m s-1]
outvar = wrf_DPA.createVariable('Vmet10', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0, :, :] = add_fill_value(viento[1][:, :], fill_value)
outvar.standard_name = 'northward_wind'
outvar.long_name = "10-m V wind component"
outvar.units = "m s-1"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"


## Rafagas
try:
        ### para la parametrización MYJ qu calcula TKE  (Kurbatova et al. 2018)
        data = getvar(wrfout, 'TKE_PBL')
        data = np.sqrt(viento[0]**2 + viento[1]**2) + 3*np.sqrt(data[0,:,:])

        outvar = wrf_DPA.createVariable('Gust10', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
        outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
        outvar.standard_name = 'wind_speed_of_gust'
        outvar.long_name = "10-m wind gust"
        outvar.units = "m s-1"
        outvar.coordinates = "XLAT XLONG"
        outvar.grid_mapping = "Lambert_Conformal"

except:
        #### Para la PRM YSU que calcula TKE pero si la velocidad de friccion (Born et al. 2012)
        U, V = getvar(wrfout, 'uvmet') #Viento 3D
        Z = getvar(wrfout, 'height_agl') #Altura sobre el nivel del suelo
        u_30m = interplevel(U, Z, 30)
        v_30m = interplevel(V, Z, 30)
        V_30m = np.sqrt(u_30m**2 + v_30m**2)
        u_fric = getvar(wrfout, 'UST')#"U* IN SIMILARITY THEORY"
        data = np.abs(V_30m) + 3.0 * 2.4 * u_fric
        outvar = wrf_DPA.createVariable('Gust10', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
        outvar[0, :, :] = add_fill_value(data[:, :], fill_value)
        outvar.standard_name = 'wind_speed_of_gust'
        outvar.long_name = "10-m wind gust"
        outvar.units = "m s-1"
        outvar.coordinates = "XLAT XLONG"
        outvar.grid_mapping = "Lambert_Conformal"


# W en puntos de masa [m s-1]
W = getvar(wrfout, "wa")
P = getvar(wrfout, "pressure")
W_int = np.full((2, nlat, nlon), np.nan)
for i,lev in enumerate(LEV):
    W_int[i, :, :] = interplevel(W, P, lev)

outvar = wrf_DPA.createVariable('W', np.float32, ('time', 'lev', 'y', 'x'), zlib = True, fill_value = fill_value)
outvar[0,:,:,:] = add_fill_value(W_int[:,:,:], fill_value)
outvar.standard_name = "upward_air_velocity"
outvar.long_name = "z-wind component"
outvar.units = "m s-1"
outvar.coordinates = "XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"


wrf_DPA.close()
wrfout.close()


end = datetime.now()
dif = end - inicio

print('Para generar ' + post_name + ' tardamos ' + str(dif.seconds) + ' segundos')


