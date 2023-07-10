## Script que calcula UH 
import glob
import sys
#import argparse
import os
import numpy as np
from netCDF4 import Dataset
from wrf import getvar, udhel, ALL_TIMES

import argparse
from datetime import datetime
from datetime import timedelta
import pyproj
from funciones_post import add_fill_value

#########
#caso='20181110_36h'
#minutos='10'
#fill_value = 1e20
#filename = sorted(glob.glob('/data/mili/workspace3/wrfutilV3.1/RUNs/deterministico_UH/WRF/00_'+minutos+'min_36h/UH/UH_2018-11-11_12_00_00'))
#path_out = '/data/mili/UH/'+caso+'/post_'+minutos+'min/' 

######
path_out = os.environ["PATHOUT_UH"]
filename = os.environ["FILEIN_UH"]
fill_value = 1e20

wrfout = Dataset(filename, "r")
#print(ifile)
fecha_ini = datetime.strptime(wrfout.START_DATE, "%Y-%m-%d_%H:%M:%S")
fecha_fc = datetime.strptime(filename[-19:], "%Y-%m-%d_%H_%M_%S")
Yi = str(fecha_ini.year)
Mi = str(fecha_ini.month).zfill(2)
Di = str(fecha_ini.day).zfill(2)
Hi = str(fecha_ini.hour).zfill(2)
hf_out = str(int((fecha_fc - fecha_ini).total_seconds()/3600)).zfill(3)

####################
#--- Obtengo variables directas del wrfout:
####################
XLAT = wrfout.variables["XLAT"]
XLONG = wrfout.variables["XLONG"]
XTIME = wrfout.variables["XTIME"]
ntime, nlat, nlon = XLAT.shape

post_dir = path_out
post_name = "model.WRF_DET_UH_4km." + Yi + Mi + Di + "_" + Hi + "0000." + hf_out + ".OPER.nc"
wrf_DPA = Dataset(post_dir + post_name, "w")

#--- Dimensiones:
wrf_DPA.createDimension("time", None) 
wrf_DPA.createDimension("y", nlat)
wrf_DPA.createDimension("x", nlon)

#--- Atributos globales:
atributos_globales = {}
atributos_globales["title"] = "Python PostProcessing for SMN WRF-ARW Deterministic "
atributos_globales["institution"] = "Servicio Meteorologico Nacional"
atributos_globales["source"] = wrfout.TITLE
atributos_globales["creation_date"] = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
atributos_globales["start_lat"] = XLAT[0,0,0]
atributos_globales["start_lon"] = XLONG[0,0,0]
atributos_globales["end_lat"] = XLAT[0,nlat-1,nlon-1]
atributos_globales["end_lon"] = XLONG[0,nlat-1,nlon-1]
atributos_globales["MAP_PROJ"] = wrfout.MAP_PROJ
atributos_globales["MAP_PROJ_CHAR"] = wrfout.MAP_PROJ_CHAR
atributos_globales["STAND_LON"] = wrfout.STAND_LON
atributos_globales["CEN_LON"] = wrfout.CEN_LON
atributos_globales["CEN_LAT"] = wrfout.CEN_LAT
atributos_globales["TRUELAT1"] = wrfout.TRUELAT1
atributos_globales["TRUELAT2"] = wrfout.TRUELAT2
atributos_globales["DX"] = wrfout.DX
atributos_globales["DY"] = wrfout.DY
atributos_globales["START_DATE"] = wrfout.START_DATE
atributos_globales["FCST_HOURS"] = hf_out
atributos_globales["VALID_DATE"] = fecha_fc.strftime("%Y-%m-%d_%H:%M:%S")
atributos_globales['Conventions'] = "CF-1.8"

wrf_DPA.setncatts(atributos_globales)
#--- Variable que define la proyeccion
lcc = pyproj.Proj(proj="lcc", lat_1 = wrfout.TRUELAT1, lat_2 = wrfout.TRUELAT2, lat_0 = wrfout.CEN_LAT, lon_0 = wrfout.CEN_LON, a=6370000, b=6370000)
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

outvar = wrf_DPA.createVariable("XTIME", np.float32, ("time"), zlib = True, fill_value = fill_value)
outvar[:] = XTIME[:]
outvar.long_name = fecha_ini.strftime("minutes since %Y-%m-%d %H:%M:%S")
outvar.units = fecha_ini.strftime("minutes since %Y-%m-%d %H:%M:%S")
outvar.standard_name = "time"
outvar.axis = "T"
outvar.calendar = "gregorian"

outvar = wrf_DPA.createVariable("x", np.float32, ("x"), zlib = True, fill_value = fill_value)
outvar[:] = x
outvar.axis = "X"
outvar.units = "m"
outvar.standard_name = "projection_x_coordinate"
outvar.long_name = "x-coordinate in projected coordinate system"

outvar = wrf_DPA.createVariable("y", np.float32, ("y"), zlib = True, fill_value = fill_value)
outvar[:] = y
outvar.axis = "Y"
outvar.units = "m"
outvar.standard_name = "projection_y_coordinate"
outvar.long_name = "y-coordinate in projected coordinate system"

### LAT, LON
#Latitud
outvar = wrf_DPA.createVariable("XLAT", np.float32, ("time", "y", "x"), zlib = True, fill_value = fill_value)
outvar[:, :, :] = XLAT[:, :, :]
outvar.standard_name = "latitude"
outvar.long_name = "Latitude"
outvar.units = "degrees_north"
#outvar.coordinates = "XLONG XLAT"

#Longitud
outvar = wrf_DPA.createVariable("XLONG", np.float32, ("time", "y", "x"), zlib = True, fill_value = fill_value)
outvar[:, :, :] = XLONG[:, :, :]
outvar.standard_name = "longitude"
outvar.long_name = "Longitude"
outvar.units = "degrees_east"
#outvar.coordinates = "XLONG XLAT"

### VARIABLE
data = getvar(wrfout,"updraft_helicity",bottom=2000.0, top=5000.0, timeidx=ALL_TIMES)
outvar = wrf_DPA.createVariable('UH', np.float32, ('time', 'y', 'x'), zlib = True, fill_value = fill_value)
if len(data.shape) == 3:
	outvar[:,:,:] = add_fill_value(data[:,:,:], fill_value) #data[:,:,:]
else:
	outvar[0,:,:] = add_fill_value(data[:,:], fill_value) #data[:,:]
outvar.standard_name = 'updraft_helicity'
outvar.long_name = "Updraft Helicity 2-5km"
outvar.units = "m2 s-2"
outvar.coordinates = "XTIME XLAT XLONG"
outvar.grid_mapping = "Lambert_Conformal"

wrf_DPA.close()
wrfout.close()


