# -*- coding: utf-8 -*-

# ------------------------------------------
# SMN - M.E. Dillon, M. Sacco 
#
# Genera un netcdf con la estructura de un wrfout
# con el perfil vertical del RMSD update del letkf de u,v,T,q
# y ademas 2 graficos operativos:
#	Campos horizontales de RMSD update
#	Perfil vertical de RMSD update
# En todos los casos se sacan los puntos de los bordes laterales
# pasando de la dimension (1249,999) a la dimension (1230,980)
#
# Activar el entorno:
# source activate wrfutil
#
# Tambien se necesita tener en el entorno cargadas las siguientes variables:
# PATHOUT_PLOT - Directorio donde se escribiran los plots
# PATHOUT      - Directorio donde se escribiran las salidas
# FILEIN_ANA   - Path absoluto del analisis medio
# FILEIN_GUES  - Path absoluto del gues medio
#
# Ejecutar con:
# python write_and_plot_smn_asim_rmsd_update.py
# ------------------------------------------

import numpy as np
import os
from wrf import getvar, interplevel
from netCDF4 import Dataset
from datetime import datetime
from datetime import timedelta

inicio = datetime.now()
print(str(datetime.now()) + " Empezamos")
path_out = os.environ["PATHOUT"]
path_out_plot = os.environ["PATHOUT_PLOT"]
aname = os.environ["FILEIN_ANA"]
guesme = os.environ["FILEIN_GUES"]

###################
# Defino la descripcion del archivo para escribir las salidas
###################
description_short = 'UPD'
description_long = 'RMSD Update (Ana-Gues) '

####################
#--- Abro los wrfout:
####################

wrf_ana = Dataset(aname, 'r')
fecha_ini = datetime.strptime(wrf_ana.START_DATE, '%Y-%m-%d_%H:%M:%S')
fecha_fc = fecha_ini
Yi = str(fecha_ini.year)
Mi = str(fecha_ini.month).zfill(2)
Di = str(fecha_ini.day).zfill(2)
Hi = str(fecha_ini.hour).zfill(2)
hf_out = str(int((fecha_fc - fecha_ini).total_seconds()/3600)).zfill(3)

wrf_gues = Dataset(guesme, 'r')

####################
#--- Creo el netcdf de salida:
####################

print(str(datetime.now()) + " Creamos el archivo de salida")

tipo = 'LEV'
post_dir = path_out
post_name = 'model.WRF_' + description_short + '_4km.' + Yi + Mi + Di + '_' + Hi + '0000.' + hf_out + '.OPER' + tipo +'.nc'
wrf_DPA = Dataset(post_dir + post_name, 'w')

#--- Dimensiones:

XLAT = getvar(wrf_ana, "XLAT")   # Latitud
XLON = getvar(wrf_ana, "XLONG")  # Longitud
ZNU = getvar(wrf_ana, "ZNU")	# Niveles verticales "eta values on half (mass) levels"
				# bottom_top(=44) es la dimension de u,v,t,q
XLAT = XLAT[10:1240,10:990]	# paso de dimension (1249,999) a dimension (1230,980) para sacar los bordes
XLON = XLON[10:1240,10:990]

XTIME = wrf_ana.variables["XTIME"]
nlev = len(ZNU)
nlat, nlon = XLAT.shape

wrf_DPA.createDimension("time", None)
wrf_DPA.createDimension("lev", nlev)

#--- Atributos globales:

atributos_globales = {}
atributos_globales['title'] = 'Python PostProcessing for SMN WRF-ARW ' + description_long + tipo
atributos_globales["institution"] = "Servicio Meteorologico Nacional"
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

#--- Atributos y Escritura de cada variable:

atributo_variable = {}

outvar = wrf_DPA.createVariable('ZNU', np.float32, ('lev'), zlib = True)
outvar[:] = ZNU[:]
outvar.long_name = "eta values on half (mass) levels"
outvar.axis = "Z"

outvar = wrf_DPA.createVariable("XTIME", np.float32, ("time"), zlib = True)
outvar[0] = XTIME[:]/60
outvar.long_name = fecha_ini.strftime("hours since %Y-%m-%d %H:%M:%S")
outvar.units = fecha_ini.strftime("hours since %Y-%m-%d %H:%M:%S")
outvar.standard_name = "time"
outvar.axis = "T"
outvar.calendar = "gregorian"


#### 
# Leemos las variables de entrada
# Calculamos el RSMD
####
print (str(datetime.now())+" Calculando el RMSD update")

# Specific Humidity:
ana = getvar(wrf_ana,"QVAPOR")*1e3
gue = getvar(wrf_gues,"QVAPOR")*1e3
dif = ana-gue	# dim lev,lat,lon
dif = dif[:,10:1240,10:990]	# me quedo con el dominio sin los bordes
ana = 0
gue = 0

rmsd_q = np.ma.sqrt( np.ma.mean(dif*dif,axis=0) )
rmsd_q_vert = np.ma.sqrt( np.ma.mean( np.ma.mean(dif*dif,axis=1), axis=1) )

# Temperature:
ana = getvar(wrf_ana,"tk")
gue = getvar(wrf_gues,"tk")
dif = ana-gue   # dim lev,lat,lon
dif = dif[:,10:1240,10:990]     # me quedo con el dominio sin los bordes
ana = 0
gue = 0

rmsd_t = np.ma.sqrt( np.ma.mean(dif*dif,axis=0) )
rmsd_t_vert = np.ma.sqrt( np.ma.mean( np.ma.mean(dif*dif,axis=1), axis=1) )

# Wind components:
viento_ana = getvar(wrf_ana, "uvmet")
Umet_ana = viento_ana[0]        # U meteorologico [m s-1]
Vmet_ana = viento_ana[1]        # V meteorologico [m s-1]
viento_gues = getvar(wrf_gues, "uvmet")
Umet_gues = viento_gues[0]        
Vmet_gues = viento_gues[1]        

udif = Umet_ana-Umet_gues
vdif = Vmet_ana-Vmet_gues
udif = udif[:,10:1240,10:990]     # me quedo con el dominio sin los bordes
vdif = vdif[:,10:1240,10:990]     # me quedo con el dominio sin los bordes
Umet_ana = 0
Vmet_ana = 0
Umet_gues = 0
Vmet_gues = 0

rmsd_u = np.ma.sqrt( np.ma.mean(udif*udif,axis=0) )
rmsd_u_vert = np.ma.sqrt( np.ma.mean( np.ma.mean(udif*udif,axis=1), axis=1) )

rmsd_v = np.ma.sqrt( np.ma.mean(vdif*vdif,axis=0) )
rmsd_v_vert = np.ma.sqrt( np.ma.mean( np.ma.mean(vdif*vdif,axis=1), axis=1) )

# Reflectivity:
DBZA = getvar(wrf_ana, "dbz")
DBZG = getvar(wrf_gues, "dbz")
dif = DBZA-DBZG
dif = dif[:,10:1240,10:990]     # me quedo con el dominio sin los bordes
DBZA = 0
DBZG = 0

rmsd_Z = np.ma.sqrt( np.ma.mean(dif*dif,axis=0) )
rmsd_Z_vert = np.ma.sqrt( np.ma.mean( np.ma.mean(dif*dif,axis=1), axis=1) )

# Geopotential:
ana = getvar(wrf_ana, "geopt")
gue = getvar(wrf_gues, "geopt")
dif = ana - gue
dif = dif[:,10:1240,10:990]     # me quedo con el dominio sin los bordes
ana = 0
gue = 0

rmsd_geopt = np.ma.sqrt( np.ma.mean(dif*dif,axis=0) )
rmsd_geopt_vert = np.ma.sqrt( np.ma.mean( np.ma.mean(dif*dif,axis=1), axis=1) )

####
# Escribimos los rmsd verticales
####

varList=[(rmsd_q_vert, "Q_update", "Specific Humidity RMSD update (Ana - Guess mean)",  "specific_humidity",    "g Kg-1") ,
         (rmsd_t_vert, "T_update", "Temperature RMSD update (Ana - Guess mean)",    "air_temperature",  "K"),
         (rmsd_u_vert, "U_update", "Zonal Wind Component RMSD update (Ana - Guess mean)",   "eastward_wind",    "m s-1"),
         (rmsd_v_vert, "V_update", "Meridional Wind Component RMSD update (Ana - Guess mean)",  "northward_wind",   "m s-1"),
	     (rmsd_Z_vert, "Z_update", "Z Reflectivity RMSD update (Ana - Guess mean)", 'equivalent_reflectivity_factor',   "dbZ"),
	     (rmsd_geopt_vert, "GEOPT_update", "Geopotential RMSD update (Ana - Guess mean)",   "geopotential_height",  "m2 s-2")]

for (var, alias, lname, stdname, units) in varList:
    print(str(datetime.now()) + " Guardando -> " + str(alias))
    outvar = wrf_DPA.createVariable(alias, np.float32, ('time', 'lev'), zlib = True)
    outvar[0, :] = var
    outvar.standard_name = stdname
    outvar.long_name = lname
    outvar.units = units

print (str(datetime.now())+" Cerrando archivo de salida")
wrf_DPA.close()
#wrf_ana.close()  # no lo cierro todavia porque necesito otra forma de leer lat lon para grafico
wrf_gues.close()

escritura = datetime.now()
dif = escritura - inicio
print('Para generar ' + post_name + ' tardamos ' + str(dif.seconds) + ' segundos')

####
# Graficamos los rmsd verticales y el campo horizontal
####
print (str(datetime.now())+" Graficando los RMSD updates")

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from PIL import Image	# para incluir el logo
import locale

lats = wrf_ana.variables['XLAT'][:] 
lats = lats[0,10:1240,10:990]
lons = wrf_ana.variables['XLONG'][:] 
lons = lons[0,10:1240,10:990]

fig_campos = 'model.WRF_' + description_short + '_4km.' + Yi + Mi + Di + '_' + Hi + '0000.' + hf_out + '_camposRMSD.png'
fig_perfiles = 'model.WRF_' + description_short + '_4km.' + Yi + Mi + Di + '_' + Hi + '0000.' + hf_out + '_perfilesRMSD.png'

# -------------------- Campo horizontal:
plt.figure(figsize = (15, 10))

locale.setlocale(locale.LC_ALL,'es_ES.UTF-8')
levels_uv = [0.5, 1, 2, 3, 5, 8, 10]
levels_t = [0.5, 1, 2, 3, 4, 5]
levels_q = [0.2, 0.4, 0.6, 0.8, 1, 1.5, 2]
levels_z = [1, 3, 5, 8, 10, 15, 20]
levels_geo = [10, 40, 80, 120, 160, 200 ]

plt.subplot(2,3,1)
m = Basemap(projection = 'mill', lat_ts = 10, llcrnrlon = -80,  \
    urcrnrlon = -50, llcrnrlat = lats.min(), urcrnrlat = -19, \
    resolution = 'l')	# con resolution='h' me da libGL error
m.shadedrelief()
m.drawstates(color = 'grey', linewidth = 0.7)
m.drawcountries(color = 'grey', linewidth = 0.7)
m.drawparallels(np.arange(-60., -5., 5.), labels = [1, 0, 1, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
m.drawmeridians(np.arange(-80., -40., 5.), labels = [1, 1, 0, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
x, y = m(lons, lats)
m.contourf(x, y, rmsd_u, levels_uv, cmap = plt.cm.nipy_spectral, extend = 'max')
cbar = m.colorbar(ticks = levels_uv, pad = "5%")
cbar.ax.tick_params(labelsize = 8)
plt.title('Comp zonal del viento U [m s-1]')

plt.subplot(2,3,2)
m = Basemap(projection = 'mill', lat_ts = 10, llcrnrlon = -80,  \
    urcrnrlon = -50, llcrnrlat = lats.min(), urcrnrlat = -19, \
    resolution = 'l')   # con resolution='h' me da libGL error
m.shadedrelief()
m.drawstates(color = 'grey', linewidth = 0.7)
m.drawcountries(color = 'grey', linewidth = 0.7)
m.drawparallels(np.arange(-60., -5., 5.), labels = [1, 0, 1, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
m.drawmeridians(np.arange(-80., -40., 5.), labels = [1, 1, 0, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
x, y = m(lons, lats)
m.contourf(x, y, rmsd_v, levels_uv, cmap = plt.cm.nipy_spectral, extend = 'max')
cbar = m.colorbar(ticks = levels_uv, pad = "5%")
cbar.ax.tick_params(labelsize = 8)
plt.title('Comp meridional del viento V [m s-1]')

plt.subplot(2,3,4)
m = Basemap(projection = 'mill', lat_ts = 10, llcrnrlon = -80,  \
    urcrnrlon = -50, llcrnrlat = lats.min(), urcrnrlat = -19, \
    resolution = 'l')   # con resolution='h' me da libGL error
m.shadedrelief()
m.drawstates(color = 'grey', linewidth = 0.7)
m.drawcountries(color = 'grey', linewidth = 0.7)
m.drawparallels(np.arange(-60., -5., 5.), labels = [1, 0, 1, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
m.drawmeridians(np.arange(-80., -40., 5.), labels = [1, 1, 0, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
x, y = m(lons, lats)
m.contourf(x, y, rmsd_t, levels_t, cmap = plt.cm.YlOrRd, extend = 'max')
cbar = m.colorbar(ticks = levels_t, pad = "5%")
cbar.ax.tick_params(labelsize = 8)
plt.title('Temperatura [K]')

plt.subplot(2,3,5)
m = Basemap(projection = 'mill', lat_ts = 10, llcrnrlon = -80,  \
    urcrnrlon = -50, llcrnrlat = lats.min(), urcrnrlat = -19, \
    resolution = 'l')   # con resolution='h' me da libGL error
m.shadedrelief()
m.drawstates(color = 'grey', linewidth = 0.7)
m.drawcountries(color = 'grey', linewidth = 0.7)
m.drawparallels(np.arange(-60., -5., 5.), labels = [1, 0, 1, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
m.drawmeridians(np.arange(-80., -40., 5.), labels = [1, 1, 0, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
x, y = m(lons, lats)
m.contourf(x, y, rmsd_q, levels_q, cmap = plt.cm.PuRd, extend = 'max')
cbar = m.colorbar(ticks = levels_q, pad = "5%")
cbar.ax.tick_params(labelsize = 8)
plt.title('Humedad especifica [g Kg-1]')

plt.subplot(2,3,3)
m = Basemap(projection = 'mill', lat_ts = 10, llcrnrlon = -80,  \
    urcrnrlon = -50, llcrnrlat = lats.min(), urcrnrlat = -19, \
    resolution = 'l')   # con resolution='h' me da libGL error
m.shadedrelief()
m.drawstates(color = 'grey', linewidth = 0.7)
m.drawcountries(color = 'grey', linewidth = 0.7)
m.drawparallels(np.arange(-60., -5., 5.), labels = [1, 0, 1, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
m.drawmeridians(np.arange(-80., -40., 5.), labels = [1, 1, 0, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
x, y = m(lons, lats)
m.contourf(x, y, rmsd_Z, levels_z, cmap = plt.cm.magma_r, extend = 'max')
cbar = m.colorbar(ticks = levels_z, pad = "5%")
cbar.ax.tick_params(labelsize = 8)
plt.title('Reflectividad [dBZ]')

plt.subplot(2,3,6)
m = Basemap(projection = 'mill', lat_ts = 10, llcrnrlon = -80,  \
    urcrnrlon = -50, llcrnrlat = lats.min(), urcrnrlat = -19, \
    resolution = 'l')   # con resolution='h' me da libGL error
m.shadedrelief()
m.drawstates(color = 'grey', linewidth = 0.7)
m.drawcountries(color = 'grey', linewidth = 0.7)
m.drawparallels(np.arange(-60., -5., 5.), labels = [1, 0, 1, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
m.drawmeridians(np.arange(-80., -40., 5.), labels = [1, 1, 0, 1], linewidth = 0.5, color = 'grey', fontsize = 7)
x, y = m(lons, lats)
m.contourf(x, y, rmsd_geopt, levels_geo, cmap = plt.cm.rainbow, extend = 'max')
cbar = m.colorbar(ticks = levels_geo, pad = "5%")
cbar.ax.tick_params(labelsize = 8)
plt.title('Altura Geopotencial [m2 s-2]')



plt.figtext(0.2, 0.95, fecha_ini.strftime('LETKF WRF - RMSD update de la media del ensamble (ANA-GUES) \n\
Válido para el %d de %B de %Y a las %H UTC'), size = 'xx-large')

os.makedirs(path_out_plot, exist_ok = True)
plt.savefig(path_out_plot + fig_campos, quality = 95, bbox_inches = 'tight')
plt.close('all')

#''' LOGO '''
mimage = Image.open(path_out_plot + fig_campos)
limage = Image.open('/data/oper/wrfutilV3.0/bin/python/LOGO_SMN.png')
# resize logo
wsize = int(min(mimage.size[0], mimage.size[1]) * 0.05)
wpercent = (wsize / float(limage.size[0]))
hsize = int((float(limage.size[1]) * float(wpercent)))

simage = limage.resize((wsize, hsize))
mbox = mimage.getbbox()
sbox = simage.getbbox()

# pongo uno en cada mapita
box = (sbox[2]+10, sbox[3]+80)	# para la derecha, para abajo (desde arriba a la izq
mimage.paste(simage, box, simage)
box = (sbox[2]+415, sbox[3]+80)
mimage.paste(simage, box, simage)
box = (sbox[2]+825, sbox[3]+80)
mimage.paste(simage, box, simage)
box = (sbox[2]+10, sbox[3]+500)  
mimage.paste(simage, box, simage)
box = (sbox[2]+415, sbox[3]+500)
mimage.paste(simage, box, simage)
box = (sbox[2]+825, sbox[3]+500)
mimage.paste(simage, box, simage)

mimage.save(path_out_plot + fig_campos)
wrf_ana.close()

plot1 = datetime.now()
dif = plot1 - escritura 
print('Para generar ' + fig_campos + ' tardamos ' + str(dif.seconds) + ' segundos')


# -------------------- Campo horizontal:
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.figure(figsize = (15, 10))

locale.setlocale(locale.LC_ALL,'es_ES.UTF-8')
vertical_levels = np.arange(1,45)	# en principio no lo paso a hPa

plt.subplot(2,3,1)
plt.plot(rmsd_u_vert,vertical_levels,'seagreen',linewidth=2)
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.ylabel('Nivel vertical del modelo', fontsize=10)
plt.xlim((0.,levels_uv[-1]))
plt.tick_params(axis='both', which='major', labelsize=9)
plt.tick_params(axis='both', which='minor', labelsize=9)
plt.title('Comp zonal del viento U [m s-1]')

plt.subplot(2,3,2)
plt.plot(rmsd_v_vert,vertical_levels,'darkcyan',linewidth=2)
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.xlim((0.,levels_uv[-1]))
plt.tick_params(axis='both', which='major', labelsize=9)
plt.tick_params(axis='both', which='minor', labelsize=9)
plt.title('Comp meridional del viento V [m s-1]')

plt.subplot(2,3,4)
plt.plot(rmsd_t_vert,vertical_levels,'red',linewidth=2)
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.ylabel('Nivel vertical del modelo', fontsize=10)
plt.xlim((0.,levels_t[-1]))
plt.tick_params(axis='both', which='major', labelsize=9)
plt.tick_params(axis='both', which='minor', labelsize=9)
plt.title('Temperatura [K]')

plt.subplot(2,3,5)
plt.plot(rmsd_q_vert,vertical_levels,'darkorchid',linewidth=2)
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.xlim((0.,levels_q[-1]))
plt.tick_params(axis='both', which='major', labelsize=9)
plt.tick_params(axis='both', which='minor', labelsize=9)
plt.title('Humedad especifica [g Kg-1]')

plt.subplot(2,3,3)
plt.plot(rmsd_Z_vert,vertical_levels,'dodgerblue',linewidth=2)
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.xlim((0.,levels_z[-1]))
plt.tick_params(axis='both', which='major', labelsize=9)
plt.tick_params(axis='both', which='minor', labelsize=9)
plt.title('Reflectividad [dBZ]')

plt.subplot(2,3,6)
plt.plot(rmsd_geopt_vert,vertical_levels,'maroon',linewidth=2)
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.xlim((0.,levels_geo[-1]))
plt.tick_params(axis='both', which='major', labelsize=9)
plt.tick_params(axis='both', which='minor', labelsize=9)
plt.title('Altura Geopotencial [m2 s-2]')


plt.figtext(0.2, 0.95, fecha_ini.strftime('LETKF WRF - RMSD update de la media del ensamble (ANA-GUES) \n\
Válido para el %d de %B de %Y a las %H UTC'), size = 'xx-large')

os.makedirs(path_out_plot, exist_ok = True)
plt.savefig(path_out_plot + fig_perfiles, quality = 95, bbox_inches = 'tight')
plt.close('all')

#''' LOGO '''
mimage = Image.open(path_out_plot + fig_perfiles)
limage = Image.open('/data/oper/wrfutilV3.0/bin/python/LOGO_SMN.png')
# resize logo
wsize = int(min(mimage.size[0], mimage.size[1]) * 0.05)
wpercent = (wsize / float(limage.size[0]))
hsize = int((float(limage.size[1]) * float(wpercent)))

simage = limage.resize((wsize, hsize))
mbox = mimage.getbbox()
sbox = simage.getbbox()

# pongo uno en cada mapita
box = (sbox[2]+300, sbox[3]+80)  # para la derecha, para abajo (desde arriba a la izq
mimage.paste(simage, box, simage)
box = (sbox[2]+710, sbox[3]+80)
mimage.paste(simage, box, simage)
box = (sbox[2]+1115, sbox[3]+80)
mimage.paste(simage, box, simage)
box = (sbox[2]+300, sbox[3]+500)
mimage.paste(simage, box, simage)
box = (sbox[2]+710, sbox[3]+500)
mimage.paste(simage, box, simage)
box = (sbox[2]+1115, sbox[3]+500)
mimage.paste(simage, box, simage)


mimage.save(path_out_plot + fig_perfiles)


plot2 = datetime.now()
dif = plot2 - plot1
print('Para generar ' + fig_perfiles + ' tardamos ' + str(dif.seconds) + ' segundos')


