#!/usr/bin/env python

# Version 1
# Para graficar las observaciones asimiladas separadas por la fuente, a partir de txt generados con

# En 10.10.23.20:
# source activate npp
# export dir=/ms270/reldata/rra_cheyenne/OBS/
# ln -sf ${dir}/20181110_12/obs_20181110_12_asimiladas.dat ./obs.dat
# ./obs_dat2txt.exe obs.dat > obs.txt
# ./graficar_obs_desdetxt.py YYYY MM DD HH Min obs.txt


import numpy as np
import os
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from PIL import Image

# no hace falta setear el path, eso quedo en el .sh, asique lo dejo en blanco
path_in = ''
path_out = ''

# Ingreso la fecha como argumentos:
parser = argparse.ArgumentParser(description='Year Month Day Hour Minute OBSASSIM')
parser.add_argument('Year',type=int)
parser.add_argument('Month',type=int)
parser.add_argument('Day',type=int)
parser.add_argument('Hour',type=int)
parser.add_argument('Minute',type=int)
parser.add_argument('OBSASSIM',type=str)
args = parser.parse_args()

Y = args.Year
M = args.Month
D = args.Day
H = args.Hour
Mi = args.Minute

YY = str(Y)
MM = str(M).zfill(2)
DD = str(D).zfill(2)
HH = str(H).zfill(2)
MiMi = str(Mi).zfill(2)

FECHA_save = YY + MM + DD + '_' +  HH + MiMi + '00' 
FECHA_title = YY + '-' + MM + '-' + DD + ' ' + HH + ':' + MiMi + 'Z'

# Abro los archivos
NAME = args.OBSASSIM
f2 = np.loadtxt(path_in + NAME)

# Defino limites para el grafico de todo el pais
lat_min = -60
lat_max = -10
lon_min = 265
lon_max = 325

# Hago la seleccion para cada tipo de dato y voy graficando
# para hacer mas eficiente el uso de memoria


fig = plt.figure(figsize=(15,10))
m = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max,projection='cyl',resolution = 'l')

################## RADAR REFLECTIVIDAD (4001 primera columna)
mascara = np.where( f2[:,0] == 4001 )[0]
obs = f2[mascara,:]
lat_obs = obs[:,2]
lon_obs = obs[:,1] + 360
num_radar = str(np.size(lon_obs))

lons, lats = m(lon_obs, lat_obs)
l1 = m.scatter(lons, lats, marker = 'o', color='lightsalmon', edgecolors='lightsalmon', s=5, label='RADAR')

################# ASCATW (20) -> no estamos asimilando, pero lo dejo preparado:cambiar colores
# mascara = np.where( f[:,6] == 20 )[0]
# obs = f2[mascara,:]
# lat_obs = obs[:,2]
# lon_obs = obs[:,1]
# num_ascat = str(np.size(lon_obs))

#lons, lats = m(lon_obs, lat_obs)
#l2 = m.scatter(lons, lats, marker = 'o', color='lightsalmon', edgecolors='salmon', s=20, label='ASCAT')

################## SATWND (4)
mascara = np.where( f2[:,6] == 4 )[0]
obs = f2[mascara,:]
lat_obs = obs[:,2]
lon_obs = obs[:,1]
num_satwnd = str(np.size(lon_obs))

lons, lats = m(lon_obs, lat_obs)
l3 = m.scatter(lons, lats, marker = 'x', color='cornflowerblue', s=20, label='GOES')

################## AIRSRT (21)
mascara = np.where( f2[:,6] == 21 )[0]
obs = f2[mascara,:]
lat_obs = obs[:,2]
lon_obs = obs[:,1]
num_airs = str(np.size(lon_obs))

lons, lats = m(lon_obs, lat_obs)
l4 = m.scatter(lons, lats, marker = 'x', color='magenta', s=20, label='SAT POLAR')

################## ADPAUT (22)
mascara = np.where( f2[:,6] == 22 )[0]
obs = f2[mascara,:]
lat_obs = obs[:,2]
lon_obs = obs[:,1]
num_adpaut = str(np.size(lon_obs))

lons, lats = m(lon_obs, lat_obs)
l5 = m.scatter(lons, lats, marker = '*', color='limegreen', s=15, label='AUTOMATICA')

################## AIRCFT (3)
mascara = np.where( f2[:,6] == 3 )[0]
obs = f2[mascara,:]
lat_obs = obs[:,2]
lon_obs = obs[:,1]
num_aircft = str(np.size(lon_obs))

lons, lats = m(lon_obs, lat_obs)
l6 = m.scatter(lons, lats, marker = '<', color='r', s=20, label='AVION')

################## ADPSFC (8)
mascara = np.where( f2[:,6] == 8 )[0]
obs = f2[mascara,:]
lat_obs = obs[:,2]
lon_obs = obs[:,1]
num_adpsfc = str(np.size(lon_obs))

lons, lats = m(lon_obs, lat_obs)
l7 = m.scatter(lons, lats, marker = 'D', color='blueviolet', edgecolors='gray', s=20, label='SUPERFICIE')

################## ADPUPA (1)
mascara = np.where( f2[:,6] == 1 )[0]
obs = f2[mascara,:]
lat_obs = obs[:,2]
lon_obs = obs[:,1]
num_adpupa = str(np.size(lon_obs))

lons, lats = m(lon_obs, lat_obs)
l8 = m.scatter(lons, lats, marker = 'o', color='gold', edgecolors='gray', s=20, label='SONDEOS')

################## SFCSHP (9)
mascara = np.where( f2[:,6] == 9 )[0]
obs = f2[mascara,:]
lat_obs = obs[:,2]
lon_obs = obs[:,1]
num_sfcshp = str(np.size(lon_obs))

lons, lats = m(lon_obs, lat_obs)
l9 = m.scatter(lons, lats, marker = '>', color='deepskyblue', edgecolors='slateblue', s=20, label='BARCOS')

####### Caracteristicas del grafico:
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.25)
m.drawparallels(np.arange(lat_min,lat_max,5.),labels=[1,0,0,0],fontsize=18)
m.drawmeridians(np.arange(lon_min,lon_max+1,5.),labels=[0,0,0,1],fontsize=18)

# incluir ascat cuando esten:
fig.legend([l1, l3, l4, l5, l6, l7, l8, l9], [l1.get_label(), l3.get_label(), l4.get_label(), l5.get_label(), l6.get_label(), l7.get_label(), l8.get_label(), l9.get_label()], fontsize=18, loc="center right")


fig.suptitle('Observaciones asimiladas el ' + FECHA_title + '\n RADAR=' + num_radar +', GOES=' + num_satwnd +', SAT POLAR=' + num_airs +', AUTOMATICA=' + num_adpaut  + '\n AVION=' + num_aircft + ', SUPERFICIE=' + num_adpsfc + ', SONDEOS=' + num_adpupa + ', BARCOS=' + num_sfcshp, fontsize=20)

fig_obs = 'model.WRF_ENS_4km.' + FECHA_save + '.000_Assimilated_Observations.jpg'
fig.savefig( path_out + fig_obs, quality = 80, bbox_inches = 'tight' )
plt.close('all')

#''' LOGO '''
mimage = Image.open(path_out + fig_obs)
limage = Image.open('/data/oper/wrfutilV3.0/bin/python/LOGO_SMN.png')
# resize logo
wsize = int(min(mimage.size[0], mimage.size[1]) * 0.08)
wpercent = (wsize / float(limage.size[0]))
hsize = int((float(limage.size[1]) * float(wpercent)))

simage = limage.resize((wsize, hsize))
mbox = mimage.getbbox()
sbox = simage.getbbox()

# pongo uno en cada mapita
box = (sbox[2]+20, sbox[3]+50)  # para la derecha, para abajo (desde arriba a la izq
mimage.paste(simage, box, simage)

mimage.save(path_out + fig_obs)







