# -*- coding: utf-8 -*-
import numpy as np
from netCDF4 import Dataset
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from cpt_convert import loadCPT
from matplotlib.colors import LinearSegmentedColormap
import datetime
import sys
import locale
from PIL import Image

locale.setlocale(locale.LC_ALL, 'es_ES.utf8')

#######################
# Variables a definir #
#######################

data_path = '/data/fcutraro/CRTM/salida_HM/'
out_path = '/data/fcutraro/CRTM/graficado/pruebas/'
aux_path = '/data/fcutraro/CRTM/graficado/aux/'

#Lista con canales que se procesaron con el CRTM
canales = ['07', '08', '09', '10', '11', '12', '13', '14', '15', '16']

#Cargo las barras de colores
cpt_topes = loadCPT(aux_path + 'smn_topes.cpt')
cmap_topes = LinearSegmentedColormap('cpt', cpt_topes)
    
cpt_vapor = loadCPT(aux_path + 'SVGAWVX_TEMP.cpt')
cmap_vapor = LinearSegmentedColormap('cpt', cpt_vapor)

#WRFOUT a procesar
filename = sys.argv[1] 

data = filename.split('.')

dia_ini = datetime.datetime.strptime(data[-3], '%Y%m%d_%H%M%S')
fcst = data[-2]


#Paso el dia de UTC a HOA
dia_ini_HOA = dia_ini + datetime.timedelta(hours = -3)
dia_prono = dia_ini + datetime.timedelta(hours = int(fcst)) + datetime.timedelta(hours = -3)

print(dia_prono)

#Leo el dato del CRTM
try:
    CRTM = Dataset(filename)
except:
    sys.exit()

lats = CRTM.variables['XLAT'][:]
lons = CRTM.variables['XLONG'][:]

################
## Projeccion ##
################

#m = Basemap(projection='lcc', width=4000*999, height=4000*1249, resolution='i',
#             lat_1=-35,lat_2=-35,lat_0=-35,lon_0=-65)
m = Basemap(projection='mill', lat_ts = 10, llcrnrlon = -80,  urcrnrlon = -50, 
            llcrnrlat = lats.min(), urcrnrlat = -19, resolution = 'i')

#Paso la latitud y longitud a coordenadas de la proyeccion
x, y = m(lons, lats)

y_size = 24
x_size = 8

#Itero sobre los canales a procesar    
for canal in canales:

  #Elijo la barra de colores en funcion del canal
  if canal in ['08', '09', '10']:
    colorbar = cmap_vapor
    vmin = -112.15
    vmax = 56.85
  else:
    colorbar = cmap_topes
    vmin = -90
    vmax = 50

  fig = plt.figure(figsize=(y_size,x_size))

  m.drawstates(color = 'white', linewidth = 0.5) #Grafica provincias
  m.drawcountries(color = 'white', linewidth = 0.5) #Grafica paises
  m.drawcoastlines(color = 'white', linewidth = 0.5) #Grafica costas
  m.drawparallels(np.arange(-60., -5., 5.), labels = [1, 0, 1, 1], linewidth = 0.5, color = 'white', fontsize = 7) #Grafica paralelos
  m.drawmeridians(np.arange(-80., -40., 5.), labels = [1, 1, 0, 1], linewidth = 0.5, color = 'white', fontsize = 7) #Grafica meridianos

  #Grafico la imagen satelital sintetica
  m.pcolormesh(x, y, CRTM.variables['CH'+str(canal)][:]-273.15, cmap = colorbar, vmin = vmin, vmax = vmax)

  #Barra de colores 
  cbar = m.colorbar(ticks = np.arange(-90,51,10), pad = "5%")
  cbar.ax.tick_params(labelsize = 6)

  #Titulo
  plt.title(dia_prono.strftime('Imagen satelital sintética GOES16 ABI Canal ' + canal + '\nVálido para el %d de %B de %Y a las %H HOA'))

  plt.figtext(0.4, 0.05, (dia_ini_HOA.strftime("Inicializado el %-d/%-m/%Y %H HOA")))

  #Archivo de salida
  fileout = out_path + 'CRTM_C' + str(canal) + dia_ini.strftime('.%Y%m%d_%H%M%S.') + fcst

  #Guardo
  plt.savefig(fileout + '.png', bbox_inches='tight')

  plt.close()

  ''' LOGO '''
  mimage = Image.open(fileout + '.png')
  limage = Image.open(aux_path + 'LOGO_SMN.png')
 
  # resize logo
  wsize = int(min(mimage.size[0], mimage.size[1]) * 0.10)
  wpercent = (wsize / float(limage.size[0]))
  hsize = int((float(limage.size[1]) * float(wpercent)))
 
  simage = limage.resize((wsize, hsize))
  mbox = mimage.getbbox()
  sbox = simage.getbbox()
 
  # right bottom corner
  box = (sbox[2]+7, sbox[3]-8)
  mimage.paste(simage, box, simage)
  mimage.save(fileout + '_SMN.png')

                        
