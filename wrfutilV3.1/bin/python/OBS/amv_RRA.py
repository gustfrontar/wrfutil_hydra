# -*- coding: utf-8 -*-
import settings as st
import numpy as np
import pandas as pd
from netCDF4 import Dataset 
import glob
import os
import sys
#import matplotlib.pyplot as plt
from scipy.io import FortranFile
import datetime 

def superobbing(datos_asimilables):
    """
    Esta funcion realiza un superobbing a los datos de Atmospheric Motion Vectors
    que se obtienen del GOES-16.
    Se promedian los datos que se encuentran dentro de regiones definidas a partir 
    de una cantidad a definir de kilometros en la horizontal y de hPa en la vertical. 
    """

    #Parametros con que defino las cajas
    dist_horizontal = 30 #Km
    dist_vertical = 25 #hPa
    RT = 6371 #Radio de la Tierra

    #Array donde se guardan los datos que pasaron por el superobbing
    datos_final = np.full([datos_asimilables.shape[0]*2, 8], np.nan)

    dist_grad = dist_horizontal/111

    #Limites de las cajas donde se promedian las observaciones
    lat_so = np.arange(lat_min, lat_max + dist_grad - 1e-6, dist_grad)
    pres_so = np.arange(10, 1050+1e-6, dist_vertical)

    #print(lat_min, lat_max, lon_min, lon_max)

    l = 0
    for i in range(lat_so.size-1):
        dist_grad = 360/(RT*np.cos(np.deg2rad(lat_so[i]))*2*np.pi/dist_horizontal)
        lon_so = np.arange(lon_min, lon_max + dist_grad - 1e-6, dist_grad) + 360

        lat_true = np.logical_and(datos_asimilables[:,1]>lat_so[i], datos_asimilables[:,1]<=lat_so[i+1])

        for j in range(lon_so.size-1):
            #Busco donde se cumple la condicion de que haya datos para las lat-lon de las cajas
            lon_true = np.logical_and(datos_asimilables[:,0]>lon_so[j], datos_asimilables[:,0]<=lon_so[j+1])

            #Busco donde se cumplen ambas condiciones
            a = np.argwhere((lat_true == True)*(lon_true == True) == True)
            if len(a)>0:
                for k in range(pres_so.size-1):
                    #Busco en que nivel de presion estan los datos
                    pres_true = np.logical_and(datos_asimilables[:,2][a]>pres_so[k], datos_asimilables[:,2][a]<=pres_so[k+1])
                    b = np.argwhere(pres_true == True)[:,0]
                    if len(b)>0:
                        #Guardo los datos promediando los vientos y asignando el valor al punto central de la caja
                        datos_final[l:l+2,5] = 7.5
                        datos_final[l:l+2,6] = 4
                        for m in range(2):
                            datos_final[l, 0] = 2819 + m
                            datos_final[l, 1] = (lon_so[j+1]+lon_so[j])/2
                            datos_final[l, 2] = (lat_so[i+1]+lat_so[i])/2
                            datos_final[l, 3] = (pres_so[k+1] + pres_so[k])/2
                            datos_final[l, 4] = np.nanmean(datos_asimilables[:, 3+m][a][b])
                            datos_final[l, 7] = np.nanstd(datos_asimilables[:, 3+m][a][b])
                            l = l+1

    return datos_final


  

####
# Parametros a editar
####

ventana = int(os.environ['OBSFREC'])*60
lat_min = -60
lat_max = -5
lon_min = -90
lon_max = -40


anio = str(sys.argv[1])
mes = str(sys.argv[2])
dia = str(sys.argv[3])
hora = str(sys.argv[4])
minuto = str(sys.argv[5])


#Genero un vector con la fecha
dia_vect = datetime.datetime.strptime(anio + mes + dia + hora + minuto + '00', '%Y%m%d%H%M%S')

data_path = st.DATOS_DIR #Path del repositorio de datos
carpetas= ['/'] #Carpeta en la que buscar los datos

#print(dia_vect)

##Primero busco los datos que vienen por PDA, si no están actualizados busca por
##GNCA
#fuente = 'PDA'
#try:
#  for file in glob.glob(data_path+carpetas[0]+'OR_ABI-L2-*-M?C*_G16_*'):
#    fecha = datetime.datetime.strptime(file.split('/')[-1].split('_')[3][1:12], '%Y%j%H%M')
#    if fecha == dia_vect:
#      raise Exception
#  data_path = '/home/ftp/gnca/DMWF/' #Path a los datos de GNCA
#  carpetas = ['DMWF-C02/', 'DMWF-C07/', 'DMWF-C08/', 'DMWF-C09/', 'DMWF-C10/',
#              'DMWF-C14/', 'DMWVF-C08/'] #Carpetas en la que buscar los datos
#  fuente = 'GNCA'
#except:
#  pass



k = 0

#For para recorrer las carpetas
for carpeta in carpetas:
#Busco en la carpeta de los datos los archivos que tengan datos en la ventana de 10 minutos
#centrada en la hora a asimilar
	files = []
	for filename in glob.glob(data_path+carpeta+'OR_ABI-L2-*-M?C*_G16_*'):
		filenamesplit=filename.split('/')[-1].split('_')
		delta = np.abs(datetime.datetime.strptime(filenamesplit[3][1:14], '%Y%j%H%M%S')-dia_vect).total_seconds()
		if delta < ventana/2:
			files.append(filename)
			fecha_img = datetime.datetime.strptime(filenamesplit[3][1:14], '%Y%j%H%M%S')
			if filenamesplit[1][10] == 'V':
				modo = filenamesplit[1][14]
			else:
				modo = filenamesplit[1][13]
	for file in files: 
			#print(file)

			#Levanto los datos a asimilar
			#datos = Dataset(data_path + carpeta + file[0])
			datos = Dataset(file)

			presion = datos.variables['pressure'][:]
			direccion = datos.variables['wind_direction'][:]
			velocidad = datos.variables['wind_speed'][:]
			lat = datos.variables['lat'][:]
			lon = datos.variables['lon'][:]

			#Flag que contiene información de si el viento es bueno o no en una posición
			flag = datos.variables['DQF'][:]	  
			datos.close()

			#Me quedo con los datos en los que flag=0 que son los que el viento esta
			#bien estimado según el algoritmo
			ind_good = np.argwhere(flag==0)

			lat = np.squeeze(lat[ind_good])
			lon = np.squeeze(lon[ind_good])
			presion = np.squeeze(presion[ind_good])
			direccion = np.squeeze(direccion[ind_good])
			velocidad = np.squeeze(velocidad[ind_good])

			#Los datos son de todo el escaneo del GOES-16 entonces tomo solamente un sector
			ind_region = np.argwhere((lat>=lat_min) & (lat<=lat_max) & (lon>=lon_min) & (lon<=lon_max)) 

			lat = np.squeeze(lat[ind_region])
			lon = np.squeeze(lon[ind_region])
			presion = np.squeeze(presion[ind_region])
			direccion = np.squeeze(direccion[ind_region])
			velocidad = np.squeeze(velocidad[ind_region])

			#Si hay 1 solo dato en el dominio, la linea de abajo que calcula la componente 
			#del viento rompe. Por el momento descarto el archivo pero habria que buscar 
			#una solucion para conservar el dato.
			try:
			
				#Descompongo la magnitud del viento en sus componentes zonal y meridional
				u = -velocidad[:] * np.sin((np.pi/180) * direccion[:])
				v = -velocidad[:] * np.cos((np.pi/180) * direccion[:])

				#Genero una matriz que almacene los datos antes de realizar el superobbing
				ndatos = int(velocidad.shape[0])

				amv = np.full((ndatos, 5), np.nan)

				for i in range(ndatos):
					amv[i, 0] = lon[i]+360 
					amv[i, 1] = lat[i]
					amv[i, 2] = presion[i]
					amv[i, 3] = u[i]
					amv[i, 4] = v[i]
			except:
				#Descompongo la magnitud del viento en sus componentes zonal y meridional
				u = -velocidad * np.sin((np.pi/180) * direccion) 
				v = -velocidad * np.cos((np.pi/180) * direccion)
                                                                                                    
				#Genero una matriz que almacene los datos antes de realizar el superobbing
				ndatos = 1
                                 
				amv = np.full((ndatos, 5), np.nan) 
                                             
				for i in range(ndatos):   
					amv[i, 0] = lon+360
					amv[i, 1] = lat
					amv[i, 2] = presion
					amv[i, 3] = u
					amv[i, 4] = v

			if k == 0:
				datos_asimilables = amv
			else:
				datos_asimilables = np.concatenate((datos_asimilables, amv),axis=0)

			k = k+1
#  elif len(file)>1:
#    #Por si al buscar los archivos se detecta más de 1
#    print('Se está tratando de procesar más de un archivo por canal')
#    continue

#Si en los archivos leídos no hay dato o en la hora de asimilación no hay ningún
#archivo se sale del script
if k == 0:
  sys.exit()


#print('Motion vector from ' + fuente)

#Elimino los datos que tienen igual latitud y longitud
datos_asimilables = pd.DataFrame(datos_asimilables)
datos_asimilables = datos_asimilables.drop_duplicates(subset = [0, 1])
datos_asimilables = datos_asimilables.values

#plt.barbs(datos_asimilables[:, 0], datos_asimilables[:, 1], datos_asimilables[:, 3], datos_asimilables[:, 4], datos_asimilables[:, 2])



#Llamo a la función que realiza el superobbing
datos_final = superobbing(datos_asimilables)

#Borro NaN que podrian haber quedado
delete = np.argwhere(np.isnan(datos_final))[:, 0]
datos_final = np.delete(datos_final, delete, axis = 0)

##############
## GUARDADO ##
##############

#Paso los datos a float32 que es el formato que lee Fortran
datos_final = np.float32(datos_final)


#Genero la fecha de la imagen del medio. Para el scaneo del globo es 900 segundos despues
#de la hora de la primera


if (dia_vect-fecha_img).total_seconds() > 0:
  if modo == '3':
    fecha_asim = dia_vect + datetime.timedelta(0,600)
  elif modo == '4':
    fecha_asim = dia_vect
  elif modo == '6':
    fecha_asim = dia_vect + datetime.timedelta(0,600)
  else:
    print('El modo ' + modo + ' no esta incluido')
else:
  if modo == '3':
    fecha_asim = dia_vect + datetime.timedelta(0,1200)
  elif modo == '4':
    fecha_asim = dia_vect + datetime.timedelta(0,600)
  elif modo == '6':
    fecha_asim = dia_vect + datetime.timedelta(0,600)
  else:
    print('El modo ' + modo + ' no esta incluido')

#Nombre del archivo de salida
fileout = st.DATOS_OUT + st.SATWND + fecha_asim.strftime('%Y%m%d%H%M%S') + '.dat'

#print(fileout)

f = FortranFile(fileout, 'w')
for i in range(datos_final.shape[0]):
  f.write_record(datos_final[i,:7])

f.close()
  
 
