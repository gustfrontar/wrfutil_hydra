# -*- coding: utf-8 -*-
"""
Script para abrir archivos con datos obs y guardar un binario para asimilar
"""
import settings as st
import numpy as np
import pandas as pd
import metpy.calc as mcalc
from metpy.units import units
import math
import sys
import os.path
from scipy.io import FortranFile
from datetime import datetime


#Errores asociados a la medición
errores = {'viento':1.4, 'temperatura':2.0,
           'humedadRelativa':10.0, 'presion':1.6}

###############################################################################

def calculo_vientoComponentes(intensidad,direccion):
    """ Calculo de las componentes u y v
    #####################################
    Keyword arguments:
    intensidad: Magnitud del viento en nudos
    direccion: Direccion del viento en grados
    """
    u = -intensidad*math.sin(math.radians(direccion))*0.514444 #Paso de nudos a m/s
    v = -intensidad*math.cos(math.radians(direccion))*0.514444 #Paso de nudos a m/s
    if abs(u) < 0.0001:
        u = 0.0
    if abs(v) < 0.0001:
        v = 0.0

    return u,v


# DEFINO ALGUNAS COSAS ANTES DE EMPEZAR
tipo = {'ADPUPA':1,'AIRCFT':3,'SATWND':4,'SATEMP':7,'ADPSFC':8,'SFCSHP':9,
'SFCBOG':10,'SYNDAT':12,'ASCATW':20,'AIRS':21}

elemento = {'vientoZonal':2819,'vientoMeridional':2820,'temperatura':3073,
'humedadEspecifica':3330,'humedadRelativa':3331,'temperaturaVirtual':3079,
'presion':14593}


###############################################################################
# 	                          APERTURA DE ARCHIVOS
###############################################################################

anio = str(sys.argv[1])
mes = str(sys.argv[2])
dia = str(sys.argv[3])
hora = str(sys.argv[4])
minuto = str(sys.argv[5])

fecha_asim = datetime.strptime(anio + mes + dia + hora + minuto + '00', '%Y%m%d%H%M%S')

filename = fecha_asim.strftime(st.DATOS_DIR + 'asm_boyas_%Y%m%d%H%M.lst')

#filename = st.DATOS_DIR + 'asm_boyas_' + hora + '.lst'

#Si el archivo no existe se sale del script
if os.path.isfile(filename) == False:
  sys.exit()



columns = ['fecha','hora','boya','latitud','longitud','temperatura',
         'temperaturaRocio','presionEstacion','presionNMM',
         'direccionViento','intensidadViento','humedadRelativa', '9']


data = pd.read_csv(filename, names = columns)

#Como el archivo contiene datos de las últimas 3 horas trabajo solo con
#los datos de la hora que se quiere asimilar
data = data[(data['hora'] == int(hora))]

if data.size==0:
  sys.exit()

#Como el archivo se va sobrescribiendo todos los días chequeo que los datos
#no sean de otro día
fecha_data = str(data['fecha'].iloc[0])

anio_data = fecha_data[:4]
mes_data = fecha_data[4:6]
dia_data = fecha_data[6:8]

if (anio!=anio_data or mes!=mes_data or dia!=dia_data):
  sys.exit()

data['temperatura'] = data['temperatura'] + 273.15


data = data.dropna(subset=['latitud', 'longitud'])

data = data.drop_duplicates()


datos_asimilacion = np.full([data.shape[0]*5,7], np.nan)

variables = ['viento', 'temperatura', 'presion', 'humedadRelativa']
i = 0
for index in data.index:
  for var in variables:
    if var == 'viento':
      intensidad = data['intensidadViento'].loc[index]
      direccion = data['direccionViento'].loc[index]

      #Si los datos asociados al viento son NaN no lo guardo
      if np.logical_and(np.isnan(direccion), np.isnan(intensidad)):
          continue

      u = np.nanmean(-intensidad*np.sin(np.deg2rad(direccion)))*0.514444
      v = np.nanmean(-intensidad*np.cos(np.deg2rad(direccion)))*0.514444


      if abs(u) < 0.1:
          u = 0.0
      if abs(v) < 0.1:
          v = 0.0

      datos_asimilacion[i,0] = elemento['vientoZonal']      #Variable
      datos_asimilacion[i,1] = data['longitud'].loc[index]  #Longitud
      datos_asimilacion[i,2] = data['latitud'].loc[index]   #Latitud
      datos_asimilacion[i,3] = 0                            #Altura
      datos_asimilacion[i,4] = np.round(u,1)                #Valor
      datos_asimilacion[i,5] = errores[var]                 #Error
      datos_asimilacion[i,6] = 9                            #Tipo de dato
      i = i+1

      datos_asimilacion[i,0] = elemento['vientoMeridional'] #Variable
      datos_asimilacion[i,1] = data['longitud'].loc[index]  #Longitud
      datos_asimilacion[i,2] = data['latitud'].loc[index]   #Latitud
      datos_asimilacion[i,3] = 0                            #Altura
      datos_asimilacion[i,4] = np.round(v,1)                #Valor
      datos_asimilacion[i,5] = errores[var]                 #Error
      datos_asimilacion[i,6] = 9                            #Tipo de dato
      i = i+1

    elif var == 'presion':

      if ~np.isnan(data['presionNMM'].loc[index]):
        tipo_pres = 'presionNMM'
      elif ~np.isnan(data['presionEstacion'].loc[index]):
        tipo_pres = 'presionEstacion'
      else:
        continue  

      datos_asimilacion[i,0] = elemento[var]                   #Variable
      datos_asimilacion[i,1] = data['longitud'].loc[index]     #Longitud
      datos_asimilacion[i,2] = data['latitud'].loc[index]      #Latitud
      datos_asimilacion[i,3] = 0                               #Altura
      datos_asimilacion[i,4] = np.round(data[tipo_pres].loc[index],1)
      datos_asimilacion[i,5] = errores[var]                    #Error
      datos_asimilacion[i,6] = 9                               #Tipo de dato
      i = i+1

    else:
      if np.isnan(data[var].loc[index]):
        continue

      datos_asimilacion[i,0] = elemento[var]                   #Variable
      datos_asimilacion[i,1] = data['longitud'].loc[index]     #Longitud
      datos_asimilacion[i,2] = data['latitud'].loc[index]      #Latitud
      datos_asimilacion[i,3] = 0                               #Altura
      datos_asimilacion[i,4] = np.round(data[var].loc[index],1)  #Valor
      datos_asimilacion[i,5] = errores[var]                    #Error
      datos_asimilacion[i,6] = 9                               #Tipo de dato
      i = i+1


#Elimino las filas de NaN que pudieran haber quedado 
delete = np.argwhere(np.isnan(datos_asimilacion))[:, 0]
datos_asimilacion = np.delete(datos_asimilacion, delete, axis = 0)


##############
## GUARDADO ##
##############


# Paso los datos a float 32 porque es el formato que lee Fortran
datos_asimilacion=np.float32(datos_asimilacion)

#fileout = st.DATOS_OUT + st.BOYAS + anio + mes + dia + hora + minuto + '00.dat'
fileout = fecha_asim.strftime(st.DATOS_OUT + st.BOYAS + '%Y%m%d%H%M00.dat')

f = FortranFile(fileout, 'w')
for i in range(datos_asimilacion.shape[0]):
  f.write_record(datos_asimilacion[i,:])

f.close()


