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
import os
from scipy.io import FortranFile
from datetime import datetime

#Errores asociados a la medición
errores = {'viento':1.4, 'temperatura':2.0,
           'humedadRelativa':10.0, 'presionNMM':1.6}

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

elemento = {'vientoZonal':82819,'vientoMeridional':82820,'temperatura':83073,
'humedadEspecifica':83330,'humedadRelativa':83331,'temperaturaVirtual':83079,
'presionNMM':14593}


###############################################################################
# 	                          APERTURA DE ARCHIVOS
###############################################################################

anio = str(sys.argv[1])
mes = str(sys.argv[2])
dia = str(sys.argv[3])
hora = str(sys.argv[4])
minuto = str(sys.argv[5])

fecha_asim = datetime.strptime(anio + mes + dia + hora + minuto + '00', '%Y%m%d%H%M%S')

filename = fecha_asim.strftime(st.DATOS_DIR + 'asm_ship_%Y%m%d%H%M.lst')

#filename = HOME + 'asm_ship_' + hora + '.lst'

#Si el archivo no existe se sale del script
if os.path.isfile(filename) == False:
  sys.exit()

#Si el archivo esta vacio salgo
if os.stat(filename).st_size == 0:
  sys.exit()


columns = ['fecha','hora','barco','latitud','longitud','temperatura',
         'temperaturaRocio','presionNMM','direccionViento','intensidadViento', '9']
#         'humedadRelativa', '9']


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
data['temperaturaRocio'] = data['temperaturaRocio'] + 273.15


data = data.dropna(subset=['latitud', 'longitud'])


datos_asimilacion = np.full([data.shape[0]*5,7], np.nan)

variables = ['viento', 'temperatura', 'presionNMM', 'temperaturaRocio']
i = 0
for index in data.index:
  for var in variables:
    if var == 'viento':
      intensidad = data['intensidadViento'].loc[index]
      direccion = data['direccionViento'].loc[index]

      #Si los datos asociados al viento son NaN no lo guardo
      if np.logical_or(np.isnan(direccion), np.isnan(intensidad)):
          continue

      u = np.nanmean(-intensidad*np.sin(np.deg2rad(direccion)))*0.514444
      v = np.nanmean(-intensidad*np.cos(np.deg2rad(direccion)))*0.514444


      if abs(u) < 0.1:
          u = 0.0
      if abs(v) < 0.1:
          v = 0.0

      datos_asimilacion[i,0] = elemento['vientoZonal']        #Variable
      datos_asimilacion[i,1] = data['longitud'].loc[index]    #Longitud
      datos_asimilacion[i,2] = data['latitud'].loc[index]     #Latitud
      datos_asimilacion[i,3] = 0                              #Altura
      datos_asimilacion[i,4] = np.round(u,1)                  #Valor
      datos_asimilacion[i,5] = errores[var]                   #Error
      datos_asimilacion[i,6] = 9                              #Tipo de dato
      i = i+1

      datos_asimilacion[i,0] = elemento['vientoMeridional']   #Variable
      datos_asimilacion[i,1] = data['longitud'].loc[index]    #Longitud
      datos_asimilacion[i,2] = data['latitud'].loc[index]     #Latitud
      datos_asimilacion[i,3] = 0                              #Altura
      datos_asimilacion[i,4] = np.round(v,1)                  #Valor
      datos_asimilacion[i,5] = errores[var]                   #Error
      datos_asimilacion[i,6] = 9                              #Tipo de dato
      i = i+1

    else:
      if np.isnan(data[var].loc[index]):
        continue
      datos_asimilacion[i,1] = data['longitud'].loc[index]     #Longitud
      datos_asimilacion[i,2] = data['latitud'].loc[index]      #Latitud
      datos_asimilacion[i,3] = 0                               #Altura

      if var == 'temperaturaRocio':
        datos_asimilacion[i,0] = elemento['humedadRelativa']   #Variable
        HR = mcalc.relative_humidity_from_dewpoint(data['temperatura'].loc[index]*units.K, data['temperaturaRocio'].loc[index]*units.K)
        datos_asimilacion[i,4] = np.round(HR*100,1)              #Valor
        datos_asimilacion[i,5] = errores['humedadRelativa']    #Error
      else:
        datos_asimilacion[i,0] = elemento[var]                 #Variable
        datos_asimilacion[i,4] = data[var].loc[index]          #Valor
        datos_asimilacion[i,5] = errores[var]                  #Error

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

#fileout = st.DATOS_OUT + st.SFCSHP + anio + mes + dia + hora + minuto + '00.dat'
fileout = fecha_asim.strftime(st.DATOS_OUT + st.SFCSHP + '%Y%m%d%H%M00.dat')

f = FortranFile(fileout, 'w')
for i in range(datos_asimilacion.shape[0]):
  f.write_record(datos_asimilacion[i,:])

f.close()


