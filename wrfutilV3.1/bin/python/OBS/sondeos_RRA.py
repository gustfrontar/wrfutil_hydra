# -*- coding: utf-8 -*-
import settings as st
import numpy as np
import pandas as pd
import sys
import metpy.calc
from metpy.units import units
from scipy.io import FortranFile
from datetime import datetime
import os


anio = str(sys.argv[1])
mes = str(sys.argv[2])
dia = str(sys.argv[3])
hora = str(sys.argv[4])
minuto = str(sys.argv[5])

#Levanto el txt con los datos
fecha_asim = datetime.strptime(anio + mes + dia + hora + minuto + '00', '%Y%m%d%H%M%S')

filename = fecha_asim.strftime(st.DATOS_DIR + 'asm_temp_%Y%m%d%H%M.lst')

#filename = st.DATOS_DIR+'asm_temp_' + hora + '.lst'

#Si el archivo no existe se sale del script
if os.path.isfile(filename) == False:
  sys.exit()


#Leo el archivo
datos = np.genfromtxt(filename, dtype=float, delimiter = ',', usecols = (2, 3, 4, 5, 6, 7, 8, 9, 10))

fecha_data = np.genfromtxt(filename, dtype=str, delimiter = ',', usecols = (0, 1))

#Como el archivo con los datos se sobrescribe todos los días hago un check 
#para no leer datos de otro día que no sea el que se quiere asimilar

dia_data = datetime.strptime((fecha_data[0, 0] + fecha_data[0, 1]).replace(" ", ""), "%d%m%Y%H%M%S")

if fecha_asim != dia_data:
  print('La fecha de los datos no se corresponde con la ingresada')
  sys.exit()


#Le pongo nombre a las columnas
columna = ['Estacion', 'Latitud', 'Longitud', 'Presion', 
           'T', 'Td', 'Altura', 'Direccion', 'Velocidad']


#Lo paso a un dataframe de pandas
datos = pd.DataFrame(datos,columns = columna)

#Guardo el número de las estaciones que reportaron dato
estaciones = pd.DataFrame(datos['Estacion'])
estaciones = estaciones.drop_duplicates()


#Errores de las observaciones
P_error_viento = np.array([10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400, 450,
                           500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000])
error_viento = np.array([2.7, 2.7, 2.7, 2.1, 2.7, 2.7, 3.0, 3.3, 3.3, 3.3, 3.0,
                         2.8, 2.5, 2.3, 2.0, 1.8, 1.6, 1.4, 1.3, 1.1, 1.1, 1.1,
                         1.1, 1.1])
error_T = 1.0
error_RH = np.array([15, 10])


#For que recorre todas las estaciones
for i in range(estaciones.shape[0]):
 
  #Me quedo con los datos que pertenecen a la estación que indica el for
  datos_estacion = datos[datos['Estacion'] == int(estaciones.iloc[i])]

  #Elimino los datos que no reportaron presión dado que es un dato que se 
  #requiere para la asimilación
  datos_estacion = datos_estacion[~np.isnan(datos_estacion['Presion'])]

  #Elimino los datos con presión repetida
  datos_estacion = datos_estacion.drop_duplicates('Presion', keep='first')

  #Ordeno los datos de mayor a menor en función de la presión
  datos_estacion = datos_estacion.sort_values('Presion', ascending=False)

  #Genero un df solo con los datos donde la T no es NaN 
  data_T = datos_estacion[~np.isnan(datos_estacion['T'])]

  #Genero otro df solo con los datos donde el viento no es NaN
  data_viento = datos_estacion[(~np.isnan(datos_estacion['Velocidad']))]

  #Genero un último df con los datos donde la humedad no es NaN
  data_humedad = datos_estacion[(~np.isnan(datos_estacion['Td']))]

  presion_humedad = np.asarray(datos_estacion[(~np.isnan(datos_estacion['Td']))])[:,5]

  RH = metpy.calc.relative_humidity_from_dewpoint(np.asarray(data_humedad)[:,5]*units.celsius,
                 np.asarray(data_humedad)[:,6]*units.celsius)*100

  #Creo una matriz para almacenar los datos a asimilar con el largo de la cantidad
  #de datos de temperatura + la cantidad de datos de viento*2 (u y v)
  datos_asimilables = np.full((data_T.shape[0] + data_viento.shape[0]*2 + 
                               data_humedad.shape[0],7), np.nan)

  #Creo un índice para ir llenando la matriz
  indice = 0

  #Cargo los datos de temperatura primero
  for j in range(data_T.shape[0]):
    datos_asimilables[indice, 0] = 3073 #Tipo de variable
    datos_asimilables[indice, 1] = data_T['Longitud'].iloc[j] #Longitud del dato
    datos_asimilables[indice, 2] = data_T['Latitud'].iloc[j] #Latitud del dato
    datos_asimilables[indice, 3] = data_T['Presion'].iloc[j] #Presion del dato
    datos_asimilables[indice, 4] = data_T['T'].iloc[j] + 273.15 #Dato
    datos_asimilables[indice, 5] = error_T #Error del dato
    datos_asimilables[indice, 6] = 1 #Tipo de dato
    indice = indice + 1

  #Cargo los datos de humedad
  for j in range(data_humedad.shape[0]):
    datos_asimilables[indice, 0] = 3331 #Tipo de variable
    datos_asimilables[indice, 1] = data_humedad['Longitud'].iloc[j] #Longitud del dato
    datos_asimilables[indice, 2] = data_humedad['Latitud'].iloc[j] #Latitud del dato
    datos_asimilables[indice, 3] = data_humedad['Presion'].iloc[j] #Presion del dato
    datos_asimilables[indice, 4] = RH[j] #Dato 
    P = np.argmin(np.abs(float(presion_humedad[j])-np.array([1000, 850])))
    datos_asimilables[indice, 5] = error_RH[P] #Error del dato
    datos_asimilables[indice, 6] = 1 #Tipo de dato
    indice = indice + 1


  #Descompongo el viento es sus componentes zonal y meridional y lo paso de nudo a m/s
  u = np.array(-np.sin(data_viento['Direccion']*np.pi/180)*data_viento['Velocidad']) * 0.514444
  v = np.array(-np.cos(data_viento['Direccion']*np.pi/180)*data_viento['Velocidad']) * 0.514444

  #Cargo los datos de viento
  for j in range(data_viento.shape[0]):

    datos_asimilables[indice:indice+2, 1] = data_viento['Longitud'].iloc[j] #Longitud del dato
    datos_asimilables[indice:indice+2, 2] = data_viento['Latitud'].iloc[j] #Latitud del dato
    datos_asimilables[indice:indice+2, 3] = data_viento['Presion'].iloc[j] #Presion del dato
    P = np.argmin(np.abs(float(data_viento['Presion'].iloc[j])-P_error_viento))
    datos_asimilables[indice:indice+2, 5] = error_viento[P] #Error del dato
    datos_asimilables[indice:indice+2, 6] = 1 #Tipo de dato
    
    for k in range(2):
      if k == 0:
        datos_asimilables[indice, 0] = 2819 #Tipo de variable
        datos_asimilables[indice, 4] = u[j] #Dato
      elif k == 1:
        datos_asimilables[indice, 0] = 2820 #Tipo de variable
        datos_asimilables[indice, 4] = v[j] #Dato
      indice = indice + 1


  if i == 0:
    datos_final = datos_asimilables
  else:
    datos_final = np.concatenate((datos_final, datos_asimilables),axis=0)


#Borro NaN que podrian haber quedado
delete = np.argwhere(np.isnan(datos_final))[:, 0]
datos_final = np.delete(datos_final, delete, axis = 0)

##############
## GUARDADO ##
##############

#Paso los datos a float32 que es el formato que lee Fortran
datos_final = np.float32(datos_final)

#fileout = st.DATOS_OUT + st.ADPUPA + anio + mes + dia + hora + minuto + '00.dat'
fileout = fecha_asim.strftime(st.DATOS_OUT + st.ADPUPA + '%Y%m%d%H%M00.dat')

f = FortranFile(fileout, 'w')
for i in range(datos_final.shape[0]):
  f.write_record(datos_final[i,:])

f.close()





