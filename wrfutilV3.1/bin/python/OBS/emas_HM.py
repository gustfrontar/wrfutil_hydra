# -*- coding: utf-8 -*-
import settings as st
import numpy as np
import pandas as pd
import sys
from datetime import datetime,timedelta
import glob
from scipy.io import FortranFile
import os

anio = sys.argv[1]
mes = sys.argv[2]
dia = sys.argv[3]
hora = sys.argv[4]
minuto = sys.argv[5]

#Frecuencia de las observaciones en segundos
ventana = int(os.environ['OBSFREC'])*60

#Erroes de cada dato
errores = {'humedad':10.0, 'presion':1.0, 'temperatura':2.0, 'viento':1.4}

#Tipo de dato
elemento = {'vientoZonal':82819,'vientoMeridional':82820,'temperatura':83073,
'humedadEspecifica':83330,'humedad':83331,'temperaturaVirtual':83079,
'presion':14593}


dia_vect = datetime.strptime(anio + mes + dia + hora + minuto + '00', '%Y%m%d%H%M%S')

#DataFrame que contiene la ubicacion de cada estacion
estaciones = pd.read_json(st.DATOS_DIR + 'estaciones.json')
estaciones = estaciones.set_index('idEstacion')

#Elimino las estaciones de las que no se tienen datos importantes
estaciones = estaciones.dropna(axis=0, subset=['altura', 'latitud', 'longitud'])


files = glob.glob(st.DATOS_DIR + 'prop*')

data = pd.DataFrame()
#DataFrame con las observaciones
for file in files:
    data = data.append(pd.read_json(file), sort = False)

#Paso el dato de fecha y hora de string a tiempo
data['fechaHora'] = pd.to_datetime(data['fechaHora'], format = '%Y-%m-%dT%H:%M:%S')

#Defino el indice de los DataFrame con el id de la estacion
data = data.set_index('idEstacion')


#Defino las hora de inicio y fin de la ventana de asimialcion
hh_inicio = datetime.strptime(anio + mes + dia + hora + minuto + '00', '%Y%m%d%H%M%S') - timedelta(seconds = ventana)
hh_final = datetime.strptime(anio + mes + dia + hora + minuto + '00', '%Y%m%d%H%M%S') + timedelta(seconds = ventana)

#Me quedo con los datos que estan dentro de la ventana de asimilacion
data = data.loc[np.logical_and(data['fechaHora'] >= hh_inicio, data['fechaHora'] < hh_final)]

#Paso el dato de temperatura a kelvin
data['temperatura'] = data['temperatura'] + 273.15

estaciones['longitud'] = estaciones['longitud'] + 360.

#De cada estacion saco el valor de U, V, HR, P y T entonces creo una array 
#de dimension cantidad de estaciones* 5 para guardar los datos a asimilar
datos_asimilacion = np.full([data.shape[0]*5,7], np.nan)


#Nombres de las variables que quiero usar para asimilar
variables = ['humedad', 'presion', 'temperatura', 'viento']
i = 0

#For que recorre a los propietarios de las estaciones
for p in range(1, len(files) + 1):
    #Genero 2 DataFrame unicamente con los datos de un propietario
    data_prop = data[data['idPropietario'] == p]
    estaciones_prop = estaciones[estaciones['idPropietario'] == p]

    #For para las estaciones de cada propietario
    for index in data_prop.index.drop_duplicates():

        #Antes habia eliminado estaciones sin datos importantes, entonces chequeo que no haya eliminado la estacion
        if index not in estaciones_prop.index:
            continue

        #For para cada variable a asimilar
        for var in variables:
            if var == 'viento': 
                intensidad = data_prop.loc[(index), 'velViento']
                direccion = data_prop.loc[(index), 'dirViento']

                #Si los datos asociados al viento son NaN no lo guardo
                if np.isnan(intensidad).all() or np.isnan(direccion).all():
                    continue

                u = np.nanmean(-intensidad*np.sin(np.deg2rad(direccion))/3.6)
                v = np.nanmean(-intensidad*np.cos(np.deg2rad(direccion))/3.6)
    
                if abs(u) < 0.1:
                    u = 0.0
                if abs(v) < 0.1:
                    v = 0.0

                datos_asimilacion[i,0] = elemento['vientoZonal']
                datos_asimilacion[i,1] = estaciones_prop.loc[(index), 'longitud']
                datos_asimilacion[i,2] = estaciones_prop.loc[(index), 'latitud']
                datos_asimilacion[i,3] = estaciones_prop.loc[(index), 'altura']
                datos_asimilacion[i,4] = np.round(u,1)
                datos_asimilacion[i,5] = errores[var]
                datos_asimilacion[i,6] = 22
                i = i+1

                datos_asimilacion[i,0] = elemento['vientoMeridional']
                datos_asimilacion[i,1] = estaciones_prop.loc[(index), 'longitud']
                datos_asimilacion[i,2] = estaciones_prop.loc[(index), 'latitud']
                datos_asimilacion[i,3] = estaciones_prop.loc[(index), 'altura']
                datos_asimilacion[i,4] = np.round(v,1)
                datos_asimilacion[i,5] = errores[var]
                datos_asimilacion[i,6] = 22
                i = i+1            
            else:
                #Si el dato a asimilar es NaN no lo guardo
                if np.isnan(data_prop.loc[(index), var]).all():
                    continue

                datos_asimilacion[i,0] = elemento[var]
                datos_asimilacion[i,1] = estaciones_prop.loc[(index), 'longitud']
                datos_asimilacion[i,2] = estaciones_prop.loc[(index), 'latitud']
                datos_asimilacion[i,3] = estaciones_prop.loc[(index), 'altura']
                datos_asimilacion[i,4] = np.round(np.nanmean(data_prop.loc[(index), var]), 1)
                datos_asimilacion[i,5] = errores[var]
                datos_asimilacion[i,6] = 22
                i = i+1


#Borro NaN que podrian haber quedado
delete = np.argwhere(np.isnan(datos_asimilacion))[:, 0]
datos_asimilacion = np.delete(datos_asimilacion, delete, axis = 0)


##########
# Guardo #
##########

# Paso a float32 porque el es formato que lee Fortran 
datos_asimilacion=np.float32(datos_asimilacion)

fileout = st.DATOS_OUT + st.ADPAUT + anio + mes + dia + hora + minuto + '00.dat'

f = FortranFile(fileout, 'w')
for i in range(datos_asimilacion.shape[0]):
  f.write_record(datos_asimilacion[i,:])

f.close()

