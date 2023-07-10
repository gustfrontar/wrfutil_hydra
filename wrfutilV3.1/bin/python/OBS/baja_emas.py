import pandas as pd
import sys
import os
import urllib.request
import socket

socket.setdefaulttimeout(30) #Defino un timeout en segundos para la descarga de archivos

#Fechas iniciales y finales de las observaciones a descargar
obs_ini = sys.argv[1]
obs_fin = sys.argv[2]

#Directorio donde se van a descargar los archivos
DIROBSREPO = os.environ['DIROBSREPO'] + 'EMAs/'

#Descargo la lista de propietarios de estaciones
try:
    propietarios = pd.read_json('http://192.168.5.213:8080/aws-api/propietaries')
except:
    print('Fallo la descarga de los propietarios de las EMAs')
    sys.exit()

#IDs de los propietarios
IDs = propietarios['id'].values

#Descargo los datos de cada propietario
for id in IDs:
    filein = 'http://192.168.5.213:8080/aws-api/observations/{}/{}/{}'.format(id, obs_ini, obs_fin)
    fileout = DIROBSREPO + 'prop{}.json'.format(id)
    try:
        urllib.request.urlretrieve(filein, fileout)
    except:
        continue


#Leo el listado de estaciones y lo guardo
filein = 'http://192.168.5.213:8080/aws-api/stations/'
fileout = DIROBSREPO + 'estaciones.json'
try:
    urllib.request.urlretrieve(filein, fileout)
except:
    print('Fallo la descarga de los Id de las EMAs')
    sys.exit()

