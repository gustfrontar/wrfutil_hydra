# -*- coding: utf-8 -*-
import settings as st
import numpy as np
from scipy.io import FortranFile
import sys
import os
import glob
from datetime import datetime, timedelta


def amdar_dpa(fecha):

    filename = fecha.strftime('{}/amdar_%d%m%Y_{}_%Mz'.format(st.DATOS_DIR, fecha.hour))

    #Si el archivo no existe salgo del script
    if os.path.isfile(filename) == False:
        return np.full(7, None)

    # Leo los datos de latitud, longitud, presion y temperatura
    datanum = np.genfromtxt(filename, skip_header = 2, usecols = (1, 2, 6, 7), encoding="latin-1")
    if datanum.size == 0 :
        return np.full(7, None)


    # Leo los datos de viento(°/kts) como string
    datastr = np.genfromtxt(filename, dtype = 'str', skip_header = 2, usecols = 8,encoding="latin-1")

    if len(datanum.shape) == 1:
        datanum = datanum[np.newaxis, :] 
        datastr = datastr[np.newaxis]


    #Defino una lista con el número asociado a cada medicion que aporta el AMDAR
    elemento = [2819, 2820, 3073]

    tipo = 3

    #Errores en las mediciones
    error_viento = 3.6
    error_temperatura = 1

    datos_amdar = np.full((datanum.shape[0]*3,7),np.nan) # El fila*3 viene de que para cada dato recibido se
                                                         # requiere guardar la temperatura y el viento zonal
                                                         # y meridional 

    index = 0

    i = 0

    while i < datanum.shape[0]:
        if datastr[i] == '///////':
            # Si el avion se encuentra en superficie no se disponen de datos de viento entonces
            # se almacenan solo datos de temperatura
            datos_amdar[index,0] = elemento[2]  #Defino el tipo de elemento
            datos_amdar[index,1] = datanum[i,1]+360 #Longitud
            datos_amdar[index,2] = datanum[i,0] #Latitud     
            datos_amdar[index,3] = datanum[i,2] #Presion
            datos_amdar[index,4] = datanum[i,3] + 273.15 #Tempteratura
            datos_amdar[index,5] = error_temperatura #Error de la temperatura
            datos_amdar[index,6] = tipo #Tipo de dato
            index = index+1
        else:
            # Si se disponen de datos de viento, almaceno entonces los 3 datos disponibles
            viento = float(datastr[i][4:7])*0.514444 #En los datos el viento esta en kts => lo paso a m/s
            direccion = np.deg2rad(float(datastr[i][0:3]))
            for j in range(3):
                datos_amdar[index,0] = elemento[j] #Elemento
                datos_amdar[index,1] = datanum[i,1]+360 #Longitud
                datos_amdar[index,2] = datanum[i,0] #Latitud
                datos_amdar[index,3] = datanum[i,2] #Presion
                if j == 0:
                    datos_amdar[index,4] = -np.sin(direccion)*viento #Viento zonal
                    datos_amdar[index,5] = error_viento #Error del viento zonal
                if j == 1:
                    datos_amdar[index,4] = -np.cos(direccion)*viento #Viento meridional
                    datos_amdar[index,5] = error_viento #Error del viento emridional
                if j == 2:
                    datos_amdar[index,4] = datanum[i,3] + 273.15 #Temperatura
                    datos_amdar[index,5] = error_temperatura #Error de la temperatura
                datos_amdar[index,6] = tipo #Tipo de dato
                index = index+1

        repetido = np.argwhere(np.sum((datanum[i,:] == datanum), axis = 1)==4)
        if len(repetido)>1:
            datanum = np.delete(datanum, repetido[1:], axis = 0)
            datastr = np.delete(datastr, repetido[1:], axis = 0)

        i += 1

    return datos_amdar


###################################
### AMDAR Aerolineas Argentinas ###
###################################


def file_correcto(file, ini, fin):
    tmp = file.split('_')[-1]
    res = []
    if (tmp>=ini and tmp<=fin):
        res.append(file)

    return res


def FL2pres(FL):
    P0 = 1013.25
    h = FL * 100 * 0.3048
    P = np.round(P0*(1 - 0.0000226 * h)**5.255)

    return P

def parser(line, ini, fin):

    elemento = [3073, 2819, 2820]
    tipo = 3
    error_viento = 3.6
    error_temperatura = 1

    #No pongo separar por split(' ') porque hay veces que hay mas de un espacio
    #entre variables y hace mas complicado todo
    data = line.split()

    #j es un indice para agarrar los casos en que hay algun espacion en blanco
    #de mas en el archivo
    #n es un indice para contar cuantas variables se extrajeron de la linea
    j = 0
    n = 0

    #A veces aparecen datos de otros tiempos que no se corresponden con los del
    #nombre del archivo entonces chequeo que los datos caigan en la ventana de 
    #asimilacion
    if fin.strftime('%d%H%M') < data[4] or ini.strftime('%d%H%M') > data[4]:
        return np.full(7, None)


    #Latitud
    try:
        signo = -1
        if data[2][-1] == 'N':
            signo = 1
        lat = int(data[2][:-1])/100 * signo
    except:
        return np.full(7, None)

    #Longitud
    try:
        if data[3][-1] == 'W':
            lon = -int(data[3][:-1])/100 + 360
        elif data[3][-1] == 'E':
            lon = int(data[3][:-1])/100
        else:
            lon = -int(data[3])/1000 + 360
    except:
        return np.full(7, None)

    #Nivel de vuelo => presion
    try:
        if data[5] == 'F':
            FL = int(data[6])
            j = j + 1
        else:
            FL = int(data[5][1:])

        presion = FL2pres(FL)
    except:
        return np.full(7, None)

    #Temperatura
    try:
        signo = 1
        if data[6+j] == 'PS' or data[6+j] == 'MS':
            if data[6+j] == 'MS':
                signo = -1
            tmp = int(data[7+j])/10
            j = j + 1
        else:
            tmp = int(data[6+j][2:])/10
            if data[6+j][:2] == 'MS':
                signo = -1

        temperatura = tmp * signo + 273.15
        n += 1
    except:
        pass

    #Viento
    try:
        tmp = data[7+j].split('=')[0].split('/')
        viento = int(tmp[1])*0.514444 #En los datos el viento esta en kts => lo paso a m/s
        if tmp[0][0] == ',':
            direccion = np.deg2rad(int(tmp[0][1:]))
        else:
            direccion = np.deg2rad(int(tmp[0]))

        u = -np.sin(direccion)*viento
        v = -np.cos(direccion)*viento
        n += 2
    except:
        pass

    #Si n es 0 es que no se pudieron obtener datos de temperatura ni viento
    if n == 0:
        return np.full(7, None)

    dato_amdar = np.full([n, 7], np.nan)

    #Datos que comprarten las 3 variables a asimilar
    for i in range(n):
        dato_amdar[i, 0] = elemento[i]
        dato_amdar[i, 1] = lon
        dato_amdar[i, 2] = lat
        dato_amdar[i, 3] = presion
        dato_amdar[i, 6] = tipo

    if n == 1: #Solo datos de temperatura
        dato_amdar[0, 4] = temperatura
        dato_amdar[0, 5] = error_temperatura
    elif n == 2: #Solo datos de viento
        dato_amdar[0, 4] = u
        dato_amdar[1, 4] = v
        dato_amdar[:, 5] = error_viento
    elif n == 3: #Datos de temperatura y viento
        dato_amdar[0, 4] = temperatura
        dato_amdar[0, 5] = error_temperatura
        dato_amdar[1, 4] = u
        dato_amdar[2, 4] = v
        dato_amdar[1:3, 5] = error_viento


    return dato_amdar



def amdar_AA(fecha):

    path = st.DATOS_DIR + '/AA_AMDAR/'

    obsfrec = int(os.environ['OBSFREC'])

    #Inicio y fin de la ventana
    ini = (fecha - timedelta(minutes = obsfrec/2))
    fin = (fecha + timedelta(minutes = obsfrec/2))

    ini_str = ini.strftime('%Y%m%d%H%M.txt')
    fin_str = fin.strftime('%Y%m%d%H%M.txt')


    #Listo todos los archivos disponibles
    files = glob.glob(fecha.strftime(path + '/*'))

    #Listo solo los archivos cuya fecha cae dentro de la ventana de asimilacion
    archivos = []
    for f in files:
        archivos += file_correcto(f, ini_str, fin_str)

    #Opera sobre los archivos para acomodar los datos a como los necesita el LETKF
    i = 0
    for file in archivos:
        with open(file, 'r') as f:
            for line in f.readlines():
                datos = parser(line, ini, fin)
                if (datos == None).all():
                    continue

                if i == 0:
                    datos_final = datos
                else:
                    datos_final = np.concatenate([datos_final, datos])
                i += 1

    if i == 0:
        return np.full(7, None)
        
    return datos_final





anio = int(sys.argv[1])
mes = int(sys.argv[2])
dia = int(sys.argv[3])
hora = int(sys.argv[4])
minuto = int(sys.argv[5])

fecha_asim = datetime(anio, mes, dia, hora, minuto)


datos_dpa = amdar_dpa(fecha_asim)
datos_AA = amdar_AA(fecha_asim)

if (datos_dpa != None).all() and (datos_AA != None).all():
    datos_amdar = np.concatenate([datos_dpa, datos_AA])
elif (datos_AA != None).all() and not (datos_dpa != None).all():
    datos_amdar = datos_AA
elif (datos_dpa != None).all() and not (datos_AA != None).all():
    datos_amdar = datos_dpa
else:
    sys.exit()

#Borro NaN que podrian haber quedado
delete = np.argwhere(np.isnan(datos_amdar))[:, 0]
datos_amdar = np.delete(datos_amdar, delete, axis = 0)

##############
## GUARDADO ##
##############

#Paso los datos a float32 que es el que lee Fortran
datos_amdar = np.float32(datos_amdar)

fileout = fecha_asim.strftime(st.DATOS_OUT + st.AIRCFT + '%Y%m%d%H%M00.dat')

f = FortranFile(fileout, 'w')
for i in range(datos_amdar.shape[0]):
  f.write_record(datos_amdar[i,:])

f.close()

