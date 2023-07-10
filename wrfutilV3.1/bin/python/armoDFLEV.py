#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03feb2020 ==> ESTE SCRIPT ARMA EL DFwrf PARA LA NUBOSIDAD 

Correr este script en el entorno verif haciendo: 
source activate verif
""" 
 #VER los valores DE LOS FALTANTES!
# ===========================================================================
# Paso 0: LIMPIO ENTORNO Y CARGO LIBRERIAS
# ===========================================================================

# Cargo librerias
import argparse
from time import time
import numpy as np
import os
import time
import pandas as pd
from netCDF4 import Dataset
import datetime
from datetime import date, timedelta
import auxiliary_functions_DF as af
import glob

def main():

   FECHA_INI = os.environ['FECHA_INI']
   CICLO = os.environ['CICLO']
   INPDIR = os.environ['DIRIN']
   OUTBASEDIR = os.environ['OUTDIR_DF']
   BASEDIR = os.environ['BASEDIR']
   PLAZOS = int(os.environ['PLAZO'])
   VARIABLES = eval(os.environ['VARIABLES'])
   NOMBRE = os.environ['NOMBRE']

   WRFUTILDIR = os.environ['WRFUTILDIR']
   DIR_TEMPLATES = WRFUTILDIR + "/templates/"

   modelo = os.environ['FILEIN_TYPE']
    
   FECHA_INI_obj = datetime.datetime.strptime(FECHA_INI+CICLO,'%Y/%m/%d%H')
   FECHA_INI = FECHA_INI_obj.strftime('%Y%m%d')
   CICLO = FECHA_INI_obj.strftime('%H')
   
   if modelo == 'WRF':
       FILEIN = glob.glob(BASEDIR + '/WRF/*/WRFOUT/wrfout_d01_*')[0]
   else:
       print("Modelos {} desconocido !!".format(modelo))
       exit(10)
   
   OUTDIR = OUTBASEDIR + FECHA_INI_obj.strftime('%Y%m%d_%H') + '0000/'+NOMBRE+'/'

   Puntos_interes = af.merge_Puntos_Interes(DIR_TEMPLATES, OUTDIR, FILEIN, modelo, estimo_ij = True)

   cols = Puntos_interes[['CIPI', 'tipo', 'ID']]
   cols.columns = ['CIPI', 'Tipo', 'ID']
   cols.set_index('CIPI', inplace = True)

   CIPI = np.sort(Puntos_interes['CIPI'].unique()) 
  
   # ===========================================================================
   # PASO 1: GENERACION DEL DATAFRAME LLENO DE NANS CON SUS INDICES ASOCIADOS
   # ===========================================================================

# 03feb abro 1 archivo para sacar los niveles y que sean las columnas:
   INPFILE = '/00/model.WRF_DET_4km.'+FECHA_INI+'_'+CICLO+'0000.000.OPERLEV.nc'  # Salidas del WRF deterministico
   filename_POST = INPDIR + INPFILE

   if not os.path.isfile(filename_POST):
      print(filename_POST)
      print('No se encontro el POST. Salgo')
      return

   f = Dataset(INPDIR+INPFILE)
   NIVELES = np.squeeze(f.variables['lev'][:]).astype(int)      # 19 niveles, el 1ero es 1000, el ultimo es 100 hPa
   Nniveles = len(NIVELES)
 
   ANIO = FECHA_INI_obj.year
   MES = FECHA_INI_obj.month
   DIA = FECHA_INI_obj.day
   
   DateIndex = pd.date_range(datetime.date(ANIO, MES, DIA), periods=1,freq = 'D')
   NPlazo = np.array(np.arange(PLAZOS + 1))
   NVariables = len(VARIABLES) # cant de variables a guardar en el dataframe
 
   # Creamos estructura de indices 
   #miindex = pd.MultiIndex.from_product([DateIndex,[CICLO],NPlazo,ciudades['Numero'],VARIABLES])
   miindex = pd.MultiIndex.from_product([DateIndex, [CICLO], NPlazo, CIPI, VARIABLES])
 
   # Creamos dataFrame vacio 
   DFwrflev = pd.DataFrame(np.full((len(miindex), Nniveles), np.nan), index = miindex, columns=NIVELES)
 
#   DFwrflev['Validez'] = DFwrflev.index[0][0]  # ==> agrega una columna llena de datos iguales al 1er indice
# 
#   for indice,row in DFwrflev.iterrows():
#      indice_fecha = indice[0]
#      indice_ciclo = timedelta(hours=int(indice[1]))
#      indice_plazo = timedelta(hours=int(indice[2]))
#      validez = indice_fecha + indice_ciclo + indice_plazo
#      DFwrflev['Validez'].loc[indice[0],indice[1],indice[2],indice[3]] = validez  # <== asigno el valor de validez a la columna 'Validez'

   DFwrflev['Validez'] = DFwrflev.index[0][0]  # ==> agrega una columna llena de datos iguales al 1er indice
   validez_ini = FECHA_INI_obj
   for ip in NPlazo.tolist():
      DFwrflev.loc[(slice(None), slice(None), ip, slice(None), slice(None)), 'Validez'] = validez_ini + timedelta(hours = ip)

   
   # Agregamos nombres a los indices del DFwrf que no tenian:
   nombres = ['Fecha','Ciclo','Plazo','CIPI','Variable']
   DFwrflev = DFwrflev.rename_axis(nombres)
 
   DFwrflev.reset_index(inplace = True) # reseteo los indices y los pone como columnas de datos
   DFwrflev = DFwrflev.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez' ,'CIPI', 'Variable']) # seteo los indices con este nuevo orden
 
 # ===========================================================================
 # Paso 2: LECTURA DE LAS SALIDAS NETCDF DE ALTURA OPERATIVAS
 # ===========================================================================
 
   # Loop para leer cada plazo del pronostico
   for ip in NPlazo:
       PLAZO = str(ip).zfill(2)
       val = datetime.datetime.strftime(FECHA_INI_obj + timedelta(hours=float(ip)) ,'%Y-%m-%d %H:00:00')
       INPFILE = '/00/model.WRF_DET_4km.'+FECHA_INI+'_'+CICLO+'0000.0'+PLAZO+'.OPERLEV.nc'           # Salidas del WRF deterministico
       filename_POST = INPDIR + INPFILE

       if not os.path.isfile(filename_POST):
          continue

       f = Dataset(INPDIR+INPFILE)
       #f.variables  # Me da el nombre de todas las variables del archivo
       lista = f.variables.keys()  # muestra la lista de variables que tiene el .nc
       CLDFRA = np.squeeze(f.variables['CLDFRA'][:]) # (19, 1249, 999) para 1 plazo fijo
       niveles = f.variables['lev'][:]	# 19 niveles [[1000.,975.,950.,925.,900.,850.,800.,750.,700.,650.,600.,550.,500.,400.,300.,250.,200.,150.,100.]
       f.close()
       for indice,row in Puntos_interes.iterrows(): # loop para todas las estaciones: 
           i = row['i']
           j = row['j']
           numero = row['CIPI']
           #print(Numero,row['Nombre']) 
           vector = CLDFRA[:,i,j].data # array con los valores de cldfra en el punto
           DFwrflev.loc[FECHA_INI, CICLO, int(PLAZO), val, numero, 'CLDFRA'] = vector # para cada fila, ponemos los valores del vector
   
   # reemplazo faltantes fill_value=9.96921e+36 por NaN:
   DFwrflev[DFwrflev >= 1e20] = np.nan

# ===========================================================================
# Paso 3: GUARDADO DEL DATAFRAME   
# ===========================================================================

   #Agrego las columnas de tipo y ID
   DFwrflev = DFwrflev.join(cols, on = 'CIPI')

   #Agrego las columnas de tipo y ID como indices
   DFwrflev.reset_index(inplace = True) # reseteo los indices y los pone como columnas de datos
   DFwrflev = DFwrflev.set_index(['Fecha', 'Ciclo', 'Plazo', 'Validez', 'CIPI', 'Tipo', 'ID', 'Variable']) # seteo los indices con este nuevo orden


   DFwrflev.to_hdf(OUTDIR + '/' + FECHA_INI + '_' + CICLO + '0000' + '_lev.h5', key = 'DFwrflev', mode = 'w')
   


if __name__ == "__main__":

   print("*************************** empieza a correr el DFwrf_LEV ****************************")

   main()
   
   print('Fin >>> Se genero el DF con los datos de cada plazo para cada estacion/ciudad')



