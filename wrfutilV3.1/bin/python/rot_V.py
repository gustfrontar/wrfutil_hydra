# Modulo para rotar el viento
# Sept 2016 - Maru Dillon, Maxi Sacco
# Julio 2021 - Modificaciones para realizar en cada corrida el archivo de los indices
# ------------------------------------------------- 

import numpy as np
import numpy.matlib
import os
from datetime import datetime
from datetime import timedelta
import pygrib
import sys

# ********* rotate_uv ************

class rotateUV():

  def __init__(self,filename):
      self.alpha=0;
      self.alphafile=filename
  def calcularAlpha(self,truelat1,stand_lon,lon2d,lat2d):
      RAD_PER_DEG = np.pi/180
      cone = np.sin(np.abs(truelat1)*RAD_PER_DEG)
      diff = lon2d - stand_lon
      for i,val in np.ndenumerate(diff):
            if val > 180 :
                  self.diff[i] = val - 360
            if val < -180 :
                  self.diff[i] = val + 360
            self.alpha = np.zeros(np.shape(lat2d))

            for i,val in np.ndenumerate(lat2d):
                if val < 0 :
                        self.alpha[i] = -diff[i] * cone * RAD_PER_DEG
                else:
                  self.alpha[i] = diff[i] * cone * RAD_PER_DEG
  def rotarUV(self,u,v):
    urot= v * np.sin(self.alpha) + u * np.cos(self.alpha)
    vrot= v * np.cos(self.alpha) - u * np.sin(self.alpha)
    return (urot, vrot)

  def guardarAlpha(self):
    np.save(self.alphafile,self.alpha)    
  def cargarAlpha(self):
    self.alpha=np.load(self.alphafile)


def rotate_uv(u,v,truelat1,stand_lon,lon2d,lat2d):

  """
  Para aplicar a u,v en coordenadas x-y de proyeccion Lambert, y rotarlas  
  a coordenadas de la Tierra W-E,N-S, es decir a las coordenadas posta de U y V
  que son las coordenadas de la proyeccion cilindrica equidistante (o sea lat-lon)
  basado en module_calc_uvmet.f90 del src del ARWpost

  Asumo que truelat1=truelat2 (sino hay que hacer otra cuenta para el cone)
  
  Keywords args:
  u,v -- vientos coord x-y 
  truelat1, stand_lon -- valores del namelist.wps
  lon2d, lat2d -- array de 2 dimensiones de longitudes y latitudes armado con el meshgrid
  
  Output: 
  (urot,vrot) 2 array separados 
  """

  RAD_PER_DEG = np.pi/180
  cone = np.sin(np.abs(truelat1)*RAD_PER_DEG)

  diff = lon2d - stand_lon

  for i,val in np.ndenumerate(diff):
        if val > 180 :
               diff[i] = val - 360
        if val < -180 :
               diff[i] = val + 360

        alpha = np.zeros(np.shape(lat2d))

  for i,val in np.ndenumerate(lat2d):
        if val < 0 :
               alpha[i] = -diff[i] * cone * RAD_PER_DEG
        else:
               alpha[i] = diff[i] * cone * RAD_PER_DEG

        urot = v * np.sin(alpha) + u * np.cos(alpha)
        vrot = v * np.cos(alpha) - u * np.sin(alpha)

  return (urot, vrot) 

def actualizarIndicesUV(gribFH):
  U = grbs.select(name='U component of wind')
  U10 = grbs.select(name='10 metre U wind component')
  U.extend(U10)
  Ui = np.array(list(map(lambda x:x.messagenumber,U)))
  Ui.sort()
  return Ui





###### Desde aca es el script #######

if __name__ == "__main__" :
  #filename="/mini_r5c/WRF/20180605_12/UPPOUT/wrfprs_d01.019"
  if len(sys.argv)==4:
    filenameIN=sys.argv[2]
    filenameOUT=sys.argv[3]
    fileAlpha=str(sys.argv[1])
  else:
    print("ERROR: debe pasar 3 argumentos: <archivo de alphas> <archivo de entrada> <archivo de salida>")
    exit(1)

  if not os.path.isfile(filenameIN):
    print(f"ERROR: el archivo {filenameIN} no existe") 
    exit(1)
  grbs = pygrib.open(filenameIN)  ## Obtenemos un File Handle del archivo grib
  # Para cargar variables usamos:
  Ui=actualizarIndicesUV(grbs)
  Ui=Ui-1
  # Abrimos el archivo grib
  grbout = open(filenameOUT,'wb')
  # Creamos el objeto rotador que contiene a la matriz alpha
  rotador=rotateUV(fileAlpha)

  if not os.path.isfile(fileAlpha):
    # Para generar las matrices lat lon hicimos:  (esto lo necesitamos 1 vez)
    # Las siguientes 6 lineas se deben descomentar si se desea generar una nueva matriz alpha
    grbs.seek(Ui[0])
    truelat1 = sys.environ["TRUELAT1"]
    stand_lon = sys.environ["STAND_LON"]
    lat, lon = grbs.read(1)[0].latlons()
    rotador.calcularAlpha(truelat1,stand_lon,lon,lat)
    rotador.guardarAlpha()
  rotador.cargarAlpha()

  inicio=0
  actual=0
  fin=grbs.messages
  grbs.seek(inicio)
  for i in Ui:
    mensajes=grbs.read(i-actual)
    for grb in mensajes:
      grbout.write(grb.tostring())
    mensajes=grbs.read(2)
    (urot,vrot) = rotador.rotarUV(mensajes[0]['values'],mensajes[1]['values'])
    mensajes[0]['values']=np.copy(urot)
    mensajes[1]['values']=np.copy(vrot)
    grbout.write(mensajes[0].tostring())
    grbout.write(mensajes[1].tostring())
    actual=i+2
  mensajes=grbs.read(fin-actual)
  for grb in mensajes:
    grbout.write(grb.tostring())
  grbout.close()
  grbs.close()
