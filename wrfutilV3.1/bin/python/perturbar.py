# -*- coding: utf-8 -*-
from netCDF4 import Dataset
import numpy as np
import glob
import os
import sys
import datetime
from multiprocessing import Pool

##########################
# version paralela
# Crear las carpetas destino por fuera de este script
# Funciona hacer links de los met_em desde DET1/WPS a /V2/originales
#########################

## Calcula la media del ensamble.
## Esto ya viene calculado, verificar si usarla lista es mas rapido
def media(ensamble,varlist):
    var_mean={}
    for iv in varlist:
        var=np.array([miem.variables[iv][:] for miem in ensamble.values()])
        var_mean[iv] = (np.mean(var,axis=0),len(var[0].shape))
    return var_mean

# Esto calcula para perturvacion respecto de la media para cada variale
def perturbacion(ensamble,varlist,media):
    out={}
    for iv in varlist:
        out[iv]=(np.array([miem.variables[iv][:]-media[iv][0] for miem in ensamble.values()]),media[iv][1])
    return out

# Dada un perturbacion, perturba un original
# y escibe el archivo de salida
def perturbar(ensamble,pert,varlist,met_det,fileout):
    for i,miem in enumerate(ensamble.values(),1):
         fileout[str(i)].setncatts(met_det.__dict__)
         for dname, the_dim in met_det.dimensions.items():
             fileout[str(i)].createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
         for v_name, varin in met_det.variables.items():
             outVar = fileout[str(i)].createVariable(v_name, varin.datatype, varin.dimensions)
             outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
             outVar[:] = varin[:]
         for iv in varlist:	
             if pert[iv][1] == 4:
                 fileout[str(i)].variables[iv][0,:,:,:] = pert[iv][0][i-1,:,:,:][np.newaxis,:,:,:] + met_det.variables[iv][:]
             elif pert[iv][1] == 3:
                 fileout[str(i)].variables[iv][0,:,:] = pert[iv][0][i-1,:,:][np.newaxis,:,:] + met_det.variables[iv][:]

# Esto es la funcion que se encarga de separar en Threads, tantos como pedacitos haya
def pert_thread(varlist,fecha,nmiembros,p,nped):
    if nped >1 :
       deterministico = Dataset(path_in + '00/met_em.d01.'+fecha.strftime('%Y-%m-%d_%H:%M:%S')+'.nc_'+"{:04d}".format(p)) # DET
    else:
       deterministico = Dataset(path_in + '00/met_em.d01.'+fecha.strftime('%Y-%m-%d_%H:%M:%S')+'.nc') 
    ensamble={}
    fileout={}
    for i in range(1,nmiembros+1):
        if nped>1 :
             ensamble[str(i)] = Dataset(glob.glob(path_in  + "{:02d}".format(i)+ '/met_em.d01.'+fecha.strftime('%Y-%m-%d_%H:%M:%S')+'.nc_'+"{:04d}".format(p))[0])
             fileout[str(i)] = Dataset(path_out + "{:02d}".format(i) + '/met_em.d01.'+fecha.strftime('%Y-%m-%d_%H:%M:%S')+'.nc_'+"{:04d}".format(p),'w') # 00_pert/ deberia estar creada
        else :
             ensamble[str(i)] = Dataset(glob.glob(path_in  + "{:02d}".format(i)+ '/met_em.d01.'+fecha.strftime('%Y-%m-%d_%H:%M:%S')+'.nc')[0])
             fileout[str(i)] = Dataset(path_out + "{:02d}".format(i) + '/met_em.d01.'+fecha.strftime('%Y-%m-%d_%H:%M:%S')+'.nc','w') # 00_pert/ deberia estar creada
    var_mean=media(ensamble,varlist)
    pert=perturbacion(ensamble,varlist,var_mean)
    perturbar(ensamble,pert,varlist,deterministico,fileout)

    for i in range(1,nmiembros+1):
        fileout[str(i)].close()
        ensamble[str(i)].close()
    deterministico.close()




# El resto de los parametros estan en Variables de Entorno de BASH
# que se configuran via el experimento.conf
fecha = datetime.datetime.strptime(os.environ['FECHA_PERT'], '%Y%m%d_%H%M%S')  # argumento de entrada => fecha de la corrida 20190424_000000
npedacitos=int(os.environ['WPSPROC'])
nmiembros=int(os.environ['MIEMBROS'])
nhoras = int(os.environ['PLAZO'])+1
path_in = os.environ['MET_ORI'] # Path donde se encuentran los Mets originales 
path_out = os.environ['MET_PERT'] # Path donde se guardaran los Mets perturbados
ARRAY_CNT=int(os.environ['ARRAYCNT']) # Cantidad de procesos pedidos
ARRAY_ID=int(os.environ['ARRAYID']) # El proceso que soy yo!
THREADS_NUM=int(os.environ['OMP_NUM_THREADS'])

# Aca tratamos de balancear lo mejor posible las horas de pronostico entre los nodos disponibles.
horas_por_proceso=int(nhoras/ARRAY_CNT)
resto=nhoras - (horas_por_proceso * ARRAY_CNT)
hora_ini=horas_por_proceso*ARRAY_ID
hora_fin=hora_ini+horas_por_proceso
if (ARRAY_ID < resto):
         hora_ini=hora_ini+ARRAY_ID
         hora_fin=hora_fin+ARRAY_ID+1
else:
         hora_ini=hora_ini+resto
         hora_fin=min(hora_fin+resto,nhoras)
	

print("SOY: "+str(ARRAY_ID)+" ; de un total de : "+str(ARRAY_CNT)+" y voy a correr de: "+str(hora_ini) +" hasta: "+str(hora_fin))

# las cuentas las hacemos sobre las variables dinamicas   

varlist = ["PRES", "SM", "ST", "GHT", "SKINTEMP",
           "ST100200", "ST040100", "ST010040", "ST000010",
           "SM100200", "SM040100", "SM010040", "SM000010",
           "PSFC", "RH", "UU", "VV", "TT", "PMSL" ]
pool = Pool(processes=THREADS_NUM)
pool.starmap(pert_thread,[(varlist,fechaHora,nmiembros,p,npedacitos) for p in range(npedacitos) for fechaHora in [fecha + datetime.timedelta(hours=1*i) for i in range(hora_ini,hora_fin)]])








