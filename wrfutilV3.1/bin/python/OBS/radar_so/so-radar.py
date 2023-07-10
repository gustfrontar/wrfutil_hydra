"""
Script para correr el superobbing de datos de radar meteorológico.

 La configuración del superobbing se hace a través de un archivo de configuración
que define variables de entorno. Estas variables se levantan con la clase BaseOptions.
Para poder iniciar la clase BaseOptions es necesario proporcionarle una carpeta con los
archivos de radar a procesar, una lista de nombres de radares a procesar, y el período
de tiempo de obtención de los datos.

"""
from utils.base_options import BaseOptions
from utils.files_to_dataframe import get_so_radar_data
from superobbing import radar_superobbing
import os
import pickle
import numpy as np
import pyart
from datetime import datetime, timedelta

if __name__ == '__main__':
    #Levantamos la configuracion del experimento.
    opt = BaseOptions()

    print('Leo los archivos')
    df = get_so_radar_data(root=opt.dataroot, opt=opt)
    print(df['file_path'])
    radar_list = sorted(list(df['file_path']))
    for corrected_radar in radar_list:
        superobbing = radar_superobbing(opt.reference_date,
                                        opt.window_length,
                                        opt.window_type,
                                        corrected_radar,
                                        opt.output_path,
                                        opt)
    os.system("rm -rf {}".format(opt.output_path+"/grid"))
