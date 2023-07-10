import os
import subprocess
import numpy as np
import pandas as pd
import datetime
from datetime import datetime, timedelta


def get_so_radar_data(root=None, opt=None):
    """
    Parametros
    ----------
    root : str
        Nombre de la carpeta que contiene a los archivos de radar.

    opt : clase BaseOptions
        Opciones que vienen de las variables de entorno.

    Devuelve
    --------
    df : DataFrame de pandas indexado por fecha.
        Las columnas del DataFrame son:
        'file_path': Path completo del archivo.
        'ext': Extension del archivo (nc, vol, H5).
        'radar_id': Nombre del radar (RMA1, RMA2,..., ANG, PAR, PER).

    """

    # Genero una lista con todos los archivos y sus path
    lista_radar = [os.path.join(path, name)
                   for path, subdirs, files in os.walk(root)
                   for name in files]

    # Compute slot inital and final time according to window type
    if isinstance(opt.reference_date, str):
       reference_date = datetime.strptime(opt.reference_date, '%Y%m%d_%H%M%S')

    if opt.window_type == 'centered':
        start_date = reference_date - timedelta(minutes=opt.window_length/2.0)
        end_date = reference_date + timedelta(minutes=opt.window_length/2.0)
    elif opt.window_type == 'backward':
        start_date = reference_date - timedelta(minutes=opt.window_length)
        end_date = reference_date
    elif opt.window_type == 'forward':
        start_date = reference_date
        end_date = reference_date + timedelta(minutes=opt.window_length)
    # Genero el DataFrame
    df = pd.DataFrame(lista_radar)
    df['file_path'] = df[0]  # genero columna file_path con todo el path

    # Get the dates
    dates = [datetime.strptime(file.split('.')[-3], '%Y%m%d_%H%M%S')
             for file in df['file_path']]
    df['dates'] = dates
    df = df.set_index('dates') #set dates as index
    mask = (df.index > start_date) & (df.index <= end_date)
    #print("DF antes de la mascara:{}".format(df))
    #df = df.loc[mask]
    df = df[mask]
    #print("DF despues de la mascara:{}".format(df))
    # Get the radar name
    radar_id = [file.split('.')[-2].split('_')[0] for file in df['file_path']]
    df['radar_id'] = radar_id

    instruments = opt.instrument_list
    print("Lista de radares:{}".format(instruments))
    df = df[df['radar_id'].isin(instruments)]

    return df
