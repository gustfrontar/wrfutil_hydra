import pyart
import numpy as np
import os
from utils.so_methods import so  
from datetime import datetime, timedelta
from .superobbing import SoOutput
import pickle
from collections import defaultdict

def get_var_names(radar,options):
    """ Busca los nombres de las variables disponibles en el radar """
    keys = [x for x in radar.fields.keys()]
    print("Las variables dentro del radar son : {}".format(keys))
    variables = options.so_vars[0]
    return list(set(keys).intersection(variables))

def get_dates(radar):
    """ Busca la fecha inicial y final del radar """
    initime = radar.time['units'].split(' ')[2]
    initial_date = datetime.strptime(initime, '%Y-%m-%dT%H:%M:%SZ')
    delta = timedelta(seconds=radar.time['data'][-1].tolist())
    end_date = initial_date + delta
    return initial_date, end_date

def datespan(startDate, endDate, delta):
    '''
    La funcion devuelve un "generator" que contiene un objecto date
    Input:
        starDate (objeto): de la clase datetime que indica la fecha inicial
        endDate (objeto): de la clase datetime que indica la fecha final
        delta (objeto): de la clase datetime que indica el intervalo temporal
    '''
    currentDate = startDate
    while currentDate < endDate:
        yield currentDate
        currentDate += delta

def get_letkf_outputs(inidate, enddate, options):
    """
    Create a list of letkf output filenames
    initime, endtime : datetime objects
    out_freq : int in seconds
    """
    freq = options.window_length*60
    ref_date =  datetime(inidate.year, inidate.month, 1, 0, 0, 0)

    # Get first possible output date
    freqtimes = freq*round((inidate-ref_date).total_seconds()/freq)
    outputini = ref_date + timedelta(seconds=freqtimes)

    # Get last possible output date
    freqtimes = freq*round((enddate-ref_date).total_seconds()/freq)
    outputend = ref_date + timedelta(seconds=freqtimes)

    # Create output list
    delta = timedelta(seconds=freq)
    output_dates = []
    #output_dates = defaultdict(list)
    for date in datespan(outputini, outputend+delta, delta):
    #    iniinterval = date - timedelta(seconds=freq/2)
    #    endinterval = date + timedelta(seconds=freq/2)
    #    output_dates[date] = [iniinterval, endinterval]
         output_dates.append(date)

    # Check radar intial and final dates are in output list
    #check_date_in_interval(inidate, output_dates[outputini][0], output_dates[outputini][1])
    #check_date_in_interval(enddate, output_dates[outputend][0], output_dates[outputend][1])

    return output_dates

def check_date_in_interval(date, lower, upper):
    if not lower <= date <= upper:
        raise ValueError('Date not in interval: ' + date)
    return True

def radar_superobbing(reference_date, time_delta, window_type, corrected_radar, outputpath, options):
    """
    Hace el superobbing a un volumen de radar
    
    Parametros
    ----------
    reference_date : str
        fecha y hora de referencia para el calculo del superobbing (slot time)
        
    time_delta : str
        ventana de tiempo donde guardo datos para asimilar en minutos
       
    window_type : str
        Tipo de vetana a emplear: 'centered', 'backward', 'forward' 

    corrected_radar: str
        Path del archivo netcdf del qc-radar

    outputpath : str
        carpeta donde voy a guardar los archivos que luego se van a asimilar
        
    options : class
        objeto SuperobbingOptions con todas las configuraciones del superobbing
        
    Salida
    ------
        un archivo con datos preparados para hacer el analisis
    """
    # Asignamos variables
    if isinstance(reference_date, str):
        reference_date = datetime.strptime(reference_date, "%Y%m%d_%H%M%S")

    # Creamos los directorios donde guardar la salida
    if outputpath == None  :
        outputpath='./'
    os.makedirs(outputpath, exist_ok=True)
    os.makedirs(outputpath + '/grid/' ,exist_ok=True)
    
    # Establecemos el objeto radar
    if isinstance(corrected_radar, str):
        radar = pyart.io.read(corrected_radar)
    else:
        radar = corrected_radar

    print("Entrando Superobbing de : {}".format(radar.metadata['instrument_name']))

    # Get variable names and radar initial and end dates
    vars_names =  get_var_names(radar, options)
    inidate, enddate = get_dates(radar)

    # Create output file list according to time_delta
    output_dates = get_letkf_outputs(inidate, enddate, options)
 
    outfile_list = [] 
    if not vars_names:
        return outfile_list   

    print(" ")
    print("El superobbing corresponde a la fecha {}, usando una ventana {}, con una longitud de {} minutos".format(reference_date, window_type, timedelta(minutes=time_delta)))

    print("El radar empieza a las {} y termina a las {}".format(inidate, enddate))

    for cdate in output_dates:

       # Get radar rays that contribute to current date
       top_second = ((cdate + timedelta(seconds=time_delta*60/2.0)) - inidate).total_seconds()
       bot_second = ((cdate - timedelta(seconds=time_delta*60/2.0)) - inidate).total_seconds()
       rays = np.squeeze(np.where(np.logical_and(radar.time['data'] >= bot_second, radar.time['data'] < top_second)))

       if rays.size > 1:
          ray_limits = [rays[0], rays[-1]]
          print("Para la fecha {} me sirven los rayos entre {}".format(cdate, ray_limits))

          # Hacer el superobbing
          so = SoOutput(radar, options.so_grid, vars_names, ray_limits, options.var_opt, weight_vars = False)

          # Check if there is an exisiting file to update the box average
          tmpfile = "{}/grid/{}_{}.pkl".format(outputpath, radar.metadata['instrument_name'], datetime.strftime(cdate, "%Y%m%d_%H%M%S"))

          if os.path.exists(tmpfile):
             print('Updating boxmean from previous file ' + tmpfile)
             tmp_so = load_object(tmpfile)
             update_box_average(tmp_so, so)

          # Escribir un archivo intermedio
          write_object(tmpfile, so.fields)

          # Escribir un archivo para el LETKF
          #if radar.metadata['instrument_name'] in ['PAR', 'PER', 'ANG']:
          #   outfile = "{}/letkf/{}_{}_{}.dat".format(outputpath, radar.metadata['instrument_name'], int(radar.range['data'].max()), cdate)
          #else:
          #   outfile = "{}/letkf/{}_{}.dat".format(outputpath, radar.metadata['instrument_name'], datetime.strftime(cdate, "%Y%m%d_%H%M%S"))
          outfile = "{}/{}_{}.dat".format(outputpath, radar.metadata['instrument_name'], datetime.strftime(cdate, "%Y%m%d_%H%M%S"))

          outfile_list.append(outfile)
          write_letkf(outfile, so, options)

    return outfile_list    


def load_object(filename):
    """ Abre un archivo pickle """
    with open(filename, 'rb') as f:
        return pickle.load(f)

def update_box_average(tmp, new):
    """ Actualizar el objeto superobing """
    for key in new.fields.keys():
        if key in tmp.keys():
            nobs_old = tmp[key]['nobs']
            nobs_new = new.fields[key]['nobs']
            nobs_tot = nobs_old + nobs_new

            for subkey in new.fields[key].keys():
                if subkey != 'nobs' and subkey != 'id' and subkey != 'error' and subkey != 'min':
                    tmp_sum = nobs_new*new.fields[key][subkey] + nobs_old*tmp[key][subkey] 
                    new.fields[key][subkey][nobs_tot > 0] = tmp_sum[nobs_tot > 0]/nobs_tot[nobs_tot > 0]
            new.fields[key]['nobs'] = nobs_tot

def write_object(filename, obj):
    #print('WRITING PICKLE FILE ' + filename)
    with open(filename, 'wb') as fileout:  # Overwrites any existing file.
        pickle.dump(obj, fileout, pickle.HIGHEST_PROTOCOL)


def write_letkf(output_path, data, options):
    """ Escribe un archivo de letkf"""
    #print("WRITING BINARY LETKF FILE {}".format(output_path))
    
    tmp1 = np.array([4])
    tmp2 = tmp1*7
    nobs = 0
    wk = np.empty(7)

    nvar = len(data.fields.keys())

    tmp_error = np.zeros(nvar)
    tmp_id = np.zeros(nvar)
    tmp_lambda =  12.0 #TODO check this value and get it from the radar object.
                       #Temporarily this value is assumed to be equal to S-BAND 
                       #an specific C-BAND operator has not been coded yet in the LETKF.


    for iv, var in enumerate(data.fields.keys()):
        if iv == 0:
           [nlev, nlat, nlon] = np.shape(data.fields[var]['data'])
           tmp_data = np.zeros((nlev,nlat,nlon,nvar))
           tmp_az = np.zeros((nlev,nlat,nlon,nvar))
           tmp_ra = np.zeros((nlev,nlat,nlon,nvar))
           tmp_el = np.zeros((nlev,nlat,nlon,nvar))
           tmp_n = np.zeros((nlev,nlat,nlon,nvar)).astype(int)

        tmp_data[:,:,:,iv] = data.fields[var]['data']
        tmp_az[:,:,:,iv] = data.fields[var]['az']
        tmp_el[:,:,:,iv] = data.fields[var]['el']
        tmp_ra[:,:,:,iv] = data.fields[var]['ra']
        tmp_n[:,:,:,iv] = data.fields[var]['nobs']
        tmp_error[iv] = data.fields[var]['error']
        tmp_id[iv] = data.fields[var]['id']

    #Filter grid points in which the number of data points is less than min_n observations
    min_n = options.min_nobs
    tmp_n[tmp_n < min_n]=0

    so.write_radar(nlon=nlon,nlat=nlat,nlev=nlev,nvar=nvar,
                   data_in=tmp_data,ndata_in=tmp_n,
                   grid_az=tmp_az,grid_el=tmp_el,grid_ra=tmp_ra,
                   error=tmp_error,ido=tmp_id,lambdar=tmp_lambda,  
                   filename=output_path,
                   radar_lon=data.radar.longitude['data'],
                   radar_lat=data.radar.latitude['data'] ,
                   radar_z=data.radar.altitude['data'] )

    #Write summary.
    for iv , var in enumerate(data.fields.keys()) :
        totalnobs=np.sum( tmp_n[:,:,:,iv] > 0 )
        if totalnobs > 0: # and False:
           print('Total observations for obstype ' + str(np.array(data.fields[var]['id'])) + ' = ' + str(totalnobs) )
           '''
           print('Max value is ' + str(np.max(tmp_data[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ])))
           print('Min value is ' + str(np.min(tmp_data[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ])))
           print('Grid properties')
           print('Max azimuth is  ',np.max(tmp_az[:,:,:,iv][ tmp_n[:,:,:,iv] > 0])) 
           print('Min azimuth is  ',np.min(tmp_az[:,:,:,iv][ tmp_n[:,:,:,iv] > 0]))
           print('Max elevation is  ',np.max(tmp_el[:,:,:,iv][ tmp_n[:,:,:,iv] > 0]))
           print('Min elevation is  ',np.min(tmp_el[:,:,:,iv][ tmp_n[:,:,:,iv] > 0]))
           print('Max range is  ',np.max(tmp_ra[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ]))
           print('Min range is  ',np.min(tmp_ra[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ]))
           '''
