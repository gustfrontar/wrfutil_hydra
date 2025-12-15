import xarray as xr
import matplotlib.pyplot as plt
import os
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import prob_utils as utils
import copy

START_INIT_DATE_STR = '20240319150000'
END_INIT_DATE_STR = '20240320060000'
INIT_FREQUENCY_HOURS = 3              
VARIABLE_NAME = 'MDBZ'                # The original variable name (e.g., PRCP, TEMP, etc.)
OUTPUT_VAR_NAME = f'Prob_{VARIABLE_NAME}' # The name of the probabilistic variable in the NetCDF
BASE_DIR = '../../../'
DATA_DIR = BASE_DIR + '/POST/DAFCST/'
RADAR_MOSAIC_DIR='/home/ra000007/a04037/data/DATOS_RADAR/mosaic_data/'
OUT_FILE=BASE_DIR + '/POST/' + OUTPUT_VAR_NAME + 'prob_scores.npz'



# Convert start/end strings to datetime objects
current_init_time = datetime.strptime(START_INIT_DATE_STR, '%Y%m%d%H%M%S')
end_init_time = datetime.strptime(END_INIT_DATE_STR, '%Y%m%d%H%M%S')
init_interval = timedelta(hours=INIT_FREQUENCY_HOURS)

# --- MAIN INITIALIZATION LOOP ---
main_auc_list = []
main_tpr_list = []
main_fpr_list = []
main_rfp_list = []
main_rof_list = []
main_rfn_list = []

while current_init_time <= end_init_time:
    auc_list = []
    tpr_list = []
    fpr_list = []
    rfp_list = []
    rof_list = []
    rfn_list = []
    
    # Format strings for the current initialization
    INIT_DATE = current_init_time.strftime('%Y%m%d%H%M%S')
    INPUT_FILENAME = DATA_DIR + f'/{INIT_DATE}/' + f'{OUTPUT_VAR_NAME}_{INIT_DATE}.nc' 
    
    print(f"\n===================================================")
    print(f"🖼️ Computing scores for Init: {INIT_DATE}")
    print(f"===================================================")

    if not os.path.exists(INPUT_FILENAME):
        print(f"File not found: {INPUT_FILENAME}. Skipping this initialization.")
        current_init_time += init_interval
        continue

    # Load the probabilistic NetCDF file
    with xr.open_dataset(INPUT_FILENAME) as ds:
       ds = utils.preprocess_wrf(ds)
       prob_data = ds[OUTPUT_VAR_NAME] # Dimensions: (time, threshold, lat, lon)
       #prob_data = prob_data.where( prob_data > 0.01 ) 
            
       # --- ITERATE THROUGH TIME AND THRESHOLDS ---
       thresholds = prob_data['threshold'].values        
       # Loop over all time steps
       for time_step in prob_data['time'].values:
          valid_time = pd.to_datetime(time_step) # Convert numpy datetime64 to pandas/python datetime
          obs_file = RADAR_MOSAIC_DIR + '/mdbz_mosaic.' + valid_time.strftime('%Y%m%d%H%M') + '00.npz'
          # Load the observed data. 
          if os.path.exists( obs_file ) :
             print('Will read the following file : ' + obs_file )
             tmp_data = np.load( obs_file )  
             mdbz = tmp_data['MDBZ']
             mdbz[ mdbz == 10.0 ] = np.nan
             obs_mask = np.isnan( mdbz ) 
             obs_lat = tmp_data['lat']
             obs_lon = tmp_data['lon']
             prob_forecast = prob_data.sel(time=time_step).squeeze().to_numpy()
             print( prob_forecast.shape )
             forecast_lon = prob_data.longitude.values
             forecast_lat = prob_data.latitude.values
                    
             # Compute ROC scores
             tmp_auc , tmp_fpr , tmp_tpr , tmp_thr = utils.compute_multiple_roc_metrics_nan_aware(
                 mdbz,
                 obs_lat,
                 obs_lon,
                 prob_forecast,
                 forecast_lon,
                 forecast_lat,
                 thresholds,
                 interpolation_method = 'linear' )

             # Compute Reliability scores
             tmp_rfp , tmp_rof , tmp_rfn , tmp_nbins , tmp_thr = utils.compute_reliability_curve(
                 mdbz,
                 obs_lat,
                 obs_lon,
                 prob_forecast,
                 forecast_lon,
                 forecast_lat,
                 thresholds,
                 n_bins = 10, # Number of probability bins (e.g., 10 for 0.1 increments)
                 interpolation_method = 'linear' )
       auc_list.append( copy.deepcopy( auc_list ) )
       tpr_list.append( copy.deepcopy( tpr_list ) )
       fpr_list.append( copy.deepcopy( fpr_list ) )
       rfp_list.append( copy.deepcopy( rfp_list ) )
       rof_list.append( copy.deepcopy( rof_list ) )
       rfn_list.append( copy.deepcopy( rfn_list ) )
     

    main_auc_list.append( copy.deepcopy( auc_list ) )
    main_tpr_list.append( copy.deepcopy( tpr_list ) )
    main_fpr_list.append( copy.deepcopy( fpr_list ) )
    main_rfp_list.append( copy.deepcopy( rfp_list ) )
    main_rof_list.append( copy.deepcopy( rof_list ) )
    main_rfn_list.append( copy.deepcopy( rfn_list ) )

    # Advance to the next initialization time
    current_init_time += init_interval

print("\nPlotting process complete.")


auc = np.array( main_auc_list )
tpr = np.array( main_tpr_list )
fpr = np.array( main_fpr_list )
rfp = np.array( main_rfp_list )
rof = np.array( main_rof_list )
rfn = np.array( main_rfn_list )

np.savez(OUT_FILE, auc=auc , tpr=tpr , fpr=fpr , rfp=rfp , rof=rof , rfn=rfn , nbins=tmp_nbins , thresholds=thresholds , times = prob_data['time'].values )

