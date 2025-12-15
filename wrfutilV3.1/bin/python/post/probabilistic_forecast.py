import xarray as xr
import numpy as np
import glob
import os
from datetime import datetime, timedelta
import prob_utils as utils
import argparse

START_INIT_DATE_STR = '20240319150000' # The first initialization time (e.g., 12 UTC March 19)
END_INIT_DATE_STR = '20240320060000'   # The last initialization time to process
INIT_FREQUENCY_HOURS = 3               # Frequency of initialization (every 3 hours)

# --- General Configuration (mostly the same) ---
BASE_DIRECTORY = '../../../POST/DAFCST/' # Parent directory containing all initialization folders
ENSEMBLE_SIZE = 60
FORECAST_LENGTH_HOURS = 24
TIME_STEP_MINUTES = 30


# --- Argument Parsing Setup ---
def parse_arguments():
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute probabilistic forecasts from WRF ensemble output files."
    )
    # Argument for the variable name (required)
    parser.add_argument(
        'variable_name', 
        type=str, 
        help='The name of the variable to process (e.g., PRCP, T2, QRAIN).'
    )
    # Argument for the thresholds list (required)
    parser.add_argument(
        'thresholds',
        type=str,
        help='A comma-separated list of thresholds (e.g., "0.1,1.0,5.0").'
    )
    return parser.parse_args()

# Execute argument parsing
args = parse_arguments()

# Set the VARIABLE_NAME
print( args )
VARIABLE_NAME = args.variable_name
OUTPUT_VARIABLE_NAME = f'Prob_{VARIABLE_NAME}'

# Set the THRESHOLDS list
try:
    # Split the string by comma, strip whitespace, and convert each part to float
    THRESHOLDS = [float(t.strip()) for t in args.thresholds.split(',')]
    if not THRESHOLDS:
        raise ValueError("Thresholds list is empty.")
    THRESHOLDS.sort() # Sorting is good practice, though not strictly required
except ValueError as e:
    print(f"Error parsing thresholds: {e}. Please ensure they are comma-separated numbers.")
    sys.exit(1) # Exit if thresholds are invalid
    
print(f"✅ Processing variable: **{VARIABLE_NAME}**")
print(f"✅ Using thresholds: **{THRESHOLDS}**")

TIME_STEPS = int(FORECAST_LENGTH_HOURS * 60 / TIME_STEP_MINUTES)

# Convert start/end strings to datetime objects
current_init_time = datetime.strptime(START_INIT_DATE_STR, '%Y%m%d%H%M%S')
end_init_time = datetime.strptime(END_INIT_DATE_STR, '%Y%m%d%H%M%S')
init_interval = timedelta(hours=INIT_FREQUENCY_HOURS)

# --- MAIN INITIALIZATION LOOP ---
while current_init_time <= end_init_time:
    
    # Format strings for the current initialization
    INIT_DATE = current_init_time.strftime('%Y%m%d%H%M%S')
    INIT_FOLDER = os.path.join(BASE_DIRECTORY, INIT_DATE)
    OUTPUT_FILENAME = INIT_FOLDER + '/' + OUTPUT_VARIABLE_NAME + f'_{INIT_DATE}.nc' 
    
    print(f"\n=======================================================")
    print(f"🚀 Starting Processing for Initialization: {INIT_DATE}")
    print(f"=======================================================")

    # List to store the final probability DataArrays for all time steps of this initialization
    all_time_steps_prob = []

    # --- FORECAST TIME STEP LOOP (Nested) ---
    for i in range(0, TIME_STEPS + 1):
        valid_time = current_init_time + timedelta(minutes=TIME_STEP_MINUTES * i)
        valid_time_str = valid_time.strftime('%Y%m%d%H%M%S')
        
        # Data loading and file globbing logic (same as previous response)
        # -----------------------------------------------------------------
        
        file_pattern = os.path.join(INIT_FOLDER, f'WRF.D{valid_time_str}.T2D.M*.nc')
        file_list = glob.glob(file_pattern)
        
        if len(file_list) != ENSEMBLE_SIZE:
            print(f"Warning: Found {len(file_list)} files for valid time {valid_time_str}, expected {ENSEMBLE_SIZE}. Skipping.")
            continue

        try:
            ensemble_data = []
            for f in file_list:
                with xr.open_dataset(f) as ds: # Use 'with' for safer file handling
                    prcp_var = ds[VARIABLE_NAME].squeeze() 
                    ensemble_data.append(prcp_var)

            prcp_array = xr.concat(ensemble_data, dim='member')

            # --- Probabilistic Calculation for Multiple Thresholds ---
            prob_fields_for_this_time = []
            
            for threshold in THRESHOLDS:
                exceedance_count = (prcp_array > threshold).sum(dim='member')
                probability_exceedance = exceedance_count / ENSEMBLE_SIZE
                
                # Add metadata (threshold, name, attributes)
                probability_exceedance = probability_exceedance.expand_dims(
                    dim={'threshold': [threshold]}
                ).copy()
                probability_exceedance.name = OUTPUT_VARIABLE_NAME
                probability_exceedance.attrs['long_name'] = f'Probability of {VARIABLE_NAME}'
                
                prob_fields_for_this_time.append(probability_exceedance)

            # Combine thresholds
            time_step_prob = xr.concat(prob_fields_for_this_time, dim='threshold')
            
            # Add the valid time coordinate
            time_step_prob = time_step_prob.expand_dims(dim={'time': [valid_time]})
            
            all_time_steps_prob.append(time_step_prob)
            
        except Exception as e:
            print(f"An error occurred while processing valid time {valid_time_str}: {e}")
            
        # -----------------------------------------------------------------

    # --- FINAL CONCATENATION AND NETCDF SAVE for the current INIT_DATE ---
    if all_time_steps_prob:
        final_prob_da = xr.concat(all_time_steps_prob, dim='time')
        final_prob_ds = final_prob_da.to_dataset(name=OUTPUT_VARIABLE_NAME)
        
        final_prob_ds.attrs['title'] = f'Probabilistic ' + VARIABLE_NAME + ' Forecast from WRF Ensemble initialized at ' + INIT_DATE
        #final_prob_ds.Prob_PRCP.attrs['units'] = '1'
        final_prob_ds[OUTPUT_VARIABLE_NAME].attrs['units'] = 'unitless' 
        final_prob_ds.to_netcdf(OUTPUT_FILENAME)
        
        print(f"\n✅ Successfully saved probabilistic forecasts for {INIT_DATE} to: **{OUTPUT_FILENAME}**")
    else:
        print(f"\n⚠️ No data was processed for {INIT_DATE}. Output file was not created.")

    # --- Advance to the next initialization time ---
    current_init_time += init_interval

print("\nProcessing complete for all initialization dates.")

