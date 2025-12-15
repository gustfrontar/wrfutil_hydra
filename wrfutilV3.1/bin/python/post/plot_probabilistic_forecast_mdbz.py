import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature # Import the feature module
from cartopy.feature import NaturalEarthFeature
import os
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import prob_utils as utils

START_INIT_DATE_STR = '20240319150000'
END_INIT_DATE_STR = '20240320060000'
INIT_FREQUENCY_HOURS = 3              
VARIABLE_NAME = 'MDBZ'                # The original variable name (e.g., PRCP, TEMP, etc.)
OUTPUT_VAR_NAME = f'Prob_{VARIABLE_NAME}' # The name of the probabilistic variable in the NetCDF
BASE_DIR = '../../../'
PLOT_DIR = BASE_DIR + f'/PLOT/PROB/{VARIABLE_NAME}/'        # Directory to save the output images
DATA_DIR = BASE_DIR + '/POST/DAFCST/'
RADAR_MOSAIC_DIR='/home/ra000007/a04037/data/DATOS_RADAR/mosaic_data/'



def plot_probability_field(data_array,obs,obs_lon,obs_lat,obs_mask, threshold, init_time, valid_time, plot_dir):
    """Generates and saves a single map plot of the probability field."""
    
    # Extract the map projection information (Assumes WRF files use standard coordinates)
    # You may need to adjust the projection setup based on the actual WRF grid info
    # For a general non-geospatial plot, you can use ax=plt.axes()
    
    fig = plt.figure(figsize=(12, 10))
    # Example: Using PlateCarree projection for simplicity. For WRF, you might need a Lambert Conformal projection.
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree()) 
    
    # Define color levels for probability (0 to 1)
    levels = np.linspace(0, 1, 11) 
    
    # Plotting the data
    im = data_array.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(), # Data coordinate system
        levels=levels,
        cmap='BuGn',
        extend='max',
        x='longitude',
        y='latitude',
        cbar_kwargs={'label': 'Probability (0-1)'}
    )
    ax.contour( obs_lat , obs_lon , obs , transform=ccrs.PlateCarree(), levels=np.linspace(-0.5,1.0,3),colors='b',linewidths=1.0)
    ax.contour( obs_lat , obs_lon , obs_mask , transform=ccrs.PlateCarree(), levels=np.linspace(-0.5,1.0,3),colors='k',linewidths=0.2)

    # Map features
    ax.coastlines(resolution='50m')
    # Add Country Borders
    # Use the 50m resolution for political borders
    ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=0.8, edgecolor='gray', zorder=5)
    # --- NEW: Add State/Province Boundaries ---
    # Use NaturalEarthFeature to load the 10m resolution state/province boundaries
    states_provinces = NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        edgecolor='darkorange',
        facecolor='none'
    )
    ax.add_feature(states_provinces, linestyle='-', linewidth=0.8, zorder=5)
    # Add other features for context (optional)
    ax.add_feature(cfeature.LAND, facecolor='#eeeeee')
    ax.add_feature(cfeature.OCEAN, facecolor='#ddddff')
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    ax.set_title(
        f'{OUTPUT_VAR_NAME} > {threshold} mm\n'
        f'Init: {init_time.strftime("%Y-%m-%d %H:%M UTC")} | Valid: {valid_time.strftime("%Y-%m-%d %H:%M UTC")}',
        fontsize=14
    )
    
    # Save the plot
    filename = f'{OUTPUT_VAR_NAME}_{init_time.strftime("%Y%m%d%H")}_valid_{valid_time.strftime("%Y%m%d%H%M")}_thr_{threshold:.1f}.png'
    
    plt.savefig(os.path.join(plot_dir, filename), bbox_inches='tight')
    plt.close(fig)
    print(f"  --> Saved {filename}")


# Create the plot directory if it doesn't exist
os.makedirs(PLOT_DIR, exist_ok=True)


# Convert start/end strings to datetime objects
current_init_time = datetime.strptime(START_INIT_DATE_STR, '%Y%m%d%H%M%S')
end_init_time = datetime.strptime(END_INIT_DATE_STR, '%Y%m%d%H%M%S')
init_interval = timedelta(hours=INIT_FREQUENCY_HOURS)

# --- MAIN INITIALIZATION LOOP ---
while current_init_time <= end_init_time:
    
    # Format strings for the current initialization
    INIT_DATE = current_init_time.strftime('%Y%m%d%H%M%S')
    INPUT_FILENAME = DATA_DIR + f'/{INIT_DATE}/' + f'{OUTPUT_VAR_NAME}_{INIT_DATE}.nc' 
    
    print(f"\n===================================================")
    print(f"🖼️ Plotting probabilities for Init: {INIT_DATE}")
    print(f"===================================================")

    if not os.path.exists(INPUT_FILENAME):
        print(f"File not found: {INPUT_FILENAME}. Skipping this initialization.")
        current_init_time += init_interval
        continue

    # Load the probabilistic NetCDF file
    with xr.open_dataset(INPUT_FILENAME) as ds:
       ds = utils.preprocess_wrf(ds)
       prob_data = ds[OUTPUT_VAR_NAME] # Dimensions: (time, threshold, lat, lon)
       prob_data = prob_data.where( prob_data > 0.01 ) 
            
       # --- ITERATE THROUGH TIME AND THRESHOLDS ---
            
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
                    
             # Loop over all thresholds
             for threshold_value in prob_data['threshold'].values:
                ob_to_plot =  ( mdbz > threshold_value ).astype(float)
                # Select the data for the current time and threshold
                # .sel() is used to select coordinates
                da_to_plot = prob_data.sel(time=time_step, threshold=threshold_value).squeeze()
                # Call the plotting function
                plot_probability_field(
                  da_to_plot,
                  ob_to_plot,obs_lat,obs_lon,obs_mask,
                  threshold_value, 
                  current_init_time, 
                  valid_time, 
                  PLOT_DIR
                )


    # Advance to the next initialization time
    current_init_time += init_interval

print("\nPlotting process complete.")


