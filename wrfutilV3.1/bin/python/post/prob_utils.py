import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score, confusion_matrix
from scipy.interpolate import griddata
from typing import List, Dict, Any, Tuple
import sys
import xarray as xr
import copy

def create_xarray_dataset(data_array: np.ndarray, lat_array: np.ndarray, lon_array: np.ndarray, var_name: str = 'my_variable') -> xr.Dataset:
    """
    Creates an xarray Dataset from 2D lat/lon arrays and a 2D data array (curvilinear grid).

    Args:
        data_array: 2D NumPy array of the data (shape: (N_y, N_x)).
        lat_array: 2D NumPy array of latitudes (shape: (N_y, N_x)).
        lon_array: 2D NumPy array of longitudes (shape: (N_y, N_x)).
        var_name: The desired name for the data variable in the Dataset.

    Returns:
        A new xarray.Dataset object.
    """
    
    N_y, N_x = data_array.shape
    
    # 1. Verify all shapes are consistent
    if not (data_array.shape == lat_array.shape and data_array.shape == lon_array.shape):
        raise ValueError(
            f"Shape mismatch: All input arrays must have the same 2D shape. "
            f"Data: {data_array.shape}, Lat: {lat_array.shape}, Lon: {lon_array.shape}"
        )

    # 2. Define Dimensions
    # The dimensions are simple integer indices (y and x) since lat/lon are not 1D.
    dims = ('y', 'x')

    # 3. Create the DataArray for the primary variable
    data_var = xr.DataArray(
        data=data_array,
        dims=dims,
        name=var_name,
        # Create default integer coordinate arrays for the dimensions
        coords={'y': np.arange(N_y), 'x': np.arange(N_x)},
        attrs={'units': 'arbitrary_units', 'long_name': 'Generated Sample Data'}
    )

    # 4. Create DataArrays for the 2D coordinates
    lon_coords = xr.DataArray(lon_array, dims=dims, name='longitude')
    lat_coords = xr.DataArray(lat_array, dims=dims, name='latitude')

    # 5. Combine into a Dataset
    ds = data_var.to_dataset()
    
    # 6. Attach the 2D coordinates to the Dataset
    ds['longitude'] = lon_coords
    ds['latitude'] = lat_coords
    
    # 7. Explicitly tell xarray which variables are the non-dimension coordinates
    ds = ds.set_coords(['longitude', 'latitude'])
    
    return ds


def preprocess_wrf(ds):
    """
    Renames the WRF coordinate variables (XLONG, XLAT) to 'longitude' and 
    'latitude' for better xarray/plotting integration.
    """
    if 'XLONG' in ds and 'XLAT' in ds:
        # Rename the coordinate variables themselves
        ds = ds.rename({'XLONG': 'longitude', 'XLAT': 'latitude'})
        
        # Ensure the dimensions are also correctly set if needed (optional, 
        # but good practice if dimensions are generic like 'x'/'y' or 'south_north'/'west_east')
        # You may need to inspect your file to confirm the dimension names. 
        # Assuming the data variable (e.g., 'PRCP') uses dimensions like 
        # (Time, south_north, west_east)
        # If your dimensions are named 'south_north' and 'west_east', 
        # they might need to be linked explicitly.
        
        # If the coordinate variables XLAT/XLONG are defined by dimensions 
        # 'south_north' and 'west_east', we set them as coordinates:
        ds = ds.set_coords(['longitude', 'latitude'])
            
    return ds



def compute_multiple_roc_metrics_nan_aware(
    observed_prcp: np.ndarray, 
    observed_lat: np.ndarray, 
    observed_lon: np.ndarray, 
    prob_forecasts: np.ndarray, 
    forecast_lon: np.ndarray,  
    forecast_lat: np.ndarray,  
    thresholds: List[float],
    interpolation_method: str = 'linear',
    custom_thresholds = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
) -> Dict[str, Any]:
    """
    Computes independent ROC metrics (AUC, FPR, TPR) for each precipitation threshold, 
    ensuring robustness against NaNs in both observed and interpolated forecast fields.

    Args:
        observed_prcp: The 2D array of observed precipitation.
        ... (coordinate and forecast arrays as before) ...
        thresholds: List of precipitation thresholds.
        interpolation_method: Scipy interpolation method.

    Returns:
        A dictionary containing lists/arrays of metrics, indexed by the precipitation threshold.
    """
    
    auc_values = []
    fpr_list = []
    tpr_list = []
    
    # Pre-calculate the observed truth NaNs (needed for masking later)
    observed_nan_mask = np.isnan(observed_prcp.ravel())
    
    # --- Define Interpolation Grid ---
    source_points = np.stack([forecast_lon.ravel(), forecast_lat.ravel()], axis=1)
    target_points = (observed_lon, observed_lat)
    
    # --- Iteration over all Precipitation Thresholds ---
    for idx, precip_threshold in enumerate(thresholds):
        
        print(f"Calculating ROC metrics for event P_obs > {precip_threshold}...")

        # 1. Select the probability field
        original_prob_field = prob_forecasts[idx]

        # 2. Interpolate the probability field onto the observed grid
        source_values = original_prob_field.ravel()
        
        try:
            interpolated_prob_field = griddata(
                source_points, 
                source_values, 
                target_points, 
                method=interpolation_method
            )
        except Exception as e:
            print(f"Error during interpolation for threshold {precip_threshold}: {e}. Skipping.")
            continue

        # 3. Define the Observed Event (The "Truth")
        # Observed event occurs if precipitation exceeds the current threshold
        observed_event = (observed_prcp > precip_threshold).ravel()
        
        # 4. Create the Combined NaN-Aware Mask
        
        # Mask 1: Points where interpolation failed (outside forecast domain)
        forecast_nan_mask = np.isnan(interpolated_prob_field.ravel())
        
        # Combined Mask: True only where BOTH the observed truth AND the interpolated 
        # forecast are valid (i.e., NOT NaN).
        valid_mask = ~(forecast_nan_mask | observed_nan_mask)
        
        if np.sum(valid_mask) <= 1:
            print(f"Warning: Insufficient valid data points ({np.sum(valid_mask)} remaining) for threshold {precip_threshold}. Skipping.")
            continue
            
        # 5. Apply the Mask
        final_observed_event = observed_event[valid_mask]
        final_forecast_probabilities = interpolated_prob_field.ravel()[valid_mask]
        
        # --- 6. Compute ROC Metrics ---
        
        # Compute FPR, TPR (The core ROC curve points)
        # y_true must contain at least one positive and one negative sample for AUC/ROC to work well.
        if len(np.unique(final_observed_event)) < 2:
            print(f"Warning: Observed event ({precip_threshold}) is singular (all event or all non-event) over valid domain. Skipping ROC for this threshold.")
            auc = np.nan
        else :
            # Compute AUC
            auc = roc_auc_score(
               y_true=final_observed_event,
               y_score=final_forecast_probabilities
            )
 
        # Manually defining thresholds
        fpr_custom = []
        tpr_custom = []

        for threshold in custom_thresholds:
           y_score = final_forecast_probabilities
           y_true  = final_observed_event
           y_pred_binary = (y_score >= threshold).astype(int)
           tn, fp, fn, tp = confusion_matrix(y_true, y_pred_binary,labels=(0,1)).ravel()
           fpr_custom.append(fp / (fp + tn))
           tpr_custom.append(tp / (tp + fn))

        auc_values.append(copy.deepcopy(auc))
        fpr_list.append(copy.deepcopy(fpr_custom))
        tpr_list.append(copy.deepcopy(tpr_custom))
        

    # --- 7. Return Results ---
    return np.array( auc_values) , np.array( fpr_list ) , np.array( tpr_list ) , np.array(thresholds) 
     #   "precipitation_thresholds": np.array(thresholds),
     #   "AUC_scores": np.array(auc_values),
     #   "FPR_arrays": fpr_list,
     #   "TPR_arrays": tpr_list,
     #   "description": "NaN-aware ROC metrics computed independently for each precipitation threshold."
     #}
    
    

def compute_dichotomous_scores(
    observed_prcp: np.ndarray, 
    observed_lat: np.ndarray, 
    observed_lon: np.ndarray, 
    prob_forecasts: np.ndarray, 
    forecast_lon: np.ndarray,  
    forecast_lat: np.ndarray,  
    thresholds: List[float],
    prob_threshold: float = 0.5, # The probability cutoff to define a "forecast event"
    interpolation_method: str = 'linear' 
) -> Dict[float, Dict[str, float]]:
    """
    Computes dichotomous verification scores (POD, FAR, CSI, etc.) for each 
    precipitation threshold after interpolating probabilistic forecasts onto the 
    observed grid.

    Args:
        observed_prcp: The 2D array of observed precipitation.
        observed_lat: The 1D or 2D array of observed latitude coordinates.
        observed_lon: The 1D or 2D array of observed longitude coordinates.
        prob_forecasts: The 3D array of probabilistic forecasts (num_thresholds, Y, X).
        forecast_lon: The 2D array of forecast longitude (XLONG).
        forecast_lat: The 2D array of forecast latitude (XLAT).
        thresholds: List of precipitation thresholds corresponding to the 
                    first dimension of prob_forecasts.
        prob_threshold: The probability cutoff used to convert the forecast 
                        probability into a binary forecast event (0.0 to 1.0).
        interpolation_method: Scipy interpolation method ('linear', 'nearest', 'cubic').

    Returns:
        A dictionary where keys are the precipitation thresholds and values 
        are dictionaries of verification scores (POD, FAR, CSI, etc.).
    """
    
    results = {}
    
    # --- Define Interpolation Grid ---
    source_points = np.stack([forecast_lon.ravel(), forecast_lat.ravel()], axis=1)
    target_points = (observed_lon, observed_lat)
    
    # --- Iteration over all Precipitation Thresholds ---
    for idx, precip_threshold in enumerate(thresholds):
        
        # 1. Select the relevant 2D probability field
        original_prob_field = prob_forecasts[idx]

        # 2. Interpolate the probability field onto the observed grid
        source_values = original_prob_field.ravel()
        
        try:
            interpolated_prob_field = griddata(
                source_points, 
                source_values, 
                target_points, 
                method=interpolation_method
            )
        except Exception as e:
            print(f"Error during interpolation for threshold {precip_threshold}: {e}")
            continue

        # 3. Define the Observed Event (Truth)
        # Observed event occurs if precipitation exceeds the current threshold
        observed_event = (observed_prcp > precip_threshold).ravel()
        
        # 4. Define the Forecast Event
        # Forecast event occurs if probability exceeds the probability cutoff
        forecast_event = (interpolated_prob_field > prob_threshold).ravel()

        # 5. Mask Invalid Points
        # Exclude points where interpolation resulted in NaN (outside forecast domain)
        valid_mask = ~np.isnan(interpolated_prob_field).ravel()
        
        if np.sum(valid_mask) == 0:
            print(f"Warning: No valid data points for threshold {precip_threshold}. Skipping.")
            continue
            
        final_observed_event = observed_event[valid_mask]
        final_forecast_event = forecast_event[valid_mask]

        # 6. Compute Contingency Table (Confusion Matrix)
        # Confusion matrix order: [[TN, FP], [FN, TP]]
        # TP: True Positives (Forecast YES, Observed YES)
        # FN: False Negatives (Forecast NO, Observed YES)
        # FP: False Positives (Forecast YES, Observed NO)
        # TN: True Negatives (Forecast NO, Observed NO)
        
        tn, fp, fn, tp = confusion_matrix(
            y_true=final_observed_event, 
            y_pred=final_forecast_event
        ).ravel()
        
        # 7. Compute Dichotomous Scores
        
        hits = tp
        false_alarms = fp
        misses = fn
        correct_negatives = tn
        
        total_events = hits + misses
        total_forecasts = hits + false_alarms
        
        # Probability of Detection (Sensitivity)
        POD = hits / total_events if total_events > 0 else np.nan
        
        # False Alarm Ratio
        FAR = false_alarms / total_forecasts if total_forecasts > 0 else np.nan
        
        # Critical Success Index (Threat Score)
        CSI = hits / (hits + false_alarms + misses) if (hits + false_alarms + misses) > 0 else np.nan
        
        # Bias (Frequency Bias)
        BIAS = total_forecasts / total_events if total_events > 0 else np.nan
        
        # Accuracy (ACC)
        ACC = (hits + correct_negatives) / final_observed_event.size

        # 8. Store Results
        results[precip_threshold] = {
            "POD": POD,
            "FAR": FAR,
            "CSI": CSI,
            "BIAS": BIAS,
            "ACC": ACC,
            "TP": hits,
            "FP": false_alarms,
            "FN": misses,
            "TN": correct_negatives,
            "N_samples": final_observed_event.size
        }
        
    return results


def compute_reliability_curve(
    observed_prcp: np.ndarray, 
    observed_lat: np.ndarray, 
    observed_lon: np.ndarray, 
    prob_forecasts: np.ndarray, 
    forecast_lon: np.ndarray,  
    forecast_lat: np.ndarray,  
    thresholds: List[float],
    n_bins: int = 10, # Number of probability bins (e.g., 10 for 0.1 increments)
    interpolation_method: str = 'linear' 
) -> Dict[str, Any]:
    """
    Computes independent Reliability Curve coordinates for the event defined 
    by *each* precipitation threshold, handling interpolation and NaNs.

    Args:
        observed_prcp: The 2D array of observed precipitation.
        ... (coordinate and forecast arrays as before) ...
        thresholds: List of precipitation thresholds.
        n_bins: Number of bins to group the forecast probabilities (default: 10).
        interpolation_method: Scipy interpolation method.

    Returns:
        A dictionary containing lists of Observed Frequency and Mean Forecast Probability 
        arrays, indexed by the precipitation threshold.
    """
    
    # Pre-calculate the observed truth NaNs (needed for masking later)
    observed_nan_mask = np.isnan(observed_prcp.ravel())
    
    # --- Define Interpolation Grid ---
    source_points = np.stack([forecast_lon.ravel(), forecast_lat.ravel()], axis=1)
    target_points = (observed_lon, observed_lat)
    
    # Initialize results lists
    mean_forecast_prob_list = []
    observed_freq_list = []
    n_forecast_prob_list = []
    
    # Define the probability bins (e.g., [0.0, 0.1, 0.2, ..., 1.0])
    prob_bins = np.linspace(0.0, 1.0, n_bins + 1)
    
    # --- Iteration over all Precipitation Thresholds ---
    for idx, precip_threshold in enumerate(thresholds):
        
        print(f"Calculating Reliability Curve for event P_obs > {precip_threshold}...")

        # 1. Select and Interpolate the probability field
        original_prob_field = prob_forecasts[idx]
        source_values = original_prob_field.ravel()
        
        try:
            interpolated_prob_field = griddata(
                source_points, 
                source_values, 
                target_points, 
                method=interpolation_method
            )
        except Exception as e:
            print(f"Error during interpolation for threshold {precip_threshold}: {e}. Skipping.")
            continue

        # 2. Define the Observed Event (The "Truth")
        # Observed event occurs if precipitation exceeds the current threshold
        observed_event = (observed_prcp > precip_threshold).astype(int).ravel() # Convert boolean to 0/1

        # 3. Create the Combined NaN-Aware Mask
        forecast_nan_mask = np.isnan(interpolated_prob_field.ravel())
        valid_mask = ~(forecast_nan_mask | observed_nan_mask)
        
        if np.sum(valid_mask) <= n_bins: # Need at least one sample per bin to be ideal
            print(f"Warning: Insufficient valid data points ({np.sum(valid_mask)} remaining) for threshold {precip_threshold}. Skipping.")
            continue
            
        # 4. Apply the Mask
        final_observed_event = observed_event[valid_mask]
        final_forecast_probabilities = interpolated_prob_field.ravel()[valid_mask]

        # --- 5. Compute Reliability Curve Coordinates ---
        
        # Calculate which bin each valid probability score falls into
        # 'digitize' returns the index of the right bin edge. Subtract 1 for the bin index [0, N_bins-1]
        bin_indices = np.digitize(final_forecast_probabilities, prob_bins[1:])
        
        # Initialize arrays for results for this threshold
        mean_prob_in_bin = np.full(n_bins, np.nan)
        observed_frequency_in_bin = np.full(n_bins, np.nan)
        n_prob_in_bin = np.full( n_bins,0)
        
        # Iterate through each bin index
        for bin_idx in range(n_bins):
            # Select all forecast/observation pairs belonging to this bin
            bin_mask = (bin_indices == bin_idx)
            
            if np.sum(bin_mask) > 0:
                # Mean Forecast Probability (X-coordinate of the reliability curve)
                mean_prob_in_bin[bin_idx] = final_forecast_probabilities[bin_mask].mean()
                
                # Observed Frequency (Y-coordinate of the reliability curve)
                # Observed Frequency = (Sum of Observed Events) / (Number of Samples in Bin)
                observed_frequency_in_bin[bin_idx] = final_observed_event[bin_mask].mean()

                n_prob_in_bin[bin_idx] = np.sum( bin_mask )
         
        mean_forecast_prob_list.append( copy.deepcopy( mean_prob_in_bin ) )
        observed_freq_list.append( copy.deepcopy( observed_frequency_in_bin ) )
        n_forecast_prob_list.append( copy.deepcopy( n_prob_in_bin ) )  
    # --- 6. Return Results ---
    
    return np.array( mean_forecast_prob_list ) , np.array( observed_freq_list ) , np.array( n_forecast_prob_list ) , n_bins , np.array( thresholds )
#        "precipitation_thresholds": np.array(thresholds),
#        "n_bins": n_bins,
#        "Mean_Forecast_Probabilities": mean_forecast_prob_list, # X-coordinates
#        "Observed_Frequencies": observed_freq_list,             # Y-coordinates
#        "description": "Reliability curve coordinates computed independently for each precipitation threshold."
#    }
