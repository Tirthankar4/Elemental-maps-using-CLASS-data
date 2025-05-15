import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from shapely.geometry import box
import rasterio
from rasterio.transform import from_origin
import math
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# Configuration Parameters
LAT_MIN, LAT_MAX = -90, 90
LON_MIN, LON_MAX = -180, 180
FINE_RESOLUTION = 0.08  # degrees per pixel

# Calculate grid dimensions
LAT_PIXELS = int((LAT_MAX - LAT_MIN) / FINE_RESOLUTION)  # 2249
LON_PIXELS = int((LON_MAX - LON_MIN) / FINE_RESOLUTION)  # 4499

# Gaussian Parameters
GAUSSIAN_SIGMA = 2  # Standard deviation in fine grid pixels

# File Paths
CSV_FILE_PATH = "your_input_file.csv"  # Replace with your actual CSV file path

# Element Configuration
ELEMENT = "Fe"  # Example element
R2_COL = f"{ELEMENT}_r2_score"  # Replace with your actual R² column name
AMPLITUDE_RATIO_COL = f"{ELEMENT}_amplitude_ratio_vs_Si"  # Replace with your actual amplitude ratio column name

# Initialize Accumulation Arrays
accumulated_data = np.zeros((LAT_PIXELS, LON_PIXELS), dtype=np.float64)
accumulated_weights = np.zeros((LAT_PIXELS, LON_PIXELS), dtype=np.float64)

# Load CSV Data
df = pd.read_csv(CSV_FILE_PATH)

# Extract Relevant Columns
latitudes = df['BORE_LAT'].values
longitudes = df['BORE_LON'].values
amplitude_ratios = df[AMPLITUDE_RATIO_COL].values
r2_values = df[R2_COL].values

# Define Gaussian Weight Function
def gaussian_weight(distance, sigma=GAUSSIAN_SIGMA):
    return np.exp(- (distance ** 2) / (2 * sigma ** 2))

# Define Function to Convert Latitude to Grid Index
def lat_to_index(lat):
    return int((lat - LAT_MIN) / FINE_RESOLUTION)

# Define Function to Convert Longitude to Grid Index
def lon_to_index(lon):
    return int((lon - LON_MIN) / FINE_RESOLUTION)

# Define Function to Compute Distance in Grid Pixels
def compute_distance(lat_center_idx, lon_center_idx, i, j):
    return math.sqrt((i - lat_center_idx) ** 2 + (j - lon_center_idx) ** 2)

# Optional: Define the transform for shapely (if needed for area calculations)
# Not strictly necessary here since we're dealing with degree-based grids
transform = from_origin(LON_MIN, LAT_MAX, FINE_RESOLUTION, FINE_RESOLUTION)

# Processing Each Patch
for idx in range(len(df)):
    lat_center = latitudes[idx]
    lon_center = longitudes[idx]
    amp_ratio = amplitude_ratios[idx]
    r2_weight = r2_values[idx]
    
    # Handle Near-Pole Condition to avoid division by zero
    margin = 0.01
    if abs(lat_center - 90) < margin:
        lat_center = 90 - margin
    elif abs(lat_center + 90) < margin:
        lat_center = -90 + margin
    
    # Compute Patch Boundaries
    delta_latitude = 0.4  # degrees
    try:
        delta_longitude = 0.4 / math.cos(math.radians(lat_center))
    except ZeroDivisionError:
        delta_longitude = 360  # At the poles, longitude spans all
    
    lat_min = lat_center - delta_latitude / 2
    lat_max = lat_center + delta_latitude / 2
    lon_min = lon_center - delta_longitude / 2
    lon_max = lon_center + delta_longitude / 2
    
    # Convert Boundaries to Grid Indices
    lat_start_idx = lat_to_index(lat_min)
    lat_end_idx = lat_to_index(lat_max) + 1  # +1 to include the upper boundary
    lon_start_idx = lon_to_index(lon_min)
    lon_end_idx = lon_to_index(lon_max) + 1
    
    # Boundary Checks
    lat_start_idx = max(lat_start_idx, 0)
    lat_end_idx = min(lat_end_idx, LAT_PIXELS)
    lon_start_idx = max(lon_start_idx, 0)
    lon_end_idx = min(lon_end_idx, LON_PIXELS)
    
    # Compute Center Grid Indices
    lat_center_idx = lat_to_index(lat_center)
    lon_center_idx = lon_to_index(lon_center)
    
    # Define Patch Geometry
    patch_box = box(lon_min, lat_min, lon_max, lat_max)
    
    # Iterate Over Overlapping Fine Grid Pixels
    for i in range(lat_start_idx, lat_end_idx):
        for j in range(lon_start_idx, lon_end_idx):
            # Define Fine Grid Pixel Geometry
            pixel_lon_min = LON_MIN + j * FINE_RESOLUTION
            pixel_lat_min = LAT_MIN + i * FINE_RESOLUTION
            pixel_lon_max = pixel_lon_min + FINE_RESOLUTION
            pixel_lat_max = pixel_lat_min + FINE_RESOLUTION
            pixel_box = box(pixel_lon_min, pixel_lat_min, pixel_lon_max, pixel_lat_max)
            
            # Compute Fractional Overlap using Shapely
            intersection = patch_box.intersection(pixel_box)
            if intersection.is_empty:
                fractional_overlap = 0.0
            else:
                overlap_area = intersection.area
                pixel_area = FINE_RESOLUTION ** 2
                fractional_overlap = overlap_area / pixel_area  # Between 0 and 1
            
            if fractional_overlap <= 0:
                continue  # No contribution
            
            # Compute Distance from Patch Center in Grid Pixels
            distance = compute_distance(lat_center_idx, lon_center_idx, i, j)
            
            # Compute Gaussian Weight
            fracOverlap = 0
            s = 2
            g_weight = gaussian_weight(distance)
            final_weight = math.sqrt(g_weight * r2_weight)
            old_weights = accumulated_weights[i, j]
            accumulated_weights[i, j] += final_weight * fractional_overlap                                                                                                                                      
            pixelVal = (amp_ratio*fracOverlap*final_weight*s**2 + accumulated_data[i, j]*old_weights)/accumulated_weights[i,j]
            accumulated_data[i,j]=pixelVal
            # accumulated_data[i, j] += amp_ratio * final_weight * fractional_overlap

            # Compute Final Weight (Geometric Mean of Gaussian and R²)
            # final_weight = math.sqrt(g_weight * r2_weight)
            
            # Accumulate Weighted Amplitude Ratio and Weights
            # accumulated_data[i, j] += amp_ratio * final_weight * fractional_overlap
            # accumulated_weights[i, j] += final_weight * fractional_overlap

# Compute Final Averaged Map
with np.errstate(divide='ignore', invalid='ignore'):
    final_map = accumulated_data / accumulated_weights
    final_map[np.isnan(final_map)] = 0  # Assign 0 to pixels with no coverage

# Generate Coordinate Arrays for Plotting
lats = LAT_MIN + (np.arange(LAT_PIXELS) + 0.5) * FINE_RESOLUTION
lons = LON_MIN + (np.arange(LON_PIXELS) + 0.5) * FINE_RESOLUTION
LON, LAT = np.meshgrid(lons, lats)

# Plot the Final Map
fig, ax = plt.subplots(figsize=(12, 6))
norm = Normalize(vmin=np.nanmin(final_map), vmax=np.nanmax(final_map), clip=True)
cmap = plt.cm.viridis
im = ax.pcolormesh(LON, LAT, final_map, cmap=cmap, norm=norm, shading='auto')
plt.colorbar(im, ax=ax, label='Amplitude ratio vs Si')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(f'Custom Drizzle-Like Weighted Average {ELEMENT} Ratio Map at 0.08° Resolution')
ax.set_aspect('equal', adjustable='box')
plt.show()
