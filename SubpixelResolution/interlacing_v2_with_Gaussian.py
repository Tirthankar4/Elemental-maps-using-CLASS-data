from re import sub
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pdb 
from tqdm import tqdm
from multiprocessing import Pool
#########################################
# Configuration
#########################################

LAT_MIN, LAT_MAX = -90, 90
LON_MIN, LON_MAX = -180, 180
LAT_PIXELS = int(1800)
LON_PIXELS = int(3600)
# NUM_PASSES = 10000  # 0th for average, 1.. for patches


# Derived resolutions (approximately 0.08 degrees per pixel)
lat_res = (LAT_MAX - LAT_MIN) / LAT_PIXELS
lon_res = (LON_MAX - LON_MIN) / LON_PIXELS

element = "Al"  # Choose your element
ratio_col = "{0}_amplitude_ratio_vs_Si".format(element)

#########################################
# Helper Functions
#########################################

def lat_to_index(lat):
    # Map latitude from [-90,90] to [0, LAT_PIXELS-1]
    return int((lat - LAT_MIN) / lat_res)

def lon_to_index(lon):
    # Map longitude from [-180,180] to [0, LON_PIXELS-1]
    return int((lon - LON_MIN) / lon_res)

def index_to_lat(idx):
    # Map latitude from [-90,90] to [0, LAT_PIXELS-1]
    return idx*lat_res + LAT_MIN

def index_to_lon(idx):
    # Map longitude from [-180,180] to [0, LON_PIXELS-1]
    return idx*lon_res + LON_MIN

def gaussian(x, sigma):
    return np.exp(-0.5*(x/sigma)**2)

def linear_drop(x, cutoff):
    if (x <= cutoff):
        value = 1 - (x / cutoff)
    else:
        value = 0
    return value

def square_root_drop(x, a):
    if (x > a):
        return 0
    return np.sqrt((a - x) / a)

def tanh_drop(x, a):
    if (x > a):
        return 0
    return np.tanh((a - x) / a)

#########################################
# Initialize Arrays
#########################################

"""
fine_array dimensions: [LAT_PIXELS, LON_PIXELS, NUM_PASSES]
Initialize with NaN. Slice 0: final average, slices 1..: individual patch values
"""
fine_array_2d = np.full((LAT_PIXELS, LON_PIXELS), np.nan, dtype=float)
# fine_array_2d = np.full((LAT_PIXELS, LON_PIXELS), 0, dtype=float)

# weight_array tracks how many patches have covered each pixel
weight_array = np.zeros((LAT_PIXELS, LON_PIXELS), dtype=float)

#########################################
# Load the CSV
#########################################

# csv_file = "Al_subset_temp.csv"  # Replace with your CSV path
csv_file = "{0}_subset.csv".format(element)  # Replace with your CSV path

df = pd.read_csv(csv_file)

latitudes = df['BORE_LAT'].values
longitudes = df['BORE_LON'].values
amplitude_ratio = df[ratio_col].values


# Gaussian parameters
R = 12.5
radius = 0.5*3475
sigma = R / 30
myGaussian = lambda x: gaussian(x, sigma=sigma)

#########################################
# Process each patch (row in CSV)
#########################################

# tqdm(range(outer_loop_iterations), desc="Outer Loop", position=0, leave=True):

for iterator_ in tqdm(range(len(amplitude_ratio)), desc="Generating high-res maps", position=0, leave=True):
    lat, lon, amp_ratio = latitudes[iterator_], longitudes[iterator_], amplitude_ratio[iterator_]
    # Calculate patch boundaries
    delta_latitude = np.rad2deg(12.5*2/3475)
    margin = 0.01

    # Handle near-pole condition
    if abs(lat - 90) < margin:
        lat = 90 - margin

    delta_longitude = delta_latitude / np.cos(np.deg2rad(lat))
    lat_min = lat - delta_latitude
    lat_max = lat + delta_latitude
    lon_min = lon - delta_longitude
    lon_max = lon + delta_longitude

    # lat_min = lat - delta_latitude/2
    # lat_max = lat + delta_latitude/2
    # lon_min = lon - delta_longitude/2
    # lon_max = lon + delta_longitude/2

    # Convert patch boundaries to indices in the fine grid
    lat_start = lat_to_index(lat_min)
    lat_end = lat_to_index(lat_max) + 1  # +1 to include the boundary pixel
    lon_start = lon_to_index(lon_min)
    lon_end = lon_to_index(lon_max) + 1

    # Bound checks
    lat_start = max(lat_start, 0)
    lat_end = min(lat_end, LAT_PIXELS)
    lon_start = max(lon_start, 0)
    lon_end = min(lon_end, LON_PIXELS)

    """
    Fill the patch area in the fine array
    For each pixel in the patch, store amp_ratio in the next available slot
    """
    for i in range(lat_start, lat_end):
        for j in range(lon_start, lon_end):
            current_count = weight_array[i, j]
            # pdb.set_trace()
            distance = radius*np.sqrt((np.deg2rad(index_to_lat(i)-lat))**2 + np.cos(np.deg2rad(lat))**2*np.deg2rad(index_to_lon(j)-lon)**2)
            # print(distance)
            # if input().lower() == "n": exit()
            # weight = linear_drop(distance, 12.5)
            # weight = square_root_drop(distance, 12.5)
            weight = tanh_drop(distance, 20)
            # weight = tanh_drop(distance, 20) * tanh_drop(abs(amp_ratio - 1.5), 0.25)
            # weight = tanh_drop(distance, 20) * tanh_drop(abs(amp_ratio - 0.5), 0.5)

            # weight = myGaussian(distance)
            fine_array_2d[i,j] = np.nansum([fine_array_2d[i,j], amp_ratio*weight], axis=0)
            weight_array[i,j] += weight

subpixelArray = fine_array_2d / weight_array


lats = np.linspace(LAT_MIN, LAT_MAX, LAT_PIXELS)
lons = np.linspace(LON_MIN, LON_MAX, LON_PIXELS)
LON, LAT = np.meshgrid(lons, lats)

for i in [subpixelArray]: #, fine_array_2d, weight_array]:
    fig, ax = plt.subplots(figsize=(12,6))
    norm = Normalize(vmin=0, vmax=1.75, clip=True)
    cmap = plt.cm.viridis
    im = ax.pcolormesh(LON, LAT, i, cmap=cmap, norm=norm, shading='auto')
    plt.colorbar(im, ax=ax, label='Amplitude ratio vs Si')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(f'Average {element} Ratio Map at 0.08° resolution')
    ax.set_aspect('equal', adjustable='box')
plt.show()
