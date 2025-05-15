import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pdb 


#########################################
# Helper Functions
#########################################

def lat_to_index(lat):
    # Map latitude from [-90,90] to [0, LAT_PIXELS-1]
    return int((lat - LAT_MIN) / lat_res)

def lon_to_index(lon):
    # Map longitude from [-180,180] to [0, LON_PIXELS-1]
    return int((lon - LON_MIN) / lon_res)



def plotThisElement(fine_array_2d, element):
    lats = np.linspace(LAT_MIN, LAT_MAX, LAT_PIXELS)
    lons = np.linspace(LON_MIN, LON_MAX, LON_PIXELS)
    LON, LAT = np.meshgrid(lons, lats)

    fig, ax = plt.subplots(figsize=(12,6))
    norm = Normalize(vmin=0, vmax=2, clip=True)
    cmap = plt.cm.viridis
    im = ax.pcolormesh(LON, LAT, fine_array_2d, cmap=cmap, norm=norm, shading='auto')
    plt.colorbar(im, ax=ax, label='Amplitude ratio vs Si')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(f'Average {element} Ratio Map at 0.08° resolution')
    ax.set_aspect('equal', adjustable='box')
    plt.show()


def getElementData(element):
    ratio_col = "{0}_amplitude_ratio_vs_Si".format(element)

    #########################################
    # Load the CSV
    #########################################

    # csv_file = "Al_subset_temp.csv"  # Replace with your CSV path
    csv_file = "{0}_subset.csv".format(element)  # Replace with your CSV path

    df = pd.read_csv(csv_file)

    latitudes = df['BORE_LAT'].values
    longitudes = df['BORE_LON'].values
    amplitude_ratio = df[ratio_col].values

    return latitudes, longitudes, amplitude_ratio



def processElementData(latitudes, longitudes, amplitude_ratio, LAT_PIXELS, LON_PIXELS):

    #########################################
    # Initialize Arrays
    #########################################
    fine_array_2d = np.full((LAT_PIXELS, LON_PIXELS), np.nan, dtype=float)

    # fill_count tracks how many patches have covered each pixel
    fill_count = np.zeros((LAT_PIXELS, LON_PIXELS), dtype=int)


    #########################################
    # Process each patch (row in CSV)
    #########################################

    for lat, lon, amp_ratio in zip(latitudes, longitudes, amplitude_ratio):
        # Calculate patch boundaries
        delta_latitude = 0.4
        margin = 0.01

        # Handle near-pole condition
        if abs(lat - 90) < margin:
            lat = 90 - margin

        delta_longitude = 0.4 / np.cos(np.deg2rad(lat))
        lat_min = lat - delta_latitude/2
        lat_max = lat + delta_latitude/2
        lon_min = lon - delta_longitude/2
        lon_max = lon + delta_longitude/2

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

        # Fill the patch area in the fine array
        # For each pixel in the patch, store amp_ratio in the next available slot
        for i in range(lat_start, lat_end):
            for j in range(lon_start, lon_end):
                current_count = fill_count[i, j]
                fine_array_2d[i,j] = np.nansum([fine_array_2d[i,j]*fill_count[i,j], amp_ratio], axis=0)/(current_count+1)
                fill_count[i, j] += 1

    return fine_array_2d


def main():

    #########################################
    # Configuration
    #########################################

    global LAT_MAX, LAT_MIN, LON_MIN, LON_MAX, LAT_PIXELS, LON_PIXELS, lat_res, lon_res
    LAT_MIN, LAT_MAX = -90, 90
    LON_MIN, LON_MAX = -180, 180
    LAT_PIXELS = 1800
    LON_PIXELS = 1800*4489//2249

    # Derived resolutions (approximately 0.08 degrees per pixel)
    lat_res = (LAT_MAX - LAT_MIN) / LAT_PIXELS
    lon_res = (LON_MAX - LON_MIN) / LON_PIXELS

    elements = "Na Mg Fe Si Ca Ti Al O".split(' ')

    for element in elements:
        latitudes, longitudes, amplitude_ratio = getElementData(element=element)
        fine_array_2d = processElementData(latitudes, longitudes, amplitude_ratio, LAT_PIXELS, LON_PIXELS)
        plotThisElement(fine_array_2d, element)


if __name__ == "__main__":
    main()
