import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pdb 
from tqdm import tqdm


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



def plotThisElement(fine_array_2d, element, v_min_max, name):
    lats = np.linspace(LAT_MIN, LAT_MAX, LAT_PIXELS)
    lons = np.linspace(LON_MIN, LON_MAX, LON_PIXELS)
    LON, LAT = np.meshgrid(lons, lats)

    fig, ax = plt.subplots(figsize=(12,6))
    norm = Normalize(vmin=v_min_max[0], vmax=v_min_max[1], clip=True)
    print("Minimum value in the array: ", np.nanmin(fine_array_2d))
    print("Maximum value in the array: ", np.nanmax(fine_array_2d))
    print("Median value in the array: ", np.nanmedian(fine_array_2d))
    print("Standard deviation of the array: ", np.nanstd(fine_array_2d))
    cmap = plt.cm.viridis
    im = ax.pcolormesh(LON, LAT, fine_array_2d, cmap=cmap, norm=norm, shading='auto')
    # plt.colorbar(im, ax=ax, label='Amplitude ratio vs Si')
    # ax.set_xlabel('Longitude')
    # ax.set_ylabel('Latitude')
    # ax.set_title(f'Average {element} Ratio Map at sub-pixel resolution')
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')
    # plt.show()
    fig.savefig(f"SubpixelImages/{element}_high_res_{name}.png", bbox_inches='tight', pad_inches=0, transparent=True, dpi=1000)
    print(f"Saved the high resolution map for {element} with {name} interlacing function.")
    return None 


def getElementData(element):
    ratio_col = "{0}_amplitude_ratio_vs_Si".format(element)

    #########################################
    # Load the CSV
    #########################################

    # csv_file = "Al_subset_temp.csv"  # Replace with your CSV path
    csv_file = "Subsets/{0}_subset.csv".format(element)  # Replace with your CSV path

    df = pd.read_csv(csv_file)

    latitudes = df['BORE_LAT'].values
    longitudes = df['BORE_LON'].values
    amplitude_ratio = df[ratio_col].values
    gmr2 = np.sqrt(df['Si_r2_score'] * df[f'{element}_r2_score'])
    return latitudes, longitudes, amplitude_ratio, gmr2



def processElementData(element, latitudes, longitudes, amplitude_ratio, gmr2, LAT_PIXELS, LON_PIXELS, interlacingParams):

    #########################################
    # Initialize Arrays
    #########################################
    """
    fine_array dimensions: [LAT_PIXELS, LON_PIXELS, NUM_PASSES]
    Initialize with NaN. Slice 0: final average, slices 1..: individual patch values
    """
    fine_array_2d = np.full((LAT_PIXELS, LON_PIXELS), np.nan, dtype=float)

    # weight_array tracks how many patches have covered each pixel
    weight_array = np.zeros((LAT_PIXELS, LON_PIXELS), dtype=float)

    # interlacing parameters
    sigma, func = interlacingParams

    #########################################
    # Process each patch (row in CSV)
    #########################################

    for iterator_ in tqdm(range(len(amplitude_ratio)), desc=f"Generating high-res maps for {element}", position=0, leave=True):
        lat, lon, amp_ratio, gmeanr2 = latitudes[iterator_], longitudes[iterator_], amplitude_ratio[iterator_], gmr2[iterator_]
        
        # Calculate patch boundaries
        delta_latitude = (180.0 / np.pi) * (12.5 / 1737.5)
        margin = 0.01

        # Handle near-pole condition
        if abs(lat - 90) < margin:
            lat = 90 - margin

        delta_longitude = delta_latitude / np.cos(np.deg2rad(lat))
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

        for i in range(lat_start, lat_end):
            for j in range(lon_start, lon_end):
                current_count = weight_array[i, j]
                distance = moonRadius*np.sqrt((np.deg2rad(index_to_lat(i)-lat))**2 + np.cos(np.deg2rad(lat))**2*np.deg2rad(index_to_lon(j)-lon)**2)
                # print(distance)
                # if input().lower() == "n": exit()
                # weight = linear_drop(distance, 12.5)
                # weight = square_root_drop(distance, 12.5)
                # weight = tanh_drop(distance, 20)
                # weight = tanh_drop(distance, 20) * tanh_drop(abs(amp_ratio - 1.5), 0.25)
                # weight = tanh_drop(distance, 20) * tanh_drop(abs(amp_ratio - 0.5), 0.5)
                weight = func(distance, sigma)# * gmeanr2

                # weight = myGaussian(distance)
                fine_array_2d[i,j] = np.nansum([fine_array_2d[i,j], amp_ratio*weight], axis=0)
                weight_array[i,j] += weight

    return fine_array_2d / weight_array


def main():

    #########################################
    # Configuration
    #########################################

    global LAT_MAX, LAT_MIN, LON_MIN, LON_MAX, LAT_PIXELS, LON_PIXELS, lat_res, lon_res
    global moonRadius, ch2UnitLength # ch2UnitLength is the unit of length for chandrayaan 2 viz. 12.5km

    LAT_MIN, LAT_MAX = -90, 90
    LON_MIN, LON_MAX = -180, 180
    LAT_PIXELS = 1800
    LON_PIXELS = 3600

    # interlacing parameters
    ch2UnitLength = 12.5
    moonRadius = 0.5*3475
    sigmaList = [ch2UnitLength / 30, ch2UnitLength, ch2UnitLength]
    intFuncList = [gaussian, linear_drop, tanh_drop]
    name_array = ["Gaussian", "Linear", "Tanh"]
    interlacingParamList = [
        [sigma, func_] for sigma, func_ in zip(sigmaList, intFuncList)
        # [sigma, lambda x: func_(x, sigma)] for sigma, func_ in zip(sigmaList, intFuncList)
    ]

    # Derived resolutions (approximately 0.08 degrees per pixel)
    lat_res = (LAT_MAX - LAT_MIN) / LAT_PIXELS
    lon_res = (LON_MAX - LON_MIN) / LON_PIXELS

    elements = "Na Mg Fe Si Ca Ti Al O".split(' ')
    # plot_elements = [False, True, False, False, False, False, False, True]
    plot_elements = [False] + [True] + [False]*4 + [True] + [False]
    v_min_max_array = [[0.6, 1.88], [0.25, 1.25], [0.2, 0.3], [0, 2], [0, 0.4], [0.015, 0.045], [0.5, 2], [0, 16]]
    for element in elements:
        if not plot_elements[elements.index(element)]:
            continue
        latitudes, longitudes, amplitude_ratio, gmr2 = getElementData(element=element)

        for params in interlacingParamList:
            subpixelArray = processElementData(element, latitudes, longitudes, amplitude_ratio, gmr2, LAT_PIXELS, LON_PIXELS, params)
            plotThisElement(subpixelArray, element, v_min_max_array[elements.index(element)],  name=f"drop={name_array[interlacingParamList.index(params)]}")

    print("High res maps generated for all elements successfully!")

    return None 


if __name__ == "__main__":
    main()
