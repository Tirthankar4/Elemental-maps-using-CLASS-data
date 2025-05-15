import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from shapely.geometry import box
from tqdm import tqdm
import math
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

#########################################
# Helper Functions
#########################################

def lat_to_index(lat):
    """
    Convert latitude to grid index.
    """
    return int((lat - LAT_MIN) / lat_res)

def lon_to_index(lon):
    """
    Convert longitude to grid index.
    """
    return int((lon - LON_MIN) / lon_res)

def index_to_lat(idx):
    """
    Convert grid index to latitude.
    """
    return idx * lat_res + LAT_MIN + lat_res / 2

def index_to_lon(idx):
    """
    Convert grid index to longitude.
    """
    return idx * lon_res + LON_MIN + lon_res / 2

def gaussian(x, sigma):
    """
    Gaussian weighting function.
    """
    return np.exp(-0.5 * (x / sigma) ** 2)

def linear_drop(x, cutoff):
    """
    Linear drop-off function.
    """
    if x <= cutoff:
        return 1 - (x / cutoff)
    else:
        return 0

def square_root_drop(x, a):
    """
    Square root drop-off function.
    """
    if x > a:
        return 0
    return math.sqrt((a - x) / a)

def tanh_drop(x, a):
    """
    Hyperbolic tangent drop-off function.
    """
    if x > a:
        return 0
    return math.tanh((a - x) / a)

def plotThisElement(fine_array_2d, element, v_min_max, name):
    """
    Plot and save the high-resolution map for a given element.
    """
    lats = np.linspace(LAT_MIN, LAT_MAX, LAT_PIXELS)
    lons = np.linspace(LON_MIN, LON_MAX, LON_PIXELS)
    LON, LAT = np.meshgrid(lons, lats)

    fig, ax = plt.subplots(figsize=(12,6))
    norm = Normalize(vmin=v_min_max[0], vmax=v_min_max[1], clip=True)
    print(f"Plotting {element} with {name} interlacing function.")
    print("Minimum value in the array: ", np.nanmin(fine_array_2d))
    print("Maximum value in the array: ", np.nanmax(fine_array_2d))
    print("Median value in the array: ", np.nanmedian(fine_array_2d))
    print("Standard deviation of the array: ", np.nanstd(fine_array_2d))
    cmap = plt.cm.viridis
    im = ax.pcolormesh(LON, LAT, fine_array_2d, cmap=cmap, norm=norm, shading='auto')
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')
    fig.savefig(f"SubpixelImages/{element}_high_res_{name}.png", bbox_inches='tight', pad_inches=0, transparent=True, dpi=1000)
    print(f"Saved the high resolution map for {element} with {name} interlacing function.")
    plt.close(fig)
    return None 

def getElementData(element):
    """
    Load and extract data for a specific element from the CSV.
    """
    ratio_col = "{0}_amplitude_ratio_vs_Si".format(element)

    #########################################
    # Load the CSV
    #########################################

    csv_file = "Subsets/{0}_subset.csv".format(element)  # Replace with your CSV path

    df = pd.read_csv(csv_file)

    latitudes = df['BORE_LAT'].values
    longitudes = df['BORE_LON'].values
    amplitude_ratio = df[ratio_col].values
    # Compute geometric mean of R² scores
    r2_col_si = 'Si_r2_score'  # Adjust if different
    r2_col_element = f'{element}_r2_score'  # Adjust if different
    gmr2 = np.sqrt(df[r2_col_si] * df[r2_col_element])
    return latitudes, longitudes, amplitude_ratio, gmr2

def processElementData(element, latitudes, longitudes, amplitude_ratio, gmr2, LAT_PIXELS, LON_PIXELS, interlacingParams):
    """
    Process data for a specific element to generate a high-resolution map.
    """
    #########################################
    # Initialize Arrays
    #########################################
    """
    fine_array dimensions: [LAT_PIXELS, LON_PIXELS]
    Initialize with zeros. We'll accumulate weighted amplitude ratios.
    """
    fine_array_2d = np.zeros((LAT_PIXELS, LON_PIXELS), dtype=float)
    accumulated_weights = np.zeros((LAT_PIXELS, LON_PIXELS), dtype=float)

    # interlacing parameters
    sigma, func = interlacingParams  # sigma: Gaussian sigma, func: drop function

    #########################################
    # Process each patch (row in CSV)
    #########################################

    for iterator_ in tqdm(range(len(amplitude_ratio)), desc=f"Processing patches for {element}", position=0, leave=True):
        lat, lon, amp_ratio, gmeanr2 = latitudes[iterator_], longitudes[iterator_], amplitude_ratio[iterator_], gmr2[iterator_]

        # Calculate patch boundaries
        delta_latitude = (180.0 / math.pi) * (12.5 / 1737.5)  # Convert km to degrees
        margin = 0.01

        # Handle near-pole condition
        if abs(lat - 90) < margin:
            lat = 90 - margin
        elif abs(lat + 90) < margin:
            lat = -90 + margin

        # Compute delta_longitude, handling poles
        if math.cos(math.radians(lat)) != 0:
            delta_longitude = delta_latitude / math.cos(math.radians(lat))
        else:
            delta_longitude = 360  # At the poles

        lat_min = lat - delta_latitude / 2
        lat_max = lat + delta_latitude / 2
        lon_min = lon - delta_longitude / 2
        lon_max = lon + delta_longitude / 2

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

        # Compute Center Grid Indices
        lat_center_idx = lat_to_index(lat)
        lon_center_idx = lon_to_index(lon)

        # Define Patch Geometry using Shapely
        patch_box = box(lon_min, lat_min, lon_max, lat_max)

        # Iterate Over Overlapping Fine Grid Pixels
        for i in range(lat_start, lat_end):
            for j in range(lon_start, lon_end):
                # Define Fine Grid Pixel Geometry
                pixel_lon_min = LON_MIN + j * lon_res
                pixel_lat_min = LAT_MIN + i * lat_res
                pixel_lon_max = pixel_lon_min + lon_res
                pixel_lat_max = pixel_lat_min + lat_res
                pixel_box = box(pixel_lon_min, pixel_lat_min, pixel_lon_max, pixel_lat_max)

                # Compute Fractional Overlap using Shapely
                intersection = patch_box.intersection(pixel_box)
                if intersection.is_empty:
                    fractional_overlap = 0.0
                else:
                    overlap_area = intersection.area
                    pixel_area = lon_res * lat_res
                    fractional_overlap = overlap_area / pixel_area  # Between 0 and 1

                if fractional_overlap <= 0:
                    continue  # No contribution

                # Compute Distance from Patch Center in Grid Pixels
                # Convert lat/lon difference to radians for accurate distance
                # Use the grid indices to compute distance in grid units                                                                                                                                                                                                                             
                # Each grid step is approximately lat_res degrees
                lat_diff_deg = index_to_lat(i) - lat
                lon_diff_deg = index_to_lon(j) - lon
                lat_diff_rad = math.radians(lat_diff_deg)
                lon_diff_rad = math.radians(lon_diff_deg)
                distance = math.sqrt((lat_diff_rad)**2 + (math.cos(math.radians(lat)) * lon_diff_rad)**2)

                # Compute Gaussian Weight
                g_weight = gaussian(distance, sigma)

                # Compute Final Weight (Geometric Mean of Gaussian and R²)
                #final_weight = math.sqrt(g_weight * gmeanr2)
                final_weight = g_weight * gmeanr2

                # Accumulate Weighted Amplitude Ratio and Weights
                fine_array_2d[i, j] += amp_ratio * final_weight * fractional_overlap
                accumulated_weights[i, j] += final_weight * fractional_overlap
                print(fractional_overlap)

    # Compute Final Averaged Map
    with np.errstate(divide='ignore', invalid='ignore'):
        final_map = fine_array_2d / accumulated_weights
        final_map[np.isnan(final_map)] = 0  # Assign 0 to pixels with no coverage

    return final_map

def main():

    #########################################
    # Configuration
    #########################################

    global LAT_MAX, LAT_MIN, LON_MIN, LON_MAX, LAT_PIXELS, LON_PIXELS, lat_res, lon_res

    LAT_MIN, LAT_MAX = -90, 90
    LON_MIN, LON_MAX = -180, 180
    LAT_PIXELS = 1800
    LON_PIXELS = 3600

    # Derived resolutions (approximately 0.08 degrees per pixel)
    lat_res = (LAT_MAX - LAT_MIN) / LAT_PIXELS
    lon_res = (LON_MAX - LON_MIN) / LON_PIXELS

    # Interlacing parameters
    # Define sigma and corresponding drop functions
    sigmaList = [2, 12.5, 12.5]  # Example sigma values; adjust as needed
    intFuncList = [gaussian, linear_drop, tanh_drop]
    name_array = ["Gaussian", "Linear", "Tanh"]
    interlacingParamList = [
        [sigma, func_] for sigma, func_ in zip(sigmaList, intFuncList)
    ]

    # Elements to process
    elements = "Na Mg Fe Si Ca Ti Al O".split(' ')
    # Define which elements to plot (index corresponds to elements list)
    plot_elements = [False, True, False, False, False, False, False, True]
    # Define v_min and v_max for plotting (adjust based on your data)
    v_min_max_array = [
        [0.6, 1.88],    # Na
        [0.25, 1.25],   # Mg
        [0.2, 0.3],     # Fe
        [0, 2],         # Si
        [0, 0.4],       # Ca
        [0.015, 0.045], # Ti
        [0.5, 2],       # Al
        [0, 16]         # O
    ]

    for idx, element in enumerate(elements):
        if not plot_elements[idx]:
            continue
        latitudes, longitudes, amplitude_ratio, gmr2 = getElementData(element=element)

        for params in interlacingParamList:
            sigma, func_ = params
            final_map = processElementData(
                element, latitudes, longitudes, amplitude_ratio, gmr2,
                LAT_PIXELS, LON_PIXELS, params
            )
            plotThisElement(
                final_map, element, v_min_max_array[idx],
                name=f"drop={name_array[interlacingParamList.index(params)]}"
            )

    print("High-res maps generated for all elements successfully!")

    return None 

if __name__ == "__main__":
    main()
