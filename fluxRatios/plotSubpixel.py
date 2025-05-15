import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from matplotlib import cm
import os

def load_and_process_data(element):
    # Import the CSV file
    file_path = f"Subsets/{element}_subset.csv"
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} was not found.")

    data = pd.read_csv(file_path)

    # Ensure the relevant headers exist
    required_columns = ['BORE_LAT', 'BORE_LONG', f'{element}_amplitude_ratio_vs_Si']
    if not all(col in data.columns for col in required_columns):
        raise ValueError(f"Missing one or more required columns: {required_columns}")

    return data

def plot_heatmap_with_variable_pixels(data, element, base_image=None):
    # Define map resolution parameters (latitudinal dependent)
    latitudes = data['BORE_LAT'].values
    longitudes = data['BORE_LONG'].values
    amplitude_ratios = data[f'{element}_amplitude_ratio_vs_Si'].values

    # Calculate pixel size based on latitude (0.4 degrees / cos(latitude))
    pixel_sizes = 0.4 / np.cos(np.radians(latitudes))

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))

    # Scatter the data onto the map with the calculated pixel size
    sc = ax.scatter(longitudes, latitudes, c=amplitude_ratios, s=pixel_sizes*50, cmap='viridis', 
                    norm=Normalize(vmin=min(amplitude_ratios), vmax=max(amplitude_ratios)), edgecolors="k", alpha=0.7)

    # Set axis labels
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(f'{element} Amplitude Ratio vs Si Heatmap')

    # Add colorbar for the amplitude ratio
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical')
    cbar.set_label(f'{element} Amplitude Ratio')

    # Handle missing data by setting it as transparent
    transparent = np.where(np.isnan(amplitude_ratios), np.nan, amplitude_ratios)

    # Overlay the base image if provided
    if base_image:
        img = plt.imread(base_image)
        ax.imshow(img, extent=[min(longitudes), max(longitudes), min(latitudes), max(latitudes)], alpha=0.3)

    # Display the map
    plt.show()

# Example usage
element = "Al"  # Replace with your actual element name
base_image_path = "path_to_base_map_image.jpg"  # Replace with your base map image path, or None if no base map

# Load the data
data = load_and_process_data(element)

# Plot the heatmap with optional base map
plot_heatmap_with_variable_pixels(data, element, base_image=base_image_path)
