import csv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
import numpy as np
from matplotlib.colors import Normalize
import pdb
import os 
from PIL import Image


############# mare, highland detection ############
from mareDetection import mare_regions, in_mare
from highlandDetection import highland_regions, in_highland
#########################################

Image.MAX_IMAGE_PIXELS = None

def parse_and_generate_subsets(figs_and_axs, csv_file, csv_subset_files, mare_subset_files, highland_subset_files):
    # Read the CSV data
    df = pd.read_csv(csv_file)
    
    # Define the columns related to each element
    elements = ['Fe', 'Ti', 'Ca', 'Al', 'Mg', 'Na', 'O', 'Si']
    
    # Create an empty dictionary to store subsets for each element
    subsets = {element: [] for element in elements}
    mare_subsets = {mare: {element: [] for element in elements} for mare in mare_regions.keys()}
    highland_subsets = {highland: {element: [] for element in elements} for highland in highland_regions.keys()}
    
    # Iterate through each unique filename (assuming filenames are unique per set)
    unique_filenames = df['filename'].unique()
    
    for filename in unique_filenames:
        file_data = df[df['filename'] == filename]  # Filter data for the current file
        # Flag to keep track of which elements have already been added for this file
        found_elements = {element: False for element in elements}
        
        # Iterate over the rows in the file
        for index, row in file_data.iterrows():
            # Check for the first row where both Si_detection and element_detection are True
            for element in elements:
                if not found_elements[element]:  # Only process if the element has not been found
                    detection_col = f'{element}_detection'  # Detection column for the element
                    if row['Si_detection'] == True and row[detection_col] == True:
                        # Add the row to the corresponding subset
                        subsets[element].append(row)
                        found_elements[element] = True  # Mark the element as found    

                        # add the row to the corresponding mare subset
                        found, mareName = in_mare(row['BORE_LAT'], row['BORE_LON'])
                        if found:
                            mare_subsets[mareName][element].append(row)
                        else:
                            pass

                        # add the row to the corresponding highland subset
                        found, highlandName = in_highland(row['BORE_LAT'], row['BORE_LON'])
                        if found:
                            highland_subsets[highlandName][element].append(row)
                        else:
                            pass


    # Convert subsets from lists of rows to DataFrames for easy inspection
    subset_dfs = {element: pd.DataFrame(subsets[element]) for element in elements}
    mare_subsets_dfs = {mare: {element: pd.DataFrame(mare_subsets[mare][element]) for element in elements} for mare in mare_regions.keys()}
    highland_subsets_dfs = {highland: {element: pd.DataFrame(highland_subsets[highland][element]) for element in elements} for highland in highland_regions.keys()}

    # Save each mare subset to a separate CSV file
    for mare, element_subset in mare_subsets_dfs.items():
        for element, mare_subset_df in element_subset.items():
            if not subset_df.empty:
                file = mare_subset_files[mare][element]
                if not os.path.isfile(file):
                    mare_subset_df.to_csv(file, mode='a', header=True, index=False)
                else:
                    mare_subset_df.to_csv(file, mode='a', header=False, index=False)
            else:
                pass
                # print(f"Subset for {element} in {mare} is empty. Skipping plot.")

    # Save each highland subset to a separate CSV file
    for highland, element_subset in highland_subsets_dfs.items():
        for element, highland_subset_df in element_subset.items():
            if not subset_df.empty:
                file = highland_subset_files[highland][element]
                if not os.path.isfile(file):
                    highland_subset_df.to_csv(file, mode='a', header=True, index=False)
                else:
                    highland_subset_df.to_csv(file, mode='a', header=False, index=False)
            else:
                pass
                # print(f"Subset for {element} in {highland} is empty. Skipping plot.")


    # Save each element subset to a separate CSV file
    for data, canvas in zip(subset_dfs.items(), figs_and_axs):
        element, subset_df = data
        fig, ax = canvas
        if not subset_df.empty:  # Check if the subset is not empty
            # file_path = f'Subsets/{element}_subset.csv'
            file_path = f'/home/kaustav_b_ph.iitr/Plotting/Subsets/{element}_subset.csv'
            if not os.path.isfile(file_path):
                subset_df.to_csv(file_path, mode='a', header=True, index=False)
            else:
                subset_df.to_csv(file_path, mode='a', header=False, index=False)
            
            # Plot the shaded region for the non-empty subset
            plot_shaded_region_on_map(fig, ax, subset_df, element)
        else:
            pass
            # print(f"Subset for {element} is empty. Skipping plot.")
        
    # print(f"Subsets generated and saved for each element: {', '.join(elements)}")


def plot_shaded_region_on_map(fig, ax, subset_df, element):
    # Extract the relevant data
    latitudes = subset_df['BORE_LAT']
    longitudes = subset_df['BORE_LON']
    amplitude_ratio = subset_df[f'{element}_amplitude_ratio_vs_Si']
    
    # Normalize the amplitude ratios for coloring
    norm = Normalize(vmin=0, vmax=2, clip=True)
    margin = 0
    # norm = Normalize(vmin=amplitude_ratio.min()-margin, vmax=amplitude_ratio.max()+margin)
    cmap = plt.cm.viridis
    
    # Create ScalarMappable for colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # To ensure it works with colorbar
    
    # Iterate over each point and plot the shaded region (0.4°x0.4°)
    for lat, lon, amp_ratio in zip(latitudes, longitudes, amplitude_ratio):
        # Calculate the boundary of the square region (corresponding to 12.5 km side square)
        delta_latitude = 0.4 

        margin = 0.01
        if abs(lat - 90) < margin:
            # if the orbiter is very close to the poles then we take latitude to be 89.99 deg 
            # to avoid division by zero error
            lat  = 90 - margin

        delta_longitude = 0.4/np.cos(np.deg2rad(lat))
        lat_min = lat - delta_latitude/2
        lat_max = lat + delta_latitude/2
        lon_min = lon - delta_longitude/2
        lon_max = lon + delta_longitude/2
        
        # Get the color corresponding to the amplitude ratio
        color = cmap(norm(amp_ratio))
        
        # Create a shaded square region (using rectangle)
        rect = patches.Rectangle((lon_min, lat_min), delta_longitude, delta_latitude, linewidth=0, edgecolor='none', facecolor=color, alpha=0.9)
        ax.add_patch(rect)
    
def overlay_on_moon_and_save(fig, element):
    # moon_image_path = 'Lunar Projections/8k_moon.jpg'  # Path to the moon image
    # moon_image_path = 'Lunar Projections/WAC_GLOBAL_E000N0000_064P.jpg'  # Path to the moon image
    moon_image_path = '/home/kaustav_b_ph.iitr/Plotting/WAC_GLOBAL_E000N0000_064P.jpg'  # Path to the moon image

    moon_image = Image.open(moon_image_path).convert('RGBA')  # Open the moon image
    
    # Save the figure to a buffer
    # fig.savefig(f'Images/{element}_temp.png', bbox_inches='tight', pad_inches=0, transparent=True, dpi=1000)
    fig.savefig(f'/home/kaustav_b_ph.iitr/Plotting/Images/{element}_temp.png', bbox_inches='tight', pad_inches=0, transparent=True, dpi=1000)

    # Load the saved figure as an image
    # fig_image = Image.open(f'Images/{element}_temp.png').convert('RGBA')
    fig_image = Image.open(f'/home/kaustav_b_ph.iitr/Plotting/Images/{element}_temp.png').convert('RGBA')
   
    # Resize the figure image to match the moon image size (if needed)
    fig_image = fig_image.resize(moon_image.size, Image.Resampling.LANCZOS)
    
    # Overlay the figure image on the moon image
    moon_image.paste(fig_image, (0, 0), fig_image)  # Paste with transparency
    
    # Convert the moon image to RGB to save it as JPEG (remove alpha channel)
    moon_image_rgb = moon_image.convert('RGB')
    
    # Save the final overlayed image as a JPEG
    moon_image_rgb.save(f'Images/{element}_overlayed.jpg', 'JPEG', dpi = (1000, 1000))
    moon_image_rgb.save(f'/home/kaustav_b_ph.iitr/Plotting/Images/{element}_overlayed.jpg', 'JPEG', dpi = (1000, 1000))

    return

# Directories containing CSV files
# csvFilesDirectories = [
#     "/home/sammy/Downloads/OutputC2019-22Complete/Output/C",
#     "/home/sammy/Downloads/OutputC2023-24Complete/Output/C",
#     "/home/sammy/Downloads/OutputM2019-22Complete/Output/M",
#     "/home/sammy/Downloads/OutputM2023-24Complete/Output/M",
#     "/home/sammy/Downloads/OutputX2023-24Complete/Output/X"
# ]
csvFilesDirectories = [
    "/home/kaustav_b_ph.iitr/Plotting/Data/OutputC2019-22Complete/Output/C",
    "/home/kaustav_b_ph.iitr/Plotting/Data/OutputC2023-24Complete/Output/C",
    "/home/kaustav_b_ph.iitr/Plotting/Data/OutputM2019-22Complete/Output/M",
    "/home/kaustav_b_ph.iitr/Plotting/Data/OutputM2023-24Complete/Output/M",
    "/home/kaustav_b_ph.iitr/Plotting/Data/OutputX2023-24Complete/Output/X"
]

elements_ = ['Fe', 'Ti', 'Ca', 'Al', 'Mg', 'Na', 'O', 'Si']
figs_and_axs = [plt.subplots(figsize=(12, 6)) for i in range(len(elements_))]

highland_subset_files = dict([(element, f"/home/kaustav_b_ph.iitr/Plotting/Subsets/{highland}_{element}_subset.csv") for highland in highland_regions.keys() for element in elements_])
mare_subset_files = dict([(element, f"/home/kaustav_b_ph.iitr/Plotting/Subsets/{mare}_{element}_subset.csv") for mare in mare_regions.keys() for element in elements_])
csv_subset_files = dict([(element, f"Subsets/{element}_subset.csv") for element in elements_])
# print(csv_subset_files)

# Collect all CSV files from the directories
csv_files = []
for directory in csvFilesDirectories:
    allFiles = os.listdir(directory)
    csv_files.extend([os.path.join(directory, f) for f in allFiles if f.endswith(".csv")])

# Process each CSV file
for csv_file in csv_files:
    try:
        print(csv_file)
        parse_and_generate_subsets(figs_and_axs, csv_file, csv_subset_files)
    except Exception as e:
        print(f"Useless File: {csv_file}, Error: {e}")



# Generate and save plots for each element
for canvas, element in zip(figs_and_axs, elements_):
    fig, ax = canvas

    # Add grid lines at 15-degree intervals
    ax.grid(True, which='both', axis='both', linestyle=':', linewidth=0.5)
    ax.set_xticks(range(-180, 185, 15))
    ax.set_yticks(range(-90, 95, 15))
    
    # Set title
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')
    fig.savefig(f'/home/kaustav_b_ph.iitr/Plotting/Images/{element}.png', bbox_inches='tight', pad_inches=0, transparent=True)
    overlay_on_moon_and_save(fig, element)