import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
import numpy as np
from matplotlib.colors import Normalize

def parse_and_generate_subsets(csv_file):
    # Read the CSV data
    df = pd.read_csv(csv_file)
    
    # Define the columns related to each element
    elements = ['Fe', 'Ti', 'Ca', 'Al', 'Mg', 'Na', 'O', 'Si']
    
    # Create an empty dictionary to store subsets for each element
    subsets = {element: [] for element in elements}
    
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

    # Convert subsets from lists of rows to DataFrames for easy inspection
    subset_dfs = {element: pd.DataFrame(subsets[element]) for element in elements}
    
    # Save each subset to a separate CSV file
    for element, subset_df in subset_dfs.items():
        if not subset_df.empty:  # Check if the subset is not empty
            subset_df.to_csv(f'{element}_subset.csv', index=False)
            # Plot the shaded region for the non-empty subset
            plot_shaded_region_on_map(subset_df, element)
        else:
            print(f"Subset for {element} is empty. Skipping plot.")
        
    print(f"Subsets generated and saved for each element: {', '.join(elements)}")


def plot_shaded_region_on_map(subset_df, element):
    # Extract the relevant data
    latitudes = subset_df['BORE_LAT']
    longitudes = subset_df['BORE_LON']
    amplitude_ratio = subset_df[f'{element}_amplitude_ratio_vs_Si']
    
    # Set up the plot
    fig, ax = plt.subplots(figsize=(12, 6))  # Explicitly create axes
    
    # Normalize the amplitude ratios for coloring
    # norm = Normalize(vmin=np.min(amplitude_ratio), vmax=np.max(amplitude_ratio))
    norm = Normalize(vmin=0, vmax=2)
    cmap = plt.cm.viridis
    
    # Create ScalarMappable for colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # To ensure it works with colorbar
    
    # Iterate over each point and plot the shaded region (0.4°x0.4°)
    for lat, lon, amp_ratio in zip(latitudes, longitudes, amplitude_ratio):

        # norm = Normalize(vmin=np.min(amp_ratio), vmax=np.max(amplitude_ratio))
        # Calculate the boundary of the square region
        lat_min = lat - 0.2
        lat_max = lat + 0.2
        lon_min = lon - 0.2
        lon_max = lon + 0.2
        
        # Get the color corresponding to the amplitude ratio
        color = cmap(norm(amp_ratio))

        # Create a shaded square region (using rectangle)
        rect = patches.Rectangle((lon_min, lat_min), 0.4, 0.4, linewidth=0, edgecolor='none', facecolor=color, alpha=0.7)
        ax.add_patch(rect)
    
    # Add colorbar for the heatmap
    cbar = plt.colorbar(sm, ax=ax, label=f'{element} Amplitude Ratio vs Si')
    
    # Add grid lines at 15-degree intervals
    ax.grid(True, which='both', axis='both', linestyle=':', linewidth=0.5)
    ax.set_xticks(range(-180, 185, 15))
    ax.set_yticks(range(-90, 95, 15))
    
    # Set title
    ax.set_title(f'{element} Region Shading with Heatmap (0.4° x 0.4°)')
    ax.set_aspect('equal', adjustable='box')
    # Show the plot
    plt.show()
    
# Example usage
csv_file = '/home/sammy/Downloads/output_64.csv'
parse_and_generate_subsets(csv_file)
