from mareDetection import MareGrid, mareList
from mareDetection import *
import pandas as pd
import glob
import os
import numpy as np

# Initialize MareGrid
mare_grid = MareGrid('marehpc.txt', 'marename_hpc.txt')

# Get list of all CSV files in Subsets directory
csv_files = glob.glob('Subsets/*_subset.csv')

# Dictionary to store results for each mare and element
mare_element_data = {mare: {} for mare in mareList}

# Process each CSV file
# ... existing code ...

# Try different encodings when reading CSV files
def read_csv_with_encoding(file_path, usecols):
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    for encoding in encodings:
        try:
            return pd.read_csv(file_path, usecols=usecols, encoding=encoding)
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError(f"Unable to read {file_path} with any of the attempted encodings")

# Update the CSV reading part in your functions

# Process each CSV file
for csv_file in csv_files:
    element_name = os.path.basename(csv_file).replace('_subset.csv', '')
    
    # Read CSV with required columns
    df = read_csv_with_encoding(csv_file, usecols=['BORE_LAT', 'BORE_LON', f'{element_name}_amplitude_ratio_vs_Si'])

    
    # For each row, determine which mare it belongs to
    for index, row in df.iterrows():
        cell = mare_grid.mareAtLocation((row['BORE_LAT'], row['BORE_LON']))
        mare_name = cell.mareName if cell.inMare else None
        
        if mare_name:  # If point is in a mare
            # Initialize list for this element in this mare if it doesn't exist
            if element_name not in mare_element_data[mare_name]:
                mare_element_data[mare_name][element_name] = []
            
            # Append the amplitude ratio
            mare_element_data[mare_name][element_name].append(
                row[f'{element_name}_amplitude_ratio_vs_Si']
            )

            
# Save results to CSV files
for mare_name in mareList:
    if mare_element_data[mare_name]:  # If mare has any data
        # Convert to DataFrame
        mare_df = pd.DataFrame(
            {elem: pd.Series(values) 
             for elem, values in mare_element_data[mare_name].items()}
        )
        
        # Save to CSV
        mare_df.to_csv(f'mare_data/{mare_name.replace(" ", "_")}_elements.csv', index=False)

# Calculate and print statistics for each mare and element
for mare_name in mareList:
    if mare_element_data[mare_name]:  # If mare has any data
        print(f"\n{mare_name}:")
        for element, values in mare_element_data[mare_name].items():
            values_array = np.array(values)
            mean = np.mean(values_array)
            std = np.std(values_array)
            print(f"  {element:15} - Mean: {mean:.3f}, Std: {std:.3f}, Count: {len(values)}")

from mareDetection import MareGrid, mareList
from mareDetection import *
import pandas as pd
import glob
import os
import numpy as np

def calculate_region_averages(min_lat, max_lat, min_lon, max_lon):
    """
    Calculate element averages for a specified latitude and longitude range
    
    Parameters:
    min_lat, max_lat: latitude range
    min_lon, max_lon: longitude range
    """
    # Get list of all CSV files in Subsets directory
    csv_files = glob.glob('Subsets/*_subset.csv')
    
    # Dictionary to store results for each element
    region_element_data = {}
    
    # Process each CSV file
    for csv_file in csv_files:
        element_name = os.path.basename(csv_file).replace('_subset.csv', '')
        
        # Read CSV with required columns
        df = pd.read_csv(csv_file, usecols=['BORE_LAT', 'BORE_LON', f'{element_name}_amplitude_ratio_vs_Si'])
        
        # Filter data for specified region
        region_df = df[
            (df['BORE_LAT'] >= min_lat) & 
            (df['BORE_LAT'] <= max_lat) & 
            (df['BORE_LON'] >= min_lon) & 
            (df['BORE_LON'] <= max_lon)
        ]
        
        # Store the amplitude ratios for this element
        region_element_data[element_name] = region_df[f'{element_name}_amplitude_ratio_vs_Si'].tolist()
    
    # Calculate and print statistics for each element
    print(f"\nStatistics for region (Lat: {min_lat} to {max_lat}, Lon: {min_lon} to {max_lon}):")
    for element, values in region_element_data.items():
        if values:  # If there's data for this element
            values_array = np.array(values)
            mean = np.mean(values_array)
            std = np.std(values_array)
            print(f"  {element:15} - Mean: {mean:.3f}, Std: {std:.3f}, Count: {len(values)}")
        else:
            print(f"  {element:15} - No data points found in specified region")
    
    return region_element_data

from mareDetection import MareGrid, mareList
from mareDetection import *
import pandas as pd
import glob
import os
import numpy as np

def get_antipodal_point(lat, lon):
    """
    Calculate the antipodal point for given latitude and longitude
    
    Parameters:
    lat: latitude of the original point
    lon: longitude of the original point
    
    Returns:
    tuple: (antipodal_lat, antipodal_lon)
    """
    antipodal_lat = -lat
    antipodal_lon = lon + 180 if lon < 0 else lon - 180
    return antipodal_lat, antipodal_lon

def calculate_antipodal_region_stats(center_lat, center_lon, region_size=30):
    """
    Calculate element statistics for a region around the antipodal point
    
    Parameters:
    center_lat: latitude of the original impact point
    center_lon: longitude of the original impact point
    region_size: size of the region to analyze in degrees (default 30)
    """
    # Calculate antipodal point
    anti_lat, anti_lon = get_antipodal_point(center_lat, center_lon)
    
    # Calculate region boundaries
    half_size = region_size / 2
    min_lat = anti_lat - half_size
    max_lat = anti_lat + half_size
    min_lon = anti_lon - half_size
    max_lon = anti_lon + half_size
    
    # Get list of all CSV files in Subsets directory
    csv_files = glob.glob('Subsets/*_subset.csv')
    
    # Dictionary to store results for each element
    region_element_data = {}
    
    print(f"\nAnalyzing antipodal region for impact at ({center_lat}°, {center_lon}°)")
    print(f"Antipodal point: ({anti_lat}°, {anti_lon}°)")
    print(f"Region size: {region_size}° x {region_size}°")
    
    # Process each CSV file
    for csv_file in csv_files:
        element_name = os.path.basename(csv_file).replace('_subset.csv', '')
        
        # Read CSV with required columns
        df = pd.read_csv(csv_file, usecols=['BORE_LAT', 'BORE_LON', f'{element_name}_amplitude_ratio_vs_Si'])
        
        # Filter data for specified region
        region_df = df[
            (df['BORE_LAT'] >= min_lat) & 
            (df['BORE_LAT'] <= max_lat) & 
            (df['BORE_LON'] >= min_lon) & 
            (df['BORE_LON'] <= max_lon)
        ]
        
        # Store the amplitude ratios for this element
        region_element_data[element_name] = region_df[f'{element_name}_amplitude_ratio_vs_Si'].tolist()
    
    # Calculate and print statistics for each element
    print(f"\nStatistics for antipodal region:")
    for element, values in region_element_data.items():
        if values:  # If there's data for this element
            values_array = np.array(values)
            mean = np.mean(values_array)
            std = np.std(values_array)
            print(f"  {element:15} - Mean: {mean:.3f}, Std: {std:.3f}, Count: {len(values)}")
        else:
            print(f"  {element:15} - No data points found in specified region")
    
    return region_element_data

# Analyze antipodal region for a specific point
data = calculate_antipodal_region_stats(0,0)  # Imbrium Basin example