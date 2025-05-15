# ############ For mare detection ############
# """
# Define the function to return true if the given point is inside the mare region
# """

# def mare1(lat, lon):
#     # Mare 1
#     mareName = "Mare 1"
#     if lat >= -20 and lat <= 20 and lon >= -20 and lon <= 20:
#         return True
#     return False, mareName

# def mare2(lat, lon):
#     # Mare 2
#     mareName = "Mare 2"
#     if lat >= 20 and lat <= 40 and lon >= -20 and lon <= 20:
#         return True
#     return False, mareName

# def mare3(lat, lon):
#     # Mare 3
#     mareName = "Mare 3"
#     if lat >= -20 and lat <= 20 and lon >= 20 and lon <= 40:
#         return True
#     return False, mareName

# def in_mare(lat, lon):
#     for mare in mare_regions:
#         inMare, mareName = mare(lat, lon)
#         if mare(lat, lon):
#             return True, mareName
#     return False, None

# # key values pairs with mare names as keys and mare functions as values
# mare_regions = {
#     "Mare 1": mare1,
#     "Mare 2": mare2,
#     "Mare 3": mare3
# }

import geopandas as gpd
from shapely.geometry import Point
import os
import pdb 
os.environ["PROJ_IGNORE_CELESTIAL_BODY"] = "YES"

def load_mare_shapefile(shapefile_path):
    """
    Load the shapefile containing the lunar mare boundaries.
    """
    # Read the shapefile into a GeoDataFrame
    gdf = gpd.read_file(shapefile_path)
    return gdf

def get_mare_name(lat, lon, gdf):
    """
    Given a latitude and longitude, return the name of the lunar mare
    if the point lies within the boundaries of any mare.
    """
    # Create a Point object from the input coordinates
    point = Point(lon, lat)
    
    # Iterate over each geometry (mare) in the GeoDataFrame
    for _, row in gdf.iterrows():
        # print(row['geometry'])
        if row['geometry'].contains(point):  # Check if point is inside the polygon (mare)
            return row['Name']  # Adjust field name if necessary
    return "No mare found at these coordinates"

# Load the shapefile
# shapefile_path = '/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/Mare/LROC_GLOBAL_MARE_180/LROC_GLOBAL_MARE_180.SHP'  # Replace with actual path to shapefile
shapefile_path = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/Mare/LROC_GLOBAL_MARE_360/LROC_GLOBAL_MARE_360.SHP"
mare_gdf = load_mare_shapefile(shapefile_path)

# Print the CRS of the shapefile
print(f"Shapefile CRS: {mare_gdf.crs}")

mare_gdf.plot()

# # Reproject to EPSG:4326 (WGS84 lat/lon)
# if mare_gdf.crs != "EPSG:4326":
#     mare_gdf = mare_gdf.to_crs(epsg=4326)
#     print("Reprojected CRS to EPSG:4326")

# Example usage:
lat = float(input("Enter the latitude: "))
lon = float(input("Enter the longitude: "))
mare_name = get_mare_name(lat, lon, mare_gdf)
print(f"Point lies in: {mare_name}")

for lat in range(-90, 91, 10):
    for lon in range(-180, 181, 10):
        try:
            mare_name = get_mare_name(lat, lon, mare_gdf)
            print(f"Point ({lat}, {lon}) lies in: {mare_name}")
        except:
            pdb.set_trace()
