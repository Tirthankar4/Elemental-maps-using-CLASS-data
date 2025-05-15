import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
import seaborn as sns
import concurrent.futures
import numpy as np
import os
import csv
import pdb 
# os.environ["PROJ_IGNORE_CELESTIAL_BODY"] = "YES"


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
        if row['geometry'].contains(point):  # Check if point is inside the polygon (mare)
            return row['MARE_NAME']  # Adjust field name if necessary
    return None


def showMarePlot(ax, mare_gdf):
    # Annotate each mare with its name
    for _, row in mare_gdf.iterrows():
        # Get the geometry (polygon) and the name of the mare
        geometry = row['geometry']
        mare_name = row['MARE_NAME']

        # Get the coordinates of the centroid of the polygon for annotation
        centroid = geometry.centroid
        x, y = centroid.x, centroid.y

        # Annotate the plot with the mare name at the centroid
        ax.text(x, y, mare_name, fontsize=8, ha='center', color='black')

    # Set the title and labels for the plot
    ax.set_title("Lunar Mare Polygon")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    # Show the plot
    plt.show()


def sortByArea(mare_gdf):
    
    # Calculate the area of each mare (ensure CRS is in a projected coordinate system)
    mare_gdf['area'] = mare_gdf['geometry'].area

    # Sort the GeoDataFrame by area in descending order
    gdf_sorted = mare_gdf.sort_values(by='area', ascending=False)

    # Get the top 7 largest mares
    top_20_mares = gdf_sorted.head(20)

    # Print the names and areas of the top 7 largest mares
    for idx, row in top_20_mares.iterrows():
        mare_name = row['MARE_NAME']
        mare_area = row['area']
        print(f"{mare_name}")
        # print(f"Mare: {mare_name}, Area: {mare_area:.2f} square degrees")
with open('mare_boolean.txt', mode='a', newline='') as file:
                writer = csv.writer(file)

# Function to check if a point lies inside any of the major mares
def check_point_in_mares(point, gdf_major_mares):
    lat, lon = point
    # Create a Point object for the current lat/lon
    p = Point(lon, lat)
    
    # Check if the point is inside any of the major mares
    for _, row in gdf_major_mares.iterrows():
        if row['geometry'].contains(p):
            return [True, row['MARE_NAME']]
    return [False, None]


# Function to process a chunk of grid points
def process_grid_chunk(chunk, gdf_major_mares):
    resultList = [check_point_in_mares(point, gdf_major_mares) for point in chunk]
    bools = [row[0] for row in resultList]
    names = [row[1] for row in resultList]
    return list(bools), list(names)


# Parallelize the task using ProcessPoolExecutor
def parallelize_grid_check(grid_points, gdf_major_mares, num_workers=4):
    # Split the grid points into chunks (for 4 cores)
    chunk_size = len(grid_points) // num_workers
    chunks = [grid_points[i:i + chunk_size] for i in range(0, len(grid_points), chunk_size)]

    # Use concurrent.futures to process the chunks in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        _ = list(executor.map(process_grid_chunk, chunks, [gdf_major_mares]*len(chunks)))
        # pdb.set_trace()

        results = [row[0] for row in _]
        nameList = [row[1] for row in _]

    # Flatten the results from the chunks
    results_, nameList_ = [], []
    results_ = [item for sublist in results for item in sublist]
    nameList_ = [item for sublist in nameList for item in sublist]
    return results_, nameList_


class MareGrid(object):

    class GridCell(object):
        def __init__(self, gridObject, location, mare=None, insideMare=False):
            self.gridObject = gridObject
            self.gridLocation = location
            self.geometricLocation = self.gridObject.__gridLocationToCoordinates__(location)
            self.mareName = mare
            self.inMare = insideMare



    def __init__(
        self, resultFile, nameFile, mareList=[
        "Oceanus Procellarum", "Mare Imbrium", "Mare Frigoris", "Mare Insularum", "Mare Tranquillitatis", 
        "Mare Serenitatis", "Mare Fecunditatis", "Mare Nubium", "Mare Crisium", 
        "Mare Humorum", "Mare Australe", "Lacus Somniorium", "Mare Cognitum", 
        "Mare Nectaris", "Cognitum East", "Mare Vaporum", "Lacus Mortis", 
        "Mare Orientale", "Mare Marginis", "Mare Ingenii"]
        ) -> None: 
        #, "Mare Crisium", "Mare Humorum", "Mare Australe", "Lacus Somniorium", "Mare Cognitum", 
        #"Mare Nectaris", "Cognitum East", "Mare Vaporum", "Lacus Mortis", 
        #"Mare Orientale", "Mare Marginis", "Mare Ingenii"
        # print(mare_gdf.shape)
        # pdb.set_trace()
        # self.gdf = mare_gdf[mare_gdf['MARE_NAME'].isin(mareList)]
        # print(self.gdf.shape)
        self.mareList = mareList
        self.__loadGrid__(resultFile, nameFile)
        # self.__setCellSize__()
        return 

    def __setCellSize__(self):
        # latitudes: -90 to 90 in degrees 
        # longitudes: -180 to 180 in degrees 
        self.latitudeResolution = 180/self.size[0]
        self.longitudeResolution = 360/self.size[1]

    def __generateGrid__(self, results, names):
        nRows, nCols = self.size
        self.grid = []

        for i in range(nRows):
            row = []
            for j in range(nCols):
                cell = self.GridCell(self, (i, j))
                # mare = None
                print(i, j)
                mare = names[i][j]
                if mare in self.mareList:
                    print(mare)
                    cell.mareName = mare
                    cell.inMare = True
                else:
                    pass
                row.append(cell)
            self.grid.append(row)

        return None

    def __coordinatesToGridLocation__(self, coordinates):
        lat, lon = coordinates
        return (int((lat+90)/self.latitudeResolution), int((lon+180)/self.longitudeResolution))

    def __gridLocationToCoordinates__(self, location):
        i, j = location
        return (-90+self.latitudeResolution*(i+0.5), -180+self.longitudeResolution*(j+0.5))

    def mareAtLocation(self, coordinates):
        i, j = self.__coordinatesToGridLocation__(coordinates)
        return self.grid[i][j]

    def __loadGrid__(self, resultFile, nameFile):
        # pdb.set_trace()
        results = []
        names = []
        with open(resultFile, mode='r', newline='') as file:
            reader = csv.reader(file)
            
            # Read the first row and convert it into a list
            results = next(reader)

        with open(nameFile, mode='r', newline='') as file:
            reader = csv.reader(file)
            
            # Read the first row and convert it into a list
            names = next(reader)
        results_ = np.array(results).reshape((nLats, nLongs))
        names_ = np.array(names).reshape((nLats, nLongs))
        # tempr_ = list(results)
        # tempn_= list(names)
        # names_ = []
        # for i in tempn_:
        #     if len(i) == 0:
        #         names_.append(None)
        #     else:
        #         names.append(i[0])
        # results = []
        # results = [tempr_[i:(i + nLats)] for i in range(0, nLats*nLongs, nLats)]
        # names = [tempn[i:i + nLats] for i in range(0, len(tempn_), nLats)]

        self.size = [nLats, nLongs]
        self.__setCellSize__()
        self.__generateGrid__(results_, names_)
        return

    # def __saveLoad__(self, resultList, nameList):
    #     # boolean
    #     try:
    #         with open('mare_boolean.txt', mode='a', newline='') as file:
    #             writer = csv.writer(file)
            
    #         # Write the list as a single row in the CSV file
    #         writer.writerow(resultList)
    #     except:
    #         with open('mare_boolean.txt', mode='a', newline='') as file:
    #             writer = csv.writer(file)
            
    #         # Write the list as a single row in the CSV file
    #         writer.writerow(resultList)

    #     # names
    #     try:
    #         with open('mare_names.txt', mode='a', newline='') as file:
    #             writer = csv.writer(file)
            
    #         # Write the list as a single row in the CSV file
    #         writer.writerow(resultList)
    #     except:
    #         with open('mare_names.txt', mode='a', newline='') as file:
    #             writer = csv.writer(file)
            
    #         # Write the list as a single row in the CSV file
    #         writer.writerow(nameList)
    #     return


def generateMareGrid(size, mare_gdf):
    return 

def heatmap(boolean_grid):

    # Create a 2D boolean grid (for example, a 10x10 grid)
    # boolean_grid = np.random.choice([True, False], size=(10, 10))

    # Create a heatmap using seaborn
    plt.figure(figsize=(8, 6))
    sns.heatmap(boolean_grid, cmap='Blues', cbar=False, square=True, annot=True, fmt="d", linewidths=0.5)

    # Set title and labels
    plt.title("Boolean Grid Heatmap", fontsize=16)
    plt.xlabel("X-axis", fontsize=12)
    plt.ylabel("Y-axis", fontsize=12)

    # Show the plot
    plt.show()


def main():

    # # Load the shapefile
    shapefile_path = "180/LROC_GLOBAL_MARE_180.SHP"
    mare_gdf = load_mare_shapefile(shapefile_path)

    axHandle = mare_gdf.plot(color='lightblue', edgecolor='black')

    sortByArea(mare_gdf=mare_gdf)
    # heatmap([[cell.inMare for cell in row] for row in mygrid.grid])
    # pdb.set_trace()

    # # List of the names of the 7 major mares
    # # Assume that 'Name' column has the mare names
    # major_mares = ['Mare Imbrium', 'Mare Serenitatis', 'Mare Crisium', 
    #                'Mare Tranquillitatis', 'Mare Nubium', 'Mare Humorum', 'Mare Fecunditatis']
    # major_mare_grid = mare_gdf[mare_gdf['MARE_NAME'].isin(major_mares)]

    # # Filter the GeoDataFrame for the 7 major mares
    # # gdf_major_mares = mare_gdf[mare_gdf['MARE_NAME'].isin(major_mares)]

    # # Create a (N x M) grid with random latitude and longitude (you can adjust the grid size as needed)
    # # Example grid size
    # longitude = np.linspace(-180, 180, nLongs)  # Longitude range from -180 to 180
    # latitude = np.linspace(-90, 90, nLats)    # Latitude range from -90 to 90

    # # Create a meshgrid for the N x M grid points (longitude, latitude)
    # lon_grid, lat_grid = np.meshgrid(longitude, latitude)

    # # Flatten the meshgrid to make it easier to parallelize
    # grid_points = np.vstack([lat_grid.ravel(), lon_grid.ravel()]).T

    # results, names = parallelize_grid_check(grid_points, major_mare_grid)
    # results = [int(value) for value in results]

    # # [results[i:i + M] nfor i in range(0, len(results), M)]
    # # pdb.set_trace()

    # with open('mare_boolean.txt', mode='w', newline='') as file:
    #     writer = csv.writer(file)
    
    #     # Write the list as a single row in the CSV file
    #     writer.writerow(results)


    # with open('mare_names.txt', mode='w', newline='') as file:
    #     writer = csv.writer(file)
    
    #     # Write the list as a single row in the CSV file
    #     writer.writerow(names)


    # print("Files saved!")


    # print("Trials.....")
    # (nLats, nLongs), 
    # mygrid = MareGrid('mare_boolean.txt', 'mare_names.txt')
    mygrid = MareGrid('marehpc.txt', 'marename_hpc.txt')

    # pdb.set_trace()
    # Reshape the results back into the shape of the grid (N, M)
    # heatmap_grid = np.array(results).reshape(nLats, nLongs)[::-1]

    # Plotting the heatmap (optional)
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.imshow(heatmap_grid, cmap='Blues', interpolation='nearest')
    # ax.colorbar(label='Inside Mare (True/False)')
    ax.set_title("Grid Points Inside Major Mares")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    plt.show()
    plt.figure()

    plt.plot()

global nLats, nLongs, mareList
nLats, nLongs = 1000, 2000
mareList=[
        "Oceanus Procellarum", "Mare Imbrium", "Mare Frigoris", "Mare Insularum", "Mare Tranquillitatis", 
        "Mare Serenitatis", "Mare Fecunditatis", "Mare Nubium", "Mare Crisium", 
        "Mare Humorum", "Mare Australe", "Lacus Somniorium", "Mare Cognitum", 
        "Mare Nectaris", "Cognitum East", "Mare Vaporum", "Lacus Mortis", 
        "Mare Orientale", "Mare Marginis", "Mare Ingenii"]

if __name__ == "__main__":
    main()