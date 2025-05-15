import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as mcolors

# Sample data: assuming you have a CSV or DataFrame with columns: latitude, longitude, and value
# Here we create an example DataFrame for illustration
# You would replace this with your actual grid data.

# Example grid data (replace with your actual data)
# data = {
#     'latitude': np.random.uniform(-90, 90, 100),
#     'longitude': np.random.uniform(-180, 180, 100),
#     'value': np.random.uniform(0, 100, 100)  # Random values as placeholders
# }

dir_path = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios"
elementList="Fe Ti Ca Si Al Mg Na O".split(" ")
elementName = elementList[4]
flare_num = 68
file = f"{dir_path}/results/{elementName}_{flare_num}.csv"

df = pd.read_csv(file)

# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 5))

# Set up the Basemap (equirectangular projection)
m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
            llcrnrlon=-180, urcrnrlon=180, resolution='i', ax=ax)

# Draw coastlines (or any other features you want)
m.drawcoastlines()

# Set up a color map (you can adjust this according to your needs)
cmap = plt.cm.viridis  # Change the colormap if needed
norm = mcolors.Normalize(vmin=df['value'].min(), vmax=df['value'].max())

# Plot each grid point on the map
x, y = m(df['longitude'].values, df['latitude'].values)
sc = m.scatter(x, y, c=df['value'].values, cmap=cmap, norm=norm, edgecolors='k', marker='s', s=50)

# Add a colorbar
plt.colorbar(sc, label='Value')

# Title and labels
plt.title("Moon Surface Grid with Value Color Mapping")

# Show the plot
plt.show()
