# Rose diagrams
# Script to make rose diagram of the basin orientation and fault strike.
# FJC 30/06/20

# import packages
import numpy as np
import matplotlib.pyplot as plt

# shapefiles
import geopandas as gpd

gdf = gpd.read_file('/home/bjdd72/san_andreas/SAF_only_basins.shp')

print(gdf)

# bin the orientations into 10 degree bins.
bin_width = 10
bin_edges = np.arange(0, 370, 10)
hist, bin_edges = np.histogram(gdf['azimuth'], bins=bin_edges)

# set up the figure
fig = plt.figure(figsize=(10,8))

# set up the fault strike subplot

# set up the basin orientation subplot
ax = fig.add_subplot(111, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title('Drainage basin orientation')

# plot the histogram as bars. The bin edges are used as the left hand coordinate for the bars, minus the last edge. The histogram counts are the bar heights.
ax.bar(np.deg2rad(bin_edges[0:-1]), hist, width=np.deg2rad(10), align='edge', edgecolor='k', color='deepskyblue', zorder=2, alpha=0.5)
plt.savefig('SAF_only_basin_orientation.png', dpi=300)
