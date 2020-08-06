# Rose diagrams
# Script to make rose diagram of the basin orientation and fault strike.
# FJC 30/06/20

# import packages
import numpy as np
import matplotlib.pyplot as plt

# shapefiles
import geopandas as gpd

# read in the basin shapefiles
data_dir = '/home/bjdd72/san_andreas/'
gdf = gpd.read_file(data_dir+'SAF_only_basins.shp')

print(gdf)

# bin the orientations into 10 degree bins.
bin_width = 5
bin_edges = np.arange(0, 360+bin_width, bin_width)
basins_hist, bin_edges = np.histogram(gdf['azimuth'], bins=bin_edges)

# read in the fault azimuth shapefile
pts = gpd.read_file(data_dir+'SanAndreasPoints.shp')
# double the strikes with directions at 180 degrees to get a mirrored rose diagram
strikes = pts['azimuth'].values
strikes_mirrored = strikes + 180
all_strikes = np.concatenate([strikes, strikes_mirrored])
strike_hist, bin_edges = np.histogram(all_strikes, bins=bin_edges)

# set up the figure
fig = plt.figure(figsize=(14,8))

# set up the fault strike subplot
ax = fig.add_subplot(121, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title('Fault strike orientation', y=1.1, fontsize=16)
ax.bar(np.deg2rad(bin_edges[0:-1]), strike_hist, width=np.deg2rad(bin_width), align='edge', edgecolor='0.2', color='0.5', zorder=2, alpha=0.5)
ax.tick_params(labelsize=12)

# set up the basin orientation subplot
ax = fig.add_subplot(122, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title('Drainage basin orientation', y=1.1, fontsize=16)
ax.tick_params(labelsize=12)

# plot the histogram as bars. The bin edges are used as the left hand coordinate for the bars, minus the last edge. The histogram counts are the bar heights.
ax.bar(np.deg2rad(bin_edges[0:-1]), basins_hist, width=np.deg2rad(bin_width), align='edge', edgecolor='0.2', color='0.5', zorder=2, alpha=0.5)
#plt.subplots_adjust(wspace=0.5)
plt.savefig('SAF_rose_diagrams.png', dpi=300)
