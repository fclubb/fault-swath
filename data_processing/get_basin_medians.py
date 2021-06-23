# get_basin_medians.py
# This script creates a shapefile of the basins along the SAF and gets the median channel gradient,
# hillslope gradient and hilltop curvature in each basin
# FJC 14/06/21

# import modules
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Point

def percentile(n):
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_


# read in the river profile CSV
data_dir = '/raid/fclubb/san_andreas/USGS_NED_10m/SAF_combined/'
fname = 'SAF_combined_10m'
df = pd.read_csv(data_dir+fname+'_profiles_SO3.csv')
df = df[df['slope'] > 0]
df.columns

# read in the hillslope metrics CSV
hs_df = pd.read_csv(data_dir+fname+'_hillslopes_SO3.csv')

# read in the hilltop metrics CSV
ht_df = pd.read_csv(data_dir+fname+'_RidgeData_SO3.csv')


# convert the river csv to a geodataframe. Remove the non-unique ID labels - these will be replaced by unique basin IDs
geometry = [Point(xy) for xy in zip(df.longitude, df.latitude)]
crs = 'epsg:4326' #http://www.spatialreference.org/ref/epsg/2263/
river_gdf = gpd.GeoDataFrame(df.drop(['latitude','longitude','basin_id','id','new_id','node'], axis=1), crs=crs, geometry=geometry)
river_gdf_clean = river_gdf[river_gdf.geometry.type == 'Point']

# convert the hillslope csv to a geodataframe. Remove the non-unique ID labels
geometry = [Point(xy) for xy in zip(hs_df.longitude_outlet, hs_df.latitude_outlet)]
hs_gdf = gpd.GeoDataFrame(hs_df, crs=crs, geometry=geometry)

# convert the hilltop csv to a geodataframe. Remove the non-unique ID labels
geometry = [Point(xy) for xy in zip(ht_df.longitude, ht_df.latitude)]
ht_gdf = gpd.GeoDataFrame(ht_df.drop(['latitude','longitude','basin_id','new_id'], axis=1), crs=crs, geometry=geometry)

# add a unique id to the basin
basin_gdf = gpd.read_file(data_dir+fname+'_basins_SO3.shp')
# convert the basin GDF to WGS84
basin_gdf = basin_gdf.to_crs('epsg:4326')
#basin_gdf = basin_gdf.drop(['basin_id'], axis=1)
basin_gdf['unique_id'] = basin_gdf.index
basin_gdf = basin_gdf[basin_gdf.geometry.type == 'Polygon']

# merge the river and basins gdf and calculate the median channel slope in each basin
join = gpd.sjoin(river_gdf, basin_gdf, how='left', op='intersects')
gr = join.groupby(['unique_id'])['slope'].agg(['median', 'std', percentile(16), percentile(84)]).rename(columns={'median': 'channel_slope_median', 'std': 'channel_slope_std', 'percentile_16': 'channel_slope_16th', 'percentile_84': 'channel_slope_84th'}).reset_index()
basin_gdf = basin_gdf.merge(gr, on='unique_id')

# now join the hillslope data
join = gpd.sjoin(basin_gdf, hs_gdf, how='left', op='contains')

# now join the hilltop data - find points within the basin and get the median curvature in each basin
join = join.drop(['index_right'], axis=1)
ht_join = gpd.sjoin(ht_gdf, join, how='left', op='within')
gr = ht_join.groupby(['unique_id'])['curvature'].agg(['median', 'std', percentile(16), percentile(84)]).rename(columns={'median': 'ht_curv_median', 'std': 'ht_curv_std', 'percentile_16': 'ht_curv_16th', 'percentile_84': 'ht_curv_84th'}).reset_index()
join = join.merge(gr, on='unique_id')
print(len(join.unique_id.unique()))

# write to shapefile
join.to_file(data_dir+fname+'_channels_plus_hilltops_by_basin_SO3.shp')




