# Read in the basin data and collate it into some catchment averaged statistics
# The basins themseleves are a raster - read this in and get the centroid of each
# to use as the latitude and longitude.
# The info is from the hilltop flow routing algorithms as a csv file. use pandas
# to read this in and then get statistics for each basin. Output to a csv which will
# be appended for the whole fault.

import pandas as pd
import numpy as np
import osgeo.gdal as gdal
import os
from os.path import exists
from osgeo import osr
from pyproj import Proj, transform

# shapefiles
#import rasterio
#from rasterio.features import shapes
from shapely.geometry import shape, Polygon, mapping, Point, LineString
#import fiona
import geopandas as gpd
import math

def azimuth(row):
    """
    azimuth between 2 shapely points from a geodataframe row
    """
    point1 = row['centroids']
    point2 = row['outlet_coords']
    angle = math.atan2(point2.x - point1.x, point2.y - point1.y)
    initial_az = math.degrees(angle)

    # Now we have the initial bearing but math.atan2 returns values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    final_az = (initial_az + 360) % 360
    #print(initial_az, final_az)
    return final_az

def get_basin_orientation(DataDirectory, basin_name):
    """
    This function takes the basin shapefile and gets the orientation of each basin.
    It calculates the minimum bounding rectangle (MBR) of each basin and then finds the
    longest axis, then uses the angle of this axis as the orientation

    FJC 24/06/20
    """

    # read in the outlet coords to a geodataframe
    df = pd.read_csv(DataDirectory+basin_name+'/'+basin_name+'_hillslopes_SO3.csv')
    geometry = [Point(xy) for xy in zip(df.longitude_outlet, df.latitude_outlet)]
    crs = 'epsg:4326' #http://www.spatialreference.org/ref/epsg/2263/
    gdf_outlets = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
    # convert to UTM10
    gdf_outlets = gdf_outlets.to_crs("EPSG:32610")

    # read in the shapefile to a geodataframe
    gdf = gpd.read_file(DataDirectory+'/'+'basin_shapefiles/'+basin_name+'_basins.shp')
    # get centroids of each polygon
    gdf['centroids'] = gdf['geometry'].centroid

    # join outlets gdf with basins gdf on basin ID
    gdf = gdf.merge(gdf_outlets, left_on='DN', right_on='basin_id')
    # rename the geometry columns
    gdf = gdf.rename(columns={'geometry_x':'geometry', 'geometry_y': 'outlet_coords'})
    print(gdf.columns)

    # calculate aziumth between centroid and outlet (outlet as 2nd point) to get orientation of the basin.
    gdf['azimuth'] = gdf.apply(azimuth, axis=1)
    final_gdf = gpd.GeoDataFrame(gdf[['basin_id', 'azimuth', 'latitude_outlet', 'longitude_outlet', 'slope_median', 'slope_16th', 'slope_84th']], geometry=gdf['geometry'], crs='EPSG:32610')
    # now save the gdf to a shapefile for checking
    final_gdf.to_file(DataDirectory+'/'+'basin_shapefiles/'+basin_name+'_basins_az.shp')


DataDirectory = '/raid/fclubb/san_andreas/SouthernSAF/'
filenames = [f.name for f in os.scandir(DataDirectory) if f.is_dir()]
for f in filenames:
    if 'SouthernSAF_' in f:
        print(f)
        get_basin_orientation(DataDirectory, f)
