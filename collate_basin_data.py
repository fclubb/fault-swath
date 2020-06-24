# Read in the basin data and collate it into some catchment averaged statistics
# The basins themseleves are a raster - read this in and get the centroid of each
# to use as the latitude and longitude.
# The info is from the hilltop flow routing algorithms as a csv file. use pandas
# to read this in and then get statistics for each basin. Output to a csv which will
# be appended for the whole fault.

import pandas as pd
import numpy as np
import osgeo.gdal as gdal
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

def azimuth(point1, point2):
    """
    azimuth between 2 shapely points.
    """
    angle = math.atan2(point2.x - point1.x, point2.y - point1.y)
    initial_az = math.degrees(angle)

    # Now we have the initial bearing but math.atan2 returns values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    final_az = (initial_az + 360) % 360
    print(initial_az, final_az)
    return final_az

def get_basin_orientation(basin_shp):
    """
    This function takes the basin shapefile and gets the orientation of each basin.
    It calculates the minimum bounding rectangle (MBR) of each basin and then finds the
    longest axis, then uses the angle of this axis as the orientation

    FJC 24/06/20
    """

    # read in the outlet coords
    outlets = pd.read_csv('/home/bjdd72/TopographicData/TopographicData/san_andreas/NorthernSAF/tile_3/tile_3_hillslopes_SO3.csv')
    points = outlets['latitude_outlet']

    print(outlets)



    # read in the shapefile to a geopandas dataframe
    gdf = gpd.read_file(basin_shp)
    # get centroids of each polygon
    gdf['centroids'] = gdf['geometry'].centroid

    # TODO
    # read outlets csv to a gdf
    # join outlets gdf with basins gdf on basin ID
    # calculate aziumth between centroid and outlet (outlet as 2nd point) to get orientation of the basin.
    print(gdf)
    # for rec in c:
    #     polygon = Polygon(shape(gdf['geometry']))
    #     # get the minimum bounding rectangle and zip coordinates into a list of point-tuples
    #     mbr_points = list(zip(*polygon.minimum_rotated_rectangle.exterior.coords.xy))
    #     print(polygon.minimum_rotated_rectangle.exterior.coords.xy)
    #
    #     # calculate the length of each side of the minimum bounding rectangle
    #     axes = [LineString((mbr_points[i], mbr_points[i+1])) for i in range(len(mbr_points) -1)]
    #     mbr_lengths = [ax.length for ax in axes]
    #     #print(mbr_lengths)
    #
    #     # get major/minor axis measurements
    #     minor_axis = min(mbr_lengths)
    #     major_axis = max(mbr_lengths)
    #     # find the index of the major axis and get the shapely line of this axis
    #     major_idx = mbr_lengths.index(max(mbr_lengths))
    #     major_ax = axes[major_idx]
    #     #print(major_ax.coords[0])
    #
    #     # now calculate the bearing from the first and last point of the line
    #     az = azimuth(Point(major_ax.coords[0]), Point(major_ax.coords[-1]))
    #     rec['azimuth'] = az
        #print(az)
    #gdf.to_file(DataDirectory+'tile_3_basins_az.shp')


DataDirectory = '/home/bjdd72/TopographicData/TopographicData/san_andreas/NorthernSAF/basin_shapefiles/'
basin_shp = 'tile_3_basins.shp'
get_basin_orientation(DataDirectory+basin_shp)
