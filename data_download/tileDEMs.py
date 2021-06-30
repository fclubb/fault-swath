# tileDEMs.py

# This script looks in the master shapefile for a series of polygons. It then looks for all the DEMs in the directory, finds which DEM
# the polygon is in, and then clips the raster to this polygon. It writes the raster to a tile sub-directory. 

# FJC 30/06/21

import osgeo.gdal as gdal
import geopandas as gpd
import os

def create_raster_from_tile(row):
  """
  loop through the dataframe and create a raster from the geometry of each tile.
  You need to have a column called "location" which tells us the raster that each
  tile corresponds to.
  This writes a temporary shapefile for each tile.
  """
  temp_shp = row.to_file(data_dir+'tile_temp.shp')
  this_tile_fname = 'tile_'+str(row['id'])
  # create a new directory for this tile
  if not os.path.isdir(data_dir+this_tile_fname):
    os.mkdir(data_dir+this_tile_fname)
  # define the destination and source files
  dest_dst = data_dir+this_tile_fname+'/tile_'+str(row['id'])+'.bil' 
  src_dst = data_dir+row['location']
  # do the clipping to a new raster for the tile
  gdal.warp(dest_dst, src_dst, cutline=data_dir+'tile_temp.shp', cropToCutline=True, format='ENVI')

# first read in the tiling shapefile
data_dir = '/raid/fclubb/san_andreas/USGS_Lidar_B3/'
tiles = gpd.read_file(data_dir+'USGS_Lidar_B3_tiles.shp')

# read in the raster footprint shapefile
rasters = gpd.read_file(data_dir+'USGS_Lidar_B3_rasters.shp')

# polygon intersection opereation to find out the raster corresponding to each tile.
tiles_plus_rasters = gpd.sjoin(tiles, rasters, how='left', op='intersects')

# this should give a new column with the raster that each tile corresponds to. now use gdal to do the tiling
tiles_plus_rasters.apply(create_raster_from_tile, axis=1)