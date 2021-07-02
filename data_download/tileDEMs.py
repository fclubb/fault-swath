# tileDEMs.py

# This script looks in the master shapefile for a series of polygons. It then looks for all the DEMs in the directory, finds which DEM
# the polygon is in, and then clips the raster to this polygon. It writes the raster to a tile sub-directory. 

# FJC 30/06/21

#import osgeo.gdal as gdal
import pandas as pd
import geopandas as gpd
import os
import numpy as np
import subprocess

def create_raster_from_tile(row):
  """
  loop through the dataframe and create a raster from the geometry of each tile.
  You need to have a column called "location" which tells us the raster that each
  tile corresponds to.
  This writes a temporary shapefile for each tile.
  """
  if(str(row['raster']) != 'nan'):
    print("TILE ID: ", row['id'])
    sub_df = tiles_plus_rasters.iloc[[row.name]]
    temp_shp = sub_df.to_file(data_dir+'tile_temp.shp')
    this_tile_fname = 'tile_'+str(int(row['id']))
    # create a new directory for this tile
    if not os.path.isdir(data_dir+this_tile_fname):
      os.mkdir(data_dir+this_tile_fname)
    # define the destination and source files
    dest_dst = data_dir+this_tile_fname+'/tile_'+str(int(row['id']))+'.bil'
    raster_name = str(row['raster'])
    src_dst = data_dir+raster_name[:len(raster_name) - 5]+'.'+row['ext']
    print(dest_dst, src_dst)
    # do the clipping to a new raster for the tile
    gdal_cmd = 'gdalwarp -cutline {}tile_temp.shp -crop_to_cutline -dstnodata -9999 -of ENVI {} {}'.format(data_dir, src_dst, dest_dst)
    #subprocess.run(gdal_cmd, check=True, shell=True)
    print(row['ext'])
    #gdal.Warp(dest_dst, src_dst, cutlineDSName=data_dir+'tile_temp.shp', cropToCutline=True, format='ENVI')

# first read in the tiling shapefile
data_dir = '/raid/fclubb/san_andreas/USGS_LidarB3/'
tiles = gpd.read_file(data_dir+'USGS_LidarB3_tiles.shp')

if not os.path.isfile(data_dir+'USGS_LidarB3_Rasters.shp'):
  if not os.path.isfile(data_dir+'USGS_LPC_CA_Sonoma_1m_DTM_mask.shp'):
    # loop through the rasters in the directory and create a footprint shapefile, then merge
    for f in os.listdir(data_dir):
      if f.endswith('.tif') or f.endswith('.bil'):
        prefix = f.split('.')[0]
        print(prefix)
        # first gdal command - create a masking raster with an alpha band
        gdal_cmd = 'gdalwarp -dstnodata 0 -dstalpha -of GTiff {} {}_mask.tif'.format(data_dir+f, data_dir+prefix)
        print(gdal_cmd)
        subprocess.run(gdal_cmd, check=True, shell=True)
        # second step - polygonise masking raster to shapefile
        gdal_cmd = 'gdal_polygonize.py {}_mask.tif -b 2 -f "ESRI Shapefile" {}_mask.shp'.format(data_dir+prefix, data_dir+prefix)
        print(gdal_cmd)
        subprocess.run(gdal_cmd, check=True, shell=True)

  # merge these to a single shapefile for all the rasters
  ogr_cmd = 'ogrmerge.py -o {}USGS_LidarB3_Rasters.shp {}*_mask.shp -single -f "ESRI Shapefile" -src_layer_field_name raster'.format(data_dir, data_dir)
  subprocess.run(ogr_cmd, check=True, shell=True)

# read in the raster footprint shapefile
rasters = gpd.read_file(data_dir+'USGS_LidarB3_Rasters.shp')

# polygon intersection opereation to find out the raster corresponding to each tile.
tiles_plus_rasters = gpd.sjoin(tiles, rasters, how='left', op='intersects')
tiles_plus_rasters.drop(columns='DN', inplace=True)
tiles_plus_rasters = tiles_plus_rasters[tiles_plus_rasters['raster'].str.len() > 3]
tiles_plus_rasters.drop_duplicates(subset='id', keep='first', inplace=True)
tiles_plus_rasters['raster'] = tiles_plus_rasters['raster'].str[:-5]
# filter the DF to only have one row for each tile (remove overlapping rasters)
#print(tiles_plus_rasters)

# loop through the directory and add a column for whether it's a bil or a tif
raster_list = []
ext_list = []
for f in os.listdir(data_dir):
  if f.endswith('.tif') or f.endswith('.bil'):
    if not 'mask' in f:
      raster_list.append(f.split('.')[0])
      ext_list.append(f.split('.')[1])

dir_rasters = pd.DataFrame({'raster': raster_list, 'ext': ext_list})
del raster_list, ext_list
#print(dir_rasters)
# merge with the tiles df
tiles_plus_rasters = tiles_plus_rasters.merge(dir_rasters, left_on='raster', right_on='raster')
#print(tiles_plus_rasters)
	

# this should give a new column with the raster that each tile corresponds to. now use gdal to do the tiling
tiles_plus_rasters.apply(create_raster_from_tile, axis=1)
