#!/usr/bin/env bash

# read in the list of tiles and make a directory for each one

file1="/raid/fclubb/san_andreas/NorthernSAF/all_tiles/NorthernSAF_basins.shp"
file2="/raid/fclubb/san_andreas/SouthernSAF/all_tiles_merged/SouthernSAF_basins.shp"


output_shp='/raid/fclubb/san_andreas/SAF_combined/SAF_all_basins.shp'
echo $output_shp

ogr2ogr -f 'ESRI Shapefile' -update -append $output_shp $file1 $file2

