#!/usr/bin/env bash

# copy the cluster driver to each directory
subdirs=$(find -mindepth 1 -maxdepth 1 -type d -name '*tile_*')
fname="river_metrics.driver"

for dir in $subdirs
do
  echo $dir
 # if [ ! -e "$dir/$fname" ]; then
  cp $fname $dir 
  #fi
  path="$PWD"
  echo $path


  dir_fname="${dir##*/}"
  # edit the lines in the driver to point to this directory and filename
  sed -i "7s#.*#read path: $path/$dir_fname/#" $dir/$fname
  sed -i "8s#.*#write path: $path/$dir_fname/#" $dir/$fname 
  sed -i "9s#.*#read fname: $dir_fname#" $dir/$fname
  sed -i "10s#.*#write fname: $dir_fname#" $dir/$fname
done
#echo $subdirs

