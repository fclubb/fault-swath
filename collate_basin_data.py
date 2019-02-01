# Read in the basin data and collate it into some catchment averaged statistics
# The basins themseleves are a raster - read this in and get the centroid of each
# to use as the latitude and longitude.
# The info is from the hilltop flow routing algorithms as a csv file. use pandas
# to read this in and then get statistics for each basin. Output to a csv which will
# be appended for the whole fault.

import pandas as pd
import numpy as np

# shapefiles
import rasterio
from rasterio.features import shapes
from shapely.geometry import shape, Polygon, mapping, Point
import fiona

def PolygoniseRaster(DataDirectory, RasterFile, OutputShapefile='polygons'):
    """
    This function takes in a raster and converts to a polygon shapefile using rasterio
    from https://gis.stackexchange.com/questions/187877/how-to-polygonize-raster-to-shapely-polygons/187883#187883?newreg=8b1f507529724a8488ce4789ba787363

    Args:
        DataDirectory (str): the data directory with the basin raster
        RasterFile (str): the name of the raster
        OutputShapefile (str): the name of the output shapefile WITHOUT EXTENSION. Default = 'polygons'

    Returns:
        Dictionary where key is the raster value and the value is a shapely polygon

    Author: FJC
    """
    # import modules

    # define the mask
    #mask = None
    raster_band = 1

    # get raster no data value
    NDV = getNoDataValue(DataDirectory+RasterFile)

    # load in the raster using rasterio
    with rasterio.open(DataDirectory+RasterFile) as src:
        image = src.read(raster_band, masked=False)

        msk = src.read_masks(1)

        results = (
        {'properties': {'raster_val': v}, 'geometry': s}
        for i, (s, v)
        in enumerate(
            shapes(image, mask=msk, transform=src.transform)))

    # define shapefile attributes
    # crs = src.crs.wkt
    # print (crs)
    crs = GetUTMEPSG(DataDirectory+RasterFile)
    schema = {'geometry': 'Polygon',
              'properties': { 'ID': 'float'}}



    # This is necessary to filter the basin results
    geoms = list(results)
    #print("Geom size is: "+str(len(geoms)))

    filtered_geoms = {}
    area_dict = {}
    for f in geoms:
        this_shape = Polygon(shape(f['geometry']))
        this_val = float(f['properties']['raster_val'])
        #print("ID is: "+str(this_val))
        this_area = this_shape.area
        if this_val in filtered_geoms.keys():
            print("Whoops. Found a repeated ID. Getting rid of the smaller one.")
            if area_dict[this_val] < this_area:
                filtered_geoms[this_val] = f
                area_dict[this_val] = this_area
                print("Found a repeated ID. Keeping the one with area of "+str(this_area))
            else:
                print("Keeping the initial ID.")
        else:
            filtered_geoms[this_val] = f
            area_dict[this_val] = this_area

    new_geoms = []
    for key,item in filtered_geoms.items():
        this_shape = Polygon(shape(item['geometry']))
        this_val = float(item['properties']['raster_val'])
        #print("ID is: "+str(this_val))
        this_area = this_shape.area
        #print("Area is: "+str(this_area))
        new_geoms.append(item)
    #print("Geom size is: "+str(len(new_geoms)))

    # transform results into shapely geometries and write to shapefile using fiona
    PolygonDict = {}
    with fiona.open(DataDirectory+OutputShapefile, 'w', crs=crs, driver='ESRI Shapefile', schema=schema) as output:
        for f in new_geoms:
            this_shape = Polygon(shape(f['geometry']))
            this_val = float(f['properties']['raster_val'])
            print("ID is: "+str(this_val))
            if this_val != NDV: # remove no data values
                output.write({'geometry': mapping(this_shape), 'properties':{'ID': this_val}})
            PolygonDict[this_val] = this_shape

    return PolygonDict

def GetBasinOutlines(DataDirectory, basins_fname):
	"""
	This function takes in the raster of basins and gets a dict of basin polygons,
	where the key is the basin ID and the value is a shapely polygon of the basin.

	Args:
		DataDirectory (str): the data directory with the basin raster
		basins_fname (str): the basin raster

	Returns:
		list of shapely polygons with the basins

	Author: FJC
	"""
	# read in the basins raster
	this_fname = basins_fname.split('.')
	print(basins_fname)
	OutputShapefile = this_fname[0]+'.shp'

	# polygonise the raster
	BasinDict = LSDMap_IO.PolygoniseRaster(DataDirectory, basins_fname, OutputShapefile)
	return BasinDict

def GetBasinCentroids(DataDirectory, basins_fname):
	"""
	This function takes in the raster of basins and returns a dict where the
	key is the basin ID and the value is the shapely point of the centroid

	Args:
		DataDirectory (str): the data directory with the basin raster
		basins_fname (str): the basin raster

	Returns:
		dict of centroid points

	Author: FJC
	"""
	# get the basin polygons
	BasinDict = GetBasinOutlines(DataDirectory, basins_fname)

	# get the centroids
	CentroidDict = {}
	for basin_key, basin in BasinDict.iteritems():
		CentroidDict[basin_key] = Point(basin.centroid)

	return CentroidDict

def CollateBasinStatistics(DataDirectory, fname_prefix, basins_fname):
    """
    This function reads in the hilltop data: "fname_prefix_HilltopData.csv"
    and gets some statistics about the hilltops for each basin in the
    basin raster. This is written to an output csv which has the basin ID,
    the latitude and longitude of the basin centroid, and the statistics.

    Args:
        DataDirectory (str): the data directory with the basin raster
        fname_prefix (str): the prefix for the DEM
        basins_fname (str): the name of the basin raster

    Returns:
       outputs csv of basin statistics

    Author: FJC
    """
    # read in the csv file with the hilltop data
    df = pd.read_csv(DataDirectory+fname_prefix+'_HilltopData.csv')

    # get the basin centroids
    CentroidDict = GetBasinCentroids(DataDirectory, basins_fname)

    # loop through each basin id
    for id, centroid in CentroidDict.items():
        # mask the dataframe for this basin
        this_df = df[df['BasinID'] == id]

        # get the latitude and longitude of the centroid
        # original coords are easting and northing so we need to convert...
        print(centroid.x, centroid.y)

        # # slope
        # this_median = this_df.S.median()
        # median_S.append(this_median)
        # S_lowerP = np.percentile(this_df.S, 16)
        # S_upperP = np.percentile(this_df.S, 84)
        # S_lower_err.append(this_median-S_lowerP)
        # S_upper_err.append(S_upperP-this_median)
        #
        # # hilltop curvature
        # this_median = abs(this_df.Cht.median())
        # median_cht.append(this_median)
        # cht_lowerP = np.percentile(this_df.Cht, 16)
        # cht_upperP = np.percentile(this_df.Cht, 84)
        # cht_lower_err.append(this_median-abs(cht_upperP)) # these are the opposite way round because
        # cht_upper_err.append(abs(cht_lowerP)-this_median) # I am inverting the axis to show positive Chts
        #
        # # hillslope length
        # this_median = this_df.Lh.median()
        # median_Lh.append(this_median)
        # Lh_lowerP = np.percentile(this_df.Lh, 16)
        # Lh_upperP = np.percentile(this_df.Lh, 84)
        # Lh_lower_err.append(this_median-Lh_lowerP)
        # Lh_upper_err.append(Lh_upperP-this_median)
        #
        # # R Star
        # this_median = this_df.R_Star.median()
        # median_Rstar.append(this_median)
        # Rstar_lowerP = np.percentile(this_df.R_Star, 16)
        # Rstar_upperP = np.percentile(this_df.R_Star, 84)
        # Rstar_lower_err.append(this_median-Rstar_lowerP)
        # Rstar_upper_err.append(Rstar_upperP-this_median)
        #
        # # E Star
        # this_median = this_df.E_Star.median()
        # median_Estar.append(this_median)
        # Estar_lowerP = np.percentile(this_df.E_Star, 16)
        # Estar_upperP = np.percentile(this_df.E_Star, 84)
        # Estar_lower_err.append(this_median-Estar_lowerP)
        # Estar_upper_err.append(Estar_upperP-this_median)
