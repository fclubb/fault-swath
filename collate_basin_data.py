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
import rasterio
from rasterio.features import shapes
from shapely.geometry import shape, Polygon, mapping, Point
import fiona

def getNoDataValue(rasterfn):
    """This gets the nodata value from the raster

    Args:
        rasterfn (str): The filename (with path and extension) of the raster

    Returns:
        float: nodatavalue; the nodata value

    Author: SMM
    """
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    return band.GetNoDataValue()

# This gets the UTM zone, if it exists
def GetUTMEPSG(FileName):
    """Uses GDAL to get the EPSG string from the raster.

    Args:
        FileName (str): The filename (with path and extension) of the raster.

    Return:
        str: The EPSG string

    Author: SMM
    """
    if exists(FileName) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + FileName + '\'')

    # see if the file exists and get the dataset
    SourceDS = gdal.Open(FileName, gdal.GA_ReadOnly)
    if SourceDS == None:
        raise Exception("Unable to read the data file")

    EPSG_string = 'NULL'

    # get the projection
    print("Let me get that projection for you")
    prj=SourceDS.GetProjection()
    srs=osr.SpatialReference(wkt=prj)

    if srs.IsProjected:
        #print("Trying projcs")
        #print(str(srs.GetAttrValue(str('PROJCS'),0)))


        print(srs.GetAttrValue(str('projcs')))
        proj_str = srs.GetAttrValue(str('projcs'))
        print("The projection string is: "+proj_str)

        print(proj_str)


        if proj_str != None:

            # extract the UTM information
            if "UTM Zone" in proj_str:
                first_split = proj_str.split(',')
                first_half = first_split[0]
                second_half = first_split[1]
                if "Northern" in second_half:
                    N_or_S = "N"
                else:
                    N_or_S = "S"
                second_split = first_half.split(' ')
                zone = second_split[2]

            elif "_Hemisphere" in proj_str:
                if "Southern" in proj_str:
                    N_or_S = "S"
                else:
                    N_or_S = "N"
                proj_split = proj_str.split('_')
                zone = proj_split[2]
                
            
            else:

                proj_split = proj_str.split('_')
                zone = proj_split[-1]

                N_or_S = zone[-1]
                zone = zone[:-1]


            # adding some logic for zones < 10
            if len(zone) < 2:
                zone = '0'+zone
                
            EPSG_string = 'epsg:'
            if N_or_S == 'S':
                EPSG_string = EPSG_string+'327'+zone
            else:
                EPSG_string = EPSG_string+'326'+zone
            print("The EPSG string is: "+EPSG_string)
    else:
        raise Exception("This is not a projected coordinate system!")

    return EPSG_string

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
	BasinDict = PolygoniseRaster(DataDirectory, basins_fname, OutputShapefile)
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
	for basin_key, basin in BasinDict.items():
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
    this_epsg = GetUTMEPSG(DataDirectory+basins_fname)
    wgs84 = Proj('+init=EPSG:4326', preserve_flags=True)
    this_proj = Proj('+init='+this_epsg)

    # lists for writing output
    lon_centroid, lat_centroid, median_S, S_lower_err, S_upper_err, median_cht, cht_lower_err, cht_upper_err, median_Lh, Lh_upper_err, Lh_lower_err, median_Rstar, Rstar_lower_err, Rstar_upper_err, median_Estar, Estar_lower_err, Estar_upper_err =  ([] for i in range(17))

    # loop through each basin id
    for id, centroid in CentroidDict.items():
        # mask the dataframe for this basin
        this_df = df[df['BasinID'] == id]

        # get the latitude and longitude of the centroid
        # original coords are easting and northing so we need to convert...
        #print(centroid.x, centroid.y)
        lon, lat = transform(this_proj, wgs84, centroid.x, centroid.y)
        #print(lon, lat)
        lon_centroid.append(lon)
        lat_centroid.append(lat)

        # slope
        this_median = this_df.S.median()
        median_S.append(this_median)
        S_lowerP = np.percentile(this_df.S, 16)
        S_upperP = np.percentile(this_df.S, 84)
        S_lower_err.append(this_median-S_lowerP)
        S_upper_err.append(S_upperP-this_median)
        
        # hilltop curvature
        this_median = abs(this_df.Cht.median())
        median_cht.append(this_median)
        cht_lowerP = np.percentile(this_df.Cht, 16)
        cht_upperP = np.percentile(this_df.Cht, 84)
        cht_lower_err.append(this_median-abs(cht_upperP)) # these are the opposite way round because
        cht_upper_err.append(abs(cht_lowerP)-this_median) # I am inverting the axis to show positive Chts
        
        # hillslope length
        this_median = this_df.Lh.median()
        median_Lh.append(this_median)
        Lh_lowerP = np.percentile(this_df.Lh, 16)
        Lh_upperP = np.percentile(this_df.Lh, 84)
        Lh_lower_err.append(this_median-Lh_lowerP)
        Lh_upper_err.append(Lh_upperP-this_median)
        
        # R Star
        this_median = this_df.R_Star.median()
        median_Rstar.append(this_median)
        Rstar_lowerP = np.percentile(this_df.R_Star, 16)
        Rstar_upperP = np.percentile(this_df.R_Star, 84)
        Rstar_lower_err.append(this_median-Rstar_lowerP)
        Rstar_upper_err.append(Rstar_upperP-this_median)
        
        # E Star
        this_median = this_df.E_Star.median()
        median_Estar.append(this_median)
        Estar_lowerP = np.percentile(this_df.E_Star, 16)
        Estar_upperP = np.percentile(this_df.E_Star, 84)
        Estar_lower_err.append(this_median-Estar_lowerP)
        Estar_upper_err.append(Estar_upperP-this_median)
    
    output_list = [('basin_id', list(CentroidDict.keys())),
                   ('longitude_centroid', lon_centroid),
                   ('latitude_centroid', lat_centroid),
                   ('S_median', median_S),
                   ('S_lower_err', S_lower_err),
                   ('S_upper_err', S_upper_err),
                   ('Lh_median', median_Lh),
                   ('Lh_lower_err', Lh_lower_err),
                   ('Lh_upper_err', Lh_upper_err),
                   ('cht_median', median_cht),
                   ('cht_lower_err', cht_lower_err),
                   ('cht_upper_err', cht_upper_err),
                   ('Rstar_median', median_Rstar),
                   ('Rstar_lower_err', Rstar_lower_err),
                   ('Rstar_upper_err', Rstar_upper_err),
                   ('Estar_median', median_Estar),
                   ('Estar_lower_err', Estar_lower_err),
                   ('Estar_upper_err', Estar_upper_err)]

    # write output to csv
    OutDF = pd.DataFrame.from_items(output_list)
    csv_outname = DataDirectory+fname_prefix+'_basin_data.csv'
    OutDF.to_csv(csv_outname,index=False)

DataDirectory = '/raid/fclubb/san_andreas/SouthernSAF/SouthernSAF_17/'
fname_prefix = 'SouthernSAF_17'
basins_fname = fname_prefix+'_basins.bil'
CollateBasinStatistics(DataDirectory, fname_prefix, basins_fname)
