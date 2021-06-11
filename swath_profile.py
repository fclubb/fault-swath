# Read in a shapefile of a fault and plot river profile morphology along it.
# FJC 26/11/18

# set backend to run on server
#import matplotlib
#matplotlib.use('Agg')

# general modules
import numpy as np
import numpy.ma as ma
import pandas as pd
import math
import matplotlib.pyplot as plt
import os
import time
from matplotlib import rcParams
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from inspect import signature

# geospatial modules
from fiona import collection
from shapely.geometry import shape, LineString, mapping
from shapely.geometry import Point as shapelyPoint
import pyproj as pyproj
from geopy.distance import distance as GeoPyDist
import geopandas as gpd

# peak detection
import peakutils
from peakutils.plot import plot as pplot

# Set up fonts for plots
label_size = 12
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = label_size
plt.rc('axes', titlesize=10)     # fontsize of the axes title

def azimuth(point1, point2):
    """
    azimuth between 2 shapely points (interval 0 - 180°, for my work)
    Only works with increasing X (point 2 must be E of point 1).
    """
    angle = math.atan2(point2.x - point1.x, point2.y - point1.y)
    return math.degrees(angle)

def deflection(row, fault_pts):
    """
    Use a row of a basin dataframe to calculate the deflection angle relative to
    the fault strike.
    """
    # get the basin orientation vector
    point1 = row['centroids']
    point2 = [row['longitude'], row['latitude']]
    basin_vec = np.array([(row['longitude'] - row['centroids'].x), (row['latitude'] - row['centroids'].y)])

    # find the nearest point along the fault to this basin.
    fault_dist = fault_pts['distance']
    fault_y = fault_pts['geometry'].y
    fault_x = fault_pts['geometry'].x
    idx = find_nearest_index(fault_dist, row['fault_dist'])
    nearest_dist = fault_pts['distance'][idx]
    #print(row['fault_dist'], nearest_dist)

    # get the fault strike vector
    # calculate nearest fault point - basin point. If negative, get the preceding point. If positive, get the following point.
    #if (nearest_dist - row['fault_dist'] < 0):
    #    fault_vec = np.array([(fault_x[idx] - fault_x[idx-1]), (fault_y[idx] - fault_y[idx-1])])
    #else:
    fault_vec = np.array([(fault_x[idx+1] - fault_x[idx]), (fault_y[idx+1] - fault_y[idx])])
    #print(fault_vec)

    # we have the two vectors. now calculate the dot product between them to find theta_d
    # theta_d = arccos(abs((basin_vec . fault_vec)/magnitude(basin_vec) * magnitude(fault_vec)|))
    theta_d = math.acos(abs(basin_vec.dot(fault_vec)/(np.linalg.norm(basin_vec)*np.linalg.norm(fault_vec))))

    # return the deflection metric in degrees
    return(math.degrees(theta_d))

def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

def find_vicenty_distance_along_line(line):
    """
    find the vicenty distance of the coords along a line
    This is a bit slow because we need to calculate the distance for
    each point on the line
    """
    coords = line.coords
    start_point = coords[0]
    temp_dist=0
    distances=[]
    for i in range (len(coords)):
        #print(coords[i][0])
        if i == 0:
            dist = GeoPyDist((coords[i][1], coords[i][0]), (start_point[1], start_point[0])).km
        else:
            dist = GeoPyDist((coords[i][1], coords[i][0]), (coords[i-1][1], coords[i-1][0])).km
        temp_dist+=dist
        distances.append(temp_dist)
    return distances

def gaussian_weighted_average(x, y, power=100., lenscale=3):
    """
    function to compute a gaussian weighted average of an array y with x
    """
    new_y = np.empty(len(y))
    #print(x)
    for i in range(0, len(x)):
        weights= np.exp(-(x-x[i])**2/lenscale)
        #print(weights)
        summe=np.sum(weights*y**power)
        #print(summe)
        new_y[i]=(summe/np.sum(weights))**(1/power)
        #print(new_y[i])

    return new_y

def get_points_along_line(DataDirectory, baseline_shapefile, output_shapefile, n=1024):
    """
    Interpolate a series of points at equal distances along an input line shapefile. Arguments that need to be supplied are:
    * DataDirectory: the directory of the input/output shapefiles
    * baseline_shapefile: the name of the input line shapefile with extension
    * n: total number of points, should be power2: n = 2^m
    * output_shapefile: the name of the output points shapefile with extension
    """

    points = []
    distances = []
    azimuths = []
    # read in the baseline shapefile
    c = collection(DataDirectory+baseline_shapefile, 'r')
    rec = c.next()
    line = LineString(shape(rec['geometry']))
    line_rvs = LineString(list(line.coords)[::-1])
    # get the coordinate system from the input shapefile
    crs = c.crs

    total_distance = line_rvs.length
    # get the spacing based on the total distance and the number of points
    dist = (total_distance/n)
    print("The total distance is", total_distance, ": returning ", n, "points at a spacing of ", dist)
    temp_dist=0
    metric_dist=0
    temp_azimuth=np.nan
    # have a point at the start of the line
    for j in range(n+1):
        point = line_rvs.interpolate(temp_dist)
        if j == 0:
            temp_metric = 0
            temp_azimuth = np.nan
        else:
            #print(list(point.coords))
            # find the distance between this point and the previous point in metres (vicenty)
            temp_metric = GeoPyDist((point.y, point.x), (points[-1].y, points[-1].x)).km
            temp_azimuth = azimuth(points[-1], point)
            print(temp_azimuth)
        metric_dist+=temp_metric
        #print(metric_dist)
        temp_dist+=dist
        points.append(shapelyPoint(point))
        distances.append(metric_dist)
        azimuths.append(temp_azimuth)


    #output schema
    schema={'geometry': 'Point', 'properties': {'distance': 'float', 'id': 'int', 'azimuth': 'float'} }

    # write the points to a shapefile
    with collection(DataDirectory+output_shapefile, 'w', crs=crs, driver='ESRI Shapefile', schema=schema) as output:
        for i in range (n+1):
            #print point
            output.write({'properties':{'distance':distances[i], 'id':i, 'azimuth':azimuths[i]},'geometry': mapping(points[i])})

    return points, distances

def get_distance_along_fault_from_points(DataDirectory, baseline_shapefile, df, output_pts_csv):
    """
    Find the distance along the fault shapefile from a DataFrame
    with lat/long coordinates
    Args:
        df: the dataframe with the points
        output_pts_csv: the csv filename to save the output dataframe with additional column of distances along the fault
    """

    # read in the baseline shapefile
    c = collection(DataDirectory+baseline_shapefile, 'r')
    rec = c.next()
    line = LineString(shape(rec['geometry']))
    line_rvs = LineString(list(line.coords)[::-1])
    # get the coordinate system from the input shapefile
    crs = c.crs
    # get the vincenty distances along the line
    line_dist = find_vicenty_distance_along_line(line_rvs)
    line_lat = [l[1] for l in line_rvs.coords]
    #print(line_lat)
    line_lon = [l[0] for l in line_rvs.coords]

    # for each point find the length along the line of the nearest point
    lon = df.longitude.values
    lat = df.latitude.values
    distances = []
    fault_normal_dists = [] # the distance away from the fault

    for i in range(len(lat)):
        if not np.isnan(lat[i]):
            #create a shapely point
            point = shapelyPoint(lon[i], lat[i])
            # find the distance of the point along the line
            dist = line_rvs.project(point)
            # shapely point on the line nearest to the initial point
            line_pt = line.interpolate(line.project(point))
            # getting the distance - difficult because we need vicenty distance
            # find the nearest point in the line distances (closest latitude)
            idx = find_nearest_index(line_lat, lat[i])
            dist = GeoPyDist((line_lat[idx], line_lon[idx]), (line_pt.y, line_pt.x)).km
            #print(dist)
            if line_lat[idx] < lat[i]:
                # line point is at a lower latitude than the point - need to find the distance between them
                # and minus that from the line distance
                dist = line_dist[idx] - dist
            else:
                dist = line_dist[idx] + dist
                #print(dist)

            distances.append(dist)
            # write the distance away from the fault for each point
            away_dist = GeoPyDist((line_pt.y, line_pt.x), (point.y, point.x)).km
            fault_normal_dists.append(away_dist)
        else:
            distances.append(np.nan)
            fault_normal_dists.append(np.nan)

    # save the distances to csv
    df['fault_dist'] = distances
    df['fault_normal_dist'] = fault_normal_dists
    df.to_csv(output_pts_csv, index=False)
    return df

def get_orthogonal_coefficients(pts):
    """
    For a list of shapely points, get the line orthogonal to each point and then save
    the parameters for the equation of the line where ax + by + c = 0
    """
    from sympy import Point, Line, N
    # loop through the points between start, stop
    coeffs = []
    #print("Getting the coefficients of each line")
    for i in range(0, len(pts)):
        # start of line
        if i == 0:
            x1 = pts[i].x
            y1 = pts[i].y
        # end of line
        if i == len(pts)-1:
            x2 = pts[i].x
            y2 = pts[i].y
        else:
            x1 = pts[i-1].x
            x2 = pts[i+1].x
            y1 = pts[i-1].y
            y2 = pts[i+1].y
        # line between pts[i-1], i+1
        l1 = Line(Point(x1, y1), Point(x2, y2))
        l2 = l1.perpendicular_line(Point(pts[i].x, pts[i].y))
        # a needs to be > 0: get the sign of a and multiply the coefficients by this.
        _ = np.sign(N(l2.coefficients[0]))
        coeffs.append([_*N(l2.coefficients[0]), _*N(l2.coefficients[1]), _*N(l2.coefficients[2])])

    return coeffs

def bisection_method(points, coeffs, distances, df, output_csv):
    """
    Use a bisectioning method (https://en.wikipedia.org/wiki/Bisection_method)
    to find the nearest point along the fault to each channel in the river cluster csv

    Args:
        points: list of shapely points along the fault
        coeffs: list of the coefficients of the line orthogonal to each point
        distances: list of the distances along the fault
        cluster_csv: name of the csv file with the cluster info
    """
    print("CHECKING SOME LENGTHS")
    # read in the csv to a pandas dataframe
    #df = pd.read_csv(cluster_csv)
    cluster_x = df.longitude.values
    cluster_y = df.latitude.values
    alpha = np.empty(len(cluster_x))

    m = int(math.log(len(points),2))
    #print(m)
    # start with the middle baseline point (q = 2^(m-1))
    i = 1
    qs = np.full(len(cluster_x), 2**(int(m-1)))
    #print(qs)
    #print(coeffs)
    while (int(m)-1-i >= -1):
        t0=time.time()
        cs = np.asarray([coeffs[q] for q in qs])
        t1=time.time()
        print("cs executed in: ", t1-t0)
        # calculate alpha: ax+bx+c=alpha for each river point (x and y from river point,
        # a, b, and c from the line)
        t0=time.time()
        alpha = np.sign(cs[:,0]*cluster_x + cs[:,1]*cluster_y + cs[:,2])
        t1=time.time()
        print("alpha executed in: ", t1-t0)
        # if alpha > 0, river point is to the right of the line. If alpha < 0, river point is
        # to the left. If alpha = 0, point is on the line.
        # Take further baseline point depending on whether alpha is positive or negative
        t0=time.time()
        if (m-1-i >= 0):
            qs = qs + alpha*2**(m-1-i)
        else:
            # last iteration
            qs = qs + (alpha-1)/2
            #qs = [q_old + int((a-1)/2) for q_old,a in zip(qs, alpha)]
        t1=time.time()
        print("qs executed in: ", t1-t0)
        print("Iteration", i, "executed in", t1-t0)
        i+=1

    # qs is now the index of the distance of the closest distance along the fault to each point in the cluster dataframe.  # first get the distance that each one represents
    fault_dists = pd.Series([distances[int(q)] for q in qs])
    # append the distances to the dataframe
    df['fault_dist'] = fault_dists.values
    #print(df)

    # use the indices to assign whether the points are east or west of the fault
    # we calculate this based on the equation
    # dir = (x - x_1)(y_2 - y_1) - (y - y_1)(x_2 - x_1)
    # if dir < 0 : point is E of line, if dir > 0 : point is W of line
    dirs = pd.Series([(cluster_x[i] - points[int(q)].x)*(points[int(q+1)].y - points[int(q)].y) - (cluster_y[i] - points[int(q)].y)*(points[int(q+1)].x - points[int(q)].x) for i, q in enumerate(qs)])
    print(dirs)
    df['direction'] = dirs.values


    df.to_csv(output_csv, index=False)
    return df

def percentile(n):
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_

def get_median_slope_in_basins(river_gdf, basin_gdf):
    """
    Function to get the median slope in each basin and return a new dataframe
    with the median and percentiles of the slope, plus latitude and longitude of the
    basin outlet
    This requires passing of two geopandas dataframes: the river gdf and the basin gdf.
    It uses a spatial join to merge them.
    ***NOTE 10/06/21 - this is deprecated and replaced with new function to make sure that
    the IDs are unique to each basin. This needs to be replaced. ***
    """

    gr = df.groupby(['basin_id'])['slope'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'median': 'slope_median', 'std': 'slope_std', 'percentile_25': 'slope_q1', 'percentile_75': 'slope_q2'}).reset_index()

    # now get the lat, long, and fault distance of the outlet
    u = df.groupby('basin_id')['distance_from_outlet'].idxmin()
    dist_df = df.loc[u, ['basin_id','latitude', 'longitude', 'fault_dist']].reset_index(drop=1)
    df_merge = pd.merge(gr, dist_df, on='basin_id')
    print(df_merge)

    return df_merge

def q1(x):
    return x.quantile(0.25)
def q2(x):
    return x.quantile(0.75)

def log_sum(x):
    return np.log(np.sum(10**x))


def process_gps_data(gps_df, threshold_record_length=5, threshold_uplift=5):
    """
    Thin the gps data to remove stations with short record lengths and
    unrealistically high uplift or subsidence rates
    """
    gps_df = gps_df[gps_df['record_length'] > threshold_record_length]
    gps_df = gps_df[gps_df['RU(mm/yr)'] < threshold_uplift]
    gps_df = gps_df[gps_df['RU(mm/yr)'] > -threshold_uplift]
    return gps_df

#--------------------------------------------------------------------------------------------------
# PLOTTING FUNCTIONS
#--------------------------------------------------------------------------------------------------

def plot_channel_slopes_along_fault(DataDirectory, fname_prefix, stream_order, median_river_shp, labels_csv, slip_rate_csv, fault_points, pcp_shp, plate_azimuth=135):
    """
    Function to make plot of channel slopes along the fault compared to various datasets
    """

    # csv with the river profiles
    river_df = gpd.read_file(median_river_shp)

    # set up a figure
    fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(10,15), sharex=True, sharey=False)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

    # plot the channel slope data

    # first, all the slopes east of the fault (direction < 0)
    east_df = river_df[river_df['direction'] < 0]
    # then all the slopes west of the fault (direction > 0)
    west_df = river_df[river_df['direction'] > 0]

    all_dfs = [east_df, west_df]
    titles = ['North American Plate (east of SAF)', 'Pacific Plate (west of SAF)']
    colors = ['r', 'b']
    figlabels = ['a', 'b', 'c', 'd', 'e']
    for i, df in enumerate(all_dfs):

        ax[i].grid(color='0.8', linestyle='--', which='both')
        ax[i].set_ylim(0,0.7)
        ax[i].text(0.04,0.85, titles[i], fontsize=14, transform=ax[i].transAxes, bbox=dict(facecolor='white'))
        ax[i].set_ylabel('Median channel\ngradient (m/m)', labelpad=10, fontsize=14)

        # now group by the fault dist and plot percentages
        gr = df.groupby(['fault_dist'])['channel_sl'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()

        # plot the medians
        ax[i].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.4', mec='0.4', c='0.4', capsize=2, alpha=0.1)

        # rolling median of channel slopes
        slopes_df = gr.sort_values(by='fault_dist')
        slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(10, center=True).median()
        print(df.columns)
        # add lat and lon information to this
        #slopes_df = pd.merge(slopes_df, df[['latitude','longitude']], on=['fault_dist'])
        #print(slopes_df)

        # create a mask for gaps in the median slopes
        these_dists = slopes_df['fault_dist'].values
        mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
        #print(mask_starts)
        mc = ma.array(slopes_df['slope_rollmedian'].values)
        mc[mask_starts] = ma.masked
        ax[i].plot(slopes_df['fault_dist'], mc, c=colors[i], zorder=100, lw=3, ls='--')

        # find and plot peaks in the rolling median
        indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.5, min_dist=30))
        #print(indexes)
        peak_dists = slopes_df['fault_dist'].iloc[indexes]
        peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
        #peak_dists = these_dists[indexes]
        #peak_slopes = mc[indexes]
        #print(peak_dists)
        #print(peak_slopes)
        #print("Channel slope peak distances: ", peak_dists.values)
        ax[i].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
        for j, txt in enumerate(list(peak_dists)):
            ax[i].annotate(str(int(txt))+' km', (list(peak_dists)[j]-10, list(peak_slopes)[j]+0.05), zorder=300, fontsize=10)

        ax[i].text(0.97, 0.9, figlabels[i], bbox=dict(facecolor='white', edgecolor='k'), transform=ax[i].transAxes)

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],0.7), xytext=(labels_dist[i], 0.8), ha='center', fontsize=14, arrowprops=dict(facecolor='k', arrowstyle="->"))

    # plot the slip rates
    slip_df = pd.read_csv(slip_rate_csv)
    ax[2].grid(color='0.8', linestyle='--', which='both')
    ax[2].axvspan(400, 580, facecolor='0.5', alpha=0.6)
    ax[2].errorbar(x=slip_df['fault_dist'], y=slip_df['slip_rate'], yerr=slip_df['slip_rate_u'], fmt='o',ms=8, marker='D', mfc='0.3', mec='k', c='k', capsize=4)
    ax[2].set_ylabel('InSAR-derived right\nlateral slip rate (mm/yr)', fontsize=14)
    #ax[2].set_ylim(0, slip_df['slip_rate'].max()+0.1)

    ax[2].set_xlim(100,1100)
    ax[2].text(0.97, 0.9, figlabels[2], bbox=dict(facecolor='white', edgecolor='k'), transform=ax[2].transAxes)

    # plot the azimuths
    pts = gpd.read_file(DataDirectory+fault_points)

    # calculate azimuth deltas. Decrease in angle = restraining bend, increase = releasing bend
    #az_deltas = pts['azimuth'].diff()
    az_deltas = pts['azimuth'] - plate_azimuth
    ax[3].plot(pts['distance'], az_deltas, c = 'k', lw = 2)
    ax[3].set_ylabel('SAF strike azimuth\nvs plate motion ($^\circ$)', fontsize=14)
    ax[3].axhspan(0, 50, alpha=0.4, color='deepskyblue')
    ax[3].axhspan(-50, 0, alpha=0.4, color='red')
    ax[3].grid(color='0.8', linestyle='--', which='both')
    ax[3].text(120, -35, 'Restraining (uplift proxy)')
    ax[3].text(120, 35, 'Releasing (subsidence proxy)')
    ax[3].set_ylim(-45,45)
    ax[3].invert_yaxis()

    ax[3].set_xlim(100,1100)
    ax[3].text(0.97, 0.9, figlabels[3], bbox=dict(facecolor='white', edgecolor='k'), transform=ax[3].transAxes)
    #plt.ylim(0,0.4)
    plt.xlabel('Distance along fault (km)', fontsize=14)

    # plot the mean annual precipitation
    pcp_gdf = gpd.read_file(pcp_shp)
    pcp_gdf = pcp_gdf[pcp_gdf['rvalue_1'] > 0]
    ax[4].grid(color='0.8', linestyle='--', which='both')
    ax[4].plot(pcp_gdf["distance"], pcp_gdf["rvalue_1"], lw=2, c='purple')
    ax[4].fill_between(pcp_gdf["distance"], pcp_gdf["rvalue_1"], 0, facecolor='purple', alpha=0.3)
    ax[4].set_ylabel('Mean annual\nprecipitation (mm yr$^{-1}$)', fontsize=14)
    ax[4].set_xlim(100,1100)
    ax[4].set_ylim(0, 2100)
    ax[4].text(0.97, 0.9, figlabels[4], bbox=dict(facecolor='white', edgecolor='k'), transform=ax[4].transAxes)

    # save figure
    plt.subplots_adjust(hspace=0.3)
    plt.tight_layout()
    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_slopes_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

    return peak_dists, peak_slopes

def plot_earthquakes_along_fault(DataDirectory, fname_prefix, eq_csv, labels_csv, peak_dists):
    """
    Make a plot of the earthquake frequency and summed magnitude along fault. Plot the channel slope peaks
    on top of this.

    Args:
        eq_csv: csv file with the earthquake data
        labels_csv: point locations of city labels to add
        peak_dists: array of distances along fault of each channel slope peak
    """

    # set up a figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,8), sharex=True, sharey=False)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    # add a grid
    ax[0].grid(color='0.8', linestyle='--', which='both')

    # read in the EQ csv to a dataframe
    eq_df = pd.read_csv(eq_csv)

    # separate by magntiude to make frequency plot
    magnitude_thr = 5
    lower_eq_df = eq_df[eq_df['Mw'] < 6]
    upper_eq_df = eq_df[eq_df['Mw'] >= 6]
    dfs = [lower_eq_df, upper_eq_df]
    colours = ['k', 'b']

    for i, df in enumerate(dfs):
        gr = df.groupby(['fault_dist'])['Mw'].agg(['count']).reset_index()
        gr['normalized_count'] = gr['count']/gr['count'].sum()
        # mc = ma.array(gr['count'].values)
        #
        # # create a mask for data gaps
        # these_dists = gr['fault_dist'].values
        # mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
        # mc[mask_starts] = ma.masked
        # plot the frequency for this magnitude
        ax[0].bar(x=gr['fault_dist'], height=gr['normalized_count'], facecolor=colours[i], zorder=100)

    # finalise the frequency plot
    ax[0].axvspan(400, 580, facecolor='0.5', alpha=0.6)
    #ax[0].set_ylim(-10,)
    ax[0].set_ylabel('Frequency', labelpad=10, fontsize=14)

    #now group to make the magnitude plot
    gr = eq_df.groupby(['fault_dist'])['Mw'].agg(['max',log_sum]).reset_index()
    # create a mask for data gaps
    these_dists = gr['fault_dist'].values
    mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],210), xytext=(labels_dist[i], 240), ha='center', fontsize=10, arrowprops=dict(facecolor='k', arrowstyle="->"))

    # plot the summed earthquake magnitude
    #mw = ma.array(gr['log_sum'].values)
    #mw[mask_starts] = ma.masked
    #ax[1].grid(color='0.8', linestyle='--', which='both')
    #ax[1].plot(gr['fault_dist'], mw, c='b',zorder=100)
    #ax[1].set_ylabel('$\Sigma$ $M_w$', labelpad=10, fontsize=14)

    # plot the maximum magnitude in each fault dist bin
    mmax = ma.array(gr['max'].values)
    mmax[mask_starts] = ma.masked
    ax[1].grid(color='0.8', linestyle='--', which='both')
    markerline, stemlines, baseline = ax[1].stem(gr['fault_dist'], mmax, linefmt='blue', use_line_collection=True)
    markerline.set_markerfacecolor('blue')
    markerline.set_markeredgecolor('k')
    stemlines.set_linewidth(0.5)
    ax[1].set_ylabel('Maximum $M_w$', labelpad=10, fontsize=14)
    ax[1].set_ylim(2,)

    # grey bar for creeping segment
    ax[1].axvspan(400, 580, facecolor='0.5', alpha=0.6)
    ax[1].set_xlabel('Distance along fault (km)', fontsize=14)

    ax[1].set_xlim(100,1100)

    # add the channel slope peaks
    bbox_props = dict(boxstyle="circle,pad=0.3", fc="white", ec="k", lw=2)
    for j, txt in enumerate(list(peak_dists)):
        if j != 0:
            print(txt)
            ax[1].annotate(str(int(j)), (list(peak_dists)[j], 8.8), zorder=300, bbox=bbox_props, annotation_clip=False, ha='center')
            #ax[0].vlines(peak_dists[j], -20, -10, colors='red')
            ax[1].vlines(peak_dists[j], 8, 8.8, colors='red')


    # save the figure
    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_EQs.png', dpi=300)
    plt.clf()

def plot_channel_slopes_vs_earthquakes(DataDirectory, fname_prefix, river_csv, eq_csv, labels_csv, peak_dists):
    """
    This function makes a plot of the channel slope data vs. earthquake magnitude and frequency along fault
    """
    # csv with the river profiles
    river_df = pd.read_csv(river_csv)
    #remove negative channel slopes
    river_df = river_df[river_df['slope'] > 0]

    # set up a figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,8), sharex=True, sharey=False)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

    # read in the EQ csv to a dataframe
    eq_df = pd.read_csv(eq_csv)

    # now make a plot of the channel slope frequencies along fault
    shallow_river_df = river_df[river_df['slope'] < 0.2]
    steep_river_df = river_df[river_df['slope'] >= 0.2]

    dfs = [lower_eq_df, upper_eq_df]
    colours = ['k', 'b']

    for i, df in enumerate(dfs):
        gr = df.groupby(['fault_dist'])['Mw'].agg(['count']).reset_index()
        gr['normalized_count'] = gr['count']/gr['count'].sum()
        # mc = ma.array(gr['count'].values)
        #
        # # create a mask for data gaps
        # these_dists = gr['fault_dist'].values
        # mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
        # mc[mask_starts] = ma.masked
        # plot the frequency for this magnitude
        ax[0].bar(x=gr['fault_dist'], height=gr['normalized_count'], facecolor=colours[i], zorder=100)


def plot_basin_orientation_along_fault(DataDirectory, fname_prefix, basins, baseline_shapefile, baseline_points, labels_csv):
    """
    This function reads in a shapefile of the basins and then plots their orientation
    compared to the fault strike
    """
    # read in the shapefile of fault points
    fault_pts = gpd.read_file(DataDirectory+baseline_points)

    # check if you have calculated basin deflection, and if not then calculate it
    basins_file = DataDirectory+fname_prefix+'_basins_deflection.shp'
    if not os.path.isfile(basins_file):
        # get the basins shapefile
        gdf = gpd.read_file(basins)
        gdf['centroids'] = gdf['geometry'].centroid
        gdf['basin_area'] = gdf['geometry'].to_crs('epsg:32610').area
        gdf = gdf.rename(columns={"latitude_o": "latitude", "longitude_": "longitude"})

        # check if you have already calculated the distance along the fault for each basin.
        output_basin_csv = DataDirectory+fname_prefix+'_basins_WGS84_dist.csv'
        if not os.path.isfile(output_basin_csv):
            # change the column names to latitue and longitude so they can be read properly
            print(gdf.columns)
            points = list(fault_pts['geometry'])
            #print(points)
            distances = fault_pts['distance']
            coeffs = get_orthogonal_coefficients(points)
            basin_df = bisection_method(points, coeffs, distances, gdf, output_basin_csv)
            print(basin_df)
        else:
            basin_df = pd.read_csv(output_basin_csv)

        # merge the basin gdf with the df
        gdf['fault_dist'] = basin_df['fault_dist']
        gdf['direction'] = basin_df['direction']
        #print(fault_pts)

        gdf['deflection'] = gdf.apply(deflection, axis=1, fault_pts=fault_pts)
        gdf = gpd.GeoDataFrame(gdf[['basin_id', 'basin_area', 'azimuth', 'deflection', 'latitude', 'longitude', 'fault_dist', 'direction']], geometry=gdf['geometry'], crs='EPSG:4326')
        gdf.to_file(DataDirectory+fname_prefix+'_basins_deflection.shp')
    else:
        gdf = gpd.read_file(basins_file)

    # now do the plotting along fault
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,6), sharex=True)
    ax = ax.ravel()
    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    colors = ['r', 'b']
    # set the default colourmap
    plt.rc('image', cmap='gray')

    # plot the channel slope data
    titles = ['North American Plate', 'Pacific Plate']
    for i, title in enumerate(titles):
        # set up each subplot
        ax[i].grid(color='0.8', linestyle='--', which='both')
        ax[i].set_ylim(0,90)
        ax[i].text(0.03,0.1, titles[i], fontsize=12, transform=ax[i].transAxes, bbox=dict(facecolor='white'))
        ax[i].set_ylabel('$\\theta_d$ ($^\circ$)')

        # first, all the basins east of the fault (direction < 0)
        if i == 0:
            this_df = gdf[gdf['direction'] < 0]
        else:
        # then all the basins west of the fault (direction > 0)
            this_df = gdf[gdf['direction'] > 0]

        # now group by the fault dist and plot percentages
        #gr2 = this_df.groupby(['fault_dist'])['deflection'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()

        # group by the fault dist and get the mean deflection at each distance weighted by drainage area (larger areas = higher weights)
        g = this_df.groupby(['fault_dist'])
        gr = g.apply(lambda x: pd.Series(np.average(x['deflection'], weights=x['basin_area'])).reset_index(name='Deflection_weighted')).reset_index()
        #gr = g.apply(lambda x: pd.Series(np.average(x['deflection'], weights=None)).reset_index(name='Deflection_weighted')).reset_index()
        area = this_df.groupby(['fault_dist'])['basin_area'].median()
        print(area)
        ax[i].scatter(x=gr['fault_dist'], y=gr['Deflection_weighted'], c='0.5', edgecolor='k', zorder=2, s=area/1000, marker='D', alpha=0.2)
        ax[i].axvspan(400, 580, facecolor='0.8', alpha=0.6, zorder=1)
        #ax[i].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=2, c='0.5', alpha=0.2, capsize=2, zorder=1)

        # rolling median
        slopes_df = gr.sort_values(by='fault_dist')
        slopes_df['rollmedian'] = slopes_df['Deflection_weighted'].rolling(10, center=True).median()
        print(slopes_df)

        # create a mask for gaps in the median slopes
        these_dists = slopes_df['fault_dist'].values
        mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
        print(mask_starts)
        mc = ma.array(slopes_df['rollmedian'].values)
        mc[mask_starts] = ma.masked
        ax[i].plot(slopes_df['fault_dist'], mc, c=colors[i], zorder=100, lw=3, ls='-')

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],90), xytext=(labels_dist[i], 110), ha='center', fontsize=10, arrowprops=dict(facecolor='k', arrowstyle="->"))

    # save the figure
    plt.xlabel('Distance along fault (km)')
    plt.savefig(DataDirectory+fname_prefix+'_basins_deflection.png', dpi=300)
    plt.clf()

def plot_channel_slopes_multiple_SO(DataDirectory, fname_prefix, labels_csv):
    """
    Read in a csv file with the channel slopes and plot compared to distance
    along the fault for multiple stream orders (1-4)
    """

    stream_orders = np.arange(1,5,1)
    print("STREAM ORDERS:", stream_orders)
    # set up a figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,8), sharex=True, sharey=True)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    colors = ['r', 'b', 'g', 'orange']

    # plot the channel slope data
    titles = ['North American Plate', 'Pacific Plate']
    for i, title in enumerate(titles):

        # set up each subplot
        ax[i].grid(color='0.8', linestyle='--', which='both')
        ax[i].set_ylim(0,0.7)
        ax[i].text(0.04,0.85, titles[i], fontsize=12, transform=ax[i].transAxes, bbox=dict(facecolor='white'))

        # now add the stream order data
        for j, so in enumerate(stream_orders):
            # csv with the river profiles
            river_df = pd.read_csv(DataDirectory+fname_prefix+"_profiles_fault_dist_SO{}.csv".format(int(so)))
            print("This stream order:", so)

            #remove negative channel slopes
            river_df = river_df[river_df['slope'] > 0]

            # first, all the slopes east of the fault (direction < 0)
            if i == 0:
                df = river_df[river_df['direction'] < 0]
            else:
            # then all the slopes west of the fault (direction > 0)
                df = river_df[river_df['direction'] > 0]

            slope_df = get_median_slope_in_basins(df)

            # now group by the fault dist and plot percentages
            gr = slope_df.groupby(['fault_dist'])['slope_median'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()
            # #print(gr)

            # plot the medians
            #ax[i].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.4', mec='0.4', c='0.4', capsize=2, alpha=0.1)

            # rolling median of channel slopes
            slopes_df = gr.sort_values(by='fault_dist')
            slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(10, center=True).median()

            # create a mask for gaps in the median slopes
            these_dists = slopes_df['fault_dist'].values
            mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
            print(mask_starts)
            mc = ma.array(slopes_df['slope_rollmedian'].values)
            mc[mask_starts] = ma.masked
            ax[i].plot(slopes_df['fault_dist'], mc, c=colors[j], zorder=100, lw=3, label="{}".format(int(so)))

            # find and plot peaks in the rolling median
            # indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.35, min_dist=50))
            # print(indexes)
            # peak_dists = slopes_df['fault_dist'].iloc[indexes]
            # peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
            # #peak_dists = these_dists[indexes]
            # #peak_slopes = mc[indexes]
            # print(peak_dists)
            # print(peak_slopes)
            # #print("Channel slope peak distances: ", peak_dists.values)
            # ax[i].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
            # for j, txt in enumerate(list(peak_dists)):
            #     ax[i].annotate(str(int(txt))+' km', (list(peak_dists)[j], list(peak_slopes)[j]+0.05), zorder=300)

        #gr.plot.scatter(x='fault_dist', y='median')
    plt.ylabel('Median channel slope (m/m)', labelpad=20)
    #ax[0].set_xlim(100,580)
    ax[0].legend(loc='upper right', title="Stream order")
    #plt.show()
    #gr.to_csv(DataDirectory+threshold_lvl+fname_prefix+'_channel_slope_fault_dist.csv')

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],0.7), xytext=(labels_dist[i], 0.8), ha='center', fontsize=10, arrowprops=dict(facecolor='k', arrowstyle="->"))

    plt.xlim(100,1100)
    #plt.ylim(0,0.4)
    plt.xlabel('Distance along fault (km)')

    # save the figure
    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_slopes_EW_multiple_SO.png', dpi=300)
    plt.clf()

def plot_slopes_vs_azimuth(DataDirectory, fname_prefix, stream_order, river_csv, fault_points, plate_azimuth=135):
    """
    Make a scatter plot of the difference in azimuth along fault vs the channel slope.
    Positive delta azimuth = releasing bend, negative delta azimuth = restraining bend
    """

    # csv with the river profiles
    river_df = pd.read_csv(river_csv)
    #remove negative channel slopes
    river_df = river_df[river_df['slope'] > 0]
    print("Got the river csv")

    # set up a figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,6))

    # get the azimuths
    pts = gpd.read_file(DataDirectory+fault_points)
    print("got the fault points")
    # print(pts)
    # calculate azimuth deltas. Decrease in angle = restraining bend, increase = releasing bend
    #pts['az_deltas'] = pts['azimuth'].diff()
    # calculate difference between azimuth and overall plate motion (NW-SE)
    pts['az_deltas'] = pts['azimuth'] - plate_azimuth
    print(pts)

    # first, all the slopes east of the fault (direction < 0)
    east_df = river_df[river_df['direction'] < 0]
    # then all the slopes west of the fault (direction > 0)
    west_df = river_df[river_df['direction'] > 0]
    all_dfs = [east_df, west_df]
    colors = ['r', 'b']
    for i,df in enumerate(all_dfs):
        # get the median slopes
        slope_df = get_median_slope_in_basins(df)

        # now group by the fault dist and plot percentages
        gr = slope_df.groupby(['fault_dist'])['slope_median'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()
        print(gr)
        # round the dataframes so they merge properly
        gr = gr.round({'fault_dist': 3})
        pts = pts.round({'distance': 3})
        new_df = pd.merge(gr, pts, left_on='fault_dist', right_on='distance', how='left')
        print(new_df)
        # print(len(az_deltas))
        # print(len(gr['median']))
        ax.scatter(new_df['az_deltas'], new_df['median'], c=colors[i])

    ax.set_xlim(-5,5)
    ax.invert_xaxis()
    ax.set_xlabel('$\Delta$ strike azimuth ($^\circ$)', fontsize=14)
    ax.set_ylabel('Median channel gradient (m/m)', fontsize=14)
    plt.savefig(DataDirectory+fname_prefix+'_slopes_azimuth_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

def plot_hillslopes_along_fault(hillslope_csv):
    """
    Read in a csv file with the hillslopes and plot compared to distance
    along the fault
    """
    # csv with the river profiles
    df = pd.read_csv(hillslope_csv)
    #remove negative channel slopes
    df = df[df['slope_median'] > 0]

    # set up a figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,8), sharex=True, sharey=True)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

    # plot the channel slope data

    # first, all the slopes east of the fault (direction < 0)
    east_df = df[df['direction'] < 0]
    # then all the slopes west of the fault (direction > 0)
    west_df = df[df['direction'] > 0]

    all_dfs = [east_df, west_df]
    titles = ['North American Plate', 'Pacific Plate']
    colors = ['r', 'b']
    for i, df in enumerate(all_dfs):

        ax[i].grid(color='0.8', linestyle='--', which='both')
        #ax[i].set_ylim(0,0.7)
        ax[i].text(0.04,0.85, titles[i], fontsize=12, transform=ax[i].transAxes, bbox=dict(facecolor='white'))
        #print(gr)

        # now group by the fault dist and plot percentages
        gr = slope_df.groupby(['fault_dist'])['slope_median'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()

        ax[i].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.4', mec='0.4', c='0.4', capsize=2, alpha=0.1)

        # rolling median of channel slopes
        slopes_df = df.sort_values(by='fault_dist')
        slopes_df['slope_rollmedian'] = slopes_df['slope_median'].rolling(50, center=True).median()

        # create a mask for gaps in the median slopes
        these_dists = slopes_df['fault_dist'].values
        mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
        print(mask_starts)
        mc = ma.array(slopes_df['slope_rollmedian'].values)
        mc[mask_starts] = ma.masked
        ax[i].plot(slopes_df['fault_dist'], mc, c=colors[i], zorder=100, lw=3, ls='--')

        # find and plot peaks in the rolling median
        indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.7, min_dist=1000))
        print(indexes)
        peak_dists = slopes_df['fault_dist'].iloc[indexes]
        peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
        print("Hillslope peak distances: ", peak_dists.values)
        ax[i].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
        for j, txt in enumerate(list(peak_dists)):
            ax[i].annotate(str(int(txt))+' km', (list(peak_dists)[j], list(peak_slopes)[j]+0.02), zorder=300)

        #gr.plot.scatter(x='fault_dist', y='median')
    plt.ylabel('Median hillslope gradient (m/m)', labelpad=20)
    #ax[0].set_xlim(100,580)
    #plt.legend(title='Cluster ID', loc='upper right')
    #plt.show()
    #gr.to_csv(DataDirectory+threshold_lvl+fname_prefix+'_channel_slope_fault_dist.csv')

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],0.7), xytext=(labels_dist[i], 0.8), ha='center', fontsize=10, arrowprops=dict(facecolor='k', arrowstyle="->"))

    plt.xlim(100,1100)
    #plt.ylim(0,0.4)
    plt.xlabel('Distance along fault (km)')

    # save the figure
    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_hillslopes_EW_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

def plot_channel_slopes_normalised(DataDirectory, fname_prefix, stream_order, median_river_shp, labels_csv):
    """
    Make a plot of the channel slopes normalised by the median hillslope gradient in each basin
    """

    # csv with the basin data
    river_df = gpd.read_file(median_river_shp)

    # set up a figure
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(18,10), sharex=True, sharey=False)
    #ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

    # get the normalised channel gradient 
    river_df['slope_norm'] = river_df['channel_sl']/river_df['slope_medi']

    # plot the channel slope data

    # first, all the slopes east of the fault (direction < 0)
    east_df = river_df[river_df['direction'] < 0]
    # then all the slopes west of the fault (direction > 0)
    west_df = river_df[river_df['direction'] > 0]

    all_dfs = [east_df, west_df]
    titles = ['North American Plate (east of SAF)', 'Pacific Plate (west of SAF)']
    #colors = ['r', 'b']
    #figlabels = ['a', 'b', 'c', 'd', 'e']
    for i, df in enumerate(all_dfs):

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # CHANNEL GRADIENTS
        ax[0][i].grid(color='0.8', linestyle='--', which='both')
        ax[0][i].set_ylim(0,1)
        ax[0][i].text(0.04,0.85, titles[i], fontsize=14, transform=ax[0][i].transAxes, bbox=dict(facecolor='white'))
        if i == 0:
            ax[0][i].set_ylabel('Median channel\ngradient, $S_c$ (m/m)', labelpad=10, fontsize=16)

        # now group by the fault dist and plot percentages
        gr = df.groupby(['fault_dist'])['channel_sl'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()

        # plot the medians
        ax[0][i].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.4', mec='0.4', c='0.4', capsize=2, alpha=0.1)

        # rolling median of channel slopes
        slopes_df = gr.sort_values(by='fault_dist')
        slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(10, center=True).median()
        print(df.columns)
        # add lat and lon information to this
        #slopes_df = pd.merge(slopes_df, df[['latitude','longitude']], on=['fault_dist'])
        #print(slopes_df)

        # create a mask for gaps in the median slopes
        these_dists = slopes_df['fault_dist'].values
        mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
        #print(mask_starts)
        mc = ma.array(slopes_df['slope_rollmedian'].values)
        mc[mask_starts] = ma.masked
        ax[0][i].plot(slopes_df['fault_dist'], mc, c='b', zorder=100, lw=3, ls='-')

        # find and plot peaks in the rolling median
        indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.5, min_dist=30))
        #print(indexes)
        peak_dists = slopes_df['fault_dist'].iloc[indexes]
        peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
        #peak_dists = these_dists[indexes]
        #peak_slopes = mc[indexes]
        #print(peak_dists)
        #print(peak_slopes)
        #print("Channel slope peak distances: ", peak_dists.values)
        ax[0][i].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
        for j, txt in enumerate(list(peak_dists)):
            ax[0][i].annotate(str(int(txt))+' km', (list(peak_dists)[j]-10, list(peak_slopes)[j]+0.05), zorder=300, fontsize=11)

        # highlight the creeping segment
        ax[0][i].axvspan(400, 580, facecolor='0.5', alpha=0.6)

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # HILLSLOPES
        ax[1][i].grid(color='0.8', linestyle='--', which='both')
        ax[1][i].set_ylim(0,1)
        if i == 0:
            ax[1][i].set_ylabel('Median hillslope\ngradient, $S_h$ (m/m)', labelpad=10, fontsize=16)

        # now group by the fault dist and plot percentages
        gr = df.groupby(['fault_dist'])['slope_medi'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()

        # plot the medians
        ax[1][i].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.4', mec='0.4', c='0.4', capsize=2, alpha=0.1)

        # rolling median of channel slopes
        slopes_df = gr.sort_values(by='fault_dist')
        slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(10, center=True).median()
        print(df.columns)
        # add lat and lon information to this
        #slopes_df = pd.merge(slopes_df, df[['latitude','longitude']], on=['fault_dist'])
        #print(slopes_df)

        # create a mask for gaps in the median slopes
        these_dists = slopes_df['fault_dist'].values
        mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
        #print(mask_starts)
        mc = ma.array(slopes_df['slope_rollmedian'].values)
        mc[mask_starts] = ma.masked
        ax[1][i].plot(slopes_df['fault_dist'], mc, c='r', zorder=100, lw=3, ls='-')

        # find and plot peaks in the rolling median
        indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.5, min_dist=30))
        #print(indexes)
        peak_dists = slopes_df['fault_dist'].iloc[indexes]
        peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
        #peak_dists = these_dists[indexes]
        #peak_slopes = mc[indexes]
        #print(peak_dists)
        #print(peak_slopes)
        #print("Channel slope peak distances: ", peak_dists.values)
        ax[1][i].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
        for j, txt in enumerate(list(peak_dists)):
            ax[1][i].annotate(str(int(txt))+' km', (list(peak_dists)[j]-10, list(peak_slopes)[j]+0.05), zorder=300, fontsize=12)

        # highlight the creeping segment
        ax[1][i].axvspan(400, 580, facecolor='0.5', alpha=0.6)

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # NORMALISED

        ax[2][i].grid(color='0.8', linestyle='--', which='both')
        ax[2][i].set_ylim(0,1)
        if i == 0:
            ax[2][i].set_ylabel('Normalised channel\ngradient ($S_c/S_h$)', labelpad=10, fontsize=16)

        # now group by the fault dist and plot percentages
        gr = df.groupby(['fault_dist'])['slope_norm'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()

        # plot the medians
        ax[2][i].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.4', mec='0.4', c='0.4', capsize=2, alpha=0.1)

        # rolling median of channel slopes
        slopes_df = gr.sort_values(by='fault_dist')
        slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(10, center=True).median()
        print(df.columns)
        # add lat and lon information to this
        #slopes_df = pd.merge(slopes_df, df[['latitude','longitude']], on=['fault_dist'])
        #print(slopes_df)

        # create a mask for gaps in the median slopes
        these_dists = slopes_df['fault_dist'].values
        mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
        #print(mask_starts)
        mc = ma.array(slopes_df['slope_rollmedian'].values)
        mc[mask_starts] = ma.masked
        ax[2][i].plot(slopes_df['fault_dist'], mc, c='purple', zorder=100, lw=3, ls='-')

        # find and plot peaks in the rolling median
        indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.5, min_dist=30))
        #print(indexes)
        peak_dists = slopes_df['fault_dist'].iloc[indexes]
        peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
        #peak_dists = these_dists[indexes]
        #peak_slopes = mc[indexes]
        #print(peak_dists)
        #print(peak_slopes)
        #print("Channel slope peak distances: ", peak_dists.values)
        ax[2][i].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
        for j, txt in enumerate(list(peak_dists)):
            ax[2][i].annotate(str(int(txt))+' km', (list(peak_dists)[j]-10, list(peak_slopes)[j]+0.05), zorder=300, fontsize=10)

        # placenames
        labels_df = pd.read_csv(labels_csv)
        labels = labels_df['Label']
        labels_dist = labels_df['fault_dist']
        for k in range(0, len(labels)):
            ax[0][i].annotate(labels[k], xy=(labels_dist[k],0.99), xytext=(labels_dist[k], 1.1), ha='center', fontsize=12, arrowprops=dict(facecolor='k', arrowstyle="->"))

        # highlight the creeping segment
        ax[2][i].axvspan(400, 580, facecolor='0.5', alpha=0.6)

    plt.xlim(100,1100)
    #plt.ylim(0,0.4)
    plt.xlabel('Distance along fault (km)', fontsize=16, labelpad=10)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.1)

    # save the figure
    plt.savefig(DataDirectory+fname_prefix+'_channels_normalised_EW_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()


def burn_lithology_to_river_df(river_csv, output_csv, lithology_raster):
    """
    read in the csv file with the river profiles and burn the corresponding
    lithology as a new column
    """
    from LSDPlottingTools import LSDMap_GDALIO as IO
    from LSDPlottingTools import LSDMap_PointTools as PT

    # csv with the river profiles
    river_df = pd.read_csv(river_csv)
    #remove negative channel slopes
    river_df = river_df[river_df['slope'] > 0]

    # read in the lithology raster
    this_raster = IO.ReadRasterArrayBlocks(lithology_raster)
    EPSG_string='epsg:32610'
    NDV, xsize, ysize, GeoT, Projection, DataType = IO.GetGeoInfo(lithology_raster)
    CellSize,XMin,XMax,YMin,YMax = IO.GetUTMMaxMin(lithology_raster)

    # get the channel slope data as a point data object
    pts = PT.LSDMap_PointData(river_csv, data_type = 'csv')
    easting, northing = pts.GetUTMEastingNorthing(EPSG_string = EPSG_string)
    lith_values = []

    for x, (i,j) in enumerate(zip(northing, easting)):
	# convert to rows and cols
        X_coordinate_shifted_origin = j - XMin;
        Y_coordinate_shifted_origin = i - YMin;

        col_point = int(X_coordinate_shifted_origin/CellSize);
        row_point = (ysize - 1) - int(round(Y_coordinate_shifted_origin/CellSize));
        # check for data at this cell
        lith_values.append(this_raster[row_point][col_point])

    # now append the lithologies to the river dataframe
    river_df['lithology'] = lith_values

    river_df.to_csv(output_csv, index=False)

def plot_lithology_shapefile(DataDirectory,lithology_shp,fault_shp):
    """
    Make a simple plot of the shapefile with the same colours as the
    slope plots
    """
    from matplotlib.colors import LinearSegmentedColormap

    sf = gpd.read_file(lithology_shp)
    print(sf.columns)

    fault = gpd.read_file(DataDirectory+fault_shp)

    # remove coastal areas (lith_code = 0)
    sf = sf[sf["lith_code"] != 0]

    liths = sf["lith_code"].unique()
    #liths = np.sort(liths)
    sorted_liths = np.sort(liths)
    print(liths)
    print(sorted_liths)
    lith_colors = ['0.25', 'orange', 'red', '#00AD49']
    cmap = LinearSegmentedColormap.from_list('mycmap', lith_colors)
    ax = sf.plot(column="lith_code", cmap=cmap, alpha=0.8, edgecolor='black', lw=0.2)
    fault.plot(ax=ax,color='k')
    #for i, l in enumerate(liths):
     #   print(lith_colors[i])
     #   this_color = ListedColormap(lith_colors[i])
     #   this_df = sf[sf["lith_code"] == l]
     #   this_df.plot(cmap=this_color)
    ax.axis('off')

    plt.savefig(DataDirectory+'CA_lithology.png', dpi=300)
    plt.clf()

def plot_slopes_with_lithology(DataDirectory, fname_prefix, river_csv, labels_csv, stream_order):
    """
    Make a plot of the channel slopes with the lithology overlain

    """

    # csv with the river profiles
    river_df = pd.read_csv(river_csv)
    #remove negative channel slopes
    river_df = river_df[river_df['slope'] > 0]

    # set up a figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,8), sharex=True, sharey=True)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)


    # set up the colours for the lithology plotting
    # index 0 = water/other, index 1 = alluvium, index 2 = sedimentary, index 3 = igneous, index 4 = metamorphic
    lith_colors = ['lightskyblue', 'gray', 'orange', 'red', '#00AD49']

    # plot the channel slope data

    # first, all the slopes east of the fault (direction < 0)
    east_df = river_df[river_df['direction'] < 0]
    # then all the slopes west of the fault (direction > 0)
    west_df = river_df[river_df['direction'] > 0]

    all_dfs = [east_df, west_df]
    titles = ['North American Plate', 'Pacific Plate']
    colors = ['r', 'b']
    for i, df in enumerate(all_dfs):

        ax[i].grid(color='0.8', linestyle='--', which='both')
        ax[i].set_ylim(0,0.7)
        ax[i].set_xlim(100,1066)
        ax[i].text(0.04,0.85, titles[i], fontsize=12, transform=ax[i].transAxes, bbox=dict(facecolor='white'))

        # now group by the fault dist and plot percentages
        gr = slope_df.groupby(['fault_dist'])['slope'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()

        # rolling median of channel slopes
        slopes_df = gr.sort_values(by='fault_dist')
        slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(5, center=True).median()

        # create a mask for gaps in the median slopes
        these_dists = slopes_df['fault_dist'].values
        mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
        print(mask_starts)
        mc = ma.array(slopes_df['slope_rollmedian'].values)
        mc[mask_starts] = ma.masked
        ax[i].plot(slopes_df['fault_dist'], mc, c='k', zorder=100, lw=3, ls='--')

        # find and plot peaks in the rolling median
        indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.35, min_dist=30))
        print(indexes)
        peak_dists = slopes_df['fault_dist'].iloc[indexes]
        peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
        print("PRINTING SLOPE DF")
        print(slopes_df)
        print("Channel slope peak distances: ", peak_dists.values)
        ax[i].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
        for j, txt in enumerate(list(peak_dists)):
            ax[i].annotate(str(int(txt))+' km', (list(peak_dists)[j], list(peak_slopes)[j]+0.05), zorder=300)

        # print the peaks to a text file


        # bin the lithology data by the fault distance
        lith_gr = df.groupby(['fault_dist'])['lithology'].agg(lambda x:x.value_counts().index[0]).reset_index()
        print(lith_gr)
        these_liths = list(lith_gr['lithology'])
        these_dists = list(lith_gr['fault_dist'])
        print(these_liths)

        # get the indices where the lithology changes
        change_idx = np.where(np.roll(these_liths,1)!=these_liths)[0]
        change_dists = lith_gr['fault_dist'].iloc[change_idx]
        print(change_idx)
        print(change_dists)

        # loop through each index where it changes and plot the vspan
        for j, idx in enumerate(change_idx):
            if j == 0:
                min_dist = 0
            else:
                min_dist = these_dists[change_idx[j-1]]
            max_dist = these_dists[idx]
            # lithology of this bar
            if idx == 0:
                color_idx = int(these_liths[0])
            else:
                color_idx = int(these_liths[idx-1])
            # now add the bar
            ax[i].axvspan(min_dist, max_dist, facecolor=lith_colors[color_idx], alpha=0.5)


    # labels for final plot
    plt.ylabel('Median channel slope (m/m)')

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],0.72), xytext=(labels_dist[i], 0.8), ha='center', fontsize=10, arrowprops=dict(facecolor='k', arrowstyle="->"))

    plt.xlim(100,1066)
    #plt.ylim(0,0.4)
    plt.xlabel('Distance along fault (km)')

    # save the figure
    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_slopes_lithology_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

def plot_lithology_deltas(DataDirectory, fname_prefix, river_csv, labels_csv, stream_order):
    """
    Plot delta of channel slopes for each lithology compared to median
    """

    # csv with the river profiles
    river_df = pd.read_csv(river_csv)
    #remove negative channel slopes
    river_df = river_df[river_df['slope'] > 0]

    # set up a figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,8), sharex=True, sharey=True)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)


    # set up the colours for the lithology plotting
    # index 0 = water/other, index 1 = alluvium, index 2 = sedimentary, index 3 = igneous, index 4 = metamorphic
    lith_colors = ['lightskyblue', '0.25', 'orange', 'red', '#00AD49']

    # plot the channel slope data

    # first remove any profiles that drain multiple lithologies
    uniform = river_df.groupby(['id']).lithology.nunique().eq(1)
    river_df = river_df[river_df['id'].isin(uniform[uniform].index)]
    #print(river_df)

    # now, get all the slopes east of the fault (direction < 0)
    east_df = river_df[river_df['direction'] < 0]
    # then all the slopes west of the fault (direction > 0)
    west_df = river_df[river_df['direction'] > 0]

    all_dfs = [east_df, west_df]
    titles = ['North American Plate (east of SAF)', 'Pacific Plate (west of SAF)']
    colors = ['r', 'b']
    for i, df in enumerate(all_dfs):

        ax[i].grid(color='0.8', linestyle='--', which='both')
        #ax[i].set_ylim(0,0.7)
        ax[i].set_xlim(100,1066)
        ax[i].text(0.04,0.85, titles[i], fontsize=14, transform=ax[i].transAxes, bbox=dict(facecolor='white'))

        # get the median line for all lithologies
        gr = slope_df.groupby(['fault_dist'])['slope'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()

        # rolling median of channel slopes
        all_liths = gr.sort_values(by='fault_dist')
        all_liths['slope_rollmedian'] = all_liths['median'].rolling(5).median()


        # Loop through each lithology
        lithologies = df['lithology'].unique()
        for lith_code in lithologies:
            if lith_code != 0:
                # mask for this one lithology
                this_df = df[df['lithology'] == lith_code]
                # now group by the fault dist and plot percentages
                gr = this_df.groupby(['fault_dist'])['slope'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()
                # rolling median of channel slopes
                slopes_df = gr.sort_values(by='fault_dist')
                print(slopes_df)
                slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(5, center=True).median()

                slopes_df['delta'] = slopes_df['slope_rollmedian'] - all_liths['slope_rollmedian']
                #print("DELTA LITH CODE", lith_code)
                #print(delta)

                # create a mask for gaps in the median slopes
                these_dists = slopes_df['fault_dist'].values
                #print(these_dists)
                #print(np.roll(these_dists,1))
                mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
                print(mask_starts)
                mc = ma.array(slopes_df['delta'].values)
                mc[mask_starts] = ma.masked
                #print(slopes_df['slope_rollmedian'])
                ax[i].plot(slopes_df['fault_dist'], mc, c=lith_colors[int(lith_code)], zorder=100, alpha=0.8, lw=3)

                # find and plot peaks in the rolling median
                #indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.35, min_dist=30))

                #peak_dists = slopes_df['fault_dist'].iloc[indexes]
                #peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
                #print("Channel slope peak distances: ", peak_dists.values)
                #ax[i].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
                #for j, txt in enumerate(list(peak_dists)):
                    #ax[i].annotate(str(int(txt))+' km', (list(peak_dists)[j], list(peak_slopes)[j]+0.05), zorder=300)

    # labels for final plot
    plt.ylabel('$\delta$S', fontsize=16, labelpad=10)

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],0.7), xytext=(labels_dist[i], 0.8), ha='center', fontsize=14, arrowprops=dict(facecolor='k', arrowstyle="->"))

    plt.xlim(100,1066)
    #plt.ylim(0,0.4)
    plt.xlabel('Distance along fault (km)', fontsize=16)

    # save the figure
    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_slopes_lith_deltas_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

def plot_channel_slopes_uniform_lithology(DataDirectory, fname_prefix, river_csv, labels_csv, stream_order):
    """
    Read in the channel slope csv file and plot the slopes
    separated by lithology
    """

    # csv with the river profiles
    river_df = pd.read_csv(river_csv)
    #remove negative channel slopes
    river_df = river_df[river_df['slope'] > 0]

    # set up a figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,8), sharex=True, sharey=True)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)


    # set up the colours for the lithology plotting
    # index 0 = water/other, index 1 = alluvium, index 2 = sedimentary, index 3 = igneous, index 4 = metamorphic
    lith_colors = ['lightskyblue', '0.25', 'orange', 'red', '#00AD49']

    # plot the channel slope data

    # first remove any profiles that drain multiple lithologies
    uniform = river_df.groupby(['id']).lithology.nunique().eq(1)
    river_df = river_df[river_df['id'].isin(uniform[uniform].index)]
    #print(river_df)

    # now, get all the slopes east of the fault (direction < 0)
    east_df = river_df[river_df['direction'] < 0]
    # then all the slopes west of the fault (direction > 0)
    west_df = river_df[river_df['direction'] > 0]

    all_dfs = [east_df, west_df]
    titles = ['North American Plate (east of SAF)', 'Pacific Plate (west of SAF)']
    colors = ['r', 'b']
    for i, df in enumerate(all_dfs):

        ax[i].grid(color='0.8', linestyle='--', which='both')
        ax[i].set_ylim(0,0.7)
        ax[i].set_xlim(100,1066)
        ax[i].text(0.04,0.85, titles[i], fontsize=14, transform=ax[i].transAxes, bbox=dict(facecolor='white'))

        # Loop through each lithology
        lithologies = df['lithology'].unique()
        for lith_code in lithologies:
            if lith_code != 0:
                # mask for this one lithology
                this_df = df[df['lithology'] == lith_code]
                # now group by the fault dist and plot percentages
                gr = this_df.groupby(['fault_dist'])['slope'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()

                # rolling median of channel slopes
                slopes_df = gr.sort_values(by='fault_dist')
                print(slopes_df)
                slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(5, center=True).median()

                # create a mask for gaps in the median slopes
                these_dists = slopes_df['fault_dist'].values
                #print(these_dists)
                #print(np.roll(these_dists,1))
                mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]
                print(mask_starts)
                mc = ma.array(slopes_df['slope_rollmedian'].values)
                mc[mask_starts] = ma.masked
                #print(slopes_df['slope_rollmedian'])
                #ax[i].plot(slopes_df['fault_dist'], mc, c=lith_colors[int(lith_code)], zorder=100, alpha=0.8, lw=1)
                #if lith_code == 2:
                ax[i].fill_between(slopes_df['fault_dist'], mc, facecolor=lith_colors[int(lith_code)], zorder=1, alpha=0.5)

                # find and plot peaks in the rolling median
                #indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.35, min_dist=30))

                #peak_dists = slopes_df['fault_dist'].iloc[indexes]
                #peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
                #print("Channel slope peak distances: ", peak_dists.values)
                #ax[i].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
                #for j, txt in enumerate(list(peak_dists)):
                    #ax[i].annotate(str(int(txt))+' km', (list(peak_dists)[j], list(peak_slopes)[j]+0.05), zorder=300)

    # labels for final plot
    plt.ylabel('Median channel slope (m/m)', fontsize=16, labelpad=10)

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],0.7), xytext=(labels_dist[i], 0.8), ha='center', fontsize=14, arrowprops=dict(facecolor='k', arrowstyle="->"))

    plt.xlim(100,1066)
    #plt.ylim(0,0.4)
    plt.xlabel('Distance along fault (km)', fontsize=16)

    # save the figure
    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_slopes_uni_lith_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

def plot_uplift_rates_along_fault_slopes(river_csv, uplift_rate_csv, gps_csv, gps_csv_filt):
    """
    Read in a csv file with slip rates along the fault and plot
    compared to distance along the shapefile
    """

    # csv with the river profiles
    river_df = pd.read_csv(river_csv)
    #remove negative channel slopes
    river_df = river_df[river_df['slope'] > 0]
    # csv with the thermochron uplift rates
    sr_df = pd.read_csv(uplift_rate_csv)
    sr_df = sr_df[np.isnan(sr_df['fault_dist']) == False]
    # csv with the gps uplift rates
    gps_df = pd.read_csv(gps_csv)
    # process the gps data to remove points with a record length less than 5 years and uplift rates > 5 or < -5 mm/yr (noise)
    gps_df = process_gps_data(gps_df, 5, 5)
    gps_df.to_csv(gps_csv_filt, index=False)

    # set up a figure
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10,10), sharex=True)

    # plot the channel slope data
    # now group by the fault dist and plot percentages
    gr = river_df.groupby(['fault_dist'])['slope'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()
    #print(gr)
    ax[0].grid(color='0.8', linestyle='--', which='both')
    ax[0].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.3', mec='0.3', c='0.4', capsize=2, alpha=0.1)

    # rolling median of channel slopes
    slopes_df = gr.sort_values(by='fault_dist')
    slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(10, center=True).median()
    ax[0].plot(slopes_df['fault_dist'], slopes_df['slope_rollmedian'], c='r', zorder=100, lw=3, ls='--')

    # find and plot peaks in the rolling median
    indexes = list(peakutils.indexes(slopes_df['slope_rollmedian'], thres=0.5, min_dist=30))
    print(indexes)
    peak_dists = slopes_df['fault_dist'].iloc[indexes]
    peak_slopes = slopes_df['slope_rollmedian'].iloc[indexes]
    print("Channel slope peak distances: ", peak_dists.values)
    ax[0].scatter(peak_dists, peak_slopes, marker="*", c='k', s=100, zorder=200)
    for i, txt in enumerate(list(peak_dists)):
        ax[0].annotate(str(int(txt))+' km', (list(peak_dists)[i], list(peak_slopes)[i]+0.05), zorder=300)

    #gr.plot.scatter(x='fault_dist', y='median')
    ax[0].set_ylabel('Median channel slope (m/m)')
    #ax[0].set_xlim(100,580)
    ax[0].set_ylim(0,0.2)
    #plt.legend(title='Cluster ID', loc='upper right')
    #plt.show()
    gr.to_csv(DataDirectory+threshold_lvl+fname_prefix+'_channel_slope_fault_dist.csv')

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],0.72), xytext=(labels_dist[i], 0.8), ha='center', fontsize=10, arrowprops=dict(facecolor='k', arrowstyle="->"))

    # plot the uplift rate data from thermochron
    aft_df = sr_df[np.isnan(sr_df['AFT(Ma)']) == False]
    ahe_df = sr_df[np.isnan(sr_df['AHe(Ma)']) == False]
    sizes = np.log(1/sr_df['fault_normal_dist']*100)*10
    for i in range(len(sizes)):
        if sizes[i] < 1:
            sizes[i] = 2
    ax[1].grid(color='0.8', linestyle='--', which='both')
    ax[1].scatter(x=aft_df['fault_dist'], y=aft_df['RU(mm/yr)'], s=sizes, marker='D', c= '0.4', edgecolors='0.2', zorder=10, alpha=0.3, label='AFT')
    ax[1].scatter(x=ahe_df['fault_dist'], y=ahe_df['RU(mm/yr)'], s=sizes, marker='s', c= '0.4', edgecolors='0.2', zorder=10, alpha=0.3, label='AHe')

    # gaussian average of uplift rate to get maxima
    sorted_df = sr_df.sort_values(by='fault_dist')
    uplift_rate = sorted_df['RU(mm/yr)'].values
    dist = sorted_df['fault_dist'].values
    # nan checking
    uplift_rate = np.array([x for i, x in enumerate(uplift_rate) if not np.isnan(dist[i])])
    dist = np.array([x for i, x in enumerate(dist) if not np.isnan(dist[i])])

    new_uplift = gaussian_weighted_average(dist, uplift_rate)
    ax[1].fill_between(dist, new_uplift, zorder=5, color='0.5',edgecolor='k', alpha=0.8, lw=2)

    # find and plot peaks in the thermochron uplift rate
    indexes = list(peakutils.indexes(new_uplift, thres=0.01, min_dist=30))
    print(indexes)
    peak_dists = [dist[i] for i in indexes]
    peak_uplifts = [new_uplift[i] for i in indexes]
    print("Thermochron peak distances: ", peak_dists)
    ax[1].scatter(peak_dists, peak_uplifts, marker="*", c='k', s=100, zorder=200)
    for i, txt in enumerate(list(peak_dists)):
        ax[1].annotate(str(int(txt))+' km', (list(peak_dists)[i]-30, list(peak_uplifts)[i]+0.7), zorder=300)

    #ax[1].set_xlabel('Distance along fault (km)')
    ax[1].set_ylabel('Rock uplift rate (mm/yr)')
    ax[1].set_yscale('log')
    ax[1].set_ylim(10**-1.6,10**1.1)
    ax[1].set_title('Thermochronometry')
    ax[1].legend(loc='upper left')

    # plot the gps uplift rate data
    ax[2].grid(color='0.8', linestyle='--', which='both')
    sizes = (np.log(1/gps_df['fault_normal_dist'] * 100)*10).values
    print("SIZES:", sizes)
    for i in range(len(sizes)):
        if sizes[i] < 1:
            sizes[i] = 2
    #ax[2].scatter(x=gps_df['fault_dist'], y=gps_df['RU(mm/yr)'], s=sizes, marker='D', c='0.4', alpha=0.3, edgecolors='0.2', zorder=10)
    ax[2].errorbar(x=gps_df['fault_dist'], y=gps_df['RU(mm/yr)'], yerr=gps_df['RU_uc'], fmt='o',ms=4, marker='D', mfc='0.3', mec='0.3', c='0.4', capsize=2, alpha=0.2)


    # rolling median of gps uplift rates
    sorted_df = gps_df.sort_values(by='fault_dist')
    sorted_df['rollmedian'] = sorted_df['RU(mm/yr)'].rolling(10, center=True).median()
    ax[2].plot(sorted_df['fault_dist'], sorted_df['rollmedian'], c='k', zorder=100, lw=3)

    # add colours for uplift and subsidence
    ax[2].axhspan(0, 5, alpha=0.4, color='red')
    ax[2].axhspan(-5, 0, alpha=0.4, color='deepskyblue')

    ax[2].set_xlabel('Distance along fault (km)')
    ax[2].set_ylabel('GPS vertical rate (mm/yr)')
    #ax[2].set_yscale('log')
    ax[2].set_title('GPS')
    ax[2].set_ylim(-5,5)

    plt.xlim(100,1100)

    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_slopes_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

    # make a plot of the thermochron uplift rate vs. median channel slope
    #print('UPLIFT RATE DIST:', dist)
    #print('UPLIFT RATE MAX:', new_uplift)
    #print('CHAN SLOPE DIST:', slopes_df['fault_dist'])
    #print('CHAN SLOPES:', slopes_df['slope_rollmedian'])

def plot_junction_angles_along_fault(junction_angle_csv, slip_rate_csv, threshold_so=3):
    """
    Make a plot of the junction angle with distance along the fault
    Args:
        threshold_so: remove any stream orders below this from analysis (default=3)
    """
    # csv with the river profiles
    angle_df = pd.read_csv(junction_angle_csv)
    # remove any angles below the threshold so
    angle_df = angle_df[angle_df['junction_stream_order'] > threshold_so]
    # csv with the InSAR slip rates
    sr_df = pd.read_csv(slip_rate_csv)

    # set up a figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,10), sharex=True)

    # plot the juncton angle data
    # now group by the fault dist and plot percentages
    gr = angle_df.groupby(['fault_dist'])['junction_angle'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()
    print(gr)
    ax[0].grid(color='0.8', linestyle='--', which='both')
    ax[0].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.3', mec='0.3', c='0.4', capsize=2, alpha=0.1)

    # rolling mean of channel slopes
    sorted_df = gr.sort_values(by='fault_dist')
    sorted_df['angle_rollmedian'] = sorted_df['median'].rolling(5, center=True).median()
    ax[0].plot(sorted_df['fault_dist'], sorted_df['angle_rollmedian'], c='r', zorder=100, lw=3, ls='--')

    #gr.plot.scatter(x='fault_dist', y='median')
    ax[0].set_ylabel('Median junction angle (deg)')
    #ax[0].set_xlim(100,580)
    #plt.legend(title='Cluster ID', loc='upper right')
    #plt.show()
    gr.to_csv(DataDirectory+threshold_lvl+fname_prefix+'_junction_angles_fault_dist.csv')

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],0.72), xytext=(labels_dist[i], 0.8), ha='center', fontsize=10, arrowprops=dict(facecolor='k', arrowstyle="->"))

    # plot the slip rate data
    ax[1].scatter(sr_df['fault_dist'], sr_df['slip_rate'])
    ax[1].set_xlabel('Distance along fault (km)')
    ax[1].set_ylabel('Slip rate')

    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_JAngles.png', dpi=300)
    plt.clf()

def plot_drainage_density_along_fault(dd_csv, uplift_rate_csv, gps_csv):
    """
    Plot drainage density and uplift along the fault
    """

    # csv with the river profiles
    dd_df = pd.read_csv(dd_csv)
    # csv with the thermochron uplift rates
    sr_df = pd.read_csv(uplift_rate_csv)
    sr_df = sr_df[np.isnan(sr_df['fault_dist']) == False]
    # csv with the gps uplift rates
    gps_df = pd.read_csv(gps_csv)

    # set up a figure
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10,10), sharex=True)

    # plot the channel slope data
    # now group by the fault dist and plot percentages
    #gr = river_df.groupby(['fault_dist'])['slope'].agg(['median', 'std', percentile(25), percentile(75)]).rename(columns={'percentile_25': 'q1', 'percentile_75': 'q2'}).reset_index()
    #print(gr)
    #ax[0].grid(color='0.8', linestyle='--', which='both')
    #ax[0].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.3', mec='0.3', c='0.4', capsize=2, alpha=0.1)

    # rolling mean of channel slopes
    #slopes_df = gr.sort_values(by='fault_dist')
    #slopes_df['slope_rollmedian'] = slopes_df['median'].rolling(10).median()
    #ax[0].plot(slopes_df['fault_dist'], slopes_df['slope_rollmedian'], c='r', zorder=100, lw=3, ls='--')

    # multiply dd values by 10*6 (m/km^2)
    dd_df['drainage_density'] = dd_df['drainage_density'] * 1e6
    # rolling median
    this_df = dd_df.sort_values(by='fault_dist')
    this_df['rollmedian'] = this_df['drainage_density'].rolling(20, center=True).median()
    ax[0].plot(this_df['fault_dist'], this_df['rollmedian'], c='r', zorder=100, lw=2, ls='--')

    ax[0].scatter(dd_df['fault_dist'], dd_df['drainage_density'],s=1, alpha=0.1, c='0.5')
    #gr.plot.scatter(x='fault_dist', y='median')
    ax[0].set_ylabel('Drainage density (m/km$^2$)')
    ax[0].set_ylim(0, 100000)
    #plt.legend(title='Cluster ID', loc='upper right')
    #plt.show()

    # placenames
    labels_df = pd.read_csv(labels_csv)
    labels = labels_df['Label']
    labels_dist = labels_df['fault_dist']
    for i in range(0, len(labels)):
        ax[0].annotate(labels[i], xy=(labels_dist[i],1e5), xytext=(labels_dist[i], 1.2e5), ha='center', fontsize=10, arrowprops=dict(facecolor='k', arrowstyle="->"))

    # plot the uplift rate data from thermochron
    sizes = np.log(1/sr_df['fault_normal_dist']*100)*10
    for i in range(len(sizes)):
        if sizes[i] < 1:
            sizes[i] = 2
    ax[1].grid(color='0.8', linestyle='--', which='both')
    ax[1].scatter(x=sr_df['fault_dist'], y=sr_df['RU(mm/yr)'], s=sizes, marker='D', c= '0.4', edgecolors='k', zorder=10)

    # gaussian average of uplift rate to get maxima
    sorted_df = sr_df.sort_values(by='fault_dist')
    uplift_rate = sorted_df['RU(mm/yr)'].values
    dist = sorted_df['fault_dist'].values
    # nan checking
    uplift_rate = np.array([x for i, x in enumerate(uplift_rate) if not np.isnan(dist[i])])
    dist = np.array([x for i, x in enumerate(dist) if not np.isnan(dist[i])])

    new_uplift = gaussian_weighted_average(dist, uplift_rate)
    ax[1].fill_between(dist, new_uplift, zorder=5, color='0.5',edgecolor='0.5', alpha=0.5)

    #ax[1].set_xlabel('Distance along fault (km)')
    ax[1].set_ylabel('Rock uplift rate (mm/yr)')
    ax[1].set_yscale('log')
    ax[1].set_ylim(10**-1.6,10**1.1)
    ax[1].set_title('Thermochronometry')

    # plot the gps uplift rate data
    ax[2].grid(color='0.8', linestyle='--', which='both')
    sizes = np.log(1/gps_df['fault_normal_dist'] * 100)*10
    for i in range(len(sizes)):
        if sizes[i] < 1:
            sizes[i] = 2
    ax[2].scatter(x=gps_df['fault_dist'], y=gps_df['RU(mm/yr)'], s=sizes, marker='D', c='0.4', edgecolors='k', zorder=10)

    # gaussiain average of uplift rate to get maxima
    sorted_df_gps = gps_df.sort_values(by='fault_dist')
    uplift_rate_gps = sorted_df_gps['RU(mm/yr)'].values
    dist_gps = sorted_df_gps['fault_dist'].values
    new_uplift_gps = gaussian_weighted_average(dist_gps, uplift_rate_gps)

    ax[2].fill_between(dist_gps, new_uplift_gps, zorder=5, color='0.5',edgecolor='0.5', alpha=0.5)

    ax[2].set_xlabel('Distance along fault (km)')
    ax[2].set_ylabel('Uplift rate (mm/yr)')
    ax[2].set_yscale('log')
    ax[2].set_title('GPS')
    #ax[2].set_ylim(0,10**1.1)

    plt.xlim(100,1100)

    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_drainage_density.png', dpi=300)
    plt.clf()

def plot_stream_length_along_fault(river_csv):
    """
    Make a plot of median stream length along the fault to check that
    slopes are not biased by changes in length
    """
    river_df = pd.read_csv(river_csv)
    # for each id get the length
    #ids = river_df.id.unique()
    #print("N channels: ", len(ids))
    gr = river_df.groupby(['id'])
    counts = gr.size().to_frame(name='counts')
    counts = counts.join(gr.agg({'fault_dist': 'mean'}).rename(columns={'fault_dist': 'fault_dist_mean'}))
    counts = counts.join(gr.agg({'drainage_area': 'max'}).rename(columns={'drainage_area': 'drainage_area_max'}))
    counts.reset_index()
    print(counts)
    print('Mean n pixels: ', counts.counts.mean())
    print('Max n pixels: ', counts.counts.max())
    plt.scatter(counts.fault_dist_mean, counts.drainage_area_max/1000000)
    plt.xlabel('Distance along fault (km)')
    plt.ylabel('Max drainage area (km$^2$)')
    plt.savefig(DataDirectory+fname_prefix+'drainage_area_fault_dist.png')
