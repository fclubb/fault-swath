# Create swath profile along an input shapefile and make some plots
# of the cluster along the swath
# FJC 26/11/18

# set backend to run on server
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import os
import time
from matplotlib import rcParams
import matplotlib.cm as cm
import matplotlib.colors as mcolors
# shapefiles
from fiona import collection
from shapely.geometry import shape, LineString, mapping
from shapely.geometry import Point as shapelyPoint
import pyproj as pyproj
from geopy.distance import distance as GeoPyDist

# Set up fonts for plots
label_size = 12
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = label_size
plt.rc('axes', titlesize=10)     # fontsize of the axes title

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
    print(x)
    for i in range(0, len(x)):
        weights= np.exp(-(x-x[i])**2/lenscale)
        #print(weights)
        summe=np.sum(weights*y**power)
        #print(summe)
        new_y[i]=(summe/np.sum(weights))**(1/power)
        #print(new_y[i])

    return new_y


def get_points_along_line(n=1024):
    """
    Interpolate a series of points at equal distances along an input line shapefile. Arguments that need to be supplied are:
    * DataDirectory: the directory of the input/output shapefiles
    * baseline_shapefile: the name of the input line shapefile with extension
    * n: total number of points, should be power2: n = 2^m
    * output_shapefile: the name of the output points shapefile with extension
    """

    points = []
    distances = []
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
    # have a point at the start of the line
    for j in range(n+1):
        point = line_rvs.interpolate(temp_dist)
        if j == 0:
            temp_metric = 0
        else:
            #print(list(point.coords))
            # find the distance between this point and the previous point in metres (vicenty)
            temp_metric = GeoPyDist((point.y, point.x), (points[-1].y, points[-1].x)).km
        metric_dist+=temp_metric
        print(metric_dist)
        temp_dist+=dist
        points.append(shapelyPoint(point))
        distances.append(metric_dist)


    #output schema
    schema={'geometry': 'Point', 'properties': {'distance': 'float', 'id': 'int'} }

    # write the points to a shapefile
    with collection(DataDirectory+output_shapefile, 'w', crs=crs, driver='ESRI Shapefile', schema=schema) as output:
        for i in range (n+1):
            #print point
            output.write({'properties':{'distance':distances[i], 'id':i },'geometry': mapping(points[i])})

    return points, distances

def get_distance_along_fault_from_points(pts_csv, output_pts_csv):
    """
    Find the distance along the fault shapefile from a DataFrame
    with lat/long coordinates
    """
    df = pd.read_csv(pts_csv)
    #df = df.apply(pd.to_numeric, errors='coerce')
    print(df)
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
            #print(line_lat)
            #print(lat[i], line_lat[idx])
            dist = GeoPyDist((line_lat[idx], line_lon[idx]), (line_pt.y, line_pt.x)).km
            #print(dist)
            if line_lat[idx] < lat[i]:
                # line point is at a lower latitude than the point - need to find the distance between them
                # and minus that from the line distance
                dist = line_dist[idx] - dist
            else:
                dist = line_dist[idx] + dist
                print(dist)

            distances.append(dist)
            # write the distance away from the fault for each point
            print(lat[i], lon[i], line_pt.y, line_pt.x, point.y, point.x)
            away_dist = GeoPyDist((line_pt.y, line_pt.x), (point.y, point.x)).km
            print('Dist from fault:', away_dist)
            fault_normal_dists.append(away_dist)
        else:
            distances.append(np.nan)
            fault_normal_dists.append(np.nan)

    # save the distances to csv
    df['fault_dist'] = distances
    df['fault_normal_dist'] = fault_normal_dists
    df.to_csv(output_pts_csv, index=False)

    #plt.scatter(distances, df['slip_rate'])
    #plt.show()

def get_channel_slope_around_each_point(pts, cluster_csv, radius=1000):
    """
    Read in a shapefile of points and get a circle with a defined radius (in metres) around each point,
    then get the median and interquartile range of the channel slope within each circle
    Write to shapefile, yo
    """
    # read in the shapefile and get a list of the lat and lons
    lon = []
    lat = []
    with collection(DataDirectory+pts, 'r') as layer:
        for element in layer:
            lon.append(element['geometry']['coordinates'][0])
            lat.append(element['geometry']['coordinates'][1])

    df = pd.read_csv(cluster_csv)

    # transform the points to projected coordinate system for tracking
    wgs84=pyproj.Proj('+init=EPSG:4326')
    UTM10N=pyproj.Proj('+init=EPSG:32610')
    xs, ys = pyproj.transform(wgs84, UTM10N, lon, lat)

    # now get a buffer around each point
    buffers = [shapelyPoint(x, y).buffer(radius) for x, y in zip(xs, ys)]

    # now get the river points that are in these buffers



def get_orthogonal_coefficients(pts):
    """
    For a list of shapely points, get the line orthogonal to each point and then save
    the parameters for the equation of the line where ax + by + c = 0
    """
    from sympy import Point, Line, N
    # loop through the points between start, stop
    coeffs = []
    print("Getting the coefficients of each line")
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

def bisection_method(points, coeffs, distances, cluster_csv, output_csv):
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
    print(len(points), len(coeffs), len(distances))
    # read in the csv to a pandas dataframe
    df = pd.read_csv(cluster_csv)
    cluster_x = df.longitude.values
    cluster_y = df.latitude.values
    alpha = np.empty(len(cluster_x))

    m = int(math.log(len(points),2))
    print(m)
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
    print(df)

    df.to_csv(output_csv, index=False)

def plot_cluster_stats_along_fault(csv):
    """
    Once you have ran the bisection method then you can read in the csv file
    and plot some statistics of the clustering with distance along the fault
    """
    df = pd.read_csv(csv)

    # now group by the fault dist and plot percentages
    gr = df.groupby(['fault_dist', 'cluster_id'])[['id']].count()
    print(gr)
    plot_df = gr.unstack('cluster_id').loc[:, 'id']
    print(plot_df)
    plot_df.plot()
    plt.xlabel('Distance along fault (km)')
    plt.ylabel('Number of river pixels')
    plt.legend(title='Cluster ID', loc='upper right')
    #plt.show()
    plt.savefig(output_fname, dpi=300)
    plt.clf()

def q1(x):
    return x.quantile(0.25)
def q2(x):
    return x.quantile(0.75)

def plot_channel_slope_along_fault(csv):
    """
    Plot the median channel slope vs. distance along the fault
    """

    df = pd.read_csv(csv)
    #remove negative channel slopes
    df = df[df['slope'] > 0]

    # now group by the fault dist and plot percentages
    f = {'median' : np.median,
         'std' : np.std,
         'q1': q1,
         'q2': q2}
    gr = df.groupby(['fault_dist'])['slope'].agg(f).reset_index()
    print(gr)
    plt.errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=5, marker='D', mfc='r', mec='k', c='0.5', capsize=2)
    #gr.plot.scatter(x='fault_dist', y='median')
    plt.xlabel('Distance along fault (km)')
    plt.ylabel('Median channel slope (m/m)')
    #plt.legend(title='Cluster ID', loc='upper right')
    #plt.show()
    plt.savefig(output_fname, dpi=300)
    plt.clf()

def plot_uplift_rates_along_fault_slopes(river_csv, uplift_rate_csv, gps_csv):
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

    # set up a figure
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10,10), sharex=True)

    # plot the channel slope data
    # now group by the fault dist and plot percentages
    f = {'median' : np.median,
         'std' : np.std,
         'q1': q1,
         'q2': q2}
    gr = river_df.groupby(['fault_dist'])['slope'].agg(f).reset_index()
    #print(gr)
    ax[0].grid(color='0.8', linestyle='--', which='both')
    ax[0].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.3', mec='0.3', c='0.4', capsize=2, alpha=0.1)

    # rolling mean of channel slopes
    sorted_df = gr.sort_values(by='fault_dist')
    sorted_df['slope_rollmedian'] = sorted_df['median'].rolling(10).median()
    ax[0].plot(sorted_df['fault_dist'], sorted_df['slope_rollmedian'], c='r', zorder=100, lw=3, ls='--')

    #gr.plot.scatter(x='fault_dist', y='median')
    ax[0].set_ylabel('Median channel slope (m/m)')
    #ax[0].set_xlim(100,580)
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

    # gaussian average of uplift rate to get maxima
    sorted_df = gps_df.sort_values(by='fault_dist')
    uplift_rate = sorted_df['RU(mm/yr)'].values
    dist = sorted_df['fault_dist'].values
    new_uplift = gaussian_weighted_average(dist, uplift_rate)

    ax[2].fill_between(dist, new_uplift, zorder=5, color='0.5',edgecolor='0.5', alpha=0.5)

    ax[2].set_xlabel('Distance along fault (km)')
    ax[2].set_ylabel('Uplift rate (mm/yr)')
    ax[2].set_yscale('log')
    ax[2].set_title('GPS')
    #ax[2].set_ylim(0,10**1.1)

    plt.xlim(100,1100)

    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_slopes_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

def plot_uplift_rates_along_fault_clusters(river_csv, uplift_rate_csv):
    """
    Read in a csv file with slip rates along the fault and plot
    compared to distance along the shapefile
    """
    # csv with the river profiles
    river_df = pd.read_csv(river_csv)
    #remove negative channel slopes
    #river_df = river_df[river_df['slope'] > 0]
    # csv with the slip rates
    sr_df = pd.read_csv(uplift_rate_csv)

    # set up a figure
    cluster_ids = river_df.cluster_id.unique()
    fig, ax = plt.subplots(nrows=len(cluster_ids)+1, ncols=1, figsize=(6,10), sharex=True)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    color_dict = dict(zip(river_df.cluster_id.unique(), river_df.colour.unique()))
    print(color_dict)
    sorted_colors = [color_dict[key] for key in sorted(color_dict)]
    sorted_colors.append('k')
    print(sorted_colors)


    # now group by the fault dist and plot percentages
    gr = river_df.groupby(['fault_dist', 'cluster_id'])[['id']].count()
    #print(gr)
    # percentages groupby
    gr_pc = gr.groupby(level=0).apply(lambda x: 100 * x / float(x.sum()))
    plot_df = gr_pc.unstack('cluster_id').loc[:, 'id']
    #plot_df = pd.concat([plot_df, sr_df], ignore_index=True)

    print(plot_df.keys)
    # get some data for plotting the clusters
    col_list = sorted(list(river_df.cluster_id.unique()))
    titles = ['Cluster ' + str(col_list[i]) for i in range(len(cluster_ids))]

    # plot the clusters along the fault
    for i in range(len(ax)-1):
        this_col = float(col_list[i])
        ax[i].plot(plot_df[this_col], color=sorted_colors[i])
        ax[i].grid(color='0.8', linestyle='--', which='both')
        if i == 2:
            ax[i].set_ylabel('% channel pixels')
        ax[i].set_title(titles[i], fontsize=10)
        ax[i].set_ylim(0,100)

    # now plot the uplift rate data
    ax[-1].grid(color='0.8', linestyle='--', which='both')
    ax[-1].scatter(x=sr_df['fault_dist'], y=sr_df['RU(mm/yr)'], s=20, marker='D', c= '0.4', edgecolors='k', zorder=10)

    # gaussian average of uplift rate to get maxima
    sorted_df = sr_df.sort_values(by='fault_dist')
    uplift_rate = sorted_df['RU(mm/yr)'].values
    dist = sorted_df['fault_dist'].values
    new_uplift = gaussian_weighted_average(dist, uplift_rate)

    ax[-1].fill_between(dist, new_uplift, zorder=5, color='0.5',edgecolor='0.5', alpha=0.7)
    ax[-1].set_ylabel('Rock uplift rate (mm/yr)')
    ax[-1].set_yscale('log')
    ax[-1].set_ylim(10**-1.6,10**2)

    # axis formatting
    plt.xlabel('Distance along fault (km)')
    plt.subplots_adjust(hspace=0.5)
    #plt.show()
    #plt.xlim(100,1100)

    #save the data
    plot_df.to_csv(DataDirectory+threshold_lvl+fname_prefix+'_fault_dist_clusters_SO{}.csv'.format(stream_order))
    plt.savefig(DataDirectory+threshold_lvl+fname_prefix+'_fault_dist_clusters_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()


def plot_dominant_cluster_along_fault_with_uplift_rate(river_csv, uplift_rate_csv):
    """
    Read in a csv file with uplift rates along the fault and plot
    compared to distance along the shapefile
    """
    # csv with the river profiles
    river_df = pd.read_csv(river_csv)
    #remove negative channel slopes
    #river_df = river_df[river_df['slope'] > 0]
    # csv with the slip rates
    sr_df = pd.read_csv(uplift_rate_csv)

    # set up a figure
    cluster_ids = river_df.cluster_id.unique()
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6,6), sharex=True)
    ax = ax.ravel()

    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')


    # now group by the fault dist and plot percentages
    gr = river_df.groupby(['fault_dist', 'cluster_id'])[['id']].count()
#    # percentages groupby
    gr_pc = gr.groupby(level=0).apply(lambda x: 100 * x / float(x.sum()))
    # max percentage in each group
   # print(gr_pc)
    plot_df = gr_pc.unstack('cluster_id').loc[:, 'id']
    print(plot_df)
    plot_df['Max'] = plot_df.idxmax(axis=1)
    print(plot_df)

    norm=mcolors.Normalize(vmin=1,vmax=8)
    ax[0].grid(color='0.8', linestyle='--', which='both')
    ax[0].scatter(plot_df.index, plot_df['Max'], c=plot_df['Max'], cmap=cm.Set1, norm=norm, s=10, marker='s',zorder=10)
    ax[0].set_ylabel('Dominant cluster ID')
#    #plot_df = pd.concat([plot_df, sr_df], ignore_index=True)
#
#    print(plot_df.keys)
#    # get some data for plotting the clusters
#    col_list = sorted(list(river_df.cluster_id.unique()))
#    titles = ['Cluster ' + str(col_list[i]) for i in range(len(cluster_ids))]
#
#    # plot the clusters along the fault
#    for i in range(len(ax)-1):
#        this_col = float(col_list[i])
#        ax[i].plot(plot_df[this_col], color=sorted_colors[i])
#        ax[i].grid(color='0.8', linestyle='--', which='both')
#        if i == 2:
#            ax[i].set_ylabel('% channel pixels')
#        ax[i].set_title(titles[i], fontsize=10)
#        ax[i].set_ylim(0,100)
#

    # now plot the uplift rate data
    ax[-1].grid(color='0.8', linestyle='--', which='both')
    ax[-1].scatter(x=sr_df['fault_dist'], y=sr_df['RU(mm/yr)'], s=20, marker='D', c= '0.4', edgecolors='k', zorder=10)

    # gaussian average of uplift rate to get maxima
    sorted_df = sr_df.sort_values(by='fault_dist')
    uplift_rate = sorted_df['RU(mm/yr)'].values
    dist = sorted_df['fault_dist'].values
    new_uplift = gaussian_weighted_average(dist, uplift_rate)

    ax[-1].fill_between(dist, new_uplift, zorder=5, color='0.5',edgecolor='0.5', alpha=0.7)
    ax[-1].set_ylabel('Rock uplift rate (mm/yr)')
    ax[-1].set_yscale('log')
    #ax[-1].set_ylim(10**-1.6,10**2)



    # axis formatting
    plt.xlabel('Distance along fault (km)')
#    plt.subplots_adjust(hspace=0.5)
#    #plt.show()
#    #plt.xlim(100,1100)

    #save the data
    plot_df.to_csv(DataDirectory+threshold_lvl+fname_prefix+'_fault_dist_main_cluster_SO{}.csv'.format(stream_order))
    plt.savefig(DataDirectory+threshold_lvl+fname_prefix+'_fault_dist_main_cluster_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

def plot_junction_angles_along_fault(junction_angle_csv, uplift_rate_csv, gps_csv):
    """
    Make a plot of the junction angle with distance along the fault
    """
    # csv with the river profiles
    angle_df = pd.read_csv(junction_angle_csv)
    # csv with the thermochron uplift rates
    sr_df = pd.read_csv(uplift_rate_csv)
    # csv with the gps uplift rates
    gps_df = pd.read_csv(gps_csv)

    # set up a figure
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10,10), sharex=True)

    # plot the juncton angle data
    # now group by the fault dist and plot percentages
    f = {'median' : np.median,
         'std' : np.std,
         'q1': q1,
         'q2': q2}
    gr = angle_df.groupby(['fault_dist'])['junction_angle'].agg(f).reset_index()
    print(gr)
    ax[0].grid(color='0.8', linestyle='--', which='both')
    ax[0].errorbar(x=gr['fault_dist'], y=gr['median'], yerr=[gr['median']-gr['q1'], gr['q2']-gr['median']], fmt='o',ms=4, marker='D', mfc='0.3', mec='0.3', c='0.4', capsize=2, alpha=0.1)

    # rolling mean of channel slopes
    sorted_df = gr.sort_values(by='fault_dist')
    sorted_df['angle_rollmedian'] = sorted_df['median'].rolling(5).median()
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

    # plot the uplift rate data from thermochron
    sizes = np.log(1/sr_df['fault_normal_dist'] * 100)*10
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
    ax[2].scatter(x=gps_df['fault_dist'], y=gps_df['RU(mm/yr)'], s=sizes, marker='D', c='0.4', edgecolors='k', zorder=10)

    # gaussian average of uplift rate to get maxima
    sorted_df = gps_df.sort_values(by='fault_dist')
    uplift_rate = sorted_df['RU(mm/yr)'].values
    dist = sorted_df['fault_dist'].values
    new_uplift = gaussian_weighted_average(dist, uplift_rate)

    ax[2].fill_between(dist, new_uplift, zorder=5, color='0.5',edgecolor='0.5', alpha=0.5)

    ax[2].set_xlabel('Distance along fault (km)')
    ax[2].set_ylabel('Uplift rate (mm/yr)')
    ax[2].set_yscale('log')
    ax[2].set_title('GPS')
    #ax[2].set_ylim(0,10**1.1)

    plt.xlim(100,1100)

    plt.savefig(DataDirectory+fname_prefix+'_fault_dist_JAngles.png', dpi=300)
    plt.clf()


# set the input parameters - will eventually change this to argparse
DataDirectory='/raid/fclubb/san_andreas/SAF_combined/SAF_only/'
threshold_lvl='threshold_2/'
fname_prefix='SAF_only'
stream_order = 3

baseline_shapefile='SanAndreasFault.shp'
output_shapefile='SanAndreasPoints.shp'

#--------------------------------------------------------------------#
# channel profiles

cluster_csv = DataDirectory+threshold_lvl+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order)
output_csv=DataDirectory+threshold_lvl+fname_prefix+'_profiles_fault_dist_S0{}.csv'.format(stream_order)
# check if the fault dist csv already exists
if not os.path.isfile(output_csv):
    points, distances = get_points_along_line(n=512)
    coeffs = get_orthogonal_coefficients(points)
    bisection_method(points, coeffs, distances, cluster_csv, output_csv)

#--------------------------------------------------------------------#
# junction angles

angles_csv = DataDirectory+fname_prefix+'_JAngles.csv'
output_angles_csv = DataDirectory+fname_prefix+'_JAngles_dist.csv'
if not os.path.isfile(output_angles_csv):
    points, distances = get_points_along_line(n=512)
    coeffs = get_orthogonal_coefficients(points)
    bisection_method(points, coeffs, distances, angles_csv, output_angles_csv)

#--------------------------------------------------------------------#
# uplift rates

uplift_rate_csv='/raid/fclubb/san_andreas/Uplift_rates/Spotila_2007.csv'
output_uplift_csv='/raid/fclubb/san_andreas/Uplift_rates/Spotila_2007_dist.csv'
if not os.path.isfile(output_uplift_csv):
    get_distance_along_fault_from_points(uplift_rate_csv, output_uplift_csv)

#--------------------------------------------------------------------#
# gps data

gps_csv='/raid/fclubb/san_andreas/Uplift_rates/gps/MIDAS_IGS08_SAF_50km.csv'
output_gps_csv='/raid/fclubb/san_andreas/Uplift_rates/gps/MIDAS_IGS08_SAF_50km_dist.csv'
if not os.path.isfile(output_gps_csv):
    get_distance_along_fault_from_points(gps_csv, output_gps_csv)

#--------------------------------------------------------------------#
# labels

labels_csv='/raid/fclubb/san_andreas/Uplift_rates/placenames.csv'
get_distance_along_fault_from_points(labels_csv, labels_csv)

plot_uplift_rates_along_fault_slopes(output_csv, output_uplift_csv, output_gps_csv)
#plot_uplift_rates_along_fault_clusters(output_csv, output_uplift_csv)
#plot_dominant_cluster_along_fault_with_uplift_rate(output_csv, output_uplift_csv)
plot_junction_angles_along_fault(output_angles_csv, output_uplift_csv, output_gps_csv)
