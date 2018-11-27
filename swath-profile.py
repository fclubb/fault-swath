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
# shapefiles
from fiona import collection
from shapely.geometry import shape, LineString, mapping
from shapely.geometry import Point as shapelyPoint


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
    temp_dist=0
    # get the spacing based on the total distance and the number of points
    dist = (total_distance/n)
    print("The total distance is", total_distance*1000, ": returning ", n, "points at a spacing of ", dist*1000)
    # have a point at the start of the line
    for j in range(n+1):
        point = line_rvs.interpolate(temp_dist)
        #print(list(point.coords))
        points.append(shapelyPoint(point))
        distances.append(temp_dist*1000)
        temp_dist+=dist

    #output schema
    schema={'geometry': 'Point', 'properties': {'distance': 'float', 'id': 'int'} }

    # write the points to a shapefile
    with collection(DataDirectory+output_shapefile, 'w', crs=crs, driver='ESRI Shapefile', schema=schema) as output:
        for i in range (n+1):
            #print point
            output.write({'properties':{'distance':distances[i], 'id':i },'geometry': mapping(points[i])})

    return points, distances

# get_channel_slope_around_each_point(pts):
#     """
#     Read in a list of shapely points and get a circle with a defined radius around each point

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

DataDirectory='/home/clubb/OneDrive/san_andreas/NorthernSAF/'
#subdirs = [x[0] for x in os.walk(DataDirectory)]
baseline_shapefile='SanAndreasFault.shp'
output_shapefile='SanAndreasPoints.shp'
points, distances = get_points_along_line(n=1024)
coeffs = get_orthogonal_coefficients(points)
cluster_csv = DataDirectory+'tile_70/threshold_2/tile_70_profiles_clustered_SO3.csv'
output_csv=DataDirectory+'tile_70/threshold_2/tile_70_profiles_fault_dist.csv'
bisection_method(points, coeffs, distances, cluster_csv, output_csv)
output_fname=DataDirectory+'tile_70/threshold_2/tile_70_profiles_fault_dist.png'
plot_cluster_stats_along_fault(output_csv)
