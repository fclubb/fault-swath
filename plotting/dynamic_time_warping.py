# dynamic time warping
# FJC 15/06/21
# This script uses dynamic time warping to look at similarity between the geomorphic metrics.

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from dtw import *

def dynamic_time_warping(DataDirectory, fname_prefix, median_basin_shp, fault_points, plate_azimuth=135):

    # csv with the basin data
    river_df = gpd.read_file(median_basin_shp)
    print(river_df['fault_dist'])

    # invert the curvature values
    river_df['ht_curv_me'] = np.abs(river_df['ht_curv_me'])


    # get the fault azimuths
    pts = gpd.read_file(DataDirectory+fault_points)

    # calculate azimuth deltas. Increase in angle = restraining bend, decrease = releasing bend
    #az_deltas = pts['azimuth'].diff()
    pts['az_deltas'] = (pts['azimuth'] - plate_azimuth)*-1
    print(pts)

    # calculate the rolling median of channel slopes
    gr = river_df.groupby(['fault_dist'])['channel_sl'].agg(['median']).reset_index()
    slopes_df = gr.sort_values(by='fault_dist')
    slopes_df['river_median'] = slopes_df['median'].rolling(10, center=True).median()

    # merge the 2 dataframes
    join = slopes_df.merge(pts, left_on='fault_dist', right_on='distance')
    print(join)

    # calculate rolling median of hillslopes
    gr = river_df.groupby(['fault_dist'])['slope_medi'].agg(['median']).reset_index()
    slopes_df = gr.sort_values(by='fault_dist')
    slopes_df['hillslope_median'] = slopes_df['median'].rolling(10, center=True).median()
    join = join.merge(slopes_df, left_on='fault_dist', right_on='fault_dist')
    print(join)

    # calculate rolling median of hilltop curvature
    # calculate rolling median of hillslopes
    gr = river_df.groupby(['fault_dist'])['ht_curv_me'].agg(['median']).reset_index()
    slopes_df = gr.sort_values(by='fault_dist')
    slopes_df['ht_curv_median'] = slopes_df['median'].rolling(10, center=True).median()
    join = join.merge(slopes_df, left_on='fault_dist', right_on='fault_dist')
    print(join)      

    #fig, ax = plt.subplots(nrows=1, ncols=1)
    # plt.plot(join['fault_dist'], join['az_deltas'])
    # plt.plot(join['fault_dist'], join['slope_rollmedian'])
    # # ax.invert_yaxis()
    # plt.show()

    # normalise each metric by its maximum value so they can be compared
    deltas = join['az_deltas'].to_numpy() 
    deltas_norm = deltas/np.max(deltas)
    channel_slopes = join['river_median'].to_numpy()
    channel_slopes_norm = channel_slopes/np.nanmax(channel_slopes)
    hillslopes = join['hillslope_median'].to_numpy()
    hillslopes_norm = hillslopes/np.nanmax(hillslopes)
    hilltops = join['ht_curv_median'].to_numpy()
    hilltops_norm = hilltops/np.nanmax(hilltops)

    # remove nans
    deltas_final = deltas_norm[~np.isnan(channel_slopes_norm)]
    channel_slopes_final = channel_slopes_norm[~np.isnan(channel_slopes_norm)]
    hillslopes_final = hillslopes_norm[~np.isnan(channel_slopes_norm)]
    hilltops_final = hilltops_norm[~np.isnan(channel_slopes_norm)]
    #print(len(channel_slopes_final))
    #print(len(deltas_final))
    plt.plot(channel_slopes_final, c='blue')
    plt.plot(hillslopes_final, c='red')
    plt.plot(hilltops_final, c='purple')
    plt.plot(deltas_norm, c='black')
    plt.show()
    
    alignment = dtw(deltas_final, channel_slopes_final, step_pattern=rabinerJuangStepPattern(6, "c"), keep_internals = True)
    print(alignment)
    #alignment.plot(type="threeway")
    alignment.plot(type="twoway", match_indices=len(channel_slopes_final))
    #print(rabinerJuangStepPattern(6,"c"))
    #rabinerJuangStepPattern(6,"c").plot()
    #alignment.dtwPlotAlignment()
    #print(dist)
    #corr_df.to_csv(DataDirectory+fname_prefix+'_pearson_correlation.csv')