# dynamic time warping
# FJC 15/06/21
# This script uses dynamic time warping to look at similarity between the geomorphic metrics.

import geopandas as gpd
import pandas as pd
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import math
from dtw import *

def dynamic_time_warping(DataDirectory, fname_prefix, median_basin_df, fault_points, plate_azimuth=135):

    # csv with the basin data
    river_df = gpd.read_file(median_basin_df)
    print(river_df.columns)
    river_df = river_df[river_df['channel_sl'] > 0.01]

    # invert the curvature values
    river_df['ht_curv_me'] = np.abs(river_df['ht_curv_me'])
    river_df = river_df.sort_values(by='fault_dist')

    # first, all the slopes east of the fault (direction < 0)
    east_df = river_df[river_df['direction'] < 0]
    # then all the slopes west of the fault (direction > 0)
    west_df = river_df[river_df['direction'] > 0]

    all_dfs = [east_df, west_df]

    # get the fault azimuths
    #pts = gpd.read_file(DataDirectory+fault_points)

    # calculate azimuth deltas. Increase in angle = restraining bend, decrease = releasing bend
    #az_deltas = pts['azimuth'].diff()
    #pts['az_deltas'] = (pts['azimuth'] - plate_azimuth)*-1
    #print(pts)
    fault_dists = river_df.fault_dist.unique()
    river_data = pd.DataFrame()
    river_data['fault_dist'] = fault_dists
    hillslope_data = []
    hilltop_data = []

    for i,df in enumerate(all_dfs):

        # calculate the rolling median of channel slopes
        gr = df.groupby(['fault_dist'])['channel_sl'].agg(['median']).reset_index()
        slopes_df = gr.sort_values(by='fault_dist')
        slopes_df['river_median'] = slopes_df['median'].rolling(10, center=True).median()

        # merge the 2 dataframes
        #join = slopes_df.merge(pts, left_on='fault_dist', right_on='distance')
        print(slopes_df.columns)
        if i == 0:
            river_data = river_data.merge(slopes_df, on='fault_dist').rename(columns={'river_median': 'river_median_east'}).drop(columns='median')
        else:
            river_data = river_data.merge(slopes_df, on='fault_dist').rename(columns={'river_median': 'river_median_west'}).drop(columns='median')

        # # calculate rolling median of hillslopes
        gr = df.groupby(['fault_dist'])['slope_medi'].agg(['median']).reset_index()
        slopes_df = gr.sort_values(by='fault_dist')
        slopes_df['hillslope_median'] = slopes_df['median'].rolling(10, center=True).median()
        if i == 0:
            river_data = river_data.merge(slopes_df, on='fault_dist').rename(columns={'hillslope_median': 'hillslope_median_east'}).drop(columns='median')
        else:
            river_data = river_data.merge(slopes_df, on='fault_dist').rename(columns={'hillslope_median': 'hillslope_median_west'}).drop(columns='median')


        # # calculate rolling median of hilltop curvature
        # # calculate rolling median of hillslopes
        gr = df.groupby(['fault_dist'])['ht_curv_me'].agg(['median']).reset_index()
        slopes_df = gr.sort_values(by='fault_dist')
        slopes_df['ht_curv_median'] = slopes_df['median'].rolling(10, center=True).median()
        if i == 0:
            river_data = river_data.merge(slopes_df, on='fault_dist').rename(columns={'ht_curv_median': 'ht_curv_median_east'}).drop(columns='median')
        else:
            river_data = river_data.merge(slopes_df, on='fault_dist').rename(columns={'ht_curv_median': 'ht_curv_median_west'}).drop(columns='median')


    print(river_data)
    # print(new_river_data['fault_dist'])

    # loop through each metric and perform DTW
    metrics = ['river_median', 'hillslope_median', 'ht_curv_median']
    # different colours for each metric
    colours = [['#2b93a1', '#1b5c65'],
               ['#EA3546', '#b91324'],
               ['#662E9B', '#411d62']]
    y_labels = ['Median channel\ngradient, $S_c$ (m/m)', 
                'Median hillslope\ngradient, $S_h$ (m/m)',
                'Median hilltop\ncurvature, $C_{ht}$ (m$^{-1}$)']

    fig, ax = plt.subplots(nrows=2, ncols=1)
    ax=ax.ravel()
    ax[0].plot(river_data['fault_dist'], river_data['river_median_east']/river_data['river_median_east'].max(), label='river')
    ax[0].plot(river_data['fault_dist'], river_data['hillslope_median_east']/river_data['hillslope_median_east'].max(), label='hillslope')
    ax[0].plot(river_data['fault_dist'], river_data['ht_curv_median_east']/river_data['ht_curv_median_east'].max(), label='hilltop')

    ax[1].plot(river_data['fault_dist'], river_data['river_median_west']/river_data['river_median_west'].max(), label='river')
    ax[1].plot(river_data['fault_dist'], river_data['hillslope_median_west']/river_data['hillslope_median_west'].max(), label='hillslope')
    ax[1].plot(river_data['fault_dist'], river_data['ht_curv_median_west']/river_data['ht_curv_median_west'].max(), label='hilltop')
    plt.legend(loc='upper left')
    plt.show()

    # for j,m in enumerate(metrics):
    #     river_data = river_data.dropna()
    #     east_river_data = river_data[m+'_east'].values
    #     west_river_data = river_data[m+'_west'].values
    #     print("METRIC: ", m)

    #     # calculate the DTW alignment object
    #     # alignment index 1 = indices of the matched points on the NORTH AMERICAN PLATE (East)
    #     # alignment index 2 = indices of the matched points on the PACIFIC PLATE (West)
    #     alignment = dtw(east_river_data, west_river_data, keep_internals=True, window_type="slantedband", window_args={'window_size': 50})

    #     # make a twoway plot of the alignment
    #     alignment.plot(type="twoway", match_indices=len(east_river_data))
    #     #plt.savefig(DataDirectory+fname_prefix+"{}_alignment.png".format(m), dpi=300)

    #     # get the indices to distances along the fault.
    #     shifted_west_x = np.empty(len(west_river_data))
    #     shifted_west_y = np.empty(len(west_river_data))
    #     river_data['alignment_dist'] = np.nan
    #     river_data['shifted_west_data'] = np.nan
    #     # save the indices of the matching points to the dataframe
    #     for i, idx in enumerate(alignment.index1):
    #         # find the fault distance of the east matching point
    #         east_dist = river_data.iloc[alignment.index1[i]].at['fault_dist']
    #         #print(east_dist)
    #         # find the fault distance of the west matching point
    #         west_dist = river_data.iloc[alignment.index2[i]].at['fault_dist']
    #         # find the channel slopes of the Pacific Plate and shift the index to the matching point on the
    #         # North American Plate
    #         west_gradient = river_data.iloc[alignment.index2[i]].at[m+'_west']
    #         river_data.loc[idx, 'alignment_dist'] = east_dist - west_dist # North American curve - Pacific curve.
    #         river_data.loc[idx, 'shifted_west_data'] = west_gradient

    #         west_data_test = west_river_data[alignment.index2[i]]
    #         shifted_west_y[idx] = west_data_test
    #         shifted_west_x[idx] = east_dist

    #     fig,ax = plt.subplots(nrows=3, ncols=1, figsize=(10,10))
    #     ax = ax.ravel()

    #     # create a mask for gaps in the median slopes
    #     these_dists = river_data['fault_dist'].values
    #     mask_starts = np.where(these_dists-np.roll(these_dists,1) > 10)[0]

    #     # mask and plot the east channel gradients
    #     east_masked = ma.array(river_data[m+'_east'].values)
    #     east_masked[mask_starts] = ma.masked
    #     ax[0].grid(color='0.8', linestyle='--', which='both')
    #     ax[0].plot(river_data['fault_dist'], east_masked, lw=3, c=colours[j][0], label='North American Plate')
    #     #ax[0].set_ylim(0,0.4)
    #     #mask and plot the west channel gradients
    #     west_masked = ma.array(river_data[m+'_west'].values)
    #     west_masked[mask_starts] = ma.masked
    #     ax[0].plot(river_data['fault_dist'], west_masked, lw=3, ls=(0, (3, 1, 1, 1)), c=colours[j][1], label='Pacific Plate')
    #     ax[0].set_ylabel(y_labels[j])
    #     ax[0].legend(loc='upper left')
    #     ax[0].set_title('Original data along strike')

    #     # plot the aligned curves
    #     shifted_west_masked = ma.array(shifted_west_y)
    #     shifted_west_masked[mask_starts] = ma.masked
    #     ax[1].grid(color='0.8', linestyle='--', which='both')
    #     #ax[1].set_ylim(0,0.4)
    #     ax[1].plot(river_data['fault_dist'], east_masked, lw=3, c=colours[j][0], label='North American Plate (fixed)')
    #     ax[1].plot(shifted_west_x, shifted_west_masked, lw=3, ls=(0, (3, 1, 1, 1)), c=colours[j][1], label='Pacific Plate (horizontal shift)')
    #     ax[1].set_ylabel(y_labels[j])
    #     ax[1].legend(loc='upper left')
    #     ax[1].set_title('Data after dynamic time warping (North American Plate fixed)')

    #     # mask and plot the offsets
    #     offset_masked = ma.array(river_data['alignment_dist'].values)
    #     offset_masked[mask_starts] = ma.masked
    #     ax[2].grid(color='0.8', linestyle='--', which='both')
    #     ax[2].plot(river_data['fault_dist'], offset_masked, lw=2, c='black')
    #     plt.xlabel('Distance along fault (km)')
    #     ax[2].set_ylabel('Horizontal offset (km)')
    #     ax[2].set_title('Local horizontal offset between North American Plate and Pacific Plate')

    #     plt.subplots_adjust(hspace=0.2)
    #     #plt.savefig(DataDirectory+fname_prefix+"_{}_dtw.png".format(m), transparent=False, dpi=300)
    #     #plt.show()





    