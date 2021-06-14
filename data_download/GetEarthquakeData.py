#!/usr/bin/env python
# coding: utf-8

# # San Andreas fault seismic analysis
#
# This is a script to download earthquake data from USGS ComCat.
# Firstly it downloads the earthquake events within the specific locations bounds and with the magnitude filters.
# Then it gets the focal mechanism for each one, which could be used to work out fault orientation.
#
# Uses the libcomcat python wrappers (https://github.com/usgs/libcomcat) and associated jupyter notebooks.

# stdlib imports
from datetime import datetime
from time import time

# Third party imports
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth
import pandas as pd
from IPython.display import display, HTML

# Local imports
from libcomcat.dataframes import get_detail_data_frame
from libcomcat.search import (get_event_by_id, get_authoritative_info, search)

# ## Get the ComCat data
# The first step is to download the general earthquake data from the ComCat database.
# We set the latitude and longitude bounds of the study area and the minimum and maximum magnitude for analysis.

# latitude and longitude bounds
bounds = [-125.2, -115.312, 32.946, 40.84]

# start/end times for analysis
stime = datetime(1900,1,1)
etime = datetime.utcnow()

# magnitude range
minmag = 3
maxmag = 4.5

# We use the `search` function to get a list of all events matching the criteria above, and then the
# `get_detail_data_frame` function. This gives a lot of information about each earthquake including the
# latitude, longitude, depth, magnitude, and focal mechanism information.

# Retrieve list of events
eventlist = search(starttime=stime,
                  endtime=etime,
                  minlatitude=bounds[2],
                  maxlatitude=bounds[3],
                  minlongitude=bounds[0],
                  maxlongitude=bounds[1],
                  minmagnitude=minmag,
                  maxmagnitude=maxmag,
                  eventtype='earthquake')
print("Number of events:", len(eventlist))

# get the detailed dataframe for these events
df = get_detail_data_frame(eventlist)

DataDirectory='/home/bjdd72/TopographicData/TopographicData/san_andreas/Earthquakes/'
df.to_csv(DataDirectory+'California_earthquakes_1900-2020_Mw3_4-5.csv',index=False)
