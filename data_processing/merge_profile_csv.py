# script to merge individual csvs into a big one for all tiles

import numpy as np
import pandas as pd
from glob import glob
import os

data_dir = '/raid/fclubb/san_andreas/NorthernSAF/'

subdirs = next(os.walk(data_dir))[1]
print(subdirs)

master_df = pd.DataFrame()

# make a new column for the new basin ID
first_id = 0
for d in subdirs:
    if 'tile_' in d:
        fname = data_dir+d+'/'+d+'_profiles_SO3.csv'
        print(fname)
        if os.path.isfile(fname):
            df = pd.read_csv(fname)

            # generate an array for the new basin ids
            basin_ids = df['basin_id'].unique()
            n_ids = len(basin_ids)
            new_ids = np.arange(first_id, first_id+n_ids, step=1)
            print(new_ids)

            # loop through and assign the new ids to the correct rows
            for i, x in enumerate(basin_ids):
                df.loc[df.basin_id == x, 'new_id'] = new_ids[i]
            
            print(df['new_id'])

            # append to master df
            master_df = master_df.append(df, ignore_index= True)
            
            # update the first id
            first_id = first_id+n_ids
            

        else:
            print('File {} does not exist'.format(fname))

print(master_df)
master_df.to_csv(data_dir+'all_tiles/'+'NorthernSAF_profiles_SO3.csv',index=False)

