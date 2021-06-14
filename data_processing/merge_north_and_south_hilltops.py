# script to merge individual csvs into a big one for all tiles

import numpy as np
import pandas as pd
from glob import glob
import os

data_dir = '/raid/fclubb/san_andreas/'
df1 = pd.read_csv(data_dir+'NorthernSAF/all_tiles/NorthernSAF_RidgeData_SO3.csv')
df2 = pd.read_csv(data_dir+'SouthernSAF/all_tiles_merged/SouthernSAF_RidgeData_SO3.csv')

master_df = pd.DataFrame()

master_df=master_df.append(df1, ignore_index=True)

df1_basins = df1.new_id.unique()
#print('N northern basins:', len(df1_basins))
df2_basins = df2.new_id.unique()
#print('N southern basins:', len(df2_basins))
new_ids = np.arange(len(df1_basins), len(df1_basins)+len(df2_basins), step=1)
print(new_ids)
for i, x in enumerate(df2_basins):
    df2.loc[df2.new_id == x, 'new_id'] = new_ids[i]
master_df=master_df.append(df2, ignore_index=True)

print(master_df)
master_df.to_csv(data_dir+'SAF_combined/all_data/'+'SAF_RidgeData_SO3.csv',index=False)

