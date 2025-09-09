import numpy as np
import pandas as pd

psm_data = pd.read_csv('../data/all_psm.csv')
#psm_data.loc[:, 'roi_name'] = psm_data['roi_name'].str.strip('"')
psm_data.loc[:, 'roi_name'] = psm_data['roi_name'].str.replace('"', '').str.strip()

psm_data.loc[:, 'Z'] = psm_data['z'] * psm_data['zsize']
psm_data.loc[:, 'X'] = psm_data['X'] * psm_data['xsize']
psm_data.loc[:, 'Y'] = psm_data['Y'] * psm_data['ysize']


rois_to_keep = ['noto', 'psm', 'som', 'n']
metadata_to_keep = ['species', 'overview_id', 'somite_count', 'xsize', 'zsize', 'psm_id']


def compute_psm_dimensions(group):

    group = group.reset_index(drop = True)

    if 'w1' in group['roi_name'].values and 'w2' in group['roi_name'].values:
        w1 = np.array(group.loc[group['roi_name'] == 'w1', ['X', 'Y', 'Z']])[0]
        w2 = np.array(group[group['roi_name'] == 'w2'][['X', 'Y', 'Z']])[0]

        psm_width = np.linalg.norm(np.linalg.norm(w1-w2))

    else:
        psm_width = np.nan

    # Filter for rows where 'roi_name' is an empty string
    filtered_rows = group[pd.isna(group['roi_name']) & (group['type'] == 'point')].reset_index(drop=True)

    if filtered_rows.shape[0]==2:
        # Extract the first and second rows if they exist
        d1 = filtered_rows.iloc[0]
        d2 = filtered_rows.iloc[1]

        d1 = d1[['X', 'Y', 'Z']].values
        d2 = d2[['X', 'Y', 'Z']].values

        depth = np.linalg.norm(d1-d2)
    else:
        depth = np.nan

    row = group.loc[0, metadata_to_keep]

    for roi in rois_to_keep:
        value = group[group['roi_name'] == roi]['length'].mean()
        row[roi] = value


    row['psm_depth'] = depth
    row['psm_width'] = psm_width

    return row

manual_meta = pd.read_csv('../../cell_counts/data/manual_metadata.csv')
manual_meta['somite_count'] = manual_meta['stage']

psm_data.reset_index(drop = True, inplace = True)

# compute the psm data.
computed_psm_data = psm_data.groupby('psm_id').apply(compute_psm_dimensions)
computed_psm_data.reset_index(drop = True, inplace = True)

# Step 1 - merge dataframes
merged_df = pd.merge(computed_psm_data, manual_meta[['psm_id', 'somite_count', 'species']], on='psm_id', how='left')

# Step 2: Overwrite the columns 'stage' and 'species' in 'metadata' with those from 'manual_meta'
merged_df['species'] = merged_df['species_y'].combine_first(
    merged_df['species_x'])
merged_df['somite_count'] = merged_df['somite_count_y'].combine_first(merged_df['somite_count_x'])

merged_df.drop(
    ['species_x', 'species_y', 'somite_count_x', 'somite_count_y'],
    axis=1)

# Compute aspect ratios
merged_df['len_width'] = merged_df['psm'] / merged_df['psm_width']
merged_df['len_depth'] = merged_df['psm'] / merged_df['psm_depth']
merged_df['width_depth'] = merged_df['psm_width'] / merged_df['psm_depth']
merged_df['mpz'] = merged_df['psm'] - merged_df['noto']
merged_df['psm_som'] = merged_df['psm'] / merged_df['som']

print(merged_df.shape)

merged_df.to_csv('../data/computed_psm_dimensions.csv')