

import numpy as np
import os
import pandas as pd
from skimage.measure import label


volumes = pd.read_csv('../data/volumes.csv')

# Remove data with large voxel sizes.
volumes = volumes[volumes['zsize'] < 1]
volumes = volumes[volumes['xsize'] < 1]

manual_metadata = pd.read_csv('../data/cell_counts_to_keep_v4.csv')
manual_metadata = manual_metadata[manual_metadata['KEEP']]
print(manual_metadata.shape)
cell_numbers = volumes[volumes['psm_id'].isin(manual_metadata['psm_id'])]
print(cell_numbers.shape)

def process(row):
    psm_id = row['psm_id']
    if f'{psm_id}.npy' in os.listdir('../data/spots/'):
        spots = np.load(f'../data/spots/{psm_id}.npy').T
        n_spots = spots.shape[1]
        row['cells_in_psm'] = n_spots
    if f'{psm_id}.npy' in os.listdir('../data/somite_spots/'):
        spots = np.load(f'../data/somite_spots/{psm_id}.npy').T
        n_spots = spots.shape[1]
        n_somites = row['num_segmented_somites']
        row['total_cells_in_somites'] = n_spots
        row['cells_per_somite'] = n_spots / n_somites

    return(row)



cell_numbers = cell_numbers.apply(process, axis = 1)


from datetime import date
today = date.today()

cell_numbers.to_csv(f'../../psm_morphometrics/data/archive/cell_numbers_{today}.csv')
cell_numbers.to_csv(f'../../psm_morphometrics/data/cell_numbers.csv')

print('Done!')