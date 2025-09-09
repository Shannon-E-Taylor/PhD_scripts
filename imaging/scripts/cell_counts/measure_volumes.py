import numpy as np
import os
import pandas as pd
from skimage.measure import label


path_to_metadata = '../data/metadata/'

all_metadata = []

fnames = os.listdir('../data/isotropic_psm/')
fnames = [i.split('.')[0] for i in fnames]


for f in os.listdir(path_to_metadata):
    if f.split('.')[0] in fnames:
        meta = pd.read_csv(f'{path_to_metadata}/{f}')
        all_metadata.append(meta)


metadata = pd.concat(all_metadata)

print(metadata.shape)
# only measure volume for good embryos
# metadata = metadata[metadata['quantify_vol']!=False]

metadata.to_csv('../data/metadata.csv')

manual_meta = pd.read_csv('../data/manual_metadata.csv')
manual_meta['somite_count'] = manual_meta['stage']


# Step 1 - merge dataframes
merged_df = pd.merge(metadata,
                     manual_meta[
                         ['psm_id', 'somite_count', 'species']],
                          on='psm_id', how='left')

# Step 2: Overwrite the columns 'stage' and 'species' in 'metadata' with those from 'manual_meta'
merged_df['species'] = merged_df['species_y'].combine_first(
    merged_df['species_x'])
merged_df['somite_count'] = merged_df['somite_count_y'].combine_first(merged_df['somite_count_x'])

merged_df = merged_df.drop(
    ['species_x', 'species_y', 'somite_count_x', 'somite_count_y'],
    axis=1)

# for now, don't measure any embryos without metadata
# merged_df = merged_df.dropna(subset='species')

def measure_volumes(row):
    psm_id = int(row['psm_id'])

    # Load the image
    img = np.load(f'../data/isotropic_psm/{psm_id}.npz')['label']

    # Measure tissue volumes
    psm_volume = np.sum(img[img == 1])  # PSM with label 1
    nt_volume = np.sum(img[img == 2]) / 2  # Neural tube with label 2
    noto_volume = np.sum(img[img == 3]) / 3  # Notocord with label 3

    # Measure somite volumes
    somite_annotation = img.copy()
    somite_annotation[somite_annotation != 4] = 0
    somite_annotation = somite_annotation / 4

    total_volume_of_somites = np.sum(somite_annotation)

    # Correct for the number of somites
    somite_annotation = label(somite_annotation)
    num_segmented_somites = np.max(somite_annotation)
    total_volume_of_somites = total_volume_of_somites / num_segmented_somites if num_segmented_somites != 0 else 0

    # Update row with calculated values
    row['psm_volume'] = psm_volume
    row['neuraltube_volume'] = nt_volume
    row['notocord_volume'] = noto_volume
    row['num_segmented_somites'] = num_segmented_somites
    row['mean_volume_of_somites'] = total_volume_of_somites

    return row


def measure_somite(row):
    psm_id = int(row['psm_id'])

    # just do nothing if we don't have a labelled somite.
    if psm_id not in somite_ids:
        return(row)

    # Load the image
    img = np.load(f'../data/isotropic_somites/{psm_id}.npz')['label']


    somite_annotation = img.copy()

    total_volume_of_somites = np.sum(somite_annotation)
    # Correct for the number of somites
    somite_annotation = label(somite_annotation)
    num_segmented_somites = np.max(somite_annotation)
    total_volume_of_somites = total_volume_of_somites / num_segmented_somites if num_segmented_somites != 0 else 0

    if row['mean_volume_of_somites'] == 0:
        row['mean_volume_of_somites'] = total_volume_of_somites
        row['num_segmented_somites'] = num_segmented_somites

    else:
        print(f'{psm_id} already has somites!!')
    return (row)

somite_ids = os.listdir('../data/isotropic_somites/')
somite_ids = [int(i.split('.')[0]) for i in somite_ids]

# Apply the function to each row in the dataframe
metadata = merged_df.apply(measure_volumes, axis=1)
metadata = metadata.apply(measure_somite, axis=1)


# get today's date
from datetime import date
today = date.today()

# save a dated copy
metadata.to_csv(f'../../psm_morphometrics/data/archive/volumes_{today}.csv')

# save the in-use copy
metadata.to_csv('../../psm_morphometrics/data/volumes.csv')

metadata.to_csv('../data/volumes.csv')

