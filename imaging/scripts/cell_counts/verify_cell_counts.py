import numpy as np

import os
import napari

import pandas as pd

import pyclesperanto as cle
cle.select_device('RTX')
cle.get_device()

metadata = pd.read_csv('../data/metadata.csv')

manual_metadata = pd.read_csv('../data/cell_counts_to_keep_v3.csv')
#manual_metadata = manual_metadata[manual_metadata['KEEP']]
print(manual_metadata.shape)
# Skip data that has already been analysed.
# metadata = metadata[~metadata['psm_id'].isin(manual_metadata['psm_id'])]
print(metadata.shape)

metadata[['xsize',  'ysize', 'zsize']] = metadata[['xsize',  'ysize', 'zsize']].astype(float)

# metadata = metadata[metadata['psm_id'].isin([1124777, 1124776, 1124773, 1113239, 1113235,
#                                              1111984, 1111980, 1103439, 1103437, 1103432,
#                                              1103430, 1103429, 1103428, 1103427, 1103426,
#                                              1103424, 1079162, 1066302, 1048645, 1029094,
#                                              1025594, 1025593,

#                                              ])]


# metadata = metadata[metadata['psm_id'].isin([
#     1067832, 1102783, 1108841, 1108843, 1108844
# ])]

print(metadata.shape)
for _, row in metadata.iterrows():
    psm_id = row['psm_id']
    if f'{psm_id}.npy' in os.listdir('../data/spots/'):
        print(_)
        print(psm_id)

        dapi = np.load(f'../data/raw_dapi/{psm_id}.npz')['dapi']
        dapi = np.moveaxis(dapi, [2], [0])
        xsize = row['xsize']
        zsize = row['zsize']

        input_image = cle.push(dapi)

        input_image = cle.scale(input_image,
                                factor_x=xsize*2,
                                factor_y=xsize*2,
                                factor_z=zsize*2,
                                interpolate=True,
                                resize = True
        )

        spots = np.load(f'../data/spots/{psm_id}.npy')
        viewer = napari.Viewer()
        viewer.add_image(cle.pull(input_image))
        viewer.add_points(spots, size = 5, face_color = 'red')
        napari.run()
    else:
        print (f'{psm_id} not found!')

