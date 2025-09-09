import numpy as np
import os
import napari

from skimage.exposure import equalize_hist, rescale_intensity

import pyclesperanto as cle
import pandas as pd
# choose my GPU
cle.select_device('RTX')
cle.get_device()

def old_spot_counting(img, label):

    input_image = cle.push(img)
    input_image = cle.scale(input_image,
                            factor_x=xsize*2,
                            factor_y=xsize*2,
                            factor_z=zsize*2,
                            interpolate=True,
                            resize = True
    )
    label = cle.scale(label,
                            factor_x=xsize*2,
                            factor_y=xsize*2,
                            factor_z=zsize*2,
                            interpolate = False,
                            resize = True)

    blurred = cle.gaussian_blur(input_image,
                                        sigma_x=sigma,  sigma_y=sigma, sigma_z=0
                                        )
    detected_maxima = cle.detect_maxima(blurred,
                                        connectivity='sphere',
                                        radius_x=radius, radius_y=radius, radius_z=radius
                                        )

    # extract PSM after everything is detected.
    detected_maxima = detected_maxima * label
    labeled_maxima = cle.label_spots(detected_maxima)
    labeled_maxima = labeled_maxima * label
    pointlist = cle.labelled_spots_to_pointlist(labeled_maxima)
    detected_spots = np.array([pointlist[2], pointlist[1], pointlist[0]]).T
    spots = detected_spots[~np.isnan(detected_spots).any(axis=1)]

    return(spots, input_image)

def count_spots(img, label):#

    input_image = cle.push(img)


    input_image = cle.scale(input_image,
                            factor_x=xsize*2,
                            factor_y=xsize*2,
                            factor_z=zsize*2,
                            interpolate=True,
                            resize = True
    )
    label = cle.scale(label,
                            factor_x=xsize*2,
                            factor_y=xsize*2,
                            factor_z=zsize*2,
                            interpolate = False,
                            resize = True)

    blurred = cle.gaussian_blur(input_image,
                                sigma_x=sigma,  sigma_y=sigma, sigma_z=0
                                )
    detected_maxima = cle.detect_maxima(blurred,
                                        connectivity='box',
                                        radius_x=radius, radius_y=radius, radius_z=radius
                                        )

    # binary dilate to merge nearby spots
    detected_maxima = cle.binary_dilate(detected_maxima,
                                        radius_x = dilate_radius,
                                        radius_y = dilate_radius,
                                        radius_z = 1,
                                        connectivity='sphere')
    detected_maxima = detected_maxima * label
    labelled_ = cle.label(detected_maxima)
    centroids = cle.reduce_labels_to_centroids(labelled_)
    # labeled_maxima = cle.label_spots(detected_maxima)

    pointlist = cle.labelled_spots_to_pointlist(centroids)
    detected_spots = np.array([pointlist[2], pointlist[1], pointlist[0]]).T
    spots = detected_spots[~np.isnan(detected_spots).any(axis=1)]

    return(spots, input_image)


overwrite = False

metadata = pd.read_csv('../data/metadata.csv')

metadata[['xsize',  'ysize', 'zsize']] = metadata[['xsize',  'ysize', 'zsize']].astype(float)

# Remove data with large voxel sizes.
metadata = metadata[metadata['zsize'] < 1]
metadata = metadata[metadata['xsize'] < 1]

to_skip = [
    1086827,  # too large
    1048642, 1048645, # poor signal in z
]

metadata = metadata[~metadata['psm_id'].isin(to_skip)]
# manual_metadata = pd.read_csv('../data/cell_counts_to_keep.csv')
# keep in overdetected  samples
# manual_metadata = manual_metadata[manual_metadata['Rationale'] != 'overdetection']
# manual_metadata = manual_metadata[manual_metadata['Rationale'] != 'Loss of signal in z']
# metadata = metadata[~metadata['psm_id'].isin(manual_metadata['psm_id'])]

print(metadata.shape)

# We will rescale all images to isotropic (0.5 * 0.5 * 0.5 voxels)
sigma = 1
radius = 3
dilate_radius = 3


for _, row in metadata.iterrows():
    # first extract metadata
    xsize = row['xsize']
    zsize = row['zsize']
    anisotropy = xsize / zsize

    psm_id = row['psm_id']
    print(f'Processing image number {psm_id}')

    if overwrite or f'{psm_id}.npy' not in os.listdir('../data/spots/'):
        try:
            dapi = np.load(f'../data/raw_dapi/{psm_id}.npz')['dapi']
            label = np.load(f'../data/raw_label/{psm_id}.npz')['label']
            # extract the PSM
            label[label!=1] = 0
        except:
            print(f'{psm_id} not found')
            continue
        # switch to ZYX order from XYZ
        dapi = np.moveaxis(dapi, [2], [0])
        label = np.moveaxis(label, [2], [0])

        img = rescale_intensity(dapi)

        spots, input_image = count_spots(img, label)

        print(spots.shape)
        np.save(f'../data/spots/{psm_id}', spots)
        # viewer = napari.Viewer()
        # viewer.add_image(cle.pull(input_image))
        # # viewer.add_image(cle.pull(blurred))
        # viewer.add_points(spots, size = 5, face_color = 'red', out_of_slice_display=True)
        # napari.run()


print('Processing somites...')

# now count cells in somite - for somites annotated within
for _, row in metadata.iterrows():
    # first extract metadata
    xsize = row['xsize']
    zsize = row['zsize']
    anisotropy = xsize / zsize

    psm_id = row['psm_id']
    if f'{psm_id}.npy' not in os.listdir('../data/somite_spots/'):
        if True:
            dapi = np.load(f'../data/raw_dapi/{psm_id}.npz')['dapi']
            label = np.load(f'../data/raw_label/{psm_id}.npz')['label']
            # check if our image has a somite label
            if 4 in label:
                label[label!=4] = 0
                print('found somite')
            else:
                # print(np.max(label))
                continue
        # except:
        #     print(f'{psm_id} not found')
        #     continue
        # switch to ZYX order from XYZ
        dapi = np.moveaxis(dapi, [2], [0])
        label = np.moveaxis(label, [2], [0])

        img = rescale_intensity(dapi)

        spots, input_image = count_spots(img, label)

        print(spots.shape)
        np.save(f'../data/somite_spots/{psm_id}', spots)
        # viewer = napari.Viewer()
        # viewer.add_image(cle.pull(input_image))
        # # viewer.add_image(cle.pull(blurred))
        # viewer.add_points(spots, size = 5, face_color = 'red', out_of_slice_display=True)
        # napari.run()



# Now count cells for images where somites have individual labels
somite_labels = os.listdir('../data/raw_somite_label/')
somite_labels = [int(i.split('.')[0]) for i in somite_labels]

# Filter metadata to exclude rows where 'psm_id' is in downloaded_dapi
somite_metadata = metadata[
    metadata['psm_id'].dropna().astype(int).isin(somite_labels)
    ]

for _, row in somite_metadata.iterrows():
    # first extract metadata
    xsize = row['xsize']
    zsize = row['zsize']
    anisotropy = xsize / zsize

    psm_id = row['psm_id']
    print(f'Processing image number {psm_id}')

    if overwrite or f'{psm_id}.npy' not in os.listdir('../data/somite_spots/'):
        try:
            dapi = np.load(f'../data/raw_dapi/{psm_id}.npz')['dapi']
            label = np.load(f'../data/raw_somite_label/{psm_id}.npz')['label']
            # extract the PSM
            label[label!=1] = 0
        except:
            print(f'{psm_id} not found')
            continue
        # switch to ZYX order from XYZ
        dapi = np.moveaxis(dapi, [2], [0])
        label = np.moveaxis(label, [2], [0])

        img = rescale_intensity(dapi)

        spots, input_image = count_spots(img, label)

        print(spots.shape)
        np.save(f'../data/somite_spots/{psm_id}', spots)




