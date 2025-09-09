# conda activate napari_env
import pandas as pd
import numpy as np

import math

# Image processing libraries
import pyclesperanto as cle
from scipy import ndimage

import napari

# choose my GPU
cle.select_device('RTX')
cle.get_device()

sigmaxy_in_um = 5


def subdivide_psm(img, pixel_sizes, nRegions):
    # img[img > 1] = 0 # delete mpz and psm labels

    # smooth image using a gaussian blur and thresholding
    sigmaxy = sigmaxy_in_um / pixel_sizes[0]
    anisotropy = pixel_sizes[2] / pixel_sizes[0]
    img = cle.gaussian_blur(
        img*1,
        sigma_x = sigmaxy,
        sigma_y = sigmaxy,
        sigma_z = sigmaxy/anisotropy)

    img = cle.threshold_otsu(img)

    img = np.rint(np.array(img))
    # extract indices in the array where the PSM exists
    where_psm = np.where(np.sum(img, axis = (0, 1)) > 0)[0]
    # split these into five regions
    indices = np.round(np.linspace(where_psm[0], where_psm[-1], nRegions+1)).astype(int)

    for val, index in enumerate(indices):
        img[:, :, index:] = img[:, :, index:] * (val + 1)

    return(img)



def get_volume_of_labels(rotated, pixel_sizes):

    xsize, ysize, zsize = pixel_sizes

    img2 = cle.scale(
        rotated,
            factor_z = zsize,
            factor_y = ysize,
            factor_x = xsize,
            interpolate = False,
            resize=True
        )

    img2 = cle.pull(img2)

    volume_dictionary = {}

    for lab in np.unique(img2):
        tmp = img2.copy()
        tmp[tmp!=lab]=0
        tmp[tmp==lab]=1
        volume_dictionary[lab] = np.sum(tmp)

    return(volume_dictionary)


def process(row):

    psm_id = row['psm_id']
    print('psm id: ', psm_id)

    # Compute angle of axis
    src0, src1, dst0, dst1 = row[['src0', 'src1', 'dst0', 'dst1']]
    angle_to_rotate =  math.atan2(
        (dst0 - src0),
        (dst1 - src1)
    )
    angle_to_rotate = 90 - np.rad2deg(angle_to_rotate)
    # if this is nonsense, just skip out of the loop
    if math.isnan(angle_to_rotate):
        return(row)
    print('rotation angle ' , angle_to_rotate)

    # Extract voxel sizes
    xsize = row['xsize']
    ysize = row['ysize']
    zsize = row['zsize']

    print(xsize, zsize)

    # Import annotated PSM,
    # and pre-process such that it is the same dimensions as the spot data
    label_img = np.rint(np.load(f'../data/raw_label/{psm_id}.npz')['label'])
    label_img = np.moveaxis(label_img, [2], [0]).astype(int)

    print(label_img.shape)
    label_img[label_img!=1] = 0
    label_img = cle.push(label_img)
    label_img = cle.scale(label_img,
                          factor_x=xsize*2,
                          factor_y=xsize*2,
                          factor_z=zsize*2,
                          interpolate = False,
                          resize=True)

    rotated_img = ndimage.rotate(label_img, angle_to_rotate, axes = (1, 2))

    # rotated_img = cle.rotate(label_img,
    #                          angle_y = angle_to_rotate,
    #                          interpolate = False, resize = True
    # )
    print('rotated_img', rotated_img.shape)

    # Now run the code to actually split the PSM in to X segments
    subdivided_psm = subdivide_psm(rotated_img, [0.5, 0.5, 0.5], nRegions = 5)
    del rotated_img

    subdivided_psm = cle.push(subdivided_psm.astype(int))
    # measure the volume of the label regions
    volume_dictionary = get_volume_of_labels(subdivided_psm, [0.5, 0.5, 0.5])

    # Now we process the spots
    spots = np.load(f'../data/spots/{psm_id}.npy').T
    spots = spots.astype(int)
    labelled_spots = label_img*0

    labelled_spots[spots[0], spots[1], spots[2]] = 1

    rotated_spots = ndimage.rotate(labelled_spots, angle_to_rotate, axes = (1, 2))

    # rotated_spots = cle.rotate(labelled_spots,
    #                            angle_y = angle_to_rotate,
    #                            interpolate = False,
    #                            resize = True)
    del labelled_spots

    rotated_spots = cle.pull(cle.label(rotated_spots))

    if View:
        viewer = napari.Viewer()
        # viewer.add_labels(label_img)
        viewer.add_labels(subdivided_psm.astype(int))
        viewer.add_labels(rotated_spots)
        napari.run()

    cellcount_dictionary = {}

    subdivided_psm = cle.pull(subdivided_psm)

    if subdivided_psm.shape != rotated_spots.shape:
        print('ERROR')
        # assert True == False

    for lab in [1, 2, 6, 24, 120]:
        tmp = rotated_spots.copy()
        # Delete all spots in regions outside the ROI
        tmp[subdivided_psm!=lab]=0
        # then count the number of unique labels
        cellcount_dictionary[lab] = len(np.unique(tmp)) - 1 # the code will count 0 as a label, so delete it
        # print(lab, len(np.unique(tmp)))

    del rotated_spots
    del subdivided_psm

    for key in [1, 2, 6, 24, 120]:
        print(cellcount_dictionary[key] / volume_dictionary[key])
        row[f'Region_{key}_volume'] = volume_dictionary[key]
        row[f'Region_{key}_cellcount'] = cellcount_dictionary[key]
        row[f'Region_{key}_density'] = cellcount_dictionary[key] / volume_dictionary[key]


    row['total_cell_number'] = spots.shape[1]

    return(row)

# NB I am importing the volumes file, as it will also have psm volume
# (important for density)
# this means that volumes MUST be computed before spot detection can be analysed
metadata = pd.read_csv('../data/volumes.csv')

# Toggle Napari viewer for validation
View = True

metadata[['xsize',  'ysize', 'zsize']] = metadata[['xsize',  'ysize', 'zsize']].astype(float)

# Remove data with large voxel sizes.
metadata = metadata[metadata['zsize'] < 1]
metadata = metadata[metadata['xsize'] < 1]

manual_metadata = pd.read_csv('../data/cell_counts_to_keep_v4.csv')
manual_metadata = manual_metadata[manual_metadata['KEEP']]
print(manual_metadata.shape)
metadata = metadata[metadata['psm_id'].isin(manual_metadata['psm_id'])]

print(metadata.shape)

# toskip = [1048651, 1103379]# causing OOM error

# metadata = metadata[~metadata['psm_id'].isin(toskip)]

# metadata[metadata['psm_id'] == 1103379].apply(process, axis=1)
metadata = metadata.apply(process, axis=1)




from datetime import date
today = date.today()


# metadata.to_csv(f'../../psm_morphometrics/data/archive/Density_across_psm_{today}.csv')
# metadata.to_csv(f'../../psm_morphometrics/data/Density_across_psm.csv')
