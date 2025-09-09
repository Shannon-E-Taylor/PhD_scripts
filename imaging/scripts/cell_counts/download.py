import numpy as np
import ezomero
from omero.gateway import BlitzGateway, MapAnnotationWrapper
from omero.rtypes import unwrap
from omero.model import LineI
from ezomero import get_image

import pandas as pd

import math
import os

# for rotating label image
from skimage.transform import resize
from scipy import ndimage



def extract_annotations(image_object):
    '''
    extracts all key-value pairs from an image.
    '''
    key_value_dict = {}
    # iterate through all annotations
    # and identify any of type MapAnnotation (aka key-value pairs)
    for ann in image_object.listAnnotations():
        if isinstance(ann, MapAnnotationWrapper):
            key_value_pairs = dict(ann.getValue())  # This returns a dictionary of (key, value) pairs
            key_value_dict = {**key_value_dict, **key_value_pairs}
    return(key_value_dict)


def get_somite_count(image_id, conn):
    '''
    This script counts all point ROIs associated with the image
    We assume that all point ROIs denote somites, so the number of somites = number of point rois
    Input: image ID to query, and connection (conn)
    Output: number of somites
    '''
    # extract ROIS from image
    roi_service = conn.getRoiService()
    result = roi_service.findByImage(image_id, None)
    rois = result.rois

    # count point rois
    somite_count = 0

    for roi in rois:
        for shape in roi.copyShapes():
            shape_type = shape.__class__.__name__.rstrip('I').lower()
            if shape_type =='point':
                somite_count+=1

    return(somite_count)


def get_metadata_for_embryo(Embryo_id, conn, ids_to_query):
    '''
    Returns a dictionary with the following metadata:
    psm_id (omero ID for the psm image)
    overview_id (omero ID for the overview)
    embryo_id (embryo name - eg. Attenborough_1)
    species (from overview)
    somite_count(number of point ROIs in overview)
    dapi_index (index of the dapi channel in psm image)

    assumes everything is correctly annotated.
    '''
    # allow both embryo_id and Embryo_id
    matching_ids = ezomero.filter_by_kv(conn, ids_to_query, key = 'Embryo_id', value = Embryo_id)
    matching_ids2 = ezomero.filter_by_kv(conn, ids_to_query, key = 'embryo_id', value = Embryo_id)
    matching_ids = matching_ids + matching_ids2

    print(matching_ids)

    metadata_dictionary = {}

    for id in matching_ids:
        image_object, _ = ezomero.get_image(conn,
                                id,
                                no_pixels=True)

        metadata = extract_annotations(image_object)

        if 'type' in metadata.keys() and metadata['type'] == 'overview':
            somite_count = get_somite_count(id, conn)
            metadata_dictionary['somite_count'] = somite_count
            metadata_dictionary['species'] = metadata['species']
            metadata_dictionary['overview_id'] = id

        if 'type' in metadata.keys() and metadata['type'] == 'psm':
            metadata_dictionary['psm_id'] = id

            # get size of the image in pixels
            metadata_dictionary['x_max'] = image_object.getSizeX()
            metadata_dictionary['y_max'] = image_object.getSizeY()
            metadata_dictionary['z_max'] = image_object.getSizeZ()

            # get size of each pixel in microns
            metadata_dictionary['xsize'] = image_object.getPixelSizeX()
            metadata_dictionary['ysize'] = image_object.getPixelSizeY()
            metadata_dictionary['zsize'] = image_object.getPixelSizeZ()

            # extract channel labels and enforce lowercase
            channel_labels = image_object.getChannelLabels()
            channel_labels = [i.lower() for i in channel_labels]
            if 'dapi' in channel_labels:
                metadata_dictionary['dapi_index'] = channel_labels.index('dapi')

    return(metadata_dictionary)

def extract_ROI_position(conn, id, roi_name):
    '''
    This code takes as input the OMERO id for a image (id),
    the ROI name to search for (roi_name),
    and an OMERO connection,
    and outputs a dictionary containing the ROI start and end points.
    '''

    roi_service = conn.getRoiService()
    result = roi_service.findByImage(id, None)

    success = False

    for roi in result.rois:
        for s in roi.copyShapes():
            if s.getTextValue() and roi_name in s.getTextValue().getValue() and type(s) == LineI:
                roi_info = s
                success = True

    if not success:
        print ('ROI name ' + str(roi_name) + ' is not in ' + str(id))
        return None

    z = unwrap(roi_info.getTheZ())

    roi_information = {}

    roi_information['z'] = z
    roi_information['src'] = [roi_info.getY1().getValue(), roi_info.getX1().getValue()]
    roi_information['dst'] = [roi_info.getY2().getValue(), roi_info.getX2().getValue()]

    return roi_information


def normalize_label_image(label_image,
                          metadata_dictionary,
                          roi_information):
    '''
    Takes the label image, resizes it to 1um/pixel,
    and then rotates along the axis defined by roi_information
    '''

    # extract pixel sizes
    xsize, ysize, zsize = [metadata_dictionary['xsize'],
                           metadata_dictionary['ysize'],
                           metadata_dictionary['zsize']
    ]

    # TODO fix this bodge
    if zsize == None:
        zsize = 1.0

    # resize the label image
    img2 = np.rint(resize(
        label_image, (
            round(label_image.shape[0]*float(xsize)),
            round(label_image.shape[1]*float(ysize)),
            round(label_image.shape[2]*float(zsize))
            ),
            preserve_range=True, anti_aliasing=False
        ))

    # rotate the image
    rads = math.atan2(
        (roi_information['dst'][0] - roi_information['src'][0]),
        (roi_information['dst'][1] - roi_information['src'][1]))

    rotated_img = ndimage.rotate(
        img2, np.rad2deg(rads),
        axes = (0, 1),
        order = 0) # order = 0 should prevent interpolation of labels

    rotated_img = np.rint(rotated_img) # Convert back to integer values
    return (rotated_img)



def download_data_for_label(label_id, conn, ids_to_query):

    label_object, label_pixels = ezomero.get_image(conn,
                                   label_id,
                                   no_pixels=True,
                                   xyzct = True)

    ann = extract_annotations(label_object)

    if 'Embryo_id' in ann.keys():
        Embryo_id = ann['Embryo_id']
    elif 'embryo_id' in ann.keys():
        Embryo_id = ann['embryo_id']
    else:
        print(f'No embryo id in {label_id}!!')
        return

    print(Embryo_id)
    metadata_dictionary = get_metadata_for_embryo(Embryo_id, conn, ids_to_query)

    if 'psm_id' not in metadata_dictionary.keys():
        print(f'Error: {label_id} has no associated psm!')
        return

    roi_data = extract_ROI_position(conn, metadata_dictionary['psm_id'], 'axis')

    # Now download the label pixels
    # This is so the code runs faster if there are errors.
    label_object, label_pixels = ezomero.get_image(conn,
                                   label_id,
                                   no_pixels=False,
                                   xyzct = True)
    label_pixels = np.squeeze(label_pixels)

    normalized_label = normalize_label_image(label_pixels,
                                       metadata_dictionary=metadata_dictionary,
                                       roi_information=roi_data)


    # download dapi channel
    psm_object, dapi_scan = get_image(conn,
                        metadata_dictionary['psm_id'],
                        start_coords=(0, 0, 0,
                                        metadata_dictionary['dapi_index'],
                                        0),
                        axis_lengths=(metadata_dictionary['x_max'],
                                        metadata_dictionary['y_max'],
                                        metadata_dictionary['z_max'],
                                        1,
                                        1),
                        xyzct = True)

    dapi_scan = np.squeeze(dapi_scan)

    psm_id = str(metadata_dictionary['psm_id'])
    metadata_dictionary['embryo_id'] = Embryo_id
    metadata_dictionary['label_id'] = label_id



    #NB all data is saved with the OMERO id for the *PSM* not the label image.
    np.savez_compressed(
        f'../data/isotropic_psm/{psm_id}',
        label = normalized_label)
    np.savez_compressed(
        f'../data/raw_label/{psm_id}',
        label = label_pixels)
    np.savez_compressed(
        f'../data/raw_dapi/{psm_id}',
        dapi = dapi_scan)

    pd.DataFrame(
        metadata_dictionary, index=[0]
        ).to_csv(f'../data/metadata/{psm_id}.csv')
    ## code to read data--
    # dapi = np.load('../data/raw_dapi/1060919.npz')['dapi']


def download_metadata_only(label_id, conn, ids_to_query):
    label_object, label_pixels = ezomero.get_image(conn,
                                   label_id,
                                   no_pixels=True,
                                   xyzct = True)

    ann = extract_annotations(label_object)
    if 'Embryo_id' in ann.keys():
        Embryo_id = ann['Embryo_id']
    elif 'embryo_id' in ann.keys():
        Embryo_id = ann['embryo_id']
    else:
        print(f'No embryo id in {label_id}')
        return
    metadata_dictionary = get_metadata_for_embryo(Embryo_id, conn, ids_to_query)

    metadata_dictionary['x_max_lab'] = label_object.getSizeX()
    metadata_dictionary['y_max_lab'] = label_object.getSizeY()
    metadata_dictionary['z_max_lab'] = label_object.getSizeZ()

    metadata_dictionary = {**metadata_dictionary, **ann}


    if 'psm_id' not in metadata_dictionary.keys():
        print(f'Error: {label_id} has no associated psm!')
        return


    psm_id = str(metadata_dictionary['psm_id'])
    metadata_dictionary['embryo_id'] = Embryo_id
    metadata_dictionary['label_id'] = label_id

    roi_information = extract_ROI_position(conn, metadata_dictionary['psm_id'], 'axis')

    if roi_information:
        metadata_dictionary['z'] = roi_information['z']
        metadata_dictionary['src0'] = roi_information['src'][0]
        metadata_dictionary['src1'] = roi_information['src'][1]
        metadata_dictionary['dst0'] = roi_information['dst'][0]
        metadata_dictionary['dst1'] = roi_information['dst'][1]

    # print(metadata_dictionary)

    pd.DataFrame(
        metadata_dictionary, index=[0]
        ).to_csv(f'../data/metadata/{psm_id}.csv')
