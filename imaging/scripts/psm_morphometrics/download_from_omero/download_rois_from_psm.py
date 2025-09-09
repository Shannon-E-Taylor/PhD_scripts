import omero
from omero.gateway import BlitzGateway
import ezomero

from datetime import date
import shutil

import pandas as pd
import numpy as np

from omero_utils import *

# get passwords etc
from config import OMEROUSER, OMEROHOST, OMEROPORT, OMEROWEB

import sys

# for testing.
import getpass
OMEROPASS = getpass.getpass()

from concurrent.futures import ThreadPoolExecutor
# Define the number of threads
num_threads = 8

conn = BlitzGateway(OMEROUSER, OMEROPASS, port=OMEROPORT, host=OMEROHOST)
ret = conn.connect()
print(ret)
assert ret  # Exit if the connection is unsuccessful



today = date.today()


script_params = {
    'Data_Type': 'Image',
    'Include_Points_Coords': False,
    'File_Name': f'../data/archive/blah.csv',
    'IDs': [id],
    'Channels': [1],
    'Export_All_Planes': False
    }




today = date.today()


script_params = {
    'Data_Type': 'Image',
    'Include_Points_Coords': False,
    'File_Name': f'../data/archive/blah.csv',
    'IDs': [id],
    'Channels': [1],
    'Export_All_Planes': False
    }



def extract_measurements_from_psm(conn, image_id):


    label_object, label_pixels = ezomero.get_image(conn,
                                   image_id,
                                   no_pixels=True,
                                   xyzct = True)

    ann = extract_annotations(label_object)

    if 'Embryo_id' in ann.keys():
        Embryo_id = ann['Embryo_id']
    elif 'embryo_id' in ann.keys():
        Embryo_id = ann['embryo_id']
    else:
        print(f'No embryo id in {image_id}!!')
        return

    image_annotations = get_metadata_for_embryo(Embryo_id, conn, all_image_ids)

    # get image
    image = conn.getObject("Image", image_id)

    # extract units for image.
    pixel_size_x = image.getPixelSizeX(units=True)
    try:
        units = pixel_size_x.getUnit()
        units_symbol = pixel_size_x.getSymbol()

        # skip any embryos with no units symbol.
        units_symbol == 'µm'

    except:
        return

    roi_data = get_export_data(conn, script_params,
                               image, units = units)

    if len(roi_data) > 0:
        roi_data = pd.DataFrame(roi_data)
    else:
        return

    for key in image_annotations:
        roi_data[key] = image_annotations[key]

    return roi_data


project_list = [
    9510,  # axial-elongation
    10516, # NMPs
    10463, # phalloidin_stains
    10954  # seg clock stains
]

all_image_ids = []
all_psm_ids = []

for proj_id in project_list:
    id_list = ezomero.get_image_ids(conn, proj_id)
    psm_ids = ezomero.filter_by_kv(conn,
                                   id_list,
                                   key = 'type',
                                   value = 'psm')
    all_psm_ids.extend(psm_ids)
    all_image_ids.extend(id_list)

all_image_ids = np.unique(all_image_ids).tolist()
all_psm_ids = np.unique(all_psm_ids).tolist()
print(len(np.unique(all_image_ids)))
print(len(np.unique(all_psm_ids)))

# Use ThreadPoolExecutor with a specified number of threads
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    # Parallelize using map, which collects results directly
    all_data = list(
        executor.map(
            lambda image_id: extract_measurements_from_psm(
                conn, int(image_id)), all_psm_ids
                ))


misc_data_df = pd.concat(all_data)
misc_data_df.to_csv(f'../data/archive/all_psm_{today}.csv')

shutil.copy(f'../data/archive/all_psm_{today}.csv', '../data/all_psm.csv')

print('psms downloaded!')


conn.close()