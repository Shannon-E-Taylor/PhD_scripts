import omero
from omero.gateway import BlitzGateway

from datetime import date
import shutil

import ezomero
import pandas as pd

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



def extract_measurements_from_overview(conn, image_id):
    # get image
    image = conn.getObject("Image", image_id)

    # extract units for image.
    pixel_size_x = image.getPixelSizeX(units=True)
    units = pixel_size_x.getUnit()
    units_symbol = pixel_size_x.getSymbol()

    # skip any embryos with no units symbol.
    if units_symbol != 'µm':
        return

    roi_data = get_export_data(conn, script_params,
                               image, units = units)

    if len(roi_data) > 0:
        roi_data = pd.DataFrame(roi_data)
    else:
        return

    # get somite count for embryo
    somite_count = get_somite_count(image_id, conn)
    # omit any point ROIs as these are unnecessary.
    roi_data = roi_data[roi_data['type']!='point']
    roi_data['somite_count'] = somite_count

    # now get image annotations
    image_annotations = extract_annotations(image)

    for key in image_annotations:
        roi_data[key] = image_annotations[key]

    return roi_data


whole_embryo_ids = ezomero.get_image_ids(conn, dataset = 26554)
print(f'There are')
print(len(whole_embryo_ids))
print('whole embryo misc to download')


# Use ThreadPoolExecutor with a specified number of threads
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    # Parallelize using map, which collects results directly
    all_data = list(
        executor.map(
            lambda image_id: extract_measurements_from_overview(
                conn, int(image_id)), whole_embryo_ids
                ))


misc_data_df = pd.concat(all_data)
misc_data_df.to_csv(f'../data/archive/misc_overviews_{today}.csv')

shutil.copy(f'../data/archive/misc_overviews_{today}.csv', '../data/misc_overviews.csv')

print('Phalloidin downloaded!')




def get_overview_ids_from_project(conn, project_id):
    # get all overviews from phalloidin dataset
    all_ids= ezomero.get_image_ids(conn, project = project_id)
    overview_ids = ezomero.filter_by_kv(conn,
                                        all_ids,
                                        key = 'type',
                                        value = 'overview')
    # skip downloaded imageids
    overview_ids = [i for i in overview_ids if i not in whole_embryo_ids]
    return overview_ids

##############
# Download phalloidin data!
##############

phalloidin_ids = get_overview_ids_from_project(conn, 10463)
print(f'There are')
print(len(phalloidin_ids))
print('whole phalloidin stained embryos to download')

# Use ThreadPoolExecutor with a specified number of threads
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    # Parallelize using map, which collects results directly
    all_data = list(
        executor.map(
            lambda image_id: extract_measurements_from_overview(
                conn, int(image_id)), phalloidin_ids
                ))


phal_data_df = pd.concat(all_data)
phal_data_df.to_csv(f'../data/archive/phalloidin_overviews_{today}.csv')

shutil.copy(f'../data/archive/phalloidin_overviews_{today}.csv', '../data/phalloidin_overviews.csv')

print('Phalloidin downloaded!')

##############
# Download NMP data!
##############

nmp_ids = get_overview_ids_from_project(conn, 10516)
print(f'There are')
print(len(nmp_ids))
print('whole nmp stained embryos to download')


all_data = []
# Use ThreadPoolExecutor with a specified number of threads
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    # Parallelize using map, which collects results directly
    all_data = list(
        executor.map(
            lambda image_id: extract_measurements_from_overview(
                conn, int(image_id)), nmp_ids
                ))


nmp_data_df = pd.concat(all_data)
nmp_data_df.to_csv(f'../data/archive/nmp_overviews_{today}.csv')

shutil.copy(f'../data/archive/nmp_overviews_{today}.csv',
            '../data/nmp_overviews.csv')

print('nmp downloaded!')

##############
# Download segclock data!
##############

seg_clock_ids = get_overview_ids_from_project(conn, 10954)
print(f'There are')
print(len(seg_clock_ids))
print('whole seg clock stained embryos to download')


# Use ThreadPoolExecutor with a specified number of threads
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    # Parallelize using map, which collects results directly
    all_data = list(
        executor.map(
            lambda image_id: extract_measurements_from_overview(
                conn, int(image_id)), seg_clock_ids
                ))


segclock_data_df = pd.concat(all_data)
segclock_data_df.to_csv(f'../data/archive/segclock_overviews_{today}.csv')

shutil.copy(f'../data/archive/segclock_overviews_{today}.csv', '../data/segclock_overviews.csv')

print('segclock downloaded!')




conn.close()

all_data = pd.concat([segclock_data_df,
                      nmp_data_df, phal_data_df,
                      misc_data_df])

all_data.to_csv(f'../data/archive/all_overviews_{today}.csv')

shutil.copy(f'../data/archive/all_overviews_{today}.csv', '../data/all_overviews.csv')
