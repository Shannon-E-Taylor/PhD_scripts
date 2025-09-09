
import numpy as np
import ezomero
from omero.gateway import BlitzGateway, MapAnnotationWrapper
from omero.rtypes import unwrap
from omero.model import LineI
from ezomero import get_image
import pandas as pd
import os

# for rotating label image
from skimage.transform import resize
from scipy import ndimage

from download import *

# get passwords etc
from config import OMEROUSER, OMEROHOST, OMEROPORT, OMEROWEB
import getpass

# If script has been run from a shell script,
# read the omero password from shell
import sys
if len(sys.argv) > 1:
    OMEROPASS = sys.argv[1]
else:
    OMEROPASS = getpass.getpass()

conn = BlitzGateway(OMEROUSER, OMEROPASS, port=OMEROPORT, host=OMEROHOST)
ret = conn.connect()
assert ret
print(ret)

####
# functions
####

############################
# get list of label images #
############################

# We have two datasets to query:
# 'hcr_labels' has older labels, for methanol-stored data
# 'phalloidin_labels' has labels for non-methanol stored data
# NB I have pre-filtered 'phalloidin_labels' to exclude data from
#        partially imaged tailbuds

hcr_labels = ezomero.get_image_ids(conn, dataset = 28503)
phalloidin_labels = ezomero.get_image_ids(conn, dataset = 31207)
nmp_labels = ezomero.get_image_ids(conn, dataset = 31263)
stellaris_labels = ezomero.get_image_ids(conn, dataset = 31408)
label_dataset_ids = hcr_labels + phalloidin_labels + nmp_labels + stellaris_labels

# remove duplicates
label_dataset_ids = np.unique(label_dataset_ids)

n_labs = len(label_dataset_ids)
print(f'There are {n_labs} total label images to analyse')


###############################################
# get list of previously downloaded label ids #
###############################################

path_to_metadata = '../data/metadata/'

# all_metadata = []

# for f in os.listdir(path_to_metadata):
#     meta = pd.read_csv(f'{path_to_metadata}/{f}')
#     all_metadata.append(meta)

# metadata = pd.concat(all_metadata)
# # label_ids_with_metadata = []
# label_ids_with_metadata = list(metadata['label_id'].dropna().astype(int).values)

# metadata.to_csv('metadata.csv')

# label_ids_to_extract_metadata = [i for i in label_dataset_ids if int(i) not in label_ids_with_metadata]#

# n_to_metadata = len(label_ids_to_extract_metadata)
# print(f'We have {n_to_metadata} images to extract metadata from')

##############################################
# get list of embryos that could have labels #
##############################################

# We need to get a list of all ids that may be associated with our images
# overview/whole embryo images (for somite counts) are probably in 'whole_embryo_scans'
# hcr and dapi stained PSMs could be in those folders
# and we need to scan the entire phalloidin stain project, as
# in this dataset images are organized by imaging session/slide, not by dataset

ax_elongation_proj = ezomero.get_image_ids(conn,
                                     project = 9510)
phalloidin_project = ezomero.get_image_ids(conn, project=10463)
nmp_project = ezomero.get_image_ids(conn, project=10516)
stellaris_project = ezomero.get_image_ids(conn, project=10955)

ids_to_query = ax_elongation_proj + phalloidin_project + nmp_project + stellaris_project
n_images = len(ids_to_query)

ids_to_query = np.unique(ids_to_query)
ids_to_query = [int(i) for i in ids_to_query]

print(f'There are {n_images} images in the whole dataset')



#####################
# Download metadata #
#####################

# for id in label_dataset_ids:
#     print(id)
#     download_metadata_only(int(id), conn, ids_to_query)

# print('finished downloading metadata')


############################
# Download microscopy data #
############################


# exclude previously downloaded data from list.
downloaded_dapi = os.listdir('../data/raw_dapi/')
downloaded_dapi = [int(i.split('.')[0]) for i in downloaded_dapi]

print(f'There are {len(downloaded_dapi)} images already downloaded')

# reload updated metadata
all_metadata = []

for f in os.listdir(path_to_metadata):
    meta = pd.read_csv(f'{path_to_metadata}/{f}')
    all_metadata.append(meta)

metadata = pd.concat(all_metadata)


# Filter metadata to exclude rows where 'psm_id' is in downloaded_dapi
with_dapi = metadata[
    metadata['psm_id'].dropna().astype(int).isin(downloaded_dapi)
    ]


# Convert the 'label_id' column to int for comparison
downloaded_labels = with_dapi['label_id'].dropna().astype(int).tolist()
print(f'There are {len(downloaded_labels)} images already downloaded!')


# Find label_dataset_ids that are not in downloaded_labels
undownloaded_data = [i for i in label_dataset_ids if int(i) not in downloaded_labels]
print(f'There are {len(undownloaded_data)} images to download')

conn.close()

for id in undownloaded_data[::-1]:
    print(id)
    conn = BlitzGateway(OMEROUSER, OMEROPASS, port=OMEROPORT, host=OMEROHOST)
    ret = conn.connect()
    download_data_for_label(int(id), conn, ids_to_query)
    conn.close()