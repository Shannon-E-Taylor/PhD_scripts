import numpy as np
import pandas as pd

import ezomero
from ezomero import get_image
from omero.gateway import BlitzGateway, MapAnnotationWrapper
from omero.rtypes import unwrap
from omero.model import LineI

# for rotating label image
from skimage.transform import resize
from scipy import ndimage

# all download
from download import *


import math
import os

# get passwords etc
from config import OMEROUSER, OMEROHOST, OMEROPORT, OMEROWEB
import getpass

OMEROPASS = getpass.getpass()

conn = BlitzGateway(OMEROUSER, OMEROPASS, port=OMEROPORT, host=OMEROHOST)
ret = conn.connect()

print(ret)

overwrite = True

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
label_dataset_ids = hcr_labels + phalloidin_labels + nmp_labels

# remove duplicates
label_dataset_ids = np.unique(label_dataset_ids)

n_labs = len(label_dataset_ids)
print(f'There are {n_labs} total label images to analyse')


###############################################
# get list of previously downloaded label ids #
###############################################


path_to_metadata = '../data/metadata/'

all_metadata = []

for f in os.listdir(path_to_metadata):
    meta = pd.read_csv(f'{path_to_metadata}/{f}')
    all_metadata.append(meta)

metadata = pd.concat(all_metadata)
# label_ids_with_metadata = []
label_ids_with_metadata = list(metadata['label_id'].dropna().astype(int).values)

metadata.to_csv('metadata.csv')

if not overwrite:
    label_ids_to_extract_metadata = [i for i in label_dataset_ids if int(i) not in label_ids_with_metadata]#
else:
    label_ids_to_extract_metadata = label_dataset_ids

n_to_metadata = len(label_ids_to_extract_metadata)
print(f'We have {n_to_metadata} images to extract metadata from')

##############################################
# get list of embryos that could have labels #
##############################################

# We need to get a list of all ids that may be associated with our images
# overview/whole embryo images (for somite counts) are probably in 'whole_embryo_scans'
# hcr and dapi stained PSMs could be in those folders
# and we need to scan the entire phalloidin stain project, as
# in this dataset images are organized by imaging session/slide, not by dataset

hcr_dataset_list = ezomero.get_image_ids(conn, dataset = 26441)
overviews = ezomero.get_image_ids(conn, dataset = 26554)
dapi_id_list = ezomero.get_image_ids(conn, dataset = 29659)
phalloidin_project = ezomero.get_image_ids(conn, project=10463)
nmp_project = ezomero.get_image_ids(conn, project=10516)

ids_to_query = hcr_dataset_list + overviews + dapi_id_list + phalloidin_project + nmp_project
n_images = len(ids_to_query)

print(f'There are {n_images} images in the whole dataset')

#####################
# Download metadata #
#####################

for idx, id in enumerate(label_ids_to_extract_metadata):
    print(f'Image {idx} of {n_to_metadata}')
    download_metadata_only(int(id), conn, ids_to_query)

print('finished downloading metadata')

conn.close()

all_metadata = []

for f in os.listdir(path_to_metadata):
    meta = pd.read_csv(f'{path_to_metadata}/{f}')
    all_metadata.append(meta)

metadata = pd.concat(all_metadata)
label_ids_with_metadata = []
label_ids_with_metadata = list(metadata['label_id'].dropna().astype(int).values)

metadata.to_csv('metadata.csv')

