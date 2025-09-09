import numpy as np 
import ezomero
from omero.gateway import BlitzGateway
from ezomero import get_image
from omero.gateway import BlitzGateway

import pandas as pd 

## run me in ezomero_env
# conda activate ezomero_env

# get passwords etc 
from config import OMEROUSER, OMEROHOST, OMEROPORT, OMEROWEB

import sys

OMEROPASS = sys.argv[1]

conn = BlitzGateway(OMEROUSER, OMEROPASS, port=OMEROPORT, host=OMEROHOST)
ret = conn.connect()
assert ret #exit if connection is unsucesful
 
whole_embryo_scan_dataset_id = 26554 # ID to the whole_embryo_scans dataset 
tbox_dataset_id = 26441

def import_metadata_from_dataset(dataset_id, dataset_name): 

    id_list = ezomero.get_image_ids(conn, dataset = dataset_id)

    callipeteras = ezomero.filter_by_kv(conn, id_list, 'species', 'A. calliptera')
    chillis = ezomero.filter_by_kv(conn, id_list, 'species', 'R. chillingali')
    zebras = ezomero.filter_by_kv(conn, id_list, 'species', 'M. zebra')

    # conn.close() # remember to close the connection! 


    list_of_metadata = []

    for id in id_list: 
        im_object, _ = get_image(conn, id, no_pixels=True) # we still want to read metadata in case this has been updated 
        xsize, ysize, zsize = [
            im_object.getPixelSizeX(), im_object.getPixelSizeY(), im_object.getPixelSizeZ()
            ]
        if id in callipeteras: 
            species = 'A. calliptera'
        elif id in chillis: 
            species = 'R. chillingali'
        elif id in zebras: 
            species = 'M. zebra'
        else: 
            species = 'not annotated'


        annotation_id = ezomero.get_map_annotation_ids(conn, 'Image', int(id))
        annotation_dictionary = {}
        for ann_id in annotation_id: 
            ann_dict = ezomero.get_map_annotation(conn, ann_id)
            annotation_dictionary.update(ann_dict)



        if species == 'A. calliptera': 
            annotation_dictionary['max_stage'] = 30 

        elif species == 'R. chillingali': 
            annotation_dictionary['max_stage'] = 38

        elif species == 'M. zebra': 
            annotation_dictionary['max_stage'] = 30 
    
        else: annotation_dictionary['max_stage'] = 0 

        # if it's a late or end stage embryo, we just say it's 30 or 38ss
        if 'stage' in annotation_dictionary.keys(): 
            if annotation_dictionary['stage'] == 'end' or annotation_dictionary['stage'] == 'late': 
                try: 
                    annotation_dictionary['stage'] = annotation_dictionary['max_stage']

                except: 
                    annotation_dictionary['stage'] = 1
                    print(id)

        dictionary = {
            'image_id': id, 'xsize': xsize, 'ysize': ysize, 'zsize': zsize, 'species': species
        }
        annotation_dictionary.update(dictionary)
        list_of_metadata.append(annotation_dictionary)

    metadata = pd.DataFrame(list_of_metadata)
    metadata.to_csv(f'../data/{dataset_name}_metadata.csv')
    print('saved')

import_metadata_from_dataset(whole_embryo_scan_dataset_id, 'whole_embryo_scan')
import_metadata_from_dataset(tbox_dataset_id, 'tbox_hcr')
import_metadata_from_dataset(29659, 'tailbud_dapi')

# np.savetxt(
#     '../calliptera_image_ids.txt', 
#     np.array([str(i) for i in callipeteras]), 
#     delimiter =", ", 
#     fmt ='% s')
# np.savetxt(
#     '../chilli_image_ids.txt', 
#     np.array([str(i) for i in chillis]), 
#     delimiter =", ", 
#     fmt ='% s')



dataset_id = 26441 # ID to tbox stains 
id_list = ezomero.get_image_ids(conn, dataset = dataset_id)
callipeteras = ezomero.filter_by_kv(conn, id_list, 'species', 'A. calliptera')
chillis = ezomero.filter_by_kv(conn, id_list, 'species', 'R. chillingali')
zebras = ezomero.filter_by_kv(conn, id_list, 'species', 'M. zebra')


print(id_list)

for id in id_list: 
    #print(id)
    file_ann_ids = ezomero.get_file_annotation_ids(conn, 'Image', id)
    #print(file_ann_ids)
    if len(file_ann_ids) > 0: # assume that each file only has one annotation 
        #print(file_ann_ids)
        tch_path = ezomero.get_file_annotation(conn, file_ann_ids[0], f'../line_profiles/')
    if id not in callipeteras and id not in chillis and id not in zebras: 
        print(str(id) + ' is not annotated')
    #print(tch_path)

conn.close()

np.savetxt(
    '../data/calliptera_image_ids_tbox.txt', 
    np.array([str(i) for i in callipeteras]), 
    delimiter =", ", 
    fmt ='% s')
np.savetxt(
    '../data/chilli_image_ids_tbox.txt', 
    np.array([str(i) for i in chillis]), 
    delimiter =", ", 
    fmt ='% s')

np.savetxt(
    '../data/zebra_image_ids_tbox.txt', 
    np.array([str(i) for i in zebras]), 
    delimiter =", ", 
    fmt ='% s')