#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
#   Copyright (C) 2018-2020 University of Dundee. All rights reserved.

#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License along
#   with this program; if not, write to the Free Software Foundation, Inc.,
#   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# ------------------------------------------------------------------------------

"""This script exports ROI intensities for selected images."""


import omero.scripts as scripts
from omero.gateway import BlitzGateway, MapAnnotationWrapper
from omero.rtypes import rlong, rint, rstring, robject, unwrap
from omero.model import RectangleI, EllipseI, LineI, PolygonI, PolylineI, \
    MaskI, LabelI, PointI
from math import sqrt, pi
import re

import ezomero

DEFAULT_FILE_NAME = "Batch_ROI_Export.csv"
INSIGHT_POINT_LIST_RE = re.compile(r'points\[([^\]]+)\]')


def log(data):
    """Handle logging or printing in one place."""
    print(data)



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


def get_metadata_for_embryo(Embryo_id,
                            conn, ids_to_query):
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


def get_export_data(conn, script_params, image, units=None):
    """Get pixel data for shapes on image and returns list of dicts."""
    # log("Image ID %s..." % image.id)

    # Get pixel size in SAME units for all images
    pixel_size_x = None
    pixel_size_y = None
    if units is not None:
        pixel_size_x = image.getPixelSizeX(units=units)
        pixel_size_x = pixel_size_x.getValue() if pixel_size_x else None
        pixel_size_y = image.getPixelSizeY(units=units)
        pixel_size_y = pixel_size_y.getValue() if pixel_size_y else None

    roi_service = conn.getRoiService()
    all_planes = script_params["Export_All_Planes"]
    include_points = script_params.get("Include_Points_Coords", False)
    size_c = image.getSizeC()
    # Channels index
    channels = script_params.get("Channels", [1])
    ch_indexes = []
    for ch in channels:
        if ch < 1 or ch > size_c:
            log("Channel index: %s out of range 1 - %s" % (ch, size_c))
        else:
            # User input is 1-based
            ch_indexes.append(ch - 1)

    ch_names = image.getChannelLabels()

    ch_names = [ch_name.replace(",", ".") for ch_name in ch_names]
    image_name = image.getName().replace(",", ".")

    result = roi_service.findByImage(image.getId(), None)

    rois = result.rois
    # Sort by ROI.id (same as in iviewer)
    rois.sort(key=lambda r: r.id.val)
    export_data = []

    for roi in rois:
        for shape in roi.copyShapes():
            label = unwrap(shape.getTextValue())
            # wrap label in double quotes in case it contains comma
            label = "" if label is None else '"%s"' % label.replace(",", ".")
            shape_type = shape.__class__.__name__.rstrip('I').lower()
            # If shape has no Z or T, we may go through all planes...
            the_z = unwrap(shape.theZ)
            z_indexes = [the_z]
            if the_z is None and all_planes:
                z_indexes = range(image.getSizeZ())
            # Same for T...
            the_t = unwrap(shape.theT)
            t_indexes = [the_t]
            if the_t is None and all_planes:
                t_indexes = range(image.getSizeT())

            # get pixel intensities
            for z in z_indexes:
                for t in t_indexes:
                    if z is None or t is None:
                        stats = None
                    else:
                        # hopefully this just skips any images where
                        # getShapeStatsRestricted fails
                        try:
                            stats = roi_service.getShapeStatsRestricted(
                                [shape.id.val], z, t, ch_indexes)
                        except:
                            stats = None
                    for c, ch_index in enumerate(ch_indexes):
                        row_data = {
                            "image_id": image.getId(),
                            "image_name": '"%s"' % image_name,
                            "roi_id": roi.id.val,
                            "shape_id": shape.id.val,
                            "type": shape_type,
                            "roi_name": label,
                            "z": z + 1 if z is not None else "",
                            "t": t + 1 if t is not None else "",
                            "channel": ch_names[ch_index],
                            "points": stats[0].pointsCount[c] if stats else "",
                            "min": stats[0].min[c] if stats else "",
                            "max": stats[0].max[c] if stats else "",
                            "sum": stats[0].sum[c] if stats else "",
                            "mean": stats[0].mean[c] if stats else "",
                            "std_dev": stats[0].stdDev[c] if stats else ""
                        }
                        add_shape_coords(shape, row_data,
                                         pixel_size_x, pixel_size_y,
                                         include_points)
                        export_data.append(row_data)
    return export_data


# well_id, well_row, well_column, well_label inserted if SPW
# Points inserted if exporting shape points string
COLUMN_NAMES = ["image_id",
                "image_name",
                "roi_id",
                "shape_id",
                "type",
                "text",
                "z",
                "t",
                "channel",
                "area",
                "length",
                "points",
                "min",
                "max",
                "sum",
                "mean",
                "std_dev",
                "X",
                "Y",
                "Width",
                "Height",
                "RadiusX",
                "RadiusY",
                "X1",
                "Y1",
                "X2",
                "Y2"]


def add_shape_coords(shape, row_data, pixel_size_x, pixel_size_y,
                     include_points=True):
    """Add shape coordinates and length or area to the row_data dict."""
    if shape.getTextValue():
        row_data['Text'] = shape.getTextValue().getValue()
    if isinstance(shape, (RectangleI, EllipseI, PointI, LabelI, MaskI)):
        row_data['X'] = shape.getX().getValue()
        row_data['Y'] = shape.getY().getValue()
    if isinstance(shape, (RectangleI, MaskI)):
        row_data['Width'] = shape.getWidth().getValue()
        row_data['Height'] = shape.getHeight().getValue()
        row_data['area'] = row_data['Width'] * row_data['Height']
    if isinstance(shape, EllipseI):
        row_data['RadiusX'] = shape.getRadiusX().getValue()
        row_data['RadiusY'] = shape.getRadiusY().getValue()
        row_data['area'] = pi * row_data['RadiusX'] * row_data['RadiusY']
    if isinstance(shape, LineI):
        row_data['X1'] = shape.getX1().getValue()
        row_data['X2'] = shape.getX2().getValue()
        row_data['Y1'] = shape.getY1().getValue()
        row_data['Y2'] = shape.getY2().getValue()
        dx = (row_data['X1'] - row_data['X2'])
        dx = dx if pixel_size_x is None else dx * pixel_size_x
        dy = (row_data['Y1'] - row_data['Y2'])
        dy = dy if pixel_size_y is None else dy * pixel_size_y
        row_data['length'] = sqrt((dx * dx) + (dy * dy))
    if isinstance(shape, (PolygonI, PolylineI)):
        point_list = shape.getPoints().getValue()
        match = INSIGHT_POINT_LIST_RE.search(point_list)
        if match is not None:
            point_list = match.group(1)
        if include_points:
            row_data['Points'] = '"%s"' % point_list
    if isinstance(shape, PolylineI):
        coords = point_list.strip(" ").split(" ")
        try:
            coords = [[float(x.strip(", ")) for x in coord.split(",", 1)]
                      for coord in coords]
        except ValueError:
            log("Invalid Polyline coords: %s" % coords)
        else:
            lengths = []
            for i in range(len(coords)-1):
                dx = (coords[i][0] - coords[i + 1][0])
                dy = (coords[i][1] - coords[i + 1][1])
                dx = dx if pixel_size_x is None else dx * pixel_size_x
                dy = dy if pixel_size_y is None else dy * pixel_size_y
                lengths.append(sqrt((dx * dx) + (dy * dy)))
            row_data['length'] = sum(lengths)
    if isinstance(shape, PolygonI):
        # https://www.mathopenref.com/coordpolygonarea.html
        coords = point_list.strip(" ").split(" ")
        try:
            coords = [[float(x.strip(", ")) for x in coord.split(",", 1)]
                      for coord in coords]
        except ValueError:
            log("Invalid Polygon coords: %s" % coords)
        else:
            total = 0
            for c in range(len(coords)):
                coord = coords[c]
                next_c = coords[(c + 1) % len(coords)]
                total += (coord[0] * next_c[1]) - (next_c[0] * coord[1])
            row_data['area'] = abs(0.5 * total)
    if 'area' in row_data and pixel_size_x and pixel_size_y:
        row_data['area'] = row_data['area'] * pixel_size_x * pixel_size_y




def get_csv_header(units_symbol):
    csv_header = ",".join(COLUMN_NAMES)
    if units_symbol is None:
        units_symbol = "pixels"
    csv_header = csv_header.replace(",length,", ",length (%s)," % units_symbol)
    csv_header = csv_header.replace(",area,", ",area (%s)," % units_symbol)
    return csv_header


