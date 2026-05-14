import numpy as np
import pandas as pd
from skimage.measure import regionprops_table, regionprops, mesh_surface_area, marching_cubes, label
# from skimage.morphology import binary_erosion
from skimage.measure._regionprops import _props_to_dict

from typing import Tuple, Any, Union, List
from pathlib import Path
import itertools
import warnings
import time

# from scipy.ndimage import maximum_position, center_of_mass
from scipy.ndimage import sum as ndi_sum
from scipy.sparse import coo_matrix

from infer_subc.core.img import apply_mask
from infer_subc.utils.batch import list_image_files, find_segmentation_tiff_files
from infer_subc.core.file_io import read_czi_image, read_tiff_image
from infer_subc.quantification.distribution import get_XY_distribution, get_Z_distribution

########################
########################
# Morphology method
#######################
#######################

def _my_props_to_dict(
    rp, label_image, intensity_image=None, properties=("label", "area", "centroid", "bbox"), extra_properties=None, spacing: Union[tuple, None] = None):
    """
    helper for get_summary_stats
    """
    if extra_properties is not None:
        properties = list(properties) + [prop.__name__ for prop in extra_properties]

    if len(rp) == 0:
        ndim = label_image.ndim
        label_image = np.zeros((3,) * ndim, dtype=int)
        label_image[(1,) * ndim] = 1
        if intensity_image is not None:
            intensity_image = np.zeros(label_image.shape + intensity_image.shape[ndim:], dtype=intensity_image.dtype)

        regions = regionprops(label_image, intensity_image=intensity_image, extra_properties=extra_properties, spacing=spacing)

        out_d = _props_to_dict(regions, properties=properties, separator="-")
        return {k: v[:0] for k, v in out_d.items()}

    return _props_to_dict(rp, properties=properties, separator="-")





### USED ### OLDER VERSION OF THE ABOVE
def get_org_morphology_3D(segmentation_img: np.ndarray, 
                           seg_name: str, 
                           intensity_img, 
                           mask: np.ndarray, 
                           scale: Union[tuple, None]=None):
    """
    Parameters
    ------------
    segmentation_img:
        a 3D (ZYX) np.ndarray of segmented objects 
    seg_name: str
        a name or nickname of the object being measured; this will be used for record keeping in the output table
    intensity_img:
        a 3D (ZYX) np.ndarray contain gray scale values from the "raw" image the segmentation is based on )single channel)
    mask:
        a 3D (ZYX) binary np.ndarray mask of the area to measure from
    scale: tuple, optional
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)


    Regionprops measurements included:
    ------------------------
    ['label',
    'centroid',
    'bbox',
    'area',
    'equivalent_diameter',
    'extent',
    'feret_diameter_max',
    'euler_number',
    'convex_area',
    'solidity',
    'axis_major_length',
    'axis_minor_length',
    'max_intensity',
    'mean_intensity',
    'min_intensity']

    Additional measurements included:
    -----------------------
    ['standard_deviation_intensity',
    'surface_area']


    Returns
    -------------
    pandas dataframe of containing regionprops measurements (columns) for each object in the segmentation image (rows) and the regionprops object
    
    """
    ###################################################
    ## MASK THE ORGANELLE OBJECTS THAT WILL BE MEASURED
    ###################################################
    # in case we sent a boolean mask (e.g. cyto, nucleus, cellmask)
    input_labels = _assert_uint16_labels(segmentation_img)

    # mask
    input_labels = apply_mask(input_labels, mask)

    ##########################################
    ## CREATE LIST OF REGIONPROPS MEASUREMENTS
    ##########################################
    # start with LABEL
    properties = ["label"]

    # add position
    properties = properties + ["centroid", "bbox"]

    # add area
    properties = properties + ["area", "equivalent_diameter"] # "num_pixels", 

    # add shape measurements
    properties = properties + ["extent", "euler_number", "solidity", "axis_major_length"] # ,"feret_diameter_max", "axis_minor_length"]

    # add intensity values (used for quality checks)
    properties = properties + ["min_intensity", "max_intensity", "mean_intensity"]

    #######################
    ## ADD EXTRA PROPERTIES
    #######################
    def standard_deviation_intensity(region, intensities):
        return np.std(intensities[region])

    extra_properties = [standard_deviation_intensity]

    ##################
    ## RUN REGIONPROPS
    ##################
    props = regionprops_table(input_labels, 
                           intensity_image=intensity_img, 
                           properties=properties,
                           extra_properties=extra_properties,
                           spacing=scale)

    props_table = pd.DataFrame(props)
    props_table.insert(0, "object", seg_name)
    props_table.rename(columns={"area": "volume"}, inplace=True)

    if scale is not None:
        round_scale = (round(scale[0], 4), round(scale[1], 4), round(scale[2], 4))
        props_table.insert(loc=2, column="scale", value=f"{round_scale}")
    else: 
        props_table.insert(loc=2, column="scale", value=f"{tuple(np.ones(segmentation_img.ndim))}") 

    ##################################################################
    ## RUN SURFACE AREA FUNCTION SEPARATELY AND APPEND THE PROPS_TABLE
    ##################################################################
    surface_area_tab = pd.DataFrame(surface_area_from_props(input_labels, props, scale))

    props_table.insert(12, "surface_area", surface_area_tab)
    props_table.insert(14, "SA_to_volume_ratio", props_table["surface_area"].div(props_table["volume"]))

    ################################################################
    ## ADD SKELETONIZATION OPTION FOR MEASURING LENGTH AND BRANCHING
    ################################################################


    return props_table


# def get_summary_stats_3D(input_labels: np.ndarray, intensity_img, mask: np.ndarray) -> Tuple[Any, Any]:
#     """collect volumentric stats from skimage.measure.regionprops
#         properties = ["label","max_intensity", "mean_intensity", "min_intensity" ,"area"->"volume" , "equivalent_diameter",
#         "centroid", "bbox","euler_number", "extent"
#         +   extra_properties = [standard_deviation_intensity]

#     Parameters
#     ------------
#     input_labels:
#         a 3d  np.ndarray image of the inferred organelle labels
#     intensity_img:
#         a 3d np.ndarray image of the florescence intensity
#     mask:
#         a 3d image containing the cellmask object (mask)

#     Returns
#     -------------
#     pandas dataframe of stats and the regionprops object
#     """

#     # in case we sent a boolean mask (e.g. cyto, nucleus, cellmask)
#     input_labels = _assert_uint16_labels(input_labels)

#     # mask
#     input_labels = apply_mask(input_labels, mask)

#     # start with LABEL
#     properties = ["label"]
#     # add intensity:
#     properties = properties + ["max_intensity", "mean_intensity", "min_intensity"]

#     # arguments must be in the specified order, matching regionprops
#     def standard_deviation_intensity(region, intensities):
#         return np.std(intensities[region])

#     extra_properties = [standard_deviation_intensity]

#     # add area
#     properties = properties + ["area", "equivalent_diameter"]
#     #  position:
#     properties = properties + ["centroid", "bbox"]  # , 'bbox', 'weighted_centroid']
#     # etc
#     properties = properties + ["euler_number", "extent"]  # only works for BIG organelles: 'convex_area','solidity',

#     rp = regionprops(input_labels, intensity_image=intensity_img, extra_properties=extra_properties)

#     props = _my_props_to_dict(
#         rp, input_labels, intensity_image=intensity_img, properties=properties, extra_properties=extra_properties
#     )

#     props["surface_area"] = surface_area_from_props(input_labels, props)
#     props_table = pd.DataFrame(props)
#     props_table.rename(columns={"area": "volume"}, inplace=True)
#     #  # ETC.  skeletonize via cellprofiler /Users/ahenrie/Projects/Imaging/CellProfiler/cellprofiler/modules/morphologicalskeleton.py
#     #         if x.volumetric:
#     #             y_data = skimage.morphology.skeletonize_3d(x_data)
#     # /Users/ahenrie/Projects/Imaging/CellProfiler/cellprofiler/modules/measureobjectskeleton.py

#     return props_table, rp


# creating a function to measure the surface area of each object. This function utilizes "marching_cubes" to generate a mesh (non-pixelated object)

### USED ###
def surface_area_from_props(labels, props, scale: Union[tuple,None]=None):
    # SurfaceArea
    surface_areas = np.zeros(len(props["label"]))
    # TODO: spacing = [1, 1, 1] # this is where we could deal with anisotropy in Z

    for index, lab in enumerate(props["label"]):
        # this seems less elegant than you might wish, given that regionprops returns a slice,
        # but we need to expand the slice out by one voxel in each direction, or surface area freaks out
        volume = labels[
            max(props["bbox-0"][index] - 1, 0) : min(props["bbox-3"][index] + 1, labels.shape[0]),
            max(props["bbox-1"][index] - 1, 0) : min(props["bbox-4"][index] + 1, labels.shape[1]),
            max(props["bbox-2"][index] - 1, 0) : min(props["bbox-5"][index] + 1, labels.shape[2]),
        ]
        volume = volume == lab
        if scale is None:
            scale=(1.0,) * labels.ndim
        verts, faces, _normals, _values = marching_cubes(
            volume,
            method="lewiner",
            spacing=scale,
            level=0,
        )
        surface_areas[index] = mesh_surface_area(verts, faces)

    return surface_areas


def _assert_uint16_labels(inp: np.ndarray) -> np.ndarray:
    """
    wrapper to enforce having the right labels
    """
    if inp.dtype == "bool" or inp.dtype == np.uint8:
        return label(inp > 0).astype(np.uint16)
    return inp

### USED ###
def get_region_morphology_3D(region_seg: np.ndarray, 
                              region_name: str,
                              intensity_img: np.ndarray, 
                              channel_names: List[str],
                              mask: np.ndarray, 
                              scale: Union[tuple, None]=None) -> Tuple[Any, Any]:
    """
    Parameters
    ------------
    region_seg:
        a 3D (ZYX) np.ndarray of segmented objects 
    region_name: str
        a name or nickname of the object being measured; this will be used for record keeping in the output table
    intensity_img:
        a 3D (ZYX) np.ndarray contain gray scale values from the "raw" image the segmentation is based on )single channel)
    mask:
        a 3D (ZYX) binary np.ndarray mask of the area to measure from
    scale: tuple, optional
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)


    Regionprops measurements:
    ------------------------
    ['label',
    'centroid',
    'bbox',
    'area',
    'equivalent_diameter',
    'extent',
    'feret_diameter_max',
    'euler_number',
    'convex_area',
    'solidity',
    'axis_major_length',
    'axis_minor_length',
    'max_intensity',
    'mean_intensity',
    'min_intensity']

    Additional measurements:
    -----------------------
    ['standard_deviation_intensity',
    'surface_area']


    Returns
    -------------
    pandas dataframe of containing regionprops measurements (columns) for each object in the segmentation image (rows) and the regionprops object

    """
    if len(channel_names) != intensity_img.shape[0]:
        ValueError("You have not provided a name for each channel in the intensity image. Make sure there is a channel name for each channel in the intensity image.")
    
    ###################################################
    ## MASK THE REGION OBJECTS THAT WILL BE MEASURED
    ###################################################
    # in case we sent a boolean mask (e.g. cyto, nucleus, cellmask)
    input_labels = _assert_uint16_labels(region_seg)

    input_labels = apply_mask(input_labels, mask)

    ##########################################
    ## CREATE LIST OF REGIONPROPS MEASUREMENTS
    ##########################################
    # start with LABEL
    properties = ["label"]
    # add position
    properties = properties + ["centroid", "bbox"]
    # add area
    properties = properties + ["area", "equivalent_diameter"] # "num_pixels", 
    # add shape measurements
    properties = properties + ["extent", "euler_number", "solidity", "axis_major_length"] # ,"feret_diameter_max", , "axis_minor_length"]
    # add intensity values (used for quality checks)
    properties = properties + ["min_intensity", "max_intensity", "mean_intensity"]

    #######################
    ## ADD EXTRA PROPERTIES
    #######################
    def standard_deviation_intensity(region, intensities):
        return np.std(intensities[region])

    extra_properties = [standard_deviation_intensity]

    ##################
    ## RUN REGIONPROPS
    ##################
    intensity_input = np.moveaxis(intensity_img, 0, -1)

    rp = regionprops(input_labels, 
                    intensity_image=intensity_input, 
                    extra_properties=extra_properties, 
                    spacing=scale)

    props = regionprops_table(label_image=input_labels, 
                              intensity_image=intensity_input, 
                              properties=properties, 
                              extra_properties=extra_properties,
                              spacing=scale)

    props_table = pd.DataFrame(props)
    props_table.insert(0, "object", region_name)
    props_table.rename(columns={"area": "volume"}, inplace=True)

    if scale is not None:
        round_scale = (round(scale[0], 4), round(scale[1], 4), round(scale[2], 4))
        props_table.insert(loc=2, column="scale", value=f"{round_scale}")
    else: 
        props_table.insert(loc=2, column="scale", value=f"{tuple(np.ones(region_seg.ndim))}") 

    rename_dict = {}
    for col in props_table.columns:
        for idx, name in enumerate(channel_names):
            if col.endswith(f"intensity-{idx}"):
                rename_dict[f"{col}"] = f"{col[:-1]}{name}_ch"

    props_table = props_table.rename(columns=rename_dict)

    ##################################################################
    ## RUN SURFACE AREA FUNCTION SEPARATELY AND APPEND THE PROPS_TABLE
    ##################################################################
    surface_area_tab = pd.DataFrame(surface_area_from_props(input_labels, props, scale))

    props_table.insert(12, "surface_area", surface_area_tab)
    props_table.insert(14, "SA_to_volume_ratio", props_table["surface_area"].div(props_table["volume"]))

    ################################################################
    ## ADD SKELETONIZATION OPTION FOR MEASURING LENGTH AND BRANCHING
    ################################################################
    #  # ETC.  skeletonize via cellprofiler /Users/ahenrie/Projects/Imaging/CellProfiler/cellprofiler/modules/morphologicalskeleton.py
    #         if x.volumetric:
    #             y_data = skimage.morphology.skeletonize_3d(x_data)
    # /Users/ahenrie/Projects/Imaging/CellProfiler/cellprofiler/modules/measureobjectskeleton.py

    return props_table

### USED ###
def get_contact_metrics_3D(a: np.ndarray,
                            a_name: str, 
                            b: np.ndarray, 
                            b_name:str, 
                            mask: np.ndarray, 
                            scale: Union[tuple, None]=None,
                            include_dist:bool=False, 
                            dist_centering_obj: Union[np.ndarray, None]=None,
                            dist_num_bins: Union[int, None]=None,
                            dist_zernike_degrees: Union[int, None]=None,
                            dist_center_on: Union[bool, None]=None,
                            dist_keep_center_as_bin: Union[bool, None]=None):
    """
    collect volumentric measurements of organelle `a` intersect organelle `b`

    Parameters
    ------------
    a: np.ndarray
        3D (ZYX) np.ndarray of one set of objects that will be assessed as part of a "contact"
    a_name: str
        the name or nickname of object a; this will be used for record keeping purposed in the output dataframe 
    b: np.ndarray
        3D (ZYX) np.ndarray of one set of objects that will be assessed as part of a "contact"
    b_name: str
        the name or nickname of object b; this will be used for record keeping purposed in the output dataframe 
    mask: np.ndarray
        3D (ZYX) binary mask of the area to measure contacts from
    include_dist:bool=False
        *optional*
        True = include the XY and Z distribution measurements of the contact sites within the masked region 
        (utilizing the functions get_XY_distribution() and get_Z_distribution() from Infer-subc)
        False = do not include distirbution measurements
    dist_centering_obj: Union[np.ndarray, None]=None
        ONLY NEEDED IF include_dist=True; if None, the center of the mask will be used
        3D (ZYX) np.ndarray containing the object to use for centering the XY distribution mask
    dist_num_bins: Union[int, None]=None
        ONLY NEEDED IF include_dist=True; if None, the default is 5
    dist_zernike_degrees: Unions[int, None]=None,
        ONLY NEEDED IF include_dist=True; if None, the zernike share measurements will not be included in the distribution
        the number of zernike degrees to include for the zernike shape descriptors
    dist_center_on: Union[bool, None]=None
        ONLY NEEDED IF include_dist=True; if None, the default is False
        True = distribute the bins from the center of the centering object
        False = distribute the bins from the edge of the centering object
    dist_keep_center_as_bin: Union[bool, None]=None
        ONLY NEEDED IF include_dist=True; if None, the default is True
        True = include the centering object area when creating the bins
        False = do not include the centering object area when creating the bins


    Regionprops measurements:
    ------------------------
    ['label',
    'centroid',
    'bbox',
    'area',
    'equivalent_diameter',
    'extent',
    'feret_diameter_max',
    'euler_number',
    'convex_area',
    'solidity',
    'axis_major_length',
    'axis_minor_length']

    Additional measurements:
    ----------------------
    ['surface_area']

    
    Returns
    -------------
    pandas dataframe of containing regionprops measurements (columns) for each overlap region between a and b (rows)
    
    """
    
    #########################
    ## CREATE OVERLAP REGIONS
    #########################
    a = _assert_uint16_labels(a)
    b = _assert_uint16_labels(b)

    a_int_b = np.logical_and(a > 0, b > 0)

    labels = label(apply_mask(a_int_b, mask)).astype("int")

    ##########################################
    ## CREATE LIST OF REGIONPROPS MEASUREMENTS
    ##########################################
    # start with LABEL
    properties = ["label"]

    # add position
    properties = properties + ["centroid", "bbox"]

    # add area
    properties = properties + ["area", "equivalent_diameter"] # "num_pixels", 

    # add shape measurements
    properties = properties + ["extent", "euler_number", "solidity", "axis_major_length", "slice"] # "feret_diameter_max", "axis_minor_length", 

    ##################
    ## RUN REGIONPROPS
    ##################
    props = regionprops_table(labels, intensity_image=None, properties=properties, extra_properties=None, spacing=scale)

    ##################################################################
    ## RUN SURFACE AREA FUNCTION SEPARATELY AND APPEND THE PROPS_TABLE
    ##################################################################
    surface_area_tab = pd.DataFrame(surface_area_from_props(labels, props, scale))

    ######################################################
    ## LIST WHICH ORGANELLES ARE INVOLVED IN THE CONTACTS
    ######################################################
    label_a = []
    index_ab = []
    label_b = []
    for index, lab in enumerate(props["label"]):
        # this seems less elegant than you might wish, given that regionprops returns a slice,
        # but we need to expand the slice out by one voxel in each direction, or surface area freaks out
        volume = labels[props["slice"][index]]
        la = a[props["slice"][index]]
        lb = b[props["slice"][index]]
        volume = volume == lab
        la = la[volume]
        lb = lb[volume]

        all_as = np.unique(la[la>0]).tolist()
        all_bs = np.unique(lb[lb>0]).tolist()
        if len(all_as) != 1:
            print(f"we have an error.  as-> {all_as}")
        if len(all_bs) != 1:
            print(f"we have an error.  bs-> {all_bs}")

        label_a.append(f"{all_as[0]}" )
        label_b.append(f"{all_bs[0]}" )
        index_ab.append(f"{all_as[0]}_{all_bs[0]}")


    ######################################################
    ## CREATE COMBINED DATAFRAME OF THE QUANTIFICATION
    ######################################################
    props_table = pd.DataFrame(props)
    props_table.drop(columns=['slice', 'label'], inplace=True)
    props_table.insert(0, 'label',value=index_ab)
    props_table.insert(0, "object", f"{a_name}X{b_name}")
    props_table.rename(columns={"area": "volume"}, inplace=True)

    props_table.insert(11, "surface_area", surface_area_tab)
    props_table.insert(13, "SA_to_volume_ratio", props_table["surface_area"].div(props_table["volume"]))

    if scale is not None:
        round_scale = (round(scale[0], 4), round(scale[1], 4), round(scale[2], 4))
        props_table.insert(loc=2, column="scale", value=f"{round_scale}")
    else: 
        props_table.insert(loc=2, column="scale", value=f"{tuple(np.ones(labels.ndim))}") 


    ######################################################
    ## optional: DISTRIBUTION OF CONTACTS MEASUREMENTS
    ######################################################
    if include_dist is True:
        XY_contact_dist, XY_bins, XY_wedges = get_XY_distribution(mask=mask, 
                                                                  obj=a_int_b,
                                                                  obj_name=f"{a_name}X{b_name}",
                                                                  centering_obj=dist_centering_obj,
                                                                  scale=scale,
                                                                  center_on=dist_center_on,
                                                                  keep_center_as_bin=dist_keep_center_as_bin,
                                                                  num_bins=dist_num_bins,
                                                                  zernike_degrees=dist_zernike_degrees)
        
        Z_contact_dist = get_Z_distribution(mask=mask,
                                            obj=a_int_b,
                                            obj_name=f"{a_name}X{b_name}",
                                            center_obj=dist_centering_obj,
                                            scale=scale)
        
        contact_dist_tab = pd.merge(XY_contact_dist, Z_contact_dist, on=["object", "scale"])

        return props_table, contact_dist_tab
    else:
        return props_table

# def get_aXb_stats_3D(a, b, mask, use_shell_a=False):
#     """
#     collect volumentric stats of `a` intersect `b`
#     """
#     properties = ["label"]  # our index to organelles
#     # add area
#     properties = properties + ["area", "equivalent_diameter"]
#     #  position:
#     properties = properties + ["centroid", "bbox"]  # ,  'weighted_centroid']
#     # etc
#     properties = properties + ["slice"]

#     # in case we sent a boolean mask (e.g. cyto, nucleus, cellmask)
#     a = _assert_uint16_labels(a)
#     b = _assert_uint16_labels(b)

#     if use_shell_a:
#         a_int_b = np.logical_and(np.logical_xor(a > 0, binary_erosion(a > 0)), b > 0)
#     else:
#         a_int_b = np.logical_and(a > 0, b > 0)

#     labels = label(apply_mask(a_int_b, mask)).astype("int")

#     props = regionprops_table(labels, intensity_image=None, properties=properties, extra_properties=None)

#     props["surface_area"] = surface_area_from_props(labels, props)
#     # there could be a bug here if there are spurious labels in the corners of the slices

#     label_a = []
#     index_ab = []
#     label_b = []
#     for index, lab in enumerate(props["label"]):
#         # this seems less elegant than you might wish, given that regionprops returns a slice,
#         # but we need to expand the slice out by one voxel in each direction, or surface area freaks out
#         volume = labels[props["slice"][index]]
#         la = a[props["slice"][index]]
#         lb = b[props["slice"][index]]
#         volume = volume == lab
#         la = la[volume]
#         lb = lb[volume]

#         all_as = np.unique(la[la>0]).tolist()
#         all_bs = np.unique(lb[lb>0]).tolist()
#         if len(all_as) != 1:
#             print(f"we have an error.  as-> {all_as}")
#         if len(all_bs) != 1:
#             print(f"we have an error.  bs-> {all_bs}")

#         label_a.append(all_as[0] )
#         label_b.append(all_bs[0] )
#         index_ab.append(f"{all_as[0]}_{all_bs[0]}")


#     props["label_a"] = label_a #[np.unique(a[s])[:1].tolist() for s in props["slice"]]
#     props["label_b"] = label_b #[np.unique(b[s])[:1].tolist() for s in props["slice"]]
#     props_table = pd.DataFrame(props)
#     props_table.rename(columns={"area": "volume"}, inplace=True)
#     props_table.drop(columns="slice", inplace=True)
#     props_table.insert(loc=0,column='label_',value=index_ab)
#     props_table.insert(loc=0,column='shell',value=use_shell_a)

#     return props_table

###################################
### DISTRIBUTIONAL STATS
###################################


# taken from cellprofiler_core.utilities.core.object
def size_similarly(labels, secondary):
    """Size the secondary matrix similarly to the labels matrix

    labels - labels matrix
    secondary - a secondary image or labels matrix which might be of
                different size.
    Return the resized secondary matrix and a mask indicating what portion
    of the secondary matrix is bogus (manufactured values).

    Either the mask is all ones or the result is a copy, so you can
    modify the output within the unmasked region w/o destroying the original.
    """
    if labels.shape[:2] == secondary.shape[:2]:
        return secondary, np.ones(secondary.shape, bool)
    if labels.shape[0] <= secondary.shape[0] and labels.shape[1] <= secondary.shape[1]:
        if secondary.ndim == 2:
            return (
                secondary[: labels.shape[0], : labels.shape[1]],
                np.ones(labels.shape, bool),
            )
        else:
            return (
                secondary[: labels.shape[0], : labels.shape[1], :],
                np.ones(labels.shape, bool),
            )

    # Some portion of the secondary matrix does not cover the labels
    result = np.zeros(
        list(labels.shape) + list(secondary.shape[2:]), secondary.dtype
    )
    i_max = min(secondary.shape[0], labels.shape[0])
    j_max = min(secondary.shape[1], labels.shape[1])
    if secondary.ndim == 2:
        result[:i_max, :j_max] = secondary[:i_max, :j_max]
    else:
        result[:i_max, :j_max, :] = secondary[:i_max, :j_max, :]
    mask = np.zeros(labels.shape, bool)
    mask[:i_max, :j_max] = 1
    return result, mask

# DEPRICATE
def get_simple_stats_3D(a, mask):
    """collect volumentric stats of `a`"""

    properties = ["label"]  # our index to organelles
    # add area
    properties = properties + ["area", "equivalent_diameter"]
    #  position:
    properties = properties + ["centroid", "bbox"]  # ,  'weighted_centroid']
    # etc
    properties = properties + ["slice"]

    # in case we sent a boolean mask (e.g. cyto, nucleus, cellmask)
    labels = _assert_uint16_labels(a)

    # props = regionprops_table(labels, intensity_image=None,
    #                             properties=properties, extra_properties=[])

    rp = regionprops(labels, intensity_image=None, extra_properties=[])
    props = _my_props_to_dict(rp, labels, intensity_image=None, properties=properties, extra_properties=None)

    stats_table = pd.DataFrame(props)
    stats_table.rename(columns={"area": "volume"}, inplace=True)

    return stats_table, rp


# untested 2D version
def get_summary_stats_2D(input_labels, intensity_img, mask):
    """collect volumentric stats"""

    # in case we sent a boolean mask (e.g. cyto, nucleus, cellmask)
    input_labels = _assert_uint16_labels(input_labels)

    # mask
    input_labels = apply_mask(input_labels, mask)

    properties = ["label"]
    # add intensity:
    properties = properties + ["max_intensity", "mean_intensity", "min_intensity"]

    # arguments must be in the specified order, matching regionprops
    def standard_deviation_intensity(region, intensities):
        return np.std(intensities[region])

    extra_properties = [standard_deviation_intensity]

    # add area
    properties = properties + ["area", "equivalent_diameter"]
    #  position:
    properties = properties + ["centroid", "bbox"]

    #  perimeter:
    properties = properties + ["perimeter", "perimeter_crofton"]

    rp = regionprops(input_labels, intensity_image=intensity_img, extra_properties=extra_properties)
    props = _my_props_to_dict(
        rp, input_labels, intensity_image=intensity_img, properties=properties, extra_properties=extra_properties
    )
    props_table = pd.DataFrame(props)
    #  # ETC.  skeletonize via cellprofiler /Users/ahenrie/Projects/Imaging/CellProfiler/cellprofiler/modules/morphologicalskeleton.py
    #             y_data = skimage.morphology.skeletonize(x_data)
    # /Users/ahenrie/Projects/Imaging/CellProfiler/cellprofiler/modules/measureobjectskeleton.py
    return props_table, rp
