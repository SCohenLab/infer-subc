from typing import List, Union
from pathlib import Path
import warnings
import time

import numpy as np
import pandas as pd
from skimage.measure import regionprops_table, marching_cubes, mesh_surface_area

from infer_subc.core.img import *
from infer_subc.organelles import * 
from infer_subc.utils.batch import list_image_files, find_segmentation_tiff_files
from infer_subc.core.file_io import read_czi_image, read_tiff_image, export_inferred_organelle
from infer_subc.quantification.batch import append_atomic_csv, load_existing_keys_csv
from infer_subc.quantification.skeletonization import create_skel, get_skeleton_metrics


def surface_area_from_props(labels: np.ndarray,
                             props: dict,
                             scale: Union[tuple, None]=None):
    """ 
    a function for getting surface area of volumetric objects

    Parameters:
    ----------
    lables:
        the segmentation np.ndarray with each object labeled a different number
    props:
        region props dictionary resulting from the _my_props_to_dict() function
    spacing:
        tuple of the dimension lengths in the same order as the dimension of your np.ndarray labels input
    """
    surface_areas = np.zeros(len(props["label"]))

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
        verts, faces, _normals, _values = marching_cubes(volume,
                                                         method="lewiner",
                                                         spacing=scale,
                                                         level=0)
        
        surface_areas[index] = mesh_surface_area(verts, faces)

    return surface_areas


# Apply skimage regionprops function to quantify the morphology of a single object type from a single cell
def get_morphology_metrics(segmentation_img: np.ndarray, 
                            seg_name: str, 
                            intensity_img: Union[np.ndarray, None],
                            intensity_ch_names: List[str],
                            channel_axis: int,
                            mask: Union[np.ndarray, None]=None, 
                            mask_name: Union[str, None]=None,
                            scale: Union[tuple, None]=None):
    """
    Parameters
    ------------
    segmentation_img:
        an np.ndarray of segmented objects 
    seg_name: str
        a name or nickname (usually the segmentation file suffix) of the object being measured; this will be used for record keeping in the output table
    intensity_img:
        a single-channel np.ndarray contain gray scale values from the "raw" image the segmentation is based on; this image should be the same shape as the segmentation file
    intensity_ch_names: List[str]
        a list of names for each channel in the intensity image; used to rename intensity measurement columns
    channel_axis:
        the index of the channel dimension axis in the intensity image
    mask: Union[np.ndarray, None]
        a binary np.ndarray mask of the area to measure from; this image should be the same shape as the segmentation file
    mask_name: Union[str, None]
        the name of the mask region being analyzed
    scale: tuple, optional
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)


    Regionprops measurements:
    ------------------------
    'label',
    'centroid',
    'bbox',
    'area',
    'equivalent_diameter',
    'extent',
    'euler_number',
    'solidity',
    'axis_major_length',
    'min_intensity',
    'max_intensity',
    'mean_intensity'

    Additional measurements:
    -----------------------
    'standard_deviation_intensity',
    'surface_area',
    'SA_to_volume_ratio`


    Returns
    -------------
    pandas dataframe of containing regionprops measurements (columns) for each object in the segmentation image (rows) and the regionprops object
    
    """
    # dealing with numerous solidity warning from regionprops
    warnings.simplefilter("ignore")

    ###################################################
    ## MASK THE ORGANELLE OBJECTS THAT WILL BE MEASURED
    ###################################################
    # apply mask to overlap image if provided
    if mask_name is None and mask is not None:
        raise ValueError("The mask_name parameter must be provided if mask is not None")
    elif mask is None and mask_name is None:
        input_labels = segmentation_img
        mask_name = "whole_image"
    else:
        input_labels = apply_mask(segmentation_img, mask)

    ##########################################
    ## CREATE LIST OF REGIONPROPS MEASUREMENTS
    ##########################################
    properties = ["label", "centroid", "bbox", "area", 
                  "equivalent_diameter", "extent", "euler_number", "solidity", "axis_major_length",
                  "min_intensity", "max_intensity", "mean_intensity"]

    #######################
    ## ADD EXTRA PROPERTIES
    #######################
    def standard_deviation_intensity(region, intensities):
        return np.std(intensities[region])

    extra_properties = [standard_deviation_intensity]

    #########################
    ## ADJUST INTENSITY IMAGE
    #########################
    if scale is None:
        scale = (1.0,) * segmentation_img.ndim

    if intensity_img is not None:
        if channel_axis == len(scale):
            pass
        else:
            intensity_img = np.moveaxis(intensity_img, channel_axis, -1)

    ##################
    ## RUN REGIONPROPS
    ##################
    props = regionprops_table(input_labels, 
                           intensity_image=intensity_img, 
                           properties=properties,
                           extra_properties=extra_properties,
                           spacing=scale)
    
    # measure the mask volume as well for easier normalization in downstream functions
    mask_vol = regionprops_table(mask,properties=["area"], spacing=scale)['area'][0]

    props_table = pd.DataFrame(props)

    ##################################################################
    ## RUN SURFACE AREA FUNCTION SEPARATELY AND APPEND THE PROPS_TABLE
    ##################################################################
    surface_area_tab = pd.DataFrame(surface_area_from_props(input_labels, props, scale))

    #############################################
    ## RENAME AND ADD ADDITIONAL METADATA COLUMNS
    #############################################
    props_table.insert(0, "object", seg_name)
    props_table.rename(columns={"area": "volume"}, inplace=True)

    if scale is not None:
        round_scale = (round(scale[0], 4), round(scale[1], 4), round(scale[2], 4))
        props_table.insert(0, column="scale", value=f"{round_scale}")
    else: 
        props_table.insert(0, column="scale", value=f"{tuple(np.ones(segmentation_img.ndim))}") 

    props_table.insert(props_table.columns.get_loc('volume') + 1, "surface_area", surface_area_tab)
    props_table.insert(props_table.columns.get_loc('surface_area') + 1, "SA_to_volume_ratio", props_table["surface_area"].div(props_table["volume"]))
    props_table.insert(0, column="mask_name", value=mask_name)
    props_table[f"{mask_name}_volume"] = mask_vol

    for col in [c for c in props_table.columns if "intensity" in c]:
        props_table.rename(columns={col:col[:-1]+intensity_ch_names[int(col[-1])]+"-ch"}, inplace=True)

    # print this statement to let user known of suppressed warnings
    if Warning: print(f"Warning(s) suppressed while quantifying {seg_name}. See 'method_morphology.ipynb' notebook for more details.")

    return props_table


# quantify the morphology of one or more organelle from one cell
def get_org_morphology(source_file_path: str,
                        list_obj_names: List[str],
                        list_obj_segs: List[np.ndarray],
                        list_intensity_img: Union[List[np.ndarray], None]=None,
                        list_region_names: Union[List[str], None]=None,
                        list_region_segs: Union[List[np.ndarray], None]=None,
                        mask_name: Union[str, None]=None,
                        scale: Union[tuple,None] = None,
                        include_skel: Union[List[str], None] = [],
                        all_skel_tab: bool = False):
    """
    Measure the amount, size, and shape of multiple organelles from a single cell

    Parameters:
    ----------
    source_file: str
        Path to the source image file. This will be used as part of the metadata information in the output table. 
        The input images are not derived from this path, but rather are provided directly as arrays in the list_obj_segs and 
        list_intensity_img variables below.
    list_obj_names: List[str]
        List of organelle names. These names should match the suffix on the segmentation image files.
    list_obj_segs: List[np.ndarray]
        List of 3D organelle segmentation arrays matching the order included in list_obj_names.
    list_intensity_img: Union[List[np.ndarray], None]
        List of 3D intensity channels from the raw image used to produce the segmentations in list_obj_segs.
        The order here should match the list_obj_segs and list_obj_names variables.
        Additional intensity channels not matching one of the segmented organelles/included in list_obj_names should not be included.
        If no intensity analysis is to be included, specify None here.
    list_region_names: Union[List[str], None]
        List of segmented region/mask names. These names should match the suffix on the segmentation image files.
        This should include:
            - a mask segmentation, such as the cell mask, for masking during all interactions analysis; else, the entire image will be 
            quantified. Only one objects per mask image will be analyzed. If there are more than one included, they will be combined 
            prior to analysis and the entire region will be quantified. If no mask is provided, the entire image will be quantified.
            - a centering object, such as the nucleus, for distribution analysis; else the center of the mask region will be used as 
            the XY distribution centering point if distribution analysis is included.
    list_region_segs: Union[List[np.ndarray], None]
        List of 3D region segmentation arrays matching the order specified in list_region_names. Specify None if no regions are provided.
    mask_name: Union[str, None]
        Name of the region to use as the mask for analysis; if not specified, the entire image will be quantified.
    scale: Union[tuple,None] = None
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)
    include_skel: Union[List[str], None]
             List of 3D organelle segmentation arrays to be skeletonized in addition to morphology metrics.
             The names should match those included in list_obj_names. If no skeletonization is to be included, specify [] (empty list).
    all_skel_tab: bool
        Whether to output all skeleton tables (True) or just the skeleton metrics table (False).

    Returns:
    ----------
    Dataframe of measurements of organelle morphology

    """
    # Validate inputs
    if not list_obj_names:
        raise ValueError("list_obj_names cannot be empty")
    if len(list_obj_names) != len(list_obj_segs):
        raise ValueError(f"Mismatch: {len(list_obj_names)} items in list_obj_names but {len(list_obj_segs)} items in list_obj_segs")
    if list_intensity_img and len(list_obj_names) != len(list_intensity_img):
        raise ValueError(f"Mismatch: {len(list_obj_names)} items in list_obj_names but {len(list_intensity_img)} items in list_intensity_img")
    if list_region_names and list_region_segs and len(list_region_names) != len(list_region_segs):
        raise ValueError(f"Mismatch: {len(list_region_names)} items in list_region_names but {len(list_region_segs)} items in list_region_segs")
    
    if isinstance(source_file_path, str): source_file_path = Path(source_file_path)
    print(f"Quantifying organelle morphology from {source_file_path.name}")

    # specify the mask image to use during quantification
    if list_region_names is None or list_region_segs is None:
        print("No regions provided. No mask will be applied before analysis.")
        mask = None
    elif mask_name is None or mask_name not in list_region_names:
        if mask_name is not None:
            raise ValueError(f"Mask '{mask_name}' not found. No mask will be applied before analysis.")
        mask = None
        mask_name = None
    else:
        mask = list_region_segs[list_region_names.index(mask_name)]
    
    # merge intensity images to create a single np.ndarray
    if list_intensity_img is None:
        intensity_img = None
        print("No intensity images provided. Morphology metrics that require intensity images will not be calculated.")
    else:
        intensity_img = np.stack(list_intensity_img, axis=0)

    # empty list to collect a morphology data for each organelle
    org_tabs = []
    skel_arr_imgs = {}
    # additional skeleton tables if include_skel and all_skel_tab are both True
    branch_tabs = []
    node_tabs = []

    # loop through the list of organelles and run the get_morphology_metrics function
    for j, target in enumerate(list_obj_names):        
        # select segmentation and if ER, ensure it is only one object
        if target == 'ER':
            org_obj = (list_obj_segs[j] > 0).astype(np.uint16)  
        else:
            org_obj = list_obj_segs[j]
        
        # run get_morphology_metrics function to output a table of measurements
        org_metrics = get_morphology_metrics(segmentation_img=org_obj, 
                                            seg_name=target,
                                            intensity_img=intensity_img, 
                                            intensity_ch_names=list_obj_names,
                                            channel_axis=0, # default to 0 because intensities are merged from list on axis 0
                                            mask=mask,
                                            mask_name=mask_name,
                                            scale=scale)
        
        # run get_skeleton_metrics function to add skeleton quantification if organelle is in include_skel list
        if target in include_skel:
            skel_arr = create_skel(org_obj)
            if np.sum(skel_arr.astype(bool))>1:
                if all_skel_tab:
                    skel_branch_table, skel_node_table, skel_metrics = get_skeleton_metrics(org_skel_arr=skel_arr,
                                                        seg_name=target, 
                                                        segmentation=org_obj,
                                                        mask=mask,
                                                        mask_name=mask_name,
                                                        scale=scale,
                                                        output_all_tables = True)
                    
                    skel_branch_table.insert(0, "branch_id", skel_branch_table.index)
                    skel_branch_table = skel_branch_table.reset_index(drop=True)
                    skel_branch_table.insert(0, "object", target)
                    skel_branch_table.insert(0, "scale", str(scale))
                    skel_branch_table.insert(0, column="mask_name", value=mask_name)
                    skel_branch_table = skel_branch_table.rename(columns={"skel_obj_id": "label"})
                    branch_tabs.append(skel_branch_table)

                    skel_node_table.insert(0, "object", target)
                    skel_node_table.insert(0, "scale", str(scale))
                    skel_node_table.insert(0, column="mask_name", value=mask_name)
                    skel_node_table = skel_node_table.rename(columns={"obj_id": "label"})
                    node_tabs.append(skel_node_table)
                else:
                    skel_metrics = get_skeleton_metrics(org_skel_arr=skel_arr,
                                                        seg_name=target, 
                                                        segmentation=org_obj,
                                                        mask=mask,
                                                        mask_name=mask_name,
                                                        scale=scale,
                                                        output_all_tables = False)
                
            else:
                print(f"Skeletonization will not be carried out for {target} because less than two voxels are present in the skeleton array")
                skel_arr = None
            skel_arr_imgs[target] = skel_arr
            # Rename skeleton metric columns to distingush from morphology metrics
            skel_cols = [col for col in skel_metrics.columns[:6]]
            for i in skel_metrics.columns[6:]:
                skel_cols.append("skel_" + i)
            skel_metrics.columns = skel_cols

            # dropping list related measurements as they do not add interpretability and cause tables to be very tall and difficult to read
            skel_drop = ['skel_branch_ids',
                'skel_brh_type_0_id',
                'skel_brh_type_1_ids',
                'skel_brh_type_2_ids',
                'skel_brh_type_3_ids',
                'skel_point_ids']

            org_metrics = pd.merge(org_metrics, skel_metrics.drop(columns=skel_drop), on=['mask_name','scale','object', 'label'], validate='one_to_one')

        # add table to list above
        org_tabs.append(org_metrics)

    # combine the lists for each organelle into one table
    final_org_tab = pd.concat(org_tabs, ignore_index=True)

    # add a new column to list the name of the image these data are derived from 
    final_org_tab.insert(loc=0,column='image_name',value=source_file_path.stem)

    if include_skel and all_skel_tab:
        final_branch_tab = pd.concat(branch_tabs)
        final_branch_tab.insert(0, "image_name", source_file_path.stem)
        final_node_tab = pd.concat(node_tabs)
        final_node_tab.insert(0, "image_name", source_file_path.stem)
        return final_org_tab, final_branch_tab, final_node_tab, skel_arr_imgs
    else:
        return final_org_tab, skel_arr_imgs


# batch process organelle morphology quantification for multiple cells from a single experiment
def batch_process_org_morph(dataset_name: str,
                             raw_path: Union[Path,str], 
                             seg_path: Union[Path,str],
                             quant_path: Union[Path, str], 
                             raw_file_type: str,
                             channel_axis: int,
                             organelle_names: List[str],
                             organelle_channels: Union[List[int], None]=None,
                             region_names: Union[List[str], None]=None,
                             mask_name: Union[str, None]=None,
                             use_scale:bool=True,
                             include_skel: Union[List[str], None] = [],
                             all_skel_tab: bool = False,
                             seg_suffix:Union[str, None]=None):
    """  
    batch process segmentation quantification (morphology, distribution, contacts); this function is currently optimized to process images from one file folder per image type (e.g., raw, segmentation)
    the output csv files are saved to the indicated quant_path folder

    Parameters:
    ----------
    dataset_name : str
        A unique string identifier for the dataset being processed. It will be included as metadata in output tables and as 
        part of the output files names. It will be used to identify if any data has already been collected for this dataset.
    raw_path: Union[Path,str]
        Path or str to the folder that contains the raw image files
    seg_path: Union[Path,str]
        Path or str to the folder that contains the segmentation tiff files
    quant_path: Union[Path, str]
        Path or str to the folder that the output datatables will be saved to
    raw_file_type: str
        File type of the raw images (e.g., "czi", "tiff")
    channel_axis : int
        Axis corresponding to the channels in the image data
    organelle_names: List[str]
        List of organelle names to analyze. These names should match the suffix on the organelle segmentation files
    organelle_channels: Union[List[int], None]=None
        List of intensity channel indices in the raw files corresponding to each organelle included in organelle_names.
        The order should match organelle_names. 
        If no intensity analysis is to be included, specify None here.
    region_names: Union[List[str], None]=None
        List of region names to analyze. Usually ['cell', 'nuc'] for cell mask and nucleus.
        If no regions are to be included, specify None here.
    mask_name: Union[str, None]=None
        Name of the region to use for segmentation (if any). This name should be included in the regions_name variable.
        If None, the entire image will be quantified.
    use_scale: bool=True
        Whether to apply scaling to the quantitative data; scaled data will be in real world units (e.g., microns) rather than pixels/voxels
    include_skel: Union[List[str], None]
        List of 3D organelle segmentation arrays to be skeletonized in addition to morphology metrics.
        The names should match those included in list_obj_names. If no skeletonization is to be included, specify [] (empty list).
    all_skel_tab: bool
        Whether to output all skeleton tables (True) or just the skeleton metrics table (False).
    seg_suffix:Union[str, None]=None
        Any additional text that is included in the segmentation tiff files between the file stem and the segmentation suffix, not including the initial "-"

    Returns:
    ----------
    None: files are saved to quant_path directly
    """
    
    start = time.time()
    count = 0

    # create path objects if inputs are strings
    if isinstance(raw_path, str): raw_path = Path(raw_path)
    if isinstance(seg_path, str): seg_path = Path(seg_path)
    if isinstance(quant_path, str): quant_path = Path(quant_path)
    
    # create directory is it doesn't exist
    if not Path.exists(quant_path):
        Path.mkdir(quant_path)
        print(f"Output file path not found. Making {quant_path}.")

    # check if any existing data is present in outfiles to skip already processed images
    unique_keys = ['dataset', 'image_name']

    morpho_path = quant_path / f"{dataset_name}_org_morphology_metrics.csv"
    existing_morpho_keys = load_existing_keys_csv(morpho_path, unique_keys)

    # reading list of files from the raw path
    img_file_list = list_image_files(raw_path, raw_file_type)
    len_file_list = len(img_file_list)

    # list of organelle segmentation and masks files to collect from each image
    segs_to_collect = organelle_names + region_names if region_names is not None else organelle_names

    # loop through list of cell analyzing each and appending the data to the empty list
    for img_f in img_file_list:
        img_start = time.time()
        count = count + 1
        # skip files that have already been processed
        if (dataset_name, img_f.stem) in existing_morpho_keys:
            print(f"Skipping {img_f.name} as it is already listed in the output file(s).")
            continue
        # process analysis for this cells
        else:
            filez = find_segmentation_tiff_files(img_f, segs_to_collect, seg_path, seg_suffix)

            # read in raw file and metadata
            img_data, meta_dict = read_czi_image(filez["raw"])

            # create intensities from raw file as list based on the channel order provided
            if organelle_channels is None:
                intensities = None
                print("No intensity channel information provided.")
            else:
                if channel_axis != 0:
                    img_data = np.moveaxis(img_data, channel_axis, 0)
                intensities = [img_data[ch] for ch in organelle_channels]

            # store organelle images as list
            organelles = [read_tiff_image(filez[org]) for org in organelle_names]

            # load regions as a list based on order in list (should match order in "masks" file)
            if region_names is None:
                regions = None
            else:
                regions = [read_tiff_image(filez[r]) for r in region_names] 

            # define the scale
            if use_scale is True:
                scale = meta_dict['scale']
            else:
                scale = None

            if include_skel and all_skel_tab:
                org_metrics, branch_table, node_table, skel_dict_arr = get_org_morphology(source_file_path=img_f,
                                                    list_obj_names=organelle_names,
                                                    list_obj_segs=organelles,
                                                    list_intensity_img=intensities,
                                                    list_region_names=region_names,
                                                    list_region_segs=regions,
                                                    mask_name=mask_name,
                                                    scale=scale,
                                                    include_skel=include_skel,
                                                    all_skel_tab=all_skel_tab)
                # save the morphology (or labels only) table data per image directly to csv
                org_metrics.insert(loc=0,column='dataset',value=dataset_name)
                append_atomic_csv(morpho_path, org_metrics)
                del org_metrics  # free up memory

                # save branch and node tables if skeletonization is included
                branch_path = quant_path / f"{dataset_name}_skeleton_branch_data.csv"
                node_path = quant_path / f"{dataset_name}_skeleton_node_data.csv"
                append_atomic_csv(branch_path, branch_table)
                append_atomic_csv(node_path, node_table)
                del branch_table, node_table  # free up memory

                org_metrics, skel_dict_arr = get_org_morphology(source_file_path=img_f,
                                                    list_obj_names=organelle_names,
                                                    list_obj_segs=organelles,
                                                    list_intensity_img=intensities, 
                                                    list_region_names=region_names,
                                                    list_region_segs=regions, 
                                                    mask_name=mask_name,
                                                    scale=scale,
                                                    include_skel=include_skel)

                # save the morphology table data per image directly to csv
                org_metrics.insert(loc=0,column='dataset',value=dataset_name)
                append_atomic_csv(morpho_path, org_metrics)
                del org_metrics  # free up memory
            
            # save the morphology table data per image directly to csv
            org_metrics.insert(loc=0,column='dataset',value=dataset_name)
            append_atomic_csv(morpho_path, org_metrics)
            del org_metrics  # free up memory

            if include_skel:
             # save the skeleton images
                for skel_name, skel_img in skel_arr_dict.items():
                    if skel_img is not None:
                        if not (Path(img_f)/f"{meta_dict['file_name'].stem}-{skel_name}-skeleton.tiff").exists():
                            export_inferred_organelle(skel_img, f"{skel_name}-skeleton", meta_dict, img_f)
                        else:
                            warnings.warn(f"Some of the skeleton images already exist for {meta_dict['file_name'].stem} in {img_f}. They will not be overwritten.", UserWarning)
                del skel_arr_dict  # free up memory

            end2 = time.time()
            print(f"Completed quantification of {meta_dict['file_name']} in {(end2-img_start)/60} mins.")
            print(f"{count}/{len_file_list} images have been processed.")
            print(f"Time elapsed: {(end2-img_start)/60} mins")

    end = time.time()
    print(f"Quantification for {count} files is COMPLETE! Files saved to '{quant_path}'.")
    print(f"It took {(end - start)/60} minutes to quantify these files.")


# summarize morphology values per organelle per cell across one or more experiments
def batch_org_morph_summary_stats(csv_path_list: List[str],
                                   out_path: str,
                                   out_prefix: str,
                                   organelle_names: List[str]):
    """" 
    csv_path_list: List[str],
        A list of path strings where .csv files to analyze are located.
    out_path: str,
        A path string where the summary data file will be output to
    out_prefix: str
        The prefix used to name the output file. An "_" will be included between this prefix and the file suffix.
    organelle_names: List[str],
        A list of organelle names used in the organelle morphology quantification (batch_process_org_morph function)
    """
    # for keeping track of dataset and file numbers
    ds_count = 0
    fl_count = 0

    ###################
    # Read in the csv files and combine them into one of each type
    ###################
    # create empty list to hold the morphology tables from different experiments
    org_tabs = []

    # loop through all of the locations listed above and find the _org_morph files; append them to the list above
    for loc in csv_path_list:
        # list all csv files in the location
        files_store = sorted(loc.glob("*.csv"))

        # find the unique datasets in this location based on the prefixes before "_org_morphology_metrics"
        prefixes = set(f.name.split("_org_morphology_metrics")[0] for f in files_store if "_org_morphology_metrics" in f.name)
        for prefix in prefixes:
            ds_count += 1
            # select only the files from this dataset
            files_subset = [f for f in files_store if f.name.startswith(prefix +"_org_morphology_metrics")]
            for file in files_subset:
                fl_count += 1
                stem = file.stem
                if "_org_morph" in stem:
                    test_orgs = pd.read_csv(file, index_col=0)
                    org_tabs.append(test_orgs)

    # combine the org_morph lists found above into one table
    org_df = pd.concat(org_tabs,axis=0, join='outer').reset_index()

    print(f"Found {fl_count} files from {ds_count} dataset(s) across {len(csv_path_list)} location(s).")


    ###################
    # summary stat group
    ###################
    group_by = ['dataset', 'image_name', 'mask_name', 'scale', 'object']
    sharedcolumns = ["SA_to_volume_ratio", "equivalent_diameter", "extent", "euler_number", "solidity", "axis_major_length"]
    
    # check if skeletonization metrics are included in the data by looking for any column names that start with "skel_"
    include_skel = any(col.startswith('skel_') for col in org_df.columns)

    # add skeleton metrics to the shared columns if they are included in the data
    if include_skel:
        skel_cols = ["skel_total_length",
                    "skel_abs_punc_count",
                    "skel_ep_count",
                    "skel_jn_count",
                    "skel_node_count",
                    "skel_brh_count"]
        skel_cols_2 = ["skel_brh_type_0_tot",
                            "skel_brh_type_1_tot",
                            "skel_brh_type_2_tot",
                                "skel_brh_type_3_tot",
                            "skel_comp_count",
                            "skel_ave_jn_deg",
                            "skel_max_deg",
                            "skel_mean_brh_str",
                            "skel_width"]
    ag_func_standard = ['mean', 'median', 'std']


    ###################
    # summarize shared measurements between org_df and contacts_df
    ###################
    tab1 = org_df[group_by + ['label']].groupby(group_by).agg(['count'])
    tab1.rename(columns={'label': 'org'}, inplace=True)
    tab2 = org_df[group_by + ['volume', 'surface_area']].groupby(group_by).agg(['sum'] + ag_func_standard)
    tab3 = org_df[group_by + sharedcolumns].groupby(group_by).agg(ag_func_standard)
    org_summary = pd.merge(tab1, tab2, 'outer', on=group_by)
    org_summary = pd.merge(org_summary, tab3, 'outer', on=group_by)
    if include_skel:
        tab4 = org_df[group_by + skel_cols].groupby(group_by).agg(['sum'] + ag_func_standard)
        tab5 = org_df[group_by + skel_cols_2].groupby(group_by).agg(ag_func_standard)
        org_summary = pd.merge(org_summary, tab4, 'outer', on=group_by)
        org_summary = pd.merge(org_summary, tab5, 'outer', on=group_by)

    # Get mask_name and corresponding volume column per group & calculate volume fraction
    mask_names = org_df.groupby(group_by)['mask_name'].first()
    mask_volume_data = org_df.groupby(group_by).first().apply(lambda row: row[f"{mask_names.loc[row.name]}_volume"], axis=1)
    org_summary.insert(org_summary.columns.get_loc(('volume', 'sum')) + 1, ('volume', 'fraction'), org_summary[('volume', 'sum')]/mask_volume_data)

    ###################
    # additional skeleton summarization inspired by mitograph measurments (if skeletonization is included)
    ###################

    if include_skel:
        # a temporary table used to calculate the average node degree, fusion score, fission score, connectivity, and heterogeneity for each organelle per cell
        # these measurements are inspired by the mitograph quantification metrics, but adapted to be more generalizable to different organelles
        skel_ff = org_df.assign(sum_jnxdeg = lambda x: x['skel_ave_jn_deg'] * x['skel_jn_count'],
                                is_punc = lambda df: (df['skel_type'] == "Punctate").astype(int),
                                is_rod = lambda df: (df['skel_type'] == "Rod").astype(int),
                                is_iso_cycle = lambda df: (df['skel_type'] == "Isolated Cycle").astype(int),
                                is_network = lambda df: (df['skel_type'] == "Network").astype(int)).groupby(group_by).agg({
            'label': 'count',
            'volume': ['mean','std'],
            'skel_total_length': ['max','sum','mean','std'],
            'skel_brh_count': ['sum','mean','std'],
            'skel_abs_punc_count': ['sum'],
            'skel_ep_count': 'sum',
            'skel_jn_count': ['sum'],
            'skel_node_count': ['sum','mean','std'],
            'sum_jnxdeg': 'sum',
            'skel_width': ['mean','std'],
            'is_punc': 'sum',
            'is_rod': 'sum',
            'is_iso_cycle': 'sum',
            'is_network': 'sum'
        }).assign(
            punctate_count = lambda df: df['is_punc'],
            rod_count = lambda df: df['is_rod'],
            iso_cycle_count = lambda df: df['is_iso_cycle'],
            network_count = lambda df: df['is_network'],
            skel_avg_node_deg = lambda df: (df['sum_jnxdeg','sum'] + df['skel_ep_count','sum']) / df['skel_node_count','sum'],
            skel_fusion_score = lambda df: (df['skel_total_length','max']/df['skel_total_length','sum']) + 
                (df['skel_total_length','sum']/df['skel_brh_count','sum']) + df['skel_avg_node_deg'],
            skel_fission_score = lambda df: (df['label','count']/df['skel_total_length','sum']) + 
                (df['skel_node_count','sum']/df['skel_total_length','sum']) + (df['skel_brh_count','sum']/df['skel_total_length','sum']),
            skel_connectivity = lambda df: df['skel_fusion_score']/df['skel_fission_score'],
            skel_heterogeneity = lambda df: (df['skel_node_count','std']/df['skel_node_count','mean']) + 
                (df['skel_brh_count','std']/df['skel_brh_count','mean']) +
                (df['skel_total_length','std']/df['skel_total_length','mean']) +
                (df['volume','std']/df['volume','mean']) +
                (df['skel_width','std']/df['skel_width','mean']) +
                df['skel_avg_node_deg']
        )

        # columns from the fusion/fission table to add to the morphology summary table
        skel_sum_metrics = ["punctate_count", "rod_count", "iso_cycle_count", "network_count",
            "skel_avg_node_deg", "skel_fusion_score", "skel_fission_score", "skel_connectivity", "skel_heterogeneity"]

        # merge these additional metrics back into the org_summary table
        org_summary = pd.merge(org_summary, skel_ff[skel_sum_metrics], on=group_by, how='outer')

    ###################
    # fill gaps & NA values
    ###################
    # Ensure all possible interactions are represented (if missing fill with NaN)
    for ind in org_summary.index.droplevel(4).unique().to_list():
        for row in organelle_names:
            if ind+(row,) not in org_summary.index:
                org_summary.loc[ind+(row,)] = np.nan

    # fill NA with 0 for specific columns
    fill_dict = {('org', 'count'): 0, 
                ('volume', 'sum'): 0,
                ('surface_area', 'sum'): 0,
                ('volume', 'fraction'): 0}
    org_summary = org_summary.fillna(value=fill_dict)

    # if (org, count) is 1, set mean, median, and std to NaN
    single_site_mask = org_summary[('org', 'count')] == 1
    for col in sharedcolumns+['volume', 'surface_area']:
        org_summary.loc[single_site_mask, (col, 'std')] = np.nan

    org_summary.sort_index(inplace=True)

    ###################
    # flatten datasheet and export
    ###################
    # export before unstacking
    if (Path(out_path) / f"{out_prefix}_per_org_morphology_summarystats.csv").exists():
        raise FileExistsError(f"CAUTION: {out_prefix}_per_org_morphology_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `quant_data_path` to continue without error.")
    else:
        org_summary.to_csv(str(out_path) + f"/{out_prefix}_per_org_morphology_summarystats.csv", mode='x')
        print(f"Exported per-organelle morphology summary statistics (before unstacking) to {out_path}/{out_prefix}_per_org_morphology_summarystats.csv")

    org_morph_final = org_summary.unstack(-1)
    org_morph_final.columns = ["_".join((col_name[1], col_name[-1], col_name[0])) for col_name in org_morph_final.columns.to_flat_index()]
    org_morph_final.columns = [col.replace('sum', 'total') for col in org_morph_final.columns]
    org_morph_final.columns = [col.replace(col, 'mask_volume') if 'mask' in col else col for col in org_morph_final.columns]
    org_morph_final = org_morph_final.loc[:, ~org_morph_final.columns.duplicated()]
    org_morph_final.reset_index(inplace=True)

    ###################
    # export summary sheets
    ###################
    if (Path(out_path) / f"{out_prefix}_organelle_morphology_summarystats.csv").exists():
        raise FileExistsError(f"CAUTION: {out_prefix}_organelle_morphology_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `quant_data_path` to continue without error.")
    else:
        org_morph_final.to_csv(str(out_path) + f"/{out_prefix}_organelle_morphology_summarystats.csv", mode='x')
        print(f"Exported organelle morphology summary statistics (after unstacking) to {out_path}/{out_prefix}_organelle_morphology_summarystats.csv")
    print(f"Organelle morphology summary is complete.")
    return org_summary