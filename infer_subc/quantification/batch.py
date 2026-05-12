import os
import tempfile
from pandas import read_csv, concat
from pathlib import Path
from typing import List, Union
import time
import warnings

import numpy as np
import pandas as pd

from infer_subc.core.file_io import export_inferred_organelle
from infer_subc.utils.batch import list_image_files, find_segmentation_tiff_files
from infer_subc.core.file_io import read_czi_image, read_tiff_image
from infer_subc.quantification.morphology import get_org_morphology
from infer_subc.quantification.interactions import get_interaction_metrics, all_combos, perorg_interactions_cnt
from infer_subc.quantification.regions import get_regions_morphology
from infer_subc.quantification.distribution import get_distribution_metrics
from infer_subc.quantification.csv_io import append_atomic_csv, load_existing_keys_csv




# batch processing function for combined organelle signature analysis
def batch_process_quantification(dataset_name: str,
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
                                seg_suffix:Union[str, None]=None,

                                # morphology settings
                                include_org_morpho:bool=True,

                                # regions settings
                                include_regions:bool=True,

                                # interaction settings
                                include_interactions:bool=True,
                                int_splitter:str="X",
                                include_inter_morpho:bool=True,
                                include_inter_degrees:bool=True,
                                export_inter_degree_imgs:bool=True,
                                export_interaction_sites:bool=True,
                                include_inter_dist:bool=True,

                                # distribution settings
                                include_org_dist:bool=True, 
                                dist_centering_obj: Union[str, None]=None,
                                dist_num_bins: Union[int, None]=5,
                                dist_center_on: Union[bool, None]=False,
                                dist_keep_center_as_bin: Union[bool, None]=True,
                                dist_zernike_degrees: Union[int, None]=9,
                                export_dist_bins_imgs: bool = True):
    """
    Batch process interaction quantification for a single dataset (e.g., images collected on the same data). 
    Morphology, distribution, and degree of interaction metrics analysis are all optionally available. 
    Interaction site segmentations and degree of interaction images can also be exported.
    
    Parameters:
    -----------
    dataset_name : str
        A unique string identifier for the dataset being processed. It will be included as metadata in output tables and as 
        part of the output files names. It will also be used to identify if any data that has already been collected for this dataset.
    raw_path : Union[Path,str]
        Path or str to the folder that contains the raw image files.
    seg_path : Union[Path,str]
        Path or str to the folder that contains the segmentation tiff files.
    quant_path : Union[Path, str]
        Path or str to the folder that the output datatables will be saved to.
    raw_file_type : str
        File type of the raw images (e.g., "czi", "tiff")
    channel_axis : int
        Axis corresponding to the channels in the image data (e.g., 0 for (C,Z,Y,X) data, 3 for (Z,Y,X,C) data)
    organelle_names : List[str]
        List of organelle names to analyze. These names should match the suffixes used to name the organelle segmentation files.
        (e.g., "mito" for mito segmentation files named "imageX-mito.tiff")
    organelle_channels : Union[List[int], None], default=None
        List of intensity channel indices from the raw files; the order should correspond to each organelle name included in organelle_names 
        parameter. If no intensity analysis is to be included, specify None here. If not specified, the default is None.
    region_names : Union[List[str], None], default=None
        List of region names to analyze. These names should match the suffixes used to name the organelle segmentation files.
        Usually, the input will be ['cell', 'nuc'] for cell mask and nucleus regions. If no regions are to be included, specify None here.
        If not specified, the default is None.
    mask_name : Union[str, None], default=None
        Name of the region to use as a mask on all images/segmentations before quantification. This name should be included in the 
        regions_name variable. If None, the entire image will be quantified. If not specified, the default is None.
    use_scale : bool, default=True
        Whether to apply scaling to the quantitative data; scaled data will be in real world units (e.g., microns) rather than pixels/voxels.
        The scale is determined from the metadata of the raw image files. If no scale information is found, all data will be in pixel/voxel units.
    seg_suffix : Union[str, None], default=None
        Any additional text that is included in the segmentation tiff files between the file stem and the segmentation suffix, 
        not including the initial dash ("-") separating the file name from the suffix. (e.g., "20230426_test-" for segmentation files named 
        "imageX-20230426_test-mito.tiff")
    include_org_morph : bool, default=True
        Whether to compute morphology metrics for each organelle object. If not specified, the default is True.
    include_regions : bool, default=True
        Whether to compute morphology metrics for each provided region. If not specified, the default is True.
    include_interactions : bool, default=True
        Whether to compute interaction site analyses. If not specified, the default is True. If `include_inter_morpho`, `include_interaction_degrees`,
        and/or `include_dist` are all False, only interaction site labels will be exported.
    int_splitter: str, default="X"
        Character used to separate organelles within the interaction site names (e.g., "X" in "mitoXlyso", the name for mitochondria and lysosome 
        interaction sites).
    include_inter_morpho : bool, default=True
        Whether to compute morphology metrics for each interaction site. If not specified, the default is True.
    include_inter_degrees : bool, default=True
        Whether to compute interaction degree analysis. If not specified, the default is True.
    export_inter_degree_imgs : bool, default=True
        Whether to export interaction degree images. If not specified, the default is True.
    export_interaction_sites : bool, default=True
        Whether to export interaction site images (including interaction site objects across the entire image; not masked). If not specified, the 
        default is True.
    include_inter_dist : bool, default=True
        Whether to compute interaction distribution metrics. If not specified, the default is True.
    include_org_dist : bool, default=True
        Whether to compute organelle distribution metrics. If not specified, the default is True.
    dist_centering_obj : str or None, default=None
        Name of the region to use to find the center of the XY region in the distribution analysis. This region should be included in the 
        list_region_names and list_region_segs variables. If None, the center of the mask, or entire image if no mask was specified, will be used as 
        the centering object. If not specified, the default is None. 
        This setting is used in both the interaction site and organelle distribution 
        analysis.
    dist_num_bins : int or None, default=5
        Number of radial bins to create in the XY distribution analysis. None is only allowed if include_inter_dist=False and include_org_dist=False. 
        This setting is used in both the interaction site and organelle distribution analysis.
    dist_center_on : bool or None, default=True
        Whether to start creation of the XY distribution bins from the center of the centering object (True) or edge (False). None is only allowed if 
        include_inter_dist=False and include_org_dist=False. 
        This setting is used in both the interaction site and organelle distribution analysis.
    dist_keep_center_as_bin : bool or None
        Whether to keep centering object as the first XY bin. None is only allowed if include_inter_dist=False and include_org_dist=False. 
        This setting is used in both the interaction site and organelle distribution analysis.
    dist_zernike_degrees : int or None, default=9
        Zernike polynomial degree for shape analysis in the XY distribution analysis. If None and include_inter_dist=True or include_org_dist=True, 
        no Zernike features will be calculated. 
        This setting is used in both the interaction site and organelle distribution analysis.
    export_dist_bins_imgs : bool, default=True
        Whether to export images of the distribution bins used in the distribution analysis. These images will be exported into a new directory: 
        quant_path / f"{dataset_name}-distribution_bins_imgs". 
        If not specified, the default is True.

    Returns:
    --------
    None
        Saves output files to the specified quantification path
    """
    #####
    # start timing & count for number of images processed
    batch_start = time.time()
    count = 0


    #####
    # check (and update if necessary) format of file paths
    if isinstance(raw_path, str): raw_path = Path(raw_path)
    if isinstance(seg_path, str): seg_path = Path(seg_path)
    if isinstance(quant_path, str): quant_path = Path(quant_path)

    # create output directory if it does not already exist
    if not Path.exists(quant_path):
        Path.mkdir(quant_path)
        print(f"Output file path not found. Making {quant_path}.")

    # list files from the raw path that will be processed below
    img_file_list = list_image_files(raw_path, raw_file_type)

    # print warning and stop if no files found
    len_file_list = len(img_file_list)
    if len_file_list==0:
        raise ValueError(f"No '{raw_file_type}' files found in {raw_path}. Please check the input path and file type.")

    # checking for any analyses
    if not any([include_org_morpho, include_regions, include_interactions, include_org_dist]):
        raise ValueError("No analyses selected. Update the \"include_*\" parameters to select at least one analysis.")

    #####
    # define output file paths & check if any existing data is present for this dataset
    # define unique keys in the output file for checking existing data
    unique_keys = ['dataset', 'image_name']

    # Build analysis configuration and load existing keys based on include_* parameters
    analyses_config = {
        'org_morpho': (include_org_morpho, f"{dataset_name}-organelle_morphology_metrics.csv"),
        'regions': (include_regions, f"{dataset_name}-regions_morphology_metrics.csv"),
        'inter_morpho': (include_interactions and include_inter_morpho, f"{dataset_name}-interactions_morphology_metrics.csv"),
        'inter_labs': (include_interactions and not include_inter_morpho, f"{dataset_name}-interactions_labels.csv"),
        'inter_degrees': (include_interactions and include_inter_degrees, f"{dataset_name}-interactions_degree_metrics.csv"),
        'inter_dist': (include_interactions and include_inter_dist, f"{dataset_name}-interactions_distribution_metrics.csv"),
        'org_dist': (include_org_dist, f"{dataset_name}-organelle_distribution_metrics.csv"),
    }

    # Load existing keys for each enabled analysis
    existing_keys_dict = {}
    for analysis_name, (is_included, filename) in analyses_config.items():
        if is_included:
            existing_keys_dict[analysis_name] = load_existing_keys_csv(quant_path / filename, unique_keys)
        else:
            existing_keys_dict[analysis_name] = None

    # Compute intersection: only skip if file is complete in ALL requested analyses
    relevant_keys = [keys for keys in existing_keys_dict.values() if keys is not None]
    existing_keys = set.intersection(*relevant_keys) if relevant_keys else set()

    # Define specific existing keys for interaction analyses
    inter_subset = {atype: existing_keys_dict[atype] for atype in ['inter_morpho', 'inter_labs', 'inter_degrees', 'inter_dist'] if analyses_config[atype][0]}
    inter_relevant_keys = [keys for keys in inter_subset.values() if keys]
    existing_inter_keys = set.intersection(*inter_relevant_keys) if inter_relevant_keys else set()

    # Define export paths for interactions outputs
    int_degree_img_path = quant_path / f"{dataset_name}-interaction_degree_images" if export_inter_degree_imgs else None
    interaction_sites_path = quant_path / f"{dataset_name}-interaction_site_segmentations" if export_interaction_sites else None
    dist_bins_path = quant_path / f"{dataset_name}-distribution_bins_imgs" if export_dist_bins_imgs else None

    #####
    # create a combined list segmentation images to collect for each image
    # if no masks are provided, only organelle segmentations will be collected
    segs_to_collect = organelle_names + region_names if region_names is not None else organelle_names

    # loop through list of images; analyze each and save the data the csv file
    for img_f in img_file_list:
        # skip files that have already been processed based on existing keys
        if (dataset_name, img_f.stem) in existing_keys:
            print(f"Skipping {img_f.name} as it is already listed in ALL of the output file(s).")
            continue

        # if not in existing keys, process analysis for this cell
        else:
            # start second timer for each image & increment count of images processed
            img_start = time.time()
            count = count+1
            
            # find file paths for raw and all segmentation files
            filez = find_segmentation_tiff_files(img_f, segs_to_collect, seg_path, seg_suffix)

            # read in raw image and metadata
            img_data, meta_dict = read_czi_image(filez["raw"])

            # format intensity channel separately from the raw file based on provided organelle channels input information
            # if not organelle channels provided, set to None (no intensity analysis will be performed)
            if organelle_channels is None:
                intensities = None
                print("No organelle channels provided; intensity analysis will be skipped.")
            else:
                if channel_axis != 0:
                    img_data = np.moveaxis(img_data, channel_axis, 0)
                intensities = [img_data[ch] for ch in organelle_channels]

            # store organelle segmentation images as list
            if organelle_names is None:
                raise ValueError("No organelle names provided; organelle segmentations are required for quantification.")
            else:
                organelles = [read_tiff_image(filez[org]) for org in organelle_names]

            # store regions segmentations as a list
            if region_names is None:
                regions = None
                mask_name = None  # if no regions provided, no mask can be applied
                print("No region names provided; regions analysis will not be calculated, no mask will be applied, and no centering object will be used in distribution analysis if specified.")
            else:
                regions = [read_tiff_image(filez[r]) for r in region_names]

            # define the scale for quantification
            if use_scale is True:
                try:
                    assert 'scale' in meta_dict.keys(), "No scale information found in image metadata."
                except AssertionError as e:
                    print("No scale information found in image metadata; all measurements will be in pixel/voxel units.")
                    scale = None
                else:
                    scale = meta_dict['scale']
                    print(f"Using pixel/voxel scale: {scale}")
            else:
                scale = None
                print(f"No scale was provided; all measurements will be in pixel/voxel units.")


            # process organelle morphology analysis, if specified
            if include_org_morpho:
                # skip files that have already been processed for this analysis
                if (dataset_name, img_f.stem) in existing_keys_dict['org_morpho']:
                    print(f"Skipping organelle morphology analysis for {img_f.name} as it is already listed in the output file(s).")
                else:
                    # get organelle morphology metrics for all included organelles
                    org_morph_metrics = get_org_morphology(source_file_path=img_f,
                                                            list_obj_names=organelle_names,
                                                            list_obj_segs=organelles,
                                                            list_intensity_img=intensities, 
                                                            list_region_names=region_names,
                                                            list_region_segs=regions, 
                                                            mask_name=mask_name,
                                                            scale=scale)
                
                    # add dataset name column to the organelle morphology data
                    org_morph_metrics.insert(loc=0,column='dataset',value=dataset_name)

                    # save the organelle morphology data for this image directly to csv
                    append_atomic_csv(quant_path / analyses_config['org_morpho'][1], org_morph_metrics)
                    del org_morph_metrics  # delete table variable to free up memory

            # process regions morphology analysis, if specified
            if include_regions and (region_names is not None):
                # skip files that have already been processed for this analysis
                if (dataset_name, img_f.stem) in existing_keys_dict['regions']:
                    print(f"Skipping regions analysis for {img_f.name} as it is already listed in the output file(s).")
                else:
                    # get region morphology metrics for all included regions
                    # intensity analysis will only be included for channels that correspond to organelles (as specified in the organelle_channels input)
                    regions_metrics = get_regions_morphology(source_file_path=img_f,
                                                            list_region_names=region_names,
                                                            list_region_segs=regions,
                                                            list_intensity_img=intensities,
                                                            list_channel_names=organelle_names,
                                                            mask_name=mask_name,
                                                            scale=scale)
                    
                    # add dataset name to regions table and save table directly to csv
                    regions_metrics.insert(loc=0,column='dataset',value=dataset_name)
                    append_atomic_csv(quant_path / analyses_config['regions'][1], regions_metrics)
                    del regions_metrics  # delete table to free up memory

            # process organelle interactions analysis, if specified
            # the same distribution analysis settings used for organelle distribution analysis will be applied to interaction site distribution analysis, if interactions distribution anlaysis included
            if include_interactions:
                # skip files that have already been processed for all of the interactions analysis
                if (dataset_name, img_f.stem) in existing_inter_keys:
                    print(f"Skipping interactions analysis for {img_f.name} as it is already listed in all of the interactions output file(s).")
                else:
                    inter_dict, inter_morph_tab, inter_dist_tab, int_degree_img, int_degree_tab, XY_bins, XY_wedges = get_interaction_metrics(source_file_path=img_f,
                                                                                                                        list_obj_names=organelle_names,
                                                                                                                        list_obj_segs=organelles,
                                                                                                                        list_intensity_img=intensities,
                                                                                                                        list_region_names=region_names,
                                                                                                                        list_region_segs=regions,
                                                                                                                        mask_name=mask_name,
                                                                                                                        splitter=int_splitter,
                                                                                                                        scale=scale,
                                                                                                                        include_morpho=include_inter_morpho,
                                                                                                                        include_inter_degrees=include_inter_degrees,
                                                                                                                        include_dist=include_inter_dist, 
                                                                                                                        dist_centering_obj=dist_centering_obj,
                                                                                                                        dist_num_bins=dist_num_bins,
                                                                                                                        dist_center_on=dist_center_on,
                                                                                                                        dist_keep_center_as_bin=dist_keep_center_as_bin,
                                                                                                                        dist_zernike_degrees=dist_zernike_degrees)

                    # if interactions morphology included, skip files that have already been processed for this analysis
                    if include_inter_morpho:
                        if (dataset_name, img_f.stem) in existing_keys_dict['inter_morpho']:
                            print(f"Skipping interactions morphology analysis for {img_f.name} as it is already listed in the output file(s).")
                        else:
                            # save interaction morphology table data per image directly to csv
                            inter_morph_tab.insert(loc=0,column='dataset',value=dataset_name)
                            append_atomic_csv(quant_path / analyses_config['inter_morpho'][1], inter_morph_tab)
                            del inter_morph_tab  # free up memory
                    # if interactions morphology, not included, interactions labels are; check for files then save if necessary
                    else:
                        if (dataset_name, img_f.stem) in existing_keys_dict['inter_labs']:
                            print(f"Skipping interactions labels analysis for {img_f.name} as it is already listed in the output file(s).")
                        else:
                            # save interaction labels table data per image directly to csv
                            inter_morph_tab.insert(loc=0,column='dataset',value=dataset_name)
                            append_atomic_csv(quant_path / analyses_config['inter_labs'][1], inter_morph_tab)
                            del inter_morph_tab  # free up memory

                    # if specified, export interaction site segmentation images
                    if export_interaction_sites:
                        inter_site_cnt=0
                        for inter_name, inter_img in inter_dict.items():
                            inter_site_cnt+=1
                            # check that the file does not already exist before exporting
                            if not (Path(interaction_sites_path)/f"{img_f.name}-{inter_name}.tiff").exists():
                                export_inferred_organelle(inter_img, f"{inter_name}", meta_dict, interaction_sites_path)
                            else:
                                if inter_site_cnt<=1:
                                    warnings.warn(f"Some of the interaction site images already exist for {img_f.name} in {interaction_sites_path}. They will not be overwritten.", UserWarning)
                    del inter_dict  # free up memory
                
                    # save the distribution table data per image directly to csv
                    if include_inter_dist:
                        # skip files that have already been processed for this analysis
                        if (dataset_name, img_f.stem) in existing_keys_dict['inter_dist']:
                            print(f"Skipping interactions distribution analysis for {img_f.name} as it is already listed in the output file(s).")
                        else:
                            # TODO: remove .astype(str) once method distribution functions have been updated
                            inter_dist_tab = inter_dist_tab.astype(str)  # ensure all data is string to avoid dtype issues
                            inter_dist_tab.insert(loc=0,column='dataset',value=dataset_name)
                            append_atomic_csv(quant_path / analyses_config['inter_dist'][1], inter_dist_tab)

                            if export_dist_bins_imgs:
                                # export XY bins and wedges as images
                                if not Path(dist_bins_path / f"{img_f.stem}-XY_bins.tiff").exists():
                                    export_inferred_organelle(XY_bins.astype(np.uint16), "XY_bins", meta_dict, dist_bins_path)
                                else:
                                    warnings.warn(f"The XY distribution bins images already exist for {img_f.stem} in {dist_bins_path}. They will not be overwritten.", UserWarning)

                                if not Path(dist_bins_path / f"{img_f.stem}-XY_wedges.tiff").exists():
                                    export_inferred_organelle(XY_wedges.astype(np.uint16), "XY_wedges", meta_dict, dist_bins_path)
                                else:
                                    warnings.warn(f"The XY distribution wedges images already exist for {img_f.stem} in {dist_bins_path}. They will not be overwritten.", UserWarning)
                    del inter_dist_tab  # free up memory
                    del XY_bins  # free up memory
                    del XY_wedges  # free up memory

                    # save the degree table data per image directly to csv
                    if include_inter_degrees:
                        # skip files that have already been processed for this analysis
                        if (dataset_name, img_f.stem) in existing_keys_dict['inter_degrees']:
                            print(f"Skipping interactions degree analysis for {img_f.name} as it is already listed in the output file(s).")
                        else:
                            int_degree_tab.insert(loc=0,column='dataset',value=dataset_name)
                            append_atomic_csv(quant_path / analyses_config['inter_degrees'][1], int_degree_tab)

                        # save the degree image
                        if export_inter_degree_imgs:
                            if not (Path(int_degree_img_path)/f"{img_f.name}-interactions_degree.tiff").exists():
                                export_inferred_organelle(int_degree_img.astype(np.uint16), "interactions_degree", meta_dict, int_degree_img_path)
                            else:
                                warnings.warn(f"The {img_f.name}-interactions_degree.tiff image already exists in {int_degree_img_path}. It will not be overwritten.")
                    del int_degree_tab  # free up memory
                    del int_degree_img  # free up memory

            # process organelle distribution analysis, if specified
            if include_org_dist:
                if (dataset_name, img_f.stem) in existing_keys_dict['org_dist']:
                    print(f"Skipping organelle distribution analysis for {img_f.name} as it is already listed in the output file(s).")
                else:
                    dist_tab, XY_bins_img, XY_wedges_img = get_distribution_metrics(source_file_path=img_f,
                                                                                    list_obj_names=organelle_names,
                                                                                    list_obj_segs=organelles, 
                                                                                    list_region_names=region_names,
                                                                                    list_region_segs=regions, 
                                                                                    mask_name=mask_name,
                                                                                    scale=scale,
                                                                                    centering_obj=dist_centering_obj,
                                                                                    num_bins=dist_num_bins,
                                                                                    center_on=dist_center_on,
                                                                                    keep_center_as_bin=dist_keep_center_as_bin,
                                                                                    zernike_degrees=dist_zernike_degrees)
                    
                    dist_tab = dist_tab.astype(str)  # ensure all data is string to avoid dtype issues

                    # save the distribution (or labels only) table data per image directly to csv
                    dist_tab.insert(loc=0,column='dataset',value=dataset_name)
                    append_atomic_csv(quant_path / analyses_config['org_dist'][1], dist_tab)
                    del dist_tab  # free up memory

                    if export_dist_bins_imgs and not include_inter_dist:  # only export distribution bins images here if interaction distribution analysis is not included, as the same distribution bins images will be used for both analyses
                        # export XY bins and wedges as images
                        if not Path(dist_bins_path / f"{img_f.stem}-XY_bins.tiff").exists():
                            export_inferred_organelle(XY_bins_img.astype(np.uint16), "XY_bins", meta_dict, dist_bins_path)
                        else:
                            warnings.warn(f"The XY distribution bins images already exist for {img_f.stem} in {dist_bins_path}. They will not be overwritten.", UserWarning)

                        if not Path(dist_bins_path / f"{img_f.stem}-XY_wedges.tiff").exists():
                            export_inferred_organelle(XY_wedges_img.astype(np.uint16), "XY_wedges", meta_dict, dist_bins_path)
                        else:
                            warnings.warn(f"The XY distribution wedges images already exist for {img_f.stem} in {dist_bins_path}. They will not be overwritten.", UserWarning)

            # end timer for single image
            end2 = time.time()
            print(f"Completed quantification of {img_f.name} in {(end2-img_start)/60} mins.")
            print(f"{count}/{len_file_list} images have been processed. \n")

    # end timer for entire batch
    batch_end = time.time()
    print(f"Quantification for {count} files is COMPLETE in {(batch_end - batch_start)/60} minutes! Files saved to '{quant_path}'.")



# summary function to combine and summarize the data collected across all images in the batch_process_quantification function
def batch_process_summarystats(out_prefix: str,
                               csv_path_list: List[str],
                               out_path: str,
                               organelle_names: List[str],
                               region_names: Union[List[str], None]=None,
                               mask_name: Union[str, None]=None,
                               splitter: Union[str, None] = "X"):
    """ 
    Batch process quantification summary statistics for all possible quantification types from multiple datasets.

    Parameters:
    -----------
    out_prefix: str
        The prefix used to name the output file. A dash ("-") will be included between this prefix and the output file name. (e.g., 
        "experiment1" will result in output files named "experiment1-organelle_morphology_summary.csv", etc.)
    csv_path_list: List[str],
        A list of path strings where .csv files to analyze are located. The csv files in each location should be a result of the 
        batch_process_quantification function. The function will search through all csv files in each location to find the relevant
        quantification data files. If there is more than one dataset in a location, all datasets will be included in the summary statistics.
    out_path: str,
        The path to the location where the summary data files will be output.
    organelle_names: List[str],
        A list of organelle names used in the interaction quantification analysis. This list should match the organelle names used in the
        batch_process_quantification function and should be consistent across all datasets being analyzed per run.
    region_names: Union[List[str], None]=None,
        A list of region names used in the interaction quantification analysis. This list should match the organelle names used in the
        batch_process_quantification function and should be consistent across all datasets being analyzed per run. If no regions were included in 
        the quantification analysis, specify None here. If not specified, the default is None.
    mask_name: Union[str, None]=None,
        The name of the mask used in the quantification analysis. If no mask was used, specify None here. If not specified, the default is None.
    splitter: Union[str, None], default="X"
        The character used to split interaction site names in the batch_process_quantification function. This should be consistent across all 
        datasets being analyzed per run. If no interaction sites were included in the quantification analysis, you can specify None here or ignore this parameter. 
        If not specified, the default is "X".
    """
    
    # validate inputs
    if not isinstance(out_prefix, str): 
        raise ValueError("Output prefix must be a string.")
    
    if not Path(out_path).exists():
        raise ValueError(f"Output path does not exist: {out_path}")
    
    for path in csv_path_list:
        if not Path(path).exists():
            raise ValueError(f"CSV path does not exist: {path}")
        
    if not isinstance(organelle_names, list) or not all(isinstance(name, str) for name in organelle_names):
        raise ValueError("Organelle names must be provided as a list of strings.")
    
    if region_names is not None:
        if not isinstance(region_names, list) or not all(isinstance(name, str) for name in region_names):
            raise ValueError("Region names must be provided as a list of strings.")



    ###############################################################
    # Read in the csv files and combine them into one of each type
    ###############################################################
    # for keeping track of dataset and file numbers
    ds_count = 0
    fl_count = 0

    # create empty list to hold the morphology tables from different experiments
    data_frames = {'org_morph': [],
                   'regions': [],
                   'int_morph': [],
                   'int_labs': [],
                   'int_dist': [],
                   'int_degree': [],
                   'org_dist': []}

    # name patterns to look for in the csv files
    metric_mapping = {"-organelle_morphology_metrics.csv": 'org_morph',
                      "-regions_morphology_metrics.csv": 'regions',
                      "-interactions_morphology_metrics.csv": 'int_morph',
                      "-interactions_labels.csv": 'int_labs',
                      "-interactions_distribution_metrics.csv": 'int_dist',
                      "-interactions_degree_metrics.csv": 'int_degree',
                      "-organelle_distribution_metrics.csv": 'org_dist'}
    metric_suffixes = tuple(metric_mapping.keys())
    

    # loop through all of the locations listed above and find the _org_morph files; append them to the list above
    for loc in csv_path_list:
        # list all csv files in the location
        files_store = sorted(loc.glob("*.csv"))

        # find the unique datasets in this location based on the dataset names
        prefixes = {f.name.rsplit("-", 1)[0] for f in files_store if f.name.endswith(metric_suffixes)}
        print(f"Found the following datasets in {loc}: {prefixes}")

        for prefix in prefixes:
            ds_count += 1
            files_subset = [f for f in files_store if f.name.startswith(tuple([prefix + item for item in metric_suffixes]))]


            # if both morphology and labels files are present, remove the labels file from the list to be processed
            if any("-interactions_morphology_metrics.csv" in f.name for f in files_subset) and any("-interactions_labels.csv" in f.name for f in files_subset):
                files_subset = [f for f in files_subset if "-interactions_labels.csv" not in f.name]

            for file in files_subset:
                fl_count += 1
                stem = file.name
                
                # determine which type of metric this file contains
                metric_type = next((key for key in metric_mapping if key in stem), None)

                # read in the file and append to the appropriate list
                if metric_type:
                    data_frames[metric_mapping[metric_type]].append(pd.read_csv(file))
                else:
                    print(f"File {stem} not recognized; skipping.")

    print(f"Found {fl_count} files from {ds_count} dataset(s) across {len(csv_path_list)} location(s).")

    # combine the org_morph lists found above into one combined table with all data
    org_morpho_df = pd.concat(data_frames['org_morph'], axis=0, join='outer', ignore_index=True, copy=False) if data_frames['org_morph'] else None
    regions_df = pd.concat(data_frames['regions'], axis=0, join='outer', ignore_index=True, copy=False) if data_frames['regions'] else None
    inter_labs_df = pd.concat(data_frames['int_labs'], axis=0, join='outer', ignore_index=True, copy=False) if data_frames['int_labs'] else None
    inter_morph_df = pd.concat(data_frames['int_morph'], axis=0, join='outer', ignore_index=True, copy=False) if data_frames['int_morph'] else None
    inter_degree_df = pd.concat(data_frames['int_degree'], axis=0, join='outer', ignore_index=True, copy=False) if data_frames['int_degree'] else None
    inter_dist_df = pd.concat(data_frames['int_dist'], axis=0, join='outer', ignore_index=True, copy=False) if data_frames['int_dist'] else None
    org_dist_df = pd.concat(data_frames['org_dist'], axis=0, join='outer', ignore_index=True, copy=False) if data_frames['org_dist'] else None


    ## TODO: break the following down into subfunctions for each analysis type - update respective notebooks/individaul functions ot reflect

    #############################################
    # Summarize organelle count & morphology data
    #############################################
    if org_morpho_df is not None:
        # summary stat group
        org_group_by = ['dataset', 'image_name', 'mask_name', 'scale', 'object']
        org_sharedcolumns = ["SA_to_volume_ratio", "equivalent_diameter", "extent", "euler_number", "solidity", "axis_major_length"] + list(org_morpho_df.filter(regex=".*intensity.*").columns)
        org_ag_func_standard = ['mean', 'median', 'std']

        # summarize shared measurements between org_df and contacts_df
        org_tab1 = org_morpho_df[org_group_by + ['label']].groupby(org_group_by).agg(['count'])
        org_tab1.rename(columns={'label': 'org'}, inplace=True)
        org_tab2 = org_morpho_df[org_group_by + ['volume', 'surface_area']].groupby(org_group_by).agg(['sum'] + org_ag_func_standard)
        org_tab3 = org_morpho_df[org_group_by + org_sharedcolumns].groupby(org_group_by).agg(org_ag_func_standard)
        org_summary = pd.merge(org_tab1, org_tab2, 'outer', on=org_group_by)
        org_summary = pd.merge(org_summary, org_tab3, 'outer', on=org_group_by)

        # Get mask_name and corresponding volume column per group & calculate volume fraction
        org_mask_names = org_morpho_df.groupby(org_group_by)['mask_name'].first()
        org_mask_volume_data = org_morpho_df.groupby(org_group_by).first().apply(lambda row: row[f"{org_mask_names.loc[row.name]}_volume"], axis=1)
        org_summary.insert(org_summary.columns.get_loc(('volume', 'sum')) + 1, ('volume', 'fraction'), org_summary[('volume', 'sum')]/org_mask_volume_data)

        # fill gaps & NA values
        # Ensure all possible interactions are represented (if missing fill with NaN)
        for ind in org_summary.index.droplevel(4).unique().to_list():
            for row in organelle_names:
                if ind+(row,) not in org_summary.index:
                    org_summary.loc[ind+(row,)] = np.nan

        # fill NA with 0 for specific columns
        org_fill_dict = {('org', 'count'): 0, 
                     ('volume', 'sum'): 0,
                     ('surface_area', 'sum'): 0,
                     ('volume', 'fraction'): 0}
        org_summary = org_summary.fillna(value=org_fill_dict)

        # if (org, count) is 1, set mean, median, and std to NaN
        org_single_site_mask = org_summary[('org', 'count')] == 1
        for col in org_sharedcolumns+['volume', 'surface_area']:
            org_summary.loc[org_single_site_mask, (col, 'std')] = np.nan

        org_summary.sort_index(inplace=True)

        # flatten datasheet and export
        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-organelle_morphology_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-organelle_morphology_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            org_summary.to_csv(str(out_path) + f"/{out_prefix}-organelle_morphology_summarystats.csv", mode='x')
            print(f"Exported per-organelle morphology summary statistics (before unstacking) to {out_path}/{out_prefix}-organelle_morphology_summarystats.csv")
        org_morph_final = org_summary.unstack(-1)
        org_morph_final.columns = ["_".join((col_name[1], col_name[-1], col_name[0])) for col_name in org_morph_final.columns.to_flat_index()]
        org_morph_final.columns = [col.replace('sum', 'total') for col in org_morph_final.columns]
        org_morph_final.columns = [col.replace(col, 'mask_volume') if 'mask' in col else col for col in org_morph_final.columns]
        org_morph_final = org_morph_final.loc[:, ~org_morph_final.columns.duplicated()]
        # org_morph_final.reset_index(inplace=True)
        # org_morph_final.set_index(['dataset', 'image_name', 'mask_name', 'scale'], inplace=True)

        final_combo_tab = org_morph_final
    else:
        final_combo_tab = pd.DataFrame()


    ####################################
    # Summarize regions morphology data
    ####################################
    if regions_df is not None:
        # summary stat group
        reg_group_by = ['dataset', 'image_name', 'mask_name', 'scale', 'object']
        reg_sharedcolumns = ["SA_to_volume_ratio", "equivalent_diameter", "extent", "euler_number", "solidity", "axis_major_length"] + list(regions_df.filter(regex=".*intensity.*").columns)
        reg_ag_func_standard = ['mean', 'median', 'std']

        # summarize morphology metrics
        reg_tab1 = regions_df[reg_group_by + ['label']].groupby(reg_group_by).agg(['count'])
        reg_tab1.rename(columns={'label': 'region'}, inplace=True)
        reg_tab2 = regions_df[reg_group_by + ['volume', 'surface_area']].groupby(reg_group_by).agg(['sum'] + reg_ag_func_standard)
        reg_tab3 = regions_df[reg_group_by + reg_sharedcolumns].groupby(reg_group_by).agg(reg_ag_func_standard)
        regions_summary = pd.merge(reg_tab1, reg_tab2, 'outer', on=reg_group_by)
        regions_summary = pd.merge(regions_summary, reg_tab3, 'outer', on=reg_group_by)

        # Get mask_name and corresponding volume column per group & calculate volume fraction
        reg_mask_names = regions_df.groupby(reg_group_by)['mask_name'].first()
        reg_mask_volume_data = regions_df.groupby(reg_group_by).first().apply(lambda row: row[f"{reg_mask_names.loc[row.name]}_volume"], axis=1)
        regions_summary.insert(regions_summary.columns.get_loc(('volume', 'sum')) + 1, ('volume', 'fraction'), regions_summary[('volume', 'sum')]/reg_mask_volume_data)

        # fill gaps & NA values
        # Ensure all possible interactions are represented (if missing fill with NaN)
        for ind in regions_summary.index.droplevel(4).unique().to_list():
            for row in region_names:
                if ind+(row,) not in regions_summary.index:
                    regions_summary.loc[ind+(row,)] = np.nan

        # fill NA with 0 for specific columns
        reg_fill_dict = {('region', 'count'): 0, 
                    ('volume', 'sum'): 0,
                    ('surface_area', 'sum'): 0,
                    ('volume', 'fraction'): 0}
        regions_summary = regions_summary.fillna(value=reg_fill_dict)

        # if (region, count) is 1, set mean, median, and std to NaN
        reg_single_site_mask = regions_summary[('region', 'count')] == 1
        for col in reg_sharedcolumns+['volume', 'surface_area']:
            regions_summary.loc[reg_single_site_mask, (col, 'std')] = np.nan

        regions_summary.sort_index(inplace=True)

        # flatten datasheet and export
        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-region_morphology_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-region_morphology_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            regions_summary.to_csv(str(out_path) + f"/{out_prefix}-region_morphology_summarystats.csv", mode='x')
            print(f"Exported per-region morphology summary statistics (before unstacking) to {out_path}/{out_prefix}-region_morphology_summarystats.csv")
        regions_morph_final = regions_summary.unstack(-1)
        regions_morph_final.columns = ["_".join((col_name[1], col_name[-1], col_name[0])) for col_name in regions_morph_final.columns.to_flat_index()]
        regions_morph_final.columns = [col.replace('sum', 'total') for col in regions_morph_final.columns]
        regions_morph_final.columns = [col.replace(col, 'mask_volume') if 'mask' in col else col for col in regions_morph_final.columns]
        regions_morph_final = regions_morph_final.loc[:, ~regions_morph_final.columns.duplicated()]
        # regions_morph_final.reset_index()
        # regions_morph_final.set_index(['dataset', 'image_name', 'mask_name', 'scale'], inplace=True)
  
        final_combo_tab = pd.concat([final_combo_tab, regions_morph_final], axis=1)
    else:
        final_combo_tab = final_combo_tab


    ################################################
    # Summarize interactions count & morphology data
    ################################################
    # list all possible interaction site combinations
    all_pos = all_combos(organelle_names, splitter)

    if inter_morph_df is not None:
        ### calculate interaction count/volume & summarize per organelle object for all interaction sites
        inter_per_org_summary = perorg_interactions_cnt(interaction_morpho_df=inter_morph_df, 
                                                   org_list=organelle_names,
                                                   splitter=splitter)

        # summarization parameters
        inter_count_vol_group_by = ["dataset", "image_name", "mask_name", "scale", "object"]
        inter_count_vol_cols = [col for col in inter_per_org_summary.columns if col.endswith(("_count", "_volume"))]
        inter_count_vol_ag_func_standard = {"num_interaction_types": ['mean', 'median', 'std']} | {col: ['sum', 'mean', 'median', 'std'] for col in inter_count_vol_cols}

        # summarize per organelle type per image
        inter_org_sum_tab = inter_per_org_summary.groupby(inter_count_vol_group_by).agg(inter_count_vol_ag_func_standard)
    
        # Ensure all possible interactions are represented (if missing fill with NaN)
        for ind in inter_org_sum_tab.index.droplevel(4).unique().to_list():
            for row in organelle_names:
                if ind+(row,) not in inter_org_sum_tab.index:
                    inter_org_sum_tab.loc[ind+(row,)] = np.nan
        inter_org_sum_tab.sort_index(inplace=True)

        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-interaction_count_volume_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-interaction_count_volume_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            inter_org_sum_tab.to_csv(str(out_path) + f"/{out_prefix}-interaction_count_volume_summarystats.csv", mode='x')
            print(f"Exported per-organelle interaction count/volume summary statistics (before unstacking) to {out_path}/{out_prefix}-interaction_count_volume_summarystats.csv")
        # unstack and format interaction count/volume summary table
        inter_count_vol_final = inter_org_sum_tab.unstack(-1)
        for col in inter_count_vol_final.columns:
            if col[0].endswith(('_count', '_volume')):
                if col[2] not in col[0]:
                    inter_count_vol_final.drop(col,axis=1, inplace=True)

        inter_count_vol_final.columns = ["_".join((col_name[1], col_name[0], "per", col_name[-1])) for col_name in inter_count_vol_final.columns.to_flat_index()]
        inter_count_vol_final.columns = [col.replace('sum', 'total') for col in inter_count_vol_final.columns]
        inter_count_vol_final.columns = [col.replace('per', 'in') if 'total' in col else col for col in inter_count_vol_final.columns]
        inter_count_vol_final.fillna(0, inplace=True)
        # inter_count_vol_final.reset_index()


        ### summarize interaction morphology per interaction site
        # summarization paramters
        inter_group_by = ["dataset", "image_name", "mask_name", "scale", "object"]
        inter_cols = ["SA_to_volume_ratio", "equivalent_diameter", "extent", "euler_number", "solidity", "axis_major_length"] + list(inter_morph_df.filter(regex=".*intensity.*").columns)
        inter_ag_func_standard = ['mean', 'median', 'std']

        # summarize counts of interaction sites per image
        tab1 = inter_morph_df[inter_group_by + ['ID']].groupby(inter_group_by).agg(['count'])
        tab1.rename(columns={'ID': 'sites'}, inplace=True)
        tab2 = inter_morph_df.copy()[inter_morph_df['in_higher_order'] == True][inter_group_by + ['ID']].groupby(inter_group_by).agg(['count'])
        tab2.rename(columns={'ID': 'sites_in_higher_order'}, inplace=True)
        tab3 = inter_morph_df.copy()[inter_morph_df['in_higher_order'] == False][inter_group_by + ['ID']].groupby(inter_group_by).agg(['count'])
        tab3.rename(columns={'ID': 'sites_not_in_higher_order'}, inplace=True)
        inter_sum_tab = pd.merge(tab1, tab2, 'outer', on=inter_group_by)
        inter_sum_tab = pd.merge(inter_sum_tab, tab3, 'outer', on=inter_group_by)

        # summarize all interaction sites
        tab4 = inter_morph_df[inter_group_by + ['volume', 'surface_area']].groupby(inter_group_by).agg(['sum'] + inter_ag_func_standard)
        tab5 = inter_morph_df[inter_group_by+inter_cols].groupby(inter_group_by).agg(inter_ag_func_standard)
        inter_sum_tab = pd.merge(inter_sum_tab, tab4, 'outer', on=inter_group_by)
        inter_sum_tab = pd.merge(inter_sum_tab, tab5, 'outer', on=inter_group_by)

        # Get mask_name and corresponding volume column per group & calculate volume fraction
        mask_names = inter_morph_df.groupby(inter_group_by)['mask_name'].first()
        mask_volume_data = inter_morph_df.groupby(inter_group_by).first().apply(lambda row: row[f"{mask_names.loc[row.name]}_volume"], axis=1)
        inter_sum_tab.insert(inter_sum_tab.columns.get_loc(('volume', 'sum')) + 1, ('volume', 'fraction'), inter_sum_tab[('volume', 'sum')]/mask_volume_data)

        # Ensure all possible interactions are represented (if missing fill with NaN)
        for ind in inter_sum_tab.index.droplevel(4).unique().to_list():
            for row in all_pos:
                if ind+(row,) not in inter_sum_tab.index:
                    inter_sum_tab.loc[ind+(row,)] = np.nan

        # fill NA with 0 for specific columns
        inter_fill_dict = {('sites', 'count'): 0, 
                    ('sites_in_higher_order', 'count'): 0, 
                    ('sites_not_in_higher_order', 'count'): 0,
                    ('volume', 'sum'): 0,
                    ('surface_area', 'sum'): 0,
                    ('volume', 'fraction'): 0}
        inter_sum_tab = inter_sum_tab.fillna(value=inter_fill_dict)

        # if (sites, count) is 1, set mean, median, and std to NaN
        inter_single_site_mask = inter_sum_tab[('sites', 'count')] == 1
        for col in inter_cols+['volume', 'surface_area']:
            inter_sum_tab.loc[inter_single_site_mask, (col, 'std')] = np.nan

        inter_sum_tab.sort_index(inplace=True)

        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-interaction_morphology_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-interaction_morphology_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            inter_sum_tab.to_csv(str(out_path) + f"/{out_prefix}-interaction_morphology_summarystats.csv", mode='x')
            print(f"Exported per-interaction morphology summary statistics (before unstacking) to {out_path}/{out_prefix}-interaction_morphology_summarystats.csv")
        # unstack and format interaction morphology summary table   
        inter_morph_final = inter_sum_tab.unstack(-1)
        inter_morph_final.columns = ["_".join((col_name[1], col_name[-1], col_name[0])) for col_name in inter_morph_final.columns.to_flat_index()]
        inter_morph_final.columns = [col.replace('sum', 'total') for col in inter_morph_final.columns]
        inter_morph_final.columns = [col.replace(col, 'mask_volume') if 'mask' in col else col for col in inter_morph_final.columns]
        inter_morph_final = inter_morph_final.loc[:, ~inter_morph_final.columns.duplicated()]

        # combine count/volume and morphology summaries
        inter_morph_final = pd.merge(inter_morph_final, inter_count_vol_final, on=["dataset", "image_name", "mask_name", "scale"])
        final_combo_tab = pd.concat([final_combo_tab, inter_morph_final], axis=1)
    else:
        inter_morph_final = None


    ###################################
    # Summarize interaction labels data
    ###################################
    if inter_labs_df is not None:
        ### summarize interaction site counts
        # define summarization paramters
        labs_group_by = ["dataset", "image_name", "mask_name", "object"]

        # summarize counts of interaction sites per image
        labs_tab1 = inter_labs_df[labs_group_by + ['ID']].groupby(labs_group_by).agg(['count'])
        labs_tab1.rename(columns={'ID': 'sites'}, inplace=True)
        labs_tab2 = inter_labs_df.copy()[inter_labs_df['in_higher_order'] == True][labs_group_by + ['ID']].groupby(labs_group_by).agg(['count'])
        labs_tab2.rename(columns={'ID': 'sites_in_higher_order'}, inplace=True)
        labs_tab3 = inter_labs_df.copy()[inter_labs_df['in_higher_order'] == False][labs_group_by + ['ID']].groupby(labs_group_by).agg(['count'])
        labs_tab3.rename(columns={'ID': 'sites_not_in_higher_order'}, inplace=True)
        labs_inter_sum_tab = pd.merge(labs_tab1, labs_tab2, 'outer', on=labs_group_by)
        labs_inter_sum_tab = pd.merge(labs_inter_sum_tab, labs_tab3, 'outer', on=labs_group_by)

        # Ensure all possible interactions (all_pos) are represented (if missing fill with NaN):
        for ind in labs_inter_sum_tab.index.droplevel(3).unique().to_list():
            for row in all_pos:
                if ind+(row,) not in labs_inter_sum_tab.index:
                    labs_inter_sum_tab.loc[ind+(row,)] = np.nan

        # fill NA with 0 for specific columns
        labs_fill_dict = {('sites', 'count'): 0, 
                    ('sites_in_higher_order', 'count'): 0, 
                    ('sites_not_in_higher_order', 'count'): 0}
        labs_inter_sum_tab = labs_inter_sum_tab.fillna(value=labs_fill_dict)

        labs_inter_sum_tab.sort_index(inplace=True)

        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-interaction_labels_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-interaction_labels_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            labs_inter_sum_tab.to_csv(str(out_path) + f"/{out_prefix}-interaction_labels_summarystats.csv", mode='x')
            print(f"Exported per-interaction count summary statistics (before unstacking) to {out_path}/{out_prefix}-interaction_labels_summarystats.csv")
            print("NOTE: The interactions counts dataframe will not be included into the combined summary stats table with the other analyses.")



    #########################################
    # Summarize interaction distribution data
    ########################################
    if org_dist_df is not None:
        combo_dist_df = org_dist_df
        if inter_dist_df is not None:
            combo_dist_df = pd.concat([combo_dist_df, inter_dist_df])
    else:
        combo_dist_df = None
        
    if combo_dist_df is not None:
        # mask name checker
        mask_name = "whole_image" if mask_name is None else mask_name

        # extract centering object metrics
        if 'XY_center_vox_cnt_perbin' in list(combo_dist_df.columns): # if there is a centering object
            nuc_dist_df = combo_dist_df[["dataset", "image_name", "mask_name", 'scale',
                                "XY_bins", "XY_center_vox_cnt_perbin", f"XY_{mask_name}_vox_cnt_perbin", "XY_center_cv_perbin",
                                "XY_wedges", "XY_center_vox_cnt_perwedge", f"XY_{mask_name}_vox_cnt_perwedge",
                                "Z_slices", "Z_center_vox_cnt", f"Z_{mask_name}_vox_cnt"]].drop_duplicates(subset=['dataset', 'image_name'])
            nuc_dist_df.columns = nuc_dist_df.columns.str.replace('center', 'obj', regex=False)
            nuc_dist_df.insert(loc=3,column='object',value='nuc')
            nuc_dist_df.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object'], inplace=True)

            # select relevant columns from dist dataset
            select_dist_df = combo_dist_df[list(nuc_dist_df.reset_index().columns)]
            select_dist_df.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object'], inplace=True)
            combo_dist_df = pd.concat([nuc_dist_df, select_dist_df], axis=0)
        else: # if there is not a centering object
            combo_dist_df.set_index(['dataset', 'image_name', "mask_name", 'scale', 'object'], inplace=True)

        # loop through each row of data and calculate histogram statistics
        hist_dfs = []
        for ind in combo_dist_df.index:
            selection = combo_dist_df.loc[[ind]].reset_index()
            bins_df = pd.DataFrame()
            wedges_df = pd.DataFrame()
            Z_df = pd.DataFrame()
            CV_df = pd.DataFrame()

            bins_df[['bins', 'masks', 'obj']] = selection[['XY_bins', f'XY_{mask_name}_vox_cnt_perbin', 'XY_obj_vox_cnt_perbin']]
            wedges_df[['bins', 'masks', 'obj']] = selection[['XY_wedges', f'XY_{mask_name}_vox_cnt_perwedge', 'XY_obj_vox_cnt_perwedge']]
            Z_df[['bins', 'masks', 'obj']] = selection[['Z_slices', f'Z_{mask_name}_vox_cnt', 'Z_obj_vox_cnt']]
            CV_df[['XY_obj_cv_perbin']] = selection[['XY_obj_cv_perbin']]

            dfs = [selection[['dataset', 'image_name', 'mask_name','scale', 'object']].reset_index()]

            # for each group of data, calculate histogram statistics
            for df, prefix in zip([bins_df, wedges_df, Z_df, CV_df], ["XY_bins_", "XY_wedges_", "Z_slices_", "CV_perbin_"]):
                if prefix != "CV_perbin_":
                    single_df = pd.DataFrame(list(zip(df["bins"].values[0][1:-1].split(", "), 
                                                    df["obj"].values[0][1:-1].split(", "), 
                                                    df["masks"].values[0][1:-1].split(", "))), columns =['bins', 'obj', 'mask']).astype(int)
                    
                    if "Z_" in prefix:
                        single_df =  single_df.drop(single_df[single_df['mask'] == 0].index)
                        single_df['bins'] = (single_df["bins"]/max(single_df.bins)*9.99).apply(np.floor)+1
                        single_df = single_df.groupby("bins").agg(['sum']).reset_index()
                        single_df.columns = ['bins',"obj","mask"]
                
                    single_df['mask_fract'] = single_df['mask']/single_df['mask'].max()
                    # single_df['obj_normed_tocell'] = (single_df["obj"]*single_df["mask_fract"]).fillna(0)
                    single_df['obj_perc_per_bin'] = (single_df["obj"] / single_df["obj"].sum())*100
                    single_df['obj_portion_normed_tobin'] = (single_df["obj_perc_per_bin"]/single_df["mask_fract"]).fillna(0)

                    sumstats_df = pd.DataFrame()

                    s = single_df['bins'].repeat(single_df['obj_portion_normed_tobin']*100)

                    sumstats_df['hist_mean']=[s.mean()]
                    sumstats_df['hist_median']=[s.median()]
                    if single_df['obj_portion_normed_tobin'].sum() != 0: sumstats_df['hist_mode']=[s.mode().iloc[0]]
                    else: sumstats_df['hist_mode']=['NaN']
                    sumstats_df['hist_min']=[s.min()]
                    sumstats_df['hist_max']=[s.max()]
                    sumstats_df['hist_range']=[s.max() - s.min()]
                    sumstats_df['hist_stdev']=[s.std()]
                    sumstats_df['hist_skew']=[s.skew()]
                    sumstats_df['hist_kurtosis']=[s.kurtosis()]
                    sumstats_df['hist_var']=[s.var()]
                    sumstats_df.columns = [prefix+col for col in sumstats_df.columns]
                    sumstats_df.reset_index(drop=True, inplace=True)

                    dfs.append(sumstats_df)
                    
                if prefix == 'CV_perbin_':
                    CV_df = pd.DataFrame(list(zip(df["XY_obj_cv_perbin"].values[0][1:-1].split(", "))), columns =['CV']).astype(float)
                    sumstats_CV_df = pd.DataFrame()
                    sumstats_CV_df['XY_bin_CV_mean'] = CV_df.mean()
                    sumstats_CV_df['XY_bin_CV_median'] = CV_df.median()
                    sumstats_CV_df['XY_bin_CV_std'] = CV_df.std()
                    sumstats_CV_df.reset_index(drop=True, inplace=True)
                    sumstats_df = pd.concat([sumstats_df, sumstats_CV_df], axis=1)

                    dfs.append(sumstats_df)
                
            combined_df = pd.concat(dfs,axis=1).drop(columns="index")
            combined_df.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object'], inplace=True)
            hist_dfs.append(combined_df)

        dist_summary = pd.concat(hist_dfs).sort_values(by=['dataset', 'image_name', 'mask_name', 'scale', 'object'])

        # Ensure all possible interactions (all_pos) are represented (if missing fill with NaN):
        for ind in dist_summary.index.droplevel(4).unique().to_list():
            for row in all_pos+organelle_names:
                if ind+(row,) not in dist_summary.index:
                    dist_summary.loc[ind+(row,)] = np.nan

        dist_summary.reset_index(inplace=True)

        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-distribution_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-distribution_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            dist_summary.to_csv(str(out_path) + f"/{out_prefix}-distribution_summarystats.csv", mode='x')
            print(f"Exported per-object (interaction or organelle) distribution summary statistics (before unstacking) to {out_path}/{out_prefix}-distribution_summarystats.csv")
        
        # unstack and format interaction distribution summary table
        dist_final = dist_summary.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object']).unstack(-1)
        dist_final.columns = ["_".join((col_name[1], col_name[0])) for col_name in dist_final.columns.to_flat_index()]
        dist_final = dist_final.reset_index()
        dist_final.set_index(['dataset', 'image_name', 'mask_name', 'scale'], inplace=True)

        # combine with previous summary table
        final_combo_tab = pd.concat([final_combo_tab, dist_final], axis=1)
    else:
        final_combo_tab = final_combo_tab


    ###################################
    # Summarize interaction degree data
    ###################################
    if inter_degree_df is not None:
        # degree_df not exported before unstacking becuase it is already summarized per image originally

        # unstack and format interaction degree summary table
        inter_degree_final = inter_degree_df.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object']).unstack(-1)
        inter_degree_final.columns = ["_".join((col_name[1], col_name[0])) for col_name in inter_degree_final.columns.to_flat_index()]
        inter_degree_final.reset_index(inplace=True)
        inter_degree_final.set_index(['dataset', 'image_name', 'mask_name', 'scale'], inplace=True)
  
        # combine with previous summary table
        final_combo_tab = pd.concat([final_combo_tab, inter_degree_final], axis=1)
    else:
        final_combo_tab = final_combo_tab

    

    ##########################
    # Export combined results
    ##########################
    if (Path(out_path) / f"{out_prefix}-combined_summarystats.csv").exists():
        raise FileExistsError(f"CAUTION: {out_prefix}-combined_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
    else:
        final_combo_tab.to_csv(str(out_path) + f"/{out_prefix}-combined_summarystats.csv", mode='x')
        print(f"Exported ALL summary statistics combined to {out_path}/{out_prefix}-combined_summarystats.csv")


    print(f"Summary statistics are complete.")
    return final_combo_tab